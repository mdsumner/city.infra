// Run with: node --test site/mapper/test/
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  parseQuery, nameMatches, resolvePart, orderSegments, wkbLines,
  mappingsToCSV, csvToMappings, nodesOf,
} from "../core.js";

const fixture = JSON.parse(readFileSync(new URL("./hobart-fixture.json", import.meta.url)));
const byName = (q) => fixture.filter((s) => nameMatches(s, q));

test("parseQuery: parts, bounds, open bounds", () => {
  assert.deepEqual(parseQuery(""), []);
  assert.deepEqual(parseQuery("Argyle"), [{ main: "Argyle", cross: [] }]);
  assert.deepEqual(parseQuery("Argyle [Federal, Burnett]"), [{ main: "Argyle", cross: ["Federal", "Burnett"] }]);
  assert.deepEqual(parseQuery("Argyle [Federal,]"), [{ main: "Argyle", cross: ["Federal", ""] }]);
  assert.deepEqual(parseQuery("Argyle [Federal]"), [{ main: "Argyle", cross: ["Federal"] }]);
  assert.deepEqual(
    parseQuery("Macquarie [Darcy, Elboden], Weld, Anglesea [Adelaide, Macquarie]").map((p) => p.main),
    ["Macquarie", "Weld", "Anglesea"]);
});

test("whole street: every matching segment, connected runs in walk order", () => {
  const main = byName("Argyle");
  const r = resolvePart(main, []);
  assert.equal(r.segs.length, main.length);
  assert.ok(r.segs.every((x) => x.inRange));
});

test("bounded range on real Hobart data: Argyle between Federal and Burnett", () => {
  const main = byName("Argyle");
  const r = resolvePart(main, [byName("Federal"), byName("Burnett")], 0);
  assert.deepEqual(r.messages, []);
  assert.ok(r.segs.length > 0);
  // the range is a single path: consecutive segments share a node
  const segs = r.segs.map((x) => x.seg);
  for (let i = 1; i < segs.length; i++) {
    const a = new Set([segs[i - 1].node0, segs[i - 1].node1]);
    assert.ok(a.has(segs[i].node0) || a.has(segs[i].node1), "path is contiguous");
  }
  // and it starts at a Federal node and ends at a Burnett node
  const fed = nodesOf(byName("Federal")), bur = nodesOf(byName("Burnett"));
  const ends = [segs[0].node0, segs[0].node1, segs.at(-1).node0, segs.at(-1).node1];
  assert.ok(ends.some((n) => fed.has(n)));
  assert.ok(ends.some((n) => bur.has(n)));
});

test("bounded range adds context segments flagged out of range", () => {
  const main = byName("Argyle");
  const r0 = resolvePart(main, [byName("Federal"), byName("Burnett")], 0);
  const r3 = resolvePart(main, [byName("Federal"), byName("Burnett")], 3);
  assert.equal(r3.segs.filter((x) => x.inRange).length, r0.segs.length);
  assert.ok(r3.segs.some((x) => !x.inRange));
});

test("open bound runs to the far end of the street", () => {
  const main = byName("Argyle");
  const r = resolvePart(main, [byName("Federal"), null], 0);
  assert.ok(r.segs.length > 0);
  assert.ok(r.segs.every((x) => x.inRange));
});

test("cross street that does not meet the main street is reported", () => {
  const main = byName("Argyle");
  const r = resolvePart(main, [byName("Liverpool Crescent"), byName("Burnett")], 0);
  assert.ok(r.messages.some((m) => /does not meet/.test(m)));
});

test("orderSegments covers all segments exactly once", () => {
  const main = byName("Macquarie");
  const o = orderSegments(main);
  assert.equal(o.length, main.length);
  assert.equal(new Set(o).size, main.length);
});

test("orderSegments handles a self-loop", () => {
  const segs = [
    { id: 1, len: 10, node0: 1, node1: 2 },
    { id: 2, len: 10, node0: 2, node1: 2 },
    { id: 3, len: 10, node0: 2, node1: 3 },
  ];
  assert.equal(orderSegments(segs).length, 3);
});

test("wkbLines decodes LineString and MultiLineString, both byte orders", () => {
  const ls = (le, pts) => {
    const b = new DataView(new ArrayBuffer(9 + pts.length * 16));
    b.setUint8(0, le ? 1 : 0); b.setUint32(1, 2, le); b.setUint32(5, pts.length, le);
    pts.forEach(([x, y], i) => { b.setFloat64(9 + i * 16, x, le); b.setFloat64(17 + i * 16, y, le); });
    return new Uint8Array(b.buffer);
  };
  assert.deepEqual(wkbLines(ls(true, [[147.3, -42.88], [147.31, -42.89]])), [[[147.3, -42.88], [147.31, -42.89]]]);
  assert.deepEqual(wkbLines(ls(false, [[1, 2], [3, 4]])), [[[1, 2], [3, 4]]]);
  const a = ls(true, [[1, 2], [3, 4]]), c = ls(true, [[5, 6], [7, 8]]);
  const head = new DataView(new ArrayBuffer(9));
  head.setUint8(0, 1); head.setUint32(1, 5, true); head.setUint32(5, 2, true);
  const multi = new Uint8Array([...new Uint8Array(head.buffer), ...a, ...c]);
  assert.deepEqual(wkbLines(multi), [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]);
});

test("CSV round trip, quoting, and reading the R app's format", () => {
  const m = { CP001: { id: "CP001", description: 'Argyle "bike" lane, north', transeg_ids: [12345, 12346],
                       snapshot: "2025-07-10", query: "Argyle [Federal, Burnett]", saved: "2026-11-07T00:00:00Z" } };
  const back = csvToMappings(mappingsToCSV(m));
  assert.deepEqual(back, m);
  const r = csvToMappings('id,description,transeg_id\nCP001,"Argyle St bike lane",12345\nCP001,"Argyle St bike lane",12346\n');
  assert.deepEqual(r.CP001.transeg_ids, [12345, 12346]);
});
