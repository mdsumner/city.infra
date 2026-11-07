// Segment mapper core: pure functions, no DOM, no network.
// Works on segments shaped like rows of tas-transport-web:
//   { id, name, sec, len, node0, node1, ... }
// Topology comes from precomputed node ids, so all matching is integer
// equality (no coordinate tolerances).

// ---- query syntax ----------------------------------------------------------
// "Argyle"                          whole street
// "Argyle [Federal, Burnett]"       Argyle between the Federal and Burnett intersections
// "Argyle [Federal,]"               from the Federal intersection to the far end
// "Macquarie [Darcy, Elboden], Weld, Anglesea [Adelaide, Macquarie]"
// Commas outside brackets separate route parts.
export function parseQuery(query) {
  const q = String(query || "").trim();
  if (!q) return [];
  const parts = [];
  let depth = 0, start = 0;
  for (let i = 0; i <= q.length; i++) {
    const ch = q[i];
    if (ch === "[") depth++;
    else if (ch === "]") depth = Math.max(0, depth - 1);
    if (i === q.length || (ch === "," && depth === 0)) {
      const chunk = q.slice(start, i).trim();
      if (chunk) parts.push(parsePart(chunk));
      start = i + 1;
    }
  }
  return parts;
}

function parsePart(chunk) {
  const m = chunk.match(/^(.+?)\s*\[(.*)\]\s*$/);
  if (!m) return { main: chunk, cross: [] };
  const cross = m[2].split(",").map((s) => s.trim());
  // "[A,]" keeps an empty second bound; "[A]" is the same as "[A,]"
  return { main: m[1].trim(), cross: cross.filter((s, i) => s || i > 0).slice(0, 2) };
}

// ---- name matching ---------------------------------------------------------
export function nameMatches(seg, q) {
  const u = String(q).toUpperCase();
  return (seg.name && seg.name.toUpperCase().includes(u)) ||
         (seg.sec && seg.sec.toUpperCase().includes(u));
}

// ---- graph helpers ---------------------------------------------------------
export function adjacency(segs) {
  const adj = new Map();
  const add = (n, s, other) => {
    if (!adj.has(n)) adj.set(n, []);
    adj.get(n).push({ seg: s, other });
  };
  for (const s of segs) {
    add(s.node0, s, s.node1);
    if (s.node1 !== s.node0) add(s.node1, s, s.node0);
  }
  return adj;
}

export function nodesOf(segs) {
  const out = new Set();
  for (const s of segs) { out.add(s.node0); out.add(s.node1); }
  return out;
}

// Dijkstra over segment length from a set of start nodes.
// Returns { dist: Map(node -> d), prev: Map(node -> {seg, from}) }.
export function dijkstra(segs, starts) {
  const adj = adjacency(segs);
  const dist = new Map(), prev = new Map(), done = new Set();
  const heap = [];
  const push = (d, n) => { heap.push([d, n]); let i = heap.length - 1;
    while (i > 0) { const p = (i - 1) >> 1; if (heap[p][0] <= heap[i][0]) break; [heap[p], heap[i]] = [heap[i], heap[p]]; i = p; } };
  const pop = () => { const top = heap[0], last = heap.pop();
    if (heap.length) { heap[0] = last; let i = 0;
      for (;;) { const l = 2 * i + 1, r = l + 1; let m = i;
        if (l < heap.length && heap[l][0] < heap[m][0]) m = l;
        if (r < heap.length && heap[r][0] < heap[m][0]) m = r;
        if (m === i) break; [heap[m], heap[i]] = [heap[i], heap[m]]; i = m; } }
    return top; };
  for (const n of starts) { dist.set(n, 0); push(0, n); }
  while (heap.length) {
    const [d, n] = pop();
    if (done.has(n)) continue;
    done.add(n);
    for (const { seg, other } of adj.get(n) || []) {
      const nd = d + (Number(seg.len) || 0);
      if (!dist.has(other) || nd < dist.get(other)) {
        dist.set(other, nd); prev.set(other, { seg, from: n }); push(nd, other);
      }
    }
  }
  return { dist, prev };
}

function walkBack(prev, end) {
  const path = [];
  let n = end;
  while (prev.has(n)) { const p = prev.get(n); path.push(p.seg); n = p.from; }
  return path.reverse();
}

// Segments of the connected component containing any of the start nodes.
export function component(segs, starts) {
  const adj = adjacency(segs);
  const seen = new Set(starts), out = new Set(), queue = [...starts];
  while (queue.length) {
    const n = queue.shift();
    for (const { seg, other } of adj.get(n) || []) {
      out.add(seg);
      if (!seen.has(other)) { seen.add(other); queue.push(other); }
    }
  }
  return segs.filter((s) => out.has(s));
}

// Order segments as a walk: start at a leaf (degree-1 node) of each
// component, depth first, components ordered by size (largest first).
export function orderSegments(segs) {
  const adj = adjacency(segs);
  const used = new Set(), out = [];
  const degree = (n) => (adj.get(n) || []).length;
  const comps = [];
  const remaining = new Set(segs);
  while (remaining.size) {
    const first = remaining.values().next().value;
    const comp = component([...remaining], [first.node0]);
    comp.forEach((s) => remaining.delete(s));
    comps.push(comp);
  }
  comps.sort((a, b) => b.length - a.length);
  for (const comp of comps) {
    const ns = [...nodesOf(comp)];
    const start = ns.find((n) => degree(n) === 1) ?? ns[0];
    const stack = [start];
    while (stack.length) {
      const n = stack.pop();
      for (const { seg, other } of [...(adj.get(n) || [])].reverse()) {
        if (used.has(seg)) continue;
        used.add(seg); out.push(seg); stack.push(other);
      }
    }
  }
  return out;
}

// ---- resolve one query part -------------------------------------------------
// main:   segments matching the main street name
// crossLists: for each cross street given, the segments matching its name
// Returns { segs: [{seg, inRange}], messages: [] } with in-range segments in
// walk order, then context segments (k hops either side) flagged inRange=false.
export function resolvePart(main, crossLists = [], context = 3) {
  const messages = [];
  if (!main.length) return { segs: [], messages: ["no segments match the street name"] };
  const mainNodes = nodesOf(main);
  const hits = crossLists.map((cl) => {
    if (cl === null) return null;   // open bound
    const common = [...nodesOf(cl)].filter((n) => mainNodes.has(n));
    return common;
  });
  const bounds = hits.filter((h) => h !== null);
  if (bounds.length === 0) {
    return { segs: orderSegments(main).map((seg) => ({ seg, inRange: true })), messages };
  }
  bounds.forEach((h, i) => { if (!h.length) messages.push(`cross street ${i + 1} does not meet the main street`); });
  const found = bounds.filter((h) => h.length);
  let range = [];
  if (found.length >= 2) {
    const { dist, prev } = dijkstra(main, found[0]);
    const end = found[1].filter((n) => dist.has(n)).sort((a, b) => dist.get(a) - dist.get(b))[0];
    if (end === undefined) {
      messages.push("the two intersections are not connected along the main street");
      return { segs: [], messages };
    }
    range = walkBack(prev, end);
  } else if (found.length === 1) {
    // one bound: from the intersection to the farthest end of the street
    const comp = component(main, found[0]);
    const adj = adjacency(comp);
    const { dist, prev } = dijkstra(comp, found[0]);
    const leaves = [...dist.keys()].filter((n) => (adj.get(n) || []).length === 1);
    const pool = leaves.length ? leaves : [...dist.keys()];
    const end = pool.sort((a, b) => dist.get(b) - dist.get(a))[0];
    range = walkBack(prev, end);
    if (!range.length) messages.push("intersection is at the end of the street");
  } else {
    messages.push("no intersections found; showing the whole street");
    return { segs: orderSegments(main).map((seg) => ({ seg, inRange: true })), messages };
  }
  const inRange = new Set(range);
  const ctx = [];
  if (context > 0) {
    const adj = adjacency(main);
    let frontier = [...nodesOf(range)], seen = new Set(frontier);
    for (let k = 0; k < context; k++) {
      const next = [];
      for (const n of frontier) for (const { seg, other } of adj.get(n) || []) {
        if (inRange.has(seg) || ctx.includes(seg)) continue;
        ctx.push(seg);
        if (!seen.has(other)) { seen.add(other); next.push(other); }
      }
      frontier = next;
    }
  }
  return { segs: [...range.map((seg) => ({ seg, inRange: true })), ...ctx.map((seg) => ({ seg, inRange: false }))], messages };
}

// ---- WKB -> coordinate arrays -------------------------------------------------
// Handles (Multi)LineString, 2D or Z (ISO or EWKB flags), either byte order.
export function wkbLines(buf) {
  const u8 = buf instanceof Uint8Array ? buf : new Uint8Array(buf);
  const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
  let o = 0;
  const lines = [];
  function geom() {
    const le = dv.getUint8(o) === 1; o += 1;
    let t = dv.getUint32(o, le); o += 4;
    let hasZ = false, hasM = false;
    if (t & 0x80000000) hasZ = true;
    if (t & 0x40000000) hasM = true;
    if (t & 0x20000000) { o += 4; }          // EWKB SRID
    t &= 0x0fffffff;
    if (t >= 3000) { hasZ = true; hasM = true; t -= 3000; }
    else if (t >= 2000) { hasM = true; t -= 2000; }
    else if (t >= 1000) { hasZ = true; t -= 1000; }
    const dim = 2 + (hasZ ? 1 : 0) + (hasM ? 1 : 0);
    if (t === 2) {
      const n = dv.getUint32(o, le); o += 4;
      const line = new Array(n);
      for (let i = 0; i < n; i++) {
        line[i] = [dv.getFloat64(o, le), dv.getFloat64(o + 8, le)];
        o += 8 * dim;
      }
      lines.push(line);
    } else if (t === 5) {
      const n = dv.getUint32(o, le); o += 4;
      for (let i = 0; i < n; i++) geom();
    } else {
      throw new Error("unsupported WKB geometry type " + t);
    }
  }
  geom();
  return lines;
}

// ---- mappings <-> CSV ---------------------------------------------------------
// One row per segment: id, description, transeg_id, snapshot, query, saved
export const CSV_FIELDS = ["id", "description", "transeg_id", "snapshot", "query", "saved"];

function csvCell(v) {
  const s = v === null || v === undefined ? "" : String(v);
  return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
}

export function mappingsToCSV(mappings) {
  const rows = [CSV_FIELDS.join(",")];
  for (const m of Object.values(mappings)) {
    for (const t of m.transeg_ids) {
      rows.push([m.id, m.description, t, m.snapshot, m.query, m.saved].map(csvCell).join(","));
    }
  }
  return rows.join("\n") + "\n";
}

export function parseCSV(text) {
  const rows = [];
  let row = [], cell = "", q = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (q) {
      if (c === '"' && text[i + 1] === '"') { cell += '"'; i++; }
      else if (c === '"') q = false;
      else cell += c;
    } else if (c === '"') q = true;
    else if (c === ",") { row.push(cell); cell = ""; }
    else if (c === "\n" || c === "\r") {
      if (c === "\r" && text[i + 1] === "\n") i++;
      row.push(cell); rows.push(row); row = []; cell = "";
    } else cell += c;
  }
  if (cell !== "" || row.length) { row.push(cell); rows.push(row); }
  return rows.filter((r) => r.length > 1 || r[0] !== "");
}

// Accepts both this format and the R app's (id, description, transeg_id).
export function csvToMappings(text) {
  const rows = parseCSV(text);
  if (!rows.length) return {};
  const head = rows[0].map((h) => h.trim().toLowerCase());
  const col = (name) => head.indexOf(name);
  const iId = col("id"), iDesc = col("description"), iT = col("transeg_id");
  if (iId < 0 || iT < 0) throw new Error("CSV needs id and transeg_id columns");
  const out = {};
  for (const r of rows.slice(1)) {
    const id = r[iId];
    if (!id) continue;
    if (!out[id]) out[id] = {
      id, description: iDesc >= 0 ? r[iDesc] : "", transeg_ids: [],
      snapshot: col("snapshot") >= 0 ? r[col("snapshot")] : "",
      query: col("query") >= 0 ? r[col("query")] : "",
      saved: col("saved") >= 0 ? r[col("saved")] : "",
    };
    const t = Number(r[iT]);
    if (Number.isFinite(t) && !out[id].transeg_ids.includes(t)) out[id].transeg_ids.push(t);
  }
  return out;
}

// Merge b into a; an id present in both takes b's version.
export function mergeMappings(a, b) {
  return Object.assign({}, a, b);
}
