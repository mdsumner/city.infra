// Segment mapper: browser app. Reads the tas-transport-web GeoParquet from
// Source Cooperative with hyparquet (pure JS, HTTP range requests, no server)
// and draws with Leaflet (global L, loaded by the page).
//
// On load it reads the small attribute columns for every segment (~10 MB:
// ids, names, length, node ids, bbox); searches then run in memory, and only
// the geometry pages for matching rows are fetched (offset index range reads).
import * as core from "./core.js";
import { asyncBufferFromUrl, cachedAsyncBuffer, parquetMetadataAsync, parquetReadObjects }
  from "https://cdn.jsdelivr.net/npm/hyparquet@1.31.5/+esm";
import { compressors } from "https://cdn.jsdelivr.net/npm/hyparquet-compressors@1.1.2/+esm";

const ATTR_COLS = ["TRANSEG_ID", "PRI_NAME", "SEC_NAME", "COMP_LEN", "TRANS_TYPE", "node0", "node1", "geometry_bbox"];
const ROW_GAP = 512;   // merge geometry reads whose row ranges are closer than this
const BASE = "https://data.source.coop/mdsumner/tas-transport-web/";
const STORE_KEY = "city.infra.segment-mapper.v1";
const AREAS = {
  hobart: { label: "Greater Hobart", bbox: [147.0, -43.2, 147.6, -42.65] },
  view: { label: "Current map view", bbox: null },
  state: { label: "All Tasmania", bbox: null },
};
const STYLE = {
  candidate: { color: "#777777", weight: 4, opacity: 0.75 },
  context: { color: "#bbbbbb", weight: 3, opacity: 0.6, dashArray: "4 4" },
  selected: { color: "#00a33a", weight: 7, opacity: 0.95 },
};
const SAVED_COLOURS = ["#2a6fdb", "#d1495b", "#8e44ad", "#e08e0b", "#00798c", "#6a994e", "#b5651d", "#c2185b"];

const $ = (id) => document.getElementById(id);
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

const state = {
  snapshot: "", url: "", file: null, metadata: null,
  segs: [],                // attribute rows, index = parquet row
  geom: new Map(),         // row -> Uint8Array WKB
  candidates: [],          // [{seg, inRange, order}]
  selected: new Set(),     // TRANSEG_IDs
  layers: new Map(),       // id -> polyline
  mappings: loadMappings(),
  query: "",
};

// ---- persistence ------------------------------------------------------------
function loadMappings() {
  try { return JSON.parse(localStorage.getItem(STORE_KEY) || "{}"); } catch { return {}; }
}
function persist() {
  try { localStorage.setItem(STORE_KEY, JSON.stringify(state.mappings)); } catch { /* private mode etc. */ }
}

// ---- status -------------------------------------------------------------------
function status(msg, kind = "info") {
  const el = $("sm-status");
  el.textContent = msg;
  el.dataset.kind = kind;
}

// ---- data -----------------------------------------------------------------------
async function initData() {
  const latest = await (await fetch(BASE + "latest.json", { cache: "no-cache" })).json();
  state.snapshot = latest.snapshot;
  state.url = BASE + latest.path + "list_transport_segments.parquet";
  $("sm-snapshot").textContent = latest.snapshot;

  state.file = cachedAsyncBuffer(await asyncBufferFromUrl({ url: state.url }));
  state.metadata = await parquetMetadataAsync(state.file, { geoparquet: false });
  const rows = await parquetReadObjects({
    file: state.file, metadata: state.metadata, compressors, geoparquet: false, columns: ATTR_COLS,
  });
  state.segs = rows.map((r, row) => ({
    row, id: Number(r.TRANSEG_ID), name: r.PRI_NAME, sec: r.SEC_NAME, len: Number(r.COMP_LEN),
    type: r.TRANS_TYPE, node0: Number(r.node0), node1: Number(r.node1), bbox: r.geometry_bbox,
  }));
  state.byId = new Map(state.segs.map((s) => [s.id, s]));
}

function inBbox(s, bbox) {
  if (!bbox) return true;
  const b = s.bbox;
  return b.xmax >= bbox[0] && b.xmin <= bbox[2] && b.ymax >= bbox[1] && b.ymin <= bbox[3];
}

function currentBbox() {
  const area = $("sm-area").value;
  if (area === "view") {
    const b = map.getBounds();
    return [b.getWest(), b.getSouth(), b.getEast(), b.getNorth()];
  }
  return AREAS[area].bbox;
}

function byName(q, bbox) {
  return state.segs.filter((s) => inBbox(s, bbox) && core.nameMatches(s, q));
}

function byIds(ids) {
  return ids.map((t) => state.byId.get(Number(t))).filter(Boolean);
}

// Fetch WKB for the given segments (by parquet row), merging nearby rows into
// one range read; results are cached.
async function loadGeometry(segs) {
  const rows = [...new Set(segs.map((s) => s.row))].filter((r) => !state.geom.has(r)).sort((a, b) => a - b);
  const ranges = [];
  for (const r of rows) {
    const last = ranges.at(-1);
    if (last && r - last[1] < ROW_GAP) last[1] = r + 1; else ranges.push([r, r + 1]);
  }
  await Promise.all(ranges.map(async ([rowStart, rowEnd]) => {
    const got = await parquetReadObjects({
      file: state.file, metadata: state.metadata, compressors, geoparquet: false, utf8: false,
      columns: ["geometry"], rowStart, rowEnd, useOffsetIndex: true,
    });
    got.forEach((g, i) => state.geom.set(rowStart + i, g.geometry));
  }));
}

// ---- search -----------------------------------------------------------------------
async function search() {
  const q = $("sm-query").value;
  const parts = core.parseQuery(q);
  if (!parts.length) { status("Type a street name, e.g. Argyle [Federal, Burnett]"); return; }
  const bbox = currentBbox();
  const ctx = $("sm-context").checked ? 3 : 0;
  status("Searching...");
  $("sm-search").disabled = true;
  try {
    const out = [], seen = new Set(), notes = [];
    for (const p of parts) {
      const main = byName(p.main, bbox);
      let cross = [];
      if (p.cross.length) {
        cross = p.cross.map((c) => (c ? byName(c, bbox) : null));
        if (cross.length === 1) cross.push(null);
      }
      const r = core.resolvePart(main, cross, ctx);
      r.messages.forEach((m) => notes.push(`${p.main}: ${m}`));
      for (const x of r.segs) if (!seen.has(x.seg.id)) { seen.add(x.seg.id); out.push(x); }
    }
    out.forEach((x, i) => (x.order = i + 1));
    await loadGeometry(out.map((x) => x.seg));
    state.candidates = out;
    state.query = q;
    state.selected = new Set();
    drawCandidates();
    renderTable();
    renderSelection();
    const n = out.filter((x) => x.inRange).length;
    status(`${n} segments in range${out.length > n ? `, ${out.length - n} context` : ""} across ${parts.length} part(s)` +
           (notes.length ? ". " + notes.join("; ") : ""), notes.length ? "warn" : "info");
  } catch (e) {
    console.error(e);
    status("Search failed: " + e.message, "error");
  } finally {
    $("sm-search").disabled = false;
  }
}

// ---- map --------------------------------------------------------------------------------
let map, candLayer, savedLayer;

function initMap() {
  map = L.map("sm-map", { preferCanvas: false }).setView([-42.88, 147.33], 14);
  L.tileLayer("https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png", {
    maxZoom: 20, subdomains: "abcd",
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>',
  }).addTo(map);
  savedLayer = L.layerGroup().addTo(map);
  candLayer = L.layerGroup().addTo(map);
  map.attributionControl.addAttribution('Network: "LIST Transport Segments" from theLIST &copy; State of Tasmania');
}

function styleFor(x) {
  if (state.selected.has(x.seg.id)) return STYLE.selected;
  return x.inRange ? STYLE.candidate : STYLE.context;
}

function latlngs(seg) {
  return core.wkbLines(state.geom.get(seg.row)).map((line) => line.map(([x, y]) => [y, x]));
}

function drawCandidates() {
  candLayer.clearLayers();
  state.layers.clear();
  const showCtx = $("sm-context").checked;
  const bounds = [];
  for (const x of state.candidates) {
    if (!x.inRange && !showCtx) continue;
    const ll = latlngs(x.seg);
    const pl = L.polyline(ll, styleFor(x)).addTo(candLayer);
    pl.bindTooltip(`#${x.order} ${esc(x.seg.name || "(unnamed)")}${x.inRange ? "" : " (context)"}<br>` +
                   `ID ${x.seg.id}, ${Math.round(x.seg.len)} m`, { sticky: true });
    pl.on("click", () => toggle(x.seg.id));
    state.layers.set(x.seg.id, { pl, x });
    ll.forEach((l) => bounds.push(...l));
  }
  if (bounds.length) map.fitBounds(L.latLngBounds(bounds).pad(0.15));
}

function restyle() {
  for (const { pl, x } of state.layers.values()) {
    pl.setStyle(styleFor(x));
    if (state.selected.has(x.seg.id)) pl.bringToFront();
  }
}

function toggle(id) {
  if (state.selected.has(id)) state.selected.delete(id); else state.selected.add(id);
  restyle();
  renderSelection();
  const cb = document.querySelector(`#sm-table input[data-id="${id}"]`);
  if (cb) cb.checked = state.selected.has(id);
}

// ---- table and selection ---------------------------------------------------------------
function renderTable() {
  const showCtx = $("sm-context").checked;
  const rows = state.candidates.filter((x) => x.inRange || showCtx).map((x) => `
    <tr class="${x.inRange ? "" : "sm-ctx"}">
      <td><input type="checkbox" data-id="${x.seg.id}" ${state.selected.has(x.seg.id) ? "checked" : ""} aria-label="select ${x.seg.id}"></td>
      <td>${x.order}</td><td>${x.seg.id}</td><td>${esc(x.seg.name)}</td><td>${esc(x.seg.sec)}</td>
      <td class="sm-num">${Math.round(x.seg.len)}</td><td>${x.inRange ? "" : "context"}</td>
    </tr>`).join("");
  $("sm-table").innerHTML = `<table class="table table-sm"><thead><tr><th></th><th>#</th><th>ID</th><th>Name</th><th>Secondary</th><th class="sm-num">m</th><th></th></tr></thead><tbody>${rows}</tbody></table>`;
  $("sm-table").querySelectorAll("input[data-id]").forEach((cb) =>
    cb.addEventListener("change", () => toggle(Number(cb.dataset.id))));
}

function renderSelection() {
  const ids = [...state.selected];
  const byId = new Map(state.candidates.map((x) => [x.seg.id, x.seg]));
  const total = ids.reduce((a, id) => a + (byId.get(id)?.len || 0), 0);
  $("sm-selection").textContent = ids.length ? `${ids.length} segments, ${Math.round(total)} m` : "None selected";
  $("sm-save").disabled = !ids.length;
}

function selectAllInRange() {
  state.candidates.filter((x) => x.inRange).forEach((x) => state.selected.add(x.seg.id));
  restyle(); renderTable(); renderSelection();
}

// ---- save, list, export ---------------------------------------------------------------
function save() {
  const id = $("sm-id").value.trim();
  if (!id) { status("Give the grouping an ID first", "warn"); return; }
  if (!state.selected.size) return;
  if (state.mappings[id] && !confirm(`Replace existing grouping "${id}"?`)) return;
  state.mappings[id] = {
    id, description: $("sm-desc").value.trim(), transeg_ids: [...state.selected],
    snapshot: state.snapshot, query: state.query, saved: new Date().toISOString().replace(/\.\d+Z$/, "Z"),
  };
  persist();
  renderSaved();
  status(`Saved ${id} with ${state.selected.size} segments`);
}

function renderSaved() {
  const ms = Object.values(state.mappings);
  $("sm-saved").innerHTML = ms.length ? ms.map((m, i) => `
    <li><span class="sm-swatch" style="background:${SAVED_COLOURS[i % SAVED_COLOURS.length]}"></span>
      <b>${esc(m.id)}</b> ${m.transeg_ids.length} segs${m.snapshot && m.snapshot !== state.snapshot ? ` <span class="sm-warn" title="made on snapshot ${esc(m.snapshot)}">(${esc(m.snapshot)})</span>` : ""}
      <br><small>${esc(m.description)}</small>
      <span class="sm-actions">
        <button class="btn btn-link btn-sm p-0" data-act="load" data-id="${esc(m.id)}">edit</button>
        <button class="btn btn-link btn-sm p-0" data-act="del" data-id="${esc(m.id)}">delete</button>
      </span></li>`).join("") : "<li><small>None yet</small></li>";
  $("sm-saved").querySelectorAll("button[data-act]").forEach((b) => b.addEventListener("click", () => {
    const m = state.mappings[b.dataset.id];
    if (b.dataset.act === "del") {
      if (confirm(`Delete grouping "${m.id}"?`)) { delete state.mappings[m.id]; persist(); renderSaved(); drawSaved(); }
    } else loadMapping(m);
  }));
}

async function loadMapping(m) {
  status(`Loading ${m.id}...`);
  const segs = byIds(m.transeg_ids);
  await loadGeometry(segs);
  state.candidates = segs.map((seg, i) => ({ seg, inRange: true, order: i + 1 }));
  state.selected = new Set(segs.map((s) => s.id));
  state.query = m.query || "";
  $("sm-query").value = m.query || "";
  $("sm-id").value = m.id;
  $("sm-desc").value = m.description || "";
  drawCandidates(); renderTable(); renderSelection();
  const missing = m.transeg_ids.length - segs.length;
  status(missing ? `${m.id}: ${missing} segment ids not in snapshot ${state.snapshot}` : `Loaded ${m.id}`, missing ? "warn" : "info");
}

async function drawSaved() {
  savedLayer.clearLayers();
  if (!$("sm-show-saved").checked) return;
  const ms = Object.values(state.mappings);
  const ids = [...new Set(ms.flatMap((m) => m.transeg_ids))];
  if (!ids.length) return;
  const found = byIds(ids);
  await loadGeometry(found);
  const segs = new Map(found.map((s) => [s.id, s]));
  ms.forEach((m, i) => {
    const colour = SAVED_COLOURS[i % SAVED_COLOURS.length];
    for (const t of m.transeg_ids) {
      const s = segs.get(Number(t));
      if (!s) continue;
      L.polyline(latlngs(s), { color: colour, weight: 5, opacity: 0.7 })
        .bindTooltip(`<b>${esc(m.id)}</b><br>${esc(m.description)}`, { sticky: true })
        .addTo(savedLayer);
    }
  });
}

function exportCSV() {
  const blob = new Blob([core.mappingsToCSV(state.mappings)], { type: "text/csv" });
  const a = Object.assign(document.createElement("a"), {
    href: URL.createObjectURL(blob), download: `segment_mappings_${new Date().toISOString().slice(0, 10)}.csv`,
  });
  document.body.appendChild(a); a.click(); a.remove();
}

async function importCSV(file) {
  try {
    const incoming = core.csvToMappings(await file.text());
    const n = Object.keys(incoming).length;
    state.mappings = core.mergeMappings(state.mappings, incoming);
    persist(); renderSaved(); drawSaved();
    status(`Imported ${n} groupings`);
  } catch (e) {
    status("Import failed: " + e.message, "error");
  }
}

// ---- wire up ------------------------------------------------------------------------
async function main() {
  initMap();
  renderSaved();
  $("sm-search").addEventListener("click", search);
  $("sm-query").addEventListener("keydown", (e) => { if (e.key === "Enter") search(); });
  $("sm-context").addEventListener("change", () => { drawCandidates(); renderTable(); });
  $("sm-all").addEventListener("click", selectAllInRange);
  $("sm-clear").addEventListener("click", () => { state.selected.clear(); restyle(); renderTable(); renderSelection(); });
  $("sm-save").addEventListener("click", save);
  $("sm-export").addEventListener("click", exportCSV);
  $("sm-import").addEventListener("change", (e) => { if (e.target.files[0]) importCSV(e.target.files[0]); e.target.value = ""; });
  $("sm-show-saved").addEventListener("change", drawSaved);
  $("sm-search").disabled = true;
  status("Loading network attributes (about 10 MB, once)...");
  try {
    await initData();
    $("sm-search").disabled = false;
    status("Ready. Try: Argyle [Federal, Burnett]");
    window.__smReady = true;
  } catch (e) {
    console.error(e);
    status("Could not load the network: " + e.message, "error");
  }
}

main();
