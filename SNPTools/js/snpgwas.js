/* =====================================================================
 *  snpgwas.js — GWAS Explorer (Manhattan plot) as a native SNPTools panel.
 *  Registers 'snpgwas'. Loads AFTER core.js + data.js (needs SNPTools, S,
 *  go(), and Data.chromLengths() for shared genome geometry).
 *
 *  Multi-dataset model: data/gwas/manifest.json lists every available
 *  trait × population × publication combination (id, trait, population,
 *  publication, file, threshP, ...). The picker at the top of the page
 *  is a faceted filter over that list — narrowing one dropdown narrows
 *  what the others can offer, so a user can never land on a combination
 *  that doesn't exist. Each combination's SNPs live in their own CSV
 *  (same format as before), fetched lazily and cached per-dataset once
 *  visited. All datasets share one genome layout (chromosome count and
 *  lengths come from Data.chromLengths(), the same table SNPVersity
 *  uses), so only the per-SNP arrays and the significance threshold
 *  change when you switch datasets — the axis geometry does not.
 *
 *  Everything below the picker (canvas rendering, zoom/pan math,
 *  binary-search region lookup, minimap, region panel, CSV export via
 *  Blob + <a download>, "Send to SNPVersity" handoff (region and/or NAM
 *  accessions, independently checkbox-gated) through
 *  window.versityRequest(...)) is unchanged from the single-dataset
 *  version — it just reads the currently active dataset instead of a
 *  hardcoded one.
 * ===================================================================== */
(function () {
  'use strict';

  const CFG = {
    manifestUrl: './data/gwas/manifest.json',
    dataBase:    './data/gwas/',
    chrCount:    10,
  };
  const FACETS = ['trait', 'traitCategory', 'population', 'publication'];
  const FILTERED_SNPS_NOTE = 'Significant markers shown in red. Non-significant SNPs shown in blue, with dark and light blue alternating by chromosome. <br> NOTE: SNPs with raw p values above 0.001 (−log₁₀ P < 3) are not shown, to save memory.';
  /* the hard floor described above, as a number — a chosen threshold looser
     than this can't be honored since those rows were never in the file */
  const FLOOR_NEGLOG = 3;

  /* ------------------------------------------------------------------ *
   *  SHARED GENOME GEOMETRY — computed once from the same chromosome    *
   *  length table SNPVersity uses, so every dataset's axis lines up.    *
   * ------------------------------------------------------------------ */
  let GEOM = null; // {chrLen, chrOffset, TOTAL_CUM}

  function ensureGeometry() {
    if (GEOM) return GEOM;
    if (typeof Data === 'undefined' || typeof Data.chromLengths !== 'function') {
      throw new Error('Data.chromLengths() is unavailable — cannot size the genome axis.');
    }
    const ref = Data.chromLengths();
    const CHR_COUNT = CFG.chrCount;
    const chrLen = new Float64Array(CHR_COUNT + 1);
    for (let c = 1; c <= CHR_COUNT; c++) chrLen[c] = ref['chr' + c] || 0;
    let totalLen = 0;
    for (let c = 1; c <= CHR_COUNT; c++) totalLen += chrLen[c];
    const GAP = (totalLen / CHR_COUNT) * 0.03;
    const chrOffset = new Float64Array(CHR_COUNT + 1);
    let cum = 0;
    for (let c = 1; c <= CHR_COUNT; c++) { chrOffset[c] = cum; cum += chrLen[c] + GAP; }
    GEOM = { chrLen: chrLen, chrOffset: chrOffset, TOTAL_CUM: cum - GAP };
    return GEOM;
  }

  /* ------------------------------------------------------------------ *
   *  MANIFEST — the catalog of available trait × population ×          *
   *  publication combinations.                                          *
   * ------------------------------------------------------------------ */
  let MANIFEST = null;
  let manifestPromise = null;
  let manifestError = null;

  function loadManifest() {
    if (MANIFEST) return Promise.resolve(MANIFEST);
    if (manifestPromise) return manifestPromise;
    manifestPromise = fetch(CFG.manifestUrl, { cache: 'no-store' }).then(function (r) {
      if (!r.ok) throw new Error('HTTP ' + r.status + ' fetching ' + CFG.manifestUrl);
      return r.json();
    }).then(function (list) {
      if (!Array.isArray(list) || !list.length) throw new Error('No datasets listed in ' + CFG.manifestUrl);
      MANIFEST = list;
      resolveSelection();
      if (activeEntry) loadDataset(activeEntry);
      return MANIFEST;
    }).catch(function (err) {
      manifestError = err;
      manifestPromise = null;
      throw err;
    });
    return manifestPromise;
  }

  /* ------------------------------------------------------------------ *
   *  FACETED SELECTION — selection[facet] narrows candidateEntries();   *
   *  optionsFor(facet) computes what's still reachable in that facet    *
   *  given the OTHER set facets, so every offered option leads to at    *
   *  least one real dataset (no dead-end combinations).                 *
   * ------------------------------------------------------------------ */
  let selection = { trait: null, traitCategory: null, measure: null, population: null, publication: null };
  let activeEntry = null;   // resolved manifest entry currently shown (or being loaded)

  function entryMatches(e, sel, excludeFacet) {
    return FACETS.every(function (f) {
      if (f === excludeFacet) return true;
      return sel[f] == null || e[f] === sel[f];
    });
  }
  function candidateEntries() {
    return MANIFEST.filter(function (e) { return entryMatches(e, selection, null); });
  }
  function optionsFor(facet) {
    const seen = new Map();
    MANIFEST.forEach(function (e) {
      if (entryMatches(e, selection, facet)) seen.set(e[facet], (seen.get(e[facet]) || 0) + 1);
    });
    return Array.from(seen.keys()).sort(function (a, b) { return a.localeCompare(b); })
      .map(function (v) { return { value: v, count: seen.get(v) }; });
  }

  /* Narrow to a unique dataset automatically; otherwise only drop the
     active dataset if it no longer matches the (broadened/changed)
     selection — broadening a filter shouldn't yank away a plot that's
     still a valid match, it should just surface the other options too. */
  function resolveSelection() {
    const cands = candidateEntries();
    if (cands.length === 1) { setActiveEntry(cands[0]); return; }
    if (activeEntry && !cands.some(function (e) { return e.id === activeEntry.id; })) {
      activeEntry = null; DATA = null; datasetError = null; currentRegion = null;
    }
  }

  function setActiveEntry(entry) {
    if (activeEntry && activeEntry.id === entry.id && (DATA || datasetPromiseId === entry.id)) return;
    activeEntry = entry;
    DATA = null; datasetError = null; currentRegion = null;
    selection = { trait: entry.trait, traitCategory: entry.traitCategory, measure: entry.measure, population: entry.population, publication: entry.publication };
    loadDataset(entry);
  }

  function onFacetChange(facet, value) {
    selection[facet] = value || null;
    resolveSelection();
    render(document.getElementById('page'));
  }

  /* ------------------------------------------------------------------ *
   *  DATASET LOAD + PARSE — fetched lazily, cached per manifest entry   *
   *  id so revisiting a dataset in the same session is instant.         *
   * ------------------------------------------------------------------ */
  let DATA = null;
  let datasetCache = new Map();   // entry.id -> parsed DATA
  let datasetPromiseId = null;    // entry.id currently in flight (guards a stale fetch from clobbering a newer pick)
  let datasetError = null;

  /* ------------------------------------------------------------------ *
   *  SIGNIFICANCE THRESHOLD — either the publication's own value        *
   *  (threshP in the manifest, as-is — this is the default) or a        *
   *  live user-provided one, computed as alpha / (SimpleM's effective   *
   *  marker count, or Bonferroni's raw marker count). Recomputed        *
   *  whenever the dataset, source, method, or alpha changes; drives     *
   *  point coloring, counts, and the plot's threshold line(s) — the     *
   *  published line always stays visible for comparison once the user  *
   *  switches to a custom value.                                       *
   * ------------------------------------------------------------------ */
  let threshSource = 'published'; // 'published' | 'custom'
  let threshMethod = 'simplem';   // 'simplem' | 'bonferroni' (only matters when threshSource==='custom')
  let threshAlpha = 0.05;
  let ACTIVE_THRESH_P = 0;
  let ACTIVE_THRESH_NEGLOG = 0;
  let ACTIVE_TOTAL_SIG = 0;

  /* ------------------------------------------------------------------ *
   *  BULK EXPORT DIALOG STATE — entirely independent of the live plot's *
   *  threshSource/threshMethod/threshAlpha above; always resets to the  *
   *  Published default each time the dialog opens, rather than          *
   *  inheriting whatever the plot currently happens to be showing.      *
   * ------------------------------------------------------------------ */
  let exportOpen = false;
  let exportBusy = false;
  let exportScope = 'filtered';   // 'all' | 'filtered' | 'selected'
  let exportThreshSource = 'published';
  let exportThreshMethod = 'simplem';
  let exportThreshAlpha = 0.05;
  let exportProgress = { phase: '', done: 0, total: 0 };
  let exportSummary = null;       // {clamped:[entries], failed:[entries]} after a run, shown until dialog closes

  /* Pure: p-value + -log10 for a given entry/source/method/alpha combo,
     with no dependency on which dataset happens to be active right now —
     reused by the live plot (against activeEntry) and by bulk export
     (against any candidate entry, cached or not). */
  function computeThresholdFor(entry, source, method, alpha) {
    let p;
    if (source === 'published') {
      p = entry.threshP;
    } else {
      const denom = method === 'bonferroni' ? entry.totalMarkers : entry.effectiveMarkers;
      p = denom ? (alpha / denom) : entry.threshP;
    }
    return { p: p, negLog: p > 0 ? -Math.log10(p) : 99 };
  }

  /* Same, but never looser than the floor every CSV is already filtered
     to — a looser choice would otherwise silently return every row in
     the file, falsely implying completeness for a cutoff that was never
     actually applied upstream. */
  function clampedThreshold(entry, source, method, alpha) {
    const t = computeThresholdFor(entry, source, method, alpha);
    const clamped = t.negLog < FLOOR_NEGLOG;
    return { negLog: clamped ? FLOOR_NEGLOG : t.negLog, clamped: clamped };
  }

  function recomputeThreshold() {
    if (!DATA || !activeEntry) return;
    const t = computeThresholdFor(activeEntry, threshSource, threshMethod, threshAlpha);
    ACTIVE_THRESH_P = t.p;
    ACTIVE_THRESH_NEGLOG = t.negLog;
    let n = 0;
    for (let i = 0; i < DATA.n; i++) if (DATA.negLogP[i] >= ACTIVE_THRESH_NEGLOG) n++;
    ACTIVE_TOTAL_SIG = n;
    if (currentRegion) {
      const sig = [];
      currentRegion.indices.forEach(function (idx) { if (DATA.negLogP[idx] >= ACTIVE_THRESH_NEGLOG) sig.push(idx); });
      currentRegion.sigIndices = sig;
    }
  }

  function publishedThreshNegLog() { return -Math.log10(activeEntry.threshP); }
  function publishedTotalSig() {
    const nl = publishedThreshNegLog();
    let n = 0;
    for (let i = 0; i < DATA.n; i++) if (DATA.negLogP[i] >= nl) n++;
    return n;
  }

  function loadDataset(entry) {
    if (datasetCache.has(entry.id)) {
      DATA = datasetCache.get(entry.id);
      view = { x0: 0, x1: ensureGeometry().TOTAL_CUM };
      threshSource = 'published'; threshMethod = 'simplem'; threshAlpha = 0.05;
      recomputeThreshold();
      rerenderIfActive();
      return;
    }
    datasetPromiseId = entry.id;
    fetch(CFG.dataBase + entry.file).then(function (r) {
      if (!r.ok) throw new Error('HTTP ' + r.status + ' fetching ' + entry.file);
      return r.text();
    }).then(function (raw) {
      ensureGeometry();
      const parsed = parseCsv(raw);
      datasetCache.set(entry.id, parsed);
      if (activeEntry && activeEntry.id === entry.id) {
        DATA = parsed;
        view = { x0: 0, x1: GEOM.TOTAL_CUM };
        threshSource = 'published'; threshMethod = 'simplem'; threshAlpha = 0.05;
        recomputeThreshold();
      }
      rerenderIfActive();
    }).catch(function (err) {
      if (activeEntry && activeEntry.id === entry.id) datasetError = err;
      rerenderIfActive();
    });
  }

  function rerenderIfActive() { if (S.tool === 'snpgwas') render(document.getElementById('page')); }

  function parseCsv(raw) {
    const lines = raw.replace(/\r/g, '').trim().split('\n');
    lines.shift(); // header: Chr,SNP,bp,A1,A2,Freq,b,se,p
    const n = lines.length;

    const chrOf = new Uint8Array(n);
    const bpOf = new Uint32Array(n);
    const negLogP = new Float32Array(n);
    const fields = new Array(n);

    for (let i = 0; i < n; i++) {
      const f = lines[i].split(',');
      fields[i] = f;
      chrOf[i] = parseInt(f[0], 10);
      bpOf[i] = parseInt(f[2], 10);
      const p = parseFloat(f[8]);
      negLogP[i] = p > 0 ? -Math.log10(p) : 99;
    }

    const CHR_COUNT = CFG.chrCount;
    const chrStart = new Int32Array(CHR_COUNT + 1).fill(-1);
    const chrEnd = new Int32Array(CHR_COUNT + 1).fill(-1);
    for (let i = 0; i < n; i++) {
      const c = chrOf[i];
      if (chrStart[c] === -1) chrStart[c] = i;
      chrEnd[c] = i + 1;
    }

    let globalMaxNegLog = 0;
    for (let i = 0; i < n; i++) if (negLogP[i] > globalMaxNegLog) globalMaxNegLog = negLogP[i];
    const Y_MAX = globalMaxNegLog * 1.10;

    return {
      n: n, chrOf: chrOf, bpOf: bpOf, negLogP: negLogP, fields: fields,
      chrStart: chrStart, chrEnd: chrEnd, Y_MAX: Y_MAX,
    };
  }

  /* Side-effect-free fetch+parse+cache for any manifest entry — unlike
     loadDataset(), never touches DATA/activeEntry/view. This is what lets
     bulk export pull in datasets the user never actually visited without
     disturbing whatever's currently on screen. */
  function ensureDatasetCached(entry) {
    if (datasetCache.has(entry.id)) return Promise.resolve(datasetCache.get(entry.id));
    return fetch(CFG.dataBase + entry.file).then(function (r) {
      if (!r.ok) throw new Error('HTTP ' + r.status + ' fetching ' + entry.file);
      return r.text();
    }).then(function (raw) {
      const parsed = parseCsv(raw);
      datasetCache.set(entry.id, parsed);
      return parsed;
    });
  }

  /* binary search helpers — dataset is contiguous per chromosome, sorted by bp */
  function lowerBound(lo, hi, bpOf, target) {
    while (lo < hi) { const mid = (lo + hi) >>> 1; if (bpOf[mid] < target) lo = mid + 1; else hi = mid; }
    return lo;
  }
  function upperBound(lo, hi, bpOf, target) {
    while (lo < hi) { const mid = (lo + hi) >>> 1; if (bpOf[mid] <= target) lo = mid + 1; else hi = mid; }
    return lo;
  }
  function cumAt(idx) { return GEOM.chrOffset[DATA.chrOf[idx]] + DATA.bpOf[idx]; }
  function globalLowerBound(target) {
    let lo = 0, hi = DATA.n;
    while (lo < hi) { const mid = (lo + hi) >>> 1; if (cumAt(mid) < target) lo = mid + 1; else hi = mid; }
    return lo;
  }
  function globalUpperBound(target) {
    let lo = 0, hi = DATA.n;
    while (lo < hi) { const mid = (lo + hi) >>> 1; if (cumAt(mid) <= target) lo = mid + 1; else hi = mid; }
    return lo;
  }
  function chrForCum(x) {
    for (let c = 1; c <= CFG.chrCount; c++) {
      if (x >= GEOM.chrOffset[c] && x <= GEOM.chrOffset[c] + GEOM.chrLen[c]) return c;
    }
    for (let c = 1; c <= CFG.chrCount; c++) if (x < GEOM.chrOffset[c]) return c > 1 ? c - 1 : c;
    return CFG.chrCount;
  }

  function fmtBp(bp) { return Math.round(bp).toLocaleString('en-US'); }
  function fmtMb(bp) { return (bp / 1e6).toFixed(2) + ' Mb'; }
  function fmtP(pStr) { const v = parseFloat(pStr); return isFinite(v) ? v.toExponential(2) : pStr; }
  function escAttr(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }
  function cap1(s) { return s.charAt(0).toUpperCase() + s.slice(1); }

  /* ------------------------------------------------------------------ *
   *  VIEW / INTERACTION STATE — module scope, persists across visits   *
   *  to this tool (so pan/zoom position survives navigating away and   *
   *  back); DOM refs are re-queried fresh every render() since #page   *
   *  is rebuilt from scratch on each visit.                            *
   * ------------------------------------------------------------------ */
  let view = null;           // {x0,x1} in cumulative-bp space
  let mode = 'pan';          // 'pan' | 'select'
  const MARGIN = { top: 50, right: 20, bottom: 34, left: 54 };
  const MIN_SPAN = 300;

  let DOM = {};
  let plotW = 0, plotH = 0, dpr = 1;
  let dragging = false, dragKind = 'pan', dragStartPx = 0, dragStartView = null, dragStartCum = 0;
  let miniDragging = false;
  let miniBins = null;

  let currentRegion = null;  // {indices, sigIndices, chrsTouched, bpStart, bpEnd, chrLabel, coordLabel, spanLabel}
  let regionFilter = 'all';
  const TABLE_CAP = 1500;

  /* table sort — sortCol null means natural (genomic position) order.
     currentTableList is the exact filtered+sorted index list currently on
     screen, kept in sync by renderRegionPanel() so "export shown SNPs"
     always matches what's actually shown. */
  let sortCol = null;
  let sortDir = 'asc';
  let currentTableList = [];

  const SORT_COLS = [
    { key: 'snp',     label: 'SNP',      get: function (idx) { return DATA.fields[idx][1]; },                              type: 'str' },
    { key: 'chr',     label: 'Chr',      get: function (idx) { return DATA.chrOf[idx]; },                                  type: 'num' },
    { key: 'bp',      label: 'BP',       get: function (idx) { return DATA.bpOf[idx]; },                                   type: 'num' },
    { key: 'a1a2',    label: 'A1/A2',    get: function (idx) { return DATA.fields[idx][3] + '/' + DATA.fields[idx][4]; },  type: 'str' },
    { key: 'freq',    label: 'Freq',     get: function (idx) { return parseFloat(DATA.fields[idx][5]); },                  type: 'num' },
    { key: 'beta',    label: 'Beta',     get: function (idx) { return parseFloat(DATA.fields[idx][6]); },                  type: 'num' },
    { key: 'se',      label: 'SE',       get: function (idx) { return parseFloat(DATA.fields[idx][7]); },                  type: 'num' },
    { key: 'p',       label: 'P',        get: function (idx) { return parseFloat(DATA.fields[idx][8]); },                  type: 'num' },
    { key: 'neglogp', label: '−log₁₀P',  get: function (idx) { return DATA.negLogP[idx]; },                                type: 'num' },
  ];

  function sortedIndices(list, key, dir) {
    const col = SORT_COLS.filter(function (c) { return c.key === key; })[0];
    if (!col) return list;
    const mul = dir === 'desc' ? -1 : 1;
    return list.slice().sort(function (a, b) {
      const va = col.get(a), vb = col.get(b);
      return col.type === 'num' ? (va - vb) * mul : String(va).localeCompare(String(vb)) * mul;
    });
  }

  function tableHeadHTML() {
    return '<tr>' + SORT_COLS.map(function (c) {
      const active = sortCol === c.key;
      const arrow = active ? (sortDir === 'asc' ? ' ▲' : ' ▼') : '';
      return '<th class="gwx-sortable' + (active ? ' sorted' : '') + '" data-sort="' + c.key + '">' + c.label + arrow + '</th>';
    }).join('') + '</tr>';
  }

  function clampView(x0, x1) {
    let span = x1 - x0;
    if (span < MIN_SPAN) span = MIN_SPAN;
    if (span > GEOM.TOTAL_CUM) span = GEOM.TOTAL_CUM;
    if (x0 < 0) { x0 = 0; x1 = span; }
    if (x1 > GEOM.TOTAL_CUM) { x1 = GEOM.TOTAL_CUM; x0 = GEOM.TOTAL_CUM - span; }
    return { x0: x0, x1: x1 };
  }
  function setView(x0, x1) {
    const v = clampView(x0, x1);
    view.x0 = v.x0; view.x1 = v.x1;
    updateReadout();
    renderPlot();
    renderMinimapViewport();
  }
  function updateReadout() {
    if (!DOM.readout) return;
    const span = view.x1 - view.x0;
    let label;
    if (view.x0 <= 0 && view.x1 >= GEOM.TOTAL_CUM) {
      label = 'Whole genome (' + fmtMb(GEOM.TOTAL_CUM) + ')';
    } else {
      const startChr = chrForCum(view.x0), endChr = chrForCum(view.x1);
      if (startChr === endChr) {
        const bpStart = Math.max(0, view.x0 - GEOM.chrOffset[startChr]);
        const bpEnd = Math.min(GEOM.chrLen[startChr], view.x1 - GEOM.chrOffset[startChr]);
        label = 'Chr ' + startChr + ': ' + fmtBp(bpStart) + '–' + fmtBp(bpEnd) + ' (' + fmtMb(span) + ')';
      } else {
        label = 'Chr ' + startChr + '–' + endChr + ' (' + fmtMb(span) + ')';
      }
    }
    DOM.readout.innerHTML = label;
  }

  /* ------------------------------------------------------------------ *
   *  CANVAS SIZING                                                      *
   * ------------------------------------------------------------------ */
  function resizeCanvases() {
    if (!DOM.plotCanvas || !document.contains(DOM.plotCanvas)) return;
    dpr = Math.max(1, window.devicePixelRatio || 1);
    const rect = DOM.plotCanvas.getBoundingClientRect();
    plotW = rect.width; plotH = rect.height;
    DOM.plotCanvas.width = Math.round(plotW * dpr);
    DOM.plotCanvas.height = Math.round(plotH * dpr);
    DOM.plotCtx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const mrect = DOM.miniCanvas.getBoundingClientRect();
    DOM.miniCanvas.width = Math.round(mrect.width * dpr);
    DOM.miniCanvas.height = Math.round(mrect.height * dpr);
    DOM.miniCtx.setTransform(dpr, 0, 0, dpr, 0, 0);

    renderPlot();
    buildMinimapBins();
    renderMinimap();
  }

  function xToPx(cumPos) { return MARGIN.left + (cumPos - view.x0) / (view.x1 - view.x0) * (plotW - MARGIN.left - MARGIN.right); }
  function yToPx(nl) { const h = plotH - MARGIN.top - MARGIN.bottom; return MARGIN.top + h - (nl / DATA.Y_MAX) * h; }
  function pxToCum(px) { return view.x0 + (px - MARGIN.left) / (plotW - MARGIN.left - MARGIN.right) * (view.x1 - view.x0); }

  function niceStep(rough) {
    const pow10 = Math.pow(10, Math.floor(Math.log10(rough)));
    const f = rough / pow10;
    const nf = f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10;
    return nf * pow10;
  }

  /* chart-only accent colors (not reused elsewhere on the site, so kept
     as plain constants rather than :root custom properties) */
  function chartColors() {
    const cs = getComputedStyle(document.documentElement);
    function v(name, fallback) { const x = cs.getPropertyValue(name).trim(); return x || fallback; }
    return {
      chrA: v('--navy-700', '#1c3360'), chrB: '#5b8bff',
      sig: '#c0392b', grid: v('--line-2', '#eef2f8'), axis: '#c9d3e4',
      panelBg: v('--white', '#fff'), muted: v('--muted', '#62718a'), ink: v('--ink', '#151b2c'),
      selectFill: 'rgba(207,138,18,.20)', selectBorder: v('--gold', '#cf8a12'),
      viewportFill: 'rgba(37,99,235,.16)', viewportBorder: v('--blue', '#2563eb'),
    };
  }

  /* ------------------------------------------------------------------ *
   *  MAIN PLOT                                                          *
   * ------------------------------------------------------------------ */
  function renderPlot() {
    if (!DOM.plotCanvas) return;
    const ctx = DOM.plotCtx;
    const col = chartColors();
    ctx.clearRect(0, 0, plotW, plotH);
    ctx.fillStyle = col.panelBg;
    ctx.fillRect(0, 0, plotW, plotH);

    const plotLeft = MARGIN.left, plotRight = plotW - MARGIN.right;
    const plotTop = MARGIN.top, plotBottom = plotH - MARGIN.bottom;
    const plotAreaH = plotBottom - plotTop;

    ctx.fillStyle = col.ink || col.muted;
    ctx.font = '600 20px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText(activeEntry.trait + ' - ' + measureShort(activeEntry.measure), (plotLeft + plotRight) / 2, 8);

    const yStep = niceStep(DATA.Y_MAX / 6);
    ctx.strokeStyle = col.grid;
    ctx.lineWidth = 1;
    ctx.fillStyle = col.muted;
    ctx.font = '11px ui-monospace, Menlo, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let yv = 0; yv <= DATA.Y_MAX; yv += yStep) {
      const py = yToPx(yv);
      ctx.beginPath(); ctx.moveTo(plotLeft, py); ctx.lineTo(plotRight, py); ctx.stroke();
      ctx.fillText(yv.toFixed(0), plotLeft - 8, py);
    }
    ctx.save();
    ctx.translate(14, plotTop + plotAreaH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.font = '11px system-ui, sans-serif';
    ctx.fillStyle = col.muted;
    ctx.fillText('−log₁₀(p)', 0, 0);
    ctx.restore();

    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    const visibleChrs = [];
    for (let c = 1; c <= CFG.chrCount; c++) {
      const segX0 = GEOM.chrOffset[c], segX1 = GEOM.chrOffset[c] + GEOM.chrLen[c];
      if (segX1 < view.x0 || segX0 > view.x1) continue;
      visibleChrs.push(c);
      const px0 = Math.max(plotLeft, xToPx(segX0)), px1 = Math.min(plotRight, xToPx(segX1));
      if (px1 - px0 > 40) {
        ctx.fillStyle = col.muted;
        ctx.font = '11px system-ui, sans-serif';
        ctx.fillText('Chr ' + c, (px0 + px1) / 2, plotBottom + 8);
      }
      if (c > 1) {
        const bx = xToPx(segX0);
        if (bx > plotLeft && bx < plotRight) {
          ctx.strokeStyle = col.grid === col.panelBg ? '#e3e9f2' : col.grid;
          ctx.beginPath(); ctx.moveTo(bx, plotTop); ctx.lineTo(bx, plotBottom); ctx.stroke();
        }
      }
    }
    if (visibleChrs.length === 1) {
      const c = visibleChrs[0];
      const bpSpan = Math.min(GEOM.chrLen[c], view.x1 - GEOM.chrOffset[c]) - Math.max(0, view.x0 - GEOM.chrOffset[c]);
      const step = niceStep(bpSpan / 6);
      const startTick = Math.ceil(Math.max(0, view.x0 - GEOM.chrOffset[c]) / step) * step;
      ctx.font = '10.5px ui-monospace, Menlo, monospace';
      ctx.fillStyle = col.muted;
      for (let bp = startTick; bp <= Math.min(GEOM.chrLen[c], view.x1 - GEOM.chrOffset[c]); bp += step) {
        const px = xToPx(GEOM.chrOffset[c] + bp);
        if (px < plotLeft || px > plotRight) continue;
        ctx.strokeStyle = col.grid;
        ctx.beginPath(); ctx.moveTo(px, plotTop); ctx.lineTo(px, plotBottom); ctx.stroke();
        const label = step >= 1e6 ? (bp / 1e6).toFixed(bp === 0 ? 0 : 1) + ' Mb' : (bp / 1e3).toFixed(0) + ' kb';
        ctx.fillText(label, px, plotBottom + 22);
      }
    }

    /* active threshold line always shows; the published one only draws
       separately when it differs from the active value — at the default
       (source==='published') they're the same line, so only one is drawn. */
    function drawThreshLine(negLog, label, color, dash) {
      const y = yToPx(negLog);
      if (y < plotTop || y > plotBottom) return;
      ctx.save();
      ctx.strokeStyle = color;
      ctx.setLineDash(dash);
      ctx.lineWidth = 1.25;
      ctx.beginPath(); ctx.moveTo(plotLeft, y); ctx.lineTo(plotRight, y); ctx.stroke();
      ctx.restore();
      ctx.fillStyle = color;
      ctx.font = '10.5px system-ui, sans-serif';
      ctx.textAlign = 'left'; ctx.textBaseline = 'bottom';
      ctx.fillText(label, plotLeft + 6, y - 3);
    }
    const pubNegLog = publishedThreshNegLog();
    const activeLabel = threshSource === 'published'
      ? ((activeEntry.thresholdMethod || 'Published') + ' threshold (−log₁₀P > ' + ACTIVE_THRESH_NEGLOG.toFixed(2) + ', α=' + activeEntry.thresholdAlpha + ')')
      : ((threshMethod === 'bonferroni' ? 'Bonferroni' : 'SimpleM') + ' threshold (−log₁₀P > ' + ACTIVE_THRESH_NEGLOG.toFixed(2) + ', α=' + threshAlpha + ')');
    drawThreshLine(ACTIVE_THRESH_NEGLOG, activeLabel, col.muted, [5, 4]);
    if (Math.abs(ACTIVE_THRESH_NEGLOG - pubNegLog) > 1e-6) {
      drawThreshLine(pubNegLog, (activeEntry.thresholdMethod || 'Published') + ' threshold (−log₁₀P > ' + pubNegLog.toFixed(2) + ', α=' + activeEntry.thresholdAlpha + ')', col.selectBorder, [2, 3]);
    }

    let pointR = 1.6;
    const spanFrac = (view.x1 - view.x0) / GEOM.TOTAL_CUM;
    if (spanFrac < 0.25) pointR = 2.3;
    if (spanFrac < 0.05) pointR = 3.0;
    if (spanFrac < 0.01) pointR = 3.8;

    for (let c = 1; c <= CFG.chrCount; c++) {
      const segX0 = GEOM.chrOffset[c], segX1 = GEOM.chrOffset[c] + GEOM.chrLen[c];
      if (segX1 < view.x0 || segX0 > view.x1) continue;
      const bpLo = Math.max(0, view.x0 - GEOM.chrOffset[c]);
      const bpHi = Math.min(GEOM.chrLen[c], view.x1 - GEOM.chrOffset[c]);
      const s = DATA.chrStart[c], e = DATA.chrEnd[c];
      if (s === -1) continue;
      const lo = lowerBound(s, e, DATA.bpOf, bpLo), hi = upperBound(s, e, DATA.bpOf, bpHi);
      const baseColor = (c % 2 === 1) ? col.chrA : col.chrB;
      for (let idx = lo; idx < hi; idx++) {
        const px = xToPx(GEOM.chrOffset[c] + DATA.bpOf[idx]);
        const py = yToPx(DATA.negLogP[idx]);
        const isSig = DATA.negLogP[idx] >= ACTIVE_THRESH_NEGLOG;
        ctx.fillStyle = isSig ? col.sig : baseColor;
        ctx.beginPath(); ctx.arc(px, py, isSig ? pointR + 0.4 : pointR, 0, 6.2832); ctx.fill();
      }
    }

    ctx.strokeStyle = '#c9d3e4';
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(plotLeft, plotBottom); ctx.lineTo(plotRight, plotBottom); ctx.stroke();
  }

  /* ------------------------------------------------------------------ *
   *  MINIMAP                                                            *
   * ------------------------------------------------------------------ */
  const MINI_BIN_COUNT = 900;

  function buildMinimapBins() {
    if (!DOM.miniCanvas) return;
    const rect = DOM.miniCanvas.getBoundingClientRect();
    const w = rect.width;
    const binCount = Math.max(200, Math.min(MINI_BIN_COUNT, Math.round(w)));
    const binW = GEOM.TOTAL_CUM / binCount;
    const bins = new Array(binCount);
    for (let b = 0; b < binCount; b++) bins[b] = { maxNL: 0, chr: 0 };
    for (let i = 0; i < DATA.n; i++) {
      const cp = GEOM.chrOffset[DATA.chrOf[i]] + DATA.bpOf[i];
      const b = Math.min(binCount - 1, Math.floor(cp / binW));
      if (DATA.negLogP[i] > bins[b].maxNL) { bins[b].maxNL = DATA.negLogP[i]; bins[b].chr = DATA.chrOf[i]; }
      else if (bins[b].chr === 0) bins[b].chr = DATA.chrOf[i];
    }
    miniBins = { bins: bins, binW: binW };
  }

  function renderMinimapBackgroundOnly() {
    const ctx = DOM.miniCtx;
    const rect = DOM.miniCanvas.getBoundingClientRect();
    const w = rect.width, h = rect.height;
    const col = chartColors();
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = col.grid; ctx.fillRect(0, 0, w, h);
    if (!miniBins) return;
    const bins = miniBins.bins;
    const padBottom = 4, usableH = h - 4 - padBottom;
    const pxPerBin = w / bins.length;
    for (let b = 0; b < bins.length; b++) {
      const bin = bins[b];
      if (bin.chr === 0) continue;
      const barH = Math.max(1, (bin.maxNL / DATA.Y_MAX) * usableH);
      const x = b * pxPerBin;
      const isSig = bin.maxNL >= ACTIVE_THRESH_NEGLOG;
      ctx.fillStyle = isSig ? col.sig : ((bin.chr % 2 === 1) ? col.chrA : col.chrB);
      ctx.fillRect(x, h - padBottom - barH, Math.max(1, pxPerBin), barH);
    }
  }

  function renderMinimap() {
    if (!DOM.miniCanvas) return;
    renderMinimapBackgroundOnly();
    renderMinimapViewport();
  }

  function renderMinimapViewport() {
    if (!DOM.miniCanvas) return;
    const ctx = DOM.miniCtx;
    const rect = DOM.miniCanvas.getBoundingClientRect();
    const w = rect.width, h = rect.height;
    const col = chartColors();
    renderMinimapBackgroundOnly();
    const x0 = (view.x0 / GEOM.TOTAL_CUM) * w;
    const x1 = (view.x1 / GEOM.TOTAL_CUM) * w;
    ctx.fillStyle = col.viewportFill;
    ctx.fillRect(x0, 0, Math.max(1, x1 - x0), h);
    ctx.strokeStyle = col.viewportBorder;
    ctx.lineWidth = 1.5;
    ctx.strokeRect(x0 + 0.75, 0.75, Math.max(1, x1 - x0) - 1.5, h - 1.5);
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x0, -3); ctx.lineTo(x0, h + 3); ctx.moveTo(x1, -3); ctx.lineTo(x1, h + 3); ctx.stroke();
  }

  /* ------------------------------------------------------------------ *
   *  INTERACTION — pan / zoom / select on the main plot                 *
   * ------------------------------------------------------------------ */
  /* Pan/Zoom mode's plain drag pans the view (unchanged) — holding Shift
     while dragging borrows the same rubber-band selection visuals as
     Select-region mode, but zooms into the span on release instead of
     opening the SNP panel. dragKind is fixed for the life of one drag so
     a mode switch or key release mid-drag can't change its behavior. */
  function onPlotMouseDown(e) {
    const rect = DOM.plotCanvas.getBoundingClientRect();
    const px = e.clientX - rect.left;
    if (px < MARGIN.left || px > plotW - MARGIN.right) return;
    dragging = true;
    dragStartPx = px;
    dragStartView = { x0: view.x0, x1: view.x1 };
    dragStartCum = pxToCum(px);
    dragKind = mode === 'select' ? 'region' : (e.shiftKey ? 'zoom' : 'pan');
    if (dragKind === 'pan') {
      DOM.plotCanvas.classList.add('dragging');
    } else {
      hideCrosshair();
      DOM.selrect.style.display = 'block';
      updateSelRectVisual(px, px);
    }
    hideTooltip();
  }

  function onWindowMouseMove(e) {
    if (!DOM.plotCanvas || !document.contains(DOM.plotCanvas)) { dragging = false; return; }
    const rect = DOM.plotCanvas.getBoundingClientRect();
    const px = e.clientX - rect.left, py = e.clientY - rect.top;
    if (dragging) {
      if (dragKind === 'pan') {
        const deltaCum = (px - dragStartPx) / (plotW - MARGIN.left - MARGIN.right) * (dragStartView.x1 - dragStartView.x0);
        setView(dragStartView.x0 - deltaCum, dragStartView.x1 - deltaCum);
      } else {
        updateSelRectVisual(dragStartPx, px);
      }
      return;
    }
    if (mode === 'pan') updateShiftCursor(e.shiftKey);
    const wantsCrosshair = mode === 'select' || (mode === 'pan' && e.shiftKey);
    if (wantsCrosshair) {
      const inPlot = px >= MARGIN.left && px <= plotW - MARGIN.right && py >= 0 && py <= plotH;
      if (inPlot) { showCrosshair(px); showCoordTooltip(px, e.clientX, e.clientY); }
      else { hideCrosshair(); hideTooltip(); }
      return;
    }
    if (px >= 0 && px <= plotW && py >= 0 && py <= plotH) handleHover(px, py, e.clientX, e.clientY);
    else hideTooltip();
  }

  function onWindowMouseUp(e) {
    if (miniDragging) miniDragging = false;
    if (!dragging) return;
    dragging = false;
    if (!DOM.plotCanvas || !document.contains(DOM.plotCanvas)) return;
    DOM.plotCanvas.classList.remove('dragging');
    if (dragKind === 'region' || dragKind === 'zoom') {
      const rect = DOM.plotCanvas.getBoundingClientRect();
      const px = Math.max(MARGIN.left, Math.min(plotW - MARGIN.right, e.clientX - rect.left));
      DOM.selrect.style.display = 'none';
      hideSelLabels();
      const cumA = dragStartCum, cumB = pxToCum(px);
      const selMin = Math.min(cumA, cumB), selMax = Math.max(cumA, cumB);
      if (Math.abs(px - dragStartPx) >= 4) {
        if (dragKind === 'region') openRegionPanel(selMin, selMax);
        else setView(selMin, selMax);
      }
    }
  }

  function onMiniWindowMouseMove(e) { if (miniDragging) jumpFromMinimap(e); }

  function updateSelRectVisual(pxA, pxB) {
    const left = Math.min(pxA, pxB), width = Math.abs(pxB - pxA);
    DOM.selrect.style.left = left + 'px';
    DOM.selrect.style.width = width + 'px';
    updateSelLabels(pxA, pxB);
  }

  function clampLabelPx(px) { return Math.max(30, Math.min(plotW - 30, px)); }

  function showCrosshair(px) {
    if (!DOM.crosshair) return;
    DOM.crosshair.style.left = px + 'px';
    DOM.crosshair.style.display = 'block';
  }
  function hideCrosshair() { if (DOM.crosshair) DOM.crosshair.style.display = 'none'; }

  /* two separate pill labels overlap into unreadable text once the
     selection is narrow on screen — merge into one centered label instead */
  const SEL_LABEL_MERGE_PX = 150;

  function updateSelLabels(pxA, pxB) {
    if (!DOM.selLabelStart) return;
    const pxMin = Math.min(pxA, pxB), pxMax = Math.max(pxA, pxB);
    const cumMin = pxToCum(pxMin), cumMax = pxToCum(pxMax);
    if (pxMax - pxMin < SEL_LABEL_MERGE_PX) {
      DOM.selLabelStart.textContent = combinedCoordLabel(cumMin, cumMax);
      DOM.selLabelStart.style.left = clampLabelPx((pxMin + pxMax) / 2) + 'px';
      DOM.selLabelStart.style.display = 'block';
      DOM.selLabelEnd.style.display = 'none';
    } else {
      DOM.selLabelStart.textContent = coordLabelAt(cumMin);
      DOM.selLabelStart.style.left = clampLabelPx(pxMin) + 'px';
      DOM.selLabelStart.style.display = 'block';
      DOM.selLabelEnd.textContent = coordLabelAt(cumMax);
      DOM.selLabelEnd.style.left = clampLabelPx(pxMax) + 'px';
      DOM.selLabelEnd.style.display = 'block';
    }
  }

  function combinedCoordLabel(cumMin, cumMax) {
    const cA = chrForCum(cumMin), cB = chrForCum(cumMax);
    const bpA = Math.max(0, Math.min(GEOM.chrLen[cA], cumMin - GEOM.chrOffset[cA]));
    const bpB = Math.max(0, Math.min(GEOM.chrLen[cB], cumMax - GEOM.chrOffset[cB]));
    return cA === cB
      ? ('Chr' + cA + ':' + fmtBp(bpA) + '–' + fmtBp(bpB))
      : ('Chr' + cA + ':' + fmtBp(bpA) + ' – Chr' + cB + ':' + fmtBp(bpB));
  }
  function hideSelLabels() {
    if (DOM.selLabelStart) DOM.selLabelStart.style.display = 'none';
    if (DOM.selLabelEnd) DOM.selLabelEnd.style.display = 'none';
  }

  function onPlotWheel(e) {
    e.preventDefault();
    const rect = DOM.plotCanvas.getBoundingClientRect();
    const px = e.clientX - rect.left;
    const cumAtMouse = pxToCum(Math.max(MARGIN.left, Math.min(plotW - MARGIN.right, px)));
    const factor = e.deltaY > 0 ? 1.18 : 1 / 1.18;
    const newSpan = (view.x1 - view.x0) * factor;
    const frac = (cumAtMouse - view.x0) / (view.x1 - view.x0);
    setView(cumAtMouse - frac * newSpan, cumAtMouse - frac * newSpan + newSpan);
  }

  function zoomBy(factor) {
    const mid = (view.x0 + view.x1) / 2;
    const newSpan = (view.x1 - view.x0) * factor;
    setView(mid - newSpan / 2, mid + newSpan / 2);
  }

  function panBy(fraction) {
    const delta = (view.x1 - view.x0) * fraction;
    setView(view.x0 + delta, view.x1 + delta);
  }

  function setMode(m) {
    mode = m;
    if (DOM.modePan) DOM.modePan.classList.toggle('on', m === 'pan');
    if (DOM.modeSelect) DOM.modeSelect.classList.toggle('on', m === 'select');
    if (DOM.plotCanvas) DOM.plotCanvas.classList.toggle('mode-select', m === 'select');
    if (m !== 'select') { hideCrosshair(); hideTooltip(); }
    updateShiftCursor();
  }

  /* Shift-to-zoom cursor feedback in Pan/Zoom mode; MouseEvent.shiftKey on
     a plain mousemove is enough, no dedicated keydown tracking needed. */
  function updateShiftCursor(shiftKey) {
    if (DOM.plotCanvas) DOM.plotCanvas.classList.toggle('shift-zoom', mode === 'pan' && !!shiftKey);
  }

  function jumpFromMinimap(e) {
    const rect = DOM.miniCanvas.getBoundingClientRect();
    const px = e.clientX - rect.left;
    const cp = (px / rect.width) * GEOM.TOTAL_CUM;
    const span = view.x1 - view.x0;
    if (DOM.chrJump) DOM.chrJump.value = '';
    setView(cp - span / 2, cp + span / 2);
  }

  /* ------------------------------------------------------------------ *
   *  HOVER TOOLTIP — binary search over genomic coordinate, correct     *
   *  and fast at every zoom level (incl. whole-genome, 60k points)      *
   * ------------------------------------------------------------------ */
  function handleHover(px, py, clientX, clientY) {
    const cumPerPx = (view.x1 - view.x0) / (plotW - MARGIN.left - MARGIN.right);
    const cumAtMouse = pxToCum(px);
    const windowCum = cumPerPx * 10;
    const lo = globalLowerBound(cumAtMouse - windowCum);
    const hi = globalUpperBound(cumAtMouse + windowCum);
    let best = -1, bestD = 100;
    for (let idx = lo; idx < hi; idx++) {
      const ppx = xToPx(cumAt(idx)), ppy = yToPx(DATA.negLogP[idx]);
      const dx = ppx - px, dy = ppy - py, d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; best = idx; }
    }
    if (best === -1) { hideTooltip(); return; }
    showTooltip(best, clientX, clientY);
  }

  function showTooltip(idx, clientX, clientY) {
    const f = DATA.fields[idx];
    const isSig = DATA.negLogP[idx] >= ACTIVE_THRESH_NEGLOG;
    const wrapRect = DOM.plotCanvas.parentElement.getBoundingClientRect();
    DOM.tooltip.innerHTML =
      '<div class="coord">Chr' + f[0] + ':' + fmtBp(parseInt(f[2], 10)) + '</div>' +
      '<div class="snp">' + f[1] + '</div>' +
      '<div class="row"><span>Alleles</span><b>' + f[3] + '/' + f[4] + '</b></div>' +
      '<div class="row"><span>Freq</span><b>' + f[5] + '</b></div>' +
      '<div class="row"><span>Beta ± SE</span><b>' + f[6] + ' ± ' + f[7] + '</b></div>' +
      '<div class="row"><span>P</span><b>' + fmtP(f[8]) + '</b></div>' +
      (isSig ? '<span class="gwx-sig-badge">significant</span>' : '');
    let left = clientX - wrapRect.left + 14;
    const top = clientY - wrapRect.top + 14;
    if (left + 214 > plotW) left = clientX - wrapRect.left - 224;
    DOM.tooltip.style.left = left + 'px';
    DOM.tooltip.style.top = top + 'px';
    DOM.tooltip.style.display = 'block';
  }
  function hideTooltip() { if (DOM.tooltip) DOM.tooltip.style.display = 'none'; }

  /* coordinate-only readout that follows the cursor in select mode, before
     the user clicks to start dragging out a region */
  function showCoordTooltip(px, clientX, clientY) {
    const wrapRect = DOM.plotCanvas.parentElement.getBoundingClientRect();
    DOM.tooltip.innerHTML = '<div class="coord">' + coordLabelAt(pxToCum(px)) + '</div>';
    let left = clientX - wrapRect.left + 14;
    const top = clientY - wrapRect.top + 14;
    if (left + 214 > plotW) left = clientX - wrapRect.left - 224;
    DOM.tooltip.style.left = left + 'px';
    DOM.tooltip.style.top = top + 'px';
    DOM.tooltip.style.display = 'block';
  }

  /* ------------------------------------------------------------------ *
   *  SNP SEARCH                                                         *
   * ------------------------------------------------------------------ */
  /* JBrowse-style coordinate: "chr6:169202719..169209873" (region) or
     "chr6:169202719" (single position); "chr" prefix optional, commas
     in numbers tolerated. */
  const COORD_RE = /^(?:chr)?(\d{1,2})\s*:\s*([\d,]+)\s*(?:\.\.\s*([\d,]+))?\s*$/i;

  function coordLabelAt(cumPos) {
    const c = chrForCum(cumPos);
    const bp = Math.max(0, Math.min(GEOM.chrLen[c], cumPos - GEOM.chrOffset[c]));
    return 'Chr' + c + ':' + fmtBp(bp);
  }

  /* a real SNP + region drawn from the loaded dataset, so the example is
     always something that will actually work */
  function searchExampleHint() {
    if (!DATA || !DATA.n) return '';
    const f = DATA.fields[Math.floor(DATA.n / 2)];
    const chr = f[0], bp = parseInt(f[2], 10);
    const lo = Math.max(0, bp - 5000), hi = bp + 5000;
    return ' Try a SNP like ' + f[1] + ', or a region like chr' + chr + ':' + lo + '..' + hi + '.';
  }

  let searchStatusTimer = null;
  function searchWarn(msg) {
    if (!DOM.searchStatus) return;
    DOM.searchStatus.textContent = msg;
    DOM.searchStatus.classList.add('show');
    clearTimeout(searchStatusTimer);
    searchStatusTimer = setTimeout(clearSearchStatus, 5000);
  }
  function clearSearchStatus() {
    clearTimeout(searchStatusTimer);
    searchStatusTimer = null;
    if (DOM.searchStatus) { DOM.searchStatus.textContent = ''; DOM.searchStatus.classList.remove('show'); }
  }

  function jumpToRegion(chrNum, bpLo, bpHi) {
    const s = DATA.chrStart[chrNum], e = DATA.chrEnd[chrNum];
    const found = s !== -1 && upperBound(s, e, DATA.bpOf, bpHi) > lowerBound(s, e, DATA.bpOf, bpLo);
    if (DOM.chrJump) DOM.chrJump.value = '';
    if (found) {
      const pad = Math.max((bpHi - bpLo) * 0.08, 200);
      setView(GEOM.chrOffset[chrNum] + Math.max(0, bpLo - pad), GEOM.chrOffset[chrNum] + Math.min(GEOM.chrLen[chrNum], bpHi + pad));
    } else {
      searchWarn('No SNPs in Chr' + chrNum + ':' + fmtBp(bpLo) + '–' + fmtBp(bpHi) + ' — showing the surrounding region.');
      setView(GEOM.chrOffset[chrNum] + Math.max(0, bpLo - 1000), GEOM.chrOffset[chrNum] + Math.min(GEOM.chrLen[chrNum], bpHi + 1000));
    }
  }

  function jumpToPosition(chrNum, bp) {
    const s = DATA.chrStart[chrNum], e = DATA.chrEnd[chrNum];
    let exact = false;
    if (s !== -1) {
      const lo = lowerBound(s, e, DATA.bpOf, bp);
      exact = lo < e && DATA.bpOf[lo] === bp;
    }
    if (DOM.chrJump) DOM.chrJump.value = '';
    if (exact) {
      const halfSpan = Math.max(MIN_SPAN, GEOM.chrLen[chrNum] * 0.01);
      const cp = GEOM.chrOffset[chrNum] + bp;
      setView(cp - halfSpan, cp + halfSpan);
    } else {
      searchWarn('No SNP at Chr' + chrNum + ':' + fmtBp(bp) + ' — showing the surrounding region.');
      setView(GEOM.chrOffset[chrNum] + Math.max(0, bp - 1000), GEOM.chrOffset[chrNum] + Math.min(GEOM.chrLen[chrNum], bp + 1000));
    }
  }

  function onSearchKeydown(e) {
    if (e.key !== 'Enter') return;
    const raw = DOM.snpSearch.value.trim();
    if (!raw) return;
    clearSearchStatus();

    const m = raw.match(COORD_RE);
    if (m) {
      const chrNum = parseInt(m[1], 10);
      if (chrNum < 1 || chrNum > CFG.chrCount) {
        searchWarn('No chromosome ' + chrNum + ' in this dataset (1–' + CFG.chrCount + ').');
        return;
      }
      const bpA = parseInt(m[2].replace(/,/g, ''), 10);
      if (m[3] != null) {
        const bpB = parseInt(m[3].replace(/,/g, ''), 10);
        jumpToRegion(chrNum, Math.min(bpA, bpB), Math.max(bpA, bpB));
      } else {
        jumpToPosition(chrNum, bpA);
      }
      return;
    }

    const qLower = raw.toLowerCase();
    for (let i = 0; i < DATA.n; i++) {
      if (DATA.fields[i][1].toLowerCase() === qLower) {
        const cp = GEOM.chrOffset[DATA.chrOf[i]] + DATA.bpOf[i];
        const halfSpan = Math.max(MIN_SPAN, GEOM.chrLen[DATA.chrOf[i]] * 0.01);
        if (DOM.chrJump) DOM.chrJump.value = '';
        setView(cp - halfSpan, cp + halfSpan);
        return;
      }
    }
    searchWarn('No SNP, position, or region matches “' + raw + '”.' + searchExampleHint());
  }

  /* ------------------------------------------------------------------ *
   *  REGION SELECTION PANEL                                             *
   * ------------------------------------------------------------------ */
  function openRegionPanel(cumMin, cumMax) {
    const indices = [], sigIndices = [], chrsTouched = [];
    for (let c = 1; c <= CFG.chrCount; c++) {
      const segX0 = GEOM.chrOffset[c], segX1 = GEOM.chrOffset[c] + GEOM.chrLen[c];
      if (segX1 < cumMin || segX0 > cumMax) continue;
      const bpLo = Math.max(0, cumMin - GEOM.chrOffset[c]);
      const bpHi = Math.min(GEOM.chrLen[c], cumMax - GEOM.chrOffset[c]);
      const s = DATA.chrStart[c], e = DATA.chrEnd[c];
      if (s === -1) continue;
      const lo = lowerBound(s, e, DATA.bpOf, bpLo), hi = upperBound(s, e, DATA.bpOf, bpHi);
      if (hi > lo) chrsTouched.push(c);
      for (let idx = lo; idx < hi; idx++) {
        indices.push(idx);
        if (DATA.negLogP[idx] >= ACTIVE_THRESH_NEGLOG) sigIndices.push(idx);
      }
    }
    const startChr = chrsTouched[0] || chrForCum(cumMin);
    const endChr = chrsTouched[chrsTouched.length - 1] || chrForCum(cumMax);
    /* cumMin/cumMax are pixel-derived, so the per-chromosome offsets come
       out fractional. Snap to whole base pairs (widen outward so the sent
       region still fully covers what was dragged) before anything stores
       or hands them off — SNPVersity expects integer coordinates. */
    const bpStart = Math.max(0, Math.floor(cumMin - GEOM.chrOffset[startChr]));
    const bpEnd = Math.max(0, Math.ceil(cumMax - GEOM.chrOffset[endChr]));
    const chrLabel = chrsTouched.length <= 1 ? ('Chr ' + startChr) : ('Chr ' + startChr + '–' + endChr);
    const coordLabel = chrsTouched.length <= 1
      ? ('Chr' + startChr + ':' + fmtBp(bpStart) + '–' + fmtBp(bpEnd))
      : ('Chr' + startChr + ':' + fmtBp(bpStart) + '–Chr' + endChr + ':' + fmtBp(bpEnd));

    currentRegion = {
      indices: indices, sigIndices: sigIndices, chrsTouched: chrsTouched,
      chr: startChr, bpStart: bpStart, bpEnd: bpEnd,
      chrLabel: chrLabel, coordLabel: coordLabel, spanLabel: fmtMb(cumMax - cumMin),
    };
    regionFilter = 'all';
    sortCol = null; sortDir = 'asc';
    if (DOM.filterAll) DOM.filterAll.classList.add('on');
    if (DOM.filterSig) DOM.filterSig.classList.remove('on');
    renderRegionPanel();
    if (DOM.panel) DOM.panel.classList.add('open');
    if (DOM.panelBackdrop) DOM.panelBackdrop.classList.add('open');
  }

  function renderRegionPanel() {
    if (!currentRegion || !DOM.regionCoord) return;
    DOM.regionCoord.textContent = currentRegion.coordLabel;
    DOM.regionSub.textContent = currentRegion.chrLabel + ' · span ' + currentRegion.spanLabel;
    DOM.statTotal.textContent = currentRegion.indices.length.toLocaleString();
    DOM.statSig.textContent = currentRegion.sigIndices.length.toLocaleString();

    let list = regionFilter === 'sig' ? currentRegion.sigIndices : currentRegion.indices;
    if (sortCol) list = sortedIndices(list, sortCol, sortDir);
    currentTableList = list;
    if (DOM.tableHead) DOM.tableHead.innerHTML = tableHeadHTML();
    DOM.tableBody.innerHTML = '';
    if (list.length === 0) {
      DOM.emptyState.style.display = 'block';
      DOM.shownHint.textContent = '';
      DOM.exportRegionBtn.disabled = true;
    } else {
      DOM.emptyState.style.display = 'none';
      DOM.exportRegionBtn.disabled = false;
      const cap = Math.min(list.length, TABLE_CAP);
      let rowsHtml = '';
      for (let k = 0; k < cap; k++) {
        const idx = list[k], f = DATA.fields[idx];
        const isSig = DATA.negLogP[idx] >= ACTIVE_THRESH_NEGLOG;
        rowsHtml += '<tr class="' + (isSig ? 'sig' : '') + '">' +
          '<td class="snpid">' + f[1] + '</td>' +
          '<td>' + f[0] + '</td>' +
          '<td>' + fmtBp(parseInt(f[2], 10)) + '</td>' +
          '<td>' + f[3] + '/' + f[4] + '</td>' +
          '<td>' + f[5] + '</td>' +
          '<td>' + f[6] + '</td>' +
          '<td>' + f[7] + '</td>' +
          '<td>' + fmtP(f[8]) + '</td>' +
          '<td>' + DATA.negLogP[idx].toFixed(2) + '</td>' +
          '</tr>';
      }
      DOM.tableBody.innerHTML = rowsHtml;
      DOM.shownHint.textContent = list.length > TABLE_CAP
        ? ('showing first ' + TABLE_CAP.toLocaleString() + ' of ' + list.length.toLocaleString() + ' — export for full list')
        : (list.length.toLocaleString() + ' shown');
    }
  }

  function closePanel() {
    if (DOM.panel) DOM.panel.classList.remove('open');
    if (DOM.panelBackdrop) DOM.panelBackdrop.classList.remove('open');
  }

  /* ------------------------------------------------------------------ *
   *  SEND TO SNPVERSITY — one combined handoff, region and NAM           *
   *  accessions each gated by their own checkbox so either or both can  *
   *  go together. A region has no meaning across chromosomes in         *
   *  SNPVersity's single-S.chr model (checkbox disabled off a multi-    *
   *  chromosome selection); accessions are only offered for NAM         *
   *  datasets, for now. The shared Handoff "Replace current selection"  *
   *  control only ever governs the accession half — a region-only send  *
   *  always adds, since there is no accession list to offer or replace  *
   *  with in that case.                                                 *
   * ------------------------------------------------------------------ */
  const NAM_FOUNDER_PROJECT_RE = /nested association mapping/i;
  const NAM_REFERENCE_PROJECT_RE = /Zm-B73-REFERENCE-NAM/i;

  function defaultVersityDatasetId() {
    return (typeof Data !== 'undefined' && Data.datasets && Data.datasets()[0]) ? Data.datasets()[0].id : '';
  }

  function namAccessionInfo() {
    if (typeof Data === 'undefined' || !Data.datasets || !Data.projectsFor || !Data.accessionsFor) return null;
    const ds = defaultVersityDatasetId();
    if (!ds) return null;
    const projects = Data.projectsFor(ds) || [];
    const wantProjIds = projects
      .filter(function (p) { return NAM_FOUNDER_PROJECT_RE.test(p.title) || NAM_REFERENCE_PROJECT_RE.test(p.title); })
      .map(function (p) { return p.id; });
    if (!wantProjIds.length) return null;
    const ids = (Data.accessionsFor(ds) || [])
      .filter(function (a) { return wantProjIds.indexOf(a.proj) !== -1; })
      .map(function (a) { return a.id; });
    if (!ids.length) return null;
    return { dataset: ds, ids: ids };
  }

  function regionSendable() { return !!currentRegion && currentRegion.chrsTouched.length === 1; }
  function accessionsOfferable() { return !!activeEntry && activeEntry.population === 'NAM' && !!namAccessionInfo(); }

  /* "Send to SNPVersity" opens a small popup with three independent
     checkboxes: region, accessions-replace, accessions-add. The two
     accession checkboxes are mutually exclusive (checking one unchecks
     the other) since they're two mutually-exclusive merge modes of the
     same action, not two different things to send. */
  function openSendPopup() {
    if (!DOM.sendPopup) return;
    const canRegion = regionSendable();
    const canAcc = accessionsOfferable();

    if (DOM.sendRegionChk) { DOM.sendRegionChk.disabled = !canRegion; DOM.sendRegionChk.checked = canRegion; }
    if (DOM.sendAccReplaceChk) { DOM.sendAccReplaceChk.disabled = !canAcc; DOM.sendAccReplaceChk.checked = false; }
    if (DOM.sendAccAddChk) { DOM.sendAccAddChk.disabled = !canAcc; DOM.sendAccAddChk.checked = false; }

    const n = (typeof S !== 'undefined' && S.selected) ? S.selected.size : 0;
    if (DOM.accCountReplace) DOM.accCountReplace.textContent = n;
    if (DOM.accCountAdd) DOM.accCountAdd.textContent = n;

    updateSendPopupState();
    DOM.sendPopup.classList.add('open');
    DOM.sendPopupBackdrop.classList.add('open');
  }

  function closeSendPopup() {
    if (DOM.sendPopup) DOM.sendPopup.classList.remove('open');
    if (DOM.sendPopupBackdrop) DOM.sendPopupBackdrop.classList.remove('open');
  }

  function onAccReplaceChange() {
    if (DOM.sendAccReplaceChk.checked && DOM.sendAccAddChk) DOM.sendAccAddChk.checked = false;
    updateSendPopupState();
  }
  function onAccAddChange() {
    if (DOM.sendAccAddChk.checked && DOM.sendAccReplaceChk) DOM.sendAccReplaceChk.checked = false;
    updateSendPopupState();
  }

  function updateSendPopupState() {
    if (!DOM.sendRegionChk) return;
    const canRegion = regionSendable();
    const wantsRegion = canRegion && DOM.sendRegionChk.checked;
    const wantsReplace = !!DOM.sendAccReplaceChk && DOM.sendAccReplaceChk.checked;
    const wantsAdd = !!DOM.sendAccAddChk && DOM.sendAccAddChk.checked;
    if (DOM.sendConfirmBtn) DOM.sendConfirmBtn.disabled = !wantsRegion && !wantsReplace && !wantsAdd;
    if (DOM.hoHint) {
      DOM.hoHint.textContent = canRegion ? ''
        : 'This region spans multiple chromosomes, so it can’t be sent to SNPVersity; SNPVersity needs a single continuous region, on one chromosome.';
    }
  }

  function confirmSendToVersity() {
    const canRegion = regionSendable();
    const wantsRegion = canRegion && !!DOM.sendRegionChk && DOM.sendRegionChk.checked;
    const wantsReplace = !!DOM.sendAccReplaceChk && DOM.sendAccReplaceChk.checked;
    const wantsAdd = !!DOM.sendAccAddChk && DOM.sendAccAddChk.checked;
    const wantsAcc = (wantsReplace || wantsAdd) && accessionsOfferable();
    if (!wantsRegion && !wantsAcc) return;

    const payload = { from: 'GWAS Explorer' };
    const noteParts = [];

    if (wantsRegion) {
      payload.chr = 'chr' + currentRegion.chr;
      payload.start = currentRegion.bpStart;
      payload.end = currentRegion.bpEnd;
      noteParts.push(activeEntry.trait + ' — ' + activeEntry.measure + ' (' + activeEntry.population + ') — ' + currentRegion.coordLabel);
    }
    if (wantsAcc) {
      const info = namAccessionInfo();
      payload.dataset = info.dataset;
      payload.accessions = info.ids;
      payload.merge = wantsReplace ? 'replace' : 'add';
      noteParts.push('NAM population accessions (B73 v5 reference + NAM founder panel)');
    } else {
      /* No accessions in this payload — never let a leftover merge choice
         apply to (and potentially wipe) the user's existing accession
         selection when we aren't offering any of our own. */
      payload.merge = 'add';
    }
    payload.note = noteParts.join(' · ');

    if (typeof Handoff !== 'undefined' && Handoff.toVersity) Handoff.toVersity(payload);
    else if (typeof window.versityRequest === 'function') window.versityRequest(payload);
  }

  /* ------------------------------------------------------------------ *
   *  CSV EXPORT — Blob + temporary <a download> (no Claude sandbox here)*
   * ------------------------------------------------------------------ */
  const CSV_HEADER = 'Chr,SNP,bp,A1,A2,Freq,b,se,p,neg_log10_p\n';

  /* # -prefixed provenance lines ahead of the header — these files are a
     one-way export (never re-parsed by this app), so a leading comment is
     safe; it just carries the dataset's identity/source out of the app
     with the data, rather than only living on-screen. */
  function csvSourceComment(entry) {
    let c = '# ' + entry.trait + ' - ' + measureShort(entry.measure) + ' · ' + entry.population + ' · ' + entry.publication + '\n';
    if (entry.publicationUrl) c += '# For more information about this dataset, see ' + entry.publicationUrl + '\n';
    return c;
  }

  function buildCsv(entry, data, indices) {
    const out = new Array(indices.length);
    for (let k = 0; k < indices.length; k++) {
      const idx = indices[k], f = data.fields[idx];
      out[k] = f[0] + ',' + f[1] + ',' + f[2] + ',' + f[3] + ',' + f[4] + ',' + f[5] + ',' + f[6] + ',' + f[7] + ',' + f[8] + ',' + data.negLogP[idx].toFixed(4);
    }
    return csvSourceComment(entry) + CSV_HEADER + out.join('\n') + '\n';
  }

  function downloadCsv(filenameBase, csvText) {
    const blob = new Blob([csvText], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = filenameBase + '.csv';
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }

  function exportAllSignificant() {
    const sigAll = [];
    for (let i = 0; i < DATA.n; i++) if (DATA.negLogP[i] >= ACTIVE_THRESH_NEGLOG) sigAll.push(i);
    downloadCsv('GWAS_' + activeEntry.id + '_significant_SNPs_' + sigAll.length, buildCsv(activeEntry, DATA, sigAll));
  }

  function exportRegion() {
    if (!currentRegion) return;
    /* currentTableList is the exact filtered+sorted list currently on
       screen, so the export always matches what's shown. */
    const label = currentRegion.chrLabel.replace(/\s+/g, '').replace(/[–:]/g, '-');
    const suffix = regionFilter === 'sig' ? 'significant' : 'all';
    downloadCsv('GWAS_' + activeEntry.id + '_region_' + label + '_' + suffix, buildCsv(activeEntry, DATA, currentTableList));
  }

  /* ------------------------------------------------------------------ *
   *  BULK EXPORT (Export results…) — zips multiple datasets' significant *
   *  SNPs at once. JSZip is loaded lazily, only when a zip is actually    *
   *  needed, rather than as a static <script> tag — unlike d3 (always    *
   *  needed by PanEffect), this is for one occasional action in one of   *
   *  many tools, not worth adding to every page load for every user.     *
   * ------------------------------------------------------------------ */
  let jszipPromise = null;
  function loadJSZip() {
    if (typeof JSZip !== 'undefined') return Promise.resolve(JSZip);
    if (jszipPromise) return jszipPromise;
    jszipPromise = new Promise(function (resolve, reject) {
      const s = document.createElement('script');
      s.src = 'https://cdn.jsdelivr.net/npm/jszip@3.10.1/dist/jszip.min.js';
      s.onload = function () { resolve(window.JSZip); };
      s.onerror = function () { jszipPromise = null; reject(new Error('Failed to load JSZip from the CDN.')); };
      document.head.appendChild(s);
    });
    return jszipPromise;
  }

  function downloadZip(files, zipName) {
    return loadJSZip().then(function (ZipCtor) {
      const zip = new ZipCtor();
      const used = new Set();
      files.forEach(function (f) {
        let name = exportFilenameFor(f.entry) + '_significant_SNPs.csv', n = 1;
        while (used.has(name)) { n++; name = exportFilenameFor(f.entry) + '_significant_SNPs_' + n + '.csv'; }
        used.add(name);
        zip.file(name, f.csv);
      });
      return zip.generateAsync({ type: 'blob', compression: 'DEFLATE' });
    }).then(function (blob) {
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = zipName + '.zip';
      document.body.appendChild(a); a.click(); document.body.removeChild(a);
      setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
    });
  }

  /* ------------------------------------------------------------------ *
   *  EXPORT DIALOG — "Export results…", opened from the dataset picker  *
   *  card (visible even before any dataset is selected, since "all" /   *
   *  "filtered" scope don't need one). renderExportPopup() regenerates  *
   *  only this popup's own subtree from state — never the whole page —  *
   *  so a bulk export's progress ticks don't disturb the plot, picker,  *
   *  or scroll position while it runs.                                  *
   * ------------------------------------------------------------------ */
  function exportProgressText() {
    if (exportProgress.phase === 'fetching') {
      return 'Fetching ' + exportProgress.done + ' of ' + exportProgress.total +
        ' dataset' + (exportProgress.total === 1 ? '' : 's') + '…';
    }
    if (exportProgress.phase === 'zipping') return 'Building zip file…';
    return '';
  }

  function exportSummaryHTML() {
    if (!exportSummary) return '';
    const parts = [];
    if (exportSummary.clamped.length) {
      parts.push('Threshold capped at p &le; 0.001 (&minus;log&#8321;&#8320;P &ge; 3) for ' +
        exportSummary.clamped.length + ' dataset' + (exportSummary.clamped.length === 1 ? '' : 's') +
        ' — results below that were never available.');
    }
    if (exportSummary.failed.length) {
      parts.push('Could not fetch: ' + exportSummary.failed.map(function (e) {
        return escAttr(e.trait + ' — ' + measureShort(e.measure));
      }).join(', ') + '.');
    }
    if (exportSummary.zipError) parts.push('Zip build failed: ' + escAttr(exportSummary.zipError));
    if (!parts.length) parts.push('Done.');
    return '<div class="gwx-export-note">' + parts.join(' ') + '</div>';
  }

  /* Which entries the current scope choice resolves to — shared by the
     threshold-validity check below and runExport() so they can never
     disagree about what's about to be exported. */
  function resolveExportEntries() {
    return exportScope === 'all' ? MANIFEST.slice()
      : exportScope === 'filtered' ? candidateEntries()
      : (activeEntry ? [activeEntry] : []);
  }

  /* True when the currently configured custom threshold is looser than
     the floor every CSV is already filtered to, for ANY entry the current
     scope would export — i.e. downloading now would silently imply
     completeness for a cutoff that was never actually applied upstream. */
  function exportThresholdInvalid() {
    if (exportThreshSource !== 'custom') return false;
    const entries = resolveExportEntries();
    if (!entries.length) return false;
    return entries.some(function (e) {
      return computeThresholdFor(e, exportThreshSource, exportThreshMethod, exportThreshAlpha).negLog < FLOOR_NEGLOG;
    });
  }

  function exportPopupBodyHTML() {
    const cands = candidateEntries();
    const allCount = MANIFEST.length;
    const canSelected = !!activeEntry;
    const busy = exportBusy;
    const invalid = exportThresholdInvalid();
    return (
      '<div class="gwx-send-popup-head">' +
        '<h3>Export results</h3>' +
        '<button class="gwx-panel-close" id="gwxExportCloseBtn" aria-label="Close"' + (busy ? ' disabled' : '') + '>&times;</button>' +
      '</div>' +
      '<div class="gwx-send-popup-body">' +
        '<div class="gwx-export-scope">' +
          '<label class="gwx-ho-check"><input type="radio" name="gwxExportScope" id="gwxExportScopeAll"' +
            (exportScope === 'all' ? ' checked' : '') + (busy ? ' disabled' : '') + '>' +
            '<span>All available results <em>(' + allCount + ' dataset' + (allCount === 1 ? '' : 's') + ')</em></span></label>' +
          '<label class="gwx-ho-check"><input type="radio" name="gwxExportScope" id="gwxExportScopeFiltered"' +
            (exportScope === 'filtered' ? ' checked' : '') + (busy ? ' disabled' : '') + '>' +
            '<span>Currently filtered <em>(' + cands.length + ' dataset' + (cands.length === 1 ? '' : 's') + ')</em></span></label>' +
          '<label class="gwx-ho-check"><input type="radio" name="gwxExportScope" id="gwxExportScopeSelected"' +
            (exportScope === 'selected' ? ' checked' : '') + (busy || !canSelected ? ' disabled' : '') + '>' +
            '<span>Currently selected' + (canSelected ? ' <em>(' + escAttr(activeEntry.trait) + ' — ' + escAttr(measureShort(activeEntry.measure)) + ')</em>' : ' <em>(pick a dataset first)</em>') + '</span></label>' +
        '</div>' +
        '<div class="gwx-thresh-bar" style="margin-top:0">' +
          '<span class="gwx-thresh-label">Significance threshold: </span>' +
          defPopoverHTML('gwxExpThreshHelpBtn', 'gwxExpThreshHelpPop', 'Significance threshold help', EXPORT_THRESH_DEF_TERMS) +
          '<div class="gwx-seg" role="group" aria-label="Export threshold source">' +
            '<button id="gwxExpThreshPublished" class="' + (exportThreshSource === 'published' ? 'on' : '') + '" type="button"' + (busy ? ' disabled' : '') + '>Published</button>' +
            '<button id="gwxExpThreshCustom" class="' + (exportThreshSource === 'custom' ? 'on' : '') + '" type="button"' + (busy ? ' disabled' : '') + '>Custom</button>' +
          '</div>' +
          '<div class="gwx-thresh-custom" id="gwxExpThreshCustomControls" style="display:' + (exportThreshSource === 'custom' ? 'flex' : 'none') + '">' +
            '<div class="gwx-seg" role="group" aria-label="Export threshold method">' +
              '<button id="gwxExpThreshSimpleM" class="' + (exportThreshMethod === 'simplem' ? 'on' : '') + '" type="button"' + (busy ? ' disabled' : '') + '>SimpleM</button>' +
              '<button id="gwxExpThreshBonferroni" class="' + (exportThreshMethod === 'bonferroni' ? 'on' : '') + '" type="button"' + (busy ? ' disabled' : '') + '>Bonferroni</button>' +
            '</div>' +
            '<label class="gwx-thresh-alpha">&alpha; <input type="text" id="gwxExpAlphaInput" value="' + escAttr(exportThreshAlpha) + '" inputmode="decimal"' + (busy ? ' disabled' : '') + ' /></label>' +
          '</div>' +
        '</div>' +
        // (exportThreshSource === 'custom' ? '<div class="gwx-export-floor-note">Threshold must be at least p &le; 0.001 (&minus;log&#8321;&#8320;P &ge; 3).</div>' : '') +
        (invalid ? '<div class="gwx-export-note gwx-export-warn">This custom threshold is less stringent than the data\'s floor — every file here only keeps p &le; 0.001 (&minus;log&#8321;&#8320;P &ge; 3). Choose <b>Published</b>, or a stricter &alpha;/method, before exporting.</div>' : '') +
        (exportProgress.phase ? '<div class="gwx-export-progress" id="gwxExportProgress">' + exportProgressText() + '</div>' : '') +
        exportSummaryHTML() +
      '</div>' +
      '<div class="gwx-send-popup-foot">' +
        '<button class="btn" id="gwxExportCancelBtn" type="button"' + (busy ? ' disabled' : '') + '>Cancel</button>' +
        '<button class="btn primary" id="gwxExportConfirmBtn" type="button"' + (busy || invalid ? ' disabled' : '') + ' title="' + (invalid ? 'Choose a valid threshold to enable export' : '') + '">' + (busy ? 'Working…' : 'Export') + '</button>' +
      '</div>'
    );
  }

  function exportDialogHTML() {
    return '<div class="gwx-send-popup-backdrop" id="gwxExportBackdrop"></div>' +
      '<div class="gwx-send-popup" id="gwxExportPopup" role="dialog" aria-label="Export results">' + exportPopupBodyHTML() + '</div>';
  }

  function renderExportPopup() {
    const popup = document.getElementById('gwxExportPopup');
    const backdrop = document.getElementById('gwxExportBackdrop');
    if (!popup || !backdrop) return;
    popup.innerHTML = exportPopupBodyHTML();
    popup.classList.toggle('open', exportOpen);
    backdrop.classList.toggle('open', exportOpen);
    wireExportPopupInner();
  }

  /* Wires everything *inside* the popup body — safe to re-run on every
     renderExportPopup() call since innerHTML just replaced those nodes.
     The backdrop and the button that opens the dialog live outside this
     subtree and are wired once from wirePicker() instead. */
  function wireExportPopupInner() {
    const closeBtn = document.getElementById('gwxExportCloseBtn');
    if (!closeBtn) return;
    closeBtn.addEventListener('click', closeExportDialog);
    const cancelBtn = document.getElementById('gwxExportCancelBtn');
    if (cancelBtn) cancelBtn.addEventListener('click', closeExportDialog);
    const confirmBtn = document.getElementById('gwxExportConfirmBtn');
    if (confirmBtn) confirmBtn.addEventListener('click', function () { runExport(); });
    wireHelpPopover('gwxExpThreshHelpBtn', 'gwxExpThreshHelpPop');

    function scopeHandler(scope) { return function () { exportScope = scope; renderExportPopup(); }; }
    const scopeAll = document.getElementById('gwxExportScopeAll');
    const scopeFiltered = document.getElementById('gwxExportScopeFiltered');
    const scopeSelected = document.getElementById('gwxExportScopeSelected');
    if (scopeAll) scopeAll.addEventListener('change', scopeHandler('all'));
    if (scopeFiltered) scopeFiltered.addEventListener('change', scopeHandler('filtered'));
    if (scopeSelected) scopeSelected.addEventListener('change', scopeHandler('selected'));

    const pub = document.getElementById('gwxExpThreshPublished');
    const custom = document.getElementById('gwxExpThreshCustom');
    const simpleM = document.getElementById('gwxExpThreshSimpleM');
    const bonferroni = document.getElementById('gwxExpThreshBonferroni');
    const alphaInput = document.getElementById('gwxExpAlphaInput');
    if (pub) pub.addEventListener('click', function () { exportThreshSource = 'published'; renderExportPopup(); });
    if (custom) custom.addEventListener('click', function () { exportThreshSource = 'custom'; renderExportPopup(); });
    if (simpleM) simpleM.addEventListener('click', function () { exportThreshMethod = 'simplem'; renderExportPopup(); });
    if (bonferroni) bonferroni.addEventListener('click', function () { exportThreshMethod = 'bonferroni'; renderExportPopup(); });
    if (alphaInput) {
      const commit = function () {
        const v = parseFloat(alphaInput.value);
        if (!isFinite(v) || v <= 0 || v >= 1) { alphaInput.value = exportThreshAlpha; return; }
        exportThreshAlpha = v;
        renderExportPopup();
      };
      alphaInput.addEventListener('change', commit);
      alphaInput.addEventListener('keydown', function (e) { if (e.key === 'Enter') { e.preventDefault(); commit(); alphaInput.blur(); } });
    }
  }

  function openExportDialog() {
    exportOpen = true;
    exportBusy = false;
    exportThreshSource = 'published'; exportThreshMethod = 'simplem'; exportThreshAlpha = 0.05;
    exportProgress = { phase: '', done: 0, total: 0 };
    exportSummary = null;
    if (exportScope === 'selected' && !activeEntry) exportScope = 'filtered';
    renderExportPopup();
  }

  function closeExportDialog() {
    if (exportBusy) return;
    exportOpen = false;
    renderExportPopup();
  }

  async function runExport() {
    if (exportBusy || exportThresholdInvalid()) return;
    const entries = resolveExportEntries();
    if (!entries.length) return;

    exportBusy = true;
    exportSummary = null;
    exportProgress = { phase: 'fetching', done: 0, total: entries.length };
    renderExportPopup();

    const failed = [];
    for (let i = 0; i < entries.length; i++) {
      if (!datasetCache.has(entries[i].id)) {
        try { await ensureDatasetCached(entries[i]); }
        catch (err) { failed.push(entries[i]); }
      }
      exportProgress = { phase: 'fetching', done: i + 1, total: entries.length };
      renderExportPopup();
    }

    const clamped = [];
    const files = [];
    entries.forEach(function (e) {
      const data = datasetCache.get(e.id);
      if (!data) return;
      const c = clampedThreshold(e, exportThreshSource, exportThreshMethod, exportThreshAlpha);
      if (c.clamped) clamped.push(e);
      const idx = [];
      for (let i = 0; i < data.n; i++) if (data.negLogP[i] >= c.negLog) idx.push(i);
      files.push({ entry: e, csv: buildCsv(e, data, idx) });
    });

    let zipError = null;
    try {
      if (exportScope === 'selected' && files.length === 1) {
        downloadCsv(exportFilenameFor(files[0].entry) + '_significant_SNPs', files[0].csv);
      } else if (files.length) {
        exportProgress = { phase: 'zipping', done: entries.length, total: entries.length };
        renderExportPopup();
        await downloadZip(files, 'GWAS_export_' + files.length + '_dataset' + (files.length === 1 ? '' : 's'));
      }
    } catch (err) {
      zipError = String(err && err.message || err);
    }

    exportBusy = false;
    exportProgress = { phase: '', done: 0, total: 0 };
    exportSummary = { clamped: clamped, failed: failed, zipError: zipError };
    renderExportPopup();
  }

  /* ------------------------------------------------------------------ *
   *  GLOBAL LISTENERS — attached once; look up "current" elements via  *
   *  the DOM ref object, which is reassigned on every visit to the tool*
   * ------------------------------------------------------------------ */
  let globalWired = false;
  function wireGlobalOnce() {
    if (globalWired) return;
    globalWired = true;
    window.addEventListener('mousemove', onWindowMouseMove);
    window.addEventListener('mousemove', onMiniWindowMouseMove);
    window.addEventListener('mouseup', onWindowMouseUp);
    window.addEventListener('resize', resizeCanvases);
    document.addEventListener('click', function () {
      document.querySelectorAll('.gwx-help-pop.show').forEach(function (p) { p.classList.remove('show'); });
      const traitPop = document.getElementById('gwxTraitDdPop');
      if (traitPop) traitPop.classList.remove('show');
    });
  }

  /* ------------------------------------------------------------------ *
   *  PER-VISIT WIRING (plot UI — only called once a dataset is loaded) *
   * ------------------------------------------------------------------ */
  function wireInteractions() {
    DOM.plotCanvas = document.getElementById('gwxPlot');
    DOM.miniCanvas = document.getElementById('gwxMinimap');
    DOM.selrect = document.getElementById('gwxSelrect');
    DOM.crosshair = document.getElementById('gwxCrosshair');
    DOM.selLabelStart = document.getElementById('gwxSelLabelStart');
    DOM.selLabelEnd = document.getElementById('gwxSelLabelEnd');
    DOM.tooltip = document.getElementById('gwxTooltip');
    DOM.readout = document.getElementById('gwxReadout');
    DOM.chrJump = document.getElementById('gwxChrJump');
    DOM.snpSearch = document.getElementById('gwxSearch');
    DOM.searchStatus = document.getElementById('gwxSearchStatus');
    DOM.modePan = document.getElementById('gwxModePan');
    DOM.modeSelect = document.getElementById('gwxModeSelect');
    DOM.modeHelpBtn = document.getElementById('gwxModeHelpBtn');
    DOM.modeHelpPop = document.getElementById('gwxModeHelpPop');
    DOM.panel = document.getElementById('gwxPanel');
    DOM.panelBackdrop = document.getElementById('gwxPanelBackdrop');
    DOM.regionCoord = document.getElementById('gwxRegionCoord');
    DOM.regionSub = document.getElementById('gwxRegionSub');
    DOM.statTotal = document.getElementById('gwxStatTotal');
    DOM.statSig = document.getElementById('gwxStatSig');
    DOM.tableHead = document.getElementById('gwxTableHead');
    DOM.tableBody = document.getElementById('gwxTableBody');
    DOM.emptyState = document.getElementById('gwxEmpty');
    DOM.shownHint = document.getElementById('gwxShownHint');
    DOM.filterAll = document.getElementById('gwxFilterAll');
    DOM.filterSig = document.getElementById('gwxFilterSig');
    DOM.exportRegionBtn = document.getElementById('gwxExportRegionBtn');
    DOM.openSendBtn = document.getElementById('gwxOpenSendBtn');
    DOM.sendPopup = document.getElementById('gwxSendPopup');
    DOM.sendPopupBackdrop = document.getElementById('gwxSendPopupBackdrop');
    DOM.sendPopupCloseBtn = document.getElementById('gwxSendPopupCloseBtn');
    DOM.sendPopupCancelBtn = document.getElementById('gwxSendPopupCancelBtn');
    DOM.sendConfirmBtn = document.getElementById('gwxSendConfirmBtn');
    DOM.sendRegionChk = document.getElementById('gwxSendRegionChk');
    DOM.sendAccReplaceChk = document.getElementById('gwxSendAccReplaceChk');
    DOM.sendAccAddChk = document.getElementById('gwxSendAccAddChk');
    DOM.accCountReplace = document.getElementById('gwxAccCountReplace');
    DOM.accCountAdd = document.getElementById('gwxAccCountAdd');
    DOM.hoHint = document.getElementById('gwxHoHint');
    DOM.plotCtx = DOM.plotCanvas.getContext('2d');
    DOM.miniCtx = DOM.miniCanvas.getContext('2d');

    setMode(mode);

    DOM.plotCanvas.addEventListener('mousedown', onPlotMouseDown);
    DOM.plotCanvas.addEventListener('wheel', onPlotWheel, { passive: false });
    DOM.plotCanvas.addEventListener('dblclick', function () { setView(0, GEOM.TOTAL_CUM); });
    DOM.miniCanvas.addEventListener('mousedown', function (e) { miniDragging = true; jumpFromMinimap(e); });

    document.getElementById('gwxZoomInBtn').addEventListener('click', function () { zoomBy(1 / 1.6); });
    document.getElementById('gwxZoomOutBtn').addEventListener('click', function () { zoomBy(1.6); });
    document.getElementById('gwxPanLeftBtn').addEventListener('click', function () { panBy(-0.4); });
    document.getElementById('gwxPanRightBtn').addEventListener('click', function () { panBy(0.4); });
    document.getElementById('gwxResetBtn').addEventListener('click', function () {
      DOM.chrJump.value = '';
      setView(0, GEOM.TOTAL_CUM);
      resetThresholdSource();
    });
    DOM.chrJump.addEventListener('change', function () {
      const c = parseInt(DOM.chrJump.value, 10);
      if (!c) { setView(0, GEOM.TOTAL_CUM); return; }
      setView(GEOM.chrOffset[c], GEOM.chrOffset[c] + GEOM.chrLen[c]);
    });
    DOM.modePan.addEventListener('click', function () { setMode('pan'); });
    DOM.modeSelect.addEventListener('click', function () { setMode('select'); });
    if (DOM.modeHelpBtn) DOM.modeHelpBtn.addEventListener('click', function (e) {
      e.stopPropagation();
      DOM.modeHelpPop.classList.toggle('show');
    });
    DOM.snpSearch.addEventListener('keydown', onSearchKeydown);
    DOM.snpSearch.addEventListener('input', clearSearchStatus);

    DOM.filterAll.addEventListener('click', function () {
      regionFilter = 'all'; DOM.filterAll.classList.add('on'); DOM.filterSig.classList.remove('on'); renderRegionPanel();
    });
    DOM.filterSig.addEventListener('click', function () {
      regionFilter = 'sig'; DOM.filterSig.classList.add('on'); DOM.filterAll.classList.remove('on'); renderRegionPanel();
    });
    DOM.tableHead.addEventListener('click', function (e) {
      const th = e.target.closest('th[data-sort]');
      if (!th) return;
      const key = th.getAttribute('data-sort');
      if (sortCol === key) sortDir = sortDir === 'asc' ? 'desc' : 'asc';
      else { sortCol = key; sortDir = 'asc'; }
      renderRegionPanel();
    });
    document.getElementById('gwxPanelCloseBtn').addEventListener('click', closePanel);
    DOM.panelBackdrop.addEventListener('click', closePanel);
    DOM.exportRegionBtn.addEventListener('click', exportRegion);
    DOM.openSendBtn.addEventListener('click', openSendPopup);
    DOM.sendPopupCloseBtn.addEventListener('click', closeSendPopup);
    DOM.sendPopupCancelBtn.addEventListener('click', closeSendPopup);
    DOM.sendPopupBackdrop.addEventListener('click', closeSendPopup);
    DOM.sendRegionChk.addEventListener('change', updateSendPopupState);
    if (DOM.sendAccReplaceChk) DOM.sendAccReplaceChk.addEventListener('change', onAccReplaceChange);
    if (DOM.sendAccAddChk) DOM.sendAccAddChk.addEventListener('change', onAccAddChange);
    DOM.sendConfirmBtn.addEventListener('click', confirmSendToVersity);
    document.getElementById('gwxExportAllBtn').addEventListener('click', exportAllSignificant);

    THRESH_BAR_SUFFIXES.forEach(wireThresholdBar);
    wireHelpPopover('gwxThreshHelpBtn', 'gwxThreshHelpPop');
    wireHelpPopover('gwxThreshHelpBtnPanel', 'gwxThreshHelpPopPanel');
    wireHelpPopover('gwxRegionHelpBtn', 'gwxRegionHelpPop');
    wireHelpPopover('gwxIntroHelpBtn', 'gwxIntroHelpPop');
    wireGlobalOnce();
  }

  /* ------------------------------------------------------------------ *
   *  THRESHOLD BAR — Published (the manifest's stored value, as-is) vs  *
   *  Custom (a live alpha / method calculation). Only re-renders the    *
   *  plot + minimap + readout (and the region panel, if one is open)    *
   *  rather than the whole page, so pan/zoom/region state isn't lost.   *
   *                                                                      *
   *  Rendered twice — once above the plot, once inside the region panel *
   *  (see thresholdBarHTML's suffix) — because the panel is a fixed,    *
   *  full-viewport modal that covers the plot copy while a region is    *
   *  open, which is exactly when a user most wants to change it. Both   *
   *  copies drive the same threshSource/threshMethod/threshAlpha state  *
   *  and are kept visually in sync by syncThresholdBarUI(), regardless  *
   *  of which copy was actually clicked.                                *
   * ------------------------------------------------------------------ */
  const THRESH_BAR_SUFFIXES = ['', 'Panel'];

  function syncThresholdBarUI() {
    const hasSimpleM = !!(activeEntry && activeEntry.effectiveMarkers);
    const hasBonferroni = !!(activeEntry && activeEntry.totalMarkers);
    THRESH_BAR_SUFFIXES.forEach(function (suffix) {
      const published = document.getElementById('gwxThreshPublished' + suffix);
      if (!published) return;
      const custom = document.getElementById('gwxThreshCustom' + suffix);
      const customControls = document.getElementById('gwxThreshCustomControls' + suffix);
      const simpleM = document.getElementById('gwxThreshSimpleM' + suffix);
      const bonferroni = document.getElementById('gwxThreshBonferroni' + suffix);
      const alphaInput = document.getElementById('gwxAlphaInput' + suffix);
      const readout = document.getElementById('gwxThreshReadout' + suffix);
      published.classList.toggle('on', threshSource === 'published');
      custom.classList.toggle('on', threshSource === 'custom');
      customControls.style.display = threshSource === 'custom' ? 'flex' : 'none';
      simpleM.classList.toggle('on', threshMethod === 'simplem');
      simpleM.disabled = !hasSimpleM;
      bonferroni.classList.toggle('on', threshMethod === 'bonferroni');
      bonferroni.disabled = !hasBonferroni;
      if (document.activeElement !== alphaInput) alphaInput.value = threshAlpha;
      if (readout) readout.innerHTML = threshReadoutText();
    });
  }

  function wireThresholdBar(suffix) {
    suffix = suffix || '';
    const published = document.getElementById('gwxThreshPublished' + suffix);
    const custom = document.getElementById('gwxThreshCustom' + suffix);
    const simpleM = document.getElementById('gwxThreshSimpleM' + suffix);
    const bonferroni = document.getElementById('gwxThreshBonferroni' + suffix);
    const alphaInput = document.getElementById('gwxAlphaInput' + suffix);
    if (!published) return;

    published.addEventListener('click', function () { threshSource = 'published'; onThresholdChange(); });
    custom.addEventListener('click', function () { threshSource = 'custom'; onThresholdChange(); });
    simpleM.addEventListener('click', function () { if (simpleM.disabled) return; threshMethod = 'simplem'; onThresholdChange(); });
    bonferroni.addEventListener('click', function () { if (bonferroni.disabled) return; threshMethod = 'bonferroni'; onThresholdChange(); });

    function commitAlpha() {
      const v = parseFloat(alphaInput.value);
      if (!isFinite(v) || v <= 0 || v >= 1) { alphaInput.value = threshAlpha; return; }
      threshAlpha = v;
      onThresholdChange();
    }
    alphaInput.addEventListener('change', commitAlpha);
    alphaInput.addEventListener('keydown', function (e) { if (e.key === 'Enter') { e.preventDefault(); commitAlpha(); alphaInput.blur(); } });
  }

  function resetThresholdSource() {
    if (threshSource === 'published' && threshMethod === 'simplem' && threshAlpha === 0.05) return;
    threshSource = 'published'; threshMethod = 'simplem'; threshAlpha = 0.05;
    onThresholdChange();
  }

  function onThresholdChange() {
    recomputeThreshold();
    renderPlot();
    renderMinimap();
    syncThresholdBarUI();
    if (currentRegion) renderRegionPanel();
  }

  /* ------------------------------------------------------------------ *
   *  DATASET PICKER — faceted trait/population/publication filter      *
   * ------------------------------------------------------------------ */
  const TRAIT_CATEGORY_ORDER = ['Plant Architecture', 'Yield Component', 'Flowering Time'];

  function facetFieldHTML(facet, label) {
    const opts = optionsFor(facet);
    const cur = selection[facet];
    let optsHtml = '<option value="">All ' + label.toLowerCase() + 's</option>';
    opts.forEach(function (o) {
      optsHtml += '<option value="' + escAttr(o.value) + '"' + (cur === o.value ? ' selected' : '') + '>' +
        escAttr(o.value) + (o.count > 1 ? ' (' + o.count + ')' : '') + '</option>';
    });
    return '<div class="field gwx-picker-field"><label>' + label + '</label>' +
      '<select id="gwxFacet' + cap1(facet) + '" data-facet="' + facet + '">' + optsHtml + '</select></div>';
  }

  /* Trait picker is a custom dropdown (not a native <select>) so that the
     category headers can be clicked directly to browse/select a whole
     category — a native <optgroup> label isn't an interactive element. */
  function traitDropdownGroups() {
    function inScope(e) {
      return (selection.population == null || e.population === selection.population) &&
        (selection.publication == null || e.publication === selection.publication);
    }
    const counts = new Map(), catOf = new Map();
    MANIFEST.forEach(function (e) {
      if (!inScope(e)) return;
      counts.set(e.trait, (counts.get(e.trait) || 0) + 1);
      catOf.set(e.trait, e.traitCategory);
    });
    const groups = new Map();
    TRAIT_CATEGORY_ORDER.forEach(function (c) { groups.set(c, []); });
    Array.from(counts.keys()).sort(function (a, b) { return a.localeCompare(b); }).forEach(function (t) {
      const cat = catOf.get(t) || 'Other';
      if (!groups.has(cat)) groups.set(cat, []);
      groups.get(cat).push({ value: t, count: counts.get(t) });
    });
    return groups;
  }

  function traitDropdownHTML() {
    const groups = traitDropdownGroups();
    const label = selection.trait || selection.traitCategory || 'All traits';
    let itemsHtml = '<button class="gwx-trait-dd-opt' + (!selection.trait && !selection.traitCategory ? ' on' : '') +
      '" data-kind="all" type="button">All traits</button>';
    groups.forEach(function (list, cat) {
      if (!list.length) return;
      const catTotal = list.reduce(function (s, o) { return s + o.count; }, 0);
      itemsHtml += '<button class="gwx-trait-dd-cat' + (selection.traitCategory === cat && !selection.trait ? ' on' : '') +
        '" data-kind="category" data-value="' + escAttr(cat) + '" type="button">' + escAttr(cat) +
        ' <span class="gwx-trait-dd-count">(' + catTotal + ')</span></button>';
      list.forEach(function (o) {
        itemsHtml += '<button class="gwx-trait-dd-opt' + (selection.trait === o.value ? ' on' : '') +
          '" data-kind="trait" data-value="' + escAttr(o.value) + '" type="button">' + escAttr(o.value) +
          (o.count > 1 ? ' <span class="gwx-trait-dd-count">(' + o.count + ')</span>' : '') + '</button>';
      });
    });
    return '<div class="field gwx-picker-field gwx-trait-dd">' +
      '<label>Trait</label>' +
      '<button class="gwx-trait-dd-btn" id="gwxTraitDdBtn" type="button">' + escAttr(label) + '</button>' +
      '<div class="gwx-trait-dd-pop" id="gwxTraitDdPop">' + itemsHtml + '</div>' +
    '</div>';
  }

  function measureShort(measure) {
    if (measure === 'Trait value (intercept)') return 'Trait';
    if (measure === 'Linear plasticity (slope)') return 'Plasticity';
    return measure;
  }

  /* human-readable, filesystem-safe identifier for a manifest entry —
     used so exported filenames are identifiable outside the app without
     needing to decode entry.id */
  function slugify(s) {
    return String(s == null ? '' : s).trim().replace(/[^\w\s-]/g, '').replace(/\s+/g, '-');
  }
  function exportFilenameFor(entry) {
    return slugify(entry.trait) + '_' + slugify(measureShort(entry.measure)) + '_' +
      slugify(entry.population) + '_' + slugify(entry.publication);
  }

  function chipListHTML(list, activeId) {
    return '<div class="gwx-picker-chips">' + list.map(function (e) {
      return '<button class="gwx-chip' + (e.id === activeId ? ' on' : '') + '" data-entry="' + escAttr(e.id) + '" type="button">' +
        escAttr(e.trait) + ' - ' + escAttr(measureShort(e.measure)) +
        '<span class="gwx-chip-sub">' + escAttr(e.population) + ' · ' + escAttr(e.publication) + '</span></button>';
    }).join('') + '</div>';
  }

  const CHIP_CAP = 20;

  /* The shared publicationUrl across a set of entries, or null if they
     don't all cite the same publication (or none have a URL at all). */
  function sharedPublicationUrl(entries) {
    if (!entries.length) return null;
    const url = entries[0].publicationUrl;
    if (!url) return null;
    for (let i = 1; i < entries.length; i++) {
      if (entries[i].publication !== entries[0].publication || entries[i].publicationUrl !== url) return null;
    }
    return url;
  }

  function doiNoteHTML(url) {
    return url ? '<div class="gwx-picker-doi">For more information about this publication, see <a href="' + escAttr(url) + '" target="_blank" rel="noopener">' + escAttr(url) + '</a></div>' : '';
  }

  function pickerHTML() {
    const cands = candidateEntries();
    let status;
    if (activeEntry) {
      const alt = cands.filter(function (e) { return e.id !== activeEntry.id; });
      status = '<div class="gwx-picker-active">Showing <b>' + escAttr(activeEntry.trait) + '</b> — ' + escAttr(activeEntry.measure) + ' · ' +
        escAttr(activeEntry.population) + ' · ' + escAttr(activeEntry.publication) +
        (alt.length ? ' · ' + alt.length + ' other match' + (alt.length === 1 ? '' : 'es') + ' for this filter' : '') +
        ' <button class="gwx-picker-clear" id="gwxPickerClear" type="button">change dataset</button>' +
        (activeEntry.publicationUrl ? doiNoteHTML(activeEntry.publicationUrl) : '') +
        '</div>' +
        (alt.length && cands.length <= CHIP_CAP ? chipListHTML(cands, activeEntry.id) : '');
    } else if (cands.length === 0) {
      status = '<div class="gwx-picker-empty">No dataset matches this combination. ' +
        '<button class="gwx-picker-clear" id="gwxPickerClear" type="button">Clear filters</button></div>';
    } else {
      status = (selection.publication ? doiNoteHTML(sharedPublicationUrl(cands)) : '') +
        '<div class="gwx-picker-empty">' + cands.length + ' datasets match — narrow the filters above' +
        (cands.length <= CHIP_CAP ? ', or pick one directly:' : '.') + '</div>' +
        (cands.length <= CHIP_CAP ? chipListHTML(cands, null) : '');
    }
    return '<div class="card pad gwx-picker">' +
      '<div class="gwx-picker-row">' +
        facetFieldHTML('publication', 'Publication') + traitDropdownHTML() + facetFieldHTML('population', 'Population') +
      '</div>' +
      status +
      '<div class="gwx-picker-export-row">' +
        '<button class="btn gwx-btn-sm" id="gwxOpenExportBtn" type="button">' + (typeof ICONS !== 'undefined' && ICONS.download ? ICONS.download : '') + ' Export results…</button>' +
      '</div>' +
    '</div>' +
    exportDialogHTML();
  }

  function wirePicker() {
    FACETS.forEach(function (facet) {
      const el = document.getElementById('gwxFacet' + cap1(facet));
      if (el) el.addEventListener('change', function () { onFacetChange(facet, el.value); });
    });
    const traitDdBtn = document.getElementById('gwxTraitDdBtn');
    const traitDdPop = document.getElementById('gwxTraitDdPop');
    if (traitDdBtn && traitDdPop) {
      traitDdBtn.addEventListener('click', function (e) {
        e.stopPropagation();
        traitDdPop.classList.toggle('show');
      });
      traitDdPop.querySelectorAll('.gwx-trait-dd-opt,.gwx-trait-dd-cat').forEach(function (btn) {
        btn.addEventListener('click', function (e) {
          e.stopPropagation();
          const kind = btn.getAttribute('data-kind');
          if (kind === 'all') { selection.trait = null; selection.traitCategory = null; }
          else if (kind === 'category') { selection.trait = null; selection.traitCategory = btn.getAttribute('data-value'); }
          else { selection.trait = btn.getAttribute('data-value'); selection.traitCategory = null; }
          resolveSelection();
          render(document.getElementById('page'));
        });
      });
    }
    const clearBtn = document.getElementById('gwxPickerClear');
    if (clearBtn) clearBtn.addEventListener('click', function () {
      selection = { trait: null, traitCategory: null, measure: null, population: null, publication: null };
      activeEntry = null; DATA = null; datasetError = null; currentRegion = null;
      render(document.getElementById('page'));
    });
    document.querySelectorAll('.gwx-chip').forEach(function (btn) {
      btn.addEventListener('click', function () {
        const entry = MANIFEST.find(function (e) { return e.id === btn.getAttribute('data-entry'); });
        if (!entry) return;
        setActiveEntry(entry);
        render(document.getElementById('page'));
      });
    });
    const retry = document.getElementById('gwxDatasetRetry');
    if (retry) retry.addEventListener('click', function () {
      datasetError = null;
      if (activeEntry) loadDataset(activeEntry);
      render(document.getElementById('page'));
    });

    const openExportBtn = document.getElementById('gwxOpenExportBtn');
    if (openExportBtn) openExportBtn.addEventListener('click', openExportDialog);
    const exportBackdrop = document.getElementById('gwxExportBackdrop');
    if (exportBackdrop) exportBackdrop.addEventListener('click', closeExportDialog);
    renderExportPopup();
  }

  /* ------------------------------------------------------------------ *
   *  MARKUP                                                             *
   * ------------------------------------------------------------------ */
  function manifestLoadingHTML() {
    return '<div class="sec"><div class="bar"></div><div><h2 style="font-size:16px">Loading available datasets…</h2>' +
      '<p>Fetching ' + escAttr(CFG.manifestUrl) + '.</p></div></div>';
  }
  function manifestErrorHTML() {
    return '<div class="card pad" style="border-color:#f3c2bd;background:#fff6f5;margin-top:14px">' +
      '<h3 style="font-family:var(--disp);color:#b42318;margin:0 0 8px">Could not load the dataset list</h3>' +
      '<p style="margin:0 0 12px;color:var(--muted);font-size:13px">' + escAttr(String(manifestError && manifestError.message || manifestError)) + '</p>' +
      '<button class="btn" id="gwxManifestRetry">Retry</button></div>';
  }
  function datasetErrorHTML() {
    return '<div class="card pad" style="border-color:#f3c2bd;background:#fff6f5;margin-top:14px">' +
      '<h3 style="font-family:var(--disp);color:#b42318;margin:0 0 8px">Dataset failed to load</h3>' +
      '<p style="margin:0 0 12px;color:var(--muted);font-size:13px">' + escAttr(String(datasetError && datasetError.message || datasetError)) + '</p>' +
      '<button class="btn" id="gwxDatasetRetry">Retry</button></div>';
  }
  function datasetLoadingHTML() {
    return '<div class="sec"><div class="bar"></div><div><h2 style="font-size:16px">Loading ' + escAttr(activeEntry.trait) + '…</h2>' +
      '<p>Fetching ' + escAttr(activeEntry.file) + '.</p></div></div>';
  }

  function formatPValue(p) {
    return p.toExponential(2).replace(/e([+-])(\d)$/, 'e$10$2');
  }

  /* Always-visible general header — shown in every state (manifest
     loading/error, picker, dataset selected), same style as SNPVersity's
     own top header. Explains the tool once; the per-dataset header below
     (introHTML) only adds the facts specific to whatever's selected, so
     the two don't repeat the same framing. */
  function generalHeaderHTML() {
    return '<div class="sec"><div class="bar"></div><div style="width:100%">' +
      '<div class="n">GWAS MANHATTAN PLOT EXPLORER </div>' +
      '<h2>Explore published GWAS results</h2>' +
      '<p>Browse Manhattan plots displaying results from curated GWAS publications across traits and populations. ' +
      'Pick a dataset below, then pan, zoom, and drag-select regions to inspect significant SNPs; ' +
      'export them as a CSV, or hand a region and/or population off to SNPVersity to keep exploring.</p>' +
    '</div></div>';
  }

  function introHTML() {
    const totalMarkersStr = activeEntry.totalMarkers != null ? activeEntry.totalMarkers.toLocaleString() : DATA.n.toLocaleString();
    return '<div class="sec"><div class="bar"></div><div style="width:100%">' +
      '<div style="display:flex;align-items:center;gap:8px"><h2 style="font-size:16px;margin:0">Dataset details</h2>' +
        defPopoverHTML('gwxIntroHelpBtn', 'gwxIntroHelpPop', 'Dataset details definitions', INTRO_DEF_TERMS) +
      '</div>' +
      '<p style="max-width:none;margin-top:6px">Trait: <b>' + escAttr(activeEntry.trait) + ' - ' + escAttr(measureShort(activeEntry.measure)) + '</b> · Population: <b>' + escAttr(activeEntry.population) + '</b> · Publication: <b>' +
      (activeEntry.publicationUrl ? '<a href="' + escAttr(activeEntry.publicationUrl) + '" target="_blank" rel="noopener">' + escAttr(activeEntry.publication) + '</a>' : escAttr(activeEntry.publication)) +
      '</b> · GWAS method: <b>' + escAttr(activeEntry.gwasMethod) + '</b><br>' +
      'Total markers: <b>' + totalMarkersStr + '</b> · Published significance threshold: <b>p &lt; ' + formatPValue(activeEntry.threshP) + '</b> (−log₁₀P &gt; ' + publishedThreshNegLog().toFixed(2) + ') (' + escAttr(activeEntry.thresholdMethod || 'SimpleM') + ', α=' + escAttr(activeEntry.thresholdAlpha) + ')· ' +
      'Significant markers: <b>' + publishedTotalSig().toLocaleString() + '</b>' +
      (activeEntry.publicationUrl ? '<br>For more information about this dataset, see <a href="' + escAttr(activeEntry.publicationUrl) + '" target="_blank" rel="noopener">' + escAttr(activeEntry.publicationUrl) + '</a>' : '') +
      '</p>' +
    '</div></div>';
  }

  function threshReadoutText() {
    const label = threshSource === 'published'
      ? ((activeEntry.thresholdMethod || 'SimpleM') + ', α=' + activeEntry.thresholdAlpha + ' (published)')
      : (threshMethod === 'bonferroni' ? 'Bonferroni' : 'SimpleM') + ', α=' + threshAlpha;
    return 'p &lt; ' + formatPValue(ACTIVE_THRESH_P) + ' (−log₁₀P &gt; ' + ACTIVE_THRESH_NEGLOG.toFixed(2) + '), ' + label + ' — <b>' +
      ACTIVE_TOTAL_SIG.toLocaleString() + '</b> significant';
  }

  /* ------------------------------------------------------------------ *
   *  DEFINITIONS POPOVERS — reuse the Help & FAQ page's own "GWAS       *
   *  Explorer" glossary (SNPHelp.definitionsFor) rather than duplicate  *
   *  the wording here, so editing snphelp.js is the one place that     *
   *  keeps every "?" button in sync. Degrades to nothing (no button)    *
   *  if SNPHelp hasn't loaded or none of the requested terms exist.     *
   * ------------------------------------------------------------------ */
  function gwasDefTerms(names) {
    if (typeof SNPHelp === 'undefined' || typeof SNPHelp.definitionsFor !== 'function') return [];
    const items = SNPHelp.definitionsFor('GWAS Explorer') || [];
    const byName = {};
    items.forEach(function (pair) { byName[pair[0]] = pair[1]; });
    return names.filter(function (n) { return byName[n]; }).map(function (n) { return [n, byName[n]]; });
  }

  /* Turns a bare URL inside already-escaped text into a clickable link,
     trimming trailing sentence punctuation (periods, commas, closing
     parens, ...) out of the link itself. */
  function linkify(escaped) {
    return escaped.replace(/https?:\/\/[^\s<]+/g, function (url) {
      let trail = '';
      while (/[.,;:)]$/.test(url)) { trail = url.slice(-1) + trail; url = url.slice(0, -1); }
      return '<a href="' + url + '" target="_blank" rel="noopener">' + url + '</a>' + trail;
    });
  }

  function defPopoverHTML(btnId, popId, ariaLabel, names) {
    const terms = gwasDefTerms(names);
    if (!terms.length) return '';
    const body = terms.map(function (t) { return '<b>' + escAttr(t[0]) + '</b>: ' + linkify(escAttr(t[1])); }).join('<br>');
    return '<div class="gwx-help-wrap">' +
      '<button class="gwx-help-btn" id="' + btnId + '" type="button" aria-label="' + escAttr(ariaLabel) + '">?</button>' +
      '<div class="gwx-help-pop" id="' + popId + '">' + body + '</div>' +
    '</div>';
  }

  function wireHelpPopover(btnId, popId) {
    const btn = document.getElementById(btnId);
    const pop = document.getElementById(popId);
    if (!btn || !pop) return;
    btn.addEventListener('click', function (e) {
      e.stopPropagation();
      document.querySelectorAll('.gwx-help-pop.show').forEach(function (p) { if (p !== pop) p.classList.remove('show'); });
      pop.classList.toggle('show');
    });
  }

  const THRESH_DEF_TERMS = ['Published significance threshold', 'Custom threshold', 'SimpleM', 'Bonferroni', 'α (alpha)', '−log₁₀P', 'Significant markers'];
  const REGION_DEF_TERMS = ['SNP', 'Chr', 'BP', 'A1/A2', 'Freq', 'Beta', 'SE', 'P', '−log₁₀P'];
  const INTRO_DEF_TERMS = ['Trait', 'Measure: Trait value (intercept)', 'Measure: Linear plasticity (slope)', 'GWAS method', 'Total markers', 'Significant markers', 'Published significance threshold'];
  const EXPORT_THRESH_DEF_TERMS = ['Published significance threshold', 'Custom threshold', 'SimpleM', 'Bonferroni', 'α (alpha)'];

  /* suffix lets the exact same control render twice with non-colliding
     ids — once above the plot, once inside the region panel (see
     THRESH_BAR_SUFFIXES) — because the panel is a fixed, full-viewport
     modal that sits on top of the plot copy and swallows its clicks
     while a region is selected, which is precisely when a user most
     wants to try a different threshold against the table. */
  function thresholdBarHTML(suffix) {
    suffix = suffix || '';
    const hasSimpleM = !!activeEntry.effectiveMarkers;
    const hasBonferroni = !!activeEntry.totalMarkers;
    return '<div class="card gwx-thresh-bar">' +
      '<span class="gwx-thresh-label">Display significance threshold: </span>' +
      defPopoverHTML('gwxThreshHelpBtn' + suffix, 'gwxThreshHelpPop' + suffix, 'Significance threshold help', THRESH_DEF_TERMS) +
      '<div class="gwx-seg" role="group" aria-label="Threshold source">' +
        '<button id="gwxThreshPublished' + suffix + '" class="' + (threshSource === 'published' ? 'on' : '') + '" type="button">Published</button>' +
        '<button id="gwxThreshCustom' + suffix + '" class="' + (threshSource === 'custom' ? 'on' : '') + '" type="button">Custom</button>' +
      '</div>' +
      '<div class="gwx-thresh-custom" id="gwxThreshCustomControls' + suffix + '" style="display:' + (threshSource === 'custom' ? 'flex' : 'none') + '">' +
        '<div class="gwx-seg" role="group" aria-label="Threshold method">' +
          '<button id="gwxThreshSimpleM' + suffix + '" class="' + (threshMethod === 'simplem' ? 'on' : '') + '" type="button"' +
            (hasSimpleM ? '' : ' disabled title="No effective marker count available for this dataset"') + '>SimpleM</button>' +
          '<button id="gwxThreshBonferroni' + suffix + '" class="' + (threshMethod === 'bonferroni' ? 'on' : '') + '" type="button"' +
            (hasBonferroni ? '' : ' disabled title="No marker count available for this dataset"') + '>Bonferroni</button>' +
        '</div>' +
        '<label class="gwx-thresh-alpha">&alpha; <input type="text" id="gwxAlphaInput' + suffix + '" value="' + escAttr(threshAlpha) + '" inputmode="decimal" /></label>' +
      '</div>' +
      '<span class="gwx-thresh-readout" id="gwxThreshReadout' + suffix + '">' + threshReadoutText() + '</span>' +
    '</div>';
  }

  function toolbarHTML() {
    const chrOpts = [];
    for (let c = 1; c <= CFG.chrCount; c++) chrOpts.push('<option value="' + c + '">Chromosome ' + c + '</option>');
    return '<div class="card gwx-toolbar">' +
      '<div class="gwx-toolbar-row">' +
        '<div class="gwx-seg" role="group" aria-label="Interaction mode">' +
          '<button id="gwxModePan" class="on" type="button">Pan / Zoom</button>' +
          '<button id="gwxModeSelect" type="button">Select region</button>' +
        '</div>' +
        '<div class="gwx-help-wrap">' +
          '<button class="gwx-help-btn" id="gwxModeHelpBtn" type="button" aria-label="Interaction mode help">?</button>' +
          '<div class="gwx-help-pop" id="gwxModeHelpPop"><b>Pan / Zoom mode</b>: Drag to pan, scroll to zoom, or hold <b>Shift</b> and drag to zoom into a region.<br><b>Select region mode</b>: Click and drag across a peak to inspect its SNPs, and to send the region and/or accessions to SNPVersity.</div>' +
        '</div>' +
        '<div class="gwx-div"></div>' +
        '<button class="btn gwx-btn-sm" id="gwxResetBtn" type="button">Reset view</button>' +
        '<div style="margin-left:auto;display:flex;align-items:center;gap:10px">' +
          '<button class="btn gwx-btn-sm" id="gwxExportAllBtn" type="button">' + (typeof ICONS !== 'undefined' && ICONS.download ? ICONS.download : '') + ' Export significant SNPs</button>' +
          '<div class="gwx-readout" id="gwxReadout">Chr 1–' + CFG.chrCount + ' · whole genome</div>' +
        '</div>' +
      '</div>' +
      '<div class="gwx-toolbar-row">' +
        '<div class="gwx-seg" role="group" aria-label="Pan">' +
          '<button class="gwx-btn-labeled" id="gwxPanLeftBtn" type="button">&larr; Pan left</button>' +
          '<button class="gwx-btn-labeled" id="gwxPanRightBtn" type="button">Pan right &rarr;</button>' +
        '</div>' +
        '<div class="gwx-seg" role="group" aria-label="Zoom">' +
          '<button class="gwx-btn-labeled" id="gwxZoomOutBtn" type="button">&minus; Zoom out</button>' +
          '<button class="gwx-btn-labeled" id="gwxZoomInBtn" type="button">Zoom in +</button>' +
        '</div>' +
        '<select id="gwxChrJump" aria-label="Jump to chromosome"><option value="">Whole genome</option>' + chrOpts.join('') + '</select>' +
        '<div class="gwx-search-wrap">' +
          '<input type="text" id="gwxSearch" placeholder="Jump to SNP or region…" title="SNP ID, chr6:169202719, or chr6:169202719..169209873" />' +
          '<div class="gwx-search-status" id="gwxSearchStatus"></div>' +
        '</div>' +
      '</div>' +
    '</div>';
  }

  function plotHTML() {
    return '<div class="card gwx-plot-wrap">' +
        '<canvas id="gwxPlot"></canvas>' +
        '<div class="gwx-selrect" id="gwxSelrect"></div>' +
        '<div class="gwx-crosshair" id="gwxCrosshair"></div>' +
        '<div class="gwx-sel-label" id="gwxSelLabelStart"></div>' +
        '<div class="gwx-sel-label" id="gwxSelLabelEnd"></div>' +
        '<div class="gwx-tooltip" id="gwxTooltip"></div>' +
      '</div>' +
      '<p class="gwx-plot-legend">' + FILTERED_SNPS_NOTE + '</p>' +
      '<div class="card gwx-minimap-wrap">' +
        '<p class="gwx-minimap-label">Genome overview</p>' +
        '<canvas id="gwxMinimap"></canvas>' +
      '</div>';
  }

  function panelHTML() {
    return '<div class="gwx-panel-backdrop" id="gwxPanelBackdrop"></div>' +
    '<aside class="gwx-panel" id="gwxPanel" aria-label="Selected region SNPs">' +
      '<div class="gwx-panel-head">' +
        '<div class="bar"></div>' +
        '<div class="gwx-panel-head-body">' +
          '<div class="gwx-panel-head-row">' +
            '<div><span class="kicker">Region</span><h2 id="gwxRegionCoord">—</h2><div class="region-sub" id="gwxRegionSub">—</div></div>' +
            '<button class="gwx-panel-close" id="gwxPanelCloseBtn" aria-label="Close">&times;</button>' +
          '</div>' +
          '<div class="gwx-panel-stats">' +
            '<div class="gwx-stat"><div class="v" id="gwxStatTotal">0</div><div class="k">SNPs in region</div></div>' +
            '<div class="gwx-stat sig"><div class="v" id="gwxStatSig">0</div><div class="k">Significant</div></div>' +
          '</div>' +
        '</div>' +
      '</div>' +
      '<div class="gwx-panel-thresh">' + thresholdBarHTML('Panel') + '</div>' +
      '<div class="gwx-panel-controls">' +
        '<div style="display:flex;align-items:center;gap:8px">' +
          '<div class="gwx-seg" role="group" aria-label="Filter">' +
            '<button id="gwxFilterAll" class="on" type="button">All SNPs</button>' +
            '<button id="gwxFilterSig" type="button">Significant only</button>' +
          '</div>' +
          defPopoverHTML('gwxRegionHelpBtn', 'gwxRegionHelpPop', 'Table column definitions', REGION_DEF_TERMS) +
        '</div>' +
        '<span class="gwx-hint" id="gwxShownHint"></span>' +
      '</div>' +
      '<div class="gwx-panel-table-wrap">' +
        '<table class="gwx-table" id="gwxTable"><thead id="gwxTableHead">' + tableHeadHTML() + '</thead><tbody id="gwxTableBody"></tbody></table>' +
        '<div class="gwx-empty" id="gwxEmpty" style="display:none">No SNPs match this filter in the selected region.</div>' +
      '</div>' +
      '<div class="gwx-panel-foot">' +
        '<div class="row">' +
          '<button class="btn" id="gwxExportRegionBtn" type="button">Export shown SNPs (CSV)</button>' +
          '<button class="btn primary" id="gwxOpenSendBtn" type="button">Send to SNPVersity →</button>' +
        '</div>' +
      '</div>' +
    '</aside>' +
    '<div class="gwx-send-popup-backdrop" id="gwxSendPopupBackdrop"></div>' +
    '<div class="gwx-send-popup" id="gwxSendPopup" role="dialog" aria-label="Send to SNPVersity">' +
      '<div class="gwx-send-popup-head">' +
        '<h3>Send to SNPVersity</h3>' +
        '<button class="gwx-panel-close" id="gwxSendPopupCloseBtn" aria-label="Close">&times;</button>' +
      '</div>' +
      '<div class="gwx-send-popup-body">' +
        '<label class="gwx-ho-check"><input type="checkbox" id="gwxSendRegionChk">' +
          '<span>Send this genomic region</span></label>' +
        (activeEntry.population === 'NAM'
          ? '<label class="gwx-ho-check"><input type="checkbox" id="gwxSendAccReplaceChk">' +
              '<span>Send accessions — replace the <b id="gwxAccCountReplace">0</b> accessions currently selected in SNPVersity</span></label>' +
            '<label class="gwx-ho-check"><input type="checkbox" id="gwxSendAccAddChk">' +
              '<span>Send accessions — keep the <b id="gwxAccCountAdd">0</b> accessions currently selected in SNPVersity and add to them</span></label>'
          : '') +
        '<div class="gwx-ho-hint" id="gwxHoHint"></div>' +
      '</div>' +
      '<div class="gwx-send-popup-foot">' +
        '<button class="btn" id="gwxSendPopupCancelBtn" type="button">Cancel</button>' +
        '<button class="btn primary" id="gwxSendConfirmBtn" type="button">Send</button>' +
      '</div>' +
    '</div>';
  }

  function styleCSS() {
    return '' +
    '.gwx-picker{margin-top:14px}' +
    '.gwx-picker-row{display:flex;gap:14px;flex-wrap:wrap}' +
    '.gwx-picker-field{flex:1 1 200px;min-width:180px;margin:0}' +
    '.gwx-picker-active{margin-top:12px;font-size:12.5px;color:var(--muted)}' +
    '.gwx-picker-doi{margin-top:4px;font-size:11.5px;color:var(--faint)}' +
    '.gwx-picker-active b{color:var(--ink)}' +
    '.gwx-picker-empty{margin-top:12px;font-size:12.5px;color:var(--muted)}' +
    '.gwx-picker-clear{appearance:none;border:none;background:transparent;color:var(--blue-600);font:inherit;font-size:12.5px;font-weight:600;cursor:pointer;padding:0;margin-left:2px;text-decoration:underline}' +
    '.gwx-picker-chips{display:flex;flex-wrap:wrap;gap:8px;margin-top:10px}' +
    '.gwx-picker-export-row{margin-top:12px;padding-top:12px;border-top:1px solid var(--line)}' +
    '.gwx-export-scope{display:flex;flex-direction:column;gap:10px}' +
    '.gwx-export-scope em{font-style:normal;color:var(--faint)}' +
    '.gwx-export-progress{font-size:12.5px;color:var(--blue-600);font-weight:600}' +
    '.gwx-export-note{font-size:12px;color:#9a6700;background:#fff8e6;border:1px solid #f0dca0;border-radius:var(--rad-sm);padding:8px 10px}' +
    '.gwx-export-warn{color:#b42318;background:#fff6f5;border-color:#f3c2bd}' +
    '.gwx-export-floor-note{font-size:11.5px;color:var(--faint)}' +
    '.gwx-chip{appearance:none;border:1px solid var(--line);background:#fff;border-radius:10px;padding:8px 12px;font:inherit;font-size:12.5px;font-weight:600;color:var(--ink);cursor:pointer;text-align:left;transition:.12s}' +
    '.gwx-chip:hover{border-color:var(--blue);background:var(--blue-50)}' +
    '.gwx-chip.on{border-color:var(--blue);background:var(--blue-50);color:var(--blue-600)}' +
    '.gwx-chip-sub{display:block;font-weight:500;font-size:11px;color:var(--muted);margin-top:2px}' +
    '.gwx-thresh-bar{display:flex;align-items:center;gap:12px;flex-wrap:wrap;padding:9px 12px;margin-top:14px}' +
    '.gwx-thresh-label{font-size:12.5px;font-weight:600;color:var(--muted)}' +
    '.gwx-thresh-custom{display:flex;align-items:center;gap:12px;flex-wrap:wrap}' +
    '.gwx-thresh-alpha{display:flex;align-items:center;gap:6px;font-size:12.5px;font-weight:600;color:var(--ink)}' +
    '.gwx-thresh-alpha input{width:70px;font:inherit;font-family:var(--mono);font-size:12.5px;padding:6px 8px;border-radius:var(--rad-sm);border:1px solid var(--line);background:#fff;color:var(--ink)}' +
    '.gwx-thresh-alpha input:focus{outline:none;border-color:var(--blue);box-shadow:0 0 0 3px var(--blue-50)}' +
    '.gwx-thresh-readout{margin-left:auto;font-family:var(--mono);font-size:12px;color:var(--blue-600);font-variant-numeric:tabular-nums}' +
    '.gwx-thresh-readout b{color:var(--ink)}' +
    '.gwx-toolbar{display:flex;flex-direction:column;gap:10px;padding:9px 12px;margin-top:14px}' +
    '.gwx-toolbar-row{display:flex;align-items:center;gap:10px;flex-wrap:wrap}' +
    '.gwx-seg{display:inline-flex;background:var(--line-2);border-radius:999px;padding:3px;gap:2px}' +
    '.gwx-seg button{appearance:none;border:none;background:transparent;color:var(--muted);font:inherit;font-size:12.5px;font-weight:600;padding:6px 13px;border-radius:999px;cursor:pointer;transition:background .12s,color .12s}' +
    '.gwx-seg button.on{background:var(--blue);color:#fff}' +
    '.gwx-toolbar select{font:inherit;font-size:12.5px;font-weight:600;padding:7px 11px;border-radius:var(--rad-sm);border:1px solid var(--line);background:#fff;color:var(--ink);width:auto;flex:0 0 auto}' +
    '.gwx-toolbar select:focus,.gwx-toolbar input:focus{outline:none;border-color:var(--blue);box-shadow:0 0 0 3px var(--blue-50)}' +
    '.gwx-toolbar input[type=text]{font:inherit;font-family:var(--mono);font-size:12.5px;padding:7px 11px;border-radius:var(--rad-sm);border:1px solid var(--line);background:#fff;color:var(--ink);width:210px}' +
    '.gwx-search-wrap{position:relative}' +
    '.gwx-search-status{position:absolute;top:100%;left:0;margin-top:6px;background:#fff8ec;border:1px solid #f3e3c2;color:#7a5a17;font-size:11.5px;line-height:1.4;padding:6px 10px;border-radius:8px;box-shadow:var(--shadow-lg);white-space:nowrap;z-index:15;display:none}' +
    '.gwx-search-status.show{display:block}' +
    '.gwx-help-wrap{position:relative;display:inline-flex}' +
    '.gwx-help-btn{width:20px;height:20px;border-radius:50%;border:1px solid var(--line);background:#fff;color:var(--muted);font:inherit;font-size:11px;font-weight:700;cursor:pointer;display:inline-flex;align-items:center;justify-content:center;padding:0;line-height:1}' +
    '.gwx-help-btn:hover{border-color:var(--blue);color:var(--blue-600)}' +
    '.gwx-help-pop{position:absolute;top:100%;left:0;margin-top:6px;background:#fff;border:1px solid var(--line);color:var(--muted);font-size:11.5px;line-height:1.45;padding:8px 10px;border-radius:8px;box-shadow:var(--shadow-lg);white-space:normal;width:280px;max-height:260px;overflow-y:auto;z-index:20;display:none}' +
    '.gwx-help-pop b{color:var(--ink)}' +
    '.gwx-help-wrap:hover .gwx-help-pop,.gwx-help-pop.show{display:block}' +
    '.gwx-trait-dd{position:relative}' +
    '.gwx-trait-dd-btn{width:100%;text-align:left;appearance:none;border:1px solid var(--line);border-radius:9px;padding:9px 11px;font:inherit;font-size:13.5px;font-weight:600;background:#fff;color:var(--ink);cursor:pointer;display:flex;align-items:center;justify-content:space-between;gap:8px}' +
    '.gwx-trait-dd-btn:hover{border-color:var(--blue)}' +
    '.gwx-trait-dd-btn::after{content:"\\25BE";color:var(--muted);font-size:11px}' +
    '.gwx-trait-dd-pop{position:absolute;top:100%;left:0;right:0;margin-top:4px;background:#fff;border:1px solid var(--line);border-radius:10px;box-shadow:var(--shadow-lg);z-index:20;display:none;max-height:360px;overflow-y:auto;padding:6px}' +
    '.gwx-trait-dd-pop.show{display:block}' +
    '.gwx-trait-dd-cat{display:block;width:100%;text-align:left;appearance:none;border:none;background:var(--blue-50);color:var(--blue-600);font:inherit;font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:.4px;padding:7px 10px;border-radius:6px;cursor:pointer;margin:8px 0 3px}' +
    '.gwx-trait-dd-cat:first-of-type{margin-top:0}' +
    '.gwx-trait-dd-cat:hover{border-color:var(--blue)}' +
    '.gwx-trait-dd-cat.on{background:var(--blue);color:#fff}' +
    '.gwx-trait-dd-opt{display:block;width:100%;text-align:left;appearance:none;border:none;background:transparent;color:var(--ink);font:inherit;font-size:13px;padding:7px 10px 7px 18px;border-radius:6px;cursor:pointer}' +
    '.gwx-trait-dd-opt:hover{background:var(--blue-50)}' +
    '.gwx-trait-dd-opt.on{background:var(--blue-50);color:var(--blue-600);font-weight:600}' +
    '.gwx-trait-dd-count{color:var(--muted);font-weight:500}' +
    '.gwx-trait-dd-cat.on .gwx-trait-dd-count{color:rgba(255,255,255,.8)}' +
    '.gwx-div{width:1px;align-self:stretch;background:var(--line);margin:2px 2px}' +
    '.gwx-readout{font-family:var(--mono);font-size:12px;color:var(--blue-600);font-weight:500;font-variant-numeric:tabular-nums;text-align:right;white-space:nowrap}' +
    '.gwx-readout b{color:var(--ink)}' +
    '.gwx-btn-sm{padding:6px 11px;font-size:12.5px}' +
    '.gwx-btn-icon{padding:6px 10px}' +
    '.gwx-btn-labeled{padding:6px 14px;white-space:nowrap}' +
    '.gwx-plot-wrap{position:relative;margin-top:14px}' +
    '#gwxPlot{display:block;width:100%;height:440px;cursor:grab;border-radius:var(--rad)}' +
    '#gwxPlot.mode-select{cursor:crosshair}' +
    '#gwxPlot.shift-zoom{cursor:crosshair}' +
    '#gwxPlot.dragging{cursor:grabbing}' +
    '.gwx-selrect{position:absolute;top:0;bottom:0;background:rgba(207,138,18,.20);border-left:1.5px solid var(--gold);border-right:1.5px solid var(--gold);pointer-events:none;display:none}' +
    /* top/bottom match MARGIN.top(16)/MARGIN.bottom(34) so the line spans exactly the plot area */
    '.gwx-crosshair{position:absolute;top:16px;bottom:34px;width:0;border-left:1.5px dashed var(--blue-600);pointer-events:none;display:none;z-index:4}' +
    '.gwx-sel-label{position:absolute;top:3px;transform:translateX(-50%);background:var(--gold);color:#fff;font-family:var(--mono);font-size:11px;font-weight:700;padding:2px 7px;border-radius:6px;white-space:nowrap;pointer-events:none;display:none;z-index:6}' +
    '.gwx-tooltip{position:absolute;pointer-events:none;background:#fff;border:1px solid var(--line);border-radius:10px;box-shadow:var(--shadow-lg);padding:9px 11px;font-size:12px;line-height:1.5;color:var(--ink);min-width:198px;z-index:5;display:none}' +
    '.gwx-tooltip .coord{font-family:var(--mono);font-weight:700;font-size:12.5px;color:var(--blue-600);margin-bottom:2px}' +
    '.gwx-tooltip .snp{font-family:var(--mono);font-weight:600;color:var(--muted);font-size:11px;margin-bottom:5px}' +
    '.gwx-tooltip .row{display:flex;justify-content:space-between;gap:14px;color:var(--muted);white-space:nowrap}' +
    '.gwx-tooltip .row b{font-family:var(--mono);color:var(--ink);font-weight:600;font-variant-numeric:tabular-nums}' +
    '.gwx-sig-badge{display:inline-block;margin-top:4px;font-size:10.5px;font-weight:700;padding:1px 7px;border-radius:20px;background:#c0392b;color:#fff}' +
    '.gwx-plot-legend{margin:8px 2px 0;font-size:10.5px;line-height:1.5;color:var(--muted);text-align:center}' +
    '.gwx-minimap-wrap{margin-top:14px;padding:10px 12px 12px}' +
    '.gwx-minimap-label{font-size:10.5px;color:var(--faint);text-transform:uppercase;letter-spacing:.06em;margin:0 0 6px;font-weight:700}' +
    '#gwxMinimap{display:block;width:100%;height:54px;cursor:pointer;border-radius:var(--rad-sm)}' +
    '.gwx-footnote{font-size:11.5px;color:var(--faint);text-align:center;margin-top:10px}' +
    '.gwx-footnote b{color:var(--muted)}' +
    /* z-index above 41 — the site's own .rail nav sits at z-index:40 and would
       otherwise paint over the left edge of a centered (as opposed to the old
       right-docked) modal */
    '.gwx-panel-backdrop{position:fixed;inset:0;background:rgba(10,15,28,.45);z-index:41;display:none}' +
    '.gwx-panel-backdrop.open{display:block}' +
    /* max-height (not a fixed height) so the modal only grows as tall as its
       content needs, up to the cap; combined with the table's own min-height
       below, this stops a tall footer from starving the table down to
       nothing on a shorter laptop viewport. overflow-y is a last-resort
       fallback — normally the table's own scrollbar is the only one that
       ever shows. */
    '.gwx-panel{position:fixed;top:50%;left:50%;width:min(1100px,94vw);max-height:min(88vh,900px);background:#fff;border:1px solid var(--line);border-radius:var(--rad);box-shadow:var(--shadow-lg);z-index:42;display:flex;flex-direction:column;overflow-y:auto;overflow-x:hidden;transform:translate(-50%,-50%) scale(.96);opacity:0;pointer-events:none;transition:transform .16s ease,opacity .16s ease}' +
    '.gwx-panel.open{transform:translate(-50%,-50%) scale(1);opacity:1;pointer-events:auto}' +
    '.gwx-panel-head{display:flex;gap:12px;padding:18px 18px 14px;border-bottom:1px solid var(--line)}' +
    '.gwx-panel-head .bar{width:4px;align-self:stretch;min-height:30px;border-radius:4px;background:var(--green);flex:none;margin-top:2px}' +
    '.gwx-panel-head-body{flex:1;min-width:0}' +
    '.gwx-panel-head-row{display:flex;align-items:flex-start;justify-content:space-between;gap:10px}' +
    '.gwx-panel .kicker{font-family:var(--disp);font-weight:700;font-size:9.5px;color:var(--green-600);letter-spacing:.12em;text-transform:uppercase;background:var(--green-50);padding:3px 8px;border-radius:5px;display:inline-block;margin-bottom:6px}' +
    '.gwx-panel h2{font-family:var(--disp);font-size:15px;font-weight:700;margin:0 0 2px}' +
    '.gwx-panel .region-sub{font-family:var(--mono);font-size:12px;color:var(--muted)}' +
    '.gwx-panel-close{appearance:none;border:none;background:transparent;color:var(--faint);font-size:18px;line-height:1;cursor:pointer;padding:4px 6px;border-radius:6px}' +
    '.gwx-panel-close:hover{background:var(--line-2);color:var(--ink)}' +
    '.gwx-panel-stats{display:flex;gap:20px;margin-top:10px}' +
    '.gwx-stat .v{font-family:var(--mono);font-size:18px;font-weight:700;color:var(--navy-700)}' +
    '.gwx-stat.sig .v{color:#c0392b}' +
    '.gwx-stat .k{font-size:10.5px;color:var(--faint);text-transform:uppercase;letter-spacing:.05em;font-weight:600}' +
    '.gwx-panel-thresh{padding:12px 18px 0}' +
    '.gwx-panel-thresh .gwx-thresh-bar{margin-top:0;padding:9px 12px}' +
    '.gwx-panel-controls{display:flex;align-items:center;justify-content:space-between;gap:10px;padding:12px 18px;border-bottom:1px solid var(--line)}' +
    '.gwx-panel-table-wrap{flex:1 1 auto;min-height:220px;overflow:auto;padding:0 0 8px}' +
    '.gwx-panel-head,.gwx-panel-controls,.gwx-panel-foot{flex:0 0 auto}' +
    '.gwx-table{width:100%;border-collapse:collapse;font-size:12px;font-family:var(--mono)}' +
    '.gwx-table th{position:sticky;top:0;background:#fff;text-align:left;font-family:var(--body);font-size:10px;text-transform:uppercase;letter-spacing:.05em;color:var(--faint);font-weight:700;padding:8px 10px;border-bottom:1px solid var(--line);white-space:nowrap}' +
    '.gwx-table th.gwx-sortable{cursor:pointer;user-select:none}' +
    '.gwx-table th.gwx-sortable:hover{color:var(--blue-600)}' +
    '.gwx-table th.sorted{color:var(--blue-600)}' +
    '.gwx-table td{padding:6px 10px;border-bottom:1px solid var(--line);font-variant-numeric:tabular-nums;color:var(--muted);white-space:nowrap}' +
    '.gwx-table td.snpid{color:var(--navy-700);font-weight:700}' +
    '.gwx-table tr.sig td.snpid::before{content:"";display:inline-block;width:6px;height:6px;border-radius:50%;background:#c0392b;margin-right:6px}' +
    '.gwx-panel-foot{padding:12px 18px 16px;border-top:1px solid var(--line);display:flex;flex-direction:column;gap:8px}' +
    '.gwx-panel-foot .row{display:flex;gap:8px}' +
    '.gwx-ho-check{display:flex;align-items:baseline;gap:8px;font-size:13px;color:var(--ink);cursor:pointer;line-height:1.4}' +
    '.gwx-ho-check input{margin:0;position:relative;top:2px;cursor:pointer;flex:0 0 auto}' +
    '.gwx-ho-check input:disabled{cursor:not-allowed}' +
    '.gwx-ho-check input:disabled ~ span{color:var(--faint)}' +
    '.gwx-ho-check b{font-family:var(--mono)}' +
    '.gwx-send-popup-backdrop{position:fixed;inset:0;background:rgba(10,15,28,.35);z-index:43;display:none}' +
    '.gwx-send-popup-backdrop.open{display:block}' +
    '.gwx-send-popup{position:fixed;top:50%;left:50%;width:min(480px,90vw);max-height:85vh;background:#fff;border:1px solid var(--line);border-radius:var(--rad);box-shadow:var(--shadow-lg);z-index:44;display:flex;flex-direction:column;overflow:hidden;transform:translate(-50%,-50%) scale(.96);opacity:0;pointer-events:none;transition:transform .16s ease,opacity .16s ease}' +
    '.gwx-send-popup.open{transform:translate(-50%,-50%) scale(1);opacity:1;pointer-events:auto}' +
    '.gwx-send-popup-head{display:flex;align-items:center;justify-content:space-between;padding:16px 18px 12px;border-bottom:1px solid var(--line)}' +
    '.gwx-send-popup-head h3{font-family:var(--disp);font-size:15px;font-weight:700;margin:0;color:var(--ink)}' +
    '.gwx-send-popup-body{padding:16px 18px;display:flex;flex-direction:column;gap:14px;overflow:auto}' +
    '.gwx-send-popup-foot{padding:12px 18px 16px;border-top:1px solid var(--line);display:flex;justify-content:flex-end;gap:8px}' +
    '.gwx-empty{padding:30px 18px;text-align:center;color:var(--faint);font-size:12.5px}' +
    '.gwx-hint{font-size:12px;color:var(--faint)}' +
    '.gwx-ho-hint{font-size:11px;color:var(--faint);min-height:14px}' +
    '@media (max-width:720px){ .gwx-readout{display:none} #gwxPlot{height:340px} }';
  }

  function datasetBodyHTML() {
    if (datasetError) return datasetErrorHTML();
    if (!activeEntry) return '';
    if (!DATA) return datasetLoadingHTML();
    return introHTML() + thresholdBarHTML() + toolbarHTML() + plotHTML() + panelHTML();
  }

  function shellHTML() {
    return '<style>' + styleCSS() + '</style>' + generalHeaderHTML() + pickerHTML() + datasetBodyHTML();
  }

  /* ------------------------------------------------------------------ *
   *  RENDER — router entry point                                       *
   * ------------------------------------------------------------------ */
  function render(page) {
    page.className = 'page fade';
    if (!MANIFEST) {
      if (manifestError) {
        page.innerHTML = generalHeaderHTML() + manifestErrorHTML();
        const retry = document.getElementById('gwxManifestRetry');
        if (retry) retry.addEventListener('click', function () { manifestError = null; render(page); });
        return;
      }
      page.innerHTML = generalHeaderHTML() + manifestLoadingHTML();
      loadManifest().then(rerenderIfActive).catch(rerenderIfActive);
      return;
    }
    page.innerHTML = shellHTML();
    wirePicker();
    if (activeEntry && DATA) {
      wireInteractions();
      updateReadout();
      resizeCanvases();
    }
  }

  SNPTools.register('snpgwas', { render: render });
})();
