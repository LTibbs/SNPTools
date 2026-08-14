/* =====================================================================
 *  snpgwas.js — GWAS Explorer (Manhattan plot) as a native SNPTools panel.
 *  Registers 'snpgwas'. Loads AFTER core.js (needs SNPTools, S, go()).
 *
 *  Ported from a standalone prototype (see GWAS_EXPLORER_HANDOFF.md) —
 *  the canvas rendering, zoom/pan math, binary-search region lookup and
 *  minimap are unchanged. What changed to fit the suite:
 *    - data is fetched from ./data/gwas/*.csv instead of embedded inline
 *    - CSV export uses a Blob + temporary <a download> instead of the
 *      Claude Artifacts sandbox's window.claude.downloads API
 *    - a "Send region to SNPVersity" button hands the selected interval
 *      off through the same window.versityRequest(...) entry point
 *      SNPFunction and friends use (see snpversity.js applyPendingRequest)
 * ===================================================================== */
(function () {
  'use strict';

  const CFG = {
    csvUrl:     './data/gwas/GCTA_LW_intercept_fixed.csv',
    threshP:    6.234414e-05,
    threshLabel:'6.234414e-05',
    traitLabel: 'LW — intercept (fixed effect)',
    chrCount:   10,
  };
  const THRESH_NEGLOG = -Math.log10(CFG.threshP);

  /* ------------------------------------------------------------------ *
   *  DATA LOAD + PARSE  (fetched once, cached for the session)          *
   * ------------------------------------------------------------------ */
  let DATA = null;
  let loadPromise = null;
  let loadError = null;

  function loadData() {
    if (DATA) return Promise.resolve(DATA);
    if (loadPromise) return loadPromise;
    loadPromise = fetch(CFG.csvUrl).then(function (r) {
      if (!r.ok) throw new Error('HTTP ' + r.status + ' fetching ' + CFG.csvUrl);
      return r.text();
    }).then(function (raw) {
      DATA = parseCsv(raw);
      view = { x0: 0, x1: DATA.TOTAL_CUM };
      return DATA;
    }).catch(function (err) {
      loadError = err;
      loadPromise = null;
      throw err;
    });
    return loadPromise;
  }

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
    const chrLen = new Float64Array(CHR_COUNT + 1);
    for (let i = 0; i < n; i++) {
      const c = chrOf[i];
      if (chrStart[c] === -1) chrStart[c] = i;
      chrEnd[c] = i + 1;
      if (bpOf[i] > chrLen[c]) chrLen[c] = bpOf[i];
    }

    let totalLen = 0;
    for (let c = 1; c <= CHR_COUNT; c++) totalLen += chrLen[c];
    const GAP = (totalLen / CHR_COUNT) * 0.03;

    const chrOffset = new Float64Array(CHR_COUNT + 1);
    let cum = 0;
    for (let c = 1; c <= CHR_COUNT; c++) { chrOffset[c] = cum; cum += chrLen[c] + GAP; }
    const TOTAL_CUM = cum - GAP;

    let globalMaxNegLog = 0;
    for (let i = 0; i < n; i++) if (negLogP[i] > globalMaxNegLog) globalMaxNegLog = negLogP[i];
    const Y_MAX = globalMaxNegLog * 1.10;

    let totalSig = 0;
    for (let i = 0; i < n; i++) if (negLogP[i] >= THRESH_NEGLOG) totalSig++;

    return {
      n: n, chrOf: chrOf, bpOf: bpOf, negLogP: negLogP, fields: fields,
      chrStart: chrStart, chrEnd: chrEnd, chrLen: chrLen, chrOffset: chrOffset,
      TOTAL_CUM: TOTAL_CUM, Y_MAX: Y_MAX, totalSig: totalSig,
    };
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
  function cumAt(idx) { return DATA.chrOffset[DATA.chrOf[idx]] + DATA.bpOf[idx]; }
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
      if (x >= DATA.chrOffset[c] && x <= DATA.chrOffset[c] + DATA.chrLen[c]) return c;
    }
    for (let c = 1; c <= CFG.chrCount; c++) if (x < DATA.chrOffset[c]) return c > 1 ? c - 1 : c;
    return CFG.chrCount;
  }

  function fmtBp(bp) { return Math.round(bp).toLocaleString('en-US'); }
  function fmtMb(bp) { return (bp / 1e6).toFixed(2) + ' Mb'; }
  function fmtP(pStr) { const v = parseFloat(pStr); return isFinite(v) ? v.toExponential(2) : pStr; }
  function escAttr(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }

  /* ------------------------------------------------------------------ *
   *  VIEW / INTERACTION STATE — module scope, persists across visits   *
   *  to this tool (so pan/zoom position survives navigating away and   *
   *  back); DOM refs are re-queried fresh every render() since #page   *
   *  is rebuilt from scratch on each visit.                            *
   * ------------------------------------------------------------------ */
  let view = null;           // {x0,x1} in cumulative-bp space
  let mode = 'pan';          // 'pan' | 'select'
  const MARGIN = { top: 16, right: 20, bottom: 34, left: 54 };
  const MIN_SPAN = 300;

  let DOM = {};
  let plotW = 0, plotH = 0, dpr = 1;
  let dragging = false, dragKind = 'pan', dragStartPx = 0, dragStartView = null, dragStartCum = 0;
  let miniDragging = false;
  let miniBins = null;

  let currentRegion = null;  // {indices, sigIndices, chrsTouched, bpStart, bpEnd, chrLabel, coordLabel, spanLabel}
  let regionFilter = 'all';
  const TABLE_CAP = 1500;

  function clampView(x0, x1) {
    let span = x1 - x0;
    if (span < MIN_SPAN) span = MIN_SPAN;
    if (span > DATA.TOTAL_CUM) span = DATA.TOTAL_CUM;
    if (x0 < 0) { x0 = 0; x1 = span; }
    if (x1 > DATA.TOTAL_CUM) { x1 = DATA.TOTAL_CUM; x0 = DATA.TOTAL_CUM - span; }
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
    if (view.x0 <= 0 && view.x1 >= DATA.TOTAL_CUM) {
      label = 'Whole genome (' + fmtMb(DATA.TOTAL_CUM) + ')';
    } else {
      const startChr = chrForCum(view.x0), endChr = chrForCum(view.x1);
      if (startChr === endChr) {
        const bpStart = Math.max(0, view.x0 - DATA.chrOffset[startChr]);
        const bpEnd = Math.min(DATA.chrLen[startChr], view.x1 - DATA.chrOffset[startChr]);
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
      panelBg: v('--white', '#fff'), muted: v('--muted', '#62718a'),
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
      const segX0 = DATA.chrOffset[c], segX1 = DATA.chrOffset[c] + DATA.chrLen[c];
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
      const bpSpan = Math.min(DATA.chrLen[c], view.x1 - DATA.chrOffset[c]) - Math.max(0, view.x0 - DATA.chrOffset[c]);
      const step = niceStep(bpSpan / 6);
      const startTick = Math.ceil(Math.max(0, view.x0 - DATA.chrOffset[c]) / step) * step;
      ctx.font = '10.5px ui-monospace, Menlo, monospace';
      ctx.fillStyle = col.muted;
      for (let bp = startTick; bp <= Math.min(DATA.chrLen[c], view.x1 - DATA.chrOffset[c]); bp += step) {
        const px = xToPx(DATA.chrOffset[c] + bp);
        if (px < plotLeft || px > plotRight) continue;
        ctx.strokeStyle = col.grid;
        ctx.beginPath(); ctx.moveTo(px, plotTop); ctx.lineTo(px, plotBottom); ctx.stroke();
        const label = step >= 1e6 ? (bp / 1e6).toFixed(bp === 0 ? 0 : 1) + ' Mb' : (bp / 1e3).toFixed(0) + ' kb';
        ctx.fillText(label, px, plotBottom + 22);
      }
    }

    const thY = yToPx(THRESH_NEGLOG);
    if (thY >= plotTop && thY <= plotBottom) {
      ctx.save();
      ctx.strokeStyle = col.muted;
      ctx.setLineDash([5, 4]);
      ctx.lineWidth = 1.25;
      ctx.beginPath(); ctx.moveTo(plotLeft, thY); ctx.lineTo(plotRight, thY); ctx.stroke();
      ctx.restore();
      ctx.fillStyle = col.muted;
      ctx.font = '10.5px system-ui, sans-serif';
      ctx.textAlign = 'left'; ctx.textBaseline = 'bottom';
      ctx.fillText('SimpleM threshold (p < ' + CFG.threshLabel + ')', plotLeft + 6, thY - 3);
    }

    let pointR = 1.6;
    const spanFrac = (view.x1 - view.x0) / DATA.TOTAL_CUM;
    if (spanFrac < 0.25) pointR = 2.3;
    if (spanFrac < 0.05) pointR = 3.0;
    if (spanFrac < 0.01) pointR = 3.8;

    for (let c = 1; c <= CFG.chrCount; c++) {
      const segX0 = DATA.chrOffset[c], segX1 = DATA.chrOffset[c] + DATA.chrLen[c];
      if (segX1 < view.x0 || segX0 > view.x1) continue;
      const bpLo = Math.max(0, view.x0 - DATA.chrOffset[c]);
      const bpHi = Math.min(DATA.chrLen[c], view.x1 - DATA.chrOffset[c]);
      const s = DATA.chrStart[c], e = DATA.chrEnd[c];
      if (s === -1) continue;
      const lo = lowerBound(s, e, DATA.bpOf, bpLo), hi = upperBound(s, e, DATA.bpOf, bpHi);
      const baseColor = (c % 2 === 1) ? col.chrA : col.chrB;
      for (let idx = lo; idx < hi; idx++) {
        const px = xToPx(DATA.chrOffset[c] + DATA.bpOf[idx]);
        const py = yToPx(DATA.negLogP[idx]);
        const isSig = DATA.negLogP[idx] >= THRESH_NEGLOG;
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
    const binW = DATA.TOTAL_CUM / binCount;
    const bins = new Array(binCount);
    for (let b = 0; b < binCount; b++) bins[b] = { maxNL: 0, chr: 0 };
    for (let i = 0; i < DATA.n; i++) {
      const cp = DATA.chrOffset[DATA.chrOf[i]] + DATA.bpOf[i];
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
      const isSig = bin.maxNL >= THRESH_NEGLOG;
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
    const x0 = (view.x0 / DATA.TOTAL_CUM) * w;
    const x1 = (view.x1 / DATA.TOTAL_CUM) * w;
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
    const bpA = Math.max(0, Math.min(DATA.chrLen[cA], cumMin - DATA.chrOffset[cA]));
    const bpB = Math.max(0, Math.min(DATA.chrLen[cB], cumMax - DATA.chrOffset[cB]));
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
    const cp = (px / rect.width) * DATA.TOTAL_CUM;
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
    const isSig = DATA.negLogP[idx] >= THRESH_NEGLOG;
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
    const bp = Math.max(0, Math.min(DATA.chrLen[c], cumPos - DATA.chrOffset[c]));
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
      setView(DATA.chrOffset[chrNum] + Math.max(0, bpLo - pad), DATA.chrOffset[chrNum] + Math.min(DATA.chrLen[chrNum], bpHi + pad));
    } else {
      searchWarn('No SNPs in Chr' + chrNum + ':' + fmtBp(bpLo) + '–' + fmtBp(bpHi) + ' — showing the surrounding region.');
      setView(DATA.chrOffset[chrNum] + Math.max(0, bpLo - 1000), DATA.chrOffset[chrNum] + Math.min(DATA.chrLen[chrNum], bpHi + 1000));
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
      const halfSpan = Math.max(MIN_SPAN, DATA.chrLen[chrNum] * 0.01);
      const cp = DATA.chrOffset[chrNum] + bp;
      setView(cp - halfSpan, cp + halfSpan);
    } else {
      searchWarn('No SNP at Chr' + chrNum + ':' + fmtBp(bp) + ' — showing the surrounding region.');
      setView(DATA.chrOffset[chrNum] + Math.max(0, bp - 1000), DATA.chrOffset[chrNum] + Math.min(DATA.chrLen[chrNum], bp + 1000));
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
        const cp = DATA.chrOffset[DATA.chrOf[i]] + DATA.bpOf[i];
        const halfSpan = Math.max(MIN_SPAN, DATA.chrLen[DATA.chrOf[i]] * 0.01);
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
      const segX0 = DATA.chrOffset[c], segX1 = DATA.chrOffset[c] + DATA.chrLen[c];
      if (segX1 < cumMin || segX0 > cumMax) continue;
      const bpLo = Math.max(0, cumMin - DATA.chrOffset[c]);
      const bpHi = Math.min(DATA.chrLen[c], cumMax - DATA.chrOffset[c]);
      const s = DATA.chrStart[c], e = DATA.chrEnd[c];
      if (s === -1) continue;
      const lo = lowerBound(s, e, DATA.bpOf, bpLo), hi = upperBound(s, e, DATA.bpOf, bpHi);
      if (hi > lo) chrsTouched.push(c);
      for (let idx = lo; idx < hi; idx++) {
        indices.push(idx);
        if (DATA.negLogP[idx] >= THRESH_NEGLOG) sigIndices.push(idx);
      }
    }
    const startChr = chrsTouched[0] || chrForCum(cumMin);
    const endChr = chrsTouched[chrsTouched.length - 1] || chrForCum(cumMax);
    const bpStart = Math.max(0, cumMin - DATA.chrOffset[startChr]);
    const bpEnd = Math.max(0, cumMax - DATA.chrOffset[endChr]);
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

    const list = regionFilter === 'sig' ? currentRegion.sigIndices : currentRegion.indices;
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
        const isSig = DATA.negLogP[idx] >= THRESH_NEGLOG;
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
    updateSendVersityState();
  }

  function closePanel() {
    if (DOM.panel) DOM.panel.classList.remove('open');
    if (DOM.panelBackdrop) DOM.panelBackdrop.classList.remove('open');
  }

  /* ------------------------------------------------------------------ *
   *  SEND REGION TO SNPVERSITY — the cross-tool handoff. A region has   *
   *  no meaning across chromosomes in SNPVersity's single-S.chr model,  *
   *  so this only enables when the selection sits on one chromosome.    *
   * ------------------------------------------------------------------ */
  function updateSendVersityState() {
    if (!DOM.sendVersityBtn) return;
    const single = !!currentRegion && currentRegion.chrsTouched.length === 1;
    DOM.sendVersityBtn.disabled = !single;
    if (DOM.hoHint) {
      DOM.hoHint.textContent = single ? ''
        : 'Select within a single chromosome to hand this region off to SNPVersity.';
    }
  }

  function sendRegionToVersity() {
    if (!currentRegion || currentRegion.chrsTouched.length !== 1) return;
    if (typeof window.versityRequest !== 'function') return;
    window.versityRequest({
      chr: 'chr' + currentRegion.chr,
      start: currentRegion.bpStart,
      end: currentRegion.bpEnd,
      from: 'GWAS Explorer',
      note: CFG.traitLabel + ' — ' + currentRegion.coordLabel,
      /* GWAS Explorer has no accession list to offer — 'add' with an empty
         accessions array leaves whatever the user already picked in
         SNPVersity untouched, instead of applyPendingRequest's default
         'replace' silently clearing it down to zero. */
      merge: 'add',
    });
  }

  /* ------------------------------------------------------------------ *
   *  CSV EXPORT — Blob + temporary <a download> (no Claude sandbox here)*
   * ------------------------------------------------------------------ */
  const CSV_HEADER = 'Chr,SNP,bp,A1,A2,Freq,b,se,p,neg_log10_p\n';

  function buildCsv(indices) {
    const out = new Array(indices.length);
    for (let k = 0; k < indices.length; k++) {
      const idx = indices[k], f = DATA.fields[idx];
      out[k] = f[0] + ',' + f[1] + ',' + f[2] + ',' + f[3] + ',' + f[4] + ',' + f[5] + ',' + f[6] + ',' + f[7] + ',' + f[8] + ',' + DATA.negLogP[idx].toFixed(4);
    }
    return CSV_HEADER + out.join('\n') + '\n';
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
    for (let i = 0; i < DATA.n; i++) if (DATA.negLogP[i] >= THRESH_NEGLOG) sigAll.push(i);
    downloadCsv('GWAS_LW_intercept_significant_SNPs_' + sigAll.length, buildCsv(sigAll));
  }

  function exportRegion() {
    if (!currentRegion) return;
    const list = regionFilter === 'sig' ? currentRegion.sigIndices : currentRegion.indices;
    const label = currentRegion.chrLabel.replace(/\s+/g, '').replace(/[–:]/g, '-');
    const suffix = regionFilter === 'sig' ? 'significant' : 'all';
    downloadCsv('GWAS_LW_intercept_region_' + label + '_' + suffix, buildCsv(list));
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
  }

  /* ------------------------------------------------------------------ *
   *  PER-VISIT WIRING                                                   *
   * ------------------------------------------------------------------ */
  function wireInteractions() {
    DOM = {
      plotCanvas: document.getElementById('gwxPlot'),
      miniCanvas: document.getElementById('gwxMinimap'),
      selrect: document.getElementById('gwxSelrect'),
      crosshair: document.getElementById('gwxCrosshair'),
      selLabelStart: document.getElementById('gwxSelLabelStart'),
      selLabelEnd: document.getElementById('gwxSelLabelEnd'),
      tooltip: document.getElementById('gwxTooltip'),
      readout: document.getElementById('gwxReadout'),
      chrJump: document.getElementById('gwxChrJump'),
      snpSearch: document.getElementById('gwxSearch'),
      searchStatus: document.getElementById('gwxSearchStatus'),
      modePan: document.getElementById('gwxModePan'),
      modeSelect: document.getElementById('gwxModeSelect'),
      panel: document.getElementById('gwxPanel'),
      panelBackdrop: document.getElementById('gwxPanelBackdrop'),
      regionCoord: document.getElementById('gwxRegionCoord'),
      regionSub: document.getElementById('gwxRegionSub'),
      statTotal: document.getElementById('gwxStatTotal'),
      statSig: document.getElementById('gwxStatSig'),
      tableBody: document.getElementById('gwxTableBody'),
      emptyState: document.getElementById('gwxEmpty'),
      shownHint: document.getElementById('gwxShownHint'),
      filterAll: document.getElementById('gwxFilterAll'),
      filterSig: document.getElementById('gwxFilterSig'),
      exportRegionBtn: document.getElementById('gwxExportRegionBtn'),
      sendVersityBtn: document.getElementById('gwxSendVersityBtn'),
      hoHint: document.getElementById('gwxHoHint'),
    };
    DOM.plotCtx = DOM.plotCanvas.getContext('2d');
    DOM.miniCtx = DOM.miniCanvas.getContext('2d');

    setMode(mode);

    DOM.plotCanvas.addEventListener('mousedown', onPlotMouseDown);
    DOM.plotCanvas.addEventListener('wheel', onPlotWheel, { passive: false });
    DOM.plotCanvas.addEventListener('dblclick', function () { setView(0, DATA.TOTAL_CUM); });
    DOM.miniCanvas.addEventListener('mousedown', function (e) { miniDragging = true; jumpFromMinimap(e); });

    document.getElementById('gwxZoomInBtn').addEventListener('click', function () { zoomBy(1 / 1.6); });
    document.getElementById('gwxZoomOutBtn').addEventListener('click', function () { zoomBy(1.6); });
    document.getElementById('gwxPanLeftBtn').addEventListener('click', function () { panBy(-0.4); });
    document.getElementById('gwxPanRightBtn').addEventListener('click', function () { panBy(0.4); });
    document.getElementById('gwxResetBtn').addEventListener('click', function () { DOM.chrJump.value = ''; setView(0, DATA.TOTAL_CUM); });
    DOM.chrJump.addEventListener('change', function () {
      const c = parseInt(DOM.chrJump.value, 10);
      if (!c) { setView(0, DATA.TOTAL_CUM); return; }
      setView(DATA.chrOffset[c], DATA.chrOffset[c] + DATA.chrLen[c]);
    });
    DOM.modePan.addEventListener('click', function () { setMode('pan'); });
    DOM.modeSelect.addEventListener('click', function () { setMode('select'); });
    DOM.snpSearch.addEventListener('keydown', onSearchKeydown);
    DOM.snpSearch.addEventListener('input', clearSearchStatus);

    DOM.filterAll.addEventListener('click', function () {
      regionFilter = 'all'; DOM.filterAll.classList.add('on'); DOM.filterSig.classList.remove('on'); renderRegionPanel();
    });
    DOM.filterSig.addEventListener('click', function () {
      regionFilter = 'sig'; DOM.filterSig.classList.add('on'); DOM.filterAll.classList.remove('on'); renderRegionPanel();
    });
    document.getElementById('gwxPanelCloseBtn').addEventListener('click', closePanel);
    DOM.panelBackdrop.addEventListener('click', closePanel);
    DOM.exportRegionBtn.addEventListener('click', exportRegion);
    DOM.sendVersityBtn.addEventListener('click', sendRegionToVersity);
    document.getElementById('gwxExportAllBtn').addEventListener('click', exportAllSignificant);

    wireGlobalOnce();
  }

  /* ------------------------------------------------------------------ *
   *  MARKUP                                                             *
   * ------------------------------------------------------------------ */
  function loadingHTML() {
    return '<div class="sec"><div class="bar"></div><div><h2 style="font-size:16px">Loading GWAS results…</h2>' +
      '<p>Fetching ' + escAttr(CFG.csvUrl) + ' — 60,000+ SNPs; this can take a moment on a slow connection.</p></div></div>';
  }
  function errorHTML(err) {
    return '<div class="card pad" style="border-color:#f3c2bd;background:#fff6f5;margin-top:14px">' +
      '<h3 style="font-family:var(--disp);color:#b42318;margin:0 0 8px">GWAS data failed to load</h3>' +
      '<p style="margin:0 0 12px;color:var(--muted);font-size:13px">' + escAttr(String(err && err.message || err)) + '</p>' +
      '<button class="btn" id="gwxRetry">Retry</button></div>';
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
      '<div class="card gwx-minimap-wrap">' +
        '<p class="gwx-minimap-label">Genome overview</p>' +
        '<canvas id="gwxMinimap"></canvas>' +
      '</div>' +
      '<p class="gwx-footnote">Drag to pan, scroll to zoom, or hold <b>Shift</b> and drag to zoom into a region. Switch to <b>Select region</b> and drag across a peak to inspect its SNPs and send the region to SNPVersity.</p>';
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
      '<div class="gwx-panel-controls">' +
        '<div class="gwx-seg" role="group" aria-label="Filter">' +
          '<button id="gwxFilterAll" class="on" type="button">All SNPs</button>' +
          '<button id="gwxFilterSig" type="button">Significant only</button>' +
        '</div>' +
        '<span class="gwx-hint" id="gwxShownHint"></span>' +
      '</div>' +
      '<div class="gwx-panel-table-wrap">' +
        '<table class="gwx-table" id="gwxTable"><thead><tr>' +
          '<th>SNP</th><th>Chr</th><th>BP</th><th>A1/A2</th><th>Freq</th><th>Beta</th><th>SE</th><th>P</th><th>−log₁₀P</th>' +
        '</tr></thead><tbody id="gwxTableBody"></tbody></table>' +
        '<div class="gwx-empty" id="gwxEmpty" style="display:none">No SNPs match this filter in the selected region.</div>' +
      '</div>' +
      '<div class="gwx-panel-foot">' +
        '<div class="row"><button class="btn" id="gwxExportRegionBtn" type="button">Export shown SNPs (CSV)</button></div>' +
        '<div class="row"><button class="btn primary" id="gwxSendVersityBtn" type="button">Send region to SNPVersity →</button></div>' +
        '<div class="gwx-ho-hint" id="gwxHoHint"></div>' +
      '</div>' +
    '</aside>';
  }

  function styleCSS() {
    return '' +
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
    '.gwx-minimap-wrap{margin-top:14px;padding:10px 12px 12px}' +
    '.gwx-minimap-label{font-size:10.5px;color:var(--faint);text-transform:uppercase;letter-spacing:.06em;margin:0 0 6px;font-weight:700}' +
    '#gwxMinimap{display:block;width:100%;height:54px;cursor:pointer;border-radius:var(--rad-sm)}' +
    '.gwx-footnote{font-size:11.5px;color:var(--faint);text-align:center;margin-top:10px}' +
    '.gwx-footnote b{color:var(--muted)}' +
    '.gwx-panel-backdrop{position:fixed;inset:0;background:rgba(10,15,28,.38);z-index:20;display:none}' +
    '.gwx-panel-backdrop.open{display:block}' +
    '.gwx-panel{position:fixed;top:0;right:0;bottom:0;width:min(480px,92vw);background:#fff;border-left:1px solid var(--line);box-shadow:var(--shadow-lg);z-index:21;display:flex;flex-direction:column;transform:translateX(100%);transition:transform .18s ease}' +
    '.gwx-panel.open{transform:translateX(0)}' +
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
    '.gwx-panel-controls{display:flex;align-items:center;justify-content:space-between;gap:10px;padding:12px 18px;border-bottom:1px solid var(--line)}' +
    '.gwx-panel-table-wrap{flex:1;overflow:auto;padding:0 0 8px}' +
    '.gwx-table{width:100%;border-collapse:collapse;font-size:12px;font-family:var(--mono)}' +
    '.gwx-table th{position:sticky;top:0;background:#fff;text-align:left;font-family:var(--body);font-size:10px;text-transform:uppercase;letter-spacing:.05em;color:var(--faint);font-weight:700;padding:8px 10px;border-bottom:1px solid var(--line)}' +
    '.gwx-table td{padding:6px 10px;border-bottom:1px solid var(--line);font-variant-numeric:tabular-nums;color:var(--muted);white-space:nowrap}' +
    '.gwx-table td.snpid{color:var(--navy-700);font-weight:700}' +
    '.gwx-table tr.sig td.snpid::before{content:"";display:inline-block;width:6px;height:6px;border-radius:50%;background:#c0392b;margin-right:6px}' +
    '.gwx-panel-foot{padding:12px 18px 16px;border-top:1px solid var(--line);display:flex;flex-direction:column;gap:8px}' +
    '.gwx-panel-foot .row{display:flex;gap:8px}' +
    '.gwx-empty{padding:30px 18px;text-align:center;color:var(--faint);font-size:12.5px}' +
    '.gwx-hint{font-size:12px;color:var(--faint)}' +
    '.gwx-ho-hint{font-size:11px;color:var(--faint);min-height:14px}' +
    '@media (max-width:720px){ .gwx-readout{display:none} #gwxPlot{height:340px} }';
  }

  function shellHTML() {
    return '<style>' + styleCSS() + '</style>' +
    '<div class="sec"><div class="bar"></div><div style="width:100%">' +
      '<div class="n">GWAS MANHATTAN EXPLORER · GCTA MLMA</div>' +
      '<h2>Scan the genome for trait-associated SNPs</h2>' +
      '<p>Trait: <b>' + CFG.traitLabel + '</b> · ' + CFG.chrCount + ' NAM chromosomes, ' + DATA.n.toLocaleString() + ' SNPs · ' +
      'significance threshold <b>p &lt; ' + CFG.threshLabel + '</b> (SimpleM) — ' + DATA.totalSig.toLocaleString() + ' significant SNPs. ' +
      'Use "Select Region" mode to drag-select a region to inspect it and download SNPS or hand them off to SNPVersity. Note: SNPs with raw p values above 0.001 (-log10 P < 3) are not shown to save memory. </p>' +
    '</div></div>' +
    toolbarHTML() +
    plotHTML() +
    panelHTML();
  }

  /* ------------------------------------------------------------------ *
   *  RENDER — router entry point                                       *
   * ------------------------------------------------------------------ */
  function render(page) {
    page.className = 'page fade';
    if (loadError) {
      page.innerHTML = errorHTML(loadError);
      const retry = document.getElementById('gwxRetry');
      if (retry) retry.addEventListener('click', function () { loadError = null; render(page); });
      return;
    }
    if (!DATA) {
      page.innerHTML = loadingHTML();
      loadData().then(function () {
        if (S.tool === 'snpgwas') render(document.getElementById('page'));
      }).catch(function () {
        if (S.tool === 'snpgwas') render(document.getElementById('page'));
      });
      return;
    }
    page.innerHTML = shellHTML();
    wireInteractions();
    updateReadout();
    resizeCanvases();
  }

  SNPTools.register('snpgwas', { render: render });
})();
