/* =====================================================================
 *  paneffect-heatmap.js — high-performance heatmap renderer.
 *
 *  WHY THIS EXISTS
 *  The original renderHeatmap()/renderHeatmapPan() build ONE <div> per
 *  cell (positions × 20 subs, or alignment_length × N genomes — tens of
 *  thousands of nodes), attach THREE listeners to each, and append a
 *  separate tooltip <div> to <body> per cell. That is the actual cost the
 *  full-protein views pay, and the per-cell tooltips leak on every slider
 *  move because innerHTML='' clears the cells but not the body tooltips.
 *
 *  This module draws the whole grid to a single <canvas> in one pass and
 *  uses ONE shared tooltip resolved by hit-testing. ~1 DOM node instead of
 *  ~50k, no per-cell listeners, no leak. It is framework-agnostic: pass in
 *  a color function (their colorScale / colorScalePan work as-is).
 *
 *  It returns a controller with .highlight(), .redraw() and .destroy() so a
 *  panel can clean up completely when the user navigates to another tool.
 * ===================================================================== */
(function (global) {
  'use strict';

  /**
   * renderGrid(container, cells, opts) -> controller
   *
   * container : DOM element to render into (canvas is appended here)
   * cells     : array of { x, y, score, ...meta }
   *               x = 1-based column (protein position / MSA column)
   *               y = 1-based row    (1..20 substitution, or genome index)
   * opts:
   *   cols        : total columns (defaults to max x in cells)
   *   rows        : total rows    (defaults to max y in cells)
   *   cellW,cellH : pixel size per cell (default 3 × 14)
   *   color       : (score) => cssColor   (pass colorScale or colorScalePan)
   *   tooltip     : (cell)  => htmlString (optional; default shows x/y/score)
   *   highlight   : { x, y } to ring on first draw (optional)
   *   dpr         : override devicePixelRatio (testing)
   */
  function renderGrid(container, cells, opts) {
    opts = opts || {};
    const color   = opts.color   || defaultColor;
    const tipHtml = opts.tooltip || defaultTooltip;
    const cellW   = opts.cellW || 3;
    const cellH   = opts.cellH || 14;

    let cols = opts.cols, rows = opts.rows;
    if (!cols || !rows) {
      let mx = 0, my = 0;
      for (const c of cells) { if (c.x > mx) mx = c.x; if (c.y > my) my = c.y; }
      cols = cols || mx; rows = rows || my;
    }

    const W = cols * cellW, H = rows * cellH;
    const dpr = opts.dpr || (global.devicePixelRatio || 1);

    /* build a fast lookup so hit-testing is O(1): key = x*10000 + y */
    const index = new Map();
    for (const c of cells) index.set(c.x * 100000 + c.y, c);

    const canvas = document.createElement('canvas');
    canvas.width  = Math.max(1, Math.round(W * dpr));
    canvas.height = Math.max(1, Math.round(H * dpr));
    canvas.style.width  = W + 'px';
    canvas.style.height = H + 'px';
    canvas.style.display = 'block';
    container.appendChild(canvas);

    const ctx = canvas.getContext('2d');
    ctx.scale(dpr, dpr);

    let highlight = opts.highlight || null;

    function draw() {
      ctx.clearRect(0, 0, W, H);
      for (const c of cells) {
        ctx.fillStyle = color(c.score, c);
        /* +0.5 width mirrors the original's slight overlap so there are no seams */
        ctx.fillRect((c.x - 1) * cellW, (c.y - 1) * cellH, cellW + 0.5, cellH);
      }
      if (highlight && highlight.y != null) {
        const hx = (highlight.x - 1) * cellW, hy = (highlight.y - 1) * cellH;
        ctx.lineWidth = 2;
        ctx.strokeStyle = '#111';
        ctx.strokeRect(hx - 1, hy - 1, cellW + 2, cellH + 2);
        ctx.strokeStyle = '#fff';
        ctx.strokeRect(hx - 2.5, hy - 2.5, cellW + 5, cellH + 5);
      } else if (highlight && highlight.x != null) {
        /* position known, but not which row (e.g. the amino-acid substitution
           couldn't be resolved) — still mark the column so the jump isn't silent */
        const hx = (highlight.x - 1) * cellW;
        ctx.lineWidth = 2;
        ctx.strokeStyle = '#111';
        ctx.strokeRect(hx - 1, -1, cellW + 2, H + 2);
        ctx.strokeStyle = '#fff';
        ctx.strokeRect(hx - 2.5, -2.5, cellW + 5, H + 5);
      }
    }
    draw();

    /* one shared tooltip for the whole grid */
    const tip = document.createElement('div');
    tip.className = 'pe-canvas-tip';
    tip.style.cssText =
      'position:fixed;z-index:9999;pointer-events:none;visibility:hidden;' +
      'background:rgba(0,0,0,.78);color:#fff;padding:5px 9px;border-radius:6px;' +
      'font:12px/1.35 Arial,sans-serif;white-space:nowrap';
    document.body.appendChild(tip);

    function cellAt(evt) {
      const r = canvas.getBoundingClientRect();
      const px = evt.clientX - r.left, py = evt.clientY - r.top;
      if (px < 0 || py < 0 || px >= W || py >= H) return null;
      const x = Math.floor(px / cellW) + 1;
      const y = Math.floor(py / cellH) + 1;
      return index.get(x * 100000 + y) || null;
    }

    function onMove(e) {
      const c = cellAt(e);
      if (!c) { tip.style.visibility = 'hidden'; return; }
      tip.innerHTML = tipHtml(c);
      tip.style.left = (e.clientX + 12) + 'px';
      tip.style.top  = (e.clientY + 12) + 'px';
      tip.style.visibility = 'visible';
    }
    function onLeave() { tip.style.visibility = 'hidden'; }

    canvas.addEventListener('mousemove', onMove);
    canvas.addEventListener('mouseleave', onLeave);

    return {
      canvas,
      redraw: draw,
      highlight(pos) { highlight = pos; draw(); },
      /* full teardown — call this when the panel is torn down */
      destroy() {
        canvas.removeEventListener('mousemove', onMove);
        canvas.removeEventListener('mouseleave', onLeave);
        if (tip.parentNode) tip.parentNode.removeChild(tip);
        if (canvas.parentNode) canvas.parentNode.removeChild(canvas);
        index.clear();
      },
    };
  }

  function defaultColor(score) {
    if (score == null || Number.isNaN(score)) return '#dddddd';
    return score > -7 ? '#3860aa' : '#db4c4d';
  }
  function defaultTooltip(c) {
    return `Position: ${c.x}<br>Row: ${c.y}<br>Score: ${c.score}`;
  }

  global.PanEffectHeatmap = { renderGrid };
})(typeof window !== 'undefined' ? window : globalThis);
