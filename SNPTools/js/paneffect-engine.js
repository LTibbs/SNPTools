/* =====================================================================
 *  paneffect-engine.js — native PanEffect engine for SNPTools.
 *
 *  PanEffect's scripts (support.js, genome.js, pan.js, dom.js, main.js) are
 *  classic scripts, so they all share ONE global scope. This engine (also a
 *  classic script) therefore reads and writes their globals directly
 *  (gene_model, main_option, gene_model_length, GN_array, …) and reassigns
 *  renderHeatmap / renderHeatmapPan to canvas versions.
 *
 *  Load order (see index.html):
 *    d3  ->  support.js -> genome.js -> pan.js -> dom.js -> main.js
 *        ->  paneffect-heatmap.js -> paneffect-engine.js -> snppaneffect.js
 *
 *  Exposes window.PanEffectEngine.render(container, opts) / .teardown().
 *    opts = { gene, esm, option, variant:{pos,wt,sub}, dataBase }
 * ===================================================================== */
(function () {
  'use strict';
  if (typeof window === 'undefined') return;

  var DATA_DIRS = ['csv','heatmap','target','query','pfam','uniprot','synonym','traits','dssp','structures'];
  var ESM_OK = ['ESM1','ESM2','ESM3'];

  var state = {
    installedFetch: false, origFetch: null,
    base: './paneffect/', variant: null,
    grids: { b73: null, pan: null },
    /* last rows handed to each canvas renderer, so a colour-scheme change
       can re-render without re-running the whole fetch pipeline */
    lastData: { b73: null, pan: null },
    /* genome.js's original updateHeatmapZoom, captured once before we wrap it */
    origUpdateHeatmapZoom: null,
  };

  /* ---------------- fetch rebaser -------------------------------------
     Rewrites only relative paths whose first segment is a known PanEffect
     data directory, so SNPTools' own requests are never touched. */
  function normalizeBase(b) {
    b = String(b || './paneffect/');
    if (b.slice(-1) !== '/') b += '/';
    return b;
  }
  function rebase(url) {
    if (typeof url !== 'string') return url;
    if (/^[a-z]+:\/\//i.test(url) || url.indexOf('//') === 0) return url; // absolute
    var u = url.replace(/^\.\//, '').replace(/^\//, '');
    var seg = u.split('/')[0].split('?')[0];
    if (DATA_DIRS.indexOf(seg) === -1) return url;
    return state.base + u;
  }
  function installFetch() {
    if (state.installedFetch) return;
    var orig = window.fetch.bind(window);
    state.origFetch = orig;
    window.fetch = function (input, init) {
      try {
        var url = (typeof input === 'string') ? input : (input && input.url);
        if (typeof url === 'string') {
          var r = rebase(url);
          if (r !== url) input = (typeof input === 'string') ? r : new Request(r, input);
        }
      } catch (e) { /* fall through to original */ }
      return orig(input, init);
    };
    state.installedFetch = true;
  }

  /* ---------------- id parsing (writes shared globals) --------------- */
  function parseId(id) {
    id = String(id || '').trim();
    if (id.indexOf('_T') >= 0)      { gene_model = id.split('_')[0]; transcript = id;                   protein = id.replace('_T', '_P'); can_flag = false; }
    else if (id.indexOf('_P') >= 0) { gene_model = id.split('_')[0]; transcript = id.replace('_P','_T'); protein = id;                     can_flag = false; }
    else                            { gene_model = id; can_flag = true; transcript = id + '_T001';        protein = id + '_P001'; }
  }

  /* ---------------- canvas heatmap renderers (override) -------------- */
  var Y_ROWS = 20;
  function subToRow(sub) {
    if (!sub) return null;
    sub = String(sub).toUpperCase();
    for (var y = 1; y <= Y_ROWS; y++) { if (numberToAminoAcid(y) === sub) return y; }
    return null;
  }

  function canvasRenderGene(data) {
    state.lastData.b73 = data;
    var host = document.getElementById('full-heatmap');
    if (!host || !window.PanEffectHeatmap) return;
    if (state.grids.b73) { state.grids.b73.destroy(); state.grids.b73 = null; }
    Array.prototype.forEach.call(host.querySelectorAll('canvas.pe-hm'), function (c) { c.remove(); });

    var cells = data.map(function (d) {
      return { x: +d.X, y: +d.Y, score: +d.Score, wt: d.WT, wgs2024: +d.WGS2024, wgs2026: +d.WGS2026 };
    });
    var vis = function (c) {
      return wgs_status || (wgs2024_status && c.wgs2024 === 1) || (wgs2026_status && c.wgs2026 === 1);
    };
    var hl = state.variant && state.variant.pos ? { x: state.variant.pos, y: subToRow(state.variant.sub) } : null;

    var grid = PanEffectHeatmap.renderGrid(host, cells, {
      cols: gene_model_length, rows: Y_ROWS,
      cellW: window_length / gene_model_length, cellH: 14,
      color: function (score, c) { return vis(c) ? THEME.colorFor('gene', score) : '#FFFFFF'; },
      tooltip: function (c) {
        return 'Position: ' + c.x + '<br>Substitution: ' + c.wt + ' &rarr; ' +
               numberToAminoAcid(c.y) + '<br>Score: ' + c.score;
      },
      highlight: hl,
    });
    grid.canvas.classList.add('pe-hm');
    state.grids.b73 = grid;
    legend(host, '.heatmap-container', 'pe-legend-gene');
  }

  function canvasRenderPan(data) {
    state.lastData.pan = data;
    var host = document.getElementById('full-heatmap-pan');
    if (!host || !window.PanEffectHeatmap) return;
    if (state.grids.pan) { state.grids.pan.destroy(); state.grids.pan = null; }
    Array.prototype.forEach.call(host.querySelectorAll('canvas.pe-hm'), function (c) { c.remove(); });

    var rows = GN_size || 0;
    var cells = data.map(function (d) {
      var wt = String(d.WT == null ? '' : d.WT).trim();
      var sub = String(d.Sub == null ? '' : d.Sub).trim();
      var raw = d.Score;
      var sc = (raw === '' || raw == null) ? NaN : +raw;
      var gap = wt.charAt(0) === '-' || sub.charAt(0) === '-' || raw === '' || raw == null || isNaN(sc);
      if (!rows && +d.Y > rows) rows = +d.Y;
      return { x: +d.X, y: +d.Y, score: sc, gap: gap, wt: d.WT, sub: d.Sub, x2: d.X2, x3: d.X3 };
    });

    var grid = PanEffectHeatmap.renderGrid(host, cells, {
      cols: alignment_length, rows: rows || 1,
      cellW: window_length / alignment_length, cellH: 5,
      color: function (score, c) { return c.gap ? THEME.gapColor() : THEME.colorFor('pan', score); },
      tooltip: function (c) {
        return 'B73 Position: ' + (c.x2 == null ? '' : c.x2) +
               '<br>Target Position: ' + (c.x3 == null ? '' : c.x3) +
               '<br>Genome: ' + (GN_array[c.y] == null ? '' : GN_array[c.y]) +
               '<br>G.M.: ' + (GM_array[c.y] == null ? '' : GM_array[c.y]) +
               '<br>Substitution: ' + c.wt + ' to ' + c.sub +
               '<br>Score: ' + (isNaN(c.score) ? '' : c.score);
      },
    });
    grid.canvas.classList.add('pe-hm');
    state.grids.pan = grid;
    var box = document.getElementById('colorBox-pan');
    if (box) { box.innerHTML = ''; box.appendChild(buildLegendBar()); }
  }

  /* ---------------- heterotic-group key (pan view) ------------------
     The pan-genome rows are one genome/accession each, tinted by
     heterotic group (support.js colorGenome). The canvas re-render drops
     the DOM key that pan.js's renderHeatmapPan used to append, so we
     rebuild it here. Colours mirror colorGenome() exactly — keep the two
     in sync if either changes. */
  var HETEROTIC_GROUPS = [
    { label: 'Stiff stalk',     color: 'black'   },
    { label: 'Mix',             color: '#666666' },
    { label: 'Non-stiff-stalk', color: '#455edd' },
    { label: 'Iodent',          color: '#a807ed' },
    { label: 'Lancaster',       color: '#da9af5' },
    { label: 'European flint',  color: '#9ae6f5' },
    { label: 'Chinese',         color: '#f50707' },
    { label: 'Tang SiPingTou',  color: '#fa7d7d' },
    { label: 'Popcorn',         color: '#ce58ce' },
    { label: 'Sweet corn',      color: 'pink'    },
    { label: 'Tropical',        color: '#30c727' },
    { label: 'PanAnd',          color: '#773510' },
    { label: 'Teosinte',        color: '#ca854c' },
    { label: 'Hi/Lo',           color: '#1B9E77' }
  ];
  function buildHeteroticLegend() {
    var wrap = document.createElement('div');
    wrap.style.display = 'inline-block';
    wrap.style.border = '1px solid #000';
    wrap.style.borderRadius = '5px';
    wrap.style.padding = '10px 14px';
    wrap.style.margin = '12px 0';
    wrap.style.fontFamily = 'Arial, sans-serif';
    wrap.style.fontSize = '13px';
    wrap.style.lineHeight = '1.55';
    wrap.style.whiteSpace = 'nowrap';   // let the box widen instead of wrapping names

    var title = document.createElement('div');
    title.innerText = 'Heterotic group';
    title.style.fontWeight = '600';
    title.style.textDecoration = 'underline';
    title.style.marginBottom = '4px';
    wrap.appendChild(title);

    HETEROTIC_GROUPS.forEach(function (g) {
      var row = document.createElement('div');
      row.style.display = 'flex';
      row.style.alignItems = 'center';

      var sw = document.createElement('span');
      sw.style.display = 'inline-block';
      sw.style.width = '12px';
      sw.style.height = '12px';
      sw.style.marginRight = '6px';
      sw.style.borderRadius = '2px';
      sw.style.background = g.color;
      sw.style.border = '1px solid rgba(0,0,0,.25)';
      sw.style.flex = '0 0 auto';
      row.appendChild(sw);

      var txt = document.createElement('span');
      txt.innerText = g.label;
      txt.style.color = g.color;
      txt.style.whiteSpace = 'nowrap';
      row.appendChild(txt);

      wrap.appendChild(row);
    });
    return wrap;
  }
  function populateHeteroticLegend() {
    var box = document.getElementById('heterotic-legend-pan');
    if (!box) return;
    /* Park the key to the right of the zoomed pan view's genome-name
       labels — those right-hand labels sit at ~x=1206px (window_length =
       1200) and are tinted by colorGenome(), so a genome's label colour
       can be read straight across to the named group here. Positioned
       relative to .heatmap-container-zoom-pan (position:relative). */
    box.style.position = 'absolute';
    box.style.top = '4px';
    box.style.left = (window_length + 220) + 'px';
    box.innerHTML = '';
    box.appendChild(buildHeteroticLegend());
  }

  /* ---------------- zoomed-region highlight overlay -------------------
     updateHeatmapZoom() (genome.js) builds the zoomed div-per-cell grid
     directly, clearing and re-populating #zoomed-heatmap on every call —
     it isn't ours to rewrite, so we wrap it instead (same technique as
     the canvas swap above): call the original to draw the cells, then
     drop a bordered overlay <div> on top of whichever cell matches
     state.variant, if that position falls inside the current 50-residue
     window. Falls back to marking the whole column when the amino-acid
     row can't be resolved, same as the full-protein canvas ring. */
  function zoomHighlightOverlay(start) {
    var host = document.getElementById('zoomed-heatmap');
    if (!host) return;
    if (!state.variant || state.variant.pos == null) return;
    var x = state.variant.pos;
    var end = start + 49;
    if (x < start || x > end) return; // variant isn't in the visible window

    var cellWidth = window_length / 50;
    var cellHeight = 20;
    var yRow = subToRow(state.variant.sub);

    var ring = document.createElement('div');
    ring.className = 'pe-zoom-ring';
    ring.style.position = 'absolute';
    ring.style.pointerEvents = 'none';
    ring.style.boxSizing = 'border-box';
    ring.style.border = '2px solid #111';
    ring.style.boxShadow = '0 0 0 2px #fff';
    ring.style.left = ((x - start) * cellWidth - 1) + 'px';
    ring.style.width = (cellWidth + 2) + 'px';
    if (yRow != null) {
      ring.style.top = (cellHeight * (yRow - 1) - 1) + 'px';
      ring.style.height = (cellHeight + 2) + 'px';
    } else {
      ring.style.top = '-1px';
      ring.style.height = (cellHeight * Y_ROWS + 2) + 'px';
    }
    host.appendChild(ring);
  }

  /* Install the wrap exactly once: capture the true genome.js original
     before ever reassigning the global, so re-navigating (render() called
     again for a different gene) re-points to the same wrapper rather than
     wrapping an already-wrapped function. */
  function installZoomHighlight() {
    if (state.origUpdateHeatmapZoom || typeof updateHeatmapZoom !== 'function') return;
    state.origUpdateHeatmapZoom = updateHeatmapZoom;
    updateHeatmapZoom = function (data, start) {
      state.origUpdateHeatmapZoom(data, start);
      try { zoomHighlightOverlay(start); } catch (e) { console.warn('[PanEffectEngine] zoom highlight failed:', e); }
    };
  }

  /* ================= colour theme (inlined) ==========================
     Deliberately inlined rather than kept in a separate file: a missing
     or 404'd <script> tag silently degraded to a no-op install(), which
     recoloured the canvas heatmaps while leaving the zoomed views on the
     old scheme and the toggle inert. One file, no load-order to get wrong.

     Two schemes:
       'classic'  — delegates to the real colorScale / colorScalePan from
                    support.js, so it is byte-identical to the original.
       'snptools' — SNPImpact's ramp from snpimpact.js scoreColor().

     Sentinels (support.js) — NOT ramp values, and both schemes honour them:
                  undefined / NaN     score === 0
       gene       white               black   -> wild-type diagonal
       pan        #DDDDDD             #00429d -> exact match to B73
     ================================================================== */
  var THEME = (function () {
    var STORE_KEY = 'snptools.pe.scheme';

    /* Confirmed against support.js: thresholds run from >0 down to <=-22
       in steps of 2, so the scored range is exactly -22..0. */
    var DOMAINS = {
      impact: { lo: -12, hi: 6 },
      gene:   { lo: -22, hi: 0 },
      pan:    { lo: -22, hi: 0 },
    };

    /* GAP matches the '#e6e6e6' hardcoded in pan.js so the canvas and
       zoomed pan views agree. */
    var GAP = '#e6e6e6';
    var NEUTRAL = {
      classic:  { gene: { blank: 'white',   zero: 'black'   },
                  pan:  { blank: '#DDDDDD', zero: '#00429d' } },
      snptools: { gene: { blank: '#ffffff', zero: '#13264a' },
                  pan:  { blank: '#DDDDDD', zero: '#2563eb' } },
    };

    var CLASSIC_COLORS = ['#00429d','#3860aa','#587fb3','#78a0b7','#9ac0b3','#c1e19e',
      '#ffff00','#ffd337','#fea447','#f1784d','#db4c4d','#bd2147','#93003a'];

    var ORIGINAL = { gene: null, pan: null };
    var scheme = 'snptools';

    function clamp01(t) { return t < 0 ? 0 : (t > 1 ? 1 : t); }
    function snpRGB(t) {
      t = clamp01(t);
      var r = t < 0.5 ? 255 : Math.round(255 * (1 - (t - 0.5) * 2));
      var g = t < 0.5 ? Math.round(255 * t * 2) : 200;
      return [r, Math.max(60, g), 60];
    }
    function snpRamp(t) { var c = snpRGB(t); return 'rgb(' + c[0] + ',' + c[1] + ',' + c[2] + ')'; }
    function snpHex(t) {
      var c = snpRGB(t), s = '#', i;
      for (i = 0; i < 3; i++) s += ('0' + c[i].toString(16)).slice(-2);
      return s;
    }
    function classicRamp(t) {
      return CLASSIC_COLORS[Math.round(clamp01(1 - t) * (CLASSIC_COLORS.length - 1))];
    }
    function tFor(kind, score) {
      var d = DOMAINS[kind] || DOMAINS.gene;
      return clamp01((+score - d.lo) / (d.hi - d.lo));
    }

    function colorFor(kind, score) {
      kind = (kind === 'pan') ? 'pan' : (kind === 'impact' ? 'impact' : 'gene');
      var nk = (kind === 'impact') ? 'gene' : kind;

      /* Classic hands the raw value straight to support.js so every input —
         including degenerate ones like NaN — behaves exactly as it always
         has. Checking sentinels first would diverge on the edges. */
      if (scheme === 'classic') {
        var orig = ORIGINAL[nk];
        if (typeof orig === 'function') {
          try { return orig(score); } catch (e) { /* fall through */ }
        }
        var nc = NEUTRAL.classic[nk];
        if (score === undefined || score === null || score === '' || isNaN(+score)) return nc.blank;
        if (+score === 0) return nc.zero;
        return classicRamp(tFor(kind, +score));
      }

      var n = NEUTRAL.snptools[nk];
      if (score === undefined || score === null || score === '' || isNaN(+score)) return n.blank;
      var v = +score;
      if (v === 0) return n.zero;
      return snpRamp(tFor(kind, v));
    }

    function stops(reverse) {
      var out = [], i;
      if (scheme === 'classic') out = CLASSIC_COLORS.slice().reverse();
      else for (i = 0; i < 9; i++) out.push(snpHex(i / 8));
      return reverse ? out.reverse() : out;
    }
    function gradientCSS(deg, reverse) {
      var s = stops(reverse), i;
      deg = (deg == null) ? 90 : deg;
      if (scheme === 'classic') {
        var parts = [], w = 100 / s.length;
        for (i = 0; i < s.length; i++) {
          parts.push(s[i] + ' ' + (i * w).toFixed(2) + '% ' + ((i + 1) * w).toFixed(2) + '%');
        }
        return 'linear-gradient(' + deg + 'deg,' + parts.join(',') + ')';
      }
      return 'linear-gradient(' + deg + 'deg,' + s.join(',') + ')';
    }

    /* support.js declares colorScale / colorScalePan as top-level function
       declarations, so they are writable properties of the global object.
       genome.js and pan.js call them BY NAME at draw time, which is why
       re-pointing them here reaches the zoomed views too. */
    function install() {
      if (!ORIGINAL.gene && typeof window.colorScale === 'function' && !window.colorScale.__peTheme) {
        ORIGINAL.gene = window.colorScale;
      }
      if (!ORIGINAL.pan && typeof window.colorScalePan === 'function' && !window.colorScalePan.__peTheme) {
        ORIGINAL.pan = window.colorScalePan;
      }
      var gene = function (score) { return colorFor('gene', score); };
      var pan  = function (score) { return colorFor('pan',  score); };
      gene.__peTheme = pan.__peTheme = true;

      window.colorScale = gene;
      window.colorScalePan = pan;
      try { colorScale = gene; colorScalePan = pan; } catch (e) { /* window write counts */ }

      return window.colorScale === gene;
    }

    function getScheme() { return scheme; }
    function setScheme(next) {
      next = (next === 'classic') ? 'classic' : 'snptools';
      scheme = next;
      try { window.localStorage.setItem(STORE_KEY, scheme); } catch (e) {}
      install();
      return scheme;
    }
    function restore() {
      try {
        var v = window.localStorage.getItem(STORE_KEY);
        if (v === 'classic' || v === 'snptools') scheme = v;
      } catch (e) {}
      return scheme;
    }
    function ready() {
      return !!(ORIGINAL.gene && ORIGINAL.pan &&
                window.colorScale && window.colorScale.__peTheme);
    }

    return {
      DOMAINS: DOMAINS, CLASSIC_COLORS: CLASSIC_COLORS,
      colorFor: colorFor, gapColor: function () { return GAP; },
      stops: stops, gradientCSS: gradientCSS,
      getScheme: getScheme, setScheme: setScheme, restore: restore,
      install: install, ready: ready,
    };
  })();
  window.PanEffectTheme = THEME;



  /* compact horizontal score legend — benign (green) on the left,
     strong effect (red) on the right, matching SNPImpact's score bar. */
  function buildLegendBar() {
    var wrap = document.createElement('div');
    wrap.className = 'pe-legend';
    var bar = document.createElement('div');
    bar.className = 'pe-legend-bar';
    /* reversed: the ramp runs red -> green, the legend reads benign -> strong */
    bar.style.background = THEME.gradientCSS(90, true);
    var labels = document.createElement('div'); labels.className = 'pe-legend-labels';
    labels.innerHTML = '<span>&gt; 0 (benign)</span><span>&minus;11</span><span>&lt; &minus;22 (strong)</span>';
    /* score 0 and no-data are sentinels in support.js, not ramp values, so
       they get their own swatches rather than being read off the gradient */
    var keys = document.createElement('div'); keys.className = 'pe-legend-keys';
    keys.innerHTML =
      '<span><i style="background:' + THEME.colorFor('gene', 0) + '"></i>wild type</span>' +
      '<span><i style="background:' + THEME.colorFor('pan', 0) + '"></i>exact match</span>' +
      '<span><i style="background:' + THEME.gapColor() + '"></i>gap / no data</span>';
    wrap.appendChild(bar); wrap.appendChild(labels); wrap.appendChild(keys);
    return wrap;
  }
  function legend(hostEl, containerSel, cls) {
    var container = hostEl.closest ? hostEl.closest(containerSel) : null;
    if (!container || !container.parentNode) return;
    var prev = container.parentNode.querySelector('.' + cls);
    if (prev) prev.remove();
    var el = buildLegendBar(); el.classList.add(cls);
    container.parentNode.insertBefore(el, container.nextSibling);
  }

  /* ---------------- colour-scheme toggle -----------------------------
     One radio pair per view (B73 and pan). Both pairs drive the same
     PanEffectTheme state and stay mirrored, so switching in one view is
     already applied when you cross to the other. */
  function schemeRow(sfx) {
    var name = 'peScheme' + sfx;
    return '' +
      '<div class="pe-scheme" id="pe-scheme' + sfx + '">' +
        '<span class="pe-scheme-lab">Colour scheme</span>' +
        '<label><input type="radio" name="' + name + '" value="snptools" checked> SNPTools</label>' +
        '<label><input type="radio" name="' + name + '" value="classic"> Classic</label>' +
      '</div>';
  }

  function syncSchemeInputs(value) {
    Array.prototype.forEach.call(
      document.querySelectorAll('.pe-scheme input[type=radio]'),
      function (r) { r.checked = (r.value === value); }
    );
  }

  function wireScheme() {
    syncSchemeInputs(THEME.getScheme());
    Array.prototype.forEach.call(
      document.querySelectorAll('.pe-scheme input[type=radio]'),
      function (r) {
        r.addEventListener('change', function () {
          if (!r.checked) return;
          THEME.setScheme(r.value);
          syncSchemeInputs(r.value);
          refreshColors();
        });
      }
    );
  }

  /* Repaint everything the scheme touches.

     The zoomed views need no colour logic of their own: updateHeatmapZoom
     (genome.js) and updateHeatmapZoomPan (pan.js) build plain <div> cells
     and call the GLOBAL colorScale / colorScalePan by name at draw time.
     Both are top-level function declarations despite being indented, so
     they are callable from here — we just re-run them with the same rows
     and the current slider position, and they pick up the new scheme.
     Dispatching 'input' on the sliders is kept as a second path in case a
     build wires the redraw differently. */
  function sliderStart(id) {
    var s = document.getElementById(id);
    var v = s ? parseInt(s.value, 10) : 1;
    return (isNaN(v) || v < 1) ? 1 : v;
  }

  function redrawZoomViews() {
    var done = { gene: false, pan: false };
    try {
      if (typeof updateHeatmapZoom === 'function' && state.lastData.b73) {
        updateHeatmapZoom(state.lastData.b73, sliderStart('zoom-slider'));
        done.gene = true;
      }
    } catch (e) { console.warn('[PanEffectEngine] zoom redraw (gene) failed:', e); }
    try {
      if (typeof updateHeatmapZoomPan === 'function' && state.lastData.pan) {
        updateHeatmapZoomPan(state.lastData.pan, sliderStart('zoom-slider-pan'));
        done.pan = true;
      }
    } catch (e) { console.warn('[PanEffectEngine] zoom redraw (pan) failed:', e); }

    /* fallback: let PanEffect redraw them the way centerOnVariant does */
    if (!done.gene || !done.pan) {
      ['zoom-slider', 'zoom-slider-pan'].forEach(function (id) {
        var s = document.getElementById(id);
        if (s) s.dispatchEvent(new Event('input'));
      });
    }
    return done;
  }

  function refreshColors() {
    /* canvas grids: re-render from the cached rows */
    if (state.lastData.b73) { try { canvasRenderGene(state.lastData.b73); } catch (e) {} }
    if (state.lastData.pan) { try { canvasRenderPan(state.lastData.pan); } catch (e) {} }
    redrawZoomViews();
    rebuildLegends();
  }

  function rebuildLegends() {
    Array.prototype.forEach.call(document.querySelectorAll('.pe-legend-bar'), function (bar) {
      bar.style.background = THEME.gradientCSS(90, true);
    });
    var box = document.getElementById('colorBox-pan');
    if (box) { box.innerHTML = ''; box.appendChild(buildLegendBar()); }
  }

  /* ---------------- view visibility --------------------------------- */
  function applyVisibility(option) {
    var sm = document.getElementById('summary');
    var b = document.getElementById('b73');
    var p = document.getElementById('pan-genome');
    if (sm) sm.classList.add('active');
    if (b) b.classList.toggle('active', option !== 'pan');
    if (p) p.classList.toggle('active', option !== 'b73');
  }

  function centerOnVariant(pos) {
    var s = document.getElementById('zoom-slider');
    if (!s) return;
    var start = Math.max(1, Math.round(pos - 25));
    var max = +s.max || start;
    s.value = Math.min(max, start);
    s.dispatchEvent(new Event('input'));
  }

  /* ---------------- WGS view default -------------------------------- */
  function applyWgs(mode) {
    var is2026 = (mode === 'maize2026');
    wgs_status = !is2026;         /* 'all variants' unless a 2026 handoff */
    wgs2024_status = false;
    wgs2026_status = is2026;
    var rAll = document.getElementById('allVariants');
    var r26 = document.getElementById('maizeWGS2026');
    if (r26) r26.checked = is2026;
    if (rAll) rAll.checked = !is2026;
  }

  /* ---------------- clean summary (override) ------------------------ */
  function cleanSummary() {
    var gnPrint = (typeof gn === 'string' && gn.trim() && gn.trim() !== 'N/A') ? ' (' + gn.trim() + ')' : '';
    var file = (typeof gene_model_file !== 'undefined' && gene_model_file)
      ? gene_model_file : ('./csv/' + currentESM + '/' + protein + '.csv');
    var csv = state.base + String(file).replace(/^\.\//, '');
    var sm = document.getElementById('summary');
    if (sm) {
      sm.innerHTML =
        '<div class="pe-sum">' +
          '<div class="pe-sum-kicker">Missense variant effects · ' + currentESM + ' · B73 v5</div>' +
          '<h2 class="pe-sum-title">' + (gene_model || '') + gnPrint + '</h2>' +
          '<a class="pe-sum-dl" href="' + csv + '" download>' +
            '<svg viewBox="0 0 24 24" width="16" height="16" fill="none" aria-hidden="true">' +
            '<path d="M12 4v10m0 0l-4-4m4 4l4-4M5 20h14" stroke="currentColor" stroke-width="2" ' +
            'stroke-linecap="round" stroke-linejoin="round"/></svg>' +
            'Download variant effects file</a>' +
        '</div>';
    }
    var r = document.getElementById('reference_gm'); if (r) r.innerHTML = 'B73 Reference View';
    var p = document.getElementById('pan_gm'); if (p) p.innerHTML = 'Pan-genome View';
  }

  /* ---------------- teardown ---------------------------------------- */
  function teardown() {
    state.lastData.b73 = state.lastData.pan = null;
    if (state.grids.b73) { try { state.grids.b73.destroy(); } catch (e) {} state.grids.b73 = null; }
    if (state.grids.pan) { try { state.grids.pan.destroy(); } catch (e) {} state.grids.pan = null; }
    /* remove any leaked tooltips from the DOM zoomed views + stray canvas tips */
    Array.prototype.forEach.call(
      document.querySelectorAll('body > .tooltip, body > .tooltipdomain, .pe-canvas-tip'),
      function (n) { n.remove(); }
    );
  }

  /* ---------------- render ------------------------------------------ */
  function render(container, opts) {
    opts = opts || {};
    return (async function () {
      teardown();

      /* Fail loudly & clearly if the ported PanEffect scripts aren't loaded.
         typeof is ReferenceError-safe for undeclared identifiers, unlike a
         bare read/write of main_id. */
      if (typeof runPanEffect !== 'function' ||
          typeof renderHeatmap !== 'function' ||
          typeof renderHeatmapPan !== 'function' ||
          typeof colorScale !== 'function' ||
          typeof numberToAminoAcid !== 'function') {
        var err = new Error(
          'PanEffect scripts not loaded. Expected ./js/pe/{support,genome,pan,dom,main}.js ' +
          '(main.js must be the edited version that defines runPanEffect). ' +
          'Check the <script src> paths in index.html and the Network tab for 404s.');
        err.code = 'PE_SCRIPTS_MISSING';
        console.error('[PanEffectEngine] ' + err.message);
        throw err;
      }

      installFetch();
      state.base = normalizeBase(opts.dataBase || state.base);
      state.variant = opts.variant || null;

      /* drive PanEffect's shared globals */
      main_id = opts.gene || '';
      main_option = opts.option || 'both';
      currentESM = ESM_OK.indexOf(opts.esm) >= 0 ? opts.esm : 'ESM2';
      parseId(opts.gene || '');

      /* swap the two full heatmaps to canvas (perf) */
      renderHeatmap = canvasRenderGene;
      renderHeatmapPan = canvasRenderPan;
      /* overlay a ring on the matching cell in the zoomed (div-based) view */
      installZoomHighlight();
      /* re-skin the d3 zoomed views: genome.js / pan.js read colorScale and
         colorScalePan at draw time, so reassigning them here is enough */
      THEME.restore();
      if (!THEME.install()) {
        console.warn('[PanEffectEngine] could not re-point colorScale/colorScalePan; ' +
          'the zoomed views will keep the classic colours.');
      }
      /* clean, minimal summary + skip the GWAS traits fetch (that lives in SNPFunction) */
      populateSummary = cleanSummary;
      loadAndDisplayTraits = function () {};

      /* inject the view skeleton */
      if (container.classList) container.classList.add('pe-native-root');
      container.innerHTML = skeleton();
      applyVisibility(main_option);

      /* external entries (SNPVersity / SNPFold) open in the MaizeGDB 2026 view */
      applyWgs(opts.wgs);

      /* run the (now callable) PanEffect pipeline */
      var ok = false;
      try { ok = await runPanEffect(); }
      catch (e) { console.error('[PanEffectEngine] pipeline error:', e); }

      wireScheme();
      populateHeteroticLegend();

      if (state.variant && state.variant.pos) {
        /* let the sliders finish wiring, then centre on the variant.
           fetchDataAndSetup() kicks off a pfam/DSSP fetch it does NOT await
           (deliberately, so the initial render isn't blocked on it), and that
           callback re-draws the zoomed number-line ticks for whatever start
           position was current when IT happens to resolve. Since it's not
           awaited, it can settle before or after this first centering pass —
           if it loses that race it can silently reset the ticks back to the
           default range even though the heatmap cells / WT line / box are
           already showing the right spot. Re-apply once more after giving
           that fetch time to finish, so the centered position is always the
           final word regardless of network timing. */
        setTimeout(function () { centerOnVariant(state.variant.pos); }, 0);
        setTimeout(function () { centerOnVariant(state.variant.pos); }, 400);
      }
      return ok;
    })();
  }

  /* ---------------- skeleton ---------------------------------------- */
  function skeleton() {
    return '' +
'<div id="loading-icon" style="display:none">Loading…</div>' +
'<div id="reference_gm_top" class="pe-gm-head"></div>' +

/* summary content is rendered by the clean populateSummary override */
'<div id="summary" class="content"></div>' +

'<div id="b73" class="content">' +
  '<span class="gene" id="reference_gm"></span>' +
  '<span id="wgs_span"><br><br>' +
    '<input type="radio" id="allVariants" name="variantEffect" value="all" checked>' +
    '<label for="allVariants">Show all variant effects</label>' +
    /* MaizeGDB 2024 option disabled for now
    '<input type="radio" id="maizeWGS" name="variantEffect" value="maize2024">' +
    '<label for="maizeWGS">MaizeGDB 2024 High Coverage variant effects</label>' +
    */
    '<input type="radio" id="maizeWGS2026" name="variantEffect" value="maize2026">' +
    '<label for="maizeWGS2026">MaizeGDB 2026 High Coverage variant effects</label><br>' +
  '</span>' +
  schemeRow('') +
  '<br>' +
  '<div class="sectionHeader">PFAM Domains</div>' +
  '<div id="pfam-wrap"><div id="pfamNumberLine" class="numberLine"></div><div id="pfamGeneModel" class="geneModel"></div></div>' +
  '<div class="sectionHeader">Secondary Structure</div>' +
  '<div id="dssp"><canvas id="proteinStructure" width="1400" height="80"></canvas></div>' +
  '<div class="sectionHeader">Variant Effects of full protein (heatmap)</div>' +
  '<div id="heatNumberLine" class="numberLine"></div>' +
  '<div class="heatmap-container"><div id="full-heatmap"><div id="highlight-box"></div></div></div>' +
  '<div class="slider-container" id="slider-container"><span id="slider"></span><span id="slider-value">1</span></div>' +
  '<div class="sectionHeader">Variant Effects of zoomed in region (heatmap)</div>' +
  '<div id="zoomNumberLine" class="numberLine"></div>' +
  '<div id="heatmap-container-zoom" class="heatmap-container-zoom"><div id="zoomed-heatmap"></div></div>' +
  '<div id="zoomWTLine" class="numberLine"></div>' +
'</div>' +

'<div id="pan-genome" class="content">' +
  '<span class="gene" id="pan_gm"></span>' +
  schemeRow('-pan') +
  '<div class="sectionHeader">PFAM Domains</div>' +
  '<div id="pfam-wrap-pan"><div id="pfamNumberLine-pan" class="numberLine"></div><div id="pfamGeneModel-pan" class="geneModel"></div></div>' +
  '<div class="sectionHeader">Secondary Structure</div>' +
  '<div id="dssp-pan"><canvas id="proteinStructure-pan" width="1400" height="80"></canvas></div>' +
  '<div class="sectionHeader" id="heatheader-pan">Heatmap of full protein</div>' +
  '<div id="heatNumberLine-pan" class="numberLine"></div>' +
  '<div id="colorBox-pan"></div>' +
  '<div class="heatmap-container-pan"><div id="full-heatmap-pan"><div id="highlight-box-pan"></div></div></div>' +
  '<div class="slider-container" id="slider-container-pan"><span id="slider-pan"></span><span id="slider-value-pan">1</span></div>' +
  '<div class="sectionHeader">Heatmap of zoomed in region</div>' +
  '<div id="zoomNumberLine-pan" class="numberLine"></div>' +
  '<div id="heatmap-container-zoom-pan" class="heatmap-container-zoom-pan"><div id="zoomed-heatmap-pan"></div><div id="heterotic-legend-pan" class="pe-heterotic-legend-wrap"></div></div>' +
  '<div id="zoomWTLine-pan" class="numberLine"></div>' +
'</div>';
  }

  /* ---------------- diagnostics --------------------------------------
     Run PanEffectEngine.diagnose() in the console if a colour change does
     not take. Each line is a precondition for the toggle working. */
  function diagnose() {
    var r = {
      themeLoaded:        !!window.PanEffectTheme,
      scheme:             THEME.getScheme(),
      originalsCaptured:  THEME.ready(),
      colorScaleRepointed:    !!(window.colorScale && window.colorScale.__peTheme),
      colorScalePanRepointed: !!(window.colorScalePan && window.colorScalePan.__peTheme),
      zoomFnGene:  typeof window.updateHeatmapZoom    === 'function',
      zoomFnPan:   typeof window.updateHeatmapZoomPan === 'function',
      radiosFound: document.querySelectorAll('.pe-scheme input[type=radio]').length,
      cachedRows:  { gene: state.lastData.b73 ? state.lastData.b73.length : 0,
                     pan:  state.lastData.pan ? state.lastData.pan.length : 0 },
      sampleGene:  { '0': THEME.colorFor('gene', 0),  '-11': THEME.colorFor('gene', -11) },
      samplePan:   { '0': THEME.colorFor('pan', 0),   '-11': THEME.colorFor('pan', -11) },
    };
    console.table(r);
    if (!r.colorScaleRepointed) console.warn('colorScale was not re-pointed — call PanEffectEngine.render() first.');
    if (!r.radiosFound) console.warn('No scheme radios in the DOM — the skeleton has not been injected yet.');
    if (!r.cachedRows.gene && !r.cachedRows.pan) console.warn('No cached rows — the data pipeline has not run.');
    return r;
  }

  window.PanEffectEngine = {
    render: render,
    teardown: teardown,
    refreshColors: refreshColors,
    setScheme: function (v) { THEME.setScheme(v); syncSchemeInputs(THEME.getScheme()); refreshColors(); },
    diagnose: diagnose,
  };
})();
