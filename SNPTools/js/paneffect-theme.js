/* =====================================================================
 *  paneffect-theme.js — switchable colour schemes for PanEffect.
 *
 *  Two schemes:
 *    'classic'  — PanEffect's original stepped blue->red ramp. Delegates
 *                 to the real colorScale / colorScalePan from support.js,
 *                 so it is pixel-identical to the original build.
 *    'snptools' — SNPImpact's continuous red->amber->green ramp, lifted
 *                 verbatim from snpimpact.js scoreColor():
 *
 *                   t = (v - lo) / (hi - lo)              clamped [0,1]
 *                   r = t < .5 ? 255 : 255 * (1 - (t-.5)*2)
 *                   g = t < .5 ? 255 * t * 2 : 200        floored at 60
 *                   b = 60
 *
 *  ---- Sentinels (see support.js) -------------------------------------
 *  Score is not purely continuous. support.js reserves two values, and
 *  BOTH schemes must honour them or the views lose information:
 *
 *              undefined / NaN        score === 0
 *    gene      white (no data)        black  -> the wild-type diagonal
 *    pan       #DDDDDD (no data)      #00429d -> exact match to B73
 *
 *  A naive ramp paints score 0 at the benign end, which erases the WT
 *  diagonal in the gene view and every exact-match cell in the pan view.
 *  The SNPTools scheme therefore keeps distinct sentinel colours.
 *
 *  ---- How the toggle reaches the zoomed views ------------------------
 *  updateHeatmapZoom (genome.js) and updateHeatmapZoomPan (pan.js) build
 *  plain <div> cells and call the GLOBAL colorScale / colorScalePan by
 *  name at draw time. So re-pointing those globals is enough — the views
 *  just need to be redrawn, which the engine does by dispatching 'input'
 *  on the zoom sliders (the same path centerOnVariant already uses).
 *
 *  Load order (classic script, shares PanEffect's global scope):
 *    d3 -> support.js -> genome.js -> pan.js -> dom.js -> main.js
 *       -> paneffect-heatmap.js -> paneffect-theme.js
 *       -> paneffect-engine.js  -> snppaneffect.js
 *
 *  paneffect-theme.js MUST load after support.js (to capture the
 *  originals) and before paneffect-engine.js.
 *
 *  Exposes window.PanEffectTheme.
 * ===================================================================== */
(function () {
  'use strict';
  if (typeof window === 'undefined') return;

  var STORE_KEY = 'snptools.pe.scheme';

  /* ---------------- score domains -------------------------------------
     Confirmed against support.js: the stepped thresholds run from >0 down
     to <=-22 in steps of 2, so the scored range is exactly -22..0. */
  var DOMAINS = {
    impact: { lo: -12, hi: 6 },   // SNPImpact / SNPVersity DNA+protein LM
    gene:   { lo: -22, hi: 0 },   // PanEffect B73 view
    pan:    { lo: -22, hi: 0 },   // PanEffect pan view
  };

  /* ---------------- sentinel + neutral colours ------------------------
     GAP matches the "#e6e6e6" hardcoded in pan.js so the canvas and
     zoomed pan views agree. */
  var NEUTRAL = {
    gap: '#e6e6e6',
    classic: {
      gene: { blank: 'white',   zero: 'black'   },
      pan:  { blank: '#DDDDDD', zero: '#00429d' },
    },
    snptools: {
      /* WT diagonal stays the darkest ink in the palette, so it reads as
         structure rather than as a score; exact-match keeps SNPTools blue */
      gene: { blank: '#ffffff', zero: '#13264a' },
      pan:  { blank: '#DDDDDD', zero: '#2563eb' },
    },
  };

  var CLASSIC_COLORS = ['#00429d','#3860aa','#587fb3','#78a0b7','#9ac0b3','#c1e19e',
    '#ffff00','#ffd337','#fea447','#f1784d','#db4c4d','#bd2147','#93003a'];

  /* the genuine support.js scales, captured before we re-point them */
  var ORIGINAL = { gene: null, pan: null };

  var scheme = 'snptools';
  var listeners = [];

  /* ---------------- ramp ------------------------------------------------ */
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
    var n = CLASSIC_COLORS.length;
    return CLASSIC_COLORS[Math.round(clamp01(1 - t) * (n - 1))];
  }
  function tFor(kind, score) {
    var d = DOMAINS[kind] || DOMAINS.gene;
    return clamp01((+score - d.lo) / (d.hi - d.lo));
  }

  /* ---------------- public colour lookup --------------------------------
     Mirrors support.js's contract exactly, including both sentinels. */
  function colorFor(kind, score) {
    kind = (kind === 'pan') ? 'pan' : (kind === 'impact' ? 'impact' : 'gene');
    var neutralKind = (kind === 'impact') ? 'gene' : kind;

    /* Classic: hand the raw value straight to support.js so every input —
       including degenerate ones like NaN — behaves exactly as it always
       has. Checking sentinels ourselves first would diverge on the edges. */
    if (scheme === 'classic') {
      var orig = ORIGINAL[neutralKind];
      if (typeof orig === 'function') {
        try { return orig(score); } catch (e) { /* fall through to the copy */ }
      }
      var nc = NEUTRAL.classic[neutralKind];
      if (score === undefined || score === null || score === '' || isNaN(+score)) return nc.blank;
      if (+score === 0) return nc.zero;
      return classicRamp(tFor(kind, +score));
    }

    var n = NEUTRAL.snptools[neutralKind];
    if (score === undefined || score === null || score === '' || isNaN(+score)) return n.blank;
    var v = +score;
    if (v === 0) return n.zero;
    return snpRamp(tFor(kind, v));
  }
  function percentFor(kind, score) { return tFor(kind, score) * 100; }
  function gapColor() { return NEUTRAL.gap; }

  /* ---------------- legend helpers -------------------------------------- */
  function stops(n, reverse) {
    var out = [], i;
    if (scheme === 'classic') {
      out = CLASSIC_COLORS.slice().reverse();      // deleterious -> benign
    } else {
      n = n || 9;
      for (i = 0; i < n; i++) out.push(snpHex(i / (n - 1)));
    }
    return reverse ? out.reverse() : out;
  }
  function gradientCSS(deg, reverse) {
    var s = stops(9, reverse), i;
    deg = (deg == null) ? 90 : deg;
    if (scheme === 'classic') {
      /* hard stops so the discrete bands stay visible */
      var parts = [], w = 100 / s.length;
      for (i = 0; i < s.length; i++) {
        parts.push(s[i] + ' ' + (i * w).toFixed(2) + '% ' + ((i + 1) * w).toFixed(2) + '%');
      }
      return 'linear-gradient(' + deg + 'deg,' + parts.join(',') + ')';
    }
    return 'linear-gradient(' + deg + 'deg,' + s.join(',') + ')';
  }

  /* ---------------- scheme switching ------------------------------------ */
  function getScheme() { return scheme; }
  function setScheme(next, opts) {
    next = (next === 'classic') ? 'classic' : 'snptools';
    if (next === scheme) return scheme;
    scheme = next;
    try { window.localStorage.setItem(STORE_KEY, scheme); } catch (e) {}
    install();
    if (!(opts && opts.silent)) {
      listeners.forEach(function (fn) { try { fn(scheme); } catch (e) {} });
    }
    return scheme;
  }
  function onChange(fn) { if (typeof fn === 'function') listeners.push(fn); }
  function restore() {
    try {
      var v = window.localStorage.getItem(STORE_KEY);
      if (v === 'classic' || v === 'snptools') scheme = v;
    } catch (e) {}
    return scheme;
  }

  /* ---------------- install over PanEffect's globals ---------------------
     support.js declares colorScale / colorScalePan as top-level function
     declarations, so they are writable properties of the global object.
     genome.js and pan.js call them by name at draw time, which means
     re-pointing here re-skins the zoomed views too. */
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

    window.colorScale    = gene;
    window.colorScalePan = pan;
    try { /* eslint-disable no-undef */
      colorScale = gene; colorScalePan = pan;
    } catch (e) { /* the window write above is the one that counts */ }

    return window.colorScale === gene;
  }

  /* true once we hold the real support.js scales — useful for a console
     sanity check: PanEffectTheme.ready() should be true after render() */
  function ready() {
    return !!(ORIGINAL.gene && ORIGINAL.pan && window.colorScale && window.colorScale.__peTheme);
  }

  window.PanEffectTheme = {
    DOMAINS: DOMAINS,
    CLASSIC_COLORS: CLASSIC_COLORS,
    colorFor: colorFor,
    gapColor: gapColor,
    percentFor: percentFor,
    stops: stops,
    gradientCSS: gradientCSS,
    getScheme: getScheme,
    setScheme: setScheme,
    onChange: onChange,
    restore: restore,
    install: install,
    ready: ready,
  };
})();
