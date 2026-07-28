/* =====================================================================
 *  snphandoff.js — shared accession-handoff plumbing.
 *
 *  One place that owns:
 *    · the session-wide "add vs replace" mode  (S.handoff.mode)
 *    · the markup + wording for the "Replace current selection" control
 *    · the delta wording ("+18 added · 6 already selected · now 52")
 *    · dispatch of a handoff payload into SNPVersity
 *
 *  Any tool that hands accessions to SNPVersity (SNPFunction, SNPFold,
 *  PanEffect, …) should drop a mount point into its markup and dispatch
 *  through Handoff.toVersity() rather than rolling its own copy:
 *
 *    <span data-ho-mount data-ho-id="fnMergeReplace"
 *          data-ho-target="SNPVersity" data-ho-dataset="nam2026"></span>
 *
 *  …then call Handoff.sync() after the container's innerHTML is written.
 *
 *  Load order: after core.js (needs the global S), before the tools.
 * ===================================================================== */
(function () {
  'use strict';

  /* ---------- session state ---------------------------------------- */
  /* S.handoff.mode:  null = auto (decided by the heuristic below)
                      'add' | 'replace' = pinned by the user            */
  function st(){
    if (typeof S === 'undefined' || !S) return null;
    if (!S.handoff) S.handoff = { mode: null };
    return S.handoff;
  }

  function defaultsFor(ds){
    try {
      if (typeof Data !== 'undefined' && typeof Data.defaultSelectionFor === 'function')
        return Data.defaultSelectionFor(ds) || [];
    } catch (e) { /* no data layer yet */ }
    return [];
  }

  /* A selection is "pristine" while it is still exactly the dataset's
     preselected founders — nothing of the user's own is in it yet, so
     replacing it destroys no work. Once they have curated it at all,
     the safe default flips to add. */
  function isPristine(){
    if (typeof S === 'undefined' || !S || !S.selected) return true;
    const def = defaultsFor(S.dataset);
    if (S.selected.size !== def.length) return false;
    for (let i = 0; i < def.length; i++) if (!S.selected.has(def[i])) return false;
    return true;
  }

  function mode(){
    const s = st();
    if (s && (s.mode === 'add' || s.mode === 'replace')) return s.mode;
    return 'add';   // default: off (merge into the current selection) until the user pins it
  }
  function setMode(m){
    const s = st(); if (!s) return;
    s.mode = (m === 'add') ? 'add' : 'replace';   // explicit choice pins it for the session
  }
  function isPinned(){ const s = st(); return !!(s && s.mode); }

  function selectedCount(){
    return (typeof S !== 'undefined' && S && S.selected) ? S.selected.size : 0;
  }

  /* Accession IDs are not portable between dataset families, so a handoff
     that switches dataset can only ever reset the selection. */
  function sameDataset(ds){
    if (ds == null || ds === '') return true;
    if (typeof S === 'undefined' || !S) return true;
    return String(ds) === String(S.dataset);
  }

  function esc(s){
    return String(s == null ? '' : s)
      .replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
  }

  /* ---------- the shared control ------------------------------------ */
  /* Same label everywhere it appears (SNPFunction, SNPFold, the upload
     panel, the arrival banner) — that repetition is what teaches it. */
  const LABEL = 'Replace current selection';

  function controlHTML(o){
    o = o || {};
    const id     = o.id || 'hoReplace';
    const target = o.target || '';
    const cross  = !sameDataset(o.dataset);
    const rep    = cross || mode() === 'replace';
    const n      = selectedCount();
    const where  = target ? ` in ${esc(target)}` : '';
    const noun   = `accession${n === 1 ? '' : 's'}`;

    const hint = cross
      ? `Different dataset — the selection${where} will be reset either way.`
      : rep ? `Replaces the ${n} ${noun} currently selected${where}.`
            : `Keeps the ${n} ${noun} currently selected${where} and adds to them.`;

    return `<label class="ho-opt${cross ? ' off' : ''}"${cross ? ' title="Accessions are not shared between datasets."' : ''}>
      <input type="checkbox" id="${esc(id)}" ${rep ? 'checked' : ''} ${cross ? 'disabled' : ''}
             onchange="Handoff.setMode(this.checked?'replace':'add');Handoff.sync()">
      <span class="ho-lbl">${LABEL}</span>
      <span class="ho-hint">${hint}</span>
    </label>`;
  }

  /* Fill every mount point on the page (or inside `root`). Safe to call
     as often as you like — it is a full re-render of the controls only. */
  function sync(root){
    injectCSS();
    const scope = root || document;
    if (!scope || !scope.querySelectorAll) return;
    scope.querySelectorAll('[data-ho-mount]').forEach(function (el) {
      el.innerHTML = controlHTML({
        id:      el.getAttribute('data-ho-id'),
        target:  el.getAttribute('data-ho-target'),
        dataset: el.getAttribute('data-ho-dataset')
      });
    });
  }

  /* ---------- shared delta wording ---------------------------------- */
  /* d = {added, already, removed, missing, total} */
  function deltaText(d){
    d = d || {};
    const bits = [];
    if (d.added)   bits.push(`<b>+${d.added}</b> added`);
    if (d.already) bits.push(`${d.already} already selected`);
    if (d.removed) bits.push(`${d.removed} replaced`);
    if (d.missing) bits.push(`${d.missing} not in this dataset`);
    if (!bits.length) bits.push('nothing changed');
    const t = d.total || 0;
    return bits.join(' · ') + ` · selection is now <b>${t}</b> accession${t === 1 ? '' : 's'}`;
  }

  /* ---------- dispatch ---------------------------------------------- */
  /* Adds the merge flag if the caller did not set one, then hands over.
     Falls back to the pre-hook contract (S.pendingVersity + go) so an
     older snpversity.js still receives the payload. */
  function toVersity(payload){
    if (!payload) return;
    const p = {};
    for (const k in payload) if (Object.prototype.hasOwnProperty.call(payload, k)) p[k] = payload[k];
    if (p.merge !== 'add' && p.merge !== 'replace') p.merge = mode();

    if (typeof window.versityRequest === 'function'){ window.versityRequest(p); return; }
    if (typeof S !== 'undefined' && S){
      S.pendingVersity = p;
      if (typeof go === 'function') go('snpversity');
    }
  }

  /* ---------- styles ------------------------------------------------ */
  function injectCSS(){
    if (document.getElementById('snphandoff-css')) return;
    const s = document.createElement('style'); s.id = 'snphandoff-css';
    s.textContent = `
      .ho-opt{display:inline-flex;align-items:baseline;gap:7px;flex-wrap:wrap;cursor:pointer;
        font-size:12px;color:var(--ink,#141922);line-height:1.5}
      .ho-opt input{margin:0;position:relative;top:1px;cursor:pointer;flex:0 0 auto}
      .ho-lbl{font-weight:600}
      .ho-hint{font-size:11.5px;color:var(--muted,#6b7789);font-weight:400}
      .ho-opt.off{cursor:default;opacity:.7}
      .ho-opt.off input{cursor:default}
      .ho-undo{margin-left:2px}
    `;
    document.head.appendChild(s);
  }

  /* ---------- public ------------------------------------------------ */
  window.Handoff = {
    LABEL: LABEL,
    mode: mode,
    setMode: setMode,
    isPinned: isPinned,
    isPristine: isPristine,
    sameDataset: sameDataset,
    selectedCount: selectedCount,
    controlHTML: controlHTML,
    sync: sync,
    deltaText: deltaText,
    toVersity: toVersity
  };
})();
