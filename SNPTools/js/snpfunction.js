/* =====================================================================
 *  snpfunction.js — Gene function & allele mining.
 *
 *  Gene-scoped view (independent of a variant region): a functional
 *  dossier (Pfam domains, size, links), the gene's variant burden across
 *  the WHOLE panel, and a damaging / knockout allele catalog listing which
 *  accessions carry each damaging allele. Pulls Data.geneFunction(gene).
 * ===================================================================== */
(function () {

  const FN = {
    gene:'Zm00001eb406050',        // default example gene model
    dataset:null, data:null, loading:false, openId:null,
    root:null,                     // persistent DOM container — survives tool switches
    loaded:false,                  // a gene's content is (being) rendered into root
    loadedGene:null,               // which gene that content is for
    annotation:undefined,          // functional-annotation record: undefined=loading, null=none, object=loaded
    goCurated:false,               // GO filter: false = all terms, true = curated only (drops predicted-only)
  };

  /* directory of per-gene functional-annotation JSON (one file per canonical model) */
  const ANN_DIR = './data/function/annotations/';

  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  function scoreColor(v){ const lo=-12,hi=6,t=Math.max(0,Math.min(1,(v-lo)/(hi-lo))); const r=t<.5?255:Math.round(255*(1-(t-.5)*2)); const g=t<.5?Math.round(255*t*2):200; return `rgb(${r},${Math.max(60,g)},60)`; }

  /* Choose black or white text from the actual score-chip background.
     Pale yellow and green backgrounds receive dark text automatically,
     while darker red/green backgrounds retain white text. */
  function contrastTextColor(background){
    const text = String(background || '').trim();
    let rgb = null;
    const hex = text.match(/^#([0-9a-f]{3}|[0-9a-f]{6})$/i);
    if (hex){
      let h = hex[1];
      if (h.length === 3) h = h.split('').map(c => c+c).join('');
      rgb = [parseInt(h.slice(0,2),16), parseInt(h.slice(2,4),16), parseInt(h.slice(4,6),16)];
    } else {
      const nums = text.match(/[0-9.]+/g);
      if (nums && nums.length >= 3) rgb = nums.slice(0,3).map(Number);
    }
    if (!rgb || rgb.some(v => !Number.isFinite(v))) return '#ffffff';
    const linear = rgb.map(v => {
      const c = Math.max(0,Math.min(255,v))/255;
      return c <= 0.04045 ? c/12.92 : Math.pow((c+0.055)/1.055,2.4);
    });
    const luminance = 0.2126*linear[0] + 0.7152*linear[1] + 0.0722*linear[2];
    return luminance > 0.179 ? '#111827' : '#ffffff';
  }
  function scoreCell(v){
    if (v==null) return '<span style="color:var(--faint)">—</span>';
    const bg = scoreColor(v);
    const fg = contrastTextColor(bg);
    return `<span class="imp-score" style="background:${bg};color:${fg} !important;text-shadow:none">${v>0?'+':''}${v.toFixed(1)}</span>`;
  }
  function prioPill(p){ return `<span class="prio ${(p||'LOW').toLowerCase()}">${p||'LOW'}</span>`; }
  function consPill(v){ return `<span class="cons ${v.consClass}">${v.consequence}</span>`; }

  /* ---------- per-allele jumps: PanEffect (effects ↗) and SNPFold (fold ↗) ----------
     Both act on a row in the damaging-allele catalog, so they resolve the gene from
     FN.data rather than taking it as an argument. */
  function isMissense(v){
    return !!v && (v.consClass==='missense' || /missense/i.test(v.consequence||''));
  }
  /* Consequence classes SNPFold has a variant record for. Synonymous alleles are
     deliberately absent — SNPFold's coding-variant table does not carry them — as is
     'other' (non-coding). Trim or extend this one list to change which rows link. */
  const FOLD_CONS = ['missense', 'lof', 'indel', 'splice'];
  function canFold(v){ return !!v && FOLD_CONS.indexOf(v.consClass) !== -1; }
  function fnGene(){ return (FN.data && FN.data.gene) || FN.gene || ''; }

  /* Protein substitution ("A123T") when the record carries one, in whichever shape.
     Returns '' when the allele is not a single amino-acid change. */
  function aaSub(v){
    if (!v) return '';
    if (nonEmptyStr(v.sub)) return String(v.sub);
    if (nonEmptyStr(v.aaRef) && nonEmptyStr(v.aaAlt) && v.resi!=null) return `${v.aaRef}${v.resi}${v.aaAlt}`;
    /* v.variant is often already the protein notation for coding alleles */
    const m = String(v.variant||'').match(/^([A-Z*])(\d+)([A-Z*])$/i);
    return m ? m[0] : '';
  }
  function nonEmptyStr(x){ return x != null && String(x) !== ''; }

  function peJump(v){
    if (!isMissense(v) || !fnGene()) return '';
    const sub = aaSub(v);
    return ` <a class="pe-jump" href="#" title="View ${esc(sub||'this substitution')} in PanEffect"
      onclick="event.stopPropagation();FUNCTION.panEffect('${esc(v.id)}');return false;">effects ↗</a>`;
  }
  function foldJump(v){
    if (!canFold(v) || !fnGene()) return '';
    return ` <a class="fold-jump" href="#" title="Show ${esc(aaSub(v)||v.variant||'this allele')} on the predicted protein structure in SNPFold"
      onclick="event.stopPropagation();FUNCTION.fold('${esc(v.id)}');return false;">fold ↗</a>`;
  }

  /* view switches — same contracts the other tools use */
  function panEffectTo(v){
    const gene = fnGene(); if (!gene) return;
    if (typeof goPanEffect !== 'function') return go('paneffect');
    const sub = isMissense(v) ? aaSub(v) : '';
    goPanEffect(gene, sub ? {variant:sub} : undefined);
  }
  /* Hands the gene *and* the site to SNPFold. SNPFold reuses an already-loaded gene,
     so clicking several of these in a row triggers no further HDF5 queries. */
  function foldTo(v){
    const gene = fnGene(); if (!gene) return;
    if (typeof goFold !== 'function') return go('snpfold');
    if (!canFold(v)) return goFold(gene);          // nothing to highlight — open the gene
    const d = FN.data || {};
    const site = { sub: aaSub(v), effect: v.consequence };
    if (d.chr != null)   site.chr  = d.chr;        // descriptor fields are optional;
    if (v.pos != null)   site.pos  = v.pos;        // include only what this record has,
    if (v.ref != null)   site.ref  = v.ref;        // so the matcher falls back cleanly
    if (v.alt != null)   site.alt  = v.alt;
    if (v.resi != null)  site.resi = v.resi;
    goFold(gene, site);
  }

  /* ---------- dataset chooser (compact; shares Data.datasets() with SNPVersity) ---------- */
  function nonEmpty(v){ return v != null && v !== ''; }
  function datasetText(dataset){
    if (dataset == null) return '';
    if (typeof dataset === 'string' || typeof dataset === 'number') return String(dataset);
    if (typeof dataset === 'object'){
      return [dataset.id, dataset.key, dataset.value, dataset.name, dataset.label, dataset.title].filter(nonEmpty).join(' ');
    }
    return String(dataset);
  }
  function fnDatasets(){
    try {
      if (typeof Data !== 'undefined' && typeof Data.datasets === 'function'){
        const ds = Data.datasets();
        if (Array.isArray(ds)) return ds;
      }
    } catch (e) { /* no dataset catalog available */ }
    return [];
  }
  function datasetId(d){
    if (d && typeof d === 'object'){
      if (nonEmpty(d.id))  return d.id;
      if (nonEmpty(d.key)) return d.key;
      return datasetText(d);
    }
    return d;
  }
  function isCurrentDataset(d){
    const a = String(datasetId(d)), b = String(datasetId(FN.dataset));
    if (a && a === b) return true;
    const na = datasetText(d).toLowerCase().replace(/[^a-z0-9]+/g, '');
    const nb = datasetText(FN.dataset).toLowerCase().replace(/[^a-z0-9]+/g, '');
    return !!na && na === nb;
  }
  function datasetLabel(d){
    const name = (d && nonEmpty(d.name)) ? String(d.name) : datasetText(d);
    let sub = (d && nonEmpty(d.sub)) ? String(d.sub) : '';
    if (sub && name.toLowerCase().includes(sub.toLowerCase())) sub = '';
    return { name: name || 'Dataset', sub };
  }
  function datasetCardsHTML(){
    const list = fnDatasets();
    if (!list.length) return '';
    return list.map(d=>{
      const id  = datasetId(d);
      const sel = isCurrentDataset(d);
      const { name, sub } = datasetLabel(d);
      return `<button type="button" class="fn-ds ${sel?'sel':''}" onclick="FUNCTION.pickDataset('${esc(id)}')">
        <span class="fn-ds-dot"></span>
        <span class="fn-ds-txt"><span class="fn-ds-name">${esc(name)}</span>${sub?`<span class="fn-ds-sub">${esc(sub)}</span>`:''}</span>
      </button>`;
    }).join('');
  }
  function datasetChooser(){
    const cards = datasetCardsHTML();
    if (!cards) return '';   // if the data layer exposes no catalog, show nothing
    return `<div class="card pad fn-ds-card" style="margin-bottom:16px">
      <div class="fn-ds-head">Dataset</div>
      <div class="fn-ds-grid" id="fnDsGrid">${cards}</div>
    </div>`;
  }
  function syncDatasetChooser(){
    const grid = FN.root && FN.root.querySelector('#fnDsGrid');
    if (grid) grid.innerHTML = datasetCardsHTML();
  }

  /* ---------- persistent container ----------
     Content is kept in FN.root and re-attached whenever the tool is shown again, so
     switching to another tool and back does not leave the page empty or force a
     re-analysis. We hold a JS reference, so #page being overwritten only detaches it. */
  function ensureRoot(page){
    if (!FN.root){
      FN.root = document.createElement('div');
      FN.root.className = 'snpfunction-root';
    }
    if (FN.root.parentNode !== page){
      page.innerHTML = '';
      page.appendChild(FN.root);
    }
  }

  /* ---------- render ---------- */
  /* Entry point called by the suite shell on navigation.
     - Arriving from another tool (SNPVersity/SNPImpact set S.functionGene) → autoload.
     - Arriving from the side-panel menu → do NOT auto-analyze; show a landing page the
       first time and simply re-show the already-loaded page on later visits. */
  function render(page){
    page = page || document.getElementById('page');
    injectCSS();
    ensureRoot(page);

    /* reflect the app-wide dataset for the picker highlight (until the user picks one) */
    if (FN.dataset == null && typeof S !== 'undefined' && S && S.dataset != null) FN.dataset = S.dataset;

    const requestedGene = (typeof S !== 'undefined' && S && S.functionGene) ? S.functionGene : null;

    if (requestedGene){
      S.functionGene = null;
      if (typeof S !== 'undefined' && S && S.functionDataset != null) FN.dataset = S.functionDataset;
      if (requestedGene !== FN.loadedGene || !FN.loaded){
        FN.gene = requestedGene; FN.data = null; FN.openId = null;
        analyzeGene();
      } else {
        reshowLoaded();
      }
      return;
    }

    if (FN.loaded){ reshowLoaded(); return; }   // preserve rendered content, no reload
    renderLanding();                            // first menu visit: no auto-analysis
  }

  /* re-show already-rendered content without re-analyzing */
  function reshowLoaded(){
    syncDatasetChooser();
    if (typeof attachTT==='function') attachTT();
  }

  /* landing: dataset chooser + gene search (prefilled with the default example),
     but no analysis until the user presses the button */
  function renderLanding(){
    FN.root.innerHTML = datasetChooser() + searchBar() + emptyState();
    if (typeof attachTT==='function') attachTT();
  }

  /* the heavy path: analyze the gene across the panel and render the full dossier */
  function analyzeGene(){
    FN.loaded = true;                // commit to the loaded view (keep across navigation)
    FN.loadedGene = FN.gene;
    FN.loading = true;
    FN.goCurated = false;            // reset GO filter for the new gene
    loadAnnotation(FN.gene);         // fetch functional annotation in parallel (independent of variant data)
    FN.root.innerHTML = datasetChooser() + searchBar() + `<div class="loading" style="padding:44px;text-align:center"><div class="spinner"></div><div>Analyzing <b>${esc(FN.gene)}</b> across the panel…</div></div>`;
    Data.geneFunction(FN.gene, FN.dataset)
      .then(d => { FN.data = d; FN.loading = false; paint(); })
      .catch(e => { FN.loading = false; FN.data = {gene:FN.gene, error:(e&&e.message)||'failed'}; paint(); });
  }

  /* Load the functional-annotation record for a gene from ANN_DIR/<gene>.json.
     Independent of the variant analysis: a missing file (404) or parse error just means
     "no annotation" and never blocks the rest of the page. We only repaint once the main
     variant data is present, so the annotation slots in without flashing the landing view. */
  function loadAnnotation(gene){
    const g = gene;
    FN.annotation = undefined;                       // mark as loading (shows a placeholder)
    fetch(ANN_DIR + encodeURIComponent(gene) + '.json', {cache:'force-cache'})
      .then(r => r.ok ? r.json() : null)
      .then(a => { if (FN.gene !== g) return; FN.annotation = a || null; if (FN.data) paint(); })
      .catch(() => { if (FN.gene !== g) return; FN.annotation = null; if (FN.data) paint(); });
  }

  /* render the loaded content from cached FN.data (used after analysis + on toggle) */
  function paint(){
    const d = FN.data;
    if (!d){ renderLanding(); return; }
    if (d.error){
      // variant analysis failed, but functional annotation may still be worth showing
      FN.root.innerHTML = datasetChooser() + searchBar()
        + notice(`Couldn’t analyze “${esc(FN.gene)}” across the panel: ${esc(d.error)}`)
        + annotationSection();
      if (typeof attachTT==='function') attachTT();
      return;
    }
    FN.root.innerHTML = datasetChooser() + searchBar() + hero(d) + dossier(d) + annotationSection() + burden(d) + catalog(d);
    if (typeof attachTT==='function') attachTT();
  }

  function searchBar(){
    return `<div class="card pad" style="margin-bottom:16px;display:flex;gap:10px;align-items:center;flex-wrap:wrap">
      <span style="font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px">Gene model</span>
      <input id="fnGeneInput" value="${esc(FN.gene||'')}" placeholder="e.g. Zm00001eb406050" spellcheck="false"
        style="flex:1;min-width:280px;border:1px solid var(--line);border-radius:9px;padding:9px 11px;font-family:var(--mono);font-size:13px"
        onkeydown="if(event.key==='Enter')FUNCTION.load()">
      <button class="btn" onclick="FUNCTION.load()">Analyze gene</button>
    </div>`;
  }
  function emptyState(){
    return `<div class="empty-state"><div class="ei">${ICONS.leaf||ICONS.star||''}</div>
      <h3>Pick a gene to mine its functional variation</h3>
      <p>SNPFunction summarizes a gene across the whole panel: its Pfam domains, the burden of
      coding variation, and a catalog of damaging / knockout alleles with the accessions that carry them.
      Choose a dataset above and a gene model, then press <b>Analyze gene</b>.</p>
      <div style="margin-top:14px"><button class="btn primary" onclick="FUNCTION.load()">Analyze ${esc(FN.gene||'')}</button></div></div>`;
  }
  function notice(html){ return `<div style="text-align:center;padding:26px;color:var(--muted);max-width:640px;margin:0 auto">${html}</div>`; }

  /* ---------- hero + dossier ---------- */
  function hero(d){
    return `<div class="sec"><div class="bar"></div><div>
      <div class="n">GENE FUNCTION &amp; ALLELE MINING</div>
      <h2 style="margin-bottom:2px">${esc(d.gene)}</h2>
      <p>Functional variation for this gene across <b>${d.nAccessions.toLocaleString()}</b> ${esc(d.datasetName)} accessions —
      domains, variant burden, and which lines carry damaging or knocked-out copies.</p>
    </div></div>`;
  }
  function dossier(d){
    const region = `${d.chr}:${(+d.start).toLocaleString()}–${(+d.end).toLocaleString()}${d.strand?` (${d.strand})`:''}`;
    const jb = `https://jbrowse.maizegdb.org/index.html?data=B73&loc=${encodeURIComponent(d.gene)}`;
    const pe = `https://maizegdb.org/effect/maize_v2/index.html?id=${encodeURIComponent(d.gene)}`;
    const mg = `https://maizegdb.org/gene_center/gene/${encodeURIComponent(d.gene)}`;
    const pg = `https://pangenome-viewer.maizegdb.org/?set=NAM&geneID=${encodeURIComponent(d.gene)}`;
    const doms = (d.domains||[]).length
      ? d.domains.map(x=>domTag(`${x.name} (${x.pfam})`)).join(' ')
      : '<span class="muted">No Pfam domains annotated</span>';
    return `<div class="card pad" style="margin-bottom:16px">
      <div class="fn-grid">
        <div><div class="fn-k">Location</div><div class="fn-v mono">${region}</div></div>
        <div><div class="fn-k">Protein</div><div class="fn-v">${d.protLen?`${d.protLen} aa`:'—'}${d.protein?` <span class="muted mono">${esc(d.protein)}</span>`:''}</div></div>
        <div><div class="fn-k">Variants in gene</div><div class="fn-v"><b>${d.nVariants.toLocaleString()}</b></div></div>
        <div style="flex:1 1 100%"><div class="fn-k">Pfam domains</div><div class="fn-v">${doms}</div></div>
      </div>
      <div style="display:flex;gap:8px;flex-wrap:wrap;margin-top:12px">
        <button class="btn" onclick="goFold('${esc(d.gene)}')">${ICONS.fold} View structure (SNPFold)</button>
        <a class="btn ghost" href="${pe}" target="_blank" rel="noopener">PanEffect ↗</a>
        <a class="btn ghost" href="${jb}" target="_blank" rel="noopener">JBrowse ↗</a>
        <a class="btn ghost" href="${mg}" target="_blank" rel="noopener">MaizeGDB ↗</a>
        <a class="btn ghost" href="${pg}" target="_blank" rel="noopener">Pangenome viewer ↗</a>
      </div>
    </div>`;
  }

  /* ---------- burden ---------- */
  function burden(d){
    const b = d.burden, bc = b.byClass;
    const sec = Data.hasSecondaryScores(d.dataset);
    const total = d.nVariants || 1;
    const barSeg = (n,cls,label)=> n? `<span class="fn-seg ${cls}" style="width:${Math.max(2,100*n/total)}%" title="${label}: ${n}"></span>`:'';
    const nss = b.nonsynSyn==null ? '∞' : b.nonsynSyn;
    const af = b.afSpectrum;
    return `<div class="card pad" style="margin-bottom:16px">
      <div class="fn-h">Variant burden</div>
      <div class="fn-stats">
        ${stat('Non-syn : syn', nss, 'coding constraint (higher = more nonsynonymous)')}
        ${stat('Exon : intron', b.exonIntron==null?'∞':b.exonIntron, 'variants in exons vs introns of the gene model')}
        ${stat('Domain-disrupting', b.domainDisrupting, 'coding variants inside a Pfam domain')}
        ${stat('Candidate KO lines', d.koLines, 'accessions homozygous for a loss-of-function allele')}
        ${stat('Mean PlantCAD1', b.meanPlantcad==null?'—':b.meanPlantcad, 'average PlantCAD1 DNA language-model score')}
        ${sec?stat('Mean PlantCAD2', b.meanPlantcad2==null?'—':b.meanPlantcad2, 'average PlantCAD2 (2026) DNA language-model score'):''}
        ${stat('Mean ESM', b.meanEsm==null?'—':b.meanEsm, 'average ESM protein language-model score')}
        ${sec?stat('Mean ESM2', b.meanEsm2==null?'—':b.meanEsm2, 'average ESM2 (2026) protein language-model score'):''}
        ${sec?stat('Mean ESM3', b.meanEsm3==null?'—':b.meanEsm3, 'average ESM3 (2026) protein language-model score'):''}
      </div>
      <div class="fn-barwrap">
        <div class="fn-bar">
          ${barSeg(bc.lof,'lof','Loss-of-function')}${barSeg(bc.splice,'splice','Splice')}${barSeg(bc.missense,'missense','Missense')}${barSeg(bc.indel,'indel','In-frame indel')}${barSeg(bc.syn,'syn','Synonymous')}${barSeg(bc.other,'other','Non-coding / other')}
        </div>
        <div class="fn-legend">
          ${leg('lof','LOF',bc.lof)} ${leg('splice','Splice',bc.splice)} ${leg('missense','Missense',bc.missense)} ${leg('indel','Indel',bc.indel)} ${leg('syn','Synon.',bc.syn)} ${leg('other','Other',bc.other)}
        </div>
      </div>
      <div class="fn-af">Allele frequency: <b>${af.rare}</b> rare (&lt;1%) · <b>${af.low}</b> low (1–5%) · <b>${af.common}</b> common (≥5%)</div>
    </div>`;
  }
  function stat(k,v,tip){ return `<div class="fn-stat" title="${esc(tip)}"><div class="fn-statv">${v}</div><div class="fn-statk">${k}</div></div>`; }
  function leg(cls,label,n){ return `<span class="fn-lg"><span class="fn-sw ${cls}"></span>${label} ${n}</span>`; }

  /* ---------- damaging / knockout catalog ---------- */
  function catalog(d){
    if (!d.damaging.length)
      return `<div class="card pad"><div class="fn-h">Damaging &amp; knockout alleles</div><div class="muted" style="padding:6px 0">No loss-of-function or high-impact damaging alleles found in this gene across the panel.</div></div>`;
    const sec = Data.hasSecondaryScores(d.dataset);
    const rows = d.damaging.map(v=>{
      const open = FN.openId===v.id;
      return `<tr class="imp-row ${open?'open':''}" onclick="FUNCTION.toggle('${v.id}')">
        <td class="c-mono c-alt" style="padding-left:11px">${v.variant}</td>
        <td>${consPill(v)}${peJump(v)}${foldJump(v)}</td>
        <td>${domTag(v.domain)}</td>
        <td class="num">${scoreCell(v.plantcad)}</td>
        ${sec?`<td class="num">${scoreCell(v.plantcad2)}</td>`:''}
        <td class="num">${scoreCell(v.esm)}</td>
        ${sec?`<td class="num">${scoreCell(v.esm2)}</td><td class="num">${scoreCell(v.esm3)}</td>`:''}
        <td>${prioPill(v.priority)}</td>
        <td class="num">${v.het}</td>
        <td class="num">${v.hom?`<b>${v.hom}</b>`:'0'}</td>
        <td class="num">${(v.af*100).toFixed(1)}%</td>
        <td class="fn-send">${(v.hom+v.het)>0
          ? `<button class="btn tiny" title="Open ${esc(d.gene)} in SNPVersity with all ${v.hom+v.het} carriers of this allele preselected"
               onclick="event.stopPropagation();FUNCTION.toVersity('${esc(v.id)}','all')">SNPVersity →</button>`
          : '<span class="muted">—</span>'}</td>
      </tr>${open?carrierRow(v,sec,d):''}`;
    }).join('');
    const anyCarrier = d.damaging.some(v=>(v.hom+v.het)>0);
    return `<div class="card pad">
      <div class="fn-h" style="display:flex;align-items:center;gap:10px;flex-wrap:wrap">Damaging &amp; knockout alleles
        <span class="muted" style="font-weight:400;font-size:12px">${d.damaging.length} alleles · click a row for carrier lines</span>
        <span style="margin-left:auto;display:flex;gap:8px;flex-wrap:wrap">
          ${anyCarrier?`<button class="btn" title="Open this gene's region in SNPVersity with every accession carrying any damaging allele preselected"
            onclick="FUNCTION.toVersity('all','all')">${ICONS.dna||''} Send all carriers to SNPVersity</button>`:''}
          <button class="btn" onclick="FUNCTION.exportCSV()">${ICONS.download||''} Export CSV</button>
        </span>
      </div>
      <div class="tbl-wrap" style="max-height:none"><table class="vcf imp">
        <thead><tr><th style="padding-left:11px">Allele</th><th>Consequence</th><th>Domain</th><th class="num">PlantCAD1</th>${sec?'<th class="num">PlantCAD2</th>':''}<th class="num">ESM</th>${sec?'<th class="num">ESM2</th><th class="num">ESM3</th>':''}<th>Priority</th><th class="num" title="heterozygous carriers">Het</th><th class="num" title="homozygous carriers">Hom</th><th class="num">AF</th><th></th></tr></thead>
        <tbody>${rows}</tbody>
      </table></div>
    </div>`;
  }
  function carrierRow(v, sec, d){
    const chip = (id,cls)=>`<span class="carrier ${cls}">${esc(id)}</span>`;
    const homs = v.carriersHom.slice(0,60).map(id=>chip(id,'hom')).join('');
    const hets = v.carriersHet.slice(0,60).map(id=>chip(id,'het')).join('');
    const send = (mode,label,n)=> n
      ? `<button class="btn tiny" onclick="event.stopPropagation();FUNCTION.toVersity('${esc(v.id)}','${mode}')">${label} (${n}) →</button>`
      : '';
    return `<tr class="fn-carriers"><td colspan="${sec?13:10}">
      <div class="fn-cwrap">
        <div><div class="fn-k">Homozygous ${v.consClass==='lof'?'(candidate knockouts)':''} · ${v.carriersHom.length}</div>
          <div class="fn-chips">${homs||'<span class="muted">none</span>'}${v.carriersHom.length>60?` <span class="muted">+${v.carriersHom.length-60} more</span>`:''}</div></div>
        <div style="margin-top:8px"><div class="fn-k">Heterozygous · ${v.carriersHet.length}</div>
          <div class="fn-chips">${hets||'<span class="muted">none</span>'}${v.carriersHet.length>60?` <span class="muted">+${v.carriersHet.length-60} more</span>`:''}</div></div>
        <div class="fn-sendrow">
          <span class="fn-k" style="margin:0">Open in SNPVersity</span>
          ${send('hom','Homozygous carriers',v.carriersHom.length)}
          ${send('het','Heterozygous carriers',v.carriersHet.length)}
          ${send('all','All carriers',v.carriersHom.length+v.carriersHet.length)}
        </div>
      </div></td></tr>`;
  }

  /* =====================================================================
   *  FUNCTIONAL ANNOTATION  (per-gene JSON from ANN_DIR)
   *  A gene-scoped functional dossier: identity + description, protein
   *  domain architecture (to-scale diagram), Gene Ontology grouped by
   *  aspect, KEGG pathways / orthology, and external cross-references.
   * ===================================================================== */

  /* small provenance/aspect palettes ------------------------------------ */
  const GO_ASPECT = {
    BP:{label:'Biological process', abbr:'BP'},
    MF:{label:'Molecular function', abbr:'MF'},
    CC:{label:'Cellular component', abbr:'CC'},
  };
  const DOM_PALETTE = ['#2f6ad0','#176c3a','#b8862b','#7a3fb0','#0e7c86','#c0362c','#3a5a99'];

  /* provenance badge for a functional_description.source string */
  function descSourceBadge(src){
    src = src || 'none';
    let cls='pred', label=src;
    if (/^MaizeGDB/i.test(src)){ cls='curated'; label='MaizeGDB · curated'; }
    else if (/^UniProt/i.test(src)){ cls='uniprot'; label='UniProt'; }
    else if (/^InterPro/i.test(src)){ cls='pred'; label='InterPro · predicted'; }
    else if (src==='none'){ cls='none'; label='no informative source'; }
    return `<span class="ann-srcbadge ${cls}" title="functional_description.source = ${esc(src)}">${esc(label)}</span>`;
  }
  /* one small badge per GO evidence source */
  function goSourceBadges(sources){
    const meta = {MaizeGDB:'curated', UniProt:'uniprot', InterPro2GO:'pred'};
    return (sources||[]).map(s=>`<span class="go-src ${meta[s]||'pred'}" title="${esc(s)}${s==='InterPro2GO'?' (predicted from domain)':''}">${esc(s)}</span>`).join('');
  }
  function isPredictedOnly(t){ const s=t.sources||[]; return s.length>0 && s.every(x=>x==='InterPro2GO'); }
  function xrefChip(text, href, title){
    return href
      ? `<a class="xref" href="${href}" target="_blank" rel="noopener" title="${esc(title||text)}">${esc(text)} <span class="xref-ext">↗</span></a>`
      : `<span class="xref" title="${esc(title||text)}">${esc(text)}</span>`;
  }

  /* ---------- top-level annotation section ---------- */
  function annotationSection(){
    const a = FN.annotation;
    if (a === undefined){
      return `<div class="card pad" style="margin-bottom:16px"><div class="fn-h">Functional annotation</div>
        <div class="ann-load"><div class="spinner sm"></div><span>Loading functional annotation…</span></div></div>`;
    }
    if (a === null){
      return `<div class="card pad" style="margin-bottom:16px"><div class="fn-h">Functional annotation</div>
        <div class="muted" style="padding:4px 0">No functional-annotation record for <span class="mono">${esc(FN.gene)}</span>.
        Records exist for the 39,756 canonical B73 v5 gene models (<span class="mono">Zm00001eb…</span>).</div></div>`;
    }
    return annHeader(a) + annDomains(a) + annGO(a) + annPathways(a) + annXrefs(a);
  }

  /* ---------- identity header: symbol, name, description, evidence ---------- */
  function annHeader(a){
    const nm = a.names||{}, fd = a.functional_description||{};
    const symbol = nm.primary_symbol || '';
    const fullName = nm.full_name || (a.functional_description && a.functional_description.text) || '';
    const aliases = (nm.aliases||[]).filter(Boolean);
    const cov = a.annotation_coverage||{};
    // "evidence present" chips — an at-a-glance summary of what's annotated
    const ev = [
      ['GO',      cov.has_go],
      ['Pfam',    cov.has_pfam],
      ['KEGG KO', cov.has_ko],
      ['Pathway', cov.has_pathway],
      ['UniProt', cov.has_uniprot],
      ['Symbol',  cov.has_symbol],
    ].map(([k,on])=>`<span class="ev-chip ${on?'on':'off'}" title="${on?'present':'not available'}">${on?'✓':'·'} ${k}</span>`).join('');
    const ko = ((a.kegg&&a.kegg.orthology)||[])[0];
    return `<div class="card pad ann-head" style="margin-bottom:16px">
      <div class="ann-idrow">
        ${symbol?`<span class="ann-symbol" title="primary gene symbol">${esc(symbol)}</span>`:''}
        <div class="ann-titlewrap">
          <div class="ann-fullname">${esc(fullName||'Uncharacterized protein')}</div>
          ${aliases.length?`<div class="ann-aliases">also: ${aliases.map(x=>`<span class="mono">${esc(x)}</span>`).join(', ')}</div>`:''}
        </div>
      </div>
      <div class="ann-desc">
        <span class="ann-desc-txt">${esc(fd.text||'Uncharacterized protein')}</span>
        ${descSourceBadge(fd.source)}
      </div>
      ${ko?`<div class="ann-ko"><span class="fn-k" style="display:inline">KEGG orthology</span> ${esc(ko.ko_id)} · ${esc(ko.name)}</div>`:''}
      <div class="ann-ev">${ev}</div>
    </div>`;
  }

  /* ---------- protein domain architecture (to-scale diagram) ---------- */
  function annDomains(a){
    const prot = a.protein||{}, L = +prot.length||0;
    const doms = (a.pfam_domains||[]).slice().sort((x,y)=>(x.start||0)-(y.start||0));
    const svg = domainSVG(L, doms);
    const list = doms.length ? doms.map((x,i)=>{
      const col = DOM_PALETTE[i%DOM_PALETTE.length];
      const span = (x.start!=null&&x.end!=null)?`${x.start}–${x.end} aa`:'—';
      const cover = (L&&x.start!=null&&x.end!=null)?` · ${Math.round(100*(x.end-x.start+1)/L)}% of protein`:'';
      return `<div class="dom-row">
        <span class="dom-swatch" style="background:${col}"></span>
        <span class="dom-name">${esc(x.pfam_name||x.pfam_id||'domain')}</span>
        <span class="mono dom-ids">${esc(x.pfam_id||'')}${x.interpro_id?` · ${esc(x.interpro_id)}`:''}</span>
        <span class="dom-span mono">${span}${cover}</span>
        ${x.evalue?`<span class="dom-ev mono" title="InterProScan E-value">E=${esc(x.evalue)}</span>`:''}
      </div>`;
    }).join('') : `<div class="muted" style="padding:4px 0">No Pfam domains annotated for this protein (Pfam-only InterProScan input).</div>`;
    return `<div class="card pad" style="margin-bottom:16px">
      <div class="fn-h" style="display:flex;align-items:baseline;gap:10px">Protein domain architecture
        <span class="muted" style="font-weight:400;font-size:12px">${prot.protein_id?esc(prot.protein_id)+' · ':''}${L?L+' aa':'length n/a'}</span>
      </div>
      <div class="dom-arch">${svg}</div>
      <div class="dom-list">${list}</div>
    </div>`;
  }

  /* SVG: backbone scaled to protein length, Pfam domains as positioned blocks */
  function domainSVG(L, doms){
    const W=1000, padL=10, padR=10, trackW=W-padL-padR;
    const sc = x => padL + (L ? (x/L)*trackW : 0);
    // greedy row packing so overlapping domains don't collide
    const rowsEnd=[]; const placed=doms.map(d=>{
      let r=0; while(rowsEnd[r]!=null && rowsEnd[r] > (d.start||0)) r++;
      rowsEnd[r]=(d.end||0); return {d,row:r};
    });
    const nRows=Math.max(1,rowsEnd.length);
    const rowH=30, gap=8, top=30, bot=26;
    const H=top+nRows*rowH+(nRows-1)*gap+bot;
    const trackMid=(row)=>top+row*(rowH+gap);
    // ruler ticks at 0, 25%, 50%, 75%, 100%
    let ticks='';
    if (L){
      for (let f=0; f<=1.0001; f+=0.25){
        const x=sc(L*f), val=Math.round(L*f);
        ticks += `<line x1="${x}" y1="${top-8}" x2="${x}" y2="${top+nRows*rowH+(nRows-1)*gap+4}" class="dom-tick"/>
                  <text x="${x}" y="${top-12}" class="dom-tick-lbl" text-anchor="${f===0?'start':f>=1?'end':'middle'}">${val}</text>`;
      }
    }
    const backbone = `<rect x="${padL}" y="${top+ (nRows*rowH+(nRows-1)*gap)/2 - 4}" width="${trackW}" height="8" rx="4" class="dom-backbone"/>`;
    const blocks = placed.map((p,i)=>{
      const d=p.d, x=sc(d.start||0), w=Math.max(3, sc(d.end||0)-sc(d.start||0));
      const y=trackMid(p.row), col=DOM_PALETTE[i%DOM_PALETTE.length];
      const label=d.pfam_name||d.pfam_id||'';
      const showText = w > label.length*7.5 + 12;   // only inline the label if it fits
      const tip=`${label}${d.pfam_id?' ('+d.pfam_id+')':''} · ${d.start}-${d.end} aa${d.evalue?' · E='+d.evalue:''}`;
      return `<g class="dom-block"><title>${esc(tip)}</title>
        <rect x="${x}" y="${y}" width="${w}" height="${rowH-6}" rx="5" fill="${col}"/>
        ${showText?`<text x="${x+w/2}" y="${y+(rowH-6)/2}" class="dom-blk-lbl" text-anchor="middle" dominant-baseline="central">${esc(label)}</text>`:''}
      </g>`;
    }).join('');
    return `<svg viewBox="0 0 ${W} ${H}" class="dom-svg" preserveAspectRatio="xMidYMid meet" role="img"
      aria-label="Protein domain architecture, ${L} amino acids, ${doms.length} Pfam domain(s)">
      ${ticks}${backbone}${blocks}
    </svg>`;
  }

  /* ---------- Gene Ontology grouped by aspect ---------- */
  function annGO(a){
    const all = a.go_terms||[];
    if (!all.length){
      return `<div class="card pad" style="margin-bottom:16px"><div class="fn-h">Gene Ontology</div>
        <div class="muted" style="padding:4px 0">No GO terms annotated for this gene.</div></div>`;
    }
    const curatedOnly = FN.goCurated;
    const terms = curatedOnly ? all.filter(t=>!isPredictedOnly(t)) : all;
    // aspect counts (of the currently shown set) for the mini summary bar
    const counts = {BP:0,MF:0,CC:0};
    terms.forEach(t=>{ if(counts[t.aspect]!=null) counts[t.aspect]++; });
    const shown = terms.length || 1;
    const summary = ['BP','MF','CC'].map(k=>counts[k]?`<span class="go-seg asp-${k}" style="width:${100*counts[k]/shown}%" title="${GO_ASPECT[k].label}: ${counts[k]}"></span>`:'').join('');
    const legend = ['BP','MF','CC'].map(k=>`<span class="go-lg"><span class="go-dot asp-${k}"></span>${GO_ASPECT[k].abbr} ${counts[k]}</span>`).join('');

    const group = k => {
      const list = terms.filter(t=>t.aspect===k);
      if (!list.length) return '';
      const chips = list.map(t=>{
        const obs = t.status && t.status!=='current';
        return `<div class="go-chip asp-${k} ${obs?'obs':''}" title="${esc(t.go_id)}${obs?' · '+esc(t.status):''}">
          <a class="go-id" href="https://amigo.geneontology.org/amigo/term/${encodeURIComponent(t.go_id)}" target="_blank" rel="noopener">${esc(t.go_id)}</a>
          <span class="go-name">${esc(t.name)}</span>
          <span class="go-srcs">${goSourceBadges(t.sources)}</span>
          ${obs?`<span class="go-obs" title="obsolete term retained/flagged">obsolete</span>`:''}
        </div>`;
      }).join('');
      return `<div class="go-group">
        <div class="go-ghead"><span class="go-dot asp-${k}"></span>${GO_ASPECT[k].label} <span class="muted">· ${list.length}</span></div>
        <div class="go-chips">${chips}</div>
      </div>`;
    };
    const nPred = all.filter(isPredictedOnly).length;
    return `<div class="card pad" style="margin-bottom:16px">
      <div class="fn-h" style="display:flex;align-items:center;gap:10px">Gene Ontology
        <span class="muted" style="font-weight:400;font-size:12px">${terms.length} of ${all.length} terms${curatedOnly?' · curated only':''}</span>
        ${nPred?`<button class="btn" style="margin-left:auto;font-size:12px;padding:6px 11px" onclick="FUNCTION.toggleGO()">${curatedOnly?'Show all evidence':'Curated only'}</button>`:''}
      </div>
      <div class="go-summary"><div class="go-bar">${summary}</div><div class="go-legend">${legend}</div></div>
      ${group('BP')}${group('MF')}${group('CC')}
      <div class="go-note muted">Source confidence: <span class="go-src curated">MaizeGDB</span>/<span class="go-src uniprot">UniProt</span> are curated; <span class="go-src pred">InterPro2GO</span> is predicted from domains.</div>
    </div>`;
  }

  /* ---------- KEGG pathways & orthology ---------- */
  function annPathways(a){
    const paths = a.pathways||[];
    const kegg = a.kegg||{};
    const kos = kegg.orthology||[];
    const genes = kegg.kegg_gene_ids||[];
    if (!paths.length && !kos.length && !genes.length){
      return `<div class="card pad" style="margin-bottom:16px"><div class="fn-h">Pathways &amp; orthology</div>
        <div class="muted" style="padding:4px 0">No KEGG pathway or orthology mapping for this gene. KEGG coverage tracks UniProt/Entrez cross-references (~5–10% of models).</div></div>`;
    }
    const koHTML = kos.length ? `<div class="kv"><div class="fn-k">KEGG orthology (KO)</div><div class="fn-chips">${
      kos.map(k=>xrefChip(`${k.ko_id} · ${k.name}`, `https://www.kegg.jp/entry/${encodeURIComponent(k.ko_id)}`, k.name)).join('')
    }</div></div>` : '';
    const pathHTML = paths.length ? `<div class="kv"><div class="fn-k">Pathways</div><div class="fn-chips">${
      paths.map(p=>xrefChip(`${p.kegg_pathway_id||''} ${p.name||''}`.trim(), p.kegg_pathway_id?`https://www.kegg.jp/pathway/${encodeURIComponent(p.kegg_pathway_id)}`:null, p.name)).join('')
    }</div></div>` : '';
    const geneHTML = genes.length ? `<div class="kv"><div class="fn-k">KEGG genes</div><div class="fn-chips">${
      genes.map(g=>xrefChip(g, `https://www.kegg.jp/entry/${encodeURIComponent(g)}`)).join('')
    }</div></div>` : '';
    return `<div class="card pad" style="margin-bottom:16px">
      <div class="fn-h">Pathways &amp; orthology</div>
      ${koHTML}${pathHTML}${geneHTML}
    </div>`;
  }

  /* ---------- external cross-references ---------- */
  function annXrefs(a){
    const x = a.cross_references||{};
    const rows = [];
    if ((x.uniprot||[]).length) rows.push(['UniProt', x.uniprot.map(u=>xrefChip(u, `https://www.uniprot.org/uniprotkb/${encodeURIComponent(u)}/entry`))]);
    if ((x.ncbi_entrez||[]).length) rows.push(['NCBI Gene', x.ncbi_entrez.map(g=>xrefChip(g, `https://www.ncbi.nlm.nih.gov/gene/${encodeURIComponent(g)}`))]);
    if ((x.gene_model_v4||[]).length) rows.push(['B73 v4 model', x.gene_model_v4.map(g=>xrefChip(g, `https://maizegdb.org/gene_center/gene/${encodeURIComponent(g)}`))]);
    if ((x.gene_model_v3||[]).length) rows.push(['B73 v3 model', x.gene_model_v3.map(g=>xrefChip(g, `https://maizegdb.org/gene_center/gene/${encodeURIComponent(g)}`))]);
    if (!rows.length) return '';
    const build = a.build_date ? `<div class="ann-build muted">Annotation build ${esc(a.build_date)} · ${esc((a.assembly||''))} ${esc((a.annotation_version||''))}</div>` : '';
    return `<div class="card pad" style="margin-bottom:16px">
      <div class="fn-h">Cross-references</div>
      ${rows.map(([k,chips])=>`<div class="kv"><div class="fn-k">${esc(k)}</div><div class="fn-chips">${chips.join('')}</div></div>`).join('')}
      ${build}
    </div>`;
  }

  /* ---------- public handlers ---------- */
  window.FUNCTION = {
    load(){ const el=document.getElementById('fnGeneInput'); if(!el)return; const g=el.value.trim(); if(!g)return;
      FN.gene=g; FN.data=null; FN.openId=null; analyzeGene(); },
    toggle(id){ FN.openId = (FN.openId===id?null:id); paint(); },
    panEffect(id){ const d=FN.data; if(!d||!d.damaging) return;
      panEffectTo(d.damaging.find(v=>String(v.id)===String(id))); },
    fold(id){ const d=FN.data; if(!d||!d.damaging) return;
      foldTo(d.damaging.find(v=>String(v.id)===String(id))); },

    /* Hand this gene's region + the accessions carrying an alternative allele
       over to SNPVersity.  alleleId 'all' = every damaging allele in the table.
       mode: 'hom' | 'het' | 'all' (homozygous, heterozygous, or both). */
    toVersity(alleleId, mode){
      const d=FN.data; if(!d||!d.damaging) return;
      mode = mode || 'all';
      const vs = alleleId==='all' ? d.damaging : d.damaging.filter(v=>String(v.id)===String(alleleId));
      if(!vs.length) return;
      const acc=new Set();
      vs.forEach(v=>{
        if(mode!=='het') (v.carriersHom||[]).forEach(a=>acc.add(a));
        if(mode!=='hom') (v.carriersHet||[]).forEach(a=>acc.add(a));
      });
      if(!acc.size){ alert('No carriers to send for this allele.'); return; }
      const what = mode==='hom' ? 'homozygous carriers'
                 : mode==='het' ? 'heterozygous carriers'
                 : 'carriers of an alternative allele';
      const note = alleleId==='all'
        ? `${what} across ${vs.length} damaging allele${vs.length>1?'s':''}`
        : `${what} of ${vs[0].variant}`;
      const payload = {
        gene:d.gene, chr:d.chr, start:+d.start, end:+d.end,
        dataset:d.dataset!=null?d.dataset:FN.dataset,
        accessions:[...acc], from:'SNPFunction', note,
        allele: alleleId==='all' ? null : vs[0].variant, mode
      };
      if (typeof window.versityRequest === 'function'){ window.versityRequest(payload); return; }
      // fallback if snpversity.js is an older build without the inbound hook
      if (typeof S !== 'undefined' && S){
        S.pendingVersity = payload;
        if (typeof go === 'function') go('snpversity');
      }
    },
    toggleGO(){ FN.goCurated = !FN.goCurated; paint(); },
    /* compact dataset chooser — selecting a dataset only records the choice + moves the
       highlight. The page is NOT re-analyzed here; the new dataset is applied on the next
       "Analyze gene" click (FUNCTION.load), which reads the current FN.dataset. */
    pickDataset(id){
      const ds = fnDatasets().find(d=>String(datasetId(d))===String(id));
      const val = ds ? datasetId(ds) : id;
      if (isCurrentDataset(ds || val)){ syncDatasetChooser(); return; }   // no change
      FN.dataset = val;
      if (typeof S !== 'undefined' && S) S.dataset = val;   // keep the whole app in sync
      syncDatasetChooser();
    },
    exportCSV(){
      const d=FN.data; if(!d||!d.damaging) return;
      const sec=Data.hasSecondaryScores(d.dataset);
      const cols=['variant','consequence','domain','plantcad']
        .concat(sec?['plantcad2']:[])
        .concat(['esm'])
        .concat(sec?['esm2','esm3']:[])
        .concat(['combined','priority','het','hom','af','homozygous_carriers','het_carriers']);
      const line=v=>{
        const base=[v.variant,v.consequence,v.domain,v.plantcad];
        if(sec) base.push(v.plantcad2);
        base.push(v.esm);
        if(sec) base.push(v.esm2, v.esm3);
        base.push(v.combined,v.priority,v.het,v.hom,v.af.toFixed(4),
          '"'+v.carriersHom.join(';')+'"','"'+v.carriersHet.join(';')+'"');
        return base.map(x=>x==null?'':x).join(',');
      };
      const csv=[cols.join(',')].concat(d.damaging.map(line)).join('\n');
      const blob=new Blob([csv],{type:'text/csv'}), u=URL.createObjectURL(blob), a=document.createElement('a');
      a.href=u; a.download=`snpfunction_${d.gene}_damaging.csv`; document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(u),1500);
    },
  };

  function injectCSS(){
    if (document.getElementById('snpfunction-css')) return;
    const s=document.createElement('style'); s.id='snpfunction-css';
    s.textContent=`
      .fn-grid{display:flex;gap:26px;flex-wrap:wrap}
      .fn-k{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:4px}
      .fn-v{font-size:14px;color:var(--ink)} .fn-v .dom-tag{margin:0 4px 4px 0;display:inline-block}
      .fn-h{font-weight:600;font-size:15px;margin-bottom:12px;color:var(--ink)}
      .fn-stats{display:flex;gap:14px;flex-wrap:wrap;margin-bottom:14px}
      .fn-stat{background:#f5f8fc;border:1px solid var(--line);border-radius:10px;padding:10px 14px;min-width:120px}
      .fn-statv{font-size:20px;font-weight:700;color:var(--ink);font-family:var(--mono)}
      .fn-statk{font-size:11px;color:var(--muted);margin-top:2px}
      .fn-bar{display:flex;height:16px;border-radius:8px;overflow:hidden;background:#eef1f5}
      .fn-seg{display:inline-block;height:100%}
      .fn-seg.lof,.fn-sw.lof{background:#c0362c}.fn-seg.splice,.fn-sw.splice{background:#b8862b}
      .fn-seg.missense,.fn-sw.missense{background:#2f6ad0}.fn-seg.indel,.fn-sw.indel{background:#176c3a}
      .fn-seg.syn,.fn-sw.syn{background:#9aa7bd}.fn-seg.other,.fn-sw.other{background:#cdd6e3}
      .fn-legend{display:flex;gap:14px;flex-wrap:wrap;margin-top:8px;font-size:12px;color:var(--muted)}
      .fn-lg{display:inline-flex;align-items:center;gap:5px} .fn-sw{width:10px;height:10px;border-radius:3px;display:inline-block}
      .fn-af{margin-top:12px;font-size:13px;color:var(--muted)} .fn-af b{color:var(--ink)}
      .fn-carriers td{background:#fbfcfe;border-bottom:1px solid var(--line)}
      .fn-cwrap{padding:10px 12px} .fn-chips{display:flex;flex-wrap:wrap;gap:5px;margin-top:5px}
      .carrier{font-family:var(--mono);font-size:11px;padding:2px 7px;border-radius:6px;border:1px solid var(--line)}
      .carrier.hom{background:#fdecea;border-color:#f0c4bd;color:#8f281c}
      .carrier.het{background:#eef4ff;border-color:#cfe0ff;color:#274b8f}
      .fn-sendrow{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-top:12px;
        border-top:1px solid var(--line);padding-top:10px}
      .btn.tiny{font-size:11px;padding:4px 9px;border-radius:7px;line-height:1.3;white-space:nowrap}
      td.fn-send{text-align:right;padding-right:10px;white-space:nowrap}
      /* per-row jump links in the Consequence column */
      table.imp .fold-jump{font-size:11px;white-space:nowrap;margin-left:4px;
        color:#176c3a;text-decoration:none;border-bottom:1px dotted #9ecdb1}
      table.imp .fold-jump:hover{color:#0f4f2a;border-bottom-style:solid}
      /* compact dataset chooser (shares Data.datasets() with SNPVersity) */
      .fn-ds-card{padding:12px 14px}
      .fn-ds-head{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:9px}
      .fn-ds-grid{display:flex;flex-wrap:wrap;gap:10px}
      .fn-ds{display:flex;align-items:center;gap:9px;text-align:left;cursor:pointer;background:#fff;
        border:1px solid var(--line);border-radius:10px;padding:9px 13px;min-width:210px;transition:border-color .12s,box-shadow .12s}
      .fn-ds:hover{border-color:#c3cee0}
      .fn-ds.sel{border-color:#9db4dd;box-shadow:0 0 0 2px rgba(47,106,208,.12);background:#f7faff}
      .fn-ds-dot{flex:0 0 auto;width:10px;height:10px;border-radius:50%;border:2px solid #c3cee0;background:#fff}
      .fn-ds.sel .fn-ds-dot{border-color:#2f6ad0;background:#2f6ad0;box-shadow:inset 0 0 0 2px #fff}
      .fn-ds-txt{display:flex;flex-direction:column;line-height:1.25;min-width:0}
      .fn-ds-name{font-weight:600;font-size:13px;color:var(--ink)}
      .fn-ds-sub{font-size:11.5px;color:var(--muted)}

      /* ===== functional annotation ===== */
      .ann-load{display:flex;align-items:center;gap:10px;color:var(--muted);padding:6px 0}
      .spinner.sm{width:16px;height:16px;border-width:2px}
      .kv{margin-top:12px} .kv:first-of-type{margin-top:4px}
      .fn-chips{display:flex;flex-wrap:wrap;gap:6px;margin-top:5px}
      /* identity header */
      .ann-idrow{display:flex;align-items:center;gap:14px;flex-wrap:wrap}
      .ann-symbol{font-family:var(--mono);font-weight:700;font-size:17px;color:#fff;background:var(--mg-blue,rgb(8,44,148));
        padding:5px 12px;border-radius:9px;letter-spacing:.02em;white-space:nowrap}
      .ann-titlewrap{min-width:0}
      .ann-fullname{font-size:17px;font-weight:600;color:var(--ink);line-height:1.2}
      .ann-aliases{font-size:12px;color:var(--muted);margin-top:2px}
      .ann-desc{margin-top:12px;display:flex;align-items:center;gap:10px;flex-wrap:wrap;font-size:13.5px;color:var(--ink)}
      .ann-ko{margin-top:8px;font-size:12.5px;color:var(--muted);font-family:var(--mono)}
      .ann-srcbadge{font:600 10.5px/1 var(--body,'Inter',sans-serif);padding:4px 8px;border-radius:20px;white-space:nowrap;text-transform:uppercase;letter-spacing:.3px}
      .ann-srcbadge.curated{background:#e7f3ec;color:#176c3a;border:1px solid #bfe0cb}
      .ann-srcbadge.uniprot{background:#eaf1fc;color:#274b8f;border:1px solid #cbdcf6}
      .ann-srcbadge.pred{background:#f2f0ea;color:#7a5b12;border:1px solid #e2dcc9}
      .ann-srcbadge.none{background:#f2f3f5;color:#7a828f;border:1px solid #e2e6ec}
      .ann-ev{display:flex;flex-wrap:wrap;gap:6px;margin-top:13px}
      .ev-chip{font:600 11px/1 var(--body,'Inter',sans-serif);padding:5px 9px;border-radius:7px;border:1px solid var(--line)}
      .ev-chip.on{background:#f0f7f2;color:#1c6b3c;border-color:#cfe6d7}
      .ev-chip.off{background:#f7f8fa;color:#aab2be}
      /* domain architecture */
      .dom-arch{margin:6px 0 4px} .dom-svg{width:100%;height:auto;display:block;overflow:visible}
      .dom-backbone{fill:#dfe4ec}
      .dom-tick{stroke:#e3e8ef;stroke-width:1}
      .dom-tick-lbl{fill:#98a1af;font-family:var(--mono);font-size:12px}
      .dom-blk-lbl{fill:#fff;font-family:var(--body,'Inter',sans-serif);font-size:12px;font-weight:600;pointer-events:none}
      .dom-block rect{transition:opacity .12s} .dom-block:hover rect{opacity:.85;cursor:default}
      .dom-list{margin-top:10px;border-top:1px solid var(--line);padding-top:10px;display:flex;flex-direction:column;gap:7px}
      .dom-row{display:flex;align-items:center;gap:9px;flex-wrap:wrap;font-size:13px}
      .dom-swatch{width:11px;height:11px;border-radius:3px;flex:0 0 auto}
      .dom-name{font-weight:600;color:var(--ink)}
      .dom-ids{font-size:11.5px;color:var(--muted)}
      .dom-span{font-size:11.5px;color:var(--muted);margin-left:auto}
      .dom-ev{font-size:11px;color:#8a6d1e;background:#f6f1e2;border:1px solid #e6dcc2;border-radius:6px;padding:1px 6px}
      /* gene ontology */
      .go-summary{margin:2px 0 14px}
      .go-bar{display:flex;height:12px;border-radius:7px;overflow:hidden;background:#eef1f5}
      .go-seg{display:inline-block;height:100%}
      .go-legend{display:flex;gap:14px;flex-wrap:wrap;margin-top:8px;font-size:12px;color:var(--muted)}
      .go-lg{display:inline-flex;align-items:center;gap:5px}
      .go-dot{width:10px;height:10px;border-radius:50%;display:inline-block}
      .go-group{margin-top:14px}
      .go-ghead{font-size:12.5px;font-weight:600;color:var(--ink);display:flex;align-items:center;gap:7px;margin-bottom:8px}
      .go-chips{display:flex;flex-wrap:wrap;gap:8px}
      .go-chip{display:inline-flex;align-items:center;gap:8px;border:1px solid var(--line);border-left-width:4px;
        border-radius:9px;padding:6px 10px;background:#fff;font-size:12.5px;max-width:100%}
      .go-chip.obs{opacity:.6}
      .go-id{font-family:var(--mono);font-size:11.5px;color:var(--muted);text-decoration:none;white-space:nowrap}
      .go-id:hover{text-decoration:underline}
      .go-name{color:var(--ink)}
      .go-srcs{display:inline-flex;gap:4px}
      .go-src{font:600 9.5px/1 var(--body,'Inter',sans-serif);padding:3px 6px;border-radius:5px;text-transform:uppercase;letter-spacing:.3px;white-space:nowrap}
      .go-src.curated{background:#e7f3ec;color:#176c3a} .go-src.uniprot{background:#eaf1fc;color:#274b8f} .go-src.pred{background:#f2f0ea;color:#7a5b12}
      .go-obs{font-size:10px;color:#a23b2c;background:#fbeae7;border-radius:5px;padding:2px 6px}
      .go-note{font-size:11.5px;margin-top:14px;display:flex;align-items:center;gap:6px;flex-wrap:wrap}
      /* aspect colors: BP green, MF blue, CC purple */
      .asp-BP{--asp:#176c3a} .asp-MF{--asp:#2f6ad0} .asp-CC{--asp:#7a3fb0}
      .go-seg.asp-BP{background:#176c3a}.go-seg.asp-MF{background:#2f6ad0}.go-seg.asp-CC{background:#7a3fb0}
      .go-dot.asp-BP{background:#176c3a}.go-dot.asp-MF{background:#2f6ad0}.go-dot.asp-CC{background:#7a3fb0}
      .go-chip.asp-BP{border-left-color:#176c3a}.go-chip.asp-MF{border-left-color:#2f6ad0}.go-chip.asp-CC{border-left-color:#7a3fb0}
      /* cross-references */
      .xref{display:inline-flex;align-items:center;gap:4px;font-family:var(--mono);font-size:12px;padding:4px 9px;border-radius:7px;
        border:1px solid var(--line);background:#f7f9fc;color:var(--ink);text-decoration:none}
      a.xref:hover{border-color:#c3cee0;background:#eef4ff}
      .xref-ext{color:var(--muted);font-family:var(--body,'Inter',sans-serif)}
      .ann-build{font-size:11px;margin-top:12px}
      @media (max-width:640px){ .dom-span{margin-left:0} }`;
    document.head.appendChild(s);
  }

  SNPTools.register('snpfunction', { render });
})();
