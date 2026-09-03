/* =====================================================================
 *  core.js — suite shell: icons, tool registry, shared state, router,
 *  nav, placeholder pages, tooltips, menus. Loaded FIRST.
 *  Tool files (snpversity.js, snpimpact.js, …) register themselves with
 *  SNPTools.register(id, { render(page){…} }).  data.js owns all data.
 * ===================================================================== */

/* ================= ICONS ================= */
const ICONS = {
  dna:'<svg viewBox="0 0 24 24" fill="none"><path d="M6 3c3 4 9 4 12 9M6 12c3 4 9 4 12 9" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"/><path d="M7 6h6M9 9h6M9 15h6M11 18h6" stroke="currentColor" stroke-width="1.4" stroke-linecap="round" opacity=".6"/></svg>',
  leaf:'<svg viewBox="0 0 24 24" fill="none"><path d="M5 19c0-8 6-14 14-14 0 8-6 14-14 14z" stroke="currentColor" stroke-width="1.7" stroke-linejoin="round"/><path d="M5 19C9 14 13 11 17 9" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"/></svg>',
  star:'<svg viewBox="0 0 24 24" fill="none"><circle cx="12" cy="12" r="9" stroke="currentColor" stroke-width="1.7"/><path d="M12 7.5l1.3 2.7 3 .4-2.2 2.1.5 3-2.6-1.4-2.6 1.4.5-3-2.2-2.1 3-.4z" fill="currentColor"/></svg>',
  func:'<svg viewBox="0 0 24 24" fill="none"><rect x="3" y="4" width="18" height="14" rx="2" stroke="currentColor" stroke-width="1.7"/><path d="M7 9h6M7 13h10" stroke="currentColor" stroke-width="1.6" stroke-linecap="round"/><circle cx="17" cy="9" r="1.4" fill="currentColor"/></svg>',
  compare:'<svg viewBox="0 0 24 24" fill="none"><rect x="3" y="5" width="7" height="14" rx="1.5" stroke="currentColor" stroke-width="1.7"/><rect x="14" y="5" width="7" height="14" rx="1.5" stroke="currentColor" stroke-width="1.7"/></svg>',
  tree:'<svg viewBox="0 0 24 24" fill="none"><path d="M5 12h5m4 0h5M10 6v12M19 6v12" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"/><circle cx="5" cy="12" r="2" stroke="currentColor" stroke-width="1.6"/><circle cx="19" cy="6" r="2" stroke="currentColor" stroke-width="1.6"/><circle cx="19" cy="18" r="2" stroke="currentColor" stroke-width="1.6"/></svg>',
  grid:'<svg viewBox="0 0 24 24" fill="none"><rect x="4" y="4" width="6" height="6" rx="1" stroke="currentColor" stroke-width="1.7"/><rect x="14" y="4" width="6" height="6" rx="1" stroke="currentColor" stroke-width="1.7"/><rect x="4" y="14" width="6" height="6" rx="1" stroke="currentColor" stroke-width="1.7"/><rect x="14" y="14" width="6" height="6" rx="1" stroke="currentColor" stroke-width="1.7"/></svg>',
  fold:'<svg viewBox="0 0 24 24" fill="none"><path d="M7 4c4 2 4 6 0 8s-4 6 0 8M17 4c-4 2-4 6 0 8s4 6 0 8" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"/></svg>',
  density:'<svg viewBox="0 0 24 24" fill="none"><path d="M4 18V9M9 18V5M14 18v-7M19 18v-4" stroke="currentColor" stroke-width="2" stroke-linecap="round"/></svg>',
  germ:'<svg viewBox="0 0 24 24" fill="none"><path d="M12 21v-7" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"/><path d="M12 14c0-3-2-5-5-5 0 3 2 5 5 5zM12 12c0-3 2-5 5-5 0 3-2 5-5 5z" stroke="currentColor" stroke-width="1.6" stroke-linejoin="round"/></svg>',
  impute:'<svg viewBox="0 0 24 24" fill="none"><path d="M20 12a8 8 0 10-2.3 5.6" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"/><path d="M20 7v5h-5" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"/></svg>',
  search:'<svg viewBox="0 0 24 24" fill="none"><circle cx="11" cy="11" r="6" stroke="currentColor" stroke-width="2"/><path d="M16 16l4 4" stroke="currentColor" stroke-width="2" stroke-linecap="round"/></svg>',
  check:'<svg viewBox="0 0 24 24" fill="none"><path d="M5 12l4 4 10-10" stroke="currentColor" stroke-width="3" stroke-linecap="round" stroke-linejoin="round"/></svg>',
  caret:'<svg viewBox="0 0 24 24" fill="none"><path d="M9 6l6 6-6 6" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/></svg>',
  download:'<svg viewBox="0 0 24 24" fill="none"><path d="M12 4v10m0 0l-4-4m4 4l4-4M5 19h14" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/></svg>',
  table:'<svg viewBox="0 0 24 24" fill="none"><rect x="3" y="4" width="18" height="16" rx="2" stroke="currentColor" stroke-width="1.8"/><path d="M3 9h18M9 9v11" stroke="currentColor" stroke-width="1.6"/></svg>',
  effect:'<svg viewBox="0 0 24 24" fill="none"><rect x="3" y="4" width="18" height="16" rx="2" stroke="currentColor" stroke-width="1.7"/><path d="M8 4v16M13 4v16M3 9h18M3 14h18" stroke="currentColor" stroke-width="1.2" opacity=".6"/></svg>',
  gwas:'<svg viewBox="0 0 24 24" fill="none"><path d="M3 20h18" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"/><path d="M3 9h18" stroke="currentColor" stroke-width="1.2" stroke-dasharray="2 2" opacity=".55"/><circle cx="6" cy="16" r="1.3" fill="currentColor"/><circle cx="9" cy="13" r="1.3" fill="currentColor"/><circle cx="12" cy="6" r="1.6" fill="currentColor"/><circle cx="15" cy="14" r="1.3" fill="currentColor"/><circle cx="18" cy="10" r="1.3" fill="currentColor"/></svg>',
};

/* ================= TOOL REGISTRY ================= */
const GROUPS = [
  {label:'Visualization & Search', tools:['snpversity','snpgwas','snptrait']},
  {label:'Explore & Analyze', tools:['snpimpact','snpfunction', 'snpfold']},
  {label:'Compare & Relate', tools:['snpcompare','snptree']},
  {label:'External plugin', tools:['paneffect']},
  //{label:'Impute & Predict', tools:['snpimpute','snpfold']},
  //{label:'Collection', tools:['snpgermplasm']},
];
const TOOLS = {
  snpgwas:{name:'GWAS Explorer', icon:'gwas', color:'#cf8a12', cat:'Visualization & Search',
    tag:'Scan the genome for trait-associated SNPs',
    desc:'Genome-wide Manhattan plot for NAM GWAS results — pan and zoom across all 10 chromosomes, search by SNP ID, and drag-select a peak to inspect its variants. Hand a selected region straight off to SNPVersity to keep exploring.',
    feats:[['gwas','Genome-wide scan','Canvas Manhattan plot across all NAM chromosomes.'],['search','Select a region','Drag a peak to list and filter its SNPs.'],['compare','Hand off region','Send the selected interval to SNPVersity.']]},
  snpversity:{name:'SNPVersity', icon:'dna', color:'#2563eb', cat:'Visualization & Search',
    tag:'Explore variation across lines and populations',
    desc:'Explore extensive variant datasets across maize accessions. Enter a genomic interval, choose accessions, and get a color-coded table plus a downloadable VCF — with allele states, effect annotations, and DNA/protein language-model scores.'},
  snptrait:{name:'SNPTrait', icon:'leaf', color:'#1f8a4c', cat:'Visualization & Search',
    tag:'Connect variation to phenotype and trait data',
    desc:'Connect genomic variation to phenotype and trait records from the National Germplasm collection. Search, sort, and filter accessions by trait values and metadata, then move selected sets straight into other SNPTools.',
    feats:[['search','Search & filter','Filter 20,000+ accessions by trait, phenotype, and metadata.'],['leaf','Trait records','Disease resistance, yield, composition, and evaluation data.'],['compare','Hand off sets','Send selected lines directly to SNPVersity and beyond.']]},
  snpimpact:{name:'SNPImpact', icon:'star', color:'#7c3aed', cat:'Explore & Analyze',
    tag:'Prioritize candidate variants using AI',
    desc:'Rank variants using AI-based allele scores together with functional annotations. Combine PlantCAD and ESM predictions with predicted effects to prioritize candidate causal variation at scale.',
    feats:[['star','Rank by impact','Order variants by combined AI score and consequence.'],['func','Score + effect','PlantCAD / ESM scores beside predicted protein effects.'],['compare','Shortlist alleles','Flag alleles for comparison and validation.']]},
  snpfunction:{name:'SNPFunction', icon:'func', color:'#2563eb', cat:'Explore & Analyze',
    tag:'Give variants biological and genomic context',
    desc:'Add biological context by integrating gene models, regulatory features, conservation, and nearby genomic evidence — interpreting variants within their local genomic region.',
    feats:[['func','Gene & domain view','Place variants on gene models and protein domains.'],['dna','Conservation layers','Overlay conservation and functional evidence tracks.'],['search','Local context','Read variants against nearby genomic features.']]},
  snpcompare:{name:'SNPCompare', icon:'compare', color:'#0e7490', cat:'Compare & Relate',
    tag:'Compare accessions, loci, and haplotypes',
    desc:'Support side-by-side comparison of selected accessions, loci, or haplotypes to reveal shared and distinguishing variants.',
    feats:[['compare','Side by side','Line up accessions or loci in one view.'],['grid','Shared vs unique','Highlight variants that agree or differ.'],['tree','Haplotype blocks','Compare haplotype structure across sets.']]},
  snptree:{name:'SNPTree', icon:'tree', color:'#15803d', cat:'Compare & Relate',
    tag:'Build local phylogenetic trees',
    desc:'Build local phylogenetic trees from regional variation data to examine evolutionary relationships, introgression, and haplotype structure.',
    feats:[['tree','Regional trees','Neighbor-joining trees from the active region.'],['dna','Introgression','Spot shared ancestry and introgressed blocks.'],['compare','Cluster lines','Group accessions by sequence similarity.']]},
  snpmatrix:{name:'SNPMatrix', icon:'grid', color:'#b45309', cat:'Compare & Relate',
    tag:'Pairwise identity and nearest neighbors',
    desc:'Calculate pairwise identity across all accessions and instantly rank the nearest genetic neighbors to a selected accession.',
    feats:[['grid','Identity matrix','Pairwise identity-by-state across the set.'],['search','Nearest neighbors','Rank closest lines to any accession.'],['density','Similarity heatmap','Read structure at a glance.']]},
  paneffect:{name:'PanEffect', icon:'effect', color:'#b45309', cat:'Explore & Analyze',
    tag:'Predict missense variant effects on proteins',
    desc:'Visualize the predicted effect of every amino-acid substitution across a protein — in B73 and across the maize pan-genome — using ESM protein language-model scores, with Pfam domains and secondary structure for context.',
    feats:[['star','Every substitution','ESM scores for all 20 substitutions at each residue.'],['grid','B73 + pan-genome','Reference view plus natural variation across the pan-genome.'],['func','Domain context','Pfam domains and predicted secondary structure overlaid.']]},
  snpimpute:{name:'SNPImpute', icon:'impute', color:'#0891b2', cat:'Impute & Predict',
    tag:'Impute sequence and function',
    desc:'Pan-genome–guided imputation that fills missing genotypes and predicts function across varying depths of sequencing.',
    feats:[['impute','Fill genotypes','Impute missing calls from the pan-genome.'],['dna','Depth aware','Works across light and deep sequencing.'],['func','Predict function','Carry functional predictions through imputation.']]},
  snpfold:{name:'SNPFold', icon:'fold', color:'#be185d', cat:'Impute & Predict',
    tag:'Link variants to protein structure',
    desc:'Link coding variants to protein structure, domain architecture, and AI-based folding or effect predictions.',
    feats:[['fold','Structure view','Map variants onto folded protein models.'],['func','Domain architecture','See which domains a change disrupts.'],['star','Effect prediction','AI predictions of structural impact.']]},
  snpdensity:{name:'SNPDensity', icon:'density', color:'#9333ea', cat:'Collection',
    tag:'SNP and INDEL density and burden',
    desc:'Measure SNP and INDEL density, burden, and distribution across genes and genomic regions to highlight mutational load, constraint, and diversification.',
    feats:[['density','Density tracks','Variant density across genes and regions.'],['star','Mutational load','Surface burden and constraint.'],['search','Diversification','Find rapidly diversifying regions.']]},
  snpgermplasm:{name:'SNPGermplasm', icon:'germ', color:'#16a34a', cat:'Collection',
    tag:'Genotype-driven collection management',
    desc:'Apply genotype-driven analytics to germplasm management — identifying redundancy, uniqueness, and priority materials for curation and deployment.',
    feats:[['germ','Curate collections','Manage accessions with genotype evidence.'],['grid','Redundancy','Find duplicate and unique materials.'],['star','Priority lines','Flag materials for deployment.']]},
};


/* ================= SHARED HELPERS ================= */
function rnd(a,b){return a+Math.random()*(b-a)}
function pick(a){return a[Math.floor(Math.random()*a.length)]}

/* ---- MaizeGDB gene model links — shared across tools ---- */
const MAIZEGDB_GENE_BASE = 'https://www.maizegdb.org/gene_center/gene/';
function maizegdbGeneURL(g){ return MAIZEGDB_GENE_BASE + encodeURIComponent(g); }
/* True only for a bare single gene model id — intergenic / boundary loci are
   encoded as ranges (contain "..") and are not single gene models. */
function isSingleGeneModel(g){ return !!g && g!=='—' && !String(g).includes('..') && !/\s/.test(g); }

/* ================= SHARED STATE =================
   Cross-tool query context lives here (region + selected accessions) so a
   selection made in one tool can flow into another. Tool-only view-state
   should live inside that tool\u2019s own file (see IMP in snpimpact.js). */
const S = {
  tool:'snpversity',
  dataset:'mgdb2026_hq',
  chr:'chr10', start:9788000, end:9826500, perPage:100,
  selected:new Set(),
  results:null, page:1,
  fImpact:'all', fEffect:'all', fMaf:0,
};

/* ================= TOOL REGISTRY ================= */
const SNPTools = { registry:{}, register(id,def){ this.registry[id]=def; } };

/* ================= RENDER NAV ================= */
function renderNav(){
  const nav=document.getElementById('nav');
  nav.innerHTML = GROUPS.map(g=>`
    <div class="navgroup">
      <div class="gl">${g.label}</div>
      ${g.tools.map(id=>{
        const t=TOOLS[id]; const active=id===S.tool?'active':'';
        const _bl={snpgwas:'new',snpversity:'updated',snptree:'new',snpcompare:'new',snpimpact:'new',snpfold:'new',snpfunction:'new',snpmatrix:'new',paneffect:'new'}[id]||'soon';
        const _bc={new:'#1f8a4c',demo:'#2563eb',soon:'#c0362c'}[_bl];
        const soon=`<span class="soon" style="color:${_bc};border-color:${_bc}">${_bl}</span>`;
        return `<button class="navitem ${active}" onclick="go('${id}')">
          <span class="ico">${ICONS[t.icon]}</span><span>${t.name}</span>${soon}</button>`;
      }).join('')}
    </div>`).join('');
  // hydrate send-menu icons
  document.querySelectorAll('[data-ico]').forEach(e=>{e.innerHTML=ICONS[e.dataset.ico]||'';e.style.display='inline-flex';e.style.width='16px';e.style.color='var(--blue-600)';});
}


/* ================= ROUTER ================= */
function go(id){
  S.tool=id; closeMenu(); toggleRail(false);
  document.getElementById('crumbTool').innerHTML='<b>'+TOOLS[id].name+'</b>';
  document.getElementById('sendBtn').style.display = id==='snpversity'?'inline-flex':'none';
  renderNav();
  window.scrollTo(0,0);
  const page=document.getElementById('page');
  page.className='page fade';
  const tool=SNPTools.registry[id];
  try{
    if(tool && typeof tool.render==='function') tool.render(page);
    else renderToolPage(id);
  }catch(err){
    console.error('['+id+'] render error:',err);
    page.innerHTML='<div class="card pad" style="border-color:#f3c2bd;background:#fff6f5;margin-top:14px">'
      +'<h3 style="font-family:var(--disp);color:#b42318;margin:0 0 8px">'+(TOOLS[id]?TOOLS[id].name:id)+' hit an error</h3>'
      +'<pre style="white-space:pre-wrap;background:#fff;border:1px solid #f0cfca;border-radius:8px;padding:12px;font-family:var(--mono);font-size:12px;color:#8a2a20">'+String(err&&err.stack||err)+'</pre></div>';
  }
}

/* ================= TOOL PLACEHOLDER PAGES ================= */
function renderToolPage(id){
  const t=TOOLS[id]; const p=document.getElementById('page'); p.className='page fade';
  const feats=t.feats||[['star','In development','This module is part of the SNPTools roadmap.'],['compare','Connected','It will share selections with the rest of the suite.'],['dna','Same data','Built on the unified variant database.']];
  p.innerHTML=`
    <div class="tool-hero">
      <div class="ti" style="background:${t.color}">${ICONS[t.icon]}</div>
      <div>
        <div class="cat">${t.cat}</div>
        <h1>${t.name}</h1>
        <p>${t.desc}</p>
      </div>
      <span class="soon-badge">IN DEVELOPMENT</span>
    </div>
    <div class="feat-grid">
      ${feats.map(f=>`<div class="feat"><div class="fi">${ICONS[f[0]]}</div><h4>${f[1]}</h4><p>${f[2]}</p></div>`).join('')}
    </div>
    <div class="mock-strip">
      <div class="ms-h">Part of the integrated SNPTools platform</div>
      <p style="margin:0;color:var(--muted);font-size:13px">Selections flow between tools — pick accessions in <a href="#" onclick="go('snpversity');return false">SNPVersity</a>, then send them here for ${t.tag.toLowerCase()}. Built on the same unified variant database, annotations, and B73 v5 coordinates.</p>
    </div>`;
}


/* ================= IBS COST MODEL =================
   Work is measured in (variants × accessions²). Two per-unit rates, calibrated
   2026-09 against real parsed results:
     SNPTree / SNPMatrix distance matrix : ~4.5 ns/unit
     SNPCompare (region) all-pairs       : ~1.2 ns/unit  (~3.8x cheaper)
   Used to (a) warn in each tool before a long synchronous compute and (b) grey
   out the SNPVersity send buttons when a query would essentially kill the tab.
   SNPImpact never touches genotypes, so it scales on variant count alone. */
const IBS_COST = {
  // SNPTree / SNPMatrix (≈4.5 ns/unit)
  workWarn:    2e9,      // ≈ 9 s  -> tool shows a "build anyway" notice
  workBlock:   1.5e10,   // ≈ 65 s -> SNPVersity greys the send button
  // SNPCompare region all-pairs (≈1.2 ns/unit)
  cmpWorkWarn: 7.5e9,    // ≈ 9 s
  cmpWorkBlock: 5e10,    // ≈ 60 s
  cmpMemWarn:  1.5e8,    // variants × accessions bytes ≈ 150 MB genotype matrix
  cmpMemBlock: 3e8,      //                             ≈ 300 MB
  // SNPMatrix render
  mtxAccWarn:  450,      // an A×A SVG heatmap past here is slow to draw and hard to read
  // SNPImpact (variant count only)
  impactSoft:  600000,   // ranking starts to lag
  impactBlock: 2000000,  // too many to rank in the browser
};
function ibsWork(v, a){ return v * a * a; }
function ibsSeconds(v, a){ return Math.round(4.5e-9 * v * a * a); }   // SNPTree / SNPMatrix
function cmpSeconds(v, a){ return Math.round(1.2e-9 * v * a * a); }   // SNPCompare region


/* ================= TOOLTIPS ================= */
const tt=document.getElementById('tt');
function attachTT(){
  document.querySelectorAll('[data-tt]').forEach(el=>{
    el.addEventListener('mouseenter',e=>{tt.innerHTML=el.dataset.tt;tt.classList.add('show');});
    el.addEventListener('mousemove',e=>{tt.style.left=(e.clientX+12)+'px';tt.style.top=(e.clientY+14)+'px';});
    el.addEventListener('mouseleave',()=>tt.classList.remove('show'));
  });
}


/* ================= MENUS / RAIL ================= */
function toggleMenu(){document.getElementById('sendMenu').classList.toggle('open');}
function closeMenu(){document.getElementById('sendMenu').classList.remove('open');}
document.addEventListener('click',e=>{if(!e.target.closest('.send-wrap'))closeMenu();});
function toggleRail(show){document.getElementById('rail').classList.toggle('show',show);document.getElementById('scrim').classList.toggle('show',show);}


/* ================= INIT (runs after every tool file has registered) ================= */
function init(){
  try{ renderNav(); go('snpversity'); }
  catch(err){
    console.error('SNPTools init error:',err);
    var p=document.getElementById('page');
    if(p) p.innerHTML='<pre style="white-space:pre-wrap;color:#b42318;padding:20px;font-family:monospace">'+String(err&&err.stack||err)+'</pre>';
  }
}
// tool <script>s are synchronous and sit before this point in the body,
// so by the time DOMContentLoaded fires they have all called SNPTools.register().
if(document.readyState==='complete') setTimeout(init,0);
else document.addEventListener('DOMContentLoaded',init);
