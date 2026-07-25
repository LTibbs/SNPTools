/* =====================================================================
 *  snpversity.js — Visualization & Search tool.
 *  Registers \u2018snpversity\u2019. All render fns are global so the inline
 *  onclick/onchange handlers in the markup resolve against them.
 * ===================================================================== */

/* reference data + chromosome geometry come from the data layer */
const DATASETS    = Data.datasets();
/* PROJECTS + ACCESSIONS depend on the chosen dataset (each .h5 family has
   its own real accession columns), so they are reassigned in selectDataset. */
let   PROJECTS    = Data.projectsFor(S.dataset);
let   ACCESSIONS  = Data.accessionsFor(S.dataset);
const GENE_MODELS = Data.geneModels();
const CHR_LEN     = Data.chromLengths();
const CENTRO      = Data.centromeres();

/* friendly default: preselect up to 12 founders (one run each) for this dataset */
Data.defaultSelectionFor(S.dataset).forEach(id=>S.selected.add(id));

/* default query region */
S.chr='chr10'; S.start=3750832; S.end=3755732;

/* ================= SNPVERSITY PAGE ================= */
function renderVersity(){
  injectVersityCSS();
  const inbound = applyPendingRequest();   // e.g. carriers handed over from SNPFunction
  const p=document.getElementById('page');
  p.className='page fade';
  p.innerHTML = `
    ${inbound ? inboundBanner(inbound) : ''}
    <div class="sec"><div class="bar"></div><div style="width:100%">
      <div class="n">VCF BUILDER & VIEWER · B73 v5</div>
      <h2>Build a variant view across maize accessions</h2>
      <p>Choose a dataset, set a genomic interval, and pick the accessions you want. SNPVersity returns a color-coded variant table and a downloadable VCF — alleles, effects, and DNA/protein language-model scores included.</p>
      <div style="margin-top:12px;display:flex;gap:8px;flex-wrap:wrap">
        <button class="btn" onclick="sendToCompare()">${ICONS.compare||ICONS.grid||''} Send to SNPCompare</button>
      </div>
    </div></div>

    <!-- DATASET -->
    <div class="sec"><div class="bar"></div><div><h2 style="font-size:16px">1 · Choose a dataset</h2></div></div>
    <div class="ds-grid" id="dsGrid"></div>

    <!-- REGION -->
    <div class="sec"><div class="bar"></div><div><h2 style="font-size:16px">2 · Select a genomic interval</h2></div></div>
    <div class="card region-card" id="regionCard"></div>

    <!-- ACCESSIONS -->
    <div class="sec"><div class="bar"></div><div><h2 style="font-size:16px">3 · Choose accessions</h2>
      <p>Search and toggle accessions, grab a random sample, or upload a list. Selected lines collect on the right.</p></div></div>
    <div class="card acc-card" id="accCard"></div>

    <!-- RUN -->
    <div class="runbar" id="runbar"></div>

    <!-- RESULTS -->
    <div id="resultsAnchor" style="margin-top:30px"></div>
  `;
  renderDatasets(); renderRegion(); renderAccPicker(); renderRunbar();
  if(inbound){
    // prefill the gene box so "Load" re-fetches the same coordinates
    const gi=document.getElementById('geneInput');
    if(gi && inbound.gene){ gi.value=inbound.gene; }
    const st=document.getElementById('geneStatus');
    if(st && inbound.gene){ st.className='status ok'; st.textContent=`Region set from ${inbound.gene} (${S.chr}:${S.start.toLocaleString()}–${S.end.toLocaleString()})`; }
    if(inbound.missing && inbound.missing.length){
      uplSay('warn', `${inbound.missing.length} carrier${inbound.missing.length>1?'s are':' is'} not present in this dataset's accession list and could not be selected.`);
      renderUplReport({matched:inbound.selected, considered:inbound.selected+inbound.missing.length,
        unmatched:inbound.missing.map((t,i)=>({line:i+1, text:t})), dupes:[], fanout:[], source:'SNPFunction', skippedHeader:null, applied:true});
    }
    document.getElementById('runbar').scrollIntoView({behavior:'smooth',block:'center'});
    return;
  }
  // returning from another tool: restore the last results instead of forcing a rebuild
  if(S.results && S.results.rows && S.results.rows.length) renderResults();
}

/* =====================================================================
 *  INBOUND HANDOFF  —  another tool (SNPFunction) hands us a gene model
 *  and a set of accessions to preselect.
 *    payload = {gene, chr, start, end, dataset, accessions:[…], from, note}
 * ===================================================================== */
window.versityRequest = function(payload){
  if(!payload) return;
  S.pendingVersity = payload;
  S.results = null;                  // the old result no longer matches the new query
  if(typeof go==='function') go('snpversity');
};

/* consume S.pendingVersity: switch dataset, set the region, resolve + select
   the accessions. Returns a summary for the banner, or null. */
function applyPendingRequest(){
  const req = S.pendingVersity; if(!req) return null;
  S.pendingVersity = null;

  // 1 · dataset (rebinds this dataset's real accession catalog)
  if(req.dataset && req.dataset!==S.dataset && DATASETS.some(d=>d.id===req.dataset)){
    S.dataset=req.dataset;
    PROJECTS=Data.projectsFor(S.dataset);
    ACCESSIONS=Data.accessionsFor(S.dataset);
  }

  // 2 · region
  if(req.chr){ S.chr=String(req.chr).startsWith('chr')?req.chr:'chr'+req.chr; }
  const flank=+req.flank||0;
  if(req.start!=null && req.end!=null){
    const lo=Math.min(+req.start,+req.end), hi=Math.max(+req.start,+req.end);
    S.start=Math.max(0,lo-flank); S.end=hi+flank;
  }

  // 3 · accessions (resolved the same way as an uploaded list, so founder
  //     names / run IDs / composite IDs all work)
  const wanted=[...new Set((req.accessions||[]).filter(Boolean).map(String))];
  const idx=buildAccIndex(); const missing=[];
  S.selected.clear();
  wanted.forEach(w=>{
    const hit=idxHit(idx,w) || (w.match(RUN_RE)?idxHit(idx,w.match(RUN_RE)[1]):null);
    if(hit) hit.forEach(id=>S.selected.add(id)); else missing.push(w);
  });
  accFilter='';
  S.page=1;
  return {gene:req.gene||'', from:req.from||'SNPFunction', note:req.note||'',
          requested:wanted.length, selected:S.selected.size, missing};
}
function inboundBanner(i){
  return `<div class="from-fn">
    <b>From ${escAttr(i.from)}</b>
    <span>${i.gene?`<span class="mono">${escAttr(i.gene)}</span> — `:''}${escAttr(i.note||'carrier accessions')}:
      <b>${i.selected}</b> of ${i.requested} accession${i.requested===1?'':'s'} preselected${i.missing.length?`, ${i.missing.length} not in this dataset`:''}.</span>
    <button class="btn" style="margin-left:auto" onclick="runQuery()">${ICONS.dna||''} Build VCF &amp; view</button>
  </div>`;
}

function renderDatasets(){
  document.getElementById('dsGrid').innerHTML = DATASETS.map(d=>`
    <div class="ds ${d.id===S.dataset?'sel':''}" onclick="selectDataset('${d.id}')">
      <div class="dot"></div>
      <div class="t">${d.name}</div>
      <div class="ref">${d.sub} · aligned to ${d.ref}</div>
      <div class="stats">
        <div class="stat"><div class="v">${d.acc}</div><div class="k">accessions</div></div>
        <div class="stat"><div class="v">${d.sites}</div><div class="k">variant sites</div></div>
      </div>
      <div class="badges">
        <span class="chiplet ${d.het?'on':''}">${d.het?'✓':'×'} heterozygous</span>
        <span class="chiplet ${d.indel?'on':''}">${d.indel?'✓':'×'} INDELs</span>
        <span class="chiplet ${d.impute?'on':''}">${d.impute?'✓':'×'} imputed</span>
      </div>
    </div>`).join('');
}
function selectDataset(id){
  S.dataset=id;
  // switch to this dataset's real accession catalog
  PROJECTS   = Data.projectsFor(id);
  ACCESSIONS = Data.accessionsFor(id);
  S.selected.clear();
  Data.defaultSelectionFor(id).forEach(x=>S.selected.add(x));
  accFilter='';
  renderDatasets(); renderAccPicker(); renderRunbar();
  const m=document.getElementById('mAcc'); if(m)m.textContent=S.selected.size;
}

/* ---- region card + genome ribbon ---- */
function renderRegion(){
  const len=CHR_LEN[S.chr]||3e8;
  document.getElementById('regionCard').innerHTML=`
    <div class="ribbon-pane">
      <div class="ribbon-head">
        <span class="chr" id="ribChr">Chromosome ${S.chr.replace('chr','')}</span>
        <span class="coord" id="ribCoord"></span>
      </div>
      <div class="ribbon">
        <div class="chrom-track"></div>
        <div class="centro" id="centro"></div>
        <div class="window" id="window"></div>
        <div class="ticks"><span>0</span><span>${(len/1e6).toFixed(0)} Mb</span></div>
      </div>
      <div class="region-meta">
        <div class="m"><div class="v" id="mSpan"></div><div class="k">interval span</div></div>
        <div class="m"><div class="v" id="mMode"></div><div class="k">view mode</div></div>
        <div class="m"><div class="v" id="mAcc"></div><div class="k">accessions</div></div>
      </div>
    </div>
    <div class="form-pane">
      <div class="field">
        <label>Chromosome</label>
        <select id="chrInput" onchange="onRegion()">
          ${Object.keys(CHR_LEN).map(c=>`<option value="${c}" ${c===S.chr?'selected':''}>Chromosome ${c.replace('chr','')}</option>`).join('')}
        </select>
      </div>
      <div class="field">
        <label>Interval (bp)</label>
        <div class="row">
          <input type="number" id="startInput" value="${S.start}" min="0" oninput="onRegion()" aria-label="start">
          <span style="align-self:center;color:var(--faint)">–</span>
          <input type="number" id="endInput" value="${S.end}" min="0" oninput="onRegion()" aria-label="end">
        </div>
      </div>
      <div class="gene-row">
        <div class="field" style="margin:0">
          <label>Gene model · auto-fill coordinates</label>
          <input type="text" id="geneInput" placeholder="Zm00001eb…" class="mono-in">
        </div>
        <button class="btn" onclick="loadGene()" style="margin-bottom:1px">Load</button>
      </div>
      <div style="display:flex;gap:14px;align-items:center">
        <button class="link-btn" onclick="exampleGene()">Use an example gene</button>
        <div class="field" style="margin:0;flex:1">
          <label>± flank (bp)</label>
          <input type="number" id="flankInput" value="0" min="0" class="mono-in">
        </div>
        <div class="field" style="margin:0;width:120px">
          <label>Loci / page</label>
          <select id="perPage" onchange="S.perPage=+this.value">
            <option>50</option><option selected>100</option><option>250</option><option>500</option>
          </select>
        </div>
      </div>
      <div class="status" id="geneStatus"></div>
    </div>`;
  drawRibbon();
}
function onRegion(){
  S.chr=document.getElementById('chrInput').value;
  S.start=+document.getElementById('startInput').value||0;
  S.end=+document.getElementById('endInput').value||0;
  drawRibbon(); renderRunbar();
  // update chr label
  document.getElementById('ribChr').textContent='Chromosome '+S.chr.replace('chr','');
}
function drawRibbon(){
  const len=CHR_LEN[S.chr]||3e8;
  const a=Math.max(0,Math.min(S.start,len)), b=Math.max(0,Math.min(S.end,len));
  const lo=Math.min(a,b), hi=Math.max(a,b);
  const w=document.getElementById('window');
  const left=(lo/len)*100, width=Math.max(((hi-lo)/len)*100,0.4);
  if(w){w.style.left=left+'%';w.style.width=width+'%';}
  const ce=document.getElementById('centro'); if(ce)ce.style.left=((CENTRO[S.chr]||.4)*100)+'%';
  const span=hi-lo;
  const fmt=n=>n>=1e6?(n/1e6).toFixed(2)+' Mb':n>=1e3?(n/1e3).toFixed(1)+' kb':n+' bp';
  document.getElementById('ribCoord').textContent=`${S.chr}:${lo.toLocaleString()}–${hi.toLocaleString()}`;
  document.getElementById('mSpan').textContent=fmt(span);
  document.getElementById('mMode').textContent= span>1e6?'VCF download':'Interactive table';
  document.getElementById('mAcc').textContent=S.selected.size;
}
async function loadGene(){
  const id=document.getElementById('geneInput').value.trim();
  const st=document.getElementById('geneStatus');
  if(!id){st.className='status err';st.textContent='Enter a gene model ID first.';return;}
  st.className='status'; st.textContent='Looking up '+id+'…';
  let g;
  try{
    g=await Data.lookupGene(id);
  }catch(err){
    st.className='status err'; st.textContent='Lookup failed: '+((err&&err.message)?err.message:err);
    console.error('Gene lookup failed:', err);
    return;
  }
  if(!g){st.className='status err';st.textContent='Gene model “'+id+'” not found.';return;}
  const flank=+document.getElementById('flankInput').value||0;
  S.chr=g.chr; S.start=Math.max(0,g.start-flank); S.end=g.end+flank;
  const chrSel=document.getElementById('chrInput');
  if(chrSel){
    if(![...chrSel.options].some(o=>o.value===g.chr)){
      chrSel.add(new Option('Chromosome '+g.chr.replace('chr',''), g.chr));
    }
    chrSel.value=g.chr;
  }
  document.getElementById('startInput').value=S.start;
  document.getElementById('endInput').value=S.end;
  document.getElementById('ribChr').textContent='Chromosome '+g.chr.replace('chr','');
  st.className='status ok'; st.textContent='Loaded '+id+' ('+g.chr+':'+g.start.toLocaleString()+'–'+g.end.toLocaleString()+')';
  drawRibbon(); renderRunbar();
}
function exampleGene(){
  const ids=Data.exampleGenes(); const id=pick(ids);
  document.getElementById('geneInput').value=id; loadGene();
}

/* ---- accession picker ---- */
let accFilter='';
const openProjects=new Set();   // project ids currently expanded (persists across re-renders)
function escAttr(s){return String(s==null?'':s).replace(/&/g,'&amp;').replace(/"/g,'&quot;').replace(/</g,'&lt;').replace(/>/g,'&gt;');}
/* REF/ALT sequences can be long indels; truncate the visible text to keep the
   table columns from stretching, while the full sequence stays available on
   hover via the data-tt tooltip. */
function alleleDisp(seq){
  const s=String(seq==null?'':seq);
  return s.length>10 ? s.slice(0,10)+'…' : s;
}
function alleleTT(seq,label){
  const s=String(seq==null?'':seq);
  return s.length>10 ? `${label}: ${s}` : label;
}
function renderAccPicker(){
  document.getElementById('accCard').innerHTML=`
    <div class="acc-top">
      <div class="acc-search">
        ${ICONS.search}
        <input type="text" id="accSearch" placeholder="Search accession, run, or founder…" oninput="accFilter=this.value.toLowerCase();renderAccList()">
      </div>
      <div class="quick">
        <button class="qbtn" onclick="randomSel(.02)">Random 2%</button>
        <button class="qbtn" onclick="randomSel(.05)">5%</button>
        <button class="qbtn" onclick="randomSel(.10)">10%</button>
        <button class="qbtn" onclick="randomSel(.25)">25%</button>
        <button class="qbtn" onclick="onePerFounder()">All unique accessions</button>
        ${Data.namFoundersFor(S.dataset).length?`<button class="qbtn" onclick="oneNAMEach()">One per NAM founder</button>`:''}
        <button class="qbtn solid" onclick="allSel(true)">Select all</button>
        <button class="qbtn" onclick="allSel(false)">Clear</button>
      </div>
    </div>
    <div class="acc-body">
      <div class="acc-list" id="accList"></div>
      <div class="acc-side">
        <div class="sh"><span class="ttl">Selected</span><span class="ct" id="selCount"></span></div>
        <div class="sel-chips" id="selChips"></div>
        <div class="upl">
          <label>Or upload a list (one accession per line)</label>
          <div class="file-row">
            <input type="file" id="fileUpload" accept=".txt,.tsv,.csv,.list,text/plain" onchange="onAccFilePicked(this)">
          </div>
          <div class="upl-name" id="uplName">No file chosen — you can also paste a list below.</div>
          <textarea id="accPaste" class="upl-paste" spellcheck="false"
            placeholder="ACC.8750_SRR12460455&#10;ACC.8782_SRR12460453&#10;SRR12460421&#10;…"
            oninput="refreshUplButtons()"></textarea>
          <div class="upl-actions">
            <button class="btn primary" id="uplLoadBtn" onclick="loadAccList()" disabled>Load accessions</button>
            <button class="btn" onclick="clearAccUpload()">Clear</button>
          </div>
          <label class="upl-opt"><input type="checkbox" id="uplReplace" checked> Replace current selection</label>
          <div class="upl-status" id="uplStatus"></div>
          <div class="upl-report" id="uplReport"></div>
        </div>
      </div>
    </div>`;
  // fresh picker (first open or dataset switch): start with the first section open
  openProjects.clear();
  if(PROJECTS[0] && PROJECTS[0].count<=250) openProjects.add(PROJECTS[0].id);
  renderAccList(); renderSelected();
}
function accMatch(a){
  if(!accFilter) return true;
  return a.id.toLowerCase().includes(accFilter)
      || (a.founder||'').toLowerCase().includes(accFilter)
      || (a.run||'').toLowerCase().includes(accFilter)
      || (a.label||'').toLowerCase().includes(accFilter);
}
function projAccIds(p){
  const ids=[]; p.groups.forEach(g=>g.accessions.forEach(a=>ids.push(a.id))); return ids;
}
function findProject(pid){ return PROJECTS.find(p=>p.id===pid); }
function accChipHTML(a){
  return `<span class="acc-chip ${S.selected.has(a.id)?'on':''}" onclick="toggleAcc('${a.id}')">
    <span class="cb">${ICONS.check}</span>
    <span>${a.run}</span><span class="founder">${a.founder}${a.reps>1?' · r'+a.rep:''}</span>
  </span>`;
}
function renderAccList(){
  const html=PROJECTS.map((p,i)=>{
    const groups=p.groups.map(g=>({name:g.name, items:g.accessions.filter(accMatch)}))
                         .filter(g=>g.items.length);
    const shown=groups.reduce((n,g)=>n+g.items.length,0);
    if(accFilter && !shown) return '';
    const selN=projAccIds(p).reduce((n,id)=>n+(S.selected.has(id)?1:0),0);
    const open=((accFilter && shown) || openProjects.has(p.id))?'open':'';
    const meta=[
      (p.bioprojects&&p.bioprojects.length)?p.bioprojects.join(', '):'',
      p.ncbiUrl?`<a href="${p.ncbiUrl}" target="_blank" rel="noopener" onclick="event.stopPropagation()">NCBI</a>`:'',
      p.referenceUrl?`<a href="${p.referenceUrl}" target="_blank" rel="noopener" onclick="event.stopPropagation()">reference</a>`:''
    ].filter(Boolean).join(' · ');
    return `<div class="proj ${open}" id="proj_${p.id}">
      <div class="proj-h" onclick="toggleProj('${p.id}')">
        <span class="caret">${ICONS.caret}</span>
        <span class="swatch" style="background:${p.color}"></span>
        <div style="min-width:0">
          <div class="pt">${p.title}</div>
          <div class="pm">${p.count} accessions${meta?' · '+meta:''}</div>
        </div>
        <span class="pc">${selN}/${p.count}</span>
      </div>
      <div class="proj-items">
        <div style="flex:0 0 100%;display:flex;gap:6px;flex-wrap:wrap;margin:2px 0 8px">
          <button class="qbtn" onclick="event.stopPropagation();selProject('${p.id}',0)">None</button>
          <button class="qbtn" onclick="event.stopPropagation();selProject('${p.id}',.25)">Random 25%</button>
          <button class="qbtn" onclick="event.stopPropagation();selProject('${p.id}',.5)">Random 50%</button>
          <button class="qbtn solid" onclick="event.stopPropagation();selProject('${p.id}',1)">Select all</button>
        </div>
        ${groups.map(g=>`
          ${g.name?`<div style="flex:0 0 100%;font-weight:600;font-size:11.5px;color:var(--muted);margin:8px 0 2px">Group: ${g.name}</div>`:''}
          <div style="flex:0 0 100%;display:flex;flex-wrap:wrap;gap:6px">${g.items.map(accChipHTML).join('')}</div>
        `).join('')}
      </div></div>`;
  }).join('');
  document.getElementById('accList').innerHTML=html||'<div class="empty" style="padding:30px;text-align:center;color:var(--faint)">No accessions match that search.</div>';
}
function toggleProj(pid){
  const el=document.getElementById('proj_'+pid); if(!el)return;
  el.classList.toggle('open');
  if(el.classList.contains('open')) openProjects.add(pid); else openProjects.delete(pid);
}
function selProject(pid,frac){
  const p=findProject(pid); if(!p)return;
  const ids=projAccIds(p);
  if(frac<=0){ ids.forEach(id=>S.selected.delete(id)); }
  else if(frac>=1){ ids.forEach(id=>S.selected.add(id)); }
  else {
    ids.forEach(id=>S.selected.delete(id));
    const shuffled=ids.slice().sort(()=>Math.random()-.5);
    const n=Math.ceil(ids.length*frac);
    for(let i=0;i<n;i++) S.selected.add(shuffled[i]);
  }
  renderAccList(); renderSelected(); renderRunbar();
}
function oneNAMEach(){
  // add one accession per tagged NAM founder (2026 dataset)
  const seen=new Set();
  ACCESSIONS.forEach(a=>{ if(a.namFounder && !seen.has(a.namFounder)){ seen.add(a.namFounder); S.selected.add(a.id); } });
  renderAccList(); renderSelected(); renderRunbar();
}
function renderSelected(){
  const arr=[...S.selected];
  const c=document.getElementById('selCount'); if(c)c.textContent=arr.length+' / '+ACCESSIONS.length;
  const box=document.getElementById('selChips'); if(!box)return;
  if(!arr.length){box.innerHTML='<div class="empty">Nothing selected yet.<br>Toggle accessions or grab a random sample.</div>';}
  else{
    box.innerHTML=arr.map(id=>{const a=ACCESSIONS.find(x=>x.id===id);
      return `<span class="sel-chip"><span class="dotc" style="display:inline-block;width:7px;height:7px;border-radius:2px;background:${a?a.projColor:'#999'}"></span>${a?a.run:id}<button onclick="toggleAcc('${id}')" aria-label="remove">×</button></span>`;
    }).join('');
  }
  const m=document.getElementById('mAcc'); if(m)m.textContent=arr.length;
}
/* =====================================================================
 *  ACCESSION LIST UPLOAD  (file or paste)
 *  Reads a plain list, resolves each entry against the current dataset's
 *  accession catalog, and reports exactly what matched / what did not.
 * ===================================================================== */
let uplFileText = null;      // text of the last successfully read file
let uplFileName = '';

const MAX_UPL_BYTES = 5 * 1024 * 1024;
const RUN_RE = /\b([SEDC]RR\d{4,})\b/i;                 // SRR / ERR / DRR / CRR run accessions
const RUN_SUFFIX_RE = /[_\s,\t-]+[SEDC]RR\d{4,}\s*$/i;  // trailing "_SRR12460455"

function normKey(s){ return String(s==null?'':s).replace(/\uFEFF/g,'').trim().replace(/^["']+|["']+$/g,'').toLowerCase(); }
function looseKey(s){ return normKey(s).replace(/[^a-z0-9]+/g,''); }

/* key → Set(accession id). Every accession is indexed under its id, run,
   founder, label and the founder_run composite (both orders), in exact and
   punctuation-insensitive form. Rebuilt on demand so it follows the dataset. */
function buildAccIndex(){
  const idx=new Map();
  const add=(k,id)=>{ if(!k) return; if(!idx.has(k)) idx.set(k,new Set()); idx.get(k).add(id); };
  ACCESSIONS.forEach(a=>{
    const keys=[a.id, a.run, a.founder, a.label];
    if(a.founder&&a.run){ keys.push(a.founder+'_'+a.run, a.run+'_'+a.founder); }
    keys.filter(Boolean).forEach(v=>{ add(normKey(v),a.id); add(looseKey(v),a.id); });
  });
  return idx;
}
function idxHit(idx,key){ return idx.get(normKey(key)) || idx.get(looseKey(key)) || null; }

/* resolve one line to accession id(s); returns {ok, ids, via} or {ok:false, text} */
function resolveAccEntry(raw, idx){
  const line=String(raw).replace(/\uFEFF/g,'').trim();
  if(!line || /^[#!]/.test(line)) return null;                 // blank or comment → skipped
  // whole line first, then each delimited field (handles csv/tsv exports)
  const toks=[line].concat(line.split(/[,\t;|]+/)).map(t=>t.trim().replace(/^["']+|["']+$/g,'')).filter(Boolean);
  for(const t of toks){ const hit=idxHit(idx,t); if(hit) return {ok:true, ids:[...hit], via:t, text:line}; }
  // fall back to the embedded run accession …
  const m=line.match(RUN_RE);
  if(m){ const hit=idxHit(idx,m[1]); if(hit) return {ok:true, ids:[...hit], via:m[1], text:line}; }
  // … then to the founder part with the run suffix stripped
  const f=line.replace(RUN_SUFFIX_RE,'').trim();
  if(f && f!==line){ const hit=idxHit(idx,f); if(hit) return {ok:true, ids:[...hit], via:f, text:line}; }
  return {ok:false, text:line};
}

function uplSay(cls,msg){ const el=document.getElementById('uplStatus'); if(el){ el.className='upl-status '+(cls||''); el.innerHTML=msg||''; } }
function refreshUplButtons(){
  const btn=document.getElementById('uplLoadBtn'); if(!btn) return;
  const pasted=(document.getElementById('accPaste')||{}).value||'';
  btn.disabled = !(uplFileText && uplFileText.trim()) && !pasted.trim();
}
function clearAccUpload(){
  uplFileText=null; uplFileName='';
  const f=document.getElementById('fileUpload'); if(f) f.value='';
  const t=document.getElementById('accPaste');  if(t) t.value='';
  const n=document.getElementById('uplName');   if(n) n.textContent='No file chosen — you can also paste a list below.';
  const r=document.getElementById('uplReport'); if(r) r.innerHTML='';
  uplSay('',''); refreshUplButtons();
}

/* file chosen → validate + read into memory (nothing is selected until "Load") */
function onAccFilePicked(input){
  const file=input.files && input.files[0];
  const nameEl=document.getElementById('uplName');
  const rep=document.getElementById('uplReport'); if(rep) rep.innerHTML='';
  uplFileText=null; uplFileName='';
  if(!file){ if(nameEl) nameEl.textContent='No file chosen — you can also paste a list below.'; uplSay('',''); refreshUplButtons(); return; }
  if(file.size>MAX_UPL_BYTES){
    if(nameEl) nameEl.textContent=file.name;
    uplSay('err',`That file is ${(file.size/1048576).toFixed(1)} MB. Please upload a plain list under 5 MB.`);
    refreshUplButtons(); return;
  }
  if(file.size===0){
    if(nameEl) nameEl.textContent=file.name;
    uplSay('err','That file is empty.'); refreshUplButtons(); return;
  }
  if(/\.(xlsx|xls|pdf|docx?|zip|gz|bam|h5|vcf)$/i.test(file.name)){
    if(nameEl) nameEl.textContent=file.name;
    uplSay('err','Unsupported file type. Use a plain text, .csv or .tsv list with one accession per line.');
    refreshUplButtons(); return;
  }
  const fr=new FileReader();
  fr.onerror=()=>{ uplSay('err','Could not read that file.'); refreshUplButtons(); };
  fr.onload=()=>{
    const txt=String(fr.result||'');
    if(/\u0000/.test(txt)){ uplSay('err','That looks like a binary file, not a text list.'); refreshUplButtons(); return; }
    uplFileText=txt; uplFileName=file.name;
    if(nameEl) nameEl.textContent=`${file.name} · ${(file.size/1024).toFixed(1)} KB · ${txt.split(/\r\n|\r|\n/).filter(l=>l.trim()).length} non-empty lines`;
    uplSay('ok','File read. Press <b>Load accessions</b> to match it against this dataset.');
    refreshUplButtons();
  };
  fr.readAsText(file);
}

/* main entry: parse whatever is available (file wins, else the textarea) */
function loadAccList(){
  const pasted=((document.getElementById('accPaste')||{}).value||'');
  const text = (uplFileText && uplFileText.trim()) ? uplFileText : pasted;
  const src  = (uplFileText && uplFileText.trim()) ? (uplFileName||'uploaded file') : 'pasted list';
  if(!text || !text.trim()){ uplSay('err','Nothing to load — choose a file or paste a list first.'); return; }
  applyAccList(text, src);
}

function applyAccList(text, source){
  const idx=buildAccIndex();
  const lines=String(text).split(/\r\n|\r|\n/);
  const matched=new Map();      // accession id → the entry that matched it
  const unmatched=[];           // {line, text}
  const dupes=[];               // entries that resolved to something already matched
  const fanout=[];              // entries that resolved to >1 accession (e.g. a founder with reps)
  let considered=0, skippedHeader=null;

  lines.forEach((raw,i)=>{
    const r=resolveAccEntry(raw, idx);
    if(r===null) return;                                   // blank / comment
    // tolerate a single header row at the top
    if(!r.ok && considered===0 && skippedHeader===null &&
       /^(accession|accessions|id|sample|sample_?id|name|run|taxa|line)\b/i.test(r.text)){
      skippedHeader=r.text; return;
    }
    considered++;
    if(!r.ok){ unmatched.push({line:i+1, text:r.text}); return; }
    if(r.ids.length>1) fanout.push({text:r.text, n:r.ids.length});
    r.ids.forEach(id=>{ if(matched.has(id)) dupes.push(r.text); else matched.set(id, r.text); });
  });

  if(!considered){ uplSay('err','No usable entries found — the list looks empty or contains only comments.'); return; }
  if(!matched.size){
    uplSay('err',`None of the ${considered} entries in <b>${escAttr(source)}</b> matched an accession in this dataset. `+
                 `Check that you picked the right dataset above, or that the file is one accession per line.`);
    renderUplReport({matched:0, considered, unmatched, dupes, fanout, source, skippedHeader, applied:false});
    return;
  }

  const replace=(document.getElementById('uplReplace')||{}).checked;
  if(replace) S.selected.clear();
  matched.forEach((_,id)=>S.selected.add(id));

  const cls = unmatched.length ? 'warn' : 'ok';
  uplSay(cls, `Loaded <b>${matched.size}</b> of ${considered} entries from <b>${escAttr(source)}</b>`+
              (unmatched.length?` · <b>${unmatched.length}</b> not recognized`:'')+
              ` · selection is now <b>${S.selected.size}</b> accessions.`);
  renderUplReport({matched:matched.size, considered, unmatched, dupes, fanout, source, skippedHeader, applied:true});

  renderAccList(); renderSelected(); renderRunbar();
}

function renderUplReport(r){
  const box=document.getElementById('uplReport'); if(!box) return;
  const bits=[];
  if(r.skippedHeader) bits.push(`<div class="upl-note">Skipped header row: <span class="mono">${escAttr(r.skippedHeader)}</span></div>`);
  if(r.dupes.length)  bits.push(`<div class="upl-note">${r.dupes.length} duplicate entr${r.dupes.length>1?'ies were':'y was'} ignored.</div>`);
  if(r.fanout.length) bits.push(`<div class="upl-note">${r.fanout.length} entr${r.fanout.length>1?'ies':'y'} matched more than one run (all replicates were selected).</div>`);
  if(r.unmatched.length){
    const show=r.unmatched.slice(0,25);
    bits.push(`<details class="upl-bad" open><summary>${r.unmatched.length} entr${r.unmatched.length>1?'ies':'y'} not found in this dataset</summary>
      <ul>${show.map(u=>`<li><span class="ln">line ${u.line}</span> <span class="mono">${escAttr(u.text)}</span></li>`).join('')}</ul>
      ${r.unmatched.length>show.length?`<div class="upl-note">…and ${r.unmatched.length-show.length} more.</div>`:''}
      <div class="upl-note">Accepted forms: SNPVersity ID, run accession (SRR/ERR/DRR/CRR…), founder name, or <span class="mono">FOUNDER_RUN</span>.</div>
    </details>`);
  }
  box.innerHTML=bits.join('');
}

/* minimal styling for the uploader (kept local so it can't clash with the suite CSS) */
function injectVersityCSS(){
  if(document.getElementById('snpversity-upl-css')) return;
  const s=document.createElement('style'); s.id='snpversity-upl-css';
  s.textContent=`
    .upl-name{font-size:11.5px;color:var(--muted);margin:6px 0 6px;word-break:break-all}
    .upl-paste{width:100%;box-sizing:border-box;min-height:64px;resize:vertical;border:1px solid var(--line);
      border-radius:8px;padding:7px 9px;font-family:var(--mono);font-size:11.5px;color:var(--ink);background:#fff}
    .upl-actions{display:flex;gap:6px;flex-wrap:wrap;margin-top:7px}
    .upl-opt{display:flex;align-items:center;gap:6px;font-size:11.5px;color:var(--muted);margin-top:8px;cursor:pointer}
    .upl-opt input{margin:0}
    .upl-status{font-size:11.5px;margin-top:8px;line-height:1.45}
    .upl-status.ok{color:#176c3a} .upl-status.err{color:#c0362c} .upl-status.warn{color:#8a6d1e}
    .upl-report{margin-top:8px}
    .upl-note{font-size:11px;color:var(--muted);margin-top:5px}
    .upl-bad{margin-top:7px;border:1px solid #f0c4bd;background:#fdf6f5;border-radius:8px;padding:7px 9px}
    .upl-bad summary{cursor:pointer;font-size:11.5px;font-weight:600;color:#8f281c}
    .upl-bad ul{margin:7px 0 0;padding-left:16px;max-height:190px;overflow:auto}
    .upl-bad li{font-size:11px;margin-bottom:3px}
    .upl-bad .ln{color:var(--faint);margin-right:5px}
    .upl-bad .mono,.upl-note .mono,.upl-name .mono{font-family:var(--mono)}
    .from-fn{background:#eef4ff;border:1px solid #cfe0ff;color:#274b8f;border-radius:9px;padding:9px 12px;
      font-size:12.5px;margin-bottom:12px;display:flex;gap:9px;align-items:center;flex-wrap:wrap}
    .from-fn .mono{font-family:var(--mono)}
    /* per-row jump links in the Effect column (matches .pe-jump) */
    .effect-cell .fold-jump{font-size:11px;white-space:nowrap;margin-left:2px;
      color:#176c3a;text-decoration:none;border-bottom:1px dotted #9ecdb1}
    .effect-cell .fold-jump:hover{color:#0f4f2a;border-bottom-style:solid}`;
  document.head.appendChild(s);
}

function toggleAcc(id){S.selected.has(id)?S.selected.delete(id):S.selected.add(id);renderAccList();renderSelected();renderRunbar();}
function allSel(on){ACCESSIONS.forEach(a=>on?S.selected.add(a.id):S.selected.delete(a.id));renderAccList();renderSelected();renderRunbar();}
function randomSel(p){S.selected.clear();const idx=[...ACCESSIONS.keys()].sort(()=>Math.random()-.5);const n=Math.ceil(ACCESSIONS.length*p);for(let i=0;i<n;i++)S.selected.add(ACCESSIONS[idx[i]].id);renderAccList();renderSelected();renderRunbar();}
function onePerFounder(){S.selected.clear();const seen=new Set();ACCESSIONS.forEach(a=>{if(!seen.has(a.founder)){seen.add(a.founder);S.selected.add(a.id);}});renderAccList();renderSelected();renderRunbar();}

/* ---- run bar ---- */
function renderRunbar(){
  const rb=document.getElementById('runbar'); if(!rb)return;
  const d=DATASETS.find(x=>x.id===S.dataset);
  const lo=Math.min(S.start,S.end),hi=Math.max(S.start,S.end);
  const span=hi-lo;
  const ready=S.selected.size>0 && span>0;
  rb.innerHTML=`
    <div class="summ">
      <b>${d.name} ${d.sub}</b> · <span class="mono">${S.chr}:${lo.toLocaleString()}–${hi.toLocaleString()}</span><br>
      <span style="color:#9fb4d6">${S.selected.size} accessions · ${(span/1e3).toFixed(1)} kb interval · ${span>1e6?'returns a VCF file':'returns an interactive table'}</span>
    </div>
    <div class="spacer"></div>
    <button class="btn primary" ${ready?'':'disabled'} onclick="runQuery()">
      ${ICONS.dna} Build VCF & view
    </button>`;
}

/* ---- run query (LIVE: processForm.php -> h5_to_vcf.py -> VCF) ---- */
function noticeCard(title, body){
  return `<div class="card pad fade" style="text-align:center">
    <h3 style="font-family:var(--disp);margin:0 0 8px">${title}</h3>
    <p style="color:var(--muted);margin:0 auto;max-width:560px">${body}</p></div>`;
}
async function runQuery(){
  const anchor=document.getElementById('resultsAnchor');
  anchor.innerHTML=`<div class="card" style="overflow:hidden"><div class="loading"><div class="spinner"></div><div>Reading HDF5 store and assembling VCF…</div></div></div>`;
  anchor.scrollIntoView({behavior:'smooth',block:'start'});
  const lo=Math.min(S.start,S.end), hi=Math.max(S.start,S.end);
  try{
    S.results = await Data.queryVariants(S.dataset, S.chr, lo, hi, [...S.selected]);
    if(S.results.wide){
      anchor.innerHTML=noticeCard('Interval too wide for the table view',
        `This interval spans ${((hi-lo)/1e6).toFixed(2)} Mb. The table view is available for regions up to 1 Mb. `+
        `<br><br><button class="btn" onclick="downloadVCF()">${ICONS.download} Download the VCF file</button>`);
      return;
    }
    if(S.results.empty || !S.results.rows.length){
      anchor.innerHTML=noticeCard('No variants in this range',
        'The query returned no variant sites for the selected accessions. Try a wider interval or a different region.');
      return;
    }
    S.page=1; renderResults();
  }catch(err){
    anchor.innerHTML=noticeCard('The query failed',
      `<span style="color:var(--muted)">${(err&&err.message?err.message:err)}</span>`);
    console.error('SNPVersity query failed:', err);
  }
}
/* download the VCF produced by the most recent query */
function downloadVCF(){
  const url=S.results&&S.results.vcfUrl; if(!url)return;
  const a=document.createElement('a');
  a.href=url; a.download=url.split('/').pop();
  document.body.appendChild(a); a.click(); a.remove();
}
function gColor(s){
  if(s===null||s===undefined) return '#e7ebf1';
  const min=-15,max=10, c=Math.max(min,Math.min(s,max));
  if(c<=1&&c>=-1) return '#aab4be';
  const r=(c-min)/(max-min);
  return `rgb(${Math.round(255*(1-r))},${Math.round(255*r)},0)`;
}

function renderResults(){
  const a=document.getElementById('resultsAnchor');
  const lo=Math.min(S.start,S.end), hi=Math.max(S.start,S.end);
  a.innerHTML=`
    <div class="sec"><div class="bar"></div>
      <div style="display:flex;align-items:center;gap:14px;flex-wrap:wrap;width:100%">
        <h2 style="font-size:16px;margin:0">Variant view · <span class="c-mono" style="color:var(--blue-600)">${S.chr}:${lo.toLocaleString()}–${hi.toLocaleString()}</span></h2>
        <div style="margin-left:auto;display:flex;gap:8px;flex-wrap:wrap">
          <button class="btn" onclick="sendToImpact()">${ICONS.star||''} Send to SNPImpact</button>
          <button class="btn" onclick="sendToCompare()">${ICONS.compare||ICONS.grid||''} Send to SNPCompare</button>
          <button class="btn" onclick="sendToTree()">${ICONS.tree} Send data to SNPTree</button>
          <button class="btn" onclick="openPangenomeRegion()">${ICONS.grid||''} Pangenome viewer ↗</button>
          <button class="btn" onclick="downloadVCF()">${ICONS.download} Download VCF</button>
        </div>
      </div>
    </div>
    <div class="result-tabs">
      <button class="rtab active" data-rt="table" onclick="switchRT('table')">${ICONS.table} Table view</button>
    </div>
    <div id="rtBody"></div>`;
  switchRT('table');
}
function switchRT(rt){
  document.querySelectorAll('.rtab').forEach(t=>t.classList.toggle('active',t.dataset.rt===rt));
  return renderTable();
}
/* hand the generated VCF matrix (+ metadata) to SNPTree for a local IBS phylogeny */
function sendToTree(){
  if(S.results&&S.results.rows&&S.results.rows.length){
    S.treeInput={
      rows:S.results.rows, accs:S.results.accs,
      chr:S.chr, start:Math.min(S.start,S.end), end:Math.max(S.start,S.end),
      dataset:S.dataset,
      datasetName:(Data.datasets().find(d=>d.id===S.dataset)||{}).name||S.dataset,
      vcfUrl:S.results.vcfUrl
    };
  }
  go('snptree');   // navigates even with no result yet (SNPTree shows a guided empty state)
}
/* hand the region's genotype matrix to SNPMatrix (IBS distance heatmap) */
function sendToMatrix(){
  if(S.results&&S.results.rows&&S.results.rows.length){
    S.matrixInput={
      rows:S.results.rows, accs:S.results.accs,
      chr:S.chr, start:Math.min(S.start,S.end), end:Math.max(S.start,S.end),
      dataset:S.dataset,
      datasetName:(Data.datasets().find(d=>d.id===S.dataset)||{}).name||S.dataset,
      vcfUrl:S.results.vcfUrl
    };
  }
  go('snpmatrix');
}
/* hand the generated VCF matrix (+ metadata) to SNPCompare for local vs global IBS */
function sendToCompare(){
  if(S.results&&S.results.rows&&S.results.rows.length){
    S.compareInput={
      rows:S.results.rows, accs:S.results.accs,
      chr:S.chr, start:Math.min(S.start,S.end), end:Math.max(S.start,S.end),
      dataset:S.dataset,
      datasetName:(Data.datasets().find(d=>d.id===S.dataset)||{}).name||S.dataset,
      vcfUrl:S.results.vcfUrl
    };
  }
  go('snpcompare');
}
/* hand the region's variants to SNPImpact for ranking (accessions irrelevant there) */
function sendToImpact(){
  if(S.results&&S.results.rows&&S.results.rows.length){
    S.impactInput={
      rows:S.results.rows, accs:S.results.accs,
      chr:S.chr, start:Math.min(S.start,S.end), end:Math.max(S.start,S.end),
      dataset:S.dataset,
      datasetName:(Data.datasets().find(d=>d.id===S.dataset)||{}).name||S.dataset,
      vcfUrl:S.results.vcfUrl
    };
  }
  go('snpimpact');
}

function setMaf(el){
  let v=parseFloat(el.value);
  if(isNaN(v)) v=0;
  v=Math.max(0, Math.min(0.5, v));      // clamp to valid MAF range
  v=Math.round(v*100)/100;              // snap to 0.01
  S.fMaf=v; S.page=1; renderTable();
}
function renderTable(){
  const {rows,accs}=S.results;
  // accession header height scales to the longest full ID so it isn't clipped
  const maxIdLen=accs.length?Math.max(...accs.map(a=>String(a.id).length)):8;
  const thH=Math.max(118, Math.min(300, Math.round(maxIdLen*6.4)+30));
  // filter
  let fr=rows.filter(r=>
    (S.fImpact==='all'||r.impact===S.fImpact) &&
    (S.fEffect==='all'||r.effect===S.fEffect) &&
    (r.maf>=S.fMaf));
  const perPage=S.perPage, pages=Math.max(1,Math.ceil(fr.length/perPage));
  if(S.page>pages)S.page=1;
  const slice=fr.slice((S.page-1)*perPage, S.page*perPage);
  const effects=['all',...new Set(rows.map(r=>r.effect))];

  const b=document.getElementById('rtBody');
  b.innerHTML=`
    <div class="filterbar">
      <div class="fld"><label>Effect impact</label>
        <select onchange="S.fImpact=this.value;S.page=1;renderTable()">
          ${['all','HIGH','MODERATE','LOW','MODIFIER'].map(v=>`<option value="${v}" ${v===S.fImpact?'selected':''}>${v==='all'?'All impacts':v}</option>`).join('')}
        </select></div>
      <div class="fld"><label>Effect type</label>
        <select onchange="S.fEffect=this.value;S.page=1;renderTable()">
          ${effects.map(v=>`<option value="${v}" ${v===S.fEffect?'selected':''}>${v==='all'?'All effects':v}</option>`).join('')}
        </select></div>
      <div class="fld"><label>Min. MAF (0.0–0.5)</label>
        <input type="number" min="0" max="0.5" step="0.01" value="${S.fMaf.toFixed(2)}" inputmode="decimal"
          title="Enter a value between 0.0 and 0.5" onchange="setMaf(this)"
          onkeydown="if(event.key==='Enter'){this.blur();}">
      </div>
    </div>
    <div class="legend">
      <span style="font-weight:600;color:var(--ink)">Genotype</span>
      <span class="li"><span class="sw" style="background:#e9f4ec"></span>Reference allele (0)</span>
      <span class="li"><span class="sw" style="background:#fdf3d2"></span>Alternative heterozygous allele (1)</span>
      <span class="li"><span class="sw" style="background:#f4a259"></span>Alternative homozygous allele (2)</span>
      <span class="li"><span class="sw" style="background:#eaeef4"></span>Missing (.)</span>
      <span style="margin-left:14px;font-weight:600;color:var(--ink)">LM score</span>
      <span class="li"><span class="sw" style="width:54px;background:linear-gradient(90deg,rgb(255,0,0),#aab4be,rgb(0,255,0))"></span>deleterious → tolerated</span>
    </div>
    ${genesPanel()}
    <div class="tbl-wrap"><table class="vcf">
      <thead><tr>
        <th>CHR</th><th>POS</th><th class="num">REF</th><th class="num">ALT</th>
        <th>Gene model</th><th>Effect</th><th>SNPEff Impact</th><th>Domain</th>
        <th class="num">MQ</th><th class="num">COMP</th><th class="num">maxR²</th><th class="num">MAF</th>
        <th class="num">PlantCAD1</th><th class="num">PlantCAD2</th><th class="num">ESM1</th><th class="num">ESM2</th><th class="num">ESM3</th>
        ${accs.map(a=>`<th class="acc-th" style="height:${thH}px" title="${escAttr((a.projTitle||a.proj||'')+' — '+a.id)}"><span class="proj-bar" style="background:${a.projColor};height:8px" title="${escAttr(a.projTitle||a.proj||'')}"></span><span class="v">${a.id}</span></th>`).join('')}
      </tr></thead>
      <tbody>
        ${slice.map(r=>rowHTML(r)).join('')}
      </tbody>
    </table></div>
    <div class="pager">
      <div class="info">Showing <span class="mono">${slice.length}</span> of <span class="mono">${fr.length}</span> loci · page <span class="mono">${S.page}/${pages}</span></div>
      <div class="pages">
        <button onclick="setPage(${Math.max(1,S.page-1)})">‹</button>
        ${pageBtns(pages)}
        <button onclick="setPage(${Math.min(pages,S.page+1)})">›</button>
      </div>
    </div>`;
  attachTT();
}
function rowHTML(r){
  const lo=r.pos-10000,hi=r.pos+10000;
  const link=`https://jbrowse.maizegdb.org/?data=B73&loc=${S.chr}:${lo}..${hi}&highlight=${S.chr}:${r.pos}..${r.pos}`;
  const effectText = String(r.effect||'');
  const isMis = /missense/i.test(effectText);
  const isSyn = /synonymous/i.test(effectText) && !/non[-_ ]?synonymous/i.test(effectText);
  const isFoldCoding = /missense|protein[_ -]?altering|non[_ -]?synonymous|stop[_ -]?gained|nonsense|frameshift|start[_ -]?lost|initiator[_ -]?codon|stop[_ -]?lost|inframe[_ -]?(insertion|deletion)/i.test(effectText);
  const hasGene = r.gene && r.gene!=='—';
  const peJump = (isMis && hasGene)
    ? ` <a class="pe-jump" href="#" title="View this substitution in PanEffect" onclick="goPanEffect('${r.gene}',{variant:'${escAttr(r.sub||'')}'});return false;">effects ↗</a>`
    : '';
  /* SNPFold: residue-level coding changes and severe LOF changes get a structure
     link. The genomic position/alleles are also passed so stop-gained records can
     still be selected when the VCF SUB field lacks a protein substitution. */
  const foldJump = (isFoldCoding && hasGene)
    ? ` <a class="fold-jump" href="#" title="Show this variant on the predicted protein structure in SNPFold" onclick="goFold('${r.gene}',{chr:'${escAttr(S.chr)}',pos:${r.pos},ref:'${escAttr(r.ref)}',alt:'${escAttr(r.alt)}',sub:'${escAttr(r.sub||'')}',effect:'${escAttr(r.effect||'')}'});return false;">fold ↗</a>`
    : '';
  const eff = (r.sub?`<span class="sub">(${r.sub})</span> ${r.effect}`:r.effect) + peJump + foldJump;
  const sc=(v)=>`<td class="score ${v===null?'na':''}" style="${v===null?'':'background:'+gColor(v)}">${v===null?'N/A':v}</td>`;
  return `<tr>
    <td class="c-mono" style="padding-left:11px">${S.chr.replace('chr','')}</td>
    <td class="c-pos"><a class="gene-link" href="${link}" target="_blank" rel="noopener">${r.pos.toLocaleString()}</a></td>
    <td class="c-allele c-ref" data-tt="${escAttr(alleleTT(r.ref,'REF allele'))}">${escAttr(alleleDisp(r.ref))}</td>
    <td class="c-allele c-alt" data-tt="${escAttr(alleleTT(r.alt,'ALT allele'))}">${escAttr(alleleDisp(r.alt))}</td>
    <td>${r.gene}</td>
    <td class="effect-cell">${eff}</td>
    <td><span class="pill ${r.impact.toLowerCase()}">${r.impact}</span></td>
    <td>${domTag(r.domain)}</td>
    <td class="num">${r.mq}</td>
    <td class="num">${r.comp}</td>
    <td class="num">${r.r2===null?'<span style="color:var(--faint)">NA</span>':r.r2}</td>
    <td class="num">${r.maf}</td>
    ${sc(r.pc1)}${sc(r.pc2)}${sc(r.esm1)}${sc(r.esm2)}${sc(r.esm3)}
    ${r.gts.map(g=>{
      const cls=g==='0/0'?'gt-00':(g==='1/1'?'gt-11':(g==='./.'?'gt-na':'gt-01'));
      const v=g==='0/0'?'0':(g==='1/1'?'2':(g==='./.'?'·':'1'));
      return `<td class="gt ${cls}" data-tt="${g}">${v}</td>`;
    }).join('')}
  </tr>`;
}
function genesPanel(){
  const rows=(S.results&&S.results.rows)||[];
  const genes=[...new Set(rows.map(r=>r.gene).filter(g=>g&&g!=='—'))].sort();
  if(!genes.length) return '';
  const jb=g=>`https://jbrowse.maizegdb.org/index.html?data=B73&loc=${encodeURIComponent(g)}`;
  const items=genes.map(g=>`<div style="display:flex;justify-content:space-between;align-items:center;gap:12px;padding:5px 2px;border-bottom:1px solid #eef1f5">
      <span class="c-mono" style="font-size:12.5px">${g}</span>
      <span style="display:flex;gap:12px;font-size:12px;white-space:nowrap">
        <a href="${jb(g)}" target="_blank" rel="noopener">JBrowse ↗</a>
        <a href="${pangenomeGeneURL(g)}" target="_blank" rel="noopener">Pangenome ↗</a>
        <a href="#" onclick="goPanEffect('${g}');return false;">PanEffect →</a>
        <a href="#" onclick="goFold('${g}');return false;">SNPFold →</a>
        <a href="#" onclick="goFunction('${g}');return false;">SNPFunction →</a>
      </span>
    </div>`).join('');
  return `<details class="card pad" style="margin-bottom:14px">
    <summary style="cursor:pointer;font-weight:600">${genes.length} gene model${genes.length>1?'s':''} in this region</summary>
    <div style="display:grid;grid-template-columns:1fr;gap:0;margin-top:10px">${items}</div>
  </details>`;
}
/* ---- MaizeGDB Pangenome viewer links (B73 v5 coordinates / gene models) ---- */
const PANGENOME_BASE = 'https://pangenome-viewer.maizegdb.org/';
/* link for a single B73 gene model */
function pangenomeGeneURL(gene, set){
  return `${PANGENOME_BASE}?set=${encodeURIComponent(set||'NAM')}&geneID=${encodeURIComponent(gene)}`;
}
/* link for a genomic interval */
function pangenomeRegionURL(chr, start, end, set){
  const lo=Math.min(start,end), hi=Math.max(start,end);
  const c=String(chr).startsWith('chr')?chr:`chr${chr}`;
  return `${PANGENOME_BASE}?set=${encodeURIComponent(set||'NAM')}&chr=${encodeURIComponent(c)}&start=${lo}&end=${hi}`;
}
/* open the current query region in the Pangenome viewer */
function openPangenomeRegion(){
  window.open(pangenomeRegionURL(S.chr,S.start,S.end),'_blank','noopener');
}
/* Jump to SNPFold for a specific gene model.
   `variant` is optional: {chr,pos,ref,alt,sub,effect} (or a plain "A123T" string).
   When present SNPFold selects that site — lollipop, table row, 3D label — after
   the gene is shown. Passing it through state (not a reload) means clicking fold
   links repeatedly on the same gene never re-queries the HDF5 store. */
function goFold(gene, variant){
  S.foldGene=gene;
  S.foldVariant = variant
    ? (typeof variant==='string' ? {sub:variant} : variant)
    : null;
  go('snpfold');
}
/* jump to SNPFunction (gene function & allele mining) for a specific gene */
function goFunction(gene, dataset){ S.functionGene=gene; S.functionDataset=dataset||S.dataset; go('snpfunction'); }
function pageBtns(pages){
  let out='';const cur=S.page;
  const show=new Set([1,pages,cur,cur-1,cur+1,cur-2,cur+2]);
  let last=0;
  for(let i=1;i<=pages;i++){
    if(!show.has(i))continue;
    if(i-last>1)out+='<span style="align-self:center;color:var(--faint);padding:0 2px">…</span>';
    out+=`<button class="${i===cur?'on':''}" onclick="setPage(${i})">${i}</button>`;
    last=i;
  }
  return out;
}
function setPage(p){S.page=p;renderTable();document.querySelector('.tbl-wrap').scrollTop=0;}

function treePlaceholder(){
  const accs=S.results.accs.slice(0,14);
  // simple radial-ish NJ mock using SVG
  const cx=300,cy=210,R=160;
  const leaves=accs.map((a,i)=>{
    const ang=(i/accs.length)*Math.PI*2-Math.PI/2;
    const x=cx+Math.cos(ang)*R, y=cy+Math.sin(ang)*R;
    const mx=cx+Math.cos(ang)*(R*.55), my=cy+Math.sin(ang)*(R*.55);
    return {a,x,y,mx,my,ang};
  });
  return `<div class="card pad fade">
    <div style="display:flex;gap:24px;flex-wrap:wrap;align-items:center">
      <svg viewBox="0 0 600 420" style="flex:1;min-width:320px;max-width:600px">
        ${leaves.map(l=>`<path d="M${cx} ${cy} Q ${l.mx} ${l.my} ${l.x} ${l.y}" stroke="${l.a.projColor}" stroke-width="1.6" fill="none" opacity=".75"/>`).join('')}
        ${leaves.map(l=>`<circle cx="${l.x}" cy="${l.y}" r="4" fill="${l.a.projColor}"/>
          <text x="${l.x+(Math.cos(l.ang)>=0?8:-8)}" y="${l.y+3}" font-size="10" font-family="IBM Plex Mono" fill="#0f1b2d" text-anchor="${Math.cos(l.ang)>=0?'start':'end'}">${l.a.founder}</text>`).join('')}
        <circle cx="${cx}" cy="${cy}" r="5" fill="#13264a"/>
      </svg>
      <div style="flex:1;min-width:240px">
        <h3 style="font-family:var(--disp);margin:0 0 8px">Neighbor-joining tree</h3>
        <p style="color:var(--muted);margin:0 0 12px">A local phylogeny built from the variants in this region across your selected accessions. In the full tool this is interactive — reroot, collapse clades, and color by bioproject.</p>
        <p style="color:var(--faint);font-size:12px;margin:0">Tree computation is powered by VCF2PopTree on the parsed VCF. Showing first ${accs.length} accessions for preview.</p>
      </div>
    </div>
  </div>`;
}
function ibsPlaceholder(){
  const accs=S.results.accs.slice(0,12);
  const n=accs.length;
  const cell=Math.min(34, Math.floor(520/Math.max(n,1)));
  let g='';
  for(let i=0;i<n;i++)for(let j=0;j<n;j++){
    const v=i===j?1:rnd(.55,.99);
    const col=`rgba(37,99,235,${(v-.5).toFixed(2)})`;
    g+=`<div data-tt="${accs[i].founder} vs ${accs[j].founder}: ${(v*100|0)}% IBS" style="width:${cell}px;height:${cell}px;background:${col};border:1px solid #fff"></div>`;
  }
  return `<div class="card pad fade">
    <h3 style="font-family:var(--disp);margin:0 0 4px">Identity-by-state matrix</h3>
    <p style="color:var(--muted);margin:0 0 16px">Pairwise sequence identity across selected accessions for this region. Darker = more similar.</p>
    <div style="display:grid;grid-template-columns:repeat(${n},${cell}px);width:max-content">${g}</div>
  </div>`;
}

/* register with the suite shell */
SNPTools.register('snpversity', { render(){ renderVersity(); } });
