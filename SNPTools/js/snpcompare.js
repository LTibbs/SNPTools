/* =====================================================================
 *  snpcompare.js — SNPCompare: identity-by-state similarity viewer.
 *
 *  Scopes
 *    • Genome-wide (global): precomputed, served per-focal-accession by
 *      ibsCompare.php  (only the focal row is sent to the browser).
 *    • This region (local): computed in-browser from a SNPVersity result.
 *    • Both: global vs local side-by-side, Δ = local − global.
 *
 *  Views
 *    • Table       — ranked list (all scopes).
 *    • Matrix      — all-pairs heatmap, split triangles, clustered order
 *                    and dendrogram (region scope only).
 *    • PCoA        — classical MDS of the distance matrix (region only).
 *    • Focal spread— 1D beeswarm of every similarity to the focal line.
 *    • Diagnostics — co-called sites vs similarity, exposes the coverage
 *                    confound directly (region only). [TEMPORARILY DISABLED —
 *                    the button and dispatch are commented out; renderDiag()
 *                    is kept intact so it can be restored in one edit.]
 *
 *  Running
 *    Normal navigation does not run an analysis on page load or when options
 *    change. Configure the focal accession, scope and filters, then press a
 *    view button (Table / Matrix / PCoA / Focal spread) to submit. An explicit
 *    tool handoff may request and autorun a view. Once a result is shown, its
 *    refinement controls (scale, order, mask, sort, legend) update it live.
 *
 *  Definitions
 *    Similarity  = fraction of co-called SNP sites (both non-missing)
 *                  that have the same genotype.
 *    Missing%    = percent of sites where at least one of the pair is
 *                  missing a genotype.
 *    Unmeasured  = fewer than `minSites` co-called sites. These are NOT
 *                  similarity 0. They are masked out of the colour ramp
 *                  and drawn in a neutral hatch so that coverage gaps
 *                  can never be mistaken for divergence.
 *
 *  Colour scaling
 *    Similarity between maize accessions lives in a narrow band near the
 *    top of 0–1, so a fixed 0–1 ramp throws away nearly all structure.
 *    Scales offered: log distance (default), percentile-clipped linear,
 *    rank/quantile, z-score (diverging) and raw 0–1 for reference. The
 *    legend is a histogram of the data with draggable endpoints.
 *
 *  Depends on core.js (SNPTools, S, go, ICONS) and data.js (Data).
 * ===================================================================== */
const SNPCompare = (function () {

  const CFG = {
    globalEndpoint : 'ibsCompare.php', // ?focal=<ID> -> {rows:[{id,similarity,missing}]}
    useDemoGlobal  : false,
    defaultDataset : 'mgdb2026_hq',
    minSitesFloor  : 20,               // never mask below this many sites
    minSitesFrac   : 0.05,             // default mask threshold = 5% of sites
  };

  const ST = {
    mode:'global', view:'table', dataset:CFG.defaultDataset, focal:null,
    sortKey:'sim', sortDir:-1,
    fSimMin:null, fSimMax:null, fMissMax:null, fProj:'all',
    /* analysis is only (re)run when a view button is pressed — nothing runs
       on page load or when options change. `ran` tracks whether a result is
       currently on screen; `ranFocal`/`ranMode` remember what produced it. */
    ran:false, ranFocal:null, ranMode:null, force:false,
    mdsLabels:'auto',      // PCoA label strategy: auto | all | focal | none
    input:null,            // local hand-off {rows, accs, chr, start, end, dataset,...}
    globalCache:{},        // focalId -> {rows, demo}
    allRows:[], gdemo:false, _meta:null, _metaDs:null,

    /* --- matrix / scale state --- */
    scale:'logd',          // logd | clip | rank | z | fixed
    lower:'cocall',        // lower triangle: cocall | sim
    order:'cluster',       // cluster | name | sim
    minSites:null,         // null = auto from CFG
    dropMiss:null,         // hide accessions whose mean missing% exceeds this
    dom:null,              // {lo,hi} user override of the colour domain
    _pairs:null, _pairsKey:null, _scale:null, _layout:null,
  };

  /* ---------------- genotype + IBS ----------------
     gts entries are 1-byte codes from Data.parseVcf: 0/1/2 = dosage, 3 = missing. */
  function dose(g){ return (g == null || g === 3) ? null : g; }

  function minSites(total){
    if(ST.minSites!=null) return ST.minSites;
    return Math.max(CFG.minSitesFloor, Math.round(CFG.minSitesFrac*(total||0)));
  }

  /* focal-vs-all, used by the table. Masked pairs return sim:null, never 0. */
  function localCompute(input, focalId){
    if(!input) return null;
    const P=allPairs(input); if(!P) return null;
    const fi=P.ids.indexOf(focalId); if(fi<0) return null;
    const thr=minSites(P.total), out=[];
    for(let o=0;o<P.n;o++){
      const k=fi*P.n+o, both=P.both[k], s=P.sim[k];
      out.push({ id:P.ids[o],
                 sim: (both>=thr && !isNaN(s)) ? s : null,
                 raw: isNaN(s)?null:s,
                 miss: P.miss[k], both });
    }
    return out;
  }

  /* full n×n IBS from a SNPVersity result, cached on the hand-off */
  function pairKey(input){
    return [input.dataset,input.chr,input.start,input.end,
            input.rows.length,input.accs.map(a=>a.id).join('|')].join('~');
  }
  function allPairs(input){
    if(!input||!input.accs||!input.rows) return null;
    const key=pairKey(input);
    if(ST._pairs && ST._pairsKey===key) return ST._pairs;

    const accs=input.accs, n=accs.length, m=input.rows.length;
    // transpose to per-accession dose vectors; -1 = missing
    const D=[]; for(let i=0;i<n;i++) D.push(new Int8Array(m));
    for(let s=0;s<m;s++){
      const g=input.rows[s].gts;
      for(let i=0;i<n;i++){ const d=dose(g[i]); D[i][s]= (d==null? -1 : d); }
    }
    const sim=new Float64Array(n*n), both=new Int32Array(n*n), miss=new Float64Array(n*n);
    for(let i=0;i<n;i++){
      const a=D[i];
      for(let j=i;j<n;j++){
        const b=D[j]; let bo=0, mt=0;
        for(let s=0;s<m;s++){ const x=a[s], y=b[s];
          if(x<0||y<0) continue; bo++; if(x===y) mt++; }
        const v = bo? mt/bo : NaN, ms = m? 100*(m-bo)/m : 0;
        sim[i*n+j]=sim[j*n+i]=v; both[i*n+j]=both[j*n+i]=bo; miss[i*n+j]=miss[j*n+i]=ms;
      }
    }
    ST._pairs={ n, total:m, ids:accs.map(a=>a.id), sim, both, miss };
    ST._pairsKey=key; ST._layout=null;
    return ST._pairs;
  }

  /* accessions kept after the "drop low-coverage lines" filter */
  function keptIdx(P){
    const idx=[];
    for(let i=0;i<P.n;i++){
      if(ST.dropMiss!=null){
        let s=0,c=0; for(let j=0;j<P.n;j++){ if(j===i)continue; s+=P.miss[i*P.n+j]; c++; }
        if(c && s/c > ST.dropMiss) continue;
      }
      idx.push(i);
    }
    return idx;
  }

  /* ---------------- colour scales ---------------- */
  /* every transform is monotonically increasing in similarity, so the hot
     end of the ramp always means "more similar" regardless of scale.     */
  const SEQ=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];        // viridis
  const DIVP=[[5,48,97],[103,169,207],[247,247,247],[239,138,98],[178,24,43]];    // blue→red
  function ramp(pal,t){
    t=Math.max(0,Math.min(1,t));
    const x=t*(pal.length-1), i=Math.min(pal.length-2,Math.floor(x)), f=x-i;
    const a=pal[i], b=pal[i+1];
    return `rgb(${Math.round(a[0]+(b[0]-a[0])*f)},${Math.round(a[1]+(b[1]-a[1])*f)},${Math.round(a[2]+(b[2]-a[2])*f)})`;
  }
  const MASK_FILL='#e6e9ef', MASK_LINE='#b6bdc9', DIAG_FILL='#39424f';

  const SCALES={
    logd : {label:'Log distance',        hint:'color ∝ −log₁₀(1 − similarity); spreads the crowded top end'},
    clip : {label:'Percentile clipped',  hint:'linear similarity, domain = 2nd–98th percentile of observed pairs'},
    rank : {label:'Rank / quantile',     hint:'color = rank among all pairs; maximal contrast, relative meaning'},
    z    : {label:'Z-score (diverging)', hint:'standard deviations from the panel mean; centred palette'},
    fixed: {label:'Raw 0–1',             hint:'the unscaled ramp, for reference'},
  };

  /* `res` = smallest distance the data can actually resolve (≈ 0.5 / median
     co-called sites). Without it, identical pairs pin the log scale at an
     arbitrary floor and squash everything else into the dark end.        */
  function buildScale(sims, res){
    const floor=Math.min(0.05, Math.max(1e-6, res||1e-4));
    const s=sims.filter(v=>v!=null && isFinite(v)).sort((a,b)=>a-b);
    const n=s.length;
    const q=p=> n? s[Math.min(n-1,Math.max(0,Math.round(p*(n-1))))] : 0;
    let mu=0; for(const v of s) mu+=v; mu = n? mu/n : 0;
    let sd=0; for(const v of s) sd+=(v-mu)*(v-mu); sd = n>1? Math.sqrt(sd/(n-1)) : 1; if(!sd) sd=1e-9;

    let tf, inv, div=false;
    switch(ST.scale){
      case 'fixed': tf=v=>v; inv=t=>t; break;
      case 'clip':  tf=v=>v; inv=t=>t; break;
      case 'rank':  tf=v=>{ let lo=0,hi=n; while(lo<hi){const mid=(lo+hi)>>1; if(s[mid]<v)lo=mid+1;else hi=mid;}
                            return n>1? lo/(n-1) : 0.5; };
                    inv=t=>q(t); break;
      case 'z':     tf=v=>(v-mu)/sd; inv=t=>mu+t*sd; div=true; break;
      default:      tf=v=>-Math.log10(Math.max(1-v,floor)); inv=t=>1-Math.pow(10,-t);
    }
    const T=s.map(tf);
    let lo,hi;
    if(ST.scale==='fixed'){ lo=0; hi=1; }
    else if(ST.scale==='rank'){ lo=0; hi=1; }
    else if(ST.scale==='z'){ const a=Math.abs(T[Math.round(0.02*(n-1))]||0), b=Math.abs(T[Math.round(0.98*(n-1))]||0);
                             hi=Math.max(a,b,1e-9); lo=-hi; }
    else { lo=T[Math.max(0,Math.round(0.02*(n-1)))]; hi=T[Math.min(n-1,Math.round(0.98*(n-1)))]; }
    if(!(hi>lo)){ hi=lo+1e-6; }
    if(ST.dom){ lo=ST.dom.lo; hi=ST.dom.hi; }

    const pal = div? DIVP : SEQ;
    const color=v=>{ if(v==null||!isFinite(v)) return MASK_FILL;
                     return ramp(pal,(tf(v)-lo)/(hi-lo)); };
    return { tf, inv, lo, hi, div, pal, color, values:T, raw:s,
             tmin: T.length? T[0]:0, tmax: T.length? T[T.length-1]:1 };
  }

  /* ---------------- clustering (average linkage, UPGMA) ---------------- */
  function cluster(P, idx){
    const k=idx.length;
    if(k<2) return { order:idx.slice(), nodes:[] };
    // distance, masked pairs imputed with the mean observed distance
    const thr=minSites(P.total);
    const d=[]; let sum=0,cnt=0;
    for(let a=0;a<k;a++){ d.push(new Float64Array(k)); }
    for(let a=0;a<k;a++) for(let b=a+1;b<k;b++){
      const kk=idx[a]*P.n+idx[b], ok=P.both[kk]>=thr && !isNaN(P.sim[kk]);
      const v= ok? 1-P.sim[kk] : NaN;
      d[a][b]=d[b][a]=v; if(ok){ sum+=v; cnt++; }
    }
    const fill = cnt? sum/cnt : 0.2;
    for(let a=0;a<k;a++) for(let b=0;b<k;b++) if(a!==b && isNaN(d[a][b])) d[a][b]=fill;

    const active=[];
    for(let a=0;a<k;a++){ active.push({leaves:[idx[a]], size:1, id:a, node:{leaf:idx[a],height:0}}); }
    const dist=new Map(); const key=(x,y)=>x<y? x+','+y : y+','+x;
    for(let a=0;a<k;a++) for(let b=a+1;b<k;b++) dist.set(key(a,b), d[a][b]);
    let next=k;
    while(active.length>1){
      let bi=0,bj=1,bd=Infinity;
      for(let a=0;a<active.length;a++) for(let b=a+1;b<active.length;b++){
        const v=dist.get(key(active[a].id,active[b].id));
        if(v<bd){ bd=v; bi=a; bj=b; }
      }
      const A=active[bi], B=active[bj];
      const merged={ id:next++, size:A.size+B.size, leaves:A.leaves.concat(B.leaves),
                     node:{ left:A.node, right:B.node, height:bd } };
      for(const C of active){ if(C===A||C===B) continue;
        const va=dist.get(key(A.id,C.id)), vb=dist.get(key(B.id,C.id));
        dist.set(key(merged.id,C.id), (va*A.size+vb*B.size)/(A.size+B.size));
      }
      active.splice(Math.max(bi,bj),1); active.splice(Math.min(bi,bj),1);
      active.push(merged);
    }
    return { order:active[0].leaves, root:active[0].node };
  }

  function layout(P){
    const key=[ST.order,ST.dropMiss,ST.minSites,ST.order==='sim'?ST.focal:''].join('~');
    if(ST._layout && ST._layout.key===key) return ST._layout;
    const idx=keptIdx(P);
    let order=idx.slice(), root=null;
    if(ST.order==='cluster'){ const c=cluster(P,idx); order=c.order; root=c.root; }
    else if(ST.order==='name'){ const meta=metaMap(ST.dataset);
      order.sort((a,b)=>((meta[P.ids[a]]||{}).name||P.ids[a]).localeCompare((meta[P.ids[b]]||{}).name||P.ids[b],undefined,{numeric:true})); }
    else if(ST.order==='sim'){ const fi=P.ids.indexOf(ST.focal);
      if(fi>=0) order.sort((a,b)=>{ const va=P.sim[fi*P.n+a], vb=P.sim[fi*P.n+b];
        return (isNaN(vb)?-1:vb)-(isNaN(va)?-1:va); }); }
    ST._layout={ key, order, root };
    return ST._layout;
  }

  /* ---------------- PCoA (classical MDS) ---------------- */
  function pcoa(P, order){
    const k=order.length; if(k<3) return null;
    const thr=minSites(P.total);
    const D2=[]; let sum=0,cnt=0;
    for(let a=0;a<k;a++) D2.push(new Float64Array(k));
    for(let a=0;a<k;a++) for(let b=a+1;b<k;b++){
      const kk=order[a]*P.n+order[b], ok=P.both[kk]>=thr && !isNaN(P.sim[kk]);
      const v= ok? 1-P.sim[kk] : NaN;
      D2[a][b]=D2[b][a]=v; if(ok){ sum+=v; cnt++; }
    }
    const fill=cnt? sum/cnt : 0.2;
    for(let a=0;a<k;a++) for(let b=0;b<k;b++) D2[a][b]= a===b?0:(isNaN(D2[a][b])?fill:D2[a][b])**2;
    // double centre
    const rm=new Float64Array(k); let gm=0;
    for(let a=0;a<k;a++){ let s=0; for(let b=0;b<k;b++) s+=D2[a][b]; rm[a]=s/k; gm+=s; }
    gm/=k*k;
    const B=[]; for(let a=0;a<k;a++){ B.push(new Float64Array(k));
      for(let b=0;b<k;b++) B[a][b]=-0.5*(D2[a][b]-rm[a]-rm[b]+gm); }
    const eig=(M)=>{
      let v=new Float64Array(k); for(let i=0;i<k;i++) v[i]=Math.sin(i*1.7+0.3);
      let lam=0;
      for(let it=0;it<300;it++){
        const w=new Float64Array(k);
        for(let a=0;a<k;a++){ let s=0; const Ma=M[a]; for(let b=0;b<k;b++) s+=Ma[b]*v[b]; w[a]=s; }
        let nrm=0; for(let i=0;i<k;i++) nrm+=w[i]*w[i]; nrm=Math.sqrt(nrm);
        if(nrm<1e-12) break;
        for(let i=0;i<k;i++) w[i]/=nrm; lam=nrm; v=w;
      }
      return {lam, v};
    };
    const e1=eig(B);
    const B2=B.map((row,a)=>{ const r=new Float64Array(k);
      for(let b=0;b<k;b++) r[b]=row[b]-e1.lam*e1.v[a]*e1.v[b]; return r; });
    const e2=eig(B2);
    let tr=0; for(let a=0;a<k;a++) tr+=B[a][a];
    const pts=order.map((gi,a)=>({ gi,
      x:e1.v[a]*Math.sqrt(Math.max(e1.lam,0)),
      y:e2.v[a]*Math.sqrt(Math.max(e2.lam,0)) }));
    return { pts, v1: tr>0? e1.lam/tr:0, v2: tr>0? e2.lam/tr:0 };
  }

  /* ---------------- metadata join ---------------- */
  function metaMap(ds){
    if(ST._meta && ST._metaDs===ds) return ST._meta;
    const projBio={};
    (Data.projectsFor(ds)||[]).forEach(p=>{ projBio[p.id]=(p.bioprojects&&p.bioprojects.length)?p.bioprojects.join(', '):(p.title||''); });
    const m={};
    (Data.accessionsFor(ds)||[]).forEach(a=>{ m[a.id]={id:a.id, name:a.founder, run:a.run, proj:a.proj, projColor:a.projColor, bio:projBio[a.proj]||''}; });
    ST._meta=m; ST._metaDs=ds; return m;
  }
  function projectOptions(ds){
    const seen={}, out=[];
    (Data.projectsFor(ds)||[]).forEach(p=>{ const b=(p.bioprojects&&p.bioprojects.length)?p.bioprojects.join(', '):(p.title||p.id);
      if(!seen[b]){seen[b]=1; out.push({id:p.id, label:b});} });
    return out;
  }
  function labelFor(id){ const m=metaMap(ST.dataset)[id]; return (m&&m.name)? m.name : id; }

  /* ---------------- global source (endpoint + demo fallback) ---------------- */
  function hash(s){ let h=2166136261>>>0; for(let i=0;i<s.length;i++){ h^=s.charCodeAt(i); h=Math.imul(h,16777619);} return h>>>0; }
  function globalDemo(ds, focalId){
    return (Data.accessionsFor(ds)||[]).map(a=>{
      if(a.id===focalId) return {id:a.id, sim:1, miss:0};
      const h=hash(focalId+'|'+a.id);
      return {id:a.id,
        sim:+(0.85+(h%1500)/1500*0.149).toFixed(4),
        miss:+(((h>>>7)%3000)/100).toFixed(2)};
    });
  }
  async function getGlobal(ds, focalId){
    if(ST.globalCache[focalId]) return ST.globalCache[focalId];
    let res;
    if(CFG.useDemoGlobal){ res={rows:globalDemo(ds,focalId), demo:true}; }
    else {
      const resp=await fetch(`${CFG.globalEndpoint}?focal=${encodeURIComponent(focalId)}`,{cache:'no-store'});
      if(!resp.ok) throw new Error('ibsCompare.php failed (HTTP '+resp.status+')');
      const raw=await resp.text(); let j;
      try{ j=JSON.parse(raw); }catch(e){ throw new Error('ibsCompare.php did not return JSON:\n'+raw.slice(0,500)); }
      const rows=(j.rows|| (j.ids? j.ids.map((id,i)=>({id, similarity:j.sim[i], missing:j.miss[i]})):[]))
        .map(r=>({id:r.id, sim:+(r.similarity!=null?r.similarity:r.sim), miss:Math.max(0,+(r.missing!=null?r.missing:r.miss))}));
      res={rows, demo:false};
    }
    ST.globalCache[focalId]=res; return res;
  }

  /* ---------------- build the combined row set ---------------- */
  async function buildRows(){
    const ds=ST.dataset, meta=metaMap(ds);
    let g=null, l=null;
    if(ST.mode==='global'||ST.mode==='both') g=await getGlobal(ds, ST.focal);
    if(ST.mode==='local' ||ST.mode==='both') l=localCompute(ST.input, ST.focal);
    const byId={};
    const add=id=>{ if(!byId[id]) byId[id]=Object.assign({id}, meta[id]||{name:id,run:id,bio:''}); return byId[id]; };
    if(g) g.rows.forEach(r=>{ const x=add(r.id); x.gsim=r.sim; x.gmiss=r.miss; });
    if(l) l.forEach(r=>{ const x=add(r.id); x.lsim=r.sim; x.lmiss=r.miss; x.lboth=r.both; });
    let rows=Object.values(byId);
    if(ST.mode==='local') rows=rows.filter(r=>r.lboth!=null || r.id===ST.focal);
    if(ST.mode==='both')  rows=rows.filter(r=>r.lboth!=null);
    rows.forEach(r=>{ r.dsim=(r.lsim!=null&&r.gsim!=null)? r.lsim-r.gsim : null;
      r.sim = ST.mode==='local'? r.lsim : r.gsim; r.miss = ST.mode==='local'? r.lmiss : r.gmiss; });
    ST.allRows=rows; ST.gdemo=!!(g&&g.demo);
  }

  /* ---------------- sort + filter ---------------- */
  function viewRows(){
    let rows=ST.allRows.slice();
    const num=v=>v==null?-Infinity:v;
    if(ST.fSimMin!=null) rows=rows.filter(r=>num(r.sim)>=ST.fSimMin);
    if(ST.fSimMax!=null) rows=rows.filter(r=>num(r.sim)<=ST.fSimMax);
    if(ST.fMissMax!=null)rows=rows.filter(r=>(r.miss==null?0:r.miss)<=ST.fMissMax);
    if(ST.fProj!=='all') rows=rows.filter(r=>r.proj===ST.fProj);
    const k=ST.sortKey, d=ST.sortDir;
    rows.sort((a,b)=>{
      let va=a[k], vb=b[k];
      if(k==='id'||k==='name'||k==='run'||k==='bio'){ va=(va||'').toString(); vb=(vb||'').toString(); return d*va.localeCompare(vb,undefined,{numeric:true}); }
      va=va==null?-Infinity:va; vb=vb==null?-Infinity:vb; return d*(va-vb);
    });
    return rows;
  }

  /* ---------------- shell ---------------- */
  function globalAvailable(){ try{ return Data.familyOf(ST.dataset)==='mgdb2026'; }catch(e){ return true; } }
  function hasLocal(){ return !!(ST.input&&ST.input.accs&&ST.input.accs.length>1); }

  function render(){
    injectCSS();
    const crumb=document.getElementById('crumbTool'); if(crumb) crumb.innerHTML='<b>SNPCompare</b>';
    const page=document.getElementById('page');
    // Optional one-shot navigation request, e.g. from SNPTree. The request
    // selects a view/scope and can run it immediately with the same local
    // SNPVersity result held in S.compareInput.
    const request=S.compareRequest||null;
    S.compareRequest=null;
    const prevDs=ST.dataset;
    if(S.compareInput){ if(S.compareInput!==ST.input){ ST._pairs=null; ST._layout=null; ST.force=false; } ST.input=S.compareInput; }
    if(ST.input){ ST.dataset=ST.input.dataset||ST.dataset; }
    if(ST.dataset!==prevDs){ ST.focal=null; ST._meta=null; ST.globalCache={}; ST.allRows=[]; ST._pairs=null; }
    const gAvail=globalAvailable();
    if(!ST.focal){
      if(ST.input && ST.input.accs && ST.input.accs.length){
        ST.focal=ST.input.accs[0].id; ST.mode = gAvail ? 'both' : 'local';
      } else {
        const a=(Data.accessionsFor(ST.dataset)||[])[0]; ST.focal=a?a.id:null; ST.mode = gAvail ? 'global' : 'local';
      }
    }
    if(!gAvail && ST.mode!=='local') ST.mode='local';
    if(ST.focal){
      const inRegion = ST.input && ST.input.accs && ST.input.accs.some(a=>a.id===ST.focal);
      const inDataset = !!(Data.accessionById && Data.accessionById(ST.focal));
      const ok = ST.input ? (inRegion || (gAvail && inDataset)) : inDataset;
      if(!ok){
        ST.focal = (ST.input && ST.input.accs && ST.input.accs.length) ? ST.input.accs[0].id
                 : ((Data.accessionsFor(ST.dataset)||[])[0]||{}).id || null;
      }
    }
    if(request){
      if(request.mode==='local' && hasLocal()) ST.mode='local';
      else if(request.mode==='both' && hasLocal() && gAvail) ST.mode='both';
      else if(request.mode==='global' && gAvail) ST.mode='global';

      if(request.view==='matrix' && hasLocal()) ST.view='matrix';
      else if(request.view==='mds' && hasLocal()) ST.view='mds';
      else if(request.view==='spread') ST.view='spread';
      else ST.view='table';
    }
    if(!hasLocal() && ST.view!=='table' && ST.view!=='spread') ST.view='table';
    page.className='page fade';
    page.innerHTML=shell();
    bindResize();
    ST.ran=false;
    // Normal navigation stays idle; an explicit handoff can open and run a
    // requested view immediately.
    if(request && request.autorun) recompute(); else showIdle();
  }

  function shell(){
    const ds=ST.dataset;
    const gAvail=globalAvailable(), local=hasLocal();
    const ids=gAvail ? (Data.accessionsFor(ds)||[]).map(a=>a.id)
                     : (local ? ST.input.accs.map(a=>a.id) : []);
    const region=local?`${ST.input.chr}:${(+ST.input.start).toLocaleString()}–${(+ST.input.end).toLocaleString()}`:null;
    const dsName=ST.input?(ST.input.datasetName||ST.dataset):ST.dataset;
    return `
    <section class="sec"><div class="bar"></div><div style="width:100%">
      <h1>SNPCompare</h1>
      <p>Rank and map accessions by identity-by-state similarity to a focal accession — genome-wide and within a queried region.</p>
    </div></section>

    <div class="card pad" style="margin-bottom:16px">
      <div style="display:flex;gap:22px;flex-wrap:wrap;align-items:flex-end">
        <div style="min-width:480px;flex:1 1 480px">
          <div class="fl-lbl">Focal accession</div>
          <div style="display:flex;gap:8px">
            <input id="cmpFocal" list="cmpFocalList" value="${esc(ST.focal||'')}" placeholder="type a SNPVersity ID…"
              oninput="SNPCompare.syncFocal(this.value)"
              onkeydown="if(event.key==='Enter'){event.preventDefault();SNPCompare.runCurrent();}"
              style="flex:1;min-width:340px;border:1px solid var(--line);border-radius:9px;padding:9px 11px;font-family:var(--mono);font-size:13px">
            <datalist id="cmpFocalList">${ids.slice(0,4000).map(i=>`<option value="${esc(i)}">`).join('')}</datalist>
          </div>
        </div>
        <div>
          <div class="fl-lbl">Scope</div>
          ${scopeBtn('global','Genome-wide',gAvail,'Genome-wide matrix available for MaizeGDB 2026 only')}
          ${scopeBtn('local','This region',local,'Send a result from SNPVersity to enable')}
          ${scopeBtn('both','Both (Δ)',gAvail&&local, !gAvail?'Genome-wide matrix available for MaizeGDB 2026 only':'Send a result from SNPVersity to enable')}
        </div>
        ${region?`<div><div class="fl-lbl">Region</div><div class="c-mono" style="color:var(--blue-600);font-size:13px;padding:8px 0">${region}</div></div>`:''}
      </div>
      ${!gAvail?`<div class="mtx-note" style="margin-top:12px">Genome-wide precomputed IBS is available for <b>MaizeGDB 2026</b> only. For <b>${esc(dsName)}</b>, use <b>This region</b> — SNPCompare computes identity-by-state live from your SNPVersity result.</div>`:''}

      ${accListHTML()}

      <div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin-top:14px">
        <div class="fl-lbl" style="margin:0 4px 0 0">Run</div>
        ${viewBtn('table','Table',true,'')}
        ${viewBtn('matrix','Matrix',local,'Send a multi-accession result from SNPVersity to enable')}
        ${viewBtn('mds','PCoA',local,'Send a multi-accession result from SNPVersity to enable')}
        ${viewBtn('spread','Focal spread',true,'')}
        <!-- Diagnostics view temporarily disabled -->
        <span id="cmpPending" style="display:none;margin-left:8px;font-size:11.5px;color:#8a6d1e">Options changed — press a view to update.</span>
      </div>

      ${local?`<div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin-top:14px">
        <span style="margin-left:auto;display:flex;gap:8px;flex-wrap:wrap">
          <button class="qbtn" onclick="SNPCompare.exportMatrixCSV()">Export matrix</button>
          <button class="qbtn" onclick="SNPCompare.toTree()">Open region in SNPTree</button>
        </span>
      </div>`:''}
    </div>

    <div id="cmpBanner"></div>
    <div id="cmpCtl"></div>
    <div id="cmpCount" style="color:var(--muted);font-size:13px;margin:0 0 8px"></div>
    <div class="card pad fade" id="cmpViewWrap" style="overflow:auto;max-height:78vh"></div>
    <div id="cmpTip" class="cmp-tip" style="display:none"></div>`;
  }
  function scopeBtn(m,label,enabled,reason){
    return `<button class="qbtn ${ST.mode===m?'solid':''}" ${enabled?'':`disabled title="${esc(reason||'Unavailable')}"`}
      style="${enabled?'':'opacity:.45;cursor:not-allowed'}" onclick="SNPCompare.setMode('${m}')">${label}</button>`;
  }
  function viewBtn(v,label,enabled,reason){
    return `<button class="qbtn ${ST.view===v?'solid':''}" ${enabled?'':`disabled title="${esc(reason||'Unavailable')}"`}
      style="${enabled?'':'opacity:.45;cursor:not-allowed'}" onclick="SNPCompare.setView('${v}')">${label}</button>`;
  }

  /* ---------------- imported-accession list ---------------- */
  /* When a result was handed over from SNPVersity, list every accession in
     the local analysis (mirrors SNPVersity's "Selected" panel). Clicking a
     chip just populates the focal field — it does not run the analysis. */
  function accChips(){
    if(!(ST.input && ST.input.accs && ST.input.accs.length)) return '';
    const meta=metaMap(ST.dataset);
    return ST.input.accs.map(a=>{
      const id=a.id, m=meta[id]||{};
      const label=m.name||a.founder||a.run||a.label||id;
      const col=m.projColor||a.projColor||'#8894a6';
      const on=id===ST.focal;
      return `<span class="cmp-acc-chip${on?' on':''}" title="${esc(id)}" onclick="SNPCompare.pickFocal('${jsStr(id)}')"><span class="cmp-acc-dot" style="background:${col}"></span>${esc(label)}</span>`;
    }).join('');
  }
  function accListHTML(){
    if(!(ST.input && ST.input.accs && ST.input.accs.length)) return '';
    const n=ST.input.accs.length;
    return `<div class="cmp-acc-wrap">
      <div class="fl-lbl" style="text-transform:none;letter-spacing:0;margin-bottom:0">
        <span style="text-transform:uppercase;letter-spacing:.4px">Accessions in this result</span>
        <span class="c-mono" style="color:var(--blue-600,#274b8f)"> · ${n}</span>
        <span style="font-weight:400;color:var(--muted)"> — click one to set it as the focal accession</span>
      </div>
      <div class="cmp-acc-box" id="cmpAccList">${accChips()}</div>
    </div>`;
  }
  /* re-highlight the chip list in place (no re-render of the whole shell) */
  function refreshAccList(){ const box=document.getElementById('cmpAccList'); if(box) box.innerHTML=accChips(); }

  /* ---------------- controls strip (scale, mask, order) ---------------- */
  function controls(){
    const el=document.getElementById('cmpCtl'); if(!el) return;
    if(ST.view==='table'){
      const projOpts=projectOptions(ST.dataset);
      el.innerHTML=`<div class="card pad" style="margin-bottom:12px">
        <div class="fl-lbl" style="margin-bottom:10px">Table filters</div>
        <div style="display:flex;gap:14px;flex-wrap:wrap;align-items:flex-end">
          <div><div class="fl-lbl">Similarity ≥</div><input class="cmp-f" id="fSimMin" type="number" step="0.01" min="0" max="1" value="${ST.fSimMin??''}" oninput="SNPCompare.setF('fSimMin',this.value)"></div>
          <div><div class="fl-lbl">Similarity ≤</div><input class="cmp-f" id="fSimMax" type="number" step="0.01" min="0" max="1" value="${ST.fSimMax??''}" oninput="SNPCompare.setF('fSimMax',this.value)"></div>
          <div><div class="fl-lbl">Missing% ≤</div><input class="cmp-f" id="fMissMax" type="number" step="1" min="0" max="100" value="${ST.fMissMax??''}" oninput="SNPCompare.setF('fMissMax',this.value)"></div>
          <div><div class="fl-lbl">Project</div>
            <select class="cmp-f" style="width:auto;min-width:230px;max-width:360px" onchange="SNPCompare.setF('fProj',this.value)">
              <option value="all">All projects</option>
              ${projOpts.map(p=>`<option value="${p.id}" ${ST.fProj===p.id?'selected':''}>${esc(p.label)}</option>`).join('')}
            </select></div>
          <button class="qbtn" onclick="SNPCompare.clearFilters()">Clear filters</button>
          <span style="margin-left:auto"><button class="qbtn" onclick="SNPCompare.exportCSV()">Export table CSV</button></span>
        </div>
      </div>`;
      return;
    }
    const P = hasLocal()? allPairs(ST.input) : null;
    const total = P? P.total : (ST.input? ST.input.rows.length : 0);
    const thr = minSites(total);
    const needScale = ST.view==='matrix'||ST.view==='spread';
    const needMask  = ST.view!=='table';
    el.innerHTML=`<div class="card pad" style="margin-bottom:12px">
      <div style="display:flex;gap:16px;flex-wrap:wrap;align-items:flex-end">
        ${needScale?`<div>
          <div class="fl-lbl">Color scale</div>
          <select class="cmp-f" style="width:auto;min-width:190px" onchange="SNPCompare.setScale(this.value)">
            ${Object.keys(SCALES).map(k=>`<option value="${k}" ${ST.scale===k?'selected':''}>${esc(SCALES[k].label)}</option>`).join('')}
          </select></div>`:''}
        ${ST.view==='matrix'?`<div>
          <div class="fl-lbl">Lower triangle</div>
          <select class="cmp-f" style="width:auto;min-width:190px" onchange="SNPCompare.setLower(this.value)">
            <option value="cocall" ${ST.lower==='cocall'?'selected':''}>Co-called sites (coverage)</option>
            <option value="sim" ${ST.lower==='sim'?'selected':''}>Similarity (mirror)</option>
          </select></div>
        <div>
          <div class="fl-lbl">Order</div>
          <select class="cmp-f" style="width:auto;min-width:170px" onchange="SNPCompare.setOrder(this.value)">
            <option value="cluster" ${ST.order==='cluster'?'selected':''}>Clustered (UPGMA)</option>
            <option value="sim" ${ST.order==='sim'?'selected':''}>Similarity to focal</option>
            <option value="name" ${ST.order==='name'?'selected':''}>Accession name</option>
          </select></div>`:''}
        ${needMask?`<div><div class="fl-lbl">Min co-called sites</div>
          <input class="cmp-f" type="number" min="0" step="10" value="${ST.minSites??thr}"
                 oninput="SNPCompare.setMinSites(this.value)" title="Pairs below this are drawn as unmeasured, not as low similarity">
          </div>
        <div><div class="fl-lbl">Hide lines with mean missing% &gt;</div>
          <input class="cmp-f" type="number" min="0" max="100" step="5" value="${ST.dropMiss??''}" placeholder="off"
                 oninput="SNPCompare.setDropMiss(this.value)"></div>`:''}
        ${ST.view==='mds'?`<div><div class="fl-lbl">Labels</div>
          <select class="cmp-f" style="width:auto;min-width:160px" onchange="SNPCompare.setMdsLabels(this.value)"
                  title="With many accessions, labels overlap. Auto hides ones that would collide; hover any point for its name.">
            <option value="auto" ${ST.mdsLabels==='auto'?'selected':''}>Auto (declutter)</option>
            <option value="all" ${ST.mdsLabels==='all'?'selected':''}>Show all</option>
            <option value="focal" ${ST.mdsLabels==='focal'?'selected':''}>Focal only</option>
            <option value="none" ${ST.mdsLabels==='none'?'selected':''}>None</option>
          </select></div>`:''}
        ${needScale?`<div style="flex:1 1 320px;min-width:300px">
          <div class="fl-lbl">Color domain — drag the handles</div>
          <div id="cmpLegend"></div></div>`:''}
        <div style="margin-left:auto;align-self:flex-end">
          <button class="qbtn" onclick="SNPCompare.saveImage()">Save ${ST.view==='matrix'?'matrix':ST.view==='mds'?'PCoA':'focal spread'} PNG</button>
        </div>
      </div>
      ${needScale?`<div class="mtx-note" style="margin:10px 0 0">${esc(SCALES[ST.scale].hint)}. Pairs with fewer than
        <b>${thr.toLocaleString()}</b> co-called sites of <b>${(total||0).toLocaleString()}</b> are shown as
        <span class="cmp-swatch"></span> unmeasured rather than as low similarity.</div>`:''}
    </div>`;
  }

  /* ---------------- legend: histogram with draggable endpoints ---------------- */
  /* Dragging must not rebuild the control strip, or the SVG under the pointer
     is destroyed mid-gesture. The handles are moved in place and only the
     view beneath is repainted.                                            */
  let legendDragging=false; // true while a handle is being dragged; blocks rebuilds below
  function drawLegend(sc){
    const host=document.getElementById('cmpLegend'); if(!host||!sc) return;
    // A drag in progress calls renderView() -> renderMatrix()/renderMDS() -> drawLegend()
    // on every pointermove. Rebuilding the SVG mid-gesture destroys the element that
    // holds pointer capture and the closure holding the active drag state, so the
    // drag silently stops after the first move. repaintLegend() (below, inside the
    // still-live closure) already keeps the handles/colours in sync during the drag,
    // so skip the rebuild entirely while dragging.
    if(legendDragging) return;
    const W=Math.max(260, host.clientWidth||320), H=64, pad=6, bins=44;
    const lo=sc.tmin, hi=sc.tmax, span=(hi-lo)||1;
    const counts=new Array(bins).fill(0);
    for(const t of sc.values){ let b=Math.floor((t-lo)/span*bins); if(b<0)b=0; if(b>=bins)b=bins-1; counts[b]++; }
    const mx=Math.max(1,...counts);
    const x=t=>pad+((t-lo)/span)*(W-2*pad);
    const bw=Math.max(1,(W-2*pad)/bins-1);
    const bars=counts.map((c,i)=>{
      const t0=lo+span*i/bins, h=Math.round((c/mx)*(H-26));
      return `<rect class="cmp-bin" data-t="${t0}" x="${x(t0).toFixed(1)}" y="${H-14-h}" width="${bw.toFixed(1)}" height="${h}"></rect>`;
    }).join('');
    const handle=(t,cls)=>`<g class="cmp-h" data-h="${cls}" style="cursor:ew-resize">
      <line x1="${x(t)}" y1="4" x2="${x(t)}" y2="${H-12}" stroke="#39424f" stroke-width="1.5"></line>
      <rect x="${(x(t)-5)}" y="4" width="10" height="12" rx="3" fill="#39424f"></rect></g>`;
    host.innerHTML=`<svg id="cmpLegendSvg" width="${W}" height="${H}" style="display:block;touch-action:none">
        ${bars}${handle(sc.lo,'lo')}${handle(sc.hi,'hi')}
        <text id="cmpLegLo" x="${pad}" y="${H-2}" font-size="10" fill="#6b7280">${esc(fmtInv(sc,sc.lo))}</text>
        <text id="cmpLegHi" x="${W-pad}" y="${H-2}" font-size="10" fill="#6b7280" text-anchor="end">${esc(fmtInv(sc,sc.hi))}</text>
      </svg>`;

    const svg=document.getElementById('cmpLegendSvg');
    const gs={ lo:svg.querySelector('[data-h="lo"]'), hi:svg.querySelector('[data-h="hi"]') };
    const bins$=Array.from(svg.querySelectorAll('.cmp-bin'));
    const txt={ lo:svg.querySelector('#cmpLegLo'), hi:svg.querySelector('#cmpLegHi') };
    function repaintLegend(){
      for(const r of bins$){ const t0=+r.getAttribute('data-t');
        r.setAttribute('fill', ramp(sc.pal,(t0-sc.lo)/(sc.hi-sc.lo)));
        r.setAttribute('opacity', (t0<sc.lo||t0>sc.hi)?0.35:1); }
      for(const k of ['lo','hi']){ const px=x(sc[k]);
        gs[k].querySelector('line').setAttribute('x1',px);
        gs[k].querySelector('line').setAttribute('x2',px);
        gs[k].querySelector('rect').setAttribute('x',px-5);
        txt[k].textContent=fmtInv(sc,sc[k]); }
    }
    repaintLegend();

    let drag=null;
    const val=e=>{ const r=svg.getBoundingClientRect();
      const px=Math.max(pad,Math.min(W-pad,e.clientX-r.left));
      return lo+((px-pad)/(W-2*pad))*span; };
    const move=e=>{ if(!drag) return;
      const v=val(e), cur={lo:sc.lo,hi:sc.hi}; cur[drag]=v;
      if(cur.hi-cur.lo < span*0.02) return;
      sc.lo=cur.lo; sc.hi=cur.hi; ST.dom={lo:cur.lo,hi:cur.hi};
      repaintLegend(); renderView();
    };
    svg.addEventListener('pointerdown',e=>{
      const g=e.target.closest('.cmp-h');
      drag = g? g.getAttribute('data-h') : (Math.abs(val(e)-sc.lo)<Math.abs(val(e)-sc.hi)?'lo':'hi');
      legendDragging=true;
      try{ svg.setPointerCapture(e.pointerId); }catch(_){}
      move(e);
    });
    svg.addEventListener('pointermove',move);
    svg.addEventListener('pointerup',()=>{ drag=null; legendDragging=false; drawLegend(sc); });
    svg.addEventListener('pointercancel',()=>{ drag=null; legendDragging=false; drawLegend(sc); });
    svg.addEventListener('dblclick',()=>{ ST.dom=null; paint(); });
  }
  function fmtInv(sc,t){ const v=sc.inv(t); return isFinite(v)? Math.max(0,Math.min(1,v)).toFixed(4) : '—'; }
  function median(arr){ if(!arr.length) return 0; const a=arr.slice().sort((x,y)=>x-y); return a[a.length>>1]; }

  /* ---------------- view dispatch ---------------- */
  function paint(){ controls(); renderView(); }
  function renderView(){
    const wrap=document.getElementById('cmpViewWrap'); if(!wrap) return;
    if(ST.view==='table')  return renderTable();
    if(ST.view==='matrix') return renderMatrix();
    if(ST.view==='mds')    return renderMDS();
    if(ST.view==='spread') return renderSpread();
    // Diagnostics view temporarily disabled — renderDiag() kept below for easy restore.
    // if(ST.view==='diag')   return renderDiag();
  }

  /* ---------------- table ---------------- */
  function cols(){
    const c=[{k:'rank',t:'#'},{k:'id',t:'Final_ID'}];
    if(ST.mode==='both'){
      c.push({k:'gsim',t:'Global sim'},{k:'lsim',t:'Local sim'},{k:'dsim',t:'Δ (local−global)'},
             {k:'gmiss',t:'Global miss%'},{k:'lmiss',t:'Local miss%'},{k:'lboth',t:'Co-called sites'});
    } else if(ST.mode==='local'){
      c.push({k:'lsim',t:'Similarity'},{k:'lmiss',t:'Missing%'},{k:'lboth',t:'Co-called sites'});
    } else {
      c.push({k:'gsim',t:'Similarity'},{k:'gmiss',t:'Missing%'});
    }
    c.push({k:'bio',t:'Project'},{k:'run',t:'SRA ID'},{k:'name',t:'Accession Name'});
    return c;
  }
  const fmtSim=v=>v==null?'—':v.toFixed(4);
  const fmtMiss=v=>v==null?'—':v.toFixed(2);
  const fmtD=v=>v==null?'—':(v>=0?'+':'')+v.toFixed(4);
  function cellVal(r,k){
    if(k==='lsim' && r.lsim==null && r.lboth!=null) return '<span class="cmp-unmeas" title="too few co-called sites to estimate">unmeasured</span>';
    if(k==='gsim'||k==='lsim') return fmtSim(r[k]);
    if(k==='gmiss'||k==='lmiss') return fmtMiss(r[k]);
    if(k==='dsim') return fmtD(r.dsim);
    if(k==='lboth') return r.lboth==null?'—':r.lboth.toLocaleString();
    return esc(r[k]==null?'':r[k]);
  }
  function renderTable(){
    const rows=viewRows();
    const wrap=document.getElementById('cmpViewWrap'); if(!wrap) return;
    const cs=cols();
    const sortKey = ST.mode==='local' && ST.sortKey==='sim' ? 'lsim' : (ST.mode!=='local'&&ST.sortKey==='sim'?'gsim':ST.sortKey);
    const head=cs.map(c=>{
      if(c.k==='rank') return `<th>#</th>`;
      const active = (c.k===sortKey) || (c.k==='gsim'&&ST.sortKey==='sim'&&ST.mode!=='local') || (c.k==='lsim'&&ST.sortKey==='sim'&&ST.mode==='local');
      const arrow = active ? (ST.sortDir<0?' ▾':' ▴') : '';
      return `<th onclick="SNPCompare.sortBy('${c.k}')" style="cursor:pointer;white-space:nowrap">${esc(c.t)}${arrow}</th>`;
    }).join('');
    const body=rows.map((r,i)=>{
      const focal=r.id===ST.focal;
      const tds=cs.map(c=> c.k==='rank'
        ? `<td class="num">${i+1}</td>`
        : `<td class="${(c.k.endsWith('sim')||c.k.endsWith('miss')||c.k==='dsim'||c.k==='lboth')?'num':''}">${cellVal(r,c.k)}</td>`).join('');
      return `<tr class="${focal?'cmp-focal':''}">${tds}</tr>`;
    }).join('');
    wrap.innerHTML=`<table class="cmp-table"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table>`;
    counter(`Showing <b>${rows.length.toLocaleString()}</b> of <b>${ST.allRows.length.toLocaleString()}</b> accessions for “${esc(ST.focal||'')}”.`);
    banner();
  }
  function counter(html){ const c=document.getElementById('cmpCount'); if(c) c.innerHTML=html; }
  function banner(){
    const ban=document.getElementById('cmpBanner'); if(!ban) return;
    ban.innerHTML = (ST.gdemo && ST.mode!=='local')
      ? `<div class="mtx-note">Genome-wide values shown here are <b>demonstration data</b>. Connect <span class="c-mono">ibsCompare.php</span> and set <span class="c-mono">useDemoGlobal=false</span> in snpcompare.js to load your precomputed IBS scores.</div>` : '';
  }

  /* ---------------- matrix ---------------- */
  function renderMatrix(){
    const wrap=document.getElementById('cmpViewWrap'); if(!wrap) return;
    const P=allPairs(ST.input);
    if(!P||P.n<2){ wrap.innerHTML=notice('Send a result with at least two accessions from SNPVersity.'); return; }
    const L=layout(P), order=L.order, k=order.length;
    if(!k){ wrap.innerHTML=notice('Every accession was removed by the missing-data filter. Raise the threshold.'); return; }
    const thr=minSites(P.total);

    // scale from measurable off-diagonal pairs only
    const vals=[], boths=[];
    for(let a=0;a<k;a++) for(let b=a+1;b<k;b++){ const kk=order[a]*P.n+order[b];
      if(P.both[kk]>=thr && !isNaN(P.sim[kk])){ vals.push(P.sim[kk]); boths.push(P.both[kk]); } }
    const sc=buildScale(vals, 0.5/Math.max(1,median(boths))); ST._scale=sc;

    // co-call ramp for the lower triangle (greys→teal, distinct from the similarity palette)
    let maxBoth=1; for(let a=0;a<k;a++) for(let b=a+1;b<k;b++) maxBoth=Math.max(maxBoth,P.both[order[a]*P.n+order[b]]);
    const COV=[[247,249,252],[203,224,232],[122,183,196],[38,120,142],[16,60,80]];

    const dendH = (ST.order==='cluster' && L.root)? 74 : 0;
    const labels=order.map(i=>labelFor(P.ids[i]));
    const cv=document.createElement('canvas');
    const ctx=cv.getContext('2d');
    ctx.font='12px system-ui, sans-serif';
    let lw=0; for(const s of labels) lw=Math.max(lw, ctx.measureText(s).width);
    const mar=Math.min(210, Math.max(90, Math.ceil(lw)+16));
    const avail=Math.max(320,(wrap.clientWidth||900)-mar-30);
    const cell=Math.max(7, Math.min(34, Math.floor(avail/k)));
    const W=mar+cell*k+16, H=mar+dendH+cell*k+16;
    const dpr=window.devicePixelRatio||1;
    cv.width=W*dpr; cv.height=H*dpr; cv.style.width=W+'px'; cv.style.height=H+'px';
    cv.id='cmpMatrixCv'; cv.style.cursor='crosshair';
    ctx.scale(dpr,dpr);
    ctx.fillStyle='#fff'; ctx.fillRect(0,0,W,H);

    const x0=mar, y0=mar+dendH;

    // dendrogram
    if(dendH){
      let maxh=0; (function walk(nd){ if(!nd)return; if(nd.height>maxh)maxh=nd.height; walk(nd.left); walk(nd.right); })(L.root);
      const pos=new Map(); order.forEach((g,i)=>pos.set(g,i));
      ctx.strokeStyle='#94a3b8'; ctx.lineWidth=1;
      const yFor=h=> (mar+dendH-6) - (maxh? (h/maxh)*(dendH-14) : 0);
      (function draw(nd){
        if(!nd) return null;
        if(nd.leaf!==undefined) return x0+(pos.get(nd.leaf)+0.5)*cell;
        const a=draw(nd.left), b=draw(nd.right), y=yFor(nd.height);
        ctx.beginPath(); ctx.moveTo(a,yFor(nd.left.height||0)); ctx.lineTo(a,y); ctx.lineTo(b,y);
        ctx.lineTo(b,yFor(nd.right.height||0)); ctx.stroke();
        return (a+b)/2;
      })(L.root);
    }

    // cells
    for(let a=0;a<k;a++) for(let b=0;b<k;b++){
      const gx=order[b], gy=order[a], kk=gy*P.n+gx;
      const px=x0+b*cell, py=y0+a*cell;
      if(a===b){ ctx.fillStyle=DIAG_FILL; ctx.fillRect(px,py,cell,cell); continue; }
      const measurable = P.both[kk]>=thr && !isNaN(P.sim[kk]);
      const upper = b>a;
      if(upper || ST.lower==='sim'){
        if(measurable){ ctx.fillStyle=sc.color(P.sim[kk]); ctx.fillRect(px,py,cell,cell); }
        else {
          ctx.fillStyle=MASK_FILL; ctx.fillRect(px,py,cell,cell);
          ctx.strokeStyle=MASK_LINE; ctx.lineWidth=1; ctx.beginPath();
          ctx.moveTo(px+1,py+cell-1); ctx.lineTo(px+cell-1,py+1); ctx.stroke();
        }
      } else {
        const t=P.both[kk]/maxBoth;
        ctx.fillStyle=ramp(COV,t); ctx.fillRect(px,py,cell,cell);
        if(!measurable){ ctx.strokeStyle=MASK_LINE; ctx.lineWidth=1; ctx.beginPath();
          ctx.moveTo(px+1,py+cell-1); ctx.lineTo(px+cell-1,py+1); ctx.stroke(); }
      }
    }
    // grid + focal highlight
    ctx.strokeStyle='rgba(255,255,255,.55)'; ctx.lineWidth=1;
    for(let i=0;i<=k;i++){ ctx.beginPath(); ctx.moveTo(x0+i*cell,y0); ctx.lineTo(x0+i*cell,y0+k*cell); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(x0,y0+i*cell); ctx.lineTo(x0+k*cell,y0+i*cell); ctx.stroke(); }
    const fpos=order.indexOf(P.ids.indexOf(ST.focal));
    if(fpos>=0){ ctx.strokeStyle='#b45309'; ctx.lineWidth=2;
      ctx.strokeRect(x0+fpos*cell,y0,cell,k*cell); ctx.strokeRect(x0,y0+fpos*cell,k*cell,cell); }

    // labels
    ctx.fillStyle='#39424f'; ctx.font=`${Math.min(12,Math.max(8,cell-1))}px system-ui, sans-serif`;
    ctx.textBaseline='middle';
    for(let i=0;i<k;i++){
      const lab=labels[i], isF=order[i]===P.ids.indexOf(ST.focal);
      ctx.fillStyle= isF? '#b45309':'#39424f';
      ctx.textAlign='right'; ctx.fillText(lab, x0-6, y0+i*cell+cell/2);
      ctx.save(); ctx.translate(x0+i*cell+cell/2, y0-6); ctx.rotate(-Math.PI/2);
      ctx.textAlign='left'; ctx.fillText(lab, 0, 0); ctx.restore();
    }

    wrap.innerHTML='';
    wrap.appendChild(cv);
    const key=document.createElement('div');
    key.className='mtx-note'; key.style.marginTop='10px';
    key.innerHTML=`<b>Upper triangle</b> similarity (${esc(SCALES[ST.scale].label)}). <b>Lower triangle</b> ${
      ST.lower==='cocall'? 'co-called sites — dark = well covered, pale = little shared data' : 'the same similarity, mirrored'
    }. <span class="cmp-swatch"></span> = fewer than ${thr.toLocaleString()} co-called sites, unmeasured. Click a cell to make that accession focal.`;
    wrap.appendChild(key);
    drawLegend(sc);

    // hover + click
    const tip=document.getElementById('cmpTip');
    cv.addEventListener('mousemove',e=>{
      const r=cv.getBoundingClientRect();
      const b=Math.floor((e.clientX-r.left-x0)/cell), a=Math.floor((e.clientY-r.top-y0)/cell);
      if(a<0||b<0||a>=k||b>=k){ tip.style.display='none'; return; }
      const gy=order[a], gx=order[b], kk=gy*P.n+gx;
      const meas=P.both[kk]>=thr && !isNaN(P.sim[kk]);
      tip.innerHTML=`<b>${esc(labelFor(P.ids[gy]))}</b> × <b>${esc(labelFor(P.ids[gx]))}</b><br>
        Similarity: ${ gy===gx? 'self' : (meas? P.sim[kk].toFixed(4) : `<i>unmeasured</i> (${isNaN(P.sim[kk])?'no co-called sites':P.sim[kk].toFixed(4)+' from too few sites'})`)}<br>
        Co-called: ${P.both[kk].toLocaleString()} / ${P.total.toLocaleString()} · missing ${P.miss[kk].toFixed(1)}%`;
      tip.style.display='block';
      tip.style.left=(e.clientX+14)+'px'; tip.style.top=(e.clientY+14)+'px';
    });
    cv.addEventListener('mouseleave',()=>{ tip.style.display='none'; });
    cv.addEventListener('click',e=>{
      const r=cv.getBoundingClientRect();
      const b=Math.floor((e.clientX-r.left-x0)/cell);
      if(b>=0&&b<k){ pick(P.ids[order[b]]); }
    });

    const dropped=P.n-k;
    counter(`<b>${k}</b> accessions × <b>${P.total.toLocaleString()}</b> sites${dropped?` · <b>${dropped}</b> hidden by the missing-data filter`:''} · region ${esc(ST.input.chr)}:${(+ST.input.start).toLocaleString()}–${(+ST.input.end).toLocaleString()}`);
    banner();
  }

  /* ---------------- PCoA map ---------------- */
  function renderMDS(){
    const wrap=document.getElementById('cmpViewWrap'); if(!wrap) return;
    const P=allPairs(ST.input);
    if(!P||P.n<3){ wrap.innerHTML=notice('PCoA needs at least three accessions in the region result.'); return; }
    const L=layout(P), res=pcoa(P,L.order);
    if(!res){ wrap.innerHTML=notice('Not enough measurable pairs to place the accessions.'); return; }
    const meta=metaMap(ST.dataset);
    const W=Math.max(420,Math.min(980,(wrap.clientWidth||820))), H=Math.round(W*0.66), pad=54;
    const xs=res.pts.map(p=>p.x), ys=res.pts.map(p=>p.y);
    const xr=[Math.min(...xs),Math.max(...xs)], yr=[Math.min(...ys),Math.max(...ys)];
    const sx=v=>pad+((v-xr[0])/((xr[1]-xr[0])||1))*(W-2*pad);
    const sy=v=>H-pad-((v-yr[0])/((yr[1]-yr[0])||1))*(H-2*pad);
    const thr=minSites(P.total);
    // pixel positions first, so labels can be placed with collision-avoidance
    const nodes=res.pts.map(p=>{
      const id=P.ids[p.gi], m=meta[id]||{}, isF=id===ST.focal;
      // mean missing% for this accession, used to size/soften unreliable points
      let s=0,c=0; for(let j=0;j<P.n;j++){ if(j===p.gi)continue; s+=P.miss[p.gi*P.n+j]; c++; }
      const mm=c? s/c : 0;
      return { id, isF, mm, px:sx(p.x), py:sy(p.y), r:isF?7:5,
               col:isF? '#b45309' : (m.projColor||'#3b6fd4'), label:labelFor(id) };
    });
    // shrink the type as the panel gets crowded, then declutter
    const fs = nodes.length>60 ? 9 : (nodes.length>36 ? 10 : 11);
    const plan = planLabels(nodes, fs, W, H);

    const circles=nodes.map(nd=>`<g class="cmp-pt" onclick="SNPCompare.pick('${jsStr(nd.id)}')" style="cursor:pointer">
        <title>${esc(nd.label)} — mean missing ${nd.mm.toFixed(1)}%</title>
        <circle cx="${nd.px.toFixed(1)}" cy="${nd.py.toFixed(1)}" r="${nd.r}" fill="${nd.col}"
                fill-opacity="${(1-Math.min(0.7,nd.mm/100)).toFixed(2)}" stroke="#fff" stroke-width="1.4"></circle></g>`).join('');
    const labels=plan.map(L=>L.show
        ? `<text x="${L.x.toFixed(1)}" y="${L.y.toFixed(1)}" font-size="${fs}" text-anchor="${L.anchor}"
             fill="${L.isF?'#b45309':'#39424f'}" font-weight="${L.isF?600:400}" style="pointer-events:none">${esc(L.label)}</text>`
        : '').join('');

    const hidden=plan.filter(L=>!L.show).length;
    const declutNote = ST.mdsLabels==='auto' && hidden
      ? ` <b>${hidden}</b> label${hidden>1?'s were':' was'} hidden to avoid overlap — hover any point for its name, or switch <b>Labels</b> to “Show all”.`
      : (ST.mdsLabels==='none' ? ' Labels are off — hover any point for its name.'
        : (ST.mdsLabels==='focal' ? ' Only the focal accession is labelled — hover any point for its name.' : ''));

    wrap.innerHTML=`<svg width="${W}" height="${H}" style="display:block;max-width:100%">
      <rect x="0" y="0" width="${W}" height="${H}" fill="#fff"></rect>
      <line x1="${pad-10}" y1="${H-pad}" x2="${W-pad+10}" y2="${H-pad}" stroke="#dde3ec"></line>
      <line x1="${pad}" y1="${pad-10}" x2="${pad}" y2="${H-pad+10}" stroke="#dde3ec"></line>
      <text x="${W/2}" y="${H-14}" font-size="11" fill="#6b7280" text-anchor="middle">PCoA 1 — ${(100*res.v1).toFixed(1)}% of variation</text>
      <text x="16" y="${H/2}" font-size="11" fill="#6b7280" text-anchor="middle" transform="rotate(-90 16 ${H/2})">PCoA 2 — ${(100*res.v2).toFixed(1)}%</text>
      ${circles}${labels}</svg>
      <div class="mtx-note" style="margin-top:10px">Classical multidimensional scaling of the distance matrix (1 − similarity).
      Distances between unmeasurable pairs are imputed with the panel mean, so a poorly covered accession drifts toward the centre
      rather than to an extreme — points are faded in proportion to their mean missing data. Click a point to make it focal.${declutNote}</div>`;
    counter(`<b>${res.pts.length}</b> accessions placed · masking below <b>${thr.toLocaleString()}</b> co-called sites.`);
    banner();
  }

  /* Greedy label placement for the PCoA map. Focal is always labelled and the
     most peripheral points get priority; each label tries right / left / above
     / below and is dropped if every candidate overlaps an already-placed one or
     leaves the panel. Modes: auto (declutter), all, focal, none.            */
  function planLabels(nodes, fs, W, H){
    const mode=ST.mdsLabels||'auto';
    const cw=fs*0.6, margin=1.5;
    const placed=[];
    const hit=(a,b)=> !(a.x2<b.x1-margin || a.x1>b.x2+margin || a.y2<b.y1-margin || a.y1>b.y2+margin);
    const candidates=nd=>{
      const w=Math.max(8, nd.label.length*cw);
      const mid=nd.py+fs*0.35;
      return [
        {anchor:'start',  x:nd.px+nd.r+4,   y:mid,           x1:nd.px+nd.r+4,     x2:nd.px+nd.r+4+w},
        {anchor:'end',    x:nd.px-nd.r-4,   y:mid,           x1:nd.px-nd.r-4-w,   x2:nd.px-nd.r-4},
        {anchor:'middle', x:nd.px,          y:nd.py-nd.r-4,  x1:nd.px-w/2,        x2:nd.px+w/2},
        {anchor:'middle', x:nd.px,          y:nd.py+nd.r+fs+2, x1:nd.px-w/2,      x2:nd.px+w/2},
      ].map(c=>({...c, y1:c.y-fs, y2:c.y+2}));
    };
    // focal first, then outermost-to-innermost
    const order=nodes.map((_,i)=>i).sort((a,b)=>{
      if(nodes[a].isF!==nodes[b].isF) return nodes[a].isF?-1:1;
      const da=Math.hypot(nodes[a].px-W/2, nodes[a].py-H/2), db=Math.hypot(nodes[b].px-W/2, nodes[b].py-H/2);
      return db-da;
    });
    const out=nodes.map(nd=>({label:nd.label, isF:nd.isF, show:false, x:0, y:0, anchor:'start'}));
    for(const i of order){
      const nd=nodes[i];
      if(mode==='none') continue;
      if(mode==='focal' && !nd.isF) continue;
      const cs=candidates(nd);
      let chosen=null;
      if(mode==='all'){ chosen=cs[0]; }
      else {
        for(const c of cs){
          const bx={x1:Math.min(c.x1,c.x2), x2:Math.max(c.x1,c.x2), y1:c.y1, y2:c.y2};
          if(bx.x1<2 || bx.x2>W-2 || bx.y1<2 || bx.y2>H-2) continue;   // keep inside the panel
          if(!placed.some(b=>hit(bx,b))){ chosen=c; placed.push(bx); break; }
        }
        if(!chosen && nd.isF){ chosen=cs[0]; }                          // never drop the focal label
      }
      if(chosen) out[i]={label:nd.label, isF:nd.isF, show:true, x:chosen.x, y:chosen.y, anchor:chosen.anchor};
    }
    return out;
  }

  /* ---------------- focal spread (beeswarm) ---------------- */
  function renderSpread(){
    const wrap=document.getElementById('cmpViewWrap'); if(!wrap) return;
    // Table filters are intentionally ignored by graphical views.
    const rows=ST.allRows.slice().filter(r=>r.id!==ST.focal);
    const useLocal = ST.mode==='local';
    const vals=rows.map(r=>useLocal? r.lsim : r.gsim).filter(v=>v!=null&&isFinite(v));
    if(!vals.length){ wrap.innerHTML=notice('No measurable similarities for this focal accession.'); return; }
    const boths=rows.map(r=>r.lboth).filter(v=>v!=null&&v>0);
    const sc=buildScale(vals, boths.length? 0.5/median(boths) : 1e-4); ST._scale=sc;
    const W=Math.max(460,Math.min(1000,(wrap.clientWidth||860))), pad=56;
    const lo=sc.lo, hi=sc.hi, span=(hi-lo)||1;
    const sx=v=>pad+Math.max(0,Math.min(1,(sc.tf(v)-lo)/span))*(W-2*pad);
    // beeswarm packing
    const items=rows.map(r=>({r, v: useLocal? r.lsim : r.gsim})).filter(o=>o.v!=null&&isFinite(o.v))
                    .sort((a,b)=>a.v-b.v);
    const lanes=[], R=5.2;
    for(const it of items){
      const x=sx(it.v); let lane=0;
      while(lanes[lane]!=null && x-lanes[lane] < 2*R) lane++;
      lanes[lane]=x; it.lane=lane;
    }
    const nl=Math.max(1,lanes.length), H=Math.min(520, 92+nl*11);
    const base=H-40;
    const dots=items.map(o=>{
      const cx=sx(o.v), cy=base-o.lane*11;
      return `<g onclick="SNPCompare.pick('${esc(o.r.id)}')" style="cursor:pointer">
        <title>${esc(labelFor(o.r.id))} — similarity ${o.v.toFixed(4)}${o.r.lboth!=null?`, ${o.r.lboth.toLocaleString()} co-called sites`:''}</title>
        <circle cx="${cx.toFixed(1)}" cy="${cy.toFixed(1)}" r="${R-0.6}" fill="${sc.color(o.v)}" stroke="#fff" stroke-width="1"></circle></g>`;
    }).join('');
    const unmeas=rows.filter(r=>useLocal && r.lsim==null).length;
    const ticks=[0,0.25,0.5,0.75,1].map(f=>{ const t=lo+span*f, x=pad+f*(W-2*pad);
      return `<g><line x1="${x}" y1="${base+8}" x2="${x}" y2="${base+13}" stroke="#c8d0dc"></line>
        <text x="${x}" y="${base+26}" font-size="10" fill="#6b7280" text-anchor="middle">${fmtInv(sc,t)}</text></g>`;}).join('');
    wrap.innerHTML=`<svg width="${W}" height="${H}" style="display:block;max-width:100%">
      <line x1="${pad}" y1="${base+8}" x2="${W-pad}" y2="${base+8}" stroke="#c8d0dc"></line>
      ${ticks}${dots}
      <text x="${W/2}" y="${H-4}" font-size="11" fill="#6b7280" text-anchor="middle">similarity to ${esc(labelFor(ST.focal||''))} (${esc(SCALES[ST.scale].label)} axis)</text>
    </svg>
    <div class="mtx-note" style="margin-top:10px">Every accession as one point on a single axis, so the crowded 0.8–1.0 band
    becomes the full width of the panel. ${unmeas?`<b>${unmeas}</b> accession${unmeas>1?'s are':' is'} not shown — too few co-called sites to estimate.`:''}
    Hover for values, click a point to make it focal.</div>`;
    drawLegend(sc);
    counter(`<b>${items.length}</b> accessions plotted for “${esc(ST.focal||'')}”.`);
    banner();
  }

  /* ---------------- diagnostics ---------------- */
  function renderDiag(){
    const wrap=document.getElementById('cmpViewWrap'); if(!wrap) return;
    const P=allPairs(ST.input);
    if(!P||P.n<2){ wrap.innerHTML=notice('Send a result with at least two accessions from SNPVersity.'); return; }
    const thr=minSites(P.total);
    const pts=[]; const fi=P.ids.indexOf(ST.focal);
    for(let a=0;a<P.n;a++) for(let b=a+1;b<P.n;b++){
      const kk=a*P.n+b; if(isNaN(P.sim[kk])) continue;
      pts.push({x:Math.max(1,P.both[kk]), y:P.sim[kk], f:(a===fi||b===fi), a, b});
    }
    if(!pts.length){ wrap.innerHTML=notice('No pair in this region shares a single called site.'); return; }
    const W=Math.max(460,Math.min(980,(wrap.clientWidth||860))), H=Math.round(W*0.55), pad=58;
    const xmin=1, xmax=Math.max(10,...pts.map(p=>p.x));
    const lx=v=>pad+(Math.log10(v)-Math.log10(xmin))/((Math.log10(xmax)-Math.log10(xmin))||1)*(W-2*pad);
    const ys=pts.map(p=>p.y), y0=Math.min(...ys), y1=Math.max(...ys);
    const ly=v=>H-pad-((v-y0)/((y1-y0)||1))*(H-2*pad);
    const dots=pts.map(p=>`<circle cx="${lx(p.x).toFixed(1)}" cy="${ly(p.y).toFixed(1)}" r="${p.f?4:2.8}"
        fill="${p.f?'#b45309':'#3b6fd4'}" fill-opacity="${p.f?0.95:0.45}"><title>${esc(labelFor(P.ids[p.a]))} × ${esc(labelFor(P.ids[p.b]))} — ${p.y.toFixed(4)} over ${p.x.toLocaleString()} sites</title></circle>`).join('');
    const decades=[]; for(let d=0; Math.pow(10,d)<=xmax; d++){ const v=Math.pow(10,d);
      decades.push(`<g><line x1="${lx(v)}" y1="${pad-8}" x2="${lx(v)}" y2="${H-pad}" stroke="#eef1f6"></line>
      <text x="${lx(v)}" y="${H-pad+16}" font-size="10" fill="#6b7280" text-anchor="middle">${v.toLocaleString()}</text></g>`); }
    wrap.innerHTML=`<svg width="${W}" height="${H}" style="display:block;max-width:100%">
      ${decades.join('')}
      <line x1="${lx(thr)}" y1="${pad-8}" x2="${lx(thr)}" y2="${H-pad}" stroke="#b45309" stroke-dasharray="4 3"></line>
      <text x="${lx(thr)+5}" y="${pad+2}" font-size="10" fill="#b45309">mask below ${thr.toLocaleString()} sites</text>
      <line x1="${pad}" y1="${H-pad}" x2="${W-pad}" y2="${H-pad}" stroke="#c8d0dc"></line>
      <line x1="${pad}" y1="${pad-8}" x2="${pad}" y2="${H-pad}" stroke="#c8d0dc"></line>
      <text x="${pad-8}" y="${ly(y1)+4}" font-size="10" fill="#6b7280" text-anchor="end">${y1.toFixed(3)}</text>
      <text x="${pad-8}" y="${ly(y0)+4}" font-size="10" fill="#6b7280" text-anchor="end">${y0.toFixed(3)}</text>
      ${dots}
      <text x="${W/2}" y="${H-16}" font-size="11" fill="#6b7280" text-anchor="middle">co-called sites (log scale)</text>
      <text x="16" y="${H/2}" font-size="11" fill="#6b7280" text-anchor="middle" transform="rotate(-90 16 ${H/2})">similarity</text>
    </svg>
    <div class="mtx-note" style="margin-top:10px">Each point is one pair. Similarity estimated from few sites fans out wildly at the
    left — that fan is the artefact that a 0–1 color ramp reads as real divergence. Pairs left of the dashed line are masked in the
    matrix and the spread. Orange points involve the focal accession.</div>`;
    counter(`<b>${pts.length.toLocaleString()}</b> pairs · <b>${pts.filter(p=>p.x<thr).length.toLocaleString()}</b> fall below the mask threshold.`);
    banner();
  }

  /* ---------------- orchestration ---------------- */
  async function recompute(){
    const wrap=document.getElementById('cmpViewWrap');
    ST.ran=false; updatePending();
    if(ST.mode!=='global' && !ST.input){ if(wrap) wrap.innerHTML=notice('Send a result from SNPVersity to enable region (local) comparison.'); return; }
    if(!ST.focal){ if(wrap) wrap.innerHTML=notice('Pick a focal accession to begin.'); return; }
    if(ST.mode!=='global' && ST.input && !ST.input.accs.some(a=>a.id===ST.focal)){
      if(wrap) wrap.innerHTML=notice(`Focal accession “${esc(ST.focal)}” isn’t in the region result. Pick one of the accessions you queried, or switch scope to Genome-wide.`); return;
    }
    if(ST.mode!=='global' && ST.input){
      const A=ST.input.accs.length, V=(ST.input.rows&&ST.input.rows.length)||0;
      if(ibsWork(V,A) > IBS_COST.cmpWorkBlock || V*A > IBS_COST.cmpMemBlock){
        if(wrap) wrap.innerHTML=notice('Selected data too large for SNPCompare - choose a smaller region, a lower-density SNP set, or fewer accessions, then re-run in SNPVersity.'); return;
      }
      if((ibsWork(V,A) > IBS_COST.cmpWorkWarn || V*A > IBS_COST.cmpMemWarn) && !ST.force){
        if(wrap) wrap.innerHTML=notice(
          `<b>${A}</b> accessions × <b>${V.toLocaleString()}</b> variants — computing all pairwise IBS will make the page unresponsive for roughly ${Math.max(1,cmpSeconds(V,A))}s.`+
          `<br><br><button class="btn" onclick="SNPCompare.forceRun()">Compute anyway</button>`); return;
      }
    }
    if(wrap) wrap.innerHTML=`<div class="loading"><div class="spinner"></div><div>Computing similarities…</div></div>`;
    try{ await buildRows(); paint(); ST.ran=true; ST.ranFocal=ST.focal; ST.ranMode=ST.mode; updatePending(); }
    catch(err){ if(wrap) wrap.innerHTML=notice('Could not load similarities: '+esc(err.message||err)); }
  }
  function forceRun(){ ST.force=true; recompute(); }
  /* idle placeholder shown before anything is run (or after options change) */
  function idleMessage(){
    if(!hasLocal() && !globalAvailable())
      return 'Send a multi-accession result from SNPVersity to enable comparison.';
    const opts = hasLocal() ? '<b>Table</b>, <b>Matrix</b>, <b>PCoA</b> or <b>Focal spread</b>' : '<b>Table</b> or <b>Focal spread</b>';
    return `Set the focal accession and scope above, then press ${opts} to run the comparison.`;
  }
  function showIdle(){
    const ctl=document.getElementById('cmpCtl');
    if(ST.view==='table') controls(); else if(ctl) ctl.innerHTML='';
    const cnt=document.getElementById('cmpCount'); if(cnt) cnt.innerHTML='';
    const ban=document.getElementById('cmpBanner'); if(ban) ban.innerHTML='';
    const wrap=document.getElementById('cmpViewWrap'); if(wrap) wrap.innerHTML=notice(idleMessage());
    updatePending();
  }
  /* subtle hint when the field/scope no longer match the displayed result */
  function updatePending(){
    const el=document.getElementById('cmpPending'); if(!el) return;
    const stale = ST.ran && (ST.focal!==ST.ranFocal || ST.mode!==ST.ranMode);
    el.style.display = stale ? 'inline' : 'none';
  }
  let _rz=null;
  function bindResize(){
    if(window._cmpResizeBound) return; window._cmpResizeBound=true;
    window.addEventListener('resize',()=>{ clearTimeout(_rz); _rz=setTimeout(()=>{
      if(document.getElementById('cmpViewWrap') && ST.view!=='table') paint(); },180); });
  }

  /* ---------------- public setters ---------------- */
  /* focal field typing: update state + chip highlight only, never run */
  function syncFocal(v){ ST.focal=(v||'').trim()||null; ST._layout=null; refreshAccList(); updatePending(); }
  /* clicking a chip in the imported-accession list: just populate the field */
  function pickFocal(id){ ST.focal=id; ST._layout=null;
    const el=document.getElementById('cmpFocal'); if(el) el.value=id; refreshAccList(); updatePending(); }
  /* run the currently selected view (Enter in the focal field, or programmatic) */
  function runCurrent(){ let v=ST.view||'table'; if((v==='matrix'||v==='mds')&&!hasLocal()) v='table'; setView(v); }
  /* legacy entry point — commit the field, then run the current view */
  function setFocalFromInput(){ const el=document.getElementById('cmpFocal'); if(el){ const v=el.value.trim(); if(v) ST.focal=v; } runCurrent(); }
  /* clicking inside a rendered plot refocuses AND re-runs that view */
  function pick(id){ ST.focal=id; ST._layout=null;
    const el=document.getElementById('cmpFocal'); if(el) el.value=id; refreshAccList(); recompute(); }
  /* scope is an option: switching it does not run — it returns to idle */
  function setMode(m){
    if((m==='local'||m==='both') && !ST.input) return;
    if((m==='global'||m==='both') && !globalAvailable()) return;
    ST.mode=m; ST.ran=false;
    document.getElementById('page').innerHTML=shell(); showIdle(); }
  /* the view buttons ARE the submit action: commit the field, then run */
  function setView(v){
    if(v!=='table' && v!=='spread' && !hasLocal()) return;
    const el=document.getElementById('cmpFocal'); if(el){ const val=el.value.trim(); if(val) ST.focal=val; }
    ST.view=v;
    if((v==='matrix'||v==='mds') && ST.mode==='global' && hasLocal()) ST.mode='local';
    document.getElementById('page').innerHTML=shell(); recompute(); }
  /* table filters live-update only the table; graphical views always use the full result */
  function setF(k,v){ ST[k]= (v===''||v==null)?null:(k==='fProj'?v:parseFloat(v)); if(k==='fProj')ST.fProj=v; if(ST.ran && ST.view==='table') renderTable(); }
  function setScale(v){ ST.scale=v; ST.dom=null; if(ST.ran) paint(); }
  function setLower(v){ ST.lower=v; if(ST.ran) paint(); }
  function setOrder(v){ ST.order=v; ST._layout=null; if(ST.ran) paint(); }
  function setMinSites(v){ ST.minSites=(v===''||v==null)?null:Math.max(0,parseInt(v,10)||0); ST._layout=null; ST.dom=null; if(ST.ran) paint(); }
  function setDropMiss(v){ ST.dropMiss=(v===''||v==null)?null:parseFloat(v); ST._layout=null; if(ST.ran) paint(); }
  function setMdsLabels(v){ ST.mdsLabels=v; if(ST.ran) paint(); }
  function clearFilters(){
    ST.fSimMin=ST.fSimMax=ST.fMissMax=null; ST.fProj='all';
    if(ST.ran) paint(); else showIdle();
  }
  function sortBy(k){ if(k==='rank'||!ST.ran)return;
    if(ST.sortKey===k){ ST.sortDir*=-1; } else { ST.sortKey=k; ST.sortDir=(k==='id'||k==='name'||k==='run'||k==='bio')?1:-1; }
    renderTable(); }
  function toTree(){ if(!ST.input)return; S.treeInput=ST.input; go('snptree'); }
  function toMatrix(){ if(!ST.input)return; S.matrixInput=ST.input; go('snpmatrix'); }

  function imageFileName(){
    const view={matrix:'distance_matrix',mds:'pcoa',spread:'focal_spread'}[ST.view]||'view';
    const region=ST.input ? `${ST.input.chr||'region'}_${ST.input.start||0}_${ST.input.end||0}` : 'genomewide';
    const focal=(ST.focal||'focal').replace(/[^\w.-]+/g,'_');
    return `snpcompare_${view}_${region}_${focal}.png`.replace(/[^\w.-]+/g,'_');
  }
  function downloadBlob(name, blob){
    const u=URL.createObjectURL(blob), a=document.createElement('a');
    a.href=u; a.download=name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(()=>URL.revokeObjectURL(u),1500);
  }
  function saveCanvasPNG(canvas, name){
    if(canvas.toBlob){
      canvas.toBlob(blob=>{ if(blob) downloadBlob(name,blob); },'image/png');
      return;
    }
    const a=document.createElement('a'); a.href=canvas.toDataURL('image/png'); a.download=name;
    document.body.appendChild(a); a.click(); a.remove();
  }
  function saveSvgPNG(svg, name){
    const w=Math.max(1,parseFloat(svg.getAttribute('width')) || (svg.viewBox&&svg.viewBox.baseVal.width) || svg.getBoundingClientRect().width);
    const h=Math.max(1,parseFloat(svg.getAttribute('height')) || (svg.viewBox&&svg.viewBox.baseVal.height) || svg.getBoundingClientRect().height);
    const clone=svg.cloneNode(true);
    clone.setAttribute('xmlns','http://www.w3.org/2000/svg');
    clone.setAttribute('width',w); clone.setAttribute('height',h);
    if(!clone.getAttribute('viewBox')) clone.setAttribute('viewBox',`0 0 ${w} ${h}`);
    const source=new XMLSerializer().serializeToString(clone);
    const u=URL.createObjectURL(new Blob([source],{type:'image/svg+xml;charset=utf-8'}));
    const img=new Image();
    img.onload=()=>{
      const scale=Math.max(2,window.devicePixelRatio||1);
      const canvas=document.createElement('canvas');
      canvas.width=Math.ceil(w*scale); canvas.height=Math.ceil(h*scale);
      const ctx=canvas.getContext('2d');
      ctx.scale(scale,scale); ctx.fillStyle='#fff'; ctx.fillRect(0,0,w,h); ctx.drawImage(img,0,0,w,h);
      URL.revokeObjectURL(u); saveCanvasPNG(canvas,name);
    };
    img.onerror=()=>{ URL.revokeObjectURL(u); };
    img.src=u;
  }
  function saveImage(){
    if(!ST.ran || ST.view==='table') return;
    const wrap=document.getElementById('cmpViewWrap'); if(!wrap) return;
    const name=imageFileName();
    if(ST.view==='matrix'){
      const canvas=document.getElementById('cmpMatrixCv'); if(canvas) saveCanvasPNG(canvas,name);
      return;
    }
    const svg=wrap.querySelector('svg'); if(svg) saveSvgPNG(svg,name);
  }

  function download(name, text){
    const blob=new Blob([text],{type:'text/csv'});
    const u=URL.createObjectURL(blob), a=document.createElement('a');
    a.href=u; a.download=name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(()=>URL.revokeObjectURL(u),1500);
  }
  function exportCSV(){
    const cs=cols(), rows=viewRows();
    const head=cs.map(c=>c.k==='rank'?'rank':c.t).join(',');
    const lines=rows.map((r,i)=>cs.map(c=>{
      if(c.k==='rank') return i+1;
      let v = c.k==='dsim'?r.dsim : r[c.k];
      if(v==null) return '';
      if(typeof v==='string' && /[",\n]/.test(v)) return '"'+v.replace(/"/g,'""')+'"';
      return v;
    }).join(','));
    download(`snpcompare_${ST.mode}_${(ST.focal||'focal').replace(/[^\w.-]/g,'_')}.csv`, [head].concat(lines).join('\n'));
  }
  /* long-format matrix export: one row per pair, with the mask flag kept explicit */
  function exportMatrixCSV(){
    const P=allPairs(ST.input); if(!P) return;
    const thr=minSites(P.total);
    const out=['a_id,b_id,a_name,b_name,similarity,distance,co_called_sites,total_sites,missing_pct,measurable'];
    for(let a=0;a<P.n;a++) for(let b=a+1;b<P.n;b++){
      const kk=a*P.n+b, s=P.sim[kk], ok=P.both[kk]>=thr && !isNaN(s);
      out.push([P.ids[a],P.ids[b],`"${labelFor(P.ids[a])}"`,`"${labelFor(P.ids[b])}"`,
        isNaN(s)?'':s.toFixed(6), isNaN(s)?'':(1-s).toFixed(6),
        P.both[kk], P.total, P.miss[kk].toFixed(3), ok?1:0].join(','));
    }
    download(`snpcompare_pairs_${(ST.input.chr||'region')}_${ST.input.start}_${ST.input.end}.csv`, out.join('\n'));
  }

  /* ---------------- helpers ---------------- */
  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  /* safe to drop inside a single-quoted inline JS string (e.g. onclick="...('X')") */
  function jsStr(s){ return String(s==null?'':s).replace(/\\/g,'\\\\').replace(/'/g,"\\'").replace(/</g,'\\u003c'); }
  function notice(html){ return `<div style="text-align:center;padding:26px;color:var(--muted);max-width:640px;margin:0 auto">${html}</div>`; }
  function injectCSS(){
    if(document.getElementById('snpcompare-css'))return;
    const s=document.createElement('style'); s.id='snpcompare-css';
    s.textContent=`
      .fl-lbl{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:7px}
      .cmp-f{border:1px solid var(--line);border-radius:9px;padding:8px 10px;font-size:13px;width:110px}
      .mtx-note{font-size:12px;color:var(--muted);background:var(--blue-50,#eef4ff);border:1px solid var(--line);border-radius:8px;padding:9px 12px;margin-bottom:12px}
      .mtx-note b{color:var(--ink)}
      .btn.solid{background:var(--green-50,#e9f4ec);border-color:#bfe5cb;color:var(--green-600,#1f8a4c)}
      .cmp-table{border-collapse:collapse;width:100%;font-size:13px}
      .cmp-table thead th{position:sticky;top:0;background:#eef1f6;border-bottom:1px solid var(--line);
        text-align:left;padding:9px 12px;font-weight:600;color:var(--ink);z-index:1}
      .cmp-table td{padding:8px 12px;border-bottom:1px solid #eef1f5;white-space:nowrap}
      .cmp-table td.num{text-align:right;font-family:var(--mono)}
      .cmp-table tbody tr:nth-child(even){background:#fafbfd}
      .cmp-table tbody tr:hover{background:#eef4ff}
      .cmp-table tr.cmp-focal{background:#fff6cf !important;font-weight:600}
      .cmp-unmeas{color:#8a93a3;font-style:italic}
      .cmp-swatch{display:inline-block;width:11px;height:11px;border-radius:2px;background:#e6e9ef;
        border:1px solid #b6bdc9;vertical-align:-1px;margin:0 2px}
      .cmp-tip{position:fixed;z-index:60;background:#232a33;color:#fff;font-size:11.5px;line-height:1.45;
        padding:7px 9px;border-radius:7px;pointer-events:none;max-width:320px;box-shadow:0 6px 18px rgba(0,0,0,.22)}
      .cmp-pt text{pointer-events:none}
      /* imported-accession list (mirrors SNPVersity's "Selected" panel, kept compact) */
      .cmp-acc-wrap{margin-top:14px}
      .cmp-acc-box{margin-top:7px;border:1px solid var(--line);border-radius:9px;background:#fafbfd;padding:8px;
        max-height:104px;overflow:auto;display:flex;flex-wrap:wrap;gap:6px;align-content:flex-start}
      .cmp-acc-chip{display:inline-flex;align-items:center;gap:6px;background:#fff;border:1px solid var(--line);
        border-radius:7px;padding:4px 9px;font-size:12px;line-height:1.1;color:var(--ink);cursor:pointer;white-space:nowrap}
      .cmp-acc-chip:hover{border-color:#9db4d6;background:#eef4ff}
      .cmp-acc-chip.on{background:#fff6cf;border-color:#e6cf6a;font-weight:600}
      .cmp-acc-dot{width:8px;height:8px;border-radius:2px;flex:0 0 auto}`;
    document.head.appendChild(s);
  }

  if(typeof SNPTools!=='undefined') SNPTools.register('snpcompare', { render });

  return { render, setFocalFromInput, syncFocal, pickFocal, runCurrent, pick, setMode, setView,
           setF, setScale, setLower, setOrder, setMinSites, setDropMiss, setMdsLabels,
           clearFilters, sortBy, toTree, toMatrix, exportCSV, exportMatrixCSV, saveImage, forceRun,
           // testing / debugging
           dose, localCompute, allPairs, buildScale, cluster, pcoa, planLabels, getGlobal, _CFG:CFG, _ST:ST };
})();
if(typeof window!=='undefined') window.SNPCompare = SNPCompare;
