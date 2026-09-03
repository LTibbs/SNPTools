/* =====================================================================
 *  snptree.js — SNPTree: local phylogeny & similarity for a SNPVersity
 *  result set. Builds an identity-by-state (IBS) distance tree from the
 *  genotype matrix already in memory (no re-parsing) and renders it with
 *  the accession metadata (bioproject color, founder, full ID).
 *
 *  Method lineage: UPGMA / Neighbour-Joining, Newick / MEGA / PHYLIP
 *  outputs — the same vocabulary as VCF2PopTree (Kumar lab), reimplemented
 *  here to run directly on the in-memory matrix and our metadata.
 *
 *  Hand-off: SNPVersity sets  S.treeInput = {rows, accs, chr, start, end,
 *  dataset, datasetName, vcfUrl}  then calls go('snptree').
 *
 *  Depends on core.js (SNPTools, S, go, ICONS) and data.js (Data).
 * ===================================================================== */
const SNPTree = (function () {

  /* ------------------------------------------------------------------ *
   *  1. IBS distance from the genotype matrix
   * ------------------------------------------------------------------ */
  // gts entries are 1-byte codes from Data.parseVcf: 0/1/2 = dosage, 3 = missing.
  function dosage(gt){ return (gt == null || gt === 3) ? null : gt; }

  // rows: [{gts:[...n]}], returns {M:n×n mean allele-difference, C:n×n shared-site counts, sites, informative}
  function ibsMatrix(rows, n){
    const M = Array.from({length:n}, () => new Float64Array(n));
    const C = Array.from({length:n}, () => new Float64Array(n));
    let informative = 0;
    const d = new Array(n);
    for (const r of rows){
      const g = r.gts; let seen0=0, seen2=0, seenH=0, nonmiss=0;
      for (let k=0;k<n;k++){ const v=dosage(g[k]); d[k]=v;
        if (v!=null){ nonmiss++; if(v===0)seen0++; else if(v===2)seen2++; else seenH++; } }
      // count as informative if it actually varies among the selected samples
      if (nonmiss>1 && ((seen0>0?1:0)+(seen2>0?1:0)+(seenH>0?1:0))>1) informative++;
      for (let i=0;i<n;i++){ const di=d[i]; if (di==null) continue;
        for (let j=i+1;j<n;j++){ const dj=d[j]; if (dj==null) continue;
          M[i][j] += Math.abs(di-dj)/2; C[i][j] += 1;
        }
      }
    }
    for (let i=0;i<n;i++) for (let j=i+1;j<n;j++){
      const c=C[i][j]; const dist = c? M[i][j]/c : 1;  // no shared sites -> max distance
      M[i][j]=M[j][i]=dist; C[j][i]=C[i][j];
    }
    return {M, C, sites: rows.length, informative};
  }

  /* ------------------------------------------------------------------ *
   *  2. Tree construction — UPGMA (rooted) & Neighbour-Joining
   * ------------------------------------------------------------------ */
  const leaf = (name,i) => ({name, idx:i, children:null, length:0, height:0, size:1});

  function upgma(M0, labels){
    const n = labels.length;
    const nodes = labels.map(leaf);
    const D = M0.map(r => Array.from(r));
    while (nodes.length > 1){
      // nearest pair
      let mi=0,mj=1,min=Infinity;
      for (let a=0;a<nodes.length;a++) for (let b=a+1;b<nodes.length;b++)
        if (D[a][b] < min){ min=D[a][b]; mi=a; mj=b; }
      const A=nodes[mi], B=nodes[mj], h=min/2;
      A.length=Math.max(0,h-A.height); B.length=Math.max(0,h-B.height);
      const U={name:null, children:[A,B], length:0, height:h, size:A.size+B.size};
      // merged distances (weighted average = UPGMA)
      const newRow=[];
      for (let k=0;k<nodes.length;k++){ if(k===mi||k===mj){newRow.push(0);continue;}
        newRow.push((A.size*D[mi][k]+B.size*D[mj][k])/(A.size+B.size)); }
      // rebuild node/matrix lists without mi,mj, then append U
      const keep=[]; for(let k=0;k<nodes.length;k++) if(k!==mi&&k!==mj) keep.push(k);
      const nn = keep.map(k=>nodes[k]); nn.push(U);
      const nD = keep.map(k=>{ const row=keep.map(k2=>D[k][k2]); row.push(newRow[k]); return row; });
      const lastRow = keep.map(k=>newRow[k]); lastRow.push(0); nD.push(lastRow);
      nodes.length=0; Array.prototype.push.apply(nodes,nn);
      D.length=0; Array.prototype.push.apply(D,nD);
    }
    return {root:nodes[0], unrooted:false};
  }

  function nj(M0, labels){
    const n = labels.length;
    let nodes = labels.map(leaf);
    let D = M0.map(r => Array.from(r));
    while (nodes.length > 2){
      const m = nodes.length;
      const r = new Array(m).fill(0);
      for (let i=0;i<m;i++){ let s=0; for(let j=0;j<m;j++) s+=D[i][j]; r[i]=s; }
      let mi=0,mj=1,min=Infinity;
      for (let i=0;i<m;i++) for (let j=i+1;j<m;j++){
        const q=(m-2)*D[i][j]-r[i]-r[j];
        if (q<min){ min=q; mi=i; mj=j; }
      }
      const A=nodes[mi], B=nodes[mj];
      const dij=D[mi][mj];
      let la=0.5*dij + (r[mi]-r[mj])/(2*(m-2));
      let lb=dij-la;
      A.length=Math.max(0,la); B.length=Math.max(0,lb);
      const U={name:null, children:[A,B], length:0, height:0, size:A.size+B.size};
      const keep=[]; for(let k=0;k<m;k++) if(k!==mi&&k!==mj) keep.push(k);
      const newDist = keep.map(k=>0.5*(D[mi][k]+D[mj][k]-dij));
      const nn = keep.map(k=>nodes[k]); nn.push(U);
      const nD = keep.map((k,ki)=>{ const row=keep.map(k2=>D[k][k2]); row.push(newDist[ki]); return row; });
      const lastRow=newDist.slice(); lastRow.push(0); nD.push(lastRow);
      nodes=nn; D=nD;
    }
    // two remain -> join under a root (arbitrary/displayed root; tree is unrooted)
    const A=nodes[0], B=nodes[1], d=D[0][1];
    A.length=Math.max(0,d/2); B.length=Math.max(0,d/2);
    return {root:{name:null, children:[A,B], length:0, height:0, size:A.size+B.size}, unrooted:true};
  }

  /* ------------------------------------------------------------------ *
   *  3. Text outputs
   * ------------------------------------------------------------------ */
  function toNewick(node){
    const rec = nd => nd.children
      ? '(' + nd.children.map(rec).join(',') + ')' + (nd.length? ':'+nd.length.toFixed(5):'')
      : safeName(nd.name) + ':' + (nd.length||0).toFixed(5);
    return rec(node) + ';';
  }
  const safeName = s => /[\s(),:;']/.test(s) ? "'" + String(s).replace(/'/g,"''") + "'" : s;

  function phylip(M, labels){
    let out = ' ' + labels.length + '\n';
    for (let i=0;i<labels.length;i++){
      let name = labels[i].slice(0,10).padEnd(10,' ');
      out += name + '  ' + M[i].map(v=>v.toFixed(5)).join(' ') + '\n';
    }
    return out;
  }
  // MEGA lower-triangular pairwise distances
  function mega(M, labels){
    let out = '#mega\n!Title: SNPTree local IBS distances;\n!Format DataType=Distance DataFormat=LowerLeft NTaxa=' + labels.length + ';\n\n';
    labels.forEach((l,i)=> out += '['+(i+1)+'] #' + l + '\n');
    out += '\n';
    for (let i=0;i<labels.length;i++){
      let row=[]; for (let j=0;j<i;j++) row.push(M[i][j].toFixed(5));
      out += '['+(i+1)+'] ' + row.join(' ') + '\n';
    }
    return out;
  }

  /* ------------------------------------------------------------------ *
   *  4. Layout
   * ------------------------------------------------------------------ */
  function collectLeaves(root){ const L=[]; (function rec(n){ n.children? n.children.forEach(rec) : L.push(n); })(root); return L; }
  function assignDepth(root){ (function rec(n,d){ n.depth=d+(n.length||0); (n.children||[]).forEach(c=>rec(c,n.depth)); })(root,0); }
  function maxDepth(root){ let m=0; (function rec(n){ if(n.depth>m)m=n.depth; (n.children||[]).forEach(rec); })(root); return m; }

  function layoutRect(root, W, H, pad){
    assignDepth(root);
    const leaves=collectLeaves(root);
    const md=maxDepth(root)||1;
    const stepY=(H-2*pad)/Math.max(1,leaves.length-1);
    leaves.forEach((l,i)=> l.y=pad+i*stepY);
    (function setY(n){ if(n.children){ n.children.forEach(setY);
      n.y=(n.children[0].y+n.children[n.children.length-1].y)/2; } })(root);
    const sx=(W-2*pad)/md;
    (function setX(n){ n.x=pad+n.depth*sx; (n.children||[]).forEach(setX); })(root);
    return {leaves, md};
  }

  function layoutRadial(root, W, H, pad){
    assignDepth(root);
    const leaves=collectLeaves(root);
    const md=maxDepth(root)||1;
    const cx=W/2, cy=H/2, R=Math.min(W,H)/2-pad;
    leaves.forEach((l,i)=> l.ang=(i/leaves.length)*Math.PI*2);
    (function setA(n){ if(n.children){ n.children.forEach(setA);
      n.ang=(n.children[0].ang+n.children[n.children.length-1].ang)/2; } })(root);
    (function setP(n){ const rad=(n.depth/md)*R; n.x=cx+Math.cos(n.ang)*rad; n.y=cy+Math.sin(n.ang)*rad;
      (n.children||[]).forEach(setP); })(root);
    return {leaves, cx, cy, R};
  }

  /* ------------------------------------------------------------------ *
   *  5. State + tool render
   * ------------------------------------------------------------------ */
  const ST = { method:'upgma', layout:'rectangular', format:'newick', colorBy:'proj', markers:'color',
               view:'tree', input:null, built:null, force:false, _pi:null };

  function metaFor(name){ return ST.metaMap ? ST.metaMap[name] : null; }

  function render(){
    const page=document.getElementById('page');
    injectCSS();
    const crumb=document.getElementById('crumbTool'); if(crumb) crumb.innerHTML='<b>SNPTree</b>';
    const inp = S.treeInput || ST.input;
    if (inp !== ST.input) ST.force = false;
    ST.input = inp;
    if (!inp || !inp.rows || !inp.rows.length){
      page.innerHTML = emptyState();
      return;
    }
    // metadata lookup by accession id (= tree label)
    ST.metaMap = {}; inp.accs.forEach(a=>{ ST.metaMap[a.id]=a; });
    ST._pi = null;   // recompute project->shape mapping for this dataset
    page.innerHTML = shell(inp);
    draw();
  }

  function emptyState(){
    return `<section class="sec"><div class="bar"></div><div>
      <h1>SNPTree</h1><p>Local phylogeny &amp; similarity from a SNPVersity result.</p></div></section>
      <div class="card pad fade" style="text-align:center">
        <div style="max-width:560px;margin:0 auto">
          <h3 style="font-family:var(--disp);margin:0 0 8px">No data yet</h3>
          <p style="color:var(--muted)">Run a query in <b>SNPVersity</b>, then use
          <b>“Send data to SNPTree”</b> in the results view. SNPTree builds an identity-by-state
          tree from the variants in that exact interval, across the accessions you selected.</p>
          <button class="btn solid" onclick="go('snpversity')">Go to SNPVersity</button>
        </div>
      </div>`;
  }

  function shell(inp){
    const n=inp.accs.length;
    const region=`${inp.chr}:${(+inp.start).toLocaleString()}–${(+inp.end).toLocaleString()}`;
    return `
    <section class="sec"><div class="bar"></div><div style="width:100%">
      <h1>SNPTree <span style="font-weight:400;color:var(--faint);font-size:14px">· powered by the VCF2PopTree method</span></h1>
      <p>Identity-by-state phylogeny for <span class="c-mono" style="color:var(--blue-600)">${region}</span> · ${inp.datasetName||inp.dataset}</p>
    </div></section>

    <div class="card pad" style="margin-bottom:16px">
      <div style="display:flex;gap:26px;flex-wrap:wrap;align-items:flex-end">
        <div><div class="fl-lbl">Construct tree</div>
          ${radio('method','upgma','UPGMA (rooted)')}${radio('method','nj','Neighbour-Joining (unrooted)')}</div>
        <div><div class="fl-lbl">Layout</div>
          ${radio('layout','rectangular','Rectangular')}${radio('layout','radial','Radial')}</div>
        <div><div class="fl-lbl">Color tips by</div>
          ${radio('colorBy','proj','Bioproject')}${radio('colorBy','none','None')}</div>
        <div><div class="fl-lbl">Markers</div>
          ${radio('markers','color','Color')}${radio('markers','shape','Shapes')}</div>
        <div><div class="fl-lbl">Text output</div>
          ${radio('format','newick','Newick')}${radio('format','mega','Pairwise (MEGA)')}${radio('format','phylip','PHYLIP')}</div>
        <div style="margin-left:auto;display:flex;gap:8px;align-items:flex-end">
          <button class="btn solid" onclick="SNPTree.draw()">Draw tree</button>
        </div>
      </div>
      <div style="display:flex;gap:8px;flex-wrap:wrap;margin-top:14px">
        <button class="qbtn" onclick="SNPTree.toCompare('table')">SNPCompare (Table view) →</button>
        <button class="qbtn" onclick="SNPTree.toCompare('matrix')">SNPCompare (Distance matrix) →</button>
        <span style="margin-left:auto"></span>
        <button class="qbtn" onclick="SNPTree.downloadSVG()">Download SVG</button>
        <button class="qbtn" onclick="SNPTree.downloadPNG()">Download PNG</button>
        <button class="qbtn" onclick="SNPTree.downloadPDF()">Download PDF</button>
      </div>
    </div>

    <div id="treeMeta" class="meta-strip"></div>
    <div class="card pad fade" id="treeStage" style="overflow:auto"></div>
    <div class="card pad" style="margin-top:16px">
      <div style="display:flex;align-items:center;gap:10px;margin-bottom:8px">
        <b id="fmtTitle">Newick</b>
        <button class="qbtn" style="margin-left:auto" onclick="SNPTree.copyText()">Copy</button>
        <button class="qbtn" onclick="SNPTree.downloadText()">Download</button>
      </div>
      <textarea id="fmtOut" readonly style="width:100%;height:150px;font-family:var(--mono);font-size:12px;
        border:1px solid var(--line);border-radius:9px;padding:10px;resize:vertical;white-space:pre;overflow:auto"></textarea>
    </div>`;
  }

  function radio(group,val,label){
    const on = ST[group]===val;
    return `<label style="display:inline-flex;align-items:center;gap:6px;margin-right:14px;cursor:pointer;font-size:13px">
      <input type="radio" name="snptree_${group}" ${on?'checked':''}
        onchange="SNPTree.setOpt('${group}','${val}')"> ${label}</label>`;
  }

  /* ----- option setters ----- */
  function setOpt(g,v){ ST[g]=v;
    if (g==='format') refreshText();     // text output can update without rebuilding
    else draw();                          // topology/layout/color need a redraw
  }
  function toCompare(view){
    if(ST.input){
      // Reuse the exact SNPVersity result that produced this tree.
      S.compareInput=ST.input;
      S.compareRequest={
        view:view==='matrix'?'matrix':'table',
        mode:'local',
        autorun:true
      };
    }
    go('snpcompare');
  }

  /* ------------------------------------------------------------------ *
   *  6. Build + draw
   * ------------------------------------------------------------------ */
  function build(){
    const inp=ST.input;
    const accs=inp.accs, n=accs.length;
    const labels=accs.map(a=>a.id);
    const {M,C,sites,informative}=ibsMatrix(inp.rows, n);
    const tree = ST.method==='nj' ? nj(M,labels) : upgma(M,labels);
    // mean shared sites across pairs
    let sp=0,cnt=0; for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){ sp+=C[i][j]; cnt++; }
    ST.built = {M, C, tree, labels, sites, informative, meanShared: cnt? sp/cnt:0, n};
    return ST.built;
  }

  function draw(){
    const inp=ST.input; if(!inp) return;
    const n=inp.accs.length;
    const stage=document.getElementById('treeStage');
    const meta=document.getElementById('treeMeta');
    if (n<3){ stage.innerHTML=notice('Need at least 3 accessions to build a tree. Select more in SNPVersity.'); meta.innerHTML=''; return; }
    const V = (inp.rows && inp.rows.length) || 0;
    if (ibsWork(V,n) > IBS_COST.workBlock){
      stage.innerHTML = notice('Selected data too large for SNPTree - choose a smaller region, a lower-density SNP set, or fewer accessions, then re-run in SNPVersity.');
      meta.innerHTML=''; return;
    }
    if (ibsWork(V,n) > IBS_COST.workWarn && !ST.force){
      stage.innerHTML = notice(
        `<b>${n}</b> accessions × <b>${V.toLocaleString()}</b> variants. Computing the distance matrix and drawing the tree will make the page unresponsive for roughly ${ibsSeconds(V,n)}s.`+
        `<br><br><button class="btn" onclick="SNPTree.forceBuild()">Build anyway</button>`);
      meta.innerHTML=''; return;
    }
    stage.innerHTML = `<div class="loading"><div class="spinner"></div><div>Computing IBS distances &amp; building ${ST.method.toUpperCase()} tree…</div></div>`;
    // let the spinner paint before the heavy sync work
    setTimeout(()=>{
      const b = build();
      meta.innerHTML = metaStrip(inp, b);
      stage.innerHTML = (ST.layout==='radial'? radialSVG(b) : rectSVG(b));
      refreshText();
    }, 20);
  }
  function forceBuild(){ ST.force=true; draw(); }

  function metaStrip(inp, b){
    const projs = uniqueProjects(inp.accs);
    const legend = (ST.colorBy==='none' && ST.markers!=='shape') ? '' :
      `<div class="legend-wrap">${projs.map(p=>{
        const col = ST.colorBy==='none' ? '#5b6b83' : p.color;
        const glyph = ST.markers==='shape'
          ? `<svg width="15" height="15" viewBox="0 0 15 15" style="vertical-align:-2px">${shapeMarkup(shapeForProjId(p.proj),7.5,7.5,5.2,col)}</svg>`
          : `<span class="sw" style="background:${col}"></span>`;
        return `<span class="li">${glyph}${esc(p.label)}</span>`;
      }).join('')}</div>`;
    return `<div class="meta-row">
      <span><b>${b.n}</b> accessions</span>
      <span><b>${b.sites.toLocaleString()}</b> variant sites</span>
      <span><b>${b.informative.toLocaleString()}</b> informative</span>
      <span>~<b>${Math.round(b.meanShared).toLocaleString()}</b> sites/pair</span>
      <span>metric: <b>IBS allele distance</b></span>
      <span>method: <b>${ST.method==='nj'?'Neighbour-Joining':'UPGMA'}</b>${b.tree.unrooted?' · unrooted':''}</span>
    </div>${legend}`;
  }
  function uniqueProjects(accs){
    const seen={}, out=[];
    accs.forEach(a=>{ const pid=a.proj||a.projTitle||'—', key=a.projTitle||a.proj||'—';
      if(!seen[key]){ seen[key]=1; out.push({label:key, color:a.projColor||'#5b6b83', proj:pid}); } });
    return out.slice(0,24);
  }
  function tipColor(a){
    if (!a || ST.colorBy==='none') return '#5b6b83';
    return a.projColor || '#5b6b83';
  }

  /* ----- project markers: color (circle) or distinct shapes ----- */
  const SHAPES=['circle','square','triangleUp','diamond','triangleDown','star','plus','hexagon','pentagon','ex'];
  function projIndex(){
    if (ST._pi) return ST._pi;
    const m={}; let i=0;
    (ST.input?ST.input.accs:[]).forEach(a=>{ const k=a.proj||a.projTitle||'—'; if(m[k]==null) m[k]=i++; });
    ST._pi=m; return m;
  }
  function shapeForProjId(pid){ const idx=projIndex()[pid]; return SHAPES[(idx||0)%SHAPES.length]; }
  function shapeForAcc(a){ return shapeForProjId(a.proj||a.projTitle||'—'); }
  function markerSVG(a,x,y,r,col){
    if (!a || ST.markers!=='shape') return `<circle cx="${x}" cy="${y}" r="${r}" fill="${col}"/>`;
    return shapeMarkup(shapeForAcc(a), x, y, r+1, col);
  }
  function regPts(cx,cy,r,n,off){ const p=[]; for(let i=0;i<n;i++){ const a=off+i*2*Math.PI/n; p.push([cx+Math.cos(a)*r, cy+Math.sin(a)*r]); } return p; }
  function starPts(cx,cy,r,n){ const p=[]; for(let i=0;i<2*n;i++){ const rr=i%2?r*0.45:r; const a=-Math.PI/2+i*Math.PI/n; p.push([cx+Math.cos(a)*rr, cy+Math.sin(a)*rr]); } return p; }
  function shapeMarkup(type,x,y,r,fill){
    const poly=pts=>`<polygon points="${pts.map(p=>p[0].toFixed(1)+','+p[1].toFixed(1)).join(' ')}" fill="${fill}"/>`;
    switch(type){
      case 'square':        return `<rect x="${x-r}" y="${y-r}" width="${2*r}" height="${2*r}" fill="${fill}"/>`;
      case 'triangleUp':    return poly([[x,y-r],[x-r,y+r],[x+r,y+r]]);
      case 'triangleDown':  return poly([[x,y+r],[x-r,y-r],[x+r,y-r]]);
      case 'diamond':       return poly([[x,y-r],[x+r,y],[x,y+r],[x-r,y]]);
      case 'star':          return poly(starPts(x,y,r,5));
      case 'pentagon':      return poly(regPts(x,y,r,5,-Math.PI/2));
      case 'hexagon':       return poly(regPts(x,y,r,6,0));
      case 'plus':          { const t=r*0.42; return `<path d="M${x-t} ${y-r}H${x+t}V${y-t}H${x+r}V${y+t}H${x+t}V${y+r}H${x-t}V${y+t}H${x-r}V${y-t}H${x-t}Z" fill="${fill}"/>`; }
      case 'ex':            { const t=r*0.42,s=r; return `<path d="M${x-s} ${y-s+t} L${x-s+t} ${y-s} L${x} ${y-t} L${x+s-t} ${y-s} L${x+s} ${y-s+t} L${x+t} ${y} L${x+s} ${y+s-t} L${x+s-t} ${y+s} L${x} ${y+t} L${x-s+t} ${y+s} L${x-s} ${y+s-t} L${x-t} ${y} Z" fill="${fill}"/>`; }
      case 'circle': default: return `<circle cx="${x}" cy="${y}" r="${r}" fill="${fill}"/>`;
    }
  }

  /* ----- SVG renderers ----- */
  function labelText(name){ const a=metaFor(name); return name; }
  function tipTitle(name){ const a=metaFor(name); if(!a) return name;
    return `${a.id}\nfounder: ${a.founder}\nproject: ${a.projTitle||a.proj||'—'}${a.group?'\ngroup: '+a.group:''}`; }

  function rectSVG(b){
    const leaves=collectLeaves(b.tree.root);
    const rowH=18, pad=24;
    const H=Math.max(240, leaves.length*rowH + pad*2);
    const labelW=Math.min(320, 30+Math.max(...b.labels.map(l=>l.length))*7);
    const W=Math.max(760, 900);
    const treeW=W-labelW-pad;
    layoutRect(b.tree.root, treeW, H, pad);
    let edges='', tips='';
    (function walk(n){
      if(n.children){
        const ys=n.children.map(c=>c.y);
        edges+=`<path d="M${n.x} ${Math.min(...ys)} L${n.x} ${Math.max(...ys)}" stroke="#9fb0c8" stroke-width="1.3" fill="none"/>`;
        n.children.forEach(c=>{ edges+=`<path d="M${n.x} ${c.y} L${c.x} ${c.y}" stroke="#9fb0c8" stroke-width="1.3" fill="none"/>`; c.__px=n.x; walk(c); });
      }
    })(b.tree.root);
    leaves.forEach(l=>{
      const a=metaFor(l.name), col=a?tipColor(a):'#5b6b83';
      tips+=`<g><title>${esc(tipTitle(l.name))}</title>`+
        `<line x1="${l.x}" y1="${l.y}" x2="${treeW+6}" y2="${l.y}" stroke="#e7ebf1" stroke-width="1" stroke-dasharray="2 3"/>`+
        `${markerSVG(a,l.x,l.y,3.4,col)}`+
        `<text x="${treeW+12}" y="${l.y+3.5}" font-size="11" font-family="var(--mono)" fill="#25324a">${esc(labelText(l.name))}</text></g>`;
    });
    return svgWrap(W,H,`${edges}${tips}`);
  }

  function radialSVG(b){
    const size=Math.max(640, Math.min(60+b.labels.length*16, 1100));
    const pad=140;
    const lay=layoutRadial(b.tree.root, size, size, pad);
    let edges='', tips='';
    (function walk(n){ (n.children||[]).forEach(c=>{ edges+=`<path d="M${n.x} ${n.y} L${c.x} ${c.y}" stroke="#9fb0c8" stroke-width="1.2" fill="none"/>`; walk(c); }); })(b.tree.root);
    collectLeaves(b.tree.root).forEach(l=>{
      const a=metaFor(l.name), col=a?tipColor(a):'#5b6b83';
      const right=Math.cos(l.ang)>=0;
      const deg=l.ang*180/Math.PI + (right?0:180);
      tips+=`<g><title>${esc(tipTitle(l.name))}</title>`+
        `${markerSVG(a,l.x,l.y,3.2,col)}`+
        `<text x="${l.x}" y="${l.y}" font-size="10.5" font-family="var(--mono)" fill="#25324a" `+
        `text-anchor="${right?'start':'end'}" transform="rotate(${deg} ${l.x} ${l.y})" dx="${right?6:-6}" dy="3">${esc(labelText(l.name))}</text></g>`;
    });
    return svgWrap(size,size,`${edges}${tips}`);
  }

  /* distance-matrix rendering now lives in SNPMatrix (snpmatrix.js) */

  function svgWrap(W,H,inner){
    return `<svg id="snptreeSVG" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" style="max-width:100%;height:auto;background:#fff">${inner}</svg>`;
  }

  /* ------------------------------------------------------------------ *
   *  7. Text output + downloads
   * ------------------------------------------------------------------ */
  function currentText(){
    const b=ST.built; if(!b) return '';
    if (ST.format==='mega')   return mega(b.M, b.labels);
    if (ST.format==='phylip') return phylip(b.M, b.labels);
    return toNewick(b.tree.root);
  }
  function refreshText(){
    const ta=document.getElementById('fmtOut'); const ti=document.getElementById('fmtTitle');
    if(!ta) return;
    if(!ST.built) build();
    ta.value=currentText();
    if(ti) ti.textContent = ST.format==='mega'?'Pairwise diversity (MEGA)':ST.format==='phylip'?'PHYLIP distance matrix':'Newick';
  }
  function copyText(){ const ta=document.getElementById('fmtOut'); if(!ta)return; ta.select();
    try{ navigator.clipboard.writeText(ta.value); }catch(e){ document.execCommand('copy'); } }
  function downloadText(){
    const ext=ST.format==='mega'?'meg':ST.format==='phylip'?'phy':'nwk';
    dl(new Blob([currentText()],{type:'text/plain'}), fname(ext));
  }

  function svgEl(){ return document.getElementById('snptreeSVG'); }
  function serializeSVG(){ const s=svgEl(); if(!s) return null;
    const clone=s.cloneNode(true);
    // inline the mono font so exports render consistently
    clone.querySelectorAll('text').forEach(t=>{ if((t.getAttribute('font-family')||'').includes('mono')) t.setAttribute('font-family','monospace'); });
    return '<?xml version="1.0" encoding="UTF-8"?>\n'+new XMLSerializer().serializeToString(clone);
  }
  function downloadSVG(){ const s=serializeSVG(); if(!s)return; dl(new Blob([s],{type:'image/svg+xml'}), fname('svg')); }

  function rasterize(cb){
    const s=svgEl(); if(!s){cb(null);return;}
    const W=+s.getAttribute('width'), H=+s.getAttribute('height');
    const scale=2;
    const img=new Image();
    const svg=serializeSVG();
    const url='data:image/svg+xml;base64,'+btoa(unescape(encodeURIComponent(svg)));
    img.onload=()=>{ const c=document.createElement('canvas'); c.width=W*scale; c.height=H*scale;
      const ctx=c.getContext('2d'); ctx.fillStyle='#fff'; ctx.fillRect(0,0,c.width,c.height);
      ctx.setTransform(scale,0,0,scale,0,0); ctx.drawImage(img,0,0); cb(c); };
    img.onerror=()=>cb(null);
    img.src=url;
  }
  function downloadPNG(){ rasterize(c=>{ if(!c)return; c.toBlob(b=>dl(b,fname('png')),'image/png'); }); }
  function downloadPDF(){
    // dependency-free: open the image in a print window sized to the tree; user saves as PDF
    rasterize(c=>{ if(!c)return; const data=c.toDataURL('image/png');
      const w=window.open('','_blank'); if(!w){ alert('Please allow pop-ups to export PDF.'); return; }
      const region=`${ST.input.chr}:${(+ST.input.start).toLocaleString()}–${(+ST.input.end).toLocaleString()}`;
      w.document.write(`<html><head><title>SNPTree ${region}</title>
        <style>@page{size:auto;margin:12mm} body{margin:0;font-family:sans-serif}
        h3{font:600 13px sans-serif;margin:0 0 6px} img{max-width:100%}</style></head>
        <body onload="setTimeout(()=>window.print(),250)">
        <h3>SNPTree · ${ST.method.toUpperCase()} · ${region} · ${ST.input.datasetName||ST.input.dataset}</h3>
        <img src="${data}"></body></html>`);
      w.document.close();
    });
  }

  /* ------------------------------------------------------------------ *
   *  helpers
   * ------------------------------------------------------------------ */
  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;'); }
  function notice(html){ return `<div style="text-align:center;padding:28px;color:var(--muted);max-width:620px;margin:0 auto">${html}</div>`; }
  function fname(ext){ const r=`${ST.input.chr}_${ST.input.start}_${ST.input.end}`; return `snptree_${ST.method}_${r}.${ext}`; }
  function dl(blob,name){ const u=URL.createObjectURL(blob); const a=document.createElement('a'); a.href=u; a.download=name;
    document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(u),1500); }

  function injectCSS(){
    if (document.getElementById('snptree-css')) return;
    const s=document.createElement('style'); s.id='snptree-css';
    s.textContent=`
      .fl-lbl{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:7px}
      .meta-strip{margin:0 0 14px}
      .meta-row{display:flex;gap:16px;flex-wrap:wrap;font-size:12.5px;color:var(--muted);padding:10px 14px;background:var(--paper);border:1px solid var(--line);border-radius:10px}
      .meta-row b{color:var(--ink)}
      .legend-wrap{display:flex;gap:14px;flex-wrap:wrap;margin-top:10px;font-size:11.5px;color:var(--muted)}
      .legend-wrap .li{display:inline-flex;align-items:center;gap:6px}
      .legend-wrap .sw{width:11px;height:11px;border-radius:3px;display:inline-block}
      .btn.solid{background:var(--green-50,#e9f4ec);border-color:#bfe5cb;color:var(--green-600,#1f8a4c)}
      .mtx-note{font-size:12px;color:var(--muted);background:var(--blue-50,#eef4ff);border:1px solid var(--line);border-radius:8px;padding:9px 12px;margin-bottom:12px;max-width:820px}
      .mtx-note b{color:var(--ink)}`;
    document.head.appendChild(s);
  }

  /* register with the suite shell */
  if (typeof SNPTools !== 'undefined') SNPTools.register('snptree', { render });

  return { render, draw, setOpt, toCompare, forceBuild, copyText, downloadText,
           downloadSVG, downloadPNG, downloadPDF,
           // exposed for testing
           ibsMatrix, upgma, nj, toNewick, phylip, mega, dosage };
})();
if (typeof window !== 'undefined') window.SNPTree = SNPTree;
