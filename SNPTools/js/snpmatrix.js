/* =====================================================================
 *  snpmatrix.js — SNPMatrix: identity-by-state distance matrix.
 *
 *  Moved out of SNPTree. Builds the pairwise IBS distance matrix from an
 *  in-memory genotype matrix (SNPVersity / SNPTree hand-off) and renders a
 *  heatmap. Adds clustered ordering, %-identity view, bioproject bars, and
 *  downloads (CSV distance matrix, PHYLIP, PNG, SVG).
 *
 *  Hand-off: S.matrixInput = {rows, accs, chr, start, end, dataset, ...}
 *  (falls back to S.treeInput). Depends on core.js (SNPTools,S,go,ICONS).
 * ===================================================================== */
const SNPMatrix = (function () {

  const ST = { input:null, order:'input', value:'distance', bars:true, force:false, built:null, _sig:null };

  /* ---- IBS distance from the genotype matrix ----
     gts entries are 1-byte codes from Data.parseVcf: 0/1/2 = dosage, 3 = missing. */
  function dosage(gt){ return (gt == null || gt === 3) ? null : gt; }
  function ibsMatrix(rows, n){
    const M = Array.from({length:n}, () => new Float64Array(n));
    const C = Array.from({length:n}, () => new Float64Array(n));
    const d = new Array(n);
    for (const r of rows){
      const g = r.gts;
      for (let k=0;k<n;k++) d[k]=dosage(g[k]);
      for (let i=0;i<n;i++){ const di=d[i]; if (di==null) continue;
        for (let j=i+1;j<n;j++){ const dj=d[j]; if (dj==null) continue;
          M[i][j] += Math.abs(di-dj)/2; C[i][j] += 1; } }
    }
    for (let i=0;i<n;i++){ M[i][i]=0; C[i][i]=rows.length;
      for (let j=i+1;j<n;j++){ const c=C[i][j]; const dist=c? M[i][j]/c : 1;
        M[i][j]=M[j][i]=dist; C[j][i]=C[i][j]; } }
    return {M, C, sites: rows.length};
  }

  /* ---- UPGMA leaf order (for clustered display) ---- */
  function clusterOrder(M){
    const n=M.length; if (n<3) return M.map((_,i)=>i);
    let nodes = M.map((_,i)=>({idx:i, children:null, size:1}));
    const D = M.map(r=>Array.from(r));
    while (nodes.length>1){
      let mi=0,mj=1,min=Infinity;
      for (let a=0;a<nodes.length;a++) for (let b=a+1;b<nodes.length;b++) if (D[a][b]<min){min=D[a][b];mi=a;mj=b;}
      const A=nodes[mi], B=nodes[mj], U={idx:null, children:[A,B], size:A.size+B.size};
      const newRow=[];
      for (let k=0;k<nodes.length;k++){ if(k===mi||k===mj){newRow.push(0);continue;}
        newRow.push((A.size*D[mi][k]+B.size*D[mj][k])/(A.size+B.size)); }
      const keep=[]; for(let k=0;k<nodes.length;k++) if(k!==mi&&k!==mj) keep.push(k);
      const nn=keep.map(k=>nodes[k]); nn.push(U);
      const nD=keep.map(k=>{ const row=keep.map(k2=>D[k][k2]); row.push(newRow[k]); return row; });
      const last=keep.map(k=>newRow[k]); last.push(0); nD.push(last);
      nodes=nn; D.length=0; Array.prototype.push.apply(D,nD);
    }
    const out=[]; (function rec(nd){ if(!nd)return; if(!nd.children){out.push(nd.idx);return;} nd.children.forEach(rec); })(nodes[0]);
    return out;
  }

  /* ---- render ---- */
  function render(){
    injectCSS();
    const crumb=document.getElementById('crumbTool'); if(crumb) crumb.innerHTML='<b>SNPMatrix</b>';
    const page=document.getElementById('page');
    if (S.matrixInput && S.matrixInput !== ST.input){ ST.input = S.matrixInput; ST.force = false; }
    else if (!ST.input && S.treeInput) ST.input = S.treeInput;
    const inp = ST.input;
    page.className='page fade';

    if (!inp || !inp.accs || inp.accs.length<2){
      page.innerHTML = shellHead() + `<div class="empty-state"><div class="ei">${ICONS.grid||ICONS.compare||''}</div>
        <h3>Send a region from SNPVersity</h3>
        <p>SNPMatrix computes the pairwise identity-by-state distance among the accessions you queried.
        Run a query in SNPVersity (pick at least two accessions) and use <b>Send to SNPMatrix</b>,
        or open it from SNPTree.</p>
        <button class="btn solid" onclick="go('snpversity')">Go to SNPVersity</button></div>`;
      return;
    }
    page.innerHTML = shellHead() + controls(inp) + `<div id="mtxStage" class="card pad">${stageInner(inp)}</div>`;
  }

  function shellHead(){
    return `<section class="sec"><div class="bar"></div><div style="width:100%">
      <h1>SNPMatrix</h1>
      <p>Pairwise identity-by-state distances across your selected accessions, as a heatmap and downloadable matrix.</p>
    </div></section>`;
  }
  function controls(inp){
    const region = `${inp.chr}:${(+inp.start).toLocaleString()}–${(+inp.end).toLocaleString()}`;
    return `<div class="card pad" style="margin-bottom:14px;display:flex;gap:20px;flex-wrap:wrap;align-items:flex-end">
      <div><div class="fl-lbl">Order</div>${seg('order','input','Input')}${seg('order','clustered','Clustered')}</div>
      <div><div class="fl-lbl">Values</div>${seg('value','distance','IBS distance')}${seg('value','identity','% identity')}</div>
      <div><div class="fl-lbl">Bioproject bars</div>${seg('bars',true,'On')}${seg('bars',false,'Off')}</div>
      <div style="margin-left:auto;display:flex;gap:8px;flex-wrap:wrap">
        <button class="qbtn" onclick="SNPMatrix.download('csv')">Download matrix (CSV)</button>
        <button class="qbtn" onclick="SNPMatrix.download('phylip')">PHYLIP</button>
        <button class="qbtn" onclick="SNPMatrix.download('png')">PNG</button>
        <button class="qbtn" onclick="SNPMatrix.download('svg')">SVG</button>
        <button class="qbtn" onclick="SNPMatrix.toTree()">Open in SNPTree</button>
      </div>
      <div style="flex:1 1 100%;color:var(--muted);font-size:12px">Region <span class="c-mono" style="color:var(--blue-600)">${region}</span> · ${inp.datasetName?inp.datasetName+' · ':''}${inp.accs.length} accessions</div>
    </div>`;
  }
  function seg(key,val,label){
    const active = String(ST[key])===String(val);
    return `<button class="qbtn ${active?'solid':''}" onclick="SNPMatrix.set('${key}','${val}')">${label}</button>`;
  }

  function build(inp){
    const sig = `${inp.chr}:${inp.start}-${inp.end}:${inp.accs.length}:${(inp.rows||[]).length}`;
    if (ST.built && ST._sig===sig) return ST.built;
    const n = inp.accs.length;
    const {M,C,sites} = ibsMatrix(inp.rows, n);
    ST.built = {M, C, sites, n, accs:inp.accs, labels:inp.accs.map(a=>a.id)};
    ST._sig = sig;
    return ST.built;
  }
  function orderedIdx(b){ return ST.order==='clustered' ? clusterOrder(b.M) : b.M.map((_,i)=>i); }

  function stageInner(inp){
    const n = inp.accs.length;
    const V = (inp.rows && inp.rows.length) || 0;
    if (ibsWork(V,n) > IBS_COST.workBlock){
      return `<div class="mtx-note">Selected data too large for SNPMatrix - choose a smaller region, a lower-density SNP set, or fewer accessions, then re-run in SNPVersity.</div>`;
    }
    if ((n > IBS_COST.mtxAccWarn || ibsWork(V,n) > IBS_COST.workWarn) && !ST.force){
      const bits = [];
      if (n > IBS_COST.mtxAccWarn) bits.push(`a ${n}×${n} heatmap is ${(n*n).toLocaleString()} cells — slow to draw and hard to read`);
      if (ibsWork(V,n) > IBS_COST.workWarn) bits.push(`the distance matrix will make the page unresponsive for roughly ${ibsSeconds(V,n)}s`);
      return `<div class="mtx-note"><b>${n}</b> accessions × <b>${V.toLocaleString()}</b> variants — ${bits.join('; ')}.
        <button class="btn" style="margin-left:8px" onclick="SNPMatrix.force()">Build anyway</button></div>`;
    }
    const b = build(inp);
    return heatmap(b);
  }

  /* ---- heatmap ---- */
  function projColor(a){ return (a && a.projColor) || '#8aa0c0'; }
  function heatmap(b){
    const n=b.n, ord=orderedIdx(b), labels=b.labels;
    const cell=Math.max(5, Math.min(26, Math.floor(760/n)));
    const barW = ST.bars ? Math.max(4, Math.min(10, cell)) : 0;
    const gap = ST.bars ? 3 : 0;
    const showLabels = n<=80 && cell>=9;
    const labPad = showLabels ? 150 : 16;
    const pad = labPad + barW + gap;
    const W = pad + n*cell + 20, H = pad + n*cell + 20;
    let cells='';
    for (let ii=0;ii<n;ii++){ const oi=ord[ii];
      for (let jj=0;jj<n;jj++){ const oj=ord[jj];
        const v=b.M[oi][oj], t=1-v;
        const col=`rgb(${Math.round(37+(1-t)*180)},${Math.round(99+(1-t)*90)},${Math.round(235-(1-t)*120)})`;
        const shown = ST.value==='identity' ? `${(t*100).toFixed(1)}% identity` : `IBS distance ${v.toFixed(3)}`;
        cells+=`<rect x="${pad+jj*cell}" y="${pad+ii*cell}" width="${cell}" height="${cell}" fill="${col}">`+
          `<title>${esc(labels[oi])} vs ${esc(labels[oj])}\n${shown} · ${(t*100).toFixed(1)}% identity · ${b.C[oi][oj]} shared sites</title></rect>`;
      }
    }
    let bars='';
    if (ST.bars){
      for (let ii=0;ii<n;ii++){ const a=b.accs[ord[ii]], col=projColor(a);
        bars+=`<rect x="${pad-barW-gap+0}" y="${pad+ii*cell}" width="${barW}" height="${cell}" fill="${col}"><title>${esc(a.projTitle||a.proj||'')}</title></rect>`;
        bars+=`<rect x="${pad+ii*cell}" y="${pad-barW-gap+0}" width="${cell}" height="${barW}" fill="${col}"><title>${esc(a.projTitle||a.proj||'')}</title></rect>`;
      }
    }
    let text='';
    if (showLabels){
      for (let ii=0;ii<n;ii++){ const L=esc(labels[ord[ii]].slice(0,22));
        text+=`<text x="${labPad-6}" y="${pad+ii*cell+cell/2+3}" font-size="9" text-anchor="end" font-family="var(--mono)" fill="#5b6b83">${L}</text>`;
        text+=`<text x="${pad+ii*cell+cell/2}" y="${labPad-6}" font-size="9" text-anchor="start" font-family="var(--mono)" fill="#5b6b83" transform="rotate(-90 ${pad+ii*cell+cell/2} ${labPad-6})">${L}</text>`;
      }
    }
    const svg = `<svg id="snpmatrixSVG" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" style="max-width:100%;height:auto;background:#fff">${bars}${cells}${text}</svg>`;
    const note = showLabels ? '' :
      `<div class="mtx-note">Labels are hidden above 80 accessions (you have <b>${n}</b>). Hover any cell for the full IDs, ${ST.value==='identity'?'% identity':'IBS distance'}, and shared-site count.${ST.order==='clustered'?' Rows/columns are ordered by UPGMA clustering.':''}</div>`;
    const meanId = matMeanIdentity(b);
    const legend = `<div class="legend-wrap" style="margin-top:10px">
        <span class="li"><span class="sw" style="background:rgb(37,99,235)"></span>identical (0 distance)</span>
        <span class="li"><span class="sw" style="background:rgb(217,189,115)"></span>divergent</span>
        <span class="li" style="margin-left:auto;color:var(--muted)">${b.sites.toLocaleString()} sites · mean pairwise identity ${(meanId*100).toFixed(1)}%</span>
      </div>${ST.bars?projLegend(b):''}`;
    return note + svg + legend;
  }
  function matMeanIdentity(b){ let s=0,c=0; for(let i=0;i<b.n;i++)for(let j=i+1;j<b.n;j++){ s+=1-b.M[i][j]; c++; } return c?s/c:1; }
  function projLegend(b){
    const seen={}, out=[];
    b.accs.forEach(a=>{ const k=a.proj||a.projTitle||''; if(k && !seen[k]){ seen[k]=1; out.push(a); } });
    if (out.length<2) return '';
    return `<div class="legend-wrap" style="margin-top:6px;flex-wrap:wrap">`+
      out.slice(0,24).map(a=>`<span class="li"><span class="sw" style="background:${projColor(a)}"></span>${esc(a.projTitle||a.proj)}</span>`).join('')+`</div>`;
  }

  /* ---- downloads ---- */
  function orderedMatrixText(sep, asIdentity){
    const b=build(ST.input), ord=orderedIdx(b), labels=b.labels;
    const head = [''].concat(ord.map(o=>labels[o])).join(sep);
    const lines=[head];
    for (let ii=0;ii<b.n;ii++){ const oi=ord[ii];
      const row=[labels[oi]];
      for (let jj=0;jj<b.n;jj++){ const v=b.M[oi][ord[jj]]; row.push((asIdentity?(1-v):v).toFixed(6)); }
      lines.push(row.join(sep));
    }
    return lines.join('\n');
  }
  function phylipText(){
    const b=build(ST.input), ord=orderedIdx(b), labels=b.labels;
    let out=' '+b.n+'\n';
    for (let ii=0;ii<b.n;ii++){ const oi=ord[ii];
      let name=labels[oi].slice(0,10).padEnd(10,' ');
      const row=[]; for (let jj=0;jj<b.n;jj++) row.push(b.M[oi][ord[jj]].toFixed(4));
      out += name+' '+row.join(' ')+'\n';
    }
    return out;
  }
  function saveText(text, name, mime){
    const blob=new Blob([text],{type:mime||'text/plain'}), u=URL.createObjectURL(blob), a=document.createElement('a');
    a.href=u; a.download=name; document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(u),1500);
  }
  function download(kind){
    if (!ST.input) return;
    const base=`snpmatrix_${(ST.input.chr||'region')}_${ST.input.start}_${ST.input.end}`;
    if (kind==='csv')    return saveText(orderedMatrixText(',', ST.value==='identity'), `${base}_${ST.value}.csv`, 'text/csv');
    if (kind==='phylip') return saveText(phylipText(), `${base}.phy`, 'text/plain');
    if (kind==='svg'){ const el=document.getElementById('snpmatrixSVG'); if(el) saveText(el.outerHTML, `${base}.svg`, 'image/svg+xml'); return; }
    if (kind==='png'){
      const el=document.getElementById('snpmatrixSVG'); if(!el) return;
      const xml=new XMLSerializer().serializeToString(el);
      const img=new Image(); const svgBlob=new Blob([xml],{type:'image/svg+xml'}), url=URL.createObjectURL(svgBlob);
      img.onload=()=>{ const vb=el.viewBox.baseVal, scale=2;
        const cv=document.createElement('canvas'); cv.width=(vb.width||el.width.baseVal.value)*scale; cv.height=(vb.height||el.height.baseVal.value)*scale;
        const cx=cv.getContext('2d'); cx.fillStyle='#fff'; cx.fillRect(0,0,cv.width,cv.height); cx.scale(scale,scale); cx.drawImage(img,0,0);
        cv.toBlob(bl=>{ const u=URL.createObjectURL(bl), a=document.createElement('a'); a.href=u; a.download=`${base}.png`; document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(u),1500); URL.revokeObjectURL(url); }); };
      img.onerror=()=>URL.revokeObjectURL(url); img.src=url;
    }
  }

  /* ---- public ---- */
  function set(k,v){ ST[k] = (v==='true')?true:(v==='false')?false:v;
    const stage=document.getElementById('mtxStage'); if(stage) stage.innerHTML=stageInner(ST.input);
    document.querySelectorAll('.page .qbtn').forEach(()=>{}); render(); }
  function force(){ ST.force=true; const stage=document.getElementById('mtxStage'); if(stage) stage.innerHTML=stageInner(ST.input); }
  function toTree(){ if(ST.input){ S.treeInput=ST.input; } go('snptree'); }

  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  function injectCSS(){
    if (document.getElementById('snpmatrix-css')) return;
    const s=document.createElement('style'); s.id='snpmatrix-css';
    s.textContent=`.fl-lbl{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:7px}
      .mtx-note{font-size:12px;color:var(--muted);background:var(--blue-50,#eef4ff);border:1px solid var(--line);border-radius:8px;padding:9px 12px;margin-bottom:12px}
      .mtx-note b{color:var(--ink)}
      .legend-wrap{display:flex;gap:14px;align-items:center;font-size:12px;color:var(--muted)}
      .legend-wrap .li{display:inline-flex;align-items:center;gap:6px}
      .legend-wrap .sw{width:12px;height:12px;border-radius:3px;display:inline-block}`;
    document.head.appendChild(s);
  }

  if (typeof SNPTools!=='undefined') SNPTools.register('snpmatrix', { render });
  return { render, set, force, download, toTree, ibsMatrix, clusterOrder, _ST:ST };
})();
if (typeof window!=='undefined') window.SNPMatrix = SNPMatrix;
