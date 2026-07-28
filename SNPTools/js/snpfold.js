/* =====================================================================
 *  snpfold.js — Structure-aware variant interpretation.
 *
 *  Integrates coding variants with predicted protein structure three ways:
 *    1. a linear "protein browser" track that aligns variants, Pfam domains,
 *       secondary structure, and per-residue pLDDT confidence by residue;
 *    2. an interactive 3D viewer (3Dmol.js, loaded on demand) colorable by
 *       confidence / domain / variant impact, with variant residues marked;
 *    3. a per-variant structural-context readout that interprets each change
 *       (domain, local confidence, secondary structure, predicted ΔΔG).
 *
 *  Structure data + PDB come from structure-<gene>.js via Data.structureFor()
 *  / Data.pdbFor(); variants from Data.queryFoldVariants().
 * ===================================================================== */
(function () {
  const THREEDMOL_URL = 'https://3Dmol.org/build/3Dmol-min.js';
  const W = 1000;   // svg user-units width (scales to container)

  const FD = {
    gene: 'Zm00001eb406050',
    struct: null, pdb: null, variants: [],
    selId: null,
    colorMode: 'plddt',     // plddt | domain | impact
    showVar: true,
    viewer: null, libState: 'idle',
    dataset: null, sec: false,  // sec = show PlantCAD2/ESM2/ESM3 (MaizeGDB 2026 only)
    structSource: null,  // 'alphafold' | 'boltz' | 'esmfold' | null — which folder the current FD.struct/FD.pdb came from
    modelPref: 'best',   // 'best' | 'alphafold' | 'boltz' | 'esmfold' — user's model choice for the next load
    carriers: null, openCarrier: null,   // pos|ref|alt -> {carriersHom,carriersHet,het,hom} (whole-panel, via geneFunction)
    locus: null,         // {chr,start,end,dataset} for this gene — used for the SNPVersity handoff
    sort: { key: null, dir: 'asc' },     // variant-table sort: column key + direction ('asc'|'desc'); null key = file order
    pendingVariant: null,// {chr,pos,ref,alt,sub,…} handed over by another tool; selected once variants are in hand
    truncation: null,    // {structureLength,maxVariantResidue,beyondCount,sourceLabel} when variants extend past the loaded model
    root: null,          // persistent DOM container — survives navigation to other tools
    loaded: false,       // a gene's heavy content is (being) rendered into root
    loadedGene: null,    // which gene that content is for
    pdbResidueMap: null, // residue number -> one-letter amino acid, built lazily from the loaded PDB
  };

  /* ---------- palettes ---------- */
  function plddtHex(b){ return b>=90?'#0053d6': b>=70?'#65cbf3': b>=50?'#ffdb13':'#ff7d45'; }
  function plddtBand(b){ return b>=90?'Very high': b>=70?'Confident': b>=50?'Low':'Very low'; }
  function impactHex(s){
    const n = finiteNumber(s);
    if (n == null) return '#8a93a3';
    const lo=-12,hi=4,t=Math.max(0,Math.min(1,(n-lo)/(hi-lo)));
    const r=t<.5?220:Math.round(220*(1-(t-.5)*2)), g=t<.5?Math.round(180*t*2):170;
    return `rgb(${r},${Math.max(50,g)},55)`;
  }
  const DOM_FILL = ['#cdeccf','#cfe0fb','#efe0fb','#fde9cc'];
  const CONS_FILL = { lof:'#d6322a', missense:'#2f5bbf', lod:'#b54708', splice:'#6d28d9', indel:'#176c3a', syn:'#8a93a3' };
  const SS_LABEL = { H:'α-helix', E:'β-strand', C:'loop / coil' };

  /* ---------- structure source (AlphaFold2 vs Boltz2) ----------
     Checked in this order — first folder that actually has the gene's file wins. */
  const STRUCT_SOURCES = [
    { key: 'alphafold', dir: './data/structures/alphafold', label: 'AlphaFold2', badge: 'AF2' },
    { key: 'boltz',      dir: './data/structures/boltz',     label: 'Boltz2',     badge: 'B2'  },
    { key: 'esmfold',    dir: './data/structures/esmfold',   label: 'ESMFold',    badge: 'ESM'  },
  ];

  /* Loads structure-<gene>.js as a <script> tag (same mechanism the app already uses to
     bring in per-gene structure files) and resolves true once it has run, or rejects if
     the file 404s / errors — which is how we tell "not in this folder" from "loaded". */
  function loadScriptOnce(src){
    return new Promise((resolve, reject) => {
      const tag = 'data-foldsrc';
      const existing = document.querySelector(`script[${tag}="${src}"]`);
      if (existing){ resolve(true); return; }
      const el = document.createElement('script');
      el.src = src;
      el.async = true;
      el.setAttribute(tag, src);
      el.onload = () => resolve(true);
      el.onerror = () => { el.remove(); reject(new Error('not found: ' + src)); };
      document.head.appendChild(el);
    });
  }

  /* Data.structureFor(gene) / Data.pdbFor(gene) is a SINGLE shared slot per gene: whichever
     structure-<gene>.js last actually *executed* is what's sitting in it. loadScriptOnce also
     skips re-running a <script> that's already in the DOM. Put those two things together and
     re-selecting a source you've already loaded before would read Data's shared slot as-is
     without re-executing that source's script — so it silently returns whatever source last
     wrote there, not the one you just asked for. In practice that source is whichever one was
     loaded last (e.g. ESMFold, if it's last in the 'best' fallback order or was picked at some
     point), and it then "sticks" no matter what you pick afterward.
     Fix: snapshot each source's struct/pdb into our OWN cache the instant its script loads,
     before any other source's script can overwrite Data's shared slot — then always read from
     this cache rather than re-reading Data.structureFor/pdbFor after the fact. */
  const structCache = Object.create(null);   // "<gene>::<sourceKey>" -> { struct, pdb }

  async function resolveStructureSource(gene, pref){
    const order = (!pref || pref === 'best') ? STRUCT_SOURCES : STRUCT_SOURCES.filter(s => s.key === pref);
    for (const src of order){
      const cacheKey = gene + '::' + src.key;
      if (structCache[cacheKey]) return { src, data: structCache[cacheKey] };  // already resolved this one — use our snapshot, not Data's shared slot
      try { await loadScriptOnce(`${src.dir}/structure-${gene}.js`); }
      catch (e){ continue; }             // no file for this gene in this folder — try next
      const struct = Data.structureFor(gene);
      if (struct){
        const data = { struct, pdb: Data.pdbFor(gene) };
        structCache[cacheKey] = data;    // snapshot now, before another source's script can overwrite Data's slot
        return { src, data };
      }
    }
    return null;
  }

  /* small "AlphaFold2" / "Boltz2" pill shown next to the gene once a model has loaded */
  function structSourceBadge(){
    const src = STRUCT_SOURCES.find(s => s.key === FD.structSource);
    if (!src) return '';
    return `<span class="fold-src-badge ${src.key}" title="Predicted structure from ${src.label}">
      <span class="fold-src-dot"></span>${src.label}</span>`;
  }

  /* language-model score pill — same gradient (red→green) as SNPImpact / SNPFunction */
  function finiteNumber(v){
    if (v == null || v === '' || v === '.' || v === 'NA' || v === 'N/A') return null;
    const n = Number(v);
    return Number.isFinite(n) ? n : null;
  }
  function scoreColor(v){
    const n = finiteNumber(v);
    if (n == null) return '#8a93a3';
    const lo=-12,hi=6,t=Math.max(0,Math.min(1,(n-lo)/(hi-lo)));
    const r=t<.5?255:Math.round(255*(1-(t-.5)*2));
    const g=t<.5?Math.round(255*t*2):200;
    return `rgb(${r},${Math.max(60,g)},60)`;
  }

  /* Choose black or white text from the actual chip background. The 0.179
     luminance cutoff is the point at which black and white have equal WCAG
     contrast; pale yellow/green chips therefore get dark text automatically. */
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
    const n = finiteNumber(v);
    if (n == null) return '<span style="color:var(--faint)">—</span>';
    const bg = scoreColor(n);
    const fg = contrastTextColor(bg);
    return `<span class="imp-score" style="background:${bg};color:${fg} !important;text-shadow:none">${n>0?'+':''}${n.toFixed(1)}</span>`;
  }
  function scoreText(v){
    const n = finiteNumber(v);
    return n == null ? '—' : `${n > 0 ? '+' : ''}${n.toFixed(1)}`;
  }

  /* Support both the established fields and common 2026 aliases. */
  const SCORE_KEYS = {
    plantcad:  ['plantcad', 'plantcad1', 'plantcad_1', 'PLANTCAD', 'PLANTCAD1'],
    plantcad2: ['plantcad2', 'plantcad_2', 'PLANTCAD2'],
    esm:       ['esm', 'esm1', 'esm_1', 'ESM', 'ESM1'],
    esm2:      ['esm2', 'esm_2', 'ESM2'],
    esm3:      ['esm3', 'esm_3', 'ESM3'],
  };
  function modelScore(v, model){
    if (!v) return null;
    const keys = SCORE_KEYS[model] || [model];
    for (const key of keys){
      if (Object.prototype.hasOwnProperty.call(v, key) && finiteNumber(v[key]) != null) return finiteNumber(v[key]);
    }
    return null;
  }
  function hasSecondaryVariantScores(variants){
    return (variants || []).some(v =>
      modelScore(v, 'plantcad2') != null ||
      modelScore(v, 'esm2') != null ||
      modelScore(v, 'esm3') != null
    );
  }
  function aiScoreSummary(v){
    const scores = [
      ['PlantCAD', modelScore(v, 'plantcad')],
      ...(FD.sec ? [['PlantCAD2', modelScore(v, 'plantcad2')]] : []),
      ['ESM1', modelScore(v, 'esm')],
      ...(FD.sec ? [
        ['ESM2', modelScore(v, 'esm2')],
        ['ESM3', modelScore(v, 'esm3')],
      ] : []),
    ];
    if (!scores.some(([, value]) => value != null)) return 'n/a';
    //return scores.map(([label, value]) => `${label} ${scoreText(value)}`).join(' · ');
    return scores.map(([label, value]) => `${label} ${scoreText(value)}`).join('<br>');
  }

  function escFold(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  function carrierKey(v){ return v.pos+'|'+v.refNt+'|'+v.altNt; }
  function carrierOf(v){ return FD.carriers ? (FD.carriers[carrierKey(v)] || null) : null; }

  /* SNPVersity annotates domains by genomic position, while SNPFold needs
     protein-coordinate intervals. Load the canonical protein-domain index
     directly and try the gene, protein/transcript IDs, and common B73-v5
     isoform aliases. Per-variant genomic domain text remains a fallback. */
  function usableDomainText(v){
    if (v == null) return null;
    const s = String(v).trim();
    return (!s || s === '—' || s === '-' || /^N\/?A$/i.test(s)) ? null : s;
  }
  function domainFromText(value){
    const text = usableDomainText(value);
    if (!text) return null;
    const m = text.match(/(PF\d{4,6})/i);
    const pfam = m ? m[1].toUpperCase() : null;
    const name = pfam ? text.replace(/\s*\(?\s*PF\d{4,6}\s*\)?\s*/i, ' ').trim() : text;
    return { name:name || text, pfam, kind:'domain', source:'variant' };
  }
  function normalizeDomainRecord(raw){
    if (!raw) return null;
    if (Array.isArray(raw)){
      const start = finiteNumber(raw[0]), end = finiteNumber(raw[1]);
      if (start == null || end == null) return null;
      return { start:Math.min(start,end), end:Math.max(start,end),
        name:String(raw[2] || raw[3] || 'Pfam domain'),
        pfam:usableDomainText(raw[3]), kind:String(raw[4] || 'domain').toLowerCase() === 'region' ? 'region' : 'domain' };
    }
    if (typeof raw !== 'object') return null;
    const start = finiteNumber(raw.start != null ? raw.start : (raw.aaStart != null ? raw.aaStart : (raw.begin != null ? raw.begin : raw.from)));
    const end   = finiteNumber(raw.end   != null ? raw.end   : (raw.aaEnd   != null ? raw.aaEnd   : (raw.stop  != null ? raw.stop  : raw.to)));
    if (start == null || end == null) return null;
    const kind = String(raw.kind || raw.type || 'domain').toLowerCase() === 'region' ? 'region' : 'domain';
    return { ...raw, start:Math.min(start,end), end:Math.max(start,end), kind,
      name:String(raw.name || raw.label || raw.domain || raw.description || 'Pfam domain'),
      pfam:usableDomainText(raw.pfam || raw.accession || raw.pfamId || raw.id) };
  }
  function addUniqueValue(list, value){
    if (value == null || value === '') return;
    const text = String(value).trim();
    if (text && !list.includes(text)) list.push(text);
  }
  function canonicalDomainKeys(fn){
    const keys = [];
    addUniqueValue(keys, FD.gene);
    if (fn){
      addUniqueValue(keys, fn.gene); addUniqueValue(keys, fn.protein);
      addUniqueValue(keys, fn.proteinId); addUniqueValue(keys, fn.transcript);
      addUniqueValue(keys, fn.transcriptId);
    }
    const st = FD.struct || {};
    [st.gene, st.id, st.protein, st.proteinId, st.protein_id,
     st.transcript, st.transcriptId, st.transcript_id, st.model].forEach(v=>addUniqueValue(keys,v));
    const gene = String(FD.gene || '').trim();
    if (gene){
      ['_P001','.P001','-P001','_P002','.P002','-P002',
       '_T001','.T001','-T001','_T002','.T002','-T002'].forEach(suffix=>addUniqueValue(keys,gene+suffix));
    }
    return keys;
  }
  async function loadCanonicalDomainRecords(fn){
    const records = [];
    if (fn && Array.isArray(fn.domains)) records.push({key:FD.gene, domains:fn.domains});
    if (typeof Data === 'undefined' || typeof Data.geneDomains !== 'function') return records;
    if (typeof Data.ensureGeneDomains === 'function'){
      try { await Data.ensureGeneDomains(); } catch (e) { /* retain structure/variant domains */ }
    }
    canonicalDomainKeys(fn).forEach(key=>{
      try {
        const gd = Data.geneDomains(key);
        if (gd && Array.isArray(gd.domains) && gd.domains.length) records.push({key, ...gd});
      } catch (e) { /* try the remaining aliases */ }
    });
    return records;
  }
  function mergeCanonicalDomains(fn, domainRecords){
    if (!FD.struct) return;
    const current = (FD.struct.domains || []).map(normalizeDomainRecord).filter(Boolean);
    const canonical = [];
    if (fn && Array.isArray(fn.domains)) canonical.push(...fn.domains);
    (domainRecords || []).forEach(record=>{
      if (record && Array.isArray(record.domains)) canonical.push(...record.domains);
    });
    const out = [], seen = new Set(), L = finiteNumber(FD.struct.length);
    current.concat(canonical.map(normalizeDomainRecord).filter(Boolean)).forEach(d => {
      if (L != null){
        d.start = Math.max(1, Math.min(L, d.start));
        d.end   = Math.max(1, Math.min(L, d.end));
        if (d.end < d.start) return;
      }
      const key = [d.start,d.end,d.name,d.pfam||'',d.kind].join('|');
      if (!seen.has(key)){ seen.add(key); out.push(d); }
    });
    out.sort((a,b)=>a.start-b.start || a.end-b.end);
    FD.struct.domains = out;
  }

  function foldVariantKey(v){
    if (!v || v.pos == null) return null;
    const ref = v.refNt != null ? v.refNt : v.ref;
    const alt = v.altNt != null ? v.altNt : v.alt;
    return [v.pos, ref == null ? '' : ref, alt == null ? '' : alt].join('|');
  }
  function severeFoldConsequence(v){
    const s = String((v && v.consequence) || '').toLowerCase();
    return (v && v.consClass === 'lof') || /stop gained|nonsense|frameshift|start lost|stop lost/.test(s);
  }
  function foldRelevant(v){ return v && (v.consClass === 'missense' || v.consClass === 'lof' || v.consClass === 'indel'); }

  const AA3_FOLD = {ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',
    LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V',TER:'*',SEC:'U'};
  /* Parse only genuine protein notation. A genomic fallback such as
     "123456 A>G" starts with a number and must never become residue 123456. */
  function aaFromVariantLabel(label){
    const raw = String(label || '').trim();
    if (!raw) return {ref:'',resi:null,alt:''};
    const hadProteinPrefix = /^p\./i.test(raw);
    const s = raw.replace(/^p\./i,'').replace(/[()\s]/g,'');
    if (!hadProteinPrefix && !/^(?:[A-Z*]|[A-Za-z]{3})\d+/.test(s)) return {ref:'',resi:null,alt:''};
    let m = s.match(/^([A-Z*])(\d+)(fs\*?|ins|del|ext\*?\d*|\?|[A-Z*]|=)?$/i);
    if (m) return {ref:(m[1]||'').toUpperCase(), resi:Number(m[2]), alt:(m[3]||'').replace(/^fs\*$/i,'fs')};
    m = s.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3}|Ter|\*|fs\*?|ins|del|ext\*?\d*|\?)?/i);
    if (m){
      const ref = AA3_FOLD[m[1].toUpperCase()] || '';
      const tok = m[3] || '';
      const alt = AA3_FOLD[tok.toUpperCase()] || (/^fs/i.test(tok)?'fs':tok);
      return {ref, resi:Number(m[2]), alt};
    }
    return {ref:'',resi:null,alt:''};
  }
  function plausibleResidue(value, maxLength){
    const n = finiteNumber(value), max = finiteNumber(maxLength);
    if (n == null || n < 1 || Math.floor(n) !== n) return null;
    if (max != null && n > max) return null;
    return n;
  }
  function modelProteinLength(model){
    const cds = modelCds(model);
    if (!cds.length) return null;
    const nt = cds.reduce((sum,iv)=>sum+(iv[1]-iv[0]+1),0);
    return nt > 0 ? Math.floor(nt/3) : null;
  }

  const AA1_FOLD = new Set('ARNDCQEGHILKMFPSTWYVUO'.split(''));
  function cleanAaLetter(value){
    const aa = String(value == null ? '' : value).trim().toUpperCase();
    return aa.length === 1 && AA1_FOLD.has(aa) ? aa : '';
  }
  function sequenceAaAt(resi){
    if (!FD.struct || resi == null) return '';
    const fields = ['sequence','seq','proteinSequence','protein_sequence','aaSequence','aa_sequence'];
    for (const key of fields){
      let seq = FD.struct[key];
      if (Array.isArray(seq)) seq = seq.join('');
      if (typeof seq !== 'string') continue;
      seq = seq.replace(/[^A-Za-z*]/g,'').toUpperCase();
      const aa = cleanAaLetter(seq.charAt(resi-1));
      if (aa) return aa;
    }
    return '';
  }
  function buildPdbResidueMap(){
    if (FD.pdbResidueMap) return FD.pdbResidueMap;
    const map = new Map();
    String(FD.pdb || '').split(/\r?\n/).forEach(line => {
      if (!/^(ATOM  |HETATM)/.test(line) || line.length < 26) return;
      const resi = Number.parseInt(line.slice(22,26).trim(),10);
      if (!Number.isFinite(resi) || map.has(resi)) return;
      const aa = AA3_FOLD[line.slice(17,20).trim().toUpperCase()] || '';
      if (aa && aa !== '*') map.set(resi, aa);
    });
    FD.pdbResidueMap = map;
    return map;
  }
  function structureAaAt(resi){
    const r = plausibleResidue(resi);
    if (r == null) return '';
    return sequenceAaAt(r) || buildPdbResidueMap().get(r) || '';
  }
  function referenceAaAt(v, resi, supplied){
    /* The loaded protein is the authority. This prevents a genomic REF base
       (A/C/G/T) from being mistaken for a one-letter amino-acid code. */
    const fromStructure = structureAaAt(resi);
    if (fromStructure) return fromStructure;
    const explicit = [v && v.aaRef, v && v.refAa, v && v.refAA,
      v && v.referenceAa, v && v.referenceAA, supplied];
    for (const value of explicit){
      const aa = cleanAaLetter(value);
      if (aa) return aa;
    }
    return '';
  }
  function modelCds(model){
    if (!model) return [];
    let cds = model.cds || model.CDS || model.cdsIntervals || model.coding ||
      (model.intervals_by_feature && model.intervals_by_feature.cds) || [];
    if (!Array.isArray(cds)) return [];
    return cds.map(iv=>{
      if (Array.isArray(iv)) return [finiteNumber(iv[0]),finiteNumber(iv[1])];
      if (iv && typeof iv === 'object') return [finiteNumber(iv.start),finiteNumber(iv.end)];
      return [null,null];
    }).filter(iv=>iv[0]!=null && iv[1]!=null).map(iv=>[Math.min(iv[0],iv[1]),Math.max(iv[0],iv[1])]);
  }
  /* Convert a one-based genomic position to a one-based amino-acid residue by
     walking the canonical CDS in transcript order. */
  function genomicToResidue(pos, model){
    pos = finiteNumber(pos);
    if (pos == null) return null;
    const strand = String((model && model.strand) || '+') === '-' ? '-' : '+';
    const cds = modelCds(model).sort((a,b)=>strand==='-' ? b[0]-a[0] : a[0]-b[0]);
    let codingBefore = 0;
    for (const [start,end] of cds){
      if (pos >= start && pos <= end){
        const within = strand === '-' ? end-pos : pos-start;
        return plausibleResidue(Math.floor((codingBefore + within)/3) + 1);
      }
      codingBefore += end-start+1;
    }
    return null;
  }
  async function loadFoldGeneModel(fn){
    if (typeof Data === 'undefined' || typeof Data.geneModelOf !== 'function') return null;
    let chr = fn && fn.chr;
    if (!chr && typeof Data.lookupGene === 'function'){
      try { const g = await Data.lookupGene(FD.gene); chr = g && g.chr; } catch (e) { /* no coordinate fallback */ }
    }
    if (!chr) return null;
    if (typeof Data.ensureGeneModels === 'function'){
      try { await Data.ensureGeneModels(chr); } catch (e) { return null; }
    }
    try { return Data.geneModelOf(chr, FD.gene); } catch (e) { return null; }
  }
  function proteinVariantLabel(v, resi, ref, alt){
    if (resi == null) return v.variant || (v.pos != null ? `${v.pos} ${v.ref||v.refNt||''}>${v.alt||v.altNt||''}` : v.consequence);
    const consequence = String(v.consequence || '').toLowerCase();
    const r = referenceAaAt(v, resi, ref);
    /* Compact one-letter display used elsewhere in SNPFold. When the loaded
       structure still cannot identify the residue, use X rather than the
       invalid HGVS-like form p.?298*. */
    const displayRef = r || 'X';
    if (/stop gained|nonsense/.test(consequence)) return `p.${displayRef}${resi}*`;
    if (/frameshift/.test(consequence)) return `p.${displayRef}${resi}fs`;
    if (/start lost/.test(consequence)) return `p.${displayRef === 'X' ? 'M' : displayRef}${resi}?`;
    if (/stop lost/.test(consequence)) return `p.*${resi}${alt && alt!=='*' ? alt : 'ext'}`;
    if (/in-frame insertion|inframe insertion/.test(consequence)) return `p.${displayRef}${resi}ins`;
    if (/in-frame deletion|inframe deletion/.test(consequence)) return `p.${displayRef}${resi}del`;
    return `p.${displayRef}${resi}${alt || '?'}`;
  }
  function resolvedResidue(v, model){
    const parsed = aaFromVariantLabel(v && v.variant);
    const proteinLength = modelProteinLength(model);
    const parsedResi = plausibleResidue(parsed.resi, proteinLength);
    const convertedResi = genomicToResidue(v && v.pos, model);
    const genomicPos = finiteNumber(v && v.pos);
    /* Explicit residue fields are a final fallback because some data sources have
       historically populated them with the genomic DNA position. Validate them
       against the canonical CDS length when available and reject an exact copy of
       the genomic coordinate. */
    const explicit = [v && v.resi, v && v.residue, v && v.aaPos, v && v.aa_position,
      v && v.proteinPos, v && v.protein_position]
      .map(x=>plausibleResidue(x, proteinLength))
      .find(x=>x!=null && (genomicPos==null || x!==genomicPos));
    return { parsed, resi: parsedResi != null ? parsedResi :
      (convertedResi != null ? convertedResi : (explicit != null ? explicit : null)) };
  }
  function normalizeFunctionVariant(v, index, model){
    if (!foldRelevant(v)) return null;
    const resolved = resolvedResidue(v, model), aa = resolved.parsed, resi = resolved.resi;
    const ref = referenceAaAt(v, resi, v.aaRef || v.refAa || aa.ref || '');
    const alt = v.aaAlt || v.altAa || aa.alt || null;
    return {
      id:'ffn'+index, gene:FD.gene, resi, ref, alt,
      variant:proteinVariantLabel(v, resi, ref, alt),
      consequence:v.consequence || (v.consClass === 'lof' ? 'Loss-of-function' : 'Coding variant'),
      consClass:v.consClass, structural:v.consClass === 'missense',
      pos:v.pos, refNt:v.refNt != null ? v.refNt : v.ref, altNt:v.altNt != null ? v.altNt : v.alt,
      impact:v.impact || v.impactLevel || null, maf:v.maf != null ? v.maf : v.af, domain:v.domain,
      plantcad:v.plantcad != null ? v.plantcad : v.pc1, plantcad2:v.plantcad2 != null ? v.plantcad2 : v.pc2,
      esm:v.esm != null ? v.esm : v.esm1, esm2:v.esm2, esm3:v.esm3, combined:v.combined,
      priority:severeFoldConsequence(v) ? 'TOP' : (v.priority || null),
    };
  }
  function normalizePrimaryVariant(v, model, index){
    const copy = {...v};
    const resolved = resolvedResidue(copy, model);
    copy.resi = resolved.resi;
    if (copy.resi != null){
      copy.ref = referenceAaAt(copy, copy.resi, copy.ref || resolved.parsed.ref || '');
      copy.alt = copy.alt || resolved.parsed.alt || null;
      if (!aaFromVariantLabel(copy.variant).resi || /^\d+\s+[ACGTN-]+>[ACGTN-]+$/i.test(String(copy.variant||'')))
        copy.variant = proteinVariantLabel(copy, copy.resi, copy.ref, copy.alt);
    }
    if (severeFoldConsequence(copy)) copy.priority = 'TOP';
    if (!copy.id) copy.id = 'fp'+index;
    return copy;
  }
  function mergeFoldVariants(primary, fn, model){
    const out = (primary || []).map((v,i) => normalizePrimaryVariant(v, model, i));
    const byKey = new Map();
    out.forEach(v => { const key = foldVariantKey(v); if (key) byKey.set(key, v); });
    ((fn && fn.variants) || []).forEach((fv, i) => {
      const normalized = normalizeFunctionVariant(fv, i, model);
      if (!normalized) return;
      const key = foldVariantKey(normalized);
      const existing = key ? byKey.get(key) : null;
      if (existing){
        if (!usableDomainText(existing.domain) && usableDomainText(normalized.domain)) existing.domain = normalized.domain;
        if (existing.resi == null && normalized.resi != null) existing.resi = normalized.resi;
        if ((!existing.variant || /^\d+\s/.test(existing.variant)) && normalized.variant) existing.variant = normalized.variant;
        if (severeFoldConsequence(normalized)) existing.priority = 'TOP';
        return;
      }
      if (normalized.resi != null || severeFoldConsequence(normalized)){
        out.push(normalized);
        if (key) byKey.set(key, normalized);
      }
    });
    return out.sort((a,b) => {
      const ar = plausibleResidue(a.resi), br = plausibleResidue(b.resi);
      if (ar == null && br == null) return (finiteNumber(a.pos)||0) - (finiteNumber(b.pos)||0);
      if (ar == null) return 1; if (br == null) return -1;
      return ar - br || ((finiteNumber(a.pos)||0) - (finiteNumber(b.pos)||0));
    });
  }

  /* Detect a likely C-terminally truncated structure by comparing the model's
     residue count with protein-residue positions present in the coding-variant
     list. We intentionally use the normalized residue positions, not genomic
     coordinates. A warning is shown only when at least one coding variant lies
     beyond the final modeled residue. */
  function detectStructureTruncation(){
    const structureLength = FD.struct ? finiteNumber(FD.struct.length) : null;
    if (structureLength == null || structureLength < 1){ FD.truncation = null; return null; }
    const beyond = (FD.variants || []).map(v=>({v, resi:finiteNumber(v.resi)}))
      .filter(x=>x.resi != null && x.resi >= 1 && Math.floor(x.resi) === x.resi && x.resi > structureLength);
    if (!beyond.length){ FD.truncation = null; return null; }
    const maxVariantResidue = Math.max(...beyond.map(x=>x.resi));
    const src = STRUCT_SOURCES.find(x=>x.key===FD.structSource);
    FD.truncation = {
      structureLength,
      maxVariantResidue,
      beyondCount:beyond.length,
      sourceLabel:(src && src.label) || 'The selected structure model',
    };
    return FD.truncation;
  }
  function truncationWarningHTML(){
    const t = FD.truncation;
    if (!t) return '';
    return `<div class="fold-trunc-warning" role="alert">
      <div class="fold-trunc-icon" aria-hidden="true">⚠</div>
      <div><b>Possible truncated protein structure</b>
        <div>${escFold(t.sourceLabel)} contains <b>${t.structureLength}</b> modeled amino acids, but the coding-variant list extends to residue <b>${t.maxVariantResidue}</b>${t.beyondCount>1?` (${t.beyondCount} variants fall beyond the model)`:''}. The structure prediction method most likely truncated this protein. Check the other structure models above for a more complete protein prediction.</div>
      </div>
    </div>`;
  }

  /* ---------- active dataset detection ---------- */
  function nonEmpty(v){ return v != null && v !== ''; }
  function datasetText(dataset){
    if (dataset == null) return '';
    if (typeof dataset === 'string' || typeof dataset === 'number') return String(dataset);
    if (typeof dataset === 'object'){
      return [dataset.id, dataset.key, dataset.value, dataset.name, dataset.label, dataset.title]
        .filter(nonEmpty).join(' ');
    }
    return String(dataset);
  }
  function datasetRepresentations(dataset){
    const out = [];
    const add = v => {
      if (!nonEmpty(v)) return;
      if (!out.some(x => x === v)) out.push(v);
    };
    add(dataset);
    if (dataset && typeof dataset === 'object'){
      add(dataset.id); add(dataset.key); add(dataset.value);
      add(dataset.name); add(dataset.label); add(dataset.title);
    }
    return out;
  }
  function resolveDataset(){
    const candidates = [];
    const add = v => { if (nonEmpty(v)) candidates.push(v); };

    if (typeof S !== 'undefined' && S){
      add(S.dataset);
      add(S.datasetId);
      add(S.datasetKey);
      add(S.selectedDataset);
      add(S.activeDataset);
      add(S.currentDataset);
      add(S.snpDataset);
      add(S.ds);
    }

    if (typeof Data !== 'undefined' && Data){
      if (typeof Data.currentDataset === 'function'){
        try { add(Data.currentDataset()); } catch (e) { /* use other candidates */ }
      } else {
        add(Data.currentDataset);
      }
      add(Data.dataset);
      add(Data.datasetId);
      add(Data.datasetKey);
      add(Data.selectedDataset);
      add(Data.activeDataset);
    }

    if (typeof window !== 'undefined'){
      add(window.CURRENT_DATASET);
      add(window.currentDataset);
      add(window.activeDataset);
    }

    add(FD.dataset);
    return candidates.length ? candidates[0] : null;
  }
  function hasSecondaryScores(dataset){
    /* Prefer the application's own helper, while trying both object and ID/name forms. */
    if (typeof Data !== 'undefined' && typeof Data.hasSecondaryScores === 'function'){
      const reps = datasetRepresentations(dataset);
      if (!reps.length) reps.push(undefined);
      for (const rep of reps){
        try {
          if (Data.hasSecondaryScores(rep) === true) return true;
        } catch (e) { /* fall through to normalized-name check */ }
      }
    }

    /* Fallback accepts labels/IDs such as MaizeGDB2026, maizegdb_2026_hq,
       or "MaizeGDB 2026 (High Coverage)". */
    const normalized = datasetText(dataset).toLowerCase().replace(/[^a-z0-9]+/g, '');
    return normalized.includes('maizegdb2026');
  }

  /* ---------- dataset chooser (compact; shares Data.datasets() with SNPVersity) ---------- */
  function foldDatasets(){
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
    const a = String(datasetId(d)), b = String(datasetId(FD.dataset));
    if (a && a === b) return true;
    const na = datasetText(d).toLowerCase().replace(/[^a-z0-9]+/g, '');
    const nb = datasetText(FD.dataset).toLowerCase().replace(/[^a-z0-9]+/g, '');
    return !!na && na === nb;
  }
  function datasetLabel(d){
    const name = (d && nonEmpty(d.name)) ? String(d.name) : datasetText(d);
    let sub = (d && nonEmpty(d.sub)) ? String(d.sub) : '';
    if (sub && name.toLowerCase().includes(sub.toLowerCase())) sub = '';
    return { name: name || 'Dataset', sub };
  }
  function datasetCardsHTML(){
    const list = foldDatasets();
    if (!list.length) return '';
    return list.map(d=>{
      const id  = datasetId(d);
      const sel = isCurrentDataset(d);
      const { name, sub } = datasetLabel(d);
      return `<button type="button" class="fold-ds ${sel?'sel':''}" onclick="FOLD.pickDataset('${escFold(id)}')">
        <span class="fold-ds-dot"></span>
        <span class="fold-ds-txt"><span class="fold-ds-name">${escFold(name)}</span>${sub?`<span class="fold-ds-sub">${escFold(sub)}</span>`:''}</span>
      </button>`;
    }).join('');
  }
  function datasetChooser(){
    const cards = datasetCardsHTML();
    if (!cards) return '';   // if the data layer exposes no catalog, show nothing
    return `<div class="card pad fold-ds-card" style="margin-bottom:16px">
      <div class="fold-ds-head">Dataset</div>
      <div class="fold-ds-grid" id="foldDsGrid">${cards}</div>
    </div>`;
  }
  function syncDatasetChooser(){
    const grid = FD.root && FD.root.querySelector('#foldDsGrid');
    if (grid) grid.innerHTML = datasetCardsHTML();
  }

  /* ---------- structural context for a residue ---------- */
  function domainAt(resi, variant){
    const ds = FD.struct.domains || [];
    const r = finiteNumber(resi);
    const structural = r == null ? null
      : (ds.find(d=>d.kind==='domain' && r>=d.start && r<=d.end)
        || ds.find(d=>r>=d.start && r<=d.end) || null);
    return structural || domainFromText(variant && variant.domain);
  }
  function ctxFor(v){
    const r = finiteNumber(v.resi);
    const i = r == null ? -1 : r-1;
    const inModel = FD.struct.plddt && i>=0 && i<FD.struct.plddt.length;
    const plddt = inModel ? FD.struct.plddt[i] : null;
    const ss = inModel && FD.struct.ss && FD.struct.ss[i] ? FD.struct.ss[i] : 'C';
    return { plddt, ss, ssLabel:SS_LABEL[ss], domain:domainAt(r, v), inModel };
  }
  function interpret(v, c){
    if (!c.inModel){
      if (finiteNumber(v.resi) == null){
        const where = v.pos != null ? ` at genomic position ${v.pos}` : '';
        const tone = v.consClass === 'lof' ? 'lof' : 'mid';
        return [tone, `${v.consequence}${where} — the VCF does not provide a protein residue, so it is listed but cannot be placed on the structure.`];
      }
      return ['mid', `Residue ${v.resi} is outside the modeled region (1–${FD.struct.length} aa) — shown in the list, but not placed on the structure.`];
    }
    const inDom = c.domain && c.domain.kind==='domain';
    const where = inDom ? (' in the '+c.domain.name+' domain')
                        : (c.domain && c.domain.kind==='region' ? (' in the '+c.domain.name) : ' in an inter-domain region');
    if (v.consClass==='lof'){
      const after = inDom ? ('within / before the '+c.domain.name+' domain') : 'in the C-terminal region';
      const what = v.consequence==='Start lost' ? 'Disrupts the start codon'
                 : v.consequence==='Stop lost'  ? 'Removes the stop codon'
                 : 'Truncates the protein';
      return ['lof', `${what} ${after} — position shown; downstream effect not modeled.`];
    }
    if (v.consClass==='indel'){
      return ['mid', `${v.consequence}${where} — in-frame length change; position shown, folding effect not modeled.`];
    }
    if (v.combined==null){
      return ['mid', `Missense${where} — no language-model score available.`];
    }
    const sev = v.combined<=-6 ? 'Likely damaging' : v.combined<=-3 ? 'Possibly damaging' : 'Likely tolerated';
    const wellStruct = c.plddt>=70;
    const tone = sev==='Likely damaging' ? 'bad' : sev==='Possibly damaging' ? 'mid' : 'ok';
    return [tone, `${sev}: ${wellStruct?'well-structured':'low-confidence'} ${c.ssLabel}${where}.`];
  }

  /* =================== RENDER =================== */
  function modelPrefRadios(){
    const opts = [{ key:'best', label:'Best' }, ...STRUCT_SOURCES.map(s=>({ key:s.key, label:s.label }))];
    return `<div class="fold-model-pref" role="radiogroup" aria-label="Structure model">
      <span style="font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px">Model</span>
      ${opts.map(o => `<label class="fold-model-opt">
        <input type="radio" name="foldModelPref" value="${o.key}" ${FD.modelPref===o.key?'checked':''}
          onchange="FOLD.setModelPref('${o.key}')"> ${o.label}
      </label>`).join('')}
    </div>`;
  }

  function searchBar(){
    return `<div class="card pad" style="margin-bottom:16px;display:flex;gap:10px;align-items:center;flex-wrap:wrap">
      <span style="font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px">Gene model</span>
      <input id="foldGeneInput" value="${FD.gene||''}" placeholder="e.g. Zm00001eb406050" spellcheck="false"
        style="flex:1;min-width:280px;border:1px solid var(--line);border-radius:9px;padding:9px 11px;font-family:var(--mono);font-size:13px"
        onkeydown="if(event.key==='Enter')FOLD.loadGene()">
      <button class="btn" onclick="FOLD.loadGene()">Load structure</button>
      ${modelPrefRadios()}
    </div>`;
  }

  /* ---------- persistent container ----------
     SNPFold's content is expensive to render (structure fetch + 3D viewer), so we
     build it once into FD.root and re-attach that same node whenever the tool is
     shown again — instead of rebuilding from scratch. Because we hold a JS reference
     to FD.root, it survives another tool overwriting #page (it is simply detached,
     not destroyed), which lets the 3D view stay warm across navigation. */
  function ensureRoot(page){
    if (!FD.root){
      FD.root = document.createElement('div');
      FD.root.className = 'snpfold-root';
    }
    if (FD.root.parentNode !== page){
      page.innerHTML = '';
      page.appendChild(FD.root);
    }
  }

  /* Entry point called by the suite shell on navigation.
     - Arriving from another tool (SNPVersity sets S.foldGene) → autoload that gene.
     - Arriving from the side-panel menu → do NOT autoload; show a landing page the
       first time, and simply re-show the already-loaded page on later visits. */
  async function render(page){
    page = page || document.getElementById('page');
    injectCSS();
    ensureRoot(page);

    /* keep FD.dataset in sync with the rest of the app (side panel / other tools) */
    const activeDataset = resolveDataset();
    if (activeDataset != null) FD.dataset = activeDataset;
    FD.sec = hasSecondaryScores(FD.dataset);

    const requestedGene = (typeof S !== 'undefined' && S && S.foldGene) ? S.foldGene : null;
    /* Optional site to highlight, e.g. the "fold ↗" link on a SNPVersity row.
       Consumed here so a stale request never leaks into a later visit. */
    const requestedVariant = (typeof S !== 'undefined' && S && S.foldVariant) ? S.foldVariant : null;
    if (requestedVariant && typeof S !== 'undefined' && S) S.foldVariant = null;
    if (requestedVariant) FD.pendingVariant = requestedVariant;

    if (requestedGene){
      /* came from another tool — honor the explicit request and autoload */
      S.foldGene = null;
      if (requestedGene !== FD.loadedGene || !FD.loaded){
        FD.gene = requestedGene;
        FD.selId = null;
        FD.openCarrier = null;
        await loadStructure();   // applies FD.pendingVariant once variants arrive
      } else {
        /* Same gene already loaded — re-show it and just move the selection.
           No structure fetch, no variant query, no HDF5 round-trip: clicking
           several "fold" links on one gene costs nothing after the first. */
        reshowLoaded();
        applyPendingVariant();
      }
      return;
    }

    /* came from the side-panel menu */
    if (FD.loaded){
      reshowLoaded();            // preserve the rendered page, don't reload
      applyPendingVariant();
      return;
    }
    renderLanding();             // first visit from the menu: no autoload
  }

  /* =================== incoming variant selection ===================
     Another tool (SNPVersity) can name a specific site alongside the gene. The
     descriptor is loose on purpose — genomic position, alleles, and/or a protein
     substitution — so we match on whatever is most specific and available. */
  function normSub(s){ return String(s==null?'':s).replace(/^p\./i,'').replace(/[()\s]/g,'').toUpperCase(); }
  function matchFoldVariant(req){
    const list = FD.variants || [];
    if (!req || !list.length) return null;

    /* 1 · an explicit SNPFold variant id wins outright */
    if (req.id != null){
      const byId = list.find(v => String(v.id) === String(req.id));
      if (byId) return byId;
    }

    /* 2 · genomic position (+ alleles when several records share a position) */
    const pos = finiteNumber(req.pos);
    if (pos != null){
      const here = list.filter(v => finiteNumber(v.pos) === pos);
      if (here.length === 1) return here[0];
      if (here.length > 1){
        const ref = String(req.ref || '').toUpperCase();
        const alt = String(req.alt || '').toUpperCase();
        const exact = here.find(v =>
          String(v.refNt != null ? v.refNt : '').toUpperCase() === ref &&
          String(v.altNt != null ? v.altNt : '').toUpperCase() === alt);
        if (exact) return exact;
        return here[0];
      }
    }

    /* 3 · protein substitution, e.g. "A123T" — compare against the printed
           label and against ref+resi+alt, since notations differ between tools */
    const sub = normSub(req.sub || req.variant);
    if (sub){
      const bySub = list.find(v =>
        normSub(v.variant) === sub ||
        normSub(`${v.ref||''}${v.resi}${v.alt||''}`) === sub);
      if (bySub) return bySub;

      /* fall back to the residue number alone, if the substitution carries one */
      const m = sub.match(/(\d+)/);
      if (m){
        const resi = Number(m[1]);
        const byResi = list.find(v => finiteNumber(v.resi) === resi);
        if (byResi) return byResi;
      }
    }
    return null;
  }

  /* Select FD.pendingVariant against the variants currently in memory. Purely a
     view operation: it re-renders the track, table, and context panel and moves
     the 3D camera — it never refetches. */
  function applyPendingVariant(){
    const req = FD.pendingVariant;
    FD.pendingVariant = null;
    if (!req) return;

    const v = matchFoldVariant(req);
    if (!v){
      const label = req.sub || req.variant || (req.pos != null ? `position ${req.pos}` : 'that variant');
      foldNotice(`${escFold(label)} is not among the coding variants SNPFold has for ${escFold(FD.gene)}.`);
      return;
    }

    FD.selId = v.id;
    refreshSelection();
    focusResidue(true);
    scrollToSelectedRow();
  }

  function scrollToSelectedRow(){
    if (!FD.root) return;
    const run = () => {
      const row = FD.root.querySelector('tr.fold-row.sel');
      if (row && typeof row.scrollIntoView === 'function'){
        row.scrollIntoView({ behavior:'smooth', block:'center' });
      }
    };
    if (typeof requestAnimationFrame === 'function') requestAnimationFrame(run); else run();
  }

  /* small dismissible bar at the top of the tool — used when a handoff cannot be honored */
  function foldNotice(html){
    if (!FD.root) return;
    const old = FD.root.querySelector('.fold-notice');
    if (old) old.remove();
    const bar = document.createElement('div');
    bar.className = 'fold-notice';
    bar.innerHTML = `<span>${html}</span><button type="button" title="Dismiss">×</button>`;
    bar.querySelector('button').onclick = () => bar.remove();
    FD.root.insertBefore(bar, FD.root.firstChild);
  }

  /* re-show the already-rendered content without refetching */
  function reshowLoaded(){
    syncDatasetChooser();
    if (FD.struct && FD.pdb){
      if (!FD.viewer) buildViewer();           // (re)init if it never got a live DOM
      else { try { FD.viewer.resize(); FD.viewer.render(); } catch (e){} }
    }
    // the SNPVersity selection may have changed while we were away — refresh the
    // "Replace current selection (N)" hint so the count on screen is never stale
    if (typeof Handoff!=='undefined') Handoff.sync(FD.root);
    if (typeof attachTT==='function') attachTT();
  }

  /* Landing shown when the tool is opened from the menu — dataset chooser + gene
     search, but no structure/variant fetch until the user asks for it. */
  function renderLanding(){
    FD.root.innerHTML = datasetChooser() + searchBar() + `
      <div class="empty-state"><div class="ei">${ICONS.fold}</div>
        <h3>Structure-aware variant interpretation</h3>
        <p>Pick a dataset above and a gene model, then press <b>Load structure</b> to map
        coding variants onto the predicted protein — a linear browser, an interactive 3D
        model, and per-variant structural context.</p>
        <div style="margin-top:14px"><button class="btn primary" onclick="FOLD.loadGene()">Load structure for ${escFold(FD.gene)}</button></div>
      </div>`;
    if (typeof attachTT==='function') attachTT();
  }

  /* The heavy path: fetch the model + variants and render the full SNPFold UI into
     FD.root. Used both for autoload-from-tool and for explicit user loads. */
  async function loadStructure(){
    FD.viewer = null;
    FD.truncation = null;
    FD.loaded = true;               // commit to the loaded view (keep across navigation)
    FD.loadedGene = FD.gene;

    const header = datasetChooser() + searchBar();

    FD.root.innerHTML = header + `<div class="loading" style="padding:48px;text-align:center">
      <div class="spinner"></div><div>Loading model for ${FD.gene}…</div></div>`;

    FD.structSource = null;
    let matched = null;
    try { matched = await resolveStructureSource(FD.gene, FD.modelPref); } catch (e){ /* no file for this gene in the checked folder(s) */ }

    FD.struct = matched ? matched.data.struct : null;
    FD.pdb    = matched ? matched.data.pdb    : null;
    FD.pdbResidueMap = null;

    if (!FD.struct || !matched){
      const checked = (!FD.modelPref || FD.modelPref === 'best')
        ? STRUCT_SOURCES
        : STRUCT_SOURCES.filter(s => s.key === FD.modelPref);
      const folderList = checked.map(s => `<span class="c-mono">data/structures/${s.key}/</span>`)
        .join(checked.length > 1 ? ', ' : '');
      const whatWasChecked = checked.length > 1
        ? `checked ${folderList} for a structure file for this gene and found none of them`
        : `checked ${folderList} for a structure file for this gene and didn't find one`;
      const modelNote = (!FD.modelPref || FD.modelPref === 'best')
        ? ''
        : ` You have <b>${STRUCT_SOURCES.find(s=>s.key===FD.modelPref)?.label}</b> selected specifically —
           switch to <b>Best</b> above to let SNPFold fall back to the other model sources.`;
      FD.root.innerHTML = datasetChooser() + searchBar() + `<div class="empty-state"><div class="ei">${ICONS.fold}</div>
        <h3>No predicted model for “${FD.gene}”</h3>
        <p>SNPFold ${whatWasChecked}. This usually means the protein wasn't modeled — most
        likely because it's part of a complex or its sequence is too long to fold.${modelNote}
        Search another gene model above, or generate its structure file and drop it in one
        of those folders.</p></div>`;
      FD.pendingVariant = null;     // nothing to select without a model
      if (typeof attachTT==='function') attachTT();
      return;
    }
    FD.structSource = matched.src.key;
    const s = FD.struct;
    const plddtValues = Array.isArray(s.plddt) ? s.plddt.map(finiteNumber).filter(v => v != null) : [];
    const meanP = plddtValues.length
      ? (plddtValues.reduce((a,b)=>a+b,0)/plddtValues.length).toFixed(0)
      : '—';

    FD.root.innerHTML = datasetChooser() + searchBar() + `<div class="loading" style="padding:48px;text-align:center">
      <div class="spinner"></div><div>Loading coding variants for ${s.gene}…</div></div>`;
    try {
      /* Passing the dataset as a second argument is backward-compatible in JavaScript:
         older one-argument implementations simply ignore it. */
      const variantsPromise = Promise.resolve(Data.queryFoldVariants(FD.gene, FD.dataset));
      const functionPromise = typeof Data.geneFunction === 'function'
        ? Promise.resolve(Data.geneFunction(FD.gene, FD.dataset)).catch(()=>null)
        : Promise.resolve(null);

      const [variants, fn] = await Promise.all([variantsPromise, functionPromise]);
      const [geneModel, domainRecords] = await Promise.all([
        loadFoldGeneModel(fn),
        loadCanonicalDomainRecords(fn),
      ]);

      /* Use canonical protein-coordinate domains even when the full-panel variant
         query fails, and convert genomic LOF positions through the canonical CDS
         before drawing the residue lollipop. */
      mergeCanonicalDomains(fn, domainRecords);
      FD.variants = mergeFoldVariants(variants, fn, geneModel);
      detectStructureTruncation();

      /* If the app state did not expose the dataset, actual returned 2026 score fields
         still enable the extra columns. */
      FD.sec = FD.sec || hasSecondaryVariantScores(FD.variants);

      FD.carriers = (fn && fn.variants)
        ? Object.fromEntries(fn.variants.map(v => [v.pos+'|'+v.ref+'|'+v.alt, v]))
        : null;

      /* keep the gene's genomic coordinates — needed to hand a region to SNPVersity */
      FD.locus = (fn && fn.chr != null)
        ? { chr:fn.chr, start:+fn.start, end:+fn.end, dataset:(fn.dataset!=null?fn.dataset:FD.dataset) }
        : null;
    }
    catch (e){
      console.error('SNPFold variant loading error', e);
      FD.variants = [];
      FD.carriers = null;
      FD.locus = null;
      FD.truncation = null;
    }

    FD.root.innerHTML = datasetChooser() + searchBar() + `
      <div class="sec"><div class="bar"></div><div>
        <div class="n">STRUCTURE-AWARE INTERPRETATION · ${STRUCT_SOURCES.find(x=>x.key===FD.structSource)?.label || 'Predicted structure'} + PlantCAD/ESM</div>
        <h2>See where variants, domains, and local structure align</h2>
        <p>Coding variants mapped onto the predicted protein. Read each change against its
        Pfam domain, secondary structure, and local model confidence (pLDDT) — in a linear
        browser and in 3D — to judge whether it strikes a structured core or a flexible loop.</p>
      </div></div>

      <div class="fold-context">
        <span class="g"><b class="mono">${s.gene}</b></span>
        ${structSourceBadge()}
        ${s.title?`<span class="dot">·</span><span>${s.title}</span>`:''}
        <span class="dot">·</span><span><b>${s.length}</b> aa</span>
        ${s.uniprot?`<span class="dot">·</span><span>UniProt <a href="https://www.uniprot.org/uniprotkb/${s.uniprot}" target="_blank" rel="noopener">${s.uniprot}</a></span>`:''}
        <span class="dot">·</span><span>mean pLDDT <b>${meanP}</b></span>
        <span class="dot">·</span><span><b>${FD.variants.length}</b> coding variants</span>
      </div>
      ${truncationWarningHTML()}

      <!-- linear protein browser -->
      <div class="sec" style="margin-top:24px"><div class="bar"></div><div><h2 style="font-size:16px">Protein browser</h2>
        <p>Variants (lollipops, height = severity) over domains, secondary structure, and pLDDT. Click a variant.</p></div></div>
      <div class="card pad">
        <div class="fold-track" id="foldTrack">${trackSVG()}</div>
        <div class="fold-legend" style="display:flex;flex-wrap:wrap;gap:10px 34px;align-items:flex-start">
          <div class="lg-group" style="display:flex;flex-direction:column;align-items:flex-start;gap:4px">
            <div class="lg-title" style="font-weight:600;color:var(--ink)">Lollipop:</div>
            <span class="lg"><span class="sw" style="background:${CONS_FILL.missense}"></span>Missense</span>
            <span class="lg"><span class="sw" style="background:${CONS_FILL.lof}"></span>Loss-of-function</span>
            <span class="lg"><span class="sw" style="background:${CONS_FILL.indel}"></span>In-frame indel</span>
          </div>
          <div class="lg-group" style="display:flex;flex-direction:column;align-items:flex-start;gap:4px">
            <div class="lg-title" style="font-weight:600;color:var(--ink)">Secondary structure:</div>
            <span class="lg"><span class="sw" style="background:#d24b6a"></span>helix</span>
            <span class="lg"><span class="sw" style="background:#2f6fd0"></span>strand</span>
            <span class="lg"><span class="sw" style="background:#cbd4e1"></span>coil</span>
          </div>
          <div class="lg-group" style="display:flex;flex-direction:column;align-items:flex-start;gap:4px">
            <div class="lg-title" style="font-weight:600;color:var(--ink)">AlphaFold confidence score (pLDDT):</div>
            <span class="lg"><span class="sw" style="background:#0053d6"></span>≥90</span>
            <span class="lg"><span class="sw" style="background:#65cbf3"></span>70–90</span>
            <span class="lg"><span class="sw" style="background:#ffdb13"></span>50–70</span>
            <span class="lg"><span class="sw" style="background:#ff7d45"></span>&lt;50</span>
          </div>
        </div>
      </div>

      <!-- 3D + context -->
      <div class="fold-main">
        <div class="card" style="overflow:hidden">
          <div class="fold-toolbar">
            <div class="seg">
              <button class="seg-b ${FD.colorMode==='plddt'?'on':''}"  onclick="FOLD.color('plddt')">Confidence</button>
              <button class="seg-b ${FD.colorMode==='domain'?'on':''}" onclick="FOLD.color('domain')">Domain</button>
              <button class="seg-b ${FD.colorMode==='impact'?'on':''}" onclick="FOLD.color('impact')">Variant impact</button>
            </div>
            <label class="chk"><input type="checkbox" ${FD.showVar?'checked':''} onchange="FOLD.toggleVar(this.checked)"> Variant residues</label>
            <span style="margin-left:auto">${structSourceBadge()}</span>
            <button class="btn ghost" onclick="FOLD.reset()">Reset view</button>
          </div>
          <div id="fold3d" class="fold-3d"></div>
        </div>
        <div class="card pad fold-ctx" id="foldCtx">${ctxHTML()}</div>
      </div>

      <!-- variant table -->
      <div class="sec" style="margin-top:24px"><div class="bar"></div>
        <div style="display:flex;align-items:center;gap:14px;flex-wrap:wrap;width:100%">
          <h2 style="font-size:16px;margin:0">Coding variants on this gene</h2>
          <span style="margin-left:auto;display:flex;gap:8px;flex-wrap:wrap">
            ${FD.gene?`<button class="btn" title="Open ${escFold(FD.gene)}'s functional dossier in SNPFunction"
              onclick="FOLD.toFunction()">${ICONS.leaf||ICONS.star||''} Send gene to SNPFunction</button>`:''}
            ${anyCarriers()?`<button class="btn"
              title="Open this gene's region in SNPVersity with every accession carrying any coding variant preselected"
              onclick="FOLD.toVersity('all','all')">${ICONS.dna||''} Send all carriers to SNPVersity</button>`:''}
          </span>
        </div>
      </div>
      ${anyCarriers()?`<div class="fn-handoff">
        <span class="fn-handoff-k">Handoff to SNPVersity</span>
        <span data-ho-mount data-ho-id="foldMergeReplace" data-ho-target="SNPVersity" data-ho-dataset="${escFold(datasetId(FD.dataset)==null?'':datasetId(FD.dataset))}"></span>
      </div>`:''}
      <div class="tbl-wrap" style="max-height:none"><table class="vcf fold-table">
        <thead id="foldTableHead">${tableHeadHTML()}</thead>
        <tbody id="foldTableBody">${sortedVariants().map(rowHTML).join('')}</tbody>
      </table></div>
    `;
    buildViewer();
    if (typeof Handoff!=='undefined') Handoff.sync(FD.root);
    if (typeof attachTT==='function') attachTT();
    applyPendingVariant();          // honor an incoming "show me this site" request
  }

  /* ---------- linear track ---------- */
  function x(r){ return ((r-1)/(FD.struct.length-1))*W; }
  function rleBands(arr, fn){
    const out=[]; let st=0;
    for (let i=1;i<=arr.length;i++){
      if (i===arr.length || fn(arr[i])!==fn(arr[st])){ out.push({from:st+1,to:i,key:fn(arr[st])}); st=i; }
    }
    return out;
  }
  function ssRuns(ss){
    const out=[]; let st=0;
    for (let i=1;i<=ss.length;i++){
      if (i===ss.length || ss[i]!==ss[st]){ out.push({from:st+1,to:i,t:ss[st]}); st=i; }
    }
    return out;
  }
  function trackSVG(){
    const s=FD.struct, N=s.length;
    const H=176, lolliTop=8, lolliH=58, base=lolliTop+lolliH;      // lollipop baseline
    const domY=base+6, domH=18, ssY=domY+domH+22, ssH=12, pY=ssY+ssH+10, pH=12;
    const maxSev=10;
    let g='';

    // pLDDT strip (RLE by band)
    rleBands(s.plddt, b=>plddtHex(b)).forEach(seg=>{
      const x1=x(seg.from), x2=x(seg.to);
      g+=`<rect x="${x1.toFixed(1)}" y="${pY}" width="${Math.max(.6,x2-x1).toFixed(1)}" height="${pH}" fill="${seg.key}"/>`;
    });
    g+=`<text x="0" y="${pY+pH+11}" class="tlab">pLDDT confidence</text>`;

    // secondary structure
    ssRuns(s.ss).forEach(r=>{
      const x1=x(r.from), x2=x(r.to), w=Math.max(.6,x2-x1);
      if (r.t==='H')      g+=`<rect x="${x1.toFixed(1)}" y="${ssY}" width="${w.toFixed(1)}" height="${ssH}" rx="3" fill="#d24b6a"/>`;
      else if (r.t==='E') g+=`<rect x="${x1.toFixed(1)}" y="${ssY+2}" width="${w.toFixed(1)}" height="${ssH-4}" fill="#2f6fd0"/>`;
      else                g+=`<rect x="${x1.toFixed(1)}" y="${ssY+ssH/2-1}" width="${w.toFixed(1)}" height="2" fill="#cbd4e1"/>`;
    });
    g+=`<text x="0" y="${ssY-4}" class="tlab">secondary structure</text>`;

    // domains
    (s.domains||[]).forEach((d,i)=>{
      const x1=x(d.start), x2=x(d.end), w=x2-x1, mid=x1+w/2;
      const label = escFold([d.name, d.pfam].filter(Boolean).join(' · '));
      if (d.kind==='region'){
        g+=`<rect x="${x1.toFixed(1)}" y="${domY+3}" width="${w.toFixed(1)}" height="${domH-6}" rx="3" fill="#eef1f6" stroke="#d3dae6" stroke-dasharray="3 2"/>`;
        g+=`<text x="${mid.toFixed(1)}" y="${domY+domH-4}" class="dlab" text-anchor="middle">${label}</text>`;
      } else {
        g+=`<rect x="${x1.toFixed(1)}" y="${domY}" width="${w.toFixed(1)}" height="${domH}" rx="5" fill="${DOM_FILL[i%DOM_FILL.length]}" stroke="#bcd0f5"/>`;
        g+=`<text x="${mid.toFixed(1)}" y="${domY+domH-5}" class="dlab" text-anchor="middle">${label}</text>`;
      }
    });
    /* If the canonical interval file has no matching key, SNPVersity's
       per-variant domain hit still appears at the converted residue as a
       narrow dashed marker. This shows a position-supported hit without
       inventing unsupported domain boundaries. */
    if (!(s.domains||[]).length){
      const seenHits = new Set();
      (FD.variants||[]).forEach(v=>{
        const resi=plausibleResidue(v.resi), text=usableDomainText(v.domain);
        if (resi==null || resi>N || !text) return;
        const key=resi+'|'+text; if (seenHits.has(key)) return; seenHits.add(key);
        const xx=x(resi), label=escFold(text);
        g+=`<g data-tt="${label} · domain hit at residue ${resi}">`;
        g+=`<rect x="${(xx-4).toFixed(1)}" y="${domY}" width="8" height="${domH}" rx="3" fill="#fff" stroke="#6d7fa7" stroke-width="1.4" stroke-dasharray="2 1"/>`;
        g+=`<text x="${xx.toFixed(1)}" y="${domY+domH+10}" class="dlab" text-anchor="middle">${label}</text></g>`;
      });
    }

    // ruler ticks
    for (let r=1;r<=N;r+=100){ const xx=x(r);
      g+=`<line x1="${xx.toFixed(1)}" y1="${base}" x2="${xx.toFixed(1)}" y2="${base+3}" stroke="#aab4c4"/>`;
      g+=`<text x="${xx.toFixed(1)}" y="${base-3}" class="rlab">${r}</text>`;
    }
    g+=`<line x1="0" y1="${base}" x2="${W}" y2="${base}" stroke="#dde3ec"/>`;

    // selection guide
    const sel = FD.variants.find(v=>v.id===FD.selId);
    const selResi = sel ? plausibleResidue(sel.resi) : null;
    if (selResi != null && selResi<=N){ const xx=x(selResi); g+=`<line x1="${xx.toFixed(1)}" y1="${lolliTop}" x2="${xx.toFixed(1)}" y2="${pY+pH}" stroke="#13264a" stroke-dasharray="2 2" opacity=".5"/>`; }

    // lollipops
    FD.variants.forEach(v=>{
      const resi = plausibleResidue(v.resi);
      if (resi == null || resi>N) return;   // only protein residues can be placed
      const xx=x(resi);
      const sev=severeFoldConsequence(v) ? maxSev
        : (v.combined!=null?Math.min(maxSev, Math.abs(v.combined)):maxSev*0.5);
      const head=base - 6 - (sev/maxSev)*(lolliH-10);
      const col=CONS_FILL[v.consClass]||'#2f5bbf';
      const on=v.id===FD.selId;
      g+=`<g class="lolli ${on?'on':''}" onclick="FOLD.select('${v.id}')" data-tt="${v.variant} · ${v.consequence} · residue ${resi}">`;
      g+=`<line x1="${xx.toFixed(1)}" y1="${base}" x2="${xx.toFixed(1)}" y2="${head.toFixed(1)}" stroke="${col}" stroke-width="${on?2:1.4}"/>`;
      g+=`<circle cx="${xx.toFixed(1)}" cy="${head.toFixed(1)}" r="${on?6:4.5}" fill="${col}" stroke="#fff" stroke-width="1.5"/>`;
      if (v.consClass==='lof') g+=`<text x="${xx.toFixed(1)}" y="${(head-8).toFixed(1)}" class="vlab" text-anchor="middle">✱</text>`;
      g+=`</g>`;
    });

    return `<svg viewBox="0 0 ${W} ${H}" class="track-svg" preserveAspectRatio="xMinYMin meet">${g}</svg>`;
  }

  /* ---------- PanEffect jump ---------- */
  function isMissense(v){ return v && v.consClass==='missense' && v.resi; }
  /* internal view switch — highlights the substitution when missense,
     otherwise just opens PanEffect on the gene */
  function panEffect(v){
    if (!FD.gene) return;
    if (typeof goPanEffect !== 'function') return go('paneffect');
    if (isMissense(v)) goPanEffect(FD.gene, {variant:{pos:v.resi, wt:v.ref||'', sub:v.alt||''}});
    else goPanEffect(FD.gene);
  }
  function peJump(v){
    if (!isMissense(v) || !FD.gene) return '';
    const sub = `${v.ref||''}${v.resi}${v.alt||''}`;
    return ` <a class="pe-jump" href="#" title="View ${escFold(sub)} in PanEffect"
      onclick="event.stopPropagation();FOLD.panEffect('${v.id}');return false;">effects ↗</a>`;
  }

  /* ---------- variant table: sortable column model ----------
     One entry per <th>, in display order. `get` returns the value the column is
     sorted on (null/undefined => always sorted to the bottom, either direction).
     `sec` marks the MaizeGDB-2026-only columns, so the header and the row markup
     stay in sync automatically. `desc1` = first click sorts high→low, which reads
     better for counts/ranks; everything else starts low→high. */
  const PRIO_RANK = { top:4, high:3, moderate:2, medium:2, low:1, modifier:0 };
  const FOLD_COLS = [
    { key:'variant',     label:'Variant',     type:'str',
      get:v => v.variant },
    { key:'consequence', label:'Consequence', type:'str',
      get:v => v.consequence },
    { key:'resi',        label:'Residue',     type:'num', num:true,
      get:v => finiteNumber(v.resi) },
    { key:'domain',      label:'Domain',      type:'str',
      get:v => { const d = ctxFor(v).domain; return d ? d.name : null; } },
    { key:'plddt',       label:'Local pLDDT', type:'num', num:true, desc1:true,
      get:v => finiteNumber(ctxFor(v).plddt) },
    { key:'ss',          label:'Structure',   type:'str',
      get:v => { const c = ctxFor(v); return c.inModel ? c.ssLabel : null; } },
    { key:'plantcad',    label:'PlantCAD',    type:'num', num:true,
      get:v => modelScore(v, 'plantcad') },
    { key:'plantcad2',   label:'PlantCAD2',   type:'num', num:true, sec:true,
      get:v => modelScore(v, 'plantcad2') },
    { key:'esm',         label:'ESM1',        type:'num', num:true,
      get:v => modelScore(v, 'esm') },
    { key:'esm2',        label:'ESM2',        type:'num', num:true, sec:true,
      get:v => modelScore(v, 'esm2') },
    { key:'esm3',        label:'ESM3',        type:'num', num:true, sec:true,
      get:v => modelScore(v, 'esm3') },
    { key:'priority',    label:'Priority',    type:'num', desc1:true,
      get:v => { const p = v.priority ? PRIO_RANK[String(v.priority).toLowerCase()] : null;
                 return p == null ? null : p; } },
    { key:'carriers',    label:'Carriers',    type:'num', num:true, desc1:true,
      get:v => { const c = carrierOf(v);
                 if (!c) return null;
                 const n = (Number(c.hom) || 0) + (Number(c.het) || 0);
                 return n === 0 ? null : n; } },
  ];
  function foldVisibleCols(){ return FOLD_COLS.filter(c => !c.sec || FD.sec); }
  function foldCol(key){ return FOLD_COLS.find(c => c.key === key) || null; }

  /* Stable sort: ties (and blanks) keep their original order, so repeated sorts
     never shuffle rows arbitrarily. */
  function sortedVariants(){
    const list = FD.variants || [];
    const col  = foldCol(FD.sort.key);
    if (!col) return list.slice();
    const dir = FD.sort.dir === 'desc' ? -1 : 1;
    return list
      .map((v, i) => ({ v, i, k: col.get(v) }))
      .sort((a, b) => {
        const aNull = a.k == null || a.k === '';
        const bNull = b.k == null || b.k === '';
        if (aNull && bNull) return a.i - b.i;
        if (aNull) return 1;          // blanks always last
        if (bNull) return -1;
        const r = col.type === 'num'
          ? (Number(a.k) - Number(b.k))
          : String(a.k).localeCompare(String(b.k), undefined, { numeric:true, sensitivity:'base' });
        return r ? r * dir : a.i - b.i;
      })
      .map(x => x.v);
  }

  function tableHeadHTML(){
    return '<tr>' + foldVisibleCols().map(c => {
      const on = FD.sort.key === c.key;
      const arrow = on ? (FD.sort.dir === 'desc' ? '▼' : '▲') : '↕';
      const tip = on
        ? `Sorted ${FD.sort.dir === 'desc' ? 'high to low' : 'low to high'} — click to reverse`
        : `Sort by ${c.label}`;
      return `<th class="fold-th${c.num ? ' num' : ''}${on ? ' sorted' : ''}" title="${escFold(tip)}"
        aria-sort="${on ? (FD.sort.dir === 'desc' ? 'descending' : 'ascending') : 'none'}"
        onclick="FOLD.sortBy('${c.key}')"><span class="fold-th-in">${escFold(c.label)}<span class="fold-ar">${arrow}</span></span></th>`;
    }).join('') + '<th class="fold-send-th"></th></tr>';
  }

  /* ---------- variant table ---------- */
  function rowHTML(v){
    const c = ctxFor(v);
    const on = v.id===FD.selId;
    const cr = carrierOf(v);
    const openC = FD.openCarrier===v.id;
    return `<tr class="fold-row ${on?'sel':''}" onclick="FOLD.select('${v.id}')">
      <td class="c-mono c-alt" style="padding-left:11px">${v.variant}</td>
      <td><span class="cons ${v.consClass}">${v.consequence}</span>${peJump(v)}</td>
      <td class="num">${finiteNumber(v.resi)==null?'<span style="color:var(--faint)">—</span>':v.resi}</td>
      <td>${c.domain ? (c.domain.kind==='domain'?`<span class="dom-tag">${c.domain.name}</span>`:`<span style="color:var(--muted);font-size:11px">${c.domain.name}</span>`) : '<span style="color:var(--faint)">—</span>'}</td>
      <td class="num">${c.plddt==null?'<span style="color:var(--faint)">—</span>':`<span class="plddt-chip" style="background:${plddtHex(c.plddt)};color:${c.plddt>=70?'#06294f':'#5c3a06'}">${c.plddt.toFixed(0)}</span>`}</td>
      <td>${c.inModel?`<span class="ss-chip ss-${c.ss}">${c.ssLabel}</span>`:'<span style="color:var(--faint)">—</span>'}</td>
      <td class="num">${scoreCell(modelScore(v, 'plantcad'))}</td>
      ${FD.sec?`<td class="num">${scoreCell(modelScore(v, 'plantcad2'))}</td>`:''}
      <td class="num">${scoreCell(modelScore(v, 'esm'))}</td>
      ${FD.sec?`<td class="num">${scoreCell(modelScore(v, 'esm2'))}</td><td class="num">${scoreCell(modelScore(v, 'esm3'))}</td>`:''}
      <td>${v.priority?`<span class="prio ${v.priority.toLowerCase()}">${v.priority}</span>`:'<span style="color:var(--faint)">—</span>'}</td>
      <td style="text-align:center">${carrierBtn(v, cr, openC)}</td>
      <td class="fold-send">${sendBtn(v, cr)}</td>
    </tr>${openC?carrierRow(v, cr):''}`;
  }
  function foldCols(){ return foldVisibleCols().length + 1; }   // colspan for the expanded carrier row (+ send column)
  /* per-variant "open in SNPVersity with these carriers preselected" */
  function sendBtn(v, cr){
    const n = cr ? ((Number(cr.hom)||0) + (Number(cr.het)||0)) : 0;
    if (!n) return '<span style="color:var(--faint)">—</span>';
    return `<button class="btn tiny" title="Open ${escFold(FD.gene)} in SNPVersity with all ${n} carrier${n>1?'s':''} of ${escFold(v.variant)} preselected"
      onclick="event.stopPropagation();FOLD.toVersity('${v.id}','all')">SNPVersity →</button>`;
  }
  function carrierBtn(v, cr, open){
    if (!cr || (cr.hom===0 && cr.het===0)) return '<span style="color:var(--faint)">—</span>';
    return `<button class="fold-cbtn ${open?'on':''}" title="Show carrier accessions (whole panel)"
      onclick="event.stopPropagation();FOLD.carriers('${v.id}')"><span class="cc hom">${cr.hom}</span><span class="cc het">${cr.het}</span></button>`;
  }
  function carrierRow(v, cr){
    if (!cr) return `<tr class="fn-carriers"><td colspan="${foldCols()}"><div class="fn-cwrap muted">No carrier data for this variant.</div></td></tr>`;
    const carriersHom = Array.isArray(cr.carriersHom) ? cr.carriersHom : [];
    const carriersHet = Array.isArray(cr.carriersHet) ? cr.carriersHet : [];
    const chip = (id,cls)=>`<span class="carrier ${cls}">${escFold(id)}</span>`;
    const homs = carriersHom.slice(0,60).map(id=>chip(id,'hom')).join('');
    const hets = carriersHet.slice(0,60).map(id=>chip(id,'het')).join('');
    const send = (mode,label,n)=> n
      ? `<button class="btn tiny" onclick="event.stopPropagation();FOLD.toVersity('${v.id}','${mode}')">${label} (${n}) →</button>`
      : '';
    return `<tr class="fn-carriers"><td colspan="${foldCols()}">
      <div class="fn-cwrap">
        <div><div class="fn-k">Homozygous ${v.consClass==='lof'?'(candidate knockouts)':''} · ${carriersHom.length}</div>
          <div class="fn-chips">${homs||'<span class="muted">none</span>'}${carriersHom.length>60?` <span class="muted">+${carriersHom.length-60} more</span>`:''}</div></div>
        <div style="margin-top:8px"><div class="fn-k">Heterozygous · ${carriersHet.length}</div>
          <div class="fn-chips">${hets||'<span class="muted">none</span>'}${carriersHet.length>60?` <span class="muted">+${carriersHet.length-60} more</span>`:''}</div></div>
        <div class="fold-sendrow">
          <span class="fn-k" style="margin:0">Open in SNPVersity</span>
          ${send('hom','Homozygous carriers',carriersHom.length)}
          ${send('het','Heterozygous carriers',carriersHet.length)}
          ${send('all','All carriers',carriersHom.length+carriersHet.length)}
        </div>
      </div></td></tr>`;
  }

  /* ---------- context panel ---------- */
  /* does any variant in the table have carrier data? drives the "send all" button */
  function anyCarriers(){
    return (FD.variants||[]).some(v=>{ const c=carrierOf(v); return c && ((Number(c.hom)||0)+(Number(c.het)||0))>0; });
  }

  /* Genomic span of the current gene, for the SNPVersity handoff.
     Preferred source is geneFunction (already fetched); falls back to fields on the
     structure record, then to an on-demand Data.lookupGene(). Cached in FD.locus. */
  async function foldLocus(){
    if (FD.locus && FD.locus.chr != null) return FD.locus;
    const s = FD.struct || {};
    if (s.chr != null && s.start != null && s.end != null){
      FD.locus = { chr:s.chr, start:+s.start, end:+s.end, dataset:datasetId(FD.dataset) };
      return FD.locus;
    }
    if (typeof Data !== 'undefined' && typeof Data.lookupGene === 'function'){
      try {
        const g = await Data.lookupGene(FD.gene);
        if (g && g.chr != null){
          FD.locus = { chr:g.chr, start:+g.start, end:+g.end, dataset:datasetId(FD.dataset) };
          return FD.locus;
        }
      } catch (e){ console.error('SNPFold: gene lookup failed', e); }
    }
    return null;
  }

  function ctxHTML(){
    const v = FD.variants.find(x=>x.id===FD.selId);
    if (!v) return `<div class="ctx-empty">
        <div class="ei2">${ICONS.fold}</div>
        <div><b>Select a variant</b><div class="muted">Click a lollipop, a table row, or a residue in 3D to see its structural context.</div></div>
      </div>`;
    const c = ctxFor(v);
    const [tone, msg] = interpret(v, c);
    return `
      <div class="ctx-head">
        <div>
          <div class="ctx-v mono">${v.variant}</div>
          <div><span class="cons ${v.consClass}">${v.consequence}</span> ${v.priority?`<span class="prio ${v.priority.toLowerCase()}">${v.priority}</span>`:''}</div>
        </div>
        <button class="btn ghost" onclick="FOLD.focus()">Focus in 3D</button>
      </div>
      <div class="interp ${tone}">${msg}</div>
      <div class="ctx-grid">
        <div class="ck"><div class="kk">Residue</div><div class="vv mono">${finiteNumber(v.resi)==null?'—':`${v.ref||''}${v.resi}${v.consClass==='missense'&&v.alt?v.alt:''}`}</div></div>
        <div class="ck"><div class="kk">Domain</div><div class="vv">${c.domain?c.domain.name:'—'}${c.domain&&c.domain.pfam&&c.domain.pfam!=='region'?` <span class="muted mono">${c.domain.pfam}</span>`:''}</div></div>
        <div class="ck"><div class="kk">Local confidence</div><div class="vv">${c.plddt==null?'<span class="muted">outside model</span>':`<span class="plddt-chip" style="background:${plddtHex(c.plddt)};color:${c.plddt>=70?'#06294f':'#5c3a06'}">${c.plddt.toFixed(0)}</span> <span class="muted">${plddtBand(c.plddt)}</span>`}</div></div>
        <div class="ck"><div class="kk">Secondary structure</div><div class="vv"><span class="ss-chip ss-${c.ss}">${c.ssLabel}</span></div></div>
        <div class="ck"><div class="kk">AI scores</div><div class="vv mono">${aiScoreSummary(v)}</div></div>
      </div>
      <div class="ctx-actions">
        ${v.consClass==='missense' && v.resi ? `<button class="btn" onclick="FOLD.panEffect('${v.id}')">${ICONS.effect||ICONS.star} PanEffect</button>` : ''}
        ${(()=>{ const cr=carrierOf(v); const n=cr?((Number(cr.hom)||0)+(Number(cr.het)||0)):0;
          return n ? `<button class="btn" title="Preselect the ${n} accession${n>1?'s':''} carrying this allele in SNPVersity"
            onclick="FOLD.toVersity('${v.id}','all')">${ICONS.dna||''} Carriers to SNPVersity</button>` : ''; })()}
        <button class="btn" onclick="go('snpimpact')">${ICONS.star} SNPImpact</button>
        <button class="btn" onclick="go('snpcompare')">${ICONS.compare} Send to SNPCompare</button>
      </div>`;
  }

  /* =================== 3D VIEWER =================== */
  function ensure3Dmol(cb){
    if (window.$3Dmol) return cb(true);
    if (FD.libState==='failed') return cb(false);
    const done = ()=>cb(!!window.$3Dmol);
    if (FD.libState==='loading'){ const t=setInterval(()=>{ if(window.$3Dmol){clearInterval(t);cb(true);} },120);
      setTimeout(()=>{clearInterval(t); if(!window.$3Dmol) cb(false);},6000); return; }
    FD.libState='loading';
    const sc=document.createElement('script'); sc.src=THREEDMOL_URL;
    sc.onload=()=>{ FD.libState='ready'; done(); };
    sc.onerror=()=>{ FD.libState='failed'; cb(false); };
    document.head.appendChild(sc);
  }
  function buildViewer(){
    const el = document.getElementById('fold3d'); if (!el) return;
    el.innerHTML = `<div class="v-loading"><div class="spinner"></div><div>Loading 3D viewer…</div></div>`;
    ensure3Dmol(ok=>{
      if (!ok || !FD.pdb){
        el.innerHTML = `<div class="fold-3d-fallback">
          <div class="ei2">${ICONS.fold}</div>
          <b>3D viewer unavailable offline</b>
          <p>The interactive structure needs the 3Dmol.js library (loaded from a CDN). The protein
          browser, domains, and structural context above work without it. To enable 3D on a
          restricted network, self-host <span class="mono">3Dmol-min.js</span> and point SNPFold at it.</p>
        </div>`;
        return;
      }
      try{
        el.innerHTML='';
        FD.viewer = $3Dmol.createViewer(el, { backgroundColor:'white' });
        FD.viewer.addModel(FD.pdb, 'pdb');
        applyStyle();
        FD.viewer.zoomTo();
        FD.viewer.render();
        if (FD.selId) focusResidue(false);
      }catch(e){
        console.error('3Dmol error', e);
        el.innerHTML = `<div class="fold-3d-fallback"><b>Could not render structure</b><p class="mono">${String(e.message||e)}</p></div>`;
      }
    });
  }
  function applyStyle(){
    if (!FD.viewer) return;
    const v=FD.viewer; v.setStyle({}, {});
    if (FD.colorMode==='plddt'){
      v.setStyle({}, { cartoon:{ colorfunc:(a)=>plddtHex(a.b) } });
    } else if (FD.colorMode==='domain'){
      v.setStyle({}, { cartoon:{ color:'#dfe5ee' } });
      (FD.struct.domains||[]).filter(d=>d.kind==='domain').forEach((d,i)=>{
        v.setStyle({ resi:d.start+'-'+d.end }, { cartoon:{ color: ['#7fae7f','#6f97d6','#a98fd0','#d6a45b'][i%4] } });
      });
    } else {
      v.setStyle({}, { cartoon:{ color:'#d8dde6' } });
    }
    if (FD.showVar){
      FD.variants.forEach(vr=>{
        const resi = plausibleResidue(vr.resi);
        if (resi==null || !FD.struct || resi>FD.struct.length) return;
        const col = FD.colorMode==='impact' ? impactHex(vr.combined) : (CONS_FILL[vr.consClass]||'#2f5bbf');
        v.addStyle({ resi }, { stick:{ radius:.18 } });
        v.addStyle({ resi }, { sphere:{ scale: vr.id===FD.selId ? 0.7 : 0.45, color: col } });
      });
    }
    v.render();
  }
  function focusResidue(animate){
    const v=FD.variants.find(x=>x.id===FD.selId); if(!FD.viewer||!v) return;
    const resi = plausibleResidue(v.resi);
    if (resi==null || !FD.struct || resi>FD.struct.length) return;
    FD.viewer.removeAllLabels();
    FD.viewer.addLabel(v.variant, { position:{ resi }, backgroundColor:'#13264a', backgroundOpacity:.9,
      fontColor:'white', fontSize:11, borderThickness:0 }, { resi });
    applyStyle();
    FD.viewer.zoomTo({ resi }, animate?500:0);
    FD.viewer.render();
  }

  /* =================== handlers =================== */
  function refreshSelection(){
    const t=document.getElementById('foldTrack'); if(t) t.innerHTML=trackSVG();
    const c=document.getElementById('foldCtx'); if(c) c.innerHTML=ctxHTML();
    refreshTable();
  }
  /* re-renders header (sort arrows) and body (sorted rows) together */
  function refreshTable(){
    const h=document.getElementById('foldTableHead'); if(h) h.innerHTML=tableHeadHTML();
    const b=document.getElementById('foldTableBody'); if(b) b.innerHTML=sortedVariants().map(rowHTML).join('');
    if (typeof attachTT==='function') attachTT();
  }
  window.FOLD = {
    select(id){ FD.selId = (FD.selId===id?null:id); refreshSelection(); if(FD.selId) focusResidue(true); else if(FD.viewer){FD.viewer.removeAllLabels();applyStyle();FD.viewer.zoomTo();FD.viewer.render();} },
    carriers(id){ FD.openCarrier = (FD.openCarrier===id?null:id); refreshTable(); },

    /* Hand this gene's region + the accessions carrying an alternative allele over to
       SNPVersity.  id 'all' = every coding variant with carriers in the table.
       mode: 'hom' | 'het' | 'all'. */
    async toVersity(id, mode){
      mode = mode || 'all';
      const vs = id==='all'
        ? (FD.variants||[]).filter(v=>{ const c=carrierOf(v); return c && ((Number(c.hom)||0)+(Number(c.het)||0))>0; })
        : (FD.variants||[]).filter(v=>String(v.id)===String(id));
      if (!vs.length){ alert('No carrier data available for this variant.'); return; }

      const acc=new Set();
      vs.forEach(v=>{
        const c=carrierOf(v); if(!c) return;
        if(mode!=='het') (c.carriersHom||[]).forEach(a=>acc.add(a));
        if(mode!=='hom') (c.carriersHet||[]).forEach(a=>acc.add(a));
      });
      if(!acc.size){ alert('No carriers to send for this allele.'); return; }

      const locus = await foldLocus();
      if(!locus){
        alert('Could not determine the genomic coordinates for '+FD.gene+', so the region cannot be set in SNPVersity.');
        return;
      }
      const what = mode==='hom' ? 'homozygous carriers'
                 : mode==='het' ? 'heterozygous carriers'
                 : 'carriers of an alternative allele';
      const note = id==='all'
        ? `${what} across ${vs.length} coding variant${vs.length>1?'s':''}`
        : `${what} of ${vs[0].variant}`;
      const payload = {
        gene:FD.gene, chr:locus.chr, start:locus.start, end:locus.end,
        dataset:(locus.dataset!=null?locus.dataset:datasetId(FD.dataset)),
        accessions:[...acc], from:'SNPFold', note,
        allele: id==='all' ? null : vs[0].variant, mode,
        merge: (typeof Handoff!=='undefined') ? Handoff.mode() : 'replace'
      };
      if (typeof Handoff!=='undefined'){ Handoff.toVersity(payload); return; }
      if (typeof window.versityRequest === 'function'){ window.versityRequest(payload); return; }
      if (typeof S !== 'undefined' && S){          // older snpversity.js build
        S.pendingVersity = payload;
        if (typeof go === 'function') go('snpversity');
      }
    },
    /* Hand this gene model over to SNPFunction for its full functional dossier
       (domains, GO terms, variant burden, damaging-allele catalog). Prefers the
       shared goFunction(gene, dataset) helper; falls back to setting the
       inbound-request fields directly for older builds without it. */
    toFunction(){
      if (!FD.gene) return;
      const ds = datasetId(FD.dataset);
      if (typeof goFunction === 'function'){ goFunction(FD.gene, ds); return; }
      if (typeof S !== 'undefined' && S){
        S.functionGene = FD.gene;
        S.functionDataset = ds;
        if (typeof go === 'function') go('snpfunction');
      }
    },
    /* Click a column header to sort; click the same header again to reverse.
       Blank cells (—) always sort to the bottom, whichever direction is active. */
    sortBy(key){
      const col = foldCol(key); if (!col) return;
      if (FD.sort.key === key) FD.sort = { key, dir: FD.sort.dir === 'asc' ? 'desc' : 'asc' };
      else                     FD.sort = { key, dir: col.desc1 ? 'desc' : 'asc' };
      refreshTable();
    },
    panEffect(id){ panEffect(FD.variants.find(v=>v.id===id)); },
    color(m){ FD.colorMode=m; document.querySelectorAll('.fold-toolbar .seg-b').forEach(b=>b.classList.remove('on'));
      const map={plddt:0,domain:1,impact:2}; const btns=document.querySelectorAll('.fold-toolbar .seg-b'); if(btns[map[m]])btns[map[m]].classList.add('on'); applyStyle(); },
    toggleVar(on){ FD.showVar=on; applyStyle(); },
    reset(){ if(FD.viewer){ FD.viewer.zoomTo(); FD.viewer.render(); } },
    focus(){ focusResidue(true); },
    loadGene(){ const el=document.getElementById('foldGeneInput'); if(!el)return;
      const g=el.value.trim(); if(!g)return; FD.gene=g; FD.selId=null; FD.openCarrier=null; loadStructure(); },
    /* Best/AlphaFold2/Boltz2/ESMFold radio — reload the current gene's structure under
       the new preference if one is already loaded, otherwise just record the choice. */
    setModelPref(pref){
      FD.modelPref = pref;
      if (FD.loaded) loadStructure();
    },
    setDataset(dataset){ FD.dataset=dataset; if(typeof S!=='undefined'&&S)S.dataset=datasetId(dataset); FD.sec=hasSecondaryScores(dataset);
      if(FD.loaded) loadStructure(); else renderLanding(); },
    /* compact dataset chooser on this page — select which dataset to use directly */
    pickDataset(id){
      const ds = foldDatasets().find(d=>String(datasetId(d))===String(id));
      const val = ds ? datasetId(ds) : id;
      if (isCurrentDataset(ds || val)) { syncDatasetChooser(); return; }   // no change
      FD.dataset = val;
      if (typeof S !== 'undefined' && S) S.dataset = val;   // keep the whole app in sync
      FD.sec = hasSecondaryScores(val);
      /* Selecting a dataset only records the choice + moves the highlight. The page
         is NOT reloaded here — the new dataset is applied when the user presses
         "Load structure" (FOLD.loadGene), which reads the current FD.dataset. */
      syncDatasetChooser();
      if (typeof Handoff!=='undefined') Handoff.sync(FD.root);
    },
  };

  /* carrier-chip + expandable-row styles — shared look with SNPFunction */
  function injectCSS(){
    if (document.getElementById('snpfold-carrier-css')) return;
    const s=document.createElement('style'); s.id='snpfold-carrier-css';
    s.textContent=`
      /* sortable variant-table headers
         Hover/active states are translucent overlays rather than opaque fills, so they
         tint whatever the .vcf header background already is and never fight its text
         colour. On a light header swap the two rgba(255,255,255,…) values for
         rgba(0,0,0,.05) / rgba(0,0,0,.08) to darken instead of lighten. */
      .fold-table th.fold-th{cursor:pointer;user-select:none;white-space:nowrap;position:relative;
        transition:background-color .12s}
      .fold-table th.fold-th:hover{background-color:rgba(55,55,55,.33)}
      .fold-table th.fold-th .fold-th-in{display:inline-flex;align-items:center;gap:5px}
      .fold-table th.fold-th.num .fold-th-in{justify-content:flex-end}
      .fold-table th.fold-th .fold-ar{font-size:9px;line-height:1;opacity:.35;transition:opacity .12s}
      .fold-table th.fold-th:hover .fold-ar{opacity:.7}
      .fold-table th.fold-th.sorted{background-color:rgba(100,100,100,.75)}
      .fold-table th.fold-th.sorted .fold-ar{opacity:1}
      /* handoff notice (e.g. a requested variant is not in this gene's table) */
      .fold-notice{display:flex;align-items:center;gap:10px;background:#fffaf0;border:1px solid #f0dcb4;
        color:#8a6d1e;border-radius:9px;padding:8px 12px;font-size:12.5px;margin-bottom:12px}
      .fold-notice button{margin-left:auto;border:0;background:none;color:inherit;font-size:16px;
        line-height:1;cursor:pointer;padding:0 2px}
      /* Persistent warning when normalized protein variants extend beyond the loaded model. */
      .fold-trunc-warning{display:flex;align-items:flex-start;gap:11px;background:#fff8e8;
        border:1px solid #e9c86c;border-left:4px solid #d99a16;color:#6f5012;border-radius:9px;
        padding:11px 14px;margin:12px 0 18px;font-size:12.5px;line-height:1.5}
      .fold-trunc-warning b{color:#62430a}
      .fold-trunc-icon{font-size:18px;line-height:1.2;flex:0 0 auto}
      .snpfold-root .prio.top{background:#fde2df;color:#9f2018;border-color:#efb8b1}
      .fold-cbtn{display:inline-flex;align-items:center;border:1px solid var(--line);border-radius:7px;overflow:hidden;background:#fff;cursor:pointer;padding:0}
      .fold-cbtn:hover{border-color:#c3cee0}
      .fold-cbtn.on{border-color:#9db4dd;box-shadow:0 0 0 2px rgba(47,106,208,.12)}
      .fold-cbtn .cc{font-family:var(--mono);font-size:11px;font-weight:600;padding:2px 7px}
      .fold-cbtn .cc.hom{background:#fdecea;color:#8f281c}
      .fold-cbtn .cc.het{background:#eef4ff;color:#274b8f}
      .fn-carriers td{background:#fbfcfe;border-bottom:1px solid var(--line)}
      .fn-cwrap{padding:10px 12px} .fn-cwrap.muted{color:var(--muted)}
      .fn-chips{display:flex;flex-wrap:wrap;gap:5px;margin-top:5px}
      .fn-k{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:4px}
      .carrier{font-family:var(--mono);font-size:11px;padding:2px 7px;border-radius:6px;border:1px solid var(--line)}
      .carrier.hom{background:#fdecea;border-color:#f0c4bd;color:#8f281c}
      .carrier.het{background:#eef4ff;border-color:#cfe0ff;color:#274b8f}
      .fold-sendrow{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-top:12px;
        border-top:1px solid var(--line);padding-top:10px}
      /* handoff mode — sits directly above the table whose buttons it governs */
      .fn-handoff{display:flex;align-items:baseline;gap:12px;flex-wrap:wrap;
        background:#f7faff;border:1px solid #dce6f5;border-radius:9px;
        padding:8px 12px;margin:12px 0}
      .fn-handoff-k{font-size:10.5px;font-weight:600;color:var(--muted);
        text-transform:uppercase;letter-spacing:.4px;flex:0 0 auto}
      .btn.tiny{font-size:11px;padding:4px 9px;border-radius:7px;line-height:1.3;white-space:nowrap}
      td.fold-send{text-align:right;padding-right:10px;white-space:nowrap}
      .fold-table th.fold-send-th{width:1%}
      /* compact dataset chooser (shares Data.datasets() with SNPVersity) */
      .fold-ds-card{padding:12px 14px}
      .fold-ds-head{font-size:10.5px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.4px;margin-bottom:9px}
      .fold-ds-grid{display:flex;flex-wrap:wrap;gap:10px}
      .fold-ds{display:flex;align-items:center;gap:9px;text-align:left;cursor:pointer;background:#fff;
        border:1px solid var(--line);border-radius:10px;padding:9px 13px;min-width:210px;transition:border-color .12s,box-shadow .12s}
      .fold-ds:hover{border-color:#c3cee0}
      .fold-ds.sel{border-color:#9db4dd;box-shadow:0 0 0 2px rgba(47,106,208,.12);background:#f7faff}
      .fold-ds-dot{flex:0 0 auto;width:10px;height:10px;border-radius:50%;border:2px solid #c3cee0;background:#fff}
      .fold-ds.sel .fold-ds-dot{border-color:#2f6ad0;background:#2f6ad0;box-shadow:inset 0 0 0 2px #fff}
      .fold-ds-txt{display:flex;flex-direction:column;line-height:1.25;min-width:0}
      .fold-ds-name{font-weight:600;font-size:13px;color:var(--ink)}
      .fold-ds-sub{font-size:11.5px;color:var(--muted)}
      /* structure-source badge (AlphaFold2 vs Boltz2 vs ESMFold) shown in the gene context line */
      .fold-src-badge{display:inline-flex;align-items:center;gap:5px;font-size:11px;font-weight:600;
        padding:2px 9px 2px 7px;border-radius:999px;letter-spacing:.2px;vertical-align:middle}
      .fold-src-dot{width:7px;height:7px;border-radius:50%;flex:0 0 auto}
      .fold-src-badge.alphafold{background:#eaf1fd;color:#1a4fa0}
      .fold-src-badge.alphafold .fold-src-dot{background:#2f6ad0}
      .fold-src-badge.boltz{background:#fdeeea;color:#a03d1a}
      .fold-src-badge.boltz .fold-src-dot{background:#d0682f}
      .fold-src-badge.esmfold{background:#eafaf0;color:#0f7a45}
      .fold-src-badge.esmfold .fold-src-dot{background:#22a35f}
      /* Best / AlphaFold2 / Boltz2 / ESMFold model picker in the search bar */
      .fold-model-pref{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-left:6px}
      .fold-model-opt{display:inline-flex;align-items:center;gap:4px;font-size:12.5px;color:var(--ink);
        cursor:pointer;white-space:nowrap}
      .fold-model-opt input{accent-color:#2f6ad0;cursor:pointer}`;
    document.head.appendChild(s);
  }

  window.addEventListener('resize', ()=>{ if(FD.viewer) try{ FD.viewer.resize(); }catch(e){} });

  SNPTools.register('snpfold', { render });
})();
