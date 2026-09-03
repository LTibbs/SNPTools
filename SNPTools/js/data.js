/* =====================================================================
 *  data.js — the DATA LAYER (LIVE / real HDF5 build).
 *
 *  queryVariants() no longer fabricates demo data: it POSTs the region +
 *  accession list to processForm.php, which runs h5_to_vcf.py against the
 *  real .h5 store, writes a VCF, and returns its path. We then fetch that
 *  VCF and parse it into the exact row shape the tools already expect.
 *
 *  Accession IDs come from accessions.real.js (window.SNP_REAL_ACCESSIONS),
 *  which holds the actual column names inside the .h5 files, so a selection
 *  maps to real HDF5 columns.
 *
 *  PFAM / protein-domain data is not in the HDF5 yet — the "Domain" column
 *  is filled with "N/A" until that separate data structure is built.
 *
 *  SNPImpact / SNPFold still use demo generators (they are not backed by
 *  these .h5 files); those bodies are left untouched below.
 *
 *  Depends on rnd() and pick() from core.js.
 * ===================================================================== */
const Data = (function () {

  /* =============================================================
   *  BACKEND CONFIG — edit these two lines if your paths differ.
   * ============================================================= */
  const CFG = {
    endpoint : 'processForm.php',  // relative to index.html
    vcfDir   : 'vcf/',             // web-served, writable; MUST match processForm.php
    geneEndpoint : 'lookupGeneModel.php',  // gene model -> coordinates
    structDir : 'data/structures/',  // per-gene structure-<gene>.js files (SNPFold)
    domainsDir : 'data/domains/by_chr/',              // per-chromosome Pfam files: <chr>.json (preferred)
    domainsUrl : 'data/domains/domains.by_chr.json',  // combined file (fallback if per-chrom absent)
    domainsGeneUrl : 'data/domains/domains.by_gene.json',  // gene -> canonical protein domains (SNPImpact detail)
    geneModelsDir : 'data/genemodels/by_chr/',        // per-chromosome exon/CDS structure (SNPImpact gene model)
    // Hard span backstop — a region wider than this is download-only regardless
    // of how few accessions are selected. The real gate is the variants ×
    // accessions "table work" budget below (TABLE_WORK_MAX / TABLE_ROWS_MAX);
    // this just catches an absurd region.
    tableMaxSpan : 20_000_000,
  };

  /* ---------------- datasets (UI cards) ----------------
   * `id`     -> sent to processForm.php as dataSet
   * `family` -> which real accession list to show + which .h5 family
   */
  const DATASETS = [
    {id:'mgdb2026_hq', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Quality',   ref:'B73 v5', acc:'2,710', sites:'98M',
     filters:['MQ ≥ 30','Coverage ≥ 50%','LD max R² > 0.5'], het:true,  indel:true,  impute:false},
    {id:'mgdb2026_hc', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Coverage',  ref:'B73 v5', acc:'2,710', sites:'290M',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:true,  indel:true,  impute:false},
  ];

  /*const DATASETS = [
    {id:'mgdb2026_hq', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Quality',   ref:'B73 v5', acc:'2,710', sites:'98M',
     filters:['MQ ≥ 30','Coverage ≥ 50%','LD max R² > 0.5'], het:true,  indel:true,  impute:false},
    {id:'mgdb2026_hc', family:'mgdb2026',     name:'MaizeGDB 2026', sub:'High Coverage',  ref:'B73 v5', acc:'2,710', sites:'290M',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:true,  indel:true,  impute:false},
    {id:'schnable2023',family:'schnable2023', name:'Schnable 2023', sub:'Imputed markers',ref:'B73 v5', acc:'1,515', sites:'12M',
     filters:['Imputed'], het:false, indel:false, impute:true},
    {id:'nam2021',     family:'nam2021',      name:'NAM 2021',      sub:'Founder panel',  ref:'B73 v5', acc:'27',    sites:'78M',
     filters:['MQ ≥ 30','Founder panel'], het:true,  indel:true,  impute:false},
    {id:'mgdb2024_hq', family:'mgdb2024',     name:'MaizeGDB 2024', sub:'High Quality',   ref:'B73 v5', acc:'1,498', sites:'83M',
     filters:['MQ ≥ 30','Coverage ≥ 50%'], het:true,  indel:true,  impute:false},
  ];
  */

  /* colors + a friendly label per family, for the accession picker group */
  const FAMILY_META = {
    mgdb2026:     {name:'MaizeGDB 2026', color:'#2563eb'},
    mgdb2024:     {name:'MaizeGDB 2024', color:'#cf8a12'},
    schnable2023: {name:'Schnable 2023', color:'#1f8a4c'},
    nam2021:      {name:'NAM 2021',      color:'#7c3aed'},
  };

  /* ---------------- accession catalog (projects -> groups -> accessions) ----------------
   * Source: window.SNP_CATALOG (compiled from data/accessions.tsv + data/projects.tsv).
   * Falls back to the older flat window.SNP_REAL_ACCESSIONS if present.
   */
  const CATALOG = (typeof window !== 'undefined' && window.SNP_CATALOG) || null;
  const REAL    = (typeof window !== 'undefined' && window.SNP_REAL_ACCESSIONS) || {};

  function familyOf(datasetId){
    const d = DATASETS.find(x => x.id === datasetId);
    return d ? d.family : 'mgdb2026';
  }
  // MaizeGDB 2026 is the only family carrying the second-generation language-model
  // scores (PlantCAD2 / ESM2). Tools call this to show those columns conditionally.
  function hasSecondaryScores(datasetId){ return familyOf(datasetId) === 'mgdb2026'; }
  function famNode(datasetId){
    return CATALOG ? CATALOG.families[familyOf(datasetId)] : null;
  }

  // Projects (bioproject sections) with metadata + groups, for the picker.
  function projectsFor(datasetId){
    const fam = famNode(datasetId);
    if (fam) return fam.projects;
    // legacy fallback: single synthetic project from the flat list
    const list = REAL[familyOf(datasetId)] || [];
    return [{id:familyOf(datasetId), title:familyOf(datasetId), bioprojects:[], color:'#2563eb',
             count:list.length, namFounders:[],
             groups:[{name:null, accessions:list}]}];
  }

  // Flat accession list (one object per accession) for search, chips, table headers.
  const _accCache = {};
  function accessionsFor(datasetId){
    const fam = familyOf(datasetId);
    if (_accCache[fam]) return _accCache[fam];
    const node = famNode(datasetId);
    let out = [];
    if (node){
      node.projects.forEach(p => p.groups.forEach(g => g.accessions.forEach(a => {
        out.push({id:a.id, run:a.run, founder:a.founder, rep:a.rep, reps:a.reps,
                  label:a.label, group:g.name, namFounder:a.namFounder,
                  proj:p.id, projColor:p.color, projTitle:p.title});
      })));
    } else {
      out = (REAL[fam] || []).map(a => ({...a}));
    }
    _accCache[fam] = out;
    return out;
  }

  function namFoundersFor(datasetId){
    const fam = famNode(datasetId);
    return fam ? (fam.namFounders || []) : [];
  }

  // Default selection: for a dataset with tagged NAM founders (2026), preselect
  // one accession per NAM founder; otherwise one per the first 12 founders.
  function defaultSelectionFor(datasetId){
    const list = accessionsFor(datasetId);
    const nam = namFoundersFor(datasetId);
    if (nam.length){
      const pick = {}, ids = [];
      for (const a of list){
        if (a.namFounder && !pick[a.namFounder]){ pick[a.namFounder] = 1; ids.push(a.id); }
      }
      return ids;
    }
    const ids = [], seen = new Set();
    for (const a of list){
      if (seen.has(a.founder)) continue;
      seen.add(a.founder); ids.push(a.id);
      if (ids.length >= 12) break;
    }
    return ids;
  }

  // union index for accessionById()
  let _byId = null;
  function idIndex(){
    if (_byId) return _byId;
    _byId = new Map();
    ['mgdb2026','mgdb2024','schnable2023','nam2021'].forEach(fam => {
      const dsid = (DATASETS.find(d => d.family === fam) || {}).id;
      if (dsid) accessionsFor(dsid).forEach(a => { if (!_byId.has(a.id)) _byId.set(a.id, a); });
    });
    return _byId;
  }
  function accessionById(id){ return idIndex().get(id) || null; }

  /* ---------------- genome geometry ---------------- */
  const GENE_MODELS = {
    'Zm00001eb374090':{chr:'chr8', start:163450112, end:163454880},
    'Zm00001eb067740':{chr:'chr2', start:21008440,  end:21013990},
    'Zm00001eb404760':{chr:'chr10',start:9821400,   end:9826110},
    'Zm00001eb404740':{chr:'chr10',start:9788220,   end:9794010},
    'Zm00001eb233650':{chr:'chr5', start:8841220,   end:8849510},
    'Zm00001eb313510':{chr:'chr7', start:174221000, end:174229800},
  };
  const CHR_LEN = {chr1:308452471,chr2:243675191,chr3:238017767,chr4:250330460,chr5:226353449,
    chr6:181357234,chr7:185808916,chr8:182411202,chr9:163004744,chr10:152435371};
  const CENTRO = {chr10:.34};

  /* ---------------- gene model -> coordinates (live) ----------------
   * Resolves a B73 v5 gene model ID to its interval via lookupGeneModel.php,
   * which reads the serialized GFF store on the server. Returns
   *   {id, chr, start, end}  or  null when the ID isn't found.
   */
  const EXAMPLE_GENES = [
    'Zm00001eb374090','Zm00001eb067740','Zm00001eb374230',
    'Zm00001eb056510','Zm00001eb233650','Zm00001eb313510',
  ];
  async function lookupGene(id){
    id = (id || '').trim();
    if (!id) return null;
    const url = `${CFG.geneEndpoint}?geneModelId=${encodeURIComponent(id)}`;
    const resp = await fetch(url, {cache:'no-store'});
    if (!resp.ok) throw new Error('Gene lookup failed (HTTP ' + resp.status + ')');
    const raw = await resp.text();
    let d;
    try { d = JSON.parse(raw); }
    catch (e) { throw new Error('lookupGeneModel.php did not return JSON:\n' + raw.slice(0, 600)); }
    // The PHP returns id:'empty' (and 0/0) when the gene isn't in the store.
    if (!d || d.id === 'empty' || d.chromosome == null) return null;
    const start = parseInt(d.start, 10), end = parseInt(d.end, 10);
    if (!Number.isFinite(start) || !Number.isFinite(end) || (start === 0 && end === 0)) return null;
    return {id, chr: d.chromosome, start, end};
  }

  /* =============================================================
   *  SNPVERSITY — LIVE query.
   *  region + accession ids -> VCF (via processForm.php) -> rows
   * ============================================================= */
  const SEVERITY = {HIGH:3, MODERATE:2, LOW:1, MODIFIER:0};

  function uniqueOutName(lo, hi){
    const ts = Date.now();
    const rnd = Math.random().toString().slice(2, 11);
    return `${CFG.vcfDir}snpv_${ts}_${rnd}_${lo}_${hi}.vcf.gz`;
  }

  /* processForm.php returns a gzip-compressed VCF (.vcf.gz). Inflate it in
     the browser; fall back to treating the bytes as plain text if they are
     not gzip-framed (a proxy that already decompressed, or a plain .vcf
     from some other caller). */
  async function readVcfResponse(resp){
    const buf = new Uint8Array(await resp.arrayBuffer());
    const gzipped = buf.length > 2 && buf[0] === 0x1f && buf[1] === 0x8b;
    if (gzipped && typeof DecompressionStream === 'function'){
      const stream = new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'));
      return await new Response(stream).text();
    }
    return new TextDecoder().decode(buf);
  }

  /* Parse the whole INFO column (';'-separated key=value tags) once per row
     into a plain object. This replaced ~20 calls/row to a helper that
     compiled a fresh RegExp each time — millions of compilations for a wide
     region. A missing key reads back as undefined, which numOrNull /
     firstNum / cleanTok already treat the same as null. */
  function parseInfo(infoStr){
    const o = {};
    if (!infoStr) return o;
    const parts = infoStr.split(';');
    for (let i = 0; i < parts.length; i++){
      const p = parts[i];
      const eq = p.indexOf('=');
      if (eq === -1) o[p] = ''; else o[p.slice(0, eq)] = p.slice(eq + 1);
    }
    return o;
  }
  const cleanTok = s => {
    if (s == null) return '';
    const tok = String(s).split(/[;,]+/)[0].replace(/_/g, ' ').trim();
    return tok === '.' ? '' : tok;
  };
  function numOrNull(v){
    if (v == null || v === '.' || v === '') return null;
    const f = parseFloat(v);
    return Number.isNaN(f) ? null : f;
  }
  function firstNum(v){                       // first of a possibly comma-listed value
    if (v == null) return null;
    return numOrNull(String(v).split(/[;,]+/)[0]);
  }

  /* A region's genotype matrix is variants × accessions — by far the biggest
     thing in memory. Store each call as a 1-byte code in an Int8Array rather
     than a string:  0 = ref/ref, 1 = het, 2 = alt/alt, 3 = missing (null).
     ~8x smaller than an Array of strings, and the downstream dosage() calls
     become the identity. GT_CODE is the fast path — a direct lookup on the
     raw VCF cell, no allocation. */
  const GT_CODE = Object.create(null);
  GT_CODE['0/0'] = 0; GT_CODE['0|0'] = 0;
  GT_CODE['1/0'] = 1; GT_CODE['0/1'] = 1; GT_CODE['1|0'] = 1; GT_CODE['0|1'] = 1;
  GT_CODE['1/1'] = 2; GT_CODE['1|1'] = 2;
  GT_CODE['./.'] = 3; GT_CODE['.|.'] = 3; GT_CODE['.'] = 3;
  function gtCode(cell){
    if (cell == null) return 3;
    const q = GT_CODE[cell];
    if (q !== undefined) return q;
    // rare: FORMAT subfields (0/1:...) or non-0/1 allele indices
    let s = cell;
    const c = s.indexOf(':'); if (c !== -1) s = s.slice(0, c);
    const q2 = GT_CODE[s]; if (q2 !== undefined) return q2;
    const m = s.split(/[/|]/);
    if (m.length < 2 || m[0] === '.' || m[1] === '.' || m[0] === '' || m[1] === '') return 3;
    return (m[0] !== '0' ? 1 : 0) + (m[1] !== '0' ? 1 : 0);
  }

  function parseVcf(text, ids){
    const rows = [];
    if (!text) return {rows, sampleCols:{}, header:[]};
    const sampleCols = {};      // sampleName -> column index in the row
    const nIds = ids.length;
    const N = text.length;
    /* Walk lines by index rather than text.split('\n') — the split would
       hold a second full copy of the payload as an array of line strings. */
    let p = 0;
    while (p < N){
      let nl = text.indexOf('\n', p);
      if (nl === -1) nl = N;
      let end = nl;
      if (end > p && text.charCodeAt(end - 1) === 13) end--;   // strip trailing \r
      if (end === p){ p = nl + 1; continue; }                  // blank line
      if (text.charCodeAt(p) === 35){                          // '#'
        if (text.charCodeAt(p + 1) !== 35){                    // '#CHROM' header row, not '##meta'
          const cols = text.slice(p + 1, end).split('\t');
          for (let i = 9; i < cols.length; i++) sampleCols[cols[i]] = i;
        }
        p = nl + 1; continue;
      }
      const t = text.slice(p, end).split('\t');
      p = nl + 1;
      if (t.length < 8) continue;
      const II = parseInfo(t[7] || '');

      // impact: take the most severe when a site lists several
      let impact = 'MODIFIER';
      const effTokens = (II.EFFECT || '').split(/[;,]+/);
      let best = -1;
      for (let k = 0; k < effTokens.length; k++){
        const key = effTokens[k].trim().toUpperCase();
        if (key in SEVERITY && SEVERITY[key] > best){ best = SEVERITY[key]; impact = key; }
      }

      // genotypes in the exact order the caller selected, as 1-byte codes
      const gts = new Int8Array(nIds);
      for (let k = 0; k < nIds; k++){
        const ci = sampleCols[ids[k]];
        gts[k] = (ci == null) ? 3 : gtCode(t[ci]);
      }

      const mq = firstNum(II.MQ), cvp = firstNum(II.CVP), pos = parseInt(t[1], 10);
      rows.push({
        pos,
        ref:    t[3],
        alt:    t[4],
        gene:   cleanTok(II.GENEMODEL) || '—',
        effect: cleanTok(II.TYPE) || 'intergenic',
        impact,
        sub:    cleanTok(II.SUB),
        domain: domainAt(t[0], pos),     // Pfam domain covering this position (or '—')
        mq:     (mq != null) ? Math.round(mq) : 'N/A',
        comp:   (cvp != null) ? cvp : 'N/A',
        r2:     firstNum(II.MAXR2),
        maf:    (firstNum(II.MAF) != null) ? firstNum(II.MAF) : 0,
        // 2026 uses plantcad1/2 + ESM1/2/3; older projects (2024/Schnable/NAM)
        // use a single DNA_SCORE (PlantCaduceus) and AA_SCORE (ESM1b) -> map to col 1.
        pc1:    numOrNull(II.plantcad1_score != null ? II.plantcad1_score : II.DNA_SCORE),
        pc2:    numOrNull(II.plantcad2_score),
        esm1:   numOrNull(II.ESM1_score != null ? II.ESM1_score : II.AA_SCORE),
        esm2:   numOrNull(II.ESM2_score),
        esm3:   numOrNull(II.ESM3_score),
        gts,
      });
    }
    return {rows, sampleCols};
  }

  // Build the accession objects the table header needs, in selection order.
  function accsFor(datasetId, ids){
    const map = new Map(accessionsFor(datasetId).map(a => [a.id, a]));
    return ids.map(id => map.get(id) || {id, run:id, founder:id, proj:familyOf(datasetId), projColor:'#8a94a6'});
  }

  // Order the submitted accessions by project (catalog order), then by accession
  // name — so each project's color renders as one contiguous block in the table.
  function sortIds(datasetId, ids){
    const order = {};
    projectsFor(datasetId).forEach((p, i) => { order[p.id] = i; });
    const map = new Map(accessionsFor(datasetId).map(a => [a.id, a]));
    const keyOf = id => {
      const a = map.get(id);
      const pi = a && (order[a.proj] != null) ? order[a.proj] : 9999;
      return {pi, name: a ? a.id : id};
    };
    return ids.slice().sort((x, y) => {
      const kx = keyOf(x), ky = keyOf(y);
      if (kx.pi !== ky.pi) return kx.pi - ky.pi;
      return kx.name.localeCompare(ky.name, undefined, {numeric:true, sensitivity:'base'});
    });
  }

  /* Beyond this much table work (variants × accessions, or raw variant count)
     the in-browser parse + render would freeze the tab for many seconds, so the
     result is download-only — no table. Chosen so a typical selection (tens of
     accessions, a few Mb) always renders, but a full panel or a genome-arm-wide
     region does not. */
  const TABLE_WORK_MAX = 4e7;
  const TABLE_ROWS_MAX  = 400000;

  /* genome-wide variant density per bp for a dataset, from its manifest `sites`
     count ('98M' etc.) over the B73 v5 assembly length — used only to *predict*
     in the run bar whether a query will come back as a table or a VCF download
     (the real decision uses the exact server count). Regional density varies,
     so this is a rough guide. */
  const GENOME_LEN = Object.keys(CHR_LEN).reduce((s, c) => s + CHR_LEN[c], 0);
  function siteCount(s){
    const m = String(s || '').match(/([\d.]+)\s*([KMG]?)/i);
    if (!m) return 0;
    return parseFloat(m[1]) * ({K:1e3, M:1e6, G:1e9}[m[2].toUpperCase()] || 1);
  }
  /* estimateResult(datasetId, spanBp, nAccessions) -> {estVariants, willDownload} */
  function estimateResult(datasetId, spanBp, nAccessions){
    const d = DATASETS.find(x => x.id === datasetId);
    const density = d ? siteCount(d.sites) / GENOME_LEN : 0;
    const estVariants = Math.round(Math.max(0, spanBp) * density);
    const willDownload = spanBp > CFG.tableMaxSpan
      || estVariants > TABLE_ROWS_MAX
      || estVariants * Math.max(1, nAccessions) > TABLE_WORK_MAX;
    return {estVariants, willDownload};
  }

  /**
   * queryVariants(dataset, chr, lo, hi, ids, opts) -> Promise<{rows, accs, chr, vcfUrl, span, wide, empty, variants}>
   * `wide` is true when the result is too big for the table — the caller should
   * offer the VCF download only. `opts.forceTable` skips the size gate (used by
   * the internal gene-scoped callers, whose regions are always tiny).
   */
  async function queryVariants(dataset, chr, lo, hi, ids, opts){
    opts = opts || {};
    const forceTable = !!opts.forceTable;
    ids = sortIds(dataset, ids);      // group by project, then accession name
    const accs = accsFor(dataset, ids);
    const span = Math.max(hi - lo, 0);
    const outName = uniqueOutName(lo, hi);

    const body = new URLSearchParams({
      start: String(lo), end: String(hi), chr: chr, dataSet: dataset,
      genotypes: JSON.stringify(ids), outName: outName,
    });
    const resp = await fetch(CFG.endpoint, {
      method: 'POST',
      headers: {'Content-Type': 'application/x-www-form-urlencoded'},
      body,
    });
    const raw = await resp.text();
    let json;
    try { json = JSON.parse(raw); }
    catch (e) { throw new Error('processForm.php did not return JSON:\n' + raw.slice(0, 1000)); }

    // A genuinely empty interval is a valid result, not a failure.
    if (json.status === 'empty'){
      return {rows:[], accs, chr, vcfUrl:null, span, wide:false, empty:true};
    }
    if (!resp.ok || json.status !== 'success'){
      if (json.output)  console.error('[processForm.php] h5_to_vcf.py output:\n' + json.output);
      if (json.command) console.error('[processForm.php] command:\n' + json.command);
      const err = new Error(json.message || ('Request failed (HTTP ' + resp.status + ')'));
      err.detail = {command: json.command, output: json.output};
      throw err;
    }

    const vcfUrl = json.outFile || outName;
    const variants = (json.variants != null && isFinite(+json.variants)) ? +json.variants : null;

    // Too big for an in-browser table — download only. The server tells us the
    // exact variant count, so this decision happens *before* the fetch + parse.
    const tooBig = span > CFG.tableMaxSpan
      || (variants != null && (variants > TABLE_ROWS_MAX || variants * ids.length > TABLE_WORK_MAX));
    if (tooBig && !forceTable){
      return {rows:[], accs, chr, vcfUrl, span, variants, wide:true, empty:false};
    }

    // Fetch the VCF, retrying briefly — right after the server writes it, a
    // same-machine dev server can 404 it for a moment (Windows file handle) or
    // lose the race against a user-triggered download of the same URL.
    let vcfResp = null;
    for (let attempt = 0; attempt < 4; attempt++){
      if (attempt) await new Promise(r => setTimeout(r, 250 * attempt));
      try { vcfResp = await fetch(vcfUrl, {cache:'no-store'}); } catch (e) { vcfResp = null; }
      if (vcfResp && vcfResp.ok) break;
    }
    if (!vcfResp || !vcfResp.ok){
      const err = new Error('The VCF was built but could not be retrieved (' + (vcfResp ? 'HTTP ' + vcfResp.status : 'network error') + '). It is still available to download.');
      err.detail = {vcfUrl};
      throw err;
    }
    const vcfText = await readVcfResponse(vcfResp);
    await ensureDomains(chr);                    // load just this chromosome's Pfam file (cached)
    const {rows} = parseVcf(vcfText, ids);
    return {rows, accs, chr, vcfUrl, span, variants, wide:false, empty: rows.length === 0};
  }

  /* =============================================================
   *  SNPImpact query  (DEMO — not backed by these .h5 files)
   * ============================================================= */
  const BASES = ['A','C','G','T'];
  const CONSEQ = [
    {t:'Loss-of-function', cls:'lof',      base:-8.0},
    {t:'Loss-of-domain',  cls:'lod',      base:-6.0},
    {t:'Splice',          cls:'splice',   base:-3.0},
    {t:'Missense',        cls:'missense', base:-1.4},
    {t:'In-frame deletion',cls:'indel',   base:-0.6},
    {t:'Synonymous',      cls:'syn',      base: 0.3},
  ];
  const DOM_NAMES = ['Kinase domain','NB-ARC','bZIP','NAC domain','WRKY','DNA-binding domain','PPR repeat','F-box'];
  function priorityFromScore(s){ return s<=-7 ? 'TOP' : s<=-4 ? 'HIGH' : s<=-1 ? 'MODERATE' : 'LOW'; }

  function queryImpact(opts){
    opts = opts || {};
    const out = [];
    for (let i=0; i<46; i++){
      const c = pick(CONSEQ);
      const plantcad = +(c.base + rnd(-2,2)).toFixed(1);
      const esm      = +(c.base*0.72 + rnd(-1.5,1.5)).toFixed(1);
      const combined = +((plantcad + esm)/2).toFixed(2);
      const hasDom   = (c.cls==='lod' || c.cls==='missense' || Math.random()<.35);
      const aa = 90 + Math.floor(Math.random()*520);
      const exons = 4 + Math.floor(Math.random()*4);
      const affectedExon = 1 + Math.floor(Math.random()*exons);
      out.push({
        id:'v'+i,
        gene:'Zm00001eb'+(100000+Math.floor(Math.random()*899999)),
        variant: c.cls==='lof'      ? 'p.'+pick(['W','Q','R','E','K'])+aa+'*'
               : c.cls==='missense' ? 'p.'+pick(['A','G','R','D','V'])+aa+pick(['R','K','L','P','S'])
               : c.cls==='lod'      ? 'Δ Exon '+affectedExon
               : c.cls==='splice'   ? 'splice-site'
               : c.cls==='indel'    ? 'deletion'
               :                      'c.'+aa+pick(BASES)+'>'+pick(BASES),
        consequence:c.t, consClass:c.cls,
        domain: hasDom ? pick(DOM_NAMES) : '—',
        plantcad, esm, combined,
        priority: priorityFromScore(combined),
        percentile: Math.max(1, Math.min(99, Math.round(50 - combined*5 + rnd(-4,4)))),
        protLen: 280 + Math.floor(Math.random()*520),
        exons, affectedExon, aa,
      });
    }
    out.sort((a,b) => a.combined - b.combined);
    return out;
  }

  /* =============================================================
   *  SNPFold — protein structure + coding variants (DEMO curated).
   * ============================================================= */
  function structureFor(gene){ return (window.SNPFOLD_STRUCT||{})[gene] || null; }
  function pdbFor(gene){ return (window.SNPFOLD_PDB||{})[gene] || null; }
  /* Lazily load js/structures/structure-<gene>.js on demand. Resolves when the
     gene's model is available; rejects if there's no file for it. */
  function ensureStructure(gene){
    gene = (gene||'').trim();
    if (!gene) return Promise.reject(new Error('no gene'));
    if ((window.SNPFOLD_STRUCT||{})[gene] || (window.SNPFOLD_PDB||{})[gene]) return Promise.resolve(true);
    return new Promise((resolve, reject)=>{
      const s = document.createElement('script');
      s.src = `${CFG.structDir}structure-${encodeURIComponent(gene)}.js`;
      s.onload  = ()=> ((window.SNPFOLD_STRUCT||{})[gene] || (window.SNPFOLD_PDB||{})[gene])
        ? resolve(true) : reject(new Error('no model data for '+gene));
      s.onerror = ()=> reject(new Error('no structure file for '+gene));
      document.head.appendChild(s);
    });
  }

  /* ---- Pfam domains by genomic position (SNPVersity Domain column / SNPImpact) ----
     Loads ONE chromosome's file on demand (data/domains/by_chr/<chr>.json), cached
     per chromosome. Falls back to a combined domains.by_chr.json if per-chrom is absent,
     and degrades to '—' if neither exists. */
  const _domByChr = {};            // chr -> sorted intervals
  const _domChrProm = {};          // chr -> in-flight promise
  let _domCombined = null, _domCombinedProm = null;
  function _loadCombinedDomains(){
    if (_domCombined) return Promise.resolve(_domCombined);
    if (_domCombinedProm) return _domCombinedProm;
    _domCombinedProm = fetch(CFG.domainsUrl, {cache:'force-cache'})
      .then(r => r.ok ? r.json() : {}).catch(() => ({}))
      .then(x => { _domCombined = x || {}; return _domCombined; });
    return _domCombinedProm;
  }
  function ensureDomains(chr){
    if (chr in _domByChr) return Promise.resolve(_domByChr[chr]);
    if (_domChrProm[chr]) return _domChrProm[chr];
    _domChrProm[chr] = fetch(CFG.domainsDir + encodeURIComponent(chr) + '.json', {cache:'force-cache'})
      .then(r => { if (!r.ok) throw 0; return r.json(); })
      .then(arr => { _domByChr[chr] = arr || []; return _domByChr[chr]; })
      .catch(() => _loadCombinedDomains().then(idx => { _domByChr[chr] = (idx && idx[chr]) || []; return _domByChr[chr]; }));
    return _domChrProm[chr];
  }
  // rows: [g_start, g_end, name, pfam, type], sorted by g_start.
  // returns "Name (PFxxxxx)" for the most specific domain covering pos, else '—'.
  function domainAt(chr, pos){
    const a = _domByChr[chr];
    if (!a || !a.length) return '—';
    let lo = 0, hi = a.length;                 // first index with g_start > pos
    while (lo < hi){ const m = (lo + hi) >> 1; if (a[m][0] <= pos) lo = m + 1; else hi = m; }
    let best = null;
    for (let i = lo - 1; i >= 0; i--){
      const iv = a[i];
      if (iv[0] < pos - 100000) break;         // domain blocks are exon-sized
      if (iv[0] <= pos && iv[1] >= pos){
        if (!best || (iv[1] - iv[0]) < (best[1] - best[0])) best = iv;   // smallest = most specific
      }
    }
    return best ? `${best[2]}${best[3] ? ` (${best[3]})` : ''}` : '—';
  }

  /* Coding-consequence classification for the structural view.
       missense -> full residue-level treatment (AA change + LM scores)
       lof      -> truncation / position marker (stop gained, frameshift, start/stop lost)
       indel    -> in-frame insertion / deletion (position marker)
       null     -> not shown (synonymous, stop_retained, splice, UTR, intron, intergenic ...) */
  function classifyConsequence(effect){
    const e = String(effect || '').toLowerCase().replace(/[\s]+/g, '_');
    if (/missense|protein_altering|non[_-]?synonymous/.test(e)) return {klass:'missense', label:'Missense',            structural:true};
    if (/stop_gained|nonsense/.test(e))                         return {klass:'lof',      label:'Stop gained',         structural:false};
    if (/frameshift/.test(e))                                   return {klass:'lof',      label:'Frameshift',          structural:false};
    if (/start_lost|initiator_codon/.test(e))                   return {klass:'lof',      label:'Start lost',          structural:false};
    if (/stop_lost/.test(e))                                    return {klass:'lof',      label:'Stop lost',           structural:false};
    if (/inframe_insertion/.test(e))                            return {klass:'indel',    label:'In-frame insertion',  structural:false};
    if (/inframe_deletion/.test(e))                             return {klass:'indel',    label:'In-frame deletion',   structural:false};
    return null;
  }

  const AA3 = {ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',
    LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V',TER:'*',SEC:'U'};
  /* Parse INFO SUB into {ref, resi, alt}. Handles 1-letter (A1V, R441*),
     3-letter (Ala1Val, Trp441Ter) and frameshift forms. null if no residue. */
  function parseSub(sub){
    if (sub == null) return null;
    const s = String(sub).trim().replace(/^p\./i, '');
    let m = s.match(/^([A-Z*])(\d+)(fs\*?|[A-Z*]|=)?$/);
    if (m) return {ref:m[1], resi:+m[2], alt:(m[3]||'').replace('fs*','fs') || null};
    m = s.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3}|Ter|\*|fs)?/);
    if (m){
      const ref = AA3[m[1].toUpperCase()] || m[1][0].toUpperCase();
      let alt = null;
      if (m[3]) alt = AA3[m[3].toUpperCase()] || (m[3] === '*' ? '*' : (/fs/i.test(m[3]) ? 'fs' : m[3][0].toUpperCase()));
      return {ref, resi:+m[2], alt};
    }
    m = s.match(/(\d+)/);
    return m ? {ref:'', resi:+m[1], alt:null} : null;
  }
  function hgvsProtein(p, cls){
    if (!p) return null;
    switch (cls.label){
      case 'Frameshift':  return 'p.' + (p.ref || '') + p.resi + 'fs';
      case 'Stop gained': return 'p.' + (p.ref || '') + p.resi + '*';
      case 'Start lost':  return 'p.' + (p.ref || 'M') + p.resi + '?';
      case 'Stop lost':   return 'p.*' + p.resi + (p.alt && p.alt !== '*' ? p.alt : 'ext');
    }
    if (cls.klass === 'indel') return 'p.' + (p.ref || '') + p.resi + (cls.label.indexOf('insertion') >= 0 ? 'ins' : 'del');
    return 'p.' + (p.ref || '') + p.resi + (p.alt || '');   // missense
  }

  /* LIVE: coding variants for one gene, from the same HDF5 -> VCF pipeline as
     SNPVersity. Residue + AA change come from INFO SUB; consequence from TYPE;
     scores from PlantCAD/ESM (pc1/esm1). Nothing fabricated. Async.
     `ids` defaults to a representative selection (site-level INFO is panel-wide,
     so this stays small); pass a broader set for an exhaustive allele catalog. */
  async function queryFoldVariants(gene, dataset, ids){
    dataset = dataset || DATASETS[0].id;
    const g = await lookupGene(gene);
    if (!g || !g.chr) return [];
    ids = ids || defaultSelectionFor(dataset);
    let res;
    try { res = await queryVariants(dataset, g.chr, g.start, g.end, ids, {forceTable:true}); }
    catch (e){ console.warn('queryFoldVariants:', e && e.message); return []; }
    const out = [];
    for (const r of (res.rows || [])){
      if (gene && r.gene && r.gene !== gene && r.gene !== '—') continue;
      const cls = classifyConsequence(r.effect);
      if (!cls) continue;
      const p = parseSub(r.sub);
      if (!p || p.resi == null) continue;                 // need a residue to place it
      const pc  = (r.pc1  != null ? r.pc1  : null);
      const esm = (r.esm1 != null ? r.esm1 : null);
      const combined = (pc != null && esm != null) ? +(((pc + esm) / 2)).toFixed(2)
                     : (pc != null ? pc : (esm != null ? esm : null));
      out.push({
        id:'f' + out.length, gene, resi:p.resi, ref:p.ref, alt:p.alt,
        variant: hgvsProtein(p, cls), consequence: cls.label, consClass: cls.klass,
        structural: cls.structural, pos:r.pos, refNt:r.ref, altNt:r.alt,
        impact: r.impact || null, maf: (r.maf != null ? r.maf : null),
        plantcad: pc, esm: esm,
        plantcad2: (r.pc2 != null ? r.pc2 : null), esm2: (r.esm2 != null ? r.esm2 : null),
        esm3: (r.esm3 != null ? r.esm3 : null),
        combined,
        priority: combined == null ? null
                : (combined <= -7 ? 'TOP' : combined <= -4 ? 'HIGH' : combined <= -1 ? 'MODERATE' : 'LOW'),
      });
    }
    out.sort((a, b) => a.resi - b.resi);
    return out;
  }

  /* ---------------- SNPImpact: rank a region's variants ---------------- */
  // broader than the fold classifier: also covers splice / UTR / intron / intergenic
  function impactClass(effect){
    const e = String(effect || '').toLowerCase().replace(/[\s]+/g, '_');
    if (/missense|protein_altering|non[_-]?synonymous/.test(e)) return {klass:'missense', label:'Missense',           severe:false};
    if (/stop_gained|nonsense/.test(e))         return {klass:'lof',    label:'Stop gained',        severe:true};
    if (/frameshift/.test(e))                    return {klass:'lof',    label:'Frameshift',         severe:true};
    if (/start_lost|initiator/.test(e))          return {klass:'lof',    label:'Start lost',         severe:true};
    if (/stop_lost/.test(e))                     return {klass:'lof',    label:'Stop lost',          severe:true};
    if (/splice_(acceptor|donor)/.test(e))       return {klass:'splice', label:'Splice site',        severe:true};
    if (/splice/.test(e))                        return {klass:'splice', label:'Splice region',      severe:false};
    if (/inframe_insertion/.test(e))             return {klass:'indel',  label:'In-frame insertion', severe:false};
    if (/inframe_deletion/.test(e))              return {klass:'indel',  label:'In-frame deletion',  severe:false};
    if (/synonymous|stop_retained/.test(e))      return {klass:'syn',    label:'Synonymous',         severe:false};
    if (/intron/.test(e))                        return {klass:'other',  label:'Intron',             severe:false};
    if (/5_prime_utr|five_prime/.test(e))        return {klass:'other',  label:'5\u2032 UTR',        severe:false};
    if (/3_prime_utr|three_prime/.test(e))       return {klass:'other',  label:'3\u2032 UTR',        severe:false};
    if (/upstream/.test(e))                      return {klass:'other',  label:'Upstream',           severe:false};
    if (/downstream/.test(e))                    return {klass:'other',  label:'Downstream',         severe:false};
    if (/intergenic/.test(e))                    return {klass:'other',  label:'Intergenic',         severe:false};
    return {klass:'other', label: (cleanTok(effect) || 'Other'), severe:false};
  }
  function impactPriority(cls, combined, level){
    if (cls.severe) return 'TOP';                          // LOF, splice donor/acceptor
    if (cls.klass === 'missense' || cls.klass === 'indel'){
      if (combined != null){
        if (combined <= -7) return 'TOP';
        if (combined <= -4) return 'HIGH';
        if (combined <= -1) return 'MODERATE';
        return 'LOW';
      }
      return level === 'HIGH' ? 'HIGH' : level === 'MODERATE' ? 'MODERATE' : 'LOW';
    }
    return level === 'HIGH' ? 'HIGH' : 'LOW';               // syn / non-coding
  }
  // rank every variant in a region (from parseVcf rows). Accessions are irrelevant here.
  function rankImpact(rows){
    const out = [];
    for (const r of (rows || [])){
      const cls = impactClass(r.effect);
      const pc = r.pc1 != null ? r.pc1 : null, esm = r.esm1 != null ? r.esm1 : null;
      const combined = (pc != null && esm != null) ? +(((pc + esm) / 2)).toFixed(2)
                     : (pc != null ? pc : (esm != null ? esm : null));
      const p = parseSub(r.sub);
      const coding = (cls.klass === 'missense' || cls.klass === 'lof' || cls.klass === 'indel');
      const variant = (p && p.resi != null && coding)
        ? hgvsProtein(p, {label: cls.label, klass: cls.klass})
        : `${r.pos} ${r.ref}>${r.alt}`;
      out.push({
        id: 'i' + out.length, gene: r.gene, pos: r.pos, ref: r.ref, alt: r.alt,
        variant, consequence: cls.label, consClass: cls.klass,
        resi: p ? p.resi : null, aaRef: p ? p.ref : null, aaAlt: p ? p.alt : null,
        domain: r.domain || '\u2014', impactLevel: r.impact || null,
        plantcad: pc, esm: esm, combined,
        plantcad2: (r.pc2 != null ? r.pc2 : null),
        pc1: r.pc1, pc2: r.pc2, esm1: r.esm1, esm2: r.esm2, esm3: r.esm3,
        maf: r.maf, r2: r.r2, mq: r.mq,
        priority: impactPriority(cls, combined, r.impact), percentile: null,
      });
    }
    // real region percentile: most deleterious (most negative combined) -> highest
    const scored = out.filter(v => v.combined != null).slice().sort((a, b) => a.combined - b.combined);
    const N = scored.length;
    scored.forEach((v, i) => { v.percentile = N ? Math.round(100 * (N - i) / N) : null; });
    return out;
  }

  /* gene -> canonical protein domains (protein coords) for the SNPImpact detail track */
  let _geneDom = null, _geneDomProm = null;
  function ensureGeneDomains(){
    if (_geneDom) return Promise.resolve(_geneDom);
    if (_geneDomProm) return _geneDomProm;
    _geneDomProm = fetch(CFG.domainsGeneUrl, {cache:'force-cache'})
      .then(r => r.ok ? r.json() : {}).catch(() => ({}))
      .then(x => { _geneDom = x || {}; return _geneDom; });
    return _geneDomProm;
  }
  function geneDomains(gene){ return (_geneDom || {})[gene] || null; }

  /* per-chromosome exon/CDS structure (canonical transcript) for SNPImpact's gene-model view */
  const _gmByChr = {}, _gmProm = {};
  function ensureGeneModels(chr){
    if (chr in _gmByChr) return Promise.resolve(_gmByChr[chr]);
    if (_gmProm[chr]) return _gmProm[chr];
    _gmProm[chr] = fetch(CFG.geneModelsDir + encodeURIComponent(chr) + '.json', {cache:'force-cache'})
      .then(r => r.ok ? r.json() : {}).catch(() => ({}))
      .then(x => { _gmByChr[chr] = x || {}; return _gmByChr[chr]; });
    return _gmProm[chr];
  }
  function geneModelOf(chr, gene){ return (_gmByChr[chr] || {})[gene] || null; }

  /* ---------- SNPFunction: gene-scoped functional dossier + allele mining ---------- */
  // gts entries are now 1-byte codes from parseVcf: 0/1/2 = dosage, 3 = missing.
  function _dose(g){ return (g == null || g === 3) ? null : g; }
  function _avg(a){ return a.length ? +(a.reduce((s,x)=>s+x,0)/a.length).toFixed(2) : null; }

  // Aggregate a gene across the WHOLE panel: which accessions carry damaging/LOF alleles.
  async function geneFunction(gene, dataset){
    dataset = dataset || DATASETS[0].id;
    const dsName = (DATASETS.find(d=>d.id===dataset)||{}).name || dataset;
    const g = await lookupGene(gene);
    if (!g || !g.chr) return {gene, dataset, datasetName:dsName, error:'No gene-model coordinates found for '+gene+'.'};
    const ids = (accessionsFor(dataset)||[]).map(a=>a.id);            // FULL panel
    await Promise.all([ensureDomains(g.chr), ensureGeneDomains(), ensureGeneModels(g.chr)]);
    let res;
    try { res = await queryVariants(dataset, g.chr, g.start, g.end, ids, {forceTable:true}); }
    catch (e){ return {gene, chr:g.chr, start:g.start, end:g.end, dataset, datasetName:dsName, error:'Variant query failed: '+(e&&e.message)}; }

    const accs = res.accs || [];
    const rows = (res.rows || []).filter(r => r.gene === gene);
    const gd = geneDomains(gene), gm = geneModelOf(g.chr, gene);

    const variants = rows.map((r, vi) => {
      const cls = impactClass(r.effect);
      const pc = r.pc1 != null ? r.pc1 : null, esm = r.esm1 != null ? r.esm1 : null;
      const pc2 = r.pc2 != null ? r.pc2 : null, esm2 = r.esm2 != null ? r.esm2 : null;
      const esm3 = r.esm3 != null ? r.esm3 : null;
      const combined = (pc != null && esm != null) ? +(((pc + esm) / 2)).toFixed(2) : (pc != null ? pc : (esm != null ? esm : null));
      let het=0, hom=0, called=0; const homIds=[], hetIds=[];
      for (let k=0;k<accs.length;k++){ const d=_dose(r.gts[k]);
        if (d==null) continue; called++;
        if (d===2){ hom++; homIds.push(accs[k].id); } else if (d===1){ het++; hetIds.push(accs[k].id); } }
      const an = 2*called, ac = het + 2*hom, af = an ? ac/an : 0;
      const p = parseSub(r.sub);
      const coding = (cls.klass==='missense'||cls.klass==='lof'||cls.klass==='indel');
      const variant = (p && p.resi!=null && coding) ? hgvsProtein(p, {label:cls.label, klass:cls.klass}) : `${r.pos} ${r.ref}>${r.alt}`;
      return {id:'fx'+vi, pos:r.pos, ref:r.ref, alt:r.alt, variant, consequence:cls.label, consClass:cls.klass, severe:cls.severe,
        domain:r.domain||'\u2014', resi:p?p.resi:null, aaRef:p?p.ref:null, aaAlt:p?p.alt:null, plantcad:pc, esm, plantcad2:pc2, esm2, esm3, combined,
        priority: impactPriority(cls, combined, r.impact),
        het, hom, af, carriersHom:homIds, carriersHet:hetIds};
    });

    const byClass = {missense:0, lof:0, splice:0, indel:0, syn:0, other:0};
    variants.forEach(v => { byClass[v.consClass] = (byClass[v.consClass]||0) + 1; });
    const nonsyn = byClass.missense + byClass.lof + byClass.indel + byClass.splice, syn = byClass.syn;
    const domainDisrupting = variants.filter(v => v.domain!=='\u2014' && (v.consClass==='missense'||v.consClass==='lof'||v.consClass==='indel')).length;
    const afSpectrum = {
      rare:   variants.filter(v => v.af>0 && v.af<0.01).length,
      low:    variants.filter(v => v.af>=0.01 && v.af<0.05).length,
      common: variants.filter(v => v.af>=0.05).length,
    };
    // exon vs intron: prefer the real gene model; else fall back to consequence
    let exonic=0, intronic=0;
    for (const v of variants){
      let inExon;
      if (gm && gm.exons) inExon = gm.exons.some(e => v.pos>=e[0] && v.pos<=e[1]);
      else inExon = (v.consClass!=='other') || /utr/i.test(v.consequence);
      if (/intron/i.test(v.consequence) || (gm && gm.exons && !inExon)) intronic++;
      else if (inExon) exonic++;
    }
    const PR = ['TOP','HIGH','MODERATE','LOW'];
    const damaging = variants
      .filter(v => v.consClass==='lof' || v.severe || (v.consClass==='missense' && v.combined!=null && v.combined<=-4))
      .sort((a,b) => (PR.indexOf(a.priority)-PR.indexOf(b.priority)) || ((a.combined==null?0:a.combined)-(b.combined==null?0:b.combined)));
    const koGenotypes = damaging.filter(v=>v.consClass==='lof').reduce((n,v)=>n+v.hom, 0);
    const koLines = new Set(); damaging.filter(v=>v.consClass==='lof').forEach(v=>v.carriersHom.forEach(id=>koLines.add(id)));

    return {
      gene, chr:g.chr, start:g.start, end:g.end, strand: gm?gm.strand:null, dataset, datasetName:dsName,
      nAccessions: accs.length, nVariants: variants.length,
      protLen: gd?gd.len:null, protein: gd?gd.protein:null, domains: gd?(gd.domains||[]):[],
      burden: { byClass, nonsyn, syn, nonsynSyn: syn ? +(nonsyn/syn).toFixed(2) : (nonsyn?null:0),
                exonic, intronic, exonIntron: intronic ? +(exonic/intronic).toFixed(2) : (exonic?null:0),
                domainDisrupting, meanPlantcad:_avg(variants.map(v=>v.plantcad).filter(x=>x!=null)),
                meanEsm:_avg(variants.map(v=>v.esm).filter(x=>x!=null)),
                meanPlantcad2:_avg(variants.map(v=>v.plantcad2).filter(x=>x!=null)),
                meanEsm2:_avg(variants.map(v=>v.esm2).filter(x=>x!=null)),
                meanEsm3:_avg(variants.map(v=>v.esm3).filter(x=>x!=null)), afSpectrum },
      damaging, koGenotypes, koLines: koLines.size, variants,
    };
  }


  const DEFAULT_DS = DATASETS[0].id;
  return {
    datasets:    () => DATASETS,
    // dataset-aware accession accessors
    projectsFor, accessionsFor, defaultSelectionFor, familyOf, hasSecondaryScores, namFoundersFor,
    // backwards-compatible defaults (first dataset)
    projects:    () => projectsFor(DEFAULT_DS),
    accessions:  () => accessionsFor(DEFAULT_DS),
    defaultSelection: () => defaultSelectionFor(DEFAULT_DS),
    geneModels:  () => GENE_MODELS,
    exampleGenes:() => EXAMPLE_GENES,
    lookupGene,             // async: gene model id -> {id, chr, start, end} | null
    chromLengths:() => CHR_LEN,
    centromeres: () => CENTRO,
    estimateResult,                        // run-bar prediction: table vs VCF download
    accessionById,
    queryVariants,          // now async (returns a Promise)
    queryImpact,
    structureFor,
    pdbFor,
    ensureStructure,
    queryFoldVariants,
    rankImpact,
    ensureGeneDomains,
    geneDomains,
    ensureGeneModels,
    geneModelOf,
    geneFunction,
  };
})();

/* Global helper: render a "Name (PFxxxxx)" domain string as a chip with the
   Pfam accession linked to InterPro. Used by SNPVersity / SNPImpact / SNPFunction. */
function pfamHref(pf){ return 'https://www.ebi.ac.uk/interpro/entry/pfam/' + pf + '/'; }
function domTag(dom){
  if (dom == null || dom === '\u2014' || dom === 'N/A' || dom === '')
    return '<span style="color:var(--faint)">\u2014</span>';
  var esc = function(s){ return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;'); };
  var m = String(dom).match(/(PF\d{4,6})/);
  if (m){
    var pf = m[1];
    var name = String(dom).replace(/\s*\(?\s*PF\d{4,6}\s*\)?\s*/, ' ').trim();
    return '<span class="dom-tag">' + esc(name) +
      ' <a href="' + pfamHref(pf) + '" target="_blank" rel="noopener" title="Pfam ' + pf +
      '" style="color:inherit;text-decoration:underline;text-decoration-style:dotted">' + pf + '</a></span>';
  }
  return '<span class="dom-tag">' + esc(dom) + '</span>';
}
