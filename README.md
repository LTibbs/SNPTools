# SNPTools — SNPVersity 2.1

> A browser-based suite of tools for exploring maize genomic variation, developed by
> [MaizeGDB](https://maizegdb.org) and funded by USDA-ARS.

SNPTools is a single-page web application that lets researchers query variant datasets
across thousands of maize accessions, annotate variants with predicted effects and
AI-based language-model scores, and hand a selection off to a family of connected
analysis tools (phylogeny, distance matrix, similarity, gene function, protein
structure, and more). All coordinates are on the **B73 v5** reference assembly.

<!-- TODO: add a screenshot or short GIF of SNPVersity here, e.g. docs/screenshot.png -->
<!-- TODO: add a live demo URL, if one is hosted -->

---

## Table of contents

- [How it works](#how-it-works)
- [Repository layout](#repository-layout)
- [The tools](#the-tools)
- [Datasets](#datasets)
- [Backend services](#backend-services)
- [Data files consumed by the frontend](#data-files-consumed-by-the-frontend)
- [Data-generation scripts](#data-generation-scripts)
- [Getting started](#getting-started)
- [Configuration](#configuration)
- [Data provenance & citation](#data-provenance--citation)
- [Acknowledgments](#acknowledgments)
- [License](#license)

---

## How it works

The frontend is plain, dependency-light JavaScript (no build step). It is organized as a
tiny shell plus a set of self-registering tool modules:

- **`core.js`** loads first and provides the shared "shell": the icon set, the tool
  registry, cross-tool shared state (`S`), the hash-free router (`go(id)`), the left-rail
  navigation, tooltips, and menus. Each tool module registers itself with
  `SNPTools.register(id, { render })`, so `core.js` never needs to know their internals.
- **`data.js`** is the single data layer. Everything that touches the variant store,
  accession catalog, gene models, domains, or protein structures goes through the `Data`
  object. Tools never fetch data directly.
- **Tool modules** (`snpversity.js`, `snpimpact.js`, …) each own one page and its
  view-state. A selection made in one tool is stashed on the shared `S` object (e.g.
   `S.treeInput`, `S.impactInput`) and picked up by the next, which is how
  "Send selection to…" works.

At runtime, a SNPVersity query POSTs a region + accession list to a small PHP backend,
which runs a Python extractor against the real HDF5 variant store, writes a VCF, and
returns its path. The frontend fetches and parses that VCF into the row shape every tool
expects, then layers on domain and gene-model annotations loaded on demand per chromosome.

```
Browser (core.js + tool modules)
        │  region + accessions
        ▼
processForm.php ──► h5_to_vcf.py ──► .h5 variant store ──► writes VCF
        │  returns VCF path
        ▼
data.js parses VCF ──► rows ──► SNPVersity table + hand-off to other tools
        │
        ├─ data/domains/by_chr/<chr>.json     (Pfam domain by position)
        ├─ data/genemodels/by_chr/<chr>.json  (exon/CDS structure)
        └─ data/structures/structure-<gene>.js (AlphaFold model, for SNPFold)
```

---

## Repository layout

The exact on-disk layout may differ from your deployment; this reflects what the code
references.

```
.
├── index.html                      # app shell + MaizeGDB masthead
├── css/
│   └── main.css                    # <!-- TODO: not included in this drop; add it -->
├── images/
│   └── maizegdb_logo.png           # <!-- TODO: add brand asset -->
├── js/
│   ├── core.js                     # shell: registry, router, nav, state
│   ├── data.js                     # data layer (HDF5 → VCF → rows, annotations)
│   ├── accessions.catalog.js       # window.SNP_CATALOG (projects → groups → accessions)
│   ├── accessions.real.js          # window.SNP_REAL_ACCESSIONS (legacy flat fallback)
│   ├── snpversity.js               # Visualization & Search (front door)
│   ├── snpimpact.js                # AI-based variant prioritization
│   ├── snpfold.js                  # protein-structure-aware interpretation
│   ├── snptree.js                  # local phylogeny (IBS / NJ / UPGMA)
│   ├── snpcompare.js               # similarity to a focal accession
│   ├── snpfunction.js              # gene function & allele mining
│   └── snphelp.js                  # Help & FAQ page
├── data/
│   ├── structures/
│       ├── alphafold/
│           └── structure-<gene>.js     # per-gene AlphaFold2 model + metadata (SNPFold)
│       ├── boltz/
│           └── structure-<gene>.js     # per-gene Boltz2 model + metadata (SNPFold)
|       └── esmfold/                 
│           └── structure-<gene>.js     # per-gene ESMFold model + metadata (SNPFold)
│   ├── domains/
│   │   ├── by_chr/<chr>.json        # Pfam domains by genomic position
│   │   ├── domains.by_protein.json  # domains in protein coords (SNPFold track)
│   │   ├── domains.by_gene.json     # canonical protein's domains per gene
│   │   └── domain_lookup.js         # binary-search helper (generated)
│   └── genemodels/
│       └── by_chr/<chr>.json        # canonical exon/CDS structure (SNPImpact)
├── vcf/                             # server-written VCF outputs (web-served, writable)
│
├── scripts/                         # data-generation pipeline (see below)
│   ├── map_domains_to_dna.py
│   ├── compile_domains.py
│   ├── compile_exons.py
│   ├── pdb_cif_to_snpfold_dssp_parallel.py
│   └── alphafold2snpfold_dssp_parallel.sh
│
└── backend/                         # <!-- TODO: PHP endpoints not included in this drop -->
    ├── processForm.php              # region+accessions → VCF (calls h5_to_vcf.py)
    ├── lookupGeneModel.php          # gene model ID → coordinates
    ├── ibsCompare.php               # genome-wide IBS similarity for a focal accession
    ├── h5_to_vcf.py                 # HDF5 extractor invoked by processForm.php
    └── build_catalog.py             # builds accessions.catalog.js from TSVs
```

> **Note:** `css/main.css`, `images/`, the backend PHP/Python endpoints
> (`processForm.php`, `lookupGeneModel.php`, `ibsCompare.php`), `h5_to_vcf.py`, and
> `build_catalog.py` are referenced by the app but were not part of this code drop. Add
> them (or update the TODOs above) before publishing.

---

## The tools

The suite ships several **live** tools plus additional **in-development** modules that
render informational placeholder pages. All live tools share the same variant database,
annotations, and B73 v5 coordinates, and selections flow freely between them.

| Tool | Status | What it does |
| --- | --- | --- |
| **SNPVersity** | Live | The front door. Choose a dataset, enter a genomic interval or a B73 v5 gene-model ID, and pick accessions. Returns a color-coded genotype table plus a downloadable VCF — allele states, predicted effects, Pfam domain per position, and DNA/protein language-model scores included. Use **Send selection to…** to hand the result to any other tool. |
| **SNPImpact** | Live | Ranks the variants in a region (independent of which accessions were chosen) by a combined AI score (PlantCAD DNA-model + ESM protein-model) together with predicted consequence and Pfam domain. Filter by consequence, priority, score, or domain, and build a shortlist of candidate causal alleles. |
| **SNPFunction** | Live | A gene-scoped dossier that is independent of any one region. Enter a gene and it summarizes the Pfam domains and protein, computes the gene's variant burden across the **whole panel**, and lists which accessions carry each damaging or knockout allele. |
| **SNPFold** | Live | Structure-aware interpretation of coding variants. Aligns variants, Pfam domains, secondary structure, and per-residue pLDDT confidence on a linear "protein browser" track; renders an interactive 3D model (3Dmol.js, loaded on demand); and gives a per-variant structural-context readout. Structure data comes from precomputed `structure-<gene>.js` files. |
| **SNPCompare** | Live | Ranks every accession by identity-by-state similarity to a chosen focal accession, in three scopes: genome-wide (precomputed, served per focal accession), this region (computed in-browser from a SNPVersity result), or both side-by-side with Δ = local − global to surface region-specific relatedness (introgression / selection). |
| **SNPTree** | Live | Builds a local phylogenetic / similarity tree from the in-memory genotype matrix using identity-by-state distances (UPGMA / Neighbor-Joining). Exports Newick / MEGA / PHYLIP. Same method vocabulary as VCF2PopTree, reimplemented to run on the in-memory matrix and MaizeGDB metadata. |
| **SNPTrait** | In development | Will connect variation to phenotype and trait records from the National Germplasm collection; filter 20,000+ accessions by trait/metadata, then hand a set to the genomic tools. |
| **SNPImpute** | In development | Pan-genome–guided imputation to fill missing genotypes and predict function across sequencing depths. |
| **SNPGermplasm** | In development | Genotype-driven collection management: redundancy, uniqueness, and priority materials. |

<!-- TODO: confirm the intended public status of SNPImpute / SNPDensity / SNPGermplasm /
     SNPTrait before release; adjust the table if any have shipped. -->

---

## Datasets

Datasets are defined in `data.js` (`DATASETS`). Each maps to a `.h5` family and its real
accession columns. Only the **MaizeGDB 2026** family carries the second-generation
language-model scores (PlantCAD2 / ESM2 / ESM3); older families expose a single
`DNA_SCORE` (PlantCaduceus) and `AA_SCORE` (ESM1b).

| Dataset ID | Name | Subset | Reference | Accessions | Sites | Filters | Het | INDELs | Imputed |
| --- | --- | --- | --- | --- | --- | --- | :-: | :-: | :-: |
| `mgdb2026_hq` | MaizeGDB 2026 | High Quality | B73 v5 | 2,710 | 98M | MQ ≥ 30, Coverage ≥ 50%, LD max R² > 0.5 | ✓ | ✓ | – |
| `mgdb2026_hc` | MaizeGDB 2026 | High Coverage | B73 v5 | 2,710 | 290M | MQ ≥ 30, Coverage ≥ 50% | ✓ | ✓ | – |
| `schnable2023` | Schnable 2023 | Imputed markers | B73 v5 | 1,515 | 12M | Imputed | – | – | ✓ |
| `nam2021` | NAM 2021 | Founder panel | B73 v5 | 27 | 78M | MQ ≥ 30, Founder panel | ✓ | ✓ | – |
| `mgdb2024_hq` | MaizeGDB 2024 | High Quality | B73 v5 | 1,498 | 83M | MQ ≥ 30, Coverage ≥ 50% | ✓ | ✓ | – |

### Accession catalog

Accessions are organized as **projects (bioprojects) → groups → accessions** in
`window.SNP_CATALOG` (`accessions.catalog.js`, auto-generated from `accessions.tsv` +
`projects.tsv`). A legacy flat list, `window.SNP_REAL_ACCESSIONS`
(`accessions.real.js`), is used as a fallback. Accession IDs are the actual column names
inside the `.h5` files, so a selection maps directly to HDF5 columns.

The underlying resequencing data spans many public bioprojects (NCBI SRA / CNCB),
including PRJNA783885, PRJCA009749, PRJNA531553, PRJNA389800, PRJNA399729, PRJNA479960,
PRJEB14212, PRJEB56320, and others; friendly names are kept in `popNames`.

<!-- TODO: confirm accession/site counts against the current .h5 build before release;
     the numbers above are the UI-card values from data.js. -->

---

## Backend services

The frontend talks to three server endpoints (paths configurable in `data.js` → `CFG`).
These were **not** part of this code drop — add them or point the config at your service.

| Endpoint | Called by | Purpose |
| --- | --- | --- |
| `processForm.php` | `Data.queryVariants` | Accepts `chr`, `start`, `end`, `dataSet`, `genotypes` (JSON list), `outName`. Runs `h5_to_vcf.py` against the `.h5` store, writes a VCF into `vcf/`, and returns JSON `{ status, outFile, … }`. |
| `lookupGeneModel.php` | `Data.lookupGene` | Resolves a B73 v5 gene-model ID (`?geneModelId=…`) to `{ chromosome, start, end }` from the serialized GFF store. |
| `ibsCompare.php` | SNPCompare (global mode) | Returns precomputed genome-wide IBS similarity rows for a focal accession (`?focal=<ID>` → `{ rows: [{ id, similarity, missing }] }`). |

VCF `INFO` tags the parser understands include: `GENEMODEL`, `TYPE`, `SUB`, `MQ`, `CVP`,
`MAXR2`, `MAF`, and score fields `plantcad1_score` / `plantcad2_score` /
`ESM1_score` / `ESM2_score` / `ESM3_score` (2026), or `DNA_SCORE` / `AA_SCORE` (older
families).

<!-- TODO: document the exact h5_to_vcf.py CLI and the .h5 store layout / filename
     convention (the code expects a "chr10"-style chromosome token in the filename). -->

---

## Data files consumed by the frontend

These static files are produced by the [data-generation scripts](#data-generation-scripts)
and loaded on demand (cached) by `data.js`:

- **`data/domains/by_chr/<chr>.json`** — Pfam domains keyed by genomic position, one entry
  per exon block, sorted by start; powers the SNPVersity "Domain" column and SNPImpact.
  Falls back to a combined `data/domains/domains.by_chr.json` if the per-chromosome file is
  absent.
- **`data/domains/domains.by_protein.json`** — domains in protein coordinates; SNPFold's
  linear track.
- **`data/domains/domains.by_gene.json`** — each gene's canonical protein domains; used as
  the SNPImpact detail track and as the canonical-transcript map for the pipeline.
- **`data/genemodels/by_chr/<chr>.json`** — canonical-transcript exon/CDS structure per
  gene; SNPImpact's gene-model view and SNPFunction's exon/intron split.
- **`data/structures/structure-<gene>.js`** — per-gene AlphaFold model plus lightweight
  metadata (`seq`, `plddt[]`, `ss[]`, `domains[]`, `length`, `uniprot`); loaded lazily by
  SNPFold. Defines `window.SNPFOLD_PDB[gene]` and `window.SNPFOLD_STRUCT[gene]`.

---

## Data-generation scripts

The `scripts/` directory builds the static annotation and structure files above. A typical
pipeline order:

### 1. `map_domains_to_dna.py`

Projects Pfam domains from **protein** coordinates onto **genomic DNA** coordinates using
the CDS structure and strand from a GFF3 (walks codons back through the CDS,
reverse-complement order on the minus strand, skipping introns). Emits one BED-like row
per contiguous genomic interval.

```bash
python3 map_domains_to_dna.py \
  --gff B73v5.gff3 \
  --domains protein_pfam_domains.tsv \
  --out domain_cds_dna_coords.tsv
```

### 2. `compile_domains.py`

Turns the two Pfam tables (protein-coordinate hits + the genomic blocks from step 1) into
the compact, view-ready lookups the frontend loads: `domains.by_chr/` (per chromosome),
`domains.by_protein.json`, `domains.by_gene.json`, and a `domain_lookup.js` binary-search
helper. Pfam versions are stripped (e.g. `PF02536.20` → `PF02536`, full id kept as
`pfam_ver`). Canonical transcripts only by default.

```bash
python3 compile_domains.py \
  --protein protein_pfam_domains.tsv \
  --dna     domain_cds_dna_coords.tsv \
  --outdir  data/domains/ \
  --minify
```

### 3. `compile_exons.py`

Emits the canonical-transcript exon/CDS structure for every gene, split per chromosome,
for SNPImpact's gene-model view. Canonical transcript is taken from a `--canonical` map
(e.g. `domains.by_gene.json`) with a fallback to the lowest `_T###` isoform (or the one
with the most coding sequence).

```bash
python3 compile_exons.py \
  --gff B73v5.gff3 \
  --outdir data/genemodels/ \
  --canonical data/domains/domains.by_gene.json
```

### 4. `pdb_cif_to_snpfold_dssp_parallel.py`

Converts an AlphaFold model (PDB or mmCIF) into a `structure-<gene>.js` file for SNPFold.
Reads per-residue pLDDT from the B-factor column, runs DSSP for secondary structure when
available (otherwise all-coil), and attaches Pfam domains in protein coordinates. CIF
parsing tries gemmi, then Biopython, then a built-in minimal parser. Supports sharding and
multiprocessing for whole-proteome runs.

```bash
python3 pdb_cif_to_snpfold_dssp_parallel.py \
  --indir  ./alphafold_pdbs/B73/ \
  --outdir ./js/structures_dssp_parallel/ \
  --domains-json ./data/domains/domains.by_protein.json \
  --canonical    ./data/domains/domains.by_gene.json \
  --workers 8 --skip-existing
```

### 5. `alphafold2snpfold_dssp_parallel.sh`

A SLURM batch wrapper that runs step 4 across a job array (parallelism = array size ×
`cpus-per-task`, e.g. 50 tasks × 8 cores = 400-way). Each task writes unique filenames into
a shared `--outdir`, and `--skip-existing` makes re-submission safe for mopping up
timed-out shards.

```bash
sbatch alphafold2snpfold_dssp_parallel.sh
```

### Also referenced (not included in this drop)

- **`h5_to_vcf.py`** — extracts a region + accession list from the `.h5` store to a VCF
  (invoked by `processForm.php`).
- **`build_catalog.py`** — builds `accessions.catalog.js` (`window.SNP_CATALOG`) from
  `data/accessions.tsv` + `data/projects.tsv`.

<!-- TODO: add h5_to_vcf.py and build_catalog.py to scripts/, or document where they live. -->

---

## Getting started

### Prerequisites

- A web server that can serve static files **and** execute the PHP endpoints
  (e.g. Apache or nginx + PHP-FPM). <!-- TODO: confirm required PHP version -->
- Python 3 on the server for `h5_to_vcf.py` and the data-generation scripts.
- Python packages for the pipeline: standard library for the domain/exon scripts; the
  structure script optionally uses `gemmi` and/or `biopython`, and DSSP (`mkdssp`) for
  secondary structure. <!-- TODO: pin exact versions / provide requirements.txt -->
- The HDF5 variant store (`.h5` files) referenced by `h5_to_vcf.py`.
  <!-- TODO: document where to obtain / stage these. -->

### Local development (frontend only)

Because the frontend is static, you can preview the UI with any static server, though
variant queries will fail without the PHP backend and `.h5` store:

```bash
# from the repository root
python3 -m http.server 8000
# then open http://localhost:8000/
```

### Full deployment

1. Place the frontend (`index.html`, `css/`, `js/`, `images/`, `data/`) under your web root.
2. Deploy the backend endpoints and point `data.js` → `CFG` at them.
3. Stage the `.h5` variant store and make `vcf/` writable by the web server.
4. Generate the static annotation/structure files with the pipeline scripts.

<!-- TODO: add concrete server config (vhost example, PHP settings, permissions). -->

---

## Configuration

Backend paths and behavior live in `data.js` under `CFG`:

| Key | Default | Purpose |
| --- | --- | --- |
| `endpoint` | `processForm.php` | Variant query endpoint. |
| `vcfDir` | `vcf/` | Web-served, writable output dir (must match `processForm.php`). |
| `geneEndpoint` | `lookupGeneModel.php` | Gene-model → coordinates lookup. |
| `structDir` | `data/structures/` | Per-gene SNPFold structure files. |
| `domainsDir` | `data/domains/by_chr/` | Per-chromosome Pfam files (preferred). |
| `domainsUrl` | `data/domains/domains.by_chr.json` | Combined domains fallback. |
| `domainsGeneUrl` | `data/domains/domains.by_gene.json` | Gene → canonical domains. |
| `geneModelsDir` | `data/genemodels/by_chr/` | Per-chromosome exon/CDS structure. |
| `tableMaxSpan` | `1_000_000` | Intervals wider than this offer a VCF download instead of an in-browser table. |

SNPCompare has its own `CFG` (`ibsCompare.php` endpoint, default dataset, and a
`useDemoGlobal` toggle for previewing without a backend).

---

## Data provenance & citation

Variation data is aggregated from public maize resequencing bioprojects (NCBI SRA and CNCB;
see `popNames` / `accessions.catalog.js` for the full list and reference DOIs). Domains are
from Pfam / InterPro; protein models are AlphaFold; secondary structure is from DSSP.
Language-model scores are from PlantCAD / PlantCaduceus (DNA) and ESM (protein).

<!-- TODO: add the preferred citation for SNPTools / SNPVersity 2.1 (paper, DOI, or
     "How to cite" text), and any dataset-specific citation requirements. -->

---

## Acknowledgments

Developed by [MaizeGDB](https://maizegdb.org). Funded by USDA-ARS.

<!-- TODO: add contributor list, grant/award numbers, and third-party library credits
     (3Dmol.js, DSSP, gemmi, Biopython, AlphaFold, Pfam/InterPro, etc.). -->

---

## License

<!-- TODO: add a license. If this is a U.S. government work, state the applicable public
     terms; otherwise add a LICENSE file and reference it here. -->
