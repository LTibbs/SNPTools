/* =====================================================================
 *  snphelp.js — Help & FAQ page for the SNPTools suite.
 *
 *  Registers itself as the 'help' tool so the router can reach it, and
 *  exposes openHelp() for the MaizeGDB masthead link (which doesn't want
 *  to depend on TOOLS[] being populated). Explains every page in the
 *  suite and every dataset defined in data.js, plus a scores/annotations
 *  glossary and a short FAQ.
 *
 *  Depends on core.js (SNPTools, S, go, ICONS, renderNav) and — for the
 *  live dataset table — data.js (Data.datasets). Degrades gracefully if
 *  Data isn't present.
 * ===================================================================== */
const SNPHelp = (function () {

  /* ---- tool descriptions, in user language, kept in step with core.js ---- */
  const LIVE = 'live', SOON = 'in development';
  const PAGES = [
    { id:'snpversity', name:'SNPVersity', icon:'dna', color:'#2563eb', status:LIVE,
      tag:'Build a variant view across accessions',
      what:'The front door of the suite and the starting point for most work. Choose a dataset, type a genomic interval (or a B73 v5 gene model ID), and pick the accessions you want. SNPVersity queries the variant store and returns a color-coded genotype table plus a downloadable VCF — allele states, predicted effects, and DNA/protein language-model scores included.',
      give:'A dataset, a region or gene, and a set of accessions.',
      get:'A genotype table and a VCF. From here, "Send selection to…" hands the same result to any other tool.' },
    { id:'snptrait', name:'SNPTrait', icon:'leaf', color:'#1f8a4c', status:SOON,
      tag:'Connect variation to traits',
      what:'Links genomic variation to phenotype and trait records from the National Germplasm collection. Search, sort, and filter accessions by trait values and metadata, then move a selected set straight into SNPVersity and the rest of the suite.',
      give:'Trait, phenotype, or metadata filters across 20,000+ accessions.',
      get:'A shortlisted set of lines you can hand off to the genomic tools.' },
    { id:'snpimpact', name:'SNPImpact', icon:'star', color:'#7c3aed', status:LIVE,
      tag:'Rank candidate variants',
      what:'Prioritizes the variants in a region regardless of which accessions you picked. It orders them by an AI-based score (PlantCAD DNA-model and ESM protein-model predictions) combined with predicted consequence and Pfam domain annotation, so likely causal changes rise to the top. Filter by consequence, priority, score, or domain, and flag a shortlist.',
      give:'A region sent from SNPVersity (accessions are ignored here).',
      get:'A ranked, filterable variant table and a shortlist of candidate alleles.' },
    { id:'snpfunction', name:'SNPFunction', icon:'func', color:'#2563eb', status:LIVE,
      tag:'Gene function & allele mining',
      what:'A gene-scoped dossier, independent of any one region. It opens with a functional-annotation record — gene symbol and names, a sourced description, protein domain architecture drawn to scale, Gene Ontology grouped by aspect, KEGG pathways and orthology, and external cross-references, each tagged by provenance (curated vs predicted). Below that it computes the gene\u2019s variant burden across the whole panel and lists which accessions carry each damaging or knockout allele.',
      give:'A gene model ID and a dataset.',
      get:'A functional-annotation dossier (identity, domains, GO, pathways, cross-references), a variant-burden breakdown, and a damaging-allele catalog with carrier lines.' },
    { id:'snpcompare', name:'SNPCompare', icon:'compare', color:'#0e7490', status:LIVE,
      tag:'Similarity to a focal accession',
      what:'Ranks every accession by identity-by-state similarity to one focal accession. Similarity is the fraction of co-called sites with the same genotype. Works genome-wide (precomputed), for the current region (computed in-browser from a SNPVersity result), or both side-by-side with a delta that surfaces region-specific relatedness such as introgression.',
      give:'A focal accession, and optionally a region handed off from SNPVersity.',
      get:'A ranked similarity table with project, SRA ID, and accession name.' },
    { id:'snptree', name:'SNPTree', icon:'tree', color:'#15803d', status:LIVE,
      tag:'Local phylogeny',
      what:'Builds a local phylogenetic tree from the genotype matrix already in memory, using identity-by-state distances (UPGMA / neighbour-joining). Useful for reading haplotype structure, introgression, and how lines cluster in a region.',
      give:'A SNPVersity result set (up to 250 accessions before a warning).',
      get:'An interactive tree, downloadable as Newick, MEGA, or PHYLIP.' },
    { id:'snpmatrix', name:'SNPMatrix', icon:'grid', color:'#b45309', status:LIVE,
      tag:'Pairwise distance matrix',
      what:'Computes the pairwise identity-by-state distance among your selected accessions and draws it as a heatmap. Reorder by input order or by clustering, switch between IBS distance and % identity, and color rows by bioproject.',
      give:'A SNPVersity result set (up to 400 accessions before a warning).',
      get:'A heatmap plus downloads: CSV distance matrix, PHYLIP, PNG, and SVG.' },
    { id:'snpimpute', name:'SNPImpute', icon:'impute', color:'#0891b2', status:SOON,
      tag:'Impute sequence & function',
      what:'Pan-genome\u2013guided imputation that fills missing genotypes and carries functional predictions through, across light and deep sequencing depths.',
      give:'A genotype set with missing calls.',
      get:'Imputed genotypes with function predictions.' },
    { id:'snpfold', name:'SNPFold', icon:'fold', color:'#be185d', status:LIVE,
      tag:'Variants on protein structure',
      what:'Maps coding variants onto predicted protein structure. A linear protein browser aligns variants with Pfam domains, secondary structure, and per-residue pLDDT confidence; an on-demand 3D viewer shows the fold colored by confidence, domain, or impact; and a per-variant readout interprets each change (domain, local confidence, secondary structure, predicted \u0394\u0394G).',
      give:'A gene with an available structure model.',
      get:'A structure-aware, per-variant interpretation of coding changes.' },
    { id:'paneffect', name:'PanEffect', icon:'effect', color:'#b45309', status:LIVE,
      tag:'Missense effects across the pan-genome',
      what:'Visualizes the predicted effect of every possible amino-acid substitution across a protein using ESM protein-language-model scores. A B73 reference view and a pan-genome view each show a full-length substitution heatmap with a zoomable window, aligned to Pfam domains and predicted secondary structure; the pan-genome view adds the natural variation seen across the maize assemblies, colored by heterotic group. Choose a gene model and an ESM model, or arrive from SNPVersity or SNPFold on a specific missense call and that substitution is highlighted in the MaizeGDB 2026 view.',
      give:'A gene model and an ESM model (ESM1 / ESM2 / ESM3); optionally a missense variant handed off from another tool.',
      get:'B73 and pan-genome substitution heatmaps with domain and secondary-structure context, plus a downloadable per-variant effects file.' },
    { id:'snpdensity', name:'SNPDensity', icon:'density', color:'#9333ea', status:SOON,
      tag:'Density & burden',
      what:'Measures SNP and INDEL density, burden, and distribution across genes and regions to highlight mutational load, constraint, and diversification.',
      give:'A region or gene set.',
      get:'Density tracks and burden summaries.' },
    { id:'snpgermplasm', name:'SNPGermplasm', icon:'germ', color:'#16a34a', status:SOON,
      tag:'Collection management',
      what:'Applies genotype-driven analytics to germplasm management, identifying redundancy, uniqueness, and priority materials for curation and deployment.',
      give:'A collection of accessions.',
      get:'Redundancy, uniqueness, and priority flags.' },
  ];

  /* ---- glossary ---- */
  const GLOSSARY = [
    ['B73 RefGen v5', 'Every coordinate, gene model, and annotation in the suite is anchored to the B73 version 5 maize reference genome.'],
    ['Genotype dosage', 'Each call is read as 0 (0/0, reference), 1 (heterozygous), or 2 (1/1, alternate homozygous). Missing calls (./.) are left out of a comparison rather than counted as a match.'],
    ['IBS distance / % identity', 'Identity-by-state compares two accessions site by site over the calls they share. Distance is the mean allele difference; % identity is 100 − distance. Used by SNPTree, SNPMatrix, and SNPCompare.'],
    ['Predicted effect & impact', 'Each variant carries a predicted consequence (missense, LOF, splice, indel, synonymous, …) rolled up into an impact tier: HIGH, MODERATE, LOW, or MODIFIER, most severe wins when a site lists several.'],
    ['PlantCAD score', 'A DNA language-model prediction of how disruptive a change is. MaizeGDB 2026 carries a second-generation PlantCAD2; older datasets carry a single DNA score.'],
    ['ESM score', 'A protein language-model prediction of the effect of an amino-acid change. MaizeGDB 2026 carries ESM2 and ESM3 alongside the original; older datasets carry a single AA score.'],
    ['MAF', 'Minor allele frequency, the frequency of the less common allele. SNPVersity can filter a region by a minimum MAF.'],
    ['Pfam domain', 'When a variant falls inside a known protein domain, that domain is shown and linked to InterPro. Domain annotation is still being loaded for some regions, where it reads as \u2014.'],
    ['VCF', 'The Variant Call Format file SNPVersity generates for your query. It is the exact matrix the other tools reuse when you send a selection.'],
    ['Gene Ontology (GO)', 'Standardized terms describing a gene product\u2019s biological process, molecular function, and cellular component. SNPFunction groups them by aspect and tags each by source: MaizeGDB and UniProt are curated, InterPro2GO is predicted from domains.'],
    ['Curated vs predicted', 'Annotation provenance. Curated evidence (MaizeGDB, UniProt) is human-reviewed; predicted evidence (InterPro / InterPro2GO) is inferred from protein domains. SNPFunction labels each description and GO term accordingly and can hide predicted-only terms.'],
    ['Heterotic group', 'A maize breeding classification (stiff-stalk, non-stiff-stalk, Iodent, Lancaster, tropical, teosinte, and others). PanEffect colors the pan-genome rows by the assembly\u2019s heterotic group.'],
  ];


  /* ---- complete definitions by tool ----
   * Keep table-header wording identical to the live tools so users can
   * search this page for the label they see on screen.
   */
  const DEFINITIONS = [
    { tool:'Shared terms', color:'#64748b', items:[
      ['Accession', 'A named maize line, sample, or sequencing run represented by one genotype column.'],
      ['Allele', 'One observed DNA state at a genomic position. REF is the B73 v5 reference allele; ALT is an alternate allele.'],
      ['Variant', 'A genomic position where at least one accession differs from the B73 v5 reference.'],
      ['SNP', 'Single-nucleotide polymorphism: a one-base substitution.'],
      ['INDEL', 'Insertion or deletion relative to the reference sequence.'],
      ['Gene model', 'The B73 v5 identifier and annotated exon, intron, CDS, and transcript structure assigned to a gene.'],
      ['Consequence / Effect', 'The predicted molecular result of a variant, such as synonymous, missense, splice-site, frameshift, stop gained, or intronic.'],
      ['Impact', 'A broad severity class assigned from the consequence: HIGH, MODERATE, LOW, or MODIFIER.'],
      ['Priority', 'The SNPTools ranking tier (TOP, HIGH, MODERATE, or LOW) produced by combining consequence, model scores, and domain context.'],
      ['Domain', 'A Pfam-annotated protein domain overlapping the affected residue; — means no loaded domain hit.'],
      ['Het', 'Heterozygous: the accession carries one reference and one alternate allele, usually displayed as 0/1 or 1/0.'],
      ['Hom', 'Alternate homozygous: the accession carries two alternate alleles, displayed as 1/1.'],
      ['Missing / ./.', 'No usable genotype call. Missing calls are excluded from pairwise similarity calculations.'],
      ['AF', 'Alternate-allele frequency: the frequency of the ALT allele among called chromosomes in the analyzed panel.'],
      ['MAF', 'Minor-allele frequency: the frequency of the less common allele, constrained to 0–0.5.'],
      ['Carrier', 'An accession with at least one copy of the alternate allele.'],
      ['Co-called sites', 'Sites where both accessions in a pair have non-missing genotype calls.'],
      ['B73 RefGen v5', 'The maize reference assembly used for coordinates, REF alleles, gene models, and annotations throughout SNPTools.'],
    ]},
    { tool:'SNPVersity', color:'#2563eb', items:[
      ['CHR', 'Reference chromosome containing the variant.'],
      ['POS', 'One-based genomic coordinate on B73 v5.'],
      ['REF', 'Reference allele in B73 v5.'],
      ['ALT', 'Alternate allele represented by the row.'],
      ['Gene model', 'B73 v5 gene model overlapping or associated with the variant.'],
      ['Effect', 'Predicted variant consequence from the annotation source.'],
      ['Impact', 'Predicted severity category: HIGH, MODERATE, LOW, or MODIFIER.'],
      ['Domain', 'Pfam protein domain overlapping the affected coding residue, when available.'],
      ['MQ', 'Mapping quality: a phred-scaled measure of confidence that reads were aligned to the correct genomic location; higher is better.'],
      ['COMP', 'Completeness: the proportion of accessions with a non-missing genotype call at that site.'],
      ['maxR²', 'Maximum linkage-disequilibrium r² used by the dataset filter or imputation-quality workflow; values closer to 1 indicate stronger correlation.'],
      ['MAF', 'Minor-allele frequency among the selected or source accessions, depending on the returned record.'],
      ['PlantCAD1 / PlantCAD2', 'DNA language-model variant scores. PlantCAD2 is available for MaizeGDB 2026 datasets. More extreme disruptive scores are prioritized according to the score convention used by the data pipeline.'],
      ['ESM1 / ESM2 / ESM3', 'Protein language-model scores for amino-acid substitutions. ESM2 and ESM3 are available for MaizeGDB 2026 datasets.'],
      ['Accession genotype columns', 'Each accession column shows its genotype at the site: 0/0 reference homozygous, 0/1 heterozygous, 1/1 alternate homozygous, or ./. missing.'],
      ['Dataset', 'A defined variant collection with its own accession panel, filters, included variant types, and score columns.'],
      ['Sites', 'Number of variant positions in the complete dataset, not necessarily the number returned by the current query.'],
      ['Imputed', 'Whether missing genotypes were statistically inferred in that dataset.'],
    ]},
    { tool:'SNPImpact', color:'#7c3aed', items:[
      ['Gene', 'Gene model associated with the candidate variant.'],
      ['Variant', 'Genomic change, generally shown as position and REF→ALT alleles.'],
      ['Consequence', 'Specific predicted molecular consequence of the change.'],
      ['Domain', 'Pfam domain containing the affected amino acid, when present.'],
      ['PlantCAD1 / PlantCAD2', 'DNA language-model scores used to estimate regulatory or sequence disruption.'],
      ['ESM / ESM2 / ESM3', 'Protein language-model scores used to estimate the effect of an amino-acid substitution.'],
      ['Priority', 'Integrated candidate tier. TOP is the strongest prioritization, followed by HIGH, MODERATE, and LOW.'],
      ['Shortlist / flag', 'A user-selected marker for retaining a candidate variant for later review or export.'],
      ['Gene-model diagram', 'A compact display of exons, introns, coding sequence, strand, and the variant position.'],
      ['Exon', 'A transcript segment retained in the mature RNA; coding portions contribute to the protein sequence.'],
      ['Intron', 'A transcribed segment removed during RNA splicing.'],
      ['CDS', 'Coding sequence: the portion of exons translated into protein.'],
      ['Imputation r²', 'Estimated squared correlation between imputed and true genotypes; values nearer 1 indicate greater confidence.'],
      ['MQ', 'Mapping quality for the variant site.'],
      ['MAF', 'Minor-allele frequency for the candidate variant.'],
    ]},
    { tool:'SNPFunction', color:'#2563eb', items:[
      ['Non-syn : syn', 'The number or ratio of nonsynonymous coding variants to synonymous coding variants. A higher value indicates more amino-acid-changing variation relative to silent variation; it is descriptive and is not by itself a formal dN/dS estimate.'],
      ['Exon : intron', 'The number or ratio of variants in annotated exons to variants in introns of the gene model. Infinity (∞) means exon variants were observed but no intron variants were counted.'],
      ['Domain-disrupting', 'Coding variants that alter an amino acid located inside an annotated Pfam domain.'],
      ['Knockout alleles', 'Alleles predicted to strongly disrupt gene function, such as frameshift, stop-gained, essential splice, or other loss-of-function changes.'],
      ['Mean PlantCAD1 / Mean PlantCAD2', 'Average DNA language-model score across the gene variants included in the burden summary.'],
      ['Mean ESM / Mean ESM2 / Mean ESM3', 'Average protein language-model score across scored amino-acid-changing variants in the gene.'],
      ['Allele', 'The specific genomic REF→ALT change represented by a damaging-allele row.'],
      ['Consequence', 'Predicted molecular effect of that allele.'],
      ['Domain', 'Pfam domain overlapping the affected residue.'],
      ['PlantCAD1 / PlantCAD2', 'DNA language-model score for the allele.'],
      ['ESM / ESM2 / ESM3', 'Protein language-model score for the resulting amino-acid change.'],
      ['Priority', 'Integrated SNPTools evidence tier for the allele.'],
      ['Het', 'Number of accessions carrying the allele heterozygously.'],
      ['Hom', 'Number of accessions carrying the allele as alternate homozygous.'],
      ['AF', 'Alternate-allele frequency across the whole analyzed panel.'],
      ['Variant burden', 'The count and composition of variants assigned to the gene across the full dataset panel.'],
      ['Damaging allele', 'An allele selected because its consequence and/or prediction scores indicate a potentially important functional effect.'],
      ['Functional annotation', 'A per-gene dossier assembled for the 39,756 canonical B73 v5 gene models: identity, description, protein domains, Gene Ontology, pathways, and cross-references.'],
      ['Gene symbol / aliases', 'The primary gene symbol and any additional names or synonyms recorded for the model.'],
      ['Description source', 'Provenance of the functional description, shown as a badge: MaizeGDB (curated), UniProt, InterPro (predicted from domains), or no informative source.'],
      ['Evidence chips', 'An at-a-glance row of what is annotated for the gene: GO, Pfam, KEGG KO, Pathway, UniProt, Symbol.'],
      ['Protein domain architecture', 'A to-scale diagram of the protein with Pfam domains as positioned blocks, plus a list giving Pfam and InterPro IDs, residue span, percent of protein covered, and InterProScan E-value.'],
      ['Gene Ontology (BP / MF / CC)', 'GO terms grouped by aspect — Biological process, Molecular function, Cellular component — each linking out to AmiGO.'],
      ['GO source (MaizeGDB / UniProt / InterPro2GO)', 'Provenance of each GO term. MaizeGDB and UniProt are curated; InterPro2GO is predicted from protein domains.'],
      ['Curated only', 'A toggle that hides GO terms supported only by prediction (InterPro2GO), leaving curated evidence.'],
      ['Obsolete term', 'A GO term whose status is no longer current; it is retained but flagged.'],
      ['KEGG orthology (KO) / Pathways', 'KEGG orthology assignments, pathway memberships, and KEGG gene IDs, where the gene cross-references to UniProt or Entrez.'],
      ['Cross-references', 'External identifiers for the gene: UniProt, NCBI Gene, and B73 v4 / v3 gene models.'],
      ['Annotation build', 'The build date, assembly, and annotation version the functional record was generated from.'],
    ]},
    { tool:'SNPCompare', color:'#0e7490', items:[
      ['#', 'Current rank after sorting and filtering.'],
      ['Project', 'BioProject or dataset project associated with the accession.'],
      ['SRA ID', 'Sequence Read Archive run identifier associated with the accession.'],
      ['Accession name', 'Human-readable line or germplasm name.'],
      ['Global sim', 'Genome-wide identity-by-state similarity between the focal accession and the comparison accession: matching genotypes divided by co-called sites.'],
      ['Local sim', 'Identity-by-state similarity calculated only from variants in the current SNPVersity region.'],
      ['Δ (local−global)', 'Local similarity minus global similarity. Positive values indicate the pair is more similar in the selected region than genome-wide; negative values indicate less similarity in the region.'],
      ['Global miss%', 'Percentage of genome-wide comparison sites where at least one accession lacks a genotype call.'],
      ['Local miss%', 'Percentage of sites in the selected region where at least one accession lacks a genotype call.'],
      ['Similarity', 'Fraction of co-called sites at which the two accessions have the same genotype.'],
      ['Missing%', 'Percentage of evaluated sites where either member of the pair is missing a genotype.'],
      ['Co-called sites', 'Number of sites in the local region with genotype calls for both accessions.'],
      ['Focal accession', 'The accession against which every other accession is ranked.'],
      ['Global', 'Precomputed genome-wide comparison scope.'],
      ['This region / Local', 'Comparison calculated from the current SNPVersity genotype matrix.'],
      ['Both', 'Side-by-side display of global and local values plus their difference.'],
    ]},
    { tool:'SNPTree', color:'#15803d', items:[
      ['IBS allele distance', 'Mean pairwise allele-dosage difference across co-called sites. Identical genotypes contribute 0, opposite homozygotes 1, and a homozygote-versus-heterozygote comparison 0.5.'],
      ['Informative sites', 'Variant sites that contain more than one observed genotype state among the selected accessions.'],
      ['Shared sites', 'Sites with non-missing calls for both accessions in a pair.'],
      ['UPGMA', 'Unweighted Pair Group Method with Arithmetic Mean, an agglomerative clustering method that assumes an ultrametric tree.'],
      ['Neighbour-Joining (NJ)', 'A distance-based tree-building method that does not require equal evolutionary rates among branches.'],
      ['Branch length', 'Distance assigned to a tree edge from the pairwise IBS distance calculation.'],
      ['Newick', 'Compact parenthetical text format for tree topology and branch lengths.'],
      ['MEGA pairwise format', 'Distance-matrix text format that can be imported by MEGA software.'],
      ['PHYLIP', 'Fixed-layout distance-matrix format used by many phylogenetic programs.'],
      ['Local phylogeny', 'A distance tree describing similarity in the selected genomic interval; it should not automatically be interpreted as a whole-genome species tree.'],
    ]},
    { tool:'SNPMatrix', color:'#b45309', items:[
      ['IBS distance', 'Pairwise mean allele-dosage difference across co-called sites; 0 means identical across compared calls and larger values indicate more difference.'],
      ['% identity', 'Similarity view derived from the distance matrix and displayed as a percentage.'],
      ['Shared-site count', 'Number of non-missing sites used for a particular pairwise matrix cell.'],
      ['Input order', 'Rows and columns remain in the same order as the accession selection.'],
      ['Clustered order', 'Rows and columns are reordered by UPGMA clustering to place similar accessions near one another.'],
      ['Bioproject bars', 'Color strips indicating the project or BioProject associated with each accession.'],
      ['Heatmap cell', 'The pairwise distance or identity value for the row accession versus the column accession.'],
      ['CSV distance matrix', 'Comma-separated pairwise IBS distance table.'],
      ['PNG / SVG', 'Raster and scalable-vector image exports of the displayed heatmap.'],
      ['PHYLIP', 'Distance-matrix export for compatible phylogenetic software.'],
    ]},
    { tool:'SNPFold', color:'#be185d', items:[
      ['Variant', 'Coding DNA change mapped to a protein residue.'],
      ['Consequence', 'Predicted molecular effect of the coding change.'],
      ['Residue', 'One-based amino-acid position in the displayed protein sequence.'],
      ['Domain', 'Pfam domain overlapping the residue.'],
      ['Local pLDDT', 'AlphaFold confidence score near the affected residue, from 0 to 100. Higher values indicate greater confidence in the local predicted structure.'],
      ['Structure', 'DSSP-style secondary-structure assignment at the residue: α-helix, β-strand, or loop/coil.'],
      ['PlantCAD / PlantCAD2', 'DNA language-model scores for the underlying nucleotide variant.'],
      ['ESM1 / ESM2 / ESM3', 'Protein language-model scores for the amino-acid substitution.'],
      ['IUPred2', 'Predicted intrinsic disorder at the residue, from 0 to 1 (IUPred2). Higher values mean the residue is more likely to be intrinsically disordered than to adopt a fixed fold.'],
      ['Anchor2', 'Predicted probability (0 to 1) that a disordered residue lies in a protein-binding region (ANCHOR2). Higher values flag disordered segments likely to become ordered upon binding a partner.'],
      ['Activity', 'Annotated functional site overlapping the residue, drawn from InterProScan member databases (for example an active site, binding site, or conserved functional feature), when present.'],
      ['Priority', 'Integrated evidence tier for the structure-mapped variant.'],
      ['Carriers', 'Number of accessions carrying the alternate allele; expanded details separate heterozygous and homozygous carriers.'],
      ['pLDDT', 'Predicted Local Distance Difference Test score from AlphaFold. Common interpretation: ≥90 very high confidence, 70–89 confident, 50–69 low, and <50 very low.'],
      ['Secondary structure', 'Local protein conformation classified as helix, strand, or coil/loop, generated from the structure model.'],
      ['ΔΔG', 'Predicted change in protein folding free energy after mutation. Positive and negative interpretations depend on the scoring convention used by the source model; magnitude reflects predicted structural effect.'],
      ['Protein browser', 'Linear alignment of protein residues, domains, secondary structure, confidence, and variant markers.'],
      ['3D viewer color: confidence', 'Colors residues by pLDDT confidence band.'],
      ['3D viewer color: domain', 'Colors residues by Pfam-domain membership.'],
      ['3D viewer color: impact', 'Colors or highlights variants according to predicted functional severity.'],
      ['3D', 'A per-row control in the variant table that force-shows that residue in the 3D structure viewer, independent of the blanket "Variant residues" toggle.'],
      ['InterProScan sites (by program)', 'Residue-level functional-site and feature annotations from InterProScan, drawn as lanes in the protein browser and grouped by the member database (program) that produced them — CDD, PIRSR, and SFLD. All program lanes are shown even when the current protein has no hits, so the annotation is discoverable.'],
      ['CDD (Conserved Domains Database)', "NCBI's Conserved Domain Database. In the InterProScan site lanes, CDD contributes conserved-domain footprints and the functional-site residues within them."],
      ['PIRSR (PIR Site Rules)', 'PIR Site Rules from the Protein Information Resource. Curated site rules propagate functional-site residues (such as active or binding sites) to proteins that match a site-rule signature.'],
      ['SFLD (Structure–Function Linkage Database)', 'The Structure–Function Linkage Database links protein structural features to specific chemical functions and contributes conserved functional-site residues, chiefly for enzyme superfamilies.'],
    ]},
    { tool:'PanEffect', color:'#b45309', items:[
      ['B73 Reference View', 'The substitution heatmap computed on the B73 v5 protein: every residue (columns) against all 20 possible amino acids (rows).'],
      ['Pan-genome View', 'The same substitution heatmap projected onto a protein multiple-sequence alignment across the maize assemblies, so natural variation and its predicted effect can be read together. Available for canonical transcripts.'],
      ['ESM score', 'A protein-language-model estimate of how tolerated an amino-acid substitution is; more negative is more disruptive.'],
      ['ESM model (ESM1 / ESM2 / ESM3)', 'The protein language model used to score substitutions, selectable in the panel.'],
      ['Substitution heatmap', 'A per-residue grid colored by predicted effect. Hover a cell to read the position, wild-type → substitution, and score.'],
      ['Zoomed region', 'A 50-residue window of the full heatmap, positioned with the slider, showing per-cell substitution letters and the wild-type residue track.'],
      ['MaizeGDB 2026 view', 'Colors the B73 heatmap using substitutions observed in the MaizeGDB 2026 High-Coverage set. Hand-offs from SNPVersity or SNPFold open in this view; "Show all variant effects" colors every possible substitution.'],
      ['PFAM Domains track', 'Pfam domains drawn to scale along the protein and linked to InterPro, shared with the SNPFold and SNPFunction views.'],
      ['Secondary structure', 'Predicted per-residue secondary structure (helix / strand / coil) drawn above the heatmap for context.'],
      ['Species / Gene model', 'In the pan-genome zoomed view, switches the row labels between assembly (species) names and their gene-model IDs.'],
      ['Heterotic group color', 'Pan-genome rows are colored by the assembly\u2019s heterotic group (stiff-stalk, non-stiff-stalk, Iodent, Lancaster, tropical, teosinte, and others).'],
      ['Variant effects file', 'The per-substitution ESM score table for the gene, downloadable from the summary.'],
    ]},
    { tool:'Dataset table', color:'#475569', items:[
      ['Dataset', 'Named variant collection and filtering configuration used for a query.'],
      ['Reference', 'Genome assembly to which reads and variants were aligned.'],
      ['Accessions', 'Number of samples or accession columns available in the dataset.'],
      ['Sites', 'Approximate or reported number of variant sites in the complete dataset.'],
      ['Filters', 'Quality and inclusion rules used to construct the dataset.'],
      ['Het', 'Whether heterozygous genotype calls are retained.'],
      ['INDELs', 'Whether insertion and deletion variants are included.'],
      ['Imputed', 'Whether missing genotype calls have been statistically inferred.'],
    ]},
  ];


  /* ---- FAQ ---- */
  const FAQ = [
    ['Where does the data come from?',
     'A query sends your region and accession list to the server, which reads the real HDF5 variant store, writes a VCF for exactly that slice, and returns it. The tools parse that VCF into the tables and matrices you see, so everything downstream is one consistent result.'],
    ['How do I move a selection between tools?',
     'Run a query in SNPVersity, then use "Send selection to…" in the top bar (or the buttons on a result). The same genotype matrix is handed to SNPImpact, SNPCompare, SNPTree, SNPMatrix, and SNPTrait without re-querying. SNPMatrix and SNPTree can also pass their set on to each other.'],
    ['Why did a wide region give me a download instead of a table?',
     'Regions larger than one million bases skip the in-browser table and return a downloadable VCF instead, so the page stays responsive. Narrow the interval to get the interactive table back.'],
    ['Why are some cells \u2014 or N/A?',
     'A \u2014 in a domain column means Pfam annotation for that position hasn\u2019t been loaded yet. Blank second-generation scores (PlantCAD2, ESM2, ESM3) mean the dataset is not MaizeGDB 2026 — only that family carries them. MQ and coverage read N/A when the source didn\u2019t record them.'],
    ['Is there a limit on how many accessions I can compare?',
     'The distance tools warn before doing very large computations in the browser: SNPTree above 250 accessions and SNPMatrix above 400. You can build anyway, it just may be slow. SNPImpact renders up to 1,500 variants at a time.'],
    ['What can I download?',
     'A VCF from SNPVersity; a CSV distance matrix, PHYLIP, PNG, and SVG from SNPMatrix; Newick, MEGA, and PHYLIP trees from SNPTree. Comparison and impact tables can be exported from their own pages.'],
    ['Which tools are ready to use now?',
     'SNPVersity, SNPImpact, SNPFunction, SNPCompare, SNPTree, SNPMatrix, SNPFold, and PanEffect are live. SNPTrait, SNPImpute, SNPDensity, and SNPGermplasm are on the roadmap and marked in development in the sidebar.'],
  ];

  function esc(s){ return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
  function ico(k){ return (typeof ICONS!=='undefined' && ICONS[k]) || ''; }

  /* ---- dataset table (live from data.js when available) ---- */
  function datasets(){
    try { if (typeof Data!=='undefined' && Data.datasets) return Data.datasets(); } catch(e){}
    return [];
  }
  function flag(on){ return on ? '<span class="hp-yes">Yes</span>' : '<span class="hp-no">\u2013</span>'; }
  function datasetTable(){
    const ds = datasets();
    if (!ds.length) return '';
    const rows = ds.map(d=>`<tr>
        <td><b>${esc(d.name)}</b><div class="hp-sub">${esc(d.sub||'')}</div></td>
        <td class="hp-mono">${esc(d.ref)}</td>
        <td class="hp-num">${esc(d.acc)}</td>
        <td class="hp-num">${esc(d.sites)}</td>
        <td>${(d.filters||[]).map(f=>`<span class="hp-chip">${esc(f)}</span>`).join(' ')}</td>
        <td class="hp-c">${flag(d.het)}</td>
        <td class="hp-c">${flag(d.indel)}</td>
        <td class="hp-c">${flag(d.impute)}</td>
      </tr>`).join('');
    return `<div class="hp-tablewrap"><table class="hp-table">
      <thead><tr>
        <th>Dataset</th><th>Reference</th><th>Accessions</th><th>Sites</th>
        <th>Filters</th><th>Het</th><th>INDELs</th><th>Imputed</th>
      </tr></thead><tbody>${rows}</tbody></table></div>
      <p class="hp-fine">Accessions are grouped into families (MaizeGDB 2026, MaizeGDB 2024, Schnable 2023, NAM 2021). Only the <b>MaizeGDB 2026</b> family carries the second-generation language-model scores (PlantCAD2, ESM2, ESM3); the others provide a single DNA score and a single protein score.</p>`;
  }

  /* ---- in-depth MaizeGDB 2026 dataset description ---- */
  const MGDB2026_PROJECTS = [
    [539,'PRJCA009749','WGS resequencing of 1,604 maize inbred lines','10.1038/s41477-022-01190-2'],
    [521,'PRJNA531553','Deep DNA resequencing of the association mapping panel','10.1038/s41588-019-0427-6'],
    [453,'PRJNA609577','Zea mays genome sequencing','10.1038/s41588-020-0671-9'],
    [340,'PRJNA783885','Maize landrace whole-genome resequencing','10.1038/s41467-022-32180-9'],
    [232,'PRJEB56320','Maize Wisconsin Diversity Panel resequencing','10.1111/tpj.16123'],
    [183,'PRJEB56320','Zea mays sequences (teosinte)','10.1038/s41588-022-01184-y'],
    [77,'PRJNA641489','Maize Nested Association Mapping (NAM)','10.1126/science.abg5289'],
    [76,'PRJEB31061','Maize Haplotype Map version 3','10.1093/gigascience/gix134'],
    [67,'PRJNA399729','Maize landraces from six highland and lowland populations','10.1093/gigascience/gix134'],
    [67,'PRJNA300309','European maize diversity','10.1371/journal.pgen.1006666'],
    [57,'PRJNA783885','Genetic diversity of Zea (teosinte)','10.1038/ng.2313'],
    [49,'PRJNA389800','Whole-genome sequencing of the maize 282 panel','10.1093/gigascience/gix134'],
    [33,'PRJNA260788','European maize genomes','10.1038/s41588-020-0671-9'],
    [7,'PRJNA479960','South American maize genome sequencing','10.1126/science.aav0207'],
    [4,'PRJEB32225','Zm-B73-REFERENCE-NAM-5.0','10.1126/science.abg5289'],
    [3,'PRJEB56265','Resequencing of three Polish maize inbred lines','10.1111/tpj.16123'],
    [1,'PRJEB61159','Coastal preceramic maize from Paredones, Peru','10.7554/eLife.83149'],
    [1,'PRJNA352392','A 5,310-year-old maize cob from the Tehuacan Valley, Mexico','10.1016/j.cub.2016.09.036']
  ];

  const MGDB2026_EFFECTS = [
    ['Intergenic',274890726,90284097],
    ["5' UTR",836712,398940],
    ['Synonymous',1409639,670709],
    ['Missense',1767459,747451],
    ['Stop-related',84653,26405],
    ['Frameshift',222401,71174],
    ['Intron',12043697,5916121],
    ['Non-coding',6777,3133],
    ["3' UTR",1231422,640830],
    ['Other',120889,46629],
    ['Total',290043644,97572347]
  ];

  const MGDB2026_CHR = [
    ['Chr1',308452471,51382911,41928591,14262099],
    ['Chr2',243675191,40457547,32643386,11008991],
    ['Chr3',238017767,39598148,32693292,10861052],
    ['Chr4',250330460,42142813,35830675,12496251],
    ['Chr5',226353449,36679188,30747268,10298954],
    ['Chr6',181357234,28914586,23742034,8001535],
    ['Chr7',185808916,30391975,24862777,8231456],
    ['Chr8',182411202,29883638,24964435,8054636],
    ['Chr9',163004744,26922259,21747152,7278487],
    ['Chr10',152435371,25319421,20883764,7078616],
    ['Total',2131846805,351692486,290043374,97572077]
  ];

  function fmtInt(n){ return Number(n).toLocaleString('en-US'); }

  function maize2026Details(){
    const projects = MGDB2026_PROJECTS.map(r=>`<tr>
      <td class="hp-num">${fmtInt(r[0])}</td>
      <td class="hp-mono">${esc(r[1])}</td>
      <td>${esc(r[2])}</td>
      <td class="hp-mono">${esc(r[3])}</td>
    </tr>`).join('');

    const effects = MGDB2026_EFFECTS.map((r,i)=>`<tr${i===MGDB2026_EFFECTS.length-1?' class="hp-total"':''}>
      <td>${esc(r[0])}</td><td class="hp-num">${fmtInt(r[1])}</td><td class="hp-num">${fmtInt(r[2])}</td>
    </tr>`).join('');

    const chr = MGDB2026_CHR.map((r,i)=>`<tr${i===MGDB2026_CHR.length-1?' class="hp-total"':''}>
      <td class="hp-mono">${esc(r[0])}</td><td class="hp-num">${fmtInt(r[1])}</td>
      <td class="hp-num">${fmtInt(r[2])}</td><td class="hp-num">${fmtInt(r[3])}</td><td class="hp-num">${fmtInt(r[4])}</td>
    </tr>`).join('');

    return `
        <br>
      <details class="hp-data-detail" open>
        <summary><b>MaizeGDB 2026 dataset: composition, processing, and scale</b><span class="hp-chev">${ico('caret')}</span></summary>
        <div class="hp-data-body">
          <p>The MaizeGDB 2026 resource combines public whole-genome resequencing from <b>2,710 maize accessions</b>, including diverse inbred lines, landraces, teosintes, association panels, NAM founders and related materials, and several historically important samples. All reads were processed through a standardized variant-calling workflow against <b>B73 RefGen_v5</b>, so coordinates, reference alleles, gene models, and downstream annotations use one common reference system.</p>

          <h3>How the two MaizeGDB 2026 datasets differ</h3>
          <div class="hp-comparegrid">
            <div class="hp-dcard"><b>High Coverage</b><p>Contains approximately 290 million loci that passed mapping-quality and genotypic-coverage requirements. This version retains broader variation and is useful when sensitivity and variant discovery are the main goals.</p></div>
            <div class="hp-dcard"><b>High Quality</b><p>Contains approximately 98 million loci. It applies the same mapping-quality and coverage filters plus an additional high-confidence linkage-disequilibrium criterion. This more conservative set is useful when specificity and confidence are priorities.</p></div>
          </div>
          <p class="hp-fine">“High Coverage” describes the broader filtered set; “High Quality” is the stricter subset. A locus count refers to a genomic variant position in the complete dataset, not the number of rows returned for a particular region or accession selection.</p>

          <h3>What is stored and displayed</h3>
          <p>The underlying resource is maintained in VCF and HDF5 forms. SNPVersity sends the selected dataset, genomic interval, and accession list to the server, which extracts the requested slice and returns a VCF. SNPTools then reuses that same genotype matrix in SNPImpact, SNPCompare, SNPTree, SNPMatrix, and related pages. Each site may include REF and ALT alleles, accession genotypes, predicted molecular consequence, gene association, mapping quality, genotype completeness, linkage-disequilibrium support, allele frequency, Pfam domain context, and DNA- or protein-language-model scores when available.</p>

          <h3>Source projects</h3>
          <p>The panel was assembled from multiple public projects rather than one experiment. This increases biological and geographic diversity, but it also means sequencing depth, library preparation, and project design can differ among accessions. The standardized alignment, calling, and filtering workflow reduces—though does not completely remove—these study-to-study differences.</p>
          <div class="hp-tablewrap"><table class="hp-table hp-compact">
            <thead><tr><th>Accessions</th><th>BioProject</th><th>Project or population</th><th>DOI</th></tr></thead>
            <tbody>${projects}</tbody>
          </table></div>

          <h3>Variant-effect composition</h3>
          <p>Most loci are intergenic because much of the maize genome lies outside annotated coding regions. Coding and gene-associated categories are much smaller but are especially important for SNPImpact, SNPFunction, and SNPFold. “Stop-related” summarizes variants annotated in the source statistics as stop effects; the exact transcript-level consequence shown in a result may be more specific.</p>
          <div class="hp-tablewrap"><table class="hp-table hp-compact">
            <thead><tr><th>Effect category</th><th>High Coverage</th><th>High Quality</th></tr></thead>
            <tbody>${effects}</tbody>
          </table></div>

          <h3>Distribution by chromosome</h3>
          <p>Variant counts broadly track chromosome length, although local diversity, repetitive sequence, mappability, selection, and the composition of the accession panel also affect density. “Raw” is the pre-filter total reported in the source summary; the two filtered columns show the successive retained sets.</p>
          <div class="hp-tablewrap"><table class="hp-table hp-compact">
            <thead><tr><th>Chromosome</th><th>Length (bp)</th><th>Raw variants</th><th>High Coverage</th><th>High Quality</th></tr></thead>
            <tbody>${chr}</tbody>
          </table></div>

          <h3>Important interpretation notes</h3>
          <ul class="hp-notes">
            <li><b>Reference-relative calls:</b> REF and ALT are defined relative to B73 v5; “alternate” does not mean rare, harmful, or derived.</li>
            <li><b>Accessions are not all independent:</b> some projects include replicates, related lines, founders, or multiple runs. Project and accession metadata should be considered when interpreting similarity.</li>
            <li><b>Missingness varies:</b> genotype completeness can differ by site and accession. Pairwise tools exclude sites where either member lacks a usable call.</li>
            <li><b>Consequence is transcript dependent:</b> one genomic variant may receive several annotations when it overlaps multiple transcripts; interfaces generally display the most severe or most relevant consequence.</li>
            <li><b>Model scores are predictions:</b> PlantCAD and ESM scores prioritize candidates but do not establish biological causality. Use them together with frequency, consequence, domain, structure, phenotype, and experimental evidence.</li>
            <li><b>Counts can differ slightly among summaries:</b> totals generated at different pipeline stages or from differently normalized records may vary by a small number of sites. The live dataset metadata and downloadable files are authoritative for an analysis.</li>
          </ul>
        </div>
      </details>`;
  }

  function statusPill(s){
    const live = s===LIVE;
    return `<span class="hp-status ${live?'on':'off'}">${live?'Live':'In development'}</span>`;
  }
  function pageCard(p){
    const live = p.status===LIVE;
    const openBtn = live
      ? `<button class="hp-open" onclick="go('${p.id}')">Open ${esc(p.name)} ${ico('caret')}</button>`
      : `<span class="hp-open muted" title="On the roadmap">Coming soon</span>`;
    return `<details class="hp-tool">
      <summary>
        <span class="hp-ti" style="background:${p.color}">${ico(p.icon)}</span>
        <span class="hp-tt"><b>${esc(p.name)}</b><span class="hp-tag">${esc(p.tag)}</span></span>
        ${statusPill(p.status)}
        <span class="hp-chev">${ico('caret')}</span>
      </summary>
      <div class="hp-tbody">
        <p>${esc(p.what)}</p>
        <div class="hp-io">
          <div><span class="hp-lbl">You give it</span>${esc(p.give)}</div>
          <div><span class="hp-lbl">You get back</span>${esc(p.get)}</div>
        </div>
        <div class="hp-actions">${openBtn}</div>
      </div>
    </details>`;
  }


  function definitionGroup(g){
    const items = g.items.map(x=>`<div class="hp-def"><dt>${esc(x[0])}</dt><dd>${esc(x[1])}</dd></div>`).join('');
    return `<details class="hp-defgroup">
      <summary><span class="hp-defdot" style="background:${g.color}"></span><b>${esc(g.tool)}</b><span class="hp-count">${g.items.length} definitions</span><span class="hp-chev">${ico('caret')}</span></summary>
      <dl class="hp-deflist">${items}</dl>
    </details>`;
  }

  function render(page){
    injectCSS();
    page = page || document.getElementById('page');
    page.className = 'page fade';
    const crumb=document.getElementById('crumbTool'); if(crumb) crumb.innerHTML='<b>Help &amp; FAQ</b>';

    const toolCards = PAGES.map(pageCard).join('');
    const gloss = GLOSSARY.map(g=>`<div class="hp-gl"><dt>${esc(g[0])}</dt><dd>${esc(g[1])}</dd></div>`).join('');
    const definitions = DEFINITIONS.map(definitionGroup).join('');
    const faq = FAQ.map(f=>`<details class="hp-faq"><summary>${esc(f[0])}<span class="hp-chev">${ico('caret')}</span></summary><div>${esc(f[1])}</div></details>`).join('');

    page.innerHTML = `
      <section class="hp-hero">
        <div class="hp-eyebrow">SNPTools · SNPVersity 2.1 · B73 RefGen v5</div>
        <h1>Help &amp; FAQ</h1>
        <p>SNPTools is an integrated suite for exploring maize sequence variation. Everything starts from a genomic
           query and flows between tools without re-running it — this page explains each tool, the datasets behind
           them, and the vocabulary you\u2019ll meet along the way.</p>
        <div class="hp-jump">
          <a href="#hp-flow">How it works</a>
          <a href="#hp-tools">The tools</a>
          <a href="#hp-data">Datasets</a>
          <a href="#hp-definitions">Definitions</a>
          <a href="#hp-gloss">Scores &amp; annotations</a>
          <a href="#hp-faq">FAQ</a>
        </div>
      </section>

      <section id="hp-flow" class="hp-sec">
        <div class="hp-h"><span class="hp-n">01</span><h2>How the suite fits together</h2></div>
        <div class="hp-flow">
          <div class="hp-step"><span class="hp-si" style="background:#2563eb">${ico('search')}</span>
            <b>Query</b><p>In SNPVersity, choose a dataset, a region or gene, and the accessions you care about.</p></div>
          <div class="hp-arrow">${ico('caret')}</div>
          <div class="hp-step"><span class="hp-si" style="background:#1f8a4c">${ico('table')}</span>
            <b>Result</b><p>You get a genotype table and a VCF for exactly that slice of the genome.</p></div>
          <div class="hp-arrow">${ico('caret')}</div>
          <div class="hp-step"><span class="hp-si" style="background:#b45309">${ico('compare')}</span>
            <b>Send onward</b><p>Hand the same matrix to any other tool with "Send selection to…" — no re-query.</p></div>
        </div>
        <p class="hp-fine">A query is run against the real variant store and returned as a VCF; every other tool reuses that
          one result, so a set you build once stays consistent as you rank it, compare it, cluster it, or draw it.</p>
      </section>

      <section id="hp-tools" class="hp-sec">
        <div class="hp-h"><span class="hp-n">02</span><h2>The tools</h2></div>
        <p class="hp-lead">Twelve tools across five stages. Eight are live today; the rest are on the roadmap and share
          the same data and coordinates. Expand any tool for what it does, what it takes, and what it returns.</p>
        <div class="hp-tools">${toolCards}</div>
      </section>

      <section id="hp-data" class="hp-sec">
        <div class="hp-h"><span class="hp-n">03</span><h2>Datasets</h2></div>
        <p class="hp-lead">Each query runs against one dataset. All are called against the B73 v5 reference; they differ
          in how they were filtered, how many accessions and sites they hold, and which score columns they carry.</p>
        ${datasetTable()}
        ${maize2026Details()}
      </section>

      <section id="hp-definitions" class="hp-sec">
        <div class="hp-h"><span class="hp-n">04</span><h2>Definitions &amp; table columns</h2></div>
        <p class="hp-lead">Definitions are grouped by tool and use the same labels shown in the interfaces. Expand a group or use your browser's find command to locate a column heading.</p>
        <div class="hp-defgroups">${definitions}</div>
      </section>

      <section id="hp-gloss" class="hp-sec">
        <div class="hp-h"><span class="hp-n">05</span><h2>Scores &amp; annotations</h2></div>
        <dl class="hp-gloss">${gloss}</dl>
      </section>

      <section id="hp-faq" class="hp-sec">
        <div class="hp-h"><span class="hp-n">06</span><h2>Frequently asked</h2></div>
        <div class="hp-faqs">${faq}</div>
      </section>

      <section class="hp-foot">
        <div>
          <b>Still stuck?</b>
          <p>SNPTools is part of MaizeGDB, the Maize Genetics and Genomics Database.</p>
        </div>
        <div class="hp-foot-btns">
          <button class="hp-open" onclick="go('snpversity')">Start in SNPVersity ${ico('caret')}</button>
          <a class="hp-open ghost" href="https://maizegdb.org" target="_blank" rel="noopener">Visit MaizeGDB ${ico('caret')}</a>
        </div>
      </section>`;

    // jump links: smooth scroll within the .scroll container
    page.querySelectorAll('.hp-jump a').forEach(a=>{
      a.addEventListener('click', e=>{
        const t=document.querySelector(a.getAttribute('href')); if(!t) return;
        e.preventDefault(); t.scrollIntoView({behavior:'smooth', block:'start'});
      });
    });
  }

  /* Entry point used by the MaizeGDB masthead link. Mirrors go() but doesn't
     assume 'help' is present in the TOOLS map, so core.js needs no changes. */
  function open(){
    if (typeof S!=='undefined') S.tool='help';
    if (typeof closeMenu==='function') closeMenu();
    if (typeof toggleRail==='function') toggleRail(false);
    const sb=document.getElementById('sendBtn'); if(sb) sb.style.display='none';
    if (typeof renderNav==='function') renderNav();
    document.querySelectorAll('.mg-help').forEach(el=>el.classList.add('on'));
    window.scrollTo(0,0);
    const sc=document.querySelector('.scroll'); if(sc) sc.scrollTop=0;
    render(document.getElementById('page'));
  }

  function injectCSS(){
    if (document.getElementById('snphelp-css')) return;
    const s=document.createElement('style'); s.id='snphelp-css';
    s.textContent = `
      .hp-hero{padding:6px 0 18px;border-bottom:1px solid var(--line,#e6e9ef);margin-bottom:26px}
      .hp-eyebrow{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--blue-600,#2563eb);margin-bottom:10px}
      .hp-hero h1{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:40px;line-height:1.05;letter-spacing:-.02em;margin:0 0 12px;color:var(--ink,#141922)}
      .hp-hero p{max-width:70ch;color:var(--muted,#5b6b83);font-size:15px;line-height:1.6;margin:0 0 18px}
      .hp-jump{display:flex;gap:8px;flex-wrap:wrap}
      .hp-jump a{font-size:12.5px;font-weight:600;color:var(--ink,#141922);text-decoration:none;padding:7px 13px;border:1px solid var(--line,#e6e9ef);border-radius:999px;background:#fff;transition:.15s}
      .hp-jump a:hover{border-color:var(--blue-600,#2563eb);color:var(--blue-600,#2563eb)}

      .hp-sec{margin:0 0 42px;scroll-margin-top:16px}
      .hp-h{display:flex;align-items:baseline;gap:12px;margin-bottom:16px}
      .hp-h h2{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:23px;letter-spacing:-.01em;margin:0;color:var(--ink,#141922)}
      .hp-n{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:12px;font-weight:600;color:#c8a24a;border:1px solid #e7d6a8;background:#fdf8ec;border-radius:6px;padding:3px 7px}
      .hp-lead{max-width:74ch;color:var(--muted,#5b6b83);font-size:14px;line-height:1.6;margin:0 0 16px}
      .hp-fine{max-width:78ch;color:var(--muted,#5b6b83);font-size:12.5px;line-height:1.6;margin:12px 0 0}
      .hp-fine b,.hp-lead b{color:var(--ink,#141922)}

      /* flow */
      .hp-flow{display:flex;align-items:stretch;gap:10px;flex-wrap:wrap}
      .hp-step{flex:1 1 210px;border:1px solid var(--line,#e6e9ef);border-radius:12px;background:#fff;padding:16px}
      .hp-step b{display:block;font-family:var(--disp,'Space Grotesk',sans-serif);font-size:15px;color:var(--ink,#141922);margin:10px 0 5px}
      .hp-step p{margin:0;font-size:13px;line-height:1.55;color:var(--muted,#5b6b83)}
      .hp-si{display:inline-flex;width:34px;height:34px;border-radius:9px;align-items:center;justify-content:center;color:#fff}
      .hp-si svg{width:19px;height:19px}
      .hp-arrow{display:flex;align-items:center;color:var(--faint,#aab4c4)}
      .hp-arrow svg{width:22px;height:22px}

      /* tool accordions */
      .hp-tools{display:flex;flex-direction:column;gap:10px}
      .hp-tool{border:1px solid var(--line,#e6e9ef);border-radius:12px;background:#fff;overflow:hidden}
      .hp-tool[open]{border-color:#cdd6e6;box-shadow:0 1px 0 rgba(20,25,34,.03)}
      .hp-tool summary{display:flex;align-items:center;gap:13px;padding:14px 16px;cursor:pointer;list-style:none}
      .hp-tool summary::-webkit-details-marker{display:none}
      .hp-ti{display:inline-flex;width:34px;height:34px;border-radius:9px;align-items:center;justify-content:center;color:#fff;flex:0 0 auto}
      .hp-ti svg{width:19px;height:19px}
      .hp-tt{display:flex;flex-direction:column;gap:2px;min-width:0}
      .hp-tt b{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:15.5px;color:var(--ink,#141922)}
      .hp-tag{font-size:12.5px;color:var(--muted,#5b6b83)}
      .hp-status{margin-left:auto;font-size:10.5px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;padding:3px 9px;border-radius:999px;border:1px solid}
      .hp-status.on{color:#1f8a4c;border-color:#bfe3cc;background:#f0faf3}
      .hp-status.off{color:#b06a12;border-color:#f0dcbb;background:#fdf6ec}
      .hp-chev{color:var(--faint,#aab4c4);transition:transform .18s;flex:0 0 auto}
      .hp-chev svg{width:18px;height:18px;transform:rotate(90deg)}
      .hp-tool[open] .hp-chev svg{transform:rotate(-90deg)}
      .hp-tbody{padding:2px 16px 18px 63px}
      .hp-tbody>p{margin:0 0 14px;font-size:13.5px;line-height:1.62;color:#3a465a;max-width:76ch}
      .hp-io{display:flex;gap:24px;flex-wrap:wrap;padding:13px 15px;background:var(--blue-50,#f2f6ff);border:1px solid #e4ecfb;border-radius:10px}
      .hp-io>div{flex:1 1 240px;font-size:13px;line-height:1.5;color:#3a465a}
      .hp-lbl{display:block;font-size:10px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:var(--blue-600,#2563eb);margin-bottom:4px}
      .hp-actions{margin-top:14px}
      .hp-open{display:inline-flex;align-items:center;gap:6px;font-family:inherit;font-size:13px;font-weight:600;color:#fff;background:var(--ink,#141922);border:1px solid var(--ink,#141922);border-radius:9px;padding:9px 15px;cursor:pointer;text-decoration:none;transition:.15s}
      .hp-open:hover{background:#000}
      .hp-open svg{width:15px;height:15px}
      .hp-open.ghost{background:#fff;color:var(--ink,#141922)}
      .hp-open.ghost:hover{border-color:var(--blue-600,#2563eb);color:var(--blue-600,#2563eb)}
      .hp-open.muted{background:#f4f6fa;border-color:var(--line,#e6e9ef);color:var(--muted,#5b6b83);cursor:default}

      /* dataset table */
      .hp-tablewrap{overflow-x:auto;border:1px solid var(--line,#e6e9ef);border-radius:12px;background:#fff}
      .hp-table{border-collapse:collapse;width:100%;font-size:13px;min-width:720px}
      .hp-table th{text-align:left;font-size:10.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted,#5b6b83);font-weight:700;padding:11px 14px;border-bottom:1px solid var(--line,#e6e9ef);background:#fafbfd;white-space:nowrap}
      .hp-table td{padding:12px 14px;border-bottom:1px solid #eef1f6;vertical-align:top;color:#3a465a}
      .hp-table tr:last-child td{border-bottom:none}
      .hp-sub{font-size:11.5px;color:var(--muted,#5b6b83);margin-top:2px;font-weight:400}
      .hp-mono{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:12px;white-space:nowrap}
      .hp-num{font-family:var(--mono,'IBM Plex Mono',monospace);font-size:12.5px;white-space:nowrap;color:var(--ink,#141922)}
      .hp-c{text-align:center}
      .hp-chip{display:inline-block;font-size:11px;color:#3a465a;background:#f1f4f9;border:1px solid #e4e9f1;border-radius:6px;padding:2px 7px;margin:2px 2px 0 0;white-space:nowrap}
      .hp-yes{color:#1f8a4c;font-weight:600}
      .hp-no{color:var(--faint,#aab4c4)}


      /* complete definitions */
      .hp-defgroups{display:flex;flex-direction:column;gap:9px}
      .hp-defgroup{border:1px solid var(--line,#e6e9ef);border-radius:11px;background:#fff;overflow:hidden}
      .hp-defgroup[open]{border-color:#cdd6e6}
      .hp-defgroup summary{display:flex;align-items:center;gap:10px;padding:13px 15px;cursor:pointer;list-style:none}
      .hp-defgroup summary::-webkit-details-marker{display:none}
      .hp-defgroup summary b{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:15px;color:var(--ink,#141922)}
      .hp-defdot{width:10px;height:10px;border-radius:50%;flex:0 0 auto}
      .hp-count{font-size:11.5px;color:var(--muted,#5b6b83);margin-left:2px}
      .hp-defgroup summary .hp-chev{margin-left:auto}
      .hp-deflist{margin:0;padding:0 15px 8px 35px;display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:0 28px}
      .hp-def{padding:11px 0;border-top:1px solid #eef1f6}
      .hp-def dt{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:13.5px;font-weight:600;color:var(--ink,#141922);margin-bottom:3px}
      .hp-def dd{margin:0;font-size:12.8px;line-height:1.55;color:var(--muted,#5b6b83)}


      /* glossary */
      .hp-gloss{margin:0;display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:2px 26px}
      .hp-gl{padding:13px 0;border-bottom:1px solid #eef1f6}
      .hp-gl dt{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:14px;font-weight:600;color:var(--ink,#141922);margin-bottom:4px}
      .hp-gl dd{margin:0;font-size:13px;line-height:1.55;color:var(--muted,#5b6b83)}

      /* faq */
      .hp-faqs{display:flex;flex-direction:column;gap:8px}
      .hp-faq{border:1px solid var(--line,#e6e9ef);border-radius:10px;background:#fff}
      .hp-faq summary{display:flex;align-items:center;gap:12px;padding:14px 16px;cursor:pointer;font-family:var(--disp,'Space Grotesk',sans-serif);font-size:14.5px;font-weight:600;color:var(--ink,#141922);list-style:none}
      .hp-faq summary::-webkit-details-marker{display:none}
      .hp-faq .hp-chev{margin-left:auto}
      .hp-faq>div{padding:0 16px 16px;font-size:13.5px;line-height:1.62;color:#3a465a;max-width:80ch}

      /* footer */
      .hp-foot{display:flex;align-items:center;gap:18px;flex-wrap:wrap;justify-content:space-between;border:1px solid var(--line,#e6e9ef);border-radius:14px;background:linear-gradient(180deg,#fff,#fbfcfe);padding:20px 22px;margin-bottom:24px}
      .hp-foot b{font-family:var(--disp,'Space Grotesk',sans-serif);font-size:16px;color:var(--ink,#141922)}
      .hp-foot p{margin:4px 0 0;font-size:13px;color:var(--muted,#5b6b83)}
      .hp-foot-btns{display:flex;gap:10px;flex-wrap:wrap}

      @media (max-width:620px){
        .hp-hero h1{font-size:31px}
        .hp-tbody{padding-left:16px}
        .hp-arrow{transform:rotate(90deg);align-self:center}
      }`;
    document.head.appendChild(s);
  }

  if (typeof SNPTools!=='undefined') SNPTools.register('help', { render });
  return { render, open };
})();
if (typeof window!=='undefined'){ window.SNPHelp = SNPHelp; window.openHelp = SNPHelp.open; }
