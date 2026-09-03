"""
h5_to_vcf.py  --  region + accession list  ->  VCF, pulled from the HDF5 store.

    python h5_to_vcf.py <db.h5> <out.vcf> <start> <end> <accessions>

<accessions> is either a path to a JSON file holding an array of accession
ids (how processForm.php passes it, to dodge Windows shell-quoting of '"'),
or a raw JSON array string (any direct caller).

The heavy lifting is vectorized: genotype codes (stored as int8) are mapped
to VCF strings through a NumPy lookup table over the whole variant x
accession block at once, and rows are assembled and written in chunks. The
per-cell Python loop it replaced was O(variants x accessions) with a
function call + dict lookup each, which is what made wide / full-panel
queries take minutes. POS is decoded with a single astype() rather than a
Python loop over every position in the file.

The INFO column is stored verbatim in the HDF5 and passed straight through,
so every annotation key -- GENEMODEL, TYPE, EFFECT, SUB, MQ, CVP, MAXR2,
MAF, plantcad1/2_score, ESM1/2/3_score -- survives unchanged.
"""
import gzip
import h5py
import os
import sys
import json
import datetime
import numpy as np

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
hdf5_file_path = sys.argv[1]
output_vcf_path = sys.argv[2]
lower_bound = int(sys.argv[3])
upper_bound = int(sys.argv[4])

arg5 = sys.argv[5]
if os.path.isfile(arg5):
    with open(arg5, 'r') as _f:
        json_string = _f.read()
else:
    json_string = arg5
genome_list = json.loads(json_string)

# HDF5 genotype code -> VCF genotype string. Anything outside 0..3 is clamped
# to 3 ('./.'), matching the old dict-with-default behaviour.
GT_LUT = np.array([b'0/0', b'1/0', b'1/1', b'./.'], dtype='S3')
FIXED_COLS = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO']
ROW_CHUNK = 20000
current_date = datetime.date.today().strftime('%Y%m%d')


def as_bytes_col(arr):
    """h5 column (object-of-bytes, 'S', or 'U') -> contiguous fixed-width 'S' array."""
    if arr.dtype.kind == 'U':
        return np.char.encode(arr)
    return arr.astype('S')


# ---------------------------------------------------------------------------
# Read the requested slice out of the HDF5 store
# ---------------------------------------------------------------------------
with h5py.File(hdf5_file_path, 'r') as hdf5_file:
    if 'POS' not in hdf5_file:
        print("No 'POS' dataset found in the file.")
        sys.exit(1)

    pos_raw = hdf5_file['POS'][:]
    try:
        pos_data = pos_raw.astype(np.int64)
    except (ValueError, TypeError):
        pos_data = np.array(
            [int(p.decode('utf-8')) if isinstance(p, bytes) else int(p) for p in pos_raw],
            dtype=np.int64,
        )

    lower_index = int(np.searchsorted(pos_data, lower_bound, side='left'))
    upper_index = int(np.searchsorted(pos_data, upper_bound, side='right'))
    n_rows = upper_index - lower_index

    if n_rows == 0:
        print("No data found in the specified position range.")
        sys.exit(0)

    fixed_cols = {c: as_bytes_col(hdf5_file[c][lower_index:upper_index]) for c in FIXED_COLS}

    # Only read columns that actually exist in this file; a stray/foreign
    # accession id in the request is filled with './.' rather than aborting.
    present = [g for g in genome_list if g in hdf5_file]
    missing = [g for g in genome_list if g not in hdf5_file]
    if missing:
        print(f"Note: {len(missing)} requested accession(s) not in {hdf5_file_path}; "
              f"filling with ./.  e.g. {missing[:5]}")

    # Genotype code matrix (n_rows, n_accessions) in request order.
    gt_codes = np.empty((n_rows, len(genome_list)), dtype=np.int8)
    for k, g in enumerate(genome_list):
        gt_codes[:, k] = hdf5_file[g][lower_index:upper_index] if g in hdf5_file else 3
    np.clip(gt_codes, 0, 3, out=gt_codes)

n_acc = len(genome_list)
ID_COL = np.full(n_rows, b'.', dtype='S1')
FILTER_COL = np.full(n_rows, b'.', dtype='S1')
FORMAT_COL = np.full(n_rows, b'GT', dtype='S2')


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
def common_info_header():
    return ["##fileformat=VCFv4.2", "##fileDate=" + current_date]


def info_defs():
    return [
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">',
        '##INFO=<ID=CVC,Number=1,Type=Integer,Description="The number of accessions that have genotype data for a particular variant">',
        '##INFO=<ID=CVP,Number=1,Type=Float,Description="The percent of accessions that have genotype data for a particular variant.">',
        '##INFO=<ID=TYPE,Number=.,Type=String,Description="The type of effect using Sequence Ontology terms">',
        '##INFO=<ID=EFFECT,Number=.,Type=String,Description="An estimation of putative impact/deleteriousness">',
        '##INFO=<ID=GENEMODEL,Number=.,Type=String,Description="The name of the gene model affected by the variant">',
        '##INFO=<ID=SUB,Number=.,Type=String,Description="The amino acid substitution for missense and non-synonymous variants">',
        '##INFO=<ID=MAXR2,Number=1,Type=Float,Description="The maximum R2 for a given loci">',
        '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor Allele Frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + genome_list),
    ]


def header_text():
    lines = []
    if "maizegdb" in hdf5_file_path:
        lines += common_info_header()
        lines.append("##source=MaizeGDB")
        lines.append("##reference=Andorf CM, Ross-Ibarra J, Seetharam AS, Hufford MB, Woodhouse MR. (2024) A unified VCF data set from nearly 1,500 diverse maize accessions and resources to explore the genomic landscape of maize. G3 Genes|Genomes|Genetics.")
        lines.append("##doi=https://doi.org/10.1101/2024.04.30.591904")
        lines += info_defs()
    if "nam" in hdf5_file_path:
        lines += common_info_header()
        lines.append("##source=MaizeGDB+NAM")
        lines.append("##reference=Hufford MB, Seetharam AS, Woodhouse MR, et al. De novo assembly, annotation, and comparative analysis of 26 diverse maize genomes. Science. 2021;373(6555):655-662.")
        lines.append("##doi=https://doi.org/10.1126/science.abg5289")
        lines += info_defs()
    if "schnable" in hdf5_file_path:
        lines += common_info_header()
        lines.append("##source=Schnable2023")
        lines.append("##reference=Grzybowski MW, Mural RV, Xu G, Turkus J, Yang J, Schnable JC. A common resequencing-based genetic marker data set for global maize diversity. Plant J. 2023;113(6):1109-1121.")
        lines.append("##doi=https://doi.org/10.1111/tpj.16123")
        lines += info_defs()
    return ("\n".join(lines) + "\n").encode() if lines else b""


# ---------------------------------------------------------------------------
# Write the VCF, one row-chunk at a time. A '.gz' output path is gzip-
# compressed transparently. compresslevel=1: the genotype text is ~90% '0/0',
# so LZ77 alone gets ~10x and level 1 is ~30x faster than the default level 9
# for a barely-larger file (a 250 MB VCF: L1 = 25 MB in ~1s, L9 = 11 MB in ~38s).
# h5_to_vcf.py called directly with a plain '.vcf' path still writes plain text.
# ---------------------------------------------------------------------------
def _open(path, mode):
    return gzip.open(path, mode, compresslevel=1) if path.endswith('.gz') else open(path, mode)
with _open(output_vcf_path, 'wb') as vcf_file:
    vcf_file.write(header_text())

    for a in range(0, n_rows, ROW_CHUNK):
        b = min(a + ROW_CHUNK, n_rows)
        m = b - a

        # 9 fixed columns, tab-joined: only 9 vectorized adds.
        line = fixed_cols['CHROM'][a:b]
        for col in (fixed_cols['POS'][a:b], ID_COL[a:b], fixed_cols['REF'][a:b],
                    fixed_cols['ALT'][a:b], fixed_cols['QUAL'][a:b], FILTER_COL[a:b],
                    fixed_cols['INFO'][a:b], FORMAT_COL[a:b]):
            line = np.char.add(np.char.add(line, b'\t'), col)

        # Genotype block: one vectorized LUT map over the whole (m, n_acc)
        # slice, prefix each cell with a tab, then reinterpret each row's
        # contiguous 4-byte cells ('\t' + 'x/y') as a single string.
        if n_acc:
            cells = np.ascontiguousarray(np.char.add(b'\t', GT_LUT[gt_codes[a:b]]))
            gt_lines = cells.view(f'S{4 * n_acc}').reshape(m)
            line = np.char.add(line, gt_lines)

        vcf_file.write(b'\n'.join(line.tolist()))
        vcf_file.write(b'\n')

print(f"variants: {n_rows}")
print(f"VCF data has been saved to {output_vcf_path}")
