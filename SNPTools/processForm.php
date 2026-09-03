<?php
/* =====================================================================
 *  processForm.php — region + accession list  ->  VCF (via h5_to_vcf.py)
 *
 *  Local testing (MAMP) vs production (Linux/Docker): nothing in this file
 *  changes between machines. The Python interpreter is resolved from the
 *  PYTHON_PATH environment variable, falling back to 'python3' on PATH.
 *  Python deps are declared in requirements.txt (pip-installed at image build).
 * ===================================================================== */

header('Content-Type: application/json');

/* Building a VCF for a wide region / large accession set legitimately takes
 * more than PHP's default 30 s max_execution_time (which on Windows counts
 * the blocking shell_exec below as wall-clock). Lift it so the request runs
 * to completion instead of dying with a Fatal error mid-stream. */
set_time_limit(0);

/* ---------------------------------------------------------------------
 *  CONFIG — provided by the environment; no machine paths committed to git
 * ------------------------------------------------------------------- */
// 1) Python interpreter that has h5py + numpy installed. Resolution order:
//      a. PYTHON_PATH environment variable (first hit wins). Set it per host,
//         not in code, so nothing here is machine-specific:
//           - Docker:  ENV PYTHON_PATH=python3   (or just rely on the fallback)
//           - MAMP:    SetEnv PYTHON_PATH /path/to/venv/bin/python  in an Apache
//                      conf / .htaccess, or export it in the shell that starts MAMP
//                      — keeps your local absolute path out of the repo.
//      b. 'python3' on the system PATH (the default). In the Docker image the
//         deps from requirements.txt are pip-installed globally, so python3 works
//         out of the box with no configuration.
$PYTHON_PATH = getenv('PYTHON_PATH');
if (!$PYTHON_PATH) {
    $PYTHON_PATH = 'python3';                     // resolved via PATH — portable across OSes
}

// 2) Where the .h5 files live, relative to this PHP file.
$VERSION_PATH = './hdf5/version3/';

// 3) Where VCFs are written (must be web-served AND writable). Matches CFG.vcfDir in data.js.
$VCF_DIR = './vcf/';

/* ---------------------------------------------------------------------
 *  INPUT
 * ------------------------------------------------------------------- */
if ($_SERVER['REQUEST_METHOD'] !== 'POST') {
    echo json_encode(array('status' => 'error', 'message' => 'POST required'));
    exit;
}

$start         = isset($_POST['start'])     ? $_POST['start']     : '';
$end           = isset($_POST['end'])       ? $_POST['end']       : '';
$chr           = isset($_POST['chr'])       ? $_POST['chr']       : '';
$dataset       = isset($_POST['dataSet'])   ? $_POST['dataSet']   : '';
$genotypesJson = isset($_POST['genotypes']) ? $_POST['genotypes'] : '[]';
$outName       = isset($_POST['outName'])   ? $_POST['outName']   : '';

// numeric interval
if (!is_numeric($start) || !is_numeric($end)) {
    echo json_encode(array('status' => 'error', 'message' => 'Invalid interval (start/end must be numeric).'));
    exit;
}
// chromosome token must look like chr10 (used verbatim in the .h5 filename)
if (!preg_match('/^chr[0-9]{1,2}$/', $chr)) {
    echo json_encode(array('status' => 'error', 'message' => "Invalid chromosome '$chr'."));
    exit;
}

/* ---------------------------------------------------------------------
 *  DATASET  ->  (family, quality)  ->  <family>_<chr>_<quality>.h5
 * ------------------------------------------------------------------- */
$ds_part0 = 'maizegdb2026';   // family
$ds_part2 = 'HQ';             // quality tier

switch ($dataset) {
    case 'mgdb2026_hq':  $ds_part0 = 'maizegdb2026'; $ds_part2 = 'HQ';     break;
    case 'mgdb2026_hc':  $ds_part0 = 'maizegdb2026'; $ds_part2 = 'HC';     break;
    case 'mgdb2024_hq':  $ds_part0 = 'maizegdb2024'; $ds_part2 = 'HQ';     break;
    case 'mgdb2024_hc':  $ds_part0 = 'maizegdb2024'; $ds_part2 = 'HC';     break;
    case 'schnable2023': $ds_part0 = 'schnable2023'; $ds_part2 = 'impute'; break;
    case 'nam2021':      // new UI sends the bare id
    case 'nam2021_hq':   $ds_part0 = 'nam2021';      $ds_part2 = 'HQ';     break;
    case 'nam2021_hc':   $ds_part0 = 'nam2021';      $ds_part2 = 'HC';     break;
    default:
        echo json_encode(array('status' => 'error', 'message' => "Unknown dataset '$dataset'."));
        exit;
}

$db_filename = $VERSION_PATH . $ds_part0 . '_' . $chr . '_' . $ds_part2 . '.h5';
if (!is_file($db_filename)) {
    echo json_encode(array('status' => 'error',
        'message' => 'HDF5 file not found: ' . $db_filename));
    exit;
}

/* ---------------------------------------------------------------------
 *  OUTPUT PATH — force it inside $VCF_DIR (no path traversal)
 * ------------------------------------------------------------------- */
if (!is_dir($VCF_DIR)) { @mkdir($VCF_DIR, 0775, true); }
if (!is_writable($VCF_DIR)) {
    echo json_encode(array('status' => 'error',
        'message' => 'VCF directory is not writable: ' . $VCF_DIR));
    exit;
}
/* Output is gzip-compressed (.vcf.gz) — genotype text compresses ~15-20x.
 * h5_to_vcf.py writes gzip transparently when the path ends in .gz, and
 * data.js inflates the response with DecompressionStream. */
$base = basename($outName ? $outName : ('snpv_' . time() . '_' . mt_rand() . '.vcf.gz'));
if (substr($base, -7) !== '.vcf.gz') {
    $base = preg_replace('/\.vcf(\.gz)?$/', '', $base) . '.vcf.gz';
}
$vcf_path = rtrim($VCF_DIR, '/') . '/' . $base;

/* ---------------------------------------------------------------------
 *  GENOTYPES — pass the selected accession IDs through as a JSON array
 * ------------------------------------------------------------------- */
$genotypesArray = json_decode($genotypesJson);
if (!is_array($genotypesArray)) { $genotypesArray = array(); }

/* The accession list is handed to Python through a sidecar JSON file, not
 * as a command-line argument. Windows' escapeshellarg() strips '"' (and
 * '%', '!') from arguments, which shreds a JSON array on the command line;
 * a file sidesteps shell quoting entirely and behaves the same on every
 * OS. h5_to_vcf.py's 5th arg is now that file's path. */
$acc_path = $vcf_path . '.acc.json';
file_put_contents($acc_path, json_encode(array_values($genotypesArray)));

/* ---------------------------------------------------------------------
 *  RUN  h5_to_vcf.py  <db> <out.vcf> <start> <end> <accessions.json>
 * ------------------------------------------------------------------- */
$command = escapeshellarg($PYTHON_PATH) . ' ' . escapeshellarg('h5_to_vcf.py') . ' '
         . escapeshellarg($db_filename) . ' '
         . escapeshellarg($vcf_path)    . ' '
         . escapeshellarg($start)       . ' '
         . escapeshellarg($end)         . ' '
         . escapeshellarg($acc_path)    . ' 2>&1';

$output = shell_exec($command);

@unlink($acc_path);

// The script writes the VCF as a side effect; success = the file now exists.
if (is_file($vcf_path)) {
    $variants = null;
    if ($output !== null && preg_match('/^variants:\s*(\d+)/m', $output, $mm)) {
        $variants = (int) $mm[1];
    }
    echo json_encode(array(
        'status'   => 'success',
        'outFile'  => $vcf_path,
        'variants' => $variants,   // exact site count — lets the client size the result before parsing
        'message'  => 'VCF written',
        'output'   => $output,
    ));
} else if ($output !== null && strpos($output, 'No data found in the specified position range') !== false) {
    // Python ran fine, the interval simply contained no variants.
    echo json_encode(array(
        'status'  => 'empty',
        'message' => 'No variants in this interval.',
        'output'  => $output,
    ));
} else {
    // Real failure (bad Python path, missing h5py/numpy, dataset key error, ...).
    echo json_encode(array(
        'status'  => 'error',
        'message' => 'No VCF produced (script error). See output.',
        'command' => $command,
        'output'  => ($output === null ? '(no output — check that $PYTHON_PATH is correct and executable)' : $output),
    ));
}
?>
