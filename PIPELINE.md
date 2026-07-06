# XMM Automated Spectral Fitting Pipeline

End-to-end pipeline for fitting X-ray spectra from the 5XMM catalogue using
BXA (Bayesian X-ray Analysis) or standard XSPEC fitting. Designed to run on
a multi-core server with HEASOFT/PyXSPEC installed.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Prerequisites](#2-prerequisites)
3. [Input Data](#3-input-data)
4. [Quick Start](#4-quick-start)
5. [Full Procedure](#5-full-procedure)
6. [CLI Reference](#6-cli-reference)
7. [Pipeline Flow](#7-pipeline-flow)
8. [Output](#8-output)
9. [Monitoring](#9-monitoring)
10. [Resuming After Interruption](#10-resuming-after-interruption)
11. [Performance and Timing](#11-performance-and-timing)
12. [Error Codes](#12-error-codes)
13. [Spectral Models](#13-spectral-models)
14. [Troubleshooting](#14-troubleshooting)

---

## 1. Overview

The pipeline processes a list of source IDs (SRCIDs) from the 5XMM stacked
catalogue through the following stages:

```
Catalog lookup  →  List spectra  →  Validate  →  Merge  →  Background check
    →  Fit (BXA or XSPEC)  →  Export results
```

For each source, the pipeline:
- Looks up associated OBS_IDs and SRC_NUMs in the stacked catalogue
- Finds extracted spectra on disk (SRSPEC files)
- Validates spectra (checks counts, background, response files)
- Merges multiple exposures per instrument (PN, MOS1, MOS2)
- Optionally checks background quality (chi-squared p-value)
- Fits using BXA (nested sampling) or Levenberg-Marquardt
- Exports summary statistics (median, 16th/84th percentiles) to FITS

---

## 2. Prerequisites

### Software

| Package    | Purpose                          | Version tested |
|------------|----------------------------------|----------------|
| HEASOFT    | XSPEC/PyXSPEC spectral fitting   | 6.35.2         |
| SAS        | XMM-Newton Science Analysis (optional for shared libs) | DR11 |
| Python 3   | Pipeline orchestration           | 3.8+           |
| BXA        | Bayesian X-ray Analysis          | `pip install bxa` |
| UltraNest  | Nested sampling engine           | `pip install ultranest` |
| astropy    | FITS I/O                         | Any recent     |
| scipy      | Chi-squared p-value computation  | Any recent     |
| matplotlib | Corner plots                     | Any recent     |
| corner     | Posterior visualization          | `pip install corner` |

### Hardware

- Multi-core server recommended (30+ cores for production runs)
- RAM: ~2 GB per worker (BXA + XSPEC)
- Disk: ~5-10 GB for results (with `--cleanup_chains`); ~200+ GB without

---

## 3. Input Data

### Spectral Data Directory

```
{DATA_DIR}/
├── {OBS_ID_1}/
│   └── product/          # or 'pps' for 4XMM
│       ├── SRSPEC0001...  # Source spectra (hex SRC_NUM)
│       ├── BGSPEC0001...  # Background spectra
│       └── *.arf          # Ancillary response files
├── {OBS_ID_2}/
│   └── product/
│       └── ...
└── ...
```

Spectrum filenames follow the pattern `SRSPEC{hex(SRC_NUM)}` (e.g.,
`SRSPEC0067` for SRC_NUM=103).

Response matrices `*.rmf` in `RESPONSES_DIR`

### Stacked Catalogue

A FITS file cross-matching SRCIDs to their constituent OBS_ID/SRC_NUM
pairs. Expected columns include SRCID, OBS_ID, SRC_NUM.

### SRCID List

A text file with one SRCID per line, listing all sources to process:

```
3000011010100001
3000011010100002
3000011010100003
...
```

---

## 4. Quick Start

### Test run (10 sources, 2 workers)

```bash
bash run_pipeline.sh --test
```

### Full production run (30 workers, disk-safe)

```bash
bash run_pipeline.sh --cleanup_chains \
    > ~/pipeline_run.log 2>&1 &
```

### Monitor progress

```bash
bash monitor_pipeline.sh /path/to/output_dir
# or auto-refresh every 5 minutes:
watch -n 300 bash monitor_pipeline.sh /path/to/output_dir
```

---

## 5. Full Procedure

### Step 1: Configure paths

Edit the configuration section at the top of `run_pipeline.sh`:

```bash
REPO_DIR="/home/mcoriat/Work/XMM/5XMM/automated_fits"
DATA_DIR="/mnt/xmmcat/5XMM_data/Spectra"
RESPONSES_DIR="/home/mcoriat/Work/XMM/5XMM/RESPONSES"
CATALOG="${REPO_DIR}/5xmm_matched_for_pipeline.fits"
OUTPUT_DIR="/home/mcoriat/Work/XMM/5XMM/pipeline_output"
SRCID_FILE="${REPO_DIR}/srcids.txt"
NWORKERS=30
MODEL="powerlaw"
SUBDIR="product"
```

### Step 2: Verify environment

The pipeline auto-initializes HEASOFT and SAS, then checks that `xspec`,
`bxa`, and `ultranest` can be imported. If any dependency is missing, it
exits with a clear error message.

### Step 3: Find testable sources (optional)

If you want to verify the pipeline with sources that actually have spectra
on disk:

```bash
python3 find_testable_srcids.py \
    /mnt/xmmcat/5XMM_data/Spectra \
    5xmm_matched_for_pipeline.fits \
    --subdir product \
    --n_test 20 \
    --output test_srcids.txt
```

Then test with:

```bash
bash run_pipeline.sh --test --srcid_file test_srcids.txt
```

### Step 4: Launch production run

```bash
bash run_pipeline.sh --cleanup_chains \
    > ~/pipeline_run.log 2>&1 &
```

Key flags:
- `--cleanup_chains`: Delete chain.fits and corner.png after extracting
  statistics, to keep disk usage minimal.
- Without `--cleanup_chains`: Full posterior samples are preserved on disk
  (~1-2 MB per successful fit).

### Step 5: Monitor

```bash
bash monitor_pipeline.sh /path/to/output_dir
```

### Step 6: Collect results

When all chunks complete, the pipeline automatically merges per-chunk
FITS files into a single `fit_results_all.fits`:

```
{OUTPUT_DIR}/fit_results_all.fits
```

If the merge step didn't run (e.g., pipeline was backgrounded with Ctrl+C),
you can trigger it manually:

```python
from astropy.table import Table, vstack
import glob, os

output_dir = "/path/to/output_dir"
files = sorted(glob.glob(os.path.join(output_dir, "fit_results_chunk_*.fits")))
tables = [Table.read(f) for f in files]
merged = vstack(tables)
merged.write(os.path.join(output_dir, "fit_results_all.fits"), overwrite=True)
print(f"Merged {len(merged)} rows from {len(files)} chunks")
```

---

## 6. CLI Reference

### run_pipeline.sh

```
Usage: bash run_pipeline.sh [OPTIONS]

Options:
  --test             Run on 10 sources with 2 workers
  --resume           Skip already-completed sources
  --nworkers N       Number of parallel workers (default: 30)
  --model NAME       Spectral model (default: powerlaw)
                     Choices: powerlaw, apec_single, blackbody, bremss
  --srcid_file F     Path to srcid list file
  --output_dir D     Output directory
  --cleanup_chains   Delete chain.fits/corner.png after extracting stats
  --no_prefit        Skip the Levenberg-Marquardt pre-fit step
                     (used for the 5XMM-DR15 production run;
                     matches the D6.2/Viitanen+25 setup)
  --help             Show help
```

### automated_fits.py

```
Usage: python3 automated_fits.py SRCID DATA_DIR SCRIPT_DIR RESPONSES_DIR
                                 OUTPUT_DIR CATALOG OUTPUT [OPTIONS]

Positional arguments:
  srcid              SRCID to fit (or 0 for batch mode)
  data_dir           Path to spectral data directory
  script_dir         Path to scripts directory
  responses_dir      Path to response matrices directory
  output_dir         Path to output directory
  catalog            Stacked catalogue FITS file
  output             Output filename (legacy, use --export_filename)

Key options:
  --srcid_file F        Text file with SRCIDs (batch mode)
  --use_bxa             Use BXA fitting (default: standard XSPEC)
  --model_name NAME     Model: powerlaw|apec_single|blackbody|bremss
  --subdir NAME         Subdirectory under OBS_ID (default: pps, use
                        'product' for 5XMM)
  --redshift Z          Redshift for redshifted models (default: 0.0)
  --cleanup_chains      Delete chains after extracting statistics
  --skip_completed      Skip sources already fitted
  --skip_bkg_check      Skip background quality check
  --export_results_fits Export results to FITS table
  --export_filename F   Output FITS filename (default: fit_results.fits)
  --no_prefit           Disable the LM pre-fit (only tightens the
                        lg10Flux prior when enabled; DR15 production
                        ran with --no_prefit)
  --no_iin              Disable the inter-instrument normalisation
                        constant (default: constant frozen at 1 for
                        pn, free [0.5, 2.0] for MOS)
  --iin_lo X            Lower bound of the free IIN constant (0.5)
  --iin_hi X            Upper bound of the free IIN constant (2.0)
  --intrinsic_flux      Measure intrinsic flux (phabs*cflux*model,
                        5XMM legacy). Default measures observed flux
                        (cflux*phabs*model, XMM2ATHENA convention)
  --use_galabs N        Fix foreground galactic absorption (0 or 1)
  --overwrite N         Overwrite existing model (default: 1)
```

### monitor_pipeline.sh

```
Usage: bash monitor_pipeline.sh [OUTPUT_DIR]

Default: /home/mcoriat/Work/XMM/5XMM/pipeline_output

Auto-refresh:
  watch -n 300 bash monitor_pipeline.sh /path/to/output_dir
```

### find_testable_srcids.py

```
Usage: python3 find_testable_srcids.py DATA_DIR CATALOG [OPTIONS]

Options:
  --subdir NAME     Subdirectory name (default: product)
  --n_test N        Number of SRCIDs to output (default: 20)
  --output FILE     Output filename (default: test_srcids.txt)
  --max_obsids N    Max OBS_IDs to scan (default: 500)
```

---

## 7. Pipeline Flow

### Per-source processing (process_one_source)

```
 1. Create source directory: {output_dir}/{srcid}/

 2. Skip check (if --skip_completed):
    a. Check for chain.fits on disk
    b. Check for SRCID in fit_results*.fits tables
    → If found: return (0, None), skip source

 3. Catalog lookup → OBS_ID / SRC_NUM pairs
    → ERROR 1 if SRCID not in catalogue
    → ERROR 2 if no OBS_ID/SRC_NUM pairs

 4. List spectra on disk (SRSPEC file matching)
    → ERROR 3 if no spectra found

 5. Validate spectra (check_spectra):
    - Open spectrum and background files
    - Verify response (RMF) and ancillary (ARF) files exist
    - Compute source/background/net counts, SNR
    - Returns 8-element tuples:
      (file, counts, bg_counts, net_counts, exposure, flag, snr, instrument)
    - Filter: keep only spectra with flag == 0
    → ERROR 4 if no valid spectra

 6. Merge spectra per instrument (merge_spectra):
    - Two instrument groups: pn and MOS (MOS1+MOS2 share one group)
    - Outputs at most ONE spectrum per group — true merging is not
      implemented; the highest-SNR spectrum of each group is selected
    - So a fit uses at most 2 spectra (pn + MOS), i.e. at most one
      free inter-instrument constant
    - Returns 9-element tuples: (above 7 fields + sp_dic)
    - sp_dic keys: SPECFILE, BACKFILE, RESPFILE, ANCRFILE
    → ERROR 5 if merge fails

 7. Background check (unless --skip_bkg_check):
    - Quick Levenberg-Marquardt fit to each spectrum
    - Compute chi-squared p-value
    - Accept if p >= 0.01, reject otherwise
    → ERROR 3 if no spectra pass background check

 8. Spectral fitting:
    BXA mode (--use_bxa):
      a. Build model and priors (get_model_and_priors)
      b. Optional LM pre-fit; only the lg10Flux prior is
         tightened (skipped with --no_prefit, as in the
         DR15 production run)
      c. Run BXA/UltraNest nested sampling
      d. Read chain.fits → compute median, p16, p84
      e. Generate corner plot
      → ERROR 6 if fit fails

    XSPEC mode (no --use_bxa):
      a. Standard Levenberg-Marquardt fitting
      → Returns (0, None)

 9. Export and cleanup:
    - Accumulate results in memory
    - Flush to FITS every 50 successful fits (crash safety)
    - Delete chain directory if --cleanup_chains
    → Return (0, results_dict)
```

### Parallel orchestration (run_pipeline.sh)

```
 1. Initialize environment (HEASOFT, SAS, check dependencies)
 2. Split SRCID file into NWORKERS chunks
 3. Record start time → .pipeline_start_time
 4. Launch NWORKERS parallel automated_fits.py processes
 5. Wait for all jobs to complete
 6. Merge fit_results_chunk_*.fits → fit_results_all.fits
```

---

## 8. Output

### Directory structure

```
{OUTPUT_DIR}/
├── .pipeline_start_time             # Unix timestamp (for monitor)
├── chunks/
│   ├── chunk_000.txt                # SRCID lists per worker
│   ├── chunk_001.txt
│   └── ...
├── chunk_logs/
│   ├── chunk_000.log                # Per-worker stdout/stderr
│   ├── chunk_001.log
│   └── ...
├── fit_results_chunk_000.fits       # Per-chunk results (updated every 50 fits)
├── fit_results_chunk_001.fits
├── ...
├── fit_results_all.fits             # Merged final results
├── batch_processing.log
└── {SRCID}/                         # Per-source directory
    ├── {SRCID}_process_log_{model}.txt
    └── {model}_{DDMMYYYY_HHMM}/     # (deleted if --cleanup_chains)
        ├── chain.fits               # Full posterior samples
        ├── corner.png               # Corner plot
        └── (ultranest run files)
```

### Results FITS table (fit_results_all.fits)

One row per successfully fitted source. Columns:

| Column              | Type  | Description                        |
|---------------------|-------|------------------------------------|
| SRCID               | int64 | Source identifier                   |
| {param}_median      | float | Posterior median                   |
| {param}_p16         | float | 16th percentile (lower 1-sigma)    |
| {param}_p84         | float | 84th percentile (upper 1-sigma)    |

Where `{param}` depends on the model (see [Spectral Models](#13-spectral-models)):
- powerlaw: `nH`, `PhoIndex`, `lg10Flux`
- apec_single: `nH`, `kT`, `lg10Flux`
- blackbody: `nH`, `kT`, `lg10Flux`
- bremss: `nH`, `kT`, `lg10Flux`

`lg10Flux` is the log10 of the unabsorbed flux in erg/cm^2/s in the
0.5-10.0 keV band (from the XSPEC `cflux` convolution model).

---

## 9. Monitoring

### monitor_pipeline.sh

Displays real-time progress of a running pipeline:

```
==============================================
 Pipeline Progress Monitor
 Output dir: /home/mcoriat/Work/XMM/5XMM/pipeline_output
 Time: Tue Feb 24 09:41:53 UTC 2026
==============================================
 Overall progress:
   Total SRCIDs:    646108
   Touched:         2526 (processed: 2526, skipped: 0)
   Successful fits: 305
   Error rate:      88.0%

 Per-chunk status:
 ─────────────────────────────────────────
   chunk_000     [RUNNING]  proc: 16    skip: 0    fits: 11   last: 300001101...
   chunk_001     [RUNNING]  proc: 54    skip: 0    fits: 11   last: 300491501...
   ...

 Running processes: 30
 Disk usage: 1023M

 Timing (elapsed: 0.9h):
   Processing rate: ~2963 sources/hour
   Fit rate:        ~358 fits/hour
   Expected fits:   ~78013 (of 646108 total)
   Est. remaining:  ~217.2 hours
==============================================
```

**Chunk statuses:**
- `RUNNING`: Process actively writing to log
- `DONE`: "Batch processing complete" found in log
- `PENDING`: Log exists but no sources processed yet
- `STOPPED?`: Log has content but process is no longer writing

**ETA calculation:** Uses the overall processing rate (all sources, including
fast failures) rather than just the fit rate. This accounts for the ~64% of
sources that lack extracted spectra and are processed in seconds (observed
fit yield in the 5XMM production runs: ~36%).

---

## 10. Resuming After Interruption

If the pipeline is interrupted (killed, crashed, server reboot):

```bash
# Resume: skip sources already in fit_results*.fits
bash run_pipeline.sh --resume --cleanup_chains \
    > ~/pipeline_run.log 2>&1 &
```

How `--resume` works:
1. At startup, each worker scans `fit_results_chunk_*.fits` in the output
   directory
2. Extracts the set of SRCIDs with existing results
3. For each source, checks:
   - Does `chain.fits` exist on disk? (works without `--cleanup_chains`)
   - Is the SRCID already in a results FITS table? (works with `--cleanup_chains`)
4. Skips matched sources immediately

**Crash safety:** Results are flushed to FITS every 50 successful fits.
At most ~50 fits worth of work is lost on any interruption.

**Killing the pipeline:**
```bash
pkill -f "automated_fits.py"       # graceful
pkill -9 -f "automated_fits.py"    # force (if needed)
```

---

## 11. Performance and Timing

### BXA fitting (--use_bxa) — observed in the 5XMM production runs

| Metric                  | Observed value             |
|-------------------------|----------------------------|
| Time per BXA fit        | typically minutes          |
| Fit yield               | ~36% of catalogue sources  |
| Disk (--cleanup_chains) | ~5-10 GB                   |
| Disk (full posteriors)  | ~200+ GB                   |
| RAM per worker          | ~1-2 GB                    |

The ~64% of sources that produce no fit mostly have no extracted
spectra on disk (ERROR 3). These fail in seconds and do not
significantly impact total runtime.

### Standard XSPEC fitting (without --use_bxa)

| Metric                  | Typical value              |
|-------------------------|----------------------------|
| Time per fit            | 5-10 seconds               |
| Speedup vs BXA          | ~20-50x                    |
| Expected total runtime  | Hours (not days)           |

### Pre-fitting (optional, off in production)

BXA has an optional pre-fit step (`_prefit_and_tighten`) that:
1. Runs a quick Levenberg-Marquardt fit with wide parameter ranges
2. Tightens the lg10Flux prior to best-fit ± 2 dex (only if LM moved it)
3. NEVER tightens nH or the shape parameter — doing so creates a
   starting-value pile-up artefact for low-S/N sources (the 5e20 cm^-2
   nH peak; see D6.2 / Viitanen+25)

Because only one parameter is tightened, the speed-up is modest.
The 5XMM-DR15 production run disabled the step entirely
(`run_pipeline.sh --no_prefit`) to match the D6.2/Viitanen+25 setup.

### Parallelism

- The SRCID list is split into N equal chunks (one per worker)
- Each worker runs as an independent Python process
- Workers share no state; parallelism is embarrassingly parallel
- Thread environment variables are set to prevent internal parallelism
  within XSPEC (`OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`)
- Recommended: `NWORKERS = total_cores - 10` (leave headroom for OS and I/O)

---

## 12. Error Codes

Returned by `process_one_source()`:

| Code | Meaning                                                     |
|------|-------------------------------------------------------------|
| 0    | Success                                                     |
| 1    | SRCID not found in catalogue                                |
| 2    | SRCID found but no OBS_ID/SRC_NUM pairs                     |
| 3    | No extracted spectra on disk, or background check failed    |
| 4    | Spectra found but none are valid (can't open, missing files)|
| 5    | Valid spectra exist but merging failed                       |
| 6    | Fit ran but failed to produce results                        |

Most sources fail with ERROR 3 (no spectra on disk). This is expected
for faint catalogue sources without pipeline-extracted spectra.

### Internal spectrum flags

Set by `check_spectra()` at tuple position [5]:

| Flag | Meaning                                        |
|------|------------------------------------------------|
| -2   | Cannot open spectral file (internal only)      |
| -1   | Cannot open background file (internal only)    |
|  0   | No issues detected                             |
|  1   | Zero or negative background counts             |
|  2   | Zero or negative source (total) or net counts  |
|  3   | Merge failed, or pre-BXA quick-fit screen failed / p < 0.01 |
|  4   | BXA solver failed (timeout, crash, no chain)   |
|  5-11| Quality warnings on a completed fit (poor GoF, PhoIndex pegged, nH ≥ 1e24 cm^-2; see automated_fits.py docstring) |

---

## 13. Spectral Models

All models use `phabs` for absorption and `cflux` for flux measurement
in the 0.5-10.0 keV band (the fit itself uses 0.3-10 keV).

By default the flux is the OBSERVED (absorbed) flux — `cflux` wraps
`phabs*emission` (XMM2ATHENA convention, Viitanen+25). With
`--intrinsic_flux` the order becomes `phabs*cflux*emission` and
lg10Flux measures the intrinsic (unabsorbed) flux (5XMM legacy).
For simultaneous multi-instrument fits a leading `constant*` is
prepended for the inter-instrument normalisation (IIN): frozen at 1
for pn, free in [0.5, 2.0] with a log-uniform prior for MOS
(disable with `--no_iin`).

Priors: nH is log-uniform (Jeffreys) over its full range; PhoIndex/kT
are uniform; lg10Flux is uniform (= log-uniform on the linear flux).

### powerlaw (`constant * cflux * phabs * zpowerlw`, Redshift frozen)

| Parameter | Range          | Description                              |
|-----------|----------------|------------------------------------------|
| nH        | 0.001 - 1000.0 | Column density (10^22 cm^-2), log-uniform prior |
| PhoIndex  | 1.0 - 3.0      | Photon index                             |
| lg10Flux  | -15.0 - -9.0   | log10 flux (erg/cm^2/s, 0.5-10 keV)      |

### apec_single (`constant * cflux * phabs * apec`)

| Parameter | Range          | Description                              |
|-----------|----------------|------------------------------------------|
| nH        | 0.001 - 1000.0 | Column density (10^22 cm^-2), log-uniform prior |
| kT        | 0.1 - 10.0     | Plasma temperature (keV)                 |
| lg10Flux  | -15.0 - -9.0   | log10 flux (erg/cm^2/s, 0.5-10 keV)      |

### blackbody (`constant * cflux * phabs * bbody`)

| Parameter | Range          | Description                              |
|-----------|----------------|------------------------------------------|
| nH        | 0.001 - 1000.0 | Column density (10^22 cm^-2), log-uniform prior |
| kT        | 0.01 - 2.0     | Temperature (keV)                        |
| lg10Flux  | -15.0 - -9.0   | log10 flux (erg/cm^2/s, 0.5-10 keV)      |

### bremss (`constant * cflux * phabs * bremss`)

| Parameter | Range          | Description                              |
|-----------|----------------|------------------------------------------|
| nH        | 0.001 - 1000.0 | Column density (10^22 cm^-2), log-uniform prior |
| kT        | 0.1 - 20.0     | Plasma temperature (keV)                 |
| lg10Flux  | -15.0 - -9.0   | log10 flux (erg/cm^2/s, 0.5-10 keV)      |

### BXA solver parameters

| Parameter            | Value                       | Notes                              |
|----------------------|-----------------------------|------------------------------------|
| n_live_points        | max(35 * n_free, 100)       | 105 for 3-parameter models         |
| evidence_tolerance   | 1.0                         | OK for parameter estimation; use 0.5 for model comparison |
| speed                | "safe"                      | BXA default step sampler           |
| Lepsilon             | 0.1                         | Likelihood tolerance               |
| pre-fit              | On by default, off with --no_prefit | Only tightens lg10Flux; DR15 production ran with --no_prefit |

---

## 14. Troubleshooting

### "No result files found to merge"

The merge step runs after all workers complete. If the pipeline was
killed before any chunk finished, no `fit_results_chunk_*.fits` files
exist. With the incremental flush (every 50 fits), partial results are
saved even if a chunk doesn't complete. Check for partial files:

```bash
ls -la /path/to/output_dir/fit_results_chunk_*.fits
```

### Pipeline runs but 0 fits produced

Check a chunk log for the actual error:

```bash
tail -100 /path/to/output_dir/chunk_logs/chunk_000.log
```

Common causes:
- `Redshift: value: None` → `--redshift` not set (fixed: defaults to 0.0)
- `'float' object is not callable` → `Fit.testStatistic()` should be
  `Fit.testStatistic` (property, not method) — fixed in current code
- All sources fail with ERROR 3 → sources in the SRCID list don't have
  spectra on disk. Use `find_testable_srcids.py` to find sources with
  spectra.

### Disk space growing too fast

Use `--cleanup_chains` to delete chain.fits and corner.png after
extracting summary statistics. This reduces per-fit disk usage from
~1-2 MB to essentially zero.

To clean up chains from a run that was started without `--cleanup_chains`:

```bash
find /path/to/output_dir -name "chain.fits" -printf '%h\n' | xargs -r rm -rf
```

### Restarting after cleanup

If you cleaned up chain files, `--resume` still works because it checks
the `fit_results_chunk_*.fits` tables (not just chain.fits on disk).

### XSPEC/HEASOFT not found

Ensure HEASOFT is properly installed and the init script path in
`run_pipeline.sh` is correct:

```bash
source /path/to/heasoft/headas-init.sh
python3 -c "import xspec; print('OK')"
```

### BXA/UltraNest import errors

```bash
pip install bxa ultranest corner
```
