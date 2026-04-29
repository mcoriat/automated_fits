#!/bin/bash
# ===========================================================
# run_validation.sh
# ===========================================================
# Runs the prefit-vs-no_prefit validation experiment on bell.
#
# This is a SCALED-DOWN version of run_pipeline.sh, designed to
# fit a small (~600 source) sample with --no_prefit so its
# results can be compared to the existing production catalogue
# (which used prefit=True).
#
# Workflow:
#   1. Select stratified sample of SRCIDs from the existing
#      production FITS catalogue.
#   2. Run automated_fits.py with --no_prefit on that sample.
#   3. Merge the output FITS chunks.
#   4. Run compare_prefit_validation.py to produce histograms
#      and summary statistics.
#
# Usage:
#   ./run_validation.sh [N_PER_BIN] [N_WORKERS]
#
#   N_PER_BIN   : sources per SNR bin (default 200)
#                 total = 3 * N_PER_BIN
#   N_WORKERS   : parallel workers (default 30)
#
# IMPORTANT: edit the PRODUCTION_CATALOG path below to point to
# your existing prefit=True catalogue before running.
# ===========================================================

set -e
set -u

# -- USER CONFIG --
PRODUCTION_CATALOG="${PRODUCTION_CATALOG:-/data/scratch/pipeline_output/fit_results_all.fits}"
DATA_DIR="${DATA_DIR:-/storage/xmmcat/5XMM_data/Spectra}"
SCRIPT_DIR="${SCRIPT_DIR:-$(pwd)}"
RESPONSES_DIR="${RESPONSES_DIR:-/home/mcoriat/XMM/RESPONSES}"
CATALOG="${CATALOG:-/home/mcoriat/XMM/5XMM/entry_products/5xmm_pps_matched_for_spec_pipe.fits}"
MODEL_NAME="${MODEL_NAME:-powerlaw}"
SUBDIR="${SUBDIR:-product}"   # 'product' or 'pps'
HEASOFT_BASE="${HEASOFT_BASE:-/home/mcoriat/Software}"

# -- ARGS --
N_PER_BIN="${1:-200}"
NWORKERS="${2:-30}"

# ===========================================================
# Initialize HEASOFT, venv, SAS (copied from run_pipeline.sh)
# ===========================================================
HEADAS_FOUND=""
for d in "${HEASOFT_BASE}"/heasoft-*/x86_64-*; do
    if [ -f "${d}/headas-init.sh" ]; then
        HEADAS_FOUND="${d}"
    fi
done
if [ -n "${HEADAS_FOUND}" ]; then
    export HEADAS="${HEADAS_FOUND}"
    echo "Initializing HEASOFT: ${HEADAS}"
    source "${HEADAS}/headas-init.sh"
else
    echo "ERROR: No HEASOFT found under ${HEASOFT_BASE}/"
    exit 1
fi

VENV_DIR="${SCRIPT_DIR}/.venv"
if [ -f "${VENV_DIR}/bin/activate" ]; then
    echo "Activating venv: ${VENV_DIR}"
    source "${VENV_DIR}/bin/activate"
else
    echo "WARNING: No venv at ${VENV_DIR}, using system Python"
fi

# SAS (optional)
if [ -z "${SAS_DIR:-}" ]; then
    for sas_candidate in \
        /home/mcoriat/Software/sas/sas-setup.sh \
        /opt/sas/sas-setup.sh; do
        if [ -f "${sas_candidate}" ]; then
            export SAS_DIR="$(dirname "${sas_candidate}")"
            set +eu; source "${sas_candidate}" 2>/dev/null; set -eu
            break
        fi
    done
fi

# Quick sanity check
echo "Checking dependencies..."
python -c "import xspec; print(f'  XSPEC:    {xspec.__file__}')" || {
    echo "ERROR: cannot import xspec"; exit 1; }
python -c "import bxa; print(f'  BXA:      {bxa.__file__}')" || {
    echo "ERROR: cannot import bxa"; exit 1; }
echo "  All OK."
echo ""

# All paths must be absolute — the export function joins
# output_dir + export_filename, and relative paths double up.
VALIDATION_DIR="$(pwd)/validation_$(date +%Y%m%d_%H%M)"
SRCID_FILE="${VALIDATION_DIR}/validation_srcids.txt"
META_FILE="${VALIDATION_DIR}/validation_srcids.meta.fits"
OUTPUT_DIR="${VALIDATION_DIR}/no_prefit_run"
PLOT_DIR="${VALIDATION_DIR}/comparison_plots"

mkdir -p "${VALIDATION_DIR}" "${OUTPUT_DIR}" "${PLOT_DIR}"

echo "==========================================================="
echo "  Validation experiment: prefit vs no_prefit"
echo "==========================================================="
echo "  Production catalogue : ${PRODUCTION_CATALOG}"
echo "  N per SNR bin        : ${N_PER_BIN}"
echo "  Workers              : ${NWORKERS}"
echo "  Output dir           : ${VALIDATION_DIR}"
echo "==========================================================="

# -- 1. Select stratified sample --
echo
echo "[1/4] Selecting stratified sample..."
python "${SCRIPT_DIR}/select_validation_sample.py" \
    --catalog "${PRODUCTION_CATALOG}" \
    --n_per_bin "${N_PER_BIN}" \
    --output "${SRCID_FILE}" \
    --metadata "${META_FILE}"
N_SRC=$(wc -l < "${SRCID_FILE}")
echo "  Selected ${N_SRC} sources"

# -- 2. Split into chunks for parallel workers --
echo
echo "[2/4] Splitting into ${NWORKERS} chunks..."
split -d -a 3 -n l/${NWORKERS} "${SRCID_FILE}" \
    "${VALIDATION_DIR}/chunk_"
ls "${VALIDATION_DIR}"/chunk_*

# -- 3. Launch workers in parallel --
echo
echo "[3/4] Launching ${NWORKERS} parallel workers (no_prefit)..."
ulimit -n 16384 || true

PIDS=()
for chunk in "${VALIDATION_DIR}"/chunk_*; do
    chunk_id=$(basename "${chunk}")
    log="${OUTPUT_DIR}/${chunk_id}.log"
    out_fits="${OUTPUT_DIR}/${chunk_id}.fits"
    (
        python "${SCRIPT_DIR}/automated_fits.py" 0 \
            "${DATA_DIR}" "${SCRIPT_DIR}" "${RESPONSES_DIR}" \
            "${OUTPUT_DIR}/${chunk_id}_results" \
            "${CATALOG}" "${out_fits}" \
            --srcid_file "${chunk}" \
            --use_bxa \
            --model_name "${MODEL_NAME}" \
            --subdir "${SUBDIR}" \
            --no_prefit \
            --export_results_fits \
            --export_filename "${out_fits}" \
            --cleanup_chains \
            > "${log}" 2>&1
    ) &
    PIDS+=($!)
    echo "  worker ${chunk_id}: PID $!"
done

echo
echo "Waiting for ${#PIDS[@]} workers (this will take HOURS — "
echo "no_prefit is slow on purpose). Monitor with:"
echo "  tail -f ${OUTPUT_DIR}/chunk_*.log"
wait "${PIDS[@]}"
echo "  All workers done."

# -- 4. Merge per-chunk FITS into single validation catalogue --
echo
echo "[4/4] Merging chunk FITS files..."
python <<EOF
from astropy.table import Table, vstack
import glob, os
fits_files = sorted(glob.glob("${OUTPUT_DIR}/chunk_*.fits"))
tables = [Table.read(f) for f in fits_files if os.path.getsize(f) > 0]
if tables:
    merged = vstack(tables)
    merged.write("${VALIDATION_DIR}/no_prefit_merged.fits",
                 overwrite=True)
    print(f"Merged {len(tables)} chunks, "
          f"{len(merged)} total sources -> "
          "${VALIDATION_DIR}/no_prefit_merged.fits")
else:
    print("WARNING: no chunk FITS files found — check logs")
EOF

# -- 5. Run comparison --
echo
echo "[5/5] Running comparison vs prefit catalogue..."
python "${SCRIPT_DIR}/compare_prefit_validation.py" \
    --prefit "${PRODUCTION_CATALOG}" \
    --no_prefit "${VALIDATION_DIR}/no_prefit_merged.fits" \
    --metadata "${META_FILE}" \
    --output_dir "${PLOT_DIR}"

echo
echo "==========================================================="
echo "  Validation complete."
echo "  Plots:   ${PLOT_DIR}/*.png"
echo "  Summary: ${PLOT_DIR}/summary.json"
echo "==========================================================="
