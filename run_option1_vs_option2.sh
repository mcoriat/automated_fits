#!/bin/bash
# ===========================================================
# run_option1_vs_option2.sh
# ===========================================================
# Test driver: compares two BXA fitting setups on the SAME
# stratified sample to decide whether to keep the pre-fit step.
#
#   Option 1: --no_prefit   (BXA only, paper-aligned)
#   Option 2: prefit ON     (current code: log-uniform nH,
#                            no nH/PhoIndex tightening; only
#                            lg10Flux is tightened from LM)
#
# The script:
#   1. Selects a stratified SRCID sample (via
#      select_validation_sample.py).
#   2. Splits it into N chunks.
#   3. Runs Option 1 across all chunks; records wall-time.
#   4. Runs Option 2 on the SAME chunks; records wall-time.
#   5. Merges per-chunk FITS for each option.
#   6. Runs compare_prefit_validation.py on
#      (Option 2 = "prefit" arg) vs (Option 1 = "no_prefit" arg).
#   7. Prints a timing summary.
#
# Usage:
#   ./run_option1_vs_option2.sh [N_PER_BIN] [N_WORKERS]
#       N_PER_BIN  : sources per SNR bin (default 200 → 600 total)
#       N_WORKERS  : parallel workers per option run (default 30)
#
# Notes:
# - The two options run SEQUENTIALLY (not concurrently) so
#   each one gets the full machine and the timing is not
#   contaminated by contention.
# - The sample is selected ONCE; both options process the same
#   SRCIDs in the same order, so the comparison is paired.
# - Edit the USER CONFIG block below before first run.
# ===========================================================

set -e
set -u

# ---- USER CONFIG ----
PRODUCTION_CATALOG="${PRODUCTION_CATALOG:-/data/scratch/pipeline_output/fit_results_all.fits}"
DATA_DIR="${DATA_DIR:-/storage/xmmcat/5XMM_data/Spectra}"
SCRIPT_DIR="${SCRIPT_DIR:-$(pwd)}"
RESPONSES_DIR="${RESPONSES_DIR:-/home/mcoriat/XMM/RESPONSES}"
CATALOG="${CATALOG:-/home/mcoriat/XMM/5XMM/entry_products/5xmm_pps_matched_for_spec_pipe.fits}"
MODEL_NAME="${MODEL_NAME:-powerlaw}"
SUBDIR="${SUBDIR:-product}"
HEASOFT_BASE="${HEASOFT_BASE:-/home/mcoriat/Software}"

# ---- ARGS ----
N_PER_BIN="${1:-200}"
NWORKERS="${2:-30}"

# ===========================================================
# Initialize HEASOFT, venv, SAS (copied verbatim from
# run_validation.sh so behaviour is identical)
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

if [ -z "${SAS_DIR:-}" ]; then
    for sas_candidate in \
        /home/mcoriat/Software/sas/sas-setup.sh \
        /opt/sas/sas-setup.sh; do
        if [ -f "${sas_candidate}" ]; then
            export SAS_DIR="$(dirname "${sas_candidate}")"
            set +eu
            source "${sas_candidate}" 2>/dev/null
            set -eu
            break
        fi
    done
fi

echo "Checking dependencies..."
python -c "import xspec; print(f'  XSPEC: {xspec.__file__}')" \
    || { echo "ERROR: cannot import xspec"; exit 1; }
python -c "import bxa; print(f'  BXA:   {bxa.__file__}')" \
    || { echo "ERROR: cannot import bxa"; exit 1; }
echo "  All OK."
echo ""

# ===========================================================
# Output layout
# ===========================================================
EXP_DIR="$(pwd)/exp_opt1_vs_opt2_$(date +%Y%m%d_%H%M)"
SRCID_FILE="${EXP_DIR}/srcids.txt"
META_FILE="${EXP_DIR}/srcids.meta.fits"
OPT1_DIR="${EXP_DIR}/option1_no_prefit"
OPT2_DIR="${EXP_DIR}/option2_prefit"
PLOT_DIR="${EXP_DIR}/comparison_plots"
TIME_LOG="${EXP_DIR}/timings.txt"

mkdir -p "${EXP_DIR}" "${OPT1_DIR}" "${OPT2_DIR}" "${PLOT_DIR}"

echo "==========================================================="
echo "  Test: Option 1 (no prefit) vs Option 2 (prefit)"
echo "==========================================================="
echo "  Production catalogue : ${PRODUCTION_CATALOG}"
echo "  N per SNR bin        : ${N_PER_BIN}  (total ~3*N)"
echo "  Workers per option   : ${NWORKERS}"
echo "  Output dir           : ${EXP_DIR}"
echo "==========================================================="

# ===========================================================
# 1. Select stratified sample (once — same for both options)
# ===========================================================
echo
echo "[1/6] Selecting stratified sample..."
python "${SCRIPT_DIR}/select_validation_sample.py" \
    --catalog "${PRODUCTION_CATALOG}" \
    --n_per_bin "${N_PER_BIN}" \
    --output "${SRCID_FILE}" \
    --metadata "${META_FILE}"
N_SRC=$(wc -l < "${SRCID_FILE}")
echo "  Selected ${N_SRC} sources"

# ===========================================================
# 2. Split into NWORKERS chunks (once — same for both options)
# ===========================================================
echo
echo "[2/6] Splitting into ${NWORKERS} chunks..."
split -d -a 3 -n l/${NWORKERS} "${SRCID_FILE}" \
    "${EXP_DIR}/chunk_"
ls "${EXP_DIR}"/chunk_* | head -3
echo "  ..."

# ===========================================================
# Helper: run all chunks for one option, time it
# ===========================================================
run_one_option() {
    local LABEL="$1"
    local OUTDIR="$2"
    local EXTRA_FLAG="$3"

    echo
    echo "==========================================================="
    echo "  Running ${LABEL} (${NWORKERS} workers)"
    echo "==========================================================="
    ulimit -n 16384 || true

    local START=$(date +%s)
    local PIDS=()
    for chunk in "${EXP_DIR}"/chunk_*; do
        local chunk_id=$(basename "${chunk}")
        local log="${OUTDIR}/${chunk_id}.log"
        local out_fits="${OUTDIR}/${chunk_id}.fits"
        (
            python "${SCRIPT_DIR}/automated_fits.py" 0 \
                "${DATA_DIR}" "${SCRIPT_DIR}" "${RESPONSES_DIR}" \
                "${OUTDIR}/${chunk_id}_results" \
                "${CATALOG}" "${out_fits}" \
                --srcid_file "${chunk}" \
                --use_bxa \
                --model_name "${MODEL_NAME}" \
                --subdir "${SUBDIR}" \
                ${EXTRA_FLAG} \
                --export_results_fits \
                --export_filename "${out_fits}" \
                --cleanup_chains \
                > "${log}" 2>&1
        ) &
        PIDS+=($!)
    done
    echo "  Workers launched: ${#PIDS[@]}"
    echo "  Monitor: tail -f ${OUTDIR}/chunk_*.log"
    wait "${PIDS[@]}"
    local END=$(date +%s)
    local ELAPSED=$((END - START))
    local PER_SRC=$(awk -v e=${ELAPSED} -v n=${N_SRC} \
        -v w=${NWORKERS} 'BEGIN { printf "%.1f", e*w/n }')

    echo "${LABEL}: total wall = ${ELAPSED}s, per-source = ${PER_SRC}s (×${NWORKERS} workers)" \
        | tee -a "${TIME_LOG}"
}

# ===========================================================
# 3. Run Option 1 (no prefit) first — fully parallel
# ===========================================================
echo
echo "[3/6] Running Option 1 (--no_prefit)..."
run_one_option "Option1_no_prefit" "${OPT1_DIR}" "--no_prefit"

# ===========================================================
# 4. Run Option 2 (prefit on, current patched code)
# ===========================================================
echo
echo "[4/6] Running Option 2 (prefit ON, log-uniform nH, only lg10Flux tightened)..."
run_one_option "Option2_prefit" "${OPT2_DIR}" ""

# ===========================================================
# 5. Merge per-chunk FITS for each option
# ===========================================================
echo
echo "[5/6] Merging per-chunk FITS..."
for d in "${OPT1_DIR}" "${OPT2_DIR}"; do
    label=$(basename "${d}")
    python <<EOF
from astropy.table import Table, vstack
import glob, os
fits_files = sorted(glob.glob("${d}/chunk_*.fits"))
tables = [Table.read(f) for f in fits_files if os.path.getsize(f) > 0]
if tables:
    merged = vstack(tables)
    out = "${d}/merged.fits"
    merged.write(out, overwrite=True)
    print(f"  ${label}: merged {len(tables)} chunks, "
          f"{len(merged)} sources -> {out}")
else:
    print(f"  WARNING (${label}): no chunk FITS files found")
EOF
done

# ===========================================================
# 6. Run comparison (reuses existing script)
# ===========================================================
echo
echo "[6/6] Running comparison..."
python "${SCRIPT_DIR}/compare_prefit_validation.py" \
    --prefit "${OPT2_DIR}/merged.fits" \
    --no_prefit "${OPT1_DIR}/merged.fits" \
    --metadata "${META_FILE}" \
    --output_dir "${PLOT_DIR}" \
    || echo "  WARNING: compare script returned non-zero — check output_dir"

# ===========================================================
# Summary
# ===========================================================
echo
echo "==========================================================="
echo "  Summary"
echo "==========================================================="
echo
echo "  Sample size  : ${N_SRC} sources"
echo "  Workers      : ${NWORKERS}"
echo
echo "  Wall-times:"
cat "${TIME_LOG}"
echo
echo "  Outputs:"
echo "    Option 1 catalog : ${OPT1_DIR}/merged.fits"
echo "    Option 2 catalog : ${OPT2_DIR}/merged.fits"
echo "    Comparison plots : ${PLOT_DIR}/"
echo "    Timing log       : ${TIME_LOG}"
echo
echo "  Decision rule:"
echo "    - If wall-times are within ~30%, go with Option 1"
echo "      (simpler, paper-aligned)."
echo "    - If Option 2 is >2× faster AND posteriors agree,"
echo "      go with Option 2 (keep the pre-fit for speed)."
echo "    - If posteriors DISAGREE significantly, investigate"
echo "      the residual lg10Flux tightening in Option 2."
echo "==========================================================="
