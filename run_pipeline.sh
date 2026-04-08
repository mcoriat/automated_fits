#!/bin/bash
# ==============================================================
#  run_pipeline.sh — Launch the automated spectral fitting
#  pipeline in parallel across multiple cores.
#
#  Usage:
#    # Test run (10 sources, 2 workers):
#    bash run_pipeline.sh --test
#
#    # Full run (all sources, 60 workers):
#    bash run_pipeline.sh
#
#    # Restart after interruption (skips completed sources):
#    bash run_pipeline.sh --resume
#
#    # Custom settings:
#    bash run_pipeline.sh --nworkers 20 --model powerlaw
#
# ==============================================================

set -euo pipefail

# Raise file descriptor limit to prevent "Too many open files"
# during long BXA fitting runs (XSPEC + ultranest + matplotlib).
# The outer fd fence in automated_fits.py should keep the
# steady-state count low, but a generous limit avoids any
# edge-case exhaustion before the fence runs.
ulimit -n 16384 2>/dev/null || ulimit -n 4096 2>/dev/null || true

# ==============================================================
# CONFIGURATION — edit these paths for your setup
# ==============================================================
REPO_DIR="/home/mcoriat/XMM/5XMM/automated_fits"
DATA_DIR="/storage/xmmcat/5XMM_data/Spectra"
RESPONSES_DIR="/home/mcoriat/XMM/RESPONSES"
CATALOG="${REPO_DIR}/5xmm_matched_for_pipeline.fits"
OUTPUT_DIR="/data/scratch/pipeline_output"
SRCID_FILE="${REPO_DIR}/srcids.txt"

# Default parallel settings
NWORKERS=60          # Number of parallel batch jobs (Lovelace: 128 threads, keep headroom)
MODEL="powerlaw"     # Spectral model: powerlaw, apec_single, blackbody, bremss
SUBDIR="product"     # Subdirectory under OBS_ID: 'product' (5XMM) or 'pps' (4XMM)
TEST_NSOURCES=10     # Number of sources for test run
CLEANUP_CHAINS=false # Delete chain.fits/corner.png after extracting statistics

# ==============================================================
# PARSE COMMAND-LINE OPTIONS
# ==============================================================
DO_TEST=false
DO_RESUME=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --test)
            DO_TEST=true
            shift
            ;;
        --resume)
            DO_RESUME=true
            shift
            ;;
        --nworkers)
            NWORKERS="$2"
            shift 2
            ;;
        --model)
            MODEL="$2"
            shift 2
            ;;
        --srcid_file)
            SRCID_FILE="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --cleanup_chains)
            CLEANUP_CHAINS=true
            shift
            ;;
        --help|-h)
            echo "Usage: bash run_pipeline.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --test           Run on ${TEST_NSOURCES} sources with 2 workers"
            echo "  --resume         Skip already-completed sources"
            echo "  --nworkers N     Number of parallel workers (default: ${NWORKERS})"
            echo "  --model NAME     Spectral model (default: ${MODEL})"
            echo "  --srcid_file F   Path to srcid list (default: ${SRCID_FILE})"
            echo "  --output_dir D   Output directory (default: ${OUTPUT_DIR})"
            echo "  --cleanup_chains Delete chain.fits/corner.png after extracting stats"
            echo "  --help           Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ==============================================================
# ENVIRONMENT SETUP
# ==============================================================
echo "=============================================="
echo " XMM Automated Spectral Fitting Pipeline"
echo "=============================================="
echo ""
echo " Started: $(date)"
echo ""

# Initialize HEASOFT (required for PyXSPEC)
if [ -z "${HEADAS:-}" ]; then
    # Auto-detect HEASOFT installation
    HEASOFT_BASE="/home/mcoriat/Software"
    HEADAS_FOUND=""
    for d in "${HEASOFT_BASE}"/heasoft-*/x86_64-*; do
        if [ -f "${d}/headas-init.sh" ]; then
            HEADAS_FOUND="${d}"
        fi
    done

    if [ -n "${HEADAS_FOUND}" ]; then
        export HEADAS="${HEADAS_FOUND}"
        echo " Initializing HEASOFT..."
        echo "   ${HEADAS}"
        source "${HEADAS}/headas-init.sh"
    else
        echo " ERROR: No HEASOFT installation found under ${HEASOFT_BASE}/"
        echo "        Expected: ${HEASOFT_BASE}/heasoft-*/x86_64-*/headas-init.sh"
        exit 1
    fi
else
    echo " HEASOFT already loaded: ${HEADAS}"
fi

# Activate Python venv (contains bxa, ultranest, astropy, etc.)
# Must be activated AFTER HEASOFT so that PYTHONPATH from
# headas-init.sh is preserved inside the venv.
VENV_DIR="${REPO_DIR}/.venv"
if [ -f "${VENV_DIR}/bin/activate" ]; then
    echo " Activating Python venv: ${VENV_DIR}"
    source "${VENV_DIR}/bin/activate"
else
    echo " WARNING: No venv found at ${VENV_DIR}"
    echo "          Assuming system Python has all dependencies."
fi

# Initialize SAS if available (not required for BXA fits, but
# some shared libraries may depend on it).
if [ -z "${SAS_DIR:-}" ]; then
    # Try common SAS locations — skip silently if not installed
    for sas_candidate in \
        /home/mcoriat/Software/sas/sas-setup.sh \
        /opt/sas/sas-setup.sh; do
        if [ -f "${sas_candidate}" ]; then
            export SAS_DIR="$(dirname "${sas_candidate}")"
            echo " Initializing SAS: ${SAS_DIR}"
            set +eu
            source "${sas_candidate}" 2>/dev/null
            set -eu
            break
        fi
    done
    if [ -z "${SAS_DIR:-}" ]; then
        echo " SAS not found (not required for BXA fits)."
    fi
else
    echo " SAS already loaded: ${SAS_DIR}"
fi

# Verify critical dependencies
echo ""
echo " Checking dependencies..."
python3 -c "import xspec; print(f'   XSPEC: {xspec.__file__}')" 2>/dev/null || {
    echo "   ERROR: Cannot import xspec. Is HEASOFT+PyXSPEC set up?"
    exit 1
}
python3 -c "import bxa; print(f'   BXA:   {bxa.__file__}')" 2>/dev/null || {
    echo "   ERROR: Cannot import bxa. Is the venv activated?"
    echo "          Expected venv at: ${VENV_DIR}"
    exit 1
}
python3 -c "import ultranest; print(f'   UltraNest: {ultranest.__file__}')" 2>/dev/null || {
    echo "   ERROR: Cannot import ultranest."
    exit 1
}
echo "   All OK."

# ==============================================================
# VERIFY PATHS
# ==============================================================
echo ""
echo " Verifying paths..."
for label_path in \
    "Repo:${REPO_DIR}" \
    "Data:${DATA_DIR}" \
    "Responses:${RESPONSES_DIR}" \
    "Catalog:${CATALOG}" \
    "SRCID file:${SRCID_FILE}"; do
    label="${label_path%%:*}"
    path="${label_path#*:}"
    if [ -e "${path}" ]; then
        echo "   ${label}: OK"
    else
        echo "   ${label}: MISSING — ${path}"
        echo "   ERROR: Fix this path before running."
        exit 1
    fi
done

# Create output directory
mkdir -p "${OUTPUT_DIR}"
echo "   Output dir: ${OUTPUT_DIR}"

# ==============================================================
# PREPARE SRCID CHUNKS
# ==============================================================
TOTAL_SRCIDS=$(wc -l < "${SRCID_FILE}" | tr -d ' ')
echo ""
echo " Total SRCIDs in file: ${TOTAL_SRCIDS}"

# For test mode, use a small subset
if [ "${DO_TEST}" = true ]; then
    NWORKERS=2
    # If the user provided a custom srcid_file, use it as-is;
    # otherwise take the first N lines from the default list.
    if [ "${TOTAL_SRCIDS}" -le "${TEST_NSOURCES}" ]; then
        echo ""
        echo " *** TEST MODE: ${TOTAL_SRCIDS} sources (from custom file), ${NWORKERS} workers ***"
        SRCID_FILE_ACTIVE="${SRCID_FILE}"
        TOTAL_SRCIDS="${TOTAL_SRCIDS}"
    else
        echo ""
        echo " *** TEST MODE: ${TEST_NSOURCES} sources, ${NWORKERS} workers ***"
        SRCID_FILE_ACTIVE="${OUTPUT_DIR}/srcids_test.txt"
        head -n "${TEST_NSOURCES}" "${SRCID_FILE}" > "${SRCID_FILE_ACTIVE}"
        TOTAL_SRCIDS="${TEST_NSOURCES}"
    fi
else
    SRCID_FILE_ACTIVE="${SRCID_FILE}"
fi

# Split SRCID file into chunks
CHUNK_DIR="${OUTPUT_DIR}/chunks"
mkdir -p "${CHUNK_DIR}"
# Clean old chunks
rm -f "${CHUNK_DIR}"/chunk_*.txt

# Calculate lines per chunk (ceiling division)
LINES_PER_CHUNK=$(( (TOTAL_SRCIDS + NWORKERS - 1) / NWORKERS ))

echo " Splitting into ${NWORKERS} chunks (~${LINES_PER_CHUNK} SRCIDs each)..."
split -l "${LINES_PER_CHUNK}" -d -a 3 \
    "${SRCID_FILE_ACTIVE}" "${CHUNK_DIR}/chunk_"

# Rename to .txt
for f in "${CHUNK_DIR}"/chunk_*; do
    if [[ ! "${f}" =~ \.txt$ ]]; then
        mv "${f}" "${f}.txt"
    fi
done

NCHUNKS=$(ls "${CHUNK_DIR}"/chunk_*.txt 2>/dev/null | wc -l)
echo " Created ${NCHUNKS} chunk files."

# ==============================================================
# BUILD SKIP FLAG
# ==============================================================
SKIP_FLAG=""
if [ "${DO_RESUME}" = true ]; then
    SKIP_FLAG="--skip_completed"
    echo ""
    echo " RESUME MODE: will skip sources already in fit_results*.fits"
fi

CLEANUP_FLAG=""
if [ "${CLEANUP_CHAINS}" = true ]; then
    CLEANUP_FLAG="--cleanup_chains"
    echo ""
    echo " CLEANUP MODE: chain.fits and corner.png deleted after stats extraction"
fi

# ==============================================================
# LAUNCH PARALLEL JOBS
# ==============================================================
echo ""
echo " Launching ${NCHUNKS} parallel batch jobs..."
echo " Model: ${MODEL}"
echo " Output: ${OUTPUT_DIR}"
echo ""

# Create a log directory for per-chunk stdout/stderr
LOG_DIR="${OUTPUT_DIR}/chunk_logs"
mkdir -p "${LOG_DIR}"

# Record start time (also write to file for monitor_pipeline.sh)
START_TIME=$(date +%s)
echo "${START_TIME}" > "${OUTPUT_DIR}/.pipeline_start_time"

# Launch all chunks in background
PIDS=()
CHUNK_IDX=0

for chunk_file in "${CHUNK_DIR}"/chunk_*.txt; do
    chunk_name=$(basename "${chunk_file}" .txt)
    n_in_chunk=$(wc -l < "${chunk_file}" | tr -d ' ')
    chunk_log="${LOG_DIR}/${chunk_name}.log"
    export_file="fit_results_${chunk_name}.fits"

    echo "   Starting ${chunk_name} (${n_in_chunk} SRCIDs) → ${chunk_log}"

    python3 -u "${REPO_DIR}/automated_fits.py" \
        0 \
        "${DATA_DIR}" \
        "${REPO_DIR}" \
        "${RESPONSES_DIR}" \
        "${OUTPUT_DIR}" \
        "${CATALOG}" \
        "dummy_output.txt" \
        --srcid_file "${chunk_file}" \
        --subdir "${SUBDIR}" \
        --use_bxa \
        --model_name "${MODEL}" \
        --export_results_fits \
        --export_filename "${export_file}" \
        ${SKIP_FLAG} \
        ${CLEANUP_FLAG} \
        > "${chunk_log}" 2>&1 &

    PIDS+=($!)
    CHUNK_IDX=$((CHUNK_IDX + 1))
done

echo ""
echo " All ${NCHUNKS} jobs launched (PIDs: ${PIDS[*]})"
echo ""
echo " Monitor progress with:"
echo "   tail -f ${LOG_DIR}/chunk_000.log"
echo "   bash ${REPO_DIR}/monitor_pipeline.sh ${OUTPUT_DIR}"
echo ""

# ==============================================================
# WAIT FOR ALL JOBS TO COMPLETE
# ==============================================================
echo " Waiting for all jobs to complete..."
echo " (This may take hours/days. If you Ctrl+C or disconnect,"
echo "  the workers continue in background but the merge step"
echo "  below will NOT run. Use recover_results.py to merge later,"
echo "  or relaunch with --resume to re-attach.)"
echo ""

# Trap SIGINT/SIGTERM: skip wait loop but still run the merge
# on whatever partial results exist.
INTERRUPTED=false
trap 'echo ""; echo " Interrupted — skipping to merge step..."; INTERRUPTED=true' INT TERM

N_DONE=0
N_FAILED=0
if [ "${INTERRUPTED}" = false ]; then
    for i in "${!PIDS[@]}"; do
        pid=${PIDS[$i]}
        chunk_file=$(ls "${CHUNK_DIR}"/chunk_*.txt | sed -n "$((i+1))p")
        chunk_name=$(basename "${chunk_file}" .txt)

        if wait "${pid}" 2>/dev/null; then
            N_DONE=$((N_DONE + 1))
            echo "   ${chunk_name} finished OK (${N_DONE}/${NCHUNKS})"
        else
            if [ "${INTERRUPTED}" = true ]; then
                break
            fi
            N_FAILED=$((N_FAILED + 1))
            echo "   ${chunk_name} FAILED (check ${LOG_DIR}/${chunk_name}.log)"
        fi
    done
fi

# Reset trap
trap - INT TERM

END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))
HOURS=$(( ELAPSED / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))

echo ""
echo "=============================================="
if [ "${INTERRUPTED}" = true ]; then
    echo " Pipeline interrupted after ${HOURS}h ${MINS}m"
    echo " Workers still running in background."
else
    echo " Pipeline complete!"
    echo " Time: ${HOURS}h ${MINS}m"
    echo " Chunks OK: ${N_DONE}  |  Failed: ${N_FAILED}"
fi
echo "=============================================="

# ==============================================================
# MERGE RESULTS
# ==============================================================
echo ""
echo " Merging fit results from all chunks..."

OUTPUT_DIR_ESCAPED="${OUTPUT_DIR}" python3 -u << 'MERGE_SCRIPT'
import os
import sys
from astropy.table import Table, vstack

output_dir = os.environ["OUTPUT_DIR_ESCAPED"]

# Find all chunk result files
result_files = sorted([
    os.path.join(output_dir, f) for f in os.listdir(output_dir)
    if f.startswith("fit_results_chunk_") and f.endswith(".fits")
])

if not result_files:
    print("   No result files found to merge.")
    print("   If workers are still running, wait for them to finish")
    print("   then run: python3 recover_results.py " + output_dir)
    sys.exit(0)

print(f"   Found {len(result_files)} result files")

tables = []
for rf in result_files:
    try:
        t = Table.read(rf)
        tables.append(t)
        print(f"   {os.path.basename(rf)}: {len(t)} rows")
    except Exception as e:
        print(f"   WARNING: Could not read {rf}: {e}")

if tables:
    merged = vstack(tables)
    merged_path = os.path.join(output_dir, "fit_results_all.fits")
    merged.write(merged_path, format='fits', overwrite=True)
    print(f"\n   Merged {len(merged)} results -> {merged_path}")
else:
    print("   No valid result tables to merge.")
MERGE_SCRIPT

echo ""
echo " Done! Results in: ${OUTPUT_DIR}/fit_results_all.fits"
echo " Finished: $(date)"
