#!/bin/bash
# ==============================================================
#  run_pipeline.sh — Launch the automated spectral fitting
#  pipeline on bell in parallel across multiple cores.
#
#  Usage:
#    # Test run (10 sources, 2 workers):
#    bash run_pipeline.sh --test
#
#    # Full run (all sources, 30 workers):
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

# ==============================================================
# CONFIGURATION — edit these paths for your setup
# ==============================================================
REPO_DIR="/home/mcoriat/Work/XMM/5XMM/automated_fits"
DATA_DIR="/mnt/xmmcat/5XMM_data/Spectra"
RESPONSES_DIR="/home/mcoriat/Work/XMM/5XMM/RESPONSES"
CATALOG="${REPO_DIR}/5xmm_matched_for_pipeline.fits"
OUTPUT_DIR="/home/mcoriat/Work/XMM/5XMM/pipeline_output"
SRCID_FILE="${REPO_DIR}/srcids.txt"

# Default parallel settings
NWORKERS=30          # Number of parallel batch jobs (keep ~10 cores free)
MODEL="powerlaw"     # Spectral model: powerlaw, apec_single, blackbody, bremss
SUBDIR="product"     # Subdirectory under OBS_ID: 'product' (5XMM) or 'pps' (4XMM)
TEST_NSOURCES=10     # Number of sources for test run

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
    export HEADAS=/home/mcoriat/Software/heasoft-6.35.2/x86_64-pc-linux-gnu-libc2.31
    if [ -f "${HEADAS}/headas-init.sh" ]; then
        echo " Initializing HEASOFT..."
        source "${HEADAS}/headas-init.sh"
    else
        echo " ERROR: HEADAS init script not found at:"
        echo "        ${HEADAS}/headas-init.sh"
        exit 1
    fi
else
    echo " HEASOFT already loaded: ${HEADAS}"
fi

# Initialize SAS (may be needed for shared libraries).
# The SAS setup script uses variables that may be unbound,
# so we temporarily relax strict mode.
if [ -z "${SAS_DIR:-}" ]; then
    export SAS_DIR=/home/filippos/SASscreeningDR11/
    export SAS_CCFPATH=/home/filippos/sasbuild/ccf/pub
    export SAS_PATH="${SAS_DIR}"
    if [ -f "${SAS_DIR}/sas-setup.sh" ]; then
        echo " Initializing SAS..."
        set +eu
        source "${SAS_DIR}/sas-setup.sh" 2>/dev/null
        set -eu
        echo " SAS initialized (or skipped gracefully)."
    else
        echo " WARNING: SAS setup not found at ${SAS_DIR}/sas-setup.sh"
        echo "          Continuing without SAS (not needed for BXA fits)."
    fi
else
    echo " SAS already loaded: ${SAS_DIR}"
fi

# Extra library path (Qt libs needed by some HEASOFT/SAS components)
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/sasbuild/local/sasbld03n/GNU_CC_CXX_9.2.0/qt-x11-free/lib/"

# Verify critical dependencies
echo ""
echo " Checking dependencies..."
python3 -c "import xspec; print(f'   XSPEC: {xspec.__file__}')" 2>/dev/null || {
    echo "   ERROR: Cannot import xspec. Is HEASOFT+PyXSPEC set up?"
    exit 1
}
python3 -c "import bxa; print(f'   BXA:   {bxa.__file__}')" 2>/dev/null || {
    echo "   ERROR: Cannot import bxa. Install with: pip install bxa"
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
    echo " RESUME MODE: will skip sources with existing chain.fits"
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

# Record start time
START_TIME=$(date +%s)

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
echo " (This may take hours/days. Safe to Ctrl+C — jobs"
echo "  continue in background. Rerun with --resume to track.)"
echo ""

N_DONE=0
N_FAILED=0
for i in "${!PIDS[@]}"; do
    pid=${PIDS[$i]}
    chunk_file=$(ls "${CHUNK_DIR}"/chunk_*.txt | sed -n "$((i+1))p")
    chunk_name=$(basename "${chunk_file}" .txt)

    if wait "${pid}"; then
        N_DONE=$((N_DONE + 1))
        echo "   ${chunk_name} finished OK (${N_DONE}/${NCHUNKS})"
    else
        N_FAILED=$((N_FAILED + 1))
        echo "   ${chunk_name} FAILED (check ${LOG_DIR}/${chunk_name}.log)"
    fi
done

END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))
HOURS=$(( ELAPSED / 3600 ))
MINS=$(( (ELAPSED % 3600) / 60 ))

echo ""
echo "=============================================="
echo " Pipeline complete!"
echo " Time: ${HOURS}h ${MINS}m"
echo " Chunks OK: ${N_DONE}  |  Failed: ${N_FAILED}"
echo "=============================================="

# ==============================================================
# MERGE RESULTS
# ==============================================================
echo ""
echo " Merging fit results from all chunks..."

python3 -u << 'MERGE_SCRIPT'
import os
import sys
from astropy.table import Table, vstack

output_dir = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "OUTPUT_DIR", ".")

# Find all chunk result files
result_files = sorted([
    os.path.join(output_dir, f) for f in os.listdir(output_dir)
    if f.startswith("fit_results_chunk_") and f.endswith(".fits")
])

if not result_files:
    print("   No result files found to merge.")
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
    merged.write(merged_path, overwrite=True)
    print(f"\n   Merged {len(merged)} results → {merged_path}")
else:
    print("   No valid result tables to merge.")
MERGE_SCRIPT

echo ""
echo " Done! Results in: ${OUTPUT_DIR}/fit_results_all.fits"
echo " Finished: $(date)"
