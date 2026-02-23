#!/bin/bash
# ==============================================================
#  monitor_pipeline.sh — Monitor the progress of a running
#  pipeline launched by run_pipeline.sh.
#
#  Usage:
#    bash monitor_pipeline.sh /path/to/output_dir
#    bash monitor_pipeline.sh              # uses default path
#
#    # Auto-refresh every 60 seconds:
#    watch -n 60 bash monitor_pipeline.sh /path/to/output_dir
#
# ==============================================================

OUTPUT_DIR="${1:-/home/mcoriat/Work/XMM/5XMM/pipeline_output}"

if [ ! -d "${OUTPUT_DIR}" ]; then
    echo "ERROR: Output directory not found: ${OUTPUT_DIR}"
    exit 1
fi

echo "=============================================="
echo " Pipeline Progress Monitor"
echo " Output dir: ${OUTPUT_DIR}"
echo " Time: $(date)"
echo "=============================================="
echo ""

# Count total SRCIDs from chunks
CHUNK_DIR="${OUTPUT_DIR}/chunks"
if [ -d "${CHUNK_DIR}" ]; then
    TOTAL=$(cat "${CHUNK_DIR}"/chunk_*.txt 2>/dev/null | wc -l | tr -d ' ')
else
    TOTAL="?"
fi

# Count source directories (each SRCID gets one)
N_STARTED=$(find "${OUTPUT_DIR}" -maxdepth 1 -type d -name '[0-9]*' 2>/dev/null | wc -l | tr -d ' ')

# Count successful fits (chain.fits files)
N_CHAINS=$(find "${OUTPUT_DIR}" -name "chain.fits" -type f 2>/dev/null | wc -l | tr -d ' ')

# Count per model
echo " Overall progress:"
echo "   Total SRCIDs:    ${TOTAL}"
echo "   Started:         ${N_STARTED}"
echo "   Chains produced: ${N_CHAINS}"
echo ""

# Per-chunk progress from log files
LOG_DIR="${OUTPUT_DIR}/chunk_logs"
if [ -d "${LOG_DIR}" ]; then
    echo " Per-chunk status:"
    echo " ─────────────────────────────────────────"
    for log_file in "${LOG_DIR}"/chunk_*.log; do
        [ -f "${log_file}" ] || continue
        chunk_name=$(basename "${log_file}" .log)

        # Count processed sources from log
        n_processed=$(grep -c "Processing SRCID" "${log_file}" 2>/dev/null || echo 0)
        n_skipped=$(grep -c "Skipping SRCID" "${log_file}" 2>/dev/null || echo 0)
        n_errors=$(grep -c "ERROR" "${log_file}" 2>/dev/null || echo 0)

        # Check if still running (look for "Batch processing complete")
        if grep -q "Batch processing complete" "${log_file}" 2>/dev/null; then
            status="DONE"
        elif [ -n "$(fuser "${log_file}" 2>/dev/null)" ]; then
            status="RUNNING"
        else
            # File not being written and not marked complete
            if [ "${n_processed}" -eq 0 ]; then
                status="PENDING"
            else
                status="STOPPED?"
            fi
        fi

        # Get last SRCID being processed
        last_srcid=$(grep "Processing SRCID" "${log_file}" 2>/dev/null | tail -1 | grep -oP 'SRCID \K[0-9]+' || echo "—")

        printf "   %-12s %8s  processed: %-6s  skipped: %-6s  errors: %-4s  last: %s\n" \
            "${chunk_name}" "[${status}]" "${n_processed}" "${n_skipped}" "${n_errors}" "${last_srcid}"
    done
    echo ""
fi

# Check for running Python processes
N_RUNNING=$(pgrep -f "automated_fits.py" 2>/dev/null | wc -l | tr -d ' ')
echo " Running processes: ${N_RUNNING}"

# Disk usage
DISK_USAGE=$(du -sh "${OUTPUT_DIR}" 2>/dev/null | cut -f1)
echo " Disk usage: ${DISK_USAGE}"

# Estimate completion
if [ "${N_CHAINS}" -gt 0 ] && [ "${TOTAL}" != "?" ]; then
    # Find earliest chain.fits to estimate start time
    FIRST_CHAIN=$(find "${OUTPUT_DIR}" -name "chain.fits" -type f -printf '%T@ %p\n' 2>/dev/null | sort -n | head -1 | cut -d' ' -f1)
    if [ -n "${FIRST_CHAIN}" ]; then
        NOW=$(date +%s)
        ELAPSED=$(echo "${NOW} - ${FIRST_CHAIN}" | bc 2>/dev/null || echo 0)
        if [ "${ELAPSED%.*}" -gt 0 ] && [ "${N_CHAINS}" -gt 0 ]; then
            RATE=$(echo "scale=2; ${N_CHAINS} / (${ELAPSED} / 3600)" | bc 2>/dev/null || echo "?")
            REMAINING=$(echo "scale=0; (${TOTAL} - ${N_STARTED}) / (${N_CHAINS} / (${ELAPSED} / 3600))" | bc 2>/dev/null || echo "?")
            echo ""
            echo " Throughput: ~${RATE} chains/hour"
            echo " Est. remaining: ~${REMAINING} hours"
        fi
    fi
fi

echo ""
echo "=============================================="
