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

# Count successful fits (chain.fits files on disk)
N_CHAINS=$(find "${OUTPUT_DIR}" -name "chain.fits" -type f 2>/dev/null | wc -l | tr -d ' ')

# Count successful fits from logs (works with --cleanup_chains)
LOG_DIR="${OUTPUT_DIR}/chunk_logs"
N_SUCCESS_LOGS=0
if [ -d "${LOG_DIR}" ]; then
    N_SUCCESS_LOGS=$(grep -rh "Fit completed successfully" "${LOG_DIR}"/ 2>/dev/null | wc -l | tr -d ' ')
fi

# Use whichever count is higher (chains on disk or from logs)
if [ "${N_SUCCESS_LOGS}" -gt "${N_CHAINS}" ]; then
    N_SUCCESS=${N_SUCCESS_LOGS}
else
    N_SUCCESS=${N_CHAINS}
fi

# Error rate
if [ "${N_STARTED}" -gt 0 ]; then
    ERROR_RATE=$(echo "scale=1; 100 * (1 - ${N_SUCCESS} / ${N_STARTED})" | bc 2>/dev/null || echo "?")
else
    ERROR_RATE="—"
fi

echo " Overall progress:"
echo "   Total SRCIDs:    ${TOTAL}"
echo "   Processed:       ${N_STARTED}"
echo "   Successful fits: ${N_SUCCESS}"
echo "   Error rate:      ${ERROR_RATE}%"
echo ""

# Per-chunk progress from log files
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
if [ "${N_STARTED}" -gt 0 ] && [ "${TOTAL}" != "?" ]; then
    # Find elapsed time from the oldest chunk log file
    FIRST_LOG=$(find "${LOG_DIR}" -name "chunk_*.log" -type f -printf '%T@\n' 2>/dev/null | sort -n | head -1)
    if [ -n "${FIRST_LOG}" ]; then
        NOW=$(date +%s)
        ELAPSED=$(echo "${NOW} - ${FIRST_LOG}" | bc 2>/dev/null || echo 0)
        ELAPSED_INT=${ELAPSED%.*}
        if [ "${ELAPSED_INT}" -gt 0 ]; then
            ELAPSED_H=$(echo "scale=1; ${ELAPSED} / 3600" | bc 2>/dev/null || echo "?")

            # Processing rate (all sources including fast failures)
            PROC_RATE=$(echo "scale=0; ${N_STARTED} / (${ELAPSED} / 3600)" | bc 2>/dev/null || echo "?")

            # Successful fit rate
            FIT_RATE=$(echo "scale=1; ${N_SUCCESS} / (${ELAPSED} / 3600)" | bc 2>/dev/null || echo "?")

            # ETA based on processing rate (accounts for error rate)
            REMAINING_SOURCES=$(echo "${TOTAL} - ${N_STARTED}" | bc 2>/dev/null || echo 0)
            ETA_H=$(echo "scale=1; ${REMAINING_SOURCES} / (${N_STARTED} / (${ELAPSED} / 3600))" | bc 2>/dev/null || echo "?")

            # Expected total successful fits
            if [ "${N_STARTED}" -gt 0 ]; then
                EXPECTED_FITS=$(echo "scale=0; ${TOTAL} * ${N_SUCCESS} / ${N_STARTED}" | bc 2>/dev/null || echo "?")
            else
                EXPECTED_FITS="?"
            fi

            echo ""
            echo " Timing (elapsed: ${ELAPSED_H}h):"
            echo "   Processing rate: ~${PROC_RATE} sources/hour"
            echo "   Fit rate:        ~${FIT_RATE} fits/hour"
            echo "   Expected fits:   ~${EXPECTED_FITS} (of ${TOTAL} total)"
            echo "   Est. remaining:  ~${ETA_H} hours"
        fi
    fi
fi

echo ""
echo "=============================================="
