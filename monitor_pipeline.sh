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

OUTPUT_DIR="${1:-/data/scratch/pipeline_output}"

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

# Log directory
LOG_DIR="${OUTPUT_DIR}/chunk_logs"

# ---------- Per-chunk + aggregate counts ----------
# Count from each chunk log individually and sum.
# This does a single grep pass per log (3 patterns via -E),
# much faster than 4 separate grep -rch on the entire log dir.
N_PROCESSED=0
N_SKIPPED=0
N_SUCCESS=0

CHUNK_LINES=""

if [ -d "${LOG_DIR}" ]; then
    for log_file in "${LOG_DIR}"/chunk_*.log; do
        [ -f "${log_file}" ] || continue
        chunk_name=$(basename "${log_file}" .log)

        # Single grep pass: count the 3 key patterns
        n_processed=$(grep -c "Processing SRCID" "${log_file}" 2>/dev/null || true)
        n_skipped=$(grep -c "Skipping SRCID" "${log_file}" 2>/dev/null || true)
        n_fits=$(grep -c "Fit completed successfully" "${log_file}" 2>/dev/null || true)

        # Accumulate totals
        N_PROCESSED=$((N_PROCESSED + n_processed))
        N_SKIPPED=$((N_SKIPPED + n_skipped))
        N_SUCCESS=$((N_SUCCESS + n_fits))

        # Check if still running (look for "Batch processing complete")
        if grep -q "Batch processing complete" "${log_file}" 2>/dev/null; then
            status="DONE"
        elif lsof "${log_file}" >/dev/null 2>&1; then
            status="RUNNING"
        else
            if [ "${n_processed}" -eq 0 ] 2>/dev/null; then
                status="PENDING"
            else
                status="STOPPED?"
            fi
        fi

        # Get last SRCID being processed
        last_srcid=$(grep "Processing SRCID" "${log_file}" 2>/dev/null | tail -1 | sed -n 's/.*SRCID \([0-9]*\).*/\1/p' || true)
        last_srcid=${last_srcid:-—}

        CHUNK_LINES="${CHUNK_LINES}$(printf "   %-12s %10s  proc: %-6s  skip: %-6s  fits: %-4s  last: %s\n" \
            "${chunk_name}" "[${status}]" "${n_processed}" "${n_skipped}" "${n_fits}" "${last_srcid}")\n"
    done
fi

# Total touched = processed + skipped
N_TOUCHED=$((N_PROCESSED + N_SKIPPED))

# Error rate (only among sources that were actually fitted, not skipped)
if [ "${N_PROCESSED}" -gt 0 ]; then
    ERROR_RATE=$(echo "scale=1; 100 - 100 * ${N_SUCCESS} / ${N_PROCESSED}" | bc 2>/dev/null || echo "?")
else
    ERROR_RATE="—"
fi

echo " Overall progress:"
echo "   Total SRCIDs:    ${TOTAL}"
echo "   Touched:         ${N_TOUCHED} (processed: ${N_PROCESSED}, skipped: ${N_SKIPPED})"
echo "   Successful fits: ${N_SUCCESS}"
echo "   Error rate:      ${ERROR_RATE}%"
echo ""

# Print per-chunk status (already collected above)
if [ -n "${CHUNK_LINES}" ]; then
    echo " Per-chunk status:"
    echo " ─────────────────────────────────────────"
    printf "%b" "${CHUNK_LINES}"
    echo ""
fi

# Check for running Python processes
N_RUNNING=$(pgrep -f "automated_fits.py" 2>/dev/null | wc -l | tr -d ' ')
echo " Running processes: ${N_RUNNING}"

# Check for result FITS files already written
N_RESULT_FILES=$(ls "${OUTPUT_DIR}"/fit_results_chunk_*.fits 2>/dev/null | wc -l | tr -d ' ')
if [ "${N_RESULT_FILES}" -gt 0 ]; then
    TOTAL_ROWS=$(python3 -c "
from astropy.table import Table
import os, sys
total = 0
for f in sorted(os.listdir('${OUTPUT_DIR}')):
    if f.startswith('fit_results_chunk_') and f.endswith('.fits'):
        try:
            total += len(Table.read(os.path.join('${OUTPUT_DIR}', f)))
        except: pass
print(total)
" 2>/dev/null || echo "?")
    echo " Result files: ${N_RESULT_FILES} chunks, ${TOTAL_ROWS} rows written"
else
    echo " Result files: none yet"
fi

# Disk usage — use df on the mount point (instant) instead of
# du -sh (walks entire tree, 30+ min on NFS)
MOUNT_POINT=$(df "${OUTPUT_DIR}" 2>/dev/null | tail -1 | awk '{print $NF}')
DISK_USAGE=$(df -h "${OUTPUT_DIR}" 2>/dev/null | tail -1 | awk '{printf "%s used / %s total (%s)", $3, $2, $5}')
echo " Disk usage: ${DISK_USAGE} (on ${MOUNT_POINT})"

# Estimate completion (only if enough data)
if [ "${N_TOUCHED}" -gt 0 ] && [ "${TOTAL}" != "?" ]; then
    # Read start time from marker file written by run_pipeline.sh
    START_FILE="${OUTPUT_DIR}/.pipeline_start_time"
    PIPELINE_START=""
    if [ -f "${START_FILE}" ]; then
        PIPELINE_START=$(cat "${START_FILE}" | tr -d '[:space:]')
    fi
    if [ -n "${PIPELINE_START}" ]; then
        NOW=$(date +%s)
        ELAPSED=$((NOW - PIPELINE_START))
        ELAPSED_INT=${ELAPSED}

        # Only show rates after at least 5 minutes of data
        if [ "${ELAPSED_INT}" -gt 300 ] 2>/dev/null; then
            # bc omits leading zero; printf adds it
            ELAPSED_H=$(printf "%.1f" "$(echo "scale=2; ${ELAPSED} / 3600" | bc 2>/dev/null)" 2>/dev/null || echo "?")

            # Rates: use N*3600/ELAPSED to avoid intermediate division by zero
            PROC_RATE=$(echo "scale=0; ${N_TOUCHED} * 3600 / ${ELAPSED}" | bc 2>/dev/null || echo "?")

            FIT_RATE=$(printf "%.1f" "$(echo "scale=2; ${N_SUCCESS} * 3600 / ${ELAPSED}" | bc 2>/dev/null)" 2>/dev/null || echo "?")

            # ETA based on total processing rate
            REMAINING_SOURCES=$(echo "${TOTAL} - ${N_TOUCHED}" | bc 2>/dev/null || echo 0)
            ETA_H=$(printf "%.1f" "$(echo "scale=2; ${REMAINING_SOURCES} * ${ELAPSED} / ${N_TOUCHED} / 3600" | bc 2>/dev/null)" 2>/dev/null || echo "?")

            # Expected total successful fits
            if [ "${N_PROCESSED}" -gt 0 ]; then
                EXPECTED_FITS=$(echo "scale=0; ${TOTAL} * ${N_SUCCESS} / ${N_PROCESSED}" | bc 2>/dev/null || echo "?")
            else
                EXPECTED_FITS="?"
            fi

            echo ""
            echo " Timing (elapsed: ${ELAPSED_H}h):"
            echo "   Processing rate: ~${PROC_RATE} sources/hour"
            echo "   Fit rate:        ~${FIT_RATE} fits/hour"
            echo "   Expected fits:   ~${EXPECTED_FITS} (of ${TOTAL} total)"
            echo "   Est. remaining:  ~${ETA_H} hours"
        else
            ELAPSED_MIN=$(( (ELAPSED_INT > 0 ? ELAPSED_INT : 0) / 60 ))
            echo ""
            echo " Timing: ${ELAPSED_MIN}min elapsed (rates shown after 5min)"
        fi
    fi
fi

echo ""
echo "=============================================="
