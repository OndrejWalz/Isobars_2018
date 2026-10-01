#!/bin/bash
###############################################################################
# run_Tower_QA.sh
#
# 1) Sets XRootD/ROOT network timeout env vars so hung/dead servers fail fast
#    instead of retrying for minutes.
# 2) Pre-checks every file in the input list with a hard per-file timeout,
#    splitting it into a "good" list (files that opened fine) and a
#    "skipped_files.txt" log (files that timed out / failed to open).
# 3) Runs RunPicoTowerTest.C in small CHUNKS of the good list, each chunk
#    under its own timeout. A chunk that hangs mid-read (a failure mode the
#    pre-check in step 2 cannot catch, since the file DID open fine) is
#    killed after ANALYSIS_TIMEOUT seconds; its files are logged to
#    skipped_files.txt instead of hanging the whole job. All surviving
#    chunk output files are merged into the final output with hadd.
#
# 
# Defaults match RunPicoTowerTest.C's own defaults.
###############################################################################

set -uo pipefail

# --------------------------- Config / arguments -----------------------------
INPUT_LIST="${1:-filelists/test_list.list}"
OUT_ROOT_FILE="${2:-Test_Tower.root}"
PYEAR="${3:-2018}"
TESTING="${4:-false}"

# Where to put the derived file lists / logs (kept alongside the input list)
LIST_DIR="$(dirname "$INPUT_LIST")"
BASE_NAME="$(basename "$INPUT_LIST")"
GOOD_LIST_DIR="./Tower_QA_good_lists"
GOOD_LIST="${GOOD_LIST_DIR}/${BASE_NAME%.*}_good.list"
SKIPPED_LOG="skipped_files.txt"
 
mkdir -p "$GOOD_LIST_DIR"

# Per-file connection/open timeout (seconds) used during the pre-check pass
PRECHECK_TIMEOUT=5

# Number of files to test in parallel during pre-check (tune to taste)
PARALLEL_JOBS=8

# --- Chunked-run settings (protects against mid-stream hangs during Step 3) --
CHUNK_SIZE=3          # files per analysis job; smaller = finer-grained isolation
                      # of a hang, but more DB/calib-table reload overhead per job
ANALYSIS_TIMEOUT=180  # seconds allowed per chunk before it's killed as "hung"

TOWER_QA_TEMP_DIR="./Tower_QA_temp"
mkdir -p "$TOWER_QA_TEMP_DIR"
WORK_DIR="$(mktemp -d "${TOWER_QA_TEMP_DIR}/.tower_qa_chunks.XXXXXX")"
CHUNK_OUTPUTS=()

cleanup() {
  rm -rf "$WORK_DIR"
}
trap cleanup EXIT

# --------------------------- Step 1: env vars --------------------------------
echo "==> [1/3] Setting XRootD/ROOT network timeout environment variables"

export XRD_CONNECTIONWINDOW=5      # max time (s) to establish a TCP connection
export XRD_STREAMTIMEOUT=5         # max time (s) to wait for a response on an open stream
export XRD_REQUESTTIMEOUT=5        # max time (s) to wait for a single request to complete
export XRD_CONNECTIONRETRY=1       # number of reconnection attempts before giving up
export XRD_TIMEOUTRESOLUTION=3     # granularity (s) at which timeouts are checked

echo "    XRD_CONNECTIONWINDOW=$XRD_CONNECTIONWINDOW"
echo "    XRD_STREAMTIMEOUT=$XRD_STREAMTIMEOUT"
echo "    XRD_REQUESTTIMEOUT=$XRD_REQUESTTIMEOUT"
echo "    XRD_CONNECTIONRETRY=$XRD_CONNECTIONRETRY"
echo "    XRD_TIMEOUTRESOLUTION=$XRD_TIMEOUTRESOLUTION"

# --------------------------- Step 2: pre-check files -------------------------
echo "==> [2/3] Pre-checking files in: $INPUT_LIST"

if [ ! -f "$INPUT_LIST" ]; then
  echo "ERROR: input list '$INPUT_LIST' not found." >&2
  exit 1
fi

: > "$GOOD_LIST"
: > "$SKIPPED_LOG"

# Function that checks a single file: opens it with TFile::Open and verifies
# it is not null / not a zombie. Exits 0 on success, 1 on failure.
check_one_file() {
  local file="$1"
  timeout "$PRECHECK_TIMEOUT" root4star -l -b -q \
    -e "gErrorIgnoreLevel = kFatal+1; TFile *f = TFile::Open(\"$file\"); Int_t ok = (f && !f->IsZombie()) ? 0 : 1; if (f) f->Close(); gSystem->Exit(ok);" \
    > /dev/null 2>&1
  return $?
}
export -f check_one_file
export PRECHECK_TIMEOUT

N_TOTAL=$(grep -c '' "$INPUT_LIST")
echo "    $N_TOTAL files to check (timeout ${PRECHECK_TIMEOUT}s each, ${PARALLEL_JOBS} in parallel)"

# Run checks in parallel, tag each result with GOOD/BAD + filename,
# then split the tagged results into the good list and the skipped log.
# Using xargs for portability (no dependency on GNU parallel).
cat "$INPUT_LIST" | xargs -P "$PARALLEL_JOBS" -I{} bash -c '
  file="{}"
  if check_one_file "$file"; then
    echo "GOOD $file"
  else
    echo "BAD  $file"
  fi
' > /tmp/precheck_results.$$.txt

grep '^GOOD ' /tmp/precheck_results.$$.txt | sed 's/^GOOD //' > "$GOOD_LIST"
grep '^BAD  ' /tmp/precheck_results.$$.txt | sed 's/^BAD  //' >> "$SKIPPED_LOG"
rm -f /tmp/precheck_results.$$.txt

N_GOOD=$(grep -c '' "$GOOD_LIST")
N_BAD_OPEN=$(grep -c '' "$SKIPPED_LOG")

echo "    Good files (opened fine): $N_GOOD"
echo "    Skipped (failed to open): $N_BAD_OPEN   (logged to $SKIPPED_LOG)"

if [ "$N_GOOD" -eq 0 ]; then
  echo "ERROR: no accessible files found, aborting before running the analysis." >&2
  exit 1
fi

# --------------------------- Step 3: chunked analysis run --------------------
echo "==> [3/3] Running RunPicoTowerTest.C in chunks of $CHUNK_SIZE file(s),"
echo "          each capped at ${ANALYSIS_TIMEOUT}s to guard against mid-stream hangs"
echo "          (a failure mode the open-only pre-check in step 2 cannot detect)."

split -l "$CHUNK_SIZE" -d -a 5 "$GOOD_LIST" "${WORK_DIR}/chunk_" --additional-suffix=.list

CHUNK_NUM=0
for chunk_file in "${WORK_DIR}"/chunk_*.list; do
  [ -e "$chunk_file" ] || continue
  CHUNK_NUM=$((CHUNK_NUM + 1))
  chunk_out="${WORK_DIR}/out_${CHUNK_NUM}.root"

  echo "    -> chunk $CHUNK_NUM ($(grep -c '' "$chunk_file") file(s))"

  if timeout "$ANALYSIS_TIMEOUT" root4star -q \
      "RunPicoTowerTest.C(\"${chunk_file}\", \"${chunk_out}\", ${PYEAR}, ${TESTING})" \
      > "${WORK_DIR}/log_${CHUNK_NUM}.txt" 2>&1; then
    if [ -s "$chunk_out" ]; then
      CHUNK_OUTPUTS+=("$chunk_out")
      echo "       OK -> $chunk_out"
    else
      echo "       WARNING: job exited 0 but no output file produced; treating chunk as failed"
      cat "$chunk_file" >> "$SKIPPED_LOG"
    fi
  else
    echo "       TIMED OUT / FAILED (see ${WORK_DIR}/log_${CHUNK_NUM}.txt) -> logging its files as skipped"
    cat "$chunk_file" >> "$SKIPPED_LOG"
  fi
done

N_CHUNK_OK=${#CHUNK_OUTPUTS[@]}
N_SKIPPED_TOTAL=$(grep -c '' "$SKIPPED_LOG" 2>/dev/null || echo 0)

echo "    Chunks completed: $N_CHUNK_OK / $CHUNK_NUM"

if [ "$N_CHUNK_OK" -eq 0 ]; then
  echo "ERROR: every chunk failed or timed out; no output produced." >&2
  echo "       Check logs under: $WORK_DIR (not cleaned up due to failure)"
  trap - EXIT   # keep the work dir around for inspection
  exit 1
fi

echo "==> Merging $N_CHUNK_OK chunk output(s) into: $OUT_ROOT_FILE"
hadd -f "$OUT_ROOT_FILE" "${CHUNK_OUTPUTS[@]}" > "${WORK_DIR}/hadd_log.txt" 2>&1 || {
  echo "ERROR: hadd failed, see ${WORK_DIR}/hadd_log.txt" >&2
  trap - EXIT
  exit 1
}

echo "==> Done. Output: $OUT_ROOT_FILE"
echo "==> Total skipped files (open failures + hung/failed chunks): $N_SKIPPED_TOTAL"
echo "    See $SKIPPED_LOG"