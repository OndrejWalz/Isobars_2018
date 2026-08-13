#!/usr/bin/env bash
# Merge DATA outputs (trees + histograms) from job_ru_2018 into one file.
#
# Default behavior:
#   - Assumes this script lives in:  <REPO>/trees
#   - Looks under:                   <REPO>/data/submit/<YYYY-MM-DD>/job_ru_2018/production
#   - Picks the newest <YYYY-MM-DD> directory
#
# Usage:
#   ./merge_data_all.sh [BASE] [OUT_DIR] [OUT_NAME]
#
# Defaults:
#   BASE     = <REPO>/data/submit/<newest_YYYY-MM-DD>
#   OUT_DIR  = <REPO>/trees
#   OUT_NAME = data_merged_Ru.root
#
# Env:
#   MIN_SIZE=5000        # ignore tiny files (bytes)
#   BATCH=100            # files per shard
#   CLEAN_SHARDS=1       # remove each good shard immediately after it's folded into the final file
#   CHECK_INTEGRITY=1    # open each SHARD (not each input file) with ROOT and skip it if it's corrupted
#   ROOT_CMD=root        # override the ROOT executable if it's not on PATH as "root"
#
# Integrity checking now happens at the SHARD level rather than per input
# file: each batch of BATCH files is hadd'd into a shard first, and only
# that shard is opened and validated with ROOT before being folded into the
# final output. This is far cheaper than spawning a ROOT process per input
# file when there are thousands of files -- with 13000 files and BATCH=100,
# that's ~130 ROOT invocations instead of 13000. The trade-off: if a shard
# is bad, we don't automatically know which single input file caused it --
# the whole shard is set aside for you to investigate (see BAD_DIR below).

set -euo pipefail

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

########################
# Locate repo + defaults
########################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SUBMIT_DIR="${REPO_DIR}/data/submit"

latest_run=""
if [[ -d "$SUBMIT_DIR" ]]; then
  latest_run="$(find "$SUBMIT_DIR" -maxdepth 1 -mindepth 1 -type d -printf '%f\n' \
                | sort | tail -n 1 || true)"
fi

if [[ -z "$latest_run" ]]; then
  echo "ERROR: No run directories found under: $SUBMIT_DIR"
  echo "       Expected something like: $SUBMIT_DIR/2025-11-28"
  exit 1
fi

DEFAULT_BASE="${SUBMIT_DIR}/${latest_run}"
#DEFAULT_BASE="${SUBMIT_DIR}/2026-06-05"
DEFAULT_OUT_DIR="${SCRIPT_DIR}"

########################
# Arguments
########################

BASE="${1:-$DEFAULT_BASE}"
OUT_DIR="${2:-$DEFAULT_OUT_DIR}"
OUT_NAME="${3:-data_merged_Ru_new.root}"

MIN_SIZE="${MIN_SIZE:-5000}"
BATCH="${BATCH:-100}"
CLEAN_SHARDS="${CLEAN_SHARDS:-1}"
CHECK_INTEGRITY="${CHECK_INTEGRITY:-1}"
ROOT_CMD="${ROOT_CMD:-root}"

PROD_DIR="$BASE/job_ru_2018/production"

mkdir -p "$OUT_DIR"

# Detect hadd flags (old ROOT may not support -k)
HADD_FLAGS="-f"
if hadd -h 2>&1 | grep -q -- '-k'; then
  HADD_FLAGS="-fk"
fi

FINAL="$OUT_DIR/$OUT_NAME"
LIST="$OUT_DIR/inlist_data.txt"
BAD_DIR="$OUT_DIR/bad_shards"
: > "$LIST"

########################
# Integrity check helper
########################

# A tiny ROOT macro used to validate each shard. Written once to a temp dir
# (its filename must match the function name, per ROOT convention).
CHECK_MACRO_DIR="$(mktemp -d)"
trap 'rm -rf "$CHECK_MACRO_DIR"' EXIT
CHECK_MACRO="$CHECK_MACRO_DIR/check_root_file.C"
cat > "$CHECK_MACRO" <<'EOF'
void check_root_file(const char *path) {
  TFile *f = TFile::Open(path);
  if (!f || f->IsZombie()) {
    printf("RESULT:BAD reason=open_failed_or_zombie\n");
    return;
  }
  if (f->TestBit(TFile::kRecovered)) {
    printf("RESULT:BAD reason=recovered_from_corruption\n");
    f->Close();
    return;
  }
  if (f->GetListOfKeys() == 0 || f->GetListOfKeys()->GetEntries() == 0) {
    printf("RESULT:BAD reason=no_keys\n");
    f->Close();
    return;
  }
  printf("RESULT:OK\n");
  f->Close();
}
EOF

# Opens a ROOT file and makes sure it's not a zombie, wasn't silently
# "recovered" from a truncated/corrupted write, and actually has content.
# Returns 0 (good) or 1 (bad). On failure, optionally prints the raw ROOT
# output so the reason (corrupt file vs. a broken check itself) is visible
# instead of silently flagging everything as bad.
check_root_file() {
  local f="$1"
  local out
  out="$("$ROOT_CMD" -l -b -q "${CHECK_MACRO}(\"${f}\")" 2>&1)"
  if grep -q 'RESULT:OK' <<<"$out"; then
    return 0
  fi
  if [[ "${CHECK_DEBUG:-0}" -eq 1 ]]; then
    log "  [debug] check output for $f:"
    sed 's/^/      /' <<<"$out"
  fi
  return 1
}

log "Merging DATA files from job_ru_2018"
log "  BASE        : $BASE"
log "  PROD_DIR    : $PROD_DIR"
log "  OUT_DIR     : $OUT_DIR"
log "  OUT_NAME    : $OUT_NAME"
log "  Newest run  : $latest_run"
log "  Integrity   : $([[ "$CHECK_INTEGRITY" -eq 1 ]] && echo "enabled (per shard)" || echo "disabled")"
echo

log "Scanning data production files..."
if [[ ! -d "$PROD_DIR" ]]; then
  log "No data production directory found (expected: $PROD_DIR)"
  exit 1
fi

if [[ "$CHECK_INTEGRITY" -eq 1 ]] && ! command -v "$ROOT_CMD" >/dev/null 2>&1; then
  log "ERROR: CHECK_INTEGRITY=1 but '$ROOT_CMD' was not found on PATH."
  log "       Set ROOT_CMD to the correct executable, or run with CHECK_INTEGRITY=0."
  exit 1
fi

if [[ "$CHECK_INTEGRITY" -eq 1 ]]; then
  # Self-test: prove the checker mechanism itself works (macro compiles,
  # TFile/RESULT sentinel is found, etc.) before trusting it on real data.
  # Otherwise a broken checker (bad ROOT version, missing symbol, whatever)
  # would silently flag every real shard as corrupted.
  #
  # NOTE: the self-test file is created via a macro FILE (same invocation
  # style as check_root_file below), not the "-e" flag -- some ROOT builds
  # (e.g. STAR ROOT) silently no-op "-e" without creating anything.
  MAKE_SELFTEST_MACRO="$CHECK_MACRO_DIR/make_selftest.C"
  cat > "$MAKE_SELFTEST_MACRO" <<'EOF'
void make_selftest(const char *path) {
  TFile f(path, "RECREATE");
  TTree t("t", "t");
  int x = 1;
  t.Branch("x", &x, "x/I");
  t.Fill();
  t.Write();
  f.Close();
}
EOF

  selftest_good="$CHECK_MACRO_DIR/selftest_good.root"
  "$ROOT_CMD" -l -b -q "${MAKE_SELFTEST_MACRO}(\"${selftest_good}\")" >/dev/null 2>&1 || true
  if [[ ! -s "$selftest_good" ]]; then
    log "ERROR: could not create a self-test ROOT file with '$ROOT_CMD' -- your ROOT setup looks broken."
    log "       Command used: $ROOT_CMD -l -b -q \"${MAKE_SELFTEST_MACRO}(\\\"${selftest_good}\\\")\""
    log "       Run that manually to see what's going wrong, or CHECK_INTEGRITY=0 to skip checking."
    exit 1
  fi
  if ! check_root_file "$selftest_good"; then
    log "ERROR: integrity checker flagged a known-good self-test file as BAD."
    log "       The checker itself is broken (not your data). Re-run with CHECK_DEBUG=1 for details:"
    CHECK_DEBUG=1 check_root_file "$selftest_good" || true
    log "       You can proceed without checking via CHECK_INTEGRITY=0."
    exit 1
  fi
  if check_root_file "$CHECK_MACRO_DIR/does_not_exist.root"; then
    log "ERROR: integrity checker flagged a nonexistent/bad file as OK -- checker is not discriminating correctly."
    exit 1
  fi
  log "Integrity checker self-test passed."
fi

# Just collect candidate files by size here -- no per-file ROOT check.
# Integrity is validated per-SHARD, after hadd, further down.
has_files=0
while IFS= read -r -d '' f; do
  sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
  if [[ "$sz" -gt "$MIN_SIZE" ]]; then
    echo "$f" >> "$LIST"
    has_files=1
  fi
done < <(find "$PROD_DIR" -type f -name '*.root' -print0 | sort -z)

if [[ "$has_files" -eq 0 ]]; then
  log "No usable ROOT files in data production directory."
  exit 1
fi

N=$(wc -l < "$LIST" | tr -d ' ')
log "Total files to merge: $N  (batch size: $BATCH)"

# Fast path for a single file: it IS the shard, so check it directly.
if [[ "$N" -eq 1 ]]; then
  src=$(head -n1 "$LIST")
  if [[ "$CHECK_INTEGRITY" -eq 1 ]] && ! check_root_file "$src"; then
    mkdir -p "$BAD_DIR"
    cp -f "$src" "$BAD_DIR/"
    rm -f "$LIST"
    log "ERROR: the only candidate file failed the integrity check; copied to $BAD_DIR"
    exit 1
  fi
  cp -f "$src" "$FINAL"
  rm -f "$LIST"
  log "Merged 1 file -> $FINAL"
  exit 0
fi

# Shard to avoid argv-too-long; old-ROOT friendly
rm -f "$OUT_DIR"/data_shard_*.root "$OUT_DIR"/inlist_data.batch.* || true
split -l "$BATCH" -d -a 4 "$LIST" "$OUT_DIR/inlist_data.batch."

shard_files=( "$OUT_DIR"/inlist_data.batch.* )
num_shards="${#shard_files[@]}"

########################
# Build + check each shard
########################

log "Creating and checking $num_shards shard(s)..."
good_shards=()
n_shard_fail_hadd=0
n_shard_fail_check=0
shard_idx=0
for sub in "${shard_files[@]}"; do
  shard=$(printf "%s/data_shard_%04d.root" "$OUT_DIR" "$shard_idx")
  files_in_shard=$(wc -l < "$sub" | tr -d ' ')
  log "  [shard $((shard_idx+1))/$num_shards] merging $files_in_shard files -> $(basename "$shard")"

  hadd_rc=0
  hadd $HADD_FLAGS "$shard" $(tr '\n' ' ' < "$sub") >/dev/null 2>"$OUT_DIR/.hadd_shard_${shard_idx}.log" || hadd_rc=$?

  if [[ "$hadd_rc" -ne 0 ]]; then
    n_shard_fail_hadd=$((n_shard_fail_hadd + 1))
    log "  [shard $((shard_idx+1))/$num_shards] hadd FAILED (rc=$hadd_rc) -- see $OUT_DIR/.hadd_shard_${shard_idx}.log"
    mkdir -p "$BAD_DIR"
    [[ -f "$shard" ]] && mv -f "$shard" "$BAD_DIR/$(basename "$shard")"
    cp -f "$sub" "$BAD_DIR/$(basename "$sub").filelist"
    shard_idx=$((shard_idx + 1))
    continue
  fi
  rm -f "$OUT_DIR/.hadd_shard_${shard_idx}.log"

  if [[ "$CHECK_INTEGRITY" -eq 1 ]]; then
    if ! check_root_file "$shard"; then
      n_shard_fail_check=$((n_shard_fail_check + 1))
      log "  [shard $((shard_idx+1))/$num_shards] FAILED integrity check -- setting aside"
      mkdir -p "$BAD_DIR"
      mv -f "$shard" "$BAD_DIR/$(basename "$shard")"
      cp -f "$sub" "$BAD_DIR/$(basename "$sub").filelist"
      shard_idx=$((shard_idx + 1))
      continue
    fi
  fi

  good_shards+=("$shard")
  shard_idx=$((shard_idx + 1))
done

n_good_shards="${#good_shards[@]}"
log "Shard summary: $n_good_shards good, $n_shard_fail_hadd failed to build, $n_shard_fail_check failed integrity check"
if [[ "$n_shard_fail_hadd" -gt 0 || "$n_shard_fail_check" -gt 0 ]]; then
  log "  Bad shards (and their input file lists) saved under: $BAD_DIR"
fi

if [[ "$n_good_shards" -eq 0 ]]; then
  log "ERROR: no good shards to merge."
  exit 1
fi

########################
# Fold good shards into FINAL
########################

log "Folding $n_good_shards good shard(s) into final output -> $FINAL"
TMP_FINAL="$OUT_DIR/.tmp_${OUT_NAME}"
i=0
for shard in "${good_shards[@]}"; do
  if [[ "$i" -eq 0 ]]; then
    # First good shard: it becomes the initial final file.
    mv -f "$shard" "$FINAL"
    log "  [merge $((i+1))/$n_good_shards] initialized $(basename "$FINAL") from $(basename "$shard")"
  else
    # Merge current FINAL + this shard into a temp file, then replace FINAL.
    rm -f "$TMP_FINAL"
    hadd $HADD_FLAGS "$TMP_FINAL" "$FINAL" "$shard" >/dev/null
    mv -f "$TMP_FINAL" "$FINAL"
    log "  [merge $((i+1))/$n_good_shards] folded $(basename "$shard") into $(basename "$FINAL")"

    if [[ "$CLEAN_SHARDS" -eq 1 ]]; then
      rm -f "$shard"
    fi
  fi
  i=$((i + 1))
done
rm -f "$TMP_FINAL" || true

# Cleanup
rm -f "$OUT_DIR"/inlist_data.txt "$OUT_DIR"/inlist_data.batch.* || true
rm -f "$OUT_DIR"/data_shard_*.root || true

log "Done. Merged $n_good_shards/$num_shards shard(s) (from $N input files) -> $FINAL"
if [[ "$n_shard_fail_hadd" -gt 0 || "$n_shard_fail_check" -gt 0 ]]; then
  log "NOTE: some shards were excluded due to failures -- check $BAD_DIR"
fi