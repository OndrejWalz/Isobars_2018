#!/usr/bin/env bash
#
# run_compare_svd_binning.sh
#
# Runs compare_svd_binning.C for one or more tags, inside the same
# Apptainer container (roounfold.sif) used for the main unfolding macro,
# so libRooUnfold is available on the library path.
#
# Usage:
#   ./run_compare_svd_binning.sh                      # uses default tag list below
#   ./run_compare_svd_binning.sh TAG1 TAG2 TAG3        # runs only the tags you pass
#
# Expects roounfold.sif to sit next to this script (same directory),
# same as your run_unfold.sh setup. Edit SIF below if it lives elsewhere.

set -euo pipefail

########################
# Paths relative to this script
########################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SIF="${SCRIPT_DIR}/roounfold.sif"
MACRO="${SCRIPT_DIR}/compare_svd_binning.C"

########################
# Tags to process
########################

# Default set of tags to loop over if none are given on the command line.
# Edit this list to match the folder names inside your response files.
DEFAULT_TAGS=(
  "R0.2_CENT_0_10_ptlead5"
  "R0.3_CENT_0_10_ptlead5"
  "R0.4_CENT_0_10_ptlead5"
  "R0.2_PERI_60_80_ptlead5"
  "R0.3_PERI_60_80_ptlead5"
  "R0.4_PERI_60_80_ptlead5"
)

if [[ $# -gt 0 ]]; then
  TAGS=("$@")
else
  TAGS=("${DEFAULT_TAGS[@]}")
fi

########################
# Checks
########################

echo "----------------------------------------"
echo "SCRIPT_DIR : $SCRIPT_DIR"
echo "SIF        : $SIF"
echo "Macro      : $MACRO"
echo "Tags       : ${#TAGS[@]}"
printf '             - %s\n' "${TAGS[@]}"
echo "----------------------------------------"

[[ -f "$SIF"   ]] || { echo "ERROR: SIF not found:   $SIF";   exit 1; }
[[ -f "$MACRO" ]] || { echo "ERROR: MACRO not found: $MACRO"; exit 1; }

FAILED_TAGS=()

########################
# Run inside container, once per tag
########################

for tag in "${TAGS[@]}"; do
  echo "=== Running compare_svd_binning(\"${tag}\") ==="

  if apptainer exec -e -B /gpfs01 "$SIF" root -l -b <<EOF
gSystem->Load("libRooUnfold");
.x ${MACRO}+("${tag}");
.q
EOF
  then
    echo "--- OK: ${tag} ---"
  else
    echo "--- FAILED: ${tag} ---" >&2
    FAILED_TAGS+=("$tag")
  fi
  echo
done

if [ "${#FAILED_TAGS[@]}" -gt 0 ]; then
  echo "[warn] the following tags failed:" >&2
  printf '    - %s\n' "${FAILED_TAGS[@]}" >&2
  exit 1
fi

echo "----------------------------------------"
echo "Done. Plots saved under: svd_binning_comparison/"
echo "----------------------------------------"