#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# SYSTEM can be given either as the first positional argument
# (bash run_qa_pdfs.sh Ru ...) or as an environment variable (SYSTEM=Ru
# bash run_qa_pdfs.sh ...). The positional argument takes precedence if
# both are given. If a positional SYSTEM is consumed, the remaining
# arguments shift down (i.e. INPUT becomes $1, OUTDIR becomes $2).
if [[ $# -ge 1 && ( "$1" == "Ru" || "$1" == "Zr" || "$1" == "test" ) ]]; then
    SYSTEM="$1"
    shift
fi

if [[ -z "${SYSTEM:-}" ]]; then
    echo "Error: SYSTEM is not set. Pass it as the first argument (Ru|Zr|test) or set \$SYSTEM. Supported systems are: Ru, Zr, test."
    exit 1
fi

# Directory where the input ROOT trees live.
TREES_DIR="${SCRIPT_DIR}/../../trees"

# 1st remaining arg: input ROOT file basename in TREES_DIR (or an
# absolute/relative path you pass in yourself, which is used as-is).
if [[ "$SYSTEM" == "Ru" ]]; then
    INPUT="${1:-data_merged_Ru.root}"
elif [[ "$SYSTEM" == "Zr" ]]; then
    INPUT="${1:-data_merged_Zr.root}"
elif [[ "$SYSTEM" == "test" ]]; then
    TREES_DIR="${SCRIPT_DIR}/../../data"
    INPUT="${1:-test_data.jetTree.root}"
else
    echo "Error: Unknown system '$SYSTEM'. Supported systems are: Ru, Zr, test."
    exit 1
fi

# Determine the species string to pass into make_qa_pdfs.C
# If test system is run, default to Zr calibration unless specified.
SPECIES_ARG="$SYSTEM"
if [[ "$SYSTEM" == "test" ]]; then
    SPECIES_ARG="Zr"
fi

# Resolve INPUT to an absolute path so it survives the later `pushd`.
# If the user passed an absolute path, or a path that already exists
# relative to the current directory, honor it as-is. Otherwise, assume
# it's a basename that lives in TREES_DIR.
if [[ "$INPUT" = /* ]]; then
    INPUT_PATH="$INPUT"
elif [[ -f "$INPUT" ]]; then
    INPUT_PATH="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
else
    INPUT_PATH="${TREES_DIR}/${INPUT}"
fi

if [[ ! -f "$INPUT_PATH" ]]; then
    echo "Error: Input file not found: $INPUT_PATH"
    exit 1
fi

MACRO="${SCRIPT_DIR}/make_qa_pdfs.C"

# Where PDFs will be created (macro uses relative paths pdf/QA, pdf/Jets)
OUTDIR="${2:-$SCRIPT_DIR}"
mkdir -p "$OUTDIR"

echo "----------------------------------------"
echo "Running QA PDF production"
echo "Script dir : $SCRIPT_DIR"
echo "System     : $SYSTEM (Species passed to ROOT: $SPECIES_ARG)"
echo "Input      : $INPUT_PATH"
echo "Macro      : $MACRO"
echo "Output dir : $OUTDIR"
echo "----------------------------------------"

# Run from OUTDIR so pdf/QA and pdf/Jets land there (not wherever you launched from)
pushd "$OUTDIR" >/dev/null

# Fix: Pass both INPUT_PATH and SPECIES_ARG explicitly into the C++ macro
root -l -b -q "${MACRO}(\"${INPUT_PATH}\", \"${SPECIES_ARG}\")"

popd >/dev/null

echo "----------------------------------------"
echo "Done. PDFs are in: ${OUTDIR}/pdf/QA and ${OUTDIR}/pdf/Jets"
echo "----------------------------------------"