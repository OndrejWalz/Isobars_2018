#!/bin/bash

# ──────────────────────────────────────────────────────────────
# run_divide.sh
# Usage:  ./run_divide.sh
# ──────────────────────────────────────────────────────────────

echo "======================================"
echo "   Ru/Zr Jet Spectra Ratio Macro"
echo "======================================"
echo ""
echo "Select operation:"
echo "  1) Divide Ru / Zr      —  MB data"
echo "  2) Divide Ru / Zr      —  HT2 data"
echo "  3) Divide MB / HT2 —  Ru"
echo "  4) Divide MB / HT2 —  Zr"
echo ""
read -p "Enter choice [1-4]: " choice

case $choice in
    1)
        DATA_TYPE="Ru/Zr MB"
        FILE_NUM="hists_data_Ru_MB.root"
        FILE_DEN="hists_data_Zr_MB.root"
        OUTPUT="ratio_Ru_over_Zr_MB.root"
        LABEL_NUM="Ru"
        LABEL_DEN="Zr"
        ;;
    2)
        DATA_TYPE="Ru/Zr HT2"
        FILE_NUM="hists_data_Ru_HT2.root"
        FILE_DEN="hists_data_Zr_HT2.root"
        OUTPUT="ratio_Ru_over_Zr_HT2.root"
        LABEL_NUM="Ru"
        LABEL_DEN="Zr"
        ;;
    3)
        DATA_TYPE="MB / HT2 Ru"
        FILE_NUM="hists_data_Ru_MB.root"
        FILE_DEN="hists_data_Ru_HT2.root"
        OUTPUT="ratio_MB_over_HT2_Ru.root"
        LABEL_NUM="MB"
        LABEL_DEN="HT2"
        ;;
    4)
        DATA_TYPE="MB / HT2 Zr"
        FILE_NUM="hists_data_Zr_MB.root"
        FILE_DEN="hists_data_Zr_HT2.root"
        OUTPUT="ratio_MB_over_HT2_Zr.root"
        LABEL_NUM="MB"
        LABEL_DEN="HT2"
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

# Derive folder name from output root filename (strip .root)
OUTPUT_DIR="${OUTPUT%.root}"

echo ""
echo "--------------------------------------"
echo " Operation    : $DATA_TYPE"
echo " Numerator    : $FILE_NUM"
echo " Denominator  : $FILE_DEN"
echo " Output file  : $OUTPUT"
echo " PDF folder   : $OUTPUT_DIR/"
echo "--------------------------------------"
echo ""

# Check input files exist
for f in "$FILE_NUM" "$FILE_DEN"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

# Create output folder for PDFs
mkdir -p "$OUTPUT_DIR"
echo "Created folder: $OUTPUT_DIR/"
echo ""

# Run the ROOT macro
root4star -l -b -q "divideJetSpectra.C(\"$FILE_NUM\", \"$FILE_DEN\", \"$OUTPUT\", \"$LABEL_NUM\", \"$LABEL_DEN\", \"$OUTPUT_DIR\")"

if [ $? -eq 0 ]; then
    echo ""
    echo "SUCCESS: Output written to $OUTPUT"
    echo "         PDFs saved in $OUTPUT_DIR/"
else
    echo ""
    echo "ERROR: ROOT macro failed. Check output above."
    exit 1
fi