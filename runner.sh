#!/usr/bin/env zsh
# Test runner for fixed MEDALT Python pipeline

# Set paths
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PACKAGE_PATH="${SCRIPT_DIR}"
# Change working dir to project root so Python can resolve the src package
cd "${SCRIPT_DIR}"
RNA_INPUT="${SCRIPT_DIR}/example/scRNA.CNV.txt"
OUTPUT_DIR="${SCRIPT_DIR}/example/outputRNA_fixed"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "Running fixed MEDALT Python pipeline"
echo "=========================================="
echo "Package path: ${PACKAGE_PATH}"
echo "Input: ${RNA_INPUT}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# Run with debugging enabled and with permutation
python -m src.SC1_py_sctree \
    -P "${PACKAGE_PATH}" \
    -I "${RNA_INPUT}" \
    -G hg19 \
    -O "${OUTPUT_DIR}" \
    -D R \
    -W 30 \
    -R F \
    --debug \
    --seed 42

echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
echo ""

# Check if analysis was successful before creating additional visualizations
if [ -f "${OUTPUT_DIR}/3_CNV.tree.txt" ]; then
    echo "Creating additional visualizations..."
    echo ""
    
    # Create visualizations using the standalone script
    python visualize_medalt.py "${OUTPUT_DIR}/3_CNV.tree.txt" \
        --lsa-file "${OUTPUT_DIR}/segmental.LSA.txt" \
        --output-path "${OUTPUT_DIR}"
    
    echo ""
    echo "=========================================="
    echo "All visualizations complete!"
    echo "=========================================="
    echo ""
    
    echo "Generated files:"
    echo "Data files:"
    echo "- Binned CNV data: ${OUTPUT_DIR}/2_scRNA.CNV_bin_30.csv"
    echo "- Tree structure: ${OUTPUT_DIR}/3_CNV.tree.txt"
    echo "- LSA results: ${OUTPUT_DIR}/segmental.LSA.txt"
    echo "- Permutations: ${OUTPUT_DIR}/permutation/"
    echo ""
    echo "Visualization files:"
    echo "- Single cell tree: ${OUTPUT_DIR}/singlecell_tree.pdf"
    echo "- LSA tree: ${OUTPUT_DIR}/LSA_tree.pdf"
    echo "- Comprehensive report: ${OUTPUT_DIR}/MEDALT_visualization_report.pdf"
    echo ""
    echo "To compare with R output, run:"
    echo "python debug_comparison.py example/outputRNA example/outputRNA_fixed"
    echo ""
    
    # Compare with actual R output
    echo "=========================================="
    echo "🔍 R-Python Results Comparison"
    echo "=========================================="
    
    R_OUTPUT_FILE="${SCRIPT_DIR}/example/outputRNAT/segmental.LSA.txt"
    PYTHON_OUTPUT_FILE="${OUTPUT_DIR}/segmental.LSA.txt"
    
    if [ -f "$R_OUTPUT_FILE" ] && [ -f "$PYTHON_OUTPUT_FILE" ]; then
        echo "Comparing Python output with actual R reference..."
        
        # Count entries
        R_COUNT=$(tail -n +2 "$R_OUTPUT_FILE" | wc -l | tr -d ' ')
        PYTHON_COUNT=$(tail -n +2 "$PYTHON_OUTPUT_FILE" | wc -l | tr -d ' ')
        
        echo "📊 Entry counts:"
        echo "  R output: ${R_COUNT} entries"
        echo "  Python output: ${PYTHON_COUNT} entries"
        
        # Show top cells from each
        echo ""
        echo "🔍 Top cells comparison:"
        echo "R top cells:"
        tail -n +2 "$R_OUTPUT_FILE" | head -5 | cut -f5 | sort -u | sed 's/^/  /'
        echo "Python top cells:"
        tail -n +2 "$PYTHON_OUTPUT_FILE" | head -5 | cut -f3 | sort -u | sed 's/^/  /'
        
        # Show top regions
        echo ""
        echo "🧬 Region comparison:"
        echo "R top regions:"
        tail -n +2 "$R_OUTPUT_FILE" | head -5 | cut -f1 | sed 's/^/  /'
        echo "Python top regions:"
        tail -n +2 "$PYTHON_OUTPUT_FILE" | head -5 | cut -f1 | sed 's/^/  /'
        
        # Check if any cells overlap
        R_CELLS=$(tail -n +2 "$R_OUTPUT_FILE" | cut -f5 | sort -u)
        PYTHON_CELLS=$(tail -n +2 "$PYTHON_OUTPUT_FILE" | cut -f3 | sort -u)
        
        echo ""
        echo "📋 Analysis Summary:"
        if [ "$R_COUNT" -eq "$PYTHON_COUNT" ]; then
            echo "✅ Same number of entries ($R_COUNT)"
        else
            echo "❌ Different entry counts: R=$R_COUNT, Python=$PYTHON_COUNT"
        fi
        
        # Check for any common cells
        COMMON_CELLS=""
        for r_cell in $R_CELLS; do
            for p_cell in $PYTHON_CELLS; do
                if [ "$r_cell" = "$p_cell" ]; then
                    COMMON_CELLS="${COMMON_CELLS} $r_cell"
                fi
            done
        done
        
        if [ -n "$COMMON_CELLS" ]; then
            echo "✅ Common cells found:$COMMON_CELLS"
        else
            echo "❌ No common cells found between R and Python outputs"
        fi
        
        echo ""
        echo "⚠️  RESULTS MISMATCH: Python and R outputs are significantly different"
        echo "🔧 Further investigation needed to match actual R implementation"
        
    else
        echo "❌ ERROR: Cannot find R reference file ($R_OUTPUT_FILE) or Python output"
        echo "Please ensure both files exist for comparison"
    fi
    
    echo ""
    echo "=========================================="
    echo "🚀 MEDALT Analysis Complete!"
    echo "=========================================="
    
else
    echo "ERROR: Analysis failed - tree file not generated"
    echo "Check the output above for error messages"
    exit 1
fi