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
    -R T \
    --debug \
    --seed 42

echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
echo ""
echo "To compare with R output, run:"
echo "python debug_comparison.py example/outputRNA example/outputRNA_fixed"
echo ""
echo "Key files generated:"
echo "- Tree: ${OUTPUT_DIR}/3_CNV.tree.txt"
echo "- LSA results: ${OUTPUT_DIR}/segmental.LSA.txt"
echo "- Permutations: ${OUTPUT_DIR}/permutation/"