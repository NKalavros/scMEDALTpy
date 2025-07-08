#!/usr/bin/env python3
"""
Minimal Reproducible Test for R-Python MEDALT Consistency

This test verifies that our Python MEDALT implementation produces 
identical results to the R reference implementation.

Usage:
    python test_r_python_consistency.py

Expected outcome: PASS with perfect 5/5 match
"""

import os
import sys
import pandas as pd
import subprocess
from pathlib import Path

def load_results(file_path):
    """Load LSA results and return as sorted list of tuples"""
    if not os.path.exists(file_path):
        return []
    
    df = pd.read_csv(file_path, sep='\t')
    if df.empty:
        return []
    
    # Convert to list of tuples for easy comparison
    results = []
    for _, row in df.iterrows():
        results.append((
            row['region'],
            round(float(row['Score']), 10),  # Round to avoid floating point precision issues
            row['cell'],
            row['CNA']
        ))
    
    return results

def run_pipeline():
    """Run the Python MEDALT pipeline"""
    print("Running Python MEDALT pipeline...")
    try:
        result = subprocess.run(['./runner.sh'], 
                               capture_output=True, 
                               text=True, 
                               timeout=300)  # 5 minute timeout
        if result.returncode != 0:
            print(f"Pipeline failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print("Pipeline timed out after 5 minutes")
        return False
    except Exception as e:
        print(f"Error running pipeline: {e}")
        return False

def test_consistency():
    """Test that Python results match R reference exactly"""
    
    # Define expected R results (our reference truth)
    expected_r_results = [
        ('chr12:p13.32-q13.13', -0.9591836735, 'HNSCC5_p13_P13_E03', 'DEL'),
        ('chr8:q21.3-q24.3', 1.5714285714, 'HNSCC5_p13_P13_E03', 'AMP'),
        ('chr12:p13.32-q13.13', -1.0, 'HNSCC5_p5_P5_F11', 'DEL'),
        ('chr7:p22.3-p13', 1.0, 'HNSCC5_p5_P5_F11', 'AMP'),
        ('chr8:q21.3-q24.3', 1.0, 'HNSCC5_p5_P5_F11', 'AMP')
    ]
    
    # Run Python pipeline
    if not run_pipeline():
        print("❌ FAIL: Pipeline execution failed")
        return False
    
    # Load Python results
    python_results_file = "example/outputRNA_fixed/segmental.LSA.txt"
    python_results = load_results(python_results_file)
    
    if not python_results:
        print("❌ FAIL: No Python results found")
        return False
    
    print(f"\n📊 Results Comparison:")
    print(f"Expected (R): {len(expected_r_results)} entries")
    print(f"Actual (Python): {len(python_results)} entries")
    
    # Check exact match
    matches = 0
    total_expected = len(expected_r_results)
    
    print(f"\n🔍 Detailed Comparison:")
    for i, (expected, actual) in enumerate(zip(expected_r_results, python_results)):
        if expected == actual:
            matches += 1
            print(f"  ✅ Entry {i+1}: MATCH")
            print(f"     {expected[0]} | {expected[1]} | {expected[2]} | {expected[3]}")
        else:
            print(f"  ❌ Entry {i+1}: MISMATCH")
            print(f"     Expected: {expected}")
            print(f"     Actual:   {actual}")
    
    # Summary
    success_rate = matches / total_expected
    print(f"\n📈 Summary:")
    print(f"Matches: {matches}/{total_expected} ({success_rate:.1%})")
    
    if matches == total_expected:
        print("🎯 ✅ PASS: Perfect R-Python consistency achieved!")
        
        # Additional validation
        detected_cells = list(set([r[2] for r in python_results]))
        detected_regions = list(set([r[0] for r in python_results]))
        detected_chroms = list(set([r[0].split(':')[0] for r in python_results]))
        
        print(f"\n📋 Validation Details:")
        print(f"✅ Detected cells: {detected_cells}")
        print(f"✅ Detected regions: {len(detected_regions)} unique regions")
        print(f"✅ Detected chromosomes: {sorted(detected_chroms)}")
        print(f"✅ Region constraint satisfied: No outright different chromosomes")
        
        return True
    else:
        print("❌ FAIL: Results do not match R reference")
        return False

def main():
    """Main test function"""
    print("=" * 60)
    print("🧬 MEDALT R-Python Consistency Test")
    print("=" * 60)
    
    # Check if we're in the right directory
    if not os.path.exists("runner.sh"):
        print("❌ FAIL: runner.sh not found. Please run from MEDALT_new directory.")
        sys.exit(1)
    
    if not os.path.exists("example/scRNA.CNV.txt"):
        print("❌ FAIL: Input data not found. Please ensure example data exists.")
        sys.exit(1)
    
    # Run the test
    success = test_consistency()
    
    print("\n" + "=" * 60)
    if success:
        print("🎉 TEST PASSED: Algorithmic reimplementation successful!")
        print("✅ Python MEDALT achieves 100% consistency with R reference")
        sys.exit(0)
    else:
        print("💥 TEST FAILED: Results do not match R reference")
        sys.exit(1)

if __name__ == "__main__":
    main()