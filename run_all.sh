#!/bin/bash
# =============================================================================
# run_all.sh — One-shot reproduction of all figures and tables
#
# Paper: "Landslide susceptibility peaks at intermediate vegetation density:
#         a global satellite analysis of critical NDVI thresholds"
# =============================================================================
set -e
cd "$(dirname "$0")"

echo "============================================================"
echo "  Landslide-Vegetation Analysis — Full Reproduction"
echo "============================================================"

# Step 1: Core breakpoint analysis (~2–4 min, bootstrap 10,000 resamples)
echo ""
echo "[1/5] Core breakpoint analysis (bootstrap CIs)..."
Rscript code/01_core_analysis.R

# Step 2: Sensitivity analyses + contamination test
echo ""
echo "[2/5] Sensitivity analyses + contamination test..."
Rscript code/02_sensitivity_and_contamination.R

# Step 3: Main-text figures
echo ""
echo "[3/5] Main figures (Fig 1 and Fig 3)..."
Rscript code/03_main_figures.R

# Step 4: Supplementary figures
echo ""
echo "[4/5] Supplementary figures (Fig S9–S12) + Table S14..."
Rscript code/04_supplementary_figures.R

# Step 5: Grid RF validation (Tables S6–S12)
echo ""
echo "[5/5] Grid RF validation..."
python3 code/05_grid_rf_validation.py

echo ""
echo "============================================================"
echo "  DONE — all outputs written to figures/ and supplementary/"
echo "============================================================"
echo ""
echo "  Key result (from models_rev/bootstrap_ci_ndvi.csv):"
Rscript -e "
  x <- read.csv('models_rev/bootstrap_ci_ndvi.csv')
  cat(sprintf('  BP1 = %.3f [%.3f, %.3f] BCa 95%%\n', x\$Estimate[1], x\$CI_Lower[1], x\$CI_Upper[1]))
  cat(sprintf('  BP2 = %.3f [%.3f, %.3f] BCa 95%%\n', x\$Estimate[2], x\$CI_Lower[2], x\$CI_Upper[2]))
"
