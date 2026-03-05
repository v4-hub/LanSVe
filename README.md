# LanSVe: Landslide–Vegetation Interaction Analysis

**A non-monotonic vegetation–landslide relationship revealed by global satellite observations: critical NDVI thresholds and pre-failure dynamics**

## Overview

This repository contains the data, code, and supplementary materials for a global analysis of landslide–vegetation interactions using satellite-derived vegetation indices. Our analysis of 12,459 rainfall-triggered landslide events (2004–2023) reveals a non-monotonic relationship between vegetation density (NDVI) and landslide susceptibility, challenging the traditional assumption that denser vegetation monotonically stabilizes slopes.

### Key Findings

- **Vegetation Paradox**: Landslide susceptibility peaks at intermediate vegetation density (NDVI 0.769–0.868), defining a Critical Transition Zone (CTZ)
- **Three Biophysical Zones**: Data-driven segmented regression identifies Instability-Dominated (IDZ), Critical Transition (CTZ), and Stability-Dominated (SDZ) zones
- **Pre-failure Dynamics**: NDVI trajectories diverge 16–32 days before slope failure, offering potential early warning signals
- **Global Robustness**: The non-monotonic pattern persists across biomes, climate zones, and slope gradients

## Repository Structure

```
├── code/
│   ├── analysis_v12.R                    # Main statistical analysis (R)
│   ├── deep_validation_analyses.py       # Deep validation analyses (Python)
│   ├── gee_extract_terrain_colab.py      # GEE terrain extraction (Colab)
│   └── gee_extract_terrain_local.py      # GEE terrain extraction (local)
├── data/
│   └── landslide_data_with_terrain.csv   # Processed landslide dataset with terrain attributes
├── figures/
│   ├── fig1.png                          # Data-driven discovery of CTZ
│   ├── fig2.png                          # Multi-scale validation
│   ├── fig3.png                          # Pre-failure vegetation dynamics
│   └── fig4.png                          # Mechanistic framework
└── supplementary/
    ├── Supplementary_Fig_*.png           # 9 supplementary figures
    ├── Supplementary_Table_*.csv         # Supplementary statistical tables
    ├── Grid_RF_*.csv                     # Random Forest model outputs
    ├── RF_*.csv                          # Feature importance results
    ├── GGIG_Catalog_Metadata.csv         # Data catalog metadata
    ├── Slope_Stratified_Zone_Distribution.csv
    └── Zone_NDVI_Change_Significance_Tests.csv
```

## Data Sources

| Dataset | Source | Period |
|---------|--------|--------|
| NASA GLC | NASA Global Landslide Catalog | 2007–2018 |
| GGIG | Global Geo-Intelligence Group | 2004–2023 |
| GFD | Global Fatal Landslide Database | 2004–2017 |
| GFL | Global Fatal Landslide Database (extension) | 2004–2016 |
| MODIS NDVI | MOD13A2.061 (250m, 16-day) | 2004–2023 |
| SRTM DEM | 30m elevation data | — |

## Requirements

### R (≥ 4.3)

Key packages: `tidyverse`, `segmented`, `mgcv`, `ranger`, `sf`, `terra`, `ggplot2`

### Python (≥ 3.9)

Key packages: `pandas`, `numpy`, `scipy`, `scikit-learn`, `ee` (Google Earth Engine)

## Usage

1. **Data Preparation**: Use `gee_extract_terrain_*.py` to extract terrain attributes via Google Earth Engine
2. **Main Analysis**: Run `analysis_v12.R` for segmented regression, validation, and figure generation
3. **Deep Validation**: Run `deep_validation_analyses.py` for additional robustness checks

## Citation

If you use this code or data, please cite:

> Shao et al. A non-monotonic vegetation–landslide relationship revealed by global satellite observations: critical NDVI thresholds and pre-failure dynamics. *One Earth* (under review).

## License

This project is licensed under the MIT License. See individual data source licenses for data usage restrictions.
