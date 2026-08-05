# LanSVe — Global satellite analysis of landslide susceptibility and vegetation density

Code and data accompanying:

> Yuan S., Wang Q., Gong Q., Gonzalez-Rodriguez M.A., Song G., Zong Y., Hao Y., Wang J., Chen Y.
> **Global satellite analysis reveals landslide susceptibility peaks at intermediate vegetation density.**
> *iScience* (accepted). Manuscript ID ISCIENCE-D-26-03863.

Lead contact: Shaoxiong Yuan (yuanshx@gdas.ac.cn), Guangzhou Institute of Geography, GDAS.

**The complete archive — code, data, main-text figures, supplementary figures S1–S10,
supplementary tables S1–S16, and metadata — is on Zenodo:**
[https://doi.org/10.5281/zenodo.21190727](https://doi.org/10.5281/zenodo.21190727)

This GitHub repository is a **code + data mirror** for convenient `git clone` and issue
tracking. Figures, tables, and other paper materials are not duplicated here — please
read them in the manuscript or download the full deposit from Zenodo.

---

## What's in this repository

```
data/
  landslide_data_with_terrain.csv    18,028 rainfall-triggered landslides (2000–2024)
                                     with pre-failure NDVI/EVI/LAI, land cover, climate
                                     zone, SRTM terrain, and QA flags.

code/
  01_core_analysis.R                 Segmented breakpoint detection + 10,000-resample
                                     BCa bootstrap CIs. Single source of truth for
                                     BP1 = 0.769 and BP2 = 0.868.

  02_sensitivity_and_contamination.R QA / bandwidth / sensor-noise sensitivity, and
                                     the compositing-contamination test that led us to
                                     retract the pre-failure "early-warning" claim.

  03_main_figures.R                  Regenerates the four main-text figures.
  04_supplementary_figures.R         Regenerates supplementary Figures S1–S10.
  05_grid_rf_validation.py           Grid (0.5°) Random Forest validation and terrain
                                     robustness checks; supplementary tables.

  06_gee_data_collection.py          Documents how vegetation and terrain layers were
                                     extracted from Google Earth Engine. Informational
                                     only — outputs are already in data/.

run_all.sh                           Runs the R and Python pipeline end-to-end.
LICENSE_CODE.txt   MIT (code)
LICENSE_DATA.txt   CC-BY-4.0 (dataset — underlying catalogs retain their own licenses)
```

Figures and supplementary tables are produced into `figures/` and `supplementary/`
subdirectories at run time (git-ignored); they are also archived on Zenodo above.

---

## Reproduction

```
R >= 4.3      segmented, randomForest, boot, mgcv, dplyr, ggplot2, patchwork, tidyr,
              scales, treemapify, cowplot, sf, rnaturalearth, ggnewscale, ggalluvial
Python >= 3.10   pandas, numpy, scipy, scikit-learn >= 1.0
```

```bash
git clone https://github.com/v4-hub/LanSVe.git
cd LanSVe
bash run_all.sh
```

Or step-by-step (must run from repo root; scripts 03/04 depend on `models_rev/` from steps 1–2):

```bash
Rscript code/01_core_analysis.R              # ~2–4 min (bootstrap)
Rscript code/02_sensitivity_and_contamination.R
Rscript code/03_main_figures.R
Rscript code/04_supplementary_figures.R
python3 code/05_grid_rf_validation.py        # ~5–10 min
```

---

## Data sources

| Dataset | Source |
|---|---|
| NASA COOLR / Global Landslide Catalog | Kirschbaum et al. (2015) |
| e-ITALICA catalog | Brunetti et al. (2025) |
| GGIG catalog | Guangzhou Institute of Geography, GDAS |
| MODIS MOD13Q1 (NDVI / EVI) | Didan (2015) |
| MODIS MOD15A2H (LAI / FPAR) | Myneni et al. (2015) |
| MODIS MCD12Q1 (land cover) | Friedl & Sulla-Menashe (2019) |
| Sentinel-2 L2A | ESA Copernicus |
| Copernicus CGLS-LC100 | Buchhorn et al. (2020) |
| Hansen Global Forest Change | Hansen et al. (2013) |
| SRTM 30 m DEM | Farr et al. (2007) |
| MERIT Hydro (TWI, upslope area) | Yamazaki et al. (2019) |
| Köppen-Geiger climate zones | Beck et al. (2018) |

Please cite the original data providers when using the derived variables.

---

## Citation

If you use this code or data, please cite both the paper and the Zenodo archive:

> Yuan S. et al. Global satellite analysis reveals landslide susceptibility peaks at
> intermediate vegetation density. *iScience* (2026).

> Yuan S. et al. LanSVe: Landslide Susceptibility and Vegetation — data and code.
> Zenodo. [https://doi.org/10.5281/zenodo.21190727](https://doi.org/10.5281/zenodo.21190727)
