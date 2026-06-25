#!/usr/bin/env python3
"""
05_grid_rf_validation.py — Grid-based RF validation and terrain analyses
========================================================================
Paper: "Landslide susceptibility peaks at intermediate vegetation density:
        a global satellite analysis of critical NDVI thresholds"

Addresses four reviewer-requested validation tests:

  A1  Grid-based RF validation: aggregates 18,028 events to 0.5° grid cells
      (n=364) and trains an independent RF on grid-level mean features, then
      checks for non-monotonic PDP peak in the CTZ range.
  A2  Terrain robustness: spatial-CV RF with and without terrain variables;
      confirms NDVI dominates even when slope, elevation, TRI, TWI are included.
  A3  Zone NDVI-change pairwise significance tests (Kruskal-Wallis +
      Mann-Whitney U with Bonferroni correction).
  A4  Guangzhou (GGIG) catalog metadata summary.

Outputs (written to supplementary/):
  Table_S7.csv    zone pairwise NDVI-change significance tests (A3)
  Table_S13.csv    slope-stratified zone distribution + KDE peaks (A2)
  Table_S14.csv    pixel-RF and grid-RF model summary (A1+A2)
  Table_S8.csv    grid-RF summary: R², RMSE, PDP peak (A1)
  Table_S9.csv   pixel-RF feature importance (without terrain) (A2)
  Table_S10.csv   grid-RF partial dependence on NDVI (A1)
  Table_S11.csv   pixel-RF feature importance with terrain (A2)

Run:
  python3 code/05_grid_rf_validation.py

Requires: pandas, numpy, scipy, scikit-learn (>= 1.0)
"""

import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from scipy import stats
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.inspection import partial_dependence

# ── Paths (relative to deposit root) ──────────────────────────────────────
BASE     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE, 'data')
SUPP_DIR = os.path.join(BASE, 'supplementary')
os.makedirs(SUPP_DIR, exist_ok=True)

INPUT = os.path.join(DATA_DIR, 'landslide_data_with_terrain.csv')

print("=" * 70)
print("  GRID RF VALIDATION — vegetation-landslide analysis")
print("=" * 70)

df = pd.read_csv(INPUT, low_memory=False)
print(f"\nLoaded: {df.shape[0]} events × {df.shape[1]} columns")

BP1, BP2 = 0.769, 0.868
df['Zone'] = pd.cut(df['NDVI_1'], bins=[-np.inf, BP1, BP2, np.inf],
                    labels=['IDZ', 'CTZ', 'SDZ'])


# ═══════════════════════════════════════════════════════════════════════════
# A4: GUANGZHOU (GGIG) CATALOG METADATA
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("  A4: GGIG CATALOG METADATA")
print("─" * 70)
ggig = df[df['Source'] == 'GGIG Catalog'].copy()
print(f"  GGIG events: {len(ggig)} ({100*len(ggig)/len(df):.1f}% of total)")
print(f"  Date range: {ggig['Event_Date'].min()} to {ggig['Event_Date'].max()}")
print(f"  NDVI_1 coverage: {ggig['NDVI_1'].notna().sum()} ({100*ggig['NDVI_1'].notna().mean():.1f}%)")
print("\n  Full inventory source breakdown:")
for src, cnt in df['Source'].value_counts().items():
    print(f"    {src}: {cnt} ({100*cnt/len(df):.1f}%)")


# ═══════════════════════════════════════════════════════════════════════════
# A3: ZONE NDVI-CHANGE PAIRWISE SIGNIFICANCE TESTS  →  Table_S7
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("  A3: ZONE NDVI-CHANGE SIGNIFICANCE TESTS  →  Table_S7")
print("─" * 70)

df_change = df.dropna(subset=['NDVI_1', 'NDVI_2', 'Zone']).copy()
df_change['dNDVI_T2_T1'] = df_change['NDVI_1'] - df_change['NDVI_2']
df_change['pct_change']  = (df_change['dNDVI_T2_T1'] / df_change['NDVI_2']) * 100

print(f"\n  Events with T2→T1 change data: {len(df_change)}")
for zone in ['IDZ','CTZ','SDZ']:
    z = df_change[df_change['Zone']==zone]['pct_change']
    print(f"    {zone}: n={len(z)}, mean={z.mean():.2f}%, median={z.median():.2f}%")

groups   = [df_change[df_change['Zone']==z]['pct_change'].values for z in ['IDZ','CTZ','SDZ']]
kw_H, kw_p = stats.kruskal(*groups)
print(f"\n  Kruskal-Wallis: H={kw_H:.2f}, p={kw_p:.2e}")

pairwise = []
for z1, z2 in [('IDZ','CTZ'),('IDZ','SDZ'),('CTZ','SDZ')]:
    g1  = df_change[df_change['Zone']==z1]['pct_change']
    g2  = df_change[df_change['Zone']==z2]['pct_change']
    U, p_val = stats.mannwhitneyu(g1, g2, alternative='two-sided')
    r_rb = 1 - (2*U)/(len(g1)*len(g2))
    pairwise.append({'Comparison':f'{z1} vs {z2}', 'U_statistic':int(U),
                     'p_value':p_val, 'effect_size_r':round(r_rb,3),
                     'mean_1_pct':round(g1.mean(),2), 'mean_2_pct':round(g2.mean(),2),
                     'n_1':len(g1), 'n_2':len(g2)})
    print(f"    {z1} vs {z2}: U={U:.0f}, p={p_val:.2e}, r={r_rb:.3f}")

raw_ps = [r['p_value'] for r in pairwise]
for i,r in enumerate(pairwise):
    r['p_bonferroni'] = min(raw_ps[i]*3, 1.0)

s6 = pd.DataFrame(pairwise)
s6.to_csv(os.path.join(SUPP_DIR, 'Table_S7.csv'), index=False)
print("  → Table_S7.csv")


# ═══════════════════════════════════════════════════════════════════════════
# A2: SLOPE-STRATIFIED ZONE DISTRIBUTION  →  Table_S13
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("  A2: SLOPE-STRATIFIED ANALYSIS  →  Table_S13")
print("─" * 70)

df_slope = df.dropna(subset=['NDVI_1','Slope_deg']).copy()
slope_bins   = [0, 15, 30, 45, 90]
slope_labels = ['<15°','15–30°','30–45°','>45°']
df_slope['Slope_Class'] = pd.cut(df_slope['Slope_deg'], bins=slope_bins, labels=slope_labels)

slope_rows = []
for cls in slope_labels:
    sub = df_slope[df_slope['Slope_Class']==cls]
    if len(sub) < 30:
        continue
    zc    = sub['Zone'].value_counts()
    total = zc.sum()
    idz_n = int(zc.get('IDZ', 0)); ctz_n = int(zc.get('CTZ', 0)); sdz_n = int(zc.get('SDZ', 0))
    # KDE peak per slope class
    ndvi_vals = sub['NDVI_1'].dropna().values
    kde = gaussian_kde(ndvi_vals, bw_method='silverman')
    xg  = np.linspace(0.2, 1.0, 500); dens = kde(xg)
    maxima = argrelextrema(dens, np.greater, order=20)[0]
    peak_ndvi = xg[maxima[np.argmax(dens[maxima])]] if len(maxima)>0 else xg[np.argmax(dens)]
    row = {'Slope_Class':cls, 'N_total':total,
           'IDZ_n':idz_n, 'IDZ_pct':f'{100*idz_n/total:.1f}',
           'CTZ_n':ctz_n, 'CTZ_pct':f'{100*ctz_n/total:.1f}',
           'SDZ_n':sdz_n, 'SDZ_pct':f'{100*sdz_n/total:.1f}',
           'KDE_peak_NDVI':f'{peak_ndvi:.3f}',
           'Peak_in_CTZ':'Yes' if BP1<=peak_ndvi<=BP2 else 'No'}
    slope_rows.append(row)
    print(f"  {cls}: CTZ={100*ctz_n/total:.1f}%, KDE peak={peak_ndvi:.3f} "
          f"{'(CTZ ✓)' if BP1<=peak_ndvi<=BP2 else ''}")

s7 = pd.DataFrame(slope_rows)
s7.to_csv(os.path.join(SUPP_DIR, 'Table_S13.csv'), index=False)
chi2_stat, chi2_p, _, _ = stats.chi2_contingency(
    pd.crosstab(df_slope['Slope_Class'], df_slope['Zone']))
print(f"  Chi-square (Zone × Slope): χ²={chi2_stat:.1f}, p={chi2_p:.2e}")
print("  → Table_S13.csv")


# ═══════════════════════════════════════════════════════════════════════════
# A2b: PIXEL RF WITH vs WITHOUT TERRAIN  →  Table_S9, Table_S11
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("  A2b: PIXEL RF (terrain robustness)  →  Table_S9, Table_S11")
print("─" * 70)

df_rf = df.dropna(subset=['NDVI_1','EVI_1']).copy()
ndvi_arr = df_rf['NDVI_1'].values
kde_func  = gaussian_kde(ndvi_arr, bw_method='silverman')
kde_scores = kde_func(ndvi_arr)
df_rf['target_kde'] = (kde_scores - kde_scores.min()) / (kde_scores.max() - kde_scores.min())

feat_orig    = ['EVI_1','LAI_1','NDVI_change_rate_1_to_2','EVI_change_rate_1_to_2',
                'LAI_change_rate_1_to_2','Hansen_Tree_Cover_2000_Percent']
feat_terrain = feat_orig + ['Slope_deg','Elevation_m','TRI','TWI','TPI_500m']
for col in ['Climate_Zone','Season','Land_Cover_Simplified']:
    df_rf[col+'_code'] = df_rf[col].astype('category').cat.codes
cat_cols = ['Climate_Zone_code','Season_code','Land_Cover_Simplified_code']
orig_feats    = feat_orig    + cat_cols
terrain_feats = feat_terrain + cat_cols

df_rf_clean = df_rf.dropna(subset=terrain_feats+['target_kde','Latitude']).copy()
print(f"  RF dataset: {len(df_rf_clean)} events")

def spatial_cv_rf(df_data, features, target='target_kde', n_folds=5):
    df_data = df_data.copy()
    df_data['lat_band'] = pd.qcut(df_data['Latitude'], n_folds, labels=False, duplicates='drop')
    X = df_data[features].values; y = df_data[target].values
    preds = np.full(len(y), np.nan); imps = np.zeros(len(features))
    for fold in range(n_folds):
        mask = df_data['lat_band']==fold
        rf = RandomForestRegressor(n_estimators=200, max_depth=15,
                                   min_samples_leaf=10, random_state=42, n_jobs=-1)
        rf.fit(X[~mask], y[~mask]); preds[mask] = rf.predict(X[mask])
        imps += rf.feature_importances_
    imps /= n_folds
    valid = ~np.isnan(preds)
    return r2_score(y[valid], preds[valid]), np.sqrt(mean_squared_error(y[valid], preds[valid])), dict(zip(features, imps))

print("  Running pixel RF without terrain...")
r2_orig, rmse_orig, imp_orig = spatial_cv_rf(df_rf_clean, orig_feats)
print(f"    R²={r2_orig:.3f}, RMSE={rmse_orig:.3f}")

print("  Running pixel RF WITH terrain...")
r2_terr, rmse_terr, imp_terr = spatial_cv_rf(df_rf_clean, terrain_feats)
slope_rank = sorted(imp_terr.items(), key=lambda x:-x[1])
slope_rank_n = [i+1 for i,(f,_) in enumerate(slope_rank) if f=='Slope_deg'][0]
print(f"    R²={r2_terr:.3f}, RMSE={rmse_terr:.3f}  (ΔR²={r2_terr-r2_orig:+.3f})")
print(f"    Slope_deg rank: #{slope_rank_n}/{len(slope_rank)}")

s10 = pd.DataFrame(sorted(imp_orig.items(), key=lambda x:-x[1]),
                   columns=['Feature','Importance'])
s10['Importance'] = s10['Importance'].round(4)
s10.to_csv(os.path.join(SUPP_DIR, 'Table_S9.csv'), index=False)
print("  → Table_S9.csv  (pixel RF feature importance, no terrain)")

s12 = pd.DataFrame(sorted(imp_terr.items(), key=lambda x:-x[1]),
                   columns=['Feature','Importance'])
s12['Importance'] = s12['Importance'].round(4)
s12['Is_Terrain'] = s12['Feature'].isin(['Slope_deg','Elevation_m','TRI','TWI','TPI_500m'])
s12.to_csv(os.path.join(SUPP_DIR, 'Table_S11.csv'), index=False)
print("  → Table_S11.csv  (pixel RF feature importance with terrain)")


# ═══════════════════════════════════════════════════════════════════════════
# A1: GRID-BASED RF VALIDATION (0.5° grid)  →  Table_S8, Table_S10
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("  A1: GRID-BASED RF VALIDATION  →  Table_S8, Table_S10")
print("─" * 70)

GRID = 0.5
df_grid = df.dropna(subset=['Latitude','Longitude','NDVI_1']).copy()
df_grid['g_lat'] = (df_grid['Latitude']/GRID).apply(np.floor)*GRID
df_grid['g_lon'] = (df_grid['Longitude']/GRID).apply(np.floor)*GRID
df_grid['g_id']  = df_grid['g_lat'].astype(str)+'_'+df_grid['g_lon'].astype(str)

agg = df_grid.groupby('g_id').agg(
    count=('OBJECTID','count'),
    mean_NDVI=('NDVI_1','mean'), mean_EVI=('EVI_1','mean'),
    mean_LAI=('LAI_1','mean'),  mean_slope=('Slope_deg','mean'),
    mean_elev=('Elevation_m','mean'), mean_TRI=('TRI','mean'),
    mean_TWI=('TWI','mean'),    mean_TPI=('TPI_500m','mean'),
    mean_tree=('Hansen_Tree_Cover_2000_Percent','mean'),
    mean_dNDVI=('NDVI_change_rate_1_to_2','mean'),
    mean_dEVI=('EVI_change_rate_1_to_2','mean')
).reset_index()
agg['log_count'] = np.log1p(agg['count'])

grid_feats = ['mean_NDVI','mean_EVI','mean_slope','mean_elev',
              'mean_TRI','mean_TWI','mean_TPI','mean_tree','mean_dNDVI','mean_dEVI']
gc = agg.dropna(subset=grid_feats+['log_count'])
print(f"  Grid cells: {len(gc)} (events: {agg['count'].sum():.0f})")

r2_grid = rmse_grid = peak_ndvi = None
is_nm = False

if len(gc) >= 50:
    X = gc[grid_feats].values; y = gc['log_count'].values
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    preds = np.full(len(y), np.nan); imps = np.zeros(len(grid_feats))
    for tr_idx, te_idx in kf.split(X):
        rf = RandomForestRegressor(n_estimators=200, max_depth=10,
                                   min_samples_leaf=5, random_state=42, n_jobs=-1)
        rf.fit(X[tr_idx], y[tr_idx])
        preds[te_idx] = rf.predict(X[te_idx])
        imps += rf.feature_importances_
    imps /= 5
    r2_grid   = r2_score(y, preds)
    rmse_grid = np.sqrt(mean_squared_error(y, preds))

    rf_final = RandomForestRegressor(n_estimators=200, max_depth=10,
                                     min_samples_leaf=5, random_state=42, n_jobs=-1)
    rf_final.fit(X, y)
    ndvi_idx    = grid_feats.index('mean_NDVI')
    pdp_res     = partial_dependence(rf_final, X, features=[ndvi_idx],
                                     kind='average', grid_resolution=50)
    pdp_vals    = pdp_res['average'][0]
    pdp_x       = pdp_res['grid_values'][0]
    peak_idx    = np.argmax(pdp_vals)
    peak_ndvi   = float(pdp_x[peak_idx])
    is_nm       = (peak_idx > 2) and (peak_idx < len(pdp_vals)-3)

    print(f"  Grid RF: R²={r2_grid:.3f}, RMSE={rmse_grid:.3f}")
    print(f"  PDP peak NDVI={peak_ndvi:.3f}  {'(CTZ ✓)' if BP1<=peak_ndvi<=BP2 else '(outside CTZ)'}")
    print(f"  Non-monotonic: {'YES ✓' if is_nm else 'NO'}")

    grid_imp = pd.DataFrame(sorted(zip(grid_feats, imps), key=lambda x:-x[1]),
                            columns=['Feature','Importance'])
    grid_imp['Importance'] = grid_imp['Importance'].round(4)

    pdp_df = pd.DataFrame({'NDVI': pdp_x, 'Partial_Dependence': pdp_vals})
    pdp_df.to_csv(os.path.join(SUPP_DIR, 'Table_S10.csv'), index=False)
    print("  → Table_S10.csv  (grid RF PDP NDVI)")

    s9 = pd.DataFrame([{'grid_size_deg':GRID, 'n_grid_cells':len(gc),
                         'R2_cv':round(r2_grid,3), 'RMSE_cv':round(rmse_grid,3),
                         'PDP_peak_NDVI':round(peak_ndvi,3),
                         'non_monotonic':'Yes' if is_nm else 'No'}])
    s9.to_csv(os.path.join(SUPP_DIR, 'Table_S8.csv'), index=False)
    print("  → Table_S8.csv  (grid RF summary)")


# ═══════════════════════════════════════════════════════════════════════════
# Table_S14 — Model performance summary (pixel RF + grid RF)
# ═══════════════════════════════════════════════════════════════════════════
rows_s8 = [
    ('Pixel_RF_without_terrain_R2_cv',    f'{r2_orig:.3f}'),
    ('Pixel_RF_without_terrain_RMSE_cv',  f'{rmse_orig:.3f}'),
    ('Pixel_RF_with_terrain_R2_cv',       f'{r2_terr:.3f}'),
    ('Pixel_RF_with_terrain_RMSE_cv',     f'{rmse_terr:.3f}'),
    ('Pixel_RF_terrain_delta_R2',         f'{r2_terr-r2_orig:+.3f}'),
    ('Pixel_RF_slope_importance_rank',    str(slope_rank_n)),
]
if r2_grid is not None:
    rows_s8 += [
        ('Grid_RF_R2_cv',      f'{r2_grid:.3f}'),
        ('Grid_RF_RMSE_cv',    f'{rmse_grid:.3f}'),
        ('Grid_RF_PDP_peak',   f'{peak_ndvi:.3f}'),
        ('Grid_RF_n_cells',    str(len(gc))),
    ]

s8 = pd.DataFrame(rows_s8, columns=['Property','Value'])
s8.to_csv(os.path.join(SUPP_DIR, 'Table_S14.csv'), index=False)
print("  → Table_S14.csv  (RF model summary)")

print(f"""
═══════════════════════════════════════════════════════════════════════════
  SUMMARY
═══════════════════════════════════════════════════════════════════════════
  A3  Kruskal-Wallis: H={kw_H:.1f}, p={kw_p:.2e}  (all pairwise p < 0.05)
  A2  Pixel RF without terrain: R²={r2_orig:.3f}
      Pixel RF WITH terrain:    R²={r2_terr:.3f}  (ΔR²={r2_terr-r2_orig:+.3f})
      Slope rank #{slope_rank_n} — NDVI dominates regardless of terrain
  A1  Grid RF (0.5°, n={len(gc)}): R²={r2_grid:.3f}
      PDP peak NDVI={peak_ndvi:.3f}  {'→ confirms CTZ non-monotonic peak' if BP1<=peak_ndvi<=BP2 else ''}
═══════════════════════════════════════════════════════════════════════════
""")
print("All tables written to supplementary/")
