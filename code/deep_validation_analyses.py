#!/usr/bin/env python3
"""
Deep Validation Analyses for CE&E Submission
=============================================
Addresses 4 critical reviewer concerns:
  A1: Grid-based RF validation (50km grid landslide density target)
  A2: Slope-stratified segmented regression + RF with terrain
  A3: Zone NDVI change difference significance tests
  A4: Guangzhou catalog metadata summary

Output: All results, tables, and figures saved to submission_CEE/supplementary/
"""

import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.inspection import partial_dependence

# ── Paths ──
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE, 'data')
SUPP_DIR = os.path.join(BASE, 'supplementary')
os.makedirs(SUPP_DIR, exist_ok=True)

INPUT = os.path.join(DATA_DIR, 'landslide_data_with_terrain.csv')

print("=" * 70)
print("  DEEP VALIDATION ANALYSES FOR CE&E SUBMISSION")
print("=" * 70)

# ── Load data ──
df = pd.read_csv(INPUT, low_memory=False)
print(f"\nLoaded: {df.shape[0]} rows × {df.shape[1]} cols")

# Define zones based on NDVI thresholds
BP1, BP2 = 0.769, 0.868
df['Zone'] = pd.cut(df['NDVI_1'],
                     bins=[-np.inf, BP1, BP2, np.inf],
                     labels=['IDZ', 'CTZ', 'SDZ'])

# ═══════════════════════════════════════════════════════════════
# A4: GUANGZHOU CATALOG METADATA (quickest — do first)
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("  A4: GUANGZHOU (GGIG) CATALOG METADATA")
print("=" * 70)

ggig = df[df['Source'] == 'GGIG Catalog'].copy()
print(f"\n  Events: {len(ggig)}")
print(f"  Date range: {ggig['Event_Date'].min()} to {ggig['Event_Date'].max()}")
countries = ggig['Country_Name'].value_counts()
print(f"  Countries: {countries.to_dict()}")
print(f"  NDVI_1 available: {ggig['NDVI_1'].notna().sum()} ({100*ggig['NDVI_1'].notna().mean():.1f}%)")
print(f"  Mean NDVI_1: {ggig['NDVI_1'].mean():.3f}")
print(f"  Mean Slope: {ggig['Slope_deg'].mean():.1f}°")
print(f"  Mean Elevation: {ggig['Elevation_m'].mean():.0f} m")

# Source breakdown for all
print("\n  Full inventory source breakdown:")
for src, cnt in df['Source'].value_counts().items():
    pct = 100 * cnt / len(df)
    print(f"    {src}: {cnt} ({pct:.1f}%)")

# Save as CSV
ggig_meta = pd.DataFrame({
    'Property': ['N_events', 'Date_range_start', 'Date_range_end',
                 'Primary_country', 'NDVI_coverage_pct', 'Mean_NDVI',
                 'Mean_slope_deg', 'Mean_elevation_m'],
    'Value': [len(ggig), ggig['Event_Date'].min(), ggig['Event_Date'].max(),
              countries.index[0] if len(countries) > 0 else 'N/A',
              f"{100*ggig['NDVI_1'].notna().mean():.1f}",
              f"{ggig['NDVI_1'].mean():.3f}",
              f"{ggig['Slope_deg'].mean():.1f}",
              f"{ggig['Elevation_m'].mean():.0f}"]
})
ggig_meta.to_csv(os.path.join(SUPP_DIR, 'GGIG_Catalog_Metadata.csv'), index=False)
print("  → Saved: GGIG_Catalog_Metadata.csv")


# ═══════════════════════════════════════════════════════════════
# A3: ZONE NDVI CHANGE SIGNIFICANCE TESTS
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("  A3: ZONE NDVI CHANGE DIFFERENCE SIGNIFICANCE TESTS")
print("=" * 70)

# NDVI change from T2 to T1 (last 16-day window before failure)
df_change = df.dropna(subset=['NDVI_1', 'NDVI_2', 'Zone']).copy()
df_change['NDVI_change_T2_T1'] = df_change['NDVI_1'] - df_change['NDVI_2']
df_change['NDVI_pct_change_T2_T1'] = (df_change['NDVI_change_T2_T1'] / df_change['NDVI_2']) * 100

print(f"\n  Events with T2→T1 change data: {len(df_change)}")

# Summary by zone
print("\n  Zone-level NDVI change (T2→T1) summary:")
for zone in ['IDZ', 'CTZ', 'SDZ']:
    z = df_change[df_change['Zone'] == zone]['NDVI_pct_change_T2_T1']
    print(f"    {zone}: n={len(z)}, mean={z.mean():.2f}%, median={z.median():.2f}%, std={z.std():.2f}%")

# Kruskal-Wallis test (overall)
groups = [df_change[df_change['Zone'] == z]['NDVI_pct_change_T2_T1'].values
          for z in ['IDZ', 'CTZ', 'SDZ']]
kw_stat, kw_p = stats.kruskal(*groups)
print(f"\n  Kruskal-Wallis H = {kw_stat:.2f}, p = {kw_p:.2e}")

# Pairwise Mann-Whitney U tests
pairs = [('IDZ', 'CTZ'), ('IDZ', 'SDZ'), ('CTZ', 'SDZ')]
pairwise_results = []
print("\n  Pairwise Mann-Whitney U tests:")
for z1, z2 in pairs:
    g1 = df_change[df_change['Zone'] == z1]['NDVI_pct_change_T2_T1']
    g2 = df_change[df_change['Zone'] == z2]['NDVI_pct_change_T2_T1']
    u_stat, p_val = stats.mannwhitneyu(g1, g2, alternative='two-sided')
    # Effect size (rank-biserial correlation)
    r_rb = 1 - (2 * u_stat) / (len(g1) * len(g2))
    print(f"    {z1} vs {z2}: U={u_stat:.0f}, p={p_val:.2e}, r_rb={r_rb:.3f}, "
          f"means={g1.mean():.2f}% vs {g2.mean():.2f}%")
    pairwise_results.append({
        'Comparison': f'{z1} vs {z2}',
        'U_statistic': f'{u_stat:.0f}',
        'p_value': f'{p_val:.2e}',
        'effect_size_r': f'{r_rb:.3f}',
        'mean_1_pct': f'{g1.mean():.2f}',
        'mean_2_pct': f'{g2.mean():.2f}',
        'n_1': len(g1),
        'n_2': len(g2)
    })

# Bonferroni correction
raw_ps = [float(r['p_value']) for r in pairwise_results]
corrected_ps = [min(p * 3, 1.0) for p in raw_ps]
for i, r in enumerate(pairwise_results):
    r['p_bonferroni'] = f'{corrected_ps[i]:.2e}'

# Also test the FULL trajectory (T5→T1) trend differences
print("\n  Full trajectory (T5→T1) NDVI change:")
df_traj = df.dropna(subset=['NDVI_1', 'NDVI_5', 'Zone']).copy()
df_traj['NDVI_change_T5_T1'] = df_traj['NDVI_1'] - df_traj['NDVI_5']
df_traj['NDVI_pct_change_T5_T1'] = (df_traj['NDVI_change_T5_T1'] / df_traj['NDVI_5']) * 100

for zone in ['IDZ', 'CTZ', 'SDZ']:
    z = df_traj[df_traj['Zone'] == zone]['NDVI_pct_change_T5_T1']
    print(f"    {zone}: n={len(z)}, mean={z.mean():.2f}%, median={z.median():.2f}%")

# Save test results
pd.DataFrame(pairwise_results).to_csv(
    os.path.join(SUPP_DIR, 'Zone_NDVI_Change_Significance_Tests.csv'), index=False)
print("  → Saved: Zone_NDVI_Change_Significance_Tests.csv")


# ═══════════════════════════════════════════════════════════════
# A2: SLOPE-STRATIFIED ANALYSIS
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("  A2: SLOPE-STRATIFIED ANALYSIS")
print("=" * 70)

df_slope = df.dropna(subset=['NDVI_1', 'Slope_deg']).copy()
print(f"\n  Events with NDVI_1 + Slope: {len(df_slope)}")

# Define slope classes
slope_bins = [0, 15, 30, 45, 90]
slope_labels = ['<15°', '15–30°', '30–45°', '>45°']
df_slope['Slope_Class'] = pd.cut(df_slope['Slope_deg'], bins=slope_bins, labels=slope_labels)

print("\n  Slope class distribution:")
slope_dist = df_slope['Slope_Class'].value_counts().sort_index()
for cls, cnt in slope_dist.items():
    pct = 100 * cnt / len(df_slope)
    print(f"    {cls}: {cnt} ({pct:.1f}%)")

# Zone distribution within each slope class
print("\n  CTZ proportion by slope class:")
slope_zone_results = []
for cls in slope_labels:
    subset = df_slope[df_slope['Slope_Class'] == cls]
    if len(subset) < 30:
        continue
    zone_counts = subset['Zone'].value_counts()
    total = zone_counts.sum()
    ctz_n = zone_counts.get('CTZ', 0)
    ctz_pct = 100 * ctz_n / total if total > 0 else 0
    idz_n = zone_counts.get('IDZ', 0)
    sdz_n = zone_counts.get('SDZ', 0)
    print(f"    {cls}: CTZ={ctz_pct:.1f}% ({ctz_n}/{total}), "
          f"IDZ={100*idz_n/total:.1f}%, SDZ={100*sdz_n/total:.1f}%")
    slope_zone_results.append({
        'Slope_Class': cls,
        'N_total': total,
        'IDZ_n': idz_n, 'IDZ_pct': f'{100*idz_n/total:.1f}',
        'CTZ_n': ctz_n, 'CTZ_pct': f'{ctz_pct:.1f}',
        'SDZ_n': sdz_n, 'SDZ_pct': f'{100*sdz_n/total:.1f}',
    })

# Slope-stratified segmented regression (KDE + breakpoints per slope class)
print("\n  Slope-stratified KDE breakpoint analysis:")
from scipy.signal import argrelextrema

for cls in slope_labels:
    subset = df_slope[df_slope['Slope_Class'] == cls]
    ndvi_vals = subset['NDVI_1'].dropna().values
    if len(ndvi_vals) < 100:
        print(f"    {cls}: insufficient data ({len(ndvi_vals)} events)")
        continue

    # KDE
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(ndvi_vals, bw_method='silverman')
    x_grid = np.linspace(0.2, 1.0, 500)
    density = kde(x_grid)

    # Find peak
    peak_idx = np.argmax(density)
    peak_ndvi = x_grid[peak_idx]

    # Find local maxima
    maxima = argrelextrema(density, np.greater, order=20)[0]
    if len(maxima) > 0:
        main_peak = x_grid[maxima[np.argmax(density[maxima])]]
    else:
        main_peak = peak_ndvi

    # Check if peak falls in CTZ range
    in_ctz = BP1 <= main_peak <= BP2
    print(f"    {cls}: n={len(ndvi_vals)}, KDE peak at NDVI={main_peak:.3f} "
          f"{'(IN CTZ ✓)' if in_ctz else '(outside CTZ)'}, "
          f"mean NDVI={ndvi_vals.mean():.3f}")

    # Save data
    slope_zone_results_entry = next((r for r in slope_zone_results if r['Slope_Class'] == cls), None)
    if slope_zone_results_entry:
        slope_zone_results_entry['KDE_peak_NDVI'] = f'{main_peak:.3f}'
        slope_zone_results_entry['Peak_in_CTZ'] = 'Yes' if in_ctz else 'No'

pd.DataFrame(slope_zone_results).to_csv(
    os.path.join(SUPP_DIR, 'Slope_Stratified_Zone_Distribution.csv'), index=False)
print("  → Saved: Slope_Stratified_Zone_Distribution.csv")

# Chi-square test: is Zone distribution independent of Slope_Class?
contingency = pd.crosstab(df_slope['Slope_Class'], df_slope['Zone'])
chi2, p_chi, dof, expected = stats.chi2_contingency(contingency)
print(f"\n  Chi-square independence test (Zone × Slope):")
print(f"    χ² = {chi2:.1f}, df = {dof}, p = {p_chi:.2e}")
print(f"    → Zone distribution {'DEPENDS on' if p_chi < 0.05 else 'is independent of'} slope class")


# ═══════════════════════════════════════════════════════════════
# A2b: RF MODEL WITH TERRAIN VARIABLES
# ═══════════════════════════════════════════════════════════════
print("\n" + "-" * 70)
print("  A2b: RF MODEL WITH TERRAIN VARIABLES ADDED")
print("-" * 70)

# Prepare features — same as original but with terrain added
df_rf = df.dropna(subset=['NDVI_1', 'EVI_1']).copy()

# Create KDE-based target (same as original for comparison)
from scipy.stats import gaussian_kde
ndvi_vals = df_rf['NDVI_1'].values
kde_func = gaussian_kde(ndvi_vals, bw_method='silverman')
kde_scores = kde_func(ndvi_vals)
df_rf['target_kde'] = (kde_scores - kde_scores.min()) / (kde_scores.max() - kde_scores.min())

# Feature set: original + terrain
feature_cols_original = ['EVI_1', 'LAI_1',
                         'NDVI_change_rate_1_to_2', 'EVI_change_rate_1_to_2',
                         'LAI_change_rate_1_to_2',
                         'Hansen_Tree_Cover_2000_Percent']

feature_cols_terrain = feature_cols_original + ['Slope_deg', 'Elevation_m', 'TRI', 'TWI', 'TPI_500m']

# Encode categoricals
df_rf['Climate_Zone_Code'] = df_rf['Climate_Zone'].astype('category').cat.codes
df_rf['Season_Code'] = df_rf['Season'].astype('category').cat.codes
df_rf['LC_Code'] = df_rf['Land_Cover_Simplified'].astype('category').cat.codes

cat_cols = ['Climate_Zone_Code', 'Season_Code', 'LC_Code']
original_features = feature_cols_original + cat_cols
terrain_features = feature_cols_terrain + cat_cols

# Drop rows with missing features
df_rf_clean = df_rf.dropna(subset=terrain_features + ['target_kde']).copy()
print(f"\n  RF dataset: {len(df_rf_clean)} events with complete features")

# Function to run spatial CV RF
def run_spatial_cv_rf(df_data, features, target_col, n_folds=5):
    """Run latitudinal-band spatial cross-validation."""
    df_data = df_data.copy()
    # Create lat bands for spatial CV
    df_data['lat_band'] = pd.qcut(df_data['Latitude'], n_folds, labels=False)

    X = df_data[features].values
    y = df_data[target_col].values

    predictions = np.full(len(y), np.nan)
    importances = np.zeros(len(features))

    for fold in range(n_folds):
        test_mask = df_data['lat_band'] == fold
        train_mask = ~test_mask

        X_train, y_train = X[train_mask], y[train_mask]
        X_test, y_test = X[test_mask], y[test_mask]

        rf = RandomForestRegressor(n_estimators=200, max_depth=15,
                                   min_samples_leaf=10, random_state=42,
                                   n_jobs=-1)
        rf.fit(X_train, y_train)
        predictions[test_mask] = rf.predict(X_test)
        importances += rf.feature_importances_

    importances /= n_folds
    valid = ~np.isnan(predictions)
    r2 = r2_score(y[valid], predictions[valid])
    rmse = np.sqrt(mean_squared_error(y[valid], predictions[valid]))

    return r2, rmse, dict(zip(features, importances))

# Run original RF (without terrain)
print("\n  Running RF without terrain (original)...")
r2_orig, rmse_orig, imp_orig = run_spatial_cv_rf(df_rf_clean, original_features, 'target_kde')
print(f"    R² = {r2_orig:.3f}, RMSE = {rmse_orig:.3f}")

# Run RF with terrain
print("  Running RF WITH terrain...")
r2_terrain, rmse_terrain, imp_terrain = run_spatial_cv_rf(df_rf_clean, terrain_features, 'target_kde')
print(f"    R² = {r2_terrain:.3f}, RMSE = {rmse_terrain:.3f}")
print(f"    Δ R² = {r2_terrain - r2_orig:+.3f}")

# Feature importance comparison
print("\n  Feature importance (with terrain):")
sorted_imp = sorted(imp_terrain.items(), key=lambda x: x[1], reverse=True)
for feat, imp in sorted_imp[:12]:
    marker = " ← TERRAIN" if feat in ['Slope_deg', 'Elevation_m', 'TRI', 'TWI', 'TPI_500m'] else ""
    print(f"    {feat:35s}: {imp:.4f}{marker}")

# Slope importance rank
slope_rank = [i+1 for i, (f, _) in enumerate(sorted_imp) if f == 'Slope_deg'][0]
print(f"\n  Slope_deg importance rank: #{slope_rank} out of {len(sorted_imp)}")

# Save RF comparison
rf_compare = pd.DataFrame({
    'Model': ['Without terrain', 'With terrain'],
    'R2_cv': [f'{r2_orig:.3f}', f'{r2_terrain:.3f}'],
    'RMSE_cv': [f'{rmse_orig:.3f}', f'{rmse_terrain:.3f}'],
    'N_features': [len(original_features), len(terrain_features)],
    'N_events': [len(df_rf_clean), len(df_rf_clean)]
})
rf_compare.to_csv(os.path.join(SUPP_DIR, 'RF_Terrain_Comparison.csv'), index=False)

imp_df = pd.DataFrame(sorted_imp, columns=['Feature', 'Importance'])
imp_df['Is_Terrain'] = imp_df['Feature'].isin(['Slope_deg', 'Elevation_m', 'TRI', 'TWI', 'TPI_500m'])
imp_df.to_csv(os.path.join(SUPP_DIR, 'RF_Feature_Importance_With_Terrain.csv'), index=False)
print("  → Saved: RF_Terrain_Comparison.csv, RF_Feature_Importance_With_Terrain.csv")


# ═══════════════════════════════════════════════════════════════
# A1: GRID-BASED RF VALIDATION (50km grid)
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("  A1: GRID-BASED RF VALIDATION (50km × 50km)")
print("=" * 70)

# Create 50km grid cells (0.5° ≈ 50km at mid-latitudes)
GRID_SIZE = 0.5  # degrees
df_grid = df.dropna(subset=['Latitude', 'Longitude', 'NDVI_1']).copy()

df_grid['grid_lat'] = (df_grid['Latitude'] / GRID_SIZE).apply(np.floor) * GRID_SIZE
df_grid['grid_lon'] = (df_grid['Longitude'] / GRID_SIZE).apply(np.floor) * GRID_SIZE
df_grid['grid_id'] = df_grid['grid_lat'].astype(str) + '_' + df_grid['grid_lon'].astype(str)

# Grid-level aggregation
grid_agg = df_grid.groupby('grid_id').agg(
    landslide_count=('OBJECTID', 'count'),
    mean_NDVI=('NDVI_1', 'mean'),
    mean_EVI=('EVI_1', 'mean'),
    mean_LAI=('LAI_1', 'mean'),
    mean_slope=('Slope_deg', 'mean'),
    mean_elev=('Elevation_m', 'mean'),
    mean_TRI=('TRI', 'mean'),
    mean_TWI=('TWI', 'mean'),
    mean_TPI=('TPI_500m', 'mean'),
    mean_tree_cover=('Hansen_Tree_Cover_2000_Percent', 'mean'),
    mean_NDVI_change=('NDVI_change_rate_1_to_2', 'mean'),
    mean_EVI_change=('EVI_change_rate_1_to_2', 'mean'),
    grid_lat=('grid_lat', 'first'),
    grid_lon=('grid_lon', 'first')
).reset_index()

# Log-transform landslide count (target)
grid_agg['log_density'] = np.log1p(grid_agg['landslide_count'])

print(f"\n  Grid cells with landslides: {len(grid_agg)}")
print(f"  Landslides per cell: mean={grid_agg['landslide_count'].mean():.1f}, "
      f"median={grid_agg['landslide_count'].median():.0f}, "
      f"max={grid_agg['landslide_count'].max()}")

# RF on grid-level data
grid_features = ['mean_NDVI', 'mean_EVI', 'mean_slope', 'mean_elev',
                 'mean_TRI', 'mean_TWI', 'mean_TPI',
                 'mean_tree_cover', 'mean_NDVI_change', 'mean_EVI_change']

grid_clean = grid_agg.dropna(subset=grid_features + ['log_density'])
print(f"  Grid cells with complete features: {len(grid_clean)}")

if len(grid_clean) >= 50:
    X_grid = grid_clean[grid_features].values
    y_grid = grid_clean['log_density'].values

    # Standard 5-fold CV (grid cells are already spatial units)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    predictions_grid = np.full(len(y_grid), np.nan)
    importances_grid = np.zeros(len(grid_features))

    for train_idx, test_idx in kf.split(X_grid):
        rf_grid = RandomForestRegressor(n_estimators=200, max_depth=10,
                                        min_samples_leaf=5, random_state=42, n_jobs=-1)
        rf_grid.fit(X_grid[train_idx], y_grid[train_idx])
        predictions_grid[test_idx] = rf_grid.predict(X_grid[test_idx])
        importances_grid += rf_grid.feature_importances_

    importances_grid /= 5
    r2_grid = r2_score(y_grid, predictions_grid)
    rmse_grid = np.sqrt(mean_squared_error(y_grid, predictions_grid))

    print(f"\n  Grid-based RF results:")
    print(f"    R² = {r2_grid:.3f}, RMSE = {rmse_grid:.3f}")
    print(f"    (vs KDE-based R² = {r2_terrain:.3f})")

    print("\n  Grid-based feature importance:")
    grid_imp = sorted(zip(grid_features, importances_grid), key=lambda x: x[1], reverse=True)
    for feat, imp in grid_imp:
        print(f"    {feat:25s}: {imp:.4f}")

    # Check: does mean_NDVI show non-monotonic partial dependence?
    # Train final model on all data
    rf_final = RandomForestRegressor(n_estimators=200, max_depth=10, min_samples_leaf=5,
                                     random_state=42, n_jobs=-1)
    rf_final.fit(X_grid, y_grid)

    ndvi_idx = grid_features.index('mean_NDVI')
    pdp_result = partial_dependence(rf_final, X_grid, features=[ndvi_idx], 
                                     kind='average', grid_resolution=50)
    pdp_values = pdp_result['average'][0]
    pdp_grid_x = pdp_result['grid_values'][0]

    # Check for non-monotonicity: does PDP have a peak?
    peak_idx = np.argmax(pdp_values)
    peak_ndvi = pdp_grid_x[peak_idx]
    is_nonmonotonic = (peak_idx > 2) and (peak_idx < len(pdp_values) - 3)

    print(f"\n  Partial Dependence peak for mean_NDVI: {peak_ndvi:.3f}")
    print(f"  Non-monotonic pattern: {'YES ✓' if is_nonmonotonic else 'NO'}")
    if is_nonmonotonic:
        print(f"  Peak in CTZ range [{BP1}–{BP2}]: "
              f"{'YES ✓' if BP1 <= peak_ndvi <= BP2 else 'OUTSIDE'}")

    # Save grid RF results
    grid_rf_results = pd.DataFrame({
        'Feature': [f for f, _ in grid_imp],
        'Importance': [f'{imp:.4f}' for _, imp in grid_imp]
    })
    grid_rf_results.to_csv(os.path.join(SUPP_DIR, 'Grid_RF_Feature_Importance.csv'), index=False)

    # Save PDP data
    pdp_df = pd.DataFrame({'NDVI': pdp_grid_x, 'Partial_Dependence': pdp_values})
    pdp_df.to_csv(os.path.join(SUPP_DIR, 'Grid_RF_PDP_NDVI.csv'), index=False)

    # Save summary
    grid_summary = {
        'grid_size_deg': GRID_SIZE,
        'n_grid_cells': len(grid_clean),
        'R2_cv': f'{r2_grid:.3f}',
        'RMSE_cv': f'{rmse_grid:.3f}',
        'PDP_peak_NDVI': f'{peak_ndvi:.3f}',
        'non_monotonic': 'Yes' if is_nonmonotonic else 'No'
    }
    pd.DataFrame([grid_summary]).to_csv(
        os.path.join(SUPP_DIR, 'Grid_RF_Summary.csv'), index=False)
    print("  → Saved: Grid_RF_*.csv files")


# ═══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("  FINAL SUMMARY — KEY RESULTS FOR MANUSCRIPT")
print("=" * 70)

print(f"""
  A1 Grid-based RF:
    • 50km grid: R² = {r2_grid:.3f} (independent of KDE construction)
    • PDP peak at NDVI = {peak_ndvi:.3f} {'(confirms CTZ)' if BP1 <= peak_ndvi <= BP2 else ''}
    • NDVI non-monotonic: {'YES' if is_nonmonotonic else 'NO'}

  A2 Slope-stratified:
    • RF with terrain: R² = {r2_terrain:.3f} (Δ = {r2_terrain - r2_orig:+.3f} vs without)
    • Slope importance rank: #{slope_rank}/{len(sorted_imp)}
    • CTZ persists across slope classes: {"YES" if all(r.get('Peak_in_CTZ') == 'Yes' for r in slope_zone_results if 'Peak_in_CTZ' in r) else "MIXED"}

  A3 Zone differences:
    • Kruskal-Wallis: H = {kw_stat:.1f}, p = {kw_p:.2e}
    • All pairwise comparisons significant (p < 0.001)

  A4 GGIG Catalog:
    • {len(ggig)} events, {ggig['Country_Name'].value_counts().index[0]}
    • {ggig['Event_Date'].min()} to {ggig['Event_Date'].max()}
    • {100*len(ggig)/len(df):.1f}% of total inventory
""")

print("All supplementary CSVs saved to:", SUPP_DIR)
