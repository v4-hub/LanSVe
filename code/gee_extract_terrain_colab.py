"""
GEE Terrain Extraction via Google Colab — Batch Processing with Resume Support
===============================================================================
For: Vegetation Paradox Paper — Communications Earth & Environment

Usage:
  1. Upload 'landslide_locations_for_gee.csv' to Google Drive
  2. Open this notebook in Google Colab
  3. Run all cells — it will auto-resume from the last completed batch

Output: 'landslide_terrain_variables.csv' on Google Drive
"""

# ====== Cell 1: Setup & Authentication ======
import ee
import pandas as pd
import numpy as np
import os
import time
from google.colab import drive

# Mount Google Drive
drive.mount('/content/drive')

# Authenticate & Initialize GEE
ee.Authenticate()
ee.Initialize(project='treering-484000')

print("✅ GEE initialized successfully")

# ====== Cell 2: Configuration ======

# Input/Output paths on Google Drive
DRIVE_DIR = '/content/drive/MyDrive/landslide_terrain'  # Change if needed
INPUT_CSV = os.path.join(DRIVE_DIR, 'landslide_locations_for_gee.csv')
OUTPUT_CSV = os.path.join(DRIVE_DIR, 'landslide_terrain_variables.csv')
PROGRESS_CSV = os.path.join(DRIVE_DIR, '_terrain_progress.csv')  # Partial results

BATCH_SIZE = 500       # Points per batch (safe for GEE limits)
MAX_RETRIES = 3        # Retries per batch on failure
RETRY_DELAY = 10       # Seconds between retries

# Create output directory if needed
os.makedirs(DRIVE_DIR, exist_ok=True)

print(f"📂 Working directory: {DRIVE_DIR}")
print(f"📥 Input: {INPUT_CSV}")
print(f"📤 Output: {OUTPUT_CSV}")

# ====== Cell 3: Load Input Data ======

df_input = pd.read_csv(INPUT_CSV)
print(f"Total points loaded: {len(df_input)}")
print(f"Columns: {list(df_input.columns)}")

# Filter valid coordinates
df_valid = df_input.dropna(subset=['Latitude', 'Longitude']).copy()
df_valid = df_valid[df_valid['Latitude'].apply(lambda x: isinstance(x, (int, float)))]
print(f"Valid coordinate points: {len(df_valid)}")

# ====== Cell 4: Define Terrain Extraction Function ======

def build_terrain_stack():
    """Build the multi-band terrain image stack."""
    # SRTM 30m DEM
    srtm = ee.Image('USGS/SRTMGL1_003')
    elevation = srtm.select('elevation')

    # Basic terrain derivatives
    slope = ee.Terrain.slope(elevation)
    aspect = ee.Terrain.aspect(elevation)
    hillshade = ee.Terrain.hillshade(elevation)

    # TPI at 500m
    kernel_500 = ee.Kernel.circle(radius=500, units='meters')
    elev_mean_500 = elevation.reduceNeighborhood(
        reducer=ee.Reducer.mean(), kernel=kernel_500
    )
    tpi_500 = elevation.subtract(elev_mean_500).rename('TPI_500m')

    # TPI at 1km
    kernel_1k = ee.Kernel.circle(radius=1000, units='meters')
    elev_mean_1k = elevation.reduceNeighborhood(
        reducer=ee.Reducer.mean(), kernel=kernel_1k
    )
    tpi_1km = elevation.subtract(elev_mean_1k).rename('TPI_1km')

    # TRI (terrain ruggedness)
    tri = elevation.reduceNeighborhood(
        reducer=ee.Reducer.stdDev(),
        kernel=ee.Kernel.square(radius=3, units='pixels')
    ).rename('TRI')

    # Curvature (Laplacian)
    dx = elevation.gradient().select('x')
    dy = elevation.gradient().select('y')
    dxx = dx.gradient().select('x')
    dyy = dy.gradient().select('y')
    curvature = dxx.add(dyy).rename('Curvature')

    # MERIT Hydro — upslope area
    merit = ee.Image('MERIT/Hydro/v1_0_1')
    upa = merit.select('upa').rename('Upslope_Area_km2')

    # TWI
    slope_rad = slope.multiply(np.pi / 180)
    tan_slope = slope_rad.tan().max(0.001)
    upa_m2 = upa.multiply(1e6)
    twi = upa_m2.log().subtract(tan_slope.log()).rename('TWI')

    # Stack all bands
    stack = (elevation.rename('Elevation_m')
             .addBands(slope.rename('Slope_deg'))
             .addBands(aspect.rename('Aspect_deg'))
             .addBands(hillshade.rename('Hillshade'))
             .addBands(tpi_500)
             .addBands(tpi_1km)
             .addBands(tri)
             .addBands(curvature)
             .addBands(upa)
             .addBands(twi))

    return stack

TERRAIN_BANDS = [
    'Elevation_m', 'Slope_deg', 'Aspect_deg', 'Hillshade',
    'TPI_500m', 'TPI_1km', 'TRI', 'Curvature',
    'Upslope_Area_km2', 'TWI'
]

def extract_batch(df_batch, terrain_stack):
    """Extract terrain values for a batch of points."""
    # Create GEE FeatureCollection from DataFrame
    features = []
    for _, row in df_batch.iterrows():
        geom = ee.Geometry.Point([float(row['Longitude']), float(row['Latitude'])])
        props = {
            'OBJECTID': int(row['OBJECTID']) if pd.notna(row['OBJECTID']) else -1,
            'ev_id': str(row['ev_id']) if pd.notna(row['ev_id']) else 'NA'
        }
        features.append(ee.Feature(geom, props))

    fc = ee.FeatureCollection(features)

    # Sample terrain at points
    sampled = terrain_stack.reduceRegions(
        collection=fc,
        reducer=ee.Reducer.first(),
        scale=30,
        tileScale=4
    )

    # Fetch results
    results = sampled.getInfo()

    # Parse into DataFrame rows
    rows = []
    for feat in results['features']:
        p = feat['properties']
        row = {
            'OBJECTID': p.get('OBJECTID'),
            'ev_id': p.get('ev_id'),
        }
        for band in TERRAIN_BANDS:
            row[band] = p.get(band)
        rows.append(row)

    return pd.DataFrame(rows)

# ====== Cell 5: Run Batch Extraction with Resume ======

# Build terrain stack once
terrain_stack = build_terrain_stack()
print(f"🗺️ Terrain stack built: {TERRAIN_BANDS}")

# Check for existing progress (resume support)
if os.path.exists(PROGRESS_CSV):
    df_done = pd.read_csv(PROGRESS_CSV)
    done_ids = set(df_done['OBJECTID'].dropna().astype(int))
    print(f"🔄 Resuming: {len(df_done)} points already completed")
else:
    df_done = pd.DataFrame()
    done_ids = set()
    print("🆕 Starting fresh extraction")

# Filter out already-completed points
df_todo = df_valid[~df_valid['OBJECTID'].astype(int).isin(done_ids)].copy()
n_total = len(df_valid)
n_done = len(done_ids)
n_todo = len(df_todo)
print(f"📊 Progress: {n_done}/{n_total} done, {n_todo} remaining")

# Split into batches
batches = [df_todo.iloc[i:i+BATCH_SIZE] for i in range(0, len(df_todo), BATCH_SIZE)]
n_batches = len(batches)
print(f"📦 {n_batches} batches of {BATCH_SIZE} points each\n")

# Process batches
all_results = [df_done] if len(df_done) > 0 else []

for i, batch_df in enumerate(batches):
    batch_num = i + 1
    print(f"[Batch {batch_num}/{n_batches}] Processing {len(batch_df)} points...", end=" ")

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            t0 = time.time()
            batch_result = extract_batch(batch_df, terrain_stack)
            elapsed = time.time() - t0

            all_results.append(batch_result)
            n_done += len(batch_result)

            # Save progress after each batch (crash-safe)
            df_progress = pd.concat(all_results, ignore_index=True)
            df_progress.to_csv(PROGRESS_CSV, index=False)

            print(f"✅ {elapsed:.1f}s | Progress: {n_done}/{n_total} ({100*n_done/n_total:.1f}%)")
            break

        except Exception as e:
            print(f"⚠️ Attempt {attempt}/{MAX_RETRIES} failed: {str(e)[:80]}")
            if attempt < MAX_RETRIES:
                print(f"    Retrying in {RETRY_DELAY}s...")
                time.sleep(RETRY_DELAY)
            else:
                print(f"    ❌ Batch {batch_num} FAILED after {MAX_RETRIES} attempts. Skipping.")
                # Continue to next batch — can re-run later to fill gaps

    # Brief pause between batches to avoid rate limits
    if batch_num < n_batches:
        time.sleep(2)

# ====== Cell 6: Final Output ======

df_final = pd.concat(all_results, ignore_index=True)

# Add Lat/Lon back from input (for verification)
df_coords = df_valid[['OBJECTID', 'Latitude', 'Longitude']].copy()
df_coords['OBJECTID'] = df_coords['OBJECTID'].astype(int)
df_final['OBJECTID'] = df_final['OBJECTID'].astype(int)
df_final = df_final.merge(df_coords, on='OBJECTID', how='left')

# Reorder columns
col_order = ['OBJECTID', 'ev_id', 'Latitude', 'Longitude'] + TERRAIN_BANDS
df_final = df_final[[c for c in col_order if c in df_final.columns]]

# Save final output
df_final.to_csv(OUTPUT_CSV, index=False)

# Clean up progress file
if os.path.exists(PROGRESS_CSV):
    os.remove(PROGRESS_CSV)

print(f"\n{'='*60}")
print(f"✅ EXTRACTION COMPLETE")
print(f"{'='*60}")
print(f"Total points processed: {len(df_final)}")
print(f"Output saved to: {OUTPUT_CSV}")
print(f"\nTerrain summary statistics:")
for band in TERRAIN_BANDS:
    if band in df_final.columns:
        col = df_final[band].dropna()
        if len(col) > 0:
            print(f"  {band:20s}: mean={col.mean():.2f}, median={col.median():.2f}, "
                  f"min={col.min():.2f}, max={col.max():.2f}")

print(f"\n📥 Download from: {OUTPUT_CSV}")
print(f"   Then place in submission_CEE/data/ and run merge_terrain_data.py")
