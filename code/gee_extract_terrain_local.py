#!/usr/bin/env python3
"""
GEE Terrain Extraction — Local Batch Processing with Resume Support
===================================================================
Extracts 10 terrain variables from SRTM 30m for 18,028 landslide points.
Runs locally using earthengine-api. Supports resume on interruption.
"""

import ee
import pandas as pd
import numpy as np
import os
import sys
import time

# ── Configuration ──
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), 'data')
INPUT_CSV = os.path.join(DATA_DIR, 'landslide_locations_for_gee.csv')
OUTPUT_CSV = os.path.join(DATA_DIR, 'landslide_terrain_variables.csv')
PROGRESS_CSV = os.path.join(DATA_DIR, '_terrain_progress.csv')

BATCH_SIZE = 500
MAX_RETRIES = 3
RETRY_DELAY = 10

# ── Initialize GEE ──
print("Initializing Google Earth Engine...")
ee.Initialize(project='treering-484000')
print("✅ GEE initialized\n")

# ── Terrain bands ──
TERRAIN_BANDS = [
    'Elevation_m', 'Slope_deg', 'Aspect_deg', 'Hillshade',
    'TPI_500m', 'TPI_1km', 'TRI', 'Curvature',
    'Upslope_Area_km2', 'TWI'
]

def build_terrain_stack():
    srtm = ee.Image('USGS/SRTMGL1_003')
    elevation = srtm.select('elevation')
    slope = ee.Terrain.slope(elevation)
    aspect = ee.Terrain.aspect(elevation)
    hillshade = ee.Terrain.hillshade(elevation)

    k500 = ee.Kernel.circle(radius=500, units='meters')
    tpi_500 = elevation.subtract(
        elevation.reduceNeighborhood(ee.Reducer.mean(), k500)
    ).rename('TPI_500m')

    k1k = ee.Kernel.circle(radius=1000, units='meters')
    tpi_1km = elevation.subtract(
        elevation.reduceNeighborhood(ee.Reducer.mean(), k1k)
    ).rename('TPI_1km')

    tri = elevation.reduceNeighborhood(
        ee.Reducer.stdDev(),
        ee.Kernel.square(radius=3, units='pixels')
    ).rename('TRI')

    dx = elevation.gradient().select('x')
    dy = elevation.gradient().select('y')
    curvature = dx.gradient().select('x').add(
        dy.gradient().select('y')
    ).rename('Curvature')

    merit = ee.Image('MERIT/Hydro/v1_0_1')
    upa = merit.select('upa').rename('Upslope_Area_km2')

    slope_rad = slope.multiply(np.pi / 180)
    tan_slope = slope_rad.tan().max(0.001)
    twi = upa.multiply(1e6).log().subtract(tan_slope.log()).rename('TWI')

    return (elevation.rename('Elevation_m')
            .addBands(slope.rename('Slope_deg'))
            .addBands(aspect.rename('Aspect_deg'))
            .addBands(hillshade.rename('Hillshade'))
            .addBands(tpi_500).addBands(tpi_1km)
            .addBands(tri).addBands(curvature)
            .addBands(upa).addBands(twi))


def extract_batch(df_batch, terrain_stack):
    features = []
    for _, row in df_batch.iterrows():
        geom = ee.Geometry.Point([float(row['Longitude']), float(row['Latitude'])])
        props = {
            'OBJECTID': int(row['OBJECTID']) if pd.notna(row['OBJECTID']) else -1,
            'ev_id': str(row['ev_id']) if pd.notna(row['ev_id']) else 'NA'
        }
        features.append(ee.Feature(geom, props))

    fc = ee.FeatureCollection(features)
    sampled = terrain_stack.reduceRegions(
        collection=fc, reducer=ee.Reducer.first(),
        scale=30, tileScale=4
    )
    results = sampled.getInfo()

    rows = []
    for feat in results['features']:
        p = feat['properties']
        row = {'OBJECTID': p.get('OBJECTID'), 'ev_id': p.get('ev_id')}
        for band in TERRAIN_BANDS:
            row[band] = p.get(band)
        rows.append(row)
    return pd.DataFrame(rows)


# ── Main ──
if __name__ == '__main__':
    # Load input
    df_input = pd.read_csv(INPUT_CSV)
    df_valid = df_input.dropna(subset=['Latitude', 'Longitude']).copy()
    print(f"Total valid points: {len(df_valid)}")

    # Build terrain stack
    terrain_stack = build_terrain_stack()
    print(f"🗺️  Terrain stack ready: {len(TERRAIN_BANDS)} bands\n")

    # Resume support
    if os.path.exists(PROGRESS_CSV):
        df_done = pd.read_csv(PROGRESS_CSV)
        done_ids = set(df_done['OBJECTID'].dropna().astype(int))
        print(f"🔄 Resuming: {len(df_done)} points already done")
    else:
        df_done = pd.DataFrame()
        done_ids = set()

    df_todo = df_valid[~df_valid['OBJECTID'].astype(int).isin(done_ids)].copy()
    n_total = len(df_valid)
    n_done = len(done_ids)
    print(f"📊 {n_done}/{n_total} done, {len(df_todo)} remaining\n")

    batches = [df_todo.iloc[i:i+BATCH_SIZE] for i in range(0, len(df_todo), BATCH_SIZE)]
    all_results = [df_done] if len(df_done) > 0 else []

    for i, batch_df in enumerate(batches):
        batch_num = i + 1
        sys.stdout.write(f"[{batch_num}/{len(batches)}] {len(batch_df)} pts... ")
        sys.stdout.flush()

        for attempt in range(1, MAX_RETRIES + 1):
            try:
                t0 = time.time()
                result = extract_batch(batch_df, terrain_stack)
                elapsed = time.time() - t0
                all_results.append(result)
                n_done += len(result)

                # Save progress
                pd.concat(all_results, ignore_index=True).to_csv(PROGRESS_CSV, index=False)

                print(f"✅ {elapsed:.1f}s | {n_done}/{n_total} ({100*n_done/n_total:.1f}%)")
                break
            except Exception as e:
                print(f"⚠️  attempt {attempt}: {str(e)[:80]}")
                if attempt < MAX_RETRIES:
                    time.sleep(RETRY_DELAY)
                else:
                    print(f"   ❌ Skipped")

        if batch_num < len(batches):
            time.sleep(1)

    # Final output
    df_final = pd.concat(all_results, ignore_index=True)
    df_coords = df_valid[['OBJECTID', 'Latitude', 'Longitude']].copy()
    df_coords['OBJECTID'] = df_coords['OBJECTID'].astype(int)
    df_final['OBJECTID'] = df_final['OBJECTID'].astype(int)
    df_final = df_final.merge(df_coords, on='OBJECTID', how='left')

    col_order = ['OBJECTID', 'ev_id', 'Latitude', 'Longitude'] + TERRAIN_BANDS
    df_final = df_final[[c for c in col_order if c in df_final.columns]]
    df_final.to_csv(OUTPUT_CSV, index=False)

    if os.path.exists(PROGRESS_CSV):
        os.remove(PROGRESS_CSV)

    print(f"\n{'='*60}")
    print(f"✅ DONE: {len(df_final)} points → {OUTPUT_CSV}")
    for band in TERRAIN_BANDS:
        if band in df_final.columns:
            col = df_final[band].dropna()
            if len(col) > 0:
                print(f"  {band:20s}: mean={col.mean():.1f}, med={col.median():.1f}, "
                      f"range=[{col.min():.1f}, {col.max():.1f}]")
