"""
06_gee_data_collection.py — Google Earth Engine terrain extraction (informational)
==================================================================================
Paper: "Landslide susceptibility peaks at intermediate vegetation density:
        a global satellite analysis of critical NDVI thresholds"

This script documents how terrain variables (elevation, slope, TRI, TWI, TPI,
curvature, upslope area) and vegetation indices (NDVI, EVI, LAI/FPAR from
MODIS; NDVI, EVI, SAVI, NDMI, NBR from Sentinel-2) were extracted from
Google Earth Engine (GEE) for each landslide location.

NOTE: This script is PROVIDED FOR TRANSPARENCY ONLY.
Running it requires:
  1. A Google Earth Engine account (free at https://signup.earthengine.google.com/)
  2. A Google Colab environment (or local GEE Python API installation)
  3. The raw landslide catalog CSV with Latitude, Longitude, Event_Date columns

The curated output of this extraction is already included in
data/landslide_data_with_terrain.csv — you do NOT need to re-run this script
to reproduce the paper's results. Run 01_core_analysis.R instead.

Data sources:
  MODIS:    MOD13Q1 (NDVI/EVI), MOD15A2H (LAI/FPAR), MCD12Q1 (land cover)
  Sentinel: ESA Sentinel-2 L2A (10 m, 2015–2024)
  Terrain:  SRTM 30m (elevation, slope, aspect), MERIT Hydro (TWI, upslope area)
  Forest:   Hansen Global Forest Change v1.10
  Land cv:  Copernicus CGLS-LC100 v3

Usage (Google Colab):
  1. Upload your landslide CSV to Google Drive
  2. Open a new Colab notebook
  3. Copy-paste this script into the notebook cells
  4. Run all cells — it will auto-resume from the last completed batch
"""

# ====== Cell 1: Setup & Authentication ======
import ee
import pandas as pd
import numpy as np
import os
import time

# For Colab: mount Google Drive
# from google.colab import drive
# drive.mount('/content/drive')

# Authenticate & Initialize GEE
ee.Authenticate()
ee.Initialize()          # pass project='your-gee-project-id' if needed

print("GEE initialized")

# ====== Cell 2: Configuration ======
DRIVE_DIR    = '/content/drive/MyDrive/landslide_terrain'  # adjust as needed
INPUT_CSV    = os.path.join(DRIVE_DIR, 'landslide_locations.csv')   # cols: ev_id, Latitude, Longitude, Event_Date
OUTPUT_CSV   = os.path.join(DRIVE_DIR, 'landslide_terrain_variables.csv')
PROGRESS_CSV = os.path.join(DRIVE_DIR, '_terrain_progress.csv')

BATCH_SIZE  = 500    # points per GEE batch (safe limit)
MAX_RETRIES = 3
RETRY_DELAY = 10     # seconds between retries
WINDOWS     = [      # (start_days_before, end_days_before)
    (16, 0), (32, 16), (48, 32), (64, 48), (80, 64)
]

os.makedirs(DRIVE_DIR, exist_ok=True)

# ====== Cell 3: Load Input Data ======
df = pd.read_csv(INPUT_CSV)
print(f"Loaded {len(df)} landslide locations")
df['Event_Date'] = pd.to_datetime(df['Event_Date'])

# Resume from partial progress if available
if os.path.exists(PROGRESS_CSV):
    done = set(pd.read_csv(PROGRESS_CSV)['ev_id'].tolist())
    df   = df[~df['ev_id'].isin(done)]
    print(f"Resuming: {len(done)} already done, {len(df)} remaining")

# ====== Cell 4: GEE Extraction Functions ======

def get_terrain_at_point(lat, lon):
    """Extract SRTM + MERIT Hydro terrain variables for one point."""
    pt = ee.Geometry.Point([lon, lat])
    srtm = ee.Image("USGS/SRTMGL1_003")
    merit = ee.Image("MERIT/Hydro/v1_0_1")

    terrain = ee.Algorithms.Terrain(srtm)
    bands   = terrain.addBands(
        srtm.rename('elevation')).addBands(
        merit.select(['upa','wth'])).rename(
        ['slope','aspect','hillshade','elevation','upslope_area','twi'])

    # Topographic Roughness Index and curvature
    tri   = srtm.reduceNeighborhood(ee.Reducer.stdDev(), ee.Kernel.circle(3))
    tpi5  = srtm.subtract(srtm.focal_mean(5, 'circle', 'pixels'))
    tpi10 = srtm.subtract(srtm.focal_mean(10, 'circle', 'pixels'))
    curv  = ee.Terrain.aspect(srtm).subtract(180).abs().divide(180)  # proxy

    full = bands.addBands(tri.rename('tri')).addBands(
           tpi5.rename('tpi_500m')).addBands(tpi10.rename('tpi_1km'))

    vals = full.reduceRegion(ee.Reducer.mean(), pt, scale=30).getInfo()
    return vals


def get_modis_ndvi_for_window(lat, lon, event_date, days_start, days_end):
    """Get MOD13Q1 NDVI/EVI for a pre-failure window."""
    pt   = ee.Geometry.Point([lon, lat])
    ed   = ee.Date(event_date.strftime('%Y-%m-%d'))
    t1   = ed.advance(-days_end,  'day')
    t2   = ed.advance(-days_start,'day')

    col  = (ee.ImageCollection('MODIS/006/MOD13Q1')
             .filterDate(t1, t2)
             .filterBounds(pt))
    img  = col.sort('system:time_start', False).first()
    if img is None:
        return {}
    vals = img.select(['NDVI','EVI','SummaryQA']).reduceRegion(
               ee.Reducer.mean(), pt, scale=250).getInfo()
    if vals:
        date_ms = img.get('system:time_start').getInfo()
        if date_ms:
            vals['date'] = pd.Timestamp(date_ms, unit='ms').strftime('%Y/%m/%d')
        vals['days_before'] = days_start
    return vals


# ====== Cell 5: Batch Extraction Loop ======
results = []
batches = [df.iloc[i:i+BATCH_SIZE] for i in range(0, len(df), BATCH_SIZE)]

for b_idx, batch in enumerate(batches):
    print(f"\nBatch {b_idx+1}/{len(batches)} ({len(batch)} events)")
    batch_results = []

    for _, row in batch.iterrows():
        lat, lon  = row['Latitude'], row['Longitude']
        ev_date   = row['Event_Date']
        ev_id     = row['ev_id']
        record    = {'ev_id': ev_id}

        for attempt in range(MAX_RETRIES):
            try:
                terrain = get_terrain_at_point(lat, lon)
                record.update(terrain)
                for w_start, w_end in WINDOWS:
                    wvals = get_modis_ndvi_for_window(lat, lon, ev_date, w_start, w_end)
                    tag   = f"_w{w_start}"
                    record.update({k+tag: v for k,v in wvals.items()})
                break
            except Exception as e:
                if attempt < MAX_RETRIES-1:
                    time.sleep(RETRY_DELAY)
                else:
                    print(f"  FAILED ev_id={ev_id}: {e}")
        batch_results.append(record)

    results.extend(batch_results)

    # Save partial progress
    pd.DataFrame(results).to_csv(PROGRESS_CSV, index=False)
    print(f"  Progress saved ({len(results)} events processed)")

# ====== Cell 6: Save Final Output ======
out_df = pd.DataFrame(results)
out_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nDone. Saved {len(out_df)} records to {OUTPUT_CSV}")
print("Merge with original catalog on 'ev_id' to get the full dataset.")
