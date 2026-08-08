#!/usr/bin/env python3
"""
Sentinel-1 SAR Flood Inundation Mapping -- Permanent Standalone Web Map
=======================================================================
WHY THIS APPROACH
  geemap.addLayer() writes GEE tile-server URLs into the HTML.  Those URLs
  expire after ~2 days, causing layers to go blank.  This script avoids that
  by downloading every layer as a GeoTIFF, converting it to a palette PNG and
  publishing that PNG as page-local data -- either a base64 data URI inside the
  HTML or a file in ASSET_DIR next to it.  Either way nothing points at a
  Google server, so the map never expires.

METHOD
  Change-detection: VV backscatter difference (post - pre, dB).
  Flooded open water causes specular reflection away from the sensor,
  producing a sharp drop in VV backscatter.

DATA      Copernicus Sentinel-1 GRD IW (via Google Earth Engine)
AOI       Loaded from AOI_GEOJSON (merged townships polygon)
Pre-flood  2026-07-27  (+/-DATE_WINDOW days)
Post-flood 2026-08-08  (+/-DATE_WINDOW days)
Output     index.html + ASSET_DIR/*.png  -- self-hosted, never-expiring

DATA AVAILABILITY (checked against the GEE catalogue on 2026-08-08)
  This AOI is imaged by three Sentinel-1 tracks, each on a 12-day repeat:

      ASCENDING  track 143 : 2026-07-03, 07-15, 07-27, (08-08)
      ASCENDING  track  41 : 2026-07-08, 07-20, 08-01
      DESCENDING track 106 : 2026-07-12, 07-24, 08-05

  2026-07-27 and 2026-08-08 are consecutive passes of ASCENDING track 143 --
  identical viewing geometry, exactly 12 days apart, 100% AOI coverage.  That
  is the ideal matched pair for change detection, hence ORBIT_PASS/REL_ORBIT
  below.  Note this differs from the previous run of this script, which used
  DESCENDING track 106 for the 07-24/08-05 pair: the track follows the dates,
  it is not a property of the AOI.

  At the time of writing, the 2026-08-08 scene had not yet been ingested into
  Earth Engine (latest acquisition over this AOI was 2026-08-05).  GEE
  ingestion typically lags acquisition by a few days.  main() runs a pre-flight
  availability check and exits with a clear message -- and a table of what IS
  available -- if the post-flood scene is still missing.  Re-run unchanged once
  it lands.

RESOLUTION
  Every layer -- flood mask, permanent water AND both SAR VV backdrops -- is
  produced and displayed at Sentinel-1's native 10 m GRD pixel spacing.  The
  flood-area statistic is reduced at 10 m too.  Nothing is downsampled.

  That has a real cost, which is why the output is no longer a single file.
  At 10 m this AOI is roughly 9,500 x 7,900 px (75 Mpx) per layer:

      flood / water mask @10 m ->   ~1 MB PNG  (binary, compresses ~100x)
      SAR VV backdrop    @10 m ->  ~56 MB PNG  (speckle, near-incompressible)

  Base64 inflates by 1.33x, so embedding two native SAR backdrops would make a
  ~150 MB HTML file -- over GitHub's hard 100 MB per-file push limit, and slow
  for a browser to even parse.  So each layer is:

    * cut into chunks of at most MAX_CHUNK_PX pixels, with all-nodata chunks
      dropped (the AOI corners), so no single image is monstrous and the
      browser can free memory per chunk; and
    * embedded as base64 if the whole layer is under EMBED_MAX_MB, otherwise
      written to ASSET_DIR/ and referenced by relative path.

  Chunks of a layer live in one Leaflet FeatureGroup, so they toggle as a unit.
  The SAR layers start hidden, so their bytes are only fetched if the user
  actually turns them on.

  Deploying therefore means publishing index.html *and* the ASSET_DIR folder.

REQUIREMENTS
    pip install earthengine-api geemap folium rasterio numpy Pillow

FIRST-TIME GEE SETUP
    1. Sign up at https://earthengine.google.com
    2. Enable the Earth Engine API in a Google Cloud project
    3. Run:  earthengine authenticate
       (or let this script call ee.Authenticate() automatically)
"""

import os
import json
import math
import base64
import shutil
import glob
from io import BytesIO

from rasterio.merge import merge as rio_merge
from rasterio.windows import Window

import numpy as np
import ee
import geemap
import folium
from folium.plugins import MiniMap, Fullscreen, MousePosition
import rasterio
from PIL import Image

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

# AOI is read from a GeoJSON file (FeatureCollection, Feature or bare geometry).
AOI_GEOJSON = "merged_townships.geojson"

PRE_FLOOD_DATE  = "2026-07-27"
POST_FLOOD_DATE = "2026-08-08"
DATE_WINDOW     = 1       # +/- days to search for imagery around each target date

FLOOD_DB_THRESH = -3.0    # VV drop (dB) below which pixel is classified flooded
SPECKLE_RADIUS  = 30      # focal-mean kernel radius in metres

# Sentinel-1 orbit pass used for BOTH pre- and post-flood composites.
# Change detection is only valid when both images share the same viewing
# geometry (same pass, ideally same relative orbit).  Mixing ASCENDING and
# DESCENDING scenes introduces a systematic backscatter bias that destroys the
# flood signal.  For this AOI/date pair, ASCENDING track 143 images both
# 2026-07-27 and 2026-08-08 (12 days apart) at 100% coverage.
ORBIT_PASS = "ASCENDING"

# Optional: pin a single Sentinel-1 relative orbit (ground track) so pre and
# post use the identical acquisition path.  Set to None to use every track of
# the chosen pass.  Track 143 is the only track covering this AOI on both of
# these dates -- descending track 106 (used for the earlier 07-24/08-05 pair)
# images this AOI on 07-24 and 08-05, neither of which is wanted here.
REL_ORBIT = 143

# Your GEE Cloud project ID -- leave "" to let GEE auto-detect from credentials
GEE_PROJECT = "gee-python-419405"

# Native Sentinel-1 GRD IW pixel spacing.  EVERY layer and the flood-area
# statistic use this -- there is no downsampling anywhere in this script.
EXPORT_SCALE = 10

# SAR VV backdrops are rendered at native resolution too (see RESOLUTION above).
# Kept as a separate name so the heavy display layers can be dropped to a
# coarser scale for a quick test run without touching the flood product.
SAR_DISPLAY_SCALE = 10

# Written to the repo root as index.html so GitHub Pages serves it directly.
# Never hand-edit index.html -- re-run this script to regenerate it.
OUTPUT_HTML = "index.html"
ASSET_DIR   = "assets"        # page-local PNGs for layers too big to embed
TEMP_DIR    = "_gee_tmp"      # deleted automatically after the map is built

# Target pixels per download tile.  GEE's getDownloadURL rejects requests over
# ~48 MB; every layer is fetched as uint8 (1 byte/px), so ~15 M px per tile
# leaves generous headroom.  Tile size in degrees is derived from this.
TILE_TARGET_PX = 15_000_000

# Max pixels in one displayed PNG chunk.  A 75 Mpx single image decodes to
# ~300 MB of RGBA in the browser and is a hard failure on mobile Safari;
# ~16 Mpx chunks decode to ~64 MB each and can be freed independently.
MAX_CHUNK_PX = 16_000_000

# Layers whose total PNG payload exceeds this are written to ASSET_DIR instead
# of being base64-embedded in the HTML.  Keeps index.html small enough to open
# instantly and to push to GitHub (hard 100 MB per-file limit).
EMBED_MAX_MB = 8.0

# The legend is collapsible so it does not obscure the map.  True = start as a
# small header pill (headline figure still visible, click to expand).
LEGEND_START_COLLAPSED = True

DEG_PER_M = 1.0 / 111_320.0   # approx degrees of latitude per metre

# ──────────────────────────────────────────────────────────────────────────────
# 0. AOI LOADING
# ──────────────────────────────────────────────────────────────────────────────

def load_aoi(path: str) -> dict:
    """
    Read an AOI from GeoJSON and return everything the rest of the script needs.

    Accepts a FeatureCollection, a single Feature, or a bare geometry so the
    AOI file can be swapped without touching code.
    """
    with open(path, encoding="utf-8") as fh:
        gj = json.load(fh)

    if gj.get("type") == "FeatureCollection":
        features = gj["features"]
    elif gj.get("type") == "Feature":
        features = [gj]
    else:                                    # bare geometry
        features = [{"type": "Feature", "geometry": gj, "properties": {}}]

    if not features:
        raise RuntimeError(f"{path} contains no features.")

    geoms = [f["geometry"] for f in features]
    props = features[0].get("properties", {}) or {}

    # Union multiple features into one analysis geometry.
    ee_geom = ee.Geometry(geoms[0])
    for g in geoms[1:]:
        ee_geom = ee_geom.union(ee.Geometry(g), maxError=1)

    # Bounding box over every coordinate, at any nesting depth.
    xs, ys = [], []

    def walk(coords):
        if coords and isinstance(coords[0], (int, float)):
            xs.append(coords[0])
            ys.append(coords[1])
        else:
            for c in coords:
                walk(c)

    for g in geoms:
        walk(g["coordinates"])

    bbox = {"west": min(xs), "east": max(xs),
            "south": min(ys), "north": max(ys)}

    # Human-readable naming, tolerant of whatever the file happens to carry.
    name = (props.get("Name_EN") or props.get("name")
            or os.path.splitext(os.path.basename(path))[0])
    townships = [t.get("TS") for t in props.get("Townships", []) if t.get("TS")]

    return {
        "ee_geom":    ee_geom,
        "geojson":    {"type": "FeatureCollection", "features": features},
        "bbox":       bbox,
        "name":       name,
        "townships":  townships,
        "name_mmr":   props.get("Name_MMR", ""),
    }

# ──────────────────────────────────────────────────────────────────────────────
# 1. GEE INITIALISATION
# ──────────────────────────────────────────────────────────────────────────────

def init_gee() -> None:
    kwargs = {"project": GEE_PROJECT} if GEE_PROJECT else {}
    try:
        ee.Initialize(**kwargs)
        print("[OK] GEE initialised.")
    except ee.EEException:
        print("[!] Not authenticated -- running ee.Authenticate() ...")
        ee.Authenticate()
        ee.Initialize(**kwargs)
        print("[OK] GEE initialised.")

# ──────────────────────────────────────────────────────────────────────────────
# 2. SENTINEL-1 RETRIEVAL
# ──────────────────────────────────────────────────────────────────────────────

def _date_range(center_date: str, window: int) -> tuple:
    """Inclusive +/-window-day range as (start, end) for ee filterDate.

    filterDate's end is EXCLUSIVE, so the end bound is advanced one extra day;
    otherwise '+/-1 day' would silently mean 'the day before and the day of'.
    """
    return (ee.Date(center_date).advance(-window, "day"),
            ee.Date(center_date).advance(window + 1, "day"))


def load_s1(aoi: ee.Geometry, center_date: str, window: int,
            orbit_pass: str = ORBIT_PASS, rel_orbit=REL_ORBIT) -> ee.Image:
    """Mean composite of Sentinel-1 IW GRD scenes within +/-window days.

    Both pre- and post-flood composites must share the same viewing geometry.
    orbit_pass fixes ASCENDING vs DESCENDING; rel_orbit (if not None) further
    pins one relative orbit / ground track so the pair uses the identical
    acquisition path -- the cleanest possible basis for change detection.
    """
    d_start, d_end = _date_range(center_date, window)

    def build(with_pass: bool) -> ee.ImageCollection:
        c = (
            ee.ImageCollection("COPERNICUS/S1_GRD")
            .filterBounds(aoi)
            .filterDate(d_start, d_end)
            .filter(ee.Filter.eq("instrumentMode", "IW"))
            .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
            .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
        )
        if with_pass:
            c = c.filter(ee.Filter.eq("orbitProperties_pass", orbit_pass))
        if rel_orbit is not None:
            c = c.filter(ee.Filter.eq("relativeOrbitNumber_start", rel_orbit))
        return c.select(["VV", "VH"])

    col = build(with_pass=True)
    n = col.size().getInfo()

    if n == 0:
        print(f"  [!] WARNING: no {orbit_pass} scenes within +/-{window}d of "
              f"{center_date}. Falling back to BOTH passes -- this mixes viewing "
              f"geometries and can bias the pre/post VV difference. Consider "
              f"widening DATE_WINDOW so a matching {orbit_pass} scene is found.")
        col = build(with_pass=False)
        n = col.size().getInfo()

    if n == 0:
        raise RuntimeError(
            f"No Sentinel-1 data within +/-{window} days of {center_date} "
            f"(pass={orbit_pass}, rel_orbit={rel_orbit}). "
            "Try widening DATE_WINDOW or setting REL_ORBIT = None."
        )

    track = f", track {rel_orbit}" if rel_orbit is not None else ""
    print(f"  [{center_date}] {n} scene(s) found{track} -> mean composite.")
    return col.mean().clip(aoi)

# ──────────────────────────────────────────────────────────────────────────────
# 2b. PRE-FLIGHT DATA AVAILABILITY CHECK
# ──────────────────────────────────────────────────────────────────────────────

def s1_availability(aoi: ee.Geometry, center_date: str, window: int,
                    orbit_pass: str, rel_orbit) -> tuple:
    """Scene count and AOI coverage fraction for one date on the chosen track."""
    d0, d1 = _date_range(center_date, window)
    col = (
        ee.ImageCollection("COPERNICUS/S1_GRD")
        .filterBounds(aoi).filterDate(d0, d1)
        .filter(ee.Filter.eq("instrumentMode", "IW"))
        .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
        .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
        .filter(ee.Filter.eq("orbitProperties_pass", orbit_pass))
    )
    if rel_orbit is not None:
        col = col.filter(ee.Filter.eq("relativeOrbitNumber_start", rel_orbit))
    col = col.select(["VV"])

    n = col.size().getInfo()
    if n == 0:
        return 0, 0.0
    frac = col.mosaic().mask().reduceRegion(
        reducer=ee.Reducer.mean(), geometry=aoi, scale=100, maxPixels=1e10
    ).get("VV")
    return n, (ee.Number(frac).getInfo() or 0.0)


def recent_acquisitions(aoi: ee.Geometry, days_back: int = 45) -> list:
    """
    Every S1 IW acquisition over the AOI in the last days_back days, as a list
    of (date, pass, track) -- one row per unique combination.

    Used only on the failure path, to tell the user which date/track pairs they
    could actually use instead of just saying "not available".
    """
    end   = ee.Date(POST_FLOOD_DATE).advance(3, "day")
    start = end.advance(-days_back, "day")
    col = (
        ee.ImageCollection("COPERNICUS/S1_GRD")
        .filterBounds(aoi).filterDate(start, end)
        .filter(ee.Filter.eq("instrumentMode", "IW"))
        .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
        .sort("system:time_start")
    )
    if col.size().getInfo() == 0:
        return []

    dated = col.map(lambda im: im.set(
        "_d", ee.Date(im.get("system:time_start")).format("YYYY-MM-dd")))
    info = ee.Dictionary({
        "d": dated.aggregate_array("_d"),
        "p": col.aggregate_array("orbitProperties_pass"),
        "t": col.aggregate_array("relativeOrbitNumber_start"),
    }).getInfo()

    seen, rows = set(), []
    for d, p, t in zip(info["d"], info["p"], info["t"]):
        key = (d, p, int(t))
        if key not in seen:
            seen.add(key)
            rows.append(key)
    return sorted(rows)


def preflight(aoi: ee.Geometry) -> None:
    """
    Verify both dates have imagery on the configured track before doing any
    heavy work.  Sentinel-1 reaches the Earth Engine catalogue a few days after
    acquisition, so a very recent post-flood date can simply not be there yet.
    Failing here -- with the real catalogue contents printed -- beats crashing
    halfway through a multi-gigabyte download run.
    """
    print(f"\n> Pre-flight: checking {ORBIT_PASS} track {REL_ORBIT} availability ...")
    missing = []
    for tag, date in (("pre ", PRE_FLOOD_DATE), ("post", POST_FLOOD_DATE)):
        n, cov = s1_availability(aoi, date, DATE_WINDOW, ORBIT_PASS, REL_ORBIT)
        state = f"{n} scene(s), {cov:.0%} of AOI" if n else "NO DATA"
        print(f"    {tag} {date} (+/-{DATE_WINDOW}d): {state}")
        if n == 0:
            missing.append((tag.strip(), date))
        elif cov < 0.5:
            print(f"        [!] only {cov:.0%} coverage -- most of the AOI "
                  f"cannot be assessed on this date.")

    if missing:
        rows = recent_acquisitions(aoi)
        table = "\n".join(
            f"         {d}   {p:<10s} track {t}"
            + ("   <- configured track" if (t == REL_ORBIT and p == ORBIT_PASS) else "")
            for d, p, t in rows
        ) or "         (nothing found)"
        names = ", ".join(f"{t}-flood {d}" for t, d in missing)
        raise SystemExit(
            f"\n[STOP] No Sentinel-1 imagery for: {names}\n"
            f"       (pass={ORBIT_PASS}, track={REL_ORBIT}, "
            f"window=+/-{DATE_WINDOW}d)\n\n"
            f"       Sentinel-1 usually appears in GEE a few days after it is\n"
            f"       acquired. If the missing date is recent, the scene is most\n"
            f"       likely still being ingested -- re-run this script unchanged\n"
            f"       in a day or two. Nothing else needs editing.\n\n"
            f"       Acquisitions over this AOI currently in the catalogue:\n"
            f"{table}\n\n"
            f"       Each track repeats every 12 days. To map a different pair,\n"
            f"       pick two dates from the SAME track above and set\n"
            f"       PRE_FLOOD_DATE / POST_FLOOD_DATE / ORBIT_PASS / REL_ORBIT\n"
            f"       to match -- never mix tracks or passes.\n"
        )
    print("    [OK] both dates available.")

# ──────────────────────────────────────────────────────────────────────────────
# 3. SPECKLE FILTERING
# ──────────────────────────────────────────────────────────────────────────────

def speckle_filter(image: ee.Image, radius: int = SPECKLE_RADIUS) -> ee.Image:
    """Boxcar (focal-mean) speckle suppression."""
    return image.focal_mean(radius=radius, kernelType="circle", units="meters")

# ──────────────────────────────────────────────────────────────────────────────
# 4. FLOOD DETECTION
# ──────────────────────────────────────────────────────────────────────────────

def detect_floods(pre: ee.Image, post: ee.Image,
                  aoi: ee.Geometry,
                  thresh: float = FLOOD_DB_THRESH):
    """
    Returns
    -------
    flood_mask  binary (1 = newly flooded)
    perm_water  binary (1 = JRC permanent/seasonal water)
    diff_db     continuous VV change image (dB)
    """
    diff_db = post.select("VV").subtract(pre.select("VV")).rename("VV_change")

    # JRC Global Surface Water v1.4 -- seasonality >= 4 months/year.
    #
    # unmask(0) is essential: the seasonality band carries NO DATA over land
    # that has never been observed as water.  Because ee .And() intersects
    # masks, feeding the raw (masked) image into .And(perm_water.Not()) deletes
    # every never-been-water pixel from the flood mask -- i.e. exactly the dry
    # land a flood would newly cover.  Unmasking to 0 makes it a filter ("this
    # pixel is not permanent water") instead of a data-availability constraint.
    perm_water = (
        ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
        .select("seasonality").gte(4).unmask(0).clip(aoi)
    )

    # HydroSHEDS terrain slope -- exclude steep pixels (unlikely to flood).
    # Left masked on purpose: the DEM has no data over sea, so its mask doubles
    # as a land mask and keeps open ocean out of the result.
    slope = ee.Terrain.slope(ee.Image("WWF/HydroSHEDS/03VFDEM")).clip(aoi)

    flood_mask = (
        diff_db.lt(thresh)          # large VV decrease
        .And(perm_water.Not())      # not already open water
        .And(slope.lt(5))           # flat terrain only
        .rename("flood")
    )
    return flood_mask, perm_water, diff_db

# ──────────────────────────────────────────────────────────────────────────────
# 5. STATISTICS
# ──────────────────────────────────────────────────────────────────────────────

def compute_flood_area_km2(flood_mask: ee.Image, aoi: ee.Geometry,
                           scale: int = EXPORT_SCALE) -> float:
    """Flood area, reduced at native resolution so nothing is generalised away."""
    area = (
        flood_mask
        .multiply(ee.Image.pixelArea())
        .reduceRegion(
            reducer=ee.Reducer.sum(), geometry=aoi,
            scale=scale, maxPixels=1e11,
        )
    )
    return ee.Number(area.get("flood")).divide(1e6).getInfo()


def compute_valid_fraction(flood_mask: ee.Image, aoi: ee.Geometry) -> float:
    """
    Fraction of the AOI that actually had usable pre AND post data.

    The VV difference is masked wherever either composite is missing, so a
    partially-imaged date silently shrinks the analysed area.  Reporting this
    keeps an unimaged gap from being read as 'no flooding detected here'.
    Reduced at a coarser scale -- this is a QA percentage, not a measurement.
    """
    frac = (
        flood_mask.mask()
        .reduceRegion(
            reducer=ee.Reducer.mean(), geometry=aoi,
            scale=60, maxPixels=1e10,
        )
    )
    return ee.Number(frac.get("flood")).getInfo() or 0.0

# ──────────────────────────────────────────────────────────────────────────────
# 6. DOWNLOAD GEE IMAGE -> LOCAL GEOTIFF
# ──────────────────────────────────────────────────────────────────────────────

def tile_deg_for(scale: int) -> float:
    """Tile edge in degrees that keeps one uint8 tile near TILE_TARGET_PX."""
    return math.sqrt(TILE_TARGET_PX) * scale * DEG_PER_M


def download_as_geotiff(image: ee.Image, name: str, aoi: ee.Geometry,
                        scale: int, bbox: dict) -> str:
    """
    Download a uint8 GEE image to TEMP_DIR/<name> (EPSG:4326).

    GEE's getDownloadURL API rejects requests larger than ~48 MB, which at 10 m
    is only a fraction of a degree.  The BBOX is therefore split into tiles
    sized by tile_deg_for(), downloaded separately, then merged with
    rasterio.merge -- so any resolution works, including native 10 m.
    """
    os.makedirs(TEMP_DIR, exist_ok=True)
    final_path = os.path.join(TEMP_DIR, name)

    if os.path.exists(final_path):
        print(f"  Using cached {name}.")
        return final_path

    tile_deg = tile_deg_for(scale)

    tile_defs, tile_paths = [], []
    lat = bbox["south"]
    while lat < bbox["north"]:
        lat_end = min(lat + tile_deg, bbox["north"])
        lon = bbox["west"]
        while lon < bbox["east"]:
            lon_end = min(lon + tile_deg, bbox["east"])
            tile_defs.append((lon, lat, lon_end, lat_end))
            lon = lon_end
        lat = lat_end

    n_tiles = len(tile_defs)
    print(f"  Downloading {name}: {n_tiles} tile(s) at {scale} m/px "
          f"({tile_deg:.3f} deg tiles) ...")

    for i, (w, s, e, n) in enumerate(tile_defs, 1):
        tile_region = ee.Geometry.Rectangle([w, s, e, n])
        tile_path   = os.path.join(TEMP_DIR, f"_tile_{i:03d}_{name}")
        tile_paths.append(tile_path)

        if not os.path.exists(tile_path):
            print(f"    Tile {i}/{n_tiles} ...")
            geemap.ee_export_image(
                image.clip(tile_region),
                filename=tile_path,
                scale=scale,
                region=tile_region,
                crs="EPSG:4326",
                file_per_band=False,
            )
            if not os.path.exists(tile_path):
                raise RuntimeError(
                    f"Tile {i}/{n_tiles} download failed for {name}.\n"
                    f"Try lowering TILE_TARGET_PX (currently {TILE_TARGET_PX})."
                )

    if n_tiles == 1:
        shutil.move(tile_paths[0], final_path)
    else:
        print(f"  Merging {n_tiles} tiles -> {name} ...")
        datasets = [rasterio.open(p) for p in tile_paths]
        mosaic, transform = rio_merge(datasets)
        profile = datasets[0].profile.copy()
        profile.update({
            "height":    mosaic.shape[1],
            "width":     mosaic.shape[2],
            "transform": transform,
        })
        for ds in datasets:
            ds.close()
        with rasterio.open(final_path, "w", **profile) as dst:
            dst.write(mosaic)
        for p in tile_paths:
            if os.path.exists(p):
                os.remove(p)

    print(f"  [OK] {name} saved.")
    return final_path

# ──────────────────────────────────────────────────────────────────────────────
# 7. GEOTIFF -> PALETTE PNG CHUNKS  (page-local -- never expires)
# ──────────────────────────────────────────────────────────────────────────────
#
# Palette ("P" mode) PNGs are used rather than RGBA.  At 10 m this AOI is
# ~75 Mpx: RGBA would allocate 300 MB per layer in memory and compress far
# worse, while a palette image is 1 byte/px and lets PNG's filters do their job.
# Index 0 is reserved for "no data / nothing here" and made transparent.

def _palette_flood() -> list:
    pal = [0, 0, 0, 215, 48, 39]          # 0 = transparent, 1 = red
    return pal + [0, 0, 0] * 254


def _palette_water() -> list:
    pal = [0, 0, 0, 33, 102, 172]         # 0 = transparent, 1 = blue
    return pal + [0, 0, 0] * 254


def _palette_gray() -> list:
    pal = []
    for i in range(256):                  # 0 = transparent, 1..255 = grays
        pal += [i, i, i]
    return pal


def _encode_png(arr: np.ndarray, palette: list) -> bytes:
    """uint8 array -> palette PNG bytes, index 0 transparent.

    frombytes("P", ...) treats the buffer as raw palette indices, so values
    survive exactly.  Going via fromarray().convert("P") would invite PIL to
    re-quantise, which for the gray ramp would silently shift DN values.
    """
    h, w = arr.shape
    img = Image.frombytes("P", (w, h), np.ascontiguousarray(arr).tobytes())
    img.putpalette(palette)
    buf = BytesIO()
    img.save(buf, format="PNG", optimize=True, transparency=0)
    return buf.getvalue()


def _chunk_grid(width: int, height: int, max_px: int = None) -> tuple:
    """Rows/cols needed so every chunk stays under max_px pixels."""
    max_dim = int(math.sqrt(MAX_CHUNK_PX if max_px is None else max_px))
    return (max(1, math.ceil(height / max_dim)),
            max(1, math.ceil(width  / max_dim)))


def clear_assets() -> None:
    """
    Drop PNGs written by a previous run so stale chunks can never be served.

    Only files matching the names this script generates are removed -- anything
    else a user has put in ASSET_DIR is left alone.
    """
    if not os.path.isdir(ASSET_DIR):
        return
    stale = glob.glob(os.path.join(ASSET_DIR, "*_[0-9][0-9]*.png"))
    for p in stale:
        os.remove(p)
    if stale:
        print(f"  Removed {len(stale)} PNG(s) from a previous run in {ASSET_DIR}/.")


def geotiff_to_layer(path: str, palette: list, key: str) -> dict:
    """
    Read a uint8 GeoTIFF and turn it into browser-ready PNG chunks.

    Returns {"pieces": [(image_ref, [[s,w],[n,e]]), ...],
             "external": bool, "mb": float}

    image_ref is either a base64 data URI (small layers, embedded in the HTML)
    or a relative path into ASSET_DIR (big layers).  Neither points at a GEE
    tile server, so neither expires.  Chunks that are entirely nodata -- the
    bbox corners outside the AOI polygon -- are dropped rather than shipped as
    megabytes of transparent pixels.
    """
    with rasterio.open(path) as src:
        W, H = src.width, src.height
        b    = src.bounds                    # EPSG:4326 from the download step
        xres = (b.right - b.left) / W
        yres = (b.top   - b.bottom) / H

        n_rows, n_cols = _chunk_grid(W, H)
        r_step = math.ceil(H / n_rows)
        c_step = math.ceil(W / n_cols)

        pieces, skipped = [], 0
        for r in range(n_rows):
            r0, r1 = r * r_step, min(H, (r + 1) * r_step)
            if r0 >= r1:
                continue
            for c in range(n_cols):
                c0, c1 = c * c_step, min(W, (c + 1) * c_step)
                if c0 >= c1:
                    continue

                arr = src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0))
                if arr.dtype != np.uint8:
                    arr = np.clip(arr, 0, 255).astype(np.uint8)
                if not arr.any():            # all nodata -- nothing to draw
                    skipped += 1
                    continue

                bounds = [[b.top - r1 * yres, b.left + c0 * xres],
                          [b.top - r0 * yres, b.left + c1 * xres]]
                pieces.append((bounds, _encode_png(arr, palette)))

    total_mb = sum(len(p) for _, p in pieces) / 1e6
    external = total_mb > EMBED_MAX_MB

    refs = []
    for i, (bounds, png) in enumerate(pieces, 1):
        if external:
            os.makedirs(ASSET_DIR, exist_ok=True)
            fname = f"{key}_{i:02d}.png"
            with open(os.path.join(ASSET_DIR, fname), "wb") as fh:
                fh.write(png)
            ref = f"{ASSET_DIR}/{fname}"     # forward slash: this is a URL
        else:
            ref = "data:image/png;base64," + base64.b64encode(png).decode()
        refs.append((ref, bounds))

    where = f"{ASSET_DIR}/" if external else "embedded"
    print(f"    {key:9s} {W}x{H} px -> {len(refs)} chunk(s), "
          f"{total_mb:.1f} MB, {where}"
          + (f" ({skipped} empty chunk(s) dropped)" if skipped else ""))

    return {"pieces": refs, "external": external, "mb": total_mb}

# ──────────────────────────────────────────────────────────────────────────────
# 8. BUILD PERMANENT STANDALONE FOLIUM MAP
# ──────────────────────────────────────────────────────────────────────────────

# 1x1 transparent PNG.  Handed to folium purely so it does not touch the disk;
# the real reference is assigned to .url immediately afterwards (see add_layer).
_BLANK_PNG = ("data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAf"
              "FcSJAAAAC0lEQVR42mNkYAAAAAYAAjCB0C8AAAAASUVORK5CYII=")


def add_layer(m: folium.Map, layer: dict, name: str,
              opacity: float, show: bool, smooth: bool = False) -> None:
    """
    Add one raster layer's chunks to the map as a single toggleable group.

    A FeatureGroup is used even for one chunk so LayerControl always shows one
    entry per layer, and Leaflet only builds the <img> elements when the group
    is switched on -- which is what keeps the hidden 10 m SAR backdrops from
    being fetched unless someone asks for them.

    .url is assigned after construction rather than passing the reference as
    `image=`.  folium's image_to_url() treats any string that is not a URL as a
    file path, opens it and base64-inlines the bytes -- which would silently
    drag every externalised PNG back into the HTML and produce the ~150 MB page
    that ASSET_DIR exists to avoid.  A relative path has no URL scheme, so the
    only way to emit it verbatim is to set it directly.

    smooth=True tags the chunks so CSS can turn off nearest-neighbour scaling.
    folium hard-codes `image-rendering:pixelated` on every image layer, which
    is right for a binary mask but makes 10 m speckle alias badly when the map
    is zoomed out past native resolution.
    """
    group = folium.FeatureGroup(name=name, show=show)
    extra = {"className": "fm-smooth"} if smooth else {}
    for ref, bounds in layer["pieces"]:
        overlay = folium.raster_layers.ImageOverlay(
            image=_BLANK_PNG, bounds=bounds, opacity=opacity, **extra,
        )
        overlay.url = ref
        overlay.add_to(group)
    group.add_to(m)


def build_map(overlays: dict, stats: dict, aoi_info: dict, output: str) -> None:
    """
    Build a folium map where every raster layer is served by the page itself --
    base64 data URIs for small layers, sibling PNGs for large ones.  No GEE
    tile links anywhere, so the map does not expire.
    """
    bbox = aoi_info["bbox"]
    cx = (bbox["south"] + bbox["north"]) / 2
    cy = (bbox["west"]  + bbox["east"])  / 2

    m = folium.Map(location=[cx, cy], zoom_start=10, tiles=None)

    # ── Basemaps ──────────────────────────────────────────────────────────────
    folium.TileLayer(
        tiles=("https://server.arcgisonline.com/ArcGIS/rest/services/"
               "World_Imagery/MapServer/tile/{z}/{y}/{x}"),
        attr="Esri World Imagery",
        name="Satellite (Esri)",
        overlay=False, control=True,
    ).add_to(m)
    folium.TileLayer(
        "OpenStreetMap", name="OpenStreetMap",
        overlay=False, control=True,
    ).add_to(m)

    # ── Raster overlays (page-local -- permanent) ─────────────────────────────
    # SAR backdrops start hidden: at 10 m they are the heavy layers, and
    # Leaflet does not fetch a group's images until it is switched on.
    add_layer(m, overlays["pre_sar"],
              f"Pre-flood VV SAR ({PRE_FLOOD_DATE}, {SAR_DISPLAY_SCALE} m)",
              opacity=0.85, show=False, smooth=True)
    add_layer(m, overlays["post_sar"],
              f"Post-flood VV SAR ({POST_FLOOD_DATE}, {SAR_DISPLAY_SCALE} m)",
              opacity=0.85, show=False, smooth=True)
    add_layer(m, overlays["water"],
              f"Permanent Water (JRC, {EXPORT_SCALE} m)",
              opacity=0.75, show=True)
    add_layer(m, overlays["flood"],
              f"Flood Inundation ({EXPORT_SCALE} m)",
              opacity=0.85, show=True)

    # ── AOI boundary (GeoJSON literal embedded in HTML) ───────────────────────
    folium.GeoJson(
        aoi_info["geojson"],
        name="AOI Boundary",
        style_function=lambda _: {"color": "yellow", "weight": 2,
                                  "fillOpacity": 0},
        tooltip=aoi_info["name"],
    ).add_to(m)

    # ── Controls ──────────────────────────────────────────────────────────────
    folium.LayerControl(collapsed=False).add_to(m)
    Fullscreen(position="topright").add_to(m)
    MiniMap(toggle_display=True).add_to(m)
    MousePosition(position="bottomleft", separator=" | ",
                  prefix="Lat / Lon:").add_to(m)

    # ── Dashboard panel ───────────────────────────────────────────────────────
    flood_km2  = stats["flood_km2"]
    aoi_km2    = stats["aoi_km2"]
    valid_frac = stats["valid_frac"]
    pct_aoi    = (flood_km2 / aoi_km2 * 100.0) if aoi_km2 else 0.0

    if valid_frac < 0.995:
        coverage_note = (
            f'<div style="color:#b35806;margin-top:3px">'
            f'&#9888; Analysed area: <b>{valid_frac:.0%}</b> of AOI &mdash; the rest is '
            f'sea, steep terrain, or had no Sentinel-1 data on one of the two '
            f'dates, so it is <b>not</b> assessed (blank &ne; flood-free).</div>'
        )
    else:
        coverage_note = (
            f'<div class="fm-dim" style="margin-top:3px">'
            f'Analysed area: {valid_frac:.0%} of AOI</div>'
        )

    ts_rows = ", ".join(aoi_info["townships"]) if aoi_info["townships"] else "&mdash;"

    sar_mb = overlays["pre_sar"]["mb"] + overlays["post_sar"]["mb"]
    sar_note = (
        f'<div class="fm-dim" style="margin-top:3px">'
        f'SAR backdrops are off by default &mdash; {sar_mb:,.0f} MB at 10 m, '
        f'loaded only when switched on.</div>'
    )

    # CSS/JS kept out of the f-string below so the braces need no escaping.
    legend_css_js = """
    <style>
      /* Two classes, so this beats folium's blanket .leaflet-image-layer
         { image-rendering: pixelated } and lets the SAR speckle resample
         smoothly when the map is zoomed out past 10 m. */
      .leaflet-image-layer.fm-smooth {
        image-rendering: auto; image-rendering: smooth;
      }
      #fm-legend {
        position:fixed; bottom:30px; right:10px; z-index:1000;
        background:rgba(255,255,255,0.94);
        border:1px solid #b8b8b8; border-radius:8px;
        box-shadow:0 2px 8px rgba(0,0,0,0.22);
        font-family:'Segoe UI',Arial,sans-serif;
        font-size:12px; line-height:1.55;
        max-width:262px; overflow:hidden;
      }
      #fm-legend-hdr {
        display:flex; align-items:center; gap:7px;
        padding:6px 10px; cursor:pointer; user-select:none; white-space:nowrap;
      }
      #fm-legend-hdr:hover { background:rgba(0,0,0,0.05); }
      #fm-caret { margin-left:auto; color:#888; font-size:11px; }
      #fm-legend-body {
        padding:8px 11px 10px; border-top:1px solid #e4e4e4;
        max-height:58vh; overflow-y:auto;
      }
      #fm-legend.fm-collapsed #fm-legend-body { display:none; }
      #fm-legend .fm-sw {
        display:inline-block; width:12px; height:12px; border-radius:2px;
        vertical-align:middle; margin-right:6px;
      }
      #fm-legend .fm-sep { height:1px; background:#e4e4e4; margin:6px 0; }
      #fm-legend .fm-dim { color:#999; }
      #fm-legend .fm-kpi {
        display:flex; justify-content:space-between; gap:10px; margin:1px 0;
      }
      #fm-legend .fm-kpi b { white-space:nowrap; }
    </style>
    <script>
      function fmToggleLegend() {
        var el = document.getElementById('fm-legend');
        el.classList.toggle('fm-collapsed');
        document.getElementById('fm-caret').textContent =
          String.fromCharCode(el.classList.contains('fm-collapsed') ? 0x25B8 : 0x25BE);
      }
    </script>
    """

    collapsed_cls = "fm-collapsed" if LEGEND_START_COLLAPSED else ""
    caret_glyph   = "&#9656;" if LEGEND_START_COLLAPSED else "&#9662;"
    track_txt     = f"track {REL_ORBIT}" if REL_ORBIT is not None else "all tracks"

    legend_html = f"""
    <div id="fm-legend" class="{collapsed_cls}">
      <div id="fm-legend-hdr" onclick="fmToggleLegend()"
           title="Click to show / hide the legend">
        <span class="fm-sw" style="background:#d73027;opacity:0.85"></span>
        <b>{flood_km2:,.1f} km&sup2; flooded</b>
        <span id="fm-caret">{caret_glyph}</span>
      </div>
      <div id="fm-legend-body">
        <b>Sentinel-1 Flood Map</b><br>
        <span class="fm-dim">{aoi_info["name"]}</span>
        <div class="fm-sep"></div>
        <div class="fm-kpi"><span>Pre-flood</span><b>{PRE_FLOOD_DATE}</b></div>
        <div class="fm-kpi"><span>Post-flood</span><b>{POST_FLOOD_DATE}</b></div>
        <div class="fm-kpi"><span>Orbit</span><b>{ORBIT_PASS.title()}, {track_txt}</b></div>
        <div class="fm-sep"></div>
        <span class="fm-sw" style="background:#d73027;opacity:0.85"></span>Flood Inundation<br>
        <span class="fm-sw" style="background:#2166ac;opacity:0.75"></span>Permanent Water (JRC)
        <div class="fm-sep"></div>
        <div class="fm-kpi"><span>Flood extent</span><b>{flood_km2:,.1f} km&sup2;</b></div>
        <div class="fm-kpi"><span>Share of AOI</span><b>{pct_aoi:.1f} %</b></div>
        <div class="fm-kpi"><span>AOI area</span><b>{aoi_km2:,.0f} km&sup2;</b></div>
        {coverage_note}
        <div class="fm-sep"></div>
        <div class="fm-kpi"><span>Resolution</span><b>{EXPORT_SCALE} m native</b></div>
        <div class="fm-kpi"><span>Applies to</span><b>every layer</b></div>
        {sar_note}
        <div class="fm-sep"></div>
        <span class="fm-dim">
          Townships: {ts_rows}<br>
          Threshold : VV &lt; {FLOOD_DB_THRESH} dB<br>
          Speckle filter : {SPECKLE_RADIUS} m radius<br>
          Permanent water excluded (JRC)
        </span>
      </div>
    </div>
    """
    m.get_root().html.add_child(folium.Element(legend_css_js + legend_html))

    m.save(output)
    print(f"[OK] Map saved -> {output}")

# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    init_gee()

    aoi_info = load_aoi(AOI_GEOJSON)
    aoi      = aoi_info["ee_geom"]
    bbox     = aoi_info["bbox"]
    aoi_km2  = aoi.area(maxError=10).divide(1e6).getInfo()

    print(f"\n> AOI: {aoi_info['name']}")
    if aoi_info["townships"]:
        print(f"  Townships: {', '.join(aoi_info['townships'])}")
    print(f"  Area: {aoi_km2:,.0f} km2   bbox "
          f"{bbox['west']:.3f},{bbox['south']:.3f} -> "
          f"{bbox['east']:.3f},{bbox['north']:.3f}")

    # Fail fast and clearly if either date's imagery is not in GEE yet.
    preflight(aoi)

    # ── GEE server-side processing ────────────────────────────────────────────
    print(f"\n> Loading S1 pre-flood  (+/-{DATE_WINDOW}d of {PRE_FLOOD_DATE}, {ORBIT_PASS}) ...")
    pre_raw  = load_s1(aoi, PRE_FLOOD_DATE,  DATE_WINDOW, ORBIT_PASS)

    print(f"> Loading S1 post-flood (+/-{DATE_WINDOW}d of {POST_FLOOD_DATE}, {ORBIT_PASS}) ...")
    post_raw = load_s1(aoi, POST_FLOOD_DATE, DATE_WINDOW, ORBIT_PASS)

    print("> Applying speckle filter ...")
    pre  = speckle_filter(pre_raw)
    post = speckle_filter(post_raw)

    print("> Detecting flooded pixels ...")
    flood_mask, perm_water, _ = detect_floods(pre, post, aoi)

    print(f"> Computing flood area at {EXPORT_SCALE} m ...")
    flood_km2 = compute_flood_area_km2(flood_mask, aoi)
    print(f"  => Estimated newly inundated area : {flood_km2:,.2f} km2 "
          f"({flood_km2 / aoi_km2 * 100:.1f}% of AOI)")

    valid_frac = compute_valid_fraction(flood_mask, aoi)
    print(f"  => AOI actually analysed          : {valid_frac:.1%}")
    if valid_frac < 0.995:
        print(f"  [!] WARNING: only {valid_frac:.0%} of the AOI was assessed. The rest "
              f"is sea, steep terrain, or lacked Sentinel-1 data on one of the two "
              f"dates. The flood figure covers only that {valid_frac:.0%}; "
              f"unassessed areas are unknown, not dry.")

    # ── Download rasters (all uint8 so 10 m stays tractable) ──────────────────
    # SAR dB is rescaled server-side to 1..255 (0 reserved as nodata): the PNG
    # is 8-bit anyway, so quantising here costs no visible fidelity but cuts
    # transfer and memory 4x versus float32.  This is a radiometric mapping for
    # display, not a spatial one -- pixel spacing stays at 10 m.
    def sar_byte(img: ee.Image) -> ee.Image:
        return (img.select("VV").unitScale(-25, 0)
                .multiply(254).add(1).clamp(1, 255).toByte())

    print("\n> Downloading rasters from GEE ...")
    flood_path = download_as_geotiff(flood_mask.toByte(), "flood.tif",
                                     aoi, EXPORT_SCALE, bbox)
    water_path = download_as_geotiff(perm_water.toByte(), "water.tif",
                                     aoi, EXPORT_SCALE, bbox)
    pre_path   = download_as_geotiff(sar_byte(pre),  "pre_vv.tif",
                                     aoi, SAR_DISPLAY_SCALE, bbox)
    post_path  = download_as_geotiff(sar_byte(post), "post_vv.tif",
                                     aoi, SAR_DISPLAY_SCALE, bbox)

    # ── Convert to palette PNG chunks (embedded or page-local files) ──────────
    print("> Encoding rasters as palette PNG chunks ...")
    clear_assets()
    overlays = {
        "flood":    geotiff_to_layer(flood_path, _palette_flood(), "flood"),
        "water":    geotiff_to_layer(water_path, _palette_water(), "water"),
        "pre_sar":  geotiff_to_layer(pre_path,   _palette_gray(),  "pre_sar"),
        "post_sar": geotiff_to_layer(post_path,  _palette_gray(),  "post_sar"),
    }

    # ── Build and save the map ────────────────────────────────────────────────
    print("> Building standalone HTML map ...")
    stats = {"flood_km2": flood_km2, "aoi_km2": aoi_km2,
             "valid_frac": valid_frac}
    build_map(overlays, stats, aoi_info, OUTPUT_HTML)

    # ── Clean up temp GeoTIFFs ────────────────────────────────────────────────
    shutil.rmtree(TEMP_DIR, ignore_errors=True)
    print("[OK] Temporary GeoTIFFs removed.")

    html_mb  = os.path.getsize(OUTPUT_HTML) / 1e6
    asset_mb = sum(v["mb"] for v in overlays.values() if v["external"])
    print(f"\n{'-' * 62}")
    print(f"  Done!  Open  '{OUTPUT_HTML}'  ({html_mb:.1f} MB) in any browser.")
    if asset_mb:
        n_files = len(glob.glob(os.path.join(ASSET_DIR, "*.png")))
        print(f"  Plus {ASSET_DIR}/ -- {n_files} PNG(s), {asset_mb:.1f} MB.")
        print(f"  Publish BOTH: the page loads its big layers from {ASSET_DIR}/.")
    print(f"  All layers at {EXPORT_SCALE} m native. No GEE URLs -- nothing expires.")
    print(f"{'-' * 62}")


if __name__ == "__main__":
    main()
