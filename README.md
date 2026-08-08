# Lemyethna Flood Inundation Map — 27 Jul → 8 Aug 2026

Sentinel-1 SAR flood inundation mapping over four merged townships in the
Ayeyarwady Region, Myanmar: **Lemyethna, Hinthada, Kyonpyaw and Yegyi**.

**Live map:** https://geonet-myanmar.github.io/lemyethna-flood-20260806/

## Status

`flood_mapping.py` is configured for the **2026-07-27 → 2026-08-08** pair on
ascending track 143, at native 10 m for every layer.

The 2026-08-08 scene was acquired on the day of writing and has **not yet been
ingested into the Earth Engine catalogue** (latest available over this AOI:
2026-08-05). The script's pre-flight check stops with a message listing what
*is* available. Sentinel-1 normally reaches GEE a few days after acquisition —
re-run `python flood_mapping.py` unchanged once it lands and it will produce
`index.html` plus `assets/` for the new pair.

Until then `index.html` remains the previous product (24 Jul → 5 Aug, descending
track 106, SAR backdrops at 30 m).

## Configuration

| | |
|---|---|
| Pre-flood image | 2026-07-27 |
| Post-flood image | 2026-08-08 |
| Orbit | Sentinel-1 **ascending, relative orbit (track) 143** |
| Scene coverage | 100 % of AOI on 2026-07-27 (post date pending ingest) |
| AOI area | 4,153 km² |
| Resolution | **10 m (native Sentinel-1 GRD) — every layer** |

### Why the track changed

The orbit is a property of the *dates*, not of the AOI. This AOI is imaged by
three tracks, each on a 12-day repeat:

| Track | Pass | Acquisitions |
|---|---|---|
| 143 | Ascending | 03 Jul, 15 Jul, **27 Jul**, **08 Aug** |
| 41 | Ascending | 08 Jul, 20 Jul, 01 Aug |
| 106 | Descending | 12 Jul, 24 Jul, 05 Aug |

27 Jul and 8 Aug are consecutive passes of **ascending track 143** — identical
viewing geometry, exactly 12 days apart, 100 % coverage. The previous run used
descending track 106 because that track carried its date pair (24 Jul / 5 Aug).

This matters: mixing ascending and descending geometry introduces a systematic
backscatter bias that destroys the flood signal (an earlier run over a different
AOI returned 0 km² for exactly this reason). Never pair dates across tracks.

## Method

Change detection on VV backscatter. Open floodwater reflects radar away from
the sensor, so newly flooded ground shows a sharp drop in VV.

1. Mean composite of Sentinel-1 IW GRD scenes within ±1 day of each date.
2. Both composites are taken from the **same orbit pass and same relative
   orbit** (ascending 143).
3. Boxcar speckle filter, 30 m radius.
4. A pixel is classified flooded where `post − pre < −3.0 dB`, **and** it is not
   JRC permanent/seasonal water, **and** terrain slope < 5°.
5. Flood area is reduced at 10 m, so nothing is generalised away.

### Data sources

- Copernicus **Sentinel-1** GRD IW (`COPERNICUS/S1_GRD`) via Google Earth Engine
- **JRC** Global Surface Water v1.4 — permanent/seasonal water mask
- **WWF HydroSHEDS** 3-arcsec void-filled DEM — slope mask (its sea NoData also
  serves as a land mask)

## Resolution

**Nothing is downsampled.** The flood mask, the permanent-water layer, the area
statistic *and* both SAR VV backdrops are produced and displayed at Sentinel-1's
native **10 m** GRD pixel spacing. `SAR_DISPLAY_SCALE` was 30 m in the previous
version; it is now 10 m.

At 10 m this AOI is roughly 9,500 × 7,900 px (75 Mpx) per layer, which is why
the output is no longer a single file:

| Layer @ 10 m | PNG payload |
|---|---|
| Flood / water mask | ~1 MB (binary, compresses ~100×) |
| SAR VV backdrop | ~56 MB each (speckle, near-incompressible) |

Base64 embedding inflates by 1.33×, so inlining two native SAR backdrops would
produce a ~150 MB `index.html` — over **GitHub's hard 100 MB per-file limit**,
and slow for a browser even to parse. So the generator now:

- **Chunks** each layer into pieces of at most `MAX_CHUNK_PX` (16 Mpx) and drops
  chunks that are entirely nodata — the bbox corners outside the AOI polygon. A
  single 75 Mpx image decodes to ~300 MB of RGBA in the browser and fails
  outright on mobile Safari; ~13 Mpx chunks decode independently and can be
  freed independently.
- **Embeds** a layer as base64 if its whole payload is under `EMBED_MAX_MB`
  (8 MB), otherwise writes it to `assets/` and references it by relative path.

Chunks of a layer share one Leaflet `FeatureGroup`, so they toggle as a unit and
appear as a single entry in the layer control. The two SAR layers start hidden,
and Leaflet does not fetch a group's images until it is switched on — so the
heavy bytes are only transferred if a user asks for them.

**Deploying therefore means publishing `index.html` *and* the `assets/` folder.**
`assets/` is deliberately not in `.gitignore`; ignoring it publishes a page whose
heavy layers 404.

## Why the page does not expire

`geemap.addLayer()` writes GEE tile-server URLs into the HTML, and those URLs
expire after about two days, leaving a blank map. Instead every raster is
downloaded, converted to a palette PNG and served by the page itself — as a
base64 data URI or as a file in `assets/`. Nothing points at a Google server.

Palette ("P" mode) PNGs are used rather than RGBA: at 1 byte/px instead of 4, a
75 Mpx layer occupies 75 MB in memory rather than 300 MB, and compresses far
better.

## Reproducing

```bash
pip install earthengine-api geemap folium rasterio numpy Pillow
earthengine authenticate          # first time only
python flood_mapping.py
```

`flood_mapping.py` writes **`index.html`** at the repo root plus **`assets/*.png`**
— together they *are* the published dashboard. Do not hand-edit them; re-run the
script.

Set `GEE_PROJECT` in the script to your own Earth Engine Cloud project. Note that
a bare `ee.Initialize()` without a project will not work here.

### Changing the AOI or dates

The AOI is read from `merged_townships.geojson` (a FeatureCollection, Feature or
bare geometry all work). Dates and orbit are set at the top of the script:

```python
AOI_GEOJSON     = "merged_townships.geojson"
PRE_FLOOD_DATE  = "2026-07-27"
POST_FLOOD_DATE = "2026-08-08"
ORBIT_PASS      = "ASCENDING"
REL_ORBIT       = 143
```

`ORBIT_PASS` / `REL_ORBIT` are **specific to this AOI and date pair** and must be
re-checked when either changes — a track that covers one AOI may miss another
entirely, and a track that carries one date pair will not carry the next. The
pre-flight check reports scene count and AOI coverage for both dates and, if
imagery is missing, prints every acquisition currently in the catalogue with its
pass and track so you can pick a valid same-track pair.

To trade resolution for a smaller page, set `SAR_DISPLAY_SCALE` back to 30 —
that affects the display backdrops only, never the flood product.

## Caveats

- **Blank ≠ flood-free.** Areas that are steep terrain, permanent water, or
  lacked data on one of the two dates are not assessed. The map reports the
  assessed fraction; unassessed areas are unknown, not dry.
- SAR change detection cannot always separate open floodwater from flooded rice
  paddy or wet soil — relevant across this delta landscape. Check the flood layer
  against the pre/post VV backdrops before treating figures as definitive.
- The −3.0 dB threshold is a widely used default, not a locally calibrated value.
- Turning on a 10 m SAR backdrop transfers tens of megabytes and is slow on a
  poor connection. That is the intended trade for native resolution.

## Files

| Path | Purpose |
|---|---|
| `index.html` | The published dashboard (generated — do not edit) |
| `assets/*.png` | Native-resolution layers too large to embed (generated) |
| `flood_mapping.py` | Full pipeline: GEE → GeoTIFF → palette PNG → folium map |
| `merged_townships.geojson` | AOI boundary (4 merged townships) |

## Licence / attribution

Contains modified Copernicus Sentinel data (2026), processed via Google Earth
Engine. JRC Global Surface Water and WWF HydroSHEDS used under their respective
terms.
