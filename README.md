# Food_Program_Accessibility_Mississauga
Who lives more than a 500 m walk from a food program? A reproducible food-access GIS pipeline for Mississauga.
# Food Program Accessibility — Mississauga

Identifies areas (and resident population) more than **500 m** from the nearest
food program, comparing two access models:

| Method | What it measures | Bias |
|---|---|---|
| **Euclidean buffer** | straight-line 500 m | optimistic — ignores street network & barriers |
| **Network service area** | 500 m along the walkable network | realistic for pedestrian access |

The "gap" = municipal boundary − union of served areas. Reporting both
quantifies the optimism bias of Euclidean accessibility — the **network/buffer
gap ratio is itself a finding**, directly in the Kwan critique of
container-based accessibility.

---

## Setup

GeoPandas/OSMnx depend on GDAL. Conda is the path of least resistance:

```bash
conda create -n foodaccess python=3.11
conda activate foodaccess
conda install -c conda-forge osmnx geopandas matplotlib
# or, with pip (wheels ship GDAL via pyogrio):
pip install -r requirements.txt
```

Targets **OSMnx 2.x** (verified against 2.1.0).

---

## Data you need to supply

### 1. Food program locations (required)

The OSM seed (`--seed-from-osm`) is a **convenience starting point only** — OSM
coverage of food banks/meal programs is incomplete and unverified. The
authoritative sources to build/validate your point file:

- **The Mississauga Food Bank** — member-agency / program directory
- **Region of Peel Open Data** — community & human-services locations
- **City of Mississauga Open Data Portal** — community centres, facilities

Provide as **CSV** (`name,lat,lon` — or a `name,address` column with
`--geocode-missing`) or **GeoJSON**.

### 2. Census dissemination areas (optional, for population)

For the population overlay, download **StatCan 2021 Dissemination Area
boundaries** clipped to Peel/Mississauga and join the DA population count (and,
optionally, a low-income count e.g. LIM-AT). Pass with `--census`, naming the
columns via `--pop-field` / `--income-field`.

> **Why population weighting matters here:** raw "area > 500 m" is misleading in
> Mississauga — Pearson Airport plus large industrial lands are huge and
> unpopulated, and would dominate any area-based statistic while housing nobody.
> Always report population, not just area.

---

## Usage

```bash
# Quick look using an OSM seed (validate the points before citing anything):
python food_access_pipeline.py --seed-from-osm --quick-map

# Real run with your own point file:
python food_access_pipeline.py --food-programs programs.csv --quick-map

# Full run with population overlay:
python food_access_pipeline.py \
    --food-programs programs.geojson \
    --census peel_da_2021.geojson --pop-field pop2021 \
    --income-field lim_at --pop-method areal \
    --quick-map

# Dasymetric refinement (population redistributed onto residential land):
python food_access_pipeline.py \
    --food-programs programs.geojson \
    --census peel_da_2021.geojson --pop-field pop2021 \
    --pop-method dasymetric \
    --residential-mask peel_residential.geojson \
    --quick-map
# (omit --residential-mask to pull OSM landuse=residential automatically;
#  add --use-buildings for a finer but heavier building-footprint surface)
```

### Outputs (`output/`)

- `food_programs_points.geojson`, `served_area_{network,buffer}.geojson`
- `gap_network.geojson`, `gap_buffer.geojson` — the underserved regions
- `summary.txt` — headline stats (area + population, both methods)
- `access_gap_map.png` — sanity-check map (with `--quick-map`)

---

## Methodology notes

- **CRS**: all distance work in EPSG:26917 (NAD83 / UTM 17N). Outputs reprojected
  to WGS84 for GeoJSON.
- **Edge effects**: the network and OSM seed are buffered `--network-buffer`
  (default 800 m) beyond the boundary, and programs *outside* Mississauga that
  serve border neighbourhoods should be included in your point file — otherwise
  boundary areas are falsely flagged as gaps.
- **Network service area**: multi-source Dijkstra (`cutoff=500 m`) over the walk
  network; reachable street segments buffered by `--edge-buffer` (25 m) and
  dissolved. This *under*-claims slightly at the frontier — the safe error
  direction for a food-desert study.
- **Population method**: `centroid` is fast/standard; `areal` apportions by
  intersection area (better for large DAs); `dasymetric` apportions only over
  **residential** land within each DA (best where DAs mix housing with
  airport/industrial/green land — the Mississauga case). They diverge most when
  served areas are thin relative to DA size — prefer `dasymetric` for headline
  numbers, fall back to `areal`, and report the sensitivity. When the
  residential mask covers a whole DA, `dasymetric` reduces exactly to `areal`.
- **Threshold**: 500 m ≈ 6–8 min walk. Stringent but defensible for small local
  food programs (vs. the 1 km / 1 mi common in supermarket studies). State it
  explicitly in your methods.

## Extending toward an iOS/MapKit front end

The GeoJSON outputs (WGS84) drop straight into MapKit `MKGeoJSONDecoder` →
`MKPolygon` overlays. Run the analysis here, ship `gap_network.geojson` +
`food_programs_points.geojson` into the Swift app as bundled or fetched assets.
