#!/usr/bin/env python3
"""
Food Program Accessibility Analysis — Mississauga
=================================================

Identifies areas (and, optionally, the resident population) located more than a
threshold walking distance (default 500 m) from the nearest food program, using
two methods:

    1. Euclidean buffer  — straight-line 500 m circles (optimistic baseline).
    2. Network service area — 500 m along the pedestrian street network
       (methodologically preferred; respects barriers like highways/rivers).

The "gap" for each method is the municipal boundary MINUS the union of served
areas. Reporting both quantifies the optimism bias of Euclidean accessibility.

Core analysis functions (compute_*, population_in_gap) depend only on
networkx / shapely / geopandas / numpy and are unit-testable without any
network I/O. OSMnx is imported lazily inside the data loaders.

Author: pipeline scaffold for X. Yan, GISpark Lab
CRS for all metric computation: EPSG:26917 (NAD83 / UTM zone 17N).
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import geopandas as gpd
import networkx as nx
from shapely.geometry import LineString, Point
from shapely.ops import unary_union

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

METRIC_CRS = "EPSG:26917"          # NAD83 / UTM 17N — metres, for Mississauga
GEO_CRS = "EPSG:4326"              # WGS84 — for GeoJSON output
DEFAULT_PLACE = "Mississauga, Ontario, Canada"
DEFAULT_THRESHOLD_M = 500.0        # accessibility cutoff (network & buffer)
DEFAULT_NETWORK_BUFFER_M = 800.0   # extend network beyond boundary (edge effect)
DEFAULT_EDGE_BUFFER_M = 25.0       # half-width of street "reach" band

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("food_access")


# --------------------------------------------------------------------------- #
# Internal graph helpers (no OSMnx dependency — keeps core testable)
# --------------------------------------------------------------------------- #

def _node_coords(G: nx.Graph) -> tuple[list, np.ndarray]:
    """Return (node_ids, Nx2 array of [x, y]) from node 'x'/'y' attributes."""
    nodes = list(G.nodes)
    xy = np.array([[G.nodes[n]["x"], G.nodes[n]["y"]] for n in nodes], dtype=float)
    return nodes, xy


def _nearest_nodes(G: nx.Graph, xs: Iterable[float], ys: Iterable[float]) -> list:
    """Nearest graph node to each (x, y). Brute-force; fine for a few hundred pts."""
    nodes, xy = _node_coords(G)
    out = []
    for x, y in zip(np.atleast_1d(xs), np.atleast_1d(ys)):
        d2 = (xy[:, 0] - x) ** 2 + (xy[:, 1] - y) ** 2
        out.append(nodes[int(np.argmin(d2))])
    return out


def _served_polygon_from_nodes(G: nx.Graph, reachable_nodes, edge_buffer: float):
    """
    Build a service-area polygon from the set of reachable nodes by buffering the
    street segments that connect them. Falls back to buffering node points if no
    edges qualify (e.g. a single isolated reachable node).
    """
    rset = set(reachable_nodes)
    geoms = []
    for u, v, data in G.edges(data=True):
        if u in rset and v in rset:
            geom = data.get("geometry")
            if geom is None:
                geom = LineString(
                    [(G.nodes[u]["x"], G.nodes[u]["y"]),
                     (G.nodes[v]["x"], G.nodes[v]["y"])]
                )
            geoms.append(geom)
    if not geoms:
        pts = [Point(G.nodes[n]["x"], G.nodes[n]["y"]) for n in reachable_nodes]
        if not pts:
            return None
        return unary_union([p.buffer(edge_buffer) for p in pts])
    return unary_union(geoms).buffer(edge_buffer)


# --------------------------------------------------------------------------- #
# Core analysis (testable, I/O-free)
# --------------------------------------------------------------------------- #

def compute_network_service_area(
    G: nx.Graph,
    programs_gdf: gpd.GeoDataFrame,
    threshold_m: float = DEFAULT_THRESHOLD_M,
    edge_buffer_m: float = DEFAULT_EDGE_BUFFER_M,
    weight: str = "length",
):
    """
    Union of all nodes reachable within `threshold_m` network metres of ANY food
    program (multi-source Dijkstra), rendered as a street-following polygon.

    Returns (served_polygon, n_reachable_nodes). Both G and programs_gdf must be
    in the same projected (metric) CRS.
    """
    if programs_gdf.empty:
        return None, 0
    xs = programs_gdf.geometry.x.to_numpy()
    ys = programs_gdf.geometry.y.to_numpy()
    sources = set(_nearest_nodes(G, xs, ys))
    lengths = nx.multi_source_dijkstra_path_length(
        G, sources, cutoff=threshold_m, weight=weight
    )
    reachable = list(lengths.keys())
    served = _served_polygon_from_nodes(G, reachable, edge_buffer_m)
    return served, len(reachable)


def compute_buffer_service_area(
    programs_gdf: gpd.GeoDataFrame, radius_m: float = DEFAULT_THRESHOLD_M
):
    """Union of straight-line `radius_m` buffers around food programs (metric CRS)."""
    if programs_gdf.empty:
        return None
    return unary_union(programs_gdf.geometry.buffer(radius_m).to_numpy())


def compute_gap(boundary_geom, served_geom):
    """Underserved region = boundary minus served area."""
    if served_geom is None:
        return boundary_geom
    return boundary_geom.difference(served_geom)


def population_in_gap(
    gap_geom,
    da_gdf: gpd.GeoDataFrame,
    pop_field: str,
    method: str = "centroid",
    income_field: Optional[str] = None,
    residential_mask: Optional[gpd.GeoDataFrame] = None,
) -> dict:
    """
    Estimate population inside the gap using Census dissemination areas (DAs).

    method='centroid'   : a DA counts as in-gap if its centroid falls in the gap
                          (fast, standard, slight boundary error).
    method='areal'      : population apportioned by the share of each DA's AREA
                          intersecting the gap (assumes uniform density per DA).
    method='dasymetric' : binary dasymetric — population is assumed uniform only
                          over RESIDENTIAL land within each DA, then apportioned
                          by the share of that residential land falling in the
                          gap. Requires `residential_mask` (e.g. OSM
                          landuse=residential). Most accurate where DAs mix
                          residential with industrial/airport/green land —
                          exactly the Mississauga case.

    Returns a dict of summary statistics.
    """
    da = da_gdf.to_crs(METRIC_CRS).copy().reset_index(drop=True)
    da["__da_id"] = np.arange(len(da))
    total_pop = float(da[pop_field].sum())
    result = {"total_population": total_pop, "method": method}

    def _frac_to_pop(frac_by_da):
        """frac_by_da: Series indexed by __da_id -> apportioned population dict."""
        frac = frac_by_da.reindex(da["__da_id"]).fillna(0.0).to_numpy()
        out = {"population_in_gap": float((da[pop_field].to_numpy() * frac).sum())}
        if income_field:
            out["low_income_in_gap"] = float((da[income_field].to_numpy() * frac).sum())
        return out

    if method == "centroid":
        in_gap = da[da.geometry.centroid.within(gap_geom)]
        result["population_in_gap"] = float(in_gap[pop_field].sum())
        if income_field:
            result["low_income_in_gap"] = float(in_gap[income_field].sum())

    elif method == "areal":
        da["__da_area"] = da.geometry.area
        gap_gdf = gpd.GeoDataFrame(geometry=[gap_geom], crs=METRIC_CRS)
        inter = gpd.overlay(da, gap_gdf, how="intersection")
        if inter.empty:
            result["population_in_gap"] = 0.0
        else:
            inter["__frac"] = inter.geometry.area / inter["__da_area"]
            frac_by_da = inter.groupby("__da_id")["__frac"].sum()
            result.update(_frac_to_pop(frac_by_da))

    elif method == "dasymetric":
        if residential_mask is None or residential_mask.empty:
            raise ValueError(
                "method='dasymetric' requires a non-empty residential_mask."
            )
        mask_geom = unary_union(residential_mask.to_crs(METRIC_CRS).geometry.to_numpy())
        mask_gdf = gpd.GeoDataFrame(geometry=[mask_geom], crs=METRIC_CRS)
        gap_gdf = gpd.GeoDataFrame(geometry=[gap_geom], crs=METRIC_CRS)

        # residential land within each DA (the dasymetric "population surface")
        res = gpd.overlay(da[["__da_id", "geometry"]], mask_gdf, how="intersection")
        if res.empty:
            # no residential land intersects any DA -> fall back to areal
            log.warning("Residential mask covers no DA area; falling back to areal.")
            return population_in_gap(gap_geom, da_gdf, pop_field, "areal", income_field)
        res["__res_area"] = res.geometry.area
        res_area_by_da = res.groupby("__da_id")["__res_area"].sum()

        # residential land within each DA AND inside the gap
        res_gap = gpd.overlay(res, gap_gdf, how="intersection")
        if res_gap.empty:
            result["population_in_gap"] = 0.0
            if income_field:
                result["low_income_in_gap"] = 0.0
        else:
            res_gap["__rg_area"] = res_gap.geometry.area
            rg_by_da = res_gap.groupby("__da_id")["__rg_area"].sum()
            frac_by_da = (rg_by_da / res_area_by_da)  # NaN where no res in gap
            result.update(_frac_to_pop(frac_by_da))
        # DAs with no residential land contribute 0 (frac=0) — they hold no pop.

    else:
        raise ValueError(
            f"Unknown method: {method!r} (use 'centroid', 'areal', or 'dasymetric')"
        )

    if total_pop > 0:
        result["pct_in_gap"] = 100.0 * result["population_in_gap"] / total_pop
    return result


# --------------------------------------------------------------------------- #
# Data loaders (OSMnx imported lazily)
# --------------------------------------------------------------------------- #

def load_boundary(place: str) -> gpd.GeoDataFrame:
    """Geocode the municipal boundary to a GeoDataFrame (WGS84)."""
    import osmnx as ox
    log.info("Geocoding boundary: %s", place)
    gdf = ox.geocode_to_gdf(place)
    return gdf.to_crs(GEO_CRS)


def load_network(boundary_gdf: gpd.GeoDataFrame, network_buffer_m: float) -> nx.MultiDiGraph:
    """
    Download the walkable street network for the boundary, buffered outward by
    `network_buffer_m` to avoid edge effects, and project to the metric CRS.
    """
    import osmnx as ox
    boundary_proj = boundary_gdf.to_crs(METRIC_CRS)
    net_poly_proj = unary_union(boundary_proj.geometry.to_numpy()).buffer(network_buffer_m)
    net_poly_geo = gpd.GeoSeries([net_poly_proj], crs=METRIC_CRS).to_crs(GEO_CRS).iloc[0]
    log.info("Downloading walk network (this can take a few minutes)…")
    G = ox.graph_from_polygon(net_poly_geo, network_type="walk")
    G = ox.project_graph(G, to_crs=METRIC_CRS)
    log.info("Network: %d nodes, %d edges", G.number_of_nodes(), G.number_of_edges())
    return G


def load_food_programs(
    path: Optional[str],
    boundary_gdf: gpd.GeoDataFrame,
    seed_from_osm: bool,
    network_buffer_m: float,
    geocode_missing: bool = False,
) -> gpd.GeoDataFrame:
    """
    Load food-program points. Priority:
      1. CSV/GeoJSON at `path`.
      2. OSM seed (`seed_from_osm`) — CANDIDATES ONLY, must be validated.
    CSV is expected to have lat/lon columns, or an 'address' column (+ geocode).
    """
    if path:
        return _load_programs_from_file(path, geocode_missing)
    if seed_from_osm:
        return _seed_programs_from_osm(boundary_gdf, network_buffer_m)
    raise ValueError(
        "No food-program data. Provide --food-programs PATH or pass --seed-from-osm."
    )


def _load_programs_from_file(path: str, geocode_missing: bool) -> gpd.GeoDataFrame:
    p = Path(path)
    if p.suffix.lower() in {".geojson", ".json", ".gpkg", ".shp"}:
        gdf = gpd.read_file(path)
        log.info("Loaded %d programs from %s", len(gdf), path)
        return gdf.to_crs(METRIC_CRS)

    import pandas as pd
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    lat_c = next((cols[c] for c in ("lat", "latitude", "y") if c in cols), None)
    lon_c = next((cols[c] for c in ("lon", "lng", "longitude", "x") if c in cols), None)

    if lat_c and lon_c:
        gdf = gpd.GeoDataFrame(
            df, geometry=gpd.points_from_xy(df[lon_c], df[lat_c]), crs=GEO_CRS
        )
    elif "address" in cols and geocode_missing:
        import osmnx as ox
        pts = []
        for addr in df[cols["address"]]:
            try:
                lat, lon = ox.geocode(str(addr))
                pts.append(Point(lon, lat))
            except Exception as e:                       # noqa: BLE001
                log.warning("Geocode failed for %r: %s", addr, e)
                pts.append(None)
            time.sleep(1.0)                              # Nominatim courtesy
        gdf = gpd.GeoDataFrame(df, geometry=pts, crs=GEO_CRS).dropna(subset=["geometry"])
    else:
        raise ValueError(
            "CSV needs lat/lon columns, or an 'address' column with --geocode-missing."
        )
    log.info("Loaded %d programs from %s", len(gdf), path)
    return gdf.to_crs(METRIC_CRS)


def _seed_programs_from_osm(
    boundary_gdf: gpd.GeoDataFrame, network_buffer_m: float
) -> gpd.GeoDataFrame:
    import osmnx as ox
    boundary_proj = boundary_gdf.to_crs(METRIC_CRS)
    poly_proj = unary_union(boundary_proj.geometry.to_numpy()).buffer(network_buffer_m)
    poly_geo = gpd.GeoSeries([poly_proj], crs=METRIC_CRS).to_crs(GEO_CRS).iloc[0]
    tags = {
        "social_facility": ["food_bank", "soup_kitchen", "food_pantry"],
        "amenity": ["social_facility", "social_centre"],
    }
    log.warning(
        "Seeding food programs from OSM — these are CANDIDATES ONLY. "
        "Validate against The Mississauga Food Bank member directory and "
        "Region of Peel / City of Mississauga open data before reporting."
    )
    feats = ox.features_from_polygon(poly_geo, tags)
    feats = feats[~feats.geometry.is_empty & feats.geometry.notna()].to_crs(METRIC_CRS)
    feats["geometry"] = feats.geometry.centroid     # points for analysis
    keep = [c for c in ("name", "social_facility", "amenity") if c in feats.columns]
    feats = feats[keep + ["geometry"]].reset_index(drop=True)
    log.info("OSM seed produced %d candidate food-program points", len(feats))
    return feats


def load_residential_mask(
    path: Optional[str],
    boundary_gdf: gpd.GeoDataFrame,
    network_buffer_m: float,
    use_buildings: bool = False,
) -> gpd.GeoDataFrame:
    """
    Load the residential land mask for dasymetric apportionment.
      1. polygon file at `path` (GeoJSON/SHP/GPKG), or
      2. OSM landuse=residential (+ optionally building footprints).
    Returns polygons in the metric CRS.
    """
    if path:
        gdf = gpd.read_file(path).to_crs(METRIC_CRS)
        gdf = gdf[gdf.geometry.geom_type.isin(["Polygon", "MultiPolygon"])]
        log.info("Loaded residential mask: %d polygons from %s", len(gdf), path)
        return gdf

    import osmnx as ox
    boundary_proj = boundary_gdf.to_crs(METRIC_CRS)
    poly_proj = unary_union(boundary_proj.geometry.to_numpy()).buffer(network_buffer_m)
    poly_geo = gpd.GeoSeries([poly_proj], crs=METRIC_CRS).to_crs(GEO_CRS).iloc[0]
    tags = {"landuse": ["residential"]}
    if use_buildings:
        tags["building"] = True   # finer surface, but heavier download
    log.info("Downloading residential land mask from OSM…")
    feats = ox.features_from_polygon(poly_geo, tags)
    feats = feats[feats.geometry.notna() & ~feats.geometry.is_empty]
    feats = feats[feats.geometry.geom_type.isin(["Polygon", "MultiPolygon"])]
    log.info("Residential mask: %d polygons", len(feats))
    return feats.to_crs(METRIC_CRS)


# --------------------------------------------------------------------------- #
# Output helpers
# --------------------------------------------------------------------------- #

def _to_geojson(geom, crs, path: Path):
    gpd.GeoSeries([geom], crs=crs).to_crs(GEO_CRS).to_file(path, driver="GeoJSON")
    log.info("Wrote %s", path)


def save_outputs(outdir: Path, **layers):
    outdir.mkdir(parents=True, exist_ok=True)
    for name, geom in layers.items():
        if geom is not None:
            _to_geojson(geom, METRIC_CRS, outdir / f"{name}.geojson")


def quick_map(outdir: Path, boundary, programs_gdf, gap_net, gap_buf):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 10))
    gpd.GeoSeries([boundary], crs=METRIC_CRS).plot(
        ax=ax, facecolor="none", edgecolor="black", linewidth=1.2
    )
    if gap_buf is not None:
        gpd.GeoSeries([gap_buf], crs=METRIC_CRS).plot(
            ax=ax, color="#fdae6b", alpha=0.45, label="Gap (Euclidean buffer)"
        )
    if gap_net is not None:
        gpd.GeoSeries([gap_net], crs=METRIC_CRS).plot(
            ax=ax, color="#d73027", alpha=0.55, label="Gap (network)"
        )
    programs_gdf.plot(ax=ax, color="#1a9850", markersize=14, label="Food programs")
    ax.set_title("Areas > 500 m from a food program — Mississauga", fontsize=13)
    ax.legend(loc="lower right")
    ax.set_axis_off()
    fig.tight_layout()
    out = outdir / "access_gap_map.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    log.info("Wrote %s", out)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #

def run(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)

    boundary_gdf = load_boundary(args.place)
    boundary_proj = unary_union(boundary_gdf.to_crs(METRIC_CRS).geometry.to_numpy())

    programs = load_food_programs(
        args.food_programs, boundary_gdf, args.seed_from_osm,
        args.network_buffer, args.geocode_missing,
    )
    if programs.empty:
        log.error("No food programs loaded — aborting.")
        sys.exit(1)

    # --- Euclidean buffer method ---
    served_buf = compute_buffer_service_area(programs, args.threshold)
    gap_buf = compute_gap(boundary_proj, served_buf).intersection(boundary_proj)

    # --- Network method ---
    G = load_network(boundary_gdf, args.network_buffer)
    served_net, n_reach = compute_network_service_area(
        G, programs, args.threshold, args.edge_buffer
    )
    # clip served area to boundary for fair area comparison
    served_net_clip = served_net.intersection(boundary_proj) if served_net else None
    gap_net = compute_gap(boundary_proj, served_net).intersection(boundary_proj)

    # --- Geometry stats ---
    total_area = boundary_proj.area
    stats = {
        "boundary_area_km2": total_area / 1e6,
        "n_food_programs": len(programs),
        "n_reachable_nodes": n_reach,
        "gap_buffer_km2": gap_buf.area / 1e6,
        "gap_network_km2": gap_net.area / 1e6,
        "gap_buffer_pct": 100 * gap_buf.area / total_area,
        "gap_network_pct": 100 * gap_net.area / total_area,
    }

    # --- Population overlay (optional) ---
    pop_stats = {}
    if args.census:
        da = gpd.read_file(args.census)
        res_mask = None
        if args.pop_method == "dasymetric":
            res_mask = load_residential_mask(
                args.residential_mask, boundary_gdf,
                args.network_buffer, args.use_buildings,
            )
        pop_net = population_in_gap(
            gap_net, da, args.pop_field, args.pop_method, args.income_field, res_mask
        )
        pop_buf = population_in_gap(
            gap_buf, da, args.pop_field, args.pop_method, args.income_field, res_mask
        )
        pop_stats = {"network": pop_net, "buffer": pop_buf}

    # --- Save ---
    save_outputs(
        outdir,
        food_programs=unary_union(programs.geometry.to_numpy()),
        served_area_buffer=served_buf,
        served_area_network=served_net_clip,
        gap_buffer=gap_buf,
        gap_network=gap_net,
    )
    programs.to_crs(GEO_CRS).to_file(outdir / "food_programs_points.geojson", driver="GeoJSON")

    if args.quick_map:
        quick_map(outdir, boundary_proj, programs.to_crs(METRIC_CRS), gap_net, gap_buf)

    # --- Report ---
    _print_report(stats, pop_stats, args, outdir)


def _print_report(stats, pop_stats, args, outdir):
    lines = [
        "",
        "=" * 64,
        "FOOD PROGRAM ACCESSIBILITY — SUMMARY",
        "=" * 64,
        f"Place                    : {args.place}",
        f"Threshold                : {args.threshold:.0f} m",
        f"Food programs            : {stats['n_food_programs']}",
        f"Municipal area           : {stats['boundary_area_km2']:.1f} km^2",
        "-" * 64,
        f"Gap (Euclidean buffer)   : {stats['gap_buffer_km2']:.2f} km^2  "
        f"({stats['gap_buffer_pct']:.1f}% of area)",
        f"Gap (network 500 m)      : {stats['gap_network_km2']:.2f} km^2  "
        f"({stats['gap_network_pct']:.1f}% of area)",
        f"Network/buffer ratio     : "
        f"{stats['gap_network_km2'] / max(stats['gap_buffer_km2'], 1e-9):.2f}x larger gap",
    ]
    if pop_stats:
        n, b = pop_stats["network"], pop_stats["buffer"]
        lines += [
            "-" * 64,
            f"Total population         : {n['total_population']:,.0f}",
            f"Pop in gap (buffer)      : {b['population_in_gap']:,.0f} "
            f"({b.get('pct_in_gap', 0):.1f}%)",
            f"Pop in gap (network)     : {n['population_in_gap']:,.0f} "
            f"({n.get('pct_in_gap', 0):.1f}%)",
        ]
        if "low_income_in_gap" in n:
            lines += [
                f"Low-income pop in gap    : {n['low_income_in_gap']:,.0f} (network)",
            ]
    lines += ["=" * 64, f"Outputs written to: {outdir.resolve()}", ""]
    report = "\n".join(lines)
    print(report)
    (Path(args.outdir) / "summary.txt").write_text(report)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Food program accessibility (buffer vs network) for Mississauga."
    )
    p.add_argument("--place", default=DEFAULT_PLACE, help="Geocodable place name.")
    p.add_argument("--food-programs", default=None,
                   help="CSV (lat/lon or address) or GeoJSON of food programs.")
    p.add_argument("--seed-from-osm", action="store_true",
                   help="Seed candidate programs from OSM (validate before use!).")
    p.add_argument("--geocode-missing", action="store_true",
                   help="Geocode CSV addresses lacking coordinates (slow).")
    p.add_argument("--census", default=None,
                   help="DA boundary file (GeoJSON/SHP/GPKG) with a population field.")
    p.add_argument("--pop-field", default="population",
                   help="Population column name in the census file.")
    p.add_argument("--income-field", default=None,
                   help="Optional low-income population column.")
    p.add_argument("--pop-method", default="centroid",
                   choices=["centroid", "areal", "dasymetric"],
                   help="Population apportionment method.")
    p.add_argument("--residential-mask", default=None,
                   help="Residential land polygon file (for --pop-method dasymetric). "
                        "If omitted, pulled from OSM landuse=residential.")
    p.add_argument("--use-buildings", action="store_true",
                   help="Add OSM building footprints to the residential mask (finer, heavier).")
    p.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD_M,
                   help="Accessibility distance threshold in metres.")
    p.add_argument("--network-buffer", type=float, default=DEFAULT_NETWORK_BUFFER_M,
                   help="Buffer beyond boundary when downloading network (edge effect).")
    p.add_argument("--edge-buffer", type=float, default=DEFAULT_EDGE_BUFFER_M,
                   help="Half-width of the street reach band for service areas.")
    p.add_argument("--outdir", default="output", help="Output directory.")
    p.add_argument("--quick-map", action="store_true", help="Render a PNG sanity map.")
    return p


def main(argv: Optional[list[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
