"""Validate core analysis on synthetic data (no OSM/census I/O needed)."""
import numpy as np
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point, box

import access_pipeline as fap

CRS = fap.METRIC_CRS

# --- Build a 21x21 grid graph, 100 m spacing => 2000x2000 m study area ---
G = nx.MultiDiGraph()
N, step = 21, 100.0
def nid(i, j): return i * N + j
for i in range(N):
    for j in range(N):
        G.add_node(nid(i, j), x=j * step, y=i * step)
for i in range(N):
    for j in range(N):
        if j + 1 < N:
            for a, b in [(nid(i, j), nid(i, j + 1)), (nid(i, j + 1), nid(i, j))]:
                G.add_edge(a, b, length=step)
        if i + 1 < N:
            for a, b in [(nid(i, j), nid(i + 1, j)), (nid(i + 1, j), nid(i, j))]:
                G.add_edge(a, b, length=step)
print(f"Grid graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Boundary = the full 2000x2000 m square
boundary = box(0, 0, 2000, 2000)

# Two food programs near corners
programs = gpd.GeoDataFrame(
    {"name": ["A", "B"]},
    geometry=[Point(200, 200), Point(1800, 1800)],
    crs=CRS,
)

# --- Network service area, 500 m ---
served_net, n_reach = fap.compute_network_service_area(G, programs, threshold_m=500)
print(f"Reachable nodes within 500 m: {n_reach}")
print(f"Network served area: {served_net.area/1e6:.3f} km^2 valid={served_net.is_valid}")

# --- Euclidean buffer service area, 500 m ---
served_buf = fap.compute_buffer_service_area(programs, radius_m=500)
print(f"Buffer served area:  {served_buf.area/1e6:.3f} km^2")

# --- Gaps ---
gap_net = fap.compute_gap(boundary, served_net).intersection(boundary)
gap_buf = fap.compute_gap(boundary, served_buf).intersection(boundary)
print(f"Gap (network): {gap_net.area/1e6:.3f} km^2 ({100*gap_net.area/boundary.area:.1f}%)")
print(f"Gap (buffer):  {gap_buf.area/1e6:.3f} km^2 ({100*gap_buf.area/boundary.area:.1f}%)")
assert gap_net.area >= gap_buf.area - 1e-6, "network gap should be >= buffer gap"

# --- Synthetic census DAs: 4x4 grid of 500m cells, pop 1000 each ---
das, pops = [], []
for i in range(4):
    for j in range(4):
        das.append(box(j*500, i*500, (j+1)*500, (i+1)*500))
        pops.append(1000)
da_gdf = gpd.GeoDataFrame({"population": pops}, geometry=das, crs=CRS)

pop_c = fap.population_in_gap(gap_net, da_gdf, "population", method="centroid")
pop_a = fap.population_in_gap(gap_net, da_gdf, "population", method="areal")
print(f"Population (centroid): {pop_c['population_in_gap']:.0f}/{pop_c['total_population']:.0f}"
      f" ({pop_c['pct_in_gap']:.1f}%)")
print(f"Population (areal):    {pop_a['population_in_gap']:.0f}/{pop_a['total_population']:.0f}"
      f" ({pop_a['pct_in_gap']:.1f}%)")
assert 0 <= pop_a["population_in_gap"] <= 16000
assert abs(pop_c["total_population"] - 16000) < 1e-6

# --- Dasymetric: residential land only in left half (x < 1000) ---
res_mask = gpd.GeoDataFrame(geometry=[box(0, 0, 1000, 2000)], crs=CRS)
pop_d = fap.population_in_gap(
    gap_net, da_gdf, "population", method="dasymetric", residential_mask=res_mask
)
print(f"Population (dasymetric): {pop_d['population_in_gap']:.0f}/{pop_d['total_population']:.0f}"
      f" ({pop_d['pct_in_gap']:.1f}%)")
assert 0 <= pop_d["population_in_gap"] <= 16000

# Sanity: if residential mask == whole area, dasymetric must equal areal
res_full = gpd.GeoDataFrame(geometry=[box(0, 0, 2000, 2000)], crs=CRS)
pop_full = fap.population_in_gap(
    gap_net, da_gdf, "population", method="dasymetric", residential_mask=res_full
)
print(f"Dasymetric(full mask): {pop_full['population_in_gap']:.1f}  "
      f"vs areal: {pop_a['population_in_gap']:.1f}")
assert abs(pop_full["population_in_gap"] - pop_a["population_in_gap"]) < 1.0, \
    "full-coverage dasymetric should equal areal"

print("\nALL CORE FUNCTIONS PASSED ✓")
