# build_gridded_abundance_by_cruise_and_taxon.py

"""
Builds a GeoLibre-renderable gridded NetCDF of ichthyoplankton abundance
from swfsc_ichthyo.nc, with cruise and taxon as independent GeoLibre
"leading axis" sliders over a shared regular lat/lon grid:
    abundance(cruise, taxon, lat, lon)

Supersedes join_occurrence_abundance_to_flat_point_netcdf.py's flat CF
"point" output -- confirmed via GeoLibre's own source
(packages/plugins/src/plugins/local-netcdf.ts) that its local NetCDF/HDF
reader only ever considers "numeric, 2-D or higher" variables renderable;
a flat 1-D point table, however CF-correctly annotated, never qualifies.
GeoLibre's actual model is a rectilinear raster with optional leading
("slider") axes -- `selector: Record<string, number>` takes more than one
axis at once, and its docs say "The Time Slider uses the same selector
path" for any non-spatial dimension, not just literally "time". So cruise
and taxon both become slider axes here, same mechanism.

Why cell assignment is by STATION IDENTITY (grid_key), never by distance-
matching raw coordinates: checked the real data first (not assumed) --
half of all net-tow rows (38857/73792) sit more than 0.1 degree from their
OWN station's median recorded position, across 192 of the 207 distinct
stations. grid_key ("st50-ln76.7") is the nominal CalCOFI grid-design
position; the recorded lat/lon is the ship's real, operationally-variable
tow position, which drifted from the nominal point by a non-trivial amount
over 70+ years. So every occurrence is assigned to a cell purely by its
net's own grid_key string -- never by snapping its raw coordinates to the
nearest grid tick -- which guarantees exactly one station per cell
regardless of how jittery the underlying tows were. Each station's own
*displayed* position on the grid uses the median (robust to outliers) of
its real recorded positions.

Grid resolution (0.05 deg) was chosen and verified, not assumed: the
tightest real gap between any two of the 207 distinct stations' median
positions is 0.1876 deg (st35-ln86.7 <-> st30-ln86.7) -- 0.05 deg gives a
~3.75x safety margin. After snapping, this script explicitly verifies
zero two stations landed in the same cell and raises if that's ever not
true, rather than trusting the resolution choice alone.

Geographic extent: the file's own global attribute
(geospatial_coverage = 0.0-54.4N, 179.8-77.2W) is misleading -- it's
dragged by a handful of wild individual-TOW outlier positions, not real
station locations. Using the robust per-station median instead, the real
207-station extent is 19.8-48.1N, -133.7 to -107.4W (padded below).
Several of the stations outside the "modern core" box carry a "_hist"
grid_key suffix (e.g. st80-ln55_hist) -- the older, wider historical
station layout referenced by the earliest cruises. All 207 are included,
per explicit direction to capture the full historical grid, not just the
modern reduced one.

Sizing: the full logical array (676 cruises x 963 taxa x ~600 x ~560)
is ~219 billion cells -- far too large to ever materialize as a dense
array, but fine as an HDF5-chunked/compressed, lazily-allocated NetCDF4
variable: small chunks (1, 1, lat_chunk, lon_chunk) mean only chunks a
real occurrence actually touches ever get written, so real disk usage
tracks the ~480K occurrence rows (fewer, after same-cell rows get summed
together), not the logical size. Written via per-cell indexed assignment,
never as one dense in-memory array.

Output also carries cruise_label(cruise)/taxon_label(taxon)/
scientific_name_label(taxon) string lookup arrays in the same file, since
GeoLibre's own axis slider only shows the bare numeric index -- a next
step (not built here) can use these to show real labels alongside the
slider.

Input:  data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io/swfsc_ichthyo.nc
Output: data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io/swfsc_ichthyo_gridded_cruise_taxon.nc
"""

import os
import collections
from datetime import datetime, timezone

import netCDF4
import numpy as np

BASE_DIR = "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io"
SRC_PATH = os.path.join(BASE_DIR, "swfsc_ichthyo.nc")
OUT_PATH = os.path.join(BASE_DIR, "swfsc_ichthyo_gridded_cruise_taxon.nc")

GRID_RES_DEG = 0.05          # verified ~3.75x under the tightest real station gap (0.1876 deg)
GRID_PAD_DEG = 0.5           # padding beyond the real station extent
DOUBLE_FILL = 9.96920996838687e36
LAT_CHUNK = 32
LON_CHUNK = 32


def to_char_array(strings, strlen):
    """Fixed-width bytes -> (n, strlen) 'S1' char array for a NetCDF char
    variable. Replaces netCDF4.stringtochar, which raises AttributeError
    under netCDF4 1.7.4 + numpy 2.5 (calcofi env)."""
    fixed = np.asarray(strings).astype(f"S{strlen}")
    return fixed.view("S1").reshape(fixed.shape + (strlen,))


def main():
    src = netCDF4.Dataset(SRC_PATH, "r")
    occ = src.groups["occurrence"]
    net = src.groups["net"]
    net_n = src.dimensions["net_n"].size

    # ---- join occurrence -> net (same verified join as the flat-point script) ----
    parent_index_1based = occ.variables["parent_index"][:]
    parent_index_0based = np.asarray(parent_index_1based, dtype=np.int64) - 1
    bad = (parent_index_0based < 0) | (parent_index_0based >= net_n)
    if bad.any():
        raise ValueError(f"{int(bad.sum())} occurrence row(s) have an out-of-range parent_index")

    net_grid_key = netCDF4.chartostring(net.variables["grid_key"][:])
    net_cruise_key = netCDF4.chartostring(net.variables["cruise_key"][:])
    net_lat = np.asarray(net.variables["latitude"][:])
    net_lon = np.asarray(net.variables["longitude"][:])

    occ_grid_key = net_grid_key[parent_index_0based]
    occ_cruise_key = net_cruise_key[parent_index_0based]
    abundance = np.asarray(occ.variables["abundance"][:])
    taxon_key = netCDF4.chartostring(occ.variables["taxon_key"][:])
    scientific_name = netCDF4.chartostring(occ.variables["scientific_name"][:])

    n_total = len(abundance)
    valid = occ_grid_key != ""
    n_excluded = int((~valid).sum())
    print(f"{n_total} occurrence rows, excluding {n_excluded} with no station (blank grid_key)")

    # ---- canonical (median) station positions, from net's own rows (not just occurrence's) ----
    by_station = collections.defaultdict(list)
    for gk, la, lo in zip(net_grid_key, net_lat, net_lon):
        if gk == "" or np.ma.is_masked(la) or np.ma.is_masked(lo):
            continue
        by_station[gk].append((la, lo))
    station_canonical = {gk: np.median(np.array(pts), axis=0) for gk, pts in by_station.items()}
    print(f"{len(station_canonical)} distinct stations with a valid position")

    src.close()

    # ---- build the regular lat/lon grid ----
    stat_lats = np.array([v[0] for v in station_canonical.values()])
    stat_lons = np.array([v[1] for v in station_canonical.values()])
    lat_min = np.floor((stat_lats.min() - GRID_PAD_DEG) / GRID_RES_DEG) * GRID_RES_DEG
    lat_max = np.ceil((stat_lats.max() + GRID_PAD_DEG) / GRID_RES_DEG) * GRID_RES_DEG
    lon_min = np.floor((stat_lons.min() - GRID_PAD_DEG) / GRID_RES_DEG) * GRID_RES_DEG
    lon_max = np.ceil((stat_lons.max() + GRID_PAD_DEG) / GRID_RES_DEG) * GRID_RES_DEG

    lat_axis = np.arange(lat_min, lat_max + GRID_RES_DEG / 2, GRID_RES_DEG)
    lon_axis = np.arange(lon_min, lon_max + GRID_RES_DEG / 2, GRID_RES_DEG)
    print(f"grid: {len(lat_axis)} lat x {len(lon_axis)} lon ticks at {GRID_RES_DEG} deg "
          f"({lat_min:.2f}-{lat_max:.2f}N, {lon_min:.2f}-{lon_max:.2f}E)")

    # snap each station to its nearest tick, and VERIFY no two stations collide --
    # a hard requirement, checked explicitly rather than trusted from the resolution margin alone.
    station_cell = {}
    for gk, (la, lo) in station_canonical.items():
        lat_idx = int(np.argmin(np.abs(lat_axis - la)))
        lon_idx = int(np.argmin(np.abs(lon_axis - lo)))
        station_cell[gk] = (lat_idx, lon_idx)

    cell_to_stations = collections.defaultdict(list)
    for gk, cell in station_cell.items():
        cell_to_stations[cell].append(gk)
    collisions = {cell: gks for cell, gks in cell_to_stations.items() if len(gks) > 1}
    if collisions:
        raise ValueError(
            f"{len(collisions)} grid cell(s) have more than one station -- "
            f"refusing to proceed. GRID_RES_DEG={GRID_RES_DEG} is too coarse. "
            f"Examples: {dict(list(collisions.items())[:5])}"
        )
    print(f"Verified: all {len(station_cell)} stations map to distinct grid cells -- no collisions")

    # ---- build cruise and taxon axes ----
    cruise_list = sorted(set(occ_cruise_key[valid].tolist()))
    taxon_list = sorted(set(taxon_key[valid].tolist()))
    cruise_index = {c: i for i, c in enumerate(cruise_list)}
    taxon_index = {t: i for i, t in enumerate(taxon_list)}
    # one representative scientific_name per taxon_key, for the label lookup
    taxon_sciname = {}
    for tk, sn in zip(taxon_key, scientific_name):
        if tk not in taxon_sciname and tk != "":
            taxon_sciname[tk] = sn
    print(f"{len(cruise_list)} cruises, {len(taxon_list)} taxa")

    # ---- accumulate into a sparse dict first (pure numpy/dict, no file I/O yet) ----
    accum = {}
    n_summed = 0
    for i in range(n_total):
        if not valid[i]:
            continue
        a = abundance[i]
        if np.ma.is_masked(a):
            continue
        gk = occ_grid_key[i]
        cell = station_cell.get(gk)
        if cell is None:
            continue
        key = (cruise_index[occ_cruise_key[i]], taxon_index[taxon_key[i]], cell[0], cell[1])
        if key in accum:
            n_summed += 1
        accum[key] = accum.get(key, 0.0) + float(a)

    print(f"{len(accum)} distinct (cruise, taxon, station) cells to write "
          f"({n_summed} occurrence rows summed into an already-populated cell)")

    # ---- write the output file ----
    out = netCDF4.Dataset(OUT_PATH, "w", format="NETCDF4")
    out.Conventions = "CF-1.10, ACDD-1.3"
    out.title = "SWFSC Ichthyoplankton abundance, gridded by cruise and taxon"
    out.summary = (
        "Specimen abundance (count) from swfsc_ichthyo.nc, gridded onto a regular "
        f"{GRID_RES_DEG} deg lat/lon mesh with cruise and taxon as independent leading "
        "axes, for GeoLibre's slider-based NetCDF rendering. Cell assignment is by "
        "station identity (net.grid_key), not by distance-matching raw tow "
        "coordinates -- see this script's own module docstring for why."
    )
    out.source = "Derived from swfsc_ichthyo.nc (CalCOFI integrated database release v2026.09.06)"
    out.history = (
        f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}: created by "
        "build_gridded_abundance_by_cruise_and_taxon.py"
    )

    out.createDimension("lat", len(lat_axis))
    out.createDimension("lon", len(lon_axis))
    out.createDimension("cruise", len(cruise_list))
    out.createDimension("taxon", len(taxon_list))
    out.createDimension("cruise_strlen", 32)
    out.createDimension("taxon_strlen", 64)
    out.createDimension("sciname_strlen", 64)

    v_lat = out.createVariable("lat", "f8", ("lat",))
    v_lat.standard_name = "latitude"
    v_lat.units = "degrees_north"
    v_lat[:] = lat_axis

    v_lon = out.createVariable("lon", "f8", ("lon",))
    v_lon.standard_name = "longitude"
    v_lon.units = "degrees_east"
    v_lon[:] = lon_axis

    v_cruise = out.createVariable("cruise", "i4", ("cruise",))
    v_cruise.long_name = "cruise index -- see cruise_label for the real cruise_key"
    v_cruise[:] = np.arange(len(cruise_list))

    v_taxon = out.createVariable("taxon", "i4", ("taxon",))
    v_taxon.long_name = "taxon index -- see taxon_label/scientific_name_label for the real identity"
    v_taxon[:] = np.arange(len(taxon_list))

    v_cruise_label = out.createVariable("cruise_label", "S1", ("cruise", "cruise_strlen"))
    v_cruise_label.long_name = "cruise_key for each cruise index"
    v_cruise_label[:] = to_char_array(cruise_list, 32)

    v_taxon_label = out.createVariable("taxon_label", "S1", ("taxon", "taxon_strlen"))
    v_taxon_label.long_name = "taxon_key for each taxon index"
    v_taxon_label[:] = to_char_array(taxon_list, 64)

    v_sciname_label = out.createVariable("scientific_name_label", "S1", ("taxon", "sciname_strlen"))
    v_sciname_label.long_name = "scientific_name for each taxon index"
    v_sciname_label[:] = to_char_array([taxon_sciname.get(t, "") for t in taxon_list], 64)

    v_abund = out.createVariable(
        "abundance", "f8", ("cruise", "taxon", "lat", "lon"),
        fill_value=DOUBLE_FILL,
        chunksizes=(1, 1, min(LAT_CHUNK, len(lat_axis)), min(LON_CHUNK, len(lon_axis))),
        zlib=True, complevel=4,
    )
    v_abund.long_name = "Specimen count per net tow, summed per (cruise, taxon, station) cell"
    v_abund.units = "count"

    # sparse, per-cell indexed writes -- never a dense in-memory array
    for i, ((ci, ti, lai, loi), val) in enumerate(accum.items()):
        v_abund[ci, ti, lai, loi] = val
        if (i + 1) % 50000 == 0:
            print(f"  written {i + 1}/{len(accum)} cells")

    out.close()
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
