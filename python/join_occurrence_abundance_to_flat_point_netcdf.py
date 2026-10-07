# join_occurrence_abundance_to_flat_point_netcdf.py

"""
Flattens swfsc_ichthyo.nc's grouped, ragged-array structure into a single,
CF-compliant "point" feature-type NetCDF that GeoLibre (or any other
generic netCDF/GIS point-rendering tool) can actually read.

Why this is needed: the source file's occurrence group stores
`abundance(occurrence_n)` with no lat/lon/time attached at that level --
those live three groups up, in `net` (latitude, longitude, time), reached
only via `occurrence.parent_index`, a 1-based explicit link. This is a
documented, deliberate design (see the source file's own global attribute
`cf_scope`: "...CF defines no feature type for this nesting."), not a bug
-- each level stores its own properties exactly once rather than
duplicating lat/lon/time into every occurrence row. But it also means no
generic tool can auto-discover coordinates for `abundance`, which is
exactly why GeoLibre reports "No renderable (2-D or higher) variables
found in the file."

What this script does:
  1. Reads occurrence (abundance, taxon_key, scientific_name, life_stage,
     depth_min_m, parent_index) and net (latitude, longitude, time).
  2. Resolves parent_index (1-based into net) to pull lat/lon/time onto
     every occurrence row -- a single-level join, since net already has
     geo-coordinates directly (no need to go further up through tow/site).
  3. Writes a flat NetCDF4 file with lat/lon/time as real CF coordinate
     variables on a single `occurrence` dimension, `abundance` carrying a
     `coordinates` attribute pointing at them, and global attribute
     `featureType = "point"` -- the standard CF Discrete Sampling Geometry
     point encoding that GeoLibre/QGIS/Panoply all know how to render
     without any custom logic.
  4. taxon_key/scientific_name/life_stage/depth_min_m ride along as plain
     ancillary variables on the same dimension, so taxa can be filtered/
     colored by directly in GeoLibre's own UI rather than needing one
     output file per taxon.

Input:  data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io/swfsc_ichthyo.nc
Output: data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io/swfsc_ichthyo_converted.nc
"""

import os
from datetime import datetime, timezone

import netCDF4
import numpy as np

BASE_DIR = "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io"
SRC_PATH = os.path.join(BASE_DIR, "swfsc_ichthyo.nc")
OUT_PATH = os.path.join(BASE_DIR, "swfsc_ichthyo_converted.nc")

# netCDF4's own default double/int fill values -- matches what the source
# file already uses for these exact variables (abundance, depth_min_m,
# order_occ/count), so the output stays consistent with the source rather
# than introducing a second, different magic number.
DOUBLE_FILL = 9.96920996838687e36
INT_FILL = -2147483647


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

    occurrence_n = src.dimensions["occurrence_n"].size
    net_n = src.dimensions["net_n"].size

    # parent_index is documented as "1-based index into the net group" --
    # convert to 0-based for numpy indexing, and validate range before
    # trusting it for a join (a silently-wrong index would put the wrong
    # tow's position under a real abundance count, a much worse failure
    # mode than just erroring out here).
    parent_index_1based = occ.variables["parent_index"][:]
    if np.ma.is_masked(parent_index_1based):
        missing = int(np.ma.getmaskarray(parent_index_1based).sum())
        if missing:
            raise ValueError(
                f"{missing} occurrence row(s) have no parent_index (net link) at all -- "
                "can't join these to a position. Not expected for this file; investigate "
                "before proceeding rather than silently dropping rows."
            )
    parent_index_0based = np.asarray(parent_index_1based, dtype=np.int64) - 1

    bad = (parent_index_0based < 0) | (parent_index_0based >= net_n)
    if bad.any():
        raise ValueError(
            f"{int(bad.sum())} occurrence row(s) have a parent_index outside "
            f"the valid net range [1, {net_n}] -- refusing to join against an "
            "out-of-range position rather than guessing."
        )

    # The actual join: every occurrence row picks up its net event's own
    # lat/lon/time by position.
    net_lat = net.variables["latitude"][:]
    net_lon = net.variables["longitude"][:]
    net_time = net.variables["time"][:]
    lat = net_lat[parent_index_0based]
    lon = net_lon[parent_index_0based]
    time = net_time[parent_index_0based]

    abundance = occ.variables["abundance"][:]
    depth_min_m = occ.variables["depth_min_m"][:]
    taxon_key = netCDF4.chartostring(occ.variables["taxon_key"][:])
    scientific_name = netCDF4.chartostring(occ.variables["scientific_name"][:])
    life_stage = netCDF4.chartostring(occ.variables["life_stage"][:])

    # Capture these before closing src -- occ/net's own Variable objects
    # (and their .shape) become invalid once the parent Dataset is closed.
    taxon_strlen = occ.variables["taxon_key"].shape[1]
    sciname_strlen = occ.variables["scientific_name"].shape[1]
    stage_strlen = occ.variables["life_stage"].shape[1]

    src.close()

    print(f"Joined {occurrence_n} occurrence rows to net-level lat/lon/time "
          f"(net group has {net_n} rows)")
    print(f"abundance range: {np.nanmin(abundance):.1f} to {np.nanmax(abundance):.1f}")
    print(f"lat range: {np.nanmin(lat):.3f} to {np.nanmax(lat):.3f}")
    print(f"lon range: {np.nanmin(lon):.3f} to {np.nanmax(lon):.3f}")
    print(f"distinct taxa: {len(set(taxon_key.tolist()))}")

    # ---- write the flat, CF point-featureType output ----
    out = netCDF4.Dataset(OUT_PATH, "w", format="NETCDF4")

    out.Conventions = "CF-1.10, ACDD-1.3"
    out.featureType = "point"
    out.title = "SWFSC Ichthyoplankton abundance by taxon -- flat point form"
    out.summary = (
        "Per-occurrence specimen abundance (count), with latitude/longitude/time "
        "resolved from the source file's net group via its parent_index link, so "
        "this file is directly renderable by generic netCDF/GIS point tools "
        "(GeoLibre, QGIS, Panoply) without needing to understand the source's "
        "nested site/tow/net/occurrence group hierarchy."
    )
    out.source = (
        "Derived from swfsc_ichthyo.nc (CalCOFI integrated database release "
        "v2026.09.06) by joining the occurrence group's own abundance/taxon_key/ "
        "scientific_name/life_stage/depth_min_m against the net group's "
        "latitude/longitude/time, via occurrence.parent_index (1-based, into net)."
    )
    out.history = (
        f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}: created by "
        "join_occurrence_abundance_to_flat_point_netcdf.py"
    )

    out.createDimension("occurrence", occurrence_n)
    out.createDimension("taxon_strlen", taxon_strlen)
    out.createDimension("sciname_strlen", sciname_strlen)
    out.createDimension("stage_strlen", stage_strlen)

    v_lat = out.createVariable("lat", "f8", ("occurrence",), fill_value=DOUBLE_FILL)
    v_lat.standard_name = "latitude"
    v_lat.long_name = "latitude"
    v_lat.units = "degrees_north"
    v_lat[:] = lat

    v_lon = out.createVariable("lon", "f8", ("occurrence",), fill_value=DOUBLE_FILL)
    v_lon.standard_name = "longitude"
    v_lon.long_name = "longitude"
    v_lon.units = "degrees_east"
    v_lon[:] = lon

    v_time = out.createVariable("time", "f8", ("occurrence",), fill_value=DOUBLE_FILL)
    v_time.standard_name = "time"
    v_time.long_name = "time"
    v_time.units = "seconds since 1970-01-01T00:00:00Z"
    v_time.calendar = "standard"
    v_time[:] = time

    v_abund = out.createVariable("abundance", "f8", ("occurrence",), fill_value=DOUBLE_FILL)
    v_abund.long_name = "Specimen count per net tow (headline occurrence; standardize via std_haul_factor in the source file's net group)"
    v_abund.units = "count"
    v_abund.coordinates = "time lat lon"
    v_abund[:] = abundance

    v_depth = out.createVariable("depth_min_m", "f8", ("occurrence",), fill_value=DOUBLE_FILL)
    v_depth.long_name = "depth min m"
    v_depth.coordinates = "time lat lon"
    v_depth[:] = depth_min_m

    v_taxon = out.createVariable("taxon_key", "S1", ("occurrence", "taxon_strlen"))
    v_taxon.long_name = "taxon key"
    v_taxon.coordinates = "time lat lon"
    v_taxon[:] = to_char_array(taxon_key, taxon_strlen)

    v_sci = out.createVariable("scientific_name", "S1", ("occurrence", "sciname_strlen"))
    v_sci.long_name = "scientific name"
    v_sci.coordinates = "time lat lon"
    v_sci[:] = to_char_array(scientific_name, sciname_strlen)

    v_stage = out.createVariable("life_stage", "S1", ("occurrence", "stage_strlen"))
    v_stage.long_name = "life stage"
    v_stage.coordinates = "time lat lon"
    v_stage[:] = to_char_array(life_stage, stage_strlen)

    out.close()

    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
