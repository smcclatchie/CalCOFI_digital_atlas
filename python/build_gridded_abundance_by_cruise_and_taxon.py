# build_gridded_abundance_by_cruise_and_taxon.py

"""
Converts the CalCOFI.io swfsc_ichthyo.nc download (grouped, ragged
site -> tow -> net -> occurrence hierarchy, which generic NetCDF tools
cannot render) into a filtered, standardised, CF-1.10 gridded NetCDF:

    abundance_per_10m2 (cruise, sampling, net_type, life_stage, taxon, lat, lon)
        oblique and vertical nets (C1, CB, CV, PV): eggs or larvae per 10 m2
    abundance_per_100m3(cruise, sampling, life_stage, taxon, lat, lon)
        Manta surface net (MT): eggs or larvae per 100 m3
    tows               (cruise, sampling, net_type, lat, lon)   oblique/vertical tows
    manta_tows         (cruise, sampling, lat, lon)             Manta tows

restricted to CalCOFI cruises and CalCOFI grid stations. Every non-spatial
dimension is a leading axis, so slider-based viewers (GeoLibre, ncview,
Panoply) can step through cruise, net type, life stage and taxon over a
shared regular lat/lon grid.

Usage (defaults are this project's own paths and release):
    python build_gridded_abundance_by_cruise_and_taxon.py \
        --src swfsc_ichthyo.nc --out swfsc_ichthyo_calcofi_gridded.nc \
        --release v2026.10.06

Dependencies: numpy, pandas, netCDF4, duckdb (reads the public CalCOFI.io
Parquet `cruise` table; no credentials needed). Internet access is also
needed for NOAA ERDDAP (Pairovet net sides, below).

CalCOFI selection (select_calcofi_nets) -- applied per net tow:
  1. Hydrographic CalCOFI cruise: the cruise also contributed CalCOFI bottle
     or CTD data (`cruise.cruise_key_datasets` contains calcofi_bottle or
     calcofi_ctd-cast). Covers CTD-only modern cruises too.
  2. Restored CalCOFI-ship cruise: no bottle/CTD data in the database, but
     the ship ran at least one cruise of type 1 AND >= 90% of the cruise's
     tows are on grid stations inside the domain. Recovers early (1950s)
     CalCOFI cruises whose hydrography is not in the database. Robust to
     the threshold: 91/90/89 cruises at 80/90/100%; the CalCOFI-ship
     cruises left out have a median of 7% in-domain tows.
  3. Domain: 18-42N (California-Oregon border), 140-105W, outside the Gulf
     of California.
A tow is kept if its cruise passes 1 or 2, the tow is in the domain, and
it has a grid_key (sits on a CalCOFI grid station). Release v2026.10.06:
609 cruises (519 + 90), 68,675 of 76,512 net tows.

Standardisation (verified against NOAA ERDDAP erdCalCOFIlrvcnt, exact
match including partially sorted samples):
  per 10 m2  = count * std_haul_factor / prop_sorted   (oblique, vertical)
  per 100 m3 = 100 * count / volume_sampled / prop_sorted   (Manta)
NOAA publishes no per-10 m2 value for Manta tows, so they are kept in their
own variable rather than mixed into abundance_per_10m2. Invertebrate
occurrences (life_stage "invert") are excluded; eggs and larvae are kept
separate on the life_stage axis, and net types on the net_type axis.

Sampling type: high-resolution and special-purpose net surveys (e.g. the
1980-86 anchovy egg surveys, 1987-91 Yellowfin and 1994 McArthur Pairovet
surveys, 2002-04 cowcod conservation area sampling, spring sardine
surveys) are kept apart from standard CalCOFI sampling on the `sampling`
axis and never averaged with it. A tow is "standard" when its line/station
has CalCOFI bottle or CTD casts (CalCOFI.io `sample` table) AND either its
cruise contributed CalCOFI bottle/CTD data, or (restored cruises with no
hydrographic data) >= 80% of the cruise's tows are at such stations. Every
other tow is "special survey". Hydrographic cruises are judged tow by tow
because several CalCOFI cruises (e.g. the April 2003-2010 spring cruises)
added survey stations to the standard pattern.

Paired nets: the protocol uses the right-hand (starboard) net of paired
samplers. The CalCOFI.io file has no net-side field. Bongo (CB) tows carry
one net each, and matching against NOAA ERDDAP (which records
net_location) shows it is the starboard net. 1,044 Pairovet (PV) tows
carry both nets; each is matched to ERDDAP's erdCalCOFIlrvcnt/eggcnt
records by cruise, line, station and standard haul factor, and the port
net is dropped. Pairs with no ERDDAP record for either net (no eggs or
larvae recorded there) keep both nets, which are then averaged. Adding a
net-side field to the CalCOFI.io net table would remove the need for this
lookup.

Stations: a station is the occupied CalCOFI line and station (`site_key`,
e.g. "091.2 032.4"), placed at the median recorded position of its kept
tows. `grid_key` is NOT used for placement: it is an area that can contain
many distinct stations (up to 260 km apart in the historical layout).
Stations are snapped to a regular GRID_RES_DEG (0.02 deg, ~2 km) grid;
the few distinct line/stations whose positions fall in the same ~2 km cell
(some have identical recorded positions) are treated as one location.

Repeat tows: when one cell was towed more than once with the same net type
on one cruise, the cell value is the MEAN standardised abundance over those
tows, with tows that did not catch the taxon counting as zero (i.e. sum of
densities / number of tows). The `tows`/`manta_tows` variables give that
tow count. A fill value in an abundance variable means the taxon was not
caught there; where the matching tow count is > 0 that is a sampled zero,
and where it is fill the cell was not sampled by that net type.

Sizing: the logical arrays are large but are written as chunked,
compressed NetCDF4 variables with 32 x 32 spatial chunks and length-1
chunks on every other axis, so only chunks holding real data are stored.
"""

import argparse
import os
from datetime import datetime, timezone

import duckdb
import netCDF4
import numpy as np
import pandas as pd

BASE_DIR = "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io"
DEFAULT_SRC = os.path.join(BASE_DIR, "swfsc_ichthyo.nc")
DEFAULT_OUT = os.path.join(BASE_DIR, "swfsc_ichthyo_calcofi_gridded.nc")
DEFAULT_RELEASE = "v2026.10.06"
PARQUET_URL = "https://storage.googleapis.com/calcofi-db/ducklake/releases/{release}/parquet/{table}.parquet"

# ---- CalCOFI selection ----
HYDRO_DATASETS = ("calcofi_bottle", "calcofi_ctd-cast")
RESTORE_MIN_IN_DOMAIN = 0.9
DOMAIN_LON = (-140.0, -105.0)
DOMAIN_LAT = (18.0, 42.0)          # 42N = California-Oregon border
# Gulf of California: east of the Baja peninsula's Gulf coast, closed across
# the mouth from Cabo San Lucas to Cabo Corrientes (IHO southern limit).
GULF_OF_CALIFORNIA = np.array([
    (-115.0, 33.0), (-115.0, 32.0), (-114.6, 30.0), (-113.2, 28.0),
    (-111.8, 26.0), (-110.3, 24.0), (-109.9, 22.9), (-105.7, 20.4),
    (-104.0, 20.4), (-104.0, 33.0),
])

# ---- standard CalCOFI sampling vs special surveys ----
SAMPLING = ["standard", "special survey"]
STANDARD_CRUISE_MIN_HYDRO_STATIONS = 0.8

# ---- standardisation ----
LIFE_STAGES = ["egg", "larva"]                 # "invert" excluded
AREAL_NET_TYPES = {                            # standardised per 10 m2
    "C1": "C1: CalCOFI 1 m ring net (oblique)",
    "CB": "CB: CalCOFI bongo (oblique)",
    "CV": "CV: CalVET (vertical)",
    "PV": "PV: Pairovet (vertical)",
}
MANTA = "MT"                                   # standardised per 100 m3
PAIROVET = "PV"
ERDDAP_PV_SIDES = [
    "https://coastwatch.pfeg.noaa.gov/erddap/tabledap/{dataset}.csv?cruise,line,station,net_location,"
    "standard_haul_factor&tow_type=%22PV%22&distinct()".format(dataset=d)
    for d in ("erdCalCOFIlrvcnt", "erdCalCOFIeggcnt")
]

# ---- gridding ----
GRID_RES_DEG = 0.02          # ~2 km
GRID_PAD_DEG = 0.5           # padding beyond the kept station extent
DOUBLE_FILL = 9.96920996838687e36
INT_FILL = -2147483647
LAT_CHUNK = 32
LON_CHUNK = 32


def to_char_array(strings, strlen):
    """Fixed-width bytes -> (n, strlen) 'S1' char array for a NetCDF char
    variable. Replaces netCDF4.stringtochar, which raises AttributeError
    under netCDF4 1.7.4 + numpy 2.5."""
    fixed = np.asarray(strings).astype(f"S{strlen}")
    return fixed.view("S1").reshape(fixed.shape + (strlen,))


def points_in_polygon(lon, lat, polygon):
    """Even-odd ray-casting test; NaN positions are outside."""
    inside = np.zeros(len(lon), dtype=bool)
    x1, y1 = polygon[-1]
    for x2, y2 in polygon:
        crosses = ((y1 > lat) != (y2 > lat)) & (lon < (x2 - x1) * (lat - y1) / (y2 - y1 + 1e-300) + x1)
        inside ^= crosses
        x1, y1 = x2, y2
    return inside


def read_cruise_table(release):
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs;")
    url = PARQUET_URL.format(release=release, table="cruise")
    return con.sql(f"SELECT cruise_key, cruise_key_datasets, ship_name FROM read_parquet('{url}')").df()


def read_hydro_sites(release):
    """Every line/station (site_key) where CalCOFI made a bottle or CTD cast."""
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs;")
    url = PARQUET_URL.format(release=release, table="sample")
    datasets = ", ".join(f"'{d}'" for d in HYDRO_DATASETS)
    return set(con.sql(f"""SELECT DISTINCT site_key FROM read_parquet('{url}')
                           WHERE dataset_key IN ({datasets}) AND site_key IS NOT NULL""").df().site_key)


def classify_sampling(tows, hydro_sites, per_cruise):
    """'standard' or 'special survey' per tow. tows needs cruise_key and
    site_key; per_cruise is select_calcofi_nets()'s per-cruise table.

    A tow is 'standard' when its line/station has CalCOFI bottle/CTD casts
    AND either its cruise contributed CalCOFI bottle/CTD data (judged tow by
    tow, so a CalCOFI cruise's extra survey stations become special survey
    while its standard stations stay standard) or, for restored cruises
    with no hydrographic data, >= STANDARD_CRUISE_MIN_HYDRO_STATIONS of the
    cruise's tows are at such stations."""
    at_hydro = tows.site_key.isin(hydro_sites)
    hydro_cruise = tows.cruise_key.map(per_cruise.hydro_cruise).fillna(False).astype(bool)
    standard_share = at_hydro.groupby(tows.cruise_key).transform("mean") >= STANDARD_CRUISE_MIN_HYDRO_STATIONS
    standard = at_hydro & (hydro_cruise | standard_share)
    return pd.Series(np.where(standard, SAMPLING[0], SAMPLING[1]), index=tows.index)


def pairovet_port_nets(nets):
    """Index of port-side nets to drop from two-net Pairovet tows, using
    NOAA ERDDAP's net_location. nets needs cruise_key, site_key, tow_type,
    tow (parent tow index) and std_haul_factor. Returns (port net index,
    number of two-net tows left unresolved)."""
    erddap = pd.concat([pd.read_csv(url, skiprows=[1]) for url in ERDDAP_PV_SIDES]).drop_duplicates()
    erddap["shf"] = erddap.standard_haul_factor.round(2)
    pairs = nets[nets.tow_type == PAIROVET]
    pairs = pairs[pairs.groupby("tow").tow.transform("size") == 2].copy()
    pairs["cruise"] = pairs.cruise_key.str[:4].astype(int) * 100 + pairs.cruise_key.str[5:7].astype(int)
    pairs["line"] = pairs.site_key.str[:5].astype(float)
    pairs["station"] = pairs.site_key.str[6:].astype(float)
    pairs["shf"] = pairs.std_haul_factor.round(2)
    pairs["side"] = (pairs.reset_index()
                          .merge(erddap[["cruise", "line", "station", "shf", "net_location"]],
                                 on=["cruise", "line", "station", "shf"], how="left")
                          .groupby("index").net_location
                          .agg(lambda x: "".join(sorted(set(x.dropna())))))
    port = []
    unresolved = 0
    for _, tow in pairs.groupby("tow"):
        sides = tow.side.tolist()
        if "S" in sides:                     # keep the starboard net
            port += tow.index[[x != "S" for x in sides]].tolist()
        elif "P" in sides:                   # port identified: the other net is starboard
            port += tow.index[[x == "P" for x in sides]].tolist()
        else:
            unresolved += 1
    return pd.Index(port), unresolved


def select_calcofi_nets(nets, cruise):
    """nets: DataFrame with cruise_key, grid_key, lat, lon (one row per net
    tow). cruise: the CalCOFI.io cruise table. Returns (selected boolean
    Series aligned to nets, per-cruise DataFrame with the decision and
    reason)."""
    d = nets[["cruise_key", "grid_key", "lat", "lon"]].merge(cruise, on="cruise_key", how="left")
    d.index = nets.index
    lat, lon = d.lat.to_numpy(), d.lon.to_numpy()
    in_domain = (
        d.lat.between(*DOMAIN_LAT) & d.lon.between(*DOMAIN_LON)
        & ~points_in_polygon(lon, lat, GULF_OF_CALIFORNIA)
    )
    on_grid = d.grid_key != ""
    hydro = d.cruise_key_datasets.fillna("").str.contains("|".join(HYDRO_DATASETS))
    calcofi_ships = set(d.loc[hydro, "ship_name"].dropna())

    per_cruise = (
        d.assign(hydro=hydro, in_domain_on_grid=in_domain & on_grid)
         .groupby("cruise_key")
         .agg(ship_name=("ship_name", "first"), datasets=("cruise_key_datasets", "first"),
              hydro_cruise=("hydro", "first"), tows=("lat", "size"),
              frac_in_domain_on_grid=("in_domain_on_grid", "mean"))
    )
    restored = (~per_cruise.hydro_cruise & per_cruise.ship_name.isin(calcofi_ships)
                & (per_cruise.frac_in_domain_on_grid >= RESTORE_MIN_IN_DOMAIN))
    per_cruise["reason"] = np.select(
        [per_cruise.hydro_cruise, restored],
        ["hydrographic CalCOFI cruise", "restored: CalCOFI ship, on-grid in domain"],
        default="dropped")

    cruise_ok = d.cruise_key.map(per_cruise.reason) != "dropped"
    return cruise_ok & in_domain & on_grid, per_cruise


def chunks_for(shape):
    """Length-1 chunks on every leading axis, 32 x 32 spatial chunks."""
    return tuple([1] * (len(shape) - 2) + [min(LAT_CHUNK, shape[-2]), min(LON_CHUNK, shape[-1])])


def write_sparse(var, table, index_cols, value_col):
    """Per-cell indexed writes from a long table -- never a dense array."""
    for row in table[index_cols + [value_col]].itertuples(index=False):
        var[tuple(row[:-1])] = row[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--src", default=DEFAULT_SRC, help="CalCOFI.io swfsc_ichthyo.nc")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output gridded NetCDF")
    ap.add_argument("--release", default=DEFAULT_RELEASE, help="CalCOFI.io release for the cruise table")
    args = ap.parse_args()

    src = netCDF4.Dataset(args.src, "r")
    occ = src.groups["occurrence"]
    net = src.groups["net"]
    net_n = src.dimensions["net_n"].size
    src_release = getattr(src, "db_release", "unknown")

    def num(var):
        return np.ma.filled(var[:].astype("f8"), np.nan)

    nets = pd.DataFrame({
        "cruise_key": netCDF4.chartostring(net.variables["cruise_key"][:]),
        "grid_key": netCDF4.chartostring(net.variables["grid_key"][:]),
        "tow_type": netCDF4.chartostring(net.variables["tow_type"][:]),
        "site_key": netCDF4.chartostring(net.variables["site_key"][:]),
        "tow": np.asarray(net.variables["parent_index"][:], dtype=np.int64) - 1,
        "lat": num(net.variables["latitude"]),
        "lon": num(net.variables["longitude"]),
        "std_haul_factor": num(net.variables["std_haul_factor"]),
        "prop_sorted": num(net.variables["prop_sorted"]),
        "volume_sampled": num(net.variables["volume_sampled"]),
    })

    # ---- CalCOFI selection, per net tow ----
    selected, per_cruise = select_calcofi_nets(nets, read_cruise_table(args.release))
    n_cruises = int((per_cruise.reason != "dropped").sum())
    print(f"selected {int(selected.sum())} of {net_n} net tows; {n_cruises} cruises")

    # Paired Pairovet nets: keep the starboard (right-hand) net.
    port, unresolved = pairovet_port_nets(nets)
    print(f"Pairovet: dropping {len(port)} port nets; {unresolved} two-net tows have no ERDDAP "
          f"side record and keep both nets (averaged)")
    selected = selected & ~nets.index.isin(port)

    # Per-net multiplier from raw count to standardised abundance.
    nets["areal"] = nets.tow_type.isin(list(AREAL_NET_TYPES))
    nets["manta"] = nets.tow_type == MANTA
    nets["factor"] = np.where(
        nets.areal, nets.std_haul_factor / nets.prop_sorted,
        np.where(nets.manta, 100.0 / (nets.volume_sampled * nets.prop_sorted), np.nan))
    usable = selected & (nets.areal | nets.manta) & np.isfinite(nets.factor) & (nets.factor > 0)
    print(f"usable tows: {int(usable.sum())} "
          f"({nets[usable].tow_type.value_counts().to_dict()}); "
          f"{int((selected & ~usable).sum())} selected tows lack a valid standardisation factor or net type")
    kept = nets[usable].copy()

    # Standard CalCOFI sampling vs high-resolution / special surveys.
    kept["sampling"] = classify_sampling(kept, read_hydro_sites(args.release), per_cruise)
    survey_cruises = kept.groupby("cruise_key").sampling.apply(lambda s: (s == SAMPLING[0]).sum() == 0)
    print(f"sampling: {kept.sampling.value_counts().to_dict()} tows; "
          f"{int(survey_cruises.sum())} cruises are entirely special survey")

    # ---- occurrences ----
    occ_df = pd.DataFrame({
        "net": np.asarray(occ.variables["parent_index"][:], dtype=np.int64) - 1,   # 1-based in source
        "count": num(occ.variables["abundance"]),
        "taxon_key": netCDF4.chartostring(occ.variables["taxon_key"][:]),
        "scientific_name": netCDF4.chartostring(occ.variables["scientific_name"][:]),
        "life_stage": netCDF4.chartostring(occ.variables["life_stage"][:]),
    })
    src.close()
    if ((occ_df.net < 0) | (occ_df.net >= net_n)).any():
        raise ValueError("occurrence rows with an out-of-range parent_index")
    print(f"{len(occ_df)} occurrence rows; life stages {occ_df.life_stage.value_counts().to_dict()}")
    occ_df = occ_df[occ_df.net.isin(kept.index) & occ_df.life_stage.isin(LIFE_STAGES)
                    & occ_df["count"].notna() & (occ_df.taxon_key != "")]
    occ_df = occ_df.join(kept[["cruise_key", "site_key", "tow_type", "sampling", "factor"]], on="net")
    occ_df["density"] = occ_df["count"] * occ_df.factor
    print(f"keeping {len(occ_df)} egg/larva occurrence rows")

    # ---- station (occupied line/station) positions: median of kept tows ----
    station_pos = kept.groupby("site_key")[["lat", "lon"]].median()

    # ---- regular lat/lon grid ----
    lat_min = np.floor((station_pos.lat.min() - GRID_PAD_DEG) / GRID_RES_DEG) * GRID_RES_DEG
    lat_max = np.ceil((station_pos.lat.max() + GRID_PAD_DEG) / GRID_RES_DEG) * GRID_RES_DEG
    lon_min = np.floor((station_pos.lon.min() - GRID_PAD_DEG) / GRID_RES_DEG) * GRID_RES_DEG
    lon_max = np.ceil((station_pos.lon.max() + GRID_PAD_DEG) / GRID_RES_DEG) * GRID_RES_DEG
    lat_axis = np.round(np.arange(lat_min, lat_max + GRID_RES_DEG / 2, GRID_RES_DEG), 4)
    lon_axis = np.round(np.arange(lon_min, lon_max + GRID_RES_DEG / 2, GRID_RES_DEG), 4)
    print(f"grid: {len(lat_axis)} lat x {len(lon_axis)} lon at {GRID_RES_DEG} deg "
          f"({lat_min:.2f} to {lat_max:.2f}N, {lon_min:.2f} to {lon_max:.2f}E), {len(station_pos)} stations")

    station_pos["lat_i"] = [int(np.argmin(np.abs(lat_axis - v))) for v in station_pos.lat]
    station_pos["lon_i"] = [int(np.argmin(np.abs(lon_axis - v))) for v in station_pos.lon]
    shared = station_pos.duplicated(["lat_i", "lon_i"], keep=False)
    print(f"{len(station_pos)} line/stations in {len(station_pos.drop_duplicates(['lat_i', 'lon_i']))} cells; "
          f"{int(shared.sum())} line/stations share a ~{GRID_RES_DEG} deg cell with another")

    # ---- axes ----
    cruise_list = sorted(kept.cruise_key.unique())
    net_list = [t for t in AREAL_NET_TYPES if (kept.tow_type == t).any()]
    taxon_list = sorted(occ_df.taxon_key.unique())
    taxon_sciname = occ_df.drop_duplicates("taxon_key").set_index("taxon_key").scientific_name
    index = {
        "cruise_i": {c: i for i, c in enumerate(cruise_list)},
        "sampling_i": {s: i for i, s in enumerate(SAMPLING)},
        "net_i": {t: i for i, t in enumerate(net_list)},
        "stage_i": {s: i for i, s in enumerate(LIFE_STAGES)},
        "taxon_i": {t: i for i, t in enumerate(taxon_list)},
    }
    print(f"{len(cruise_list)} cruises, net types {net_list} + {MANTA}, {len(taxon_list)} taxa")

    def add_indices(df):
        df = df.join(station_pos[["lat_i", "lon_i"]], on="site_key")
        df["cruise_i"] = df.cruise_key.map(index["cruise_i"])
        df["sampling_i"] = df.sampling.map(index["sampling_i"])
        if "tow_type" in df:
            df["net_i"] = df.tow_type.map(index["net_i"])
        if "life_stage" in df:
            df["stage_i"] = df.life_stage.map(index["stage_i"])
        if "taxon_key" in df:
            df["taxon_i"] = df.taxon_key.map(index["taxon_i"])
        return df

    # Tow (net) counts per (cruise, sampling, net type, cell) -- the denominator
    # for the mean. net_i = -1 marks Manta tows, which have their own variables.
    cell_keys = ["cruise_i", "sampling_i", "net_i", "lat_i", "lon_i"]
    kept_ix = add_indices(kept)
    kept_ix["net_i"] = kept_ix.net_i.fillna(-1).astype(int)
    tow_counts = kept_ix.groupby(cell_keys).size().rename("n_tows").reset_index()
    print(f"repeat tows: {int((tow_counts.n_tows > 1).sum())} of {len(tow_counts)} "
          f"(cruise, sampling, net type, cell) combinations have more than one tow")

    # Mean density over tows = summed density / number of tows (absences count as zero).
    occ_ix = add_indices(occ_df)
    occ_ix["net_i"] = occ_ix.net_i.fillna(-1).astype(int)
    cells = (occ_ix.groupby(cell_keys + ["stage_i", "taxon_i"]).density.sum()
                   .rename("density_sum").reset_index()
                   .merge(tow_counts, on=cell_keys))
    cells["mean_density"] = cells.density_sum / cells.n_tows
    areal_cells = cells[cells.net_i >= 0]
    manta_cells = cells[cells.net_i < 0]
    areal_tows = tow_counts[tow_counts.net_i >= 0]
    manta_tows = tow_counts[tow_counts.net_i < 0]
    print(f"{len(areal_cells)} per-10 m2 cells, {len(manta_cells)} per-100 m3 (Manta) cells")

    # ---- write ----
    out = netCDF4.Dataset(args.out, "w", format="NETCDF4")
    out.Conventions = "CF-1.10, ACDD-1.3"
    out.title = ("CalCOFI ichthyoplankton abundance, gridded by cruise, sampling type, net type, "
                 "life stage and taxon")
    out.summary = (
        "Fish egg and larva abundance from the CalCOFI.io swfsc_ichthyo download, restricted "
        "to CalCOFI cruises and CalCOFI grid stations (south of 42N, outside the Gulf of "
        "California) and gridded onto a regular "
        f"{GRID_RES_DEG} deg lat/lon mesh. Oblique and vertical tows are standardised to "
        "eggs or larvae per 10 m2 (count * std_haul_factor / prop_sorted); Manta surface tows "
        "to eggs or larvae per 100 m3 (100 * count / volume_sampled / prop_sorted). Repeat "
        "tows of one net type at one station on one cruise are averaged, absences counting "
        "as zero. Paired samplers use the starboard (right-hand) net. Stations are the "
        "occupied CalCOFI line/station, placed at the median position of their tows. "
        "Standard CalCOFI sampling and high-resolution / special surveys are kept apart on "
        "the sampling axis and never averaged together."
    )
    out.source = (f"swfsc_ichthyo.nc (CalCOFI.io release {src_release}); cruise selection "
                  f"from the CalCOFI.io cruise table (release {args.release})")
    out.references = ("Smith, P.E. and Richardson, S.L. (1977) Standard techniques for pelagic fish "
                      "egg and larva surveys. FAO Fish. Tech. Pap. 175; https://calcofi.io")
    out.history = (f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}: created by "
                   "build_gridded_abundance_by_cruise_and_taxon.py")
    out.geospatial_lat_min, out.geospatial_lat_max = float(lat_axis[0]), float(lat_axis[-1])
    out.geospatial_lon_min, out.geospatial_lon_max = float(lon_axis[0]), float(lon_axis[-1])

    for name, size in [("lat", len(lat_axis)), ("lon", len(lon_axis)), ("cruise", len(cruise_list)),
                       ("sampling", len(SAMPLING)), ("sampling_strlen", 16), ("net_type", len(net_list)), ("life_stage", len(LIFE_STAGES)),
                       ("taxon", len(taxon_list)), ("cruise_strlen", 32), ("net_strlen", 48),
                       ("stage_strlen", 8), ("taxon_strlen", 64), ("sciname_strlen", 64)]:
        out.createDimension(name, size)

    v = out.createVariable("lat", "f8", ("lat",))
    v.standard_name, v.units = "latitude", "degrees_north"
    v[:] = lat_axis
    v = out.createVariable("lon", "f8", ("lon",))
    v.standard_name, v.units = "longitude", "degrees_east"
    v[:] = lon_axis

    for dim, labels, strlen_dim, strlen, label_name, desc in [
        ("cruise", cruise_list, "cruise_strlen", 32, "cruise_label", "cruise_key"),
        ("sampling", SAMPLING, "sampling_strlen", 16, "sampling_label", "sampling type"),
        ("net_type", [AREAL_NET_TYPES[t] for t in net_list], "net_strlen", 48, "net_type_label", "net type"),
        ("life_stage", LIFE_STAGES, "stage_strlen", 8, "life_stage_label", "life stage"),
        ("taxon", taxon_list, "taxon_strlen", 64, "taxon_label", "taxon_key"),
    ]:
        v = out.createVariable(dim, "i4", (dim,))
        v.long_name = f"{dim} index -- see {label_name}"
        v[:] = np.arange(len(labels))
        v = out.createVariable(label_name, "S1", (dim, strlen_dim))
        v.long_name = f"{desc} for each {dim} index"
        v[:] = to_char_array(labels, strlen)
    v = out.createVariable("scientific_name_label", "S1", ("taxon", "sciname_strlen"))
    v.long_name = "scientific_name for each taxon index"
    v[:] = to_char_array([taxon_sciname.get(t, "") for t in taxon_list], 64)

    areal_dims = ("cruise", "sampling", "net_type", "life_stage", "taxon", "lat", "lon")
    v = out.createVariable("abundance_per_10m2", "f8", areal_dims, fill_value=DOUBLE_FILL,
                           chunksizes=chunks_for([out.dimensions[d].size for d in areal_dims]),
                           zlib=True, complevel=4)
    v.long_name = ("Eggs or larvae per 10 m2 of sea surface (oblique and vertical tows), mean over "
                   "tows at the station; fill = not caught (see tows)")
    v.units = "count per 10 m2"
    write_sparse(v, areal_cells, ["cruise_i", "sampling_i", "net_i", "stage_i", "taxon_i", "lat_i", "lon_i"],
                 "mean_density")

    manta_dims = ("cruise", "sampling", "life_stage", "taxon", "lat", "lon")
    v = out.createVariable("abundance_per_100m3", "f8", manta_dims, fill_value=DOUBLE_FILL,
                           chunksizes=chunks_for([out.dimensions[d].size for d in manta_dims]),
                           zlib=True, complevel=4)
    v.long_name = ("Eggs or larvae per 100 m3 filtered (Manta surface tows), mean over tows at the "
                   "station; fill = not caught (see manta_tows)")
    v.units = "count per 100 m3"
    write_sparse(v, manta_cells, ["cruise_i", "sampling_i", "stage_i", "taxon_i", "lat_i", "lon_i"], "mean_density")

    tow_dims = ("cruise", "sampling", "net_type", "lat", "lon")
    v = out.createVariable("tows", "i4", tow_dims, fill_value=INT_FILL,
                           chunksizes=chunks_for([out.dimensions[d].size for d in tow_dims]),
                           zlib=True, complevel=4)
    v.long_name = "Number of oblique/vertical tows at the station (fill = not sampled)"
    v.units = "1"
    write_sparse(v, areal_tows, ["cruise_i", "sampling_i", "net_i", "lat_i", "lon_i"], "n_tows")

    manta_tow_dims = ("cruise", "sampling", "lat", "lon")
    v = out.createVariable("manta_tows", "i4", manta_tow_dims, fill_value=INT_FILL,
                           chunksizes=chunks_for([out.dimensions[d].size for d in manta_tow_dims]),
                           zlib=True, complevel=4)
    v.long_name = "Number of Manta tows at the station (fill = not sampled)"
    v.units = "1"
    write_sparse(v, manta_tows, ["cruise_i", "sampling_i", "lat_i", "lon_i"], "n_tows")

    out.close()
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
