# build_station_geoparquet.py

"""
Station-level GeoParquet tables of CalCOFI ichthyoplankton abundance, for
vector GIS and web tools (GeoLibre, QGIS, DuckDB, MapLibre/deck.gl):

  calcofi_ichthyo_stations_catch.parquet  one point per cruise x sampling
      type x net type x life stage x taxon x station, non-zero catches only,
      with standardised abundance and taxon names
  calcofi_ichthyo_stations_tows.parquet   one point per station occupation
      (cruise x sampling type x net type x station) with the tow count, so
      a map can show where a station was sampled and nothing was caught

Every selection and standardisation rule is
build_gridded_abundance_by_cruise_and_taxon.load_selected_occurrences(), so
these tables and the gridded NetCDF always agree. The differences are in
the geometry only: each point is the station's median recorded position
(no grid snapping), and repeat tows are averaged per station rather than
per grid cell (mean over tows, absences counting as zero).

Usage (defaults are this project's own paths and the current release):
    python build_station_geoparquet.py --src swfsc_ichthyo.nc --out-dir . --release latest

Dependencies: numpy, pandas, netCDF4, duckdb (spatial extension, installed
on first use, writes GeoParquet 1.0 metadata).
"""

import argparse
import os

import duckdb
import numpy as np
import pandas as pd

import build_gridded_abundance_by_cruise_and_taxon as build

CATCH_NAME = "calcofi_ichthyo_stations_catch.parquet"
TOWS_NAME = "calcofi_ichthyo_stations_tows.parquet"
UNITS = {True: "count per 100 m3", False: "count per 10 m2"}   # keyed by "is Manta"


def read_taxon_table(release):
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs;")
    url = build.PARQUET_URL.format(release=release, table="taxon")
    return con.sql(f"""SELECT taxon_key, common_name, rank, family, order_taxon, class
                       FROM read_parquet('{url}')""").df().set_index("taxon_key")


def write_geoparquet(df, path, sort_by):
    """Write df (with lon/lat columns) as GeoParquet, point geometry last."""
    con = duckdb.connect()
    con.sql("INSTALL spatial; LOAD spatial;")
    con.register("df", df)
    order = ", ".join(sort_by)
    con.sql(f"""COPY (SELECT *, ST_Point(lon, lat) AS geometry FROM df ORDER BY {order})
                TO '{path}' (FORMAT parquet, COMPRESSION zstd)""")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--src", default=build.DEFAULT_SRC, help="CalCOFI.io swfsc_ichthyo.nc")
    ap.add_argument("--out-dir", default=build.BASE_DIR, help="directory for the two GeoParquet files")
    ap.add_argument("--release", default=build.DEFAULT_RELEASE,
                    help="CalCOFI.io release for the cruise/sample/taxon tables ('latest' = current)")
    args = ap.parse_args()
    args.release = build.resolve_release(args.release)
    print(f"CalCOFI.io release {args.release}")

    kept, occ_df, per_cruise, src_release = build.load_selected_occurrences(args.src, args.release)

    # ---- stations: occupied line/station at the median recorded position ----
    stations = kept.groupby("site_key")[["lat", "lon"]].median()
    stations["line"] = stations.index.str[:5].astype(float)
    stations["station"] = stations.index.str[6:].astype(float)

    # ---- station occupations (the denominator for the mean) ----
    occ_keys = ["cruise_key", "sampling", "tow_type", "site_key"]
    tows = (kept.groupby(occ_keys)
                .agg(n_tows=("tow", "size"), time=("time", "min"))
                .reset_index())
    tows["date"] = pd.to_datetime(tows.time, unit="s", utc=True).dt.date
    # Plain ISO string too: GeoLibre's Time Slider binding detects ISO date
    # strings, epoch numbers and bare years, not necessarily Parquet DATE.
    tows["date_iso"] = tows.date.astype(str)
    tows["year"] = tows.cruise_key.str[:4].astype(int)
    tows["month"] = tows.cruise_key.str[5:7].astype(int)
    tows["ship_name"] = tows.cruise_key.map(per_cruise.ship_name)
    tows["net_description"] = tows.tow_type.map({**build.AREAL_NET_TYPES, build.MANTA: "MT: Manta (surface)"})
    tows = tows.drop(columns="time").join(stations, on="site_key").rename(columns={"tow_type": "net_type"})
    print(f"{len(stations)} stations, {len(tows)} station occupations")

    # ---- non-zero catches: mean standardised abundance over the occupation's tows ----
    catch = (occ_df.groupby(occ_keys + ["life_stage", "taxon_key"])
                   .agg(density_sum=("density", "sum"), scientific_name=("scientific_name", "first"))
                   .reset_index()
                   .rename(columns={"tow_type": "net_type"})
                   .merge(tows, on=["cruise_key", "sampling", "net_type", "site_key"]))
    catch["abundance"] = catch.density_sum / catch.n_tows
    catch["log10_abundance"] = np.log10(catch.abundance)   # skewed: for graduated styling
    catch["units"] = (catch.net_type == build.MANTA).map(UNITS)
    catch = catch.drop(columns="density_sum").join(read_taxon_table(args.release), on="taxon_key")
    print(f"{len(catch)} non-zero catch rows, {catch.taxon_key.nunique()} taxa "
          f"({catch.drop_duplicates('taxon_key').common_name.notna().sum()} with a common name)")

    column_order = ["cruise_key", "date", "date_iso", "year", "month", "ship_name", "sampling", "net_type",
                    "net_description", "site_key", "line", "station", "lat", "lon", "n_tows"]
    catch = catch[column_order[:8] + ["life_stage", "taxon_key", "scientific_name", "common_name",
                                      "rank", "family", "order_taxon", "class", "abundance",
                                      "log10_abundance", "units"]
                  + column_order[8:]]
    tows = tows[column_order]

    os.makedirs(args.out_dir, exist_ok=True)
    catch_path = os.path.join(args.out_dir, CATCH_NAME)
    tows_path = os.path.join(args.out_dir, TOWS_NAME)
    # Sorted by taxon so a browser query for one species touches few row groups.
    write_geoparquet(catch, catch_path, ["taxon_key", "cruise_key", "site_key"])
    write_geoparquet(tows, tows_path, ["cruise_key", "site_key"])
    for p in (catch_path, tows_path):
        print(f"wrote {p} ({os.path.getsize(p) / 1e6:.1f} MB)")
    print(f"source: swfsc_ichthyo.nc (CalCOFI.io release {src_release}); "
          f"cruise/sample/taxon tables release {args.release}")


if __name__ == "__main__":
    main()
