# build_bottle_standard_depths.py

"""
Temperature, salinity, density (sigma-theta), dissolved oxygen (ml/L),
nitrate, silicate and phosphate (umol/L) at the CalCOFI standard depths, one row per
cast and depth, from the CalCOFI bottle database, for the hydrographic atlas.

Source: the CalCOFI.org database CSVs (194903-202105_Bottle.csv and _Cast.csv,
https://calcofi.org/downloads/database/), NOT CalCOFI.io's calcofi_bottle.nc:
the NetCDF drops the record indicator (RecInd) and the quality flags, and both
are needed here.

Standard depths: 0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400,
500 m. At each, a cast has one record: RecInd 7 (interpolated onto the
standard depth), RecInd 3 (a bottle observed exactly there) or RecInd 5.
RecInd 6 only ever appears as a second copy beside a 3 and is dropped.

Quality flags (T_qual, S_qual, SThtaq; O_qual, NO3q, SiO3qu follow the same
pattern), as read from the data's own pattern
(flag x RecInd): 9 = missing (never has a value); 8 = questionable (only on
observed bottles) -> excluded; 6 = interpolated/estimated (almost all on the
RecInd 7 standard-depth records) -> kept; blank -> kept; the rare other codes
(~0.1%, including hundreds of stray S_qual values) -> excluded. The flags miss a
few gross errors (an unflagged T = 99 fill value interpolated into its
neighbours, salinity overwritten by sigma-theta), so a wide range check
follows; density is dropped wherever its temperature or salinity is.

Domain: the same CalCOFI domain as the ichthyoplankton products
(build_gridded_abundance_by_cruise_and_taxon.DOMAIN_LAT / DOMAIN_LON, outside
the Gulf of California), by the cast's recorded position.

Season: winter = Dec-Feb, spring = Mar-May, summer = Jun-Aug, autumn =
Sep-Nov. December belongs to the FOLLOWING year's winter (season_year).

Output: data/CalCOFI_bottle_data/calcofi_bottle_standard_depths.parquet

Usage:  python build_bottle_standard_depths.py [--src-dir DIR] [--out FILE]
"""

import argparse
import os

import duckdb
import pandas as pd

import build_gridded_abundance_by_cruise_and_taxon as build

SRC_DIR = os.path.join(build.REPO_DIR, "data", "CalCOFI_bottle_data")
OUT = os.path.join(SRC_DIR, "calcofi_bottle_standard_depths.parquet")
STANDARD_DEPTHS = [0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500]
KEEP_FLAGS = {6}                     # plus blank (NaN)
T_RANGE = (-2, 35)                   # degC; unflagged fill values (99) and bad interpolations sit outside
S_RANGE = (20, 35.5)                 # wide enough for real nearshore plume values (~28)
# Oxygen (ml/L), nitrate and silicate (umol/L): the same flag rules; ranges
# catch unflagged gross errors (e.g. an oxygen value of -200) and allow small
# negative nutrient readings at the detection limit.
O2_RANGE = (0, 12)
NO3_RANGE = (-0.5, 60)
SIO3_RANGE = (-0.5, 200)
PO4_RANGE = (-0.1, 6)
DYNHT_RANGE = (0, 3)                 # dynamic metres
# CalCOFI's reported dynamic height (R_DYNHT, dyn m): at each depth, the
# dynamic height from the surface down to it, so at 500 m it is the 0/500 dbar
# dynamic height. It has no quality flag.
SEASONS = {12: "winter", 1: "winter", 2: "winter", 3: "spring", 4: "spring", 5: "spring",
           6: "summer", 7: "summer", 8: "summer", 9: "autumn", 10: "autumn", 11: "autumn"}
VARIABLES = {"temperature": ("T_degC", "T_qual"), "salinity": ("Salnty", "S_qual"),
             "sigma_theta": ("STheta", "SThtaq"), "oxygen": ("O2ml_L", "O_qual"),
             "nitrate": ("NO3uM", "NO3q"), "silicate": ("SiO3uM", "SiO3qu"),
             "phosphate": ("PO4uM", "PO4q"), "dynamic_height": ("R_DYNHT", None)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--src-dir", default=SRC_DIR, help="folder holding the CalCOFI.org bottle/cast CSVs")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    bottle = pd.read_csv(os.path.join(args.src_dir, "194903-202105_Bottle.csv"), encoding="latin-1", low_memory=False,
                         usecols=["Cst_Cnt", "Depthm", "RecInd"] + [c for pair in VARIABLES.values() for c in pair if c])
    cast = pd.read_csv(os.path.join(args.src_dir, "194903-202105_Cast.csv"), encoding="latin-1", low_memory=False,
                       usecols=["Cst_Cnt", "Cruise_ID", "Cruise", "Sta_ID", "Date", "Year", "Month",
                                "Lat_Dec", "Lon_Dec", "Ship_Name", "Ship_Code"])
    print(f"{len(bottle)} bottle records, {len(cast)} casts")

    # ---- one record per cast and standard depth ----
    b = bottle[bottle.Depthm.isin(STANDARD_DEPTHS) & bottle.RecInd.isin([3, 5, 7])].copy()
    dup = b.duplicated(["Cst_Cnt", "Depthm"], keep=False)
    print(f"standard-depth records (RecInd 3/5/7): {len(b)}; duplicate cast/depth pairs: {int(dup.sum())}")
    b = b.drop_duplicates(["Cst_Cnt", "Depthm"], keep="first")

    # ---- quality flags ----
    out = b[["Cst_Cnt", "Depthm", "RecInd"]].rename(columns={"Depthm": "depth_m", "RecInd": "rec_ind"}).copy()
    for name, (value, flag) in VARIABLES.items():
        ok = pd.Series(True, index=b.index) if flag is None else (b[flag].isna() | b[flag].isin(KEEP_FLAGS))
        out[name] = b[value].where(ok)
        print(f"  {name}: {int(b[value].notna().sum())} values, {int((b[value].notna() & ~ok).sum())} excluded by flag")

    # ---- gross range check: errors the flags miss ----
    bad_t = ~out.temperature.between(*T_RANGE) & out.temperature.notna()
    bad_s = (~out.salinity.between(*S_RANGE) & out.salinity.notna()) | (out.salinity == out.sigma_theta)
    out.loc[bad_t, "temperature"] = None
    out.loc[bad_s, "salinity"] = None
    out.loc[bad_t | bad_s, "sigma_theta"] = None
    print(f"  range check: {int(bad_t.sum())} temperature, {int(bad_s.sum())} salinity values excluded (before the domain filter)")
    for name, rng in [("oxygen", O2_RANGE), ("nitrate", NO3_RANGE), ("silicate", SIO3_RANGE), ("phosphate", PO4_RANGE),
                      ("dynamic_height", DYNHT_RANGE)]:
        bad = ~out[name].between(*rng) & out[name].notna()
        out.loc[bad, name] = None
        print(f"  range check: {int(bad.sum())} {name} values excluded")
    out = out.dropna(subset=list(VARIABLES), how="all")

    # ---- cast metadata, domain, season ----
    out = out.merge(cast, on="Cst_Cnt", how="left")
    lon, lat = out.Lon_Dec.to_numpy(float), out.Lat_Dec.to_numpy(float)
    in_domain = (out.Lat_Dec.between(*build.DOMAIN_LAT) & out.Lon_Dec.between(*build.DOMAIN_LON)
                 & ~build.points_in_polygon(lon, lat, build.GULF_OF_CALIFORNIA))
    print(f"casts in domain: {out[in_domain].Cst_Cnt.nunique()} of {out.Cst_Cnt.nunique()}")
    out = out[in_domain].copy()
    out["season"] = out.Month.map(SEASONS)
    out["season_year"] = out.Year + (out.Month == 12)
    out["line"] = out.Sta_ID.str[:5].astype(float)
    out["station"] = out.Sta_ID.str[6:].astype(float)
    out = out.rename(columns={"Cst_Cnt": "cast", "Cruise_ID": "cruise_id", "Cruise": "cruise_yyyymm",
                              "Sta_ID": "site_key", "Date": "date", "Year": "year", "Month": "month",
                              "Lat_Dec": "lat", "Lon_Dec": "lon", "Ship_Name": "ship_name", "Ship_Code": "ship_code"})
    out["date"] = pd.to_datetime(out.date, format="%m/%d/%Y").dt.strftime("%Y-%m-%d")
    cols = ["cast", "cruise_id", "cruise_yyyymm", "ship_name", "ship_code", "date", "year", "month", "season",
            "season_year", "site_key", "line", "station", "lat", "lon", "depth_m", "rec_ind"] + list(VARIABLES)
    out = out[cols].sort_values(["date", "cast", "depth_m"])

    con = duckdb.connect()
    con.register("out", out)
    con.sql(f"COPY out TO '{args.out}' (FORMAT parquet, COMPRESSION zstd)")
    print(f"wrote {args.out}: {len(out)} rows, {out.cast.nunique()} casts, {out.site_key.nunique()} stations, "
          f"{out.year.min()}-{out.year.max()} ({os.path.getsize(args.out) / 1e6:.1f} MB)")
    print(out.groupby("season").cast.nunique().rename("casts").to_dict())


if __name__ == "__main__":
    main()
