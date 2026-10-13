# build_carbonate_standard_depths.py

"""
Dissolved inorganic carbon (DIC) and total alkalinity (umol/kg) at the
CalCOFI standard depths, for the hydrographic atlas, from CalCOFI.io's
CalCOFI DIC dataset (calcofi_dic, 1983-2021: ~3,900 bottle samples at ~50
stations, sparse before 2008).

Each DIC profile (one cruise, station and day) is interpolated to the
standard depths as the CTD casts' bottle nutrients are (build_ctd_standard_
depths.at_standard_depths): a sample within 0.5 m of a standard depth is
taken as it is, otherwise linear between the samples either side, across
gaps up to MAX_GAP_M (the samples below 100 m are ~60-200 m apart), never
extrapolated; 0 m = a sample within SURFACE_MAX_M.

Profiles: the samples at one position (rounded to 0.01 deg) in one month.
Since ~2008 the dataset records the station's nominal position and a
placeholder date (the 20th of the cruise month), so each profile is attached
to the nearest bottle cast, and the nearest CTD cast, within MATCH_DEG and
MATCH_DAYS (nearest in space, then in time), keyed by that source's cast.

Outputs (data/CalCOFI_carbonate/): calcofi_{bottle,ctd}_carbonate.parquet,
one row per cast and standard depth (cast, depth_m, dic, alkalinity).

Usage:  python build_carbonate_standard_depths.py [--release latest]
"""

import argparse
import os

import duckdb
import numpy as np
import pandas as pd

import build_gridded_abundance_by_cruise_and_taxon as build
from build_bottle_standard_depths import OUT as BOTTLE_TABLE
from build_ctd_standard_depths import DB, OUT as CTD_TABLE, at_standard_depths

OUT_DIR = os.path.join(build.REPO_DIR, "data", "CalCOFI_carbonate")
MAX_GAP_M = 210.0
SURFACE_MAX_M = 5.0
MATCH_DEG = 0.08                     # ~9 km: nominal vs actual station positions
MATCH_DAYS = 20
RANGES = {"dic": (1800, 2500), "alkalinity": (2100, 2500)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--release", default="latest")
    args = ap.parse_args()
    release = build.resolve_release(args.release)
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs")
    src = f"{DB.format(release=release)}/obs/dataset_key=calcofi_dic/data_0.parquet"
    o = con.sql(f"""SELECT strftime(datetime, '%Y-%m') AS ym, strftime(datetime, '%Y-%m-%d') AS day,
                           round(latitude, 2) AS lat, round(longitude, 2) AS lon,
                           depth_min_m AS z, measurement_type AS t, measurement_value AS v
                    FROM '{src}' WHERE measurement_type IN ('dic', 'alkalinity') AND measurement_value IS NOT NULL""").df()
    for t, (lo, hi) in RANGES.items():
        bad = (o.t == t) & ~o.v.between(lo, hi)
        print(f"{t}: {int((o.t == t).sum())} values, {int(bad.sum())} outside {lo}-{hi} left out")
        o = o[~bad]
    print(f"CalCOFI.io release {release}: {o.groupby(['ym', 'lat', 'lon']).ngroups} DIC profiles")

    rows = []
    for (ym, lat, lon), g in o.groupby(["ym", "lat", "lon"]):
        rec = {"day": pd.Timestamp(g.day.min()), "lat": lat, "lon": lon}
        for t in ("dic", "alkalinity"):
            p = g[g.t == t].groupby("z", as_index=False).v.mean().sort_values("z")
            rec[t] = at_standard_depths(p.z.to_numpy(float), p.v.to_numpy(float), max_gap=MAX_GAP_M, surface_max=SURFACE_MAX_M)
        rows.append(rec)
    profiles = pd.DataFrame(rows)

    os.makedirs(OUT_DIR, exist_ok=True)
    from build_ctd_standard_depths import TARGETS
    for source, table in [("bottle", BOTTLE_TABLE), ("ctd", CTD_TABLE)]:
        casts = con.sql(f"""SELECT DISTINCT "cast", date, lat, lon FROM '{table}'""").df()
        casts["date"] = pd.to_datetime(casts.date)
        picks = []
        for r in profiles.itertuples():
            c = casts[(casts.date - r.day).abs() <= pd.Timedelta(days=MATCH_DAYS)]
            d = np.hypot(c.lat - r.lat, c.lon - r.lon)
            c = c.assign(d=d, dt=(c.date - r.day).abs())[d < MATCH_DEG]
            if len(c):
                picks.append((r.Index, c.sort_values(["d", "dt"]).iloc[0]["cast"]))
        m = profiles.loc[[i for i, _ in picks]].assign(cast=[k for _, k in picks]).drop_duplicates("cast")
        out = []
        for r in m.itertuples():
            out.append(pd.DataFrame({"cast": r.cast, "depth_m": TARGETS.astype(int), "dic": r.dic, "alkalinity": r.alkalinity}))
        out = pd.concat(out, ignore_index=True).dropna(subset=["dic", "alkalinity"], how="all")
        path = os.path.join(OUT_DIR, f"calcofi_{source}_carbonate.parquet")
        con.register("out", out)
        con.sql(f"COPY out TO '{path}' (FORMAT parquet, COMPRESSION zstd)")
        con.unregister("out")
        print(f"{source}: {len(m)} of {len(profiles)} profiles matched to a cast; {out['cast'].nunique()} casts, "
              f"{int(out.dic.notna().sum())} DIC and {int(out.alkalinity.notna().sum())} alkalinity values -> {path}")


if __name__ == "__main__":
    main()
