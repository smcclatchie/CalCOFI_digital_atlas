# build_isopycnal_properties.py

"""
Temperature, salinity and dissolved oxygen on the sigma-theta = 25.8 and 26.4
surfaces, one row per CTD cast and surface, for the hydrographic atlas's Water
masses view, plus the depth of each surface.

Source: CalCOFI.io's thinned CTD downcasts (calcofi_ctd-cast; 1 m bins
thinned to every 10 m plus the points that keep the profile's shape), not the
standard depths: these surfaces lie at ~100-400 m, where the standard depths
are 50-100 m apart. Variables as in build_ctd_standard_depths.py: two-sensor
mean temperature; bottle-corrected salinity (else uncorrected, 2025-26
preliminary casts); sensor oxygen, station-corrected first. Sigma-theta is
computed (EOS-80, eos80.py) on each cast's own profile from that temperature
and salinity.

Each surface is taken at its first crossing going down the profile; the
depth and the three properties are interpolated linearly between the two
profile points bracketing it (never across a gap wider than MAX_GAP_M). A
cast that doesn't reach the surface, or starts below it, has no value.

Output: data/CalCOFI_isopycnal/calcofi_ctd_isopycnals.parquet

Usage:  python build_isopycnal_properties.py [--release latest]
"""

import argparse
import os

import duckdb
import numpy as np
import pandas as pd

import build_gridded_abundance_by_cruise_and_taxon as build
from build_bottle_standard_depths import O2_RANGE, S_RANGE, SEASONS, T_RANGE
from build_ctd_standard_depths import BAD_FLAGS, DB, OXYGEN_PREFERENCE, SALINITY_PREFERENCE
from eos80 import pressure_from_depth, sigma_theta

OUT = os.path.join(build.REPO_DIR, "data", "CalCOFI_isopycnal", "calcofi_ctd_isopycnals.parquet")
SURFACES = [25.8, 26.4]
MAX_GAP_M = 30.0


def first_crossing(z, sig, level):
    """Index i such that the profile first reaches `level` between z[i-1] and z[i] (or None)."""
    if not len(z) or sig[0] >= level:
        return None                                # starts denser: the surface is above the cast
    hit = np.nonzero(sig >= level)[0]
    if not len(hit):
        return None
    i = hit[0]
    return i if z[i] - z[i - 1] <= MAX_GAP_M else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--release", default="latest")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    release = build.resolve_release(args.release)
    src = f"{DB.format(release=release)}/obs/dataset_key=calcofi_ctd-cast/data_0.parquet"
    print(f"CalCOFI.io release {release}: {src}")
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs")
    types = ", ".join(f"'{t}'" for t in ["temperature_ave"] + SALINITY_PREFERENCE + OXYGEN_PREFERENCE)
    obs = con.sql(f"""
        SELECT sample_key, cruise_key, grid_key, latitude AS lat, longitude AS lon, datetime,
               depth_min_m AS z, measurement_type AS type, measurement_value AS v
        FROM '{src}'
        WHERE measurement_type IN ({types}) AND sample_key LIKE '%d'
          AND measurement_value IS NOT NULL
          AND (measurement_qual IS NULL OR measurement_qual NOT IN {BAD_FLAGS})
    """).df()
    print(f"{len(obs)} downcast values, {obs.sample_key.nunique()} casts")

    # One profile per cast: depth x (temperature, salinity, oxygen), each the preferred variant.
    wide = obs.pivot_table(index=["sample_key", "z"], columns="type", values="v", aggfunc="mean").reset_index()

    def preferred(df, names):
        out = pd.Series(np.nan, index=df.index)
        source = pd.Series(None, index=df.index, dtype=object)
        for n in names:                                            # first variant the cast has at all
            if n not in df:
                continue
            has = df.groupby("sample_key")[n].transform(lambda s: s.notna().any())
            take = has & source.isna()
            out[take] = df.loc[take, n]
            source[take] = n
        return out, source

    wide["temperature"] = wide["temperature_ave"]
    wide["salinity"], wide["salinity_source"] = preferred(wide, SALINITY_PREFERENCE)
    wide["oxygen"], wide["oxygen_source"] = preferred(wide, OXYGEN_PREFERENCE)
    for name, rng in [("temperature", T_RANGE), ("salinity", S_RANGE), ("oxygen", O2_RANGE)]:
        wide.loc[~wide[name].between(*rng), name] = np.nan

    meta = obs.groupby("sample_key").agg(cruise_key=("cruise_key", "first"), grid_key=("grid_key", "first"),
                                         lat=("lat", "first"), lon=("lon", "first"), datetime=("datetime", "first"))
    wide = wide.merge(meta[["lat"]], left_on="sample_key", right_index=True)
    wide["sigma_theta"] = sigma_theta(wide.temperature, wide.salinity, pressure_from_depth(wide.z, wide.lat))

    rows = []
    for key, g in wide.dropna(subset=["sigma_theta"]).sort_values("z").groupby("sample_key", sort=False):
        z, sig = g.z.to_numpy(float), g.sigma_theta.to_numpy(float)
        for level in SURFACES:
            i = first_crossing(z, sig, level)
            if i is None:
                continue
            f = (level - sig[i - 1]) / (sig[i] - sig[i - 1])
            interp = lambda col: g[col].iloc[i - 1] + f * (g[col].iloc[i] - g[col].iloc[i - 1])
            rows.append({"sample_key": key, "sigma_theta": level, "depth_m": z[i - 1] + f * (z[i] - z[i - 1]),
                         "temperature": interp("temperature"), "salinity": interp("salinity"), "oxygen": interp("oxygen"),
                         "salinity_source": g.salinity_source.iloc[i], "oxygen_source": g.oxygen_source.iloc[i]})
    out = pd.DataFrame(rows).merge(meta, left_on="sample_key", right_index=True)

    in_domain = (out.lat.between(*build.DOMAIN_LAT) & out.lon.between(*build.DOMAIN_LON)
                 & ~build.points_in_polygon(out.lon.to_numpy(float), out.lat.to_numpy(float), build.GULF_OF_CALIFORNIA))
    out = out[in_domain].copy()
    dt = pd.to_datetime(out.datetime)
    out["date"] = dt.dt.strftime("%Y-%m-%d")
    out["year"], out["month"] = dt.dt.year, dt.dt.month
    out["season"] = out.month.map(SEASONS)
    out["season_year"] = out.year + (out.month == 12)
    gk = out.grid_key.str.extract(r"^st(?P<station>[-\d.]+)-ln(?P<line>[-\d.]+)$").astype(float)
    out["line"], out["station"] = gk.line, gk.station
    out["site_key"] = [f"{l:05.1f} {s:05.1f}" if np.isfinite(l) and np.isfinite(s) else None
                       for l, s in zip(out.line, out.station)]
    out["cast"] = out.sample_key.str.split(":").str[-1]
    cols = ["cast", "cruise_key", "date", "year", "month", "season", "season_year", "site_key", "line", "station",
            "lat", "lon", "sigma_theta", "depth_m", "temperature", "salinity", "oxygen", "salinity_source", "oxygen_source"]
    out = out[cols].sort_values(["date", "cast", "sigma_theta"])

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    con.register("out", out)
    con.sql(f"COPY out TO '{args.out}' (FORMAT parquet, COMPRESSION zstd)")
    print(f"wrote {args.out}: {len(out)} rows, {out.cast.nunique()} casts, {out.year.min()}-{out.year.max()}")
    for level, g in out.groupby("sigma_theta"):
        print(f"  sigma-theta {level}: {len(g)} casts; depth median {g.depth_m.median():.0f} m "
              f"(1-99%: {g.depth_m.quantile(.01):.0f}-{g.depth_m.quantile(.99):.0f}); "
              f"T {g.temperature.median():.2f} C, S {g.salinity.median():.3f}, O2 {g.oxygen.median():.2f} ml/L")


if __name__ == "__main__":
    main()
