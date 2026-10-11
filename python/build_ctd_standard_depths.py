# build_ctd_standard_depths.py

"""
Temperature, salinity and density (sigma-theta) at the CalCOFI standard
depths, one row per CTD cast and depth, from CalCOFI.io's CTD cast dataset
(calcofi_ctd-cast, 1993 onward), for the hydrographic atlas. The output has
the same columns as build_bottle_standard_depths.py's, so the atlas treats
both sources alike (they are never mixed).

Source: the CalCOFI.io database's thinned CTD table
(obs, dataset_key = calcofi_ctd-cast): 1 m-binned downcast profiles thinned
to every 10 m plus the points needed to keep the profile's shape. Values at
the standard depths are taken where a bin sits on the depth (within 0.5 m)
and otherwise interpolated linearly between the bins either side, never
across a gap wider than MAX_GAP_M and never extrapolated. The 0 m value is
the shallowest bin within SURFACE_MAX_M of the surface (the CTD's top bin is
usually 1-3 m). Upcasts (14 of ~9,700 casts) are ignored.

Variables:
  temperature  temperature_ave (mean of the two sensors)
  salinity     salinity_ave_corr (bottle-corrected, mean of the two sensors),
               else salinity_1_corr, else salinity_1 (uncorrected, mostly
               preliminary recent cruises); salinity_source records which
  sigma_theta  recomputed (EOS-80, eos80.py) from the temperature and the
               salinity chosen, as in the bottle database's STheta; the CTD
               table's own sigma_theta_1 uses sensor 1's uncorrected salinity
Quality flags 8 (questionable) and 9 (missing) are excluded; the same wide
range check as the bottle data follows.

Domain, seasons and station keys follow build_bottle_standard_depths.py.

Output: data/CalCOFI_CTD/calcofi_ctd_standard_depths.parquet

Usage:  python build_ctd_standard_depths.py [--release latest] [--out FILE]
"""

import argparse
import os

import duckdb
import numpy as np
import pandas as pd

import build_gridded_abundance_by_cruise_and_taxon as build
from build_bottle_standard_depths import S_RANGE, SEASONS, STANDARD_DEPTHS, T_RANGE
from eos80 import pressure_from_depth, sigma_theta

OUT = os.path.join(build.REPO_DIR, "data", "CalCOFI_CTD", "calcofi_ctd_standard_depths.parquet")
DB = "https://storage.googleapis.com/calcofi-db/ducklake/releases/{release}/parquet"
SALINITY_PREFERENCE = ["salinity_ave_corr", "salinity_1_corr", "salinity_1"]
TYPES = ["temperature_ave"] + SALINITY_PREFERENCE
BAD_FLAGS = ("8", "9")
SURFACE_MAX_M = 3.0
MAX_GAP_M = 20.0
TARGETS = np.array(STANDARD_DEPTHS, dtype=float)


def at_standard_depths(z, v):
    """Profile (depth z ascending, value v) -> values at TARGETS (NaN where not supported)."""
    out = np.full(len(TARGETS), np.nan)
    if not len(z):
        return out
    if z[0] <= SURFACE_MAX_M:
        out[0] = v[0]
    for k, target in enumerate(TARGETS[1:], start=1):
        i = np.searchsorted(z, target)
        if i < len(z) and abs(z[i] - target) <= 0.5:
            out[k] = v[i]
        elif i > 0 and abs(z[i - 1] - target) <= 0.5:
            out[k] = v[i - 1]
        elif 0 < i < len(z) and z[i] - z[i - 1] <= MAX_GAP_M:
            out[k] = v[i - 1] + (v[i] - v[i - 1]) * (target - z[i - 1]) / (z[i] - z[i - 1])
    return out


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
    types = ", ".join(f"'{t}'" for t in TYPES)
    obs = con.sql(f"""
        SELECT sample_key, cruise_key, grid_key, latitude AS lat, longitude AS lon, datetime,
               depth_min_m AS z, measurement_type AS type, measurement_value AS v
        FROM '{src}'
        WHERE measurement_type IN ({types}) AND sample_key LIKE '%d'
          AND measurement_value IS NOT NULL
          AND (measurement_qual IS NULL OR measurement_qual NOT IN {BAD_FLAGS})
    """).df()
    print(f"{len(obs)} downcast values, {obs.sample_key.nunique()} casts")

    casts = obs.groupby("sample_key").agg(cruise_key=("cruise_key", "first"), grid_key=("grid_key", "first"),
                                          lat=("lat", "first"), lon=("lon", "first"), datetime=("datetime", "first"))
    profiles = {k: g for k, g in obs.groupby(["sample_key", "type"])}
    rows = []
    for key in casts.index:
        def profile(t):
            g = profiles.get((key, t))
            if g is None:
                return None
            g = g.groupby("z", as_index=False).v.mean().sort_values("z")
            return at_standard_depths(g.z.to_numpy(), g.v.to_numpy())
        temp = profile("temperature_ave")
        sal, source = None, None
        for t in SALINITY_PREFERENCE:
            sal = profile(t)
            if sal is not None and np.isfinite(sal).any():
                source = t
                break
        if temp is None:
            temp = np.full(len(TARGETS), np.nan)
        if source is None:
            sal = np.full(len(TARGETS), np.nan)
        rows.append(pd.DataFrame({"sample_key": key, "depth_m": STANDARD_DEPTHS, "temperature": temp,
                                  "salinity": sal, "salinity_source": source}))
    out = pd.concat(rows, ignore_index=True).merge(casts, left_on="sample_key", right_index=True)
    out = out.dropna(subset=["temperature", "salinity"], how="all")
    print(f"standard-depth rows: {len(out)}; salinity source by cast:",
          out.groupby("salinity_source").sample_key.nunique().to_dict())

    # ---- range check, density ----
    bad_t = ~out.temperature.between(*T_RANGE) & out.temperature.notna()
    bad_s = ~out.salinity.between(*S_RANGE) & out.salinity.notna()
    out.loc[bad_t, "temperature"] = np.nan
    out.loc[bad_s, "salinity"] = np.nan
    print(f"  range check: {int(bad_t.sum())} temperature, {int(bad_s.sum())} salinity values excluded")
    out["sigma_theta"] = sigma_theta(out.temperature, out.salinity, pressure_from_depth(out.depth_m, out.lat))
    out = out.dropna(subset=["temperature", "salinity", "sigma_theta"], how="all")

    # ---- domain, season, station keys ----
    in_domain = (out.lat.between(*build.DOMAIN_LAT) & out.lon.between(*build.DOMAIN_LON)
                 & ~build.points_in_polygon(out.lon.to_numpy(float), out.lat.to_numpy(float), build.GULF_OF_CALIFORNIA))
    print(f"casts in domain: {out[in_domain].sample_key.nunique()} of {out.sample_key.nunique()}")
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
            "lat", "lon", "depth_m", "temperature", "salinity", "sigma_theta", "salinity_source"]
    out = out[cols].sort_values(["date", "cast", "depth_m"])

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    con.register("out", out)
    con.sql(f"COPY out TO '{args.out}' (FORMAT parquet, COMPRESSION zstd)")
    print(f"wrote {args.out}: {len(out)} rows, {out.cast.nunique()} casts, {out.site_key.nunique()} stations, "
          f"{out.year.min()}-{out.year.max()} ({os.path.getsize(args.out) / 1e6:.1f} MB)")
    print(out.groupby("season").cast.nunique().rename("casts").to_dict())


if __name__ == "__main__":
    main()
