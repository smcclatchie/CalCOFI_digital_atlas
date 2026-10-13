# build_derived_variables.py

"""
Derived variables for the hydrographic atlas's Derived variables view, per
source (bottle, CTD -- never mixed), from the standard-depth tables made by
build_bottle_standard_depths.py and build_ctd_standard_depths.py, plus
CalCOFI's own CTD-derived values.

One value per cast (a depth, m):
  z_sigma25         depth of the sigma-theta = 25.0 surface: its first
                    crossing down the standard-depth profile, interpolated
                    linearly; 0 where the surface water is already denser
                    (the surface outcrops); missing if the profile never
                    reaches 25.0
  nitricline_depth  depth of the steepest nitrate gradient: the midpoint of
                    the standard-depth interval with the largest increase of
                    nitrate with depth (dNO3/dz > 0), consecutive standard
                    depths only
  dyn_height        the 0/500 dbar dynamic height (dyn m): CalCOFI's own
                    dynamic height at 500 m (bottle R_DYNHT; CTD dynamic_height,
                    both accumulated from the surface), not recomputed; casts
                    reaching 500 m only
  chl_max_depth     CTD only: CalCOFI's chl_max_depth (the depth of the 3 m
                    running-mean chlorophyll maximum), not recomputed
  mld               CTD only: CalCOFI's headline mixed layer depth,
                    mld_sigma_theta_002 (sigma-theta 0.02 kg/m3 above its
                    10 m value), not recomputed
  chl_integrated    CTD only: CalCOFI's chlorophyll a integrated over the top
                    200 m (or to the cast bottom), mg/m2, not recomputed
  secchi            bottle only: the Secchi disk depth (m) recorded with the
                    cast (CalCOFI.org Cast.csv, 1969-2021, daylight stations)
  pp_integrated     bottle only: 14C primary production integrated over the
                    sampled depths (mg C m-2 per half light-day), from the
                    mean of the replicate assimilations (MeanAs) at the
                    observed light-level bottles (RecInd 3, flag 8 out):
                    trapezoids between bottles, the shallowest value held up
                    to the surface, nothing below the deepest bottle (the
                    incubations reach ~1% light, so little production is
                    missed); casts with at least 3 such bottles

At each standard depth / interval:
  spiciness         TEOS-10 spiciness referenced to 0 dbar (gsw.spiciness0,
                    kg/m3) from temperature and salinity at the standard depth
  n2                Brunt-Vaisala frequency squared (gsw.Nsquared), 1e-4 s^-2,
                    between consecutive standard depths; stored at the upper
                    depth of each interval (0 = 0-10 m, ..., 400 = 400-500 m)

Outputs (data/CalCOFI_derived/):
  calcofi_{bottle,ctd}_derived_casts.parquet   one row per cast
  calcofi_{bottle,ctd}_derived_depths.parquet  one row per cast and depth

Usage:  python build_derived_variables.py [--release latest]
"""

import argparse
import os

import duckdb
import gsw
import numpy as np
import pandas as pd

import build_gridded_abundance_by_cruise_and_taxon as build
from build_bottle_standard_depths import OUT as BOTTLE_TABLE
from build_bottle_standard_depths import SRC_DIR as BOTTLE_DIR
from build_ctd_standard_depths import DB, OUT as CTD_TABLE

OUT_DIR = os.path.join(build.REPO_DIR, "data", "CalCOFI_derived")
SIGMA_SURFACE = 25.0
CALCOFI_CTD_DERIVED = {"chl_max_depth": "chl_max_depth", "mld": "mld_sigma_theta_002", "chl_integrated": "chl_integrated"}
CAST_META = ["cast", "date", "year", "month", "season", "season_year", "site_key", "line", "station", "lat", "lon"]


def isopycnal_depth(z, sig, level=SIGMA_SURFACE):
    """First depth where the sigma-theta profile reaches `level` (0 if the surface already does)."""
    ok = np.isfinite(sig)
    z, sig = z[ok], sig[ok]
    if not len(z):
        return np.nan
    if sig[0] >= level:
        return 0.0 if z[0] <= 10 else np.nan     # outcrops only if the profile starts near the surface
    for i in range(1, len(z)):
        if sig[i] >= level:
            return z[i - 1] + (level - sig[i - 1]) * (z[i] - z[i - 1]) / (sig[i] - sig[i - 1])
    return np.nan


def steepest_gradient_depth(z, no3, standard):
    """Midpoint of the consecutive-standard-depth interval with the largest positive dNO3/dz."""
    ok = np.isfinite(no3)
    z, no3 = z[ok], no3[ok]
    best, depth = 0.0, np.nan
    for i in range(1, len(z)):
        if standard.index(z[i]) != standard.index(z[i - 1]) + 1:
            continue                                  # a missing standard depth in between
        g = (no3[i] - no3[i - 1]) / (z[i] - z[i - 1])
        if g > best:
            best, depth = g, (z[i] + z[i - 1]) / 2
    return depth


def derive(df, standard):
    casts, depths = [], []
    for cast, g in df.sort_values("depth_m").groupby("cast", sort=False):
        z = g.depth_m.to_numpy(float)
        meta = g.iloc[0]
        row = {c: meta[c] for c in CAST_META}
        row["z_sigma25"] = isopycnal_depth(z, g.sigma_theta.to_numpy(float))
        row["nitricline_depth"] = steepest_gradient_depth(z, g.nitrate.to_numpy(float), standard)
        at500 = g.loc[g.depth_m == 500, "dynamic_height"]
        row["dyn_height"] = float(at500.iloc[0]) if len(at500) and pd.notna(at500.iloc[0]) else np.nan
        casts.append(row)

        t, sp = g.temperature.to_numpy(float), g.salinity.to_numpy(float)
        p = gsw.p_from_z(-z, meta.lat)
        sa = gsw.SA_from_SP(sp, p, meta.lon, meta.lat)
        ct = gsw.CT_from_t(sa, t, p)
        spice = gsw.spiciness0(sa, ct)
        n2 = np.full(len(z), np.nan)
        ok = np.isfinite(sa) & np.isfinite(ct)
        zi = z[ok]
        if ok.sum() >= 2:
            n2_mid, _ = gsw.Nsquared(sa[ok], ct[ok], p[ok], meta.lat)
            for k in range(len(zi) - 1):
                if standard.index(zi[k + 1]) == standard.index(zi[k]) + 1:      # consecutive standard depths
                    n2[np.where(z == zi[k])[0][0]] = n2_mid[k] * 1e4
        for k in range(len(z)):
            depths.append({"cast": cast, "depth_m": int(z[k]), "spiciness": spice[k], "n2": n2[k]})
    return pd.DataFrame(casts), pd.DataFrame(depths)


PP_MIN_BOTTLES = 3


def bottle_primary_production():
    """Integrated 14C primary production per bottle cast (see the module docstring), keyed by cast."""
    b = pd.read_csv(os.path.join(BOTTLE_DIR, "194903-202105_Bottle.csv"), encoding="latin-1", low_memory=False,
                    usecols=["Cst_Cnt", "Depthm", "RecInd", "MeanAs", "MeanAq"])
    b = b[(b.RecInd == 3) & b.MeanAs.notna() & (b.MeanAq != 8)].sort_values(["Cst_Cnt", "Depthm"])
    out = {}
    for cast, g in b.groupby("Cst_Cnt"):
        g = g.drop_duplicates("Depthm")
        if len(g) < PP_MIN_BOTTLES:
            continue
        z, v = g.Depthm.to_numpy(float), g.MeanAs.to_numpy(float).clip(min=0)
        out[cast] = float(z[0] * v[0] + np.trapezoid(v, z))
    return pd.Series(out, name="pp_integrated")


def calcofi_ctd_derived(release):
    """CalCOFI's own per-cast CTD values (sample_measurement), keyed by our cast id."""
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs")
    types = ", ".join(f"'{t}'" for t in CALCOFI_CTD_DERIVED.values())
    m = con.sql(f"""SELECT sample_key, measurement_type, measurement_value
                    FROM '{DB.format(release=release)}/sample_measurement.parquet'
                    WHERE dataset_key = 'calcofi_ctd-derived' AND measurement_type IN ({types})
                      AND (measurement_qual IS NULL OR measurement_qual NOT IN ('8', '9'))""").df()
    m["cast"] = m.sample_key.str.split(":").str[-1]
    wide = m.pivot_table(index="cast", columns="measurement_type", values="measurement_value", aggfunc="first")
    return wide.rename(columns={v: k for k, v in CALCOFI_CTD_DERIVED.items()}).reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--release", default="latest")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    con = duckdb.connect()
    standard = sorted(con.sql(f"SELECT DISTINCT depth_m FROM '{BOTTLE_TABLE}'").df().depth_m.astype(int).tolist())
    for source, table in [("bottle", BOTTLE_TABLE), ("ctd", CTD_TABLE)]:
        df = con.sql(f"SELECT * FROM '{table}' WHERE site_key IS NOT NULL").df()
        casts, depths = derive(df, standard)
        if source == "bottle":
            sec = pd.read_csv(os.path.join(BOTTLE_DIR, "194903-202105_Cast.csv"), encoding="latin-1", low_memory=False,
                              usecols=["Cst_Cnt", "Secchi"]).dropna().rename(columns={"Cst_Cnt": "cast", "Secchi": "secchi"})
            casts = casts.merge(sec, on="cast", how="left")
            print(f"Secchi depth: {int(casts.secchi.notna().sum())} casts, median {casts.secchi.median():.0f} m")
            pp = bottle_primary_production()
            casts = casts.merge(pp, left_on="cast", right_index=True, how="left")
            print(f"14C integrated primary production: {int(casts.pp_integrated.notna().sum())} casts, "
                  f"median {casts.pp_integrated.median():.0f} mg C m-2 per half light-day")
        if source == "ctd":
            release = build.resolve_release(args.release)
            cd = calcofi_ctd_derived(release)
            casts = casts.merge(cd, on="cast", how="left")
            print(f"CalCOFI CTD-derived values (release {release}): "
                  f"chl_max_depth {int(casts.chl_max_depth.notna().sum())}, mld {int(casts.mld.notna().sum())}, "
                  f"chl_integrated {int(casts.chl_integrated.notna().sum())} casts")
        for name, out in [("casts", casts), ("depths", depths)]:
            path = os.path.join(OUT_DIR, f"calcofi_{source}_derived_{name}.parquet")
            con.register("out", out)
            con.sql(f"COPY out TO '{path}' (FORMAT parquet, COMPRESSION zstd)")
            con.unregister("out")
        print(f"{source}: {len(casts)} casts; z_sigma25 {int(casts.z_sigma25.notna().sum())} "
              f"(outcropping {int((casts.z_sigma25 == 0).sum())}), nitricline {int(casts.nitricline_depth.notna().sum())}, "
              f"0/500 dynamic height {int(casts.dyn_height.notna().sum())} (median {casts.dyn_height.median():.3f} dyn m); "
              f"spiciness {int(depths.spiciness.notna().sum())}, N2 {int(depths.n2.notna().sum())} values")


if __name__ == "__main__":
    main()
