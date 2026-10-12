# build_water_mass_fractions.py

"""
Water-mass fractions -- Pacific Subarctic (PSA), Pacific Equatorial (PEW) and
North Pacific Central (NPCW) water -- at every cast and CalCOFI standard depth
from 0 to 400 m, for the hydrographic atlas's Water masses view.

Method: the extended Optimum Multiparameter (OMP) analysis of Bograd,
Schroeder & Jacox (2019, GRL 46:6690-6698), as implemented in the Fishing
Maps repository's calcofi_water_mass_omp.py (recreational_fishing_high_
resolution), there solved at three reference stations in winter; here the
same water types and solver are applied to every cast and standard depth:
  - 6 water types (upper/lower end-members of PSA, PEW, NPCW): T, S, O2, PO4
    and SiO4 from that analysis's Table 2 end-members; NO3 end-members from
    its regressions of bottle nitrate on temperature at the reference
    stations (CalCOFI.io v2026.10.08, values below)
  - 6 parameters (T, S, O2 [umol/L], PO4, NO3, SiO4) + mass conservation,
    each property row normalised across the water types (Tomczak & Large
    1989), the mass row weighted like temperature, plus a phosphate
    remineralisation term deltaP (Redfield O/P 170, N/P 16, Si/P 18);
    non-negative least squares; fractions renormalised to sum to 1
  - PSA, PEW, NPCW = the sum of each water mass's two end-members
Only rows with all six properties are solved. "upper ocean (0-400 db)", as
in that analysis.

Inputs: the standard-depth tables (build_bottle_standard_depths.py,
build_ctd_standard_depths.py; CTD nutrients are the bottles on the CTD casts).
Outputs: data/CalCOFI_water_masses/calcofi_{bottle,ctd}_water_mass_fractions.parquet

Usage:  python build_water_mass_fractions.py
"""

import os

import duckdb
import numpy as np
import pandas as pd
from scipy.optimize import nnls

import build_gridded_abundance_by_cruise_and_taxon as build
from build_bottle_standard_depths import OUT as BOTTLE_TABLE
from build_ctd_standard_depths import OUT as CTD_TABLE

OUT_DIR = os.path.join(build.REPO_DIR, "data", "CalCOFI_water_masses")
MAX_DEPTH_M = 400
ML_L_TO_UMOL_L = 44.661
R_O_P, R_N_P, R_SI_P = 170.0, 16.0, 18.0
# Table 2 end-members (O2 in mL/L), as in calcofi_water_mass_omp.py
TABLE2_ENDPOINTS = {
    "PSA":  {"T_lo": 8.87, "T_hi": 11.91, "S_lo": 33.59, "S_hi": 33.46, "O2_lo_mll": 3.58, "O2_hi_mll": 5.44,
             "PO4_lo": 1.72, "PO4_hi": 0.90, "SiO4_lo": 23.20, "SiO4_hi": 6.91},
    "PEW":  {"T_lo": 7.66, "T_hi": 10.26, "S_lo": 34.32, "S_hi": 34.07, "O2_lo_mll": 0.70, "O2_hi_mll": 2.39,
             "PO4_lo": 2.86, "PO4_hi": 2.00, "SiO4_lo": 55.10, "SiO4_hi": 25.98},
    "NPCW": {"T_lo": 6.25, "T_hi": 9.37, "S_lo": 34.13, "S_hi": 33.92, "O2_lo_mll": 1.21, "O2_hi_mll": 3.66,
             "PO4_lo": 2.81, "PO4_hi": 1.65, "SiO4_lo": 63.11, "SiO4_hi": 20.78},
}
# NO3 end-members (lower, upper), umol/L: calcofi_water_mass_omp.derive_no3_endpoints
# on CalCOFI.io v2026.10.08 (PSA n=174 r=-0.82; PEW n=20 r=-0.96; NPCW n=251 r=-0.91)
NO3_ENDPOINTS = {"PSA": (22.8565, 8.0305), "PEW": (34.7704, 24.4490), "NPCW": (37.6931, 22.1450)}
PARAMS = ["T", "S", "O2", "PO4", "NO3", "SiO4"]
WATER_MASSES = ["PSA", "PEW", "NPCW"]
WT_NAMES = [f"{wm}_{b}" for wm in WATER_MASSES for b in ("lower", "upper")]
COLUMNS = {"T": "temperature", "S": "salinity", "O2": "oxygen", "PO4": "phosphate", "NO3": "nitrate", "SiO4": "silicate"}


def water_types():
    wt = {}
    for wm in WATER_MASSES:
        e, (no3_lo, no3_hi) = TABLE2_ENDPOINTS[wm], NO3_ENDPOINTS[wm]
        wt[f"{wm}_lower"] = {"T": e["T_lo"], "S": e["S_lo"], "O2": e["O2_lo_mll"] * ML_L_TO_UMOL_L,
                             "PO4": e["PO4_lo"], "NO3": no3_lo, "SiO4": e["SiO4_lo"]}
        wt[f"{wm}_upper"] = {"T": e["T_hi"], "S": e["S_hi"], "O2": e["O2_hi_mll"] * ML_L_TO_UMOL_L,
                             "PO4": e["PO4_hi"], "NO3": no3_hi, "SiO4": e["SiO4_hi"]}
    return wt


def normalized_g(wt):
    raw = np.array([[wt[n][p] for n in WT_NAMES] for p in PARAMS])
    means, stds = raw.mean(axis=1, keepdims=True), raw.std(axis=1, keepdims=True)
    g = np.hstack([(raw - means) / stds, np.array([0.0, 0.0, -R_O_P, 1.0, R_N_P, R_SI_P]).reshape(-1, 1) / stds])
    t_std = stds[PARAMS.index("T"), 0]
    g = np.vstack([g, np.array([[1.0] * len(WT_NAMES) + [0.0]]) / t_std])
    return g, means.ravel(), stds.ravel(), t_std


def solve(df):
    g, means, stds, t_std = normalized_g(water_types())
    obs = np.column_stack([df[COLUMNS[p]].to_numpy(float) for p in PARAMS])
    obs[:, PARAMS.index("O2")] *= ML_L_TO_UMOL_L
    d = np.column_stack([(obs - means) / stds, np.full(len(obs), 1.0 / t_std)])
    out = np.full((len(obs), len(WATER_MASSES) + 2), np.nan)
    for i, row in enumerate(d):
        x, resid = nnls(g, row)
        frac, total = x[:-1], x[:-1].sum()
        if total > 0:
            frac = frac / total
            out[i, :3] = [frac[2 * k] + frac[2 * k + 1] for k in range(3)]
            out[i, 3], out[i, 4] = x[-1], resid
    return pd.DataFrame(out, columns=["PSA", "PEW", "NPCW", "delta_p", "residual"], index=df.index)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    con = duckdb.connect()
    for source, table in [("bottle", BOTTLE_TABLE), ("ctd", CTD_TABLE)]:
        cols = ", ".join(COLUMNS.values())
        df = con.sql(f"""SELECT "cast", date, year, month, season, season_year, site_key, line, station, lat, lon,
                                depth_m, {cols}
                         FROM '{table}' WHERE depth_m <= {MAX_DEPTH_M} AND site_key IS NOT NULL""").df()
        df = df.dropna(subset=list(COLUMNS.values())).reset_index(drop=True)
        out = pd.concat([df, solve(df)], axis=1)
        path = os.path.join(OUT_DIR, f"calcofi_{source}_water_mass_fractions.parquet")
        con.register("out", out)
        con.sql(f"COPY out TO '{path}' (FORMAT parquet, COMPRESSION zstd)")
        con.unregister("out")
        means = out[["PSA", "PEW", "NPCW"]].mean().round(3).to_dict()
        print(f"{source}: {len(out)} cast-depths, {out['cast'].nunique()} casts, {out.year.min()}-{out.year.max()}; "
              f"mean fractions {means}; median residual {out.residual.median():.3f}")


if __name__ == "__main__":
    main()
