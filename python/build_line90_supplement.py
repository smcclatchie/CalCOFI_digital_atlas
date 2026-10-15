# build_line90_supplement.py

"""
Line 90 spring larval counts newer than the CalCOFI.io release, as station
tables in the same layout as build_station_geoparquet.py's, for the
ichthyoplankton atlas.

Source: data/CalCOFI_ichthyoplankton/larval_line 90 through 2025_from_
Andrew.csv -- one row per line 90 station occupation on the spring cruises
1951-2025 (S_C cruise YYMM, S_SC ship code, S_L line, S_S station), one
column per taxon, larvae per 10 m2 (the standard oblique tow: CalCOFI ring
net to 1977, bongo since). Where it overlaps the CalCOFI.io tables, ~90% of
station values agree exactly; the rest are months with more than one tow or
net at a station.

Only the survey months the CalCOFI.io tables don't have on line 90 are kept
(2023-04, 2024-04 and 2025-04 as of the 2026 release), so a newer release
supersedes them automatically.

Taxa: a column is kept only where it is one atlas taxon -- its name (with
SPELLING fixes) matches a CalCOFI.io taxon, and on the overlap it agrees
with that taxon's values (>= MIN_AGREE of station values within 5%) at least
as well as with the sum of every record under that name. That second test
drops the lumped columns: several "<genus or family>1" columns (e.g.
Vinciguerria1, Gonostomatidae1) are every larva of that group without a
column of its own, not the genus- or family-level record, so they would be
wrong as either. Others (e.g. Sebastes1 = the Sebastes spp. record) pass.

Exception (user, 2026-10-16): Vinciguerria1 stands in for Vinciguerria
lucetia (SUBSTITUTE), so the atlas's warm-water group, which needs
V. lucetia, also reaches 2025 -- almost all Vinciguerria larvae are
V. lucetia, and Vinciguerria1 matches the summed Vinciguerria records on
97% of the overlap. Superseded once CalCOFI.io holds these months.

Outputs (data/CalCOFI_ichthyoplankton/line90_supplement/):
  line90_supplement_catch.parquet  as calcofi_ichthyo_stations_catch.parquet
                                   (no geometry): non-zero catches only
  line90_supplement_tows.parquet   as calcofi_ichthyo_stations_tows.parquet
  line90_supplement_taxa.csv       the taxa reported (taxon_key, column):
                                   for any other taxon these months are
                                   unsampled, not zero

Usage:  python build_line90_supplement.py
"""

import os

import duckdb
import numpy as np
import pandas as pd

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ICHTHYO_DIR = os.path.join(REPO_DIR, "data", "CalCOFI_ichthyoplankton")
SRC = os.path.join(ICHTHYO_DIR, "larval_line 90 through 2025_from_Andrew.csv")
DB_DIR = os.path.join(ICHTHYO_DIR, "ichthyoplankton_from_calCOFI.io")
OUT_DIR = os.path.join(ICHTHYO_DIR, "line90_supplement")
SPELLING = {"Lampanyctus ritterii": "Lampanyctus ritteri"}
SUBSTITUTE = {"Vinciguerria1": "Vinciguerria lucetia"}   # column -> taxon it stands in for (see above)
SHIPS = {"SH": ("BELL M. SHIMADA", "3322"), "RL": ("REUBEN LASKER", "33UD")}
MIN_AGREE = 0.85
ID_COLUMNS = ["S_C", "S_SC", "S_L", "S_S", "year", "season"]


def survey_month(code):
    yy, mm = divmod(int(code), 100)
    return f"{(1900 if yy > 40 else 2000) + yy}-{mm:02d}"


def agreement(a, b):
    """Share of station values within 5% (or 0.05) of each other."""
    j = a.to_frame("a").join(b.rename("b"), how="left").fillna({"b": 0}).dropna()
    return float(((j.a - j.b).abs() <= 0.05 * j[["a", "b"]].max(axis=1) + 0.05).mean())


def main():
    src = pd.read_csv(SRC, low_memory=False).drop(columns="Unnamed: 0")
    src["period"] = src.S_C.map(survey_month)
    con = duckdb.connect()
    db = con.sql(f"""SELECT strftime(date, '%Y-%m') AS d, substr(cruise_key, 1, 7) AS period, station, taxon_key,
                            scientific_name, family, order_taxon, abundance
                     FROM '{DB_DIR}/calcofi_ichthyo_stations_catch.parquet'
                     WHERE line = 90 AND life_stage = 'larva' AND sampling = 'standard' AND net_type IN ('C1', 'CB')""").df()
    db_months = set(con.sql(f"SELECT DISTINCT substr(cruise_key, 1, 7) FROM '{DB_DIR}/calcofi_ichthyo_stations_tows.parquet' "
                            "WHERE line = 90").df().iloc[:, 0])
    taxa = con.sql(f"SELECT DISTINCT scientific_name, taxon_key, common_name, rank, family, order_taxon, class "
                   f"FROM '{DB_DIR}/calcofi_ichthyo_stations_catch.parquet' WHERE life_stage = 'larva'").df()
    by_name = taxa.drop_duplicates("scientific_name").set_index("scientific_name")

    # ---- which columns are single atlas taxa ----
    overlap = src[src.period.isin(db_months)].set_index(["period", "S_S"])
    overlap.index.names = ["period", "station"]
    kept, report = {}, []
    for col in [c for c in src.columns if c not in ID_COLUMNS + ["period"]]:
        if col in SUBSTITUTE:
            kept[col] = SUBSTITUTE[col]
            report.append((col, by_name.loc[SUBSTITUTE[col], "taxon_key"], np.nan, np.nan, f"kept as {SUBSTITUTE[col]}"))
            continue
        name = SPELLING.get(col.replace(".", " "), col.replace(".", " ")).rstrip("1")
        if name not in by_name.index:
            report.append((col, "-", np.nan, np.nan, "no CalCOFI.io taxon of that name"))
            continue
        exact = db[db.scientific_name == name].groupby(["period", "station"]).abundance.mean()
        first = db.scientific_name.str.split().str[0]
        under = db[(db.scientific_name == name) | (first == name) | (db.family == name) | (db.order_taxon == name)]
        lump = under.groupby(["period", "station", "scientific_name"]).abundance.mean().groupby(["period", "station"]).sum()
        a_exact, a_lump = agreement(overlap[col], exact), agreement(overlap[col], lump)
        ok = a_exact >= MIN_AGREE and a_exact >= a_lump
        report.append((col, by_name.loc[name, "taxon_key"], round(a_exact, 2), round(a_lump, 2),
                       "kept" if ok else "dropped: lumped group" if a_lump > a_exact else "dropped: disagrees"))
        if ok:
            kept[col] = name
    report = pd.DataFrame(report, columns=["column", "taxon_key", "agree_taxon", "agree_all_under_name", "decision"])
    print(report[~report.decision.str.startswith("kept")].to_string(index=False))
    print(report[report.decision.str.startswith("kept as")].to_string(index=False))
    print(f"{len(kept)} of {len(report)} columns kept")

    # ---- the new survey months ----
    new = src[~src.period.isin(db_months)].copy()
    print("new survey months:", new.groupby("period").size().to_dict())
    unknown = set(new.S_SC) - set(SHIPS)
    if unknown:
        raise SystemExit(f"ship codes without a name: {unknown}")
    stations = con.sql(f"SELECT DISTINCT station, lat, lon FROM '{DB_DIR}/calcofi_ichthyo_stations_tows.parquet' "
                       "WHERE line = 90").df().groupby("station")[["lat", "lon"]].median()
    new["ship_name"] = new.S_SC.map(lambda s: SHIPS[s][0])
    new["cruise_key"] = new.period + "-" + new.S_SC.map(lambda s: SHIPS[s][1])
    new["line"], new["station"] = 90.0, new.S_S.astype(float)
    missing = set(new.station) - set(stations.index)
    if missing:
        raise SystemExit(f"line 90 stations with no known position: {missing}")
    new = new.join(stations, on="station")
    new["site_key"] = new.apply(lambda r: f"{r.line:05.1f} {r.station:05.1f}", axis=1)
    new["year"] = new.period.str[:4].astype(int)
    new["month"] = new.period.str[5:].astype(int)
    tows = pd.DataFrame({
        "cruise_key": new.cruise_key, "date": pd.NaT, "date_iso": None, "year": new.year, "month": new.month,
        "ship_name": new.ship_name, "sampling": "standard", "net_type": "CB",
        "net_description": "CB: CalCOFI bongo (oblique)", "site_key": new.site_key, "line": new.line,
        "station": new.station, "lat": new.lat, "lon": new.lon, "n_tows": 1})

    long = new.melt(id_vars=["cruise_key", "site_key"], value_vars=list(kept), var_name="column", value_name="abundance")
    long = long[long.abundance > 0].merge(tows.drop(columns="n_tows"), on=["cruise_key", "site_key"])
    long["scientific_name"] = long.column.map(kept)
    long = long.join(by_name[["taxon_key", "common_name", "rank", "family", "order_taxon", "class"]], on="scientific_name")
    long["life_stage"], long["units"], long["n_tows"] = "larva", "count per 10 m2", 1
    long["log10_abundance"] = np.log10(long.abundance)
    catch_cols = ["cruise_key", "date", "date_iso", "year", "month", "ship_name", "sampling", "net_type", "life_stage",
                  "taxon_key", "scientific_name", "common_name", "rank", "family", "order_taxon", "class", "abundance",
                  "log10_abundance", "units", "net_description", "site_key", "line", "station", "lat", "lon", "n_tows"]
    catch = long[catch_cols]

    os.makedirs(OUT_DIR, exist_ok=True)
    for df, name in [(catch, "catch"), (tows, "tows")]:
        con.register("t", df)
        con.sql(f"COPY t TO '{OUT_DIR}/line90_supplement_{name}.parquet' (FORMAT parquet, COMPRESSION zstd)")
        con.unregister("t")
    (pd.DataFrame({"taxon_key": [by_name.loc[n, "taxon_key"] for n in kept.values()], "column": list(kept)})
       .to_csv(f"{OUT_DIR}/line90_supplement_taxa.csv", index=False))
    report.to_csv(f"{OUT_DIR}/line90_supplement_column_check.csv", index=False)
    print(f"{len(tows)} tows, {len(catch)} non-zero catches, {len(kept)} taxa -> {OUT_DIR}")


if __name__ == "__main__":
    main()
