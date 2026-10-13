# build_zooplankton_tables.py

"""
Zooplankton tables for the atlas's Zooplankton view, from the CalCOFI.io
database, with the same CalCOFI selection as the ichthyoplankton products
(build_gridded_abundance_by_cruise_and_taxon.select_calcofi_nets and
classify_sampling: CalCOFI cruises, the CalCOFI domain, CalCOFI grid
stations, standard vs special-survey sampling). One row per tow:

  plankton_volume  swfsc_ichthyo nets with displacement volumes: total and
                   small plankton volume (mL) per 1,000 m3 filtered
                   (volume / volume_sampled * 1000), by net type (C1 ring
                   net to 1978, CB bongo after: kept apart, as for the
                   fish -- the net change is a known break in the series)
  krill            cce-lter_euphausiids (Brinton & Townsend), abundance per
                   m2 by species and life stage. Sub-stages (furcilia F1-F7,
                   calyptopis C1-C3) are summed into their stage: they are
                   recorded instead of the stage, never beside it (1 of
                   5,740 tow x species cases has both). "damaged", "larvae"
                   and "metanauplius" count only in the all-stages total.
  zooscan          cce-lter_zooscan (PRPOOS net, ZooScan imaging), per
                   image class: abundance (per m2), carbon biomass (mg C
                   m-2), mean Feret diameter (mm), mean individual carbon
                   (ug C)

Stations and cruises: the krill and ZooScan tows were taken on CalCOFI
cruises, often in the same hauls as the ichthyoplankton nets, but their own
keys are rough: the krill file rounds line numbers (63, 67, 73, 77, 83, 87,
93 for 63.3, 66.7, 73.3, 76.7, 83.3, 86.7, 93.3; 2,112 tows) and sometimes
station numbers, and ZooScan tows since 2016 have no cruise key (421). Each
tow therefore takes the station and cruise of the ichthyoplankton net on the
same date within MATCH_DEG, or failing a net the CTD cast (which also gives
the ZooScan tows without one their cruise), so all three views agree; unmatched tows keep their own keys, the rounded lines
mapped to the standard ones (ROUNDED_LINES).

CUFES (swfsc_cufes, Continuous Underway Fish Egg Sampler, 1996-2022): pumped
samples along the ship's track, eggs per m3 of anchovy, sardine, jack
mackerel, hake, market squid and other fish. Not at stations, so each sample
keeps its own position; selected by cruise (a CalCOFI cruise, as above) and
the CalCOFI domain, its sampling "standard" on a hydrographic CalCOFI cruise
and "special survey" otherwise (e.g. the sardine survey cruises).

Outputs (data/CalCOFI_zooplankton/):
  plankton_volume_tows.parquet   one row per net tow, both volumes
  krill_tows.parquet, krill_values.parquet      tows; non-zero values per
                                                tow x taxon x stage
  zooscan_tows.parquet, zooscan_values.parquet  tows; values per tow x
                                                class x measurement
  cufes_samples.parquet, cufes_values.parquet   samples; non-zero eggs per m3
                                                per sample x taxon

Usage:  python build_zooplankton_tables.py [--release latest]
"""

import argparse
import os

import duckdb
import numpy as np
import pandas as pd

import build_gridded_abundance_by_cruise_and_taxon as build

OUT_DIR = os.path.join(build.REPO_DIR, "data", "CalCOFI_zooplankton")
KRILL_STAGES = ["adult", "juvenile", "furcilia", "calyptopis"]       # sub-stages roll up; others: all-stages only
ZOOSCAN_MEASUREMENTS = ["zooscan_abundance", "zooscan_biomass_carbon", "zooscan_feret_diameter", "zooscan_carbon_individual"]
TOW_COLUMNS = ["sample_key", "cruise_key", "ship_name", "site_key", "line", "station", "lat", "lon", "date", "year", "month",
               "season", "season_year", "sampling"]
MATCH_DEG = 0.05                                                     # ~5 km
ROUNDED_LINES = {"063.0": "063.3", "067.0": "066.7", "073.0": "073.3", "077.0": "076.7",
                 "083.0": "083.3", "087.0": "086.7", "093.0": "093.3"}
SEASONS = {12: "winter", 1: "winter", 2: "winter", 3: "spring", 4: "spring", 5: "spring",
           6: "summer", 7: "summer", 8: "summer", 9: "autumn", 10: "autumn", 11: "autumn"}


def match_ichthyo(tows, nets, label):
    """Take each tow's site and cruise keys from the ichthyoplankton net (or
    CTD cast) on the same date within MATCH_DEG (see the module docstring)."""
    t = tows.copy()
    t["day"] = pd.to_datetime(t.datetime).dt.strftime("%Y-%m-%d")
    m = t[["sample_key", "day", "latitude", "longitude"]].merge(nets, on="day", suffixes=("", "_n"))
    m["d"] = np.hypot(m.latitude - m.latitude_n, m.longitude - m.longitude_n)
    m = m[m.d < MATCH_DEG]
    best = m.loc[m.groupby("sample_key").d.idxmin()].set_index("sample_key")
    hit = t.sample_key.isin(best.index)
    t.loc[hit, "site_key"] = t.loc[hit, "sample_key"].map(best.site_key_n)
    t.loc[hit, "cruise_key"] = t.loc[hit, "sample_key"].map(best.cruise_key_n)
    line = t.site_key.str[:5]
    t.loc[~hit & line.isin(list(ROUNDED_LINES)), "site_key"] = (
        line.map(ROUNDED_LINES) + t.site_key.str[5:])[~hit & line.isin(list(ROUNDED_LINES))]
    print(f"{label}: {int(hit.sum())} of {len(t)} tows matched to an ichthyoplankton net or CTD cast")
    return t.drop(columns="day")


def select(samples, release, label):
    """CalCOFI selection and sampling type, as the ichthyoplankton products."""
    s = samples.rename(columns={"latitude": "lat", "longitude": "lon"}).copy()
    s["grid_key"] = s.grid_key.fillna("")
    cruise = build.read_cruise_table(release)
    keep, per_cruise = build.select_calcofi_nets(s, cruise)
    s = s[keep & s.site_key.notna()].copy()
    hydro_sites = build.read_hydro_sites(release)
    s["sampling"] = build.classify_sampling(s, hydro_sites, per_cruise)
    # A tow with no cruise key (ZooScan, 2016 onward, unmatched: dated a day off
    # its cast) at a standard station is standard sampling.
    s.loc[s.cruise_key.isna() & s.site_key.isin(hydro_sites), "sampling"] = build.SAMPLING[0]
    dt = pd.to_datetime(s.datetime)
    s["date"] = dt.dt.strftime("%Y-%m-%d")
    s["year"], s["month"] = dt.dt.year, dt.dt.month
    s["season"] = s.month.map(SEASONS)
    s["season_year"] = s.year + (s.month == 12)
    s["ship_name"] = s.cruise_key.map(cruise.set_index("cruise_key").ship_name)
    s["line"] = s.site_key.str[:5].astype(float)
    s["station"] = s.site_key.str[6:].astype(float)
    print(f"{label}: {int(keep.sum())} of {len(keep)} tows selected; {s.sampling.value_counts().to_dict()}")
    return s


def write(df, name):
    path = os.path.join(OUT_DIR, name)
    con = duckdb.connect()
    con.register("df", df)
    con.sql(f"COPY df TO '{path}' (FORMAT parquet, COMPRESSION zstd)")
    print(f"  wrote {name}: {len(df)} rows ({os.path.getsize(path) / 1e6:.1f} MB)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--release", default="latest")
    args = ap.parse_args()
    release = build.resolve_release(args.release)
    url = lambda t: build.PARQUET_URL.format(release=release, table=t)
    obs = lambda dk: url(f"obs/dataset_key={dk}/data_0")
    print(f"CalCOFI.io release {release}")
    os.makedirs(OUT_DIR, exist_ok=True)
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs")
    sample_cols = "sample_key, cruise_key, grid_key, site_key, latitude, longitude, datetime, tow_type"

    # ---- plankton displacement volume (ichthyoplankton nets) ----
    v = con.sql(f"""
        SELECT s.{sample_cols.replace(', ', ', s.')},
               max(CASE WHEN m.measurement_type = 'total_plankton_biomass' THEN m.measurement_value END) AS total_ml,
               max(CASE WHEN m.measurement_type = 'small_plankton_biomass' THEN m.measurement_value END) AS small_ml,
               max(CASE WHEN m.measurement_type = 'volume_sampled' THEN m.measurement_value END) AS volume_m3
        FROM '{url("sample")}' s JOIN '{url("sample_measurement")}' m USING (sample_key)
        WHERE s.dataset_key = 'swfsc_ichthyo'
          AND m.measurement_type IN ('total_plankton_biomass', 'small_plankton_biomass', 'volume_sampled')
        GROUP BY ALL
        HAVING total_ml IS NOT NULL OR small_ml IS NOT NULL""").df()
    v = select(v, release, "plankton volume")
    v = v[v.volume_m3 > 0].copy()
    v["total_ml_per_1000m3"] = v.total_ml / v.volume_m3 * 1000
    v["small_ml_per_1000m3"] = v.small_ml / v.volume_m3 * 1000
    v = v.rename(columns={"tow_type": "net_type"})
    print(f"  nets: {v.net_type.value_counts().to_dict()}; total volume per 1,000 m3 median "
          f"{v.total_ml_per_1000m3.median():.0f} (C1 {v[v.net_type == 'C1'].total_ml_per_1000m3.median():.0f}, "
          f"CB {v[v.net_type == 'CB'].total_ml_per_1000m3.median():.0f})")
    write(v[TOW_COLUMNS + ["net_type", "volume_m3", "total_ml", "small_ml", "total_ml_per_1000m3", "small_ml_per_1000m3"]]
          .sort_values(["date", "sample_key"]), "plankton_volume_tows.parquet")

    # ---- names: the taxon table, else the dataset's own (ZooScan classes) ----
    names = con.sql(f"""
        SELECT dt.dataset_key, dt.taxon_key, coalesce(t.scientific_name, dt.ds_scientific_name) AS scientific_name,
               coalesce(t.common_name, dt.ds_common_name) AS common_name, dt.ds_scientific_name, t.rank
        FROM '{url("dataset_taxon")}' dt LEFT JOIN '{url("taxon")}' t USING (taxon_key)
        WHERE dt.dataset_key IN ('cce-lter_euphausiids', 'cce-lter_zooscan')""").df()

    nets = con.sql(f"""SELECT strftime(datetime, '%Y-%m-%d') AS day, latitude AS latitude_n, longitude AS longitude_n,
                              site_key AS site_key_n, cruise_key AS cruise_key_n
                       FROM '{url("sample")}' WHERE (dataset_key = 'swfsc_ichthyo' AND sample_type = 'net'
                                                      OR dataset_key = 'calcofi_ctd-cast')
                         AND site_key IS NOT NULL AND cruise_key IS NOT NULL""").df()     # nets, else CTD casts
    for dk, label in [("cce-lter_euphausiids", "krill"), ("cce-lter_zooscan", "zooscan")]:
        tows = con.sql(f"SELECT {sample_cols} FROM '{url('sample')}' WHERE dataset_key = '{dk}'").df()
        tows = select(match_ichthyo(tows, nets, label), release, label)
        o = con.sql(f"""SELECT sample_key, taxon_key, life_stage, measurement_type, measurement_value AS value
                        FROM '{obs(dk)}' WHERE measurement_value IS NOT NULL""").df()
        o = o[o.sample_key.isin(tows.sample_key)]
        nm = names[names.dataset_key == dk].drop_duplicates("taxon_key").set_index("taxon_key")
        if label == "krill":
            o["stage"] = o.life_stage.str.split().str[0]
            staged = o[o.stage.isin(KRILL_STAGES)].groupby(["sample_key", "taxon_key", "stage"]).value.sum().reset_index()
            total = o.groupby(["sample_key", "taxon_key"]).value.sum().reset_index().assign(stage="all")
            vals = pd.concat([staged, total], ignore_index=True)
            vals = vals[vals.value > 0]
            vals["measurement"] = "abundance_per_m2"
        else:
            vals = o[o.measurement_type.isin(ZOOSCAN_MEASUREMENTS)].rename(columns={"measurement_type": "measurement"})
            vals["stage"] = None
            vals = vals[["sample_key", "taxon_key", "stage", "measurement", "value"]]
        vals["scientific_name"] = vals.taxon_key.map(nm.scientific_name)
        vals["common_name"] = vals.taxon_key.map(nm.common_name)
        print(f"  {label}: {vals.taxon_key.nunique()} taxa; {vals.groupby(['measurement']).size().to_dict()}")
        write(tows[TOW_COLUMNS].sort_values(["date", "sample_key"]), f"{label}_tows.parquet")
        write(vals.sort_values(["sample_key", "taxon_key"]), f"{label}_values.parquet")


    # ---- CUFES: underway samples, each at its own position ----
    s = con.sql(f"""SELECT sample_key, cruise_key, grid_key, latitude, longitude, datetime
                    FROM '{url("sample")}' WHERE dataset_key = 'swfsc_cufes'""").df()
    s = s.rename(columns={"latitude": "lat", "longitude": "lon"})
    s["grid_key"] = "on track"                      # select by cruise and domain only
    cruise = build.read_cruise_table(release)
    keep, per_cruise = build.select_calcofi_nets(s, cruise)
    s = s[keep].copy()
    s["sampling"] = np.where(s.cruise_key.map(per_cruise.hydro_cruise).fillna(False).astype(bool),
                             build.SAMPLING[0], build.SAMPLING[1])
    dt = pd.to_datetime(s.datetime)
    s["date"] = dt.dt.strftime("%Y-%m-%d")
    s["year"], s["month"] = dt.dt.year, dt.dt.month
    s["season"] = s.month.map(SEASONS)
    s["season_year"] = s.year + (s.month == 12)
    s["ship_name"] = s.cruise_key.map(cruise.set_index("cruise_key").ship_name)
    print(f"cufes: {len(s)} samples selected on {s.cruise_key.nunique()} cruises; {s.sampling.value_counts().to_dict()}")
    o = con.sql(f"""SELECT sample_key, taxon_key, measurement_value AS value FROM '{obs("swfsc_cufes")}'
                    WHERE measurement_type = 'egg_concentration' AND measurement_value > 0""").df()
    o = o[o.sample_key.isin(s.sample_key)]
    nm = con.sql(f"SELECT taxon_key, scientific_name, common_name FROM '{url('taxon')}'").df().set_index("taxon_key")
    o["scientific_name"] = o.taxon_key.map(nm.scientific_name)
    o["common_name"] = o.taxon_key.map(nm.common_name)
    print(f"  cufes: {o.groupby('scientific_name').size().to_dict()} non-zero samples")
    write(s[["sample_key", "cruise_key", "ship_name", "lat", "lon", "date", "year", "month", "season", "season_year",
             "sampling"]].sort_values(["date", "sample_key"]), "cufes_samples.parquet")
    write(o.sort_values(["sample_key", "taxon_key"]), "cufes_values.parquet")


if __name__ == "__main__":
    main()
