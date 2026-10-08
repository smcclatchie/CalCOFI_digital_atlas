# explore_cruise_based_station_map.py
"""
Map of how the full swfsc_ichthyo.nc net-tow dataset is filtered to the
CalCOFI selection, coloured by the CalCOFI.io `grid` table's station
pattern. Writes three separate maps, 300 dpi for print: (1) every tow,
kept (navy) vs dropped (amber), at the full extent of the dropped tows, (2) the CalCOFI selection (fitted to the kept tows), and (3) the
same selection split into standard CalCOFI sampling vs high-resolution /
special surveys.

The selection and the standard/special-survey split are
build_gridded_abundance_by_cruise_and_taxon.py's select_calcofi_nets() and
classify_sampling(), so this figure and the gridded NetCDF always agree --
see that script's docstring for the rules. (The NetCDF additionally keeps
only the starboard net of paired Pairovet tows; both nets sit at the same
position, so the map is unaffected.)

Station pattern: swfsc_ichthyo.nc (release v2026.09.06) uses the previous
grid's keys; the current release renamed/redrew the grid, so keys are
mapped through `grid_crosswalk` (each old key -> its best-overlapping new
station). Tows with no grid_key at all are "not on grid".

Also writes the per-cruise decision (with reason) as CSV.
"""
import os

import duckdb
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from build_gridded_abundance_by_cruise_and_taxon import (
    DEFAULT_RELEASE, DEFAULT_SRC, PARQUET_URL, SAMPLING, classify_sampling, read_cruise_table,
    read_hydro_sites, select_calcofi_nets,
)

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root
OUT_DIR = f"{BASE}/figures/exploration"
os.makedirs(OUT_DIR, exist_ok=True)       # figures/ is gitignored, absent in a fresh clone
OUT_NAMES = ["station_map_1_kept_and_dropped.png", "station_map_2_selected_by_station_pattern.png",
             "station_map_3_selected_by_sampling_type.png"]
OUT_CSV = f"{BASE}/data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io/calcofi_cruise_selection.csv"
FRAME_PAD_DEG = 1.0
FIG_WIDTH_IN = 13

# Fixed categorical order (validated: CVD dE >= 9); marker shape is the
# secondary encoding.
PATTERNS = [
    ("standard", "Standard", "#2a78d6", "o"),
    ("extended", "Extended", "#eb6834", "s"),
    ("historical", "Historical", "#1baf7a", "^"),
    ("not on grid", "Not on grid", "#eda100", "D"),
]
# Sampling-type panel: separate hues from the station-pattern panels so the
# two encodings are never confused (validated: CVD and contrast pass).
SAMPLING_STYLE = {
    SAMPLING[0]: ("Standard CalCOFI sampling", "#184f95", "o"),
    SAMPLING[1]: ("High-resolution / special survey", "#d55181", "s"),
}

con = duckdb.connect()
con.sql("INSTALL httpfs; LOAD httpfs;")
crosswalk = con.sql(f"""
    SELECT prev_grid_key AS grid_key, arg_max(grid_key, prev_frac) AS current_grid_key
    FROM read_parquet('{PARQUET_URL.format(release=DEFAULT_RELEASE, table="grid_crosswalk")}')
    GROUP BY prev_grid_key""").df()
grid = con.sql(f"""SELECT grid_key AS current_grid_key, pattern
                   FROM read_parquet('{PARQUET_URL.format(release=DEFAULT_RELEASE, table="grid")}')""").df()

net = netCDF4.Dataset(DEFAULT_SRC)["net"]
tows = pd.DataFrame({
    "cruise_key": netCDF4.chartostring(net["cruise_key"][:]),
    "grid_key": netCDF4.chartostring(net["grid_key"][:]),
    "site_key": netCDF4.chartostring(net["site_key"][:]),
    "lat": np.ma.filled(net["latitude"][:].astype("f8"), np.nan),
    "lon": np.ma.filled(net["longitude"][:].astype("f8"), np.nan),
})
tows["selected"], per_cruise = select_calcofi_nets(tows, read_cruise_table(DEFAULT_RELEASE))
per_cruise.to_csv(OUT_CSV)

tows = (tows.join(per_cruise.reason, on="cruise_key")
            .merge(crosswalk, on="grid_key", how="left")
            .merge(grid, on="current_grid_key", how="left"))
tows["pattern"] = tows.pattern.fillna("not on grid")
sel, drop = tows[tows.selected], tows[~tows.selected].copy()
sel = sel.assign(sampling=classify_sampling(sel, read_hydro_sites(DEFAULT_RELEASE), per_cruise))
n_survey_cruises = int((sel.groupby("cruise_key").sampling.agg(lambda s: (s == SAMPLING[0]).sum()) == 0).sum())
n_hydro = sel[sel.reason == "hydrographic CalCOFI cruise"].cruise_key.nunique()
n_rest = sel[sel.reason != "hydrographic CalCOFI cruise"].cruise_key.nunique()

# Selected panel fitted to the kept tows; dropped panel spans every dropped
# tow, rounded out to 5 degrees.
selected_extent = [np.floor(sel.lon.min() - FRAME_PAD_DEG), np.ceil(sel.lon.max() + FRAME_PAD_DEG),
                   np.floor(sel.lat.min() - FRAME_PAD_DEG), np.ceil(sel.lat.max() + FRAME_PAD_DEG)]
dropped_extent = [max(-180, 5 * np.floor(drop.lon.min() / 5)), 5 * np.ceil(drop.lon.max() / 5),
                  5 * np.floor(drop.lat.min() / 5), 5 * np.ceil(drop.lat.max() / 5)]

panels = [
    (tows, dropped_extent, 6,
     f"All {len(tows):,} tows: {len(sel):,} kept, {len(drop):,} dropped "
     f"(other surveys, Gulf of California, north of 42°N, not on grid)"),
    (sel, selected_extent, 4,
     f"Selected: {n_hydro + n_rest} CalCOFI cruises ({n_hydro} with hydrographic data + {n_rest} "
     f"filtered by CalCOFI ship and station pattern)\n"
     f"{len(sel):,} tows on grid stations; south of 42°N, outside the Gulf of California"),
    (sel, selected_extent, 4,
     f"Selected tows by sampling type: standard CalCOFI stations vs high-resolution / special surveys\n"
     f"{n_survey_cruises} cruises are entirely special survey (e.g. 1980–86 anchovy, 2002–04 cowcod)"),
]
release_note = f"SWFSC ichthyoplankton net tows, CalCOFI.io {DEFAULT_RELEASE}"
for i, ((d, extent, size, title), name) in enumerate(zip(panels, OUT_NAMES)):
    aspect = (extent[3] - extent[2]) / (extent[1] - extent[0])
    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_WIDTH_IN * aspect + 1.4))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="#e8e6df", zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor="#6b6a64")
    ax.add_feature(cfeature.BORDERS, linewidth=0.4, edgecolor="#9a9890")
    if i == 0:
        legend_title = "Selection"
        # Navy/amber: blue-yellow contrast survives common colour-vision
        # deficiencies (validated: worst-case CVD dE 30, contrast 3:1).
        groups = [(d[d.selected], f"Kept, {sel.cruise_key.nunique()} cruises", "#184f95", "o"),
                  (d[~d.selected], f"Dropped, {drop.cruise_key.nunique()} cruises", "#c98500", "x")]
    elif i == 1:
        legend_title = "Station pattern"
        groups = [(d[d.pattern == key], label, color, marker) for key, label, color, marker in PATTERNS]
    else:
        legend_title = "Sampling type"
        groups = [(d[d.sampling == key], *SAMPLING_STYLE[key]) for key in SAMPLING]
    # In the sampling map, special surveys are drawn underneath so the
    # standard station grid stays visible; legend order is unchanged.
    # Kept tows on top of dropped ones (dropped tows at kept stations would
    # otherwise hide the kept pattern); special surveys underneath standard.
    zorders = {0: [3, 2], 2: [3, 2]}.get(i, [2] * len(groups))
    for (p, label, color, marker), z in zip(groups, zorders):
        if p.empty:
            continue
        ax.scatter(p.lon, p.lat, s=size, alpha=0.4, color=color, marker=marker,
                   linewidths=0.6 if marker == "x" else 0,
                   transform=ccrs.PlateCarree(), label=f"{label} ({len(p):,} tows)", zorder=z)
    ax.set_title(f"{title}\n{release_note}", fontsize=11)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="#bdbbb3", alpha=0.7)
    gl.top_labels = gl.right_labels = False
    leg = ax.legend(title=legend_title, loc="lower left", markerscale=3, fontsize=10, framealpha=0.9)
    for h in leg.legend_handles:
        h.set_alpha(1)
    # Leave room on the left for gridliner labels, which tight_layout ignores.
    fig.tight_layout(rect=[0.05, 0.0, 1.0, 1.0])
    fig.savefig(f"{OUT_DIR}/{name}", dpi=300, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print("wrote", f"{OUT_DIR}/{name}")
print("wrote", OUT_CSV)
print(per_cruise.reason.value_counts())
print("tows selected", len(sel), "dropped", len(drop), "extents", selected_extent, dropped_extent)
print("sampling:", sel.sampling.value_counts().to_dict(), "; entirely special-survey cruises:", n_survey_cruises)
