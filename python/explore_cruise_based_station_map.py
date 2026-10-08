# explore_cruise_based_station_map.py
"""
Exploratory map of how the full swfsc_ichthyo.nc net-tow dataset is
filtered to the CalCOFI selection, coloured by the CalCOFI.io `grid`
table's station pattern.

Selection rule (applied per tow; all inputs from CalCOFI.io parquet):
  1. Hydrographic CalCOFI cruise: the cruise also contributed CalCOFI
     bottle or CTD data (`cruise.cruise_key_datasets`).
  2. Restored early/other CalCOFI-ship cruise: the cruise has no bottle/CTD
     data in the database, BUT its ship ran at least one cruise of type 1
     AND >= 90% of its tows sit on grid stations inside the domain below.
     (Robust: 91/90/89 cruises restored at an 80/90/100% threshold; the
     CalCOFI-ship cruises left out have a median of 7% in-domain tows.)
  3. Domain: south of the California-Oregon border (42N), outside the Gulf
     of California, within DOMAIN_BOX.
A tow is selected if its cruise passes 1 or 2 AND the tow itself is in
the domain AND on a grid station (has a grid_key). Everything else is
dropped. The left panel's frame is then fitted to the selected tows.

Station pattern: swfsc_ichthyo.nc (release v2026.09.06) uses the previous
grid's keys; the current release renamed/redrew the grid, so keys are
mapped through `grid_crosswalk` (each old key -> its best-overlapping new
station). Tows with no grid_key at all are "not on grid".

Also writes the per-cruise selection (with reason) as CSV for the
gridded-NetCDF build step.
"""
import duckdb
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
import cartopy.crs as ccrs
import cartopy.feature as cfeature

RELEASE = "v2026.10.06"
PARQUET = f"https://storage.googleapis.com/calcofi-db/ducklake/releases/{RELEASE}/parquet"
BASE = "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas"
SRC = f"{BASE}/data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io/swfsc_ichthyo.nc"
OUT = f"{BASE}/figures/exploration/cruise_based_station_map.png"
OUT_CSV = f"{BASE}/data/CalCOFI_ichthyoplankton/ichthyoplankton_from_calCOFI.io/calcofi_cruise_selection.csv"
DOMAIN_BOX = [-140, -105, 18, 50]   # lon_min, lon_max, lat_min, lat_max
FRAME_PAD_DEG = 1.0
CA_OR_BORDER_LAT = 42.0
RESTORE_MIN_IN_DOMAIN = 0.9

# Gulf of California: east of the Baja peninsula's Gulf coast, closed across
# the mouth from Cabo San Lucas to Cabo Corrientes (IHO southern limit).
GULF_OF_CALIFORNIA = Path([
    (-115.0, 33.0), (-115.0, 32.0), (-114.6, 30.0), (-113.2, 28.0),
    (-111.8, 26.0), (-110.3, 24.0), (-109.9, 22.9), (-105.7, 20.4),
    (-104.0, 20.4), (-104.0, 33.0),
])

# Fixed categorical order (validated: CVD dE >= 9); marker shape is the
# secondary encoding.
PATTERNS = [
    ("standard", "Standard", "#2a78d6", "o"),
    ("extended", "Extended", "#eb6834", "s"),
    ("historical", "Historical", "#1baf7a", "^"),
    ("not on grid", "Not on grid", "#eda100", "D"),
]

con = duckdb.connect()
con.sql("INSTALL httpfs; LOAD httpfs;")
cruise = con.sql(f"""SELECT cruise_key, cruise_key_datasets, ship_name
                     FROM read_parquet('{PARQUET}/cruise.parquet')""").df()
crosswalk = con.sql(f"""
    SELECT prev_grid_key, arg_max(grid_key, prev_frac) AS grid_key
    FROM read_parquet('{PARQUET}/grid_crosswalk.parquet') GROUP BY prev_grid_key""").df()
grid = con.sql(f"SELECT grid_key, pattern FROM read_parquet('{PARQUET}/grid.parquet')").df()

net = netCDF4.Dataset(SRC)["net"]
tows = pd.DataFrame({
    "cruise_key": netCDF4.chartostring(net["cruise_key"][:]),
    "prev_grid_key": netCDF4.chartostring(net["grid_key"][:]),
    "lat": np.ma.filled(net["latitude"][:], np.nan),
    "lon": np.ma.filled(net["longitude"][:], np.nan),
})
tows = (tows.merge(cruise, on="cruise_key", how="left")
            .merge(crosswalk, on="prev_grid_key", how="left")
            .merge(grid, on="grid_key", how="left"))
tows["pattern"] = tows.pattern.fillna("not on grid")

lon_min, lon_max, lat_min, lat_max = DOMAIN_BOX
tows["in_domain"] = (
    tows.lat.between(lat_min, CA_OR_BORDER_LAT) & tows.lon.between(lon_min, lon_max)
    & ~GULF_OF_CALIFORNIA.contains_points(tows[["lon", "lat"]].to_numpy())
)
tows["hydro_cruise"] = tows.cruise_key_datasets.fillna("").str.contains("calcofi_bottle|calcofi_ctd-cast")

calcofi_ships = set(tows.loc[tows.hydro_cruise, "ship_name"].dropna())
per_cruise = tows.groupby("cruise_key").agg(
    ship_name=("ship_name", "first"), datasets=("cruise_key_datasets", "first"),
    hydro_cruise=("hydro_cruise", "first"), tows=("lat", "size"),
    frac_in_domain_on_grid=("in_domain", lambda s: (s & (tows.loc[s.index, "prev_grid_key"] != "")).mean()),
)
restored = (~per_cruise.hydro_cruise & per_cruise.ship_name.isin(calcofi_ships)
            & (per_cruise.frac_in_domain_on_grid >= RESTORE_MIN_IN_DOMAIN))
per_cruise["reason"] = np.select(
    [per_cruise.hydro_cruise, restored],
    ["hydrographic CalCOFI cruise", "restored: CalCOFI ship, on-grid in domain"],
    default="dropped")
per_cruise.to_csv(OUT_CSV)

tows = tows.join(per_cruise.reason, on="cruise_key")
tows["selected"] = (tows.reason != "dropped") & tows.in_domain & (tows.pattern != "not on grid")
sel, drop = tows[tows.selected], tows[~tows.selected]
n_hydro = sel[sel.reason == "hydrographic CalCOFI cruise"].cruise_key.nunique()
n_rest = sel[sel.reason != "hydrographic CalCOFI cruise"].cruise_key.nunique()

# Left panel fitted to the selected tows; right panel spans every dropped
# tow, rounded out to 5 degrees.
selected_extent = [np.floor(sel.lon.min() - FRAME_PAD_DEG), np.ceil(sel.lon.max() + FRAME_PAD_DEG),
                   np.floor(sel.lat.min() - FRAME_PAD_DEG), np.ceil(sel.lat.max() + FRAME_PAD_DEG)]
dropped_extent = [max(-180, 5 * np.floor(drop.lon.min() / 5)), 5 * np.ceil(drop.lon.max() / 5),
                  5 * np.floor(drop.lat.min() / 5), 5 * np.ceil(drop.lat.max() / 5)]

panels = [
    (sel, selected_extent, 4,
     f"Selected: {n_hydro + n_rest} CalCOFI cruises ({n_hydro} with hydrographic data + {n_rest} filtered by CalCOFI ship and station pattern)\n"
     f"{len(sel):,} tows on grid stations; south of 42°N, outside the Gulf of California"),
    (drop, dropped_extent, 6,
     f"Dropped: {drop.cruise_key.nunique()} cruises with tows outside the selection "
     f"(other surveys, Gulf of California, north of 42°N, not on grid)\n{len(drop):,} tows"),
]
widths = [(e[1] - e[0]) / (e[3] - e[2]) for _, e, _, _ in panels]
fig = plt.figure(figsize=(21, 8.5))
gs = fig.add_gridspec(1, 2, width_ratios=widths)
for i, (d, extent, size, title) in enumerate(panels):
    ax = fig.add_subplot(gs[0, i], projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="#e8e6df", zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor="#6b6a64")
    ax.add_feature(cfeature.BORDERS, linewidth=0.4, edgecolor="#9a9890")
    for key, label, color, marker in PATTERNS:
        p = d[d.pattern == key]
        if p.empty:
            continue
        ax.scatter(p.lon, p.lat, s=size, alpha=0.4, color=color, marker=marker, linewidths=0,
                   transform=ccrs.PlateCarree(), label=f"{label} ({len(p):,})")
    ax.set_title(title, fontsize=10.5)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="#bdbbb3", alpha=0.7)
    gl.top_labels = gl.right_labels = False
    leg = ax.legend(title="Station pattern", loc="lower left", markerscale=3, fontsize=9, framealpha=0.9)
    for h in leg.legend_handles:
        h.set_alpha(1)
fig.suptitle(f"SWFSC ichthyoplankton net tows: full dataset filtered to the CalCOFI selection "
             f"(CalCOFI.io {RELEASE})", fontsize=13)
fig.tight_layout()
fig.savefig(OUT, dpi=130)
print("wrote", OUT, "and", OUT_CSV)
print(per_cruise.reason.value_counts())
print("tows selected", len(sel), "dropped", len(drop), "extents", selected_extent, dropped_extent)
