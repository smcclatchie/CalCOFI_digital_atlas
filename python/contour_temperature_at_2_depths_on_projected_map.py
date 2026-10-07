# contour_temperature_at_2_depths_on_projected_map.py

import numpy as np
import xarray as xr
import pandas as pd

import holoviews as hv
from holoviews import opts

import geoviews as gv
import geoviews.feature as gf
from geoviews import opts

from cartopy import crs as ccrs

import plotly.graph_objects as go

hv.extension('bokeh')
# hv.extension('plotly')
# hv.extension('matplotlib')

# load 2 depths from one cruise of the dataset
ds1 = pd.read_csv('/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/merged_2020_cruise_depth_50m.csv')

ds2 = pd.read_csv('/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/merged_2020_cruise_depth_100m.csv')

# # plot station locations
stations = ds1[['Lon_Dec', 'Lat_Dec']].drop_duplicates()
# stations = ds1[['Lon_Dec', 'Lat_Dec', 'R_Depth']].drop_duplicates()
# stations_hv = hv.Points(stations, kdims=['Lon_Dec', 'Lat_Dec'])
# stations_hv.opts(opts.Points(size=5, color='blue', tools=['hover']))

lon = stations['Lon_Dec']
lat = stations['Lat_Dec']
# depth = stations['R_Depth']
temperature_50m = ds1['R_TEMP']
temperature_100m = ds2['R_TEMP']

lon2d, lat2d = np.meshgrid(lon, lat)

z1 = np.sin(np.deg2rad(lon2d)) * np.cos(np.deg2rad(lat2d))
z2 = np.cos(np.deg2rad(lon2d)) * np.sin(np.deg2rad(lat2d))

dataset1 = gv.Dataset((lon, lat, z1), kdims=[lon, lat], vdims=[temperature_50m])
dataset2 = gv.Dataset((lon, lat, z2), kdims=[lon, lat], vdims=[temperature_100m])

contour1 = dataset1.to(gv.Contours, kdims=[lon, lat], vdims=[temperature_50m]).opts(
    cmap='Magma', alpha=0.8, colorbar=True
)

contour2 = dataset2.to(gv.Contours, kdims=[lon, lat], vdims=[temperature_100m]).opts(
    cmap='Blues', alpha=0.6, colorbar=True
)

stacked_map = contour1 * contour2
# stacked_map = gf.ocean * gf.coastline * contour1 * contour2

# stacked_map.opts(
#     projection=ccrs.PlateCarree(), 
#     global_extent=False, 
#     width=600, 
#     height=600


