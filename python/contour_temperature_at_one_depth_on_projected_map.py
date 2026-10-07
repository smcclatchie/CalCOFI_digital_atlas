# contour_temperature_at_one_depth_on_projected_map.py

import numpy as np
import xarray as xr
import pandas as pd

import holoviews as hv
from holoviews import opts

import geoviews as gv
import geoviews.feature as gf
from geoviews import opts

from cartopy import crs

import plotly.graph_objects as go

hv.extension('bokeh')
# hv.extension('plotly')
# hv.extension('matplotlib')

# load one cruise of the dataset
ds = pd.read_csv('/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/merged_2020_cruise_depth.csv')
print(ds.head())

# # plot station locations
stations = ds[['Lon_Dec', 'Lat_Dec', 'R_Depth']].drop_duplicates()
# stations_hv = hv.Points(stations, kdims=['Lon_Dec', 'Lat_Dec'])
# stations_hv.opts(opts.Points(size=5, color='blue', tools=['hover']))

lon = stations['Lon_Dec']
lat = stations['Lat_Dec']
depth = stations['R_Depth']
temperature_50m = ds['R_TEMP']

lon2d, lat2d = np.meshgrid(lon, lat)

dataset1 = gv.Dataset((lon, lat, temperature_50m), kdims=['Longitude', 'Latitude'], vdims=['Temperature'])

# plot the data as a 2D contour plot
# go.Figure(data = 
                # go.Contour(x=lon, y=lat, z=temperature, colorscale='Viridis'))
# fig.show()  
