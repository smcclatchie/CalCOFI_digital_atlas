# test_holoviews.py

import numpy as np

import xarray as xr

import holoviews as hv
from holoviews import opts
import plotly.graph_objects as go

import pandas as pd
import xarray as xr

hv.extension('plotly')
# hv.extension('matplotlib')

# load one cruise of the dataset
ds = pd.read_csv('/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/merged_2020_cruise_depth.csv')
print(ds.head())

# # plot station locations
stations = ds[['Lon_Dec', 'Lat_Dec', 'R_Depth']].drop_duplicates()
# stations_hv = hv.Points(stations, kdims=['Lon_Dec', 'Lat_Dec'])
# stations_hv.opts(opts.Points(size=5, color='blue', tools=['hover']))

# decalare dedependent variables
# vdims = [('R_TEMP', 'Temperature (°C)')]
# cruise3D = hv.Dataset(ds, ['Lon_Dec', 'Lat_Dec'], vdims=vdims)
lon = stations['Lon_Dec']
lat = stations['Lat_Dec']
depth = stations['R_Depth']
temperature = ds['R_TEMP']

# plot the data as a 2D contour plot
go.Figure(data = 
                go.Contour(x=lon, y=lat, z=temperature, colorscale='Viridis'))
# fig.show()  
