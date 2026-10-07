# image_gridded_MUR_SST.py

import numpy as np
import xarray as xr
import holoviews as hv
import geoviews as gv
import geoviews.feature as gf

from cartopy import crs
from geoviews import opts

gv.extension('matplotlib')

gv.output(size=150)

xr_ensemble = xr.open_dataset('/data_7TB/mnt/data/dynamic_data/projects/projects2022/recreational_fishing_high_resolution/data/regional_downloads/latest_MUR_SST_Southwest_US.nc').load()
xr_ensemble

# NOTE: The daily MUR-SST variable is called analysed_sst, but the monthly MUR-SST variable is called sst.
kdims = ['time', 'longitude', 'latitude']
vdims = ['sst']

xr_dataset = gv.Dataset(xr_ensemble, kdims=kdims, vdims=vdims)

xr_dataset.to(gv.Image, ['longitude', 'latitude'])

max_sst = xr_dataset.range('sst')[1]
print(max_sst)
# use a new variable name for the  new assigned vartiable
sst_scaled = xr_dataset.redim.range(max_sst=(10, max_sst)).to(gv.Image,['longitude', 'latitude'])

sst_scaled.opts(opts.Image(projection=crs.Mercator(), cmap='Viridis', colorbar=True))