# prepare_larval_plot.py

""" 
Script to load larval data, write geotiff metadata, and plot a map with labels and colorbar.

-- Loads data from netcdf
-- Uses rioxarray to write the geotiff metadata and appy a projection
-- Saves a geotiff

Output
A tif file  "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/" + region + "_larval_reference.tif"
and 
../figures/larval_plot_"+ region+ ".tif"
"""

import os
from sys import argv

import pickle
from datetime import datetime, timedelta
from pytz import timezone

import cartopy
import cartopy.feature as cfeature
import cartopy.crs as ccrs

import colorcet as cc
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import rioxarray as rxr
import xarray as xr

# Pass in the region list variable from shell script
REGION_STR = argv[1]
REGION_LIST = REGION_STR.split(",")
DAYADJ = argv[2]

for region in REGION_LIST:

    DAYADJ = int(DAYADJ)

    ######################################
    # load data from a netCDF data file using xarray
    ######################################

    DATADIR = "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/"
    f = "latest_MUR_SST_" + region + ".nc"
    d = DATADIR + f

    ds = xr.open_dataset(d)
    xr_sst = ds["analysed_sst"]
    # load the region directory 
    import region_dictionary
    region_bounds = region_dictionary.get_region_bounds(region)
    ######################################
    # Add geotiff metadata to a reference file
    ######################################

    # remove the time coordinate  to match the number of dimensions (2) that chlor_a has.
    xr_sst = xr_sst.squeeze()

    # Create the reference geotiff

    lat = ds["latitude"]
    lon = ds["longitude"]

    xr_sst.rio.write_crs("epsg:4326", inplace=True)
    xr_sst.rio.to_raster(
        "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/"
        + region
        + "_sst_reference.tif")

    # copy the sst reference tif file to a separate directory 
    # so these region files can be used for Sobel gradient calculations
    os.system("cp /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/" + region + "_sst_reference.tif /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_for_sobel_calculation/" + region + "_sst_reference.tif")


    ######################################
    # load the 4-day regional max/min SST calculated from the data
    ######################################

    # Import the saved max and min SST from all 4 days combined
    with open("sst_range.bin", "rb") as f: # "rb" because we want to read in binary mode
        sst_range = pickle.load(f)

    min_4_day_sst = sst_range[0]
    max_4_day_sst = sst_range[1]

    # convert SST range to degrees F if region is in California or Baja
    if region_bounds["super_region"] == "Southwest_US" or region_bounds["super_region"] == "Baja_Mexico_only":
        min_4_day_sstF = (min_4_day_sst * 1.8) + 32
        max_4_day_sstF = (max_4_day_sst * 1.8) + 32
        # Reassign sst range
        min_4_day_sst = min_4_day_sstF
        max_4_day_sst = max_4_day_sstF

    print("%%%%%%%%%%%%%%%%%%%%%")
    print("min_4 day sst", min_4_day_sst)
    print("max_4 day sst", max_4_day_sst)
    print("%%%%%%%%%%%%%%%%%%%%%")

    ######################################
    # plot the map
    ######################################
   
    # Load the reference geotiff
    sst_3d = rxr.open_rasterio(
        "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/"
        + region
        + "_sst_reference.tif",
        masked=True)
    sst = sst_3d.squeeze()

    # convert SST to degrees F if region is in California or Baja
    if region_bounds["super_region"] == "Southwest_US" or region_bounds["super_region"] == "Baja_Mexico_only":
        sstF = (sst * 1.8) + 32
        # Reassign sst
        sst = sstF
    
    
    # Note: the projection is ultimately is set with GDAL warp.
    ax = plt.axes()
    # ax = plt.axes(projection=ccrs.Mercator()) # this is EPSG:4326

    # # NOTE The cloud-free GHRSST product is acquired by the satellite at 08:00h UTC and  uploaded to the Copernicus Marine database by 11:00h UTC. 

    # Add 1 day to  file date (in GMT) for NZ. Cal/ Baja keeps the same day as GMT date. 
 
    out = ds["time"].dt.strftime("%Y %m %d %H %M %S")
    
    dl = np.fromstring(out.item(0),dtype=int, sep=" ").tolist()
    utc = timezone("UTC")
    date_object = utc.localize(datetime(dl[0],dl[1],dl[2],dl[3],dl[4],dl[5]))

    # Time zone conversions
    # NOTE: All map region timezones in region_dictionary.py are set to GMT since the data times are UTC
    tz = timezone(region_bounds["time_zone"])
    local_time = date_object.astimezone(tz)

    # NOTE: I've removed the timestamp so just the day shows in the label
    # The \n in 'string1\nstring2' should put a carriage return between the strings

    if region_bounds["super_region"] == "New_Zealand":
        label_time =  local_time + timedelta(days=1)
        label = datetime.strftime(label_time, "%e %b") + " SST" + r" ($^\circ$$C$)"  
    elif region_bounds["super_region"] == "Southwest_US": 
        label = datetime.strftime(local_time, "%e %b") + " SST" + r" ($^\circ$$F$)"  
    elif region_bounds["super_region"] == "Baja_Mexico_only": 
        label = datetime.strftime(local_time, "%e %b") + " SST" + r" ($^\circ$$F$)"  
    elif region_bounds["super_region"] == "Central_America": 
        label_time =  local_time + timedelta(days=1)
        label = datetime.strftime(label_time, "%e %b") + " SST" + r" ($^\circ$$C$)"  
    elif region_bounds["super_region"] == "Mauritius_Reunion": 
        label_time =  local_time + timedelta(days=1)
        label = datetime.strftime(label_time, "%e %b") + " SST" + r" ($^\circ$$C$)"  
    elif region_bounds["super_region"] == "Western_Australia": 
        label_time =  local_time + timedelta(days=1)
        label = datetime.strftime(label_time, "%e %b") + " SST" + r" ($^\circ$$C$)"  
    elif region_bounds["super_region"] == "Eastern_Australia": 
        label_time =  local_time + timedelta(days=1)
        label = datetime.strftime(label_time, "%e %b") + " SST" + r" ($^\circ$$C$)"  
     
    FULL_LABEL = str(label)   

    # Colorbar settings 
    # label positioning
    X =region_bounds["colorbar_position_X"]
    Y =region_bounds["colorbar_position_Y"]

    ax.text(
        X,
        Y,
        FULL_LABEL,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        color="white",
        bbox={"facecolor": "black", "alpha": 1.0, "pad": 2},
    )

    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    color_map = cc.m_rainbow_bgyrm_35_85_c69  # a colorcet colormap
    # color_map =cc.cm["gouldian"] # a colorcet colormap
    # color_map ="cmc.lapaz" # a crameri colormap, reversed
    # color_map ="cmc.imola" # a crameri colormap, reversed
    # color_map ="cmc.romaO_r" # a crameri colormap, reversed
    # color_map ="cmc.vikO" # a crameri colormap, reversed
    # color_map ="cmc.roma_r" # a crameri colormap, reversed
    # color_map ="cmc.hawaii_r" # a crameri colormap, reversed
    # color_map ="plasma"
    # color_map ="rainbow"
    # color_map ="viridis"
    # color_map ="magma"
    # color_map = "PuBu_r"
    # color_map ="cmo.thermal"
    # color_map ="cmo.ice"
    # color_map ="cmo.dense"

    # Set contour levels based on the ranges of SST values
    upper_contour_limit = np.ceil(max_4_day_sst)
    lower_contour_limit = np.floor(min_4_day_sst)
    # contour interval set in the region directory for best effect
    c_interval = region_bounds["SST_contour_interval"]
    contour_levels = np.arange(
        lower_contour_limit, upper_contour_limit, c_interval
    ).tolist()

    # set the color range based on the max and min for the region before subsetting by day so that all of the colourbars in the 4 day plots will be the same. You don't want to do this for the pooled regions because you will lose detail.

    pt_sst = plt.pcolormesh(
        lon,
        lat,
        sst,
        shading="auto", 
        cmap=color_map,
        alpha=0.85,
        vmin=min_4_day_sst,
        vmax=max_4_day_sst,
        zorder=0)

    # add contours for "temperature breaks" -------------------
    CS = plt.contour(
        lon, lat, sst, levels=contour_levels, linewidths=0.5, colors="white")
    ax.clabel(
        CS,
        CS.levels,
        inline=True,
        fontsize=5)

    
    # load bathymetric data
    bathy = xr.open_dataset("/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/bathymetry/"+region+"_pygmt_earth_relief_bathymetry.nc")
    depth = bathy["elevation"]
    lat = bathy["lat"]
    lon = bathy["lon"]
    
    # Set bathymetric contour levels based on the ranges of depth values
    bt_upper_contour_limit = 0
    bt_lower_contour_limit = -3200

    # unevenly-spaced contour intervals --------------------- 
 
    if all([region !="North_Island_NZ", region !="South_Island_NZ", region !="Mauritius", region !="Reunion", region !="Ilwaco_Westport_WA", region !="Rottnest_Abrolhos", region !="Exmouth"]):
        bt_c_interval = [-3000,-2500,-2000,-1800,-1600,-1400,-1200,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100]
        # was [-3000,-2500,-2000,-1800,-1600,-1400,-1200,-1000,-800,-600,-400,-200,-150,-100,-50]
        topo_c_interval = [0,4000]
    elif all([region == "North_Island_NZ", region =="Soudan_Banks"]):
        bt_c_interval = [-4000,-3000,-2000,-1000,-800,-600,-400,-200,-100,-50,0]
        topo_c_interval = [0,4000]
    else:
        bt_c_interval = [-4000,-3000,-2000,-1000,0]
        topo_c_interval = [0,4000]
        
    # evenly-spaced contour intervals ---------------------
    # # contour interval set in the region directory for best effect
    # bt_c_interval = region_bounds["topo_contour_interval"]
    # bt_contour_levels = np.arange(
    #     bt_lower_contour_limit, bt_upper_contour_limit, bt_c_interval
    # ).tolist()
    # ----------------------------------------------------

    BT = plt.contour(
        lon, lat, depth, levels=bt_c_interval, linewidths=0.9, linestyles="solid", colors="black")
        # was linewidths=0.75
    ax.clabel(
        BT,
        BT.levels,
        inline=True,
        fontsize=5) # was 4
    # # second label to create 2 labels per line
    # ax.clabel(
    #     BT,
    #     BT.levels,
    #     inline=True,
    #     fontsize=4)

    Topo = plt.contourf(
        lon, lat, depth, levels=topo_c_interval, linestyles="solid", colors="black", zorder=2)
    #artists with higher zorder are drawn on top


    # NOTE: The reason I put the colorbar and label inside
    # the plot is because placing them outside screws up the overlay  with the bathymetry.

    # sets the size  of the background box that contains the colorbar
    # proportion of percentage of the frame width
    cbbox = inset_axes(
        ax, width='16%', height='40%', loc=region_bounds["colorbar_position"])
        
    cbbox.tick_params(
        axis="both",
        left=False,
        top=False,
        right=False,
        bottom=False,
        labelleft=False,
        labeltop=False,
        labelright=False,
        labelbottom=False,
    )
    cbbox.set_facecolor([1,1,1,1])
   
   # colorbar size in the background box
   # e.g height spans 90% of the height of the background box
    cbaxes = inset_axes(cbbox, "30%", "90%", loc=region_bounds["colorbar_position"])
   
    cb=plt.colorbar(pt_sst, cax=cbaxes) #make colorbar
    
    # adjusts which side of the colorbar the labels are on 
    if region_bounds["tick_location"] == "left":
        cbaxes.yaxis.tick_left()

    cb.outline.set_edgecolor('black')
    cb.outline.set_linewidth(2)

    cb.ax.tick_params(labelsize=7, colors="black")

    plt.clim(lower_contour_limit, upper_contour_limit)

    cb.solids.set(alpha=1)
    # Save the plot in /figures.

    # GDAL translate operates on this file to first georeference it
    # using the corners from the reference geotiff file in /data,
    # and then to warp it to a spherical surface using the projection.
    # The final product is the /figures/[regional filename]_proj.tif

    plt.savefig("/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/SST_plot_"+ region+ ".tif",bbox_inches="tight", pad_inches=0, dpi=600)
