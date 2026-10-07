#!/bin/bash

# master_shell_script_CalCOFI_larvae.sh

# USAGE:
# First download the data file:
# sh ./download_CalCOFI_larval_data.sh

# then from the shell directory:
# ./master_shell_script_CalCOFI_larvae

##############################################
# NOTE:
# You need to first activate one of the environments with 
# Python 3.n installed.
# e.g. cd ../shell
# conda activate geo_env
#############################################

# # pass in the arguments
# superregion=$1

# # Error trapping
# # test if the required data file downloaded
# if [[ $1 == "Southwest_US" ]]
# then
#     echo $1
#     test -f "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/latest_MUR_SST_Southwest_US.nc";
#     if [ $? -ne 0 ]; then
#         echo "Missing data file: No cloud-free SST data for Southwest US from ERDDAP"
#         cp /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/latest_MUR_SST_Southwest_US_NASA_final.nc /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/latest_MUR_SST_Southwest_US.nc
#     fi
# fi

# Delete previous figures and data
##################################################
rm /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/*.*
rm /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/*.*

# Create the SST maps
########################################################
conda run -n calcofi python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/prepare_larval_reference_geotiff_and_plot.py" $regions $dayadj

# NOTE:
# Output here is a tif file
# ../figures/SST_plot_"+ region + ".tif"

# Enhance the SST map 
####################################################
conda run -n image python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/batch_sharpen_SST_maps.py" $regions

# Use image magick to change resolution of enhanced files from 72dpi to 600dpi
# Use '2>/dev/null' to remove warning about tif tags
# Keep this or your SST images will be low resolution!
mogrify -set density 600 "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/SST_plot_"$regions"_clahe.tif" 

# Calculate the Sobel gradients
###############################################
conda run -n wbt python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/sobel_gradient.py" $regions

# Create the Sobel map
###############################################
conda run -n geo_env python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/sobel_gradient_plot.py" $regions
# conda run -n geo_env python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/prepare_sobel_reference_geotiff_and_plot.py" $regions $dayadj

# NOTE:
# Output here is a tif file
# ../figures_for_sobel_calculation/" + region + "_sobel_gradient.tif"

# Geolocate and project the SST  maps
#######################################################################
conda run -n geo_env python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/georeference_SST_map_using_GDAL.py" $regions

# NOTE:
# Output here is a tif file
# ../figures/SST_plot_"+region+"_proj.tif"

# echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
# echo "      Finished georeference_SST_map_using_GDAL.py"
# echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

# # add 10 minute grid to the map writing to the same file
# not working/ commented
# ####################################
# conda run -n pygmt python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/gridlines_for_SST_navigation_map.py" $regions

# Geolocate and project the Sobel maps
#######################################################################
conda run -n geo_env python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/georeference_Sobel_map_using_GDAL.py" $regions

# NOTE:
# Output here is a tif file
# ../figures/Sobel_plot_"+region+"_proj.tif"

# echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
# echo "      Finished georeference_Sobel_map_using_GDAL.py"
# echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

mv "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/SST_plot_"$regions"_proj.tif"  "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/"$regions"_overlay_proj.tif" 

mv "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/Sobel_plot_"$regions"_proj.tif"  "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/"$regions"_sobel_overlay_proj.tif"         

# NOTE: Since there is now no need to warp as second time (it is commented out above), just copy figures/"$regions"_overlay_proj.tif to figures/"$regions"_overlay_proj_print.tif
# so the montage will work

# SST 
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/"$regions"_overlay_proj.tif" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/"$regions"_overlay_proj_print.tif"
#Sobel gradient
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/"$regions"_sobel_overlay_proj.tif" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/"$regions"_sobel_overlay_proj_print.tif"        

conda run -n pygmt python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/plot_the_geotiff_with_pygmt.py" $regions

# output is
# ../figures_Avenza/"+region+"_overlay_proj_printable.tif"

conda run -n pygmt python "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/python/plot_the_sobel_geotiff_with_pygmt.py" $regions

# output is
# ../figures_Avenza/"+region+"_sobel_overlay_proj_printable.tif"

# Copy just the navigation geotiff to the figures Avenza directory (the printable map tiff file is already there)
#########################################################################
# SST
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/"$regions"_overlay_proj.tif" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_overlay_proj.tif"
# Sobel gradient
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/"$regions"_sobel_overlay_proj.tif" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_sobel_overlay_proj.tif"

# Convert the printable map tiff file to a png file for smaller downloads from website
#####################################################################
# Use '2>/dev/null' to remove warning about tif tags
cd /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza

mogrify -format png $regions"*_printable.tif" 2>/dev/null

# Move the printable map png file to  /recreational_fishing_high_resolution/figures_website_png_downloads for downloads and to serve as the source for each montage panel
#####################################################################
# SST

mv "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_overlay_proj_printable.png" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_website_png_downloads/"$regions"_overlay_proj_printable.png"
# Sobel gradient
mv "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_sobel_overlay_proj_printable.png" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_website_png_downloads/"$regions"_sobel_overlay_proj_printable.png"

# Create each montage panel from the printable map png file, incrementing on each pass of the loop
#####################################################  
# SST
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_website_png_downloads/"$regions"_overlay_proj_printable.png"  "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_d"$montage_map".png"

rm "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_website_png_downloads/"$regions"_overlay_proj_printable.png"

#Sobel gradient
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_website_png_downloads/"$regions"_sobel_overlay_proj_printable.png"  "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_e"$montage_map".png"

rm "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_website_png_downloads/"$regions"_sobel_overlay_proj_printable.png"

# copy the tif file needed for Leaflet
cp /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/$regions"_overlay_proj.tif" /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/leaflet_map_files/$regions"_sst_overlay_proj_day_"$dayadj".tif"

cp /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures/$regions"_sobel_overlay_proj.tif" /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/leaflet_map_files/$regions"_sst_fronts_overlay_proj_day_"$dayadj".tif"

((montage_map=montage_map+1))

# Delete the 4 montage panel files created in the loop
# best to keep them for trouble-shooting
# rm /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/*.png

done

# Create the SST and Sobel montages from the 4 montage panel files created in the loop
#########################################
cd /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage || exit

#SST
montage $regions"_d1.png" $regions"_d3.png" $regions"_d2.png" $regions"_d4.png" -geometry +2+2 $regions"_montage_SST.png"
#Sobel gradient
montage $regions"_e1.png" $regions"_e3.png" $regions"_e2.png" $regions"_e4.png" -geometry +2+2 $regions"_montage_Sobel_gradient.png"

# # animation code --------------------------------

# #############################################
# # Copy the SST and Sobel gradient daily files to the directory used to assemble the SST/ gradient animation 
# #############################################
# cp $regions"_d"* $regions"_e"* "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/sst_video/"

# ############################
# # create the SST/ gradient animation for this region
# # need to rename the files first so they will be in the correct order
# cd /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/sst_video/

# mv $regions"_d1.png" $regions"_1.png"
# mv $regions"_e1.png" $regions"_2.png"
# mv $regions"_d2.png" $regions"_3.png"
# mv $regions"_e2.png" $regions"_4.png"
# mv $regions"_d3.png" $regions"_5.png"
# mv $regions"_e3.png" $regions"_6.png"
# mv $regions"_d4.png" $regions"_7.png"
# mv $regions"_e4.png" $regions"_8.png"

# ############################
# # create the animation file
# # control the resolution using scale to reduce file size
# # use scale=-1:1000 to keep original aspect ratio
# /usr/bin/ffmpeg -framerate 0.4 -pattern_type glob -i '*.png' -y -vf scale=-1:1000 -vcodec libx264 -acodec aac $regions"_SST_animation.mp4"

#############################################
# Copy the SST montage, SST latest day, and the SST navigation file to the directory used to update the Ecwid store
#############################################

# # Copy the SST/ Sobel animation file
# cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/sst_video/"$regions"_SST_animation.mp4" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/latest_files_for_subscription/"$regions"_SST_animation.mp4"

# # be sure to cleanup the png files before re-running ffmpeg for the next region or the .avi file will use all of the remaining png files
# rm /data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/sst_video/*.png

# # end animation code --------------------------------

#############################################
# Reduce file sizes using image magick resize to 20%
# overwrite original file name
#############################################
# Resize the SST montage file
convert "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_montage_SST.png" -resize 50% "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_montage_SST.png"

# Resize the latest day of SST
convert "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_d4.png" -resize 50% "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_d4.png"

# Resize the Sobel gradient montage file
convert "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_montage_Sobel_gradient.png" -resize 50% "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_montage_Sobel_gradient.png"

# Resize the latest day of Sobel gradient
convert "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_e4.png" -resize 50% "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_e4.png"

# #############################################
# # Reduce navigation file sizes using image magick resize to 20%
# # DO NOT overwrite original file
# # It is needed to transfer the tag data from the original geotiff to the resized geotiff
# #############################################
# # Resize and add tags to the SST navigation file
# convert "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_overlay_proj.tif" -resize 50% "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_overlay_proj_resized.tif"

# exiftool -TagsFromFile "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_overlay_proj.tif" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_overlay_proj_resized.tif"

# # Resize and add tags to the Sobel gradient navigation file
# convert "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_sobel_overlay_proj.tif" -resize 50% "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_sobel_overlay_proj_resized.tif"

# exiftool -TagsFromFile "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_sobel_overlay_proj.tif" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_sobel_overlay_proj_resized.tif"


#############################################
# Copy the SST & Sobel montage to the directory used to update the Ecwid store 
#############################################

# Copy the SST montage file
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_montage_SST.png" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/latest_files_for_subscription/"$regions"_SST_last_4_days.png"

# Copy the Sobel gradient montage file
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_montage_Sobel_gradient.png" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/latest_files_for_subscription/"$regions"_SST_Breaks_last_4_days.png"

# Copy the latest day of SST
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_d4.png" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/latest_files_for_subscription/"$regions"_SST_to_print.png"

# Copy the latest day of Sobel gradient
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_animation/montage/"$regions"_e4.png" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/latest_files_for_subscription/"$regions"_SST_Breaks_to_print.png"

# Copy the SST navigation file
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_overlay_proj.tif" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/latest_files_for_subscription/"$regions"_SST_navigation_file.tif"

# Copy the Sobel gradient navigation file
cp "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/figures_Avenza/"$regions"_sobel_overlay_proj.tif" "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/latest_files_for_subscription/"$regions"_SST_Breaks_navigation_file.tif"

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "      Finished cleanup and copying files"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"