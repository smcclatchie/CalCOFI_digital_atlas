#!/bin/bash
# download_regional_datasets_SST.sh

# ##########################################################
# USAGE (from ./shell):
# sh ./download_regional_datasets_SST.sh "Southwest_US"
# sh ./download_regional_datasets_SST.sh "Baja_Mexico_only"

# where the argument is the broader region containing the sub-region maps

# pass in the arguments
superregion=$1
# day_start=$2

# Delete previous data
##################################################

Rscript "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/R/download_GHRSST_MUR-JPL-L4-GLOB-v4.1_daily_Regional.R" --args $superregion
# if [ $? -ne 0 ]; then
#     echo "ERROR: Download failed for $superregion after all retries"
#     echo "ERDDAP SST download failed for region '$superregion' on $(hostname) at $(date). Check the cron log for details." \
#         | mail -s "Fishing Maps: Download Failed - $superregion" smcclatchie@fishingmaps.info
#     exit 1
# fi
echo "Superregion = "$1

