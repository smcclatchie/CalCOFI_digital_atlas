#!/bin/bash
# download_CalCOFI_larval_data.sh

##########################################################
# USAGE (from ./shell):
# sh ./download_CalCOFI_larval_data.sh 

# where the argument is the broader region containing the sub-region maps

# pass in the arguments
superregion=$1
# day_start=$2

# Delete previous data
##################################################
# rm `/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/`CalCOFI_larval_$superregion*.nc

Rscript "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/R/download_CalCOFI_larval_data.R"

# NOTE: NZ download takes about 80 sec, so allow for that in pause
# then use a loop with a test to determine if the file needs re-downloading. Allow 90 sec for each repeat download and increase the pause by 5 seconds for each iteration of the loop. Roy recommends 5-10 repeats when accessing data from a script.

test -f "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/CalCOFI_larval.nc"
if [ $? -ne 0 ] 
then
    for pause in 90 95 
    do
        echo "Missing data file. Retrying in "$pause" seconds ..."
        # wait the specified number of seconds to retry download
        sleep $pause
        Rscript "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/R/download_CalCOFI_larval_data.R"
        # test  if download was successful (exit status=0)
        test -f "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/CalCOFI_larval.nc"
        if [ $? -eq 0 ] 
        then
            break
        fi
    done
fi

