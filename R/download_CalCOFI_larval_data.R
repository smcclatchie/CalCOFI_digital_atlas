### download_CalCOFI_larval_data.R

# Just use one start and end date and adjust the labels in prepare ... reference geotiff and plot.py
# Day selection checked
# start_date <- paste(Sys.Date() - 1000, "T09:00:00Z", sep='')
# end_date <- paste(Sys.Date() - 2, "T09:00:00Z", sep='')

# print("-------------------------------------------------")
# print(start_date)
# print(end_date)
# print("-------------------------------------------------")
 
# download CalCOFI larval data as csv
# NOTE: Download all years, seasons, nets from ERDDAP 
# and subset later using R

# remove second header line of CSV file
# 2 do

##############################################
url <- "https://upwell.pfeg.noaa.gov/erddap/tabledap/erdCalCOFIlrvcnt.csv?cruise%2Ctow_type%2Clatitude%2Clongitude%2Cline%2Cstation%2Cscientific_name%2Ccalcofi_species_code%2Clarvae_10m2&time%3E=1951-01-09T18%3A16%3A00Z&time%3C=2023-01-25T18%3A45%3A00Z"

filename <- "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/CalCOFI_larval.csv"
download.file(url,filename)
    
print("############ end R script (data download) ##############")