### download GHRSST, MUR-JPL-L4-GLOB-v4.1 daily

# library(httr2)

args <- commandArgs(trailingOnly = TRUE)

input <- strsplit(args, split=" ")
region <- input[[2]]

setwd("/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/R")

# Just use one start and end date and adjust the labels in prepare ... reference geotiff and plot.py
# Day selection checked
start_date <- "2025-06-16T00:00:00Z"
end_date <- "2026-05-16T00:00:00Z"

print("-------------------------------------------------")
print(start_date)
print(end_date)
print("-------------------------------------------------")

# Download one ERDDAP URL to a file.
# Aborts with a clear message if the response is not NetCDF (e.g. an ERDDAP HTML error page).
# Retries up to 3 times with a 30-second wait on HTTP 429/503 (ERDDAP busy responses).
# download_file <- function(url, filename, region_name) {
#   resp <- tryCatch(
#     request(url) |>
#       req_timeout(600) |>
#       req_retry(max_tries = 3, backoff = ~ 30) |>
#       req_perform(),
#     error = function(e) stop(paste("Download failed for", region_name, ":", conditionMessage(e)))
#   )

#   content_type <- resp_content_type(resp)
#   if (!grepl("netcdf|octet-stream", content_type, ignore.case = TRUE)) {
#     stop(paste(
#       "Download for", region_name, "returned unexpected content-type:", content_type,
#       "- likely an ERDDAP error page, not a NetCDF file."
#     ))
#   }

#   writeBin(resp_body_raw(resp), filename)
#   cat("Downloaded", region_name, "->", basename(filename),
#       "| content-type:", content_type, "\n")
# }

# download MUR GHRSST data (Use monthly, not daily)
##############################################

# Southwest US
# Use Celcius dataset and do Fahrenheit conversion manually because ERDDAP jplMURSST41F database updates later than jplMURSST41
if ("Southwest_US" %in% region) {
    # url <-https://upwell.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.nc?sst%5B(2026-01-16):1:(2026-05-16T00:00:00Z)%5D%5B(28):1:(49)%5D%5B(-128):1:(-115)%5D,mask%5B(2026-01-16):1:(2026-05-16T00:00:00Z)%5D%5B(28):1:(49)%5D%5B(-128):1:(-115)%5D
    url <- paste("https://upwell.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.nc?sst%5B(",end_date,"):1:(",start_date,")%5D%5B(28.0):1:(49.0)%5D%5B(-128.0):1:(-115.0)%5D", sep="")
    print(url)
    filename <- "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/latest_MUR_SST_Southwest_US.nc"
    download.file(url, filename)
    }
# end Southwest US

# Baja_Mexico_only
# Use Celcius dataset and do Fahrenheit conversion manually because ERDDAP jplMURSST41F database updates later than jplMURSST41
if ("Baja_Mexico_only" %in% region) {
    url <- paste("https://upwell.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.nc?sst%5B(",end_date,"):1:(",start_date,")%5D%5B(18.0):1:(31.0)%5D%5B(-118.0):1:(-100.0)%5D", sep="")
    print(url)
    filename <- "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/latest_MUR_SST_Baja_Mexico_only.nc"
    download.file(url, filename)
    }
# end Baja_Mexico_only

# Northwest US
if ("Northwest_US" %in% region) {
    url <- paste("https://upwell.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.nc?sst%5B(",end_date,"):1:(",start_date,")%5D%5B(39.0):1:(49.0)%5D%5B(-128.0):1:(-122.0)%5D", sep="")
    print(url)
    filename <- "/data_7TB/mnt/data/dynamic_data/projects/projects2026/CalCOFI_digital_atlas/data/regional_downloads/latest_MUR_SST_Northwest_US.nc"
    download.file(url, filename)
    }
# end Northwest US


# download MUR GHRSST anomaly data
##############################################


print("############ end R script (data download) ##############")
