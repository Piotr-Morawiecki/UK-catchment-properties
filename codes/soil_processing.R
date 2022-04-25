## ---------------------------
##
## Script name: soil_processing.R
##
## Purpose of script: Extracting MvG model parameters for UK catchments
##
## Author: Piotr Morawiecki
##
## Date Created: 2021-08-20
##
## Copyright (c) Piotr Morawiecki, 2021
## Email: pwm27@bath.ac.uk
##
## ---------------------------
##
## Notes:
##   
##  Script requires two input datasets:
##    - data/NRFA_catchment_boundaries/nrfa_public_ex20200925.shp
##      (more information: https://nrfa.ceh.ac.uk/content/catchment-boundary-and-areas)
##    - data/ESDAC_soil_conductivity/...
##      (can be requested here: https://esdac.jrc.ec.europa.eu/content/3d-soil-hydraulic-database-europe-1-km-and-250-m-resolution)
##
##  Script extractes mean value of all Mualem-Van Genuchten model
##  (K_s, alpha, theta_s, theta_r, n) on different depths for UK catchments.
##  The results are exported to 'output/soil_processing_summary.csv' file.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Directory to the input datasets (see section 'Notes' in the header above)

boundaries_file <- 'data//NRFA_catchment_boundaries//nrfa_public_ex20200925.shp'
ESDAC_dataset_directory <- 'data//ESDAC_soil_conductivity//'

# File to which the results are saved

output_file <- "output//soil_summary.csv"

## ---------------------------

# Installing required libraries

libraries_required <- c('rgdal', 'rgeos', 'sp', 'raster')
lapply(libraries_required, function(x) {
  if (!require(x, character.only=T)) {install.packages(x);require(x)}
})

# Read datasets with catchment boundaries

catchments <- readOGR(dsn = boundaries_file, stringsAsFactors = F)

# Remove ring self-intersections, which may rarely occur in the input datasets

catchments <- gBuffer(catchments, byid=TRUE, width=0)

# Transform catchment boundaries to the coordinate system used in ESDAC dataset 

catchments <- spTransform(catchments, CRS('+proj=longlat +datum=WGS84 +no_defs'))

# Extract UK catchment geographical extent

UK_extent <- extent(catchments)

# List of parameters to be extracted. Available parameters are:
# 'FC', 'HCC_alp', 'HCC_K0', 'HCC_L', 'HCC_m', 'HCC_n', 'HCC_thr', 'HCC_ths',
# 'KS', 'MRC_alp', 'MRC_m', 'MRC_n', 'MRC_mthr', 'MRC_ths', 'THS', 'WP')
# Description of paramters is available in the ESDAC dataset documentation.
# Here we listMualem-Van Genuchten model (K_s, alpha, theta_s, theta_r, n).

data_types <- c('KS', 'HCC_thr', 'HCC_ths', 'HCC_alp', 'HCC_n')

# Depths from which data are extracted. All possible depths are listed.

depths <- c(0, 5, 15, 30, 60, 100, 200)

# Conversion of units from the ESDAC dataset to SI units

unit_converstion <- c(1/(100*100*24*3600), 1/10000, 1/10000, 100/10000, 1/10000)

# Initialize dataframe to record extracted values for each catchment

summary <- data.frame(id=c(), depth=c(), KS=c(), HCC_thr=c(),
                      HCC_ths=c(), HCC_alp=c(), HCC_n=c())

# start file count and a clock to track computation time

files_done <- 0
files_total <- length(depths) * length(data_types)
start_time <- Sys.time()

# There is a different dataset for each depth and data_type.
# Outer loop iterates over all depths, and inner over data types.

for (i in 1:length(depths)) {
  depth <- depths[i]
  
  # A new row is initiated, which will then be appended to the main summary.
  
  summary_part <- data.frame(id=catchments$STATION, depth=depth)
  
  # For each depth iterate over all data types.
  
  for (j in 1:length(data_types)) {
    type <- data_types[j]
    
    # Log the current progess
    
    print(paste('Processed dataset: ', type, ' (depth: ', depth, ')', sep=''), quote=FALSE)
    
    # Download raster representing values of given data type at given depth
    
    raster_name <- paste(ESDAC_dataset_directory, type, '_sl', i, '.tif', sep='')
    imported_raster <- raster(raster_name)
    
    # Crop map of the whole Europe from file to area of UK catchments 
    
    imported_raster <- crop(imported_raster, UK_extent)
    
    # Extract mean value of imported_raster for each UK catchment
    
    mean_value <- unlist(lapply(extract(imported_raster, catchments), mean))
    
    # Convert obtained mean_value to SI units
    
    mean_value <- mean_value * unit_converstion[j] 
    
    # Add value to the new row of the summary
    
    summary_part[,type] <- mean_value
    
    # Estimate remaining computation time and print the computation progress
    
    files_done <- files_done + 1
    time_passed <- as.numeric(difftime(Sys.time(), start_time), units='hours')
    remaining_time <- round(time_passed / files_done * (files_total - files_done), digits=2)
    print(paste('Remaining time:', remaining_time,'hours'), quote=FALSE)
  }
  
  # Add new row to the main summary data frame
  
  summary <- rbind(summary, summary_part)
  
  # Write the current version of summary to the output_file.
  # If the script is interrupted one can read this csv file and continue
  # running for loop starting from the last i value of the inner loop.
  
  write.csv(summary, output_file, row.names=FALSE)
}

# Sort rows of the summary by catchment id

summary <- summary[order(summary$id),]

# Write final version of the summary to output_file

#write.csv(summary, output_file, row.names=FALSE)
