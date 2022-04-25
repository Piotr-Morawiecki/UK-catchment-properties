## ---------------------------
##
## Script name: river_processing_2.R
##
## Purpose of script: Extracting length of all streams for UK catchments
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
##    - data/OS_VectorMap/...
##      (available here: https://www.ordnancesurvey.co.uk/business-government/products/vectormap-district)
##
##  Script extracts total length of all streams from OS_VectorMap database.
##  The results are exported to 'output/river_processing_2_summary.csv' file.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Directory to the input datasets (see section 'Notes' in the header above)

boundaries_file <- 'data//NRFA_catchment_boundaries//nrfa_public_ex20200925.shp'
vectorMap_directory <- 'data//OS_VectorMap//'

# File to which the results are saved

output_file <- "output//river_processing_2_summary.csv"

# Each tile of Ordnance Survey National Grid will be divided into
# n_subtiles x n_subtiles subtiles in order to reduce the total processing time

n_subtiles <- 5

## ---------------------------

# Installing required libraries

libraries_required <- c('rgdal', 'rgeos', 'sp', 'raster')
lapply(libraries_required, function(x) {
  if (!require(x, character.only=T)) {install.packages(x);require(x)}
})

# ------------------------
# | Function definitions |
# ------------------------

# Function try_readOGR attepts to read shapefile for specified region and dataset
# region - two letter code of Ordnance Survey National Grid
# dataset - 'SurfaceWater_Line' for small rivers and channels (lines)
#           'SurfaceWater_Area' for major rivers and lakes (polygons)
# crs - desired coordinate reference system
# Given shapefile is returned. If given shapefile is not available returns NA.

try_readOGR <- function(region, dataset, crs) {
  
  # Construct filename
  
  filename <- paste(vectorMap_directory, region, '//', region, '_', dataset,
                    '.shp', sep='')
  
  #If file does not exist return NA
  
  if (!file.exists(filename)) return(NA)
  
  # Otherwise read given shapefile
  
  data <- readOGR(dsn = filename, stringsAsFactors = F, verbose=FALSE)
  
  # Transform shapefile object to the desired coordinate system
  
  data <- spTransform(data, crs)
  
  # Return the shapefile object
  
  return(data)
}

# Function compute_river_length for given polyline x computes its length

compute_river_length = function(x) {
  
  # coordinates of line breakpoints are extracted
  
  x = x@Lines[[1]]@coords
  
  # If there are only two points length of the stream is computed from
  # Pythagorean Theorem
  
  if (nrow(x) == 2) return(sqrt(sum((x[1:(nrow(x)-1),] - x[2:nrow(x),]) ** 2)))
  
  # Otherwise the length of each section of stream is calculated from
  # Pythagorean Theorem and then all sections are summed.
  
  return(sum(sqrt(apply((x[1:(nrow(x)-1),] - x[2:nrow(x),]) ** 2, 1, sum))))
}

# Read file with catchment boundaries

catchments <- readOGR(dsn = boundaries_file, stringsAsFactors = F)

# Remove ring self-intersections, which may rarely occur in the input datasets

catchments <- gBuffer(catchments, byid=TRUE, width=0)

# Read list of all regions from vectorMap_directory

regions <- list.files(vectorMap_directory)

# Data types to be analyzed: surface water bodies both in form of lines (small
# stream and channels) and polygons (large rivers and lakes)

tile_types <- c('SurfaceWater_Line', 'SurfaceWater_Area')

# Initialise summary data frame in which total length of the streams for each
# catchment will be recorded

summary <- data.frame(id=catchments$STATION, length_B=0)

# For loop iterates over all regions, for each of them finding the total length
# of streams in given region overlapping with each catchment boundary.

for (i in 1:length(regions)) {
  region <- regions[i]
  
  # Report current progress
  
  print(paste('Region ', region, ' started (', i, '/', length(regions),').', sep=''), quote=FALSE)
  print('---------------------------', quote=FALSE)
  
  # For each data type different shapefile is read and processed
  
  for (type in tile_types) {
    
    # Import river data
    
    rivers <- try_readOGR(region, type, catchments@proj4string)
    
    # If such dataset does not exist proceed to the next iteration
    
    if (class(rivers) == "logical") next
    
    # Find extent of all rivers in a given tile
    
    tile_extent <- extent(rivers)
    xmin <- tile_extent@xmin
    xmax <- tile_extent@xmax
    ymin <- tile_extent@ymin
    ymax <- tile_extent@ymax
    
    # Find catchments which overlap with given regions
    
    catchments_in_tile <- crop(catchments, tile_extent)
    
    # If no catchments overlap proceed to the next iteration
    
    if (length(catchments_in_tile) == 0 ) next
    
    # Region is divided into n_subtiles x n_subtiles subtiles.
    # For each subtile (with coordinates x, y) rivers inside this subtile are
    # extracted and intersected with catchments. boundaries.
    # This reduces number of rivers that are intersected and therefore reduces
    # computational time.
    
    for (x in 1:n_subtiles){
      for (y in 1:n_subtiles) {
        
        # Report current progress
        
        print(paste('Subtile processed: ', n_subtiles * (x-1) + y, '/', n_subtiles**2, sep=''))
        
        # Compute subtile boundaries
        
        xmin2 <- xmin + (xmax - xmin) * (x-1) / n_subtiles
        ymin2 <- ymin + (ymax - ymin) * (y-1) / n_subtiles
        xmax2 <- xmin + (xmax - xmin) * x / n_subtiles
        ymax2 <- ymin + (ymax - ymin) * y / n_subtiles
        tile_extent <- extent(xmin2, xmax2, ymin2, ymax2)
        
        # Find catchments that may overlap with given subtile boundaries
        # (in order to reduce number of investigated catchments)
        
        catchments_subset <- crop(catchments_in_tile, tile_extent)
        
        # If there are no such catchments proceed to the next iteration
        
        if (length(catchments_subset) == 0 ) next
        
        # Find catchments that do overlap with given subtile boundaries
        
        rivers_in_tile <- crop(rivers, tile_extent)
        
        # If there are no such catchments proceed to the next iteration
        
        if (length(rivers_in_tile) == 0 ) next
        
        # Iterate over all catchments to find total length of all streams
        # in this subtile with given catchment boundaries
        
        for (j in 1:length(catchments_subset)) {
          catchment <- catchments_subset[j,]
          
          # Find row in summary data frame corresponding to the given catchment
          
          row_id <- which(summary$id==catchment$STATION)
          
          # Find and select rivers overlapping with the given catchment
          
          rivers_subset <- over(rivers_in_tile, catchment)
          rivers_subset <- rivers_in_tile[!is.na(rivers_subset$STATION),]
          
          # If there are no such rivers proceed to the next catchment
          
          if (length(rivers_subset)==0) next
          
          # Intersect extracted rivers with the catchment boundaries
          
          rivers_subset <- intersect(rivers_subset, catchment)
          
          # If the rivers are line objects calculate their length.
          # If they are polygons than additionally divide their perimeter by 2
          # (to take into account that river banks would be counted separately).
          
          if (type == 'SurfaceWater_Line') {
            lengths <- do.call(sum, lapply(rivers_subset@lines, function(x)
              compute_river_length(x@Lines[[1]]@coords)))
          } else {
            lengths <- do.call(rbind, lapply(rivers_subset@polygons, function(x)
              compute_river_length(x@Polygons[[1]]@coords))) / 2
          }
          
          # Add computed length to summary data frame
          
          summary$length_B[row_id] <- summary$length_B[row_id] + lengths
        }
      }
    }
  }
  
  # The summary is saved in csv file. If calculations are interrupted one can
  # continue script by reading this csv file and running for loop from the
  # last i value.
  
  write.csv(summary, 'output//river_processing_2_summary.csv', row.names = FALSE)
}

# Save the final version of the summary to the output_file.

write.csv(summary, 'output//river_processing_2_summary.csv', row.names = FALSE)
