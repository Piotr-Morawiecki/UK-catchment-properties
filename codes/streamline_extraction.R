## ---------------------------
##
## Script name: streamline_extraction.R
##
## Purpose of script: creating datasets of streamlines from DTM data for further 
##                    processing by streamline_processing.R script
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
##   data//GADM_UK_boundary//gadm36_GBR_0.shp
##  Script requires three input datasets:
##    - data/GADM_UK_boundary/gadm36_GBR_0.shp
##      (available here: https://biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_GBR_shp.zip)
##    - data/OS_VectorMap/...
##      (available here: https://www.ordnancesurvey.co.uk/business-government/products/vectormap-district)
##    - data/DIGIMAP_elevation/...
##      (available here: https://www.ordnancesurvey.co.uk/business-government/products/terrain-50)
##
##  Script for each tile of UK dtm raster finds streamline following the
##  steepest descent path. For each streamline its start point, end point,
##  length, elevation difference and type of object at the outlet is recorded.
##
##  The results are exported to 'output/streamflow_samples' directory.
##  Then they are used by streamline_processing.R script to estimate width of
##  UK catchments and mean value of elevation gradient along the hillslopes.
##
##  More information can be found in the appendix of the report.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Directory to the input datasets (see section 'Notes' in the header above)

UK_boudnary_file <- "data//GADM_UK_boundary//gadm36_GBR_0.shp"
vectorMap_directory <- 'data//OS_VectorMap//'
dtm_directory <-'data//DIGIMAP_elevation//'

# File to which the results are saved

output_directory <- "output//streamlines//"

## ---------------------------

# Installing required libraries

libraries_required <- c('rgdal', 'rgeos', 'raster')
lapply(libraries_required, function(x) {
  if (!require(x, character.only=T)) {install.packages(x);require(x)}
})

# ------------------------
# | Function definitions |
# ------------------------

# Function returnNA returns NA (used in tryCatch function later)

returnNA <- function(cond) return(NA)

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

# Function update_surface_water takes raster surface_water and replaces all
# tiles overlapping with boundary_shp with value id. This is used to construct
# surface_water raster containing information about the types of surface water
# localized in each tile.
# tile_extent is geographical extent from which objects of boundary_shp can
# overlap with the surface_water raster

update_surface_water <- function(surface_water, boundary_shp, tile_extent, id) {
  
  # Try cropping boundary_shp shape by provided tile_extent.
  # If the boundary_shp is empty then NA is returned.
  
  boundary <- tryCatch(crop(boundary_shp, tile_extent), error=returnNA)
  
  # Check if boundary is empty or NA (in case of empty boundary_shp argument)
  
  if (class(boundary)!="logical" & class(boundary)!="NULL") {
    
    # if boundary is made of polygons convert it into lines
    
    if (class(boundary) == "SpatialPolygonsDataFrame") {
      boundary <- as(boundary, 'SpatialLinesDataFrame')
    }
    
    # Replace surface_water values which overlap with the boundary with id value 
    
    surface_water[!is.na(mask(surface_water, boundary))] <- id
  }
  
  # Function retuns updated surface_water raster
  
  return(surface_water)
}

# Read shapefile with UK boundary

UK_boundary <- readOGR(dsn = UK_boudnary_file, stringsAsFactors = F)

# Find list of all dtm raster files in dtm_directory

dtm_files <- list.files(dtm_directory, pattern = "\\.asc$", recursive=TRUE)

# Read coordinate reference system from the first dtm raster file

dtm_crs <- raster(paste(dtm_directory, dtm_files[1], sep=''))@crs

# Transform UK boundary to dtm coordinate system

UK_boundary <- spTransform(UK_boundary, dtm_crs)

# Read list of all Ordnance Survey National Grid region codes from vectorMap_directory

regions <- list.files(vectorMap_directory)

# n_processed variable is used to keep track of number of tiles, which has been processed
# n_tiles is the total number of tiles to be processed

n_processed <- 0
n_tiles <- length(dtm_files)

# start_time is a reference time used for estimating expected computation time

start_time <- Sys.time()

# For loop iterates over all regions to find all streamlines for given region 

for (i in 1:length(regions)) {
  region <- regions[i]
  
  # Report the current computation progress
  
  print(paste('Region ', region, ' started (', i, '/', length(regions),').', sep=''), quote=FALSE)
  print('---------------------------', quote=FALSE)
  
  # Construct directory with dtm files for the given region
  
  directory <- paste(dtm_directory, tolower(region), sep='')
  
  # Find all dtm raster files in this directory
  
  dtm_files <- list.files(directory, pattern = "\\.asc$")
  
  # If no dtm files were found proceed to the next region 
  
  if (length(dtm_files) == 0) next
  
  # Construct data frame for storing streamline information
  # x, y                - starting point coordinates
  # x_outlet, y_outlet  - outlet coordinates
  # outlet_type         - id of the outlet type; possible values are:
  #                         1 - tidal boundary,
  #                         2 - large channel (polygon), 
  #                         3 - small channel (line),
  #                         4 - tile boundary
  #                         5 - local elevation minimum (without surface water)
  # distance            - length of the streamline
  # dtm_point           - evelation of the starting point
  # dtm_outlet          - elevation fo the outlet
  
  samples <- data.frame(x=c(), y=c(), x_outlet=c(), y_outlet=c(), outlet_type=c(),
                        distance=c(), dtm_point=c(), dtm_outlet=c())
  
  # For given region read surface water objects, icluding:
  # - small streams and channels stored as lines (rivers_line)
  # - rivers and lakes stored as polygons (rivers_area)
  # - tidal boundaries (e.g. coastline, bays) stored as lines (tidal_boudary)
  
  rivers_line <- try_readOGR(region, 'SurfaceWater_Line', dtm_crs)
  rivers_area <- try_readOGR(region, 'SurfaceWater_Area', dtm_crs)
  tidal_boundary <- try_readOGR(region, 'TidalBoundary', dtm_crs)
  
  # For each dtm file in a given region streamlines starting form and ending
  # at a given raster (or its boundary) are recorded.
  
  for (j in 1:length(dtm_files)) {
    
    # Report current progress
    
    print(paste('Region', i, '/', length(regions)), quote=FALSE)
    print(paste('DTM', j, '/', length(dtm_files)), quote=FALSE)
    
    # Read given dtm file
    
    dtm <- raster(paste(directory, dtm_files[j], sep='//'))
    
    # Extract geographical extent of the raster
    
    tile_extent <- extent(dtm)
    
    # find terrain gradient from dtm (it points to 1 out of 8 neighboring tiles
    # located at the direction of the steepest descent). Find more information
    # in terrain function documentation.
    
    gradient <- terrain(dtm, opt='flowdir')
    
    # The gradient data frame is used to find dx and dy at each point,
    # where dx shows x-component of the steepest descent direction (-1, 0 or 1),
    # and dy shows its y-component (-1, 0, 1)
    
    dx <- gradient
    dx[,] <- 0
    dy <- dx
    dx[gradient%in%c(128, 1, 2)] <- 1
    dx[gradient%in%c(8, 16, 32)] <- -1
    dy[gradient%in%c(2, 4, 8)] <- -1
    dy[gradient%in%c(32, 64, 128)] <- 1
    
    # dx and dy are converted to matrices
    
    dx_array <- as.matrix(dx)
    dy_array <- as.matrix(dy)
    
    # Find locations of local elevation minimums (points with elevation lower
    # than elevation of their 8 neighbours). See focal function documentation
    # for more details. Output is converted to matrix form.
    
    local_minimum <-  focal(dtm, w=matrix(1,nrow=3,ncol=3),
                            fun=function(x, ...) min(x) == x[5],
                            pad=TRUE, padValue=FALSE, na.rm=TRUE)
    local_minimum <- as.matrix(local_minimum)
    
    # surface_water raster is created. It has the same size as dtm raster,
    # but stores information about the type of the terrain at given point.
    # It uses the same values as outlet_type parameter in samples data frame
    # defined before - see their descriptions in the comment above.
    # Initially it has values 0.
    
    surface_water <- dtm
    values(surface_water) <- 0
    
    # Values of tiles overlapping with rivers_line are assigned to 3
    surface_water <- update_surface_water(surface_water, rivers_line, tile_extent, 3)
    
    # Values of tiles overlapping with rivers_area are assigned to 2
    surface_water <- update_surface_water(surface_water, rivers_area, tile_extent, 2)
    
    # Values of tiles overlapping with tidal_boundary are assigned to 1
    surface_water <- update_surface_water(surface_water, tidal_boundary, tile_extent, 1)
    
    # Values of tiles laying outside the UK_boundary are assigned to 1
    surface_water[is.na(mask(surface_water, UK_boundary))] <- 1
    
    # Plotting option allowing to plot surface_water raster and surface water
    # shapefiles (see example in the report).
    plotting <- false
    if (plotting) {
      par(mfrow=c(1,2))
      plot(dtm, main='Digital terrain model with rivers')
      plot(crop(rivers_line, tile_extent),add=TRUE)
      plot(crop(rivers_area, tile_extent),add=TRUE)
      plot(surface_water, main='Surface water raster')
    }
    
    # Convert surface_water to matrix
    
    surface_water <- as.matrix(surface_water)
    
    # Assign value 4 at dtm raster boundaries
    
    surface_water[1,] <- 4
    surface_water[,1] <- 4
    surface_water[nrow(surface_water),] <- 4
    surface_water[,ncol(surface_water)] <- 4
    
    # Assign value 5 at local elevation minima without surface water
    
    surface_water[local_minimum==1 & surface_water==0] <- 5
    
    # Data frame data is initiated to stores information about streamlines.
    # Points with default value (0) assigned to them in the surface_water
    # matrix are starting points of constructed streamlines.
    # Initially the endpoint is identical to start point, however then
    # iterative it is updated following steepest descent direction unless
    # reaching surface water, local minimum or catchment boundary. 
    
    data <- which(surface_water==0, arr.ind = T)
    data <- data.frame(x=data[,1], y=data[,2], x_outlet=data[,1], y_outlet=data[,2],
                       outlet_type=numeric(nrow(data)), distance=numeric(nrow(data)))
    
    # data_final data frame is initiatied. It will store the streamlines, which
    # already reached an outlet and don't need to be further processed.
    
    data_final <- data.frame(x=c(), y=c(), x_outlet=c(), y_next=c(), outlet_type=c(), distance=c())
    
    # Streamlines are updated unless all are transferred to data_final.
    
    while (nrow(data) > 0) {
      
      # Find value of dx and dy for all streamline ends
      
      dx <- as.numeric(Map(function(i, j) dx_array[i, j], data$x_outlet, data$y_outlet))
      dy <- as.numeric(Map(function(i, j) dy_array[i, j], data$x_outlet, data$y_outlet))
      
      # Update outlest coordinates following steepest descent (dx, dy)
      
      data$x_outlet <- data$x_outlet - dy
      data$y_outlet <- data$y_outlet + dx
      
      # Update distance by adding length of displacement vector
      
      data$distance <- data$distance + sqrt(dx^2 + dy^2)
      
      # Update outlet_type value by checking value of surface_water at new
      # outlet coordinates
      
      data$outlet_type <- as.numeric(Map(function(i, j) surface_water[i, j], data$x_outlet, data$y_outlet))
      
      # Add streamlines, which reached non-default outlet_type
      
      data_final <- rbind(data_final, data[data$outlet_type>0,])
      
      # Only streamlines which did not reach non-default outlet_type are left
      
      data <- data[data$outlet_type==0,]
    }
    
    # dtm is conveted to matrix form and used to find elevation value at
    # start and end point of each streamline
    
    dtm <- as.matrix(dtm)
    data_final$dtm_point <- as.numeric(Map(function(i, j) dtm[i, j], data_final$x, data_final$y))
    data_final$dtm_outlet <- as.numeric(Map(function(i, j) dtm[i, j], data_final$x_outlet, data_final$y_outlet))
    
    # The distance is multiplied by raster resolution (50 meters)
    
    data_final$distance <- data_final$distance * 50
    
    # x and y coordinates of start and end points are converted to values
    # corresponding to easting and northing as defined by dtm crs
    
    data_final$x <- tile_extent@xmin + 50 * data_final$x + 25
    data_final$y <- tile_extent@ymin + 50 * data_final$y + 25
    data_final$x_outlet <- tile_extent@xmin + 50 * data_final$x_outlet + 25
    data_final$y_outlet <- tile_extent@ymin + 50 * data_final$y_outlet + 25
    
    # Append streamlines to samples data frame
    
    samples <- rbind(samples, data_final)
    
    # Write the current samples data frame to .csv file
    
    output_file <- paste(output_directory, 'samples_', region, '.csv', sep='')
    write.csv(samples, output_file, row.names = FALSE)
    
    # update number of tiles processed
    
    n_processed <- n_processed + 1
    
    # check time taken by the computations and express in minutes
    time_taken <- Sys.time() - start_time
    time_taken <- as.numeric(time_taken, units = "mins")
    
    # estimate the remaining time by assuming that remaining tiles will on
    # average takes the same amount of time per tile
    
    time_remaining <- time_taken / n_processed * (n_tiles - n_processed)
    
    # Report the progress, time passed and time remaining
    
    print(paste('Progress:', n_processed, '/', n_tiles), quote=FALSE)
    print(paste('Time passed:', time_taken, 'minutes.'), quote=FALSE)
    print(paste('Time remaining:', time_remaining, 'minutes.'), quote=FALSE)
    print('---------------------------', quote=FALSE)
  }
}