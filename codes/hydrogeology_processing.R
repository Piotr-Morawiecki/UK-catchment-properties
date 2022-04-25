## ---------------------------
##
## Script name: hydrogeology_processing.R
##
## Purpose of script: Extracting hydrogeologic properties for UK catchments
##
## Author: Piotr Morawiecki
##
## Date Created: 2021-08-19
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
##    - data/DIGIMAP_hydrogeology/HydrogeologyUK_IoM_v5.shp
##      (available here: https://www.bgs.ac.uk/datasets/hydrogeology-625k/)
##
##  Script extracts area of each catchment which is classified to
##  different aquifier character (CHARACTER) and flow mechanisms (FLOW_MECHA)
##  classes according to the hydrogeology-625k dataset. The results are
##  exported to 'output/hydrogeology_summary.csv' file.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Directory to the input datasets (see section 'Notes' in the header above)

boundaries_file <- 'data//NRFA_catchment_boundaries//nrfa_public_ex20200925.shp'
hydrogeology_file <- 'data//DIGIMAP_hydrogeology//HydrogeologyUK_IoM_v5.shp'

# File to which the results are saved

output_file <- "output//hydrogeology_summary.csv"

# Entire UK will be divided into nx x ny tiles. For each tile intersection of
# hydrogeology_data and catchment_boundaries polygons will be found
# independently for faster processing.

nx <- 20
ny <- 30
n_parts <- nx * ny

## ---------------------------

# installing required libraries

libraries_required <- c('rgdal', 'rgeos', 'sf', 'raster', 'reshape2', 'plyr')
lapply(libraries_required, function(x) {
  if (!require(x, character.only=T)) {install.packages(x);require(x)}
})

# Read input datasets

catchment_boundaries <- readOGR(dsn = boundaries_file, stringsAsFactors = F)
hydrogeology_data <- readOGR(dsn = hydrogeology_file, stringsAsFactors = F)

# Remove ring self-intersections, which may rarely occur in the input datasets

catchment_boundaries <- gBuffer(catchment_boundaries, byid=TRUE, width=0)
hydrogeology_data <- gBuffer(hydrogeology_data, byid=TRUE, width=0)

# Not available flow mechanisms are replaced with an 'Unknown' label

hydrogeology_data@data$FLOW_MECHA[is.na(hydrogeology_data@data$FLOW_MECHA)] <- 'Unknown'

# All unique values of aquifer character and flow mechanism are extracted

character_types <- unique(hydrogeology_data@data$CHARACTER)
mechanism_types <- unique(hydrogeology_data@data$FLOW_MECHA)

# Extent of each polygon forming hydrogeology_data and catchment_boundaries
# datasets are found. They will be used to narrow down the number of potentially
# overlapping polygons.

hydrogeology_extent <- sapply(1:nrow(hydrogeology_data), function(x) extent(hydrogeology_data[x,]))
hydrogeology_xmin <- unlist(lapply(hydrogeology_extent, function(x) x@xmin))
hydrogeology_xmax <- unlist(lapply(hydrogeology_extent, function(x) x@xmax))
hydrogeology_ymin <- unlist(lapply(hydrogeology_extent, function(x) x@ymin))
hydrogeology_ymax <- unlist(lapply(hydrogeology_extent, function(x) x@ymax))

catchment_extent <- sapply(1:nrow(catchment_boundaries), function(x) extent(catchment_boundaries[x,]))
catchment_xmin <- unlist(lapply(catchment_extent, function(x) x@xmin))
catchment_xmax <- unlist(lapply(catchment_extent, function(x) x@xmax))
catchment_ymin <- unlist(lapply(catchment_extent, function(x) x@ymin))
catchment_ymax <- unlist(lapply(catchment_extent, function(x) x@ymax))

# Find extent of all UK catchments

UK_extent <- extent(catchment_boundaries)

# Initially empty table is initiated with columns corresponding to all aquifer
# character and flow mechanism types appearing in the hydrogeology dataset.

colnames <- c("STATION", character_types, mechanism_types)
summary <- data.frame(matrix(ncol = length(colnames), nrow = 0))
colnames(summary) <- colnames

# x and y denote coordinates of the tile currently being processed
for (y in 1:ny) {
  for (x in 1:nx) {
    
    # Print information about the progress
    
    part <- (y - 1) * nx + x
    print(paste('Processing part ', part, "/", n_parts, sep=''))
    
    # Extent of the tile is found
    # (UK_extent@xmax - UK_extent@xmin) / nx is width of each tile
    # (UK_extent@ymax - UK_extent@ymin) / nx is height of each tile
    
    xmin <- UK_extent@xmin + (UK_extent@xmax - UK_extent@xmin) * (x-1) / nx
    ymin <- UK_extent@ymin + (UK_extent@ymax - UK_extent@ymin) * (y-1) / ny
    xmax <- UK_extent@xmin + (UK_extent@xmax - UK_extent@xmin) * x / nx
    ymax <- UK_extent@ymin + (UK_extent@ymax - UK_extent@ymin) * y / ny
    tile_extent <- extent(xmin, xmax, ymin, ymax)
    
    # catchment_subset stores catchment_boundaries polygons, which may
    # overlap with the given tile.
    
    catchment_subset <- catchment_boundaries[
      catchment_xmin < tile_extent@xmax & catchment_xmax > tile_extent@xmin &
        catchment_ymin < tile_extent@ymax & catchment_ymax > tile_extent@ymin, ]
    
    # If no catchments appear in this tile we proceed to the next tile
    
    if (length(catchment_subset) == 0) next
    
    # Selected catchments are cropped to the tile boundaries
    
    catchment_subset <- crop(catchment_subset, tile_extent)
    
    # If no catchments are left after cropping we proceed to the next tile
    
    if (length(catchment_subset) == 0) next
    
    # tile extent is further narrowed down to the extent of selected catchments
    
    tile_extent <- extent(catchment_subset)
    
    # hydrogeology_data_subset stores hydrology_data polygons at least partially
    # overlapping the given tile.
    
    hydrogeology_data_subset <- hydrogeology_data[
                                    hydrogeology_xmin < tile_extent@xmax &
                                    hydrogeology_xmax > tile_extent@xmin &
                                    hydrogeology_ymin < tile_extent@ymax &
                                    hydrogeology_ymax > tile_extent@ymin, ]
    
    # If no polygons appear in this tile we proceed to the next tile
    
    if (length(hydrogeology_data_subset) == 0 ) next
    
    # Intersection between the hydrogeology and catchment polygons is found 
    
    intersections <- intersect(hydrogeology_data_subset, catchment_subset)
    
    # area of each intersection is calculated
    
    intersections$area <- area(intersections)
    
    # data frames recording the total area of intersections for each catchment
    # (indentified with STATION value) and flow machanism and aquifer character
    # are created
    
    mechanism_agg <- aggregate(area ~ STATION + FLOW_MECHA, data=intersections, FUN=sum)
    character_agg <- aggregate(area ~ STATION + CHARACTER, data=intersections, FUN=sum)
    
    mechanism_agg <- dcast(mechanism_agg, STATION ~ FLOW_MECHA, value.var='area')
    character_agg <- dcast(character_agg, STATION ~ CHARACTER, value.var='area')
    
    # obtained data frames are added to the summary table
    
    summary <- rbind.fill(summary, merge(mechanism_agg, character_agg))
    summary[is.na(summary)] <- 0
  }
}

# all records in the summary are aggegated in order to find the total area of
# each catchment classified to each flow mechanism and aquifer character

summary <- aggregate(. ~ STATION, summary, sum)

# STATION identifier's name is changed to 'id'

colnames(summary)[1] <- 'id'

# the summary is recorded in the output_file

write.csv(summary, output_file, row.names=FALSE)
