## ---------------------------
##
## Script name: aquifer_processing.R
##
## Purpose of script: estimating mean depth of aquifier for UK catchments
##
## Author: Piotr Morawiecki
##
## Date Created: 2021-08-21
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
##    - data/BGS_Boreholes/BGS_UK3D_v2015_3D_Cross_Sections.shp
##      (available here: https://www.bgs.ac.uk/datasets/uk3d/)
##
##  Script uses 3D profiles obtained from borehole measurements to
##  estimate mean depth of different aquifer layers. The method is described
##  in detail in the report.
##  Results are exported to 'output//aquifer_summary.csv' file.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Directory to the input datasets (see section 'Notes' in the header above)

boundaries_file <- 'data//NRFA_catchment_boundaries//nrfa_public_ex20200925.shp'
boreholes_file <- "data/BGS_Boreholes/BGS_UK3D_v2015_3D_Cross_Sections.shp"

# File to which the results are saved

output_file <- "output//aquifer_summary.csv"

## ---------------------------

# Installing required libraries

libraries_required <- c('rgdal', 'sf', 'raster', 'sp')
lapply(libraries_required, function(x) {
  if (!require(x, character.only=T)) {install.packages(x);require(x)}
})

# Function triangle_area calculates area of a triangle with coordinates given
# by 3x2 matrix x. The formula is derived from the formula for the cross product
# (see https://math.stackexchange.com/questions/128991/how-to-calculate-the-area-of-a-3d-triangle)

triangle_area <- function(x) {
  ab <- x[1,] - x[2,]
  ac <- x[1,] - x[3,]
  area <- 0.5 * sqrt((ab[2]*ac[3]-ab[3]*ac[2])^2 + (ab[3]*ac[1]-ab[1]*ac[3])^2 + (ab[1]*ac[2]-ab[1]*ac[2])^2)
  return(area)
}

# Function total_area computes total area of polygon x composed of triangular
# polygons. x should be a list containing a single list of triangular polygons
# (this is format used in the input UK3D dataset)

total_area <- function(x) sum(unlist(lapply(x, function(x) triangle_area(x[[1]]))))

# Import input UK3D dataset

profiles <- st_read(boreholes_file)

# Compute area of each polygons

profiles[,"area"] <- unlist(lapply(profiles$geometry, function(x) total_area(x)))

# ---------------------------------------
# | Computing mean depth over entire UK |
# ---------------------------------------

# Find total area occupied by each aquifer type. Categories are:
# AQUIFER_D - describes given aquifer zone as either 'Principal', 'SecondaryA',
#             'SecondaryB', 'SecondaryU', 'Unproductive' or 'NULL' (unknown);
#             More information about each type are provided in the report
#             and the UK3D documentation
# HYDRO_D   - described productivity of given aquifer and ranges from
#             highly productive aquifers with intergranular flow
#             to aquifers with esentially no groundwater

area_summary_A <- aggregate(area ~ AQUIFER_D, profiles, sum)
area_summary_AH <- aggregate(area ~ AQUIFER_D + HYDRO_D, profiles, sum)

# According to UK3D documentation profiles have total length of around 20,000 km
# (here we convert it into meters)

cross_section_length <- 20e3 * 1e3

# Mean depth is estimated by dividing total polygon area corresponding to
# given category by the total profile length

area_summary_A$area <- area_summary_A$area / cross_section_length
area_summary_AH$area <- area_summary_AH$area / cross_section_length

# Colname is changed from area to mean depth

colnames(area_summary_A)[2] <- "mean depth"
colnames(area_summary_AH)[3] <- "mean depth"

# Obtained depth are displayed

area_summary_A
area_summary_AH

# ------------------------------------------------------
# | Computing mean depth over individual UK catchments |
# ------------------------------------------------------

# Define width of polygons used to estimate profile length
# (more details can be found in the report)

d <- 15

# Read file with catchment boundaries

catchments <- readOGR(dsn = boundaries_file, stringsAsFactors = F)

# Remove ring self-intersections, which may rarely occur in the input datasets

catchments <- gBuffer(catchments, byid=TRUE, width=0)

# Function convert to line takes a 3D (potentially tilted) triangular polygon x
# divides it in to parts and projects it into xy-plane according to the method
# described in the report.
# At the output list of the following objects is obtained:
# line1, line2  - spatial line describing the projection of the first
#                 and second triangle on the xy-plane respectively
# base1, base2  - length of above projections
# dz            - height of the triangle at the midpoint
#                 (the same for both triangles)
# x_mid, y_mid  - coordinates of the midpoint
# polygon       - polygon of width 2d with centerline at the base projection
#                 used to calculate the profile length

convert_to_line <- function(x, d) {
  
  # Polygons coordinates are obtained
  
  x <- x[[1]]
  
  # Find extent of the polygon in x and y direction
  
  values_range <-apply(x[,1:2], 2, function(x) diff(range(x)))
  
  # If the x and y coordinates of all points on the polygon are identical
  # return NA
  
  if (sum(values_range) == 0) return(NA)
  
  # Check whether extent is higher in x (dir=0) or y-direction (dir=2)
  # It is found to guarantee that given coordinate change along the polygon
  
  dir <- which.max(values_range)
  
  # Find the location of the start and end point along this dimension in x table
  
  start <- which.min(x[1:3,dir])[1]
  end <- which.max(x[1:3,dir])[1]
  
  # Find the location of the remaining midpoint
  
  mid <- setdiff(1:3, c(start, end))
  
  # Extract midpoint coordinates
  
  x_mid <- x[mid,1]
  y_mid <- x[mid,2]
  
  # Find the position along z-axis of the point laying on the line joining
  # start and end point, and located at the same xy position as the midpoint
  
  z_mid <- (x[end,3] * (x[mid,dir] - x[start,dir]) +
    x[start,3] * (x[end,dir] - x[mid,dir])) / (x[end,dir] - x[start,dir])
  
  # Find vertical extent of the triangle at the mid point
  
  dz <- abs(z_mid - x[mid, 3])
  
  # Construct Line objectsdescribing the projection of the first and second
  # triangle on the xy-plane respectively
  
  line1 <- Line(x[c(start, mid), 1:2])
  line2 <- Line(x[c(end, mid), 1:2])
  
  # Find length of these projections
  
  base1 <- sqrt(sum(diff(x[c(start, mid), 1:2])^2))
  base2 <- sqrt(sum(diff(x[c(end, mid), 1:2])^2))
  
  # Find coordinates of the start and end points
  
  r1 <- x[start, 1:2]
  r2 <- x[end, 1:2]
  
  # Find direction in which line joining start and end points is heading
  
  phi <- atan2(r1[2] - r2[2], r1[1] - r2[1])
  
  # Find shift vector of length d perpendicular to r1 - r2 vector
  
  shift <- d * c(sin(phi), -cos(phi))
  
  # Find coordinates of the polygon of width 2d with centerline at the base
  # projection
  
  coords <- rbind(r1+shift, r1-shift, r2-shift, r2+shift, r1+shift)
  
  # Construct a polygon object
  
  polygon <- Polygon(coords)
  
  # Return the results in the form of a list
  
  return(list(line1, line2, base1, base2, dz, x_mid, y_mid, polygon))
}

# Apply above function on all polygons forming 3D profiles

profile_structure <- unlist(lapply(profiles$geometry, function(x) lapply(x,
                            function(x) convert_to_line(x,d))), recursive = FALSE)

# All polygons with no horizontal extent (and therefore no area) are removed
# from the dataset

invalid_entries <- unlist(lapply(profile_structure, function(x) is.na(x[1])))
profile_structure <- profile_structure[!invalid_entries]

# All information returned by convert_to_line function results are extracted
# to separate variables

profiles <- unlist(c(lapply(profile_structure, function(x) x[[1]]),
                     lapply(profile_structure, function(x) x[[2]])))
base1_values <- unlist(lapply(profile_structure, function(x) x[[3]]))
base2_values <- unlist(lapply(profile_structure, function(x) x[[4]]))
dz_values <- unlist(lapply(profile_structure, function(x) x[[5]]))
mid_points <- data.frame(x=unlist(lapply(profile_structure, function(x) x[[6]])),
                         y=unlist(lapply(profile_structure, function(x) x[[7]])))
polygons <- unlist(c(lapply(profile_structure, function(x) x[[8]])))

# profiles are grouped together to form Lines object, which is then converted
# to a SpatialLines object

profiles <- sapply(1:length(profiles), function(x) Lines(list(profiles[[x]]), as.character(x)))
profiles <- SpatialLines(profiles, proj4string=catchments@proj4string)

# polygons are grouped together to form Polygons object, which is then converted
# to a SpatialPolygons object

polygons <- sapply(1:length(polygons), function(x) Polygons(list(polygons[[x]]), as.character(x)))
polygons <- SpatialPolygons(polygons, proj4string=catchments@proj4string)

# Number of polygons belonging to each profile in profiles data fream is found

number_of_polygons <- lapply(profiles$geometry, function(x) length(x))

# Here profiles are converted to a Spatial Lines Data Frame

# Initially data frame consisting of aquifer types (AQUIFER_D) is constructed

profiles_sldf <- data.frame(AQUIFER_D=profiles$AQUIFER_D)

# Each entry is duplicated by the number of polygons belonging to given profile
# (each triangular polygon will have a separate row in profiles_sldf)

profiles_sldf <- profiles_sldf[rep(1:nrow(profiles_sldf), number_of_polygons), ]

# Remove invalid entires (polygons with no horizontal extent)

profiles_sldf <- profiles_sldf[!invalid_entries,]

# Add dz parameter returned by convert_to_line to the data frame

profiles_sldf[,"dz"] <- dz_values

# Each triangle is divided into two triangles, so the rows are duplicated

profiles_sldf <- rbind(profiles_sldf, profiles_sldf)

# Add base length returned by convert_to_line to the data frame and use it
# to calculate polygon's area (from equation for area of triangle)

profiles_sldf[,"base"] <- c(base1_values, base2_values)
profiles_sldf[,"area"] <- (profiles_sldf$base * profiles_sldf$dz) / 2

# Assign an id to each profile

profiles_sldf[,"profile_id"] <- 1:nrow(profiles_sldf)

# Name rows with consecutive integer numbers 

row.names(profiles_sldf) <- as.character(1:nrow(profiles_sldf))

# Form a SpatialLinesDataFrame by adding data frame constructed above
# to the SpatialLines object constructed above

profiles <- SpatialLinesDataFrame(profiles, profiles_sldf)

# Duplicate mid_points (so each midpoint corresponds to one triangular polygon)

mid_points <- rbind(mid_points, mid_points)

# Form SpatialPoints object out of the mid points

mid_points <- SpatialPoints(mid_points, proj4string=catchments@proj4string)

# The constructed variables are save to a back up file. One can run the
# code from this load.image line if the computations had to be interrupted

save(list=c("profiles", "mid_points", "polygons", "profile_structure"),
     file='aquifer_processing.RData')
load.image(file='aquifer_processing.RData')

# Find number of polygons

n_polygons <- length(polygons)

# Entire UK will be divided into nx x ny tiles.
# By investigating intersection between polygons and catchments belonging to
# a single time at a time, it is possible to limit number of polygons being
# intersected, and therefore reduce computation time. Catchments belonging to
# the given tile are picked based on their centroids.

nx <- 4
ny <- 8

# Find catchments centroids

centroid <- coordinates(catchments)

# Find extent of all centroids

UK_extent <- extent(centroid)

# Initialize summary data frame in which length of profiles and mean depth
# of different types of aquifer is found.

summary <- data.frame(id=catchments$STATION, length=0)
for (type in unique(profiles$AQUIFER_D)) {
  summary[, type] <- 0
}

# Iterate over all tiles

for (x in 1:nx) {
  for (y in 1:ny) {
    
    # Report current progress
    
    print(paste('Processing tile ', (x - 1) * ny + y, '/', nx * ny), quote=FALSE)
    
    # Find tile's geographical extent
    
    xmin <- UK_extent@xmin + (x - 1) * (UK_extent@xmax - UK_extent@xmin) / nx
    xmax <- UK_extent@xmin + x * (UK_extent@xmax - UK_extent@xmin) / nx
    ymin <- UK_extent@ymin + (y - 1) * (UK_extent@ymax - UK_extent@ymin) / ny
    ymax <- UK_extent@ymin + y * (UK_extent@ymax - UK_extent@ymin) / ny
    
    # Find catchments, which centroids belong to the given tile  
    
    ids <- as.numeric(which(xmin <= centroid[,1] & centroid[,1] <= xmax &
                            ymin <= centroid[,2] & centroid[,2] <= ymax))
    
    # If no catchments belong to this tile proceed to the next tile
    
    if (sum(ids) == 0) next
    
    # Extract subset of catchments belonging to the given tile
    
    catchments_in_tile <- catchments[ids,]
    
    # Find geogeaphical extent of these catchments
    
    tile_extent <- extent(catchments_in_tile)
    
    # Reduce number of profiles to only those within extent computed above
    
    profiles_in_tile <- crop(profiles, tile_extent)
    
    # For each catchment in the tile mean depth of each aquifer type is computed
    
    for (i in 1:nrow(catchments_in_tile)) {
      catchment <- catchments_in_tile[i,]
      
      # Report the current progress
      
      print(paste('Processing catchment', i, '/', nrow(catchments_in_tile)), quote=FALSE)
      
      # Check which profiles overlap with the given catchment
      
      valid_profiles_ids <- over(profiles_in_tile, catchment)
      valid_profiles_ids <- !is.na(valid_profiles_ids$STATION)
      
      # If no profiles do overlap continue to the next catchment
      
      if (sum(valid_profiles_ids)==0) next
      
      # Extract the profiles
      
      profiles_subset <- profiles_in_tile[valid_profiles_ids,]
      
      # Use ids of extracted profiles to find ids of the polygons, corresponding
      # to these profiles (i.e. polygons which overlap with the catchment)
      
      valid_profiles_ids <- profiles$profile_id %in% profiles_subset$profile_id
      valid_polygons_ids <- valid_profiles_ids[1:n_polygons] |
        valid_profiles_ids[(n_polygons+1):(2*n_polygons)]
      polygons_subset <- polygons[valid_polygons_ids,]
      
      # intersect polygons and profiles with the catchment boundaries
      
      polygons_subset <- intersect(polygons_subset, catchment)
      profiles_subset <- intersect(profiles_subset, catchment)
      
      # Estimate the profile length by calculating union of all polygons and
      # dividing it by 2d (see report for the justification)
      
      polygon_union <- st_union(st_as_sf(polygons_subset))
      profiles_total_length <- as.numeric(st_area(polygon_union) / (2 * d))
      
      # Find ids of profiles intersecting with the catchment boundary
      # (only part of their area contribute to estimating mean depth)
      
      boundary_profiles <- over(profiles_subset, as(catchment, 'SpatialLinesDataFrame'))
      boundary_profiles <- !is.na(boundary_profiles$STATION)
      
      # Extract mid points belonging to boundary profiles and check whether
      # they are located outside the catchment (in this case different formula
      # for area is used)
      
      mid_points_subset <- profiles_subset[boundary_profiles,]$profile_id
      mid_points_subset <- mid_points[mid_points_subset,]
      mid_points_outside <- over(mid_points_subset, catchment)
      mid_points_outside <- is.na(mid_points_outside$STATION)
      
      # Calculate new base of rectangle (length of projection on xy plane)
      # and compare it to the original base before applying intersection
      
      base_new <- unlist(lapply(profiles_subset[boundary_profiles,]@lines, LinesLength))
      base_new <- base_new / profiles_subset[boundary_profiles,]$base
      
      # Check area of each polygon located at the boundary. In case of polygons
      # laying on the boundary, compute fraction of the area located inside
      # the catchment. and overwrite its previous value in the profile_subset.
      
      area_new <- profiles_subset[boundary_profiles,]$area
      area_new[mid_points_outside] <- area_new[mid_points_outside] *
                                      base_new[mid_points_outside] ^ 2
      area_new[!mid_points_outside] <- area_new[!mid_points_outside] *
                                       (1 - (1 - base_new[!mid_points_outside])^2)
      profiles_subset$area[boundary_profiles] <- area_new
      
      # Find total area of the profiles belonging to each AQUIFER_D class
      
      aquifer_summary <- aggregate(area ~ AQUIFER_D, profiles_subset, sum)
      
      # Check location of the catchment in summary data frame
      
      row_id <- which(summary$id == catchment$STATION)
      
      # Add total profile length to the data frame
      
      summary[row_id, 'length'] <- summary[row_id, 'length'] + profiles_total_length
      
      # Add total polygon area of each AQUIFER_D class to the data frame
      
      for (j in 1:nrow(aquifer_summary)) {
        aquifer_type <- aquifer_summary$AQUIFER_D[j]
        summary[row_id, aquifer_type] <- summary[row_id, aquifer_type] + aquifer_summary$area[j]
      }
    }
  }
}

# Find mean depth of each AQUIFER_D class by dividing corresponding
# total profile area by total profile length

# If profile length is zero the set mean depth to NA

summary[,unique(profiles$AQUIFER_D)] <- summary[,unique(profiles$AQUIFER_D)] / summary$length
summary[summary$length==0, unique(profiles$AQUIFER_D)] <- NA

# Export the results to output_file

write.csv(summary, output_file, row.names=FALSE)
