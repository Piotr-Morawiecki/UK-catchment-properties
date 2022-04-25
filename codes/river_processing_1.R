## ---------------------------
##
## Script name: river_processing_1.R
##
## Purpose of script: Extracting main rivers length and their mean elevation
##                    gradient for UK catchments
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
##  Script requires three input datasets:
##    - data/NRFA_catchment_boundaries/nrfa_public_ex20200925.shp
##      (more information: https://nrfa.ceh.ac.uk/content/catchment-boundary-and-areas)
##    - data/DIGIMAP_river_network/WatercourseLink.shp
##      (available here: https://www.ordnancesurvey.co.uk/business-government/products/open-map-rivers)
##    - data/DIGIMAP_elevation/...
##      (available here: https://www.ordnancesurvey.co.uk/business-government/products/terrain-50)
##
##  Script extracts total length of all rivers from WatercourseLink.shp file
##  and mean value of elevation gradient along these rivers for UK catchments.
##  The results are exported to 'output/river_processing_summary.csv' file.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Directory to the input datasets (see section 'Notes' in the header above)

boundaries_file <- 'data//NRFA_catchment_boundaries//nrfa_public_ex20200925.shp'
river_network_filename <- 'data//DIGIMAP_river_network//WatercourseLink.shp'
dtm_directory <- 'data//DIGIMAP_elevation//'

# File to which the results are saved

output_file <- "output//river_processing_summary.csv"

## ---------------------------

# Installing required libraries

libraries_required <- c('rgdal', 'rgeos', 'sp', 'raster', 'plyr', 'igraph')
lapply(libraries_required, function(x) {
  if (!require(x, character.only=T)) {install.packages(x);require(x)}
})

# Read input datasets

catchments <- readOGR(dsn = boundaries_file, stringsAsFactors = F)
rivers <- readOGR(dsn = river_network_filename, stringsAsFactors = F)

# Remove ring self-intersections, which may rarely occur in the input datasets

catchments <- gBuffer(catchments, byid=TRUE, width=0)

# find each stream start and end point coordinates

start <- as.data.frame(do.call(rbind, lapply(rivers@lines, function(x) x@Lines[[1]]@coords[1,])))
end <- as.data.frame(do.call(rbind, lapply(rivers@lines, function(x) x@Lines[[1]]@coords[nrow(x@Lines[[1]]@coords),])))

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

# Calculate length of each stream in rivers dataset

lengths <- as.data.frame(do.call(rbind, lapply(rivers@lines, compute_river_length)))

# Save length, start point and end point coordiantes in the data frame.
# Column dtm1 and dtm2 will be filled with dtm value at start and end point,
# respectivelly.

stream_data <- data.frame(length=lengths, x1=start, x2=end, dtm1=NA, dtm2=NA)

# Rename columns

colnames(stream_data) <- c('length', 'x1', 'y1', 'x2', 'y2', 'dtm1', 'dtm2')

# Read list of regions from the dtm_directory (ignoring 'docs' file if present)

regions <- list.files(dtm_directory)
regions <- regions[regions!='docs']

# function inside_extent checks whether given point is inside a given extent

inside_extent <- function(x, y, extent) x >= extent@xmin & x <= extent@xmax & y >= extent@ymin & y <= extent@ymax 

# For loop for each region finds dtm of start and end points of rivers
# laying within a given region

for (i in 1:length(regions)) {
  region <- regions[i]
  
  # Log the current progress
  
  print(paste('Region ', region, ' started (', i, '/', length(regions),').', sep=''), quote=FALSE)
  print('---------------------------', quote=FALSE)
  
  # Construct directory with DTM datasets for the given region
  
  directory <- paste(dtm_directory, tolower(region), sep='')
  
  # List all DTM_files for the given region
  
  dtm_files <- list.files(directory, pattern = "\\.asc$")
  
  # If no files were found proceed to the next region
  
  if (length(dtm_files) == 0) next
  
  # Iterate over all dtm files for the given region
  
  for (j in 1:length(dtm_files)) {
    print(paste('DTM', j, '/', length(dtm_files)), quote=FALSE)
    
    # import given dtm raster
    
    dtm <- raster(paste(directory, dtm_files[j], sep='//'))
    
    # check start points inside the raster 
    
    points_inside <- inside_extent(stream_data$x1, stream_data$y1, extent(dtm))
    
    # extract dtm value at these points and save in stream_data data frame
    
    stream_data$dtm1[points_inside] <- extract(dtm, stream_data[points_inside,c("x1", "y1")])
    
    # Repeat the last two steps with the endpoints
    
    points_inside <- inside_extent(stream_data$x2, stream_data$y2, extent(dtm))
    stream_data$dtm2[points_inside] <- extract(dtm, stream_data[points_inside,c("x2", "y2")])
  }
}

# stream_data is saved as csv file. If the script is interrupted one can
# continue from this line without a need to recompute the previous for loop.

write.csv(stream_data, 'output//river_gradient_data.csv', row.names = FALSE)
stream_data <- read.csv('output//river_gradient_data.csv')

# Elevation difference is calculated as difference between the start and end point

stream_data$elevation_diff <- abs(stream_data$dtm1 - stream_data$dtm2)

# TYpical value of gradient along the river can be estimated by dividing the
# total elevation difference by the length of all streams. Only streams with
# both dtm1 and dtm2 available are used to calculated this sum.

valid <- is.finite(stream_data$dtm1) & is.finite(stream_data$dtm2)
mean_gradient <- sum(stream_data$elevation_diff[valid]) / sum(stream_data$length[valid])
print(paste('Mean gradient along river:', mean_gradient), quote=FALSE)

# The stream_data are saved in rivers SpatialLinesDataFrame

rivers@data <- cbind(rivers@data, stream_data, by='id')

# Initialise summary data frame for extaracted total stream length (length_main)
# and mean gradient along river (gradient_parallel) for each catchment

summary <- data.frame(id=c(),
                      length_main=c(),
                      max_river_length=c(),
                      avg_segment_length=c(),
                      gradient_parallel=c())

# The for loop iterates over all catchments to find above parameters

for (i in 1:nrow(catchments)) {
  catchment <- catchments[i,]
  
  # Current progress is reported
  
  print(paste("Catchment", i, "/", nrow(catchments)))
  
  # Rivers overlapping with given catchments are found
  
  rivers_subset <- rivers[!is.na(over(rivers, catchment)$STATION),]
  
  # If no rivers were found set parameters in summary data frame to NA
  
  if (length(rivers_subset) == 0) {
    summary <- rbind(summary, data.frame(id=catchment$STATION,length_A=NA, gradient_parallel=NA))
    next
  }
  
  # Intersect overlapping rivers with the catchment boundaries 
  
  rivers_subset <- intersect(rivers_subset, catchment)
  
  # If no rivers are overlapping set parameters in summary data frame to NA
  
  if (length(rivers_subset) == 0) {
    summary <- rbind(summary, data.frame(id=catchment$STATION,length_A=NA, gradient_parallel=NA))
    next
  }
  
  # Now we calculate four catchment parameters, three measures of catchment
  # length and mean elevation gradient along the streams
  
  # 1. COMPUTING THE TOTAL LENGTH OF ALL STREAMS L_MAIN
  
  # Compute length of all streams inside the catchment
  # Note that streams passing through boundary (usually only one per catchment)
  # will be shorter than before intersection with catchment boundary
  
  lengths <- lapply(rivers_subset@lines, compute_river_length)
  lengths <- unlist(as.data.frame(do.call(rbind, lengths)))
  
  # Sum above lengths to find total stream length over the catchment
  
  length_main <- sum(lengths)
  
  # 2. COMPUTING THE LENGTH OF THE LONGEST STREAM L_MAX
  
  # Find all unique node ids (they are represented with strings)
  
  nodes <- unique(c(rivers_subset$startNode, rivers_subset$endNode))
  
  # Then convert them into natural numbers: 1, 2, 3, ... 
  
  start_node <- as.numeric(mapvalues(rivers_subset$startNode, nodes,
                                     1:length(nodes), warn_missing = FALSE))
  
  end_node <- as.numeric(mapvalues(rivers_subset$endNode, nodes,
                                   1:length(nodes), warn_missing = FALSE))
  
  # Represent river network in form of an adjacency matrix
  # (first create it, then assigned stream lengths between appropriate nodes)
  
  adjacency_matrix <- matrix(0, nrow=length(nodes), ncol=length(nodes))
  adjacency_matrix[cbind(start_node, end_node)] <- lengths
  
  # Construct graph object from adjacency matrix and find its diameter;
  # This way we find longest distance between the nodes, and since it is
  # a directed graph it corresponds to the longest stream in our catchment 
  
  graph <- graph_from_adjacency_matrix(adjacency_matrix, weighted=TRUE)
  max_river_length <- diameter(graph)
  
  # 3. COMPUTING THE MEAN LENGTH BETWEEN TRIBUTARIES L_TRIB
  
  # Remove segments of the stream, which are dead ends, from the subset,
  # so that only segments between river tributaries are left.
  
  rivers_subset_trib <- rivers_subset[start_node %in% end_node,]
  
  # Simplify the segments using Douglas-Peuker algorithm with 100m
  # tolerance; this way we get rid of small river meanders leaving only
  # bends longer than 100 meters.
  
  rivers_subset_trib <- gSimplify(rivers_subset_trib, tol=100)
  
  # compute lengths of all segments and find their mean value,
  # i.e. mean distance between tributaries
  
  lengths <- lapply(temp2@lines, compute_river_length)
  avg_segment_length <- mean(unlist(as.data.frame(do.call(rbind, lengths))))
  
  # Calculate mean gradient as sum of elevation difference over stream length
  
  gradient <- sum(rivers_subset$elevation_diff) / sum(rivers_subset$length)
  
  # Obtained estimates are added to the summary data frame
  
  summary <- rbind(summary, data.frame(id=catchment$STATION,
                                       length_main=length_main,
                                       length_max=max_river_length,
                                       length_trib=avg_segment_length,
                                       gradient_parallel=gradient))
  
  # The summary is saved in csv file. If calculations are interrupted one can
  # continue script by reading this csv file and running for loop from the
  # last i value.
  
  write.csv(summary, output_file, row.names = FALSE)
}

# Save the final version of the summary to the output_file.

write.csv(summary, output_file, row.names = FALSE)
