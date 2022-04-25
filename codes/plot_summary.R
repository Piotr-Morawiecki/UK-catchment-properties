## ---------------------------
##
## Script name: plot_summary.R
##
## Purpose of script: Plotting maps of catchment parameters
##
## Author: Piotr Morawiecki
##
## Date Created: 2022-04-25
##
## Copyright (c) Piotr Morawiecki, 2022
## Email: pwm27@bath.ac.uk
##
## ---------------------------
##
## Notes:
##  
##  Script requires one input datasets:
##    - output/merged_summary.csv generated with codes/merge_all.R
##
##  Sript takes parameters values saved in 'merged_summary.csv' file
##  and for selected subset of parameters plots maps showing their value
##  across UK. The plots are saved to 'figures/' directory.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Directory to the input datasets (see section 'Notes' in the header above)

summary_filename <- "output/merged_summary.csv"
catchment_boundaries_filename <- "data//NRFA_catchment_boundaries//nrfa_public_ex20200925.shp"

# Directory to which the results are saved

output_directory <- "figures//"

# List of parameters to be plotted - names have to agree with columns names
# as specified in summary_filename

parameters <- c('rainfall', 'runoff', # 'evapotransmission',
                'density_A',
                'density_B', 'width_B', 'gradient_parallel',
                'gradient_perpendicular', 'aquifer_thickness',
                'hydraulic_conductivity_A', 'MvG_alpha', 'MvG_thetaR',
                'MvG_thetaS', 'MvG_n', 'Manning_n_hillslope')

# List of corresponding titles for the plots

titles <- c('Rainfall intensity', 'Mean runoff', # 'Mean evapotransmission',
            'Main streams density', 'All streams density', 'Catchment width',
            'Gradient along hillslope', 'Gradient along river',
            'Aquifer thickness', 'Hydraulic conductivity',
            expression(paste("MvG", alpha, "parameter")),
            'Residual water content', 'Saturated water content',
            'MvG n parameter', 'Manning\'s n for the surface' )

## ---------------------------

# Installing required libraries

libraries_required <- c('ggplot2', 'gridExtra', 'latex2exp', 'sp', 'rgdal', 'rgeos', 'raster')
lapply(libraries_required, function(x) {
  if (!require(x, character.only=T)) {install.packages(x);require(x)}
})

# Read summary with catchments parameters
summary <- read.csv(summary_filename)

# Read catchments boundaries from the file
catchment_boundaries <- readOGR(dsn = catchment_boundaries_filename, stringsAsFactors = F)

# Remove ring self-intersections from polygons
catchment_boundaries <- gBuffer(catchment_boundaries, byid=TRUE, width=0)

# add summary data frame to the catchment_boundaries spatial polygons data frame
catchment_boundaries@data <- merge(data.frame(id=catchment_boundaries$STATION),
                                   summary, all.x=TRUE, sort=FALSE)

# remove catchments with unspecified area
catchment_boundaries <- catchment_boundaries[!is.na(catchment_boundaries$area),]

# simplify catchment boundaries using Douglas-Peuker algorithm with 100m
# tolerance. It significantly reduces the output figures size
simple_boundaries <- gSimplify(catchment_boundaries, tol=100, topologyPreserve=TRUE)

# Create new spatial polygons data frame, this time with simplified boundaries
simple_boundaries <- SpatialPolygonsDataFrame(simple_boundaries, catchment_boundaries@data, match.ID=FALSE)

# Download UK boundary from GADM database of global administrative boundaries
UK_boundary <- getData('GADM', country='United Kingdom', level=0)

# Transform UK_boundary coordinate system to be consistent with catchment_boundaries 
UK_boundary <- spTransform(UK_boundary, crs(catchment_boundaries))

# Convert simple_boundaries spatial polygons to a data frame
# (this is the format that ggplot requires for plotting spatial datasets)
fortified_data <- fortify(simple_boundaries, region = "id")

# Add all polygons data to the constructed data frame 
fortified_data <- merge(fortified_data, simple_boundaries@data, by = "id")

# Reorder data in order of increasing area - this way smaller catchments
# are plotted on top and won't be covered by large catchments
fortified_data <- fortified_data[rev(order(fortified_data$area)),]
fortified_data$group <- ordered(fortified_data$group, unique(fortified_data$group))

# In a simple way UK_boundary is converted into data frame
UK_boundary <- fortify(UK_boundary)

# plot_map function takes two input parameters:
#   parameter: name of parameter to be plotted (as it appears in fortified_data)
#   title: title of the graph in form of a string
# function returns map showing values of given parameter for all UK catchments

plot_map <- function(parameter, title) {
  
  # Copy fortified_data data frame adding chosen parameter to plotted_variable
  plotted_data <- fortified_data
  plotted_data[, 'plotted_variable'] <- plotted_data[,parameter]
  
  # Catchments for which given parameter is not specified are removed
  # from the dataset, and hence are not plotted.
  plotted_data <- plotted_data[!is.na(plotted_data$plotted_variable),]
  
  # The plotted variable is cut into five intervals, each of which
  # will be plotted with different color
  plotted_data$plotted_variable <- cut(plotted_data$plotted_variable, breaks=5)
  
  # ggplot is used to produce a plot with UK boundary and catchments outlines
  # colored appropriately on top of it
  map <- ggplot(plotted_data,
                aes(x=long, y=lat, group=group, fill=plotted_variable)) + 
    geom_polygon(data=UK_boundary, fill='white', col='black') +
    geom_polygon(col=NA) + ggtitle(title) + theme_void() + coord_fixed() +
    xlim(0, 660e3) + ylim(0, 1200e3) +
    theme_minimal() +
    scale_fill_brewer(palette="GnBu", name = '', na.value="white") +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_blank())
  
  # map is returned, so that it can be later plotted or saved to a file
  return(map)
}

# plot_map function is applied to all parameters specified in the beginning of the script
plots <- lapply(1:length(parameters), function(i) plot_map(parameters[i], titles[i]))

# Plot grid with all plots on the screen 
do.call(grid.arrange,c(plots, ncol=7, nrow=2))

# Save grid of all plots into a pdf and png files 
g <- do.call(arrangeGrob,c(plots, ncol=7, nrow=2))
ggsave(paste(output_directory,'UK_map_grid.pdf',sep=''), g, width=5.31, height=4, units = "cm")
ggsave(paste(output_directory,'UK_map_grid.png',sep=''), g, width=29.7, height=21, units = "cm")

# To remove from repo (TikZ version)

plot_map_blank <- function(parameter, title) {
  plotted_data <- fortified_data
  plotted_data[, 'plotted_variable'] <- plotted_data[,parameter]
  plotted_data <- plotted_data[!is.na(plotted_data$plotted_variable),]
  plotted_data$plotted_variable <- cut(plotted_data$plotted_variable, breaks=5)
  map <- ggplot(plotted_data,
                aes(x=long, y=lat, group=group, fill=plotted_variable)) + 
    geom_polygon(data=UK_boundary, fill='white', col='black') +
    geom_polygon(col=NA) +
    coord_fixed() +
    xlim(0, 660e3) + ylim(0, 1200e3) +
    theme_nothing() +
    scale_fill_brewer(palette="GnBu", name = '', na.value="white")
  return(map)
}

for (i in 1:length(parameters)) {
  g <- plot_map_blank(parameters[i], titles[i])
  ggsave(paste("output//UK_map_grid_plot", i, ".png", sep=''), g, width=29.7, height=21, units = "cm")
}
