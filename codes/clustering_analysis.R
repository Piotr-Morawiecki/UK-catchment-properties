## ---------------------------
##
## Script name: clustering_analysis.R
##
## Purpose of script: postprocessing using correlation, cluster analysis and PCA 
##
## Author: Piotr Morawiecki
##
## Date Created: 2022-04-25
##
## Copyright (c) Piotr Morawiecki, 2021
## Email: pwm27@bath.ac.uk
##
## ---------------------------
##
## Notes:
##  
##  Script requires two input datasets:
##    - output/merged_summary.csv generated with codes/merge_all.R
##    - data/NRFA_catchment_boundaries/nrfa_public_ex20200925.shp
##      (more information: https://nrfa.ceh.ac.uk/content/catchment-boundary-and-areas)
##
##  Script uses merged summary with all catchment marameters to perform standard
##  statistical analysis, including:
##    - correlation analysis, showing correlation between all variables,
##    - cluster analysis, showing main types of catchments in UK,
##    - principal component analysis (PCA) highlighting key differences
##      between the clusters in terms of the most varying combinations of parameters
##
##  Multiple output figures and files are generated.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Directory to the input datasets (see section 'Notes' in the header above)

summary_file <- 'output/merged_summary.csv'
catchment_filename <- "data//NRFA_catchment_boundaries//nrfa_public_ex20200925.shp"

# Specify number of clusters to be generated

n_clusters <- 3

# Specify which variables should be used for clustering and PCA;
# Label should match with column names in summary_file

variables <- c("rainfall", "runoff", "Manning_n_hillslope", "width_B",
               "gradient_perpendicular", "gradient_parallel",
               "aquifer_thickness", "hydraulic_conductivity_A", "MvG_alpha",
               "MvG_thetaS", "MvG_n")

# Names of output files

corrplot_file <- 'figures//corr_plot.pdf'
clustering_file <- 'figures//clustering_results.pdf'
corr_file <- 'output//corr_matrix.csv'
PCA_file <- 'output//pca_components.csv'
PCA_diagram_file <- 'figures//PCA_diagram.pdf'

## ---------------------------

# Installing required libraries

libraries_required <- c('corrplot', 'RColorBrewer', 'rgdal', 'rgeos', 'raster', 'colorspace')
lapply(libraries_required, function(x) {
  if (!require(x, character.only=T)) {install.packages(x);require(x)}
})

# Read summary with catchment parameters

summary = read.csv(summary_file)


# PART 1. CORRELATION ANALYSIS

# Calculate correlation matrix and save it

M <- cor(summary[,2:ncol(summary)], use='pairwise.complete.obs')
write.csv(M, corr_file)

# Represent it graphically and save it

pdf(corrplot_file)
par(mfrow=c(1,1), mar = c(5.1,4.1,4.1,2.1))
corrplot(M, type="upper", order="hclust", col=brewer.pal(n=8, name="RdYlBu"))
dev.off()


# PART 2. CLUSTER ANALYSIS

# Create reduced summary consisting only chosen set of variables

summary_reduced <- summary[,c('id', variables)]

# Remove rows with incomplete set of variables

summary_reduced <- summary_reduced[apply(summary_reduced, 1, function(x) sum(is.na(x))==0),]

# Use k-means clustering algorithm to assign catchments to given number of clusters

fit <- kmeans(scale(summary_reduced[,variables]), n_clusters)

# Print mean value of parameters belonging to each cluster 

print(aggregate(summary_reduced[,variables],by=list(fit$cluster),FUN=mean))

# Add cluster id to the summary_reduced data frame

summary_reduced[,'cluster'] <- fit$cluster

# Calculate number of catchments belonging to each cluster

n_catchments <- sapply(1:n_clusters, function(x) sum(fit$cluster==x))

# Plot clusters to pdf file

pdf(clustering_file, width=10, height=8)
par(mfrow=c(3,4), mar = c(3,3,1,2))

# Firstly column plot with number of catchments in each cluster

barplot(n_catchments~cluster,
        data.frame(cluster=1:n_clusters, n_catchments=n_catchments),
        xlab="", ylab="number of catchments", main="(A) Cluster size")

# Then boxplots with distribution of each variable

i <- 1
for (var in variables) {
  i <- i + 1
  boxplot(as.formula(paste(var,'cluster',sep='~')), data=summary_reduced,
          xlab="", ylab=var, main=paste('(',LETTERS[i],') ', var, sep=''))
}

dev.off()


# PART 3. PRINCIPAL COMPONENT ANALYSIS (PCA)

# Perform PCA on parameters used for clustering (after normalization)

pca_components <- prcomp(summary_reduced[, variables], scale = TRUE)

# Add PCA components to reduced summary forming summary_pca data frame

summary_pca <- data.frame(summary_reduced, pca_components$x)

# Extract weighs defining each principal components and write them to a csv file   

pca_components <- pca_components$rotation
write.csv(pca_components, 'pca_components.csv')

# Now we want to plot a UK map showing location of clusters.

# Read shapefile with catchment boundaries

catchment_boundaries <- readOGR(dsn = catchment_filename, stringsAsFactors = F)

# Removes ring self-intersections from polygons

catchment_boundaries <- gBuffer(catchment_boundaries, byid=TRUE, width=0)

# Order catchment boundaries by id number and merge with PCA summary

catchment_boundaries <- catchment_boundaries[order(catchment_boundaries$STATION),]
catchment_boundaries@data <- merge(data.frame(id=catchment_boundaries$STATION),
                                   summary_pca[,c('id', 'cluster', 'PC1', 'PC2')],
                                   all.x=TRUE, sort=TRUE)

# simplify catchment boundaries using Douglas-Peuker algorithm with 100m
# tolerance. It significantly reduces the output figures size

simple_boundaries <- gSimplify(catchment_boundaries, tol=100, topologyPreserve=TRUE)

# Create new spatial polygons data frame, this time with simplified boundaries

simple_boundaries <- SpatialPolygonsDataFrame(simple_boundaries, catchment_boundaries@data, match.ID=FALSE)

# Select only boundaries, which has been assigned to any cluster

simple_boundaries <- simple_boundaries[!is.na(simple_boundaries$cluster),]

# Download UK boundary from GADM database of global administrative boundaries

UK_boundary <- getData('GADM', country='United Kingdom', level=0)

# Transform UK_boundary coordinate system to be consistent with catchment_boundaries 

UK_boundary <- spTransform(UK_boundary, crs(catchment_boundaries))

# Also simplify it using Douglas-Peuker algorithm with 100m tolerance

UK_boundary <- gSimplify(UK_boundary, tol=100, topologyPreserve=TRUE)

# Define color palette - we use Brewer palette Green-Blue

colors_palette <- brewer.pal(n_clusters, "GnBu")

# Now we create 2x2 figure.

pdf(PCA_diagram_file, width=10, height=10)
par(mfrow=c(2,2))

# In the top row we show location of all catchment on PC1 vs PC2 diagram and UK map.

plot(summary_pca$PC1, summary_pca$PC2, col=colors_palette[summary_pca$cluster],
     pch=16, xlab='1st principal component', ylab='2nd principal component')
legend(3,3,lapply(1:n_clusters, function(x) paste('Cluster', x)),
       col=colors_palette,pch=16)

colors <- colors_palette[simple_boundaries$cluster]
plot(UK_boundary, xlim=c(0,656e3), ylim=c(5e3, 1296e3))
plot(simple_boundaries, col=colors, border=darken(colors, amount=0.3), add=TRUE)

# In the second row we show with distribution of PC1 and PC2 values across UK

colors_palette <- brewer.pal(7, "GnBu")
colors <- cut(simple_boundaries$PC1, breaks=seq(-6,8,2), labels=1:7)
colors <- colors_palette[colors]
plot(UK_boundary, xlim=c(0,656e3), ylim=c(5e3, 1296e3))
plot(simple_boundaries, col=colors, border=darken(colors, amount=0.3), add=TRUE)

colors_palette <- brewer.pal(5, "GnBu")
colors <- cut(simple_boundaries$PC2, breaks=seq(-6,4,2), labels=1:5)
colors <- colors_palette[colors]
plot(UK_boundary, xlim=c(0,656e3), ylim=c(5e3, 1296e3))
plot(simple_boundaries, col=colors, border=darken(colors, amount=0.3), add=TRUE)

dev.off()

# Version with a proper legend is available in our paper (see main directory).