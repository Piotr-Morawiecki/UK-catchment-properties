## ---------------------------
##
## Script name: merge_all.R
##
## Purpose of script: Merging summary with physical parameters
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
##  Script requires the following input files:
##
##    - output/NRFA_summary.csv
##      (generated using NRFA_processing.R script)
##
##    - output/streamline_summary.csv
##      (generated using streamline_processing.R script)
##
##    - output/river_summary.CSV
##      (generated using river_processing.R script)
##
##    - output/river_2_summary.CSV
##      (generated using river_processing_2.R script)
##
##    - output/aquifer_summary.csv
##      (generated using aquifer_processing.R script)
##
##    - output/soil_summary.csv
##      (generated using soil_processing.R script)
##
##    - output/permeability_summary.csv
##      (generated using permeability_processing.R script)
##
##    - output/hydrogeology_summary.csv
##      (generated using hydrogeology_processing.R script)
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Read all input files (see section 'Notes' in the header above)

NRFA_summary <- read.csv('output//NRFA_summary.csv')
#streamline_summary <- read.csv('output//streamline_summary.csv')

#streamline_summary$width_B[is.na(streamline_summary$gradient_perpendicular)] <- NA
streamline_summary <- read.csv('output//streamline_summary_corrected.csv')

longest_river_summary <- read.csv('output//longest_river_processing_summary.CSV')
river_summary <- read.csv('output//river_summary.CSV')
river2_summary <- read.csv('output//river_2_summary.CSV')
aquifer_summary <- read.csv('output//aquifer_summary.csv')
soil_summary <- read.csv('output//soil_summary.csv')
aquifer_properties <- read.csv('output//permeability_summary.csv')
hydrogeology_summary <- read.csv('output//hydrogeology_summary.csv')

## ---------------------------

# For soil_summary aquifer_thickness is calculated by adding all depth of
# selected aquifer_types. Here 'Principal' and 'SecondaryA' are picked as
# two most conductive zones.
# The mean depth is calculated as long as the length of profiles from which
# data was extracted per unit area is higher than given threshold value.

# Define threshold value

profile_threshold <- 1e-5

# Define aquifer types which contribute to aquifer depth

aquifer_types <- c('Principal', 'SecondaryA')

# Area of catchment from NRFA summary is added to aquifer_summary

aquifer_summary <- merge(aquifer_summary, NRFA_summary[,c('id','area')])

# Length of profiles per unit area is calulcated for each catchment

aquifer_summary$lengthPerArea <- aquifer_summary$length / aquifer_summary$area

# Check which entries are below threshold

entries_below_threshold <- aquifer_summary$lengthPerArea < profile_threshold

# Calculate total depth of chosen layer and overwrite aquifer_summary with them

thickness <- apply(aquifer_summary[,aquifer_types], 1, sum)
aquifer_summary <- data.frame(id=aquifer_summary$id,
                              aquifer_thickness=thickness)

# Assign NA to entires below the threshold level

aquifer_summary$aquifer_thickness[entries_below_threshold] <- NA

## ---------------------------

# From soil_summary dataset only entries from a given depth are extracted

# Choose depth from which entries are recorded

depth <- 30

# extract rows with given depth value

soil_summary <- soil_summary[soil_summary$depth==depth,]

# Construct data frame with new parameter names

soil_summary <- data.frame(id=soil_summary$catchment,
                           hydraulic_conductivity_A=soil_summary$KS,
                           MvG_alpha=soil_summary$HCC_alp,
                           MvG_thetaR=soil_summary$HCC_thr,
                           MvG_thetaS=soil_summary$HCC_ths,
                           MvG_n=soil_summary$HCC_n)

## ---------------------------

# Merge all datasets into a single data frame

summary <- Reduce(function(x,y) merge(x = x, y = y, by = "id", all = TRUE),
                  list(NRFA_summary, streamline_summary, river_summary,
                       river2_summary, aquifer_summary, soil_summary,
                       longest_river_summary))

# Estimate catchment width by dividing area (originally from NRFA_summary) by
# total stream length (originally from river2_summary)

summary[,'width_A'] <- summary$area / summary$length_B

# If total stream length is 0 set NA as catchment width

summary$width_A[summary$length_B==0] <- NA



# url <- "https://nrfaapps.ceh.ac.uk/nrfa/ws/station-info?station=*&format=html&fields=all"
# catchments <- htmltab(doc = url, which = 1, encoding = "UTF-8")
# temp <- data.frame(id=catchments$id,
#                    x=catchments$easting,
#                    y=catchments$northing)
# temp <- merge(summary, temp)
# temp <- temp[is.na(summary$aquifer_thickness),]
# plot(temp$x, temp$y)



# Remove all intermediate summaries from the variables

rm(list=setdiff(ls(), c("summary","aquifer_properties","hydrogeology_summary")))

write.csv(summary, 'output/merged_summary.csv', row.names=FALSE)

## ---------------------------

library(stringr)
str_replace_all(string, pattern, replacement)

format_number <- function(x, sf) {
  if (x > 1000) x <- round(x, digits=sf-1-floor(log10(x)))
  x <- format(x, digits=3)
  x <- paste('$', x, '$', sep='')
  return(x)
}

for (var in colnames(summary)) {
  values <- summary[,var]
  values <- values[is.finite(values)]
  q25 <- format_number(quantile(values, probs=0.25), sf=3)
  med <- format_number(median(values), sf=3)
  q75 <- format_number(quantile(values, probs=0.75), sf=3)
  print(paste(var, q25, med, q75, sep=' & '),quote=FALSE)
}
summary_statistics <- summary(summary)
summary_statistics

area <- apply(aquifer_properties[,2:6], 1, sum)
aquifer_properties[,2:6] <- aquifer_properties[,2:6] / area
aquifer_properties[,7:11] <- aquifer_properties[,7:11] / area

area <- apply(hydrogeology_summary[,2:7], 1, sum)
hydrogeology_summary[,2:7] <- hydrogeology_summary[,2:7] / area
hydrogeology_summary[,8:10] <- hydrogeology_summary[,8:10] / area

hydrogeology_summary$Rocks.with.essentially.no.groundwater


temp <- merge(summary, hydrogeology_summary)


plot(temp$gradient_perpendicular, temp$width_B, log='xy', pch=16)
highlighted = temp$Rocks.with.essentially.no.groundwater>0.9
points(temp$gradient_perpendicular[highlighted], temp$width_B[highlighted],
     log='xy', col='red', pch=16)


temp$param <- temp$gradient_perpendicular / temp$width_B
temp$param <- temp$gradient_perpendicular * 1 * 1e-4 / (temp$width_B * temp$runoff)

mean(temp$gradient_perpendicular[highlighted], na.rm=TRUE)
mean(temp$width_B[highlighted], na.rm=TRUE)
mean(temp$runoff[highlighted], na.rm=TRUE)
mean(temp$param[highlighted], na.rm=TRUE)
summary(temp[highlighted,])

(1e-4*4e-2)/(6e2*8e-9)
