## ---------------------------
##
## Script name: NRFA_processing.R
##
## Purpose of script: Extracting physical parameters for UK catchments
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
##  Script requires web acces to download NRFA dataset with catchment and
##  gauging station information.
##
##  From this dataset mean rainfall, runoff (river flow per unit area) and
##  mean evapotranspiration (mean rainfall - mean runoff) for each catchment
##  is calculated and exported to 'output//NRFA_summary.csv' file.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# File to which the results are saved

output_file <- 'output//NRFA_summary.csv'

# Set an appropriate working directory

if (!require('htmltab', character.only=T)) {install.packages('htmltab');require('htmltab')}

# download dataset from the following url

url <- "https://nrfaapps.ceh.ac.uk/nrfa/ws/station-info?station=*&format=html&fields=all"
catchments <- htmltab(doc = url, which = 1, encoding = "UTF-8")

# extaract information about area covered with different surface types

terrain_data <- data.frame(arable = catchments$'lcm2007-arable-horticultural',
                           grassland = catchments$'lcm2007-grassland',
                           mountain = catchments$'lcm2007-mountain-heath-bog',
                           urban = catchments$'lcm2007-urban',
                           woodland = catchments$'lcm2007-woodland')

# convert dataset entries from character to numeric value

for (i in 1:ncol(terrain_data)) {
  terrain_data[,i] <- as.numeric(terrain_data[,i])
}

# replace missing entires with zeros

terrain_data[is.na(terrain_data)] <- 0

# fractional area covered by each type of terrain is caluclated

terrain_data <- terrain_data / apply(terrain_data, 1, sum)

# Mean Manning's n value is estimated as weighted average of Manning's n
# value for each type of the terrain. There values are taken from:
# https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/HecRAS/NEDC/lectures/docs/Manning%92s%20n-values%20for%20Kansas%20Dam%20Breach%20Analyses%20-%20Adopted%20071216.pdf

Manning_n  <- 0.035 * terrain_data$arable + 0.035 * terrain_data$grassland +
  0.025 * terrain_data$mountain + 0.15 * terrain_data$urban + 0.16 * terrain_data$woodland

# The catchment area is converted from km^2 to m^2

area_m2 <- 1e6 * as.numeric(catchments$`catchment-area`)

# Mean rainfall is taken from data from 1961-1990 and is converted from mm/year to m/s

rainfall_ms <- 1e-3 * as.numeric(catchments$`saar-1961-1990`) / (365.24 * 24 * 3600)

# If rainfall is negative set its value to NA

rainfall_ms[rainfall_ms < 0] <- NA

# Mean runoff per unit area [m/s] is calculated by dividing mean flow [m^3/s]
# by the catchment's area [m^2]

runoff_ms <- as.numeric(catchments$`gdf-mean-flow`) / area_m2

# Mean evapotranspiration is estimated as difference between mean rainfall and
# mean runoff

evap_ms <- rainfall_ms - runoff_ms

# Estimate channel dimensions from mean gauged flow using power laws derived in
# "A Study of the bank-full discharges of rivers in England and Wales"
# by M. Nixon (1959)

channel_width <- 0.911 * as.numeric(catchments$`gdf-mean-flow`) ^ (1/2)
channel_depth <- 0.545 * as.numeric(catchments$`gdf-mean-flow`) ^ (1/3)

# All parameters obtained are saved in a single data frame

summary <- data.frame(id=catchments$id,
                      area=area_m2,
                      rainfall=rainfall_ms,
                      runoff=runoff_ms,
                      evapotransmission=evap_ms,
                      Manning_n_hillslope=Manning_n,
                      channel_width=channel_width,
                      channel_depth=channel_depth)

# Results are saved in output_file

write.csv(summary, output_file, row.names=FALSE)