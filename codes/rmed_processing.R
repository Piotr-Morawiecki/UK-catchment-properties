## ---------------------------
##
## Script name: rmed_processing.R
##
## Purpose of script: Extracting median of annual peak rainfall for UK catchments
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
##  Script requires web acces to download NRFA datasets with daily gauged rainfall.
##
##  From this dataset historical daily rainfall data is read, processes
##  and exported to 'output//RMED_summary.csv' file.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Set an appropriate working directory

setwd('C://Users//pwm27//Desktop//PhD_floods//codes//data_extraction')

# File to which the results are saved

output_file <- 'output//Rmed_summary.csv'

# Set an appropriate working directory

if (!require('rjson', character.only=T)) {install.packages('htmltab');require('rjson')}

# download list of all catchments

link <- "https://nrfaapps.ceh.ac.uk/nrfa/ws/station-ids?format=json-object"
ids <- fromJSON(paste(readLines(link), collapse=""))[[1]]

# initialize data frame for saving Rmed values

Rmed_data <- data.frame(id=ids, Rmed=NA)

# start a clock to track computation time

start_time <- Sys.time()

# For each catchment Rmed is estimated

for (i in 1:length(ids)) {
  
  # For a given catchment download the rainfall time series
  
  id <- ids[i]
  rain_link <- "https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series?format=nrfa-csv&data-type=cdr&station="
  link <- paste(rain_link, id, sep='')
  rain <- read.csv(url(link))
  
  # Check if file includes any data (data are recorded starting from 19th row)
  
  if (nrow(rain) >= 19) {
    
    # Remove the file prefix (first 18 rows)
    
    rain <- rain[19:nrow(rain),1:2]
    
    # Name the columns and reformat data into a date and numeric format
    
    colnames(rain) <- c('date', 'daily_rain')
    rain$date <- as.Date(rain$date, "%Y-%m-%d")
    rain$daily_rain <- as.numeric(rain$daily_rain)
    
    # Extract year from each date
    
    rain[,"year"] <- format(rain$date, format="%Y")
    
    # For each year find maximum daily rainfall
    
    rain_max <- aggregate(daily_rain~year, rain, max)
    
    # Calculate the median of maximum annual rainfall (Rmed)
    
    Rmed_data$Rmed[i] <- median(rain_max$daily_rain)
  }
  
  # Estimate remaining computation time and print the computation progress
  
  time_passed <- as.numeric(difftime(Sys.time(), start_time), units='mins')
  remaining_time <- round(time_passed / i * (length(ids) - i), digits=2)
  print(paste('Catchment', i, '/', length(ids),'anaysed.'), quote=FALSE)
  print(paste('Remaining time:', remaining_time,'minutes.'), quote=FALSE)
}

# Save constructed data frame in the output file

write.csv(Rmed_data, output_file)