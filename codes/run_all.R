## ---------------------------
##
## Script name: run_all.R
##
## Purpose of script: Running other scripts in the repository
##
## Author: Piotr Morawiecki
##
## Date Created: 2022-04-24
##
## Copyright (c) Piotr Morawiecki, 2022
## Email: pwm27@bath.ac.uk
##
## ---------------------------
##
## Notes:
##   
##  Before running the codes make sure that all required input data
##  are available (see README for details). Aslo some scripts require
##  web acces to download NRFA data.
##
##  Scripts processing scripts can be run independently.
##  Postprocessing scripts require datasets produced by all processing scripts,
##  however their precomputed versions are already included in the repository.
##  This allows the postprocessing scripts to be run without going through
##  processing scripts, some of which take very long to compute.
##
## ---------------------------

# Set an appropriate working directory (works only when run from RStudio)
# Otherwise replace with setwd('full path to UK-catchment-properties' directory)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Set the codes you want to rerun to TRUE.

# Preprocessing script
run_streamline_extraction <- FALSE    # WARNING: computations take a lot of time

# Processing codes
run_aquifer_processing <- FALSE
run_hydrogeology_processing <- FALSE
run_land_cover_processing <- FALSE
run_NRFA_processing <- FALSE
run_river_processing_1 <- FALSE
run_river_processing_2 <- FALSE       # WARNING: computations take a lot of time
run_rmed_processing <- FALSE
run_soil_processing <- FALSE
run_streamline_processing <- FALSE

# Postprocessing scripts
run_merge_all <- TRUE
run_cluster_analysis <- TRUE
run_summary_plotting <- TRUE

# Selected codes are run
if (run_streamline_extraction)   source("codes//streamline_extraction.R", local=TRUE)
if (run_aquifer_processing)      source("codes//aquifer_processing.R", local=TRUE)
if (run_hydrogeology_processing) source("codes//hydrogeology_processing.R", local=TRUE)
if (run_land_cover_processing)   source("codes//land_cover_processing.R", local=TRUE)
if (run_NRFA_processing)         source("codes//NRFA_processing.R", local=TRUE)
if (run_river_processing_1)      source("codes//river_processing_1.R", local=TRUE)
if (run_river_processing_2)      source("codes//river_processing_2.R", local=TRUE)
if (run_rmed_processing)         source("codes//rmed_processing.R", local=TRUE)
if (run_soil_processing)         source("codes//soil_processing.R", local=TRUE)
if (run_streamline_processing)   source("codes//streamline_processing.R", local=TRUE)
if (run_summary_plotting)        source("codes//run_summary_plotting.R", local=TRUE)

if (run_merge_all) source("codes//merge_all.R", local=TRUE)
if (run_cluster_analysis) source("codes//cluster_analysis.R", local=TRUE)
if (run_summary_plotting) source("codes//summary_plotting.R", local=TRUE)

# Output by individual scripts is available in: ./output (all spreadsheets)
# and .//figures (all figures from postprocessing)