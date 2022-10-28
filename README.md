# Parameter estimation of catchment properties in the United Kingdom
by Piotr Morawiecki and Philippe H. Trinh, Departament of Mathematical Sciences, University of Bath (pwm27@bath.ac.uk)

:warning: **This repository is still under construction.**

## About the repository

This GitHub repository include codes allowing to extract and process physical parameters describing catchment properties in UK for use in rainfall-runoff modelling. It is based on available datasets shared by their developers - these need to be downloaded before running the codes. A single dataset `output/merge_summary.csv`, which includes values of extracted parameters assigned to each catchment. A detailed description of processing algorithms together with discussion of their results is available in paper "*On the development and analysis of mathematical models of catchments. Part 1. Parameter estimation of catchment properties*" [[1]](#1). Preprint is available inside the repository (`parameter_estimation_of_catchment_properties.pdf`).

Futher in this readme we describe:
* required input datasets (and how they should be added),
* the code structure and dependencies between them,
* parameters included in the output summary.

## Required input datasets

Below we present information about all input dataset included in `data` directory. 

**Important:** Due to copyrights the first datasets is **not** included in the repository. It need to be downloaded manually following instructions in *accessing data* row. However it is only required by `soil_processing.R` script - other scripts can be performed without downloading the dataset. 

### 3D Soil Hydraulic Database of Europe at 250 m resolution

Author(s) | Brigitta Tóth, Melanie Weynants, László Pásztor, Tomislav Hengl
--- | --- 
Description | Data include soil properties for the whole Europe. They include as hydraulic conductivity and parameters from Mualem-van Genuchten (MvG) model at seven different depths up to 2 meters. The data was obtained using European pedotransfer functions[[2]](#2)
Accessing data | The dataset need to be requested via European Soil Data Centre website. Follow this link https://esdac.jrc.ec.europa.eu/content/3d-soil-hydraulic-database-europe-1-km-and-250-m-resolution and fill in the request form.
Data location | Place the file in: `data/ESDAC_soil_conductivity`
Content | The dataset should include following rasters:<ul><li>`HCC_K0_sl*.tif` - hydraulic</li><li>`HCC_alp_sl*.tif` - MvG α parameter</li><li>`HCC_n_sl*.tif` - MvG n parameter</li><li>`HCC_thr_sl*.tif` - residual water content</li><li>`HCC_ths_sl*.tif` - saturated water content</li></ul>
Dataset size | 15.1 GB (it can be reduced by deleting unused files)

### OS VectorMap District

Author(s) | Ordnance Survey
--- | --- 
Description | OS VectorMap GIS data contains data for use in district level mapping, and includes spatial datasets describing buildings, road, railway and energy infrastratcures, woodlands, surface water bodies etc. In our framework we only use shapefiles describing surface water bodies.  
Accessing data | Data is available at https://osdatahub.os.uk/downloads/open/VectorMapDistrict, and was already uploaded to our repository under Open Government Licence for public sector information.Data location | `data/OS_VectorMap_District`
Content | Dataset is divided into subdirectories representing individual National Grid Reference regions. From the large original dataset we use only the following shapefiles (as long as they are defined for given region):<ul><li>`**_SurfaceWater_Area.shp` - incl. wide rivers, lakes, reserviors etc.</li><li>`**_SurfaceWater_Line.shp` - incl. small rivers, channels and streams</li><li>`**_TidalBoundary.shp` - incl. coastline</li></ul>
Dataset size | 1.68 GB

### OS Open Rivers

Author(s) | Ordnance Survey
--- | --- 
Description | OS Open Rivers GIS data contains over 144,000 km of water bodies and watercourses map data. These include freshwater rivers, tidal estuaries and canals.
Accessing data | Data is available at https://osdatahub.os.uk/downloads/open/OpenRivers, and was already uploaded to our repository under Open Government Licence for public sector information.
Data location | `data/OS_Open_Rivers`
Content | Original data set includes two shapefiles, from which only 'WatercourseLink.shp' is used in our framework
Dataset size | 244 MB

### OS Terrain 50


Author(s) | Ordnance Survey
--- | --- 
Description | OS Terrain 50 contains digital terrain model (DTM) data describing elevation over entire UK in resolution of 50m.
Accessing data | Data is available at https://osdatahub.os.uk/downloads/open/Terrain50, and was already uploaded to our repository under Open Government Licence for public sector information.
Data location | `data/OS_Terrain_50`
Content | Dataset is divided into subdirectories representing individual National Grid Reference regions. Each region is divided into 10x10 subtiles. DTM raster for each of them is included in a separate .asc format.
Dataset size | 576 MB

### National River Flow Archive (NRFA)

Author(s) | UK Centre for Ecology & Hydrology 
--- | ---
Description | NRFA consist a wide range of UK catchment descriptors (including mean rainfall, peak flows, land cover etc.) as well as time series of Gauged Daily Flow data (measured at the gauging stations) and Catchment Daily Rainfall [[3]](#3).
Accessing data | Data from NRFA do not need to be manually downloaded. It is automatically obtained by relevant R scripts using using NRFA API. The catchment boundaries are already included in `data/NRFA_catchment_boundaries` directory.

### UK3D

Author(s) | British Geological Survey 
--- | --- 
Description | UK3D is a national-scale network of intersecting cross-sections (also known as a *fence diagram*) of complex rocks and structures that make up the UK landmass [[4]](#4). Their properties include the produvitivity of the aquifer, which allows us to estimate its depth. 
Accessing data | Data is available at https://www.bgs.ac.uk/datasets/uk3d/, and was already uploaded to our repository under Open Government Licence for public sector information.
Data location | `data/BGS_UK3D`
Content | Original data set includes three shapefiles, from which only 'BGS_UK3D_v2015_3D_Cross_Sections.shp' is used in our framework
Dataset size | 57.8 MB

## Code structure

The structure of the code is presented in the figure below. The processing codes are divided into files based on the input dataset they are using, so that not all datasets need to be downloaded to rerun given part of the code. The datasets required by each script are listed at the beginning of a given script. Results of `streamline_extraction.R` were pregenerated, since it takes very long time to compute them.

All the processing codes generate summaries in .csv format, which are then merged by `merge_summary.R` code to a single file, `merge_summary.csv`. The postprocessing files allow to perform statistical analysis of this summary (including correlation, clustering and PCA) and generate maps showing distribution of parameter values over UK.

![Image cannot be uploaded](https://people.bath.ac.uk/pwm27/code_structure.svg)

## Output variables

Table below lists all output variables and their description.

parameter | symbol | unit | description | extraction code
--- | --- | --- | --- | --- 
id | - | - | catchment's id number as specified in the NRFA database; the detailed description of all catchments can be found in https://nrfaapps.ceh.ac.uk/nrfa/ws/station-info?station=%a&format=html&fields=all | - |
catchment_width_B | <img src="https://render.githubusercontent.com/render/math?math=L_x^\text{stream}"> | [m] | catchment width defined as an average length of a streamline | `streamline_processing.R`
... | ... | ... | ... | ... 
... | ... | ... | ... | ... 
... | ... | ... | ... | ... 
... | ... | ... | ... | ... 
... | ... | ... | ... | ... 

## Licence

## References
<a id="1">[1]</a> 
P. W. Morawiecki, P. H. Trinh. "On the development and analysis of mathematical models of catchments. Part 1. Parameter estimation of catchment properties" (2022)

<a id="2">[2]</a> 
T., Brigitta, et al. "3D soil hydraulic database of Europe at 250 m resolution." Hydrological Processes 31.14 (2017): 2662-2666.

<a id="3">[3]</a> 
M. J. Fry, O. Swain. "Hydrological data management systems within a national river flow archive." British Hydrological Society (2010).

<a id="4">[4]</a> 
C. N. Waters, et al. "The construction of a bedrock geology model for the UK: UK3D_v2015." British Geological Survey (2016)
