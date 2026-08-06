# Otego Creek Brook Trout Habitat Modeling
This repository holds files and code for modeling Brook Trout habitat in Otego Creek, NY.

# Directories
`GIS/` holds all georeferenced files necessary for running compiling networks and data, and running spatial stream network models. Directories containing files used in the analysis are described in the list below: 

* `LSN3/` Directory containing pre-processed landscape network (LSN) using the SSNbler package following An Introduction to ‘SSNbler’. The pre-processing included fine-scale correction of topologically validated stream flow network files available through The National Stream Internet Project and the tutorial available through the SSNbler GitHub repository.

* `LSN4/` A working LSN directory where further processing within SSN2_eels_covariates.R can read/write. Online file versions are directory placeholders and will be overwritten in that code currently.

* `mid_atlantic_preds/` 1-km spaced prediction points along the flow network subsetted from the Prediction Points dataset available through The National Stream Internet Project. Within this directory, preds_elevation was clipped to the upper Susquehanna River watershed using the aoi shapefile in upper_susquehanna/

* `sites/` Spatial datasets for American eel observations within the upper Susquehanna River. observed/obs_elevations.shp is used in SSN2_eels_covariates.R to model occurrence.

`HSI/` Habitat suitability index field data for analysis.

`Literature/` Directory for pertinent documents related to data collection or analysis.

`results/` Directory for figures generated from the analysis.

`otego_habitat_ssn2.R` The R script used to compile observed site data, prediction site data, and the LSN into a spatial stream network (SSN) object and conduct SSN regression analysis for American eel in the upper Susquehanna River watershed.

