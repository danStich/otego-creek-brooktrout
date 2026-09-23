# Otego Creek Brook Trout Habitat Modeling
This repository holds files and code for modeling Brook Trout (*Salvelinus fontinalis*) habitat in Otego Creek, NY. Data for this project were collected by Braeden Victory, Caitlin Brislin, and Dan Stich through the SUNY Oneonta Biological Field Station summer internship program in summer 2026. All work was done in coordination with New York State Department of Environmental Conservation. 

# Directories
`GIS/` holds all georeferenced files necessary for running compiling networks and data, and running spatial stream network models. Directories containing files used in the analysis are described in the list below: 

* `LSN3/` Directory containing pre-processed landscape network (LSN) using the [`SSNbler` package](https://pet221.github.io/SSNbler/index.html) following [An Introduction to ‘SSNbler’](https://pet221.github.io/SSNbler/articles/introduction.html). The pre-processing included fine-scale correction of topologically validated stream flow network files available through [The National Stream Internet Project](https://research.fs.usda.gov/rmrs/projects/national-stream-internet#download-data) and the tutorial available through the [SSNbler GitHub repository](https://github.com/pet221/SSNbler/tree/main/inst/tutorials).

* `LSN4/` A working LSN directory where further processing within `otego_habitat_ssn2.R` can read/write. Online file versions are directory placeholders and will be overwritten when the code in that file is run.

* `mid_atlantic_preds/` 1-km spaced prediction points along the flow network subsetted from the Prediction Points dataset available through [The National Stream Internet Project](https://research.fs.usda.gov/rmrs/projects/national-stream-internet#download-data). Within this directory, `otego_preds.*` shapefiles are NSIP points clipped to the Otego Creek, NY watershed. Elevations were sampled from USGS digital elevation products in `otego_preds_elevation.*` shape files. The `otego_preds_muahaha.*` shape files are working fine-scale prediction points added separately from the NSIP products for finer resolution predictions in the watershed.

* `sites/` Spatial datasets for Brook Trout habitat sampling sites and HSI variables created from the data file in `HSI/`.

* **other spatial files** additional subdirectories contain cartographic, municipal, and other GIS files used for study planning, execution, and creation of study site maps.

`HSI/` holds Brook Trout habitat suitability index (HSI) data and derived calculations for Otego Creek, NY during summer 2026.

`Literature/` Directory for pertinent documents related to data collection or analysis.

`results/` Directory for figures generated from the analysis.

`otego_habitat_ssn2.R` The R script used to compile observed site data, prediction site data, and the LSN into a spatial stream network (SSN) object and conduct SSN regression analysis for American eel in the upper Susquehanna River watershed.

