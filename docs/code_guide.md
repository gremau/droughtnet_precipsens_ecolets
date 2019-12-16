# Code guide

Descriptions of raw data extraction and processing R code used in this DroughtNet project. All scripts will be found in the `bin` or `jupyter` directory of this repository. Data files derived from many of these scripts will be saved to the DATADIR. See the [Data Guide](data_guide.md) for details on source data (raw) and derived data.

**To point the important scripts here to a new DATADIR, edit the paths in [`create_ushcn_files.R`](../bin/create_ushcn_files.R) and [`common_stats_functions`](../jupyter/common_stats_functions.R).**

## USHCN/GHCN scripts

* create_ushcn_files.R - the master script for processing raw USHCN data to monthly and annual data, calculating SPEI and derived metrics, and processing raw NDVI data to match the USHCN data. Basically it makes all the files currently needed for data analysis/figures/writing.
* xhcn_data_reduction.R - A bunch of functions are contained here for processing USHCN and MODIS data. They are called from create_ushch_files.R 
* xhcn_calculate_values.R - Again, several functions for making calculations from met data.


## MODIS NDVI scripts

* extract_modisNDVI_ushcn.R - this is the main script that extracts NDVI timeseries for each USHCN site from the Spruce et al. MODIS dataset. It outputs 'MODISoutput_ushcn.csv'. Should be run on a server with R, with access to the 16 years of gridded MODIS data (see Data Guide).

## Google EarthEngine scripts for MOD17 and EVI

All present in [`DGS_modis`](../bin/DGS_modis/). Load these into Google EarthEngine, along with the USHCN_10km_polygons.shp file in DATADIR, and these scripts should extract long-term data from each of the MODIS data products.

## Data analysis in jupyter notebooks

These notebooks contain data analysis, statistics, and figure creation code, and associated R functions (common_stats_functions.R). Documentation of processing, analysis, and figure creation is in the notebooks themselves. The notebooks will render with full output at this GitHub repository, or at <https://nbviewer.jupyter.org/>.
