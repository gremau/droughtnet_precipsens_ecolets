# Data Guide

The data directory (referred to as DATADIR here) for this project contains both raw and derived data. Most of the raw data come from public sources, which are linked in the sections below. We have also tried to describe the derived data products in this document and with appropriate metadata files (in .yaml format). Where possible we describe which scripts produce the derived files, but see the [Code guide](code_guide.md) for more details.

## DATADIR/USHCN_raw - raw USHCN data

We used USHCN long-term climate data in the analysis. USHCN version 2.5 data files were downloaded from the NOAA National Climatic Data Center (NCDC) ftp server at <ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/v2.5>. Files for individual climate variables (precip, temperature, etc.) come as tar.gz archives.
* The USHCN station inventory was also used and is from the same ftp server (<ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/v2.5/ushcn-v2.5-stations.txt>).
* For the Maurer et al. _Ecology Letters_ article we downloaded data on 21 May 2017, and the original archive files are included in this directory and can be extracted and used to reproduce the study.

## DATADIR/USHCN_derived - processed and derived USHCN data

Raw USHCN data was processed and reduced in several ways. The resulting derived data files are in this directory. These files include:

* monthly_ushcn_allstations.csv - monthly data in csv format processed directly from the raw USHCN data.
* monthly_ushcn_spei_allstations.csv - monthly USHCN data with SPEI and other related metrics calculated and added as columns.
* annual_ushcn_spei_allstations.csv - annual climate metrics and SPEI indices.
* allyr_ushcn_calcs_allstations.csv - station info, data coverage metrics, and long-term trend calculations for 1218 USHCN stations
* USHCN_10km_polygons.shp - A shapefile with 10km square polygons centered around each USHCN station location, which defines the MODIS averaging area for each.

All tabular files above (not the shapefile) have a .yaml metadata file describing columns and origin of the file (which script, etc.). Scripts that produced these files will be found in [`../bin/`](../bin/), or see the [Code guide](code_guide.md).

## DATADIR/MODIS - derived MODIS data

This directory has several files with derived data from MODIS (for the USHCN sites)

* MODISoutput_ushcn.csv - 8-day modis NDVI averages for the 10km area around each USHCN site from 2000-2016, derived direcly from the raw dataset.
* annual_ushcn_NDVI.csv - annually calculated statistics for NDVI derived from the 8-day dataset above.
* annual_ushcn_eeMODIS - annual EVI, NPP, and GPP data for the 10km area around each USHCN site from 2000-2016.

All tabular files above have a .yaml metadata file describing columns and origin of the file (which script, etc.).

Raw NDVI data come from the [Spruce et al. 2015](https://doi.org/10.3334/ORNLDAAC/1299) dataset. To reproduce the derived data in this directory, 2000-2016 .ncdf4 files can be downloaded and stored on a server that has available storage capacity and processing time. The gridded data can then be processed to NDVI timeseries with the [`extract_modisNDVI_ushcn.R`](../bin/extract_modisNDVI_ushcn.R) script (`bin` directory, produces 'MODISoutput_ushcn.csv) and reduced to annual values when ['create_ushcn_files.R'](../bin/create_USHCN_files.R) is run.

Raw EVI, NPP, and GPP datasets have their own DOIs and are available on Google Earth Engine.

* EVI: <https://dx.doi.org/10.3334/ORNLDAAC/1299> (EE: MODIS/006/MOD13Q1)
* GPP and NPP: <https://doi.org/10.1002/rse2.74> (EE: UMT/NTSG/v2/MODIS/GPP and UMT/NTSG/v2/MODIS/NPP)

Code to produce the annual_ushcn_eeMODIS.csv derived file is found in the [`DGS_modis`](../bin/DGS_modis/) directory (in `bin`).

## DATADIR/sensitivity_analysis - data products from statistical analyses

All files here are outputs from statistical analysis of MODIS and climate data together. Generally they come from the jupyter notebooks in [`../jupyter/`](../jupyter/) and are well described there.

## Other

The DATADIR may also contain:

* The NEON domain shapefile (NEONDomains.zip) - See the jupyter notebooks in `../jupyter/` for examples of how to assign domains to USHCN station
* Huxman - raw and derived data for examining the production-precipitation sensitivity relationships at the sites from Huxman et al. 2004.
