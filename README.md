# DroughtNetUNM - ecolets-figshare branch

Files for the UNM/ASU DroughtNet group, mostly consisting of data management and analysis code. This group is studying ecosystem responses to climatic variability in the continental U.S. (mostly).

## In this repository:

* _bin/_ - Scripts for raw data extraction and processing (MODIS, USHCN). See [Code Guide](docs/code_guide.md)
* _docs/_ - documentation; also see the table of contents below
* _figs/_ - figures output by analysis scripts (jupyter notebooks, mostly)
* _jupyter/_ - Jupyter notebooks containing data analysis, statistics, figure creation, and associated R functions (common_stats_functions.R). Documentation of processing, analysis, and figure creation is in the notebooks themselves. Also see [Code Guide](docs/code_guide.md).

## Source and derived data

The raw data used in this study are mostly public datasets and can be downloaded as needed. We have a shared directory of the raw data, and derived data produced in data reduction and analysis (DATADIR). A subset of the raw and derived data are collected in the Figshare collection associated with this project:

* https://doi.org/10.6084/m9.figshare.c.4780313

See the Data Guide document below for more detailed information on locating source or derived data. **To point data processing and analysis scripts to a new DATADIR, edit the paths in [`create_ushcn_files.R`](../bin/create_ushcn_files.R) and [`common_stats_functions`](../jupyter/common_stats_functions.R).**


## Documentation

Table of contents:

* [Code Guide](docs/code_guide.md) - how to run various data processing and analysis scripts
* [Data Guide](docs/data_guide.md) - guide to data sources and derived data files in DATADIR
* [File metadata](docs/metadata.md) - guide to processed data file names, column names, etc
* [Statistics](docs/statistics.md) - explanation of statistical models and outputs 

## Associated publications

