# This file does data reduction and calculations for the DroughtNet UNM
# study (and subsequent paper in Ecology Letters).
# Raw data files come from the NOAA USHCN network and MODIS (multiple products)

# First get some functions we wrote
source('xhcn_data_reduction.R')
source('xhcn_calculate_values.R')

# *** Path to DATADIR ***
# This should point to your copy of the Droughtnet precip sensitivity data
# All other output paths will be relative to this.
dn_path <- '~/GD_gregmaurer/droughtnet_precipsens_data/'

# USHCN paths
# This directory is created upon extraction of the 3 archives
raw_ushcnPath <- paste0(dn_path, 'USHCN_raw/', 'ushcn.v2.5.5.20170521/')
# This file is found at ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/v2.5
ushcnInventory <- paste0(dn_path, 'USHCN_raw/',
                         'ushcn-v2.5-stations.txt')
# Where to put derived USHCN data
ushcnPath <- paste0(dn_path, 'USHCN_derived/')

# Make a monthly ushcn file
# Takes raw USHCN data and inventory and creates a monthly file
# named "monthly_ushcn_allstations.csv"
monthlyfname <- monthly_ushcn_file(raw_ushcnPath, ushcnInventory, ushcnPath)

# Calculate SPEI values for above monthly file and create station calcs file
fnames <- monthly_xhcn_calcs(monthlyfname, ushcnInventory, ushcnPath,
                             hcn_network='ushcn', name='allstations')
monthlySPEIfname <- fnames[[1]]
stnCalcsFname <- fnames[[2]]

# Calculate full century SPEI and SPEI CV trends with the monthly SPEI file
# above - these are added to the station calcs file
spei_trend_calcs(monthlySPEIfname, stnCalcsFname, ushcnPath,
                   hcn_network='ushcn', name='allstations')

# Calculate annual SPEI values for the above monthly SPEI file
annual_calc_file(monthlySPEIfname, ushcnPath,
                   hcn_network='ushcn', name='allstations')

# Set path for MODIS data
modisPath <- paste0(dn_path,'MODIS/')

# This file is in the droughtnet data directory,
# It is output from "extract_modisNDVI_ushcn.R"
rawNDVIfname <- paste0(modisPath, 'MODISoutput_ushcn.csv')

# Calculate annual NDVI and zNDVI values for each USHCN site file
annual_NDVI_file(rawNDVIfname, modisPath, hcn_network='ushcn')

# These files are output from google earth engine using the 
# included js scripts (bin/DGS_modis)
rawEVIfname <- '~/GD_gregmaurer/ee_output/modEVI_ushcn_2.csv'
rawNPPfname <- '~/GD_gregmaurer/ee_output/modNPP_ushcn_1.csv'
rawGPPfname <- '~/GD_gregmaurer/ee_output/modGPP_ushcn_1.csv'
# Make annual MODIS EVI, GPP, and NPP data file
annual_eeMODIS_file(rawEVIfname, rawNPPfname, rawGPPfname, modisPath,
		    hcn_network='ushcn')
