library('readr')
library('dplyr')
library('rgdal')
library('ggplot2')
library('sp')
library('forcats')

# Functions in this file are used to generate dataframes, transform datasets,
# and analyze data for the droughtnet project. They are sourced and used in
# several of the jupyter notebook analyses (ecosystem sensitivity, data,
# long-term climate stats, some periodic metadata checks, etc.)

# *** Path to DATADIR ***
# This should point to your copy of the Droughtnet precip sensitivity data
# All other output paths will be relative to this.
dn_path <- '~/GD_gregmaurer/droughtnet_precipsens_data/'

# Annual MODIS NDVI and zNDVI data, including full-year and growing-season
# means, sums, and other stats. One value per year per station since 2000.
annNDVI_fname <- paste0(dn_path, 'MODIS/annual_ushcn_NDVI.csv')
annNDVI_fname_Hux <- paste0(dn_path, 'Huxman/annual_huxman_NDVI.csv')

# Annual and water-year precip, temperature, SPEI, SPI, and other climate stats
# for study locations. One value per year per station since late 1800s
annSPEI_fname <- paste0(dn_path,
                        'USHCN_derived/annual_ushcn_spei_allstations.csv')
annSPEI_fname_Hux <- paste0(dn_path, 'Huxman/annual_huxman_spei.csv')

# Study location aggregate summary data including location, elevation, data 
# stats, (percent missing, etc) and climate trends (SPEI, SPI, T, P) for each
# station. 
USHCN_fname <- paste0(dn_path,
                      'USHCN_derived/allyr_ushcn_calcs_allstations.csv')
Huxman_fname <- paste0(dn_path, 'Huxman/allyr_huxman_calcs.csv')

# Neon domains shapefile path
neonDom_fname <- paste0(dn_path, 'NEONdomains/')

# Annual MODIS EVI, NPP and GPP data, including full-year and growing-season
# means, sums, and other stats. One value per year per station since 2000.
annMODIS_fname <- paste0(dn_path, 'MODIS/annual_ushcn_eeMODIS.csv')


sensitivity_stats_df <- function(locs='ushcn'){
    # This function opens the NDVI, SPEI/climate, and study location stats
    # data files on the server and merges them into one dataframe that
    # can be used for statistical analysis - namely calculating ecosystem
    # sensitivity statistics.
    #
    # It also transforms the data and adds things to the finished dataframe:
    #
    # 1. Assigns each study location to a NEON domain using a shapefile
    # 2. Trims out climate statistics earlier than the year 2000
    # 3. Removes a few rows that are not assigned to NEON domains (not sure why)
    #
    # The resulting dataframe is used in several jupyter notebooks
    if (locs=='ushcn') {
        annNDVI <- annNDVI_fname
        annSPEI <- annSPEI_fname
        annMODIS <- annMODIS_fname
        locs_fname <- USHCN_fname
    } else if (locs=='huxman') {
        annNDVI <- annNDVI_fname_Hux
        annSPEI <- annSPEI_fname_Hux
        locs_fname <- Huxman_fname
    }
    
    print(paste('Loading', annNDVI))
    annNDVI <- read_csv(annNDVI)
    print(paste('Loading', annSPEI))
    annSPEI <- read_csv(annSPEI)
    #print(paste('Loading', annMODIS))
    #annMODIS <- read_csv(annMODIS)
    print(paste('Loading', locs_fname))
    locCalcs <- read_csv(locs_fname)
    
    # Assign NEON domain names to each station on the calcs USHCN dataframe
    # First open shapefile and reproject to geographic WGS84
    print('Assigning NEON domain names to locCalcs dataframe...')
    neonDomains <- readOGR(neonDom_fname,"NEON_Domains")
    neonDomains <- spTransform(neonDomains, CRS("+proj=longlat +datum=WGS84"))
    # Exclude three domains outside the continental U.S.
    neonDomains<- subset(neonDomains, DomainID!=18 & DomainID!=19 & 
                         DomainID!=20)
    # Generate a points dataset from locCalcs
    ushcnpoints <- SpatialPoints(locCalcs[,c('longitude', 'latitude')],
                                 proj4string=CRS("+proj=longlat +datum=WGS84"))
    # over assigns a neon Domain name to each ushcn point based on which
    # polygon it is contained in
    pt.in.poly <- over(ushcnpoints, neonDomains, returnlist=TRUE)
    locCalcs$neonDomainName <- pt.in.poly$DomainName
    # This should remove unused factor levels (Tundra, Taiga, Pacific Tropical)
    locCalcs$neonDomainName <- factor(locCalcs$neonDomainName)
    
    if (locs=='ushcn') {
        print(paste('Loading', annMODIS))
        annMODIS <- read_csv(annMODIS)
        # Subset SPEI to year 2000+, then merge SPEI, NDVI, and USHCN 
        # dataframes into one
        print('Subsetting SPEI (>1999) and merging three dataframes...')
        annSPEIs <- subset(annSPEI, year > 1999)
        df <- merge(annSPEIs, annNDVI, by=c('stationid', 'year'))
        # Merge in MODIS, but add nans for year 2000 (starts 2001)
        df <- merge(df, annMODIS, by=c('stationid', 'year'), all.x=TRUE)
        df <- merge(df, locCalcs, by='stationid')
        } else if (locs=='huxman') {
        print('Subsetting SPEI (>1999) and merging three dataframes...')
        annSPEIs <- subset(annSPEI, year > 1999)
        df <- merge(annSPEIs, annNDVI, by=c('stationid', 'year'))
        df <- merge(df, locCalcs, by='stationid')
        }
    
    # For some reason some rows have no neon domain
    # and are assigned NA. Remove them
    rows_1 = dim(df)[1]
    df <- subset(df, !is.na(neonDomainName))
    print(paste('Removing', rows_1-dim(df)[1], 'rows missing domain names...'))
    
    return(df)
    }

longterm_climate_df <- function(){
    # This function opens the SPEI/climate and USHCN stats
    # data files on the server and merges them into one dataframe that
    # can be used for statistical analysis - namely calculating long-term
    # averages for climate, aridity indices, etc.
    #
    # It also transforms the data and adds things to the finished dataframe:
    #
    # 1. Assigns each USHCN station to a NEON domain using a shapefile
    # 2. Removes some rows that are not assigned to NEON domains (not sure why)
    #
    # The resulting dataframe is used in several jupyter notebooks
    print(paste('Loading', annSPEI_fname))
    annSPEI <- read_csv(annSPEI_fname)
    print(paste('Loading', USHCN_fname))
    locCalcs <- read_csv(USHCN_fname)
    
    # Assign NEON domain names to each station on the calcs USHCN dataframe
    # First open shapefile and reproject to geographic WGS84
    print('Assigning NEON domain names to locCalcs dataframe...')
    neonDomains <- readOGR(neonDom_fname,"NEON_Domains")
    neonDomains <- spTransform(neonDomains, CRS("+proj=longlat +datum=WGS84"))
    # Exclude three domains outside the continental U.S.
    neonDomains<- subset(neonDomains, DomainID!=18 & DomainID!=19 & 
                         DomainID!=20)
    # Generate a points dataset from locCalcs
    ushcnpoints <- SpatialPoints(locCalcs[,c('longitude', 'latitude')],
                                 proj4string=CRS("+proj=longlat +datum=WGS84"))
    # over assigns a neon Domain name to each ushcn point based on which
    # polygon it is contained in
    pt.in.poly <- over(ushcnpoints, neonDomains, returnlist=TRUE)
    locCalcs$neonDomainName <- pt.in.poly$DomainName
    # This should remove unused factor levels (Tundra, Taiga, Pacific Tropical)
    locCalcs$neonDomainName <- factor(locCalcs$neonDomainName)
    
    # Subset SPEI to year 2000+, then merge SPEI, NDVI, and USHCN 
    # dataframes into one
    print('Merging two dataframes...')
    df <- merge(annSPEI, locCalcs, by='stationid')
    
    # For some reason some rows have no neon domain
    # and are assigned NA. Remove them
    rows_1 = dim(df)[1]
    df <- subset(df, !is.na(neonDomainName))
    print(paste('Removing', rows_1-dim(df)[1], 'rows missing domain names...'))
    
    return(df)
    }

get_locCalcs_and_stats <- function(){
    # Load the locCalcs dataset and print some statistics about it
    print(paste('Loading', USHCN_fname))
    df <- read_csv(USHCN_fname)
    # First, number of sites and missing data stats
    nsites <- dim(df)[1]
    print(paste(nsites, 'USHCN sites total'))
    
    print(paste('Maximum percent of record missing (all stations) =',
                max(df$pctmissing)))
    print(paste('Maximum years of data missing (all stations) =',
                (max(df$pctmissing)/100)* max(df$nyears)))
    print(paste('Minimum number of years in record (all stations) =', 
                min(df$nyears)))
    # Now number of spei trends
    print(paste(sum(is.na(df$spei12mo_trend)),
                'sites have 12mo SPEI trends assigned NA'))
    print(paste(sum(df$spei12mo_trend_sig > 0.05, na.rm=T),
                'non-significant 12mo SPEI trends (alpha=0.05)'))
    print(paste(sum(df$spei12mo_trend_sig > (0.05/(nsites - 1)), na.rm=T),
                paste0('non-significant 12mo SPEI trends (alpha=0.05,',
                       'Bonferroni corrected)')))
    # Now number of spei CV trends
    print(paste(sum(is.na(df$spei12mo_cv5yr_trend)),
                'sites have 12mo SPEI CV trends assigned NA'))
    print(paste(sum(df$spei12mo_cv5yr_trend_sig > 0.05, na.rm=T),
                'non-significant 12mo SPEI CV trends (alpha=0.05)'))
    print(paste(sum(df$spei12mo_cv5yr_trend_sig > (0.05/(nsites-1)), na.rm=T),
                paste0('non-significant 12mo SPEI CV trends (alpha=0.05,',
                       'Bonferroni corrected)')))
    return(df)
    }   

get_annSPEI_and_stats <- function(){
    # Load the annSPEI dataset and print some statistics about it
    print(paste('Loading', annSPEI_fname))
    annSPEI <- read_csv(annSPEI_fname)
    print(paste('Earliest/latest year in USHCN (all stations) =',
                min(annSPEI$year), '/', max(annSPEI$year)))
    # Max year by station
    maxsta <- annSPEI %>% group_by(stationid) %>% summarise(Value = max(year))
    # Min year by station
    minsta <- annSPEI %>% group_by(stationid) %>% summarise(Value = min(year))
    print(paste('Latest station record starting year =', max(minsta$Value)))
    print(paste('Earliest station record ending year =', min(maxsta$Value)))
    print(paste('Shortest USHCN record is',
                min(maxsta$Value)-max(minsta$Value),'years'))
    return(annSPEI)
    }

rename_domains_for_plotting <- function(df){
    df <- df %>% mutate(neonDomainName = fct_recode(neonDomainName,
                "S Rockies/Col. Plat."="Southern Rockies / Colorado Plateau",
                "Appalachians/C. Plat." = "Appalachians / Cumberland Plateau"))
    return(df)
}

# Function for extracting fixed & random effects, and residuals from lme model
get_fe_2re <- function(model, modelname, random1, random2){
    # Get fixed and random1 effects (domain)
    model_fe <- data.frame(rbind(model[['coefficients']][['fixed']],
                      model[['coefficients']][['random']][[random1]]))
    colnames(model_fe) <- c(paste0(modelname, '_Intcpt'),
                            paste0(modelname, '_Slope'))
    model_fe[random1] <- c('FixedEffects',
                           row.names(model_fe)[2:length(row.names(model_fe))])
    # Get random2 effects (station)
    model_re <- data.frame(model[['coefficients']][['random']][[random2]])
    colnames(model_re) <- c(paste0(modelname, '_randomIntcpt'),
                            paste0(modelname, '_randomSlope'))
    # Get domain and station names from the random effects
    model_re[random2] <- paste0('U',do.call(
        'rbind', strsplit(as.character(rownames(model_re)),
                          '/U',fixed=TRUE))[,2])
    # Calculate and merge in the rmse for each station
    # There are 3 different residuals calculated for lme, not sure what is 
    # best. see here:
    # https://www.rdocumentation.org/packages/nlme/versions/3.1-141/topics/residuals.lme
    rmse <- function(x) sqrt(mean(x^2))
    # This is an "unbiased" rmse dividing by df
    rmse.ub <- function(x) sqrt(sum(x^2)/(length(x)-3))
    rmse_raw <- lapply(residuals(model, level=2, asList=T), rmse)
    #rmse_pear <- lapply(residuals(model, level=2, asList=T, type='pearson'),
    #                    rmse)
    #rmse_norm <- lapply(residuals(model, level=2, asList=T, type='normalized'),
    #                    rmse)
    rmse_ub <- lapply(residuals(model, level=2, asList=T), rmse.ub)
    #rmse_ubnorm <- lapply(residuals(model,level=2, asList=T,type='normalized'),
    #                      rmse.ub)
    # Make and rename a dataframe
    rmses <- data.frame(rmse_raw=unlist(rmse_raw),#rmse_pear=unlist(rmse_pear),
                        #rmse_norm=unlist(rmse_norm),
                        rmse_ub=unlist(rmse_ub))
                        #rmse_ubnorm=unlist(rmse_ubnorm))
    colnames(rmses) <- c(paste0(modelname, '_rmse_raw'),
                         #paste0(modelname, '_rmse_pear'),
                         #paste0(modelname, '_rmse_norm'),
                         paste0(modelname, '_rmse_ub'))
                         #paste0(modelname, '_rmse_ubnorm'))
    rmses[random2] <- paste0('U', do.call(
        'rbind',strsplit(as.character(rownames(rmses)),
                         '/U',fixed=TRUE))[,2])
    # To get the normalized RMSE, first get y values
    yval <- fitted(model, level=2) + residuals(model, level=2)
    # This gives a named vector, the next creates a named list
    yval.list <- split(unname(yval),names(yval))
    # Calculate the range in y values
    yval.range <- function(x) max(x) - min(x)
    yval_range <- lapply(yval.list, yval.range)
    ranges <- data.frame(yval_range=unlist(yval_range))
    colnames(ranges) <- c(paste0(modelname, '_yval_range'))
    ranges[random2] <- paste0('U',do.call(
        'rbind',strsplit(as.character(rownames(ranges)),
                         '/U',fixed=TRUE))[,2])
    # Now merge ranges and rmses into one dataframe
    model_re <- merge(model_re, ranges, by='stationid', all.x=T)
    model_re <-  merge(model_re, rmses, by='stationid', all.x=T)
    # Calculate the normalized RMSE (using unbiased)
    model_re[paste0(modelname,'_nrmse_ub')] <- (
        model_re[paste0(modelname,'_rmse_ub')]/model_re[paste0(modelname,
                                                               '_yval_range')])
    return(list(model_fe, model_re))
}


#### Mapping functions ####


## Create some background layers

# US map boundary
map_base <- map_data('usa')# change to states for state boundaries
map_base_aes <- aes(x=long, y=lat, group=group)
#map_base_poly <- geom_polygon(data=map_base, map_base_aes)
us_border <- geom_polygon(data=map_base, map_base_aes, color="black",
                          fill=NA, size=0.25)

#Load shapefile of neon domains
neonDomains <- readOGR(dsn=paste0(dn_path, 'NEONdomains/'),
                       layer="NEON_Domains")
shp <- spTransform(neonDomains, CRS("+proj=longlat +datum=WGS84"))
domains <- fortify(shp)
# Remove Hawaii and Alaska
domains$id <- as.numeric(domains$id) -1
domains_poly <- subset(domains,
                       !is.na(id) & id!=18 & id!=19 & id!=20
                       & lat<50 & lat>20)


# adapted from: http://www.howtobuildsoftware.com/index.php/how-do/bniC/r-ggplot2-geospatial-interpolation-smoothing-out-ggplot2-map
krige_area <- function(df,datvar,latvar,longvar,map_base_data=map_data('usa')){
    library(data.table)
    library(ggplot2)
    library(automap)
    library(plyr)
    # Data munging
    coord_vars <- c(latvar,longvar)
    df_sub <- df[,c(latvar,longvar,datvar)]
    data_vars <- setdiff(colnames(df_sub), coord_vars)
    sp_points <- SpatialPoints(df_sub[,coord_vars])
    sp_df <- SpatialPointsDataFrame(sp_points, df_sub[,data_vars,drop=FALSE])

    # Create a fine grid
    pixels_per_side <- 200
    bottom.left <- apply(sp_points@coords,2,min)
    top.right <- apply(sp_points@coords,2,max)
    margin <- abs((top.right-bottom.left))/10
    bottom.left <- bottom.left-margin
    top.right <- top.right+margin
    pixel.size <- abs(top.right-bottom.left)/pixels_per_side
    g <- GridTopology(cellcentre.offset=bottom.left,
             cellsize=pixel.size,
             cells.dim=c(pixels_per_side,pixels_per_side))
    # Rename columns in map base data to match latvar and longvar
    colnames(map_base_data)[match(c("long","lat"),
                                  colnames(map_base_data))] <- c(longvar,latvar)
    # Loop through each region in map data and make a polygon list and fit
    # points within the polygons
    polys <- function(x) {
      pg <- unique(x$region)
      print(pg)
      Polygons(list(Polygon(x[,c(latvar,longvar)])),ID=pg)
    }
    mbd_pg <- SpatialPolygons(dlply(map_base_data, .(region), polys))
    grid_points <- SpatialPoints(g)
    in_points <- !is.na(over(grid_points,mbd_pg))
    fit_points <- SpatialPoints(as.data.frame(grid_points)[in_points,])
    head(sp_df)
    # Do kriging
    krig <- autoKrige(eval(parse(text=datvar))~1, sp_df, new_data=fit_points)
    interp_data <- as.data.frame(krig$krige_output)
    colnames(interp_data) <- c(latvar,longvar, paste(datvar, '_pred', sep=''),
                           paste(datvar, '_var', sep=''),
                               paste(datvar, '_stdev', sep=''))

    return(interp_data)
}



# Map an interpolated surface
plot_interp <- function(df_interp, pred_var_df, scale_expr, nonan_points,
                        lowCol='#003366',midCol='white',highCol='red',
                        midpt=0, nbin=20, plot_contours=T, plot_stations=T){
    if(!exists('us_border')){
        # Create borders for map
        map_base <- map_data('usa')# change to states for state boundaries
        map_base_aes <- aes(x=long, y=lat, group=group)
        us_border <- geom_polygon(data=map_base, map_base_aes,
                                  color="black", fill=NA, size=0.25)
    }
    if(!exists('domains_poly')){
        #Load shapefile of neon domains
        neonDomains <- readOGR(dsn=paste0(dn_path, 'NEONdomains/'),
                               layer="NEON_Domains")
        shp <- spTransform(neonDomains, CRS("+proj=longlat +datum=WGS84"))
        domains <- fortify(shp)
        # Remove Hawaii and Alaska
        domains$id <- as.numeric(domains$id) -1
        domains_poly <- subset(domains,
                               !is.na(id) & id!=18 & id!=19 & id!=20 & 
                               lat<50 & lat>20)
    }
    # Now create the figure
    fig1 <- ggplot(data=df_interp, aes(x=longitude, y=latitude)) + 
      geom_tile(aes_string(fill=pred_var_df),color=NA) +
      scale_fill_gradient2(low=lowCol,mid=midCol,high=highCol,
                           # name='GPPint\nSensitivity\nto Precip.',
                           midpoint=midpt) + 
      #mean(interp_NdviPrcp_sens$ndvi_prcp_randomSlope2_pred)) +
      us_border + coord_map(projection='mercator') +
      labs(fill = scale_expr)
    # Optional things to plot
    if(plot_contours){
        fig1 <- fig1 + stat_contour(aes_string(z=pred_var_df), bins=nbin,
                                    color="#999999")}
    if(plot_stations){
        fig1 <- fig1 + geom_point(data=nonan_points,color="black",
                                  alpha=0.5, size=0.1)}
    
    fig1 <- fig1 +
      # Comment out to remove neon boundaries
      geom_polygon(aes(x=long, y=lat, group=id), data=domains_poly,
                   colour='black', fill='black', alpha=0, size=.25) +
      ylab('Latitude') + xlab('Longitude') +
      theme_minimal() + theme(legend.position = c(0.9, 0.30))

    return(fig1)
}
