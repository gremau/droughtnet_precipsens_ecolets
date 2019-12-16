library('readr')
library('dplyr')
library('lubridate')
library('reshape2')


ushcn_wide_to_long <- function(df, varname){
    # Take wide dataframe, add id to each row, melt 12 month rows to 1 month,
    # sort by id, then remove the id
    df_melt <- df %>% mutate(id=seq_len(n())) %>% 
        melt(id.vars=c('stationid', 'year', 'id'),
             measure.vars=c('jan', 'feb', 'mar', 'apr', 'may', 'jun',
                            'jul', 'aug', 'sep', 'oct', 'nov', 'dec')) %>%
    arrange(stationid, id, variable) %>% select(-id)
    colnames(df_melt) <- c('stationid', 'year', 'month', varname)
    dates <- as.Date(paste(as.character(df_melt$year),'-',
                           as.character(df_melt$month), '-01',
                           sep=''), format='%Y-%b-%d')
    df_melt$date <- dates + (days_in_month(dates) - 1)
    df_melt <- select(df_melt, -year, -month)
    df_melt <- df_melt[c('stationid','date', varname)]
    #print(head(df_melt))
    return(df_melt)
}


load_ushcn_inventory <- function(fname){
    # Helper function to load USHCN station inventory file
    print(paste('Loading ushcn inventory:', fname, '...'))
    widths <- c(11, 1, 8, 1, 9, 2, 5, 1, 2, 1, 30, 1, 6, 1, 6, 1, 6, 1, 2)
    cnames <- c('stationid', 'd1', 'latitude', 'd2', 'longitude', 'd3',
                'elev', 'd4', 'state', 'd5', 'name', 'd6',
                'flag1', 'd7', 'flag2', 'd8', 'flag3', 'd9',
                'utcoffset')
    inv <- read_fwf(fname, fwf_widths(widths,col_names=cnames))
    inv <- inv[c('stationid', 'latitude', 'longitude', 'elev', 'state',
                 'name', 'flag1', 'flag2', 'flag3', 'utcoffset')]
    print('Done.')
    return(inv)
}


monthly_ushcn_file <- function(ushcn_path, inventoryfile, outpath){
    # Load inventory
    inv <- load_ushcn_inventory(inventoryfile)

    # Get list of sites
    sites <- unique(inv$stationid)
    sites_good <- sites
    print(paste(length(sites_good), 'sites from inventory will be read.' ))

    widths <- c(11, 1, 4, rep(c(6, 1, 1, 1), times=12))
    cnames <- c('stationid', 'discard','year',
                'jan', 'dmflag1', 'qcflag1', 'dsflag1',
                'feb', 'dmflag2', 'qcflag2', 'dsflag2',
                'mar', 'dmflag3', 'qcflag3', 'dsflag3',
                'apr', 'dmflag4', 'qcflag4', 'dsflag4',
                'may', 'dmflag5', 'qcflag5', 'dsflag5',
                'jun', 'dmflag6', 'qcflag6', 'dsflag6',
                'jul', 'dmflag7', 'qcflag7', 'dsflag7',
                'aug', 'dmflag8', 'qcflag8', 'dsflag8',
                'sep', 'dmflag9', 'qcflag9', 'dsflag9',
                'oct', 'dmflag10', 'qcflag10', 'dsflag10',
                'nov', 'dmflag11', 'qcflag11', 'dsflag11',
                'dec', 'dmflag12', 'qcflag12', 'dsflag12')
    
    count <- 1
    for (s in sites_good){
        print(paste('Reading data for site', s))
        prcp_file <- paste(ushcn_path, s, '.FLs.52j.prcp', sep='')
        tmax_file <- paste(ushcn_path, s, '.FLs.52j.tmax', sep='')
        tmin_file <- paste(ushcn_path, s, '.FLs.52j.tmin', sep='')
        tavg_file <- paste(ushcn_path, s, '.FLs.52j.tavg', sep='')
        # Read and add columns, convert missing values to NA
        prcp <- read_fwf(file=prcp_file, na=c("","NA","-9999"),
                        fwf_widths(widths, col_names=cnames))
        tmax <- read_fwf(file=tmax_file, na=c("","NA","-9999"),
                        fwf_widths(widths, col_names=cnames))
        tmin <- read_fwf(file=tmin_file, na=c("","NA","-9999"),
                        fwf_widths(widths, col_names=cnames))
        tavg <- read_fwf(file=tavg_file, na=c("","NA","-9999"),
                        fwf_widths(widths, col_names=cnames))
        # Wide to long format, add dates, etc 
        prcp_m <- ushcn_wide_to_long(prcp, 'prcp')
        tavg_m <- ushcn_wide_to_long(tavg, 'tavg')
        tmax_m <- ushcn_wide_to_long(tmax, 'tmax')
        tmin_m <- ushcn_wide_to_long(tmin, 'tmin')
        
        joined <- left_join(prcp_m, tavg_m, by=c('stationid', 'date')) %>%
                left_join(., tmax_m, by=c('stationid', 'date')) %>%
                left_join(., tmin_m, by=c('stationid', 'date'))

        # Convert units
        joined$prcp <- joined$prcp/10
        joined$tavg <- joined$tavg/100
        joined$tmax <- joined$tmax/100
        joined$tmin <- joined$tmin/100
         
        #print(head(joined))
        if (count==1){
            ushcn_df <- joined
        } else {
            ushcn_df <- rbind(ushcn_df, joined)
        }
        count <- count + 1
    }
    outfile <- paste(outpath,  'monthly_ushcn_allstations.csv', sep='')
    print(paste('Writing monthly output file:', outfile))
    write_csv(ushcn_df, outfile)
    return(outfile)
}

ghcnd_to_monthly_file <- function(datafile, coveragefile, outpath, name,
                                  return_df=FALSE){
    # Function to take daily GHCN data (downloaded with rnoaa), impose a little
    # filtering by coverage, add tavg where missing, and then reduce the data
    # to monthly frequency.
    library('xts')
    # Load data
    print(paste('Loading', datafile, '...'))
    df <- read_csv(datafile, na=c('NA','','-9999'),
                   col_types=cols(id = col_character(),
                                  date = col_date(format = ""),
                                  prcp = col_integer(),
                                  tavg = col_integer(),
                                  tmax = col_integer(),
                                  tmin = col_integer()))
    print('Done.')
    # Load coverage file
    print(paste('Loading', coveragefile, '...'))
    df.cov <- read_csv(coveragefile)
    print('Done.')
    
    # Does datafile contain tavg?
    tavgflag <- 'tavg'%in%colnames(df)

    # Get list of sites in datafile and coverage file
    sites_df <- unique(df$id)
    sites_cov <- df.cov$id
    print(paste('There are ', length(sites_df), 'sites in daily datafile.' ))
    ### FILTER sites with coverage less than 90%
    prcp.lowcov <- df.cov$prcp < 0.9
    tmax.lowcov <- df.cov$tmax < 0.9
    tmin.lowcov <- df.cov$tmin < 0.9
    # Coverage logic varies depending on if tavg is present
    if (tavgflag){
        maxminlow <- tmax.lowcov | tmin.lowcov
        lowcoverage <- prcp.lowcov | (maxminlow & df.cov$tavg < 0.9)
    } else {
        lowcoverage <- prcp.lowcov | tmax.lowcov | tmin.lowcov
    }
    sites_goodcov <- sites_cov[!lowcoverage]
    print(paste('There are ', length(sites_goodcov), '/', length(sites_cov),
                'with coverage > 90%.'))
    sites_sel <- sites_df[sites_df %in% sites_goodcov]
    print(paste(length(sites_sel),
                'of these will be pulled from daily datafile'))

    count <- 1
    for (s in sites_sel){
        # Subset the data by site s
        df.site <- subset(df, id==s)
        print(paste('Daily-monthly conversion for site ', s))
        # Reduce to monthly, note that NA values are removed without knowing
        # how many are missing - potential source of bias, but only for T
        df.site.xts <- xts(df.site[,3:ncol(df.site)],
                           order.by=as.Date(df.site$date, "%Y-%m-%d"))
        df.site.prcp <- apply.monthly(df.site.xts$prcp, sum)
        df.site.tmax <- apply.monthly(df.site.xts$tmax, mean, na.rm=TRUE)
        df.site.tmin <- apply.monthly(df.site.xts$tmin, mean, na.rm=TRUE)
        # If we have a tavg column use that and backfill NA values with the
        # average of tmax/tmin
        if (tavgflag){
            df.site.tavg <- apply.monthly(df.site.xts$tavg, mean, na.rm=TRUE)
            df.site.mmtavg <- (df.site.tmax + df.site.tmin)/2
            tavgNA <- is.na(df.site.tavg)
            df.site.tavg[tavgNA] <- df.site.mmtavg[tavgNA]
        } else {
            df.site.tavg <- (df.site.tmax + df.site.tmin)/2
        }        
        colnames(df.site.tavg) <- 'tavg'
        df.site.mon <- merge(df.site.prcp, df.site.tavg, df.site.tmax,
                             df.site.tmin)
        newdf <- data.frame(date=index(df.site.mon), coredata(df.site.mon))
        newdf$stationid <- s
        newdf <- newdf[c('stationid','date','prcp','tavg','tmax','tmin')]
        # Convert units
        newdf$prcp <- newdf$prcp/10
        newdf$tavg <- newdf$tavg/10
        newdf$tmax <- newdf$tmax/10
        newdf$tmin <- newdf$tmin/10
        if (count < 2){
            df.site.mon.rbind <- newdf
        } else {
            df.site.mon.rbind <- rbind(df.site.mon.rbind, newdf)
        }
        count <- count + 1   
    }
    print(paste('Writing output file ', outpath,
                'monthly_ghcn_', name, '.csv', sep=''))
    write_csv(df.site.mon.rbind, paste(outpath, 'monthly_ghcn_',
                                       name, '.csv', sep=''))
    print('Done.')
    if (return_df){
        return(df.site.mon.rbind)
    }
}


ghcnd_to_monthly_teamloop <- function(concat_outpath){
    # Loop through all the team's files they downloaded from noaa,
    # filter sites, calculate monthly values (with ghcnd_to_monthly), 
    # and write a concatenated file for all sites
    datafiles <- c('~/sftp/data/GHCN/scott/droughtnet_GHCN_data_scott.csv',
                   '~/sftp/data/GHCN/alesia/droughtnet_GHCN_data_alesia.csv',
                   '~/sftp/data/GHCN/renee/droughtnet_GHCN_data_renee.csv',
                   '~/sftp/data/GHCN/greg/droughtnet_GHCN_data_greg_1.csv',
                   '~/sftp/data/GHCN/greg/droughtnet_GHCN_data_greg_2.csv',
                   '~/sftp/data/GHCN/greg/droughtnet_GHCN_data_greg_3.csv')
    covfiles <- c('~/sftp/data/GHCN/scott/droughtnet_GHCN_coverage_scott.csv',
                  '~/sftp/data/GHCN/alesia/droughtnet_GHCN_coverage_alesia.csv',
                  '~/sftp/data/GHCN/renee/droughtnet_GHCN_coverage_renee.csv',
                  '~/sftp/data/GHCN/greg/droughtnet_GHCN_coverage_greg_1.csv',
                  '~/sftp/data/GHCN/greg/droughtnet_GHCN_coverage_greg_2.csv',
                  '~/sftp/data/GHCN/greg/droughtnet_GHCN_coverage_greg_3.csv')
    outpath <- c('~/sftp/data/GHCN/scott/', '~/sftp/data/GHCN/alesia/',
                 '~/sftp/data/GHCN/renee/', '~/sftp/data/GHCN/greg/',
                 '~/sftp/data/GHCN/greg/', '~/sftp/data/GHCN/greg/')
    names <- c('scott', 'alesia', 'renee', 'greg_1', 'greg_2', 'greg_3')
    #concat_outpath <- outpath

    for (i in 1:length(datafiles)){
        dfi <- ghcnd_to_monthly_file(datafiles[i], covfiles[i], outpath[i],
                                     names[i], return_df=TRUE)
        if (i==1){dfc <- dfi} else {dfc <- rbind(dfc, dfi)}
    }
    write_csv(dfc, paste(concat_outpath, 'monthly_ghcn_selstations.csv',
                         sep=''))
}

monthly_huxman_file <- function(huxclim_path, inventoryfile, outpath){
    # Load inventory
    inv <- read_csv(inventoryfile,
		    col_types = cols(inconus = col_logical()))

    # Get list of sites in CONUS
    sites <- inv$stationid[inv$inconus]
    sites_good <- sites
    print(paste(length(sites_good), 'sites from inventory will be read.' ))
    
    count <- 1
    for (s in sites_good){
	prismfile <- paste(huxclim_path,
		'PRISM_ppt_tmin_tmean_tmax_provisional_4km_190001_201901_',
		s, '.csv', sep='')
    	prism <- read_csv(prismfile, skip=10)
	# Clean columns and dates
	#names(prism) <- sub(" \\(.+\\)","", names(prism))
	names(prism) <- c('date','prcp','tmin','tavg','tmax')
	prism$stationid <- s
	#prism$Date2 <- as.Date(as.yearmon(prism$Date), frac = 1)
	prism$date <- ymd(paste0(prism$date, '-01')) + months(1) - days(1)

        print(paste('Reading PRISM data for site', s))        
        print(head(prism))

        if (count==1){
            hux_df <- prism
        } else {
            hux_df <- rbind(hux_df, prism)
        }
        count <- count + 1
    }
    outfile <- paste(outpath,  'monthly_climate_huxman.csv', sep='')
    print(paste('Writing monthly output file:', outfile))
    write_csv(hux_df, outfile)
}

annual_calc_file <- function(monthly_calcfile, outpath, hcn_network=NULL,
                        name=NULL){
    if (is.null(hcn_network)){stop('Please enter HCN network')}
    if (is.null(name)){stop('Please enter a file name string')}
    # Load file
    print(paste('Loading', monthly_calcfile, '...'))
    df <- read_csv(monthly_calcfile,
                   col_types=cols(stationid = col_character(),
                                  date = col_date(format = ""),
                                  prcp = col_double(),
                                  tavg = col_double()
                                  #tmax = col_double(),
                                  #tmin = col_double()
                                  ))
    print('Done.')
    #modify name to keep file naming convention
    if (name!=''){name <- paste0('_', name)}
    # Change date to POSIXct
    df$dates <- with(df,as.POSIXct(date,format="%Y-%m-%d"))
    
    # Make new variables, year and month and wyear
    df <- mutate(df,
           month=as.numeric(format(dates,"%m")),
           year=as.numeric(format(dates,"%Y")),
           wyear=year)
    # Create a water year column (add 1 to Oct-Dec)
    df$wyear[df$month > 9] <- df$year[df$month > 9] + 1

    # Use ddply to reduce to annual data
    ## According to year
    by_stationyear <- group_by(df, stationid, year)
    df_ann <- summarise(by_stationyear,
                        tavg_mean=mean(tavg, na.rm=T),
                        #meantmax=mean(tmax, na.rm=T),
                        #meantmin=mean(tmin, na.rm=T),
                        prcp_sum=sum(prcp),
                        pet_mean=mean(pet, na.rm=T),
                        spei3mo_mean=mean(spei3mo, na.rm=T),
                        spei6mo_mean=mean(spei6mo, na.rm=T),
                        spei9mo_mean=mean(spei9mo, na.rm=T),
                        spei12mo_mean=mean(spei12mo, na.rm=T),
                        spei18mo_mean=mean(spei18mo, na.rm=T),
                        spei24mo_mean=mean(spei24mo, na.rm=T),
                        spei12mo_december=last(spei12mo),
			spi12mo_mean=mean(spi12mo, na.rm=T),
			spi12mo_december=last(spi12mo))
    #print(head(df_ann))
    ## According to water year
    by_stationWyear <- group_by(df, stationid, wyear)
    df_WYann <- summarise(by_stationWyear,
                          prcp_wysum=sum(prcp),
                          tavg_wymean=mean(tavg, na.rm=T),
                          #wymeantmax=mean(tmax, na.rm=T),
                          #wymeantmin=mean(tmin, na.rm=T),
                          pet_wymean=mean(pet, na.rm=T),
                          spei3mo_wymean=mean(spei3mo, na.rm=T),
                          spei6mo_wymean=mean(spei6mo, na.rm=T),
                          spei9mo_wymean=mean(spei9mo, na.rm=T),
                          spei12mo_wymean=mean(spei12mo, na.rm=T),
                          spei18mo_wymean=mean(spei18mo, na.rm=T),
                          spei24mo_wymean=mean(spei24mo, na.rm=T),
                          spei6mo_september=last(spei6mo),
                          spei12mo_september=last(spei12mo),
			  spi12mo_wymean=mean(spi12mo, na.rm=T),
			  spi12mo_september=last(spi12mo))

    df_ann <- merge(df_ann, df_WYann, by.x=c("stationid", "year"),
                    by.y=c("stationid", "wyear"))
    
    out_fname <- paste(outpath, 'annual_', hcn_network, '_spei',
                       name, '.csv', sep='')
    print(paste('Writing', out_fname, '...'))
    write_csv(df_ann, out_fname)
    print('Done.')
}

annual_NDVI_file <- function(raw_NDVIfile, outpath,
                             ndviCol=quo(mean.MODIS.NDVI),
                             doyCol=quo(MODIS.end.DOY),
                             hcn_network=NULL){
    if (is.null(hcn_network)){stop('Please enter HCN network')}
    # Load file
    print(paste('Loading', raw_NDVIfile, '...'))
    rawNDVI <- read_csv(raw_NDVIfile)
    print('Done.')

    # Get date and month column
    rawNDVI$date <- as.Date(rawNDVI[[quo_name(doyCol)]] + 4,
                            origin=as.Date('1999-12-31'))
    rawNDVI$month <- as.numeric(format(rawNDVI$date, '%m'))

    # Create a monthly dataframe
    # Have to do funny symbol unquoting for the column variables
    rawNDVI_m <- rawNDVI %>% group_by(stationid, year, month) %>% 
        summarise(ndvi_mean=mean(!!ndviCol),ndvi_max=max(!!ndviCol),
                  ndvi_min=min(!!ndviCol))

    # Calculate z score of monthly mean NDVI for each site
    sites <- unique(rawNDVI_m$stationid)

    for (s in sites){
        siteNDVI <- subset(rawNDVI_m, stationid==s)
        siteNDVI_m <- mean(siteNDVI$ndvi_mean, na.rm=T)
        siteNDVI_sd <- sd(siteNDVI$ndvi_mean, na.rm=T)
        siteNDVI$zNDVI<- (siteNDVI$ndvi_mean - siteNDVI_m)/siteNDVI_sd
        if (s==sites[1]){
            newNDVI_m <- siteNDVI
        }
        else{
            newNDVI_m <- rbind(newNDVI_m, siteNDVI)
        }
    }

    # Create annual dataframe from the monthly one
    newNDVI_ann <- newNDVI_m %>% group_by(stationid, year) %>% 
        summarise(ndvi_amean=mean(ndvi_mean), ndvi_amax=max(ndvi_max),
                  ndvi_amin=min(ndvi_min), ndvi_sum=sum(ndvi_mean),
                  zndvi_mean=mean(zNDVI), zndvi_max=max(zNDVI),
                  zndvi_min=min(zNDVI),zndvi_sum=sum(zNDVI))

    newNDVI_ann <- rename(newNDVI_ann, ndvi_mean=ndvi_amean,
                          ndvi_max=ndvi_amax, ndvi_min=ndvi_amin)

    # Create a growing season dataframe from the monthly one
    newNDVI_gs <- newNDVI_m %>% group_by(stationid, year) %>% 
        filter(month>=4 & month<=9) %>%
        summarise(ndvi_gsmean=mean(ndvi_mean), ndvi_gssum=sum(ndvi_mean),
                  zndvi_gsmean=mean(zNDVI), zndvi_gssum=sum(zNDVI))

    # Merge annual and growing season dataframes
    newNDVI_ann <- merge(newNDVI_ann, newNDVI_gs, by=c('stationid', 'year'))
    
    # Write file
    out_fname <- paste0(outpath, 'annual_', hcn_network,
                        '_NDVI', '.csv')
    print(paste('Writing', out_fname, '...'))
    write_csv(newNDVI_ann, out_fname)
    print('Done.')
}

annual_eeMODIS_file <- function(raw_EVIfile, raw_NPPfile, raw_GPPfile, 
				    outpath, eviCol=quo(mean),
				    nppCol=quo(NPP), gppCol=quo(GPP),
				    dateCol=quo(date),
				    hcn_network=NULL){
    if (is.null(hcn_network)){stop('Please enter HCN network')}
    # Load files
    print(paste('Loading', raw_EVIfile, '...'))
    rawEVI <- read_csv(raw_EVIfile)
    print(paste('Loading', raw_NPPfile, '...'))
    rawNPP <- read_csv(raw_NPPfile)
    print(paste('Loading', raw_GPPfile, '...'))
    rawGPP <- read_csv(raw_GPPfile)

    print('Done.')
    
    rawNPP$date <- as.POSIXct(rawNPP[[quo_name(dateCol)]], '%Y-%m-%dT%H:%M:%S')
    rawNPP$year <- as.numeric(format(rawNPP$date, '%Y'))

    newNPP_ann <- rawNPP %>% rename(stationid=id, npp_qc=QC,
				    npp_sum=annualNPP) %>%
        select(-"system:index", -".geo")

    # Function to calculate monthly and annual dataframes from
    # raw MODIS dataframes
    daily_to_ann <- function(raw, varCol, dateCol, freq=8, varMult=1){
        # Get date and month column
        raw$date <- as.POSIXct(raw[[quo_name(dateCol)]], '%Y-%m-%dT%H:%M:%S')
        raw$year <- as.numeric(format(raw$date, '%Y'))
        raw$month <- as.numeric(format(raw$date, '%m'))
	raw[[quo_name(varCol)]] <- raw[[quo_name(varCol)]] * varMult

        # Create a monthly MODIS dataframe
        # Have to do funny symbol unquoting for the column variables
        raw_m <- raw %>% rename(stationid=id) %>%
          select(-"system:index", -".geo") %>%
          group_by(stationid, year, month) %>% 
          summarise(date=max(!!dateCol), mean=mean(!!varCol),
                    max=max(!!varCol),min=min(!!varCol))
                
        # Get # days in month and integrate mean gpp for months
        raw_m['ndays'] <- days_in_month(raw_m$date)
        raw_m['mint'] <- raw_m$mean/freq * raw_m$ndays

        # Calculate z score of monthly integrated MODIS var for each site
        sites <- unique(raw_m$stationid)

        for (s in sites){
            siteVar <- subset(raw_m, stationid==s)
            siteVar_m <- mean(siteVar$mint, na.rm=T)
            siteVar_sd <- sd(siteVar$mint, na.rm=T)
            siteVar$z <- (siteVar$mint - siteVar_m)/siteVar_sd
            if (s==sites[1]){
                newVar_m <- siteVar
            }
            else{
                newVar_m <- rbind(newVar_m, siteVar)
            }
        }

        # Create annual dataframe from the monthly one
        newVar_ann <- newVar_m %>% group_by(stationid, year) %>% 
	    summarise(amean=mean(mean), amax=max(max),
		      amin=min(min), sum=sum(mint),
		      z_mean=mean(z), z_max=max(z),
                      z_min=min(z),z_sum=sum(z))

        # Create a growing season dataframe from the monthly one
        newVar_gs <- newVar_m %>% group_by(stationid, year) %>% 
            filter(month>=4 & month<=9) %>%
            summarise(gsmean=mean(mean), gssum=sum(mean),
                      z_gsmean=mean(z), z_gssum=sum(z))

	# Merge annual and growing season dataframes, then return
	out <- merge(newVar_ann, newVar_gs, by=c('stationid', 'year'))
	return(out)
    }

    newGPP_ann <- daily_to_ann(rawGPP, gppCol, dateCol) %>%
	rename(gpp_mean=amean, gpp_max=amax, gpp_min=amin, gpp_sum=sum,
	       zgpp_mean=z_mean, zgpp_max=z_max, zgpp_min=z_min, zgpp_sum=z_sum,
	       gpp_gsmean=gsmean, gpp_gssum=gssum, zgpp_gsmean=z_gsmean,
	       zgpp_gssum=z_gssum)


    newEVI_ann <- daily_to_ann(rawEVI, eviCol, dateCol, freq=16,
			       varMult=0.0001) %>%
	rename(evi_mean=amean, evi_max=amax, evi_min=amin, evi_sum=sum,
	       zevi_mean=z_mean, zevi_max=z_max, zevi_min=z_min, zevi_sum=z_sum,
	       evi_gsmean=gsmean, evi_gssum=gssum, zevi_gsmean=z_gsmean,
	       zevi_gssum=z_gssum)

    # Merge annual EVI, NPP, GPP, dataframes
    newMOD_ann <- merge(newEVI_ann, newGPP_ann, by=c('stationid', 'year'))
    newMOD_ann <- merge(newMOD_ann, newNPP_ann, by=c('stationid', 'year'))
    
    # Write file
    out_fname <- paste0(outpath, 'annual_', hcn_network,
                        '_eeMODIS', '.csv')
    print(paste('Writing', out_fname, '...'))
    write_csv(newMOD_ann, out_fname)
    print('Done.')
}

