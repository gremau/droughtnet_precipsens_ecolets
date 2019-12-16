library('readr')
library('dplyr')
library('lubridate')
library('reshape2')
library('SPEI')
library('xts')
source('xhcn_data_reduction.R')

get_spei <- function(cwdiff, freq=12, int_per=6, plot=TRUE){
    # Helper function to calculate spei
    # This comes from my climate tools repo on github

    # Get start date
    startmon <- as.yearmon(index(cwdiff[1]))
    startyr <- floor(as.numeric(startmon))
    startmon <- as.numeric( format( startmon, '%m'))

    # Get spei for that integration period
    spei_int <- spei(ts(cwdiff, frequency=freq, start=c(startyr, startmon)),
                     int_per, na.rm=T)

    # Check for invalid values
    values <- spei_int$fitted
    if (sum(is.na(values)) > int_per-1){
        browser()
    }
    values[!is.finite(values)] <- NA
    if (sum(is.na(values)) > (int_per - 1)){
        print('WARNING!!! - there are invalid values in the SPEI series')
    }
    if (plot){
        plot(spei_int)
    }
    return(spei_int)
}

get_spi <- function(prcp, freq=12, int_per=6, plot=TRUE){
    # Helper function to calculate spi
    # This comes from my climate tools repo on github

    # Get start date
    startmon <- as.yearmon(index(prcp[1]))
    startyr <- floor(as.numeric(startmon))
    startmon <- as.numeric( format( startmon, '%m'))

    # Get spi for that integration period
    spi_int <- spei(ts(prcp, frequency=freq, start=c(startyr, startmon)),
                     int_per, na.rm=T)

    # Check for invalid values
    values <- spi_int$fitted
    if (sum(is.na(values)) > int_per-1){
        browser()
    }
    values[!is.finite(values)] <- NA
    if (sum(is.na(values)) > (int_per - 1)){
        print('WARNING!!! - there are invalid values in the SPI series')
    }
    if (plot){
        plot(spi_int)
    }
    return(spi_int)
}

get_cv <- function(spei_xts, window_yr){
    # Rolling CV calculation
    # Specify window width
    wid <- window_yr*12
    # Calculate adjusted spei (>0 for calculating CV) and get its windowed mean
    spei_adj <- spei_xts - min(spei_xts, na.rm=TRUE)
    spei_adj_mean <- rollapply(spei_adj, wid, mean, align='right', fill=NA)

    # Standard deviation of adjusted data
    spei_adj_sd <- rollapply(spei_adj, wid, sd, align='right', fill=NA)
    # CV
    spei_adj_cv <- spei_adj_sd/spei_adj_mean
    return(spei_adj_cv)
}

monthly_xhcn_calcs <- function(monthly_file, inventoryfile, outpath,
                               hcn_network=NULL, name=NULL){
    # Calculate monthly SPEI and CV of SPEI values for each site
    # Calculate length of record and percent missing data per site
    # Output monthly_spei file and ushcn_allstations_calcs file

    # Check optional inputs
    if (is.null(hcn_network)){stop('Please enter HCN network')}
    if (is.null(name)){stop('Please enter a file name string')}
    # Create ignore file
    if (hcn_network=='ghcn'){
        ignore <- c('USC00102667','USC00187230','USC00480432')
    } else {
        ignore <- c('None')
    }

    # Load data
    print(paste('Loading', monthly_file, '...'))
    df <- read_csv(monthly_file,
                   col_types=cols(stationid = col_character(),
                                  date = col_date(format = ""),
                                  prcp = col_double(),
                                  tavg = col_double(),
                                  tmax = col_double(),
                                  tmin = col_double()))
    print('Done.')
    #modify name to keep file naming convention
    if (name!=''){name <- paste0('_', name)}
    # Load the inventory of sites
    if (hcn_network=='ghcn'){
        print('Fetch GHCN inventory...')
        df.inv <- read_csv(inventoryfile)
    } else if (hcn_network=='huxman') {
	print('Fetch Huxman inventory...')
	df.inv <- read_csv(inventoryfile,
		    col_types = cols(inconus = col_logical()))
	df.inv <- df.inv[df.inv$inconus,]
    } else if (hcn_network=='ushcn') {
        print('Fetch USHCN inventory...')
        df.inv <- load_ushcn_inventory(inventoryfile)
    } else { stop('No HCN network given!')
    }
    # Add number of years and missing data columns
    df.inv$nyears <- NA
    df.inv$pctmissing <- NA
    print('Done.')
    
    # Get list of sites
    sites <- unique(df$stationid)
    print(paste('There are ', length(sites), 'sites in monthly datafile.' ))
    sites_good <- sites[!(sites%in%ignore)]
    print(paste(length(sites)-length(sites_good), 'sites are being ignored.' ))

    count <- 1
    for (s in sites_good){
        # Subset the data by site s
        df.site <- subset(df, stationid==s)
        # If tavg is NA, try to fill with average of 
        # tmax and tmin
        masktavg <- is.na(df.site$tavg)
        mmtavg <- (df.site$tmax + df.site$tmin)/2
        df.site$tavg[masktavg] <- mmtavg[masktavg]

        # Get latitude from the inventory
        lat <- df.inv$latitude[df.inv$stationid==s]
        print(paste('Calculating SPEI for site', s, 'at latitude', lat))
        # Create XTS object from precip and temp columns
        df.site.xts <- xts(df.site[,c('prcp', 'tavg')],
                           order.by=as.Date(df.site$date, "%Y-%m-%d"))
        # Remove leading and trailing rows with NAs
        df.site.xts <- na.trim(df.site.xts, is.na='any')
        masktavg <- is.na(df.site.xts$tavg)
        maskprcp <- is.na(df.site.xts$prcp)
        # Print some stats, add to inventory, fill missing values
        nyears <- length(maskprcp)/12
        pctmissing <- 100 * (sum(masktavg | maskprcp)/length(maskprcp))
        print(paste0(as.character(nyears), ' years of met data'))
        print(paste0(as.character(sum(masktavg | maskprcp)), ' months missing'))
        print(paste0(as.character(pctmissing), ' percent missing'))
        df.inv[df.inv$stationid==s, 'nyears'] <- nyears
        df.inv[df.inv$stationid==s, 'pctmissing'] <- pctmissing
        df.site.xts$tavg <- na.approx(df.site.xts$tavg)
        df.site.xts$prcp <- na.approx(df.site.xts$prcp)
        # Note T already converted (from deg/100) in monthly data
        df.site.xts$pet <- NA
        #mask <- is.na(df.site.xts$tavg)
        petin <- thornthwaite(df.site.xts$tavg, lat)
        df.site.xts[, 'pet'] <- petin

        # Calculate climatic water differential
        # Note P already converted (from mm/10) in monthly data
        cwdiff <- df.site.xts$prcp - df.site.xts$pet
        # And SPEI at 3, 6, 9, 12, 18, and 24mo
        spei3mo <- get_spei(cwdiff, int_per=3, plot=FALSE)
        df.site.xts$spei3mo <- xts(as.vector(spei3mo$fitted),
                                   order.by=index(cwdiff))
        spei6mo <- get_spei(cwdiff, plot=FALSE)
        df.site.xts$spei6mo <- xts(as.vector(spei6mo$fitted),
                                   order.by=index(cwdiff))
        spei9mo <- get_spei(cwdiff, int_per=9, plot=FALSE)
        df.site.xts$spei9mo <- xts(as.vector(spei9mo$fitted),
                                   order.by=index(cwdiff))
        spei12mo <- get_spei(cwdiff, int_per=12, plot=FALSE)
        df.site.xts$spei12mo <- xts(as.vector(spei12mo$fitted), 
                                    order.by=index(cwdiff))
        spei18mo <- get_spei(cwdiff, int_per=18, plot=FALSE)
        df.site.xts$spei18mo <- xts(as.vector(spei18mo$fitted), 
                                    order.by=index(cwdiff))
        spei24mo <- get_spei(cwdiff, int_per=24, plot=FALSE)
        df.site.xts$spei24mo <- xts(as.vector(spei24mo$fitted), 
                                    order.by=index(cwdiff))
	# And SPI at 12 mo
	spi12mo <- get_spi(df.site.xts$prcp, int_per=12, plot=FALSE)
        df.site.xts$spi12mo <- xts(as.vector(spi12mo$fitted), 
                                    order.by=index(df.site.xts$prcp))
        # Get 5 and 10 year rolling coefficient of variation
	df.site.xts$tavg_cv5yr <- get_cv(df.site.xts$tavg, 5)
        df.site.xts$tavg_cv10yr <- get_cv(df.site.xts$tavg, 10)
	df.site.xts$prcp_cv5yr <- get_cv(df.site.xts$prcp, 5)
        df.site.xts$prcp_cv10yr <- get_cv(df.site.xts$prcp, 10)
        df.site.xts$spei12mo_cv5yr <- get_cv(df.site.xts$spei12mo, 5)
        df.site.xts$spei12mo_cv10yr <- get_cv(df.site.xts$spei12mo, 10)
	df.site.xts$spi12mo_cv5yr <- get_cv(df.site.xts$spi12mo, 5)
        df.site.xts$spi12mo_cv10yr <- get_cv(df.site.xts$spi12mo, 10)
        # Build new dataframe
        newdf <- data.frame(date=index(df.site.xts), coredata(df.site.xts))
        newdf$stationid <- s
        newdf <- newdf[,c(ncol(newdf), 1, 3:ncol(newdf)-1)]
        if (count < 2){
            df.site.xts.rbind <- newdf
        } else {
            df.site.xts.rbind <- rbind(df.site.xts.rbind, newdf)
        }
        count <- count + 1   
    }
    # Write out the monthly file
    out_fname1 <- paste0(outpath, 'monthly_', hcn_network, '_spei', name, '.csv')
    print(paste('Writing output file', out_fname1))
    #write.zoo(df.sites.xts.rbind, out_fname, sep=',')
    write_csv(df.site.xts.rbind, out_fname1)

    out_fname2 <- paste0(outpath, 'allyr_', hcn_network, '_calcs', name, '.csv')
    print(paste('Writing output file', out_fname2))
    #write.zoo(df.sites.xts.rbind, out_fname, sep=',')
    write_csv(df.inv, out_fname2)
    return(list(out_fname1, out_fname2))

}

spei_trend_calcs <- function(spei_mon_file, stncalcsfile, outpath,
                             hcn_network=NULL, name=NULL){
    # Check optional inputs
    if (is.null(hcn_network)){stop('Please enter HCN network')}
    if (is.null(name)){stop('Please enter a file name string')}
    library('forecast')
    # Read in data
    speidat <- read_csv(spei_mon_file)

    # Load the station_calcs file ( comes from monthly_xhcn_calcs() )
    if (hcn_network=='ghcn'){
        print('Fetch GHCN station data...')
        inv <- read_csv(stncalcsfile)
    } else if (hcn_network=='huxman') {
	print('Fetch Huxman station data...')
	inv <- read_csv(stncalcsfile)
    } else if (hcn_network=='ushcn') {
        print('Fetch USHCN station_Data...')
        inv <- read_csv(stncalcsfile)
    } else { stop('No HCN network given!')
    }
    print('Done.')
    inv$tavg_trend <- NA
    inv$tavg_trend_sig <- NA
    inv$tavg_cv5yr_trend <- NA
    inv$tavg_cv5yr_trend_sig <- NA
    inv$prcp_trend <- NA
    inv$prcp_trend_sig <- NA
    inv$prcp_cv5yr_trend <- NA
    inv$prcp_cv5yr_trend_sig <- NA
    inv$spei12mo_trend <- NA
    inv$spei12mo_trend_sig <- NA
    inv$spei12mo_cv5yr_trend <- NA
    inv$spei12mo_cv5yr_trend_sig <- NA
    inv$spi12mo_trend <- NA
    inv$spi12mo_trend_sig <- NA
    inv$spi12mo_cv5yr_trend <- NA
    inv$spi12mo_cv5yr_trend_sig <- NA
    sites <- unique(speidat$stationid)
    for (i in sites){
        # Subset site
        site.df <- subset(speidat, stationid==i)
    	# Get inventory row
    	site.inv <- subset(inv, stationid==i)
        # If site is less than 85 years flag site
        shortflag <- site.inv$nyears < 100
        # If more than 15% of data are missing flag site (trend will be NA)
        missingflag <- site.inv$pctmissing >= 16
        if (!(shortflag | missingflag)) {
	    # XTS object
            site.df.xts <- xts(site.df[,3:ncol(site.df)],
                            order.by=as.Date(site.df$date, "%Y-%m-%d"))
	    # TAVG trend
	    minyr <- format(index(site.df.xts)[1], format = "%Y")
            minmo <- format(index(site.df.xts)[1], format = "%m")
            ## Sometimes the spei calculator adds some Inf values - change to NA
            #site.df.xts$tavg[is.infinite(site.df.xts$spei12mo)] <- NA
            tavg_ts <- ts(site.df.xts$tavg[1:length(site.df.xts$tavg)],
                          frequency=12, start=c(as.numeric(minyr),
                                                as.numeric(minmo)))
            tavg.fit <- tslm(tavg_ts ~ trend)
            # Get fit statistics and put in inventory
            invsite <- inv$stationid==i
            inv$tavg_trend[invsite] <- coefficients(tavg.fit)[2]
            inv$tavg_trend_sig[invsite] <- 
                summary(tavg.fit)$coefficients['trend',4]
	    
	    # PRCP trend
            ## Sometimes the spei calculator adds some Inf values - change to NA
            #site.df.xts$tavg[is.infinite(site.df.xts$spei12mo)] <- NA
            prcp_ts <- ts(site.df.xts$prcp[1:length(site.df.xts$prcp)],
                          frequency=12, start=c(as.numeric(minyr),
                                                as.numeric(minmo)))
            prcp.fit <- tslm(prcp_ts ~ trend)
            # Get fit statistics and put in inventory
            inv$prcp_trend[invsite] <- coefficients(prcp.fit)[2]
            inv$prcp_trend_sig[invsite] <- 
                summary(prcp.fit)$coefficients['trend',4]
	
	    # 12mo SPEI trend
            minyr <- format(index(site.df.xts)[12], format = "%Y")
            minmo <- format(index(site.df.xts)[12], format = "%m")
            # Sometimes the spei calculator adds some Inf values - change to NA
            site.df.xts$spei12mo[is.infinite(site.df.xts$spei12mo)] <- NA
            spei_ts <- ts(site.df.xts$spei12mo[12:length(site.df.xts$spei12mo)],
                          frequency=12, start=c(as.numeric(minyr),
                                                as.numeric(minmo)))
            spei.fit <- tslm(spei_ts ~ trend)
            # Get fit statistics and put in inventory
            inv$spei12mo_trend[invsite] <- coefficients(spei.fit)[2]
            inv$spei12mo_trend_sig[invsite] <- 
                summary(spei.fit)$coefficients['trend',4]

	    # 12mo SPI trend
            # Sometimes the spei calculator adds some Inf values - change to NA
            site.df.xts$spi12mo[is.infinite(site.df.xts$spi12mo)] <- NA
            spi_ts <- ts(site.df.xts$spi12mo[12:length(site.df.xts$spi12mo)],
                          frequency=12, start=c(as.numeric(minyr),
                                                as.numeric(minmo)))
            spi.fit <- tslm(spi_ts ~ trend)
            # Get fit statistics and put in inventory
            inv$spi12mo_trend[invsite] <- coefficients(spi.fit)[2]
            inv$spi12mo_trend_sig[invsite] <- 
                summary(spi.fit)$coefficients['trend',4]

            # CV5 trends
            minyr <- format(index(site.df.xts)[6*12], format = "%Y")
            minmo <- format(index(site.df.xts)[6*12], format = "%m")
            # TAVG and PRCP
	    tavgcv_ts <- ts(site.df.xts$tavg_cv5yr[
                        (6*12):length(site.df.xts$tavg_cv5yr)],
                frequency=12, start=c(as.numeric(minyr),
                                      as.numeric(minmo)))
	    tavgcv.fit <- tslm(tavgcv_ts ~ trend)
            # Get fit statistics and put in inventory
            inv$tavg_cv5yr_trend[invsite] <- coefficients(tavgcv.fit)[2]
            inv$tavg_cv5yr_trend_sig[invsite] <- 
                summary(tavgcv.fit)$coefficients['trend',4]

	    prcpcv_ts <- ts(site.df.xts$prcp_cv5yr[
                        (6*12):length(site.df.xts$prcp_cv5yr)],
                frequency=12, start=c(as.numeric(minyr),
                                      as.numeric(minmo)))
	    prcpcv.fit <- tslm(prcpcv_ts ~ trend)
            # Get fit statistics and put in inventory
            inv$prcp_cv5yr_trend[invsite] <- coefficients(prcpcv.fit)[2]
            inv$prcp_cv5yr_trend_sig[invsite] <- 
                summary(prcpcv.fit)$coefficients['trend',4]

	    # SPEI CV
	    # Sometimes the spei calculator adds some Inf values - change to NA
            site.df.xts$spei12mo_cv5yr[is.infinite(site.df.xts$spei12mo_cv5yr)
                                      ] <- NA
            speicv_ts <- ts(site.df.xts$spei12mo_cv5yr[
                        (6*12):length(site.df.xts$spei12mo_cv5yr)],
                frequency=12, start=c(as.numeric(minyr),
                                      as.numeric(minmo)))
            tryCatch(
                {speicv.fit <- tslm(speicv_ts ~ trend)
                # Get fit statistics and put in inventory
                inv$spei12mo_cv5yr_trend[invsite] <- coefficients(speicv.fit)[2]
                inv$spei12mo_cv5yr_trend_sig[invsite] <- 
                    summary(speicv.fit)$coefficients['trend',4]
                },
                error=function(e) {
                    print('Missing values in CV')
                    # Choose a return value in case of error
                    inv$spei12mo_cv5yr_trend[invsite] <- NA
                    inv$spei12mo_cv5yr_trend_sig[invsite] <- NA
                }
	    )
	    # SPI CV
	    # Sometimes the spei calculator adds some Inf values - change to NA
            site.df.xts$spi12mo_cv5yr[is.infinite(site.df.xts$spi12mo_cv5yr)
                                      ] <- NA
            spicv_ts <- ts(site.df.xts$spi12mo_cv5yr[
                        (6*12):length(site.df.xts$spi12mo_cv5yr)],
                frequency=12, start=c(as.numeric(minyr),
                                      as.numeric(minmo)))
            tryCatch(
                {spicv.fit <- tslm(spicv_ts ~ trend)
                # Get fit statistics and put in inventory
                inv$spi12mo_cv5yr_trend[invsite] <- coefficients(spicv.fit)[2]
                inv$spi12mo_cv5yr_trend_sig[invsite] <- 
                    summary(spicv.fit)$coefficients['trend',4]
                },
                error=function(e) {
                    print('Missing values in CV')
                    # Choose a return value in case of error
                    inv$spi12mo_cv5yr_trend[invsite] <- NA
                    inv$spi12mo_cv5yr_trend_sig[invsite] <- NA
                }
            )
        
        } else {
            print(paste('Site', i, 'has <85 years or >50% NA, skipping'))
            invsite <- inv$stationid==i
            inv$spei12mo_trend[invsite] <- NA
            inv$spei12mo_trend_sig[invsite] <- NA
        }
    }
    out_fname <- stncalcsfile
    print(paste('Writing output file', out_fname))
    #write.zoo(df.sites.xts.rbind, out_fname, sep=',')
    write_csv(inv, paste(out_fname, sep=''))
}
    
