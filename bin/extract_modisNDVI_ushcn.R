### UNM DroughtNet MODIS NDVI
# Formerly named AJ_DroughtNetMODIS

# This script downloads and processes MODIS NDVI data from North America and arranges it into a user-friendly format
# Current version written by Alesia

# Install required packages
#install.packages("ncdf4")
#install.packages("raster")

# Load required libraries
library(ncdf4)
library(raster)

# Write .csv header
header <- c("stationid", "latitude", "longitude", "year", "MODIS.end.DOY", "mean.MODIS.NDVI")
write.table(t(header), "C:/DroughtNet/MODIS_data/MODISoutput_ushcn.csv", row.names=F, col.names=F, sep=",")

# Get coordinates from met stations
ushcn_coords <- read.csv("C:/inetpub/ftproot/ushcn_allstations_calcs.csv", sep=',', header=T, strip.white=T)

# Read in MODIS data one year at a time
# Put in a loop so we read all years
allMODISfiles <- paste0('C:/inetpub/ftproot/DNdata/',
                        list.files('C:/inetpub/ftproot/DNdata/', 
                                   pattern='.unaccum.nc4'))

for (yr in 1:length(allMODISfiles)) {
  
  filename <- allMODISfiles[yr]
  MODIS.USA <- nc_open(filename)
  # Change depending on filename - select year from filename
  year <- unlist(substr(strsplit(filename, split="[.]")[[1]][2], 2, 5))
  
  #MODIS.USA
  # Get julian dates for each 8 day composite 
  MODIS.USA.time <- ncvar_get(MODIS.USA, 'time_bnds')
  
  # Read in latitude and longitude matrices, upper left pixel coordinate
  MODIS.latitudes <- ncvar_get(MODIS.USA, 'lat')
  MODIS.longitudes <- ncvar_get(MODIS.USA, 'lon')
  
  # Create coordinate matrix for selecting subsets,
  # Breaks up the search area into 27 blocks, making the code run faster
  lat.long.clip <- data.frame(cluster=NA, start.x=NA, start.y=NA, count.x=NA, count.y=NA)
  lat.long.clip[1,] <- c("a", 1, 10735, 2759-1, 13076-10735)
  lat.long.clip[2,] <- c("b", 2150, 10109, 5775-2150, 13076-10109)
  lat.long.clip[3,] <- c("c", 5468, 9899, 8853-5468, 12494-9899)
  lat.long.clip[4,] <- c("d", 8852, 9899, 12236-8852, 12494-9899)
  lat.long.clip[5,] <- c("e", 11930, 10109, 15554-11930, 13076-10109)
  lat.long.clip[6,] <- c("f", 15400, 10735, 20250-15400, 14277-10735)
  lat.long.clip[7,] <- c("g", 1, 8387, 2151-1, 11385-8387)
  lat.long.clip[8,] <- c("h", 1579, 7222, 5469-1579, 10736-7222)
  lat.long.clip[9,] <- c("i", 5182, 7222, 8853-5182, 10110-7222)
  lat.long.clip[10,] <- c("j", 8852, 7500, 12522-8852, 10110-7500)
  lat.long.clip[11,] <- c("k", 12235, 7722, 16225-12235, 10736-7722)
  lat.long.clip[12,] <- c("l", 15553, 8387, 20250-15553, 12437-8387)
  lat.long.clip[13,] <- c("m", 1, 6036, 1580-1, 8830-6036)
  lat.long.clip[14,] <- c("n", 1048, 5338, 5183-1048, 8388-5338)
  lat.long.clip[15,] <- c("o", 4916, 5106, 8853-4916, 7501-5106)
  lat.long.clip[16,] <- c("p", 8852, 5106, 12788-8852, 7723-5106)
  lat.long.clip[17,] <- c("q", 12521, 5338, 16656-12521, 8388-5338)
  lat.long.clip[18,] <- c("r", 16124, 6036, 20250-16124, 9492-6036)
  lat.long.clip[19,] <- c("s", 1, 2962, 4917-1, 6311-2962)
  lat.long.clip[20,] <- c("t", 4670, 2720, 8853-4670, 5339-2720)
  lat.long.clip[21,] <- c("u", 8852, 2720, 13034-8852, 5339-2720)
  lat.long.clip[22,] <- c("v", 12787, 2962, 17147-12787, 6037-2962)
  lat.long.clip[23,] <- c("w", 16655, 3686, 20250-16655, 7143-3686)
  lat.long.clip[24,] <- c("x", 4392, 1, 8853-4392, 2963-1)
  lat.long.clip[25,] <- c("y", 8852, 1, 13312-8852, 2963-1)
  lat.long.clip[26,] <- c("z", 13033, 1, 18000-13033, 3687-1)
  lat.long.clip[27,] <- c("unk", 1, 1, 20249, 14276)
  lat.long.clip$start.x <- as.numeric(lat.long.clip$start.x)
  lat.long.clip$start.y <- as.numeric(lat.long.clip$start.y)
  lat.long.clip$count.x <- as.numeric(lat.long.clip$count.x)
  lat.long.clip$count.y <- as.numeric(lat.long.clip$count.y)
  
  # Loop through each site
  for (site.num in 1:length(ushcn_coords$latitude)) {
    #length(ghcn_coords$latitude)
    print(paste("Year:", year, "Site being processed:", site.num, "/", length(ushcn_coords$latitude))) 
    
    # Extract lat and long from each met station
    site.id <- as.character(ushcn_coords$stationid[site.num])
    lat <- ushcn_coords$latitude[site.num]
    long <- ushcn_coords$longitude[site.num]
    
    # Find MODIS index values of desired lat/long
    # Select the search area from lat.long.clip
    if (lat>50) next
    if (long>-66) next
    if (lat>45 & long<(-120)) {cluster = 1}
    if (lat>45 & long<(-110) & long>(-120)) {cluster = 2}
    if (lat>45 & long<(-100) & long>(-110)) {cluster = 3}
    if (lat>45 & long<(-90) & long>(-100)) {cluster = 4}
    if (lat>45 & long<(-80) & long>(-90)) {cluster = 5}
    if (lat>45 & long>(-80)) {cluster = 6}
    if (lat>40 & lat<45 & long>(-120)) {cluster = 7}
    if (lat>40 & lat<45 & long<(-110) & long>(-120)) {cluster = 8}
    if (lat>40 & lat<45 & long<(-100) & long>(-110)) {cluster = 9}
    if (lat>40 & lat<45 & long<(-90) & long>(-100)) {cluster = 10}
    if (lat>40 & lat<45 & long<(-80) & long>(-90)) {cluster = 11}
    if (lat>40 & lat<45 & long>(-80)) {cluster = 12}
    if (lat>35 & lat<40 & long<(-120)) {cluster = 13}
    if (lat>35 & lat<40 & long<(-110) & long>(-120)) {cluster = 14}
    if (lat>35 & lat<40 & long<(-100) & long>(-110)) {cluster = 15}
    if (lat>35 & lat<40 & long<(-90) & long>(-100)) {cluster = 16}
    if (lat>35 & lat<40 & long<(-80) & long>(-90)) {cluster = 17}
    if (lat>35 & lat<40 & long>(-80)) {cluster = 18}
    if (lat>30 & lat<35 & long<(-110)) {cluster = 19}
    if (lat>30 & lat<35 & long<(-100) & long>(-110)) {cluster = 20}
    if (lat>30 & lat<35 & long<(-90) & long>(-100)) {cluster = 21}
    if (lat>30 & lat<35 & long<(-80) & long>(-90)) {cluster = 22}
    if (lat>30 & lat<35 & long>(-80)) {cluster = 23}
    if (lat<30 & long<(-100)) {cluster = 24}
    if (lat<30 & long<(-90) & long>(-100)) {cluster = 25}
    if (lat<30 & long>(-90)) {cluster = 26}
    
    # Get MODIS latitude and then longitude for the cluster
    MODIS.clip.lat <- ncvar_get(MODIS.USA, 'lat',
                                start=c(lat.long.clip$start.x[cluster], lat.long.clip$start.y[cluster]),
                                count=c(lat.long.clip$count.x[cluster], lat.long.clip$count.y[cluster]))
    
    MODIS.clip.long <- ncvar_get(MODIS.USA, 'lon',
                                 start=c(lat.long.clip$start.x[cluster], lat.long.clip$start.y[cluster]),
                                 count=c(lat.long.clip$count.x[cluster], lat.long.clip$count.y[cluster]))
    # Calculate distance from station lat/lon to MODIS corner
    lat.long.dif <- abs(MODIS.clip.lat - lat) + abs(abs(MODIS.clip.long) - abs(long))
    # Find minimum lat and long distance which should correspond to the corner
    # of the pixel containing the station
    lat.long.indices <- data.frame(which(
      lat.long.dif == min(lat.long.dif, na.rm=T), arr.ind=T))
    
    # Subset NDVI raster to 10km x 10km region around met station for each date
    # 41 by 41 pixels, 46 aquisitions per year
    MODIS.NDVI.clip <- ncvar_get(MODIS.USA, 'NDVI',
                                 start=c(lat.long.clip$start.x[cluster] + lat.long.indices$row[1]-20, 
                                         lat.long.clip$start.y[cluster] + lat.long.indices$col[1]-20, 1),
                                 count=c(41,41,46))
    # Make a timeseries of the mean of each acquisitions
    for (time.num in 1:46) {
      # Mask for landcover here? then calculate stats
      mean.NDVI <- mean(MODIS.NDVI.clip[,,time.num])
      #image(1:41, 1:41, MODIS.NDVI.clip[,,time.num], col = (brewer.pal(10, "BrBG")))
      
      # Write out new line of data
      newline <- c(site.id, lat, long, year, MODIS.USA.time[1,time.num], mean.NDVI)
      write.table(t(newline),"C:/DroughtNet/MODIS_data/MODISoutput_ushcn.csv", 
                  append=T, row.names=F, col.names=F, sep=",")
    }
  }
}




