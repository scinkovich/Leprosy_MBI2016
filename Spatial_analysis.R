################################################################################## Model of Ebola virus disease (EVD) outbreak in Liberia to compare projections ## made with and without accounting for asymptomatic, non-infectious infections.
##
## US_CAN Institutes Epidemiology Summer School
## Leprosy Group
## Mentors: Jane Heffernan & Julien Arino
##
## Stephanie Cinkovich, June 27, 2016

## Make sure you install these all first using install.packages("")
library(dismo)
library(rgbif)
library(maps)
library(ggplot2)
library(spatial.tools)
library(rgdal)
library(raster)
library(rWBclimate)
library(scales) 
library(gridExtra)
library(rgeos)
library(maptools)
library(plyr)
library(Rmisc)
library(ggmap)
library(RColorBrewer)
library(microbenchmark)
library(compiler)

## Grabbing species occurrence data from Global Biodiversity Information Facility (GBIF)
## using the search terms of the Nine-banded armadillo scientific name
key <- name_suggest(q='Dasypus novemcinctus', rank='species')$key[1]
Dv <- occ_search(taxonKey=key, return='data')
gbifmap(Dv) #map the species presence on worldmap
locs <- subset(Dv, select = c("country", "decimalLatitude", "decimalLongitude")) #take out only location data
head(locs)  # prove that the correct headings are there
locs <- subset(locs, locs$decimalLatitude < 90) #gets rid of coordinates with errors
coordinates(locs) <- c("decimalLongitude", "decimalLatitude")  #set spatial coordinates
plot(locs) #plots on a coordinate grid
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
proj4string(locs) <- crs.geo  # define projection system of our data
summary(locs)

# load in US states outline
us <- getData("GADM", country="USA", level=1)

## Make pseudo-absences - only if you want to run a species distribution model (SDM)
#mask <- raster("tmean_fixed.grd")
#set.seed(1963)
#bg <- randomPoints(mask, 500 )
#par( mfrow = c( 1, 2 ), oma = c( 0, 0, 2, 0 ) )
#plot(!is.na(mask), legend=FALSE)
#points(bg, cex=0.5)
# now we repeat the sampling, but limit the area of sampling using a spatial extent
#e <- extent(-130, -65, 40, 55) 
#bg2 <- randomPoints(mask, 100, ext=e)
#plot(!is.na(mask), legend=FALSE)
#plot(e, add=TRUE, col='red')
#points(bg2, cex=0.5)
#title("Generating pseudo-absences", outer=TRUE)


## Climate: I used WorldClim or BioClim to download the annual precip and mean temperature
## >>> The smaller the arc-minutes, the higher the resolution of the data (30s is best)
## Soil: I used 40-100cm Soil Moisture data from Global Land Data Assimilation System 2 (GLDAS-2)
## Land Cover: I used 0.5km MODIS-based Global Land Cover Climatology data from USGS LCI
## Make sure all of your rasters are in your working directory folder

## Only necessary to run the first time in order to write out your working raster
usext <- c(-130, -65, 20, 55) #sets the extent of our maps to the US only
temp <- raster(paste(getwd(), "/MeanTempCurrent/tmean1.bil", sep = "")) #mean temp in January
temp.c <- crop(temp, usext) #crops the raster based on the US extent above

precip <- raster(paste(getwd(), "/BioClimCurrent/bio12.bil", sep = "")) #annual precip
precip.c <- crop(precip, usext)

tempvar <- raster(paste(getwd(), "/BioClimCurrent/bio4.bil", sep="")) #temperature seasonality
tempvar.c <- crop(tempvar, usext)

soil <- raster("soil.grd") #soil moisture
soil <- crop(soil, usext)
soil.c <- resample(soil, tmean1.c, method='bilinear') #change resolution to match other rasters

cover <- raster("LCType.tif") #Land cover type
rcover <- resample(cover, precip.c, method="ngb") #change dimensions to match other rasters

## Write all of those cropped rasters as new files for easy read-in for future work 
writeRaster(cover, filename="landcover_us.grd")
writeRaster(tmean1.c, filename="temp_us.grd")
writeRaster(precip.c, filename="precip_us.grd")
writeRaster(tempvar.c, filename="vartemp_us.grd")
writeRaster(soil.c, filename="soil_us.grd")

## Stack layers and extract their values @ armadillo points 
predictors.c <- stack(temp, tempvar, precip, land, soil)
presvals.c <- extract(predictors.c, locs)
presvals.c <- na.omit(presvals.c)#because there are points outside of the US
backgr.c <- randomPoints(predictors.c, 485)
absvals.c <- extract(predictors.c, backgr.c)
pb.c <- c(rep(1, nrow(presvals.c)), rep(0, nrow(absvals.c)))
sdmdata.c <- data.frame(cbind(pb.c, rbind(presvals.c, absvals.c)))
sdmdata.c[,'LCType'] = as.factor(sdmdata.c[,'LCType'])
head(sdmdata.c)
pairs(sdmdata.c, cex=0.1, fig=TRUE) #check for layer correlation

## Getting the proportions of Land Cover @ the armadillo pres points
dat = as.data.frame(count(sdmdata.c, "LCType")) #make a dataframe of the frequency
dat$prop=dat$freq/nrow(sdmdata.c) #add a column of proportions
dat$risk=round(dat$prop/max(dat$prop), 3) #add column for the risk of each factor
#use these values when assigning new values under "Raster Calculations"

#**************************************
# VARIOUS PLOTS 
#*************************************
## Model prediction of armadillo presence using SDM approach
m1.c <- glm(pb.c ~ tmean1 + LCType + bio12 + bio4 +Soil.moisture, data=sdmdata.c)
summary(m1.c)
p.c <- predict(predictors.c, m1.c)
dev.off()
plot(p.c, xlim=c(-130,-65), ylim=c(22,50))
plot(us, add=T)
dnov <- readShapePoly("Dnov_range/species_6290.shp")
plot(dnov, add=TRUE, density=c(5,10,15),  border=TRUE,xlim=c(-130, -65),ylim=c(20, 55))
title('Armadillo suitable habitats - Current')

## START HERE AFTER FIRST TIME RUNNING THE CODE ABOVE
## Reading in the saved rasters 
precip <- raster("precip_us.grd")
temp <- raster("temp_us.grd")
land <- raster("landcover_us.grd")
tempvar <- raster("vartemp_us.grd")
soil <- raster("soil_us.grd")

## Making a USA Leprosy Choropleth (# cases per state by color)
## You need to download a shapefile of the US - available online need all files in same folder
## but you only read in the ".shp" file

set.seed(8000) #makes sure that every time you run this you will get same outcome
states <- readShapeSpatial("states_21basic/states.shp") #read in the shapefile for the US
cases <- read.csv("Leprosy2004_2014.csv") #I entered case data from pdf online into a csv by state
head(states) #allows you to see what the column names are and what they contain
states.shp <- fortify(states, region = "STATE_ABBR") #makes the shapefile stable using the unique region
states.shp <- subset(states.shp, id != "AK") #exclude Alaska
states.shp.final <- subset(states.shp, id != "HI") #exclude Hawaii
states.data <- merge(states.shp.final, cases, by=("id"), all.x=TRUE) #merge case data with shapefile
#make sure the csv file has same "id" column so it is easy to merge
final.plot<-states.data[order(states.data$order), ] #puts the state data in order - better graphing
ggplot() + #map of leprosy cases from 2004-2014
  geom_polygon(data = final.plot, 
               aes(x=long, y=lat, group = group, fill = cases), 
               color = "black", size = 0.2) + 
  coord_map()+
  scale_fill_distiller(name="Cases", palette = "YlGn", breaks = pretty_breaks(n =5))+
  theme_nothing(legend = TRUE)+
  labs(title="Cumulative Cases of Leprosy from 2004-2014")

## READ IN SOIL DATA AND CONVERT TO RASTER
## This is how you read in ".nc4" file type... it is terrible. Avoid if possible.
library(RNetCDF)
library(ncdf4)
library(ncdf)
fname <- file.choose()
nc<-nc_open(fname) 
TAS <- raster(fname, varname="SoilMoi40_100cm_inst")
writeRaster(TAS, filename="soil.grd", bandorder='BIL', overwrite=TRUE)

#********************************************************
# CALCULATIONS ON RASTER DATA FOR RISK INDICES
#*******************************************************
usext <- c(-130, -65, 20, 55) #extent for US map again
precip <- raster("precip_us.grd")
temp <- raster("temp_us.grd")
land <- raster("landcover_us.grd")
tempvar <- raster("vartemp_us.grd")
soil <- raster("soil_us.grd")

## Annual Precipitation 
rain_rules <- function(r) {
  rr <- r #call in the raster as my matrix
  rr[rr < 340] <- 0 #any value < 340 will be assigned a 0 in new raster
  rr[rr < 380 & rr > 339] <- 0.5 #any value >339 and <380 is assigned 0.5 in new raster
  rr[rr >= 380] <- 1 #any value > or = to 380 will be assigned a 1
  r <- rr #replace the original raster with new values
  r #print or return back the new raster
}
rain_proc <- calc(precip, rain_rules) #calculate the function above on my "precip" data
RV3 <- rain_proc #save the new raster as my random variable 3 (RV3)

## Soil Moisture
soilMoist <- function(r){
  rr <- r
  out <- r/max(r, na.rm=TRUE)
  return(out)
}
soil_proc = calc(soil, soilMoist)
#above = take each value and divide by the maximum value of the raster, ignoring NA values
RV4 <- soil_proc #save the new raster as my random variable 4 (RV4)

## Mean Temperature in January
temp_rules <- function(r) {
  rr <- r
  rr/10 #convert to degree C 
  (rr*1.8)+32 #convert to degree F
  rr[rr < 14] <- 0 
  rr[rr < 29 & rr > 13] <- 0.25
  rr[rr < 33 & rr > 28] <- 0.5
  rr[rr > 32] <- 1
  r <- rr
  return(r)
}
temp_proc <- calc(temp, temp_rules) #might need to use temp not temp.1
RV1 <- temp_proc

## Seasonal Variation in Temperature
seasonal <- function(r){
  rr <- r/max(r, na.rm=TRUE)
  out <- 1-rr
  return(out)
}
tempvar_proc <- calc(tempvar, seasonal)
#above = 1 minus each value divided by max temp in raster, ignoring the NA values
#so you will be assigning rasters with lowest value the highest risk
RV5 <- tempvar_proc

## Land Cover 
## See earlier proportion analysis for values*
## The land cover data was categorical - that's why the raster "math" is long and ugly
pink_fairy_armadillos_rule = function(nine_banded_armadillos_are_loosers) {
  result = nine_banded_armadillos_are_loosers
  result[which(nine_banded_armadillos_are_loosers == 0)] = 0.158  #water
  result[which(nine_banded_armadillos_are_loosers == 1)] = 0.168 #evergreen needle leaf forest
  result[which(nine_banded_armadillos_are_loosers == 2)] = 0.015 #evergreen broadleaf forest
  result[which(nine_banded_armadillos_are_loosers == 3)] = 0 #deciduous needle leaf forest
  result[which(nine_banded_armadillos_are_loosers == 4)] = 0.125 #deciduous broadleaf forest
  result[which(nine_banded_armadillos_are_loosers == 5)] = 0.231 #mixed forest
  result[which(nine_banded_armadillos_are_loosers == 6)] = 0 #closed shrublands
  result[which(nine_banded_armadillos_are_loosers == 7)] = 0.271 #open shrublands
  result[which(nine_banded_armadillos_are_loosers == 8)] = 0.454 #woody savanna
  result[which(nine_banded_armadillos_are_loosers == 9)] = 0.015 #savannas
  result[which(nine_banded_armadillos_are_loosers == 10)] = 1 #grasslands
  result[which(nine_banded_armadillos_are_loosers == 11)] = 0.081 #permanent wetland
  result[which(nine_banded_armadillos_are_loosers == 12)] = 0.264 #croplands
  result[which(nine_banded_armadillos_are_loosers == 13)] = 0.341 #urban/built-up
  result[which(nine_banded_armadillos_are_loosers == 14)] = 0.238 #cropland/natural veg
  result[which(nine_banded_armadillos_are_loosers == 15)] = 0 #snow/ice
  result[which(nine_banded_armadillos_are_loosers == 16)] = 0.018 #barren/sparsely vegetated
  return(result)
}
land_processed = calc(land,pink_fairy_armadillos_rule)
RV2 <- land_processed

#***********************RISK INDICES PRODUCT*********************
RV.stack <- stack(RV1,RV2,RV3,RV4,RV5) #stack all RV's so they are in one dataset
flext <- c(-88, -78, 25, 34) #the extent for only Florida
RI <- calc(RV.stack, fun=function(x){(x[[1]]+x[[2]]+x[[3]]+x[[4]]+x[[5]])/5}) #take average of all RV's as final risk
RI.fl <- crop(RI, flext) #crop the final risk index for only Florida

#****************PLOTTING RISK MAPS****************************
fx.ggplot<-function(RV, title) { #define a function for all maps to be generated
  RV.p <- rasterToPoints(RV)
  df.RV <- data.frame(RV.p)
  colnames(df.RV) <- c('Longitude', 'Latitude', 'RV')
  ggplot(df.RV, aes(x=Longitude, y=Latitude)) +
      geom_raster(aes(fill=RV)) +
      ggtitle(title)+
      theme_bw() +
      ditch_the_axes + #defined below
      coord_equal() +
      scale_fill_distiller('Risk', limits=c(0,1), palette = 'Spectral')
}
ditch_the_axes <- theme( #function for you to get rid of lat/long axes
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank()
)

## Applying the rasterToPoints and mapping function on all variables
plot1 <- fx.ggplot(RV1, title="Mean Temp - Jan")
plot2 <- fx.ggplot(RV2, title="Land Cover")
plot3 <- fx.ggplot(RV3, title="Annual Precip")
plot4 <- fx.ggplot(RV4, title="Soil Moisture - Jul")
plot5 <- fx.ggplot(RV5, title="Temp Variation")
plot.final <- fx.ggplot(RI, title="Overall Risk") #the full risk index (RI) plot
plot.FL <- fx.ggplot(RI.fl, title="Overall Risk - FL") #the full RI - focused on FL
multiplot(plot1, plot2, plot3,plot4, plot5, cols=3) #all indiv. RI plots together