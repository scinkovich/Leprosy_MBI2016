# Leprosy_MBI2016
Leprosy project for the US-CAN Institutes Epidemiology Summer School

Species distribution modeling (SDM) relies on presence and absence data of a species to make accurate estimations from the environmental factors. Though there are a few work-arounds for generating pseudo-absence data, it would not be sufficient in the case of armadillos. Armadillos continue to pop up in areas not previously known to support populations and because of this our presence data for armadillos is not complete or reliable. Using a spatially-explicit habitat suitablity risk index model takes out the assumptions that must be made when running a SDM.  

##### Temperature for armadillos  
  + Regions with mean January temperatures of -2 deg C limits armadillo movement.  
  + The ability of yearlings to survive winter climatic extremes may limit the establishment of a permanent range at the northern boundaries [@Udvardy1969, @Taulman1996].  
  + The optimal temperature for in vivo growth of _M. leprae_ has been shown to be 33 deg C [@Shepard1965], therefore the hotter temperatures are receiving the highest risk as a two-fold effect of good conditions for the armadillos and bacterium.  
  + Data was retrieved online from WorldClim.  
  + **NOTE**: The temperatures in the raster "temp_us" are reported as deg C *10.     

&nbsp;  
  
##### Land type use by Armadillos 
  + Armadillos inhabit dense shady cover, such as brush, woodland or pine forests [@IFAS]. Armadillos are often considered a pest species because they often burrow under households or driveways [@IFAS].   
  + Land type measures the risk of transmission for _M. leprae_ by ranking the proportion of land cover used by armadillos. Using species data from GBIF, known locations of armadillos were plotted onto a land cover raster. After extracting the values (categorical) for land type at each plotted location, the highest proportions of armadillos were found in grasslands (0.29). This cover type was assigned the highest risk of 1.0 and the proportion of all other locations were divided by the max proportion to obtain values of risk. 
  + Data was 2010 data retrieved onine from MODIS-based Global Land Cover Climatology data from USGS LCI.    
  
&nbsp;   

##### Annual precipitation  
  + Armadillo range in the United States may be primarily restricted to regions receiving at least 380 mm yearly precipitation [@Taulman1996], though there have been armadillos reported in areas with about 340 mm of yearly rain.  
  + Since armadillos are only bound by the lower limits of rain per year, those above the 380mm threshold had a 1.0 probability of being present - based on annual precipiation alone.  
  + *NOTE*: Raster reports precipitation in mm.  

&nbsp;  

##### Soil moisture    
  + Soil texture is a factor in the animal's habitat selection and often burrow up to 100cm into the soil [@Taulman1996].
  + They prefer sandy or loam soils that are relatively easy to excavate [@IFAS].  
  + Wet soil is optimal habitat for the invertebrates (such as worms) that armadillos predominantly feed on.  
  + 40-100 cm soil moisture data were retrieved online from Global Land Data Assimilation System 2 (GLDAS-2).  
  + *NOTE*: The raster is reported moisture in $\frac{kg}{m^2}$.

&nbsp;   

##### Temperature Seasonality   
  + The idea behind using the temperature seasonality (sd*100) is that if the area is a high risk area for armadillos (in any of the other variables) and has a small annual temperature variation, then that area is more stable habitat.  
  + This more stable habitat would, in theory, be more likely to host larger armadillo populations and thus increase the risk of Leprosy to human populations.    

&nbsp;  

##### Final Risk Index (RI)  
  + The overall risk in an area averages the risk over all variables at that point.  

