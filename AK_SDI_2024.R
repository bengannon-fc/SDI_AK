####################################################################################################
#### Suppression Difficulty Index, Alaska 2024
#### Author: Ben Gannon (benjamin.gannon@usda.gov)
#### Date Created: 04/01/2024
#### Last Modified: 05/18/2024
####################################################################################################
# This script will model Suppression Difficulty Index (SDI) as defined in Rodriguez y Silva et al. 
# 2020 with modifications for use in Alaska. This script was adapted from the earlier work of Matt 
# Panunto, Quresh Latif, and Pyrologix.
#
# This script assumes that the streets and trails shapefiles have already been clipped to the study
# area.
#
# The neighborhood weights matrix for the focal mean analyses took some experimentation to match
# the ArcGIS results. Know that the Arc NbrCircle tool does not behave as you might assume - it 
# usually makes a larger neighborhood than if starting from the center pixel.
# The biggest differences between the R and arcpy implementations are associated with the
# accessibility sub-index (distance from roads) and the penetrability sub-index (trail density). In
# general, differences are small enough (+/- 10) that there should be no change in interpretation.
####################################################################################################
#-> Set working directory
setwd('F:/AK_SDI/SDI')
####################################################################################################

###########################################START MESSAGE############################################
cat('Suppression Difficulty Index, Alaska 2024\n',sep='')
cat('Started at: ',as.character(Sys.time()),'\n\n',sep='')
cat('Messages, Errors, and Warnings (if they exist):\n')
####################################################################################################

############################################START SET UP############################################

#-> Load packages
.libPaths(.libPaths()[1]) # Specify libPath if there are multiple
packages <- c('terra')
for(package in packages){
	if(suppressMessages(!require(package,character.only=T))){
		install.packages(package,repos='https://repo.miserver.it.umich.edu/cran/')
		suppressMessages(library(package,character.only=T))
	}
}

#-> Maximize raster processing speed
terraOptions(memmax=10^9)

#-> Grab data and parameters from settings file
settings <- read.csv('SDI_inputs.csv',header=T)
inDEM <- settings[settings$Variable == 'DEM','Value'] # DEM raster (meters)
inFBFM40 <- settings[settings$Variable == 'FBFM40','Value'] # fire behavior fuel model raster
inFL <- settings[settings$Variable == 'FL','Value'] # FlamMap flame length raster (meters)
inHUA <- settings[settings$Variable == 'HUA','Value'] # FlamMap heat per unit area raster (kJ/m2)
inSA <- settings[settings$Variable == 'SA','Value'] # Study area shapefile for clipping
streets <- settings[settings$Variable == 'streets','Value'] # Streets/roads shapefile
trails <- settings[settings$Variable == 'trails','Value'] # Trails shapefile
opath <- settings[settings$Variable == 'opath','Value'] # Raster output path
oname <- settings[settings$Variable == 'oname','Value'] # Raster output name prefix

#############################################END SET UP#############################################	

###########################################START ANALYSIS###########################################

#-> Load study area extent
SA <- vect(inSA)

#-> Clip rasters to study area and create raster objects
FL <- crop(rast(inFL),SA)
HUA <- crop(rast(inHUA),SA)
DEM <- crop(rast(inDEM),FL) # Use FL raster to mitigate any extent differences
FBFM40 <- crop(rast(inFBFM40),FL)

#-> Standardize CRS info to avoid warnings (does not change projection)
crs(FBFM40) <- crs(DEM)
crs(FL) <- crs(DEM)
crs(HUA) <- crs(DEM)

#-> ENERGY BEHAVIOR SUB-INDEX
cat('#-> Calculating Energy Behavior Sub Index\n')

# Reclassify Flammap flamelength raster values based on Table 2 of Rodriguez y Silva et al
# Convert flame length units to meters before reclassifying, if necessary
FL_rtab <- data.frame(Low=c(-1,0.5,1,1.5,2,2.5,3,3.5,4,4.5),
                      High=c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,1000),
					  Val=c(1,2,3,4,5,6,7,8,9,10))
FL_recl <- classify(FL,FL_rtab)
cat('Reclassified flame lengths to assigned index values\n')

# Reclassify Flammap heat per unit area raster values based on Table 2 of Rodriguez y Silva et al
# Convert heat per unit area raster from kJ/m2 to kcal/m2
HUA_kcal <- HUA*0.2388459
HUA_rtab <- data.frame(Low=c(-1,380,1265,1415,1610,1905,2190,4500,6630,8000),
                       High=c(380,1265,1415,1610,1905,2190,4500,6630,8000,1000000),
					   Val=c(1,2,3,4,5,6,7,8,9,10))
HUA_recl <- classify(HUA_kcal,HUA_rtab)
cat('Reclassified heat per unit area to assigned index values\n')

# Calculate Energy Behavior raster
# The focal statistics smoothing operation replaces the previous propFuel multplication factor
eb_rast <- (2*FL_recl*HUA_recl)/(FL_recl + HUA_recl)
cat('Created energy behavior raster\n')
# Smooth, but set non-burnable pixels to value of 1
# This matches the Pyrologix "NbrCircle(56.41896, 'MAP')" methods
wm <- focalMat(DEM,60,'circle') # Create weights matrix to define neighborhood
wm[wm>0] <- 1; wm[wm==0] <- NA
eb_smooth <- focal(eb_rast,w=wm,fun='mean',na.policy='all',na.rm=T)
eb_smooth[FBFM40 < 100] <- 1
writeRaster(eb_smooth,paste0(opath,'/',oname,'_energy_behavior.tif'),overwrite=T)
cat('Smoothed energy behavior raster\n\n')

#-> ACCESSIBILITY SUB INDEX
# Run Euclidean Distance on the streets layer
# Classify areas within 100 m of a road as a 10, decreasing to 1 where the nearest road is >900 m
cat('#-> Calculating accessibility sub index\n')
streets_lyr <- vect(streets)
streets_rast <- rasterize(streets_lyr,DEM,touches=T)
roaddist <- distance(streets_rast)
cat('Calculated distance to roads\n')
ACC_rtab <- data.frame(Low=c(-1,100,200,300,400,500,600,700,800,900),
                       High=c(100,200,300,400,500,600,700,800,900,1000000),
					   Val=c(10,9,8,7,6,5,4,3,2,1))
accessibility <- classify(roaddist,ACC_rtab)
writeRaster(accessibility,paste0(opath,'/',oname,'_accessibility.tif'),overwrite=T)
cat('Reclassified distance to roads raster to accessibility index value\n\n')

#-> MOBILITY SUB-INDEX
# Skipping, setting mobility = 1 in calculation

#-> PENETRABILITY SUB-INDEX
# Alaska edits:
# 1. drop aspect 
# 2. customize fuel control values
cat('#-> Calculating penetrability sub index\n')

# Calculate Percent Slope
slope_deg <- terrain(DEM,v='slope',neighbors=8,unit='degrees')
slope <- tan((pi/180)*slope_deg)*100
slope[slope > 500] <- 500
cat('Calculated percent slope\n')

# Convert Percent Slope Raster into assigned values
SLP_rtab <- data.frame(Low=c(0,6,11,16,21,26,31,36,41,46),
                       High=c(6,11,16,21,26,31,36,41,46,1000),
					   Val=c(10,9,8,7,6,5,4,3,2,1))
slope_class <- classify(slope,SLP_rtab,include.lowest=T) # Python uses lowest
cat('Converted percent slope into assigned values\n')

# Reclassify fuels into assigned RTC values using reclass table
# This will convert the fuel types into a fuel control difficulty weight originally based on 
# estimated fireline production rates for a 20-person crew. This table was modified for use
# in Alaska by their fire modeling and analysis committee. 
RTC_rtab <- data.frame(FBFM40=c(-9999,
                                91,92,93,98,99,
                                101,102,103,104,105,106,107,108,109,
                                121,122,123,124,
								141,142,143,144,145,146,147,148,149,
								161,162,163,164,165,
								181,182,183,184,185,186,187,188,189,
								201,202,203,204),
                       Val=c(NA,
					         NA,NA,NA,NA,NA,
					         10,8,8,7,8,10,3,3,3,
					         10,9,7,3,
                             5,6,5,4,4,4,4,2,2,
                             8,6,3,2,2,
                             8,8,8,4,6,7,2,8,8,
                             3,2,1,1))
fuel_cntrl <- classify(FBFM40,RTC_rtab)
cat('Reclassified fuels into assigned fireline production values\n')

# Load trails layer
trails_lyr <- vect(trails)
cat('Loaded trails layer\n')

# Calculate length of trails within 1 ha (56.41896m radius) moving window.
# Two options for the calculations are provided. The first is the closest to the ESRI calculations 
# (line length within a circle around the cell center), but it is very slow for large study areas.
# The second option is a reasonable approximation suitable for use with large study areas.
#trail_len <- rasterizeGeom(trails_lyr,DEM,fun='length',unit='m') # Option 1
#trail_len_perha <- focal(trail_len,w=wm,fun='sum',na.policy='all',na.rm=T)
trail_len <- rasterize(trails_lyr,DEM,touches=T,field=30,background=0) # Option 2
trail_len_perha <- focal(trail_len,w=wm,fun='sum',na.policy='all',na.rm=T)
writeRaster(trail_len_perha,paste0(opath,'/',oname,'_trail_len_perha.tif'),overwrite=T)
cat('Calculated length of trails within 1 ha\n')

# Convert pre-suppression trails raster into assigned values
TRL_rtab <- data.frame(Low=c(0,10,20,30,40,50,60,70,80,90),
                       High=c(10,20,30,40,50,60,70,80,90,10000),
					   Val=c(1,2,3,4,5,6,7,8,9,10))
trail_len_class <- classify(trail_len_perha,TRL_rtab,include.lowest=T) # Python uses lowest
cat('Reclassified road/trail lengths to assigned values\n')

# Create raster for Penetrability Sub Index
# The focal statistics smoothing operation replaces the previous propFuel multplication factor
# Alaska edit: divided by 4 instead of 5 because aspect component dropped
penetrability <- (slope_class + fuel_cntrl + 2*trail_len_class)/4
cat('Calculated penetrability index\n')
pen_smooth <- focal(penetrability,w=wm,fun='mean',na.policy='all',na.rm=T)
writeRaster(pen_smooth,paste0(opath,'/',oname,'_penetrability.tif'),overwrite=T)
cat('Smoothed penetrability index\n\n')

#-> FIRELINE OPENING/CREATION SUB-INDEX
cat('#-> Calculating fireline opening/creation sub index\n')

# Create Slope Adjustment Raster for hand work
HND_rtab <- data.frame(Low=c(0,16,31,46),
                       High=c(16,31,46,1000),
					   Val=c(1,0.8,0.6,0.5))
slope_adjust_hand <- classify(slope,HND_rtab,include.lowest=T) # Python uses lowest
cat('Calculated slope adjustment raster for hand work\n')

# Create Slope Adjustment Raster for machine work
MCN_rtab <- data.frame(Low=c(0,16,21,26,31,36),
                       High=c(16,21,26,31,36,1000),
					   Val=c(1,0.8,0.7,0.6,0.5,0))	   
slope_adjust_mach <- classify(slope,MCN_rtab,include.lowest=T) # Python uses lowest
cat('Calculated slope adjustment raster for machine work\n')

# Calculate Fireline Opening Sub Index
firelineopening <- fuel_cntrl*slope_adjust_hand + fuel_cntrl*slope_adjust_mach
writeRaster(firelineopening,paste0(opath,'/',oname,'_firelineopening.tif'),overwrite=T)
cat('Calculated fireline opening index\n\n')

#-> Kit O'Connor (KO) UNIVERSAL SLOPE MOBILITY HAZARD ADJUSTMENTS
cat('#-> Calculating KO universal slope mobility hazard adjustment\n')

# Create slope hazard raster
SHZ_rtab <- data.frame(Low=c(0,10,20,30,40,50),
                       High=c(10,20,30,40,50,1000),
					   Val=c(0.01,0.1,0.2,0.4,0.6,0.8))	   
slope_haz <- classify(slope,SHZ_rtab,include.lowest=T) # Python uses lowest
cat('Calculated slope mobility hazard\n\n')

#-> SDI CALCULATION
cat('#-> Calculating SDI from sub-indices\n')

# Create SDI raster without slope modification
# mobility = 1, 2020 publication version
SDI <- eb_smooth/(accessibility + 1 + pen_smooth + firelineopening)
SDI[is.na(SDI) & (FBFM40 < 100)] <- 0 # Setting non-burnable pixels to zero
cat('Calculated SDI without slope modification\n')

# Optional: save denominator for inspection
sdi_denom <- accessibility + 1 + pen_smooth + firelineopening
writeRaster(sdi_denom,paste0(opath,'/',oname,'_denominator.tif'),overwrite=T)

# Create SDI raster with universal slope modification
SDI_SLP <- SDI + slope_haz
cat('Calculated SDI with slope modification\n')

# Set water to zero SDI
SDI_SLP_H2O <- ifel(FBFM40==98,0,SDI_SLP)
cat('Set water to zero SDI\n')
    
# Combine using base model for roads & trails, and slope modification for everything else
trails_rds <- rbind(streets_lyr,trails_lyr)
trails_rds_rast <- rasterize(trails_rds,DEM,touches=T)
fSDI <- ifel(is.na(trails_rds_rast),SDI_SLP_H2O,SDI)
cat('Combined base and slope adjusted SDI for roaded/trailed or unroaded/untrailed areas\n\n')

# Convert to integer
# This is a true round, but ESRI rounds down - use floor if it is necessary to match
fSDI100 <- round(fSDI*100,0)
cat('Converted to integer\n')

# Add attribute table with default symbology
colcodes <- c('#3DA1D1','#9AC5B4','#D7EA90','#F6D865','#F88F3E','#F0261C')
rgbvals <- col2rgb(colcodes)
colbreaks <- c(-1,10,20,40,70,100,350)
classnames <- c('0-10 (lowest difficulty)','10-20','20-40','40-70','70-100',
                '100-350 (highest difficulty)')
rat <- data.frame(Value=as.numeric(unlist(unique(fSDI100))),Classname=NA,
                  Red=NA,Green=NA,Blue=NA,Alpha=255)
for(i in 1:nrow(rat)){
	cind <- cut(rat$Value[i],colbreaks,labels=F)
	rat$Classname[i] <- classnames[cind]
	rat$Red[i] <- rgbvals[1,cind]
	rat$Green[i] <- rgbvals[2,cind]
	rat$Blue[i] <- rgbvals[3,cind]
}
levels(fSDI100)[[1]] <- rat
cat('Added raster attribute table with default symbology\n')

# Save raster	
writeRaster(fSDI100,paste0(opath,'/',oname,'_SDI.tif'),overwrite=T)
cat('Saved raster\n\n')

############################################END ANALYSIS############################################

####################################################################################################
cat('\nFinished at: ',as.character(Sys.time()),'\n\n',sep='')
############################################END MESSAGE#############################################
