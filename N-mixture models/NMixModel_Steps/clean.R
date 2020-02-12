# -----------------------
# Supplementary Material for Deith and Brodie 2020; "Predicting defaunation - accurately mapping bushmeat hunting pressure over large areas"
#    doi: 10.1098/rspb.2019-2677
# -----------------------

# clean.R:
#   Clean and modify datasets from load.R for analysis of circuit-theoretic map
#   The output of this script is an unmarked data frame that can be used for p-count Nmix model fitting
#   First created by Mairin Deith on Aug 26, 2017
#   Last modified by Mairin Deith on Oct 5, 2019

# Clean environment ------------------------------------------------------------

rm(list = ls())
setwd("~/Documents/GradSchool/PhD/Research/CircuitTheory_Borneo/PRSB_Revision2/Nmix_Models/Code/")

# Load Libraries and data ------------------------------------------------------

library(tidyverse)
library(sf)
library(sp)
library(raster)
library(stars)
library(viridis)
library(unmarked)

load("loadedData.RData")

# Camera trap station data -----------------------------------------------------
colnames(humanCap) <- c('rowID','station','year','month','day','hour','minute')

humanCap$timeStamp <- ISOdatetime(humanCap$year, humanCap$month, humanCap$day, humanCap$hour, humanCap$minute, 0)
humanCap$date <- ISOdate(humanCap$year, humanCap$month, humanCap$day)
humanCap$timeDiff_sec <- NA

# Identify independent captures ------------------------------------------------

for (i in 1:nrow(humanCap)){
  if(i==1){
    # First row of observations: set arbitrary time since last capture at 2000 seconds
    humanCap$timeDiff_sec[i]=2000
  } else {
    if(humanCap$station[i]==humanCap$station[i-1]){
      # If the station is the same as the last in the list, calculate the time between adjacent captures
      humanCap$timeDiff_sec[i]=as.numeric(difftime(humanCap$timeStamp[i],humanCap$timeStamp[i-1],units='secs'))
    } else {
      # If the station is unique,is not the same, set the arbitrary 2000 limits
      humanCap$timeDiff_sec[i]=2000
    }
  }
}

# Select only observations which are more than 30 minutes apart
humanCap$Independent <- FALSE
humanCap$Independent[which(humanCap$timeDiff_sec>=1800)] = TRUE

# indHumanCap contains only independent data
indHumanCap <- humanCap[which(humanCap$timeDiff_sec>=1800),]
indHumanCap$station <- as.character(indHumanCap$station)
nrow(humanCap) # 19488 total captures
nrow(indHumanCap) # 1589 independent observations

# Summary of human captures
sumHumanCap <- data.frame(humanCap %>%
                            group_by(station, Independent) %>%
                            summarize(humans = length(unique(timeStamp)))
)

# Only captures that are >1800s apart (30 minutes)
sumHumanCap_sub <- data.frame(indHumanCap %>%
                                group_by(station) %>%
                                summarize(humans = length(unique(timeStamp)))
)

# Identify independent sampling days -------------------------------------------
colnames(ctEffort_o) <- c('station','year','month','day','hours')
# Remove effort data from Danum, where there no detection data are available
# Note: these data were used for another study and human captures were not
#   available.

ctEffort_o <- ctEffort_o %>%
  dplyr::filter(!str_detect(station, "ag"))

ctEffort_o$date <- ISOdate(ctEffort_o$year,ctEffort_o$month, ctEffort_o$day)

# List of all camera traps that "exerted" some effort
camTraps <- as.character(unique(ctEffort_o$station))

longestEffort <- ctEffort_o %>%
  group_by(station) %>%
  summarize(n.trapping.days = length(unique(date)))

maxDays <- max(longestEffort$n.trapping.days)
minDays <- min(longestEffort$n.trapping.days)

# Day-by-day matrix of per-camera trap image collection for all sites (not only the 76 that observed hunters)
# Observations:
y.array <- as.data.frame(matrix(nrow=length(camTraps), ncol = max(maxDays)))
for(col in 1:ncol(y.array)){
  colnames(y.array)[col] <- paste0('day',col)
}
rownames(y.array) <- camTraps

# Trapping effort:
camTrapHours <- as.data.frame(matrix(nrow=length(camTraps), ncol = max(maxDays)))
for(col in 1:ncol(camTrapHours)){
  colnames(camTrapHours)[col] <- paste0('day',col)
}

rownames(camTrapHours) <- camTraps

# Seasons as detection covariate -----------------------------------------------
# NE Monsoon - November to February (extending to Oct-March)
# SW Monsoon - May to August (extending to April-Sept)

# 1 = NE monsoon
# 0 = SW monsoon

season <- as.data.frame(matrix(nrow=length(camTraps), ncol = max(maxDays)))
for(col in 1:ncol(season)){
  colnames(season)[col] <- paste0('day',col)
}
rownames(season) <- camTraps

ne.monsoon <- c(1:3, 10:12)
sw.monsoon <- c(4:9)

# Populate capture y-array with detection covars -------------------------------
for(i in 1:nrow(y.array)){
  # Create subset dataframes for each camera trap station
  camTrapEffortSub <- ctEffort_o[which(ctEffort_o$station==camTraps[[i]]),]
  capturesSub <- indHumanCap[which(indHumanCap$station==camTraps[[i]]),]
  # How many days of active capture time were recorded?
  capturesLength <- length(unique(camTrapEffortSub$date))
  # Fill in temporary '0' for day-by-day captures
  y.array[i,1:capturesLength] <- 0
  camTrapHours[i,1:nrow(camTrapEffortSub)] <- camTrapEffortSub$hours
  season.tmp <- (camTrapEffortSub$month %in% ne.monsoon)*1
  season[i,1:length(season.tmp)] <- season.tmp
  for(h in 1:nrow(camTrapEffortSub)){
    # Some stations have >1 camera, and therefore apparently have 2x the effort.
    # This should be ignored.
    if(camTrapHours[i,h]>24){
      camTrapHours[i,h] <- 24
    }
  }
  # Skip filling in the y.array if there are no human captures at that station
  if(nrow(capturesSub)==0){
    next
  } else {
    # For each human image, add '1' to the date that matches
    for(c in 1:nrow(capturesSub)){
      captureDateInteger <- which(camTrapEffortSub$date==capturesSub$date[c])
      y.array[i,captureDateInteger] <- y.array[i,captureDateInteger]+1
    }
  }
}

# Remove ctCovar data from Danum, where we do not have records of human captures

ctCovars <- ctCovars_in %>%
  dplyr::select(-starts_with("field")) %>%
  dplyr::filter(!str_detect(station, "ag")) %>%
  # Only include camera traps that were not stolen or malfunctioning
  dplyr::filter(station %in% row.names(y.array))

# Spatial points dataframe of camera trap locations and local covariates
ctDF <- SpatialPointsDataFrame(coords=ctCovars[,c('longitude','latitude')],
                               proj4string = CRS("+proj=longlat +ellps=WGS84"),
                               data=ctCovars
)

# Trail as a detection covariate -----------------------------------------------

trailObs <- camTrapHours
for(i in 1:nrow(trailObs)){
  trail.tmp <- as.numeric(ctDF@data$trail[i])-1
  trail.obs.tmp <- as.numeric(trailObs[i,])*0+trail.tmp
  trailObs[i,] <- as.factor(trail.obs.tmp)
}

# Extract accessibility values at each site  -----------------------------------
# Raster Extract: Accessibility measures  --------------------------------------
ctDF@data$CTaccess <- raster::extract(x=CTaccess,y=SpatialPoints(ctDF@coords), buffer=90, fun=mean, na.rm=T)

# Raster Extract: Others --------------------------------------------------
GRUMPS_urbanRural <- '/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/urban_rural_access_50k/acc_50k.tif'
URaccess.global <- raster(GRUMPS_urbanRural, colname='URAccess', proj4string = CRS("+proj=longlat +ellps=WGS84"))
URaccess <- crop(URaccess.global, extent(109.5379, 119.2997, 0.855001, 7.380556))

rm(URaccess.global)

ctDF@data$URAccessibility <- raster::extract(x=URaccess,y=SpatialPoints(ctDF@coords),buffer=90, fun=mean, na.rm=T)
ctDF@data$roadDist <- raster::extract(x=distanceRast,y=SpatialPoints(ctDF@coords),buffer=90, fun=mean, na.rm=T)*111.320 # convert to km

popdistance_VS_Rast <- raster("~/Documents/GradSchool/PhD/Research/CircuitTheory_Borneo/PRSB_Revision2/distanceToPopOrPlantation/gazetteer_VS_growdist_plantation_m_Oct4.asc")
crs(popplantdistance_VS_Rast) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ctDF@data$popPlant_VS_Dist <- raster::extract(x=popplantdistance_VS_Rast,y=SpatialPoints(ctDF@coords),buffer=90,fun=mean, na.rm=T)/1000 # convert to km

# Transform accessibility ------------------------------------------------------
# load("CTDF_Data_Oct1_2019.RData")
ctDF@data$logCTaccess <- log10(ctDF@data$CTaccess+0.00001)

ctDF@data$logURAccess <- log10(ctDF@data$URAccessibility+0.00001)
ctDF@data$logRDAccess <- log10(ctDF@data$roadDist+0.00001)
ctDF@data$logPPlantD_VS_Access <- log10(ctDF@data$popPlant_VS_Dist+0.00001)

# Standardize covariates (mean=0, SD=1)-----------------------------------------
# These are factors
not.standardized <- c('trail', 'state', 'protection', 'logged')
standardize.this <- colnames(ctDF@data)[c(7,11:ncol(ctDF@data))]

standardize.index <- ctDF@data[, standardize.this]
no.standardize.index <- ctDF@data[, not.standardized]

standardizedCovars <- cbind(no.standardize.index,
                            apply(X = standardize.index, MARGIN = 2, FUN = scale)
)

rawCovars <- cbind(no.standardize.index,
                   standardize.index)

# unmarkedFrame creation -------------------------------------------------------

model.uframe <- unmarkedFramePCount(y=y.array,
                                    siteCovs= standardizedCovars,
                                    obsCovs = list(hours=camTrapHours,
                                                   season=season,
                                                   trail=trailObs))

# Create unstandardized model.uframe -------------------------------------------
model.uframe.raw <- unmarkedFramePCount(y=y.array,
                                        siteCovs= rawCovars,
                                        obsCovs = list(hours=camTrapHours,
                                                       season=season,
                                                       trail=trailObs))

save(list("model.uframe", "model.uframe.raw"), "model.uframe.data.RData")
