# -----------------------
# Supplementary Material for Deith and Brodie 2020; "Predicting defaunation - accurately mapping bushmeat hunting pressure over large areas"
#    doi: 10.1098/rspb.2019-2677
# -----------------------

### workflow.R:
# This Rscript prepares data for and fits N-mixture models of hunter site usage
#   Contains 3 steps of analysis - data loading, cleaning/preparation, and model fitting - into one script.
# Created by Mairin Deith in Spring 2018 for manuscript:
#   Deith and Brodie 'Predicting defaunation – accurately mapping bushmeat hunting effort over large areas'

### Note:
# Data included in the Supplementary Materials of this publication are
# already formatted for the creation of an unmarkedFrame. This Rscript
# provides code to incorporate raw data into an N-mixture model format.
# The siteCovs, obsCovs, and ObsToY objects are presented in different sheets of
# the file "DeithBrodie_ProcRSB_NMixtureModelData.xlsx"
# These can be used to directly create the model.uframe object with:
#
#     model.uframe <- unmarkedFramePCount(
#                                    obsToY= [OBS TO Y MATRIX],
#                                    siteCovs = [SITE SPECIFIC COVARIATESS DATAFRAME],
#                                    obsCovs = [OBSERVATION COVARIATES DATAFRAME]
#                       )

### Step 1: Load data

# Previously in load.R script:
#   Load all datasets required for analysis of circuit-theoretic map
#   First created by Mairin Deith on Aug 26, 2017
#   Last modified by Mairin Deith on August 12, 2019

#  Clean environment -----------------------------------------------------------
rm(list = ls())

#  Load libraries for all stages -----------------------------------------------
library(raster)
library(sf)
library(tidyverse)
library(stars)
library(sp)
library(viridis)
library(unmarked)

# Camera trap station data -----------------------------------------------------
# Images of hunters:
humanCapturesFile <- file.path("CTData/HumanEncounters.csv")
humanCap <- read_csv(humanCapturesFile, col_types = cols())

# Environmental covars
ctCovarsFile <- file.path("CTData/Updated_site_covariates_protected.csv")
ctCovars_in <- read_csv(ctCovarsFile, col_types = cols(
  trail = col_factor(),
  logged = col_factor(),
  protection = col_factor()
))

# Trapping effort (duration) at each station
ctEffortFile <- file.path("CTData/Updated sampling effort.csv")
ctEffort_in <- read_csv(ctEffortFile, col_types = cols())

ctEffort_o <- ctEffort_in %>%
  arrange(station)

# Spatial data -----------------------------------------------------------------

### Administrative data
state.poly.sabswk <- read_sf("InputSpatialData/AdminBoundaries/MSY_adm/MSY_Borneo_Merged.shp")
state.poly.sabah <- read_sf("InputSpatialData/AdminBoundaries/MSY_adm/MSY_Borneo.shp")
state.poly.swk <- read_sf("InputSpatialData/AdminBoundaries/MSY_adm/MSY_Borneo_swk.shp")

### Accessibility map files
# CTaccess
base_ct_path <- '/home/mairin/Documents/GradSchool/Research/CircuitTheory_Borneo/PRSB_Revision2/Revised_CTMapCreation/GFlowOutputs/ClusterFlexOutputs/AllMaps'
ct_f <- file.path(base_ct_path, 'GrassSummation_VillageSabah_Impassable.asc')
CTaccess <- raster(ct_f, colname='CTaccess', proj4string = CRS("+proj=longlat +ellps=WGS84"))

### URaccess
GRUMPS_urbanRural <- '/home/mairin/Documents/GradSchool/PhD/Research/CircuitTheory_Borneo/urban_rural_access_50k/acc_50k.tif'
URaccess.global <- raster(GRUMPS_urbanRural, colname='URAccess', proj4string = CRS("+proj=longlat +ellps=WGS84"))
URaccess <- crop(URaccess.global, extent(109.5379, 119.2997, 0.855001, 7.380556))

rm(URaccess.global)

### Distance to the nearest road
distance_file <- '/home/mairin/Documents/GradSchool/PhD/Research/CircuitTheory_Borneo/InputSpatialData/BorneoRoadDensity/clippedRoadDistance.tif'
distanceRast <- raster(distance_file,colname='roadDist', proj4string = CRS("+proj=longlat +ellps=WGS84"))

### Distance to EITHER nearest village OR plantation
pop_plant_file_f <- '/home/mairin/Documents/GradSchool/PhD/Research/CircuitTheory_Borneo/PRSB_Revision2/distanceToPopOrPlantation/rgrowdist_rasterizedGaveauPlantationsGazetteerSettlements_km.tif'
popplantdistanceRast <- raster(pop_plant_file_f, colname='popPlantDist', proj4string = CRS("+proj=longlat +ellps=WGS84"))

### Save intermediate outputs
save.image("loadedData.RData")

### Step 2: Clean data
#   Clean and modify loaded datasets to create data frames in Nmix model format
#   The output of this section is an unmarked data frame that can be used for p-count Nmix model fitting

### If starting from loaded data, uncomment this:
# load("loadedData.RData")

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

### Step 3: Generate list of candidate model formulae and fit N-mix models
# modelfitting.R
# Code to fit data in a prepared unmarked data frame to
#   validate various accessibility maps including the circuit-theoretic approach
#   Designed to be run on a multi-core cluster with parallel computing abilities
# First created by Mairin Deith on Aug 26, 2017
# Last modified by Mairin Deith on Oct 21, 2019

# Initialize parallel and load libraries in each core --------------------------
req.libraries <- c("foreach", "doParallel")
for(rl in c(req.libraries)){
  if(rl %in% rownames(installed.packages())==FALSE){
    # For the cluster, manually install VGAM
    install.packages(rl, repos = "https://cloud.r-project.org")
  } else {
    print(paste("Already installed:",rl))
  }
  library(rl, character.only = T)
}

cores <- detectCores()
cl <- makeCluster(cores[1]/2)
registerDoParallel(cl)

print('...cluster registered.')

#  Import libraries --------------------------------------------------------
# Check to see which libraries are installed already:

req.libraries <- c("devtools","unmarked", "AICcmodavg")

for(rl in c(req.libraries, "VGAM")){
  if(rl %in% rownames(installed.packages())==FALSE){
    # For the cluster, manually install VGAM
    if(rl=="VGAM"){
      install_version("VGAM", version = "1.0-3")
    } else {
      install.packages(rl, repos = "https://cloud.r-project.org")
    }
  } else {
    print(paste("Already installed:",rl))
  }
  library(rl, character.only = T)
}

### If need to load data, uncomment:
# load("~/Documents/GradSchool/PhD/Research/CircuitTheory_Borneo/PRSB_Revision2/Nmix_Models/Code/Oct.21.19cleanedData.RData")

# for-loop to automatically generate model formulae -----------------------

### Non-accessibility covariates
keep.covars <- c('protection', 'logged', 'elevation')

### Accessibility covariates

# tct = transformed CT access
presence.covars.tct <- c('logCTAccess', keep.covars)
# trd = transformed roadDist
presence.covars.trd <- c('logMDAccess', keep.covars)
# tur = transformed urban-rural
presence.covars.tur <- c('logURAccess', keep.covars)
# uppd = untransformed distance to the nearest village OR plantation
presence.covars.uppd <- c('popPlant_VS_Dist', keep.covars)

### Polynomial models
presence.covars.poly.tct <- c('logCTAccess + I(logCTAccess^2)', keep.covars)
presence.covars.poly.trd <- c('logMDAccess + I(logMDAccess^2)', keep.covars)
presence.covars.poly.tur <- c('logURAccess + I(logURAccess^2)', keep.covars)
presence.covars.poly.uppd <- c('popPlant_VS_Dist + I(popPlant_VS_Dist^2)', keep.covars)

# Observation covariates - hours not included as this was not significant
obs.covars <- c('season','trail')

# Let's pull in the old CT covariates here as well:
presence.covars <- list(
                        presence.covars.tct,
                        presence.covars.tur,
                        presence.covars.trd,
                        presence.covars.uppd,
                        presence.covars.poly.tct,
                        presence.covars.poly.trd,
                        presence.covars.poly.tur,
                        presence.covars.poly.uppd
)

for(access in 1:length(presence.covars)){
  access_covar <- presence.covars[[access]][1]
  model.list <- list()

  # Linear models (no interacting terms)

  presence.covars.tmp <- presence.covars[[access]]
  for(i in 0:length(presence.covars.tmp)){
    # Include cases with zero presence covariates
    if(i==0){
      pres.sub.tmp <- '1'
    } else {
      pres.sub.tmp <- apply(combn(presence.covars.tmp, i), 2, paste0, collapse="+")
    }
    for(j in 0:length(obs.covars)){
      # Include cases with zero observation covariates
      if(j==0){
        obs.sub.tmp <- '1'
      } else {
        obs.sub.tmp <- apply(combn(obs.covars, j), 2, paste0, collapse="+")
      }
      model.list <- append(model.list, paste0("~",obs.sub.tmp," ~",pres.sub.tmp))
    }
  model.list.u <- unique(model.list)
  formula.list <- lapply(unique(model.list.u), formula)
  fit.model.names <- list(model.list.u)

  fit.model.list <- foreach(f=1:length(formula.list), .packages="unmarked") %dopar% {
      fitmodel.tmp <- pcount(formula.list[[f]], model.uframe, K=50)
      model.name.tmp <- paste0("f",f,".pcount.k50")
      fitmodel.tmp
    }
    # names(fit.model.list) <- fit.model.names
    # save("fit.model.names", file = paste0(today, "_", access_covar, "_untransformed_modelNames.Rdata"))
    # save("fit.model.list", file = paste0(today, "_", access_covar, "_untransformed_modelFitList.Rdata"))

    # AIC ranking -------------------------------------------------------------

    aicPoisson <- aictab(fit.model.list)
    aicPoissonDF <- as.data.frame(aicPoisson)
    rownames(aicPoissonDF) <- 1:nrow(aicPoissonDF)
    aicPoissonDF$ModelID <- as.numeric(gsub(aicPoissonDF$Modnames,
                                            pattern = "Mod", replacement = ""))
    aicPoissonDF$formula <- NaN
    for(i in 1:nrow(aicPoissonDF)){
      aicPoissonDF$formula[[i]] <- as.character(fit.model.names[[1]][[aicPoissonDF$ModelID[[i]]]])
    }
    ### Save model rankings for that accessibility measure
    write.csv(aicPoissonDF, paste0(today, "_", access_covar, "_UntransformedModelFittingResults_AIC.csv"))
    # BIC ranking -------------------------------------------------------------
    bicPoisson <- bictab(fit.model.list)
    bicPoissonDF <- as.data.frame(bicPoisson)
    rownames(bicPoissonDF) <- 1:nrow(bicPoissonDF)
    bicPoissonDF$ModelID <- as.numeric(gsub(bicPoissonDF$Modnames,
                                            pattern = "Mod", replacement = ""))
    bicPoissonDF$formula <- NaN
    for(i in 1:nrow(bicPoissonDF)){
      bicPoissonDF$formula[[i]] <- as.character(fit.model.names[[1]][[bicPoissonDF$ModelID[[i]]]])
    }
    ### Save model rankings for that accessibility measure
    write.csv(bicPoissonDF, paste0(today, "_", access_covar, "_UntransformedModelFittingResults_BIC.csv"))
  }
  stopCluster(cl)
}
