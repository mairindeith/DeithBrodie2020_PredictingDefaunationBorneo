# -----------------------
# Supplementary Material for Deith and Brodie 2020; "Predicting defaunation - accurately mapping bushmeat hunting pressure over large areas"
#    doi: 10.1098/rspb.2019-2677
# -----------------------

# load.R:
#   Load all datasets required for analysis of circuit-theoretic map
#   First created by Mairin Deith on Aug 26, 2017
#   Last modified by Mairin Deith on Oct 4, 2019

#  Clean environment --------------------------------------------------------

rm(list = ls())
setwd("~/Documents/GradSchool/Research/CircuitTheory_Borneo/")

#  Libraries --------------------------------------------------------
library(raster)
library(sf)
library(tidyverse)
library(stars)

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

