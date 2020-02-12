# Proposed R structure:

# 1: load.R - load original data
# 2: clean.R - data pre-processing
# 3: functions.R - all functions needed to run the analysis
# 4: main.R - this uses source() to call the above and any other analysis scripts

setwd("~/Documents/GradSchool/Research/CircuitTheory_Borneo/GFlowMapAnalysis/Code/PRSB_Revision2/")

# Load helper functions
source("functions.R", echo = F, print.eval = F)

# Load original data
source("load.R", echo = F, print.eval = F)

# Clean data and prepare unmarked dataframe for model fitting
# May take some time due to raster extraction
source("clean.R", echo = F, print.eval = F)

### Only has to be run once
# source("exploration.R")

# Should be run on cluster, this script takes a very long time
#   to fit
source("modelfitting.R")

#   Plot model predictions from fitted N-mixture models
source("plot.R")
