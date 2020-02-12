# -----------------------
# Supplementary Material for Deith and Brodie 2020; "Predicting defaunation - accurately mapping bushmeat hunting pressure over large areas"
#    doi: 10.1098/rspb.2019-2677
# -----------------------

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
