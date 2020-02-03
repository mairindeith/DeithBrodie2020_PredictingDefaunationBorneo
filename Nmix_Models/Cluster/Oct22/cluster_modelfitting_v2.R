# Code to fit data in a prepared unmarked data frame to 
#   validate various accessibility maps including the circuit-theoretic approach
#
# Code written by M. Deith beginning Aug 26, 2017
# Last edited by M. Deith on Feb 26, 2019

### Differences between this and the 2017 version of the script:
#   Remove all data and manipulation of mammal captures - human only
#   Revision from line 183 to include season as a detection covariate
#   Revision from line 300 to include trail as a detection covariate
#   Remove model fitting and creation of iterative models in new script, '03.02.2019_ProcBModelFitting.R'

### Begin script
#   Clear the environment

rm(list = ls())
setwd("~/GFlowMapAnalysis/PRSB_Revision2/Oct22")

# Initialize parallel -----------------------------------------------------
req.libraries <- c("foreach", "doParallel")

for(rl in c(req.libraries)){
  if(rl %in% rownames(installed.packages())==FALSE){
    # For the cluster, manually install VGAM
    install.packages(rl,
                     lib = "../GFlowMapAnalysis",
                     repos = "https://cloud.r-project.org")
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

req.libraries <- c("devtools","unmarked", "AICcmodavg", #"dplyr",
                   "foreach", "doParallel")

for(rl in c(req.libraries, "VGAM")){
  if(rl %in% rownames(installed.packages())==FALSE){
    # For the cluster, manually install VGAM
    if(rl=="VGAM"){
      install_version("VGAM", version = "1.0-3", repos = "https://cloud.r-project.org")
    } else {
      install.packages(rl,
                       lib = "../GFlowMapAnalysis",
                       repos = "https://cloud.r-project.org")
    }
  } else {
    print(paste("Already installed:",rl))
  }
  library(rl, character.only = T)
}

# Load data frame ---------------------------------------------------------
load('Oct.22.19_modeluframes_v2.RData')
# load('14.03.2019.UnmarkedDataframe.Rdata')
print('...dataframe loaded')

# for-loop to automatically generate model formulae -----------------------

keep.covars <- c('protection', 'logged', 'elevation')

# uct = untransformed CT access
presence.covars.uct.allp <- c('CTaccess_GazAll_pass', keep.covars)
presence.covars.uct.alli <- c('CTaccess_GazAll_impass', keep.covars)

presence.covars.uct.vsp <- c('CTaccess_GazVS_pass', keep.covars)
presence.covars.uct.vsi <- c('CTaccess_GazVS_impass', keep.covars)

presence.covars.uct.vsop <- c('CTaccess_GazVSO_pass', keep.covars)
presence.covars.uct.vsoi <- c('CTaccess_GazVSO_impass', keep.covars)

presence.covars.uct.wallp <- c('CTaccess_Walk_GazAll_pass', keep.covars)
presence.covars.uct.walli <- c('CTaccess_Walk_GazAll_impass', keep.covars)

presence.covars.uct.wvsp <- c('CTaccess_Walk_GazVS_pass', keep.covars)
presence.covars.uct.wvsi <- c('CTaccess_Walk_GazVS_impass', keep.covars)

presence.covars.uct.wvsop <- c('CTaccess_Walk_GazVSO_pass', keep.covars)
presence.covars.uct.wvsoi <- c('CTaccess_Walk_GazVSO_impass', keep.covars)

# umd = untransformed minDist
presence.covars.umd <- c('minDist', keep.covars)
# uur = untransformed urban-rural
presence.covars.uur <- c('URAccessibility', keep.covars)
# upd = untransformed distance to the nearest village
presence.covars.upd <- c('popDist', keep.covars)
# uppd = untransformed distance to the nearest village OR plantation
presence.covars.uppd <- c('popPlantDist', keep.covars)

presence.covars.upd.vs <- c('pop_VS_Dist', keep.covars)
# uppd = untransformed distance to the nearest village OR plantation
presence.covars.uppd.vs <- c('popPlant_VS_Dist', keep.covars)


# uct = untransformed CT access
presence.covars.tct.allp <- c('logCTaccess_GazAll_pass', keep.covars)
presence.covars.tct.alli <- c('logCTaccess_GazAll_impass', keep.covars)

presence.covars.tct.vsp <- c('logCTaccess_GazVS_pass', keep.covars)
presence.covars.tct.vsi <- c('logCTaccess_GazVS_impass', keep.covars)

presence.covars.tct.vsop <- c('logCTaccess_GazVSO_pass', keep.covars)
presence.covars.tct.vsoi <- c('logCTaccess_GazVSO_impass', keep.covars)

presence.covars.tct.wallp <- c('logCTaccess_Walk_GazAll_pass', keep.covars)
presence.covars.tct.walli <- c('logCTaccess_Walk_GazAll_impass', keep.covars)

presence.covars.tct.wvsp <- c('logCTaccess_Walk_GazVS_pass', keep.covars)
presence.covars.tct.wvsi <- c('logCTaccess_Walk_GazVS_impass', keep.covars)

presence.covars.tct.wvsop <- c('logCTaccess_Walk_GazVSO_pass', keep.covars)
presence.covars.tct.wvsoi <- c('logCTaccess_Walk_GazVSO_impass', keep.covars)

# umd = untransformed minDist
presence.covars.tmd <- c('logMDAccess', keep.covars)
# uur = untransformed urban-rural
presence.covars.tur <- c('logURAccess', keep.covars)
# upd = untransformed distance to the nearest village
presence.covars.tpd <- c('logPDAccess', keep.covars)
# uppd = untransformed distance to the nearest village OR plantation
presence.covars.tppd <- c('logPPlantDAccess', keep.covars)
# upd = untransformed distance to the nearest village
presence.covars.tpd.vs <- c('logPD_VS_Access', keep.covars)
# uppd = untransformed distance to the nearest village OR plantation
presence.covars.tppd.vs <- c('logPPlantD_VS_Access', keep.covars)

# hours is not significant, so only include these observation covariates:
obs.covars <- c('season','trail')

presence.covars <- list(
  # presence.covars.uct.allp,
  # presence.covars.uct.alli,
  presence.covars.uct.vsp,
  presence.covars.uct.vsi,
  presence.covars.uct.vsop,
  presence.covars.uct.vsoi,
  # presence.covars.uct.wallp,
  # presence.covars.uct.walli,
  # presence.covars.uct.wvsp,
  # presence.covars.uct.wvsi,
  # presence.covars.uct.wvsop,
  # presence.covars.uct.wvsoi,
  # presence.covars.tct.allp, # Transformed
  # presence.covars.tct.alli,
  presence.covars.tct.vsp,
  presence.covars.tct.vsi,
  presence.covars.tct.vsop,
  presence.covars.tct.vsoi,
  # presence.covars.tct.wallp,
  # presence.covars.tct.walli,
  # presence.covars.tct.wvsp,
  # presence.covars.tct.wvsi,
  # presence.covars.tct.wvsop,
  # presence.covars.tct.wvsoi,
  # presence.covars.uct.old, # Old CT maps
  presence.covars.uur,
  presence.covars.umd, 
  presence.covars.upd, 
  # presence.covars.uppd, #,
  presence.covars.upd.vs, 
  presence.covars.uppd.vs, #,
  # presence.covars.tct.old,
  presence.covars.tur,
  presence.covars.tmd, 
  # presence.covars.tpd, 
  # presence.covars.tppd,
  presence.covars.tpd.vs, 
  presence.covars.tppd.vs
)

# model.count <- 1
# model.list <- list()

today <- format(Sys.Date(), "%b_%d_%y")

for(access in 1:length(presence.covars)){
    model.count <- 1
    model.list <- list()
    print(access)  
    presence.covars.tmp <- presence.covars[[access]]
    for(i in 0:length(presence.covars.tmp)){
        # Include cases with zero presence covariates 
        if(i==0){   
            next
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
        model.count <- model.count+1
        }
     }
  print("Model list generated")
  formula.list <- lapply(model.list, formula)
  fit.model.names <- list()
  fit.model.list <- foreach(f=1:length(formula.list), .packages="unmarked") %dopar% { 
    fitmodel.tmp <- pcount(formula.list[[f]], model.uframe, K=100) #, K=250)
    model.name.tmp <- paste0("m",f,paste0(as.vector(gsub("~","",formula.list[[f]])),collapse="."),sep=".")
    print(paste("Model fit success:",model.name.tmp))
    assign(model.name.tmp, fitmodel.tmp)
    fit.model.names <- append(fit.model.names, model.name.tmp)
    fitmodel.tmp
    }
    saveRDS(fit.model.list, file = paste0(today, "_", access, "_fitlist.RData"))
    # AICc and BICc tables:
    aicPoisson <- aictab(fit.model.list)
    aicPoissonDF <- as.data.frame(aicPoisson)
    rownames(aicPoissonDF) <- 1:nrow(aicPoissonDF)
    aicPoissonDF$ModelID <- substring(aicPoissonDF$Modnames, first=4)
    aicPoissonDF$ranking <- NaN
    aicPoissonDF$formula <- NaN

    for(i in 1:nrow(aicPoissonDF)){
      modID <- as.numeric(aicPoissonDF$ModelID[i])
      aicPoissonDF$ranking[i] <- rownames(aicPoissonDF)[i]
      aicPoissonDF$formula[i] <- paste0("m",modID,paste0(as.vector(gsub("~","",fit.model.list[[modID]]@formula)),
                                                 collapse="."),
                                    sep="")
    }
    saveRDS(aicPoissonDF, file = paste0(today, "_", access, "_aicTable.RData"))
    write.csv(aicPoissonDF, file = paste0(today, "_", access, "_aicTable.csv"))

    # BIC ranking -------------------------------------------------------------

    bicPoisson <- bictab(fit.model.list)
    bicPoissonDF <- as.data.frame(bicPoisson)
    rownames(bicPoissonDF) <- 1:nrow(bicPoissonDF)
    bicPoissonDF$ModelID <- substring(bicPoissonDF$Modnames, first=4)
    bicPoissonDF$ranking <- NaN
    bicPoissonDF$formula <- NaN
    
    for(i in 1:nrow(bicPoissonDF)){
      modID <- as.numeric(bicPoissonDF$ModelID[i])
      bicPoissonDF$ranking[i] <- rownames(bicPoissonDF)[i]
      bicPoissonDF$formula[i] <- paste0("m",modID,paste0(as.vector(gsub("~","",fit.model.list[[modID]]@formula)),
                                                     collapse="."),
                                    sep="")
    }
    saveRDS(bicPoissonDF, file = paste0(today, "_", access, "_bicTable.RData"))
    write.csv(bicPoissonDF, file = paste0(today, "_", access, "_bicTable.csv"))
    rm(fit.model.list)
    gc()
}

print("COMPLETE")
