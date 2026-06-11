# load required libraries
library(fastDummies)
require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)

source("fun_hp_baseline.R")

# read modeling data
db <- read.csv("../4modelingDataBats/data_for_modeling.csv")

# assuming first 31 columns correspond to the species
species <- names(db)[1:31]

# baseline predictors (uncorrelated, pearson < 0.6)
predictors <- c(
  "bio1", "bio2", "bio4", "bio8", "bio12", "bio15",
  "prop_trees", "prop_shrubs", "prop_grasslands",
  "prop_crops", "prop_built", "prop_bare_soil",
  "prop_water", "prop_wetlands", "elevation",    
  "samp_bias"
)

# iterate through species and run evaluation
for(i in 1:length(species)) {

    current_spp <- species[i]
    message(paste("testing hyperparameters combinations for species", current_spp, "....", i, "of", length(species)))
    
    blocks_h <- read.csv(paste0("./blocks/blocksCV_", current_spp, ".csv"))
    
    # bind data and block validation information
    db2 <- cbind(db[, current_spp, drop=FALSE], db[, predictors], blocks_h)
    db2 <- na.omit(db2)
  
    # run hyperparameter evaluation testing values that avoid overfitting
    bartHO_Eq(
        records = db2[, current_spp],
        explainVars = db2[, predictors],
        folds = db2[, c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")],
        name = current_spp,
        ntree = c(200),
        k = c(2),
        power = c(2),
        base = c(0.95)
    )     
}