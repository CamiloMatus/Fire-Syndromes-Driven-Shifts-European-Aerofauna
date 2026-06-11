# load required libraries
library(fastDummies)
require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)

source("fun_hp_pyro.R")

# read modeling data
db <- read.csv("../04modelingDataBats/data_for_modeling.csv")

species <- names(db)[1:31]

db$Pyrome_Consenso <- as.factor(db$Pyrome_Consenso)

# generate dummy variables for pyromes
env_vars_dummy <- dummy_cols(db[, c("Pyrome_Consenso")], 
                             select_columns = "Pyrome_Consenso", 
                             remove_selected_columns = TRUE)

# combine predictors list
predictors <- c(
    "bio1", "bio2", "bio4", "bio8", "bio12", "bio15",
    "prop_trees", "prop_shrubs", "prop_grasslands",
    "prop_crops", "prop_built", "prop_bare_soil",
    "prop_water", "prop_wetlands", "elevation",    
    "Pyrome_Consenso_1", "Pyrome_Consenso_2", "Pyrome_Consenso_3", 
    "Pyrome_Consenso_4", "Pyrome_Consenso_5", "Pyrome_Consenso_6", 
    "Pyrome_Consenso_7", "Pyrome_Consenso_8", "Pyrome_Consenso_9", 
    "Pyrome_Consenso_10", "Pyrome_Consenso_11", "samp_bias"
)

# bind everything properly
full_db <- cbind(db, env_vars_dummy)

# iterate through species and run evaluation
for(i in 1:length(species)) {

    current_spp <- species[i]
    message(paste("testing hyperparameters combinations for species", current_spp, "....", i, "of", length(species)))
    
    blocks_h <- read.csv(paste0("./blocks/blocksCV_", current_spp, ".csv"))
    
    db2 <- cbind(full_db[, current_spp, drop=FALSE], full_db[, predictors], blocks_h)
    db2 <- na.omit(db2)
  
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