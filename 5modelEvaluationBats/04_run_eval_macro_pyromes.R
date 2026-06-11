# load required libraries
library(fastDummies)
require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)

source("fun_hp_macro_pyro.R")

# read modeling data
db <- read.csv("../4modelingDataBats/data_for_modeling.csv")

species <- names(db)[1:31]

db$macro_pyro <- as.factor(db$macro_pyro)

# generate dummy variables for macro pyromes
env_vars_dummy <- dummy_cols(db[, c("macro_pyro")], 
                             select_columns = "macro_pyro", 
                             remove_selected_columns = TRUE)

# combine predictors list
predictors <- c(
    "bio1", "bio2", "bio4", "bio8", "bio12", "bio15",
    "prop_trees", "prop_shrubs", "prop_grasslands",
    "prop_crops", "prop_built", "prop_bare_soil",
    "prop_water", "prop_wetlands", "elevation",    
    "macro_pyro_M1", "macro_pyro_M2", "macro_pyro_M3", 
    "macro_pyro_M4", "macro_pyro_M5", "samp_bias"
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