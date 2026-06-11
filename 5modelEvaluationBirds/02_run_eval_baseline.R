# load required libraries
require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)

# source baseline evaluation function
source("fun_hp_baseline.R")

# read modeling data
db <- read.csv("../4modelingDataBirds/data_for_modeling.csv")

species <- names(db)[3:296]

# update predictors to match standard english names
predictors <- c("bio1", "bio2", "bio4", "bio8", "bio12", "bio15", "elevation",
                "prop_trees", "prop_shrubs", "prop_grasslands", 
                "prop_crops", "prop_built", "prop_bare_soil", 
                "prop_water", "prop_wetlands", "samp_bias")

env_vars <- db[, predictors]

# calculate correlation matrix for baseline predictors
tab_cor_final <- round(cor(env_vars[, predictors]), 2)
tab_cor_final_df <- data.frame(tab_cor_final)
write.csv(tab_cor_final_df, "tab_cor_predictors_baseline.csv")

# iterate through species and run evaluation
for(i in 1:length(species)) {
    
    current_spp <- species[i]
    
    blocks_h <- read.csv(paste0("./blocks/blocksCV_", current_spp, ".csv"))
    db2 <- cbind(db, blocks_h)
    db2 <- na.omit(db2)
    
    # testing fixed parameters for baseline reference
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