# load required libraries
require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)
library(fastDummies)
library(correlation)

# source macro pyrome evaluation function
source("fun_hp_macro_pyro.R")

# read modeling data
db <- read.csv("../04modelingDataBirds/data_for_modeling.csv")

species <- names(db)[3:296]

env_vars <- db[, c("bio1", "bio2", "bio4", "bio8", "bio12", "bio15", "elevation",
                   "prop_trees", "prop_shrubs", "prop_grasslands", 
                   "prop_crops", "prop_built", "prop_bare_soil", 
                   "prop_water", "prop_wetlands", "macro_pyro", "samp_bias")]

env_vars$macro_pyro <- as.factor(env_vars$macro_pyro)

# calculate correlations including categorical macro pyrome variable
results <- correlation(env_vars, include_factors = TRUE, method = "auto")
filtered_results <- results %>% filter(str_detect(Parameter1, "macro_pyro") | str_detect(Parameter2, "macro_pyro"))

write.csv(data.frame(filtered_results), "tab_cor_predictors_macro_pyro.csv", row.names = FALSE)

# create dummy variables and remove baseline macro_pyro M0
env_vars_dummy <- dummy_cols(env_vars, 
                             select_columns = "macro_pyro", 
                             remove_selected_columns = TRUE)

env_vars_dummy <- env_vars_dummy[, !names(env_vars_dummy) %in% "macro_pyro_M0"]

# iterate through species and run evaluation
for(i in 1:length(species)) {

    current_spp <- species[i]
    
    blocks_h <- read.csv(paste0("./blocks/blocksCV_", current_spp, ".csv"))
    db2 <- cbind(db, blocks_h)
    db2 <- na.omit(db2)
  
    bartHO_Eq(
        records = db2[, current_spp],
        explainVars = env_vars_dummy,
        folds = db2[, c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")],
        name = current_spp,
        ntree = c(200),
        k = c(2),
        power = c(2),
        base = c(0.95)
    )       
}