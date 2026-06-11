# load required libraries
library(fastDummies)
require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)

# helper function to detect columns with only zeros
only_zeros <- function(x) {
    all(x == 0, na.rm = TRUE)
}

# create output directory if it does not exist
dir.create("./var_imp", showWarnings = FALSE)

# load modeling data
db <- read.csv("../04modelingDataBirds/data_for_modeling.csv") |> na.omit()
species_list <- names(db)[3:296]

# select environment variables and convert macro_pyro to factor
base_predictors <- c("bio1", "bio2", "bio4", "bio8", "bio12", "bio15", "elevation",
                     "prop_trees", "prop_shrubs", "prop_grasslands",
                     "prop_crops", "prop_built", "prop_bare_soil",
                     "prop_water", "prop_wetlands", "macro_pyro", "samp_bias")

env_vars <- db[, base_predictors]
env_vars$macro_pyro <- as.factor(env_vars$macro_pyro)

# create dummy variables for macro pyromes
env_vars_dummy <- dummy_cols(env_vars, 
                             select_columns = "macro_pyro", 
                             remove_selected_columns = TRUE)

db <- cbind(db, env_vars_dummy[, c("macro_pyro_M1", "macro_pyro_M2", 
                                   "macro_pyro_M3", "macro_pyro_M4", "macro_pyro_M5")])

# evaluate variable importance per species
for (spp in species_list) {
    
    names_ev <- c("bio1", "bio2", "bio4", "bio8", "bio12", "bio15", "elevation",
                  "prop_trees", "prop_shrubs", "prop_grasslands",
                  "prop_crops", "prop_built", "prop_bare_soil",
                  "prop_water", "prop_wetlands", "macro_pyro_M1",
                  "macro_pyro_M2", "macro_pyro_M3", "macro_pyro_M4", 
                  "macro_pyro_M5", "samp_bias")

    message(paste("working on", spp, "-", which(species_list == spp), "of", length(species_list)))
    
    # condition for species with few records
    if(sum(db[, spp]) < 150) { nrep <- 10 } else { nrep <- 1 }
    
    for (replication in 1:nrep) {
        
        # load random pseudo-absences used in model evaluation
        rp_path <- paste0("../05modelEvaluationBirds/random_pseudo_absences_macro_pyro/randomPseudoAbs_Replication_", replication, "_", spp, ".csv")
        rp <- read.csv(rp_path)
        
        # subset data for model fitting
        dat_absences <- db[row.names(db) %in% rp$x, ]
        dat_presences <- db[db[, spp] == 1, ]
        dat_h <- rbind(dat_absences, dat_presences)

        # clean columns that contain only zeros to prevent bart failures
        cols_all_zeros <- names(which(sapply(dat_h, only_zeros)))
        dat_h <- dat_h[, !names(dat_h) %in% cols_all_zeros]
        names_ev <- names_ev[!names_ev %in% cols_all_zeros]
        
        # evaluate variable importance
        vi_h <- varimp.diag(dat_h[, names_ev], dat_h[, spp], iter = 2)
        
        out_filename <- paste0("./var_imp/var_imp_macro_pyro_", spp, "_rep_", replication, ".csv")
        write.csv(vi_h$data, file = out_filename, row.names = FALSE)
    }
}