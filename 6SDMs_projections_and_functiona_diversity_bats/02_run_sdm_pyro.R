# load required libraries
library(fastDummies)
require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)

only_zeros <- function(x) {
    all(x == 0, na.rm = TRUE)
}

dir.create("./predictions", showWarnings = FALSE)

# load modeling data
db <- read.csv("../4modelingDataBats/data_for_modeling.csv") |> na.omit()
species_list <- names(db)[1:31]

db$Pyrome_Consenso <- as.factor(db$Pyrome_Consenso)

env_vars_dummy <- dummy_cols(db, 
                             select_columns = "Pyrome_Consenso", 
                             remove_selected_columns = TRUE)

# update main db with dummy pyromes
db <- env_vars_dummy

# generate predictions by species
for(spp in species_list) {

    names_ev <- c("bio1", "bio2", "bio4", "bio8", "bio12", "bio15", "elevation",
                  "prop_trees", "prop_shrubs", "prop_grasslands",
                  "prop_crops", "prop_built", "prop_bare_soil",
                  "prop_water", "prop_wetlands",
                  "Pyrome_Consenso_1", "Pyrome_Consenso_2", "Pyrome_Consenso_3", 
                  "Pyrome_Consenso_4", "Pyrome_Consenso_5", "Pyrome_Consenso_6", 
                  "Pyrome_Consenso_7", "Pyrome_Consenso_8", "Pyrome_Consenso_9", 
                  "Pyrome_Consenso_10", "Pyrome_Consenso_11", "samp_bias")

    message(paste("working on", spp, "-", which(species_list == spp), "of", length(species_list)))

    if(sum(db[, spp]) < 250) { nrep <- 10 } else { nrep <- 1 }
  
    for (replication in 1:nrep) {
        
        rp_path <- paste0("../05modelEvaluationBats/random_pseudo_absences_pyro/randomPseudoAbs_Replication_", replication, "_", spp, ".csv")
        rp <- read.csv(rp_path)
        
        dat_1 <- db[row.names(db) %in% rp$x, ]
        dat_2 <- db[db[, spp] == 1, ]
        dat_h <- rbind(dat_1, dat_2)

        cols_all_zeros <- names(which(sapply(dat_h, only_zeros)))
        dat_h <- dat_h[, !names(dat_h) %in% cols_all_zeros]
        names_ev <- names_ev[!names_ev %in% cols_all_zeros]
        
        mod <- bart2(
            formula = dat_h[, names_ev],
            k = 2, power = 2, base = 0.95, n.tree = 200, n.chains = 1,
            data = dat_h[, spp], 
            keepTrees = TRUE, verbose = FALSE
        )

        s_list <- list() 
        for (e in 1:length(names_ev)) {
            r_h <- raster(nrow = nrow(db), ncol = 1)
            if(names_ev[e] != "samp_bias") {
                values(r_h) <- db[, names_ev[e]]
            } else {
                values(r_h) <- rep(0, nrow(db))
            }
            s_list[[e]] <- r_h 
        } 

        s_stack <- stack(s_list)
        names(s_stack) <- names_ev
        
        pred1 <- predict2.bart(object = mod, x.layers = s_stack, quiet = TRUE) |> values() |> as.numeric()
        gc()
        
        if(replication == 1) { pred_h <- pred1 } else { pred_h <- cbind(pred_h, pred1) }
    }

    if(nrep == 1) { pred_mean <- pred_h } else {
        pred_mean <- rowMeans(pred_h, na.rm = TRUE)
    }

    d_h <- data.frame(1:nrow(db), db[, spp], pred_mean)
    th <- optimal.thresholds(d_h, opt.methods = 3, threshold = 1000)[2] |> as.numeric()
    
    d_h$bin <- 0
    names(d_h)[1] <- "id"
    names(d_h)[2] <- "obs_record"
    names(d_h)[3] <- paste("prob", spp, sep = ".")
    names(d_h)[4] <- paste("bin", spp, sep = ".")
    
    d_h[d_h[, 3] >= th, 4] <- 1
    d_h$x <- db$x
    d_h$y <- db$y
    
    out_file <- paste0("./predictions/predictions_pyro_", spp, ".csv")
    write.csv(d_h, file = out_file, row.names = FALSE)
}

# correction for estimates = 1 in pyromes where the species had no observations
for(u in species_list) {

    d_con <- db[db[, u] == 1, c(u, "Pyrome_Consenso_1", "Pyrome_Consenso_2",
                                "Pyrome_Consenso_3", "Pyrome_Consenso_4", "Pyrome_Consenso_5",
                                "Pyrome_Consenso_6", "Pyrome_Consenso_7", "Pyrome_Consenso_8",
                                "Pyrome_Consenso_9", "Pyrome_Consenso_10", "Pyrome_Consenso_11")]

    if(TRUE %in% sapply(d_con, only_zeros) == TRUE) {
        
        d_h <- read.csv(paste0("./predictions/predictions_pyro_", u, ".csv"))
        write.csv(d_h, paste0("./predictions/predictions_pyro_backup_", u, ".csv"), row.names = FALSE)
        
        vars_h <- names(which(sapply(d_con, only_zeros)))
        
        for(v in 1:length(vars_h)) {
            d_h2 <- cbind(d_h, db[, vars_h[v]])
            d_h[row.names(d_h2[d_h2[, 7] == 1, ]), 4] <- 0
        }
        
        write.csv(d_h, file = paste0("./predictions/predictions_pyro_", u, ".csv"), row.names = FALSE)
    }
}