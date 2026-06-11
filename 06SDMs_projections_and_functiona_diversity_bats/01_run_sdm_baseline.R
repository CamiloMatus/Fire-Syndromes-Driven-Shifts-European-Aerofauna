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
db <- read.csv("../04modelingDataBats/data_for_modeling.csv") |> na.omit()
species_list <- names(db)[1:31]

# generate predictions by species
for(spp in species_list) {
  
    names_ev <- c("bio1", "bio2", "bio4", "bio8", "bio12", "bio15", "elevation",
                  "prop_trees", "prop_shrubs", "prop_grasslands",
                  "prop_crops", "prop_built", "prop_bare_soil",
                  "prop_water", "prop_wetlands", "samp_bias")

    message(paste("working on", spp, "-", which(species_list == spp), "of", length(species_list)))

    if(sum(db[, spp]) < 250) { nrep <- 10 } else { nrep <- 1 }
  
    for (replication in 1:nrep) {
        
        # load random pseudo-absences used in model evaluation
        rp_path <- paste0("../05modelEvaluationBats/random_pseudo_absences_baseline/randomPseudoAbs_Replication_", replication, "_", spp, ".csv")
        rp <- read.csv(rp_path)
        
        dat_1 <- db[row.names(db) %in% rp$x, ]
        dat_2 <- db[db[, spp] == 1, ]
        dat_h <- rbind(dat_1, dat_2)

        # filter out variables with only zero values
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
                # samprate equal to 0 for sampling bias mitigation
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
    
    # compute optimal threshold
    th <- optimal.thresholds(d_h, opt.methods = 3, threshold = 1000)[2] |> as.numeric()
    
    d_h$bin <- 0
    names(d_h)[1] <- "id"
    names(d_h)[2] <- "obs_record"
    names(d_h)[3] <- paste("prob", spp, sep = ".")
    names(d_h)[4] <- paste("bin", spp, sep = ".")
    
    d_h[d_h[, 3] >= th, 4] <- 1
    d_h$x <- db$x
    d_h$y <- db$y
    
    out_file <- paste0("./predictions/predictions_baseline_", spp, ".csv")
    write.csv(d_h, file = out_file, row.names = FALSE)
}