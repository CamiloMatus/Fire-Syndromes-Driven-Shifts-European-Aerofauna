bartHO_Eq <- function(records, explainVars, folds, ntree=200, k=2, power=2, base=0.95, name=NULL) {

  only_zeros <- function(x) { all(x == 0, na.rm = TRUE) }
  set.seed(123)
  if(sum(records[records == 1]) < 250) { nrep <- 10 } else { nrep <- 1 }

  if(!dir.exists("random_pseudo_absences_pyro")) dir.create("random_pseudo_absences_pyro")
  if(!dir.exists("hyper_params_eval_pyro")) dir.create("hyper_params_eval_pyro")
    
  if(is.null(name)) { name <- "speciesX" }  

  HyptoTest <- tidyr::expand_grid(ntree, k, power, base) %>% data.frame()

  dat.h00 <- data.frame(records, explainVars, folds)
  for (r in 1:nrep) {
    rp.out <- c()
    for (f in 1:ncol(folds)) {
      n.p <- nrow(dat.h00[dat.h00$records == 1 & dat.h00[, ((ncol(explainVars) + 1) + f)] == FALSE, ])
      rp.h <- sample(rownames(dat.h00[dat.h00$records == 0 & dat.h00[, ((ncol(explainVars) + 1) + f)] == FALSE, ]), size = n.p)
      rp.out <- c(rp.out, rp.h)
    }
    assign(paste0("randomPoint_Replication_", r), rp.out)
    write.csv(get(paste0("randomPoint_Replication_", r)),
              file = paste0("random_pseudo_absences_pyro/randomPseudoAbs_Replication_", r, "_", name, ".csv"))  
  }

  for (h in 1:nrow(HyptoTest)) {
    auc <- c(); boyce <- c(); tss <- c(); kappa <- c(); sens <- c(); spec <- c()

    for (run in 1:ncol(folds)) {  
      message(paste("Testing hyperparameter in Fold", run))
      
      dat.h0 <- data.frame(records, explainVars, folds[,run])    

      if(nrep != 1) {
        for (replication in 1:nrep) {
          d.h <- dat.h0[row.names(dat.h0) %in% get(paste0("randomPoint_Replication_", replication)), ]
          if(replication == 1) { dat.h.RP.meanPred <- d.h } else { dat.h.RP.meanPred <- rbind(dat.h.RP.meanPred, d.h) }
        }
      }
          
      pred.h <- c()
      for (replication in 1:nrep) {
        dat.h.0 <- rbind(dat.h0[dat.h0$records == 1, ], dat.h0[row.names(dat.h0) %in% get(paste0("randomPoint_Replication_", replication)), ])
        
        if(nrow(dat.h.0[dat.h.0[, "folds...run."] == FALSE, ]) != 0 | nrow(dat.h.0[dat.h.0[, "folds...run."] == FALSE, ]) == 1) {
          if(!any(sapply(dat.h.0[dat.h.0[, "folds...run."] == TRUE, names(explainVars)], only_zeros))) {

            mod <- bart(x.train = dat.h.0[dat.h.0[, "folds...run."] == TRUE, names(explainVars)], 
                        y.train = dat.h.0[dat.h.0[, "folds...run."] == TRUE, "records"],
                        k = as.numeric(paste(HyptoTest[h,"k"])), power = as.numeric(paste(HyptoTest[h,"power"])),
                        base = as.numeric(paste(HyptoTest[h,"base"])), ntree = as.numeric(paste(HyptoTest[h,"ntree"])),
                        nchain = 1, keeptrees = TRUE, verbose = FALSE)
           
            if(nrep == 1) {
              s <- list() 
              for (e in 1:length(explainVars)) {
                r.h <- raster(nrow = nrow(dat.h.0[dat.h.0[, "folds...run."] == FALSE, ]), ncol = 1)
                values(r.h) <- dat.h.0[dat.h.0[, "folds...run."] == FALSE, names(explainVars)[e]]
                s[[e]] <- r.h
              } 
              s <- stack(s); names(s) <- names(explainVars)
            } else {
              s <- list() 
              for (e in 1:length(explainVars)) {
                dat.h <- rbind(dat.h0[dat.h0$records == 1, ], dat.h.RP.meanPred)
                r.h <- raster(nrow = nrow(dat.h[dat.h[, "folds...run."] == FALSE, ]), ncol = 1)
                values(r.h) <- dat.h[dat.h[, "folds...run."] == FALSE, names(explainVars)[e]]
                s[[e]] <- r.h
              } 
              s <- stack(s); names(s) <- names(explainVars)
            }
            
            s[["samp_bias"]] <- 0
            pred <- predict2.bart(object = mod, x.layers = s) %>% values() %>% as.numeric()
            if(replication == 1) { pred.h <- pred } else { pred.h <- cbind(pred.h, pred) }
          }
        }
      }
        
      if(nrow(dat.h.0[dat.h.0[, "folds...run."] == FALSE, ]) != 0 | nrow(dat.h.0[dat.h.0[, "folds...run."] == FALSE, ]) == 2) {
        if (length(pred.h) == 0) {
          auc.h <- NA; boyce.h <- NA; tss.h <- NA; kappa.h <- NA; sens.h <- NA; spec.h <- NA
        } else {
          if(nrep == 1) { pred.out <- pred.h } else { pred.out <- rowMeans(pred.h, na.rm = TRUE) }
          
          if(nrep != 1) {
             obs.h <- dat.h[dat.h[, "folds...run."] == FALSE, "records"]
          } else {
             dat.h <- rbind(dat.h0[dat.h0$records == 1, ], dat.h0[row.names(dat.h0) %in% get(paste0("randomPoint_Replication_", replication)), ])
             obs.h <- dat.h[dat.h[, "folds...run."] == FALSE, "records"]
          }
             
          auc.h <- auc(DATA = data.frame(1:length(obs.h), obs.h, pred.out))[1] %>% as.numeric()
          dat.h3 <- data.frame(obs.h, pred.out)
          boyce.h <- ecospat.boyce(fit = pred.out, obs = dat.h3[dat.h3[, 1] == 1, 2], PEplot = FALSE)$cor %>% as.numeric()
          
          if(length(unique(obs.h)) > 1 && !all(is.na(pred.out))) {
            thresh_seq <- seq(min(pred.out, na.rm=TRUE), max(pred.out, na.rm=TRUE), length.out = 100)
            tss_vec <- numeric(100)
            for(t_idx in 1:100) {
              t_val <- thresh_seq[t_idx]
              tp <- sum(obs.h == 1 & pred.out >= t_val, na.rm=TRUE); tn <- sum(obs.h == 0 & pred.out < t_val, na.rm=TRUE)
              fp <- sum(obs.h == 0 & pred.out >= t_val, na.rm=TRUE); fn <- sum(obs.h == 1 & pred.out < t_val, na.rm=TRUE)
              sens_t <- tp / (tp + fn); spec_t <- tn / (tn + fp)
              tss_vec[t_idx] <- (sens_t + spec_t) - 1
            }
            opt_t_idx <- which.max(tss_vec)
            tss.h <- tss_vec[opt_t_idx]; t_val <- thresh_seq[opt_t_idx]
            tp <- sum(obs.h == 1 & pred.out >= t_val, na.rm=TRUE); tn <- sum(obs.h == 0 & pred.out < t_val, na.rm=TRUE)
            fp <- sum(obs.h == 0 & pred.out >= t_val, na.rm=TRUE); fn <- sum(obs.h == 1 & pred.out < t_val, na.rm=TRUE)
            sens.h <- tp / (tp + fn); spec.h <- tn / (tn + fp)
            po <- (tp + tn) / (tp + tn + fp + fn)
            pe <- ((tp + fn) * (tp + fp) + (fp + tn) * (fn + tn)) / ((tp + tn + fp + fn)^2)
            kappa.h <- (po - pe) / (1 - pe)
          } else {
            tss.h <- NA; kappa.h <- NA; sens.h <- NA; spec.h <- NA
          }
        } 
      } else {
        auc.h <- NA; boyce.h <- NA; tss.h <- NA; kappa.h <- NA; sens.h <- NA; spec.h <- NA
      }
      auc <- c(auc, auc.h); boyce <- c(boyce, boyce.h); tss <- c(tss, tss.h)
      kappa <- c(kappa, kappa.h); sens <- c(sens, sens.h); spec <- c(spec, spec.h)
    }   
    
    dat.out <- data.frame(auc, boyce, tss, kappa, sens, spec)
    dat.out$ntree <- HyptoTest[h, "ntree"]; dat.out$k <- HyptoTest[h, "k"]
    dat.out$power <- HyptoTest[h, "power"]; dat.out$base <- HyptoTest[h, "base"]
    dat.out$species <- name; dat.out$Fold <- 1:ncol(folds)
    
    if(h == 1) { dat.out.f <- dat.out } else { dat.out.f <- rbind(dat.out.f, dat.out) }
  }
  
  dat.out.f2 <- dat.out.f
  dat.out.f2$hp <- paste(dat.out.f2$ntree, dat.out.f2$k, dat.out.f2$power, dat.out.f2$base)
  
  BOYCE <- tapply(dat.out.f2$boyce, dat.out.f2$hp, mean, na.rm = TRUE)
  AUC <- tapply(dat.out.f2$auc, dat.out.f2$hp, mean, na.rm = TRUE)
  TSS <- tapply(dat.out.f2$tss, dat.out.f2$hp, mean, na.rm = TRUE)
  KAPPA <- tapply(dat.out.f2$kappa, dat.out.f2$hp, mean, na.rm = TRUE)
  SENS <- tapply(dat.out.f2$sens, dat.out.f2$hp, mean, na.rm = TRUE)
  SPEC <- tapply(dat.out.f2$spec, dat.out.f2$hp, mean, na.rm = TRUE)
  
  out <- data.frame(BOYCE, AUC, TSS, KAPPA, SENS, SPEC)
  out$hp <- rownames(out)
  
  out2 <- merge(x = out, y = dat.out.f2[, c("ntree", "k", "power", "base", "hp")], by = "hp", all.y = FALSE)
  out3 <- distinct(out2, .keep_all = TRUE)
  
  write.csv(out3, file = paste0("hyper_params_eval_pyro/hyperEval_", name, "_pyro.csv"), row.names = FALSE)
}