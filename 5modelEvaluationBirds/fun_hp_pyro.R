# this function is used to test hyperparameters of the bart algorithm. 
# the same number of pseudo absences (random) and presences are used 
# as recommended for machine learning techniques based on decision trees. 
# also, mean predictions are made based on 10 repetitions 
# (each repetition with random pseudo absences), 
# in order to achieve the best predictive performances as recommended.
# both recommendations come from the study of:
# https://doi.org/10.1111/j.2041-210X.2011.00172.x

bartHO_Eq <- function(records, explainVars, folds, ntree=200, k=2, power=2, base=0.95, name=NULL) {

  only_zeros <- function(x) {
    all(x == 0, na.rm = TRUE)
  }

  set.seed(123)

  if(sum(records[records==1]) < 150) { nrep <- 10 } else { nrep <- 1 }

  if(file.exists("random_pseudo_absences_pyro") == FALSE) {
    dir.create(file.path(getwd(), "random_pseudo_absences_pyro"))
  }

  if(file.exists("hyper_params_eval_pyro") == FALSE) {
    dir.create(file.path(getwd(), "hyper_params_eval_pyro"))
  }
    
  if(is.null(name) == TRUE) { name <- "speciesX" }  

  HyptoTest <- tidyr::expand_grid(ntree, k, power, base) %>% data.frame()

  dat.h00 <- data.frame(records, explainVars, folds)
  for (r in 1:nrep) {
    rp.out <- c()
    for (f in 1:ncol(folds)) {
      
      n.p <- nrow(dat.h00[dat.h00$records==1 & dat.h00[,((ncol(explainVars)+1)+f)]==FALSE,])
      rp.h <- sample(rownames(dat.h00[dat.h00$records==0 & dat.h00[,((ncol(explainVars)+1)+f)]==FALSE,]),
                     size = n.p)
      rp.out <- c(rp.out, rp.h)
    }
    
    assign(paste("randomPoint_Replication_", r, sep = ""), rp.out)
    write.csv(get(paste("randomPoint_Replication_", r, sep = "")),
              file = paste("random_pseudo_absences_pyro/randomPseudoAbs_Replication_", r, "_", name, ".csv", sep = ""))  
  }

  for (h in 1:nrow(HyptoTest)) {
    auc <- c()
    boyce <- c()  
    tss <- c()
    kappa <- c()
    sens <- c()
    spec <- c()

    for (run in 1:ncol(folds)) {  
      message(paste("Testing hyperparameter in Fold", run))
      message(paste("ntrees =", HyptoTest[h,"ntree"],
                    "; k =", HyptoTest[h,"k"], "; base =",
                    HyptoTest[h,"base"], "; power =",
                    HyptoTest[h,"power"],
                    "...", "combination", h, "of", nrow(HyptoTest), sep=" "))

      dat.h0 <- data.frame(records, explainVars, folds[,run])    

      if(nrep != 1) {
        for (replication in 1:nrep) {
          d.h <- dat.h0[row.names(dat.h0) %in% get(paste("randomPoint_Replication_", replication, sep = "")),]
          if(replication == 1) { dat.h.RP.meanPred <- d.h } else { dat.h.RP.meanPred <- rbind(dat.h.RP.meanPred, d.h) }
        }
      }
          
      pred.h <- c()
      for (replication in 1:nrep) {
        if(nrep != 1) {
          message(paste("Replication", replication, "of", nrep))
        } else { NULL }

        dat.h.0 <- rbind(dat.h0[dat.h0$records==1,],
                         dat.h0[row.names(dat.h0) %in% get(paste("randomPoint_Replication_", replication, sep = "")),])
        
        if(nrow(dat.h.0[dat.h.0[,"folds...run."]==FALSE,]) != 0 | nrow(dat.h.0[dat.h.0[,"folds...run."]==FALSE,]) == 1) {
        
          if(any(sapply(dat.h.0[dat.h.0[,"folds...run."]==TRUE, names(explainVars)], only_zeros)) != TRUE) {

            mod <- bart(x.train = dat.h.0[dat.h.0[,"folds...run."]==TRUE, names(explainVars)], 
                        y.train = dat.h.0[dat.h.0[,"folds...run."]==TRUE, "records"],
                        k = paste(HyptoTest[h,"k"]) %>% as.numeric(),
                        power = paste(HyptoTest[h,"power"]) %>% as.numeric(),
                        base = paste(HyptoTest[h,"base"]) %>% as.numeric(),
                        ntree = paste(HyptoTest[h,"ntree"]) %>% as.numeric(),
                        nchain = 1, 
                        keeptrees = TRUE, verbose = F)
           
            message("var contribution")
            print(varimp(mod))

            if(nrep == 1) {
              s <- list() 
              for (e in 1:length(explainVars)) {
                r.h <- raster(nrow=nrow(dat.h.0[dat.h.0[,"folds...run."]==FALSE,]), ncol=1)
                values(r.h) <- dat.h.0[dat.h.0[,"folds...run."]==FALSE, names(explainVars)[e]]
                s[[e]] <- r.h
              } 
              s <- stack(s)
              names(s) <- names(explainVars)
            } else {
            
              s <- list() 
              for (e in 1:length(explainVars)) {
                dat.h <- rbind(dat.h0[dat.h0$records==1,],
                               dat.h.RP.meanPred)
                r.h <- raster(nrow=nrow(dat.h[dat.h[,"folds...run."]==FALSE,]), ncol=1)
                values(r.h) <- dat.h[dat.h[,"folds...run."]==FALSE, names(explainVars)[e]]
                s[[e]] <- r.h
              } 
              s <- stack(s)
              names(s) <- names(explainVars)
            }
            
            s[["samp_bias"]] <- 0

            pred <- predict2.bart(object = mod, x.layers = s) %>% values() %>% as.numeric()
            
            if(replication == 1) { pred.h <- pred } else { pred.h <- cbind(pred.h, pred) }

          } else { NULL }  
        } else { NULL }
      }
        
      if(nrow(dat.h.0[dat.h.0[,"folds...run."]==FALSE,]) != 0 | nrow(dat.h.0[dat.h.0[,"folds...run."]==FALSE,]) == 2) {
        
        if (length(pred.h) == 0) {
          auc.h <- NA; boyce.h <- NA; tss.h <- NA; kappa.h <- NA; sens.h <- NA; spec.h <- NA
        } else {
          
          if(nrep == 1) { pred.out <- pred.h } else {
            pred.out <- rowMeans(pred.h, na.rm = T)
          }
          
          if(nrep != 1) {
             obs.h <- dat.h[dat.h[,"folds...run."]==FALSE, "records"]
          } else {
             dat.h <- rbind(dat.h0[dat.h0$records==1,],
                            dat.h0[row.names(dat.h0) %in% get(paste("randomPoint_Replication_", replication, sep = "")),])
             obs.h <- dat.h[dat.h[,"folds...run."]==FALSE, "records"]
          }
             
          auc.h <- auc(DATA = data.frame(1:length(obs.h), obs.h, pred.out))[1] %>% as.numeric()
          dat.h3 <- data.frame(obs.h, pred.out)
          boyce.h <- ecospat.boyce(fit = pred.out, obs = dat.h3[dat.h3[,1]==1, 2], PEplot = F)$cor %>% as.numeric()
          
          if(length(unique(obs.h)) > 1 && !all(is.na(pred.out))) {
            thresh_seq <- seq(min(pred.out, na.rm=T), max(pred.out, na.rm=T), length.out = 100)
            tss_vec <- numeric(100)
            
            for(t_idx in 1:100) {
              t_val <- thresh_seq[t_idx]
              tp <- sum(obs.h == 1 & pred.out >= t_val, na.rm=T)
              tn <- sum(obs.h == 0 & pred.out < t_val, na.rm=T)
              fp <- sum(obs.h == 0 & pred.out >= t_val, na.rm=T)
              fn <- sum(obs.h == 1 & pred.out < t_val, na.rm=T)
              
              sens_t <- tp / (tp + fn)
              spec_t <- tn / (tn + fp)
              tss_vec[t_idx] <- (sens_t + spec_t) - 1
            }
            
            opt_t_idx <- which.max(tss_vec)
            tss.h <- tss_vec[opt_t_idx]
            
            t_val <- thresh_seq[opt_t_idx]
            tp <- sum(obs.h == 1 & pred.out >= t_val, na.rm=T)
            tn <- sum(obs.h == 0 & pred.out < t_val, na.rm=T)
            fp <- sum(obs.h == 0 & pred.out >= t_val, na.rm=T)
            fn <- sum(obs.h == 1 & pred.out < t_val, na.rm=T)
            
            sens.h <- tp / (tp + fn)
            spec.h <- tn / (tn + fp)
            
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
      
      auc <- c(auc, auc.h)
      boyce <- c(boyce, boyce.h)
      tss <- c(tss, tss.h)
      kappa <- c(kappa, kappa.h)
      sens <- c(sens, sens.h)
      spec <- c(spec, spec.h)
      
      if(nrep != 1) { message(paste("Results of mean predictions")) }
      message(paste("auc = ", round(auc.h, 2)))
      message(paste("boyce = ", round(boyce.h, 2)))
      message(paste("tss = ", round(tss.h, 2)))
      message(paste("kappa = ", round(kappa.h, 2)))
      message(paste("sens = ", round(sens.h, 2)))
      message(paste("spec = ", round(spec.h, 2)))
    }   
    
    if(h == 1) {
      dat.out <- data.frame(auc, boyce, tss, kappa, sens, spec)
      dat.out$ntree <- HyptoTest[h,"ntree"]
      dat.out$k <- HyptoTest[h,"k"]
      dat.out$power <- HyptoTest[h,"power"]
      dat.out$base <- HyptoTest[h,"base"]
      dat.out$species <- name
      dat.out$Fold <- 1:ncol(folds)
      dat.out.f <- dat.out  
    } else {
      dat.out <- data.frame(auc, boyce, tss, kappa, sens, spec)
      dat.out$ntree <- HyptoTest[h,"ntree"]
      dat.out$k <- HyptoTest[h,"k"]
      dat.out$power <- HyptoTest[h,"power"]
      dat.out$base <- HyptoTest[h,"base"]
      dat.out$species <- name
      dat.out$Fold <- 1:ncol(folds)
      dat.out.f <- rbind(dat.out.f, dat.out)  
    }
  }
  
  dat.out.f2 <- dat.out.f
  dat.out.f2$hp <- paste(dat.out.f2$ntree, dat.out.f2$k, dat.out.f2$power, dat.out.f2$base)
  
  BOYCE <- tapply(dat.out.f2$boyce, dat.out.f2$hp, mean, na.rm = T)
  BOYCE.SD <- tapply(dat.out.f2$boyce, dat.out.f2$hp, sd, na.rm = T)
  AUC <- tapply(dat.out.f2$auc, dat.out.f2$hp, mean, na.rm = T)
  AUC.SD <- tapply(dat.out.f2$auc, dat.out.f2$hp, sd, na.rm = T)
  
  TSS <- tapply(dat.out.f2$tss, dat.out.f2$hp, mean, na.rm = T)
  TSS.SD <- tapply(dat.out.f2$tss, dat.out.f2$hp, sd, na.rm = T)
  KAPPA <- tapply(dat.out.f2$kappa, dat.out.f2$hp, mean, na.rm = T)
  KAPPA.SD <- tapply(dat.out.f2$kappa, dat.out.f2$hp, sd, na.rm = T)
  SENS <- tapply(dat.out.f2$sens, dat.out.f2$hp, mean, na.rm = T)
  SPEC <- tapply(dat.out.f2$spec, dat.out.f2$hp, mean, na.rm = T)
  
  out <- data.frame(BOYCE, BOYCE.SD, AUC, AUC.SD, TSS, TSS.SD, KAPPA, KAPPA.SD, SENS, SPEC)
  out$hp <- rownames(out)
  
  out2 <- merge(x = out, y = dat.out.f2[, c(7:10, 13)], by = "hp", all.y = F)
  out3 <- distinct(out2, .keep_all = TRUE)
  
  write.csv(out3[, 2:15], file = paste("hyper_params_eval_pyro/hyperEval", name, "_pyro.csv", sep = ""), row.names = F)
  write.csv(out3[, 2:15], file = paste("hyper_params_eval_pyro/hyperEval", name, "_pyro.csv", sep = ""), row.names = F)
}