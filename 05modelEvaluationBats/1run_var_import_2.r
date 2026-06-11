library(fastDummies)
require(embarcadero)
require(PresenceAbsence)
require(sf)
require(dplyr)
require(ecospat)
require(foreign)



solo_ceros <- function(x) {
  all(x == 0, na.rm = TRUE)
}


#reading data
db <- read.csv("../2make_modeling_data/datForModeling.csv")%>%na.omit()

species <- names(db)[3:33]

# envVars<-db[,names(db)[c(297:323,335,337)]]

envVars<- db[,names(db)[c(34:64)]]


envVars$Pyrome_Consenso <- as.factor(envVars$Pyrome_Consenso)

names(envVars)

envVars_dummy <- dummy_cols(envVars, 
                            select_columns = "Pyrome_Consenso", 
                            remove_selected_columns = TRUE)

db<-cbind(db,envVars_dummy[,c("Pyrome_Consenso_1","Pyrome_Consenso_2",
                             "Pyrome_Consenso_3","Pyrome_Consenso_4","Pyrome_Consenso_5",
                             "Pyrome_Consenso_6","Pyrome_Consenso_7","Pyrome_Consenso_8",
                             "Pyrome_Consenso_9","Pyrome_Consenso_10","Pyrome_Consenso_11")])

#predicctions by species
for(i in species){
  


namesEV <- c("bio1","bio2", "bio4", "bio8", "bio12", "bio15", "elevacion",
                             "prop_arboles", "prop_matorrales", "prop_pastizales",
                             "prop_cultivos", "prop_construido", "prop_suelo_desnudo",
                             "prop_agua", "prop_humedales","Pyrome_Consenso_1","Pyrome_Consenso_2",
                             "Pyrome_Consenso_3","Pyrome_Consenso_4","Pyrome_Consenso_5",
                             "Pyrome_Consenso_6","Pyrome_Consenso_7","Pyrome_Consenso_8",
                             "Pyrome_Consenso_9","Pyrome_Consenso_10","Pyrome_Consenso_11", 
                              "sampBias")
                              
  message(paste("working on ",i,"   ",which(species==i ),"of",length(species)))
  
  
  #condition for species with few records
  if(sum(db[,i])<150){nrep <- 10}else{nrep <- 1}
  
  for (replication in 1:nrep) {
    if(nrep!=1){
      message(paste("Fitting for replication",replication,"of",nrep))}else{NULL}
    
    #random pseudo-absences used in model evaluation
    rp <- read.csv(paste("../3sdmEval_pyrosStudy/randomPseudo-Absences/randomPseudoAbs_Replication_",replication,"_",i,".csv",sep = ""))
    
    #data with pseudo-absences used in model evaluation
    dat.1 <- db[row.names(db)%in%rp$x,]
    
    #data with only presences records
    dat.2 <- db[db[,i]==1,]
    
    #data for model fitting
    dat.h <- rbind(dat.1,dat.2)

    dat.h<-dat.h[,names(dat.h)[!names(dat.h)%in%   names( sapply(dat.h,solo_ceros)[sapply(dat.h,solo_ceros)])]]
 
    namesEV<-namesEV[namesEV%in%names(sapply(dat.h,solo_ceros))]
    
    #evaluating variable importance
    message("evaluating variable importance")
    vi.h<-varimp.diag(dat.h[,namesEV],dat.h[,i],iter=2)
    
    write.csv(vi.h$data,file=paste0("./varImp/","varImp_Pyromes_",i,"_rep_",replication,".csv"),row.names=F)
    }
}