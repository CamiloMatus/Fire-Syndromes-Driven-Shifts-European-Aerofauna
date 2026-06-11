# load required libraries
require(tidyverse)
require(ade4)
require(mFD)
require(stringr)
require(foreach)
require(doParallel)

db_ref <- read.csv("../4modelingDataBirds/data_for_modeling.csv") |> na.omit()
species_list <- names(db_ref)[3:296]

# load and filter traits dataset
df_traits <- read.csv("trait_bird_Europe_data.csv")
df_traits <- df_traits[df_traits$Scientific_Name_ID %in% species_list, ]
df_traits$Scientific_Name_ID <- str_replace(df_traits$Scientific_Name_ID, "_", " ")

df_traits$Diet.Category <- as.factor(df_traits$Diet.Category)
df_traits$Phen.MigrantStatus <- as.factor(df_traits$Phen.MigrantStatus)
df_traits$Ecol.Habitat <- as.factor(df_traits$Ecol.Habitat)
df_traits$Ecol.TrophicNiche <- as.factor(df_traits$Ecol.TrophicNiche)
df_traits$Ecol.PrimaryLifestyle  <- as.factor(df_traits$Ecol.PrimaryLifestyle)

trait_type <- lapply(df_traits, class)[2:ncol(df_traits)] |> data.frame() |> t()
trait_name <- names(df_traits)[2:ncol(df_traits)]
fuzzy_name <- substr(trait_name, 1, 4)

# weight of traits
dw <- data.frame((1 / length(unique(fuzzy_name))) / table(fuzzy_name))

df3 <- data.frame(trait_name, trait_type, fuzzy_name)
row.names(df3) <- 1:nrow(df3)
df4 <- merge(df3, dw, by = "fuzzy_name")
df5 <- df4[, 2:4]
names(df5)[3] <- "trait_weight"

df5[df5$trait_type == "numeric", "trait_type"] <- "Q"
df5[df5$trait_type == "factor", "trait_type"] <- "N"
df5[df5$trait_type == "integer", "trait_type"] <- "Q"

# compute functional distance
row.names(df_traits) <- df_traits$Scientific_Name_ID
dist_trait <- funct.dist(
    sp_tr = df_traits[, df5$trait_name],
    tr_cat = df5,
    metric = "gower",
    stop_if_NA = FALSE
)

fspaces_quality <- mFD::quality.fspaces(
    sp_dist = dist_trait,
    maxdim_pcoa = 20,
    deviation_weighting = "absolute",
    fdist_scaling = FALSE,
    fdendro = "average"
)

# save quality metrics
dir.create("./final_tabs", showWarnings = FALSE)
sink("./final_tabs/quality_fs_baseline.txt")
print(fspaces_quality$quality_fspaces)
sink()

# retrieving principal coordinates
pco <- dudi.pco(dist_trait, scann = FALSE)

sink("./final_tabs/summary_pco_baseline.txt")
inertia.dudi(pco)
sink()

names_pco <- row.names(pco$tab)
sp_faxes <- as.matrix(pco$tab)
row.names(sp_faxes) <- row.names(pco$tab)

# making df estimates in each geographic grid
files <- list.files("./predictions")
files <- files[str_detect(files, "predictions_baseline_")]

for(i in 1:length(files)) {
    f_h <- files[i]
    d_h <- read.csv(paste0("./predictions/", f_h))
    c_h <- data.frame(d_h[, 4])
    
    if(i == 1) { out <- c_h } else { out <- cbind(out, c_h) }
}

nf <- str_remove(files, "predictions_baseline_")
nf2 <- str_remove(nf, "\\.csv")
nf3 <- str_replace(nf2, "_", " ")
names(out) <- nf3

db <- out

# register parallel cluster
cl <- makeCluster(8)
registerDoParallel(cl)

db$sp_richn <- NA
db$fdis <- NA
db$fmpd <- NA
db$fnnd <- NA
db$feve <- NA
db$fric <- NA
db$fdiv <- NA
db$fori <- NA
db$fspe <- NA

# processing functional diversity predictions
output_file <- "func_div_predictions_baseline_incremental.csv"
chunk_size <- 500 

write.csv(db[0, ], file = output_file, row.names = FALSE)

row_chunks <- split(1:nrow(db), ceiling(seq_along(1:nrow(db)) / chunk_size))

for (i in 1:length(row_chunks)) {
    
    current_rows <- row_chunks[[i]]
    message("processing batch ", i, " of ", length(row_chunks))
    
    out_chunk <- foreach(j = current_rows, .packages = "dplyr", .combine = "rbind") %dopar% {
        
        a <- rep(0, 294)
        a <- data.frame(a)
        row.names(a) <- row.names(pco$tab)
        names_h <- names(colSums(db[j, 1:294])[colSums(db[j, 1:294]) > 0])
        a[rownames(a) %in% names_h, "a"] <- 1
        
        if (length(names_h) > 7) {
            alpha_fd_indices <- mFD::alpha.fd.multidim(
                sp_faxes_coord = sp_faxes[, 1:7],
                asb_sp_w = t(a),
                ind_vect = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", "fspe"),
                scaling = TRUE, check_input = TRUE, details_returned = TRUE, verbose = FALSE
            )
            
            result_row <- db[j, ]
            result_row$sp_richn <- length(names_h)
            result_row$fdis <- alpha_fd_indices$functional_diversity_indices$fdis
            result_row$fmpd <- alpha_fd_indices$functional_diversity_indices$fmpd
            result_row$fnnd <- alpha_fd_indices$functional_diversity_indices$fnnd
            result_row$feve <- alpha_fd_indices$functional_diversity_indices$feve
            result_row$fric <- alpha_fd_indices$functional_diversity_indices$fric
            result_row$fdiv <- alpha_fd_indices$functional_diversity_indices$fdiv
            result_row$fori <- alpha_fd_indices$functional_diversity_indices$fori
            result_row$fspe <- alpha_fd_indices$functional_diversity_indices$fspe
            
        } else {
            result_row <- db[j, ]
            result_row$sp_richn <- length(names_h)
        }
        
        result_row 
    } 
    
    write.table(out_chunk, file = output_file, sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
}

stopCluster(cl)