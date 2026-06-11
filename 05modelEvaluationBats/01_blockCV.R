# load required libraries
require(blockCV)
require(sf)
require(doParallel)
require(tidyverse)
require(terra)

# set up parallel processing cluster
cl <- makeCluster(5)
registerDoParallel(cl)

# load modeling data and reference raster from previous step
modeling_data <- read.csv("../04modelingDataBats/data_for_modeling.csv")
ref_raster <- rast("../04modelingDataBats/data_for_modeling.tif")

# convert to spatial vector
# assuming first 31 columns are species based on the bats modeling data structure
data_vect <- vect(modeling_data, geom = c("x", "y"), crs = crs(ref_raster))
species_list <- names(data_vect)[1:31]

# convert to sf object for blockcv compatibility
data_sf <- st_as_sf(data_vect)

# create output directory if it does not exist
dir.create("./blocks", showWarnings = FALSE)

# run spatial cross-validation blocking in parallel for each species
foreach(sp_name = species_list, .packages = "blockCV", .combine = "c") %dopar% {
    
    spatial_blocks <- cv_spatial(
        x = data_sf[, sp_name],
        column = sp_name,
        size = 200000,
        k = 5,
        plot = FALSE,
        selection = "systematic"
    )

    blocks_df <- data.frame(spatial_blocks$biomod_table)
    out_filename <- paste0("./blocks/blocksCV_", sp_name, ".csv")
    write.csv(blocks_df, file = out_filename, row.names = FALSE)
    
    return(NULL)
}

# stop cluster
stopCluster(cl)