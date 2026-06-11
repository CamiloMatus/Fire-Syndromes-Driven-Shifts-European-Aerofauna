# load required libraries
require(sampbias)
require(terra)
require(stringr)

# load bat data and filtered raster
bat_data_raw <- read.csv("bat_records_europe.csv")
filtered_raster <- rast("raster_presence_bats_at_least_25_records.tif")

# extract valid species from raster
valid_species <- str_replace(names(filtered_raster), "_", " ")

# filter original data to match valid species
bat_data_filtered <- bat_data_raw[bat_data_raw$species %in% valid_species, ]

# format columns for sampbias requirements
bat_data_final <- bat_data_filtered[, c("species", "decimalLongitude", "decimalLatitude")]

# calculate and project sampling bias
bias_calc <- calculate_bias(x = bat_data_final, res = 0.2)
bias_proj <- project_bias(bias_calc)

# save raster with sampling rate
# extracting the 4th layer which contains the total bias weight
raster::writeRaster(bias_proj[[4]], file = "samp_bias_bats.tif", overwrite = TRUE)