# load required libraries
library(sampbias)
library(terra)
library(data.table)

# load filtered bird records
bird_data <- fread("filtered_bird_records_atLeast250px.csv")

# taxonomy correction map for species with multiple names
correction_map <- list(
    "Corvus corone" = c("Corvus corone", "Corvus cornix"),
    "Passer domesticus" = c("Passer domesticus", "Passer hispaniolensis", "Passer italiae"),
    "Phylloscopus collybita" = c("Phylloscopus collybita", "Phylloscopus ibericus"),
    "Picus viridis" = c("Picus viridis", "Picus sharpei"),
    "Curruca hortensis" = c("Curruca hortensis", "Curruca crassirostris"),
    "Iduna pallida" = c("Iduna opaca", "Iduna pallida")
)

# create a lookup table for replacements
updates_list <- list()
for (correct_name in names(correction_map)) {
    old_names <- correction_map[[correct_name]]
    old_names_to_replace <- setdiff(old_names, correct_name)
    
    if (length(old_names_to_replace) > 0) {
        updates_list[[correct_name]] <- data.table(
            old_name = old_names_to_replace,
            new_name = correct_name
        )
    }
}

corrections_dt <- rbindlist(updates_list)

# update names directly in the data.table using an update join
bird_data[corrections_dt, on = .(SCIENTIFIC_NAME = old_name), SCIENTIFIC_NAME := i.new_name]

# format columns for sampbias package requirements
names(bird_data)[names(bird_data) == "LONGITUDE"] <- "decimalLongitude"
names(bird_data)[names(bird_data) == "LATITUDE"] <- "decimalLatitude"
names(bird_data)[names(bird_data) == "SCIENTIFIC_NAME"] <- "species"

# reorder columns to species, decimalLongitude, decimalLatitude
bird_data <- bird_data[, .(species, decimalLongitude, decimalLatitude)]

# calculate and project sampling bias
bias_calc <- calculate_bias(x = bird_data, res = 0.2)
bias_proj <- project_bias(bias_calc)

# save raster with sampling rate
# extracting the 4th layer which contains the total bias weight
raster::writeRaster(bias_proj[[4]], file = "samp_bias.tif", overwrite = TRUE)