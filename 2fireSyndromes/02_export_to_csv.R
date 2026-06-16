# load required library
library(terra)

# set up file paths
output_folder <- "final_results"
input_raster_path <- file.path(output_folder, "fire_regimes.tif")
output_csv_path <- file.path(output_folder, "fire_vars_data.csv")

# load the raster file
if (!file.exists(input_raster_path)) stop("[error] raster file not found.")
regimes_raster <- rast(input_raster_path)

# convert raster to data.frame including xy coordinates
# na.rm = true removes empty pixels
regimes_df <- as.data.frame(regimes_raster, xy = TRUE, na.rm = TRUE)

# save data as csv without row names
write.csv(regimes_df, output_csv_path, row.names = FALSE)