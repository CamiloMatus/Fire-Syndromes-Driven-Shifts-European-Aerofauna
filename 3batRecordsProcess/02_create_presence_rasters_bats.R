# load required libraries
library(terra)
library(data.table)
library(dplyr)

# file paths configuration
template_raster_path <- "../baselinePredictors/biosTopos_europe_5min.tif"
study_area_path <- "../baselinePredictors/europaPaises.gpkg"
bat_data_path <- "bat_records_europe.csv"

out_all_records <- "raster_presence_bats_all_records.tif"
out_filtered_records <- "raster_presence_bats_at_least_25_records.tif"

# threshold of minimum presence pixels to keep a species
pixel_threshold <- 25 

# load spatial base elements
template_raster <- rast(template_raster_path)
study_area_vect <- vect(study_area_path)
study_area_proj <- project(study_area_vect, crs(template_raster))

# load bat data
bat_dt_raw <- fread(bat_data_path)

# clean and standardize columns
bat_dt <- bat_dt_raw |>
    select(
        scientific_name = species,
        latitude = decimalLatitude,
        longitude = decimalLongitude
    ) |>
    filter(
        !is.na(latitude) & !is.na(longitude), 
        !is.na(scientific_name),            
        scientific_name != ""                  
    )
  
species_to_model <- unique(bat_dt$scientific_name)

# define simple rasterization function
create_simple_presence_raster <- function(species_df, template_r, study_area_p) {
    
    final_df <- species_df |> mutate(presence = 1)
    
    if (nrow(final_df) == 0) {
        empty_r <- rast(template_r)
        values(empty_r) <- 0
        return(crop(mask(empty_r, study_area_p), study_area_p))
    }
    
    points_vect <- vect(final_df, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    points_proj <- project(points_vect, crs(template_r))
    
    presence_r <- rasterize(x = points_proj, y = template_r, field = "presence", fun = "max")
    presence_r[is.na(presence_r)] <- 0
    
    return(crop(mask(presence_r, study_area_p), study_area_p))
}

# process each species
species_rasters_list <- list()

for (current_species in species_to_model) {
    
    current_data <- bat_dt[scientific_name == current_species]
    if (nrow(current_data) == 0) next
    
    result_raster <- create_simple_presence_raster(
        species_df = current_data,
        template_r = template_raster,
        study_area_p = study_area_proj
    )
    
    layer_name <- gsub(" ", "_", current_species)
    names(result_raster) <- layer_name
    species_rasters_list[[layer_name]] <- result_raster
}

# stack and save all records
if(length(species_rasters_list) > 0) {
    multi_raster <- rast(species_rasters_list)
    writeRaster(multi_raster, out_all_records, overwrite = TRUE)
} else {
    stop("[error] no valid rasters generated.")
}

# filter by pixel count
pixels_count <- global(multi_raster, fun = "sum", na.rm = TRUE)
pixels_count$species <- rownames(pixels_count)

species_to_keep <- subset(pixels_count, sum >= pixel_threshold)$species

if(length(species_to_keep) > 0) {
    filtered_raster <- multi_raster[[species_to_keep]]
    writeRaster(filtered_raster, out_filtered_records, overwrite = TRUE)
}