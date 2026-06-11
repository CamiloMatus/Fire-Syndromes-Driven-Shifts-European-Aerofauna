# load required libraries
library(sf)
library(dplyr)
library(terra)
library(lubridate)
library(stringr)
library(circular)

# helper functions
sd_robust <- function(x, na.rm = TRUE) { 
    if (na.rm) x <- x[!is.na(x)]
    if (length(x) < 2) return(0)
    return(sd(x, na.rm = na.rm)) 
}

# initial setup and output directories
output_folder <- "final_results"
dir.create(output_folder, showWarnings = FALSE)
terraOptions(tempdir = file.path(output_folder, "terra_temp"))

main_folder <- "./cicatrices_2024/"
# reference raster comes from baselinePredictors process
ref_raster_path <- "../baselinePredictors/biosTopos_europe_5min.tif" 

out_file_pyromes <- file.path(output_folder, "fire_regimes.tif")
out_file_sdm <- file.path(output_folder, "fire_regimes_sdm_imputed.tif")
out_shp_filtered <- file.path(output_folder, "fires_unique_filtered.gpkg")

europe_projection <- "EPSG:3035"

# load and process fire data
shapefile_list <- list.files(main_folder, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
if (length(shapefile_list) == 0) stop("no .shp files found.")

all_fires_sf <- lapply(shapefile_list, function(shp_path) {
    base_name <- basename(shp_path)
    start_year <- NA
    end_year <- NA
    
    match1 <- str_match(base_name, "(\\d{4})_to_(\\d{4})")
    match2 <- str_match(base_name, "to(\\d{4})")
    
    if (!is.na(match1[1,1])) { 
        start_year <- as.numeric(match1[1,2])
        end_year <- as.numeric(match1[1,3]) 
    } else if (!is.na(match2[1,1])) { 
        start_year <- 2001
        end_year <- as.numeric(match2[1,2]) 
    } else { 
        return(NULL) 
    }
    
    obs_period <- end_year - start_year + 1
    shp_data <- tryCatch(st_read(shp_path, quiet = TRUE), error = function(e) NULL)
    
    if (is.null(shp_data) || nrow(shp_data) == 0) return(NULL)
    
    cols_to_convert <- c("ig_day", "ig_year", "event_dur", "tot_pix", "tot_ar_km2", "fsr_px_dy", "fsr_km2_dy", "tot_perim")
    shp_data <- shp_data |> 
        mutate(
            across(any_of(cols_to_convert), as.numeric), 
            ig_date = as_date(ig_date), 
            ig_month = month(ig_date), 
            ig_jday = yday(ig_date), 
            obs_period = obs_period
        )
    return(shp_data)
}) |> bind_rows()

# deduplicate fires
all_fires_deduplicated_sf <- all_fires_sf |> 
    mutate(x_key = round(ig_utm_x, -2), y_key = round(ig_utm_y, -2)) |> 
    group_by(ig_date, x_key, y_key) |> 
    slice_max(order_by = obs_period, n = 1, with_ties = FALSE) |> 
    ungroup() |> 
    mutate(unique_fire_id = 1:n())

# assign country names by spatial location
countries_sf <- read_sf("../baselinePredictors/europaPaises.gpkg")
countries_projected_sf <- st_transform(countries_sf, st_crs(all_fires_deduplicated_sf))
fires_centroids <- st_centroid(all_fires_deduplicated_sf)
fires_with_country <- st_join(fires_centroids, countries_projected_sf, join = st_intersects)

all_fires_deduplicated_sf$country <- fires_with_country$country_name
st_write(all_fires_deduplicated_sf, out_shp_filtered, append = FALSE)

# prepare final raster and data
r <- rast(ref_raster_path)
r_p <- project(r, europe_projection)[[1]]

all_fires_projected_sf <- st_transform(all_fires_deduplicated_sf, crs = europe_projection)
all_fires_projected_sf <- all_fires_projected_sf |> 
    mutate(shape_pa = tot_perim / (tot_ar_km2 + 1e-6))

# calculate fire regimes
cell_id_r <- r_p
values(cell_id_r) <- 1:ncell(r_p)
names(cell_id_r) <- "cell_id"

fire_cell_links <- terra::extract(cell_id_r, all_fires_projected_sf, list = FALSE, ID = TRUE)
fire_attributes <- all_fires_projected_sf |> st_drop_geometry() |> mutate(ID = 1:n())
fire_data_for_summary <- as_tibble(fire_cell_links) |> 
    left_join(fire_attributes, by = "ID") |> 
    filter(!is.na(cell_id))

cell_summary_event <- fire_data_for_summary |> 
    group_by(cell_id) |> 
    summarise(
        fire_years = n_distinct(ig_year, na.rm = TRUE), 
        obs_period = max(obs_period, na.rm = TRUE),
        mean_size = mean(tot_ar_km2, na.rm = TRUE), 
        sd_size = sd_robust(tot_ar_km2, na.rm = TRUE),
        max_size = max(tot_ar_km2, na.rm = TRUE), 
        sum_area = sum(tot_ar_km2, na.rm = TRUE),
        mean_dur = mean(event_dur, na.rm = TRUE), 
        sd_dur = sd_robust(event_dur, na.rm = TRUE),
        mean_exp_rate = mean(fsr_km2_dy, na.rm = TRUE), 
        sd_exp_rate = sd_robust(fsr_km2_dy, na.rm = TRUE),
        mean_shape = mean(shape_pa, na.rm = TRUE), 
        sd_shape = sd_robust(shape_pa, na.rm = TRUE),
        season_conc = {
            circ_days <- circular((ig_jday/365.25)*360, units="degrees")
            as.numeric(rho.circular(circ_days, na.rm=TRUE))
        },
        .groups = 'drop'
    )

annual_summary <- fire_data_for_summary |> 
    group_by(cell_id, ig_year) |> 
    summarise(annual_ba = sum(tot_ar_km2, na.rm = TRUE), .groups = 'drop')

cell_summary_cv_ba <- annual_summary |> 
    group_by(cell_id) |> 
    summarise(
        cv_ba = sd_robust(annual_ba, na.rm = TRUE) / (mean(annual_ba, na.rm = TRUE) + 1e-9), 
        .groups = 'drop'
    )

# join and calculate final metrics
cell_summary_final <- left_join(cell_summary_event, cell_summary_cv_ba, by = "cell_id") |> 
    mutate(
        burn_prob = fire_years / obs_period, 
        fri = 1 / (burn_prob + 1e-9),
        maab = sum_area / obs_period,
        cv_size = sd_size / (mean_size + 1e-9),
        cv_dur = sd_dur / (mean_dur + 1e-9),
        cv_exp_rate = sd_exp_rate / (mean_exp_rate + 1e-9),
        cv_shape = sd_shape / (mean_shape + 1e-9)
    ) |> 
    filter(!is.na(cell_id))

# create and save raster for pyromes
template_raster <- r_p
values(template_raster) <- NA
final_raster_list <- list()

layer_names <- c("burn_prob", "fri", "maab", "cv_ba", 
                 "mean_size", "cv_size", "max_size", 
                 "mean_dur", "cv_dur",
                 "mean_exp_rate", "cv_exp_rate",
                 "mean_shape", "cv_shape",
                 "season_conc")

for (layer_name in layer_names) {
    r_out <- template_raster
    if (layer_name %in% names(cell_summary_final)) { 
        r_out[cell_summary_final$cell_id] <- cell_summary_final[[layer_name]] 
    }
    names(r_out) <- layer_name
    final_raster_list[[layer_name]] <- r_out
}

pyromes_raster <- rast(final_raster_list)

# impute, crop and save raster for sdm
sdm_raster <- pyromes_raster
vars_to_zero <- c("burn_prob", "maab", "cv_ba", "mean_size", "cv_size", "max_size", "mean_dur", "cv_dur", "mean_exp_rate", "cv_exp_rate", "mean_shape", "cv_shape", "season_conc")

for (var in vars_to_zero) { 
    sdm_raster[[var]] <- ifel(is.na(sdm_raster[[var]]), 0, sdm_raster[[var]]) 
}
sdm_raster[["fri"]] <- ifel(is.na(sdm_raster[["fri"]]), 999, sdm_raster[["fri"]])

study_area <- read_sf("../baselinePredictors/europaPaises.gpkg")
projected_study_area <- st_transform(study_area, crs(sdm_raster))

final_pyromes_raster <- crop(pyromes_raster, projected_study_area, mask = TRUE)
final_sdm_raster <- crop(sdm_raster, projected_study_area, mask = TRUE)

writeRaster(final_pyromes_raster, out_file_pyromes, overwrite = TRUE, datatype = "FLT4S")
writeRaster(final_sdm_raster, out_file_sdm, overwrite = TRUE, datatype = "FLT4S")

terra::tmpFiles(remove = TRUE)