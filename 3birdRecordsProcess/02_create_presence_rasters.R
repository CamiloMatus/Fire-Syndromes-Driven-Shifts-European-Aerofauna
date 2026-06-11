# load required libraries
library(terra)
library(data.table)
library(dplyr)
library(lubridate)

# set temp directory for terra to avoid filling up main drive
terraOptions(tempdir = "temp_terra_dir")

# file paths configuration
template_raster_path <- "../baselinePredictors/biosTopos_europe_5min.tif"
study_area_path <- "../baselinePredictors/europaPaises.gpkg"
ebird_data_path <- "ebd_europe_prefiltered_clean.tsv.gz"
phenology_table_path <- "tabFeno.csv"

out_all_records <- "raster_presence_final_allRecords.tif"
out_filtered_250 <- "raster_presence_final_atLeast250Records.tif"
out_final_corrected <- "records_species_final.tif"

# load spatial base elements
template_raster <- rast(template_raster_path)
study_area_vect <- vect(study_area_path)
study_area_proj <- project(study_area_vect, crs(template_raster))

# load and process phenology table
pheno_tab <- fread(phenology_table_path)

# map spanish months from csv to numeric values
month_map <- c(
    "Enero" = 1, "Febrero" = 2, "Marzo" = 3, "Abril" = 4, 
    "Mayo" = 5, "Junio" = 6, "Julio" = 7, "Agosto" = 8, 
    "Septiembre" = 9, "Octubre" = 10, "Noviembre" = 11, "Diciembre" = 12
)

pheno_proc <- pheno_tab |>
    mutate(
        start_month_num = month_map[`Mes de Inicio Cría`],
        end_month_num = month_map[`Mes de Fin Cría`]
    ) |>
    select(
        scientific_name = `Scientific_Name`,
        classification = `Classification`,
        start_month = start_month_num,
        end_month = end_month_num,
        obs_mode = `Modo de Observación Principal`
    )

species_to_model <- pheno_proc |>
    filter(classification != "Excluir") |>
    pull(scientific_name)

aerial_species <- pheno_proc |>
    filter(obs_mode == "Aéreo") |>
    pull(scientific_name)

# load and filter ebird data (keep only 2015 onwards)
ebird_dt <- fread(
    cmd = paste("gzip -dc", shQuote(ebird_data_path)),
    sep = "\t",
    encoding = "UTF-8"
)

ebird_dt <- ebird_dt[year(as.Date(`OBSERVATION DATE`)) >= 2015]

# define rasterization function
create_presence_raster <- function(species_df, sci_name, species_rules, aerial_list, template_r, study_area_p) {
    
    # filter by sampling quality
    filtered_df <- species_df |>
        filter(`ALL SPECIES REPORTED` == 1,
               is.na(`DURATION MINUTES`) | `DURATION MINUTES` <= 300,
               is.na(`EFFORT DISTANCE KM`) | `EFFORT DISTANCE KM` <= 5)
    
    # filter by flight behavior (except aerial species)
    if (!(sci_name %in% aerial_list)) {
        filtered_df <- filtered_df |>
            filter(is.na(`BEHAVIOR CODE`) | !grepl("F", `BEHAVIOR CODE`))
    }
    
    # filter by breeding season if applicable
    if (species_rules$classification == "Filtrar_Cria" && !is.na(species_rules$start_month)) {
        filter_months <- seq(species_rules$start_month, species_rules$end_month)
        filtered_df <- filtered_df |>
            mutate(obs_month = month(as.Date(`OBSERVATION DATE`))) |>
            filter(obs_month %in% filter_months)
    }
    
    # prepare for rasterization
    final_df <- filtered_df |>
        filter(!is.na(LATITUDE) & !is.na(LONGITUDE)) |>
        mutate(presence = 1)
    
    if (nrow(final_df) == 0) {
        empty_r <- rast(template_r)
        values(empty_r) <- 0
        return(crop(mask(empty_r, study_area_p), study_area_p))
    }
    
    # rasterize points
    points_vect <- vect(final_df, geom = c("LONGITUDE", "LATITUDE"), crs = "EPSG:4326")
    points_proj <- project(points_vect, crs(template_r))
    
    presence_r <- rasterize(x = points_proj, y = template_r, field = "presence", fun = "max")
    presence_r[is.na(presence_r)] <- 0
    
    return(crop(mask(presence_r, study_area_p), study_area_p))
}

# main loop to process each species
species_rasters_list <- list()

for (current_species in species_to_model) {
    
    current_rules <- pheno_proc[scientific_name == current_species]
    current_data <- ebird_dt[`SCIENTIFIC NAME` == current_species]
    
    if (nrow(current_data) == 0) next
    
    result_raster <- create_presence_raster(
        species_df = current_data,
        sci_name = current_species,
        species_rules = current_rules,
        aerial_list = aerial_species,
        template_r = template_raster,
        study_area_p = study_area_proj
    )
    
    if (nlyr(result_raster) == 1) {
        layer_name <- gsub(" ", "_", current_species)
        names(result_raster) <- layer_name
        species_rasters_list[[layer_name]] <- result_raster
    }
}

# stack all valid rasters
if (length(species_rasters_list) == 0) stop("[error] no valid rasters generated.")
multi_raster <- rast(species_rasters_list)
writeRaster(multi_raster, out_all_records, overwrite = TRUE)

# filter layers with more than 250 presence pixels
pixels_count <- global(multi_raster, fun = "sum", na.rm = TRUE)
pixels_count$species <- rownames(pixels_count)

species_to_keep <- subset(pixels_count, sum >= 250)$species
if(length(species_to_keep) == 0) stop("[error] no species passed the 250 records threshold.")

filtered_raster <- multi_raster[[species_to_keep]]
writeRaster(filtered_raster, out_filtered_250, overwrite = TRUE)

# taxonomy correction map
correction_map <- list(
    "Corvus_corone" = c("Corvus_corone", "Corvus_cornix"),
    "Passer_domesticus" = c("Passer_domesticus", "Passer_hispaniolensis", "Passer_italiae"),
    "Phylloscopus_collybita" = c("Phylloscopus_collybita", "Phylloscopus_ibericus"),
    "Picus_viridis" = c("Picus_viridis", "Picus_sharpei"),
    "Curruca_hortensis" = c("Curruca_hortensis", "Curruca_crassirostris"),
    "Iduna_pallida" = c("Iduna_opaca", "Iduna_pallida")
)

# merge duplicate layers
merged_layers <- list()
all_merged_names <- unlist(correction_map)
available_layers <- names(filtered_raster)

for (correct_name in names(correction_map)) {
    layers_to_merge <- correction_map[[correct_name]]
    valid_layers_to_merge <- intersect(layers_to_merge, available_layers)
    
    if (length(valid_layers_to_merge) > 0) {
        subset_raster <- filtered_raster[[valid_layers_to_merge]]
        merged_layer <- max(subset_raster, na.rm = TRUE)
        names(merged_layer) <- correct_name
        merged_layers[[correct_name]] <- merged_layer
    }
}

# combine untouched layers with merged layers
merged_raster_stack <- rast(merged_layers)
untouched_names <- setdiff(available_layers, all_merged_names)
untouched_raster_stack <- filtered_raster[[untouched_names]]

final_corrected_raster <- c(merged_raster_stack, untouched_raster_stack)

# save final output
writeRaster(final_corrected_raster, out_final_corrected, overwrite = TRUE)

# cleanup temp files
terra::tmpFiles(remove = TRUE)