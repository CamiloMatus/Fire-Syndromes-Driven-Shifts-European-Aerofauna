# load required libraries
library(terra)
library(data.table)
library(dplyr)
library(lubridate)

# file paths configuration
template_raster_path <- "../baselinePredictors/biosTopos_europe_5min.tif"
study_area_path <- "../baselinePredictors/europaPaises.gpkg"
ebird_data_path <- "ebd_europe_prefiltered_clean.tsv.gz"
phenology_table_path <- "tabFeno.csv"
output_csv_path <- "filtered_bird_records_atLeast250px.csv"

template_raster <- rast(template_raster_path)

# load study area and ensure project matches coordinates for ebird data compatibility
study_area_vect_original <- vect(study_area_path)
study_area_vect <- project(study_area_vect_original, "EPSG:4326")

# load and process reference tables
pheno_tab <- fread(phenology_table_path)

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

ebird_dt <- fread(ebird_data_path)

# filter observations for recent years only
ebird_dt <- ebird_dt[year(as.Date(`OBSERVATION DATE`)) >= 2015]

# define filtering function for species records
filter_species_records <- function(species_df, sci_name, species_rules, aerial_list) {
    filtered_df <- species_df |>
        filter(`ALL SPECIES REPORTED` == 1,
               is.na(`DURATION MINUTES`) | `DURATION MINUTES` <= 300,
               is.na(`EFFORT DISTANCE KM`) | `EFFORT DISTANCE KM` <= 5)
    
    if (!(sci_name %in= aerial_list)) {
        filtered_df <- filtered_df |>
            filter(is.na(`BEHAVIOR CODE`) | !grepl("F", `BEHAVIOR CODE`))
    }
    
    if (species_rules$classification == "Filtrar_Cria" && !is.na(species_rules$start_month)) {
        filter_months <- seq(species_rules$start_month, species_rules$end_month)
        filtered_df <- filtered_df |>
            mutate(obs_month = month(as.Date(`OBSERVATION DATE`))) |>
            filter(obs_month %in% filter_months)
    }
    
    final_df <- filtered_df |>
        filter(!is.na(LATITUDE) & !is.na(LONGITUDE))
    
    return(final_df)
}

# process each species records and spatial matching
filtered_data_list <- list() 
pixel_counts <- c()           

for (current_species in species_to_model) {
    
    current_rules <- pheno_proc[scientific_name == current_species]
    current_data <- ebird_dt[`SCIENTIFIC NAME` == current_species]
    
    if (nrow(current_data) == 0) next
    
    result_df <- filter_species_records(
        species_df = current_data,
        sci_name = current_species,
        species_rules = current_rules,
        aerial_list = aerial_species
    )
    
    if (nrow(result_df) > 0) {
        points_vect <- vect(result_df, geom = c("LONGITUDE", "LATITUDE"), crs = "EPSG:4326")
        points_in_area <- points_vect[study_area_vect, ]
        
        if (length(points_in_area) > 0) {
            final_result_df <- as.data.frame(points_in_area, geom = "XY") |>
                rename(LONGITUDE = x, LATITUDE = y)
        } else {
            final_result_df <- data.frame()
        }
    } else {
        final_result_df <- data.frame()
    }
    
    if (nrow(final_result_df) > 0) {
        points_vect_count <- vect(final_result_df, geom = c("LONGITUDE", "LATITUDE"), crs = "EPSG:4326")
        points_proj_count <- project(points_vect_count, crs(template_raster))
        
        cells_idx <- terra::cells(template_raster, points_proj_count)
        unique_pixels_count <- length(unique(cells_idx[, "cell"]))
        
        pixel_counts[current_species] <- unique_pixels_count
        filtered_data_list[[current_species]] <- final_result_df
    } else {
        pixel_counts[current_species] <- 0
    }
}

# filter species by pixel count threshold and export dataset
species_to_keep <- names(pixel_counts[pixel_counts >= 250])

if(length(species_to_keep) > 0) {
    final_dfs_list <- filtered_data_list[species_to_keep]
    complete_final_data <- rbindlist(final_dfs_list, fill = TRUE)
    
    final_csv_data <- complete_final_data |>
        select(
            SCIENTIFIC_NAME = `SCIENTIFIC NAME`,
            LATITUDE,
            LONGITUDE
        )

    fwrite(final_csv_data, output_csv_path)
}