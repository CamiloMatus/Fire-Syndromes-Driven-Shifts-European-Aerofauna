# environment setup
library(terra)
library(sf)
library(stringr)

# quality parameter: minimum 5% land to consider cell valid
land_threshold <- 0.05 

if (!dir.exists("temp.terra")) dir.create("temp.terra")
terraOptions(tempdir = "temp.terra", memfrac = 0.5, progress = 0)

template_file <- "biosTopos_europe_5min.tif"
if(!file.exists(template_file)) stop("[ERROR] template file not found.")
template_raster <- rast(template_file)
raster_crs <- crs(template_raster)

# land mask preparation with 50m negative buffer
if(!file.exists("europaPaises.gpkg")) stop("[ERROR] europaPaises.gpkg not found.")
countries_raw_sf <- read_sf("europaPaises.gpkg")

# apply 50m coastal retraction via terra for speed
# 1. unify geometries to avoid inland border gaps
# 2. project to european metric standard epsg:3035
# 3. apply 50m negative buffer
countries_vect <- vect(countries_raw_sf)
continental_mask <- countries_vect |> 
    aggregate() |>              
    project("EPSG:3035") |>     
    buffer(width = -50) |>      
    project(raster_crs)          

# tiles preparation
tile_paths <- list.files("./LUC_tilesEuropa/Europa", pattern = "\\.tif$", full.names = TRUE)
tile_paths <- tile_paths[!str_detect(tile_paths, "aux")]
global_vrt <- vrt(tile_paths)

total_tiles <- length(tile_paths)
land_classes <- list(
    prop_trees = 10, prop_shrubs = 20, prop_grasslands = 30, 
    prop_crops = 40, prop_built = 50, prop_bare_soil = 60, 
    prop_water = 80, prop_wetlands = 90, prop_mangroves = 95
)

dir.create("tiles_5min_processed", showWarnings = FALSE)

# process each tile
for (i in 1:total_tiles) {
    base_name <- basename(tile_paths[i])
    temp_out_path <- file.path("tiles_5min_processed", paste0("proc_", base_name))
    
    if (file.exists(temp_out_path)) next
    
    # load buffered zone
    tile_ext <- ext(rast(tile_paths[i]))
    zone_raster <- crop(global_vrt, tile_ext + 0.1)
    
    # crop mask to tile area
    zone_mask <- crop(continental_mask, tile_ext + 0.1)
    
    # skip if tile is purely ocean
    if(nrow(zone_mask) == 0) next
    
    # apply mask (ocean becomes NA)
    zone_raster <- mask(zone_raster, zone_mask)
    
    # calculate land availability for quality filter
    land_10m <- !is.na(zone_raster)
    land_prop_5min <- project(land_10m, template_raster, method = "average")
    
    # calculate proportions
    tile_classes_list <- list()
    for (c_idx in 1:length(land_classes)) {
        c_name <- names(land_classes)[c_idx]
        c_code <- land_classes[[c_idx]]
        
        # average ignores NA representing ocean to calculate true proportion over land
        final_prop <- project(zone_raster == c_code, template_raster, method = "average")
        
        # filter: discard cell if it has very little land
        final_prop[land_prop_5min < land_threshold] <- NA
        
        names(final_prop) <- c_name
        tile_classes_list[[c_name]] <- final_prop
    }
    
    # save processed tile
    tile_result <- crop(rast(tile_classes_list), tile_ext)
    writeRaster(tile_result, temp_out_path, overwrite = TRUE)
    
    rm(zone_raster, tile_classes_list, tile_result, zone_mask, land_10m, land_prop_5min)
    gc()
}

# final mosaic and crop
final_tiles <- list.files("tiles_5min_processed", full.names = TRUE, pattern = ".tif")
final_mosaic <- mosaic(sprc(final_tiles))

# sync crs of original vector for aesthetic cropping
final_countries_vect <- project(countries_vect, final_mosaic)
out_raster <- crop(final_mosaic, final_countries_vect, mask = TRUE)

output_filename <- "land_proportions_europe_5min_final.tif"
writeRaster(out_raster, output_filename, overwrite = TRUE)