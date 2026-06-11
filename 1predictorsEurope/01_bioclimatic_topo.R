# load required libraries
require(geodata)
require(terra)
require(spatialEco)
require(elevatr)
require(sf)

# set up the area of interest (europe excluding russia)
europe_country_codes <- c(
    "ALB", "AND", "AUT", "BLR", "BEL", "BIH", "BGR", "HRV", "CYP", "CZE", 
    "DNK", "EST", "FIN", "FRA", "DEU", "GRC", "HUN", "ISL", "IRL", "ITA",
    "LVA", "LIE", "LTU", "LUX", "MLT", "MDA", "MCO", "MNE", "NLD",
    "MKD", "NOR", "POL", "PRT", "ROU", "SMR", "SRB", "SVK", "SVN", 
    "ESP", "SWE", "CHE", "UKR", "GBR", "VAT" 
)

# download and validate polygon geometries
countries_vect <- gadm(country = europe_country_codes, level = 0, path = getwd())
countries_sf <- st_as_sf(countries_vect) |> st_make_valid()

# set resolutions to process in minutes
resolutions <- c(2.5, 5, 10)

# main loop for each resolution
for (res in resolutions) {

    # get and prep bioclimatic variables
    bios <- worldclim_global(var = "bio", res = res, path = getwd())
    bios_cropped <- crop(bios, countries_sf, mask = TRUE)

    # get elevation data (adjusting zoom level based on resolution)
    zoom_level <- ifelse(res == 2.5, 6, ifelse(res == 5, 5, 4))
    high_res_elev_raster <- get_elev_raster(locations = countries_sf, z = zoom_level, src = "aws")
    
    # convert, resample and calculate topo variables
    high_res_elev_terra <- rast(high_res_elev_raster)
    elev_resampled <- resample(high_res_elev_terra, bios_cropped)

    slope_aspect <- terrain(elev_resampled, v = c("slope", "aspect"), unit = "degrees")
    tpi_raster <- tpi(elev_resampled)

    # create, mask and save final raster stack
    unmasked_stack <- c(bios_cropped, elev_resampled, slope_aspect, tpi_raster)
    names(unmasked_stack) <- c(paste0("bio", 1:19), "elevation", "slope", "aspect", "tpi")

    final_mask <- bios_cropped[[1]]
    final_stack <- mask(unmasked_stack, final_mask)

    # save output file
    output_filename <- paste0("biosTopos_europe_", res, "min.tif")
    writeRaster(final_stack, output_filename, overwrite = TRUE)
    print(final_stack)
}

# extra step: clip 5min resolution to custom europe borders geopackage
p_europe <- read_sf("europaPaises.gpkg")
r_h <- rast("biosTopos_europe_5min.tif")

p_europe <- st_transform(p_europe, crs(r_h))
r_out <- crop(r_h, p_europe, mask = TRUE)

writeRaster(r_out, "biosTopos_Europe_final.tif", overwrite = TRUE)