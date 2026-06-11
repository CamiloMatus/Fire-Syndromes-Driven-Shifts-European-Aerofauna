# load required libraries
require(terra)

# load raster layers from previous steps
r0 <- rast("../3batRecordsProcess/raster_presence_bats_at_least_25_records.tif")
r1 <- rast("../1predictorsEurope/biosTopos_europe_5min.tif")
r2 <- rast("../1predictorsEurope/land_proportions_europe_5min_final.tif")
r3 <- rast("../2fireSyndromes/db_pyromes_final.tif") 
r4 <- rast("../3batRecordsProcess/samp_bias_bats.tif")

# define template raster
template_raster <- r1

# ---------------------------------------------------------
# adjust sampling bias data (r4)
# ---------------------------------------------------------

if (!same.crs(r4, template_raster)) {
    r4 <- project(r4, template_raster)
}

# resample continuous bias data using bilinear method
r4_resampled <- resample(r4, template_raster, method = "bilinear")
r4_resampled[!is.finite(r4_resampled)] <- NA

# ensure template matches its own grid correctly
if (!same.crs(r1, template_raster) || any(res(r1) != res(template_raster)) || !ext(r1) == ext(template_raster)) {
    r1 <- resample(project(r1, template_raster), template_raster, method = "bilinear")
}

# create mask from template and fill bias nas with 0 where template has data
mask_r1 <- !is.na(r1[[1]])
r4_covered <- ifel(is.na(r4_resampled) & mask_r1, 0, r4_resampled)

# ---------------------------------------------------------
# adjust pyromes data (r3)
# ---------------------------------------------------------

# separate continuous and categorical layers
continuous_names <- c("Prob_Pyrome_1", "Prob_Pyrome_2", "Prob_Pyrome_3", "Prob_Pyrome_4", 
                      "Prob_Pyrome_5", "Prob_Pyrome_6", "Prob_Pyrome_7", "Prob_Pyrome_8", 
                      "Prob_Pyrome_9", "Prob_Pyrome_10", "Prob_Pyrome_11")
categorical_names <- c("Pyrome_Consenso", "macro_pyro")

r3_continuous <- subset(r3, continuous_names)
r3_categorical <- subset(r3, categorical_names)

# align continuous part using bilinear and clamp to keep probabilities between 0 and 1
r3_cont_temp <- project(r3_continuous, template_raster, method = "bilinear")
r3_cont_aligned <- clamp(r3_cont_temp, lower = 0, upper = 1)

# align categorical part using nearest neighbor
r3_cat_aligned <- project(r3_categorical, template_raster, method = "near")

# merge aligned pyromes layers
r3_aligned_full <- c(r3_cont_aligned, r3_cat_aligned)

# ---------------------------------------------------------
# create final stack
# ---------------------------------------------------------

r0_aligned <- extend(r0, template_raster)
r2_extended <- extend(r2, r0_aligned)

final_stack <- c(r0_aligned, r1, r2_extended, r3_aligned_full, r4_covered)
writeRaster(final_stack, "data_for_modeling.tif", overwrite = TRUE)

# ---------------------------------------------------------
# convert to dataframe and clean up
# ---------------------------------------------------------

full_df <- as.data.frame(final_stack, xy = TRUE)

# filter valid pixels based on climate and land cover presence
df_valid_climate <- full_df[!is.na(full_df$bio1), ]
df_valid_pixels <- df_valid_climate[!is.na(df_valid_climate$prop_bare_soil), ]

# fill NA pyrome probabilities and consensus with 0
prob_cols <- continuous_names
for (col in prob_cols) {
    df_valid_pixels[is.na(df_valid_pixels[[col]]), col] <- 0
}
df_valid_pixels[is.na(df_valid_pixels$Pyrome_Consenso), "Pyrome_Consenso"] <- 0

# clean macro pyrome categories
df_valid_pixels$macro_pyro <- as.character(df_valid_pixels$macro_pyro)
df_valid_pixels[is.na(df_valid_pixels$macro_pyro), "macro_pyro"] <- "M0"

# rename sampling bias column
names(df_valid_pixels)[ncol(df_valid_pixels)] <- "samp_bias"

# extract species names dynamically from r0 to select columns safely
species_names <- names(r0)
predictor_names <- c("x", "y", names(r1), names(r2_extended), continuous_names, "Pyrome_Consenso", "macro_pyro", "samp_bias")

final_modeling_data <- df_valid_pixels[, c(species_names, predictor_names)]
final_modeling_data <- na.omit(final_modeling_data)

write.csv(final_modeling_data, "data_for_modeling.csv", row.names = FALSE)