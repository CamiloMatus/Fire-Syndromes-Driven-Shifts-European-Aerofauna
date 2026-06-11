# load required libraries
require(stringr)
library(dplyr)
require(ggplot2)

# read and compile data
file_list <- list.files("./var_imp")
rep1_files <- file_list[str_detect(file_list, "rep_1.csv")] 

macro_files <- rep1_files[str_detect(rep1_files, "var_imp_macro_pyro")]
pyro_files <- rep1_files[str_detect(rep1_files, "var_imp_pyro")]

# function to extract species name and process data
process_var_imp <- function(files, model_name) {
    temp_df <- data.frame()
    
    if(length(files) > 0) {
        for(i in 1:length(files)) {
            d_h <- read.csv(paste0("./var_imp/", files[i]))
            
            # extract species name removing prefixes and suffixes
            spp_name <- files[i] |>
                str_remove("^var_imp_") |>
                str_remove("macro_pyro_|pyro_") |>
                str_remove("_rep_1\\.csv$")
            
            m_h <- d_h |>
                group_by(names) |>             
                summarise(avg_imp = mean(imp), .groups = 'drop')  
            
            m_h$spp <- spp_name 
            m_h$model <- model_name
            
            temp_df <- bind_rows(temp_df, m_h)
        }
    }
    return(temp_df)
}

# compile both models
d_macro <- process_var_imp(macro_files, "Macro-Pyromes extended")
d_pyro <- process_var_imp(pyro_files, "Pyromes extended")

# merge all into a single dataframe
d_out <- bind_rows(d_macro, d_pyro)

# classification and clean names
d_out <- d_out |>
    mutate(
        predictor_group = case_when(
            grepl("^bio", names) ~ "Climate",
            grepl("^prop", names) ~ "Land use",
            grepl("^elevation", names) ~ "Topo.",
            grepl("^macro_pyro|^Pyrome_Consenso", names) ~ "Fire Syndromes",
            names == "samp_bias" ~ "Bias"
        ),
        predictor_group = factor(predictor_group, levels = c("Climate", "Topo.", "Land use", "Fire Syndromes", "Bias"))
    )

# clean x-axis names dynamically
d_out <- d_out |>
    mutate(names_clean = case_when(
        names == "elevation" ~ "Elevation",
        names == "prop_trees" ~ "Trees",
        names == "prop_shrubs" ~ "Shrubs",
        names == "prop_grasslands" ~ "Grasslands",
        names == "prop_crops" ~ "Croplands",
        names == "prop_built" ~ "Built-up",
        names == "prop_bare_soil" ~ "Bare soil",
        names == "prop_water" ~ "Water",
        names == "prop_wetlands" ~ "Wetlands",
        names == "samp_bias" ~ "Samp. Rate",
        grepl("^macro_pyro_", names) ~ sub("macro_pyro_", "", names), 
        grepl("^Pyrome_Consenso_", names) ~ sub("Pyrome_Consenso_", "P", names), 
        TRUE ~ names 
    ))

# define global order for all possible variables
x_order <- c("bio1", "bio2", "bio4", "bio8", "bio12", "bio15", "Elevation",
             "Trees", "Shrubs", "Grasslands", "Croplands", "Built-up", "Bare soil", "Water", "Wetlands",
             "M1", "M2", "M3", "M4", "M5", 
             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11",
             "Samp. Rate")

# apply order
d_out <- d_out |>
    mutate(
        names_clean = factor(names_clean, levels = x_order),
        model = factor(model, levels = c("Macro-Pyromes extended", "Pyromes extended"))
    )

# combined plot with shared axis
png("var_imp_combined_plot_bats.png", width = 12, height = 7, units = 'in', res = 300)

ggplot(d_out, aes(x = names_clean, y = avg_imp * 100, fill = model)) +
    geom_boxplot(color = "black", outlier.shape = 21, outlier.size = 1.2, alpha = 0.8, 
                 position = position_dodge2(preserve = "single")) +
    facet_grid(~ predictor_group, scales = "free_x", space = "free_x") +
    scale_fill_manual(
        values = c(
            "Macro-Pyromes extended" = "#4C9BE8", 
            "Pyromes extended" = "#FFA500"        
        )
    ) +
    labs(
        x = "Predictor variable",
        y = "Relative importance (%)",
        fill = "SDMs type", 
        title = "B) Bats "
    ) +
    theme_classic(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(color = "black", size = 18),
        title = element_text(color = "black", size = 18),
        legend.position = "bottom", 
        legend.title = element_text(color = "black"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
        strip.text = element_text(color = "black", size = 13), 
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black", linewidth = 0.8),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.spacing = unit(0.2, "lines")
    )

dev.off()

# save summary tables
dir.create("./final_tabs", showWarnings = FALSE)
write.csv(d_out, "./final_tabs/var_imp_bats.csv", row.names = FALSE)

df_summary <- d_out |>
    group_by(model, names_clean) |>
    summarise(
        min_val = min(avg_imp, na.rm = TRUE),
        mean_val = mean(avg_imp, na.rm = TRUE),
        median_val = median(avg_imp, na.rm = TRUE),
        max_val = max(avg_imp, na.rm = TRUE),
        sd_val = sd(avg_imp, na.rm = TRUE),
        .groups = "drop"
    ) 

write.csv(df_summary, "./final_tabs/var_imp_summary_bats.csv", row.names = FALSE)