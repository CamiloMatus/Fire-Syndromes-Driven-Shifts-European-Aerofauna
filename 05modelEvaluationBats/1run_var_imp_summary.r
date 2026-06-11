require(stringr)
library(dplyr)
require(ggplot2)

# --- 1. LECTURA Y COMPILACIÓN DE DATOS ---
l <- list.files("./varImp")
l2 <- l[str_detect(l, "rep_1.csv")] 

l_macro <- l2[str_detect(l2, "varImp_MacroPyros")]
l_pyro <- l2[str_detect(l2, "varImp_Pyromes")]

# Función modificada para extraer el nombre de la especie
procesar_varimp <- function(archivos, nombre_modelo) {
  d_temp <- data.frame()
  if(length(archivos) > 0) {
    for(i in 1:length(archivos)){
      d.h <- read.csv(paste0("./varImp/", archivos[i]))
      
      # Extraemos el nombre de la especie eliminando los prefijos y el sufijo
      # Se asume un formato como: varImp_MacroPyros_Nombre_Especie_rep_1.csv
      spp_name <- archivos[i] %>%
        str_remove("^varImp_") %>%
        str_remove("MacroPyros_|Pyromes_") %>%
        str_remove("_rep_1\\.csv$")
      
      m.h <- d.h %>%
        group_by(names) %>%             
        summarise(promedio_imp = mean(imp), .groups = 'drop')  
      
      m.h$spp <- spp_name # Guardamos el nombre de la especie extraído
      m.h$Model <- nombre_modelo
      
      d_temp <- bind_rows(d_temp, m.h)
    }
  }
  return(d_temp)
}

# Compilamos ambos modelos
cat("Procesando Macro Pyromes...\n")
d_macro <- procesar_varimp(l_macro, "Macro-Pyromes extended")

cat("Procesando Pyromes...\n")
d_pyro <- procesar_varimp(l_pyro, "Pyromes extended")

# Unimos todo en un solo dataframe
d.out <- bind_rows(d_macro, d_pyro)

# --- 2. CLASIFICACIÓN Y LIMPIEZA DE NOMBRES ---
# Creamos la columna de categoría (Grupo)
d.out <- d.out %>%
  mutate(
    Predictor_Group = case_when(
      grepl("^bio", names) ~ "Climate",
      grepl("^prop", names) ~ "Land use",
      grepl("^elevacion", names) ~ "Topo.",
      grepl("^macroPyro|^Pyrome_Consenso", names) ~ "Fire Syndromes",
      names == "sampBias" ~ "Bias"
    ),
    # Fijamos el orden lógico de los paneles con los nuevos nombres
    Predictor_Group = factor(Predictor_Group, levels = c("Climate", "Topo.", "Land use", "Fire Syndromes", "Bias"))
  )

# Limpiamos los nombres del eje X dinámicamente
d.out <- d.out %>%
  mutate(names_clean = case_when(
    names == "elevacion" ~ "Elevation",
    names == "prop_arboles" ~ "Trees",
    names == "prop_matorrales" ~ "Shrubs",
    names == "prop_pastizales" ~ "Grasslands",
    names == "prop_cultivos" ~ "Croplands",
    names == "prop_construido" ~ "Built-up",
    names == "prop_suelo_desnudo" ~ "Bare soil",
    names == "prop_agua" ~ "Water",
    names == "prop_humedales" ~ "Wetlands",
    names == "sampBias" ~ "Samp. Rate",
    grepl("^macroPyro_", names) ~ sub("macroPyro_", "", names), 
    grepl("^Pyrome_Consenso_", names) ~ sub("Pyrome_Consenso_", "P", names), 
    TRUE ~ names 
  ))

# Definimos el orden global de todas las posibles variables
orden_x <- c("bio1", "bio2", "bio4", "bio8", "bio12", "bio15", "Elevation",
             "Trees", "Shrubs", "Grasslands", "Croplands", "Built-up", "Bare soil", "Water", "Wetlands",
             "M1", "M2", "M3", "M4", "M5", 
             "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11",
             "Samp. Rate")

# Aplicamos el orden
d.out <- d.out %>%
  mutate(
    names_clean = factor(names_clean, levels = orden_x),
    Model = factor(Model, levels = c("Macro-Pyromes extended", "Pyromes extended"))
  )

# --- 3. GRÁFICO FINAL (UN SOLO EJE COMPARTIDO) ---
png("varImp_Combined_English_Bats.png", width = 12, height = 7, units = 'in', res = 300)

ggplot(d.out, aes(x = names_clean, y = promedio_imp * 100, fill = Model)) +
  # position_dodge2(preserve = "single") permite colocar cajas lado a lado sin estirar las que están solas
  geom_boxplot(color = "black", outlier.shape = 21, outlier.size = 1.2, alpha = 0.8, 
               position = position_dodge2(preserve = "single")) +
  
  # facet_grid divide el eje X en las categorías (Climate, Land Use, etc.)
  facet_grid(~ Predictor_Group, scales = "free_x", space = "free_x") +
  
  # Colores para diferenciar los modelos
  scale_fill_manual(
    values = c(
      "Macro-Pyromes extended" = "#4C9BE8", # Azul
      "Pyromes extended" = "#FFA500"        # Naranja
    )
  ) +
  
  labs(
    x = "Predictor variable",
    y = "Relative importance (%)",
    fill = "SDMs type", # Cambio solicitado
    title = "B) Bats "
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    # Tipografía sin bold
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black", size = 18),
    title = element_text(color = "black", size = 18),
    legend.position = "bottom", # Leyenda movida abajo
    legend.title = element_text(color = "black"),
    
    # Estilo de las cabeceras de los paneles
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(color = "black", size = 13), 
    
    # Bordes y marcas
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    
    # Ajuste fino: un pequeño espacio entre los paneles
    panel.spacing = unit(0.2, "lines")
  )

dev.off()


head(d.out)

write.csv(d.out,"./tabs_Final/varImp_bats.csv",row.names = F)


df_summary_2 <- d.out %>%
  group_by(Model,names_clean) %>%
  summarise(
    Min_val = min(promedio_imp, na.rm = TRUE),
    Mean_val = mean(promedio_imp, na.rm = TRUE),
    Median_val = median(promedio_imp, na.rm = TRUE),
    Max_val = max(promedio_imp, na.rm = TRUE),
    SD_val = sd(promedio_imp, na.rm = TRUE),
    .groups = "drop"
  ) 

write.csv(df_summary_2,"./tabs_Final/varImp_summary_bats.csv",row.names = F)
