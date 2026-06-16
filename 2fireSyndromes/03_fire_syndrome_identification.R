# load required libraries
library(mclust)
library(tidyverse)
library(caret)
library(clue)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(factoextra)

# data preparation
db <- read_csv("final_results/fire_vars_data.csv")

features <- db |> 
    select(-x, -y) |>
    na.omit()

# remove near zero variance columns
cols_to_remove <- nearZeroVar(features)
if(length(cols_to_remove) > 0) {
    features <- features[, -cols_to_remove]
}

scaled_features <- scale(features)

# run gmm model evaluation
gmm_model <- Mclust(scaled_features, G = 2:35, verbose = FALSE)
save(gmm_model, file = "gmm_model.RData")

png(file.path("bic_vs_clusters.png"), width = 1200, height = 1200, res = 100)
plot(gmm_model, what = "BIC")
abline(v = 11)
dev.off()

bic_values <- data.frame(gmm_model$BIC[1:34,])
bic_col <- bic_values$VEV[1:28]
clusters_seq <- 2:29
bic_df <- data.frame(clusters = clusters_seq, bic = bic_col)

# elbow analysis for cluster selection
k_start <- 2
k_elbow <- 11
k_peak <- 29

bic_start <- bic_df |> filter(clusters == k_start) |> pull(bic)
bic_elbow <- bic_df |> filter(clusters == k_elbow) |> pull(bic)
bic_peak <- bic_df |> filter(clusters == k_peak) |> pull(bic)

gain_phase1 <- bic_elbow - bic_start
steps_phase1 <- k_elbow - k_start
avg_gain_phase1 <- gain_phase1 / steps_phase1

gain_phase2 <- bic_peak - bic_elbow
steps_phase2 <- k_peak - k_elbow
avg_gain_phase2 <- gain_phase2 / steps_phase2

# save elbow analysis report to text file quietly
sink("analysis_elbow_clusters.txt")
cat("====================================================\n")
cat("QUANTITATIVE INFLECTION POINT ANALYSIS (ELBOW)\n")
cat("====================================================\n\n")
cat("--- Phase 1: Growth up to k =", k_elbow, "---\n")
cat("Total BIC Gain:", round(gain_phase1), "\n")
cat("Added clusters:", steps_phase1, "\n")
cat("AVG GAIN PER CLUSTER:", round(avg_gain_phase1), "\n\n")
cat("--- Phase 2: Growth from k =", k_elbow, "to k =", k_peak, "---\n")
cat("Total BIC Gain:", round(gain_phase2), "\n")
cat("Added clusters:", steps_phase2, "\n")
cat("AVG GAIN PER CLUSTER:", round(avg_gain_phase2), "\n\n")
cat("--- Verdict ---\n")
cat("Improvement ratio (Phase 1 vs Phase 2):", round(avg_gain_phase1 / avg_gain_phase2, 1), "times\n")
sink()

# consensus clustering preparation
k_final <- 11
n_runs <- 100

set.seed(1)
ref_model <- Mclust(scaled_features, G = k_final, verbose = FALSE)
ref_partition <- as.cl_partition(ref_model$classification)

aligned_probs_list <- list()

# 100 runs loop with real-time alignment
for (i in 1:n_runs) {
    set.seed(i)
    
    current_model <- tryCatch({
        Mclust(scaled_features, G = k_final, verbose = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(current_model)) {
        current_partition <- as.cl_partition(current_model$classification)
        
        agreement_table <- table(cl_class_ids(ref_partition), 
                                 cl_class_ids(current_partition))
        
        p <- solve_LSAP(agreement_table, maximum = TRUE)
        
        aligned_prob <- current_model$z[, p, drop = FALSE]
        aligned_probs_list[[i]] <- aligned_prob
    }
}

# average probabilities across runs
probs_array <- array(unlist(aligned_probs_list),
                     dim = c(nrow(aligned_probs_list[[1]]),
                             ncol(aligned_probs_list[[1]]),
                             length(aligned_probs_list)))

consensus_probs <- apply(probs_array, c(1, 2), mean)
colnames(consensus_probs) <- paste0("prob_pyrome_", 1:k_final)

consensus_classification <- apply(consensus_probs, 1, which.max)
max_probability <- apply(consensus_probs, 1, max)

db_consensus <- db |>
    mutate(
        pyrome_consensus = consensus_classification,
        prob_consensus = max_probability
    ) |>
    bind_cols(as.data.frame(consensus_probs))

write.csv(db_consensus, "db_pyromes_consensus.csv", row.names = FALSE)

# pca analysis for plotting
pca_result <- prcomp(scaled_features, center = TRUE, scale. = TRUE)

full_plot_data <- as.data.frame(pca_result$x[, 1:2]) |>
    mutate(
        cluster = as.factor(db_consensus$pyrome_consensus),
        max_prob = db_consensus$prob_consensus
    )

# filter high certainty points (>= 70%)
filtered_plot_data <- full_plot_data |>
    filter(max_prob >= 0.70)

# pca visualization
n_clusters <- length(unique(filtered_plot_data$cluster))
get_palette <- colorRampPalette(brewer.pal(min(n_clusters, 12), "Paired"))

png(file.path("pca_clusters_classic.png"), width = 7, height = 7, res = 300, units = "in")
ggplot(filtered_plot_data, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_d3(palette = "category20") +
    labs(
        x = "PC1 (37.3%)",
        y = "PC2 (23.6%)",
        color = "Clusters (Pyromes)"
    ) +
    theme_classic(base_size = 14) +
    theme(
        axis.text = element_text(color = "black"),
        legend.position = "right",
        axis.line = element_line(linewidth = 0.8, color = "black"),
        plot.margin = margin(10, 10, 10, 10)
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
dev.off()

# hierarchical clustering of centroids
data_with_clusters <- as.data.frame(scaled_features) |>
    mutate(cluster = as.factor(db_consensus$pyrome_consensus))

centroids <- data_with_clusters |>
    group_by(cluster) |>
    summarise(across(everything(), mean), .groups = 'drop')

cluster_labels <- paste0("P", centroids$cluster)
centroids_matrix <- centroids |> select(-cluster) |> as.matrix()
rownames(centroids_matrix) <- cluster_labels

dist_centroids <- dist(centroids_matrix, method = "euclidean")
hclust_centroids <- hclust(dist_centroids, method = "ward.D2")

# custom dendrogram theme and relations
theme_dendro <- theme_classic(base_family = "serif", base_size = 15) +
    theme(
        axis.line.y = element_line(linewidth = 1),
        plot.title = element_text(face = "bold")
    )

# define macro pyromes directly on the dataset
db_consensus$macro_pyro <- "M"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 8] <- "M1"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 10] <- "M1"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 6] <- "M2"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 9] <- "M2"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 11] <- "M3"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 1] <- "M4"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 5] <- "M4"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 2] <- "M5"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 7] <- "M5"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 3] <- "M5"
db_consensus$macro_pyro[db_consensus$pyrome_consensus == 4] <- "M5"

# build relations table for palette
pyro_relations <- db_consensus |> select(pyrome_consensus, macro_pyro) |> distinct()

base_colors <- c("Greens", "Reds", "BuPu", "YlOrBr", "Blues")
final_palette <- c()

unique_macros <- unique(pyro_relations$macro_pyro)

for (i in seq_along(unique_macros)) {
    m <- unique_macros[i]
    pyros_in_m <- pyro_relations |> filter(macro_pyro == m) |> pull(pyrome_consensus)
    n <- length(pyros_in_m)
    
    if(n == 1) {
        ramp_colors <- brewer.pal(9, base_colors[i])[6]
    } else if(n == 2) {
        ramp_colors <- colorRampPalette(brewer.pal(9, base_colors[i])[3:4])(n)
    } else if(n == 3) {
        ramp_colors <- colorRampPalette(brewer.pal(9, base_colors[i])[6:8])(n)
    } else if(n >= 4) {
        ramp_colors <- colorRampPalette(brewer.pal(9, base_colors[i])[3:6])(n)
    }
    
    names(ramp_colors) <- pyros_in_m
    final_palette <- c(final_palette, ramp_colors)
}

# dendrogram plot
p1 <- fviz_dend(
    hclust_centroids, k = 11, 
    cex = 1.3, fill = "red", k_colors = final_palette,
    color_labels_by_k = FALSE, ggtheme = theme_dendro,
    horiz = FALSE, main = "Pyrome hierarchical structure", 
    lwd = 1.5, ylab = "Height (Distance between centroids)",
    xlab = "Pyromes (Clusters)"
) + 
ylim(c(-0.5, 10.5)) + 
theme(
    plot.title = element_text(size = 22, hjust = 0, margin = margin(b = 15)),
    axis.title.x = element_text(size = 20, margin = margin(t = 12)),
    axis.title.y = element_text(size = 20, margin = margin(r = 12))
) +
geom_hline(yintercept = 4, color = "gray50", linetype = "dashed", linewidth = 0.5)

p1$layers[[2]]$aes_params$angle <- 45
p1$layers[[2]]$aes_params$colour <- "gray20"

dir.create("figs_paper", showWarnings = FALSE)
png(file.path("figs_paper/dendrogram_clusters.png"), width = 10, height = 7, units = 'in', res = 300)
print(p1)
dev.off()

# save final dataset with macro pyromes
write.csv(db_consensus, "db_pyromes_final.csv", row.names = FALSE)