library(factoextra)
library(factoextra)
library(dplyr)
library(ggrepel)

set.seed(123)
sample_idx <- sample(1:nrow(raw.env.data), 1000)  # sample points

# Select environmental variables for PCA
env_pca <- raw.env.data %>%
  select(mean_moist_pct, mean_air, mean_T1, mean_T2, vegetation_height)

# Scale and run PCA
env_matrix <- scale(env_pca)
pca_res <- prcomp(env_matrix, center = TRUE, scale. = TRUE)

# Extract PC scores and subset to sampled points
pc_scores <- as.data.frame(pca_res$x) %>%
  bind_cols(raw.env.data %>% select(Elevation)) %>%
  slice(sample_idx)  # keep only sampled points

# Identify unique elevations among sampled points
unique_elevs <- pc_scores %>%
  group_by(Elevation) %>%
  slice(1)  # first occurrence of each Elevation

### Plotting ####
loadings <- as.data.frame(pca_res$rotation) %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(Variable = recode(Variable,
                           mean_moist_pct = "Soil moisture",
                           mean_T1 = "Soil temp",
                           mean_T2 = "Ground temp",
                           mean_air = "Air temp",
                           vegetation_height = "VegHeight"
  ))

# --- Create biplot manually for full control ---
pca_plot <- ggplot() +
  # points
  geom_point(data = pc_scores,
             aes(x = PC1, y = PC2),
             color = "grey30", alpha = 0.5, size = 1) +
  # arrows
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 0.8) +
  # arrow labels
  geom_text_repel(data = loadings,
                  aes(x = PC1 * 3.2, y = PC2 * 3.2, label = Variable),
                  color = "black", size = 4, fontface = "italic") +
  # elevation labels
  geom_text_repel(data = unique_elevs,
                  aes(x = PC1, y = PC2, label = Elevation),
                  color = "black", size = 4,
                  box.padding = 0.3,
                  segment.color = "grey50") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  labs(
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

pca_plot

# --- Extract PC scores and rejoin back to main dataframe ---
pc_scores <- as.data.frame(pca_res$x)

raw.env.data_pca <- raw.env.data %>%
  bind_cols(pc_scores)
