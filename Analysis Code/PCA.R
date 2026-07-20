library(factoextra)
library(factoextra)
library(dplyr)
library(ggrepel)

country_colors <- c(
  "Norway"       = "#0072B2",
  "South Africa" = "#E69F00"
)

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
  bind_cols(raw.env.data %>% select(Elevation, Country)) %>% # include colors
  slice(sample_idx)  # keep only sampled points

# Identify unique elevations among sampled points
unique_elevs <- pc_scores %>%
  group_by(Elevation) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Elev_label = if_else(Elevation < 1500,
                              paste0("N", round(Elevation)),
                              paste0("SA", round(Elevation))))

### Plotting ####
loadings <- as.data.frame(pca_res$rotation) %>%
  tibble::rownames_to_column("Variable") %>%
  mutate(Variable = dplyr::recode(
    Variable,
    "mean_moist_pct"    = "VWC",
    "mean_T1"           = "Soil\ntemperature",
    "mean_T2"           = "Surface\ntemperature",
    "mean_air"          = "Air\ntemperature",
    "vegetation_height" = "Vegetation height"
  ))


# --- Create biplot manually for full control ---
arrow_scale <- 4.5
loadings_repel  <- loadings %>% filter(Variable != "Vegetation height")
loadings_manual <- loadings %>% filter(Variable == "Vegetation height")

pca_plot <- ggplot() +
  geom_point(data = pc_scores,
             aes(x = PC1, y = PC2, color = Country),
             alpha = 0.5, size = 1) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text_repel(data = loadings_repel,
                  aes(x = PC1 * (arrow_scale + 0.3), y = PC2 * (arrow_scale + 0.3), label = Variable),
                  color = "black", size = 5, lineheight = 0.85,
                  box.padding = 0.15, point.padding = 0.1,
                  force = 0.5, max.overlaps = Inf,
                  segment.color = "grey40", segment.size = 0.3,
                  seed = 42) +
  # Vegetation height: manual placement, pulled left/down away from N700
  geom_text(data = loadings_manual,
            aes(x = PC1 * (arrow_scale + 0.3), y = PC2 * (arrow_scale + 0.3), label = Variable),
            color = "black", size = 5, lineheight = 0.85,
            hjust = 0, nudge_x = 0.1, nudge_y = -0.2) +
  geom_text_repel(data = unique_elevs,
                  aes(x = PC1, y = PC2, label = Elev_label, color = Country),
                  size = 5,
                  box.padding = 0.4, point.padding = 0.2,
                  force = 3, max.overlaps = Inf,
                  segment.color = "grey50",
                  seed = 42) +
  scale_color_manual(values = country_colors, name = "") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  theme_classic(base_size = 18) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

pca_plot

# --- Extract PC scores and rejoin back to main dataframe ---
pc_scores <- as.data.frame(pca_res$x)

raw.env.data_pca <- raw.env.data %>%
  bind_cols(pc_scores)

