library(factoextra)
library(factoextra)
library(dplyr)
library(ggrepel)

country_colors <- c(
  "Norway"       = "#0072B2",
  "South Africa" = "#E69F00"
)

set.seed(123) # stm: I don't think this is necessary, but leaving so as to not break anything (far) downstream

# Use one row per site; the five env variables are site-level means (DataPrep.R),
# so distinct() collapses point-level replication and to give sites equal weight
env_site <- raw.env.data %>%
  distinct(Elevation, Country,
           mean_moist_pct, mean_air, mean_T1, mean_T2, vegetation_height)

env_pca <- env_site %>%
  select(mean_moist_pct, mean_air, mean_T1, mean_T2, vegetation_height)

# prcomp centers and scales once (drop the separate scale() call)
pca_res <- prcomp(env_pca, center = TRUE, scale. = TRUE)

# Pin PC1 orientation to cool-to-warm (positive air-temperature loading)
if (pca_res$rotation["mean_air", "PC1"] < 0) {
  pca_res$rotation[, "PC1"] <- -pca_res$rotation[, "PC1"]
  pca_res$x[, "PC1"]        <- -pca_res$x[, "PC1"]
}

# Pin PC2 orientation to dry-to-wet (positive VWC loading)
if (pca_res$rotation["mean_moist_pct", "PC2"] < 0) {
  pca_res$rotation[, "PC2"] <- -pca_res$rotation[, "PC2"]
  pca_res$x[, "PC2"]        <- -pca_res$x[, "PC2"]
}

# One score row per site
pc_scores <- as.data.frame(pca_res$x) %>%
  bind_cols(env_site %>% select(Elevation, Country))

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

# Correlation-biplot scaling: stretch each loading component by its PC's sd so
# arrow reach on each axis matches the point spread on that axis
sdev_pc <- pca_res$sdev[1:2]
loadings <- loadings %>%
  mutate(PC1 = PC1 * sdev_pc[1],
         PC2 = PC2 * sdev_pc[2])

# --- Data-driven arrow scale (PC2/VWC is the binding axis) ---
label_frac  <- 0.95
score_rng   <- c(diff(range(pc_scores$PC1)), diff(range(pc_scores$PC2)))
load_rng    <- c(diff(range(loadings$PC1)),  diff(range(loadings$PC2)))
arrow_scale <- label_frac * min(score_rng / load_rng)
message(sprintf("Resolved arrow_scale = %.2f", arrow_scale))

# --- Manual loading-label placement (deterministic; no leader segments) ---
# Anchor is the arrow tip; nx/ny nudge in PC-score units; tweak in ~0.05 steps.
lab_pos <- tibble::tribble(
  ~Variable,               ~hjust, ~vjust,  ~nx,    ~ny,
  "VWC",                    0,      0.5,     0.10,   0.00,   # right of tip
  "Soil\ntemperature",      0,      0,       0.08,   0.12,   # up-right, above the line
  "Surface\ntemperature",   0,      1,       0.08,  -0.12,   # down-right, below the line
  "Air\ntemperature",       0,      1,       0.05,  -0.12,   # below its tip
  "Vegetation height",      1,      0.5,    -0.10,  -0.05    # left of the downward arrow
)
loadings_lab <- loadings %>%
  left_join(lab_pos, by = "Variable") %>%
  mutate(xlab = PC1 * arrow_scale + nx,
         ylab = PC2 * arrow_scale + ny)

pca_plot <- ggplot() +
  geom_point(data = pc_scores,
             aes(x = PC1, y = PC2, color = Country),
             alpha = 0.5, size = 1) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = loadings_lab,
            aes(x = xlab, y = ylab, label = Variable),
            hjust = loadings_lab$hjust, vjust = loadings_lab$vjust,
            color = "black", size = 5, lineheight = 0.85) +
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
    x = sprintf("PC1 (%.1f%%)", summary(pca_res)$importance[2, 1] * 100),
    y = sprintf("PC2 (%.1f%%)", summary(pca_res)$importance[2, 2] * 100)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

pca_plot

# --- Rejoin site-level PC scores to point-level data by Elevation ---
site_scores <- env_site %>%
  select(Elevation) %>%
  bind_cols(as.data.frame(pca_res$x))

raw.env.data_pca <- raw.env.data %>%
  left_join(site_scores, by = "Elevation")

