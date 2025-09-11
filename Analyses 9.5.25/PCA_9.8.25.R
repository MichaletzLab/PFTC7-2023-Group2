library(factoextra)

# Select the environmental variables for PCA
env_pca <- raw.env.data %>%
  select(mean_moist_pct, mean_air, mean_T1, mean_T2)

# Scale and run PCA
env_matrix <- scale(env_pca)

pca_res <- prcomp(env_matrix, center = TRUE, scale. = TRUE)

# Check summary
summary(pca_res)

# Scree plot
fviz_eig(pca_res)

# Biplot (color points by country for extra clarity)
fviz_pca_biplot(
  pca_res,
  repel = TRUE,
  col.var = "blue"#,
  #col.ind = raw.env.data$Country
)

# --- Extract PC scores and rejoin back to your main dataframe ---
pc_scores <- as.data.frame(pca_res$x)

raw.env.data_pca <- raw.env.data %>%
  bind_cols(pc_scores)

# Quick check
head(raw.env.data_pca)
