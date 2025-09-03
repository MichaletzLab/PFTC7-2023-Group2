# Soil moisture
moist_pca <- moist_by_elev %>%
  select(Elevation, mean_moist_pct)

# Air temperature
air_pca <- air_by_elev %>%
  select(Elevation, mean_air)

# Soil T1
T1_pca <- T1_by_elev %>%
  select(Elevation, mean_T1)

# Soil T2
T2_pca <- T2_by_elev %>%
  select(Elevation, mean_T2)

env_pca <- moist_pca %>%
  left_join(air_pca, by = "Elevation") %>%
  left_join(T1_pca, by = "Elevation") %>%
  left_join(T2_pca, by = "Elevation")

# Check the structure
str(env_pca)

env_matrix <- env_pca %>%
  select(-Elevation) %>%
  scale() 

pca_res <- prcomp(env_matrix, center = TRUE, scale. = TRUE)
summary(pca_res)

library(factoextra)

# Scree plot
fviz_eig(pca_res)

# Biplot (sites + environmental variables)
fviz_pca_biplot(pca_res,
                repel = TRUE,
                col.var = "blue",   # arrows for variables
                col.ind = "red")    # points for sites
