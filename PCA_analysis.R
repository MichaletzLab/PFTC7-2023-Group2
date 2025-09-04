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

#Extract PC1 as a predictor variable 
# PC scores for each site (rows = sites, columns = PCs)
pc_scores <- as.data.frame(pca_res$x)

# Add back the Elevation/site ID
pc_scores <- env_pca %>%
  select(Elevation) %>%
  bind_cols(pc_scores)

head(pc_scores)
dat_with_PC1 <- raw.dat %>%
  left_join(pc_scores %>% select(Elevation, PC1), by = "Elevation")

gam_A <- gam(A ~ s(Tleaf, k=5) + s(PC1, k=5), data = dat_with_PC1)
summary(gam_A)
gam_Ab <- gam(A ~ s(Tleaf, k=5) + s(mean_air, k=5), data = dat.airTemp)
summary(gam_Ab)





# Now the same but include mins and maxes--------------------------------------
## --- Soil moisture ----------------------------------------------------------
tomst_by_elev <- tomst %>%
  group_by(Elevation) %>%
  summarise(
    mean_moist_pct = mean(moist_vol, na.rm = TRUE),
    min_moist_pct  = min(moist_vol, na.rm = TRUE),
    max_moist_pct  = max(moist_vol, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "tomst")

norway_by_elev <- norwayTemps %>%
  filter(str_to_lower(climate_variable) %in% c("soil_moisture", "soil moisture")) %>%
  mutate(value = pmax(value, 0), value_pct = value * 100) %>%
  group_by(Elevation) %>%
  summarise(
    mean_moist_pct = mean(value_pct, na.rm = TRUE),
    min_moist_pct  = min(value_pct, na.rm = TRUE),
    max_moist_pct  = max(value_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "norway")

moist_by_elev <- bind_rows(tomst_by_elev, norway_by_elev) %>%
  group_by(Elevation) %>%
  summarise(
    mean_moist_pct = mean(mean_moist_pct, na.rm = TRUE),
    min_moist_pct  = mean(min_moist_pct,  na.rm = TRUE),
    max_moist_pct  = mean(max_moist_pct,  na.rm = TRUE),
    .groups = "drop"
  )

## --- Air temperature --------------------------------------------------------
tomst_air_by_elev <- tomst %>%
  group_by(Elevation) %>%
  summarise(
    mean_air = mean(AirTemp, na.rm = TRUE),
    min_air  = min(AirTemp, na.rm = TRUE),
    max_air  = max(AirTemp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "tomst")

norway_air_by_elev <- norwayTemps %>%
  filter(str_to_lower(climate_variable) %in% c("air_temperature", "air temp", "air temp.")) %>%
  group_by(Elevation) %>%
  summarise(
    mean_air = mean(value, na.rm = TRUE),
    min_air  = min(value, na.rm = TRUE),
    max_air  = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "norway")

air_by_elev <- bind_rows(tomst_air_by_elev, norway_air_by_elev) %>%
  group_by(Elevation) %>%
  summarise(
    mean_air = mean(mean_air, na.rm = TRUE),
    min_air  = mean(min_air,  na.rm = TRUE),
    max_air  = mean(max_air,  na.rm = TRUE),
    .groups = "drop"
  )

## --- Soil T1 
tomst_T1_by_elev <- tomst %>%
  group_by(Elevation) %>%
  summarise(
    mean_T1 = mean(T1, na.rm = TRUE),
    min_T1  = min(T1, na.rm = TRUE),
    max_T1  = max(T1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "tomst")

norway_soil_by_elev <- norwayTemps %>%
  filter(str_to_lower(climate_variable) %in% c("soil_temperature", "soil temp")) %>%
  group_by(Elevation) %>%
  summarise(
    mean_T1 = mean(value, na.rm = TRUE),
    min_T1  = min(value, na.rm = TRUE),
    max_T1  = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "norway")

T1_by_elev <- bind_rows(tomst_T1_by_elev, norway_soil_by_elev) %>%
  group_by(Elevation) %>%
  summarise(
    mean_T1 = mean(mean_T1, na.rm = TRUE),
    min_T1  = mean(min_T1,  na.rm = TRUE),
    max_T1  = mean(max_T1,  na.rm = TRUE),
    .groups = "drop"
  )

## --- Ground T2 
tomst_T2_by_elev <- tomst %>%
  group_by(Elevation) %>%
  summarise(
    mean_T2 = mean(T2, na.rm = TRUE),
    min_T2  = min(T2, na.rm = TRUE),
    max_T2  = max(T2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "tomst")

norway_ground_by_elev <- norwayTemps %>%
  filter(str_to_lower(climate_variable) %in% c("ground_temperature", "ground temp")) %>%
  group_by(Elevation) %>%
  summarise(
    mean_T2 = mean(value, na.rm = TRUE),
    min_T2  = min(value, na.rm = TRUE),
    max_T2  = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "norway")

T2_by_elev <- bind_rows(tomst_T2_by_elev, norway_ground_by_elev) %>%
  group_by(Elevation) %>%
  summarise(
    mean_T2 = mean(mean_T2, na.rm = TRUE),
    min_T2  = mean(min_T2,  na.rm = TRUE),
    max_T2  = mean(max_T2,  na.rm = TRUE),
    .groups = "drop"
  )


## --- Build PCA input 
moist_pca <- moist_by_elev %>%
  select(Elevation, mean_moist_pct, min_moist_pct, max_moist_pct)

air_pca <- air_by_elev %>%
  select(Elevation, mean_air, min_air, max_air)

T1_pca <- T1_by_elev %>%
  select(Elevation, mean_T1, min_T1, max_T1)

T2_pca <- T2_by_elev %>%
  select(Elevation, mean_T2, min_T2, max_T2)

env_pca <- moist_pca %>%
  left_join(air_pca, by = "Elevation") %>%
  left_join(T1_pca, by = "Elevation") %>%
  left_join(T2_pca, by = "Elevation")

## --- PCA --------------------------------------------------------------------
env_matrix <- env_pca %>%
  #select(-Elevation) %>%
  scale()

pca_res <- prcomp(env_matrix, center = TRUE, scale. = TRUE)

summary(pca_res)

library(factoextra)

# Scree plot
fviz_eig(pca_res)

# Biplot
fviz_pca_biplot(
  pca_res,
  repel = TRUE,
  col.var = "blue",   # variable arrows
  col.ind = "red"     # site points
)

# Run on physiological data as well: -------------------------------------------
library(dplyr)

# Join environmental summaries onto your physiological data
dat_pca <- raw.dat %>%
  left_join(air_by_elev %>% select(Elevation, mean_air), by = "Elevation") %>%
  left_join(T1_by_elev %>% select(Elevation, mean_T1), by = "Elevation") %>%
  left_join(T2_by_elev %>% select(Elevation, mean_T2), by = "Elevation") %>%
  left_join(moist_by_elev %>% select(Elevation, mean_moist_pct), by = "Elevation") %>%
  select(A, E, gsw, mean_air, mean_T1, mean_T2, mean_moist_pct, Elevation)
# Scale numeric variables
pca_matrix <- dat_pca %>% select(-Elevation) %>% scale()  # optionally include Elevation

pca_res <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)

summary(pca_res)

# Scree plot
library(factoextra)
fviz_eig(pca_res)

# Biplot: points colored by species (if you want to see species patterns)
dat_pca$Species <- raw.dat$Species  # add back species for plotting

fviz_pca_biplot(pca_res,
                repel = TRUE,
                col.var = "blue",   # arrows for variables
                col.ind = dat_pca$Species)    # points colored by species

