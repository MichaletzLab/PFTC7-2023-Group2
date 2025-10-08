veg.height.SA <- read.csv("data/ii_PFTC7_clean_elevationgradient_community_structure_2023.csv")
veg.height.N <- read.csv("data/PFTC6_ElevationalGradient_vegetation_height_2019_2022.csv")
veg_summary <- veg.height.SA %>%
  rename(Elevation = elevation_m_asl) %>%
  filter(variable == "vegetation_height") %>%  # ✅ fixed name
  group_by(Elevation) %>%
  summarise(
    vegetation_height = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

#Norway:
elev_map <- c(
  Vik = 469,
  Hog = 700,
  Joa = 920,
  Lia = 1290
)

veg_summary_N <- veg.height.N %>%
  mutate(site_prefix = substr(siteID, 1, 3)) %>%
  mutate(Elevation = elev_map[site_prefix]) %>%
  group_by(Elevation) %>%
  summarise(vegetation_height = mean(vegetation_height, na.rm = TRUE)) %>%
  arrange(Elevation)

veg_summary_N

combined_summary <- bind_rows(veg_summary, veg_summary_N)
# Join veg height into your environmental data
raw.env.data <- raw.env.data %>%
  left_join(combined_summary, by = "Elevation")



#summary of GAM without veg height:
summary(gam_mod.A_species) #from document GAM_Elevation_9.8.25.R

#summary of GAM with veg height
veg.GAM <- gam(A ~ s(Tleaf, k = 3) +
                 s(Elevation, k = 3) +
                 s(vegetation_height, k=3) +
                 s(mean_T2, k = 3) +
                 s(mean_moist_pct, k = 3) +
                 Species +
                 ti(Tleaf, Elevation, k = 3) +
                 ti(Tleaf, vegetation_height, k = 3) +
                 ti(Tleaf, mean_T2, k = 3) +
                 ti(Tleaf, mean_moist_pct, k = 3) +
                 s(Tleaf, by = Species, k = 3),
               data = raw.env.data,
               method = "REML"
)
summary(veg.GAM)

AIC(gam_mod.A_species, veg.GAM)# Vegetation height seems important to include but we do not have it for Norway...

