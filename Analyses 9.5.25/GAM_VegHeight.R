veg.height <- read.csv("data/ii_PFTC7_clean_elevationgradient_community_structure_2023.csv")
#Note that there is no veg height for Norway, just S. Africa.
veg_summary <- veg.height %>%
  rename(Elevation = elevation_m_asl) %>%
  filter(variable == "variable_height") %>%
  group_by(Elevation) %>%
  summarise(vegetation_height = mean(value, na.rm = TRUE),
            .groups = "drop")

raw.env.data <- raw.env.data %>%
  left_join(veg_summary, by = "Elevation")  

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