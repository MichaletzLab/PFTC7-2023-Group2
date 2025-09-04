species_lookup_raw <- tibble::tribble(
  ~raw,                  ~clean,
  "Achillea millefolium", "Achillea millefolium",
  "Agrostis capillaris",  "Agrostis capillaris",
  "Alchemilla alpina",    "Alchemilla alpina",
  "Vaccinium vitis-idaea","Vaccinium vitis-idaea",
  "dimorphotheca_jucunda","Dimorphotheca jucunda",
  "eucomis_cf_humilis",   "Eucomis bicolor",
  "helichrysum_ecklonis", "Helichrysum ecklonis",
  "helichrysum_nudifolium","Helichrysum nudifolium",
  "helichrysum_pallidum", "Helichrysum pallidum",
  "helichrysum_pilosellum","Helichrysum piloselum",
  "senecio_glaberrimus",  "Senecio glaberrimus",
  "senecio_tall",         "Senecio tall"
)
species_lookup_weib <- tibble::tribble(
  ~raw,                  ~clean,
  "Achillea", "Achillea millefolium",
  "Agrostis",  "Agrostis capillaris",
  "Alchemilla",    "Alchemilla alpina",
  "Vaccinium","Vaccinium vitis-idaea",
  "dimorphotheca_jucunda","Dimorphotheca jucunda",
  "eucomis_cf_humilis",   "Eucomis bicolor",
  "helichrysum_ecklonis", "Helichrysum ecklonis",
  "helichrysum_nudifolium","Helichrysum nudifolium",
  "helichrysum_pallidum", "Helichrysum pallidum",
  "helichrysum_pilosellum","Helichrysum piloselum",
  "senecio_glaberrimus",  "Senecio glaberrimus",
  "senecio_tall",         "Senecio tall"
)

raw.dat <- raw.dat %>%
  left_join(species_lookup_raw, by = c("Species" = "raw")) %>%
  mutate(Species = coalesce(clean, Species)) %>%
  select(-clean)
weibull.discarded <- weibull.discarded %>%
  left_join(species_lookup_weib, by = c("Species" = "raw")) %>%
  mutate(Species = coalesce(clean, Species)) %>%
  select(-clean)


mean_Topt <- weibull.discarded %>%
  group_by(Species) %>%
  summarise(mean_Topt = mean(T_opt, na.rm = TRUE))
mean_iWUE <- raw.dat %>%
  mutate(iWUE = A / gsw) %>%
  group_by(Species) %>%
  summarise(mean_iWUE = mean(iWUE, na.rm = TRUE))
species_summary <- mean_Topt %>%
  left_join(mean_iWUE, by = "Species")


ggplot(species_summary, aes(x = mean_Topt, y = mean_iWUE, label = Species)) +
  geom_point(size = 3, color = "darkblue") +
  geom_text(vjust = -0.5, size = 3) +
  labs(
    x = expression(Mean ~ T[opt] ~ (degree*C)),
    y = expression(Mean ~ iWUE ~ (mu*mol~CO[2]~mol^{-1}~H[2]*O))
  ) +
  theme_classic()


# Now group by curveID
iWUE_curveID <- raw.dat %>%
  mutate(iWUE = A / gsw) %>%
  group_by(curveID) %>%
  summarise(mean_iWUE = mean(iWUE, na.rm = TRUE))
iWUE_curve <- weibull.discarded %>%
  left_join(iWUE_curveID, by = "curveID")
iWUE_curve <- iWUE_curve%>%
  mutate(Species = coalesce(Species,species))

# Fit linear model
lm_fit <- lm(mean_iWUE ~ T_opt, data = iWUE_curve)

# Extract slope and p-value
lm_summary <- summary(lm_fit)
slope <- coef(lm_summary)[2, 1]
pval <- coef(lm_summary)[2, 4]

# Format nicely
pval_text <- paste0("Slope = ", round(slope, 3),
                    "\nP = ", signif(pval, 3))

# Plot with annotation
ggplot(iWUE_curve, aes(x = T_opt, y = mean_iWUE, color = Species)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(
    x = expression(T[opt] ~ (degree*C)),
    y = expression(Mean ~ iWUE ~ (mu*mol~CO[2]~mol^{-1}~H[2]*O))
  ) +
  annotate("text",
           x = min(iWUE_curve$T_opt, na.rm = TRUE) + 1,   # adjust placement
           y = max(iWUE_curve$mean_iWUE, na.rm = TRUE), 
           label = pval_text,
           hjust = 0,
           size = 4) +
  theme_classic()


