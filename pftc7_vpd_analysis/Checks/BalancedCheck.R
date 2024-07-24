#This is to confirm the two step method gets the same results as if we run the fits simultaneously.

presence_absence <- dat %>%
  distinct(site, species) %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = species, values_from = presence, values_fill = list(presence = 0))

#Three balanced panels possible:
  #subset1: Achillea m, Alchemilla a, and vaccinium @ site 6,7,8,9
  #subset2: ecklonis and pallidum @ site 1,4
  #subset3: pallidum and senecio @ site 2,5

#Make the three balanced panel datasets:
# Recode site and species variables
discard.weibull <- read_csv("outputs/discard.weibull.SANW.csv")%>%
  mutate(country = case_when(curveID > 1000~"Norway",curveID < 1000~"SAfrica"))

subset1 <- discard.weibull %>%
  filter(country == "Norway" & Species != "Agrostis capillaris") %>%
  mutate(
    site = factor(site, levels = c(6, 7, 8, 9)),
    Species = factor(Species, levels = c("Achillea millefolium", "Alchemilla alpina", "Vaccinium vitis-idaea"))
  )
species_site_counts <- subset1 %>%
  group_by(Species, site) %>%
  summarise(count = n()) %>%
  ungroup()
## This is irreparably unbalanced.....


subset2 <- discard.weibull %>%
  filter(country == "SAfrica" & 
           site %in% c(1, 4) & 
           Species %in% c("Helichrysum ecklonis", "Helichrysum pallidum"))%>%
  mutate(site = factor(site, levels = c(1,4)))

balanced_subset2 <- subset2 %>%
  group_by(Species, site) %>%
  slice_sample(n = 3) %>%
  ungroup()


## Balance subset3::::
rows_to_delete <- discard.weibull %>%
  filter(country == "SAfrica" & 
           site %in% c(2, 5) & 
           Species == "Senecio glaberrimus") %>%
  mutate(site = factor(site, levels = c(2, 5)))

# Select 2 random rows from site 2
rows_to_delete_site2 <- rows_to_delete %>%
  filter(site == 2) %>%
  sample_n(2)

# Select 1 random row from site 5
rows_to_delete_site5 <- rows_to_delete %>%
  filter(site == 5) %>%
  sample_n(2)

# Combine rows to delete
rows_to_delete_final <- bind_rows(rows_to_delete_site2, rows_to_delete_site5)

discard.weibull <- discard.weibull %>%
  mutate(site = as.character(site))

rows_to_delete_final <- rows_to_delete_final %>%
  mutate(site = as.character(site))

# Create a custom identifier for each row in rows_to_delete_final
rows_to_delete_final <- rows_to_delete_final %>%
  rowwise() %>%
  mutate(row_id = paste(Species, site, T_opt, sep = "_")) %>%
  ungroup()

# Create a custom identifier for each row in discard.weibull
discard.weibull <- discard.weibull %>%
  rowwise() %>%
  mutate(row_id = paste(Species, site, T_opt, sep = "_")) %>%
  ungroup()

# Remove these rows from the original dataset using anti_join with the custom identifier
subset3 <- discard.weibull %>%
  anti_join(rows_to_delete_final, by = "row_id") %>%
  select(-row_id)  # Remove the custom identifier column

# Re-filter for the balanced dataset
subset3 <- subset3 %>%
  filter(country == "SAfrica" & 
           site %in% c(2, 5) & 
           Species %in% c("Senecio glaberrimus", "Helichrysum pallidum")) %>%
  mutate(site = factor(site, levels = c(2, 5)))

#Run some test linear regressions to see what I am working with:
summary(lmer(T_opt ~ Elevation + (1 | Species), data = subset2)) ## Shows no significance of species nor elevation
summary(lm(T_opt~factor(Species)+factor(Elevation), data=balanced_subset2)) ## Shows no significance of species
summary(lmer(T_opt~(1|Species)+Elevation, data=subset3))

summary(lmer(T_opt~Elevation + (1|country), data=discard.weibull))
summary(lm(T_opt~Species*country, data=discard.weibull))
