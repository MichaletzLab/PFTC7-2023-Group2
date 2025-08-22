#### Make VPD binned dataset ####
# make VPD bins (you can adjust number of bins with breaks = )
dat_VPD <- at.subset3 %>%
  mutate(vpd_bin = cut(VPDleaf,
                       breaks = pretty(VPDleaf, n = 10),   # ~10 bins
                       include.lowest = TRUE))%>%
  filter(gsw<0.6, gsw>-1)%>%
  filter(!is.na(Tair), !is.na(A), !is.na(Species))

# count number of obs per bin
vpd_counts <- dat_VPD %>%
  group_by(vpd_bin) %>%
  summarise(n = n(), .groups = "drop")

# create nice labels with ranges and counts
# (example: "0.3–1.2 (247)")
vpd_labels <- vpd_counts %>%
  mutate(label = paste0(as.character(vpd_bin), " (", n, ")")) %>%
  select(vpd_bin, label)

# join labels back to dataframe
dat_VPD <- dat_VPD %>%
  left_join(vpd_labels, by = "vpd_bin")

#### New binned VPD plot 8.21.25 ####
# plot like panel a: Photosynthesis vs AirTemp, colored by VPD bin
ggplot(dat_VPD, aes(x = Tleaf, y = A, color = label)) +
  geom_point(alpha = 0.01, size = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_viridis_d(name = "VPD (kPa)\n(n obs)") +
  labs(x = "Leaf temperature (°C)",
       y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1}))) +
  theme_classic() +
  theme(legend.position = "right") #+
  #facet_wrap(~Elevation) #expot as png 700x700

#### Plot the slope of the VPD categorized lines above ####
# Step 1: fit lm(A ~ Tleaf) within each Elevation × VPD bin
slopes <- dat_VPD %>%
  group_by(Elevation, vpd_bin, label) %>%
  do(broom::tidy(lm(A ~ Tleaf, data = .))) %>%
  ungroup() %>%
  filter(term == "Tleaf") %>%
  rename(slope = estimate)

# Step 2: get numeric midpoint of each VPD bin
# extract from factor levels like "(0.3,1.2]"
slopes <- slopes %>%
  mutate(vpd_mid = sapply(strsplit(as.character(vpd_bin), ","), function(x) {
    as.numeric(gsub("[^0-9\\.]", "", x[1])) + 
      (as.numeric(gsub("[^0-9\\.]", "", x[2])) - as.numeric(gsub("[^0-9\\.]", "", x[1]))) / 2
  }))

# Step 3: plot slopes vs VPD midpoint
ggplot(slopes, aes(x = vpd_mid, y = slope)) +
  #geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "VPD (kPa)",
       y = "Slope of A vs Tleaf",
       color = "Elevation") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~as.factor(Elevation)) # export as PNG 700x700


#### Plot just AT categorized by elevation ####
ggplot(dat_VPD2, aes(x = Tleaf, y = A, color = factor(Elevation))) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~Elevation) +
  labs(x = "Leaf temperature (°C)",
       y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
       color = "Elevation") +
  scale_color_viridis_d(name = "Elevation (masl)", direction = -1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")

#### Calculate slopes of A ~ Tleaf for each elevation ####
slopes_elev <- dat_VPD2 %>%
  group_by(Elevation) %>%
  do(tidy(lm(A ~ Tleaf, data = .))) %>%
  filter(term == "Tleaf") %>%
  select(Elevation, slope = estimate)

#### Plot slopes vs elevation ####
ggplot(slopes_elev, aes(x = factor(Elevation), y = slope)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Elevation", y = "Slope of A vs Tleaf") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
#### Plot the slopes of AT vs. elevation ####
#### Photosynthesis vs. Tleaf for categories of soil moisture ####
#### Fix up the tomst data for both countries ####
#Norway load data
norwayTemps=read.csv("data/PFTC6_microclimate_2022.csv")
norwayTemps <- norwayTemps %>%
  mutate(datetime = ymd_hms(datetime))%>%
  filter(
    (month(datetime) == 7 & day(datetime) >= 23) | (month(datetime) == 8 & day(datetime) <= 1),
    hour(datetime) >= 8 & hour(datetime) <= 18)%>%
  mutate(
    site = case_when(
      destSiteID == "Vikesland" ~ 6,
      destSiteID == "Hogsete" ~ 7,
      destSiteID == "Joasete" ~ 8,
      destSiteID == "Liahovden" ~ 9,
      TRUE ~ NA_integer_
    )
  )%>%
  mutate(
    Elevation = case_when(
      destSiteID == "Vikesland" ~ 469,
      destSiteID == "Hogsete" ~ 700,
      destSiteID == "Joasete" ~ 920,
      destSiteID == "Liahovden" ~ 1290,
      TRUE ~ NA_integer_
    ))
#SA load data:
tomst = read_csv("data/PFTC7_Tomst_Data.csv")
tomst <- tomst %>%
  mutate(Elevation = case_when(
    site == 1 ~ 2000,
    site == 2 ~ 2200,
    site == 3 ~ 2400,
    site == 4 ~ 2600,
    site == 5 ~ 2800))%>%
  mutate(Aspect = case_when(
    aspect == "east" ~ "E",
    aspect == "west" ~ "W"))%>%
  mutate(datetime = as.POSIXct(datetime, tz = "Europe/London"))%>%
  filter(!rowSums(is.na(.)))%>%
  filter(format(datetime, "%H:%M:%S") >= "08:00:00" &
           format(datetime, "%H:%M:%S") <= "18:00:00")%>%
  rename(AirTemp = T3)%>%
  mutate(site=as.character(site))%>%
  mutate(hourmin_str = format(datetime, "%Y-%m-%d %H:%M:%S"))

tomst_clean <- tomst %>%
  group_by(Elevation) %>%
  summarise(mean_moist = mean(moist, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(source = "tomst")
norway_clean <- norwayTemps %>%
  filter(climate_variable == "soil_moisture") %>%
  group_by(Elevation) %>%
  summarise(mean_moist = mean(value, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(source = "norwayTemps")
soil_moisture_summary <- bind_rows(tomst_clean, norway_clean)
soil_moisture_summary <- soil_moisture_summary %>%
  mutate(moist_cat = cut(mean_moist,
                         breaks = quantile(mean_moist, probs = c(0, 1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9, 1), na.rm = TRUE),
                         labels = c("SoilMoisture 0-11","SoilMoisture 11-22", "SoilMoisture 22-33", "SoilMoisture 33-44", "SoilMoisture 44-55", "SoilMoisture 55-66", "SoilMoisture 66-77", "SoilMoisture 77-88", "SoilMoisture 88-100"),
                         include.lowest = TRUE))
dat_VPD2 <- dat_VPD %>%
  left_join(soil_moisture_summary %>% select(Elevation, moist_cat), 
            by = "Elevation")
#### Make A vs. Tleaf plot at different soil moisture contents ####
annot_df <- dat_VPD2 %>%
  group_by(moist_cat) %>%
  summarise(Elevation_value = unique(Elevation), .groups = "drop") %>%
  # choose coordinates for the text in each panel
  mutate(x = min(dat_VPD2$Tleaf, na.rm = TRUE) + 10,
         y = max(dat_VPD2$A, na.rm = TRUE) - 1,
         label = paste0(Elevation_value,"masl"))

# Plot with annotations
ggplot(dat_VPD2, aes(x = Tleaf, y = A, color = moist_cat)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~moist_cat) +
  geom_text(data = annot_df, aes(x = x, y = y, label = label), 
            color = "black", inherit.aes = FALSE) +
  labs(x = "Leaf temperature (°C)",
       y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1}))) +
  scale_color_viridis_d(name = "Soil moisture level", direction = -1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")
slopes <- dat_VPD2 %>%
  group_by(moist_cat) %>%
  do(tidy(lm(A ~ Tleaf, data = .))) %>%
  filter(term == "Tleaf") %>%
  select(moist_cat, slope = estimate)
ggplot(slopes, aes(x = moist_cat, y = slope)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x="", y = "Slope of A vs Tleaf")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
