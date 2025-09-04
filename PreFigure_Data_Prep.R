library(rlang)
# Arrange Data -----------------------------------------------------------------
raw.dat <- read.csv('data/raw.discardHooks_data.csv')
at.subset3 <- raw.dat
raw.dat <- raw.dat%>%
  filter(gsw<0.6, gsw>-1)%>%
  filter(!is.na(Tair), !is.na(A), !is.na(Species))
weibull.discarded <- read.csv('weibull.discard.hooks.SANW.csv')
weibull.discarded <- weibull.discarded %>%
  mutate(Species = coalesce(Species, species))

TOPT = mean(weibull.discarded$T_opt)
hot.dat <- raw.dat %>%
  filter(Tleaf > TOPT)

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

#create a mean moisture summary-------------------------------------------------
tomst_by_elev <- tomst %>%
  group_by(Elevation) %>%
  summarise(mean_moist_pct = mean(moist_vol, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "tomst")

norway_by_elev <- norwayTemps %>%
  filter(str_to_lower(climate_variable) %in% c("soil_moisture", "soil moisture")) %>%
  mutate(
    value = pmax(value, 0),         # clamp negatives to 0
    value_pct = value * 100         # convert fraction → %
  ) %>%
  group_by(Elevation) %>%
  summarise(mean_moist_pct = mean(value_pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "norway")

## --- 3) Combine, then average across sources at the same Elevation (optional) ---
moist_by_elev_source <- bind_rows(tomst_by_elev, norway_by_elev)

# If you want ONE mean per Elevation (averaging TOMST + Norway where both exist):
moist_by_elev <- moist_by_elev_source %>%
  group_by(Elevation) %>%
  summarise(mean_moist_pct = mean(mean_moist_pct, na.rm = TRUE), .groups = "drop")

## --- 4) Classify Low / High using global median across elevations ---
## --- 4) Classify Low / Medium / High using tertiles across elevations ---
moist_by_elev <- moist_by_elev %>%
  mutate(
    moist_cat = cut(mean_moist_pct,
                    breaks = quantile(mean_moist_pct, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                    labels = c("Low","Medium","High"),
                    include.lowest = TRUE)
  )

## --- 5) Join onto your raw.dat so every row gets its elevation’s moisture class ---
dat_environ <- raw.dat %>%
  left_join(moist_by_elev, by = "Elevation") %>%
  filter(!is.na(moist_cat)) %>%
  mutate(moist_cat = factor(moist_cat, levels = c("Low","Medium","High")))


#create a mean temperature summary ---------------------------------------------
# ==============================
# 1) Air Temperature (already done)
# ==============================
# TOMST: AirTemp
tomst_air_by_elev <- tomst %>%
  group_by(Elevation) %>%
  summarise(mean_air = mean(AirTemp, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "tomst")

# Norway: air_temperature
norway_air_by_elev <- norwayTemps %>%
  filter(str_to_lower(climate_variable) %in% c("air_temperature", "air temp", "air temp.")) %>%
  group_by(Elevation) %>%
  summarise(mean_air = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "norway")

# Combine & classify
air_by_elev_source <- bind_rows(tomst_air_by_elev, norway_air_by_elev)

air_by_elev <- air_by_elev_source %>%
  group_by(Elevation) %>%
  summarise(mean_air = mean(mean_air, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    air_cat = cut(mean_air,
                  breaks = quantile(mean_air, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                  labels = c("Low", "Medium", "High"),
                  include.lowest = TRUE)
  )

dat.airTemp <- raw.dat %>%
  left_join(air_by_elev, by = "Elevation") %>%
  filter(!is.na(air_cat)) %>%
  mutate(air_cat = factor(air_cat, levels = c("Low", "Medium", "High")))

# ==============================
# 2) T1 / Soil Temperature
# ==============================
tomst_T1_by_elev <- tomst %>%
  group_by(Elevation) %>%
  summarise(mean_T1 = mean(T1, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "tomst")

norway_soil_by_elev <- norwayTemps %>%
  filter(str_to_lower(climate_variable) %in% c("soil_temperature", "soil temp")) %>%
  group_by(Elevation) %>%
  summarise(mean_T1 = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "norway")

T1_by_elev_source <- bind_rows(tomst_T1_by_elev, norway_soil_by_elev)

T1_by_elev <- T1_by_elev_source %>%
  group_by(Elevation) %>%
  summarise(mean_T1 = mean(mean_T1, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    T1_cat = cut(mean_T1,
                 breaks = quantile(mean_T1, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                 labels = c("Low", "Medium", "High"),
                 include.lowest = TRUE)
  )

dat.T1Temp <- raw.dat %>%
  left_join(T1_by_elev, by = "Elevation") %>%
  filter(!is.na(T1_cat)) %>%
  mutate(T1_cat = factor(T1_cat, levels = c("Low", "Medium", "High")))

# ==============================
# 3) T2 / Ground Temperature
# ==============================
tomst_T2_by_elev <- tomst %>%
  group_by(Elevation) %>%
  summarise(mean_T2 = mean(T2, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "tomst")

norway_ground_by_elev <- norwayTemps %>%
  filter(str_to_lower(climate_variable) %in% c("ground_temperature", "ground temp")) %>%
  group_by(Elevation) %>%
  summarise(mean_T2 = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "norway")

T2_by_elev_source <- bind_rows(tomst_T2_by_elev, norway_ground_by_elev)

T2_by_elev <- T2_by_elev_source %>%
  group_by(Elevation) %>%
  summarise(mean_T2 = mean(mean_T2, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    T2_cat = cut(mean_T2,
                 breaks = quantile(mean_T2, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                 labels = c("Low", "Medium", "High"),
                 include.lowest = TRUE)
  )

dat.T2Temp <- raw.dat %>%
  left_join(T2_by_elev, by = "Elevation") %>%
  filter(!is.na(T2_cat)) %>%
  mutate(T2_cat = factor(T2_cat, levels = c("Low", "Medium", "High")))

# ==============================
# Quick check
# ==============================
table(dat.airTemp$air_cat)
table(dat.T1Temp$T1_cat)
table(dat.T2Temp$T2_cat)
