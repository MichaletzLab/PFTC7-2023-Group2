# Purpose: Check correlations
# Plots:  plots of topt, breadth, Ea and Ed vs. air temperature as well as supplementary figures
# Dependencies: 
# Outputs: "TempSummary.csv"


###########      NORWAY        ###############
norwayTemps=read.csv("data/PFTC6_microclimate_2022.csv")

norwayTemps <- norwayTemps %>%
  mutate(datetime = ymd_hms(datetime))

# Apply the filtering
norwayTemps <- norwayTemps %>%
  filter(
    (month(datetime) == 7 & day(datetime) >= 23) | (month(datetime) == 8 & day(datetime) <= 1),
    hour(datetime) >= 8 & hour(datetime) <= 18,
    climate_variable == "air_temperature")%>%
  rename(AirTemp = value)%>%
  mutate(
    site = case_when(
      destSiteID == "Vikesland" ~ 6,
      destSiteID == "Hogsete" ~ 7,
      destSiteID == "Joasete" ~ 8,
      destSiteID == "Liahovden" ~ 9,
      TRUE ~ NA_integer_
    )
  )
# Summary dataframe of average daily temperatures for the time span of the Faster curves
summaryNorwayTemps = norwayTemps%>%
  group_by(site) %>%
  summarize(AirTemp = mean(AirTemp, na.rm = TRUE))


###########      S. Africa        ###############
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

# Summary dataframe of average daily temperatures for the time span of the Faster curves
summarySAfricaTemps_aspect <- tomst %>%
  group_by(site, aspect) %>%
  summarise(mean_tomst.Tair = mean(AirTemp, na.rm = TRUE))

summarySAfricaTemps <- tomst %>%
  group_by(site) %>%
  summarise(AirTemp = mean(AirTemp, na.rm = TRUE))
summarySAfricaTemps$site = as.numeric(summarySAfricaTemps$site)

###########      Combine files       ###############
TempSummary = bind_rows(summarySAfricaTemps,summaryNorwayTemps)
write.csv(TempSummary,"outputs/TempSummary.csv")

#####################################
#Check SAfrica data to see if site2 temperatures have high predictive power of T_opt and other parametners
#Using site2 since it is nearly the same location as measurements occured.

# Calculate average daily air temperature at site 2
temp_summary <- tomst %>%
  filter(site == "2") %>%
  group_by(date) %>%
  summarise(avg_air_temp = mean(AirTemp, na.rm = TRUE))

discard.weibull.2temps <- read_csv("outputs/discard.weibull.SANW.csv") %>%
  filter(site < 6) %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))  # Ensure Date is in Date format

# Join the average daily air temperature for site 2 with discard.weibull
discard.weibull.2temps <- discard.weibull.2temps %>%
  left_join(temp_summary, by = c("Date" = "date")) %>%
  rename(site2Temp = avg_air_temp)%>%
  filter(!is.na(site2Temp))

summary(lm(T_opt~site2Temp, data=discard.weibull.2temps))

