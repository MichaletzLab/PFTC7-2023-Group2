library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(assertr)
library(lubridate)
library(tidyr)

setwd("data")

#5E
#7:21-18:14 (only start at 8:30 (problem before that)) -- Data filtered to start after 8:30
flir5E <- read.csv('5E.4.24.24.csv')
data.start.5E = ISOdatetime(2023, 12, 09, 07, 21, 00)
filter_raw.5E <- names(flir5E) %in% str_subset(names(flir5E), "time_raw_[[:digit:]]")
flir5E <- flir5E[!filter_raw.5E] %>% select(-c("file"))
flir5E$datetime = data.start.5E + flir5E$time_elapsed
flir5E <- flir5E%>%
  mutate(site = "5")%>%
  mutate(Elevation = "2800")%>%
  mutate(Aspect = "E")%>%
  filter(hour(datetime) > 9 | (hour(datetime) == 9 & minute(datetime) >= 30))
#__Done

#1W
#7:06-18:08
flir1W <- read.csv("1W.4.24.24.csv")
data.start.1W = ISOdatetime(2023, 12, 12, 07, 06, 00)
filter_raw.1W <- names(flir1W) %in% str_subset(names(flir1W), "time_raw_[[:digit:]]")
flir1W <- flir1W[!filter_raw.1W] %>% select(-c("file"))
flir1W$datetime = data.start.1W + flir1W$time_elapsed
flir1W <- flir1W %>%
  mutate(site = "1") %>%
  mutate(Elevation = "2000")%>%
  mutate(Aspect = "W")%>%
  filter(hour(datetime) > 8 | (hour(datetime) == 8 & minute(datetime) >= 30))

#5W
#flir5W <- read.csv('5W_4.24.24.csv')
flir5W <- read.csv('5W.5.14.24.csv')
data.start.5W = ISOdatetime(2023, 12, 08, 09, 34, 00)
filter_raw.5W <- names(flir5W) %in% str_subset(names(flir5W), "time_raw_[[:digit:]]")
flir5W <- flir5W[!filter_raw.5W] %>% select(-c("file"))
flir5W$datetime = data.start.5W + flir5W$time_elapsed
flir5W <- flir5W %>%
  mutate(site = "5") %>%
  mutate(Elevation = "2800")%>%
  mutate(Aspect = "W")%>%
  filter(hour(datetime) > 11 | (hour(datetime) == 11 & minute(datetime) >= 20))

#1E
flir1EA <- read.csv("1Ea.2.24.24.csv")
data.start.1EA = ISOdatetime(2023, 12, 14, 09, 26, 00)
filter_raw.1EA <- names(flir1EA) %in% str_subset(names(flir1EA), "time_raw_[[:digit:]]")
flir1EA <- flir1EA[!filter_raw.1EA] %>% select(-c("file"))
flir1EA$datetime = data.start.1EA + flir1EA$time_elapsed
flir1EA <- flir1EA %>%
  mutate(site = "1") %>%
  mutate(Elevation = "2000")%>%
  mutate(Aspect = "E")

flir1EB <- read.csv("1Eb.4.24.24.csv")
data.start.1EB = ISOdatetime(2023, 12, 14, 13, 34, 00) ### Just fix this and then stick them together 
filter_raw.1EB <- names(flir1EB) %in% str_subset(names(flir1EB), "time_raw_[[:digit:]]")
flir1EB <- flir1EB[!filter_raw.1EB] %>% select(-c("file"))
flir1EB$datetime = data.start.1EB + flir1EB$time_elapsed
flir1EB <- flir1EB %>%
  mutate(site = "1") %>%
  mutate(Elevation = "2000")%>%
  mutate(Aspect = "E")

flir1E <- rbind(flir1EA, flir1EB)


all_columns <- unique(c(names(flir1E), names(flir1W), names(flir5E),names(flir5W)))

for (col in all_columns) {
  if (!col %in% names(flir1E)) flir1E[[col]] <- NA
  if (!col %in% names(flir1W)) flir1W[[col]] <- NA
  if (!col %in% names(flir5E)) flir5E[[col]] <- NA
  if (!col %in% names(flir5W)) flir5W[[col]] <- NA
}

FLIRdata <- rbind(flir1E, flir1W, flir5E, flir5W)





weibulls <- read_csv("weibull_results.csv")
as.numeric(weibulls$curveID)
at.meta$curveID = as.numeric(at.meta$curveID)
weib.meta <-left_join(weibulls, at.meta, by = 'curveID')


meantemp.1E <- flir1E %>%
  summarise(Tveg = mean(thermal_mean),
            sdTveg=mean(thermal_sd),
            Tair = mean(temp_atm))%>%
  mutate(site = "1") %>%
  mutate(Elevation = 2000)%>%
  mutate(Aspect = "E")

meantemp.1W <- flir1W %>%
  summarise(Tveg = mean(thermal_mean),
            sdTveg=mean(thermal_sd),
            Tair = mean(temp_atm))%>%
  mutate(site = "1") %>%
  mutate(Elevation = 2000)%>%
  mutate(Aspect = "W")

meantemp.5E <- flir5E %>%
  summarise(Tveg = mean(thermal_mean),
            sdTveg=mean(thermal_sd),
            Tair = mean(temp_atm))%>%
  mutate(site = "5")%>%
  mutate(Elevation = 2800)%>%
  mutate(Aspect = "E")

meantemp.5W <- flir5W %>%
  summarise(Tveg = mean(thermal_mean),
            sdTveg=mean(thermal_sd),
            Tair = mean(temp_atm))%>%
  mutate(site = "5")%>%
  mutate(Elevation = 2800)%>%
  mutate(Aspect = "W")
meantemps = bind_rows(meantemp.1E,meantemp.1W, meantemp.5E, meantemp.5W)


dat.4sites <- left_join(weib.meta, meantemps, by = c("site", "Aspect", "Elevation"))

#join weib.meta with tomst instead!
weib.tom <- left_join(weib.meta, mean_temp.tomst, by = 'site')
