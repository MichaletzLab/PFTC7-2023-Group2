library(ggplot2)
library(dplyr)
library(fuzzyjoin)
library(lmerTest)
library(MuMIn)
library(car)
library(lubridate)

# Extract time component from datetime
FLIRdata$time_of_day <- format(FLIRdata$datetime, "%H:%M:%S")
#Remove temp_atm values for 1E
FLIRdata = FLIRdata%>%
  mutate(temp_atm = case_when(
    site == 1 & Aspect == "E" ~ NA_real_,
    TRUE ~ temp_atm))


tomst = read_csv("PFTC7_Tomst_Data.csv")
tomstF <- tomst %>%
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
  filter((as.Date(datetime) == as.Date("2023-12-12") & site==1 & Aspect =="W")|
           (as.Date(datetime) == as.Date("2023-12-14") & site==1 & Aspect =="E")|
           (as.Date(datetime) == as.Date("2023-12-08") & site==5 & Aspect =="W")|
           (as.Date(datetime) == as.Date("2023-12-09") & site==5 & Aspect =="E"))%>%
  filter(format(datetime, "%H:%M:%S") >= "08:00:00" &
           format(datetime, "%H:%M:%S") <= "19:00:00")%>%
  rename(tomst.Tair = T3)%>%
  mutate(site=as.character(site))%>%
  mutate(hourmin_str = format(datetime, "%Y-%m-%d %H:%M:%S"))


#testplot2 = tomstF%>%
#  filter(site==1, Aspect=="E")%>%
#  ggplot(aes(y=tomst.Tair, x=datetime)) +geom_point()
#testplot2

ave.day.tomst <- tomstF %>%
  group_by(site, aspect) %>%
  summarise(mean_tomst.Tair = mean(tomst.Tair, na.rm = TRUE))

FLIRt <- FLIRdata %>%
  mutate(datetime = as.POSIXct(datetime))%>%
  group_by(site, Aspect) %>%
  mutate(hourmin = as.POSIXct(round(as.numeric(datetime) / 600) * 600, origin = "1970-01-01")) %>%
  group_by(hourmin) %>%
  summarise(across(where(is.numeric), mean,na.rm = TRUE),
            across(where(is.factor), first),
            across(where(is.character), first)) %>%
  rename(datetime = hourmin)%>%
  ungroup()%>%
  mutate(Elevation = as.numeric(Elevation))%>%
  mutate(site=as.character(site))%>%
  mutate(hourmin_str = format(datetime, "%Y-%m-%d %H:%M:%S"))

FLIRtomst <- left_join(FLIRt, tomstF, by = c("hourmin_str", "site", "Aspect", "Elevation"))%>%
  select(-datetime.x, -datetime.y, -aspect, -plot)%>%
  mutate(temp_atm = case_when(
    site == 1 & Aspect == "E" ~ NA_real_,TRUE ~ temp_atm))%>%
  rename(flir.Tleaf = thermal_mean)%>%
  rename(flir.Tair = temp_atm)%>% 
  mutate(flir.Tair = flir.Tair - 273.15)%>%
  mutate(flir.Tleaf = flir.Tleaf - 273.15)%>%
  mutate(hourmin_str = ymd_hms(hourmin_str)) %>%
  mutate(TOD = format(hourmin_str, "%H:%M:%S"))%>%
  mutate(TOD_seconds = as.numeric(strptime(TOD, format = "%H:%M:%S")) - as.numeric(strptime("00:00:00", format = "%H:%M:%S")))



#testplot3 = tomstF%>%
#  filter(site==1, Aspect=="E")%>%
#  ggplot(aes(y=tomst.Tair, x=hourmin_str)) +geom_point()
#testplot3


#FLIRtomst <- FLIRtomst %>%
#  group_by(hourmin) %>%
#  summarise_all(~ if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))

################################
################################
FLIR.dt = FLIRt%>%
  mutate(datetime = hourmin_str)
tomst = tomst%>%
  mutate(datetime = as.character(datetime))%>%
  mutate(site=as.character(site))%>%
  rename(Aspect=aspect)%>%
  mutate(Elevation = case_when(
    site == 1 ~ 2000,
    site == 2 ~ 2200,
    site == 3 ~ 2400,
    site == 4 ~ 2600,
    site == 5 ~ 2800))%>%
  mutate(Aspect = case_when(
    Aspect == "east" ~ "E",
    Aspect == "west" ~ "W"))

tomstFLIR <- left_join(tomst, FLIR.dt, by = c("datetime", "site","Aspect", "Elevation"))

tomstFLIR <- tomstFLIR %>%
  rename(flir.Tleaf = thermal_mean)%>%
  rename(flir.Tair = temp_atm)%>% 
  rename(tomst.Tair = T3)%>%
  mutate(flir.Tair = case_when(
    site == 1 & Aspect == "E" ~ NA_real_,
    TRUE ~ flir.Tair))%>%
  mutate(flir.Tair = flir.Tair - 273.15)%>%
  mutate(flir.Tleaf = flir.Tleaf - 273.15)%>%
  mutate(hourmin_str = ymd_hms(hourmin_str)) %>%
  mutate(datetime = ymd_hms(datetime)) %>%
  filter(format(datetime, "%H:%M:%S") >= "08:30:00" &
           format(datetime, "%H:%M:%S") <= "18:20:00")%>%
  mutate(TOD = format(datetime, "%H:%M:%S"))%>%
  mutate(TOD_seconds = as.numeric(strptime(TOD, format = "%H:%M:%S")) - as.numeric(strptime("00:00:00", format = "%H:%M:%S")))

 
# Create a combined variable for Elevation and Aspect
FLIRdata <- FLIRdata %>%
  mutate(Elev_Aspect = interaction(Elevation, Aspect, sep = " - "))%>%
  mutate(time_of_day = as.POSIXct(time_of_day, format = "%H:%M"))
FLIRtomst <- FLIRtomst %>%
  mutate(Elev_Aspect = interaction(Elevation, Aspect, sep = " - "))%>%
  mutate(time_of_day = as.POSIXct(hourmin_str, format = "%H:%M"))

# Define custom colors for each combination of Elevation and Aspect
custom_colors <- c("2000 - E" = "darkblue", "2000 - W" = "#0077BB",
                   "2800 - E" = "#882255", "2800 - W" = "#AA4499")

# part of FIGURE 3
tleaf.time <- FLIRdata %>%
  ggplot(aes(x = as.POSIXct(time_of_day), y = thermal_mean - 273.15, color = Elev_Aspect)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = custom_colors, name = "Elevation - Aspect", 
                     labels = c("2000m East", "2800m East","2000m West", "2800m West")) +
  xlab("Time of day") + ylab("Leaf temperature (˚C)") +
  theme_classic() +
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
  ) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours")  # Set breaks every 2 hours with specific format

tleaf.time


#plot tleaf vs. tair for both sites
tleaf.tair <- FLIRdata %>%
  ggplot(aes(x = temp_atm - 273.15, y = thermal_mean - 273.15, color = Elev_Aspect)) +
  geom_point(size = 2.5, alpha = 0.6) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_manual(values = custom_colors, name = "Elevation - Aspect", 
                     labels = c("2000m East", "2800m East", "2000m West","2800m West")) +
  xlab("FLIR air temperature (˚C)") + ylab("Leaf temperature (˚C)") +
  theme_classic() +
  coord_cartesian(xlim = c(26, 47), ylim = c(26, 47)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black", aes(group = 1)) +  # Using color for regression line
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

tleaf.tair


#Make a lil function for formating the time in the legend:
format_time <- function(x) {
  hours <- floor(x)
  minutes <- (x - hours) * 60
  sprintf("%02d:%02d", hours, round(minutes))
}

#Make the FLIR faceted plot with fading color for time:
tleaf.tair.ft <- FLIRdata %>%
  mutate(time_of_day = hour(datetime) + minute(datetime) / 60 + second(datetime) / 3600)%>%
  filter(Elev_Aspect != "2000 - E") %>%
  ggplot(aes(x = temp_atm - 273.15, y = thermal_mean - 273.15, color = time_of_day)) +
  geom_point(size = 2.5, alpha = 0.6) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_gradient(low = "blue", high = "red", name = "Time of Day", labels = format_time) +
  xlab("Air temperature (˚C)") + ylab("Leaf temperature (˚C)") +
  theme_classic() +
  coord_cartesian(xlim = c(26, 47), ylim = c(26, 47)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black", aes(group = 1)) +  # Using color for regression line
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  facet_wrap(~ Elev_Aspect)

tleaf.tair.ft


####Now the same 2 but with tomst air data
#plot tleaf vs. tair for both sites
tleaf.tairt <- FLIRtomst %>%
  ggplot(aes(x = tomst.Tair, y = flir.Tleaf, color = Elev_Aspect)) +
  geom_point(size = 2.5, alpha = 0.6) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_manual(values = custom_colors, name = "Elevation - Aspect", 
                     labels = c("2000m East", "2000m West", "2800m East", "2800m West")) +
  xlab("Tomst air temperature (˚C)") + ylab("Leaf temperature (˚C)") +
  theme_classic() +
  coord_cartesian(xlim = c(10, 47), ylim = c(10, 47)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black", aes(group = 1)) +  # Using color for regression line
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

tleaf.tairt

##Now make the faceted version:
tleaf.tair.ft.t <- FLIRtomst %>%
  mutate(time_of_dayz = hour(time_of_day) + minute(time_of_day) / 60 + second(time_of_day) / 3600)%>%
  ggplot(aes(x = tomst.Tair, y = flir.Tleaf, color = time_of_dayz)) +
  geom_point(size = 2.5, alpha = 0.6) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_gradient(low = "blue", high = "red", name = "Time of Day", labels = format_time) +
  xlab("Air temperature (˚C)") + ylab("Leaf temperature (˚C)") +
  theme_classic() +
  coord_cartesian(xlim = c(10, 47), ylim = c(10, 47)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black", aes(group = 1)) +  # Using color for regression line
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  facet_wrap(~ Elev_Aspect)

tleaf.tair.ft.t


figure3=ggarrange(tleaf.tair.ft.t+rremove("ylab"), tleaf.time+rremove("ylab"), nrow=2, ncol=1, heights=c(2,1), labels=c("A","B"))
(figure3= annotate_figure(figure3, left = textGrob("Leaf temperature (˚C)", rot = 90, vjust = 1, gp = gpar(cex = 1.3))))
#gotta manually add the p-values...
#1E: p = 0.003
#1W: p = 1.57e-21
#5E: p = 0.021
#5W: p = 1.37e-11 

#same but with tomst air instead of flir air
model.all2 <- lm(flir.Tleaf ~ tomst.Tair, data = FLIRtomst)
summary(model.all2)
slope_se <- summary(model.all2)$coefficients["tomst.Tair", "Std. Error"]
t_value <- coef(model.all2)["tomst.Tair"] / slope_se
df <- summary(model.all2)$df[2]  # Residual degrees of freedom
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)

FLIRtomst.1E <- FLIRtomst%>%
  filter(site==1, Aspect == "E")
model.all2.1E <- lm(flir.Tleaf ~ tomst.Tair, data = FLIRtomst.1E)
summary(model.all2.1E)
slope <- coef(model.all2.1E)["tomst.Tair"]
slope_se <- summary(model.all2.1E)$coefficients["tomst.Tair", "Std. Error"]
t_value <- (slope - 1) / slope_se
df <- summary(model.all2.1E)$df[2] 
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)

FLIRtomst.1W <- FLIRtomst%>%
  filter(site==1, Aspect == "W")
model.all2.1W <- lm(flir.Tleaf ~ tomst.Tair, data = FLIRtomst.1W)
summary(model.all2.1W)
slope <- coef(model.all2.1W)["tomst.Tair"]
slope_se <- summary(model.all2.1W)$coefficients["tomst.Tair", "Std. Error"]
t_value <- (slope - 1) / slope_se
df <- summary(model.all2.1W)$df[2] 
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)

FLIRtomst.5E <- FLIRtomst%>%
  filter(site==5, Aspect == "E")
model.all2.5E <- lm(flir.Tleaf ~ tomst.Tair, data = FLIRtomst.5E)
summary(model.all2.5E)
slope <- coef(model.all2.5E)["tomst.Tair"]
slope_se <- summary(model.all2.5E)$coefficients["tomst.Tair", "Std. Error"]
t_value <- (slope - 1) / slope_se
df <- summary(model.all2.5E)$df[2] 
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)

FLIRtomst.5W <- FLIRtomst%>%
  filter(site==5, Aspect == "W")
model.all2.5W <- lm(flir.Tleaf ~ tomst.Tair, data = FLIRtomst.5W)
summary(model.all2.5W)
slope <- coef(model.all2.5W)["tomst.Tair"]
slope_se <- summary(model.all2.5W)$coefficients["tomst.Tair", "Std. Error"]
t_value <- (slope - 1) / slope_se
df <- summary(model.all2.5W)$df[2] 
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)




filt.fltm = FLIRtomst%>%
  filter(site==5, Aspect=="W")
summary(lm(flir.Tleaf~tomst.Tair, data=filt.fltm))

# Linear regression of tleaf vs tair for all data - t-test to see if slope is less than 1. 
model.all <- lm(thermal_mean ~ temp_atm, data = FLIRdata)
summary(model.all)
slope_se <- summary(model.all)$coefficients["temp_atm", "Std. Error"]
t_value <- coef(model.all)["temp_atm"] / slope_se
df <- summary(model.all)$df[2]  # Residual degrees of freedom
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)

### Now for each plot individually:
FLIRdata.1W <- FLIRdata%>%
  filter(site==1, Aspect == "W")
mod.1W <- lm(thermal_mean ~ temp_atm, data = FLIRdata.1W)
summary(mod.1W)
slope <- coef(mod.1W)["temp_atm"]
slope_se <- summary(mod.1W)$coefficients["temp_atm", "Std. Error"]
t_value <- (slope - 1) / slope_se
df <- summary(mod.1W)$df[2] 
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)

FLIRdata.5W <- FLIRdata%>%
  filter(site==5, Aspect == "W")
mod.5W <- lm(thermal_mean ~ temp_atm, data = FLIRdata.5W)
summary(mod.5W)
slope <- coef(mod.5W)["temp_atm"]
slope_se <- summary(mod.5W)$coefficients["temp_atm", "Std. Error"]
t_value <- (slope - 1) / slope_se
df <- summary(mod.5W)$df[2] 
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)

FLIRdata.5E <- FLIRdata%>%
  filter(site==5, Aspect == "E")
mod.5E <- lm(thermal_mean ~ temp_atm, data = FLIRdata.5E)
summary(mod.5E)
slope <- coef(mod.5E)["temp_atm"]
slope_se <- summary(mod.5E)$coefficients["temp_atm", "Std. Error"]
t_value <- (slope - 1) / slope_se
df <- summary(mod.5E)$df[2] 
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)




################################################################
#Plot and perform regression on tomst / FLIR data for predictions
flir.tom.mod = lm(flir.Tleaf~tomst.Tair, data=FLIRtomst) #Rsq = 0.27
summary(flir.tom.mod)

flir.tom.mod1 = lm(flir.Tleaf~tomst.Tair + T1 + T2, data=FLIRtomst) 
summary(flir.tom.mod1)
vif(flir.tom.mod1)# shows that we need to get rid of T2 becuase it is over 5

flir.tom.mod2 = lm(flir.Tleaf~tomst.Tair + T1, data=FLIRtomst) #Rsq = 0.28
summary(flir.tom.mod2)
vif(flir.tom.mod2) #Good


FLIRtomst$TOD_seconds = as.numeric(FLIRtomst$TOD_seconds)
tomstFLIR$TOD_seconds = as.numeric(tomstFLIR$TOD_seconds)
FLIRtomst$site = as.character(FLIRtomst$site)
tomstFLIR$site = as.character(tomstFLIR$site)

flir.tom.mod3 = lmer(flir.Tleaf~tomst.Tair*TOD_seconds+ (1|site/Aspect), data=FLIRtomst) #Rsq = 0.53
r.squaredGLMM(flir.tom.mod3)
summary(flir.tom.mod3)
############################
flir.tom.mod5 = lmer(flir.Tleaf~tomst.Tair*TOD_seconds*Elevation + (1|Aspect), data=FLIRtomst) #Rsq = 0.52
r.squaredGLMM(flir.tom.mod5)
summary(flir.tom.mod5)
############################
############################
############################

tomstFLIR.preds <- tomstFLIR %>%
  filter(!is.na(TOD_seconds))%>%
  mutate(predicted.Tleaf = predict(flir.tom.mod5, newdata = .))


plot.fl.tm <- ggplot(tomstFLIR.preds, aes(x = tomst.Tair)) +
  geom_point(aes(y = predicted.Tleaf, color = "Predicted"), alpha = 0.3) +
  geom_point(aes(y = flir.Tleaf, color = "Measured"), alpha = 0.6) +
  geom_smooth(aes(y = predicted.Tleaf, color = "Predicted"), color = "#007766", method = 'lm', se = TRUE) +
  geom_smooth(aes(y = flir.Tleaf, color = "Measured"),color="#220066", method = 'lm', se = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  coord_cartesian(xlim = c(10, 50), ylim = c(10, 50)) +
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)) +
  labs(x = "Tomst air temperature (˚C)", y = "Leaf temperature (˚C)", color = "Tleaf Source") +
  theme_classic() +
  scale_color_manual(values = c("Measured" = "#332299", "Predicted" = "#44AA99"))

plot.fl.tm

##### See if the predicted Tleaf vs Tair is still homeothermy
model.preds <- lm(predicted.Tleaf ~ tomst.Tair, data = tomstFLIR.preds)
summary(model.preds)
slope <- coef(model.preds)["tomst.Tair"]
slope_se <- summary(model.preds)$coefficients["tomst.Tair", "Std. Error"]
t_value <- (slope - 1) / slope_se
df <- summary(model.preds)$df[2] 
p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
print(p_value)
  
################################
#plot of air temperature vs. elevation (boxplot)
(ggplot(tomstFLIR, aes(x=Elevation, y=tomst.Tair, group = Elevation, color=Elevation)) +
    geom_boxplot() + 
    theme_classic()+
    xlab("Elevation (m)")+ ylab("Air Temperature (˚C)")+labs(color = "Elevation (m)")
  )
(summary(lm(tomst.Tair~Elevation, data=tomstFLIR)))


#Air temperatures from daytime only (8:30-6:20) but every day that tomst was out (12/07/2023 – 12/15/2023)
#FIGURE 1: Air temperature vs Elevation
x_breaks <- c(2000, 2200, 2400, 2600, 2800)
airelevation.plot <- ggplot(tomstFLIR.preds, aes(x = Elevation, y = tomst.Tair)) +
  geom_point(position = position_jitter(width = 100, height = 0), size = 2, alpha = 0.2, color="#009988")+
  theme_classic() +
  xlab("Elevation (m)") + 
  ylab("Air temperature (˚C)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 17),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey")+
  scale_x_continuous(breaks = x_breaks)

model.tairelev <- lm(tomst.Tair ~ Elevation, data = tomstFLIR.preds)
p_value <- summary(model.tairelev)$coefficients[2,4]
formatted_p_value <- sprintf("p-value < 2e-16 ***")
(airelev.plot <- airelevation.plot +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))

#################################


