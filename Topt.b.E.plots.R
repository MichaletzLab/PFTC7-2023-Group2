############## Load Data: ##############
discard.weibull <- read_csv("outputs/discard.weibull.SANW.csv")%>%
  left_join(TempSummary, by="site")
discard.weibull$curveID = as.numeric(discard.weibull$curveID)
discard.weibull= discard.weibull%>%
  mutate(country = case_when(curveID > 1000~"Norway",curveID < 1000~"SAfrica"))
discard.schoolfield <- read.csv("outputs/discard.schoolfield.SANW.csv")%>%
  mutate(country = case_when(curveID > 1000~"Norway",curveID < 1000~"SAfrica"))%>%
  left_join(TempSummary, by="site")
topts.gam.add = read.csv("gam_fits_add.csv")%>%
  mutate(country = case_when(curveid > 1000~"Norway",curveid < 1000~"SAfrica"))%>%
  rename(curveID=curveid)%>%
  select(curveID, tleaf_cond_breadth)
discard.weibull <- left_join(discard.weibull, topts.gam.add, by="curveID")



############## Elevation vs. Topt - Weibull - All data ############## 
elev.asp.plot <- ggplot(discard.weibull, aes(x = Elevation, y = T_opt, color = AirTemp)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black") +
  geom_point(position = position_jitter(width = 100, height = 0), size = 2) +
  scale_color_gradient(low = "blue", high = "red", name = "Air Temperature (°C)") +
  labs(shape = "Aspect") +
  theme_classic() +
  xlab("Elevation (m)") + 
  ylab("Topt (˚C)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
# Calculate the linear regression model to extract p-value
model.toptelev <- lm(T_opt ~ Elevation, data = discard.weibull)
p_value <- summary(model.toptelev)$coefficients[2,4]
# Format the p-value
formatted_p_value <- sprintf("p-value = %.3f", p_value)
# Annotate the plot with the p-value
(elev.asp.plot <- elev.asp.plot +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))



############## Elevation vs. Topt - Weibull - facet by country ############## 
elev.asp.plot <- ggplot(discard.weibull, aes(x = Elevation, y = T_opt)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black") +
  geom_point(position = position_jitter(width = 50, height = 0), size = 2, color = "black") +
  labs(shape = "Aspect") +
  theme_classic() +
  xlab("Elevation (m)") + 
  ylab("Topt (˚C)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  facet_wrap(~ country, scales = "free")  # Facet by country

# Calculate linear regression models and extract p-values for each facet
p_values <- discard.weibull %>%
  group_by(country) %>%
  summarise(p_value = summary(lm(T_opt ~ Elevation))$coefficients[2, "Pr(>|t|)"])

# Print the plot
elev.asp.plot



############## Air temperature vs. Topt - Weibull - All data ############## 
AT.plot <- ggplot(discard.weibull, aes(x = AirTemp, y = T_opt, color = Elevation)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black") +
  geom_point(position = position_jitter(width = 2, height = 0), size = 2) +
  scale_color_gradient(low = "orange", high = "forestgreen", name = "Elevation (m)") +
  labs(shape = "Aspect") +
  theme_classic() +
  xlab("Air Temperature (˚C)") + 
  ylab("Topt (˚C)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
# Calculate the linear regression model to extract p-value
model.toptAT <- lm(T_opt ~ AirTemp, data = discard.weibull)
p_value <- summary(model.toptAT)$coefficients[2,4]
# Format the p-value
formatted_p_value <- sprintf("p-value = %.3f", p_value)
# Annotate the plot with the p-value
(AT.plot <- AT.plot +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))


############## Elevation vs. breadth - Weibull - All data ############## 
breadth.elev <- discard.weibull %>%
  ggplot(aes(x = Elevation, y = getbreadth_90, color=AirTemp)) +
  geom_point(size = 2.5, 
             position = position_jitter(width = 100, height = 0)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black") +
  scale_color_gradient(low = "blue", high = "red", name = "Air Temperature (°C)") +
  xlab("Elevation (m)") + ylab("Thermal breadth (˚C)") +
  theme_classic() +
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

model.b.elev <- lm(getbreadth_90 ~ Elevation, data = discard.weibull)
p_value <- summary(model.b.elev)$coefficients[2,4]
# Format the p-value
formatted_p_value <- sprintf("p-value = %.3f", p_value)
# Annotate the plot with the p-value
(breadth.elev <- breadth.elev +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))



############## Elevation vs. breadth - GAM - All data ############## 
breadth.elev.gam <- discard.weibull %>%
  ggplot(aes(x = Elevation, y = tleaf_cond_breadth)) +
  geom_point(size = 2.5, 
             position = position_jitter(width = 100, height = 0)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black") +
  xlab("Elevation (m)") + ylab("Thermal breadth (˚C)") +
  theme_classic() +
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

model.b.elev <- lm(tleaf_cond_breadth ~ Elevation, data = discard.weibull)
p_value <- summary(model.b.elev)$coefficients[2,4]
# Format the p-value
formatted_p_value <- sprintf("p-value = %.3f", p_value)
# Annotate the plot with the p-value
(breadth.elev.gam <- breadth.elev.gam +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))


############## Elevation vs. breadth - Weibull - facet by country ############## 
elev.b.plot <- ggplot(discard.weibull, aes(x = Elevation, y = getbreadth_90)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black") +
  geom_point(position = position_jitter(width = 50, height = 0), size = 2, color = "black") +
  labs(shape = "Aspect") +
  theme_classic() +
  xlab("Elevation (m)") + 
  ylab("Breadth (˚C)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  facet_wrap(~ country, scales = "free")  # Facet by country

# Calculate linear regression models and extract p-values for each facet
p_values <- discard.weibull %>%
  group_by(country) %>%
  summarise(p_value = summary(lm(getbreadth_90 ~ Elevation))$coefficients[2, "Pr(>|t|)"])

# Print the plot
elev.b.plot



############## Air temperature vs. breadth - Weibull - All data ############## 
breadth.AT <- discard.weibull %>%
  ggplot(aes(x = AirTemp, y = getbreadth_90, color=Elevation)) +
  geom_point(size = 2.5, 
             position = position_jitter(width = 2, height = 0)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black") +
  scale_color_gradient(low = "orange", high = "forestgreen", name = "Elevation (m)") +
  xlab("Air temperature (˚C)") + ylab("Thermal breadth (˚C)") +
  theme_classic() +
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

model.breadth.AT <- lm(getbreadth_90 ~ AirTemp, data = discard.weibull)
p_value <- summary(model.breadth.AT)$coefficients[2,4]
# Format the p-value
formatted_p_value <- sprintf("p-value = %.3f", p_value)
# Annotate the plot with the p-value
(breadth.AT <- breadth.AT +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))



############## Elevation vs. E - Schoolfield - All data ############## 
E.elev <- discard.schoolfield %>%
  ggplot(aes(x = Elevation, y = E, color=AirTemp)) +
  geom_point(size = 2.5, 
             position = position_jitter(width = 100, height = 0)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black") +
  scale_color_gradient(low = "blue", high = "red", name = "Air Temperature (°C)") +
  xlab("Elevation (m)") + ylab("Activation energy") +
  theme_classic() +
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

model.E.elev <- lm(E ~ Elevation, data = discard.schoolfield)
p_value <- summary(model.E.elev)$coefficients[2,4]
# Format the p-value
formatted_p_value <- sprintf("p-value = %.3f", p_value)
# Annotate the plot with the p-value
(E.elev <- E.elev +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))



############## Elevation vs. E - Schoolfield - Facet by country ############## 
elev.E.plot <- ggplot(discard.schoolfield, aes(x = Elevation, y = E)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black") +
  geom_point(position = position_jitter(width = 50, height = 0), size = 2, color = "black") +
  labs(shape = "Aspect") +
  theme_classic() +
  xlab("Elevation (m)") + 
  ylab("Activation energy") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  facet_wrap(~ country, scales = "free")  # Facet by country

# Calculate linear regression models and extract p-values for each facet
p_values <- discard.schoolfield %>%
  group_by(country) %>%
  summarise(p_value = summary(lm(E ~ Elevation))$coefficients[2, "Pr(>|t|)"])

# Print the plot
elev.E.plot




############## Air temperature vs. E - Schoolfield - All data ############## 
E.AT <- discard.schoolfield %>%
  ggplot(aes(x = AirTemp, y = E, color=Elevation)) +
  geom_point(size = 2.5, 
             position = position_jitter(width = 2, height = 0)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black") +
  scale_color_gradient(low = "orange", high = "forestgreen", name = "Elevation (m)") +
  xlab("Air temperature (˚C)") + ylab("Activation energy") +
  theme_classic() +
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

model.E.AT <- lm(E ~ AirTemp, data = discard.schoolfield)
p_value <- summary(model.E.AT)$coefficients[2,4]
# Format the p-value
formatted_p_value <- sprintf("p-value = %.3f", p_value)
# Annotate the plot with the p-value
(E.AT <- E.AT +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))



############## Elevation vs. E_D - Schoolfield - All data ############## 
Ed.elev <- discard.schoolfield %>%
  ggplot(aes(x = Elevation, y = E_D, color=AirTemp)) +
  geom_point(size = 2.5, 
             position = position_jitter(width = 100, height = 0)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black") +
  scale_color_gradient(low = "blue", high = "red", name = "Air Temperature (°C)") +
  xlab("Elevation (m)") + ylab("Deactivation energy") +
  theme_classic() +
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

model.Ed.elev <- lm(E_D ~ Elevation, data = discard.schoolfield)
p_value <- summary(model.Ed.elev)$coefficients[2,4]
# Format the p-value
formatted_p_value <- sprintf("p-value = %.3f", p_value)
# Annotate the plot with the p-value
(Ed.elev <- Ed.elev +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))


############## Elevation vs. E_D - Schoolfield - Facet by country ############## 
elev.Ed.plot <- ggplot(discard.schoolfield, aes(x = Elevation, y = E_D)) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black") +
  geom_point(position = position_jitter(width = 50, height = 0), size = 2, color = "black") +
  labs(shape = "Aspect") +
  theme_classic() +
  xlab("Elevation (m)") + 
  ylab("Deactivation energy") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  facet_wrap(~ country, scales = "free")  # Facet by country

# Calculate linear regression models and extract p-values for each facet
p_values <- discard.schoolfield %>%
  group_by(country) %>%
  summarise(p_value = summary(lm(E_D ~ Elevation))$coefficients[2, "Pr(>|t|)"])

# Print the plot
elev.Ed.plot



############## Air temperature vs. E_D - Schoolfield - All data ############## 
Ed.AT <- discard.schoolfield %>%
  ggplot(aes(x = AirTemp, y = E_D, color=Elevation)) +
  geom_point(size = 2.5, 
             position = position_jitter(width = 2, height = 0)) +
  geom_smooth(method = 'lm', se = TRUE, color = "black") +
  scale_color_gradient(low = "orange", high = "forestgreen", name = "Elevation (m)") +
  xlab("Elevation (m)") + ylab("Deactivation energy") +
  theme_classic() +
  theme(
    text = element_text(size = 17),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

model.Ed.AT <- lm(E_D ~ AirTemp, data = discard.schoolfield)
p_value <- summary(model.Ed.AT)$coefficients[2,4]
# Format the p-value
formatted_p_value <- sprintf("p-value = %.3f", p_value)
# Annotate the plot with the p-value
(Ed.AT <- Ed.AT +
    annotate("text", x = Inf, y = Inf, label = formatted_p_value, 
             hjust = 1.1, vjust = 25, size = 5))
