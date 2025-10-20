school = read.csv("schoolfield.discard.hooks.SANW.csv")
weibull = read.csv("weibull.discard.hooks.SANW.csv")

school <- school %>%
  rename(Species.y=taxon)%>%
  rename(Species.x=Species)%>%
  rename(Individual.y = turf_number)%>%
  rename(Individual.x = Individual)%>%
  rename(T_opt_school = T_opt)%>%
  rename(T_opt_SE_school = T_opt_SE)%>%
  rename(AIC_school = AIC)%>%
  rename(r_sq_school = r_sq)

weibull <- weibull %>%
  rename(Species.y=taxon)%>%
  rename(Species.x=Species)%>%
  rename(Individual.y = turf_number)%>%
  rename(Individual.x = Individual)%>%
  rename(T_opt_weib = T_opt)%>%
  rename(T_opt_SE_weib = T_opt_SE)%>%
  rename(AIC_weib = AIC)%>%
  rename(r_sq_weib = r_sq)

# clean up data sets
fix_joined_cols <- function(df) {
  base_names <- unique(gsub("\\.(x|y)$", "", names(df)))
  out <- lapply(base_names, function(b) {
    xcol <- paste0(b, ".x")
    ycol <- paste0(b, ".y")
    if (xcol %in% names(df) & ycol %in% names(df)) {
      coalesce(df[[xcol]], df[[ycol]])
    } else if (xcol %in% names(df)) {
      df[[xcol]]
    } else if (ycol %in% names(df)) {
      df[[ycol]]
    } else {
      df[[b]]
    }})
  
  out <- as.data.frame(out)
  names(out) <- base_names
  out}

#Run function to clean data sets
schoolfield_params <- fix_joined_cols(school)
weibull_params <- fix_joined_cols(weibull)

#Join data sets and add a PC1 column corresponding to each site
parameter_dat <- left_join(schoolfield_params, weibull_params, by=c("curveID","Elevation", "site","Species","Individual", "Aspect","BarcodeLeaf","BarcodeCutout"))
PC1_df <- raw.env.data_pca %>%
  group_by(curveID) %>%
  summarise(PC1 = mean(PC1, na.rm = TRUE))

ThermTraits.dat <- left_join(parameter_dat, PC1_df, by="curveID")
#Run a gam with y=thermal traits, x=PC1
summary(mod.T_opt_sch <- gam(T_opt_school ~ s(PC1, k=3) + ###
                       Species,
                      data = ThermTraits.dat,
                      method="REML"))
summary(mod.T_opt_weib <- gam(T_opt_weib ~ s(PC1, k=3) +
                               Species,
                             data = ThermTraits.dat,
                             method="REML"))
summary(mod.Ea <- gam(E ~ s(PC1, k=3) +  ###
                               Species,
                             data = ThermTraits.dat,
                             method="REML"))
summary(mod.Ed <- gam(E_D ~ s(PC1, k=3) + ###
                        Species,
                      data = ThermTraits.dat,
                      method="REML"))
summary(mod.breadth <- gam(breadth_95 ~ s(PC1, k=3) +
                        Species,
                      data = ThermTraits.dat,
                      method="REML"))
summary(mod.getbreadth <- gam(getbreadth_90 ~ s(PC1, k=3) + ###
                       Species,
                     data = ThermTraits.dat,
                     method="REML"))
##None are significant! Interesting... Small problem though because this loses 
  #the error terms that I have stored (SE). I think this will overstate our 
  #confidence though since the error is understated which means that since  
  #our results are insignificant this is OK.

#Visualize:
School_Topt_Plot <- ggplot(ThermTraits.dat, aes(x = PC1, y = T_opt_school, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = expression(T[opt]~"(°C)"),
       color = "Species",
       title="A")
Ea_Plot <- ggplot(ThermTraits.dat, aes(x = PC1, y = E, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = expression(E[a]~"(eV)"),
       color = "Species")
Ed_Plot <- ggplot(ThermTraits.dat, aes(x = PC1, y = E_D, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = expression(E[d]~"(eV)"),
       color = "Species")
breadth_Plot <- ggplot(ThermTraits.dat, aes(x = PC1, y = breadth_95, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = "breadth (°C)",
       color = "Species")

A.Therm.plot <- ggarrange(School_Topt_Plot, breadth_Plot, Ea_Plot, Ed_Plot, nrow=2, ncol=2, common.legend = TRUE, labels = c("A","B","C","D"),legend="right")
