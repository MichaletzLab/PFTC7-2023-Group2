school.iWUE = read.csv("schoolfield.iWUE.discard.hooks.SANW.csv")

school.iWUE <- school.iWUE %>%
  rename(Species.y=taxon)%>%
  rename(Species.x=Species)%>%
  rename(Individual.y = turf_number)%>%
  rename(Individual.x = Individual)%>%
  rename(T_opt_school = T_opt)%>%
  rename(T_opt_SE_school = T_opt_SE)%>%
  rename(AIC_school = AIC)%>%
  rename(r_sq_school = r_sq)

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
schoolfield_params.i <- fix_joined_cols(school.iWUE)

#Join data sets and add a PC1 column corresponding to each site
parameter_dat.i <- schoolfield_params.i
PC1_df <- raw.env.data_pca %>%
  group_by(curveID) %>%
  summarise(PC1 = mean(PC1, na.rm = TRUE))

ThermTraits.dat.i <- left_join(parameter_dat.i, PC1_df, by="curveID")
#Run a gam with y=thermal traits, x=PC1
summary(mod.T_opt_sch <- gam(T_opt_school ~ s(PC1, k=3) + ###
                               Species,
                             data = ThermTraits.dat.i,
                             method="REML"))
summary(mod.Ea <- gam(E ~ s(PC1, k=3) +  ###
                        Species,
                      data = ThermTraits.dat.i,
                      method="REML"))
summary(mod.Ed <- gam(E_D ~ s(PC1, k=3) + ###
                        Species,
                      data = ThermTraits.dat.i,
                      method="REML"))
summary(mod.breadth <- gam(breadth_95 ~ s(PC1, k=3) +
                             Species,
                           data = ThermTraits.dat.i,
                           method="REML"))
##None are significant! Interestini... Small problem though because this loses 
#the error terms that I have stored (SE). I think this will overstate our 
#confidence though since the error is understated which means that since  
#our results are insignificant this is OK.

#Visualize:
i.School_Topt_Plot <- ggplot(ThermTraits.dat.i, aes(x = PC1, y = T_opt_school, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = expression(T[opt]~"(°C)"),
       color = "Species",
       title="iWUE (A/gsw)")
i.Ea_Plot <- ggplot(ThermTraits.dat.i, aes(x = PC1, y = E, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = expression(E[a]~"(eV)"),
       color = "Species")
i.Ed_Plot <- ggplot(ThermTraits.dat.i, aes(x = PC1, y = E_D, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = expression(E[d]~"(eV)"),
       color = "Species")
i.breadth_Plot <- ggplot(ThermTraits.dat.i, aes(x = PC1, y = breadth_95, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = "breadth (°C)",
       color = "Species")

iWUE.Therm.plot <- ggarrange(i.School_Topt_Plot, i.breadth_Plot, i.Ea_Plot, i.Ed_Plot, nrow=2, ncol=2, common.legend = TRUE, labels = c("A","B","C","D"),legend="right")
