# Purpose: Check correlations
# Plots:  plots of topt, breadth, Ea and Ed vs. air temperature as well as supplementary figures
# Dependencies:
# Outputs:

## Plot VPD vs. Tleaf==================
model.vpd.t <- gam(vpdl ~ s(tleaf, bs = 'ts', k = 5), data = dat)

# Extract the R-squared value
rsq <- summary(model.vpd.t)$r.sq
formatted_rsq1 <- sprintf("r² = %.3f", rsq)

# Create the plot
(p_out <- dat %>% 
  ggplot(aes(tleaf, vpdl)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = 'gam',
              formula = y ~ s(x, bs = 'ts', k = 5),
              color = "#cf0000") + 
  labs(x = "Leaf temperature (°C)",
       y = "VPD (kPa)") + 
  theme_classic() +
  annotate("text", x = 20, y = Inf, label = formatted_rsq1, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "italic"))

## Plot A vs gsw==================
(a_out <- dat %>% 
  ggplot(aes(cond,photo))+
  geom_point(alpha=0.1)+
  geom_smooth(method='gam',
              formula = y~s(x,bs='ts',k=5),
              color = "#cf0000") + 
  labs(x = "Conductance",
       y = "Photosynthesis") + 
  theme_classic()+
  annotate("text", x = 20, y = Inf, label = formatted_rsq1, 
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "italic"))


## Plot tleaf vs gsw==================
model.gsw.t <- gam(cond ~ s(tleaf, bs = 'ts', k = 5), data = dat)

# Extract the R-squared value
rsq <- summary(model.gsw.t)$r.sq
formatted_rsq1 <- sprintf("r² = %.3f", rsq)

# Create the plot
(p_out <- dat %>% 
    ggplot(aes(tleaf, cond)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method="lm", se=TRUE) + 
    coord_cartesian(ylim = c(0, 0.56)) +
    labs(x = "Leaf temperature (°C)",
         y = expression(paste("Conductance (µmol ", CO[2], " ", m^-2, s^-1, ")"))) + 
    theme_classic() +
    annotate("text", x = 20, y = Inf, label = formatted_rsq1, 
             hjust = 1.1, vjust = 1.1, size = 5, fontface = "italic"))


## Plot m vs. max gs.==================
# Calc reference gs ===================
ref_gs <- dat %>% 
  filter(vpdl <= 1.4) %>%
  # # filter(vpdl %between% c(0.8,1.2)) %>%
  group_by(species) %>% 
  summarize(gs_q90 = quantile(cond[vpdl <= 1.4],0.9, na.rm=T),
            gs_max = max(cond[vpdl <= 1.4], na.rm=T),
            gs_sd = sd(cond,na.rm=T),
            nobs = n()) %>% 
  ungroup() %>% 
  mutate(gsref = ifelse(is.na(gs_q90)==T, yes = gs_max, no = gs_q90)) %>% 
  as.data.table()


gs_dat <- merge(dat,ref_gs,by='species')

fn_m <- function(tmp_dat){
  gam(I(cond-gsref) ~ 0 + 
        log(vpdl), 
      data=tmp_dat,
      method='REML') %>% 
    broom::tidy(parametric = T)
  
}
gs_dat[,(fn_m(.SD)), by=.(site)] %>% 
  fwrite(., 
         file = "pftc7_vpd_analysis/data.files/fit_m_gs-vpd-sens_site.csv")
gs_dat[,(fn_m(.SD)), by=.(species)] %>% 
  fwrite(., 
         file = "pftc7_vpd_analysis/data.files/fit_m_gs-vpd-sens_species.csv")
gs_dat[,(fn_m(.SD)), by=.(elevation,species)] %>% 
  fwrite(., 
         file = "pftc7_vpd_analysis/data.files/fit_m_gs-vpd-sens_elevation_species.csv")
gs_dat[,(fn_m(.SD)), by=.(site,species)] %>% 
  fwrite(., 
         file = "pftc7_vpd_analysis/data.files/fit_m_gs-vpd-sens_site_species.csv")


fits_species <- gs_dat[,(fn_m(.SD)), by=.(species)]
fits_species %>% 
  select(-term) %>% 
  rename(m = estimate) %>% 
  mutate(m = -m) %>% 
  # kableExtra::kable(format = "pipe", digits = c(1,2,3,3,3)) %>% 
  gt::gt() %>% 
  gt::fmt_number(columns = m,
                 decimals = 3) %>% 
  gt::fmt_number(columns = c(std.error,statistic,p.value),
                 decimals = 3) %>% 
  gt::gtsave(., filename="pftc7_vpd_analysis/data.files/fit_m_gs-vpd-sens_species.docx")

# Plot m ~ gsref =====================

merge(fits_species,ref_gs,by='species') %>% 
  gam(I(-estimate)~gsref,data=.) %>% 
  summary()
merge(fits_species,ref_gs,by='species') %>% 
  ggplot(aes(gsref,-estimate))+
  geom_point(pch=21)+
  geom_smooth(method='lm', 
              formula = y~x,
              color = 'black',
              alpha =0.5) + 
  labs(x = expression(paste(g[s][" ref."])),
       y = expression(paste(italic(m))))+
  fn_theme()
ggsave(filename = "pftc7_vpd_analysis/figures/gs_sensitivity_m_gsref.png",
       width = 12,
       height = 10,
       units='cm',
       dpi = 600)
merged_data <- merge(fits_species, ref_gs, by = 'species')
summary(lm(I(-estimate) ~ gsref, data = merged_data))

########################
site.spec.m.dat = read.csv("pftc7_vpd_analysis/data.files/fit_m_gs-vpd-sens_elevation_species.csv")
ggplot(data = site.spec.m.dat, aes(x = elevation, y = estimate)) +
  geom_point(aes(color = species)) + # Keep species colors for points
  geom_smooth(method = "loess", se = TRUE, linetype = "solid", size = 1, color = "black") + # Add single curved line
  theme_classic() +
  ylab("m") + 
  xlab("Elevation (m)")
summary(lm(estimate~elevation+species, data=site.spec.m.dat))

########################
#Check OLS regression with just Tleaf, then just VPD, then both to ensure concurvity/multicollinearity is not an issue
########################
summary(lm(photo~species+elevation+tleaf,data=dat)) #-0.081 ----> 0.23
summary(lm(photo~species+elevation+cond,data=dat)) #28.49
summary(lm(photo~species+elevation+cond+tleaf,data=dat))
summary(lm(photo~species+elevation+vpdl,data=dat)) #-0.93 -----> -2.79
summary(lm(photo~species+elevation+vpdl+tleaf,data=dat))
vif(modmulti)

summary(glm(photo~species+elevation+poly(tleaf,2),data=dat,family = gaussian(link = "identity"))) 
summary(glm(photo~species+elevation+poly(vpdl,2),data=dat,family = gaussian(link = "identity")))
summary(lm(photo~species+elevation+poly(vpdl,2)+poly(tleaf,2),data=dat,family = gaussian(link = "identity")))

########################
