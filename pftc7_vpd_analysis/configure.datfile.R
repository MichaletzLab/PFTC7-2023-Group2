# Purpose: Load and organize data for modelling preparation
# Plots:  
# Dependencies:
# Outputs: "modelling_dat.csv"

pacman::p_load(tidyverse,data.table,dplyr,tidyverse,ggplot2,car,stringr,nls.multstart,nlstools,
               mgcv,mgcViz,data.table,pacman,lubridate,gt,cowplot,grid,gridExtra,lme4,knitr,
               gratia,kableExtra,patchwork) 

# Data import =================================
dat <- read_csv("outputs/raw.discardHooks_data.csv")

weibull.topts.disc <- read_csv("weibull.discard.hooks.SANW.csv")
weibull.topts.disc <- weibull.topts.disc %>%
  select(curveID, T_opt,getbreadth_90)%>%
  rename(breadth_95 = getbreadth_90)

school.topts.disc <- read_csv("schoolfield.discard.hooks.SANW.csv")
school.topts.disc <- school.topts.disc %>%
  select(curveID, T_opt)


Faster_Key = read.csv("data/Faster_Key.csv")%>%
  rename(curveID = Obs)
norway_key = read.csv("data/Norway.Key.csv")
dat <- left_join(dat, Faster_Key, by = "curveID")
dat<-left_join(dat,norway_key,by="curveID")
dat <- left_join(dat, weibull.topts.disc, by = "curveID")
dat <- left_join(dat, school.topts.disc, by = "curveID")
dat <- dat %>%
  select(Date.y,Date.measured.x, curveID,site.x,site.y, Elevation.masl,Elevation.y, T_opt.x,breadth_95,T_opt.y,Species.y,taxon, Individual.y, BarcodeLeaf.y, hhmmss, Tleaf, A, gsw, Ci, Ca, Emm, VPDleaf, Tair, CO2_r, CO2_s)%>%
  mutate(Elevation = ifelse(is.na(Elevation.y), Elevation.masl, Elevation.y)) %>%
  mutate(site.y = site.y + 5)%>%
  mutate(site = ifelse(is.na(site.x), site.y, site.x))%>%
  mutate(Species = ifelse(is.na(Species.y), taxon, Species.y))%>%
  mutate(Date = ifelse(is.na(Date.y), Date.measured.x,Date.y))%>%
  select(Date, curveID,site, Elevation, T_opt.x,breadth_95,T_opt.y,Species, Individual.y, BarcodeLeaf.y, hhmmss, Tleaf, A, gsw, Ci, Ca, Emm, VPDleaf, Tair, CO2_r, CO2_s)


dat <- as.data.table(dat)
dat[is.na(Species), Species := "hypoxis_costata"]
# Rename columns in Dat to remain consistent with program script
names(dat) <- c(
  "date", "curveID","site", "elevation","weib.Topt","weib.breadth","School.Topt", "species", "rep", "leaf", "time", "tleaf", "photo",
  "cond", "ci", "ca", "trmmol", "vpdl", "tair", "co2r", "co2s"
)

# add this because it's important for filtering
dat$cica <- dat$ci / dat$ca

names(dat) <- tolower(names(dat))
dat[,`:=`(species = str_replace(species,"[:space:]"," "))]
dat <- dat %>% 
  filter(!is.na(vpdl)) %>% 
  filter(cica %between% c(0.05,1.2))%>% 
  filter(abs(cond) < 0.8) 

dat <- dat %>% mutate(fspecies = factor(species),
                      fsite = factor(site))

vec_tleaf <- range(dat$tleaf,na.rm=T)
vec_vpdl <- range(dat$vpdl,na.rm=T)
vec_cond <- range(dat$cond,na.rm=T)
dat[, `:=`(
  log_vpdl = log(vpdl), 
  log_cond = log(cond), 
  tleaf_c = scale(tleaf, scale = FALSE)[, 1], 
  vpdl_c = scale(vpdl, scale = FALSE)[, 1],
  cond_c = scale(cond, scale = FALSE)[, 1]
)]

dat[,`:=`(mod_group = paste0(site))]
dat[, species := str_trim(species)]
dat[, species := tolower(species)]

ref_dat <- dat[,.(photo_u = mean(photo,na.rm=T), 
                  photo_sd = sd(photo)),by=species]
dat <- merge(dat, ref_dat, by = 'species', all.x = TRUE)
dat[,`:=`(photo_rel = (photo - photo_u)/photo_u, 
          photo_z = (photo - photo_u)/photo_sd)]

vec_mod_group <- unique(dat$mod_group)

species_summary <- dat %>%
  group_by(mod_group) %>%
  summarize(unique_species = n_distinct(fspecies))
dat <- dat %>% mutate(mg = factor(paste(species,rep,leaf)))%>%
  mutate(country = case_when(curveid > 1000~"Norway",curveid < 1000~"SAfrica"))

write.csv(dat, "modelling_dat.csv")

###############
###############
## Theme ============
fn_theme <- function (base_size = 12, base_family = "", base_line_size = base_size/22, 
                      base_rect_size = base_size/24){
  half_line <- base_size/2
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.text = element_text(colour = "black", size = rel(1)),
          axis.text.x = element_text(margin = margin(5,0,0,0)), 
          axis.text.y = element_text(margin = margin(0,5,0,0)),
          axis.ticks = element_line(colour = "black", 
                                    linewidth = rel(0.5)),
          axis.ticks.length = unit(-0.1,'cm'),
          panel.border = element_rect(fill = NA, colour = "black", 
                                      linewidth = rel(1)), panel.grid = element_line(colour = "black"), 
          panel.grid.major = element_blank(), #element_line(linewidth = rel(0.1)), 
          panel.grid.minor = element_blank(), #element_line(linewidth = rel(0.05)), 
          strip.background = element_rect(fill = "transparent",
                                          color = 'transparent'), 
          strip.text = element_text(colour = "black", size = rel(0.9), 
                                    margin = margin(0.8 * half_line, 0.8 * half_line, 
                                                    0.8 * half_line, 0.8 * half_line)), 
          legend.position = c(1,1),
          legend.justification = c(1,1),
          # legend.direction = 'vertical',
          legend.background =  element_rect(fill='transparent', 
                                            color='transparent'),
          legend.text = element_text(size = rel(0.9),lineheight = 1),
          legend.spacing.y = unit(0.01, 'cm'),
          complete = TRUE)
}
# schoolfield.fit = read.csv("outputs/discard.schoolfield.SANW.csv")%>%
#   rename(curveid = curveID)%>%
#   select(curveid, T_opt,J_ref, E, E_D,breadth)%>%
#   rename(school_breadth = breadth)%>%
#   rename(T_opt.s = T_opt)
# weibull.fit = read.csv("outputs/discard.weibull.SANW.csv")%>%
#   rename(curveid = curveID)%>%
#   rename(T_opt.w = T_opt)
# # Merge dat1 with schoolfield.fits based on curveid
# dats <- dat %>%
#   left_join(schoolfield.fit, by = "curveid")%>%
#   left_join(weibull.fit, by = "curveid")
# 
# # Calculate the predicted values using the Schoolfield model
# dats <- dats %>%
#   mutate(predicted = schoolfield(tleaf, J_ref, E, E_D, T_opt.s))