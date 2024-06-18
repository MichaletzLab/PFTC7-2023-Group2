# Libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(nls.multstart)
library(nlstools)
library(data.table)
library(mgcv)
library(mgcViz)
library(gratia)
library(patchwork)
library(pacman)
library(data.table)
library(lubridate)

pacman::p_load(tidyverse,data.table,
               mgcv,mgcViz,
               gratia,
               patchwork) 

# Data import =================================
dat <- read_csv("outputs/raw.discardHooks_data.csv")

weibull.topts.disc <- read_csv("outputs/discard.weibull.SANW.csv")
weibull.topts.disc <- weibull.topts.disc %>%
  select(curveID, T_opt)

#Data for hooks ==============================
weibull.cuthooks.raw <- read_csv("outputs/rawData.at.subset2.SANW.csv")
weibull.cuthooks <- weibull.cuthooks.raw %>%
  select(curveID)
weibull.cuthooks <- distinct(weibull.cuthooks)

onlyhooks <- anti_join(weibull.cuthooks, weibull.topts.disc, by = "curveID") %>%
  select(curveID)

raw.hooks <- semi_join(weibull.cuthooks.raw, onlyhooks, by = "curveID")

#dat = raw.hooks
#=============================================

Faster_Key = read.csv("data/Faster_Key.csv")%>%
  rename(curveID = Obs)
norway_key = read.csv("data/Norway.Key.csv")
dat <- left_join(dat, Faster_Key, by = "curveID")
dat<-left_join(dat,norway_key,by="curveID")
dat <- left_join(dat, weibull.topts.disc, by = "curveID")
dat <- dat %>%
  select(Date.y,Date.measured.x, curveID,site.x,site.y, Elevation.masl,Elevation.y, T_opt,Species.y,taxon, Individual.y, BarcodeLeaf.y, hhmmss, Tleaf, A, gsw, Ci, Ca, Emm, VPDleaf, Tair, CO2_r, CO2_s)%>%
  mutate(Elevation = ifelse(is.na(Elevation.y), Elevation.masl, Elevation.y)) %>%
  mutate(site.y = site.y + 5)%>%
  mutate(site = ifelse(is.na(site.x), site.y, site.x))%>%
  mutate(Species = ifelse(is.na(Species.y), taxon, Species.y))%>%
  mutate(Date = ifelse(is.na(Date.y), Date.measured.x,Date.y))%>%
  select(Date, curveID,site, Elevation, T_opt,Species, Individual.y, BarcodeLeaf.y, hhmmss, Tleaf, A, gsw, Ci, Ca, Emm, VPDleaf, Tair, CO2_r, CO2_s)


dat <- as.data.table(dat)
dat[is.na(Species), Species := "hypoxis_costata"]
# Rename columns in Dat to remain consistent with program script
names(dat) <- c(
  "date", "curveID","site", "elevation","weib.Topt", "species", "rep", "leaf", "time", "tleaf", "photo",
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
dat <- dat %>% mutate(mg = factor(paste(species,rep,leaf)))
