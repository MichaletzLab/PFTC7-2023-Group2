# Objective: Fit GAM (General Additive Models)
# to LI-COR measurements from PFTC7

# Date: 2024.06.06

# ====================================================
# The following code for fitting GAM models
# was adapted from code written by Martijn Slot.
# Source:
#   Slot M, Winter K. 2017.
#   New Phyt doi: 10.1111/nph.14469
#
# ====================================================

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

pacman::p_load(tidyverse,data.table,mgcv,mgcViz,
               gratia,
               patchwork, 
               nls.multstart)

# Data import =================================
dat <- read_delim("./raw.discardHooks_data.csv", delim=",") %>% 
  setDT()

weibull.topts.disc <- read.csv("discard.weibull.SANW.csv")
weibull.topts.disc <- weibull.topts.disc %>%
  select(curveID, T_opt)

#Data for hooks ==============================
weibull.cuthooks <- read.csv("weibull.cut.hooks.csv")
weibull.cuthooks <- weibull.cuthooks %>%
  select(curveID, T_opt)

onlyhooks <- anti_join(weibull.cuthooks, weibull.topts.disc, by = "curveID") %>%
  select(curveID)

hooks <- read.csv("rawData.at.subset2.csv")
raw.hooks <- semi_join(hooks, onlyhooks, by = "curveID")

dat = hooks
#=============================================
Faster_Key = read.csv("Faster_Key.csv")%>%
  rename(curveID = Obs)
norway_key = read.csv("trait.data.with.area.csv")
dat <- left_join(dat, Faster_Key, by = "curveID")
dat<-left_join(dat,norway_key,by="curveID")
dat <- left_join(dat, weibull.topts.disc, by = "curveID")
dat <- dat %>%
  select(Date.y, curveID,site.x,site.y, Elevation.masl,Elevation.y, T_opt,Species.y,taxon, Individual.y, BarcodeLeaf.y, hhmmss, Tleaf, A, gsw, Ci, Ca, Emm, VPDleaf, Tair, CO2_r, CO2_s)%>%
  mutate(Elevation = ifelse(is.na(Elevation.y), Elevation.masl, Elevation.y)) %>%
  mutate(site.y = site.y + 5)%>%
  mutate(site = ifelse(is.na(site.x), site.y, site.x))%>%
  mutate(Species = ifelse(is.na(Species.y), taxon, Species.y))%>%
  select(Date.y, curveID,site, Elevation, T_opt,Species, Individual.y, BarcodeLeaf.y, hhmmss, Tleaf, A, gsw, Ci, Ca, Emm, VPDleaf, Tair, CO2_r, CO2_s)


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
                      #fform = factor(form))

vec_tleaf <- range(dat$tleaf,na.rm=T)
vec_vpdl <- range(dat$vpdl,na.rm=T)

dat$vpdl = as.numeric(dat$vpdl)
dat$tleaf = as.numeric(dat$tleaf)
setDT(dat)

# Create new columns
dat[, `:=`(
  log_vpdl = log(vpdl), 
  tleaf_c = scale(tleaf, scale = FALSE)[, 1], 
  vpdl_c = scale(vpdl, scale = FALSE)[, 1]
)]

dat[,`:=`(mod_group = paste0(site))]

ref_dat <- dat[,.(photo_u = mean(photo,na.rm=T), 
                  photo_sd = sd(photo)),by=species]
dat <- merge(dat,ref_dat,by='species')
dat[,`:=`(photo_rel = (photo - photo_u)/photo_u, 
          photo_z = (photo - photo_u)/photo_sd)]

# fdoy= factor(doy)
dat <- dat %>% mutate(mg = factor(paste(species,rep,leaf)))

# # Fit a 2D kernel density ===============================
# mat <- cbind(dat$tleaf,dat$vpdl)
# bins <- KernSmooth::bkde2D(mat, bandwidth = c(0.25,0.25), gridsize = c(50L, 50L))
# tmp <- cbind(data.table(tleaf=rep(bins$x1,length(bins$x2)),
#                         vpdl=sort(rep(bins$x2,length(bins$x1)))),
#              reshape2::melt(bins$fhat,value.name = 'p'))
# 
# # Fit a 2D GAM to the kernel density smooth =============
# pm <- gam(p ~ te(tleaf,vpdl),
#           family = betar,
#           data=tmp)
# pdat <- data.table(vpdl = seq(0.5,7.23,length.out=300)) %>% 
#   expand_grid(., tleaf = seq(20,50,length.out=30)) %>% 
#   mutate(log_vpdl = log(vpdl)) %>% 
#   mutate(p = predict(pm, type='response',newdata=.)) %>% 
#   filter(p >= 0.005) %>% # can change 
#   mutate(fsite = "PNM") %>% 
#   mutate(fspecies = levels(dat$fspecies)[20]) %>% 
#   setDT()
# pdat[,`:=`(log_vpdl = log(vpdl), 
#            tleaf_c = scale(tleaf,scale=F)[,1],
#            vpdl_c = scale(vpdl,scale=F)[,1])]



# Sub models ===============================================
vec_mod_group <- unique(dat$mod_group)


vv <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group){
    ## only factors and random effects that notably improved AIC
    tmp_dat <- dat[mod_group ==sel_group]
    o <- gam(photo ~ 
               s(fspecies,bs='re') + 
               s(vpdl,bs='ts',k=5) + 
               s(tleaf,bs='ts',k=5), 
             data=tmp_dat,
             select=T,
             method = "REML")
    o$mod_group <- sel_group
    # o$form <- unique(tmp_dat$form)
    o$site <- unique(tmp_dat$site)
    return(o)
  })

vv2 <- vv %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(tleaf)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')
vv2$R2 %>% summary


vv3 <- vv %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(vpdl)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')
vv3$R2 %>% summary


ww <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group){
    tmp_dat <- dat[mod_group ==sel_group]
    o <- gam(photo ~ 
               s(fspecies,bs='re') + 
               s(tleaf,bs='ts',k=5), 
             data=tmp_dat,
             select=T,
             method = "REML")
    o$mod_group <- sel_group
    o$form <- unique(tmp_dat$form)
    o$site <- unique(tmp_dat$site)
    return(o)
  })

ww2 <- ww %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(tleaf)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = "m_t")
ww2$R2 %>% summary

zz <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group){
    tmp_dat <- dat[mod_group ==sel_group]
    o <- gam(photo ~ 
               s(fspecies,bs='re') + 
               s(vpdl,bs='ts',k=5),
             data=tmp_dat,
             select=T,
             method = "REML")
    o$mod_group <- sel_group
    o$form <- unique(tmp_dat$form)
    o$site <- unique(tmp_dat$site)
    return(o)
  })

zz2 <- zz %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(vpdl)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_v')
zz2$R2 %>% summary


# Extract Topt ==================
fig_dat_topt <- 
  rbindlist(list(vv2,ww2),fill=T) %>% 
  group_by(mod,site) %>% 
  filter(.estimate == max(.estimate)) %>% 
  filter(.estimate < 41) %>% 
  ungroup() %>% 
  select(mod,site,tleaf) %>% 
  rename(Topt = tleaf)

avg_weib_topt <- dat %>%
  group_by(site) %>%
  summarize(avg_weib_topt = mean(weib.topt, na.rm = TRUE)) %>%
  ungroup()

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


## Temperature Response Figure =============
p_photo <- rbindlist(list(vv2,ww2),fill=T) %>% 
  ggplot(aes(tleaf, .estimate,color=mod,fill=mod))+
  geom_hline(aes(yintercept = 0), 
             col = 'grey40',
             lwd=0.5,
             lty=1) + 
  geom_vline(data=fig_dat_topt,
             show.legend = F,
             aes(xintercept = Topt,
                 color = mod),
             lty = 2) +
  # geom_hline(aes(yintercept=0),
  #            col='grey',
  #            lty=2) + 
  geom_ribbon(aes(tleaf,
                  ymin=.estimate-1.96*.se,
                  ymax=.estimate+1.96*.se, 
                  fill = mod),
              color=NA,
              alpha=0.25)+
  geom_line(lwd=1)+
  scale_color_manual(
    values = c("m_t" = "#cf0000",
               "m_tv" = "navy"),
    labels = c("m_t" = "Temperature",
               "m_tv" = "Temperature + VPD")
  ) +
  scale_fill_manual(
    values = c("m_t" = "#cf0000",
               "m_tv" = "navy"),
    labels = c("m_t" = "Temperature",
               "m_tv" = "Temperature + VPD")
  ) +
  scale_x_continuous(breaks = seq(5, 40, 5), 
                     expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(color='GAM',
       fill='GAM', 
       y = expression(paste("(Partial) effect on Photosynthesis (µmol ",
                            CO[2]," ",
                            m**-2, s**-1,")")),
       # x = expression(paste(T[leaf]~"(°C)"))
       x = "Leaf temperature (°C)") +
  facet_wrap(~factor(site,levels=c("6","7","8","9","1", "2", "3", "4", "5"),
                     labels=c("Norway 1","Norway 2", "Norway 3","Norway 4","S.Africa 1", "S.Africa 2", "S.Africa 3", "S.Africa 4", "S.Africa 5")), 
             scales = 'free_y',
             ncol = 1,
             shrink = T) +
  fn_theme()+
  geom_vline(data = avg_weib_topt, show.legend = FALSE, aes(xintercept = avg_weib_topt), color = 'grey', lty = 2) +
  theme(legend.position = 'bottom', 
        strip.text = element_text(face='bold',
                                  hjust=0)); p_photo



## Conductance models and figures ------------------------------------------
vv <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group){
    ## only factors and random effects that notably improved AIC
    tmp_dat <- dat[mod_group ==sel_group]
    o <- gam(cond ~ 
               s(fspecies,bs='re') + 
               s(vpdl,bs='ts',k=5) + 
               s(tleaf,bs='ts',k=5), 
             data=tmp_dat,
             select=T,
             method = "REML")
    o$mod_group <- sel_group
    # o$form <- unique(tmp_dat$form)
    o$site <- unique(tmp_dat$site)
    return(o)
  })

vv2 <- vv %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(tleaf)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')
vv2$R2 %>% summary


vv3 <- vv %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(vpdl)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')
vv3$R2 %>% summary


ww <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group){
    tmp_dat <- dat[mod_group ==sel_group]
    o <- gam(cond ~ 
               s(fspecies,bs='re') + 
               s(tleaf,bs='ts',k=5), 
             data=tmp_dat,
             select=T,
             method = "REML")
    o$mod_group <- sel_group
    o$form <- unique(tmp_dat$form)
    o$site <- unique(tmp_dat$site)
    return(o)
  })

ww2 <- ww %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(tleaf)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = "m_t")
ww2$R2 %>% summary


zz <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group){
    tmp_dat <- dat[mod_group ==sel_group]
    o <- gam(cond ~ 
               s(fspecies,bs='re') + 
               s(vpdl,bs='ts',k=5),
             data=tmp_dat,
             select=T,
             method = "REML")
    o$mod_group <- sel_group
    o$form <- unique(tmp_dat$form)
    o$site <- unique(tmp_dat$site)
    return(o)
  })

zz2 <- zz %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(vpdl)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_v')
zz2$R2 %>% summary


## Make Figure =============
p_cond <- rbindlist(list(vv2,ww2),fill=T) %>% 
  ggplot(aes(tleaf, .estimate,color=mod,fill=mod))+
  geom_hline(aes(yintercept = 0), 
             col = 'grey40',
             lwd=0.5,
             lty=1) + 
  # geom_vline(data=fig_dat_topt,
  #            show.legend = F,
  #            aes(xintercept = Topt,
  #                color = mod),
  #            lty = 2) +
  # geom_hline(aes(yintercept=0),
  #            col='grey',
  #            lty=2) + 
  geom_ribbon(aes(tleaf,
                  ymin=.estimate-1.96*.se,
                  ymax=.estimate+1.96*.se, 
                  fill = mod),
              color=NA,
              alpha=0.25)+
  geom_line(lwd=1)+
  scale_color_manual(
    values = c("m_t" = "#cf0000",
               "m_tv" = "navy"),
    labels = c("m_t" = "Temperature",
               "m_tv" = "Temperature + VPD")
  ) +
  scale_fill_manual(
    values = c("m_t" = "#cf0000",
               "m_tv" = "navy"),
    labels = c("m_t" = "Temperature",
               "m_tv" = "Temperature + VPD")
  ) +
  scale_x_continuous(breaks = seq(5, 40, 5), 
                     expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), 
                     position = 'right') + 
  labs(color='GAM',
       fill='GAM', 
       y = expression(paste("(Partial) effect on Conductance (mol ",
                            # CO[2]," ",
                            m**-2, s**-1,")")),
       # x = expression(paste(T[leaf]~"(°C)"))
       x = "Leaf temperature (°C)") +
  facet_wrap(~factor(as.factor(site), levels=c("6","7","8","9","1", "2", "3", "4", "5"),
                     labels=c("469 m","700","920","1290","2000 m", "2200 m", "2400 m", "2600 m", "2800 m")
                     # labels=c("San Lorenzo","Parque Metropolitano")
  ), 
  scales = 'free_y',
  ncol = 1,
  shrink = T) +
  fn_theme() +
  theme(legend.position = 'bottom'); p_cond


p_out <- (p_photo|p_cond) + 
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')&
  theme(legend.position = 'bottom')

#Make sure to call it figures/HOOKSvpd_and_tleaf","_temperature-effect_ 
ggsave(p_out, 
       filename=
         paste0("figures/vpd_and_tleaf","_temperature-effect_",
                Sys.Date(),
                ".png"),
       device = grDevices::png,
       width=12*2,
       height=24,
       units='cm',
       scale = 0.8,
       dpi=600)

