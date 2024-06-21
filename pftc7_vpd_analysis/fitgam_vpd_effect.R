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

#Start by running the file: configure.datfile.R 

# Sub models ===============================================
vv <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group){
    tmp_dat <- dat[mod_group ==sel_group]
    o <- gam(photo ~ 
               s(fspecies, bs='re') +
               s(vpdl,k=5,bs='ts') + 
               s(tleaf,k=5,bs='ts'), 
             data=tmp_dat,
             select=T,
             method = "REML")
    o$mod_group <- sel_group
    o$site <- unique(tmp_dat$site)
    return(o)
  })

vv2 <- vv %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    # form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,freq=F, 
                             select='s(vpdl', 
                             partial_match = T) %>% 
      # mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist()


ww <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group){
    tmp_dat <- dat[mod_group ==sel_group]
    o <- gam(cond ~ 
               s(fspecies, bs='re') + 
               s(vpdl,k=5,bs='ts') + 
               s(tleaf,k=5,bs='ts'), 
             data=tmp_dat,
             select=T,
             method = "REML")
    o$mod_group <- sel_group
    o$site <- unique(tmp_dat$site)
    return(o)
  })
ww2 <- ww %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,
                             unconditional=T,freq=F,
                             select='s(vpdl',
                             partial_match = T) %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist()



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

## Make Figure =============
#006BA4
#FF800E
# cols4all::c4a_gui()
colors_norway <- c("plum1", "violet", "violetred", "violetred4") 
colors_safrica <- c("lightblue1", "skyblue1", "dodgerblue","blue", "blue4")  
colors_combined <- c(colors_norway, colors_safrica)
p1 <- vv2 %>% 
  mutate(site = factor(site,
                       ordered = T,
                       levels=c("6","7","8","9","1", "2", "3", "4", "5"),
                       labels=c("Norway 1","Norway 2", "Norway 3","Norway 4","S.Africa 1", "S.Africa 2", "S.Africa 3", "S.Africa 4", "S.Africa 5"))) %>% 
  ggplot(aes(vpdl, .estimate,color=site,fill=site))+
  geom_hline(aes(yintercept = 0), 
             col = 'grey40',
             lwd=0.5,
             lty=1) + 
  geom_ribbon(aes(vpdl,
                  ymin=.estimate-1.96*.se,
                  ymax=.estimate+1.96*.se, 
                  fill = site),
              color=NA,
              alpha=0.25)+
  geom_line(lwd=1)  +
  scale_color_manual(values = colors_combined) +
  scale_fill_manual(values = colors_combined) +
  labs(color='Site',
       fill='Site', 
       y = expression(paste("(Partial) effect on Photosynthesis (umol ",
                            CO[2]," ",
                            m**-2, s**-1,")")),
       # x = expression(paste(T[leaf]~"(°C)"))
       x = "VPD (kPa)") +
  fn_theme() +
  theme(legend.position = 'none'); p1

p2 <- ww2 %>% 
  mutate(site = factor(site,
                       ordered = T,
                       levels=c("6","7","8","9","1", "2", "3", "4", "5"),
                       labels=c("469 m","700","920","1290","2000 m", "2200 m", "2400 m", "2600 m", "2800 m"))) %>% 
  ggplot(aes(vpdl, .estimate,color=site,fill=site))+
  geom_hline(aes(yintercept = 0), 
             col = 'grey40',
             lwd=0.5,
             lty=1) + 
  geom_ribbon(aes(vpdl,
                  ymin=.estimate-1.96*.se,
                  ymax=.estimate+1.96*.se, 
                  fill = site),
              color=NA,
              alpha=0.25)+
  geom_line(lwd=1)  +
  labs(color='Site',
       fill='Site', 
       y = expression(paste("(Partial) effect on Conductance (mol ",
                            # CO[2]," ",
                            m**-2, s**-1,")")),
       # x = expression(paste(T[leaf]~"(°C)"))
       x = "VPD (kPa)") +
  scale_color_manual(values = colors_combined) +
  scale_fill_manual(values = colors_combined) +
  fn_theme() + theme(legend.position = "bottom"); p2

p_out <- (p1|p2) + plot_annotation(tag_levels = 'a', 
                                   tag_prefix = '(',
                                   tag_suffix = ')')
p_out



ggsave(p_out,
       filename=
         paste0("pftc7_vpd_analysis/figures/plot_photo_gs_site_tv_with-re_vpd-effect","_",
                Sys.Date(),
                ".png"),
       width=22,
       height=12,
       units='cm',
       device = grDevices::png,
       scale = 1,
       dpi=600)
