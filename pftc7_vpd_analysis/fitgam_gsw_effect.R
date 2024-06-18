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

# Fit GAM models for Conductance (cond) with smoothing terms for fspecies, and tleaf
ww <- vec_mod_group %>% 
  lapply(function(sel_group) {
    tmp_dat <- dat[mod_group == sel_group]
    o <- tryCatch(
      {
        gam(photo ~ 
              s(fspecies, bs = 're') + 
              s(cond, k = 5, bs = 'ts') + 
              s(tleaf, k = 5, bs = 'ts'), 
            data = tmp_dat,
            select = TRUE,
            method = "REML")
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (is.null(o)) return(NULL)
    o$mod_group <- sel_group
    o$site <- unique(tmp_dat$site)
    return(o)
  })

# Remove failed models
ww <- ww[!sapply(ww, is.null)]

# Extract smooth estimates
ww2 <- ww %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    site <- x$site
    smooth_estimates <- tryCatch(
      {
        gratia::smooth_estimates(x, unconditional = TRUE, freq = FALSE, select = 's(cond', partial_match = TRUE)
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (is.null(smooth_estimates)) return(NULL)
    smooth_estimates %>%
      mutate(species = unique(x$model$species),
             R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist(fill = TRUE)


# Plotting data
colors_norway <- c("plum1", "violet", "violetred", "violetred4") 
colors_safrica <- c("lightblue1", "skyblue1", "dodgerblue","blue", "blue4")  
colors_combined <- c(colors_norway, colors_safrica)

p2 <- ww2 %>% 
  mutate(site = factor(site,
                       ordered = TRUE,
                       levels = c("6", "7", "8", "9", "1", "2", "3", "4", "5"),
                       labels = c("469 m", "700", "920", "1290", "2000 m", "2200 m", "2400 m", "2600 m", "2800 m"))) %>% 
  ggplot(aes(cond, .estimate, color = site, fill = site)) +
  geom_hline(aes(yintercept = 0), 
             col = 'grey40',
             lwd = 0.5,
             lty = 1) + 
  geom_ribbon(aes(cond,
                  ymin = .estimate - 1.96 * .se,
                  ymax = .estimate + 1.96 * .se, 
                  fill = site),
              color = NA,
              alpha = 0.25) +
  geom_line(lwd = 1) +
  labs(color = 'Site',
       fill = 'Site', 
       y = expression(paste("(Partial) effect on photosynthesis (mol ",
                            m^ -2, s^-1, ")")),
       x = "Conductance") +
  scale_color_manual(values = colors_combined) +
  scale_fill_manual(values = colors_combined) +
  theme_classic() +  # Add this line if `fn_theme()` is a custom theme function
  theme(legend.position = "bottom")
p2


ggsave(p2,
       filename=
         paste0("pftc7_vpd_analysis/figures/EffectOfConductanceOnA","_",
                Sys.Date(),
                ".png"),
       width=18,
       height=12,
       units='cm',
       device = grDevices::png,
       scale = 1,
       dpi=600)

################################################
################################################
# Multipanel site plot ===============================================
vv2 <- ww %>% 
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


vv3 <- ww %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(cond)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')
vv3$R2 %>% summary


qq <- vec_mod_group %>% 
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

qq2 <- qq %>% 
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
               s(cond,bs='ts',k=5),
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
    gratia::smooth_estimates(x,unconditional=T,smooth='s(cond)') %>% 
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
  rbindlist(list(vv2,qq2),fill=T) %>% 
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

## Temperature Response Figure =============
site_levels <- c("S.Africa 1","Norway 1", "S.Africa 2","Norway 2", "S.Africa 3","Norway 3","S.Africa 4", "Norway 4","S.Africa 5")
p_photo <- rbindlist(list(vv2,qq2),fill=T)
p_photo_plot <- ggplot(p_photo, aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_point() +
  geom_line(lwd = 1) +
  theme_classic() +
  geom_hline(aes(yintercept = 0), col = 'grey40', lwd = 0.5, lty = 1) + 
  geom_vline(data = fig_dat_topt, show.legend = FALSE, aes(xintercept = Topt, color = mod), lty = 2) +
  geom_ribbon(aes(ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = mod), color = NA, alpha = 0.25) +
  scale_color_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy"), labels = c("m_t" = "Temperature", "m_tv" = "Temperature + Conductance")) +
  scale_fill_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy"), labels = c("m_t" = "Temperature", "m_tv" = "Temperature + Conductance")) +
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(color = 'GAM', fill = 'GAM', y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m**-2, s**-1,")")), x = "Leaf temperature (°C)") +
  facet_wrap(~ factor(site, levels = c("1", "6", "2", "7", "3", "8", "4", "9", "5"),
                      labels = site_levels),
             scales = 'free_y', ncol = 2, shrink = TRUE) +
  theme(legend.position = 'bottom', strip.text = element_text(face='bold', hjust = 0)) +
  geom_vline(data = avg_weib_topt, show.legend = FALSE, aes(xintercept = avg_weib_topt), color = 'grey', lty = 2) +
  theme(legend.position = 'bottom', strip.text = element_text(face='bold', hjust = 0))
p_photo_plot

##############################
##############################
#Species facet plot:
##############################
# Fit the models and include species
ww <- vec_mod_group %>% 
  lapply(function(sel_group) {
    tmp_dat <- dat[mod_group == sel_group]
    o <- tryCatch(
      {
        gam(photo ~ 
              s(fspecies, bs = 're') + 
              s(cond, k = 5, bs = 'ts') + 
              s(tleaf, k = 5, bs = 'ts'), 
            data = tmp_dat,
            select = TRUE,
            method = "REML")
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (is.null(o)) return(NULL)
    o$mod_group <- sel_group
    o$species <- unique(tmp_dat$species)[1] # Ensure single species
    return(o)
  })

# Remove failed models
ww <- ww[!sapply(ww, is.null)]

# Extract smooth estimates for `cond`
ww2 <- ww %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    smooth_estimates <- tryCatch(
      {
        gratia::smooth_estimates(x, unconditional = TRUE, freq = FALSE, select = 's(cond)', partial_match = TRUE)
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (is.null(smooth_estimates)) return(NULL)
    smooth_estimates %>%
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist(fill = TRUE)

# Extract smooth estimates for `tleaf`
vv2 <- ww %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    gratia::smooth_estimates(x, unconditional = TRUE, smooth = 's(tleaf)') %>% 
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')

# Extract smooth estimates for `cond` again, if needed
vv3 <- ww %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    gratia::smooth_estimates(x, unconditional = TRUE, smooth = 's(cond)') %>% 
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')

# Repeat for qq and zz models
qq <- vec_mod_group %>% 
  lapply(function(sel_group) {
    tmp_dat <- dat[mod_group == sel_group]
    o <- gam(photo ~ 
               s(fspecies, bs = 're') + 
               s(tleaf, bs = 'ts', k = 5), 
             data = tmp_dat,
             select = TRUE,
             method = "REML")
    o$mod_group <- sel_group
    o$form <- unique(tmp_dat$form)
    o$species <- unique(tmp_dat$species)[1] # Ensure single species
    return(o)
  })

qq2 <- qq %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    gratia::smooth_estimates(x, unconditional = TRUE, smooth = 's(tleaf)') %>% 
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = "m_t")

zz <- vec_mod_group %>% 
  lapply(function(sel_group) {
    tmp_dat <- dat[mod_group == sel_group]
    o <- gam(photo ~ 
               s(fspecies, bs = 're') + 
               s(cond, bs = 'ts', k = 5),
             data = tmp_dat,
             select = TRUE,
             method = "REML")
    o$mod_group <- sel_group
    o$form <- unique(tmp_dat$form)
    o$species <- unique(tmp_dat$species)[1] # Ensure single species
    return(o)
  })

zz2 <- zz %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    gratia::smooth_estimates(x, unconditional = TRUE, smooth = 's(cond)') %>% 
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_v')
# Extract Topt
fig_dat_topt <- 
  rbindlist(list(vv2, qq2), fill = TRUE) %>% 
  group_by(mod, species) %>% 
  filter(.estimate == max(.estimate)) %>% 
  filter(.estimate < 41) %>% 
  ungroup() %>% 
  select(mod, species, tleaf) %>% 
  rename(Topt = tleaf)

avg_weib_topt <- dat %>%
  group_by(species) %>%
  summarize(avg_weib_topt = mean(weib.topt, na.rm = TRUE)) %>%
  ungroup()

spec_photo <- rbindlist(list(vv2,qq2),fill=T)

# Create the updated plot
spec_photo_plot <- ggplot(spec_photo, aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_point() +
  geom_line(lwd = 1) +
  theme_classic() +
  geom_hline(aes(yintercept = 0), col = 'grey40', lwd = 0.5, lty = 1) + 
  geom_vline(data = fig_dat_topt, show.legend = FALSE, aes(xintercept = Topt, color = mod), lty = 2) +
  geom_ribbon(aes(ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = mod), color = NA, alpha = 0.25) +
  scale_color_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy"), labels = c("m_t" = "Temperature", "m_tv" = "Temperature + Conductance")) +
  scale_fill_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy"), labels = c("m_t" = "Temperature", "m_tv" = "Temperature + Conductance")) +
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  labs(color = 'GAM', fill = 'GAM', y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m**-2, s**-1, ")")), x = "Leaf temperature (°C)") +
  facet_wrap(~ species, scales = 'free_y', ncol = 2, shrink = TRUE) +
  theme(legend.position = 'bottom', strip.text = element_text(face = 'bold', hjust = 0)) +
  geom_vline(data = avg_weib_topt, show.legend = FALSE, aes(xintercept = avg_weib_topt), color = 'grey', lty = 2) +
  theme(legend.position = 'bottom', strip.text = element_text(face = 'bold', hjust = 0))

spec_photo_plot
