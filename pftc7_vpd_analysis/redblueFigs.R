########################################################
#Make multi-panel version for site (can comment in or out the species fixed effect)
########################################################

vec_mod_group <- unique(dat$mod_group)

# Model 1: Temperature + VPD
vvt <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[mod_group == sel_group]
  o <- gam(photo ~ factor(species)+
             #s(fspecies, bs='re') + 
             s(vpdl, bs='ts', k=5) + 
             s(tleaf, bs='ts', k=5), 
           data=tmp_dat,
           select=TRUE,
           method="REML")
  o$mod_group <- sel_group
  o$site <- unique(tmp_dat$site)
  return(o)
})

vv2t <- lapply(vvt, function(x) {
  smooth_estimates(x, unconditional=TRUE, smooth='s(tleaf)') %>%
    mutate(mod="m_tv", R2=summary(x)$r.sq, site=x$site)
}) %>% rbindlist()

# Model 2: Temperature
wwt <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[mod_group == sel_group]
  o <- gam(photo ~ factor(species)+
             #s(fspecies, bs='re') + 
             s(tleaf, bs='ts', k=5), 
           data=tmp_dat,
           select=TRUE,
           method="REML")
  o$mod_group <- sel_group
  o$site <- unique(tmp_dat$site)
  return(o)
})

ww2t <- lapply(wwt, function(x) {
  smooth_estimates(x, unconditional=TRUE, smooth='s(tleaf)') %>%
    mutate(mod="m_t", R2=summary(x)$r.sq, site=x$site)
}) %>% rbindlist()

# Model 3: Temperature + Cond
ww <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[mod_group == sel_group]
  o <- gam(photo ~ 
             s(fspecies, bs='re') + 
             s(cond, bs='ts', k=5) + 
             s(tleaf, bs='ts', k=5), 
           data=tmp_dat,
           select=TRUE,
           method="REML")
  o$mod_group <- sel_group
  o$site <- unique(tmp_dat$site)
  return(o)
})

vv2 <- lapply(ww, function(x) {
  smooth_estimates(x, unconditional=TRUE, smooth='s(tleaf)') %>%
    mutate(mod="m_tc", R2=summary(x)$r.sq, site=x$site)
}) %>% rbindlist()

# Combine All Data
combined_data <- rbindlist(list(vv2t, ww2t, vv2), fill=TRUE)

# Extract Topt
fig_dat_topt <- combined_data %>% 
  group_by(mod, site) %>% 
  filter(.estimate == max(.estimate)) %>% 
  filter(.estimate < 41) %>% 
  ungroup() %>% 
  select(mod, site, tleaf) %>% 
  rename(Topt = tleaf)

# Calculate Average Weibull Topt
avg_weib_topt <- dat %>%
  group_by(site) %>%
  summarize(avg_weib_topt = mean(weib.topt, na.rm = TRUE)) %>%
  ungroup()

# Define site levels
site_levels <- c("Norway 1", "Norway 2", "Norway 3", "Norway 4", 
                 "S.Africa 1", "S.Africa 2", "S.Africa 3", "S.Africa 4", "S.Africa 5")

# Plot
p <- combined_data %>%
  ggplot(aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_point() +
  geom_line(lwd = 1) +
  geom_hline(aes(yintercept = 0), col = 'grey40', lwd = 0.5, lty = 1) +
  geom_vline(data = fig_dat_topt, aes(xintercept = Topt, color = mod), lty = 2, show.legend = FALSE) +
  geom_ribbon(aes(ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = mod), color = NA, alpha = 0.25) +
  scale_color_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy", "m_tc" = "green"),
                     labels = c("m_t" = "Temperature", "m_tv" = "Temperature + VPD", "m_tc" = "Temperature + Conductance")) +
  scale_fill_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy", "m_tc" = "green"),
                    labels = c("m_t" = "Temperature", "m_tv" = "Temperature + VPD", "m_tc" = "Temperature + Conductance")) +
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(color = 'GAM', fill = 'GAM', 
       y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^2, s^-1, ")")), 
       x = "Leaf temperature (°C)") +
  facet_wrap(~ factor(site, levels = unique(dat$site), labels = site_levels),
             scales = 'free_y', ncol = 1, shrink = TRUE) +
  theme_classic() +
  theme(legend.position = 'bottom', strip.text = element_text(face = 'bold', hjust = 0)) +
  geom_vline(data = avg_weib_topt, aes(xintercept = avg_weib_topt), color = 'grey', lty = 2, show.legend = FALSE)

# Save the Plot
ggsave(p, 
       filename = paste0("pftc7_vpd_analysis/figures/threeGAM_site_fixSpec", Sys.Date(), ".png"),
       device = grDevices::png,
       width = 20,
       height = 35,
       units = 'cm',
       scale = 0.8,
       dpi = 600)



########################################################
#Make one panel version
########################################################
# Define unique model groups
vec_mod_group <- unique(dat$mod_group)

# Model 1: Temperature + VPD
vvt <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[mod_group == sel_group]
  o <- gam(photo ~ 
             #s(fspecies, bs='re') + 
             s(vpdl, bs='ts', k=5) + 
             s(tleaf, bs='ts', k=5), 
           data=tmp_dat,
           select=TRUE,
           method="REML")
  o$mod_group <- sel_group
  return(o)
})

vv2t <- lapply(vvt, function(x) {
  smooth_estimates(x, unconditional=TRUE, smooth='s(tleaf)') %>%
    mutate(mod="m_tv", R2=summary(x)$r.sq)
}) %>% rbindlist()

# Model 2: Temperature
wwt <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[mod_group == sel_group]
  o <- gam(photo ~ 
             s(tleaf, bs='ts', k=5), 
           data=tmp_dat,
           select=TRUE,
           method="REML")
  o$mod_group <- sel_group
  return(o)
})

ww2t <- lapply(wwt, function(x) {
  smooth_estimates(x, unconditional=TRUE, smooth='s(tleaf)') %>%
    mutate(mod="m_t", R2=summary(x)$r.sq)
}) %>% rbindlist()

# Model 3: Temperature + Cond
ww <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[mod_group == sel_group]
  o <- gam(photo ~ 
             #s(fspecies, bs='re') + 
             s(cond, bs='ts', k=5) + 
             s(tleaf, bs='ts', k=5), 
           data=tmp_dat,
           select=TRUE,
           method="REML")
  o$mod_group <- sel_group
  return(o)
})

vv2 <- lapply(ww, function(x) {
  smooth_estimates(x, unconditional=TRUE, smooth='s(tleaf)') %>%
    mutate(mod="m_tc", R2=summary(x)$r.sq)
}) %>% rbindlist()

# Combine All Data
combined_data <- rbindlist(list(vv2t, ww2t, vv2), fill=TRUE)

# Extract Topt
fig_dat_topt <- combined_data %>% 
  group_by(mod) %>% 
  filter(.estimate == max(.estimate)) %>% 
  filter(.estimate < 41) %>% 
  ungroup() %>% 
  select(mod, tleaf) %>% 
  rename(Topt = tleaf)

# Calculate Average Weibull Topt
avg_weib_topt <- dat %>%
  summarize(avg_weib_topt = mean(weib.topt, na.rm = TRUE)) %>%
  ungroup()

# Ensure data is ordered by tleaf
combined_data <- combined_data %>%
  arrange(tleaf)

# Plot without faceting
p_single <- combined_data %>%
  ggplot(aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"), se = TRUE, size = 1, alpha = 0.25) +
  geom_hline(aes(yintercept = 0), col = 'grey40', size = 0.5, linetype = 1) +
  geom_vline(data = fig_dat_topt, aes(xintercept = Topt, color = mod), linetype = 2, show.legend = FALSE) +
  scale_color_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy", "m_tc" = "green"),
                     labels = c("m_t" = "Temperature", "m_tv" = "Temperature + VPD", "m_tc" = "Temperature + Conductance")) +
  scale_fill_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy", "m_tc" = "green"),
                    labels = c("m_t" = "Temperature", "m_tv" = "Temperature + VPD", "m_tc" = "Temperature + Conductance")) +
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(color = 'GAM', fill = 'GAM', 
       y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")), 
       x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = 'bottom') +
  geom_vline(data = avg_weib_topt, aes(xintercept = avg_weib_topt), color = 'grey', linetype = 2, show.legend = FALSE)

# Save the Plot
ggsave(p_single, 
       filename = paste0("pftc7_vpd_analysis/figures/threeGAM", Sys.Date(), ".png"),
       device = grDevices::png,
       width = 20,
       height = 20,
       units = 'cm',
       scale = 0.8,
       dpi = 600)




###############################################################
#Make one panel version with a line for schoolfield + cond too
###############################################################
# Define data and start parameters
vec_mod_group <- unique(dat$mod_group)
start_params <- c(J_ref = 35, E = 0.5, E_D = 1.5, T_opt = 20)
# Fit Schoolfield GAM
ssmod <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  fit_results <- fit_schoolfield_gam("tleaf", "photo", "cond", T_ref = 25, start_params = start_params, data = tmp_dat)
  return(list(mod_group = sel_group, gam_fit = fit_results$gam_fit))
})

# Extract smooth estimates for Schoolfield GAM
ssmod_sm <- lapply(ssmod, function(x) {
  smooth_estimates(x$gam_fit, unconditional = TRUE,smooth = 's(tleaf)', partial_match = TRUE) %>%
    mutate(mod = "ss", R2 = summary(x$gam_fit)$r.sq)
}) %>% rbindlist()



# Fit Schoolfield GAM
ssmod.tleaf <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  fit_results <- fit_schoolfield_gam("tleaf", "photo", T_ref = 25, start_params = start_params, data = tmp_dat)
  return(list(mod_group = sel_group, gam_fit = fit_results$gam_fit))
})


# Extract smooth estimates for Schoolfield GAM
ssmod_sm_tleaf <- lapply(ssmod.tleaf, function(x) {
  smooth_estimates(x$gam_fit, unconditional = TRUE,smooth = 's(tleaf)', partial_match = TRUE) %>%
    mutate(mod = "sst", R2 = summary(x$gam_fit)$r.sq)
}) %>% rbindlist()



# Model 1: Temperature + VPD
vvt <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  o <- gam(photo ~ 
             s(vpdl, bs='ts', k=5) + 
             s(tleaf, bs='ts', k=5), 
           data = tmp_dat,
           select = TRUE,
           method = "REML")
  o$mod_group <- sel_group
  return(o)
})

vv2t <- lapply(vvt, function(x) {
  smooth_estimates(x, unconditional = TRUE, smooth = 's(tleaf)', partial_match = TRUE) %>%
    mutate(mod = "m_tv", R2 = summary(x)$r.sq)
}) %>% rbindlist()

# Model 2: Temperature
wwt <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  o <- gam(photo ~ 
             s(tleaf, bs='ts', k=5), 
           data = tmp_dat,
           select = TRUE,
           method = "REML")
  o$mod_group <- sel_group
  return(o)
})

ww2t <- lapply(wwt, function(x) {
  smooth_estimates(x, unconditional = TRUE, smooth = 's(tleaf)', partial_match = TRUE) %>%
    mutate(mod = "m_t", R2 = summary(x)$r.sq)
}) %>% rbindlist()

# Model 3: Temperature + Cond
ww <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  o <- gam(photo ~ 
             s(cond, bs='ts', k=5) + 
             s(tleaf, bs='ts', k=5), 
           data = tmp_dat,
           select = TRUE,
           method = "REML")
  o$mod_group <- sel_group
  return(o)
})

vv2 <- lapply(ww, function(x) {
  smooth_estimates(x, unconditional = TRUE, smooth = 's(tleaf)', partial_match = TRUE) %>%
    mutate(mod = "m_tc", R2 = summary(x)$r.sq)
}) %>% rbindlist()

# Combine all data
combined_data <- rbindlist(list(vv2t, ww2t, vv2, ssmod_sm, ssmod_sm_tleaf), fill = TRUE)

# Function to find local maxima using LOESS
find_local_maxima <- function(data) {
  loess_fit <- loess(.estimate ~ tleaf, data = data, span = 0.3)
  smoothed <- predict(loess_fit)
  peaks <- which(diff(sign(diff(smoothed))) == -2) + 1
  if (length(peaks) == 0) return(NA)
  return(data$tleaf[peaks[which.max(smoothed[peaks])]])
}

# Extract Topt
fig_dat_topt <- combined_data %>% 
  group_by(mod) %>% 
  nest() %>%
  mutate(Topt = map_dbl(data, find_local_maxima)) %>%
  select(mod, Topt) %>%
  ungroup()
#fig_dat_topt <- fig_dat_topt %>%
#  mutate(Topt = ifelse(mod == "sst", 26.6, Topt))

# Calculate Average Weibull Topt
avg_weib_topt <- dat %>%
  summarize(avg_weib_topt = mean(weib.topt, na.rm = TRUE)) %>%
  ungroup()

# Ensure data is ordered by tleaf
combined_data <- combined_data %>%
  group_by(mod)%>%
  arrange(tleaf)

# Plot without faceting
p_single <- combined_data %>%
  ggplot(aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"), se = TRUE, size = 1, alpha = 0.25) +
  geom_hline(aes(yintercept = 0), col = 'grey40', size = 0.5, linetype = 1) +
  geom_vline(data = fig_dat_topt, aes(xintercept = Topt, color = mod), linetype = 2, show.legend = FALSE, size=1) +
  scale_color_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy", "m_tc" = "green", "sst" = "orange","ss" = "maroon"),
                     labels = c("m_t" = "Temperature", "m_tv" = "Temperature + VPD", "m_tc" = "Temperature + gsw", "sst" = "Schoolfield + Tleaf","ss" = "Schoolfield + Tleaf + gsw")) +
  scale_fill_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy", "m_tc" = "green", "sst" = "orange","ss" = "maroon"),
                    labels = c("m_t" = "Temperature", "m_tv" = "Temperature + VPD", "m_tc" = "Temperature + gsw", "sst" = "Schoolfield + Tleaf","ss" = "Schoolfield + Tleaf + gsw")) +
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-20, 15))+
  labs(color = 'GAM', fill = 'GAM', 
       y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")), 
       x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = 'bottom') +
  geom_vline(data = avg_weib_topt, aes(xintercept = avg_weib_topt), color = 'grey', linetype = 2, show.legend = FALSE, size=1)

# Save the Plot
ggsave(p_single, 
       filename = paste0("pftc7_vpd_analysis/figures/fiveGAM", Sys.Date(), ".png"),
       device = grDevices::png,
       width = 40,
       height = 40,
       units = 'cm',
       scale = 0.8,
       dpi = 600)

#############################################
#Plot above but incrementally for presentations
#############################################
(redblueplot1 <- combined_data %>%
  filter(mod == "m_t") %>%
  ggplot(aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"), se = TRUE, size = 1, alpha = 0.25) +
  geom_hline(aes(yintercept = 0), col = 'grey40', size = 0.5, linetype = 1) +
  geom_vline(data = fig_dat_topt %>% filter(mod == "m_t"), aes(xintercept = Topt, color = mod), linetype = 2, show.legend = FALSE, size=1) +
  scale_color_manual(values = c("m_t" = "#cf0000"),labels = c("m_t" = "Temperature"))+
  scale_fill_manual(values = c("m_t" = "#cf0000"),labels = c("m_t" = "Temperature"))+
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-20, 15))+
  labs(color = 'GAM', fill = 'GAM', 
       y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")), 
       x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = 'bottom') +
  geom_vline(data = avg_weib_topt, aes(xintercept = avg_weib_topt), color = 'grey', linetype = 2, show.legend = FALSE, size=1))

(redblueplot2 <- combined_data %>%
    filter(mod %in% c("m_t", "m_tc")) %>%
    ggplot(aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"), se = TRUE, size = 1, alpha = 0.25) +
    geom_hline(aes(yintercept = 0), col = 'grey40', size = 0.5, linetype = 1) +
    geom_vline(data = fig_dat_topt %>% filter(mod %in% c("m_t", "m_tc")), aes(xintercept = Topt, color = mod), linetype = 2, show.legend = FALSE, size=1) +
    scale_color_manual(values = c("m_t" = "#cf0000", "m_tc" = "green"),
                       labels = c("m_t" = "Temperature", "m_tc" = "Temperature + gsw")) +
    scale_fill_manual(values = c("m_t" = "#cf0000", "m_tc" = "green"),
                      labels = c("m_t" = "Temperature", "m_tc" = "Temperature + gsw")) +
    scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(-20, 15)) +
    labs(color = 'GAM', fill = 'GAM', 
         y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")), 
         x = "Leaf temperature (°C)") +
    theme_classic() +
    theme(legend.position = 'bottom') +
    geom_vline(data = avg_weib_topt, aes(xintercept = avg_weib_topt), color = 'grey', linetype = 2, show.legend = FALSE, size=1))

(redblueplot3 <- combined_data %>%
    filter(mod %in% c("m_t", "m_tc", "m_tv")) %>%
    ggplot(aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"), se = TRUE, size = 1, alpha = 0.25) +
    geom_hline(aes(yintercept = 0), col = 'grey40', size = 0.5, linetype = 1) +
    geom_vline(data = fig_dat_topt %>% filter(mod %in% c("m_t", "m_tc", "m_tv")), aes(xintercept = Topt, color = mod), linetype = 2, show.legend = FALSE, size=1) +
    scale_color_manual(values = c("m_t" = "#cf0000", "m_tc" = "green", "m_tv" = "navy"),
                       labels = c("m_t" = "Temperature", "m_tc" = "Temperature + gsw", "m_tv" = "Temperature + VPD")) +
    scale_fill_manual(values = c("m_t" = "#cf0000", "m_tc" = "green", "m_tv" = "navy"),
                      labels = c("m_t" = "Temperature", "m_tc" = "Temperature + gsw", "m_tv" = "Temperature + VPD")) +
    scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(-20, 15)) +
    labs(color = 'GAM', fill = 'GAM', 
         y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")), 
         x = "Leaf temperature (°C)") +
    theme_classic() +
    theme(legend.position = 'bottom') +
    geom_vline(data = avg_weib_topt, aes(xintercept = avg_weib_topt), color = 'grey', linetype = 2, show.legend = FALSE, size=1))


(redblueplot4 <- combined_data %>%
    filter(mod %in% c("m_t", "m_tc", "m_tv", "sst")) %>%
    ggplot(aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"), se = TRUE, size = 1, alpha = 0.25) +
    geom_hline(aes(yintercept = 0), col = 'grey40', size = 0.5, linetype = 1) +
    geom_vline(data = fig_dat_topt %>% filter(mod %in% c("m_t", "m_tc", "m_tv", "sst")), aes(xintercept = Topt, color = mod), linetype = 2, show.legend = FALSE, size=1) +
    scale_color_manual(values = c("m_t" = "#cf0000", "m_tc" = "green", "m_tv" = "navy", "sst" = "orange"),
                       labels = c("m_t" = "Temperature", "m_tc" = "Temperature + gsw", "m_tv" = "Temperature + VPD", "sst" = "Schoolfield + Tleaf")) +
    scale_fill_manual(values = c("m_t" = "#cf0000", "m_tc" = "green", "m_tv" = "navy", "sst" = "orange"),
                      labels = c("m_t" = "Temperature", "m_tc" = "Temperature + gsw", "m_tv" = "Temperature + VPD", "sst" = "Schoolfield + Tleaf")) +
    scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(-20, 15)) +
    labs(color = 'GAM', fill = 'GAM', 
         y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")), 
         x = "Leaf temperature (°C)") +
    theme_classic() +
    theme(legend.position = 'bottom') +
    geom_vline(data = avg_weib_topt, aes(xintercept = avg_weib_topt), color = 'grey', linetype = 2, show.legend = FALSE, size=1))


##Then plot the full one at the start of this section::



#############################################
#Plot one panel fig but with cond on x axis:
#############################################
# Model 3: Temperature + Cond
# Extract smooth estimates for Schoolfield GAM
ssmod_sm2 <- lapply(ssmod, function(x) {
  smooth_estimates(x$gam_fit, unconditional = TRUE,smooth = 's(cond)', partial_match = TRUE) %>%
    mutate(mod = "ss", R2 = summary(x$gam_fit)$r.sq)
}) %>% rbindlist()


TCmod <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  o <- gam(photo ~ 
             s(cond, bs='ts', k=5) + 
             s(tleaf, bs='ts', k=5), 
           data = tmp_dat,
           select = TRUE,
           method = "REML")
  o$mod_group <- sel_group
  return(o)
})

TCmod_est <- lapply(TCmod, function(x) {
  smooth_estimates(x, unconditional = TRUE, smooth = 's(cond)', partial_match = TRUE) %>%
    mutate(mod = "m_tc", R2 = summary(x)$r.sq)
}) %>% rbindlist()
combo_data <- rbindlist(list(TCmod_est,ssmod_sm2), fill = TRUE)


Cmod <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  o <- gam(photo ~ 
             s(cond, bs='ts', k=5), 
           data = tmp_dat,
           select = TRUE,
           method = "REML")
  o$mod_group <- sel_group
  return(o)
})

Cmod_est <- lapply(Cmod, function(x) {
  smooth_estimates(x, unconditional = TRUE, smooth = 's(cond)', partial_match = TRUE) %>%
    mutate(mod = "m_c", R2 = summary(x)$r.sq)
}) %>% rbindlist()
combo_data <- rbindlist(list(TCmod_est,ssmod_sm2, Cmod_est), fill = TRUE)

# Plot without faceting
p_condplot <- combo_data %>%
  ggplot(aes(x = cond, y = .estimate, color = mod, fill = mod)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"), se = TRUE, size = 1, alpha = 0.25) +
  geom_hline(aes(yintercept = 0), col = 'grey40', size = 0.5, linetype = 1) +
  scale_color_manual(values = c("m_tc" = "green", "ss" = "maroon", "m_c"="purple4"),
                     labels = c("m_tc" = "Temperature + gsw", "ss" = "Schoolfield + Tleaf + gsw", "m_c" = "gsw")) +
  scale_fill_manual(values = c("m_tc" = "green", "ss" = "maroon", "m_c"="purple4"),
                    labels = c("m_tc" = "Temperature + gsw", "ss" = "Schoolfield + Tleaf + gsw", "m_c" = "gsw")) +
  labs(color = 'GAM', fill = 'GAM', 
       y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")), 
       x = "gsw") +
  theme_classic() +
  theme(legend.position = 'bottom') 

# Save the Plot
ggsave(p_condplot, 
       filename = paste0("pftc7_vpd_analysis/figures/partialeffets_cond", Sys.Date(), ".png"),
       device = grDevices::png,
       width = 40,
       height = 40,
       units = 'cm',
       scale = 0.8,
       dpi = 600)
