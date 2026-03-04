# --- Fit 3 GAMs across mod_groups ---
vec_mod_group <- unique(dat$mod_group)

# Model 1: Temperature + VPD
vvt <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  o <- gam(photo ~ 
             s(vpdl, bs='ts', k=5) + 
             s(tleaf, bs='ts', k=5), 
           data = tmp_dat, select = TRUE, method = "REML")
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
           data = tmp_dat, select = TRUE, method = "REML")
  o$mod_group <- sel_group
  return(o)
})

ww2t <- lapply(wwt, function(x) {
  smooth_estimates(x, unconditional = TRUE, smooth = 's(tleaf)', partial_match = TRUE) %>%
    mutate(mod = "m_t", R2 = summary(x)$r.sq)
}) %>% rbindlist()

# Model 3: Temperature + Conductance
ww <- lapply(vec_mod_group, function(sel_group) {
  tmp_dat <- dat[dat$mod_group == sel_group, ]
  o <- gam(photo ~ 
             s(cond, bs='ts', k=5) + 
             s(tleaf, bs='ts', k=5), 
           data = tmp_dat, select = TRUE, method = "REML")
  o$mod_group <- sel_group
  return(o)
})

vv2 <- lapply(ww, function(x) {
  smooth_estimates(x, unconditional = TRUE, smooth = 's(tleaf)', partial_match = TRUE) %>%
    mutate(mod = "m_tc", R2 = summary(x)$r.sq)
}) %>% rbindlist()

# --- Combine, find Topt, compute avg Schoolfield Topt ---
combined_data <- rbindlist(list(vv2t, ww2t, vv2), fill = TRUE) %>%
  group_by(mod) %>%
  arrange(tleaf)

find_local_maxima <- function(data) {
  loess_fit <- loess(.estimate ~ tleaf, data = data, span = 0.3)
  smoothed <- predict(loess_fit)
  peaks <- which(diff(sign(diff(smoothed))) == -2) + 1
  if (length(peaks) == 0) return(NA)
  return(data$tleaf[peaks[which.max(smoothed[peaks])]])
}

fig_dat_topt <- combined_data %>% 
  group_by(mod) %>% 
  nest() %>%
  mutate(Topt = map_dbl(data, find_local_maxima)) %>%
  select(mod, Topt) %>%
  ungroup()

avg_school_topt <- dat %>%
  summarize(avg_school_topt = mean(school.topt, na.rm = TRUE))

# --- Plot ---
redblueplot3 <- combined_data %>%
  ggplot(aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ts"), se = TRUE, size = 1, alpha = 0.25) +
  geom_hline(aes(yintercept = 0), col = 'grey40', size = 0.5, linetype = 1) +
  geom_vline(data = fig_dat_topt, aes(xintercept = Topt, color = mod), linetype = 2, show.legend = FALSE, size = 1) +
  geom_vline(data = avg_school_topt, aes(xintercept = avg_school_topt), color = 'grey2', linetype = 2, show.legend = FALSE, size = 1) +
  scale_color_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy", "m_tc" = "green"),
                     labels = c("m_t" = "Temperature", "m_tv" = "Temperature + VPD", "m_tc" = "Temperature + gsw")) +
  scale_fill_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy", "m_tc" = "green"),
                    labels = c("m_t" = "Temperature", "m_tv" = "Temperature + VPD", "m_tc" = "Temperature + gsw")) +
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-15, 9)) +
  labs(color = 'GAM', fill = 'GAM', 
       y = expression(paste("Corrected photosynthesis (µmol ", CO[2], " ", m^-2, s^-1, ")")), 
       x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = 'bottom')


redblueplot3 <- redblueplot3 +
  coord_cartesian(ylim = c(-15, 9), clip = "off") +
  theme(plot.margin = margin(5.5, 90, 5.5, 5.5, "pt")) +
  annotate("text", x = 41, y = 0, label = "Concurvity = 100%",
           hjust = 0, size = 3.5) +
  annotate("text", x = 41, y = -1.7, label = "Concurvity = 98.7%",
           hjust = 0, size = 3.5) +
  annotate("text", x = 41, y = -7, label = "Concurvity = 0%",
           hjust = 0, size = 3.5)

redblueplot3
