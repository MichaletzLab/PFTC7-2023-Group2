raw.dat <- read.csv('data/raw.discardHooks_data.csv')
raw.dat <- raw.dat%>%
  filter(gsw<0.6, gsw>-1)%>%
  filter(!is.na(Tair), !is.na(A), !is.na(Species))
#### A vs. Tleaf by elevation ####
ggplot(raw.dat, aes(x = Tleaf, y = A, color=Species)) +
  geom_point(alpha = 0.01) +
  labs(x = "Leaf temperature (°C)",
       y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
       color = "Species") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(
    data = raw.dat,
    aes(x = Tair, y = A),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3,
    inherit.aes = FALSE) +
  theme_classic() +
  facet_wrap(~Elevation) +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))

#### gs vs. Tleaf by elevation ####
ggplot(raw.dat, aes(x = Tleaf, y = gsw, color=Species)) +
  geom_point(alpha = 0.01) +
  labs(x = "Leaf temperature (°C)",
       y = expression(gsw~(mol~m^{-2}~s^{-1})),
       color = "Species") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(
    data = raw.dat,
    aes(x = Tair, y = gsw),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3,
    inherit.aes = FALSE) +
  theme_classic() +
  facet_wrap(~Elevation) +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))

#### A vs. gs by elevation ####
ggplot(raw.dat, aes(x = gsw, y = A, color=Species)) +
  geom_point(alpha = 0.01) +
  labs(y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
       x = expression(gsw~(mol~m^{-2}~s^{-1})),
       color = "Species") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(
    data = raw.dat,
    aes(x = gsw, y = A),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3,
    inherit.aes = FALSE) +
  theme_classic() +
  facet_wrap(~Elevation) +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))

#### A vs. Tleaf by species ####
ggplot(raw.dat, aes(x = Tleaf, y = A, color=as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  labs(x = "Leaf temperature (°C)",
       y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
       color = "Elevation (masl)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(
    data = raw.dat,
    aes(x = Tair, y = A),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3,
    inherit.aes = FALSE) +
  theme_classic() +
  facet_wrap(~Species) +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))

#### gs vs. Tleaf by species ####
ggplot(raw.dat, aes(x = Tleaf, y = gsw, color=as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  labs(x = "Leaf temperature (°C)",
       y = expression(gsw~(mol~m^{-2}~s^{-1})),
       color = "Elevation (masl)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(
    data = raw.dat,
    aes(x = Tair, y = gsw),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3,
    inherit.aes = FALSE) +
  theme_classic() +
  facet_wrap(~Species) +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))

#### A vs. gs by species ####
ggplot(raw.dat, aes(x = gsw, y = A, color=as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  labs(y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
       x = expression(gsw~(mol~m^{-2}~s^{-1})),
       color = "Elevation (masl)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(
    data = raw.dat,
    aes(x = gsw, y = A),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3,
    inherit.aes = FALSE) +
  theme_classic() +
  facet_wrap(~Species) +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))

#### use dat.VPD from VPD_SW_8.22.25.R for A vs. gsw by VPD ####
ggplot(dat_VPD, aes(x = gsw, y = A, color = label)) +
  geom_point(alpha = 0.05, size = 1) +
  #geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_viridis_d(
    name = "VPD (kPa)\n(n obs)", 
    direction = -1,
    guide = guide_legend(
      override.aes = list(alpha = 1, size = 3)  # full color, bigger points in legend
    )
  ) +
  labs(x = expression(gsw~(mol~m^{-2}~s^{-1})),
       y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1}))) +
  theme_classic() +
  theme(legend.position = "right")#+
  #facet_wrap(~Elevation) #change to Species or Elevation or comment out

ggplot(dat_VPD, aes(x = Tleaf, y = gsw, color = label)) +
  geom_point(alpha = 0.01, size = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_viridis_d(name = "VPD (kPa)\n(n obs)", direction=-1) +
  labs(x = "Leaf temperature (°C)",
       y = expression(gsw~(mol~m^{-2}~s^{-1}))) +
  theme_classic() +
  theme(legend.position = "right") + 
  facet_wrap(~Elevation)

#### A vs. gsw by aspect ####
ggplot(raw.dat, aes(x = gsw, y = A, color=as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  labs(y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
       x = expression(gsw~(mol~m^{-2}~s^{-1})),
       color = "Elevation") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(
    data = raw.dat,
    aes(x = gsw, y = A),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3,
    inherit.aes = FALSE) +
  theme_classic() +
  facet_wrap(~Aspect) +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))

# iWUE: A/gsw vs. Tleaf ####
#Just A:
A.Tleaf.plot <- ggplot(raw.dat, aes(x = Tleaf, y = A, color = as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  labs(
    x = "Leaf temperature (°C)",
    y= expression(A~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    color = "Elevation (masl)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 20.66579, linetype="dashed")+
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3
  ) +
  theme_classic() +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))
#Just gsw:
gsw.Tleaf.plot <- ggplot(raw.dat, aes(x = Tleaf, y = gsw, color = as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  labs(
    x = "Leaf temperature (°C)",
    y= expression(gsw~(mol~m^{-2}~s^{-1})),
    color = "Elevation (masl)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 20.66579, linetype="dashed")+
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3
  ) +
  theme_classic() +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))
#Calculate average Topt and plot as vline below:
mean(discard.weibull$T_opt)
#iWUE
iWUE.Tleaf.plot <-ggplot(raw.dat, aes(x = Tleaf, y = A/gsw, color = as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  labs(
    x = "Leaf temperature (°C)",
    y = expression(iWUE~(mu*mol~CO[2]~m^{-2}~s^{-1}~ "/" ~mol~H[2]*O~m^{-2}~s^{-1})),
    color = "Elevation (masl)"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 20.66579, linetype="dashed")+
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    color = "black",
    linewidth = 1,
    alpha = 0.3
  ) +
  theme_classic() +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))

A.Tleaf.plot <- A.Tleaf.plot + theme(legend.position = "none")
gsw.Tleaf.plot <- gsw.Tleaf.plot + theme(legend.position = "none")

a.gsw.plot <- ggarrange(A.Tleaf.plot, gsw.Tleaf.plot, nrow = 2, ncol = 1, labels = c("A","B"))
ggarrange(a.gsw.plot, iWUE.Tleaf.plot, nrow = 1, ncol = 2,common.legend = TRUE,legend = "right", labels = c("", "C"))

# Make Tleaf binned dataset ####
# make Tleaf bins based on quantiles
dat_Tleaf <- at.subset3 %>%
  mutate(Tleaf_bin = ntile(Tleaf, 12)) %>%
  filter(gsw < 0.6, gsw > -1) %>%
  filter(!is.na(Tair), !is.na(A), !is.na(Species))
# summarize counts per bin
Tleaf_counts <- dat_Tleaf %>%
  group_by(Tleaf_bin) %>%
  summarise(
    n = n(),
    min_val = min(Tleaf, na.rm = TRUE),
    max_val = max(Tleaf, na.rm = TRUE),
    .groups = "drop")
# create nice labels with actual Tleaf ranges + counts
Tleaf_labels <- Tleaf_counts %>%
  mutate(label = paste0(round(min_val, 2), "–", round(max_val, 2), " (", n, ")")) %>%
  select(Tleaf_bin, label)
# join labels back to dataframe
dat_Tleaf <- dat_Tleaf %>%
  left_join(Tleaf_labels, by = "Tleaf_bin")

# Make a plot of gsw and A vs. Tleaf ####
# Calculate scaling to match ranges of A and gsw
range_A <- range(dat_Tleaf$A, na.rm = TRUE)
range_gsw <- range(dat_Tleaf$gsw, na.rm = TRUE)

# Define a linear transformation: scale gsw to A scale
scale_gsw_to_A <- function(x) {
  (x - range_gsw[1]) / diff(range_gsw) * diff(range_A) + range_A[1]
}

scale_A_to_gsw <- function(x) {
  (x - range_A[1]) / diff(range_A) * diff(range_gsw) + range_gsw[1]
}

# Fit GAM for Photosynthesis
gam_A <- gam(A ~ s(Tleaf, k = 10), data = dat_Tleaf)

# Fit GAM for gsw
gam_gsw <- gam(gsw ~ s(Tleaf, k = 10), data = dat_Tleaf)

# Generate prediction grids
pred_grid <- data.frame(Tleaf = seq(min(dat_Tleaf$Tleaf, na.rm=TRUE),
                                    max(dat_Tleaf$Tleaf, na.rm=TRUE),
                                    length.out = 200))

pred_grid$A_pred <- predict(gam_A, newdata = pred_grid)
pred_grid$gsw_pred <- predict(gam_gsw, newdata = pred_grid)

# Transform gsw to A scale for plotting
pred_grid$gsw_pred_scaled <- scale_gsw_to_A(pred_grid$gsw_pred)

# Plot
ggplot(dat_Tleaf, aes(x = Tleaf)) +
  # Points
  geom_point(aes(y = A, color = "Photosynthesis"), alpha = 0.05, size = 2) +
  geom_point(aes(y = scale_gsw_to_A(gsw), color = "gsw"), alpha = 0.01, size = 2) +
  # GAM lines
  geom_line(data = pred_grid, aes(y = A_pred, color = "Photosynthesis"), size = 1.2) +
  geom_line(data = pred_grid, aes(y = gsw_pred_scaled, color = "gsw"), size = 1.2) +
  scale_y_continuous(
    name = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    sec.axis = sec_axis(~scale_A_to_gsw(.), name = expression(g[sw]~(mol~m^{-2}~s^{-1})))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis" = "green3", "gsw" = "lightblue")
  ) +
  labs(x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

## plot linear for categories of Tleaf:
ggplot(dat_Tleaf, aes(x = Tleaf)) +
  # Points
  geom_point(aes(y = A, color = "Photosynthesis"), alpha = 0.05, size = 2) +
  geom_point(aes(y = scale_gsw_to_A(gsw), color = "gsw"), alpha = 0.01, size = 2) +
  # Linear regression lines for each Tleaf bin
  geom_smooth(aes(y = A, group = Tleaf_bin, color = "Photosynthesis"), 
              method = "lm", se = FALSE, size = 1.2) +
  geom_smooth(aes(y = scale_gsw_to_A(gsw), group = Tleaf_bin, color = "gsw"), 
              method = "lm", se = FALSE, size = 1.2) +
  # Axes and scaling
  scale_y_continuous(
    name = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    sec.axis = sec_axis(~scale_A_to_gsw(.), 
                        name = expression(g[sw]~(mol~m^{-2}~s^{-1})))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis" = "green3", "gsw" = "lightblue")
  ) +
  labs(x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )
# Plot A and E vs Tleaf ####
# Define scaling functions between A and E
scale_E_to_A <- function(E) {
  # put E on A scale
  rng_A <- range(dat_Tleaf$A, na.rm = TRUE)
  rng_E <- range(dat_Tleaf$E, na.rm = TRUE)
  (E - rng_E[1]) / diff(rng_E) * diff(rng_A) + rng_A[1]
}

scale_A_to_E <- function(A_scaled) {
  rng_A <- range(dat_Tleaf$A, na.rm = TRUE)
  rng_E <- range(dat_Tleaf$E, na.rm = TRUE)
  (A_scaled - rng_A[1]) / diff(rng_A) * diff(rng_E) + rng_E[1]
}

# Fit GAMs
gam_A <- gam(A ~ s(Tleaf, k = 10), data = dat_Tleaf)
gam_E <- gam(E ~ s(Tleaf, k = 10), data = dat_Tleaf)

# Prediction grid
pred_grid <- data.frame(
  Tleaf = seq(min(dat_Tleaf$Tleaf, na.rm=TRUE),
              max(dat_Tleaf$Tleaf, na.rm=TRUE),
              length.out = 200)
)

# Predict with standard errors
pred_A <- predict(gam_A, newdata = pred_grid, se.fit = TRUE)
pred_E <- predict(gam_E, newdata = pred_grid, se.fit = TRUE)

pred_grid$A_pred <- pred_A$fit
pred_grid$A_se <- pred_A$se.fit

pred_grid$E_pred <- pred_E$fit
pred_grid$E_se <- pred_E$se.fit
pred_grid$E_pred_scaled <- scale_E_to_A(pred_grid$E_pred)
pred_grid$E_se_scaled <- pred_grid$E_se / diff(range(dat_Tleaf$E, na.rm=TRUE)) * diff(range(dat_Tleaf$A, na.rm=TRUE))  # scale SE the same way as E

# Plot with shaded confidence bands
ggplot(dat_Tleaf, aes(x = Tleaf)) +
  geom_point(aes(y = A, color = "Photosynthesis"), alpha = 0.01, size = 2) +
  geom_point(aes(y = scale_E_to_A(E), color = "Transpiration"), alpha = 0.01, size = 2) +
  # Confidence bands
  geom_ribbon(data = pred_grid, aes(ymin = A_pred - A_se, ymax = A_pred + A_se), 
              fill = "green3", alpha = 0.4) +
  geom_ribbon(data = pred_grid, aes(ymin = E_pred_scaled - E_se_scaled, ymax = E_pred_scaled + E_se_scaled),
              fill = "orange", alpha = 0.4) +
  # GAM lines
  geom_line(data = pred_grid, aes(y = A_pred, color = "Photosynthesis"), size = 1.2) +
  geom_line(data = pred_grid, aes(y = E_pred_scaled, color = "Transpiration"), size = 1.2) +
  scale_y_continuous(
    name = expression(A~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    sec.axis = sec_axis(~scale_A_to_E(.), name = expression(E~(mmol~m^{-2}~s^{-1})))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis" = "green3", "Transpiration" = "orange")
  ) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  labs(x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = "right")

# Make the same as above but faceted by Elevation ####
# --- Define scaling functions (unchanged) ---
scale_E_to_A <- function(E, dat) {
  rng_A <- range(dat$A, na.rm = TRUE)
  rng_E <- range(dat$E, na.rm = TRUE)
  (E - rng_E[1]) / diff(rng_E) * diff(rng_A) + rng_A[1]
}

scale_A_to_E <- function(A_scaled, dat) {
  rng_A <- range(dat$A, na.rm = TRUE)
  rng_E <- range(dat$E, na.rm = TRUE)
  (A_scaled - rng_A[1]) / diff(rng_A) * diff(rng_E) + rng_E[1]
}

# --- Fit GAMs with Elevation as grouping factor ---
gam_A <- gam(A ~ s(Tleaf, by = Elevation) + Elevation, data = dat_Tleaf)
gam_E <- gam(E ~ s(Tleaf, by = Elevation) + Elevation, data = dat_Tleaf)

# --- Prediction grid across Elevation levels ---
pred_grid <- dat_Tleaf %>%
  group_by(Elevation) %>%
  summarise(Tleaf_seq = list(seq(min(Tleaf, na.rm=TRUE),
                                 max(Tleaf, na.rm=TRUE),
                                 length.out = 200)), .groups = "drop") %>%
  tidyr::unnest(Tleaf_seq) %>%
  rename(Tleaf = Tleaf_seq)

# --- Predictions ---
pred_A <- predict(gam_A, newdata = pred_grid, se.fit = TRUE)
pred_E <- predict(gam_E, newdata = pred_grid, se.fit = TRUE)

pred_grid$A_pred <- pred_A$fit
pred_grid$A_se <- pred_A$se.fit
pred_grid$E_pred <- pred_E$fit
pred_grid$E_se <- pred_E$se.fit

# scale E predictions into A units for plotting on same axis
pred_grid <- pred_grid %>%
  group_by(Elevation) %>%
  mutate(E_pred_scaled = scale_E_to_A(E_pred, dat_Tleaf[dat_Tleaf$Elevation == first(Elevation), ]),
         E_se_scaled = E_se / diff(range(dat_Tleaf$E[dat_Tleaf$Elevation == first(Elevation)], na.rm=TRUE)) *
           diff(range(dat_Tleaf$A[dat_Tleaf$Elevation == first(Elevation)], na.rm=TRUE))) %>%
  ungroup()

# --- Plot with facets ---
ggplot(dat_Tleaf, aes(x = Tleaf)) +
  geom_point(aes(y = A, color = "Photosynthesis"), alpha = 0.01, size = 1.5) +
  geom_point(aes(y = scale_E_to_A(E, dat_Tleaf), color = "Transpiration"), alpha = 0.01, size = 1.5) +
  geom_ribbon(data = pred_grid,
              aes(ymin = A_pred - A_se, ymax = A_pred + A_se, group = Elevation),
              fill = "green3", alpha = 0.3) +
  geom_ribbon(data = pred_grid,
              aes(ymin = E_pred_scaled - E_se_scaled, ymax = E_pred_scaled + E_se_scaled, group = Elevation),
              fill = "orange", alpha = 0.3) +
  geom_line(data = pred_grid, aes(y = A_pred, color = "Photosynthesis", group = Elevation), size = 1.2) +
  geom_line(data = pred_grid, aes(y = E_pred_scaled, color = "Transpiration", group = Elevation), size = 1.2) +
  facet_wrap(~Elevation) +
  scale_y_continuous(
    name = expression(A~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    sec.axis = sec_axis(~scale_A_to_E(., dat_Tleaf), 
                        name = expression(E~(mmol~m^{-2}~s^{-1})))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis" = "green3", "Transpiration" = "orange")
  ) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  labs(x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = "right")

# Plot A and E vs. Tleaf by Species ####
# --- Define scaling functions (unchanged) ---
scale_E_to_A <- function(E, dat) {
  rng_A <- range(dat$A, na.rm = TRUE)
  rng_E <- range(dat$E, na.rm = TRUE)
  (E - rng_E[1]) / diff(rng_E) * diff(rng_A) + rng_A[1]
}

scale_A_to_E <- function(A_scaled, dat) {
  rng_A <- range(dat$A, na.rm = TRUE)
  rng_E <- range(dat$E, na.rm = TRUE)
  (A_scaled - rng_A[1]) / diff(rng_A) * diff(rng_E) + rng_E[1]
}

# --- Fit GAMs with Species as grouping factor ---
dat_Tleaf$Species <- factor(dat_Tleaf$Species)
gam_A <- gam(A ~ s(Tleaf, Species, bs = "fs") + Species, data = dat_Tleaf)
gam_E <- gam(E ~ s(Tleaf, Species, bs = "fs") + Species, data = dat_Tleaf)

# --- Prediction grid across Species levels ---
pred_grid <- dat_Tleaf %>%
  group_by(Species) %>%
  summarise(Tleaf_seq = list(seq(min(Tleaf, na.rm=TRUE),
                                 max(Tleaf, na.rm=TRUE),
                                 length.out = 200)), .groups = "drop") %>%
  tidyr::unnest(Tleaf_seq) %>%
  rename(Tleaf = Tleaf_seq)

# --- Predictions ---
pred_A <- predict(gam_A, newdata = pred_grid, se.fit = TRUE)
pred_E <- predict(gam_E, newdata = pred_grid, se.fit = TRUE)

pred_grid$A_pred <- pred_A$fit
pred_grid$A_se <- pred_A$se.fit
pred_grid$E_pred <- pred_E$fit
pred_grid$E_se <- pred_E$se.fit

# scale E predictions into A units for plotting on same axis
pred_grid <- pred_grid %>%
  group_by(Species) %>%
  mutate(E_pred_scaled = scale_E_to_A(E_pred, dat_Tleaf[dat_Tleaf$Species == first(Species), ]),
         E_se_scaled = E_se / diff(range(dat_Tleaf$E[dat_Tleaf$Species == first(Species)], na.rm=TRUE)) *
           diff(range(dat_Tleaf$A[dat_Tleaf$Species == first(Species)], na.rm=TRUE))) %>%
  ungroup()

# --- Plot with facets by Species ---
ggplot(dat_Tleaf, aes(x = Tleaf)) +
  geom_point(aes(y = A, color = "Photosynthesis"), alpha = 0.01, size = 1.5) +
  geom_point(aes(y = scale_E_to_A(E, dat_Tleaf), color = "Transpiration"), alpha = 0.01, size = 1.5) +
  geom_ribbon(data = pred_grid,
              aes(ymin = A_pred - A_se, ymax = A_pred + A_se, group = Species),
              fill = "green3", alpha = 0.3) +
  geom_ribbon(data = pred_grid,
              aes(ymin = E_pred_scaled - E_se_scaled, ymax = E_pred_scaled + E_se_scaled, group = Species),
              fill = "orange", alpha = 0.3) +
  geom_line(data = pred_grid, aes(y = A_pred, color = "Photosynthesis", group = Species), size = 1.2) +
  geom_line(data = pred_grid, aes(y = E_pred_scaled, color = "Transpiration", group = Species), size = 1.2) +
  facet_wrap(~Species) +
  scale_y_continuous(
    name = expression(A~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    sec.axis = sec_axis(~scale_A_to_E(., dat_Tleaf), 
                        name = expression(E~(mmol~m^{-2}~s^{-1})))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis" = "green3", "Transpiration" = "orange")
  ) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  labs(x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = "right")

#
# Statistical test of A vs. E at high temps: ####
AE_long <- at.subset3 %>%
  pivot_longer(cols = c(A, E), names_to = "Variable", values_to = "value")

# Now repeat the test
highT_AE <- subset(AE_long, Tleaf > 35 & Variable %in% c("A", "E"))
lm_AE <- lm(value ~ Variable * Tleaf, data = highT_AE)
anova(lm_AE)
summary(lm_AE)
#Sub-test to confirm E slope is not different than zero:
dat35 <- at.subset3%>%
  filter(Tleaf>35)
summary(lm(E~Tleaf, data=dat35))

# Slope for A:
slope_A <- coef(lm_AE)["Tleaf"]  # -0.733
# Slope for E:
slope_E <- coef(lm_AE)["Tleaf"] + coef(lm_AE)["VariableE:Tleaf"]  # -0.73262 + 0.73231 ≈ -0.0003

#annotations:
slope_A <- -0.7326
p_A <- "<<0.001"

slope_E <- -0.0003
p_E <- "<<0.001"

highT_AE <- highT_AE %>%
  mutate(Variable = recode(Variable,
                           "A" = "Photosynthesis (µmol m⁻² s⁻¹)",
                           "E" = "Transpiration (mol m⁻² s⁻¹)"))

high.AE.plot <- ggplot(highT_AE %>% filter(value < 9), aes(x = Tleaf, y = value, color = Variable)) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Leaf temperature (°C)", y = "Value") +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis (µmol m⁻² s⁻¹)" = "green3", "Transpiration (mol m⁻² s⁻¹)" = "orange")
  ) +
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  annotate("text", x = 36.5, y = 3.75, 
           label = paste0("Slope = ", round(slope_A, 2), "\n", "p ", p_A),
           color = "green3", hjust = 0) +
  annotate("text", x = 35, y = 0.5, 
           label = paste0("Slope ≈ ", round(slope_E, 4), "\n", "p ", p_E),
           color = "orange", hjust = 0) +
  theme_classic() +
  ggforce::facet_zoom(ylim = c(0, 0.01))

# Make two insets instead: 
A.inset <- ggplot(at.subset3 %>% filter(Tleaf > 35), aes(x = Tleaf, y = A)) +
  geom_smooth(method = "lm", se = TRUE, color="green3") +
  labs(
    x = "Leaf temperature (°C)", y = "A (µmol m⁻² s⁻¹)") +
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  annotate("text", x = 36, y = 1, 
           label = paste0("Slope = ", round(slope_A, 2), "\n", "p ", p_A),
           color = "green3", hjust = 0) +
  theme_classic()
A.inset
E.inset <- ggplot(at.subset3 %>% filter(Tleaf > 35), aes(x = Tleaf, y = E)) +
  geom_smooth(method = "lm", se = TRUE, color="orange") +
  labs(
    x = "Leaf temperature (°C)", y = "E (mol m⁻² s⁻¹)") +
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  annotate("text", x = 36, y = 0.001, 
           label = paste0("Slope ≈ ", round(slope_E, 4), "\n", "p ", p_E),
           color = "orange", hjust = 0) + ylim(0,0.005)+
  theme_classic()

# Combine AE vs. Tleaf plot with the high temp AE plot: ####
second.row = ggarrange(A.inset, E.inset, nrow=1, ncol=2, labels=c("B", "C"))
ggarrange(A.E.plot, second.row, nrow=2, ncol=1, labels = c("A",""), heights=c(1,0.8))


#Plot A and E vs. Elevation color####
elev_colors <- c("#FDE725", "#F68D45", "#FF0025", "#E849A1","#CC33CC", "#CC99FF", "#99CCFF", "#6699FF", "#000066")
as.factor(raw.dat$Elevation)
# Slopes for A
slopes_A <- raw.dat %>%
  filter(Tleaf > 35) %>%
  group_by(Elevation) %>%
  group_modify(~ {
    mod <- lm(A ~ Tleaf, data = .x)
    tidy(mod) %>%
      filter(term == "Tleaf")
  }) %>%
  mutate(label = paste0("slope=", round(estimate, 2),
                        ", p=", signif(p.value, 2)),
         x_pos = Inf,   # place on right edge
         y_pos = Inf)   # place on top

# Slopes for E
slopes_E <- raw.dat %>%
  filter(Tleaf > 35) %>%
  group_by(Elevation) %>%
  group_modify(~ {
    mod <- lm(E ~ Tleaf, data = .x)
    tidy(mod) %>%
      filter(term == "Tleaf")
  }) %>%
  mutate(label = paste0("slope=", round(estimate, 4),
                        ", p=", signif(p.value, 2)),
         x_pos = Inf,
         y_pos = Inf)

# Plot A
pA <- ggplot(raw.dat %>% filter(Tleaf > 35),
             aes(x = Tleaf, y = A, color = as.factor(Elevation))) +
  geom_point(alpha = 0.4, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = elev_colors) +
  geom_text(data = slopes_A,
            aes(x = x_pos, y = y_pos, label = label, color = as.factor(Elevation)),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.1, size = 3.5) +
  labs(x = "Leaf temperature (°C)", y = "A (µmol m⁻² s⁻¹)", color = "Elevation") +
  theme_classic() +
  facet_wrap(~Elevation)

# Plot E
pE <- ggplot(raw.dat %>% filter(Tleaf > 35),
             aes(x = Tleaf, y = E, color = as.factor(Elevation))) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = elev_colors) +
  geom_text(data = slopes_E,
            aes(x = x_pos, y = y_pos, label = label, color = as.factor(Elevation)),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.1, size = 3.5) +
  labs(x = "Leaf temperature (°C)", y = "E (mol m⁻² s⁻¹)", color = "Elevation") +
  theme_classic() +
  facet_wrap(~Elevation)

pA
pE

# A and then E by Tleaf for each Species ####


###############    START HERE TOMORROW !! I had to go to 30+ instead of 35 figure out if this messes things up.... AND THEN GO BACK TO THE ELEVATION VERSIONS AND MAKE SURE THEY LOOK OK AND ARENT WEIRD...



# --- Define 12 colors for 12 species ---
species_colors <- c(
  "#FDE725", "#F68D45", "#FF0025", "#E849A1",
  "#CC33CC", "#CC99FF", "#99CCFF", "#6699FF",
  "#000066", "#006633", "#66CC66", "#999999"
)

# Make sure Species is a factor
raw.dat$Species <- as.factor(raw.dat$Species)

# --- Slopes for A by species ---
slopes_A <- raw.dat %>%
  filter(Tleaf > 30) %>%
  group_by(Species) %>%
  group_modify(~ {
    mod <- lm(A ~ Tleaf, data = .x)
    tidy(mod) %>% filter(term == "Tleaf")
  }) %>%
  mutate(label = paste0("slope=", round(estimate, 2),
                        ", p=", signif(p.value, 2)),
         x_pos = Inf,   # place on right edge
         y_pos = Inf)   # place on top

# --- Slopes for E by species ---
slopes_E <- raw.dat %>%
  filter(Tleaf > 30) %>%
  group_by(Species) %>%
  group_modify(~ {
    mod <- lm(E ~ Tleaf, data = .x)
    tidy(mod) %>% filter(term == "Tleaf")
  }) %>%
  mutate(label = paste0("slope=", round(estimate, 4),
                        ", p=", signif(p.value, 2)),
         x_pos = Inf,
         y_pos = Inf)

# --- Plot A by species ---
pA <- ggplot(raw.dat %>% filter(Tleaf > 30),
             aes(x = Tleaf, y = A, color = Species)) +
  geom_point(alpha = 0.05, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = species_colors) +
  geom_text(data = slopes_A,
            aes(x = x_pos, y = y_pos, label = label, color = Species),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.1, size = 3.5) +
  labs(x = "Leaf temperature (°C)", 
       y = expression(A~(mu*mol~m^{-2}~s^{-1})), 
       color = "Species") +
  theme_classic() +
  facet_wrap(~Species)

# --- Plot E by species ---
pE <- ggplot(raw.dat %>% filter(Tleaf > 30),
             aes(x = Tleaf, y = E, color = Species)) +
  geom_point(alpha = 0.05, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = species_colors) +
  geom_text(data = slopes_E,
            aes(x = x_pos, y = y_pos, label = label, color = Species),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.1, size = 3.5) +
  labs(x = "Leaf temperature (°C)", 
       y = expression(E~(mmol~m^{-2}~s^{-1})), 
       color = "Species") +
  theme_classic() +
  facet_wrap(~Species)

# Print plots
pA
pE

########## Cuticular conductance stuff






# Cuticular conductance analyses ---------------------------------------------
#### plot A and gcw (cuticular conductance)
# Relationship: gtw = gsw + gcw + gbw
# words relationship: total conductance to water vapor = stomatal + cuticular + boundary
gsw.dat <- at.subset3 %>%
  mutate(gcw = gtw - gsw - gbw) %>%
  filter(is.finite(gcw))%>%
  filter(gcw>-5, gcw<0)

# Scaling functions with safety check
scale_gcw_to_A <- function(gcw) {
  rng_A <- range(gsw.dat$A, na.rm = TRUE)
  rng_gcw <- range(gsw.dat$gcw, na.rm = TRUE)
  if (diff(rng_gcw) == 0) return(rep(mean(rng_A, na.rm = TRUE), length(gcw)))
  (gcw - rng_gcw[1]) / diff(rng_gcw) * diff(rng_A) + rng_A[1]
}

scale_A_to_gcw <- function(A_scaled) {
  rng_A <- range(gsw.dat$A, na.rm = TRUE)
  rng_gcw <- range(gsw.dat$gcw, na.rm = TRUE)
  if (diff(rng_A) == 0) return(rep(mean(rng_gcw, na.rm = TRUE), length(A_scaled)))
  (A_scaled - rng_A[1]) / diff(rng_A) * diff(rng_gcw) + rng_gcw[1]
}

# Fit GAMs
gam_A <- gam(A ~ s(Tleaf, k = 10), data = gsw.dat)
gam_gcw <- gam(gcw ~ s(Tleaf, k = 10), data = gsw.dat)

# Prediction grid
pred_grid <- data.frame(
  Tleaf = seq(min(gsw.dat$Tleaf, na.rm=TRUE),
              max(gsw.dat$Tleaf, na.rm=TRUE),
              length.out = 200)
)

pred_grid$A_pred <- predict(gam_A, newdata = pred_grid)
pred_grid$gcw_pred <- predict(gam_gcw, newdata = pred_grid)
pred_grid$gcw_pred_scaled <- scale_gcw_to_A(pred_grid$gcw_pred)

# Plot
ggplot(gsw.dat, aes(x = Tleaf)) +
  # Points
  geom_point(aes(y = A, color = "Photosynthesis"), alpha = 0.01, size = 2) +
  geom_point(aes(y = scale_gcw_to_A(gcw), color = "gcw"), alpha = 0.02, size = 2) +
  # GAM lines
  geom_line(data = pred_grid, aes(y = A_pred, color = "Photosynthesis"), size = 1.2) +
  geom_line(data = pred_grid, aes(y = gcw_pred_scaled, color = "gcw"), size = 1.2) +
  # Axes
  scale_y_continuous(
    name = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    sec.axis = sec_axis(~scale_A_to_gcw(.), 
                        name = expression(g[cw]~(mol~m^{-2}~s^{-1})))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis" = "green3", "gcw" = "orange")
  ) +
  labs(x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

#### Plot conductances
# Gather conductances
gsw_long <- gsw.dat %>%
  select(Tleaf, gsw, gcw, gbw, gtw) %>%
  pivot_longer(cols = c(gsw, gcw, gbw, gtw), names_to = "Conductance", values_to = "value")

# Shift gcw and gbw for visualization
gsw_long_shifted <- gsw_long %>%
  mutate(value_shifted = case_when(
    Conductance == "gcw" ~ value + 3,   # shift gcw up
    Conductance == "gbw" ~ value - 3,   # shift gbw down
    TRUE ~ value                        # keep gsw and gtw as-is
  ))%>%
  filter(value_shifted<1,value_shifted>-1)

# Fit GAMs and predict with SE
gam_fits_shifted <- gsw_long_shifted %>%
  group_by(Conductance) %>%
  group_modify(~{
    gam_fit <- mgcv::gam(value_shifted ~ s(Tleaf, k = 10), data = .x)
    pred_grid <- data.frame(
      Tleaf = seq(min(.x$Tleaf), max(.x$Tleaf), length.out = 200)
    )
    pred <- predict(gam_fit, newdata = pred_grid, se.fit = TRUE)
    pred_grid$value_shifted <- pred$fit
    pred_grid$se <- pred$se.fit
    pred_grid$Conductance <- unique(.x$Conductance)
    pred_grid
  })

ggplot(gsw_long_shifted, aes(x = Tleaf)) +
  geom_point(aes(y = value_shifted, color = Conductance), 
             alpha = 0.005, size = 1.5) +
  geom_ribbon(data = gam_fits_shifted, 
              aes(x = Tleaf, 
                  ymin = value_shifted - 2*se, 
                  ymax = value_shifted + 2*se, 
                  fill = Conductance), 
              alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = gam_fits_shifted, 
            aes(x = Tleaf, y = value_shifted, color = Conductance), 
            size = 1.2, inherit.aes = FALSE) +
  labs(x = "Leaf temperature (°C)", 
       y = "Conductance (shifted for visualization)") +
  scale_color_manual(values = c("gsw" = "dodgerblue", 
                                "gcw" = "orange", 
                                "gtw" = "darkblue", 
                                "gbw" = "red")) +
  scale_fill_manual(values = c("gsw" = "dodgerblue", 
                               "gcw" = "orange", 
                               "gtw" = "darkblue", 
                               "gbw" = "red")) +
  theme_classic() +
  theme(legend.position = "right")


#### Statistical test of conductance behavior:
highT <- gsw_long %>% filter(Tleaf > quantile(Tleaf, 0.75))

lm_high <- lm(value ~ Conductance + Tleaf*Conductance, data = highT)
anova(lm_high)
summary(lm_high)


# Statistical test for decoupling by Elevation ####
df <- dat_Tleaf %>%
  mutate(
    Elevation = factor(Elevation, levels = sort(unique(Elevation))),   # ordered facets/terms
    E_eps = ifelse(E <= 0, NA, E),                                    # avoid /0 and nonpositive E
    D = A / E_eps,
    logD = log(D)
  )

# focus on the hot range where decoupling is hypothesized (adjust as needed)
hot <- df %>% filter(Tleaf >= 35)

# Linear model: log(A/E) ~ Tleaf * Elevation
m_decouple <- lm(logD ~ Tleaf * Elevation, data = hot)
summary(m_decouple)
anova(m_decouple)

# Brief extraction of per-elevation slopes (Tleaf effect within each elevation)
emtrends(m_decouple, ~ Elevation, var = "Tleaf")  # slope of log(A/E) vs Tleaf by elevation

# Now try for Species ####
df <- dat_Tleaf %>%
  mutate(
    Species = factor(Species, levels = sort(unique(Species))),  # ordered facets/terms
    E_eps   = ifelse(E <= 0, NA, E),                            # avoid /0 and nonpositive E
    D       = A / E_eps,
    logD    = log(D)
  )

# focus on the hot range where decoupling is hypothesized (adjust threshold as needed)
hot <- df %>% filter(Tleaf >= 35)

# Linear model: log(A/E) ~ Tleaf * Species
m_decouple_sp <- lm(logD ~ Tleaf * Species, data = hot)

summary(m_decouple_sp)

# ANOVA
anova(m_decouple_sp)

# slope of log(A/E) vs Tleaf by species
emtrends(m_decouple_sp, ~ Species, var = "Tleaf")
