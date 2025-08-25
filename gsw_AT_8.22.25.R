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


#### Make Tleaf binned dataset ####
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

#### Make a plot of gsw and A vs. Tleaf ####
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
#### Plot A and E vs Tleaf ####
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

pred_grid$A_pred <- predict(gam_A, newdata = pred_grid)
pred_grid$E_pred <- predict(gam_E, newdata = pred_grid)
pred_grid$E_pred_scaled <- scale_E_to_A(pred_grid$E_pred)

# Plot
ggplot(dat_Tleaf, aes(x = Tleaf)) +
  # Points
  geom_point(aes(y = A, color = "Photosynthesis"), alpha = 0.05, size = 2) +
  geom_point(aes(y = scale_E_to_A(E), color = "Transpiration"), alpha = 0.01, size = 2) +
  # GAM lines
  geom_line(data = pred_grid, aes(y = A_pred, color = "Photosynthesis"), size = 1.2) +
  geom_line(data = pred_grid, aes(y = E_pred_scaled, color = "Transpiration"), size = 1.2) +
  # Axes
  scale_y_continuous(
    name = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    sec.axis = sec_axis(~scale_A_to_E(.), 
                        name = expression(Transpiration~(mmol~m^{-2}~s^{-1})))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis" = "green3", "Transpiration" = "orange")
  ) +
  labs(x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )



#### plot A and gcw (cuticular conductance) ####
# Relationship: gtw = gsw + gcw + gbw
### words relationship: total conductance to water vapor = stomatal + cuticular + boundary
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


#### Plot conductances ####
# Pivot longer to make one column for conductance type
gsw_long <- gsw.dat %>%
  select(Tleaf, gsw, gcw, gtw, gbw) %>%
  pivot_longer(cols = c(gsw, gcw, gtw, gbw),
               names_to = "Conductance",
               values_to = "value")

# Add shifted values for plotting
gsw_long_shifted <- gsw_long %>%
  mutate(value_shifted = case_when(
    Conductance == "gcw" ~ value + 3,   # shift up
    Conductance == "gbw" ~ value - 3,   # shift down
    TRUE ~ value                        # leave gsw & gtw as-is
  ))%>%
  filter(value_shifted < 1, value_shifted >-1)

# Refit GAMs with shifted values
gam_fits_shifted <- gsw_long_shifted %>%
  group_by(Conductance) %>%
  group_modify(~ {
    gam_fit <- gam(value_shifted ~ s(Tleaf, k = 10), data = .x)
    pred_grid <- data.frame(
      Tleaf = seq(min(.x$Tleaf, na.rm = TRUE),
                  max(.x$Tleaf, na.rm = TRUE),
                  length.out = 200)
    )
    pred_grid$value_shifted <- predict(gam_fit, newdata = pred_grid)
    pred_grid$Conductance <- unique(.x$Conductance)
    pred_grid
  })

# Plot shifted
ggplot(gsw_long_shifted, aes(x = Tleaf, y = value_shifted, color = Conductance)) +
  geom_point(alpha = 0.01, size = 1) +
  geom_line(data = gam_fits_shifted, aes(y = value_shifted, color = Conductance), size = 1.2) +
  labs(
    x = "Leaf temperature (°C)",
    y = expression(Conductance~(mol~m^{-2}~s^{-1})~"(shifted for clarity)")
  ) +
  scale_color_manual(
    values = c(
      gsw = "blue",
      gcw = "orange",
      gtw = "black",
      gbw = "purple"
    )
  ) +
  theme_classic() +
  theme(legend.position = "top")
