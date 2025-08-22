raw.dat <- read.csv('data/raw.discardHooks_data.csv')
raw.dat <- raw.dat%>%
  filter(gsw<0.6, gsw>-1)%>%
  filter(!is.na(Tair), !is.na(A), !is.na(Species))
#### A vs. Tair by elevation ####
ggplot(raw.dat, aes(x = Tair, y = A, color=Species)) +
  geom_point(alpha = 0.01) +
  labs(x = "Air temperature (°C)",
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

#### gs vs. Tair by elevation ####
ggplot(raw.dat, aes(x = Tair, y = gsw, color=Species)) +
  geom_point(alpha = 0.01) +
  labs(x = "Air temperature (°C)",
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

#### A vs. Tair by species ####
ggplot(raw.dat, aes(x = Tair, y = A, color=as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  labs(x = "Air temperature (°C)",
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

#### gs vs. Tair by species ####
ggplot(raw.dat, aes(x = Tair, y = gsw, color=as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  labs(x = "Air temperature (°C)",
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
  geom_point(alpha = 0.01, size = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_viridis_d(name = "VPD (kPa)\n(n obs)", direction=-1) +
  labs(x = expression(gsw~(mol~m^{-2}~s^{-1})),
       y = expression(Photosynthesis~(mu*mol~CO[2]~m^{-2}~s^{-1}))) +
  theme_classic() +
  theme(legend.position = "right") +
  facet_wrap(~Elevation) #change to Species or Elevation or comment out

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
