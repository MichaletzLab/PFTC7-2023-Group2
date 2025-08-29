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
TOPT = mean(weibull.discarded$T_opt)
#iWUE
gam_iWUE <- mgcv::gam(I(A/gsw) ~ s(Tleaf, bs = "cs"), data = raw.dat)
pred_grid <- data.frame(Tleaf = seq(min(raw.dat$Tleaf, na.rm = TRUE),
                                    max(raw.dat$Tleaf, na.rm = TRUE),
                                    length.out = 500))
pred_grid$y_pred <- predict(gam_iWUE, newdata = pred_grid)
max_point <- pred_grid[which.max(pred_grid$y_pred), ] # This is Topt for iWUE

iWUE.Tleaf.plot <- ggplot(raw.dat, aes(x = Tleaf, y = A/gsw, color = as.factor(Elevation))) +
  geom_point(alpha = 0.01) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = TOPT, linetype = "dashed") +
  geom_vline(xintercept = max_point$Tleaf, linetype = "dotted", color = "red", size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              se = TRUE, color = "black", linewidth = 1, alpha = 0.3) +
  labs(
    x = "Leaf temperature (°C)",
    y = expression(iWUE~(mu*mol~CO[2]~m^{-2}~s^{-1}~ "/" ~mol~H[2]*O~m^{-2}~s^{-1})),
    color = "Elevation (masl)"
  ) +
  theme_classic() +
  scale_color_viridis_d(guide = guide_legend(override.aes = list(alpha = 1)))

iWUE.Tleaf.plot


A.Tleaf.plot <- A.Tleaf.plot + theme(legend.position = "none")
gsw.Tleaf.plot <- gsw.Tleaf.plot + theme(legend.position = "none")

a.gsw.plot <- ggarrange(A.Tleaf.plot, gsw.Tleaf.plot, nrow = 2, ncol = 1, labels = c("A","B"))
ggarrange(a.gsw.plot, iWUE.Tleaf.plot, nrow = 1, ncol = 2,common.legend = TRUE,legend = "right", labels = c("", "C"))
