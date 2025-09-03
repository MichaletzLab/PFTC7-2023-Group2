# --- 1) Fit GAM with separate smooths for each T1_cat ---
gam_iWUE_T1 <- mgcv::gam(I(A/gsw) ~ s(Tleaf, by = T1_cat, bs = "cs") + T1_cat, 
                         data = dat.T1Temp)

# --- 2) Prediction grid ---
pred_grid_T1 <- expand.grid(
  Tleaf = seq(min(dat.T1Temp$Tleaf, na.rm = TRUE),
              max(dat.T1Temp$Tleaf, na.rm = TRUE),
              length.out = 500),
  T1_cat = levels(dat.T1Temp$T1_cat)
)

pred_T1 <- predict(gam_iWUE_T1, newdata = pred_grid_T1, se.fit = TRUE)

pred_grid_T1 <- pred_grid_T1 %>%
  mutate(
    y_pred = pred_T1$fit,
    se = pred_T1$se.fit,
    ymin = y_pred - 1.96 * se,
    ymax = y_pred + 1.96 * se
  )

# --- Topt per bin ---
topt_by_T1 <- pred_grid_T1 %>%
  group_by(T1_cat) %>%
  slice_max(order_by = y_pred, n = 1, with_ties = FALSE) %>%
  select(T1_cat, Tleaf, y_pred)

# --- Overall GAM (global smooth) ---
gam_T1_all <- mgcv::gam(I(A/gsw) ~ s(Tleaf, bs = "cs"), data = dat.T1Temp)
pred_all_T1 <- data.frame(Tleaf = seq(min(dat.T1Temp$Tleaf, na.rm = TRUE),
                                      max(dat.T1Temp$Tleaf, na.rm = TRUE),
                                      length.out = 500))
pred_all_T1$y_pred <- predict(gam_T1_all, newdata = pred_all_T1)
overall_topt_T1 <- pred_all_T1[which.max(pred_all_T1$y_pred), "Tleaf"]

# --- Define colors ---
T1_colors <- c("Low" = "grey30",
               "Medium" = "lightblue3",
               "High" = "#3333CC")

# --- 3) Filter for display ---
dat_T1_plot <- dat.T1Temp %>% filter((A/gsw) < 400)
pred_grid_T1_plot <- pred_grid_T1 %>% filter(y_pred < 400)

# --- 4) Plot ---
iWUE.Tleaf.T1.plot <- ggplot() +
  geom_point(data = dat_T1_plot,
             aes(x = Tleaf, y = A/gsw, color = T1_cat),
             alpha = 0.01) +
  geom_ribbon(data = pred_grid_T1_plot,
              aes(x = Tleaf, ymin = ymin, ymax = ymax, fill = T1_cat),
              alpha = 0.6) +
  geom_line(data = pred_grid_T1_plot,
            aes(x = Tleaf, y = y_pred, color = T1_cat),
            size = 1.2) +
  geom_vline(xintercept = overall_topt_T1, linetype = "dashed", size = 1) +
  geom_vline(data = topt_by_T1,
             aes(xintercept = Tleaf, color = T1_cat),
             linetype = "dotted", size = 1) +
  labs(
    x = "Leaf temperature (°C)",
    y = expression(iWUE~(mu*mol~CO[2]~mol^{-1}~H[2]*O)),
    color = "Soil temperature (T1)",
    fill = "Soil temperature (T1)"
  ) +
  scale_color_manual(values = T1_colors) +
  scale_fill_manual(values = T1_colors) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 400))

iWUE.Tleaf.T1.plot
