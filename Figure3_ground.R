# --- Ensure T2_cat is factor ---
dat.T2Temp <- dat.T2Temp %>%
  filter(!is.na(T2_cat)) %>%
  mutate(T2_cat = factor(T2_cat, levels = c("Low","Medium","High")))

# --- Fit GAMs ---
gam_A   <- gam(A   ~ T2_cat + s(Tleaf, by = T2_cat, k = 10), data = dat.T2Temp)
gam_E   <- gam(E   ~ T2_cat + s(Tleaf, by = T2_cat, k = 10), data = dat.T2Temp)
gam_gsw <- gam(gsw ~ T2_cat + s(Tleaf, by = T2_cat, k = 10), data = dat.T2Temp)

# --- Prediction grid ---
pred_grid <- dat.T2Temp %>%
  group_by(T2_cat) %>%
  summarise(Tleaf = seq(min(Tleaf, na.rm=TRUE),
                        max(Tleaf, na.rm=TRUE), length.out = 100), .groups = "drop")

# Predict with SE
pred_A   <- predict(gam_A, newdata = pred_grid, se.fit = TRUE)
pred_E   <- predict(gam_E, newdata = pred_grid, se.fit = TRUE)
pred_gsw <- predict(gam_gsw, newdata = pred_grid, se.fit = TRUE)

pred_grid <- pred_grid %>%
  mutate(A_fit = pred_A$fit,   A_upper = pred_A$fit + pred_A$se.fit,   A_lower = pred_A$fit - pred_A$se.fit,
         E_fit = pred_E$fit,   E_upper = pred_E$fit + pred_E$se.fit,   E_lower = pred_E$fit - pred_E$se.fit,
         gsw_fit = pred_gsw$fit, gsw_upper = pred_gsw$fit + pred_gsw$se.fit, gsw_lower = pred_gsw$fit - pred_gsw$se.fit)

# Colors
temp_colors <- c("Low" = "lightblue3", "Medium" = "orange3", "High" = "#3333CC")

# Plot function
plot_gam <- function(df, pred_grid, yvar, fitvar, lowervar, uppervar, catvar, ylab){
  ggplot() +
    geom_point(data = df, aes(x = Tleaf, y = .data[[yvar]], color = .data[[catvar]]), alpha = 0.01) +
    geom_ribbon(data = pred_grid, aes(x = Tleaf, ymin = .data[[lowervar]], ymax = .data[[uppervar]], fill = .data[[catvar]]),
                alpha = 0.2, inherit.aes = FALSE) +
    geom_line(data = pred_grid, aes(x = Tleaf, y = .data[[fitvar]], color = .data[[catvar]]), size = 1.2) +
    scale_color_manual(values = temp_colors) +
    scale_fill_manual(values = temp_colors) +
    labs(x = "Leaf Temperature (°C)", y = ylab, color = "Average Ground Temp", fill = "Average Ground Temp") +
    theme_classic()
}

pA    <- plot_gam(dat.T2Temp, pred_grid, "A", "A_fit", "A_lower", "A_upper", "T2_cat", expression(A~(mu*mol~CO[2]~m^-2~s^-1)))
pE    <- plot_gam(dat.T2Temp, pred_grid, "E", "E_fit", "E_lower", "E_upper", "T2_cat", expression(E~(mmol~m^-2~s^-1)))
pgsw  <- plot_gam(dat.T2Temp, pred_grid, "gsw", "gsw_fit", "gsw_lower", "gsw_upper", "T2_cat", expression(g[s][w]~(mol~m^-2~s^-1)))

# --- USO model ---
g1_values <- dat.T2Temp %>% group_by(T2_cat) %>% summarise(g1=fit_g1(cur_data()), .groups="drop")

dat.T2Temp <- dat.T2Temp %>% left_join(g1_values, by="T2_cat") %>%
  mutate(USO_x = 1.6*(1 + g1/sqrt(VPDleaf))*A/Ca)

annot_text <- paste0("g1 Low = ", round(g1_values$g1[g1_values$T2_cat=="Low"],2),
                     "\n g1 Med = ", round(g1_values$g1[g1_values$T2_cat=="Medium"],2),
                     "\n g1 High = ", round(g1_values$g1[g1_values$T2_cat=="High"],2))

pUSO <- ggplot(dat.T2Temp, aes(x = USO_x, y = gsw, color = T2_cat)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method="lm", se=TRUE) +
  scale_color_manual(values = temp_colors, name="Ground Temp") +
  labs(x = expression(1.6*(1+g[1]/sqrt(VPDleaf))*A/C[a]),
       y = expression(g[s][w]~(mol~m^-2~s^-1))) +
  theme_classic() +
  annotate("text", x=min(dat.T2Temp$USO_x, na.rm=TRUE),
           y=max(dat.T2Temp$gsw, na.rm=TRUE),
           label=annot_text, hjust=0, vjust=1, size=4, color="black")


ggarrange(pA, pE, pgsw, pUSO, nrow = 2, ncol = 2, common.legend = TRUE, labels = c("A","B","C","D"))

lm_USO <- lm(gsw ~ USO_x * T2_cat, data = dat.T2Temp)
summary(lm_USO)
