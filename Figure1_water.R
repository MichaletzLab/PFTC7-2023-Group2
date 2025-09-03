cat_colors <- c("Low" = "grey40", "Medium" = "lightblue3", "High" = "#3333CC")

# --- Safety: ensure factor levels & drop NAs ---
dat_moist <- dat_environ %>%
  filter(!is.na(moist_cat)) %>%
  mutate(moist_cat = factor(moist_cat, levels = c("Low","Medium","High")))


# --- Fit GAMs with factor-smooth interactions ---
gam_A_moist   <- gam(A   ~ moist_cat + s(Tleaf, by = moist_cat, k = 10), data = dat_moist)
gam_E_moist   <- gam(E   ~ moist_cat + s(Tleaf, by = moist_cat, k = 10), data = dat_moist)
gam_gsw_moist <- gam(gsw ~ moist_cat + s(Tleaf, by = moist_cat, k = 10), data = dat_moist)

# --- Prediction grid (full data) ---
pred_grid_moist <- dat_moist %>%
  group_by(moist_cat) %>%
  summarise(Tleaf = seq(min(Tleaf, na.rm=TRUE), max(Tleaf, na.rm=TRUE), length.out = 200),
            .groups = "drop")

# Helper to add fit + CI from a GAM into pred_grid
add_preds <- function(pred_df, gam_obj, resp_name){
  pr <- predict(gam_obj, newdata = pred_df, se.fit = TRUE, type = "response")
  pred_df %>%
    mutate(!!paste0(resp_name,"_fit")   := pr$fit,
           !!paste0(resp_name,"_lower") := pr$fit - pr$se.fit,
           !!paste0(resp_name,"_upper") := pr$fit + pr$se.fit)
}

pred_full_moist <- pred_grid_moist %>%
  add_preds(gam_A_moist,   "A") %>%
  add_preds(gam_E_moist,   "E") %>%
  add_preds(gam_gsw_moist, "gsw")

# --- Plot helper (full data) ---
plot_gam3 <- function(df_points, df_pred, cat_var, yvar, fit, lo, hi, ylab, color_name){
  ggplot() +
    geom_point(data = df_points, aes(x = Tleaf, y = .data[[yvar]], color = .data[[cat_var]]), alpha = 0.01) +
    geom_ribbon(data = df_pred, aes(x = Tleaf, ymin = .data[[lo]], ymax = .data[[hi]], fill = .data[[cat_var]]),
                alpha = 0.2) +
    geom_line(data = df_pred, aes(x = Tleaf, y = .data[[fit]], color = .data[[cat_var]]), size = 1.2) +
    scale_color_manual(values = cat_colors, name = color_name) +
    scale_fill_manual(values = cat_colors,  name = color_name) +
    labs(x = "Leaf temperature (°C)", y = ylab) +
    theme_classic()
}

pA_full_moist <- plot_gam3(dat_moist, pred_full_moist, "moist_cat", "A",   "A_fit",   "A_lower",   "A_upper",
                           expression(A~(mu*mol~CO[2]~m^-2~s^-1)), "Soil moisture")
pE_full_moist <- plot_gam3(dat_moist, pred_full_moist, "moist_cat", "E",   "E_fit",   "E_lower",   "E_upper",
                           expression(E~(mmol~m^-2~s^-1)),          "Soil moisture")
pg_full_moist <- plot_gam3(dat_moist, pred_full_moist, "moist_cat", "gsw", "gsw_fit", "gsw_lower", "gsw_upper",
                           expression(g[s][w]~(mol~m^-2~s^-1)),     "Soil moisture")

fig1_moist <- ggarrange(pA_full_moist, pE_full_moist, pg_full_moist,
                        nrow = 1, ncol = 3, common.legend = TRUE, legend = "top",
                        labels = c("A","B","C"))

# --- HOT subset (> TOPT) & refit for derivatives ---
hot_moist <- dat_moist %>% filter(Tleaf > TOPT)

gam_A_hot_moist   <- gam(A   ~ moist_cat + s(Tleaf, by = moist_cat, k = 10), data = hot_moist)
gam_E_hot_moist   <- gam(E   ~ moist_cat + s(Tleaf, by = moist_cat, k = 10), data = hot_moist)
gam_gsw_hot_moist <- gam(gsw ~ moist_cat + s(Tleaf, by = moist_cat, k = 10), data = hot_moist)

pred_hot_grid <- hot_moist %>%
  group_by(moist_cat) %>%
  summarise(Tleaf = seq(min(Tleaf, na.rm=TRUE), max(Tleaf, na.rm=TRUE), length.out = 200),
            .groups = "drop")

pred_hot_moist <- pred_hot_grid %>%
  add_preds(gam_A_hot_moist,   "A") %>%
  add_preds(gam_E_hot_moist,   "E") %>%
  add_preds(gam_gsw_hot_moist, "gsw")

# Derivatives per category helper
deriv_by_cat <- function(gam_obj, data_ref, cat_var){
  levs <- levels(data_ref[[cat_var]])
  out <- lapply(levs, function(z){
    nd <- data_ref %>%
      filter(.data[[cat_var]] == z) %>%
      select(Tleaf, all_of(cat_var)) %>%
      distinct() %>%
      arrange(Tleaf)
    if(nrow(nd) == 0) return(NULL)
    # compute derivative for the smooth for this level
    dd <- derivatives(gam_obj,
                      term = paste0("s(Tleaf):", cat_var, z),
                      newdata = nd,
                      partial_match = TRUE)
    dd[[cat_var]] <- z
    dd
  })
  bind_rows(out)
}

derA_moist   <- deriv_by_cat(gam_A_hot_moist,   pred_hot_grid, "moist_cat")
derE_moist   <- deriv_by_cat(gam_E_hot_moist,   pred_hot_grid, "moist_cat")
dergsw_moist <- deriv_by_cat(gam_gsw_hot_moist, pred_hot_grid, "moist_cat")

# Plot helpers (hot data + derivatives)
plot_hot <- function(df_points, df_pred, cat_var, yvar, fit, lo, hi, ylab, color_name){
  ggplot() +
    geom_point(data = df_points, aes(x = Tleaf, y = .data[[yvar]], color = .data[[cat_var]]), alpha = 0.02) +
    geom_ribbon(data = df_pred, aes(x = Tleaf, ymin = .data[[lo]], ymax = .data[[hi]], fill = .data[[cat_var]]),
                alpha = 0.2) +
    geom_line(data = df_pred, aes(x = Tleaf, y = .data[[fit]], color = .data[[cat_var]]), size = 1.2) +
    scale_color_manual(values = cat_colors, name = color_name) +
    scale_fill_manual(values = cat_colors,  name = color_name) +
    labs(x = "Leaf temperature (°C)", y = ylab) +
    theme_classic()
}

plot_deriv <- function(deriv_df, cat_var, ylab, color_name){
  ggplot(deriv_df, aes(x = Tleaf, y = .derivative, color = .data[[cat_var]])) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .data[[cat_var]]), alpha = 0.15, color = NA) +
    geom_line(size = 1.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = cat_colors, name = color_name) +
    scale_fill_manual(values = cat_colors,  name = color_name) +
    labs(x = "Leaf temperature (°C)", y = ylab) +
    theme_classic()
}

pA_hot   <- plot_hot(hot_moist, pred_hot_moist, "moist_cat", "A",   "A_fit",   "A_lower",   "A_upper",
                     expression(A~(mu*mol~CO[2]~m^-2~s^-1)), "Soil moisture")
pE_hot   <- plot_hot(hot_moist, pred_hot_moist, "moist_cat", "E",   "E_fit",   "E_lower",   "E_upper",
                     expression(E~(mmol~m^-2~s^-1)),          "Soil moisture")
pg_hot   <- plot_hot(hot_moist, pred_hot_moist, "moist_cat", "gsw", "gsw_fit", "gsw_lower", "gsw_upper",
                     expression(g[s][w]~(mol~m^-2~s^-1)),     "Soil moisture")

pA_d1    <- plot_deriv(derA_moist,   "moist_cat",
                       expression(frac(dA,dT[leaf])~~"(" * mu*mol~CO[2]~m^-2~s^-1~degree*C^-1 * ")"),
                       "Soil moisture")
pE_d1    <- plot_deriv(derE_moist,   "moist_cat",
                       expression(frac(dE,dT[leaf])~~"(" * mmol~m^-2~s^-1~degree*C^-1 * ")"),
                       "Soil moisture")
pgsw_d1  <- plot_deriv(dergsw_moist, "moist_cat",
                       expression(frac(dg[s][w],dT[leaf])~~"(" * mol~m^-2~s^-1~degree*C^-1 * ")"),
                       "Soil moisture")

fig2_moist <- ggarrange(pA_hot, pE_hot, pg_hot, pA_d1, pE_d1, pgsw_d1,
                        nrow = 2, ncol = 3, common.legend = TRUE, legend = "top",
                        labels = c("A","B","C","D","E","F"))

# Show (or save) figures:
fig1_moist
fig2_moist
