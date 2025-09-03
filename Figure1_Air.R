cat_colors <- c("Low" = "grey40", "Medium" = "lightblue3", "High" = "#3333CC")
gam_A_air   <- gam(A   ~ air_cat + s(Tleaf, by = air_cat, k = 10), data = dat.airTemp)
gam_E_air   <- gam(E   ~ air_cat + s(Tleaf, by = air_cat, k = 10), data = dat.airTemp)
gam_gsw_air <- gam(gsw ~ air_cat + s(Tleaf, by = air_cat, k = 10), data = dat.airTemp)

pred_grid_air <- dat.airTemp %>%
  group_by(air_cat) %>%
  summarise(Tleaf = seq(min(Tleaf, na.rm=TRUE), max(Tleaf, na.rm=TRUE), length.out = 200),
            .groups = "drop")

add_preds <- function(pred_df, gam_obj, resp_name){
  pr <- predict(gam_obj, newdata = pred_df, se.fit = TRUE, type = "response")
  pred_df %>%
    mutate(!!paste0(resp_name,"_fit")   := pr$fit,
           !!paste0(resp_name,"_lower") := pr$fit - pr$se.fit,
           !!paste0(resp_name,"_upper") := pr$fit + pr$se.fit)
}

pred_full_air <- pred_grid_air %>%
  add_preds(gam_A_air,   "A") %>%
  add_preds(gam_E_air,   "E") %>%
  add_preds(gam_gsw_air, "gsw")

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

pA_full_air <- plot_gam3(dat.airTemp, pred_full_air, "air_cat", "A",   "A_fit",   "A_lower",   "A_upper",
                         expression(A~(mu*mol~CO[2]~m^-2~s^-1)), "Air temperature")
pE_full_air <- plot_gam3(dat.airTemp, pred_full_air, "air_cat", "E",   "E_fit",   "E_lower",   "E_upper",
                         expression(E~(mmol~m^-2~s^-1)),          "Air temperature")
pg_full_air <- plot_gam3(dat.airTemp, pred_full_air, "air_cat", "gsw", "gsw_fit", "gsw_lower", "gsw_upper",
                         expression(g[s][w]~(mol~m^-2~s^-1)),     "Air temperature")

fig1_air <- ggarrange(pA_full_air, pE_full_air, pg_full_air,
                      nrow = 1, ncol = 3, common.legend = TRUE, legend = "top",
                      labels = c("A","B","C"))

# HOT subset
hot_air <- dat.airTemp %>% filter(Tleaf > TOPT)

gam_A_hot_air   <- gam(A   ~ air_cat + s(Tleaf, by = air_cat, k = 10), data = hot_air)
gam_E_hot_air   <- gam(E   ~ air_cat + s(Tleaf, by = air_cat, k = 10), data = hot_air)
gam_gsw_hot_air <- gam(gsw ~ air_cat + s(Tleaf, by = air_cat, k = 10), data = hot_air)

pred_hot_grid_air <- hot_air %>%
  group_by(air_cat) %>%
  summarise(Tleaf = seq(min(Tleaf, na.rm=TRUE), max(Tleaf, na.rm=TRUE), length.out = 200),
            .groups = "drop")

pred_hot_air <- pred_hot_grid_air %>%
  add_preds(gam_A_hot_air,   "A") %>%
  add_preds(gam_E_hot_air,   "E") %>%
  add_preds(gam_gsw_hot_air, "gsw")

deriv_by_cat <- function(gam_obj, data_ref, cat_var){
  levs <- levels(data_ref[[cat_var]])
  out <- lapply(levs, function(z){
    nd <- data_ref %>% filter(.data[[cat_var]] == z) %>% select(Tleaf, all_of(cat_var)) %>% distinct() %>% arrange(Tleaf)
    if(nrow(nd) == 0) return(NULL)
    dd <- derivatives(gam_obj, term = paste0("s(Tleaf):", cat_var, z), newdata = nd, partial_match = TRUE)
    dd[[cat_var]] <- z
    dd
  })
  bind_rows(out)
}

derA_air   <- deriv_by_cat(gam_A_hot_air,   pred_hot_grid_air, "air_cat")
derE_air   <- deriv_by_cat(gam_E_hot_air,   pred_hot_grid_air, "air_cat")
dergsw_air <- deriv_by_cat(gam_gsw_hot_air, pred_hot_grid_air, "air_cat")

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

pA_hot_air  <- plot_hot(hot_air, pred_hot_air, "air_cat", "A",   "A_fit",   "A_lower",   "A_upper",
                        expression(A~(mu*mol~CO[2]~m^-2~s^-1)), "Air temperature")
pE_hot_air  <- plot_hot(hot_air, pred_hot_air, "air_cat", "E",   "E_fit",   "E_lower",   "E_upper",
                        expression(E~(mmol~m^-2~s^-1)),          "Air temperature")
pg_hot_air  <- plot_hot(hot_air, pred_hot_air, "air_cat", "gsw", "gsw_fit", "gsw_lower", "gsw_upper",
                        expression(g[s][w]~(mol~m^-2~s^-1)),     "Air temperature")

pA_d1_air   <- plot_deriv(derA_air,   "air_cat",
                          expression(frac(dA,dT[leaf])~~"(" * mu*mol~CO[2]~m^-2~s^-1~degree*C^-1 * ")"),
                          "Air temperature")
pE_d1_air   <- plot_deriv(derE_air,   "air_cat",
                          expression(frac(dE,dT[leaf])~~"(" * mmol~m^-2~s^-1~degree*C^-1 * ")"),
                          "Air temperature")
pgsw_d1_air <- plot_deriv(dergsw_air, "air_cat",
                          expression(frac(dg[s][w],dT[leaf])~~"(" * mol~m^-2~s^-1~degree*C^-1 * ")"),
                          "Air temperature")

fig1_air
fig2_air <- ggarrange(pA_hot_air, pE_hot_air, pg_hot_air, pA_d1_air, pE_d1_air, pgsw_d1_air,
                      nrow = 2, ncol = 3, common.legend = TRUE, legend = "top",
                      labels = c("A","B","C","D","E","F"))
fig2_air
