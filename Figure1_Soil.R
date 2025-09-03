cat_colors <- c("Low" = "grey40", "Medium" = "lightblue3", "High" = "#3333CC")
gam_A_T1   <- gam(A   ~ T1_cat + s(Tleaf, by = T1_cat, k = 10), data = dat.T1Temp)
gam_E_T1   <- gam(E   ~ T1_cat + s(Tleaf, by = T1_cat, k = 10), data = dat.T1Temp)
gam_gsw_T1 <- gam(gsw ~ T1_cat + s(Tleaf, by = T1_cat, k = 10), data = dat.T1Temp)

pred_grid_T1 <- dat.T1Temp %>%
  group_by(T1_cat) %>%
  summarise(Tleaf = seq(min(Tleaf, na.rm=TRUE), max(Tleaf, na.rm=TRUE), length.out = 200),
            .groups = "drop")

add_preds <- function(pred_df, gam_obj, resp_name){
  pr <- predict(gam_obj, newdata = pred_df, se.fit = TRUE, type = "response")
  pred_df %>%
    mutate(!!paste0(resp_name,"_fit")   := pr$fit,
           !!paste0(resp_name,"_lower") := pr$fit - pr$se.fit,
           !!paste0(resp_name,"_upper") := pr$fit + pr$se.fit)
}

pred_full_T1 <- pred_grid_T1 %>%
  add_preds(gam_A_T1,   "A") %>%
  add_preds(gam_E_T1,   "E") %>%
  add_preds(gam_gsw_T1, "gsw")

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

pA_full_T1 <- plot_gam3(dat.T1Temp, pred_full_T1, "T1_cat", "A",   "A_fit",   "A_lower",   "A_upper",
                        expression(A~(mu*mol~CO[2]~m^-2~s^-1)), "Soil temperature (T1)")
pE_full_T1 <- plot_gam3(dat.T1Temp, pred_full_T1, "T1_cat", "E",   "E_fit",   "E_lower",   "E_upper",
                        expression(E~(mmol~m^-2~s^-1)),          "Soil temperature (T1)")
pg_full_T1 <- plot_gam3(dat.T1Temp, pred_full_T1, "T1_cat", "gsw", "gsw_fit", "gsw_lower", "gsw_upper",
                        expression(g[s][w]~(mol~m^-2~s^-1)),     "Soil temperature (T1)")

fig1_T1 <- ggarrange(pA_full_T1, pE_full_T1, pg_full_T1,
                     nrow = 1, ncol = 3, common.legend = TRUE, legend = "top",
                     labels = c("A","B","C"))

# HOT subset
hot_T1 <- dat.T1Temp %>% filter(Tleaf > TOPT)

gam_A_hot_T1   <- gam(A   ~ T1_cat + s(Tleaf, by = T1_cat, k = 10), data = hot_T1)
gam_E_hot_T1   <- gam(E   ~ T1_cat + s(Tleaf, by = T1_cat, k = 10), data = hot_T1)
gam_gsw_hot_T1 <- gam(gsw ~ T1_cat + s(Tleaf, by = T1_cat, k = 10), data = hot_T1)

pred_hot_grid_T1 <- hot_T1 %>%
  group_by(T1_cat) %>%
  summarise(Tleaf = seq(min(Tleaf, na.rm=TRUE), max(Tleaf, na.rm=TRUE), length.out = 200),
            .groups = "drop")

pred_hot_T1 <- pred_hot_grid_T1 %>%
  add_preds(gam_A_hot_T1,   "A") %>%
  add_preds(gam_E_hot_T1,   "E") %>%
  add_preds(gam_gsw_hot_T1, "gsw")

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

derA_T1   <- deriv_by_cat(gam_A_hot_T1,   pred_hot_grid_T1, "T1_cat")
derE_T1   <- deriv_by_cat(gam_E_hot_T1,   pred_hot_grid_T1, "T1_cat")
dergsw_T1 <- deriv_by_cat(gam_gsw_hot_T1, pred_hot_grid_T1, "T1_cat")

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

pA_hot_T1  <- plot_hot(hot_T1, pred_hot_T1, "T1_cat", "A",   "A_fit",   "A_lower",   "A_upper",
                       expression(A~(mu*mol~CO[2]~m^-2~s^-1)), "Soil temperature (T1)")
pE_hot_T1  <- plot_hot(hot_T1, pred_hot_T1, "T1_cat", "E",   "E_fit",   "E_lower",   "E_upper",
                       expression(E~(mmol~m^-2~s^-1)),          "Soil temperature (T1)")
pg_hot_T1  <- plot_hot(hot_T1, pred_hot_T1, "T1_cat", "gsw", "gsw_fit", "gsw_lower", "gsw_upper",
                       expression(g[s][w]~(mol~m^-2~s^-1)),     "Soil temperature (T1)")

pA_d1_T1   <- plot_deriv(derA_T1,   "T1_cat",
                         expression(frac(dA,dT[leaf])~~"(" * mu*mol~CO[2]~m^-2~s^-1~degree*C^-1 * ")"),
                         "Soil temperature (T1)")
pE_d1_T1   <- plot_deriv(derE_T1,   "T1_cat",
                         expression(frac(dE,dT[leaf])~~"(" * mmol~m^-2~s^-1~degree*C^-1 * ")"),
                         "Soil temperature (T1)")
pgsw_d1_T1 <- plot_deriv(dergsw_T1, "T1_cat",
                         expression(frac(dg[s][w],dT[leaf])~~"(" * mol~m^-2~s^-1~degree*C^-1 * ")"),
                         "Soil temperature (T1)")

fig1_T1
fig2_T1 <- ggarrange(pA_hot_T1, pE_hot_T1, pg_hot_T1, pA_d1_T1, pE_d1_T1, pgsw_d1_T1,
                     nrow = 2, ncol = 3, common.legend = TRUE, legend = "top",
                     labels = c("A","B","C","D","E","F"))
fig2_T1
