raw.env.data_pca <- raw.env.data_pca %>%
  mutate(PC1_group = as.factor(round(PC1, 4))) %>%
  mutate(WUE = A / E)

# --- Prediction helper for all 9 PC1 levels ---
make_pred_PC1_sites <- function(model, response, n = 200) {
  species_levels <- unique(model$model$Species)
  placeholder_species <- species_levels[1]
  
  # Use all unique PC1 values (one per site)
  pc1_unique <- sort(unique(model$model$PC1))
  
  pred_list <- lapply(pc1_unique, function(pc1_val) {
    newdat <- expand.grid(
      Tleaf = seq(min(model$model$Tleaf, na.rm = TRUE),
                  max(model$model$Tleaf, na.rm = TRUE),
                  length.out = n),
      PC1 = pc1_val,
      Species = placeholder_species,
      curveID = model$model$curveID[1]
    )
    
    X <- predict(model, newdata = newdat, type = "lpmatrix")
    coefs <- coef(model)
    species_cols <- grep("^Species", colnames(X))
    if (length(species_cols) > 0)
      X[, species_cols] <- rowMeans(X[, species_cols, drop = FALSE])
    
    newdat$fit <- as.vector(X %*% coefs)
    newdat$response <- response
    newdat$PC1_group <- as.factor(round(pc1_val, 4))
    newdat
  })
  
  bind_rows(pred_list)
}

# --- Prediction data for each response ---
pred_A   <- make_pred_PC1_sites(gam_mod_A_PC1, "A")
pred_E   <- make_pred_PC1_sites(gam_mod_E_PC1, "E")
pred_gsw <- make_pred_PC1_sites(gam_mod_gsw_PC1, "gsw")
pred_iW  <- make_pred_PC1_sites(gam_mod_iWUE_PC1, "iWUE")
pred_eW  <- make_pred_PC1_sites(gam_mod_WUE_PC1, "WUE")  # NEW

# --- Find optimal temperature (Tleaf) for each PC1 level ---
find_opt_Tleaf <- function(pred_df) {
  pred_df %>%
    group_by(PC1_group) %>%
    slice_max(order_by = fit, n = 1, with_ties = FALSE) %>%
    select(PC1_group, Tleaf, fit)
}

# --- Plotting function for PC1 (9 unique lines) ---
plot_top_PC1_sites <- function(pred_df, raw_df, response_name, ymax = NULL){
  y_lab <- switch(
    response_name,
    "A"    = expression(italic(A)~"(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"    = expression(italic(E)~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw"  = expression(italic(g)[sw]~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "iWUE" = expression((iWUE)~"(" * mu * mol ~ mol^-1 * ")"),
    "WUE"  = expression((WUE)~"(" * mu * mol ~ mol^-1 * ")"),
    response_name
  )
  
  if(!is.null(ymax)){
    raw_df <- raw_df %>% filter(.data[[response_name]] < ymax)
    pred_df <- pred_df %>% filter(fit < ymax)
  }
  
  # --- find maxima
  optima <- find_opt_Tleaf(pred_df)
  
  ggplot() +
    geom_point(data = raw_df,
               aes(x = Tleaf, y = .data[[response_name]], color = PC1_group),
               alpha = 0.002, size = 1) +
    geom_line(data = pred_df,
              aes(x = Tleaf, y = fit, color = PC1_group),
              linewidth = 1) +
    #geom_vline(data = optima,
    #           aes(xintercept = Tleaf, color = PC1_group),
    #           linetype = "dashed", linewidth = 0.8) +
    scale_color_viridis_d(option = "viridis", name = "PC1 score") +
    theme_classic(base_size = 24) +
    labs(x = expression(Leaf~temperature~(degree*C)),
         y = y_lab)
}

# --- Create top-row plots ---
pA <- plot_top_PC1_sites(pred_A, raw.env.data_pca, "A")
pE <- plot_top_PC1_sites(pred_E, raw.env.data_pca, "E")
pG <- plot_top_PC1_sites(pred_gsw, raw.env.data_pca, "gsw")
pI <- plot_top_PC1_sites(pred_iW, raw.env.data_pca, "iWUE", ymax = 300)
pEw <- plot_top_PC1_sites(pred_eW, raw.env.data_pca, "WUE", ymax = 10000)

# =========================================================
# --- Derivative Plots WITH 95% Confidence Bands
# =========================================================

deriv_with_ci_overall <- function(model, n = 200, alpha = 0.05) {
  
  species_levels <- unique(model$model$Species)
  placeholder_species <- species_levels[1]
  
  Tseq <- seq(
    min(model$model$Tleaf, na.rm = TRUE),
    max(model$model$Tleaf, na.rm = TRUE),
    length.out = n
  )
  
  newdat <- expand.grid(
    Tleaf  = Tseq,
    PC1    = mean(model$model$PC1, na.rm = TRUE),
    Species = placeholder_species,
    curveID = model$model$curveID[1]
  )
  
  eps <- diff(range(Tseq)) / n
  
  newdat_hi <- newdat
  newdat_lo <- newdat
  newdat_hi$Tleaf <- newdat$Tleaf + eps
  newdat_lo$Tleaf <- newdat$Tleaf - eps
  
  X_hi <- predict(model, newdat_hi, type = "lpmatrix")
  X_lo <- predict(model, newdat_lo, type = "lpmatrix")
  
  # Finite-difference derivative of the linear predictor
  Xd <- (X_hi - X_lo) / (2 * eps)
  
  # Average Species columns if present (matches your prediction logic)
  species_cols <- grep("^Species", colnames(Xd))
  if (length(species_cols) > 0) {
    Xd[, species_cols] <- rowMeans(Xd[, species_cols, drop = FALSE])
  }
  
  beta <- coef(model)
  Vb   <- vcov(model)
  
  deriv <- as.vector(Xd %*% beta)
  se    <- sqrt(rowSums((Xd %*% Vb) * Xd))
  
  crit <- qnorm(1 - alpha / 2)
  
  data.frame(
    Tleaf = Tseq,
    deriv = deriv,
    se    = se,
    lower = deriv - crit * se,
    upper = deriv + crit * se
  )
}

plot_deriv_overall <- function(deriv_df, response_name){
  
  y_lab <- switch(
    response_name,
    "A"    = bquote(frac(d*italic(A), d*italic(T)[leaf])),
    "E"    = bquote(frac(d*italic(E), d*italic(T)[leaf])),
    "gsw"  = bquote(frac(d*italic(g)[sw], d*italic(T)[leaf])),
    "iWUE" = bquote(frac(d*(iWUE), d*italic(T)[leaf])),
    "WUE"  = bquote(frac(d*(WUE), d*italic(T)[leaf])),
    bquote(frac(d*italic(.(response_name)), d*italic(T)[leaf]))
  )
  
  ggplot(deriv_df, aes(x = Tleaf)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "grey70", alpha = 0.4
    ) +
    geom_line(
      aes(y = deriv),
      color = "black", linewidth = 1.2
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic(base_size = 24) +
    labs(
      x = expression(Leaf~temperature~(degree*C)),
      y = y_lab
    )
}

# --- Compute derivatives with CIs for each response ---
deriv_A   <- deriv_with_ci_overall(gam_mod_A_PC1)
deriv_E   <- deriv_with_ci_overall(gam_mod_E_PC1)
deriv_gsw <- deriv_with_ci_overall(gam_mod_gsw_PC1)
deriv_iW  <- deriv_with_ci_overall(gam_mod_iWUE_PC1)
deriv_eW  <- deriv_with_ci_overall(gam_mod_WUE_PC1)

# --- Create derivative plots ---
dA  <- plot_deriv_overall(deriv_A, "A")
dE  <- plot_deriv_overall(deriv_E, "E")
dG  <- plot_deriv_overall(deriv_gsw, "gsw")
dI  <- plot_deriv_overall(deriv_iW, "iWUE")
dEw <- plot_deriv_overall(deriv_eW, "WUE")

# =========================================================
# --- Combine All Plots (5 on top, 5 derivatives below)
# =========================================================

final_plot <- ggarrange(
  ggarrange(
    pA, pE, pG, pEw, pI,
    ncol = 3, nrow = 2,
    labels = c("A", "B", "C", "D", "E"),
    font.label = list(size = 28, face = "bold"),
    common.legend = TRUE,
    legend = "right"
  ),
  ggarrange(
    dA, dE, dG, dEw, dI,
    ncol = 3, nrow = 2,
    labels = c("F", "G", "H", "I", "J"),
    font.label = list(size = 28, face = "bold")
  ),
  ncol = 1, nrow = 2,
  heights = c(1, 1)
)
ggsave(
  filename = "PC1_GAM_derivatives.png",
  plot = final_plot,
  width  = 20,
  height = 20,
  units  = "in",
  dpi    = 300,
  limitsize = FALSE
)



A_optima_table <- pred_A %>%
  group_by(PC1_group) %>%
  slice_max(order_by = fit, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(
    Topt_A = Tleaf,
    A_max  = fit
  ) %>%
  arrange(as.numeric(as.character(PC1_group)))
