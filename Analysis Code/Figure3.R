# libraries
library(mgcv)
library(ggpubr)
library(cowplot)

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
  
  optima <- find_opt_Tleaf(pred_df)
  
  ggplot() +
    geom_point(data = raw_df,
               aes(x = Tleaf, y = .data[[response_name]]),
               color = "grey95", alpha = 0.06, size = 1) +
    geom_line(data = pred_df,
              aes(x = Tleaf, y = fit, group = PC1_group),
              color = "white", linewidth = 2.2) +
    geom_line(data = pred_df,
              aes(x = Tleaf, y = fit, color = PC1, group = PC1_group),
              linewidth = 1.4) +
    scale_color_viridis_c(option = "cividis", begin = 0.15, end = 1, 
                          name = "PC1 score", guide = guide_colourbar()) +
    theme_classic(base_size = 24) +
    theme(
      aspect.ratio      = 1,
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA)
    ) +
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

plot_deriv_overall <- function(deriv_df, response_name, shade_range = NULL){
  
  y_lab <- switch(
    response_name,
    "A"    = bquote(frac(d*italic(A), d*italic(T)[leaf])),
    "E"    = bquote(frac(d*italic(E), d*italic(T)[leaf])),
    "gsw"  = bquote(frac(d*italic(g)[sw], d*italic(T)[leaf])),
    "iWUE" = bquote(frac(d*(iWUE), d*italic(T)[leaf])),
    "WUE"  = bquote(frac(d*(WUE), d*italic(T)[leaf])),
    bquote(frac(d*italic(.(response_name)), d*italic(T)[leaf]))
  )
  
  p <- ggplot(deriv_df, aes(x = Tleaf))
  
  # Decoupling-window band: pale warm neutral, drawn first (behind everything)
  if (!is.null(shade_range)) {
    p <- p + annotate("rect",
                      xmin = shade_range[1], xmax = shade_range[2],
                      ymin = -Inf, ymax = Inf,
                      fill = "#D6C7EA", alpha = 0.5, color = NA)
  }
  
  p +
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
    theme(
      aspect.ratio      = 1,
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA)
    ) +
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

# --- Decoupling window ---
# Finds first Tleaf where the upper CI crosses 0
find_upper_crossing <- function(deriv_df) {
  df <- deriv_df[order(deriv_df$Tleaf), ]
  idx <- which(df$upper[-nrow(df)] >= 0 & df$upper[-1] < 0)
  if (length(idx) == 0) {
    warning("No zero-crossing found in upper CI bound.")
    return(NA_real_)
  }
  i <- idx[1]
  T1 <- df$Tleaf[i];     U1 <- df$upper[i]
  T2 <- df$Tleaf[i + 1]; U2 <- df$upper[i + 1]
  T1 + (0 - U1) * (T2 - T1) / (U2 - U1)
}

# Lower bound: onset of significant A decline (dA/dTleaf CI -> entirely negative)
Tleaf_lower <- find_upper_crossing(deriv_A)
# Upper bound: onset of significant E decline (dE/dTleaf CI -> entirely negative)
Tleaf_upper <- find_upper_crossing(deriv_E)

message(sprintf("Decoupling window: %.2f to %.2f deg C", Tleaf_lower, Tleaf_upper))

# --- Rescale panel G for a narrower y-axis to better align with other panels
# Data are scaled by 10^4 (multiplier included in axis title)
scale_factor_E <- 1e4
deriv_E_scaled <- deriv_E %>%
  mutate(deriv = deriv * scale_factor_E,
         lower = lower * scale_factor_E,
         upper = upper * scale_factor_E)

# --- Create derivative plots ---
dA  <- plot_deriv_overall(deriv_A, "A", shade_range = c(Tleaf_lower, Tleaf_upper))
dE  <- plot_deriv_overall(deriv_E_scaled, "E", shade_range = c(Tleaf_lower, Tleaf_upper)) +
  labs(y = bquote(frac(d*italic(E), d*italic(T)[leaf]) ~ "(x" * 10^-4 * ")"))
dG  <- plot_deriv_overall(deriv_gsw, "gsw")
dI  <- plot_deriv_overall(deriv_iW, "iWUE")
dEw <- plot_deriv_overall(deriv_eW, "WUE")

# =========================================================
# --- Combine All Plots (5 on top, 5 derivatives below)
# =========================================================

# Shared x-axis: only show axis title on bottom panels (E and J)
# pI and dI keep their x-axis titles
pA  <- pA  + theme(axis.title.x = element_blank())
pE  <- pE  + theme(axis.title.x = element_blank())
pG  <- pG  + theme(axis.title.x = element_blank())
pEw <- pEw + theme(axis.title.x = element_blank())
dA  <- dA  + theme(axis.title.x = element_blank())
dE  <- dE  + theme(axis.title.x = element_blank())
dG  <- dG  + theme(axis.title.x = element_blank())
dEw <- dEw + theme(axis.title.x = element_blank())

# Extract the shared legend from pA before stripping legends for the grid
legend_pc1 <- get_legend(pA)

pA_nl  <- pA  + theme(legend.position = "none")
pE_nl  <- pE  + theme(legend.position = "none")
pG_nl  <- pG  + theme(legend.position = "none")
pEw_nl <- pEw + theme(legend.position = "none")
pI_nl  <- pI  + theme(legend.position = "none")

# align = "hv", axis = "tblr" equalizes each panel's actual plotting-region
# margins (not just its outer box) across the whole grid -- this is what
# makes every panel's x-axis start/end at the same screen position, so the
# decoupling band renders identically in F and G.
panel_grid <- plot_grid(
  pA_nl, dA,
  pE_nl, dE,
  pG_nl, dG,
  pEw_nl, dEw,
  pI_nl, dI,
  ncol  = 2,
  align = "hv",
  axis  = "tblr",
  labels = c("A", "F", "B", "G", "C", "H", "D", "I", "E", "J"),
  label_size = 28,
  label_fontface = "bold"
)

final_plot <- plot_grid(
  panel_grid, legend_pc1,
  ncol = 2,
  rel_widths = c(1, 0.12)
)

ggsave(
  filename = "Figure3.png",
  plot = final_plot,
  width  = 16,
  height = 28,
  units  = "in",
  dpi    = 300,
  bg     = "white",
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
