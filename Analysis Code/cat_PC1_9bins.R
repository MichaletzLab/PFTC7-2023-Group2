raw.env.data_pca <- raw.env.data_pca %>%
  mutate(PC1_group = as.factor(round(PC1, 4)))  # rounding avoids floating-point issues


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
      Species = placeholder_species
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
    "A"    = expression(A~"(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"    = expression(E~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw"  = expression(g[sw]~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "iWUE" = expression(iWUE~"(" * mu * mol ~ mol^-1 * ")"),
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
               alpha = 0.01, size = 1) +
    geom_line(data = pred_df,
              aes(x = Tleaf, y = fit, color = PC1_group),
              linewidth = 1) +
    geom_vline(data = optima,
               aes(xintercept = Tleaf, color = PC1_group),
               linetype = "dashed", linewidth = 0.8) +
    scale_color_viridis_d(option = "turbo", name = "PC1 Level (Site)") +
    theme_classic() +
    labs(x = expression(Leaf~Temperature~(degree*C)),
         y = y_lab)
}


# --- Create top-row plots ---
pA <- plot_top_PC1_sites(pred_A, raw.env.data_pca, "A")
pE <- plot_top_PC1_sites(pred_E, raw.env.data_pca, "E")
pG <- plot_top_PC1_sites(pred_gsw, raw.env.data_pca, "gsw")
pI <- plot_top_PC1_sites(pred_iW, raw.env.data_pca, "iWUE", ymax = 300)


# --- Arrange all panels ---
bin9Plot <- ggarrange(
  pA, pE, pG, pI,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  common.legend = TRUE,
  legend = "right"
)

bin9Plot


A_optima <- find_opt_Tleaf(pred_A) %>%
  mutate(PC1 = as.numeric(as.character(PC1_group)))

# --- Run regression: Topt (°C) vs. PC1 ---
Aopt_PC1_mod <- lm(Tleaf ~ PC1, data = A_optima) #Tleaf here is the Topt
summary(Aopt_PC1_mod)

# --- Optional: visualize relationship ---

ggplot(A_optima, aes(x = PC1, y = Tleaf)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_classic(base_size = 14) +
  labs(
    x = "PC1",
    y = expression(T[opt]~"(" * degree * C * ")")
  )















###############Delete the above if this works ########################3
raw.env.data_pca <- raw.env.data_pca %>%
  mutate(PC1_group = as.factor(round(PC1, 4)))  # rounding avoids floating-point issues


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
      Species = placeholder_species
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
    "A"    = expression(A~"(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"    = expression(E~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw"  = expression(g[sw]~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "iWUE" = expression(iWUE~"(" * mu * mol ~ mol^-1 * ")"),
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
               alpha = 0.01, size = 1) +
    geom_line(data = pred_df,
              aes(x = Tleaf, y = fit, color = PC1_group),
              linewidth = 1) +
    geom_vline(data = optima,
               aes(xintercept = Tleaf, color = PC1_group),
               linetype = "dashed", linewidth = 0.8) +
    scale_color_viridis_d(option = "turbo", name = "PC1 Level (Site)") +
    theme_classic() +
    labs(x = expression(Leaf~Temperature~(degree*C)),
         y = y_lab)
}


# --- Create top-row plots ---
pA <- plot_top_PC1_sites(pred_A, raw.env.data_pca, "A")
pE <- plot_top_PC1_sites(pred_E, raw.env.data_pca, "E")
pG <- plot_top_PC1_sites(pred_gsw, raw.env.data_pca, "gsw")
pI <- plot_top_PC1_sites(pred_iW, raw.env.data_pca, "iWUE", ymax = 300)


# --- Regression of A optima vs PC1 ---
A_optima <- find_opt_Tleaf(pred_A) %>%
  mutate(PC1 = as.numeric(as.character(PC1_group)))

Aopt_PC1_mod <- lm(Tleaf ~ PC1, data = A_optima)
summary(Aopt_PC1_mod)

# --- Optional: visualize regression ---
ggplot(A_optima, aes(x = PC1, y = Tleaf)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_classic(base_size = 14) +
  labs(x = "PC1", y = expression(T[opt]~"(" * degree * C * ")"))


# =========================================================
# --- ADDITIONAL SECTION: First Derivative Plots (Overall)
# =========================================================

# --- Function to compute overall derivative wrt Tleaf ---
finite_diff_deriv_overall <- function(model, n = 200) {
  species_levels <- unique(model$model$Species)
  placeholder_species <- species_levels[1]
  
  newdat <- expand.grid(
    Tleaf = seq(min(model$model$Tleaf, na.rm = TRUE),
                max(model$model$Tleaf, na.rm = TRUE),
                length.out = n),
    PC1 = mean(model$model$PC1, na.rm = TRUE),
    Species = placeholder_species
  )
  
  eps <- (max(newdat$Tleaf) - min(newdat$Tleaf)) / n
  newdat_hi <- newdat
  newdat_lo <- newdat
  newdat_hi$Tleaf <- newdat$Tleaf + eps
  newdat_lo$Tleaf <- newdat$Tleaf - eps
  
  pred_hi <- predict(model, newdata = newdat_hi, type = "response")
  pred_lo <- predict(model, newdata = newdat_lo, type = "response")
  
  data.frame(
    Tleaf = newdat$Tleaf,
    deriv = (pred_hi - pred_lo) / (2 * eps)
  )
}

# --- Function to plot derivative ---
plot_deriv_overall <- function(deriv_df, response_name){
  y_lab <- bquote(frac(d~.(response_name), d~T[leaf]))
  
  ggplot(deriv_df, aes(x = Tleaf, y = deriv)) +
    geom_line(color = "black", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    labs(x = expression(Leaf~Temperature~(degree*C)), y = y_lab)
}

# --- Compute derivatives for each response ---
deriv_A   <- finite_diff_deriv_overall(gam_mod_A_PC1)
deriv_E   <- finite_diff_deriv_overall(gam_mod_E_PC1)
deriv_gsw <- finite_diff_deriv_overall(gam_mod_gsw_PC1)
deriv_iW  <- finite_diff_deriv_overall(gam_mod_iWUE_PC1)

# --- Create derivative plots ---
dA <- plot_deriv_overall(deriv_A, "A")
dE <- plot_deriv_overall(deriv_E, "E")
dG <- plot_deriv_overall(deriv_gsw, "gsw")
dI <- plot_deriv_overall(deriv_iW, "iWUE")

# --- Arrange all panels: top = 9-line plots, bottom = derivatives ---
final_plot <- ggarrange(
  ggarrange(pA, pE, pG, pI, ncol = 2, nrow = 2,
            labels = c("A", "B", "C", "D"),
            common.legend = TRUE, legend = "right"),
  ggarrange(dA, dE, dG, dI, ncol = 2, nrow = 2,
            labels = c("E", "F", "G", "H")),
  ncol = 1, nrow = 2, heights = c(2, 1)
)

final_plot
