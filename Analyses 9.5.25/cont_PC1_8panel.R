# ---- Generic prediction function (avg Species intercepts) ----
make_pred_PC1_avg_species <- function(model, response, n = 200) {
  # Pick any valid Species level from the model data
  species_levels <- unique(model$model$Species)
  placeholder_species <- species_levels[1]
  
  newdat <- expand.grid(
    Tleaf = seq(min(model$model$Tleaf, na.rm=TRUE),
                max(model$model$Tleaf, na.rm=TRUE), length.out = n),
    PC1 = mean(model$model$PC1, na.rm=TRUE),
    Species = placeholder_species
  )
  
  X <- predict(model, newdata = newdat, type = "lpmatrix")
  coefs <- coef(model)
  species_cols <- grep("^Species", colnames(X))
  
  # Average species columns to remove bias
  if(length(species_cols) > 0) X[, species_cols] <- rowMeans(X[, species_cols, drop=FALSE])
  
  newdat$fit <- as.vector(X %*% coefs)
  newdat$se  <- sqrt(rowSums((X %*% vcov(model)) * X))
  newdat$response <- response
  newdat
}

# ---- Finite-difference derivative function ----
finite_diff_deriv_PC1 <- function(model, xvar = "Tleaf", n = 200, alpha = 0.05) {
  mf_names <- names(model$model)
  response <- mf_names[1]
  preds <- mf_names[mf_names != response]
  
  # Calculate means for all numeric predictors
  means <- lapply(preds, function(v) mean(model$model[[v]], na.rm = TRUE))
  grid <- as.data.frame(means)[rep(1, n), , drop = FALSE]
  grid[[xvar]] <- seq(
    min(model$model[[xvar]], na.rm = TRUE),
    max(model$model[[xvar]], na.rm = TRUE),
    length.out = n
  )
  
  # Ensure PC1 exists in grid
  if (!"PC1" %in% names(grid) && "PC1" %in% names(model$model)) {
    grid$PC1 <- mean(model$model$PC1, na.rm = TRUE)
  }

  # Add Species factor if present
  if (!"Species" %in% names(grid) && "Species" %in% names(model$model)) {
    grid$Species <- factor(levels(model$model$Species)[1],
                           levels = levels(model$model$Species))
  }
  
  # Build prediction matrix
  X <- predict(model, newdata = grid, type = "lpmatrix")
  coefs <- coef(model)
  
  # Average out species dummy variables
  species_cols <- grep("^Species", colnames(X))
  if (length(species_cols) > 0) {
    X[, species_cols] <- rowMeans(X[, species_cols, drop = FALSE])
  }
  
  # Finite difference
  dx <- diff(grid[[xvar]])
  step <- dx[1]
  dX <- diff(X) / step
  deriv <- dX %*% coefs
  se <- sqrt(rowSums((dX %*% vcov(model)) * dX))
  
  crit <- qnorm(1 - alpha / 2)
  data.frame(
    x = grid[[xvar]][-1],
    derivative = as.vector(deriv),
    se = as.vector(se),
    lower = deriv - crit * se,
    upper = deriv + crit * se
  )
}



# ---- Top-row plotting function ----
plot_top_PC1 <- function(pred_df, raw_df, response_name){
  x_max <- find_xmax(pred_df)
  y_lab <- switch(
    response_name,
    "A"    = expression(A ~ "(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"    = expression(E ~ "(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw"  = expression(g[sw] ~ "(" * mol ~ m^-2 ~ s^-1 * ")"),
    "iWUE" = expression(iWUE ~ "(" * mu * mol ~ mol^-1 * ")"),
    response_name
  )
  
  ggplot(raw_df, aes(x = Tleaf, y = .data[[response_name]], color = PC1)) +
    geom_point(alpha=0.005, size=1.5) +
    scale_color_viridis_c(option="C") +
    geom_line(data=pred_df, aes(x=Tleaf, y=fit), color="black", linewidth=1) +
    theme_classic() +
    geom_vline(xintercept = x_max, linetype = "dashed", color = "black") +
    labs(x = expression(Leaf~Temperature~(degree*C)),
         y = y_lab,
         color = "PC1")
}

# ---- Derivative plotting function ----
plot_deriv_PC1 <- function(deriv_df, response_name){
  deriv_df <- deriv_df %>%
    mutate(sig = lower*upper>0,
           dir = ifelse(derivative>0,"Positive","Negative"))
  
  y_lab <- switch(
    response_name,
    "A"    = expression(dA/dTleaf~"("*mu*"mol"~m^-2~s^-1~degree*C^-1*")"),
    "E"    = expression(dE/dTleaf~"(" * mol~m^-2~s^-1~degree*C^-1 * ")"),
    "gsw"  = expression(dg[sw]/dTleaf~"(" * mol~m^-2~s^-1~degree*C^-1 * ")"),
    "iWUE" = expression(iWUE/dTleaf),
    paste0("d", response_name, "/dTleaf")
  )
  
  ggplot(deriv_df, aes(x = x, y = derivative)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill="grey50") +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_line(color="black") +
    geom_point(data=filter(deriv_df, sig), aes(color=dir), size=2) +
    scale_color_manual(values=c("Negative"="red3","Positive"="lightblue3")) +
    labs(x = expression(Leaf~Temperature~(degree*C)),
         y = y_lab) +
    theme_classic() +
    theme(legend.position="right")
}

# ---- Predictions + derivatives ----
pred_A   <- make_pred_PC1_avg_species(gam_mod_A_PC1, "A")
pred_E   <- make_pred_PC1_avg_species(gam_mod_E_PC1, "E")
pred_gsw <- make_pred_PC1_avg_species(gam_mod_gsw_PC1, "gsw")
pred_iW  <- make_pred_PC1_avg_species(gam_mod_iWUE_PC1, "iWUE")

deriv_A   <- finite_diff_deriv_PC1(gam_mod_A_PC1)
deriv_E   <- finite_diff_deriv_PC1(gam_mod_E_PC1)
deriv_gsw <- finite_diff_deriv_PC1(gam_mod_gsw_PC1)
deriv_iW  <- finite_diff_deriv_PC1(gam_mod_iWUE_PC1)

# ---- Make plots ----
pA  <- plot_top_PC1(pred_A, raw.env.data_pca, "A")
pE  <- plot_top_PC1(pred_E, raw.env.data_pca, "E")
pG  <- plot_top_PC1(pred_gsw, raw.env.data_pca, "gsw")
pI  <- plot_top_PC1(pred_iW, raw.env.data_pca, "iWUE")

dA  <- plot_deriv_PC1(deriv_A, "A")
dE  <- plot_deriv_PC1(deriv_E, "E")
dG  <- plot_deriv_PC1(deriv_gsw, "gsw")
dI  <- plot_deriv_PC1(deriv_iW, "iWUE")

# ---- Arrange plots with labels (2 rows, 4 columns) ----
final_plot <- ggarrange(pA, pE, pG, pI,
                        dA, dE, dG, dI,
                        ncol=4, nrow=2,
                        labels = c("A","B","C","D","E","F","G","H"),
                        label.x=0.05, label.y=0.95,
                        common.legend=TRUE, legend="right")

final_plot
