# --- Fit GAMs ---
gam_mod.A <- gam(A ~ s(Tleaf, k = 3) +
                   s(Elevation, k = 3) +
                   s(mean_T2, k = 3) +
                   s(mean_moist_pct, k = 3) +
                   ti(Tleaf, Elevation, k = 3) +
                   ti(Tleaf, mean_T2, k = 3) +
                   ti(Tleaf, mean_moist_pct, k = 3),
                 data = raw.env.data,
                 method = "REML")
summary(gam_mod.A)

gam_mod.E <- gam(E ~ s(Tleaf, k = 3) +
                   s(Elevation, k = 3) +
                   s(mean_T2, k = 3) +
                   s(mean_moist_pct, k = 3) +
                   ti(Tleaf, Elevation, k = 3) +
                   ti(Tleaf, mean_T2, k = 3) +
                   ti(Tleaf, mean_moist_pct, k = 3),
                 data = raw.env.data,
                 method = "REML")

gam_mod.gsw <- gam(gsw ~ s(Tleaf, k = 3) +
                     s(Elevation, k = 3) +
                     s(mean_T2, k = 3) +
                     s(mean_moist_pct, k = 3) +
                     ti(Tleaf, Elevation, k = 3) +
                     ti(Tleaf, mean_T2, k = 3) +
                     ti(Tleaf, mean_moist_pct, k = 3),
                   data = raw.env.data,
                   method = "REML")
make_pred_df <- function(model, response, data) {
  newdat <- data.frame(
    Tleaf = seq(min(data$Tleaf, na.rm=TRUE),
                max(data$Tleaf, na.rm=TRUE), length.out = 200),
    Elevation      = mean(data$Elevation, na.rm=TRUE),
    mean_T2        = mean(data$mean_T2, na.rm=TRUE),
    mean_moist_pct = mean(data$mean_moist_pct, na.rm=TRUE)
  )
  p <- predict(model, newdata = newdat, se.fit = TRUE)
  newdat$fit <- p$fit
  newdat$se  <- p$se.fit
  newdat$response <- response
  newdat
}

finite_diff_deriv <- function(gam_model, xvar, length.out = 200, alpha = 0.05) {
  # Find all predictors used in the GAM
  all_vars <- attr(terms(gam_model), "term.labels")
  
  # Hold all other predictors at their means
  means <- sapply(gam_model$model[all_vars], mean, na.rm = TRUE)
  
  # Build grid
  grid <- as.data.frame(as.list(means))
  grid <- grid[rep(1, length.out), , drop = FALSE]
  grid[[xvar]] <- seq(
    min(gam_model$model[[xvar]], na.rm = TRUE),
    max(gam_model$model[[xvar]], na.rm = TRUE),
    length.out = length.out
  )
  
  # Predict and compute finite differences
  X <- predict(gam_model, newdata = grid, type = "lpmatrix")
  coefs <- coef(gam_model)
  Vb <- vcov(gam_model)
  
  dX <- diff(X) / diff(grid[[xvar]][1:2])
  deriv <- dX %*% coefs
  se <- sqrt(rowSums((dX %*% Vb) * dX))
  
  crit <- qnorm(1 - alpha/2)
  lower <- deriv - crit * se
  upper <- deriv + crit * se
  
  # Use setNames to name the first column dynamically
  out <- data.frame(
    x = grid[[xvar]][-1],
    derivative = as.vector(deriv),
    se = as.vector(se),
    lower = as.vector(lower),
    upper = as.vector(upper)
  )
  names(out)[1] <- xvar
  out
}


deriv_A <- finite_diff_deriv(gam_mod.A, "Tleaf")
deriv_E <- finite_diff_deriv(gam_mod.E, "Tleaf")
deriv_gsw <- finite_diff_deriv(gam_mod.gsw, "Tleaf")

pred_A   <- make_pred_df(gam_mod.A,   "A",   raw.env.data)
pred_E   <- make_pred_df(gam_mod.E,   "E",   raw.env.data)
pred_gsw <- make_pred_df(gam_mod.gsw, "gsw", raw.env.data)


# ---- Helper: find x at max y ----
find_xmax <- function(pred_df) {
  pred_df$Tleaf[which.max(pred_df$fit)]
}

# ---- Function for top-row response plots ----
plot_response <- function(pred_df, raw_df, response_name) {
  # Get x at max y (Topt)
  x_max <- find_xmax(pred_df)
  
  # Standard labels for Licor 6800 variables
  y_lab <- switch(
    response_name,
    "A"   = expression(A ~ "(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"   = expression(E ~ "(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw" = expression(g[sw] ~ "(" * mol ~ m^-2 ~ s^-1 * ")"),
    response_name
  )
  
  ggplot(pred_df, aes(x = Tleaf, y = fit)) +
    geom_point(
      data = raw_df,
      aes(x = Tleaf, y = .data[[response_name]]),
      inherit.aes = FALSE, alpha = 0.02, size = 0.5
    ) +
    geom_line(color = "forestgreen", linewidth = 1) +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                fill = "forestgreen", alpha = 0.2) +
    geom_vline(xintercept = x_max, linetype = "dashed", color = "black") +
    theme_classic() +
    labs(x = expression(Leaf~Temperature~(degree*C)), y = y_lab)
}

# ---- Function for derivative plots ----
plot_deriv <- function(deriv_df, response_name) {
  xvar <- names(deriv_df)[1]
  
  # Determine slope direction for coloring
  deriv_df <- deriv_df %>%
    mutate(
      sig = lower * upper > 0,              # CI excludes 0
      dir = ifelse(derivative > 0, "Positive", "Negative")
    )
  
  y_lab <- switch(
    response_name,
    "A"   = expression(dA/dTleaf~"("*mu*"mol"~m^-2~s^-1~degree*C^-1*")"),
    "E"   = expression(dE/dTleaf~"(" * mol~m^-2~s^-1~degree*C^-1 * ")"),
    "gsw" = expression(dg[sw]/dTleaf~"(" * mol~m^-2~s^-1~degree*C^-1 * ")"),
    paste0("d", response_name, "/dTleaf")
  )
  
  ggplot(deriv_df, aes(x = .data[[xvar]], y = derivative)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.2, fill = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(color = "black") +
    geom_point(
      data = filter(deriv_df, sig),
      aes(color = dir), size = 2
    ) +
    scale_color_manual(
      name = "Slope direction",
      values = c("Negative" = "red3", "Positive" = "lightblue3")
    ) +
    labs(x = expression(Leaf~Temperature~(degree*C)), y = y_lab) +
    theme_classic() +
    theme(legend.position = "right")
}

# ---- Create the plots ----
pA  <- plot_response(pred_A, raw.env.data, "A")
pE  <- plot_response(pred_E, raw.env.data, "E")
pG  <- plot_response(pred_gsw, raw.env.data, "gsw")

dA  <- plot_deriv(deriv_A, "A")
dE  <- plot_deriv(deriv_E, "E")
dG  <- plot_deriv(deriv_gsw, "gsw")

# ---- Combine panels ----
(pA | pE | pG) /
  (dA | dE | dG)
