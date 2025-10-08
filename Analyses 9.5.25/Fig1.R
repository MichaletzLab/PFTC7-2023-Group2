# Four models:
gam_mod.A <- gam(A ~ s(Tleaf, k = 3) +
                              s(Elevation, k = 3) +
                              s(vegetation_height, k=3) +
                              s(mean_T2, k = 3) +
                              s(mean_moist_pct, k = 3) +
                              Species +
                              ti(Tleaf, Elevation, k = 3) +
                              ti(Tleaf, vegetation_height, k = 3) +
                              ti(Tleaf, mean_T2, k = 3) +
                              ti(Tleaf, mean_moist_pct, k = 3),
                            data = raw.env.data,
                            method = "REML"
)
summary(gam_mod.A)
#gam_mod.A <- veg.GAM
gam_mod.E <- gam(E ~ s(Tleaf, k = 3) +
                              s(Elevation, k = 3) +
                              s(vegetation_height, k=3) +
                              s(mean_T2, k = 3) +
                              s(mean_moist_pct, k = 3) +
                              Species +
                              ti(Tleaf, Elevation, k = 3) +
                              ti(Tleaf, vegetation_height, k = 3) +
                              ti(Tleaf, mean_T2, k = 3) +
                              ti(Tleaf, mean_moist_pct, k = 3),
                            data = raw.env.data,
                            method = "REML"
)
summary(gam_mod.E)

gam_mod.gsw <- gam(gsw ~ s(Tleaf, k = 3) +
                                s(Elevation, k = 3) +
                                s(vegetation_height, k=3) +
                                s(mean_T2, k = 3) +
                                s(mean_moist_pct, k = 3) +
                                Species +
                                ti(Tleaf, Elevation, k = 3) +
                                ti(Tleaf, vegetation_height, k = 3) +
                                ti(Tleaf, mean_T2, k = 3) +
                                ti(Tleaf, mean_moist_pct, k = 3),
                              data = raw.env.data,
                              method = "REML"
)
summary(gam_mod.gsw)

# iWUE model (same terms)
gam_mod.iWUE <- gam(iWUE ~ s(Tleaf, k = 3) +
                                 s(Elevation, k = 3) +
                                 s(vegetation_height, k=3) +
                                 s(mean_T2, k = 3) +
                                 s(mean_moist_pct, k = 3) +
                                 Species +
                                 ti(Tleaf, Elevation, k = 3) +
                                 ti(Tleaf, vegetation_height, k = 3) +
                                 ti(Tleaf, mean_T2, k = 3) +
                                 ti(Tleaf, mean_moist_pct, k = 3),
                               data = raw.env.data,
                               method = "REML"
)
summary(gam_mod.iWUE)


library(ggplot2)
library(mgcv)
library(dplyr)
library(patchwork)

# ---- Updated prediction function (average Species intercepts) ----
make_pred_df_avg_species <- function(model, response, xvar = "Tleaf", n = 200) {
  mf_names <- names(model$model)
  preds <- mf_names[mf_names != response]
  if (!xvar %in% preds) stop("xvar not found among predictors in the model.")
  
  # compute means for numeric predictors only
  means <- lapply(preds, function(var) {
    if(is.numeric(model$model[[var]])) mean(model$model[[var]], na.rm = TRUE)
    else NA
  })
  names(means) <- preds
  
  # build newdat: replicate means and vary xvar
  newdat <- as.data.frame(means)[rep(1, n), , drop = FALSE]
  newdat[[xvar]] <- seq(min(model$model[[xvar]], na.rm = TRUE),
                        max(model$model[[xvar]], na.rm = TRUE),
                        length.out = n)
  
  # Add a valid factor for Species (any level)
  newdat$Species <- factor(levels(model$model$Species)[1],
                           levels = levels(model$model$Species))
  
  # predict using lpmatrix
  X <- predict(model, newdata = newdat, type = "lpmatrix")
  coefs <- coef(model)
  
  # average Species intercepts
  species_cols <- grep("^Species", colnames(X))
  if(length(species_cols) > 0) {
    X[, species_cols] <- rowMeans(X[, species_cols, drop = FALSE])
  }
  
  fit <- X %*% coefs
  se <- sqrt(rowSums((X %*% vcov(model)) * X))
  
  newdat$fit <- as.vector(fit)
  newdat$se  <- as.vector(se)
  newdat$response <- response
  newdat
}

# ---- Updated derivative function (average Species intercepts) ----
finite_diff_deriv_avg_species <- function(gam_model, xvar, length.out = 200, alpha = 0.05) {
  mf_names <- names(gam_model$model)
  response <- mf_names[1]
  preds <- mf_names[mf_names != response]
  if (!xvar %in% preds) stop("xvar not found among predictors in the model.")
  
  # compute means for numeric predictors
  means <- lapply(preds, function(v) {
    if(is.numeric(gam_model$model[[v]])) mean(gam_model$model[[v]], na.rm = TRUE) else NA
  })
  names(means) <- preds
  
  # build grid
  grid <- as.data.frame(means)[rep(1, length.out), , drop = FALSE]
  grid[[xvar]] <- seq(min(gam_model$model[[xvar]], na.rm = TRUE),
                      max(gam_model$model[[xvar]], na.rm = TRUE),
                      length.out = length.out)
  
  # Add valid Species factor
  grid$Species <- factor(levels(gam_model$model$Species)[1],
                         levels = levels(gam_model$model$Species))
  
  # lpmatrix
  X <- predict(gam_model, newdata = grid, type = "lpmatrix")
  coefs <- coef(gam_model)
  Vb <- vcov(gam_model)
  
  # average Species intercepts
  species_cols <- grep("^Species", colnames(X))
  if(length(species_cols) > 0) {
    X[, species_cols] <- rowMeans(X[, species_cols, drop = FALSE])
  }
  
  dx <- diff(grid[[xvar]])
  step <- dx[1]
  dX <- diff(X) / step
  deriv <- dX %*% coefs
  se <- sqrt(rowSums((dX %*% Vb) * dX))
  
  crit <- qnorm(1 - alpha/2)
  lower <- deriv - crit * se
  upper <- deriv + crit * se
  
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

# ---- Make predictions + derivatives for all models ----
pred_A   <- make_pred_df_avg_species(gam_mod.A,   "A",   "Tleaf")
pred_E   <- make_pred_df_avg_species(gam_mod.E,   "E",   "Tleaf")
pred_gsw <- make_pred_df_avg_species(gam_mod.gsw, "gsw", "Tleaf")
pred_iW  <- make_pred_df_avg_species(gam_mod.iWUE,"iWUE","Tleaf")

deriv_A   <- finite_diff_deriv_avg_species(gam_mod.A,   "Tleaf")
deriv_E   <- finite_diff_deriv_avg_species(gam_mod.E,   "Tleaf")
deriv_gsw <- finite_diff_deriv_avg_species(gam_mod.gsw, "Tleaf")
deriv_iW  <- finite_diff_deriv_avg_species(gam_mod.iWUE,"Tleaf")

# ---- Helper: find Topt ----
find_xmax <- function(pred_df) {
  pred_df$Tleaf[which.max(pred_df$fit)]
}

# ---- Response plot (top row) ----
plot_response <- function(pred_df, raw_df, response_name) {
  x_max <- find_xmax(pred_df)
  y_lab <- switch(
    response_name,
    "A"   = expression(A ~ "(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"   = expression(E ~ "(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw" = expression(g[sw] ~ "(" * mol ~ m^-2 ~ s^-1 * ")"),
    "iWUE"= expression(iWUE ~ "(" * mu * mol ~ mol^-1 * ")"),
    response_name
  )
  
  ggplot(pred_df, aes(x = Tleaf, y = fit)) +
    geom_point(
      data = raw_df,
      aes(x = Tleaf, y = .data[[response_name]]),
      inherit.aes = FALSE, alpha = 0.02, size = 0.5
    ) +
    geom_line(color = "green4", linewidth = 1) +
    geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se),
                fill = "green", alpha = 0.2) +
    geom_vline(xintercept = x_max, linetype = "dashed", color = "black") +
    theme_classic() +
    labs(x = expression(Leaf~Temperature~(degree*C)), y = y_lab)
}

# ---- Derivative plot (bottom row) ----
plot_deriv <- function(deriv_df, response_name) {
  xvar <- names(deriv_df)[1]
  
  deriv_df <- deriv_df %>%
    mutate(sig = lower * upper > 0,
           dir = ifelse(derivative > 0, "Positive", "Negative"))
  
  y_lab <- switch(
    response_name,
    "A"   = expression(dA/dTleaf~"("*mu*"mol"~m^-2~s^-1~degree*C^-1*")"),
    "E"   = expression(dE/dTleaf~"(" * mol~m^-2~s^-1~degree*C^-1 * ")"),
    "gsw" = expression(dg[sw]/dTleaf~"(" * mol~m^-2~s^-1~degree*C^-1 * ")"),
    "iWUE"= expression(diWUE/dTleaf),
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

# ---- Create all plots ----
pA  <- plot_response(pred_A, raw.env.data, "A")
pE  <- plot_response(pred_E, raw.env.data, "E")
pG  <- plot_response(pred_gsw, raw.env.data, "gsw")
pI  <- plot_response(pred_iW, raw.env.data, "iWUE")

dA  <- plot_deriv(deriv_A, "A")
dE  <- plot_deriv(deriv_E, "E")
dG  <- plot_deriv(deriv_gsw, "gsw")
dI  <- plot_deriv(deriv_iW, "iWUE")

all_plots <- list(pA, pE, pG, pI,
                  dA, dE, dG, dI)

# Arrange in 2 rows, 4 columns, add labels
final_plot <- ggarrange(plotlist = all_plots,
                        ncol = 4, nrow = 2,
                        labels = c("A", "B", "C", "D",
                                   "E", "F", "G", "H"),
                        label.x = 0.05, label.y = 0.95,  # position of labels
                        common.legend = TRUE, legend = "right")

final_plot
