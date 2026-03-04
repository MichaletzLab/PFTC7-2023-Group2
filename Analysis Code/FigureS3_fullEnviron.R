### --- Prediction helpers --- ###

make_pred_full <- function(model, response, n = 200) {
  
  placeholder_species <- model$model$Species[1]
  placeholder_curveID <- model$model$curveID[1]
  
  newdat <- data.frame(
    Tleaf = seq(min(model$model$Tleaf), max(model$model$Tleaf), length.out = n),
    mean_moist_pct    = mean(model$model$mean_moist_pct),
    mean_T1           = mean(model$model$mean_T1),
    mean_T2           = mean(model$model$mean_T2),
    mean_air          = mean(model$model$mean_air),
    vegetation_height = mean(model$model$vegetation_height),
    Species = placeholder_species,
    curveID           = placeholder_curveID
  )
  
  X <- predict(model, newdata = newdat, type = "lpmatrix")
  coefs <- coef(model)
  Vb <- vcov(model)
  
  species_cols <- grep("^Species", colnames(X))
  if (length(species_cols) > 0)
    X[, species_cols] <- rowMeans(X[, species_cols, drop = FALSE])
  
  curveID_cols <- grep("curveID", colnames(X))
  if (length(curveID_cols) > 0)
    X[, curveID_cols] <- 0
  
  newdat$fit   <- as.vector(X %*% coefs)
  newdat$se    <- sqrt(rowSums((X %*% Vb) * X))
  newdat$upper <- newdat$fit + 1.96 * newdat$se
  newdat$lower <- newdat$fit - 1.96 * newdat$se
  newdat$response <- response
  newdat
}


make_pred_full_PC <- function(model, response, n = 200) {
  
  placeholder_species <- model$model$Species[1]
  placeholder_curveID <- model$model$curveID[1]
  
  newdat <- data.frame(
    Tleaf = seq(min(model$model$Tleaf), max(model$model$Tleaf), length.out = n),
    PC1   = mean(model$model$PC1),
    PC2   = mean(model$model$PC2),
    PC3   = mean(model$model$PC3),
    PC4   = mean(model$model$PC4),
    PC5   = mean(model$model$PC5),
    Species = placeholder_species,
    curveID = placeholder_curveID
  )
  
  X <- predict(model, newdata = newdat, type = "lpmatrix")
  coefs <- coef(model)
  Vb <- vcov(model)
  
  species_cols <- grep("^Species", colnames(X))
  if (length(species_cols) > 0)
    X[, species_cols] <- rowMeans(X[, species_cols, drop = FALSE])
  
  curveID_cols <- grep("curveID", colnames(X))
  if (length(curveID_cols) > 0)
    X[, curveID_cols] <- 0
  
  newdat$fit   <- as.vector(X %*% coefs)
  newdat$se    <- sqrt(rowSums((X %*% Vb) * X))
  newdat$upper <- newdat$fit + 1.96 * newdat$se
  newdat$lower <- newdat$fit - 1.96 * newdat$se
  newdat$response <- response
  newdat
}


make_pred_PC1 <- function(model, response, n = 200) {
  
  placeholder_species <- model$model$Species[1]
  placeholder_curveID <- model$model$curveID[1]
  
  newdat <- data.frame(
    Tleaf = seq(min(model$model$Tleaf), max(model$model$Tleaf), length.out = n),
    PC1   = mean(model$model$PC1),
    Species = placeholder_species,
    curveID = placeholder_curveID
  )
  
  X <- predict(model, newdata = newdat, type = "lpmatrix")
  coefs <- coef(model)
  Vb <- vcov(model)
  
  species_cols <- grep("^Species", colnames(X))
  if (length(species_cols) > 0)
    X[, species_cols] <- rowMeans(X[, species_cols, drop = FALSE])
  
  curveID_cols <- grep("curveID", colnames(X))
  if (length(curveID_cols) > 0)
    X[, curveID_cols] <- 0
  
  newdat$fit   <- as.vector(X %*% coefs)
  newdat$se    <- sqrt(rowSums((X %*% Vb) * X))
  newdat$upper <- newdat$fit + 1.96 * newdat$se
  newdat$lower <- newdat$fit - 1.96 * newdat$se
  newdat$response <- response
  newdat
}


### --- Generate prediction data --- ###

pred_A    <- make_pred_full(gam_mod_full_environ_A, "A")
pred_E    <- make_pred_full(gam_mod_full_environ_E, "E")
pred_gsw  <- make_pred_full(gam_mod_full_environ_gsw, "gsw")
pred_iW   <- make_pred_full(gam_mod_full_environ_iWUE, "iWUE")
pred_eW   <- make_pred_full(gam_mod_full_environ_WUE, "WUE")

pred_A_pc   <- make_pred_full_PC(gam_mod_full_PC_A, "A")
pred_E_pc   <- make_pred_full_PC(gam_mod_full_PC_E, "E")
pred_gsw_pc <- make_pred_full_PC(gam_mod_full_PC_gsw, "gsw")
pred_iW_pc  <- make_pred_full_PC(gam_mod_full_PC_iWUE, "iWUE")
pred_eW_pc  <- make_pred_full_PC(gam_mod_full_PC_WUE, "WUE")

pred_A_pc1   <- make_pred_PC1(gam_mod_A_PC1, "A")
pred_E_pc1   <- make_pred_PC1(gam_mod_E_PC1, "E")
pred_gsw_pc1 <- make_pred_PC1(gam_mod_gsw_PC1, "gsw")
pred_iW_pc1  <- make_pred_PC1(gam_mod_iWUE_PC1, "iWUE")
pred_eW_pc1  <- make_pred_PC1(gam_mod_WUE_PC1, "WUE")


### --- Plot function --- ###

plot_environ_compare <- function(pred_env, pred_pc, pred_pc1, raw_df, response_name, ymax = NULL) {
  
  y_lab <- switch(
    response_name,
    "A"    = expression(italic(A)~"(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"    = expression(italic(E)~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw"  = expression(italic(g)[italic(sw)]~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "iWUE" = expression(iWUE~"(" * mu * mol ~ mol^-1 * ")"),
    "WUE"  = expression(WUE~"(" * mu * mol ~ mol^-1 * ")"),
    response_name
  )
  
  if (!is.null(ymax)) {
    raw_df   <- raw_df   %>% filter(.data[[response_name]] < ymax)
    pred_env <- pred_env %>% filter(fit < ymax)
    pred_pc  <- pred_pc  %>% filter(fit < ymax)
    pred_pc1 <- pred_pc1 %>% filter(fit < ymax)
  }
  
  ggplot() +
    geom_point(data = raw_df,
               aes(x = Tleaf, y = .data[[response_name]]),
               alpha = 0.005, size = 1) +
    # Confidence ribbons
    geom_ribbon(data = pred_env, aes(x = Tleaf, ymin = lower, ymax = upper),
                fill = "#D81B60", alpha = 0.3) +
    geom_ribbon(data = pred_pc, aes(x = Tleaf, ymin = lower, ymax = upper),
                fill = "#1E88E5", alpha = 0.3) +
    geom_ribbon(data = pred_pc1, aes(x = Tleaf, ymin = lower, ymax = upper),
                fill = "#FFC107", alpha = 0.3) +
    # Fitted lines
    geom_line(data = pred_env, aes(x = Tleaf, y = fit, color = "Full Env"), linewidth = 1.2) +
    geom_line(data = pred_pc,  aes(x = Tleaf, y = fit, color = "Full PC"), linewidth = 1.2, linetype = "dashed") +
    geom_line(data = pred_pc1, aes(x = Tleaf, y = fit, color = "PC1-only"), linewidth = 1.2, linetype = "dotdash") +
    scale_color_manual(values = c("Full Env" = "#D81B60", "Full PC" = "#1E88E5", "PC1-only" = "#FFC107")) +
    theme_classic() +
    labs(x = expression(Leaf~temperature~(degree*C)),
         y = y_lab,
         color = "")
}


### --- Make plots --- ###

pA  <- plot_environ_compare(pred_A,  pred_A_pc,  pred_A_pc1,  raw.env.data_pca, "A")
pE  <- plot_environ_compare(pred_E,  pred_E_pc,  pred_E_pc1,  raw.env.data_pca, "E")
pG  <- plot_environ_compare(pred_gsw, pred_gsw_pc, pred_gsw_pc1, raw.env.data_pca, "gsw")
pI  <- plot_environ_compare(pred_iW, pred_iW_pc, pred_iW_pc1, raw.env.data_pca, "iWUE", ymax = 300)
pEw <- plot_environ_compare(pred_eW, pred_eW_pc, pred_eW_pc1, raw.env.data_pca, "WUE", ymax = 10000)


### --- Arrange plots --- ###
ggarrange(pA, pE, pG, pEw, pI,
          labels = c("A","B","C","D","E"),
          common.legend = TRUE, legend = "right")
