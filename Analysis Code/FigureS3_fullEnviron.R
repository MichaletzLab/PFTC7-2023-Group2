make_pred_full <- function(model, response, n = 200) {
  
  placeholder_species <- model$model$Species[1]
  placeholder_curveID <- 9  # omit random effect for population-level prediction
  
  newdat <- data.frame(
    Tleaf = seq(min(model$model$Tleaf), max(model$model$Tleaf), length.out = n),
    mean_moist_pct    = mean(model$model$mean_moist_pct),
    mean_T1           = mean(model$model$mean_T1),
    mean_T2           = mean(model$model$mean_T2),
    mean_air          = mean(model$model$mean_air),
    vegetation_height = mean(model$model$vegetation_height),
    Species           = factor(placeholder_species, levels = levels(model$model$Species)),
    curveID           = placeholder_curveID
  )
  
  # Predict with SEs
  preds <- predict(model, newdata = newdat, type = "response", se.fit = TRUE, exclude = "s(curveID)")
  
  newdat$fit <- preds$fit
  newdat$se  <- preds$se.fit
  newdat$upper <- newdat$fit + 1.96 * newdat$se
  newdat$lower <- newdat$fit - 1.96 * newdat$se
  newdat$response <- response
  newdat
}

make_pred_full_PC <- function(model, response, n = 200) {
  
  placeholder_species <- model$model$Species[1]
  placeholder_curveID <- 9  # omit random effect for population-level prediction
  
  newdat <- data.frame(
    Tleaf = seq(min(model$model$Tleaf), max(model$model$Tleaf), length.out = n),
    PC1   = mean(model$model$PC1),
    PC2   = mean(model$model$PC2),
    PC3   = mean(model$model$PC3),
    PC4   = mean(model$model$PC4),
    PC5   = mean(model$model$PC5),
    Species   = factor(placeholder_species, levels = levels(model$model$Species)),
    curveID   = placeholder_curveID
  )
  
  # Predict with SEs
  preds <- predict(model, newdata = newdat, type = "response", se.fit = TRUE, exclude = "s(curveID)")
  
  newdat$fit <- preds$fit
  newdat$se  <- preds$se.fit
  newdat$upper <- newdat$fit + 1.96 * newdat$se
  newdat$lower <- newdat$fit - 1.96 * newdat$se
  newdat$response <- response
  newdat
}


pred_A   <- make_pred_full(gam_mod_full_environ_A, "A")
pred_E   <- make_pred_full(gam_mod_full_environ_E, "E")
pred_gsw <- make_pred_full(gam_mod_full_environ_gsw, "gsw")
pred_iW  <- make_pred_full(gam_mod_full_environ_iWUE, "iWUE")
pred_eW <- make_pred_full(gam_mod_full_environ_WUE, "WUE")

pred_A_pc   <- make_pred_full_PC(gam_mod_full_PC_A, "A")
pred_E_pc  <- make_pred_full_PC(gam_mod_full_PC_E, "E")
pred_gsw_pc <- make_pred_full_PC(gam_mod_full_PC_gsw, "gsw")
pred_iW_pc  <- make_pred_full_PC(gam_mod_full_PC_iWUE, "iWUE")
pred_eW_pc <- make_pred_full_PC(gam_mod_full_PC_WUE, "WUE")




plot_environ_compare <- function(pred_env, pred_pc, raw_df, response_name, ymax = NULL) {
  
  y_lab <- switch(
    response_name,
    "A"    = expression(A~"(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"    = expression(E~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw"  = expression(g[sw]~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "iWUE" = expression(iWUE~"(" * mu * mol ~ mol^-1 * ")"),
    "WUE"  = expression(WUE~"(" * mu * mol ~ mol^-1 * ")"),
    response_name
  )
  
  if (!is.null(ymax)) {
    raw_df   <- raw_df   %>% filter(.data[[response_name]] < ymax)
    pred_env <- pred_env %>% filter(fit < ymax)
    pred_pc  <- pred_pc  %>% filter(fit < ymax)
  }
  
  ggplot() +
    geom_point(data = raw_df,
               aes(x = Tleaf, y = .data[[response_name]]),
               alpha = 0.005, size = 1) +
    # Confidence ribbons
    geom_ribbon(data = pred_env, aes(x = Tleaf, ymin = lower, ymax = upper),
                fill = "blue", alpha = 0.4) +
    geom_ribbon(data = pred_pc, aes(x = Tleaf, ymin = lower, ymax = upper),
                fill = "red", alpha = 0.4) +
    # Fitted lines
    geom_line(data = pred_env, aes(x = Tleaf, y = fit, color = "Full Env"), linewidth = 1.2) +
    geom_line(data = pred_pc,  aes(x = Tleaf, y = fit, color = "PC GAM"), linewidth = 1.2, linetype = "dashed") +
    scale_color_manual(values = c("Full Env" = "blue", "PC GAM" = "red")) +
    theme_classic() +
    labs(x = expression(Leaf~Temperature~(degree*C)),
         y = y_lab,
         color = "")
}


# --- Create top-row plots ---
pA  <- plot_environ_compare(pred_A,  pred_A_pc,  raw.env.data_pca, "A")
pE  <- plot_environ_compare(pred_E,  pred_E_pc,  raw.env.data_pca, "E")
pG  <- plot_environ_compare(pred_gsw, pred_gsw_pc, raw.env.data_pca, "gsw")
pI  <- plot_environ_compare(pred_iW, pred_iW_pc, raw.env.data_pca, "iWUE", ymax = 300)
pEw <- plot_environ_compare(pred_eW, pred_eW_pc, raw.env.data_pca, "WUE", ymax = 10000)


ggarrange(pA, pE, pG, pEw, pI, labels = c("A","B","C","D","E"))
