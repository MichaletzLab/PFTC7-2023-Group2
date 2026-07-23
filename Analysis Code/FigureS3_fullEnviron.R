library(ggpubr)

### --- Prediction helpers --- ###
#
# Figure S3 compares three approaches for characterizing environmental 
# variation. All three arms are fitted without a site term (see note in 
# Models.R), so predictions have to be made explicitly for each site then 
# averaged. We don't just want to predict at the covariate means, because the 
# smooths are nonlinear so a prediction using the mean environment isn't the 
# mean prediction across environments, and each arm would end up evaluated at 
# its own average site. And those means would be taken over measurement points
# rather than sites, which weights sites by how heavily they were sampled (same
# bias that was removed from PCA.R in revisions)
#
# Averaging the lpmatrix over the 9 sites with equal weight puts all three
# arms on the same population, and matches deriv_with_ci_overall() in Figure3.R.

make_pred_site_avg <- function(model, response, env_vars, n = 200) {
  
  placeholder_species <- model$model$Species[1]
  placeholder_curveID <- model$model$curveID[1]
  
  # One row per site. The environmental variables are site level, so their
  # unique combinations recover the nine sites.
  sites <- unique(model$model[, env_vars, drop = FALSE])
  
  Tseq <- seq(min(model$model$Tleaf), max(model$model$Tleaf), length.out = n)
  
  X_list <- lapply(seq_len(nrow(sites)), function(k) {
    nd <- data.frame(Tleaf   = Tseq,
                     sites[rep(k, n), , drop = FALSE],
                     Species = placeholder_species,
                     curveID = placeholder_curveID,
                     row.names = NULL)
    predict(model, newdata = nd, type = "lpmatrix")
  })
  
  X <- Reduce(`+`, X_list) / length(X_list)   # equal weight per site
  
  coefs <- coef(model)
  Vb    <- vcov(model)
  
  species_cols <- grep("^Species", colnames(X))
  if (length(species_cols) > 0)
    X[, species_cols] <- rowMeans(X[, species_cols, drop = FALSE])
  
  curveID_cols <- grep("curveID", colnames(X))
  if (length(curveID_cols) > 0)
    X[, curveID_cols] <- 0
  
  fit <- as.vector(X %*% coefs)
  se  <- sqrt(rowSums((X %*% Vb) * X))
  
  data.frame(Tleaf = Tseq, fit = fit, se = se,
             upper = fit + 1.96 * se, lower = fit - 1.96 * se,
             response = response)
}

make_pred_full <- function(model, response, n = 200)
  make_pred_site_avg(model, response,
                     c("mean_moist_pct", "mean_T1", "mean_T2",
                       "mean_air", "vegetation_height"), n)

make_pred_full_PC <- function(model, response, n = 200)
  make_pred_site_avg(model, response, c("PC1", "PC2", "PC3", "PC4", "PC5"), n)

make_pred_PC1 <- function(model, response, n = 200)
  make_pred_site_avg(model, response, "PC1", n)


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

pred_A_pc1   <- make_pred_PC1(gam_mod_PC1_only_A, "A")
pred_E_pc1   <- make_pred_PC1(gam_mod_PC1_only_E, "E")
pred_gsw_pc1 <- make_pred_PC1(gam_mod_PC1_only_gsw, "gsw")
pred_iW_pc1  <- make_pred_PC1(gam_mod_PC1_only_iWUE, "iWUE")
pred_eW_pc1  <- make_pred_PC1(gam_mod_PC1_only_WUE, "WUE")


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
    coord_cartesian(ylim = if (is.null(ymax)) NULL else c(NA, ymax)) +
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
FigS3 <- ggarrange(pA, pE, pG, pEw, pI,
                   ncol = 3, nrow = 2,
                   labels = c("A","B","C","D","E"),
                   common.legend = TRUE, legend = "right")

FigS3

# Output PNG
ggsave("FigureS3.png", plot = FigS3, width = 13, height = 8, units = "in",
       dpi = 300, bg = "white", limitsize = FALSE)

# Output PDF
ggsave("FigureS3.pdf", plot = FigS3, width = 13, height = 8, device = cairo_pdf)

# Quantify differences
arm_diff <- function(pe, pp, p1, nm, peaked = FALSE) {
  mx  <- max(abs(pe$fit - pp$fit), abs(pe$fit - p1$fit), abs(pp$fit - p1$fit))
  rng <- diff(range(c(pe$fit, pp$fit, p1$fit)))
  if (peaked)
    cat(sprintf("%-5s  peak T: env %.1f, PC %.1f, PC1 %.1f degC | max|diff| %.3g (%.0f%% of range)\n",
                nm, pe$Tleaf[which.max(pe$fit)], pp$Tleaf[which.max(pp$fit)],
                p1$Tleaf[which.max(p1$fit)], mx, 100 * mx / rng))
  else
    cat(sprintf("%-5s  max|diff| %.3g (%.0f%% of range)\n", nm, mx, 100 * mx / rng))
}

arm_diff(pred_A,   pred_A_pc,   pred_A_pc1,   "A",    peaked = TRUE)
arm_diff(pred_iW,  pred_iW_pc,  pred_iW_pc1,  "iWUE", peaked = TRUE)
arm_diff(pred_E,   pred_E_pc,   pred_E_pc1,   "E")
arm_diff(pred_gsw, pred_gsw_pc, pred_gsw_pc1, "gsw")
arm_diff(pred_eW,  pred_eW_pc,  pred_eW_pc1,  "WUE")

pk <- function(d) d$Tleaf[which.max(d$fit)]
a <- c(pk(pred_A),  pk(pred_A_pc),  pk(pred_A_pc1))
i <- c(pk(pred_iW), pk(pred_iW_pc), pk(pred_iW_pc1))
cat(sprintf("A     %s  max diff %.3f\n", paste(sprintf("%.3f", a), collapse = ", "), diff(range(a))))
cat(sprintf("iWUE  %s  max diff %.3f\n", paste(sprintf("%.3f", i), collapse = ", "), diff(range(i))))
cat(sprintf("grid spacing %.3f degC\n", diff(pred_A$Tleaf[1:2])))


