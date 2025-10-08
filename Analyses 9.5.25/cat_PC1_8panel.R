raw.env.data_pca <- raw.env.data_pca %>%
  mutate(PC1_cat = cut(PC1,
                       breaks = quantile(PC1, probs=c(0,1/3,2/3,1), na.rm=TRUE),
                       labels = c("Low","Medium","High"),
                       include.lowest=TRUE))

# --- Prediction helper for three PC1 levels ---
make_pred_PC1_levels <- function(model, response, n=200) {
  species_levels <- unique(model$model$Species)
  placeholder_species <- species_levels[1]
  
  pc1_vals <- model$model$PC1
  qs <- quantile(pc1_vals, probs=c(0,1/3,2/3,1), na.rm=TRUE)
  pc1_means <- c(
    Low    = mean(pc1_vals[pc1_vals >= qs[1] & pc1_vals <= qs[2]]),
    Medium = mean(pc1_vals[pc1_vals >  qs[2] & pc1_vals <= qs[3]]),
    High   = mean(pc1_vals[pc1_vals >  qs[3] & pc1_vals <= qs[4]])
  )
  
  pred_list <- lapply(names(pc1_means), function(lvl) {
    newdat <- expand.grid(
      Tleaf = seq(min(model$model$Tleaf,na.rm=TRUE),
                  max(model$model$Tleaf,na.rm=TRUE), length.out=n),
      PC1 = pc1_means[lvl],
      Species = placeholder_species
    )
    X <- predict(model, newdata=newdat, type="lpmatrix")
    coefs <- coef(model)
    species_cols <- grep("^Species", colnames(X))
    if(length(species_cols)>0) X[,species_cols] <- rowMeans(X[,species_cols,drop=FALSE])
    newdat$fit <- as.vector(X %*% coefs)
    newdat$response <- response
    newdat$PC1_level <- lvl
    newdat
  })
  
  bind_rows(pred_list)
}

# --- Prediction data ---
pred_A   <- make_pred_PC1_levels(gam_mod_A_PC1, "A")
pred_E   <- make_pred_PC1_levels(gam_mod_E_PC1, "E")
pred_gsw <- make_pred_PC1_levels(gam_mod_gsw_PC1, "gsw")
pred_iW  <- make_pred_PC1_levels(gam_mod_iWUE_PC1,"iWUE")

# --- Derivatives (averaging PC1 and Species) ---
deriv_A   <- finite_diff_deriv_PC1(gam_mod_A_PC1)
deriv_E   <- finite_diff_deriv_PC1(gam_mod_E_PC1)
deriv_gsw <- finite_diff_deriv_PC1(gam_mod_gsw_PC1)
deriv_iW  <- finite_diff_deriv_PC1(gam_mod_iWUE_PC1)

# --- Optimal temperature vertical lines ---
find_opt_Tleaf <- function(pred_df) {
  pred_df %>%
    group_by(PC1_level) %>%
    slice_max(order_by = fit, n = 1, with_ties = FALSE) %>%
    select(PC1_level, Tleaf, fit)
}


# --- Plotting function for PC1 categories ---
plot_top_PC1_levels <- function(pred_df, raw_df, response_name, ymax = NULL){
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
    geom_point(data=raw_df, aes(x=Tleaf,y=.data[[response_name]],color=PC1_cat),
               alpha=0.005,size=1) +
    geom_line(data=pred_df, aes(x=Tleaf,y=fit,color=PC1_level),linewidth=1) +
    # vertical dashed lines at max
    geom_vline(data=optima, aes(xintercept=Tleaf, color=PC1_level),
               linetype="dashed", linewidth=0.8) +
    scale_color_manual(values = c("Low"="#0000CD", "Medium"="#A2CD5A", "High"="#CD5B45")) +
    theme_classic() +
    labs(x=expression(Leaf~Temperature~(degree*C)),
         y=y_lab,
         color="PC1 level")
}



# --- Create top-row plots ---
pA <- plot_top_PC1_levels(pred_A, raw.env.data_pca,"A")
pE <- plot_top_PC1_levels(pred_E, raw.env.data_pca,"E")
pG <- plot_top_PC1_levels(pred_gsw, raw.env.data_pca,"gsw")
pI <- plot_top_PC1_levels(pred_iW, raw.env.data_pca,"iWUE", ymax=300)
# --- Derivative plots (from your previous function) ---
dA <- plot_deriv_PC1(deriv_A,"A")
dE <- plot_deriv_PC1(deriv_E,"E")
dG <- plot_deriv_PC1(deriv_gsw,"gsw")
dI <- plot_deriv_PC1(deriv_iW,"iWUE")

# --- Arrange all panels ---
final_plot <- ggarrange(pA,pE,pG,pI,
                        dA,dE,dG,dI,
                        ncol=4,nrow=2,
                        labels=c("A","B","C","D","E","F","G","H"),
                        common.legend=TRUE,legend="right")
final_plot

