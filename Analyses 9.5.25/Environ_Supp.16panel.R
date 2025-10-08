# ----------------------------- #
# 16-panel GAM predictions plot #
# ----------------------------- #

library(mgcv)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ----------------------------- #
# Helper: split environmental variable into 3 categories
split_env_var <- function(data, var){
  qs <- quantile(data[[var]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  data[[paste0(var,"_cat")]] <- cut(data[[var]],
                                    breaks = qs,
                                    include.lowest = TRUE,
                                    labels = c("Low","Medium","High"))
  data
}

# Split all env vars
raw.env.data <- raw.env.data %>%
  split_env_var("Elevation") %>%
  split_env_var("vegetation_height") %>%
  split_env_var("mean_T2") %>%
  split_env_var("mean_moist_pct")

# ----------------------------- #
# Prediction function for 3 categories per env variable
make_pred_env_levels <- function(model, response, xvar = "Tleaf", envvar, n = 200){
  species_levels <- unique(model$model$Species)
  placeholder_species <- species_levels[1]
  
  env_vals <- model$model[[envvar]]
  qs <- quantile(env_vals, probs = c(0,1/3,2/3,1), na.rm = TRUE)
  env_means <- c(
    Low    = mean(env_vals[env_vals >= qs[1] & env_vals <= qs[2]]),
    Medium = mean(env_vals[env_vals >  qs[2] & env_vals <= qs[3]]),
    High   = mean(env_vals[env_vals >  qs[3] & env_vals <= qs[4]])
  )
  
  pred_list <- lapply(names(env_means), function(lvl){
    newdat <- expand.grid(
      Tleaf = seq(min(model$model$Tleaf, na.rm=TRUE),
                  max(model$model$Tleaf, na.rm=TRUE), length.out = n),
      Species = placeholder_species
    )
    
    # Set numeric predictors to mean, except envvar
    num_preds <- names(model$model)[sapply(model$model, is.numeric) & names(model$model) != "Tleaf"]
    for(np in num_preds){
      if(np == envvar) newdat[[np]] <- env_means[lvl]
      else newdat[[np]] <- mean(model$model[[np]], na.rm = TRUE)
    }
    
    # lpmatrix prediction
    X <- predict(model, newdata = newdat, type="lpmatrix")
    coefs <- coef(model)
    species_cols <- grep("^Species", colnames(X))
    if(length(species_cols) > 0) X[, species_cols] <- rowMeans(X[, species_cols, drop=FALSE])
    
    newdat$fit <- as.vector(X %*% coefs)
    newdat$se <- sqrt(rowSums((X %*% vcov(model)) * X))
    newdat$env_level <- lvl
    newdat$response <- response
    newdat$envvar <- envvar
    newdat
  })
  
  bind_rows(pred_list)
}

# ----------------------------- #
# Response plotting function
plot_response_env <- function(pred_df, raw_df, response_name, envvar){
  y_lab <- switch(
    response_name,
    "A"    = expression(A~"(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"    = expression(E~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw"  = expression(g[sw]~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "iWUE" = expression(iWUE~"(" * mu * mol ~ mol^-1 * ")"),
    response_name
  )
  
  ggplot() +
    geom_point(data=raw_df, 
               aes(x=Tleaf, y=.data[[response_name]], color=.data[[paste0(envvar,"_cat")]]),
               alpha=0.005, size=1) +
    geom_line(data=pred_df, aes(x=Tleaf, y=fit, color=env_level), linewidth=1) +
    scale_color_manual(values=c("Low"="#0000CD", "Medium"="#A2CD5A", "High"="#CD5B45")) +
    theme_classic() +
    labs(x=expression(Leaf~Temperature~(degree*C)),
         y=y_lab,
         color=paste(envvar,"level"))
}

# ----------------------------- #
# Models: replace these with your fitted GAMs
# gam_mod.A, gam_mod.E, gam_mod.gsw, gam_mod.iWUE
responses <- list(A=gam_mod.A, E=gam_mod.E, gsw=gam_mod.gsw, iWUE=gam_mod.iWUE)
env_vars <- c("Elevation","vegetation_height","mean_T2","mean_moist_pct")

# ----------------------------- #
# Make predictions and plots
pred_list <- list()
plot_list <- list()
for(resp in names(responses)){
  for(env in env_vars){
    pred_df <- make_pred_env_levels(responses[[resp]], resp, "Tleaf", env)
    pred_list[[paste(resp,env,sep="_")]] <- pred_df
    
    plt <- plot_response_env(pred_df, raw.env.data, resp, env)
    plot_list[[paste(resp,env,sep="_")]] <- plt
  }
}

# ----------------------------- #
# ----------------------------- #
# ----------------------------- #
# Arrange plots in columns by env_var WITH labels
col_plots <- list()


for (i in seq_along(env_vars)) {
  env <- env_vars[i]
  
  # Figure out which 4 letters go in this column
  lab_seq <- LETTERS[seq(i, 16, by = 4)]  # rowwise slice, e.g. c("A","E","I","M")
  
  col_plots[[env]] <- ggarrange(
    plot_list[[paste("A", env, sep = "_")]],
    plot_list[[paste("E", env, sep = "_")]],
    plot_list[[paste("gsw", env, sep = "_")]],
    plot_list[[paste("iWUE", env, sep = "_")]],
    ncol = 1, nrow = 4,
    labels = lab_seq,               # give A/E/I/M to column 1, B/F/J/N to col 2, etc
    common.legend = TRUE, legend = "bottom"
  )
}

# Combine the 4 labeled columns, no extra labels
final_plot <- ggarrange(
  col_plots$Elevation,
  col_plots$vegetation_height,
  col_plots$mean_T2,
  col_plots$mean_moist_pct,
  ncol = 4, nrow = 1,
  common.legend = FALSE
)



# Display
final_plot

# Save as PNG
ggsave("16_panel_plot_by_env.png",
       plot = final_plot,
       width = 20, height = 12, units = "in", dpi = 300)
