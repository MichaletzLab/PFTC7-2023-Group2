# Panels A-E: GAM predictions for each site, restricted to leaf temperatures 
# where at least 4 curves from that site were measured (can be changed below).
# Panels F-J: First derivatives averaged across sites, restricted to temperatures
# with common support (the intersection of the nine site windows).
# Models are fitted to all 117 curves in Models.R. The restriction applies to
# prediction and display only - no fitted values change.
# Run after DataPrep.R and PCA.R (will source Models.R if needed).


if (!exists("raw.env.data_pca"))
  stop("raw.env.data_pca not found. Run DataPrep.R then PCA.R first.")

# Load cached models (or run Models.R if absent)
if (file.exists("outputs/gams_fig3.rds")) {
  gams_fig3 <- readRDS("outputs/gams_fig3.rds")
  cat("Loaded cached Figure 3 GAMs.\n")
} else {
  source("Analysis Code/Models.R")
  gams_fig3 <- list(A = gam_mod_A_PC1, E = gam_mod_E_PC1,
                    gsw = gam_mod_gsw_PC1, iWUE = gam_mod_iWUE_PC1,
                    WUE = gam_mod_WUE_PC1)
  dir.create("outputs", showWarnings = FALSE)
  saveRDS(gams_fig3, "outputs/gams_fig3.rds")
}

gam_mod_A_PC1    <- gams_fig3$A
gam_mod_E_PC1    <- gams_fig3$E
gam_mod_gsw_PC1  <- gams_fig3$gsw
gam_mod_iWUE_PC1 <- gams_fig3$iWUE
gam_mod_WUE_PC1  <- gams_fig3$WUE


# libraries
library(dplyr)
library(mgcv)
library(ggplot2)
library(cowplot)

N_CURVES <- 4              # min curves per site for a plotted region
XLIM_F3  <- c(5, 40)       # shared x-axis, all ten panels

raw.env.data_pca <- raw.env.data_pca %>%
  mutate(PC1_group = as.factor(round(PC1, 4)),
         WUE = A / E)

# ---- Support windows --------------------------------------------------------
support_interval <- function(cr, N, step = 0.05) {
  g   <- seq(min(cr$Tmin), max(cr$Tmax), by = step)
  cnt <- vapply(g, function(t) sum(cr$Tmin <= t & cr$Tmax >= t), numeric(1))
  ok  <- cnt >= N
  if (!any(ok)) return(data.frame(lo = NA_real_, hi = NA_real_))
  r      <- rle(ok)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  keep   <- which(r$values)
  if (length(keep) > 1) warning("Non-contiguous support; using longest run.")
  best <- keep[which.max(r$lengths[keep])]
  data.frame(lo = g[starts[best]], hi = g[ends[best]])
}

site_windows <- raw.env.data_pca %>%
  mutate(Site = factor(Elevation)) %>%
  group_by(Site, curveID) %>%
  summarise(Tmin = min(Tleaf, na.rm = TRUE),
            Tmax = max(Tleaf, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site) %>%
  group_modify(~ support_interval(.x, N_CURVES)) %>%
  ungroup()

common_support <- c(max(site_windows$lo), min(site_windows$hi))

cat(sprintf("N = %d curves. Common support: %.2f to %.2f degC\n",
            N_CURVES, common_support[1], common_support[2]))
print(as.data.frame(site_windows), digits = 4)

# ---- Clipped site predictions ----------------------------------------------
make_pred_PC1_sites <- function(model, response, n = 200,
                                windows = site_windows, clip = TRUE) {
  placeholder_species <- unique(model$model$Species)[1]
  placeholder_curve   <- model$model$curveID[1]
  key  <- unique(model$model[, c("PC1", "Site")])
  beta <- coef(model)
  Tall <- range(model$model$Tleaf, na.rm = TRUE)
  
  bind_rows(lapply(seq_len(nrow(key)), function(k) {
    sw  <- windows[as.character(windows$Site) == as.character(key$Site[k]), ]
    rng <- if (clip && nrow(sw) == 1) c(sw$lo, sw$hi) else Tall
    
    newdat <- expand.grid(
      Tleaf   = seq(rng[1], rng[2], length.out = n),
      PC1     = key$PC1[k], Species = placeholder_species,
      curveID = placeholder_curve, Site = key$Site[k])
    
    X <- predict(model, newdata = newdat, type = "lpmatrix")
    sc <- grep("^Species", colnames(X))
    if (length(sc) > 0) X[, sc] <- rowMeans(X[, sc, drop = FALSE])
    cc <- grep("curveID", colnames(X))
    if (length(cc) > 0) X[, cc] <- 0
    
    newdat$fit       <- as.vector(X %*% beta)
    newdat$response  <- response
    newdat$PC1_group <- as.factor(round(key$PC1[k], 4))
    newdat
  }))
}

# ---- Panels A-E -------------------------------------------------------------
plot_top_PC1_sites <- function(pred_df, raw_df, response_name, ymax = NULL) {
  y_lab <- switch(
    response_name,
    "A"    = expression(italic(A)~"(" * mu * mol ~ m^-2 ~ s^-1 * ")"),
    "E"    = expression(italic(E)~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "gsw"  = expression(italic(g[sw])~"(" * mol ~ m^-2 ~ s^-1 * ")"),
    "iWUE" = expression("iWUE"~"(" * mu * mol ~ mol^-1 * ")"),
    "WUE"  = expression("WUE"~"(" * mu * mol ~ mol^-1 * ")"),
    response_name)
  
  ggplot() +
    geom_point(data = raw_df, aes(x = Tleaf, y = .data[[response_name]]),
               color = "grey95", alpha = 0.06, size = 1) +
    geom_line(data = pred_df, aes(x = Tleaf, y = fit, group = PC1_group),
              color = "white", linewidth = 2.2) +
    geom_line(data = pred_df,
              aes(x = Tleaf, y = fit, color = PC1, group = PC1_group),
              linewidth = 1.4) +
    scale_color_viridis_c(option = "cividis", begin = 0.15, end = 1,
                          name = "PC1 score",
                          guide = guide_colourbar(barheight = 20, barwidth = 1.5)) +
    coord_cartesian(xlim = XLIM_F3,
                    ylim = if (is.null(ymax)) NULL else c(NA, ymax)) +
    theme_classic(base_size = 24) +
    theme(aspect.ratio = 1,
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          legend.text  = element_text(size = 16),
          legend.title = element_text(size = 18)) +
    labs(x = expression(Leaf~temperature~(degree*C)), y = y_lab)
}

# ---- Derivatives on common support -----------------------------------------
deriv_with_ci_overall <- function(model, n = 200, alpha = 0.05,
                                  t_range = common_support) {
  placeholder_species <- unique(model$model$Species)[1]
  placeholder_curve   <- model$model$curveID[1]
  key  <- unique(model$model[, c("PC1", "Site")])
  Tseq <- seq(t_range[1], t_range[2], length.out = n)
  eps  <- diff(range(Tseq)) / n
  
  Xd <- Reduce(`+`, lapply(seq_len(nrow(key)), function(k) {
    nd <- expand.grid(Tleaf = Tseq, PC1 = key$PC1[k],
                      Species = placeholder_species,
                      curveID = placeholder_curve, Site = key$Site[k])
    nd_hi <- nd; nd_hi$Tleaf <- nd$Tleaf + eps
    nd_lo <- nd; nd_lo$Tleaf <- nd$Tleaf - eps
    (predict(model, nd_hi, type = "lpmatrix") -
        predict(model, nd_lo, type = "lpmatrix")) / (2 * eps)
  })) / nrow(key)
  
  sc <- grep("^Species", colnames(Xd))
  if (length(sc) > 0) Xd[, sc] <- rowMeans(Xd[, sc, drop = FALSE])
  cc <- grep("curveID", colnames(Xd))
  if (length(cc) > 0) Xd[, cc] <- 0
  
  beta  <- coef(model); Vb <- vcov(model)
  deriv <- as.vector(Xd %*% beta)
  se    <- sqrt(rowSums((Xd %*% Vb) * Xd))
  crit  <- qnorm(1 - alpha / 2)
  data.frame(Tleaf = Tseq, deriv = deriv, se = se,
             lower = deriv - crit * se, upper = deriv + crit * se)
}

# ---- Panels F-J -------------------------------------------------------------
plot_deriv_overall <- function(deriv_df, response_name, shade_range = NULL) {
  y_lab <- switch(
    response_name,
    "A"    = bquote(d*italic(A) / d*italic(T[leaf]) ~ "(" * mu*mol ~ m^-2 ~ s^-1 ~ degree*"C"^-1 * ")"),
    "E"    = bquote(d*italic(E) / d*italic(T[leaf]) ~ "(" * mol ~ m^-2 ~ s^-1 ~ degree*"C"^-1 * ")"),
    "gsw"  = bquote(d*italic(g[sw]) / d*italic(T[leaf]) ~ "(" * mol ~ m^-2 ~ s^-1 ~ degree*"C"^-1 * ")"),
    "iWUE" = bquote(d~"iWUE" / d*italic(T[leaf]) ~ "(" * mu*mol ~ mol^-1 ~ degree*"C"^-1 * ")"),
    "WUE"  = bquote(d~"WUE" / d*italic(T[leaf]) ~ "(" * mu*mol ~ mol^-1 ~ degree*"C"^-1 * ")"),
    bquote(d*italic(.(response_name)) / d*italic(T[leaf])))
  
  p <- ggplot(deriv_df, aes(x = Tleaf))
  if (!is.null(shade_range))
    p <- p + annotate("rect", xmin = shade_range[1], xmax = shade_range[2],
                      ymin = -Inf, ymax = Inf,
                      fill = "#D6C7EA", alpha = 0.5, color = NA)
  p +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.4) +
    geom_line(aes(y = deriv), color = "black", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_cartesian(xlim = XLIM_F3) +
    theme_classic(base_size = 24) +
    theme(aspect.ratio = 1,
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA)) +
    labs(x = expression(Leaf~temperature~(degree*C)), y = y_lab)
}

# ---- Compute ----------------------------------------------------------------
pred_A   <- make_pred_PC1_sites(gam_mod_A_PC1,    "A")
pred_E   <- make_pred_PC1_sites(gam_mod_E_PC1,    "E")
pred_gsw <- make_pred_PC1_sites(gam_mod_gsw_PC1,  "gsw")
pred_iW  <- make_pred_PC1_sites(gam_mod_iWUE_PC1, "iWUE")
pred_eW  <- make_pred_PC1_sites(gam_mod_WUE_PC1,  "WUE")

deriv_A   <- deriv_with_ci_overall(gam_mod_A_PC1)
deriv_E   <- deriv_with_ci_overall(gam_mod_E_PC1)
deriv_gsw <- deriv_with_ci_overall(gam_mod_gsw_PC1)
deriv_iW  <- deriv_with_ci_overall(gam_mod_iWUE_PC1)
deriv_eW  <- deriv_with_ci_overall(gam_mod_WUE_PC1)

# ---- Decoupling window ------------------------------------------------------
find_upper_crossing <- function(deriv_df) {
  df  <- deriv_df[order(deriv_df$Tleaf), ]
  idx <- which(df$upper[-nrow(df)] >= 0 & df$upper[-1] < 0)
  if (length(idx) == 0) { warning("No zero-crossing in upper CI."); return(NA_real_) }
  i  <- idx[1]
  T1 <- df$Tleaf[i];     U1 <- df$upper[i]
  T2 <- df$Tleaf[i + 1]; U2 <- df$upper[i + 1]
  T1 + (0 - U1) * (T2 - T1) / (U2 - U1)
}

Tleaf_lower <- find_upper_crossing(deriv_A)
Tleaf_upper <- common_support[2] # Upper bound set by underlying data support
stopifnot(!is.na(Tleaf_lower), Tleaf_lower < Tleaf_upper)


cat(sprintf("\nDecoupling window: %.2f to %.2f degC\n", Tleaf_lower, Tleaf_upper))
cat(sprintf("dA/dTleaf max upper CI above onset: %.4g  (want < 0)\n",
            max(deriv_A$upper[deriv_A$Tleaf >= Tleaf_lower])))
cat(sprintf("dE/dTleaf min lower CI over window: %.4g  (want > 0)\n",
            min(deriv_E$lower)))
cat(sprintf("dgsw/dTleaf max upper CI: %.4g  (want < 0)\n", max(deriv_gsw$upper)))

for (nm in c("WUE", "iWUE")) {
  d <- if (nm == "WUE") deriv_eW else deriv_iW
  i <- which(diff(sign(d$deriv)) != 0)
  cat(sprintf("%s derivative sign change within window at: %s\n", nm,
              if (length(i) == 0) "none"
              else paste(sprintf("%.2f", d$Tleaf[i]), collapse = ", ")))
}

# ---- Assemble ---------------------------------------------------------------
scale_factor_E <- 1e4
deriv_E_scaled <- deriv_E %>%
  mutate(deriv = deriv * scale_factor_E,
         lower = lower * scale_factor_E,
         upper = upper * scale_factor_E)

pA  <- plot_top_PC1_sites(pred_A,   raw.env.data_pca, "A")
pE  <- plot_top_PC1_sites(pred_E,   raw.env.data_pca, "E")
pG  <- plot_top_PC1_sites(pred_gsw, raw.env.data_pca, "gsw")
pEw <- plot_top_PC1_sites(pred_eW,  raw.env.data_pca, "WUE",  ymax = 10000)
pI  <- plot_top_PC1_sites(pred_iW,  raw.env.data_pca, "iWUE", ymax = 300)

dA  <- plot_deriv_overall(deriv_A, "A", shade_range = c(Tleaf_lower, Tleaf_upper))
dE  <- plot_deriv_overall(deriv_E_scaled, "E",
                          shade_range = c(Tleaf_lower, Tleaf_upper)) +
  labs(y = bquote(d*italic(E) / d*italic(T[leaf]) ~
                    "(" * 10^-4 ~ mol ~ m^-2 ~ s^-1 ~ degree*"C"^-1 * ")"))
dG  <- plot_deriv_overall(deriv_gsw, "gsw")
dEw <- plot_deriv_overall(deriv_eW,  "WUE")
dI  <- plot_deriv_overall(deriv_iW,  "iWUE")

for (nm in c("pA","pE","pG","pEw","dA","dE","dG","dEw"))
  assign(nm, get(nm) + theme(axis.title.x = element_blank()))

legend_pc1 <- get_legend(pA)

panel_grid <- plot_grid(
  pA  + theme(legend.position = "none"), dA,
  pE  + theme(legend.position = "none"), dE,
  pG  + theme(legend.position = "none"), dG,
  pEw + theme(legend.position = "none"), dEw,
  pI  + theme(legend.position = "none"), dI,
  ncol = 2, align = "hv", axis = "tblr",
  labels = c("A","F","B","G","C","H","D","I","E","J"),
  label_size = 28, label_fontface = "bold")

final_plot <- plot_grid(panel_grid, legend_pc1, ncol = 2,
                        rel_widths = c(1, 0.12))

# Output PNG
ggsave("Figure3.png", plot = final_plot, width = 16, height = 28,
       units = "in", dpi = 300, bg = "white", limitsize = FALSE)

# Output PDF
ggsave("Figure3.pdf", plot = final_plot, width = 16, height = 28,
       units = "in", bg = "white", limitsize = FALSE, device = cairo_pdf)

A_optima_table <- pred_A %>%
  group_by(PC1_group) %>%
  slice_max(order_by = fit, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(Topt_A = Tleaf, A_max = fit) %>%
  arrange(as.numeric(as.character(PC1_group)))
print(as.data.frame(A_optima_table), digits = 4)