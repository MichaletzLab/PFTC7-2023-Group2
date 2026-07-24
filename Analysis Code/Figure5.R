# Figure 5 / Figure S5: temperature sensitivities and reference values for gas
# exchange variables with monotonic temperature responses, vs a principal
# component axis.
#
# The axis and output names are parameterized. Run this file directly for
# Figure 5 (PC1). FigureS5.R sets PC_AXIS and FIG_TAG, then sources this file.

# ---- axis and output configuration -----------------------------------------
if (!exists("PC_AXIS")) PC_AXIS <- "PC1"
if (!exists("FIG_TAG")) FIG_TAG <- "Figure5"
message("Monotonic response figure: axis = ", PC_AXIS,
        ", output tag = ", FIG_TAG)

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(cowplot)

source("Analysis Code/Fig_helpers.R")

# Create WUE = A/E, with Tleaf centered at 25 C (consistent w/ Fig. 4 SS).
# The axis column is renamed to a neutral "PC" so that every formula,
# coefficient name, and aesthetic below is identical for PC1 and PC2.
dat5 <- raw.env.data_pca %>%
  mutate(WUE = A / E, Tc = Tleaf - 25, Site = factor(Elevation),
         PC = .data[[PC_AXIS]])

sp_levels <- species_level_set()   # canonical level set, shared across figures
pc_lim    <- range(dat5$PC, na.rm = TRUE)
pc_seq    <- seq(pc_lim[1], pc_lim[2], length.out = 100)
sp_labs   <- species_labels(sp_levels)

# Define symbols and units
u_flux       <- quote(mol ~ m^-2 ~ s^-1)                   # E, gsw : level
u_flux_slope <- quote(mol ~ m^-2 ~ s^-1 ~ degree*"C"^-1)   # E, gsw : slope
u_wue        <- quote(mu*mol ~ mol^-1)                     # WUE    : level
u_wue_slope  <- quote(mu*mol ~ mol^-1 ~ degree*"C"^-1)     # WUE    : slope

var_sym <- list(E = quote(italic(E)), gsw = quote(italic(g[sw])), WUE = quote(~"WUE"))
sym_level <- list(E   = quote(italic(E[25])),
                  gsw = quote(italic(g["sw,25"])),
                  WUE = quote("WUE"[25]))
var_ulv <- list(E = u_flux,       gsw = u_flux,       WUE = u_wue)       # level units
var_usl <- list(E = u_flux_slope, gsw = u_flux_slope, WUE = u_wue_slope) # slope units

# Fit one variable (full site structure, fall back to site intercept if singular)
# and return per-curve OLS estimates, population trends, p-values, and df.
fit_var <- function(col, data) {
  f_full <- as.formula(paste0(col,
                              " ~ Tc * PC + (1 + Tc || Site) + (1 + Tc || Species/curveID)"))
  mod <- suppressMessages(lmer(f_full, data = data))
  structure_used <- "site intercept + slope"
  if (isSingular(mod)) {
    f_red <- as.formula(paste0(col,
                               " ~ Tc * PC + (1 | Site) + (1 + Tc || Species/curveID)"))
    mod <- suppressMessages(lmer(f_red, data = data))
    structure_used <- if (isSingular(mod)) "site intercept only (still singular)"
    else                 "site intercept only"
  }
  pc <- do.call(rbind, lapply(split(data, data$curveID), function(d) {
    m <- lm(d[[col]] ~ d$Tc)
    data.frame(curveID = d$curveID[1], Species = d$Species[1], PC = d$PC[1],
               intercept = unname(coef(m)[1]), slope = unname(coef(m)[2]))
  }))
  pc$Species <- factor(pc$Species, levels = sp_levels)
  sc <- summary(mod)$coefficients
  list(percurve    = pc,
       trend_slope = fe_trend(mod, "Tc", "Tc:PC", pc_seq),
       trend_int   = fe_trend(mod, "(Intercept)", "PC", pc_seq),
       p_slope = sc["Tc:PC", "Pr(>|t|)"], df_slope = sc["Tc:PC", "df"],
       p_int   = sc["PC", "Pr(>|t|)"],    df_int   = sc["PC", "df"],
       structure = structure_used)
}

# Panels: per-curve points + population trend + ribbon (v = "E","gsw","WUE")
slope_panel <- function(res, v) {
  sym <- var_sym[[v]]; u <- var_usl[[v]]
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30", linewidth = 0.7) +
    geom_point(data = res$percurve, aes(x = PC, y = slope, color = Species),
               size = 2, alpha = 0.8) +
    geom_ribbon(data = res$trend_slope, aes(x = x, ymin = lwr, ymax = upr),
                fill = "grey70", alpha = 0.4) +
    geom_line(data = res$trend_slope, aes(x = x, y = fit), color = "black",
              linewidth = 0.7, linetype = p_linetype(res$p_slope)) +
    annotate("text", x = -Inf, y = Inf, label = p_fmt(res$p_slope), parse = TRUE,
             hjust = -0.20, vjust = 1.5, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_color_discrete(drop = FALSE, labels = sp_labs) +
    coord_cartesian(xlim = pc_lim) +
    labs(x = paste0(PC_AXIS, " (dimensionless)"),
         y = bquote(d * .(sym) / d * italic(T[leaf]) ~ "(" * .(u) * ")"),
         color = "Species",
         title = var_sym[[v]]) +
    theme_classic(base_size = 12)
}

int_panel <- function(res, v) {
  sym <- sym_level[[v]]; u <- var_ulv[[v]]
  ggplot() +
    geom_point(data = res$percurve, aes(x = PC, y = intercept, color = Species),
               size = 2, alpha = 0.8) +
    geom_ribbon(data = res$trend_int, aes(x = x, ymin = lwr, ymax = upr),
                fill = "grey70", alpha = 0.4) +
    geom_line(data = res$trend_int, aes(x = x, y = fit), color = "black",
              linewidth = 0.7, linetype = p_linetype(res$p_int)) +
    annotate("text", x = -Inf, y = Inf, label = p_fmt(res$p_int), parse = TRUE,
             hjust = -0.20, vjust = 1.5, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_color_discrete(drop = FALSE, labels = sp_labs) +
    coord_cartesian(xlim = pc_lim) +
    labs(x = paste0(PC_AXIS, " (dimensionless)"),
         y = bquote(.(sym) ~ "(" * .(u) * ")"),
         color = "Species") +
    theme_classic(base_size = 12)
}

# Do fitting
resE <- fit_var("E",   dat5)
resG <- fit_var("gsw", dat5)
resW <- fit_var("WUE", dat5)

# Report p and effective df
cat(sprintf("\n---- %s mixed models (Tleaf centered at 25 C, axis = %s) ----\n",
            FIG_TAG, PC_AXIS))
cat(sprintf("%-5s  %-24s %-24s  %s\n", "var",
            sprintf("slope (Tc:%s)", PC_AXIS),
            sprintf("level at 25 (%s)", PC_AXIS),
            "site structure"))
for (nm in c("E","gsw","WUE")) {
  r <- switch(nm, E = resE, gsw = resG, WUE = resW)
  cat(sprintf("%-5s  p=%.3g (df=%.1f)      p=%.3g (df=%.1f)   %s\n",
              nm, r$p_slope, r$df_slope, r$p_int, r$df_int, r$structure))
}

panels <- list(
  slope_panel(resE, "E"), slope_panel(resG, "gsw"), slope_panel(resW, "WUE"),
  int_panel(resE, "E"),   int_panel(resG, "gsw"),   int_panel(resW, "WUE")
)
legend <- get_legend(panels[[1]] + theme(legend.position = "right", legend.text = element_text()))
panels_noleg <- lapply(panels, function(p) p + theme(legend.position = "none"))
grid6 <- plot_grid(plotlist = panels_noleg, ncol = 3, nrow = 2, align = "hv",
                   labels = "AUTO", label_x = 0, label_y = 1.02)
Fig5 <- plot_grid(grid6, legend, ncol = 2, rel_widths = c(1, 0.28))

# Output PNG
ggsave(paste0(FIG_TAG, ".png"), plot = Fig5, width = 13, height = 8, units = "in",
       dpi = 300, bg = "white", limitsize = FALSE)

# Output PDF
ggsave(paste0(FIG_TAG, ".pdf"), plot = Fig5, width = 13, height = 8, device = cairo_pdf)