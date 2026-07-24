# Figure 4 / Figure S4: Sharpe-Schoolfield parameters vs a principal component axis.
# Plot for A and iWUE, which model selection indicated are best characterized by
# a peaked/unimodal model. Then, regress fitted parameters on the PC axis.
#
# The axis and output names are parameterized. Run this file directly for
# Figure 4 (PC1). FigureS4.R sets PC_AXIS and FIG_TAG, then sources this file.

# ---- axis and output configuration -----------------------------------------
if (!exists("PC_AXIS")) PC_AXIS <- "PC1"
if (!exists("FIG_TAG")) FIG_TAG <- "Figure4"
message("Sharpe-Schoolfield parameter figure: axis = ", PC_AXIS,
        ", output tag = ", FIG_TAG)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(lmerTest)   # Satterthwaite p-values for lmer

source("Analysis Code/Fig_helpers.R")

fig4 <- read.csv("figure4_ss_parameters.csv", stringsAsFactors = FALSE)

# PC score (mean per curve) and Species from the point-level PCA object.
# The axis column is renamed to a neutral "PC" so that every formula,
# coefficient name, and aesthetic below is identical for PC1 and PC2.
meta <- raw.env.data_pca %>%
  group_by(curveID) %>%
  summarise(PC = mean(.data[[PC_AXIS]], na.rm = TRUE),
            Species = first(Species),
            Elevation = first(Elevation), .groups = "drop")

fig4 <- fig4 %>% left_join(meta, by = "curveID") %>% filter(interior)   # resolved optima only

fig4 <- fig4 %>% filter(!is.na(PC))  # remove any curves without a PC score (should only be curve 2)
fig4$Site <- factor(fig4$Elevation)  # site = elevation, for the site random intercept
fig4$Species <- factor(fig4$Species, levels = species_level_set())  # canonical level set

# ============================================================================
# Species-structure comparison for the PC slope.
# For each variable x parameter, fit three specifications and report the
# PC slope, its standard error, and its p-value:
#   random : lmerTest::lmer(y ~ PC + (1 | Species))  species random intercept
#   fixed  : lm(y ~ PC + Species)                    species fixed effects
#   naive  : lm(y ~ PC)                              no species; reference only
# random and fixed are the two candidate primary tests. naive is shown only
# to expose species-turnover artifacts and is not a candidate.
# ============================================================================

compare_params <- c("Topt", "Theta", "Ea", "Ed")
compare_rows   <- list()

for (v in c("A", "iWUE")) {
  for (p in compare_params) {
    d  <- fig4[fig4$variable == v, ]
    dd <- data.frame(y = d[[p]], PC = d$PC, Species = droplevels(d$Species))
    dd <- dd[stats::complete.cases(dd), ]
    n_curve <- nrow(dd)
    n_spp   <- nlevels(droplevels(dd$Species))
    
    m_ran <- suppressWarnings(lmerTest::lmer(y ~ PC + (1 | Species), data = dd))
    m_fix <- if (n_spp >= 2) lm(y ~ PC + Species, data = dd) else lm(y ~ PC, data = dd)
    m_nai <- lm(y ~ PC, data = dd)
    
    cr <- summary(m_ran)$coefficients
    cf <- summary(m_fix)$coefficients
    cn <- summary(m_nai)$coefficients
    
    compare_rows[[paste(v, p)]] <- data.frame(
      pc_axis      = PC_AXIS,
      variable     = v,
      parameter    = p,
      n_curve      = n_curve,
      n_species    = n_spp,
      slope_random = cr["PC", "Estimate"],
      se_random    = cr["PC", "Std. Error"],
      p_random     = cr["PC", "Pr(>|t|)"],
      slope_fixed  = cf["PC", "Estimate"],
      se_fixed     = cf["PC", "Std. Error"],
      p_fixed      = cf["PC", "Pr(>|t|)"],
      p_naive      = cn["PC", "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }
}

fig4_compare <- do.call(rbind, compare_rows)
rownames(fig4_compare) <- NULL

fig4_print <- fig4_compare
is_num <- sapply(fig4_print, is.numeric)
fig4_print[is_num] <- lapply(fig4_print[is_num], function(x) signif(x, 4))
cat(sprintf("\n==== %s: %s slope, random-intercept vs fixed-effects species ====\n",
            FIG_TAG, PC_AXIS))
print(fig4_print, row.names = FALSE)

compare_csv <- paste0(FIG_TAG, "_species_structure_comparison.csv")
write.csv(fig4_compare, compare_csv, row.names = FALSE)
cat("\nWrote:", compare_csv, "\n")

params <- c("Topt", "Theta", "Ea", "Ed")

# Fit the PRIMARY model once per variable x parameter:
# site + species random intercepts, falling back to species only where site is singular
models <- list()
for (v in c("A", "iWUE")) for (p in params) {
  d  <- fig4[fig4$variable == v, ]
  dd <- data.frame(y = d[[p]], PC = d$PC, Species = d$Species, Site = d$Site)
  m  <- suppressWarnings(lmerTest::lmer(y ~ PC + (1 | Site) + (1 | Species), data = dd))
  if (lme4::isSingular(m))
    m <- suppressWarnings(lmerTest::lmer(y ~ PC + (1 | Species), data = dd))
  models[[paste(v, p)]] <- m
}

for (nm in names(models))
  cat(nm, ":", if (grepl("Site", paste(deparse(formula(models[[nm]])), collapse = " "))) "site+species" else "species only", "\n")

pc_pval <- function(model) summary(model)$coefficients["PC", "Pr(>|t|)"]

# ---- Plot: line/ribbon come from the fixed effects of the primary model ----
yl <- list(Topt  = expression(italic(T[opt]) ~ "(" * degree * "C)"),
           Theta = expression(italic("\u03B8") ~ "(" * degree * "C)"),
           Ea    = expression(italic(E[a]) ~ "(eV)"),
           Ed    = expression(italic(E[d]) ~ "(eV)"))
titles <- list(A = expression(italic(A)), iWUE = "iWUE")

sp_levels <- levels(fig4$Species)
sp_labs   <- species_labels(sp_levels)

panel <- function(v, p) {
  d  <- fig4[fig4$variable == v, ]
  m  <- models[[paste(v, p)]]
  pv <- pc_pval(m)
  pr <- fe_trend(m, "(Intercept)", "PC",
                 seq(min(d$PC, na.rm = TRUE), max(d$PC, na.rm = TRUE), length.out = 100))
  ggplot(d, aes(x = PC, y = .data[[p]])) +
    geom_point(aes(color = Species), size = 2, alpha = 0.8) +
    geom_ribbon(data = pr, aes(x = x, ymin = lwr, ymax = upr),
                inherit.aes = FALSE, fill = "grey70", alpha = 0.4) +
    geom_line(data = pr, aes(x = x, y = fit), inherit.aes = FALSE,
              color = "black", linewidth = 0.6, linetype = p_linetype(pv)) +
    annotate("text", x = -Inf, y = Inf, label = p_fmt(pv), parse = TRUE,
             hjust = -0.20, vjust = 1.5, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_color_discrete(labels = sp_labs, drop = FALSE) +
    theme_classic(base_size = 12) + theme(legend.text = element_text()) +
    labs(x = paste0(PC_AXIS, " (dimensionless)"), y = yl[[p]], color = "Species",
         title = if (p == "Topt") titles[[v]] else NULL)
}

plots <- list()
for (p in params) { plots[[length(plots)+1]] <- panel("A", p)
plots[[length(plots)+1]] <- panel("iWUE", p) }
Fig4 <- ggarrange(plotlist = plots, ncol = 2, nrow = 4,
                  common.legend = TRUE, legend = "right", labels = "AUTO")

# output PNG
ggsave(
  filename = paste0(FIG_TAG, ".png"),
  plot = Fig4,
  width  = 8,
  height = 11,
  units  = "in",
  dpi    = 300,
  bg     = "white",
  limitsize = FALSE
)

# output PDF
ggsave(paste0(FIG_TAG, ".pdf"), Fig4, width = 8, height = 11, device = cairo_pdf)