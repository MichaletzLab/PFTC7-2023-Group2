library(lme4)
library(lmerTest)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(tibble)
library(ggtext)
library(ggpubr)

# --- Fit model
mod <- lmer(E ~ Tleaf * PC1 + (1 + Tleaf || Species/curveID),
            data = raw.env.data_pca)
summary(mod)
# --- Extract fixed effects
coefs <- fixef(mod)
b1 <- coefs["Tleaf"]
b3 <- coefs["Tleaf:PC1"]
b0 <- coefs["(Intercept)"]
b_PC1 <- coefs["PC1"]
# --- Get species-level random effects
species_ranefs <- ranef(mod)$Species %>%
  as.data.frame() %>%
  rename(intercept_r = `(Intercept)`, slope_r = Tleaf) %>%
  rownames_to_column("Species")
# --- Get full range of PC1
pc1_vals <- unique(raw.env.data_pca$PC1)
# --- Compute per-species slope lines
slope_data <- expand.grid(
  PC1 = pc1_vals,
  Species = unique(species_ranefs$Species)
) %>%
  left_join(species_ranefs, by = "Species") %>%
  mutate(
    slope_Tleaf = (b1 + slope_r) + b3 * PC1
  )
# --- Compute per-species intercept lines
intercept_data <- expand.grid(
  PC1 = pc1_vals,
  Species = unique(species_ranefs$Species)
) %>%
  left_join(species_ranefs, by = "Species") %>%
  mutate(
    intercept = (b0 + intercept_r) + b_PC1 * PC1  # intercept changes with PC1
  )
# --- Extract p-values
pval_slope <- summary(mod)$coefficients["Tleaf:PC1", "Pr(>|t|)"]
pval_intercept <- summary(mod)$coefficients["PC1", "Pr(>|t|)"]

# Build plotmath labels: italic binomial, roman "var."
species_labels <- function(levels) {
  out <- lapply(levels, function(s) {
    if (grepl(" var. ", s, fixed = TRUE)) {
      parts <- strsplit(s, " var. ", fixed = TRUE)[[1]]
      g <- strsplit(parts[1], " ")[[1]]     # genus + species
      bquote(italic(.(g[1]))~italic(.(g[2]))~"var."~italic(.(parts[2])))
    } else {
      g <- strsplit(s, " ")[[1]]
      if (length(g) == 2) bquote(italic(.(g[1]))~italic(.(g[2])))
      else                bquote(italic(.(s)))   # safety for odd names
    }
  })
  do.call(expression, out)
}

sp_levels <- levels(factor(slope_data$Species))
sp_labs   <- species_labels(sp_levels)


# --- Slope plot
p_slope <- ggplot(slope_data, aes(x = PC1, y = slope_Tleaf, color = Species)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "PC1 (dimensionless)",
    y = expression("Slope of " * italic(E) ~ "vs." ~ italic(T)[leaf]),
    color = "Species"
  ) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("p = ", signif(pval_slope, 3)),
           hjust = 1.1, vjust = 8, size = 4.5) +
  scale_color_discrete(labels = sp_labs) +
  theme_classic(base_size = 13) +
  theme(legend.text = element_text())

# --- Intercept plot
p_intercept <- ggplot(intercept_data, aes(x = PC1, y = intercept, color = Species)) +
  geom_line() +
  labs(
    x = "PC1 (dimensionless)",
    y = expression("Intercept of " * italic(E) ~ "vs." ~ italic(T)[leaf]),
    color = "Species"
  ) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("p = ", signif(pval_intercept, 3)),
           hjust = 1.1, vjust = 1.5, size = 4.5) +
  scale_color_discrete(labels = sp_labs) +
  theme_classic(base_size = 13) +
  theme(legend.text = element_text())

# --- Combine with labels A and B
library(cowplot)
# Extract legend from one plot
legend <- get_legend(
  p_slope + theme(legend.position = "right")
)
# Remove legends from individual plots
p_slope_noleg <- p_slope + theme(legend.position = "none")
p_intercept_noleg <- p_intercept + theme(legend.position = "none")
# Combine plots and legend manually
combined <- plot_grid(
  p_slope_noleg, p_intercept_noleg,
  labels = c("A", "B"),
  ncol = 1, align = "v",
  label_x = 0,   # move right, away from y-axis
  label_y = 1.02   # move slightly above the panel
)
# Add shared legend
Fig5 <- plot_grid(combined, legend, rel_widths = c(1, 0.5))
#700 x 500

# output PNG
ggsave("Figure5.png", plot = Fig5,
       width = 9, height = 5, units = "in",
       dpi = 300, bg = "white", limitsize = FALSE)

# output PDF
ggsave("Figure5.pdf", plot = Fig5,
       width = 9, height = 5, device = cairo_pdf)
