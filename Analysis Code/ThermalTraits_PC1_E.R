library(lme4)
library(lmerTest)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(tibble)
library(ggpubr)

# --- Fit model
mod <- lmer(E ~ Tleaf * PC1 + (1 + Tleaf | Species/curveID), data = raw.env.data_pca)
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

# --- Slope plot
p_slope <- ggplot(slope_data, aes(x = PC1, y = slope_Tleaf, color = Species)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "PC1 (dimensionless)",
    y = expression(paste("Slope of ", E, " vs. ", T[leaf])),
    color = "Species"
  ) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("p = ", signif(pval_slope, 3)),
           hjust = 1.1, vjust = 8, size = 4.5) +
  theme_classic(base_size = 13)

# --- Intercept plot
p_intercept <- ggplot(intercept_data, aes(x = PC1, y = intercept, color = Species)) +
  geom_line() +
  labs(
    x = "PC1 (dimensionless)",
    y = expression(paste("Intercept of ", E, " vs. ", T[leaf])),
    color = "Species"
  ) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("p = ", signif(pval_intercept, 3)),
           hjust = 1.1, vjust = 1.5, size = 4.5) +
  theme_classic(base_size = 13)

# --- Combine with labels A and B
ggarrange(p_slope, p_intercept,
          labels = c("A", "B"), common.legend=TRUE, legend="right",label.x = 0.1,
          ncol = 1, nrow = 2)
