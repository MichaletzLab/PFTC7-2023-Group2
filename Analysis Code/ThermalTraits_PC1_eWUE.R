school.eWUE = read.csv("schoolfield.eWUE.discard.hooks.SANW.csv")

school.eWUE <- school.eWUE %>%
  rename(Species.y=taxon)%>%
  rename(Species.x=Species)%>%
  rename(Individual.y = turf_number)%>%
  rename(Individual.x = Individual)%>%
  rename(T_opt_school = T_opt)%>%
  rename(T_opt_SE_school = T_opt_SE)%>%
  rename(AIC_school = AIC)%>%
  rename(r_sq_school = r_sq)

# clean up data sets
fix_joined_cols <- function(df) {
  base_names <- unique(gsub("\\.(x|y)$", "", names(df)))
  out <- lapply(base_names, function(b) {
    xcol <- paste0(b, ".x")
    ycol <- paste0(b, ".y")
    if (xcol %in% names(df) & ycol %in% names(df)) {
      coalesce(df[[xcol]], df[[ycol]])
    } else if (xcol %in% names(df)) {
      df[[xcol]]
    } else if (ycol %in% names(df)) {
      df[[ycol]]
    } else {
      df[[b]]
    }})
  
  out <- as.data.frame(out)
  names(out) <- base_names
  out}

#Run function to clean data sets
schoolfield_params.e <- fix_joined_cols(school.eWUE)

#Join data sets and add a PC1 column corresponding to each site
parameter_dat.e <- schoolfield_params.e
PC1_df <- raw.env.data_pca %>%
  group_by(curveID) %>%
  summarise(PC1 = mean(PC1, na.rm = TRUE))

ThermTraits.dat.e <- left_join(parameter_dat.e, PC1_df, by="curveID")
#Run a gam with y=thermal traits, x=PC1
summary(mod.T_opt_sch <- gam(T_opt_school ~ s(PC1, k=3) + ###
                               Species,
                             data = ThermTraits.dat.e,
                             method="REML"))
summary(mod.Ea <- gam(E ~ s(PC1, k=3) +  ###
                        Species,
                      data = ThermTraits.dat.e,
                      method="REML"))
summary(mod.Ed <- gam(E_D ~ s(PC1, k=3) + ###
                        Species,
                      data = ThermTraits.dat.e,
                      method="REML"))
summary(mod.breadth <- gam(breadth_95 ~ s(PC1, k=3) +
                             Species,
                           data = ThermTraits.dat.e,
                           method="REML"))
##None are significant! Interestine... Small problem though because this loses 
#the error terms that I have stored (SE). I think this will overstate our 
#confidence though since the error is understated which means that since  
#our results are insignificant this is OK.

#Visualize:
e.School_Topt_Plot <- ggplot(ThermTraits.dat.e, aes(x = PC1, y = T_opt_school, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = expression(T[opt]~"(°C)"),
       color = "Species",
       title="iWUE (A/E)")
e.Ea_Plot <- ggplot(ThermTraits.dat.e, aes(x = PC1, y = E, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = expression(E[a]~"(eV)"),
       color = "Species")
e.Ed_Plot <- ggplot(ThermTraits.dat.e, aes(x = PC1, y = E_D, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = expression(E[d]~"(eV)"),
       color = "Species")
e.breadth_Plot <- ggplot(ThermTraits.dat.e, aes(x = PC1, y = breadth_95, color = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(x = "PC1",
       y = "breadth (°C)",
       color = "Species")

eWUE.Therm.plot <- ggarrange(e.School_Topt_Plot, e.breadth_Plot, e.Ea_Plot, e.Ed_Plot, nrow=2, ncol=2, common.legend = TRUE, labels = c("A","B","C","D"),legend="right")




# also run a linear like for E
library(lme4)
library(lmerTest)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(tibble)
library(ggpubr)

# --- Fit model
mod <- lmer(eWUE ~ Tleaf * PC1 + (1 + Tleaf | Species/curveID), data = raw.env.data_pca)
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
peW_slope <- ggplot(slope_data, aes(x = PC1, y = slope_Tleaf, color = Species)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "PC1 (Environmental Gradient)",
    y = expression(paste("Slope of ", "A/E", " vs. ", T[leaf])),
    color = "Species"
  ) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Interaction p = ", signif(pval_slope, 3)),
           hjust = 1.1, vjust = 1.5, size = 4.5) +
  theme_classic(base_size = 13)

# --- Intercept plot
peW_intercept <- ggplot(intercept_data, aes(x = PC1, y = intercept, color = Species)) +
  geom_point(alpha = 0.7) +
  labs(
    x = "PC1 (Environmental Gradient)",
    y = expression(paste("Intercept of ", "A/E", " vs. ", T[leaf])),
    color = "Species"
  ) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("p(PC1) = ", signif(pval_intercept, 3)),
           hjust = 1.1, vjust = 1.5, size = 4.5) +
  theme_classic(base_size = 13)

# --- Combine with labels A and B
ggarrange(peW_slope, peW_intercept,
          labels = c("A", "B"), common.legend=TRUE, legend="right",label.x = 0.1,
          ncol = 1, nrow = 2)
