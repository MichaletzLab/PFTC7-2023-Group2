# --- Temperature palette (kept as-is) ---
temp_colors <- c("Soil temp"   = "#b05a3c",  # warm earthy orange
                 "Ground temp" = "#9ab973",  # muted olive green
                 "Air temp"    = "#648fff")  # cool blue

country_colors <- c("S. Africa" = "#377F6C",  # muted teal
                    "Norway"    = "#984464")  # plum red

# --- Temperature plot ---
temp.plot <- env_all_long %>%
  filter(Variable %in% c("Soil temp", "Ground temp", "Air temp")) %>%
  ggplot(aes(x = Elevation, y = Value,
             color = Variable,
             shape = Country)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(aes(linetype = Country),  # <-- added this
              method = "lm", se = TRUE, linewidth = 1, alpha = 0.15) +
  scale_color_manual(values = temp_colors, name = "Temperature class") +
  scale_linetype_manual(values = c("Norway" = "solid", "S. Africa" = "dashed"), 
                        name = "Country") +
  scale_shape_manual(values = c("Norway" = 16, "S. Africa" = 17), 
                     name = "Country") +
  theme_classic(base_size = 14) +
  labs(x = "Elevation (m a.s.l.)", y = "Temperature (°C)") +
  theme(legend.position = "right",
        legend.box = "vertical")

# --- Moisture vs. Elevation ---
moist.plot <- env_all_long %>%
  filter(Variable == "Soil moisture") %>%
  ggplot(aes(x = Elevation, y = Value,
             color = Country, shape = Country)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(aes(linetype = Country),  # <-- added this
              method = "lm", se = TRUE, linewidth = 1, alpha = 0.15) +
  scale_color_manual(values = country_colors, name = "Country") +
  scale_linetype_manual(values = c("Norway" = "solid", "S. Africa" = "dashed"), 
                        name = "Country") +
  scale_shape_manual(values = c("Norway" = 16, "S. Africa" = 17), 
                     name = "Country") +
  theme_classic(base_size = 14) +
  labs(x = "Elevation (m a.s.l.)", y = "Soil moisture (%)") +
  theme(legend.position = "right")

# --- Vegetation height vs. Elevation ---
height.plot <- ggplot(raw.env.data_pca,
                      aes(x = Elevation, y = vegetation_height,
                          color = Country, shape = Country)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(aes(linetype = Country),  # <-- added this
              method = "lm", se = TRUE, linewidth = 1, alpha = 0.15) +
  scale_color_manual(values = country_colors, name = "Country") +
  scale_linetype_manual(values = c("Norway" = "solid", "S. Africa" = "dashed"),
                        name = "Country") +
  scale_shape_manual(values = c("Norway" = 16, "S. Africa" = 17),
                     name = "Country") +
  theme_classic(base_size = 14) +
  labs(x = "Elevation (m a.s.l.)", y = "Vegetation height (cm)") +
  theme(legend.position = "none")

# --- PC1 vs. Elevation ---
pc.plot <- ggplot(raw.env.data_pca,
                  aes(x = Elevation, y = PC1,
                      color = Country, shape = Country)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(aes(linetype = Country),  # <-- added this
              method = "lm", se = TRUE, linewidth = 1, alpha = 0.15) +
  scale_color_manual(values = country_colors, name = "Country") +
  scale_linetype_manual(values = c("Norway" = "solid", "S. Africa" = "dashed"),
                        name = "Country") +
  scale_shape_manual(values = c("Norway" = 16, "S. Africa" = 17),
                     name = "Country") +
  theme_classic(base_size = 14) +
  labs(x = "Elevation (m a.s.l.)", y = "PC1 (dimensionless)") +
  theme(legend.position = "none")


# --- Extract clean legends ---
legend_country <- get_legend(
  moist.plot + theme(legend.position = "bottom",
                     legend.box = "horizontal",
                     legend.title = element_text(size = 12))
)

legend_temp <- get_legend(
  temp.plot +
    guides(linetype = "none", shape = "none") +  # show only color legend
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.title = element_text(size = 12))
)

# --- Remove legends from main plots ---
temp.plot_noleg   <- temp.plot   + theme(legend.position = "none")
moist.plot_noleg  <- moist.plot  + theme(legend.position = "none")
height.plot_noleg <- height.plot + theme(legend.position = "none")
pc.plot_noleg     <- pc.plot     + theme(legend.position = "none")

# --- Combine 4-panel ---
plot.4way <- ggarrange(
  moist.plot_noleg, temp.plot_noleg,
  height.plot_noleg, pc.plot_noleg,
  nrow = 2, ncol = 2,
  labels = c("A", "B", "C", "D"),
  label.x = 0.11,   # move letters slightly right
  label.y = 0.98,
  font.label = list(size = 14, face = "bold")
)

# --- Combine with PCA and legends ---
final_plot <- plot_grid(
  plot.4way, pca_plot,
  nrow = 1, rel_widths = c(1.3, 1),
  labels = c("", "E"), label_x = 0.05
)

# --- Add separate legends (non-overlapping) ---
final_with_legends <- plot_grid(
  final_plot,
  legend_temp,
  legend_country,
  ncol = 1,
  rel_heights = c(1, 0.08, 0.08)
)

final_with_legends # 990 x 550

###########################################
## P-values for text: ####
###########################################
summary(lm(PC1 ~ Elevation * Country, data = raw.env.data_pca))
summary(lm(vegetation_height ~ Elevation * Country, data = raw.env.data_pca))
summary(lm(Value ~ Elevation * Country,
           data = env_all_long %>% filter(Variable == "Soil moisture")))

summary(lm(Value ~ Elevation * Country * Variable,
           data = env_all_long %>% 
             filter(Variable %in% c("Soil temp", "Ground temp", "Air temp"))))















library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

# --- Make figures with PC1 as x ----
# --- Temperature plot (keep the legend here) ---
temp.plot.pc1 <- env_all_long %>%
  filter(Variable %in% c("Soil temp", "Ground temp", "Air temp")) %>%
  ggplot(aes(x = PC1, y = Value, color = Variable)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(aes(group = Variable, color = Variable),
              method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = temp_colors, name = "Temperature class") +
  theme_classic(base_size = 14) +
  labs(x = "PC1 (dimensionless)", y = "Temperature (°C)") +
  theme(legend.position = "bottom",  # <-- Keep legend here
        legend.box = "horizontal")

# --- Moisture plot ---
moist.plot.pc1 <- env_all_long %>%
  filter(Variable == "Soil moisture") %>%
  ggplot(aes(x = PC1, y = Value)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "black") +
  theme_classic(base_size = 14) +
  labs(x = "PC1 (dimensionless)", y = "Soil moisture (%)") +
  theme(legend.position = "none")

# --- Vegetation height plot ---
height.plot.pc1 <- ggplot(raw.env.data_pca,
                          aes(x = PC1, y = vegetation_height)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.15, color = "black") +
  theme_classic(base_size = 14) +
  labs(x = "PC1 (dimensionless)", y = "Vegetation height (cm)") +
  theme(legend.position = "none")

# --- PC1 vs Elevation ---
pc.plot.pc1 <- ggplot(raw.env.data_pca,
                      aes(y = Elevation, x = PC1)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.15, color = "black") +
  theme_classic(base_size = 14) +
  labs(x = "PC1 (dimensionless)", y = "Elevation (m a.s.l.)") +
  theme(legend.position = "none")

# --- Extract the temperature legend BEFORE removing it from the plot ---
legend_temp <- get_legend(temp.plot.pc1)

# --- Remove legend from main plot now ---
temp.plot_noleg <- temp.plot.pc1 + theme(legend.position = "none")

# --- Combine the four PC1 plots into 2x2 grid ---
plot.4way.pc1 <- ggarrange(
  moist.plot.pc1, temp.plot_noleg,
  height.plot.pc1, pc.plot.pc1,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  label.x = 0.11,   # move letters slightly right
  label.y = 0.98,
  font.label = list(size = 14, face = "bold")
)

# --- Combine with PCA plot (E) ---
main_plot <- plot_grid(
  plot.4way.pc1, pca_plot,
  nrow = 1,
  rel_widths = c(1.2, 1),
  labels = c("", "E"),
  label_size = 14,
  align = "h"
)

# --- Add the legend below everything ---
final_plot <- plot_grid(
  main_plot,
  legend_temp,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# --- Display final figure ---
final_plot

#P-vals
# --- P-values for each plotted linear model ---

# 1. Temperature models (one for each temperature variable)
temp_models <- env_all_long %>%
  filter(Variable %in% c("Soil temp", "Ground temp", "Air temp")) %>%
  group_by(Variable) %>%
  do(model = lm(Value ~ PC1, data = .)) %>%
  mutate(summary = purrr::map(model, summary),
         slope = purrr::map_dbl(summary, ~ coef(.x)[2, 1]),
         p_value = purrr::map_dbl(summary, ~ coef(.x)[2, 4])) %>%
  select(Variable, slope, p_value)

temp_models
# 2. Moisture model
moist_model <- lm(Value ~ PC1,
                  data = env_all_long %>% filter(Variable == "Soil moisture"))
summary(moist_model)
# 3. Vegetation height model
height_model <- lm(vegetation_height ~ PC1, data = raw.env.data_pca)
summary(height_model)
# 4. Elevation vs PC1 model
elev_model <- lm(PC1 ~ Elevation, data = raw.env.data_pca)
summary(elev_model)

