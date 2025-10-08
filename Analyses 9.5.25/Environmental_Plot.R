# --- Temperature palette (kept as-is) ---
temp_colors <- c("Soil temp" = "#4A2125",   # brownish red
                 "Ground temp" = "#4DAF4A", # medium green
                 "Air temp"   = "#377EB8")  # medium blue

# --- Complementary palette for other plots ---
# Harmonized with the above but distinct
country_colors <- c("S. Africa" = "#007A4D",   # warm golden orange
                    "Norway"    = "#8B0000")   # cool sky blue

#This one uses all the data
temp.plot <- env_all_long %>%
  filter(Variable %in% c("Soil temp", "Ground temp", "Air temp")) %>%
  ggplot(aes(x = Elevation, y = Value,
             color = Variable,
             linetype = Country)) +
  geom_point(aes(shape = Country), alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.25) +
  scale_color_manual(values = temp_colors) +
  theme_classic(base_size = 14) +
  labs(x = "Elevation (m a.s.l.)",
       y = "Temperature (°C)",
       color = "Temperature class",
       shape = "Country",
       linetype = "Country") +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

# --- Moisture vs. Elevation ---
#But this one uses all the data
moist.plot <- env_all_long %>%
  filter(Variable == "Soil moisture") %>%
  ggplot(aes(x = Elevation, y = Value, color = Country)) +
  geom_point(aes(shape = Country), alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.25) +
  scale_color_manual(values = country_colors) +
  theme_classic(base_size = 14) +
  labs(x = "Elevation (m a.s.l.)",
       y = "Soil moisture (%)",
       color = "Country",
       shape = "Country") +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )


# --- Vegetation height vs. Elevation ---
height.plot <- ggplot(raw.env.data_pca, aes(x = Elevation, y = vegetation_height, color = Country)) +
  geom_point(aes(shape = Country), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.25) +
  scale_color_manual(values = country_colors) +
  theme_classic(base_size = 14) +
  labs(x = "Elevation (masl)", y = "Vegetation height (cm)", color = "Country", shape = "Country")

# --- PC1 vs. Elevation ---
pc.plot <- ggplot(raw.env.data_pca, aes(x = Elevation, y = PC1, color = Country)) +
  geom_point(aes(shape = Country), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.25) +
  scale_color_manual(values = country_colors) +
  theme_classic(base_size = 14) +
  labs(x = "Elevation (masl)", y = "PC1", color = "Country", shape = "Country")

# --- Extract legends ---
legend_country <- get_legend(
  moist.plot + theme(legend.position = "bottom")
)

legend_temp <- get_legend(
  temp.plot + guides(shape = "none") + theme(legend.position = "bottom")
)

# --- Remove legends from main plots ---
moist.plot_noleg  <- moist.plot  + theme(legend.position = "none")
temp.plot_noleg   <- temp.plot   + theme(legend.position = "none")
height.plot_noleg <- height.plot + theme(legend.position = "none")
pc.plot_noleg     <- pc.plot     + theme(legend.position = "none")

# --- Combine 4-panel ---
plot.4way <- ggarrange(
  moist.plot_noleg, temp.plot_noleg, height.plot_noleg, pc.plot_noleg,
  nrow = 2, ncol = 2,
  labels = c("A", "B", "C", "D"),
  label.x = 0.05, label.y = 0.98
)

# --- Combine with PCA and legends ---
final_plot <- plot_grid(
  plot.4way, pca_plot,
  nrow = 1, rel_widths = c(1.3, 1),
  labels = c("", "E"), label_x = 0.05
)

final_with_legends <- plot_grid(
  final_plot,
  legend_country,
  legend_temp,
  ncol = 1,
  rel_heights = c(1, 0.08, 0.1)
)

final_with_legends