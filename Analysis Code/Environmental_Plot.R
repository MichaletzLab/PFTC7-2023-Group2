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
