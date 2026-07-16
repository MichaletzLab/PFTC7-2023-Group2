(curve_plot <- raw.env.data %>% 
  filter(curveID == 1068)%>%
  ggplot(aes(Tleaf, A)) +
  geom_point(alpha = 0.5, size = 1.8, shape = 16, color = "chartreuse4") +
  labs(x = "Leaf temperature (°C)",
       y = expression(italic(A)~"(" * mu * mol ~ m^-2 ~ s^-1 * ")")) + 
  theme_classic(base_size = 16))
#1 is great (n=299)

unique(raw.env.data$curveID)
count_rows <- raw.env.data%>%
  filter(curveID == 1068)%>%
  nrow()
count_rows

# get species of plotted curve
raw.env.data %>%
  filter(curveID == 1068) %>%
  distinct(Species)

# Save file
ggsave(
  filename = "Figure2.png",
  plot     = curve_plot,
  width    = 6,
  height   = 5,
  units    = "in",
  dpi      = 300,
  limitsize = FALSE
)
