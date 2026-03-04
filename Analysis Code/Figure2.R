(curve_plot <- raw.env.data %>% 
  filter(curveID == 1068)%>%
  ggplot(aes(Tleaf, A)) +
  geom_point() +
  labs(x = "Leaf temperature (°C)",
       y = expression(A~"(" * mu * mol ~ m^-2 ~ s^-1 * ")")) + 
  theme_classic())
#1 is great (n=299)

unique(raw.env.data$curveID)
count_rows <- raw.env.data%>%
  filter(curveID == 1068)%>%
  nrow()
count_rows
