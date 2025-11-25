(curve_plot <- raw.env.data %>% 
  filter(curveID == 1)%>%
  ggplot(aes(Tleaf, A)) +
  geom_point() +
  labs(x = "Leaf temperature (°C)",
       y = expression(A~"(" * mu * mol ~ m^-2 ~ s^-1 * ")")) + 
  theme_classic())
#1,5,