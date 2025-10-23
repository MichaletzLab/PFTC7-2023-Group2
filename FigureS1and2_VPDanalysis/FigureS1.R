## Plot VPD vs. Tleaf==================
model.vpd.t <- gam(vpdl ~ s(tleaf, bs = 'ts', k = 5), data = dat)

# Extract the R-squared value
rsq <- summary(model.vpd.t)$r.sq
formatted_rsq1 <- sprintf("r² = %.3f", rsq)

# Create the plot
(p_out <- dat %>% 
    ggplot(aes(tleaf, vpdl)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = 'gam',
                formula = y ~ s(x, bs = 'ts', k = 5),
                color = "#cf0000") + 
    labs(x = "Leaf temperature (°C)",
         y = "VPD (kPa)") + 
    theme_classic() +
    annotate("text", x = 20, y = Inf, label = formatted_rsq1, 
             hjust = 1.1, vjust = 1.1, size = 5, fontface = "italic"))
