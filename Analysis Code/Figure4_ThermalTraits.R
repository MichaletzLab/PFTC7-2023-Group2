#First run all four ThermalTraits_PC1 scripts
library(cowplot)

legend_plot <- School_Topt_Plot +
  theme(
    legend.position = "right",
    legend.text = element_markdown(size = 15),
    legend.title = element_text(size = 15)
  )

shared_legend <- get_legend(legend_plot)

plots_no_legends <- list(
  School_Topt_Plot, g.School_Topt_Plot, e.School_Topt_Plot, i.School_Topt_Plot,
  breadth_Plot, g.breadth_Plot, e.breadth_Plot, i.breadth_Plot,
  Ea_Plot, g.Ea_Plot, e.Ea_Plot, i.Ea_Plot,
  Ed_Plot, g.Ed_Plot, e.Ed_Plot, i.Ed_Plot
) |> 
  lapply(\(p) p + theme(legend.position = "none"))

plot_grid_combined <- plot_grid(
  plotlist = plots_no_legends,
  ncol = 4, nrow = 4,
  align = "hv"
)

final_plot <- plot_grid(
  plot_grid_combined, shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.22)
)

final_plot
#export as 1400 x 1000
