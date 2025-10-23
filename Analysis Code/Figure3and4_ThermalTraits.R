#First run all four ThermalTraits_PC1 scripts
ggarrange(School_Topt_Plot, g.School_Topt_Plot, e.School_Topt_Plot,
          i.School_Topt_Plot, breadth_Plot, g.breadth_Plot, 
          e.breadth_Plot, i.breadth_Plot, Ea_Plot, g.Ea_Plot, 
          e.Ea_Plot, i.Ea_Plot, Ed_Plot, g.Ed_Plot, e.Ed_Plot,
          i.Ed_Plot, nrow=4, ncol=4, common.legend = TRUE, 
          labels = c("A","B","C","D","E","F","G","H","I","J","K",
                     "L","M","N","O","P"),legend="right")
#export as 1180 x 900
