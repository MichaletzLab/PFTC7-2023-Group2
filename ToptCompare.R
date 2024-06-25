topts.gam.add = read.csv("gam.fits.add.topts.csv")
topts.gam.raw = read.csv("gam.fits.raw.topts.csv")
topts.glm.raw = read.csv("glm.fits.topts.csv")

#Additive Gam mods 
t.test(topts.gam.add$tleaf_cond_opt, topts.gam.add$weib_Topt,paired=TRUE)
t.test(topts.gam.add$tleaf_opt, topts.gam.add$weib_Topt,paired=TRUE)

#glm mods 
t.test(topts.glm.raw$topt_temp_cond, topts.glm.raw$topt_weib,paired=TRUE)
t.test(topts.glm.raw$topt_temp, topts.glm.raw$topt_weib,paired=TRUE)

#Raw Gam 
t.test(topts.gam.raw$topt_temp_cond, topts.gam.raw$topt_weib,paired=TRUE)
t.test(topts.gam.raw$topt_temp, topts.gam.raw$topt_weib,paired=TRUE)








# Function to create scatter plot with regression line and p-value
create_scatter_plot <- function(x, y, xlab, ylab, main, paired_t_test) {
  # Fit linear model
  model <- lm(y ~ x)
  
  # Create a data frame for ggplot
  data <- data.frame(x = x, y = y)
  
  # Create the scatter plot
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", col = "red") +
    geom_abline(intercept = 0, slope = 1, col = "blue", linetype = "dashed") +
    labs(x = xlab, y = ylab, title = main) +
    annotate("text", x = Inf, y = Inf, label = sprintf("p = %.3f", paired_t_test$p.value), 
             hjust = 1.1, vjust = 2, size = 4, color = "black")
  
  return(p)
}

# Define the pairs of variables and their labels
pairs <- list(
  list(x = topts.gam.add$tleaf_cond_opt, y = topts.gam.add$weib_Topt, xlab = "tleaf_cond_opt", ylab = "weib_Topt", main = "Additive Gam - tleaf_cond_opt vs weib_Topt"),
  list(x = topts.gam.add$tleaf_opt, y = topts.gam.add$weib_Topt, xlab = "tleaf_opt", ylab = "weib_Topt", main = "Additive Gam - tleaf_opt vs weib_Topt"),
  list(x = topts.glm.raw$topt_temp_cond, y = topts.glm.raw$topt_weib, xlab = "topt_temp_cond", ylab = "topt_weib", main = "GLM - topt_temp_cond vs topt_weib"),
  list(x = topts.glm.raw$topt_temp, y = topts.glm.raw$topt_weib, xlab = "topt_temp", ylab = "topt_weib", main = "GLM - topt_temp vs topt_weib"),
  list(x = topts.gam.raw$topt_temp_cond, y = topts.gam.raw$topt_weib, xlab = "topt_temp_cond", ylab = "topt_weib", main = "Raw Gam - topt_temp_cond vs topt_weib"),
  list(x = topts.gam.raw$topt_temp, y = topts.gam.raw$topt_weib, xlab = "topt_temp", ylab = "topt_weib", main = "Raw Gam - topt_temp vs topt_weib")
)

# Perform paired t-tests and create scatter plots
plots <- lapply(pairs, function(pair) {
  paired_t_test <- t.test(pair$x, pair$y, paired = TRUE)
  create_scatter_plot(pair$x, pair$y, pair$xlab, pair$ylab, pair$main, paired_t_test)
})

# Arrange and display the plots
library(gridExtra)
do.call(grid.arrange, c(plots, ncol = 2))
