topts.gam.add = read.csv("gam_fits_add.csv") #effects gam for tleaf and tleaf+cond 
topts.gam.raw = read.csv("gam_fits_raw.csv") #gam for tleaf and tleaf+cond
topts.glm.raw = read.csv("glm.fits.topts.csv") #glm for tleaf and tleaf+cond
schoolfield.fit = read.csv("outputs/discard.schoolfield.SANW.csv")
weibull.fit = read.csv("outputs/discard.weibull.SANW.csv")


t.test(topts.gam.add$tleaf_cond_breadth, topts.gam.add$weib_breadth,paired=TRUE)

t.test(topts.gam.add$tleaf_cond_opt, topts.gam.add$weib_Topt,paired=TRUE)


# Assuming you have the data stored in topts.gam.add
# Create a data frame for plotting
plot_data <- data.frame(
  tleaf_cond_opt = topts.gam.add$tleaf_cond_opt,
  weib_Topt = topts.gam.add$weib_Topt
)

ggplot(plot_data, aes(x = weib_Topt, y = tleaf_cond_opt)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(x = "Weibull T_opt Values", y = "GAM + Conditional T_opt Values") +
  theme_classic()


# Create a data frame for plotting
plot_data_breadth <- data.frame(
  tleaf_cond_breadth = topts.gam.add$tleaf_cond_breadth,
  weib_breadth = topts.gam.add$weib_breadth
)

# Plot using ggplot2
ggplot(plot_data_breadth, aes(x = weib_breadth, y = tleaf_cond_breadth)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(x = "Weibull breadth (˚C)", y = "GAM(tleaf + gsw) breadth (˚C)") +
  theme_classic() +
  xlim(0, 16) +
  ylim(0, 16)







t.test(topts.gam.add$tleaf_opt, topts.gam.add$weib_Topt,paired=TRUE)
t.test(topts.gam.add$tleaf_cond_opt, topts.gam.add$weib_Topt,paired=TRUE)
t.test(topts.gam.add$tleaf_vpdl_opt, topts.gam.add$weib_Topt,paired=TRUE) ## This one is significant
t.test(topts.gam.add$tleaf_cond_vpdl_opt, topts.gam.add$weib_Topt,paired=TRUE) ## Also significant


aic.tleaf = sum(topts.gam.add$tleaf_AIC)
aic.tleaf.cond= sum(topts.gam.add$tleaf_cond_AIC)
aic.tleaf.vpdl=sum(topts.gam.add$tleaf_vpdl_AIC)
aic.tleaf.cond.vpdl=sum(topts.gam.add$tleaf_cond_vpdl_AIC)
concurv.tleaf = mean(topts.gam.add$tleaf_concurvity)
concurv.tleaf.cond= mean(topts.gam.add$tleaf_cond_concurvity)
concurv.tleaf.vpdl=mean(topts.gam.add$tleaf_vpdl_concurvity)
concurv.tleaf.cond.vpdl=mean(topts.gam.add$tleaf_cond_vpdl_concurvity)

raic.tleaf = sum(topts.gam.raw$AIC_tleaf)
raic.tleaf.cond= sum(topts.gam.raw$AIC_tleaf_cond)
raic.tleaf.vpdl=sum(topts.gam.raw$AIC_tleaf_vpdl)
raic.tleaf.cond.vpdl=sum(topts.gam.raw$AIC_tleaf_vpdl_cond)
rconcurv.tleaf = mean(topts.gam.raw$concurvity_tleaf)
rconcurv.tleaf.cond= mean(topts.gam.raw$concurvity_tleaf_cond)
rconcurv.tleaf.vpdl=mean(topts.gam.raw$concurvity_tleaf_vpdl)
rconcurv.tleaf.cond.vpdl=mean(topts.gam.raw$concurvity_tleaf_vpdl_cond)

weibull.AIC = sum(weibull.fit$AIC)

summary_table <- data.frame(
  Model = c("tleaf", "tleaf + cond", "tleaf + vpdl", "tleaf + cond + vpdl"),
  Weibull.AIC = c(weibull.AIC,NA,NA,NA),
  AIC_add = c(aic.tleaf, aic.tleaf.cond, aic.tleaf.vpdl, aic.tleaf.cond.vpdl),
  Concurvity_add = c(concurv.tleaf, concurv.tleaf.cond, concurv.tleaf.vpdl, concurv.tleaf.cond.vpdl),
  AIC_raw = c(raic.tleaf, raic.tleaf.cond, raic.tleaf.vpdl, raic.tleaf.cond.vpdl),
  Concurvity_raw = c(rconcurv.tleaf, rconcurv.tleaf.cond, rconcurv.tleaf.vpdl, rconcurv.tleaf.cond.vpdl)
)

# Print the summary table
print(summary_table)





#Additive Gam mods 
t.test(topts.gam.add$tleaf_cond_opt, topts.gam.add$weib_Topt,paired=TRUE) #not significant
t.test(topts.gam.add$tleaf_opt, topts.gam.add$weib_Topt,paired=TRUE) #not significant
t.test(topts.gam.add$tleaf_cond_breadth, topts.gam.add$weib_breadth,paired=TRUE)
t.test(topts.gam.add$tleaf_breadth, topts.gam.add$weib_breadth,paired=TRUE)

#glm mods 
t.test(topts.glm.raw$topt_temp_cond, topts.glm.raw$topt_weib,paired=TRUE)#not significant
t.test(topts.glm.raw$topt_temp, topts.glm.raw$topt_weib,paired=TRUE)#not significant
t.test(topts.glm.raw$breadth_temp_cond, topts.glm.raw$breadth_weib,paired=TRUE)#not significant
t.test(topts.glm.raw$breadth_temp, topts.glm.raw$breadth_weib,paired=TRUE)

#Raw Gam 
t.test(topts.gam.raw$topt_temp_cond, topts.gam.raw$topt_weib,paired=TRUE)#not significant
t.test(topts.gam.raw$topt_temp, topts.gam.raw$topt_weib,paired=TRUE)#not significant
t.test(topts.gam.raw$breadth_temp_cond, topts.gam.raw$breadth_weib,paired=TRUE)
t.test(topts.gam.raw$breadth_temp, topts.gam.raw$breadth_weib,paired=TRUE)#not significant





# Function to create scatter plot with regression line and p-value
create_scatter_plot <- function(x, y, xlab, ylab, main, paired_t_test) {
  model <- lm(y ~ x)
  data <- data.frame(x = x, y = y)
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", col = "red") +
    geom_abline(intercept = 0, slope = 1, col = "blue", linetype = "dashed") +
    labs(x = xlab, y = ylab, title = main) +
    annotate("text", x = Inf, y = Inf, label = sprintf("p = %.3f", paired_t_test$p.value), 
             hjust = 1.1, vjust = 2, size = 4, color = "black")
  return(p)
}


############# Topt Plot #############
pairs <- list(
  list(x = topts.gam.add$tleaf_cond_opt, y = topts.gam.add$weib_Topt, xlab = "tleaf_cond_opt", ylab = "weib_Topt", main = "Additive Gam - tleaf_cond_opt vs weib_Topt"),
  list(x = topts.gam.add$tleaf_opt, y = topts.gam.add$weib_Topt, xlab = "tleaf_opt", ylab = "weib_Topt", main = "Additive Gam - tleaf_opt vs weib_Topt"),
  list(x = topts.glm.raw$topt_temp_cond, y = topts.glm.raw$topt_weib, xlab = "topt_temp_cond", ylab = "topt_weib", main = "GLM - topt_temp_cond vs topt_weib"),
  list(x = topts.glm.raw$topt_temp, y = topts.glm.raw$topt_weib, xlab = "topt_temp", ylab = "topt_weib", main = "GLM - topt_temp vs topt_weib"),
  list(x = topts.gam.raw$topt_temp_cond, y = topts.gam.raw$topt_weib, xlab = "topt_temp_cond", ylab = "topt_weib", main = "Raw Gam - topt_temp_cond vs topt_weib"),
  list(x = topts.gam.raw$topt_temp, y = topts.gam.raw$topt_weib, xlab = "topt_temp", ylab = "topt_weib", main = "Raw Gam - topt_temp vs topt_weib")
)
plots <- lapply(pairs, function(pair) {
  paired_t_test <- t.test(pair$x, pair$y, paired = TRUE)
  create_scatter_plot(pair$x, pair$y, pair$xlab, pair$ylab, pair$main, paired_t_test)
})
do.call(grid.arrange, c(plots, ncol = 2))




############# Breadth Plot ############# 
pairs <- list(
  list(x = topts.gam.add$tleaf_cond_breadth, y = topts.gam.add$weib_breadth, xlab = "tleaf_cond_breadth", ylab = "weib_breadth", main = "Additive Gam - tleaf_cond_breadth vs weib_breadth"),
  list(x = topts.gam.add$tleaf_breadth, y = topts.gam.add$weib_breadth, xlab = "tleaf_breadth", ylab = "weib_breadth", main = "Additive Gam - tleaf_breadth vs weib_breadth"),
  list(x = topts.glm.raw$breadth_temp_cond, y = topts.glm.raw$breadth_weib, xlab = "breadth_temp_cond", ylab = "breadth_weib", main = "GLM - breadth_temp_cond vs breadth_weib"),
  list(x = topts.glm.raw$breadth_temp, y = topts.glm.raw$breadth_weib, xlab = "breadth_temp", ylab = "breadth_weib", main = "GLM - breadth_temp vs breadth_weib"),
  list(x = topts.gam.raw$breadth_temp_cond, y = topts.gam.raw$breadth_weib, xlab = "breadth_temp_cond", ylab = "breadth_weib", main = "Raw Gam - breadth_temp_cond vs breadth_weib"),
  list(x = topts.gam.raw$breadth_temp, y = topts.gam.raw$breadth_weib, xlab = "breadth_temp", ylab = "breadth_weib", main = "Raw Gam - breadth_temp vs breadth_weib")
)
plots <- lapply(pairs, function(pair) {
  paired_t_test <- t.test(pair$x, pair$y, paired = TRUE)
  create_scatter_plot(pair$x, pair$y, pair$xlab, pair$ylab, pair$main, paired_t_test)
})
do.call(grid.arrange, c(plots, ncol = 2))



ggplot(aes(x=vpdl,y=photo))+geom_point
