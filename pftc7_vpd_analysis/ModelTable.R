schoolfield.fit = read.csv("outputs/discard.schoolfield.SANW.csv")
weibull.fit = read.csv("outputs/discard.weibull.SANW.csv")
MJCschoolfield.fit = read.csv("MJCschoolfield.fits.csv")

weibull.rsq = mean(weibull.fit$r_sq)
weibull.AIC = sum(weibull.fit$AIC)
schoolfield.rsq = mean(schoolfield.fit$r_sq)
schoolfield.AIC = sum(schoolfield.fit$AIC)
MJCschoolfield.rsq = mean(MJCschoolfield.fit$r_sq)
MJCschoolfield.AIC = sum(MJCschoolfield.fit$AIC)

summary_table <- data.frame(
  model = c("weibull", "schoolfield","schoolfield x vpd"),
  stat.method = c("NLS","NLS","NLS"),
  sum_AIC = c(weibull.AIC, schoolfield.AIC, MJCschoolfield.AIC),
  mean_rsq = c(weibull.rsq, schoolfield.rsq, MJCschoolfield.rsq)
)
summary_table
