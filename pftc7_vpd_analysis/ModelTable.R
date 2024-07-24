schoolfield.fit = read.csv("outputs/discard.schoolfield.SANW.csv")
weibull.fit = read.csv("outputs/discard.weibull.SANW.csv")
MJCschoolfield.fit = read.csv("MJCschoolfield.fits.csv")
gam.fits = read.csv("GAM.compare.csv")

conc.ave.0=mean(gam.fits$concurvity_full0.para)
conc.ave.1=mean(gam.fits$concurvity_full1.para)
conc.ave.2=mean(gam.fits$concurvity_full2.para)
conc.ave.3=mean(gam.fits$concurvity_full3.para)
conc.ave.4=mean(gam.fits$concurvity_full4.para)
conc.ave.5=mean(gam.fits$concurvity_full5.para)
conc.ave.6=mean(gam.fits$concurvity_full6.para)

sconc.ave.0=mean(gam.fits$concurvity_full0.s.tleaf.)
sconc.ave.1=mean(gam.fits$concurvity_full1.s.cond.)
sconc.ave.2=mean(gam.fits$concurvity_full2.s.cond.)
sconc.ave.20=mean(gam.fits$concurvity_full2.s.vpdl.)
sconc.ave.3=mean(gam.fits$concurvity_full3.s.tleaf.)
sconc.ave.4=mean(gam.fits$concurvity_full4.s.tleaf.)
sconc.ave.5=mean(gam.fits$concurvity_full5.s.tleaf.)
sconc.ave.50=mean(gam.fits$concurvity_full5.s.vpdl.)
sconc.ave.6=mean(gam.fits$concurvity_full6.s.tleaf.)
sconc.ave.60=mean(gam.fits$concurvity_full6.s.cond.)

weibull.rsq = mean(weibull.fit$r_sq)
weibull.AIC = sum(weibull.fit$AIC)
schoolfield.rsq = mean(schoolfield.fit$r_sq)
schoolfield.AIC = sum(schoolfield.fit$AIC)
MJCschoolfield.rsq = mean(MJCschoolfield.fit$r_sq)
MJCschoolfield.AIC = sum(MJCschoolfield.fit$AIC)
mod0.rsq = mean(gam.fits$rsq0)
mod0.AIC = sum(gam.fits$AIC0)
mod1.rsq = mean(gam.fits$rsq1)
mod1.AIC = sum(gam.fits$AIC1)
mod2.rsq = mean(gam.fits$rsq2)
mod2.AIC = sum(gam.fits$AIC2)
mod3.rsq = mean(gam.fits$rsq3)
mod3.AIC = sum(gam.fits$AIC3)
mod4.rsq = mean(gam.fits$rsq4)
mod4.AIC = sum(gam.fits$AIC4)
mod5.rsq = mean(gam.fits$rsq5)
mod5.AIC = sum(gam.fits$AIC5)
mod6.rsq = mean(gam.fits$rsq6)
mod6.AIC = sum(gam.fits$AIC6)

matrix.a <- lm(photo~ tleaf + vpdl, data = dats)
matrix.b <- lm(photo~ tleaf + cond, data = dats)
matrix.c <- lm(photo~ tleaf + tair, data = dats)
# Calculate VIF
(vif_valuesa <- vif(lm(matrix.a)))
(vif_valuesb <- vif(lm(matrix.b)))
(vif_valuesc <- vif(lm(matrix.c)))

summary_table <- data.frame(
  model = c("Weibull", "Schoolfield","Schoolfield x VPD","SS + s(tleaf)","SS + s(gsw)","SS + s(gsw) + s(vpd)","SSxVPD + s(tleaf)","s(tleaf)","s(tleaf) + s(vpd)","s(tleaf) + s(gsw)"),
  stat.method = c("NLS","NLS","NLS","GAM","GAM","GAM","GAM","GAM","GAM","GAM"),
  sum_AIC = c(weibull.AIC, schoolfield.AIC, MJCschoolfield.AIC,mod0.AIC,mod1.AIC,mod2.AIC,mod3.AIC,mod4.AIC,mod5.AIC,mod6.AIC),
  mean_rsq = c(weibull.rsq, schoolfield.rsq, MJCschoolfield.rsq,mod0.rsq,mod1.rsq,mod2.rsq,mod3.rsq,mod4.rsq,mod5.rsq,mod6.rsq),
  full_concurvity = c("NA","NA","NA",conc.ave.0,conc.ave.1,conc.ave.2,conc.ave.3,conc.ave.4,conc.ave.5,conc.ave.6),
  sm_concurvity_tleaf = c("NA","NA","NA",sconc.ave.0,"NA","NA",sconc.ave.3,sconc.ave.4,sconc.ave.5,sconc.ave.6),
  sm_concurvity_gsw = c("NA","NA","NA","NA",sconc.ave.1,sconc.ave.2,"NA","NA","NA",sconc.ave.60),
  sm_concurvity_vpd = c("NA","NA","NA","NA","NA",sconc.ave.20,"NA","NA",sconc.ave.50,"NA")
)

ESA.summary.table <- data.frame(
  model = c("Weibull", "Schoolfield", "s(tleaf)", "s(tleaf) + s(vpd)", "s(tleaf) + s(gsw)"),
  stat.method = c("NLS", "NLS", "GAM", "GAM", "GAM"),
  sum_AIC = c(weibull.AIC, schoolfield.AIC, mod4.AIC, mod5.AIC, mod6.AIC),
  mean_rsq = c(weibull.rsq, schoolfield.rsq, mod4.rsq, mod5.rsq, mod6.rsq),
  full_concurvity = c("NA", "NA", conc.ave.4, conc.ave.5, conc.ave.6),
  sm_concurvity_tleaf = c("NA", "NA", sconc.ave.4, sconc.ave.5, sconc.ave.6),
  sm_concurvity_gsw = c("NA", "NA", "NA", "NA", sconc.ave.60),
  sm_concurvity_vpd = c("NA", "NA", "NA", sconc.ave.50, "NA")
)

ESA.summary.table <- ESA.summary.table %>%
  select(-full_concurvity) %>%
  mutate(across(c(sum_AIC, mean_rsq, sm_concurvity_tleaf, sm_concurvity_gsw, sm_concurvity_vpd),
                ~ ifelse(is.na(as.numeric(.)), NA, sprintf("%.3f", as.numeric(.)))))

# Create the publication-ready table
kable(ESA.summary.table, 
      col.names = c("Model", "Method", "Sum AIC", "Mean R²", 
                    "Smooth Concurvity (Tleaf)", "Smooth Concurvity (Gsw)", "Smooth Concurvity (VPD)"),
      format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                full_width = F, 
                position = "left") %>%
  column_spec(3:7, width = "3cm") %>%
  add_header_above(c(" " = 2, "Statistics" = 5)) %>%
  kable_styling(latex_options = c("striped", "scale_down"))





# Extract necessary values from the model fit
new_model_AIC <- GAM.SS.tleaf.cond.fit$AIC
new_model_concurvity_full <- mean(GAM.SS.tleaf.cond.fit$concurvity["worst", ])
new_model_concurvity_sm_tleaf <- GAM.SS.tleaf.cond.fit$concurvity["worst", "s(tleaf)"]
new_model_concurvity_sm_cond <- GAM.SS.tleaf.cond.fit$concurvity["worst", "s(cond)"]

# Compute mean R-squared value (placeholder since it's not directly given by GAM model)
# Note: You might need to calculate R-squared based on your data and model predictions
new_model_rsq <- summary(GAM.SS.tleaf.cond.fit$gam_fit)$r.sq

# Create a new row with the model information
new_row <- data.frame(
  model = "gamSS + s(tleaf) + s(gsw)",
  stat.method = "GAM",
  sum_AIC = new_model_AIC,
  mean_rsq = new_model_rsq,
  full_concurvity = new_model_concurvity_full,
  sm_concurvity_tleaf = new_model_concurvity_sm_tleaf,
  sm_concurvity_gsw = new_model_concurvity_sm_cond,
  sm_concurvity_vpd = "NA"
)

# Append the new row to the summary table
summary_table <- rbind(summary_table, new_row)

# View the updated summary table
print(summary_table)

