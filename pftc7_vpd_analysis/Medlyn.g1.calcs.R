med.dat = dat%>%
  group_by(curveid)%>%
  mutate(medlyn=(1.6*photo/(ca/sqrt(vpdl))))
  

F=as.formula(y~x)
medlyn.fig= med.dat%>%
    ggplot(aes(x =medlyn, y=cond, color=factor(site))) + geom_point(shape=18)+
    theme_classic()+theme(axis.text.x=element_text(size=15))+
    xlab(expression(1.6*A/(Ca*sqrt(VPD))))+ylab(expression('g'[s]* ' (mol / (m'^2*'s)'*''))+
    stat_smooth(aes(fill=factor(site),color=factor(site)),method="lm",formula=F)+
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")
medlyn.fig
# Calculate the slopes, standard errors
#med.dat$elevation <- factor(med.dat$elevation, levels = c("469", "700", "920", "1290", "2000","2200","2400","2600","2800"))
slopes <- med.dat %>%
  group_by(elevation) %>%
  do({
    model <- lm(cond ~ medlyn, data = .)
    slope <- coef(model)[[2]]
    se_slope <- summary(model)$coefficients[2, "Std. Error"]
    data.frame(slope = slope, se_slope = se_slope)
  })
specific_breaks <- c(469, 700, 920, 1290, 2000, 2200, 2400, 2600, 2800)
# Plot the g1 values
g1.plot=(ggplot(slopes, aes(x = elevation, y = slope))) +
            geom_point(position = position_dodge(width = 0.75), size = 3) +
            geom_errorbar(aes(ymin = slope - se_slope, ymax = slope + se_slope),width = 0,
                          position = position_dodge(width = 0.75)) +
            scale_x_continuous(breaks = specific_breaks) +
            theme_classic() +theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 12),
                                   axis.title.y = element_text(size = 15),legend.title= element_blank()) +
            ylab(expression('g'[1])) + xlab("Elevation(masl)")
g1.plot

