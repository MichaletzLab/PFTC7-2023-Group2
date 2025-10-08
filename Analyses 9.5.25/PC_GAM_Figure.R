# --- Fit GAMs ---
gam_mod_all_Elev_A <- gam(A ~ s(Tleaf, k=3) +
                            s(Elevation, k=3) +
                            s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) +
                            ti(Tleaf, Elevation, k=3) +
                            ti(Tleaf, PC1, k=3) + ti(Tleaf, PC2, k=3) +
                            ti(Tleaf, PC3, k=3),
                          data = raw.env.data_pca,
                          method="REML")
gam_mod_all_Elev_E <- gam(E ~ s(Tleaf, k=3) +
                            s(Elevation, k=3) +
                            s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) +
                            ti(Tleaf, Elevation, k=3) +
                            ti(Tleaf, PC1, k=3) + ti(Tleaf, PC2, k=3) +
                            ti(Tleaf, PC3, k=3),
                          data = raw.env.data_pca,
                          method="REML")
gam_mod_all_Elev_gsw <- gam(gsw ~ s(Tleaf, k=3) +
                              s(Elevation, k=3) +
                              s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) +
                              ti(Tleaf, Elevation, k=3) +
                              ti(Tleaf, PC1, k=3) + ti(Tleaf, PC2, k=3) +
                              ti(Tleaf, PC3, k=3),
                            data = raw.env.data_pca,
                            method="REML")

# --- Helper for predictions ---
make_pred_df <- function(model, response, split_PC1=FALSE) {
  Tseq <- seq(min(raw.env.data_pca$Tleaf,na.rm=TRUE),
              max(raw.env.data_pca$Tleaf,na.rm=TRUE), length.out=200)
  base_df <- data.frame(
    Tleaf = Tseq,
    Elevation = median(raw.env.data_pca$Elevation, na.rm=TRUE),
    PC1 = mean(raw.env.data_pca$PC1, na.rm=TRUE),
    PC2 = mean(raw.env.data_pca$PC2, na.rm=TRUE),
    PC3 = mean(raw.env.data_pca$PC3, na.rm=TRUE)
  )
  if (!split_PC1) {
    pred <- predict(model, newdata=base_df, se.fit=TRUE)
    base_df$fit <- pred$fit; base_df$se <- pred$se.fit
    base_df$response <- response; base_df$level <- "All"
    return(base_df)
  } else {
    qs <- quantile(raw.env.data_pca$PC1, probs=c(0.1,0.5,0.9), na.rm=TRUE)
    names(qs) <- c("Low","Medium","High")
    out <- lapply(names(qs), function(lvl){
      newdf <- base_df
      newdf$PC1 <- qs[lvl]
      pred <- predict(model, newdata=newdf, se.fit=TRUE)
      newdf$fit <- pred$fit; newdf$se <- pred$se.fit
      newdf$response <- response; newdf$level <- lvl
      newdf
    })
    bind_rows(out)
  }
}

# --- Prediction data ---
pred_A_all <- make_pred_df(gam_mod_all_Elev_A,"A",FALSE)
pred_E_all <- make_pred_df(gam_mod_all_Elev_E,"E",FALSE)
pred_gsw_all <- make_pred_df(gam_mod_all_Elev_gsw,"gsw",FALSE)

pred_A_split <- make_pred_df(gam_mod_all_Elev_A,"A",TRUE)
pred_E_split <- make_pred_df(gam_mod_all_Elev_E,"E",TRUE)
pred_gsw_split <- make_pred_df(gam_mod_all_Elev_gsw,"gsw",TRUE)

# --- Function to format units for y-axis ---
response_units <- function(resp){
  switch(resp,
         "A"   = expression(A~(mu*mol~m^-2~s^-1)),
         "E"   = expression(E~(mmol~m^-2~s^-1)),
         "gsw" = expression(g[sw]~(mol~m^-2~s^-1)))
}

# --- Plotting ---
plot_unsplit <- function(pred, response){
  ggplot(pred, aes(x=Tleaf,y=fit))+
    geom_point(data=raw.env.data_pca,
               aes(x=Tleaf,y=.data[[response]]),
               inherit.aes=FALSE,alpha=0.02,size=0.5)+
    geom_line(color="darkgreen",linewidth=1)+
    geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se),
                fill="lightgreen",alpha=0.3)+
    theme_classic()+
    labs(x="Tleaf (°C)", y=response_units(response),
         title=NULL)
}

plot_split <- function(pred, response){
  # enforce factor order Low-Medium-High
  pred$level <- factor(pred$level, levels=c("Low","Medium","High"))
  
  # classify raw data
  qs <- quantile(raw.env.data_pca$PC1, probs=c(0.1,0.5,0.9), na.rm=TRUE)
  raw_split <- raw.env.data_pca %>%
    mutate(level = cut(PC1,
                       breaks=c(-Inf,qs[1],qs[2],qs[3],Inf),
                       labels=c("Low","Medium","High",NA),
                       include.lowest=TRUE)) %>%
    filter(!is.na(level))
  
  ggplot(pred, aes(x=Tleaf,y=fit,color=level,fill=level))+
    geom_point(data=raw_split,
               aes(x=Tleaf,y=.data[[response]],color=level),
               inherit.aes=FALSE,alpha=0.02,size=0.5)+
    geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=level),alpha=0.2,colour=NA,show.legend=FALSE)+
    geom_line(linewidth=1)+
    scale_color_manual(values=c(Low="lightblue",Medium="dodgerblue3",High="navy"),
                       breaks=c("Low","Medium","High"))+
    guides(fill="none")+  # hide separate fill legend
    theme_classic()+
    labs(x="Tleaf (°C)", y=response_units(response),
         color="PC1 level",
         title=NULL)
}

# --- Panels ---
pA_all  <- plot_unsplit(pred_A_all,"A")
pE_all  <- plot_unsplit(pred_E_all,"E")
pg_all  <- plot_unsplit(pred_gsw_all,"gsw")

pA_split <- plot_split(pred_A_split,"A")
pE_split <- plot_split(pred_E_split,"E")
pg_split <- plot_split(pred_gsw_split,"gsw")

PC.top <- ggarrange(pA_all, pE_all,pg_all, nrow = 1, labels=c("A","B","C"))
PC.bottom <- ggarrange(
  pA_split, pE_split, pg_split,
  nrow = 1,
  labels = c("D","E","F"),
  common.legend = TRUE,
  legend = "right",
  widths = c(1, 1, 1, 0.3)  # 3 plots get most of the width, legend gets less
)
ggarrange(PC.top, PC.bottom, nrow=2, ncol=1)


