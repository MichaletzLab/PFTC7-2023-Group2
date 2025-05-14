# Code to read LI-6800 files -- set working directory to the folder of licor files (not excel)

library(tidyverse)
library(readr)
library(nls.multstart)

read_6800 <- function(filename) {
  
  raw_input = readLines(filename)                # Read in raw datafile
  data_start = grep("\\[Data\\]", raw_input) + 2 # Find where data starts (at "[Data]" line)
  data_end = grep("\\[Header\\]", raw_input) - 1 # Find where data ends (at "[Header]" line)
  data_end = c(data_end, length(raw_input))[-1]  # Add data end at end of file
  
  data_end = data_end[!(data_end < data_start[1])]
  
  compiled_data = c() # Initiate blank holder
  
  # Loop over each data set in the file
  for (i in 1:length(data_start)) {
    trimmed = raw_input[data_start[i]:data_end[i]]      # Grab data
    trimmed = trimmed[-2]                               # Remove units line
    current_data = read.csv(text = trimmed, sep = "\t") # Convert to dataframe
    compiled_data = rbind(compiled_data, current_data)  # Merge together into one 
  }
  return(compiled_data)
}
folder <- list.files()

lrc_mean_test_results <- vector(mode = "integer", length = length(folder))

for ( i in 1:length(folder)) {
  lrc <- read_6800(folder[i])
  
  PARlrc<-lrc$Qin #PAR (aka PPFD or Q)
  photolrc<-lrc$A #net photosynthetic rate (Anet)
  
  curvelrc<-data.frame(PARlrc,photolrc)
  curvelrc # *inspect raw data and check notebook (data reasonable or need edited/discarded?)
  
  par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
  plot(PARlrc,photolrc,xlab="", ylab="", ylim=c(-2,max(photolrc)+2),cex.lab=1.2,cex.axis=1.5,cex=2)
  mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1.5)
  mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
  
  # --- Nonlinear least squares regression (non-rectangular hyperbola).  4 parameter model: Amax (max gross photosytnthetic rate), Rd (dark respiration), AQY (apparent quantum yield), Theta (curvature parameter, dimensionless) ---
  # Another option is to fit AQY (initial slope) and Rd (y-intercept) separately using linear regression on data points that are not light-saturated, then use those model fits in the non-linear model to parameterize Amax and the curve parameter (theta). However, this requires user to subjectively decide which points are not light saturated (initial linear portion of curve). 
  # For more or Rd estimation see protocol text. 
  # Depending on data, quantile regression can be implemented through nlrq()
  
  #curve.nlslrc = nls(photolrc ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,start=list(Am=(max(photolrc)-min(photolrc)),AQY=0.05,Rd=-min(photolrc),theta=0.75)) 
  
  dat.df = data.frame(photolrc, PARlrc)
  
  curve.nlslrc = nls_multstart(photolrc ~ (1/(2*theta))*(AQY*PARlrc+Am-sqrt((AQY*PARlrc+Am)^2-4*AQY*theta*Am*PARlrc))-Rd,
                               iter = 1000,
                               data = dat.df,
                               start_lower = list(Am=(max(photolrc)-min(photolrc))*0.49,AQY=0.05*0.49,Rd=abs(min(photolrc))*-2,theta=0.75*0.49),
                               start_upper = list(Am=(max(photolrc)-min(photolrc))*2.1,AQY=0.05*2.1, Rd=abs(min(photolrc))*2,theta=0.75*2.1),
                               #lower = c(theta = 0, AQY=0, Am = 0, Rd=-10),
                               #upper = c(theta = 1, AQY = 10, Am = 10, Rd = 5),
                               supp_errors = "Y")
  
  #                     start=list(Am=(max(photolrc)-min(photolrc)),AQY=0.05,Rd=-min(photolrc),theta=0.75)) 
  
  
  summary(curve.nlslrc) #summary of model fit
  
  # ---Graph raw data with modeled curve---
  
  par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
  plot(PARlrc,photolrc,xlab="", ylab="", ylim=c(-2,max(photolrc)+2),cex.lab=1.2,cex.axis=1.5,cex=2)
  mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=2)
  mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2,cex=2)
  curve((1/(2*summary(curve.nlslrc)$coef["theta","Estimate"]))*(summary(curve.nlslrc)$coef["AQY","Estimate"]*x+summary(curve.nlslrc)$coef["Am","Estimate"]-sqrt((summary(curve.nlslrc)$coef["AQY","Estimate"]*x+summary(curve.nlslrc)$coef["Am","Estimate"])^2-4*summary(curve.nlslrc)$coef["AQY","Estimate"]*summary(curve.nlslrc)$coef["theta","Estimate"]*summary(curve.nlslrc)$coef["Am","Estimate"]*x))-summary(curve.nlslrc)$coef["Rd","Estimate"],lwd=2,col="blue",add=T)
  
  # ---Solve for light saturation point (LSP), PPFD where 75% of Amax is achieved (75% is arbitrary - cutoff could be changed) 
  #Change the value in this function to do 75% or 80% as needed
  x<-function(x) {(1/(2*summary(curve.nlslrc)$coef["theta","Estimate"]))*(summary(curve.nlslrc)$coef["AQY","Estimate"]*x+summary(curve.nlslrc)$coef["Am","Estimate"]-sqrt((summary(curve.nlslrc)$coef["AQY","Estimate"]*x+summary(curve.nlslrc)$coef["Am","Estimate"])^2-4*summary(curve.nlslrc)$coef["AQY","Estimate"]*summary(curve.nlslrc)$coef["theta","Estimate"]*summary(curve.nlslrc)$coef["Am","Estimate"]*x))-summary(curve.nlslrc)$coef["Rd","Estimate"]-(0.8*summary(curve.nlslrc)$coef["Am","Estimate"])+0.8*(summary(curve.nlslrc)$coef["Rd","Estimate"])}
  
  lrc_mean_test_results[i] <- uniroot(x,c(0,5000))$root
  view(lrc_mean_test_results)
  
  LSP_mean <- mean(lrc_mean_test_results)
  view(LSP_mean)  
} #LSP
