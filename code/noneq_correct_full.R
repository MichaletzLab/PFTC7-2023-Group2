# Function for computing dynamic assimilation using the full dynamic equation
# explained in the Saathof and Welles paper, and in the new LI-6800 manual

# Arguments:
# Data is a dataframe with LI-6800 data
# dt1 is the calibration "tuning" constant from licor for REF
# dt2 is the calibration "tuning" constant from licor for SAMPLE
# aV is the effective volume tuning constant
noneq_correct_full = function(data, dt1_c, dt2_c, aV_c, dt1_h, dt2_h, aV_h) {
  
  # 
  return_data = c()
  
  # Iterate over curves in dataset
  for (i in unique(data$curveID)) {
    
    # Grab first curve
    non.eq.AT = subset(data, curveID == i)
    
    # Get length
    len = dim(non.eq.AT)[1]
    
    # Get CO2 concentrations
    co = non.eq.AT$CO2_s/1e6 #[2:n] #umol/mol; convert to mol/mol
    ce = non.eq.AT$CO2_r/1e6 #[2:n] #umol/mol; convert to mol/mol
    
    # Get H2O concentrations
    wo = non.eq.AT$H2O_s/1e3 #[2:n] #mmol/mol; convert to mol/mol
    we = non.eq.AT$H2O_r/1e3 #[2:n] #mmol/mol; convert to mol/mol
    
    # 1. Get Crd and Csd at each t - effective "dry" CO2 concentration
    Crd = ce/(1-we)*1e6 # factor of 1e6 converts back to umol/mol
    Csd = co/(1-wo)*1e6
    
    Wr = non.eq.AT$H2O_r # mmol/mol
    Ws = non.eq.AT$H2O_s # mmol/mol
    
    # 2. Interpolate Crd at -dt1
      # Find the two points which straddle time-dt1 in time
      # If one of them matches exactly, use that one
      # Otherwise, interpolate between values
    
    time = non.eq.AT$time
    past_time = time-dt1_c
    
    Crd_dt1 = rep(NA, len)
    
    for (i in 1:len) {
      # check if past time < beginning of time; if so, Crd_dt1 = NA
      if (past_time[i] < time[1]) {
        Crd_dt1[i] = NA
        next
      }
      
      # check if past time matches exactly one time; if so, set Crd_dt1 to that value
      if (sum(past_time[i] == time) == 1) {
        Crd_dt1[i] = Crd[which(past_time[i] == time)]
        next
      }
      
      # otherwise, find the two straddling points in time
      diff = time-past_time[i]
      next_obs = which(min(diff[diff > 0]) == diff)
      last_obs = which(max(diff[diff < 0]) == diff)
      
      time_diff = past_time[i]-time[last_obs]
      slope = (Crd[next_obs]-Crd[last_obs])/(time[next_obs]-time[last_obs])
      Crd_last = Crd[last_obs]
      
      Crd_dt1[i] = Crd_last + slope*time_diff
      
    }
    
    # 3. Compute dCsd/dt at each t
    dCsc_dt = (Csd[2:len]-Csd[1:len-1])/(time[2:len]-time[1:len-1])
    dCsc_dt = c(NA, dCsc_dt)
    
    dCsc_dt_dt2 = rep(NA, len)
    
    past_time = time-dt2_c
    
    # 4. Interpolate dCsd/dt at -dt2
    for (i in 1:len) {
      # check if past time < beginning of time; if so, Crd_dt1 = NA
      if (past_time[i] < time[1]) {
        dCsc_dt_dt2[i] = NA
        next
      }
      
      # check if past time matches exactly one time; if so, set Crd_dt1 to that value
      if (sum(past_time[i] == time) == 1) {
        dCsc_dt_dt2[i] = dCsc_dt[which(past_time[i] == time)]
        next
      }
      
      # otherwise, find the two straddling points in time
      diff = time-past_time[i]
      next_obs = which(min(diff[diff > 0]) == diff)
      last_obs = which(max(diff[diff < 0]) == diff)
      
      time_diff = past_time[i]-time[last_obs]
      slope = (dCsc_dt[next_obs]-dCsc_dt[last_obs])/(time[next_obs]-time[last_obs])
      dCsc_dt_last = dCsc_dt[last_obs]
      
      dCsc_dt_dt2[i] = dCsc_dt_last + slope*time_diff
      
    }
    
    # 4.1 Interpolate Wr at -dt1
    past_time = time-dt1_h
    Wr_dt1 = rep(NA, len)
    
    for (i in 1:len) {
      # check if past time < beginning of time; if so, Crd_dt1 = NA
      if (past_time[i] < time[1]) {
        Wr_dt1[i] = NA
        next
      }
      # check if past time matches exactly one time; if so, set Crd_dt1 to that value
      if (sum(past_time[i] == time) == 1) {
        Wr_dt1[i] = Wr[which(past_time[i] == time)]
        next
      }
      
      # otherwise, find the two straddling points in time
      diff = time-past_time[i]
      next_obs = which(min(diff[diff > 0]) == diff)
      last_obs = which(max(diff[diff < 0]) == diff)
      
      time_diff = past_time[i]-time[last_obs]
      slope = (Wr[next_obs]-Wr[last_obs])/(time[next_obs]-time[last_obs])
      Wr_last = Wr[last_obs]
      
      Wr_dt1[i] = Wr_last + slope*time_diff
    }
    
    # 4.2 Interpolate dWs/dt at -dt2
    # Compute dWs/dt at each t
    dWs_dt = (Ws[2:len]-Ws[1:len-1])/(time[2:len]-time[1:len-1])
    dWs_dt = c(NA, dWs_dt)
    
    dWs_dt_dt2 = rep(NA, len)
    past_time = time-dt2_h
    
    # Interpolate dWs/dt at -dt2
    for (i in 1:len) {
      # check if past time < beginning of time; if so, Crd_dt1 = NA
      if (past_time[i] < time[1]) {
        dWs_dt_dt2[i] = NA
        next
      }
      # check if past time matches exactly one time; if so, set Crd_dt1 to that value
      if (sum(past_time[i] == time) == 1) {
        dWs_dt_dt2[i] = dWs_dt[which(past_time[i] == time)]
        next
      }
      
      # otherwise, find the two straddling points in time
      diff = time-past_time[i]
      next_obs = which(min(diff[diff > 0]) == diff)
      last_obs = which(max(diff[diff < 0]) == diff)
      
      time_diff = past_time[i]-time[last_obs]
      slope = (dWs_dt[next_obs]-dWs_dt[last_obs])/(time[next_obs]-time[last_obs])
      dWs_dt_last = dWs_dt[last_obs]
      
      dWs_dt_dt2[i] = dWs_dt_last + slope*time_diff
      
    }
    
    # 5. Finally, compute Adyn and Edyn
    
    S = non.eq.AT$S #[2:n] # cm^2
    Tair = non.eq.AT$Tair
    Flow = non.eq.AT$Flow
    Pa = non.eq.AT$Pa
    #H2OR = non.eq.AT$H2O_r
    
    Adyn = (Crd_dt1 - Csd - Pa*1000/(8.314*(Tair+273.15)) * aV_c/Flow * dCsc_dt_dt2) * Flow/(100*S) * (1000-Wr)/1000
    Edyn = (1/1000)*(Ws - Wr_dt1 + Pa*1000/(8.314*(Tair+273.15)) * aV_h/Flow * dWs_dt_dt2) * Flow/(100*S) * 1000/(1000-Ws)
    # Note the additional leading factor of 1/1000 here is to convert to mol/m2s
    
    non.eq.AT$A = Adyn
    non.eq.AT$E = Edyn
    
    # Also need to recalculate gsw and Ci based on new E
    non.eq.AT = non.eq.AT %>% mutate(gtw = E*(1000-(1000*0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd))/(Pa+`ΔPcham`)+H2O_s)/2)/(1000*0.61365*exp(17.502*TleafCnd/(240.97+TleafCnd))/(Pa+`ΔPcham`)-H2O_s),
                       gsw = 2 / ((1/gtw - 1/gbw) + sign(gtw)*sqrt((1/gtw - 1/gbw)^2 + 4*K/((K+1)^2)*(2/(gtw*gbw) - 1/(gbw^2)))),
                       gtc = 1 / ((K+1)/(gsw/1.6)+1/(gbw/1.37)) + K/((K+1)/(gsw/1.6) + K/(gbw/1.37)),
                       Ci = ((gtc-E/2)*Ca-A)/(gtc+E/2))
    
    return_data = bind_rows(return_data, non.eq.AT)
  }
  
  # Drop any with A = NA or E = NA
  return_data = subset(return_data, !is.na(A))
  return_data = subset(return_data, !is.na(E))
  
  return(return_data)
  
}
