# Use mgcv analysis package to ...
# a) Quantify interaction effects
# b) Test for Concurvity
# ====================================================
# The following code for fitting GAM models
# was adapted from code written by Martijn Slot.
# Source:
#   Slot M, Winter K. 2017.
#   New Phyt doi: 10.1111/nph.14469
#
# ====================================================

# Fit GAM with interaction terms for data across sites
gams_across_sites <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group) {
    tmp_dat <- dat[mod_group == sel_group]
    o <- gam(photo ~ 
               s(fspecies, bs = 're') + 
               s(vpdl, bs = 'ts', k = 5) + 
               s(tleaf, bs = 'ts', k = 5),
               # TODO: Spline fit for interaction effects
               # ti(tleaf, vpdl, bs = 'ts', k = 5),
             data = tmp_dat,
             select = TRUE,
             method = "REML")
    o$mod_group <- sel_group
    o$site <- unique(tmp_dat$site)
    return(o)
  })

# Function to extract R2 and interaction effects
# For a specified smoothing term (e.g. vpdl or tleaf)
extract_effects <- function(model_list, sm_term) {
  model_list %>% 
    lapply(., FUN = function(x) {
      mod_r2 <- summary(x)$r.sq
      mod_group <- x$mod_group
      site <- x$site
      gratia::smooth_estimates(x, unconditional = TRUE, smooth = sm_term) %>% 
        mutate(species = x$species) %>% 
        mutate(R2 = mod_r2,
               mod_group = mod_group,
               site = site)
    }) %>% 
    rbindlist() %>% 
    mutate(mod = 'm_tv')
}

# Extract effects for temperature, VPD, and their interaction
R2_temp_effects <- extract_effects(gams_across_sites, 's(tleaf)')
R2_cica_effects <- extract_effects(gams_across_sites, 's(vpdl)')
#R2_interaction_effects <- extract_effects(gams_across_sites, 'ti(tleaf,vpdl)')

# Summarize R2 values
summary(R2_temp_effects$R2)
summary(R2_cica_effects$R2)
# summary(R2_interaction_effects$R2)

# Concurvity analysis for each fitted GAM
concurvity_results <- gams_across_sites %>% 
  lapply(., FUN = function(x) {
    mod_group <- x$mod_group
    site <- x$site
    conc <- concurvity(x)
    list(mod_group = mod_group, site = site, concurvity = conc)
  })

# Report concurvity results for the model parameters
for (mod_test in concurvity_results) {
  print(mod_test$concurvity)
}
