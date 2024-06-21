##############################
#NOT WORKING YET.......

##############################
#Species facet plot:
##############################
# Fit the models and include species
ww <- vec_mod_group %>% 
  lapply(function(sel_group) {
    tmp_dat <- dat[mod_group == sel_group]
    o <- tryCatch(
      {
        gam(photo ~ 
              s(fspecies, bs = 're') + 
              s(cond, k = 5, bs = 'ts') + 
              s(tleaf, k = 5, bs = 'ts'), 
            data = tmp_dat,
            select = TRUE,
            method = "REML")
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (is.null(o)) return(NULL)
    o$mod_group <- sel_group
    o$species <- unique(tmp_dat$species)[1] # Ensure single species
    return(o)
  })

# Remove failed models
ww <- ww[!sapply(ww, is.null)]

# Extract smooth estimates for `cond`
ww2 <- ww %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    smooth_estimates <- tryCatch(
      {
        gratia::smooth_estimates(x, unconditional = TRUE, freq = FALSE, select = 's(cond)', partial_match = TRUE)
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (is.null(smooth_estimates)) return(NULL)
    smooth_estimates %>%
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist(fill = TRUE)


# Extract smooth estimates for `tleaf`
vv2 <- ww %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    gratia::smooth_estimates(x, unconditional = TRUE, smooth = 's(tleaf)') %>% 
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')

# Extract smooth estimates for `cond` again, if needed
vv3 <- ww %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    gratia::smooth_estimates(x, unconditional = TRUE, smooth = 's(cond)') %>% 
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')

# Repeat for qq and zz models
qq <- vec_mod_group %>% 
  lapply(function(sel_group) {
    tmp_dat <- dat[mod_group == sel_group]
    o <- gam(photo ~ 
               s(fspecies, bs = 're') + 
               s(tleaf, bs = 'ts', k = 5), 
             data = tmp_dat,
             select = TRUE,
             method = "REML")
    o$mod_group <- sel_group
    o$form <- unique(tmp_dat$form)
    o$species <- unique(tmp_dat$species)[1] # Ensure single species
    return(o)
  })

qq2 <- qq %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    gratia::smooth_estimates(x, unconditional = TRUE, smooth = 's(tleaf)') %>% 
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = "m_t")

zz <- vec_mod_group %>% 
  lapply(function(sel_group) {
    tmp_dat <- dat[mod_group == sel_group]
    o <- gam(photo ~ 
               s(fspecies, bs = 're') + 
               s(cond, bs = 'ts', k = 5),
             data = tmp_dat,
             select = TRUE,
             method = "REML")
    o$mod_group <- sel_group
    o$form <- unique(tmp_dat$form)
    o$species <- unique(tmp_dat$species)[1] # Ensure single species
    return(o)
  })

zz2 <- zz %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    species <- x$species # Ensure single species
    gratia::smooth_estimates(x, unconditional = TRUE, smooth = 's(cond)') %>% 
      mutate(species = species,
             R2 = mod_r2,
             mod_group = mod_group)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_v')
# Extract Topt
fig_dat_topt <- 
  rbindlist(list(vv2, qq2), fill = TRUE) %>% 
  group_by(mod, species) %>% 
  filter(.estimate == max(.estimate)) %>% 
  filter(.estimate < 41) %>% 
  ungroup() %>% 
  select(mod, species, tleaf) %>% 
  rename(Topt = tleaf)

avg_weib_topt <- dat %>%
  group_by(species) %>%
  summarize(avg_weib_topt = mean(weib.topt, na.rm = TRUE)) %>%
  ungroup()

spec_photo <- rbindlist(list(vv2,qq2),fill=T)

# Create the updated plot
spec_photo_plot <- ggplot(spec_photo, aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_point() +
  geom_line(lwd = 1) +
  theme_classic() +
  geom_hline(aes(yintercept = 0), col = 'grey40', lwd = 0.5, lty = 1) + 
  geom_vline(data = fig_dat_topt, show.legend = FALSE, aes(xintercept = Topt, color = mod), lty = 2) +
  geom_ribbon(aes(ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = mod), color = NA, alpha = 0.25) +
  scale_color_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy"), labels = c("m_t" = "Temperature", "m_tv" = "Temperature + Conductance")) +
  scale_fill_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy"), labels = c("m_t" = "Temperature", "m_tv" = "Temperature + Conductance")) +
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  labs(color = 'GAM', fill = 'GAM', y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m**-2, s**-1, ")")), x = "Leaf temperature (°C)") +
  facet_wrap(~ species, scales = 'free_y', ncol = 2, shrink = TRUE) +
  theme(legend.position = 'bottom', strip.text = element_text(face = 'bold', hjust = 0)) +
  geom_vline(data = avg_weib_topt, show.legend = FALSE, aes(xintercept = avg_weib_topt), color = 'grey', lty = 2) +
  theme(legend.position = 'bottom', strip.text = element_text(face = 'bold', hjust = 0))

spec_photo_plot
