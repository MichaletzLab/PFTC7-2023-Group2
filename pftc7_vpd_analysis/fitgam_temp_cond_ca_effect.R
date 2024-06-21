# Multipanel site plot ===============================================
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
    o$site <- unique(tmp_dat$site)
    return(o)
  })

# Remove failed models
ww <- ww[!sapply(ww, is.null)]

# Extract smooth estimates
ww2 <- ww %>% 
  lapply(function(x) {
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    site <- x$site
    smooth_estimates <- tryCatch(
      {
        gratia::smooth_estimates(x, unconditional = TRUE, freq = FALSE, select = 's(cond', partial_match = TRUE)
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (is.null(smooth_estimates)) return(NULL)
    smooth_estimates %>%
      mutate(species = unique(x$model$species),
             R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist(fill = TRUE)
vv2 <- ww %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(tleaf)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')
vv2$R2 %>% summary


vv3 <- ww %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(cond)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tv')
vv3$R2 %>% summary


qq <- vec_mod_group %>% 
  lapply(., FUN = function(sel_group){
    tmp_dat <- dat[mod_group ==sel_group]
    o <- gam(photo ~ 
               s(fspecies,bs='re') + 
               s(tleaf,bs='ts',k=5), 
             data=tmp_dat,
             select=T,
             method = "REML")
    o$mod_group <- sel_group
    o$form <- unique(tmp_dat$form)
    o$site <- unique(tmp_dat$site)
    return(o)
  })

qq2 <- qq %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(tleaf)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = "m_t")
ww2$R2 %>% summary


gam_medlyn <- vec_mod_group %>% 
  lapply(function(sel_group) {
    tmp_dat <- dat[mod_group == sel_group]
    o <- tryCatch(
      {
        gam(photo ~ 
              s(fspecies, bs = 're') + 
              s(cond, k = 5, bs = 'ts') + 
              s(cica, k = 5, bs = 'ts') +
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
    o$site <- unique(tmp_dat$site)
    return(o)
  })
gam_medlyn_sm <- gam_medlyn %>% 
  lapply(., FUN = function(x){
    mod_r2 <- summary(x)$r.sq
    mod_group <- x$mod_group
    form <- x$form
    site <- x$site
    gratia::smooth_estimates(x,unconditional=T,smooth='s(tleaf)') %>% 
      mutate(species = x$species) %>% 
      mutate(R2 = mod_r2,
             mod_group = mod_group,
             site = site)
  }) %>% 
  rbindlist() %>% 
  mutate(mod = 'm_tgc')

# Extract Topt ==================
fig_dat_topt <- 
  rbindlist(list(vv2,qq2,gam_medlyn_sm),fill=T) %>% 
  group_by(mod,site) %>% 
  filter(.estimate == max(.estimate)) %>% 
  filter(.estimate < 41) %>% 
  ungroup() %>% 
  select(mod,site,tleaf) %>% 
  rename(Topt = tleaf)

avg_weib_topt <- dat %>%
  group_by(site) %>%
  summarize(avg_weib_topt = mean(weib.topt, na.rm = TRUE)) %>%
  ungroup()

## Temperature Response Figure =============
site_levels <- c("S.Africa 1","Norway 1", "S.Africa 2","Norway 2", "S.Africa 3","Norway 3","S.Africa 4", "Norway 4","S.Africa 5")
p_photo <- rbindlist(list(vv2,qq2),fill=T)
p_photo_plot <- ggplot(p_photo, aes(x = tleaf, y = .estimate, color = mod, fill = mod)) +
  geom_point() +
  geom_line(lwd = 1) +
  theme_classic() +
  geom_hline(aes(yintercept = 0), col = 'grey40', lwd = 0.5, lty = 1) + 
  geom_vline(data = fig_dat_topt, show.legend = FALSE, aes(xintercept = Topt, color = mod), lty = 2) +
  geom_ribbon(aes(ymin = .estimate - 1.96 * .se, ymax = .estimate + 1.96 * .se, fill = mod), color = NA, alpha = 0.25) +
  scale_color_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy"), labels = c("m_t" = "Temperature", "m_tv" = "Temperature + Conductance", "m_tgc" = "Temperature + Conductance + Ci/Ca")) +
  scale_fill_manual(values = c("m_t" = "#cf0000", "m_tv" = "navy"), labels = c("m_t" = "Temperature", "m_tv" = "Temperature + Conductance", "m_tgc" = "Temperature + Conductance + Ca")) +
  scale_x_continuous(breaks = seq(5, 40, 5), expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(color = 'GAM', fill = 'GAM', y = expression(paste("(Partial) effect on Photosynthesis (µmol ", CO[2], " ", m**-2, s**-1,")")), x = "Leaf temperature (°C)") +
  facet_wrap(~ factor(site, levels = c("1", "6", "2", "7", "3", "8", "4", "9", "5"),
                      labels = site_levels),
             scales = 'free_y', ncol = 2, shrink = TRUE) +
  theme(legend.position = 'bottom', strip.text = element_text(face='bold', hjust = 0)) +
  geom_vline(data = avg_weib_topt, show.legend = FALSE, aes(xintercept = avg_weib_topt), color = 'grey', lty = 2) +
  theme(legend.position = 'bottom', strip.text = element_text(face='bold', hjust = 0))
p_photo_plot

