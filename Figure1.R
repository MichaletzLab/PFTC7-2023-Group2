# Make Figure 1 ----------------------------------------------------------------
# --- Helper functions for scaling A and E onto same axis ---
scale_E_to_A <- function(E) {
  rng_A <- range(raw.dat$A, na.rm = TRUE)
  rng_E <- range(raw.dat$E, na.rm = TRUE)
  (E - rng_E[1]) / diff(rng_E) * diff(rng_A) + rng_A[1]
}
scale_A_to_E <- function(A_scaled) {
  rng_A <- range(raw.dat$A, na.rm = TRUE)
  rng_E <- range(raw.dat$E, na.rm = TRUE)
  (A_scaled - rng_A[1]) / diff(rng_A) * diff(rng_E) + rng_E[1]
}

# --- Clean data: remove NAs ---
raw.dat <- raw.dat %>%
  filter(!is.na(E), !is.na(A), !is.na(gsw))

# --- Pivot A and E into long format ---
plot_dat <- raw.dat %>%
  pivot_longer(cols = c(A, E), names_to = "Flux", values_to = "Value") %>%
  mutate(
    Flux = recode(Flux,
                  A = "Photosynthesis",
                  E = "Transpiration"),
    Value_scaled = ifelse(Flux == "Transpiration",
                          scale_E_to_A(Value),  # put E on A axis
                          Value)
  )

# --- Fit GAMs for A and E separately ---
gam_A <- gam(A ~ s(Tleaf, k = 10), data = raw.dat)
gam_E <- gam(E ~ s(Tleaf, k = 10), data = raw.dat)

# Prediction grid
pred_grid <- data.frame(Tleaf = seq(min(raw.dat$Tleaf, na.rm=TRUE),
                                    max(raw.dat$Tleaf, na.rm=TRUE),
                                    length.out = 200))

pred_A <- predict(gam_A, newdata = pred_grid, se.fit = TRUE)
pred_E <- predict(gam_E, newdata = pred_grid, se.fit = TRUE)

# Put predictions into long format too
pred_grid_long <- bind_rows(
  tibble(
    Tleaf = pred_grid$Tleaf,
    Flux = "Photosynthesis",
    Pred_scaled = pred_A$fit,
    SE_scaled   = pred_A$se.fit
  ),
  tibble(
    Tleaf = pred_grid$Tleaf,
    Flux = "Transpiration",
    Pred_scaled = scale_E_to_A(pred_E$fit),
    SE_scaled   = pred_E$se.fit / diff(range(raw.dat$E, na.rm=TRUE)) *
      diff(range(raw.dat$A, na.rm=TRUE))
  )
)

# --- Plot combined A & E ---
A.E.full.plot <- ggplot() +
  # Raw points
  geom_point(data = plot_dat,
             aes(x = Tleaf, y = Value_scaled, color = Flux),
             alpha = 0.005, size = 2) +
  # Ribbons
  geom_ribbon(data = pred_grid_long,
              aes(x = Tleaf,
                  ymin = Pred_scaled - SE_scaled,
                  ymax = Pred_scaled + SE_scaled,
                  fill = Flux),
              alpha = 0.3, inherit.aes = FALSE) +
  # Lines
  geom_line(data = pred_grid_long,
            aes(x = Tleaf, y = Pred_scaled, color = Flux),
            size = 1.2) +
  scale_y_continuous(
    name = expression(A~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    sec.axis = sec_axis(~scale_A_to_E(.),
                        name = expression(E~(mmol~m^{-2}~s^{-1})))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis" = "green4", "Transpiration" = "blue3")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Photosynthesis" = "green4", "Transpiration" = "blue3")
  ) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  geom_vline(xintercept=TOPT, linetype="dashed", color="black") +
  labs(x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = "top")

# Plot gsw vs. Tleaf
gam_gsw <- gam(gsw ~ s(Tleaf, k = 10), data = raw.dat)

pred_grid <- data.frame(Tleaf = seq(min(raw.dat$Tleaf, na.rm=TRUE),
                                    max(raw.dat$Tleaf, na.rm=TRUE),
                                    length.out = 200))

pred_gsw <- predict(gam_gsw, newdata = pred_grid, se.fit = TRUE)
pred_grid$gsw_pred <- pred_gsw$fit
pred_grid$gsw_se <- pred_gsw$se.fit

gsw.plot <- ggplot(raw.dat, aes(x = Tleaf, y = gsw)) +
  geom_point(alpha = 0.01, size = 2, color = "turquoise") +
  geom_line(data = pred_grid, aes(x = Tleaf, y = gsw_pred), color = "turquoise", size = 1.2, inherit.aes = FALSE) +
  geom_ribbon(data = pred_grid, aes(x = Tleaf, ymin = gsw_pred - gsw_se, ymax = gsw_pred + gsw_se),
              fill = "turquoise", alpha = 0.3, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  geom_vline(xintercept=TOPT, linetype="dashed", color="black")+
  labs(x = "Leaf temperature (°C)",
       y = expression(gsw ~ (mol~m^-2~s^-1))) +
  theme_classic()

ggarrange(A.E.full.plot, gsw.plot, nrow=1, ncol=2, labels=c("A","B"), widths = c(4,3))

# Plot A and E vs. Tleaf at high temperatures (above Topt)
# --- Clean hot.dat (remove NAs) ---
hot.dat <- hot.dat %>%
  filter(!is.na(E), !is.na(A), !is.na(gsw))

# --- Pivot A and E into long format ---
plot_hot <- hot.dat %>%
  pivot_longer(cols = c(A, E), names_to = "Flux", values_to = "Value") %>%
  mutate(
    Flux = recode(Flux,
                  A = "Photosynthesis",
                  E = "Transpiration"),
    Value_scaled = ifelse(Flux == "Transpiration",
                          scale_E_to_A(Value),  # put E on A axis
                          Value)
  )

# --- Fit GAMs for A and E on hot subset ---
gam_A_hot <- gam(A ~ s(Tleaf, k = 10), data = hot.dat)
gam_E_hot <- gam(E ~ s(Tleaf, k = 10), data = hot.dat)

# Prediction grid for hot subset
pred_grid_hot <- data.frame(Tleaf = seq(min(hot.dat$Tleaf, na.rm=TRUE),
                                        max(hot.dat$Tleaf, na.rm=TRUE),
                                        length.out = 200))

pred_A_hot <- predict(gam_A_hot, newdata = pred_grid_hot, se.fit = TRUE)
pred_E_hot <- predict(gam_E_hot, newdata = pred_grid_hot, se.fit = TRUE)

# Put predictions into long format
pred_hot_long <- bind_rows(
  tibble(
    Tleaf = pred_grid_hot$Tleaf,
    Flux = "Photosynthesis",
    Pred_scaled = pred_A_hot$fit,
    SE_scaled   = pred_A_hot$se.fit
  ),
  tibble(
    Tleaf = pred_grid_hot$Tleaf,
    Flux = "Transpiration",
    Pred_scaled = scale_E_to_A(pred_E_hot$fit),
    SE_scaled   = pred_E_hot$se.fit / diff(range(hot.dat$E, na.rm=TRUE)) *
      diff(range(hot.dat$A, na.rm=TRUE))
  )
)

# --- Plot hot subset ---
A.E.hot.plot <- ggplot() +
  # Raw points
  geom_point(data = plot_hot,
             aes(x = Tleaf, y = Value_scaled, color = Flux),
             alpha = 0.005, size = 2) +
  # Ribbons
  geom_ribbon(data = pred_hot_long,
              aes(x = Tleaf,
                  ymin = Pred_scaled - SE_scaled,
                  ymax = Pred_scaled + SE_scaled,
                  fill = Flux),
              alpha = 0.3, inherit.aes = FALSE) +
  # Lines
  geom_line(data = pred_hot_long,
            aes(x = Tleaf, y = Pred_scaled, color = Flux),
            size = 1.2) +
  scale_y_continuous(
    name = expression(A~(mu*mol~CO[2]~m^{-2}~s^{-1})),
    sec.axis = sec_axis(~scale_A_to_E(.),
                        name = expression(E~(mmol~m^{-2}~s^{-1})))
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Photosynthesis" = "green4", "Transpiration" = "blue3")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Photosynthesis" = "green4", "Transpiration" = "blue3")
  ) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  geom_vline(xintercept=TOPT, linetype="dashed", color="black") +
  labs(x = "Leaf temperature (°C)") +
  theme_classic() +
  theme(legend.position = "top")



# Fit GAM for gsw
gam_gsw <- gam(gsw ~ s(Tleaf, k = 10), data = hot.dat)

pred_grid <- data.frame(Tleaf = seq(min(hot.dat$Tleaf, na.rm=TRUE),
                                    max(hot.dat$Tleaf, na.rm=TRUE),
                                    length.out = 200))

pred_gsw <- predict(gam_gsw, newdata = pred_grid, se.fit = TRUE)
pred_grid$gsw_pred <- pred_gsw$fit
pred_grid$gsw_se <- pred_gsw$se.fit

# Plot gsw vs. Tleaf
gsw.plot.hot <- ggplot(hot.dat, aes(x = Tleaf, y = gsw)) +
  geom_point(alpha = 0.01, size = 2, color = "turquoise") +
  geom_line(data = pred_grid, aes(x = Tleaf, y = gsw_pred), color = "turquoise", size = 1.2, inherit.aes = FALSE) +
  geom_ribbon(data = pred_grid, aes(x = Tleaf, ymin = gsw_pred - gsw_se, ymax = gsw_pred + gsw_se),
              fill = "turquoise", alpha = 0.3, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  geom_vline(xintercept=TOPT, linetype="dashed", color="black") +
  labs(x = "Leaf temperature (°C)",
       y = expression(gsw ~ (mol~m^-2~s^-1))) +
  theme_classic()

# Arrange plots side by side
hot.plot <- ggarrange(A.E.hot.plot, gsw.plot.hot, nrow = 1, ncol = 2, labels = c("A", "B"), widths = c(4,3))


deriv_A <- derivatives(gam_A, select = "s(Tleaf)", partial_match = TRUE)
deriv_E <- derivatives(gam_E, select = "s(Tleaf)", partial_match = TRUE)

ensure_ci_and_flag <- function(d){
  # find derivative column
  der_col <- intersect(names(d), c("derivative", "estimate", "deriv"))
  if (length(der_col) == 0) stop("No derivative column found")
  der_col <- der_col[1]
  
  # find se column
  se_col <- intersect(names(d), c("se", "se.derivative", "std.error"))
  if (length(se_col) == 0) stop("No standard error column found")
  se_col <- se_col[1]
  
  # add crit if missing
  if (!("crit" %in% names(d))) d$crit <- qnorm(0.975)
  
  # rename into standard names for consistency
  d <- d %>%
    rename(derivative = !!der_col, se = !!se_col)
  
  # create CI if missing
  if (!("lower" %in% names(d) && "upper" %in% names(d))) {
    d <- d %>%
      mutate(
        lower = derivative - crit * se,
        upper = derivative + crit * se
      )
  }
  
  d %>%
    mutate(
      sig = (lower > 0) | (upper < 0),
      dir = case_when(
        upper < 0 ~ "neg",
        lower > 0 ~ "pos",
        TRUE ~ "ns"
      )
    )
}
deriv_A <- deriv_A %>%
  mutate(
    sig = (.lower_ci > 0) | (.upper_ci < 0),
    dir = case_when(
      .upper_ci < 0 ~ "Negative",
      .lower_ci > 0 ~ "Positive",
      TRUE ~ "ns"))
deriv_E <- deriv_E %>%
  mutate(
    sig = (.lower_ci > 0) | (.upper_ci < 0),
    dir = case_when(
      .upper_ci < 0 ~ "Negative",
      .lower_ci > 0 ~ "Positive",
      TRUE ~ "ns"))
as.factor(deriv_A$dir)
as.factor(deriv_E$dir)
derivA_summary <- deriv_A %>%
  summarise(
    mean_derivative = mean(.derivative),
    mean_se = sqrt(sum(.se^2) / n()^2),
    n = n()) %>%
  mutate(
    t_value = mean_derivative / mean_se,
    p_value = 2 * (1 - pnorm(abs(t_value))),   # two-sided test
    p_value_fmt = format.pval(p_value, digits = 6, eps = .Machine$double.eps),
    ci_lower = mean_derivative - 1.96 * mean_se,
    ci_upper = mean_derivative + 1.96 * mean_se)
derivA_summary

derivE_summary <- deriv_E %>%
  summarise(
    mean_derivative = mean(.derivative),
    mean_se = sqrt(sum(.se^2) / n()^2),
    n = n()) %>%
  mutate(
    t_value = mean_derivative / mean_se,
    p_value = 2 * (1 - pnorm(abs(t_value))),   # two-sided test
    p_value_fmt = format.pval(p_value, digits = 6, eps = .Machine$double.eps),
    ci_lower = mean_derivative - 1.96 * mean_se,
    ci_upper = mean_derivative + 1.96 * mean_se)
derivE_summary

# Add dummy row for Positive slope to fix the legend
deriv_A_plot_data <- deriv_A %>%
  bind_rows(tibble(Tleaf = NA,.derivative = NA,.lower_ci = NA,.upper_ci = NA,sig = TRUE,dir = factor("Positive", levels = c("Negative","Positive"))))

deriv_A.plot <- ggplot(deriv_A_plot_data, aes(x = Tleaf, y = .derivative)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = filter(deriv_A_plot_data, sig), aes(color = dir), size = 2) +
  scale_color_manual(
    name = "Slope direction",
    values = c("Negative" = "red3", "Positive" = "lightblue3")
  ) +
  labs(
    x = "Leaf temperature (°C)",
    y = expression(A * "'" ~ "(" * mu*mol~CO[2]~m^-2~s^-1~degree*C^-1 * ")")
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title.y = element_text(size = 10)
  ) +
  annotate(
    "text",
    x = min(deriv_A$Tleaf, na.rm = TRUE) + 1,
    y = -0.75, 
    label = paste0("Mean slope = ", round(derivA_summary$mean_derivative, 3),
                   "\n95% CI: [", round(derivA_summary$ci_lower, 3),
                   ", ", round(derivA_summary$ci_upper, 3), "]"),
    hjust = 0, size = 4
  )

# Derivative plot for E
deriv_E.plot <- ggplot(deriv_E, aes(x = Tleaf, y = .derivative)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = filter(deriv_E, sig), aes(color = dir), size = 2) +
  scale_color_manual(
    name = "Slope direction",
    values = c("Negative" = "red3", "Positive" = "lightblue3"),
    drop = FALSE
  ) +
  labs(
    x = "Leaf temperature (°C)",
    y = expression(E * "'" ~ "(" * mmol~m^-2~s^-1~degree*C^-1 * ")")
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title.y = element_text(size = 10)
  ) +
  annotate(
    "text",
    x = min(deriv_E$Tleaf, na.rm = TRUE) + 1,
    y = -0.0003, 
    label = paste0(
      "Mean slope = ", format(round(derivE_summary$mean_derivative, 5), scientific = FALSE),
      "\n95% CI: [", format(round(derivE_summary$ci_lower, 5), scientific = FALSE),
      ", ", format(round(derivE_summary$ci_upper, 5), scientific = FALSE), "]"
    ),
    hjust = 0, size = 4
  )

deriv_gsw <- derivatives(gam_gsw, type = "central") %>%
  mutate(sig = (.lower_ci > 0) | (.upper_ci < 0),
         dir = case_when(
           .upper_ci < 0 ~ "Negative",
           .lower_ci > 0 ~ "Positive",
           TRUE ~ "ns"
         ))

# Summary for mean slope
deriv_gsw_summary <- deriv_gsw %>%
  summarise(mean_derivative = mean(.derivative),
            mean_se = sqrt(sum(.se^2) / n())) %>%
  mutate(
    t_value = mean_derivative / mean_se,
    p_value = 2 * (1 - pnorm(abs(t_value))),
    ci_lower = mean_derivative - 1.96 * mean_se,
    ci_upper = mean_derivative + 1.96 * mean_se
  )

# Plot 2: Derivative of gsw vs Tleaf
gsw_deriv.plot <- ggplot(deriv_gsw, aes(x = Tleaf, y = .derivative)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = filter(deriv_gsw, sig), aes(color = dir), size = 2) +
  scale_color_manual(name = "Slope direction",
                     values = c("Negative" = "red3", "Positive" = "lightblue3")) +
  labs(x = "Leaf temperature (°C)",
       y = expression(gsw * "'" ~ "(" * mol~m^-2~s^-1~degree*C^-1 * ")"))+
  theme_classic() +
  theme(legend.position = "right") +
  annotate(
    "text",
    x = min(deriv_gsw$Tleaf, na.rm = TRUE) + 1,
    y = -0.015,
    label = paste0("Mean slope = ", signif(deriv_gsw_summary$mean_derivative, 4),
                   "\n95% CI: [", signif(deriv_gsw_summary$ci_lower, 4),
                   ", ", signif(deriv_gsw_summary$ci_upper, 4), "]"),
    hjust = 0,
    size = 4
  )

hot.plot
ggarrange(deriv_A.plot, deriv_E.plot, gsw_deriv.plot, nrow=1, common.legend = TRUE,labels = c("C", "D","E"),legend = "right")
