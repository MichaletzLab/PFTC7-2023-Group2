# Stage 1: Shape classification (peaked/unimodal or monotonic) of gas exchange 
# temperature responses, for each curve. Uses an interior-optimum gate, then
# AIC within peaked or monotonic classes.
#
# The gate is used because with such high density data for each curve, AIC 
# favors flexibility, and the Sharpe-Schoolfield (SS) can exhibit a monotonic 
# response by placing T_opt outside the range of observed T_leaf values, such
# that SS is selected by AIC even for monotonic responses. The definition of 
# "peaked" is a resolved interior optimum, so we gate peaked vs monotonic on 
# whether the estiamted SS T_opt falls within [T_min + MARGIN, T_max - MARGIN]. 
# AIC is then used only within a class.
#
# Two criteria, in sequence, per curve:
#   (1) GATE : SS T_opt interior?  yes -> peaked ; no -> monotonic
#   (2) AIC  : monotonic -> linear vs exp ; peaked -> SS vs Arroyo
# SS is retained as the peaked estimator for interpretability (explicit T_opt,
# E_a, E_d, Theta) regardless of the SS-vs-Arroyo AIC outcome; that comparison
# exists to test the Methods claim that SS is preferred over Arroyo.
#
# Models (all natural scale, Gaussian likelihood, anchored at T_ref = 25 C):
#   linear (2p, monotonic): y = b0 + b1*(Tleaf - 25)
#   exp    (2p, monotonic): y = a*exp(b*(Tleaf - 25))
#   SS     (4p, peaked)   : Sharpe-Schoolfield, explicit T_opt (pipeline form)
#   Arroyo (3p, peaked)   : Arroyo et al. 2022 transition-state model, Eq. 6/7
#           with alpha = 0, reparameterized to the centered T_opt-explicit form
#           from their SI (Eq. S34): y = Y_opt*(T_opt/T)^(-a)*exp(a*(T_opt/T-1)),
#           T in Kelvin, a < 0 for a peak. Centering avoids the under/overflow of
#           the raw 1/T form and yields an explicit T_opt. (alpha is absorbed
#           into the free exponent, so Eq. 6 vs 7 changes neither fit nor AIC.)
#
# Variables: A, gsw, E, iWUE (=A/gsw), eWUE (i.e., WUE=A/E; named eWUE to follow
# earlier code).
# ============================================================================

suppressMessages({ library(dplyr); library(nls.multstart) })

# ----------------------------- user settings --------------------------------
INPUT_CSV    <- "data/raw.discardHooks_data.csv"
OUT_PERCURVE <- file.path("outputs", "model_selection_percurve.csv")
OUT_CLASS    <- file.path("outputs", "model_selection_classification.csv")
OUT_WITHIN   <- file.path("outputs", "model_selection_withinfamily.csv")
T_REF        <- 25     # anchor for centering (C)
N_ITER       <- 500    # nls_multstart random starts (drop to 250 if slow)
MARGIN       <- 4      # interior-optimum margin (C); classifications invariant over 2-5 C
PEAK_FRAC    <- 0.5    # variable is called peaked if >= this fraction of curves are peaked
MIN_POINTS   <- 6
VARS         <- c("A", "gsw", "E", "iWUE", "eWUE")   # eWUE = WUE as per earlier code
set.seed(1)
dir.create("outputs", showWarnings = FALSE)

# ----------------------------- model functions ------------------------------
mod_exp <- function(Tleaf, a, b, T_ref = T_REF) a * exp(b * (Tleaf - T_ref))

mod_ss <- function(Tleaf, J_ref, E, E_D, T_opt, T_ref = T_REF) {
  k <- 8.62e-5; Tk <- Tleaf + 273.15; Trefk <- T_ref + 273.15; Toptk <- T_opt + 273.15
  J_ref * exp(E * (1/(k*Trefk) - 1/(k*Tk))) /
    (1 + E/(E_D - E) * exp((E_D/k) * (1/Toptk - 1/Tk)))
}

# Arroyo, centered on its own T_opt (a < 0 gives a peak). No T_ref needed.
mod_arroyo <- function(Tleaf, Y_opt, T_opt, a) {
  Tk <- Tleaf + 273.15; Toptk <- T_opt + 273.15
  y <- Toptk / Tk
  Y_opt * y^(-a) * exp(a * (y - 1))
}

# ----------------------------- record helper --------------------------------
# One standard row per fit. T_opt and 'interior' are populated only for models
# that carry an explicit T_opt parameter (SS, Arroyo).
make_record <- function(curveID, variable, model, family, fit, y, pred, Tmin, Tmax) {
  base <- data.frame(curveID = curveID, variable = variable, model = model,
                     family = family, n = length(y), Tmin = Tmin, Tmax = Tmax,
                     k = NA_real_, logLik = NA_real_, AIC = NA_real_, r_sq = NA_real_,
                     converged = FALSE, T_opt = NA_real_, interior = NA,
                     p1 = NA_real_, p2 = NA_real_, p3 = NA_real_, p4 = NA_real_,
                     stringsAsFactors = FALSE)
  if (is.null(fit)) return(base)
  cf   <- coef(fit)
  ll   <- as.numeric(logLik(fit)); kk <- attr(logLik(fit), "df"); aic <- AIC(fit)
  rss  <- sum((y - pred)^2); tss <- sum((y - mean(y))^2)
  topt <- if ("T_opt" %in% names(cf)) as.numeric(cf["T_opt"]) else NA_real_
  base$k <- kk; base$logLik <- ll; base$AIC <- aic
  base$r_sq <- if (tss > 0) 1 - rss/tss else NA_real_
  base$converged <- TRUE
  base$T_opt <- topt
  base$interior <- if (!is.na(topt)) (topt >= Tmin + MARGIN & topt <= Tmax - MARGIN) else NA
  if (length(cf) >= 1) base$p1 <- cf[[1]]
  if (length(cf) >= 2) base$p2 <- cf[[2]]
  if (length(cf) >= 3) base$p3 <- cf[[3]]
  if (length(cf) >= 4) base$p4 <- cf[[4]]
  base
}

# ----------------------------- per-curve fitter -----------------------------
fit_one_curve <- function(curveID, variable, Tleaf, y) {
  ok <- is.finite(Tleaf) & is.finite(y); Tleaf <- Tleaf[ok]; y <- y[ok]
  Tmin <- if (length(Tleaf)) min(Tleaf) else NA_real_
  Tmax <- if (length(Tleaf)) max(Tleaf) else NA_real_
  fam  <- list(c("linear","monotonic"), c("exp","monotonic"),
               c("SS","peaked"), c("arroyo","peaked"))
  if (length(y) < MIN_POINTS)
    return(do.call(rbind, lapply(fam, function(mf)
      make_record(curveID, variable, mf[1], mf[2], NULL, y, NULL, Tmin, Tmax))))
  d <- data.frame(Tleaf = Tleaf, y = y); ymax <- max(abs(y))

  lin <- tryCatch({
    f <- lm(y ~ I(Tleaf - T_REF), data = d)
    make_record(curveID, variable, "linear", "monotonic", f, y, fitted(f), Tmin, Tmax)
  }, error = function(e) make_record(curveID, variable, "linear", "monotonic", NULL, y, NULL, Tmin, Tmax))

  ex <- tryCatch({
    f <- nls_multstart(y ~ mod_exp(Tleaf, a, b), data = d, iter = N_ITER,
           start_lower = c(a = -ymax, b = -0.5),
           start_upper = c(a = ymax + 1e-6, b = 0.5),
           supp_errors = 'Y', na.action = na.omit)
    if (is.null(f)) make_record(curveID, variable, "exp", "monotonic", NULL, y, NULL, Tmin, Tmax)
    else make_record(curveID, variable, "exp", "monotonic", f, y, fitted(f), Tmin, Tmax)
  }, error = function(e) make_record(curveID, variable, "exp", "monotonic", NULL, y, NULL, Tmin, Tmax))

  ss <- tryCatch({
    f <- nls_multstart(y ~ mod_ss(Tleaf, J_ref, E, E_D, T_opt), data = d, iter = N_ITER,
           start_lower = c(J_ref = 0, E = 0, E_D = 0.2, T_opt = 0),
           start_upper = c(J_ref = 2*ymax + 1e-6, E = 2, E_D = 5, T_opt = 50),
           lower = c(J_ref = -2*ymax, E = 0, E_D = 0, T_opt = 0),
           supp_errors = 'Y', na.action = na.omit)
    if (is.null(f)) make_record(curveID, variable, "SS", "peaked", NULL, y, NULL, Tmin, Tmax)
    else make_record(curveID, variable, "SS", "peaked", f, y, fitted(f), Tmin, Tmax)
  }, error = function(e) make_record(curveID, variable, "SS", "peaked", NULL, y, NULL, Tmin, Tmax))

  ar <- tryCatch({
    f <- nls_multstart(y ~ mod_arroyo(Tleaf, Y_opt, T_opt, a), data = d, iter = N_ITER,
           start_lower = c(Y_opt = 0, T_opt = Tmin, a = -3000),
           start_upper = c(Y_opt = 2*ymax + 1e-6, T_opt = Tmax, a = -1),
           supp_errors = 'Y', na.action = na.omit)
    if (is.null(f)) make_record(curveID, variable, "arroyo", "peaked", NULL, y, NULL, Tmin, Tmax)
    else make_record(curveID, variable, "arroyo", "peaked", f, y, fitted(f), Tmin, Tmax)
  }, error = function(e) make_record(curveID, variable, "arroyo", "peaked", NULL, y, NULL, Tmin, Tmax))

  rbind(lin, ex, ss, ar)
}

# ----------------------------- run contest ----------------------------------
message("Reading: ", INPUT_CSV)
dat <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)

# Whole-curve data quality exclusions
dat <- dat %>% filter(curveID != 63) # exclude curve 63 (has impossible gsw below ~29 C)
valid_curves <- dat %>%
  group_by(curveID) %>%
  summarise(has_valid = any(gsw < 0.6 & gsw > -1, na.rm = TRUE), .groups = "drop") %>% # exclude curves with no gsw value bewteen -1 and 0.6 (curveID 2; following DataPrep.R)
  filter(has_valid) %>%
  pull(curveID)
dat <- dat %>% filter(curveID %in% valid_curves)

dat <- dat %>% mutate(iWUE = A / gsw, eWUE = A / E)
stopifnot(all(c("Tleaf", "curveID", VARS) %in% names(dat)))

records <- list()
for (v in VARS) {
  message("Fitting variable: ", v)
  for (cid in sort(unique(dat$curveID))) {
    cd <- dat[dat$curveID == cid, ]
    records[[length(records) + 1]] <- fit_one_curve(cid, v, cd$Tleaf, cd[[v]])
  }
}
percurve <- do.call(rbind, records); rownames(percurve) <- NULL
write.csv(percurve, OUT_PERCURVE, row.names = FALSE)
message("Wrote per-curve fits: ", OUT_PERCURVE)

# ---------------------- Figure 4 SS-parameter export (A, iWUE) ----------------------
# Interpretable SS parameters for the PEAKED variables only, from the corrected
# (data-adaptive-bound) fits above. Theta = width of the fitted SS curve at
# THETA_THRESH of its maximum. Set THETA_THRESH to match the manuscript's Theta.
THETA_THRESH <- 0.95   # 0.95 matches fit_mod_schoolfield; use 0.80 if that is your Theta

ss_breadth <- function(J_ref, E, E_D, T_opt, Tmin, Tmax, thresh = THETA_THRESH, n = 2001) {
  if (any(is.na(c(J_ref, E, E_D, T_opt, Tmin, Tmax)))) return(NA_real_)
  Tg <- seq(Tmin, Tmax, length.out = n)
  y  <- mod_ss(Tg, J_ref, E, E_D, T_opt)          # T_ref defaults to T_REF (25)
  if (!any(is.finite(y))) return(NA_real_)
  ipk <- which.max(y); target <- thresh * y[ipk]
  yasc <- y[1:ipk];          Tasc <- Tg[1:ipk]
  ydes <- y[ipk:length(Tg)]; Tdes <- Tg[ipk:length(Tg)]
  Tlo <- if (min(yasc) < target) approx(yasc, Tasc, xout = target)$y else NA_real_
  Thi <- if (min(ydes) < target) approx(ydes, Tdes, xout = target)$y else NA_real_
  if (is.na(Tlo) || is.na(Thi)) NA_real_ else Thi - Tlo
}

sscols <- percurve[percurve$model == "SS" &
                     percurve$variable %in% c("A", "iWUE") &
                     percurve$converged, ]
fig4 <- data.frame(
  curveID  = sscols$curveID, variable = sscols$variable,
  J_ref = sscols$p1, Ea = sscols$p2, Ed = sscols$p3, Topt = sscols$T_opt,  # p1..p4 = J_ref,E,E_D,T_opt
  Tmin = sscols$Tmin, Tmax = sscols$Tmax, interior = sscols$interior,
  stringsAsFactors = FALSE)
fig4$Theta <- mapply(ss_breadth, fig4$J_ref, fig4$Ea, fig4$Ed, fig4$Topt, fig4$Tmin, fig4$Tmax)
out_ss <- file.path("outputs", "figure4_ss_parameters.csv")
write.csv(fig4, out_ss, row.names = FALSE)
message("Wrote Figure 4 SS parameters (A, iWUE): ", out_ss)

# --------------------- within-family AIC comparison helper ------------------
# Restricted to `ids` and to curves where BOTH `models` converged.
within_compare <- function(sub, ids, models) {
  s <- sub[sub$curveID %in% ids & sub$model %in% models,
           c("curveID", "model", "AIC", "converged")]
  ok <- s %>% group_by(curveID) %>%
    summarise(both = sum(converged) == length(models), .groups = "drop") %>%
    filter(both)
  s <- s[s$curveID %in% ok$curveID, ]
  n <- length(unique(s$curveID))
  if (n == 0)
    return(data.frame(model = models, n_curves = 0L, sumAIC = NA_real_,
                      mean_dAIC = NA_real_, tally_wins = NA_integer_, mean_weight = NA_real_))
  s <- s %>% group_by(curveID) %>%
    mutate(dAIC = AIC - min(AIC),
           wt   = exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC)),
           win  = AIC == min(AIC)) %>% ungroup()
  out <- s %>% group_by(model) %>%
    summarise(sumAIC = sum(AIC), mean_weight = mean(wt),
              tally_wins = sum(win), .groups = "drop")
  out$n_curves  <- n
  out$mean_dAIC <- (out$sumAIC - min(out$sumAIC)) / n
  as.data.frame(out[, c("model", "n_curves", "sumAIC", "mean_dAIC", "tally_wins", "mean_weight")])
}

# --------------------------- per-variable summary ---------------------------
class_rows <- list(); within_rows <- list()
for (v in VARS) {
  sub <- percurve[percurve$variable == v, ]
  ss  <- sub[sub$model == "SS", c("curveID", "converged", "interior")]
  ss$peaked <- ss$converged & !is.na(ss$interior) & ss$interior

  n_total   <- length(unique(sub$curveID))
  n_ss_conv <- sum(ss$converged)
  n_peaked  <- sum(ss$peaked)
  frac      <- if (n_ss_conv > 0) n_peaked / n_ss_conv else NA_real_
  verdict   <- if (!is.na(frac) && frac >= PEAK_FRAC) "peaked" else "monotonic"
  conv <- tapply(sub$converged, sub$model, sum)

  peaked_ids  <- ss$curveID[ss$peaked]
  monoton_ids <- setdiff(unique(sub$curveID), peaked_ids)

  class_rows[[v]] <- data.frame(
    variable = v, n_total = n_total, n_SS_conv = n_ss_conv, n_peaked = n_peaked,
    frac_peaked = round(frac, 3), verdict = verdict,
    conv_linear = conv["linear"], conv_exp = conv["exp"],
    conv_SS = conv["SS"], conv_arroyo = conv["arroyo"],
    stringsAsFactors = FALSE)

  wp <- within_compare(sub, peaked_ids,  c("SS", "arroyo"))   # peaked family
  wm <- within_compare(sub, monoton_ids, c("linear", "exp"))  # monotonic family
  wp$variable <- v; wp$family <- "peaked (SS vs Arroyo)"
  wm$variable <- v; wm$family <- "monotonic (linear vs exp)"
  within_rows[[v]] <- rbind(wp, wm)
}
class_tbl  <- do.call(rbind, class_rows);  rownames(class_tbl)  <- NULL
within_tbl <- do.call(rbind, within_rows); rownames(within_tbl) <- NULL
write.csv(class_tbl,  OUT_CLASS,  row.names = FALSE)
write.csv(within_tbl, OUT_WITHIN, row.names = FALSE)
message("Wrote classification: ", OUT_CLASS)
message("Wrote within-family comparisons: ", OUT_WITHIN)

# ------------------------------- report -------------------------------------
cat("\n===================== SHAPE CLASSIFICATION (interior-optimum gate) =====================\n")
cat(sprintf("(eWUE = manuscript WUE = A/E; margin = %g C; peaked if >= %g of curves interior)\n\n",
            MARGIN, PEAK_FRAC))
print(class_tbl[, c("variable","n_total","n_SS_conv","n_peaked","frac_peaked","verdict")],
      row.names = FALSE)

cat("\nSS convergence and other models (n converged of ", nrow(dat[!duplicated(dat$curveID), ]),
    " curves):\n", sep = "")
print(class_tbl[, c("variable","conv_linear","conv_exp","conv_SS","conv_arroyo")], row.names = FALSE)

cat("\n===================== WITHIN-FAMILY AIC (only where both models converged) =============\n")
wt <- within_tbl
wt$sumAIC     <- round(wt$sumAIC, 1)
wt$mean_dAIC  <- round(wt$mean_dAIC, 2)
wt$mean_weight<- round(wt$mean_weight, 3)
for (v in VARS) {
  cat(sprintf("\n---- %s ----\n", v))
  print(wt[wt$variable == v, c("family","model","n_curves","sumAIC","mean_dAIC","tally_wins","mean_weight")],
        row.names = FALSE)
}
cat("\nExpected: A and iWUE peaked; E, gsw, eWUE(WUE) monotonic.\n")
cat("Peaked-family rows test the Methods claim (does SS beat Arroyo among peaked curves?).\n")
