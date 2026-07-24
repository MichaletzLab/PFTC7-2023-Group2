# Shared helpers for Figs. 4, 5, S4, and S5.
# Source after PCA.R has run, so that raw.env.data_pca exists.
#
#   ALPHA, p_linetype(), p_fmt()  significance styling and panel annotation
#   species_labels()              plotmath labels for the species legend
#   species_level_set()           canonical species factor levels
#   fe_trend()                    fixed-effect trend and 95% CI

# ---- significance conventions ----------------------------------------------
ALPHA <- 0.05

p_linetype <- function(p) if (!is.na(p) && p < ALPHA) "solid" else "dashed"

p_fmt <- function(p) {
  if (is.na(p)) return("italic(p) == NA")
  if (p < 0.001) return('italic(p) < "0.001"')
  sprintf('italic(p) == "%s"',
          formatC(p, format = "f", digits = if (p < 0.01) 3 else 2))
}

# ---- species legend --------------------------------------------------------
# Build plotmath labels: italic binomial, roman "var."
species_labels <- function(levels) {
  out <- lapply(levels, function(s) {
    if (grepl(" var. ", s, fixed = TRUE)) {
      parts <- strsplit(s, " var. ", fixed = TRUE)[[1]]
      g <- strsplit(parts[1], " ")[[1]]     # genus + species
      bquote(italic(.(g[1]))~italic(.(g[2]))~"var."~italic(.(parts[2])))
    } else {
      g <- strsplit(s, " ")[[1]]
      if (length(g) == 2) bquote(italic(.(g[1]))~italic(.(g[2])))
      else                bquote(italic(.(s)))   # safety for odd names
    }
  })
  do.call(expression, out)
}

# Canonical species levels, set in DataPrep.R and carried through PCA.R.
# Used with drop = FALSE so colors are identical across Figures 4, 5, S4, S5
# regardless of which curves each panel retains.
species_level_set <- function() {
  if (!exists("raw.env.data_pca", envir = globalenv()))
    stop("raw.env.data_pca not found. Run DataPrep.R and PCA.R first.")
  sp <- get("raw.env.data_pca", envir = globalenv())$Species
  if (!is.factor(sp)) stop("raw.env.data_pca$Species is not a factor.")
  levels(sp)
}

# ---- fixed-effect trend ----------------------------------------------------
# Trend and 95% CI for a + b * x, where a and b are two named fixed effects.
#   Figure 4  a = "(Intercept)", b = PC_AXIS        parameter vs PC
#   Figure 5  a = "Tc",          b = "Tc:<PC_AXIS>" temperature sensitivity vs PC
#             a = "(Intercept)", b = PC_AXIS        value at Tref vs PC
# The returned x column is named "x" and carries no axis identity, so the same
# function serves any predictor. Coefficients are selected by name, not position.
fe_trend <- function(mod, a_name, b_name, xseq) {
  fe <- lme4::fixef(mod)
  V  <- as.matrix(vcov(mod))
  missing_terms <- setdiff(c(a_name, b_name), names(fe))
  if (length(missing_terms))
    stop("fixed effect(s) not found in model: ",
         paste(missing_terms, collapse = ", "))
  a <- fe[[a_name]]; b <- fe[[b_name]]
  Vaa <- V[a_name, a_name]; Vbb <- V[b_name, b_name]; Vab <- V[a_name, b_name]
  fit <- a + b * xseq
  se  <- sqrt(Vaa + xseq^2 * Vbb + 2 * xseq * Vab)
  data.frame(x = xseq, fit = fit, lwr = fit - 1.96 * se, upr = fit + 1.96 * se)
}