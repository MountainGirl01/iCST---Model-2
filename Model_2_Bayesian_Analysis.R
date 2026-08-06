# ------------------------------------------------------------------
# 1. Load data and shorten name for convenience
# ------------------------------------------------------------------

df <- THESIS.COMBINED.Complete.case.26.week.follow.up.data

# ------------------------------------------------------------------
# 2. Fix miscoded missing values
# ------------------------------------------------------------------
# BASELINE_CQCPR_20 had 0s that were actually missing-value codes,
# not genuine scores (confirmed against jamovi: 3 missing cases,
# jamovi mean = 59.4, SD = 6.46). Recode 0 -> NA before use.

df$BASELINE_CQCPR_20[df$BASELINE_CQCPR_20 == 0] <- NA

# ------------------------------------------------------------------
# 3. Filter to complete cases on the variables needed for the model
# ------------------------------------------------------------------
# Matches model 1's complete-case approach.

df_complete <- df[complete.cases(df[, c("ADAS_20_FU2",
                                        "c_BASELINE_ADAScog",
                                        "Randomisation",
                                        "BASELINE_CQCPR_20")]), ]

nrow(df_complete)   # sample size after filtering (258)

# ------------------------------------------------------------------
# 4. Standardise BASELINE_CQCPR_20 (QPRC) into a z-score
# ------------------------------------------------------------------

> # c_BASELINE_ADAScog is left as centered-only (raw ADAS-cog points),
  > # per earlier decision - it's a covariate, not part of the interaction,
  > # and keeping it in points preserves clinical interpretability.
  > 
  > df_complete <- df_complete %>%
    +     mutate(z_BASELINE_CQCPR_20 = as.numeric(scale(BASELINE_CQCPR_20)))
  > 
    > # Sanity checks
    > mean(df_complete$BASELINE_CQCPR_20, na.rm = TRUE)      # ~59.38, matches jamovi's 59.4
  [1] 59.38339
  > sd(df_complete$BASELINE_CQCPR_20, na.rm = TRUE)         # ~6.46, matches jamovi's 6.46
  [1] 6.46475
  > mean(df_complete$z_BASELINE_CQCPR_20, na.rm = TRUE)     # ~0
  [1] -3.256439e-16
  > sd(df_complete$z_BASELINE_CQCPR_20, na.rm = TRUE)       # 1
  [1] 1
  > 

  ## ------------------------------------------------------------------
## Model 2: Prior specification
## ------------------------------------------------------------------
## Coefficient names below assume Randomisation is coded with a
## reference level (e.g. "TAU") and one treatment level - check the
## exact name via get_prior() first and adjust `coef = "..."` to match.

get_prior(adas_formula, data = df_complete)

# Run this and find the exact name of the Randomisation coefficient,
# e.g. it might show as "RandomisationiCST" or similar - substitute
# that exact string into `coef = ` below wherever you see
# "RandomisationiCST".

# ------------------------------------------------------------------
# PRIMARY PRIORS
# ------------------------------------------------------------------

priors_primary <- c(
  # Beta 0: Intercept
  prior(normal(20, 7), class = "Intercept"),
  
  # Beta 1: Baseline ADAS-cog (retained from Model 1)
  prior(normal(0.5, 0.3), class = "b", coef = "c_BASELINE_ADAScog"),
  
  # Beta 2: Treatment effect (primary estimate, from CST literature)
  prior(normal(-1.92, 2), class = "b", coef = "RandomisationiCST"),
  
  # Beta 3: Main effect of QPRC
  prior(normal(0, 2), class = "b", coef = "z_BASELINE_CQCPR_20"),
  
  # Beta 4: Interaction (Group x QPRC)
  prior(normal(0, 2), class = "b",
        coef = "RandomisationiCST:z_BASELINE_CQCPR_20"),
  
  # Sigma: residual SD (same as Model 1)
  prior(exponential(0.1), class = "sigma")
)

# Check the primary priors are valid against the model/data
validate_prior(priors_primary, adas_formula, data = df_complete)

# ------------------------------------------------------------------
# SENSITIVITY ANALYSIS: OPTIMISTIC TREATMENT PRIOR
# ------------------------------------------------------------------

priors_optimistic <- priors_primary
priors_optimistic[priors_optimistic$coef == "RandomisationiCST" &
                    priors_optimistic$class == "b", "prior"] <- "normal(-2.9, 2)"

# ------------------------------------------------------------------
# SENSITIVITY ANALYSIS: PESSIMISTIC TREATMENT PRIOR
# ------------------------------------------------------------------

priors_pessimistic <- priors_primary
priors_pessimistic[priors_pessimistic$coef == "RandomisationiCST" &
                     priors_pessimistic$class == "b", "prior"] <- "normal(-0.5, 2)"

# ------------------------------------------------------------------
# FIT: PRIMARY MODEL
# ------------------------------------------------------------------

fit_primary <- brm(
  formula = adas_formula,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_primary,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 1234,
  control = list(adapt_delta = 0.95)
)

summary(fit_primary)

# ------------------------------------------------------------------
# FIT: SENSITIVITY MODELS
# ------------------------------------------------------------------

fit_optimistic <- brm(
  formula = adas_formula,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_optimistic,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 1234,
  control = list(adapt_delta = 0.95)
)

fit_pessimistic <- brm(
  formula = adas_formula,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_pessimistic,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 1234,
  control = list(adapt_delta = 0.95)
)

# Compare treatment effect estimates across all three priors
summary(fit_primary)$fixed["RandomisationiCST", ]
summary(fit_optimistic)$fixed["RandomisationiCST", ]
summary(fit_pessimistic)$fixed["RandomisationiCST", ]