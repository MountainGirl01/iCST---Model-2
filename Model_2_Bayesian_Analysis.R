# ------------------------------------------------------------------
# Load data and shorten name for convenience
# ------------------------------------------------------------------

df <- THESIS.COMBINED.Complete.case.26.week.follow.up.data

# ------------------------------------------------------------------
# Fix miscoded missing values
# ------------------------------------------------------------------
# BASELINE_CQCPR_20 had 0s that were actually missing-value codes,
# not genuine scores (confirmed against jamovi: 3 missing cases,
# jamovi mean = 59.4, SD = 6.46). Recode 0 -> NA before use.

df$BASELINE_CQCPR_20[df$BASELINE_CQCPR_20 == 0] <- NA

# ------------------------------------------------------------------
# Filter to complete cases on the variables needed for the model
# ------------------------------------------------------------------
# Matches model 1's complete-case approach.

df_complete <- df[complete.cases(df[, c("ADAS_20_FU2",
                                        "c_BASELINE_ADAScog",
                                        "Randomisation",
                                        "BASELINE_CQCPR_20")]), ]

nrow(df_complete)   # sample size after filtering (258)

# ------------------------------------------------------------------
# Standardise BASELINE_CQCPR_20 (QPRC) into a z-score
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

# ------------------------------------------------------------------
# Interaction term: estimate and posterior probabilities
# ------------------------------------------------------------------

summary(fit_primary)$fixed["RandomisationiCST:z_BASELINE_CQCPR_20", ]

hypothesis(fit_primary, "RandomisationiCST:z_BASELINE_CQCPR_20 > 0")
hypothesis(fit_primary, "RandomisationiCST:z_BASELINE_CQCPR_20 < 0")

# ------------------------------------------------------------------
# DIAGNOSTIC PLOTS
# ------------------------------------------------------------------

# Trace plots for convergence
plot(fit_primary)

# Rhat at a glance
mcmc_rhat(rhat(fit_primary))

# Posterior predictive check - does the model reproduce observed data?
pp_check(fit_primary, ndraws = 100)
pp_check(fit_primary, type = "stat", stat = "mean")

# ------------------------------------------------------------------
# COEFFICIENT / FOREST PLOT
# ------------------------------------------------------------------

fit_primary %>%
  gather_draws(b_c_BASELINE_ADAScog, b_RandomisationiCST,
               b_z_BASELINE_CQCPR_20, `b_RandomisationiCST:z_BASELINE_CQCPR_20`) %>%
  ggplot(aes(y = .variable, x = .value)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Estimate", y = NULL,
       title = "Posterior distributions of Model 2 fixed effects") +
  theme_minimal()

# ------------------------------------------------------------------
# TREATMENT EFFECT POSTERIOR (single parameter highlight)
# ------------------------------------------------------------------

mcmc_areas(fit_primary, pars = "b_RandomisationiCST", prob = 0.95) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Posterior distribution: iCST treatment effect")

# ------------------------------------------------------------------
# INTERACTION / MODERATION PLOT
# ------------------------------------------------------------------
# This is the key plot answering the research question: does QPRC
# moderate the iCST treatment effect on ADAS-cog at 26 weeks?

conditional_effects(fit_primary, effects = "z_BASELINE_CQCPR_20:Randomisation")

# ------------------------------------------------------------------
# SENSITIVITY ANALYSIS COMPARISON PLOT
# ------------------------------------------------------------------
# Overlays the treatment effect posterior across all three priors,
# showing how much the conclusion depends on prior choice.

draws_primary <- as_draws_df(fit_primary) %>%
  select(b_RandomisationiCST) %>%
  mutate(prior_type = "Primary")

draws_optimistic <- as_draws_df(fit_optimistic) %>%
  select(b_RandomisationiCST) %>%
  mutate(prior_type = "Optimistic")

draws_pessimistic <- as_draws_df(fit_pessimistic) %>%
  select(b_RandomisationiCST) %>%
  mutate(prior_type = "Pessimistic")

all_draws <- bind_rows(draws_primary, draws_optimistic, draws_pessimistic)

ggplot(all_draws, aes(x = b_RandomisationiCST, fill = prior_type)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Treatment effect (iCST vs TAU Control)", y = "Density",
       fill = "Prior specification",
       title = "Sensitivity analysis: treatment effect across prior specifications") +
  theme_minimal()

# ------------------------------------------------------------------
# Save fitted models
# ------------------------------------------------------------------

saveRDS(fit_primary, "fit_primary.rds")
saveRDS(fit_optimistic, "fit_optimistic.rds")
saveRDS(fit_pessimistic, "fit_pessimistic.rds")

