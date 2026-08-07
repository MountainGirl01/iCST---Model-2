# ------------------------------------------------------------------
# Install Packages
# ------------------------------------------------------------------

install.packages(c(
  "brms", "tidyverse", "bayesplot", "tidybayes", "loo", "usethis",
  "performance", "car", "Hmisc", "corrplot"
))

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
  # ------------------------------------------------------------------
# Multicollinearity check
# ------------------------------------------------------------------
# VIF checks whether predictors are too strongly correlated with each
# other, which would inflate standard errors / posterior uncertainty.
# Rule of thumb: VIF > 5 = concerning, VIF > 10 = serious problem.

# Option A: performance package - works directly on brms objects
check_collinearity(fit_primary)

# Option B: car::vif - requires an equivalent frequentist lm() model
# (brms doesn't have its own vif(); this proxy model is only used to
# check the design matrix, not to draw inferential conclusions from)
lm_check <- lm(ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation * z_BASELINE_CQCPR_20,
               data = df_complete)
vif(lm_check)

# Note: for models with an interaction term, VIF on the interaction
# itself will often look inflated by construction (this is expected
# and not usually a concern) - focus mainly on the VIFs for the
# main effects: c_BASELINE_ADAScog, Randomisation, z_BASELINE_CQCPR_20

# Simple correlation check between the two continuous predictors
cor(df_complete$c_BASELINE_ADAScog, df_complete$z_BASELINE_CQCPR_20,
    use = "complete.obs")

# ------------------------------------------------------------------
# Linearity Check
# ------------------------------------------------------------------
# Checks whether the relationship between each continuous predictor
# and the outcome is reasonably linear (as assumed by the model),
# rather than curved/non-linear.

# --- 2a. Residuals vs fitted values plot ---
# Should show a random scatter around 0, with no obvious curve/pattern.

df_complete$fitted_vals <- fitted(fit_primary)[, "Estimate"]
df_complete$resid_vals  <- residuals(fit_primary)[, "Estimate"]

ggplot(df_complete, aes(x = fitted_vals, y = resid_vals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(x = "Fitted values", y = "Residuals",
       title = "Residuals vs Fitted - check for curvature") +
  theme_minimal()

# --- 2b. Residuals vs each continuous predictor ---
# A curved smoothed line (rather than flat) suggests a non-linear
# relationship between that predictor and the outcome.

ggplot(df_complete, aes(x = c_BASELINE_ADAScog, y = resid_vals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(x = "Baseline ADAS-cog (centered)", y = "Residuals",
       title = "Residuals vs Baseline ADAS-cog") +
  theme_minimal()

ggplot(df_complete, aes(x = z_BASELINE_CQCPR_20, y = resid_vals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(x = "QPRC (z-score)", y = "Residuals",
       title = "Residuals vs QPRC") +
  theme_minimal()

# --- 2c. performance package - all-in-one diagnostic panel ---
# Note: performance::check_model() is built primarily for frequentist
# models; for brms it will use posterior predictive draws where
# possible. Can be slow for large models - optional.

# check_model(fit_primary)

# --- 2d. Observed vs predicted scatter (overall model fit check) ---
ggplot(df_complete, aes(x = fitted_vals, y = ADAS_20_FU2)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Predicted ADAS_20_FU2", y = "Observed ADAS_20_FU2",
       title = "Observed vs Predicted - points should hug the diagonal") +
  theme_minimal()

# ------------------------------------------------------------------
# Check Correlations
# ------------------------------------------------------------------

# Select just the continuous variables relevant to the model
cor_vars <- df_complete[, c("ADAS_20_FU2", "c_BASELINE_ADAScog", "z_BASELINE_CQCPR_20")]

# Basic correlation matrix
cor_matrix <- cor(cor_vars, use = "complete.obs")
cor_matrix

# Round for readability
round(cor_matrix, 3)

# With p-values (tests whether each correlation is significantly different from 0)
# install.packages("Hmisc")   # if not already installed
library(Hmisc)
cor_results <- rcorr(as.matrix(cor_vars))
cor_results$r   # correlation coefficients
cor_results$P   # p-values

# Visual correlation matrix (heatmap-style)
library(corrplot)
corrplot(cor_matrix, method = "number", type = "upper",
         tl.col = "black", tl.srt = 45)
  
## ------------------------------------------------------------------
## Prior specification
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
# Primary Priors
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
# Optimistic Model
# ------------------------------------------------------------------

priors_optimistic <- priors_primary
priors_optimistic[priors_optimistic$coef == "RandomisationiCST" &
                    priors_optimistic$class == "b", "prior"] <- "normal(-2.9, 2)"

# ------------------------------------------------------------------
# Pessimistic Model
# ------------------------------------------------------------------

priors_pessimistic <- priors_primary
priors_pessimistic[priors_pessimistic$coef == "RandomisationiCST" &
                     priors_pessimistic$class == "b", "prior"] <- "normal(-0.5, 2)"

# ------------------------------------------------------------------
# Primary Model
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
# Sensitivity Models
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


# ============================================================
# Run predictive checks
# ============================================================

p_prior <- pp_check(prior_check_original, ndraws = 100) +
  labs(title = "Prior Predictive Check",
       x = "ADAS-Cog score",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "right")

p_posterior <- pp_check(fit_primary, ndraws = 100) +
  labs(title = "Posterior Predictive Check",
       x = "ADAS-Cog score",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "right")

library(patchwork)

combined_plot <- (p_prior | p_posterior) +
  plot_annotation(
    title = "Prior and Posterior Predictive Checks",
    theme = theme(plot.title = element_text(face = "italic", size = 14))
  )

combined_plot

# ============================================================
# Overlaid Posterior Distributions with Probability Labels
# ============================================================

install.packages("ggrepel")
library(ggrepel)

overlay_plot <- ggplot(all_draws, aes(x = b_RandomisationiCST, fill = prior_type,
                                      color = prior_type)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel(
    data = prob_labels,
    aes(x = mean_est, y = label_y, label = label, color = prior_type),
    inherit.aes = FALSE,
    size = 3.5, fontface = "bold",
    show.legend = FALSE,
    box.padding = 0.5,
    max.overlaps = Inf,
    min.segment.length = Inf   # this removes the connector lines entirely
  ) +
  labs(
    title = "Posterior Distributions of Treatment Effect by Prior Specification",
    subtitle = "iCST vs TAU Control on ADAS-Cog at 26 weeks",
    x = "Treatment effect (\u03b2)",
    y = "Density",
    fill = "Prior specification",
    color = "Prior specification"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")

overlay_plot


# ------------------------------------------------------------------
# Interaction term: estimate and posterior probabilities
# ------------------------------------------------------------------

summary(fit_primary)$fixed["RandomisationiCST:z_BASELINE_CQCPR_20", ]

hypothesis(fit_primary, "RandomisationiCST:z_BASELINE_CQCPR_20 > 0")
hypothesis(fit_primary, "RandomisationiCST:z_BASELINE_CQCPR_20 < 0")

# ------------------------------------------------------------------
# Trace Plots
# ------------------------------------------------------------------

plot(fit_primary)

# ------------------------------------------------------------------
# Forest Plot
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
# Interaction/ Moderation Plot
# ------------------------------------------------------------------

conditional_effects(fit_primary, effects = "z_BASELINE_CQCPR_20:Randomisation")

# ------------------------------------------------------------------
# Save fitted models
# ------------------------------------------------------------------

saveRDS(fit_primary, "fit_primary.rds")
saveRDS(fit_optimistic, "fit_optimistic.rds")
saveRDS(fit_pessimistic, "fit_pessimistic.rds")

# ------------------------------------------------------------------
# Summary table of results
# ------------------------------------------------------------------

reporting_gt <- reporting_table %>%
  gt() %>%
  tab_header(
    title = "Model 2: Posterior Summary of Fixed Effects",
    subtitle = "ADAS-Cog at 26 weeks ~ Baseline ADAS-Cog + Randomisation x QPRC"
  ) %>%
  cols_label(
    Parameter = "Parameter",
    Estimate = "Estimate",
    Est.Error = "SE",
    `l-95% CI` = "95% CI (lower)",
    `u-95% CI` = "95% CI (upper)",
    Rhat = "Rhat",
    Bulk_ESS = "Bulk ESS",
    Tail_ESS = "Tail ESS",
    Posterior_Probability = "Posterior Probability"
  ) %>%
  fmt_number(columns = c(Estimate, Est.Error, `l-95% CI`, `u-95% CI`, Rhat,
                         Posterior_Probability),
             decimals = 3) %>%
  sub_missing(columns = Posterior_Probability, missing_text = "\u2014") %>%
  tab_source_note(source_note = paste0("N = ", nrow(df_complete),
                                       " complete cases. Posterior probability reflects direction ",
                                       "consistent with the estimate's sign."))

reporting_gt

