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

## ------------------------------------------------------------------
## Model 2 Categorical: Interaction using 3-Category QCPR (Round-Number Split)
## Below 55 (n=74) | 55 to 62 (n=97) | Above 62 (n=87)
## ------------------------------------------------------------------

library(brms)
library(dplyr)

# ------------------------------------------------------------------
# 1. Confirm the categorical variable and set reference level
# ------------------------------------------------------------------
# Reference level = "Below 55", so both other coefficients are
# interpreted relative to the lowest QCPR group.

df_complete$QCPR_3cat_round <- factor(
  df_complete$QCPR_3cat_round,
  levels = c("Below 55", "55 to 62", "Above 62")
)

table(df_complete$QCPR_3cat_round)

# ------------------------------------------------------------------
# 2. Model formula
# ------------------------------------------------------------------

adas_formula_3cat <- bf(
  ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation * QCPR_3cat_round
)

get_prior(adas_formula_3cat, data = df_complete)
# Check exact coefficient names before finalising priors - expect:
# RandomisationiCST, QCPR_3cat_round55to62, QCPR_3cat_roundAbove62,
# RandomisationiCST:QCPR_3cat_round55to62,
# RandomisationiCST:QCPR_3cat_roundAbove62

# ------------------------------------------------------------------
# 3. Priors
# ------------------------------------------------------------------
# Same core structure as the primary model. QCPR category and
# interaction priors are weakly informative (wider SD) given there
# is no strong prior literature for this specific categorical split.

priors_3cat <- c(
  prior(normal(20, 7),   class = "Intercept"),
  prior(normal(0.5, 0.3), class = "b", coef = "c_BASELINE_ADAScog"),
  prior(normal(-1.92, 2), class = "b", coef = "RandomisationiCST"),
  prior(normal(0, 5),     class = "b", coef = "QCPR_3cat_round55to62"),
  prior(normal(0, 5),     class = "b", coef = "QCPR_3cat_roundAbove62"),
  prior(normal(0, 5),     class = "b",
        coef = "RandomisationiCST:QCPR_3cat_round55to62"),
  prior(normal(0, 5),     class = "b",
        coef = "RandomisationiCST:QCPR_3cat_roundAbove62"),
  prior(exponential(0.1), class = "sigma")
)

# NOTE: confirm these coefficient names exactly match your
# get_prior() output above before running - brms naming for
# factor levels can vary slightly (e.g. spaces/punctuation removed).

# ------------------------------------------------------------------
# 4. Fit the model
# ------------------------------------------------------------------

fit_3cat <- brm(
  formula = adas_formula_3cat,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_3cat,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 1234,
  control = list(adapt_delta = 0.95)
)

summary(fit_3cat)

# ------------------------------------------------------------------
# 5. Interaction terms: estimates, CIs, posterior probabilities
# ------------------------------------------------------------------

summary(fit_3cat)$fixed["RandomisationiCST:QCPR_3cat_round55to62", ]
summary(fit_3cat)$fixed["RandomisationiCST:QCPR_3cat_roundAbove62", ]

hypothesis(fit_3cat, "RandomisationiCST:QCPR_3cat_round55to62 > 0")
hypothesis(fit_3cat, "RandomisationiCST:QCPR_3cat_roundAbove62 > 0")

# ------------------------------------------------------------------
# 6. Diagnostics
# ------------------------------------------------------------------

plot(fit_3cat)
brms::pp_check(fit_3cat, ndraws = 100)
cat("Max Rhat:", round(max(brms::rhat(fit_3cat)), 4), "\n")

divergences_3cat <- sum(subset(brms::nuts_params(fit_3cat),
                               Parameter == "divergent__")$Value)
cat("Divergent transitions:", divergences_3cat, "\n")

# ------------------------------------------------------------------
# 7. Interaction plot (conditional effects)
# ------------------------------------------------------------------

conditional_effects(fit_3cat, effects = "QCPR_3cat_round:Randomisation")

# ------------------------------------------------------------------
# 8. Slope-style interaction plot
# ------------------------------------------------------------------

library(ggplot2)

new_data_3cat <- expand.grid(
  Randomisation = levels(df_complete$Randomisation),
  QCPR_3cat_round = levels(df_complete$QCPR_3cat_round),
  c_BASELINE_ADAScog = 0
)

preds_3cat <- fitted(fit_3cat, newdata = new_data_3cat, summary = TRUE,
                     re_formula = NA)

plot_data_3cat <- cbind(new_data_3cat, preds_3cat)

interaction_slope_3cat <- ggplot(plot_data_3cat,
                                 aes(x = QCPR_3cat_round, y = Estimate, color = Randomisation,
                                     group = Randomisation)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0.1) +
  labs(
    title = "Interaction: Treatment x QCPR Category (3-group split)",
    subtitle = "Predicted ADAS-Cog at 26 weeks (with 95% credible intervals)",
    x = "Baseline QCPR Category",
    y = "Predicted ADAS-Cog at 26 weeks",
    color = "Randomisation"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

interaction_slope_3cat

ggsave("interaction_slope_plot_3cat.png", interaction_slope_3cat,
       width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------
# UPDATED MODEL 2 - 2 PRIORS
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# UPDATE ZSCORE
# ------------------------------------------------------------------
library(dplyr)

df_complete <- df_complete %>%
  mutate(z_BASELINE_CQCPR_20 = as.numeric(scale(BASELINE_CQCPR_20)))

mean(df_complete$z_BASELINE_CQCPR_20, na.rm = TRUE)
sd(df_complete$z_BASELINE_CQCPR_20, na.rm = TRUE)

df_complete[df_complete$BASELINE_CQCPR_20 == 69, c("BASELINE_CQCPR_20", "z_BASELINE_CQCPR_20")]

## ------------------------------------------------------------------
## Model 2: Prior Specification
## ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation * z_BASELINE_CQCPR_20
## Matches Model 1's finalized structure:
##   - Outcome truncated to the plausible ADAS-Cog scale range
##   - Sigma: half-normal, calibrated to the observed residual SD
##   - Two-prior structure: Informative vs Skeptical (treatment effect)
## ------------------------------------------------------------------

library(brms)

# ------------------------------------------------------------------
# 1. Confirm the observed SD of the outcome, to calibrate sigma
# ------------------------------------------------------------------

sd(df_complete$ADAS_20_FU2, na.rm = TRUE)
range(df_complete$ADAS_20_FU2, na.rm = TRUE)

# ------------------------------------------------------------------
# 2. Confirm Randomisation factor levels (TAU Control = reference)
# ------------------------------------------------------------------

df_complete$Randomisation <- factor(df_complete$Randomisation,
                                    levels = c("TAU Control", "iCST"))
levels(df_complete$Randomisation)

# ------------------------------------------------------------------
# 3. Model formula with truncated outcome
# ------------------------------------------------------------------
# Bounds match Model 1: 0 (true floor) to 60 (safely above observed max)

adas_formula_m2_trunc <- bf(
  ADAS_20_FU2 | trunc(lb = 0, ub = 60) ~
    c_BASELINE_ADAScog + Randomisation * z_BASELINE_CQCPR_20
)

get_prior(adas_formula_m2_trunc, data = df_complete)
# Confirm exact coefficient names before finalising priors - expect:
# c_BASELINE_ADAScog, RandomisationiCST, z_BASELINE_CQCPR_20,
# RandomisationiCST:z_BASELINE_CQCPR_20

# ------------------------------------------------------------------
# 4. INFORMATIVE PRIORS
# ------------------------------------------------------------------
# Intercept, baseline, and sigma match Model 1's finalized values
# (same outcome variable, same population). Treatment effect prior
# is literature-derived, matching Model 1's informative specification.
# QPRC main effect and interaction remain weakly informative
# (centred at 0) in both prior specifications, since there is no
# strong prior literature estimating these effects.

priors_informative_m2 <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-1.92, 2), class = b, coef = RandomisationiCST),
  prior(normal(0, 2),     class = b, coef = z_BASELINE_CQCPR_20),
  prior(normal(0, 2),     class = b,
        coef = "RandomisationiCST:z_BASELINE_CQCPR_20"),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

# ------------------------------------------------------------------
# 5. SKEPTICAL PRIORS
# ------------------------------------------------------------------
# Identical structure, except the treatment effect is centred at 0
# (no assumed effect), matching Model 1's skeptical specification.

priors_skeptical_m2 <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(0, 2),     class = b, coef = RandomisationiCST),
  prior(normal(0, 2),     class = b, coef = z_BASELINE_CQCPR_20),
  prior(normal(0, 2),     class = b,
        coef = "RandomisationiCST:z_BASELINE_CQCPR_20"),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

# ------------------------------------------------------------------
# 6. Validate both prior sets against the model/data
# ------------------------------------------------------------------

validate_prior(priors_informative_m2, adas_formula_m2_trunc, data = df_complete)
validate_prior(priors_skeptical_m2, adas_formula_m2_trunc, data = df_complete)

## ------------------------------------------------------------------
## Model 2: Fit Informative and Skeptical Models
## Matches Model 1's settings: chains=4, iter=4000, warmup=2000,
## adapt_delta=0.95
## ------------------------------------------------------------------

library(brms)

# ------------------------------------------------------------------
# 1. FIT: INFORMATIVE MODEL
# ------------------------------------------------------------------

fit_informative_m2 <- brm(
  formula = adas_formula_m2_trunc,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_informative_m2,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model2_informative_final"
)

# ------------------------------------------------------------------
# 2. FIT: SKEPTICAL MODEL
# ------------------------------------------------------------------

fit_skeptical_m2 <- brm(
  formula = adas_formula_m2_trunc,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_skeptical_m2,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model2_skeptical_final"
)

# ------------------------------------------------------------------
# 3. SUMMARY - both models
# ------------------------------------------------------------------

summary(fit_informative_m2)
summary(fit_skeptical_m2)

# ------------------------------------------------------------------
# 4. QUICK CONVERGENCE CHECK
# ------------------------------------------------------------------

cat("Max Rhat (Informative):", round(max(brms::rhat(fit_informative_m2)), 4), "\n")
cat("Max Rhat (Skeptical):", round(max(brms::rhat(fit_skeptical_m2)), 4), "\n")

div_inf_m2 <- sum(subset(brms::nuts_params(fit_informative_m2),
                         Parameter == "divergent__")$Value)
div_skep_m2 <- sum(subset(brms::nuts_params(fit_skeptical_m2),
                          Parameter == "divergent__")$Value)

cat("Divergent transitions (Informative):", div_inf_m2, "\n")
cat("Divergent transitions (Skeptical):", div_skep_m2, "\n")

library(brms)

# ------------------------------------------------------------------
# TREATMENT EFFECT: posterior probabilities, both priors
# ------------------------------------------------------------------

hypothesis(fit_informative_m2, "RandomisationiCST < 0")
hypothesis(fit_skeptical_m2, "RandomisationiCST < 0")

# ------------------------------------------------------------------
# INTERACTION TERM: posterior probabilities, both priors
# ------------------------------------------------------------------

hypothesis(fit_informative_m2, "RandomisationiCST:z_BASELINE_CQCPR_20 > 0")
hypothesis(fit_skeptical_m2, "RandomisationiCST:z_BASELINE_CQCPR_20 > 0")

# ------------------------------------------------------------------
# QUICK SUMMARY TABLE - both effects, both priors
# ------------------------------------------------------------------

draws_treat_inf <- as_draws_df(fit_informative_m2)$b_RandomisationiCST
draws_treat_skep <- as_draws_df(fit_skeptical_m2)$b_RandomisationiCST
draws_int_inf <- as_draws_df(fit_informative_m2)$`b_RandomisationiCST:z_BASELINE_CQCPR_20`
draws_int_skep <- as_draws_df(fit_skeptical_m2)$`b_RandomisationiCST:z_BASELINE_CQCPR_20`

prob_summary_m2 <- data.frame(
  Model = c("Informative", "Skeptical"),
  P_treatment_benefit = c(mean(draws_treat_inf < 0), mean(draws_treat_skep < 0)),
  P_interaction_positive = c(mean(draws_int_inf > 0), mean(draws_int_skep > 0))
)

prob_summary_m2$P_treatment_benefit <- round(prob_summary_m2$P_treatment_benefit, 3)
prob_summary_m2$P_interaction_positive <- round(prob_summary_m2$P_interaction_positive, 3)

print(prob_summary_m2)

library(brms)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------
# 1. Fit prior-only models (sample_prior = "only")
# ------------------------------------------------------------------

prior_check_informative_m2 <- brm(
  formula = adas_formula_m2_trunc,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_informative_m2,
  sample_prior = "only",
  chains  = 4,
  iter    = 1000,
  seed    = 42
)

prior_check_skeptical_m2 <- brm(
  formula = adas_formula_m2_trunc,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_skeptical_m2,
  sample_prior = "only",
  chains  = 4,
  iter    = 1000,
  seed    = 42
)

# ------------------------------------------------------------------
# 2. INFORMATIVE - side by side
# ------------------------------------------------------------------

p_prior_inf_m2 <- brms::pp_check(prior_check_informative_m2, ndraws = 100) +
  labs(title = "Prior Predictive Check",
       x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_inf_m2 <- brms::pp_check(fit_informative_m2, ndraws = 100) +
  labs(title = "Posterior Predictive Check",
       x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_informative_m2 <- p_prior_inf_m2 | p_post_inf_m2

combined_informative_m2

ggsave("model2_ppc_informative.png", combined_informative_m2,
       width = 10, height = 5, dpi = 300)

# ------------------------------------------------------------------
# 3. SKEPTICAL - side by side
# ------------------------------------------------------------------

p_prior_skep_m2 <- brms::pp_check(prior_check_skeptical_m2, ndraws = 100) +
  labs(title = "Prior Predictive Check",
       x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_skep_m2 <- brms::pp_check(fit_skeptical_m2, ndraws = 100) +
  labs(title = "Posterior Predictive Check",
       x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_skeptical_m2 <- p_prior_skep_m2 | p_post_skep_m2

combined_skeptical_m2

ggsave("model2_ppc_skeptical.png", combined_skeptical_m2,
       width = 10, height = 5, dpi = 300)


> # ------------------------------------------------------------------
> # Interaction Plot
> # ------------------------------------------------------------------

conditional_effects(fit_informative_m2, effects = "z_BASELINE_CQCPR_20:Randomisation")


> # ------------------------------------------------------------------
> # Table of Results
> # ------------------------------------------------------------------

## ------------------------------------------------------------------
## Model 2: Comprehensive Results Table
## All parameters (Intercept, Baseline, Treatment, QPRC, Interaction,
## Sigma), full diagnostics (Est.Error, CrI, Rhat, ESS), and posterior
## probabilities where relevant - both Informative and Skeptical priors.
## ------------------------------------------------------------------

library(dplyr)
library(brms)
library(gt)

# ------------------------------------------------------------------
# 1. Function to pull all parameters + diagnostics from a fitted model
# ------------------------------------------------------------------

extract_all_params_m2 <- function(fit, model_label) {
  fixed_df <- as.data.frame(summary(fit)$fixed)
  fixed_df$Parameter <- rownames(fixed_df)
  
  sigma_df <- as.data.frame(summary(fit)$spec_pars)
  sigma_df$Parameter <- rownames(sigma_df)
  
  combined <- bind_rows(fixed_df, sigma_df)
  combined$Model <- model_label
  
  combined %>%
    select(Model, Parameter, Estimate, Est.Error, `l-95% CI`, `u-95% CI`,
           Rhat, Bulk_ESS, Tail_ESS)
}

# ------------------------------------------------------------------
# 2. Extract from both models
# ------------------------------------------------------------------

params_inf_m2 <- extract_all_params_m2(fit_informative_m2, "Informative")
params_skep_m2 <- extract_all_params_m2(fit_skeptical_m2, "Skeptical")

full_table_m2 <- bind_rows(params_inf_m2, params_skep_m2)

# ------------------------------------------------------------------
# 3. Clean parameter labels and set display order
# ------------------------------------------------------------------

full_table_m2$Parameter <- dplyr::recode(full_table_m2$Parameter,
                                         "Intercept" = "Intercept",
                                         "c_BASELINE_ADAScog" = "Baseline",
                                         "RandomisationiCST" = "Treatment",
                                         "z_BASELINE_CQCPR_20" = "QPRC (z-score)",
                                         "RandomisationiCST:z_BASELINE_CQCPR_20" = "Treatment x QPRC",
                                         "sigma" = "Sigma"
)

param_order <- c("Intercept", "Baseline", "QPRC (z-score)", "Sigma",
                 "Treatment", "Treatment x QPRC")
full_table_m2$Parameter <- factor(full_table_m2$Parameter, levels = param_order)
full_table_m2$Model <- factor(full_table_m2$Model, levels = c("Informative", "Skeptical"))

# ------------------------------------------------------------------
# 4. Add posterior probability column (treatment & interaction only)
# ------------------------------------------------------------------

draws_treat_inf <- as_draws_df(fit_informative_m2)$b_RandomisationiCST
draws_treat_skep <- as_draws_df(fit_skeptical_m2)$b_RandomisationiCST
draws_int_inf <- as_draws_df(fit_informative_m2)$`b_RandomisationiCST:z_BASELINE_CQCPR_20`
draws_int_skep <- as_draws_df(fit_skeptical_m2)$`b_RandomisationiCST:z_BASELINE_CQCPR_20`

prob_lookup <- tibble::tibble(
  Model = c("Informative", "Informative", "Skeptical", "Skeptical"),
  Parameter = c("Treatment", "Treatment x QPRC", "Treatment", "Treatment x QPRC"),
  P_direction = c(
    mean(draws_treat_inf < 0),    # P(benefit)
    mean(draws_int_inf > 0),      # P(positive)
    mean(draws_treat_skep < 0),
    mean(draws_int_skep > 0)
  )
)

full_table_m2 <- full_table_m2 %>%
  left_join(prob_lookup, by = c("Model", "Parameter")) %>%
  arrange(Model, Parameter) %>%
  mutate(across(c(Estimate, Est.Error, `l-95% CI`, `u-95% CI`, Rhat, P_direction), ~round(.x, 3)),
         across(c(Bulk_ESS, Tail_ESS), ~round(.x, 0)))

print(full_table_m2)

# ------------------------------------------------------------------
# 5. Save raw CSV backup
# ------------------------------------------------------------------

write.csv(full_table_m2, "model2_comprehensive_results_table.csv", row.names = FALSE)

# ------------------------------------------------------------------
# 6. Render as a formatted gt table
# ------------------------------------------------------------------

full_gt_m2 <- full_table_m2 %>%
  rename(`Est. Error` = Est.Error,
         `Lower CrI` = `l-95% CI`,
         `Upper CrI` = `u-95% CI`,
         `Bulk ESS` = Bulk_ESS,
         `Tail ESS` = Tail_ESS,
         `R-hat` = Rhat,
         `Posterior Probability` = P_direction) %>%
  group_by(Model) %>%
  gt() %>%
  tab_header(
    title = "Model 2: Comprehensive Posterior Results",
    subtitle = "ADAS-Cog at 26 weeks ~ Baseline ADAS-Cog + Randomisation x QPRC"
  ) %>%
  fmt_number(columns = c(Estimate, `Est. Error`, `Lower CrI`, `Upper CrI`,
                         `R-hat`, `Posterior Probability`), decimals = 3) %>%
  fmt_number(columns = c(`Bulk ESS`, `Tail ESS`), decimals = 0) %>%
  sub_missing(columns = `Posterior Probability`, missing_text = "\u2014") %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = Parameter) %>%
  tab_source_note(
    source_note = paste0("N = 258 complete cases. Posterior Probability = P(Treatment < 0) ",
                         "for Treatment, P(Treatment x QPRC > 0) for the interaction term. ",
                         "Not applicable (\u2014) for Intercept, Baseline, QPRC main effect, and Sigma.")
  )

full_gt_m2

gtsave(full_gt_m2, "model2_comprehensive_results_table.png")

# ------------------------------------------------------------------
# Formatted Interaction Slope
# ------------------------------------------------------------------

library(ggplot2)

interaction_plot_m2 <- plot(conditional_effects(fit_informative_m2,
                                                effects = "z_BASELINE_CQCPR_20:Randomisation"),
                            points = FALSE)[[1]] +
  labs(
    title = "Interaction: Treatment x Relationship Quality (QPRC)",
    subtitle = "Predicted ADAS-Cog at 26 weeks by treatment group across QPRC (z-score)",
    x = "QPRC (z-score)",
    y = "Predicted ADAS-Cog at 26 weeks"
  ) +
  theme_minimal()

interaction_plot_m2

ggsave("model2_interaction_plot_final.png", interaction_plot_m2,
       width = 9, height = 6, dpi = 300)