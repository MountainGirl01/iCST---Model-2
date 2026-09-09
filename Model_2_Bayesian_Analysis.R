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

# ------------------------------------------------------------------
# Categorical Model 01.09.26
# ------------------------------------------------------------------

## ============================================================================
## Model 2 (Categorical Sensitivity Analysis): QCPR Median Split
## ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation * QCPR_median_group
## Complementary sensitivity analysis to the primary continuous QCPR model.
## ============================================================================

# ---- LIBRARIES ---------------------------------------------------------
library(brms)
library(dplyr)
library(ggplot2)
library(gt)
library(tidybayes)

# ============================================================================
# PART 1: DATA PREPARATION
# ============================================================================

# ---- Step 1: Calculate the median of BASELINE_CQCPR_20 ---------------------

qcpr_median <- median(df_complete$BASELINE_CQCPR_20, na.rm = TRUE)
qcpr_median

# ---- Step 2: Create the categorical (median-split) variable ----------------
# "Below median" is set as the reference level.

df_complete <- df_complete %>%
  mutate(
    QCPR_median_group = ifelse(BASELINE_CQCPR_20 <= qcpr_median,
                               "Below median", "Above median"),
    QCPR_median_group = factor(QCPR_median_group,
                               levels = c("Below median", "Above median"))
  )

table(df_complete$QCPR_median_group)

# ---- Step 3: Confirm Randomisation reference level (TAU Control) -----------

df_complete$Randomisation <- factor(df_complete$Randomisation,
                                    levels = c("TAU Control", "iCST"))
levels(df_complete$Randomisation)

# ============================================================================
# PART 2: MODEL SPECIFICATION AND FITTING
# ============================================================================

# ---- Formula (truncated outcome, matching primary model) -------------------

adas_formula_m2_median <- bf(
  ADAS_20_FU2 | trunc(lb = 0, ub = 60) ~
    c_BASELINE_ADAScog + Randomisation * QCPR_median_group
)

get_prior(adas_formula_m2_median, data = df_complete)

# ---- Priors (matching primary model structure) ------------------------------

priors_informative_median <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-1.92, 2), class = b, coef = RandomisationiCST),
  prior(normal(0, 5),     class = b, coef = "QCPR_median_groupAbovemedian"),
  prior(normal(0, 5),     class = b,
        coef = "RandomisationiCST:QCPR_median_groupAbovemedian"),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

priors_skeptical_median <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(0, 2),     class = b, coef = RandomisationiCST),
  prior(normal(0, 5),     class = b, coef = "QCPR_median_groupAbovemedian"),
  prior(normal(0, 5),     class = b,
        coef = "RandomisationiCST:QCPR_median_groupAbovemedian"),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

validate_prior(priors_informative_median, adas_formula_m2_median, data = df_complete)
validate_prior(priors_skeptical_median, adas_formula_m2_median, data = df_complete)

# ---- Fit both models ---------------------------------------------------------

fit_informative_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_informative_median,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model2_median_informative"
)

fit_skeptical_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_skeptical_median,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model2_median_skeptical"
)

summary(fit_informative_median)
summary(fit_skeptical_median)

# ---- Convergence check --------------------------------------------------------

cat("Max Rhat (Informative):", round(max(brms::rhat(fit_informative_median)), 4), "\n")
cat("Max Rhat (Skeptical):", round(max(brms::rhat(fit_skeptical_median)), 4), "\n")

# ============================================================================
# PART 3: RESULTS TABLE
# ============================================================================

extract_all_params_median <- function(fit, model_label) {
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

params_inf_median <- extract_all_params_median(fit_informative_median, "Informative")
params_skep_median <- extract_all_params_median(fit_skeptical_median, "Skeptical")

full_table_median <- bind_rows(params_inf_median, params_skep_median)

full_table_median$Parameter <- dplyr::recode(full_table_median$Parameter,
                                             "Intercept" = "Intercept",
                                             "c_BASELINE_ADAScog" = "Baseline",
                                             "RandomisationiCST" = "Treatment",
                                             "QCPR_median_groupAbovemedian" = "QCPR Group (Above median)",
                                             "RandomisationiCST:QCPR_median_groupAbovemedian" = "Treatment x QCPR Group",
                                             "sigma" = "Sigma"
)

param_order <- c("Intercept", "Baseline", "QCPR Group (Above median)", "Sigma",
                 "Treatment", "Treatment x QCPR Group")
full_table_median$Parameter <- factor(full_table_median$Parameter, levels = param_order)
full_table_median$Model <- factor(full_table_median$Model, levels = c("Informative", "Skeptical"))

draws_treat_inf_med <- as_draws_df(fit_informative_median)$b_RandomisationiCST
draws_treat_skep_med <- as_draws_df(fit_skeptical_median)$b_RandomisationiCST
draws_int_inf_med <- as_draws_df(fit_informative_median)$`b_RandomisationiCST:QCPR_median_groupAbovemedian`
draws_int_skep_med <- as_draws_df(fit_skeptical_median)$`b_RandomisationiCST:QCPR_median_groupAbovemedian`

prob_lookup_median <- tibble::tibble(
  Model = c("Informative", "Informative", "Skeptical", "Skeptical"),
  Parameter = c("Treatment", "Treatment x QCPR Group", "Treatment", "Treatment x QCPR Group"),
  P_direction = c(
    mean(draws_treat_inf_med < 0),
    mean(draws_int_inf_med > 0),
    mean(draws_treat_skep_med < 0),
    mean(draws_int_skep_med > 0)
  )
)

full_table_median <- full_table_median %>%
  left_join(prob_lookup_median, by = c("Model", "Parameter")) %>%
  arrange(Model, match(Parameter, param_order)) %>%
  mutate(across(c(Estimate, Est.Error, `l-95% CI`, `u-95% CI`, Rhat, P_direction), ~round(.x, 3)),
         across(c(Bulk_ESS, Tail_ESS), ~round(.x, 0)))

print(full_table_median)

write.csv(full_table_median, "model2_median_split_results_table.csv", row.names = FALSE)

full_gt_median <- full_table_median %>%
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
    title = "Model 2 (Median-Split Sensitivity Analysis): Posterior Results",
    subtitle = "ADAS-Cog at 26 weeks ~ Baseline ADAS-Cog + Randomisation x QCPR Group"
  ) %>%
  fmt_number(columns = c(Estimate, `Est. Error`, `Lower CrI`, `Upper CrI`,
                         `R-hat`, `Posterior Probability`), decimals = 3) %>%
  fmt_number(columns = c(`Bulk ESS`, `Tail ESS`), decimals = 0) %>%
  sub_missing(columns = `Posterior Probability`, missing_text = "\u2014") %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = Parameter) %>%
  tab_source_note(
    source_note = paste0("N = 258 complete cases. QCPR median = ",
                         round(qcpr_median, 0),
                         ". Reference categories: TAU Control, Below median.")
  )

full_gt_median

gtsave(full_gt_median, "model2_median_split_results_table.png")

# ============================================================================
# PART 4: INTERACTION BAR CHART
# ============================================================================

new_data_median <- expand.grid(
  Randomisation = levels(df_complete$Randomisation),
  QCPR_median_group = levels(df_complete$QCPR_median_group),
  c_BASELINE_ADAScog = 0
)

preds_median <- fitted(fit_informative_median, newdata = new_data_median,
                       summary = TRUE, re_formula = NA)

plot_data_median <- cbind(new_data_median, preds_median)

print(plot_data_median)

interaction_bar_median <- ggplot(plot_data_median,
                                 aes(x = QCPR_median_group, y = Estimate, fill = Randomisation)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                position = position_dodge(width = 0.7), width = 0.15, linewidth = 0.6) +
  scale_fill_manual(values = c("TAU Control" = "#F08080", "iCST" = "#4ECDC4")) +
  labs(
    title = "Interaction: Treatment x Relationship Quality (Median Split)",
    subtitle = "Predicted ADAS-Cog at 26 weeks (with 95% credible intervals)",
    x = "QCPR Group",
    y = "Predicted ADAS-Cog at 26 weeks",
    fill = "Randomisation"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, max(plot_data_median$Q97.5) * 1.1))

interaction_bar_median

ggsave("model2_interaction_barchart_median.png", interaction_bar_median,
       width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 5: DISTRIBUTION HISTOGRAM (WITH MEDIAN CUTOFF MARKED)
# ============================================================================

qcpr_histogram_median <- ggplot(df_complete, aes(x = BASELINE_CQCPR_20)) +
  geom_histogram(binwidth = 2, fill = "steelblue", color = "white", alpha = 0.85) +
  stat_bin(binwidth = 2, geom = "text", aes(label = after_stat(count)),
           vjust = -0.5, size = 3) +
  geom_vline(xintercept = qcpr_median, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "Distribution of Baseline QCPR Scores (Carer-Reported)",
    subtitle = paste0("Dashed line marks the sample median (", qcpr_median, ")"),
    x = "QCPR Score (0-70)",
    y = "Number of participants"
  ) +
  theme_minimal()

qcpr_histogram_median

ggsave("qcpr_histogram_median_split.png", qcpr_histogram_median,
       width = 9, height = 6, dpi = 300)

## ------------------------------------------------------------------
> ## Secondary Table: Change in Relationship Quality (QCPR) by Group
  > ## Standalone table - separate from Table 1 (baseline characteristics).
  > ## Reports descriptive change scores and the between-group Welch's
  > ## t-test comparing amount of change between iCST and TAU Control.
  > ## ------------------------------------------------------------------
> 
  > library(dplyr)
> library(gt)
> 
  > # ------------------------------------------------------------------
> # 1. Calculate change score
  > # ------------------------------------------------------------------
> 
  > dat$qcpr_change <- dat$CQCPR_20_FU2 - dat$BASELINE_CQCPR_20
> 
  > tau_only <- dat[dat$Randomisation == "TAU Control", ]
> icst_only <- dat[dat$Randomisation == "iCST", ]
> 
  > # ------------------------------------------------------------------
> # 2. Run the Welch's t-test
  > # ------------------------------------------------------------------
> 
  > qcpr_ttest <- t.test(qcpr_change ~ Randomisation, data = dat)
> qcpr_ttest

## ============================================================================
## Model 2 (Categorical, Median-Split): QCPR Group Moderation Analysis
## ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation * QCPR_median_group
## Four participants with extreme ADAS-Cog change scores excluded prior
## to model fitting (identified via inspection of the change score
## distribution, which showed right-skew).
## ============================================================================

# ---- LIBRARIES ---------------------------------------------------------
library(brms)
library(dplyr)
library(ggplot2)
library(gt)
library(tidybayes)

# ============================================================================
# PART 1: DATA PREPARATION
# ============================================================================

# ---- Load and clean the data (if not already in session) -------------------

df <- read.csv("~/Downloads/THESIS COMBINED Complete case 26 week follow up data.csv",
               header = TRUE)

df$BASELINE_CQCPR_20[df$BASELINE_CQCPR_20 == 0] <- NA

df_complete <- df[complete.cases(df[, c("ADAS_20_FU2",
                                        "c_BASELINE_ADAScog",
                                        "Randomisation",
                                        "BASELINE_CQCPR_20")]), ]

nrow(df_complete)   # 258 before outlier exclusion

# ---- Exclude the four extreme change-score cases ---------------------------
# qPIN 14201, 17018, 13044, 11051

df_complete <- df_complete[!(df_complete$qPIN %in% c(14201, 17018, 13044, 11051)), ]

nrow(df_complete)   # should be 254

# ---- Set factor levels -----------------------------------------------------

df_complete$Randomisation <- factor(df_complete$Randomisation,
                                    levels = c("TAU Control", "iCST"))
table(df_complete$Randomisation)

# ---- Calculate the median and create the categorical QCPR group -----------

qcpr_median <- median(df_complete$BASELINE_CQCPR_20, na.rm = TRUE)
qcpr_median

df_complete <- df_complete %>%
  mutate(
    QCPR_median_group = ifelse(BASELINE_CQCPR_20 <= qcpr_median,
                               "Below median", "Above median"),
    QCPR_median_group = factor(QCPR_median_group,
                               levels = c("Below median", "Above median"))
  )

table(df_complete$QCPR_median_group)

# ============================================================================
# PART 2: MODEL SPECIFICATION
# ============================================================================

adas_formula_m2_median <- bf(
  ADAS_20_FU2 | trunc(lb = 0, ub = 60) ~
    c_BASELINE_ADAScog + Randomisation * QCPR_median_group
)

get_prior(adas_formula_m2_median, data = df_complete)

priors_informative_median <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-1.92, 2), class = b, coef = RandomisationiCST),
  prior(normal(0, 5),     class = b, coef = "QCPR_median_groupAbovemedian"),
  prior(normal(0, 5),     class = b,
        coef = "RandomisationiCST:QCPR_median_groupAbovemedian"),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

priors_skeptical_median <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(0, 2),     class = b, coef = RandomisationiCST),
  prior(normal(0, 5),     class = b, coef = "QCPR_median_groupAbovemedian"),
  prior(normal(0, 5),     class = b,
        coef = "RandomisationiCST:QCPR_median_groupAbovemedian"),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

validate_prior(priors_informative_median, adas_formula_m2_median, data = df_complete)
validate_prior(priors_skeptical_median, adas_formula_m2_median, data = df_complete)

# ============================================================================
# PART 3: FIT BOTH MODELS
# ============================================================================

fit_informative_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_informative_median,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model2_median_informative_outlier_excluded"
)

fit_skeptical_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_skeptical_median,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model2_median_skeptical_outlier_excluded"
)

summary(fit_informative_median)
summary(fit_skeptical_median)

# ============================================================================
# PART 4: CONVERGENCE CHECK
# ============================================================================

cat("Max Rhat (Informative):", round(max(brms::rhat(fit_informative_median)), 4), "\n")
cat("Max Rhat (Skeptical):", round(max(brms::rhat(fit_skeptical_median)), 4), "\n")

# ============================================================================
# PART 5: POSTERIOR PROBABILITIES
# ============================================================================

hypothesis(fit_informative_median, "RandomisationiCST < 0")
hypothesis(fit_skeptical_median, "RandomisationiCST < 0")

hypothesis(fit_informative_median, "RandomisationiCST:QCPR_median_groupAbovemedian > 0")
hypothesis(fit_skeptical_median, "RandomisationiCST:QCPR_median_groupAbovemedian > 0")

# ============================================================================
# PART 6: PRIOR AND POSTERIOR PREDICTIVE CHECKS
# ============================================================================

library(patchwork)

prior_check_informative_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_informative_median,
  sample_prior = "only",
  chains  = 4,
  iter    = 1000,
  seed    = 42
)

prior_check_skeptical_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_skeptical_median,
  sample_prior = "only",
  chains  = 4,
  iter    = 1000,
  seed    = 42
)

p_prior_inf <- brms::pp_check(prior_check_informative_median, ndraws = 100) +
  labs(title = "Prior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_inf <- brms::pp_check(fit_informative_median, ndraws = 100) +
  labs(title = "Posterior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_informative <- p_prior_inf | p_post_inf
combined_informative
ggsave("model2_median_ppc_informative.png", combined_informative, width = 10, height = 5, dpi = 300)

p_prior_skep <- brms::pp_check(prior_check_skeptical_median, ndraws = 100) +
  labs(title = "Prior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_skep <- brms::pp_check(fit_skeptical_median, ndraws = 100) +
  labs(title = "Posterior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_skeptical <- p_prior_skep | p_post_skep
combined_skeptical
ggsave("model2_median_ppc_skeptical.png", combined_skeptical, width = 10, height = 5, dpi = 300)

# ============================================================================
# PART 7: INTERACTION BAR CHART (iCST first, TAU Control second)
# ============================================================================

new_data_median <- expand.grid(
  Randomisation = levels(df_complete$Randomisation),
  QCPR_median_group = levels(df_complete$QCPR_median_group),
  c_BASELINE_ADAScog = 0
)

preds_median <- fitted(fit_informative_median, newdata = new_data_median,
                       summary = TRUE, re_formula = NA)

plot_data_median <- cbind(new_data_median, preds_median)
plot_data_median$Randomisation <- factor(plot_data_median$Randomisation,
                                         levels = c("iCST", "TAU Control"))

print(plot_data_median)

interaction_bar_median <- ggplot(plot_data_median,
                                 aes(x = QCPR_median_group, y = Estimate, fill = Randomisation)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                position = position_dodge(width = 0.7), width = 0.15, linewidth = 0.6) +
  scale_fill_manual(values = c("iCST" = "#4ECDC4", "TAU Control" = "#F08080")) +
  labs(
    title = "Interaction: Treatment x Relationship Quality (Median Split)",
    subtitle = "Predicted ADAS-Cog at 26 weeks (with 95% credible intervals)",
    x = "QCPR Group",
    y = "Predicted ADAS-Cog at 26 weeks",
    fill = "Randomisation"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, max(plot_data_median$Q97.5) * 1.1))

interaction_bar_median
ggsave("model2_interaction_barchart_median.png", interaction_bar_median,
       width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 8: COMPREHENSIVE RESULTS TABLE
# ============================================================================

extract_all_params_median <- function(fit, model_label) {
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

params_inf_median <- extract_all_params_median(fit_informative_median, "Informative")
params_skep_median <- extract_all_params_median(fit_skeptical_median, "Skeptical")

full_table_median <- bind_rows(params_inf_median, params_skep_median)

full_table_median$Parameter <- dplyr::recode(full_table_median$Parameter,
                                             "Intercept" = "Intercept",
                                             "c_BASELINE_ADAScog" = "Baseline",
                                             "RandomisationiCST" = "Treatment",
                                             "QCPR_median_groupAbovemedian" = "QCPR Group (Above median)",
                                             "RandomisationiCST:QCPR_median_groupAbovemedian" = "Treatment x QCPR Group",
                                             "sigma" = "Sigma"
)

param_order <- c("Intercept", "Baseline", "QCPR Group (Above median)", "Sigma",
                 "Treatment", "Treatment x QCPR Group")
full_table_median$Parameter <- factor(full_table_median$Parameter, levels = param_order)
full_table_median$Model <- factor(full_table_median$Model, levels = c("Informative", "Skeptical"))

draws_treat_inf_med <- as_draws_df(fit_informative_median)$b_RandomisationiCST
draws_treat_skep_med <- as_draws_df(fit_skeptical_median)$b_RandomisationiCST
draws_int_inf_med <- as_draws_df(fit_informative_median)$`b_RandomisationiCST:QCPR_median_groupAbovemedian`
draws_int_skep_med <- as_draws_df(fit_skeptical_median)$`b_RandomisationiCST:QCPR_median_groupAbovemedian`

prob_lookup_median <- tibble::tibble(
  Model = c("Informative", "Informative", "Skeptical", "Skeptical"),
  Parameter = c("Treatment", "Treatment x QCPR Group", "Treatment", "Treatment x QCPR Group"),
  P_direction = c(
    mean(draws_treat_inf_med < 0),
    mean(draws_int_inf_med > 0),
    mean(draws_treat_skep_med < 0),
    mean(draws_int_skep_med > 0)
  )
)

full_table_median <- full_table_median %>%
  left_join(prob_lookup_median, by = c("Model", "Parameter")) %>%
  arrange(Model, match(Parameter, param_order)) %>%
  mutate(across(c(Estimate, Est.Error, `l-95% CI`, `u-95% CI`, Rhat, P_direction), ~round(.x, 3)),
         across(c(Bulk_ESS, Tail_ESS), ~round(.x, 0)))

print(full_table_median)

write.csv(full_table_median, "model2_median_split_results_table.csv", row.names = FALSE)

full_gt_median <- full_table_median %>%
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
    title = "Model 2 (Median-Split): Posterior Results",
    subtitle = "ADAS-Cog at 26 weeks ~ Baseline ADAS-Cog + Randomisation x QCPR Group"
  ) %>%
  fmt_number(columns = c(Estimate, `Est. Error`, `Lower CrI`, `Upper CrI`,
                         `R-hat`, `Posterior Probability`), decimals = 3) %>%
  fmt_number(columns = c(`Bulk ESS`, `Tail ESS`), decimals = 0) %>%
  sub_missing(columns = `Posterior Probability`, missing_text = "\u2014") %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = Parameter) %>%
  tab_source_note(
    source_note = paste0("N = 254 complete cases. QCPR median = ", round(qcpr_median, 0),
                         ". Reference categories: TAU Control, Below median.")
  )

full_gt_median
gtsave(full_gt_median, "model2_median_split_results_table.png")

## ------------------------------------------------------------------
## Table 1: Baseline Characteristics of Person with Dementia and Caregiver
## N=254 sample - 4 participants with extreme ADAS-Cog change scores
## excluded prior to analysis.
## ------------------------------------------------------------------

library(gtsummary)
library(dplyr)

# ------------------------------------------------------------------
# 0. Confirm sample size and group split
# ------------------------------------------------------------------

nrow(df_complete)   # should be 254
table(df_complete$Randomisation)

# Reorder so Treatment (iCST) displays first
dat_table1 <- df_complete
dat_table1$Randomisation <- factor(dat_table1$Randomisation,
                                   levels = c("iCST", "TAU Control"))

# ------------------------------------------------------------------
# 1. Build the baseline dataset (reusing confirmed codings)
# ------------------------------------------------------------------

table1_data <- dat_table1 %>%
  transmute(
    Randomisation = Randomisation,
    
    # --- Person with dementia ---
    Person_Female = ifelse(qRelGenP == 1, "Male", "Female"),
    Person_Age = cRelAgeP,
    Person_ADAScog = BASELINE_ADAS_20,
    Person_White = ifelse(cRelEthP %in% c(1, 2, 3, 4), "White", "Other"),
    Person_Spouse = ifelse(qRelationship == 1, 1, 0),
    Person_Child = ifelse(qRelationship == 3, 1, 0),
    Person_QCPR = BASELINE_PQCPR_20,
    
    # --- Caregiver ---
    Carer_Female = ifelse(qRelGenC == 1, "Male", "Female"),
    Carer_Age = cRelAgeC,
    Carer_White = ifelse(cRelEthC %in% c(1, 2, 3, 4), "White", "Other"),
    Carer_QCPR = BASELINE_CQCPR_20
  )

# ------------------------------------------------------------------
# 2. Build separate Person and Caregiver summary tables
# ------------------------------------------------------------------

person_vars <- table1_data %>%
  select(Randomisation, Person_Female, Person_Age, Person_ADAScog,
         Person_White, Person_Spouse, Person_Child, Person_QCPR)

carer_vars <- table1_data %>%
  select(Randomisation, Carer_Female, Carer_Age, Carer_White, Carer_QCPR)

tbl_person <- person_vars %>%
  tbl_summary(
    by = Randomisation,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1,
    label = list(
      Person_Female ~ "Female gender, n (%)",
      Person_Age ~ "Age, years: mean (sd)",
      Person_ADAScog ~ "ADAS-Cog: mean (sd)",
      Person_White ~ "Ethnicity: white, n (%)",
      Person_Spouse ~ "Relationship status: spouse, n (%)",
      Person_Child ~ "Relationship status: child, n (%)",
      Person_QCPR ~ "QCPR: mean (sd)"
    ),
    missing = "no"
  )

tbl_carer <- carer_vars %>%
  tbl_summary(
    by = Randomisation,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1,
    label = list(
      Carer_Female ~ "Female gender, n (%)",
      Carer_Age ~ "Age, years: mean (sd)",
      Carer_White ~ "Ethnicity: white, n (%)",
      Carer_QCPR ~ "QCPR: mean (sd)"
    ),
    missing = "no"
  )

# ------------------------------------------------------------------
# 3. Stack into final Table 1
# ------------------------------------------------------------------

final_table1 <- tbl_stack(
  list(tbl_person, tbl_carer),
  group_header = c("Person with dementia", "Caregiver")
) %>%
  modify_header(label = "**Participant**") %>%
  modify_caption("**Table 1. Baseline Characteristics of Person with Dementia and Caregiver**") %>%
  bold_labels()

final_table1

# ------------------------------------------------------------------
# 4. Save
# ------------------------------------------------------------------

final_table1 %>%
  as_gt() %>%
  gt::gtsave("table1_baseline_characteristics_N254.png")

## ============================================================================
## Model 2 (Categorical, Median-Split, N=256): QCPR Group Moderation Analysis
## ADAS_20_FU2 ~ c_BASELINE_ADAScog + Randomisation * QCPR_median_group
## Uses raw BASELINE_CQCPR_20 scores (not z-scored) to define the median split.
## Strict 3 SD outlier rule applied: only qPIN 14201 and 17018 excluded.
## ============================================================================

# ---- LIBRARIES ---------------------------------------------------------
library(brms)
library(dplyr)
library(ggplot2)
library(gt)
library(tidybayes)
library(patchwork)

# ============================================================================
# PART 1: DATA PREPARATION
# ============================================================================

df <- read.csv("~/Downloads/THESIS COMBINED Complete case 26 week follow up data.csv",
               header = TRUE)

# ---- Fix miscoded missing values ----------------------------------------

df$BASELINE_CQCPR_20[df$BASELINE_CQCPR_20 == 0] <- NA

# ---- Filter to complete cases --------------------------------------------

df_complete <- df[complete.cases(df[, c("ADAS_20_FU2",
                                        "c_BASELINE_ADAScog",
                                        "Randomisation",
                                        "BASELINE_CQCPR_20")]), ]

nrow(df_complete)   # should be 258

# ---- Exclude ONLY the 2 cases exceeding a strict 3 SD threshold -----------
# qPIN 14201 (z=4.15) and 17018 (z=4.33)

df_complete <- df_complete[!(df_complete$qPIN %in% c(14201, 17018)), ]

nrow(df_complete)   # should now be 256

# ---- Set factor levels -----------------------------------------------------

df_complete$Randomisation <- factor(df_complete$Randomisation,
                                    levels = c("TAU Control", "iCST"))
table(df_complete$Randomisation)

# ---- Calculate the median (on raw BASELINE_CQCPR_20) and create the
# categorical group -----------------------------------------------------

qcpr_median <- median(df_complete$BASELINE_CQCPR_20, na.rm = TRUE)
qcpr_median

df_complete <- df_complete %>%
  mutate(
    QCPR_median_group = ifelse(BASELINE_CQCPR_20 <= qcpr_median,
                               "Below median", "Above median"),
    QCPR_median_group = factor(QCPR_median_group,
                               levels = c("Below median", "Above median"))
  )

table(df_complete$QCPR_median_group)

# ============================================================================
# PART 2: MODEL SPECIFICATION
# ============================================================================

adas_formula_m2_median <- bf(
  ADAS_20_FU2 | trunc(lb = 0, ub = 60) ~
    c_BASELINE_ADAScog + Randomisation * QCPR_median_group
)

get_prior(adas_formula_m2_median, data = df_complete)

priors_informative_median <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(-1.92, 2), class = b, coef = RandomisationiCST),
  prior(normal(0, 5),     class = b, coef = "QCPR_median_groupAbovemedian"),
  prior(normal(0, 5),     class = b,
        coef = "RandomisationiCST:QCPR_median_groupAbovemedian"),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

priors_skeptical_median <- c(
  prior(normal(20, 4),    class = Intercept),
  prior(normal(0.5, 0.3), class = b, coef = c_BASELINE_ADAScog),
  prior(normal(0, 2),     class = b, coef = RandomisationiCST),
  prior(normal(0, 5),     class = b, coef = "QCPR_median_groupAbovemedian"),
  prior(normal(0, 5),     class = b,
        coef = "RandomisationiCST:QCPR_median_groupAbovemedian"),
  prior(normal(9, 3),     class = sigma, lb = 0)
)

validate_prior(priors_informative_median, adas_formula_m2_median, data = df_complete)
validate_prior(priors_skeptical_median, adas_formula_m2_median, data = df_complete)

# ============================================================================
# PART 3: FIT BOTH MODELS
# ============================================================================

fit_informative_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_informative_median,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model2_median_informative_N256"
)

fit_skeptical_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_skeptical_median,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "icst_model2_median_skeptical_N256"
)

summary(fit_informative_median)
summary(fit_skeptical_median)

# ============================================================================
# PART 4: CONVERGENCE CHECK
# ============================================================================

cat("Max Rhat (Informative):", round(max(brms::rhat(fit_informative_median)), 4), "\n")
cat("Max Rhat (Skeptical):", round(max(brms::rhat(fit_skeptical_median)), 4), "\n")

# ============================================================================
# PART 5: POSTERIOR PROBABILITIES
# ============================================================================

hypothesis(fit_informative_median, "RandomisationiCST < 0")
hypothesis(fit_skeptical_median, "RandomisationiCST < 0")

hypothesis(fit_informative_median, "RandomisationiCST:QCPR_median_groupAbovemedian > 0")
hypothesis(fit_skeptical_median, "RandomisationiCST:QCPR_median_groupAbovemedian > 0")

# ============================================================================
# PART 6: PRIOR AND POSTERIOR PREDICTIVE CHECKS
# ============================================================================

prior_check_informative_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_informative_median,
  sample_prior = "only",
  chains  = 4,
  iter    = 1000,
  seed    = 42
)

prior_check_skeptical_median <- brm(
  formula = adas_formula_m2_median,
  data    = df_complete,
  family  = gaussian(),
  prior   = priors_skeptical_median,
  sample_prior = "only",
  chains  = 4,
  iter    = 1000,
  seed    = 42
)

p_prior_inf <- brms::pp_check(prior_check_informative_median, ndraws = 100) +
  labs(title = "Prior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_inf <- brms::pp_check(fit_informative_median, ndraws = 100) +
  labs(title = "Posterior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_informative <- p_prior_inf | p_post_inf
combined_informative
ggsave("model2_median_ppc_informative_N256.png", combined_informative, width = 10, height = 5, dpi = 300)

p_prior_skep <- brms::pp_check(prior_check_skeptical_median, ndraws = 100) +
  labs(title = "Prior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

p_post_skep <- brms::pp_check(fit_skeptical_median, ndraws = 100) +
  labs(title = "Posterior Predictive Check", x = "ADAS-Cog score", y = "Density") +
  theme_minimal()

combined_skeptical <- p_prior_skep | p_post_skep
combined_skeptical
ggsave("model2_median_ppc_skeptical_N256.png", combined_skeptical, width = 10, height = 5, dpi = 300)

# ============================================================================
# PART 7: INTERACTION BAR CHART (iCST first, TAU Control second)
# ============================================================================

new_data_median <- expand.grid(
  Randomisation = levels(df_complete$Randomisation),
  QCPR_median_group = levels(df_complete$QCPR_median_group),
  c_BASELINE_ADAScog = 0
)

preds_median <- fitted(fit_informative_median, newdata = new_data_median,
                       summary = TRUE, re_formula = NA)

plot_data_median <- cbind(new_data_median, preds_median)
plot_data_median$Randomisation <- factor(plot_data_median$Randomisation,
                                         levels = c("iCST", "TAU Control"))

print(plot_data_median)

interaction_bar_median <- ggplot(plot_data_median,
                                 aes(x = QCPR_median_group, y = Estimate, fill = Randomisation)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                position = position_dodge(width = 0.7), width = 0.15, linewidth = 0.6) +
  scale_fill_manual(values = c("iCST" = "#4ECDC4", "TAU Control" = "#F08080")) +
  labs(
    title = "Interaction: Treatment x Relationship Quality (Median Split)",
    subtitle = "Predicted ADAS-Cog at 26 weeks (with 95% credible intervals), N=256",
    x = "QCPR Group",
    y = "Predicted ADAS-Cog at 26 weeks",
    fill = "Randomisation"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, max(plot_data_median$Q97.5) * 1.1))

interaction_bar_median
ggsave("model2_interaction_barchart_median_N256.png", interaction_bar_median,
       width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 8: COMPREHENSIVE RESULTS TABLE
# ============================================================================

extract_all_params_median <- function(fit, model_label) {
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

params_inf_median <- extract_all_params_median(fit_informative_median, "Informative")
params_skep_median <- extract_all_params_median(fit_skeptical_median, "Skeptical")

full_table_median <- bind_rows(params_inf_median, params_skep_median)

full_table_median$Parameter <- dplyr::recode(full_table_median$Parameter,
                                             "Intercept" = "Intercept",
                                             "c_BASELINE_ADAScog" = "Baseline",
                                             "RandomisationiCST" = "Treatment",
                                             "QCPR_median_groupAbovemedian" = "QCPR Group (Above median)",
                                             "RandomisationiCST:QCPR_median_groupAbovemedian" = "Treatment x QCPR Group",
                                             "sigma" = "Sigma"
)

param_order <- c("Intercept", "Baseline", "QCPR Group (Above median)", "Sigma",
                 "Treatment", "Treatment x QCPR Group")
full_table_median$Parameter <- factor(full_table_median$Parameter, levels = param_order)
full_table_median$Model <- factor(full_table_median$Model, levels = c("Informative", "Skeptical"))

draws_treat_inf_med <- as_draws_df(fit_informative_median)$b_RandomisationiCST
draws_treat_skep_med <- as_draws_df(fit_skeptical_median)$b_RandomisationiCST
draws_int_inf_med <- as_draws_df(fit_informative_median)$`b_RandomisationiCST:QCPR_median_groupAbovemedian`
draws_int_skep_med <- as_draws_df(fit_skeptical_median)$`b_RandomisationiCST:QCPR_median_groupAbovemedian`

prob_lookup_median <- tibble::tibble(
  Model = c("Informative", "Informative", "Skeptical", "Skeptical"),
  Parameter = c("Treatment", "Treatment x QCPR Group", "Treatment", "Treatment x QCPR Group"),
  P_direction = c(
    mean(draws_treat_inf_med < 0),
    mean(draws_int_inf_med > 0),
    mean(draws_treat_skep_med < 0),
    mean(draws_int_skep_med > 0)
  )
)

full_table_median <- full_table_median %>%
  left_join(prob_lookup_median, by = c("Model", "Parameter")) %>%
  arrange(Model, match(Parameter, param_order)) %>%
  mutate(across(c(Estimate, Est.Error, `l-95% CI`, `u-95% CI`, Rhat, P_direction), ~round(.x, 3)),
         across(c(Bulk_ESS, Tail_ESS), ~round(.x, 0)))

print(full_table_median)

write.csv(full_table_median, "model2_median_split_results_table_N256.csv", row.names = FALSE)

full_gt_median <- full_table_median %>%
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
    title = "Model 2 (Median-Split): Posterior Results",
    subtitle = "ADAS-Cog at 26 weeks ~ Baseline ADAS-Cog + Randomisation x QCPR Group (N = 256)"
  ) %>%
  fmt_number(columns = c(Estimate, `Est. Error`, `Lower CrI`, `Upper CrI`,
                         `R-hat`, `Posterior Probability`), decimals = 3) %>%
  fmt_number(columns = c(`Bulk ESS`, `Tail ESS`), decimals = 0) %>%
  sub_missing(columns = `Posterior Probability`, missing_text = "\u2014") %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = Parameter) %>%
  tab_source_note(
    source_note = paste0("N = 256 (strict 3 SD outlier exclusion rule; qPIN 14201 and 17018 excluded). ",
                         "QCPR median = ", round(qcpr_median, 0),
                         ". Reference categories: TAU Control, Below median.")
  )

full_gt_median
gtsave(full_gt_median, "model2_median_split_results_table_N256.png")

library(bayesplot)
library(ggplot2)

trace_plot_median <- mcmc_trace(fit_informative_median,
                                pars = c("b_Intercept", "b_c_BASELINE_ADAScog",
                                         "b_RandomisationiCST",
                                         "b_QCPR_median_groupAbovemedian",
                                         "b_RandomisationiCST:QCPR_median_groupAbovemedian",
                                         "sigma"),
                                facet_args = list(labeller = as_labeller(c(
                                  "b_Intercept" = "Intercept",
                                  "b_c_BASELINE_ADAScog" = "Baseline ADAS-Cog",
                                  "b_RandomisationiCST" = "Treatment",
                                  "b_QCPR_median_groupAbovemedian" = "QCPR Group (Above median)",
                                  "b_RandomisationiCST:QCPR_median_groupAbovemedian" = "Treatment x QCPR Group",
                                  "sigma" = "Sigma"
                                )))) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 11, colour = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

trace_plot_median

ggsave("model2_median_traceplots_informative_N256.png", trace_plot_median,
       width = 10, height = 10, dpi = 300)

## ------------------------------------------------------------------
## Model 2 (Median-Split, N=256): Full Results Table
## Est. Error column removed; formatted in Arial 11.
## ------------------------------------------------------------------

library(dplyr)
library(brms)
library(gt)

extract_all_params_median <- function(fit, model_label) {
  fixed_df <- as.data.frame(summary(fit)$fixed)
  fixed_df$Parameter <- rownames(fixed_df)
  
  sigma_df <- as.data.frame(summary(fit)$spec_pars)
  sigma_df$Parameter <- rownames(sigma_df)
  
  combined <- bind_rows(fixed_df, sigma_df)
  combined$Model <- model_label
  
  combined %>%
    select(Model, Parameter, Estimate, `l-95% CI`, `u-95% CI`,
           Rhat, Bulk_ESS, Tail_ESS)
}

params_inf_median <- extract_all_params_median(fit_informative_median, "Informative")
params_skep_median <- extract_all_params_median(fit_skeptical_median, "Skeptical")

full_table_median <- bind_rows(params_inf_median, params_skep_median)

full_table_median$Parameter <- dplyr::recode(full_table_median$Parameter,
                                             "Intercept" = "Intercept",
                                             "c_BASELINE_ADAScog" = "Baseline",
                                             "RandomisationiCST" = "Treatment",
                                             "QCPR_median_groupAbovemedian" = "QCPR Group (Above median)",
                                             "RandomisationiCST:QCPR_median_groupAbovemedian" = "Treatment x QCPR Group",
                                             "sigma" = "Sigma"
)

param_order <- c("Intercept", "Baseline", "QCPR Group (Above median)", "Sigma",
                 "Treatment", "Treatment x QCPR Group")
full_table_median$Parameter <- factor(full_table_median$Parameter, levels = param_order)
full_table_median$Model <- factor(full_table_median$Model, levels = c("Informative", "Skeptical"))

draws_treat_inf_med <- as_draws_df(fit_informative_median)$b_RandomisationiCST
draws_treat_skep_med <- as_draws_df(fit_skeptical_median)$b_RandomisationiCST
draws_int_inf_med <- as_draws_df(fit_informative_median)$`b_RandomisationiCST:QCPR_median_groupAbovemedian`
draws_int_skep_med <- as_draws_df(fit_skeptical_median)$`b_RandomisationiCST:QCPR_median_groupAbovemedian`

prob_lookup_median <- tibble::tibble(
  Model = c("Informative", "Informative", "Skeptical", "Skeptical"),
  Parameter = c("Treatment", "Treatment x QCPR Group", "Treatment", "Treatment x QCPR Group"),
  P_direction = c(
    mean(draws_treat_inf_med < 0),
    mean(draws_int_inf_med > 0),
    mean(draws_treat_skep_med < 0),
    mean(draws_int_skep_med > 0)
  )
)

full_table_median <- full_table_median %>%
  left_join(prob_lookup_median, by = c("Model", "Parameter")) %>%
  arrange(Model, match(Parameter, param_order)) %>%
  mutate(across(c(Estimate, `l-95% CI`, `u-95% CI`, Rhat, P_direction), ~round(.x, 2)),
         across(c(Bulk_ESS, Tail_ESS), ~round(.x, 0)))

print(full_table_median)

write.csv(full_table_median, "model2_median_split_results_table_N256.csv", row.names = FALSE)

full_gt_median <- full_table_median %>%
  rename(`Lower CrI` = `l-95% CI`,
         `Upper CrI` = `u-95% CI`,
         `Bulk ESS` = Bulk_ESS,
         `Tail ESS` = Tail_ESS,
         `R-hat` = Rhat,
         `Posterior Probability` = P_direction) %>%
  group_by(Model) %>%
  gt() %>%
  tab_header(
    title = "Model 2 (Median-Split): Posterior Results",
    subtitle = "ADAS-Cog at 26 weeks ~ Baseline ADAS-Cog + Randomisation x QCPR Group (N = 256)"
  ) %>%
  fmt_number(columns = c(Estimate, `Lower CrI`, `Upper CrI`,
                         `R-hat`, `Posterior Probability`), decimals = 2) %>%
  fmt_number(columns = c(`Bulk ESS`, `Tail ESS`), decimals = 0) %>%
  sub_missing(columns = `Posterior Probability`, missing_text = "\u2014") %>%
  cols_align(align = "center", columns = everything()) %>%
  cols_align(align = "left", columns = Parameter) %>%
  tab_source_note(
    source_note = paste0("N = 256 (strict 3 SD outlier exclusion rule; qPIN 14201 and 17018 excluded). ",
                         "QCPR median = ", round(qcpr_median, 0),
                         ". Reference categories: TAU Control, Below median.")
  ) %>%
  opt_table_font(font = "Arial") %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = px(13),
    heading.subtitle.font.size = px(11),
    source_notes.font.size = px(9)
  )

full_gt_median
gtsave(full_gt_median, "model2_median_split_results_table_N256.png")

> # ------------------------------------------------------------------
> # 2. Run the Welch's t-test
  > # ------------------------------------------------------------------
> 
  > qcpr_ttest <- t.test(qcpr_change ~ Randomisation, data = dat)
> qcpr_ttest
