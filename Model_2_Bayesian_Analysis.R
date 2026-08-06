# ------------------------------------------------------------------
> # 4. Standardise BASELINE_CQCPR_20 (QPRC) into a z-score
  > # ------------------------------------------------------------------
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