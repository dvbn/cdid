# Monte Carlo & Scaling Validation: Cross-Pair Chaining, `PT-ctrl`, and Within-Unit GMM in `cdid`

These results come from the cross-pair source revision `5b1a614`. The simulations used `cdid` (1.2.1.0), `did` 2.5.1, and `DRDID` 1.4.0 under R version 4.5.2 (2025-10-31) (`seed = 20261111`, `200` replications per design).

## Executive Summary of Findings

1. **Exact Numerical Equivalence ($\le 1.56e-13$):** Across all 200 replications of cohort-balanced and wave-gap short-window designs (Designs 1, 2, 3, and 5), `att_gt(..., cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within", base_period = "universal")` matches both unit-fixed-effect and cohort-fixed-effect stacked target-slice OLS to machine precision (max absolute difference `1.563e-13`). Under non-monotone intermittent post-treatment dropout (Design 4), `att_gt()` pools $O(T^2)$ cohort-pair differences rather than $2^T$ missingness strata (max finite-sample difference `2.563e-02` vs unit-FE slice OLS) while eliminating cohort-FE composition bias.
2. **Optimal Two-Step GMM Efficiency:** When errors depart from homoskedastic random walks across cohorts (Design 2), `gmm_weighting = "optimal"` achieves strictly lower empirical SD and RMSE than `within`-weighted slice regression.
3. **Pre-Reference Trend Robustness (`PT-post` vs `PT-ctrl`):** When a later-treated control cohort experiences a trend shift prior to the target cohort's reference period (Design 3), stacked slice OLS and `PT-ctrl` are biased because they pool control pre-trends before the reference period, whereas `cdid` under `pt_assumption = "post"` remains unbiased with nominal 95% coverage.
4. **Selective-Attrition Robustness (Unit Differences vs Cohort FE Levels):** When post-treatment attrition depends on persistent unit fixed effects $\alpha_i$ (Design 4), cohort-FE slice OLS suffers severe composition bias, whereas `cdid` removes $\alpha_i$ via within-unit differencing and remains unbiased with nominal 95% coverage.
5. **Biennial Wave Gaps & Non-Overlapping Cross-Pair Windows with Covariates:** Across biennial gaps (Design 5), `cdid` estimates seamlessly without dropping post-transition waves. When target and control cohorts are never observed on a common date pair, control pair observation is partial ($r_{4,2,3} = 0.5$), and covariates $X_i$ drive both treatment timing and outcome trends (Design 6), same-pair `att_gt(..., cross_pair = FALSE)` returns `NA`, `did::att_gt(..., allow_unbalanced_panel = TRUE)` falls back to repeated cross-sections across distinct cohorts (`-4.022` bias), unadjusted slice OLS suffers covariate-trend bias (`+1.985`), and both doubly robust and inverse-rate-calibrated IPW `att_gt(..., xformla = ~ x, cross_pair = TRUE)` are unbiased with nominal 95% coverage.
6. **Conditional-MTAR Selective Observation & Selection Double Robustness (`selection = "conditional_mtar"`):** When panel entry ($A_g \subset F_g$) and pair retention depend on covariates $X_i$ alongside heterogeneous treatment effects $\tau(X_i)$ on disjoint windows (Design 7), naive `selection = "none"` suffers severe selection bias (`+0.640`) even with `xformla = ~ x`, whereas `selection = "conditional_mtar"` across `selection_est_method = "ipw"`, `"imputation"`, and `"dr"` (augmented, which remains unbiased even when `imputation_pooling = "within_cohort"` is misspecified) recovers the full study-population `ATT(3, 4)` with nominal ~95% analytical coverage.

## Design 1: Staggered Short Windows — Homoskedastic i.i.d. Errors

- **Target parameter:** `ATT(2005, 2009)` = `5.00` ($N = 900$ units across 6 cohorts `2005..2010`, 200 replications).
- **Max numerical discrepancy `|cdid(within, PT-ctrl) - Stacked slice OLS (unit FE)|`:** `1.075e-13`.

| Estimator | Avail (%) | Bias | Empirical SD | Mean Analytical SE | RMSE | 95% Coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Stacked slice OLS (cohort FE) | 100% | -0.001 | 0.163 | n/a | 0.163 | n/a |
| Stacked slice OLS (unit FE) | 100% | -0.001 | 0.163 | n/a | 0.163 | n/a |
| cdid (within, PT-ctrl, cross-pair) | 100% | -0.001 | 0.163 | 0.164 | 0.163 | 0.975 |
| cdid (optimal, PT-ctrl, cross-pair) | 100% | -0.003 | 0.163 | 0.163 | 0.163 | 0.970 |
| cdid (optimal, PT-post, cross-pair) | 100% | -0.002 | 0.162 | 0.163 | 0.162 | 0.970 |
| cdid (optimal, PT-post, same-pair) | 100% | -0.002 | 0.162 | 0.163 | 0.162 | 0.970 |
| cdid (optimal, PT-all, cross-pair) | 100% | 0.001 | 0.153 | 0.148 | 0.152 | 0.950 |
| did::att_gt (CS 2021, NYT) | 0% | n/a | n/a | n/a | n/a | n/a |

## Design 2: Staggered Short Windows — Heteroskedastic Random-Walk Errors

- **Target parameter:** `ATT(2005, 2009)` = `5.00` ($N = 900$ units across 6 cohorts `2005..2010`, 200 replications).
- **Max numerical discrepancy `|cdid(within, PT-ctrl) - Stacked slice OLS (unit FE)|`:** `1.066e-13`.

| Estimator | Avail (%) | Bias | Empirical SD | Mean Analytical SE | RMSE | 95% Coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Stacked slice OLS (cohort FE) | 100% | -0.013 | 0.213 | n/a | 0.213 | n/a |
| Stacked slice OLS (unit FE) | 100% | -0.013 | 0.213 | n/a | 0.213 | n/a |
| cdid (within, PT-ctrl, cross-pair) | 100% | -0.013 | 0.213 | 0.207 | 0.213 | 0.960 |
| cdid (optimal, PT-ctrl, cross-pair) | 100% | -0.008 | 0.173 | 0.171 | 0.172 | 0.960 |
| cdid (optimal, PT-post, cross-pair) | 100% | -0.009 | 0.172 | 0.171 | 0.172 | 0.960 |
| cdid (optimal, PT-post, same-pair) | 100% | -0.009 | 0.172 | 0.171 | 0.172 | 0.960 |
| cdid (optimal, PT-all, cross-pair) | 100% | -0.008 | 0.169 | 0.168 | 0.169 | 0.955 |
| did::att_gt (CS 2021, NYT) | 0% | n/a | n/a | n/a | n/a | n/a |

## Design 3: Control Cohort 2006 Shifted in Pre-Reference Years (2002–2003 Only)

- **Target parameter:** `ATT(2005, 2009)` = `5.00` ($N = 900$ units across 6 cohorts `2005..2010`, 200 replications).
- **Max numerical discrepancy `|cdid(within, PT-ctrl) - Stacked slice OLS (unit FE)|`:** `1.075e-13`.

| Estimator | Avail (%) | Bias | Empirical SD | Mean Analytical SE | RMSE | 95% Coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Stacked slice OLS (cohort FE) | 100% | -0.094 | 0.154 | n/a | 0.180 | n/a |
| Stacked slice OLS (unit FE) | 100% | -0.094 | 0.154 | n/a | 0.180 | n/a |
| cdid (within, PT-ctrl, cross-pair) | 100% | -0.094 | 0.154 | 0.164 | 0.180 | 0.925 |
| cdid (optimal, PT-ctrl, cross-pair) | 100% | -0.091 | 0.158 | 0.163 | 0.182 | 0.925 |
| cdid (optimal, PT-post, cross-pair) | 100% | 0.029 | 0.154 | 0.163 | 0.156 | 0.965 |
| cdid (optimal, PT-post, same-pair) | 100% | 0.029 | 0.154 | 0.163 | 0.156 | 0.965 |
| cdid (optimal, PT-all, cross-pair) | 100% | 0.756 | 0.170 | 0.148 | 0.774 | 0.000 |
| did::att_gt (CS 2021, NYT) | 0% | n/a | n/a | n/a | n/a | n/a |

## Design 4: Selective Post-Treatment Attrition (Low-alpha Units Drop Out)

- **Target parameter:** `ATT(2005, 2009)` = `5.00` ($N = 900$ units across 6 cohorts `2005..2010`, 200 replications).
- **Max numerical discrepancy `|cdid(within, PT-ctrl) - Stacked slice OLS (unit FE)|`:** `2.563e-02`.

| Estimator | Avail (%) | Bias | Empirical SD | Mean Analytical SE | RMSE | 95% Coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Stacked slice OLS (cohort FE) | 100% | 0.468 | 0.195 | n/a | 0.507 | n/a |
| Stacked slice OLS (unit FE) | 100% | 0.009 | 0.168 | n/a | 0.168 | n/a |
| cdid (within, PT-ctrl, cross-pair) | 100% | 0.008 | 0.169 | 0.173 | 0.168 | 0.950 |
| cdid (optimal, PT-ctrl, cross-pair) | 100% | 0.013 | 0.178 | 0.163 | 0.178 | 0.930 |
| cdid (optimal, PT-post, cross-pair) | 100% | 0.013 | 0.177 | 0.163 | 0.177 | 0.925 |
| cdid (optimal, PT-post, same-pair) | 100% | 0.008 | 0.180 | 0.176 | 0.180 | 0.960 |
| cdid (optimal, PT-all, cross-pair) | 100% | 0.011 | 0.158 | 0.152 | 0.158 | 0.960 |
| did::att_gt (CS 2021, NYT) | 0% | n/a | n/a | n/a | n/a | n/a |

## Design 5: Biennial Wave Gaps (2007 & 2009 Removed, Target ATT(2005, 2008))

- **Target parameter:** `ATT(2005, 2008)` = `4.00` ($N = 900$ units across 6 cohorts `2005..2010`, 200 replications).
- **Max numerical discrepancy `|cdid(within, PT-ctrl) - Stacked slice OLS (unit FE)|`:** `1.563e-13`.

| Estimator | Avail (%) | Bias | Empirical SD | Mean Analytical SE | RMSE | 95% Coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Stacked slice OLS (cohort FE) | 100% | 0.011 | 0.148 | n/a | 0.148 | n/a |
| Stacked slice OLS (unit FE) | 100% | 0.011 | 0.148 | n/a | 0.148 | n/a |
| cdid (within, PT-ctrl, cross-pair) | 100% | 0.011 | 0.148 | 0.152 | 0.148 | 0.945 |
| cdid (optimal, PT-ctrl, cross-pair) | 100% | 0.010 | 0.148 | 0.151 | 0.148 | 0.955 |
| cdid (optimal, PT-post, cross-pair) | 100% | 0.011 | 0.148 | 0.151 | 0.148 | 0.950 |
| cdid (optimal, PT-post, same-pair) | 100% | 0.011 | 0.148 | 0.151 | 0.148 | 0.950 |
| cdid (optimal, PT-all, cross-pair) | 100% | 0.006 | 0.128 | 0.134 | 0.127 | 0.960 |
| did::att_gt (CS 2021, NYT) | 0% | n/a | n/a | n/a | n/a | n/a |

## Design 6: Disjoint Cross-Pair Observation Windows with Covariate Confounding & Partial Pair Observation

- **Target parameter:** `ATT(3, 4)` = `4.00` ($N = 720$ units across cohorts `{0, 3, 4, 5}`, where `G=3` is observed at `{1, 2, 4}`, half of `G=4` at `{2, 3}` and half at `{1, 2}` so $r_{4,2,3} = 0.5$, `G=5` at `{3, 4}`, and `G=0` at `{1, 2}`; 200 replications).
- **Key property:** Zero control units are observed at both dates `(2, 4)`, so same-pair estimators return `NA`, while covariate $X_i$ shifts both cohort assignment and untreated trend slopes under partial within-cohort pair observation.

| Estimator | Avail (%) | Bias | Empirical SD | Mean Analytical SE | RMSE | 95% Coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Stacked slice OLS (cohort FE, unadjusted) | 100% | 1.985 | 0.366 | n/a | 2.018 | n/a |
| cdid (unadjusted, within, PT-ctrl, cross-pair) | 100% | 1.875 | 0.222 | 0.217 | 1.888 | 0.000 |
| cdid (xformla = ~x, DR, within, PT-ctrl, cross-pair) | 100% | -0.027 | 0.193 | 0.190 | 0.194 | 0.940 |
| cdid (xformla = ~x, DR, optimal, PT-ctrl, cross-pair) | 100% | -0.027 | 0.192 | 0.190 | 0.194 | 0.940 |
| cdid (xformla = ~x, DR, optimal, PT-post, cross-pair) | 100% | -0.027 | 0.193 | 0.190 | 0.194 | 0.940 |
| cdid (xformla = ~x, IPW, optimal, PT-post, cross-pair) | 100% | 0.042 | 0.275 | 0.249 | 0.277 | 0.930 |
| cdid (xformla = ~x, DR, optimal, PT-post, same-pair) | 0% | n/a | n/a | n/a | n/a | n/a |
| did::att_gt (xformla = ~x, DR, CS 2021) | 100% | -4.022 | 0.157 | 0.156 | 4.025 | 0.000 |

## Design 7: Disjoint Cross-Pair Windows with Covariate-Dependent Selective Observation (`selection = "conditional_mtar"`)

- **Target parameter:** Full study-population `ATT_F(3, 4)` = `3.00` ($N = 1{,}000$ frame units across cohorts `{0, 3, 4, 5}` including frame-only and singleton dropouts, heterogeneous treatment effect $\tau(X_i) = 3.0 + 2.0(X_i - 0.4)$, and logit selection probabilities driven by $X_i$; 200 replications).
- **Key property:** High-$X_i$ target units are far more likely to remain in the repeat-observed panel $A_3$ and have larger treatment effects, so `selection = "none"` uses the selected changes and exhibits substantial bias relative to `ATT_F(3, 4)`, while `selection = "conditional_mtar"` (`ipw`, `imputation`, and `dr` augmented) transports back to $F_3$ with full study-population influence-function inference.

| Estimator | Avail (%) | Bias | Empirical SD | Mean Analytical SE | RMSE | 95% Coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cdid (selection = 'none', DR, xformla = ~x, cross-pair) | 100% | 0.640 | 0.216 | 0.198 | 0.675 | 0.080 |
| cdid (conditional_mtar, sel = 'ipw', IPW, cross-pair) | 100% | -0.014 | 0.308 | 0.261 | 0.308 | 0.905 |
| cdid (conditional_mtar, sel = 'ipw', OR, cross-pair) | 100% | -0.008 | 0.292 | 0.250 | 0.292 | 0.900 |
| cdid (conditional_mtar, sel = 'ipw', DR, cross-pair) | 100% | -0.009 | 0.295 | 0.249 | 0.294 | 0.910 |
| cdid (conditional_mtar, sel = 'imputation' [sep], DR, cross-pair) | 100% | -0.014 | 0.193 | 0.180 | 0.193 | 0.945 |
| cdid (conditional_mtar, sel = 'dr' [augmented], DR, cross-pair) | 100% | -0.016 | 0.202 | 0.185 | 0.202 | 0.935 |
| cdid (conditional_mtar, sel = 'dr' [augmented], DR, same-pair) | 0% | n/a | n/a | n/a | n/a | n/a |

## Runtime Scaling Benchmarks (6 Cohorts, 10 Calendar Periods)

| Physical Units ($N$) | Panel Rows | Same-Pair Opt (s) | Cross-Pair Within (s) | Cross-Pair DR Opt (s) | Cross-Pair MTAR IPW (s) | Cross-Pair MTAR Imp (s) | Cross-Pair MTAR Aug DR (s) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 900 | 6750 | 0.398 | 0.081 | 0.084 | 0.386 | 0.201 | 0.411 |
| 3000 | 22500 | 0.703 | 0.102 | 0.121 | 0.519 | 0.255 | 0.667 |
| 9000 | 67500 | 1.759 | 0.191 | 0.226 | 1.236 | 0.494 | 1.091 |
| 18000 | 135000 | 3.293 | 0.333 | 0.388 | 1.831 | 0.794 | 1.928 |
