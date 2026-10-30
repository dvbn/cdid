# Treatment adjustment: IPW, OR and DR

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

`est_method` decides how each treatment cohort is made comparable with its controls. The same choice applies to every short comparison in the chain. It is separate from the choice of how to handle missing outcomes, which has its own arguments and is described in [Handling missing outcomes](missing-outcomes.md).

## What is adjusted

A chained comparison involves a target cohort g, a control group and two dates s and t at which units are observed. The control group consists of units never treated or units not yet treated, according to `control_group`. With baseline observation assumptions or observation weighting, the outcome changes come from units observed at both dates. Pure imputation and augmented correction instead reconstruct actual changes over the study population before treatment adjustment. A comparison contrasts adjusted cohort and control changes; it can represent a treatment-effect increment or an unrestricted pre-reference differential. The package combines the valid catalog comparisons into the effects ATT(g,t) with the generalized method of moments (GMM); [Parallel trends](parallel-trends.md) explains which comparisons are used.

Treatment adjustment matters when the cohort and its controls differ in characteristics that are related to their outcome changes. Those characteristics enter through the covariate formula `xformla`.

All three methods rely on the same assumptions:

- parallel trends, conditional on the covariates in `xformla` when it has any ([Parallel trends](parallel-trends.md));
- no anticipation beyond the periods allowed by `anticipation`;
- overlap: in every comparison, the controls include units comparable to the cohort, so the probability of belonging to the cohort stays below one;
- a valid way of handling missing outcomes ([Handling missing outcomes](missing-outcomes.md));
- covariates that are fixed within each unit or read at a reference date ([Covariates](#covariates)).

## Treatment weighting (IPW)

`est_method = "ipw"` uses inverse probability weighting. For each distinct risk set, a logit model of membership in the target cohort, against membership in the control group, uses the design of `xformla` and the weights of the requested observation route. Baseline and observation-weighted routes use observed pairs; imputation and augmented routes use the system of reconstructed changes for all study units. Controls are weighted by their fitted odds so that their covariate distribution resembles the cohort's, and the weights within each group are normalized to sum to one before the mean changes are compared. Identical models are fitted once and shared by the comparisons that use them. This chained IPW estimator follows the original chained DiD paper.

What must be correct is the logit model of cohort membership in the population that the comparison represents, after any adjustment for missing outcomes. With units not yet treated as controls, the control group mixes several cohorts, and its odds can be nonlinear in the covariates even when the odds of each cohort are linear.

## Outcome regression (OR)

`est_method = "reg"` models untreated outcome changes. A single weighted least squares regression is fitted to the outcome changes of control units, pooling every comparison that the fit uses, including the comparisons reported before treatment:

```
change(i, e) = alpha(b) + L(e) * z(i)' beta + error      for control units i
```

- `e = (s, t)` is a date pair and `change(i, e)` is the outcome change of unit i from s to t.
- `z(i)` holds the covariates of `xformla`, without the intercept.
- `beta` is one slope vector shared by every comparison and every period.
- `L(e)` is the number of steps from s to t in the stored calendar, the sorted dates present in the data. It is computed from the positions of the two dates, never by subtracting their labels.
- `alpha(b)` is a local intercept for each block b, a distinct combination of date pair and control group. Blocks have intercepts when `xformla` has an intercept.

For a baseline or observation-weighted comparison of cohort g, the estimate is the normalized weighted target-cohort mean of the observed change minus its predicted untreated change, `alpha(b) + L(e) * z(i)' beta`. With imputation or augmented correction, the same treatment adjustment uses reconstructed changes. Treatment OR itself fits no treatment propensity model or actual-treated-change regression; the latter belongs to the separate observation-reconstruction step when requested.

The stored calendar sets the step counts. If the data contain the dates 2001, 2003 and 2004, then:

| Date pair | Positions in the stored calendar | `L(e)` |
| --- | --- | --- |
| 2001 to 2003 | 1 to 2 | 1 |
| 2003 to 2004 | 2 to 3 | 1 |
| 2001 to 2004 | 1 to 3 | 2 |

This specification has three practical consequences.

- A covariate shifts untreated changes by the same amount per calendar step in every period and in every comparison. There is no separate slope for each period. If the influence of a covariate on untreated trends changes over time, the model is misspecified.
- A comparison over two steps has a covariate slope twice that of a comparison over one step, so the regression of short and long changes is internally consistent.
- The shared slope is a restriction on the outcome model. It adds no parallel trends assumption.

The covariate formula determines the fitted model:

| `xformla` | Fitted model |
| --- | --- |
| `NULL` or `~ 1` | a local intercept for each block and no slopes |
| `~ x1 + x2` | local intercepts and shared slopes on `L(e) * x1` and `L(e) * x2` |
| `~ 0 + x1` | a shared slope on `L(e) * x1` and no local intercepts |
| `~ 0` | refused, because `xformla` needs at least one column |

Before fitting, the package looks for columns of the design that are exact combinations of other columns. Such a column is removed only when predictions for every relevant control and treated unit stay the same. Near dependencies that cannot be removed at arithmetic roundoff are retained for the separate rank check. An exact control null direction that fails to reconstruct required treated predictions is refused, as is a pooled design that remains rank deficient. Removing a column only changes the parameterization: the fitted predictions are unchanged, and nothing is assumed about the coefficient of the removed column. The original, retained and removed columns are recorded in `fit$chained$outcome_design` and the route-specific records of `fit$chained$nuisance_diagnostics`. Inspect them with `str()` rather than assuming one common `basis` subfield.

Baseline OR with static covariates needs the design on repeat observers. Dated covariates are read at the route's required reference rows. Imputation and augmented correction need their outcome design on every relevant row of the complete unit table.

This is a difference from `did` only when coefficient restrictions or fitting procedures differ. The same formula notation can describe a shared slope here and separately fitted slopes there; if the complete specification and fitted procedure coincide, OR/DR labels alone cannot explain different estimates. The [balanced-panel decision guide](../from-did.md#why-can-balanced-results-differ-and-which-should-i-use) separates this issue from additional NYT moments and explains the conditional tradeoffs.

## Doubly robust treatment estimation (DR)

`est_method = "dr"`, the default, combines the two models. The pooled control regression predicts each unit's untreated change, and IPW then compares the remaining changes (observed change minus prediction) of the cohort and its weighted controls, with separately normalized weights. For each comparison this gives

```
DR = IPW - (mean design of the cohort - weighted mean design of the controls)' beta
```

where the design of a unit is `L(e) * z(i)`. The local intercept cancels because both sets of weights sum to one. With an intercept only formula, DR therefore equals IPW exactly, including its influence function.

**What DR protects.** Suppose missing outcomes are handled validly and parallel trends, no anticipation, overlap and the regularity conditions hold. Then each comparison is consistent if either the treatment propensity model or the pooled outcome regression is correct. Two qualifications apply. The regression is shared, so it must be correct across all of its pooled blocks: a wrong specification in one block can move the shared slope. The guarantee is a large sample property. DR estimates are not unbiased in finite samples, and when both models are wrong they are generally biased, although the errors can cancel in special cases.

| Treatment propensity model | Outcome regression | IPW | OR | DR |
| --- | --- | --- | --- | --- |
| correct | correct | consistent | consistent | consistent |
| correct | misspecified | consistent | generally biased | consistent |
| misspecified | correct | generally biased | consistent | consistent |
| misspecified | misspecified | generally biased | generally biased | generally biased |

Every entry assumes that the remaining assumptions hold. In particular, treatment DR does not by itself protect against a wrong model for missing outcomes: if observed changes are reweighted or imputed incorrectly, a correct treatment model does not repair the estimate.

## Change imputation uses a different regression

When missing outcomes are handled by change imputation (`selection_est_method = "imputation"`) or by the augmented correction (`selection_est_method = "dr"`), the package fits a second regression. It predicts actual changes, including any treatment effects, separately for each cohort and date pair, and averages the predictions over the study population. That model and its pooling argument `imputation_pooling` are described in [Handling missing outcomes](missing-outcomes.md#change-imputation). The counterfactual OR on this page is a different model, although both use `xformla`. It always shares one slope vector across its blocks and has no pooling argument.

With change imputation, the treatment methods work on the reconstructed changes: IPW fits the treatment logits with the corresponding weights and no counterfactual regression, OR fits the shared counterfactual regression to the imputed changes of the controls and no treatment logit, and DR fits both.

## Covariates

`xformla` is a one sided formula, such as `~ x1 + x2`. Covariates that are constant within each unit (and agree with the unit table when you supply one) are used as unit characteristics.

Covariates that change within a unit, or a formula that uses the time variable, are read at a reference date. Every comparison used for the effects of cohort g reads each unit's values at the cohort's reference period, and a comparison reported before treatment reads them at its earlier date. A required value that is missing is refused with its unit, date and variable; values are never carried forward. [Changing covariates and weights](../reference/options.md#changing-covariates-and-weights) gives the complete rules.

## Custom IPW functions

`est_method` can also be a function with `did`'s panel IPW interface, `f(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...)`. It must return a list with a finite scalar `ATT` and a finite vector or one-column matrix `att.inf.func` with one entry per supplied row, in the supplied order. Chained mode calls it on the units of each comparison, with weights that already include the cohort's response rate, and passes any extra arguments given through `...`.

With `selection = "conditional_mtar"` and observation weighting, the function must also return `log_weight_gradient`, the derivative of `ATT` with respect to the logarithm of each weight; requests without it are refused. The function must implement the IPW criterion and its complete influence function correctly. Return-shape and finiteness checks cannot establish that scientific validity. Custom functions are not available with change imputation or the augmented correction, and `conditional_did_pretest()` returns a result with status `"not_applicable"` for them.

## Example: adjust for a covariate in a complete panel

In this simulated design a fixed covariate `x` raises both the probability of treatment and the untreated trend, so parallel trends hold only conditional on `x`.


```r
set.seed(20261025)

n_units <- 2000
n_periods <- 6

x <- rnorm(n_units)
first_treat <- ifelse(runif(n_units) < plogis(0.8 * x), 3, 0)
unit_effect <- rnorm(n_units)

cov_panel <- data.frame(
  id = rep(seq_len(n_units), each = n_periods),
  period = rep(seq_len(n_periods), times = n_units),
  first_treat = rep(first_treat, each = n_periods),
  x = rep(x, each = n_periods)
)

treated_now <- cov_panel$first_treat > 0 &
  cov_panel$period >= cov_panel$first_treat
effect <- ifelse(treated_now,
                 1 + 0.25 * (cov_panel$period - cov_panel$first_treat), 0)

# Untreated outcomes rise by 0.5 + 0.3 * x per period.
cov_panel$y <- rep(unit_effect, each = n_periods) +
  (0.5 + 0.3 * cov_panel$x) * cov_panel$period +
  effect +
  rnorm(nrow(cov_panel))

fit_unadjusted <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  data = cov_panel, control_group = "nevertreated",
  est_method = "ipw", bstrap = FALSE, cband = FALSE
)
fit_ipw <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = cov_panel, control_group = "nevertreated",
  est_method = "ipw", bstrap = FALSE, cband = FALSE
)
fit_or <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = cov_panel, control_group = "nevertreated",
  est_method = "reg", bstrap = FALSE, cband = FALSE
)
fit_dr <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = cov_panel, control_group = "nevertreated",
  est_method = "dr", bstrap = FALSE, cband = FALSE
)

fits <- list(unadjusted = fit_unadjusted, ipw = fit_ipw, or = fit_or, dr = fit_dr)
simple <- lapply(fits, cdid::aggte, type = "simple")
data.frame(
  fit = names(simple),
  estimate = vapply(simple, function(a) a$overall.att, numeric(1)),
  std.error = vapply(simple, function(a) a$overall.se, numeric(1))
)
```

Output:

```text
fit estimate  std.error
unadjusted unadjusted 1.846435 0.05929326
ipw               ipw 1.348377 0.06046122
or                 or 1.324931 0.05151023
dr                 dr 1.326676 0.05714047
```

`bstrap = FALSE, cband = FALSE` requests analytical standard errors with pointwise intervals, which keeps the comparison simple; `aggte()` inherits these settings.

**What the design implies.** The average effect over periods 3 to 6 is 1.375, as in [Get started](../get-started.md). Cohort 3 has larger values of `x` on average, so its untreated outcomes rise faster than those of the units never treated. The unadjusted fit is therefore expected to overstate the effects, by more at later periods. With `xformla = ~ x` the treatment logit is correctly specified, because the probability of treatment is a logit in `x`, and so is the pooled regression, because the untreated change over a pair of dates is `0.5 * L(e) + 0.3 * L(e) * x`: a block intercept plus a shared slope on `L(e) * x`. IPW, OR and DR therefore all target the effects, and their estimates will differ somewhat in a finite sample.

A variation shows what DR protects. Replace `0.3 * cov_panel$x` by `0.3 * cov_panel$x^2` in the outcome and keep `xformla = ~ x`. The outcome regression is then misspecified while the treatment logit stays correct, so the design implies that OR is generally biased while IPW and DR still target the effects.

## Where the result records the choice

| Location | Contents |
| --- | --- |
| `fit$DIDparams$est_method`, `fit$DIDparams$xformla` | the requested method and the formula with its environment |
| `fit$chained$estimator` | the estimator label: `"chained IPW-GMM"`, `"chained outcome regression/control-mean GMM"` or `"chained doubly robust GMM"`, with a suffix for change imputation or the augmented correction |
| `fit$chained$selection` | for OR and DR, the outcome and treatment model descriptions, the step convention, the assumptions and model counts such as `n_regression_fits`, `n_regression_blocks`, `n_local_intercepts`, `n_shared_slopes` and `n_treatment_fits` |
| `fit$chained$nuisance_diagnostics` | for OR and DR, the shared slopes, block intercepts, fitted comparisons and step convention |
| `fit$chained$outcome_design` | original, retained and removed design columns |

`tidy()` and `ggdid()` show the estimator label for OR and DR fits and for the imputation and augmented routes.

---

**Reading path, step 7 of 20** · Previous: [Prepare your panel and unit table](../help/prepare-data.md) · Next: [Choose IPW, OR or DR in an incomplete panel](../examples/treatment-choice.md)

<details>
<summary>All steps of the reading path</summary>

| Step | Page |
| --- | --- |
| 1 | [Get started](../get-started.md) |
| 2 | [Moving from did](../from-did.md) |
| 3 | [Incomplete panels](../examples/incomplete-panels.md) |
| 4 | [Identify effects when the observed units change composition](../examples/identification.md) |
| 5 | [How to compare precision](../examples/precision.md) |
| 6 | [Prepare your panel and unit table](../help/prepare-data.md) |
| 7 | **Treatment adjustment** (this page) |
| 8 | [Choose IPW, OR or DR in an incomplete panel](../examples/treatment-choice.md) |
| 9 | [Handling missing outcomes](missing-outcomes.md) |
| 10 | [Sampling information for the observation model](sampling-information.md) |
| 11 | [Choose an observation model](../examples/observation-models.md) |
| 12 | [Impute actual changes](../examples/change-imputation.md) |
| 13 | [Parallel trends](parallel-trends.md) |
| 14 | [Add PT-all restrictions](../examples/pt-all-comparison.md) |
| 15 | [GMM weighting and covariance regularization](gmm.md) |
| 16 | [Robustness](../examples/robustness.md) |
| 17 | [Inference](../examples/inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
