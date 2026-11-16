# Impute actual changes for missing outcomes

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage. OR/DR calls shown here explicitly retain pooled slopes; the current default uses separate comparison regressions.

Observation weighting reweights the changes of units observed at both dates. Change imputation instead predicts the actual change of every unit from its covariates, and the augmented correction combines prediction with weighting. This example compares the three corrections on the panel of [Choose an observation model for missing outcomes](observation-models.md), where every model they need is correctly specified by construction. [Handling missing outcomes](../methods/missing-outcomes.md#change-imputation) defines the methods.



## Prerequisite

Run the code under [Generate the panel](observation-models.md#generate-the-panel) first. This page uses its `sim`, `panel` and `units` exactly as generated there, with the same missing outcomes, and draws no new data. The fit `fit_none` from that page is an optional baseline only when fitted to these same objects: its assumption fails in this design, so it offers no protection.

```r
# Objects from the first example; this page draws no new data.
stopifnot(
  is.data.frame(sim), is.data.frame(panel), is.data.frame(units),
  nrow(units) == 2000, nrow(panel) == 2000 * 6,
  identical(!is.na(panel$y), sim$observed)
)
```

## Two outcome regressions

Change imputation fits a model `j` of actual changes. It is not the counterfactual regression `m` of treatment OR and DR.

| | Actual change model `j` | Counterfactual regression `m` |
| --- | --- | --- |
| Layer | handling missing outcomes: `selection_est_method = "imputation"` or `"dr"` | treatment adjustment: `est_method = "reg"` or `"dr"` |
| Predicts | the actual change of a cohort over a pair of dates, effects included | the untreated change |
| Fitted to | units of the cohort observed at both dates | control changes, reconstructed under imputation and the augmented correction |
| Slopes | shared as `imputation_pooling` declares | separate by comparison by default, or shared across blocks with `outcome_pooling = "pooled"` |

Both read `xformla`. Here `~ x` defines `j` and both treatment models of `est_method = "dr"`, the treatment logit `p` and `m`. `selection_xformla` defines only the observation model `q`.

## How change imputation works

For each cohort k and each pair of dates e that the estimator needs, `j` is fitted by weighted least squares to the units of cohort k observed at both dates:

```
j(k, e)(x) = alpha(k, e) + L(e) * x * beta(r)
```

`L(e)` counts the steps between the two dates in the stored calendar, so `L(e) = t - s` for dates 1 to 6. `alpha(k, e)` is an intercept for each cohort and pair, present because `~ x` has one, and `beta(r)` is the slope of pooling group r. The default `imputation_pooling = "separate"` gives each cohort and pair its own intercept and slope.

The fitted model is evaluated for every unit of cohort k in `units`, including units whose change was observed, and the average of these predictions replaces the observed mean. You fill in nothing; `y` stays `NA`. A unit observed at one date supplies no change. Its prediction is the whole change predicted from `x`; the method does not condition on its one observed outcome or impute a missing level from it. This is one imputation strategy, prediction and averaging of changes, not imputation of missing outcomes in general.

Because `j` predicts actual changes, it must represent how the effect varies with `x`. For a pair that spans the start of treatment, the mean change of cohort 3 has slope `0.3 * L(e) + 0.5` on `x`; for every other pair of either cohort the slope is `0.3 * L(e)`. Thus the duration scaled coefficient `beta(r)` must be `0.3 + 0.5 / L(e)` for treatment crossing pairs of cohort 3 and `0.3` otherwise. With binary `x` and separate cells, `j` is correct. The mean untreated change, `0.1 * L(e) + 0.3 * L(e) * x`, makes `m` correct with a block intercept and one slope on `L(e) * x`, and `p` is correct because the cohort is drawn independently of `x`. These are algebraic properties of the design, not fit results or guarantees.

## Three corrected fits

The three fits keep treatment DR and every setting of the first example; only the correction changes.

```r
fit_weighting <- cdid::att_gt(
  outcome_pooling = "pooled",
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x, data = panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  pt_assumption = "post", base_period = "varying",
  selection = "conditional_mtar",
  selection_est_method = "ipw",
  selection_model = "independent_period_logit",
  selection_pooling = "separate",
  selection_universe = "F",
  bstrap = FALSE, cband = FALSE
)

fit_imputation <- cdid::att_gt(
  outcome_pooling = "pooled",
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x, data = panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  pt_assumption = "post", base_period = "varying",
  selection = "conditional_mtar",
  selection_est_method = "imputation",
  selection_universe = "F",
  imputation_pooling = "separate",
  bstrap = FALSE, cband = FALSE
)

fit_augmented <- cdid::att_gt(
  outcome_pooling = "pooled",
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x, data = panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  pt_assumption = "post", base_period = "varying",
  selection = "conditional_mtar",
  selection_est_method = "dr",
  selection_model = "independent_period_logit",
  selection_pooling = "separate",
  selection_universe = "F",
  imputation_pooling = "separate",
  bstrap = FALSE, cband = FALSE
)
```

`fit_weighting` repeats the specification of `fit_date`. `fit_imputation` passes no `selection_model`, `selection_pooling` or `selection_xformla`: under pure imputation they are inactive and never evaluated, and no observation probability is fitted. `fit_augmented` uses both groups of arguments, and both models are active.

## The augmented correction

For each needed cell c, a cohort and pair, every unit of the cohort receives the reconstructed change

```
Z(c) = j(c)(x) + S(e) * { change(e) - j(c)(x) } / max(q(c)(x), 1e-6)
```

where `S(e)` is 1 for units observed at both dates and 0 otherwise, and `q(c)` is the fitted probability of observing both dates, here a product of two fitted date probabilities. For cohort 3, `j(c)` is the model of its actual change, effects included, not the counterfactual regression `m`; `m` is fitted afterward to the reconstructed changes of the controls. The correction term is not renormalized. Since version 1.1.0.0 its inverse observation denominator is floored at `1e-6` ([Probability safeguards](../methods/missing-outcomes.md#probability-safeguards)). The result is one estimator, not an average of the other two.

The examples keep `est_method = "dr"` so only the missing-outcome adjustment changes. You can instead combine any of the three adjustments with treatment `"ipw"` or `"reg"`. The augmented correction still fits both `q` and `j`; treatment IPW additionally fits only `p`, and treatment OR only `m`. Their treatment-model requirements differ as shown in [Correctness conditions across both layers](robustness.md#correctness-conditions-across-both-layers).

## Which models must be correct

| Fit | Missing outcome layer | Treatment layer |
| --- | --- | --- |
| `fit_weighting` | `q` correct in every needed cell, with independence of the two dates | `p` or `m` correct |
| `fit_imputation` | `j` correct in every needed cell | `p` or `m` correct |
| `fit_augmented` | for each needed cell, `q` or `j` correct, possibly different across cells | `p` or `m` correct |

Every fit also needs conditional mean trends missing at random given a scientifically sufficient conditioning set W, here `x`; PT-post and no anticipation; positive cohort masses, positive observation probabilities where `q` is fitted and positive training mass with prediction support where `j` is; valid moments for every comparison used; and regular finite fits of every active model, alongside the usual large sample regularity conditions. For the augmented route, W includes the covariates in `xformla` and any variables only in `selection_xformla`; a correct `j` must represent the actual-change mean on that W. Correctness must hold across the cells sharing the relevant slopes, both endpoints of an active date model and all blocks of a shared counterfactual regression. Different predictor lists need separately sufficient treatment and observation restrictions. Pure imputation needs a correct `j` even when `p` and `m` are both correct: treatment DR does not generally repair a wrong imputer. These conditions are sufficient, not necessary, and they concern large sample behavior. There is no general finite sample unbiasedness guarantee, and no fit is guaranteed the lowest standard error. Here every model is correct, so all three fits are expected to target the effects. [Compare model robustness](robustness.md#correctness-conditions-across-both-layers) gives the complete conditions.

## Three kinds of slope sharing

`imputation_pooling` shares slopes of `j`, and `selection_pooling` shares slopes of `q`; under the augmented correction they act independently. The counterfactual regression `m` uses separate comparison slopes by default; `outcome_pooling = "pooled"` shares them across blocks. Its pooling choice is independent of the other two.

| `imputation_pooling` | Slopes of `j` | Here |
| --- | --- | --- |
| `"separate"` (default) | one per cohort and pair | correct |
| `"within_cohort"` | one per cohort across its pairs | misspecified: pairs of cohort 3 that span the start of treatment have the extra slope 0.5 on `x` |
| `"across_cohorts"` | one for all cells | misspecified for the same reason |

Each cohort and pair keeps its intercept whenever `xformla` has one. Shared slopes are a restriction you choose, and they can be wrong, for example when the effect varies with covariates, as here. Pooling is not changed automatically after a failure. A pooling group whose reduced design remains rank deficient is refused; exact aliases may first be removed only when the required predictions are preserved.

## Read the results

```r
fits <- list(weighting = fit_weighting, imputation = fit_imputation,
             augmented = fit_augmented)
# Optional: add this baseline only if fit_none was fitted to these same
# panel/units objects in the preceding page, not an unrelated session object.
# fits <- c(list(none = fit_none), fits)

simple <- lapply(fits, cdid::aggte, type = "simple")
data.frame(
  handling = names(simple),
  estimate = vapply(simple, function(a) a$overall.att, numeric(1)),
  std.error = vapply(simple, function(a) a$overall.se, numeric(1))
)

# Known value of the simple aggregate for these generated units.
1.375 + 0.5 * mean(units$x[units$g == 3])

fit_imputation$chained$selection[c("method", "imputation_pooling")]
fit_augmented$chained$selection[c("method", "selection_model", "selection_pooling",
                                  "imputation_pooling", "selection_formula_inherited")]
```

<details>
<summary>Console output</summary>

```text
handling estimate  std.error
weighting   weighting 1.623850 0.03817623
imputation imputation 1.628628 0.03800096
augmented   augmented 1.628628 0.03800096
[1] 1.626773
$method
[1] "imputation"

$imputation_pooling
[1] "separate"

$method
[1] "augmented"

$selection_model
[1] "independent_period_logit"

$selection_pooling
[1] "separate"

$imputation_pooling
[1] "separate"

$selection_formula_inherited
[1] TRUE
```

</details>

With one treated cohort, the simple aggregate is the average of ATT(3, 3) to ATT(3, 6), `1.375 + 0.5 * E[x | g = 3]`, which is 1.625 in the population. The known value uses the generated `x` of cohort 3; it is a target, not a fitted aggregate, and it involves no estimated cohort shares. `method` names the route; `selection_formula_inherited` reports whether `selection_xformla` was inherited from `xformla`, which these calls leave at `NULL`. The nuisance diagnostics of `fit_imputation` list the imputation cells and pooling groups, the counterfactual regression and the treatment logits, with no observation probability model. `summary()` adds a line describing either route.

## Refusals and limits

- Pure imputation fits no observation logit, so no such logit can cause a refusal, and it can run when every unit of a cell is observed. Its own models must still pass their support and rank checks, and the treatment logits and counterfactual regression keep theirs.
- The augmented correction fits both missing-outcome models and the models active under the chosen treatment method. A failure in any active model is refused, including an active observation cell with its own free intercept in which every unit, or none, is observed. Nothing falls back to another correction.
- A route's label does not make every cell available. A cell that the comparisons cannot support is reported as unavailable with its reason.
- Version 1.1.0.0 applies fixed treatment-propensity caps and inverse-observation floors without dropping observations; binding safeguards can alter the original target and robustness claims ([Probability safeguards](../methods/missing-outcomes.md#probability-safeguards)). A multiplier bootstrap, when requested, reuses the influence functions and refits nothing. Missing outcomes or changes are not randomly drawn: this is not multiple imputation and uses no combination by Rubin's rules. The multiplier bootstrap does draw random signs for inference.
- Standard errors come from complete influence functions with one row for every unit of `units`, including units observed once or never. They account for the estimated imputation, observation, treatment and counterfactual models, the response rates and the cohort shares used. With `compute_inffunc = FALSE` a fit returns point estimates only, and `aggte()` refuses it.

The argument rules are in [Options reference](../reference/options.md#handling-missing-outcomes), the counterfactual regression in [Treatment adjustment](../methods/treatment.md#change-imputation-uses-a-different-regression), and the reported uncertainty in [Results reference](../reference/results.md#uncertainty).

---

**Reading path, step 12 of 22** · Previous: [Choose an observation model](observation-models.md) · Next: [Parallel trends](../methods/parallel-trends.md)

<details>
<summary>All steps of the reading path</summary>

| Step | Page |
| --- | --- |
| 1 | [Get started](../get-started.md) |
| 2 | [Moving from did](../from-did.md) |
| 3 | [Incomplete panels](incomplete-panels.md) |
| 4 | [Identify effects when the observed units change composition](identification.md) |
| 5 | [How to compare precision](precision.md) |
| 6 | [Prepare your panel and unit table](../help/prepare-data.md) |
| 7 | [Treatment adjustment](../methods/treatment.md) |
| 8 | [Choose IPW, OR or DR in an incomplete panel](treatment-choice.md) |
| 9 | [Handling missing outcomes](../methods/missing-outcomes.md) |
| 10 | [Sampling information for the observation model](../methods/sampling-information.md) |
| 11 | [Choose an observation model](observation-models.md) |
| 12 | **Impute actual changes** (this page) |
| 13 | [Parallel trends](../methods/parallel-trends.md) |
| 14 | [Add PT-all restrictions](pt-all-comparison.md) |
| 15 | [GMM weighting and covariance regularization](../methods/gmm.md) |
| 16 | [Robustness](robustness.md) |
| 17 | [Inference](inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |
| 21 | [Cross-pair chaining and control parallel trends](../methods/cross-pair-chaining.md) |
| 22 | [Short pre-treatment windows](short-pretreatment-windows.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
