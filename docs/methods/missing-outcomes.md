# Handling missing outcomes

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

An observed outcome change is available only when the same unit is observed at both dates. Imputation and augmented correction also reconstruct changes on the full study population. The choices on this page decide what must be true for those observed changes to stand in for the study population, and which models adjust for differences between the units observed and the population. The arguments involved begin with `selection`, together with `imputation_pooling` and the unit table `frame`.

The choice is separate from treatment adjustment (`est_method`, see [Treatment adjustment](treatment.md)). Every combination below keeps the same target: the average treatment effects of the original study population.

## The study population and the unit table

The study population is the set of units whose effects the analysis describes. The unit table, supplied as `frame`, lists it: one row per study unit, matched to `data` by `idname`, including units observed once or never. Give it the same column names as `data` for the identifier, the cohort and the covariates, and for the weights and cluster labels when you use them.

Supply positive sampling weights and covariates wherever the chosen route reads them, and cluster labels for every study unit if inference is clustered. Static attributes may be supplied in the unit table; changing values belong in dated `data` rows, with an explicit `fix_weights` reference rule for changing weights. Do not place a conflicting static copy in `frame`. Under `selection = "conditional_mtar"` every row also needs a known cohort. Where a setting allows an unknown cohort for a row present only in the unit table, aggregates that depend on cohort shares are reported as unavailable. Rows that appear only in the unit table keep their cohort, weight and cluster label, so they count in the cohort shares and in the standard errors. No row of the unit table is ever dropped.

When you omit `frame`, the units present in `data` form the study population, including units with a single observed outcome. A unit with no row in `data` is then outside the study population. If part of your population was never observed, supply it through `frame`. The package reads a missing outcome at a date as nonresponse at that date, which is correct only if your data record every observation that exists.

| Kind of unit | What it contributes |
| --- | --- |
| Observed at two or more dates | changes for usable catalog pairs whose endpoints are observed, plus its cohort, required covariates and observation pattern |
| Observed at one date | no observed change; retained population information and every active model/rate influence remain in the calculation across all study units |
| Present only in the unit table | the same as a unit observed at one date |

When inference is requested, the influence function retains every row of the unit table. The reported number of units `n` counts them all. Inactive models are not fitted, and `compute_inffunc = FALSE` returns no standard errors.

## Which changes are available

A change is available for a unit when the unit is observed at both of its dates, and the dates need not be adjacent. This table shows four observation groups over three periods:

| Observation group | Period 1 | Period 2 | Period 3 | Changes it supplies |
| --- | --- | --- | --- | --- |
| Complete observations | observed | observed | observed | 1 to 2, 2 to 3, 1 to 3 |
| First two dates | observed | observed | missing | 1 to 2 |
| Last two dates | missing | observed | observed | 2 to 3 |
| A gap between dates | observed | missing | observed | 1 to 3 |

Under the assumptions below, the average change from period 1 to period 3 can be estimated directly from the units observed at periods 1 and 3, or by adding the estimated changes from 1 to 2 and from 2 to 3, each from its own observed units. Chained DiD combines the usable catalog comparisons by GMM, accounting for shared units and fitted nuisance parameters. Graph connectivity alone is insufficient: control support, valid observation reconstruction and cancellation of unrestricted trends also matter. A requested cell unsupported by the retained catalog is unavailable with a reason; model and numerical failures instead raise explicit errors.

## The assumption about missing outcomes

`selection` states the assumption.

| `selection` | Meaning |
| --- | --- |
| `"none"` (default) | No observation model. Within a cohort, the units observed in a given pair are assumed representative of that cohort with respect to the mean outcome change over that pair. |
| `"conditional_mtar"` | Conditional mean trends missing at random, with an explicit adjustment for observation that depends on covariates. |
| `"corollary1"` | Legacy name for `"conditional_mtar"` with `selection_model = "pair_logit"`. |

### Without a correction

With `selection = "none"`, within a cohort the units observed in a given pair must be representative of that cohort with respect to the mean outcome change over that pair. The estimator scales the observed changes by the cohort's response rate for the pair. The population restriction is

```text
E[ Y(t) - Y(s) | cohort, observed at s and t ] = E[ Y(t) - Y(s) | cohort ]
```

This is equality of population expectations, not a requirement that each unit have the same change or that realized sample means agree exactly. The distributions of changes may differ. Observation may depend on a unit fixed effect removed by differencing, provided that dependence does not select the mean change. If the fixed effect predicts trends or heterogeneous treatment effects, selecting on it can also select changes. Independence is sufficient for the displayed mean restriction but is stronger than required. No observation model is fitted, and a nondefault `selection_model` is refused under this setting.

### Conditional mean trends missing at random

With `selection = "conditional_mtar"`, the package assumes that for every cohort, including the units never treated, and for every pair of dates (s, t) that the estimator actually uses, including the comparisons reported before treatment,

```
E[ Y(t) - Y(s) | X, cohort, observed at s and t ]  =  E[ Y(t) - Y(s) | X, cohort ]
```

and, for observation weighting and augmented correction, that the probability of observing both dates conditional on X and cohort is positive. Pure imputation instead needs observed training units and sufficient model/prediction support in each needed cell. The expectations refer to the study population described by the unit table and its weights.

Several points follow from this statement.

- **Actual changes of treated units are covered.** For a treated cohort, `Y(t) - Y(s)` is the actual change, which includes the treatment effect. An assumption about untreated outcomes alone does not cover these changes.
- **Every pair actually used is covered.** The restriction must hold for each date pair the estimator uses; holding it for adjacent dates alone is not enough.
- **Means, not distributions.** Only the conditional mean change of observed units must equal that of the population; the full distributions may differ.
- **What observation may depend on.** Observation may depend on covariates and a fixed unit component removed by differencing. Any further dependence is admissible only if the displayed conditional-mean equality remains valid; full conditional independence is a stronger sufficient condition.
- **The conditioning set needs justification.** For the augmented sufficient conditions below, W includes the `xformla` covariates and any variables only in `selection_xformla`. Observation weighting likewise needs valid selection and treatment restrictions for the scientific conditioning variables; different predictor lists require separately sufficient restrictions. Supplying formulas does not establish their sufficiency.

These are population conditional expectations, not exact finite-sample equalities or restrictions on each unit's change. Conditional representativeness still needs a scientifically sufficient conditioning set, the relevant models and positivity or prediction support for the selected adjustment. It is a substantive identification assumption. A well fitting observation model, balanced covariates or a plot of observed missingness cannot establish it. [Investigate missingness and data structure](../help/investigate-missingness.md) separates useful observable checks from these untestable conditions.

## How to adjust

With `selection = "conditional_mtar"`, `selection_est_method` chooses the adjustment.

| `selection_est_method` | Name in this guide | Models fitted for missing outcomes |
| --- | --- | --- |
| `"ipw"` (default) | observation weighting | probability of observing both dates |
| `"imputation"` | change imputation | the actual change of each cohort and date pair |
| `"dr"` | augmented correction | both, with `est_method = "dr"` only |

### Observation weighting

Observation weighting weights each observed change by the inverse of its fitted probability of being observed at both dates, combined with the cohort's repeat rate as the sampling adjustment. It is distinct from treatment weighting, which concerns membership in the cohort. `selection_model` chooses the probability model:

- `"pair_logit"` (default) models the probability of observing both dates directly: one logit per cohort and pair of dates for the event that both are observed. Because it models the pair as a whole, it allows observation at the two dates to be dependent, as with dropout.
- `"independent_period_logit"` models observation at each date separately and multiplies the two date probabilities. This requires an additional assumption: given the cohort and the covariates, observation at one date must be independent of observation at the other. The assumption fails, for example, when units that leave the panel do not return. This model requires `selection_universe = "F"`.

Both models are logits. A fit with separation, a boundary solution or deficient rank is refused. Since version 1.1.0.0, inverse observation denominators use the fixed `1e-6` floor described under [Probability safeguards](#probability-safeguards); no units or rows are trimmed.

### Change imputation

Change imputation predicts actual changes and averages the predictions over the study population. For each cohort k and each date pair e that the estimator needs (for the effects and for both `base_period` reporting families), a linear model of the actual change is fitted by weighted least squares to the units of cohort k observed at both dates, with their fixed weights:

```
j(k, e)(X) = alpha(k, e) + L(e) * z(X)' beta(r)
```

Here `L(e)` is the number of steps between the two dates in the stored calendar, `alpha(k, e)` is a local intercept for the cell when `xformla` has an intercept, and the slope vector `beta(r)` belongs to the pooling group r that `imputation_pooling` defines. Each fitted model is then evaluated for every unit of cohort k in the unit table, including the units whose change was observed, and the average of these predictions replaces the observed mean. No observation probability is modeled, so `selection_model`, `selection_pooling` and `selection_xformla` are inactive and never evaluated.

Change imputation requires `selection = "conditional_mtar"`, `selection_universe = "F"`, `selection_rho = NULL`, known cohorts and complete `xformla` covariates on every row of the unit table, and `est_method` equal to `"ipw"`, `"reg"` or `"dr"`. It also runs on complete panels, where an all-success observation logit with a free cell intercept has no finite interior root.

This model is different from the counterfactual outcome regression of treatment OR. It predicts the actual change of each cohort, treatment effects included, while the counterfactual regression predicts untreated changes from controls.

### Augmented correction

The augmented correction combines both models in one estimator. For each needed cell c, a cohort and date pair, every unit of the cohort receives the reconstructed change

```
Z(c) = j(c)(X) + S(e) * { change(e) - j(c)(X) } / max(q(c)(W), 1e-6)
```

where `j(c)` is the change imputation model, `q(c)` is the fitted probability of observing both dates, `S(e)` is 1 for units observed at both dates and 0 otherwise, and W contains the covariates of both formulas. In words, the package imputes each unit's change and then adds back the prediction errors of the observed units, weighted by their inverse observation probabilities. The correction term is not renormalized; version 1.1.0.0 floors its inverse observation denominator at `1e-6`. The result is a single estimator: in finite samples it generally differs from both the weighting and the imputation estimates, and it is not their average.

It requires `est_method = "dr"`; requests with `"ipw"`, `"reg"` or a custom function are refused before the data are read. The declaration requirements of change imputation apply, and complete covariates of both formulas are needed on every row of the unit table.

Each model must be estimable. An active pair cell, or endpoint-date cell for a date model, that is all observed or all unobserved cannot fit a logit with its own free intercept and is refused. Pooling slopes does not remove those cell intercepts. Pure imputation avoids that unused observation-logit obstruction but must still satisfy its own model/support checks. Failures of the imputation, treatment or counterfactual models are refused as well. Nothing falls back to change imputation or to observation weighting.

## Probability safeguards

Version 1.1.0.0 caps package-owned treatment propensities at `1 - 1e-6` and floors inverse observation denominators at `1e-6`, without trimming units, rows or comparisons. For independent date logits the floor applies to the final pair product, not to each endpoint. Positive empirical pair-rate denominators use the same floor. Raw probabilities remain in likelihood scores and information matrices; influence functions differentiate the weights actually used. Rank, separation, convergence, boundary and support refusals still apply.

When a safeguard binds with positive population probability, the clipped functional can differ from the original full-population ATT and can alter the usual robustness branches. The original-target consistency statements below require inactive clipping in population or additional conditions preserving the target; ordinary influence-function inference also needs negligible threshold-crossing remainder. Clipping does not establish identification or nominal coverage. See the probability-clipping section of `?cdid::att_gt` for diagnostics and the exact cutoff convention.

## Which models must be correct

All three adjustments need conditional mean trends missing at random for every needed cell, parallel trends and no anticipation, and the condition of the treatment adjustment: a correct treatment propensity model for `est_method = "ipw"`, a correct shared counterfactual regression for `"reg"`, and either of the two for `"dr"`. Observation weighting and augmented correction also require positive observation probabilities; pure imputation needs observed training mass and sufficient prediction support rather than a fitted probability model. In addition:

| Adjustment | What must also be correct |
| --- | --- |
| Observation weighting | the observation model for every needed cell; with `"independent_period_logit"`, also the independence of the two dates given cohort and covariates |
| Change imputation | the linear mean of the actual change within each pooling group |
| Augmented correction | for each needed cell, either the imputation model or the observation model |

For the augmented correction combined with treatment DR, the sufficient condition is therefore: (treatment propensity model correct or counterfactual regression correct) and, for every needed cell, (observation model correct or imputation model correct). The correct alternative may differ from one cell to another. These are sufficient conditions for the primitive links with several qualifications. ATT/GMM recovery additionally needs valid moments for every used catalog link and the accepted identification and regularity conditions; incompatible link pseudo-targets can make the limit depend on the GMM weights.

- Pooled slopes, the two dates of a date model and the shared counterfactual regression link cells together, so an error in another cell, block or pooling group can move a shared coefficient.
- If variables that appear only in `selection_xformla` shift the mean change, an imputation model built from `xformla` cannot be correct on its own, and the observation model must be. Adding observation covariates is therefore not automatically harmless.
- When both alternatives of a cell are wrong, or both treatment models are wrong, the estimate is generally biased.
- The results are large sample properties. There is no finite sample unbiasedness, and small fitted probabilities inflate the variance.

## Pooling: how much the models share

`selection_pooling` applies to the observation models and `imputation_pooling` to the change imputation models. Both take the same three values.

| Value | Observation models (`selection_pooling`) | Change imputation (`imputation_pooling`) |
| --- | --- | --- |
| `"separate"` (default) | separate slopes for each cell, a date pair or a date | separate slopes for each cohort and date pair |
| `"within_cohort"` | slopes shared within each cohort, across its pairs or dates | slopes shared within each cohort, across its pairs |
| `"across_cohorts"` | one slope vector shared by all cells | one slope vector shared by all cells |
| Intercepts | separate for each cell when `selection_xformla` has an intercept | separate for each cohort and pair when `xformla` has an intercept |

Shared slopes are an explicit model restriction that you choose. They can be wrong, for example when the changes of a treated cohort respond to covariates differently across date pairs. The package never pools automatically, and a pooling group whose design is rank deficient is refused. `selection_pooling` changes only the slopes of the observation models; the shared slope and local intercepts of the counterfactual outcome regression are unaffected. Under the augmented correction the two pooling arguments act independently.

## Two covariate formulas

| | `xformla` | `selection_xformla` |
| --- | --- | --- |
| Used by | the treatment logit, the counterfactual regression and change imputation | the observation probability model only |
| Default | `NULL`, an intercept only | `NULL`, which inherits `xformla` |
| Under change imputation | used | inactive and never evaluated |
| Under the augmented correction | used | defines the observation design; its additional variables never enter the outcome or treatment models |

Both are one sided formulas of covariates that are fixed within units or read at a reference date. Intercept only and no intercept forms are supported. [Options](../reference/options.md#handling-missing-outcomes) lists the remaining details.

## Supported combinations

With `selection = "conditional_mtar"`, seven built-in combinations of treatment adjustment and adjustment for missing outcomes are supported:

| `est_method` | Observation weighting | Change imputation | Augmented correction |
| --- | --- | --- | --- |
| `"ipw"` | supported | supported | refused |
| `"reg"` | supported | supported | refused |
| `"dr"` | supported | supported | supported |

`selection = "none"` works with all three methods and with a custom IPW function. A custom function also works with observation weighting when it returns `log_weight_gradient`; it is refused with change imputation and the augmented correction.

## Advanced: fitting the observation model on part of the unit table

By default the observation model is fitted on the whole unit table (`selection_universe = "F"`). For observation weighting with `selection_model = "pair_logit"`, it can instead be fitted on the units observed at least once (`"B"`) or at least twice (`"A"`). These counts describe outcome observations only; they do not define structural observation windows. The study population and the standard errors still use the whole unit table.

Fitting on a subset requires a stated justification in `selection_justification`. With `selection_rho = NULL`, the justification asserts that, within each cohort, the fitting subset has the same covariate distribution as the cohort's study population. Otherwise `selection_rho` names a column of fixed, known, positive, normalized density ratios `rho = f(X | G, U) / f(X | G)` under the declared weighted population law, with the fitting-subset density in the numerator; the package neither estimates nor renormalizes it. The package records the justification, and no fit diagnostic can confirm it. A correct logit on the whole unit table does not imply a correct logit on a subset; the two are different modeling assumptions with no ranking of efficiency. [Sampling information for the observation model](sampling-information.md) explains these fitting sets and known transport declarations.

## Example: observation that depends on a covariate

In this simulated design, a fixed covariate `x` shifts the untreated trend and the treatment effect, and outcomes with larger `x` are observed more often. Observation at each date is independent of observation at the other dates given `x`.


```r
set.seed(20261026)

n_units <- 2000
n_periods <- 6

x <- rnorm(n_units)
first_treat <- sample(c(0, 3), size = n_units, replace = TRUE)
unit_effect <- rnorm(n_units)

full_panel <- data.frame(
  id = rep(seq_len(n_units), each = n_periods),
  period = rep(seq_len(n_periods), times = n_units),
  first_treat = rep(first_treat, each = n_periods),
  x = rep(x, each = n_periods)
)

treated_now <- full_panel$first_treat > 0 &
  full_panel$period >= full_panel$first_treat
# The effect is larger for units with larger x.
effect <- ifelse(
  treated_now,
  1 + 0.25 * (full_panel$period - full_panel$first_treat) + 0.5 * full_panel$x,
  0
)
full_panel$y <- rep(unit_effect, each = n_periods) +
  (0.5 + 0.3 * full_panel$x) * full_panel$period +
  effect +
  rnorm(nrow(full_panel))

# Each outcome is observed independently across dates given x,
# more often for larger x.
observed <- runif(nrow(full_panel)) < plogis(1.5 + full_panel$x)
obs_panel <- full_panel[observed, ]

# Unit table: every study unit, including units observed once or never.
units <- data.frame(id = seq_len(n_units), first_treat = first_treat, x = x)

# Number of units by number of observed outcomes (0 to 6).
table(table(factor(obs_panel$id, levels = units$id)))
```

Output:

```text
0   1   2   3   4   5   6 
 13  56 109 213 393 540 676
```

The four fits below hold the data, the unit table, the covariate formula and the treatment adjustment (`est_method = "dr"`) fixed and change only the handling of missing outcomes.

```r
fit_none <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = obs_panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  selection = "none",
  bstrap = FALSE, cband = FALSE
)

fit_weighting <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = obs_panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  selection = "conditional_mtar",
  selection_est_method = "ipw",
  selection_model = "independent_period_logit",
  selection_pooling = "separate",
  bstrap = FALSE, cband = FALSE
)

fit_imputation <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = obs_panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  selection = "conditional_mtar",
  selection_est_method = "imputation",
  imputation_pooling = "separate",
  bstrap = FALSE, cband = FALSE
)

fit_augmented <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = obs_panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  selection = "conditional_mtar",
  selection_est_method = "dr",
  selection_model = "independent_period_logit",
  selection_pooling = "separate",
  imputation_pooling = "separate",
  bstrap = FALSE, cband = FALSE
)

fits <- list(none = fit_none, weighting = fit_weighting,
             imputation = fit_imputation, augmented = fit_augmented)
simple <- lapply(fits, cdid::aggte, type = "simple")
data.frame(
  handling = names(simple),
  estimate = vapply(simple, function(a) a$overall.att, numeric(1)),
  std.error = vapply(simple, function(a) a$overall.se, numeric(1))
)

# Value of the simple aggregate for these 2000 simulated units.
1.375 + 0.5 * mean(units$x[units$first_treat == 3])

fit_augmented$chained$selection[c("method", "selection_model", "selection_pooling",
                                  "imputation_pooling", "selection_formula_inherited")]
```

Output:

```text
handling estimate  std.error
none             none 1.545866 0.06311893
weighting   weighting 1.272803 0.06791702
imputation imputation 1.314577 0.06505744
augmented   augmented 1.271517 0.06743981
[1] 1.344878
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

**What the design implies.** In the population, the simple aggregate equals 1.375, because `x` has mean zero in cohort 3; the last calculation gives its value for these 2000 units. Both treatment models are correct here, since treatment is assigned at random and untreated changes are linear in `L(e) * x`. The four fits differ only in their handling of missing outcomes.

- `fit_none`: the assumption of `selection = "none"` fails, because observation depends on `x`, which also shifts the treated changes. Its estimates are generally biased for the population effects.
- `fit_weighting`: the date model is correctly specified, because each date's observation is a logit in `x` and dates are independent given `x`. A pair logit linear in `x` is generally misspecified here: the product of the two date logits is generally not itself logit-linear in `x`.
- `fit_imputation`: with `"separate"`, each cohort and pair has its own slope, which this design needs: for cohort 3, the change over a pair that spans the start of treatment depends on `x` differently from other pairs. Slopes shared within the cohort, `"within_cohort"` or `"across_cohorts"`, would be misspecified.
- `fit_augmented`: both alternatives are correct in every cell.

The last three fits are therefore expected to target the effects, with different finite sample estimates and standard errors.

## Where the result records the choice

| Location | Contents |
| --- | --- |
| `fit$chained$selection` | the adjustment: `method` is `"none"`, `"corollary1"` for observation weighting (its legacy record name, with `canonical_assumption = "conditional_mtar"`), `"imputation"` or `"augmented"`; for the corrections, also the observation model, the pooling choices, the fitting set, the transport, the justification, the assumptions and counts of fitted models; for the augmented correction, `selection_formula_inherited` reports whether `selection_xformla` was inherited |
| `fit$chained$nuisance_diagnostics` | the fitted cells, the pooling groups with their ranks, slopes and removed columns, the observation probability blocks (`records$selection$blocks` for the augmented correction), the counterfactual regression, the treatment logits and the support of each cell |
| `fit$chained$outcome_design` | the outcome design, and for the augmented correction both designs |
| `summary(fit)` | an `Observation correction` line naming the adjustment and its model and pooling choices |

---

**Reading path, step 9 of 20** · Previous: [Choose IPW, OR or DR in an incomplete panel](../examples/treatment-choice.md) · Next: [Sampling information for the observation model](sampling-information.md)

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
| 7 | [Treatment adjustment](treatment.md) |
| 8 | [Choose IPW, OR or DR in an incomplete panel](../examples/treatment-choice.md) |
| 9 | **Handling missing outcomes** (this page) |
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
