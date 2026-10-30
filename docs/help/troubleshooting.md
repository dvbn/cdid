# Diagnose errors and unavailable results

A chained fit can stop with an error, report effects as unavailable, or complete with surprising results. Each calls for a different check, so first tell them apart.

| What you see | What it usually means | First check |
| --- | --- | --- |
| An error before estimation starts | A refused request or input | The condition class |
| A cell or aggregate reported as `NA` with a reason | The available information cannot support that target | Cell `support.status`/`support.reason`, or aggregate `reason` |
| An error during the fit that names a model, covariance or solve | A computation failed on this sample | The named component and the data that feed it |
| A completed fit with unexpected estimates or standard errors | Assumptions, overlap, sample size or estimator choices | The observation pattern and the inference record |

The package never drops units or comparisons, pools models or substitutes another method to make a fit succeed. Version 1.1.0.0 applies fixed probability safeguards to finite weights, but these do not repair failed nuisance fits or support violations ([Probability safeguards](../methods/missing-outcomes.md#probability-safeguards)). A change you make in response, such as a different formula, pooling choice, comparison group or parallel trends assumption, changes the analysis and its assumptions. Justify it from the application, not from whether the new fit runs.

The inspection snippets use your own data and fitted result, with the names in [Prepare your panel and unit table](prepare-data.md): a long panel `panel_data` with columns `id`, `time`, `g` and `y`, the unit table `units`, and a completed chained result `fit`.

[Investigate missingness and data structure](investigate-missingness.md) provides a self-contained example of study-population coverage, selected-pair counts and the fitted diagnostics. Start there when the main question is what the observation pattern supplies.

## Read the condition class

Typed chained refusals and fit failures inherit class `cdid_error`. Keep the condition and read it instead of retrying with other settings:

```r
res <- tryCatch(
  cdid::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
               data = panel_data, frame = units, bstrap = FALSE, cband = FALSE),
  cdid_error = function(e) e
)
if (inherits(res, "cdid_error")) {
  print(class(res))
  cat(conditionMessage(res), "\n")
}
```

| Class | Cause | Where to look |
| --- | --- | --- |
| `cdid_mode_error`, `cdid_chained_panel` | An invalid `mode` or `panel`, or chained mode without `panel = TRUE` | [Repeated cross sections](../from-did.md#repeated-cross-sections) |
| `cdid_chained_balanced` | `allow_unbalanced_panel = FALSE` in chained mode | `did`'s balanced subset needs `mode = "ordinary"` |
| `cdid_chained_unavailable` | An unknown `est_method` string, or `selection_est_method = "dr"` without `est_method = "dr"` | [Supported combinations](../reference/options.md#supported-combinations) |
| `cdid_chained_request_*` | Other unsupported requests or inputs, such as an invalid `pt_start` | The argument named in the message, in the [options reference](../reference/options.md) |
| `cdid_chained_request_parallel_error` | `pl = TRUE, cores = k` beyond `control$workers` | [Resources and parallel draws](#resources-and-parallel-draws) |
| `cdid_chained_public_control_error`, `cdid_chained_public_accuracy_deferred` | An invalid `control` list, or the entry `relative_accuracy`, which is not provided | [The control list](../reference/options.md#the-control-list) |
| `cdid_custom_ipw_error` | Malformed or nonfinite output of a custom IPW function | [Treatment adjustment](../reference/options.md#treatment-adjustment) |
| `cdid_chained_request_no_inffunc_error` | `aggte()` on a fit with `compute_inffunc = FALSE` | That fit has point estimates only |
| `cdid_untested_upstream` | `did` older than 2.5.1 or `DRDID` older than 1.3.0 | The library your session loads |

An upstream-version refusal can occur while loading the package, before the call inside this handler starts; inspect the loaded library versions separately. R's own dependency check can refuse an older `did` first, with its own message instead of this class. Later versions of `did` and `DRDID` load but have not been tested with this release.

Repeated unit and date pairs, dates that are not finite numbers and cohorts that change within a unit are easiest to find before fitting, with the [base R checks](prepare-data.md#concise-base-r-data-checklist-and-examples). No count of observed outcomes reveals missing outcomes recoded as 0, so check how the data were built.

## Inspect a completed fit

Several chained records are nested lists whose layout depends on the estimator. Rely only on entries named in the [results reference](../reference/results.md#the-chained-record) or the methods pages, and read the others with `str()` first.

```r
tab <- cdid::tidy(fit)
tab[tab$available %in% FALSE, c("target.id", "role", "support.status", "support.reason")]
cdid::glance(fit)[, c("n.requested", "n.available", "inference", "clustered")]

fit$chained$inference$kind     # "analytic", "multiplier_bootstrap_iqr_se" or "not_computed"
fit$chained$pretest$status     # status of the unconditional pretrend test
fit$chained$selection$method   # handling of missing outcomes
str(fit$chained$covariance_handling, max.level = 1)
str(fit$chained$outcome_design, max.level = 1)
str(fit$chained$nuisance_diagnostics, max.level = 1)
```

## An effect is unavailable

An unavailable cell has `NA`, never zero, as its estimate, standard error and influence function column, with a support status and reason. Either no chain of observed changes connects the needed dates for the cohort or its controls, or the parallel trends choice does not cover a needed step. The status is relative to the comparisons available to the estimator; it is not a finding that the effect is unidentified in the population.

```r
# Choose a few relevant pairs in your stored calendar, rather than expanding
# all within-unit pairs. The unit table supplies the full denominators.
obs <- panel_data[!is.na(panel_data$y), c("id", "time", "g")]
selected_pairs <- data.frame(start = c(1, 2), end = c(2, 3))
pair_support <- do.call(rbind, lapply(seq_len(nrow(selected_pairs)), function(j) {
  s <- selected_pairs$start[j]
  t <- selected_pairs$end[j]
  both <- intersect(obs$id[obs$time == s], obs$id[obs$time == t])
  do.call(rbind, lapply(sort(unique(units$g)), function(g) {
    take <- units$g == g
    data.frame(start = s, end = t, cohort = g, n_roster = sum(take),
               n_both = sum(units$id[take] %in% both),
               share_both = mean(units$id[take] %in% both))
  }))
}))
pair_support

sort(unique(panel_data$time))  # a date with no rows at all is never reported
```

| Finding | Implication |
| --- | --- |
| No units of the cohort, or no eligible controls, observed at both dates of a needed pair | That pair supplies no observed comparison. Check for valid alternative paths; a route through earlier dates works only when the unrestricted trends cancel ([Comparisons through earlier dates](../methods/parallel-trends.md#comparisons-through-earlier-dates)). |
| A recovery path needs an additional trend restriction before the reference period | PT-post does not supply it. PT-all may cover the represented step from a `pt_start` chosen before seeing results, if that restriction is credible. False restrictions can bias estimates. |
| A date of the study calendar has no rows | Keep its rows with `y = NA` so that the date and its cells are reported. |

An aggregate that needs an unavailable cell is `NA`; `cdid::tidy(agg)` shows `available`, `n.intended`, `n.used`, `excluded.ids` and `reason`. `na.rm = TRUE` instead reports a restricted target over the available cells, marked in `restricted`, which is a different estimand.

## A logit model has no finite solution

Treatment propensity logits are fitted for `est_method = "ipw"` and `"dr"`; observation logits only under `selection = "conditional_mtar"` with `selection_est_method = "ipw"` or `"dr"`. Logits do not drop aliased columns. A logit fit is refused, with its component named, in these cases:

| Cause | Check |
| --- | --- |
| A covariate or factor level perfectly separates treated from comparison units in a risk set | The covariate by cohort on the actual fitting support: observed pairs for baseline/observation weighting, all study units after reconstructing missing changes for imputation/augmentation |
| A covariate perfectly separates observed from unobserved units in an observation cell | The covariate by pair or endpoint-date observation status within the cohort, according to `selection_model` |
| An active observation cell is all observed or all unobserved and has its own free intercept | Pair observation by cohort for `"pair_logit"`; observation at each endpoint for `"independent_period_logit"` |
| The design is rank deficient | Duplicated or constant columns in the formula |

```r
# Units of the unit table observed at both dates 2 and 3, by cohort
ids_both <- intersect(obs$id[obs$time == 2], obs$id[obs$time == 3])
table(cohort = units$g,
      observed_both = factor(units$id %in% ids_both, levels = c(FALSE, TRUE)))

# Optional: if region is a static factor covariate supplied in the unit table.
# Use the actual fitting support of your route when checking separation.
if ("region" %in% names(units)) {
  with(units[units$id %in% ids_both, ], table(region, cohort = g))
}
```

A zero in either column of the first table marks an all observed or all unobserved cell, and pooling slopes leaves its intercept in place. The fixed `1e-6` inverse-probability floor does not repair this refusal. No observations are trimmed, observing every unit of a cell does not make its probability a known 1, and known probabilities are not supported. The separation check suffices to show that no finite solution exists, but it is not a complete test for separation. A refusal is not a conclusion about identification: it can reflect a lack of overlap or positivity in the population, or only a small sample.

| Observation model setting | Requirement |
| --- | --- |
| `selection_universe = "F"` (default) | Fitted on the whole unit table. |
| `"B"` or `"A"` | Only with observation weighting and `selection_model = "pair_logit"`, plus a `selection_justification`. With `selection_rho = NULL` it asserts that the fitting set has the cohort's covariate distribution; otherwise `selection_rho` names known normalized density ratios. No diagnostic can confirm it, and the population and standard errors still use the whole unit table. |
| `selection_model = "independent_period_logit"` | `"F"`, and observation at the two dates independent given cohort and covariates, which generally fails under permanent dropout. |

Change imputation fits no observation logit, but replaces it with a linear model of actual changes that must be correct within each pooling group. The augmented correction refuses an all observed or all unobserved active observation cell with its own free intercept. Pure imputation avoids that observation-logit obstruction, but must still satisfy its own training, model and prediction-support checks. Pooling and formula changes restrict or change the model. Each is a different assumption, to be judged on its own terms ([Which models must be correct](../methods/missing-outcomes.md#which-models-must-be-correct)).

## An outcome regression or imputation design is refused

Outcome regression (OR), used by `est_method = "reg"` and `"dr"`, and change imputation fit linear models. Exact collinearity is preprocessed; unsupported predictions are refused.

| Situation | What the package does |
| --- | --- |
| Exactly collinear columns on the fitting rows | Pivoted QR (tolerance `1e-7`) proposes an alias, accepted only if the retained columns reconstruct every needed prediction row to arithmetic roundoff, `64 * .Machine$double.eps` times the original column count. The removal is recorded. |
| A direction that varies on treated prediction rows but not among control fitting rows, such as a level present only among treated units | Refused; its slope is never set to zero. |
| A comparison block that is rank deficient on its own | Supported when the pooled control design identifies the shared slopes. |
| Nearly collinear columns that fail the alias check | Retained. The fit is refused if the reduced design is rank deficient or badly conditioned overall; no further columns, observations or moments are dropped. |
| An imputation pooling group whose reduced design remains rank deficient | Refused; nothing is pooled automatically. |

`fit$chained$outcome_design` and `fit$chained$nuisance_diagnostics` record the designs and their preprocessing. The layout depends on the route: inspect the relevant pooled-control or imputation-group record with `str()` rather than assuming one common `basis` subfield. A removed alias reparameterizes redundant columns; it is not an assumption that a coefficient is zero. For a refused fit, tabulate the covariate by cohort as in the logit check.

## Missing covariates or weights versus missing outcomes

| Missing value | How the package reads it |
| --- | --- |
| An outcome | `y = NA` is nonresponse at that date, and the row stays. Never recode it as 0 or carry an earlier value forward. |
| A covariate that changes within units | Read at each cohort's reference date, or at the earlier date of a comparison before treatment, from the actual rows, including rows with `y = NA`. A missing required value is refused with its unit, date, variable and purpose, so keep that row with its covariates. Nothing is carried forward, and a unit table value cannot override changing values in `data`. |
| A sampling weight that changes within units | `fix_weights = NULL` and `"varying"` are refused; choose `"first_period"` or `"base_period"` from the design of the weights. |
| A treated unit's weight at the first stored date, in aggregates of such fits | An aggregate that depends on cohort shares is `NA`, with the unit identifiers. One-cohort bins and a one-bin group overall summary cancel these shares exactly. Units never treated do not need the mass. |

Which rows need covariates depends on the route ([Static versus changing covariates and sampling weights](prepare-data.md#static-versus-changing-covariates-and-sampling-weights)). Values are needed only where an estimate reads them: no route requires every unit to have an observed outcome at the reference date, rows at dates that no comparison reads, or a unit table that is a census. Under `selection = "none"`, rows present only in the unit table may have an unknown cohort, which makes aggregates that depend on cohort shares unavailable. Observation weighting, change imputation and the augmented correction need known cohorts on every row; the last two also need complete covariates there.

## Covariance regularization

Under `gmm_weighting = "optimal"`, the default, the generalized method of moments (GMM) weights comparisons by the inverse of their estimated covariance. When that covariance needs regularization after exact reductions, `gmm_covariance = "tikhonov"` adds a ridge proportional to its trace. At exactly zero covariance a separate least-squares rule applies, with the actual zero sandwich variance. The control entries `condition_trigger` (default `1e-10`) and `ridge_relative` (default `1e-4`) set the trigger and size, and `fit$chained$covariance_handling` records its use. The ridge changes only the weighting: it supplies no identifying information, cannot make an unavailable cell available, and says nothing by itself about identification. Standard errors still use the sandwich variance of the complete influence functions, identity weighting is never substituted, and nonfinite or materially indefinite covariance inputs are refused without clipping eigenvalues or enlarging the ridge.

## Pretrend tests are unavailable

| Test | Unavailable when | Record |
| --- | --- | --- |
| Unconditional Wald test, `fit$W` and `fit$Wpval` | A needed nonreference preperiod cell is unavailable or the full joint covariance fails the Cholesky/solve and conditioning checks; with no such cell it is not applicable, and with `compute_inffunc = FALSE` it is not computed | `fit$chained$pretest$status`: `"computed"`, `"unavailable"`, `"not_applicable"` or `"not_computed"`, with its reason |
| Conditional test, `cdid::conditional_did_pretest()` | A needed cell is unavailable; a covariate needed on the grid is missing or has a level absent from the fitted design; no cell before treatment exists; clustered inference has fewer than two clusters; or the estimated process is exactly zero | `pre$chained$status` and `pre$chained$reason` |

No test on an available subset is substituted. The conditional test is `"not_applicable"` for a custom treatment function. The conditional grid has one point for every unit of the unit table, so one row with a missing `xformla` covariate can make the test unavailable while every effect is estimated.

For the next template, `x1` is an application covariate supplied on the route's required supports and on the full testing grid; replace the formula with your justified design.

```r
set.seed(1)  # an available conditional test uses multiplier draws
pre <- cdid::conditional_did_pretest(
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x1, data = panel_data, frame = units
)
pre$chained$status
pre$chained$reason
```

An unavailable test is not evidence about parallel trends, and a test that does not reject does not establish them. Under PT-all, cells before `pt_start` can reject without contradicting the assumption ([Pretrend tests under each choice](../methods/parallel-trends.md#pretrend-tests-under-each-choice)).

## Standard errors look unusual

Uncertainty rests on large sample approximations: enough independent units or clusters, overlap, the model conditions of the chosen method and regularity conditions. A completed fit does not establish them, and no finite sample, few cluster or nominal coverage guarantee is made ([What the uncertainty assumes](../reference/results.md#what-the-uncertainty-assumes)).

| Symptom | Likely source | Check |
| --- | --- | --- |
| Large standard errors for some cells | Few units observed at both dates of the pairs that reach them | The pair table above |
| Large standard errors under treatment or observation weighting | Weak overlap: fitted probabilities near 0 or 1 | `str(fit$chained$nuisance_diagnostics, max.level = 2)` |
| Clustered standard errors with few clusters | The approximation needs many independent clusters | The clustering record in `str(fit$chained$inference)` |
| Multiplier results that change between runs | New random signs | `set.seed()` before the call; workers do not change the draws |

An unusual standard error alone is not evidence of a numerical error. Investigate as a possible software problem an unexplained difference where chained and `did` point estimates are expected to agree: `est_method = "ipw"`, `selection = "none"`, a complete panel, never treated controls, no covariates, each unit's weight fixed across comparisons, PT-post, and the same `anticipation` and `base_period`.

## Results differ from did

| Setting | Relation to `did` |
| --- | --- |
| The matched case above | Point estimates agree up to floating point rounding. |
| OR or doubly robust (DR) treatment adjustment | Can differ when fitted coefficient restrictions or procedures differ: chained OR shares one slope vector across comparison blocks. The method label alone does not imply a difference. |
| Units not yet treated as controls, even in a balanced panel | Can differ: chained estimation keeps additional valid control pool moments for the same date pair. |
| PT-all | Need not agree; `did` has no such option. |
| Incomplete panel | Different observations and assumptions: `did`'s balanced subset, or its repeated cross section estimators with `allow_unbalanced_panel = TRUE`. |

Match every setting, and compare point estimates, standard errors and aggregates separately. Outside the matched case, a difference can reflect these choices and is not by itself evidence of an error ([Why can balanced results differ](../from-did.md#why-can-balanced-results-differ-and-which-should-i-use)). [Identify effects when the observed units change composition](../examples/identification.md) shows an incomplete panel in which the two routes target different quantities.

## Resources and parallel draws

| `control` entry | Default | Meaning |
| --- | --- | --- |
| `numeric_workspace_bytes` | 1073741824 (1 GiB) | Admission limit for the counted numeric workspace; work is split into chunks of `chunk_rows` rows, chosen automatically |
| `counted_work_limit` | `Inf` | No estimated work cap; a finite positive whole number opts in |
| `workers` | 1 | Worker processes for multiplier draws; never raised automatically |

A fit whose conservative payload or work estimate exceeds these limits is refused before the stage runs, and nothing is dropped or relaxed. The limits do not bound the memory of the R process or CPU time, and the package does not detect the memory or speed of your machine. A larger `numeric_workspace_bytes` admits a larger counted payload; whether the machine can hold it, plus memory the count omits, is yours to judge. A finite `counted_work_limit` should be a deliberate budget: the estimate is conservative, so a tight value can refuse fits your machine could complete. `fit$chained$public$counted_work_limit` and `fit$chained$public$work_cap` record the public cap and the internal ceiling.

Parallel multiplier draws need a matching declaration of workers:

```r
fit <- cdid::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
                    data = panel_data, frame = units,
                    pl = TRUE, cores = 4, control = list(workers = 4))
```

Estimation stays serial, results equal those of `cores = 1`, worker buffers count against the workspace, and later `aggte()` calls reuse the workers. `compute_inffunc = FALSE` returns the same point estimates without influence functions, standard errors, bands or pretrend test, and `aggte()` refuses such a fit.

---

**Reading path, step 20 of 20** · Previous: [Results](../reference/results.md) · Last step; back to the [guide home](../README.md)

<details>
<summary>All steps of the reading path</summary>

| Step | Page |
| --- | --- |
| 1 | [Get started](../get-started.md) |
| 2 | [Moving from did](../from-did.md) |
| 3 | [Incomplete panels](../examples/incomplete-panels.md) |
| 4 | [Identify effects when the observed units change composition](../examples/identification.md) |
| 5 | [How to compare precision](../examples/precision.md) |
| 6 | [Prepare your panel and unit table](prepare-data.md) |
| 7 | [Treatment adjustment](../methods/treatment.md) |
| 8 | [Choose IPW, OR or DR in an incomplete panel](../examples/treatment-choice.md) |
| 9 | [Handling missing outcomes](../methods/missing-outcomes.md) |
| 10 | [Sampling information for the observation model](../methods/sampling-information.md) |
| 11 | [Choose an observation model](../examples/observation-models.md) |
| 12 | [Impute actual changes](../examples/change-imputation.md) |
| 13 | [Parallel trends](../methods/parallel-trends.md) |
| 14 | [Add PT-all restrictions](../examples/pt-all-comparison.md) |
| 15 | [GMM weighting and covariance regularization](../methods/gmm.md) |
| 16 | [Robustness](../examples/robustness.md) |
| 17 | [Inference](../examples/inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | **Diagnose errors and unavailable results** (this page) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
