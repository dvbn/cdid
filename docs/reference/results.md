# Results reference

This page describes what the fitted objects contain, how to read their tables and plots, which uncertainty each output reports, and how to export results. It covers only functions and fields that exist in the current source. Additional tables that have been proposed but not implemented are listed at the end as future work.

## Result classes and functions

| Function | Result in chained mode | Methods |
| --- | --- | --- |
| `att_gt()` | `cdid_MP`, which inherits `did`'s `MP` | `print()`, `summary()`, `tidy()`, `glance()`, `ggdid()`, and `aggte()` for aggregation |
| `aggte()` | `cdid_AGGTEobj`, which inherits `did`'s `AGGTEobj` | `print()`, `summary()`, `tidy()`, `glance()`, `ggdid()` (except for the simple type) |
| `conditional_did_pretest()` | `cdid_MP.TEST`, which inherits `did`'s `MP.TEST` | `print()`, `summary()` |

In ordinary mode each function returns `did`'s own object, unchanged. For each chained class, `print()` and `summary()` give the same display. The `att_gt()` and `aggte()` displays show `did`'s summary table under a heading that names the chained estimator, without `did`'s citation line, followed by a short chained record. For `att_gt()` the record gives the available and requested cells with the reasons for unavailable ones, the units of the unit table, the correction for missing outcomes, the parallel trends assumption and `base_period`, the GMM weighting and any covariance regularization, the inference actually performed and the pretrend test or the reason it is unavailable. For `aggte()` it gives the available and intended aggregates, every restricted or unavailable target with its reason, the bins used and excluded by the overall summary, the parallel trends assumption and the inference, and, when they apply, nondefault event-time selectors, unknown cohort shares or shares taken as masses at the first stored date. The heading shows the package versions and both paper references once. Interval family codes, the numerical policy and resource profile remain in the result without being printed, as listed below. The conditional-test display reports its CvM test and status.

## Group time results

### Fields inherited from did

| Field | Contents in a chained result |
| --- | --- |
| `group`, `t` | cohort and period of each reported cell |
| `att` | estimates; `NA` for unavailable cells, never zero filled; exactly zero for the reference cell under `base_period = "universal"` |
| `se` | standard errors, analytical or from multiplier draws according to the inference actually performed; `NA` for unavailable cells |
| `c` | critical value: the pointwise normal value `qnorm(1 - alp / 2)` or a simultaneous multiplier quantile |
| `V_analytical` | estimated asymptotic covariance of `sqrt(n) * (estimate - target)`, equal to `n` times the estimated covariance of the cell estimates; analytical standard errors are `sqrt(diag(V_analytical) / n)` |
| `inffunc` | influence function matrix with one row per unit of the unit table (row names are the unit identifiers) and one column per cell; columns of unavailable cells are `NA` |
| `n` | number of units in the unit table |
| `W`, `Wpval` | unconditional pretrend Wald statistic and p value; `NULL` when the test is not computed |
| `alp`, `DIDparams` | significance level and settings |

### The chained record

`fit$chained` holds everything specific to chained estimation.

For a worked inspection of unit table coverage, missingness, selected-pair counts and these fitted records, see [Investigate missingness and data structure](../help/investigate-missingness.md).

| Element | Contents |
| --- | --- |
| `estimator`, `references`, `versions` | estimator label, the two references, and the versions of `cdid` and `did` with the tested `did` version |
| `target_id`, `target_meta`, `role`, `available` | identifiers and descriptions of the reported cells; `role` is `"pre"`, `"post"` or `"reference"`; `available` is `TRUE` or `FALSE` |
| `support_status`, `support_reason`, `support`, `support_ledger` | why each cell is or is not available, with the supporting records |
| `deterministic_reference` | cells that are exact zeros by construction, which have no test statistic |
| `z_stat`, `p_value`, `analytical_se`, `ci_lower`, `ci_upper`, `cband_lower`, `cband_upper` | statistics, analytical standard errors, pointwise intervals and band limits for each cell, named by `target_id` |
| `band_family`, `band_families` | the interval family of the reported cells; values beginning with `"simultaneous"` denote a simultaneous band, values beginning with `"pointwise"` pointwise intervals |
| `inference` | the inference actually performed: `kind` (`"analytic"`, `"multiplier_bootstrap_iqr_se"` or `"not_computed"`), the effective `bstrap` and `cband`, the number of draws, clustering and number of clusters, `n_total`, the random number record and the requested flags |
| `pretest` | status and details of the unconditional pretrend test |
| `frame` | the unit table as used: `id`, `cohort`, `weight`, `cluster` and `eligible` |
| `base_period`, `pt` | the reporting choice and the parallel trends record |
| `covariance_handling` | the GMM weighting request and the covariance regularization record |
| `selection` | the handling of missing outcomes and model descriptions. Fit counts are nested in `selection$models` and depend on the active route. Observation weighting also records its fitting universe and active cells. |
| `nuisance_diagnostics`, `outcome_design` | the fitted treatment, outcome, observation and imputation models, and the design columns |
| `numerical_certificate` | a numerical certificate when one is produced; the current numerical policy uses conventional computation without per fit error certification |
| `request`, `public` | the honored computation flags with the admitted multiplier workers, and the resource profile with the numerical policy |
| `basis` | optional separate contrasts, kept apart from the cells |

### Settings in DIDparams

`fit$DIDparams` keeps `did`'s settings as requested, including `est_method`, `control_group`, `anticipation`, `base_period`, `xformla`, `bstrap`, `cband`, `biters`, `clustervars` and `compute_inffunc`. `fit$DIDparams$cdid` adds the mode, the estimator label, the versions, the requested and effective inference flags and the number of draws actually made. It also records `data_retained = FALSE`: a chained result does not store the panel data, so it must be aggregated with `cdid::aggte()`.

## Tables from tidy and glance

`tidy()` of a group time result returns one row per reported cell, with `did`'s columns (including `estimate`, `std.error`, `conf.low` and `conf.high`) and the following chained columns or updates:

| Column | Meaning |
| --- | --- |
| `statistic`, `p.value` | test of a zero effect for the cell; `NA` for exact zero reference cells |
| `target.id`, `role` | cell identifier and role (`"pre"`, `"post"` or `"reference"`) |
| `available`, `support.status`, `support.reason` | availability and its reason |
| `certified` | whether a numerical certificate covers the cell; `NA` when no certificate was produced, the usual case |
| `estimator` | estimator label, present for OR and DR fits and for change imputation and the augmented correction |
| `band.family` | interval family |

The `statistic` and `p.value` columns update the inherited values. P-values remain pointwise even when `conf.low` and `conf.high` use a simultaneous critical value; the inherited `point.conf.low` and `point.conf.high` are always pointwise normal intervals.

`glance()` of a group time result returns one row with `did`'s columns plus `estimator`, `n.requested`, `n.available`, `band.family`, `inference` and `clustered`.

`tidy()` of an aggregated result updates the inherited `statistic` and `p.value` and adds `label`, `available`, `restricted`, `n.intended`, `n.used`, `excluded.ids`, `reason`, `certified`, `band.family` and, as above, `estimator`. Its rows depend on the type: the group type starts with the overall summary (`label` equal to `"overall"`) followed by one row per cohort; the simple type has one row; the dynamic and calendar types have one row per event time or period and no overall row. `glance()` of an aggregated result adds `estimator`, `status`, `na.rm`, `n.bins`, `n.bins.available`, `overall.status`, `band.family` and `inference`.

## Unavailable effects

A cell is unavailable when the comparisons available to the estimator cannot support it, for example because no chain of observed changes connects the needed dates for the cohort or its controls, or because the parallel trends choice does not cover a needed step. Its estimate, standard error, influence function column and covariance entries are `NA`, never zero, and `support_status` and `support_reason` give the reason. The summary lists unavailable cells. An unavailable cell reflects the information in the available comparisons; whether the effect is identified in the population is a separate question.

Failures of a model fit, a covariance calculation or a numerical solve are different: the fit stops with an error that names the problem, and nothing is dropped or relaxed to obtain a result.

Unavailable cells also affect other outputs. A simultaneous band covers only the available cells. The pretrend test is unavailable when any of its cells is. An aggregate that needs an unavailable cell is `NA` with the missing cells named, unless `na.rm = TRUE` requests the restricted target over the available cells, which is then marked in `restricted` and `excluded.ids`.

## Aggregated results

### Four aggregation types

`aggte()` follows `did`'s definitions. Cohort shares are estimated from the unit table, and the standard errors include their uncertainty.

| Type | Separate estimates (`egt`, `att.egt`) | Overall summary (`overall.att`) |
| --- | --- | --- |
| `"simple"` | none | average of the intended effects after treatment, each weighted by its cohort's share; missing cells require explicit `na.rm = TRUE` for a restricted target |
| `"group"` | one per cohort: the average of the cohort's effects after treatment | average of the cohort estimates, weighted by cohort shares |
| `"dynamic"` | one per event time `e = t - g`, including event times before treatment | average of the estimates for event times from 0 onward |
| `"calendar"` | one per calendar period: the share weighted average over the cohorts treated by then | average of the calendar period estimates |

With several cohorts the overall summaries answer different questions. The simple type gives more weight to cohorts observed for more periods after treatment, the group type weights each cohort by its share whatever the number of its periods, and the dynamic and calendar types weight event times or periods equally. With one treated cohort, as in [Get started](../get-started.md#aggregate-the-effects), they coincide in the population. For fits with changing covariates or weights, the cohort shares are masses at the first stored date, as described in [Changing covariates and weights](options.md#changing-covariates-and-weights).

Share uncertainty enters only when the target depends on shares. For one-cohort bins and a one-bin group overall summary the ratio cancels shares exactly, so missing first-date masses do not by themselves make those targets unavailable.

### Fields of an aggregated result

The fields inherited from `did` are `overall.att`, `overall.se`, `type`, `egt`, `att.egt`, `se.egt`, `crit.val.egt`, `inf.function`, `call` and `DIDparams`, and for the dynamic type `min_e`, `max_e` and `balance_e`. For the simple type the separate estimates are `NULL` and `overall.att` holds the estimate. `inf.function` holds complete influence functions with one row per unit of the unit table, named as in `did`: `simple.att`; `selective.inf.func.g` and `selective.inf.func`; `dynamic.inf.func.e` and `dynamic.inf.func`; `calendar.inf.func.t` and `calendar.inf.func`.

`agg$chained` adds:

| Element | Contents |
| --- | --- |
| `type`, `status`, `na_rm`, `reason`, `target_definition` | what was aggregated and whether it is available |
| `cell` | one row per aggregate bin, including the scalar simple target: `label`, `key`, `available`, `restricted`, `share_dependent`, `n_intended`, `n_used`, `intended_ids`, `used_ids`, `excluded_ids`, `reason` |
| `cell_target_ids`, `cell_used_target_ids` | the cells each estimate intended to use and actually used |
| `overall` | the overall summary: `status`, `reason`, `intended_bins`, `used_bins`, `excluded_bins`, `weights`, `inference_kind`, `n_draws`, `ci_lower`, `ci_upper` and related records |
| `shares_known`, `cohort_shares`, `weights_theta` | the cohort shares and the weights of the cells |
| `z_stat`, `p_value`, `ci_lower`, `ci_upper`, `cband_lower`, `cband_upper` | statistics and limits of the separate estimates |
| `se_fixed_shares_only`, `if_fixed_full`, `if_share_full`, `decomposition` | standard errors treating the shares as fixed, kept as a diagnostic, and the parts of the influence function |
| `band_family`, `inference`, `requested` | the inference performed and the effective request, with whether each setting was inherited or explicit |
| `selectors` | how `min_e`, `max_e` and `balance_e` were applied |
| `share_law`, `cohort_mass_complete` | present for fits whose shares are first date masses |

## Uncertainty

### Which interval each output reports

The defaults are `bstrap = TRUE`, `cband = TRUE`, `biters = 1000` and `alp = 0.05`. The behavior depends on the output, following `did` 2.5.1:

| `bstrap` | `cband` | Cells of `att_gt()` | `aggte()` simple | Separate estimates of `aggte()` group, dynamic and calendar |
| --- | --- | --- | --- | --- |
| `FALSE` | `FALSE` | analytical standard errors, pointwise normal intervals | the same | the same |
| `TRUE` | `FALSE` | multiplier standard errors, pointwise normal intervals | the same | the same |
| `TRUE` | `TRUE` | multiplier standard errors, simultaneous band over the available cells | multiplier standard error, pointwise normal interval | multiplier standard errors, simultaneous band over the available estimates |
| `FALSE` | `TRUE` | analytical standard errors, pointwise normal intervals, no draws | analytical standard error, pointwise normal interval, no draws | analytical standard errors with a simultaneous critical value from multiplier draws, and `did`'s warning |

These flag conventions describe supported nondegenerate families. Empty or exactly zero families may require no draws; read the effective flags and draw count in the result rather than assuming `biters` draws ran. Separate interval families, including any retained basis contrasts, do not form one simultaneous band over their union.

The overall summaries of the group, dynamic and calendar types always have pointwise normal intervals; their standard errors are multiplier or analytical according to the effective `bstrap`. `aggte()` inherits the settings of the fit unless you override them.

### Multiplier draws

The multiplier bootstrap reuses the estimated influence functions; it does not refit the estimator. Each draw multiplies the influence function rows of each unit, or the row sums of each cluster, by a random sign of plus or minus one (Rademacher multipliers). A multiplier standard error is the interquartile range of the draws divided by the interquartile range of the standard normal distribution, on the scale of the estimate. A simultaneous critical value is the `1 - alp` empirical quantile of the largest absolute standardized draw across the covered estimates, and the band is the estimate plus or minus that value times the standard error. In the `bstrap = FALSE, cband = TRUE` case of the aggregate bins, the draws are standardized by their own interquartile range scale and the resulting critical value is applied to the analytical standard errors. Pointwise intervals use the normal critical value.

The signs come from your session's random number stream, so `set.seed()` before `att_gt()` or `aggte()` reproduces the draws; there is no fixed seed. Parallel workers only multiply the signs drawn in your session, so the results equal those of `cores = 1`.

### Analytical standard errors and clustering

Analytical standard errors come from the complete influence functions, with one row for every unit of the unit table, including units without an observed change. They account for the estimation of the treatment, outcome, observation and imputation models, the response rates and the cohort shares used by the chosen method.

With `clustervars`, influence function rows are summed within clusters, including clusters without any unit that contributes a change, and the multiplier signs are drawn by cluster. Clustering changes the standard errors and bands; the point estimates stay the same, because the GMM weights are based on units. A cluster scheme that `aggte()` cannot honor is refused.

### What the uncertainty assumes

The reported uncertainty rests on large sample approximations: enough independent units or clusters, overlap, the model conditions of the chosen methods and regularity conditions. A completed fit does not establish these conditions, and no finite sample, few cluster or nominal coverage guarantee is made. Computation is conventional, with conditioning and residual checks, and no per fit error certificate is produced.

## Point estimates only

With `compute_inffunc = FALSE`, `att_gt()` returns the same point estimates as with `TRUE`, but `inffunc`, `V_analytical`, `W` and `Wpval` are `NULL`, `se` is `NA`, `c` is the pointwise normal value, the recorded `bstrap` and `cband` are `FALSE`, `fit$chained$inference$kind` is `"not_computed"` and the pretrend test status is `"not_computed"`. `aggte()` refuses such a result.

## Pretrend tests

### Unconditional Wald test

`att_gt()` tests jointly all cells before treatment (role `"pre"`) of the chosen `base_period`, excluding exact zero reference cells, with `did`'s Wald statistic. `W` and `Wpval` are set only when every one of these cells is available and the full covariance passes the existing Cholesky/solve and conditioning policy. Otherwise they are `NULL` and `fit$chained$pretest` records the status (`"unavailable"`, `"not_applicable"` or `"not_computed"`) with its reason. A test on the available subset is never substituted. When the test is computed, `fit$chained$pretest` has `status = "computed"`, the tested `cells`, the degrees of freedom `df`, `W` and `Wpval`.

A test that does not reject does not establish parallel trends. Under PT-all, cells before the chosen start can reject without contradicting the assumption, as explained in [Parallel trends](../methods/parallel-trends.md#pretrend-tests-under-each-choice).

### Conditional test

`conditional_did_pretest()` tests a conditional null: for every cell before treatment of the chosen `base_period` (reference zeros excluded), the mean of the cell's comparison given the `xformla` design row of the unit at the earlier endpoint of its auxiliary date pair is zero. This can be weaker than conditioning on every relevant variable. The statistic is the unstandardized Cramér von Mises (CvM) statistic over a grid with one point for every unit of the unit table. Its multiplier draws use one matrix of Rademacher signs shared by all cells and grid points, at the unit or cluster level. The p value is the share of draws that exceed the statistic, and the critical value is their empirical type-1 `1 - alp` quantile.

| Field | Contents |
| --- | --- |
| `CvM` | test statistic |
| `CvMb` | the statistic in each multiplier draw |
| `CvMcval`, `CvMpval` | critical value and p value |
| `KS`, `KSb`, `KScval`, `KSpval` | `NULL`; this test is not computed |
| `clustervars`, `xformla` | settings of the test |
| `chained` | the null, `status`, `reason`, base period, tested cells with their estimates, route, grid, statistic and multiplier conventions, clusters, random number record, execution and resources |

`status` is `"available"`, `"unavailable"` (with `NA` statistics and a reason) or `"not_applicable"` for a custom `est_method` function. The test is unavailable when a needed cell is unavailable, when a covariate needed on the grid is missing or has a level absent from the fitted design, when there is no cell before treatment to test, when clustered inference has fewer than two clusters, or when the estimated process is exactly zero everywhere. Estimated effects are never affected. The conditional test is separate from `W` and `Wpval`. Its calibration requires uniform process and multiplier regularity conditions, including a growing number of independent clusters for clustered inference.

## Plots

`ggdid()` returns `did`'s `ggplot2` plots with intervals that use the result's critical values. For a group time result, the plot is built on the grid of cohorts and periods; `attr(p, "cdid_presentation_grid")` lists each grid point with its `group`, `time`, `target_id`, `att`, `se`, `status` (`"available"`, `"unavailable"` or `"unrequested"`), `support_status` and `support_reason`. For an aggregated result, `attr(p, "cdid_presentation_bins")` lists `key`, `label`, `att`, `se`, `status` (`"available"`, `"restricted"` or `"unavailable"`) and `reason`. As in `did`, `ggdid()` refuses the simple type; it has no plot. For OR and DR fits and for change imputation and the augmented correction, the plot caption shows the estimator label.

## Saving and exporting

- **Tables.** `tidy()` and `glance()` return data frames that `write.csv()` exports, with unavailable values as `NA`. The overall summary of the group type is the first row of its tidy table; for the dynamic and calendar types use `overall.att` and `overall.se`, as in [Get started](../get-started.md#aggregate-the-effects).
- **Plots.** Save with `ggplot2::ggsave()`.
- **Fitted objects.** `saveRDS()` keeps the complete result, including influence functions, for later inspection with the same package versions; `fit$chained$versions` records them. The panel data are not part of the result.

## Planned result tables

The following capabilities have been proposed for a later version. They require further authoring and, in most cases, new inspection functions whose names and contents are not yet decided. No such function exists now.

| Proposed capability | Intended content | Current state |
| --- | --- | --- |
| Table of two period comparisons | cohort, start and end dates, comparison group, numbers of usable pairs, adjustment methods, comparison estimate, availability and reason; comparisons used for effects kept apart from comparisons reported before treatment | Stored structures differ across estimators; no documented accessor. |
| Combination table for one effect | each comparison's coefficient in the reconstruction of a chosen effect and its contribution, reproducing the reported estimate | No documented accessor for the combination coefficients. |
| Aggregation weights table | effect identifiers, weights, included and excluded cells and the resulting target | Partly available in `agg$chained$cell`, `agg$chained$overall$weights`, `agg$chained$cohort_shares` and `agg$chained$weights_theta`. |
| Consistent model diagnostics | observation probabilities and coefficients, treatment and outcome fits, overlap summaries | Detailed nested records exist in `fit$chained$nuisance_diagnostics`. |

These tables will need to keep three kinds of weights apart. Observation weights adjust for missing outcomes. Combination coefficients reconstruct an effect from comparisons; they can be negative and are not population shares. Aggregation weights combine effects into summaries.

---

**Reading path, step 19 of 20** · Previous: [Options](options.md) · Next: [Diagnose errors and unavailable results](../help/troubleshooting.md)

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
| 18 | [Options](options.md) |
| 19 | **Results** (this page) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
