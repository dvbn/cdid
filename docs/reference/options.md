# Options reference

This page is the guide's single definition of each argument: its default, accepted values, when it is active and how it interacts with other choices. The help pages `?cdid::att_gt`, `?cdid::aggte` and `?cdid::conditional_did_pretest` contain the complete text. The articles link here instead of repeating these definitions.

## Signatures

These signature blocks describe the interface; they are not stand-alone example calls. Worked examples elsewhere in the guide show the corresponding output.

```r
cdid::att_gt(
  yname, tname, idname = NULL, gname, xformla = NULL, data,
  panel = TRUE, allow_unbalanced_panel = NULL,
  control_group = c("nevertreated", "notyettreated"), anticipation = 0,
  weightsname = NULL, fix_weights = NULL, alp = 0.05, bstrap = TRUE,
  cband = TRUE, biters = 1000, clustervars = NULL, est_method = "dr",
  base_period = "varying", faster_mode = TRUE, print_details = FALSE,
  pl = FALSE, cores = 1, compute_inffunc = TRUE, ...,
  mode = NULL, frame = NULL, selection = "none",
  selection_model = "pair_logit", selection_pooling = "separate",
  selection_xformla = NULL, selection_universe = "F", selection_rho = NULL,
  selection_justification = NULL, selection_est_method = "ipw",
  imputation_pooling = "separate", gmm_weighting = c("optimal", "identity"),
  gmm_covariance = "tikhonov", pt_assumption = c("post", "all"),
  pt_start = NULL, control = NULL
)

cdid::aggte(
  MP, type = "group", balance_e = NULL, min_e = -Inf, max_e = Inf,
  na.rm = FALSE, bstrap = NULL, biters = NULL, cband = NULL, alp = NULL,
  clustervars = NULL
)

cdid::conditional_did_pretest(
  yname, tname, idname = NULL, gname, xformla = NULL, data,
  panel = TRUE, allow_unbalanced_panel = NULL,
  control_group = c("nevertreated", "notyettreated"), weightsname = NULL,
  alp = 0.05, bstrap = TRUE, cband = TRUE, biters = 1000, clustervars = NULL,
  est_method = "ipw", print_details = FALSE, pl = FALSE, cores = 1,
  anticipation = 0, base_period = "varying", fix_weights = NULL,
  mode = NULL, frame = NULL, selection = "none",
  selection_model = "pair_logit", selection_pooling = "separate",
  selection_xformla = NULL, selection_universe = "F", selection_rho = NULL,
  selection_justification = NULL, selection_est_method = "ipw",
  imputation_pooling = "separate", control = NULL
)
```

In `att_gt()`, arguments before `...` retain the names and order of `did` 2.5.1. The deliberate default difference is `allow_unbalanced_panel = NULL`, a mode-dependent sentinel that maps to `did`'s `FALSE` in ordinary mode. The arguments after `...` are used only in chained mode; ordinary mode neither evaluates nor forwards them.

## Data, timing and comparison group

| Argument | Default | Accepted values | Notes for chained mode |
| --- | --- | --- | --- |
| `yname` | required | name of a numeric outcome column | Observed values must be finite. Keep missing outcomes as `NA`. |
| `tname` | required | name of a finite numeric date column | The stored calendar is the sorted set of dates present in the data. OR/imputation durations count positions in it; event time and anticipation cutoffs use numeric date values. |
| `idname` | `NULL` | name of the unit identifier column | Required, because units are aligned by identity. |
| `gname` | required | name of the column holding each unit's first treatment period | Known, constant within each panel unit. Use a finite positive treatment date on the time scale of `tname`, or 0 for never treated. The treatment date need not be a stored row. |
| `data` | required | a data frame with one row per unit and date | Each ID/date pair must be unique. A balanced panel and an observed reference-date outcome for every unit are not required. See [Prepare your panel](../help/prepare-data.md) and [Investigate missingness](../help/investigate-missingness.md). |
| `panel` | `TRUE` | `TRUE` or `FALSE` | Chained mode requires `TRUE`. With `mode = NULL`, `FALSE` selects ordinary mode. |
| `allow_unbalanced_panel` | `NULL` | `NULL`, `TRUE` or `FALSE` | `NULL` retains the full unit table, including units with no observed change. Units with repeated observations supply changes. `FALSE` is refused. In ordinary mode `NULL` means `did`'s default `FALSE`. |
| `control_group` | `"nevertreated"` | `"nevertreated"` or `"notyettreated"` | |
| `anticipation` | 0 | finite nonnegative whole number | Sets each cohort's reference period: the latest stored date d with `d + anticipation < g`. |
| `weightsname` | `NULL` | name of a column of positive sampling weights | Weights that change within a unit need `fix_weights`, see [Changing covariates and weights](#changing-covariates-and-weights). |
| `fix_weights` | `NULL` | `NULL`, `"first_period"` or `"base_period"` | `"varying"` is refused. |

The unit table `frame` is described with the [chained arguments](#chained-arguments) and in [Handling missing outcomes](../methods/missing-outcomes.md#the-study-population-and-the-unit-table).

## Treatment adjustment

| Argument | Default | Accepted values | Notes for chained mode |
| --- | --- | --- | --- |
| `est_method` | `"dr"` | `"dr"`, `"ipw"`, `"reg"` or a custom IPW function | Any other string is refused before the data are read; nothing is substituted. |
| `xformla` | `NULL` | one sided formula, such as `~ x1 + x2` | The design must have at least one column: `~ 1` and `~ 0 + x` are accepted, `~ 0` is refused. Covariates must be fixed within units or are read at reference dates. |
| `...` | | | Passed to a custom `est_method` function; refused when `est_method` is a string. |

A custom function has `did`'s panel IPW interface `f(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...)` and returns a list with a finite scalar `ATT` and a finite vector or one-column matrix `att.inf.func` with one entry per supplied row, in that order. Under observation weighting with `selection = "conditional_mtar"`, it must also return `log_weight_gradient`. [Treatment adjustment](../methods/treatment.md#custom-ipw-functions) explains these requirements and the models behind the three string methods.

## Chained arguments

| Argument | Default | Accepted values | Active when |
| --- | --- | --- | --- |
| `mode` | `NULL` | `NULL`, `"chained"` or `"ordinary"` | Always. `NULL` selects chained mode for `panel = TRUE` and ordinary mode otherwise. |
| `frame` | `NULL` | the unit table: a data frame with one row per study unit, matched to `data` by `idname` | Chained mode. Without it, the units in `data` form the study population. |
| `selection` | `"none"` | `"none"`, `"conditional_mtar"`, or the legacy name `"corollary1"` | Chained mode. `"corollary1"` means `"conditional_mtar"` with `selection_model = "pair_logit"`. |
| `selection_est_method` | `"ipw"` | `"ipw"` (observation weighting), `"imputation"` (change imputation) or `"dr"` (augmented correction) | `"imputation"` and `"dr"` require `selection = "conditional_mtar"`; `"dr"` also requires `est_method = "dr"`. |
| `selection_model` | `"pair_logit"` | `"pair_logit"` or `"independent_period_logit"` | An active conditional-MTAR correction with `selection_est_method` equal to `"ipw"` or `"dr"`. Never evaluated under `"imputation"`. A nondefault value with `selection = "none"` is refused. |
| `selection_pooling` | `"separate"` | `"separate"`, `"within_cohort"` or `"across_cohorts"` | As for `selection_model`. |
| `selection_xformla` | `NULL`, which inherits `xformla` | one sided formula | An observation model is fitted (`"ipw"` or `"dr"`). Never evaluated under `"imputation"`. |
| `selection_universe` | `"F"` | `"F"` (whole unit table), `"B"` (units observed at least once) or `"A"` (units observed at least twice) | `"B"` and `"A"` only with observation weighting and `selection_model = "pair_logit"`. |
| `selection_rho` | `NULL` | `NULL` or the name of a unit table column of known normalized density ratios | Only with `"B"` or `"A"`. |
| `selection_justification` | `NULL` | a nonempty character string | Required with `"B"` or `"A"`. |
| `imputation_pooling` | `"separate"` | `"separate"`, `"within_cohort"` or `"across_cohorts"` | `selection_est_method` is `"imputation"` or `"dr"`; otherwise never evaluated. |
| `gmm_weighting` | `"optimal"` | `"optimal"`, its alias `"two_step"`, or `"identity"` | Chained mode. |
| `gmm_covariance` | `"tikhonov"` | `"tikhonov"` | Validated in chained mode; the handler is used under `gmm_weighting = "optimal"`. |
| `pt_assumption` | `"post"` | `"post"` or `"all"` | Chained mode. |
| `pt_start` | `NULL` | one finite numeric date of the stored calendar | Required under `"all"`; a non-`NULL` value is refused under `"post"`. |
| `control` | `NULL` | a named list, see [The control list](#the-control-list) | Chained mode. |

## Handling missing outcomes

[Handling missing outcomes](../methods/missing-outcomes.md) explains the assumptions and models. The interface rules are these.

- **Formulas.** `xformla` and `selection_xformla` are one sided formulas whose `stats::terms` intercept and assignment information is respected, including intercept only (`~ 1`) and no intercept (`~ 0 + x`) forms. Only `selection_xformla` accepts a design with no columns (`~ 0`); its observation model then has no parameters and uses the fixed probability 0.5 for each pair, or 0.25 for a product of two dates.
- **Observation weighting** (`selection_est_method = "ipw"`) fits `selection_model` with `selection_pooling` on the fitting sample selected by `selection_universe`. `"independent_period_logit"` requires `selection_universe = "F"` and the independence of the two dates given cohort and covariates.
- **Change imputation** (`"imputation"`) requires `selection = "conditional_mtar"`, `selection_universe = "F"`, `selection_rho = NULL`, known cohorts and complete `xformla` covariates on every row of the unit table, and `est_method` equal to `"ipw"`, `"reg"` or `"dr"`. It fits no observation model.
- **Augmented correction** (`"dr"`) has the requirements of change imputation, requires `est_method = "dr"`, needs complete covariates of both formulas on every row of the unit table and uses every observation model and imputation argument. The two pooling arguments are independent.
- **Not supported:** known pair probabilities (including a design probability of one), estimated or unnormalized transport ratios, and probit models.

## Parallel trends and reporting

| Argument | Default | Meaning |
| --- | --- | --- |
| `pt_assumption` | `"post"` | `"post"`: parallel trends from each cohort's reference period onward. `"all"`: the same restrictions plus restrictions on untreated trends from `pt_start` up to each reference period. |
| `pt_start` | `NULL` | Under `"all"`, one finite numeric date that equals a stored calendar date exactly, common to all cohorts and chosen before looking at results. Omitting it, `NULL` and dates outside the calendar are errors. |
| `base_period` | `"varying"` | Reporting of cells before treatment, with `did`'s meaning: `"varying"` compares each stored period after the first with the preceding stored period; `"universal"` compares every preperiod with the reference period, including its exact-zero reference cell. Independent of `pt_assumption`. |
| `anticipation` | 0 | A nonnegative whole number used in the numeric-date cutoff `d + anticipation < g`; it moves the reference period back. This is not the stored-index duration of the outcome models. |

[Parallel trends](../methods/parallel-trends.md) explains PT-post, PT-all and the case where PT-all adds no restriction.

## GMM weighting and covariance regularization

Chained mode combines the available comparisons by GMM. `gmm_weighting = "optimal"`, the default, uses feasible two step weighting based on the estimated covariance of the comparisons; `"two_step"` is an alias for it. `gmm_weighting = "identity"` requests unweighted estimation explicitly. Identity weighting is never substituted automatically.

Before the covariance is checked, the package applies exact reductions that remove proven redundancies among the comparisons. When the reduced covariance is well conditioned, it is used as it is. When it is singular or nearly singular, judged by the ratio of its smallest to its largest eigenvalue, or when it is known to be rank deficient or its Cholesky factorization fails, the selected handler `gmm_covariance = "tikhonov"` adds a ridge proportional to the trace of the reduced covariance before it is inverted for the weights. At exactly zero covariance the criterion is full-catalog least squares with the actual zero sandwich variance. Invalid or materially indefinite covariance inputs still refuse. No comparison, unit or identified contrast is dropped. `"tikhonov"` is the only accepted value of `gmm_covariance`. The weighting request, whether the handler was invoked and why, its tuning and the conditioning are recorded in `fit$chained$covariance_handling`. The advanced entries `condition_trigger` and `ridge_relative` of the control list belong to this step.

Covariance regularization affects only how comparisons are weighted. It supplies no identifying information: an effect that the comparisons cannot identify remains unavailable. Standard errors and bands always use the complete influence functions and the sandwich variance, never the inverse of the regularized matrix. Feasible weighting does not make the estimates unbiased in finite samples.

## Inference and computation

| Argument | Default | Meaning in chained mode |
| --- | --- | --- |
| `alp` | 0.05 | Significance level of intervals and bands. |
| `bstrap` | `TRUE` | Standard errors from multiplier bootstrap draws instead of analytical ones. |
| `cband` | `TRUE` | Simultaneous band over the available cells where requested. |
| `biters` | 1000 | Number of multiplier draws. |
| `clustervars` | `NULL` | At most one cluster variable in addition to `idname`, with a label for every unit of the unit table. Clustering changes standard errors and bands; the point estimates stay the same. |
| `compute_inffunc` | `TRUE` | `FALSE` returns point estimates only, which equal those of `TRUE`: no influence functions, standard errors, bands or pretrend test, and `aggte()` refuses the result. |
| `pl`, `cores` | `FALSE`, 1 | `pl = TRUE, cores = k` computes the multiplier draw products on `k` worker processes, which must be declared with `control = list(workers = k)`; otherwise the request is refused before the data are read. Signs are drawn in your session, so the results equal those of `cores = 1`. Estimation itself stays serial. |
| `print_details` | `FALSE` | `TRUE` reports progress through messages. |
| `faster_mode` | `TRUE` | Accepted with `did`'s meaning; chained mode has one computation engine. |

The exact combinations of `bstrap` and `cband` at each output are listed in the [results reference](results.md#uncertainty).

## The control list

`control` is an advanced list used only in chained mode. `NULL` and an empty list give the defaults below, which the result records.

| Entry | Default | Meaning |
| --- | --- | --- |
| `numeric_workspace_bytes` | 1073741824 (1 GiB) | Admission limit for the counted numeric workspace of a fit. |
| `workers` | 1 | Number of worker processes that `pl = TRUE, cores = k` may use for multiplier draws. It is never raised automatically, and the workers' buffers count against the same workspace. |
| `counted_work_limit` | `Inf` | The estimated work cap is off by default. A finite positive whole number refuses, before a stage runs, a fit whose conservative work estimate exceeds it; an explicit `Inf` also turns it off. Zero, negative, `NA` and `NaN` values, numbers that are not whole, values out of range and nonnumeric values are refused. |
| `chunk_rows` | chosen automatically | Number of rows per processing chunk. |
| `condition_trigger` | `1e-10` | Finite, strictly between zero and one; used by the covariance regularization step. |
| `ridge_relative` | `1e-4` | Finite and positive; used by the covariance regularization step. |

Unknown entries are refused. The entry `relative_accuracy` is refused with class `cdid_chained_public_accuracy_deferred`, because a guaranteed accuracy for each fit is not provided. The fitted result records the public cap and the internal arithmetic ceiling separately, in `fit$chained$public$counted_work_limit` and `fit$chained$public$work_cap`.

Work estimates are conservative diagnostics, not measured run time. These controls limit counted numeric storage and estimated work; they do not limit the memory of the R process, CPU time or the work of a whole call, and the package does not detect the memory or speed of your machine.

## Supported combinations

With `selection = "conditional_mtar"`, seven built-in combinations of treatment adjustment and adjustment for missing outcomes are supported. The qualified custom IPW extension is listed separately in the last row:

| `est_method` | `selection_est_method = "ipw"` | `"imputation"` | `"dr"` |
| --- | --- | --- | --- |
| `"ipw"` | supported | supported | refused |
| `"reg"` | supported | supported | refused |
| `"dr"` | supported | supported | supported |
| custom IPW function | supported with `log_weight_gradient` | refused | refused |

With `selection = "none"`, `selection_est_method` stays at `"ipw"` and all of `"ipw"`, `"reg"`, `"dr"` and a custom function are supported. `selection_universe = "B"` or `"A"` is available only in the first column, with `selection_model = "pair_logit"`. Refused requests stop with an error before estimation; no other combination is substituted.

## Changing covariates and weights

Covariates and sampling weights that are constant within each unit and agree with the unit table pass the usual checks and are used as unit characteristics. For such calls, `fix_weights = NULL`, `"first_period"` and `"base_period"` give identical results.

**Weights that change within a unit** require a fixed reference:

- `"first_period"`: every comparison uses each unit's weight at the first stored date.
- `"base_period"`: every unit in the comparisons for cohort g, controls included, uses its weight at the cohort's reference period.
- `NULL` would mean weights specific to each date, as in `did`, and is refused.

**Covariates that change within a unit**, or a formula that uses the time variable, are read at reference dates, following `did`'s convention of using covariate values from the base period of each comparison. Every comparison used for the effects of cohort g reads each unit's values at the cohort's reference period, and a comparison reported before treatment reads them at its earlier date. Values come only from the actual unit and date rows, including rows whose outcome is missing, and are required only where an estimate reads them. A missing required value is refused with its unit, date, variable and purpose; values are never carried forward or substituted, and a unit table value cannot override values in `data` that change. This handling is available with every chained method, with or without corrections for missing outcomes.

**Aggregation of these fits** weights cohorts by their masses at the first stored date, and the standard errors include the uncertainty of these shares. Units never treated need no weight at that date. A share-dependent aggregate is unavailable when a required treated-cohort first-date mass is incomplete, with the unit identifiers recorded. One-cohort bins and a one-bin group overall summary cancel the shares exactly and do not require that mass. With `"base_period"`, each cell keeps the weights of its own cohort's comparisons, so an aggregate is a declared mixture of cohort effects and does not describe a single weighted population.

## Aggregation arguments

| Argument | Default | Meaning for chained results |
| --- | --- | --- |
| `MP` | required | A result of `att_gt()`. Chained results are aggregated from their stored state without refitting; ordinary results go to `did::aggte()` unchanged. |
| `type` | `"group"` | `"simple"`, `"group"`, `"dynamic"` or `"calendar"`. |
| `balance_e`, `min_e`, `max_e` | `NULL`, `-Inf`, `Inf` | Event time selectors with `did`'s conventions. For the calendar type they are ignored with `did`'s warning. The result records how they were applied. |
| `na.rm` | `FALSE` | `FALSE`: an aggregate that needs an unavailable cell is `NA`, and the missing cells are named. `TRUE`: the aggregate is computed over the available cells and labelled as a restricted target. |
| `bstrap`, `biters`, `cband`, `alp`, `clustervars` | `NULL` | `NULL` inherits the value used by `att_gt()`. An explicit value overrides it with `did`'s meaning. A cluster override that cannot be honored is refused. |

Selector details matter: simple and group aggregation use only `max_e`, measured on the sorted union of calendar and cohort dates; `min_e` and `balance_e` have no effect there. Dynamic aggregation uses `e = t - g` in original date units and applies its window and cohort-balance rules. Calendar aggregation ignores all three selectors with a warning. These selectors change the reported target, not the fitted comparison catalog.

A result fitted with `compute_inffunc = FALSE` is refused with class `cdid_chained_request_no_inffunc_error`. Multiplier draws use your session's random number stream, so `set.seed()` before the call reproduces them, and workers admitted by the fit are reused.

## Conditional pretrend test arguments

`conditional_did_pretest()` keeps `did`'s 19 arguments in order and appends the chained settings it needs. For panel data, `mode = NULL` runs the chained conditional test; `mode = "ordinary"` and repeated cross sections go to `did::conditional_did_pretest()`.

| Argument | Meaning in chained mode |
| --- | --- |
| `est_method` | `"ipw"` (default), `"reg"` or `"dr"`. A custom function returns a result with status `"not_applicable"`. |
| `xformla` | Its design columns are the coordinates of the test; variables that appear only in `selection_xformla` are not. |
| `bstrap`, `cband` | Validated only. The test always uses multiplier draws and reports no standard errors. |
| `biters`, `alp` | Number of multiplier draws and level of the reported critical value. |
| `pl`, `cores` | Parallel draws when `control$workers` declares them; results equal those of `cores = 1`. |
| `anticipation`, `base_period`, `fix_weights`, `frame`, the `selection` arguments, `imputation_pooling`, `control` | The chained settings of `att_gt()`, with the same meanings; they select the fitted route whose comparisons before treatment are tested. |

GMM weighting and the PT-all settings are not arguments of the test.

## Refusals

Refusals are errors that inherit class `cdid_error`. The most common classes are:

| Class | Cause |
| --- | --- |
| `cdid_mode_error` | invalid `mode`, or invalid `panel` while resolving `mode = NULL` |
| `cdid_chained_panel` | chained mode with `panel` not equal to `TRUE` |
| `cdid_chained_balanced` | explicit `allow_unbalanced_panel = FALSE` in chained mode |
| `cdid_chained_unavailable` | an unsupported request, such as an unknown `est_method` string (also `cdid_chained_request_est_method`) or the augmented correction with `est_method` other than `"dr"` (also `cdid_chained_request_unsupported_selection`) |
| `cdid_chained_request_*` | other unsupported requests or inputs, for example invalid `pt_start` values |
| `cdid_chained_request_parallel_error` | `pl = TRUE, cores = k` beyond the declared `control$workers` |
| `cdid_chained_public_control_error` | an invalid `control` list |
| `cdid_custom_ipw_error` | malformed or nonfinite output of a custom IPW function |
| `cdid_chained_request_no_inffunc_error` | `aggte()` on a result fitted with `compute_inffunc = FALSE` |
| `cdid_untested_upstream` | loading the package with `did` older than 2.5.1 or `DRDID` older than 1.3.0 |

Fits can also stop because a numeric workspace or work limit would be exceeded, or because a model fit, covariance or solve is invalid. Such failures concern the computation; they are not conclusions about whether an effect is identified in the population. Nothing is dropped or relaxed to make a fit succeed.

## Tested versions and dependencies

The package requires `did` 2.5.1 or later and `DRDID` 1.3.0 or later, and refuses to load with earlier versions. It has been tested with `did` 2.5.1 and `DRDID` 1.3.0. Later versions are accepted, but it is not guaranteed that everything works with them or that results are unchanged; a result records the `did` version it was computed with in `fit$chained$versions`. The metadata declares R 4.1.0 or later and imports `did`, `DRDID` and the base package `parallel`; that declaration is not a completed minimum-version or platform-support validation. `ggplot2`, which `did` uses for its plots, is needed for `ggdid()`. The [installation](../../README.md#installation) instructions use `remotes::install_github("dvbn/cdid")`; the optional [pinned validation setup](../../validation/SETUP.md#optional-pinned-validation-setup) installs exactly the tested versions.

---

**Reading path, step 18 of 20** · Previous: [Inference](../examples/inference.md) · Next: [Results](results.md)

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
| 18 | **Options** (this page) |
| 19 | [Results](results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
