# Impute actual changes for missing outcomes

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

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
| Slopes | shared as `imputation_pooling` declares | one vector for all blocks, with no pooling argument |

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

## Which models must be correct

| Fit | Missing outcome layer | Treatment layer |
| --- | --- | --- |
| `fit_weighting` | `q` correct in every needed cell, with independence of the two dates | `p` or `m` correct |
| `fit_imputation` | `j` correct in every needed cell | `p` or `m` correct |
| `fit_augmented` | for each needed cell, `q` or `j` correct, possibly different across cells | `p` or `m` correct |

Every fit also needs conditional mean trends missing at random given a scientifically sufficient conditioning set W, here `x`; PT-post and no anticipation; positive cohort masses, positive observation probabilities where `q` is fitted and positive training mass with prediction support where `j` is; valid moments for every comparison used; and regular finite fits of every active model, alongside the usual large sample regularity conditions. For the augmented route, W includes the covariates in `xformla` and any variables only in `selection_xformla`; a correct `j` must represent the actual-change mean on that W. Correctness must hold across the cells sharing the relevant slopes, both endpoints of an active date model and all blocks of a shared counterfactual regression. Different predictor lists need separately sufficient treatment and observation restrictions. Pure imputation needs a correct `j` even when `p` and `m` are both correct: treatment DR does not generally repair a wrong imputer. These conditions are sufficient, not necessary, and they concern large sample behavior. There is no general finite sample unbiasedness guarantee, and no fit is guaranteed the lowest standard error. Here every model is correct, so all three fits are expected to target the effects. [Compare model robustness](robustness.md#correctness-conditions-across-both-layers) gives the complete conditions.

## Three kinds of slope sharing

`imputation_pooling` shares slopes of `j`, and `selection_pooling` shares slopes of `q`; under the augmented correction they act independently. The counterfactual regression `m` always shares one slope vector and has no pooling argument.

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
nuisance <- fit_imputation$chained$nuisance_diagnostics
nuisance$fitted_catalog$link_id <- NULL  # opaque row identifiers, not needed here
str(nuisance, max.level = 2)
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

List of 26
 $ status                     : chr "available"
 $ method                     : chr "imputation"
 $ selection_est_method       : chr "imputation"
 $ est_method                 : chr "dr"
 $ imputation_pooling         : chr "separate"
 $ imputation_model           : chr "actual-change imputation: observed pair changes on S_e = 1 with fixed unit weights w, local cohort/pair interce"| __truncated__
 $ counterfactual_model       : chr "shared one-period counterfactual control outcome-change regression on the imputed synthetic control changes (lo"| __truncated__
 $ treatment_model            : chr "risk-set treatment propensity logit with synthetic weights w*a_k on F"
 $ fit_structure              : chr "one WLS fit per imputation pooling group (30 group(s), 30 cell(s)), one joint shared counterfactual OR fit, 15 "| __truncated__
 $ duration_convention        : chr "stored-calendar transition count L_e = i_t - i_s"
 $ selection                  : chr "conditional_mtar"
 $ selection_route            :List of 25
  ..$ method                : chr "imputation"
  ..$ canonical_assumption  : chr "conditional_mtar"
  ..$ requested_selection   : chr "conditional_mtar"
  ..$ legacy_alias          : logi FALSE
  ..$ selection_est_method  : chr "imputation"
  ..$ est_method            : chr "dr"
  ..$ imputation_pooling    : chr "separate"
  ..$ observation_model     : chr "none: no observation-probability model is fitted"
  ..$ inactive_arguments    : chr [1:3] "selection_model" "selection_pooling" "selection_xformla"
  ..$ selection_model       : chr NA
  ..$ selection_pooling     : chr NA
  ..$ fitting_universe      : chr "F"
  ..$ population            : chr "supplied representative frame"
  ..$ formula               :Class 'formula'  language ~x
  .. .. ..- attr(*, ".Environment")=<environment: 0xca0505b60> 
  ..$ transport             : chr "rho=1"
  ..$ justification         : NULL
  ..$ imputation_model      : chr "actual-change imputation: observed pair changes on S_e = 1 with fixed unit weights w, local cohort/pair interce"| __truncated__
  ..$ counterfactual_model  : chr "shared one-period counterfactual control outcome-change regression on the imputed synthetic control changes (lo"| __truncated__
  ..$ treatment_model       : chr "risk-set treatment propensity logit with synthetic weights w*a_k on F"
  ..$ duration_convention   : chr "stored-calendar transition count L_e = i_t - i_s"
  ..$ empirical_repeat_rates: Named num [1:2] 0.915 0.932
  .. ..- attr(*, "names")= chr [1:2] "3" "Inf"
  ..$ assumptions           : chr "Each needed actual cohort/pair cell requires actual-pair conditional mean trend-MAR on F, E[Delta_e | G=k, X, S"| __truncated__
  ..$ if_scale              : chr "phi_F on ALL N roster rows; core psi=(n_A/N)*phi_F with 1/n_A covariance normalization"
  ..$ unsupported           : chr [1:5] "selection_universe B/A" "transport columns (rho)" "corrected custom IPW callbacks" "augmented selection_est_method = \"dr\"" ...
  ..$ models                :List of 17
 $ inactive_selection         :List of 9
  ..$ selection_est_method    : chr "imputation"
  ..$ active                  : logi FALSE
  ..$ observation_model_fitted: logi FALSE
  ..$ selection_formula_forced: logi FALSE
  ..$ selection_model         : chr NA
  ..$ selection_pooling       : chr NA
  ..$ selection_xformla       : NULL
  ..$ universe                : chr "F"
  ..$ note                    : chr "Pure imputation fits no observation-probability model; selection_model, selection_pooling, and selection_xforml"| __truncated__
 $ outcome_recipe             :List of 9
  ..$ columns        : chr [1:2] "(Intercept)" "x"
  ..$ intercept      : logi TRUE
  ..$ intercept_index: int 1
  ..$ slope_index    : int 2
  ..$ slope_columns  : chr "x"
  ..$ p              : int 1
  ..$ duration       : chr "L_e = t_calendar_index - s_calendar_index on the stored modeled calendar"
  ..$ design         : chr "d_e = L_e * z(X) with fixed-unit z"
  ..$ key            : chr "fixed-X-summed|intercept:TRUE|z:x"
  ..- attr(*, "class")= chr "cdid_or_pooled_recipe"
 $ treatment_recipe           :List of 3
  ..$ columns: chr [1:2] "(Intercept)" "x"
  ..$ index  : int [1:2] 1 2
  ..$ p      : int 2
  ..- attr(*, "class")= chr "cdid_dr_treatment_recipe"
 $ cells                      :'data.frame':	30 obs. of  11 variables:
  ..$ cell_id         : int [1:30] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ name            : chr [1:30] "k0.e12" "k0.e23" "k0.e34" "k0.e45" ...
  ..$ cohort          : num [1:30] Inf Inf Inf Inf Inf ...
  ..$ s               : num [1:30] 1 2 3 4 5 1 2 3 4 1 ...
  ..$ t               : num [1:30] 2 3 4 5 6 3 4 5 6 4 ...
  ..$ duration_steps  : int [1:30] 1 1 1 1 1 2 2 2 2 3 ...
  ..$ group_id        : int [1:30] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ n_observed_pairs: int [1:30] 235 265 325 346 352 265 294 336 346 268 ...
  ..$ mass            : num [1:30] 0.117 0.133 0.163 0.173 0.176 ...
  ..$ local_intercept : num [1:30] 0.1284 0.1714 0.0492 0.0445 0.0812 ...
  ..$ repeat_rate     : num [1:30] 0.932 0.932 0.932 0.932 0.932 ...
 $ groups                     :List of 30
  ..$ group0 :List of 12
  ..$ group1 :List of 12
  ..$ group2 :List of 12
  ..$ group3 :List of 12
  ..$ group4 :List of 12
  ..$ group5 :List of 12
  ..$ group6 :List of 12
  ..$ group7 :List of 12
  ..$ group8 :List of 12
  ..$ group9 :List of 12
  ..$ group10:List of 12
  ..$ group11:List of 12
  ..$ group12:List of 12
  ..$ group13:List of 12
  ..$ group14:List of 12
  ..$ group15:List of 12
  ..$ group16:List of 12
  ..$ group17:List of 12
  ..$ group18:List of 12
  ..$ group19:List of 12
  ..$ group20:List of 12
  ..$ group21:List of 12
  ..$ group22:List of 12
  ..$ group23:List of 12
  ..$ group24:List of 12
  ..$ group25:List of 12
  ..$ group26:List of 12
  ..$ group27:List of 12
  ..$ group28:List of 12
  ..$ group29:List of 12
 $ counterfactual_or          :List of 9
  ..$ rank                  : int 1
  ..$ kept_columns          : chr "x"
  ..$ shared_slopes         : Named num 0.292
  .. ..- attr(*, "names")= chr "x"
  ..$ has_local_intercepts  : logi TRUE
  ..$ local_block_intercepts: Named num [1:15] 0.0862 0.1023 0.1029 0.1029 0.0913 ...
  .. ..- attr(*, "names")= chr [1:15] "e12:K={0}" "e23:K={0}" "e34:K={0}" "e45:K={0}" ...
  ..$ blocks                :'data.frame':	15 obs. of  7 variables:
  ..$ aliased_columns       : chr(0) 
  ..$ rss                   : num 23.3
  ..$ objective             : num 0.0117
 $ treatment                  :List of 15
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
  ..$ :List of 8
 $ coefficients               :List of 6
  ..$ imputation_slopes          :List of 30
  ..$ imputation_local_intercepts: Named num [1:30] 0.1284 0.1714 0.0492 0.0445 0.0812 ...
  .. ..- attr(*, "names")= chr [1:30] "k0.e12" "k0.e23" "k0.e34" "k0.e45" ...
  ..$ has_local_intercepts       : logi TRUE
  ..$ shared_slopes              : Named num 0.292
  .. ..- attr(*, "names")= chr "x"
  ..$ local_intercepts           : Named num [1:15] 0.0862 0.1023 0.1029 0.1029 0.0913 ...
  .. ..- attr(*, "names")= chr [1:15] "e12:K={0}" "e23:K={0}" "e34:K={0}" "e45:K={0}" ...
  ..$ treatment                  :List of 15
 $ fitted_catalog             :'data.frame':	16 obs. of  16 variables:
  ..$ index          : int [1:16] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ role           : chr [1:16] "catalog" "catalog" "catalog" "catalog" ...
  ..$ cohort         : num [1:16] 3 3 3 3 3 3 3 3 3 3 ...
  ..$ controls       : chr [1:16] "0" "0" "0" "0" ...
  ..$ s              : num [1:16] 1 1 1 1 1 2 2 2 2 3 ...
  ..$ t              : num [1:16] 2 3 4 5 6 3 4 5 6 4 ...
  ..$ duration_steps : int [1:16] 1 2 3 4 5 1 2 3 4 1 ...
  ..$ target_cell    : int [1:16] 16 21 25 28 30 17 22 26 29 18 ...
  ..$ or_block       : int [1:16] 1 6 10 13 15 2 7 11 14 3 ...
  ..$ treatment_model: int [1:16] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ mu             : num [1:16] -0.0386 1.1789 1.5004 1.6939 1.9489 ...
  ..$ mu_ipw         : num [1:16] -0.0386 1.1789 1.5004 1.6939 1.9489 ...
  ..$ mu_reg         : num [1:16] -0.0389 1.1788 1.5002 1.694 1.9485 ...
  ..$ mu_dr          : num [1:16] -0.0386 1.1789 1.5004 1.6939 1.9489 ...
  ..$ d_plus         : num [1:16] 0.452 0.452 0.452 0.452 0.452 ...
  ..$ d_minus        : num [1:16] 0.451 0.451 0.451 0.451 0.451 ...
 $ support                    :List of 11
  ..$ n_full_frame_units              : int 2000
  ..$ n_pair_eligible_A_units         : int 1847
  ..$ n_frame_only_units              : int 0
  ..$ design_source                   : chr "X_frame"
  ..$ design_columns                  : chr [1:2] "(Intercept)" "x"
  ..$ design_rows                     : int 2000
  ..$ prediction_rows                 : chr "every F row of each needed cohort (prediction standardization, not fill-only)"
  ..$ observation_probability_models  : int 0
  ..$ observation_probabilities_fitted: logi FALSE
  ..$ if_rows                         : chr "complete original-F rows, including off-A and frame-only units"
  ..$ scale                           : chr "core psi = (n_A/N) phi_F with covariance divisor n_A; clusters sum original-unit rows"
 $ accounting                 :List of 34
  ..$ owned_entries                : num 201289
  ..$ character_bytes              : num 98538
  ..$ borrowed_entries             : num 1405746
  ..$ retained_external_entries    : num 1254436
  ..$ design_owned_entries         : num 6443
  ..$ design_borrowed_entries      : num 0
  ..$ external_work_ops            : num 26296576
  ..$ retained_bound               : num 1611262
  ..$ fit_scratch_bound            : num 1025792
  ..$ work_ops_bound               : num 8.81e+09
  ..$ automatic_chunk_rows         : logi TRUE
  ..$ n_cells                      : int 30
  ..$ n_groups                     : int 30
  ..$ n_imputation_fits            : int 30
  ..$ n_imputation_local_intercepts: int 30
  ..$ n_imputation_slopes_by_group : int [1:30] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ n_or_fits                    : int 1
  ..$ n_or_blocks                  : int 15
  ..$ n_or_local_intercepts        : int 15
  ..$ n_or_shared_slopes           : int 1
  ..$ n_treatment_models           : int 15
  ..$ n_treatment_fits             : int 15
  ..$ n_treatment_columns          : int 2
  ..$ n_selection_fits             : int 0
  ..$ n_selection_cells            : int 0
  ..$ n_selection_blocks           : int 0
  ..$ n_selection_parameters       : int 0
  ..$ observation_logit_models     : int 0
  ..$ full_F_masks_per_model       : int 0
  ..$ full_F_prediction_matrices   : int 0
  ..$ full_F_by_link_IF_retained   : logi FALSE
  ..$ inverse_matrices_retained    : logi FALSE
  ..$ dense_dummy_regression       : logi FALSE
  ..$ peak_rss_measured            : logi FALSE
 $ frame_counts               :List of 2
  ..$ n_full_frame_units     : int 2000
  ..$ n_pair_eligible_A_units: int 1847
 $ auxiliary_descriptor_status:'data.frame':	1 obs. of  7 variables:
  ..$ index : int 1
  ..$ g     : num 3
  ..$ s     : num 1
  ..$ t     : num 2
  ..$ K     : chr "{Inf}"
  ..$ status: chr "supported"
  ..$ reason: chr NA
 $ interpretation             : chr "Pure actual-change imputation: the imputation fits use observed pair changes only, and the counterfactual OR (r"| __truncated__
```

</details>

With one treated cohort, the simple aggregate is the average of ATT(3, 3) to ATT(3, 6), `1.375 + 0.5 * E[x | g = 3]`, which is 1.625 in the population. The known value uses the generated `x` of cohort 3; it is a target, not a fitted aggregate, and it involves no estimated cohort shares. `method` names the route; `selection_formula_inherited` reports whether `selection_xformla` was inherited from `xformla`, which these calls leave at `NULL`. The nuisance diagnostics of `fit_imputation` list the imputation cells and pooling groups, the counterfactual regression and the treatment logits, with no observation probability model. `summary()` adds a line describing either route.

## Refusals and limits

- Pure imputation fits no observation logit, so no such logit can cause a refusal, and it can run when every unit of a cell is observed. Its own models must still pass their support and rank checks, and the treatment logits and counterfactual regression keep theirs.
- The augmented correction fits every model and is refused when any of them fails, including an active observation cell with its own free intercept in which every unit, or none, is observed. Nothing falls back to another correction.
- A route's label does not make every cell available. A cell that the comparisons cannot support is reported as unavailable with its reason.
- Version 1.1.0.0 applies fixed treatment-propensity caps and inverse-observation floors without dropping observations; binding safeguards can alter the original target and robustness claims ([Probability safeguards](../methods/missing-outcomes.md#probability-safeguards)). A multiplier bootstrap, when requested, reuses the influence functions and refits nothing. Missing outcomes or changes are not randomly drawn: this is not multiple imputation and uses no combination by Rubin's rules. The multiplier bootstrap does draw random signs for inference.
- Standard errors come from complete influence functions with one row for every unit of `units`, including units observed once or never. They account for the estimated imputation, observation, treatment and counterfactual models, the response rates and the cohort shares used. With `compute_inffunc = FALSE` a fit returns point estimates only, and `aggte()` refuses it.

The argument rules are in [Options reference](../reference/options.md#handling-missing-outcomes), the counterfactual regression in [Treatment adjustment](../methods/treatment.md#change-imputation-uses-a-different-regression), and the reported uncertainty in [Results reference](../reference/results.md#uncertainty).

---

**Reading path, step 12 of 20** · Previous: [Choose an observation model](observation-models.md) · Next: [Parallel trends](../methods/parallel-trends.md)

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

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
