# Prepare your panel and unit table

Preparing longitudinal data carefully before calling `cdid::att_gt()` prevents common input errors and ensures that the fitted model targets your intended study population. This page explains how to format the long panel table (`data`) and the table of study units (`frame`), how `cdid` distinguishes units with repeated observations from units with one or no observed outcomes, how reference dates work when covariates or sampling weights change over time, and how to diagnose data and model issues.

For an introduction to how chained difference in differences combines outcome changes across observation windows, see [Estimate treatment effects from incomplete panels](../examples/incomplete-panels.md), [Get started](../get-started.md), and [Coming from `did`](../from-did.md).

## Core requirements for `data` and `frame`

In chained mode, `cdid::att_gt()` reads a long panel data frame (`data`) and an optional table of study units (`frame`, referred to after first introduction as the unit table).

For the simplest unweighted specification without covariates or clustering, the essential columns are a unit identifier, a finite numeric date, a numeric outcome with `NA` for missing values, and the unit's known first-treatment cohort (constant within unit, with `0` for never treated). Each unit/date pair must be unique. A balanced panel and an observed outcome at the reference date for every unit are not required. Valid columns, or two observed outcomes per unit, do not by themselves guarantee that the requested effects are estimable.

| Requirement | Column or argument | Rule in `cdid` |
| --- | --- | --- |
| Long panel rows | `data` | One row per unit and calendar date observation. |
| Unique unit and date pairs | `idname`, `tname` | `idname` is required in chained mode because observations are aligned by unit identity. Every `(idname, tname)` pair in `data` must be unique. |
| Stored ordered numeric dates | `tname` | Calendar dates must be finite numeric values. `cdid` sorts the distinct stored dates $t_1 < t_2 < \dots < t_T$ and measures duration in shared outcome regression and change imputation models by the stored calendar index gap $L_e = i_t - i_s$ between ordered calendar positions $i_s < i_t$, rather than raw subtraction $t - s$ of date labels. Under `pt_assumption = "all"`, `pt_start` must match an exact finite date in this stored calendar. |
| Constant treatment cohort | `gname` | Each unit's first treatment period must be constant across all of its rows in `data` and agree with `frame`. Units never treated must be coded as `0` (`gname == 0`). Treated cohorts use a finite positive treatment date in the same numeric scale as `tname`; it need not itself be a stored observation date. |
| Numeric outcomes, missing values kept as `NA` | `yname` | Observed values must be finite numeric outcomes. Keep missing outcomes as `NA`. Never substitute `0` or carry forward prior outcome values in place of a missing outcome. |
| Unit table of the study population | `frame` | Optional data frame with one row per study unit, covering every ID in `data` with the same ID type. Supply required cohort, weight, cluster and covariate attributes here when they are not available from panel rows; changing values are read from dated rows in `data`. At most one cluster variable beyond `idname` is accepted. |

### Never impute zero for missing outcomes

A frequent data preparation mistake is replacing missing outcomes (`NA`) with `0` before passing `data` to an estimator. In chained DiD, replacing a missing outcome $Y_{i,t}$ with `0` tells the software that unit $i$ was actually observed at date $t$ with an outcome equal to zero. Doing so corrupts both the observation indicators $S_{i,e}$ for every date pair touching date $t$ and the within unit changes $\Delta Y_{s,t} = Y_t - Y_s$. Always retain `NA` for unobserved outcomes.

## The unit table (`frame`) and unit observation roles

### Study population versus census

The unit table (`frame`) defines the study population $F$ whose treatment effects your analysis estimates. It lists the study units and need not be a census of a wider external population. Any claim that these units represent a wider population needs its own sampling justification. When `frame = NULL` (the default), the unique units present in `data` (including units with only one observed outcome or all `NA` outcomes) define $F$, and no wider external register is presumed. You may omit `frame` when every intended study unit appears in `data`. Units omitted entirely from `data` cannot be recovered from its rows. To retain such zero-observation units in the target population, supply a complete unit table through `frame`.

Required fields depend on the route and the unit's role. Panel units need known, constant cohorts. Baseline estimation can retain frame-only units with unknown cohort when the frame omits `gname`, but cohort-share-dependent aggregates are then unavailable. Conditional selection, pure imputation and augmented observation correction require known cohorts on every F row. Supply weights, cluster labels and covariates wherever that route reads them; the package does not infer missing attributes.

### Roles of units with zero, one, or repeated observations

Across the full unit table $F$ (of size $N = |F|$), study units fall into three nested sets based on how many nonmissing outcomes they have in `data`:
- **$A$ (units with repeated observations, at least two observed outcomes):** A unit must be observed at two or more dates to form a within unit outcome change $\Delta Y_{s,t} = Y_t - Y_s$. Units in $A$ contribute outcome changes to the specific date pairs $\{s, t\}$ where both endpoints $s$ and $t$ are observed. Having at least two observed outcomes makes a unit eligible to contribute changes, though it does not mean the unit is observed at every date pair.
- **$B \setminus A$ (units with one observed outcome, also called singletons):** A unit observed at only one date cannot form a within unit outcome change, because its single outcome level still contains the unobserved unit fixed effect $\alpha_i$.
- **$F \setminus B$ (units with no observed outcomes):** Units present in `frame` (or present in `data` with `NA` outcomes at all dates) supply no observed outcomes at all.

In chained mode, `cdid` **never drops** units with one or zero observed outcomes from the study population $F$:
1. **Population shares and rates:** Every unit remains in the sample size $N = |F|$ and the appropriate population mass and rate accounting. Cohort shares require known cohort membership; unknown frame-only cohorts make share-dependent aggregates unavailable.
2. **Models for handling missing outcomes:** Under `selection = "conditional_mtar"`, observation probability models (`selection_est_method = "ipw"` or `"dr"`) fitted on `selection_universe = "F"` evaluate every unit in $F$ (including nonobservers), while change imputation (`selection_est_method = "imputation"` or `"dr"`) fits actual change regressions on observed pairs and averages predicted changes across **all** rows of $F$ in each cohort.
3. **Uncertainty and clustering:** When inference is requested, `inffunc` has one row per F unit and carries every active nuisance and rate component, retaining clusters with no A units. Components absent from a chosen method are not fitted.

When `selection_model = "pair_logit"` is used with `selection_est_method = "ipw"`, `selection_universe` can be set to `"F"` (default, all study units), `"B"` (units with at least one observed outcome), or `"A"` (units with at least two observed outcomes). Choosing `"B"` or `"A"` fits the pair observation logit on that subset and requires a nonempty `selection_justification` string (plus an optional column name in `selection_rho` containing known normalized transport ratios when $U \in \{A, B\}$ is not conditionally independent of covariates given cohort). Fitting the observation logit on `"B"` or `"A"` does not change the target parameter away from the full study population $F$; it relies on the declared transport assumption, which fit diagnostics cannot verify. Under `selection_model = "independent_period_logit"`, `selection_est_method = "imputation"`, or `selection_est_method = "dr"`, `selection_universe` must be `"F"` with `selection_rho = NULL`. For details, see [Handling missing outcomes](../methods/missing-outcomes.md) and [Option reference](../reference/options.md).

## Panel mode selection and `allow_unbalanced_panel`

`cdid::att_gt()` supports both chained estimation (`mode = "chained"`) and ordinary delegation to `did::att_gt()` (`mode = "ordinary"`). Both modes can analyze balanced or incomplete panels according to their respective options and assumptions.

- **Default mode (`mode = NULL`):** `cdid` evaluates `panel` first. When `panel = TRUE` (the default), `mode = NULL` selects `"chained"`. When `panel = FALSE` (genuine repeated cross sections), `mode = NULL` selects `"ordinary"` and forwards the call to `did::att_gt()`. Requesting `mode = "chained"` with `panel = FALSE` raises a typed error (`"cdid_chained_panel"`).
- **`allow_unbalanced_panel` in chained mode versus ordinary mode:**
  - In **chained mode**, both `allow_unbalanced_panel = NULL` (the default) and `allow_unbalanced_panel = TRUE` retain all supplied units in $F$ and use all units with repeated observations ($A$), whether the panel is balanced or incomplete.
  - Passing an explicit `allow_unbalanced_panel = FALSE` in chained mode is refused with an error (`"cdid_chained_balanced"`) that points to `mode = "ordinary"`. `cdid` never silently drops incomplete units to switch the target parameter to complete cases.
  - In **ordinary mode** (`mode = "ordinary"`), `allow_unbalanced_panel = NULL` resolves to `did`'s default `FALSE` (balancing the panel to complete cases), whereas `allow_unbalanced_panel = TRUE` forwards `TRUE` to `did::att_gt()`. See [Coming from `did`](../from-did.md) for a direct comparison.

## Static versus changing covariates and sampling weights

### Static covariates and weights

When covariates (`xformla`, `selection_xformla`) and sampling weights (`weightsname`) are constant within each unit across all rows in `data` and agree with `frame`, `fix_weights = NULL`, `fix_weights = "first_period"`, and `fix_weights = "base_period"` produce identical estimates. Required covariate coverage still depends on the route: for example, baseline IPW with a static design requires that design for all study units, while baseline OR/DR require it on repeat observers.

### Reference dates for covariates that change within a unit

When covariates in `xformla` or `selection_xformla` change across dates within a unit (or when a formula uses `tname`), `cdid` treats them as reference date covariates:
- **Posttreatment comparisons:** Every posttreatment comparison for target cohort $g$ (including intermediate bridges and comparison units) evaluates each unit's covariates at cohort $g$'s baseline reference date

$$b_g = \max\{d \in \{t_1,\ldots,t_T\}: d + \text{anticipation} < g\}$$

This cutoff compares numeric date values, using `anticipation` in that scale. It is distinct from the stored-index gap $L_e$ used by outcome regression and imputation. If there is no admissible reference date, the relevant target has no reference support.

- **Auxiliary preperiod comparisons:** Preperiod reports evaluate covariates at the earlier endpoint of the date pair.
- **Read from actual `(id, date)` rows, including rows with `y = NA`:** Reference date covariate values come directly from the `(idname, tname)` rows in `data`; a static value in `frame` cannot override changing values in `data`. Crucially, a unit's covariate row at `(id, b_g)` may be needed even when the unit's outcome `y` at date $b_g$ is missing (`NA`). For instance, if a unit in cohort $g = 3$ (with $b_3 = 2$ under `anticipation = 0`) has observed outcomes only at $t = 3$ and $t = 4$, its covariate values at `(id, time = 2)` are still required to form the reference design for cohort $g = 3$ comparisons. Keep the `(id, time = 2)` row in `data` with its covariate values and `y = NA`.
- **No blanket full history or baseline outcome requirement:** You do not need observed outcomes at $b_g$ for every unit, nor do you need covariate rows at dates that no active comparison or model reads. Covariate values are required only at the `(id, date)` anchors read by the requested estimator: on observed pairs for baseline treatment adjustment (`selection = "none"`), on every unit in $F$ of each cohort for change imputation and doubly robust adjustment for missing outcomes (`selection_est_method = "imputation"` or `"dr"`), and on the selection universe for observation probability models. If a required reference value is missing, `cdid` refuses the call and reports the unit ID, date, variable, and support context rather than carrying values forward across time.

### Changing sampling weights and `fix_weights`

When sampling weights change within a unit, `fix_weights = NULL` is refused in chained mode. The explicit `"varying"` option is refused in chained mode even for static weights. You must specify one of two explicit reference weight rules:
- `fix_weights = "first_period"`: Every comparison and rate denominator uses each unit's weight $w_{t_1}$ at the first stored calendar date $t_1$.
- `fix_weights = "base_period"`: Every comparison in target cohort $g$'s context (controls included) uses each unit's weight at cohort $g$'s reference date $b_g$.

**Aggregation shares (`aggte`) under reference inputs:** When aggregating reference input fits with `cdid::aggte()`, treated cohorts are weighted by their first calendar date masses $w_{t_1}$, and uncertainty in those cohort shares is included in the influence functions. Units never treated (`g == 0`) do not need $w_{t_1}$ for cohort aggregation shares. A share-dependent aggregate is unavailable when a required treated-cohort first-date mass is incomplete, with the affected unit IDs recorded; one-cohort bins and a one-bin group overall summary cancel shares exactly and do not require that mass. Note that under `fix_weights = "base_period"`, each cohort cell uses its own reference weight $w_{b_g}$, so an aggregate across cohorts is a declared mixture of cohort specific weighted effects rather than the effect under one common weight distribution.

## Outcome formula (`xformla`) versus observation formula (`selection_xformla`)

`cdid::att_gt()` separates the covariate formula for treatment and outcome change models from the covariate formula for the observation probability model:
- `xformla` specifies the one sided covariate design for treatment adjustment (`est_method = "dr"`, `"reg"`, or `"ipw"`, see [Treatment adjustment](../methods/treatment.md)) and, when active, change imputation (`selection_est_method = "imputation"` or `"dr"`). When covariates are used, `xformla` requires a nonempty design (such as `~1`, `~x1 + x2`, or `~0 + x1`).
- `selection_xformla` specifies the one sided covariate design for the observation probability model under `selection = "conditional_mtar"` (`selection_model = "pair_logit"` or `"independent_period_logit"`). When `selection_xformla = NULL` (the default), it inherits `xformla`.
- Under pure change imputation (`selection_est_method = "imputation"`), no observation probability model is fitted, so `selection_xformla`, `selection_model`, and `selection_pooling` are inactive and never evaluated.
- Under doubly robust adjustment for missing outcomes (`selection_est_method = "dr"`, which requires `est_method = "dr"`), `xformla` builds the imputation, counterfactual outcome regression, and treatment propensity designs on all rows of $F$, while `selection_xformla` builds the observation probability design on all rows of $F$. Variables included only in `selection_xformla` do not enter the outcome or treatment models; however, conditional mean MTAR must hold given the combined covariate set, and if variables present only in `selection_xformla` shift the outcome change mean, an `xformla` only imputer is not sufficient on its own. See [Handling missing outcomes](../methods/missing-outcomes.md) and [Option reference](../reference/options.md).

## Diagnosing separation, collinearity, and covariance issues

First inspect who is observed, at which dates, and in the particular pairs your analysis needs. [Investigate missingness and data structure](investigate-missingness.md) gives executable base R checks and shows how to inspect the fitted support, model, parallel-trends and inference records. These checks assess observable structure. They cannot establish parallel trends or the representativeness of missing changes.

When preparing covariates and inspecting fit messages, keep three distinct mechanisms separate:

### 1. Separation and boundary probabilities in logit models

Treatment propensity logits (`est_method = "ipw"` or `"dr"`) and active observation probability logits (`selection = "conditional_mtar"` with `selection_est_method = "ipw"` or `"dr"`) require interior probabilities and appropriate support. The default `selection = "none"` fits no observation logit.
- If a covariate or factor level perfectly separates treated units from comparison units in an active risk set, or separates observed units from missing units in an active observation cell, no finite interior maximum likelihood root exists and `cdid` refuses the fit with a typed error.
- In particular, if every unit in an active observation cell is observed or every unit is unobserved and that cell has a free intercept, the logit has no finite interior root (pooling slopes does not remove cell intercepts) and `cdid` refuses the separated fit. Observing all units in a sample cell does not establish a known design probability of 1 when a logit observation model is requested.

### 2. Collinearity preprocessing versus unsupported covariate predictions

In built in outcome regression (`est_method = "reg"` or `"dr"`) and change imputation (`selection_est_method = "imputation"` or `"dr"`), `cdid` uses standard pivoted QR decomposition (tolerance `1e-7`) to propose aliased (collinear) columns on the fitting design:
- **Prediction row reconstruction check:** A proposed alias is accepted only after reconstruction at arithmetic roundoff (`64 * .Machine$double.eps` times the original column count) on the required prediction rows. Counterfactual regression checks its control/treated supports; each imputation pooling group checks its own observed training rows and full-F cohort prediction rows.
- **Unsupported treated covariate directions:** If a predictor is constant or collinear among the control fitting rows (for example, a factor level that appears only among treated units) so that it varies on treated prediction rows in a null direction of the control design, the retained control columns cannot reconstruct the treated prediction rows. `cdid` refuses that unsupported prediction direction rather than silently setting its slope to zero. Near-collinear columns that fail exact-alias verification are retained and remain subject to the regression rank/conditioning check.
- **Pooled control identification across blocks:** Because the built in counterfactual outcome regression fits one shared slope vector $\beta$ on the duration scaled covariates $L_e z_i$ together with local block intercepts $\alpha_b$, a predictor that is locally rank deficient within one comparison block remains supported when the pooled control design across active blocks identifies the shared slope and reconstructs all prediction rows.
- Dropping a truly redundant aliased column is a numerical reparameterization of collinear predictors; it is not an identifying assumption that an underlying structural effect is zero, and it does not relax general ATT identification requirements. Treatment and observation logits do not drop aliased columns and refuse directly on rank deficiency.

### 3. Covariance regularization versus ATT identification

In feasible two step GMM (`gmm_weighting = "optimal"`, the default), comparison moments can have singular or nearly singular covariance. Exact structural equivalence reductions come first. The default `gmm_covariance = "tikhonov"` handles the remaining covariance problem using the recorded conditioning and rank rules, including `control$condition_trigger` (`1e-10` by default). Exactly zero covariance follows the existing full-catalog least-squares convention with zero sandwich variance.

Covariance regularization only stabilizes the weighting matrix used to combine already identified comparison moments; it **cannot** create ATT identification when date links, cohort support, or identifying assumptions are missing. Furthermore, analytical and multiplier bootstrap standard errors use the estimated sandwich covariance built from the complete influence functions, never the inverse of the regularized weighting matrix (see [Results and uncertainty](../reference/results.md)).

## Resource and workspace controls (`control`)

Chained mode accepts an optional named `control` list for internal numeric workspace and preflight work limits (see [Option reference](../reference/options.md)):
- **Configurable numeric workspace (`numeric_workspace_bytes`):** Defaults to `1073741824` (1 GiB) of counted numeric array payload, with automatic memory aware chunking controlled by `chunk_rows`.
- **Opt in estimated work cap (`counted_work_limit`):** Defaults to `Inf` (no estimated work cap). You can opt in by setting `counted_work_limit` to a finite positive whole number (such as `1e12`), which refuses a stage before it runs if its conservative preflight work estimate exceeds your budget.
- **Parallel multiplier bootstrap workers (`workers`):** Defaults to `1`. To run multiplier bootstrap draws in parallel with `pl = TRUE, cores = k`, you must also declare `control = list(workers = k)`, and all worker buffers count against `numeric_workspace_bytes`.
- **Scope of resource controls:** `numeric_workspace_bytes` and `counted_work_limit` track internal numeric array allocations and conservative preflight operation counts; they do not cap total R process resident memory (RSS), R heap usage, or wall clock CPU time, and `cdid` does not automatically detect machine RAM.

## Concise base R data checklist and examples

> [!IMPORTANT]
The checks below are templates for your `panel_data` and `units` tables. The [observation-model example](../examples/observation-models.md) shows a complete simulation, fitted results and output; the column names below illustrate how to prepare your own data.

### 1. Verify long panel structure, unit table, and missing outcome coding

```r
# Data checks for your panel and unit table.
# Suppose `panel_data` has columns: id, time, g, y, x1, z_sel, w
# and `units` (the unit table `frame`) has one row per study unit.

# 1. Check unique unit IDs in `units` and unique (id, time) pairs in `panel_data`.
stopifnot(
  identical(class(units$id), class(panel_data$id)),
  !anyNA(units$id),
  !anyDuplicated(units$id),
  all(panel_data$id %in% units$id),
  !anyNA(panel_data$id),
  !anyNA(panel_data$time),
  !anyDuplicated(panel_data[c("id", "time")])
)

# 2. Check that calendar dates are finite numeric values and inspect the sorted calendar.
stopifnot(is.numeric(panel_data$time), all(is.finite(panel_data$time)))
calendar_dates <- sort(unique(panel_data$time))

# 3. Check that cohort `g` is finite, constant within each unit, agrees with `units`,
# and is nonnegative (0 means never treated). This checklist assumes known cohorts
# for every study unit; that is convenient here, not a universal baseline requirement.
unit_g_lookup <- setNames(units$g, units$id)
stopifnot(
  is.numeric(units$g),
  all(is.finite(units$g)),
  all(units$g >= 0),
  all(panel_data$g == unit_g_lookup[as.character(panel_data$id)])
)

# 4. Count observed (not NA) outcomes per unit across the full unit table `units`.
# Counts do not reveal whether someone previously replaced missing outcomes with 0.
obs_counts <- tabulate(
  match(panel_data$id[!is.na(panel_data$y)], units$id),
  nbins = nrow(units)
)
n_total_F    <- nrow(units)          # Full study population F
n_at_least_1 <- sum(obs_counts >= 1) # Set B: at least one observed outcome
n_repeated_A <- sum(obs_counts >= 2) # Set A: units with repeated observations
n_singletons <- sum(obs_counts == 1) # Set B \ A: one observed outcome
n_zero_obs   <- sum(obs_counts == 0) # Set F \ B: no observed outcomes
```

### 2. Check reference date rows when covariates or weights change over time

For this template, `x1` and `w` vary by date and are supplied only in `panel_data`; keep their dated rows even when `y` is missing. Do not also place a conflicting static copy in `frame`. The inspection below does not assert that every F unit needs a row at every cohort baseline.

```r
# Check and call for changing covariates and weights.
# This example uses zero anticipation: b_g is the last stored date before g.
# For nonzero anticipation the cutoff is d + anticipation < g on numeric dates.
anticipation <- 0
treated_cohorts <- sort(unique(units$g[units$g > 0]))
base_dates <- vapply(treated_cohorts, function(g_val) {
  eligible <- calendar_dates[calendar_dates < g_val]
  if (length(eligible) == 0L) NA_real_ else max(eligible)
}, numeric(1))
names(base_dates) <- as.character(treated_cohorts)
base_dates

# Inspect rows at the reference dates, retaining y = NA where appropriate.
# This is an inspection aid, not an all-F baseline-row requirement. The package
# checks exact required IDs/dates by route, pair support, and recipe, including
# auxiliary preperiod recipes; it reports missing required values explicitly.
reference_rows <- panel_data[panel_data$time %in% base_dates, ]
head(reference_rows[c("id", "time", "g", "y", "x1", "w")])

# Keep varying x1 and w in panel_data, not as dateless unit-table attributes.
frame_ref <- units[c("id", "g")]

# Fit chained DiD with reference date covariates and explicit base period weights.
fit_ref <- cdid::att_gt(
  yname       = "y",
  tname       = "time",
  idname      = "id",
  gname       = "g",
  xformla     = ~x1,
  data        = panel_data,
  weightsname = "w",
  fix_weights = "base_period",
  est_method  = "dr",
  frame       = frame_ref,
  bstrap      = FALSE,
  cband       = FALSE
)
```

### 3. Separate outcome and observation formulas with optional resource controls

This is an independent, unweighted template with static covariates. It assumes `x1` and `z_sel` are available consistently in the panel or unit table on the route's required supports and the chosen observation logits have finite interior roots. Distinct predictor lists require conditional-MTAR and treatment-model restrictions sufficient for their scientific conditioning sets; adding `z_sel` is not a substitute for stating those assumptions.

```r
# Call template with separate xformla and selection_xformla and an optional work cap.
fit_sel <- cdid::att_gt(
  yname                = "y",
  tname                = "time",
  idname               = "id",
  gname                = "g",
  xformla              = ~x1,
  selection            = "conditional_mtar",
  selection_est_method = "ipw",
  selection_model      = "pair_logit",
  selection_pooling    = "within_cohort",
  selection_xformla    = ~x1 + z_sel,
  selection_universe   = "F",
  data                 = panel_data,
  frame                = units,
  est_method           = "dr",
  bstrap               = FALSE,
  cband                = FALSE,
  control              = list(
    counted_work_limit      = 1e12,
    numeric_workspace_bytes = 1073741824
  )
)
```

## Where to go next

- To see how chained DiD links pair comparisons across rotating, overlapping, and gapped panels, read [Estimate treatment effects from incomplete panels](../examples/incomplete-panels.md).
- To choose between observation weighting, change imputation, and doubly robust adjustment for missing outcomes, see [Handling missing outcomes](../methods/missing-outcomes.md).
- To compare treatment adjustment models (`est_method = "dr"`, `"reg"`, `"ipw"`), see [Treatment adjustment](../methods/treatment.md).
- To compare `PT-post` and `PT-all` (`pt_assumption` and `pt_start`) and reference period reporting (`base_period`), see [Parallel trends assumptions](../methods/parallel-trends.md).
- For complete parameter definitions and output structures, see [Option reference](../reference/options.md) and [Results and uncertainty](../reference/results.md).

---

**Reading path, step 6 of 20** · Previous: [How to compare precision](../examples/precision.md) · Next: [Treatment adjustment](../methods/treatment.md)

<details>
<summary>All steps of the reading path</summary>

| Step | Page |
| --- | --- |
| 1 | [Get started](../get-started.md) |
| 2 | [Moving from did](../from-did.md) |
| 3 | [Incomplete panels](../examples/incomplete-panels.md) |
| 4 | [Identify effects when the observed units change composition](../examples/identification.md) |
| 5 | [How to compare precision](../examples/precision.md) |
| 6 | **Prepare your panel and unit table** (this page) |
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
| 20 | [Diagnose errors and unavailable results](troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
