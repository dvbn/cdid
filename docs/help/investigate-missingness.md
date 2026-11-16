# Investigate missingness and data structure

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

Before choosing an adjustment, establish which study units are in your data, when their outcomes are observed and which changes they supply. Then inspect the fitted support and model records. These checks reveal observation patterns and model problems. They cannot prove that missing changes are representative in mean or that untreated trends are parallel.

## What data do I need?

For a simple chained specification, supply long data with these four columns:

| Column | Required information |
| --- | --- |
| Unit ID | A nonmissing identifier that links repeated observations of the same unit |
| Time | A finite numeric date, on a common scale across units |
| Outcome | A numeric outcome, finite when observed and `NA` when missing |
| First treatment cohort | A known, constant first-treatment date for each unit, with `0` for never treated, on the same numeric scale as time |

Each `(ID, time)` pair must be unique. The treatment date itself need not appear as an observed row. A balanced panel and an observed baseline outcome for every unit are not required. Keep a date in the stored calendar through rows with `y = NA` when you intend it to be part of that calendar. An entirely omitted date is not recovered automatically.

The optional `frame` is a unit table with one row per intended study unit. You can omit it when every intended unit appears somewhere in `data`, including units whose outcomes are all `NA`. A unit absent from all rows cannot be recovered from the observed panel: a complete unit table is needed to include such zero-observation units. The package cannot infer their number or their attributes from the observed units.

Additional requirements follow from your request. Covariates, positive sampling weights and cluster labels must be supplied on the units and dates read by the active models. If covariates or weights change over time, actual dated reference rows may be needed even where `y = NA`. Values are not carried forward. [Prepare your panel and unit table](prepare-data.md) explains the reference rules, model-specific coverage and the limited baseline handling of unknown frame-only cohorts. Supplying known cohorts for all study units makes population counts and cohort-share aggregation interpretable and is required by conditional observation adjustments.

Correct columns and two observations per unit are not an identification guarantee. Each needed comparison also requires cohort and eligible-control support, overlap, valid observation reconstruction, parallel trends and the relevant nuisance-model conditions. The target must be estimable from the retained catalog under those assumptions.

## A small observation-pattern example

The next blocks are self-contained base R code. They describe twelve study units in cohorts 0, 3 and 5, including units observed once or never. The artificial outcomes only make the input concrete. This tiny table illustrates checks and is not a recommended estimation sample.

```r
unit_table <- data.frame(
  id = 1:12,
  first_treat = rep(c(0, 3, 5), each = 4),
  x = rep(c(-1, 0, 1, 2), times = 3)
)

# Within each cohort: one unit with no observed outcome, one singleton,
# one unit observed at {1, 2}, and one observed at {2, 4, 5}.
panel_check <- data.frame(
  id = c(2, 3, 3, 4, 4, 4, 6, 7, 7, 8, 8, 8,
         10, 11, 11, 12, 12, 12),
  period = rep(c(2, 1, 2, 2, 4, 5), times = 3)
)
pos <- match(panel_check$id, unit_table$id)
panel_check$first_treat <- unit_table$first_treat[pos]
panel_check$x <- unit_table$x[pos]
panel_check$y <- panel_check$id / 10 + panel_check$period / 5 +
  as.numeric(panel_check$first_treat > 0 &
               panel_check$period >= panel_check$first_treat)

# Retain date 3 in the stored calendar, without declaring an observed outcome.
# Units 5 and 9 have no panel rows and are known only through the unit table.
panel_check <- rbind(panel_check, data.frame(
  id = 1L, period = 3, first_treat = 0, x = -1, y = NA_real_
))
panel_check <- panel_check[order(panel_check$id, panel_check$period), ]
rownames(panel_check) <- NULL
```

### Validate identity, dates and cohorts

For your own data, replace `panel_check` and `unit_table` with your long panel and known unit table using these column names, or adapt the names in the checks. If no separate unit table is needed, first form it from unique unit attributes in the panel and verify that each unit has only one cohort. Do not construct a supposedly complete unit table from observed IDs when the target also includes units absent from the panel.

```r
stopifnot(
  !anyNA(unit_table$id), !anyDuplicated(unit_table$id),
  !anyNA(panel_check$id),
  !anyDuplicated(panel_check[c("id", "period")]),
  is.numeric(panel_check$period), all(is.finite(panel_check$period)),
  is.numeric(panel_check$y),
  all(is.na(panel_check$y) | is.finite(panel_check$y)),
  is.numeric(unit_table$first_treat),
  all(is.finite(unit_table$first_treat)),
  all(unit_table$first_treat >= 0),
  is.numeric(panel_check$first_treat),
  all(is.finite(panel_check$first_treat))
)
pos <- match(panel_check$id, unit_table$id)
stopifnot(!anyNA(pos))
stopifnot(all(panel_check$first_treat == unit_table$first_treat[pos]))
sort(unique(panel_check$period))
```

Output:

```text
[1] 1 2 3 4 5
```

The last check verifies both constancy of each panel unit's cohort and agreement with the unit table. It does not verify the substantive accuracy of the recorded treatment dates. Check your source records as well.

### Count outcomes using all study units as the denominator

```r
observed <- !is.na(panel_check$y)
n_outcomes <- tabulate(
  match(panel_check$id[observed], unit_table$id),
  nbins = nrow(unit_table)
)
cohorts <- sort(unique(unit_table$first_treat))

coverage <- do.call(rbind, lapply(cohorts, function(g) {
  take <- unit_table$first_treat == g
  z <- n_outcomes[take]
  data.frame(cohort = g, n_study_units = sum(take),
             n_zero = sum(z == 0), n_singleton = sum(z == 1),
             n_repeated = sum(z >= 2),
             share_repeated = mean(z >= 2))
}))
coverage

date_coverage <- do.call(rbind, lapply(sort(unique(panel_check$period)), function(t) {
  seen <- unit_table$id %in%
    panel_check$id[panel_check$period == t & observed]
  do.call(rbind, lapply(cohorts, function(g) {
    take <- unit_table$first_treat == g
    data.frame(cohort = g, period = t, n_study_units = sum(take),
               n_observed = sum(seen[take]), share_observed = mean(seen[take]))
  }))
}))
date_coverage
```

Output:

```text
cohort n_study_units n_zero n_singleton n_repeated share_repeated
1      0        4      1           1          2            0.5
2      3        4      1           1          2            0.5
3      5        4      1           1          2            0.5
   cohort period n_study_units n_observed share_observed
1       0      1        4          1           0.25
2       3      1        4          1           0.25
3       5      1        4          1           0.25
4       0      2        4          3           0.75
5       3      2        4          3           0.75
6       5      2        4          3           0.75
7       0      3        4          0           0.00
8       3      3        4          0           0.00
9       5      3        4          0           0.00
10      0      4        4          1           0.25
11      3      4        4          1           0.25
12      5      4        4          1           0.25
13      0      5        4          1           0.25
14      3      5        4          1           0.25
15      5      5        4          1           0.25
```

Zeros and singletons stay in the denominators. Units with two or more observed outcomes supply changes only for the pairs at which both endpoints are observed. A high repeat-observation rate can coexist with no support for a particular pair. If you use sampling weights, distinguish these raw counts from the corresponding weighted population masses under your declared reference-weight rule.

### Investigate selected pairs

```r
pair_counts <- function(s, t) {
  stopifnot(s < t)
  at_s <- unit_table$id %in%
    panel_check$id[panel_check$period == s & observed]
  at_t <- unit_table$id %in%
    panel_check$id[panel_check$period == t & observed]
  do.call(rbind, lapply(cohorts, function(g) {
    take <- unit_table$first_treat == g
    data.frame(cohort = g, start = s, end = t, n_study_units = sum(take),
               n_start = sum(at_s[take]), n_end = sum(at_t[take]),
               n_both = sum(at_s[take] & at_t[take]),
               share_both = mean(at_s[take] & at_t[take]))
  }))
}
selected_pairs <- rbind(pair_counts(1, 2), pair_counts(2, 4), pair_counts(1, 4))
selected_pairs
```

Output:

```text
cohort start end n_study_units n_start n_end n_both share_both
1      0     1   2        4       1     3      1       0.25
2      3     1   2        4       1     3      1       0.25
3      5     1   2        4       1     3      1       0.25
4      0     2   4        4       3     1      1       0.25
5      3     2   4        4       3     1      1       0.25
6      5     2   4        4       3     1      1       0.25
7      0     1   4        4       1     1      0       0.00
8      3     1   4        4       1     1      0       0.00
9      5     1   4        4       1     1      0       0.00
```

This table asks who can supply each change. It does not reproduce the package catalog or select a valid comparison group. For example, for cohort 3 with zero anticipation and reference date 2, cohort 5 remains untreated through date 4 and can be eligible for a not-yet-treated comparison ending at 4. It cannot supply untreated controls at date 5. Never-treated controls are coded 0. Actual catalog retention also depends on the comparison's reference, control pool, pair support and the requested models.

Start with pairs relevant to your scientific targets and suspected gaps. On large panels, keep the long data and count selected pairs or dates, as above. Do not build a dense units-by-dates matrix or every possible pair merely to explore missingness. The small date-by-cohort output table is distinct from a full unit-calendar expansion. After fitting, use the package's support records to investigate which retained comparisons support each requested target.

### Compare observed covariate composition

```r
at_2 <- unit_table$id %in%
  panel_check$id[panel_check$period == 2 & observed]
at_4 <- unit_table$id %in%
  panel_check$id[panel_check$period == 4 & observed]
composition <- do.call(rbind, lapply(cohorts, function(g) {
  take <- unit_table$first_treat == g
  pair <- take & at_2 & at_4
  data.frame(cohort = g, mean_x_study = mean(unit_table$x[take]),
             n_pair = sum(pair),
             mean_x_pair = if (any(pair)) mean(unit_table$x[pair]) else NA_real_)
}))
composition
```

Output:

```text
cohort mean_x_study n_pair mean_x_pair
1      0           0.5      1           2
2      3           0.5      1           2
3      5           0.5      1           2
```

Here `x` is known for the study population and is constant within unit. In your application, inspect relevant covariates at the reference date actually read by the model, with the intended weights where needed. Compare distributions across all study units and observed pairs, not only a pooled mean. Look for absent factor levels, limited overlap and cells with no observed or no missing units. These can expose implausible specifications or a logit separation problem. Similar observed covariate distributions still do not establish representativeness of the unobserved changes.

## Read the diagnostics from a fitted result

If the fit warns about missing reference covariates, review
`fit$chained$reference_eligibility$targets` before interpreting the effects.
This table identifies the cohorts and reference dates to inspect in your data.
Use a covariate-only row for any independently known reference value.
[Reference covariate diagnostics](../reference/results.md#reference-covariate-exclusions).

Run the [rotating-panel worked example](../examples/incomplete-panels.md#worked-example-rotating-consecutive-panel-in-base-r) first to create `fit_rot`. The next block inspects that result without fitting anything again. It shows the study units and the support for each reported effect.

```r
fit <- fit_rot
summary(fit)
cdid::glance(fit)
cells <- cdid::tidy(fit)
cells[, c("target.id", "role", "estimate", "std.error", "available",
          "support.status", "support.reason")]

head(fit$chained$frame)
with(fit$chained$frame, table(cohort, eligible, useNA = "ifany"))
```

<details>
<summary>Console output</summary>

```text
cdid chained difference-in-differences (chained IPW-GMM): cdid_MP
  cdid 1.1.0.0 result built on did 2.5.1 (tested: did 2.5.1); did's MP methods tabulate the estimates below.
  Chained estimators: Bellego C., Benatia D. and V. Dortet-Bernadet (2025). "The Chained Difference-in-Differences." Journal of Econometrics 248, 105783. doi:10.1016/j.jeconom.2024.105783.
  Ordinary estimators: Callaway B. and P. H. C. Sant'Anna (2021). "Difference-in-Differences with Multiple Time Periods." Journal of Econometrics 225(2), 200-230. doi:10.1016/j.jeconom.2020.12.001.

Call:
cdid::att_gt(yname = "y", tname = "time", idname = "id", gname = "g", 
    xformla = NULL, data = panel_data, panel = TRUE, allow_unbalanced_panel = TRUE, 
    control_group = "nevertreated", est_method = "ipw", selection = "none", 
    pt_assumption = "post", bstrap = FALSE, cband = FALSE, frame = units)

Group-Time Average Treatment Effects:
 Group Time ATT(g,t) Std. Error [95% Pointwise  Conf. Band]  
     3    2   0.0469     0.0823         -0.1143      0.2081  
     3    3   0.8946     0.0922          0.7139      1.0753 *
     3    4   1.2016     0.1247          0.9573      1.4460 *
     3    5   1.3501     0.1573          1.0419      1.6583 *
---
Signif. codes: `*' confidence band does not cover 0

P-value for pre-test of parallel trends assumption:  0.568617658368036
Control Group:  Never Treated,  Anticipation Periods:  0
Estimation Method:  Inverse Probability Weighting

Chained support: 4 of 4 requested cells available.
  Units: 1000 in the unit table, 0 observed at all 5 dates, 1000 with at least two observed outcomes.
  Observation correction: none (selection = "none").
  Parallel trends: PT-post; base_period = "varying".
  GMM weighting: optimal (two-step); covariance regularization: not applied (moment covariance well conditioned).
  Inference: analytic standard errors, unit-level; pointwise intervals.
  Parallel-trends pre-test: W = 0.325001, df = 1, p = 0.5686 (Wald test of the pretrend cell).
  nobs ngroup ntime control.group est.method       estimator n.requested n.available
1 1000      1     5  nevertreated        ipw chained IPW-GMM           4           4
                          band.family inference clustered
1 pointwise_normal_over_supported_ids  analytic     FALSE
     target.id role   estimate  std.error available            support.status
1 pre(g=3,t=2)  pre 0.04689223 0.08225430      TRUE supported_direct_pretrend
2 ATT(g=3,t=3) post 0.89461117 0.09218505      TRUE                 estimable
3 ATT(g=3,t=4) post 1.20164289 0.12468476      TRUE                 estimable
4 ATT(g=3,t=5) post 1.35007838 0.15725715      TRUE                 estimable
                                                                        support.reason
1                                 direct auxiliary comparison 1->2 with controls {Inf}
2 level anchored to the cohort reference through retained comparison rows (prop:graph)
3 level anchored to the cohort reference through retained comparison rows (prop:graph)
4 level anchored to the cohort reference through retained comparison rows (prop:graph)
  id cohort weight cluster eligible n_obs frame_only
1  1    Inf      1      NA     TRUE     2      FALSE
2  2    Inf      1      NA     TRUE     2      FALSE
3  3      3      1      NA     TRUE     2      FALSE
4  4      3      1      NA     TRUE     2      FALSE
5  5    Inf      1      NA     TRUE     2      FALSE
6  6      3      1      NA     TRUE     2      FALSE
      eligible
cohort TRUE
   3    494
   Inf  506
```

</details>

| Record | What to investigate |
| --- | --- |
| `frame` | The full unit table, with `id`, `cohort`, `weight`, `cluster` and `eligible`. Its cohort representation uses `Inf` for never treated. `eligible` marks repeated-outcome units, not units guaranteed to support every pair. |
| `support_status`, `support_reason`, `support`, `support_ledger` | Which requested cells the retained catalog recovers and the reasons for unavailable cells. `tidy()` gives a convenient cell-level table. |
| `selection` and `selection$models` | The observation adjustment and its active model counts. Counts such as `n_selection_fits` are nested in `models` when that route fits them. They are not universal top-level fields. |
| `nuisance_diagnostics` and `outcome_design` | Active treatment, observation, counterfactual and imputation designs and fit diagnostics. Their detailed nesting depends on the route. Inspect the relevant record rather than assuming every method fits the same models. |
| `pt` | The maintained parallel-trends choice. Under PT-all, `pt$by_cohort` gives reference dates and the represented added restrictions, including no-extra cases. |
| `covariance_handling` | The GMM weighting and whether covariance regularization was used, with its reason and conditioning record. |
| `inference` and `pretest` | What uncertainty and pretrend tests were actually computed, with draw, cluster and availability records. |

The example uses `selection = "none"`, so it fits no observation model. To inspect fitted observation-model counts, first run the [missing-outcomes example](../methods/missing-outcomes.md#example-observation-that-depends-on-a-covariate), which creates `fit_weighting`, then inspect:

```r
str(fit_weighting$chained$selection$models, max.level = 1)
```

<details>
<summary>Console output</summary>

```text
List of 10
 $ n_selection_fits      : int 12
 $ n_selection_cells     : int 12
 $ n_selection_blocks    : int 12
 $ n_selection_parameters: int 24
 $ n_regression_fits     : int 1
 $ n_regression_blocks   : int 15
 $ n_local_intercepts    : int 15
 $ n_shared_slopes       : int 1
 $ n_supports            : int 30
 $ n_treatment_fits      : int 15
```

</details>

With inference computed, the influence-function matrix retains one row for every study unit, including units with no observed change, and carries the uncertainty of the active models and rates. Do not subset it to repeat observers when reconstructing uncertainty. With `compute_inffunc = FALSE`, no influence-function matrix or standard errors are returned.

## Interpret findings without overstating them

| Finding | What it says and what it does not establish |
| --- | --- |
| Dropout, intermittent gaps or changing observed covariates | Observable features that help assess the design and choose scientifically justified models. They do not prove mean missing-at-random or parallel trends. |
| Observed-pair units differ in outcome levels from the study population | This can matter for level comparisons. Chained DiD allows observation to depend on a unit fixed effect when the relevant mean changes remain representative. Differencing alone does not make this assumption true. |
| A cell is unavailable in `tidy(fit)` | With the chosen comparisons and assumptions, this fit cannot recover the effect. |
| A required logit is separated, or a model lacks rank or prediction support | Check which covariate levels appear in the fitting and prediction samples. In a logit, check whether predictors perfectly separate the two outcomes. |
| Moment covariance is singular or regularized | Comparisons share sampling variation. Covariance regularization stabilizes their weighting and cannot create target identification. See the two [analytical examples](../methods/gmm.md#covariance-singularity-versus-parameter-identification). |
| Standard errors or a pretrend p value are reported | The computation was available. Valid inference still needs its large-sample assumptions, and few independent clusters can make the approximation unreliable. A nonrejection does not establish parallel trends. |

Within a cohort, mean representativeness equates the population expected outcome change for units observed in a pair with the cohort's population expected change over that pair. It requires neither identical unit-level changes nor exact sample-mean equality. Missingness can depend on unit fixed effects while this holds. If those fixed effects predict trends or treatment-effect heterogeneity, the condition can fail. The conditional version additionally needs a scientifically sufficient conditioning set, appropriate models and the positivity or prediction-support conditions of the chosen adjustment. [Handling missing outcomes](../methods/missing-outcomes.md#the-assumption-about-missing-outcomes) states these conditions.

Continue with [Prepare your panel and unit table](prepare-data.md), [Incomplete panels](../examples/incomplete-panels.md) or [Results reference](../reference/results.md). For refused fits, see [Diagnose errors and unavailable results](troubleshooting.md).

[Guide home](../README.md) · [Package overview](../../README.md)
