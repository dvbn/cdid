# Add PT-all restrictions

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

In chained difference in differences (`cdid`, extending Bellégo, Benatia and Dortet-Bernadet 2025 and Callaway and Sant'Anna 2021), `pt_assumption` selects the parallel trends restrictions imposed in generalized method of moments (GMM) estimation ([Parallel trends](../methods/parallel-trends.md)). The default `pt_assumption = "post"` (`PT-post`) assumes conditional parallel trends from each cohort's reference period onward. Setting `pt_assumption = "all"` (`PT-all`) with an explicit common calendar start date `pt_start` keeps every `PT-post` condition and each cohort's original reference period, and adds zero differential restrictions on represented untreated preperiod steps between `pt_start` and each cohort's reference period.

This worked example reuses the exact unit table (`units`), incomplete panel (`panel`), complete latent outcome vector (`y_complete`), and observation mask (`observed`) created in [Choose IPW, OR or DR in an incomplete panel](treatment-choice.md). Start a fresh R session, run that script, then run this example immediately: other pages of this guide use the names `units`, `panel` and `observed` for different objects.

## Question

Holding `panel`, `units`, `xformla = ~ x`, `control_group = "nevertreated"`, `est_method = "dr"`, `selection = "none"`, `base_period = "varying"`, and `bstrap = FALSE, cband = FALSE` fixed, what changes when we switch from `pt_assumption = "post"` to `pt_assumption = "all"` with `pt_start = 1`? First, how do the two fits compare under the baseline law where conditional parallel trends hold across all dates $1$ to $6$? Second, what happens when we add a cohort specific pre baseline untreated pulse at period $1$ that violates `PT-all` from period $1$ while leaving `PT-post` and all posttreatment causal effects unchanged?

## Data

### Prerequisite objects and pre baseline pulse modification

Run [Choose IPW, OR or DR in an incomplete panel](treatment-choice.md#code) first to create `units` ($N = 2000$ units with `id`, `first_treat` in $\{0, 3, 5\}$, and `x`), `panel` ($N \times T = 12000$ rows across stored calendar dates $t \in \{1, \dots, 6\}$ with `NA` in `y` where unobserved), `y_complete` (unmasked latent outcomes), and `observed` (logical observation mask). We compare `PT-post` and `PT-all` (`pt_start = 1`) on two outcome vectors sharing `units`, latent shocks, and `observed`:

1. **Baseline panel (`panel`):** Untreated outcomes follow $Y_{i,t}(0) = \alpha_i + 0.1 t + 0.3 X_i t + \varepsilon_{i,t}$ at all dates $t \in \{1, \dots, 6\}$. With `anticipation = 0`, cohort $3$ has reference period $2$ and cohort $5$ has reference period $4$. Conditional parallel trends hold across every date pair, so both `PT-post` and `PT-all` (`pt_start = 1`) hold in the population.
2. **Pre baseline pulse panel (`pulse_panel`):** Before applying the original mask `observed`, we add $+0.5 \times \mathbf{1}\{G_i > 0, t = 1\}$ to `y_complete`, so $Y_{i,t}^{\text{pulse}}(0) = Y_{i,t}(0) + 0.5 \times \mathbf{1}\{G_i > 0, t = 1\}$ and $Y_{i,t}^{\text{pulse}} = Y_{i,t} + 0.5 \times \mathbf{1}\{G_i > 0, t = 1\}$. Applying `ifelse(observed, y_pulse_complete, NA_real_)` preserves every `NA` entry and the exact observation pattern of `panel`.

This modification is a cohort specific pre baseline untreated pulse, **not** anticipation of treatment:

- At every date $t \in \{1, \dots, 6\}$, the potential outcome difference $Y_{i,t}^{\text{pulse}}(g) - Y_{i,t}^{\text{pulse}}(0)$ remains $0$ for $t < g$ and $1 + 0.25(t - g)$ for $t \ge g$, so `anticipation = 0` holds and the posttreatment target remains $\text{ATT}(g, t) = 1 + 0.25(t - g)$.
- Because the pulse acts only at $t = 1$, untreated changes $Y_{i,t}^{\text{pulse}}(0) - Y_{i,s}^{\text{pulse}}(0)$ for all pairs $2 \le s < t$ equal those of the baseline law. Since cohort reference periods are $2$ and $4$, `PT-post` still holds completely.
- From $t = 1$ to $t = 2$, however, the conditional mean untreated change for cohorts $3$ and $5$ differs by $-0.5$ from that of never treated controls ($G_i = 0$). Thus `PT-all` with `pt_start = 1` fails.

## Assumptions

The four fits vary only `data`, `pt_assumption`, and `pt_start` ([Options](../reference/options.md)).

| Argument | Setting | Purpose |
| --- | --- | --- |
| `data`, `frame` | `panel` or `pulse_panel`, `units` | $N \times T$ panel with `NA` outcomes and fixed $N$ row unit table ([Options](../reference/options.md#chained-arguments)). |
| `xformla`, `control_group`, `est_method`, `selection` | `~ x`, `"nevertreated"`, `"dr"`, `"none"` | Doubly robust treatment adjustment and representative observed pair changes ([Treatment adjustment](../methods/treatment.md)). |
| `pt_assumption`, `pt_start` | `"post"` (`NULL`) or `"all"` (`1`) | Switches between `PT-post` and `PT-all` from common calendar start date $1$ ([Parallel trends](../methods/parallel-trends.md)). |
| `base_period`, `bstrap`, `cband` | `"varying"`, `FALSE`, `FALSE` | Consecutive unrestricted preperiod reports and analytical standard errors ([Options](../reference/options.md#inference-and-computation)). |

### How `PT-post` and `PT-all` use the comparison catalog

1. **Earlier bridges under `PT-post` have no artificial endpoint cutoff:** Under `PT-post`, two date comparisons touching dates before a cohort's reference period remain in the catalog whenever their unrestricted pre reference differentials cancel across a bridge ([Parallel trends](../methods/parallel-trends.md#comparisons-through-earlier-dates)). For cohort $3$ (reference period $2$), combining $(1, 3)$ and $(1, 2)$ cancels the period $1$ to $2$ differential because $E[Y_3 - Y_1] - E[Y_2 - Y_1] = E[Y_3 - Y_2]$ in the fixed cohort and control populations. This is a population identity: each observed pair must represent the same populations, with valid control/model support, and its fitted estimates need not satisfy an exact realized equality.
2. **Added restrictions under `PT-all` (`pt_start = 1`):** `PT-all` retains every `PT-post` condition and each cohort's original reference period ($2$ for cohort $3$, $4$ for cohort $5$), and constrains represented operator specific untreated pre reference steps from `pt_start = 1` up to each cohort's reference period to zero ($1$ to $2$ for cohort $3$, and steps from $1$ up to $4$ for cohort $5$). `pt_start` must be one finite numeric date in the stored calendar common to all cohorts (omitting it or passing `NULL` under `"all"`, or passing a non `NULL` value under `"post"`, raises an error).
3. **GMM moment restrictions and separate `base_period`:** `PT-all` adds restrictions to the existing chained GMM system over the same primitive comparison catalog and control pools. It retains their complete unit influence-function construction, but changes the combination: the final target influence functions, covariance and standard errors can change. It does not switch to a cross sectional estimator or introduce a new weighting rule. Meanwhile, `base_period` separately controls how unrestricted preperiod comparisons (`role == "pre"`) are reported.

## Code



```r
# Objects from the treatment-choice example. exists() searches only the current
# environment, so base::units() cannot satisfy the check.
needed <- c("units", "panel", "y_complete", "observed")
found <- vapply(needed, exists, logical(1), envir = environment(), inherits = FALSE)
if (!all(found)) {
  stop("Run 'Choose IPW, OR or DR in an incomplete panel' first, in a fresh R session. ",
       "Not found: ", paste(needed[!found], collapse = ", "))
}
if (!is.data.frame(units) || !all(c("id", "first_treat", "x") %in% names(units)) ||
    !is.data.frame(panel) ||
    !all(c("id", "period", "first_treat", "x", "y") %in% names(panel)) ||
    !is.numeric(panel$y) || !is.numeric(y_complete) || !is.logical(observed) ||
    length(y_complete) != nrow(panel) || !identical(!is.na(panel$y), observed)) {
  stop("units, panel, y_complete and observed do not match 'Choose IPW, OR or DR ",
       "in an incomplete panel'. Start a fresh R session, run that example, ",
       "then run this one.")
}

# 1. Baseline panel fits: PT-post vs PT-all with explicit pt_start = 1.
fit_pt_post <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = panel, frame = units,
  control_group = "nevertreated", est_method = "dr", selection = "none",
  pt_assumption = "post", base_period = "varying", bstrap = FALSE, cband = FALSE
)
fit_pt_all <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = panel, frame = units,
  control_group = "nevertreated", est_method = "dr", selection = "none",
  pt_assumption = "all", pt_start = 1, base_period = "varying", bstrap = FALSE, cband = FALSE
)

utils::str(fit_pt_post$chained$pt)
utils::str(fit_pt_all$chained$pt)
fit_pt_all$chained$pt[c(
  "pt_assumption", "pt_start", "n_extra_calendar_dates",
  "n_restricted_pre_parameters", "no_extra_restrictions"
)]
fit_pt_all$chained$pt$by_cohort

cols <- c("target.id", "group", "time", "role", "available", "estimate", "std.error")
baseline_table <- merge(
  cdid::tidy(fit_pt_post)[, cols],
  cdid::tidy(fit_pt_all)[, c("target.id", "estimate", "std.error")],
  by = "target.id", suffixes = c("_pt_post", "_pt_all")
)
baseline_table$population_truth <- ifelse(
  baseline_table$role == "post", 1 + 0.25 * (baseline_table$time - baseline_table$group), 0
)
baseline_table[baseline_table$role == "post", ]
baseline_table[baseline_table$role == "pre", ]

list(
  pt_post = fit_pt_post$chained$pretest[c("status", "reason", "df", "W", "Wpval")],
  pt_all = fit_pt_all$chained$pretest[c("status", "reason", "df", "W", "Wpval")]
)

# 2. Pre baseline untreated pulse: add +0.5 * 1{G > 0, period == 1} to y_complete
# BEFORE applying the original `observed` mask.
y_pulse_complete <- y_complete + ifelse(panel$first_treat > 0L & panel$period == 1L, 0.5, 0)
pulse_panel <- panel
pulse_panel$y <- ifelse(observed, y_pulse_complete, NA_real_)

fit_pulse_post <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = pulse_panel, frame = units,
  control_group = "nevertreated", est_method = "dr", selection = "none",
  pt_assumption = "post", base_period = "varying", bstrap = FALSE, cband = FALSE
)
fit_pulse_all <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = pulse_panel, frame = units,
  control_group = "nevertreated", est_method = "dr", selection = "none",
  pt_assumption = "all", pt_start = 1, base_period = "varying", bstrap = FALSE, cband = FALSE
)

pulse_table <- merge(
  cdid::tidy(fit_pulse_post)[, cols],
  cdid::tidy(fit_pulse_all)[, c("target.id", "estimate", "std.error")],
  by = "target.id", suffixes = c("_pt_post", "_pt_all")
)
pulse_table[pulse_table$role == "post", ]
pulse_table[pulse_table$role == "pre", ]

list(
  pulse_pt_post = fit_pulse_post$chained$pretest[c("status", "reason", "df", "W", "Wpval")],
  pulse_pt_all = fit_pulse_all$chained$pretest[c("status", "reason", "df", "W", "Wpval")]
)
```

<details>
<summary>Console output</summary>

```text
List of 2
 $ pt_assumption: chr "post"
 $ pt_start     : NULL
List of 14
 $ pt_assumption                : chr "all"
 $ pt_start                     : num 1
 $ start_selection              : chr "fixed exogenously"
 $ by_cohort                    :'data.frame':	2 obs. of  10 variables:
  ..$ g                          : num [1:2] 3 5
  ..$ ref                        : num [1:2] 2 4
  ..$ pt_start                   : num [1:2] 1 1
  ..$ baseline_status            : chr [1:2] "stored_reference" "stored_reference"
  ..$ n_extra_calendar_dates     : int [1:2] 1 3
  ..$ n_extra_represented_dates  : int [1:2] 1 3
  ..$ n_restricted_pre_parameters: int [1:2] 1 3
  ..$ n_zero_loading_links       : int [1:2] 1 6
  ..$ no_extra_calendar_window   : logi [1:2] FALSE FALSE
  ..$ no_extra_restrictions      : logi [1:2] FALSE FALSE
 $ n_extra_calendar_dates       : int 4
 $ n_restricted_pre_parameters  : int 4
 $ n_zero_loading_links         : int 7
 $ no_extra_calendar_window     : logi FALSE
 $ no_extra_restrictions        : logi FALSE
 $ restriction_is_identity      : logi FALSE
 $ free_post_param_indices      : int [1:6] 2 3 4 5 9 10
 $ restricted_post_param_indices: int [1:4] 1 6 7 8
 $ assumptions                  : chr "Common start fixed exogenously. Baseline PT-post assumptions are retained at every qualified causal date, inclu"| __truncated__
 $ pretest_interpretation       : chr "Unrestricted direct pre-comparisons retain the full requested nonreference pretrend null. Pre-start contrasts c"| __truncated__
$pt_assumption
[1] "all"

$pt_start
[1] 1

$n_extra_calendar_dates
[1] 4

$n_restricted_pre_parameters
[1] 4

$no_extra_restrictions
[1] FALSE

  g ref pt_start  baseline_status n_extra_calendar_dates n_extra_represented_dates
1 3   2        1 stored_reference                      1                         1
2 5   4        1 stored_reference                      3                         3
  n_restricted_pre_parameters n_zero_loading_links no_extra_calendar_window no_extra_restrictions
1                           1                    1                    FALSE                 FALSE
2                           3                    6                    FALSE                 FALSE
     target.id group time role available estimate_pt_post std.error_pt_post estimate_pt_all
1 ATT(g=3,t=3)     3    3 post      TRUE         1.007945        0.05946479       0.9608242
2 ATT(g=3,t=4)     3    4 post      TRUE         1.298627        0.06084339       1.2861154
3 ATT(g=3,t=5)     3    5 post      TRUE         1.528495        0.05908227       1.5149354
4 ATT(g=3,t=6)     3    6 post      TRUE         1.728774        0.06184471       1.7112627
5 ATT(g=5,t=5)     5    5 post      TRUE         1.053470        0.06544368       1.0475218
6 ATT(g=5,t=6)     5    6 post      TRUE         1.143902        0.06712190       1.1287508
  std.error_pt_all population_truth
1       0.04654280             1.00
2       0.04694743             1.25
3       0.04989078             1.50
4       0.05269931             1.75
5       0.05129701             1.00
6       0.05333004             1.25
      target.id group time role available estimate_pt_post std.error_pt_post estimate_pt_all
7  pre(g=3,t=2)     3    2  pre      TRUE     0.0733146863        0.07110089    0.0733146863
8  pre(g=5,t=2)     5    2  pre      TRUE    -0.0003672074        0.07605146   -0.0003672074
9  pre(g=5,t=3)     5    3  pre      TRUE     0.1281262590        0.07029130    0.1281262590
10 pre(g=5,t=4)     5    4  pre      TRUE    -0.0634728529        0.08537601   -0.0634728529
   std.error_pt_all population_truth
7        0.07110089                0
8        0.07605146                0
9        0.07029130                0
10       0.08537601                0
$pt_post
$pt_post$status
[1] "computed"

$pt_post$<NA>
NULL

$pt_post$df
[1] 4

$pt_post$W
[1] 4.768994

$pt_post$Wpval
[1] 0.3118316


$pt_all
$pt_all$status
[1] "computed"

$pt_all$<NA>
NULL

$pt_all$df
[1] 4

$pt_all$W
[1] 4.768994

$pt_all$Wpval
[1] 0.3118316


     target.id group time role available estimate_pt_post std.error_pt_post estimate_pt_all
1 ATT(g=3,t=3)     3    3 post      TRUE         1.007945        0.05946479       0.8167734
2 ATT(g=3,t=4)     3    4 post      TRUE         1.298627        0.06084339       1.1234910
3 ATT(g=3,t=5)     3    5 post      TRUE         1.528495        0.05908227       1.3036359
4 ATT(g=3,t=6)     3    6 post      TRUE         1.728774        0.06184471       1.5002007
5 ATT(g=5,t=5)     5    5 post      TRUE         1.053470        0.06544368       0.9324565
6 ATT(g=5,t=6)     5    6 post      TRUE         1.143902        0.06712190       1.0072939
  std.error_pt_all
1       0.04654280
2       0.04694743
3       0.04989078
4       0.05269931
5       0.05129701
6       0.05333004
      target.id group time role available estimate_pt_post std.error_pt_post estimate_pt_all
7  pre(g=3,t=2)     3    2  pre      TRUE      -0.42668531        0.07110089     -0.42668531
8  pre(g=5,t=2)     5    2  pre      TRUE      -0.50036721        0.07605146     -0.50036721
9  pre(g=5,t=3)     5    3  pre      TRUE       0.12812626        0.07029130      0.12812626
10 pre(g=5,t=4)     5    4  pre      TRUE      -0.06347285        0.08537601     -0.06347285
   std.error_pt_all
7        0.07110089
8        0.07605146
9        0.07029130
10       0.08537601
$pulse_pt_post
$pulse_pt_post$status
[1] "computed"

$pulse_pt_post$<NA>
NULL

$pulse_pt_post$df
[1] 4

$pulse_pt_post$W
[1] 57.4808

$pulse_pt_post$Wpval
[1] 9.807315e-12


$pulse_pt_all
$pulse_pt_all$status
[1] "computed"

$pulse_pt_all$<NA>
NULL

$pulse_pt_all$df
[1] 4

$pulse_pt_all$W
[1] 57.4808

$pulse_pt_all$Wpval
[1] 9.807315e-12
```

</details>

## Result

As documented in [Parallel trends](../methods/parallel-trends.md#where-the-result-records-the-choice), `fit$chained$pt` records the parallel trends configuration:

- `fit_pt_post$chained$pt` records `pt_assumption = "post"` and `pt_start = NULL`.
- `fit_pt_all$chained$pt` records `pt_assumption = "all"`, `pt_start = 1`, `n_extra_calendar_dates`, `n_restricted_pre_parameters`, `no_extra_restrictions`, the assumption text, and `by_cohort` (with columns `g`, `ref`, `pt_start`, `n_extra_calendar_dates`, `n_extra_represented_dates`, `n_restricted_pre_parameters`, `no_extra_calendar_window`, and `no_extra_restrictions` for `g = 3, ref = 2` and `g = 5, ref = 4`).
- In `baseline_table` and `pulse_table`, rows with `role == "post"` report the six posttreatment effects $(g, t) \in \{(3, 3), (3, 4), (3, 5), (3, 6), (5, 5), (5, 6)\}$ with population truth $1 + 0.25(t - g)$, while rows with `role == "pre"` report the consecutive unrestricted preperiod comparisons $(3, 2), (5, 2), (5, 3), (5, 4)$.

## Interpretation

### Baseline law versus pre baseline pulse law

On `panel`, both `fit_pt_post` and `fit_pt_all` impose true population restrictions and, under the usual identification and large sample regularity conditions, are consistent for $\text{ATT}(g, t) = 1 + 0.25(t - g)$. Their finite sample point estimates and standard errors can differ because `fit_pt_all` imposes pre reference zero differential restrictions from period $1$ up to each cohort's reference period in the two step GMM solve. A difference can reflect those restrictions and is not by itself evidence of a coding error; an unexplained difference still merits investigation. A smaller standard error is not guaranteed in every cell or finite draw ([How to compare precision](precision.md#restrictions-can-change-information-and-robustness)).

On `pulse_panel`, the $+0.5$ untreated pulse at period $1$ leaves all untreated changes from period $2$ onward unchanged. Because cohort reference periods are $2$ and $4$ and any bridge through period $1$ under `PT-post` cancels the unrestricted period $1$ to $2$ differential, `fit_pulse_post` remains consistent for $\text{ATT}(g, t) = 1 + 0.25(t - g)$. By contrast, `fit_pulse_all` (`pt_start = 1`) constrains the period $1$ to $2$ untreated differential to zero when its population value is $-0.5$. Those false restrictions can bias the posttreatment estimates through the joint GMM combination, although neither bias in every cell nor a pretrend test rejection is guaranteed in a particular finite draw.

### Unrestricted preperiod reports, pretrend availability, and choosing `pt_start`

1. **Preperiod cells (`role == "pre"`) are never zeroed out by `PT-all`:** Under both `PT-post` and `PT-all`, reported preperiod cells are direct, unrestricted comparisons, and restricted zero pre parameters never supply the unconditional Wald test (`W`, `Wpval`).
2. **Full pretrend test availability can be `"unavailable"`:** `att_gt()` computes `W` and `Wpval` only when every declared nonreference `"pre"` cell is available and the joint covariance passes its Cholesky and conditioning check; with failed support or covariance checks, `W` and `Wpval` are `NULL` and `fit$chained$pretest$status` is `"unavailable"` with its reason. In other settings, no declared preperiod comparison gives `"not_applicable"`, while `compute_inffunc = FALSE` gives `"not_computed"` ([Results](../reference/results.md#pretrend-tests)).
3. **Never choose `pt_start` to fit a pretest:** Select `pt_start` substantively from domain knowledge before viewing estimates, standard errors, or pretrend tests ([Parallel trends](../methods/parallel-trends.md#choosing-the-start-date)). A single script or single draw does not establish coverage or efficiency; repeated sample evidence remains pending.

## Next step

- Read [Parallel trends](../methods/parallel-trends.md), [Results](../reference/results.md), [Inference](inference.md), and [How to compare precision](precision.md), or return to [Choose IPW, OR or DR in an incomplete panel](treatment-choice.md).

## References

- Bellégo, C., Benatia, D. and Dortet-Bernadet, V. (2025). The Chained Difference-in-Differences. *Journal of Econometrics*. https://doi.org/10.1016/j.jeconom.2024.105783
- Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with Multiple Time Periods. *Journal of Econometrics* 225(2). https://doi.org/10.1016/j.jeconom.2020.12.001

---

**Reading path, step 14 of 20** · Previous: [Parallel trends](../methods/parallel-trends.md) · Next: [GMM weighting and covariance regularization](../methods/gmm.md)

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
| 12 | [Impute actual changes](change-imputation.md) |
| 13 | [Parallel trends](../methods/parallel-trends.md) |
| 14 | **Add PT-all restrictions** (this page) |
| 15 | [GMM weighting and covariance regularization](../methods/gmm.md) |
| 16 | [Robustness](robustness.md) |
| 17 | [Inference](inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
