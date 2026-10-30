# Moving from did

If you already use `did`, most of `cdid` will look familiar. This page explains what stays the same, how the two estimation modes differ, how each package treats an incomplete panel and when chained estimates should agree with those of `did`.

## What stays the same

- **Functions.** `att_gt()`, `aggte()`, `ggdid()`, `tidy()`, `glance()` and `conditional_did_pretest()`. The generics `ggdid()`, `tidy()` and `glance()` are exported unchanged from `did`; `cdid` registers methods only for its own result classes, so `did` results keep `did`'s methods.
- **Arguments.** `att_gt()` retains the upstream argument names and order, with `allow_unbalanced_panel = NULL` as a mode-dependent sentinel in place of `did`'s `FALSE` default. `aggte()` keeps the signature of `did` 2.5.1. `conditional_did_pretest()` retains the 19 upstream arguments in order, with the same unbalanced-panel sentinel. The arguments that only chained mode uses come after `...` in `att_gt()` and after `did`'s arguments in `conditional_did_pretest()`.
- **Defaults.** `att_gt()` uses `est_method = "dr"`, and `conditional_did_pretest()` uses `est_method = "ipw"`, as in `did`. Units never treated have `gname` equal to 0.
- **Results.** Chained results inherit `did`'s classes (`MP` for `att_gt()`, `AGGTEobj` for `aggte()`), and `did`'s own methods tabulate and plot the estimates. `cdid` adds a short chained record to the display and a `chained` element to the object.

Qualify calls as `cdid::att_gt()` and `cdid::aggte()`. A later `library(did)` then cannot change which function runs.

## Two estimation modes

| `mode` | When it applies | What runs | Returned object |
| --- | --- | --- | --- |
| `NULL` (default) | resolved from `panel` | chained mode when `panel = TRUE`, ordinary mode when `panel = FALSE` | as below |
| `"chained"` | panel data | chained DiD | `cdid_MP`, which inherits `MP` |
| `"ordinary"` | any data `did` accepts | `did::att_gt()`, unchanged | `did`'s `MP` |

In ordinary mode, `cdid` forwards only the arguments you supplied, so `did` applies its own defaults. Each argument is evaluated once, in your environment, when `did` needs it. The one exception is `panel`: with `mode = NULL`, `cdid` evaluates it first to choose the mode. The chained arguments are neither evaluated nor forwarded. `allow_unbalanced_panel = NULL` becomes `did`'s default `FALSE`. The delegated call is recorded in `DIDparams$call`, so printed summaries show it. `aggte()` passes ordinary results to `did::aggte()` unchanged, and `conditional_did_pretest(mode = "ordinary")` runs `did::conditional_did_pretest()`.

Aggregate chained results with `cdid::aggte()`. A chained result stores the state that aggregation needs but not the panel data, so `did::aggte()` cannot aggregate it.

## Repeated cross sections

With `panel = FALSE`, the default mode is ordinary, and `did`'s repeated cross section (RCS) estimators run unchanged. Chained mode needs panel data and refuses `panel = FALSE`. When different units are observed at each date there are no changes within units to chain, so the RCS approach is the available one, under its own sampling assumptions.

## Incomplete panels in did and in chained mode

`did` already accepts unbalanced panels and offers two treatments of them through `allow_unbalanced_panel`. Chained mode uses a third.

| Approach | Observations used | What each comparison uses |
| --- | --- | --- |
| `did` with `allow_unbalanced_panel = FALSE`, its default (also ordinary mode with `NULL`) | rows with missing values in the variables used are dropped, then only units observed at every date are kept | changes within units of this balanced subset |
| `did` with `allow_unbalanced_panel = TRUE` | rows with complete, finite values of the required variables | `did`'s repeated cross section estimators, which compare outcome levels among the units observed at each date |
| `cdid` chained mode | every unit of the unit table; units with at least two observed outcomes contribute changes | changes between two observed dates of the same unit, combined across dates |

These descriptions refer to `did` 2.5.1, the version `cdid` is tested with. In chained mode, `allow_unbalanced_panel = NULL` (the default) keeps every eligible unit, including in balanced panels, and an explicit `FALSE` is refused with an error that points to `mode = "ordinary"`.

The choice matters for three reasons.

1. **The balanced subset can describe a different population.** If observation depends on characteristics related to how outcomes evolve, the units observed at every date differ from the study population. When effects vary across units, an analysis of that subset estimates effects for the subset.
2. **Comparisons of levels need conditions on levels.** Callaway and Sant'Anna's repeated cross section setting assumes that the sample at each date represents the population at that date and that the distribution of cohorts and covariates is stable across dates. When observation depends on a fixed unit component, for example when units with persistently low outcomes leave the panel, the observed levels at later dates stop representing the population. Changes within units remove that component, so chained DiD can still rely on the observed changes.
3. **Units with one observed outcome form no change.** In chained mode they remain part of the study population through the unit table and still contribute their cohort, covariates and observation information. These units do not enter the regression fits for outcome changes, because computing a difference requires at least two observations per unit.

Each approach needs its own assumptions. [Handling missing outcomes](methods/missing-outcomes.md) states the chained assumptions about observed changes.

## When chained and did estimates should agree

| Setting | Point estimates of chained and ordinary mode |
| --- | --- |
| Baseline `est_method = "ipw"`, `selection = "none"`, complete panel, units never treated as controls, no covariates, each unit's weight fixed across comparisons, PT-post, the same `anticipation` and `base_period` | Expected to agree up to floating point rounding, because consecutive changes add up to the long comparisons `did` uses. [Get started](get-started.md#optional-compare-with-did) gives the code. |
| OR or DR with different fitted specifications or procedures | Can differ. Chained OR imposes one slope vector shared across comparison blocks, while `did` ordinarily fits separate comparison regressions. The method label alone is not a reason for a difference; identical formula notation need not impose identical coefficient restrictions. |
| Covariates with `est_method = "ipw"` | This guide makes no general claim. Compare the two fits directly. |
| Units not yet treated as controls | The chained catalog can retain additional fixed-horizon control-pool moments for the same date pair and combine them jointly. This differs from using only the cell's ordinary control pool; see the concrete example below. |
| `pt_assumption = "all"` | Need not agree. The added restrictions can change estimates with complete outcomes, and `did` has no such option. |
| Incomplete panel | The approaches use different observations and assumptions, as in the table above. |
| Sampling weights that change within units | Chained mode requires a fixed reference weight through `fix_weights`, described in [Changing covariates and weights](reference/options.md#changing-covariates-and-weights). |


## Why can balanced results differ, and which should I use?

With the same complete sample, balance removes missing-outcome differences. A fair comparison must also match the fitted specification and procedure: coefficient restrictions, nuisance fitting, moments and weighting. If that full procedure is identical, the point estimates should agree up to numerical rounding. The labels OR or DR, or the same `xformla` notation, do not themselves establish a difference or an equality.

| Matched balanced panel setting | Why point estimates can agree or differ | Implication for the choice |
| --- | --- | --- |
| Baseline built-in IPW, never treated controls, no covariates, PT-post and each unit's weight fixed across comparisons | The common cohort and control means telescope, so the chained combination reduces to the long comparisons used by `did`, up to floating point rounding. | This is a useful common-case check, not a claim of equality for every method. Match the input, weights, anticipation and other settings before comparing. |
| Treatment OR or DR when slope restrictions or fitting differ | The shared slope and local intercept fit differs from separate comparison regressions unless their specifications and fitted procedure have been constrained to coincide. Even different nuisance fits need not change the final estimate if their correction cancels. | If covariate effects on untreated changes vary across dates and the shared slope restriction lacks justification, separate comparison regressions are more flexible. If pooling is credible, sharing information may improve precision, but no uniform efficiency advantage is established. |
| Not-yet-treated controls | The chained catalog retains original comparisons and distinct later fixed-horizon reuses. Multiple valid control-pool moments can be combined for one effect even with complete outcomes and no covariates. | Under the full NYT parallel-trends and overlap conditions, these are valid additional moments. They may add information, but no universal finite-sample or efficiency advantage is claimed. |
| PT-all | Chained estimation adds represented untreated-trend restrictions from the chosen common start; `did` has no corresponding PT-all option. | Use these stronger restrictions only when credible. They may add identification or precision; false restrictions can bias the estimates. |

### A balanced NYT example

Consider a regular calendar, zero anticipation, cohorts first treated in periods 3 and 5, and units never treated (`G = 0`). Every outcome is observed, with no covariates and equal weights. For ATT(3,3), the reference date is 2.

| Comparison of the change from 2 to 3 | Ordinary cell and chained catalog |
| --- | --- |
| Cohort 3 against cohorts 5 and 0 together | The ordinary `did` cell uses this pool; the chained catalog retains it too. |
| Cohort 3 against cohort 0 alone | The chained catalog also retains this same date pair with the control pool fixed at a later horizon, such as 5, by which time cohort 5 is treated. |

Under the full not-yet-treated parallel-trends assumption, the target cohort's untreated changes agree with each eligible later control tail, not only the largest tail at the cell's own date. Both comparisons therefore have population expectation ATT(3,3). This implication uses the same full NYT assumption at every later cutoff; it adds no parallel-trends restriction before reference date 2 and is not PT-all. A weaker assumption asserted only for the combined `{5, 0}` pool at period 3 would not justify the never-treated-only reuse.

The two sample control averages need not coincide. `cdid` combines these and the other valid moments jointly by GMM, whereas the ordinary cell uses its own comparison. The additional moment restrictions can therefore change finite sample point estimates even in this balanced example. Joint weighting can also affect other targets through covariance with these moments. This identifies the source of a possible difference; it does not rank efficiency or establish the assumptions from the data.

For treatment DR, the correctness alternative remains: under valid observation reconstruction, parallel trends, overlap and the other regularity conditions, a correct treatment propensity model can protect a misspecified pooled counterfactual regression, or a correct pooled regression can protect a misspecified propensity model. Pooling must be correct across the blocks that share its slope for the regression branch to apply. DR is not protection against invalid parallel trends or invalid handling of missing outcomes.

Point estimates, standard errors and multiplier bands are separate comparisons. Agreement of point estimates does not establish agreement of influence functions or analytical standard errors. Different multiplier draws can also change reported bootstrap standard errors and critical values without changing the fitted effects; a reproducible inference comparison must align its settings and random number procedure.

## Differences in the chained interface

| Topic | Chained mode |
| --- | --- |
| `idname` | Required, because units are aligned by identity. |
| `clustervars` | At most one cluster variable in addition to `idname`. |
| `allow_unbalanced_panel` | `NULL` keeps every eligible unit; `FALSE` is refused. |
| `est_method` | `"dr"` (default), `"ipw"`, `"reg"` or a custom IPW function with `did`'s panel interface. Any other string is refused before the data are read, and nothing is substituted. |
| `...` | Passed to a custom `est_method` function; refused when `est_method` is a string. |
| `fix_weights` | `NULL`, `"first_period"` or `"base_period"`; `"varying"` is refused. |
| `compute_inffunc = FALSE` | Point estimates only, as in `did`: no influence functions, standard errors, pretrend test or aggregation. |
| `pl = TRUE`, `cores = k` | Multiplier draws on `k` worker processes, which you must declare with `control = list(workers = k)`. Results equal those of `cores = 1`. |
| `faster_mode` | Accepted with its `did` meaning; chained mode has one computation engine. |
| `print_details = TRUE` | Progress messages, which `suppressMessages()` silences. |

Chained mode adds arguments for the unit table, missing outcomes, parallel trends and GMM weighting. The [options reference](reference/options.md#chained-arguments) lists them.

---

**Reading path, step 2 of 20** · Previous: [Get started](get-started.md) · Next: [Incomplete panels](examples/incomplete-panels.md)

<details>
<summary>All steps of the reading path</summary>

| Step | Page |
| --- | --- |
| 1 | [Get started](get-started.md) |
| 2 | **Moving from did** (this page) |
| 3 | [Incomplete panels](examples/incomplete-panels.md) |
| 4 | [Identify effects when the observed units change composition](examples/identification.md) |
| 5 | [How to compare precision](examples/precision.md) |
| 6 | [Prepare your panel and unit table](help/prepare-data.md) |
| 7 | [Treatment adjustment](methods/treatment.md) |
| 8 | [Choose IPW, OR or DR in an incomplete panel](examples/treatment-choice.md) |
| 9 | [Handling missing outcomes](methods/missing-outcomes.md) |
| 10 | [Sampling information for the observation model](methods/sampling-information.md) |
| 11 | [Choose an observation model](examples/observation-models.md) |
| 12 | [Impute actual changes](examples/change-imputation.md) |
| 13 | [Parallel trends](methods/parallel-trends.md) |
| 14 | [Add PT-all restrictions](examples/pt-all-comparison.md) |
| 15 | [GMM weighting and covariance regularization](methods/gmm.md) |
| 16 | [Robustness](examples/robustness.md) |
| 17 | [Inference](examples/inference.md) |
| 18 | [Options](reference/options.md) |
| 19 | [Results](reference/results.md) |
| 20 | [Diagnose errors and unavailable results](help/troubleshooting.md) |

[Guide home](README.md) · [Package overview](../README.md)

</details>
