# Parallel trends: PT-post and PT-all

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

`pt_assumption` chooses which parallel trends restrictions the estimator imposes. With `"post"`, the default, untreated outcomes of each cohort and its controls are assumed to move in parallel from the cohort's reference period onward. With `"all"`, every one of those restrictions is kept and further restrictions are added between a common start date `pt_start` and each cohort's reference period. How effects before treatment are reported is a separate choice, made with `base_period`.

## Reference period and anticipation

The reference period of cohort g is the latest date d in the stored calendar (the sorted dates present in the data) such that `d + anticipation < g`. This cutoff uses numeric date values; it is distinct from the stored-calendar index gap used for OR/imputation durations. Every effect of the cohort is measured relative to it. Dates after the reference period and before g are anticipation periods, in which effects may already be nonzero.

| Cohort g | `anticipation` | Reference period with dates 1, 2, ..., 6 | Anticipation periods |
| --- | --- | --- | --- |
| 3 | 0 | 2 | none |
| 3 | 1 | 1 | 2 |
| 5 | 0 | 4 | none |

## PT-post

Under PT-post, for each cohort g and each qualified causal date used by its catalog, the average untreated outcome of the cohort would have changed from the reference period to t by the same amount as that of its comparison group, conditional on the covariates of `xformla` when it has any. The comparison group consists of units never treated or units not yet treated, according to `control_group`. Nothing is assumed about trends before the reference period.

There is no endpoint cutoff that discards comparisons touching earlier dates. Their unrestricted pre-reference contributions must cancel in an identified target, as in the next section. An effect that could be reached only by assuming parallel trends before the reference period is reported as unavailable under PT-post.

## How the estimator uses the comparisons

The reference period anchors the target. It need not be observed for every unit, and a useful comparison need not start there. The package works with a joint system rather than selecting one path through the dates:

| Step | What enters the calculation | What it requires |
| --- | --- | --- |
| 1. Form comparisons | All valid ordered-pair cohort-minus-eligible-control moments in the retained catalog, including longer spans and valid fixed-horizon control-pool reuses | Observed pairs in the cohort and eligible controls, overlap and the support needed by active models |
| 2. Reconstruct population changes | The treatment adjustment and, when requested, the observation weighting, change imputation or augmented correction | Each comparison must reconstruct its stated cohort and control populations under the same fixed covariate and weight reference law whenever cancellation is used |
| 3. Build the PT model | Population comparison means $d = A\theta$, with treatment effects and unrestricted earlier trend differentials in $\theta$ | PT-post leaves the pre-reference differentials free. PT-all sets only its represented in-window differentials to zero. |
| 4. Recover supported targets | GMM combines the estimated comparisons, accounting for their joint covariance and fitted-model uncertainty | A target contrast must be recoverable from the rows of the maintained model. Covariance regularization is a separate step. |

In $d = A\theta$, a row of $A$ records which effects or free trend differentials enter one comparison, with their signs. A target is estimable when some combination of these rows isolates its contrast. This is the row-space condition. Under PT-post, every unrestricted term in that combination must cancel. A visible path in an observation graph alone does not establish this, especially when different comparisons use different not-yet-treated control pools.

The sample vector $\hat d$ estimates $d$. It need not satisfy the population equations exactly, and GMM uses the whole retained system to combine the sample information. Exact structural reductions preserve all original comparisons, their full metric and their complete influence functions. They do not replace the system with a smaller identification forest. [GMM weighting](gmm.md#exact-structural-reductions-and-the-full-comparison-metric) explains this distinction.

## Comparisons through earlier dates

Consider stored dates 1, 2, 3, 4 and 5, cohort $g = 3$, zero anticipation and reference period $b_3 = 2$. Use the same never-treated control pool and the same reference population law in both comparisons. Suppose units are observed in windows $\{1,2\}$ or $\{1,4\}$, with no unit observed at both 2 and 4. Dates 3 and 5 can remain in the stored calendar through rows whose outcomes are `NA`.

| Observation group (cohort 3 and controls alike) | Period 1 | Period 2 | Period 3 | Period 4 | Period 5 | Change it supplies |
| --- | --- | --- | --- | --- | --- | --- |
| First group | observed | observed | missing | missing | missing | 1 to 2 |
| Second group | observed | missing | missing | observed | missing | 1 to 4 |

Let $a_2-a_1$ denote the unrestricted untreated cohort-minus-control trend differential from 1 to 2. Write $d_{st}$ for the population comparison from $s$ to $t$. PT-post restricts untreated trends from 2 onward, so

$$
d_{14} = \mathrm{ATT}(3,4) + a_2-a_1, \qquad d_{12} = a_2-a_1.
$$

The two-row model and its target combination are

$$
\begin{pmatrix}d_{14}\\d_{12}\end{pmatrix} =
\begin{pmatrix}1&1\\0&1\end{pmatrix}
\begin{pmatrix}\mathrm{ATT}(3,4)\\a_2-a_1\end{pmatrix},
\qquad d_{14}-d_{12}=\mathrm{ATT}(3,4).
$$

Subtraction cancels the earlier differential without assuming that it is zero. This is a population identity for the same cohort and control populations. The two windows can contain different units, whose changes need only be representative in mean over their respective pairs. Their fitted comparisons need not telescope exactly in the sample. Catalog support, overlap and the observation/model assumptions remain required.

If the $\{1,2\}$ window is absent and the only available comparison is $d_{14}$, the effect and $a_2-a_1$ cannot be separated under this PT-post model. PT-all with `pt_start = 1` sets the earlier differential to zero, which makes the effect recoverable from $d_{14}$ if its other assumptions hold. A false restriction biases that recovery.

### Why the control pool must match

With `control_group = "notyettreated"`, the largest eligible control set can differ from one pair to another. A comparison over $\{1,2\}$ may use cohorts that are already treated by date 4, while a comparison over $\{1,4\}$ must exclude them. Their unrestricted early differentials are then different parameters and need not cancel.

The catalog therefore retains valid common, fixed-horizon control-pool reevaluations alongside every original valid comparison. For an earlier pair, a fixed-horizon reuse evaluates it with controls that remain eligible through the specified later horizon. This can provide the same control pool and reference law for the subtraction above. Only the retained, qualified combinations justify cancellation. The package does not assume that any connected sequence of different control pools supplies the needed equality, and clean control dates alone do not establish parallel trends or valid observation reconstruction.

## PT-all

Set `pt_assumption = "all"` and give `pt_start`, one finite numeric date that matches a date of the stored calendar exactly and is common to all cohorts. PT-all then:

- keeps every PT-post restriction, including at treated dates earlier than the start;
- keeps each cohort's original reference period;
- restricts the represented pre-reference trend differentials for each cohort and its catalog control operators between `pt_start` and the original reference period to zero;
- leaves trends before `pt_start` and anticipation effects unrestricted;
- uses the same comparisons, comparison groups and combination step as PT-post, and keeps every assumption about sampling, missing outcomes, overlap and models.

The added restrictions can make previously unavailable effects estimable, as in the example above. When the effects are already available, valid and informative restrictions can also improve precision. Neither gain is guaranteed: a particular standard error can move in either direction, and restrictions that are false bias the estimates. "All" refers to the dates from the chosen start onward, not to every historical date.

Algebraically, PT-all writes $\theta = H\eta$, where $H$ keeps the free parameters and sets only those represented in-window pre-reference differentials to zero. The model becomes $d = AH\eta$. The comparison catalog, population references and GMM combination remain the same. Earlier differentials stay free, and anticipation contributions remain free causal effects rather than imposed zeros. PT-all changes the restrictions within this estimator. It does not activate a separate efficiency estimator.

A comparison between two restricted pre-reference dates can then have a zero row in $AH$: its model mean is zero for every $\eta$. The observed comparison is still noisy and may be correlated with comparisons that load on the treatment effects. Under covariance-weighted joint GMM, that correlation can inform the effect combination. Such zero-loading rows remain in the moment system rather than being discarded.

`pt_start` is required under PT-all: omitting it, passing `NULL` or giving a date that is not in the stored calendar is an error, and no rounding or default date is used. Under PT-post, a non `NULL` `pt_start` is refused.

## When PT-all adds nothing

If no cohort has a represented step between `pt_start` and its reference period, PT-all adds no restriction. The fit still runs, the model is the PT-post model, and `summary()` reports "No additional restrictions". This happens, for example, when `pt_start` equals the reference period of every cohort, or when anticipation moves the reference period back to the start. The result records the count of added restrictions for each cohort, so you can see where the assumption has content.

## Choosing the start date

Choose `pt_start` from knowledge of the application, for example the earliest period from which untreated trends are believed to be parallel, before looking at estimates, standard errors or pretrend tests. A start selected by searching for precise or significant results invalidates the interpretation of the restrictions.

## Reporting with base_period

`base_period` controls how cells before treatment are reported, with `did`'s meaning and default. It does not change the parallel trends assumption, and `pt_assumption` does not change the reporting.

| `base_period` | Cells before treatment | Cells from treatment onward |
| --- | --- | --- |
| `"varying"` (default) | each stored period after the first compared with the preceding stored period | each period compared with the reference period |
| `"universal"` | each period compared with the reference period; the cell at the reference period is an exact zero | each period compared with the reference period |

Under both values and under both parallel trends choices, the cells before treatment are direct comparisons with no restriction imposed. Restricted zeros of PT-all never appear among them.

## Pretrend tests under each choice

The fit reports a Wald test of all declared nonreference pretrend cells of the chosen `base_period`, stored as `W` and `Wpval`. Exact reference zeros are identities and are excluded. The test is computed only when every tested cell is available and the full covariance passes the existing Cholesky/solve and conditioning policy; otherwise the test is unavailable and the reason is recorded. A test on a subset of the cells is never substituted.

Under PT-all, a cell that involves dates before `pt_start` can reject without contradicting PT-all from the chosen start, and cells in anticipation periods need not be zero. Descriptive tags mark each comparison before treatment as lying before the start, inside the restricted window or in the anticipation periods. An in-window tag implies a zero contrast only under that comparison's actual control operator and selection/trend premises; an auxiliary not-yet-treated pool absent from the catalog assumptions is not covered by its date tag alone. A test that does not reject does not establish parallel trends, and restricted zeros are never used as evidence for them.

`conditional_did_pretest()` runs a separate conditional test, described in the [results reference](../reference/results.md#pretrend-tests). Its arguments do not include the PT-all settings.

## Example

This example reuses `sim_panel` from [Get started](../get-started.md#prepare-the-panel); run that code first. Cohort 3 has reference period 2, and the simulation satisfies parallel trends at every date, so a restriction from period 1 is valid.


```r
fit_post <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  data = sim_panel, control_group = "nevertreated", est_method = "ipw",
  pt_assumption = "post",
  bstrap = FALSE, cband = FALSE
)

fit_all <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  data = sim_panel, control_group = "nevertreated", est_method = "ipw",
  pt_assumption = "all", pt_start = 1,
  bstrap = FALSE, cband = FALSE
)

fit_all_start2 <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  data = sim_panel, control_group = "nevertreated", est_method = "ipw",
  pt_assumption = "all", pt_start = 2,
  bstrap = FALSE, cband = FALSE
)

fit_all_universal <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  data = sim_panel, control_group = "nevertreated", est_method = "ipw",
  pt_assumption = "all", pt_start = 1, base_period = "universal",
  bstrap = FALSE, cband = FALSE
)

fit_all$chained$pt$by_cohort
fit_all_start2$chained$pt$no_extra_restrictions
summary(fit_all_start2)

compare <- merge(
  data.frame(group = fit_post$group, time = fit_post$t,
             post = fit_post$att, post_se = fit_post$se),
  data.frame(group = fit_all$group, time = fit_all$t,
             all = fit_all$att, all_se = fit_all$se),
  by = c("group", "time")
)
compare

cdid::tidy(fit_all_universal)[, c("target.id", "role", "estimate", "std.error")]
```

<details>
<summary>Console output</summary>

```text
g ref pt_start  baseline_status n_extra_calendar_dates n_extra_represented_dates
1 3   2        1 stored_reference                      1                         1
  n_restricted_pre_parameters n_zero_loading_links no_extra_calendar_window no_extra_restrictions
1                           1                    1                    FALSE                 FALSE
[1] TRUE

cdid chained difference-in-differences (chained IPW-GMM): cdid_MP
  cdid 1.1.0.0 result built on did 2.5.1 (tested: did 2.5.1); did's MP methods tabulate the estimates below.
  Chained estimators: Bellego C., Benatia D. and V. Dortet-Bernadet (2025). "The Chained Difference-in-Differences." Journal of Econometrics 248, 105783. doi:10.1016/j.jeconom.2024.105783.
  Ordinary estimators: Callaway B. and P. H. C. Sant'Anna (2021). "Difference-in-Differences with Multiple Time Periods." Journal of Econometrics 225(2), 200-230. doi:10.1016/j.jeconom.2020.12.001.

Call:
cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "first_treat", 
    data = sim_panel, control_group = "nevertreated", est_method = "ipw", 
    pt_assumption = "all", pt_start = 2, bstrap = FALSE, cband = FALSE)

Group-Time Average Treatment Effects:
 Group Time ATT(g,t) Std. Error [95% Pointwise  Conf. Band]  
     3    2  -0.0583     0.0637         -0.1832      0.0665  
     3    3   1.0096     0.0632          0.8858      1.1334 *
     3    4   1.3274     0.0616          1.2066      1.4482 *
     3    5   1.4837     0.0630          1.3603      1.6072 *
     3    6   1.8478     0.0630          1.7243      1.9713 *
---
Signif. codes: `*' confidence band does not cover 0

P-value for pre-test of parallel trends assumption:  0.359889309273502
Control Group:  Never Treated,  Anticipation Periods:  0
Estimation Method:  Inverse Probability Weighting

Chained support: 5 of 5 requested cells available.
  Units: 2000 in the unit table, 2000 observed at all 6 dates, 2000 with at least two observed outcomes.
  Observation correction: none (selection = "none").
  Parallel trends: PT-all from common stored date 2, in addition to every PT-post restriction; base_period = "varying".
    No additional restrictions (PT-post-equivalent model).
    Cohort 3: reference period 2; stored window dates 0; represented restrictions 0.
  GMM weighting: optimal (two-step); covariance regularization: not applied (moment covariance well conditioned).
  Inference: analytic standard errors, unit-level; pointwise intervals.
  Parallel-trends pre-test: W = 0.838279, df = 1, p = 0.3599 (Wald test of the pretrend cell).
  group time        post    post_se         all     all_se
1     3    2 -0.05831968 0.06369725 -0.05831968 0.06369725
2     3    3  1.00962257 0.06316956  0.98062848 0.05465855
3     3    4  1.32743963 0.06163755  1.29986751 0.05378014
4     3    5  1.48374276 0.06297661  1.45503400 0.05461554
5     3    6  1.84780436 0.06303371  1.81769972 0.05377841
           target.id      role   estimate  std.error
1       pre(g=3,t=1)       pre 0.05831968 0.06369725
2 reference(g=3,t=2) reference 0.00000000 0.00000000
3       ATT(g=3,t=3)      post 0.98062848 0.05465855
4       ATT(g=3,t=4)      post 1.29986751 0.05378014
5       ATT(g=3,t=5)      post 1.45503400 0.05461554
6       ATT(g=3,t=6)      post 1.81769972 0.05377841
```

</details>

**What the design implies.**

- `fit_all`: the stored calendar has one date, period 1, between the start and the reference period of cohort 3, so the record for cohort 3 should show one window date and the restriction it represents. The restriction is true in this simulation. The effects from period 3 onward may differ from those of `fit_post`, and so may their standard errors, in either direction for a particular draw.
- `fit_all_start2`: the start equals the reference period, so no restriction is added. `no_extra_restrictions` should be `TRUE`, the summary should report "No additional restrictions", and the estimates should equal those of `fit_post`.
- `fit_all_universal`: only the reporting of cells before treatment changes. The cell for period 1 is now measured relative to period 2, the cell at period 2 is an exact zero with role `"reference"`, and the effects from period 3 onward should equal those of `fit_all`.

## Where the result records the choice

| Location | Contents |
| --- | --- |
| `fit$chained$pt` | `pt_assumption` and `pt_start`; under PT-all also the totals `n_extra_calendar_dates` and `n_restricted_pre_parameters`, the flag `no_extra_restrictions`, the assumption text and the table `by_cohort` |
| `fit$chained$pt$by_cohort` | one row per cohort: `g`, its reference period `ref`, `pt_start`, `n_extra_calendar_dates`, `n_extra_represented_dates`, `n_restricted_pre_parameters`, `no_extra_calendar_window` and `no_extra_restrictions` |
| `fit$chained$base_period`, `fit$DIDparams$base_period` | the reporting choice |
| `summary(fit)` | the parallel trends assumption and `base_period`; under PT-all also the common start, the number of added restrictions (or "No additional restrictions") and a line for each cohort |

---

**Reading path, step 13 of 20** · Previous: [Impute actual changes](../examples/change-imputation.md) · Next: [Add PT-all restrictions](../examples/pt-all-comparison.md)

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
| 9 | [Handling missing outcomes](missing-outcomes.md) |
| 10 | [Sampling information for the observation model](sampling-information.md) |
| 11 | [Choose an observation model](../examples/observation-models.md) |
| 12 | [Impute actual changes](../examples/change-imputation.md) |
| 13 | **Parallel trends** (this page) |
| 14 | [Add PT-all restrictions](../examples/pt-all-comparison.md) |
| 15 | [GMM weighting and covariance regularization](gmm.md) |
| 16 | [Robustness](../examples/robustness.md) |
| 17 | [Inference](../examples/inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
