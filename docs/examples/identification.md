# Identify effects when the observed units change composition

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

This example builds a small rotating panel in which the treated units observed at each date have a different mix of a fixed, unobserved unit component. Outcome changes within the same unit still identify the treatment effects, while comparisons of mean outcome levels across dates do not. The reasoning is the argument about mean pair changes in [Estimate treatment effects from incomplete panels](incomplete-panels.md#telescoping-consecutive-windows-in-layout-2), applied to one concrete design. The example introduces no new estimator, option or assumption.



## Question

A study follows units over five dates. Units in cohort 3 are first treated at date 3, and the others are never treated. Each unit is observed at exactly two consecutive dates. Among treated units, those with a high fixed component tend to be observed early and those with a low one late. Write ATT(3,t) for the average treatment effect of cohort 3 at date t. Can the changes observed within units identify ATT(3,3), ATT(3,4) and ATT(3,5) when the observed treated units differ in composition from date to date?

## The design

- **Cohort.** $G$ is 0 or 3, each with probability 1/2.
- **Fixed component.** $\alpha$ is $-1$ or $+1$, each with probability 1/2, independently of $G$. The estimators never see $\alpha$, and no covariate is supplied.
- **Untreated outcomes.** $Y_t(0) = \alpha + 0.1\,t + \varepsilon_t$. The $\varepsilon_t$ have mean zero, are independent across dates and are independent of $G$, $\alpha$ and observation.
- **Effect.** Let $\tau_t = 1$ when $G = 3$ and $t \ge 3$, and $\tau_t = 0$ otherwise. The actual outcome is $Y_t = Y_t(0) + \tau_t$, so every treated unit gains exactly 1 from date 3 onward and nothing earlier. ATT(3,3), ATT(3,4) and ATT(3,5) all equal 1, for the simulated units and for the population the design describes.
- **Observation.** Each unit is observed in one window $W$ of two consecutive dates, (1,2), (2,3), (3,4) or (4,5), drawn given $G$ and $\alpha$ with the probabilities in Table 1.

**Table 1. Window probabilities given cohort and fixed component (population algebra)**

| Units | Share of units | (1,2) | (2,3) | (3,4) | (4,5) |
| --- | --- | --- | --- | --- | --- |
| $G = 0$, either value of $\alpha$ | 1/2 | 0.25 | 0.25 | 0.25 | 0.25 |
| $G = 3$, $\alpha = +1$ | 1/4 | 0.40 | 0.30 | 0.20 | 0.10 |
| $G = 3$, $\alpha = -1$ | 1/4 | 0.10 | 0.20 | 0.30 | 0.40 |
| $G = 3$, both values together | 1/2 | 0.25 | 0.25 | 0.25 | 0.25 |

Each window contains treated and never treated units with both values of $\alpha$ with positive probability, and each window holds a quarter of each cohort. No unit is observed at all five dates: every unit has exactly two observed outcomes.

Let $R_t = 1$ when date $t$ lies in the unit's window. Write $p_t^{+}$ and $p_t^{-}$ for the probability that $R_t = 1$ in cohort 3 when $\alpha = +1$ and when $\alpha = -1$; each is the sum of the probabilities of the windows that contain $t$. Because the two values of $\alpha$ have equal mass in cohort 3,

$$m_t = E[\alpha \mid G = 3, R_t = 1] = \frac{p_t^{+} - p_t^{-}}{p_t^{+} + p_t^{-}}.$$

**Table 2. Observation and the mean fixed component by date (population algebra)**

| Quantity | $t = 1$ | $t = 2$ | $t = 3$ | $t = 4$ | $t = 5$ |
| --- | --- | --- | --- | --- | --- |
| $p_t^{+}$ | 0.4 | 0.7 | 0.5 | 0.3 | 0.1 |
| $p_t^{-}$ | 0.1 | 0.3 | 0.5 | 0.7 | 0.4 |
| $P(R_t = 1 \mid G = 3) = P(R_t = 1 \mid G = 0)$ | 0.25 | 0.5 | 0.5 | 0.5 | 0.25 |
| $m_t = E[\alpha \mid G = 3, R_t = 1]$ | 0.6 | 0.4 | 0 | $-0.4$ | $-0.6$ |
| $E[\alpha \mid G = 0, R_t = 1]$ | 0 | 0 | 0 | 0 | 0 |

The two cohorts are observed with the same probability at every date, so cohort 3 makes up half of the observed units at each date. What changes is the mix within cohort 3: its observed units have a mean $\alpha$ of 0.4 at date 2 and $-0.4$ at date 4.

## Why changes within units identify the effects

For a unit observed in window $(s, s+1)$, the fixed component cancels from its change:

$$Y_{s+1} - Y_s = 0.1 + (\tau_{s+1} - \tau_s) + \varepsilon_{s+1} - \varepsilon_s.$$

The change depends on the cohort, the dates and the noise, while the window depends only on the cohort and $\alpha$. Within each cohort, the units observed at both dates of a pair therefore have the cohort's mean change,

$$E[Y_{s+1} - Y_s \mid G, W = (s, s+1)] = E[Y_{s+1} - Y_s \mid G],$$

even though they do not represent the cohort in levels. This mean condition is what `selection = "none"` requires, as described in [Handling missing outcomes](../methods/missing-outcomes.md#without-a-correction). Observation depends on $\alpha$, but $\alpha$ never enters a change, so no model of observation is needed.

Write $D_{s,s+1}$ for the mean change of cohort 3 minus the mean change of the never treated units. Every mean change equals 0.1, except the change of cohort 3 from date 2 to date 3, which is 1.1. Hence $D_{1,2} = 0$, $D_{2,3} = 1$, $D_{3,4} = 0$ and $D_{4,5} = 0$.

Untreated outcomes satisfy $Y_t(0) - Y_2(0) = 0.1\,(t - 2) + \varepsilon_t - \varepsilon_2$ in both cohorts, so PT-post from reference date 2 holds, and there is no anticipation. The first link then identifies the effect, and the later links identify its increments:

$$D_{2,3} = \text{ATT}(3,3) = 1, \qquad D_{3,4} = \text{ATT}(3,4) - \text{ATT}(3,3) = 0, \qquad D_{4,5} = \text{ATT}(3,5) - \text{ATT}(3,4) = 0.$$

Summing the links gives $\text{ATT}(3,4) = D_{2,3} + D_{3,4} = 1$ and $\text{ATT}(3,5) = D_{2,3} + D_{3,4} + D_{4,5} = 1$. The four windows are the only pairs anyone is observed at, so chained difference in differences (chained DiD) can only combine these consecutive comparisons. The cell at date 2, reported before treatment, uses $D_{1,2}$.

## Why comparisons of levels do not

The intercept-only repeated cross section (RCS) IPW estimator used below compares normalized means of outcome levels in each cohort and date arm. Its constant propensity odds cancel within each normalized control mean. In this design the population propensity is 1/2, away from the trimming threshold, so the population comparison for cell (3,t) with reference date 2 is

$$\bigl(E[Y_t \mid G = 3, R_t = 1] - E[Y_2 \mid G = 3, R_2 = 1]\bigr) - \bigl(E[Y_t \mid G = 0, R_t = 1] - E[Y_2 \mid G = 0, R_2 = 1]\bigr).$$

The common trend cancels, and the observed never treated units have a mean $\alpha$ of zero at every date, so for $t \ge 3$ this contrast equals $1 + m_t - m_2$, with $m_2 = 0.4$.

**Table 3. What each comparison targets (population algebra)**

| Cell | Known value | Chained links | Comparison of observed levels |
| --- | --- | --- | --- |
| Date 2, before treatment (dates 1 and 2) | 0 | $D_{1,2} = 0$ | $m_2 - m_1 = 0.4 - 0.6 = -0.2$ |
| ATT(3,3) | 1 | $D_{2,3} = 1$ | $1 + (0 - 0.4) = 0.6$ |
| ATT(3,4) | 1 | $1 + 0 = 1$ | $1 + (-0.4 - 0.4) = 0.2$ |
| ATT(3,5) | 1 | $1 + 0 + 0 = 1$ | $1 + (-0.6 - 0.4) = 0$ |

Callaway and Sant'Anna's repeated cross section setting assumes that the sample at each date represents the population at that date and that the distribution of cohorts and covariates is stable across dates. Here population parallel trends hold and the cohort shares are stable, but the first condition fails for cohort 3. Its observed units come mostly from $\alpha = +1$ at early dates and mostly from $\alpha = -1$ at late dates, so each level contrast adds a change in composition to the effect.

Fixed differences between units do not cause this on their own. The never treated units take the same two values of $\alpha$, but their windows do not depend on it, so their observed levels represent their cohort at every date. If the windows of cohort 3 were also drawn independently of $\alpha$, as in the rotating example of [Estimate treatment effects from incomplete panels](incomplete-panels.md), the observed levels would represent both cohorts, and the RCS approach would target the same effects under its assumptions. `did` handles unbalanced panels with that approach; whether it answers the question depends on its sampling conditions.

## Assumptions in this example

| Item | In this example |
| --- | --- |
| Whose effect | ATT(3,t) for t = 3, 4, 5 in the study population of 2000 simulated units listed in the unit table `units` |
| Comparison group | Units never treated, `g == 0` |
| Anticipation | None, so the reference period of cohort 3 is date 2 |
| Parallel trends | PT-post from date 2, true in the design |
| Observed changes | `selection = "none"`: within each cohort, units observed at both dates of a pair have the cohort's mean change; true because observation depends only on the cohort and $\alpha$ |
| Treatment adjustment | Treatment weighting by inverse probability weighting (IPW) with an intercept only, `xformla = ~ 1`; adequate because parallel trends hold without covariates |
| Overlap and support | Every window has positive probability in both cohorts; the script stops if a realized window lacks either cohort |
| Unit table | One row per unit with `id` and `g`; no covariates or weights |

## Code



The script simulates the design, fits chained DiD and fits `did` with the same upstream arguments. For this genuinely unbalanced panel, with `panel = TRUE` and `allow_unbalanced_panel = TRUE`, `did` 2.5.1 keeps rows with complete values of the variables it uses and applies its RCS point estimator, comparing levels among the units observed at each date. It retains the unit identity distinction for unbalanced-panel inference; using an RCS point construction does not mean it discards identities ([Moving from did](../from-did.md#incomplete-panels-in-did-and-in-chained-mode)). That is the route compared here. In chained mode, `allow_unbalanced_panel = TRUE` has the same meaning as the default `NULL`.

```r
set.seed(20261024)

N <- 2000L
T_dates <- 5L

# Cohort and fixed component: independent draws, each value with probability 1/2.
# alpha only generates the data; neither estimator receives it.
g     <- sample(c(0, 3), size = N, replace = TRUE)
alpha <- sample(c(-1, 1), size = N, replace = TRUE)

# Window probabilities of Table 1; window k covers dates k and k + 1.
p_window <- rbind(
  never      = c(0.25, 0.25, 0.25, 0.25),
  treat_high = c(0.40, 0.30, 0.20, 0.10),
  treat_low  = c(0.10, 0.20, 0.30, 0.40)
)
type <- ifelse(g == 0, "never", ifelse(alpha > 0, "treat_high", "treat_low"))
window_start <- integer(N)
for (k in rownames(p_window)) {
  in_k <- which(type == k)
  window_start[in_k] <- sample.int(4L, size = length(in_k), replace = TRUE,
                                   prob = p_window[k, ])
}

# Stop, without redrawing, if a realized window lacks either cohort.
stopifnot(all(table(factor(g, levels = c(0, 3)),
                    factor(window_start, levels = 1:4)) > 0))

# Full grid of potential outcomes for every unit and date, with independent noise.
full_grid <- data.frame(
  id   = rep(seq_len(N), each = T_dates),
  time = rep(seq_len(T_dates), times = N)
)
full_grid$g     <- g[full_grid$id]
full_grid$alpha <- alpha[full_grid$id]
full_grid$noise <- stats::rnorm(nrow(full_grid), sd = 0.5)
full_grid$y0    <- full_grid$alpha + 0.1 * full_grid$time + full_grid$noise
full_grid$y1    <- full_grid$y0 + 1
full_grid$y     <- full_grid$alpha + 0.1 * full_grid$time +
  1 * (full_grid$g > 0 & full_grid$time >= full_grid$g) + full_grid$noise

# Keep every unit and date row; the outcome is NA outside the unit's window.
in_window <- full_grid$time == window_start[full_grid$id] |
  full_grid$time == window_start[full_grid$id] + 1L
panel_data <- full_grid[, c("id", "time", "g")]
panel_data$y <- ifelse(in_window, full_grid$y, NA_real_)
stopifnot(all(tapply(!is.na(panel_data$y), panel_data$id, sum) == 2))

# Unit table: one row per study unit, with its identifier and cohort only.
units <- data.frame(id = seq_len(N), g = g)

# Chained DiD: IPW with an intercept only, never treated controls, PT-post,
# and no model of observation.
fit_chained <- cdid::att_gt(
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ 1, data = panel_data,
  panel = TRUE, allow_unbalanced_panel = TRUE,
  control_group = "nevertreated", anticipation = 0,
  est_method = "ipw", base_period = "varying",
  bstrap = FALSE, cband = FALSE,
  frame = units, selection = "none", pt_assumption = "post"
)

# did's RCS route: the same upstream arguments, no unit table, no chained arguments.
fit_rcs <- did::att_gt(
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ 1, data = panel_data,
  panel = TRUE, allow_unbalanced_panel = TRUE,
  control_group = "nevertreated", anticipation = 0,
  est_method = "ipw", base_period = "varying",
  bstrap = FALSE, cband = FALSE
)

summary(fit_chained)
cdid::tidy(fit_chained)[, c("target.id", "role", "estimate", "std.error",
                            "available", "support.status")]
agg_chained <- cdid::aggte(fit_chained, type = "dynamic")
summary(agg_chained)

summary(fit_rcs)
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
    xformla = ~1, data = panel_data, panel = TRUE, allow_unbalanced_panel = TRUE, 
    control_group = "nevertreated", anticipation = 0, est_method = "ipw", 
    base_period = "varying", bstrap = FALSE, cband = FALSE, frame = units, 
    selection = "none", pt_assumption = "post")

Group-Time Average Treatment Effects:
 Group Time ATT(g,t) Std. Error [95% Pointwise  Conf. Band]  
     3    2   0.0129     0.0640         -0.1126      0.1383  
     3    3   0.8972     0.0656          0.7686      1.0258 *
     3    4   0.8883     0.0894          0.7130      1.0635 *
     3    5   0.8679     0.1098          0.6527      1.0830 *
---
Signif. codes: `*' confidence band does not cover 0

P-value for pre-test of parallel trends assumption:  0.840754265674953
Control Group:  Never Treated,  Anticipation Periods:  0
Estimation Method:  Inverse Probability Weighting

Chained support: 4 of 4 requested cells available.
  Units: 2000 in the unit table, 0 observed at all 5 dates, 2000 with at least two observed outcomes.
  Observation correction: none (selection = "none").
  Parallel trends: PT-post; base_period = "varying".
  GMM weighting: optimal (two-step); covariance regularization: not applied (moment covariance well conditioned).
  Inference: analytic standard errors, unit-level; pointwise intervals.
  Parallel-trends pre-test: W = 0.0403724, df = 1, p = 0.8408 (Wald test of the pretrend cell).
     target.id role   estimate  std.error available            support.status
1 pre(g=3,t=2)  pre 0.01286068 0.06400615      TRUE supported_direct_pretrend
2 ATT(g=3,t=3) post 0.89721369 0.06559575      TRUE                 estimable
3 ATT(g=3,t=4) post 0.88827272 0.08942485      TRUE                 estimable
4 ATT(g=3,t=5) post 0.86786563 0.10978623      TRUE                 estimable

cdid chained difference-in-differences (chained IPW-GMM): cdid_AGGTEobj
  cdid 1.1.0.0 result built on did 2.5.1 (tested: did 2.5.1); did's AGGTEobj methods tabulate the estimates below.
  Chained estimators: Bellego C., Benatia D. and V. Dortet-Bernadet (2025). "The Chained Difference-in-Differences." Journal of Econometrics 248, 105783. doi:10.1016/j.jeconom.2024.105783.
  Ordinary estimators: Callaway B. and P. H. C. Sant'Anna (2021). "Difference-in-Differences with Multiple Time Periods." Journal of Econometrics 225(2), 200-230. doi:10.1016/j.jeconom.2020.12.001.

Call:
cdid::aggte(fit_chained, type = "dynamic")

Overall summary of ATT's based on event-study/dynamic aggregation:  
    ATT    Std. Error     [ 95%  Conf. Int.]  
 0.8845          0.08     0.7277      1.0412 *


Dynamic Effects:
 Event time Estimate Std. Error [95% Pointwise  Conf. Band]  
         -1   0.0129     0.0640         -0.1126      0.1383  
          0   0.8972     0.0656          0.7686      1.0258 *
          1   0.8883     0.0894          0.7130      1.0635 *
          2   0.8679     0.1098          0.6527      1.0830 *
---
Signif. codes: `*' confidence band does not cover 0

Control Group:  Never Treated,  Anticipation Periods:  0
Estimation Method:  Inverse Probability Weighting

Chained dynamic aggregation: 4 of 4 intended aggregate(s) available; na.rm = FALSE.
  Overall summary: available over its 3 intended bins; pointwise interval.
  Parallel trends: PT-post; base_period = "varying".
  Inference: analytic standard errors, unit-level; pointwise intervals.

Call:
did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g", 
    xformla = ~1, data = panel_data, panel = TRUE, allow_unbalanced_panel = TRUE, 
    control_group = "nevertreated", anticipation = 0, bstrap = FALSE, 
    cband = FALSE, est_method = "ipw", base_period = "varying")

Reference: Callaway, Brantly and Pedro H.C. Sant'Anna.  "Difference-in-Differences with Multiple Time Periods." Journal of Econometrics, Vol. 225, No. 2, pp. 200-230, 2021. <https://doi.org/10.1016/j.jeconom.2020.12.001>, <https://arxiv.org/abs/1803.09015> 

Group-Time Average Treatment Effects:
 Group Time ATT(g,t) Std. Error [95% Pointwise  Conf. Band]  
     3    2  -0.1196     0.0819         -0.2801      0.0408  
     3    3   0.6840     0.0789          0.5293      0.8387 *
     3    4   0.2791     0.0962          0.0906      0.4676 *
     3    5  -0.0625     0.1148         -0.2875      0.1625  
---
Signif. codes: `*' confidence band does not cover 0

P-value for pre-test of parallel trends assumption:  0.14382
Control Group:  Never Treated,  Anticipation Periods:  0
Estimation Method:  Inverse Probability Weighting
```

</details>

The ordinary RCS fit reports that it drops 6,000 rows with missing outcomes. This is expected here: it compares the observed levels, whereas the chained fit uses the available within-unit changes and the supplied unit table.

Optionally, put both sets of cells beside the known values of Table 3:

```r
known <- data.frame(
  group = 3, time = 2:5,
  known_value    = c(0, 1, 1, 1),
  level_contrast = c(-0.2, 0.6, 0.2, 0)
)
estimates <- merge(
  data.frame(group = fit_chained$group, time = fit_chained$t,
             chained = fit_chained$att, chained_se = fit_chained$se),
  data.frame(group = fit_rcs$group, time = fit_rcs$t,
             rcs = fit_rcs$att, rcs_se = fit_rcs$se),
  by = c("group", "time")
)
merge(known, estimates, by = c("group", "time"))
```

Output:

```text
group time known_value level_contrast    chained chained_se         rcs     rcs_se
1     3    2           0           -0.2 0.01286068 0.06400615 -0.11964856 0.08185521
2     3    3           1            0.6 0.89721369 0.06559575  0.68400574 0.07892854
3     3    4           1            0.2 0.88827272 0.08942485  0.27907751 0.09616598
4     3    5           1            0.0 0.86786563 0.10978623 -0.06246468 0.11479458
```
| Argument | Setting | Definition |
| --- | --- | --- |
| `est_method`, `xformla` | `"ipw"`, `~ 1` | [Treatment adjustment](../reference/options.md#treatment-adjustment) |
| `control_group`, `anticipation`, `panel`, `allow_unbalanced_panel` | `"nevertreated"`, 0, `TRUE`, `TRUE` | [Data, timing and comparison group](../reference/options.md#data-timing-and-comparison-group) |
| `base_period`, `pt_assumption` | `"varying"`, `"post"` | [Parallel trends and reporting](../reference/options.md#parallel-trends-and-reporting) |
| `frame`, `selection` | `units`, `"none"` | [Chained arguments](../reference/options.md#chained-arguments) |
| `bstrap`, `cband` | `FALSE`, `FALSE` | [Inference and computation](../reference/options.md#inference-and-computation) |

The fitted result records these choices: `fit_chained$chained$selection$method` for the handling of missing outcomes and `fit_chained$chained$pt` for the parallel trends choice.

## Known values and future estimates

The design fixes the known values in Table 3: 1 for each effect and 0 for the comparison before treatment. Under the assumptions above, the cells of `fit_chained` target these values. With one treated cohort, `agg_chained` reports the same cells by event time, so event times 0, 1 and 2 have known value 1 and so does their overall average. The cells of `fit_rcs` target the level contrasts 0.6, 0.2 and 0 after treatment and $-0.2$ before it.

Estimates from a future run need not equal these population values. The algebra does not establish fit completion or the validity of reported uncertainty; [What the uncertainty assumes](../reference/results.md#what-the-uncertainty-assumes) states its separate large sample requirements. In particular, an RCS standard error cannot repair the composition shift relative to the causal ATT. One draw does not measure bias or coverage, and a comparison of precision for a common target would not be justified by these two outputs alone.

A fit can still fail in a particular sample, for example when a model has no finite solution, and a cell can be unavailable, reported as `NA` with a reason ([Unavailable effects](../reference/results.md#unavailable-effects)). Keep any such error or status as it is: the script does not retry with other settings, redraw the sample or drop units. [Diagnose errors and unavailable results](../help/troubleshooting.md) describes the checks.

## Why there is no complete case comparison

`did`'s default, `allow_unbalanced_panel = FALSE`, keeps only units observed at every date. Every unit here has exactly two observed outcomes, so that subset is empty and no complete case estimate exists. The script therefore contains no such call. In chained mode, an explicit `allow_unbalanced_panel = FALSE` is refused with class `cdid_chained_balanced`.

## What the example shows

The construction shows one way in which the two routes can target different quantities from the same data. It does not show that chained DiD is preferable in general. In other designs both routes can be valid, and either can be more precise. Chained DiD with `selection = "none"` can fail when observation changes the cohort mean of an observed pair change, for example through selection on noise or heterogeneous effects. Dependence alone does not establish that failure: the required mean representativeness may still hold. If $\alpha$ were measured and supplied as a covariate, models that condition on it would become available, each with its own conditions, and the comparison could come out differently. Here $\alpha$ is unobserved and the formula is `~ 1`.

## Next step

When observation depends on measured characteristics that also shift the changes, [Handling missing outcomes](../methods/missing-outcomes.md) explains the conditional assumption and the three ways to adjust for it.


## Output from the example

| Period | Target truth | Chained estimate | SE | RCS estimate | SE | RCS population contrast |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 0.000 | 0.013 | 0.064 | -0.120 | 0.082 | -0.200 |
| 3 | 1.000 | 0.897 | 0.066 | 0.684 | 0.079 | 0.600 |
| 4 | 1.000 | 0.888 | 0.089 | 0.279 | 0.096 | 0.200 |
| 5 | 1.000 | 0.868 | 0.110 | -0.062 | 0.115 | 0.000 |

Each unit has exactly two observed outcomes. The post-treatment truth is 1 at every date; the RCS population contrasts are 0.6, 0.2 and 0.0. The estimates above are one draw. `did` warned that it removed 6000 rows with missing outcomes; that warning is expected for this RCS call and was retained.

---

**Reading path, step 4 of 20** · Previous: [Incomplete panels](incomplete-panels.md) · Next: [How to compare precision](precision.md)

<details>
<summary>All steps of the reading path</summary>

| Step | Page |
| --- | --- |
| 1 | [Get started](../get-started.md) |
| 2 | [Moving from did](../from-did.md) |
| 3 | [Incomplete panels](incomplete-panels.md) |
| 4 | **Identify effects when the observed units change composition** (this page) |
| 5 | [How to compare precision](precision.md) |
| 6 | [Prepare your panel and unit table](../help/prepare-data.md) |
| 7 | [Treatment adjustment](../methods/treatment.md) |
| 8 | [Choose IPW, OR or DR in an incomplete panel](treatment-choice.md) |
| 9 | [Handling missing outcomes](../methods/missing-outcomes.md) |
| 10 | [Sampling information for the observation model](../methods/sampling-information.md) |
| 11 | [Choose an observation model](observation-models.md) |
| 12 | [Impute actual changes](change-imputation.md) |
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
