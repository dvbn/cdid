# Estimate treatment effects from incomplete panels

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

Many longitudinal studies track the same units over time but do not observe every unit at every calendar date. Survey rotation designs observe each wave for a limited window, administrative registers miss periods during temporary gaps, and panel studies experience late entry or dropout. `cdid` implements chained difference in differences (Bellégo, Benatia and Dortet-Bernadet 2025), extending the group and time framework of Callaway and Sant'Anna (2021) and the workflow of `did` so you can estimate treatment effects from incomplete panels without restricting the analysis to complete cases or discarding within unit links.

This guide shows how chained DiD combines short outcome changes across observation patterns, states the assumptions required for those combinations to identify the target effect, and walks through a complete rotating panel script in base R. For an initial balanced panel walkthrough, see [Get started](../get-started.md) and [Coming from `did`](../from-did.md).

For your own panel, [Prepare your panel and unit table](../help/prepare-data.md) states the minimum inputs, and [Investigate missingness and data structure](../help/investigate-missingness.md) shows counts and selected-pair checks that include units observed once or never.

## Question

Suppose a study population contains units first treated at period $G = 3$ (treatment cohort $G = 3$) and units never treated ($G = 0$) across five calendar dates $t \in \{1, 2, 3, 4, 5\}$. With no treatment anticipation (`anticipation = 0`), period $t = 2$ is the baseline reference period for cohort $G = 3$. How can we identify and estimate the average treatment effect for cohort $G = 3$ at period $t = 5$, denoted $\text{ATT}(3, 5)$, when no single unit in the study population is observed at both the reference period $t = 2$ and the target period $t = 5$?

## Usable pair changes versus raw level comparisons

An incomplete panel differs fundamentally from genuine repeated cross sections. When repeated observations on the same unit are linked by `idname`, any unit observed at two dates $s < t$ provides an observed outcome change

$$\Delta Y_{s,t} = Y_t - Y_s$$

Taking the difference between period $t$ and period $s$ within the same unit removes any additive unobserved individual component $\alpha_i$ that is constant over time.

Treating an incomplete panel as if it were repeated cross sections compares raw outcome levels across dates without linking observations from the same unit. When entry, dropout, or rotation changes the composition of unobserved fixed effects across dates, raw level samples at different dates can fail to represent the same study population even when within unit changes satisfy the assumptions for chained DiD. Repeated cross section estimators require additional representativeness conditions on outcome levels and stable sample composition across dates that assumptions about outcome changes do not supply.

At the same time, restricting a panel to complete cases (units observed at all five dates) discards every unit with missing outcomes and fails entirely in a rotating design where no unit is observed in every period. Chained DiD instead forms comparisons from units observed at both dates of each usable pair $\{s, t\}$ and chains those comparisons together across time. Because a unit must have at least two observed outcomes to supply a change, units with repeated observations contribute outcome differences, while units with one or no observed outcomes remain in the unit table (`frame`) to represent the study population.

## Four observation patterns across five dates

To see how observation patterns determine which pair comparisons are available, consider five calendar dates $t = 1, 2, 3, 4, 5$, one treated cohort first treated at $t = 3$ ($G = 3$), and a comparison group of units never treated ($G = 0$).

**Legend for the four observation tables:**
- Each row represents an **observation group** of study units that mixes **both** treated units in cohort $G = 3$ and units never treated in $G = 0$.
- `Observed` means the outcome $Y_{i,t}$ is recorded for units in that observation group at date $t$.
- `Missing` means the outcome $Y_{i,t}$ is missing (`NA`) for units in that observation group at date $t$.

### Layout 1: Balanced complete panel

In a balanced panel, a single observation group covering the entire study population is observed at all five dates.

| Observation group (each row mixes cohort `G = 3` and never treated controls `G = 0`) | Observed dates | `t = 1` | `t = 2` | `t = 3` | `t = 4` | `t = 5` |
| --- | --- | --- | --- | --- | --- | --- |
| Group 1 (`G = 3` and `G = 0` units) | `{1, 2, 3, 4, 5}` | Observed | Observed | Observed | Observed | Observed |

Every unit is observed at both endpoints of every date pair $\{s, t\}$. Direct comparisons from the reference period $t = 2$ to $t = 3$, $t = 4$, and $t = 5$ are all directly available on the full study population, alongside every adjacent and multiperiod pair.

### Layout 2: Rotating consecutive windows

In a two period rotating panel, the study population is partitioned into four observation groups. Each unit is observed for exactly one window of two consecutive dates: $\{1, 2\}$, $\{2, 3\}$, $\{3, 4\}$, or $\{4, 5\}$.

| Observation group (each row mixes cohort `G = 3` and never treated controls `G = 0`) | Observed dates | `t = 1` | `t = 2` | `t = 3` | `t = 4` | `t = 5` |
| --- | --- | --- | --- | --- | --- | --- |
| Window group 1 (`G = 3` and `G = 0` units) | `{1, 2}` | Observed | Observed | Missing | Missing | Missing |
| Window group 2 (`G = 3` and `G = 0` units) | `{2, 3}` | Missing | Observed | Observed | Missing | Missing |
| Window group 3 (`G = 3` and `G = 0` units) | `{3, 4}` | Missing | Missing | Observed | Observed | Missing |
| Window group 4 (`G = 3` and `G = 0` units) | `{4, 5}` | Missing | Missing | Missing | Observed | Observed |

No unit is observed at both $t = 2$ and $t = 4$, and no unit is observed at both $t = 2$ and $t = 5$. Ordinary two period difference in differences between the reference date $t = 2$ and $t = 5$ has zero units observed at both dates. However, Window group 2 observes $\{2, 3\}$, Window group 3 observes $\{3, 4\}$, and Window group 4 observes $\{4, 5\}$, forming a consecutive chain from $t = 2$ to $t = 5$.

### Layout 3: General overlapping pairs

Incomplete panels in practice rarely follow a single consecutive chain. Units may enter early, experience intermittent gaps, and reenter later, generating both adjacent pairs and longer spans.

| Observation group (each row mixes cohort `G = 3` and never treated controls `G = 0`) | Observed dates | `t = 1` | `t = 2` | `t = 3` | `t = 4` | `t = 5` |
| --- | --- | --- | --- | --- | --- | --- |
| Pair group 1 (`G = 3` and `G = 0` units) | `{1, 2}` | Observed | Observed | Missing | Missing | Missing |
| Pair group 2 (`G = 3` and `G = 0` units) | `{1, 4}` | Observed | Missing | Missing | Observed | Missing |
| Pair group 3 (`G = 3` and `G = 0` units) | `{2, 3}` | Missing | Observed | Observed | Missing | Missing |
| Pair group 4 (`G = 3` and `G = 0` units) | `{3, 5}` | Missing | Missing | Observed | Missing | Observed |
| Pair group 5 (`G = 3` and `G = 0` units) | `{4, 5}` | Missing | Missing | Missing | Observed | Observed |

Here multiple routes connect the reference period $t = 2$ to $t = 5$:
- A two step posttreatment route combines $\{2, 3\}$ and $\{3, 5\}$.
- An alternative route bridges backward from $t = 2$ to $t = 1$ using $\{1, 2\}$, jumps to $t = 4$ using $\{1, 4\}$, and steps to $t = 5$ using $\{4, 5\}$.

`cdid` does not discard observations outside a single handpicked path or collapse the comparison catalog to one chain. Instead, the package retains all valid pair comparisons in the catalog and combines the available routes through generalized method of moments (GMM), accounting for the correlation across comparisons when units or estimated model parameters are shared.

### Layout 4: Gapped panel with an earlier bridge and a missing link

In a gapped panel, the adjacent transition $\{2, 3\}$ from the reference period to the first treatment period may not be directly observed, while earlier pairs $\{1, 2\}$ and $\{1, 3\}$ and a later span $\{3, 5\}$ are observed.

| Observation group (each row mixes cohort `G = 3` and never treated controls `G = 0`) | Observed dates | `t = 1` | `t = 2` | `t = 3` | `t = 4` | `t = 5` |
| --- | --- | --- | --- | --- | --- | --- |
| Pair group 1 (`G = 3` and `G = 0` units) | `{1, 2}` | Observed | Observed | Missing | Missing | Missing |
| Pair group 2 (`G = 3` and `G = 0` units) | `{1, 3}` | Observed | Missing | Observed | Missing | Missing |
| Pair group 3 (`G = 3` and `G = 0` units) | `{3, 5}` | Missing | Missing | Observed | Missing | Observed |

This layout illustrates two important features of chained DiD:
1. Even though no unit is observed at both $t = 2$ and $t = 3$, the earlier pairs $\{1, 2\}$ and $\{1, 3\}$ form a valid bridge that recovers the comparison between $t = 2$ and $t = 3$, which then combines with $\{3, 5\}$ to identify $\text{ATT}(3, 5)$.
2. No observed pair links $t = 4$ to the reference period. These comparisons therefore do not recover $\text{ATT}(3, 4)$ under the stated baseline restrictions. If date 4 is retained in `data` with `NA` outcomes, the requested cell is unavailable; if every date-4 row is omitted, 4 is absent from the stored calendar and no date-4 cell is returned.

## How chained differences identify treatment effects

### Population difference in mean changes

For any ordered date pair $s < t$, let $D_{s,t}$ denote the population mean outcome change between dates $s$ and $t$ among treated cohort $G = 3$ minus the population mean outcome change between dates $s$ and $t$ among the same never treated comparison population ($G = 0$):

$$D_{s,t} = E[Y_t - Y_s \mid G = 3] - E[Y_t - Y_s \mid G = 0]$$

Notice that $D_{s,t}$ is defined on the entire study population of cohort $G = 3$ and never treated controls $G = 0$, not on changing subpopulations across dates.

### Parallel trends from the reference period (PT-post)

Under parallel trends from each cohort's reference period (`PT-post`, selected by `pt_assumption = "post"`), cohort $G = 3$ with zero anticipation (`anticipation = 0`) has reference period $t = 2$. `PT-post` restricts untreated potential outcome trends from period $t = 2$ onward:

$$E[Y_t(0) - Y_2(0) \mid G = 3] = E[Y_t(0) - Y_2(0) \mid G = 0] \quad \text{for } t \in \{3, 4, 5\}$$

No parallel trends restriction is assumed before period $t = 2$ under `PT-post`. In particular, the untreated trend differential between $t = 1$ and $t = 2$,

$$D_{1,2} = E[Y_2(0) - Y_1(0) \mid G = 3] - E[Y_2(0) - Y_1(0) \mid G = 0]$$

is unrestricted and need not equal zero.

### Telescoping consecutive windows in Layout 2

Under `PT-post` and no anticipation before $t = 3$, the adjacent population differences from $t = 2$ onward satisfy three population identities:

1. Between $t = 2$ and $t = 3$, parallel trends between $t = 2$ and $t = 3$ implies that untreated changes cancel, so

$$D_{2,3} = \text{ATT}(3, 3)$$

2. Between $t = 3$ and $t = 4$, writing $Y_4(0) - Y_3(0) = \bigl(Y_4(0) - Y_2(0)\bigr) - \bigl(Y_3(0) - Y_2(0)\bigr)$ shows that `PT-post` implies $E[Y_4(0) - Y_3(0) \mid G = 3] = E[Y_4(0) - Y_3(0) \mid G = 0]$. The untreated change between $t = 3$ and $t = 4$ cancels, leaving the increment in treatment effects:

$$D_{3,4} = \text{ATT}(3, 4) - \text{ATT}(3, 3)$$

3. Between $t = 4$ and $t = 5$, by the same logic,

$$D_{4,5} = \text{ATT}(3, 5) - \text{ATT}(3, 4)$$

Summing the consecutive population differences telescopes the intermediate effects:

$$\text{ATT}(3, 4) = D_{2,3} + D_{3,4}$$

$$\text{ATT}(3, 5) = D_{2,3} + D_{3,4} + D_{4,5}$$

Under baseline mean missingness restrictions (`selection = "none"`), within each cohort the units observed in each of Window groups 2 ($\{2, 3\}$), 3 ($\{3, 4\}$), and 4 ($\{4, 5\}$) must be representative of that cohort with respect to its mean change over that pair. The comparisons then reconstruct $D_{2,3}$, $D_{3,4}$, and $D_{4,5}$ for the **same** cohort $G = 3$ and never treated $G = 0$ populations. Changing composition that changes these means can invalidate that reconstruction. These equations equate population expectations. They require neither identical unit-level changes nor exact cancellation among estimated sample comparisons, and they promise neither finite sample unbiasedness nor exact nominal coverage.

### Earlier bridge cancellation in Layout 4 without PT-all

Now consider Layout 4, where only $\{1, 2\}$, $\{1, 3\}$, and $\{3, 5\}$ are observed. Because $\{2, 3\}$ is not observed directly, can we still identify $\text{ATT}(3, 5)$ under `PT-post` without assuming parallel trends between $t = 1$ and $t = 2$?

Yes. In the study population, the identity $Y_3 - Y_1 = (Y_3 - Y_2) + (Y_2 - Y_1)$ holds for every unit. Because $D_{1,3}$ and $D_{1,2}$ are defined as expectations over the **exact same** cohort $G = 3$ and never treated $G = 0$ populations, linearity of expectation over that fixed population gives

$$D_{1,3} - D_{1,2} = \bigl(E[Y_3 - Y_1 \mid G = 3] - E[Y_3 - Y_1 \mid G = 0]\bigr) - \bigl(E[Y_2 - Y_1 \mid G = 3] - E[Y_2 - Y_1 \mid G = 0]\bigr)$$

$$= E[Y_3 - Y_2 \mid G = 3] - E[Y_3 - Y_2 \mid G = 0] = D_{2,3}$$

In terms of potential outcomes, $D_{1,3} = \text{ATT}(3, 3) + D_{1,2}$ under `PT-post` (since the untreated trend between $t = 2$ and $t = 3$ is zero while the earlier trend differential $D_{1,2}$ between $t = 1$ and $t = 2$ is unrestricted). Subtracting $D_{1,2}$ removes the unrestricted prebaseline trend differential $D_{1,2}$ exactly:

$$D_{1,3} - D_{1,2} = \text{ATT}(3, 3)$$

Combining this bridge with $D_{3,5} = \text{ATT}(3, 5) - \text{ATT}(3, 3)$ (which holds under `PT-post` from $t = 2$ onward) yields

$$D_{1,3} - D_{1,2} + D_{3,5} = \text{ATT}(3, 5)$$

This cancellation works because both $D_{1,3}$ and $D_{1,2}$ represent the **same** cohort $G = 3$ and control $G = 0$ study populations. `cdid` under `PT-post` does not impose an artificial endpoint cutoff at $t = 2$ that throws away pairs touching $t = 1$, and using the bridge $D_{1,3} - D_{1,2}$ does not invoke additional parallel trends restrictions (`PT-all`, `pt_assumption = "all"`). Instead, `PT-post` retains earlier bridges whenever their prebaseline trend differentials cancel across comparisons representing the same population. By contrast, if $\{1, 2\}$ were also missing so only $\{1, 3\}$ and $\{3, 5\}$ were observed, $D_{1,2}$ would not cancel under `PT-post`; identifying $\text{ATT}(3, 5)$ from $\{1, 3\}$ and $\{3, 5\}$ alone needs an additional identifying restriction. One sufficient option is `PT-all` starting at $t = 1$ (`pt_assumption = "all", pt_start = 1`), which restricts $D_{1,2} = 0$. For a complete comparison of `PT-post` and `PT-all`, see [Parallel trends assumptions](../methods/parallel-trends.md).

### Missing links, necessary support, and unavailable effects

Under the baseline model illustrated here, a target needs catalog comparisons whose unrestricted trend terms cancel. If a requested cell cannot be recovered, `cdid` records catalog-relative unavailability and its reason in `result$chained`; this is not a blanket population nonidentification verdict.

In Layout 4, the observed pairs do not recover date 4. If `{3, 5}` is also removed, the earlier bridge still recovers the date-3 effect, but these remaining pairs do not recover dates 4 or 5. A date is reported only if it remains in the stored calendar: keep its rows with `y = NA` when it is part of the intended study calendar, rather than assuming an omitted date will appear as an `NA` result.

Keep two practical distinctions in mind when inspecting unavailable effects or fit errors:
1. **Observation graph connectivity alone is not sufficient.** For the baseline comparisons here, a path must connect the target and reference dates, but connectivity by itself does not guarantee identification or estimation. Each used pair $\{s, t\}$ must also have sufficient sample support (both cohort $G = 3$ units and eligible control units observed at both $s$ and $t$), covariate overlap, and valid observation selection and parallel trends assumptions so that the observed pair represents the target and control populations.
2. **Not every model or solver refusal means lack of population identification.** When a fit stops because a treatment propensity logit or observation probability logit encounters separation or boundary probabilities, or when an outcome regression design is rank deficient on the realized sample, `cdid` raises a typed error describing the failed model component. That failure reflects the fitted model and sample support, not a mathematical proof that the population effect is unidentified.

## Assumptions for the rotating panel example

Before running the code below, review the assumptions required for this baseline rotating panel analysis:
- **Study population and unit table:** The target parameter is $\text{ATT}(3, t)$ for $t \in \{3, 4, 5\}$ in the study population represented by the 1,000 units in `units` (`frame`).
- **Comparison group and no anticipation:** Units never treated (`g == 0`, `control_group = "nevertreated"`) serve as the comparison group, and treatment has no effect before period $G = 3$ (`anticipation = 0`), placing the reference period for cohort $G = 3$ at $t = 2$.
- **Parallel trends from the reference period (`PT-post`):** From $t = 2$ onward, mean untreated potential outcome changes for cohort $G = 3$ equal those of the never treated population (`pt_assumption = "post"`).
- **Representative pair observation (`selection = "none"`):** In this simulation, each unit is assigned to one consecutive two period window $\{1, 2\}$, $\{2, 3\}$, $\{3, 4\}$, or $\{4, 5\}$ uniformly at random, independent of cohort membership and potential outcomes. Consequently, within both $G = 3$ and $G = 0$, the units observed in each window $\{s, t\}$ have the same mean outcome change $E[Y_t - Y_s \mid G]$ as the full cohort population. When observation patterns depend on measured covariates, see [Handling missing outcomes](../methods/missing-outcomes.md).
- **Positive support:** Each window has positive population probability in both cohorts. The realized sample must also contain treated and never treated units in every required window; the script checks this and stops if it fails, without redrawing.

## Worked example: rotating consecutive panel in base R



### Simulated truth in the data generating process

In the simulation script below:
- There are $N = 1000$ units in the study population and $T = 5$ calendar dates ($t = 1, 2, 3, 4, 5$).
- Each unit belongs either to treatment cohort $G = 3$ or to the never treated comparison group $G = 0$.
- Untreated outcomes equal an individual normal fixed effect $\alpha_i \sim \mathcal{N}(0, 1)$, a common calendar trend, plus normal noise $\varepsilon_{i,t} \sim \mathcal{N}(0, 0.5^2)$.
- For units in cohort $G = 3$, the treatment effect at posttreatment dates $t \ge 3$ follows the formula

$$\tau(3, t) = 1 + 0.25 \times (t - 3)$$

This formula implies the following true population effects:
- At $t = 3$ (event time $e = t - G = 0$): $\text{ATT}(3, 3) = 1 + 0.25 \times (3 - 3) = 1.00$.
- At $t = 4$ (event time $e = t - G = 1$): $\text{ATT}(3, 4) = 1 + 0.25 \times (4 - 3) = 1.25$.
- At $t = 5$ (event time $e = t - G = 2$): $\text{ATT}(3, 5) = 1 + 0.25 \times (5 - 3) = 1.50$.
- Dynamic average across the three posttreatment event times ($e \in \{0, 1, 2\}$): $(1.00 + 1.25 + 1.50) / 3 = 1.25$.
- At $t = 1$ (preperiod comparison between $t = 1$ and $t = 2$): true pretrend differential equals $0.00$.

### R code

```r
set.seed(2026)

# 1. Create the table of study units (one row per unit in the study population).
N <- 1000L
T_dates <- 5L

units <- data.frame(
  id = seq_len(N),
  g  = ifelse(stats::runif(N) < 0.5, 3, 0)
)

# 2. Allocate each unit independently and uniformly to one consecutive two date
# window: {1, 2}, {2, 3}, {3, 4}, or {4, 5}.
window_start <- sample.int(4L, size = N, replace = TRUE)
stopifnot(all(table(factor(units$g, levels = c(0, 3)),
                    factor(window_start, levels = 1:4)) > 0))
alpha_i <- stats::rnorm(N, mean = 0, sd = 1)

# 3. Expand to full long panel rows (N * T rows) and mask outcomes outside each
# unit's assigned two consecutive dates with NA.
panel_data <- expand.grid(
  time = seq_len(T_dates),
  id   = units$id
)
panel_data <- panel_data[order(panel_data$id, panel_data$time), c("id", "time")]
rownames(panel_data) <- NULL

panel_data$g <- units$g[panel_data$id]

# Generate full potential outcomes with normal individual effect and normal noise.
eps_it <- stats::rnorm(nrow(panel_data), mean = 0, sd = 0.5)
treat_effect <- ifelse(
  panel_data$g == 3 & panel_data$time >= 3,
  1 + 0.25 * (panel_data$time - panel_data$g),
  0
)
y_full <- alpha_i[panel_data$id] + 0.2 * panel_data$time + treat_effect + eps_it

# Retain y only on the two consecutive dates {w, w + 1}; set y = NA elsewhere.
in_window <- (panel_data$time == window_start[panel_data$id]) |
             (panel_data$time == window_start[panel_data$id] + 1L)
panel_data$y <- ifelse(in_window, y_full, NA_real_)

# 4. Fit chained DiD with inverse probability weighting (IPW), PT-post,
# never treated controls, no covariate formula, and analytic standard errors.
fit_rot <- cdid::att_gt(
  yname                  = "y",
  tname                  = "time",
  idname                 = "id",
  gname                  = "g",
  xformla                = NULL,
  data                   = panel_data,
  panel                  = TRUE,
  allow_unbalanced_panel = TRUE,
  control_group          = "nevertreated",
  est_method             = "ipw",
  selection              = "none",
  pt_assumption          = "post",
  bstrap                 = FALSE,
  cband                  = FALSE,
  frame                  = units
)

# 5. Inspect group and time effects and dynamic event time aggregation.
summary(fit_rot)
agg_dyn <- cdid::aggte(fit_rot, type = "dynamic")
summary(agg_dyn)
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

cdid chained difference-in-differences (chained IPW-GMM): cdid_AGGTEobj
  cdid 1.1.0.0 result built on did 2.5.1 (tested: did 2.5.1); did's AGGTEobj methods tabulate the estimates below.
  Chained estimators: Bellego C., Benatia D. and V. Dortet-Bernadet (2025). "The Chained Difference-in-Differences." Journal of Econometrics 248, 105783. doi:10.1016/j.jeconom.2024.105783.
  Ordinary estimators: Callaway B. and P. H. C. Sant'Anna (2021). "Difference-in-Differences with Multiple Time Periods." Journal of Econometrics 225(2), 200-230. doi:10.1016/j.jeconom.2020.12.001.

Call:
cdid::aggte(fit_rot, type = "dynamic")

Overall summary of ATT's based on event-study/dynamic aggregation:  
    ATT    Std. Error     [ 95%  Conf. Int.]  
 1.1488        0.1125     0.9283      1.3692 *


Dynamic Effects:
 Event time Estimate Std. Error [95% Pointwise  Conf. Band]  
         -1   0.0469     0.0823         -0.1143      0.2081  
          0   0.8946     0.0922          0.7139      1.0753 *
          1   1.2016     0.1247          0.9573      1.4460 *
          2   1.3501     0.1573          1.0419      1.6583 *
---
Signif. codes: `*' confidence band does not cover 0

Control Group:  Never Treated,  Anticipation Periods:  0
Estimation Method:  Inverse Probability Weighting

Chained dynamic aggregation: 4 of 4 intended aggregate(s) available; na.rm = FALSE.
  Overall summary: available over its 3 intended bins; pointwise interval.
  Parallel trends: PT-post; base_period = "varying".
  Inference: analytic standard errors, unit-level; pointwise intervals.
```

</details>

### Key options used in this call

| Parameter | Setting used | Purpose and reference link |
| --- | --- | --- |
| `panel`, `allow_unbalanced_panel` | `TRUE`, `TRUE` | Selects chained mode for panel data and retains all units with repeated observations across incomplete windows ([Option reference](../reference/options.md)). |
| `control_group` | `"nevertreated"` | Uses units never treated (`g == 0`) as the comparison group ([Option reference](../reference/options.md)). |
| `est_method` | `"ipw"` | Uses treatment weighting (inverse probability weighting, IPW) for treatment adjustment ([Treatment adjustment](../methods/treatment.md)). |
| `selection` | `"none"` | Uses baseline mean missingness assumptions without fitting an observation probability or change imputation model ([Handling missing outcomes](../methods/missing-outcomes.md)). |
| `pt_assumption` | `"post"` | Imposes parallel trends from each cohort's reference period (`PT-post`) without restricting trends before $t = 2$ ([Parallel trends](../methods/parallel-trends.md)). |
| `bstrap`, `cband` | `FALSE`, `FALSE` | Computes analytical standard errors and pointwise normal confidence intervals for `att_gt()` ([Results and uncertainty](../reference/results.md)). |
| `frame` | `units` | Supplies the complete unit table (`N = 1000` units) representing the study population ([Prepare your data](../help/prepare-data.md)). |

## Interpretation

Even though every unit in `panel_data` has three missing outcomes (`y = NA`) and no unit is observed at both $t = 2$ and $t = 4$ or at both $t = 2$ and $t = 5$, `cdid::att_gt()` chains the consecutive window comparisons $\{2, 3\}$, $\{3, 4\}$, and $\{4, 5\}$ to estimate $\text{ATT}(3, 3)$, $\text{ATT}(3, 4)$, and $\text{ATT}(3, 5)$, while using $\{1, 2\}$ to report the preperiod comparison before treatment. Passing `fit_rot` to `cdid::aggte(fit_rot, type = "dynamic")` aggregates these effects by event time $e = t - G$ without refitting the underlying models, using the full unit table influence functions stored in `fit_rot` to account for the covariance across chained links. You can inspect the recorded method settings and support details in `fit_rot$chained` (see [Results reference](../reference/results.md)).

## Next steps

- To check your own longitudinal dataset, structure the unit table (`frame`), and handle reference dates for covariates or weights that change over time, read [Prepare your panel and unit table](../help/prepare-data.md).
- When observation patterns depend on measured characteristics, see [Handling missing outcomes](../methods/missing-outcomes.md) to adjust via observation weighting (`selection_est_method = "ipw"`), change imputation (`selection_est_method = "imputation"`), or doubly robust adjustment for missing outcomes (`selection_est_method = "dr"`).
- To compare treatment weighting (`est_method = "ipw"`), outcome regression (`est_method = "reg"`), and doubly robust treatment estimation (`est_method = "dr"`), see [Treatment adjustment](../methods/treatment.md).


## Returned effects

| target.id | estimate | std.error | conf.low | conf.high | available |
| --- | --- | --- | --- | --- | --- |
| pre(g=3,t=2) | 0.0469 | 0.0823 | -0.1143 | 0.2081 | TRUE |
| ATT(g=3,t=3) | 0.8946 | 0.0922 | 0.7139 | 1.0753 | TRUE |
| ATT(g=3,t=4) | 1.2016 | 0.1247 | 0.9573 | 1.4460 | TRUE |
| ATT(g=3,t=5) | 1.3501 | 0.1573 | 1.0419 | 1.6583 | TRUE |

---

**Reading path, step 3 of 20** · Previous: [Moving from did](../from-did.md) · Next: [Identify effects when the observed units change composition](identification.md)

<details>
<summary>All steps of the reading path</summary>

| Step | Page |
| --- | --- |
| 1 | [Get started](../get-started.md) |
| 2 | [Moving from did](../from-did.md) |
| 3 | **Incomplete panels** (this page) |
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
| 14 | [Add PT-all restrictions](pt-all-comparison.md) |
| 15 | [GMM weighting and covariance regularization](../methods/gmm.md) |
| 16 | [Robustness](robustness.md) |
| 17 | [Inference](inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
