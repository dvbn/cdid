# Cross-Pair Chaining, Control Parallel Trends (`PT-ctrl`), and Within-Unit GMM Weighting

If the target cohort and its controls are observed over different intervals, they may have no date pair in common. Cross-pair chaining uses changes in each group separately and links them through calendar dates where the trend assumptions apply. Set `cross_pair = TRUE` to allow those comparisons and use `chain_report()` to see how they identify an effect.

The [proof of identification and regression equivalence](../theory/cross_pair_and_control_trends.tex) gives the mathematical conditions. The [simulation examples](../outputs/cross_pair_monte_carlo_results.md) show how the choices matter in different observation patterns.


---

A simple extension of the original chained DiD idea links successive control cohorts. Balla-Elliott and Norwich proposed this approach for cohorts with short pre-treatment windows. Their stacked cohort fixed-effects regression reproduces the calendar trend by linking those windows. The benchmark below shows how to reproduce its point estimates in `cdid` when the cohorts share observation dates within each cohort, and when the within-unit version applies. [Cohort-Chained DiD, footnote 2 and Section 2](https://pdfs.dballaelliott.com/ccdid.pdf) · [Original chained DiD paper](https://doi.org/10.1016/j.jeconom.2024.105783)

## Options for chaining

| Dimension | Default `cdid::att_gt()` | Stacked unit-FE slice benchmark in `cdid::att_gt()` | When to use each setting |
| --- | --- | --- | --- |
| **Moment pairing (`cross_pair`)** | `FALSE`: same-pair DiD moments $(s, t)$ where both target $g$ and clean controls $K$ are observed at $(s, t)$ | `TRUE`: separate target change moments $\Delta Y_g(s, t)$ and clean control change moments $\Delta Y_k(s, t)$ | Use `cross_pair = TRUE` when target and clean controls lack shared date pairs; either `FALSE` or `TRUE` works on contiguous windows |
| **Parallel trends (`pt_assumption`)** | `"post"`: parallel trends from each cohort's reference date $b_g$ onward | `"ctrl"`: adds parallel pre-trends among clean control cohorts $k, k' > g$ from `pt_start` onward | Prefer `"post"` unless pre-reference control pooling is specifically needed for disconnected pre-periods |
| **GMM metric (`gmm_weighting`)** | `"optimal"`: feasible two-step GMM using empirical unit influence functions | `"within"`: analytical within-unit covariance of changes under uncorrelated level shocks | Use `"optimal"` to estimate the moment covariance from the data. Use `"within"` for analytical weighting. |
| **Fixed effects / attrition** | Within-unit differencing $\Delta Y_{i,s,t}$ (cancels unit $\alpha_i$) | Within-unit differencing $\Delta Y_{i,s,t}$ (matches stacked unit-FE slice OLS under nested/monotone attrition, and cohort-FE slice OLS on cohort-balanced windows) | Differencing eliminates the fixed unit level. Representative changes and the other identifying assumptions are still required. |

For the stacked target-slice benchmark, match the target and clean-control sample for each slice, set `selection = "none"` and `xformla = NULL`, and use the same unit weights. Fit with `control_group = "notyettreated"`, `base_period = "universal"`, `cross_pair = TRUE`, `pt_assumption = "ctrl"`, `pt_start = NULL` and `gmm_weighting = "within"`. With nested observation windows within each cohort, the available point estimates equal the unit-fixed-effect regression coefficients. With a common window for all units in each cohort, they also equal the cohort-fixed-effect coefficients. The installed tests check these cases against independent regressions. Inference in `cdid` uses its original-unit influence functions.

---

## 2. Identification: Contiguous Windows vs. Disjoint Cross-Pair Windows

### 2.1 Contiguous short pre-treatment windows need only `PT-post`

In staggered designs with short pre-treatment windows (for example, cohorts `2005..2010` each observed from $g - 4$ onward with no never-treated group), each control cohort's clean untreated dates form a contiguous interval of stored calendar dates.

Whenever $b_g$ and $t$ lie in the same connected component of the cohort-date observation graph under contiguous cohort windows, every consecutive step $(u_{h-1}, u_h)$ between $b_g$ and $t$ is covered by at least one cohort $k > g$ still untreated at $u_h$. That cohort belongs to the clean control pool $K_{\mathrm{orig}}$ for $(u_{h-1}, u_h)$, so **standard `cdid::att_gt(..., control_group = "notyettreated")` with `cross_pair = FALSE` and `pt_assumption = "post"` already identifies the exact same post-treatment effects**—including long-horizon cells where two-period `did::att_gt()` returns `NA`—and chains directly across biennial or irregular wave spacing.

### 2.2 Non-overlapping observation pairs (`cross_pair = TRUE`)

In rotating or irregular panels where a target cohort and its clean controls are never observed on the *same* pair of dates—for example, target cohort $G = 3$ observed at $\{1, 2, 4\}$, control cohort $G = 4$ at $\{2, 3\}$, and control cohort $G = 5$ at $\{3, 4\}$—no control unit has an observed change over $(2, 4)$, and no target unit is observed at date $3$. Same-pair DiD has zero causal links.

Setting `cross_pair = TRUE` constructs separate within-unit change moments for the target cohort ($\gamma_g(t) - \gamma_g(s)$) and for each clean control cohort ($\lambda_{g,k}(t) - \lambda_{g,k}(s)$), anchored at $\gamma_g(b_g) = \lambda_{g,k}(b_g) = 0$. The solver then recovers
$$
\mathrm{ATT}(3, 4) = \mathbb{E}[Y_4 - Y_2 \mid G = 3] - \mathbb{E}[Y_3 - Y_2 \mid G = 4] - \mathbb{E}[Y_4 - Y_3 \mid G = 5]
$$
Use `pt_assumption = "post"`, `"ctrl"` or `"all"` as appropriate for the comparison. Cross-pair fits support IPW, OR and DR treatment adjustment and reference-date covariates and weights. The outcome regression can use separate or pooled slopes.

With `selection = "conditional_mtar"`, choose observation weighting, actual-change imputation or augmentation. The full study frame supports both probability models where probabilities are used. The options `selection_universe = "A"` and `"B"` are for observation weighting with a direct pair logit. They require the [transport assumptions and inputs](missing-outcomes.md#advanced-fitting-the-observation-model-on-part-of-the-unit-table). Independent-period logits, change imputation and augmentation use `selection_universe = "F"`.

The rows of `fit$inffunc` follow the original study units. The influence function carries the fitted nuisance contributions, including contributions from units with fewer than two observed outcomes when they enter an active model or a cohort share. Analytical, clustered and multiplier inference use those contributions for cell effects and summaries.

---

## 3. Parallel-Trends Choices: `"post"`, `"ctrl"`, and `"all"`

For a target cohort $g$ with reference date $b_g$:

- **`pt_assumption = "post"` (default in `att_gt()`):** Assumes conditional parallel trends between cohort $g$ and clean control cohorts from $b_g$ onward. Pre-reference differentials are left free per control pool.
- **`pt_assumption = "ctrl"`:** Keeps `PT-post` and assumes clean control cohorts $k, k' > g$ share identical untreated trends from `pt_start` onward (`pt_start = NULL` defaults to the earliest stored calendar date). Target cohort $g$'s own pre-reference drift remains unrestricted, so pre-treatment leads remain genuine falsification checks of $g$'s pre-trends against the pooled controls.
- **`pt_assumption = "all"` (requires explicit `pt_start`):** Restricts both control cohorts and target cohort $g$ to share common untreated trends from `pt_start` onward, setting in-window pre-reference differentials to zero.

Jointly across all target cohorts, standard not-yet-treated `PT-post` for the earliest treated cohort $g_{\min}$ already implies equal untreated trends among clean cohorts from $b_{g_{\min}}$ onward; thus `PT-ctrl`'s additional empirical restriction is concentrated prior to $b_{g_{\min}}$. If a control cohort experiences an idiosyncratic trend shift before $b_g$, `PT-ctrl` (and stacked slice regression) transmits bias into post-treatment estimates, while PT-post leaves that shift unrestricted.

---

## 4. Inspecting Identification Chains with `chain_report()`

`chain_report()` shows which effects need chaining, a set of comparisons that identifies each effect, and the weights used to combine the available comparisons of a `cdid_MP` fit:

```r
fit <- att_gt(
  yname = "y", tname = "period", idname = "id", gname = "G",
  data = panel_df,
  control_group = "notyettreated", base_period = "universal",
  cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
  bstrap = FALSE
)

# Print the comparisons that identify effects for cohort 2005 (set include_pre = TRUE to include pre-treatment cells)
cr <- chain_report(fit, cohort = 2005)
print(cr)

# Tables of effects, identifying comparisons and GMM weights
as.data.frame(cr, what = "summary")
as.data.frame(cr, what = "witness")
as.data.frame(cr, what = "gmm_weights")
```

For every estimable cell, the `contribution` column (`gmm_weight * link_estimate`) in `as.data.frame(cr, what = "gmm_weights")` sums to the reported `ATT(g, t)` estimate up to floating-point rounding.

---

**Reading path, step 21 of 22** · Previous: [Diagnose errors and unavailable results](../help/troubleshooting.md) · Next: [Short pre-treatment windows](../examples/short-pretreatment-windows.md)

[Guide home](../README.md) · [Package overview](../../README.md)
