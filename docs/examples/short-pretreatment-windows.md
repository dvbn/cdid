# Short Pre-Treatment Windows and Cross-Pair Chaining

Longitudinal panels with staggered adoption often observe each unit only over a limited window around treatment—for instance, administrative registers tracking workers over fixed windows around an event, rotating household surveys, or longitudinal surveys that switch from annual to biennial waves. When no never-treated group is available, standard two-period not-yet-treated estimators (`did::att_gt()`) run out of single-cohort comparison groups at long horizons: by the time early-treated cohort $g$ reaches $g + 4$, later-treated cohorts observed back at $g - 1$ have already entered treatment, while cohorts still untreated at $g + 4$ were not yet observed at $g - 1$.

This guide shows how `cdid::att_gt()` recovers long-horizon effects in short-window designs across both same-pair and cross-pair chained modes, how to inspect the underlying control chains with `chain_report()`, and how to choose among `pt_assumption = "post"` vs `"ctrl"` and `gmm_weighting = "optimal"` vs `"within"`.

---

## 1. Simulating a Staggered Short-Window Panel (No Never-Treated Units)

Consider years `2001` through `2010` and six adoption cohorts `2005, 2006, ..., 2010`, with no never-treated units. Each cohort $g$ enters the panel four years prior to treatment ($\max(2001, g - 4)$ through `2010`):

```r
library(cdid)
set.seed(20261111)

cohorts <- 2005:2010
rows <- list()
uid <- 0L
for (g in cohorts) {
  win_t <- seq.int(max(2001L, g - 4L), 2010L)
  for (i in seq_len(100L)) {
    uid <- uid + 1L
    alpha_i <- rnorm(1, mean = 0.5 * (g - 2005), sd = 1.2)
    for (tv in win_t) {
      tau <- if (tv >= g) 1.0 + 1.0 * (tv - g) else 0.0 # ATT(2005, 2009) = 5.0
      y_it <- alpha_i + 0.4 * (tv - 2001) + tau + rnorm(1)
      rows[[length(rows) + 1L]] <- data.frame(id = uid, period = tv, G = g, y = y_it)
    }
  }
}
df_short <- do.call(rbind, rows)
```

Why does two-period `did::att_gt(..., control_group = "notyettreated")` return `NA` for `ATT(2005, 2008)` and `ATT(2005, 2009)`?
- Cohort `2005`'s reference period is `2004`.
- Only cohorts `2009` and `2010` remain untreated in `2008`, and only cohort `2010` remains untreated in `2009`.
- Cohort `2009` enters the panel in `2005` ($2009 - 4$), and cohort `2010` enters in `2006` ($2010 - 4$)—neither is observed at `2004`.

---

## 2. Analytical Within-Unit Weighting vs. Optimal Chained GMM

In this design, all units in a cohort share the same dates and neither fit adjusts for covariates or observation probabilities. The within-weighted fit matches the available point estimates from stacked target-slice unit- or cohort-fixed-effect OLS on the same clean-control samples. The second fit uses PT-post and feasible two-step GMM weighting.

```r
# 1. Analytical within-unit covariance weighting under PT-ctrl and universal base
fit_within <- att_gt(
  yname = "y", tname = "period", idname = "id", gname = "G",
  control_group = "notyettreated",
  base_period = "universal",
  cross_pair = TRUE,
  pt_assumption = "ctrl",
  gmm_weighting = "within",
  data = df_short
)
agg_within <- aggte(fit_within, type = "dynamic", na.rm = TRUE)
summary(agg_within)

# 2. Optimal two-step GMM under weaker PT-post (no pre-reference control pooling)
fit_post_opt <- att_gt(
  yname = "y", tname = "period", idname = "id", gname = "G",
  control_group = "notyettreated",
  pt_assumption = "post",
  gmm_weighting = "optimal",
  cross_pair = TRUE,
  data = df_short
)
```

Both `fit_within` and `fit_post_opt` (as well as default same-pair `att_gt(..., control_group = "notyettreated")` on contiguous windows) recover all **15 post-treatment effects**, including `ATT(2005, 2008)` ($e = 3$) and `ATT(2005, 2009)` ($e = 4$). Only `t = 2010` is `NA` because cohort `2010` is already treated in `2010`, leaving zero untreated units in the population.

---

## 3. Auditing the Control Chain with `chain_report()`

To see *how* `ATT(2005, 2009)` bridges `2004` to `2009` without any control cohort observed at both dates, call `chain_report()`:

```r
cr <- chain_report(fit_within, cohort = 2005, time = c(2005, 2008, 2009))
print(cr)
```

`chain_report()` flags `requires_chaining = FALSE` for `ATT(2005, 2005)` (directly supported by clean cohorts observed across `2004` and `2005`) and `requires_chaining = TRUE` for `ATT(2005, 2008)` and `ATT(2005, 2009)`, printing a short path of eligible control comparisons and showing the GMM weights via `as.data.frame(cr, what = "gmm_weights")`.

---

## 4. Practical Recommendations

1. **Keep biennial or irregular waves in the data:** `cdid` chains over arbitrary stored calendar dates—so transitions from annual to biennial survey waves do not require dropping post-transition periods.
2. **Always retain unit identifiers (`idname`):** When panel attrition correlates with persistent outcome levels $\alpha_i$, cohort-fixed-effect regressions confound changing cohort composition with treatment effects, while within-unit differencing cancels those persistent levels. The pair-representativeness and parallel-trends assumptions still apply.
3. **Prefer `gmm_weighting = "optimal"` and `pt_assumption = "post"` in empirical work:** Use `gmm_weighting = "within", pt_assumption = "ctrl"` when benchmarking directly against stacked target-slice regressions, and use `gmm_weighting = "optimal"` and `pt_assumption = "post"` to gain two-step GMM efficiency under serial correlation/heteroskedasticity while avoiding bias from pre-reference control trend shifts.
4. **Adjust for covariates and non-random missingness with full `cdid` generality:** Baseline and reference-date covariate differences (`xformla`, `fix_weights`) across early and late adopters are handled via built-in `"ipw"`, `"reg"`, or `"dr"` (`outcome_pooling = "separate"` or `"pooled"`), and covariate-dependent attrition (`selection = "conditional_mtar"`) is supported across `"ipw"`, `"imputation"`, and `"dr"` (`selection_est_method`) with all observation-logit and imputation pooling modes—propagating exact full-roster ($N = |F|$) M-estimation nuisance derivatives into unit- and cluster-robust standard errors, simultaneous confidence bands, and `aggte()` aggregations.

---

**Reading path, step 22 of 22** · Previous: [Cross-pair chaining and control parallel trends](../methods/cross-pair-chaining.md) · Last step; back to the [guide home](../README.md)

[Guide home](../README.md) · [Package overview](../../README.md)
