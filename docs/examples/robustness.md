# Compare model robustness across treatment and missing outcome adjustments

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

When covariates influence which units enter treatment, how untreated outcomes evolve, and which outcomes are recorded across dates, chained difference in differences (`cdid`, extending Bellégo, Benatia and Dortet-Bernadet 2025 and Callaway and Sant'Anna 2021) separates estimation into two distinct modeling layers: **treatment adjustment** (`est_method`) and **handling missing outcomes** (`selection` and `selection_est_method`).

This article explains the four models across the two layers, summarizes the large sample consistency conditions for each supported combination, and walks through concrete misspecification scenarios using the base R script below. For full mathematical definitions, see [Treatment adjustment](../methods/treatment.md) and [Handling missing outcomes](../methods/missing-outcomes.md).

## Question

Suppose an incomplete longitudinal panel tracks a study population listed in a unit table (`frame`) with a baseline covariate `x` that shifts both treatment assignment and outcome trajectories, alongside an additional variable `z` that shifts whether outcomes are observed. How do the treatment adjustment models and the missing outcome adjustment models interact in `cdid`, which estimators remain consistent when one model is misspecified, and when does pooling slopes or supplying extra selection covariates invalidate a doubly robust branch?

## Data

### Study population and four labeled models

Consider a study population of $N = 2000$ units in the unit table `units` (`frame`) across $T = 6$ calendar dates $t \in \{1, 2, 3, 4, 5, 6\}$, with treatment cohort $G = 3$ (`first_treat == 3`), never treated controls ($G = 0$, `first_treat == 0`), and static unit covariates `x` and `z`. In the panel dataset `obs_panel`, some outcomes are missing, so only a subset of units is observed at both endpoints of any ordered date pair $e = (s, t)$. In the script, `x` and `z` are independent Uniform(-1, 1) variables. The stated treatment and date-observation logits are then bounded away from zero and one in the population; this does not guarantee a successful fit in every finite draw.

To recover average treatment effects for the study population in `units`, `cdid` uses up to four labeled models across two layers:

| Symbol | Model name | Layer and active options | What it models and where it is fitted |
| --- | --- | --- | --- |
| $p$ | Treatment propensity model | Treatment adjustment (`est_method = "ipw"` or `"dr"`, using `xformla`) | For each distinct risk set, a logit model of membership in target cohort $g$ versus `control_group`, weighted by the chosen missing outcome route. |
| $m$ | Counterfactual untreated change regression | Treatment adjustment (`est_method = "reg"` or `"dr"`, using `xformla`) | One pooled weighted least squares regression (`n_regression_fits = 1`) fitted on control units across all deduplicated comparison blocks: $m_b(X) = \alpha(b) + L(e) z(X)' \beta$, with local block intercepts $\alpha(b)$ and one shared slope vector $\beta$. |
| $q$ | Observation pair probability model | Handling missing outcomes (`selection = "conditional_mtar"` with `selection_est_method = "ipw"` or `"dr"`, using `selection_xformla`) | For each needed cohort and date pair cell $c = (k, e)$, a logit model $q_c(W)$ for observing both dates of pair $e$, either directly (`"pair_logit"`) or as the product of two date logits (`"independent_period_logit"`). |
| $j$ | Actual cohort change imputer | Handling missing outcomes (`selection = "conditional_mtar"` with `selection_est_method = "imputation"` or `"dr"`, using `xformla`) | For each needed cohort and date pair cell $c = (k, e)$, a weighted least squares linear model $j_c(X) = \alpha(k, e) + L(e) z(X)' \beta(r)$ fitted to observed changes of cohort $k$, **treatment effects included**, and evaluated across every unit of cohort $k$ in `frame`. |

### How the two outcome regressions differ

Both $m$ and $j$ use `xformla` and scale covariate slopes by $L(e)$, the number of calendar steps between the positions of dates $s$ and $t$ in the stored calendar (never raw subtraction of date labels), but they serve different layers:

1. **Counterfactual untreated change regression ($m$):** Fitted only to control units across all comparison blocks (including preperiod comparisons). It always estimates one global slope vector $\beta$ across all blocks and periods, with one local intercept $\alpha(b)$ per date pair and control pool block $b$ whenever `xformla` includes an intercept (`~ 0 + x` omits local block intercepts; `~ 1` fits block intercepts with zero slopes; `~ 0` is refused). It has no pooling option.
2. **Actual cohort change imputer ($j$):** Fitted within each cohort $k$ (both treated and control cohorts) on units observed at both endpoints of pair $e$ ($S_e = 1$), using fixed unit weights. For a treated cohort after treatment starts, $j_c(X)$ predicts the **actual** outcome change $\Delta Y_e = Y_t - Y_s$, **including treatment effects**. It keeps a local intercept $\alpha(k, e)$ per cell when `xformla` has an intercept and shares its slope vector $\beta(r)$ according to `imputation_pooling`: `"separate"` (default, separate slopes per cohort and pair cell), `"within_cohort"` (one slope vector per cohort across its pairs), or `"across_cohorts"` (one slope vector across all cells).

The observation probability model $q$ keeps separate cell intercepts when `selection_xformla` includes an intercept and pools slopes via `selection_pooling` (`"separate"`, `"within_cohort"`, or `"across_cohorts"`), independently of `imputation_pooling`. Only `selection_xformla` accepts an empty design (`~ 0`), fixing each pair probability to $0.5$ (or $0.25$ for a two date product).

### Four routes for handling missing outcomes

- **Without an observation correction (`selection = "none"`):** Fits neither $q$ nor $j$. Observed pair changes within each cohort are scaled by the cohort's empirical pair response rate and passed to treatment adjustment ($p$, $m$, or both).
- **Observation weighting (`selection = "conditional_mtar", selection_est_method = "ipw"`):** Fits $q_c(W)$ only and weights observed pair changes by inverse fitted pair probabilities combined with cohort repeat rates before running $p$, $m$, or both.
- **Change imputation (`selection = "conditional_mtar", selection_est_method = "imputation"`):** Fits $j_c(X)$ on observed pairs and averages the predictions $j_c(X)$ over **every** unit of cohort $k$ in the unit table (`selection_universe = "F"`, `selection_rho = NULL`), including units whose change was observed. No observation model $q$ is fitted (`selection_model`, `selection_pooling`, and `selection_xformla` are inactive and never evaluated). Treatment adjustment ($p$, $m$, or both) runs on these imputed changes.
- **Augmented correction (`selection = "conditional_mtar", selection_est_method = "dr"`):** Requires `est_method = "dr"` (requests with `"ipw"`, `"reg"`, or a custom function are refused before the data are read) on `selection_universe = "F"` (`selection_rho = NULL`). For each needed cell $c = (k, e)$, every unit of cohort $k$ in `frame` receives the reconstructed change

$$Z_c = j_c(X) + S_e \frac{\Delta Y_e - j_c(X)}{\max(q_c(W),10^{-6})}$$

where $S_e = 1$ for units observed at both dates of pair $e$ and $0$ otherwise, and $W$ is the union of the covariates in `xformla` and `selection_xformla`. The residual adjustment is not renormalized. Version 1.1.0.0 floors its inverse observation denominator at $10^{-6}$ ([Probability safeguards](../methods/missing-outcomes.md#probability-safeguards)). Treatment doubly robust estimation ($p$ and $m$) runs on $Z_c$, forming a single estimator that generally differs in finite samples from both observation weighting and change imputation and is not their average.

## Assumptions

### Shared identification and regularity conditions

Every consistency result below is a qualified large sample property requiring the following conditions and, since version 1.1.0.0, inactive probability safeguards in population or additional conditions preserving the original ATT. When clipping binds, the usual robustness branches need not hold; influence-function inference also requires negligible threshold-crossing remainder ([Probability safeguards](../methods/missing-outcomes.md#probability-safeguards)).

- **Scientifically sufficient conditioning:** Writing formulas in `xformla` or `selection_xformla` does not prove sufficiency; the covariates must justify parallel trends and missing outcome assumptions in the population.
- **Parallel trends and no anticipation:** Conditional on `xformla`, untreated potential outcome changes satisfy `PT-post` (or `PT-all` when selected) with no anticipation beyond `anticipation` ([Parallel trends](../methods/parallel-trends.md)).
- **Valid missing outcome condition on every used pair:** Under `selection = "none"`, observed pair changes must represent each cohort's study population mean for every used pair. Under `selection = "conditional_mtar"`, conditional mean trends missing at random must hold for **every** cohort $k$ (including actual changes of treated cohorts) and **every** date pair $(s, t)$ used by the catalog, including preperiod comparisons:

$$E[Y_t - Y_s \mid W, G = k, S_{s,t} = 1] = E[Y_t - Y_s \mid W, G = k]$$

- **Overlap, support, and valid GMM moments:** Treatment propensities stay strictly below one; observation weighting and augmented correction require $q_c(W) > 0$; change imputation requires observed training mass and prediction support. Active logits need finite regular interior roots and their required design rank; active WLS models need identified, well-conditioned fitting directions and valid required predictions, and all catalog links combined by GMM must satisfy valid moment conditions for the common target parameter (incompatible link pseudo targets make the limit depend on the GMM weights).

### Correctness conditions across both layers

| Missing outcome route | `est_method = "ipw"` | `est_method = "reg"` | `est_method = "dr"` |
| --- | --- | --- | --- |
| `selection = "none"` (no $q$ or $j$) | $p$ correct | $m$ correct across all blocks | $p$ correct **or** $m$ correct |
| `selection = "conditional_mtar"`, `selection_est_method = "ipw"` (fits $q$) | $q$ correct in every cell **and** $p$ correct | $q$ correct in every cell **and** $m$ correct | $q$ correct in every cell **and** ($p$ correct **or** $m$ correct) |
| `selection = "conditional_mtar"`, `selection_est_method = "imputation"` (fits $j$) | $j$ correct in every cell/group **and** $p$ correct | $j$ correct in every cell/group **and** $m$ correct | $j$ correct in every cell/group **and** ($p$ correct **or** $m$ correct) |
| `selection = "conditional_mtar"`, `selection_est_method = "dr"` (fits $q$ and $j$) | Refused | Refused | For every needed cell $c$, ($q_c$ correct **or** $j_c$ correct), **and** ($p$ correct **or** $m$ correct) |

Four essential qualifications govern this table:

1. **Treatment DR supplies no general correction for a wrong imputer:** Under pure imputation, the sufficient condition still requires correct $j$ within each needed cell/pooling group. Choosing treatment DR provides the $p$-or-$m$ alternative after reconstruction; it does not generally repair incorrect reconstructed changes. A wrong $j$ can therefore bias the result even when $p$ and $m$ are correct. These conditions are sufficient, not necessary: special projection, aggregation or error cancellations can preserve an effect despite a misspecified nuisance.
2. **Coupling through shared fits:** Because $m$ shares one slope vector $\beta$ across all blocks, and because `imputation_pooling`, `selection_pooling`, and `"independent_period_logit"` share parameters across cells or dates, misspecification in one cell or block can shift a shared slope everywhere it is used.
3. **Common conditioning set $W$ and covariates only in `selection_xformla`:** Under the augmented correction, `xformla` builds $j$, $m$, and $p$ from covariates $X$, whereas `selection_xformla` builds $q$ from $W = (X, Z)$. Variables $Z$ only in `selection_xformla` never enter $j$, $m$, or $p$. If $Z$ shifts the conditional mean of actual changes given W but is absent from `xformla`, $j_c(X)$ cannot generally represent that W-conditional mean. The stated imputer-correct sufficient branch then fails; correct $q_c(W)$ together with valid treatment adjustment supplies the alternative branch. This does not prove that every effect is biased when both nuisance specifications fail. Adding selection covariates is not automatically harmless.
4. **All active fit failures are retained without fallback:** Active logits must pass their rank/interior-root checks, and active regressions must pass their own rank/conditioning and prediction-support checks, even when the alternative doubly robust model is correct. A completely observed or completely unobserved cell with a free intercept in $q$ has no finite interior logit root and stops with an error (and pooling slopes via `selection_pooling` retains cell intercepts, so it does not prevent separation). Pure imputation runs on complete cells because $q$ is inactive, but the augmented correction refuses whenever $q$, $j$, $p$, or $m$ fails. Fixed probability safeguards apply to finite fitted weights; they do not rescue failed fits, and nothing falls back automatically.

## Code



```r
set.seed(20261027)

n_units <- 2000L
n_periods <- 6L

# Bounded continuous covariates keep the stated population logits away from 0/1.
x <- stats::runif(n_units, min = -1, max = 1)
z <- stats::runif(n_units, min = -1, max = 1)
first_treat <- ifelse(stats::runif(n_units) < stats::plogis(0.6 * x), 3L, 0L)
unit_effect <- stats::rnorm(n_units)

units <- data.frame(id = seq_len(n_units), first_treat = first_treat, x = x, z = z)

full_panel <- data.frame(
  id          = rep(units$id, each = n_periods),
  period      = rep(seq_len(n_periods), times = n_units),
  first_treat = rep(units$first_treat, each = n_periods),
  x           = rep(units$x, each = n_periods),
  z           = rep(units$z, each = n_periods)
)

treated_now <- full_panel$first_treat > 0L & full_panel$period >= full_panel$first_treat
effect <- ifelse(
  treated_now,
  1 + 0.25 * (full_panel$period - full_panel$first_treat) + 0.5 * full_panel$x,
  0
)

full_panel$y <- rep(unit_effect, each = n_periods) +
  (0.4 + 0.3 * full_panel$x) * full_panel$period +
  effect +
  stats::rnorm(nrow(full_panel), mean = 0, sd = 0.5)

# This is scenario 1 only; later scenario changes are discussed, not executed.
observed <- stats::runif(nrow(full_panel)) <
  stats::plogis(1.2 + 0.7 * full_panel$x + 0.6 * full_panel$z)
obs_panel <- full_panel[observed, ]

fit_none_dr <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = obs_panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  selection = "none", bstrap = FALSE, cband = FALSE
)

fit_obs_or <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, selection_xformla = ~ x + z, data = obs_panel, frame = units,
  control_group = "nevertreated", est_method = "reg",
  selection = "conditional_mtar", selection_est_method = "ipw",
  selection_model = "independent_period_logit", selection_pooling = "separate",
  bstrap = FALSE, cband = FALSE
)

fit_obs_dr <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, selection_xformla = ~ x + z, data = obs_panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  selection = "conditional_mtar", selection_est_method = "ipw",
  selection_model = "independent_period_logit", selection_pooling = "separate",
  bstrap = FALSE, cband = FALSE
)

fit_imp_dr <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, data = obs_panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  selection = "conditional_mtar", selection_est_method = "imputation",
  imputation_pooling = "separate", bstrap = FALSE, cband = FALSE
)

fit_aug_dr <- cdid::att_gt(
  yname = "y", tname = "period", idname = "id", gname = "first_treat",
  xformla = ~ x, selection_xformla = ~ x + z, data = obs_panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  selection = "conditional_mtar", selection_est_method = "dr",
  selection_model = "independent_period_logit", selection_pooling = "separate",
  imputation_pooling = "separate", bstrap = FALSE, cband = FALSE
)

fits <- list(none_dr = fit_none_dr, obs_or = fit_obs_or, obs_dr = fit_obs_dr,
             imp_dr = fit_imp_dr, aug_dr = fit_aug_dr)
simple_summaries <- lapply(fits, cdid::aggte, type = "simple")
```

Table derived from the five returned `simple_summaries` objects (the code above assigns them without printing):

```text
method estimate  std.error
  None + treatment DR 1.436272 0.03399377
 Observation IPW + OR 1.381341 0.03226414
 Observation IPW + DR 1.391160 0.03292001
      Imputation + DR 1.377307 0.03253647
  Observation DR + DR 1.388394 0.03291157
```

## Interpretation

### Concrete misspecification scenarios

The code above specifies only scenario 1. Scenarios 2 to 5 are hypothetical modifications, not additional fitted objects or executed checks. Hold the formulas fixed as indicated; in scenario 4 change only the active imputation pooling options. “Consistent” below names a sufficient large sample branch, conditional on all stated assumptions and regular fitted limits. “No guarantee” means those sufficient conditions no longer establish consistency; it does not prove bias for a particular effect or exclude cancellation.

| Scenario | Correct vs misspecified models | `fit_none_dr` | `fit_obs_or` | `fit_obs_dr` | `fit_imp_dr` | `fit_aug_dr` |
| --- | --- | --- | --- | --- | --- | --- |
| **1. Baseline script:** `z` shifts observation only; linear changes in `x`; logit $p(x)$ and date logit $q(x, z)$ | $p$, $m$, $q$, $j$ all correct | Generally biased | Consistent | Consistent | Consistent | Consistent |
| **2. Nonlinear untreated trend:** replace `0.3 * x` by `0.3 * x^2` in outcomes | $p$, $q$ correct; linear $m$, $j$ misspecified | No guarantee | No guarantee | Consistent | No guarantee | Consistent |
| **3. Wrong observation model $q$:** independent date observation uses `z^2`, but fitted date logits keep `~ x + z` | $j$ correct; imputation/augmented treatment models valid on W; fitted $q$ generally misspecified | No guarantee | No guarantee | No guarantee | Consistent | Consistent if the active $q$ fit has a regular interior limit |
| **4. Mispooled imputer $j$:** set `imputation_pooling = "within_cohort"` or `"across_cohorts"` | $q$ correct and the original counterfactual model representable; imposed $j$ pooling misspecified | No guarantee | Consistent | Consistent | No guarantee | Consistent |
| **5. `z` also shifts outcomes:** add `0.4 * z * period` while keeping `xformla = ~ x` | $q$ and $p$ correct on W; $j$ and $m$ omit the outcome-shifting `z` | No guarantee | No guarantee from the stated W-based branch | Consistent through $q$ and $p$ | No guarantee | Consistent through $q$ and $p$ |

- **Scenario 1 (Baseline script):** Because `x` raises both the observation probability and the treatment effect, `fit_none_dr` (`selection = "none"`) is generally biased. Because `z` is independent of `(x, first_treat)` and affects only observation, $E[\Delta Y_e \mid x, z, G = k] = E[\Delta Y_e \mid x, G = k]$ is linear in $L(e) x$ within every cell, so $p$, $m$, $q$, and $j$ are all correctly specified.
- **Scenario 2 (Nonlinear untreated trend in `x`):** Replacing $0.3 x$ with $0.3 x^2$ misspecifies both linear regressions $m$ and $j$ (a nonlinear term is not harmless without a special cancellation). For pure imputation, treatment DR alone no longer supplies a consistency guarantee once $j$ is wrong; a particular effect could still benefit from a special cancellation. Both `fit_obs_dr` and `fit_aug_dr` remain consistent because $q_c(x, z)$ and $p(x)$ are correct.
- **Scenario 3 (Misspecified observation model $q$):** Replace the observation logit by `1.2 + 0.7 * x + 0.6 * z^2`, keeping independent date observation and fitted `selection_xformla = ~ x + z`. This preserves conditional mean-MTAR: observation uses only the covariates and independent draws, while `z` still affects no changes given `x`. The fitted $q$ is generally misspecified. Observation weighting then lacks its stated sufficient guarantee; wrong $q$ need not bias every effect. Correct separate-cell $j$ and valid treatment adjustment preserve the imputation and augmented branches, provided every active fit remains regular. Permanent dropout is not implemented here and would require its own stated missingness, support and date-dependence conditions.
The treatment propensity in a reconstructed route is a synthetic weighted-cohort logit, not necessarily the original assignment probability. With correct reconstruction in this law, cohort-rate rescaling changes its intercept while the `~ x` design remains adequate. When an observation model or imputer is wrong, do not infer correctness of every downstream fitted nuisance merely because the original assignment or untreated-outcome equation was simple.

- **Scenario 4 (Mispooled imputer $j$):** Cohort 3's actual change slope on `x` is $0.3 L(e)$ before treatment but shifts by $0.5$ across pairs spanning period 3. A single duration-scaled slope within the treated cohort or across cohorts cannot generally represent both crossing and noncrossing pairs. This misspecifies $j$ and removes the pure-imputation guarantee; it does not misspecify the separate shared counterfactual regression $m$. Correct $q$ preserves the augmented branch under its remaining assumptions.
- **Scenario 5 (Extra selection covariate `z` shifts outcome changes):** When `z` directly shifts outcome changes but is omitted from `xformla = ~ x`, $j_c(x)$ is misspecified on $W = (x, z)$. Pure imputation therefore has no guarantee from its stated $j$ branch. Here `z` is independent of treatment conditional on `x`, so $p$ remains correctly specified on W even though its formula uses only `x`; together with correct $q$, this supplies the observation-DR and augmented branches. A different treatment law would require checking that conclusion again. No claim of inevitable OR or imputation bias follows from omitted-variable labels alone.

### What a repeated sample robustness check requires

Evaluating finite sample performance across these specifications requires a separate repeated sample Monte Carlo study that computes the exact population truth from the data generating law, runs many independent replicates per design, and reports bias, root mean squared error (RMSE), empirical interval coverage with Monte Carlo standard errors, and the exact count of failed or refused fits. No repeated sample campaign has been executed for this article.

## Next steps

- Read [Treatment adjustment](../methods/treatment.md) and [Handling missing outcomes](../methods/missing-outcomes.md) for full method details.
- Continue to [Clustered inference, confidence bands, and pretrend tests](inference.md), [Options reference](../reference/options.md), and [Results reference](../reference/results.md).

---

**Reading path, step 16 of 20** · Previous: [GMM weighting and covariance regularization](../methods/gmm.md) · Next: [Inference](inference.md)

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
| 14 | [Add PT-all restrictions](pt-all-comparison.md) |
| 15 | [GMM weighting and covariance regularization](../methods/gmm.md) |
| 16 | **Robustness** (this page) |
| 17 | [Inference](inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
