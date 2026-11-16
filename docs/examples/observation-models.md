# Choose an observation model for missing outcomes

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage. OR/DR calls shown here explicitly retain pooled slopes; the current default uses separate comparison regressions.

When outcomes are missing at some dates, chained DiD uses the changes of units observed at both dates of each pair. You decide how those observed changes stand in for the study population: with no observation equation, with a logit for each date, or with a logit for each pair of dates, and with more or fewer shared slopes. This example makes that choice on one small simulated panel in which the right answer is known by construction. [Handling missing outcomes](../methods/missing-outcomes.md) defines the assumptions and models.



## The illustration

The panel has 2000 units and six dates, `time` from 1 to 6. Each unit has a binary covariate `x`, equal to 1 with probability 0.5 and fixed over time, and a cohort `g` drawn independently of `x`, the unit component and outcome noise: 0 (never treated) or 3, with equal probability. With a unit component `alpha` from a standard normal distribution and independent normal noise with standard deviation 0.5, the untreated outcome is

```
Y0(t) = alpha + 0.1 * t + 0.3 * x * t + noise
```

From date 3 onward, units of cohort 3 add the effect `1 + 0.25 * (t - 3) + 0.5 * x`. In both cohorts, each outcome is observed with probability

```
q_t(x) = plogis(-0.3 + 0.1 * (t - 1) + 0.4 * x)
```

Given `x`, observation is independent across dates and independent of potential outcomes, noise and cohort. This does not assert unconditional independence from outcomes, which also depend on `x`. All units have the same sampling weight, so no `weightsname` is used.

## Generate the panel

```r
set.seed(20261031)

n_units <- 2000
n_periods <- 6

x <- rbinom(n_units, size = 1, prob = 0.5)
g <- sample(c(0, 3), size = n_units, replace = TRUE)
alpha <- rnorm(n_units)

sim <- data.frame(
  id = rep(seq_len(n_units), each = n_periods),
  time = rep(seq_len(n_periods), times = n_units),
  g = rep(g, each = n_periods),
  x = rep(x, each = n_periods)
)

# Untreated outcome, then the actual outcome with the effect from date 3 on.
y_untreated <- rep(alpha, each = n_periods) +
  0.1 * sim$time + 0.3 * sim$x * sim$time +
  rnorm(nrow(sim), mean = 0, sd = 0.5)
treated_now <- sim$g == 3 & sim$time >= 3
sim$y_complete <- y_untreated +
  treated_now * (1 + 0.25 * (sim$time - 3) + 0.5 * sim$x)

# Given x, the observation mask is independent across dates
# and independent of potential outcomes, noise and cohort.
sim$observed <- runif(nrow(sim)) <
  plogis(-0.3 + 0.1 * (sim$time - 1) + 0.4 * sim$x)

# Analysis data: every unit and date row, NA where the outcome is missing.
panel <- sim[c("id", "time", "g", "x")]
panel$y <- ifelse(sim$observed, sim$y_complete, NA)

# Unit table: one row per study unit, including units observed once or never.
units <- data.frame(id = seq_len(n_units), g = g, x = x)

# Number of units by number of observed outcomes, 0 to 6.
table(tapply(!is.na(panel$y), panel$id, sum))
```

Output:

```text
0   1   2   3   4   5   6 
 27 126 411 585 514 283  54
```

`sim` stores the latent complete outcome `y_complete` and the original mask `observed` for checks; the fits read only `panel` and `units`. `panel` keeps every unit and date row and codes a missing outcome as `NA`, never as a filled value ([Prepare your panel and unit table](../help/prepare-data.md)), and `units` lists every unit, including units observed once or never.

## The target

The effects are those of cohort 3 in the study population that `units` lists:

```
ATT(3, t) = 1 + 0.25 * (t - 3) + 0.5 * E[x | g = 3]      for t = 3, 4, 5, 6
```

In the population, `x` has mean 0.5 in both cohorts, so ATT(3, t) is `1.25 + 0.25 * e` at event time `e = t - 3`. For the 2000 generated units, `E[x | g = 3]` becomes `mean(units$x[units$g == 3])`, so their ATT(3, 3) is generally not exactly 1.25. With one treated cohort, the dynamic aggregate at event time `e` has the value of ATT(3, 3 + e).

## Three comparable fits

The fits share the data, the unit table, `xformla = ~ x`, never treated controls, treatment DR, PT-post and `base_period = "varying"`. They keep `compute_inffunc = TRUE`, the default, so their analytical standard errors come from complete influence functions. Only the handling of missing outcomes changes.

```r
fit_none <- cdid::att_gt(
  outcome_pooling = "pooled",
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x, data = panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  pt_assumption = "post", base_period = "varying",
  selection = "none",
  bstrap = FALSE, cband = FALSE
)

fit_date <- cdid::att_gt(
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x, data = panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  pt_assumption = "post", base_period = "varying",
  selection = "conditional_mtar",
  selection_est_method = "ipw",
  selection_model = "independent_period_logit",
  selection_pooling = "separate",
  selection_universe = "F",
  bstrap = FALSE, cband = FALSE
)

fit_pair <- cdid::att_gt(
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x, data = panel, frame = units,
  control_group = "nevertreated", est_method = "dr",
  pt_assumption = "post", base_period = "varying",
  selection = "conditional_mtar",
  selection_est_method = "ipw",
  selection_model = "pair_logit",
  selection_pooling = "separate",
  selection_universe = "F",
  bstrap = FALSE, cband = FALSE
)
```

`fit_none` fits no observation model. `fit_date` and `fit_pair` use the canonical names `selection = "conditional_mtar"` and `selection_est_method = "ipw"`, and both state `selection_pooling = "separate"`, the default for either model.

## What each choice assumes

| Choice | What must hold | Here |
| --- | --- | --- |
| `selection = "none"` | within each cohort, units observed at both dates of a pair have the cohort's mean change | fails |
| `selection = "conditional_mtar"` | the same equality given `x` and cohort, for every cohort and every pair the estimator uses | holds |
| date logits, in addition | observation at the two dates is independent given cohort and `x`, and each date logit is correct | holds |
| pair logits, in addition | the logit of observing both dates is correct in every needed cell | holds with separate cells |

**Representative means.** The baseline condition is weaker than independence of observation and changes within a cohort, and it needs no model. Independence given `x` does not imply it. Here units with `x = 1` are observed more often, so they make up a larger share of the units observed at both dates than of their cohort, and changes depend on `x`: through the trend `0.3 * x * t` and, for cohort 3, through the effect `0.5 * x`. The observed mean changes are therefore not representative, in particular for the treated cohort, and `fit_none` is generally biased for the population effects. Treatment adjustment on `x` does not correct this ([Treatment adjustment](../methods/treatment.md#doubly-robust-treatment-estimation-dr)).

**Conditional mean trends missing at random.** The corrected fits assume the mean equality given `x` and cohort, which covers the actual changes of cohort 3, effects included. It holds here because the design satisfies a stronger condition: given `x`, the mask is independent of the cohort and of every potential outcome. Such conditional independence is sufficient, not required. No fitted observation model can establish either condition; in an application it rests on knowledge of how outcomes went missing.

**Independence across dates.** A correct logit for each date implies nothing about how observation at one date relates to observation at another. The date model adds independence of the two dates, given cohort and `x`, as a separate assumption and divides by `q_s(x) * q_t(x)`, which in this design is exactly the probability of observing both dates. The pair model estimates that probability directly, so the two fits estimate their denominators differently. Their population pair probabilities agree under this design, but the fitted probabilities need not agree in a sample. Persistent dropout usually violates the date model, because a unit that leaves does not return.

**Flexibility.** Pair logits allow any dependence between the two dates but need more cells and slopes: each cohort has up to fifteen pairs of dates and only six dates. Date logits need fewer parameters and the extra independence assumption. Neither is automatically better; choose from how observation arises in your data, not from which fit looks more precise.

## How much the logits share

`selection_pooling` restricts only the slopes. Each cell, a cohort and date for date logits or a cohort and pair for pair logits, keeps its own intercept whenever `selection_xformla` has an intercept, as the inherited `~ x` does.

| `selection_pooling` | Date logits | Pair logits |
| --- | --- | --- |
| `"separate"` (default) | a slope for each cohort and date: correct here | a slope for each cohort and pair: correct here, saturated on binary `x` |
| `"within_cohort"` | one slope per cohort across its dates: correct here | one slope per cohort across its pairs: generally misspecified |
| `"across_cohorts"` | one slope for all cells: correct here | one slope for all cells: generally misspecified |

These entries are algebraic properties of this design, not fit results or guarantees. With binary `x`, a logit with an intercept and a slope reproduces any two probabilities, so separate pair cells are saturated; with a covariate taking more values, the logit of a product of two date probabilities is generally not linear in it. Even here the implied pair slope changes from pair to pair, so shared pair slopes are generally misspecified under this date law. The pooled date logits are correct only because this particular law gives every date and both cohorts the slope 0.4. No pooling choice is guaranteed the smaller standard error. Pooling is a restriction you choose: the package never pools automatically and never changes the model after a failure.

## Formula, fitting set and unit table

- `selection_xformla` can be set independently of `xformla`; its default `NULL` inherits `xformla`, here `~ x`. Covariates added only to the observation model need their own scientific justification: different predictor lists require restrictions sufficient for each model, and supplying a formula does not establish them ([Two covariate formulas](../methods/missing-outcomes.md#two-covariate-formulas)).
- `selection_universe = "F"`, the default, fits the observation model on the whole unit table, and date logits require it. The subsets `"A"` and `"B"` exist only for pair logits under observation weighting, with the transport conditions in [Advanced: fitting the observation model on part of the unit table](../methods/missing-outcomes.md#advanced-fitting-the-observation-model-on-part-of-the-unit-table).
- Every unit of `units` stays in the analysis. A unit observed once or never supplies no change, but it enters the date or pair logits fitted on the whole unit table, counts in the response rates and keeps its row in the influence function.

## Read the results

```r
fits <- list(none = fit_none, date = fit_date, pair = fit_pair)
simple <- lapply(fits, cdid::aggte, type = "simple")
data.frame(
  fit = names(simple),
  estimate = vapply(simple, function(a) a$overall.att, numeric(1)),
  std.error = vapply(simple, function(a) a$overall.se, numeric(1))
)

# Known values for these generated units, computed from the design.
# ATT(3, t) for t = 3 to 6, also the dynamic values for e = 0 to 3:
1 + 0.25 * (3:6 - 3) + 0.5 * mean(units$x[units$g == 3])
# The simple aggregate, the average of those four effects:
1.375 + 0.5 * mean(units$x[units$g == 3])

```

<details>
<summary>Console output</summary>

```text
fit estimate  std.error
none none 1.675595 0.03914203
date date 1.623850 0.03817623
pair pair 1.628628 0.03800096
[1] 1.251773 1.501773 1.751773 2.001773
[1] 1.626773
```

</details>

The simple aggregate averages ATT(3, 3) to ATT(3, 6); its population value is 1.625. The lines computed from the design are known targets for these units, not estimates, and are kept apart from the fitted aggregates. Under the required identification and large sample regularity conditions, both corrected fits target these effects. Their finite sample estimates and standard errors can differ; neither is guaranteed the smaller standard error. `fit_none` is generally biased.

`fit$chained$selection` records the adjustment, the observation model, the pooling choice, the fitting set, the assumptions and the counts of fitted models. For observation weighting its `method` entry keeps the legacy record name `"corollary1"`, with `canonical_assumption = "conditional_mtar"`; the requests above use the canonical names. `fit$chained$nuisance_diagnostics` holds the fitted cells, the pooling groups with their ranks and slopes, the counterfactual regression, the treatment logits and the support of each cell; see [Where the result records the choice](../methods/missing-outcomes.md#where-the-result-records-the-choice).

## When a fit stops

Each active logit needs a finite interior fit of full rank. A date or pair cell in which every unit is observed, or none is, cannot fit a logit with its own free intercept, so the fit stops with an error; pooled slopes keep those intercepts and do not avoid it. Separation, boundary solutions and deficient rank are refused in the same way. Version 1.1.0.0 floors finite inverse observation denominators at `1e-6`; it does not retry failed fits or replace another model. When clipping binds, original-target robustness and inference need further conditions ([Probability safeguards](../methods/missing-outcomes.md#probability-safeguards)). The date probabilities of this design are bounded away from zero and one, but that does not guarantee a successful fit in every draw.

## Related pages

[Impute actual changes for missing outcomes](change-imputation.md) reuses this panel for change imputation and the augmented correction. [Compare model robustness across treatment and missing outcome adjustments](robustness.md) works through misspecified models in both layers, and [Clustered inference, confidence bands, and pretrend tests](inference.md) and the [results reference](../reference/results.md#uncertainty) explain the reported uncertainty.

---

**Reading path, step 11 of 22** · Previous: [Sampling information for the observation model](../methods/sampling-information.md) · Next: [Impute actual changes](change-imputation.md)

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
| 11 | **Choose an observation model** (this page) |
| 12 | [Impute actual changes](change-imputation.md) |
| 13 | [Parallel trends](../methods/parallel-trends.md) |
| 14 | [Add PT-all restrictions](pt-all-comparison.md) |
| 15 | [GMM weighting and covariance regularization](../methods/gmm.md) |
| 16 | [Robustness](robustness.md) |
| 17 | [Inference](inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |
| 21 | [Cross-pair chaining and control parallel trends](../methods/cross-pair-chaining.md) |
| 22 | [Short pre-treatment windows](short-pretreatment-windows.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
