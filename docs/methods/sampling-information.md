# Sampling information for the observation model

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

With `selection_model = "pair_logit"`, observation weighting fits one logit per cohort and pair of dates for the probability that a unit is observed at both dates. By default each of these pair logits is fitted on all of the cohort's rows in the unit table. This advanced page explains the smaller fitting sets that `selection_universe` offers, what `selection_rho` and `selection_justification` must then declare, and why a smaller fitting set is a different model rather than a simplification. It expands the advanced section of [Handling missing outcomes](missing-outcomes.md#advanced-fitting-the-observation-model-on-part-of-the-unit-table), which you should read first, and introduces no new option.

The declared target remains the average treatment effects of the study population listed by the unit table; changing the fitting set supplies no validity guarantee.

The identity and hand table below use covariates fixed within units and fixed unit sampling weights. Write W for the full scientific conditioning set needed by the analysis, including relevant treatment, outcome and observation covariates. An observation formula may use fewer variables only when its model remains correct for the required conditioning set. The existing [reference-date covariate and weight options](../reference/options.md#changing-covariates-and-weights) remain available; this fixed-covariate derivation alone supplies no blanket transport guarantee for changing reference contexts.

## Three sets of units

`selection_universe` names one of three nested sets of rows of the unit table:

- **F**, the whole unit table: the study population, one row per unit, including units observed once or never.
- **B**, the units with at least one observed outcome.
- **A**, the units with at least two observed outcomes.

The counts are of observed outcomes, that is, nonmissing values of the outcome column at the dates of the stored calendar. They do not describe structural observation windows: a unit scheduled for two interviews that answers only one belongs to B but not to A. A row of `data` kept for its covariates with a missing outcome adds nothing to the count, while an outcome replaced by zero counts as observed, one more reason to [keep missing outcomes as `NA`](../help/prepare-data.md#never-impute-zero-for-missing-outcomes).

If you omit `frame`, the units present in `data` form F, including units observed once and units whose outcomes are all missing. A unit that was never observed and has no row in `data` is then outside F, and no fitting set brings it back. When every unit in `data` has an observed outcome, F and B are then the same set.

These sets answer two different questions. The first is which units supply observed changes. A change needs two observed outcomes of the same unit, at both dates of a pair, so every observed change comes from A, whatever the fitting set. The second is the question `selection_universe` answers: which units enter the pair logits, as successes (observed at both dates of the pair) or failures (not observed at both). Every success belongs to A, so the three fitting sets differ only in the failures that the logits see. For a pair e, the rows of a cohort play these roles in the cohort's logit:

| Row of the unit table | Observed changes | Under `"F"` | Under `"B"` | Under `"A"` |
| --- | --- | --- | --- | --- |
| observed at both dates of e | the change over e, and over its other observed pairs | success | success | success |
| observed at two or more dates, but not at both dates of e | changes over other pairs only | failure | failure | failure |
| observed at one date | none | failure | failure | not a row |
| no observed outcome, including rows present only in the unit table | none | failure | not a row | not a row |

Whatever its role in the logits, every row keeps its cohort, weight and cluster label and counts in the cohort shares, the cohort rates and the standard errors.

## What the fitting set changes

`selection_universe` chooses the fitting set U, the rows on which each cohort's pair logits are fitted. The observation weights then combine these logits with cohort rates and, under `"B"` or `"A"`, with the declaration described below. The choice does not change the study population, the target effects or the cohort shares, and the influence function keeps one row for every row of the unit table, all counted in the reported `n`. The package determines B and A itself from the outcomes in `data`. Do not trim the unit table or `data` to U: a unit removed from the unit table, or from `data` when you omit `frame`, leaves the study population, and outcomes removed from `data` change which units count as observed.

For the fixed covariates used here, the fitting set also decides where the observation design is read. Under `"F"` they are read on every row of the unit table. Under `"B"` or `"A"`, the observation design and the column named by `selection_rho` are evaluated only on rows of U, so their values may be missing outside U as long as no other active part of the fit reads them there. Under `selection = "conditional_mtar"`, every row of the unit table still needs a known cohort, and a weight and a cluster label when you use them, and the treatment and outcome models keep their own covariate requirements. A missing or nonfinite observation covariate on a row of U is refused with an error. The package does not infer missing attributes, drop rows, or change the fitting set, the model or the target to make a fit succeed.

## Where the choice is available

| Route | Arguments | Fitting set | Models for missing outcomes |
| --- | --- | --- | --- |
| Observation weighting with the pair logit | `selection = "conditional_mtar"`, `selection_est_method = "ipw"`, `selection_model = "pair_logit"` | `"F"` (default), `"B"` or `"A"` | one pair logit per cohort and pair, fitted on the cohort's rows in U |
| Observation weighting with the product of date probabilities | the same, with `selection_model = "independent_period_logit"` | `"F"` only | date logits fitted on the whole unit table |
| Change imputation | `selection_est_method = "imputation"` | `"F"` only, with `selection_rho = NULL` | no observation model; each change model is trained on the units observed at both dates of its pair and evaluated on every row of its cohort in the unit table |
| Augmented correction | `selection_est_method = "dr"`, with `est_method = "dr"` | `"F"` only, with `selection_rho = NULL` | the observation model, fitted on every row of the cohort in the unit table, and the change models of change imputation |

In the first row, the treatment adjustment can be any that observation weighting supports. `"B"` and `"A"` are available only there; other requests are refused with an error before estimation, and no other combination is substituted.

Pooling and the formulas work the same way under every fitting set. `selection_pooling` restricts only the slopes of the pair logits: `"separate"` (the default) gives each cell, a cohort and pair, its own slopes, `"within_cohort"` shares slopes across a cohort's pairs, and `"across_cohorts"` shares one slope vector across all cells. Each cell keeps its own intercept when the observation formula has one. `selection_xformla`, or `xformla` when it is inherited, defines the covariates of the pair logits, and `xformla` alone defines the treatment and outcome models. Pooling and the fitting set are independent choices.

Only the pair logit of observation weighting has a fitting set choice. B and A describe actual outcome counts, not an initial observation schedule or a model for singleton outcome levels. Their representativeness or transport declaration and the mean-change conditions below still need justification. Known pair probabilities, including a design probability of one, are not supported.

## What you declare with B or A

With `"F"`, the fitting set is the study population itself, the ratio defined below equals 1, and `selection_rho` must be `NULL`; supplying a column is refused.

With `"B"` or `"A"`, the logits describe units of U only, and you must state how U relates to the study population. `selection_justification` must be a nonempty character string, and requests without one are refused. What it asserts depends on `selection_rho`:

| `selection_rho` | What the justification asserts |
| --- | --- |
| `NULL` | Within each cohort, the conditioning variables W have the same weighted distribution in U as in the study population, so the ratio below equals 1. Equivalently, membership in U does not depend on these covariates within the cohort. |
| the name of a column of the unit table | The column holds the known ratio defined below, and the justification states where it comes from. |

Representativeness is not automatically true: whether a unit is observed at least once or twice can depend on the same covariates that drive observation of each pair. A transport column is correct only if its source is. The package records the justification as written, and no fit diagnostic can confirm either assertion, because the logits are fitted on U alone and describe observation within U.

For known transport, the column holds, on every row of U,

```
rho_U(W, G)  =  f_P(W | G, U) / f_P(W | G)
```

where W is the full scientific conditioning set, G the cohort and `f_P` for a density under the declared weighted population law: the study population described by the unit table and its fixed sampling weights, the same law that defines the target effects. The density of the fitting set is in the numerator. The column must be numeric, with finite positive values on U, fixed and known rather than estimated from your data, and normalized: as a ratio of two densities, `rho_U` averages to one over the cohort's population. The validator checks that the named column is numeric and its values on U are finite and positive. It cannot establish that those values are the known normalized density ratio; that remains your declaration. It neither estimates nor renormalizes the column. Estimated or unnormalized ratios are not supported. Do not estimate a ratio from your sample, or rescale one until its sample average is one, to meet these requirements: the package would treat the result as known, and the standard errors would omit its estimation uncertainty.

The declaration and the logit for observation within U are separate assumptions, and either can hold while the other fails. The illustration below has a correct logit and a false declaration of representativeness.

## The identity behind the weights

Write `S(e) = 1` when a unit is observed at both dates of pair e. Under `"B"` or `"A"`, each pair logit models

```
q_U(e, W)  =  P(S(e) = 1 | W, G, U)
```

the probability of observing the pair among the cohort's units in U, which is a different probability for A and for B. Let `pi_U(G) = P(U | G)` be the share of the cohort's study population that belongs to U. A unit with `S(e) = 1` has two observed outcomes, so it belongs to both A and B. Under the positivity that weighting requires, conditioning on membership in U and applying Bayes' rule to `P(U | W, G)` gives the identity

```
P(S(e) = 1 | W, G)  =  pi_U(G) * rho_U(W, G) * q_U(e, W)
```

The left side is the probability of observing the pair in the cohort's study population. Under `"F"` the pair logit models it directly; under `"B"` or `"A"` the declaration transports the fitted probability through the right side. Within a cohort, the observation adjustment is proportional to its inverse; total weights also include the fixed unit sampling weights and cohort rate normalizations. The required positivity means that every covariate value of the cohort's population must also be possible in U, so that `rho_U` is positive, and that the pair must have a positive probability of being observed there, so that `q_U` is positive.

`pi_U(G)` is a population quantity in this identity, not an additional public input. The package computes the response and fitting-set rate quantities it needs from the declared unit table, weights and outcome counts; F remains the population and inference frame. Complete influence functions account for estimated rates and logits. The fixed, known ratio contributes no estimation uncertainty.

The identity follows from the definitions alone. It shows how the denominator is assembled, and it does not by itself establish conditional mean trends missing at random, parallel trends or valid treatment models for the pairs the estimator uses.

## An illustration with probabilities chosen by hand

The probabilities below are chosen by hand for one cohort, one pair of dates and a single fixed conditioning variable X with two values, so W = X here. They illustrate the identity; they are not fitted values, package output or a simulation design.

| | X = 0 | X = 1 |
| --- | --- | --- |
| Share of the cohort's study population | 0.5 | 0.5 |
| Probability of belonging to A | 0.2 | 0.6 |
| Share among the cohort's units in A | 0.25 | 0.75 |
| `rho_A(X, G)` | 0.5 | 1.5 |
| `q_A(e, X)` for the pair of dates 1 and 2 | 0.5 | 0.5 |
| Probability of observing the pair in the study population | 0.1 | 0.3 |

The cohort's share in A is `pi_A(G) = 0.5 * 0.2 + 0.5 * 0.6 = 0.4`, so the shares within A are `0.5 * 0.2 / 0.4 = 0.25` and `0.5 * 0.6 / 0.4 = 0.75`, and the ratios are `0.25 / 0.5 = 0.5` and `0.75 / 0.5 = 1.5`. The identity reproduces the last row: `0.4 * 0.5 * 0.5 = 0.1` and `0.4 * 1.5 * 0.5 = 0.3`. These values fit, for example, three stored dates at which each unit of A, whatever its X, is observed only at dates 1 and 2 with probability 0.5, only at dates 1 and 3 with probability 0.25, and only at dates 2 and 3 with probability 0.25.

Within A the pair probability is the same for both values of X, yet in the population it is three times larger for X = 1: a constant `q_U` does not imply a constant population probability. A logit on A describes `q_A` exactly here, and that says nothing about the ratio. Declaring `selection_rho = NULL` would assert shares of 0.5 and 0.5 within A, but they are 0.25 and 0.75. The identity would then give 0.2 for both values, every observed pair would receive the same observation adjustment, and the weighted pairs would keep the covariate shares of the observed pairs, 0.25 and 0.75, instead of the population's 0.5 and 0.5. With the correct column, the observation adjustment factors are in the proportion of 1 / 0.1 to 1 / 0.3 and restore equal shares under the declared weighted population law. A ratio of 1 is justified only when the covariate distributions match, never because the logit for `q_U` fits well.

## Why A is not always simpler

A smaller fitting set changes the probability being modeled, not only the number of rows. Because `q_U(e, W)` is `P(S(e) = 1 | W, G)` divided by `P(U | W, G)`, a logit that is correct on the whole unit table is generally not correct on A or B, and the reverse.

A smaller fitting set can also leave a cell without failures, since moving from F to B to A removes only failures. With two stored dates, A is exactly the set of units observed at both dates, so every row of the A logit is a success and `q_A = 1`. For a nonempty A fitting set this probability is one by definition, but a logit with its own free cell intercept has no finite interior root. The requested fit is refused; pooling slopes keeps the intercepts. This special identity does not provide a supported known-probability interface. In other settings an all-success sample alone does not prove a population probability of one. An all-failure cell with its own free intercept is also refused. Separation, boundary solutions and rank deficiency are refused under every fitting set, and no observations are trimmed. Since version 1.1.0.0, finite inverse observation denominators use a fixed `1e-6` floor; this safeguard does not repair a failed logit or support violation ([Probability safeguards](missing-outcomes.md#probability-safeguards)). The whole unit table faces the same checks, so a fit on F is not guaranteed to succeed either.

Fewer fitting rows guarantee neither a faster fit nor smaller standard errors, and the package does not rank the fitting sets by efficiency or flexibility. A sufficient route under any fitting set retains conditional mean trends missing at random for every pair used given a scientifically sufficient W, positivity, a correct pair logit on U for each needed cell, valid representativeness or known transport, the treatment adjustment conditions, parallel trends, no anticipation and large sample regularity. These are separate premises: a correct logit does not prove transport, and correct transport does not prove either the logit or the outcome conditions. Sharing slopes requires correctness across the cells that share them, and all active fits and catalog links must meet their own regularity and support requirements.

## Declaring the choice in a call

Each list below is an argument template, rather than a stand-alone fit. For executed fits and their output, see [Choose an observation model](../examples/observation-models.md).

Supply one of these sets of arguments in your `cdid::att_gt()` call, together with the data, unit table, formulas and treatment adjustment of your analysis.

```r
# Pair logits fitted on the whole unit table, the default.
observation_args_F <- list(
  selection            = "conditional_mtar",
  selection_est_method = "ipw",
  selection_model      = "pair_logit",
  selection_universe   = "F",
  selection_rho        = NULL
)

# Pair logits fitted on units with at least two observed
# outcomes, declared representative within each cohort.
observation_args_A <- list(
  selection               = "conditional_mtar",
  selection_est_method    = "ipw",
  selection_model         = "pair_logit",
  selection_universe      = "A",
  selection_rho           = NULL,
  selection_justification = paste(
    "Whether a unit has at least two observed outcomes depends only on its",
    "randomly assigned revisit schedule, which is independent of the",
    "scientific conditioning variables within each cohort, so within each cohort the",
    "units observed at least twice have the weighted covariate distribution",
    "of the study population."
  )
)

# Presumes that the unit table passed as `frame` already has a
# column `rho` holding the fixed, known, normalized ratio
# f_P(W | G, A) / f_P(W | G) on every row of A.
observation_args_A_transport <- list(
  selection               = "conditional_mtar",
  selection_est_method    = "ipw",
  selection_model         = "pair_logit",
  selection_universe      = "A",
  selection_rho           = "rho",
  selection_justification = paste(
    "Column rho is the known ratio f_P(W | G, A) / f_P(W | G) under the",
    "weighted population law, taken from the design documentation of the",
    "revisit sample; it is not estimated from these data."
  )
)
```

After a fit with observation weighting, the selection record shows what you declared:

```r
# Requires `fit`, a result of cdid::att_gt() fitted with
# selection = "conditional_mtar" and selection_est_method = "ipw".
str(fit$chained$selection)
```

<details>
<summary>Console output</summary>

```text
List of 23
 $ method                : chr "corollary1"
 $ canonical_assumption  : chr "conditional_mtar"
 $ requested_selection   : chr "conditional_mtar"
 $ selection_model       : chr "independent_period_logit"
 $ selection_pooling     : chr "separate"
 $ link                  : chr "logit"
 $ fitting_universe      : chr "F"
 $ population            : chr "supplied representative frame"
 $ formula               :Class 'formula'  language ~x
  .. ..- attr(*, ".Environment")=<environment: 0x8287fbb60> 
 $ transport             : chr "rho=1"
 $ justification         : NULL
 $ active_cells          :'data.frame':	12 obs. of  5 variables:
  ..$ cell_id : int [1:12] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ kind    : chr [1:12] "date" "date" "date" "date" ...
  ..$ cohort  : num [1:12] 3 3 3 3 3 ...
  ..$ label   : chr [1:12] "k1_t1" "k1_t2" "k1_t3" "k1_t4" ...
  ..$ block_id: int [1:12] 1 2 3 4 5 6 7 8 9 10 ...
 $ n_active_cells        : int 12
 $ n_selection_blocks    : int 12
 $ n_selection_parameters: int 24
 $ empirical_repeat_rates: Named num [1:2] 0.963 0.968
  ..- attr(*, "names")= chr [1:2] "3" "Inf"
 $ assumptions           : chr "Each active pair/cohort requires population conditional mean trend-MAR on F, positivity, used-pair conditional "| __truncated__
 $ if_scale              : chr "phi_F on ALL N roster rows; core psi=(n_A/N)*phi_F with 1/n_A covariance normalization"
 $ unsupported           : chr [1:4] "known pair probabilities (including design p=1)" "custom IPW callbacks without an explicit log_weight_gradient" "estimated or unnormalized transport" "probit"
 $ models                :List of 10
  ..$ n_selection_fits      : int 12
  ..$ n_selection_cells     : int 12
  ..$ n_selection_blocks    : int 12
  ..$ n_selection_parameters: int 24
  ..$ n_regression_fits     : int 1
  ..$ n_regression_blocks   : int 15
  ..$ n_local_intercepts    : int 15
  ..$ n_shared_slopes       : int 1
  ..$ n_supports            : int 30
  ..$ n_treatment_fits      : int 15
 $ outcome_model         : chr "shared one-period outcome-change control-mean regression"
 $ duration_convention   : chr "stored-calendar transition count L_e = i_t - i_s"
 $ treatment_model       : chr "risk-set treatment propensity logit"
```

</details>

Its entries include the observation model, the pooling, the fitting set, the transport, your justification, the stated assumptions and the counts of fitted models. For observation weighting the `method` entry is `"corollary1"`, the legacy record name, with `canonical_assumption = "conditional_mtar"`.

## Related pages

- [Handling missing outcomes](missing-outcomes.md#advanced-fitting-the-observation-model-on-part-of-the-unit-table): the assumptions and models of observation weighting, and the advanced section this page expands.
- [Options reference](../reference/options.md#handling-missing-outcomes): the definitions of `selection_universe`, `selection_rho` and `selection_justification`, and the supported combinations.
- [Prepare your panel and unit table](../help/prepare-data.md#the-unit-table-frame-and-unit-observation-roles): building the unit table and the roles of units with zero, one or repeated observations.
- [Compare model robustness across treatment and missing outcome adjustments](../examples/robustness.md): which models must be correct for each combination.
- [Results reference](../reference/results.md#the-chained-record): the contents of `fit$chained`.

---

**Reading path, step 10 of 20** · Previous: [Handling missing outcomes](missing-outcomes.md) · Next: [Choose an observation model](../examples/observation-models.md)

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
| 10 | **Sampling information for the observation model** (this page) |
| 11 | [Choose an observation model](../examples/observation-models.md) |
| 12 | [Impute actual changes](../examples/change-imputation.md) |
| 13 | [Parallel trends](parallel-trends.md) |
| 14 | [Add PT-all restrictions](../examples/pt-all-comparison.md) |
| 15 | [GMM weighting and covariance regularization](gmm.md) |
| 16 | [Robustness](../examples/robustness.md) |
| 17 | [Inference](../examples/inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
