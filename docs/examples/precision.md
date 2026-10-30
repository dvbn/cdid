# How to compare precision

A smaller standard error does not by itself establish a better estimator. Before comparing chained and ordinary mode, control groups, parallel trends choices or corrections for missing outcomes, match the target population and parameter, then state the assumptions that justify each fit. Balance alone decides neither validity nor a precision ranking. [Moving from did](../from-did.md#why-can-balanced-results-differ-and-which-should-i-use) explains how to choose between the two estimation modes.

This page explains covariance and the package's precision diagnostics. It reports no simulation results, executed examples or timing measurements. Executed walkthroughs elsewhere in this guide show individual fits; they are not repeated-sample precision comparisons.

## More comparisons do not mean independent information

Chained estimation combines two-date comparisons. If a vector of link estimates has covariance matrix $\Sigma$ and its coefficients $a$ are fixed, the variance of the linear combination is $a'\Sigma a$. The off-diagonal covariances matter as much as the individual link variances. This is an exact identity for fixed coefficients. Feasible GMM estimates its weights from the same data: inserting fitted, data-dependent coefficients does not give an exact unconditional finite-sample variance identity. Under the appropriate GMM regularity conditions, a first-order calculation uses the limiting combination and complete influence functions. The package's inference follows its complete influence functions, rather than treating fitted weights as an exact finite-sample variance formula.

For the same unit, consecutive changes telescope:

$$
(Y_{i,4}-Y_{i,3})+(Y_{i,3}-Y_{i,2})=Y_{i,4}-Y_{i,2}.
$$

The shared intermediate outcome cancels. If time-specific measurement errors are mutually independent, their contribution to the covariance of these adjacent changes is minus the variance of the shared error. The total covariance need not be negative: serial dependence, heterogeneous trends and other shared components also matter. With different observed units or weights at each link, estimated link averages need not telescope.

In the balanced common reduction case, chained cell estimates agree with `did` up to rounding: built-in IPW, no covariates, no correction for missing outcomes, never treated controls, weights fixed across comparisons, PT-post, and matching anticipation and reporting settings. Splitting the same units' long changes into consecutive steps adds comparison rows without independent information. This point-estimate reduction does not establish agreement of every standard error, band or aggregate. See [the precise comparison](../from-did.md#when-chained-and-did-estimates-should-agree).

## What GMM weighting and standard errors measure

Overlapping comparison paths can make the moment covariance singular even when an ATT is supported. Lack of support is a different issue: if the available comparisons and trend restrictions cannot recover a cell, weighting cannot make it available. Read `fit$chained$support_status` and `support_reason`, or the corresponding columns of `cdid::tidy(fit)`.

With `gmm_weighting = "optimal"`, feasible two-step GMM uses an estimated moment covariance. The package removes proven exact redundancies before checking that covariance. The selected `gmm_covariance = "tikhonov"` handler adds a trace-proportional ridge when the reduced covariance needs regularization; exactly zero covariance has a separate least-squares rule. Regularization affects weighting, supplies no identifying information and guarantees no finite-sample precision improvement. `fit$chained$covariance_handling` records what happened; the [options reference](../reference/options.md#gmm-weighting-and-covariance-regularization) gives the rules.

The inverse weighting matrix is not the reported estimator covariance. With influence functions computed:

- `fit$inffunc` has one row per unit in the unit table, including units without an observed change, and one column per reported cell.
- `fit$V_analytical` estimates the covariance of `sqrt(n) * (estimate - target)`: it is **n times** the estimated covariance of the cell estimates. Analytical standard errors are `sqrt(diag(fit$V_analytical) / fit$n)`, also recorded in `fit$chained$analytical_se`.
- `fit$se` and `std.error` in `cdid::tidy(fit)` use analytical or multiplier standard errors according to the inference actually performed. They need not equal the analytical values when `bstrap = TRUE`.

The complete influence functions account for fitted nuisances and, in aggregation, estimated cohort shares. Clustering sums their rows within clusters for inference; it does not change the unit-based GMM weights or point estimates. Read the effective flags and draws in `fit$chained$inference`, and compare matching pointwise intervals or matching simultaneous families. [Inference](inference.md) and the [results reference](../reference/results.md#uncertainty) explain the endpoints and `cdid::aggte()` outputs.

## Incomplete panels: match the population and sampling conditions

Complete-case analysis, repeated cross sections and chained changes each have settings where they are valid. Their standard errors answer a common precision question only when their targets and assumptions align.

| Approach | What to check before comparing precision |
| --- | --- |
| Complete-case balancing | It uses units observed at every date. Their ATT can differ from the full population ATT when observation relates to heterogeneous effects. The targets can coincide under suitable representativeness conditions; a rotating panel may contain no complete cases. |
| Outcome-level RCS comparisons | Each date's sample must satisfy the level-sampling conditions. Changing composition in an unobserved fixed component can invalidate that comparison even when observed within-unit changes remain representative. Fixed differences alone do not invalidate representative RCS sampling. |
| Chained changes | The observed or reconstructed pair changes must represent the declared population under the chosen assumptions. Overlapping pairs can identify effects without any complete unit; scarce bridges can still make estimates imprecise. |

In `did` 2.5.1, a genuinely unbalanced `panel = TRUE, allow_unbalanced_panel = TRUE` call uses RCS estimation routines internally while retaining its distinction from true repeated cross sections through unit identities. Ordinary mode forwards this behavior. [Moving from did](../from-did.md#incomplete-panels-in-did-and-in-chained-mode) describes the routes, and [Handling missing outcomes](../methods/missing-outcomes.md) states the chained conditions.

## Restrictions can change information and robustness

| Choice | What can change, holding the population target fixed |
| --- | --- |
| Never treated versus not yet treated controls | Under full NYT parallel trends and overlap, the chained catalog can retain additional fixed-horizon control-pool moments for the same pair. Joint GMM can change estimates and uncertainty even with complete outcomes. More valid moments do not establish a universal finite-sample gain. |
| PT-post versus PT-all | Credible extra untreated-trend restrictions from the chosen `pt_start` can add identification or precision. False restrictions can bias estimates; neither bias nor a gain follows from the option label alone. PT-all may add no restriction in a particular calendar. `base_period` separately controls preperiod reporting. |
| Treatment OR or DR | A shared counterfactual change slope is a restriction across comparison blocks. Credible pooling can share information; unjustified pooling can invalidate the OR branch. DR retains its treatment-propensity correctness alternative under the other required conditions. |
| Observation or imputation pooling | `selection_pooling` shares coefficients of observation-probability models; `imputation_pooling` shares regressions of actual cohort changes, including treatment effects. Neither argument controls the shared counterfactual slope of treatment OR/DR. Their distinct restrictions require justification. |

Small pair probabilities can amplify observation weights; rare links and weak covariate support can worsen uncertainty. Neither a large unit table nor a short chain alone determines precision. [Robustness](robustness.md) separates the model-correctness alternatives, and [Parallel trends](../methods/parallel-trends.md) explains which restrictions are active.

## Read a fit before making a ranking

Use `fit$n`, `fit$chained$frame`, `fit$group` and `fit$t` to check the population and cells. Inspect availability through `cdid::tidy(fit)` and `n.requested`/`n.available` in `cdid::glance(fit)`. Compare `DIDparams`, `fit$chained$pt`, `selection`, `covariance_handling` and `inference` for the restrictions and uncertainty actually used. An unavailable effect has no precision ranking.

A single fit's estimated standard error differs from the empirical standard deviation across repeated samples. Bias, RMSE and interval coverage additionally compare estimates with the population truth. A fair repeated-sample evaluation pairs the same inputs and observation masks, matches interval families, retains failed and unavailable attempts, and reports denominators and Monte Carlo uncertainty. A result in one setting cannot establish a universal ranking.

Resource diagnostics are separate. The admitted workspace is recorded in `fit$chained$public$profile$numeric_workspace_bytes`, in bytes. `fit$chained$public$counted_work_limit` is an optional cap on conservative counted-work estimates, default `Inf`; these work units are not memory bytes or measured seconds. Neither is a process-RSS, CPU-time or whole-call resource guarantee. Measured runtime and memory require measurements, with timeouts reported as censored attempts.

The [roadmap](../../README.md#roadmap) states the planned efficiency direction.

---

**Reading path, step 5 of 20** · Previous: [Identify effects when the observed units change composition](identification.md) · Next: [Prepare your panel and unit table](../help/prepare-data.md)

<details>
<summary>All steps of the reading path</summary>

| Step | Page |
| --- | --- |
| 1 | [Get started](../get-started.md) |
| 2 | [Moving from did](../from-did.md) |
| 3 | [Incomplete panels](incomplete-panels.md) |
| 4 | [Identify effects when the observed units change composition](identification.md) |
| 5 | **How to compare precision** (this page) |
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
