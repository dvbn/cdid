# Chained GMM weighting and covariance regularization

Results below were produced on October 29, 2026 with `cdid` 1.1.0.0, `did` 2.5.1 and `DRDID` 1.3.0 under R 4.5.2, using the displayed seeds where applicable. These are individual examples; they do not establish repeated-sample coverage.

Chained estimation combines two-period comparison moments across the retained catalog to recover group-time treatment effects. It separately reports direct, unrestricted pre-period comparisons ([Bellego, Benatia, and Dortet-Bernadet 2025](https://doi.org/10.1016/j.jeconom.2024.105783); [Callaway and Sant'Anna 2021](https://doi.org/10.1016/j.jeconom.2020.12.001)). This page explains how the generalized method of moments (GMM) criterion is formed, how exact structural reductions and Tikhonov covariance regularization work under feasible two-step weighting, how moment-covariance singularity differs from parameter identification, and how point-estimate weighting is kept separate from sandwich and multiplier inference.

## Moment conditions and the quadratic objective

Let $\hat{d}$ denote the $M \times 1$ vector of primitive two-period comparison estimates in the retained comparison catalog. Each entry of $\hat{d}$ compares outcome changes between an ordered pair of stored calendar dates for a target cohort against its admissible control pool (never-treated or not-yet-treated units, including any valid fixed-horizon control-pool reuses in the catalog), constructed using the selected treatment adjustment (`est_method`) and missing-outcome correction (`selection` and `selection_est_method`).

Let $\theta$ denote the GMM model parameters: post-treatment effects and any pre-period trend differentials left free by the chosen parallel-trends model. PT-all restricts represented pre-reference differentials; these model parameters are distinct from the direct, unrestricted pre-period cells reported for the pretrend test. Let $d$ denote the vector of population comparison means and $A$ the known $M \times p$ design mapping $\theta$ to them. Under sufficient parallel-trends, observation-reconstruction, overlap and nuisance-model conditions, with the required regularity,

$$
d = A\theta, \qquad \hat{d} \xrightarrow{p} d
$$

This is a population moment restriction and a large-sample statement, not a claim that nuisance-adjusted sample comparisons are unbiased. Given a symmetric positive-semidefinite $M \times M$ weighting matrix $W$, the chained GMM estimator minimizes the quadratic form

$$
Q_W(\theta) = (\hat{d} - A\theta)' W (\hat{d} - A\theta),
$$

and recovers each supported post-treatment target $\tau_j = c_j'\theta$ as a linear contrast of the fitted parameter vector $\hat{\theta}$.

Several properties follow directly from this formulation:

- **Sufficient conditions for recovery.** Validity of every retained comparison under the maintained parallel-trends, observation-selection, overlap and nuisance-model conditions, identification of $c_j$ from the row space of $A$, and standard GMM regularity are sufficient for consistency. They are not necessary in every special case: misspecification can cancel in a target contrast. Such cancellation is not a general robustness guarantee.
- **What $W$ changes when moments are valid.** When the full moment conditions hold and $W$ converges in probability to a positive-definite limit $W_0$ that preserves the rank of $A$, changing $W$ can change the linear combination of overidentifying comparisons used to reconstruct $\tau_j$ and its asymptotic and finite-sample precision. It does not change the population target.
- **No finite-sample unbiasedness or arbitrary-$W$ guarantee.** Feasible two-step GMM estimates $W$ from the same sample used to form $\hat{d}$. Inserting a data-dependent weighting matrix does not establish exact finite-sample unbiasedness, and asymptotic validity does not hold for arbitrary unstable or data-adaptive choices of $W$.
- **Sensitivity to incompatible moments.** If the population limits of the retained comparisons cannot all be represented by one $A\theta$ (for example after invalid trend or observation restrictions), the moments are incompatible. Misspecification alone need not imply incompatibility or bias in every contrast. Under the relevant convergence conditions, the limiting criterion defines a weighted pseudo-target that generally depends on the limiting weights.
- **Moment-covariance ridge is not a parameter penalty.** When the estimated moment covariance is regularized before inversion to form $W$, the positive ridge is added to the *moment covariance matrix*, not to $A'WA$ and not as an $L_2$ penalty on $\theta$. It imposes no shrinkage target on ATT and supplies no identifying restrictions on $\theta$. A particular estimate can nevertheless move closer to or farther from zero when its comparison weights change.

## Public weighting and covariance options

Two arguments of `cdid::att_gt()` configure the GMM criterion in chained mode, supplemented by two optional tuning entries in `control`:

| Argument or `control` entry | Default | Accepted values | Meaning |
| --- | --- | --- | --- |
| `gmm_weighting` | `"optimal"` | `"optimal"`, `"two_step"`, `"identity"` | `"optimal"` selects feasible two-step GMM using the unit-level moment covariance; `"two_step"` is an explicit alias for `"optimal"`; `"identity"` selects explicit unweighted GMM in the full catalog metric. |
| `gmm_covariance` | `"tikhonov"` | `"tikhonov"` | Singular or near-singular moment-covariance handler under feasible two-step weighting. `"moore_penrose"` is deferred and refused. |
| `control$condition_trigger` | `1e-10` | finite numeric scalar in $(0, 1)$ | Reciprocal spectral condition threshold $\tau_{\mathrm{cond}}$ below which `"tikhonov"` covariance regularization is triggered. |
| `control$ridge_relative` | `1e-4` | finite positive numeric scalar | Relative ridge multiplier $\alpha_{\mathrm{ridge}}$ applied to the full-catalog trace scale when regularization is triggered. |

Supplying `"two_step"` selects the same feasible two-step estimator as `"optimal"` and records in `fit$chained$covariance_handling` that the alias was used. Passing `gmm_covariance = "moore_penrose"` is refused by the request validator with class `"cdid_chained_request_input_error"`. Unsupported cutoff entries in `control`, such as `mp_cutoff`, instead raise `"cdid_chained_public_control_error"`. Both inherit `"cdid_error"`; Moore-Penrose pseudoinverse weighting is deferred and has no supported interface.

## Exact structural reductions and the full comparison metric

Before examining the sample moment covariance, `cdid` applies proved exact structural reductions when proved linear identities reconstruct the full catalog from retained comparisons (for example, on baseline built-in IPW routes where distinct candidate rows share the exact same cohort, control pool, and observation mask).

Suppose the full catalog contains $M$ comparisons and proved exact equivalences express the $M \times 1$ residual vector in terms of $L \le M$ retained comparison rows through a sparse $M \times L$ reconstruction matrix $E$:

$$
\hat{d} - A\theta = E(\hat{d}_r - A_r\theta), \qquad G = E'E = C'C,
$$

where $G$ is the symmetric positive-definite $L \times L$ Gram matrix of the reduction and $C$ is its upper-triangular Cholesky factor. When no structural reduction applies, $L = M$ and $G = C = I_M$.

Preserving $G = E'E$ is essential for three reasons:

1. **Exact preservation of the original $M$-comparison objective.** Any full-catalog weighting matrix becomes $E'WE$ on the $L$ retained rows. In particular, under `gmm_weighting = "identity"`, the full-catalog unweighted objective is
   $$
   (\hat{d} - A\theta)' I_M (\hat{d} - A\theta) = (\hat{d}_r - A_r\theta)' E'E (\hat{d}_r - A_r\theta) = \|C(\hat{d}_r - A_r\theta)\|_2^2.
   $$
   The reduced solver therefore whitens the retained system by $C$ so that `"identity"` still minimizes the exact sum of squares over all $M$ original catalog comparisons, rather than replacing $G$ by a naive $I_L$ in reduced coordinates.
2. **No arbitrary row dropping or spanning-tree subsetting.** Correlation alone does not justify removing comparison rows; an exact reduction needs a proved reconstruction identity, which can involve linear combinations rather than identical rows. Discarding correlated rows, deduplicating rows without preserving $G = E'E$, or restricting the catalog to a minimal identification forest would alter the estimator and is never performed.
3. **Algebraic identity is not extra sample information.** When two catalog entries evaluate the exact same sample difference on the same units, or when consecutive changes of the same balanced units telescope to a long change, the extra row repeats existing sample variation rather than adding independent information. Exact reduction removes structural duplication without changing the criterion.

## Covariance conditioning and Tikhonov regularization

Under `gmm_weighting = "optimal"` (or `"two_step"`), let $\widehat{V}_r$ denote the $L \times L$ comparison covariance on the solver's weighting scale, formed from complete unit-level link influence functions. This unit-IF/asymptotic normalization need not equal the estimated covariance of $\hat d_r$ itself. A common positive rescaling of $\widehat V_r$ rescales its trace-proportional ridge equally and multiplies the criterion by a constant, leaving its minimizer unchanged. To respect the original $M$-comparison geometry when $L < M$, conditioning and regularization are defined in the $C$-transformed coordinates where the full-catalog identity metric $G = C'C$ becomes $I_L$. Define the transformed $L \times L$ symmetric covariance matrix

$$
B = C \widehat{V}_r C'.
$$

Adding an isotropic ridge $\lambda I_M$ in the original $M$-comparison coordinates corresponds in the reduced transformed coordinates to replacing $B$ by $B + \lambda I_L = C(\widehat{V}_r + \lambda G^{-1})C'$. The effective reduced weighting matrix is therefore

$$
W_r(\lambda) = C'(B + \lambda I_L)^{-1} C = \left(\widehat{V}_r + \lambda G^{-1}\right)^{-1},
$$

which reduces to $(\widehat{V}_r + \lambda I_M)^{-1}$ when $G = I_M$, and otherwise scales the ridge inversely with $G$ so that identical original comparisons receive the same metric weight regardless of how exact reduction grouped them.

### How the handler chooses a branch

After checking that $\widehat{V}_r$ is finite and exactly symmetric and forming a symmetric finite $B$, the estimator computes a symmetric eigendecomposition of $B$ with minimum and maximum eigenvalues $d_{\min}(B)$ and $d_{\max}(B)$, checks that the eigendecomposition residual lies within floating-point roundoff bounds, and evaluates the reciprocal spectral condition ratio

$$
\rho = \frac{\max(d_{\min}(B), 0)}{d_{\max}(B)}
$$

(with $\rho = 0$ if $B = 0$). Materially negative eigenvalues or diagonal entries abort the fit; small negative eigenvalues consistent with floating-point roundoff are accepted only when $\widehat{V}_r$ carries authenticated unit-IF Gram construction provenance. This authenticates the covariance construction for the roundoff diagnosis; it is not a per-fit numerical error certificate.

If the initial eigensystem exceeds the reconstruction residual bound, the implementation attempts one QR reorthogonalization with Rayleigh values. It accepts the refined eigensystem only when it meets the same bound, `100 * L * .Machine$double.eps * ||B||_F`. The covariance matrix, definiteness checks and ridge policy are unchanged. `fit$chained$covariance_handling$conditioning$eigensystem_refinement` records whether this refinement was used and the original residual. A failed refinement still refuses the fit.

The weighting step then follows one of four branches, recorded in `fit$chained$covariance_handling`:

| Branch | Condition | Effective reduced weight and solve |
| --- | --- | --- |
| **Well-conditioned unregularized solve** | $B$ is positive definite, $\rho \ge \tau_{\mathrm{cond}}$, retained dimension $L$ does not exceed the sample rank bound from contributing unit rows, and ordinary Cholesky factorization of $\widehat V_r$ and the subsequent solve pass their checks | $W_r(0) = C' B^{-1} C = \widehat{V}_r^{-1}$, implemented by Cholesky whitening of $\widehat V_r$ and pivoted QR factorization on the whitened design without invoking the singular-covariance handler |
| **Tikhonov covariance ridge** | $B \ne 0$ and $\rho < \tau_{\mathrm{cond}}$, $d_{\min}(B) = 0$, $d_{\min}(B) < 0$ within authenticated Gram roundoff, $L$ exceeds the sample rank upper bound, or unregularized Cholesky of $\widehat V_r$ fails after valid covariance diagnosis | $W_r(\lambda) = C'(B + \lambda I_L)^{-1}C$ with $s = \mathrm{tr}(B)/M$ and $\lambda = \alpha_{\mathrm{ridge}} s$, solved by Cholesky factorization of $B + \lambda I_L$ and pivoted QR |
| **Exact-zero covariance limit** | $\widehat{V}_r = 0$ and $B = 0$ identically | Separate scale-free full-catalog least-squares rule with $\lambda = 0$ and reduced weight $G = C'C$ (whitened by $C$), paired with the actual zero sandwich covariance and zero-scale inference diagnostics |
| **Refusal on invalid covariance or failed solve** | $\widehat{V}_r$ or $B$ is nonfinite or materially indefinite, or the Cholesky or whitened-QR factorization fails its conditioning or residual bounds | Call stops with an error inheriting `"cdid_error"`; no estimate is returned |

Note three implementation details of this rule:

- **Trace scale uses the full catalog dimension $M$.** In the Tikhonov branch, the scale is $s = \mathrm{tr}(B) / M$, where $M$ is the full catalog comparison count and $\mathrm{tr}(B) = \mathrm{tr}(\widehat{V}_r G)$ equals the trace of the full $M \times M$ moment covariance.
- **Stable factorizations rather than explicit matrix inversion.** Neither the unregularized branch nor the Tikhonov branch forms an explicit dense inverse $B^{-1}$ or $(B + \lambda I_L)^{-1}$. The unregularized branch whitens with the Cholesky factor of $\widehat V_r$; the ridge branch uses that of $B+\lambda I_L$. Both then apply pivoted QR to the whitened design, and verify factor conditioning and relative residuals before accepting the solution.
- **No silent fallbacks or ad hoc repairs.** The default values `condition_trigger = 1e-10` and `ridge_relative = 1e-4` are numerical implementation settings, not universal optimal tuning parameters and not thresholds for parameter identification. If $B + \lambda I_L$ or the whitened QR factor still fails its conditioning guards, `cdid` refuses the fit rather than increasing $\lambda$ iteratively, clipping eigenvalues, adding random jitter, dropping comparisons, or switching automatically to `gmm_weighting = "identity"`.

## Covariance singularity versus parameter identification

Singular or near-singular moment covariance $\Sigma = \mathrm{Cov}(\hat{d})$ and non-identification of a target parameter $\tau = c'\theta$ are distinct mathematical phenomena. Moment covariance is a property of the stochastic errors across comparison rows. Parameter identification is a property of the design matrix $A$ and the maintained population restrictions $d = A\theta$.

The two analytical examples below illustrate why neither implies the other. Here $\hat d$ denotes noisy comparisons and $d = E[\hat d]$ their population means. These are algebraic examples, not executed R outputs.

### Example 1: Singular moment covariance with an identified parameter

Suppose a scalar treatment effect $\theta$ is measured by two comparison moments that carry the exact same zero-mean shock $u$ with $\mathrm{Var}(u) = \sigma^2 > 0$:

$$
\hat d_1 = \theta + u, \qquad \hat d_2 = \theta + u, \qquad A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}, \qquad \Sigma = \sigma^2 \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}.
$$

Here the design matrix $A$ has column rank $1$, so $\theta = E[\hat d_1] = E[\hat d_2]$ is identified. Yet the $2 \times 2$ covariance matrix $\Sigma$ has rank $1$ and determinant $0$ because $\hat d_1 - \hat d_2 = 0$ with probability one. In chained estimation, a singular or near-singular covariance of this kind can arise when duplicate or telescoping comparison paths span the same dates. If the duplication is an exact structural equivalence, reduction collapses the pair to $L = 1$ with $G = E'E = 2$ and positive scalar variance. If this duplicated pair is retained without a structural reduction, Tikhonov regularization replaces $\Sigma$ with $\Sigma + \lambda I_2$, which has eigenvalues $2\sigma^2 + \lambda$ and $\lambda > 0$, allowing the GMM step to form a well-defined convex combination of $\hat d_1$ and $\hat d_2$. Because both moments share the same shock $u$, the second moment adds no independent information, and the estimator variance remains $\sigma^2$.

### Example 2: Nonsingular moment covariance with an unidentified parameter

Suppose instead that two parameters $\theta = (\theta_1, \theta_2)'$ (for instance, an unrestricted pre-period trend differential $\theta_1$ and a post-treatment effect $\theta_2$) enter a single available comparison moment with $\mathrm{Var}(u) = \sigma^2 > 0$:

$$
\hat d_1 = \theta_1 + \theta_2 + u, \qquad A = \begin{pmatrix} 1 & 1 \end{pmatrix}, \qquad \Sigma = \begin{pmatrix} \sigma^2 \end{pmatrix}.
$$

Here the $1 \times 1$ moment covariance $\Sigma = (\sigma^2)$ is strictly positive definite and well conditioned, but $A$ has rank $1 < 2$. The moment identifies the sum $\theta_1 + \theta_2$, while the individual contrast $\tau = \theta_2 = (0, 1)\theta$ is not in the row space of $A$. Adding a covariance ridge $\lambda > 0$ changes the scalar weight from $1/\sigma^2$ to $1/(\sigma^2 + \lambda)$, so the normal matrix

$$
A'(\Sigma + \lambda)^{-1}A = \frac{1}{\sigma^2 + \lambda} \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}
$$

still has rank $1$ and cannot be inverted for $\theta_1$ and $\theta_2$ separately. Because covariance regularization acts on $\Sigma$ rather than adding rows to $A$ or penalizing $\theta$, it never turns an unidentified parameter into an estimable one.

### Catalog support, model refusals, and precision rankings

In `cdid`, post-treatment target availability is determined from the catalog design $A$ before and independently of covariance regularization. Reported pre-period cells retain their separate direct-comparison support rules:

- **Unsupported catalog cells (`NA`).** When a requested target contrast $c_j'\theta$ does not lie in the row space of $A$ under the retained comparison catalog and the chosen `pt_assumption`, the cell is marked unavailable (`available = FALSE`), its estimate and standard error are set to `NA`, and `fit$chained$support_status` and `fit$chained$support_reason` record why the catalog cannot recover it. This is a statement about support in the retained comparison catalog, not a claim that the effect is unidentifiable in the population under every possible assumption.
- **Computation and nuisance refusals (errors).** Separation or unresolved rank/support failure in a required model after its allowed preprocessing, or a covariance or QR factorization that fails the applicable checks, stops the call with an error rather than returning `NA` for affected cells or dropping comparisons.
- **No uniform precision dominance.** With independent units, valid moments, consistent covariance estimation, standard regularity and a nonsingular limiting covariance after exact reductions, the unregularized optimal branch has the usual asymptotic efficiency property within the retained moment set. That statement does not automatically extend to a ridge branch or to dependent units. Covariance ridge changes the weighting to obtain numerical conditioning when $B$ is singular or ill-conditioned, and `"identity"` uses fixed full-catalog weights. In finite samples, estimation noise in $\widehat{V}_r$ can offset theoretical efficiency gains, so neither `"optimal"`, nor a particular `ridge_relative`, nor `"identity"` is guaranteed to produce smaller finite-sample variance or mean squared error in every setting.

## Sandwich inference, multiplier bootstrap, and clustering

The regularized weighting matrix $W_r(\lambda)$ is used only to compute the point estimates $\hat{\theta}$ and the resulting linear target loadings. It is **never** inverted to report parameter variance or standard errors.

### Complete unit-level influence functions and `V_analytical`

For a supported post-treatment target, the GMM solver computes the realized loading vector on the retained comparison moments and composes it with their complete influence functions. Separately reported pre-period cells use their direct auxiliary comparisons, not restricted GMM pretrend parameters. The reporting system retains complete original-unit influence functions across all $n$ rows of the unit table (`fit$n` and `fit$chained$frame`). These rows include frame-only units and units with fewer than two observed outcomes, and carry the uncertainty of every active nuisance channel: cohort and repeat-observation rates, observation-selection logits, actual-change imputation fits, treatment propensity logits, shared counterfactual outcome regressions, and (in `cdid::aggte()`) estimated cohort shares.

From the resulting $n \times J$ influence-function matrix (`fit$inffunc`), the package forms the sandwich covariance matrix:

- `fit$V_analytical` estimates the asymptotic covariance of $\sqrt{n}(\hat{\tau} - \tau)$ and equals **$n$ times** the estimated covariance matrix of the cell estimates.
- Analytical standard errors are `sqrt(diag(fit$V_analytical) / fit$n)`, also stored in `fit$chained$analytical_se`.
- When bootstrap inference is actually performed, `fit$se` and `std.error` in `cdid::tidy(fit)` use multiplier interquartile-range standard errors on the estimate scale; otherwise, with inference computed, they use the analytical values. Requested flags alone do not establish what ran: inspect `fit$chained$inference` for the effective flags, kind and draws. Empty or exactly zero families may need no draws.
- With `compute_inffunc = FALSE`, the covariance work needed for GMM weighting is retained, so `fit$att` equals that under `TRUE`. Returned `inffunc`, `V_analytical`, `W` and `Wpval` are `NULL`, `se` is a numeric vector of `NA`, effective `bstrap` and `cband` are `FALSE`, and the inference kind is `"not_computed"`. `aggte()` refuses this point-only result.

For `att_gt()` cells, `bstrap = FALSE, cband = TRUE` gives analytical standard errors and pointwise normal intervals, with no draws. Group, dynamic and calendar aggregate bins instead retain a multiplier band with analytical standard errors and a warning under those flags. Overall summaries and simple aggregation use pointwise normal intervals. The [results reference](../reference/results.md#uncertainty) specifies each endpoint; separate families do not form one joint band.

### Clustering and asymptotic regularity

When `clustervars` supplies a cluster variable (at most one variable in addition to `idname`), rows of `fit$inffunc` are summed within clusters before computing `fit$V_analytical` and drawing cluster-level Rademacher signs in the multiplier bootstrap. Clustering affects inference only: point estimates `fit$att` and the GMM weighting matrix $\widehat{V}_r$ remain unit-based. Consequently, `"optimal"` weighting is unit-level feasible two-step GMM and is not claimed to be cluster-optimal.

Treating feasible two-step or regularized GMM weights as fixed in the first-order influence function relies on standard large-sample regularity: valid moment conditions, asymptotically linear complete link estimators, consistent covariance estimation, and a stable limiting loading branch (either a well-conditioned positive-definite limit covariance or a fixed limiting regularization mapping). The package claims no finite-sample cancellation theorem, no uniform inference validity near covariance singularity boundaries or with few clusters, and no per-fit numerical error certificate (passing `control = list(relative_accuracy = ...)` is refused with class `"cdid_chained_public_accuracy_deferred"`). For the full statistical discussion of standard errors, simultaneous confidence bands, and pretrend tests, see [How to compare precision](../examples/precision.md), [Inference](../examples/inference.md), and the [results reference](../reference/results.md#uncertainty).

## Inspecting GMM weighting, tuning, and resource controls

Changing `gmm_weighting` or overriding `condition_trigger` and `ridge_relative` in `control` is an explicit estimator choice recorded in the returned object. It is not an automatic fix for an unidentified effect, a violated parallel-trends restriction, or an invalid missing-outcome or treatment specification.

The snippet below reuses `sim_panel` from [Get started](../get-started.md#prepare-the-panel); run that preparation step first.


```r
# Helper comparing default feasible two-step GMM, explicit identity weighting,
# and an explicit covariance-regularization / resource control override on the
# same prepared panel (and optional unit frame).
fit_gmm_variants <- function(data, frame = NULL) {
  fit_opt <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "first_treat",
    data = data, frame = frame,
    control_group = "nevertreated", est_method = "dr",
    gmm_weighting = "optimal", gmm_covariance = "tikhonov",
    bstrap = FALSE, cband = FALSE
  )

  # "two_step" is a documented alias for "optimal"
  fit_alias <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "first_treat",
    data = data, frame = frame,
    control_group = "nevertreated", est_method = "dr",
    gmm_weighting = "two_step", gmm_covariance = "tikhonov",
    bstrap = FALSE, cband = FALSE
  )

  # Explicit unweighted GMM in the full catalog metric
  fit_id <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "first_treat",
    data = data, frame = frame,
    control_group = "nevertreated", est_method = "dr",
    gmm_weighting = "identity",
    bstrap = FALSE, cband = FALSE
  )

  # Explicit covariance-ridge tuning and optional resource limits
  fit_tuned <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "first_treat",
    data = data, frame = frame,
    control_group = "nevertreated", est_method = "dr",
    gmm_weighting = "optimal", gmm_covariance = "tikhonov",
    control = list(
      condition_trigger = 1e-8,
      ridge_relative = 1e-3,
      numeric_workspace_bytes = 1073741824,
      counted_work_limit = Inf
    ),
    bstrap = FALSE, cband = FALSE
  )

  list(optimal = fit_opt, two_step = fit_alias, identity = fit_id, tuned = fit_tuned)
}

fits <- fit_gmm_variants(sim_panel)

# Inspect the recorded GMM weighting, support, inference, and public profile records
covariance <- fits$optimal$chained$covariance_handling
covariance$target_availability$target_id <- NULL
str(covariance)
str(fits$optimal$chained$support)
str(fits$optimal$chained$inference)
profile <- fits$optimal$chained$public
profile$decisions <- NULL
profile$profile_status <- NULL
profile$numerical_policy$decision <- NULL
profile$numerical_policy$relative_accuracy_status <- NULL
str(profile)

# Compare cell estimates and analytical standard errors
cdid::tidy(fits$optimal)[, c("target.id", "role", "available", "estimate", "std.error")]
cdid::tidy(fits$identity)[, c("target.id", "role", "available", "estimate", "std.error")]
cdid::glance(fits$optimal)
```

<details>
<summary>Console output</summary>

```text
List of 28
 $ weighting_requested   : chr "optimal"
 $ weighting_canonical   : chr "optimal"
 $ alias_used            : logi FALSE
 $ handler_selected      : chr "tikhonov"
 $ handler_invoked       : logi FALSE
 $ handler_used          : chr "none"
 $ trigger_reason        : chr "well_conditioned_pd"
 $ dimensions            :List of 9
  ..$ M_full                   : int 15
  ..$ L_retained               : int 5
  ..$ p_parameters             : int 5
  ..$ k_design_rank            : int 5
  ..$ n_targets                : int 5
  ..$ n_estimable              : int 5
  ..$ M_positive_rows          : int 2000
  ..$ sample_rank_upper_bound  : int 1999
  ..$ exceeds_sample_rank_bound: logi FALSE
 $ gram_provenance       : chr "sparse_structural_reduction_E_t_E"
 $ gram_is_identity      : logi FALSE
 $ gram_is_diagonal      : logi FALSE
 $ gram_n_triplets       : int 20
 $ kept_indices          : int [1:5] 1 2 3 4 5
 $ condition_trigger     : num 1e-10
 $ ridge_relative        : num 1e-04
 $ scale                 : num 7.89
 $ trace_B               : num 118
 $ lambda                : num NA
 $ covariance_rank       : int 5
 $ covariance_rank_rule  : chr "computed eigenvalues above b + r; numerical diagnostic only"
 $ weighted_design_rank  : int 5
 $ conditioning          :List of 13
  ..$ rho                          : num 0.891
  ..$ rho_raw                      : num 0.891
  ..$ d_min                        : num 22.3
  ..$ d_max                        : num 25.1
  ..$ eigenvalues_B                : num [1:5] 25.1 24.7 23.4 23 22.3
  ..$ norm_B_frobenius             : num 53
  ..$ b_roundoff_bound             : num 5.88e-12
  ..$ r_eigen_residual             : num 2.49e-14
  ..$ eigensystem_refinement       :List of 2
  .. ..$ method           : chr "none"
  .. ..$ original_residual: num 2.49e-14
  ..$ diagnosed_numerical_negative : logi FALSE
  ..$ negative_eigenvalue_diagnosis: chr "nonnegative_spectrum"
  ..$ nonpositive_eigenvalue_count : int 0
  ..$ gram_provenance_certified    : logi TRUE
 $ zero_covariance       : logi FALSE
 $ zero_covariance_branch: chr "not_applicable"
 $ inference_claim       : chr "complete_unit_if_sandwich_conditional_on_stable_loading"
 $ target_availability   :'data.frame':	5 obs. of  2 variables:
  ..$ estimable: logi [1:5] TRUE TRUE TRUE TRUE TRUE
  ..$ status   : chr [1:5] "estimable" "estimable" "estimable" "estimable" ...
 $ tuning_overrides      : chr(0) 
 $ scope_note            : chr "Regularized criterion covariance is used only for the GMM weighting matrix; inference replays the complete unit"| __truncated__
List of 6
 $ status                : chr "all_requested_supported"
 $ available_ids         : chr [1:5] "pre(g=3,t=2)" "ATT(g=3,t=3)" "ATT(g=3,t=4)" "ATT(g=3,t=5)" ...
 $ unavailable_ids       : chr(0) 
 $ n_available           : int 5
 $ n_requested           : int 5
 $ whole_request_complete: logi TRUE
List of 16
 $ kind                    : chr "analytic"
 $ bstrap                  : logi FALSE
 $ cband                   : logi FALSE
 $ biters                  : int NA
 $ clustered               : logi FALSE
 $ n_clusters              : int 2000
 $ n_zero_eligible_clusters: int 0
 $ n_total                 : int 2000
 $ rng                     : NULL
 $ multiplier              : NULL
 $ requested               :List of 3
  ..$ bstrap: logi FALSE
  ..$ cband : logi FALSE
  ..$ biters: num 1000
 $ zero_law                : logi FALSE
 $ flag_law                : chr "explicit_bstrap_and_cband"
 $ released_false_true     : logi FALSE
 $ endpoint                : chr "cells"
 $ draws                   : logi FALSE
List of 21
 $ entry                              : chr "cdid::att_gt(mode = \"chained\", est_method = \"dr\", control = ...)"
 $ profile                            :List of 3
  ..$ numeric_workspace_bytes: num 1.07e+09
  ..$ counted_work_limit     : num Inf
  ..$ workers                : num 1
 $ profile_source                     : chr "NULL: resolved to the defaults (optional work cap disabled; memory/worker defaults)"
 $ overrides                          : chr(0) 
 $ numeric_workspace_entries          : num 1.34e+08
 $ counted_work_limit                 : num Inf
 $ work_cap                           :List of 6
  ..$ status                        : chr "disabled"
  ..$ user_limit                    : num Inf
  ..$ source                        : chr "default: no user work budget"
  ..$ effective_private_max_work_ops: num 9.01e+15
  ..$ internal_arithmetic_ceiling   : num 9.01e+15
  ..$ meaning                       : chr "no user work budget; the private stage controls carry the exact-count arithmetic ceiling 2^53 - 1 so that count"| __truncated__
 $ chunk_rows                         : num 256
 $ chunk_rule                         : chr "min(n_frame, 256, floor(workspace_entries / (2048 (links_bound + p + 1))))"
 $ dimensions                         :List of 22
  ..$ n_frame              : num 2000
  ..$ n_panel_units        : num 2000
  ..$ n_rows               : num 12000
  ..$ n_dates              : num 6
  ..$ n_cohorts            : num 2
  ..$ n_targets            : num 1
  ..$ n_covariates         : num 1
  ..$ biters               : num 1000
  ..$ n_clusters           : num 2000
  ..$ n_aggregation_types  : num 4
  ..$ pairs_per_target     : num 15
  ..$ R_full               : num 4
  ..$ candidates_bound     : num 105
  ..$ links_bound          : num 30
  ..$ Q_bound              : num 12
  ..$ params_bound         : num 61
  ..$ J_bound              : num 82
  ..$ mask_entries_bound   : num 60000
  ..$ pool_entries_bound   : num 60
  ..$ E_entries_bound      : num 180
  ..$ materialized_if_bound: num 164000
  ..$ scope                : chr "actual counts of the one prepared input and finite derived upper bounds implied by the accepted catalog/design "| __truncated__
 $ derived_bounds                     :'data.frame':	16 obs. of  2 variables:
  ..$ bound: chr [1:16] "owner/fitted/point/problem/link/covariance max_retained/matrix/output entries" "work limits (max_work_ops, catalog max_work): user work cap disabled; exact-count arithmetic ceiling 2^53 - 1" "max_observer_units" "catalog max_candidates" ...
  ..$ value: num [1:16] 1.34e+08 9.01e+15 2.00e+03 1.05e+02 3.00e+01 ...
 $ glm                                :List of 2
  ..$ epsilon: num 1e-12
  ..$ maxit  : int 50
 $ unused_private_certificate_settings:List of 3
  ..$ taylor_degree: int 12
  ..$ max_abs_eta  : num 50
  ..$ max_squarings: int 16
 $ gmm_guard                          :List of 8
  ..$ max_orthogonality_defect: num 1e-08
  ..$ min_certified_lambda_rel: num 1e-10
  ..$ max_factor_residual_rel : num 1e-10
  ..$ min_rcond_R             : num 1e-10
  ..$ min_rcond_U             : num 1e-10
  ..$ min_sigma_rel_A         : num 1e-10
  ..$ max_kkt_residual_rel    : num 1e-10
  ..$ max_loading_residual_rel: num 1e-10
 $ gmm_weighting                      :List of 3
  ..$ requested : chr "optimal"
  ..$ canonical : chr "optimal"
  ..$ alias_used: logi FALSE
 $ covariance_strategy                :List of 7
  ..$ method             : chr "tikhonov"
  ..$ condition_trigger  : num 1e-10
  ..$ ridge_relative     : num 1e-04
  ..$ overrides          : chr(0) 
  ..$ weighting_requested: chr "optimal"
  ..$ weighting_canonical: chr "optimal"
  ..$ alias_used         : logi FALSE
  ..- attr(*, "class")= chr "cdid_gmm_covariance_strategy"
 $ numerical_policy                   :List of 5
  ..$ mode                                    : chr "conventional"
  ..$ complete_per_fit_certification_requested: logi FALSE
  ..$ complete_per_fit_certification_performed: logi FALSE
  ..$ override_recorded                       : logi FALSE
  ..$ rule                                    : chr "stable Cholesky/QR with existing conditioning and factor residual policies; complete IF/inference; residual rat"| __truncated__
 $ admission                          :List of 8
  ..$ numeric_payload_bytes      : num 1.07e+09
  ..$ numeric_payload_entries    : num 1.34e+08
  ..$ counted_work_limit         : num Inf
  ..$ work_cap_status            : chr "disabled"
  ..$ private_max_work_ops       : num 9.01e+15
  ..$ scope                      : chr "bounds the reachable numeric payload the accepted preflights count; the estimated-work cap is optional and, whe"| __truncated__
  ..$ whole_process_RSS_guarantee: logi FALSE
  ..$ measurement_status         : chr "no whole-call resource guarantee or performance claim"
 $ execution                          :List of 3
  ..$ declared_workers  : num 1
  ..$ multiplier_workers: num 1
  ..$ scope             : chr "estimation is serial; only multiplier-bootstrap raw-draw products may use worker processes, whose whole buffers"| __truncated__
 $ witness                            : chr "not requested/performed on the conventional-computation public path"
 $ retention                          : chr "no observation, moment or record is dropped on refusal; the accepted pipeline's typed bound errors stop the call"
     target.id role available    estimate  std.error
1 pre(g=3,t=2)  pre      TRUE -0.05831968 0.06369725
2 ATT(g=3,t=3) post      TRUE  1.00962257 0.06316956
3 ATT(g=3,t=4) post      TRUE  1.32743963 0.06163755
4 ATT(g=3,t=5) post      TRUE  1.48374276 0.06297661
5 ATT(g=3,t=6) post      TRUE  1.84780436 0.06303371
     target.id role available    estimate  std.error
1 pre(g=3,t=2)  pre      TRUE -0.05831968 0.06369725
2 ATT(g=3,t=3) post      TRUE  1.00962257 0.06316956
3 ATT(g=3,t=4) post      TRUE  1.32743963 0.06163755
4 ATT(g=3,t=5) post      TRUE  1.48374276 0.06297661
5 ATT(g=3,t=6) post      TRUE  1.84780436 0.06303371
  nobs ngroup ntime control.group est.method                 estimator n.requested n.available
1 2000      1     6  nevertreated         dr chained doubly robust GMM           5           5
                          band.family inference clustered
1 pointwise_normal_over_supported_ids  analytic     FALSE
```

</details>

### Resource controls versus statistical and numerical limits

The advanced `control` list also carries resource-admission settings that are recorded in `fit$chained$public` and kept distinct from covariance conditioning and statistical inference:

- **Counted numeric workspace (`numeric_workspace_bytes`).** Defaults to `1073741824` (1 GiB, recorded in `fit$chained$public$profile$numeric_workspace_bytes`) and must be a finite whole number of bytes at least `8` (up to $2^{53} - 1$). Preflight checks convert this byte budget into binary64 array entries (`floor(numeric_workspace_bytes / 8)`) and refuse stages whose conservative reachable numeric buffer envelope exceeds it.
- **Optional counted-work limit (`counted_work_limit`).** Defaults to `Inf` (disabled). Supplying a finite positive whole number (up to $2^{53} - 1$) opts into a preflight cap on conservative counted-work estimates, refusing a stage before execution if its estimated work units exceed the limit; supplying `Inf` explicitly opts out. When disabled (`Inf`), internal stage constructors still use the exact binary64 integer ceiling $2^{53} - 1$ to protect integer count arithmetic from overflow; `fit$chained$public$counted_work_limit` and `fit$chained$public$work_cap` record the user setting and the internal arithmetic ceiling separately.
- **Streamed chunking and bootstrap workers (`chunk_rows`, `workers`).** `chunk_rows` overrides the automatic unit-chunk size used by streamed producers, and `workers` (default `1`) declares how many worker processes `pl = TRUE, cores = k` may use during multiplier-bootstrap draw products. GMM estimation and covariance handling always run serially, and worker buffers are charged against the same `numeric_workspace_bytes` limit.
- **Not process memory or runtime.** Counted workspace bytes measure bounded numeric array payloads, and counted work units are conservative operation counts rather than bytes or seconds. Neither setting caps R process resident set size (RSS), interpreter heap overhead, or CPU and wall-clock time, and `cdid` makes no timing, scaling, or whole-process resource guarantee.
- **Separate failure modes.** Exceeding a workspace or opt-in work limit aborts the call before the stage runs without dropping units or moments. Mathematical singularity of $B$ is handled by exact reduction and `"tikhonov"` covariance regularization during point estimation. Neither a completed covariance-ridge solve nor an admitted resource profile establishes that the sample has enough units or independent clusters for accurate asymptotic normal or multiplier-bootstrap inference, or that the full joint covariance of pre-period cells is well conditioned enough to compute the unconditional pretrend Wald statistic `W`.

## Related pages

- [Options reference: GMM weighting and covariance regularization](../reference/options.md#gmm-weighting-and-covariance-regularization) and [The control list](../reference/options.md#the-control-list)
- [Results reference: The chained record](../reference/results.md#the-chained-record) and [Uncertainty](../reference/results.md#uncertainty)
- [How to compare precision](../examples/precision.md)
- [Parallel trends: PT-post and PT-all](parallel-trends.md)
- [Moving from did](../from-did.md)

---

**Reading path, step 15 of 20** · Previous: [Add PT-all restrictions](../examples/pt-all-comparison.md) · Next: [Robustness](../examples/robustness.md)

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
| 13 | [Parallel trends](parallel-trends.md) |
| 14 | [Add PT-all restrictions](../examples/pt-all-comparison.md) |
| 15 | **GMM weighting and covariance regularization** (this page) |
| 16 | [Robustness](../examples/robustness.md) |
| 17 | [Inference](../examples/inference.md) |
| 18 | [Options](../reference/options.md) |
| 19 | [Results](../reference/results.md) |
| 20 | [Diagnose errors and unavailable results](../help/troubleshooting.md) |

[Guide home](../README.md) · [Package overview](../../README.md)

</details>
