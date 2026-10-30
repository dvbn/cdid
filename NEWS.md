# cdid 1.1.0.0

* Construct internally estimated moment covariance matrices symmetrically, avoiding platform-dependent refusals caused by rounding in diagonal blocks.
* Preserve the original coefficient bits, including signed zero, when binding optional numerical diagnostics.
* Detect exactly duplicated and all-zero predictors before fitting treatment or observation logits, with a rank-deficiency diagnostic.

* Refuse a full pretrend Wald test when its restriction count exceeds the exact centered covariance rank bound from the number of independent unit or cluster rows. Tiny positive Cholesky pivots from rounding cannot authorize an impossible full-rank test. Point estimates and their other inference remain available.

* Avoid a rare refusal of a well-conditioned covariance when a host symmetric eigensolver loses orthogonality. One QR/Rayleigh refinement is accepted only under the existing reconstruction residual bound and is recorded in covariance diagnostics. Covariance inputs, definiteness checks and ridge settings are unchanged.


* Faster computation through batching, exact reuse of shared nuisance and regression intermediates, and bounded phase-local influence-function (IF) storage across baseline and missing-outcome (selection, imputation, and augmented) routes, without claiming that every case improved.
* Default probability safeguard using `did`/`DRDID`-compatible values (`1e-6`) with no trimming, no row exclusion, and no new public tuning knob: treatment propensities have an upper cap of `1 - 1e-6` on the `did`/`DRDID` scale, inverse-observation probabilities (and inverted positive empirical pair-rate denominators) have a lower floor of `1e-6`, and for independent date logits (`selection_model = "independent_period_logit"`) the `1e-6` lower floor applies to the final pair product.
* Raw likelihood scores and information matrices remain evaluated at raw fitted probabilities, while weight derivatives in the influence functions respect inactive clipped regions (retaining derivatives where clipping is inactive and zeroing direct weight-derivative contributions where a cap or floor binds).
* Route diagnostics report clipping activity and weight summaries, while genuine support and nuisance-estimation failures continue to be reported. Original-target and inference claims are qualified when clipping binds: clipping does not guarantee identification of the unclipped target, asymptotic robustness, or nominal coverage when active.

# cdid 1.0.0.0

This version replaces cdid 0.1.1 with a new implementation built on the `did`
package. No code from 0.1.1 is included (see `inst/PROVENANCE.md`). Scripts
written for 0.1.1 should be checked against the current help pages. Version
1.0.0.0 has not been submitted to CRAN; install it from GitHub with
`remotes::install_github("dvbn/cdid")`. The optional pinned setup described in
`validation/SETUP.md` installs exactly the tested versions of did and DRDID.

## Interface

* `att_gt()` and `aggte()` keep the upstream argument order of did 2.5.1.
  `att_gt()` adds `mode` and chained-specific arguments. Its
  `allow_unbalanced_panel = NULL` default permits incomplete histories in
  chained mode and becomes did's `FALSE` in ordinary mode; explicitly setting
  it to `FALSE` is refused in chained mode. Other upstream defaults are
  retained. `aggte()` dispatches on the result class without a `mode` argument.
  `ggdid()`, `tidy()` and `glance()` are re-exported from did, with methods for
  chained results.
* `mode = "ordinary"` passes the call unchanged to did. It is also the default
  for repeated cross sections. An explicit `allow_unbalanced_panel = NULL` is
  translated to did's `FALSE`.
* Chained estimation is the default for panel data. It uses the outcome changes
  between the dates at which each unit is observed and combines them by GMM.
  Units with incomplete histories are retained, and effects the comparisons
  cannot recover are reported as unavailable with a reason.
* `frame` supplies the unit table of the study population, including units with
  one or no observed outcome.

## Estimation

* Treatment adjustment: chained doubly robust estimation (`est_method = "dr"`,
  the default), chained outcome regression (`"reg"`) and chained IPW (`"ipw"`).
  OR and DR fit one pooled regression of the controls' outcome changes with
  local block intercepts. Custom did-style IPW callbacks are accepted and may
  return the optional `log_weight_gradient`.
* Missing outcomes: `selection = "conditional_mtar"` (legacy alias
  `"corollary1"`) corrects for missingness that depends on covariates.
  `selection_est_method` chooses observation weighting (`"ipw"`, default),
  change imputation (`"imputation"`) or the augmented correction (`"dr"`, with
  `est_method = "dr"` only). Related arguments: `selection_model`
  (`"pair_logit"` or `"independent_period_logit"`), `selection_pooling`,
  `imputation_pooling`, `selection_xformla`, `selection_universe`,
  `selection_rho` and `selection_justification`.
* Parallel trends: `pt_assumption = "post"` (default) or `"all"` with an
  explicit common `pt_start`.
* GMM: `gmm_weighting = "optimal"` (alias `"two_step"`, default) or
  `"identity"`, with Tikhonov handling of a singular or near-singular moment
  covariance (`gmm_covariance = "tikhonov"`).
* Chained `fix_weights` accepts `NULL`, `"first_period"` and `"base_period"`;
  covariates and sampling weights that change within a unit are supported with
  an explicit reference date.
* `conditional_did_pretest()` runs a chained conditional pretrend test for
  panels by default and delegates to did in ordinary mode.
* Results keep did's `MP` and `AGGTEobj` fields, with complete influence
  functions, analytic or multiplier inference, clustering on at most one
  variable beyond `idname`, and a `chained` record of support, models and
  settings. `aggte()` aggregates chained results without refitting.
* `print()` and `summary()` of chained results show did's table under a
  heading naming the chained estimator, followed by a short chained record:
  available cells and the reasons for unavailable ones, the unit table, the
  correction for missing outcomes, parallel trends, GMM weighting and
  regularization, the inference performed and the pretrend test. Aggregates
  also list restricted or unavailable targets and the bins excluded from the
  overall summary. Each paper reference is shown once in the heading; detailed
  technical records stay in the result.
* `compute_inffunc = FALSE`, `print_details`, and `pl`/`cores` (with workers
  declared in `control$workers`) follow did's meaning. `control` sets numeric
  workspace, workers, chunk rows, an optional estimated-work cap and the
  covariance regularization constants.

## Qualifications

* The package requires did 2.5.1 or later and DRDID 1.3.0 or later, and
  refuses to load with earlier versions. It has been tested with did 2.5.1 and
  DRDID 1.3.0; later versions are accepted, but it is not guaranteed that
  everything works with them or that results are unchanged.
* Chained OR and DR, change imputation, the augmented correction and PT-all
  extend the published chained DiD paper. Each method requires the assumptions
  stated in its help page and in the guide.

# cdid 0.1.1

* Previously released version ([deprecated repository](https://github.com/joelcuerrier/cdid)), published on CRAN.
