# Cohort-by-date fixed effects

Use this preset for an absorbing treatment and a common linear untreated time path with unit fixed effects. Each treated cohort can have a different effect at each post-treatment date. The fit uses changes across each unit's observed dates, including changes that span gaps.

## Fit the model

Supply a long panel with one row per unit and date, an outcome, and a cohort column. Use the adoption date as the cohort code and zero for never-treated units. Include at least one never-treated cohort, use consecutive positive integer study dates, and place every adoption after the first study date. You can leave rows with missing outcomes in the panel. An optional unit table supplies units with no observed outcomes and keeps the original influence-function rows.

```r
fit <- cdid::att_gt(
  yname = "y", tname = "t", idname = "id", gname = "g",
  data = panel, est_method = "reg", base_period = "universal",
  control_group = "nevertreated",
  bstrap = FALSE, cband = FALSE,
  control = list(preset = "cohort_time_fe")
)
cdid::tidy(fit)
```

With never-treated controls, the pre-treatment period immediately preceding adoption is the reference. Other pre-treatment cells have their own estimated deviations. A joint test uses the covariance of the full declared set of pre-treatment deviations when it is available.

For the common untreated path across both never-treated and not-yet-treated cohorts, use `control_group = "notyettreated"`, `pt_assumption = "all"` and `pt_start = min(panel$t)`. That model imposes the common path on all untreated cells from the first study date. It reports post-treatment effects. Its pre-treatment restrictions are part of the fitted model.

Use zero anticipation and unit sampling weights. The example uses `selection = "none"`. The usual `clustervars` and multiplier inference options apply. The fit reports support for each requested cell. Inspect the [missing-outcomes guide](missing-outcomes.md) when observed changes may be selective. The cohort-by-date model needs the observed changes to satisfy its within-unit mean restriction. With this setting the preset fits the model to observed changes directly.

## Adjust for covariates

Set `xformla = ~ x1 + x2` to allow the untreated calendar path and each
cohort-date effect to vary with time-invariant covariates. Use
`outcome_pooling = "separate"`. The fitted design has one set of common
calendar slopes and separate cohort-date deviations.

Supply the covariates for every unit in the study frame, including units
observed once or never. The reported treatment effect averages the fitted
contrast over the corresponding cohort in that frame. A unit without an
observed outcome change contributes to that covariate average and its
influence function. The solve uses connected components with requested
contrasts and removes one coefficient origin in each of them. Their remaining
covariate design must have sufficient rank. The fit reports rank deficiency
or numerical singularity in that system. Disconnected cells remain unavailable.

## Treatment weighting and doubly robust adjustment

Set `est_method = "ipw"` to reweight clean comparison units to each
treated cohort's covariate distribution. Set `est_method = "dr"` to combine
a shared untreated change model with reweighted residuals. Both methods
estimate the average untreated change separately for each cohort, over
its entire study frame. If the same conditional trend is steeper at high
values of a covariate, cohorts with more high-covariate units can have
larger untreated changes.

Propensity models use the target and donor cohorts in the study frame
and include an intercept. They are reused across eligible pairs of
dates. A never-treated fit uses only the never-treated donor. A
not-yet-treated fit can also use each other cohort for changes ending
before its adoption. The DR outcome nuisance fits a shared linear
calendar path to clean changes. The default `gmm_weighting = "optimal"`
uses the covariance of the original-unit moment scores. Choose
`gmm_weighting = "identity"` for equal weight on each retained moment.
The exact R and Stata regression matches below use `est_method = "reg"`.

These weighting and DR fits require support for the treated-cohort
average and overlap with the clean donors used to estimate it. With
`selection = "none"`, a sufficient condition is that the observations
behind each cohort and date pair preserve its covariate distribution
and conditional mean change. The observed actual changes of the
treated cohort must represent its full-cohort changes. Given those
conditions and the stated parallel trends, DR is consistent when
either the treatment propensity or the untreated change model is
correct. The propensity-only guarantee requires the weighting odds
to remain correct after the package's probability cap. See the
[missing-outcomes guide](missing-outcomes.md).

## Adjust for selective observation

Use `selection = "conditional_mtar"` when observed within-unit changes are
representative of cohort changes conditional on the covariates in the relevant
selection or imputation model. Choose the treatment method with `est_method`
and the observation correction with `selection_est_method`.

| Observation correction | Response for each unit and supported date pair |
| --- | --- |
| `"ipw"` | The observed change divided by its probability of being observed, or zero when either outcome is missing. |
| `"imputation"` | The fitted conditional mean change for that cohort and pair. |
| `"dr"` | The fitted conditional mean change plus the observed residual divided by its observation probability. A unit with a missing change contributes the fitted mean. |

With `"pair_logit"`, the observation model predicts whether both outcomes
are observed. With `"independent_period_logit"`, it predicts observation at
each date and multiplies the probabilities. The latter requires conditional
independence of the two observation indicators. Weighting and augmentation
use `selection_xformla`, which can differ from the treatment and imputation
formula `xformla`. The independent-period model and both imputation routes
use the full study frame (`selection_universe = "F"`). Direct-pair weighting
can also use `"A"` or `"B"` with a justified transport ratio. The full-frame
probability is then the conditional fitted probability times the known ratio
and the share of the cohort in the fitting sample. See the
[missing-outcomes guide](missing-outcomes.md) for these assumptions and options.

For example, augment the observation correction and use DR for treatment:

```r
fit_dr <- cdid::att_gt(
  yname = "y", tname = "t", idname = "id", gname = "g",
  data = panel, xformla = ~ x1, est_method = "dr",
  selection = "conditional_mtar", selection_est_method = "dr",
  selection_model = "pair_logit", selection_xformla = ~ x1 + x2,
  base_period = "universal", control_group = "nevertreated",
  bstrap = FALSE, cband = FALSE, control = list(preset = "cohort_time_fe")
)
```

The adjusted responses enter the original cohort and calendar restrictions.
The regression combines them with the within-unit least-squares geometry.
Treatment IPW and DR average the adjusted responses over the full cohort
frame, including missing-change rows. In treatment DR, the calendar model
uses the adjusted clean changes and treatment propensities use the target
and donor cohorts. The inference includes estimating the active observation,
imputation, outcome and propensity models. A missing adjustment or weak
overlap can still invalidate the effect. A date pair with no observed change
in a cohort supplies no fitted comparison there. Corrected fits use the
connected comparisons that identify a requested contrast. Treatment DR
also uses the clean pairs needed to fit its common calendar changes.
All supported links in those components enter the fit.

## Match R etwfe and Stata jwdid

The following commands match the cell specification in the example above. In R, explicitly supply the unit identifier and the never-treated cohort reference.

```r
ref <- etwfe::etwfe(
  y ~ 1, tvar = t, gvar = g, ivar = id, gref = 0,
  data = panel, cgroup = "never", vcov = ~id,
  ssc = fixest::ssc(K.adj = FALSE, G.adj = FALSE)
)
```

```stata
jwdid y, ivar(id) tvar(t) gvar(g) never cluster(id)
```

For the not-yet-treated comparison, set `cgroup = "notyet"` in R and omit `never` in Stata. Match the same input rows, calendar, controls and cohort-date coefficients. On incomplete panels, verify that each command uses the same observations. For example, `jwdid` can exclude treated units whose first observed outcome is already treated. Cohort and date labels still refer to the same adoption and study periods.

With `bstrap = FALSE`, cdid uses the original-unit sandwich without a degrees-of-freedom multiplier. Pass an equivalent cluster to all commands. The R call above requests that normalization explicitly. When Stata fits with `reghdfe`, divide the stored `e(V)` by the finite-sample factor
`(e(N) - 1) / (e(N) - e(df_m) - e(df_a) - (e(df_a_nested) > 0)) * e(N_clust) / (e(N_clust) - 1)`
to compare covariance on the same scale. The corresponding standard-error conversion divides by the square root of the factor. Check critical values separately when comparing intervals.

These comparisons concern matched cohort-date effects and their joint covariance. On incomplete panels, `cdid::aggte()` weights cohorts using the study-unit frame. The reference post-estimation routines can weight the observed treated cells. Compare those summaries only after aligning the weights and whether the weights contribute to sampling uncertainty.

[etwfe model reference](https://grantmcdermott.com/etwfe/reference/etwfe.html) · [jwdid model reference](https://friosavila.github.io/app_metrics/app_metrics11.html)

## How the changes are combined

For the regression fit, eliminating a unit intercept gives the same least-squares objective as forming consecutive differences along that unit's observed dates and using the inverse of their joint covariance under uncorrelated level errors with equal variance. For example, a unit observed in dates 1, 3 and 4 contributes changes from 1 to 3 and 3 to 4. The changes share the outcome at date 3. The fit accounts for that dependence when combining them. Inference then uses the empirical unit or cluster scores.

The regression solver reuses parts of the within-unit calculation for units in the same cohort with the same observation dates. One factorization supplies all available cohort-date estimates and their scores. The observed changes must connect each cell to the corresponding calendar path. Covariate adjustment also requires enough variation in the relevant covariates to identify the fitted contrast. Treatment IPW and DR instead combine observed changes and cohort-standardized
untreated changes by GMM. Each link enters once for each supported comparison
cohort, and the default weighting accounts for the covariance across those
moments. See the [derivation](../theory/cohort_time_fe_linear.tex).

## Verified comparisons

### Intercept-only model

The intercept-only fits were compared on the same complete and incomplete synthetic panels with 128, 2,048 and 10,000 units and 6, 12 and 20 dates, respectively. The 24 cases cover both control groups and unit or larger cluster levels. The incomplete design has non-nested observation dates. Every treated unit has a pre-treatment observation in these comparisons.

Across the 24 cases, the largest absolute ATT difference is $1.9 \times 10^{-10}$ against R etwfe 0.6.2 and $3.3 \times 10^{-12}$ against Stata jwdid. With the matched sandwich convention, the largest absolute standard-error difference is $2.8 \times 10^{-8}$ against etwfe at its default fixed-effect tolerance. It falls to $3.1 \times 10^{-12}$ when etwfe uses `fixef.tol = 3e-12`. After the documented variance conversion, the Stata standard-error difference is below $7.0 \times 10^{-13}$.

These are median seconds over three runs at 10,000 units and 20 dates with one computing thread and clustering by unit. Timed calls include the regression and reported coefficient covariance. The R benchmark uses the default fixed-effect tolerance. Stata uses `tolerance(1e-12)`. Peak resident memory measures the program process, including package startup and reading the input.

| Panel | Controls | cdid seconds | R etwfe seconds | Stata jwdid seconds | cdid MiB | R etwfe MiB | Stata jwdid MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Complete | Never treated | 0.049 | 1.516 | 5.441 | 501 | 1,100 | 1,162 |
| Complete | Not yet treated | 0.037 | 1.245 | 2.739 | 478 | 958 | 829 |
| Incomplete | Never treated | 0.045 | 1.084 | 4.756 | 533 | 730 | 899 |
| Incomplete | Not yet treated | 0.036 | 0.842 | 2.385 | 482 | 653 | 553 |

The largest panel has five treated cohorts. For panels with many more distinct cohort-date cells, the current solver's dense factorization can become the largest cost. Recovering standard errors also requires retaining the requested unit-by-effect influence functions.

### Time-invariant covariates

For the covariate comparison use `y ~ x` in R and `jwdid y x, xasis ...`
in Stata. The coefficient on a cohort-date indicator can refer to different
covariate values in the two programs. Evaluate the fitted intercept and
slope contrast at the same cohort covariate mean in each program. To match
a covariance that holds that mean fixed, compare the coefficient score
at that common evaluation point. The default cdid covariance also includes
sampling variation in the estimated cohort mean.

In nine matched cases with one covariate, the largest absolute ATT
differences were $2.9 \times 10^{-12}$ against R and $2.7 \times 10^{-13}$
against Stata. The largest differences in the fixed-evaluation joint
covariances were $7.8 \times 10^{-14}$ and $2.7 \times 10^{-13}$,
respectively. At 10,000 units and 20 dates, the incomplete not-yet-treated
fit took 0.035 seconds in cdid, 3.941 in R etwfe and 6.198 in Stata jwdid
using the matching fixed-effect tolerances and covariance conventions.
The complete and incomplete 5,000-unit, 10-date fits with four covariates
took 0.017 to 0.026 seconds including dynamic aggregation, with peak
process memory between 440 and 451 MiB on this machine.

## Match the simple summary

To compare with `etwfe::emfx(ref, type = "simple", compress = FALSE)` or Stata's `estat simple`, weight cohort-date effects by the number of observed treated rows in those cells. The following computes the same linear contrast and holds the observed cell weights fixed for inference.

```r
post <- !is.na(panel$y) & panel$g > 0 & panel$t >= panel$g
cell <- paste(panel$g[post], panel$t[post], sep = ":")
fit_cell <- paste(fit$group, fit$t, sep = ":")
position <- match(cell, fit_cell)
stopifnot(length(position) > 0L, !anyNA(position))
w <- tabulate(position, nbins = length(fit_cell)) / length(position)
keep <- which(w > 0)
stopifnot(all(is.finite(fit$att[keep])))
V <- fit$V_analytical[keep, keep, drop = FALSE] / fit$n
c(estimate = sum(w[keep] * fit$att[keep]),
  std.error = sqrt(drop(crossprod(w[keep], V %*% w[keep]))))
```

On the 128-unit validation panels, these simple summaries agree across both control pools and both observation designs. The point estimates agree to numerical precision. The maximum absolute standard-error gap is $1.6 \times 10^{-8}$ against `emfx()` and below $10^{-12}$ against Stata after converting its covariance to the same normalization.
