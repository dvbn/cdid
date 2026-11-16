\name{att_gt}
\alias{att_gt}
\title{Difference-in-differences estimates for panel data}
\description{
Estimate average treatment effects by cohort and date or, for a binary
treatment that can turn on and off, by switching date and direction. For
panels, the default method combines changes measured within units across
the dates at which they are observed. Ordinary mode uses \pkg{did}.
}
\usage{
att_gt(
  yname,
  tname,
  idname = NULL,
  gname,
  xformla = NULL,
  data,
  panel = TRUE,
  allow_unbalanced_panel = NULL,
  control_group = c("nevertreated", "notyettreated"),
  anticipation = 0,
  weightsname = NULL,
  fix_weights = NULL,
  alp = 0.05,
  bstrap = TRUE,
  cband = TRUE,
  biters = 1000,
  clustervars = NULL,
  est_method = "dr",
  base_period = "varying",
  faster_mode = TRUE,
  print_details = FALSE,
  pl = FALSE,
  cores = 1,
  compute_inffunc = TRUE,
  ...,
  mode = NULL,
  frame = NULL,
  selection = "none",
  selection_model = "pair_logit",
  selection_pooling = "separate",
  selection_xformla = NULL,
  selection_universe = "F",
  selection_rho = NULL,
  selection_justification = NULL,
  selection_est_method = "ipw",
  imputation_pooling = "separate",
  gmm_weighting = c("optimal", "identity", "within"),
  gmm_covariance = "tikhonov",
  pt_assumption = c("post", "all", "ctrl"),
  pt_start = NULL,
  cross_pair = FALSE,
  control = NULL,
  outcome_pooling = "separate",
  dname = NULL,
  imputation_scale = "per_period",
  switchers = "both",
  switch_times = NULL,
  missing_reference = "exclude"
)
}
\arguments{
\item{yname}{Name of the numeric outcome column. Use \code{NA} for missing outcomes.}
\item{tname}{Name of the numeric date column, on the same scale as the treatment dates.}
\item{idname}{Name of the unit identifier column. Required in chained mode.}
\item{gname}{Name of the first-treatment-date column. Use zero for never-treated units. Supply either \code{gname} or \code{dname}.}
\item{xformla}{One-sided covariate formula, such as \code{~ age + income}.
  The default \code{NULL} uses an intercept only.}
\item{data}{A data frame in long format. A panel may have at most one row per unit and date.}
\item{panel}{Whether the data follow the same units over time.
  Use \code{FALSE} for repeated cross sections.}
\item{allow_unbalanced_panel}{In ordinary mode, \code{NULL} uses \pkg{did}'s
  default of balancing the panel. Chained mode uses incomplete panels
  with \code{NULL} or \code{TRUE}. To require a balanced sample, use
  \code{mode = "ordinary"} and \code{allow_unbalanced_panel = FALSE}.}
\item{control_group}{For cohort models, use never-treated units (\code{"nevertreated"}) or
  units not yet treated at the comparison date (\code{"notyettreated"}).}
\item{anticipation}{Number of periods before adoption in which treatment
  may affect outcomes. The reference date is the last date before that window.}
\item{weightsname}{Name of a column of positive sampling weights for cohort effects.}
\item{fix_weights}{For weights that vary over time, use \code{"first_period"}
  or \code{"base_period"} in chained mode. These use the first study date
  or the cohort's reference date, respectively. The default is suitable
  for constant weights.}
\item{alp}{Significance level for confidence intervals and bands.}
\item{bstrap}{Whether to use the multiplier bootstrap for inference.}
\item{cband}{Whether to report simultaneous confidence bands when \code{bstrap = TRUE}.}
\item{biters}{Number of bootstrap draws.}
\item{clustervars}{Names of cluster variables. Chained mode accepts at most
  one cluster variable in addition to the unit identifier. Clustering affects
  inference.}
\item{est_method}{Treatment adjustment. Use \code{"dr"} for doubly robust
  estimation, \code{"reg"} for outcome regression, or \code{"ipw"} for inverse
  probability weighting. Chained cohort estimation also accepts a custom
  IPW function, described below.}
\item{base_period}{Use \code{"varying"} for adjacent-date comparisons before
  treatment or \code{"universal"} for comparisons to the cohort reference
  date. Both use the same reference date for post-treatment effects.
  Leave at its default for switching models.}
\item{faster_mode}{Passed to \pkg{did} in ordinary mode.}
\item{print_details}{Whether to print information about the progress of the fit.}
\item{pl, cores}{Use \code{pl = TRUE} and set \code{cores} to run bootstrap
  calculations on multiple workers. In chained mode, also set
  \code{control = list(workers = cores)}.}
\item{compute_inffunc}{Set to \code{FALSE} for point estimates only. Use
  \code{TRUE} to obtain standard errors, pretrend tests and
  \code{\link{aggte}()} summaries.}
\item{...}{Additional arguments to a custom estimation function.}
\item{mode}{Use \code{"chained"} for within-unit comparisons or
  \code{"ordinary"} to call \code{did::att_gt()}. By default, panels use
  chained mode and repeated cross sections use ordinary mode.}
\item{frame}{For cohort models, an optional study-unit table with one row per unit, matched to
  \code{data} by \code{idname}. Include units with one or no observed outcomes
  if they belong to the study population. Without this table, the units in
  \code{data} define the population.}
\item{selection}{Use \code{"none"} when units observed at both dates represent
  their cohort's mean change. Use \code{"conditional_mtar"} when this
  representativeness holds conditional on covariates. In cohort models,
  \code{"corollary1"} remains an alias for the latter with direct pair logits.}
\item{selection_model}{Observation-probability model. \code{"pair_logit"} fits
  the probability of observing both dates directly. \code{"independent_period_logit"}
  multiplies two date-specific probabilities and assumes conditional independence
  of observation at the two dates. This argument applies to observation
  weighting and augmentation.}
\item{selection_pooling}{How observation models share slopes. For cohort
  models use \code{"separate"}, \code{"within_cohort"} or \code{"across_cohorts"}.
  For switching models use \code{"separate"}, \code{"within_group"} or
  \code{"pooled"}. Within-group sharing keeps switch date, direction and
  switcher or stayer status fixed. Each pair or date has its own intercept
  when the formula includes one.}
\item{selection_xformla}{One-sided formula for observation probabilities.
  The default uses \code{xformla}. Imputation uses \code{xformla}.}
\item{selection_universe}{Sample for the observation model. \code{"F"} uses
  all study units. Direct pair weighting also permits \code{"B"} (at least
  one outcome) or \code{"A"} (at least two). The target remains the full study
  population.}
\item{selection_rho}{For \code{"A"} or \code{"B"}, an optional column name containing the known, normalized density ratio
  of covariates in the fitting sample to those in the study population.
  Supply it in the unit table for cohort models or in \code{data} for
  switching models. The ratio is defined within cohort or within the
  corresponding switcher and stayer groups.
  Values must be positive. \code{NULL} uses a ratio of one.}
\item{selection_justification}{Required explanation when fitting on \code{"A"}
  or \code{"B"}. Explain why the fitting sample is representative within the cohort or switcher and stayer group, or how
  the known density ratios were obtained.}
\item{selection_est_method}{With \code{selection = "conditional_mtar"}, choose
  \code{"ipw"} for observation weighting, \code{"imputation"} to predict actual
  changes, or \code{"dr"} to combine these approaches. The last two require
  a built-in treatment method, the full study population (\code{"F"}), and
  complete required covariates for that population. Observation weighting
  and augmentation need an estimable probability model and overlap.}
\item{imputation_pooling}{How regressions of actual changes share slopes.
  Uses the same choices as \code{selection_pooling}, independently of it.}
\item{gmm_weighting}{How comparisons are combined. \code{"optimal"} (also
  \code{"two_step"}) estimates weights using the covariance of the comparisons.
  \code{"identity"} uses equal weights. \code{"within"} uses the covariance of
  within-unit changes implied by uncorrelated level shocks and requires a
  built-in cohort method.}
\item{gmm_covariance}{Use \code{"tikhonov"} to regularize singular or nearly
  singular comparison covariance matrices when needed.}
\item{pt_assumption}{Parallel trends assumption for cohorts. \code{"post"}
  starts at each cohort's reference date. \code{"ctrl"} adds earlier common
  trends among eligible controls. \code{"all"} also includes the target
  cohorts in those earlier restrictions. This is separate from \code{base_period}.}
\item{pt_start}{Common calendar date from which to apply the extra restrictions.
  Required for \code{"all"}, optional for \code{"ctrl"} (defaulting to the
  first study date), and unused with \code{"post"}.}
\item{cross_pair}{With a built-in cohort method, set to \code{TRUE} to use target
  and control changes measured over different date pairs when the assumed
  trends link those changes.}
\item{control}{Named list of advanced settings. Use
  \code{list(keep_nuisance = TRUE)} to retain fitted model coefficients and
  standard errors, or \code{list(preset = "cohort_time_fe")} for the linear
  cohort-by-date fixed-effects specification described below. Other options
  include \code{numeric_workspace_bytes} (default one GiB for numerical
  arrays), \code{workers} (default one), \code{chunk_rows},
  \code{counted_work_limit} (default \code{Inf}), \code{condition_trigger}
  (default \code{1e-10}) and \code{ridge_relative} (default \code{1e-4}).}
\item{outcome_pooling}{Use \code{"separate"} for comparison-specific outcome
  slopes or \code{"pooled"} for shared slopes with local intercepts when
  the formula includes an intercept. For cohort models, pooled slopes scale
  with the number of calendar steps in a comparison. For switching models,
  a named character vector or list can assign transition cells to groups
  that share slopes. Applies to \code{est_method = "reg"} and \code{"dr"}.}
\item{dname}{Name of the period-specific treatment indicator, coded 0 or 1.
  With this argument, effects are estimated at entry or exit under a static
  treatment model with no carryover. Supply a known treatment status for
  every unit and study date, including dates with missing outcomes.}
\item{imputation_scale}{For change imputation in a switching model,
  \code{"per_period"} scales covariate slopes by the number of calendar
  steps in the change. \code{"per_change"} leaves slopes unscaled.}
\item{switchers}{In a switching model, use \code{"entry"} for changes from 0 to
  1, \code{"exit"} for changes from 1 to 0, or \code{"both"}.}
\item{switch_times}{In a switching model, dates at which to estimate effects.
  The default uses dates with switchers in the requested directions.
  A requested date and direction with no switchers is reported as unavailable.}
\item{missing_reference}{For chained cohort models with time-varying
  covariates, how to handle a missing value required at the reference date.
  \code{"exclude"} (the default) warns and limits the affected cohort's
  comparisons to units with the required reference covariates, including
  those in the observation formula when it is used.
  \code{"error"} stops and identifies the missing inputs. Ordinary mode
  does not use this argument.}
}
\details{
Each observed change uses two outcomes from the same unit. Units with fewer
than two outcomes can still belong to the study population. Keep them in
\code{data} or supply their information in \code{frame}. Cohort and population
weights then describe the study units you specified. Chaining can connect
different sets of units that supply changes over overlapping dates, subject
to the comparison-group support and parallel trends assumptions.

For time-varying covariates, post-treatment cohort comparisons use covariates
at the cohort's reference date. Comparisons reported before treatment use
covariates from the earlier comparison date. Keep the rows providing these
values even if their outcomes are missing. When sampling weights vary, supply
the required weights at the chosen reference date. For \code{\link{aggte}()},
cohort shares use first-date sampling weights.

When a required reference covariate is missing, the default restricts the
affected cohort's comparisons to units with the required covariates at its
reference date. A unit can still contribute to another cohort whose reference
covariates it has. Supply independently known values on the appropriate
unit-date row, using \code{NA} for a missing outcome. Covariates that are
constant within a unit are read as usual. If treated units are excluded, the
cohort effect describes the remaining treated units, and \code{\link{aggte}()}
uses those units to calculate the cohort's first-date weight. Influence
functions retain a row for every original unit. Effects without sufficient
support and pretrend comparisons with missing required inputs are reported as
unavailable. See \code{fit$chained$reference_eligibility} for the affected
cohorts, dates and counts. Missing sampling weights at a required reference
date still cause an error. If the observation model is fitted on a subset,
its transport assumptions and any supplied density ratios apply within each
cohort's eligible sample.

The treatment method and the choice of missing-outcome adjustment are separate.
Outcome regression for treatment predicts the untreated change. Change
imputation predicts the actual change within its group. With augmentation,
the package also uses observation probabilities to reweight the residuals of
that imputation. With treatment \code{"dr"} and observation \code{"dr"}, either
model in each pair can be misspecified if the other is correct and the
remaining identifying conditions hold. Extreme probability weights are
bounded. Check overlap and clipping before interpreting a fit with large
weights, since clipping can affect the target and its robustness properties.

A switching model compares entrants with units staying untreated and exiters
with units staying treated. It assumes that past treatment has no effect
beyond current treatment status. Exit effects are reported in the
treated-minus-untreated direction. This route supports all three built-in
treatment methods and missing-outcome adjustments. The cohort-specific trend
and cross-pair options apply to \code{gname} models.

The linear fixed-effects specification uses
\code{control = list(preset = "cohort_time_fe")}, consecutive positive integer
dates, unit weights, zero anticipation, \code{base_period = "universal"} and a
never-treated cohort. Use \code{pt_assumption = "post"} for never-treated
controls, or \code{pt_assumption = "all"} with \code{pt_start} at the first
date for not-yet-treated controls. Any covariates must be time-invariant and
available for all study units. Regression uses within-unit least squares.
IPW and DR use the common conditional untreated path evaluated in each cohort.
Leave \code{outcome_pooling} and \code{cross_pair} at their defaults.

A custom chained cohort estimator follows the panel IPW interface
\code{f(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...)}. Return a list
with a scalar \code{ATT} and \code{att.inf.func} in the supplied row order.
For fitted observation weights, also return \code{log_weight_gradient}, the
derivative of the estimate with respect to each log weight, including refits
inside the function. See the
\href{https://github.com/dvbn/cdid/blob/main/docs/methods/treatment.md#custom-ipw-functions}{treatment adjustment guide}
for the callback requirements.
}
\value{
Ordinary mode returns the \code{MP} object from \pkg{did}. Chained
cohort estimation returns a \code{cdid_MP} object with cohort-date estimates,
uncertainty and comparison support. The joint pretrend test is reported when
the full set of requested comparisons and its covariance are usable.

Switching models return a \code{cdid_static_att_gt} object with effects by date
and direction and a switcher-weighted estimate in \code{$aggregate}.
Use \code{summary()}, \code{\link{tidy}()} or \code{\link{ggdid}()} to inspect
either type of result. Use \code{\link{aggte}()} to aggregate cohort effects.
}
\seealso{
\code{\link{aggte}()}, \code{\link{chain_report}()},
\code{\link{nuisance_models}()}, \code{\link{conditional_did_pretest}()},
\code{vignette("cdid", package = "cdid")},
\href{https://github.com/dvbn/cdid/tree/main/docs}{User guide}
}
\references{
Bellego, C., Benatia, D. and Dortet-Bernadet, V. (2025). The Chained
Difference-in-Differences. \emph{Journal of Econometrics} 248, 105783.
\doi{10.1016/j.jeconom.2024.105783}

Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with
Multiple Time Periods. \emph{Journal of Econometrics} 225(2), 200--230.
\doi{10.1016/j.jeconom.2020.12.001}

de Chaisemartin, C., D'Haultfoeuille, X., Pasquier, F., Sow, D. and
Vazquez-Bare, G. (2026). Difference-in-Differences Estimators for Treatments
Continuously Distributed at Every Period. Working paper, arXiv:2201.06898v7
for the static binary treatment model.
}
\examples{
data("mpdta", package = "did")
fit <- att_gt(
  yname = "lemp", tname = "year", idname = "countyreal",
  gname = "first.treat", data = mpdta, bstrap = FALSE, cband = FALSE
)
summary(fit)
aggte(fit, type = "dynamic")
}
