\name{conditional_did_pretest}
\alias{conditional_did_pretest}
\alias{print.cdid_MP.TEST}
\alias{summary.cdid_MP.TEST}
\title{Test for parallel trends before treatment}
\description{
Test whether pre-treatment differences in trends have zero mean conditional
on the supplied covariates. Panels use the chained test by default.
Use \code{mode = "ordinary"} or repeated cross sections to use the test
from \pkg{did}. The default treatment adjustment is \code{"ipw"}.
}
\usage{
conditional_did_pretest(
  yname,
  tname,
  idname = NULL,
  gname,
  xformla = NULL,
  data,
  panel = TRUE,
  allow_unbalanced_panel = NULL,
  control_group = c("nevertreated", "notyettreated"),
  weightsname = NULL,
  alp = 0.05,
  bstrap = TRUE,
  cband = TRUE,
  biters = 1000,
  clustervars = NULL,
  est_method = "ipw",
  print_details = FALSE,
  pl = FALSE,
  cores = 1,
  anticipation = 0,
  base_period = "varying",
  fix_weights = NULL,
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
  control = NULL,
  outcome_pooling = "separate"
)

\method{print}{cdid_MP.TEST}(x, ...)

\method{summary}{cdid_MP.TEST}(object, ...)
}
\arguments{
\item{yname, tname, idname, gname, xformla, data, panel, control_group, weightsname, clustervars, print_details}{As in
  \code{did::conditional_did_pretest()} and \code{att_gt()}. In chained
  mode the test coordinates are the columns of \code{model.matrix(xformla)}.
  Variables that appear only in \code{selection_xformla} are not test
  coordinates.}

\item{allow_unbalanced_panel}{\code{NULL} or \code{TRUE} allows incomplete panels
  in chained mode. Use \code{mode = "ordinary"} to balance the panel.
  In ordinary mode, \code{NULL} uses \pkg{did}'s default of \code{FALSE}.}

\item{alp}{Significance level of the reported critical value.}

\item{bstrap, cband}{Validated only. As in did, the conditional test always
  uses multiplier draws and reports no standard errors.}

\item{biters}{Number of multiplier draws.}

\item{est_method}{Treatment adjustment of the tested comparisons:
  \code{"ipw"} (default), \code{"reg"} or \code{"dr"} in chained mode. A
  custom function returns a result with status \code{"not_applicable"},
  because no conditional process is defined for arbitrary functions.}

\item{pl, cores}{Multiplier draws may use \code{cores} worker processes when
  \code{pl = TRUE} and \code{control$workers} declares them. The calling
  session draws the same signs in the same order for every worker count.
  On Windows, parallel evaluation can change the last digits of the
  computed statistics through floating-point rounding.}

\item{anticipation, base_period, fix_weights, mode, frame, selection, selection_model, selection_pooling, selection_xformla, selection_universe, selection_rho, selection_justification, selection_est_method, imputation_pooling, outcome_pooling, control}{The
  chained settings of \code{att_gt()}, with the same meanings. They select
  the fitted route whose pre-treatment comparisons are tested. The
  \code{selection} arguments describe the adjustment for missing outcomes
  and are separate from the treatment adjustment \code{est_method}. GMM
  weighting, cross-pair chaining and the \code{pt_assumption} settings are
  not test arguments.}

\item{x, object}{A \code{cdid_MP.TEST} result.}

\item{...}{Unused.}
}
\details{
The null is that every requested pre-treatment comparison has zero mean
conditional on the covariate design from \code{xformla}. The comparison
dates follow \code{base_period}. Reference zeros are excluded. The test uses
outcome changes for IPW and outcome residuals for regression or doubly
robust adjustment. When a missing-outcome adjustment is requested, the
comparison uses the corresponding adjusted changes.

The test evaluates the comparisons over the covariate values of all
supplied study units, including units with one or no observed outcomes.
Changing covariates are read at the reference date of the comparison.
If a required covariate is missing or a prediction is unsupported, the
test reports why it is unavailable.

The Cramer-von Mises statistic combines these comparisons. Its p-value
and critical value come from multiplier draws over units or the requested
clusters. Use \code{set.seed()} to reproduce the draws. Clustered testing
requires many independent clusters. The unconditional Wald test is
reported separately by \code{\link{att_gt}()}.
}
\value{
In chained mode, an object of class \code{c("cdid_MP.TEST", "MP.TEST")}
with the did fields \code{CvM}, \code{CvMb}, \code{CvMcval}, \code{CvMpval},
\code{clustervars} and \code{xformla}. The Kolmogorov-Smirnov fields are
\code{NULL}. The \code{chained} element records the null, status and reason,
base period, tested cells with their scalar values, route, grid rule and
columns, statistic and multiplier conventions, clusters, random-number
record, execution and resource accounting. Unavailable tests have status
\code{"unavailable"} (or \code{"not_applicable"} for custom functions) and
\code{NA} statistics. In ordinary mode, did's \code{MP.TEST} object.
}
\seealso{\code{\link{att_gt}}, \code{\link{aggte}}, \code{\link{reexports}}, \code{vignette("cdid", package = "cdid")}}
\examples{
\dontrun{
pre <- cdid::conditional_did_pretest(yname = "y", tname = "period", idname = "id",
  gname = "G", xformla = ~x, data = panel_data, control_group = "notyettreated",
  biters = 199)
summary(pre)
pre$chained$cells
}
}
