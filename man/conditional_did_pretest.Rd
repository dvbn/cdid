\name{conditional_did_pretest}
\alias{conditional_did_pretest}
\alias{print.cdid_MP.TEST}
\alias{summary.cdid_MP.TEST}
\title{Conditional pretrend test (chained by default for panels)}
\description{
\code{conditional_did_pretest()} keeps the 19 formals of
\code{did::conditional_did_pretest()} (did 2.5.1) in order and appends the
existing chained reference and selection settings it needs. For panel data the
default \code{mode = NULL} runs the chained conditional pretrend test, matching
\code{att_gt()}; \code{mode = "ordinary"} (and repeated cross sections) forwards
the upstream arguments to \code{did::conditional_did_pretest()} unchanged, with
the \code{NULL} unbalanced-panel sentinel mapped to did's \code{FALSE}. The
appended chained arguments stay lazy and are neither evaluated nor forwarded in
ordinary mode. As in did, the default \code{est_method} is \code{"ipw"},
independently of \code{att_gt()}'s \code{"dr"} default.
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
  control = NULL
)

\method{print}{cdid_MP.TEST}(x, ...)

\method{summary}{cdid_MP.TEST}(object, ...)
}
\arguments{
\item{yname, tname, idname, gname, xformla, data, panel, control_group, weightsname, clustervars, print_details}{As in
\code{did::conditional_did_pretest()} and \code{att_gt()}. In chained mode the
test coordinates are the \code{model.matrix(xformla)} columns; variables that
appear only in \code{selection_xformla} are not test coordinates.}
\item{allow_unbalanced_panel}{\code{NULL} (default) keeps every eligible unit in
chained mode (as \code{att_gt()}); \code{FALSE} is refused in chained mode. In
ordinary mode \code{NULL} means did's \code{FALSE}.}
\item{alp}{Significance level of the reported critical value.}
\item{bstrap, cband}{Validated only. As in did, the conditional test always uses
multiplier draws and reports no standard errors.}
\item{biters}{Number of multiplier draws.}
\item{est_method}{\code{"ipw"} (default), \code{"reg"} or \code{"dr"} in chained
mode. A custom function returns an explicit not-applicable result: no
conditional process is claimed for arbitrary callbacks.}
\item{pl, cores}{Multiplier draws may use \code{cores} worker processes when
\code{pl = TRUE} and \code{control$workers} declares them; signs are drawn in
the calling session, so results equal \code{cores = 1}.}
\item{anticipation, base_period, fix_weights, mode, frame, selection, selection_model, selection_pooling, selection_xformla, selection_universe, selection_rho, selection_justification, selection_est_method, imputation_pooling, control}{The
existing chained settings of \code{att_gt()} with the same meanings; they select
the fitted route whose auxiliary comparisons are tested. GMM weighting and
PT-all settings are not test arguments.}
\item{x, object}{A \code{cdid_MP.TEST} result.}
\item{...}{Unused.}
}
\details{
Null: for every declared pre-treatment cell \eqn{c} of the selected
\code{base_period} family (the pretrend test cells of \code{att_gt()}, reference zeros
excluded), \eqn{E_g[\delta_c \mid z_c] = 0}, where \eqn{z_c} is the
\code{model.matrix(xformla)} row of the unit at the cell's own reference date
(the earlier endpoint of its auxiliary comparison). This can be weaker than
conditioning on all scientific variables. The process
\eqn{J_c(u) = \sigma_c E_g[1\{z_c \le u\} \delta_c]} is the accepted auxiliary
link of the fitted union catalog with the grid indicator inserted into its
numerator integrands: IPW uses the outcome change, OR/DR the link residual and
pure imputation/augmented correction the imputed or augmented response. The
shared nuisance fits (rates, selection, propensity, imputation and the single
global outcome-regression slope) are reused unchanged; only their loadings
depend on \eqn{u}, and the complete original-unit influence function includes
all of their forcing (not \eqn{h_u} times the scalar influence function).
The links use the probability clipping of \code{att_gt()} (section
Probability clipping): capped treatment odds and floored observation
denominators enter the term weights, and the loading integrands carry the same
derivative masks. When clipping binds with positive population probability, the
tested null concerns the clipped auxiliary functional.

Grid: one row per supplied original unit (frame order; singletons, units without
observed outcomes and frame-only units included; duplicate rows keep their
mass), evaluated through the route's frozen terms, factor levels, contrasts and
transform state; grid-only rows never refit them. Dated covariates are read at
each cell's recipe date. A required grid covariate that is missing, or a level
absent from the frozen design, makes the complete test unavailable with the
unit ids (no unit is deleted and no available-only subset is tested); estimated
ATTs are unaffected. The statistic is the unstandardized
\eqn{CvM = (N/M) \sum_c \sum_m J_c(u_m)^2} with \eqn{M = N}; multiplier draws
use one Rademacher sign matrix shared by all cells and grid points, drawn from
the caller's RNG stream in the existing batch order, at the unit level or over
the full original-unit cluster alignment (clusters without influence included).
The p-value is the share of draws exceeding \eqn{CvM}; the critical value is
their type-1 quantile. A globally degenerate process (every cluster process row exactly zero)
is unavailable. Grid, indicator, process and multiplier work is streamed in
bounded chunks and charged against the existing resource caps; no
original-unit by grid influence array is formed.

Calibration requires uniform process and multiplier-bootstrap regularity
conditions, including a growing number of independent clusters for clustered
inference; scalar ATT asymptotic normality alone does not establish them, and
finite checks establish algebra, not size or power. This is not the pretrend Wald
test: the unconditional \code{W}/\code{Wpval} of \code{att_gt()} is separate and
unchanged.
}
\value{
In chained mode an object of class \code{c("cdid_MP.TEST", "MP.TEST")} with the
familiar fields \code{CvM}, \code{CvMb}, \code{CvMcval}, \code{CvMpval},
\code{clustervars} and \code{xformla}; the KS fields are \code{NULL}. The
\code{chained} element records the null, status and reason, base period, tested
cells (with their scalar values), route, grid rule and columns, statistic and
multiplier conventions, clusters, RNG provenance, execution and resource
accounting. Unavailable tests have status \code{"unavailable"} (or
\code{"not_applicable"} for custom estimators) and \code{NA} statistics. In
ordinary mode did's \code{MP.TEST} object.
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
