\name{att_gt}
\alias{att_gt}
\title{Group-time average treatment effects: chained or ordinary mode}
\description{
\code{att_gt()} has the interface of \code{did::att_gt()} (did 2.5.1) plus a
\code{mode} argument, the arguments that configure chained estimation and one
chained-only advanced \code{control} list.

\code{mode = "ordinary"} passes the call unchanged to \code{did::att_gt()}
and returns its \code{MP} object. For panel data, the default
\code{mode = NULL} selects chained estimation. Chained estimation runs for
built-in doubly robust \code{est_method = "dr"} (the default), explicit
built-in outcome regression \code{est_method = "reg"}, built-in
\code{est_method = "ipw"} (chained IPW-GMM, Bellego, Benatia and
Dortet-Bernadet 2025) or a did-style custom IPW panel callback function
\code{f(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...)} returning a
list with finite scalar \code{ATT} and finite length-\eqn{m_\ell} vector (or
1-column matrix) \code{att.inf.func} in the supplied row order. Unsupported
\code{est_method} strings are refused with an explanatory error before the data
is touched; built-in IPW is never substituted for them. For repeated cross
sections (\code{panel = FALSE}) the default resolves to ordinary mode.
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
  gmm_weighting = c("optimal", "identity"),
  gmm_covariance = "tikhonov",
  pt_assumption = c("post", "all"),
  pt_start = NULL,
  control = NULL
)
}
\arguments{
\item{yname, tname, idname, gname, xformla, data, panel, control_group,
  anticipation, weightsname, fix_weights, alp, bstrap, cband, biters,
  clustervars, est_method, base_period, faster_mode, print_details, pl, cores,
  compute_inffunc}{As in \code{did::att_gt()}, with the same names, order and
  defaults. In ordinary mode, arguments are forwarded only if you supplied
  them, so did applies its own defaults. Each expression is evaluated once,
  in your environment, when did needs it. The exception is \code{panel}:
  with \code{mode = NULL}, cdid evaluates it first to choose the mode. The
  default \code{est_method = "dr"} is kept. In chained mode \code{idname} is
  required (units are aligned by identity) and at most one cluster variable
  beyond \code{idname} is accepted. \code{compute_inffunc = FALSE} returns
  point estimates only, as in did: \code{bstrap} and \code{cband} become
  \code{FALSE}, \code{inffunc}, \code{V_analytical}, \code{W} and
  \code{Wpval} are \code{NULL}, \code{se} is \code{NA}, \code{c} is the
  pointwise normal constant and \code{aggte()} refuses the result; the
  feasible two-step/Tikhonov GMM weights still use the complete internal link
  influence functions, so every estimate equals the \code{compute_inffunc =
  TRUE} estimate. \code{pl = TRUE} with \code{cores = k} runs the multiplier
  bootstrap's draw products on \code{k} worker processes (forked on
  Unix-alikes, a PSOCK cluster on Windows; base \code{parallel}), which must
  be declared through \code{control = list(workers = k)}; otherwise the
  request is refused before the data are read. Signs are drawn by your session
  in the serial batch order, so draws, standard errors and bands equal
  \code{cores = 1}, later \code{aggte()} calls inherit the workers, and
  estimation itself stays serial. \code{print_details = TRUE} reports stage
  and comparison progress through messages (\code{suppressMessages()}
  silences it); \code{FALSE} keeps the fit quiet. \code{faster_mode} keeps
  its documented meaning: chained mode has one execution engine. The
  public never-treated code is \code{gname == 0}. Chained \code{fix_weights}
  accepts \code{NULL}, \code{"first_period"} and \code{"base_period"}
  (\code{"varying"} is refused); see \sQuote{Reference dates and changing
  weights} for covariates and sampling weights that change within a unit.}

\item{allow_unbalanced_panel}{\code{NULL} (default), \code{TRUE} or
  \code{FALSE}. In ordinary mode, \code{NULL} means did's default
  \code{FALSE} (balance the panel). Chained estimation retains every eligible
  unit (at least two observed outcomes), including balanced panels, so an
  explicit \code{FALSE} in chained mode is an error that points to
  \code{mode = "ordinary"}.}

\item{...}{Passed to a custom \code{est_method} function in ordinary mode (as
  in did) and in chained mode when \code{est_method} is a function. When
  \code{est_method} is a character string (such as \code{"ipw"}) in chained
  mode, non-empty \code{...} is refused.}

\item{mode}{\code{NULL} (default), \code{"chained"} or \code{"ordinary"}.
  \code{NULL} selects \code{"chained"} when \code{panel = TRUE} and
  \code{"ordinary"} for repeated cross sections.}

\item{selection_universe}{Selection logit fitting sample: \code{"F"} (default,
all supplied representative units), \code{"B"} (at least one observed outcome),
or \code{"A"} (at least two). \code{"B"} and \code{"A"} are supported only for
\code{selection_model = "pair_logit"}; \code{selection_model} = \code{"independent_period_logit"}
requires \code{selection_universe = "F"}. These outcome counts do not specify structural
observation windows. The entire F remains the population and inference frame.}
\item{selection_rho}{NULL for rho=1, or a column name in the unit frame containing
fixed positive KNOWN normalized rho = f_P(X|G,U)/f_P(X|G) on U. Only supported
for A/B under \code{selection_model = "pair_logit"}; no estimation or self-normalization. Missing X/rho outside U is allowed
and never evaluated. Necessary covariates must be complete on U.}
\item{selection_justification}{A nonempty explanation required for A/B. With
NULL rho this explicitly asserts U is representative within cohort under the
declared weighted population law (U independent of X conditional on G).
Otherwise it declares provenance of the known normalized transport. The
package records this assertion; fit diagnostics cannot establish it.}
\item{selection_model}{Observation-selection probability model when selection
correction is active: \code{"pair_logit"} (default, direct logit for observing
both endpoints of an ordered pair) or\cr
\code{"independent_period_logit"}
(period-observation logit with pair probability formed as the product of the two
endpoint date probabilities on \code{selection_universe = "F"}). Ignored when
\code{selection = "none"} unless a non-default value is supplied without an
active selection assumption, which is refused. Inactive and never evaluated
under \code{selection_est_method = "imputation"}; active under
\code{selection_est_method = "dr"}, where it defines the observation-probability
arm of the augmented correction.}
\item{selection_pooling}{Slope-pooling mode across active selection cells for
both \code{"pair_logit"} and \code{"independent_period_logit"}:
\code{"separate"} (default, cell-specific slopes), \code{"within_cohort"}
(one shared slope vector per cohort across its active pairs or dates, with
cell-specific intercepts when the selection formula includes an intercept), or
\code{"across_cohorts"} (one shared slope vector across all active cells, with
cell-specific intercepts when the selection formula includes an intercept).
Pooling restricts selection slopes only; shared outcome-regression slopes and
local scientific-block intercepts are unchanged. Inactive and never evaluated
under \code{selection_est_method = "imputation"}, which uses
\code{imputation_pooling} instead. Under \code{selection_est_method = "dr"} it
pools the observation-probability slopes and is independent of
\code{imputation_pooling}.}
\item{selection_est_method}{Chained-only observation-correction method, separate
from the treatment comparison method \code{est_method}: \code{"ipw"} (default,
the existing probability-weighted correction, unchanged), \code{"imputation"}
(pure actual-change imputation; see the section below) or \code{"dr"} (augmented
observation correction combining both; see the section below). \code{"imputation"}
requires \code{selection = "conditional_mtar"} (or
the legacy assumption alias \code{"corollary1"}), \code{selection_universe = "F"}
with \code{selection_rho = NULL}, known cohorts and complete \code{xformla}
covariates on every F row, and a built-in \code{est_method}; it fits no
observation-probability model, so \code{selection_model},
\code{selection_pooling} and \code{selection_xformla} are inactive and never
evaluated. \code{"dr"} has the same declaration requirements, is implemented
with \code{est_method = "dr"} only (\code{"ipw"}, \code{"reg"} and custom
callbacks are refused before the data are evaluated; nothing is substituted)
and uses every probability-model and imputation argument. Ordinary mode neither
evaluates nor forwards it.}
\item{imputation_pooling}{Slope pooling of the actual-change imputation model,
used only when \code{selection_est_method = "imputation"} or \code{"dr"}
(otherwise it is never evaluated): \code{"separate"} (default; one slope vector per actual cohort/pair
cell), \code{"within_cohort"} (one slope vector per cohort across its needed
pairs) or \code{"across_cohorts"} (one slope vector across all needed cells).
Local cohort/pair intercepts are kept whenever \code{xformla} has an intercept.
Pooling is an explicit model restriction, never an automatic fallback; it can
be wrong, for example when treated cohorts' dynamics differ across pairs.
Ordinary mode neither evaluates nor forwards it.}
\item{frame, selection, selection_xformla, gmm_weighting, gmm_covariance, pt_assumption, pt_start}{
  Chained-only: an optional unit frame with population or sampling information
  (one row per population unit, matched to \code{data} by \code{idname};
  frame-only units keep their cohort, weight and cluster for population shares
  and influence-function rows), the observation-selection assumption and its
  covariate formula (\code{selection = "none"} by default, or
  \code{"conditional_mtar"} for conditional mean trend-MAR selection correction,
  with \code{"corollary1"} accepted as the legacy direct-pair alias for
  \code{selection = "conditional_mtar"} with \code{selection_model = "pair_logit"};
  \code{selection_xformla} is inactive and never evaluated under
  \code{selection_est_method = "imputation"}; under \code{"dr"} it defines only
  the observation-probability design, \code{NULL} inheriting \code{xformla}), the
  GMM weighting (\code{"optimal"} by default for feasible two-step weighting,
  with \code{"two_step"} accepted as an explicit alias for \code{"optimal"},
  and \code{"identity"} for explicit unweighted estimation), the singular or
  near-singular moment-covariance handler under feasible two-step weighting
  (\code{gmm_covariance = "tikhonov"} by default for covariance-ridge Tikhonov
  regularization; Moore-Penrose is deferred and is rejected), and the separately selected parallel-trends assumption
  (\code{"post"} by default; \code{"all"} adds untreated-trend restrictions
  from an explicit common \code{pt_start}). Under \code{"all"}, \code{pt_start}
  must be one finite numeric date exactly in the stored calendar, chosen
  exogenously; omission, NULL and non-calendar dates are errors. Under chained
  \code{"post"}, a non-NULL start is inapplicable and refused. Ordinary mode delegates before chained-only
  validation: these arguments are neither evaluated nor forwarded to did.
  \code{pt_assumption} is separate from \code{base_period}, whose did meaning
  and default (\code{"varying"}) are unchanged.}

\item{control}{Chained-only advanced list, \code{NULL} by default. \code{NULL}
  (and an empty list) resolves to the adopted defaults recorded in the result:
  \code{numeric_workspace_bytes = 1073741824} (1 GiB), \code{workers = 1}, and
  no user estimated-work budget (\code{counted_work_limit = Inf}). The
  estimated-work cap is opt-in: set
  \code{counted_work_limit} to a finite positive whole number (exact binary64
  integer) to refuse, before the stage runs, a fit whose conservative preflight
  work ESTIMATE exceeds it; an explicit \code{Inf} also opts out. Zero,
  negative, \code{NA}/\code{NaN}, non-integral, out-of-range and non-numeric
  values are refused. Work estimates are diagnostic and conservative, not
  measured runtime. Internally a disabled cap hands the finite private
  constructors the exact-count arithmetic ceiling \code{2^53 - 1}; the result
  records the public value and that ceiling separately
  (\code{chained$public$counted_work_limit}, \code{chained$public$work_cap}).
  These controls admit counted numeric payload and work estimates; they do
  not bound process RSS, R heap, CPU time or whole-call work, and cdid does not
  detect machine RAM or time. A plain named list may also override
  \code{numeric_workspace_bytes}, \code{workers} (the number of multiplier
  bootstrap worker processes \code{pl = TRUE, cores = k} may use; the
  default 1 is never raised automatically and their whole buffers count
  against the same numeric workspace), \code{chunk_rows},
  \code{condition_trigger} (default \code{1e-10}, finite and strictly between
  zero and one), and \code{ridge_relative} (default \code{1e-4}, finite and positive).
  Unknown entries and objects are refused. The universal per-fit error
  guarantee and its \code{relative_accuracy} override are deferred; supplying
  that entry raises \code{cdid_chained_public_accuracy_deferred}. Ordinary mode never evaluates it.}
}
\value{
In ordinary mode, the \code{MP} object returned by \code{did::att_gt()},
unchanged. did records the delegated call (for example
\code{did::att_gt(yname = yname, ...)}) as \code{DIDparams$call}, so
printed summaries show that call.

In chained mode with \code{est_method = "dr"}, \code{est_method = "reg"}, \code{est_method = "ipw"} or a custom IPW callback
function, a \code{cdid_MP} object (inheriting \code{MP}, built with did's
constructors) whose did fields hold the reported chained cells: direct
did-style pretrend comparisons for every declared treated cohort and the
ATT(g,t) post cells, with complete full-frame influence functions, the analytic
(or multiplier IQR) standard errors and the critical value of the requested
interval law. For a custom \code{est_method} callback, each qualified link
invokes \code{est_method(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...)}
on the observed risk set in canonical increasing A-unit order with calibrated
weights \code{a_i = w_i / max(r_G, 1e-6)}, composes cohort-rate estimation uncertainty via
\code{b_k = -sum(psi[G == k]) / (m * r_k)} (zero when \code{r_k <= 1e-6}) and
\code{phi = (n / m) * 1[observer] * psi + rho \%*\% b} without adding a duplicate
propensity or pre-GMM cohort-share term. A callback may also return the optional
\code{log_weight_gradient}: a finite length-\eqn{m_\ell} vector (or
1-column matrix) in the supplied row order with
\eqn{g_j = \partial ATT / \partial \log} \code{i.weights[j]}, the total derivative
with every internal refit included (ATT units, never centred). When supplied the
rate channel uses \code{b_k = -sum(g[G == k]) / r_k} instead of the identity
rule (the two are never added). Under \code{selection = "conditional_mtar"}
(\code{selection_est_method = "ipw"}, every selection model, pooling and
universe) callbacks receive the selection-corrected weights \code{w a_U / (rho max(p_U, 1e-6))} and must
supply this gradient; \code{b} then loads the selection rates and logit
coefficients and \code{phi = (N / m) * 1[observer] * psi + xi_u \%*\% b} on every
original row. Without it such requests are refused: \eqn{g = \psi/m} needs
prerequisites that are not verified and \eqn{g} is never inferred. Shape checks
do not establish the callback's scientific validity. Callbacks must reimplement the existing
IPW criterion with the same weighted propensity equations, normalized outcome
means and complete fixed-rate influence function, including propensity and
normalization terms. The IF must have valid mean-zero semantics and preserve
the supplied row order; finite return shapes do not establish correctness.
The bridge retains the callback IF unchanged, records its mean diagnostically,
and imposes no centering cutoff. The existing point/design/complete-IF forest
identities apply after outcome-free verification of equal cohort/control pools
and observation masks, including balanced panels. Broader arbitrary-callback
support requires a separate public policy. Arbitrary
workspace inside the opaque callback is outside package-owned array accounting.
Cells the retained comparison catalog cannot support are \code{NA} with a
support reason; they are not population non-identification claims. The
\code{chained} element records the support ledger, the GMM weighting and
covariance-handling record (\code{covariance_handling}), the record that complete
per-fit error certification was not requested or performed, the pretrend
test (\code{W}/\code{Wpval} are set only when the full declared test is
computable; otherwise an explicit unavailable or not-applicable status and
reason), the resolved resource profile with its status, the actual input
dimensions and derived finite bounds, the resource accounting of reachable
numeric payload (not process memory or CPU), and the state later aggregation
needs so that \code{\link{aggte}} never refits. The resource accounting
records the exact auxiliary reporting plan (\code{auxiliary_plan}: the union
of the varying and universal pretrend comparisons that are fitted, the
reported pretrend columns of the chosen \code{base_period}, reference-zero and
unavailable cells), the mask registry's charged scans and retained records
(\code{mask_registry}) and, on the built-in IPW reduction route, whether the
registry's sealed mask-equivalence certificate replaced streamed mask
comparisons (\code{early_structural$mask_evidence}). The inherited flags as
honoured, the admitted multiplier workers and their route are recorded in
\code{chained$request} (\code{execution}); \code{DIDparams} carries the
requested \code{print_details}, \code{pl}, \code{cores} and
\code{compute_inffunc}.

Refusals are explicit errors inheriting \code{"cdid_error"}:
\itemize{
  \item \code{"cdid_chained_unavailable"}, also inheriting
    \code{"cdid_chained_request_est_method"}, for unsupported string chained requests.
  \item \code{"cdid_chained_unavailable"}, also inheriting
    \code{"cdid_chained_request_unsupported_selection"}, for the not implemented
    \code{selection_est_method = "dr"} with \code{est_method = "ipw"} or \code{"reg"}.
  \item \code{"cdid_custom_ipw_error"} for malformed or nonfinite custom callback outputs.
  \item \code{"cdid_chained_panel"} if \code{panel} is not \code{TRUE}.
  \item \code{"cdid_chained_balanced"} for an explicit
    \code{allow_unbalanced_panel = FALSE}.
  \item \code{"cdid_chained_public_control_error"} for an invalid \code{control}.
  \item \code{"cdid_chained_request_parallel_error"} for
    \code{pl = TRUE, cores = k} beyond the declared \code{control$workers}.
  \item \code{"cdid_chained_request_*"} for unsupported requests or inputs.
  \item The accepted pipeline's typed bound errors when the input's conservative
    numeric-payload or work envelope exceeds the profile (nothing is dropped or relaxed),
    and genuine invalid nuisance, nonfinite, covariance or solve errors.
}
Exact signed-coordinate separation is a sufficient finite-root
obstruction; it is not a complete separation
classifier or a population identification conclusion. An invalid \code{mode}, or an
invalid \code{panel} when resolving \code{mode = NULL}, raises
\code{"cdid_mode_error"}.
}
\details{
PT-all retains every baseline PT-post condition, including at qualified causal
dates before the chosen start. Each cohort keeps its original admissible
reference date. Added equalities constrain represented operator-specific
untreated pre-period differentials from the common start up to that reference;
earlier differentials and anticipation effects stay free. The same comparison
catalog, control operators, observation pools and complete influence functions
are used. Sampling, selection, overlap and nuisance assumptions remain required.
If there are no represented extra restrictions globally, the fit is allowed and
summary/print reports ``No additional restrictions'' (the PT-post-equivalent
model). Stored window dates and represented restrictions are counted separately
by cohort in \code{result$chained$pt}; no cohorts are dropped.

Both \code{base_period} values retain did's direct, unrestricted pre-period
reporting. Restricted zero pre-parameters never supply the full pretrend Wald test.
Descriptive pre-start, within-window and anticipation tags accompany the test
record. Pre-start contrasts may reject without contradicting PT-all from the
chosen start; anticipation effects need not be zero. An in-window auxiliary
comparison is implied zero only under its actual control operator's PT and
selection premises; a date tag alone does not cover an auxiliary NYT pool absent
from the catalog assumptions. There is no new subset or overidentification test.

GMM weighting and covariance handling: feasible two-step GMM weighting
(\code{gmm_weighting = "optimal"}, or alias \code{"two_step"}) is the default,
and \code{gmm_weighting = "identity"} is the explicit unweighted option.
Proved exact structural equivalence reductions are applied before covariance
conditioning; when structural singularity is removed by exact reduction alone
and the reduced transformed covariance \eqn{B = C \widehat{V}_r C'} (with
\eqn{G = E'E = C'C} from the sparse reduction reconstruction \eqn{E}) is well
conditioned (\eqn{\rho \ge \tau_{\mathrm{cond}}}), ordinary Cholesky-whitened
QR is used without invoking a singular-covariance handler. When \eqn{B} is
singular or near-singular (\eqn{\rho < \tau_{\mathrm{cond}}}), the selected
handler \code{gmm_covariance = "tikhonov"} uses covariance ridge with
\eqn{\lambda = \alpha_{\mathrm{ridge}} \operatorname{tr}(B) / M} and
effective reduced weight \eqn{C'(B + \lambda I_L)^{-1}C}, where \eqn{M}
is the full catalog size. The trigger uses the reciprocal spectral condition
\eqn{d_{\min}(B)/d_{\max}(B)}, rather than a 1-norm condition estimate.
Known unit-covariance rank deficiency or failed ordinary Cholesky on valid PSD
covariance also activates ridge. Identity is never substituted automatically.
At exactly zero covariance the criterion is full-catalog least squares (reduced
weight \eqn{G} up to scale), with the actual zero sandwich and the existing
zero-scale inference diagnostics. No moments, units or identified contrasts
are dropped. Internally formed complete-unit-IF Gram roundoff is diagnosed
under the existing eigen-residual rule; no eigenvalues are clipped and no penalty
is increased after failure. Nonfinite, materially indefinite and invalid
provenance inputs still refuse. Tuning, activation, full/retained dimensions,
scale, conditioning and provenance are recorded in the result.
Analytical and multiplier-bootstrap inference always
use the complete target influence functions and the true sandwich variance,
never the inverse regularized criterion matrix. Clusters affect inference only;
weights and tuning use unit IF covariance. Feasible weighting does not imply
finite-sample unbiasedness. Inference requires the complete-link ALR/CLT,
covariance consistency and a stable limiting loading/automatic branch; no
uniform boundary or zero-population-trace inference guarantee is claimed.

Numerical policy: conventional Cholesky/QR computation uses the existing
conditioning and factor-residual policies. Relative solve residuals remain
diagnostic, including at zero dual solutions. Complete per-fit interval error
certification is not requested or performed; no universal 1e-6-of-SE guarantee
is claimed. Its public override is explicitly deferred.

Pretrend test: the full set of declared nonreference pretrend
comparisons of the selected \code{base_period} reporting is tested with the
did-style Wald statistic when every comparison is available and the joint
covariance admits the full Cholesky/solve under the existing conditioning policy; otherwise the test is unavailable
with its reason and the cells are returned unchanged. It is the unconditional
test; the conditional pretrend test is
\code{\link{conditional_did_pretest}()}.

The package refuses to load with did older than 2.5.1 or DRDID older than
1.3.0; it has been tested with did 2.5.1 and DRDID 1.3.0, and later versions
are accepted but untested. Qualify calls as \code{cdid::att_gt()} so that a
later \code{library(did)} cannot change which function runs.

Input and result inspection: in chained mode, \code{data} is a long panel with
unique \code{(idname, tname)} rows, finite numeric \code{tname}, numeric
\code{yname} (finite when observed and kept as \code{NA} when missing, never
replaced with \code{0}), and constant unit cohort \code{gname} (\code{0} for
never-treated units). No unit requires a complete history or an observed outcome
at the cohort baseline date. The optional unit table \code{frame} defines the
study population \eqn{F} of size \eqn{N = |F|}, retaining units observed at two
or more dates (\eqn{A}), singletons (one observed outcome) and zero-outcome
units in population counts, cohort shares, active selection or imputation models
and influence functions for all study units without dropping rows. Under the default
\code{selection = "none"}, units observed at both endpoints of each used date
pair are assumed representative of their cohort in the study population with
respect to the mean outcome change over that pair. Use \code{summary()},
\code{\link{tidy}}() (columns \code{target.id}, \code{role}, \code{available},
\code{support.status}, \code{support.reason}), \code{\link{glance}}(),
\code{\link{aggte}}() and \code{\link{ggdid}}() to inspect fitted results; see
\code{\link{reexports}} and \code{vignette("cdid", package = "cdid")}.
}
\references{
Bellego, C., Benatia, D. and Dortet-Bernadet, V. (2025). The Chained
Difference-in-Differences. \emph{Journal of Econometrics} 248, 105783.
\doi{10.1016/j.jeconom.2024.105783}

Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with
Multiple Time Periods. \emph{Journal of Econometrics} 225(2), 200--230.
\doi{10.1016/j.jeconom.2020.12.001}
}
\section{Reference dates and changing weights (chained mode)}{
Calls whose covariates and sampling weights pass the existing static checks
(constant within unit and in agreement with \code{frame}) run unchanged, and
\code{fix_weights = NULL}, \code{"first_period"} and \code{"base_period"} give
identical results. Sampling weights that change within a unit require
\code{"first_period"} (every context uses \eqn{w} at the first stored calendar
date \eqn{t_1}) or \code{"base_period"} (every unit of target cohort \eqn{g}'s
context, controls included, uses \eqn{w} at \eqn{b_g = \max\{d: d +
\code{anticipation} < g\}}); \code{NULL} would mean did's period-specific cell
weights and is refused. With unit-level covariates and \eqn{w_{t_1}} present
on every study unit, \code{"first_period"} gives the existing static estimator on
\eqn{w_{t_1}}, on every route.

Covariates that change within a unit (or a formula use of \code{tname}) are
reference-date covariates, supported on every chained route: baseline IPW,
custom did-style callbacks, OR and DR (\code{selection = "none"}); conditional
selection with IPW, OR or DR (\code{pair_logit} on \code{selection_universe}
F, B or A with known \code{selection_rho}, \code{independent_period_logit} on
F, every \code{selection_pooling}); pure actual-change imputation with IPW, OR
or DR (every \code{imputation_pooling}); and the augmented observation
correction (\code{selection_est_method = "dr"}). Every post comparison of
target \eqn{g} (bridges and controls included) uses each unit's values at
\eqn{b_g}; an auxiliary pre-period report uses the earlier endpoint of its
date pair. Values come only from the actual (id, date) rows, including rows
with a missing outcome, and are required only where an estimate reads them:
the outcome/treatment design on observed pairs (for every study unit in each
cohort for pure imputation and the augmented correction), the selection design
in the selection fitting sample, and \eqn{W_g} on each rate or mass denominator. A
missing required value is refused with its id, date, variable and support
(never carried forward or substituted). Each formula, its factor levels and
transform state are evaluated once on the canonical required (date, id) rows
of all post and both auxiliary report families; \code{selection_xformla} keeps
its own design. Rates, repeat rates, selection cells, treatment logits,
imputation cells and shared-OR blocks are keyed by recipe (design date,
weight anchor) and fitted once; the shared OR keeps one global slope on
\eqn{L_e z}, and selection or imputation pooling never crosses recipes. Custom
callbacks receive the recipe design and \eqn{W_g / r_G}. Every influence
function has one row per original unit with all rate, selection, imputation,
treatment and outcome-regression channels. A frame value cannot override data
values that change. On these fits the dated baseline IPW and custom routes use
the general GMM criterion over every catalog comparison rather than the static
structural reduction.

For these reference-input fits, \code{aggte()} weights cohorts by their
first-calendar masses \eqn{H(I_g w_{t_1})}, with share uncertainty in the
influence functions; never-treated units need no \eqn{w_{t_1}}, and an aggregate
whose treated cohort lacks \eqn{w_{t_1}} for some unit is \code{NA} with the ids.
Each cell keeps its context weight; with \code{"base_period"} an aggregate is a
declared mixture of cohort effects, not the ATT of one weighted population.
}

\seealso{\code{\link{aggte}}, \code{\link{conditional_did_pretest}}, \code{\link{tidy}}, \code{\link{glance}}, \code{\link{ggdid}}, \code{vignette("cdid", package = "cdid")}}
\examples{
mpdta <- did::mpdta
set.seed(1)
out <- cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                    gname = "first.treat", xformla = ~lpop, data = mpdta,
                    mode = "ordinary", bstrap = FALSE)
summary(out)
cdid::aggte(out, type = "dynamic")

# For panel data the default mode is "chained" and default est_method is "dr".
# Unsupported est_method strings are refused rather than replaced by IPW:
res <- try(cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                        gname = "first.treat", data = mpdta,
                        est_method = "stdipw"), silent = TRUE)
cat(res)

\donttest{
# Chained doubly robust (default est_method = "dr") and IPW-GMM:
ch_dr <- cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                      gname = "first.treat", xformla = ~lpop, data = mpdta,
                      bstrap = FALSE, cband = FALSE)
summary(ch_dr)
cdid::aggte(ch_dr, type = "dynamic")

# Chained IPW-GMM with analytic inference and the default controls (no user
# estimated-work cap; 1 GiB counted numeric workspace; one worker).
ch <- cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                   gname = "first.treat", xformla = ~lpop, data = mpdta,
                   est_method = "ipw", gmm_weighting = "optimal",
                   gmm_covariance = "tikhonov", bstrap = FALSE, cband = FALSE)
summary(ch)
cdid::aggte(ch, type = "dynamic")   # no refit

# OPTIONAL finite estimated-work cap (refuses before a stage whose conservative
# work estimate exceeds it; not a CPU/time guarantee) and the existing
# configurable counted numeric workspace (not process RSS):
ch_capped <- cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                          gname = "first.treat", xformla = ~lpop, data = mpdta,
                          est_method = "ipw", bstrap = FALSE, cband = FALSE,
                          control = list(counted_work_limit = 1e12,
                                         numeric_workspace_bytes = 1073741824))

# Explicit identity weighting (two_step is an alias for optimal):
ch_id <- cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                      gname = "first.treat", xformla = ~lpop, data = mpdta,
                      est_method = "ipw", gmm_weighting = "identity",
                      bstrap = FALSE, cband = FALSE)

# Pure actual-change imputation as the observation correction (no
# observation-probability model; selection_model etc. are inactive):
ch_imp <- cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                       gname = "first.treat", xformla = ~lpop, data = mpdta,
                       est_method = "dr", selection = "conditional_mtar",
                       selection_est_method = "imputation",
                       imputation_pooling = "separate",
                       bstrap = FALSE, cband = FALSE)
ch_imp$chained$selection$imputation_pooling

# Augmented observation correction (selection_est_method = "dr" with
# est_method = "dr"): impute the actual pair changes and correct the observed
# residuals by fitted observation probabilities. Its observation logits need
# every needed pair to be partially observed, so this illustration removes a
# deterministic 1-in-7 pattern of county-years (mpdta itself is balanced).
mp_obs <- mpdta[(mpdta$countyreal + mpdta$year) \%\% 7 != 0, ]
ch_aug <- cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                       gname = "first.treat", xformla = ~lpop, data = mp_obs,
                       est_method = "dr", selection = "conditional_mtar",
                       selection_est_method = "dr",
                       imputation_pooling = "within_cohort",
                       selection_model = "pair_logit",
                       selection_pooling = "within_cohort",
                       bstrap = FALSE, cband = FALSE)
ch_aug$chained$selection[c("method", "selection_model", "selection_pooling",
                           "imputation_pooling", "selection_formula_inherited")]
}
}

\section{Conditional MTAR selection correction}{
Use \code{selection="conditional_mtar"} (or the legacy direct-pair alias
\code{selection="corollary1"}) with built-in \code{est_method="dr"},
\code{est_method="ipw"}, or \code{est_method="reg"}.
\code{selection_xformla=NULL} inherits \code{xformla}; both formulas are
one-sided, fixed unit-level designs preserving \code{stats::terms},
\code{assign} and intercept metadata, including intercept-only (\code{~1}),
no-intercept (\code{~0 + x}) specifications. Only \code{selection_xformla}
supports a zero-column design (\code{~0}, fixed cell probability 0.5 and date
product 0.25 with zero parameters); \code{xformla} requires a nonempty design.
\code{selection_model="pair_logit"} (default) fits the joint observation indicator
for each active ordered pair and cohort; \code{selection_model} = \code{"independent_period_logit"}
fits period-observation logits on \code{selection_universe="F"} and forms each
pair probability as the product of its two endpoint date probabilities without
A/B multiplication. Both models support \code{selection_pooling="separate"}
(default), \code{"within_cohort"} and \code{"across_cohorts"} over active union-catalog
cells. The legacy alias \code{selection="corollary1"} is restricted to
\code{selection_model="pair_logit"}. Identical actual models are shared across
NT/NYT pools and direct pre-comparisons, and reused date or pair scores are
summed by original unit before covariance computation. This is logit only.

F is the supplied representative unit table, in frame order. When frame is
omitted, the units listed in the supplied panel data, including singletons, define F;
no wider population register is presumed. External linked frame-only rows
may have zero outcomes but need known cohort, positive fixed weight (if
weighted), necessary covariates and cluster labels. An absent outcome is a
known pair nonresponse only under the user's linkage/observation coverage.
No frame rows are dropped. Singletons supply no differences.

For U=F, the combined correction is a/p. Retained U=A/B uses its own
conditional logit p_U and repeat rate a_U with b=a_U/(rho*p_U), the fitted
probability (or date product) floored at 1e-6 as described under Probability
clipping. IPW uses b in
the treatment score and both empirical outcome means; outcome regression uses
b in the control-only weighted regression and the target residual mean. The a_U denominator is
U; for A, a_U=1 and its rate influence is zero. Known rho contributes no
nuisance influence. Its normalization is a scientific input, not a sample
renormalization. F-logit correctness does not imply U-logit correctness;
these are different modeling assumptions, without an efficiency ranking.

Every actual catalog pair, including pre-comparisons, needs its own population
conditional mean trend-MAR restriction. Adjacent-pair MAR is insufficient.
Selection and, for IPW, synthetic treatment-odds specifications must be correct
conditional on all scientific X; OR requires the correct control-mean regression; different predictor lists require separately
sufficient restrictions. PT/no anticipation, weighted-law validity, positivity,
finite regular full-rank interior roots and complete-score second moments
also remain application assumptions. No trimming, imputation or new nuisance
regularizer is used; the probability clipping is described in its own section.

Complete influence functions replay every F row, including nonobserver
selection/rate scores. Public n is N=|F|. The GMM core uses psi=(n_A/N)*phi_F
on all F rows with its original n_A covariance scale. Population shares,
share/ATT cross covariance and all cluster rows are retained, including
clusters with no A units. Clustering changes inference only. Existing GMM
asymptotic loadings apply; they are not exact case derivatives of refitted
two-step covariance weights. Full requested pretrend nulls retain their
computed result or explicit unavailable reason.

Known pair probabilities (including design p=1), custom-IPW callbacks without
an explicit \code{log_weight_gradient}, estimated/unnormalized transport and
probit are unsupported.
All-success data never establish a design probability: fitted logit refuses
separation/boundary/rank failures. Numerical certificates are not provided
for this branch. Finite tests and the demo do not establish inference coverage.
}

\section{Control-mean outcome regression and doubly robust estimation}{
Explicit \code{est_method="reg"} uses the built-in panel shared one-period
outcome-change control-mean regression with \code{selection="none"} or explicit
conditional-MTAR (\code{selection="conditional_mtar"} / legacy \code{"corollary1"})
selection. Default or explicit \code{est_method="dr"} composes the treatment-DR
extension of this shared one-period outcome-change control-mean regression and
residualized Hajek IPW over deduplicated risk-set treatment propensity logits
from the same public \code{xformla}, keeping the full un-reduced treatment and
selection designs while applying OR alias preprocessing only to the outcome
regression recipe. For \code{est_method="reg"}, no treatment propensity or treated outcome
regression is fitted. Across all qualified catalog and auxiliary comparison
links, a single joint weighted least-squares regression (\code{n_regression_fits = 1})
fits one shared slope vector \eqn{\beta} on the duration-scaled covariate design
\eqn{d_{e,i} = L_e z_i} together with one local intercept \eqn{\alpha_b} per
unique pair/control-pool block \eqn{b} (\code{n_regression_blocks}) when
\code{xformla} includes an intercept. Here \eqn{L_e = i_t - i_s} is the
stored-calendar transition count between ordered calendar positions
\eqn{i_s < i_t}, never raw subtraction \eqn{t - s} of arbitrary numeric date
labels. Baseline weights are the fixed unit weight divided by that cohort's
actual-pair response rate among repeat observers; Corollary1 uses the same
reader, logit and known transport declaration described above. Each target link
is the normalized weighted target-cohort residual mean
\eqn{\hat{\mu}_\ell = \sum_{i \in S_\ell} q_{\ell,i}(\Delta Y_{e,i} - \hat{\alpha}_{b(\ell)} - L_{e(\ell)} z_i'\hat{\beta}) / \sum_{i \in S_\ell} q_{\ell,i}}.
All original comparisons and fixed-horizon moments enter the existing GMM,
sandwich, aggregation, multiplier and unrestricted pretest paths.

Actual \code{stats::terms} and \code{assign} metadata determine whether the
outcome formula includes an intercept and map retained columns into the outcome
recipe. Intercept-only (\code{~1}) fits local block intercepts with zero shared
slopes; no-intercept (\code{~0 + x}) fits raw uncentered shared slopes without
local block intercepts (\code{n_local_intercepts = 0}); an explicitly named
\code{`(Intercept)`} variable is preserved as an ordinary shared slope when the
formula has no intercept. When a constant slope column is aliased with a
formula intercept, the intercept is retained so the model preserves local block
intercepts plus duration-scaled slopes.

A common outcome-design basis is frozen before fitting or differentiation.
Standard R \code{lm.fit} pivoted QR (default tolerance \code{1e-7}) proposes
aliases. Each proposal is checked on every repeat-observer control and treated
prediction row. The maximum absolute reconstruction discrepancy is recorded;
the maximum scaled discrepancy is checked at arithmetic roundoff
(\code{64 * .Machine$double.eps} times the original column count).
This is a numerical diagnostic, not an exact-arithmetic or effect-error
certificate. A proposal with nontrivial reconstruction discrepancy retains that
original column. Selection predictor columns are unchanged. In the joint
pooled control fit, any remaining pooled-control null direction is accepted
only if it reconstructs all active control and treated prediction rows to
arithmetic roundoff; a locally rank-deficient control block remains supported
when the pooled control design across blocks identifies the shared slopes. The
compact pooled control fit separately checks its final reduced QR relative
diagonal against the recorded R regression tolerance. This criterion refuses
globally rank-deficient or near-collinear reduced fits without dropping further
columns, observations or moments. Original, retained and aliased names and the
frozen reconstruction maps are recorded in \code{chained$outcome_design} and
\code{chained$nuisance_diagnostics$basis}; source formula and its environment
remain in \code{DIDparams$xformla}. A dropped coordinate is not an identifying
assumption that its original coefficient is zero.

Outcome covariates need only be complete on repeat observers. Baseline uses a
reader for all study units without fitting any selection model; Corollary1 retains
all F identities, including off-A nuisance effects. Every streamed IF link
batch includes global shared-slope dependencies on all F rows, and cluster
inference sums original unit rows. Public diagnostics (\code{chained$nuisance_diagnostics})
report the shared slopes, block local intercepts, fitted catalog, duration
convention, global QR factor (retained once) and compact per-block records
(describing local intercept blocks of the single joint fit, not independent
regression fits), charging the complete diagnostic copy against the existing
finite resource limits. Control-only support/rank and conditioning failures are
computation failures, not population ATT nonidentification claims. No public OR
tolerance override, automatic fallback or numerical accuracy certificate is supplied.
}

\section{Pure actual-change imputation}{
With \code{selection = "conditional_mtar"} and \code{selection_est_method = "imputation"}
the observation correction predicts instead of weighting. For every actual
cohort/pair cell \eqn{(k, e)} needed by the union of the post catalog and both
base-period auxiliary catalogs, the actual-change model
\eqn{j_{k,e}(X) = \alpha_{k,e} + L_e z(X)'\beta_{r(k,e)}} is fitted by weighted
least squares to the observed pair changes (\eqn{S_e = 1}) of cohort \eqn{k}
with the fixed unit weights (no observation probabilities), with one local
intercept per cell when \code{xformla} has an intercept and slopes shared as
\code{imputation_pooling} declares. Each fitted \eqn{j_{k,e}} is evaluated on
every F row of cohort \eqn{k} (prediction standardization, not a fill-in of the
unobserved changes only), with the empirical repeat rates \eqn{a_k} on F.

The treatment method then uses these reconstructions: \code{est_method = "ipw"}
fits the risk-set treatment propensity logits with synthetic weights
\eqn{w a_k} and no counterfactual outcome regression; \code{est_method = "reg"}
fits the shared one-period counterfactual control outcome-change regression
(local block intercepts when declared by \code{xformla}, with slopes shared
across blocks on the retained basis) to the imputed synthetic control
changes and no treatment logit; \code{est_method = "dr"} fits both. The
counterfactual regression is a second model, distinct from the actual-change
imputation model, although both use \code{xformla}.

\code{xformla} is evaluated on every F row, including frame-only, singleton and
zero-observation units; no common repeat-observer alias reduction is applied.
Each imputation pooling group and the shared counterfactual regression apply
their own QR alias/support preprocessing on their fitting rows with
prediction-row alias verification; the treatment logit keeps the raw columns
with the existing rank and separation refusals. Intercept-only and no-intercept
formulas are supported; empty formulas are refused. A rank-deficient pooling
group refuses; nothing is pooled automatically.

Sufficient (not necessary) conditions: actual-pair conditional mean trend-MAR
on F for every needed cell, positive cohort and cell masses, a correct linear
actual-change mean within each pooling group, and in addition a correct
synthetic treatment logit (\code{"ipw"}), a correct shared counterfactual
control-mean regression across all blocks (\code{"reg"}) or either of the two
(\code{"dr"}), with PT/no anticipation and the fixed-weight population law.
\code{"separate"} is the flexible default; \code{"within_cohort"} and
\code{"across_cohorts"} are explicit restrictions that can be misspecified, for
example under heterogeneous treated dynamics. No observation-probability model,
date-independence factorization or observed-level conditioning is used, and no
refusal arises merely because an unused observation logit would be separated
(for example on a balanced panel).

Complete original-unit influence functions on all F rows, including off-A and
frame-only units, carry the repeat-rate and imputation paths and the active
treatment/counterfactual nuisance paths into the existing GMM, aggregation, clustering,
multiplier-band and pretest consumers with core scaling
\eqn{\psi = (n_A/N)\phi_F}. \code{chained$selection} reports
\code{method = "imputation"}, both method axes, the pooling and the inactive
arguments; \code{chained$nuisance_diagnostics} reports the cells, pooling groups
(ranks, slopes, aliases), counterfactual regression, treatment logits, fitted
catalog and study-population support with zero observation-probability models;
\code{chained$outcome_design} records the raw full-F design. B/A fitting
samples, transport columns and custom IPW callbacks are unavailable; the
augmented correction that adds an observation-probability arm is described in
the next section. No numerical certificate or coverage claim is made.
}

\section{Augmented observation correction}{
With \code{selection = "conditional_mtar"} (or the legacy direct-pair alias
\code{"corollary1"}), \code{selection_est_method = "dr"} and
\code{est_method = "dr"}, every needed actual cohort/pair cell \eqn{c = (k, e)}
is reconstructed by
\deqn{\widehat Z_c = \widehat j_c(X) + S_e \{\Delta_e - \widehat j_c(X)\} / \widehat q_c(W),}{Zhat_c = jhat_c(X) + S_e (Delta_e - jhat_c(X)) / qhat_c(W),}
combining the actual-change imputation \eqn{\widehat j_c}{jhat_c} of the previous section
(observed pairs, fixed unit weights, \code{imputation_pooling}) with the fitted
observation probability \eqn{\widehat q_c}{qhat_c} of the conditional-MTAR section
(\code{selection_model}, \code{selection_pooling}, fitted on every F row of the
cohort). \eqn{\Delta_e} and \eqn{1/\widehat q_c}{1/qhat_c} are evaluated only on observed
pairs (\eqn{S_e = 1}); the residual term is not renormalized, \eqn{\widehat q}{qhat} enters
only through the floored inverse \eqn{1/\max(\widehat q, 10^{-6})}{1/max(qhat, 1e-6)}
(see Probability clipping), and no fill-only rule is used. The repeat rates \eqn{a_k}, the risk-set
treatment propensity logits on the unit masses \eqn{w a_k}, the shared
counterfactual control-mean regression (now fitted to \eqn{\widehat Z}{Zhat} over every
deduplicated block) and both normalizers are those of the pure-imputation DR
route. It is one estimator, not an average of the weighted and imputed estimates,
and it has no general finite equality with either.

Two formulas, two designs: \code{xformla} builds the imputation, counterfactual
regression and treatment designs, evaluated on every F row; \code{selection_xformla}
builds the observation-probability design only (\code{NULL} inherits
\code{xformla}), also on every F row. A variable that appears only in
\code{selection_xformla} never enters the outcome or treatment models. Known
cohorts and complete covariates of both formulas are required on every F row,
including frame-only, singleton and zero-observation units. The two pooling
arguments are independent; pair and date models, all three poolings of each
arm, intercept-only and no-intercept formulas are supported with the existing
alias and rank policies.

Every active fit must meet its existing regularity, even where a correctness
branch below would make it redundant: a cell whose pairs (pair model) or dates
(date model) are all observed or all unobserved, with its own free intercept,
is refused through the existing separation/boundary refusals, although pure
imputation would run on the same data. Imputation, treatment-logit and
counterfactual-regression failures also refuse. Nothing falls back to pure
imputation or to probability weighting.

Sufficient (not necessary) conditions for each primitive link: conditional mean
trend-MAR for every needed cell given the common conditioning set
\eqn{W} (the \code{xformla} covariates and any variables only in
\code{selection_xformla}), with positive observation probability; for each needed
cell EITHER the imputation arm (a correct linear actual-change mean given
\eqn{W}, represented by the \code{xformla} design, for every cell of its imputation
pooling group) OR the observation arm (a correct observation model given
\eqn{W} for every cell of its selection block; for \code{"independent_period_logit"}
correct date logits at both endpoint blocks and conditional independence of
the two used dates given cohort and \eqn{W}); the arms may differ across
cells; and in addition EITHER a correct synthetic treatment propensity logit OR
a correct shared counterfactual control-mean regression across all deduplicated
blocks (each of whose cells then needs one arm), with PT/no anticipation and the
fixed-weight population law. Pooling, date endpoints and the shared regression
couple cells, so a reconstruction error in another group, block or regression
block can move a shared coefficient. Correctness of an \code{xformla}-only
imputer alone does not establish reconstruction when variables only in
\code{selection_xformla} shift the outcome-change mean; a correct observation
model given \eqn{W} does. Adding selection covariates is therefore not
automatically harmless. When both arms of a needed cell are wrong, or both the
propensity and the regression are wrong, the link is generally biased.
Recovered joint ATT/GMM coordinates additionally need valid catalog moments for
every used link and the accepted identification premises; with incompatible
per-link pseudo-targets the limit depends on the weighting. There is no
finite-sample unbiasedness; small fitted probabilities inflate the variance.

Complete original-unit influence functions on all F rows carry the
observation-probability, imputation, repeat-rate, treatment-propensity and
shared-regression paths and both normalizers into the existing GMM, aggregation,
clustering, multiplier-band and pretest consumers with core scaling
\eqn{\psi = (n_A/N)\phi_F}. \code{chained$selection} reports
\code{method = "augmented"}, both method axes, \code{selection_model},
\code{selection_pooling}, \code{imputation_pooling}, the selection formula,
the assumption and inference text and every fit count. The
\code{selection_formula_inherited} field flags the \code{NULL} fallback;
\code{chained$nuisance_diagnostics} reports the cells (unit counts,
probability-weighted and discrepancy masses, imputed and augmented totals), the
imputation groups, the observation-probability blocks
(\code{records$selection$blocks}), the counterfactual regression on
\eqn{\widehat Z}{Zhat}, the treatment logits, the fitted catalog and study-population support;
\code{chained$outcome_design} records both designs. B/A fitting samples,
transport columns, custom IPW callbacks, SMAR/Corollary 2, singleton outcome
levels, cross-fitting/ML nuisances and efficiency claims are outside this route.
No numerical certificate or coverage claim is made.
}

\section{Probability clipping (chained mode)}{
Package-owned chained inverse-probability weights clip fitted probabilities
with the epsilon \eqn{10^{-6}}{1e-6} of did 2.5.1 and DRDID 1.3.0. Every risk-set
treatment propensity is capped above, \eqn{p_c = \min(p, 1 - 10^{-6})}{p_c = min(p, 1 - 1e-6)}, and the
consumed control odds is \eqn{p_c/(1 - p_c)}{p_c/(1 - p_c)} (at most about \eqn{10^6}{1e6}). The cap is
one-sided, and the upstream control trimming (\code{trim.level = 0.995}) is not
applied: no row, unit or comparison is dropped. Every inverse observation
denominator is floored at \eqn{10^{-6}}{1e-6}: a fitted direct pair probability
\eqn{p_U}{p_U}, the final product \eqn{p_s p_t}{p_s * p_t} of the two endpoint date probabilities
(the endpoints themselves are not floored) and a strictly positive empirical pair
observation rate. The multiplier \eqn{a_U/(\rho \max(p_U, 10^{-6}))}{a_U/(rho max(p_U, 1e-6))} keeps the
empirical repeat rate \eqn{a_U}{a_U} and the known \eqn{\rho}{rho}, so it is bounded by
\eqn{a_U/(\rho 10^{-6})}{a_U/(rho 1e-6)} rather than by \eqn{10^6}{1e6}. The augmented correction uses the
same floored inverse. Pure imputation has no observation denominator; its
treatment odds are capped. Outcome and imputation predictions are not
probabilities and are not clipped. No argument controls this; ordinary and
repeated cross-section mode delegate to did unchanged.

Every logit score and information matrix uses raw probabilities. Observation
weights in corrected treatment and outcome fits use the floored denominators;
the raw empirical rates and their influence functions are kept. Rank, QR, separation, convergence, boundary and
nonfinite refusals remain, as do zero-support, invalid-rate and nonpositive
normalizer refusals; clipping does not rescue a failed fit. Finite predictions
replayed on other rows are evaluated stably and clipped rather than overflowing.
Influence functions differentiate the clipped weights: a capped odds or a floored
denominator has no derivative in its probability coefficients away from the
cutoff. A value exactly at a cutoff is a kink; the strict activity convention
assigns it a zero probability derivative, not an ordinary smooth derivative. The
repeat-rate numerator, every normalizer and all other channels are retained.

Ordinary influence-function inference for the clipped estimator also needs a
negligible threshold-crossing remainder, for example no probability mass of the
fitted probabilities at either cutoff. When clipping binds with positive
population probability, the estimand is the clipped functional, which can differ
from the original ATT target and can change the usual observation or treatment
robustness branches. Original-target claims require clipping to be inactive in
population or additional conditions preserving that target; no coverage of the
original target follows from software checks.
The GMM argument that omits estimated-weight derivatives assumes compatible
comparison moments; clipped moments that are incompatible remove that premise.
Every fitted treatment model, selection-weighted outcome-regression support and
augmented cell records a compact summary of its raw and clipped probability or
rate ranges and the numbers of capped or floored rows among the route diagnostics.
Treatment summaries include the largest fitting weight; selection-weighted
outcome-regression supports include the largest WLS weight, and augmented cells
include the largest observed correction weight \eqn{w a_U/\max(q, 10^{-6})}{w a_U/max(q, 1e-6)}.
}
