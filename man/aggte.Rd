\name{aggte}
\alias{aggte}
\title{Aggregate group-time average treatment effects}
\description{
Summarize fitted cohort-date effects into an overall effect, cohort averages,
an event study or calendar-date averages. With a switching model, use the
switcher-weighted estimate already stored in \code{fit$aggregate}.
}
\usage{
aggte(
  MP,
  type = "group",
  balance_e = NULL,
  min_e = -Inf,
  max_e = Inf,
  na.rm = FALSE,
  bstrap = NULL,
  biters = NULL,
  cband = NULL,
  alp = NULL,
  clustervars = NULL
)
}
\arguments{
\item{MP}{An \code{MP} or \code{cdid_MP} object returned by
  \code{\link{att_gt}}.}

\item{type}{Use \code{"simple"} for an overall average, \code{"group"} for
  cohort averages, \code{"dynamic"} for event time, or \code{"calendar"} for
  calendar-time averages.}

\item{balance_e, min_e, max_e}{Event-time selectors with did's meaning.
  Simple and group aggregation use only \code{max_e}. Dynamic aggregation
  uses event time \code{t - g} in original date units with did's window and
  cohort-balance rules. Calendar aggregation ignores all three with did's
  warning. These selectors change the reported target, not the fitted
  comparisons.}

\item{na.rm}{\code{FALSE} (default): an aggregate that needs an unavailable
  cell is \code{NA}, and the missing cells are named. \code{TRUE}: the
  aggregate is computed over the available cells and labelled as a restricted
  target.}

\item{bstrap, biters, cband, alp, clustervars}{Inference settings. For
  chained results, \code{NULL} inherits the value used by \code{att_gt()}.
  An explicit value overrides it with did's meaning. In particular
  \code{bstrap = FALSE, cband = TRUE} gives analytic standard errors with a
  multiplier band for group, dynamic and calendar bins, together with did's
  warning. Simple aggregation has no band. A cluster override that cannot be
  honored is refused rather than silently changed. For chained results the
  multiplier signs are drawn in batches of the fit's processing chunk size,
  so reproducing the draws requires the same package version, seed and RNG
  kind, the same data in the same row order, and the same
  \code{control$numeric_workspace_bytes} (and \code{control$chunk_rows}, if
  set) in the fit.}
}
\details{
For chained results, aggregation reuses the influence functions and cohort
shares stored by \code{att_gt()}. Cohorts are weighted by their population
shares in the study population. For fits with reference-date covariates or
weights that change within a unit, cohorts are weighted by their masses at
the first stored date, and the uncertainty of those shares enters the
standard errors. With \code{fix_weights = "base_period"} each cell keeps the
weights of its own cohort's comparisons, so an aggregate is a declared mixture
of cohort effects rather than the effect in one weighted population. See
\code{\link{att_gt}} for the reference-date rules.

Multiplier draws use your session's random number stream. Calling
\code{set.seed()} before \code{aggte()} reproduces them. There is no built-in
seed. Draws reuse the worker count requested by the fit's \code{pl} and
\code{cores} request and give the same results as a serial run.

A result fitted with \code{compute_inffunc = FALSE} has no influence
functions. Refit with \code{compute_inffunc = TRUE} to obtain summaries.
}
\value{
For ordinary results, the \code{AGGTEobj} returned by \code{did::aggte()}.
Its \code{call} element records the delegated call.

For chained results, a \code{cdid_AGGTEobj} (inheriting \code{AGGTEobj})
with the overall estimate and its standard error (\code{overall.att},
\code{overall.se}), the bin estimates, full-frame influence functions that
include the cohort-share uncertainty, the requested interval type, and the
record of intended, used and excluded cells. Inspect it with
\code{summary()}, \code{\link{tidy}}() (columns \code{label},
\code{available}, \code{restricted}, \code{n.intended}, \code{n.used},
\code{excluded.ids}, \code{reason}), \code{\link{glance}}() and
\code{\link{ggdid}}(). Plots exist for \code{"group"}, \code{"dynamic"} and
\code{"calendar"} aggregates. \code{"simple"} has no plot.
}
\seealso{\code{\link{att_gt}}, \code{\link{conditional_did_pretest}}, \code{\link{tidy}}, \code{\link{glance}}, \code{\link{ggdid}}, \code{vignette("cdid", package = "cdid")}}
\examples{
set.seed(1)
out <- cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                    gname = "first.treat", data = did::mpdta,
                    mode = "ordinary", bstrap = FALSE, cband = FALSE)
cdid::aggte(out, type = "simple")
cdid::aggte(out, type = "group")
}
