\name{aggte}
\alias{aggte}
\title{Aggregate group-time average treatment effects}
\description{
\code{aggte()} has the signature of \code{did::aggte()} (did 2.5.1). Results
from ordinary mode are aggregated by \code{did::aggte()} unchanged, with only
the arguments you supply forwarded. Chained results (\code{cdid_MP} objects
from \code{\link{att_gt}} with \code{est_method = "dr"}, \code{est_method = "ipw"},
\code{est_method = "reg"}, or a custom IPW callback) are aggregated from
their retained state without refitting the estimator, with did's selector and
inheritance conventions.
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
\item{MP, type, balance_e, min_e, max_e, na.rm, bstrap, biters, cband, alp,
  clustervars}{As in \code{did::aggte()}: an \code{MP} object from
  \code{\link{att_gt}}, the aggregation type (\code{"simple"},
  \code{"group"}, \code{"dynamic"} or \code{"calendar"}), event-time
  restrictions, missing-value handling and inference settings. For chained
  results, \code{NULL} inherits the estimation request's \code{bstrap},
  \code{biters}, \code{cband}, \code{alp} and \code{clustervars}; explicit
  values override them with their released did 2.5.1 meaning (including
  \code{bstrap = FALSE, cband = TRUE}: analytic standard errors with a
  multiplier band for group/dynamic/calendar bins and did's warning; simple
  aggregation has no bin band). Cluster overrides that cannot be honoured are
  refused instead of silently changed. \code{na.rm = TRUE} requests the
  explicitly labelled restricted target over the supported cells; with
  \code{na.rm = FALSE} an aggregate needing an unavailable cell is \code{NA}
  with the missing cells named. As in did, a result fitted with
  \code{compute_inffunc = FALSE} has no influence functions and is refused
  (\code{"cdid_chained_request_no_inffunc_error"}); multiplier draws of a
  chained result reuse the workers admitted by its \code{pl}/\code{cores}
  request, with identical results.}
}
\value{For ordinary results, the \code{AGGTEobj} returned by
\code{did::aggte()}; its \code{call} element records the delegated call. For
chained results, a \code{cdid_AGGTEobj} (inheriting \code{AGGTEobj}) with the
population cohort-share aggregation (first-calendar cohort masses for
reference-date or changing-weight fits on every chained route; see \code{\link{att_gt}}), complete
full-frame influence functions including the share uncertainty, the requested interval law, the record
that complete per-fit error certification was not requested or performed, and
the intended/used/excluded cell records. Multiplier draws use the caller's RNG stream
(\code{set.seed()} before the call reproduces them); no fixed seed exists. Use
\code{summary()}, \code{\link{tidy}}() (columns \code{label}, \code{available},
\code{restricted}, \code{n.intended}, \code{n.used}, \code{excluded.ids},
\code{reason}), \code{\link{glance}}(), \code{overall.att}/\code{overall.se} and
\code{\link{ggdid}}() (for \code{"group"}, \code{"dynamic"} and
\code{"calendar"} types; \code{"simple"} has no plot) to inspect aggregated
results.}
\seealso{\code{\link{att_gt}}, \code{\link{conditional_did_pretest}}, \code{\link{tidy}}, \code{\link{glance}}, \code{\link{ggdid}}, \code{vignette("cdid", package = "cdid")}}
\examples{
set.seed(1)
out <- cdid::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                    gname = "first.treat", data = did::mpdta,
                    mode = "ordinary", bstrap = FALSE)
cdid::aggte(out, type = "simple")
cdid::aggte(out, type = "group")
}
