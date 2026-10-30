\name{reexports}
\alias{reexports}
\alias{ggdid}
\alias{tidy}
\alias{glance}
\title{Functions re-exported from did}
\description{
These did functions are re-exported unchanged so that the usual workflow
works with \code{cdid::} qualification. See their documentation in the did
package: \code{did::ggdid}, \code{did::tidy} and \code{did::glance} (the
\code{generics} methods for \code{MP} and \code{AGGTEobj} objects). For chained
results (\code{cdid_MP} and \code{cdid_AGGTEobj}), \pkg{cdid} registers S3
methods for \code{print()}, \code{summary()}, \code{tidy()}, \code{glance()} and
\code{ggdid()} that extend upstream output with chained support and inference
metadata. \code{conditional_did_pretest()} is cdid's own wrapper (see
\code{\link{conditional_did_pretest}}): chained by default for panels, and
\pkg{did}'s \code{conditional_did_pretest()} unchanged in ordinary mode.
}
\details{
For chained results:
\itemize{
  \item \code{print(x)} calls \code{summary(x)}. Both print the estimate table
    with chained support, parallel-trends and inference information, then
    return the original object invisibly. For \code{cdid_MP}, the summary also
    reports unit counts, observation correction, GMM weighting and pretrend
    test status.
  \item \code{tidy(x)} returns one row per reported cell (\code{cdid_MP}) or
    aggregate bin (\code{cdid_AGGTEobj}); \code{type = "group"} also includes
    an overall row. It keeps \pkg{did}'s columns
    (\code{estimate}, \code{std.error}, \code{conf.low}, \code{conf.high},
    \code{point.conf.low}, \code{point.conf.high}), updating pointwise
    \code{statistic} and \code{p.value}, and adding \code{target.id},
    \code{role} (\code{"pre"}, \code{"post"}, \code{"reference"}),
    \code{available}, \code{support.status}, \code{support.reason},
    \code{certified} and \code{band.family} for \code{cdid_MP}, or \code{label},
    \code{available}, \code{restricted}, \code{n.intended}, \code{n.used},
    \code{excluded.ids}, \code{reason}, \code{certified} and \code{band.family}
    for \code{cdid_AGGTEobj}.
  \item \code{glance(x)} returns a one-row data frame adding \code{estimator},
    \code{n.requested}, \code{n.available}, \code{band.family}, \code{inference}
    and \code{clustered} for \code{cdid_MP}, or \code{estimator}, \code{status},
    \code{na.rm}, \code{n.bins}, \code{n.bins.available}, \code{overall.status},
    \code{band.family} and \code{inference} for \code{cdid_AGGTEobj}.
  \item \code{ggdid(x, ...)} returns a \pkg{ggplot2} plot of \code{cdid_MP}
    cells or \code{cdid_AGGTEobj} \code{"group"}, \code{"dynamic"} and
    \code{"calendar"} aggregates (\code{"simple"} has no plot) and attaches the
    plotted grid or bin status table in \code{attr(p, "cdid_presentation_grid")}
    or \code{attr(p, "cdid_presentation_bins")}.
}
}
\seealso{\code{\link{att_gt}}, \code{\link{aggte}}, \code{\link{conditional_did_pretest}}, \code{vignette("cdid", package = "cdid")}}
