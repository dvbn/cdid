\name{reexports}
\alias{reexports}
\alias{ggdid}
\alias{tidy}
\alias{glance}
\title{Functions re-exported from did}
\description{
Display results from \code{\link{att_gt}()} and \code{\link{aggte}()}.
These functions use the same interface as in \pkg{did}.
}
\usage{
ggdid(object, ...)
tidy(x, ...)
glance(x, ...)
}
\arguments{
  \item{object, x}{A result returned by \code{\link{att_gt}()} or
    \code{\link{aggte}()}.}
  \item{...}{Arguments passed to the method for the supplied result.}
}
\details{
\describe{
  \item{\code{tidy(x)}}{Returns a table of effects. Columns include
    \code{estimate}, \code{std.error}, \code{conf.low} and \code{conf.high}.
    For chained fits, \code{available} and the support or aggregation
    columns explain which effects the data support.}
  \item{\code{glance(x)}}{Returns a one-row summary of the fit, with
    information such as the estimation method, the available effects and
    the type of inference.}
  \item{\code{ggdid(object, ...)}}{Plots the effects from a cohort or switching
    fit or a cohort, event-time or calendar-time summary. Use
    \code{tidy(x)} to view a simple aggregate as a table.}
}
Use \code{summary(x)} for a printed table and a description of the fit.
}
\seealso{\code{\link{att_gt}}, \code{\link{aggte}}, \code{\link{conditional_did_pretest}}, \code{vignette("cdid", package = "cdid")}}
