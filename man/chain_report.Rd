\name{chain_report}
\alias{chain_report}
\alias{print.cdid_chain_report}
\alias{as.data.frame.cdid_chain_report}
\title{Inspect the comparisons used in a chained fit}
\description{
See the date pairs and control groups that identify a cohort-date effect and
the weights used to combine those comparisons. Supply a chained cohort fit
from \code{\link{att_gt}()}.
}
\usage{
chain_report(
  object,
  cohort = NULL,
  time = NULL,
  include_pre = FALSE,
  tol = 1e-10,
  group = NULL,
  t = NULL
)

\method{print}{cdid_chain_report}(
  x,
  digits = 4,
  show_gmm_weights = FALSE,
  max_cells = 30L,
  ...
)

\method{as.data.frame}{cdid_chain_report}(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("summary", "witness", "gmm_weights"),
  ...
)
}
\arguments{
\item{object}{A \code{cdid_MP} result produced by \code{\link{att_gt}} in
  chained mode.}

\item{cohort, group}{Optional numeric vector of treatment cohorts to keep.
  Supply at most one of \code{cohort} or \code{group}.}

\item{time, t}{Optional numeric vector of calendar periods to keep. Supply
  at most one of \code{time} or \code{t}.}

\item{include_pre}{Logical scalar (\code{FALSE} by default). When
  \code{FALSE} and no period filter is given, only post-treatment cells are
  included.}

\item{tol}{Nonnegative numeric threshold (\code{1e-10} by default). GMM
  moment weights at or below it in absolute value are omitted from the active
  weight table.}

\item{x}{A \code{cdid_chain_report} object.}

\item{digits}{Number of decimal places for printed estimates and weights.}

\item{show_gmm_weights}{Logical. If \code{TRUE}, the table of nonzero GMM
  moment weights is printed beneath each cell.}

\item{max_cells}{Maximum number of cells printed before truncating.}

\item{row.names, optional}{As in \code{as.data.frame()}.}

\item{what}{Which table to extract: \code{"summary"} (one row per cell),
  \code{"witness"} (the identification path of each cell: its date pairs and
  control pools in order), or \code{"gmm_weights"} (all nonzero GMM weights,
  whose \code{contribution} column sums to the reported estimate).}

\item{...}{Additional arguments (unused).}
}
\details{
A cell identified in one step compares the cohort's reference period with
the target date directly. A cell that needs chaining passes through
intermediate dates, each step using the units observed at both dates of
that step. The \code{"witness"} table shows one such path. The
\code{"gmm_weights"} table shows how the fit actually combined all retained
moments, which generally involves more comparisons than the single path.
Weights reflect the \code{gmm_weighting} choice recorded in the fit.
}
\value{
\code{chain_report()} returns an object of class \code{"cdid_chain_report"}
with elements \code{summary} (one row per cell), \code{chains} (a named list
holding the \code{witness} and \code{gmm_weights} data frames of each cell)
and \code{meta}. \code{as.data.frame()} returns the requested table as a
\code{data.frame}.
}
\seealso{\code{\link{att_gt}}, \code{\link{nuisance_models}}}
\examples{
set.seed(1)
rows <- list()
uid <- 0L
for (g in 2005:2008) {
  for (i in 1:20) {
    uid <- uid + 1L
    a_i <- rnorm(1)
    for (tv in max(2002L, g - 3L):2008L) {
      rows[[length(rows) + 1L]] <- data.frame(
        id = uid, period = tv, G = g, y = a_i + (tv >= g) * 2 + rnorm(1, sd = 0.2)
      )
    }
  }
}
fit <- att_gt(
  yname = "y", tname = "period", idname = "id", gname = "G",
  data = do.call(rbind, rows),
  control_group = "notyettreated", base_period = "universal",
  cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
  bstrap = FALSE, cband = FALSE
)
cr <- chain_report(fit, cohort = 2005, time = 2007)
print(cr, show_gmm_weights = TRUE)
as.data.frame(cr, what = "gmm_weights")
}
