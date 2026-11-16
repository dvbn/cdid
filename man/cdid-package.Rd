\docType{package}
\name{cdid-package}
\alias{cdid}
\alias{cdid-package}
\title{cdid: The Chained Difference-in-Differences}
\description{
Estimate difference-in-differences effects with panels in which some outcomes
are missing. Chained estimation combines changes measured within units
across the dates when they are observed. Use \code{mode = "ordinary"} for
\pkg{did}'s estimator or for repeated cross sections.
}
\details{
Use \code{gname} to estimate effects by first-treatment cohort and date
when treatment is absorbing. Use \code{dname} for effects at entry and exit
when a binary treatment can change and past treatment has no carryover.

Choose the treatment adjustment with \code{est_method}. If the units
observed at both dates represent their cohort with respect to mean changes
only after conditioning on covariates, use
\code{selection = "conditional_mtar"}. Choose observation weighting,
imputation of actual changes or their augmented combination with
\code{selection_est_method}.

\itemize{
  \item \code{\link{att_gt}()} fits the model from a long data frame and,
    for cohort models, an optional study-unit table.
  \item \code{\link{aggte}()} summarizes cohort effects. A switching fit
    already contains the switcher-weighted estimate in \code{fit$aggregate}.
  \item \code{\link{chain_report}()} shows which comparisons identify
    cohort effects and how they are weighted.
  \item \code{\link{nuisance_models}()} displays fitted equations and
    coefficient tables when requested during estimation.
  \item \code{\link{conditional_did_pretest}()} tests pre-treatment trends
    conditional on the model covariates.
  \item \code{\link{tidy}()}, \code{\link{glance}()} and \code{\link{ggdid}()}
    return tables, summaries and plots.
}

Open \code{vignette("cdid", package = "cdid")} for a hands-on introduction.
The \href{https://github.com/dvbn/cdid/tree/main/docs}{user guide}
covers preparing data, choosing a model and inspecting results.
Use \code{citation("cdid")} for references.
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
(revised 2026-07-22). Static binary specialization, used on the \code{dname}
route.
}
\seealso{\code{\link{att_gt}}, \code{\link{aggte}}, \code{\link{chain_report}}, \code{\link{conditional_did_pretest}}, \code{\link{nuisance_models}}, \code{\link{reexports}}, \code{vignette("cdid", package = "cdid")}}
