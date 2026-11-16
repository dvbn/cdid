\name{nuisance_models}
\alias{nuisance_models}
\title{Inspect fitted treatment and missing-outcome models}
\description{
Show the fitted equations, coefficient estimates and standard errors for the
models used by \code{\link{att_gt}()}. Set
\code{control = list(keep_nuisance = TRUE)} when fitting to retain the model
information. The report uses the fitted models and does not fit them again.
}
\usage{
nuisance_models(fit, component = NULL, block = NULL)
}
\arguments{
\item{fit}{A result of \code{\link{att_gt}()}.}

\item{component}{\code{NULL} (a compact index of all components) or one of
  \code{"selection"}, \code{"outcome"}, \code{"treatment"},
  \code{"imputation"}.}

\item{block}{\code{NULL}, or one fitted block of the component: its stable
  identifier (\code{"S1"}, \code{"S2"}, \ldots, in block order), its integer
  position or its label (for example \code{"within_cohort[g=2004]"}).
  Supplying \code{block} alone selects \code{component = "selection"}.}
}
\details{
Call \code{nuisance_models(fit)} to see which models are available.
Select a component, then a block when several equations were fitted.
For example, a model with separate slopes can have one block for each
date pair and cohort.

\describe{
  \item{\code{"selection"}}{Logit for observing a date or both outcomes
    in a date pair. With date-specific observation probabilities, the pair probability is the
    product of the two fitted date probabilities.}
  \item{\code{"imputation"}}{Regression of observed actual changes.
    The estimates are used to predict missing changes within the relevant
    cohort or switching group.}
  \item{\code{"treatment"}}{Logit for treatment status within the units
    eligible for the comparison.}
  \item{\code{"outcome"}}{Regression for the counterfactual outcome change.
    In the linear cohort-by-date model, the report gives coefficients
    relative to the stated reference levels.}
}

The estimation settings determine which models are fitted. A doubly robust
treatment fit with augmented observation correction can report all four
models. Some specifications, including custom estimation functions and some
models with changing covariates or weights, give a reason when a report is
unavailable.

Coefficient standard errors use unit contributions from the fitted model,
combined within clusters when clustering was requested. The table can have
missing standard errors for a point-only fit, insufficient cluster support,
or insufficient numerical workspace. See the reported reason.
}
\value{
With \code{component = NULL}, an index with model availability and fitted
blocks. With a component, an object with fitted-block and coefficient tables.
The coefficient table includes \code{term}, \code{estimate},
\code{std_error}, \code{z_value} and \code{p_value}. Selecting a block returns
its equation, fitting sample, weights and coefficient table. Use
\code{print()} to display the result. A component without a report includes
its status and reason.
}
\seealso{\code{\link{att_gt}}, \code{\link{chain_report}}}
\examples{
\dontrun{
fit <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g",
              xformla = ~ x + region, data = panel, weightsname = "w",
              est_method = "dr", selection = "conditional_mtar",
              selection_model = "independent_period_logit",
              selection_pooling = "within_cohort",
              control = list(keep_nuisance = TRUE))
nuisance_models(fit)                       # index and template equation
nuisance_models(fit, "selection")          # every block's table
nuisance_models(fit, block = "S1")         # one block
nuisance_models(fit, "outcome")            # status for this specification

# Augmented correction with treatment DR: every fitted role.
fit_aug <- att_gt(yname = "y", tname = "t", idname = "id", gname = "g",
                  xformla = ~ x + region, data = panel, weightsname = "w",
                  est_method = "dr", selection = "conditional_mtar",
                  selection_est_method = "dr", outcome_pooling = "pooled",
                  control = list(keep_nuisance = TRUE))
nuisance_models(fit_aug, "selection", block = "S1")
nuisance_models(fit_aug, "outcome")              # pooled fit O1
nuisance_models(fit_aug, "treatment", block = "T1")
nuisance_models(fit_aug, "imputation")
}
}
