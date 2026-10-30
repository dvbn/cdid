\docType{package}
\name{cdid-package}
\alias{cdid}
\alias{cdid-package}
\title{cdid: The Chained Difference-in-Differences}
\description{
cdid 1.1.0.0, an update of the package after version 0.1.1, is a small layer on
the did package. It keeps did's \code{att_gt()} and \code{aggte()} interface and
adds an estimation mode for panel data with missing outcomes.

\itemize{
  \item \code{mode = "ordinary"} delegates unchanged to \code{did::att_gt()}
    and \code{did::aggte()}. This is also the default for repeated cross
    sections.
  \item Chained estimation (Bellego, Benatia and Dortet-Bernadet, 2025)
    is the default for panel data. Built-in doubly robust \code{est_method="dr"}
    (the default), \code{est_method="ipw"}, and control-mean outcome regression
    \code{est_method="reg"} are supported, with baseline (\code{selection = "none"})
    or conditional-MTAR (\code{"conditional_mtar"} / legacy \code{"corollary1"})
    missing-outcome handling via observation weighting (\code{selection_est_method = "ipw"}),
    change imputation (\code{"imputation"}) or the augmented correction (\code{"dr"}).
    Custom panel IPW callbacks support baseline selection and, when they return
    \code{log_weight_gradient}, conditional-MTAR observation weighting
    (\code{selection_est_method = "ipw"}).
}

cdid requires did 2.5.1 or later and DRDID 1.3.0 or later, and has been
tested with did 2.5.1 and DRDID 1.3.0. Later versions are accepted without a
guarantee that everything works with them or that results are unchanged. The
optional pinned setup in \file{validation/SETUP.md} of the source repository
installs exactly the tested versions and verifies their library paths. Use
\code{citation("cdid")} for references.
}
\details{
Quick navigation for offline help:
\itemize{
  \item \code{\link{att_gt}}: estimate group-time average treatment effects from
    long panel data (\code{data}) and an optional study-population unit table
    (\code{frame}). Returns a \code{cdid_MP} object in chained mode.
  \item \code{\link{aggte}}: aggregate a fitted result into \code{"simple"},
    \code{"group"}, \code{"dynamic"} or \code{"calendar"} summaries without
    refitting (\code{cdid_AGGTEobj} in chained mode).
  \item \code{\link{conditional_did_pretest}}: run the chained conditional
    pretrend test on panel data (or \pkg{did}'s \code{conditional_did_pretest()} in
    ordinary mode).
  \item \code{\link{reexports}}: \code{\link{tidy}}, \code{\link{glance}} and
    \code{\link{ggdid}} generics re-exported from \pkg{did}, with chained S3
    methods that report cell and bin availability and support reasons.
}

Run \code{vignette("cdid", package = "cdid")} for a self-contained introductory
walkthrough. Extended documentation on panel data preparation, investigating
missingness and observation models is available at
\url{https://github.com/dvbn/cdid}.
}
\seealso{\code{\link{att_gt}}, \code{\link{aggte}}, \code{\link{conditional_did_pretest}}, \code{\link{reexports}}, \code{vignette("cdid", package = "cdid")}}
