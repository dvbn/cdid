# Stage beside the existing packaged tests; run against the patched installation.
test_that("compact raw beta0 retains supplied bits with unchanged arithmetic", {
  if (!requireNamespace("cdid", quietly = TRUE))
    stop("Actual installed cdid required; no source fallback.", call. = FALSE)
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns))
  fn <- function(name) {
    f <- get(name, envir = ns, inherits = FALSE)
    stopifnot(is.function(f), identical(environment(f), ns))
    f
  }
  # Same panel24 design and explicit diagnostic controls as point-view tests.
  panel <- utils::read.csv(testthat::test_path("fixtures", "fitted_panel24.csv"),
                           stringsAsFactors = FALSE)
  units <- panel[!duplicated(panel$id), ]
  Z <- cbind(one = 1, sign = as.double(units$sign))
  oc <- fn("ipw_owner_control")(8, 100000, 500000, 5000000, 100000)
  cc <- fn("logit_certificate_control")(18, 40, 16, 8, 500000, 100000000)
  nc <- fn("ipw_nuisance_problem_control")(8, 100000, 500000, 100000000)
  po <- fn("construct_ipw_panel_owner")(panel, "id", "time", "y", "g", control = oc)
  fo <- fn("construct_ipw_frame_owner")(po, units$id, Z, rep(1, nrow(units)), control = oc)
  nv <- fn("construct_ipw_nuisance_view")(fo, 3, Inf, 1, 2, control = oc)
  beta <- c(one = -as.double(0L), sign = as.double(0L))
  positive <- as.double(beta); positive[positive == 0L] <- as.double(0L)
  F <- diag(2L)
  make <- function(b) fn("construct_certified_ipw_nuisance_problem_view")(
    nv, b, F, 1e-6, cc, nc, oc)
  negative_problem <- make(beta)
  positive_problem <- make(positive)
  expect_false(identical(as.double(beta), positive, num.eq = FALSE))
  expect_true(identical(negative_problem$raw$beta0, as.double(beta), num.eq = FALSE))
  expect_true(identical(positive_problem$raw$beta0, positive, num.eq = FALSE))
  expect_true(identical(negative_problem$certificate$provenance$beta0, positive,
                        num.eq = FALSE))
  expect_true(fn("ipw_point_view_identical")(negative_problem$raw$beta0, as.double(beta)))
  expect_false(fn("ipw_point_view_identical")(negative_problem$raw$beta0, positive))
  for (field in c("certificate", "nuisance_enclosures", "arithmetic_refusal",
                  "controls", "preflight", "work", "status"))
    expect_true(identical(negative_problem[[field]], positive_problem[[field]], num.eq = FALSE),
                info = field)
  for (problem in list(negative_problem, positive_problem)) {
    fn("validate_certified_ipw_nuisance_problem")(problem, require_certified = FALSE)
    expect_true(identical(problem$frame_ref, fo$ref, num.eq = FALSE))
    expect_true(identical(problem$panel_ref, po$ref, num.eq = FALSE))
    expect_true(identical(problem$nuisance_ref, nv$ref, num.eq = FALSE))
  }
})
