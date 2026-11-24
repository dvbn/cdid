# The copy into the returned object of a static OR/DR fit without memory phases is charged before estimation,
# with the same check and count as the charge made after estimation.
copy_panel <- function() {
  set.seed(41)
  n <- 80L
  d <- expand.grid(id = seq_len(n), t = 1:6)
  d$g <- rep(c(0, 3, 5), length.out = n)[d$id]
  d$x <- rnorm(n)[d$id]
  d$y <- d$x + .3 * d$t + rnorm(nrow(d)) + (d$g > 0 & d$t >= d$g)
  d$y[runif(nrow(d)) < .15] <- NA_real_
  d
}

copy_fit <- function(m, mib) {
  suppressMessages(att_gt("y", "t", "id", "g", data = m, xformla = ~x, est_method = "reg",
    outcome_pooling = "separate", control_group = "nevertreated", gmm_weighting = "identity",
    bstrap = FALSE, cband = FALSE, control = list(numeric_workspace_bytes = mib * 2^20)))
}

test_that("the copy into the returned object is refused before estimation with the late count", {
  m <- copy_panel()
  ns <- asNamespace("cdid")
  entered <- new.env()
  entered$n <- 0L
  count <- function() entered$n <- entered$n + 1L
  suppressMessages(trace("cp04_stream_covariance", where = ns, print = FALSE, tracer = bquote(.(count)())))
  on.exit(suppressMessages(untrace("cp04_stream_covariance", where = ns)), add = TRUE)
  e <- tryCatch(copy_fit(m, 10), error = identity)
  expect_s3_class(e, "cdid_ipw_owner_bound_error")
  expect_identical(e$resource$check, "OR public diagnostic copy concurrent payload")
  expect_identical(e$resource$requested, 2013711)
  expect_identical(entered$n, 0L)
  expect_match(conditionMessage(e), "Increase the workspace limit.", fixed = TRUE)
  entered$n <- 0L
  f <- copy_fit(m, 16)
  expect_true(length(f$att) > 0L)
  expect_gt(entered$n, 0L)
})

copy_fit_dr_within <- function(m, mib) {
  suppressMessages(att_gt("y", "t", "id", "g", data = m, xformla = ~x, est_method = "dr",
    outcome_pooling = "separate", control_group = "nevertreated", gmm_weighting = "within",
    bstrap = FALSE, cband = FALSE, control = list(numeric_workspace_bytes = mib * 2^20)))
}

test_that("a memory-phase fit is not charged by the precheck", {
  # The copy of this DR `within` fit is counted from the returned result after estimation. The precheck's
  # formula would count 3,424,895 entries, more than the 25 MiB limit, while the copy needs about 3.16 million.
  m <- copy_panel()
  f <- copy_fit_dr_within(m, 25)
  expect_true(length(f$att) > 0L)
  expect_false(is.null(f$chained$resource_accounting$memory_phases))
})

test_that("the copy is still checked after estimation where the precheck does not run", {
  m <- copy_panel()
  ns <- asNamespace("cdid")
  entered <- new.env()
  entered$n <- 0L
  count <- function() entered$n <- entered$n + 1L
  suppressMessages(trace("cp04_stream_covariance", where = ns, print = FALSE, tracer = bquote(.(count)())))
  on.exit(suppressMessages(untrace("cp04_stream_covariance", where = ns)), add = TRUE)
  e <- tryCatch(copy_fit_dr_within(m, 24), error = identity)
  expect_s3_class(e, "cdid_ipw_owner_bound_error")
  expect_identical(e$resource$check, "DR public diagnostic copy concurrent payload")
  expect_gt(entered$n, 0L)
})
