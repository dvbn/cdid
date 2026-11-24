# Workspace accounting changes must leave public results unchanged. Fixture: the seeded 80-unit, six-date panel
# with cohorts 3, 5 and never and 15% missing outcomes.
parity_panel <- function() {
  set.seed(41)
  n <- 80L
  d <- expand.grid(id = seq_len(n), t = 1:6)
  d$g <- rep(c(0, 3, 5), length.out = n)[d$id]
  d$x <- rnorm(n)[d$id]
  d$y <- d$x + .3 * d$t + rnorm(nrow(d)) + (d$g > 0 & d$t >= d$g)
  m <- d
  m$y[runif(nrow(m)) < .15] <- NA_real_
  m$cl <- (m$id - 1L) %/% 4L + 1L
  m
}

test_that("workspace limit changes multiplier batches but not estimates or the final RNG state", {
  m <- parity_panel()
  fit_at <- function(mib) {
    set.seed(101)
    f <- suppressMessages(att_gt("y", "t", "id", "g", data = m, xformla = ~x, est_method = "dr",
      outcome_pooling = "separate", control_group = "nevertreated", bstrap = TRUE, cband = TRUE,
      biters = 199L, gmm_weighting = "within", control = list(numeric_workspace_bytes = mib * 2^20)))
    list(fit = f, rng = .Random.seed)
  }
  a <- fit_at(32)
  b <- fit_at(64)
  expect_identical(a$fit$chained$inference$rng$batch_rows, 22L)
  expect_identical(b$fit$chained$inference$rng$batch_rows, 44L)
  expect_identical(a$fit$att, b$fit$att)
  expect_identical(as.matrix(a$fit$inffunc), as.matrix(b$fit$inffunc))
  expect_false(identical(a$fit$se, b$fit$se))
  expect_identical(a$rng, b$rng)
})

expect_parity_diagnostic_shape <- function(actual, reference, info) {
  expect_identical(typeof(actual), typeof(reference), info = info)
  expect_identical(length(actual), length(reference), info = info)
  expect_identical(attributes(actual), attributes(reference), info = info)
  expect_identical(is.finite(actual), is.finite(reference), info = info)
  expect_identical(actual[!is.finite(actual)], reference[!is.finite(reference)], info = info)
  expect_identical(actual == 0, reference == 0, info = info)
}

expect_parity_covariance_handling <- function(actual, reference, key) {
  fields <- list(
    c("scale", "trace_B"),
    c("rho", "rho_raw", "d_min", "d_max", "eigenvalues_B",
      "norm_B_frobenius", "b_roundoff_bound")
  )
  observed <- list(actual, actual$conditioning)
  archived <- list(reference, reference$conditioning)
  for (index in seq_along(fields)) {
    for (name in fields[[index]]) {
      info <- paste(key, name)
      x <- observed[[index]][[name]]
      y <- archived[[index]][[name]]
      expect_parity_diagnostic_shape(x, y, info)
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
      expect_true(all(abs(x - y) <= 1e-12 * pmax(1, abs(x), abs(y))), info = info)
      observed[[index]][[name]] <- archived[[index]][[name]]
    }
  }

  conditioning <- actual$conditioning
  residual <- conditioning$r_eigen_residual
  original <- conditioning$eigensystem_refinement$original_residual
  expect_parity_diagnostic_shape(residual, reference$conditioning$r_eigen_residual, key)
  expect_parity_diagnostic_shape(
    original, reference$conditioning$eigensystem_refinement$original_residual, key
  )
  # These fixtures require no eigenvector refinement. The same residual is
  # recorded twice and must satisfy the solver's existing reconstruction bound.
  expect_identical(conditioning$eigensystem_refinement$method, "none", info = key)
  expect_identical(original, residual, info = key)
  bound <- 100.0 * as.double(actual$dimensions$L_retained) *
    .Machine$double.eps * conditioning$norm_B_frobenius
  expect_identical(conditioning$b_roundoff_bound, bound, info = key)
  expect_true(is.finite(residual) && residual >= 0 && residual <= bound, info = key)
  observed[[2]]$r_eigen_residual <- reference$conditioning$r_eigen_residual
  observed[[2]]$eigensystem_refinement$original_residual <-
    reference$conditioning$eigensystem_refinement$original_residual

  observed[[1]]$conditioning <- observed[[2]]
  expect_identical(observed[[1]], reference, info = key)
}

test_that("clustered within fits match the 811ed24 reference across treatment and observation methods", {
  ref <- readRDS(test_path("fixtures", "r4_parity12_811ed24.rds"))
  m <- parity_panel()
  for (est in c("ipw", "reg", "dr")) for (adj in c("none", "ipw", "imputation", "dr")) {
    key <- paste(est, adj, sep = "_")
    args <- list(yname = "y", tname = "t", idname = "id", gname = "g", data = m, xformla = ~x, est_method = est,
      control_group = "nevertreated", gmm_weighting = "within", clustervars = "cl", bstrap = TRUE, cband = TRUE,
      biters = 199L)
    if (adj != "none") {
      args$selection <- "conditional_mtar"
      args$selection_est_method <- adj
      args$selection_xformla <- ~x
    }
    set.seed(101)
    rng_header <- .Random.seed[1L]
    f <- suppressWarnings(suppressMessages(do.call(att_gt, args)))
    r <- ref[[key]]
    IF <- as.matrix(f$inffunc)
    expect_equal(f$att, r$att, tolerance = 1e-10, info = key)
    expect_equal(f$se, r$se, tolerance = 1e-10, info = key)
    expect_equal(f$V_analytical, r$V, tolerance = 1e-10, info = key)
    expect_equal(f$c, r$c, tolerance = 1e-10, info = key)
    expect_equal(crossprod(IF), r$if_cross, tolerance = 1e-10, info = key)
    expect_equal(IF[1, ], r$if_row1, tolerance = 1e-10, info = key)
    expect_identical(f$chained$available, r$available, info = key)
    expect_parity_covariance_handling(f$chained$covariance_handling, r$covariance_handling, key)
    # R-devel includes the binomial RNG kind in the header. The live
    # header must stay unchanged and the reference state still matches exactly.
    expect_identical(.Random.seed[1L], rng_header, info = key)
    expect_identical(.Random.seed[-1L], r$rng[-1L], info = key)
  }
})

test_that("built-in cohort fits with omitted weighting use within", {
  m <- parity_panel()
  for (est in c("ipw", "reg", "dr")) {
    args <- list(yname = "y", tname = "t", idname = "id", gname = "g",
                 data = m, xformla = ~x, est_method = est,
                 bstrap = FALSE, cband = FALSE)
    omitted <- suppressMessages(do.call(att_gt, args))
    args$gmm_weighting <- "within"
    within <- suppressMessages(do.call(att_gt, args))
    expect_identical(omitted$att, within$att, info = est)
    expect_identical(omitted$inffunc, within$inffunc, info = est)
    expect_identical(omitted$V_analytical, within$V_analytical, info = est)
    expect_identical(omitted$chained$covariance_handling,
                     within$chained$covariance_handling, info = est)
    expect_identical(omitted$chained$covariance_handling$weighting_canonical,
                     "within", info = est)
  }
})
