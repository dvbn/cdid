# Separate conditional-process composition. Existing saved process fixtures
# explicitly retain pooled outcomes; these checks reuse their raw finite data.
opcpt_ns <- function(name) getFromNamespace(name, "cdid")
opcpt_data <- function(name) utils::read.csv(testthat::test_path("fixtures", "cpt", name))
opcpt_args <- function(data, overrides = list()) {
  a <- list(yname = "y", tname = "time", idname = "id", gname = "G", xformla = ~x + z,
    data = data, panel = TRUE, allow_unbalanced_panel = NULL, control_group = "notyettreated", weightsname = "w",
    alp = 0.05, bstrap = TRUE, cband = TRUE, biters = 9L, clustervars = NULL, est_method = "reg",
    print_details = FALSE, pl = FALSE, cores = 1, anticipation = 0, base_period = "varying", fix_weights = NULL,
    mode = NULL, frame = NULL, selection = "none", selection_model = "pair_logit", selection_pooling = "separate",
    selection_xformla = NULL, selection_universe = "F", selection_rho = NULL, selection_justification = NULL,
    selection_est_method = "ipw", imputation_pooling = "separate", control = NULL, keep_rows = TRUE,
    fixed_signs = matrix(rep(c(1, -1), length.out = 9L * length(unique(data$id))), 9L))
  for (nm in names(overrides)) a[nm] <- list(overrides[[nm]])
  a
}
opcpt_run <- function(data, overrides = list()) do.call(opcpt_ns("chained_conditional_did_pretest"), opcpt_args(data, overrides))
opcpt_index <- function(state, id) {
  for (d in state$aux_plan$selections[[state$base_period]]) {
    key <- paste0("pre(g=", opcpt_ns("fmt_design_num")(d$g), ",t=", opcpt_ns("fmt_design_num")(d$t), ")")
    if (identical(key, id)) return(d$index)
  }
  stop("Missing auxiliary descriptor")
}
opcpt_backend_rows <- function(state, q) {
  b <- state$registry$cp04; j <- state$registry$cp04_aux_index[q]
  out <- numeric(state$n_total)
  for (s in seq.int(1L, state$n_total, by = 7L)) {
    rr <- seq.int(s, min(state$n_total, s + 6L))
    out[rr] <- opcpt_ns("chained_full_if_chunk")(b, rr, j, "frame")$phi[, 1L]
  }
  out
}

testthat::test_that("separate conditional extractors recover complete link IFs on fixed union fits", {
  A <- opcpt_data("fA.csv"); B <- opcpt_data("fB.csv"); Bd <- opcpt_data("fBd.csv")
  pair <- list(selection = "conditional_mtar", selection_xformla = ~sx + sz)
  runs <- list(
    baseline_or = list(A, list()), baseline_dr = list(A, list(est_method = "dr")),
    legacy_or = list(B, pair),
    general_dr = list(B, c(pair, list(est_method = "dr", selection_model = "independent_period_logit", selection_pooling = "within_cohort"))),
    imputation_reg = list(B, list(selection = "conditional_mtar", selection_est_method = "imputation", imputation_pooling = "across_cohorts")),
    dated_or = list(Bd, list(fix_weights = "base_period", base_period = "universal")),
    dated_imputation_dr = list(Bd, list(est_method = "dr", fix_weights = "base_period", selection = "conditional_mtar",
      selection_est_method = "imputation", imputation_pooling = "within_cohort")),
    dated_augmented_dr = list(Bd, c(pair, list(est_method = "dr", fix_weights = "base_period", selection_est_method = "dr",
      selection_model = "independent_period_logit", selection_pooling = "within_cohort", imputation_pooling = "within_cohort"))))
  for (nm in names(runs)) {
    rr <- runs[[nm]]
    result <- opcpt_run(rr[[1L]], c(rr[[2L]], list(outcome_pooling = "separate")))
    testthat::expect_identical(result$chained$status, "available", info = nm)
    testthat::expect_identical(result$chained$request$outcome_pooling, "separate", info = nm)
    state <- attr(result, "state"); backend <- state$registry$cp04
    testthat::expect_identical(backend$outcome_pooling, "separate", info = nm)
    frozen <- serialize(backend, NULL)
    for (cell in attr(result, "cells")) {
      q <- opcpt_index(state, cell$id)
      testthat::expect_equal(opcpt_ns("cpt_cell_scalar_rows")(cell), cell$sigma * opcpt_backend_rows(state, q),
        tolerance = 1e-9, info = paste(nm, cell$id))
      testthat::expect_equal(cell$sigma * cell$mu1, cell$estimate_scalar, tolerance = 1e-10)
      testthat::expect_equal(cell$denominator_check$supplied, cell$denominator_check$accumulated, tolerance = 1e-12)
      slopes <- grep("^(or_beta|ob):", cell$coords, value = TRUE)
      testthat::expect_true(all(grepl("^(or_beta|ob):[0-9]+:[0-9]+$", slopes)), info = nm)
      # Replaying distinct row chunks reads the same sealed union nuisance fit.
      sel <- list(id = cell$id, sign = cell$sigma)
      extractor <- opcpt_ns("cpt_cell_extractor")(state, q, sel, 5L)
      testthat::expect_identical(extractor$spec$coords, cell$coords)
      for (k in seq_len(extractor$n_chunks)) extractor$pieces(k)
    }
    testthat::expect_identical(serialize(backend, NULL), frozen, info = nm)
  }
})

testthat::test_that("conditional omission equals explicit separate with fixed signs and bounded chunks", {
  d <- opcpt_data("fA.csv")
  omitted <- opcpt_run(d)
  explicit <- opcpt_run(d, list(outcome_pooling = "separate", control = list(chunk_rows = 7L)))
  testthat::expect_equal(omitted$CvM, explicit$CvM, tolerance = 1e-10)
  testthat::expect_equal(omitted$CvMb, explicit$CvMb, tolerance = 1e-10)
  for (j in seq_along(attr(omitted, "cells"))) {
    left <- attr(omitted, "cells")[[j]]; right <- attr(explicit, "cells")[[j]]
    testthat::expect_equal(opcpt_ns("cpt_cell_process_rows")(left), opcpt_ns("cpt_cell_process_rows")(right), tolerance = 1e-9)
  }
})

# Independent raw-data frequency oracle for a frozen auxiliary comparison and
# frozen grid. Every evaluation refits cohort observation rates and control WLS
# from raw outcomes, regressors and weights. It reads fit identities only, never
# fitted coefficients, nuisance IFs, loading operators or production predictors.
opcpt_raw_or <- function(data, ids, cohorts, pair, controls, target, threshold, freq) {
  u <- data[match(ids, data$id), , drop = FALSE]
  G <- u$G; G[G == 0] <- Inf
  ci <- match(G, cohorts)
  y_at <- function(t) {
    r <- data[data$time == t, , drop = FALSE]
    r$y[match(ids, r$id)]
  }
  yl <- y_at(pair$s); yr <- y_at(pair$t)
  S <- is.finite(yl) & is.finite(yr); dy <- yr - yl
  nobs <- vapply(ids, function(id) sum(is.finite(data$y[data$id == id])), integer(1L))
  A <- nobs >= 2L
  rate <- vapply(seq_along(cohorts), function(k) {
    ii <- which(ci == k)
    den <- sum(freq[ii] * u$w[ii] * A[ii])
    if (den == 0) return(0)
    sum(freq[ii] * u$w[ii] * S[ii]) / den
  }, numeric(1L))
  q <- freq * u$w / pmax(rate[ci], 1e-6)
  dates <- sort(unique(data$time)); L <- match(pair$t, dates) - match(pair$s, dates)
  X <- cbind(1, L * as.matrix(u[c("x", "z")]))
  control_rows <- which(S & ci %in% controls)
  beta <- stats::lm.wfit(X[control_rows, , drop = FALSE], dy[control_rows], q[control_rows])$coefficients
  stopifnot(all(is.finite(beta)))
  tt <- which(S & ci == target)
  h <- as.double(u$x[tt] <= threshold[2L] & u$z[tt] <= threshold[3L])
  # The conditional process restricts the numerator, keeping its original
  # treated denominator. All N original-unit frequencies remain in each fit.
  sum(q[tt] * h * (dy[tt] - as.vector(X[tt, , drop = FALSE] %*% beta))) / sum(q[tt])
}

testthat::test_that("nonconstant separate conditional OR loadings match full raw-frequency refits", {
  d <- opcpt_data("fA.csv"); result <- opcpt_run(d, list(outcome_pooling = "separate"))
  state <- attr(result, "state"); b <- state$registry$cp04
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE); st <- reader$state
  cell <- attr(result, "cells")[[1L]]; q <- opcpt_index(state, cell$id)
  link <- b$links[[state$registry$cp04_aux_index[q]]]
  pool <- b$regression_specs[[link$regression]]
  support <- b$support_specs[[link$support_id]]
  testthat::expect_identical(opcpt_ns("or_outcome_block_fit")(b$outcome_groups, link$regression)$p, 2L)
  raw <- function(freq, threshold) cell$sigma * opcpt_raw_or(d, st$id, st$cohorts, support$pair,
    pool$cohorts, support$cohort, threshold, freq)
  # Select two nonconstant orthants of the frozen original-unit grid.
  grid <- cell$zd; all_x <- d$x[match(st$id, d$id)]; all_z <- d$z[match(st$id, d$id)]
  mass <- vapply(seq_len(nrow(grid)), function(j) sum(all_x <= grid[j, 2L] & all_z <= grid[j, 3L]), integer(1L))
  js <- head(which(mass > 0L & mass < length(st$id)), 2L)
  testthat::expect_length(js, 2L)
  expected_rows <- opcpt_ns("cpt_cell_process_rows")(cell)
  N <- length(st$id); eps <- 1e-6
  probes <- unique(c(1L, which(st$cohort_index == support$cohort)[c(2L, 5L)],
    which(st$cohort_index %in% pool$cohorts)[c(1L, 5L)], N))
  for (j in js) {
    testthat::expect_equal(raw(rep(1, N), grid[j, ]), cell$J[j], tolerance = 1e-10)
    for (i in probes) {
      direction <- N * as.double(seq_len(N) == i) - 1
      fd <- (raw(1 + eps * direction, grid[j, ]) - raw(1 - eps * direction, grid[j, ])) / (2 * eps)
      testthat::expect_equal(fd, expected_rows[i, j], tolerance = 1e-5, info = paste("row", i, "grid", j))
    }
  }
  H <- outer(all_x, grid[js, 2L], "<=") & outer(all_z, grid[js, 3L], "<=")
  shortcut <- H * opcpt_ns("cpt_cell_scalar_rows")(cell)
  testthat::expect_gt(max(abs(expected_rows[, js, drop = FALSE] - shortcut)), 1e-5)
})
