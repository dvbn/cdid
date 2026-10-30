# Installed namespace integration of the accepted source prototype. No source
# production functions or globals are used; every binding is checked on lookup.
cp04_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}
cp04_reference <- new.env(parent = asNamespace("stats"))
invisible(capture.output(sys.source(testthat::test_path("reference", "cp04_finite_reference.R"), cp04_reference)))
cp04_control <- cp04_ns("ipw_owner_control")(max_chunk_rows = 23L,
  max_retained_entries = 1e7, max_matrix_entries = 2e7, max_work_ops = 1e9, max_output_entries = 1e7)
cp04_glm_control <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
cp04_evidence <- new.env(parent = emptyenv())
cp04_fixture <- function(order = seq_len(cp04_reference$N), Z = cp04_reference$Z,
                         V = cp04_reference$V, w = cp04_reference$w) {
  ref <- cp04_reference
  g <- c(3, 5, Inf)[ref$dat$k]
  long <- data.frame(id = rep(ref$ids, each = 4L), time = rep(1:4, ref$N),
    y = as.vector(t(ref$Y)), g = rep(g, each = 4L))
  # All but one entirely outcome-free identity are linked frame-only rows;
  # one explicit all-NA panel unit also exercises the owner's zero-observer map.
  zero <- which(ref$A == 0 & rowSums(ref$obs) == 0)
  long <- long[!long$id %in% zero[-1L], ]
  panel <- cp04_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", cp04_control)
  frame <- cp04_ns("construct_ipw_frame_owner")(panel, ref$ids[order], V[order, , drop = FALSE],
    weights = w[order], control = cp04_control, frame_g = g[order], allow_frame_only = TRUE)
  reader <- cp04_ns("construct_cp04_frame_reader")(frame, ref$ids[order], Z[order, , drop = FALSE],
    g[order], w[order], cp04_control)
  make <- function(s, t, K, id) cp04_ns("construct_ipw_link_view")(frame, s, t, g = 3, K = K,
    control = cp04_control, link_id = id,
    functional_key = paste0("cmpcat:reference-F-law:g1:K", if (length(K) == 1L) 2L else 1L))
  views <- list(make(1, 2, c(5, Inf), "NYT12"), make(2, 4, c(5, Inf), "NYT24"),
                make(1, 2, Inf, "NT12"), make(1, 2, c(Inf, 5), "horizon12"))
  list(panel = panel, frame = frame, reader = reader, views = views)
}

cp04_collect <- function(fit, chunk = 13L, batch = 4L, details = TRUE, scale = "frame") {
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  phi <- matrix(NA_real_, reader$state$N, length(fit$links))
  chunks <- list()
  info <- cp04_ns("replay_cp04_if")(fit, chunk, batch, function(b) {
    phi[b$rows, b$links] <<- b$phi
    if (details) chunks[[length(chunks) + 1L]] <<- b
  }, scale, details)
  list(phi = phi, chunks = chunks, info = info)
}

cp04_reference_stack <- function(fit, chunks) {
  rows <- list()
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  st <- reader$state
  for (b in chunks) {
    dd <- b$details; m <- dd$models[1:3]
    co <- function(name) do.call(cbind, lapply(m, function(x) x$coordinates[, name]))
    # pi/r are diagnostics only. Their rows are obtained independently from the
    # empirical marginal statistics, never appended to the production link IF.
    S <- cp04_reference$S[b$rows, 1]
    I <- cp04_reference$I[b$rows, , drop = FALSE]
    w <- st$w[b$rows]; A <- as.double(st$a_index[b$rows] > 0L)
    pi <- vapply(fit$selections[1:3], `[[`, numeric(1), "pi")
    r <- vapply(fit$selections[1:3], `[[`, numeric(1), "r")
    xpi <- sweep(w * I * sweep(matrix(S, length(S), 3), 2, pi, "-"), 2, st$mass_F, "/")
    xr <- sweep(w * I * (S - A * rep(r, each = length(S))), 2, st$mass_F * st$a, "/")
    rows[[length(rows) + 1L]] <- cbind(dd$a, do.call(cbind, dd$alpha),
      do.call(cbind, lapply(m, `[[`, "beta")), co("T"), co("M"), co("dp"), co("dm"),
      dd$eta, dd$tau, dd$lambda, xr, xpi,
      do.call(cbind, lapply(m, function(x) x$normalizers_A[, 1])),
      do.call(cbind, lapply(m, function(x) x$normalizers_A[, 2])))
  }
  do.call(rbind, rows)
}

testthat::test_that("actual full-F owners match independent point and complete case-refit rows", {
  fx <- cp04_fixture()
  testthat::expect_s3_class(fx$panel, "cdid_ipw_panel_owner")
  testthat::expect_true(cp04_ns("ipw_owner_seal_ok")(fx$panel, "panel_owner"))
  testthat::expect_gt(fx$frame$n_frame_only_units, 0)
  testthat::expect_gt(fx$panel$n_singleton_units, 0)
  testthat::expect_gt(fx$panel$n_unobserved_units, 0)
  fit <- cp04_ns("fit_cp04_corollary1")(fx$reader, fx$views, cp04_control, cp04_glm_control)
  testthat::expect_equal(fit$accounting$n_selection_fits, 6L)
  testthat::expect_equal(fit$accounting$n_treatment_fits, 3L)
  testthat::expect_identical(fit$links[[1]]$model, fit$links[[4]]$model)
  testthat::expect_identical(fit$links[[1]]$functional_key, fit$links[[4]]$functional_key)
  testthat::expect_true(startsWith(fit$links[[1]]$functional_key,
    "cmpcat:CP04-F-weighted-cohort-pair-logit-a/p-etaQG|"))
  testthat::expect_identical(fit$links[[1]]$source_functional_key, fx$views[[1]]$functional_key)
  testthat::expect_identical(fit$treatment_specs[[1]]$selection_ids[c(1, 3)],
                              fit$treatment_specs[[3]]$selection_ids)
  got <- cp04_collect(fit)
  stack <- cp04_reference_stack(fit, got$chunks)
  point <- c(fx$reader$state$a, unlist(lapply(fit$selections, `[[`, "coefficients")),
    unlist(lapply(fit$treatments, `[[`, "coefficients")),
    vapply(fit$treatments, function(x) x$means["T"], numeric(1)),
    vapply(fit$treatments, function(x) x$means["M"], numeric(1)),
    vapply(fit$treatments, function(x) x$means["dp"], numeric(1)),
    vapply(fit$treatments, function(x) x$means["dm"], numeric(1)))
  testthat::expect_lt(max(abs(point - cp04_reference$f$vector[1:33])), cp04_reference$tol)
  testthat::expect_equal(dim(stack), c(cp04_reference$N, 52L))
  testthat::expect_lt(max(abs(stack - cp04_reference$last)), cp04_reference$tol)
  fdphi <- cp04_reference$last[, 22:24] - cp04_reference$last[, 25:27]
  testthat::expect_lt(max(abs(got$phi[, 1:3] - fdphi)), cp04_reference$tol)
  testthat::expect_lt(max(abs(colMeans(stack))), cp04_reference$root_tol)
  testthat::expect_gt(diff(range(fx$reader$state$a)), .2)
  testthat::expect_gt(max(abs(got$phi[cp04_reference$A == 0, ])), .1)
  testthat::expect_gt(max(abs(got$phi[cp04_reference$S[, 1] == 0, 1])), .1)
  testthat::expect_identical(got$phi[, 1], got$phi[, 4])
  # Both count-normalized point summaries and factors remain on declared scales.
  lambda <- fx$reader$state$lambda
  core <- cp04_collect(fit, details = FALSE, scale = "core")
  testthat::expect_equal(core$phi, lambda * got$phi, tolerance = 1e-14)
  testthat::expect_equal(core$info$normalization_count, sum(cp04_reference$A))
  for (j in 1:3) {
    testthat::expect_equal(unname(fit$treatments[[j]]$means[c("dp", "dm")] / lambda),
      c(cp04_reference$f$dp[j], cp04_reference$f$dm[j]) / cp04_reference$f$lambda,
      tolerance = 1e-12)
  }
  cp04_evidence$fit <- fit; cp04_evidence$fixture <- fx
  cp04_evidence$stack <- stack; cp04_evidence$rows <- got
  cp04_evidence$errors <- c(point = max(abs(point - cp04_reference$f$vector[1:33])),
    all_52 = max(abs(stack - cp04_reference$last)), link = max(abs(got$phi[, 1:3] - fdphi)),
    normalizers = max(abs(stack[, c(28:33, 47:52)] - cp04_reference$last[, c(28:33, 47:52)])))
})

testthat::test_that("cross derivative omissions, wrong control sign and zero nonobserver rows fail", {
  fit <- cp04_evidence$fit; got <- cp04_evidence$rows
  gaps <- c(alpha_to_beta = 0, a_to_beta = 0, control_sign = 0)
  for (b in got$chunks) for (j in 1:3) {
    m <- b$details$models[[j]]; f <- fit$treatments[[j]]
    ia <- seq.int(1L, ncol(m$upstream), by = 3L)
    ial <- setdiff(seq_len(ncol(m$upstream)), ia)
    vb <- f$partials$beta[2, ]
    omit <- function(ix) {
      delta <- t(cp04_ns("cp04_solve")(f$H_chol, t(m$upstream[, ix, drop = FALSE] %*%
        t(f$cross[, ix, drop = FALSE]))))
      max(abs(drop(delta %*% vb)))
    }
    gaps["alpha_to_beta"] <- max(gaps["alpha_to_beta"], omit(ial))
    gaps["a_to_beta"] <- max(gaps["a_to_beta"], omit(ia))
    gaps["control_sign"] <- max(gaps["control_sign"], max(abs(2 * m$coordinates[, "M"])))
  }
  zeroed <- got$phi; zeroed[cp04_reference$A == 0, ] <- 0
  gaps <- c(gaps, zero_A0 = max(abs(zeroed - got$phi)),
    missing_N_over_nA = max(abs(got$phi - cp04_reference$f$lambda * got$phi)))
  testthat::expect_true(all(gaps > 100 * cp04_reference$tol))
  # Compare the deliberately wrong rows against actual FD, not just each other.
  testthat::expect_gt(max(abs(zeroed[, 1:3] -
    (cp04_reference$last[, 22:24] - cp04_reference$last[, 25:27]))), .1)
  cp04_evidence$negative_control_gaps <- gaps
})

testthat::test_that("bounded batches, exact reuse, and supplied identity order are preserved", {
  fit <- cp04_evidence$fit
  reference <- cp04_evidence$rows$phi
  for (chunk in c(1L, 7L, 23L)) for (batch in c(1L, 2L, 4L)) {
    got <- cp04_collect(fit, chunk, batch, FALSE)
    testthat::expect_equal(got$phi, reference, tolerance = 1e-13)
    testthat::expect_equal(got$info$batches,
      ceiling(cp04_reference$N / chunk) * ceiling(4 / batch))
  }
  b <- cp04_ns("evaluate_cp04_if_chunk")(fit, 1:13, c(1, 3, 4), details = TRUE)
  testthat::expect_equal(b$diagnostics$n_selection_score_evaluations, 3L)
  testthat::expect_equal(b$diagnostics$n_treatment_score_evaluations, 2L)
  ord <- c(seq.int(2L, cp04_reference$N, by = 2L), seq.int(1L, cp04_reference$N, by = 2L))
  fx <- cp04_fixture(ord)
  shuffled <- cp04_ns("fit_cp04_corollary1")(fx$reader, fx$views, cp04_control, cp04_glm_control)
  testthat::expect_equal(cp04_collect(shuffled, details = FALSE)$phi,
    reference[ord, ], tolerance = 1e-11)
  # A distinct reader is refitted, even with the same labels and dimensions.
  zz <- cp04_reference$Z; zz[, 2] <- zz[, 2] + cp04_reference$dat$h / 10
  fx2 <- cp04_fixture(Z = zz)
  distinct <- cp04_ns("fit_cp04_corollary1")(fx2$reader, fx2$views, cp04_control, cp04_glm_control)
  testthat::expect_equal(distinct$accounting$n_selection_fits, 6L)
  testthat::expect_gt(max(abs(cp04_collect(distinct, details = FALSE)$phi - reference)), 1e-3)
  fx3 <- cp04_fixture(V = matrix(1, cp04_reference$N, 1L))
  differing <- cp04_ns("fit_cp04_corollary1")(fx3$reader, fx3$views, cp04_control, cp04_glm_control)
  testthat::expect_length(differing$treatments[[1]]$coefficients, 1L)
  testthat::expect_length(differing$selections[[1]]$coefficients, 2L)
  for (s in fit$selections) {
    testthat::expect_equal(names(s), c("coefficients", "H_chol", "diagnostics", "pi", "r"))
    testthat::expect_equal(dim(s$H_chol), c(2L, 2L))
  }
  testthat::expect_false(fit$accounting$full_F_by_link_IF_retained)
  testthat::expect_false(fit$accounting$inverse_matrices_retained)
  testthat::expect_equal(fit$accounting$full_F_masks_per_model, 0L)
  testthat::expect_lt(fit$accounting$owned_entries + fit$accounting$borrowed_entries,
    fit$accounting$retained_bound)
  testthat::expect_lt(fit$accounting$work_ops_bound, 1e7)
  # Deterministic 2x permuted replication preserves point and row-mapped IFs.
  ref <- cp04_reference; N2 <- 2L * ref$N; orig2 <- rep(seq_len(ref$N), 2L)
  ids2 <- seq_len(N2); g2 <- c(3, 5, Inf)[ref$dat$k][orig2]; w2 <- ref$w[orig2]
  Z2 <- ref$Z[orig2, , drop = FALSE]; V2 <- ref$V[orig2, , drop = FALSE]; Y2 <- ref$Y[orig2, , drop = FALSE]
  zero0 <- which(ref$A == 0 & rowSums(ref$obs) == 0)
  fmask <- rep(FALSE, ref$N); fmask[zero0[-1L]] <- TRUE; p_units <- which(!fmask[orig2])
  long2 <- data.frame(id = rep(ids2[p_units], each = 4L), time = rep(1:4, length(p_units)),
    y = as.vector(t(Y2[p_units, , drop = FALSE])), g = rep(g2[p_units], each = 4L))
  ord2 <- c(seq.int(2L, N2, by = 2L), seq.int(1L, N2, by = 2L))
  panel_rep <- cp04_ns("construct_ipw_panel_owner")(long2, "id", "time", "y", "g", cp04_control)
  frame_rep <- cp04_ns("construct_ipw_frame_owner")(panel_rep, ids2[ord2], V2[ord2, , drop = FALSE],
    weights = w2[ord2], control = cp04_control, frame_g = g2[ord2], allow_frame_only = TRUE)
  reader_rep <- cp04_ns("construct_cp04_frame_reader")(frame_rep, ids2[ord2], Z2[ord2, , drop = FALSE],
    g2[ord2], w2[ord2], cp04_control)
  make_rep <- function(s, t, K, id) cp04_ns("construct_ipw_link_view")(frame_rep, s, t, g = 3, K = K,
    control = cp04_control, link_id = id,
    functional_key = paste0("cmpcat:reference-F-law:g1:K", if (length(K) == 1L) 2L else 1L))
  views_rep <- list(make_rep(1, 2, c(5, Inf), "NYT12"), make_rep(2, 4, c(5, Inf), "NYT24"),
    make_rep(1, 2, Inf, "NT12"), make_rep(1, 2, c(Inf, 5), "horizon12"))
  fit_rep <- cp04_ns("fit_cp04_corollary1")(reader_rep, views_rep, cp04_control, cp04_glm_control)
  got_rep <- cp04_collect(fit_rep, details = FALSE)
  testthat::expect_equal(got_rep$phi, reference[orig2[ord2], ], tolerance = 1e-11)
})

testthat::test_that("all-frame covariance and population share cross covariance have accepted scaling", {
  ref <- cp04_reference; fit <- cp04_evidence$fit
  phi <- cp04_evidence$rows$phi[, 1:3]; lambda <- ref$f$lambda; nA <- sum(ref$A)
  psi <- lambda * phi
  unit_v <- crossprod(scale(phi, center = TRUE, scale = FALSE)) / ref$N^2
  unit_core <- crossprod(scale(psi, center = TRUE, scale = FALSE)) / nA
  testthat::expect_equal(unit_core, nA * unit_v, tolerance = 1e-13)
  clusters <- 1L + ref$ids %% 9L
  grouped <- rowsum(psi, clusters)
  cluster_v <- crossprod(grouped) / nA^2
  testthat::expect_equal(cluster_v, crossprod(rowsum(phi, clusters)) / ref$N^2, tolerance = 1e-13)
  dropped <- crossprod(rowsum(psi[ref$A == 1, ], clusters[ref$A == 1])) / nA^2
  testthat::expect_gt(max(abs(dropped - cluster_v)), 1e-4)
  tau <- cp04_evidence$stack[, 37:39]
  testthat::expect_gt(max(abs(crossprod(phi, tau))), .1)
  testthat::expect_equal(cp04_collect(fit, details = FALSE)$phi[, 1:3], phi)
  # Distinct clusters containing exclusively A=0 rows also remain in sums.
  cl <- ifelse(ref$A == 1, 1L + ref$ids %% 7L, 20L + ref$ids %% 3L)
  testthat::expect_gt(sum(rowSums(abs(rowsum(psi, cl)))[as.numeric(rownames(rowsum(psi, cl))) >= 20] > 0), 0)
  cp04_evidence$covariance <- list(unit = unit_v, Omega_A = unit_core,
    cluster = cluster_v, dropped_A0_gap = max(abs(dropped - cluster_v)))
})

testthat::test_that("rank, separation, nonfinite, unsupported and resource inputs refuse", {
  fit <- cp04_evidence$fit; fx <- cp04_evidence$fixture
  X <- cbind(1, rep(c(-1, 1), each = 4)); y <- rep(c(0, 1, 0, 1), 2)
  testthat::expect_error(cp04_ns("cp04_logit")(cbind(X, X[, 2]), y, rep(1, 8), cp04_glm_control, 8, list()),
    class = "cdid_cp04_rank_error")
  testthat::expect_error(cp04_ns("cp04_logit")(X, rep(1, 8), rep(1, 8), cp04_glm_control, 8, list()),
    class = "cdid_cp04_separation_error")
  testthat::expect_error(cp04_ns("cp04_logit")(X, rep(0, 8), rep(1, 8), cp04_glm_control, 8, list()),
    class = "cdid_cp04_separation_error")
  testthat::expect_error(cp04_ns("cp04_logit")(X, y, c(Inf, rep(1, 7)), cp04_glm_control, 8, list()),
    class = "cdid_cp04_numeric_error")
  testthat::expect_error(cp04_ns("cp04_logit")(X, rep(.5, 8), rep(1, 8), cp04_glm_control, 8, list()),
    class = "cdid_cp04_unsupported_error")
  testthat::expect_error(cp04_ns("cp04_logit")(X, as.double(X[, 2] > 0), rep(1, 8), cp04_glm_control, 8, list()),
    class = "cdid_cp04_separation_error")
  badgc <- cp04_glm_control; badgc$maxit <- 1L
  testthat::expect_error(cp04_ns("fit_cp04_corollary1")(fx$reader, fx$views, cp04_control, badgc),
    class = "cdid_cp04_convergence_error")
  zz <- cp04_reference$Z; zz[cp04_reference$A == 0, 2] <- NA_real_
  testthat::expect_error(cp04_ns("construct_cp04_frame_reader")(fx$frame, cp04_reference$ids, zz,
    c(3, 5, Inf)[cp04_reference$dat$k], cp04_reference$w, cp04_control),
    class = "cdid_cp04_numeric_error")
  ww <- cp04_reference$w; ww[1] <- ww[1] + 1
  testthat::expect_error(cp04_ns("construct_cp04_frame_reader")(fx$frame, cp04_reference$ids, cp04_reference$Z,
    c(3, 5, Inf)[cp04_reference$dat$k], ww, cp04_control), class = "cdid_cp04_identity_error")
  testthat::expect_error(cp04_ns("fit_cp04_corollary1")(fx$reader, list(list(pair = c(1, 2))),
    cp04_control, cp04_glm_control), class = "cdid_ipw_owner_integrity_error")
  nv24 <- cp04_ns("construct_ipw_nuisance_view")(fx$frame, 2, 4, g = 3, K = c(5, Inf),
    control = cp04_control, nuisance_id = "nv24")
  lv_cross <- cp04_ns("construct_ipw_link_view")(fx$frame, 1, 2, nuisance_view = nv24,
    control = cp04_control, verify_nuisance_mask = FALSE, link_id = "cross",
    functional_key = "cmpcat:reference-F-law:g1:K1")
  testthat::expect_error(cp04_ns("fit_cp04_corollary1")(fx$reader, list(lv_cross), cp04_control, cp04_glm_control),
    class = "cdid_cp04_unsupported_error")
  testthat::expect_error(cp04_ns("evaluate_cp04_if_chunk")(fit, 1:24, 1), class = "cdid_cp04_input_error")
  testthat::expect_error(cp04_ns("evaluate_cp04_if_chunk")(fit, 1:3, 1, scale = "observer"), class = "cdid_cp04_input_error")
  testthat::expect_error(cp04_ns("fit_cp04_corollary1")(fx$reader, fx$views, cp04_control, cp04_glm_control,
    known_probabilities = 1), "unused argument")
  bad <- fit; bad$selections[[1]]$coefficients[1] <- 0
  testthat::expect_error(cp04_ns("evaluate_cp04_if_chunk")(bad, 1:3, 1), class = "cdid_cp04_integrity_error")
  bad_r <- fx$reader; bad_r$ref <- NULL; bad_r$panel_ref <- new.env(parent = emptyenv())
  bad_r <- cp04_ns("ipw_owner_seal")(bad_r, "cp04_reader")
  testthat::expect_error(cp04_ns("validate_cp04_reader")(bad_r), class = "cdid_cp04_integrity_error")
  low <- cp04_control; low$max_matrix_entries <- 100
  testthat::expect_error(cp04_ns("fit_cp04_corollary1")(fx$reader, fx$views, low, cp04_glm_control),
    class = "cdid_ipw_owner_bound_error")
  low <- cp04_control; low$max_work_ops <- 100
  testthat::expect_error(cp04_ns("fit_cp04_corollary1")(fx$reader, fx$views, low, cp04_glm_control),
    class = "cdid_ipw_owner_bound_error")
  low <- cp04_control; low$max_output_entries <- 10
  testthat::expect_error(cp04_ns("fit_cp04_corollary1")(fx$reader, fx$views, low, cp04_glm_control),
    class = "cdid_ipw_owner_bound_error")
  # Rebind ONLY explicit resource limits on real fitted models to test replay
  # preflight. No fake owner, numerical state, fit, or IF expectation is used.
  limited <- function(field) {
    x <- fit; x$ref <- NULL; x$control[[field]] <- 1
    cp04_ns("ipw_owner_seal")(x, "cp04_backend")
  }
  for (field in c("max_output_entries", "max_work_ops", "max_matrix_entries", "max_retained_entries")) {
    called <- FALSE
    testthat::expect_error(cp04_ns("replay_cp04_if")(limited(field), 7L, 2L, function(b) called <<- TRUE),
      class = "cdid_ipw_owner_bound_error")
    testthat::expect_false(called)
  }
  testthat::expect_equal(unname(cp04_ns("cp04_payload")(fit)["entries"]), fit$accounting$owned_entries)
  testthat::expect_equal(unname(cp04_ns("cp04_payload")(fx$reader)["entries"]), fx$reader$accounting$owned_entries)
})

testthat::test_that("actual all-success F logits refuse and only observed-risk outcomes are read", {
  # Actual native owners: cohort 3 has S12=1 on its entire declared F risk set.
  id <- 1:16; g <- rep(c(3, Inf), each = 8)
  long <- data.frame(id = rep(id, each = 3L), time = rep(1:3, 16),
    y = rep(c(1, 2, 3), 16), g = rep(g, each = 3L))
  long$y[long$id %in% 13:16 & long$time == 2] <- NA_real_
  panel <- cp04_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", cp04_control)
  X <- matrix(1, 16, 1L)
  frame <- cp04_ns("construct_ipw_frame_owner")(panel, id, X, rep(1, 16), cp04_control, frame_g = g)
  reader <- cp04_ns("construct_cp04_frame_reader")(frame, id, X, g, rep(1, 16), cp04_control)
  link <- cp04_ns("construct_ipw_link_view")(frame, 1, 2, g = 3, K = Inf, control = cp04_control,
    link_id = "all-success", functional_key = "cmpcat:reference-F-law:g1:K1")
  testthat::expect_error(cp04_ns("fit_cp04_corollary1")(reader, list(link), cp04_control, cp04_glm_control),
    class = "cdid_cp04_separation_error")
  # Overflowing differences in a cohort outside NT are never queried. In NYT
  # those same observed endpoints must refuse rather than silently omit rows.
  fx <- cp04_evidence$fixture
  p <- fx$panel; p$ref <- NULL
  target <- p$units$g[p$observations$unit_index] == 5
  p$observations$y[target & p$observations$time == 1] <- 1e308
  p$observations$y[target & p$observations$time == 2] <- -1e308
  # Reconstruct actual owners from the same long records, preserving missingness.
  oo <- p$observations
  data <- data.frame(id = p$units$id[oo$unit_index], time = oo$time, y = oo$y,
    g = p$units$g[oo$unit_index])
  panel2 <- cp04_ns("construct_ipw_panel_owner")(data, "id", "time", "y", "g", cp04_control)
  ref <- cp04_reference; gg <- c(3, 5, Inf)[ref$dat$k]
  frame2 <- cp04_ns("construct_ipw_frame_owner")(panel2, ref$ids, ref$V, ref$w, cp04_control,
    frame_g = gg, allow_frame_only = TRUE)
  reader2 <- cp04_ns("construct_cp04_frame_reader")(frame2, ref$ids, ref$Z, gg, ref$w, cp04_control)
  make <- function(K) cp04_ns("construct_ipw_link_view")(frame2, 1, 2, g = 3, K = K,
    control = cp04_control, link_id = "risk", functional_key = paste0("cmpcat:reference-F-law:g1:K",
      if (length(K) == 1L) 2L else 1L))
  nt <- cp04_ns("fit_cp04_corollary1")(reader2, list(make(Inf)), cp04_control, cp04_glm_control)
  testthat::expect_equal(cp04_collect(nt, batch = 1L, details = FALSE)$phi[, 1],
    cp04_evidence$rows$phi[, 3], tolerance = 1e-12)
  testthat::expect_error(cp04_ns("fit_cp04_corollary1")(reader2, list(make(c(5, Inf))),
    cp04_control, cp04_glm_control), class = "cdid_cp04_numeric_error")
})

cp04_u_reference <- new.env(parent = asNamespace("stats"))
sys.source(testthat::test_path("reference", "cp04_retained_u_reference.R"), cp04_u_reference)

cp04_known_rho <- function(ref, U) {
  # A fixed supplied normalized ratio; it is NOT reestimated in case refits.
  rho <- rep(NA_real_, ref$N)
  for (k in 1:3) for (x in -1:1) {
    kk <- ref$dat$k == k; xx <- ref$dat$x == x
    rho[U & kk & xx] <- (sum(ref$w[U & kk & xx])/sum(ref$w[U & kk])) /
      (sum(ref$w[kk & xx])/sum(ref$w[kk]))
  }
  rho
}

cp04_u_fit <- function(U, rho, universe, Z = cp04_reference$Z) {
  ref <- cp04_reference; fx <- cp04_evidence$fixture
  reader <- cp04_ns("construct_cp04_frame_reader")(fx$frame, ref$ids, Z,
    c(3,5,Inf)[ref$dat$k], ref$w, cp04_control, U, rho, universe,
    "Known normalized population-law transport supplied by design")
  cp04_ns("fit_cp04_corollary1")(reader, fx$views, cp04_control, cp04_glm_control)
}

cp04_first33 <- function(fit, chunks) {
  do.call(rbind, lapply(chunks, function(b) {
    d <- b$details; m <- d$models[1:3]
    co <- function(name) do.call(cbind, lapply(m, function(x) x$coordinates[,name]))
    cbind(d$a, do.call(cbind,d$alpha), do.call(cbind,lapply(m,`[[`,"beta")),
      co("T"),co("M"),co("dp"),co("dm"))
  }))
}

testthat::test_that("retained A/B use actual U scores, normalized known transport, and complete derivatives", {
  ref <- cp04_reference
  for (universe in c("B", "A")) {
    U <- if (universe == "A") ref$A == 1 else rowSums(ref$obs) > 0
    rho <- cp04_known_rho(ref,U)
    zz <- ref$Z; zz[!U,] <- NA_real_
    fit <- cp04_u_fit(U,rho,universe,zz)
    got <- cp04_collect(fit)
    actual <- cp04_first33(fit,got$chunks)
    fd <- matrix(0,ref$N,33L); h <- 5e-5
    for (i in seq_len(ref$N)) {
      plus <- minus <- rep(1,ref$N); plus[i] <- 1+h; minus[i] <- 1-h
      fd[i,] <- ref$N * (cp04_u_reference$cp04_finite_u(ref,U,rho,plus)$vector -
        cp04_u_reference$cp04_finite_u(ref,U,rho,minus)$vector)/(2*h)
    }
    baseline <- cp04_u_reference$cp04_finite_u(ref,U,rho)
    testthat::expect_lt(max(abs(actual-fd)),ref$tol)
    testthat::expect_equal(unname(vapply(fit$treatments,function(x) x$means["T"]-x$means["M"],0.0)),
      baseline$phi_point,tolerance=1e-11)
    testthat::expect_equal(got$phi[,1:3],actual[,22:24]-actual[,25:27],tolerance=1e-12)
    testthat::expect_equal(got$phi[!U,],matrix(0,sum(!U),4L),ignore_attr=TRUE,tolerance=1e-12)
    testthat::expect_lt(max(abs(colMeans(actual))),ref$root_tol)
    if (universe == "A") {
      testthat::expect_identical(unname(baseline$a),rep(1,3L))
      testthat::expect_equal(actual[,1:3],matrix(0,ref$N,3L),ignore_attr=TRUE)
    }
    # Cohort-dependent renormalization is not an innocuous transport change.
    badrho <- rho * c(1,1.5,.7)[ref$dat$k]
    bad <- cp04_u_reference$cp04_finite_u(ref,U,badrho)
    testthat::expect_gt(max(abs(bad$phi_point-baseline$phi_point)),1e-4)
    # F and U logits are different finite fits, even under rho=1 declarations.
    Uone <- cp04_u_fit(U,rep(1,ref$N),universe,zz)
    testthat::expect_gt(max(abs(cp04_collect(Uone,details=FALSE)$phi-cp04_evidence$rows$phi)),1e-3)
  }
})

testthat::test_that("saturated conditional-U logits transport correctly without assuming F logit equality", {
  ref <- cp04_reference; U <- ref$A == 1; rho <- cp04_known_rho(ref,U)
  Z <- cbind(1,ref$dat$x == 0,ref$dat$x == 1) * 1.0
  fitF <- cp04_u_fit(rep(TRUE,ref$N),rep(1,ref$N),"F",Z)
  fitU <- cp04_u_fit(U,rho,"A",Z)
  # Saturation gives exact fitted cell probabilities in their OWN universe.
  f <- cp04_u_reference$cp04_finite_u(ref,U,rho,Z=Z)
  testthat::expect_equal(vapply(fitU$treatments,function(x) unname(x$means["T"]-x$means["M"]),0.0),f$phi_point,tolerance=1e-11)
  testthat::expect_equal(vapply(fitU$treatments,function(x) unname(x$means["T"]-x$means["M"]),0.0),
    vapply(fitF$treatments,function(x) unname(x$means["T"]-x$means["M"]),0.0),tolerance=1e-10)
  testthat::expect_gt(max(abs(fitU$selections[[1]]$coefficients-fitF$selections[[1]]$coefficients)),.1)
})

cp04_public_fixture <- function() {
  ref <- cp04_reference
  frame <- data.frame(id=ref$ids,g=c(3,4,0)[ref$dat$k],x=ref$dat$x,w=ref$w,
    cluster=ifelse(ref$A==0,1000L+ref$ids,1L+ref$ids%%7L))
  data <- data.frame(id=rep(frame$id,each=4L),time=rep(1:4,ref$N),y=as.vector(t(ref$Y)),
    g=rep(frame$g,each=4L),x=rep(frame$x,each=4L),w=rep(frame$w,each=4L))
  data <- data[!data$id %in% frame$id[rowSums(ref$obs)==0],]
  list(data=data,frame=frame)
}

cp04_public_call <- function(fx=cp04_public_fixture(), ...) {
  args <- list(yname="y",tname="time",idname="id",gname="g",data=fx$data,frame=fx$frame,
    xformla=~x,selection="corollary1",est_method="ipw",weightsname="w",
    control_group="notyettreated",bstrap=FALSE,cband=FALSE,biters=19L,gmm_weighting="identity")
  dots <- list(...); args[names(dots)] <- dots
  do.call(cdid::att_gt,args)
}

testthat::test_that("public installed roster, cluster, aggregation and pretrend bindings use all F rows", {
  fx <- cp04_public_fixture(); mp <- cp04_public_call(fx); cl <- cp04_public_call(fx,clustervars="cluster")
  testthat::expect_s3_class(mp,"cdid_MP")
  testthat::expect_identical(mp$n,cp04_reference$N)
  testthat::expect_identical(rownames(mp$inffunc),as.character(fx$frame$id))
  testthat::expect_true(any(mp$chained$frame$frame_only))
  testthat::expect_true(any(mp$chained$frame$n_obs==1L))
  testthat::expect_gt(max(abs(mp$inffunc[cp04_reference$A==0,])),.1)
  testthat::expect_equal(mp$att,cl$att,tolerance=1e-12)
  testthat::expect_equal(mp$inffunc,cl$inffunc,tolerance=1e-12)
  grouped <- rowsum(cl$inffunc,fx$frame$cluster)
  centered <- sweep(grouped,2,colMeans(grouped),"-")
  testthat::expect_equal(unname(cl$V_analytical),unname(crossprod(centered)/mp$n),tolerance=1e-12)
  testthat::expect_equal(unname(cl$se),unname(sqrt(diag(cl$V_analytical)/mp$n)),tolerance=1e-12)
  omitted <- rowsum(cl$inffunc[cp04_reference$A==1,,drop=FALSE],fx$frame$cluster[cp04_reference$A==1])
  testthat::expect_gt(max(abs(crossprod(scale(omitted,scale=FALSE))/mp$n-cl$V_analytical)),1e-5)
  testthat::expect_true(cl$chained$pretest$status %in% c("computed","unavailable"))
  if (cl$chained$pretest$status=="unavailable") testthat::expect_true(nzchar(cl$chained$pretest$reason))
  ag <- cdid::aggte(cl,type="simple")
  testthat::expect_s3_class(ag,"cdid_AGGTEobj")
  testthat::expect_true(is.finite(ag$overall.att))
  # Installed public aggregation agrees with the accepted full-F share/IF path.
  out <- cp04_ns("aggregate_chained_ipw_targets")(theta_hat=cl$att,target_if_full=cl$inffunc,
    target_meta=cl$chained$target_meta,full_frame_g=cl$chained$frame$cohort,
    full_frame_w=cl$chained$frame$weight,type="simple",cluster_full=fx$frame$cluster,
    bstrap=FALSE,cband=FALSE,inference_accuracy=NULL,eligible_full=cp04_reference$A==1,available=cl$chained$available,
    multiplier_budget=cl$chained$aggregation_state$multiplier_budget,
    calendar=1:4,cohort_grid=c(3,4,Inf))
  testthat::expect_equal(ag$overall.att,out$overall$estimate,tolerance=1e-12)
  # Input permutation is an identity change only.
  fx$frame <- fx$frame[rev(seq_len(nrow(fx$frame))),]
  fx$data <- fx$data[rev(seq_len(nrow(fx$data))),]
  perm <- cp04_public_call(fx)
  testthat::expect_equal(perm$att,mp$att,tolerance=1e-10)
  testthat::expect_equal(unname(perm$inffunc[match(cp04_reference$ids,fx$frame$id),]),unname(mp$inffunc),tolerance=1e-9)
  testthat::expect_identical(mp$chained$selection$fitting_universe,"F")
  testthat::expect_identical(mp$chained$selection$link,"logit")
})

testthat::test_that("installed GMM core uses complete F replay with original core covariance scale", {
  fx <- cp04_public_fixture()
  prep <- cp04_ns("chained_request_prepare_inputs")(fx$data,fx$frame,"y","time","id","g",~x,
    weightsname="w",selection="corollary1")
  resolved <- cp04_ns("chained_public_resolve_control")(NULL,gmm_weighting="two_step")
  dims <- cp04_ns("chained_public_dimensions")(prep,biters=19L)
  ctrl <- cp04_ns("chained_public_pipeline_control")(dims,resolved,gmm_weighting="two_step")$control
  run <- function(weighting="two_step",pt="post",start=NULL,cg="notyettreated",bp="varying",ant=0L)
    cp04_ns("run_chained_ipw_pipeline")(data=prep$data_internal,idname="id",tname="time",yname="y",gname="g",
      control=ctrl,frame_id=prep$frame_ids,frame_Z=prep$frame_Z,frame_g=prep$frame_g,frame_weights=prep$frame_w,
      cp04=prep$cp04,allow_frame_only=TRUE,target_cohorts=prep$target_cohorts,run_numerical_transfer=FALSE,
      weighting=weighting,control_group=cg,base_period=bp,pt_assumption=pt,pt_start=start,anticipation=ant,
      aggregations=character(),bstrap=FALSE,cband=FALSE,biters=19L)
  res <- run(); b <- res$nuisance_registry$cp04
  testthat::expect_true(any(res$metadata$catalog$rows$original))
  testthat::expect_true(any(res$metadata$catalog$rows$n_horizons>0L))
  testthat::expect_identical(res$reduction_summary$n_links_used,res$reduction_summary$n_links_all)
  ant <- run(ant=1L)
  testthat::expect_equal(ant$controls$anticipation,1L)
  ids <- vapply(res$point_links_used,`[[`,integer(1),"cp04_index")
  phi <- matrix(0,cp04_reference$N,length(ids))
  cp04_ns("cp04_replay_links")(b,ids,function(x) phi[x$rows,] <<- x$phi,scale="frame")
  lambda <- sum(cp04_reference$A)/cp04_reference$N
  testthat::expect_equal(unname(res$covariance_stream$V),crossprod(scale(lambda*phi,scale=FALSE))/sum(cp04_reference$A),tolerance=1e-12)
  loading <- res$gmm_fit$internal_gauge$Lambda
  testthat::expect_equal(unname(res$full_frame_transport$Phi_full),unname(phi %*% loading),tolerance=1e-10)
  testthat::expect_true(res$full_frame_transport$complete_frame_rows)
  testthat::expect_equal(nrow(res$gmm_fit$if_projection$materialized_target_if),cp04_reference$N)
  # Public normalizers and cross terms against independent all-F link replay.
  testthat::expect_equal(res$gmm_fit$if_projection$materialized_target_if/lambda,
    res$full_frame_transport$Phi_full,tolerance=1e-12)
  nt <- run(cg="nevertreated"); testthat::expect_gt(nrow(nt$metadata$catalog$rows),0)
  all <- run(pt="all",start=1)
  vac <- run(pt="all",start=4)
  testthat::expect_equal(res$reporting$estimates,vac$reporting$estimates,tolerance=1e-10)
  testthat::expect_equal(res$reporting$if_full,vac$reporting$if_full,tolerance=1e-9)
  testthat::expect_gt(all$design_used$rank_diagnostics$n_restricted_pre_parameters,0)
  uni <- run(bp="universal")
  post <- res$reporting$target_meta$role=="post"; upost <- uni$reporting$target_meta$role=="post"
  testthat::expect_equal(res$reporting$estimates[post],uni$reporting$estimates[upost],tolerance=1e-12)
  testthat::expect_equal(res$reporting$if_full[,post],uni$reporting$if_full[,upost],tolerance=1e-12)
})

testthat::test_that("retained-U public missing-X and unsupported requests refuse truthfully", {
  fx <- cp04_public_fixture(); off <- cp04_reference$A==0
  fx$frame$x[off] <- NA_real_; fx$data$x[fx$data$id %in% fx$frame$id[off]] <- NA_real_
  mp <- cp04_public_call(fx,selection_universe="A",selection_justification="Repeat roster representative within cohort under the target law")
  testthat::expect_identical(mp$chained$selection$fitting_universe,"A")
  testthat::expect_equal(mp$inffunc[off,],matrix(0,sum(off),ncol(mp$inffunc)),ignore_attr=TRUE,tolerance=1e-12)
  testthat::expect_error(cp04_public_call(fx),class="cdid_chained_request_covariate_error")
  fx$frame$rho <- cp04_known_rho(cp04_reference,!off)
  transported <- cp04_public_call(fx,selection_universe="A",selection_rho="rho",
    selection_justification="Known normalized population-law ratio from the design")
  testthat::expect_true(all(is.finite(transported$att)))
  testthat::expect_error(cp04_public_call(selection_universe="A"),"selection_justification")
  # A valid weighted callback without D64's total weight gradient must refuse
  # after its output is validated, rather than being refused at the scope gate.
  callback_without_gradient <- function(y1, y0, D, covariates, i.weights, inffunc) {
    dy <- y1 - y0; w <- i.weights
    wt <- sum(w * D); wc <- sum(w * (1 - D))
    mt <- sum(w * D * dy) / wt; mc <- sum(w * (1 - D) * dy) / wc
    list(ATT = mt - mc, att.inf.func = length(D) * w *
      (D * (dy - mt) / wt - (1 - D) * (dy - mc) / wc))
  }
  e <- testthat::expect_error(cp04_public_call(est_method = callback_without_gradient),
    class = "cdid_custom_ipw_derivative_error")
  testthat::expect_match(conditionMessage(e), "log_weight_gradient", fixed = TRUE)
  testthat::expect_error(cp04_public_call(selection_link="probit"),"Unknown argument")
  testthat::expect_error(cp04_public_call(known_probabilities=1),"Unknown argument")
  testthat::expect_error(cp04_public_call(selection_universe="B",selection_rho="rho"),"selection_justification")
})

testthat::test_that("complete F cluster multipliers and aggregation share cross terms match independent rows", {
  fx <- cp04_public_fixture()
  prep <- cp04_ns("chained_request_prepare_inputs")(fx$data,fx$frame,"y","time","id","g",~x,
    weightsname="w",clustervars="cluster",selection="corollary1")
  dims <- cp04_ns("chained_public_dimensions")(prep,biters=19L)
  ctrl <- cp04_ns("chained_public_pipeline_control")(dims,
    cp04_ns("chained_public_resolve_control")(NULL,gmm_weighting="identity"),gmm_weighting="identity")$control
  res <- cp04_ns("run_chained_ipw_pipeline")(data=prep$data_internal,idname="id",tname="time",yname="y",gname="g",
    control=ctrl,frame_id=prep$frame_ids,frame_Z=prep$frame_Z,frame_g=prep$frame_g,frame_weights=prep$frame_w,
    cp04=prep$cp04,allow_frame_only=TRUE,target_cohorts=prep$target_cohorts,run_numerical_transfer=FALSE,
    weighting="identity",control_group="notyettreated",cluster_vec=prep$cluster_full,
    aggregations=character(),bstrap=FALSE,cband=FALSE,biters=19L)
  rm <- res$reporting; N <- prep$n_frame; groups <- unique(prep$cluster_full); G <- length(groups)
  sg <- ifelse(outer(1:19,seq_len(G),function(i,j) (i*j+i+floor(j/3))%%7)<3,-1,1)*1.0
  inference <- cp04_ns("compute_chained_ipw_inference")(rm$estimates,rm$if_full,
    cluster_full=prep$cluster_full,alp=.05,cband=TRUE,biters=19L,multiplier_budget=ctrl$multiplier_budget,
    inference_accuracy=NULL,fixed_signs=sg,eligible_full=cp04_reference$A==1,bstrap=TRUE)
  raw <- rowsum(rm$if_full,match(prep$cluster_full,groups),reorder=FALSE)
  center <- sweep(raw,2,colMeans(raw),"-")
  testthat::expect_equal(unname(inference$vcov),unname(crossprod(center)/N^2),tolerance=1e-12)
  testthat::expect_gt(inference$n_zero_eligible_clusters,0)
  draws <- sg %*% center / N
  scale_iqr <- apply(draws,2,function(z) (stats::quantile(z,.75,type=1)-stats::quantile(z,.25,type=1))/
    (stats::qnorm(.75)-stats::qnorm(.25)))
  testthat::expect_equal(unname(inference$se),unname(scale_iqr),tolerance=1e-12)
  critical <- as.numeric(stats::quantile(apply(abs(sweep(draws,2,scale_iqr,"/")),1,max),.95,type=1))
  testthat::expect_equal(inference$cband_crit_val,critical,tolerance=1e-12)
  agg <- cp04_ns("aggregate_chained_ipw_targets")(rm$estimates,rm$if_full,rm$target_meta,
    prep$frame_g,prep$frame_w,type="simple",cluster_full=prep$cluster_full,
    multiplier_budget=ctrl$multiplier_budget,inference_accuracy=NULL,bstrap=FALSE,cband=FALSE,
    eligible_full=cp04_reference$A==1,available=rm$available,calendar=1:4,cohort_grid=c(3,4,Inf))
  post <- which(rm$target_meta$time>=rm$target_meta$g)
  cellg <- rm$target_meta$g[post]; th <- rm$estimates[post]
  mass <- vapply(cellg,function(g) mean(prep$frame_w*(prep$frame_g==g)),0.0)
  denominator <- sum(mass); point <- sum(mass*th)/denominator
  fixed <- drop(rm$if_full[,post,drop=FALSE] %*% (mass/denominator))
  share <- numeric(N)
  for (g in unique(cellg)) share <- share + prep$frame_w*(prep$frame_g==g)*sum(th[cellg==g]-point)/denominator
  testthat::expect_equal(agg$overall$estimate,unname(point),tolerance=1e-12)
  testthat::expect_equal(unname(agg$overall$if_full),fixed+share,tolerance=1e-11)
  testthat::expect_gt(max(abs(share)),1e-5)
  testthat::expect_gt(abs(sum(fixed*share)),1e-5)
  testthat::expect_equal(agg$overall$se,sqrt(sum(scale(rowsum(fixed+share,prep$cluster_full),scale=FALSE)^2))/N,tolerance=1e-12)
  testthat::expect_false(isTRUE(all.equal(agg$overall$if_full,fixed,tolerance=1e-10)))
  # Mandatory resource refusal precedes every callback, with complete row bounds.
  small <- b <- res$nuisance_registry$cp04
  small$ref <- NULL; small$control$max_output_entries <- 1
  small <- cp04_ns("ipw_owner_seal")(small,"cp04_backend")
  called <- FALSE
  testthat::expect_error(cp04_ns("cp04_replay_links")(small,1L,function(x) called<<-TRUE),class="cdid_ipw_owner_bound_error")
  testthat::expect_false(called)
})

testthat::test_that("panel-declared F and retained-U risk validation preserve roster rows", {
  fx <- cp04_public_fixture()
  # Frame is optional: the supplied panel roster, including singletons, declares F.
  pr <- cp04_ns("chained_request_prepare_inputs")(fx$data,NULL,"y","time","id","g",~x,
    weightsname="w",selection="corollary1")
  testthat::expect_identical(pr$n_frame,length(unique(fx$data$id)))
  testthat::expect_true(any(pr$units$n_obs==1L))
  testthat::expect_true(all(pr$cp04$U))
  for (universe in c("A","B")) {
    U <- if (universe=="A") cp04_reference$A==1 else rowSums(cp04_reference$obs)>0
    ff <- fx; ff$frame$x[!U] <- NA_real_
    ff$data$x[ff$data$id %in% ff$frame$id[!U]] <- NA_real_
    ff$frame$rho <- cp04_known_rho(cp04_reference,U)
    prepared <- cp04_ns("chained_request_prepare_inputs")(ff$data,ff$frame,"y","time","id","g",~x,
      weightsname="w",selection="corollary1",selection_universe=universe,selection_rho="rho",
      selection_justification="Known normalized design ratio")
    testthat::expect_identical(prepared$cp04$U,U)
    testthat::expect_true(all(is.na(prepared$cp04$Z[!U,])))
    testthat::expect_true(all(is.na(prepared$cp04$rho[!U])))
    testthat::expect_identical(prepared$n_frame,cp04_reference$N)
    ff$frame$x[which(U)[1L]] <- NA_real_
    testthat::expect_error(cp04_ns("chained_request_prepare_inputs")(ff$data,ff$frame,"y","time","id","g",~x,
      selection="corollary1",selection_universe=universe,
      selection_justification="Known normalized design ratio"),class="cdid_chained_request_covariate_error")
  }
  sep <- fx; sep$frame$x_trt <- sep$frame$x
  sep$frame$x_trt[cp04_reference$A==0] <- NA_real_
  sep$data$x_trt <- sep$frame$x_trt[match(sep$data$id,sep$frame$id)]
  sep_prep <- cp04_ns("chained_request_prepare_inputs")(sep$data,sep$frame,"y","time","id","g",~x,
    weightsname="w",selection="corollary1",selection_xformla=~1)
  sep_tr <- cp04_ns("chained_request_prepare_inputs")(sep$data,sep$frame,"y","time","id","g",~x_trt,
    weightsname="w",selection="corollary1",selection_xformla=~x)
  testthat::expect_true(all(is.na(sep_tr$frame_Z[cp04_reference$A==0,])))
  testthat::expect_true(all(is.finite(sep_tr$frame_Z[cp04_reference$A==1,])))
  testthat::expect_true(all(is.finite(sep_tr$cp04$Z)))
  testthat::expect_equal(unname(sep_tr$frame_Z[cp04_reference$A==1,]),unname(sep_prep$frame_Z[cp04_reference$A==1,]))
  sep$frame$x_trt[which(cp04_reference$A==1)[1L]] <- NA_real_
  testthat::expect_error(cp04_ns("chained_request_prepare_inputs")(sep$data,sep$frame,"y","time","id","g",~x_trt,
    weightsname="w",selection="corollary1",selection_xformla=~x),class="cdid_chained_request_covariate_error")
  bad <- fx; bad$frame$g <- NULL
  testthat::expect_error(cp04_ns("chained_request_prepare_inputs")(bad$data,bad$frame,"y","time","id","g",~x,
    selection="corollary1"),"known cohorts")
})


testthat::test_that("retained-U containment and known-boundary events refuse without invention", {
  ref <- cp04_reference; fx <- cp04_evidence$fixture
  U <- rep(TRUE,ref$N); U[which(ref$A==1)[1L]] <- FALSE
  testthat::expect_error(cp04_ns("construct_cp04_frame_reader")(fx$frame,ref$ids,ref$Z,
    c(3,5,Inf)[ref$dat$k],ref$w,cp04_control,U,rep(1,ref$N),"B"),class="cdid_cp04_support_error")
  rho <- rep(1,ref$N); rho[1] <- 0
  testthat::expect_error(cp04_ns("construct_cp04_frame_reader")(fx$frame,ref$ids,ref$Z,
    c(3,5,Inf)[ref$dat$k],ref$w,cp04_control,rep(TRUE,ref$N),rho,"F"),class="cdid_cp04_numeric_error")
  helpers <- new.env(parent=baseenv())
  sys.source(system.file("development","cp04_demo_fixture.R",package="cdid",mustWork=TRUE),helpers)
  d <- helpers$cp04_demo_fixture()
  testthat::expect_error(cdid::att_gt(yname="y",tname="period",idname="id",gname="G",
    data=d$panel,frame=d$frame,xformla=~x,selection="corollary1",selection_universe="A",
    selection_justification="Example declaration only; fitted all-success event must refuse",est_method="ipw",
    bstrap=FALSE,cband=FALSE),class="cdid_cp04_separation_error")
})


testthat::test_that("unavailable public pretrend retains the full requested hypothesis", {
  fx <- cp04_public_fixture()
  fx$data$y[fx$data$g==3 & fx$data$time==1] <- NA_real_
  mp <- cp04_public_call(fx)
  testthat::expect_identical(mp$chained$pretest$status,"unavailable")
  testthat::expect_match(mp$chained$pretest$reason,"full declared null")
  testthat::expect_true(any(!mp$chained$available[mp$chained$role=="pre"]))
})
