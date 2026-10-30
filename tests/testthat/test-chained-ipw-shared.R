# Substantive tests for the private CP02 shared-nuisance IPW layer and its
# structural class certificate (gmm_computation.tex Prop. 5 / Cor. 6).
# All glm.control values, comparison tolerances and eigenvalue ratios below are
# test-only fixture settings; they define no production threshold or policy.

if (!("cdid" %in% loadedNamespaces())) {
  stop("Load the candidate cdid namespace before running this test file.", call. = FALSE)
}

shared_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)
shared_bounds <- list(max_links = 50, max_classes = 10, max_class_dates = 10,
                      max_E_entries = 200, max_mask_entries = 1e6)

# Risk units (cohorts 4, 7, Inf) observe all of dates 1..6 or at most one of
# them (pair non-observers, kept A-eligible by dates 7..8). Cohort 9 is outside
# every declared risk set and observes arbitrary dates.
shared_fixture <- function(seed = 20261014, n = 180L) {
  set.seed(seed)
  G <- sample(c(4, 7, Inf, 9), n, TRUE, prob = c(0.3, 0.25, 0.3, 0.15))
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  Z <- cbind(intercept = 1, x1 = x1, x2 = x2)
  w <- runif(n, 0.5, 2)
  obs <- matrix(FALSE, n, 8L)
  for (i in seq_len(n)) {
    if (G[i] == 9) {
      obs[i, sample(1:8, sample(2:6, 1L))] <- TRUE
    } else if (runif(1L) < 0.75) {
      obs[i, 1:6] <- TRUE
      obs[i, 7:8] <- runif(2L) < 0.5
    } else {
      k <- sample(0:6, 1L)
      if (k > 0L) obs[i, k] <- TRUE
      obs[i, 7:8] <- TRUE
    }
  }
  Y <- rnorm(n) + matrix(rnorm(n * 8L), n, 8L) + outer(x1, (1:8) / 4) +
    outer(G == 4, (1:8 >= 4) * 1.5) + outer(x2 * (G == 7), rep(0.3, 8L))
  panel <- data.frame(id = rep(seq_len(n), each = 8L), time = rep(as.double(1:8), n),
                      y = as.vector(t(ifelse(obs, Y, NA_real_))), g = rep(G, each = 8L))
  pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z, weights = w)
  list(n = n, G = G, Z = Z, w = w, obs = obs, Y = Y, panel = panel, pidx = pidx, ctx = ctx)
}

pairs_A <- list(c(1, 2), c(2, 3), c(1, 3), c(3, 4), c(2, 4), c(1, 4))
pairs_B <- list(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6), c(4, 6))
pair_name <- function(p) paste0(p[1L], p[2L])

declare_all <- function(ses, nu, g, K, pairs, fk) {
  lapply(pairs, function(p) cdid:::declare_ipw_shared_link(
    ses, nu, g, K, p[1L], p[2L], shared_ctrl, link_id = paste0(fk, "_", pair_name(p)), functional_key = fk))
}

design_for <- function(pairs, fks, g = 4, ref = 3) {
  cdid:::build_chained_design(data.frame(
    g = g, s = vapply(pairs, `[`, numeric(1L), 1L), t = vapply(pairs, `[`, numeric(1L), 2L),
    ref = ref, functional_key = fks, stringsAsFactors = FALSE))
}

load_original_kernel <- function() {
  path <- testthat::test_path("reference", "chained_ipw_link_original_6550dc.R")
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required to bind the frozen kernel.")
  expect_identical(digest::digest(file = path, algo = "sha256"),
                   "6550dc592451cb7e9b7f45573d083fa14cacc1f3a09421a3330d31997183b0b5")
  env <- new.env(parent = baseenv())
  env$cdid_abort <- cdid:::cdid_abort
  ipw_reference_source(path, envir = env)
  env
}

load_potential_reference <- function() {
  env <- new.env(parent = baseenv())
  sys.source(testthat::test_path("reference", "cp02_shared_ipw_potential_reference.R"), envir = env)
  env
}

strip_shared <- function(f) {
  f <- unclass(f)
  f[setdiff(names(f), c("context_ref", "shared", "nuisance_ref", "declaration_ref", "link_ref"))]
}

# Scoped spy: private copies of the namespace functions whose enclosure is a
# child environment holding counting wrappers. The namespace is never edited.
make_spy <- function() {
  ns <- asNamespace("cdid")
  spy <- new.env(parent = ns)
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get(nm, envir = ns, inherits = FALSE)
    if (is.function(f) && identical(environment(f), ns)) {
      environment(f) <- spy
      assign(nm, f, envir = spy)
    }
  }
  counts <- new.env()
  counts$glm <- 0L
  counts$outcomes <- 0L
  counts$outcome_pairs <- list()
  real_glm <- get("ipw_glm_fit_internal", envir = ns)
  real_out <- get("ipw_pair_outcomes", envir = ns)
  spy$ipw_glm_fit_internal <- function(...) {
    counts$glm <- counts$glm + 1L
    real_glm(...)
  }
  spy$ipw_pair_outcomes <- function(context, hit) {
    counts$outcomes <- counts$outcomes + 1L
    if (length(hit$obs_s)) {
      tm <- context$observations$time
      counts$outcome_pairs[[length(counts$outcome_pairs) + 1L]] <- c(tm[hit$obs_s[1L]], tm[hit$obs_t[1L]])
    }
    real_out(context, hit)
  }
  list(fn = spy, counts = counts)
}

test_that("shared fits and complete IFs equal the immutable 6550dc kernel plus ONE matched-start adaptation bitwise and the independent potential reference", {
  fx <- shared_fixture()
  orig <- load_original_kernel()
  ref <- load_potential_reference()
  ses <- cdid:::construct_ipw_shared_session(fx$ctx, max_nuisance = 2, max_links = 20, max_unit_entries = 1e5)
  nu <- cdid:::construct_ipw_nuisance(ses, 4, c(7, Inf), 2, 5, shared_ctrl)
  decl <- declare_all(ses, nu, 4, c(7, Inf), pairs_B, "opA")
  fits <- lapply(decl, function(d) cdid:::fit_ipw_link_shared(ses, d))
  names(fits) <- vapply(pairs_B, pair_name, character(1L))

  ctx_o <- orig$construct_ipw_context(fx$pidx, id = seq_len(fx$n), Z = fx$Z, weights = fx$w)
  # The frozen kernel accepts three fields; native-only fields must be defaults.
  legacy3 <- shared_ctrl[c("epsilon", "maxit", "trace")]
  expect_identical(do.call(stats::glm.control, legacy3), shared_ctrl)
  ofits <- lapply(pairs_B, function(p) orig$fit_ipw_link(
    ctx_o, 4, c(7, Inf), p[1L], p[2L], legacy3, link_id = paste0("opA_", pair_name(p)), functional_key = "opA"))
  names(ofits) <- names(fits)
  for (k in seq_along(fits)) {
    sf <- strip_shared(fits[[k]])
    so <- strip_shared(ofits[[k]])
    expect_identical(names(sf$diagnostics), c(names(so$diagnostics), "clipping"))
    expect_identical(sf$diagnostics[names(so$diagnostics)], so$diagnostics)
    expect_true(is.list(sf$diagnostics$clipping) && identical(sf$diagnostics$clipping$eps, 1e-6))
    expect_identical(as.integer(sf$diagnostics$clipping$p_n), fits[[k]]$counts$n_observed_risk_units)
    expect_identical(as.integer(sf$diagnostics$clipping$r_n), fits[[k]]$cohorts$n_cohorts)
    sf$diagnostics <- sf$diagnostics[names(so$diagnostics)]
    expect_identical(sf, so)
  }

  for (cs in c(1L, 17L, fx$n)) {
    m <- cdid:::materialize_ipw_if(fx$ctx, fits, max_entries = 6 * fx$n, chunk_size = cs, include_components = TRUE)
    mo <- orig$materialize_ipw_if(ctx_o, ofits, max_entries = 6 * fx$n, chunk_size = cs, include_components = TRUE)
    expect_identical(m$phi_matrix, mo$phi_matrix)
    expect_identical(m$components, mo$components)
    expect_identical(m$id, mo$id)
  }

  # Independent reference: own Newton logit, accepted weighted-pair IF, potentials.
  risk <- fx$G %in% c(4, 7, Inf)
  S <- as.numeric(risk & fx$obs[, 2L] & fx$obs[, 5L])
  rn <- ref$ref_shared_nuisance(S, fx$G, fx$Z, fx$w, 4, c(7, Inf))
  expect_lt(max(abs(rn$beta - nu$stage$beta)), 1e-9)
  expect_lt(max(abs(rn$r - nu$stage$rates)), 1e-14)
  psi <- lapply(1:6, function(u) ref$ref_date_potential(rn, fx$Y[, u]))
  for (k in seq_along(pairs_B)) {
    p <- pairs_B[[k]]
    delta <- ifelse(S == 1, fx$Y[, p[2L]] - fx$Y[, p[1L]], 0)
    rl <- ref$ref_link_if(rn, delta)
    expect_lt(abs(fits[[k]]$estimate - rl$estimate), 1e-10)
    expect_lt(max(abs(m$phi_matrix[, k] - rl$phi)), 1e-8)
    # Prop. 5: mu_l = lambda_t - lambda_s and phi_l = psi_t - psi_s for every unit.
    expect_lt(abs(fits[[k]]$estimate - (psi[[p[2L]]]$lambda - psi[[p[1L]]]$lambda)), 1e-10)
    expect_lt(max(abs(m$phi_matrix[, k] - (psi[[p[2L]]]$psi - psi[[p[1L]]]$psi))), 1e-8)
  }
  # Outside-risk units: exact zeros; risk-set pair non-observers: nonzero nuisance terms.
  expect_gt(sum(!risk), 10L)
  expect_identical(max(abs(m$phi_matrix[!risk, ])), 0)
  nonobs <- risk & S == 0
  expect_identical(sum(nonobs), nu$mask$n_pair_nonobservers)
  expect_gt(sum(nonobs), 20L)
  expect_gt(min(apply(abs(m$phi_matrix[nonobs, ]), 2L, max)), 1e-3)
  expect_identical(max(abs(m$components$treated[nonobs, ])), 0)
  expect_identical(max(abs(m$components$control[nonobs, ])), 0)
  # Treated/control denominators and all cohort rate/Jacobian columns are retained.
  expect_identical(fits[[1L]]$means$sum_a, nu$stage$sum_a)
  expect_identical(dim(fits[[1L]]$partials$A_beta_r), c(3L, 3L))
  expect_true(all(abs(fits[[1L]]$partials$A_beta_r[, "g"]) > 0))
  # Outside-risk units still change the link's own pair count, never the mask.
  expect_gt(length(unique(vapply(fits, function(f) f$counts$n_observed_pair_units, integer(1L)))), 1L)
})

test_that("balanced multi-date classes: covariance rank equals class rank, canonical forests span disconnected components, E self-test", {
  fx <- shared_fixture()
  ses <- cdid:::construct_ipw_shared_session(fx$ctx, max_nuisance = 4, max_links = 40, max_unit_entries = 1e5)
  nuA <- cdid:::construct_ipw_nuisance(ses, 4, c(7, Inf), 1, 2, shared_ctrl)
  nuB <- cdid:::construct_ipw_nuisance(ses, 4, Inf, 1, 2, shared_ctrl)
  decl <- c(declare_all(ses, nuA, 4, c(7, Inf), pairs_A, "opA"),
            declare_all(ses, nuB, 4, Inf, pairs_B, "opB"))
  des <- design_for(c(pairs_A, pairs_B), rep(c("opA", "opB"), each = 6L))

  # Production order: certify outcome-free declarations, then fit retained rows only.
  sp <- make_spy()
  cert <- sp$fn$certify_ipw_structural_classes(ses, decl, des, 1:12, "optimal", shared_bounds)
  expect_identical(sp$counts$outcomes, 0L)
  expect_identical(sp$counts$glm, 0L)
  expect_identical(cert$link_kind, rep("shared_declaration", 12L))
  expect_identical(length(cert$classes), 2L)
  cA <- cert$classes[[1L]]
  cB <- cert$classes[[2L]]
  expect_identical(cA$dates, c(1, 2, 3, 4))
  expect_identical(cA$n_components, 1L)
  expect_identical(cA$cycle_rank, 3L)
  expect_identical(cA$forest_link_index, c(1L, 3L, 6L))          # (1,2),(1,3),(1,4): BFS from date 1
  expect_identical(cA$dropped_link_index, c(2L, 5L, 4L))         # canonical (2,3),(2,4),(3,4)
  expect_identical(cB$n_components, 2L)
  expect_identical(cB$roots, c(1, 4))
  expect_identical(cB$cycle_rank, 2L)
  expect_identical(cB$forest_link_index, c(7L, 9L, 10L, 12L))
  expect_identical(cB$dropped_link_index, c(8L, 11L))
  expect_identical(cert$dropped_index, c(2L, 4L, 5L, 8L, 11L))
  expect_identical(cA$n_outside_risk_a_units + cA$n_risk_a_units, fx$n)
  expect_identical(cA$n_pair_nonobservers, nuA$mask$n_pair_nonobservers)
  # Dropped (2,4) = -(1,2) + (1,4): sparse signed forest path.
  row24 <- which(cert$E$i == 5L)
  expect_identical(cert$E$j[row24], c(1L, 6L))
  expect_identical(cert$E$x[row24], c(-1L, 1L))

  retained <- lapply(decl[cert$retained_index], function(d) sp$fn$fit_ipw_link_shared(ses, d))
  seen <- unique(sp$counts$outcome_pairs)
  ret_pairs <- lapply(c(pairs_A, pairs_B)[cert$retained_index], as.double)
  expect_true(all(vapply(seen, function(p) any(vapply(ret_pairs, identical, logical(1L), p)), logical(1L))))
  expect_identical(sp$counts$glm, 0L)

  # Regression self-test only (not the certificate): evaluate dropped rows too.
  all_fits <- lapply(decl, function(d) cdid:::fit_ipw_link_shared(ses, d))
  expect_identical(all_fits[cert$retained_index], retained)
  cert_fits <- cdid:::certify_ipw_structural_classes(ses, all_fits, des, 1:12, "optimal", shared_bounds)
  expect_identical(cert_fits[c("retained_index", "dropped_index", "E")], cert[c("retained_index", "dropped_index", "E")])
  expect_identical(lapply(cert_fits$classes, function(z) z[-1L]), lapply(cert$classes, function(z) z[-1L]))
  est <- vapply(all_fits, `[[`, numeric(1L), "estimate")
  expect_lt(max(abs(cdid:::apply_ipw_reduction_map(cert, est[cert$retained_index]) - est)), 1e-12)
  phi <- cdid:::materialize_ipw_if(fx$ctx, all_fits, max_entries = 12 * fx$n, chunk_size = 23)$phi_matrix
  rec <- t(cdid:::apply_ipw_reduction_map(cert, t(phi[, cert$retained_index])))
  expect_lt(max(abs(rec - phi)), 1e-11)
  Wfull <- cdid:::materialize_chained_design(des, 1e6)$W
  expect_identical(cdid:::apply_ipw_reduction_map(cert, Wfull[cert$retained_index, , drop = FALSE]),
                   unname(Wfull * 1.0))

  cov_all <- cdid:::stream_ipw_covariance(fx$ctx, all_fits, 31)$covariance
  ev <- eigen(cov_all, symmetric = TRUE, only.values = TRUE)$values
  expect_identical(sum(ev / ev[1L] < 1e-12), 5L)              # total cycle rank 3 + 2
  cov_r <- cdid:::stream_ipw_covariance(fx$ctx, retained, 31)$covariance
  evr <- eigen(cov_r, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(evr) / max(evr), 1e-4)
  expect_match(cert$retained_covariance, "not_certified")

  expect_invisible(cdid:::verify_ipw_structural_certificate(cert, ses, decl, des, 1:12))
  bad <- cert
  bad$E$x[1L] <- -bad$E$x[1L]
  expect_error(cdid:::verify_ipw_structural_certificate(bad, ses, decl, des, 1:12),
               class = "cdid_ipw_certificate_mismatch_error")
  expect_error(cdid:::certify_ipw_structural_classes(ses, decl, des, 1:12, "identity", shared_bounds),
               class = "cdid_ipw_reduction_weighting_error")
})

test_that("one propensity solve per nuisance object, verified by a scoped spy; construction never reads outcomes", {
  fx <- shared_fixture()
  sp <- make_spy()
  ses <- sp$fn$construct_ipw_shared_session(fx$ctx, 3, 20, 1e5)
  nu <- sp$fn$construct_ipw_nuisance(ses, 4, c(7, Inf), 3, 6, shared_ctrl)
  expect_identical(sp$counts$glm, 1L)
  expect_identical(sp$counts$outcomes, 0L)
  decl <- lapply(pairs_B, function(p) sp$fn$declare_ipw_shared_link(ses, nu, 4, c(7, Inf), p[1L], p[2L], shared_ctrl))
  expect_identical(sp$counts$outcomes, 0L)
  fits <- lapply(decl, function(d) sp$fn$fit_ipw_link_shared(ses, d))
  expect_identical(sp$counts$glm, 1L)
  expect_identical(sp$counts$outcomes, 6L)
  for (f in fits) {
    expect_identical(f$nuisance_ref, nu$ref)
    expect_identical(f$factors$H_chol, nu$hessian$H_chol)
  }
  again <- sp$fn$fit_ipw_link_shared(ses, decl[[2L]])
  expect_identical(again, fits[[2L]])
  expect_identical(sp$counts$outcomes, 6L)
  # The ordinary path is a separate intentional path: one solve per call.
  for (p in pairs_B) sp$fn$fit_ipw_link(fx$ctx, 4, c(7, Inf), p[1L], p[2L], shared_ctrl)
  expect_identical(sp$counts$glm, 7L)
  expect_identical(environment(cdid:::ipw_glm_fit_internal), asNamespace("cdid"))

  # Outcome independence: poisoned outcomes leave the nuisance object unchanged.
  pidx_bad <- fx$pidx
  pidx_bad$observations$y[] <- NaN
  ctx_bad <- cdid:::construct_ipw_context(pidx_bad, id = seq_len(fx$n), Z = fx$Z, weights = fx$w)
  ses_bad <- cdid:::construct_ipw_shared_session(ctx_bad, 1, 2, 1e5)
  nu_bad <- cdid:::construct_ipw_nuisance(ses_bad, 4, c(7, Inf), 3, 6, shared_ctrl)
  for (fld in c("mask", "stage", "hessian", "support", "control", "K")) expect_identical(nu_bad[[fld]], nu[[fld]])
  d_bad <- cdid:::declare_ipw_shared_link(ses_bad, nu_bad, 4, c(7, Inf), 1, 2, shared_ctrl)
  expect_error(cdid:::fit_ipw_link_shared(ses_bad, d_bad), class = "cdid_ipw_outcome_error")
})

test_that("reuse refuses unequal masks, K order, control, cohort, foreign sessions/contexts and mutated objects", {
  fx <- shared_fixture()
  ses <- cdid:::construct_ipw_shared_session(fx$ctx, 3, 20, 1e5)
  nu <- cdid:::construct_ipw_nuisance(ses, 4, c(7, Inf), 1, 2, shared_ctrl)
  # Pairs touching dates 7/8 have a different risk-set mask (non-observers observe them).
  expect_error(cdid:::declare_ipw_shared_link(ses, nu, 4, c(7, Inf), 6, 7, shared_ctrl),
               class = "cdid_ipw_shared_mask_error")
  expect_error(cdid:::declare_ipw_shared_link(ses, nu, 4, c(7, Inf), 7, 8, shared_ctrl),
               class = "cdid_ipw_shared_mask_error")
  expect_error(cdid:::declare_ipw_shared_link(ses, nu, 4, c(Inf, 7), 2, 3, shared_ctrl),
               class = "cdid_ipw_shared_incompatible_error")
  expect_error(cdid:::declare_ipw_shared_link(ses, nu, 4, c(7, Inf), 2, 3,
                                              stats::glm.control(epsilon = 1e-10, maxit = 100)),
               class = "cdid_ipw_shared_incompatible_error")
  expect_error(cdid:::declare_ipw_shared_link(ses, nu, 7, c(4, Inf), 2, 3, shared_ctrl),
               class = "cdid_ipw_shared_incompatible_error")
  expect_identical(length(ses$ref$declaration_refs), 0L)

  # Foreign session/context with mutated covariates or weights.
  Z2 <- fx$Z
  Z2[5L, "x1"] <- Z2[5L, "x1"] + 1
  ctx_z <- cdid:::construct_ipw_context(fx$pidx, id = seq_len(fx$n), Z = Z2, weights = fx$w)
  ses_z <- cdid:::construct_ipw_shared_session(ctx_z, 2, 5, 1e5)
  expect_error(cdid:::declare_ipw_shared_link(ses_z, nu, 4, c(7, Inf), 2, 3, shared_ctrl),
               class = "cdid_ipw_nuisance_error")
  ctx_same <- cdid:::construct_ipw_context(fx$pidx, id = seq_len(fx$n), Z = fx$Z, weights = fx$w)
  ses_same <- cdid:::construct_ipw_shared_session(ctx_same, 2, 5, 1e5)
  expect_error(cdid:::declare_ipw_shared_link(ses_same, nu, 4, c(7, Inf), 2, 3, shared_ctrl),
               class = "cdid_ipw_nuisance_error")
  ctx_mut <- fx$ctx
  ctx_mut$Z_a[3L, 2L] <- 0
  expect_error(cdid:::construct_ipw_shared_session(ctx_mut, 1, 1, 1e5), class = "cdid_ipw_context_mismatch_error")
  ctx_mut <- fx$ctx
  ctx_mut$w_a[3L] <- 1
  expect_error(cdid:::construct_ipw_shared_session(ctx_mut, 1, 1, 1e5), class = "cdid_ipw_context_mismatch_error")
  expect_error(assign("context", ctx_z, envir = ses$ref))

  # Mutated or label-only objects are never reused.
  for (mut in list(function(x) { x$mask$observed_a_index <- x$mask$observed_a_index[-1L]; x },
                   function(x) { x$stage$beta[1L] <- x$stage$beta[1L] + 1e-9; x },
                   function(x) { x$serial <- 2L; x })) {
    expect_error(cdid:::declare_ipw_shared_link(ses, mut(nu), 4, c(7, Inf), 2, 3, shared_ctrl),
                 class = "cdid_ipw_nuisance_error")
  }
  expect_error(cdid:::declare_ipw_shared_link(ses, "nuisance_1", 4, c(7, Inf), 2, 3, shared_ctrl),
               class = "cdid_ipw_nuisance_error")
  fake <- structure(list(serial = 1L, mask = nu$mask, certified = TRUE), class = "cdid_ipw_nuisance")
  fake_env <- new.env()
  fake$ref <- fake_env
  assign("nuisance", fake, envir = fake_env)
  lockEnvironment(fake_env, bindings = TRUE)
  expect_error(cdid:::declare_ipw_shared_link(ses, fake, 4, c(7, Inf), 2, 3, shared_ctrl),
               class = "cdid_ipw_nuisance_error")

  d <- cdid:::declare_ipw_shared_link(ses, nu, 4, c(7, Inf), 2, 3, shared_ctrl, link_id = "l23")
  d_mut <- d
  d_mut$pair$t <- 4
  expect_error(cdid:::fit_ipw_link_shared(ses, d_mut), class = "cdid_ipw_declaration_error")
  expect_error(cdid:::fit_ipw_link_shared(ses_same, d), class = "cdid_ipw_declaration_error")
  f <- cdid:::fit_ipw_link_shared(ses, d)
  f_mut <- f
  f_mut$loadings$loading_beta <- -f_mut$loadings$loading_beta
  expect_error(cdid:::materialize_ipw_if(fx$ctx, f_mut, fx$n, 10), class = "cdid_ipw_shared_fit_error")
  expect_error(cdid:::replay_ipw_if(fx$ctx, list(f, f_mut), 10, function(b) NULL), class = "cdid_ipw_shared_fit_error")
  expect_error(cdid:::materialize_ipw_if(ctx_same, f, fx$n, 10), class = "cdid_ipw_context_mismatch_error")
})

test_that("equal counts, rates, covariates and coefficients with different unit IDs split; distinct numerically identical fits are never certified", {
  # Twins p_i (observe 1,2) and q_i (observe 3,4) are adjacent with identical
  # cohort, covariates and weight; core units observe all four dates.
  set.seed(4401)
  m_tw <- 30L
  n_core <- 60L
  G_core <- sample(c(3, 5, Inf), n_core, TRUE)
  G_tw <- sample(c(3, 5, Inf), m_tw, TRUE)
  x_core <- rnorm(n_core)
  x_tw <- rnorm(m_tw)
  G <- c(G_core, rep(G_tw, each = 2L))
  x <- c(x_core, rep(x_tw, each = 2L))
  n <- length(G)
  obs <- matrix(FALSE, n, 4L)
  obs[seq_len(n_core), ] <- TRUE
  tw <- n_core + seq_len(2L * m_tw)
  obs[tw[c(TRUE, FALSE)], 1:2] <- TRUE
  obs[tw[c(FALSE, TRUE)], 3:4] <- TRUE
  Y <- matrix(rnorm(n * 4L), n, 4L) + outer(x, 1:4 / 3)
  panel <- data.frame(id = rep(seq_len(n), each = 4L), time = rep(as.double(1:4), n),
                      y = as.vector(t(ifelse(obs, Y, NA_real_))), g = rep(G, each = 4L))
  pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = cbind(1, x))
  ses <- cdid:::construct_ipw_shared_session(ctx, 4, 20, 1e5)
  nu12 <- cdid:::construct_ipw_nuisance(ses, 3, c(5, Inf), 1, 2, shared_ctrl)
  nu34 <- cdid:::construct_ipw_nuisance(ses, 3, c(5, Inf), 3, 4, shared_ctrl)
  expect_identical(nu12$stage$pair_count, nu34$stage$pair_count)
  expect_identical(nu12$stage$rates, nu34$stage$rates)
  expect_identical(nu12$stage$beta, nu34$stage$beta)
  expect_identical(nu12$hessian, nu34$hessian)
  expect_false(identical(nu12$mask$observed_a_index, nu34$mask$observed_a_index))
  expect_error(cdid:::declare_ipw_shared_link(ses, nu12, 3, c(5, Inf), 3, 4, shared_ctrl),
               class = "cdid_ipw_shared_mask_error")
  d12 <- cdid:::declare_ipw_shared_link(ses, nu12, 3, c(5, Inf), 1, 2, shared_ctrl)
  d34 <- cdid:::declare_ipw_shared_link(ses, nu34, 3, c(5, Inf), 3, 4, shared_ctrl)
  f12 <- cdid:::fit_ipw_link_shared(ses, d12)
  f34 <- cdid:::fit_ipw_link_shared(ses, d34)
  phi <- cdid:::materialize_ipw_if(ctx, list(f12, f34), 2 * n, 25)$phi_matrix
  expect_gt(max(abs(phi[tw[c(TRUE, FALSE)], 1L] - phi[tw[c(FALSE, TRUE)], 1L])), 1e-3)

  # Two separately fitted but bitwise identical nuisance objects on one mask.
  fx <- shared_fixture()
  s2 <- cdid:::construct_ipw_shared_session(fx$ctx, 4, 20, 1e5)
  nuA <- cdid:::construct_ipw_nuisance(s2, 4, c(7, Inf), 1, 2, shared_ctrl)
  nuB <- cdid:::construct_ipw_nuisance(s2, 4, c(7, Inf), 1, 2, shared_ctrl)
  for (fld in c("mask", "stage", "hessian")) expect_identical(nuA[[fld]], nuB[[fld]])
  dl <- list(cdid:::declare_ipw_shared_link(s2, nuA, 4, c(7, Inf), 1, 2, shared_ctrl, functional_key = "opA"),
             cdid:::declare_ipw_shared_link(s2, nuA, 4, c(7, Inf), 2, 3, shared_ctrl, functional_key = "opA"),
             cdid:::declare_ipw_shared_link(s2, nuB, 4, c(7, Inf), 1, 3, shared_ctrl, functional_key = "opA"))
  des <- design_for(pairs_A[1:3], "opA")
  # Never one class: the cycle across distinct objects is refused, not reduced.
  expect_error(cdid:::certify_ipw_structural_classes(s2, dl, des, 1:3, "optimal", shared_bounds),
               class = "cdid_ipw_reduction_singular_error")
  cert <- cdid:::certify_ipw_structural_classes(s2, dl[c(1L, 3L)], design_for(pairs_A[c(1, 3)], "opA"), 1:2,
                                                "optimal", shared_bounds)
  expect_identical(cert$dropped_index, integer(0L))
  expect_identical(length(cert$classes), 2L)
  expect_identical(vapply(cert$classes, `[[`, integer(1L), "cycle_rank"), c(0L, 0L))
  # Ordinary refits are unclassed and always retained; their exact relations
  # along a date cycle are likewise refused rather than certified.
  ofits <- lapply(pairs_A[1:3], function(p) cdid:::fit_ipw_link(fx$ctx, 4, c(7, Inf), p[1L], p[2L], shared_ctrl))
  for (k in 1:3) {
    sf <- cdid:::fit_ipw_link_shared(s2, dl[[k]])
    expect_identical(ofits[[k]][c("estimate", "coefficients", "loadings", "factors")],
                     unclass(sf)[c("estimate", "coefficients", "loadings", "factors")])
  }
  expect_error(cdid:::certify_ipw_structural_classes(s2, ofits, des, 1:3, "optimal", shared_bounds),
               class = "cdid_ipw_reduction_singular_error")
  cert_o <- cdid:::certify_ipw_structural_classes(s2, ofits[1:2], design_for(pairs_A[1:2], "opA"), 1:2,
                                                  "optimal", shared_bounds)
  expect_identical(cert_o$retained_index, 1:2)
  expect_identical(cert_o$link_kind, rep("ordinary_unclassed", 2L))
})

test_that("intrinsic vertex-map, operator, mapping, duplicate and cross-operator cycle violations refuse", {
  fx <- shared_fixture()
  ses <- cdid:::construct_ipw_shared_session(fx$ctx, 4, 40, 1e5)
  nu <- cdid:::construct_ipw_nuisance(ses, 4, c(7, Inf), 1, 2, shared_ctrl)
  dA <- declare_all(ses, nu, 4, c(7, Inf), pairs_A[1:3], "opA")
  ok <- cdid:::certify_ipw_structural_classes(ses, dA, design_for(pairs_A[1:3], "opA"), 1:3,
                                              "optimal", shared_bounds)
  # Wrong mapping to rows with other dates/operators.
  expect_error(cdid:::certify_ipw_structural_classes(ses, dA, design_for(pairs_A[1:3], "opA"), c(2L, 1L, 3L),
                                                     "optimal", shared_bounds),
               class = "cdid_ipw_class_vertex_error")
  expect_error(cdid:::certify_ipw_structural_classes(ses, dA, design_for(pairs_A[1:3], "opB"), 1:3,
                                                     "optimal", shared_bounds),
               class = "cdid_ipw_class_operator_error")
  # Mutated design: date-1 endpoint of (1,3) re-pointed to the other operator's pre vertex.
  dA3 <- design_for(pairs_A[1:3], "opA")
  dB_pairs <- design_for(c(pairs_A[1:3], list(c(1, 2))), c(rep("opA", 3L), "opB"))
  mut <- dB_pairs
  mut$links <- mut$links[1:3, ]
  mut$incidence$n_links <- 3L
  vB <- which(mut$vertices$time == 1 & mut$vertices$functional_key %in% "opB")
  mut$links$start_vertex_index[3L] <- vB
  expect_error(cdid:::certify_ipw_structural_classes(ses, dA, mut, 1:3, "optimal", shared_bounds),
               class = "cdid_ipw_class_vertex_error")
  mut2 <- dA3
  mut2$links$end_vertex_index[1L] <- mut2$links$end_vertex_index[3L]
  expect_error(cdid:::certify_ipw_structural_classes(ses, dA, mut2, 1:3, "optimal", shared_bounds),
               class = "cdid_ipw_class_vertex_error")
  # A duplicated vertex with the same cohort/date/operator passes the intrinsic
  # endpoint check, but one date then has two vertices inside the class.
  mut3 <- dA3
  v1 <- mut3$links$start_vertex_index[1L]
  mut3$vertices <- rbind(mut3$vertices, mut3$vertices[v1, ])
  mut3$links$start_vertex_index[3L] <- nrow(mut3$vertices)
  expect_error(cdid:::certify_ipw_structural_classes(ses, dA, mut3, 1:3, "optimal", shared_bounds),
               class = "cdid_ipw_class_vertex_error")
  # Duplicate (s, t) inside one class (only constructible from a mutated design).
  dup_A <- cdid:::declare_ipw_shared_link(ses, nu, 4, c(7, Inf), 1, 2, shared_ctrl, functional_key = "opA")
  mut4 <- design_for(pairs_A[c(1, 3)], "opA")
  mut4$links$t[2L] <- 2
  mut4$links$end_vertex_index[2L] <- mut4$links$end_vertex_index[1L]
  expect_error(cdid:::certify_ipw_structural_classes(ses, list(dA[[1L]], dup_A), mut4, 1:2, "optimal", shared_bounds),
               class = "cdid_ipw_class_duplicate_error")
  # Cohort mismatch between link and design row.
  expect_error(cdid:::certify_ipw_structural_classes(ses, dA, design_for(pairs_A[1:3], "opA", g = 5, ref = 3), 1:3,
                                                     "optimal", shared_bounds),
               class = "cdid_ipw_class_vertex_error")
  dup <- cdid:::declare_ipw_shared_link(ses, nu, 4, c(7, Inf), 1, 2, shared_ctrl, functional_key = "opB")
  # Same nuisance, split across operators (1,2),(2,3) opA and (1,3) opB: the
  # retained rows contain an uncertified IF cycle -> refuse.
  d13B <- cdid:::declare_ipw_shared_link(ses, nu, 4, c(7, Inf), 1, 3, shared_ctrl, functional_key = "opB")
  expect_error(cdid:::certify_ipw_structural_classes(ses, list(dA[[1L]], dA[[2L]], d13B),
                                                     design_for(pairs_A[1:3], c("opA", "opA", "opB")), 1:3,
                                                     "optimal", shared_bounds),
               class = "cdid_ipw_reduction_singular_error")
  expect_error(cdid:::certify_ipw_structural_classes(ses, list(dA[[1L]], dup),
                                                     design_for(pairs_A[c(1, 1)], c("opA", "opB")), 1:2,
                                                     "optimal", shared_bounds),
               class = "cdid_ipw_reduction_singular_error")
  expect_identical(ok$dropped_index, 2L)
})

test_that("mean-normalization, control-sign and rate-column mutations fail references; forged registered objects fail integrity", {
  fx <- shared_fixture()
  ref <- load_potential_reference()
  ses <- cdid:::construct_ipw_shared_session(fx$ctx, 4, 20, 1e5)
  nu <- cdid:::construct_ipw_nuisance(ses, 4, c(7, Inf), 1, 2, shared_ctrl)
  d13 <- cdid:::declare_ipw_shared_link(ses, nu, 4, c(7, Inf), 1, 3, shared_ctrl)
  f13 <- cdid:::fit_ipw_link_shared(ses, d13)
  risk <- fx$G %in% c(4, 7, Inf)
  S <- as.numeric(risk & fx$obs[, 1L] & fx$obs[, 2L])
  rl <- ref$ref_link_if(ref$ref_shared_nuisance(S, fx$G, fx$Z, fx$w, 4, c(7, Inf)),
                        ifelse(S == 1, fx$Y[, 3L] - fx$Y[, 1L], 0))
  good <- cdid:::materialize_ipw_if(fx$ctx, f13, fx$n, 13)$phi
  expect_lt(max(abs(good - rl$phi)), 1e-8)

  # Unsealed ordinary copies carry each mutation through the existing replay.
  plain <- strip_shared(f13)
  plain$context_ref <- fx$ctx$ref
  class(plain) <- "cdid_ipw_link_fit"
  reload <- function(x, A) {
    x$partials$A_beta_r <- A
    x$loadings$loading_rate <- as.double(x$partials$G_r + drop(crossprod(A, x$loadings$loading_beta)))
    x$loadings$rate_coef_s0 <- -(x$masses$rates / x$masses$q) * x$loadings$loading_rate
    x$loadings$rate_coef_s1 <- ((1 - x$masses$rates) / x$masses$q) * x$loadings$loading_rate
    x
  }
  expect_identical(cdid:::materialize_ipw_if(fx$ctx, plain, fx$n, 13)$phi, good)
  norm_bad <- plain
  norm_bad$means$mean_a <- norm_bad$means$sum_a / f13$counts$n_observed_risk_units
  sign_bad <- plain
  sign_bad$means$mean_b <- -sign_bad$means$mean_b
  A0 <- plain$partials$A_beta_r
  A0[, "g"] <- 0
  rate_bad <- reload(plain, A0)
  A1 <- plain$partials$A_beta_r
  A1[, "k_1"] <- 0
  rate_bad_k <- reload(plain, A1)
  for (bad in list(norm_bad, sign_bad, rate_bad, rate_bad_k)) {
    expect_gt(max(abs(cdid:::materialize_ipw_if(fx$ctx, bad, fx$n, 13)$phi - rl$phi)), 1e-3)
  }

  # Forged objects re-sealed and registered by hand still fail recomputation.
  bnd <- shared_bounds
  forge <- function(edit) {
    x <- unclass(nu)
    x$ref <- NULL
    x <- edit(x)
    x$serial <- length(ses$ref$nuisance_refs) + 1L
    x <- cdid:::ipw_seal(structure(x, class = "cdid_ipw_nuisance"), "nuisance")
    ses$ref$nuisance_refs[[x$serial]] <- x$ref
    # Keep cache accounting consistent so this adversary reaches the separate
    # recomputation check rather than the earlier storage-integrity refusal.
    ses$ref$unit_entries <- ses$ref$unit_entries + x$unit_entries
    x
  }
  for (edit in list(function(x) { x$stage$beta[2L] <- x$stage$beta[2L] + 1e-6; x },
                    function(x) { x$stage$mean_a <- x$stage$sum_a; x },
                    function(x) { x$stage$rates[2L] <- x$stage$rates[2L] * (1 + 1e-12); x })) {
    fz <- forge(edit)
    dz <- cdid:::declare_ipw_shared_link(ses, fz, 4, c(7, Inf), 1, 2, shared_ctrl)
    expect_error(cdid:::certify_ipw_structural_classes(ses, list(dz), design_for(pairs_A[1], "opA"), 1L,
                                                       "optimal", bnd),
                 class = "cdid_ipw_nuisance_integrity_error")
  }
})

test_that("explicit allocation, cache, malformed-input and overflow guards", {
  fx <- shared_fixture()
  for (bad in list(NA, -1, 1.5, 2^53, "10", c(1, 2), Inf)) {
    expect_error(cdid:::construct_ipw_shared_session(fx$ctx, bad, 1, 1e5), class = "cdid_ipw_bound_error")
    expect_error(cdid:::construct_ipw_shared_session(fx$ctx, 1, 1, bad), class = "cdid_ipw_bound_error")
  }
  expect_error(cdid:::construct_ipw_shared_session(fx$ctx, 1, 1), class = "cdid_ipw_bound_error")
  sp <- make_spy()
  tiny <- sp$fn$construct_ipw_shared_session(fx$ctx, 2, 2, 100)
  expect_error(sp$fn$construct_ipw_nuisance(tiny, 4, c(7, Inf), 1, 2, shared_ctrl), class = "cdid_ipw_cache_bound_error")
  expect_identical(sp$counts$glm, 0L)
  ses <- cdid:::construct_ipw_shared_session(fx$ctx, 1, 2, 1e5)
  nu <- cdid:::construct_ipw_nuisance(ses, 4, c(7, Inf), 1, 2, shared_ctrl)
  expect_identical(ses$ref$unit_entries, nu$unit_entries)
  expect_error(cdid:::construct_ipw_nuisance(ses, 4, Inf, 1, 2, shared_ctrl), class = "cdid_ipw_cache_bound_error")
  d <- declare_all(ses, nu, 4, c(7, Inf), pairs_A[1:2], "opA")
  expect_error(cdid:::declare_ipw_shared_link(ses, nu, 4, c(7, Inf), 1, 3, shared_ctrl), class = "cdid_ipw_cache_bound_error")
  expect_error(cdid:::construct_ipw_nuisance(list(ref = ses$ref), 4, c(7, Inf), 1, 2, shared_ctrl),
               class = "cdid_ipw_session_error")
  broken <- ses
  broken$extra <- 1
  expect_error(cdid:::declare_ipw_shared_link(broken, nu, 4, c(7, Inf), 1, 3, shared_ctrl), class = "cdid_ipw_session_error")

  s3 <- cdid:::construct_ipw_shared_session(fx$ctx, 2, 20, 1e5)
  nu3 <- cdid:::construct_ipw_nuisance(s3, 4, c(7, Inf), 1, 2, shared_ctrl)
  dl <- declare_all(s3, nu3, 4, c(7, Inf), pairs_A, "opA")
  des <- design_for(pairs_A, "opA")
  certify <- function(b, links = dl, design = des, idx = 1:6, weighting = "optimal") {
    cdid:::certify_ipw_structural_classes(s3, links, design, idx, weighting, b)
  }
  expect_identical(length(certify(shared_bounds)$E$i), 6L)
  for (nm in names(shared_bounds)) {
    b <- shared_bounds
    b[[nm]] <- NULL
    expect_error(certify(b), class = "cdid_ipw_bound_error")
    b <- shared_bounds
    b[[nm]] <- NA
    expect_error(certify(b), class = "cdid_ipw_bound_error")
    b[[nm]] <- 2^53
    expect_error(certify(b), class = "cdid_ipw_bound_error")
  }
  expect_error(certify(c(shared_bounds, list(extra = 1))), class = "cdid_ipw_bound_error")
  tight <- list(max_links = 5, max_classes = 0, max_class_dates = 3, max_E_entries = 5, max_mask_entries = 6 * fx$n - 1)
  for (nm in names(tight)) {
    b <- shared_bounds
    b[[nm]] <- if (nm == "max_classes") 1 else tight[[nm]]
    if (nm == "max_classes") {
      expect_error(certify(b, links = c(dl, declare_all(s3, cdid:::construct_ipw_nuisance(s3, 4, Inf, 1, 2, shared_ctrl),
                                                         4, Inf, pairs_A[1], "opB")),
                           design = design_for(c(pairs_A, pairs_A[1]), c(rep("opA", 6L), "opB")), idx = 1:7),
                   class = "cdid_ipw_bound_error")
    } else {
      expect_error(certify(b), class = "cdid_ipw_bound_error")
    }
  }
  expect_error(certify(shared_bounds, weighting = c("optimal", "identity")), class = "cdid_ipw_input_error")
  expect_error(certify(shared_bounds, weighting = "pseudo_inverse"), class = "cdid_ipw_input_error")
  expect_error(certify(shared_bounds, links = list()), class = "cdid_ipw_input_error")
  expect_error(certify(shared_bounds, links = c(dl[1:5], list(list(nuisance = "nu_1", certified = TRUE)))),
               class = "cdid_ipw_input_error")
  expect_error(certify(shared_bounds, idx = c(1:5, 5L)), class = "cdid_ipw_design_error")
  expect_error(certify(shared_bounds, design = unclass(des)), class = "cdid_ipw_design_error")
  expect_error(certify(shared_bounds, design = design_for(pairs_A[1:5], "opA")), class = "cdid_ipw_design_error")
  other_ctx <- cdid:::construct_ipw_context(fx$pidx, id = seq_len(fx$n), Z = fx$Z, weights = fx$w)
  expect_error(certify(shared_bounds, links = c(dl[1:5], list(cdid:::fit_ipw_link(other_ctx, 4, c(7, Inf), 1, 4, shared_ctrl)))),
               class = "cdid_ipw_context_mismatch_error")
  cert <- certify(shared_bounds)
  expect_error(cdid:::apply_ipw_reduction_map(cert, c(1, 2)), class = "cdid_ipw_input_error")
  expect_error(cdid:::apply_ipw_reduction_map(unclass(cert), c(1, 2, 3)), class = "cdid_ipw_certificate_error")

  # Huge outcomes: the outcome-free nuisance is unchanged, replay stays finite,
  # and the existing covariance overflow guard still refuses for shared fits.
  panel <- fx$panel
  panel$y <- panel$y * 1e155
  ctx_h <- cdid:::construct_ipw_context(cdid:::construct_observed_pairs(panel, "id", "time", "y", "g"),
                                        id = seq_len(fx$n), Z = fx$Z, weights = fx$w)
  sh <- cdid:::construct_ipw_shared_session(ctx_h, 1, 2, 1e5)
  nuh <- cdid:::construct_ipw_nuisance(sh, 4, c(7, Inf), 1, 2, shared_ctrl)
  expect_identical(nuh$stage, nu$stage)
  fh <- cdid:::fit_ipw_link_shared(sh, cdid:::declare_ipw_shared_link(sh, nuh, 4, c(7, Inf), 1, 3, shared_ctrl))
  expect_true(all(is.finite(cdid:::materialize_ipw_if(ctx_h, fh, fx$n, 7)$phi)))
  expect_error(cdid:::stream_ipw_covariance(ctx_h, fh, 7), class = "cdid_ipw_numeric_error")
})

test_that("accepted two-pair fixture: shared path matches the immutable root plus ONE matched-start adaptation; its non-balanced pair refuses reuse", {
  root_path <- normalizePath(testthat::test_path("reference", "weighted_pair_ipw_derivative_reference.R"),
                             winslash = "/", mustWork = TRUE)
  root <- new.env(parent = globalenv())
  capture.output(ipw_reference_source(root_path, envir = root))
  n <- root$n
  panel <- data.frame(id = rep(seq_len(n), each = 4L), time = rep(as.double(1:4), n),
                      y = as.vector(t(ifelse(root$observed, root$Y, NA_real_))), g = rep(root$G, each = 4L))
  ctx <- cdid:::construct_ipw_context(cdid:::construct_observed_pairs(panel, "id", "time", "y", "g"),
                                      id = seq_len(n), Z = root$Z)
  ses <- cdid:::construct_ipw_shared_session(ctx, 2, 4, 1e5)
  nu12 <- cdid:::construct_ipw_nuisance(ses, 3, c(4, Inf), 1, 2, shared_ctrl)
  f12 <- cdid:::fit_ipw_link_shared(ses, cdid:::declare_ipw_shared_link(ses, nu12, 3, c(4, Inf), 1, 2, shared_ctrl))
  expect_lt(abs(f12$estimate - root$fits[[1L]]$value), 1e-14)
  phi <- cdid:::materialize_ipw_if(ctx, f12, n, 13)$phi
  expect_lt(max(abs(phi - root$phi[, 1L])), 1e-13)
  expect_lt(max(abs(phi - root$d2[, 1L])), 1e-6)
  expect_gt(max(abs(root$fits[[1L]]$wrong - root$d2[, 1L])), 0.1)
  # Pair (1,3) observes a different risk-set: reuse is refused, not refitted.
  expect_error(cdid:::declare_ipw_shared_link(ses, nu12, 3, c(4, Inf), 1, 3, shared_ctrl),
               class = "cdid_ipw_shared_mask_error")
  expect_identical(length(ses$ref$nuisance_refs), 1L)
})
