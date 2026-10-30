# Source-only deterministic reference; run from the workspace root:
# R_LIBS_USER=/opt/tools/r-libs Rscript --vanilla tools/cp02_custom_ipw_reference.R
# No package loading, installation, random draws, production expectation maker,
# or fitted whole-pipeline execution. The two-field callback is independently
# implemented below. Its composition is checked against the full equation stack
# and only then against supplied source helpers.
options(digits = 16, warn = 2)
out <- "verification"
dir.create(out, showWarnings = FALSE)
counts <- list(callback_calls = 0L, newton_calls = 0L, newton_iterations = 0L,
               glm_calls = 0L, gmm_calls = 0L, transport_calls = 0L)
checks <- data.frame(regime = character(), check = character(), error = numeric(),
                     tolerance = numeric(), pass = logical())
failures <- character()
diagnostics <- list()
save_checkpoint <- function() {
  write.csv(checks, file.path(out, "checks.csv"), row.names = FALSE)
  dput(list(counts = counts, failures = failures, diagnostics = diagnostics),
       file.path(out, "results.R"))
}
assert_close <- function(regime, label, x, y, tol = 1e-10) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)))
  err <- max(abs(x - y), 0)
  checks <<- rbind(checks, data.frame(regime = regime, check = label,
                  error = err, tolerance = tol, pass = err <= tol))
  if (err > tol) stop(paste(regime, label, "AR-EXACT failure:", err))
}
detect <- function(regime, label, x, y) {
  err <- max(abs(x - y), 0)
  checks <<- rbind(checks, data.frame(regime = regime, check = label,
                   error = err, tolerance = 1e-10, pass = err > 1e-10))
  if (!is.finite(err) || err <= 1e-10) stop(paste("Non-discriminating mutation", label))
}

# Fixture-only Newton, same 50 iterations / 1e-13 step criterion as the supplied
# shared-IPW independent Newton reference. A failure stops; there is no retry.
newton <- function(X, d, w) {
  counts$newton_calls <<- counts$newton_calls + 1L
  b <- numeric(ncol(X))
  for (iteration in 1:50) {
    prob <- stats::plogis(drop(X %*% b))
    step <- solve(crossprod(X, X * (w * prob * (1 - prob))),
                  colSums(X * (w * (d - prob))))
    b <- b + step
    counts$newton_iterations <<- counts$newton_iterations + 1L
    if (max(abs(step)) < 1e-13) return(b)
  }
  stop("Independent Newton failed at iteration 50")
}

# Complex-step differentiation of equations, not of a production IF. Holomorphic
# exp/ratios permit this without subtraction error. No finite-difference threshold
# is substituted for AR-EXACT; all point/IF/covariance checks retain 1e-10.
jacobian <- function(theta, equations) {
  vapply(seq_along(theta), function(j) {
    z <- as.complex(theta); z[j] <- z[j] + 1i * 1e-20
    Im(colMeans(equations(z))) / 1e-20
  }, numeric(length(theta)))
}

# Ordinary weighted pair callback, with all nuisance equations using i.weights.
# Its IF treats the supplied calibration weights as fixed. Return exactly the
# released panel callback's two required fields, in the input row order.
custom_ipw <- function(y1, y0, D, covariates, i.weights, inffunc, ...) {
  counts$callback_calls <<- counts$callback_calls + 1L
  stopifnot(isTRUE(inffunc))
  X <- covariates; d <- D; w <- i.weights; dy <- y1 - y0
  b <- newton(X, d, w); p <- ncol(X); n <- length(d)
  odds <- exp(drop(X %*% b))
  a <- w * d; c <- w * (1 - d) * odds
  T <- sum(a * dy) / sum(a); M <- sum(c * dy) / sum(c)
  eq <- function(th) {
    eta <- drop(X %*% th[seq_len(p)])
    cbind(X * (w * (d - 1 / (1 + exp(-eta)))),
          a * (dy - th[p + 1L]),
          w * (1 - d) * exp(eta) * (dy - th[p + 2L]))
  }
  theta <- c(b, T, M)
  xi <- t(solve(jacobian(theta, eq), -t(eq(theta))))
  list(ATT = T - M, att.inf.func = xi[, p + 1L] - xi[, p + 2L])
}

# Smallest proposed composition. All cohorts/rates and row identities belong to
# the caller. Only the callback's two fields enter this formula. This is valid
# under the weight-response identity proved in the accompanying TeX source.
compose <- function(ans, observer, gi, w, S, rate, n) {
  stopifnot(is.list(ans), is.numeric(ans$ATT), length(ans$ATT) == 1L,
            is.finite(ans$ATT), is.numeric(ans$att.inf.func),
            (is.null(dim(ans$att.inf.func)) || identical(dim(ans$att.inf.func), c(length(observer), 1L))))
  psi <- as.numeric(ans$att.inf.func)
  expected_obs <- which(as.logical(S) & !is.na(gi))
  stopifnot(length(psi) == length(observer), all(is.finite(psi)),
            length(gi) == n, length(w) == n, length(S) == n,
            all(is.finite(w) & w > 0), all(is.finite(rate) & rate > 0 & rate <= 1),
            !anyDuplicated(observer), all(observer %in% seq_len(n)),
            setequal(observer, expected_obs))
  m <- length(observer); I <- outer(gi, seq_along(rate), "==") * 1
  I[is.na(I)] <- 0
  stopifnot(m >= 2L, all(colSums(I) >= 1), all(tabulate(gi[observer], length(rate)) >= 1L))
  q <- colSums(I * w) / n
  rho <- I * w * (S - rep(rate, each = n)) / rep(q, each = n)
  sensitivity <- vapply(seq_along(rate), function(k)
    -sum(psi[gi[observer] == k]) / (m * rate[k]), numeric(1))
  base <- numeric(n); base[observer] <- n / m * psi
  list(estimate = ans$ATT, phi = base + drop(rho %*% sensitivity),
       base = base, rho = rho, sensitivity = sensitivity)
}

# Independent complete stack: rate rows, logit rows and TWO normalized mean
# equations. Differentiate and solve the joint Jacobian, without composition
# sensitivities, a supplied IF, or a supplied fit.
full_stack <- function(G, S, dy, X, w, cohorts) {
  n <- length(G); k <- length(cohorts); p <- ncol(X)
  I <- outer(G, cohorts, "==") * 1; R <- rowSums(I); D <- I[, 1L]
  obs <- which(S & R == 1)
  r <- colSums(I * w * S) / colSums(I * w)
  gi <- match(G, cohorts); rr <- rep(1, n); rr[!is.na(gi)] <- r[gi[!is.na(gi)]]
  b <- newton(X[obs, , drop = FALSE], D[obs], w[obs] / rr[obs])
  a <- w * D * S / rr; c <- w * (R - D) * S * exp(drop(X %*% b)) / rr
  T <- sum(a * dy) / sum(a); M <- sum(c * dy) / sum(c)
  eq <- function(th) {
    rate <- th[seq_len(k)]; bg <- th[k + seq_len(p)]
    rG <- rep(1 + 0 * th[1L], n)
    rG[!is.na(gi)] <- rate[gi[!is.na(gi)]]
    eta <- drop(X %*% bg); pw <- w * R * S / rG
    cbind(I * w * (S - rep(rate, each = n)),
          X * (pw * (D - 1 / (1 + exp(-eta)))),
          w * D * S / rG * (dy - th[k + p + 1L]),
          w * (R - D) * S * exp(eta) / rG * (dy - th[k + p + 2L]))
  }
  theta <- c(r, b, T, M); J <- jacobian(theta, eq)
  xi <- t(solve(J, -t(eq(theta))))
  rest <- k + seq_len(p + 2L)
  response <- -solve(J[rest, rest, drop = FALSE], J[rest, seq_len(k), drop = FALSE])
  list(estimate = T - M, phi = xi[, k + p + 1L] - xi[, k + p + 2L],
       rate_if = xi[, seq_len(k), drop = FALSE], J = J,
       sensitivity = response[p + 1L, ] - response[p + 2L, ],
       residual = max(abs(colMeans(eq(theta)))), beta = b)
}

main <- function() {
  # Source actual supplied definitions into a distinct comparator environment.
  # No mocks: cdid_abort is supplied by forward.R. Source-time definitions do
  # not call did/DRDID or load the native DLL. Only conventional base-R helpers
  # are invoked. In particular no whole public/package fit is attempted.
  prod <- new.env(parent = globalenv())
  for (f in c("forward", "chained_pairs", "chained_design", "chained_comparison_catalog",
              "chained_ipw_link", "chained_ipw_shared", "chained_gmm_core",
              "chained_ipw_pipeline", "chained_public_control"))
    sys.source(file.path("cdid/R", paste0(f, ".R")), envir = prod)
  # Count actual glm calls at the supplied solver's single call site, without
  # changing its body or any inputs. The tracing expression is diagnostic only.
  trace("glm.fit", tracer = quote(counts$glm_calls <<- counts$glm_calls + 1L),
        print = FALSE, where = asNamespace("stats"))
  on.exit(untrace("glm.fit", where = asNamespace("stats")), add = TRUE)
  panel <- read.csv("reference/fixtures/panel24.csv")
  ids <- sort(unique(panel$id)); n <- length(ids); calendar <- sort(unique(panel$time))
  first <- match(ids, panel$id); G0 <- panel$g[first]; sign <- panel$sign[first]
  Y <- matrix(NA_real_, n, length(calendar))
  Y[cbind(match(panel$id, ids), match(panel$time, calendar))] <- panel$y
  observed <- !is.na(Y)
  # Fixed, declared variants, chosen before fitting. The supplied CSV is never
  # edited. NYT maps original control IDs 4,8,16 to cohort 4, providing two
  # supported control cohorts and D13's smaller fixed-horizon pools. All outcome
  # values and missingness stay exactly as supplied. Weights are deterministic
  # positive mechanical case weights; this is not a survey-sampling validation.
  regimes <- list(
    NT_unit = list(G = G0, w = rep(1, n), X = cbind(intercept = rep(1, n)), cg = "nevertreated"),
    NT_weighted_covariate = list(G = G0, w = 1 + (ids %% 5) / 4,
        X = cbind(intercept = 1, sign = sign), cg = "nevertreated"),
    NYT_unit = list(G = replace(G0, ids %in% c(4, 8, 16), 4),
        w = rep(1, n), X = cbind(intercept = rep(1, n)), cg = "notyettreated"),
    NYT_weighted = list(G = replace(G0, ids %in% c(4, 8, 16), 4),
        w = 1 + (ids %% 5) / 4, X = cbind(intercept = rep(1, n)), cg = "notyettreated"))
  for (regime in names(regimes)) {
    fx <- regimes[[regime]]; G <- fx$G; w <- fx$w; X <- fx$X
    dat <- panel; dat$g <- G[match(dat$id, ids)]
    pair <- prod$construct_observed_pairs(dat, "id", "time", "y", "g")
    ctx <- prod$construct_ipw_context(pair, ids, X, w)
    stopifnot(identical(ctx$a_id, ids), ctx$n_a_units == n)
    # Enumerate support independently, then compare all D10/D13 catalog rows.
    occ <- do.call(rbind, lapply(combn(calendar, 2, simplify = FALSE), function(e)
      data.frame(cohort = sort(unique(G)), s = e[1], t = e[2], count = vapply(sort(unique(G)),
        function(g) sum(G == g & observed[, e[1]] & observed[, e[2]]), numeric(1)))))
    cat <- prod$build_comparison_catalog(calendar, sort(unique(G)), 3, 0, occ,
      fx$cg, "fixed-baseline-reference", NULL,
      list(max_candidates = 1000, max_rows = 1000, max_pool_entries = 10000,
           max_scratch = 1e7, max_work = 1e8))
    rows <- cat$rows
    pool <- function(j) cat$spec$cohorts[cat$pools$entries[seq.int(
      cat$pools$start[j], length.out = cat$pools$size[j])]]
    keys <- vapply(seq_len(nrow(rows)), function(j)
      paste(rows$s[j], rows$t[j], paste(pool(rows$pool_id[j]), collapse = ":")), character(1))
    expected_keys <- character()
    for (e in combn(calendar, 2, simplify = FALSE)) {
      S <- observed[, e[1]] & observed[, e[2]]
      if (!any(G == 3 & S)) next
      available <- sort(unique(G[G != 3 & S]))
      horizons <- if (fx$cg == "notyettreated") unique(c(e[2], calendar[calendar >= max(3, e[2])])) else e[2]
      for (h in horizons) {
        K <- available[is.infinite(available) | (fx$cg == "notyettreated" & available > 3 & available > h)]
        if (length(K)) expected_keys <- c(expected_keys, paste(e[1], e[2], paste(K, collapse = ":")))
      }
    }
    stopifnot(setequal(keys, unique(expected_keys)), !anyDuplicated(keys))
    L <- nrow(rows); composed <- stack <- builtin <- matrix(0, n, L)
    mus <- refmu <- builtinmu <- numeric(L); max_rate <- 0; raw_point_errors <- numeric()
    rate_composed <- rate_stack <- vector("list", L)
    z_perm <- NULL
    for (j in seq_len(L)) {
      K <- pool(rows$pool_id[j]); cohorts <- c(3, K); gi <- match(G, cohorts)
      s <- rows$s[j]; t <- rows$t[j]; S <- observed[, s] & observed[, t]
      obs <- which(S & !is.na(gi)); dy <- numeric(n); dy[obs] <- Y[obs, t] - Y[obs, s]
      I <- outer(G, cohorts, "==") * 1
      rate <- colSums(I * w * S) / colSums(I * w)
      # Deliberately permute pair input order; caller keeps the exact ID map.
      obs <- rev(obs)
      ans <- custom_ipw(Y[obs, t], Y[obs, s], as.numeric(G[obs] == 3),
                        X[obs, , drop = FALSE], w[obs] / rate[gi[obs]], TRUE)
      z <- compose(ans, obs, gi, w, S, rate, n)
      ref <- full_stack(G, S, dy, X, w, cohorts)
      fit <- prod$fit_ipw_link(ctx, 3, K, s, t,
        stats::glm.control(epsilon = 1e-12, maxit = 50, trace = FALSE))
      # Stream production IF only for comparison. Independent expectations above
      # are already calculated before these helper calls.
      prod$replay_ipw_if(ctx, fit, 7, function(batch) builtin[batch$a_index, j] <<- batch$phi[, 1L])
      composed[, j] <- z$phi; stack[, j] <- ref$phi
      mus[j] <- z$estimate; refmu[j] <- ref$estimate; builtinmu[j] <- fit$estimate
      assert_close(regime, paste0("rate IF ", j), z$rho, ref$rate_if)
      assert_close(regime, paste0("rate sensitivity ", j), z$sensitivity, ref$sensitivity)
      rate_composed[[j]] <- z$rho; rate_stack[[j]] <- ref$rate_if
      assert_close(regime, paste0("inverse-rate cohort mass ", j), colSums(I * w * S / rep(
        replace(rep(1, n), !is.na(gi), rate[gi[!is.na(gi)]]), length(cohorts))) / n,
        colSums(I * w) / n)
      max_rate <- max(max_rate, abs(z$phi - z$base))
      if (j == 1L && fx$cg == "notyettreated") {
        ans_perm <- ans
        in_g <- which(G[obs] == 3)
        ans_perm$att.inf.func[in_g] <- rev(ans_perm$att.inf.func[in_g])
        z_perm <- compose(ans_perm, obs, gi, w, S, rate, n)
        assert_close(regime, "within-cohort permutation preserves rate sensitivity", z_perm$sensitivity, ref$sensitivity)
        detect(regime, "within-cohort observer permutation IF", z_perm$phi, ref$phi)
      }
      if (fx$cg == "notyettreated" && length(K) > 1L) {
        detect(regime, paste0("omit estimated-rate IF ", j), z$base, ref$phi)
        raw <- custom_ipw(Y[obs, t], Y[obs, s], as.numeric(G[obs] == 3), X[obs, , drop = FALSE], w[obs], TRUE)
        raw_point_errors <- c(raw_point_errors, abs(raw$ATT - ref$estimate))
      }
    }
    assert_close(regime, "point vs full stack", mus, refmu)
    assert_close(regime, "point vs supplied helper", mus, builtinmu)
    assert_close(regime, "complete joint IF vs full stack", composed, stack)
    assert_close(regime, "complete joint IF vs supplied helper", composed, builtin)
    if (length(raw_point_errors)) detect(regime, "omit point calibration across pooled rows", raw_point_errors, rep(0, length(raw_point_errors)))
    omega <- function(phi) crossprod(sweep(phi, 2, colMeans(phi))) / nrow(phi)
    assert_close(regime, "full joint rate/moment covariance", omega(cbind(composed, do.call(cbind, rate_composed))),
      omega(cbind(stack, do.call(cbind, rate_stack))))
    V <- omega(stack); Vprod <- omega(builtin)
    assert_close(regime, "joint covariance vs supplied helper", omega(composed), Vprod)
    detect(regime, "wrong row identity", replace(composed, seq_len(n), rev(composed[, 1L])), stack)
    wrong <- composed; wrong[, 1L] <- rev(wrong[, 1L])
    detect(regime, "wrong row cross-covariance", omega(wrong), V)
    if (fx$cg == "notyettreated") {
      wrong_within <- composed; wrong_within[, 1L] <- z_perm$phi
      detect(regime, "within-cohort observer permutation cross-covariance", omega(wrong_within), V)
    }
    # Nontrivial full-frame map, two padding rows. These are deterministic
    # transport probes, not added fitting units or new population observations.
    N <- n + 2L; frame_ids <- c(-101, rev(ids), -102); map <- match(ids, frame_ids)
    embedded <- matrix(0, N, L); embedded[map, ] <- N / n * composed
    counts$transport_calls <<- counts$transport_calls + 1L
    transported <- prod$transport_if_to_full_frame(stack, N, map)$Phi_full
    assert_close(regime, "sample-root IF centering", c(colMeans(composed), colMeans(embedded)), rep(0, 2L * L))
    assert_close(regime, "full-frame transport", embedded, transported)
    assert_close(regime, "full-frame variance normalization", omega(embedded) / N, V / n)
    detect(regime, "omit full-frame scale", embedded * n / N, transported)
    # Same complete D10/D13 system, all rows retained. The supplied conventional
    # source GMM path is supported without DLL; no native/installed call is made.
    dq <- data.frame(g = 3, t = 3:5)
    design <- prod$build_comparison_catalog_design(cat, dq, 3:5)
    W <- prod$materialize_chained_design(design, 1e5)$W
    stopifnot(nrow(W) == L, ncol(W) > 0L, qr(W)$rank == ncol(W))
    # Independent saddle-point equations, no matrix inverse. With PD V and
    # identified W this computes exactly the same fixed-system GLS map.
    p <- ncol(W); KKT <- rbind(cbind(V, W), cbind(t(W), matrix(0, p, p)))
    sol <- tryCatch(solve(KKT, rbind(diag(L), matrix(0, p, L))), error = identity)
    prov <- prod$gmm_covariance_provenance("unit", n, n, n, TRUE, FALSE, "none", 1 / n,
                                        "baseline-reference", "complete-stack")
    replay <- prod$gmm_if_replay_spec(function(consume_batch) consume_batch(list(unit_id = ids, phi = builtin)),
      "unit", FALSE, n, n, n, TRUE, FALSE, "none", 1 / n, "baseline-reference", "complete-stack",
      expected_unit_ids = ids)
    carried <- prod$chained_public_carried_settings()$gmm_guard
    accuracy <- do.call(prod$gmm_accuracy_control, c(list(method = "conventional_chol_qr",
      target_error_unit = "abs", max_target_error_abs = prod$chained_public_staging_envelope(),
      max_target_error_se = NA_real_), carried, list(require_matching_scale = TRUE,
      max_loading_error_abs = prod$chained_public_staging_envelope())))
    bounds <- prod$gmm_bounds_control(1e7, n, 100, n, TRUE, 1e5)
    C <- diag(p); colnames(C) <- paste0("parameter_", seq_len(p))
    counts$gmm_calls <<- counts$gmm_calls + 1L
    fitg <- tryCatch(prod$solve_chained_gmm_core(design, mus, Vprod, prov, C, "optimal", replay,
      accuracy, bounds, "lowest", request_param_cov = TRUE, materialize_target_if = TRUE), error = identity)
    gmm_result <- list()
    if (inherits(sol, "error") || inherits(fitg, "error")) {
      # Keep genuine singular/conditioning refusals. No alternate weighting,
      # moment deletion, pinv, ridge, or more permissive criterion is tried.
      gmm_result <- list(status = "unresolved finite-fixture covariance/solve refusal (not exact/population singularity)",
        independent_error = if (inherits(sol, "error")) conditionMessage(sol) else NULL,
        comparator_error = if (inherits(fitg, "error")) conditionMessage(fitg) else NULL)
    } else {
    B <- sol[L + seq_len(p), , drop = FALSE]
    theta <- drop(B %*% refmu); psi <- stack %*% t(B)
    Sigma <- omega(psi) / n
    assert_close(regime, "GMM point", theta, fitg$targets$estimate)
    assert_close(regime, "GMM joint IF", psi, fitg$if_projection$materialized_target_if)
    # With normalization_factor=1/n the private core returns the covariance
    # of the sqrt(n) linear term. Sampling SE divides it by sqrt(n), as in the
    # pipeline's downstream inference. Do not compare these different scales.
    assert_close(regime, "GMM sampling SE", sqrt(diag(Sigma)), fitg$targets$se_sandwich / sqrt(n))
    assert_close(regime, "GMM core asymptotic covariance", Sigma * n, fitg$if_projection$sandwich_cov)
    assert_close(regime, "GMM covariance", Sigma,
      crossprod(sweep(fitg$if_projection$materialized_target_if, 2,
                      colMeans(fitg$if_projection$materialized_target_if))) / n^2)
    zcrit <- stats::qnorm(.975)
    assert_close(regime, "fixed pointwise CI", cbind(theta - zcrit * sqrt(diag(Sigma)), theta + zcrit * sqrt(diag(Sigma))),
      cbind(fitg$targets$estimate - zcrit * fitg$targets$se_sandwich / sqrt(n),
            fitg$targets$estimate + zcrit * fitg$targets$se_sandwich / sqrt(n)))
    gmm_result <- list(status = "fixed-system source comparison passed", points = theta, SE = sqrt(diag(Sigma)))
    }
    diagnostics[[regime]] <<- list(n_units = n, n_links = L, n_original = sum(rows$original),
      n_horizon_only = sum(!rows$original), pair_keys = keys, max_rate_IF = max_rate,
      rcond_covariance = rcond(V), kappa_KKT = kappa(KKT),
      max_abs_shared_covariance = max(abs(V[row(V) != col(V)])),
      GMM = gmm_result)
    save_checkpoint()
  }
  # Focused pooled-NYT non-saturated covariate edge check: verifies nonzero
  # indirect score-rate sensitivity L' H^{-1} G_beta, nonzero treated-rate
  # sensitivity b_{ell g}, failure of unweighted-propensity callbacks under
  # point calibration and rate sensitivity, and exact shared-nuisance IF alignment.
  G_e <- replace(G0, ids %in% c(4, 8, 16), 4)
  w_e <- 1 + (ids %% 5) / 4
  X_e <- cbind(intercept = 1, x = sign + (ids %% 3 - 1) / 4)
  dat_e <- panel; dat_e$g <- G_e[match(dat_e$id, ids)]
  pair_e <- prod$construct_observed_pairs(dat_e, "id", "time", "y", "g")
  ctx_e <- prod$construct_ipw_context(pair_e, ids, X_e, w_e)
  K_e <- c(4, Inf); cohorts_e <- c(3, K_e); gi_e <- match(G_e, cohorts_e)
  S_e <- observed[, 1] & observed[, 2]
  obs_e <- rev(which(S_e & !is.na(gi_e)))
  dy_e <- numeric(n); dy_e[obs_e] <- Y[obs_e, 2] - Y[obs_e, 1]
  I_e <- outer(G_e, cohorts_e, "==") * 1
  rate_e <- colSums(I_e * w_e * S_e) / colSums(I_e * w_e)
  ans_e <- custom_ipw(Y[obs_e, 2], Y[obs_e, 1], as.numeric(G_e[obs_e] == 3),
                      X_e[obs_e, , drop = FALSE], w_e[obs_e] / rate_e[gi_e[obs_e]], TRUE)
  z_e <- compose(ans_e, obs_e, gi_e, w_e, S_e, rate_e, n)
  z_mat <- compose(list(ATT = ans_e$ATT, att.inf.func = matrix(ans_e$att.inf.func, ncol = 1L)),
                   obs_e, gi_e, w_e, S_e, rate_e, n)
  ref_e <- full_stack(G_e, S_e, dy_e, X_e, w_e, cohorts_e)
  ctrl_e <- stats::glm.control(epsilon = 1e-12, maxit = 50, trace = FALSE)
  fit_e <- prod$fit_ipw_link(ctx_e, 3, K_e, 1, 2, ctrl_e)
  phi_builtin_e <- numeric(n)
  prod$replay_ipw_if(ctx_e, fit_e, 7, function(batch) phi_builtin_e[batch$a_index] <<- batch$phi[, 1L])
  sess_e <- prod$construct_ipw_shared_session(ctx_e, max_nuisance = 4, max_links = 4, max_unit_entries = 1000)
  nu_e <- prod$construct_ipw_nuisance(sess_e, 3, K_e, 1, 2, ctrl_e)
  decl_e <- prod$declare_ipw_shared_link(sess_e, nu_e, 3, K_e, 1, 2, ctrl_e)
  fit_sh_e <- prod$fit_ipw_link_shared(sess_e, decl_e)
  phi_shared_e <- numeric(n)
  prod$replay_ipw_if(ctx_e, fit_sh_e, 7, function(batch) phi_shared_e[batch$a_index] <<- batch$phi[, 1L])
  mask_err <- tryCatch({prod$declare_ipw_shared_link(sess_e, nu_e, 3, K_e, 2, 3, ctrl_e); NULL}, error = identity)
  stopifnot(inherits(mask_err, "cdid_ipw_shared_mask_error"))
  assert_close("NYT_pooled_covariate", "point vs full stack", z_e$estimate, ref_e$estimate)
  assert_close("NYT_pooled_covariate", "point vs shared helper", z_e$estimate, fit_sh_e$estimate)
  assert_close("NYT_pooled_covariate", "rate sensitivity", z_e$sensitivity, ref_e$sensitivity)
  assert_close("NYT_pooled_covariate", "complete IF vs full stack", z_e$phi, ref_e$phi)
  assert_close("NYT_pooled_covariate", "complete IF vs direct helper", z_e$phi, phi_builtin_e)
  assert_close("NYT_pooled_covariate", "complete IF vs shared helper", z_e$phi, phi_shared_e)
  assert_close("NYT_pooled_covariate", "1-column matrix IF acceptance", z_mat$phi, z_e$phi)
  k_e <- length(cohorts_e); p_e <- ncol(X_e)
  G_r_sens <- -ref_e$J[k_e + p_e + 1L, seq_len(k_e)] / ref_e$J[k_e + p_e + 1L, k_e + p_e + 1L] +
               ref_e$J[k_e + p_e + 2L, seq_len(k_e)] / ref_e$J[k_e + p_e + 2L, k_e + p_e + 2L]
  indirect_sens <- -drop((ref_e$J[k_e + p_e + 2L, k_e + seq_len(p_e)] / ref_e$J[k_e + p_e + 2L, k_e + p_e + 2L]) %*%
                    solve(ref_e$J[k_e + seq_len(p_e), k_e + seq_len(p_e)], ref_e$J[k_e + seq_len(p_e), seq_len(k_e)]))
  assert_close("NYT_pooled_covariate", "direct + indirect rate sensitivity decomposition",
               G_r_sens + indirect_sens, ref_e$sensitivity)
  detect("NYT_pooled_covariate", "omit indirect propensity-rate sensitivity", G_r_sens, ref_e$sensitivity)
  detect("NYT_pooled_covariate", "omit treated-rate sensitivity",
         c(0, z_e$sensitivity[-1L]), ref_e$sensitivity)
  b_unw <- newton(X_e[obs_e, , drop = FALSE], as.numeric(G_e[obs_e] == 3), rep(1, length(obs_e)))
  w_cal <- w_e[obs_e] / rate_e[gi_e[obs_e]]; d_e <- as.numeric(G_e[obs_e] == 3); dy_obs <- Y[obs_e, 2] - Y[obs_e, 1]
  a_u <- w_cal * d_e; c_u <- w_cal * (1 - d_e) * exp(drop(X_e[obs_e, , drop = FALSE] %*% b_unw))
  T_u <- sum(a_u * dy_obs) / sum(a_u); M_u <- sum(c_u * dy_obs) / sum(c_u)
  eq_u <- function(th) {
    eta <- drop(X_e[obs_e, , drop = FALSE] %*% th[seq_len(p_e)])
    cbind(X_e[obs_e, , drop = FALSE] * (d_e - 1 / (1 + exp(-eta))),
          a_u * (dy_obs - th[p_e + 1L]),
          w_cal * (1 - d_e) * exp(eta) * (dy_obs - th[p_e + 2L]))
  }
  xi_u <- t(solve(jacobian(c(b_unw, T_u, M_u), eq_u), -t(eq_u(c(b_unw, T_u, M_u)))))
  z_u <- compose(list(ATT = T_u - M_u, att.inf.func = xi_u[, p_e + 1L] - xi_u[, p_e + 2L]),
                 obs_e, gi_e, w_e, S_e, rate_e, n)
  detect("NYT_pooled_covariate", "unweighted-propensity point calibration failure", z_u$estimate, ref_e$estimate)
  detect("NYT_pooled_covariate", "unweighted-propensity rate-sensitivity failure", z_u$sensitivity, ref_e$sensitivity)
  diagnostics$NYT_pooled_covariate <<- list(
    estimate = z_e$estimate, sensitivity = z_e$sensitivity, G_r_only_sensitivity = G_r_sens,
    indirect_propensity_rate_sensitivity = indirect_sens,
    shared_mask_mismatch_refused = conditionMessage(mask_err))
  # Finite/type/shape/pool-membership rejection probes. Field shape cannot
  # validate semantics or discover a dishonest unlabelled row permutation.
  good <- list(ATT = 0, att.inf.func = c(0, 0))
  rejected <- 0L
  for (bad in list(list(ATT = NA_real_, att.inf.func = c(0, 0)),
                   list(ATT = 0, att.inf.func = 0),
                   list(ATT = 0, att.inf.func = c(Inf, 0)))) {
    err <- tryCatch({compose(bad, 1:2, c(1, 2), c(1, 1), c(1, 1), c(1, 1), 2); NULL}, error = identity)
    stopifnot(inherits(err, "error")); rejected <- rejected + 1L
  }
  for (bad_call in list(
    quote(compose(good, 1:2, c(1, 1, 2), c(1, 1, 1), c(1, 1, 0), c(1, 0.5), 3)),
    quote(compose(list(ATT = 0, att.inf.func = c(0, 0, 0)), 1:3, c(1, 2, NA), c(1, 1, 1), c(1, 1, 1), c(1, 1), 3))
  )) {
    err <- tryCatch({eval(bad_call); NULL}, error = identity)
    stopifnot(inherits(err, "error")); rejected <- rejected + 1L
  }
  diagnostics$negative_outputs <<- list(rejected = rejected,
    kinds = c("nonfinite ATT", "wrong IF length", "nonfinite IF",
              "zero observed-pair cohort membership", "nonobserver included in observer"))
}
tryCatch({main(); save_checkpoint(); cat("PASS:", nrow(checks), "checks; AR-EXACT tolerance 1e-10.\n"); print(counts)},
 error = function(e) {failures <<- c(failures, conditionMessage(e)); save_checkpoint();
   message("FAIL (retained): ", conditionMessage(e)); quit(status = 1L)})
