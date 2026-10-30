# Probability clipping of chained inverse-probability weights (did 2.5.1 /
# DRDID 1.3.0 epsilon 1e-6): treatment propensities capped at 1 - 1e-6,
# inverse observation denominators (fitted pair probabilities, final date
# products, positive empirical pair rates) floored at 1e-6.
#
# The oracles below are written from the clipped functional itself:
# closed-form roots of saturated binary-covariate logit cells and
# stats::glm.fit for the pooled date model. They never call the package
# clipping helpers. Package points and complete influence functions are
# compared with these oracles and with complete-refit case-weight central
# differences, on both sides of each kink. Each finite difference perturbs
# weights well inside one side of its cutoff. A value exactly at a cutoff is
# a kink, so the package convention there (counted as capped or floored, with
# zero derivative) is checked on the helpers only. Finite fixtures check
# algebra and wiring, not coverage.

clip_eps <- 1e-6
clip_cap <- 1 - 1e-6
clip_cap_odds <- (1 - 1e-6) / (1 - (1 - 1e-6))

clip_skip <- function() {
  testthat::skip_if_not(isNamespaceLoaded("cdid") || requireNamespace("cdid", quietly = TRUE),
    "cdid namespace required")
}

clip_quasi <- function(X, y, w) {
  stats::glm.fit(X, y, weights = w, family = stats::quasibinomial(),
    control = stats::glm.control(epsilon = 1e-14, maxit = 100))$coefficients
}

# ------------------------------------------------------------------------------
# Helper contract
# ------------------------------------------------------------------------------

test_that("clip constants, stable raw probabilities, capped odds and strict masks", {
  clip_skip()
  expect_identical(cdid:::ipw_clip_eps, clip_eps)
  expect_identical(cdid:::ipw_clip_cap, clip_cap)
  # Exact binary64 cap odds cap / (1 - cap); 1 - cap is exact here, and cap / 1e-6 differs.
  expect_identical(cdid:::ipw_clip_cap_odds, clip_cap / (1 - clip_cap))
  expect_false(identical(cdid:::ipw_clip_cap_odds, clip_cap / 1e-6))
  expect_lt(abs((clip_cap / 1e-6) / cdid:::ipw_clip_cap_odds - 1 - 2.8755660985637216e-11), 1e-14)

  eta <- c(-800, -40, -2, 0, 3, 13.8, 13.82, 25, 36.8, 40, 40.5, 700, 800)
  lv <- cdid:::ipw_clip_logit(eta)
  lo <- eta <= 40
  e <- exp(eta[lo])
  expect_identical(lv$prob[lo], e / (1 + e))             # original arithmetic
  expect_identical(lv$prob[!lo], rep(1, sum(!lo)))       # stable beyond 40
  p_exp <- c(e / (1 + e), rep(1, sum(!lo)))
  act <- p_exp < clip_cap
  expect_identical(lv$active, as.double(act))
  expect_identical(lv$odds[act], exp(eta[act]))
  expect_true(all(lv$odds[!act] == clip_cap_odds))
  expect_true(all(is.finite(lv$odds)))                   # no overflow at eta = 800
  expect_identical(cdid:::ipw_logit_values(cbind(1, eta), c(0, 1))$odds, lv$odds)
  nf <- cdid:::ipw_clip_logit(c(Inf, NaN, NA_real_))
  expect_identical(nf$odds, c(Inf, NaN, NA_real_))        # preserve refusal values, not NA coercion
  expect_true(all(is.na(nf$prob)))
  neginf <- cdid:::ipw_clip_logit(-Inf)
  expect_identical(neginf$eta, -Inf)                    # raw predictor guard still sees -Inf
  expect_identical(neginf$odds, 0)

  # Strict masks: equality is a kink, counted as capped or floored.
  expect_identical(cdid:::ipw_clip_odds_active(c(clip_cap - 2^-53, clip_cap, 1)), c(1, 0, 0))
  q <- c(0, 5e-7, clip_eps, clip_eps * (1 + 2^-50), 2e-6, 1)
  expect_identical(cdid:::ipw_clip_floor(q), c(clip_eps, clip_eps, clip_eps, q[4:6]))
  expect_identical(cdid:::ipw_clip_floor_active(q), c(0, 0, 0, 1, 1, 1))
})

# ------------------------------------------------------------------------------
# CP02 baseline link: capped control odds and floored empirical rates
# ------------------------------------------------------------------------------

# Saturated binary-x treatment logit on a 1-2 pair; cohort 4 carries the rate
# kink and the x = 1 cell the treatment-cap kink. "bound": r_4 ~ 5e-8 and
# p(x = 1) ~ 1 - 5.5e-9; "active": r_4 ~ 2e-6 and p(x = 1) ~ 1 - 2.0e-6.
clip_cp02_fixture <- function(side) {
  b <- identical(side, "bound")
  u <- data.frame(
    id = 1:18,
    g = c(3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, Inf, Inf, Inf, Inf, Inf, Inf, 5),
    x = c(0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0),
    w = c(1.0, 1.5, 0.7, 1.2, 0.8, 1.1, 1.0,
      if (b) 1e-7 else 4e-6, if (b) 1e-14 else 4e-12, 1.0, 1.0,
      1.0, 2.0, 1.0, if (b) 1e-9 else 1e-6, if (b) 2e-9 else 1e-6, 1.5, 1.0),
    pattern = c("123", "12", "13", "123", "12", "23", "2", "12", "12", "13", "23",
      "123", "12", "13", "12", "123", "23", "123"),
    stringsAsFactors = FALSE)
  times <- 1:3
  obs <- unname(t(vapply(u$pattern, function(p) vapply(times, function(s) grepl(as.character(s), p), logical(1)),
    logical(3))))
  Y <- outer(u$id, times, function(i, s) round(cos(1.7 * i + s) + 0.25 * s * (u$g[i] == 3), 6))
  panel <- data.frame(id = rep(u$id, each = 3L), time = rep(as.double(times), nrow(u)),
    y = as.vector(t(ifelse(obs, Y, NA_real_))), g = rep(u$g, each = 3L))
  list(units = u, obs = obs, Y = Y, panel = panel)
}

clip_cp02_oracle <- function(fx, cw = rep(1, nrow(fx$units))) {
  u <- fx$units
  w <- u$w * cw
  A <- rowSums(fx$obs) >= 2
  S <- fx$obs[, 1] & fx$obs[, 2]
  coh <- c(3, 4, Inf)
  k <- match(u$g, coh)
  r <- vapply(seq_along(coh), function(j) sum(w[A & S & k %in% j]) / sum(w[A & k %in% j]), numeric(1))
  rc <- pmax(r, clip_eps)
  risk <- A & S & !is.na(k)
  fw <- numeric(nrow(u)); fw[risk] <- w[risk] / rc[k[risk]]
  D <- as.double(u$g == 3)
  dy <- ifelse(S, fx$Y[, 2] - fx$Y[, 1], 0)
  p <- cell_odds <- numeric(2)
  for (v in 0:1) {
    mt <- sum(fw[risk & u$x == v & D == 1]); mc <- sum(fw[risk & u$x == v & D == 0])
    p[v + 1] <- mt / (mt + mc)
    cell_odds[v + 1] <- if (p[v + 1] < clip_cap) mt / mc else clip_cap_odds
  }
  odds <- cell_odds[u$x + 1]
  a <- fw * D * risk; bw <- fw * (1 - D) * risk * odds
  list(mu = sum(a * dy) / sum(a) - sum(bw * dy) / sum(bw), r = r, rc = rc, p = p, fw = fw,
    b = bw, D = D, dy = dy, risk = risk, k = k, w = w, C = sum(bw * dy) / sum(bw), sum_b = sum(bw))
}

clip_cp02_fit <- function(fx) {
  u <- fx$units
  pidx <- cdid:::construct_observed_pairs(fx$panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(pidx, id = u$id, Z = cbind(intercept = 1, x = u$x), weights = u$w)
  fit <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 2,
    control = stats::glm.control(epsilon = 1e-12, maxit = 100), link_id = "L12")
  list(ctx = ctx, fit = fit)
}

clip_cp02_fd <- function(fx, ids, n_a, h = 1e-6) {
  pos <- match(ids, fx$units$id)
  vapply(pos, function(i) {
    cp <- cm <- rep(1, nrow(fx$units)); cp[i] <- 1 + h; cm[i] <- 1 - h
    n_a * (clip_cp02_oracle(fx, cp)$mu - clip_cp02_oracle(fx, cm)$mu) / (2 * h)
  }, numeric(1))
}

# Rebuild the replay loadings of a fit from substituted partials.
clip_cp02_reload <- function(fit, G_beta, G_r, A_beta_r) {
  f <- fit
  lb <- drop(solve(fit$factors$H, G_beta))
  lr <- as.double(G_r + drop(crossprod(A_beta_r, lb)))
  f$loadings$loading_beta <- lb
  f$loadings$loading_rate <- lr
  f$loadings$rate_coef_s0 <- -(fit$masses$rates / fit$masses$q) * lr
  f$loadings$rate_coef_s1 <- ((1 - fit$masses$rates) / fit$masses$q) * lr
  f
}

test_that("CP02 link: capped control odds and a floored rate match the clipped oracle, IF and mutations", {
  clip_skip()
  fx <- clip_cp02_fixture("bound")
  run <- clip_cp02_fit(fx); ctx <- run$ctx; fit <- run$fit
  or <- clip_cp02_oracle(fx)
  expect_lt(or$r[2], clip_eps)                       # floored empirical rate
  expect_gt(or$p[2], clip_cap)                       # capped treatment cell
  expect_lt(or$p[2], 1 - 1e-12)
  expect_lt(max(abs(unname(fit$masses$rates) / or$r - 1)), 1e-12)   # raw rates kept
  expect_lt(abs(fit$estimate - or$mu), 1e-8)

  clip <- fit$diagnostics$clipping
  expect_identical(unlist(clip[c("p_capped", "p_capped_ctrl", "r_floored")], use.names = FALSE), c(5, 3, 1))
  expect_identical(clip[["r_min"]], min(unname(fit$masses$rates)))

  n_a <- ctx$n_a_units
  mat <- cdid:::materialize_ipw_if(ctx, fit, max_entries = n_a, chunk_size = 5L, include_components = TRUE)
  fd <- clip_cp02_fd(fx, mat$id, n_a)
  expect_lt(max(abs(mat$phi - fd)), 1e-6)
  expect_lt(abs(mean(mat$phi)), 1e-9)
  out <- mat$id == 18
  expect_identical(max(abs(mat$phi[out])), 0)       # outside-risk A unit
  one <- cdid:::materialize_ipw_if(ctx, fit, max_entries = n_a, chunk_size = 1L)
  expect_lt(max(abs(one$phi - mat$phi)), 1e-13)     # streamed chunks agree

  # Raw MLE score: the x = 1 score vanishes with raw probabilities, not with
  # clipped probabilities substituted.
  Zr <- cbind(1, fx$units$x)
  pr <- stats::plogis(drop(Zr %*% unname(fit$coefficients)))
  cell <- or$risk & fx$units$x == 1
  raw <- sum(or$fw[cell] * (or$D[cell] - pr[cell]))
  sub <- sum(or$fw[cell] * (or$D[cell] - pmin(pr[cell], clip_cap)))
  expect_gt(abs(sub), 1e-6)
  expect_lt(abs(raw), 1e-3 * abs(sub))

  # Mutation 1: the odds derivative of capped controls kept.
  ctrl <- or$risk & or$D == 0
  Gb_un <- -colSums((or$b * (or$dy - or$C))[ctrl] * Zr[ctrl, , drop = FALSE]) / or$sum_b
  m1 <- clip_cp02_reload(fit, Gb_un, fit$partials$G_r, fit$partials$A_beta_r)
  phi1 <- cdid:::materialize_ipw_if(ctx, m1, max_entries = n_a, chunk_size = 5L)$phi
  expect_gt(max(abs(phi1 - fd)), 5e-4)
  # Mutation 2: the derivative of the floored rate kept.
  j <- 2L; idx <- or$risk & or$k %in% j
  A_un <- fit$partials$A_beta_r
  A_un[, j] <- -colSums((or$w * (or$D - pr) / or$rc[j]^2)[idx] * Zr[idx, , drop = FALSE]) / n_a
  Gr_un <- fit$partials$G_r
  Gr_un[j] <- sum((or$w * or$b / or$fw * (or$dy - or$C))[idx]) / (or$rc[j]^2 * or$sum_b)
  m2 <- clip_cp02_reload(fit, fit$partials$G_beta, Gr_un, A_un)
  phi2 <- cdid:::materialize_ipw_if(ctx, m2, max_entries = n_a, chunk_size = 5L)$phi
  expect_gt(max(abs(phi2 - fd)), 1e-3)

  # Custom callbacks receive the floored package weights; the package-owned
  # rate derivative of the floored cohort is zero.
  seen <- new.env(parent = emptyenv())
  cb <- function(y1, y0, D, covariates, i.weights, inffunc, ...) {
    seen$w <- i.weights; seen$y0 <- y0
    m <- length(y1)
    list(ATT = 0.25, att.inf.func = rep(0, m), log_weight_gradient = seq_len(m) / m)
  }
  cfit <- cdid:::fit_custom_ipw_link(ctx, 3, c(4, Inf), 1, 2, custom_est_method = cb, link_id = "C12")
  pos <- match(seen$y0, fx$Y[, 1])
  expect_false(anyNA(pos))
  expect_lt(max(abs(seen$w / or$fw[pos] - 1)), 1e-12)
  gk <- or$k[pos]; gg <- seq_along(pos) / length(pos)
  b_exp <- vapply(1:3, function(jj) if (or$r[jj] > clip_eps) -sum(gg[gk == jj]) / or$r[jj] else 0, numeric(1))
  expect_lt(max(abs(unname(cfit$b_rate) - b_exp)), 1e-12)
  expect_identical(unname(cfit$b_rate[2]), 0)
})

test_that("CP02 link: just-active odds and rate keep their derivatives", {
  clip_skip()
  fx <- clip_cp02_fixture("active")
  run <- clip_cp02_fit(fx); ctx <- run$ctx; fit <- run$fit
  or <- clip_cp02_oracle(fx)
  expect_gt(or$r[2], clip_eps); expect_lt(or$r[2], 3e-6)
  expect_lt(or$p[2], clip_cap); expect_gt(or$p[2], 1 - 3e-6)
  expect_lt(abs(fit$estimate - or$mu), 1e-8)
  expect_identical(unlist(fit$diagnostics$clipping[c("p_capped", "r_floored")], use.names = FALSE), c(0, 0))
  n_a <- ctx$n_a_units
  mat <- cdid:::materialize_ipw_if(ctx, fit, max_entries = n_a, chunk_size = 4L)
  expect_lt(max(abs(mat$phi - clip_cp02_fd(fx, mat$id, n_a))), 1e-6)
  expect_gt(abs(fit$partials$G_r[2]), 0)            # active rate path present
})

# ------------------------------------------------------------------------------
# Selection floors: final date product (shared endpoint block) and direct pair
# ------------------------------------------------------------------------------

# 32 weighted atoms whose within-cohort date-logit law is exact. Cohort 1 at
# x = 1 has endpoint probabilities 4e-4 and 2e-3 (both above eps) with date
# product 8e-7 (below eps); its other cells and cohort 0 are active.
clip_sel_atoms <- function() {
  grid <- expand.grid(R3 = c(0, 1), R2 = c(0, 1), R1 = c(0, 1), x = c(-1, 1), cohort = c(0L, 1L))
  K <- as.integer(grid$cohort); X <- as.double(grid$x)
  R <- cbind(as.double(grid$R1), as.double(grid$R2), as.double(grid$R3))
  alpha <- rbind(c(-0.5, 0.3, 0.6), stats::qlogis(c(4e-4, 2e-3, 0.4)) + 3)
  beta <- c(0.4, -3)
  P <- stats::plogis(alpha[K + 1L, , drop = FALSE] + beta[K + 1L] * X)
  mass <- matrix(c(0.21, 0.29, 0.31, 0.19), 2L, byrow = TRUE)
  W <- length(K) * mass[cbind(K + 1L, ifelse(X > 0, 2L, 1L))] * apply(P^R * (1 - P)^(1 - R), 1L, prod)
  specs <- list()
  for (k in 1:2) for (e in list(c(1, 2), c(1, 3), c(2, 3)))
    specs[[length(specs) + 1L]] <- list(cohort = k, pair = list(s = e[1], t = e[2],
      s_calendar_index = as.integer(e[1]), t_calendar_index = as.integer(e[2])))
  list(K = K, X = X, R = R, W = W, A = as.double(rowSums(R) >= 2), N = length(K), specs = specs,
    Z = cdid:::selection_model_design(~ 1 + x, data.frame(x = X)),
    probes = list(cohort = c(0, 0, 1, 1), Z = matrix(c(-1, 1, -1, 1), ncol = 1L)))
}

# Package target order: probes, then that cohort's pair cohorts.
clip_sel_targets <- function(at, fit) {
  do.call(rbind, lapply(seq_along(at$probes$cohort), function(pi) {
    k <- match(at$probes$cohort[pi], c(0, 1))
    qs <- which(vapply(fit$pair_cohorts, function(pc) pc$cohort == k, logical(1)))
    data.frame(q = qs, k = k, x = at$probes$Z[pi, 1],
      s = vapply(fit$pair_cohorts[qs], function(pc) as.double(pc$pair$s), numeric(1)),
      t = vapply(fit$pair_cohorts[qs], function(pc) as.double(pc$pair$t), numeric(1)))
  }))
}

clip_sel_oracle <- function(at, tg, model, cw = rep(1, at$N)) {
  w <- at$W * cw
  a <- vapply(0:1, function(k) sum(w[at$K == k & at$A > 0]) / sum(w[at$K == k]), numeric(1))
  if (identical(model, "date")) {
    n <- at$N; rows <- rep(seq_len(n), 3L); tt <- rep(1:3, each = n)
    Xd <- matrix(0, 3L * n, 8L)
    Xd[cbind(seq_len(3L * n), at$K[rows] * 3L + tt)] <- 1
    Xd[cbind(seq_len(3L * n), 7L + at$K[rows])] <- at$X[rows]
    co <- clip_quasi(Xd, as.vector(at$R), w[rows])
    pd <- function(k, s, x) stats::plogis(co[(k - 1) * 3 + s] + co[6 + k] * x)
    pp <- pd(tg$k, tg$s, tg$x) * pd(tg$k, tg$t, tg$x)
  } else {
    pp <- vapply(seq_len(nrow(tg)), function(j) {
      cell <- at$K == tg$k[j] - 1L & at$X == tg$x[j]
      S <- at$R[, tg$s[j]] * at$R[, tg$t[j]]
      sum(w[cell] * S[cell]) / sum(w[cell])
    }, numeric(1))
  }
  list(pp = pp, a = a[tg$k], v = a[tg$k] / pmax(pp, clip_eps))
}

clip_sel_fd <- function(at, tg, model, what, h = 1e-5) {
  t(vapply(seq_len(at$N), function(i) {
    cp <- cm <- rep(1, at$N); cp[i] <- 1 + h; cm[i] <- 1 - h
    at$N * (clip_sel_oracle(at, tg, model, cp)[[what]] - clip_sel_oracle(at, tg, model, cm)[[what]]) / (2 * h)
  }, numeric(nrow(tg))))
}

clip_rel <- function(x, y) max(abs(x - y)) / max(1, max(abs(y)))

test_that("date-product and direct-pair floors: multipliers, raw probability and multiplier IFs", {
  clip_skip()
  at <- clip_sel_atoms()
  gc <- list(epsilon = 1e-15, maxit = 60L, trace = FALSE)
  S_list <- list("1:2" = at$R[, 1] * at$R[, 2], "1:3" = at$R[, 1] * at$R[, 3], "2:3" = at$R[, 2] * at$R[, 3])
  for (model in c("date", "pair")) {
    fit <- if (identical(model, "date")) {
      cdid:::evaluate_selection_table_if(Z = at$Z, cohort_index = at$K, w = at$W, A = at$A,
        cohorts = c(0L, 1L), all_times = c(1, 2, 3), date_responses = at$R, pair_specs = at$specs,
        selection_model = "independent_period_logit", selection_pooling = "within_cohort",
        glm_control = gc, probes = at$probes, chunk_rows = 7L, details = TRUE)
    } else {
      cdid:::evaluate_selection_table_if(Z = at$Z, cohort_index = at$K, w = at$W, A = at$A,
        cohorts = c(0L, 1L), all_times = c(1, 2, 3), pair_responses = S_list, pair_specs = at$specs,
        selection_model = "pair_logit", selection_pooling = "separate",
        glm_control = gc, probes = at$probes, chunk_rows = 7L, details = TRUE)
    }
    tg <- clip_sel_targets(at, fit)
    or <- clip_sel_oracle(at, tg, model)
    fl <- which(tg$k == 2L & tg$x == 1 & tg$s == 1 & tg$t == 2)
    expect_length(fl, 1L)
    expect_lt(or$pp[fl], clip_eps)
    expect_true(all(or$pp[-fl] > clip_eps))
    expect_lt(clip_rel(unname(fit$pair_probabilities), or$pp), 1e-8)   # raw products
    expect_lt(clip_rel(unname(fit$multipliers), or$v), 1e-8)
    expect_identical(unname(fit$multipliers[fl]), unname(fit$a[2]) / clip_eps)

    fd_v <- clip_sel_fd(at, tg, model, "v")
    fd_p <- clip_sel_fd(at, tg, model, "pp")
    expect_lt(clip_rel(unname(fit$multiplier_if), fd_v), 1e-6)
    expect_lt(clip_rel(unname(fit$pair_prob_if), fd_p), 1e-6)
    # The floored multiplier keeps only its repeat-rate path.
    expect_lt(clip_rel(fit$multiplier_if[, fl], fit$a_if[, 2] / clip_eps), 1e-12)
    # Mutation: the unmasked floored-multiplier derivative.
    a2 <- unname(fit$a[2])
    mut <- fit$a_if[, 2] / clip_eps - (a2 / clip_eps^2) * fit$pair_prob_if[, fl]
    expect_gt(clip_rel(mut, fd_v[, fl]), 1e-2)
  }
})

test_that("selection forcing carries the final pair mask on both endpoints of a shared block", {
  clip_skip()
  at <- clip_sel_atoms()
  gc <- list(epsilon = 1e-15, maxit = 60L, trace = FALSE)
  fit <- cdid:::evaluate_selection_table_if(Z = at$Z, cohort_index = at$K, w = at$W, A = at$A,
    cohorts = c(0L, 1L), all_times = c(1, 2, 3), date_responses = at$R, pair_specs = at$specs,
    selection_model = "independent_period_logit", selection_pooling = "within_cohort",
    glm_control = gc, probes = at$probes, chunk_rows = 7L, details = TRUE)
  q <- which(vapply(fit$pair_cohorts, function(pc) pc$cohort == 2L && pc$pair$s == 1 && pc$pair$t == 2,
    logical(1)))
  expect_length(q, 1L)
  pc <- fit$pair_cohorts[[q]]
  expect_identical(pc$s_block_id, pc$t_block_id)
  rows <- which(at$K == 1L & at$R[, 1] == 1 & at$R[, 2] == 1)
  ow <- rep(1, length(rows))
  acc <- cdid:::accumulate_selection_table_moments(fit, at$Z, q, rows, ow, sign = -1)
  n <- at$N; rr <- rep(seq_len(n), 3L); tt <- rep(1:3, each = n)
  Xd <- matrix(0, 3L * n, 8L)
  Xd[cbind(seq_len(3L * n), at$K[rr] * 3L + tt)] <- 1
  Xd[cbind(seq_len(3L * n), 7L + at$K[rr])] <- at$X[rr]
  co <- clip_quasi(Xd, as.vector(at$R), at$W[rr])
  ps <- stats::plogis(co[4] + co[8] * at$X[rows]); pt <- stats::plogis(co[5] + co[8] * at$X[rows])
  act <- as.double(ps * pt > clip_eps)
  expect_identical(act, as.double(at$X[rows] < 0))
  bl <- acc$blocks[[pc$s_block_id]]
  expect_lt(abs(bl$t_alpha[pc$s_local_b] + sum(ow * (1 - ps) * act) / n), 1e-9)
  expect_lt(abs(bl$t_alpha[pc$t_local_b] + sum(ow * (1 - pt) * act) / n), 1e-9)
  tb <- -sum(ow * ((1 - ps) + (1 - pt)) * act * at$X[rows]) / n
  expect_lt(abs(bl$t_beta - tb), 1e-9)
  tb_unmasked <- -sum(ow * ((1 - ps) + (1 - pt)) * at$X[rows]) / n
  expect_gt(abs(bl$t_beta - tb_unmasked), 1e-2)
  expect_lt(abs(acc$t_a[2] - sum(ow / fit$a[2]) / n), 1e-12)   # a_U path unmasked
})

# ------------------------------------------------------------------------------
# Selection-corrected IPW: direct floor and treatment cap in one link
# ------------------------------------------------------------------------------

# Selection uses z, treatment uses x. Cohort Inf at z = 0 has pair
# probability 6e-7 (floored; its observers have small weights, so they stay
# comparable inside their treatment cell). Treatment cell x = 1 has
# p ~ 1 - 2.2e-8 (capped) and x = 0 has p ~ 0.46 (active). Rows off A
# (singletons, a zero-observation row) are kept.
clip_cor_fixture <- function() {
  u <- data.frame(id = 1:18,
    g = c(rep(3, 8), rep(Inf, 10)),
    x = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0),
    z = c(0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1),
    w = c(rep(1, 8), 1e-6, 2e-6, 4, 1e-8, 2e-8, 1e-8, 1e-8, 1, 1.5, 1),
    pattern = c("123", "12", "13", "1", "12", "123", "23", "", "12", "12", "13", "12", "123", "1", "23", "2",
      "12", "13"),
    stringsAsFactors = FALSE)
  obs <- unname(t(vapply(u$pattern, function(p) vapply(1:3, function(s) grepl(as.character(s), p), logical(1)),
    logical(3))))
  Y <- outer(u$id, 1:3, function(i, s) round(0.1 * i + 0.37 * sin(2.3 * i + s) + 0.5 * s * (u$g[i] == 3), 6))
  list(units = u, obs = obs, Y = Y)
}

clip_cor_oracle <- function(fx, cw = rep(1, nrow(fx$units))) {
  u <- fx$units; w <- u$w * cw
  A <- rowSums(fx$obs) >= 2; S <- fx$obs[, 1] & fx$obs[, 2]
  kk <- match(u$g, c(3, Inf))
  a <- vapply(1:2, function(k) sum(w[kk == k & A]) / sum(w[kk == k]), numeric(1))
  q <- numeric(nrow(u))
  for (k in 1:2) for (v in 0:1) {
    cell <- kk == k & u$z == v
    q[cell] <- sum(w[cell & S]) / sum(w[cell])
  }
  fw <- ifelse(S, w * a[kk] / pmax(q, clip_eps), 0)
  D <- as.double(u$g == 3)
  dy <- ifelse(S, fx$Y[, 2] - fx$Y[, 1], 0)
  p <- cell_odds <- numeric(2)
  for (v in 0:1) {
    mt <- sum(fw[S & u$x == v & D == 1]); mc <- sum(fw[S & u$x == v & D == 0])
    p[v + 1] <- mt / (mt + mc)
    cell_odds[v + 1] <- if (p[v + 1] < clip_cap) mt / mc else clip_cap_odds
  }
  wm <- fw * (1 - D) * cell_odds[u$x + 1]; wp <- fw * D
  M <- sum(wm * dy) / sum(wm)
  list(mu = sum(wp * dy) / sum(wp) - M, q = q, p = p, fw = fw, D = D, S = S, dy = dy, wm = wm, M = M)
}

test_that("selection-corrected IPW: floored pair weight and capped odds match the oracle IF; masks matter", {
  clip_skip()
  fx <- clip_cor_fixture(); u <- fx$units; N <- nrow(u)
  res <- cdid:::fit_cp04_selection_table(Z = cdid:::selection_model_design(~ 1 + z, data.frame(z = u$z)),
    V = cbind(1, u$x), Y = fx$Y, obs = fx$obs * 1, g = u$g, w = u$w, calendar = c(1, 2, 3),
    links = list(list(g = 3, K = Inf, s = 1, t = 2, link_id = "L12")),
    selection_model = "pair_logit", selection_pooling = "separate")
  or <- clip_cor_oracle(fx)
  expect_lt(max(or$q[u$g == Inf & u$z == 0]), clip_eps)
  expect_true(all(or$q[or$S & !(u$g == Inf & u$z == 0)] > clip_eps))
  expect_gt(or$p[2], clip_cap); expect_lt(or$p[1], clip_cap)
  tr <- res$treatments[[1]]
  expect_lt(abs(unname(res$estimates[1]) - or$mu), 1e-8)

  h <- 1e-6
  fd <- vapply(seq_len(N), function(i) {
    cp <- cm <- rep(1, N); cp[i] <- 1 + h; cm[i] <- 1 - h
    N * (clip_cor_oracle(fx, cp)$mu - clip_cor_oracle(fx, cm)$mu) / (2 * h)
  }, numeric(1))
  expect_lt(max(abs(res$phi[, 1] - fd)), 1e-6)
  expect_gt(abs(fd[16]), 1e-4)       # off-A singleton moves a_U while its cell floor binds

  # Raw treatment score in the capped cell vanishes only with raw probabilities.
  X <- cbind(1, u$x); pr <- stats::plogis(drop(X %*% unname(tr$coefficients)))
  cell <- or$S & u$x == 1
  raw <- sum(or$fw[cell] * (or$D[cell] - pr[cell]))
  sub <- sum(or$fw[cell] * (or$D[cell] - pmin(pr[cell], clip_cap)))
  expect_gt(abs(sub), 1e-7)
  expect_lt(abs(raw), 1e-3 * abs(sub))

  # Mutation: the odds derivative of the capped controls kept in V.
  capc <- or$S & or$D == 0 & u$x == 1
  dV2 <- colSums((or$wm * (or$dy - or$M))[capc] * X[capc, , drop = FALSE]) / sum(or$wm)
  mut <- res$phi[, 1] - drop(tr$beta_if %*% dV2)
  expect_gt(max(abs(mut - fd)), 1e-3)
})

# ------------------------------------------------------------------------------
# Conditional-pretest loading integrand
# ------------------------------------------------------------------------------

test_that("conditional two-arm loading uses the active-odds mask on propensity coordinates", {
  clip_skip()
  X <- matrix(c(1, 1, 0.5, -1), 2L, 2L)
  lam_up <- matrix(c(0.1, 0.2), 2L, 1L)
  pc <- cdid:::cpt_two_arm_piece(rows = 1:3, rr = c(1L, 3L), D = c(1, 0), fw = c(2, 4), X = X,
    odds = c(1.5, clip_cap_odds), active = c(1, 0), rho = c(0.3, -0.7), lam_up = lam_up,
    xi = matrix(0, 3L, 3L), D_T = 0.7, D_M = 1.3, n_up = 1L, p_trt = 2L)
  nu2 <- c(0, 0, 4 * clip_cap_odds)
  lamM <- matrix(0, 3L, 3L); lamM[c(1, 3), 1] <- lam_up; lamM[1, 2:3] <- X[1, ]
  lamT <- matrix(0, 3L, 3L); lamT[c(1, 3), 1] <- lam_up
  rho <- c(0.3, 0, -0.7); nu1 <- c(2, 0, 0)
  P <- (nu1 * rho / 0.7) * lamT - (nu2 * rho / 1.3) * lamM
  expect_lt(max(abs(pc$P - P)), 1e-12 * max(1, abs(P)))
  expect_identical(pc$blam[, 2], colSums(nu2 * lamM))
  expect_identical(unname(pc$blam[2:3, 2]), c(0, 0))   # capped row: no propensity loading
})
