# CP08 O14/O15 (D63/D64/D66): the chained conditional pretrend test and the
# optional custom-callback log-weight gradient. Deterministic checks against the
# accepted independent references (saved values; nothing is refit or rerun from
# them), plus public-path RNG, chunk, worker, unavailability and resource checks.
# Finite fixtures establish formulas and wiring, not size or power.

cpt_path <- function(...) test_path("fixtures", "cpt", ...)
cpt_cache <- new.env(parent = emptyenv())
cpt_rds <- function(f) {
  if (!exists(f, envir = cpt_cache, inherits = FALSE)) assign(f, readRDS(cpt_path(f)), envir = cpt_cache)
  get(f, envir = cpt_cache, inherits = FALSE)
}
cpt_csv <- function(f) utils::read.csv(cpt_path(f))
cpt_signs <- function(x) do.call(rbind, lapply(x, function(r) as.double(unlist(r))))
cpt_key <- function(m) apply(m, 1L, function(r) paste(sprintf("%a", r + 0), collapse = "|"))

cpt_cases <- list(
  A1 = list(fx = "A", csv = "fA.csv", rds = "full.rds", args = list(), routes = c(ipw = "ipw", reg = "reg", dr = "dr")),
  B4 = list(fx = "B", csv = "fB.csv", rds = "full.rds", args = list(selection = "conditional_mtar",
    selection_model = "pair_logit", selection_pooling = "across_cohorts", selection_xformla = ~sx + sz),
    routes = c(ipw = "ipw", reg = "reg", dr = "dr")),
  B8 = list(fx = "B", csv = "fBd.csv", rds = "full.rds", args = list(selection = "conditional_mtar",
    selection_model = "independent_period_logit", selection_pooling = "within_cohort", selection_xformla = ~sx + sz,
    fix_weights = "base_period"), routes = c(ipw = "ipw", reg = "reg", dr = "dr")),
  B2 = list(fx = "B", csv = "fB.csv", rds = "imp.rds", args = list(selection = "conditional_mtar",
    selection_est_method = "imputation", imputation_pooling = "within_cohort"),
    routes = c(ipw = "imp_ipw", reg = "imp_reg", dr = "imp_dr")),
  B3 = list(fx = "B", csv = "fB.csv", rds = "imp.rds", args = list(selection = "conditional_mtar",
    selection_est_method = "dr", imputation_pooling = "within_cohort", selection_model = "independent_period_logit",
    selection_pooling = "within_cohort", selection_xformla = ~sx + sz), routes = c(dr = "aug_dr")))

cpt_defaults <- function(data) list(yname = "y", tname = "time", idname = "id", gname = "G", xformla = ~x + z,
  data = data, panel = TRUE, allow_unbalanced_panel = NULL, control_group = "notyettreated", weightsname = "w",
  alp = 0.05, bstrap = TRUE, cband = TRUE, biters = 100, clustervars = NULL, est_method = "ipw", print_details = FALSE,
  pl = FALSE, cores = 1, anticipation = 0, base_period = "varying", fix_weights = NULL, mode = NULL, frame = NULL,
  selection = "none", selection_model = "pair_logit", selection_pooling = "separate", selection_xformla = NULL,
  selection_universe = "F", selection_rho = NULL, selection_justification = NULL, selection_est_method = "ipw",
  imputation_pooling = "separate", control = NULL)

# Private executor: caller-fixed signs (private CvM arithmetic only) and kept rows.
cpt_fixed <- function(data, overrides, signs) {
  a <- cpt_defaults(data)
  for (nm in names(overrides)) a[nm] <- list(overrides[[nm]])
  a$biters <- nrow(signs); a$fixed_signs <- signs; a$keep_rows <- TRUE
  do.call(chained_conditional_did_pretest, a)
}

cpt_case_data <- function(case) {
  cs <- cpt_cases[[case]]
  d <- cpt_csv(cs$csv)
  fx <- cpt_rds("priv.rds")$fixtures$fixtures[[cs$fx]]
  d$cl <- as.integer(unlist(fx$clusters))[match(d$id, sort(unique(d$id)))]
  list(data = d, fx = fx, cs = cs)
}

cpt_cell_key <- function(id, bp) {
  gt <- regmatches(id, regexec("g=([0-9.]+),t=([0-9.]+)", id))[[1L]]
  list(arr = paste0(bp, "__pre_", if (bp == "varying") "var" else "uni", "_g", gt[2L], "_t", gt[3L], "__"),
    cell = paste0(bp, "|pre_", if (bp == "varying") "var" else "uni", "_g", gt[2L], "_t", gt[3L]))
}

# Max |deviation| of the reference mutations (R3) from the correct rows.
cpt_mutation_devs <- function(cell) {
  r <- cell$rows
  H <- cpt_indicator(cell$zd[r$map, , drop = FALSE], cell$zd)
  phi <- cpt_cell_process_rows(cell); s1 <- cpt_cell_scalar_rows(cell)
  term <- function(sel) if (any(sel)) max(abs(cell$sigma * r$xi[, sel, drop = FALSE] %*% cell$Lmat[sel, , drop = FALSE])) else 0
  c(h_times_scalar_IF = max(abs(phi - H * s1)),
    scalar_loadings = max(abs(cell$sigma * r$xi %*% (cell$Lmat - cell$L1))),
    dropped_xi_alpha = term(grepl("^or_alpha:", cell$coords)),
    dropped_xi_alpha_OR = term(grepl("^oa:", cell$coords)),
    dropped_rate_selection_loading = term(grepl("^(up|rate):", cell$coords)))
}

for (case in names(cpt_cases)) local({
  case <- case
  test_that(paste0("conditional process, draws and CvM reproduce the accepted ", case, " references"), {
    cd <- cpt_case_data(case); cs <- cd$cs; ref <- cpt_rds(cs$rds)
    base <- ref$reference$cases[[case]]$base; arr <- ref$arrays[[case]]
    S_unit <- cpt_signs(cd$fx$unit_rademacher_signs); S_cl <- cpt_signs(cd$fx$rademacher_signs)
    u_cl <- unique(cd$data$cl[!duplicated(cd$data$id)])
    for (est in names(cs$routes)) for (bp in c("varying", "universal")) {
      rn <- cs$routes[[est]]; info <- paste(case, rn, bp)
      ov <- c(cs$args, list(est_method = est, base_period = bp))
      r <- cpt_fixed(cd$data, ov, S_unit)
      expect_s3_class(r, c("cdid_MP.TEST", "MP.TEST"))
      expect_identical(r$chained$status, "available", info = info)
      expect_null(r$KS); expect_null(r$KSb); expect_null(r$KSpval)
      gold <- base$cvm[[paste0(bp, "|", rn)]]
      expect_equal(r$CvM, gold$unit$CvM, tolerance = 1e-10, info = info)
      expect_lt(max(abs(r$CvMb - unlist(gold$unit$CvMb))), 1e-9)
      expect_identical(r$CvMpval, gold$unit$CvMpval, info = info)
      expect_equal(r$CvMcval, gold$unit$CvMcval, tolerance = 1e-10, info = info)
      cells <- attr(r, "cells")
      expect_identical(length(cells), 3L)
      for (cell in cells) {
        k <- cpt_cell_key(cell$id, bp)
        rk <- base$cells[[k$cell]]
        gpts <- cpt_signs(base$grids[[rk$recipe]]$distinct_points)
        idx <- match(cpt_key(gpts), cpt_key(cell$zd))
        expect_false(anyNA(idx), info = info)
        expect_equal(sum(cell$mult), sum(!duplicated(cd$data$id)))
        expect_lt(max(abs(cpt_cell_process_rows(cell)[, idx] - arr[[paste0(k$arr, rn, "__phi_distinct")]])), 1e-10)
        expect_lt(max(abs(cell$J[cell$rows$map] - arr[[paste0(k$arr, rn, "__J_full")]])), 1e-12)
        expect_lt(max(abs(cpt_cell_scalar_rows(cell) - arr[[paste0(k$arr, rn, "__phi_att")]])), 1e-10)
        # R3 mutations (h * scalar IF, scalar loadings, dropped coordinates) are
        # distinguishable from the correct process exactly as in the reference.
        mut <- unlist(rk$routes[[rn]]$mutations_max_abs_dev_from_correct_IF)
        mine <- cpt_mutation_devs(cell)
        for (tag in names(mut)) expect_lt(abs(unname(mine[tag]) - unname(mut[tag])), 1e-8 * max(1, abs(unname(mut[tag]))))
      }
      if (identical(case, "A1") || identical(est, names(cs$routes)[1L])) {
        rc <- cpt_fixed(cd$data, c(ov, list(clustervars = c("id", "cl"))), S_cl[, u_cl + 1L, drop = FALSE])
        expect_identical(rc$chained$inference$n_clusters, length(u_cl))
        expect_equal(rc$CvM, gold$cluster$CvM, tolerance = 1e-10, info = info)
        expect_lt(max(abs(rc$CvMb - unlist(gold$cluster$CvMb))), 1e-9)
        expect_identical(rc$CvMpval, gold$cluster$CvMpval, info = info)
        expect_equal(rc$CvMcval, gold$cluster$CvMcval, tolerance = 1e-10, info = info)
      }
    }
  })
})

test_that("zero-IF clusters keep their sign columns and grid multiplicities are retained", {
  cd <- cpt_case_data("A1")
  zc <- cd$fx$zero_if_cluster_id
  u_cl <- unique(cd$data$cl[!duplicated(cd$data$id)])
  S_cl <- cpt_signs(cd$fx$rademacher_signs)[, u_cl + 1L, drop = FALSE]
  rc <- cpt_fixed(cd$data, list(clustervars = c("id", "cl")), S_cl)
  units <- sort(unique(cd$data$id))
  zero_units <- which(as.integer(unlist(cd$fx$clusters)) == zc)
  for (cell in attr(rc, "cells")) {
    expect_identical(max(abs(cpt_cell_process_rows(cell)[zero_units, , drop = FALSE])), 0)
    expect_identical(nrow(cell$AggX), length(u_cl))
  }
  expect_identical(rc$chained$inference$n_clusters, 12L)
  expect_identical(rc$chained$status, "available")
})

# Backend complete IF row (frame scale, unoriented) of auxiliary descriptor q.
cpt_backend_rows <- function(state, q) {
  N <- state$n_total
  if (is.null(state$cp04_reader)) {
    pl <- state$auxiliary$point_links[[q]]
    A <- numeric(state$n_a_units)
    replay_ipw_point_if_view(list(pl), state$control$point_control, callback_schema = TRUE,
      callback = function(b) A[b$a_index] <<- b$phi[, 1L])
    out <- numeric(N); out[state$a_to_full] <- (N / state$n_a_units) * A
    return(out)
  }
  reg <- state$registry; j <- reg$cp04_aux_index[q]; out <- numeric(N)
  B <- reg$cp04$control$max_chunk_rows
  for (s0 in seq.int(1L, N, by = B)) {
    rr <- seq.int(s0, min(N, s0 + B - 1L))
    out[rr] <- chained_full_if_chunk(reg$cp04, rr, j, "frame")$phi[, 1L]
  }
  out
}

cpt_q_of <- function(state, id) {
  for (d in state$aux_plan$selections[[state$base_period]])
    if (identical(paste0("pre(g=", fmt_design_num(d$g), ",t=", fmt_design_num(d$t), ")"), id)) return(d$index)
  NA_integer_
}

test_that("on routes without process references the h = 1 process equals the backend's own complete IF", {
  B <- cpt_csv("fB.csv"); Bd <- cpt_csv("fBd.csv")
  sel <- list(selection = "conditional_mtar", selection_model = "pair_logit", selection_pooling = "separate",
    selection_xformla = ~sx + sz)
  runs <- list(
    legacy_pair_logit_ipw = list(B, c(sel, list(est_method = "ipw")), "cdid_cp04_corollary1"),
    legacy_pair_logit_or = list(B, c(sel, list(est_method = "reg")), "cdid_or_corollary1"),
    dated_baseline_ipw = list(Bd, list(est_method = "ipw", fix_weights = "base_period"), "cdid_reference_ipw_backend"),
    dated_baseline_or = list(Bd, list(est_method = "reg", fix_weights = "base_period", base_period = "universal"), "cdid_or_baseline"),
    dated_baseline_dr = list(Bd, list(est_method = "dr", fix_weights = "base_period"), "cdid_dr_baseline"),
    dated_imputation = list(Bd, list(est_method = "dr", fix_weights = "base_period", selection = "conditional_mtar",
      selection_est_method = "imputation", imputation_pooling = "within_cohort"), "cdid_imputation_backend"),
    dated_augmented = list(Bd, list(est_method = "dr", fix_weights = "base_period", selection = "conditional_mtar",
      selection_est_method = "dr", selection_model = "independent_period_logit", selection_pooling = "within_cohort",
      selection_xformla = ~sx + sz), "cdid_augmented_backend"))
  signs <- matrix(rep(c(1, -1), length.out = 19L * 192L), 19L, 192L)
  for (nm in names(runs)) {
    rr <- runs[[nm]]
    r <- cpt_fixed(rr[[1L]], rr[[2L]], signs)
    expect_identical(r$chained$status, "available", info = nm)
    st <- attr(r, "state")
    expect_true(inherits(st$registry$cp04, rr[[3L]]), info = nm)
    for (cell in attr(r, "cells")) {
      q <- cpt_q_of(st, cell$id)
      expect_lt(max(abs(cpt_cell_scalar_rows(cell) - cell$sigma * cpt_backend_rows(st, q))), 1e-9)
      expect_equal(cell$sigma * cell$mu1, cell$estimate_scalar, tolerance = 1e-10, info = paste(nm, cell$id))
      expect_lt(max(abs(cell$denominator_check$supplied - cell$denominator_check$accumulated)), 1e-12)
    }
    if (nm %in% c("dated_baseline_ipw", "dated_imputation")) expect_identical(r$chained$grid$columns[[1L]], c("(Intercept)", "x", "z"))
  }
})

cpt_synthetic <- function(delta) {
  n <- 40L; id <- seq_len(n); G <- rep(c(3, 0), each = 20L); x <- rep(rep(c(-1, 1), each = 10L), 2L)
  e <- rep(c(-0.5, 0.5), 20L)
  d <- data.frame(id = rep(id, each = 3L), time = rep(1:3, n), G = rep(G, each = 3L), x = rep(x, each = 3L))
  y2 <- 1 + delta * x * (G == 3) + e
  d$y <- 0
  d$y[d$time == 2L] <- y2
  d$y[d$time == 3L] <- y2 + 0.5 * (G == 3)
  d
}

test_that("a conditional alternative with zero unconditional mean is detected by the process, not the scalar", {
  for (delta in c(0, 1)) {
    d <- cpt_synthetic(delta)
    signs <- matrix(rep(c(1, -1, -1, 1), length.out = 9L * 40L), 9L, 40L)
    a <- list(xformla = ~x, weightsname = NULL, control_group = "nevertreated")
    r <- cpt_fixed(d, a, signs)
    cell <- attr(r, "cells")[[1L]]
    expect_lt(abs(r$chained$cells$estimate), 1e-10)
    lo <- which(cell$zd[, 2L] == -1)
    expect_lt(max(abs(cell$J[lo] + delta / 2)), 1e-10)
    expect_lt(abs(r$CvM - 5 * delta^2), 1e-10)
    expect_identical(r$chained$status, "available")
  }
})

# Fitted union registry and auxiliary links of a request (pipeline state).
cpt_state <- function(data, overrides) {
  a <- cpt_defaults(data); for (nm in names(overrides)) a[nm] <- list(overrides[[nm]])
  a$biters <- 9; a$fixed_signs <- matrix(1, 9L, length(unique(data$id))); a$keep_rows <- TRUE
  attr(do.call(chained_conditional_did_pretest, a), "state")
}

test_that("missing grid covariates or unsupported levels make the complete test unavailable; estimates remain", {
  d <- cpt_csv("fA.csv")
  d <- d[d$G != 1, ]
  d$y[d$id == 5 & d$time > 1] <- NA
  d$x[d$id == 5] <- NA
  r <- conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G", xformla = ~x + z, data = d,
    control_group = "notyettreated", weightsname = "w", est_method = "reg", biters = 9)
  expect_s3_class(r, "cdid_MP.TEST")
  expect_identical(r$chained$status, "unavailable")
  expect_true(is.na(r$CvM)); expect_null(r$CvMb)
  expect_match(r$chained$reason, "id 5", fixed = TRUE)
  expect_identical(unique(r$chained$grid$missing$id), 5L)
  # The point estimates of the same request remain available (A-row designs).
  st <- cpt_state(d, list(est_method = "reg"))
  expect_true(all(is.finite(vapply(st$registry$point_links, `[[`, numeric(1L), "estimate"))))
  d2 <- cpt_csv("fA.csv"); d2 <- d2[d2$G != 1, ]
  d2$grp <- ifelse(d2$id %% 2L == 0L, "a", "b")
  d2$y[d2$id == 5 & d2$time > 1] <- NA
  d2$grp[d2$id == 5] <- "c"
  r2 <- conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G", xformla = ~x + z + grp,
    data = d2, control_group = "notyettreated", weightsname = "w", est_method = "dr", biters = 9)
  expect_identical(r2$chained$status, "unavailable")
  expect_identical(r2$chained$grid$unsupported$value, "c")
  expect_match(r2$chained$reason, "unsupported levels", fixed = TRUE)
})

test_that("caller RNG, grid/draw chunks and declared workers give the same draws", {
  d <- cpt_csv("fA.csv")
  call_it <- function(...) conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G",
    xformla = ~x + z, data = d, control_group = "notyettreated", weightsname = "w", biters = 120, ...)
  set.seed(11); r1 <- call_it(); u1 <- stats::runif(1L)
  set.seed(11); r2 <- call_it()
  expect_identical(r1$CvMb, r2$CvMb)
  expect_identical(r1$chained$inference$rng$mode, "caller_rng_stream_Rademacher")
  set.seed(12); r3 <- call_it()
  expect_false(identical(r1$CvMb, r3$CvMb))
  # One common sign stream: batches in draw order, column-major in a batch.
  N <- 48L; br <- r1$chained$inference$rng$batch_rows
  set.seed(11)
  S <- do.call(rbind, lapply(split(seq_len(120L), ceiling(seq_len(120L) / br)), function(ii)
    matrix(sample(c(-1.0, 1.0), length(ii) * N, replace = TRUE), length(ii), N)))
  expect_identical(stats::runif(1L), u1)
  a <- cpt_defaults(d)
  r4 <- cpt_fixed(d, list(est_method = "ipw"), S)
  expect_equal(r4$CvMb, r1$CvMb, tolerance = 1e-13)
  expect_identical(r4$CvM, r1$CvM)
  # Smaller extraction/grid chunks change no statistic or draw for given signs.
  r5 <- cpt_fixed(d, list(control = list(chunk_rows = 5)), S)
  expect_equal(r5$CvM, r4$CvM, tolerance = 1e-12)
  expect_equal(r5$CvMb, r4$CvMb, tolerance = 1e-12)
  # Declared multiplier workers draw the same signs in the calling process.
  set.seed(11); r6 <- call_it(pl = TRUE, cores = 2, control = list(workers = 2))
  expect_identical(r6$CvMb, r1$CvMb)
  expect_identical(r6$chained$execution$workers, 2)
  expect_error(call_it(pl = TRUE, cores = 2), class = "cdid_chained_request_parallel_error")
})

test_that("the conditional test charges its buffers against the existing caps", {
  dims <- list(list(Md = 50, p = 3, R = 2, Dl = 6, chunk = 50))
  ctrl <- list(owner_control = list(max_retained_entries = 1e6, max_matrix_entries = 1e6, max_work_ops = 1e12))
  mb <- list(max_draws = 100, max_workspace_cells = 1e6, max_workspace_bytes = 8e6)
  ok <- cpt_preflight(dims, 50, 50, 100, 50, 50, 50, ctrl, mb, 0, 0, 1)
  expect_true(ok$retained_entries > 0 && ok$work_ops > 0)
  expect_error(cpt_preflight(dims, 50, 50, 101, 50, 50, 50, ctrl, mb, 0, 0, 1), class = "cdid_conditional_pretest_bound_error")
  small <- ctrl; small$owner_control$max_retained_entries <- 1000
  expect_error(cpt_preflight(dims, 50, 50, 100, 50, 50, 50, small, mb, 0, 0, 1), class = "cdid_ipw_owner_bound_error")
  slow <- ctrl; slow$owner_control$max_work_ops <- 1000
  expect_error(cpt_preflight(dims, 50, 50, 100, 50, 50, 50, slow, mb, 0, 0, 1), class = "cdid_ipw_owner_bound_error")
})

test_that("conditional_did_pretest keeps did's formals, appends the chained ones and resolves the mode lazily", {
  f <- formals(conditional_did_pretest)
  expect_identical(names(f)[seq_along(upstream_conditional_did_pretest_formals)], upstream_conditional_did_pretest_formals)
  expect_identical(names(f)[-seq_along(upstream_conditional_did_pretest_formals)], cdid_conditional_pretest_chained_formals)
  expect_identical(f$est_method, "ipw"); expect_null(f$allow_unbalanced_panel); expect_identical(f$base_period, "varying")
  # Custom estimators have no claimed conditional process: not applicable, no data forced.
  r <- conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G",
    data = stop("data must not be evaluated"), est_method = function(...) NULL)
  expect_identical(r$chained$status, "not_applicable")
  expect_error(conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G",
    data = stop("data must not be evaluated"), est_method = "stdipw"), class = "cdid_chained_unavailable")
  expect_error(conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G",
    data = stop("data must not be evaluated"), allow_unbalanced_panel = FALSE), class = "cdid_chained_balanced")
})

test_that("ordinary mode delegates to did::conditional_did_pretest without evaluating chained arguments", {
  skip_if_not_installed("did")
  d <- cpt_csv("fA.csv"); d <- d[d$G != 1, ]
  run <- function(f, ...) {
    set.seed(3)
    out <- NULL
    utils::capture.output(out <- suppressMessages(f(yname = "y", tname = "time", idname = "id", gname = "G",
      xformla = ~x, data = d, biters = 29, ...)))
    out
  }
  ours <- run(conditional_did_pretest, mode = "ordinary", selection = stop("chained argument evaluated"),
    control = stop("chained argument evaluated"), base_period = stop("chained argument evaluated"))
  ref <- run(did::conditional_did_pretest)
  expect_identical(class(ours), class(ref))
  expect_equal(ours$CvM, ref$CvM); expect_equal(ours$CvMb, ref$CvMb); expect_equal(ours$CvMpval, ref$CvMpval)
})

# ------------------------------------------------------------------------------
# O15: optional explicit total log-weight gradient of custom callbacks
# ------------------------------------------------------------------------------
cpt_wlogit <- function(X, y, w) {
  th <- rep(0, ncol(X))
  for (it in 1:200) {
    p <- stats::plogis(drop(X %*% th))
    st <- drop(solve(crossprod(X, (w * p * (1 - p)) * X), crossprod(X, w * (y - p))))
    th <- th + st
    if (max(abs(st)) < 1e-14) break
  }
  th
}
# White-box weighted Hajek/logit callbacks of the reference: the propensity is
# weighted by i.weights^ps_power (1: std_ipw_did_panel type, H1 holds; 0 and
# 1/2: H1-violating witnesses). psi is the m-normalized fixed-weight IF, g the
# exact total derivative d ATT / d log i.weights (internal refit included).
cpt_callback <- function(ps_power, give_g = TRUE) function(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...) {
  X <- covariates; dY <- y1 - y0; iw <- i.weights; m <- length(D); psw <- iw^ps_power
  gam <- cpt_wlogit(X, D, psw); eta <- drop(X %*% gam); p <- stats::plogis(eta); o <- exp(eta)
  wt <- iw * D; wc <- iw * (1 - D) * o; S1 <- sum(wt); S0 <- sum(wc)
  Rp <- sum(wt * dY) / S1; Rm <- sum(wc * dY) / S0
  direct <- wt * (dY - Rp) / S1 - wc * (dY - Rm) / S0
  lever <- drop(solve(crossprod(X, (psw * p * (1 - p)) * X), crossprod(X, wc * (dY - Rm)) / S0))
  sdir <- psw * (D - p) * drop(X %*% lever)
  out <- list(ATT = Rp - Rm, att.inf.func = m * (direct - sdir))
  if (give_g) out$log_weight_gradient <- direct - ps_power * sdir
  out
}

# Fitted pipeline (no did result binding): reported pretrend IF columns.
cpt_pipeline <- function(data, est, frame = NULL, sel = list(), bp = "varying") {
  g <- function(nm, d) if (is.null(sel[[nm]])) d else sel[[nm]]
  gate <- chained_request_scope_gate(mode = NULL, panel = TRUE, allow_unbalanced_panel = NULL, est_method = est,
    selection = g("selection", "none"), selection_model = g("selection_model", "pair_logit"),
    selection_pooling = g("selection_pooling", "separate"), selection_xformla = sel$selection_xformla,
    selection_universe = g("selection_universe", "F"), selection_rho = sel$selection_rho,
    selection_justification = sel$selection_justification, selection_est_method = "ipw", imputation_pooling = "separate",
    pt_assumption = "post", pt_start = NULL, gmm_weighting = c("optimal", "identity"), gmm_covariance = "tikhonov",
    fix_weights = NULL, compute_inffunc = TRUE, pl = FALSE, cores = 1, faster_mode = TRUE, print_details = FALSE,
    control_group = "notyettreated", base_period = bp, anticipation = 0, alp = 0.05, bstrap = FALSE, cband = FALSE,
    biters = 100, dots = list(), entry = "public")
  resolved <- chained_public_resolve_control(NULL, gmm_weighting = gate$gmm_weighting_requested,
    gmm_covariance = gate$gmm_covariance)
  prep <- chained_request_prepare_inputs(data = data, frame = frame, yname = "y", tname = "time", idname = "id",
    gname = "G", xformla = ~x + z, weightsname = "w", clustervars = NULL, anticipation = 0, target_cohorts = NULL,
    selection = gate$selection, selection_model = g("selection_model", "pair_logit"),
    selection_pooling = g("selection_pooling", "separate"), selection_xformla = sel$selection_xformla,
    selection_universe = g("selection_universe", "F"), selection_rho = sel$selection_rho,
    selection_justification = sel$selection_justification, outcome_regression = identical(est, "reg"),
    or_resource = resolved$profile, doubly_robust = identical(est, "dr"), selection_est_method = "ipw",
    imputation_pooling = NULL, imputation_est_method = NULL, fix_weights = NULL)
  dims <- chained_public_dimensions(prep, biters = 100, clustered = FALSE)
  produced <- chained_public_pipeline_control(dims, resolved, gmm_weighting = gate$gmm_weighting_requested,
    gmm_covariance = gate$gmm_covariance, workers = 1)
  res <- run_chained_ipw_pipeline(data = prep$data_internal, idname = "id", tname = "time", yname = "y", gname = "G",
    control = produced$control, xformla = prep$xformla, cluster_vec = prep$cluster_full,
    target_cohorts = prep$target_cohorts, anticipation = 0L, control_group = "notyettreated", base_period = bp,
    weighting = gate$gmm_weighting_requested, pt_assumption = "post", run_numerical_transfer = FALSE,
    root_witness = NULL, aggregations = character(0L), cband = FALSE, alp = 0.05, biters = 100L, btype = "rademacher",
    seed = NULL, cp04 = prep$cp04, frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w,
    frame_g = prep$frame_g_arg, allow_frame_only = prep$n_frame_only > 0L, bstrap = FALSE, caller_rng = TRUE,
    custom_est_method = if (is.function(est)) est else NULL, custom_dots = list(),
    gmm_covariance = gate$gmm_covariance, compute_inffunc = TRUE)
  list(res = res, col = function(id) {
    j <- match(id, colnames(res$reporting$if_full))
    list(phi = unname(res$reporting$if_full[, j]), est = unname(res$reporting$estimates[j]))
  })
}

cpt_b5_setup <- function() {
  B <- cpt_csv("fB.csv")
  fx <- cpt_rds("priv.rds")$fixtures$fixtures$B
  fr <- B[!duplicated(B$id), c("id", "G", "x", "z", "sx", "sz", "w")]
  fr$rhoB <- as.double(unlist(fx$rho$B))[match(fr$id, sort(unique(B$id)))]
  list(data = B, frame = fr, sel = list(selection = "conditional_mtar", selection_model = "pair_logit",
    selection_pooling = "separate", selection_universe = "B", selection_rho = "rhoB", selection_xformla = ~sx + sz,
    selection_justification = "fixture B known normalized transport ratio"))
}

test_that("O15: explicit-gradient callback rows reproduce the accepted callback references", {
  priv <- cpt_rds("priv.rds"); B <- cpt_csv("fB.csv"); b5 <- cpt_b5_setup()
  names_cb <- c("ipw_builtin_equivalent", "h1_violating_witness", "h1_violating_sqrt_ps_witness"); pw <- c(1, 0, 0.5)
  cfgs <- list(B7_nyt_baseline_pair_rates = list(data = B, frame = NULL, sel = list()),
    B5_nyt_cor1_pair_B_separate = b5)
  for (cfg in names(cfgs)) {
    arr <- priv$arrays[[paste0(cfg, "|base")]]; refs <- priv$reference_values$configs[[cfg]]$callbacks
    for (bp in c("varying", "universal")) {
      id <- if (bp == "varying") "pre(g=3,t=2)" else "pre(g=4,t=1)"
      key <- if (bp == "varying") "varying__pre_var_g3_t2" else "universal__pre_uni_g4_t1"
      rkey <- if (bp == "varying") "varying|pre_var_g3_t2" else "universal|pre_uni_g4_t1"
      sgn <- if (bp == "varying") 1 else -1
      for (k in seq_along(pw)) {
        fit <- cpt_pipeline(cfgs[[cfg]]$data, cpt_callback(pw[k]), cfgs[[cfg]]$frame, cfgs[[cfg]]$sel, bp)
        cc <- fit$col(id)
        info <- paste(cfg, bp, names_cb[k])
        expect_lt(max(abs(sgn * cc$phi - arr[[paste0(key, "__cb-", names_cb[k], "__explicit")]])), 1e-10)
        expect_equal(sgn * cc$est, refs[[paste0(rkey, "|cb:", names_cb[k])]]$T_unoriented_at_c1, tolerance = 1e-12, info = info)
        # Explicit g with psi/m added again (double rate channel) is a mutation that must fail.
        dbl <- arr[[paste0(key, "__cb-", names_cb[k], "__mut-double_rate_channel_explicit_plus_identity")]]
        if (max(abs(dbl - arr[[paste0(key, "__cb-", names_cb[k], "__explicit")]])) > 1e-6)
          expect_gt(max(abs(sgn * cc$phi - dbl)), 1e-6)
      }
      # Without g: the baseline keeps the existing identity rule g = psi/m; the
      # estimated-observation route refuses (g is never inferred).
      if (identical(cfg, "B7_nyt_baseline_pair_rates")) {
        fit0 <- cpt_pipeline(B, cpt_callback(0, FALSE), NULL, list(), bp)
        expect_lt(max(abs(sgn * fit0$col(id)$phi -
          arr[[paste0(key, "__cb-h1_violating_witness__mut-identity_helper_psi_over_m")]])), 1e-10)
      } else {
        expect_error(cpt_pipeline(b5$data, cpt_callback(0, FALSE), b5$frame, b5$sel, bp),
          class = "cdid_custom_ipw_derivative_error")
      }
    }
  }
})

test_that("O15: a built-in-equivalent callback reproduces built-in CP04 rows on the production registry", {
  B <- cpt_csv("fB.csv")
  sel <- list(selection = "conditional_mtar", selection_model = "pair_logit", selection_pooling = "across_cohorts",
    selection_xformla = ~sx + sz)
  a <- cpt_pipeline(B, "ipw", NULL, sel); b <- cpt_pipeline(B, cpt_callback(1), NULL, sel)
  for (id in c("pre(g=3,t=2)", "pre(g=4,t=2)", "pre(g=4,t=3)")) {
    expect_lt(max(abs(a$col(id)$phi - b$col(id)$phi)), 1e-9)
    expect_equal(a$col(id)$est, b$col(id)$est, tolerance = 1e-12)
  }
  expect_true(inherits(b$res$nuisance_registry$cp04, "cdid_cp04_custom_backend"))
  expect_identical(b$res$nuisance_registry$cp04$treatments[[1L]]$derivative_source, "explicit_log_weight_gradient")
})

test_that("O15: full-F finite perturbations, including non-observers, match the composed rows", {
  # Sampling-weight perturbations of a unit change every full-F estimating
  # equation it enters. A non-observer enters only through the selection/rate
  # upstream (b' xi_u, any callback); an observer's i.weight also changes, which
  # equals the callback's fixed-weight case perturbation psi/m exactly under H1.
  b5 <- cpt_b5_setup(); d <- b5$data
  pat <- tapply(!is.na(d$y), list(d$id, d$time), any)
  obs12 <- pat[, "1"] & pat[, "2"]
  ids <- as.integer(rownames(pat)); G <- d$G[match(ids, d$id)]
  non_obs <- ids[G %in% c(3, 4, 0) & !obs12 & rowSums(pat) >= 1][1:2]
  obs <- ids[G == 0 & obs12][1L]
  N <- length(ids); h <- 1e-5
  for (case in list(list(cb = cpt_callback(0.5), units = non_obs), list(cb = cpt_callback(1), units = c(non_obs[1L], obs)))) {
    base <- cpt_pipeline(d, case$cb, b5$frame, b5$sel)$col("pre(g=3,t=2)")
    for (i in case$units) {
      T_at <- function(eps) {
        dd <- d; dd$w[dd$id == i] <- dd$w[dd$id == i] * (1 + eps)
        ff <- b5$frame; ff$w[ff$id == i] <- ff$w[ff$id == i] * (1 + eps)
        cpt_pipeline(dd, case$cb, ff, b5$sel)$col("pre(g=3,t=2)")$est
      }
      fd <- N * (T_at(h) - T_at(-h)) / (2 * h)
      expect_equal(base$phi[match(i, ids)], fd, tolerance = 1e-5, info = paste("unit", i))
      expect_true(abs(fd) > 1e-6)
    }
  }
})

test_that("O15: the optional gradient is validated before any use", {
  ok <- list(ATT = 1, att.inf.func = c(1, -1, 0))
  expect_null(custom_ipw_validate_callback_output(ok, 3)$g)
  expect_identical(custom_ipw_validate_callback_output(c(ok, list(log_weight_gradient = matrix(c(1, 2, 3), 3L))), 3)$g, c(1, 2, 3))
  expect_error(custom_ipw_validate_callback_output(c(ok, list(log_weight_gradient = c(1, 2))), 3),
    class = "cdid_custom_ipw_dimension_error")
  expect_error(custom_ipw_validate_callback_output(c(ok, list(log_weight_gradient = c(1, NA, 3))), 3),
    class = "cdid_custom_ipw_output_error")
  expect_error(custom_ipw_validate_callback_output(c(ok, list(log_weight_gradient = "a")), 3),
    class = "cdid_custom_ipw_output_error")
  val <- custom_ipw_validate_callback_output(c(ok, list(log_weight_gradient = c(0.5, 0.25, 1))), 3)
  # b_k = -sum_{j in k} g_j / r_k (explicit) versus -sum psi_j / (m r_k) (identity).
  expect_equal(custom_ipw_rate_sensitivity(val, c(1L, 2L, 2L), c(0.5, 0.8), 2L, 3L), c(-1, -1.25 / 0.8))
  expect_equal(custom_ipw_rate_sensitivity(custom_ipw_validate_callback_output(ok, 3), c(1L, 2L, 2L), c(0.5, 0.8), 2L, 3L),
    c(-1 / 1.5, 1 / 2.4))
})


test_that("global degeneracy uses cluster process rows rather than realized cancelling signs", {
  d <- cpt_synthetic(0)
  S <- matrix(rep(c(1, -1, -1, 1), length.out = 9L * 40L), 9L, 40L)
  r <- cpt_fixed(d, list(xformla = ~x, weightsname = NULL, control_group = "nevertreated"), S)
  expect_identical(r$chained$status, "available")
  # These signs cancel a nonzero process. Summation/QR roundoff need not
  # produce bitwise zero; the squared statistic is on the eps^2 scale.
  cancellation_bound <- 100 * 40 * .Machine$double.eps^2
  expect_lt(max(abs(r$CvMb)), cancellation_bound)
  expect_lt(abs(r$CvM), cancellation_bound)
  cell <- attr(r, "cells")[[1L]]
  expect_equal(crossprod(cpt_cell_process_rows(cell)), matrix(c(20, 20, 20, 40), 2L), tolerance = 1e-12)
  expect_true(cpt_process_nonzero(list(cell), 40L, 1L, 3L))
  zero <- cell
  zero$A_v[] <- 0; zero$AggX[] <- 0; zero$AggK[] <- 0
  expect_false(cpt_process_nonzero(list(zero), 40L, 1L, 3L))
})

test_that("NULL and intercept-only public grids retain a nondegenerate scalar process", {
  d <- cpt_csv("fA.csv")
  call_it <- function(fm) {
    set.seed(7)
    conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G", data = d,
      xformla = fm, control_group = "notyettreated", weightsname = "w", biters = 19)
  }
  a <- call_it(NULL); b <- call_it(~1)
  expect_identical(a$chained$status, "available")
  expect_identical(a$CvM, b$CvM); expect_identical(a$CvMb, b$CvMb)
  expect_true(all(a$chained$cells$n_distinct_grid_points == 1L))
  expect_true(any(a$CvMb > 0))
  expect_equal(a$CvM, 48 * sum(a$chained$cells$estimate^2), tolerance = 1e-12)
})

test_that("conditional bounds retain full grids, diagnostic rows and parallel cell copies", {
  dims <- list(list(Md = 20, p = 10, R = 2, Dl = 6, chunk = 5, rows_retained = 100 * 9))
  ctrl <- list(owner_control = list(max_retained_entries = 1e7, max_matrix_entries = 1e7, max_work_ops = 1e12))
  mb <- list(max_draws = 100, max_workspace_cells = 1e7, max_workspace_bytes = 8e7)
  a <- cpt_preflight(dims, 100, 100, 20, 10, 5, 5, ctrl, mb, 0, 0, 1)
  b <- cpt_preflight(dims, 100, 100, 20, 10, 5, 5, ctrl, mb, 0, 0, 2)
  cells <- 20 * (2 + 6 + 10 + 3) + 100 * (2 + 6) + 300 + 6 * 2 + 64 + 900
  expect_gte(b$multiplier_entries - a$multiplier_entries, 4 * cells)
  no_rows <- dims; no_rows[[1L]]$rows_retained <- 0
  z <- cpt_preflight(no_rows, 100, 100, 20, 10, 5, 5, ctrl, mb, 0, 0, 1)
  expect_equal(a$retained_entries - z$retained_entries, 900)
  mb$max_workspace_cells <- b$multiplier_entries - 1
  expect_error(cpt_preflight(dims, 100, 100, 20, 10, 5, 5, ctrl, mb, 0, 0, 2),
    class = "cdid_conditional_pretest_bound_error")
})


test_that("transformed factors use frozen model-frame levels and preserve unsupported ids", {
  d <- cpt_csv("fBd.csv"); d$x <- as.double(d$x > stats::median(d$x))
  r <- conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G", data = d,
    xformla = ~factor(x) + z, control_group = "notyettreated", weightsname = "w", biters = 19,
    est_method = "dr", fix_weights = "base_period", selection = "conditional_mtar",
    selection_model = "independent_period_logit", selection_pooling = "within_cohort", selection_xformla = ~sx + sz)
  expect_identical(r$chained$status, "available")
  expect_true(all(vapply(r$chained$grid$columns, function(x) "factor(x)1" %in% x, logical(1L))))
  d <- cpt_csv("fA.csv"); d <- d[d$G != 1, ]
  d$y[d$id == 5 & d$time > 1] <- NA
  d$grp <- ifelse(d$id %% 2 == 0, "a", "b"); d$grp[d$id == 5] <- "c"
  call_it <- function(d, fm) conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G",
    data = d, xformla = fm, control_group = "notyettreated", weightsname = "w", biters = 19, est_method = "reg")
  bad <- call_it(d, ~x + z + factor(grp))
  expect_identical(bad$chained$status, "unavailable")
  expect_identical(unique(bad$chained$grid$unsupported$id), 5L)
  expect_identical(unique(bad$chained$grid$unsupported$variable), "factor(grp)")
  expect_match(bad$chained$reason, "id 5", fixed = TRUE)
  # A declared unused factor column is already in the actual frozen design;
  # preserving it is distinct from introducing a new level during prediction.
  d$grp <- factor(d$grp)
  declared <- call_it(d, ~x + z + grp)
  expect_identical(declared$chained$status, "available")
  expect_true(all(vapply(declared$chained$grid$columns, function(x) "grpc" %in% x, logical(1L))))
})

test_that("non-fitting data-sourced grid rows require exact reference-date values", {
  d <- cpt_csv("fA.csv"); d <- d[d$G != 1, ]
  d$y[d$id == 5 & d$time > 1] <- NA
  d$x[d$id == 5] <- 1; d$x[d$id == 5 & d$time == 1] <- NA
  r <- conditional_did_pretest(yname = "y", tname = "time", idname = "id", gname = "G", data = d,
    xformla = ~x + z, control_group = "notyettreated", weightsname = "w", biters = 19, est_method = "reg")
  expect_identical(r$chained$status, "unavailable")
  expect_identical(unique(r$chained$grid$missing$id), 5L)
  expect_identical(unique(r$chained$grid$missing$date), 1)
  src <- list(frame_supplied = TRUE, frame = data.frame(x = c(4, 5)), frame_ids = 1:2,
    unit = c(1, 1, 2, 2), tt = c(1, 2, 1, 2), data = data.frame(x = c(4, 4, NA, 5)))
  expect_identical(cpt_grid_static_raw(src, "x", 1, c(TRUE, FALSE))$df$x, c(4, 5))
})
