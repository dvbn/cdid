# Bounded regression suite executed unconditionally by tests/testthat.R during
# standard R CMD check. Covers all supported public estimators, options,
# inference modes, ordinary-mode delegation, and refusal gates against
# accepted independent finite-sample references and frozen oracle fixtures.
# The complete 172-file development and Monte Carlo suite remains packaged in
# tests/testthat/ and runs unfiltered (filter = NULL) in CI and via
# validation/run_full_tests.R.

cran_ref_env <- new.env(parent = globalenv())
sys.source(testthat::test_path("reference", "cp03_public_tiny_reference.R"), envir = cran_ref_env)
sys.source(testthat::test_path("reference", "cp06_public_or_reference.R"), envir = cran_ref_env)
sys.source(testthat::test_path("reference", "cp06_pooled_public_or_reference.R"), envir = cran_ref_env)
sys.source(testthat::test_path("reference", "selection_public_reference.R"), envir = cran_ref_env)
sys.source(testthat::test_path("reference", "cp07_public_dr_reference.R"), envir = cran_ref_env)
sys.source(testthat::test_path("reference", "imputation_public_reference.R"), envir = cran_ref_env)
sys.source(testthat::test_path("reference", "d60_or_reference.R"), envir = cran_ref_env)
sys.source(testthat::test_path("reference", "d60_routes_reference.R"), envir = cran_ref_env)

cran_cp06_raw <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_baseline", "raw.csv"))
cran_cp06_raw$cohort <- as.numeric(replace(cran_cp06_raw$cohort, cran_cp06_raw$cohort == "NT", "Inf"))
cran_cp06_trans <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_corollary1", "transport.csv"))
cran_cp06_trans$cohort <- as.numeric(replace(cran_cp06_trans$cohort, cran_cp06_trans$cohort == "NT", "Inf"))

cran_cp06_fixture <- function(route = "none",
                              selection = if (identical(route, "none")) "none" else "conditional_mtar",
                              selection_model = "pair_logit",
                              selection_pooling = "separate",
                              est_method = "dr",
                              outcome_pooling = "pooled",
                              xformla = ~x,
                              selection_xformla = ~x) {
  d <- cran_cp06_raw
  N <- nrow(d)
  dates <- 1:4
  Y <- as.matrix(d[paste0("Y", dates)])
  ids <- paste0("unit-", d$id)
  g_mapped <- ifelse(is.infinite(d$cohort), 0, d$cohort)
  panel <- data.frame(
    id = rep(ids, each = length(dates)),
    time = rep(dates, N),
    y = as.vector(t(Y)),
    g = rep(g_mapped, each = length(dates)),
    stringsAsFactors = FALSE
  )
  # Retain one zero-outcome unit in the long panel as all-NA and make the
  # remaining zero-outcome units frame-only; singleton units (n_obs == 1) stay in panel.
  zero <- which(rowSums(!is.na(Y)) == 0L)
  panel <- panel[!panel$id %in% ids[zero[-1L]], , drop = FALSE]
  frame <- data.frame(
    id = ids,
    g = g_mapped,
    x = d$x,
    w = d$w,
    cl = paste0("cl-", (d$id - 1L) %/% 3L),
    stringsAsFactors = FALSE
  )
  if (route %in% c("B", "A")) {
    tb <- cran_cp06_trans[cran_cp06_trans$route == route, ]
    frame$rho <- tb$rho[match(paste(d$cohort, d$x), paste(tb$cohort, tb$x))]
  }
  args <- list(
    outcome_pooling = outcome_pooling,
    yname = "y", tname = "time", idname = "id", gname = "g",
    data = panel, frame = frame, xformla = xformla, weightsname = "w",
    est_method = est_method, control_group = "nevertreated",
    bstrap = FALSE, cband = FALSE, gmm_weighting = "identity"
  )
  if (!identical(route, "none")) {
    args$selection <- selection
    args$selection_model <- selection_model
    args$selection_pooling <- selection_pooling
    args$selection_xformla <- selection_xformla
    args$selection_universe <- route
    if (route %in% c("B", "A")) {
      args$selection_rho <- "rho"
      args$selection_justification <- "frozen population-law normalized transport ratio"
    }
  }
  list(args = args, frame = frame, data = panel, raw = d)
}

# Independent per-block WLS + logit reference for outcome_pooling = "separate"
# on the 144-unit CP06 fixture (base R stats::lm.wfit and stats::glm.fit only).
cran_separate_eval <- function(raw, w = raw$w, est_method = c("reg", "dr"),
                               route = "none", selection_model = "pair_logit",
                               selection_pooling = "separate",
                               Q_slope = cbind(x = raw$x),
                               X_slope = cbind(x = raw$x),
                               V_treat = cbind(`(Intercept)` = 1, x = raw$x),
                               dates = 1:4, control_group = "nevertreated",
                               anticipation = 0L, pt_assumption = "post",
                               pt_start = NULL, base_period = "varying",
                               cat_info = NULL) {
  est_method <- match.arg(est_method)
  if (is.null(cat_info)) {
    cat_info <- cran_ref_env$cp06_pooled_public_catalog(
      raw, dates = dates, control_group = control_group,
      anticipation = anticipation, pt_assumption = pt_assumption,
      pt_start = pt_start, base_period = base_period
    )
  }
  cal <- cat_info$calendar
  links_df <- cat_info$links_df
  N <- nrow(raw)
  idx_time <- function(tt) match(tt, cal)
  R_mat <- !is.na(as.matrix(raw[paste0("Y", seq_along(cal))]))
  A_vec <- as.numeric(rowSums(R_mat) >= 2L)
  B_vec <- as.numeric(rowSums(R_mat) >= 1L)
  cohorts_all <- cat_info$cohorts_all
  K_len <- length(cohorts_all)
  mass_F <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k]), numeric(1))
  mass_A <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k & A_vec == 1]), numeric(1))
  tau <- mass_F / sum(mass_F)
  names(tau) <- as.character(cohorts_all)

  if (identical(route, "none")) {
    get_q_weights <- function(rows, is, it, kv) {
      pair_rows <- raw$cohort == kv & R_mat[, is] & R_mat[, it]
      rate <- sum(w[pair_rows]) / mass_A[match(kv, cohorts_all)]
      w[rows] / rate
    }
  } else {
    U_vec <- if (route == "F") rep(1, N) else if (route == "B") B_vec else A_vec
    rho_vec <- if (route %in% c("B", "A")) raw$rho else rep(1, N)
    mass_U <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k & U_vec == 1]), numeric(1))
    a_U <- ifelse(mass_U > 0, mass_A / mass_U, 0)
    names(a_U) <- as.character(cohorts_all)
    pc_list <- list()
    for (r in seq_len(nrow(links_df))) {
      is <- links_df$is[r]; it <- links_df$it[r]
      ks <- c(links_df$g[r], sort(as.numeric(strsplit(links_df$controls[r], "|", fixed = TRUE)[[1L]])))
      for (kv in ks) {
        pc <- list(k_idx = match(kv, cohorts_all), cohort = kv, is = is, it = it)
        if (!any(vapply(pc_list, identical, logical(1), pc))) pc_list[[length(pc_list) + 1L]] <- pc
      }
    }
    pc_list <- pc_list[order(
      vapply(pc_list, `[[`, integer(1), "k_idx"),
      vapply(pc_list, `[[`, integer(1), "is"),
      vapply(pc_list, `[[`, integer(1), "it")
    )]
    cells <- list()
    if (identical(selection_model, "independent_period_logit")) {
      for (k_idx in seq_len(K_len)) for (t_idx in seq_along(cal)) {
        used <- any(vapply(pc_list, function(pc) pc$k_idx == k_idx && (pc$is == t_idx || pc$it == t_idx), logical(1)))
        if (used) cells[[length(cells) + 1L]] <- list(kind = "date", k_idx = k_idx, cohort = cohorts_all[k_idx], t_idx = t_idx)
      }
    } else {
      up_pairs <- unique(data.frame(
        is = vapply(pc_list, `[[`, integer(1), "is"),
        it = vapply(pc_list, `[[`, integer(1), "it")
      ))
      up_pairs <- up_pairs[order(up_pairs$is, up_pairs$it), , drop = FALSE]
      for (p_i in seq_len(nrow(up_pairs))) {
        is <- up_pairs$is[p_i]; it <- up_pairs$it[p_i]
        for (k_idx in seq_len(K_len)) {
          used <- any(vapply(pc_list, function(pc) pc$k_idx == k_idx && pc$is == is && pc$it == it, logical(1)))
          if (used) cells[[length(cells) + 1L]] <- list(kind = "pair", k_idx = k_idx, cohort = cohorts_all[k_idx], is = is, it = it)
        }
      }
    }
    C_cells <- length(cells)
    block_groups <- if (identical(selection_pooling, "separate")) {
      lapply(seq_len(C_cells), function(cid) cid)
    } else if (identical(selection_pooling, "within_cohort")) {
      active_k <- sort(unique(vapply(cells, `[[`, integer(1), "k_idx")))
      lapply(active_k, function(k_idx) which(vapply(cells, `[[`, integer(1), "k_idx") == k_idx))
    } else {
      list(seq_len(C_cells))
    }
    Z_sl <- as.matrix(Q_slope)
    d_s <- ncol(Z_sl)
    cell_prob <- matrix(0, N, C_cells)
    for (m in seq_along(block_groups)) {
      cids <- block_groups[[m]]; B_m <- length(cids)
      X_rows <- list(); y_vec <- numeric(0); w_vec <- numeric(0)
      for (b in seq_len(B_m)) {
        cid <- cids[b]; cell <- cells[[cid]]
        ik <- which(raw$cohort == cell$cohort & U_vec == 1)
        y_c <- if (identical(cell$kind, "date")) as.double(R_mat[ik, cell$t_idx]) else as.double(R_mat[ik, cell$is] & R_mat[ik, cell$it])
        E_b <- matrix(0, length(ik), B_m); E_b[, b] <- 1.0
        X_rows[[b]] <- cbind(E_b, Z_sl[ik, , drop = FALSE])
        y_vec <- c(y_vec, y_c); w_vec <- c(w_vec, w[ik])
      }
      coef_blk <- cran_ref_env$cp07_public_dr_fit_logit(do.call(rbind, X_rows), y_vec, w_vec)
      al_m <- coef_blk[seq_len(B_m)]
      be_m <- if (d_s > 0L) coef_blk[B_m + seq_len(d_s)] else numeric(0)
      for (b in seq_len(B_m)) {
        cid <- cids[b]
        cell_prob[, cid] <- stats::plogis(al_m[b] + if (d_s > 0L) drop(Z_sl %*% be_m) else 0.0)
      }
    }
    get_pair_prob <- function(is, it, kv) {
      k_idx <- match(kv, cohorts_all)
      if (identical(selection_model, "independent_period_logit")) {
        cid_s <- which(vapply(cells, function(c) c$k_idx == k_idx && c$t_idx == is, logical(1)))
        cid_t <- which(vapply(cells, function(c) c$k_idx == k_idx && c$t_idx == it, logical(1)))
        cell_prob[, cid_s] * cell_prob[, cid_t]
      } else {
        cid_p <- which(vapply(cells, function(c) c$k_idx == k_idx && c$is == is && c$it == it, logical(1)))
        cell_prob[, cid_p]
      }
    }
    get_q_weights <- function(rows, is, it, kv) {
      pp <- get_pair_prob(is, it, kv)[rows]
      w[rows] * a_U[as.character(kv)] / (rho_vec[rows] * pp)
    }
  }

  # Separate outcome regression per unique (pair, control_pool) block
  upairs <- unique(data.frame(
    s = links_df$s, t = links_df$t,
    is = idx_time(links_df$s), it = idx_time(links_df$t),
    stringsAsFactors = FALSE
  ))
  upairs <- upairs[order(upairs$is, upairs$it), , drop = FALSE]
  upairs$pair_id <- seq_len(nrow(upairs)); upairs$L <- upairs$it - upairs$is
  link_pair_id <- match(paste0(links_df$s, "_", links_df$t), paste0(upairs$s, "_", upairs$t))
  block_df <- unique(data.frame(pair_id = link_pair_id, controls = links_df$controls, stringsAsFactors = FALSE))
  pool_sort_key <- vapply(block_df$controls, function(cs) {
    ks <- match(as.numeric(strsplit(cs, "|", fixed = TRUE)[[1L]]), cohorts_all)
    paste(sprintf("%012d", ks), collapse = "|")
  }, character(1))
  block_df <- block_df[order(sprintf("%012d|%s", block_df$pair_id, pool_sort_key), method = "radix"), , drop = FALSE]
  G <- nrow(block_df)
  link_block_id <- match(paste0(link_pair_id, "|", links_df$controls), paste0(block_df$pair_id, "|", block_df$controls))
  p_or <- ncol(X_slope)
  alpha <- numeric(G)
  beta_mat <- matrix(0, G, p_or)
  for (b in seq_len(G)) {
    pid <- block_df$pair_id[b]
    is <- upairs$is[pid]; it <- upairs$it[pid]; L_b <- upairs$L[pid]
    ctrl_cohs <- as.numeric(strsplit(block_df$controls[b], "|", fixed = TRUE)[[1L]])
    XX <- NULL; YY <- numeric(0); WW <- numeric(0)
    for (kv in ctrl_cohs) {
      rows <- which(raw$cohort == kv & R_mat[, is] & R_mat[, it])
      qw <- get_q_weights(rows, is, it, kv)
      XX <- rbind(XX, cbind(1, L_b * X_slope[rows, , drop = FALSE]))
      YY <- c(YY, raw[[paste0("Y", it)]][rows] - raw[[paste0("Y", is)]][rows])
      WW <- c(WW, qw)
    }
    wfit <- stats::lm.wfit(x = XX, y = YY, w = WW, tol = 1e-12)
    coef_b <- unname(wfit$coefficients)
    alpha[b] <- coef_b[1L]
    beta_mat[b, ] <- coef_b[-1L]
  }

  L_links <- nrow(links_df)
  mu <- numeric(L_links)
  for (l in seq_len(L_links)) {
    b <- link_block_id[l]; pid <- link_pair_id[l]
    is <- upairs$is[pid]; it <- upairs$it[pid]; L_l <- upairs$L[pid]
    g_val <- links_df$g[l]
    if (identical(est_method, "reg")) {
      rows <- which(raw$cohort == g_val & R_mat[, is] & R_mat[, it])
      qw <- get_q_weights(rows, is, it, g_val)
      dy <- raw[[paste0("Y", it)]][rows] - raw[[paste0("Y", is)]][rows]
      pred <- alpha[b] + as.vector(L_l * X_slope[rows, , drop = FALSE] %*% beta_mat[b, ])
      mu[l] <- sum(qw * (dy - pred)) / sum(qw)
    } else {
      ctrl_cohs <- sort(as.numeric(strsplit(links_df$controls[l], "|", fixed = TRUE)[[1L]]))
      cohs <- c(g_val, ctrl_cohs)
      obs_rows <- which(raw$cohort %in% cohs & R_mat[, is] & R_mat[, it])
      D <- as.double(raw$cohort[obs_rows] == g_val)
      dy <- raw[[paste0("Y", it)]][obs_rows] - raw[[paste0("Y", is)]][obs_rows]
      pred <- alpha[b] + as.vector(L_l * X_slope[obs_rows, , drop = FALSE] %*% beta_mat[b, ])
      res <- dy - pred
      fw <- numeric(length(obs_rows))
      for (kv in cohs) {
        jj <- which(raw$cohort[obs_rows] == kv)
        fw[jj] <- get_q_weights(obs_rows[jj], is, it, kv)
      }
      X_tr <- V_treat[obs_rows, , drop = FALSE]
      b_tr <- cran_ref_env$cp07_public_dr_fit_logit(X_tr, D, fw)
      odds <- exp(drop(X_tr %*% b_tr))
      wp <- fw * D; wm <- fw * (1 - D) * odds
      mu[l] <- sum(wp * res) / sum(wp) - sum(wm * res) / sum(wm)
    }
  }
  cells_out <- cran_ref_env$cp06_pooled_public_cells_from_mu(
    cat_info, mu, anticipation = anticipation,
    pt_assumption = pt_assumption, pt_start = pt_start, base_period = base_period
  )
  agg_overall <- cran_ref_env$imputation_public_aggregate(cells_out, tau)
  list(mu = mu, cells = cells_out, att = cells_out$att, agg_overall = agg_overall, cat_info = cat_info)
}

cran_separate_reference <- function(raw, est_method = "reg", route = "none",
                                    selection_model = "pair_logit",
                                    selection_pooling = "separate",
                                    dates = 1:4, control_group = "nevertreated",
                                    anticipation = 0L, pt_assumption = "post",
                                    pt_start = NULL, base_period = "varying",
                                    eps = 5e-5) {
  cat_info <- cran_ref_env$cp06_pooled_public_catalog(
    raw, dates = dates, control_group = control_group,
    anticipation = anticipation, pt_assumption = pt_assumption,
    pt_start = pt_start, base_period = base_period
  )
  ev0 <- cran_separate_eval(
    raw, w = raw$w, est_method = est_method, route = route,
    selection_model = selection_model, selection_pooling = selection_pooling,
    dates = dates, control_group = control_group, anticipation = anticipation,
    pt_assumption = pt_assumption, pt_start = pt_start, base_period = base_period,
    cat_info = cat_info
  )
  N <- nrow(raw)
  R <- length(ev0$att)
  avail <- !is.na(ev0$att)
  avail_agg <- !is.na(ev0$agg_overall)
  inffunc <- matrix(NA_real_, N, R)
  inffunc[, avail] <- 0
  phi_agg_overall <- matrix(NA_real_, N, length(ev0$agg_overall), dimnames = list(NULL, names(ev0$agg_overall)))
  phi_agg_overall[, avail_agg] <- 0
  for (i in seq_len(N)) {
    wp <- raw$w; wp[i] <- raw$w[i] * (1 + eps)
    wm <- raw$w; wm[i] <- raw$w[i] * (1 - eps)
    up <- cran_separate_eval(
      raw, w = wp, est_method = est_method, route = route,
      selection_model = selection_model, selection_pooling = selection_pooling,
      dates = dates, control_group = control_group, anticipation = anticipation,
      pt_assumption = pt_assumption, pt_start = pt_start, base_period = base_period,
      cat_info = cat_info
    )
    um <- cran_separate_eval(
      raw, w = wm, est_method = est_method, route = route,
      selection_model = selection_model, selection_pooling = selection_pooling,
      dates = dates, control_group = control_group, anticipation = anticipation,
      pt_assumption = pt_assumption, pt_start = pt_start, base_period = base_period,
      cat_info = cat_info
    )
    inffunc[i, avail] <- N * (up$att[avail] - um$att[avail]) / (2 * eps)
    if (any(avail_agg)) {
      phi_agg_overall[i, avail_agg] <- N * (up$agg_overall[avail_agg] - um$agg_overall[avail_agg]) / (2 * eps)
    }
  }
  V_analytical <- matrix(NA_real_, R, R)
  if (any(avail)) {
    centered <- scale(inffunc[, avail, drop = FALSE], center = TRUE, scale = FALSE)
    V_analytical[avail, avail] <- crossprod(centered) / N
  }
  c(ev0, list(inffunc = inffunc, phi_agg_overall = phi_agg_overall, V_analytical = V_analytical))
}

cran_compare_mp_to_ref <- function(fit, ref, frame_ids, tol_att = 1e-11, tol_if = 1e-7,
                                   tol_v = 1e-7, tol_se = 1e-8, check_agg = TRUE) {
  N <- length(frame_ids)
  expect_s3_class(fit, "cdid_MP")
  expect_equal(fit$n, N)
  pk <- paste0(fit$group, ":", fit$t)
  rk <- paste0(ref$cells$g, ":", ref$cells$time)
  j <- match(rk, pk)
  expect_identical(length(pk), length(rk))
  expect_false(anyNA(j))
  IF <- as.matrix(fit$inffunc)
  expect_identical(nrow(IF), N)
  expect_setequal(rownames(IF), frame_ids)
  ri <- match(frame_ids, rownames(IF))
  av <- !is.na(ref$att)
  expect_identical(av, !is.na(fit$att[j]))
  a <- which(av)
  expect_lt(max(abs(fit$att[j[a]] - ref$att[a])), tol_att)
  expect_lt(max(abs(IF[ri, j[a], drop = FALSE] - ref$inffunc[, a, drop = FALSE])), tol_if)
  V <- as.matrix(fit$V_analytical)
  expect_equal(unname(V[j[a], j[a], drop = FALSE]), unname(ref$V_analytical[a, a, drop = FALSE]), tolerance = tol_v)
  nonref <- a[ref$cells$role[a] != "reference"]
  expected_se <- sqrt(diag(ref$V_analytical)[nonref] / N)
  expect_equal(unname(fit$se[j[nonref]]), unname(expected_se), tolerance = tol_se)
  structural <- a[ref$cells$role[a] == "reference"]
  if (length(structural)) {
    expect_identical(unname(fit$att[j[structural]]), rep(0, length(structural)))
    expect_true(all(is.na(fit$se[j[structural]]) | fit$se[j[structural]] == 0))
    expect_identical(max(abs(IF[, j[structural], drop = FALSE])), 0)
  }
  if (check_agg) {
    for (tp in c("simple", "group", "dynamic", "calendar")) {
      ag <- cdid::aggte(fit, type = tp, bstrap = FALSE, cband = FALSE)
      expect_lt(abs(ag$overall.att - unname(ref$agg_overall[tp])), tol_att)
      ifv <- switch(tp,
        simple = ag$inf.function$simple.att,
        group = ag$inf.function$selective.inf.func,
        dynamic = ag$inf.function$dynamic.inf.func,
        calendar = ag$inf.function$calendar.inf.func
      )
      ifv <- as.matrix(ifv)
      rows <- if (!is.null(rownames(ifv))) rownames(ifv) else rownames(IF)
      expect_identical(dim(ifv), c(N, 1L))
      expect_lt(max(abs(ifv[match(frame_ids, rows), 1L] - ref$phi_agg_overall[, tp])), tol_if)
      # Frozen imputation/augmented references include their cluster operator.
      # Other references use centered original-unit rows without clustering.
      se_ref <- if (!is.null(ref$se_agg)) unname(ref$se_agg[tp]) else
        sqrt(sum((ref$phi_agg_overall[, tp] - mean(ref$phi_agg_overall[, tp]))^2)) / N
      expect_equal(ag$overall.se, se_ref, tolerance = tol_se)
    }
  }
}

# ------------------------------------------------------------------------------
# 1. Public interface, native registration, ordinary delegation, refusal gates
# ------------------------------------------------------------------------------

test_that("exported interface, native entry points, ordinary delegation, and scope gates hold", {
  exports <- sort(getNamespaceExports("cdid"))
  expect_identical(exports, sort(c("att_gt", "aggte", "conditional_did_pretest", "ggdid", "tidy", "glance")))

  dlls <- getLoadedDLLs()
  expect_true("cdid" %in% names(dlls))
  reg <- getDLLRegisteredRoutines("cdid")
  expect_true("C_gmm_native_status" %in% names(reg[[".Call"]]))
  expect_true("C_gmm_iv_dot" %in% names(reg[[".Call"]]))
  native_status <- .Call(reg[[".Call"]][["C_gmm_native_status"]])
  expect_identical(native_status$backend, "native_c_ieee754_binary64")
  expect_true(native_status$c_registration)
  expect_true(native_status$force_symbols)

  # The parity helper requires captured value/conditions, not raw fit objects.
  p_bal <- parity_panel()
  ord_mp <- with_conditions(cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G", xformla = ~x,
    data = p_bal, bstrap = FALSE, cband = FALSE, mode = "ordinary",
    selection = stop("chained argument must not be evaluated in ordinary mode"),
    imputation_pooling = stop("chained argument must not be evaluated in ordinary mode"),
    outcome_pooling = stop("chained argument must not be evaluated in ordinary mode")
  ))
  ref_mp <- with_conditions(did::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G", xformla = ~x,
    data = p_bal, bstrap = FALSE, cband = FALSE
  ))
  expect_false(inherits(ord_mp$value, "cdid_MP"))
  expect_same_result(ord_mp, ref_mp, mp_fields)
  expect_equal(cdid::tidy(ord_mp$value), did::tidy(ref_mp$value))
  expect_equal(cdid::glance(ord_mp$value), did::glance(ref_mp$value))

  for (tp in c("simple", "group", "dynamic", "calendar")) {
    ord_ag <- with_conditions(cdid::aggte(ord_mp$value, type = tp, bstrap = FALSE, cband = FALSE))
    ref_ag <- with_conditions(did::aggte(ref_mp$value, type = tp, bstrap = FALSE, cband = FALSE))
    expect_same_result(ord_ag, ref_ag, aggte_fields)
  }

  # Repeated cross-sections (panel = FALSE) default to ordinary mode
  rc_dat <- parity_rc()
  rc_mp <- with_conditions(cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G", xformla = ~x,
    data = rc_dat, panel = FALSE, bstrap = FALSE, cband = FALSE
  ))
  rc_ref <- with_conditions(did::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G", xformla = ~x,
    data = rc_dat, panel = FALSE, bstrap = FALSE, cband = FALSE
  ))
  expect_same_result(rc_mp, rc_ref, mp_fields)

  # Refusal gates must fire before forcing the data promise
  unforced <- list(
    yname = "y", tname = "time", idname = "id", gname = "g",
    data = quote(stop("data promise must not be forced"))
  )
  refuse_before_data <- function(cls, ...) {
    err <- expect_error(do.call(cdid::att_gt, utils::modifyList(unforced, list(...))), class = cls)
    expect_no_match(conditionMessage(err), "data promise must not be forced")
    invisible(err)
  }
  refuse_before_data("cdid_mode_error", mode = "invalid_mode")
  refuse_before_data("cdid_chained_panel", panel = FALSE, mode = "chained")
  refuse_before_data("cdid_chained_balanced", allow_unbalanced_panel = FALSE, mode = "chained")
  refuse_before_data("cdid_chained_request_est_method", est_method = "stdipw")
  refuse_before_data("cdid_chained_request_unsupported_selection", selection = "invalid_sel")
  refuse_before_data("cdid_chained_request_unsupported_selection", selection_est_method = "dr", selection = "none")
  refuse_before_data("cdid_chained_request_input_error", selection = "conditional_mtar",
                     selection_universe = "B")
  refuse_before_data("cdid_chained_request_outcome_pooling", outcome_pooling = "invalid_pool")
  refuse_before_data("cdid_chained_request_input_error", gmm_weighting = "invalid_gmm")
  refuse_before_data("cdid_chained_request_input_error", gmm_covariance = "invalid_cov")
  refuse_before_data("cdid_chained_request_input_error", pt_assumption = "invalid_pt")
  refuse_before_data("cdid_chained_request_unsupported_weights", fix_weights = "invalid_anchor")
  refuse_before_data("cdid_chained_public_accuracy_deferred", control = list(relative_accuracy = 1e-6))
})

# ------------------------------------------------------------------------------
# 2. Weighted 2-period panel coincidence with upstream DRDID and separate/pooled
# ------------------------------------------------------------------------------

test_that("weighted 2-period balanced panel matches upstream DRDID across IPW, OR, and DR", {
  # In a 2-period balanced panel with one treated cohort (g = 2) and never-treated
  # controls (g = 0), there is a single comparison block (G = 1, L = 1) and unit
  # pair observation rate (r_k = 1). Chained IPW, OR, and DR specifications
  # coincide with DRDID::std_ipw_did_panel, DRDID::reg_did_panel, and DRDID::drdid_panel,
  # and outcome_pooling = "separate" vs "pooled" are identical.
  d <- cran_cp06_raw
  keep <- which(!is.na(d$Y1) & !is.na(d$Y2))
  d2 <- d[keep, , drop = FALSE]
  N <- nrow(d2)
  D <- as.numeric(is.finite(d2$cohort))
  g2 <- ifelse(D == 1, 2, 0)
  y0 <- d2$Y1
  y1 <- d2$Y2
  w <- d2$w
  X <- cbind(`(Intercept)` = 1, x = d2$x)
  panel2 <- data.frame(
    id = rep(seq_len(N), each = 2L),
    time = rep(c(1, 2), N),
    y = as.vector(rbind(y0, y1)),
    g = rep(g2, each = 2L),
    x = rep(d2$x, each = 2L),
    w = rep(w, each = 2L)
  )

  drdid_ipw <- DRDID::std_ipw_did_panel(y1 = y1, y0 = y0, D = D, covariates = X, i.weights = w, boot = FALSE, inffunc = TRUE)
  drdid_reg <- DRDID::reg_did_panel(y1 = y1, y0 = y0, D = D, covariates = X, i.weights = w, boot = FALSE, inffunc = TRUE)
  drdid_dr <- DRDID::drdid_panel(y1 = y1, y0 = y0, D = D, covariates = X, i.weights = w, boot = FALSE, inffunc = TRUE)

  refs <- list(ipw = drdid_ipw, reg = drdid_reg, dr = drdid_dr)
  for (em in c("ipw", "reg", "dr")) {
    mp_sep <- cdid::att_gt(
      yname = "y", tname = "time", idname = "id", gname = "g",
      xformla = ~x, weightsname = "w", data = panel2,
      est_method = em, outcome_pooling = "separate",
      bstrap = FALSE, cband = FALSE
    )
    ref <- refs[[em]]
    expect_equal(mp_sep$att[1L], as.numeric(ref$ATT), tolerance = 1e-11, info = em)
    phi_ref <- as.numeric(ref$att.inf.func)
    expect_equal(as.numeric(mp_sep$inffunc[, 1L]), phi_ref, tolerance = 1e-10, info = em)
    v_ref <- sum((phi_ref - mean(phi_ref))^2) / N
    expect_equal(as.numeric(mp_sep$V_analytical[1L, 1L]), v_ref, tolerance = 1e-10, info = em)
    expect_equal(mp_sep$se[1L], sqrt(v_ref / N), tolerance = 1e-10, info = em)

    if (em %in% c("reg", "dr")) {
      mp_pool <- cdid::att_gt(
        yname = "y", tname = "time", idname = "id", gname = "g",
        xformla = ~x, weightsname = "w", data = panel2,
        est_method = em, outcome_pooling = "pooled",
        bstrap = FALSE, cband = FALSE
      )
      expect_equal(mp_sep$att, mp_pool$att, tolerance = 1e-12, info = em)
      expect_equal(mp_sep$inffunc, mp_pool$inffunc, tolerance = 1e-12, info = em)
      expect_equal(mp_sep$V_analytical, mp_pool$V_analytical, tolerance = 1e-12, info = em)
      expect_equal(mp_sep$se, mp_pool$se, tolerance = 1e-12, info = em)
    }
  }
})

# ------------------------------------------------------------------------------
# 3. GMM weighting ("identity", default "optimal" two-step, active Tikhonov)
#    and PT-all vs PT-post with zero-outcome and singleton units
# ------------------------------------------------------------------------------

test_that("GMM weighting, active Tikhonov ridge, PT-all, and singleton/frame-only IF scaling match cp03 oracle", {
  d <- cran_ref_env$cp03_public_tiny_data()
  for (mode in c("identity", "two_step", "ridge")) {
    ref <- cran_ref_env$cp03_public_tiny_reference(mode)
    gw <- if (identical(mode, "identity")) "identity" else "optimal"
    ctrl <- if (identical(mode, "ridge")) list(condition_trigger = 0.99, ridge_relative = 1e-4) else NULL
    mp <- cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g", data = d,
      est_method = "ipw", gmm_weighting = gw, gmm_covariance = "tikhonov",
      pt_assumption = "all", pt_start = 2, base_period = "universal",
      control = ctrl, bstrap = FALSE, cband = FALSE
    )
    post <- which(mp$chained$role == "post")
    ref_cell <- which(mp$chained$role == "reference")
    expect_length(post, 1L)
    expect_length(ref_cell, 1L)
    expect_equal(mp$att[post], ref$point[2L], tolerance = 1e-10, info = mode)
    expect_equal(as.numeric(mp$inffunc[, post]), as.numeric(ref$IF[, 2L]), tolerance = 1e-10, info = mode)
    expect_equal(mp$V_analytical[post, post] / 8, ref$covariance[2L, 2L], tolerance = 1e-10, info = mode)
    expect_equal(mp$se[post], sqrt(ref$covariance[2L, 2L]), tolerance = 1e-10, info = mode)
    # Exact structural reference zero at the universal reference period
    expect_identical(mp$att[ref_cell], 0)
    expect_identical(max(abs(mp$inffunc[, ref_cell])), 0)
    if (identical(mode, "ridge")) {
      expect_true(mp$chained$covariance_handling$handler_invoked)
      expect_identical(mp$chained$covariance_handling$handler_used, "tikhonov")
      expect_equal(mp$chained$covariance_handling$lambda, ref$lambda, tolerance = 1e-10)
    } else {
      expect_false(isTRUE(mp$chained$covariance_handling$handler_invoked))
    }
  }

  # Adding one singleton unit (id = 9, n_obs = 1) and one frame-only unit (id = 10, n_obs = 0)
  # preserves point estimates, sets zero IF rows on non-eligible units for cell ATTs,
  # and scales original-unit IF rows by N_new / N_old = 10 / 8.
  one <- d[1L, , drop = FALSE]; one$id <- 9L
  panel10 <- rbind(d, one)
  frame10 <- data.frame(id = 1:10, g = c(rep(4, 4), rep(0, 4), 4, 0), cluster = c(rep(1:4, each = 2), 5, 5))
  mp10 <- cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g",
    data = panel10, frame = frame10, est_method = "ipw",
    pt_assumption = "all", pt_start = 2, bstrap = FALSE, cband = FALSE
  )
  mp10_cl <- cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g",
    data = panel10, frame = frame10, clustervars = "cluster", est_method = "ipw",
    pt_assumption = "all", pt_start = 2, bstrap = FALSE, cband = FALSE
  )
  expect_identical(mp10$n, 10L)
  expect_identical(sum(mp10$chained$frame$n_obs == 1L), 1L)
  expect_true(mp10$chained$frame$frame_only[10L])
  expect_identical(mp10$att, mp10_cl$att)
  expect_identical(mp10$inffunc, mp10_cl$inffunc)
  post10 <- which(mp10$chained$role == "post")
  ref_two <- cran_ref_env$cp03_public_tiny_reference("two_step")
  expect_equal(unname(mp10$inffunc[9:10, post10]), c(0, 0), tolerance = 0)
  expect_equal(unname(mp10$inffunc[1:8, post10]), ref_two$IF[, 2L] * (10 / 8), tolerance = 1e-10)
})

# ------------------------------------------------------------------------------
# 4. Selection IPW, separate/pooled OR, and DR across baseline and conditional MTAR
# ------------------------------------------------------------------------------

test_that("selection IPW matches sp_ref_public_reference across logit models, poolings, universes, and control groups", {
  # Case 1: pair_logit, within_cohort, universe = "B" with known rho, NYT, PT-all, universal base
  fx_b <- cran_cp06_fixture(
    route = "B", selection = "conditional_mtar",
    selection_model = "pair_logit", selection_pooling = "within_cohort",
    est_method = "ipw"
  )
  args_b <- utils::modifyList(fx_b$args, list(
    control_group = "notyettreated", pt_assumption = "all", pt_start = 1, base_period = "universal"
  ))
  fit_b <- do.call(cdid::att_gt, args_b)
  raw_b <- fx_b$raw; raw_b$rho <- fx_b$frame$rho
  ref_b <- cran_ref_env$sp_ref_public_reference(
    raw_b, est_method = "ipw", selection_model = "pair_logit",
    selection_pooling = "within_cohort", route = "B",
    control_group = "notyettreated", pt_assumption = "all", pt_start = 1, base_period = "universal"
  )
  cran_compare_mp_to_ref(fit_b, ref_b, fx_b$frame$id)

  # Case 2: independent_period_logit, across_cohorts, universe = "F", NT, varying base
  fx_f <- cran_cp06_fixture(
    route = "F", selection = "conditional_mtar",
    selection_model = "independent_period_logit", selection_pooling = "across_cohorts",
    est_method = "ipw"
  )
  fit_f <- do.call(cdid::att_gt, fx_f$args)
  ref_f <- cran_ref_env$sp_ref_public_reference(
    fx_f$raw, est_method = "ipw", selection_model = "independent_period_logit",
    selection_pooling = "across_cohorts", route = "F",
    control_group = "nevertreated", pt_assumption = "post", base_period = "varying"
  )
  cran_compare_mp_to_ref(fit_f, ref_f, fx_f$frame$id)

  # Case 3: legacy alias selection = "corollary1" on universe = "A" matches "conditional_mtar"
  fx_a_cm <- cran_cp06_fixture(
    route = "A", selection = "conditional_mtar",
    selection_model = "pair_logit", selection_pooling = "separate",
    est_method = "ipw"
  )
  fx_a_c1 <- cran_cp06_fixture(
    route = "A", selection = "corollary1",
    selection_model = "pair_logit", selection_pooling = "separate",
    est_method = "ipw"
  )
  fit_a_cm <- do.call(cdid::att_gt, fx_a_cm$args)
  fit_a_c1 <- do.call(cdid::att_gt, fx_a_c1$args)
  expect_equal(fit_a_cm$att, fit_a_c1$att, tolerance = 1e-12)
  expect_equal(fit_a_cm$inffunc, fit_a_c1$inffunc, tolerance = 1e-12)
  expect_equal(fit_a_cm$V_analytical, fit_a_c1$V_analytical, tolerance = 1e-12)
})

test_that("pooled and separate outcome regression (OR) match independent references across baseline and conditional MTAR", {
  # Pooled OR under baseline MTAR with NYT, anticipation = 1, universal base
  fx_or_pool <- cran_cp06_fixture(route = "none", est_method = "reg", outcome_pooling = "pooled")
  args_or_pool <- utils::modifyList(fx_or_pool$args, list(
    control_group = "notyettreated", anticipation = 1L, base_period = "universal"
  ))
  fit_or_pool <- do.call(cdid::att_gt, args_or_pool)
  ref_or_pool <- cran_ref_env$cp06_pooled_public_or_reference(
    fx_or_pool$raw, route = "none", control_group = "notyettreated",
    anticipation = 1L, base_period = "universal"
  )
  # Current pretrends are direct auxiliary comparisons at every t < g, including
  # the anticipation window ref_g < t < g (R/chained_ipw_pipeline.R:18-56).
  # The legacy oracle's cells_from_mu uses a GMM coordinate for that window
  # (reference/cp06_pooled_public_or_reference.R:173-175). Its independent mu
  # and full-unit phi_mu already contain the direct auxiliary comparisons.
  # Map ALL nonreference pre cells from those primitives, with no fitted output.
  ci <- ref_or_pool$cat_info
  for (r in which(ref_or_pool$cells$role == "pre")) {
    g <- ref_or_pool$cells$g[r]; u <- ref_or_pool$cells$time[r]
    base <- max(ci$calendar[ci$calendar + 1L < g])
    hit <- which(ci$aux_df$g == g & ci$aux_df$s == min(u, base) &
                   ci$aux_df$t == max(u, base))
    expect_length(hit, 1L)
    stopifnot(length(hit) == 1L)
    link <- nrow(ci$cat_df) + hit
    sign <- if (u < base) -1 else 1
    ref_or_pool$att[r] <- ref_or_pool$cells$att[r] <- sign * ref_or_pool$mu[link]
    ref_or_pool$inffunc[, r] <- sign * ref_or_pool$phi_mu[, link]
  }
  centered <- scale(ref_or_pool$inffunc, center = TRUE, scale = FALSE)
  ref_or_pool$V_analytical <- crossprod(centered) / nrow(fx_or_pool$raw)
  # Check every cell, full-144-row IF, all cross-covariances, and SEs.
  # This oracle's aggregate maps explicitly assume anticipation = 0.
  cran_compare_mp_to_ref(fit_or_pool, ref_or_pool, fx_or_pool$frame$id, check_agg = FALSE)

  # Reuse the fitted OR object for every exported reporting method, as in
  # test-chained-or-public.R:534-544; plotting constructs an object only.
  expect_match(paste(capture.output(print(fit_or_pool)), collapse = "\n"), "outcome regression/control-mean", fixed = TRUE)
  expect_match(paste(capture.output(summary(fit_or_pool)), collapse = "\n"), "outcome regression/control-mean", fixed = TRUE)
  expect_true(all(cdid::tidy(fit_or_pool)$estimator == fit_or_pool$chained$estimator))
  expect_identical(cdid::glance(fit_or_pool)$estimator, fit_or_pool$chained$estimator)
  expect_identical(cdid::ggdid(fit_or_pool)$labels$caption, fit_or_pool$chained$estimator)

  # Pooled OR under conditional MTAR (pair_logit, across_cohorts, universe = "B") with aggregations
  fx_or_sel <- cran_cp06_fixture(
    route = "B", selection = "conditional_mtar",
    selection_model = "pair_logit", selection_pooling = "across_cohorts",
    est_method = "reg", outcome_pooling = "pooled"
  )
  fit_or_sel <- do.call(cdid::att_gt, fx_or_sel$args)
  raw_or_b <- fx_or_sel$raw; raw_or_b$rho <- fx_or_sel$frame$rho
  ref_or_sel <- cran_ref_env$sp_ref_public_reference(
    raw_or_b, est_method = "reg", selection_model = "pair_logit",
    selection_pooling = "across_cohorts", route = "B", control_group = "nevertreated"
  )
  cran_compare_mp_to_ref(fit_or_sel, ref_or_sel, fx_or_sel$frame$id)

  # Separate OR (default outcome_pooling = "separate") under baseline MTAR and conditional MTAR
  fx_or_sep <- cran_cp06_fixture(route = "none", est_method = "reg", outcome_pooling = "separate")
  args_or_sep <- utils::modifyList(fx_or_sep$args, list(control_group = "notyettreated"))
  fit_or_sep <- do.call(cdid::att_gt, args_or_sep)
  ref_or_sep <- cran_separate_reference(fx_or_sep$raw, est_method = "reg", route = "none", control_group = "notyettreated")
  cran_compare_mp_to_ref(fit_or_sep, ref_or_sep, fx_or_sep$frame$id)

  # Verify that separate and pooled OR differ when multiple comparison blocks have distinct slopes
  fit_or_pool_nyt <- do.call(cdid::att_gt, utils::modifyList(args_or_sep, list(outcome_pooling = "pooled")))
  expect_gt(max(abs(fit_or_sep$att - fit_or_pool_nyt$att), na.rm = TRUE), 1e-4)

  # Separate OR under conditional MTAR (independent_period_logit, within_cohort, universe = "F")
  fx_or_sep_sel <- cran_cp06_fixture(
    route = "F", selection = "conditional_mtar",
    selection_model = "independent_period_logit", selection_pooling = "within_cohort",
    est_method = "reg", outcome_pooling = "separate"
  )
  fit_or_sep_sel <- do.call(cdid::att_gt, fx_or_sep_sel$args)
  ref_or_sep_sel <- cran_separate_reference(
    fx_or_sep_sel$raw, est_method = "reg", route = "F",
    selection_model = "independent_period_logit", selection_pooling = "within_cohort"
  )
  cran_compare_mp_to_ref(fit_or_sep_sel, ref_or_sep_sel, fx_or_sep_sel$frame$id)
})

test_that("pooled and separate doubly robust (DR) estimators match independent references across baseline and conditional MTAR", {
  # Pooled DR under baseline MTAR (NYT, PT-all, universal base)
  fx_dr_base <- cran_cp06_fixture(route = "none", est_method = "dr", outcome_pooling = "pooled")
  args_dr_base <- utils::modifyList(fx_dr_base$args, list(
    control_group = "notyettreated", pt_assumption = "all", pt_start = 2, base_period = "universal"
  ))
  fit_dr_base <- do.call(cdid::att_gt, args_dr_base)
  ref_dr_base <- cran_ref_env$cp07_public_dr_reference(
    fx_dr_base$raw, route = "none", control_group = "notyettreated",
    pt_assumption = "all", pt_start = 2, base_period = "universal"
  )
  cran_compare_mp_to_ref(fit_dr_base, ref_dr_base, fx_dr_base$frame$id)

  # Pooled DR under conditional MTAR (pair_logit, within_cohort, universe = "B")
  fx_dr_b <- cran_cp06_fixture(
    route = "B", selection = "conditional_mtar",
    selection_model = "pair_logit", selection_pooling = "within_cohort",
    est_method = "dr", outcome_pooling = "pooled"
  )
  fit_dr_b <- do.call(cdid::att_gt, fx_dr_b$args)
  raw_dr_b <- fx_dr_b$raw; raw_dr_b$rho <- fx_dr_b$frame$rho
  ref_dr_b <- cran_ref_env$cp07_public_dr_reference(
    raw_dr_b, route = "B", selection_model = "pair_logit",
    selection_pooling = "within_cohort", control_group = "nevertreated"
  )
  cran_compare_mp_to_ref(fit_dr_b, ref_dr_b, fx_dr_b$frame$id)

  # Separate DR (default outcome_pooling = "separate") under conditional MTAR
  fx_dr_sep <- cran_cp06_fixture(
    route = "F", selection = "conditional_mtar",
    selection_model = "independent_period_logit", selection_pooling = "across_cohorts",
    est_method = "dr", outcome_pooling = "separate"
  )
  # Exercise both public defaults against the explicit independent DR target.
  fx_dr_sep$args$est_method <- NULL
  fx_dr_sep$args$outcome_pooling <- NULL
  fit_dr_sep <- do.call(cdid::att_gt, fx_dr_sep$args)
  expect_identical(fit_dr_sep$DIDparams$est_method, "dr")
  expect_identical(fit_dr_sep$DIDparams$outcome_pooling, "separate")
  ref_dr_sep <- cran_separate_reference(
    fx_dr_sep$raw, est_method = "dr", route = "F",
    selection_model = "independent_period_logit", selection_pooling = "across_cohorts"
  )
  cran_compare_mp_to_ref(fit_dr_sep, ref_dr_sep, fx_dr_sep$frame$id)

  # Separate DR with chunk_rows override and intercept-only DR == IPW identity
  fit_dr_chunk <- do.call(cdid::att_gt, utils::modifyList(fx_dr_sep$args, list(control = list(chunk_rows = 5L))))
  expect_equal(fit_dr_sep$att, fit_dr_chunk$att, tolerance = 1e-12)
  expect_equal(fit_dr_sep$inffunc, fit_dr_chunk$inffunc, tolerance = 1e-12)

  fit_dr_int <- do.call(cdid::att_gt, utils::modifyList(fx_dr_base$args, list(
    outcome_pooling = "separate", xformla = ~1
  )))
  fit_ipw_int <- do.call(cdid::att_gt, utils::modifyList(fx_dr_base$args, list(
    est_method = "ipw", xformla = ~1
  )))
  expect_equal(fit_dr_int$att, fit_ipw_int$att, tolerance = 1e-12)
  expect_equal(fit_dr_int$inffunc, fit_ipw_int$inffunc, tolerance = 1e-12)
})

# ------------------------------------------------------------------------------
# 5. Pure change imputation and frozen augmented DR oracles (aug_public_v1)
# ------------------------------------------------------------------------------

test_that("change imputation and augmented DR match independent analytic and frozen aug_public_v1 oracles", {
  aug_dir <- testthat::test_path("fixtures", "aug_public_v1")
  aug_raw <- cran_ref_env$imputation_public_raw(utils::read.csv(file.path(aug_dir, "raw.csv")))
  # Frozen aug_public_v1 oracles use this exact selection-only covariate
  # (test-chained-augmented-public.R:27); changing V changes the fitted model.
  aug_raw$V <- ((5 * aug_raw$id) %% 11 - 5) / 5

  # 5A. Pure change imputation (selection_est_method = "imputation") across ipw, reg, dr and poolings
  fx_imp <- cran_ref_env$imputation_public_fixture(aug_raw, dates = c(1, 2, 3, 5))
  args_imp_base <- utils::modifyList(fx_imp$args, list(
    imputation_pooling = "within_cohort",
    control_group = "notyettreated",
    pt_assumption = "all", pt_start = 2, base_period = "universal",
    clustervars = "cl"
  ))
  bp_imp <- cran_ref_env$imputation_public_blueprint(
    aug_raw, imputation_pooling = "within_cohort",
    dates = c(1, 2, 3, 5), control_group = "notyettreated",
    pt_assumption = "all", pt_start = 2, base_period = "universal"
  )
  ref_imp_all <- cran_ref_env$imputation_public_reference(
    aug_raw, est_method = c("ipw", "reg", "dr"),
    imputation_pooling = "within_cohort",
    dates = c(1, 2, 3, 5), control_group = "notyettreated",
    pt_assumption = "all", pt_start = 2, base_period = "universal",
    cluster = fx_imp$frame$cl, blueprint = bp_imp
  )
  for (em in c("ipw", "reg", "dr")) {
    fit_imp <- do.call(cdid::att_gt, utils::modifyList(args_imp_base, list(est_method = em)))
    ref_em <- ref_imp_all$results[[em]]
    cran_compare_mp_to_ref(fit_imp, ref_em, fx_imp$frame$id, tol_att = 1e-11, tol_if = 1e-9, tol_v = 1e-9, tol_se = 1e-9)
  }

  # Also exercise imputation_pooling = "across_cohorts" with default outcome_pooling = "separate"
  fx_imp_sep <- cran_ref_env$imputation_public_fixture(aug_raw, dates = 1:4)
  fit_imp_sep <- do.call(cdid::att_gt, utils::modifyList(fx_imp_sep$args, list(
    outcome_pooling = "separate", imputation_pooling = "across_cohorts", est_method = "dr"
  )))
  expect_s3_class(fit_imp_sep, "cdid_MP")
  expect_identical(fit_imp_sep$DIDparams$outcome_pooling, "separate")
  expect_identical(nrow(fit_imp_sep$inffunc), 144L)

  # 5B. Frozen independent augmented DR oracles (selection_est_method = "dr")
  aug_configs <- list(
    pair_sep_jsep_NT_NULL = list(model = "pair_logit", qpool = "separate", jpool = "separate", cg = "nevertreated", q = "NULL"),
    pair_within_jacross_NYT_xV = list(model = "pair_logit", qpool = "within_cohort", jpool = "across_cohorts", cg = "notyettreated", q = "xV"),
    date_across_jacross_NYT_NULL = list(model = "independent_period_logit", qpool = "across_cohorts", jpool = "across_cohorts", cg = "notyettreated", q = "NULL"),
    pair_irregular_all_universal_cluster_xV = list(model = "pair_logit", qpool = "across_cohorts", jpool = "within_cohort", cg = "notyettreated", q = "xV", irregular = TRUE)
  )
  for (cfg_name in names(aug_configs)) {
    cfg <- aug_configs[[cfg_name]]
    dates <- if (isTRUE(cfg$irregular)) c(1, 2, 3, 5) else 1:4
    fx <- cran_ref_env$imputation_public_fixture(aug_raw, dates = dates)
    fx$frame$V <- aug_raw$V
    args <- fx$args
    args$outcome_pooling <- "pooled"
    args$frame <- fx$frame
    args$selection <- "conditional_mtar"
    args$selection_est_method <- "dr"
    args$est_method <- "dr"
    args$imputation_pooling <- cfg$jpool
    args$selection_model <- cfg$model
    args$selection_pooling <- cfg$qpool
    args$control_group <- cfg$cg
    sel_form <- switch(cfg$q, `NULL` = NULL, xV = ~x + V, V = ~V, intercept_only = ~1)
    if (!is.null(sel_form)) args$selection_xformla <- sel_form
    if (isTRUE(cfg$irregular)) {
      args$pt_assumption <- "all"
      args$pt_start <- 2
      args$base_period <- "universal"
      args$clustervars <- "cl"
    }
    fit_aug <- do.call(cdid::att_gt, args)
    ref_aug <- readRDS(file.path(aug_dir, "expected", paste0(cfg_name, ".rds")))$results$dr
    cran_compare_mp_to_ref(fit_aug, ref_aug, fx$frame$id, tol_att = 1e-11, tol_if = 3e-7, tol_v = 1e-7, tol_se = 1e-8)
  }

  # Augmented observation route with est_method = "ipw", "reg", and separate outcome pooling
  fx_aug4 <- cran_ref_env$imputation_public_fixture(aug_raw, dates = 1:4)
  fx_aug4$frame$V <- aug_raw$V
  for (em in c("ipw", "reg", "dr")) {
    fit_em <- do.call(cdid::att_gt, utils::modifyList(fx_aug4$args, list(
      frame = fx_aug4$frame, selection = "conditional_mtar", selection_est_method = "dr",
      est_method = em, outcome_pooling = "separate", selection_xformla = ~x + V
    )))
    expect_s3_class(fit_em, "cdid_MP")
    expect_identical(nrow(fit_em$inffunc), 144L)
    expect_true(all(is.finite(fit_em$att[fit_em$chained$available])))
  }
})

# ------------------------------------------------------------------------------
# 6. Time-varying covariates and changing sampling weights (fix_weights)
# ------------------------------------------------------------------------------

test_that("time-varying covariates and fix_weights anchors match independent d60_routes_reference", {
  d_cpt <- utils::read.csv(testthat::test_path("fixtures", "cpt", "fBd.csv"))
  # Construct positive time-varying sampling weights
  d_cpt$w <- d_cpt$w * (1 + 0.08 * match(d_cpt$time, sort(unique(d_cpt$time))) + 0.03 * (d_cpt$id %% 3L))
  ids_cpt <- as.character(sort(unique(d_cpt$id)))
  N_cpt <- length(ids_cpt)

  for (fw in c("first_period", "base_period")) {
    rt <- list(
      kind = "corollary1", est = "dr",
      sel_model = "independent_period_logit", sel_pooling = "within_cohort", universe = "F"
    )
    S <- cran_ref_env$d60r_setup(
      d_cpt, route = rt, formula = ~x + z, sformula = ~sx + sz,
      fix_weights = fw, control_group = "notyettreated"
    )
    ref_d60 <- cran_ref_env$d60r_if(S, eps = 1e-4)
    mp_d60 <- cdid::att_gt(
      yname = "y", tname = "time", idname = "id", gname = "G",
      xformla = ~x + z, selection = "conditional_mtar",
      selection_model = "independent_period_logit", selection_pooling = "within_cohort",
      selection_xformla = ~sx + sz, weightsname = "w", fix_weights = fw,
      est_method = "dr", outcome_pooling = "pooled",
      control_group = "notyettreated", gmm_weighting = "identity",
      data = d_cpt, bstrap = FALSE, cband = FALSE
    )
    post_idx <- which(mp_d60$chained$role == "post")
    pk <- paste0(mp_d60$group[post_idx], ":", mp_d60$t[post_idx])
    rk <- paste0(ref_d60$base$cells$g, ":", ref_d60$base$cells$t)
    m_idx <- match(rk, pk)
    expect_false(anyNA(m_idx))
    expect_lt(max(abs(mp_d60$att[post_idx[m_idx]] - ref_d60$base$cells$att)), 1e-10)
    ri <- match(as.character(S$ids), rownames(mp_d60$inffunc))
    # Retain the original D60 criterion (test-chained-reference-inputs.R:48-51).
    expect_lte(max(abs(mp_d60$inffunc[ri, post_idx[m_idx], drop = FALSE] - ref_d60$phi_cells)),
               3e-7 * max(1, max(abs(ref_d60$phi_cells))))
    V_ref <- cran_ref_env$imputation_public_covariance(ref_d60$phi_cells)
    expect_equal(unname(mp_d60$V_analytical[post_idx[m_idx], post_idx[m_idx], drop = FALSE]),
                 unname(V_ref), tolerance = 1e-7)
    expect_equal(unname(mp_d60$se[post_idx[m_idx]]), sqrt(diag(V_ref) / N_cpt), tolerance = 1e-8)

    for (tp in c("simple", "group", "dynamic", "calendar")) {
      ag <- cdid::aggte(mp_d60, type = tp, bstrap = FALSE, cband = FALSE)
      expect_lt(abs(ag$overall.att - unname(ref_d60$base$agg[tp])), 1e-10)
      ifv <- switch(tp,
        simple = ag$inf.function$simple.att,
        group = ag$inf.function$selective.inf.func,
        dynamic = ag$inf.function$dynamic.inf.func,
        calendar = ag$inf.function$calendar.inf.func
      )
      phi_ref <- ref_d60$phi_agg[, match(tp, c("simple", "group", "dynamic", "calendar"))]
      expect_lte(max(abs(as.matrix(ifv)[ri, 1L] - phi_ref)), 3e-7 * max(1, max(abs(phi_ref))))
      expect_equal(ag$overall.se, sqrt(sum((phi_ref - mean(phi_ref))^2)) / N_cpt, tolerance = 1e-8)
    }
  }

  # Refusal when weights change within a unit and fix_weights is not specified
  expect_error(
    cdid::att_gt(
      yname = "y", tname = "time", idname = "id", gname = "G",
      xformla = ~x + z, weightsname = "w", data = d_cpt, bstrap = FALSE
    ),
    class = "cdid_chained_request_unsupported_weights"
  )
})

# ------------------------------------------------------------------------------
# 7. Multiplier bootstrap, simultaneous confidence bands, pretrend Wald test,
#    conditional_did_pretest() CvM process, custom IPW callback, and refusals
# ------------------------------------------------------------------------------

test_that("multiplier bootstrap, cband, RNG discipline, pretrend Wald test, conditional_did_pretest, and custom callbacks hold", {
  fx <- cran_cp06_fixture(route = "none", est_method = "dr", outcome_pooling = "separate")

  # bstrap = FALSE does not consume caller RNG
  set.seed(20261106L)
  seed_before <- .Random.seed
  mp_analytic <- do.call(cdid::att_gt, utils::modifyList(fx$args, list(
    control_group = "notyettreated", bstrap = FALSE, cband = TRUE
  )))
  expect_identical(.Random.seed, seed_before)

  # Unconditional pretrend Wald test W and Wpval match exact quadratic form on pretrend cells
  expect_identical(mp_analytic$chained$pretest$status, "computed")
  pre_cells <- which(mp_analytic$chained$role == "pre" & mp_analytic$chained$available)
  expect_gt(length(pre_cells), 0L)
  N <- mp_analytic$n
  theta_pre <- mp_analytic$att[pre_cells]
  V_pre <- mp_analytic$V_analytical[pre_cells, pre_cells, drop = FALSE]
  W_expected <- as.double(N * crossprod(theta_pre, solve(V_pre, theta_pre)))
  expect_equal(mp_analytic$W, W_expected, tolerance = 1e-12)
  expect_equal(mp_analytic$Wpval, stats::pchisq(W_expected, df = length(pre_cells), lower.tail = FALSE), tolerance = 1e-12)

  # Multiplier bootstrap + simultaneous confidence bands with clustervars and seed reproducibility
  boot_args <- utils::modifyList(fx$args, list(
    control_group = "notyettreated", clustervars = "cl",
    bstrap = TRUE, cband = TRUE, biters = 49L
  ))
  set.seed(901L)
  mp_boot1 <- do.call(cdid::att_gt, boot_args)
  set.seed(901L)
  mp_boot2 <- do.call(cdid::att_gt, boot_args)
  expect_identical(mp_boot1$se, mp_boot2$se)
  expect_identical(mp_boot1$c, mp_boot2$c)
  expect_gte(mp_boot1$c, stats::qnorm(0.975))

  set.seed(902L)
  ag_dyn1 <- cdid::aggte(mp_boot1, type = "dynamic", bstrap = TRUE, cband = TRUE, biters = 49L)
  set.seed(902L)
  ag_dyn2 <- cdid::aggte(mp_boot1, type = "dynamic", bstrap = TRUE, cband = TRUE, biters = 49L)
  expect_identical(ag_dyn1$overall.se, ag_dyn2$overall.se)
  expect_identical(ag_dyn1$se.egt, ag_dyn2$se.egt)
  expect_identical(ag_dyn1$crit.val.egt, ag_dyn2$crit.val.egt)

  # Genuine cell unavailability and pretest unavailability when a pre/post cell is unobserved
  fx_unavail <- fx
  fx_unavail$args$control_group <- "notyettreated"
  fx_unavail$args$data$y[fx_unavail$args$data$g == 4 & fx_unavail$args$data$time == 2] <- NA_real_
  mp_unavail <- do.call(cdid::att_gt, fx_unavail$args)
  expect_true(any(!mp_unavail$chained$available))
  expect_true(all(is.na(mp_unavail$att[!mp_unavail$chained$available])))
  expect_true(all(is.na(mp_unavail$se[!mp_unavail$chained$available])))
  expect_identical(mp_unavail$chained$pretest$status, "unavailable")
  expect_null(mp_unavail$W)
  expect_null(mp_unavail$Wpval)

  # conditional_did_pretest() against frozen independent CPT fixtures (baseline and imputation)
  cpt_dir <- testthat::test_path("fixtures", "cpt")
  cpt_full <- readRDS(file.path(cpt_dir, "full.rds"))
  cpt_imp <- readRDS(file.path(cpt_dir, "imp.rds"))
  cpt_priv <- readRDS(file.path(cpt_dir, "priv.rds"))
  d_fB <- utils::read.csv(file.path(cpt_dir, "fB.csv"))
  signs_B <- do.call(rbind, lapply(cpt_priv$fixtures$fixtures$B$unit_rademacher_signs, function(r) as.double(unlist(r))))
  chained_cpt_fn <- getFromNamespace("chained_conditional_did_pretest", "cdid")

  cpt_b4 <- chained_cpt_fn(
    outcome_pooling = "pooled", yname = "y", tname = "time", idname = "id", gname = "G",
    xformla = ~x + z, data = d_fB, panel = TRUE, allow_unbalanced_panel = NULL,
    control_group = "notyettreated", weightsname = "w", alp = 0.05,
    bstrap = TRUE, cband = TRUE, biters = nrow(signs_B), clustervars = NULL,
    est_method = "dr", print_details = FALSE, pl = FALSE, cores = 1,
    anticipation = 0, base_period = "varying", fix_weights = NULL, mode = NULL, frame = NULL,
    selection = "conditional_mtar", selection_model = "pair_logit",
    selection_pooling = "across_cohorts", selection_xformla = ~sx + sz,
    selection_universe = "F", selection_rho = NULL, selection_justification = NULL,
    selection_est_method = "ipw", imputation_pooling = "separate",
    control = NULL, fixed_signs = signs_B
  )
  gold_b4 <- cpt_full$reference$cases$B4$base$cvm[["varying|dr"]]$unit
  expect_s3_class(cpt_b4, c("cdid_MP.TEST", "MP.TEST"))
  expect_equal(cpt_b4$CvM, gold_b4$CvM, tolerance = 1e-10)
  expect_lt(max(abs(cpt_b4$CvMb - unlist(gold_b4$CvMb))), 1e-9)
  expect_identical(cpt_b4$CvMpval, gold_b4$CvMpval)

  cpt_b2 <- chained_cpt_fn(
    outcome_pooling = "pooled", yname = "y", tname = "time", idname = "id", gname = "G",
    xformla = ~x + z, data = d_fB, panel = TRUE, allow_unbalanced_panel = NULL,
    control_group = "notyettreated", weightsname = "w", alp = 0.05,
    bstrap = TRUE, cband = TRUE, biters = nrow(signs_B), clustervars = NULL,
    est_method = "dr", print_details = FALSE, pl = FALSE, cores = 1,
    anticipation = 0, base_period = "varying", fix_weights = NULL, mode = NULL, frame = NULL,
    selection = "conditional_mtar", selection_model = "pair_logit",
    selection_pooling = "separate", selection_xformla = NULL,
    selection_universe = "F", selection_rho = NULL, selection_justification = NULL,
    selection_est_method = "imputation", imputation_pooling = "within_cohort",
    control = NULL, fixed_signs = signs_B
  )
  gold_b2 <- cpt_imp$reference$cases$B2$base$cvm[["varying|imp_dr"]]$unit
  expect_equal(cpt_b2$CvM, gold_b2$CvM, tolerance = 1e-10)
  expect_lt(max(abs(cpt_b2$CvMb - unlist(gold_b2$CvMb))), 1e-9)
  expect_identical(cpt_b2$CvMpval, gold_b2$CvMpval)

  # Public conditional_did_pretest() call with caller RNG reproducibility
  set.seed(77L)
  pub_cpt1 <- cdid::conditional_did_pretest(
    yname = "y", tname = "time", idname = "id", gname = "G",
    xformla = ~x + z, data = d_fB, control_group = "notyettreated",
    weightsname = "w", est_method = "dr", biters = 19L
  )
  set.seed(77L)
  pub_cpt2 <- cdid::conditional_did_pretest(
    yname = "y", tname = "time", idname = "id", gname = "G",
    xformla = ~x + z, data = d_fB, control_group = "notyettreated",
    weightsname = "w", est_method = "dr", biters = 19L
  )
  expect_identical(pub_cpt1$chained$status, "available")
  expect_identical(pub_cpt1$CvM, pub_cpt2$CvM)
  expect_identical(pub_cpt1$CvMb, pub_cpt2$CvMb)

  # Preserve the accepted ordinary conditional-pretest delegation witness
  # (test-chained-conditional-pretest.R:345-359), including lazy chained args.
  d_fA <- utils::read.csv(file.path(cpt_dir, "fA.csv"))
  d_fA <- d_fA[d_fA$G != 1, ]
  ordinary_cpt <- function(f, ...) {
    set.seed(3L)
    out <- NULL
    invisible(capture.output(out <- suppressMessages(f(
      yname = "y", tname = "time", idname = "id", gname = "G",
      xformla = ~x, data = d_fA, biters = 29L, ...
    ))))
    out
  }
  ord_cpt <- ordinary_cpt(cdid::conditional_did_pretest, mode = "ordinary",
    selection = stop("chained argument evaluated"), control = stop("chained argument evaluated"),
    base_period = stop("chained argument evaluated"))
  ref_cpt <- ordinary_cpt(did::conditional_did_pretest)
  expect_identical(class(ord_cpt), class(ref_cpt))
  expect_equal(ord_cpt$CvM, ref_cpt$CvM)
  expect_equal(ord_cpt$CvMb, ref_cpt$CvMb)
  expect_equal(ord_cpt$CvMpval, ref_cpt$CvMpval)

  # Custom IPW callback with and without explicit log_weight_gradient
  custom_hajek_cb <- function(give_grad = TRUE) {
    function(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...) {
      X <- covariates; dY <- y1 - y0; iw <- i.weights; m <- length(D)
      fit <- suppressWarnings(stats::glm.fit(
        x = X, y = D, weights = iw, family = stats::binomial(),
        control = stats::glm.control(epsilon = 1e-13, maxit = 60L)
      ))
      eta <- drop(X %*% fit$coefficients); p <- stats::plogis(eta); o <- exp(eta)
      wt <- iw * D; wc <- iw * (1 - D) * o; S1 <- sum(wt); S0 <- sum(wc)
      Rp <- sum(wt * dY) / S1; Rm <- sum(wc * dY) / S0
      direct <- wt * (dY - Rp) / S1 - wc * (dY - Rm) / S0
      lever <- drop(solve(crossprod(X, (iw * p * (1 - p)) * X), crossprod(X, wc * (dY - Rm)) / S0))
      sdir <- iw * (D - p) * drop(X %*% lever)
      out <- list(ATT = Rp - Rm, att.inf.func = m * (direct - sdir))
      if (give_grad) out$log_weight_gradient <- direct - sdir
      out
    }
  }
  fx_cb <- cran_cp06_fixture(
    route = "F", selection = "conditional_mtar",
    selection_model = "pair_logit", selection_pooling = "across_cohorts",
    est_method = "ipw"
  )
  mp_builtin_ipw <- do.call(cdid::att_gt, fx_cb$args)
  mp_custom_ipw <- do.call(cdid::att_gt, utils::modifyList(fx_cb$args, list(est_method = custom_hajek_cb(TRUE))))
  expect_equal(mp_custom_ipw$att, mp_builtin_ipw$att, tolerance = 1e-11)
  expect_equal(mp_custom_ipw$inffunc, mp_builtin_ipw$inffunc, tolerance = 1e-9)

  # Under conditional_mtar, omitting log_weight_gradient in a custom callback is refused
  expect_error(
    do.call(cdid::att_gt, utils::modifyList(fx_cb$args, list(est_method = custom_hajek_cb(FALSE)))),
    class = "cdid_custom_ipw_derivative_error"
  )
})
