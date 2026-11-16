# ------------------------------------------------------------------------------
# Comprehensive tests for cross-pair chained DiD, PT-ctrl, within-unit GMM
# weighting, spanning-forest Gram preservation, and chain_report():
# 1. Exact machine-precision equivalence to stacked target-slice unit-FE OLS
#    (and cohort-FE slice OLS on cohort-balanced windows)
# 2. Strict cross-pair identification when target and controls share zero pairs,
#    plus empty-catalog / zero-estimable-target safety
# 3. Cross-pair IPW / OR / DR covariate adjustment, exact M-estimation IFs,
#    and WLS alias / logit separation checks
# 4. Spanning-forest Gram matrix G = E'E invariance (including non-clique masks)
#    and analytical weight-overlap within-unit covariance
# 5. Selective-attrition robustness vs cohort-FE slice regression
# 6. chain_report() S3 diagnostics (including include_pre = TRUE, unreduced
#    direct-link detection, and empty catalogs) and point-only / frame contracts
# ------------------------------------------------------------------------------

fit_reference_slice_ols <- function(df, fe = c("unit", "cohort"),
                                    weightsname = NULL, anticipation = 0L,
                                    control_group = c("notyettreated", "nevertreated")) {
  fe <- match.arg(fe)
  control_group <- match.arg(control_group)
  all_g <- sort(unique(df$G[df$G > 0]))
  all_t <- sort(unique(df$period))
  out <- list()

  for (g in all_g) {
    adm_pre <- all_t[all_t < g - anticipation]
    if (!length(adm_pre)) next
    ref_t <- max(adm_pre)
    ctrl_cond <- if (control_group == "notyettreated") {
      (df$G == 0) | (df$G > g & df$period + anticipation < df$G)
    } else {
      (df$G == 0)
    }
    sub <- df[(df$G == g) | ctrl_cond, , drop = FALSE]
    if (!any(sub$G != g)) next
    sub$D_slice <- as.integer(sub$G == g)
    t_sub <- sort(unique(sub$period))
    if (!(ref_t %in% sub$period[sub$D_slice == 1L])) next

    sub$period_f <- stats::relevel(factor(sub$period, levels = t_sub), ref = as.character(ref_t))
    sub$fe_f <- if (fe == "unit") factor(sub$id) else factor(sub$G)
    w <- if (is.null(weightsname)) rep(1.0, nrow(sub)) else as.double(sub[[weightsname]])

    t_g <- setdiff(sort(unique(sub$period[sub$D_slice == 1L])), ref_t)
    X_inter <- vapply(t_g, function(tv) as.numeric(sub$D_slice == 1L & sub$period == tv),
                      numeric(nrow(sub)))
    if (!is.matrix(X_inter)) X_inter <- matrix(X_inter, nrow = nrow(sub))
    colnames(X_inter) <- paste0("att_", g, "_", t_g)
    fit <- stats::lm(sub$y ~ 0 + sub$fe_f + sub$period_f + X_inter, weights = w)
    cf <- stats::coef(fit)
    est <- unname(cf[paste0("X_interatt_", g, "_", t_g)])
    out[[length(out) + 1L]] <- data.frame(
      group = rep(g, length(t_g) + 1L),
      t = c(ref_t, t_g),
      att_ols = c(0.0, est),
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, out)
  res[order(res$group, res$t), , drop = FALSE]
}

make_short_window_panel <- function(seed = 20261111L,
                                    n_per_cohort = c(30L, 45L, 35L, 50L, 40L, 38L)) {
  set.seed(seed)
  cohorts <- 2005:2010
  rows <- list()
  uid <- 0L
  for (ci in seq_along(cohorts)) {
    g <- cohorts[ci]
    n_c <- n_per_cohort[ci]
    win_t <- seq.int(max(2001L, g - 4L), 2010L)
    for (i in seq_len(n_c)) {
      uid <- uid + 1L
      alpha_i <- stats::rnorm(1L, mean = 0.7 * (g - 2005), sd = 1.1)
      w_i <- stats::runif(1L, 0.4, 2.4)
      cl_i <- (uid - 1L) %/% 3L + 1L
      for (tv in win_t) {
        tau <- if (tv >= g) 2.0 + 0.5 * (tv - g) + 0.15 * (g - 2005) else 0.0
        y_it <- alpha_i + 0.35 * (tv - 2001) + tau + stats::rnorm(1L, sd = 0.45)
        rows[[length(rows) + 1L]] <- data.frame(
          id = uid, cluster = cl_i, period = tv, G = g, y = y_it, w = w_i
        )
      }
    }
  }
  do.call(rbind, rows)
}

test_that("reduced same-pair fits preserve unsupported targets in short windows", {
  fit <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    data = make_short_window_panel(),
    control_group = "notyettreated", est_method = "ipw",
    cross_pair = FALSE, bstrap = FALSE, cband = FALSE
  )
  last <- fit$t == 2010 & fit$group < 2010
  expect_equal(sum(last), 5L)
  expect_true(all(is.na(fit$att[last])))
  expect_true(all(fit$chained$support_status[last] == "isolated_unsupported"))
  expect_true(any(is.finite(fit$att[fit$t >= fit$group & fit$t < 2010])))
})

test_that("cross-pair PT-ctrl within-unit GMM reproduces stacked target-slice OLS to machine precision across balanced, weighted, gap, and monotone-attrition panels", {
  df <- make_short_window_panel()

  # 1. Balanced cohort windows: matches BOTH unit-FE and cohort-FE slice OLS across all 40 cells
  mp <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df, bstrap = FALSE, cband = FALSE
  )
  ols_u <- fit_reference_slice_ols(df, fe = "unit")
  ols_c <- fit_reference_slice_ols(df, fe = "cohort")
  cmp <- merge(
    merge(data.frame(group = mp$group, t = mp$t, att_cdid = mp$att), ols_u, by = c("group", "t")),
    ols_c, by = c("group", "t"), suffixes = c("_unit", "_cohort")
  )
  expect_equal(nrow(cmp), 40L)
  expect_identical(is.na(cmp$att_cdid), is.na(cmp$att_ols_unit))
  expect_identical(is.na(cmp$att_cdid), is.na(cmp$att_ols_cohort))
  id_mask <- !is.na(cmp$att_cdid)
  expect_equal(sum(id_mask), 30L)
  expect_equal(sum(id_mask & cmp$t >= cmp$group), 15L)
  expect_lt(max(abs(cmp$att_cdid[id_mask] - cmp$att_ols_unit[id_mask])), 1e-12)
  expect_lt(max(abs(cmp$att_cdid[id_mask] - cmp$att_ols_cohort[id_mask])), 1e-12)

  # 2. Dynamic event-study cohort-mass weighted aggregation
  agg_dyn <- cdid::aggte(mp, type = "dynamic", na.rm = TRUE, bstrap = FALSE)
  n_by_g <- table(df$G[!duplicated(df$id)])
  ols_post_pre <- cmp
  ols_post_pre$e <- ols_post_pre$t - ols_post_pre$group
  ols_post_pre$n_g <- as.double(n_by_g[as.character(ols_post_pre$group)])
  e_grid <- agg_dyn$egt
  dyn_manual <- vapply(e_grid, function(ev) {
    sub_e <- ols_post_pre[ols_post_pre$e == ev & is.finite(ols_post_pre$att_ols_unit), , drop = FALSE]
    if (!nrow(sub_e)) return(NA_real_)
    sum(sub_e$n_g * sub_e$att_ols_unit) / sum(sub_e$n_g)
  }, numeric(1L))
  expect_identical(is.na(agg_dyn$att.egt), is.na(dyn_manual))
  ok_e <- !is.na(dyn_manual)
  expect_lt(max(abs(agg_dyn$att.egt[ok_e] - dyn_manual[ok_e])), 1e-12)

  # 3. Survey-weighted WLS equivalence
  mp_w <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    weightsname = "w", control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df, bstrap = FALSE, cband = FALSE
  )
  ols_uw <- fit_reference_slice_ols(df, fe = "unit", weightsname = "w")
  cmp_w <- merge(data.frame(group = mp_w$group, t = mp_w$t, att_cdid = mp_w$att), ols_uw, by = c("group", "t"))
  expect_identical(is.na(cmp_w$att_cdid), is.na(cmp_w$att_ols))
  expect_lt(max(abs(cmp_w$att_cdid[!is.na(cmp_w$att_cdid)] - cmp_w$att_ols[!is.na(cmp_w$att_ols)])), 1e-12)

  # 4. Monotone unit attrition inside cohorts: matches unit-FE slice OLS to machine precision
  df_unbal <- df[!(df$id %% 4L == 0L & df$period >= 2007L), , drop = FALSE]
  mp_u <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df_unbal, bstrap = FALSE, cband = FALSE
  )
  ols_unbal <- fit_reference_slice_ols(df_unbal, fe = "unit")
  cmp_u <- merge(data.frame(group = mp_u$group, t = mp_u$t, att_cdid = mp_u$att), ols_unbal, by = c("group", "t"))
  expect_identical(is.na(cmp_u$att_cdid), is.na(cmp_u$att_ols))
  expect_lt(max(abs(cmp_u$att_cdid[!is.na(cmp_u$att_cdid)] - cmp_u$att_ols[!is.na(cmp_u$att_ols)])), 1e-12)

  # 5. Biennial wave gaps (2007 and 2009 absent) & anticipation = 1
  df_gap <- df[!(df$period %in% c(2007L, 2009L)), , drop = FALSE]
  mp_gap <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df_gap, bstrap = FALSE, cband = FALSE
  )
  ols_gap <- fit_reference_slice_ols(df_gap, fe = "unit")
  cmp_gap <- merge(data.frame(group = mp_gap$group, t = mp_gap$t, att_cdid = mp_gap$att),
                   ols_gap[is.finite(ols_gap$att_ols), ], by = c("group", "t"))
  expect_equal(nrow(cmp_gap), 20L)
  expect_lt(max(abs(cmp_gap$att_cdid - cmp_gap$att_ols)), 1e-12)

  mp_ant <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    anticipation = 1, control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df, bstrap = FALSE, cband = FALSE
  )
  ols_ant <- fit_reference_slice_ols(df, fe = "unit", anticipation = 1L)
  cmp_ant <- merge(data.frame(group = mp_ant$group, t = mp_ant$t, att_cdid = mp_ant$att),
                   ols_ant[is.finite(ols_ant$att_ols), ], by = c("group", "t"))
  expect_lt(max(abs(cmp_ant$att_cdid - cmp_ant$att_ols)), 1e-12)

  # 6. Varying vs universal base_period leaves all post-treatment ATTs and SEs identical
  mp_var <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", base_period = "varying",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df, bstrap = FALSE, cband = FALSE
  )
  post_u <- mp$t >= mp$group
  post_v <- mp_var$t >= mp_var$group
  expect_equal(mp$att[post_u], mp_var$att[post_v], tolerance = 1e-14)
  expect_equal(mp$se[post_u], mp_var$se[post_v], tolerance = 1e-14)
})

test_that("cross_pair = TRUE strictly expands identification when target and controls share zero common date pairs, and handles empty or disconnected catalogs cleanly", {
  set.seed(101)
  rows <- list()
  uid <- 0L
  # Target G=3 observed at {1, 2, 4}; G=4 at {2, 3}; G=5 at {3, 4}; G=0 at {1, 2}
  spec <- list(
    list(g = 3L, times = c(1L, 2L, 4L)),
    list(g = 4L, times = c(2L, 3L)),
    list(g = 5L, times = c(3L, 4L)),
    list(g = 0L, times = c(1L, 2L))
  )
  for (s in spec) {
    for (i in seq_len(120L)) {
      uid <- uid + 1L
      a_i <- 0.8 * s$g + stats::rnorm(1L, sd = 0.2)
      for (tv in s$times) {
        eff <- if (s$g == 3L && tv >= 3L) 3.5 else 0.0
        rows[[length(rows) + 1L]] <- data.frame(
          id = uid, period = tv, G = s$g, y = a_i + 0.6 * tv + eff + stats::rnorm(1L, sd = 0.05)
        )
      }
    }
  }
  df_disjoint <- do.call(rbind, rows)

  mp_same <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", cross_pair = FALSE, data = df_disjoint, bstrap = FALSE
  )
  i_same <- which(mp_same$group == 3 & mp_same$t == 4)
  expect_true(is.na(mp_same$att[i_same]))

  for (pt in c("post", "ctrl", "all")) {
    mp_cp <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      control_group = "notyettreated", pt_assumption = pt,
      pt_start = if (pt == "all") 1 else NULL,
      cross_pair = TRUE, data = df_disjoint, bstrap = FALSE
    )
    i_cp <- which(mp_cp$group == 3 & mp_cp$t == 4)
    expect_true(is.finite(mp_cp$att[i_cp]))
    expect_lt(abs(mp_cp$att[i_cp] - 3.5), 0.05)
  }

  mp_capture <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", pt_assumption = "all", pt_start = 1,
    cross_pair = TRUE, data = df_disjoint, bstrap = FALSE,
    control = list(keep_nuisance = TRUE)
  )
  for (field in c("att", "se", "inffunc")) {
    expect_identical(mp_capture[[field]], mp_cp[[field]])
  }
  report <- cdid::nuisance_models(mp_capture)
  expect_s3_class(report, "cdid_nuisance_index")
  expect_identical(report$components$status, rep("not_applicable", 4L))
  expect_identical(cdid::nuisance_models(mp_cp)$code, "not_requested")

  # Disconnected control graph: drop G=4 so gap 2->3 has zero clean control links
  df_broken <- df_disjoint[df_disjoint$G != 4L, , drop = FALSE]
  mp_broken <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", pt_assumption = "ctrl",
    cross_pair = TRUE, data = df_broken, bstrap = FALSE
  )
  i_br <- which(mp_broken$group == 3 & mp_broken$t == 4)
  expect_true(is.na(mp_broken$att[i_br]))

  # Zero estimable targets / empty-catalog support-only fit (e.g. nevertreated when G=0 observed at single date only)
  df_zero <- rbind(
    df_disjoint[df_disjoint$G == 3L & df_disjoint$period %in% c(2L, 4L), , drop = FALSE],
    df_disjoint[df_disjoint$G == 0L & df_disjoint$period == 1L, , drop = FALSE]
  )
  mp_zero_cp <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "nevertreated", cross_pair = TRUE,
    data = df_zero, bstrap = FALSE, cband = FALSE
  )
  expect_true(all(is.na(mp_zero_cp$att[mp_zero_cp$t >= mp_zero_cp$group])))
  cr_zero_cp <- cdid::chain_report(mp_zero_cp, include_pre = TRUE)
  expect_s3_class(cr_zero_cp, "cdid_chain_report")

  # Empty link catalog in same-pair mode: chain_report() succeeds cleanly without crashing
  df_empty_sp <- rbind(
    data.frame(id = 1:5, period = 1L, G = 3L, y = rnorm(5)),
    data.frame(id = 1:5, period = 2L, G = 3L, y = rnorm(5)),
    data.frame(id = 11:15, period = 3L, G = 0L, y = rnorm(5))
  )
  mp_empty_sp <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "nevertreated", cross_pair = FALSE,
    data = df_empty_sp, bstrap = FALSE, cband = FALSE
  )
  cr_empty <- cdid::chain_report(mp_empty_sp, include_pre = TRUE)
  expect_s3_class(cr_empty, "cdid_chain_report")
  expect_equal(cr_empty$meta$n_moments, 0L)
  expect_equal(nrow(as.data.frame(cr_empty, what = "witness")), 0L)
  expect_equal(nrow(as.data.frame(cr_empty, what = "gmm_weights")), 0L)
})

test_that("cross_pair covariate adjustment (ipw, reg, dr) has exact M-estimation IFs, rank/separation guards, clustering, and bootstrap bands", {
  set.seed(2026)
  n_per <- 35L
  rows <- list()
  uid <- 0L
  spec <- list(
    list(g = 3L, times = c(1L, 2L, 3L, 4L)),
    list(g = 4L, times = c(2L, 3L)),
    list(g = 5L, times = c(3L, 4L)),
    list(g = 0L, times = c(1L, 2L))
  )
  for (s in spec) {
    for (i in seq_len(n_per)) {
      uid <- uid + 1L
      x1 <- stats::rnorm(1L, mean = if (s$g == 3L) 0.35 else -0.15)
      cl <- (uid - 1L) %/% 2L + 1L
      a_i <- 0.5 * s$g + 0.9 * x1 + stats::rnorm(1L, sd = 0.25)
      for (tv in s$times) {
        eff <- if (s$g == 3L && tv >= 3L) 2.2 + 0.4 * (tv - 3L) else 0.0
        y <- a_i + (0.4 + 0.5 * x1) * tv + eff + stats::rnorm(1L, sd = 0.2)
        rows[[length(rows) + 1L]] <- data.frame(
          id = uid, cluster = cl, period = tv, G = s$g, x1 = x1, w = 1.0, y = y
        )
      }
    }
  }
  df_cov <- do.call(rbind, rows)
  unit_ids <- unique(df_cov$id)
  N <- length(unit_ids)

  # 1. Verify analytical unit IF against finite-difference unit-weight perturbation
  eps <- 1e-6
  set.seed(77)
  dir_u <- stats::rnorm(N)
  dir_u <- dir_u - mean(dir_u)

  for (em in c("ipw", "reg", "dr")) {
    mp0 <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1, est_method = em, control_group = "notyettreated",
      gmm_weighting = "identity", pt_assumption = "ctrl", cross_pair = TRUE,
      data = df_cov, bstrap = FALSE, cband = FALSE
    )
    df_plus <- df_cov
    df_minus <- df_cov
    m_id <- match(df_cov$id, unit_ids)
    df_plus$w <- 1.0 + eps * dir_u[m_id]
    df_minus$w <- 1.0 - eps * dir_u[m_id]

    mp_p <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1, est_method = em, weightsname = "w", control_group = "notyettreated",
      gmm_weighting = "identity", pt_assumption = "ctrl", cross_pair = TRUE,
      data = df_plus, bstrap = FALSE, cband = FALSE
    )
    mp_m <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1, est_method = em, weightsname = "w", control_group = "notyettreated",
      gmm_weighting = "identity", pt_assumption = "ctrl", cross_pair = TRUE,
      data = df_minus, bstrap = FALSE, cband = FALSE
    )

    fd_deriv <- (mp_p$att - mp_m$att) / (2.0 * eps)
    if_mat <- as.matrix(mp0$inffunc)
    an_deriv <- as.double(crossprod(dir_u, if_mat)) / as.double(N)
    ok <- is.finite(mp0$att) & mp0$att != 0
    expect_lt(max(abs(fd_deriv[ok] - an_deriv[ok])), 1e-5)
  }

  # 2. Global exact collinearity works in OR when alias holds identically on target & controls
  df_alias_ok <- df_cov
  df_alias_ok$x2 <- 2.0 * df_alias_ok$x1
  mp_alias_ok <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = ~ x1 + x2, est_method = "reg", control_group = "notyettreated",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df_alias_ok, bstrap = FALSE, cband = FALSE
  )
  mp_single_x <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = ~ x1, est_method = "reg", control_group = "notyettreated",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df_cov, bstrap = FALSE, cband = FALSE
  )
  expect_equal(mp_alias_ok$att, mp_single_x$att, tolerance = 1e-12)

  # 3. Support-specific rank deficiency (x2 = 0 on controls, x2 != 0 on target G=3) raises cdid_or_wls_rank_error
  df_bad_or <- df_cov
  df_bad_or$x2 <- ifelse(df_bad_or$G == 3L, df_bad_or$x1 + 1.0, 0.0)
  expect_error(
    cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1 + x2, est_method = "reg", control_group = "notyettreated",
      cross_pair = TRUE, data = df_bad_or, bstrap = FALSE
    ),
    class = "cdid_or_wls_rank_error"
  )

  # 4. Complete covariate separation in propensity logit raises cdid_ipw_fit_error
  df_sep <- df_cov
  df_sep$x_sep <- ifelse(df_sep$G == 3L, 10.0, -10.0)
  expect_error(
    cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x_sep, est_method = "ipw", control_group = "notyettreated",
      cross_pair = TRUE, data = df_sep, bstrap = FALSE
    ),
    class = "cdid_ipw_fit_error"
  )

  # 5. Cluster-robust SE + multiplier bootstrap simultaneous confidence bands + workflow S3 methods
  set.seed(99)
  mp_boot <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = ~ x1, est_method = "dr", clustervars = "cluster",
    control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "optimal",
    data = df_cov, bstrap = TRUE, cband = TRUE, biters = 200
  )
  expect_s3_class(mp_boot, "cdid_MP")
  expect_true(is.finite(mp_boot$c))
  expect_gt(mp_boot$c, stats::qnorm(0.975) * 0.9)
  td <- cdid::tidy(mp_boot)
  gl <- cdid::glance(mp_boot)
  expect_equal(nrow(td), length(mp_boot$att))
  expect_equal(nrow(gl), 1L)
  for (tp in c("simple", "group", "dynamic", "calendar")) {
    ag <- cdid::aggte(mp_boot, type = tp, na.rm = TRUE, biters = 200)
    expect_s3_class(ag, "AGGTEobj")
    expect_true(is.finite(ag$overall.att))
  }
})

test_that("mask spanning-forest Gram G = E'E is exact on non-clique graphs and preserves Tikhonov-regularized GMM invariance", {
  ns <- asNamespace("cdid")
  build_fg <- get("build_mask_forest_and_gram", envir = ns)

  # 1. Non-clique 4-cycle on dates {1, 2, 3, 4}: edges (1,2), (2,3), (3,4), (1,4)
  # Missing chords (1,3) and (2,4) occur when extra units are observed at {1,3} and {2,4}.
  s_cyc <- c(1, 2, 3, 1)
  t_cyc <- c(2, 3, 4, 4)
  fg <- build_fg(s_cyc, t_cyc, rep("grp1", 4L))
  expect_equal(fg$kept, c(1L, 2L, 3L))
  # Dropped edge (1,4) = (1,2) + (2,3) + (3,4), so E = rbind(I_3, c(1, 1, 1))
  E_exact <- rbind(diag(3L), c(1.0, 1.0, 1.0))
  expect_equal(fg$G, crossprod(E_exact), tolerance = 1e-15)

  # 2. Exact Tikhonov invariance: reduced + G = E'E vs unreduced full catalog
  set.seed(5150)
  rows <- list()
  uid <- 0L
  for (g in c(0L, 4L, 5L, 6L)) {
    for (i in seq_len(35L)) {
      uid <- uid + 1L
      a_i <- stats::rnorm(1L)
      w_i <- stats::runif(1L, 0.5, 2.0)
      for (tv in 1:6) {
        eff <- if (g > 0L && tv >= g) 1.8 + 0.3 * (tv - g) else 0.0
        rows[[length(rows) + 1L]] <- data.frame(
          id = uid, period = tv, G = g, w = w_i,
          y = a_i + 0.25 * tv + eff + stats::rnorm(1L, sd = 0.3)
        )
      }
    }
  }
  df_bal <- do.call(rbind, rows)

  # Force covariance-ridge Tikhonov activation via condition_trigger = 0.99 and compare
  # spanning-forest reduced fit against unreduced fit on cross-pair!
  run_cp <- get("run_chained_cross_pair_fit", envir = ns)
  resolve_ctrl <- get("chained_public_resolve_control", envir = ns)
  prep_inputs <- get("chained_request_prepare_inputs", envir = ns)
  pub_dims <- get("chained_public_dimensions", envir = ns)
  pub_ctrl <- get("chained_public_pipeline_control", envir = ns)

  res_ctrl <- resolve_ctrl(
    list(condition_trigger = 0.99, ridge_relative = 1e-2),
    gmm_weighting = "within", gmm_covariance = "tikhonov"
  )
  prep <- prep_inputs(
    data = df_bal, yname = "y", tname = "period", idname = "id",
    gname = "G", xformla = ~ 1, weightsname = "w"
  )
  dims <- pub_dims(prep, biters = 100, clustered = FALSE)
  produced <- pub_ctrl(
    dims, res_ctrl, gmm_weighting = "within", gmm_covariance = "tikhonov", workers = 1L
  )

  res_red <- run_cp(
    data = prep$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = produced$control, xformla = prep$xformla, cluster_vec = prep$cluster_full,
    target_cohorts = prep$target_cohorts, anticipation = 0L,
    control_group = "notyettreated", base_period = "universal",
    weighting = "within", pt_assumption = "ctrl", pt_start = NULL,
    est_method = "ipw", cband = FALSE, alp = 0.05, biters = 100L,
    cp04 = prep$cp04, frame_id = prep$frame_ids, frame_Z = prep$frame_Z,
    frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = FALSE, bstrap = FALSE, gmm_covariance = "tikhonov",
    reduce_structural_redundancy = TRUE, prepared = prep$prepared
  )
  res_full <- run_cp(
    data = prep$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = produced$control, xformla = prep$xformla, cluster_vec = prep$cluster_full,
    target_cohorts = prep$target_cohorts, anticipation = 0L,
    control_group = "notyettreated", base_period = "universal",
    weighting = "within", pt_assumption = "ctrl", pt_start = NULL,
    est_method = "ipw", cband = FALSE, alp = 0.05, biters = 100L,
    cp04 = prep$cp04, frame_id = prep$frame_ids, frame_Z = prep$frame_Z,
    frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = FALSE, bstrap = FALSE, gmm_covariance = "tikhonov",
    reduce_structural_redundancy = FALSE, prepared = prep$prepared
  )
  expect_true(res_red$gmm_fit$covariance_handling$handler_invoked)
  expect_true(res_full$gmm_fit$covariance_handling$handler_invoked)
  expect_false(res_red$gmm_fit$covariance_handling$gram_is_identity)
  expect_equal(res_red$reporting$target_meta$estimate,
               res_full$reporting$target_meta$estimate, tolerance = 1e-12)

  # Same-pair vs cross-pair agreement on balanced panel + partial pt_start
  mp_same_within_dr <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G", weightsname = "w",
    est_method = "dr", control_group = "notyettreated", base_period = "universal",
    pt_assumption = "ctrl", gmm_weighting = "within", cross_pair = FALSE,
    data = df_bal, bstrap = FALSE, cband = FALSE
  )
  mp_cp_within <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G", weightsname = "w",
    control_group = "notyettreated", base_period = "universal",
    pt_assumption = "ctrl", gmm_weighting = "within", cross_pair = TRUE,
    data = df_bal, bstrap = FALSE, cband = FALSE
  )
  post_mask <- mp_cp_within$t >= mp_cp_within$group
  expect_lt(max(abs(mp_same_within_dr$att[post_mask] - mp_cp_within$att[post_mask])), 1e-12)

  mp_partial_same <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", pt_assumption = "ctrl", pt_start = 2,
    cross_pair = FALSE, data = df_bal, bstrap = FALSE, cband = FALSE
  )
  mp_partial_cp <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", pt_assumption = "ctrl", pt_start = 2,
    cross_pair = TRUE, data = df_bal, bstrap = FALSE, cband = FALSE
  )
  expect_equal(mp_partial_same$chained$pt$pt_start, 2)
  expect_equal(mp_partial_cp$chained$pt$pt_start, 2)
  expect_false(anyNA(mp_partial_same$att))
  expect_false(anyNA(mp_partial_cp$att))
})

test_that("selective post-treatment attrition biases cohort-FE slice OLS while leaving within-unit cross-pair and same-pair GMM unbiased", {
  set.seed(888)
  rows <- list()
  uid <- 0L
  for (g in c(2005L, 2007L, 2009L)) {
    for (i in seq_len(400L)) {
      uid <- uid + 1L
      a_i <- stats::rnorm(1L, sd = 2.0)
      for (tv in 2003:2008) {
        if (g == 2005L && tv >= 2006L && a_i < 0) next
        eff <- if (tv >= g) 5.0 else 0.0
        rows[[length(rows) + 1L]] <- data.frame(
          id = uid, period = tv, G = g, y = a_i + 0.5 * (tv - 2003L) + eff + stats::rnorm(1L, sd = 0.05)
        )
      }
    }
  }
  df_attr <- do.call(rbind, rows)
  ols_coh <- fit_reference_slice_ols(df_attr, fe = "cohort")
  mp_cp <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df_attr, bstrap = FALSE, cband = FALSE
  )
  est_coh <- ols_coh$att_ols[ols_coh$group == 2005 & ols_coh$t == 2008]
  est_cdid <- mp_cp$att[mp_cp$group == 2005 & mp_cp$t == 2008]
  expect_gt(abs(est_coh - 5.0), 1.0)
  expect_lt(abs(est_cdid - 5.0), 0.05)
})

test_that("chain_report S3 diagnostics, include_pre = TRUE witnesses/weights, unreduced direct-link detection, frame-only units, and compute_inffunc = FALSE work cleanly", {
  df <- make_short_window_panel(seed = 123L)
  mp_cp <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    data = df, bstrap = FALSE, cband = FALSE
  )
  mp_cp_var <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", base_period = "varying",
    cross_pair = TRUE, pt_assumption = "post", gmm_weighting = "optimal",
    data = df, bstrap = FALSE, cband = FALSE
  )
  mp_sp <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", base_period = "universal",
    cross_pair = FALSE, pt_assumption = "post",
    data = df, bstrap = FALSE, cband = FALSE
  )

  # 1. Direct vs chained support detection on unreduced catalog:
  # Even though spanning-forest reduction keeps consecutive 1-step edges (2004->2005->2006),
  # ATT(2005, 2006) IS directly supported in the unreduced catalog by cohorts 2007..2010 observed at (2004, 2006)!
  cr_cp <- cdid::chain_report(mp_cp, cohort = 2005, time = c(2005, 2006, 2009))
  expect_s3_class(cr_cp, "cdid_chain_report")
  out_txt <- utils::capture.output(print(cr_cp, show_gmm_weights = TRUE))
  expect_true(any(grepl("ATT\\(g=2005,t=2009\\)", out_txt)))
  expect_true(any(grepl("gmm_weight", out_txt)))

  s_cp <- as.data.frame(cr_cp, what = "summary")
  w_cp <- as.data.frame(cr_cp, what = "witness")
  g_cp <- as.data.frame(cr_cp, what = "gmm_weights")
  expect_equal(nrow(s_cp), 3L)
  expect_false(s_cp$requires_chaining[s_cp$time == 2005])
  expect_false(s_cp$requires_chaining[s_cp$time == 2006])
  expect_true(s_cp$requires_chaining[s_cp$time == 2009])
  expect_gt(nrow(w_cp), 2L)
  for (tid in s_cp$target_id) {
    est_val <- s_cp$estimate[s_cp$target_id == tid]
    gw_sum <- sum(g_cp$contribution[g_cp$target_id == tid])
    expect_lt(abs(est_val - gw_sum), 1e-11)
  }

  # 2. include_pre = TRUE populates witnesses and exact GMM moment weights in BOTH cross-pair & same-pair!
  for (fit_obj in list(mp_cp, mp_cp_var, mp_sp)) {
    cr_pre <- cdid::chain_report(fit_obj, cohort = 2009, include_pre = TRUE)
    s_pre <- as.data.frame(cr_pre, what = "summary")
    w_pre <- as.data.frame(cr_pre, what = "witness")
    g_pre <- as.data.frame(cr_pre, what = "gmm_weights")
    pre_rows <- s_pre[s_pre$time < s_pre$g & s_pre$available & abs(s_pre$estimate) > 1e-12, , drop = FALSE]
    expect_gt(nrow(pre_rows), 0L)
    for (tid in pre_rows$target_id) {
      expect_gt(sum(w_pre$target_id == tid), 0L)
      gw_sub <- g_pre[g_pre$target_id == tid, , drop = FALSE]
      expect_gt(nrow(gw_sub), 0L)
      expect_lt(abs(sum(gw_sub$contribution) - pre_rows$estimate[pre_rows$target_id == tid]), 1e-10)
    }
  }

  # 3. Explicit unit frame with frame-only units + point-only compute_inffunc = FALSE contract
  u_df <- df[!duplicated(df$id), c("id", "G", "w", "cluster"), drop = FALSE]
  frame_extra <- rbind(
    u_df,
    data.frame(id = c(999901L, 999902L), G = c(2005L, 2008L), w = c(1.2, 0.8), cluster = c(500L, 501L))
  )
  mp_frame <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    frame = frame_extra, weightsname = "w", clustervars = "cluster",
    control_group = "notyettreated", cross_pair = TRUE, pt_assumption = "ctrl",
    data = df, bstrap = FALSE, cband = FALSE
  )
  expect_equal(mp_frame$n, nrow(frame_extra))
  expect_equal(nrow(mp_frame$inffunc), nrow(frame_extra))

  mp_no_if <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", cross_pair = TRUE, pt_assumption = "ctrl",
    gmm_weighting = "within", compute_inffunc = FALSE,
    data = df, bstrap = FALSE, cband = FALSE
  )
  expect_null(mp_no_if$inffunc)
  mp_with_if <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", cross_pair = TRUE, pt_assumption = "ctrl",
    gmm_weighting = "within", compute_inffunc = TRUE,
    data = df, bstrap = FALSE, cband = FALSE
  )
  expect_equal(mp_no_if$att, mp_with_if$att, tolerance = 1e-14)
  expect_true(all(is.na(mp_no_if$se)))
  expect_error(cdid::aggte(mp_no_if, type = "simple", na.rm = TRUE),
               class = "cdid_chained_request_no_inffunc_error")

  # 4. Scope-gate refusals (custom est_method callbacks remain refused on cross_pair and within)
  custom_cb <- function(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...) {
    n <- length(y1)
    list(ATT = 0.0, att.inf.func = numeric(n))
  }
  expect_error(
    cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                 est_method = custom_cb, cross_pair = TRUE, data = df),
    class = "cdid_chained_request_input_error"
  )
  expect_error(
    cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                 est_method = custom_cb, gmm_weighting = "within", data = df),
    class = "cdid_chained_request_input_error"
  )
  expect_error(
    cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                 pt_assumption = "ctrl", pt_start = 1999, data = df),
    class = "cdid_error"
  )
  expect_error(
    cdid::chain_report(list(), cohort = 2005),
    class = "cdid_chain_report_input_error"
  )
  expect_error(
    cdid::chain_report(mp_cp, cohort = 2005, group = 2005),
    class = "cdid_chain_report_input_error"
  )
})

test_that("inverse pair-rate propensity calibration, rate-channel IF, rank/alias guards, and unqualified direct-link chain_report detection hold exactly", {
  # 1. Exact point-estimate, SE, and unit-by-unit IF equivalence between cross_pair = TRUE and
  #    cross_pair = FALSE when r_g = 1 and |K| = 1 under unbalanced control pair rate r_k < 1,
  #    non-uniform weights, singletons (N > n_A), and xformla with or without intercept
  set.seed(20261112)
  n_per <- 30L
  ids <- seq_len(2L * n_per)
  g_u <- rep(c(4L, 6L), each = n_per)
  x1_u <- stats::rnorm(length(ids), mean = ifelse(g_u == 4L, 0.35, -0.15), sd = 0.6)
  x2_u <- stats::runif(length(ids), 0.4, 1.6)
  w_u <- stats::runif(length(ids), 0.6, 2.1)

  rows <- list()
  for (i in seq_along(ids)) {
    # Cohort 4 completely observed on {1,2,3,4} (r_4 = 1); cohort 6 split across
    # {1,2,3} (r_6 = 1/3 < 1) and {3,4} (r_6 = 2/3 < 1) with zero cross-stratum cycles
    t_obs <- if (g_u[i] == 4L) {
      c(1L, 2L, 3L, 4L)
    } else if (i %% 3L == 0L) {
      c(1L, 2L, 3L)
    } else {
      c(3L, 4L)
    }
    a_i <- stats::rnorm(1L)
    for (tv in t_obs) {
      yv <- a_i + 0.5 * tv + 0.8 * x1_u[i] * tv - 0.3 * x2_u[i] * (tv^1.2) +
        ifelse(tv >= g_u[i], 3.0, 0.0) + stats::rnorm(1L, sd = 0.35)
      rows[[length(rows) + 1L]] <- data.frame(
        id = ids[i], period = tv, G = g_u[i], x1 = x1_u[i], x2 = x2_u[i], w = w_u[i], y = yv
      )
    }
  }
  # Append singleton units (observed at a single date only, so N > n_A)
  df_rk <- rbind(
    do.call(rbind, rows),
    data.frame(
      id = 61:64, period = c(2L, 3L, 3L, 4L), G = c(4L, 4L, 6L, 6L),
      x1 = stats::rnorm(4L), x2 = stats::runif(4L, 0.4, 1.6), w = c(0.8, 1.4, 1.1, 1.7), y = stats::rnorm(4L)
    )
  )

  for (xf in list(~ x1 + x2, ~ 0 + x1 + x2)) {
    for (em in c("ipw", "reg", "dr")) {
      for (wt in c("identity", "within")) {
        sp_fit <- cdid::att_gt(
          yname = "y", tname = "period", idname = "id", gname = "G",
          xformla = xf, weightsname = "w", est_method = em,
          control_group = "notyettreated", base_period = "universal",
          cross_pair = FALSE, pt_assumption = "post", gmm_weighting = wt,
          data = df_rk, bstrap = FALSE, cband = FALSE
        )
        cp_fit <- cdid::att_gt(
          yname = "y", tname = "period", idname = "id", gname = "G",
          xformla = xf, weightsname = "w", est_method = em,
          control_group = "notyettreated", base_period = "universal",
          cross_pair = TRUE, pt_assumption = "post", gmm_weighting = wt,
          data = df_rk, bstrap = FALSE, cband = FALSE
        )
        expect_equal(cp_fit$att, sp_fit$att, tolerance = 1e-11)
        expect_equal(cp_fit$se, sp_fit$se, tolerance = 1e-11)
        expect_equal(unname(as.matrix(cp_fit$inffunc)), unname(as.matrix(sp_fit$inffunc)), tolerance = 1e-10)
      }
    }
  }

  # 2. Finite-difference case-weight derivative check under r_k < 1, ~ 0 + x1 + x2, singletons, and aggte()
  u_ids <- as.integer(rownames(cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = ~ 0 + x1 + x2, weightsname = "w", est_method = "dr",
    control_group = "notyettreated", cross_pair = TRUE, pt_assumption = "ctrl",
    gmm_weighting = "identity", data = df_rk, bstrap = FALSE, cband = FALSE
  )$inffunc))
  N_tot <- length(u_ids)
  set.seed(404)
  dir_w <- stats::rnorm(N_tot)
  dir_w <- dir_w - mean(dir_w)
  eps_w <- 1e-6
  m_idx <- match(df_rk$id, u_ids)
  df_p <- df_rk; df_p$w <- df_rk$w * (1.0 + eps_w * dir_w[m_idx])
  df_m <- df_rk; df_m$w <- df_rk$w * (1.0 - eps_w * dir_w[m_idx])

  for (em in c("ipw", "reg", "dr")) {
    f0 <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ 0 + x1 + x2, weightsname = "w", est_method = em,
      control_group = "notyettreated", base_period = "universal",
      cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "identity",
      data = df_rk, bstrap = FALSE, cband = FALSE
    )
    fp <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ 0 + x1 + x2, weightsname = "w", est_method = em,
      control_group = "notyettreated", base_period = "universal",
      cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "identity",
      data = df_p, bstrap = FALSE, cband = FALSE
    )
    fm <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ 0 + x1 + x2, weightsname = "w", est_method = em,
      control_group = "notyettreated", base_period = "universal",
      cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "identity",
      data = df_m, bstrap = FALSE, cband = FALSE
    )
    av <- which(is.finite(f0$att) & f0$chained$target_meta$role != "reference")
    fd_att <- (fp$att[av] - fm$att[av]) / (2.0 * eps_w)
    an_att <- as.double(crossprod(dir_w, as.matrix(f0$inffunc)[, av, drop = FALSE])) / as.double(N_tot)
    expect_lt(max(abs(fd_att - an_att)), 1e-5)

    ag0 <- cdid::aggte(f0, type = "simple", na.rm = TRUE, bstrap = FALSE, cband = FALSE)
    agp <- cdid::aggte(fp, type = "simple", na.rm = TRUE, bstrap = FALSE, cband = FALSE)
    agm <- cdid::aggte(fm, type = "simple", na.rm = TRUE, bstrap = FALSE, cband = FALSE)
    fd_ag <- (agp$overall.att - agm$overall.att) / (2.0 * eps_w)
    an_ag <- sum(dir_w * as.double(ag0$inf.function$simple.att)) / as.double(N_tot)
    expect_lt(abs(fd_ag - an_ag), 1e-5)
  }

  # 3. Pre-GLM exact alias refusal on IPW/DR and near-collinear non-exact alias refusal on OR
  df_dup <- df_rk
  df_dup$x_dup <- 2.0 * df_dup$x1
  for (em_ipw in c("ipw", "dr")) {
    expect_error(
      cdid::att_gt(
        yname = "y", tname = "period", idname = "id", gname = "G",
        xformla = ~ x1 + x_dup, est_method = em_ipw,
        control_group = "notyettreated", cross_pair = TRUE, data = df_dup, bstrap = FALSE
      ),
      class = "cdid_ipw_rank_error"
    )
  }
  df_near <- df_rk
  df_near$x_near <- 2.0 * df_near$x1 + 1e-9 * df_near$x2
  expect_error(
    cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1 + x_near, est_method = "reg",
      control_group = "notyettreated", cross_pair = TRUE, data = df_near, bstrap = FALSE
    ),
    class = "cdid_or_wls_rank_error"
  )

  # 4. chain_report() requires_chaining = TRUE when no qualified direct (ref, t) link exists
  set.seed(303)
  rows_uq <- list()
  uid_uq <- 0L
  for (spec_u in list(
    list(g = 3L, times = c(1L, 2L, 3L, 4L)),
    list(g = 4L, times = c(1L, 2L, 3L, 4L)),
    list(g = 5L, times = c(1L, 3L, 4L))
  )) {
    for (i in seq_len(20L)) {
      uid_uq <- uid_uq + 1L
      a_i <- stats::rnorm(1L)
      for (tv in spec_u$times) {
        yv <- a_i + 0.4 * tv + ifelse(tv >= spec_u$g, 2.0, 0.0) + stats::rnorm(1L, sd = 0.1)
        rows_uq[[length(rows_uq) + 1L]] <- data.frame(id = uid_uq, period = tv, G = spec_u$g, y = yv)
      }
    }
  }
  df_uq <- do.call(rbind, rows_uq)
  fit_uq <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    control_group = "notyettreated", base_period = "universal",
    cross_pair = FALSE, pt_assumption = "post", data = df_uq, bstrap = FALSE, cband = FALSE
  )
  al_uq <- fit_uq$chained$chain_info$all_links
  expect_false(any(al_uq$g == 3 & al_uq$s == 2 & al_uq$t == 4 & al_uq$t_qualified))
  cr_uq <- cdid::chain_report(fit_uq, cohort = 3, time = 4)
  sum_uq <- as.data.frame(cr_uq, what = "summary")
  expect_true(sum_uq$available)
  expect_true(sum_uq$requires_chaining)

  # Even if an unqualified (t_qualified == FALSE) row spans (2, 4), requires_chaining remains TRUE
  fit_uq_mod <- fit_uq
  fit_uq_mod$chained$chain_info$all_links <- rbind(
    fit_uq_mod$chained$chain_info$all_links,
    data.frame(g = 3, s = 2, t = 4, link_role = "did_pair", control_pool = "{4}", t_qualified = FALSE, stringsAsFactors = FALSE)
  )
  cr_uq_mod <- cdid::chain_report(fit_uq_mod, cohort = 3, time = 4)
  expect_true(as.data.frame(cr_uq_mod, what = "summary")$requires_chaining)
})

test_that("fractional dates, anticipation, partial-window PT-ctrl, and restrict_chained_design_pt_ctrl work cleanly", {
  ns <- asNamespace("cdid")

  # 1. Fractional calendar dates & cohorts + weights
  set.seed(401)
  rows_fr <- list()
  uid_fr <- 0L
  for (g in c(0, 2003.5, 2004.5)) {
    for (i in seq_len(20L)) {
      uid_fr <- uid_fr + 1L
      x_i <- stats::rnorm(1L)
      w_i <- 0.5 + stats::runif(1L)
      a_i <- stats::rnorm(1L)
      for (tv in c(2001.5, 2002.5, 2003.5, 2004.5)) {
        tau <- if (g > 0 && tv >= g) 2.5 else 0.0
        yv <- a_i + 0.4 * (tv - 2001.5) * (1 + 0.3 * x_i) + tau + stats::rnorm(1L, sd = 0.15)
        rows_fr[[length(rows_fr) + 1L]] <- data.frame(id = uid_fr, period = tv, G = g, x = x_i, w = w_i, y = yv)
      }
    }
  }
  df_fr <- do.call(rbind, rows_fr)
  f_fr_cp <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = ~ x, weightsname = "w", data = df_fr,
    control_group = "nevertreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within",
    bstrap = FALSE, cband = FALSE
  )
  f_fr_sp <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = ~ x, weightsname = "w", data = df_fr,
    control_group = "nevertreated", base_period = "universal",
    cross_pair = FALSE, pt_assumption = "ctrl", gmm_weighting = "within",
    bstrap = FALSE, cband = FALSE
  )
  expect_lt(max(abs(f_fr_cp$att - f_fr_sp$att)), 1e-11)
  expect_lt(max(abs(f_fr_cp$se - f_fr_sp$se)), 1e-11)

  # 2. Anticipation > 0 matches same-pair identically under nevertreated
  set.seed(402)
  rows_ant <- list()
  uid_ant <- 0L
  for (g in c(0, 5, 7)) {
    for (i in seq_len(25L)) {
      uid_ant <- uid_ant + 1L
      a_i <- stats::rnorm(1L)
      for (tv in 1:7) {
        tau <- if (g > 0 && tv >= g) 3.0 else if (g > 0 && tv == g - 1L) 0.8 else 0.0
        rows_ant[[length(rows_ant) + 1L]] <- data.frame(
          id = uid_ant, period = tv, G = g, y = a_i + 0.5 * tv + tau + stats::rnorm(1L, sd = 0.1)
        )
      }
    }
  }
  df_ant <- do.call(rbind, rows_ant)
  for (bp in c("universal", "varying")) {
    f_ant_cp <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      data = df_ant, anticipation = 1, control_group = "nevertreated", base_period = bp,
      cross_pair = TRUE, pt_assumption = "post", gmm_weighting = "within",
      bstrap = FALSE, cband = FALSE
    )
    f_ant_sp <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      data = df_ant, anticipation = 1, control_group = "nevertreated", base_period = bp,
      cross_pair = FALSE, pt_assumption = "post", gmm_weighting = "within",
      bstrap = FALSE, cband = FALSE
    )
    expect_lt(max(abs(f_ant_cp$att - f_ant_sp$att)), 1e-12)
    expect_lt(max(abs(f_ant_cp$se - f_ant_sp$se)), 1e-12)
  }

  # 3. Partial-window PT-ctrl (pt_start = 2 > min(calendar) = 1) eliminates pre-pt_start control shock bias
  set.seed(403)
  rows_pw <- list()
  uid_pw <- 0L
  for (sp in list(
    list(g = 4, win = c(3, 4), sh = 0.0),
    list(g = 6, win = c(1, 2, 3, 4, 5), sh = 10.0),
    list(g = 7, win = c(1, 2, 4, 5), sh = 0.0)
  )) {
    for (i in seq_len(30L)) {
      uid_pw <- uid_pw + 1L
      for (tv in sp$win) {
        yv <- 0.3 * tv + (if (sp$g == 4 && tv >= 4) 5.0 else 0.0) + (if (tv == 1) sp$sh else 0.0) + stats::rnorm(1L, sd = 1e-5)
        rows_pw[[length(rows_pw) + 1L]] <- data.frame(id = uid_pw, period = tv, G = sp$g, y = yv)
      }
    }
  }
  df_pw <- do.call(rbind, rows_pw)
  fit_pw1 <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    data = df_pw, control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", pt_start = NULL,
    gmm_weighting = "within", bstrap = FALSE, cband = FALSE
  )
  fit_pw2 <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    data = df_pw, control_group = "notyettreated", base_period = "universal",
    cross_pair = TRUE, pt_assumption = "ctrl", pt_start = 2,
    gmm_weighting = "within", bstrap = FALSE, cband = FALSE
  )
  idx_pw44 <- which(fit_pw1$group == 4 & fit_pw1$t == 4)
  expect_gt(abs(fit_pw1$att[idx_pw44] - 5.0), 0.5)
  expect_lt(abs(fit_pw2$att[idx_pw44] - 5.0), 1e-4)

  # Also check pt_operator_covered is FALSE on unavailable pre-trend cells (cohort 4 unobserved at t=1,2)
  tm_pw1 <- fit_pw1$chained$target_meta
  unav_pre <- which(tm_pw1$g == 4 & tm_pw1$role == "pre" & !tm_pw1$available)
  expect_gt(length(unav_pre), 0L)
  expect_true(all(!tm_pw1$pt_operator_covered[unav_pre]))

  # 4. restrict_chained_design_pt_ctrl matches direct PT-ctrl design construction
  links_in <- data.frame(
    g = c(4, 4, 4, 4), s = c(1, 1, 2, 3), t = c(2, 2, 3, 4), ref = c(3, 3, 3, 3),
    control_pool = c("{5,6}", "{6}", "{5,6}", "{5,6}"),
    functional_key = c("p56", "p6", "p56", "p56"),
    s_qualified = c(FALSE, FALSE, FALSE, TRUE),
    t_qualified = c(FALSE, FALSE, TRUE, TRUE),
    s_coordinate_role = c("pre_nuisance", "pre_nuisance", "pre_nuisance", "reference"),
    t_coordinate_role = c("pre_nuisance", "pre_nuisance", "reference", "post_att"),
    stringsAsFactors = FALSE
  )
  attr(links_in, "calendar_spec") <- ns$validate_calendar_spec(1:6, c(4, 5, 6), 0, pt_assumption = "post")
  d_post <- ns$build_chained_design(links_in)
  d_restr <- ns$restrict_chained_design_pt_ctrl(d_post, pt_start = NULL)
  links_c <- links_in
  links_c$s_coordinate_role[1:3] <- "pre_pt_ctrl_shared"
  links_c$t_coordinate_role[1:2] <- "pre_pt_ctrl_shared"
  attr(links_c, "calendar_spec") <- ns$validate_calendar_spec(1:6, c(4, 5, 6), 0, pt_assumption = "ctrl", pt_start = 1)
  d_direct <- ns$build_chained_design(links_c)
  expect_identical(d_restr$incidence, d_direct$incidence)
  expect_identical(d_restr$pt_ctrl_spec, d_direct$pt_ctrl_spec)
  expect_identical(d_restr$pt_ctrl_spec$n_merged_pre_parameters, 2L)
})


test_that("cross_pair and within support all IPW/OR/DR, outcome_pooling, conditional_mtar selection regimes, reference inputs, and full-roster inference", {
  # Build a compact rich panel with F \ B (never observed), B \ A (singletons), and A (repeat observers)
  set.seed(20261120)
  cohs <- c(4, 5, Inf)
  times <- 1:5
  n_per <- 36L
  N <- length(cohs) * n_per
  id <- seq_len(N)
  g_true <- rep(cohs, each = n_per)
  x1 <- stats::rnorm(N, mean = ifelse(g_true == 4, 0.25, ifelse(g_true == 5, -0.1, -0.15)), sd = 0.75)
  x2 <- stats::rnorm(N, sd = 0.75)
  v_sel <- stats::rnorm(N, sd = 0.65)
  w_f <- stats::runif(N, 0.6, 1.8)
  clust <- sample.int(18L, N, replace = TRUE)
  alpha_i <- 0.5 * x1 - 0.3 * x2 + stats::rnorm(N, sd = 0.45)

  S_mat <- matrix(FALSE, nrow = N, ncol = length(times))
  for (ci in seq_along(cohs)) {
    idx <- which(g_true == cohs[ci])
    idx_fb <- idx[1:4]
    idx_ba <- idx[5:8]
    idx_a  <- idx[9:n_per]
    for (i in idx_ba) S_mat[i, sample.int(length(times), 1L)] <- TRUE
    for (i in idx_a) {
      pr <- stats::plogis(0.75 + 0.3 * x1[i] - 0.25 * v_sel[i])
      s_row <- stats::runif(length(times)) < pr
      if (sum(s_row) < 2L) s_row[sample.int(length(times), 2L)] <- TRUE
      S_mat[i, ] <- s_row
    }
  }

  frame_df <- data.frame(
    id = id,
    g = ifelse(is.infinite(g_true), 0, g_true),
    x1 = x1, x2 = x2, x_dup = 2.0 * x1 - 0.5 * x2, v_sel = v_sel,
    w = w_f, clust = clust,
    stringsAsFactors = FALSE
  )
  in_B <- rowSums(S_mat) >= 1L
  in_A <- rowSums(S_mat) >= 2L

  rho_B <- rep(NA_real_, N)
  rho_A <- rep(NA_real_, N)
  raw_b <- exp(0.18 * x1)
  raw_a <- exp(-0.15 * x2)
  for (gv in unique(frame_df$g)) {
    mb <- frame_df$g == gv & in_B
    ma <- frame_df$g == gv & in_A
    mf <- frame_df$g == gv
    Wb <- sum(w_f[mb]); Wa <- sum(w_f[ma]); Wf <- sum(w_f[mf])
    rho_B[mb] <- raw_b[mb] * sum(w_f[mb] / raw_b[mb]) / Wb
    rho_A[ma] <- raw_a[ma] * sum(w_f[ma] / raw_a[ma]) / Wa
    # Note: cdid checks sum_{i in U_k} w_i / rho_i / sum_{i in U_k} w_i == 1
  }
  frame_df$rho_B <- rho_B
  frame_df$rho_A <- rho_A

  p_rows <- list()
  for (i in which(in_B)) {
    for (ti in which(S_mat[i, ])) {
      tt <- times[ti]
      trt <- if (is.finite(g_true[i]) && tt >= g_true[i]) (2.0 + 0.4 * (tt - g_true[i])) else 0.0
      y_it <- alpha_i[i] + 0.3 * tt + 0.35 * x1[i] * tt - 0.2 * x2[i] * tt + trt + stats::rnorm(1L, sd = 0.3)
      p_rows[[length(p_rows) + 1L]] <- data.frame(
        id = id[i], t = tt, g = frame_df$g[i], y = y_it,
        x1 = x1[i], x2 = x2[i], x_dup = frame_df$x_dup[i], v_sel = v_sel[i],
        w = w_f[i], clust = clust[i], stringsAsFactors = FALSE
      )
    }
  }
  panel_df <- do.call(rbind, p_rows)
  rownames(panel_df) <- NULL

  check_fd <- function(fit_fn, eta_dir = NULL, wt = "identity", rho_col = NULL, u_mask = NULL) {
    if (is.null(eta_dir)) {
      eta_dir <- stats::rnorm(N)
      eta_dir <- eta_dir - mean(eta_dir)
    }
    if (!is.null(rho_col) && !is.null(u_mask)) {
      # Project perturbation orthogonal to the fixed-rho cohort normalization constraint
      for (gv in unique(frame_df$g)) {
        mk <- which(frame_df$g == gv & u_mask)
        c_vec <- frame_df$w[mk] * (1.0 / frame_df[[rho_col]][mk] - 1.0)
        if (sum(c_vec^2) > 0) {
          eta_dir[mk] <- eta_dir[mk] - c_vec * (sum(c_vec * eta_dir[mk]) / sum(c_vec^2))
        }
      }
    }
    f0 <- fit_fn(panel_df, frame_df, wt)
    expect_identical(nrow(f0$inffunc), N)
    eps <- 1e-6
    f_p <- frame_df; f_p$w <- frame_df$w * (1.0 + eps * eta_dir)
    f_m <- frame_df; f_m$w <- frame_df$w * (1.0 - eps * eta_dir)
    p_p <- panel_df; p_p$w <- f_p$w[match(panel_df$id, f_p$id)]
    p_m <- panel_df; p_m$w <- f_m$w[match(panel_df$id, f_m$id)]
    fp <- fit_fn(p_p, f_p, wt)
    fm <- fit_fn(p_m, f_m, wt)
    av <- which(is.finite(f0$att) & f0$chained$target_meta$role != "reference")
    expect_gt(length(av), 0L)
    fd_att <- (fp$att[av] - fm$att[av]) / (2.0 * eps)
    an_att <- as.double(crossprod(eta_dir, as.matrix(f0$inffunc)[, av, drop = FALSE])) / as.double(N)
    expect_lt(max(abs(fd_att - an_att)), 1e-5)
    f0
  }

  # 1. Cross-pair baseline (`selection = "none"`): pooled vs separate OR/DR + collinear OR
  for (op in c("pooled", "separate")) {
    for (em in c("reg", "dr")) {
      f_id <- check_fd(function(p, f, wt) {
        cdid::att_gt(
          yname = "y", tname = "t", idname = "id", gname = "g",
          xformla = ~ x1 + x2, data = p, frame = f, weightsname = "w",
          control_group = "notyettreated", est_method = em, outcome_pooling = op,
          cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = wt,
          bstrap = FALSE, cband = FALSE
        )
      })
      f_opt <- cdid::att_gt(
        yname = "y", tname = "t", idname = "id", gname = "g",
        xformla = ~ x1 + x2, data = panel_df, frame = frame_df, weightsname = "w",
        control_group = "notyettreated", est_method = em, outcome_pooling = op,
        cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "optimal",
        bstrap = FALSE, cband = FALSE
      )
      expect_true(any(is.finite(f_opt$att)))
      expect_true(any(is.finite(f_opt$se) & f_opt$se > 0))
    }
  }

  # Exact structural collinearity (`x_dup = 2*x1 - 0.5*x2`) is handled cleanly by OR
  f_col_or <- cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g",
    xformla = ~ x1 + x2 + x_dup, data = panel_df, frame = frame_df, weightsname = "w",
    control_group = "notyettreated", est_method = "reg", outcome_pooling = "pooled",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "optimal",
    bstrap = FALSE, cband = FALSE
  )
  expect_true(any(is.finite(f_col_or$att)))

  # 2. Cross-pair `selection = "conditional_mtar"` across IPW, imputation, and DR (augmented)
  # 2a. Observation weighting (`selection_est_method = "ipw"`), including off-A-only perturbation & universes B/A
  eta_off_a <- numeric(N)
  eta_off_a[!in_A] <- seq_len(sum(!in_A)) - mean(seq_len(sum(!in_A)))
  f_sel_ipw <- check_fd(function(p, f, wt) {
    cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g",
      xformla = ~ x1 + x2, selection_xformla = ~ x1 + v_sel,
      data = p, frame = f, weightsname = "w",
      control_group = "notyettreated", est_method = "dr", outcome_pooling = "pooled",
      selection = "conditional_mtar", selection_est_method = "ipw",
      selection_model = "independent_period_logit", selection_pooling = "within_cohort",
      cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = wt,
      bstrap = FALSE, cband = FALSE
    )
  }, eta_dir = eta_off_a)
  expect_gt(max(abs(as.matrix(f_sel_ipw$inffunc)[!in_A, ]), na.rm = TRUE), 1e-3)

  # Sub-universes B and A with known normalized density ratios (`rho_B`, `rho_A`) and legacy `"corollary1"`
  for (u_cfg in list(
    list(sel = "conditional_mtar", u = "B", rho = "rho_B", mask = in_B, sp = "separate"),
    list(sel = "corollary1", u = "A", rho = "rho_A", mask = in_A, sp = "across_cohorts")
  )) {
    check_fd(function(p, f, wt) {
      cdid::att_gt(
        yname = "y", tname = "t", idname = "id", gname = "g",
        xformla = ~ 0 + x1 + x2, selection_xformla = ~ x1 + v_sel,
        data = p, frame = f, weightsname = "w",
        control_group = "notyettreated", est_method = "dr", outcome_pooling = "separate",
        selection = u_cfg$sel, selection_est_method = "ipw",
        selection_model = "pair_logit", selection_pooling = u_cfg$sp,
        selection_universe = u_cfg$u, selection_rho = u_cfg$rho,
        selection_justification = "Simulation transport test",
        cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = wt,
        bstrap = FALSE, cband = FALSE
      )
    }, rho_col = u_cfg$rho, u_mask = u_cfg$mask)
  }

  # 2b. Pure actual-change imputation (`selection_est_method = "imputation"`) across pooling modes & IPW/OR/DR
  for (cfg_imp in list(
    list(em = "ipw", op = "pooled", ip = "separate", xf = ~ x1 + x2, off_a = TRUE),
    list(em = "reg", op = "pooled", ip = "within_cohort", xf = ~ x1 + x2, off_a = TRUE),
    list(em = "dr", op = "separate", ip = "within_cohort", xf = ~ 0 + x1 + x2, off_a = TRUE),
    list(em = "reg", op = "separate", ip = "across_cohorts", xf = ~ x1 + x2, off_a = FALSE)
  )) {
    f_imp <- check_fd(function(p, f, wt) {
      cdid::att_gt(
        yname = "y", tname = "t", idname = "id", gname = "g",
        xformla = cfg_imp$xf, data = p, frame = f, weightsname = "w",
        control_group = "notyettreated", est_method = cfg_imp$em, outcome_pooling = cfg_imp$op,
        selection = "conditional_mtar", selection_est_method = "imputation",
        imputation_pooling = cfg_imp$ip,
        cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = wt,
        bstrap = FALSE, cband = FALSE
      )
    })
    if (cfg_imp$off_a) {
      expect_gt(max(abs(as.matrix(f_imp$inffunc)[!in_A, ]), na.rm = TRUE), 1e-3)
    }
  }

  # 2c. Augmented missing-outcome correction (`selection_est_method = "dr"`) across IPW/OR/DR
  for (cfg_aug in list(
    list(em = "ipw", op = "pooled", ip = "separate", sm = "pair_logit", sp = "separate"),
    list(em = "reg", op = "pooled", ip = "within_cohort", sm = "pair_logit", sp = "within_cohort"),
    list(em = "dr", op = "separate", ip = "across_cohorts", sm = "independent_period_logit", sp = "across_cohorts")
  )) {
    f_aug <- check_fd(function(p, f, wt) {
      cdid::att_gt(
        yname = "y", tname = "t", idname = "id", gname = "g",
        xformla = ~ x1 + x2, selection_xformla = ~ x1 + v_sel,
        data = p, frame = f, weightsname = "w",
        control_group = "notyettreated", est_method = cfg_aug$em, outcome_pooling = cfg_aug$op,
        selection = "conditional_mtar", selection_est_method = "dr",
        imputation_pooling = cfg_aug$ip,
        selection_model = cfg_aug$sm, selection_pooling = cfg_aug$sp,
        cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = wt,
        bstrap = FALSE, cband = FALSE
      )
    })
    expect_gt(max(abs(as.matrix(f_aug$inffunc)[!in_A, ]), na.rm = TRUE), 1e-3)
  }

  # 3. `gmm_weighting = "within"` on both cross_pair = FALSE and cross_pair = TRUE across all selection modes
  for (xp_flag in c(FALSE, TRUE)) {
    for (sem_w in c("none", "ipw", "imputation", "dr")) {
      f_w <- cdid::att_gt(
        yname = "y", tname = "t", idname = "id", gname = "g",
        xformla = ~ x1 + x2, data = panel_df, frame = frame_df, weightsname = "w",
        control_group = "notyettreated", est_method = "dr", outcome_pooling = "pooled",
        selection = if (sem_w == "none") "none" else "conditional_mtar",
        selection_est_method = if (sem_w == "none") "ipw" else sem_w,
        cross_pair = xp_flag, pt_assumption = "ctrl", gmm_weighting = "within",
        bstrap = FALSE, cband = FALSE
      )
      expect_identical(nrow(f_w$inffunc), N)
      expect_true(any(is.finite(f_w$att)))
      expect_true(any(is.finite(f_w$se) & f_w$se > 0))
    }
  }

  # 4. Time-varying reference-date covariates & anchored varying weights (`fix_weights`) + `aggte()` IF FD
  tv_rows <- list()
  for (i in seq_len(N)) {
    obs_y <- S_mat[i, ]
    if (sum(obs_y) < 2L) obs_y[1:2] <- TRUE
    for (ti in seq_along(times)) {
      tt <- times[ti]
      xtv <- x1[i] + 0.12 * tt + 0.08 * cos(i + tt)
      wtv <- w_f[i] * (1.0 + 0.07 * tt)
      trt <- if (is.finite(g_true[i]) && tt >= g_true[i]) 2.2 else 0.0
      y_it <- if (obs_y[ti]) (alpha_i[i] + 0.25 * tt + 0.3 * xtv + trt + stats::rnorm(1L, sd = 0.25)) else NA_real_
      tv_rows[[length(tv_rows) + 1L]] <- data.frame(
        id = id[i], t = tt, g = frame_df$g[i], y = y_it, xtv = xtv, wtv = wtv, stringsAsFactors = FALSE
      )
    }
  }
  tv_df <- do.call(rbind, tv_rows)
  eta_tv <- stats::rnorm(N); eta_tv <- eta_tv - mean(eta_tv)
  eps_tv <- 1e-6
  tv_p <- tv_df; tv_p$wtv <- tv_df$wtv * (1.0 + eps_tv * eta_tv[tv_df$id])
  tv_m <- tv_df; tv_m$wtv <- tv_df$wtv * (1.0 - eps_tv * eta_tv[tv_df$id])

  for (fw in c("first_period", "base_period")) {
    for (sem_tv in c("none", "dr")) {
      run_tv_fit <- function(d_in) {
        cdid::att_gt(
          yname = "y", tname = "t", idname = "id", gname = "g",
          xformla = ~ xtv, data = d_in, weightsname = "wtv", fix_weights = fw,
          control_group = "notyettreated", est_method = "dr", outcome_pooling = "pooled",
          selection = if (sem_tv == "none") "none" else "conditional_mtar",
          selection_est_method = if (sem_tv == "none") "ipw" else sem_tv,
          cross_pair = TRUE, pt_assumption = "post", gmm_weighting = "identity",
          bstrap = FALSE, cband = FALSE
        )
      }
      ft0 <- run_tv_fit(tv_df); ftp <- run_tv_fit(tv_p); ftm <- run_tv_fit(tv_m)
      av_t <- which(is.finite(ft0$att) & ft0$chained$target_meta$role != "reference")
      fd_t <- (ftp$att[av_t] - ftm$att[av_t]) / (2.0 * eps_tv)
      an_t <- as.double(crossprod(eta_tv, as.matrix(ft0$inffunc)[, av_t, drop = FALSE])) / as.double(N)
      expect_lt(max(abs(fd_t - an_t)), 1e-5)

      ag0 <- cdid::aggte(ft0, type = "dynamic", na.rm = TRUE, bstrap = FALSE, cband = FALSE)
      agp <- cdid::aggte(ftp, type = "dynamic", na.rm = TRUE, bstrap = FALSE, cband = FALSE)
      agm <- cdid::aggte(ftm, type = "dynamic", na.rm = TRUE, bstrap = FALSE, cband = FALSE)
      fd_ag <- (agp$overall.att - agm$overall.att) / (2.0 * eps_tv)
      an_ag <- sum(eta_tv * as.double(ag0$inf.function$dynamic.inf.func)) / as.double(N)
      expect_lt(abs(fd_ag - an_ag), 1e-5)
    }
  }

  # 5. Cluster-robust parallel multiplier bootstrap, aggte(all 4 types), tidy/glance/ggdid/chain_report
  expect_error(
    cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g",
      xformla = ~ x1 + x2, data = panel_df, frame = frame_df,
      weightsname = "w", clustervars = "clust",
      control_group = "notyettreated", est_method = "dr",
      selection = "conditional_mtar", selection_est_method = "dr",
      cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "optimal",
      bstrap = TRUE, cband = TRUE, biters = 99L, pl = TRUE, cores = 2L
    ),
    class = "cdid_chained_request_parallel_error"
  )

  fit_par <- cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g",
    xformla = ~ x1 + x2, data = panel_df, frame = frame_df,
    weightsname = "w", clustervars = "clust",
    control_group = "notyettreated", est_method = "dr", outcome_pooling = "pooled",
    selection = "conditional_mtar", selection_est_method = "dr",
    imputation_pooling = "within_cohort", selection_pooling = "within_cohort",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "optimal",
    bstrap = TRUE, cband = TRUE, biters = 99L, pl = TRUE, cores = 2L,
    control = list(workers = 2L)
  )
  expect_true(is.finite(fit_par$c))
  for (tp in c("simple", "group", "dynamic", "calendar")) {
    ag <- cdid::aggte(fit_par, type = tp, na.rm = TRUE, bstrap = TRUE, cband = TRUE, biters = 99L)
    expect_true(is.finite(ag$overall.att))
    expect_gt(ag$overall.se, 0)
  }
  expect_identical(nrow(cdid::tidy(fit_par)), length(fit_par$att))
  expect_identical(nrow(cdid::glance(fit_par)), 1L)
  expect_s3_class(cdid::ggdid(fit_par), "ggplot")
  expect_s3_class(cdid::chain_report(fit_par), "cdid_chain_report")

  model <- function(capture) cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g",
    xformla = ~ x1 + x2, data = panel_df, frame = frame_df, weightsname = "w",
    clustervars = "clust", control_group = "notyettreated", est_method = "dr",
    outcome_pooling = "pooled", selection = "conditional_mtar", selection_est_method = "dr",
    selection_model = "independent_period_logit", selection_xformla = ~v_sel,
    imputation_pooling = "within_cohort", selection_pooling = "within_cohort",
    cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "identity",
    bstrap = FALSE, cband = FALSE, control = list(keep_nuisance = capture)
  )
  with_reports <- model(TRUE)
  without_reports <- model(FALSE)
  for (field in c("att", "se", "inffunc", "V_analytical")) {
    expect_identical(with_reports[[field]], without_reports[[field]])
  }
  report <- cdid::nuisance_models(with_reports)
  expect_identical(report$components$status, rep("available", 4L))
  for (role in report$components$component) {
    block <- cdid::nuisance_models(with_reports, role, block = 1)
    expect_true(nzchar(block$equation))
    expect_true(nrow(block$coefficients) > 0L)
    expect_true(all(is.finite(block$coefficients[, "Estimate"])))
    expect_true(any(is.finite(block$coefficients[, "Std. Error"])))
  }
})

test_that("cross_pair conditional_mtar handles anticipation, PT-all, partial support, and compute_inffunc = FALSE identically", {
  set.seed(20261115)
  G_ed <- c(rep(4L, 35L), rep(5L, 30L), rep(6L, 30L), rep(7L, 30L), rep(0L, 25L), rep(8L, 20L))
  N_ed <- length(G_ed)
  x1_ed <- stats::rnorm(N_ed, mean = 0.2 * (G_ed == 4L) - 0.2 * (G_ed == 5L))
  v1_ed <- stats::rnorm(N_ed)
  cl_ed <- rep(seq_len(20L), length.out = N_ed)
  frame_ed <- data.frame(id = seq_len(N_ed), G = G_ed, x1 = x1_ed, v1 = v1_ed, cl = cl_ed)

  rows_ed <- list()
  for (i in seq_len(N_ed)) {
    gi <- G_ed[i]
    pos_in_g <- sum(G_ed[seq_len(i)] == gi)
    if (pos_in_g <= 4L) next
    ts <- if (gi != 8L && pos_in_g <= 6L) {
      1L
    } else if (gi == 4L) {
      c(1L, 2L, 3L, 5L)
    } else if (gi == 5L) {
      c(2L, 3L, 4L)
    } else if (gi == 6L) {
      c(3L, 4L, 5L)
    } else if (gi == 7L) {
      c(4L, 5L, 6L)
    } else if (gi == 0L) {
      c(1L, 2L, 3L)
    } else {
      c(8L, 9L)
    }
    yi <- 0.5 * gi + (0.4 + 0.5 * x1_ed[i]) * ts +
      ifelse(gi > 0L & ts >= gi - 1L, 2.5, 0.0) + stats::rnorm(length(ts), sd = 0.35)
    rows_ed[[length(rows_ed) + 1L]] <- data.frame(id = i, period = ts, G = gi, y = yi)
  }
  panel_ed <- do.call(rbind, rows_ed)

  # 1. anticipation = 1 + partial support + conditional_mtar across all selection regimes
  for (sm in c("none", "ipw", "imputation", "dr")) {
    sel_arg <- if (sm == "none") "none" else "conditional_mtar"
    sem_arg <- if (sm == "none") "ipw" else sm
    mp_ant <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1, selection_xformla = if (sm %in% c("ipw", "dr")) ~ x1 + v1 else NULL,
      data = panel_ed, frame = frame_ed,
      control_group = "notyettreated", anticipation = 1L,
      cross_pair = TRUE, pt_assumption = "ctrl",
      est_method = "dr", selection = sel_arg, selection_est_method = sem_arg,
      clustervars = "cl", bstrap = FALSE, cband = FALSE
    )
    idx_4_5 <- which(mp_ant$group == 4 & mp_ant$t == 5)
    expect_length(idx_4_5, 1L)
    expect_true(is.finite(mp_ant$att[idx_4_5]))
    expect_true(is.finite(mp_ant$se[idx_4_5]))
    expect_true(all(is.na(mp_ant$att[mp_ant$group == 8])))

    ag_d <- cdid::aggte(mp_ant, type = "dynamic", na.rm = TRUE, min_e = -1, max_e = 1, bstrap = FALSE)
    expect_true(is.finite(ag_d$overall.att))
    expect_gt(ag_d$overall.se, 0)
  }

  # 2. pt_assumption = "all" + cross_pair = TRUE + conditional_mtar
  for (sm in c("ipw", "imputation", "dr")) {
    mp_all <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1, selection_xformla = if (sm %in% c("ipw", "dr")) ~ x1 + v1 else NULL,
      data = panel_ed, frame = frame_ed,
      control_group = "notyettreated", cross_pair = TRUE,
      pt_assumption = "all", pt_start = 1,
      est_method = "dr", selection = "conditional_mtar", selection_est_method = sm,
      bstrap = FALSE, cband = FALSE
    )
    idx_4_5 <- which(mp_all$group == 4 & mp_all$t == 5)
    expect_true(is.finite(mp_all$att[idx_4_5]))
    expect_gt(mp_all$se[idx_4_5], 0)
  }

  # 3. compute_inffunc = FALSE exact point-estimate parity across weightings & selection regimes
  for (wt in c("optimal", "within", "identity")) {
    for (sm in c("none", "ipw", "imputation", "dr")) {
      sel_arg <- if (sm == "none") "none" else "conditional_mtar"
      sem_arg <- if (sm == "none") "ipw" else sm
      m_t <- cdid::att_gt(
        yname = "y", tname = "period", idname = "id", gname = "G",
        xformla = ~ x1, data = panel_ed, frame = frame_ed,
        control_group = "notyettreated", cross_pair = TRUE, pt_assumption = "ctrl",
        gmm_weighting = wt, est_method = "dr",
        selection = sel_arg, selection_est_method = sem_arg,
        compute_inffunc = TRUE, bstrap = FALSE, cband = FALSE
      )
      m_f <- cdid::att_gt(
        yname = "y", tname = "period", idname = "id", gname = "G",
        xformla = ~ x1, data = panel_ed, frame = frame_ed,
        control_group = "notyettreated", cross_pair = TRUE, pt_assumption = "ctrl",
        gmm_weighting = wt, est_method = "dr",
        selection = sel_arg, selection_est_method = sem_arg,
        compute_inffunc = FALSE
      )
      expect_identical(is.na(m_t$att), is.na(m_f$att))
      expect_lt(max(abs(m_t$att - m_f$att), na.rm = TRUE), 1e-14)
      expect_null(m_f$inffunc)
      expect_true(all(is.na(m_f$se)))
    }
  }
})


test_that("within weighting supports fix_weights = base_period across estimators, pooling, and conditional_mtar regimes", {
  ns <- asNamespace("cdid")

  # 1a. Balanced nevertreated with static covariates and time-varying weights:
  #     exact all-cell parity between cross_pair = FALSE and cross_pair = TRUE across separate & pooled
  set.seed(20261125)
  rows_bp_stat <- list()
  rows_bp_tvx <- list()
  uid_bp <- 0L
  for (gv in c(3L, 4L, 0L)) {
    for (u in seq_len(28L)) {
      uid_bp <- uid_bp + 1L
      x_base <- stats::rnorm(1L)
      w_base <- stats::runif(1L, 0.7, 1.8)
      a_u <- stats::rnorm(1L, sd = 0.4)
      for (tv in 1:4) {
        xt <- x_base + 0.12 * tv + stats::rnorm(1L, sd = 0.08)
        wt <- w_base * (1 + 0.05 * tv * ((uid_bp %% 3L) + 1L))
        yt <- a_u + 0.35 * tv + 0.4 * xt +
          (if (gv > 0L && tv >= gv) 1.75 else 0.0) + stats::rnorm(1L, sd = 0.3)
        rows_bp_stat[[length(rows_bp_stat) + 1L]] <- data.frame(
          id = uid_bp, period = tv, G = gv, x1 = x_base, w = wt, y = yt
        )
        rows_bp_tvx[[length(rows_bp_tvx) + 1L]] <- data.frame(
          id = uid_bp, period = tv, G = gv, x1 = xt, w = wt, y = yt
        )
      }
    }
  }
  df_bp_stat <- do.call(rbind, rows_bp_stat)
  df_bp_tvx <- do.call(rbind, rows_bp_tvx)

  for (em in c("ipw", "reg", "dr")) {
    op_vec <- if (em == "ipw") "separate" else c("separate", "pooled")
    for (op in op_vec) {
      f_sp <- cdid::att_gt(
        yname = "y", tname = "period", idname = "id", gname = "G",
        xformla = ~ x1, weightsname = "w", fix_weights = "base_period",
        data = df_bp_stat, control_group = "nevertreated",
        est_method = em, outcome_pooling = op,
        cross_pair = FALSE, gmm_weighting = "within",
        bstrap = FALSE, cband = FALSE
      )
      f_cp <- cdid::att_gt(
        yname = "y", tname = "period", idname = "id", gname = "G",
        xformla = ~ x1, weightsname = "w", fix_weights = "base_period",
        data = df_bp_stat, control_group = "nevertreated",
        est_method = em, outcome_pooling = op,
        cross_pair = TRUE, gmm_weighting = "within",
        bstrap = FALSE, cband = FALSE
      )
      expect_lt(max(abs(f_sp$att - f_cp$att)), 1e-11)
      expect_lt(max(abs(f_sp$se - f_cp$se)), 1e-11)
      expect_lt(max(abs(as.matrix(f_sp$inffunc) - as.matrix(f_cp$inffunc))), 1e-10)
    }
  }

  # 1b. With time-varying covariates + time-varying weights, separate estimators match
  #     identically on all post-treatment cells
  for (em in c("ipw", "reg", "dr")) {
    f_sp_tv <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1, weightsname = "w", fix_weights = "base_period",
      data = df_bp_tvx, control_group = "nevertreated",
      est_method = em, outcome_pooling = "separate",
      cross_pair = FALSE, gmm_weighting = "within",
      bstrap = FALSE, cband = FALSE
    )
    f_cp_tv <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1, weightsname = "w", fix_weights = "base_period",
      data = df_bp_tvx, control_group = "nevertreated",
      est_method = em, outcome_pooling = "separate",
      cross_pair = TRUE, gmm_weighting = "within",
      bstrap = FALSE, cband = FALSE
    )
    post_idx <- which(f_sp_tv$t >= f_sp_tv$group)
    expect_lt(max(abs(f_sp_tv$att[post_idx] - f_cp_tv$att[post_idx])), 1e-11)
    expect_lt(max(abs(f_sp_tv$se[post_idx] - f_cp_tv$se[post_idx])), 1e-11)
    expect_lt(max(abs(as.matrix(f_sp_tv$inffunc[, post_idx]) - as.matrix(f_cp_tv$inffunc[, post_idx]))), 1e-10)
  }

  # 2. Unbalanced notyettreated with F \ B and B \ A units under fix_weights = base_period + gmm_weighting = within
  set.seed(20261126)
  cohs_u <- c(3L, 4L, 0L)
  n_per_u <- 28L
  N_u <- length(cohs_u) * n_per_u
  g_u <- rep(cohs_u, each = n_per_u)
  x0_u <- stats::rnorm(N_u, sd = 0.7)
  v0_u <- stats::rnorm(N_u, sd = 0.6)
  w0_u <- stats::runif(N_u, 0.7, 1.8)
  frame_u <- data.frame(id = seq_len(N_u), G = g_u, x1 = x0_u, v1 = v0_u)

  rows_u <- list()
  for (ci in seq_along(cohs_u)) {
    idx <- which(g_u == cohs_u[ci])
    for (i in idx) {
      s_row <- if (i %in% idx[1:3]) {
        rep(FALSE, 4L)
      } else if (i %in% idx[4:6]) {
        c(FALSE, TRUE, FALSE, FALSE)
      } else {
        sr <- stats::runif(4L) < 0.78
        sr[c(2L, 3L)] <- TRUE
        sr
      }
      for (tv in 1:4) {
        wt <- w0_u[i] * (1 + 0.04 * tv * ((i %% 3L) + 1L))
        yt <- if (s_row[tv]) {
          0.35 * tv + 0.45 * x0_u[i] +
            (if (g_u[i] > 0L && tv >= g_u[i]) 1.5 else 0.0) + stats::rnorm(1L, sd = 0.3)
        } else {
          NA_real_
        }
        rows_u[[length(rows_u) + 1L]] <- data.frame(
          id = i, period = tv, G = g_u[i], x1 = x0_u[i], v1 = v0_u[i], w = wt, y = yt
        )
      }
    }
  }
  panel_u <- do.call(rbind, rows_u)
  eta_u <- stats::rnorm(N_u); eta_u <- eta_u - mean(eta_u)
  eps_u <- 1e-6
  panel_u_p <- panel_u; panel_u_p$w <- panel_u$w * (1.0 + eps_u * eta_u[panel_u$id])
  panel_u_m <- panel_u; panel_u_m$w <- panel_u$w * (1.0 - eps_u * eta_u[panel_u$id])

  for (sm in c("none", "ipw", "imputation", "dr")) {
    sel_arg <- if (sm == "none") "none" else "conditional_mtar"
    sem_arg <- if (sm == "none") "ipw" else sm
    for (cp in c(FALSE, TRUE)) {
      run_u_fit <- function(d_in, wt) {
        cdid::att_gt(
          yname = "y", tname = "period", idname = "id", gname = "G",
          xformla = ~ x1, selection_xformla = if (sm %in% c("ipw", "dr")) ~ x1 + v1 else NULL,
          weightsname = "w", fix_weights = "base_period",
          data = d_in, frame = frame_u,
          control_group = "notyettreated", cross_pair = cp, pt_assumption = "post",
          est_method = "dr", selection = sel_arg, selection_est_method = sem_arg,
          gmm_weighting = wt,
          bstrap = FALSE, cband = FALSE
        )
      }
      fit_pub <- run_u_fit(panel_u, "within")
      fin_idx <- which(is.finite(fit_pub$att) & fit_pub$chained$target_meta$role != "reference")
      expect_gt(length(fin_idx), 0L)
      expect_true(all(is.finite(fit_pub$se[fin_idx]) & fit_pub$se[fin_idx] > 0))
      expect_equal(nrow(fit_pub$inffunc), N_u)
      expect_lt(max(abs(colMeans(as.matrix(fit_pub$inffunc[, fin_idx, drop = FALSE])))), 1e-12)
      ag_w <- cdid::aggte(fit_pub, type = "dynamic", na.rm = TRUE, bstrap = FALSE, cband = FALSE)
      expect_true(is.finite(ag_w$overall.att))
      expect_gt(ag_w$overall.se, 0)

      fit_id <- run_u_fit(panel_u, "identity")
      fit_p <- run_u_fit(panel_u_p, "identity")
      fit_m <- run_u_fit(panel_u_m, "identity")
      fd_u <- (fit_p$att[fin_idx] - fit_m$att[fin_idx]) / (2.0 * eps_u)
      an_u <- as.double(crossprod(eta_u, as.matrix(fit_id$inffunc)[, fin_idx, drop = FALSE])) / as.double(N_u)
      expect_lt(max(abs(fd_u - an_u)), 1e-5)
    }
  }

  # Explicit control$chunk_rows agrees with automatic chunk sizing on pooled selection + pooled OR
  for (sm_ch in c("ipw", "dr")) {
    fit_auto_ch <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1, selection_xformla = ~ x1 + v1,
      data = panel_u, frame = frame_u,
      control_group = "notyettreated", cross_pair = TRUE, pt_assumption = "ctrl",
      est_method = "dr", outcome_pooling = "pooled",
      selection = "conditional_mtar", selection_est_method = sm_ch,
      selection_pooling = "within_cohort", imputation_pooling = "within_cohort",
      bstrap = FALSE, cband = FALSE
    )
    fit_ch16 <- cdid::att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x1, selection_xformla = ~ x1 + v1,
      data = panel_u, frame = frame_u,
      control_group = "notyettreated", cross_pair = TRUE, pt_assumption = "ctrl",
      est_method = "dr", outcome_pooling = "pooled",
      selection = "conditional_mtar", selection_est_method = sm_ch,
      selection_pooling = "within_cohort", imputation_pooling = "within_cohort",
      control = list(chunk_rows = 16L),
      bstrap = FALSE, cband = FALSE
    )
    expect_true(isTRUE(fit_auto_ch$chained$nuisance_diagnostics$accounting$automatic_chunk_rows))
    expect_false(isTRUE(fit_ch16$chained$nuisance_diagnostics$accounting$automatic_chunk_rows))
    expect_lt(max(abs(fit_auto_ch$att - fit_ch16$att), na.rm = TRUE), 1e-11)
    expect_lt(max(abs(fit_auto_ch$se - fit_ch16$se), na.rm = TRUE), 1e-11)
  }
})

test_that("within-unit covariance preserves exact n_A scale invariance under fix_weights = base_period", {
  set.seed(20261112L)
  N <- 120L
  TT <- 6L
  cohs <- rep(c(4, 5, 6, 0), each = N / 4L)
  x0 <- stats::rnorm(N)
  df_full <- do.call(rbind, lapply(seq_len(TT), function(tt) {
    data.frame(
      id = seq_len(N), t = tt, g = cohs, x0 = x0,
      w_tv = ifelse(tt == TT, 2.5, 2.0),
      y = stats::rnorm(N) + ifelse(cohs > 0 & tt >= cohs, 1.5, 0.0) + 0.3 * stats::rnorm(N)
    )
  }))
  df_sub <- df_full
  df_sub$y[sample.int(nrow(df_full), round(0.12 * nrow(df_full)))] <- NA_real_
  frame_u <- df_full[df_full$t == 1L, c("id", "g", "x0")]

  ns <- asNamespace("cdid")
  run_raw_cp_v <- function(wname, fw_mode, sm_mode, op_mode = "separate", w_scale = 1.0) {
    d_in <- df_sub
    if (!is.null(wname) && w_scale != 1.0) d_in[[wname]] <- d_in[[wname]] * w_scale
    sel_arg <- if (sm_mode == "none") "none" else "conditional_mtar"
    sem_arg <- if (sm_mode == "none") "ipw" else sm_mode
    fr_arg <- if (sm_mode == "none") NULL else frame_u
    prep <- ns$chained_request_prepare_inputs(
      data = d_in, frame = fr_arg, yname = "y", tname = "t", idname = "id", gname = "g",
      xformla = ~ x0, weightsname = wname, fix_weights = fw_mode,
      doubly_robust = (sm_mode %in% c("none", "ipw")), outcome_pooling = op_mode,
      selection = sel_arg, selection_est_method = sem_arg,
      imputation_pooling = if (sm_mode %in% c("imputation", "dr")) "separate" else NULL,
      imputation_est_method = if (sm_mode == "imputation") "dr" else NULL,
      augmented_est_method = if (sm_mode == "dr") "dr" else NULL,
      cross_pair = TRUE
    )
    dims <- ns$chained_public_dimensions(prep, biters = 10L, clustered = FALSE)
    res_c <- ns$chained_public_resolve_control(NULL, gmm_weighting = "within", gmm_covariance = "tikhonov")
    prod_c <- ns$chained_public_pipeline_control(dims, res_c, gmm_weighting = "within", gmm_covariance = "tikhonov", workers = 1L)
    fit <- ns$run_chained_cross_pair_fit(
      data = prep$data_internal, idname = "id", tname = "t", yname = "y", gname = "g",
      control = prod_c$control, xformla = prep$xformla, cluster_vec = prep$cluster_full,
      target_cohorts = prep$target_cohorts, anticipation = 0L,
      control_group = "notyettreated", base_period = "varying",
      weighting = "within", pt_assumption = "all", pt_start = 1,
      est_method = "dr", cband = FALSE, alp = 0.05, biters = 10L,
      cp04 = prep$cp04, frame_id = prep$frame_ids, frame_Z = prep$frame_Z,
      frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
      allow_frame_only = prep$n_frame_only > 0L, bstrap = FALSE,
      gmm_covariance = "tikhonov", prepared = prep$prepared
    )
    fit$covariance_stream$V
  }

  for (cp_mode in c(FALSE, TRUE)) {
    for (sm_mode in c("none", "ipw", "imputation", "dr")) {
      sel_arg <- if (sm_mode == "none") "none" else "conditional_mtar"
      sem_arg <- if (sm_mode == "none") "ipw" else sm_mode
      fr_arg <- if (sm_mode == "none") NULL else frame_u

      f_unw <- cdid::att_gt(
        yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x0,
        data = df_sub, frame = fr_arg, weightsname = NULL,
        control_group = "notyettreated", allow_unbalanced_panel = TRUE,
        pt_assumption = "all", pt_start = 1, gmm_weighting = "within",
        cross_pair = cp_mode, selection = sel_arg, selection_est_method = sem_arg,
        bstrap = FALSE, cband = FALSE
      )
      f_fp <- cdid::att_gt(
        yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x0,
        data = df_sub, frame = fr_arg, weightsname = "w_tv", fix_weights = "first_period",
        control_group = "notyettreated", allow_unbalanced_panel = TRUE,
        pt_assumption = "all", pt_start = 1, gmm_weighting = "within",
        cross_pair = cp_mode, selection = sel_arg, selection_est_method = sem_arg,
        bstrap = FALSE, cband = FALSE
      )
      f_bp <- cdid::att_gt(
        yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x0,
        data = df_sub, frame = fr_arg, weightsname = "w_tv", fix_weights = "base_period",
        control_group = "notyettreated", allow_unbalanced_panel = TRUE,
        pt_assumption = "all", pt_start = 1, gmm_weighting = "within",
        cross_pair = cp_mode, selection = sel_arg, selection_est_method = sem_arg,
        bstrap = FALSE, cband = FALSE
      )

      expect_lt(max(abs(f_unw$att - f_fp$att), na.rm = TRUE), 1e-12)
      expect_lt(max(abs(f_fp$att - f_bp$att), na.rm = TRUE), 1e-12)
      expect_lt(max(abs(f_fp$se - f_bp$se), na.rm = TRUE), 1e-12)
      expect_lt(max(abs(as.matrix(f_fp$inffunc) - as.matrix(f_bp$inffunc)), na.rm = TRUE), 1e-11)

      if (cp_mode) {
        V_unw <- run_raw_cp_v(NULL, NULL, sm_mode, "separate")
        V_fp <- run_raw_cp_v("w_tv", "first_period", sm_mode, "separate")
        V_bp <- run_raw_cp_v("w_tv", "base_period", sm_mode, "separate")
        expect_lt(max(abs(V_unw - V_fp)), 1e-12)
        expect_lt(max(abs(V_fp - V_bp)), 1e-12)

        V_bp_p1 <- run_raw_cp_v("w_tv", "base_period", sm_mode, "pooled", 1.0)
        V_bp_p7 <- run_raw_cp_v("w_tv", "base_period", sm_mode, "pooled", 7.5)
        expect_lt(max(abs(V_bp_p1 - V_bp_p7)), 1e-12)
      }
    }
  }
})

test_that("causal anticipation cells use caller pt_assumption and chain_report covers pre and anticipation cells", {
  # Panel: t in 1..5, cohorts g = 4 (target, anticipation = 1 => ref b_4 = 2, causal anticipation t = 3),
  # g = 5, g = 0. Pre-reference shock on control cohort g = 5 at t = 1 ONLY (before b_4 = 2).
  N_per <- 30L
  cohs <- rep(c(4, 5, 0), each = N_per)
  N <- length(cohs)
  df_ant <- do.call(rbind, lapply(1:5, function(tt) {
    # True causal effects for g = 4: tau(4,3) = 1.75 (anticipation), tau(4,4) = 3.0, tau(4,5) = 4.25
    # True causal effect for g = 5: tau(5,4) = 0.8 (anticipation), tau(5,5) = 2.0
    eff <- ifelse(cohs == 4 & tt == 3, 1.75,
           ifelse(cohs == 4 & tt == 4, 3.00,
           ifelse(cohs == 4 & tt == 5, 4.25,
           ifelse(cohs == 5 & tt == 4, 0.80,
           ifelse(cohs == 5 & tt == 5, 2.00, 0.0)))))
    # Pre-reference shock on g = 5 at t = 1 (before b_4 = 2): must NOT bias g = 4 under PT-post or PT-ctrl(pt_start = 2)
    shock <- ifelse(cohs == 5 & tt == 1, 9.0, 0.0)
    # Unit fixed effect + common time trend + tiny noise for non-degenerate covariance
    u_fe <- (seq_len(N) %% 7L) * 0.2
    data.frame(
      id = seq_len(N), t = tt, g = cohs,
      y = u_fe + 0.5 * tt + eff + shock
    )
  }))

  for (bp in c("varying", "universal")) {
    fit_post <- cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g",
      data = df_ant, control_group = "notyettreated", anticipation = 1,
      base_period = bp, pt_assumption = "post", cross_pair = TRUE,
      gmm_weighting = "within", bstrap = FALSE, cband = FALSE
    )
    idx_43 <- which(fit_post$group == 4 & fit_post$t == 3)
    idx_44 <- which(fit_post$group == 4 & fit_post$t == 4)
    expect_equal(length(idx_43), 1L)
    expect_lt(abs(fit_post$att[idx_43] - 1.75), 1e-10)
    expect_lt(abs(fit_post$att[idx_44] - 3.00), 1e-10)

    fit_ctrl2 <- cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g",
      data = df_ant, control_group = "notyettreated", anticipation = 1,
      base_period = bp, pt_assumption = "ctrl", pt_start = 2, cross_pair = TRUE,
      gmm_weighting = "within", bstrap = FALSE, cband = FALSE
    )
    expect_lt(abs(fit_ctrl2$att[which(fit_ctrl2$group == 4 & fit_ctrl2$t == 3)] - 1.75), 1e-10)
    expect_lt(abs(fit_ctrl2$att[which(fit_ctrl2$group == 4 & fit_ctrl2$t == 4)] - 3.00), 1e-10)

    rep_ant <- cdid::chain_report(fit_ctrl2, group = 4, t = 3)
    expect_s3_class(rep_ant, "cdid_chain_report")
    expect_equal(nrow(rep_ant$summary), 1L)
    expect_true(isTRUE(rep_ant$summary$available[1L]))
    expect_gt(nrow(as.data.frame(rep_ant)), 0L)
  }
})

test_that("exact numerical Jacobian oracle matches V_within for dated covariates + non-uniform weights and pooled OR + varying weights", {
  set.seed(20261113)
  N <- 36L
  TT <- 4L
  cohs <- rep(c(3, 4, 0), each = N / 3L)
  w_static <- rep(c(1.0, 2.5, 0.7, 3.2), length.out = N)
  x_base <- stats::rnorm(N)
  df <- do.call(rbind, lapply(seq_len(TT), function(tt) {
    data.frame(
      id = seq_len(N), t = tt, g = cohs,
      x_tv = x_base + 0.35 * tt * ((seq_len(N) %% 3L) - 1),
      w_st = w_static,
      w_var = w_static * (1 + 0.4 * tt),
      y = stats::rnorm(N) + ifelse(cohs > 0 & tt >= cohs, 1.8, 0.0)
    )
  }))
  # Introduce mild unbalancedness so repeat-observer subset is non-trivial
  df$y[c(2L, 17L, 29L)] <- NA_real_

  ns <- asNamespace("cdid")
  check_fd_within_oracle <- function(d_in, xform, wname, fw_mode, op_mode, est_m) {
    run_once <- function(d_cur) {
      prep <- ns$chained_request_prepare_inputs(
        data = d_cur, frame = NULL, yname = "y", tname = "t", idname = "id", gname = "g",
        xformla = xform, weightsname = wname, fix_weights = fw_mode,
        doubly_robust = identical(est_m, "dr"), outcome_pooling = op_mode,
        outcome_regression = identical(est_m, "reg"),
        selection = "none", selection_est_method = "ipw", cross_pair = TRUE
      )
      dims <- ns$chained_public_dimensions(prep, biters = 10L, clustered = FALSE)
      res_c <- ns$chained_public_resolve_control(NULL, gmm_weighting = "within", gmm_covariance = "tikhonov")
      prod_c <- ns$chained_public_pipeline_control(dims, res_c, gmm_weighting = "within", gmm_covariance = "tikhonov", workers = 1L)
      fit <- ns$run_chained_cross_pair_fit(
        data = prep$data_internal, idname = "id", tname = "t", yname = "y", gname = "g",
        control = prod_c$control, xformla = prep$xformla, cluster_vec = prep$cluster_full,
        target_cohorts = prep$target_cohorts, anticipation = 0L,
        control_group = "notyettreated", base_period = "varying",
        weighting = "within", pt_assumption = "post", pt_start = NULL,
        est_method = est_m, cband = FALSE, alp = 0.05, biters = 10L,
        cp04 = prep$cp04, frame_id = prep$frame_ids, frame_Z = prep$frame_Z,
        frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
        allow_frame_only = FALSE, bstrap = FALSE,
        gmm_covariance = "tikhonov", prepared = prep$prepared
      )
      list(fit = fit, prep = prep)
    }

    base_out <- run_once(d_in)
    m0 <- base_out$fit$covariance_stream$mhat
    V_code <- base_out$fit$covariance_stream$V
    links_used <- base_out$fit$design_used$links
    L <- length(m0)

    # Map rows of d_in on repeat-observer units A
    obs_counts <- tapply(!is.na(d_in$y), d_in$id, sum)
    a_ids <- as.integer(names(obs_counts)[obs_counts >= 2L])
    n_a <- length(a_ids)

    # Extract per-slice weights on A
    w_by_g <- lapply(sort(unique(links_used$g)), function(gg) {
      b_g <- gg - 1L
      t_w <- if (identical(fw_mode, "base_period")) b_g else 1L
      sub_w <- d_in[d_in$t == t_w, c("id", wname)]
      sub_w[[wname]][match(a_ids, sub_w$id)]
    })
    names(w_by_g) <- as.character(sort(unique(links_used$g)))

    # With fixed support, covariates and weights, pre-GMM moments are affine
    # in the outcomes. Use a dyadic step to reduce cancellation in the oracle.
    eps <- 2^-13
    U_num <- array(0.0, dim = c(n_a, TT, L))
    for (ii in seq_len(n_a)) {
      uid <- a_ids[ii]
      for (tt in seq_len(TT)) {
        r_idx <- which(d_in$id == uid & d_in$t == tt)
        if (!length(r_idx) || is.na(d_in$y[r_idx])) next
        d_plus <- d_in; d_plus$y[r_idx] <- d_in$y[r_idx] + eps
        d_minus <- d_in; d_minus$y[r_idx] <- d_in$y[r_idx] - eps
        plus <- run_once(d_plus)$fit
        minus <- run_once(d_minus)$fit
        stopifnot(identical(plus$design_used$links, links_used),
                  identical(minus$design_used$links, links_used))
        mp <- plus$covariance_stream$mhat
        mm <- minus$covariance_stream$mhat
        U_num[ii, tt, ] <- (mp - mm) / (2 * eps)
      }
    }

    V_oracle <- matrix(0.0, L, L)
    for (g_chr in names(w_by_g)) {
      gg <- as.numeric(g_chr)
      idx_g <- which(links_used$g == gg)
      w_g <- w_by_g[[g_chr]]
      pos <- w_g > 0
      # Active A_r units: under dated covariates, units with finite x_tv at b_g; here all A have finite x_tv
      scale_i <- numeric(n_a)
      scale_i[pos] <- (sum(pos) * mean(w_g[pos])) / w_g[pos]
      for (tt in seq_len(TT)) {
        U_sub <- U_num[, tt, idx_g, drop = FALSE]
        dim(U_sub) <- c(n_a, length(idx_g))
        V_oracle[idx_g, idx_g] <- V_oracle[idx_g, idx_g] + crossprod(U_sub * sqrt(scale_i))
      }
    }
    expect_lt(max(abs(V_code - V_oracle)), 1e-8)
  }

  # Case (a): Time-varying covariates + static non-uniform weights (DR, separate)
  check_fd_within_oracle(df, ~ x_tv, "w_st", "first_period", "separate", "dr")
  # Case (b): Pooled OR + base_period varying weights across multiple slices
  check_fd_within_oracle(df, ~ x_tv, "w_var", "base_period", "pooled", "reg")
})

test_that("PT-ctrl populates pretest metadata and fast cross-pair fits expose nuisance diagnostics", {
  set.seed(909)
  N <- 120L
  cohs <- rep(c(4, 5, 0), each = N / 3L)
  x0 <- stats::rnorm(N)
  df <- do.call(rbind, lapply(1:5, function(tt) {
    data.frame(
      id = seq_len(N), t = tt, g = cohs, x0 = x0,
      y = stats::rnorm(N) + 0.3 * x0 * tt + ifelse(cohs > 0 & tt >= cohs, 1.5, 0.0)
    )
  }))

  for (cp_flag in c(FALSE, TRUE)) {
    fit_ctrl <- cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x0,
      data = df, control_group = "notyettreated",
      pt_assumption = "ctrl", pt_start = 2, cross_pair = cp_flag,
      bstrap = FALSE, cband = FALSE
    )
    expect_true(!is.null(fit_ctrl$chained$pretest$interpretation))
    expect_s3_class(fit_ctrl$chained$pretest$window_tags, "data.frame")
    expect_gt(nrow(fit_ctrl$chained$pretest$window_tags), 0L)
    expect_true(all(c("pre_pt_start", "in_pt_window") %in% fit_ctrl$chained$pretest$window_tags$pt_window_tag))
    expect_equal(fit_ctrl$chained$nuisance_diagnostics$status, "available")
  }

  fit_fast_nocov <- cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g",
    data = df, control_group = "notyettreated",
    pt_assumption = "post", cross_pair = TRUE, est_method = "ipw",
    bstrap = FALSE, cband = FALSE
  )
  expect_equal(fit_fast_nocov$chained$nuisance_diagnostics$status, "available")
})

test_that("compute_inffunc = FALSE preserves exact point estimates across cross-pair weighting and selection modes", {
  set.seed(404)
  N <- 120L
  cohs <- rep(c(3, 4, 0), each = N / 3L)
  x0 <- stats::rnorm(N)
  df_full <- do.call(rbind, lapply(1:4, function(tt) {
    data.frame(
      id = seq_len(N), t = tt, g = cohs, x0 = x0,
      y = stats::rnorm(N) + 0.25 * x0 * tt + ifelse(cohs > 0 & tt >= cohs, 2.0, 0.0)
    )
  }))
  df_sub <- df_full
  df_sub$y[sample.int(nrow(df_full), round(0.15 * nrow(df_full)))] <- NA_real_
  fr <- df_full[df_full$t == 1L, c("id", "g", "x0")]

  for (wt in c("within", "identity", "optimal")) {
    f_if <- cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x0,
      data = df_sub, control_group = "notyettreated", allow_unbalanced_panel = TRUE,
      cross_pair = TRUE, gmm_weighting = wt, est_method = "dr",
      bstrap = FALSE, cband = FALSE, compute_inffunc = TRUE
    )
    f_noif <- cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x0,
      data = df_sub, control_group = "notyettreated", allow_unbalanced_panel = TRUE,
      cross_pair = TRUE, gmm_weighting = wt, est_method = "dr",
      bstrap = FALSE, cband = FALSE, compute_inffunc = FALSE
    )
    expect_null(f_noif$inffunc)
    expect_lt(max(abs(f_if$att - f_noif$att), na.rm = TRUE), 1e-12)
  }

  for (sem in c("ipw", "imputation", "dr")) {
    f_if <- cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x0,
      data = df_sub, frame = fr, control_group = "notyettreated", allow_unbalanced_panel = TRUE,
      cross_pair = TRUE, gmm_weighting = "within", est_method = "dr",
      selection = "conditional_mtar", selection_est_method = sem,
      bstrap = FALSE, cband = FALSE, compute_inffunc = TRUE
    )
    f_noif <- cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x0,
      data = df_sub, frame = fr, control_group = "notyettreated", allow_unbalanced_panel = TRUE,
      cross_pair = TRUE, gmm_weighting = "within", est_method = "dr",
      selection = "conditional_mtar", selection_est_method = sem,
      bstrap = FALSE, cband = FALSE, compute_inffunc = FALSE
    )
    expect_null(f_noif$inffunc)
    expect_lt(max(abs(f_if$att - f_noif$att), na.rm = TRUE), 1e-12)
  }
})

test_that("chain_report, as.data.frame, and print enforce typed cdid_chain_report_input_error validation", {
  df <- make_short_window_panel(n_per_cohort = rep(15L, 6L), seed = 77)
  fit <- cdid::att_gt(
    yname = "y", tname = "period", idname = "id", gname = "G",
    data = df, control_group = "notyettreated", allow_unbalanced_panel = TRUE,
    cross_pair = TRUE, gmm_weighting = "within", bstrap = FALSE, cband = FALSE
  )
  rep_ok <- cdid::chain_report(fit)

  expect_error(cdid::chain_report(list()), class = "cdid_chain_report_input_error")
  expect_error(cdid::chain_report(fit, group = NaN), class = "cdid_chain_report_input_error")
  expect_error(cdid::chain_report(fit, t = c(2005, NA_real_)), class = "cdid_chain_report_input_error")
  expect_error(cdid::chain_report(fit, group = as.Date("2020-01-01")), class = "cdid_chain_report_input_error")
  expect_error(cdid::chain_report(fit, tol = -1e-6), class = "cdid_chain_report_input_error")
  expect_error(cdid::chain_report(fit, tol = Inf), class = "cdid_chain_report_input_error")
  expect_error(cdid::chain_report(fit, tol = c(1e-8, 1e-6)), class = "cdid_chain_report_input_error")

  expect_error(as.data.frame(rep_ok, what = NA_character_), class = "cdid_chain_report_input_error")
  expect_error(as.data.frame(rep_ok, what = c("summary", "witness")), class = "cdid_chain_report_input_error")
  expect_error(as.data.frame(rep_ok, what = "bogus"), class = "cdid_chain_report_input_error")
  expect_error(print(rep_ok, max_cells = -1L), class = "cdid_chain_report_input_error")
  expect_error(print(rep_ok, max_cells = Inf), class = "cdid_chain_report_input_error")
  expect_error(print(rep_ok, digits = -1L), class = "cdid_chain_report_input_error")
  expect_error(print(rep_ok, digits = NaN), class = "cdid_chain_report_input_error")
})

test_that("pooled cross-pair reporting admits its score workspace before preparing residuals", {
  n <- 32L
  b <- 40L
  x <- seq_len(n) - mean(seq_len(n))
  X <- cbind("(Intercept)" = 1, x = x)
  core <- list2env(list(
    N = n, L = b, sel_mode = "none", est_method = "reg", outcome_pooling = "pooled",
    compute_if = TRUE, need_or = TRUE, n_or_blks = b, n_trt_fits = 0L,
    rate_list = list(), ctrl = list(owner_control = list()),
    recipe_or = list(intercept = TRUE, slope_index = 2L, columns = colnames(X),
                     slope_columns = "x", intercept_index = 1L),
    st = list(cohorts = Inf), all_times = seq_len(b + 1L), link_blk_id = seq_len(b),
    blk_fits = lapply(seq_len(b), function(j) list(group_id = 1L, support_ids = j,
      keep = 1L, alpha = 0, beta = 1, L = j, mean_d = 0, weight_sum = n)),
    sup_list = lapply(seq_len(b), function(j) list(Z_obs = X, cell = list(rows_obs = seq_len(n),
      fw = rep(1, n), dy = j * x, ci = 1L, cs = 1L, ct = j + 1L)))
  ))
  sink <- cross_pair_nuisance_sink(NULL, seq_len(n), seq_len(n), list())
  requested <- 0
  testthat::with_mocked_bindings({
    cross_pair_nuisance_weighted_collect(sink, core)
  }, cross_pair_nuisance_admit = function(s, ctrl, live, scratch, work, label) {
    requested <<- scratch
    "Synthetic workspace refusal."
  }, cross_pair_nuisance_or_piece = function(...) {
    stop("Score inputs were prepared before workspace admission.")
  }, .package = "cdid")
  rec <- sink$outcome[[1L]]
  expect_gt(requested, n * b * (b + 1))
  expect_identical(rec$se$code, "budget")
  expect_true(is.na(rec$check))
  expect_identical(tail(rec$estimate, 1L), 1)
})
