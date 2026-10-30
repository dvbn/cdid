# Host tests for the private chained result adapters (R/chained_did_results.R).
#
# Requires the pinned did 2.5.1 (constructors + registered methods) and the
# cdid namespace (installed package or source harness). Fixtures below are
# small, deterministic and hand-built: no fit, DGP, RNG, bootstrap or Monte
# Carlo. Where a fixture needs a standard error, it is formed from the fixture's
# own fixed IF numbers by the pipeline's documented definition (centered cluster
# sums, crossprod / n^2, Omega = n * vcov) so the adapters' consistency checks
# have something true to check; the adapters themselves compute nothing.

# ------------------------------------------------------------------------------
# Fixture helpers
# ------------------------------------------------------------------------------

cl_center <- function(X, cl) {
  S <- rowsum(as.matrix(X), cl)
  sweep(S, 2L, colMeans(S), "-")
}
cl_se <- function(X, cl, n) sqrt(colSums(cl_center(X, cl)^2)) / n
cl_cov <- function(X, cl, n) crossprod(cl_center(X, cl)) / n^2

fx_ids <- c(101L, 102L, 103L, 104L, 105L, 106L, 107L, 108L)
fx_cluster <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 5L)       # unit 108 alone in cluster 5
fx_eligible <- c(rep(TRUE, 7L), FALSE)                # unit 108 is frame-only
fx_cohort <- c(Inf, 3, 3, 4, 4, Inf, 3, Inf)          # full-frame cohorts (Inf = never treated)
fx_n <- 8L

# Reported cells: deliberately NOT group-major, sparse, with one exact-zero
# reference cell and one unavailable cell.
fx_target_meta <- function() {
  data.frame(
    g = c(4, 3, 3, 3, 4), time = c(4, 1, 3, 4, 5),
    role = c("post", "pre", "post", "post", "post"),
    target_id = c("ATT(g=4,t=4)", "pre(g=3,t=1)", "ATT(g=3,t=3)", "ATT(g=3,t=4)", "ATT(g=4,t=5)"),
    support_status = c("estimable", "deterministic_reference_zero", "estimable",
                       "isolated_unsupported", "estimable"),
    support_reason = c("level anchored to the cohort reference",
                       "universal base-period reference cell: exact zero with zero influence function",
                       "level anchored to the cohort reference",
                       "post vertex isolated in the retained catalog",
                       "level anchored to the cohort reference"),
    available = c(TRUE, TRUE, TRUE, FALSE, TRUE), stringsAsFactors = FALSE)
}

fx_if <- function() {
  tid <- fx_target_meta()$target_id
  IF <- matrix(NA_real_, fx_n, 5L, dimnames = list(NULL, tid))
  IF[, 1L] <- c(2, -1, 3, -2, 1, -3, 0, 0)
  IF[, 2L] <- 0
  IF[, 3L] <- c(-1, 2, -2, 1, 3, -3, 0, 0)
  IF[, 5L] <- c(1, 1, -1, -1, 2, -2, 0, 0)
  IF
}
fx_est <- c(1.5, 0, 2.0, NA, 3.0)

mk_inference_block <- function(est, se, an_se, V, n, crit, cl, mult, cert, labels) {
  names(est) <- names(se) <- names(an_se) <- labels
  z <- ifelse(!is.na(se) & se > 0, est / se, NA_real_)
  names(z) <- labels
  av <- !is.na(est)
  list(status = if (all(av)) "all_requested_supported" else if (any(av)) "supported_subvector" else "empty_supported_selection",
       estimates = est, se = se, analytical_se = an_se, z_stat = z, p_value = 2 * stats::pnorm(-abs(z)),
       ci_lower = est - stats::qnorm(0.975) * se, ci_upper = est + stats::qnorm(0.975) * se,
       cband_crit_val = crit, cband_lower = est - crit * se, cband_upper = est + crit * se,
       vcov = V, Omega = V * n, cluster_if = NULL, raw_cluster_sums = NULL,
       n_total = n, n_clusters = length(unique(cl)),
       n_zero_eligible_clusters = sum(!unique(cl) %in% cl[fx_eligible]),
       clustered = TRUE, multiplier = mult, rng = NULL,
       support = list(status = "supported_subvector", available_ids = labels[av], unavailable_ids = labels[!av],
                      n_available = sum(av), n_requested = length(av), whole_request_complete = all(av)),
       numerical_certificate = cert)
}

mk_reporting <- function(bootstrap = FALSE) {
  tm <- fx_target_meta(); tid <- tm$target_id; av <- tm$available
  IF <- fx_if(); est <- fx_est; n <- fx_n; cl <- fx_cluster
  V <- matrix(NA_real_, 5L, 5L, dimnames = list(tid, tid))
  V[av, av] <- cl_cov(IF[, av], cl, n)
  an_se <- rep(NA_real_, 5L); an_se[av] <- sqrt(pmax(0, diag(V)[av]))
  se <- an_se; crit <- stats::qnorm(0.975); mult <- NULL
  family <- "pointwise_normal_over_certified_ids"
  if (bootstrap) {
    pos <- which(av & an_se > 0)
    se[pos] <- an_se[pos] * 1.1
    crit <- 2.4
    mult <- list(accepted = TRUE, crit_val = 2.4, bootstrap_se = se[pos],
                 raw_draws = matrix(0, 3L, length(pos)), metadata = list(n_draws = 99L))
    family <- "simultaneous_multiplier_sup_t_over_certified_ids"
  }
  cert <- list(certified = TRUE, certified_target_ids = tid[av], unavailable_target_ids = tid[!av],
               band_family = family, whole_request_complete = FALSE, scope = "test fixture")
  inf <- mk_inference_block(est, se, an_se, V, n, crit, cl, mult, cert, tid)
  c(list(base_period = "universal", target_meta = tm, L_map = NULL),
    inf[c("estimates", "se", "analytical_se", "numerical_certificate", "z_stat", "p_value",
          "ci_lower", "ci_upper", "cband_crit_val", "cband_lower", "cband_upper", "vcov")],
    list(Omega_a = V[av, av] * 7, Omega_a_enclosure = NULL, if_full = IF, available = av,
         support = inf$support, inference = inf))
}

fx_units <- function(cohort = fx_cohort, cluster = fx_cluster, eligible = fx_eligible) {
  u <- data.frame(id = fx_ids, cohort = cohort, weight = 1, cluster = cluster, eligible = eligible)
  u[rev(seq_len(nrow(u))), , drop = FALSE]                 # permuted row order on purpose
}

fx_params <- function(...) {
  p <- list(yname = "y", tname = "period", idname = "id", gname = "G", xformla = ~1,
            control_group = "nevertreated", anticipation = 0, alp = 0.05, tlist = 1:5,
            glist = c(3, 4), est_method = "ipw", base_period = "universal", clustervars = "cluster",
            bstrap = FALSE, cband = FALSE)
  utils::modifyList(p, list(...))
}
fx_call <- quote(att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = d, mode = "chained"))

mk_mp <- function(reporting = mk_reporting(), units = fx_units(), params = fx_params(), ...) {
  chained_did_MP(reporting, fx_ids, fx_call, params, units = units, ...)
}

# Aggregation fixtures ---------------------------------------------------------

fx_p <- c(`3` = 3 / 8, `4` = 2 / 8)                       # full-frame cohort shares, unit weights
fx_P <- cbind(`3` = as.numeric(fx_cohort == 3) - fx_p[["3"]],
              `4` = as.numeric(fx_cohort == 4) - fx_p[["4"]])

mk_cert <- function(labels, available, family) {
  list(certified = any(available), certified_target_ids = labels[available],
       unavailable_target_ids = labels[!available],
       band_family = if (!any(available)) "none" else family, whole_request_complete = all(available),
       scope = "test fixture")
}

mk_agg <- function(type, na_rm, keys, labels, target_ids, used_ids, est, fixed_if, share_if,
                   available, restricted, share_dependent, reasons, overall, bootstrap = FALSE) {
  n <- fx_n; cl <- fx_cluster; M <- length(labels)
  total <- fixed_if + share_if
  colnames(total) <- colnames(fixed_if) <- colnames(share_if) <- labels
  se <- se_fixed <- rep(NA_real_, M)
  if (any(available)) {
    se[available] <- cl_se(total[, available, drop = FALSE], cl, n)
    se_fixed[available] <- cl_se(fixed_if[, available, drop = FALSE], cl, n)
  }
  an_se <- se; crit <- if (any(available)) stats::qnorm(0.975) else NA_real_; mult <- NULL
  family <- "pointwise_normal_over_certified_ids"
  if (bootstrap) {
    pos <- which(available & an_se > 0)
    se[pos] <- an_se[pos] * 1.05; crit <- 2.3
    mult <- list(accepted = TRUE, crit_val = 2.3, bootstrap_se = se[pos], metadata = list(n_draws = 77L))
    family <- "simultaneous_multiplier_sup_t_over_certified_ids"
  }
  V <- matrix(NA_real_, M, M, dimnames = list(labels, labels))
  if (any(available)) V[available, available] <- cl_cov(total[, available, drop = FALSE], cl, n)
  cert <- mk_cert(labels, available, family)
  inf <- mk_inference_block(est, se, an_se, V, n, crit, cl, mult, cert, labels)
  names(se_fixed) <- labels
  excluded <- mapply(function(a, b) paste(setdiff(a, b), collapse = "; "), target_ids, used_ids)
  cs <- data.frame(label = labels, available = available, restricted = restricted,
                   share_dependent = share_dependent, n_intended = lengths(target_ids),
                   n_used = lengths(used_ids), share_denominator = NA_real_,
                   intended_ids = vapply(target_ids, paste, character(1L), collapse = "; "),
                   used_ids = vapply(used_ids, paste, character(1L), collapse = "; "),
                   excluded_ids = unname(excluded), reason = reasons, stringsAsFactors = FALSE)
  ov_cert <- list(certified = identical(overall$status, "available"),
                  certified_target_ids = if (identical(overall$status, "available")) "overall" else character(0L),
                  band_family = if (identical(overall$status, "available")) "pointwise_normal_over_certified_ids" else "none",
                  whole_request_complete = identical(overall$status, "available"))
  ov_se <- if (identical(overall$status, "available")) cl_se(matrix(overall$if_full, ncol = 1L), cl, n) else NA_real_
  overall_full <- c(list(estimate = overall$estimate, se = unname(ov_se),
                         ci_lower = overall$estimate - stats::qnorm(0.975) * unname(ov_se),
                         ci_upper = overall$estimate + stats::qnorm(0.975) * unname(ov_se),
                         if_full = overall$if_full, status = overall$status,
                         reason = if (is.null(overall$reason)) NA_character_ else overall$reason),
                    overall[c("intended_bins", "used_bins", "excluded_bins", "weights", "share_gradient",
                              "share_dependent", "share_denominator")],
                    list(numerical_certificate = ov_cert))
  structure(list(
    type = type, na_rm = na_rm,
    status = if (all(available) && identical(overall$status, "available")) "available" else if (any(available)) "partially_available" else "unavailable",
    cell_keys = keys, cell_labels = labels, cell_target_ids = target_ids, cell_used_target_ids = used_ids,
    cell_status = cs, target_definition = if (na_rm) "explicit restricted targets (test fixture)" else "intended targets (test fixture)",
    shares_known = TRUE, estimates = inf$estimates, se = inf$se, se_fixed_shares_only = se_fixed,
    ci_lower = inf$ci_lower, ci_upper = inf$ci_upper, cband_crit_val = crit,
    cband_lower = inf$cband_lower, cband_upper = inf$cband_upper,
    if_full = total, if_fixed_full = fixed_if, if_share_full = share_if,
    weights_theta = NULL, gradient_shares = NULL, cohort_shares = fx_p, cohort_share_enclosure = NULL,
    if_shares_full = fx_P, if_shares_enclosure = NULL, fixed_share_numerical_certificate = NULL,
    weight_enclosure = NULL, share_gradient_enclosure = NULL, cell_estimate_enclosure = NULL,
    decomposition_enclosure = NULL, decomposition = NULL,
    inference = inf, numerical_certificate = cert, overall = overall_full),
    class = "cdid_ipw_aggregation_result")
}

fx_dynamic <- function(bootstrap = FALSE) {
  IF <- fx_if(); n <- fx_n
  # event 0: ATT(3,3)=2.0 (p3) and ATT(4,4)=1.5 (p4); shares renormalized over {3,4}
  w0 <- fx_p / sum(fx_p); e0 <- sum(w0 * c(2.0, 1.5))
  fixed <- cbind(0, w0[["3"]] * IF[, 3L] + w0[["4"]] * IF[, 1L], IF[, 5L])
  D0 <- c(`3` = (2.0 - e0) / sum(fx_p), `4` = (1.5 - e0) / sum(fx_p))
  share <- cbind(0, as.numeric(fx_P %*% D0), 0)
  total <- fixed + share
  ov_if <- 0.5 * total[, 2L] + 0.5 * total[, 3L]
  mk_agg("dynamic", TRUE, keys = c(-2, 0, 1), labels = c("event_-2", "event_0", "event_1"),
         target_ids = list("pre(g=3,t=1)", c("ATT(g=3,t=3)", "ATT(g=4,t=4)"), c("ATT(g=3,t=4)", "ATT(g=4,t=5)")),
         used_ids = list("pre(g=3,t=1)", c("ATT(g=3,t=3)", "ATT(g=4,t=4)"), "ATT(g=4,t=5)"),
         est = c(0, e0, 3.0), fixed_if = fixed, share_if = share,
         available = c(TRUE, TRUE, TRUE), restricted = c(FALSE, FALSE, TRUE), share_dependent = c(FALSE, TRUE, FALSE),
         reasons = c("intended_target", "intended_target",
                     "restricted_target (na.rm = TRUE): weights renormalized over the supported subset"),
         overall = list(estimate = 0.5 * (e0 + 3.0), if_full = ov_if, status = "available", reason = NA_character_,
                        intended_bins = c("event_0", "event_1"), used_bins = c("event_0", "event_1"),
                        excluded_bins = character(0L), weights = c(`event_-2` = 0, event_0 = 0.5, event_1 = 0.5),
                        share_gradient = c(`3` = 0, `4` = 0), share_dependent = FALSE, share_denominator = NA_real_),
         bootstrap = bootstrap)
}

fx_group <- function() {
  IF <- fx_if()
  # na.rm = TRUE: group_3 restricted to ATT(3,3)=2.0; group_4 = mean(1.5, 3.0) = 2.25 (fixed weights)
  fixed <- cbind(IF[, 3L], 0.5 * IF[, 1L] + 0.5 * IF[, 5L])
  share <- matrix(0, fx_n, 2L)
  S_g <- sum(fx_p); ow <- fx_p / S_g; ov <- sum(ow * c(2.0, 2.25))
  og <- c(`3` = (2.0 - ov) / S_g, `4` = (2.25 - ov) / S_g)
  ov_if <- as.numeric(fixed %*% ow) + as.numeric(fx_P %*% og)
  mk_agg("group", TRUE, keys = c("3", "4"), labels = c("group_3", "group_4"),
         target_ids = list(c("ATT(g=3,t=3)", "ATT(g=3,t=4)"), c("ATT(g=4,t=4)", "ATT(g=4,t=5)")),
         used_ids = list("ATT(g=3,t=3)", c("ATT(g=4,t=4)", "ATT(g=4,t=5)")),
         est = c(2.0, 2.25), fixed_if = fixed, share_if = share,
         available = c(TRUE, TRUE), restricted = c(TRUE, FALSE), share_dependent = c(FALSE, FALSE),
         reasons = c("restricted_target (na.rm = TRUE): weights renormalized over the supported subset", "intended_target"),
         overall = list(estimate = ov, if_full = ov_if, status = "available", reason = NA_character_,
                        intended_bins = c("group_3", "group_4"), used_bins = c("group_3", "group_4"),
                        excluded_bins = character(0L), weights = stats::setNames(ow, c("group_3", "group_4")),
                        share_gradient = og, share_dependent = TRUE, share_denominator = S_g))
}

fx_calendar <- function() {
  IF <- fx_if()
  # na.rm = FALSE: calendar_4 needs ATT(3,4) (unavailable) -> NA; overall unavailable
  fixed <- cbind(IF[, 3L], NA_real_, IF[, 5L])
  share <- cbind(rep(0, fx_n), NA_real_, rep(0, fx_n))
  mk_agg("calendar", FALSE, keys = c(3, 4, 5), labels = c("calendar_3", "calendar_4", "calendar_5"),
         target_ids = list("ATT(g=3,t=3)", c("ATT(g=3,t=4)", "ATT(g=4,t=4)"), "ATT(g=4,t=5)"),
         used_ids = list("ATT(g=3,t=3)", character(0L), "ATT(g=4,t=5)"),
         est = c(2.0, NA, 3.0), fixed_if = fixed, share_if = share,
         available = c(TRUE, FALSE, TRUE), restricted = c(FALSE, FALSE, FALSE), share_dependent = c(FALSE, TRUE, FALSE),
         reasons = c("intended_target", "missing_targets (na.rm = FALSE): ATT(g=3,t=4) [isolated_unsupported]", "intended_target"),
         overall = list(estimate = NA_real_, if_full = rep(NA_real_, fx_n), status = "unavailable",
                        reason = "missing_inner_aggregates (na.rm = FALSE): calendar_4",
                        intended_bins = c("calendar_3", "calendar_4", "calendar_5"), used_bins = character(0L),
                        excluded_bins = c("calendar_3", "calendar_4", "calendar_5"),
                        weights = c(calendar_3 = 0, calendar_4 = 0, calendar_5 = 0),
                        share_gradient = c(`3` = 0, `4` = 0), share_dependent = FALSE, share_denominator = NA_real_))
}

fx_simple <- function() {
  IF <- fx_if()
  den <- fx_p[["3"]] + 2 * fx_p[["4"]]
  w <- c(fx_p[["3"]], fx_p[["4"]], fx_p[["4"]]) / den          # ATT(3,3), ATT(4,4), ATT(4,5)
  est <- sum(w * c(2.0, 1.5, 3.0))
  fixed <- matrix(w[1] * IF[, 3L] + w[2] * IF[, 1L] + w[3] * IF[, 5L], ncol = 1L)
  D <- c(`3` = (2.0 - est) / den, `4` = ((1.5 - est) + (3.0 - est)) / den)
  share <- matrix(as.numeric(fx_P %*% D), ncol = 1L)
  mk_agg("simple", TRUE, keys = "simple", labels = "ATT_simple",
         target_ids = list(c("ATT(g=3,t=3)", "ATT(g=3,t=4)", "ATT(g=4,t=4)", "ATT(g=4,t=5)")),
         used_ids = list(c("ATT(g=3,t=3)", "ATT(g=4,t=4)", "ATT(g=4,t=5)")),
         est = est, fixed_if = fixed, share_if = share,
         available = TRUE, restricted = TRUE, share_dependent = TRUE,
         reasons = "restricted_target (na.rm = TRUE): weights renormalized over the supported subset",
         overall = list(estimate = est, if_full = as.numeric(fixed + share), status = "available", reason = NA_character_,
                        intended_bins = "ATT_simple", used_bins = "ATT_simple", excluded_bins = character(0L),
                        weights = c(ATT_simple = 1), share_gradient = c(`3` = 0, `4` = 0),
                        share_dependent = FALSE, share_denominator = NA_real_))
}

# Shape produced by ipw_pipeline_unavailable_aggregation(): no post cell declared.
fx_unavailable_group <- function() {
  n <- fx_n; labels <- c("group_3", "group_4"); na_v <- stats::setNames(rep(NA_real_, 2L), labels)
  inf <- list(status = "empty_supported_selection", reason = "no_declared_post_cells", estimates = na_v, se = na_v,
              analytical_se = na_v, z_stat = na_v, p_value = na_v, ci_lower = na_v, ci_upper = na_v,
              cband_crit_val = NA_real_, cband_lower = na_v, cband_upper = na_v,
              vcov = matrix(NA_real_, 2L, 2L, dimnames = list(labels, labels)),
              Omega = matrix(NA_real_, 2L, 2L, dimnames = list(labels, labels)),
              n_total = n, n_clusters = 5L, n_zero_eligible_clusters = 1L, clustered = TRUE,
              multiplier = NULL, rng = NULL,
              support = list(status = "empty_supported_selection", available_ids = character(0L),
                             unavailable_ids = labels, n_available = 0L, n_requested = 2L, whole_request_complete = FALSE),
              numerical_certificate = list(certified = FALSE, status = "no_supported_targets",
                                           certified_target_ids = character(0L), unavailable_target_ids = labels,
                                           band_family = "none", whole_request_complete = FALSE, scope = "nothing certified"))
  structure(list(type = "group", na_rm = FALSE, status = "unavailable", reason = "no_declared_post_cells",
                 cell_keys = c("3", "4"), cell_labels = labels,
                 cell_target_ids = list("ATT(g=3,t=3)", "ATT(g=4,t=4)"),
                 cell_status = data.frame(label = labels, available = FALSE, restricted = FALSE,
                                          reason = "no_declared_post_cells", stringsAsFactors = FALSE),
                 estimates = na_v, se = na_v, se_fixed_shares_only = na_v, ci_lower = na_v, ci_upper = na_v,
                 cband_crit_val = NA_real_, cband_lower = na_v, cband_upper = na_v,
                 if_full = matrix(NA_real_, n, 2L, dimnames = list(NULL, labels)),
                 if_fixed_full = matrix(NA_real_, n, 2L), if_share_full = matrix(NA_real_, n, 2L),
                 weights_theta = NULL, gradient_shares = NULL, cohort_shares = NULL, shares_known = TRUE,
                 cohort_share_enclosure = NULL, if_shares_full = NULL, if_shares_enclosure = NULL,
                 fixed_share_numerical_certificate = NULL, decomposition_enclosure = NULL, decomposition = NULL,
                 inference = inf, numerical_certificate = inf$numerical_certificate,
                 overall = list(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                                if_full = rep(NA_real_, n), status = "unavailable", reason = "no_declared_post_cells",
                                numerical_certificate = inf$numerical_certificate)),
            class = "cdid_ipw_aggregation_result")
}

upstream_tidy_mp_cols <- c("term", "group", "time", "estimate", "std.error", "statistic", "p.value",
                           "conf.low", "conf.high", "point.conf.low", "point.conf.high")

# ------------------------------------------------------------------------------
# cdid_MP construction
# ------------------------------------------------------------------------------

test_that("cdid_MP maps reported cells in exact order onto did::MP fields", {
  mp <- mk_mp()
  tm <- fx_target_meta(); tid <- tm$target_id; IF <- fx_if()
  expect_s3_class(mp, c("cdid_MP", "MP"), exact = TRUE)
  expect_true(inherits(mp, "MP"))
  expect_identical(mp$group, tm$g)
  expect_identical(mp$t, tm$time)
  expect_identical(mp$att, unname(fx_est))
  expect_identical(mp$n, fx_n)
  expect_identical(mp$call, fx_call)
  expect_identical(mp$DIDparams$call, fx_call)
  expect_null(mp$aggte)
  # V_analytical is Omega = n * vcov, on did's scale se = sqrt(diag(V)/n)
  rep <- mk_reporting()
  expect_identical(unname(mp$V_analytical), unname(rep$inference$Omega))
  expect_false(isTRUE(all.equal(unname(mp$V_analytical), unname(rep$vcov))))
  expect_false(isTRUE(all.equal(dim(mp$V_analytical), dim(rep$Omega_a))))
  av <- tm$available
  expect_equal(mp$se[av], unname(sqrt(diag(mp$V_analytical)[av] / mp$n)))
  expect_true(all(is.na(mp$V_analytical[!av, ])) && all(is.na(mp$V_analytical[, !av])))
  expect_identical(dimnames(mp$V_analytical), list(tid, tid))
  # IF: full-frame rows identified by the caller's IDs, NA column retained
  expect_identical(dim(mp$inffunc), c(fx_n, 5L))
  expect_identical(rownames(mp$inffunc), as.character(fx_ids))
  expect_identical(colnames(mp$inffunc), tid)
  expect_identical(unname(mp$inffunc[, av]), unname(IF[, av]))
  expect_true(all(is.na(mp$inffunc[, !av])))
  expect_identical(unname(mp$inffunc[8L, av]), c(0, 0, 0, 0))   # frame-only unit keeps its zero row
  expect_equal(mp$c, stats::qnorm(0.975))
  expect_identical(mp$alp, 0.05)
  expect_null(mp$W); expect_null(mp$Wpval)
  expect_identical(mp$chained$pretest$status, "not_computed")
})

test_that("cdid_MP DIDparams is honest and carries both count families", {
  mp <- mk_mp()
  dp <- mp$DIDparams
  expect_s3_class(dp, "DIDparams")
  expect_identical(dp$est_method, "ipw")
  expect_identical(dp$control_group, "nevertreated")
  expect_identical(dp$base_period, "universal")
  expect_identical(dp$anticipation, 0)
  expect_false(dp$bstrap); expect_false(dp$cband)
  expect_identical(dp$biters, formals(did::DIDparams)$biters)
  expect_false(dp$faster_mode)
  expect_true(dp$panel); expect_false(dp$true_repeated_cross_sections)
  expect_null(dp$data)
  expect_identical(dp$n, fx_n); expect_identical(dp$id_count, fx_n)
  expect_identical(dp$nG, 2L); expect_identical(dp$treated_groups_count, 2L)
  expect_identical(dp$nT, 5L); expect_identical(dp$time_periods_count, 5L)
  expect_identical(dp$tlist, as.double(1:5)); expect_identical(dp$time_periods, as.double(1:5))
  expect_identical(dp$glist, c(3, 4)); expect_identical(dp$treated_groups, c(3, 4))
  expect_identical(dp$clustervars, "cluster")
  expect_identical(dp$cluster_vector, fx_cluster)       # IF-row order, not the permuted input order
  expect_identical(dp$cluster_vector_var, "cluster")
  expect_identical(dp$cdid$mode, "chained")
  expect_false(dp$cdid$data_retained)
  expect_identical(dp$cdid$versions$did_tested, "2.5.1")
  expect_identical(nobs(mp), fx_n)
})

test_that("unit metadata is aligned by identity from permuted input and checked against inference", {
  mp <- mk_mp()
  fr <- mp$chained$frame
  expect_identical(fr$id, fx_ids)
  expect_identical(fr$cohort, fx_cohort)
  expect_identical(fr$cluster, fx_cluster)
  expect_identical(fr$eligible, fx_eligible)
  expect_identical(fr$weight, rep(1, fx_n))
  # frame-only unit retained with its (known) cohort and zero-eligible cluster
  expect_false(fr$eligible[fr$id == 108L])
  expect_identical(mp$chained$inference$n_zero_eligible_clusters, 1L)
  expect_identical(mp$chained$inference$n_clusters, 5L)
  # unknown frame-only cohort stays NA (never inferred)
  mp2 <- mk_mp(units = fx_units(cohort = replace(fx_cohort, 8L, NA)))
  expect_true(is.na(mp2$chained$frame$cohort[8L]))
  # metadata without any optional columns is accepted
  mp3 <- chained_did_MP(mk_reporting(), fx_ids, fx_call, fx_params(),
                        units = data.frame(id = rev(fx_ids), cluster = rev(fx_cluster)))
  expect_identical(mp3$chained$frame$cluster, fx_cluster)
  expect_true(all(is.na(mp3$chained$frame$cohort)))
})

test_that("cdid_MP refuses unfaithful inputs instead of repairing them", {
  rep <- mk_reporting()
  # permuted estimates relative to target_meta
  bad <- rep; bad$estimates <- bad$estimates[c(2, 1, 3, 4, 5)]
  expect_error(mk_mp(reporting = bad), class = "cdid_result_adapter_error")
  # zero-filled unavailable IF column
  bad <- rep; bad$if_full[, 4L] <- 0
  expect_error(mk_mp(reporting = bad), "never zero-filled", class = "cdid_result_adapter_error")
  # covariance not on the n scale
  bad <- rep; bad$inference$Omega <- bad$inference$vcov
  expect_error(mk_mp(reporting = bad), class = "cdid_result_adapter_error")
  # Omega_a (eligible-frame) in place of Omega: wrong dimension
  bad <- rep; bad$inference$Omega <- bad$Omega_a; bad$inference$vcov <- NULL
  expect_error(mk_mp(reporting = bad), class = "cdid_result_adapter_error")
  # duplicate / short frame IDs
  expect_error(chained_did_MP(rep, c(fx_ids[-8L], 101L), fx_call, fx_params(), units = fx_units()),
               "frame_ids", class = "cdid_result_adapter_error")
  expect_error(chained_did_MP(rep, fx_ids[-8L], fx_call, fx_params(), units = fx_units()),
               class = "cdid_result_adapter_error")
  # IF rownames that disagree with the explicit frame IDs
  bad <- rep; rownames(bad$if_full) <- as.character(rev(fx_ids))
  expect_error(mk_mp(reporting = bad), "rownames", class = "cdid_result_adapter_error")
  # clustered inference without the aligned cluster labels / wrong cluster count
  expect_error(mk_mp(units = fx_units()[, c("id", "cohort")]), "clustered", class = "cdid_result_adapter_error")
  expect_error(mk_mp(units = fx_units(cluster = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L))),
               "n_clusters", class = "cdid_result_adapter_error")
  # eligibility that contradicts the reported zero-eligible cluster count
  expect_error(mk_mp(units = fx_units(eligible = rep(TRUE, 8L))), "zero-eligible", class = "cdid_result_adapter_error")
  # missing / inconsistent explicit metadata
  expect_error(mk_mp(params = fx_params(tlist = NULL)), "explicit caller metadata", class = "cdid_result_adapter_error")
  expect_error(mk_mp(params = fx_params(glist = 4)), "glist", class = "cdid_result_adapter_error")
  expect_error(mk_mp(params = fx_params(alp = 0.10)), "critical value", class = "cdid_result_adapter_error")
  expect_error(mk_mp(params = fx_params(clustervars = NULL)), "clustervars", class = "cdid_result_adapter_error")
  expect_error(mk_mp(params = fx_params(panel = FALSE)), "panel", class = "cdid_result_adapter_error")
  # a unit with -Inf cohort is not a cohort
  expect_error(mk_mp(units = fx_units(cohort = replace(fx_cohort, 1L, -Inf))), "cohort", class = "cdid_result_adapter_error")
})

test_that("DIDparams preserves requests while chained inference records the computation", {
  mp <- mk_mp()
  expect_identical(mp$chained$inference$kind, "analytic")
  expect_identical(mp$chained$band_family, "pointwise_normal_over_certified_ids")
  expect_identical(mp$chained$band_family_source, "numerical_certificate")
  # analytic input + requested bootstrap -> refusal, not a silent bootstrap claim
  expect_error(mk_mp(params = fx_params(bstrap = TRUE)), "analytic", class = "cdid_result_adapter_error")
  expect_error(mk_mp(params = fx_params(cband = TRUE)), "band", class = "cdid_result_adapter_error")
  # consistent requested flags are accepted and recorded
  mp2 <- mk_mp(params = fx_params(bstrap = FALSE, cband = FALSE, biters = 500L))
  expect_false(mp2$DIDparams$bstrap)
  expect_identical(mp2$DIDparams$cdid$requested_bstrap, FALSE)
  expect_identical(mp2$DIDparams$cdid$requested_biters, 500L)
  expect_identical(mp2$DIDparams$biters, 500L)
  expect_true(is.na(mp2$chained$inference$biters))
  # Requested iteration count is preserved separately from the actual draws.
  mpb <- mk_mp(reporting = mk_reporting(bootstrap = TRUE),
                params = fx_params(bstrap = TRUE, cband = TRUE, biters = 999L))
  expect_true(mpb$DIDparams$bstrap); expect_true(mpb$DIDparams$cband)
  expect_identical(mpb$DIDparams$biters, 999L)
  expect_identical(mpb$chained$inference$biters, 99L)
  expect_identical(mpb$chained$inference$kind, "multiplier_bootstrap_iqr_se")
  expect_identical(mpb$c, 2.4)
  expect_null(mpb$chained$inference$multiplier$raw_draws)
  expect_error(mk_mp(reporting = mk_reporting(bootstrap = TRUE), params = fx_params(bstrap = FALSE)),
               "bootstrap", class = "cdid_result_adapter_error")
  # a certificate family that contradicts the multiplier presence is refused
  bad <- mk_reporting(); bad$numerical_certificate$band_family <- "simultaneous_multiplier_sup_t_over_certified_ids"
  expect_error(mk_mp(reporting = bad), "band family", class = "cdid_result_adapter_error")
  # no certificate: family derived from the multiplier presence, labelled as such
  nocert <- mk_reporting(); nocert$numerical_certificate <- NULL; nocert$inference$numerical_certificate <- NULL
  mpn <- mk_mp(reporting = nocert)
  expect_identical(mpn$chained$band_family_source, "inference_multiplier_presence")
  expect_identical(mpn$chained$band_family, "pointwise_normal_over_certified_ids")
})

test_that("released bstrap = FALSE, cband = TRUE endpoint laws are accepted only when the pipeline recorded them", {
  cells_law <- "released_bstrap_false_cband_true_cells_analytic_pointwise_no_draws"
  bins_law <- "released_bstrap_false_cband_true_aggte_bins_analytic_se_multiplier_sup_t_band"
  # att_gt endpoint: analytic pointwise inference recorded under the cells law
  rp <- mk_reporting(); rp$inference$bootstrap <- FALSE; rp$inference$simultaneous <- FALSE; rp$inference$flag_law <- cells_law
  mp <- mk_mp(reporting = rp, params = fx_params(bstrap = FALSE, cband = TRUE, biters = 250L))
  expect_false(mp$DIDparams$bstrap); expect_true(mp$DIDparams$cband); expect_identical(mp$DIDparams$biters, 250L)
  expect_identical(mp$DIDparams$cdid$requested_cband, TRUE); expect_identical(mp$DIDparams$cdid$effective_cband, FALSE)
  expect_identical(mp$DIDparams$cdid$flag_law, cells_law)
  expect_false(mp$chained$inference$bstrap); expect_false(mp$chained$inference$cband); expect_identical(mp$chained$inference$kind, "analytic")
  expect_true(mp$chained$inference$released_false_true); expect_identical(mp$chained$inference$endpoint, "cells")
  expect_true(is.na(mp$chained$inference$biters)); expect_identical(mp$c, stats::qnorm(0.975))
  txt <- paste(utils::capture.output(summary(mp)), collapse = "\n")
  expect_match(txt, "Pointwise", fixed = TRUE); expect_no_match(txt, "Simult.", fixed = TRUE)
  expect_match(txt, "released did 2.5.1 att_gt", fixed = TRUE)
  # the same record without the law is still the generic contradiction (refused)
  rp2 <- rp; rp2$inference$flag_law <- "explicit_bstrap_and_cband"
  expect_error(mk_mp(reporting = rp2, params = fx_params(bstrap = FALSE, cband = TRUE)), "band", class = "cdid_result_adapter_error")
  rp3 <- rp; rp3$inference$flag_law <- NULL
  expect_error(mk_mp(reporting = rp3, params = fx_params(bstrap = FALSE, cband = TRUE)), "band", class = "cdid_result_adapter_error")
  # the bins law is not an att_gt-endpoint law: refused for reported cells
  rp4 <- mk_reporting(bootstrap = TRUE); rp4$inference$bootstrap <- FALSE; rp4$inference$simultaneous <- TRUE; rp4$inference$flag_law <- bins_law
  rp4$se <- rp4$analytical_se; rp4$inference$se <- rp4$analytical_se
  expect_error(mk_mp(reporting = rp4, params = fx_params(bstrap = FALSE, cband = TRUE)), "released", class = "cdid_result_adapter_error")
  # bootstrap draws with analytic SEs but no recorded law: refused (SE kind never follows draw presence silently)
  rp5 <- mk_reporting(bootstrap = TRUE); rp5$inference$bootstrap <- FALSE
  expect_error(mk_mp(reporting = rp5, params = fx_params(bstrap = FALSE, cband = TRUE)), "analytic without", class = "cdid_result_adapter_error")
  # the cells law cannot claim bootstrap SEs or a band either
  rp6 <- mk_reporting(bootstrap = TRUE); rp6$inference$flag_law <- cells_law
  expect_error(mk_mp(reporting = rp6, params = fx_params(bstrap = FALSE, cband = TRUE)), class = "cdid_result_adapter_error")
  # aggte-bin endpoint: analytic SEs with a multiplier sup-t constant recorded under the bins law
  a <- fx_dynamic(bootstrap = TRUE)
  a$inference$bootstrap <- FALSE; a$inference$simultaneous <- TRUE; a$inference$flag_law <- bins_law
  an <- a$inference$analytical_se; a$se <- an; a$inference$se <- an
  a$cband_lower <- a$inference$cband_lower <- a$estimates - 2.3 * an; a$cband_upper <- a$inference$cband_upper <- a$estimates + 2.3 * an
  base_mp <- mk_mp(params = fx_params(bstrap = FALSE, cband = TRUE), reporting = { r0 <- mk_reporting(); r0$inference$flag_law <- cells_law; r0$inference$bootstrap <- FALSE; r0$inference$simultaneous <- FALSE; r0 })
  req <- list(bstrap = FALSE, cband = TRUE, biters = 77L, alp = 0.05, clustervars = NULL,
              sources = c(bstrap = "inherited", cband = "inherited", biters = "explicit", alp = "inherited", clustervars = "inherited"))
  ab <- chained_did_AGGTEobj(a, base_mp, fx_call, requested = req)
  expect_false(ab$DIDparams$bstrap); expect_true(ab$DIDparams$cband)
  expect_identical(ab$chained$inference$kind, "analytic"); expect_true(ab$chained$inference$cband)
  expect_identical(ab$chained$inference$biters, 77L); expect_identical(ab$crit.val.egt, 2.3)
  expect_identical(ab$se.egt, unname(an)); expect_true(ab$chained$inference$released_false_true)
  expect_identical(ab$chained$overall$inference_kind, "analytic")
  txt <- paste(utils::capture.output(summary(ab)), collapse = "\n")
  expect_match(txt, "Simult.", fixed = TRUE); expect_match(txt, "analytic standard errors", fixed = TRUE)
  expect_match(txt, "released did 2.5.1 aggte", fixed = TRUE)
  # simple endpoint never carries the bins law
  s <- fx_simple(); s$inference$flag_law <- bins_law; s$inference$bootstrap <- FALSE; s$inference$simultaneous <- FALSE
  expect_error(chained_did_AGGTEobj(s, base_mp, fx_call, requested = req), class = "cdid_result_adapter_error")
  s2 <- fx_simple(); s2$inference$flag_law <- cells_law; s2$inference$bootstrap <- FALSE; s2$inference$simultaneous <- FALSE
  sb <- chained_did_AGGTEobj(s2, base_mp, fx_call, requested = req)
  expect_true(sb$DIDparams$cband); expect_false(sb$chained$inference$cband); expect_identical(sb$chained$inference$kind, "analytic")
  # requested TRUE/TRUE cannot borrow the bins law to hide analytic SEs
  req_tt <- req; req_tt$bstrap <- TRUE
  expect_error(chained_did_AGGTEobj(a, base_mp, fx_call, requested = req_tt), "bootstrap standard errors cannot be claimed", class = "cdid_result_adapter_error")
})

test_that("zero and empty laws retain upstream bootstrap requests without requiring draws", {
  zero <- mk_reporting()
  zero$target_meta <- zero$target_meta[2L, , drop = FALSE]
  for (nm in c("estimates", "se", "analytical_se", "z_stat", "p_value", "ci_lower", "ci_upper",
               "cband_lower", "cband_upper", "available")) zero[[nm]] <- zero[[nm]][2L]
  zero$if_full <- zero$if_full[, 2L, drop = FALSE]
  zero$inference$Omega <- zero$inference$Omega[2L, 2L, drop = FALSE]
  zero$inference$vcov <- zero$inference$vcov[2L, 2L, drop = FALSE]
  zero$cband_crit_val <- 0
  zero$numerical_certificate$certified_target_ids <- zero$target_meta$target_id
  zero$numerical_certificate$unavailable_target_ids <- character(0L)
  zero$numerical_certificate$band_family <- "deterministic_exact_zero_over_certified_ids"
  # Omitted flags/count use DIDparams defaults; the frozen pipeline needs no draws.
  params <- fx_params(bstrap = NULL, cband = NULL)
  mp <- mk_mp(reporting = zero, params = params)
  expect_true(mp$DIDparams$bstrap); expect_true(mp$DIDparams$cband)
  expect_identical(mp$DIDparams$biters, formals(did::DIDparams)$biters)
  expect_true(mp$chained$inference$zero_law)
  expect_false(mp$chained$inference$bstrap); expect_null(mp$chained$inference$multiplier)
  expect_identical(mp$c, 0)
  expect_true(all(is.na(tidy(mp)$statistic)))
  txt <- paste(utils::capture.output(summary(mp)), collapse = "\n")
  expect_match(txt, "Pointwise", fixed = TRUE); expect_no_match(txt, "Simult.", fixed = TRUE)
  expect_identical(summary(mp), mp)
  empty <- zero
  for (nm in c("estimates", "se", "analytical_se", "z_stat", "p_value", "ci_lower", "ci_upper",
               "cband_lower", "cband_upper")) empty[[nm]][] <- NA_real_
  empty$if_full[] <- empty$inference$Omega[] <- empty$inference$vcov[] <- NA_real_
  empty$available[] <- empty$target_meta$available[] <- FALSE
  empty$cband_crit_val <- NA_real_
  empty$numerical_certificate$certified_target_ids <- character(0L)
  empty$numerical_certificate$band_family <- "none"
  empty$inference$status <- "empty_supported_selection"
  e <- mk_mp(reporting = empty, params = params)
  expect_true(e$DIDparams$bstrap); expect_true(e$DIDparams$cband)
  expect_true(e$chained$inference$zero_law)
  expect_true(all(is.na(tidy(e)$estimate)))
  expect_error(mk_mp(params = params), "analytic", class = "cdid_result_adapter_error")
})

test_that("basis contrasts, separate band families and a supplied pretest are carried apart from cells", {
  n <- fx_n; tid <- fx_target_meta()$target_id
  bmeta <- data.frame(basis_id = c("ATT(g=3,t=4)-ATT(g=3,t=5)"), component_id = 2L, stringsAsFactors = FALSE)
  basis <- list(meta = bmeta, estimates = c(`ATT(g=3,t=4)-ATT(g=3,t=5)` = 0.4),
                se = c(`ATT(g=3,t=4)-ATT(g=3,t=5)` = 0.2), if_full = matrix(c(1, -1, 0, 0, 0, 0, 0, 0), n, 1L),
                cross_vcov_with_reported_cells = matrix(NA_real_, 5L, 1L, dimnames = list(tid, "ATT(g=3,t=4)-ATT(g=3,t=5)")),
                numerical_certificate = list(certified = TRUE, certified_target_ids = "ATT(g=3,t=4)-ATT(g=3,t=5)",
                                             band_family = "pointwise_normal_over_certified_ids"))
  families <- list(reported_cells = list(certified_target_ids = tid[-4L], band_family = "pointwise_normal_over_certified_ids"),
                   post_basis = list(certified_target_ids = "ATT(g=3,t=4)-ATT(g=3,t=5)", band_family = "pointwise_normal_over_certified_ids"))
  pre <- list(W = 1.234, Wpval = 0.267, cells = "pre(g=3,t=1)", df = 1, description = "fixture: supplied externally")
  mp <- mk_mp(basis = basis, band_families = families, pretest = pre)
  expect_length(mp$group, 5L)                                   # no contrast inserted as a cell
  expect_identical(mp$chained$basis$meta$basis_id, bmeta$basis_id)
  expect_identical(rownames(mp$chained$basis$if_full), as.character(fx_ids))
  expect_identical(names(mp$chained$band_families), c("reported_cells", "post_basis"))
  expect_identical(mp$W, 1.234); expect_identical(mp$Wpval, 0.267)
  expect_identical(mp$chained$pretest$status, "supplied_precomputed")
  # a contrast ID colliding with a cell ID, or a cross-covariance in another order, is refused
  bad <- basis; bad$meta$basis_id <- "ATT(g=3,t=3)"; names(bad$estimates) <- "ATT(g=3,t=3)"
  expect_error(mk_mp(basis = bad), "never cohort/date cells", class = "cdid_result_adapter_error")
  bad <- basis; rownames(bad$cross_vcov_with_reported_cells) <- rev(tid)
  expect_error(mk_mp(basis = bad), class = "cdid_result_adapter_error")
  bad <- basis; rownames(bad$if_full) <- as.character(rev(fx_ids))
  expect_error(mk_mp(basis = bad), "rownames", class = "cdid_result_adapter_error")
  aligned <- basis; rownames(aligned$if_full) <- as.character(fx_ids)
  expect_identical(mk_mp(basis = aligned)$chained$basis$if_full, aligned$if_full)
  # an incomplete pretest record is refused (absence stays explicit instead)
  expect_error(mk_mp(pretest = list(W = 1.2, Wpval = 0.3)), "pretest", class = "cdid_result_adapter_error")
  expect_null(mk_mp()$W)
})

# ------------------------------------------------------------------------------
# cdid_MP methods
# ------------------------------------------------------------------------------

test_that("tidy.cdid_MP keeps did's columns, uses the NA convention and appends support metadata", {
  mp <- mk_mp()
  td <- tidy(mp)
  expect_s3_class(td, "data.frame")
  expect_identical(names(td)[seq_along(upstream_tidy_mp_cols)], upstream_tidy_mp_cols)
  expect_identical(names(td)[-seq_along(upstream_tidy_mp_cols)],
                   c("target.id", "role", "available", "support.status", "support.reason", "certified", "band.family"))
  expect_identical(nrow(td), 5L)
  expect_identical(td$term, paste0("ATT(", mp$group, ",", mp$t, ")"))
  expect_identical(td$estimate, mp$att)
  # exact-zero reference cell: statistic/p-value NA (not NaN), band collapses to 0
  expect_true(is.na(td$statistic[2L]) && !is.nan(td$statistic[2L]))
  expect_true(is.na(td$p.value[2L]) && !is.nan(td$p.value[2L]))
  expect_identical(td$conf.low[2L], 0); expect_identical(td$conf.high[2L], 0)
  expect_false(any(is.nan(td$statistic) | is.nan(td$p.value) | is.infinite(td$statistic), na.rm = TRUE))
  # unavailable cell: NA throughout, reason retained
  expect_true(is.na(td$estimate[4L]) && is.na(td$conf.low[4L]) && is.na(td$p.value[4L]))
  expect_false(td$available[4L]); expect_identical(td$support.status[4L], "isolated_unsupported")
  expect_false(td$certified[4L]); expect_true(all(td$certified[-4L]))
  # available cells reproduce did's arithmetic
  av <- c(1L, 3L, 5L)
  expect_equal(td$statistic[av], mp$att[av] / mp$se[av])
  expect_equal(td$conf.low[av], mp$att[av] - mp$c * mp$se[av])
  expect_identical(td$target.id, fx_target_meta()$target_id)
  expect_identical(td$role, c("post", "pre", "post", "post", "post"))
  expect_true(all(td$band.family == "pointwise_normal_over_certified_ids"))
})

test_that("glance.cdid_MP extends did's glance with chained diagnostics", {
  mp <- mk_mp()
  g <- glance(mp)
  expect_identical(names(g)[1:5], c("nobs", "ngroup", "ntime", "control.group", "est.method"))
  expect_identical(g$nobs, fx_n); expect_identical(g$ngroup, 2L); expect_identical(g$ntime, 5L)
  expect_identical(g$est.method, "ipw"); expect_identical(g$control.group, "nevertreated")
  expect_identical(g$estimator, "chained IPW-GMM")
  expect_identical(g$n.requested, 5L); expect_identical(g$n.available, 4L)
  expect_identical(g$inference, "analytic"); expect_true(g$clustered)

  mp_dr <- mk_mp(params = fx_params(est_method = "dr"))
  expect_identical(mp_dr$chained$estimator, "chained doubly robust GMM")
  expect_identical(glance(mp_dr)$estimator, "chained doubly robust GMM")
  expect_identical(unique(tidy(mp_dr)$estimator), "chained doubly robust GMM")
  ag_dr <- chained_did_AGGTEobj(fx_dynamic(), mp_dr, quote(aggte(out, type = "dynamic", na.rm = TRUE)))
  expect_identical(glance(ag_dr)$estimator, "chained doubly robust GMM")
  expect_identical(unique(tidy(ag_dr)$estimator), "chained doubly robust GMM")
})

test_that("summary/print.cdid_MP identify chained GMM and the D5 state concisely", {
  mp <- mk_mp()
  out <- paste(utils::capture.output(summary(mp)), collapse = "\n")
  expect_match(out, "cdid chained difference-in-differences (chained IPW-GMM):", fixed = TRUE)
  # Both references appear once in the heading; did's duplicate citation is removed.
  expect_match(out, "Chained estimators: Bellego C., Benatia D. and V. Dortet-Bernadet (2025).", fixed = TRUE)
  expect_match(out, "Journal of Econometrics 248, 105783.", fixed = TRUE)
  expect_match(out, "Ordinary estimators: Callaway", fixed = TRUE)
  expect_no_match(out, "Reference:", fixed = TRUE)
  expect_identical(names(mp$chained$references), c("estimator", "cdid", "did"))
  expect_match(out, "Group-Time Average Treatment Effects", fixed = TRUE)   # did's table
  expect_match(out, "4 of 5 requested cells available (unavailable cells are NA, not zero effects)", fixed = TRUE)
  expect_match(out, "ATT(g=3,t=4) [isolated_unsupported]", fixed = TRUE)
  expect_match(out, "Exact-zero reference cells (no test statistic): 1.", fixed = TRUE)
  expect_match(out, "Units: 8 in the unit table, 7 with at least two observed outcomes.", fixed = TRUE)
  expect_match(out, "analytic standard errors, clustered on 5 clusters (1 without eligible units); pointwise intervals.", fixed = TRUE)
  expect_match(out, "pre-test: not computed", fixed = TRUE)
  expect_match(out, "Pointwise", fixed = TRUE)
  expect_no_match(out, "Simult.", fixed = TRUE)
  # Internal identifiers and records are not part of the normal display.
  for (internal in c("pointwise_normal_over_certified_ids", "numerical certificate", "Band family",
                     "Numerical policy", "Resource profile", "full-frame"))
    expect_no_match(out, internal, fixed = TRUE)
  expect_identical(utils::capture.output(print(mp)), utils::capture.output(summary(mp)))
  outb <- paste(utils::capture.output(summary(mk_mp(reporting = mk_reporting(bootstrap = TRUE),
                                                     params = fx_params(bstrap = TRUE, cband = TRUE)))), collapse = "\n")
  expect_match(outb, "Simult.", fixed = TRUE)
  expect_match(outb, "multiplier-bootstrap IQR standard errors (99 draws)", fixed = TRUE)
  expect_match(outb, "simultaneous band over the available cells.", fixed = TRUE)
})

test_that("the chained display keeps every line of did's own summary except its citation", {
  for (obj in list(mk_mp(), chained_did_AGGTEobj(fx_dynamic(), mk_mp(), quote(aggte(out, type = "dynamic"))))) {
    shown <- obj
    shown$DIDparams <- chained_result_presentation_params(obj)
    up <- utils::capture.output(utils::getS3method("summary", class(obj)[2L], envir = asNamespace("did"))(shown))
    cite <- grep("^Reference: ", up)
    expect_true(length(cite) <= 1L)
    kept <- if (length(cite)) up[-cite] else up
    kept <- kept[nzchar(trimws(kept))]
    ours <- utils::capture.output(summary(obj))
    expect_length(grep("^Reference: ", ours), 0L)
    pos <- match(kept, ours)
    expect_false(anyNA(pos))
    expect_false(is.unsorted(pos, strictly = TRUE))
  }
  # The citation filter removes only the citation and keeps one separating blank line.
  mp_lines <- c("", "Call:", "att_gt(y)", "", "Reference: Callaway ...", "", "Group-Time Average Treatment Effects:", "x")
  expect_identical(chained_result_drop_citation(mp_lines),
                   c("", "Call:", "att_gt(y)", "", "Group-Time Average Treatment Effects:", "x"))
  agg_lines <- c("", "Call:", "aggte(f)", "", "Reference: Callaway ...", "", "", "Overall summary:", "x")
  expect_identical(chained_result_drop_citation(agg_lines), c("", "Call:", "aggte(f)", "", "Overall summary:", "x"))
  expect_identical(chained_result_drop_citation(c("Reference: Callaway ...", "", "", "Overall summary:")),
                   c("", "Overall summary:"))
  expect_identical(chained_result_drop_citation(c("Call:", "f()", "Table")), c("Call:", "f()", "Table"))
})

test_that("the heading reports loaded and tested versions", {
  mp <- mk_mp()
  mp$chained$versions <- list(cdid = "1.0.0.0", did = "2.5.1", did_tested = "2.5.1")
  out <- paste(utils::capture.output(summary(mp)), collapse = "\n")
  expect_match(out, "cdid 1.0.0.0 result built on did 2.5.1 (tested: did 2.5.1)", fixed = TRUE)
  mp$chained$versions$did <- "2.6.0"
  out <- paste(utils::capture.output(summary(mp)), collapse = "\n")
  expect_match(out, "cdid 1.0.0.0 result built on did 2.6.0 (tested: did 2.5.1)", fixed = TRUE)
  expect_identical(summary(mp), mp)
  a <- chained_did_AGGTEobj(fx_dynamic(), mp, quote(aggte(out, type = "dynamic")))
  expect_match(paste(utils::capture.output(summary(a)), collapse = "\n"),
               "cdid 1.0.0.0 result built on did 2.6.0 (tested: did 2.5.1)", fixed = TRUE)
})

test_that("summary.cdid_MP reports units, observation correction, PT, GMM weighting and the pretest from stored records", {
  mp <- mk_mp()
  mp$chained$calendar <- as.double(1:5)
  mp$chained$frame$n_obs <- c(5L, 5L, 3L, 2L, 5L, 1L, 2L, 0L)
  mp$chained$selection <- list(method = "corollary1", canonical_assumption = "conditional_mtar",
                               requested_selection = "conditional_mtar", selection_model = "pair_logit",
                               selection_pooling = "separate", fitting_universe = "F")
  mp$chained$pt <- list(pt_assumption = "post", pt_start = NULL)
  mp$chained$covariance_handling <- list(weighting_requested = "two_step", weighting_canonical = "optimal",
                                         handler_selected = "tikhonov", handler_invoked = FALSE, handler_used = "none",
                                         trigger_reason = "well_conditioned_pd", lambda = NA_real_)
  txt_of <- function(x) paste(utils::capture.output(summary(x)), collapse = "\n")
  out <- txt_of(mp)
  expect_match(out, "Units: 8 in the unit table, 3 observed at all 5 dates, 6 with at least two observed outcomes.", fixed = TRUE)
  expect_match(out, paste0("Observation correction: observation weighting under selection = \"conditional_mtar\" ",
                           "(selection_model = \"pair_logit\", selection_pooling = \"separate\")."), fixed = TRUE)
  expect_match(out, "Parallel trends: PT-post; base_period = \"universal\".", fixed = TRUE)
  expect_match(out, paste0("GMM weighting: optimal (two-step); covariance regularization: not applied ",
                           "(moment covariance well conditioned)."), fixed = TRUE)
  mp$chained$selection <- list(method = "imputation", canonical_assumption = "conditional_mtar",
                               requested_selection = "corollary1", imputation_pooling = "within_cohort")
  mp$chained$covariance_handling <- utils::modifyList(mp$chained$covariance_handling,
    list(handler_invoked = TRUE, handler_used = "tikhonov", trigger_reason = "spectral_condition_below_trigger", lambda = 1.5e-06))
  mp$chained$pt <- list(pt_assumption = "all", pt_start = 1, no_extra_restrictions = FALSE, n_restricted_pre_parameters = 2L,
                        by_cohort = data.frame(g = c(3, 4), ref = c(2, 3), n_extra_calendar_dates = c(1L, 2L),
                                               n_restricted_pre_parameters = c(1L, 1L)))
  out <- txt_of(mp)
  expect_match(out, paste0("Observation correction: change imputation under selection = \"corollary1\" ",
                           "(imputation_pooling = \"within_cohort\"); no observation probability model is fitted."), fixed = TRUE)
  expect_match(out, paste0("covariance regularization: Tikhonov ridge applied with lambda = 1.5e-06 ",
                           "(moment covariance near singular)."), fixed = TRUE)
  expect_match(out, "Parallel trends: PT-all from common stored date 1, in addition to every PT-post restriction; base_period = \"universal\".", fixed = TRUE)
  expect_match(out, "Added pre-period restrictions: 2.", fixed = TRUE)
  expect_match(out, "Cohort 4: reference period 3; stored window dates 2; represented restrictions 1.", fixed = TRUE)
  mp$chained$pt$no_extra_restrictions <- TRUE
  mp$chained$selection <- list(method = "none", canonical_assumption = "none")
  mp$chained$covariance_handling$weighting_canonical <- "identity"
  out <- txt_of(mp)
  expect_match(out, "No additional restrictions (PT-post-equivalent model).", fixed = TRUE)
  expect_match(out, "Observation correction: none (selection = \"none\").", fixed = TRUE)
  expect_match(out, "GMM weighting: identity; covariance regularization: not applicable.", fixed = TRUE)
  # Pretrend test: the statistic and its scope are printed; tested cell IDs are not;
  # an unavailable test keeps its full reason.
  pre <- list(W = 1.234, Wpval = 0.267, cells = "pre(g=3,t=1)", df = 1L, description = "fixture: supplied externally",
              status = "computed")
  out <- txt_of(mk_mp(pretest = pre))
  expect_match(out, "Parallel-trends pre-test: W = 1.234, df = 1, p = 0.267 (Wald test of the pretrend cell).", fixed = TRUE)
  mp$chained$pretest <- list(status = "computed", W = 7.5, Wpval = 0.0576, df = 3L,
                             cells = c("pre(g=3,t=1)", "pre(g=4,t=2)", "pre(g=4,t=3)"),
                             structural_reference_cells = "reference(g=3,t=2)")
  out <- txt_of(mp)
  expect_match(out, paste0("Parallel-trends pre-test: W = 7.5, df = 3, p = 0.0576 (Wald test of all 3 pretrend cells ",
                           "jointly; 1 exact reference zero excluded)."), fixed = TRUE)
  expect_no_match(out, "pre(g=4,t=2)", fixed = TRUE)
  why <- "declared pretrend comparison pre(g=3,t=1) is unavailable (isolated_unsupported); no available-only subset is substituted"
  mp$chained$pretest <- list(status = "unavailable", reason = why)
  expect_match(txt_of(mp), paste0("Parallel-trends pre-test: unavailable; ", why, "."), fixed = TRUE)
  expect_identical(summary(mp), mp)
})

test_that("ggdid.cdid_MP pairs sparse, reordered rows by their actual (g, t)", {
  mp <- mk_mp()
  grid <- chained_mp_presentation_grid(mp)
  expect_identical(nrow(grid), 8L)                             # 2 groups x 4 reported dates
  expect_identical(grid$group, rep(c(3, 4), each = 4L))
  expect_identical(grid$time, rep(c(1, 3, 4, 5), times = 2L))
  pick <- function(g, t) grid[grid$group == g & grid$time == t, ]
  expect_identical(pick(3, 3)$att, 2.0); expect_identical(pick(3, 3)$status, "available")
  expect_identical(pick(4, 4)$att, 1.5); expect_identical(pick(4, 4)$target_id, "ATT(g=4,t=4)")
  expect_identical(pick(4, 5)$att, 3.0)
  expect_identical(pick(3, 1)$att, 0); expect_identical(pick(3, 1)$se, 0)
  expect_identical(pick(3, 4)$status, "unavailable"); expect_true(is.na(pick(3, 4)$att))
  expect_identical(pick(3, 4)$support_status, "isolated_unsupported")
  expect_identical(pick(4, 1)$status, "unrequested"); expect_identical(pick(4, 3)$status, "unrequested")
  expect_identical(pick(3, 5)$status, "unrequested")
  expect_true(all(is.na(grid$target_id[grid$status == "unrequested"])))
  expect_identical(sum(grid$status == "available"), 4L)
  p <- ggdid(mp)
  expect_true(inherits(p, "ggplot"))
  expect_identical(attr(p, "cdid_presentation_grid"), grid)
  pd <- p$data
  expect_identical(nrow(pd), 8L)
  expect_identical(pd$att[pd$group == 4 & as.character(pd$year) == "4"], 1.5)
  expect_identical(pd$att[pd$group == 3 & as.character(pd$year) == "3"], 2.0)
  expect_identical(pd$att[pd$group == 4 & as.character(pd$year) == "5"], 3.0)
  expect_true(is.na(pd$att[pd$group == 3 & as.character(pd$year) == "4"]))
  expect_true(is.na(pd$att[pd$group == 4 & as.character(pd$year) == "3"]))
  expect_identical(as.character(pd$post[pd$group == 3 & as.character(pd$year) == "1"]), "0")
  p3 <- ggdid(mp, group = 3, title = "Cohort")
  expect_true(inherits(p3, "ggplot"))
  expect_true(all(p3$data$group == 3))
})

# ------------------------------------------------------------------------------
# cdid_AGGTEobj construction
# ------------------------------------------------------------------------------

test_that("dynamic aggregation maps complete share IFs under upstream field names", {
  mp <- mk_mp(); ag <- fx_dynamic()
  acall <- quote(aggte(out, type = "dynamic", na.rm = TRUE))
  a <- chained_did_AGGTEobj(ag, mp, acall)
  expect_s3_class(a, c("cdid_AGGTEobj", "AGGTEobj"), exact = TRUE)
  expect_identical(a$type, "dynamic")
  expect_identical(a$egt, c(-2, 0, 1))
  expect_identical(a$att.egt, unname(ag$estimates))
  expect_identical(a$se.egt, unname(ag$se))
  expect_equal(a$crit.val.egt, stats::qnorm(0.975))
  expect_identical(a$overall.att, unname(ag$overall$estimate))
  expect_identical(a$overall.se, unname(ag$overall$se))
  expect_identical(names(a$inf.function), c("dynamic.inf.func.e", "dynamic.inf.func"))
  expect_identical(dim(a$inf.function$dynamic.inf.func.e), c(fx_n, 3L))
  expect_identical(dimnames(a$inf.function$dynamic.inf.func.e), list(as.character(fx_ids), ag$cell_labels))
  # complete IF (fixed + share), not the fixed-share-only IF
  expect_identical(unname(a$inf.function$dynamic.inf.func.e), unname(ag$if_full))
  expect_false(isTRUE(all.equal(unname(a$inf.function$dynamic.inf.func.e), unname(ag$if_fixed_full))))
  expect_true(any(ag$if_share_full[, 2L] != 0))
  expect_identical(unname(a$inf.function$dynamic.inf.func), ag$overall$if_full)
  expect_identical(names(a$inf.function$dynamic.inf.func), as.character(fx_ids))
  expect_identical(a$min_e, -Inf); expect_identical(a$max_e, Inf); expect_null(a$balance_e)
  expect_identical(a$call, acall)
  expect_identical(a$DIDparams$call, fx_call)
  expect_false(a$DIDparams$bstrap); expect_false(a$DIDparams$cband)
  ch <- a$chained
  expect_true(ch$na_rm)
  expect_identical(ch$cell$restricted, c(FALSE, FALSE, TRUE))
  expect_identical(ch$cell$excluded_ids[3L], "ATT(g=3,t=4)")
  expect_identical(ch$cell$n_intended, c(1L, 2L, 2L)); expect_identical(ch$cell$n_used, c(1L, 2L, 1L))
  expect_identical(ch$overall$used_bins, c("event_0", "event_1"))
  expect_identical(ch$overall$band_family, "pointwise_normal_over_certified_ids")
  expect_identical(ch$se_fixed_shares_only, ag$se_fixed_shares_only)
  expect_false(ch$selectors$applied)
  expect_identical(nobs(a), fx_n)
})

test_that("group, calendar and simple aggregations follow their type-specific conventions", {
  mp <- mk_mp()
  g <- chained_did_AGGTEobj(fx_group(), mp, quote(aggte(out, type = "group", na.rm = TRUE)))
  expect_identical(g$egt, c(3, 4))
  expect_identical(names(g$inf.function), c("selective.inf.func.g", "selective.inf.func"))
  expect_identical(unname(g$inf.function$selective.inf.func.g), unname(fx_group()$if_full))
  expect_identical(unname(g$inf.function$selective.inf.func), fx_group()$overall$if_full)
  expect_true(g$chained$overall$share_dependent)
  expect_null(g$min_e); expect_null(g$max_e); expect_null(g$balance_e)
  expect_identical(g$overall.att, fx_group()$overall$estimate)
  expect_identical(g$chained$cell$restricted, c(TRUE, FALSE))

  cf <- fx_calendar()
  k <- chained_did_AGGTEobj(cf, mp, quote(aggte(out, type = "calendar")))
  expect_identical(k$egt, c(3, 4, 5))
  expect_identical(names(k$inf.function), c("calendar.inf.func.t", "calendar.inf.func"))
  expect_true(is.na(k$att.egt[2L]) && is.na(k$se.egt[2L]))
  expect_true(all(is.na(k$inf.function$calendar.inf.func.t[, 2L])))
  expect_true(is.na(k$overall.att) && is.na(k$overall.se))
  expect_true(all(is.na(k$inf.function$calendar.inf.func)))
  expect_identical(k$chained$overall$status, "unavailable")
  expect_match(k$chained$overall$reason, "missing_inner_aggregates")
  expect_match(k$chained$cell$reason[2L], "ATT\\(g=3,t=4\\)")
  expect_identical(k$chained$status, "partially_available")

  sf <- fx_simple()
  s <- chained_did_AGGTEobj(sf, mp, quote(aggte(out, type = "simple", na.rm = TRUE)))
  expect_identical(s$type, "simple")
  expect_null(s$egt); expect_null(s$att.egt); expect_null(s$se.egt); expect_null(s$crit.val.egt)
  expect_identical(s$overall.att, unname(sf$estimates))
  expect_identical(s$overall.se, unname(sf$se))
  expect_identical(names(s$inf.function), "simple.att")
  expect_identical(unname(s$inf.function$simple.att), unname(sf$if_full[, 1L]))
  expect_true(s$chained$cell$restricted)
  expect_identical(s$chained$cell$excluded_ids, "ATT(g=3,t=4)")

  u <- chained_did_AGGTEobj(fx_unavailable_group(), mp, quote(aggte(out, type = "group")))
  expect_identical(u$egt, c(3, 4))
  expect_true(all(is.na(u$att.egt)) && is.na(u$overall.att) && is.na(u$crit.val.egt))
  expect_identical(nrow(u$chained$cell), length(u$egt))
  expect_identical(u$chained$cell$n_intended, c(1L, 1L))
  expect_identical(u$chained$status, "unavailable")
  expect_identical(u$chained$band_family, "none")
  expect_identical(u$chained$cell$reason, rep("no_declared_post_cells", 2L))
})

test_that("aggregation adapter refuses selectors it cannot have applied and unfaithful inputs", {
  mp <- mk_mp(); ag <- fx_dynamic()
  expect_error(chained_did_AGGTEobj(ag, mp, quote(aggte(out)), min_e = -1), "selector", class = "cdid_result_adapter_error")
  expect_error(chained_did_AGGTEobj(ag, mp, quote(aggte(out)), max_e = 2), "selector", class = "cdid_result_adapter_error")
  expect_error(chained_did_AGGTEobj(ag, mp, quote(aggte(out)), balance_e = 1), "selector", class = "cdid_result_adapter_error")
  bad <- ag; bad$if_full <- bad$if_fixed_full[-1L, , drop = FALSE]
  expect_error(chained_did_AGGTEobj(bad, mp, quote(aggte(out))), "n_total x M", class = "cdid_result_adapter_error")
  bad <- ag; rownames(bad$if_full) <- as.character(rev(fx_ids))
  expect_error(chained_did_AGGTEobj(bad, mp, quote(aggte(out))), "rownames", class = "cdid_result_adapter_error")
  bad <- ag; names(bad$overall$if_full) <- as.character(rev(fx_ids))
  expect_error(chained_did_AGGTEobj(bad, mp, quote(aggte(out))), "names", class = "cdid_result_adapter_error")
  aligned <- ag; rownames(aligned$if_full) <- names(aligned$overall$if_full) <- as.character(fx_ids)
  expect_identical(chained_did_AGGTEobj(aligned, mp, quote(aggte(out)))$inf.function$dynamic.inf.func.e,
                   aligned$if_full)
  bad <- ag; bad$cell_target_ids[[2L]] <- c("ATT(g=3,t=3)", "ATT(g=9,t=9)")
  expect_error(chained_did_AGGTEobj(bad, mp, quote(aggte(out))), "reported targets", class = "cdid_result_adapter_error")
  bad <- fx_calendar(); bad$if_full[, 2L] <- 0
  expect_error(chained_did_AGGTEobj(bad, mp, quote(aggte(out))), "never zero-filled", class = "cdid_result_adapter_error")
  plain <- did::MP(group = 3, t = 3, att = 1, V_analytical = NULL, se = 1, c = 1.96, inffunc = NULL, n = 8L)
  expect_error(chained_did_AGGTEobj(ag, plain, quote(aggte(out))), "cdid_MP", class = "cdid_result_adapter_error")
  expect_error(chained_did_AGGTEobj(unclass(ag), mp, quote(aggte(out))), class = "cdid_result_adapter_error")
  # Aggregation retains inherited parameter fields and records actual draws separately.
  ab <- chained_did_AGGTEobj(fx_dynamic(bootstrap = TRUE), mp, quote(aggte(out, type = "dynamic")))
  expect_false(ab$DIDparams$bstrap); expect_false(ab$DIDparams$cband)
  expect_identical(ab$DIDparams$biters, formals(did::DIDparams)$biters)
  expect_true(ab$chained$inference$bstrap); expect_true(ab$chained$inference$cband)
  expect_identical(ab$chained$inference$biters, 77L)
  expect_identical(ab$DIDparams$cdid$actual_draws, 77L)
  expect_match(paste(utils::capture.output(summary(ab)), collapse = "\n"), "Simult.", fixed = TRUE)
  expect_identical(summary(ab), ab)
  expect_identical(ab$crit.val.egt, 2.3)
  expect_false(mp$DIDparams$bstrap)
})

# ------------------------------------------------------------------------------
# cdid_AGGTEobj methods
# ------------------------------------------------------------------------------

test_that("tidy.cdid_AGGTEobj keeps did's rows/columns per type with NA convention and metadata", {
  mp <- mk_mp()
  d <- tidy(chained_did_AGGTEobj(fx_dynamic(), mp, quote(aggte(out, type = "dynamic"))))
  expect_identical(names(d)[1:11], c("type", "term", "event.time", "estimate", "std.error", "statistic", "p.value",
                                     "conf.low", "conf.high", "point.conf.low", "point.conf.high"))
  expect_identical(nrow(d), 3L)                              # no overall row for dynamic (did convention)
  expect_true(is.na(d$statistic[1L]) && !is.nan(d$statistic[1L]))   # reference-only bin: se 0 -> NA, not NaN
  expect_true(is.na(d$p.value[1L]))
  expect_equal(d$statistic[2:3], d$estimate[2:3] / d$std.error[2:3])
  expect_identical(d$label, c("event_-2", "event_0", "event_1"))
  expect_identical(d$restricted, c(FALSE, FALSE, TRUE)); expect_identical(d$excluded.ids[3L], "ATT(g=3,t=4)")
  expect_true(all(d$certified)); expect_true(all(d$band.family == "pointwise_normal_over_certified_ids"))

  g <- tidy(chained_did_AGGTEobj(fx_group(), mp, quote(aggte(out, type = "group"))))
  expect_identical(nrow(g), 3L)                              # Average + 2 groups (did convention)
  expect_identical(g$term[1L], "ATT(Average)"); expect_identical(g$group, c("Average", "3", "4"))
  expect_identical(g$label, c("overall", "group_3", "group_4"))
  expect_equal(g$statistic, g$estimate / g$std.error)
  expect_identical(g$n.used, c(2L, 1L, 2L))

  k <- tidy(chained_did_AGGTEobj(fx_calendar(), mp, quote(aggte(out, type = "calendar"))))
  expect_identical(nrow(k), 3L); expect_identical(k$time, c(3, 4, 5))
  expect_true(is.na(k$estimate[2L]) && is.na(k$statistic[2L]) && is.na(k$conf.low[2L]))
  expect_false(k$available[2L]); expect_false(k$certified[2L]); expect_true(k$certified[1L])
  expect_false(any(is.nan(k$statistic), na.rm = TRUE))

  s <- tidy(chained_did_AGGTEobj(fx_simple(), mp, quote(aggte(out, type = "simple"))))
  expect_identical(nrow(s), 1L); expect_identical(s$term, "ATT(simple average)")
  expect_true(s$restricted); expect_identical(s$excluded.ids, "ATT(g=3,t=4)")

  u <- tidy(chained_did_AGGTEobj(fx_unavailable_group(), mp, quote(aggte(out, type = "group"))))
  expect_identical(nrow(u), 3L); expect_true(all(is.na(u$estimate)))
  expect_false(any(is.nan(u$statistic), na.rm = TRUE))
})

test_that("glance/summary/print/ggdid for cdid_AGGTEobj reuse did and annotate the chained state", {
  mp <- mk_mp()
  a <- chained_did_AGGTEobj(fx_dynamic(), mp, quote(aggte(out, type = "dynamic", na.rm = TRUE)))
  gl <- glance(a)
  expect_identical(names(gl)[1:6], c("type", "nobs", "ngroup", "ntime", "control.group", "est.method"))
  expect_identical(gl$nobs, fx_n); expect_identical(gl$type, "dynamic")
  expect_identical(gl$n.bins, 3L); expect_identical(gl$n.bins.available, 3L)
  expect_identical(gl$overall.status, "available"); expect_true(gl$na.rm)
  out <- paste(utils::capture.output(summary(a)), collapse = "\n")
  expect_match(out, "cdid chained difference-in-differences (chained IPW-GMM):", fixed = TRUE)
  expect_match(out, "Ordinary estimators: Callaway", fixed = TRUE); expect_match(out, "Chained estimators: Bellego", fixed = TRUE)
  expect_match(out, "event-study/dynamic aggregation", fixed = TRUE)   # did's summary
  expect_match(out, "3 of 3 intended aggregate(s) available; na.rm = TRUE.", fixed = TRUE)
  expect_match(out, "event_1 (restricted target)", fixed = TRUE)
  expect_match(out, "Overall summary: available over its 2 intended bins; pointwise interval.", fixed = TRUE)
  expect_match(out, "Inference: analytic standard errors, clustered on 5 clusters (1 without eligible units); pointwise intervals.", fixed = TRUE)
  # Default selectors and shares estimated from the unit table are not printed.
  for (internal in c("Selectors", "Cohort shares", "Band family", "numerical certificate"))
    expect_no_match(out, internal, fixed = TRUE)
  expect_identical(utils::capture.output(print(a)), utils::capture.output(summary(a)))
  expect_identical(summary(a), a)
  outk <- paste(utils::capture.output(summary(chained_did_AGGTEobj(fx_calendar(), mp, quote(aggte(out, type = "calendar"))))), collapse = "\n")
  expect_match(outk, "calendar_4 (unavailable): missing_targets", fixed = TRUE)
  expect_match(outk, "Overall summary: unavailable: missing_inner_aggregates", fixed = TRUE)
  outs <- paste(utils::capture.output(summary(chained_did_AGGTEobj(fx_simple(), mp, quote(aggte(out, type = "simple"))))), collapse = "\n")
  expect_match(outs, "ATT_simple (restricted target)", fixed = TRUE)

  pd <- ggdid(a)
  expect_true(inherits(pd, "ggplot"))
  bins <- attr(pd, "cdid_presentation_bins")
  expect_identical(bins$key, c(-2, 0, 1)); expect_identical(bins$status, c("available", "available", "restricted"))
  expect_identical(pd$data$att, a$att.egt)
  expect_true(inherits(ggdid(chained_did_AGGTEobj(fx_group(), mp, quote(aggte(out, type = "group")))), "ggplot"))
  pk <- ggdid(chained_did_AGGTEobj(fx_calendar(), mp, quote(aggte(out, type = "calendar"))))
  expect_true(inherits(pk, "ggplot"))
  expect_identical(attr(pk, "cdid_presentation_bins")$status, c("available", "unavailable", "available"))
  # did has no simple-aggregation plot; the upstream refusal is kept
  expect_error(ggdid(chained_did_AGGTEobj(fx_simple(), mp, quote(aggte(out, type = "simple")))),
               "Plot method not available")
})

test_that("summary.cdid_AGGTEobj discloses omitted bins, nondefault selectors and share caveats", {
  mp <- mk_mp()
  txt_of <- function(x) paste(utils::capture.output(summary(x)), collapse = "\n")
  a <- chained_did_AGGTEobj(fx_dynamic(), mp, quote(aggte(out, type = "dynamic", na.rm = TRUE)))
  a$chained$overall$used_bins <- "event_0"
  a$chained$overall$excluded_bins <- "event_1"
  expect_match(txt_of(a), "Overall summary: available over 1 of its 2 intended bins (excluded: event_1); pointwise interval.",
               fixed = TRUE)
  a$chained$selectors$nondefault <- TRUE
  a$chained$selectors$min_e <- 0
  a$chained$selectors$max_e <- 1
  expect_match(txt_of(a), "Selectors: min_e = 0, max_e = 1, balance_e = NULL, applied as in did 2.5.1.", fixed = TRUE)
  a$chained$shares_known <- FALSE
  expect_match(txt_of(a), "Cohort shares: unknown for some units of the unit table; share-dependent targets are unavailable.",
               fixed = TRUE)
  a$chained$shares_known <- TRUE
  a$chained$share_law <- "first_calendar_masses"
  expect_match(txt_of(a), "Cohort shares: cohort masses at the first stored date; standard errors include their uncertainty.",
               fixed = TRUE)
  k <- chained_did_AGGTEobj(fx_calendar(), mp, quote(aggte(out, type = "calendar")))
  k$chained$selectors$nondefault <- TRUE
  k$chained$selectors$ignored <- TRUE
  k$chained$selectors$max_e <- 2
  outk <- txt_of(k)
  expect_match(outk, "Selectors: min_e = -Inf, max_e = 2, balance_e = NULL, ignored for calendar aggregation as in did 2.5.1.",
               fixed = TRUE)
  expect_match(outk, "calendar_4 (unavailable): missing_targets (na.rm = FALSE): ATT(g=3,t=4) [isolated_unsupported]", fixed = TRUE)
  # Effective inference of the aggregate bins, including the multiplier draw count.
  ab <- chained_did_AGGTEobj(fx_dynamic(bootstrap = TRUE), mp, quote(aggte(out, type = "dynamic")))
  expect_match(txt_of(ab), paste0("Inference: multiplier-bootstrap IQR standard errors (77 draws), clustered on 5 clusters ",
                                  "(1 without eligible units); simultaneous band over the available estimates."), fixed = TRUE)
})

# ------------------------------------------------------------------------------
# Ordinary did objects are unaffected
# ------------------------------------------------------------------------------

test_that("ordinary MP/AGGTEobj dispatch is unchanged and no base-class method is registered by cdid", {
  ns <- asNamespace("cdid")
  for (m in c("summary", "print", "tidy", "glance", "ggdid", "nobs")) {
    expect_false(exists(paste0(m, ".MP"), envir = ns, inherits = FALSE), info = m)
    expect_false(exists(paste0(m, ".AGGTEobj"), envir = ns, inherits = FALSE), info = m)
  }
  for (m in c("summary", "print", "tidy", "glance", "ggdid")) {
    expect_identical(environment(utils::getS3method(m, "MP")), asNamespace("did"), info = m)
    expect_identical(environment(utils::getS3method(m, "AGGTEobj")), asNamespace("did"), info = m)
    expect_identical(environment(utils::getS3method(m, "cdid_MP")), ns, info = m)
    expect_identical(environment(utils::getS3method(m, "cdid_AGGTEobj")), ns, info = m)
  }
  dp <- did::DIDparams(yname = "y", tname = "t", idname = "id", gname = "g", data = NULL,
                       control_group = "nevertreated", bstrap = FALSE, cband = FALSE, est_method = "dr",
                       true_repeated_cross_sections = FALSE, n = 8L, nG = 2L, nT = 3L, tlist = 1:3, glist = 2:3,
                       call = quote(att_gt()))
  plain <- did::MP(group = c(2, 2, 3), t = c(2, 3, 3), att = c(1, 0, 2), V_analytical = diag(3),
                   se = c(1, 0, 1), c = 1.96, inffunc = matrix(0, 8L, 3L), n = 8L, DIDparams = dp)
  expect_identical(class(plain), "MP")
  td <- tidy(plain)
  expect_identical(names(td), upstream_tidy_mp_cols)          # no appended columns, upstream NaN kept
  expect_true(is.nan(td$statistic[2L]))
  expect_identical(names(glance(plain)), c("nobs", "ngroup", "ntime", "control.group", "est.method"))
  out <- paste(utils::capture.output(summary(plain)), collapse = "\n")
  expect_no_match(out, "chained", fixed = TRUE)
  # did's own display of its own objects is untouched, including its citation line.
  expect_identical(utils::capture.output(summary(plain)),
                   utils::capture.output(utils::getS3method("summary", "MP", envir = asNamespace("did"))(plain)))
  # The unchanged upstream method rebuilds a 2 x 2 rectangle from these
  # sparse three rows and refuses. Preserve that ordinary-object behavior.
  expect_error(ggdid(plain), "replacement has 3 rows, data has 4", fixed = TRUE)
  pa <- did::AGGTEobj(overall.att = 1, overall.se = 0.5, type = "dynamic", egt = c(-1, 0), att.egt = c(0, 1),
                      se.egt = c(0.2, 0.3), crit.val.egt = 1.96, DIDparams = dp)
  expect_identical(class(pa), "AGGTEobj")
  expect_identical(nrow(tidy(pa)), 2L)
  expect_identical(ncol(tidy(pa)), 11L)
  expect_no_match(paste(utils::capture.output(summary(pa)), collapse = "\n"), "chained", fixed = TRUE)
  expect_identical(utils::capture.output(summary(pa)),
                   utils::capture.output(utils::getS3method("summary", "AGGTEobj", envir = asNamespace("did"))(pa)))
})

# ------------------------------------------------------------------------------
# Real-shape smoke case: already computed N200 fields (no refit)
# ------------------------------------------------------------------------------

test_that("N200 already computed fields map to cdid_MP / cdid_AGGTEobj and the methods run", {
  fx <- readRDS(testthat::test_path("fixtures", "cp02_result_adapter_N200.rds"))
  rep <- fx$reporting; fr <- fx$frame
  # independent identity evidence: the share IF is positive exactly for frame cohort 3 rows
  expect_identical(which(fx$aggregations$dynamic$if_shares_full[, "3"] > 0), which(fr$G == 3))
  expect_true(all(fr$eligible_ge2))
  expect_false(rep$inference$clustered)
  # the caller normalizes the frame's public never-treated code 0 to internal Inf explicitly
  units <- data.frame(id = fr$id, cohort = ifelse(fr$G == 0, Inf, fr$G), weight = fr$sample_weight,
                      eligible = fr$eligible_ge2)
  units <- units[order(-units$id), ]                            # permuted on purpose
  params <- list(yname = "y", tname = "period", idname = "id", gname = "G", xformla = ~1,
                 control_group = "nevertreated", anticipation = 0, alp = 0.05, tlist = fx$calendar,
                 glist = sort(unique(units$cohort[is.finite(units$cohort)])), est_method = "ipw",
                 base_period = rep$base_period, bstrap = FALSE, cband = FALSE)
  call <- quote(att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = d, mode = "chained"))
  mp <- chained_did_MP(rep, fr$id, call, params, units = units, support = fx$support,
                       basis = fx$basis, band_families = NULL)
  expect_s3_class(mp, c("cdid_MP", "MP"), exact = TRUE)
  expect_identical(mp$n, 200L)
  expect_identical(mp$group, rep$target_meta$g); expect_identical(mp$t, rep$target_meta$time)
  expect_identical(rownames(mp$inffunc), as.character(fr$id))
  expect_identical(unname(mp$V_analytical), unname(rep$inference$Omega))
  expect_equal(unname(mp$V_analytical), unname(rep$vcov) * 200)
  expect_equal(mp$se, unname(sqrt(diag(mp$V_analytical) / mp$n)))
  expect_identical(mp$chained$frame$id, fr$id)
  expect_identical(mp$chained$frame$cohort, ifelse(fr$G == 0, Inf, fr$G))
  expect_identical(mp$DIDparams$glist, c(3, 4, 5)); expect_identical(mp$DIDparams$nT, 5L)
  expect_false(mp$DIDparams$bstrap)
  td <- tidy(mp)
  expect_identical(nrow(td), 4L)
  expect_identical(td$role, c("pre", "post", "post", "post"))
  expect_equal(td$statistic, unname(rep$z_stat)); expect_equal(td$p.value, unname(rep$p_value))
  expect_identical(glance(mp)$n.available, 4L)
  expect_match(paste(utils::capture.output(summary(mp)), collapse = "\n"), "4 of 4 requested cells available", fixed = TRUE)
  p <- ggdid(mp)
  expect_true(inherits(p, "ggplot"))
  expect_identical(attr(p, "cdid_presentation_grid")$status, rep("available", 4L))
  for (ty in c("simple", "group", "dynamic", "calendar")) {
    a <- chained_did_AGGTEobj(fx$aggregations[[ty]], mp, as.call(list(quote(aggte), quote(out), type = ty)))
    expect_s3_class(a, c("cdid_AGGTEobj", "AGGTEobj"), exact = TRUE)
    expect_identical(a$type, ty)
    expect_true(all(vapply(a$inf.function, NROW, integer(1L)) == 200L))
    expect_equal(a$overall.att, unname(fx$aggregations[[ty]]$overall$estimate))
    tda <- tidy(a)
    expect_true(is.data.frame(tda) && nrow(tda) >= 1L)
    expect_false(any(is.nan(tda$statistic), na.rm = TRUE))
    expect_true(is.data.frame(glance(a)))
    expect_true(length(utils::capture.output(summary(a))) > 5L)
    if (ty != "simple") expect_true(inherits(ggdid(a), "ggplot"))
  }
  dyn <- chained_did_AGGTEobj(fx$aggregations$dynamic, mp, quote(aggte(out, type = "dynamic")))
  expect_identical(dyn$egt, c(-1, 0, 1, 2))
  expect_identical(names(dyn$inf.function), c("dynamic.inf.func.e", "dynamic.inf.func"))
  expect_identical(unname(dyn$inf.function$dynamic.inf.func.e), unname(fx$aggregations$dynamic$if_full))
  grp <- chained_did_AGGTEobj(fx$aggregations$group, mp, quote(aggte(out, type = "group")))
  expect_identical(grp$egt, 3)
  expect_identical(nrow(tidy(grp)), 2L)
})


test_that("unavailable summaries preserve substantive reasons and omit missing reason suffixes", {
  obj <- mk_mp()
  txt <- function(x) paste(capture.output(summary(x)), collapse = "\n")
  expect_match(txt(obj), "ATT(g=3,t=4) [isolated_unsupported]: post vertex isolated in the retained catalog", fixed = TRUE)
  for (reason in c(NA_character_, "", "   ")) {
    obj$chained$support_reason[!obj$chained$available] <- reason
    output <- txt(obj)
    expect_match(output, "ATT(g=3,t=4) [isolated_unsupported]", fixed = TRUE)
    expect_no_match(output, "[isolated_unsupported]:", fixed = TRUE)
  }
})

test_that("stored worker and certificate diagnostics are retained without normal printing", {
  obj <- mk_mp()
  obj$chained$request <- list(execution = list(workers = 2L, route = "psock"))
  obj$chained$public <- list(numerical_policy = list(decision = "D20", relative_accuracy = 1e-6,
                                                   override_recorded = FALSE))
  saved <- obj
  for (method in list(summary, print)) {
    output <- paste(capture.output(result <- withVisible(method(obj))), collapse = "\n")
    expect_false(result$visible)
    expect_identical(result$value, saved)
    expect_identical(obj, saved)
    for (needle in c("Multiplier bootstrap workers:", "Numerical acceptance:", "certified error", "Resource profile:"))
      expect_no_match(output, needle, fixed = TRUE)
  }
})

test_that("summary and print invisibly return unchanged cell and aggregate objects", {
  mp <- mk_mp()
  agg <- chained_did_AGGTEobj(fx_dynamic(), mp, quote(aggte(out, type = "dynamic", na.rm = TRUE)))
  for (obj in list(mp, agg)) {
    saved <- obj
    summary_text <- capture.output(s <- withVisible(summary(obj)))
    print_text <- capture.output(p <- withVisible(print(obj)))
    expect_false(s$visible)
    expect_false(p$visible)
    expect_identical(s$value, saved)
    expect_identical(p$value, saved)
    expect_identical(obj, saved)
    expect_identical(summary_text, print_text)
  }
})


test_that("each paper citation appears exactly once before the call", {
  mp <- mk_mp()
  ag <- chained_did_AGGTEobj(fx_dynamic(), mp, quote(aggte(out, type = "dynamic")))
  for (obj in list(mp, ag)) {
    lines <- capture.output(summary(obj))
    txt <- paste(lines, collapse = "\n")
    for (doi in c("10.1016/j.jeconom.2024.105783", "10.1016/j.jeconom.2020.12.001")) {
      expect_length(regmatches(txt, gregexpr(doi, txt, fixed = TRUE))[[1L]], 1L)
    }
    expect_lt(grep("Ordinary estimators:", lines), grep("^Call:", lines))
    expect_length(grep("^Reference:", lines), 0L)
  }
})
