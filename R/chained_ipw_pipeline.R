# Common direct-precomparison request construction; no model fitting here.
# One union of both base-period reporting plans: descriptors are the distinct
# (g, s, t, ordered K) comparisons of the varying and universal families in
# first-use order; every selection keeps its family, sign, reference-zero or
# unavailable status. CP08 scaling action 3 builds it once, before the
# covariance/replay/report reservations, and reuses it in every consumer. The
# exact identity of a descriptor is its binary64 values (hex keys, -0 = 0,
# as identical()); each clean-control set is evaluated once per date.
ipw_auxiliary_spec <- function(ledger, calendar, cohorts, anticipation, control_group) {
  cal <- sort(unique(as.double(calendar))); coh <- sort(unique(as.double(cohorts))); ant <- as.double(anticipation)
  groups <- sort(unique(as.double(ledger$g)))
  refs <- vapply(groups, function(g) tryCatch(calendar_reference_date(g, cal, ant), cdid_error = function(e) NA_real_), numeric(1L))
  descriptors <- list(); selections <- list(varying = list(), universal = list())
  push <- function(bp, rec) selections[[bp]][[length(selections[[bp]]) + 1L]] <<- rec
  seen <- new.env(hash = TRUE, parent = emptyenv())
  clean_at <- new.env(hash = TRUE, parent = emptyenv())
  universe_txt <- NULL
  for (bp in c("varying", "universal")) for (gi in seq_along(groups)) {
    g <- groups[gi]; ref <- refs[gi]
    times <- if (bp == "varying") cal[-1L] else cal
    for (t in times[times < g]) {
      if (is.na(ref)) {
        push(bp, list(g = g, t = t, ref = NA_real_, from = NA_real_, index = NA_integer_, sign = NA_real_, K = numeric(),
          status = "unavailable", reason = paste0("no_admissible_reference: no stored calendar date t < g - anticipation for cohort g=", fmt_design_num(g))))
        next
      }
      from <- if (bp == "varying") cal[match(t, cal) - 1L] else ref
      if (bp == "universal" && t == ref) {
        push(bp, list(g = g, t = t, ref = ref, from = from, index = 0L, sign = 0.0, K = numeric(), status = "reference_zero", reason = NA_character_))
        next
      }
      K <- if (identical(control_group, "nevertreated")) coh[is.infinite(coh)] else {
        u <- max(t, from); uk <- sprintf("%a", u + 0)
        cl <- clean_at[[uk]]
        if (is.null(cl)) {
          cl <- vapply(coh, function(k) exact_cmp_sum(u, ant, k) < 0L, logical(1L))
          assign(uk, cl, envir = clean_at)
        }
        coh[coh != g & cl]
      }
      if (!length(K)) {
        if (is.null(universe_txt)) universe_txt <- paste(fmt_design_num(coh), collapse = ",")
        push(bp, list(g = g, t = t, ref = ref, from = from, index = NA_integer_, sign = NA_real_, K = numeric(),
          status = "unavailable", reason = paste0("no_clean_control_cohort: no cohort of the universe {",
            universe_txt, "} other than g=", fmt_design_num(g), " is untreated at date ",
            fmt_design_num(max(t, from)), " with anticipation ", fmt_design_num(ant), " (comparison dates ", fmt_design_num(min(t, from)), "->", fmt_design_num(max(t, from)), ")")))
        next
      }
      d <- list(g = g, s = min(t, from), t = max(t, from), K = as.double(K))
      key <- paste(sprintf("%a", c(d$g, d$s, d$t, d$K) + 0), collapse = ",")
      hit <- seen[[key]]
      if (is.null(hit)) {
        descriptors[[length(descriptors) + 1L]] <- d; hit <- length(descriptors)
        assign(key, hit, envir = seen)
      }
      push(bp, list(g = g, t = t, ref = ref, from = from, index = hit, sign = if (t > from) 1.0 else -1.0, K = as.double(K), status = "pending", reason = NA_character_))
    }
  }
  list(descriptors = descriptors, selections = selections)
}

# Prospective envelope of `ipw_auxiliary_spec` from the request's actual
# dimensions: J target cohorts, T dates and C cohorts give at most S =
# J (2T - 1) selections of both families (ten fields plus at most C controls
# each) and at most S descriptors (four fields plus C controls); every
# selection may carry one reason naming the cohort universe (at most 25 C + 256
# bytes) and one transient hex key (24 (C + 4) bytes); one cleanliness vector
# per date. Strings use the pipeline payload convention (cp04_payload: half a
# slot per byte) plus eight header slots each. Charged before building.
ipw_auxiliary_spec_envelope <- function(n_targets, n_dates, n_cohorts) {
  J <- as.double(n_targets); T <- as.double(n_dates); C <- as.double(n_cohorts)
  S <- ipw_owner_safe_product(J, max(0, 2 * T - 1), label = "auxiliary selections of both families")
  reason_slots <- ipw_owner_safe_product(S, ceiling((25 * C + 256) / 2) + 8, label = "auxiliary unavailable reasons")
  key_slots <- ipw_owner_safe_product(S, ceiling(24 * (C + 4) / 2) + 8, label = "auxiliary descriptor keys")
  retained <- ipw_owner_safe_sum(c(ipw_owner_safe_product(S, 2 * C + 14, label = "auxiliary selections and descriptors"),
    reason_slots, ipw_owner_safe_product(T, C + 8, label = "cleanliness vectors")), label = "auxiliary plan payload")
  scratch <- ipw_owner_safe_sum(c(key_slots, 64 * (C + 8)), label = "auxiliary plan scratch")
  work <- ipw_owner_safe_sum(c(ipw_owner_safe_product(64, ipw_owner_safe_sum(c(ipw_owner_safe_product(J, T), ipw_owner_safe_product(T, C)))),
    ipw_owner_safe_product(16, ipw_owner_safe_sum(c(ipw_owner_safe_product(S, C + 8), reason_slots, key_slots)))),
    label = "auxiliary plan work")
  list(selections_bound = S, retained_entries = retained, scratch_entries = scratch, work_ops = work)
}

# ==============================================================================
# Chained IPW-GMM CP02 Pipeline Integration
# ==============================================================================
# Connects:
#   1. Sparse panel & covariates/weights frame ownership (`R/chained_ipw_owner.R`)
#   2. D10/D13 comparison catalog & structural target qualification (`R/chained_pipeline_metadata.R`)
#   3. Shared fitted nuisances with chunked zero-outcome mask verification (`R/chained_ipw_fitted_view.R`)
#   4. Compact point link summaries & two-pass streamed moment covariance (`R/chained_ipw_point_view.R`)
#   5. Compact structural equivalence class certification & redundancy reduction (`R/chained_ipw_shared.R`)
#   6. Chained GMM solver & numerical guard (`R/chained_gmm_core.R`)
#   7. Complete outward interval numerical transfer:
#      - Nuisance root certification (`bind_fitted_nuisance_to_problem_view`)
#      - Link moment/loading/IF enclosure (`bind_point_link_to_enclosure_view`)
#      - Moment covariance transfer (`certify_ipw_covariance_transfer_view`)
#      - Theorem 12 (`thm:gmm`) changing-weight GMM target/loading/SE/covariance transfer
#   8. Full-frame IF transport, `base_period` reporting map, degree-0 cohort-share
#      aggregation with cross terms (`docs/theory/cp02_population_aggregation_if.tex`)
#   9. Analytical and one-way cluster inference + multiplier bootstrap (`R/chained_multiplier.R`)
# ==============================================================================

cdid_ipw_pipeline_abort <- function(message, class = character(), diagnostics = NULL) {
  stop(structure(
    class = unique(c(class, "cdid_ipw_pipeline_error", "cdid_error", "error", "condition")),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

# Inherited did `print_details`: TRUE reports stage and comparison progress
# through message(); FALSE returns NULL and nothing is emitted. Reporting reads
# only counts already computed; it changes no data, RNG, estimate or bound.
chained_progress_emitter <- function(print_details) {
  if (!isTRUE(print_details)) return(NULL)
  function(...) message("[cdid] ", ...)
}

# Per-comparison progress of a loop over M rows, at most about ten lines.
chained_progress_rows <- function(progress, what) {
  if (is.null(progress)) return(NULL)
  function(m, M, n_views = NULL) {
    step <- max(1, ceiling(M / 10))
    if (m == M || m %% step == 0)
      progress(what, ": ", m, " of ", M, " comparisons", if (!is.null(n_views)) paste0(" (", n_views, " distinct nuisance views)") else "")
  }
}

# compute_inffunc = FALSE (did 2.5.1 convention): point estimates only. The
# reported cells keep their estimates and availability; every standard error,
# interval and band is NA, no covariance, multiplier draw or influence function
# is formed, and the critical value is the pointwise normal constant.
ipw_pipeline_point_only_inference <- function(estimates, available, n, cluster_full, eligible_full, alp) {
  ids <- names(estimates); J <- length(estimates)
  est <- stats::setNames(ifelse(available, as.double(estimates), NA_real_), ids)
  na_v <- stats::setNames(rep(NA_real_, J), ids)
  u_cl <- if (is.null(cluster_full)) seq_len(n) else unique(cluster_full)
  map <- if (is.null(cluster_full)) seq_len(n) else match(cluster_full, u_cl)
  if (is.null(eligible_full)) eligible_full <- rep(TRUE, n)
  list(status = "inference_not_computed",
    reason = "compute_inffunc = FALSE: point estimates only (no influence functions, standard errors, bands, pre-test or aggregation)",
    estimates = est, se = na_v, analytical_se = na_v, z_stat = na_v, p_value = na_v, ci_lower = na_v, ci_upper = na_v,
    cband_crit_val = stats::qnorm(1 - alp / 2), cband_lower = na_v, cband_upper = na_v,
    vcov = NULL, Omega = NULL, n_total = n, n_clusters = length(u_cl),
    n_zero_eligible_clusters = sum(!seq_along(u_cl) %in% unique(map[eligible_full])),
    clustered = !is.null(cluster_full), multiplier = NULL, rng = NULL,
    bootstrap = FALSE, simultaneous = FALSE, flag_law = "compute_inffunc_false_point_estimates_only",
    band_family = "none", n_draws = NA_integer_,
    support = list(status = if (all(available)) "all_requested_supported" else "supported_subvector",
      available_ids = ids[available], unavailable_ids = ids[!available], n_available = sum(available),
      n_requested = J, whole_request_complete = all(available)),
    numerical_certificate = NULL)
}

ipw_pipeline_control <- function(
  owner_control,
  fitted_control,
  point_control,
  glm_control,
  certificate_control,
  problem_control,
  link_control,
  covariance_control,
  class_bounds,
  gmm_bounds,
  gmm_accuracy,
  catalog_bounds,
  multiplier_budget,
  inference_accuracy
) {
  if (missing(owner_control) || missing(fitted_control) || missing(point_control) ||
      missing(glm_control) || missing(certificate_control) || missing(problem_control) ||
      missing(link_control) || missing(covariance_control) || missing(class_bounds) ||
      missing(gmm_bounds) || missing(gmm_accuracy) || missing(catalog_bounds) ||
      missing(multiplier_budget) || missing(inference_accuracy)) {
    cdid_ipw_pipeline_abort(
      "All fourteen control structures (`owner_control`, `fitted_control`, `point_control`, `glm_control`, `certificate_control`, `problem_control`, `link_control`, `covariance_control`, `class_bounds`, `gmm_bounds`, `gmm_accuracy`, `catalog_bounds`, `multiplier_budget`, `inference_accuracy`) must be supplied explicitly; no hidden defaults are adopted.",
      c("cdid_ipw_pipeline_control_error", "cdid_ipw_control_error")
    )
  }
  structure(
    list(
      owner_control = validate_ipw_owner_control(owner_control),
      fitted_control = validate_ipw_fitted_view_control(fitted_control),
      point_control = validate_ipw_point_view_control(point_control),
      glm_control = validate_ipw_glm_control(glm_control),
      certificate_control = validate_logit_certificate_control(certificate_control),
      problem_control = validate_ipw_nuisance_problem_control(problem_control),
      link_control = validate_ipw_link_enclosure_control(link_control),
      covariance_control = validate_ipw_covariance_transfer_control(covariance_control),
      class_bounds = validate_ipw_certificate_bounds(class_bounds),
      gmm_bounds = validate_gmm_bounds_control(gmm_bounds),
      gmm_accuracy = gmm_accuracy,
      catalog_bounds = validate_catalog_bounds(catalog_bounds),
      multiplier_budget = validate_multiplier_budget(multiplier_budget),
      inference_accuracy = validate_ipw_pipeline_inference_accuracy(inference_accuracy)
    ),
    class = "cdid_ipw_pipeline_control"
  )
}

validate_ipw_pipeline_control <- function(control, weighting = "optimal") {
  if (missing(control) || is.null(control) || !is.list(control)) {
    cdid_ipw_pipeline_abort(
      "Explicit `control` constructed via `ipw_pipeline_control(...)` is required.",
      c("cdid_ipw_pipeline_control_error", "cdid_ipw_control_error")
    )
  }
  req <- c(
    "owner_control", "fitted_control", "point_control", "glm_control",
    "certificate_control", "problem_control", "link_control", "covariance_control",
    "class_bounds", "gmm_bounds", "gmm_accuracy", "catalog_bounds",
    "multiplier_budget", "inference_accuracy"
  )
  if (!all(req %in% names(control)) || anyDuplicated(names(control)) ||
      length(setdiff(names(control), req)) > 0L) {
    cdid_ipw_pipeline_abort(
      paste0("`control` must contain exactly: ", paste(req, collapse = ", "), "."),
      c("cdid_ipw_pipeline_control_error", "cdid_ipw_control_error")
    )
  }
  acc_val <- validate_gmm_accuracy_control(control$gmm_accuracy, weighting)
  ipw_pipeline_control(
    owner_control = control$owner_control,
    fitted_control = control$fitted_control,
    point_control = control$point_control,
    glm_control = control$glm_control,
    certificate_control = control$certificate_control,
    problem_control = control$problem_control,
    link_control = control$link_control,
    covariance_control = control$covariance_control,
    class_bounds = control$class_bounds,
    gmm_bounds = control$gmm_bounds,
    gmm_accuracy = acc_val,
    catalog_bounds = control$catalog_bounds,
    multiplier_budget = control$multiplier_budget,
    inference_accuracy = control$inference_accuracy
  )
}

# Small composition helpers use the accepted outward binary64 primitives. Bounds
# concern reachable numeric payload, not R heap or process RSS.
#
# Two inference-accuracy forms are accepted (October 16, 2026, D20 binding):
#   (a) the existing private ABSOLUTE form: exactly the five positive finite
#       absolute budgets (unchanged behaviour for every existing caller);
#   (b) the public SCALE-AWARE form (D20): `max_relative_error_lower_se`, a
#       positive finite level (public default 1e-6, explicit override recorded by
#       the caller), optionally together with the five absolute fields. Under (b)
#       every requested estimate, SE and interval endpoint must have a certified
#       error at most level x (certified positive LOWER endpoint of the reported
#       SE interval); covariance entries use the product of the two lower SE
#       scales; multiplier raw draws use their own dimension n/sqrt(G) x the
#       analytic lower SE scale (diagnostic only); known structural reference
#       operators with singleton-zero input enclosures are identities and must carry exactly zero error, never a ratio.
#       A non-positive or unresolved scale refuses. Absolute fields supplied
#       together with the level are still enforced (never silently weakened).
ipw_pipeline_absolute_accuracy_names <- c("max_estimate_error_abs", "max_se_error_abs", "max_vcov_error_abs",
                                          "max_multiplier_error_abs", "max_band_error_abs")
ipw_pipeline_relative_accuracy_name <- "max_relative_error_lower_se"

validate_ipw_pipeline_inference_accuracy <- function(x) {
  if (is.null(x)) return(NULL) # D23: no per-fit error certification requested.
  req <- ipw_pipeline_absolute_accuracy_names
  rel <- ipw_pipeline_relative_accuracy_name
  if (!is.list(x) || is.object(x) || is.null(names(x)) || anyDuplicated(names(x)) ||
      !all(names(x) %in% c(req, rel))) {
    cdid_ipw_pipeline_abort(paste0("Explicit inference-error budgets are required: either the five positive finite absolute ",
      "budgets (", paste(req, collapse = ", "), ") or the scale-aware level `", rel, "` (optionally with absolute budgets)."),
      "cdid_ipw_pipeline_control_error")
  }
  has_rel <- rel %in% names(x)
  abs_present <- intersect(req, names(x))
  if (!has_rel && !identical(sort(abs_present), sort(req))) {
    cdid_ipw_pipeline_abort("Explicit positive finite inference-error budgets are required.",
                            "cdid_ipw_pipeline_control_error")
  }
  if (any(!vapply(x[abs_present], gmm_scalar_numeric, logical(1L), positive = TRUE))) {
    cdid_ipw_pipeline_abort("Explicit positive finite inference-error budgets are required.",
                            "cdid_ipw_pipeline_control_error")
  }
  if (has_rel && (!gmm_scalar_numeric(x[[rel]], positive = TRUE) || x[[rel]] >= 1)) {
    cdid_ipw_pipeline_abort(paste0("`", rel, "` must be a positive finite level below 1 (D20 public default 1e-6)."),
                            "cdid_ipw_pipeline_control_error")
  }
  x
}

ipw_pipeline_accuracy_is_relative <- function(accuracy) {
  is.list(accuracy) && ipw_pipeline_relative_accuracy_name %in% names(accuracy)
}

# Directed interval matrix product on the native backend (October 16, 2026,
# measured critical path: 98% of the sampled B1000 post-fit time was this
# product inside sign_batch_fn). The kernel `C_gmm_iv_mat_product` evaluates the
# historical scalar recurrence in the same i/j/k order with the same accepted
# binary64 primitives (gmm_iv_mul_c, gmm_iv_add_c), the same zero shortcuts and
# the same typed failures; the frozen pure-R recurrence is retained only as an
# independent test reference. No bound is added here: the producer preflight
# already limits these dimensions, and the kernel keeps the exact-integer
# allocation guards of the other native entries.
ipw_pipeline_interval_product <- function(xlo, xhi, ylo, yhi = ylo) {
  gmm_iv_mat_product(xlo, xhi, ylo, yhi)
}

ipw_pipeline_interval_mean <- function(lo, hi) {
  ans <- matrix(0.0, 2L, ncol(lo))
  for (j in seq_len(ncol(lo))) {
    z <- gmm_iv(0.0)
    for (i in seq_len(nrow(lo))) z <- gmm_iv_add(z, c(lo[i,j], hi[i,j]))
    ans[,j] <- gmm_iv_div(z, gmm_iv(as.double(nrow(lo))))
  }
  ans
}

ipw_pipeline_interval_center <- function(lo, hi) {
  mu <- ipw_pipeline_interval_mean(lo, hi)
  for (j in seq_len(ncol(lo))) for (i in seq_len(nrow(lo))) {
    z <- gmm_iv_sub(c(lo[i,j], hi[i,j]), mu[,j])
    lo[i,j] <- z[1L]; hi[i,j] <- z[2L]
  }
  list(lo = lo, hi = hi)
}

ipw_pipeline_interval_around <- function(x, error) {
  lo <- hi <- as.matrix(x)
  for (j in seq_len(ncol(lo))) for (i in seq_len(nrow(lo))) {
    z <- gmm_iv_add(gmm_iv(lo[i,j]), c(-error[j], error[j]))
    lo[i,j] <- z[1L]; hi[i,j] <- z[2L]
  }
  list(lo = lo, hi = hi)
}

ipw_pipeline_interval_distance <- function(x, lo, hi) {
  out <- numeric(length(x))
  for (k in seq_along(x)) {
    a <- gmm_iv_sub(gmm_iv(as.double(x[k])), c(lo[k], hi[k]))
    out[k] <- max(abs(a))
  }
  dim(out) <- dim(x)
  out
}

# Theorem 12 supplies RMS errors. A row's error is at most sqrt(n_A) times
# its RMS error. This conservative conversion permits the existing finite
# interval arithmetic to follow every reporting/aggregation/cluster/sign map.
ipw_pipeline_target_enclosure <- function(transfer, gmm_fit, n_total, a_to_full) {
  gt <- transfer$gmm_transfer
  n_a <- length(a_to_full)
  e <- vapply(gt$R_c_cert, function(z)
    gmm_iv_mul(gmm_iv_sqrt(gmm_iv(as.double(n_a))), gmm_iv(z))[2L], numeric(1L))
  # Complete empirical IFs have exact sample mean zero at the certified score
  # root: normalized outcome/rate terms center exactly and the score sums to
  # zero (weighted-pair IF proof). Theorem 12's B is centered, so first center
  # the stored projected rows with directed arithmetic. An RMS bound alone
  # must never be applied to an unconstrained uncentered row offset.
  stored <- gmm_fit$if_projection$materialized_target_if
  rows <- ipw_pipeline_interval_center(stored,stored)
  for (j in seq_len(ncol(stored))) for(i in seq_len(nrow(stored))) {
    z <- gmm_iv_add(c(rows$lo[i,j],rows$hi[i,j]),c(-e[j],e[j]))
    rows$lo[i,j] <- z[1L];rows$hi[i,j] <- z[2L]
  }
  lo <- hi <- matrix(0.0, n_total, ncol(rows$lo))
  scale <- gmm_iv_div(gmm_iv(as.double(n_total)), gmm_iv(as.double(n_a)))
  for (j in seq_len(ncol(lo))) for (i in seq_len(n_a)) {
    z <- gmm_iv_mul(c(rows$lo[i,j], rows$hi[i,j]), scale)
    lo[a_to_full[i],j] <- z[1L]; hi[a_to_full[i],j] <- z[2L]
  }
  # D5: the compact error vector is aligned with the estimable (loading) columns
  # only; the full target table also carries NA rows for unavailable levels.
  est_compact <- gmm_fit$targets$estimate[gmm_fit$targets$estimable]
  if (length(est_compact) != length(gt$total_target_error_cert) || length(est_compact) != ncol(stored) || anyNA(est_compact)) {
    cdid_ipw_pipeline_abort("Compact target estimates, error bounds and stored IF columns are misaligned.", "cdid_ipw_pipeline_integrity_error")
  }
  theta <- ipw_pipeline_interval_around(matrix(est_compact, nrow = 1L),
                                         gt$total_target_error_cert)
  list(theta = theta, rows = list(lo = lo, hi = hi))
}

ipw_pipeline_inference_intervals <- function(enclosure, cluster_full) {
  lo <- enclosure$rows$lo; hi <- enclosure$rows$hi
  n <- nrow(lo); J <- ncol(lo)
  map <- if (is.null(cluster_full)) seq_len(n) else match(cluster_full, unique(cluster_full))
  G <- max(map)
  clo <- chi <- matrix(0.0, G, J)
  for (i in seq_len(n)) for (j in seq_len(J)) {
    z <- gmm_iv_add(c(clo[map[i],j], chi[map[i],j]), c(lo[i,j], hi[i,j]))
    clo[map[i],j] <- z[1L]; chi[map[i],j] <- z[2L]
  }
  centered <- ipw_pipeline_interval_center(clo, chi)
  # The retained nonempty GMM can have zero requested loading columns (D5).
  # Its target covariance is the exact empty 0x0 matrix; no kernel arithmetic
  # exists for that product. The native entry's positive-dimension admission
  # remains unchanged, and all preceding nuisance/GMM/transfer checks still run.
  V <- if (J == 0L) list(lo = matrix(0.0, 0L, 0L), hi = matrix(0.0, 0L, 0L)) else
    ipw_pipeline_interval_product(t(centered$lo), t(centered$hi), centered$lo, centered$hi)
  denom <- gmm_iv_mul(gmm_iv(as.double(n)), gmm_iv(as.double(n)))
  for (k in seq_along(V$lo)) {
    z <- gmm_iv_div(c(V$lo[k], V$hi[k]), denom)
    V$lo[k] <- z[1L]; V$hi[k] <- z[2L]
  }
  se <- matrix(0.0, 2L, J)
  for (j in seq_len(J)) se[,j] <- gmm_iv_sqrt(c(max(0.0, V$lo[j,j]), max(0.0, V$hi[j,j])))
  list(theta = enclosure$theta, rows = enclosure$rows, cluster = centered, vcov = V, se = se)
}

ipw_pipeline_enforce_inference <- function(cert, estimates, se, vcov, accuracy,
                                            band_lower, band_upper) {
  err_est <- ipw_pipeline_interval_distance(estimates, cert$theta$lo, cert$theta$hi)
  err_se <- ipw_pipeline_interval_distance(se, cert$se[1L,], cert$se[2L,])
  err_v <- ipw_pipeline_interval_distance(vcov, cert$vcov$lo, cert$vcov$hi)
  err_band <- max(ipw_pipeline_interval_distance(band_lower, cert$band$lo, cert$band$hi),
                  ipw_pipeline_interval_distance(band_upper, cert$band_upper$lo, cert$band_upper$hi))
  vals <- c(max_estimate_error_abs = max(err_est), max_se_error_abs = max(err_se),
            max_vcov_error_abs = max(err_v), max_multiplier_error_abs = cert$multiplier_error,
            max_band_error_abs = err_band)
  abs_names <- intersect(ipw_pipeline_absolute_accuracy_names, names(accuracy))
  if (any(!is.finite(vals)) || any(vals[abs_names] > unlist(accuracy[abs_names]))) {
    cdid_ipw_pipeline_abort("Reported inference exceeds explicit numerical-error budgets.",
                            "cdid_ipw_pipeline_accuracy_error", list(errors = vals, budgets = accuracy))
  }
  cert$errors <- vals
  if (ipw_pipeline_accuracy_is_relative(accuracy)) {
    cert <- ipw_pipeline_enforce_d20(cert, estimates, se, vcov, accuracy[[ipw_pipeline_relative_accuracy_name]],
      err_est = err_est, err_se = err_se, err_v = err_v, band_lower = band_lower, band_upper = band_upper)
  }
  cert$certified <- TRUE
  cert$scope <- paste0("Conditional finite arithmetic through reported rows, covariance, SE, fixed/generated signs and bands; no sampling or MC acceptance.",
    if (!is.null(cert$d20)) paste0(" D20: every certified error is at most ", format(cert$d20$level, digits = 6),
      " times the certified lower statistical SE scale of its own quantity (exact structural zeros excluded as identities).") else "")
  cert
}

# D20 (approved October 16, 2026): scale-aware acceptance of the already
# computed conditional error enclosures. Nothing is recomputed, widened or
# narrowed here: the inputs are the same certified boxes the absolute check
# reads, so the private absolute certificate is bound to the public policy
# without weakening it.
#   * scale s_j = lower endpoint of the certified interval of the REPORTED SE
#     (analytic cluster-sum interval, or the multiplier IQR interval under
#     bstrap = TRUE), i.e. the certified positive lower statistical SE scale;
#   * requested estimate, SE, pointwise CI endpoints and band endpoints of
#     target j: certified error <= level x s_j;
#   * covariance entry (j, k): certified error <= level x s_j x s_k (products of
#     the two corresponding lower scales; dimension of a covariance);
#   * multiplier raw-draw enclosures/errors propagate into the requested IQR-SE
#     and band boxes; their own-scale ratios are diagnostics, not a new public
#     intermediate acceptance threshold;
#   * exact structural reference zeros (zero estimate and identically zero IF
#     column) are identities: their errors must be exactly zero and no scale is
#     inverted; every other zero, negative or unresolved scale refuses.
# All ratios use the accepted outward primitives (gmm_iv_ratio_upper, gmm_iv_mul).
# A refusal certifies only that the conservative error bound could not be shown
# to meet the level; it is not a statement about the estimate's value, the
# sampling distribution or the identification of any ATT.
ipw_pipeline_enforce_d20 <- function(cert, estimates, se, vcov, level, err_est, err_se, err_v, band_lower, band_upper) {
  J <- length(estimates)
  ids <- names(estimates); if (is.null(ids)) ids <- paste0("t", seq_len(J))
  deterministic <- if (!is.null(cert$deterministic)) as.integer(cert$deterministic) else
    which(as.double(se) == 0 & as.double(estimates) == 0)
  # The public exception needs BOTH the authoritative reference operator and
  # exact singleton-zero input enclosures. A data-specific zero is insufficient.
  identity <- vapply(seq_len(J), function(j) {
    !is.null(cert$rows) && all(cert$rows$lo[, j] == 0) && all(cert$rows$hi[, j] == 0) &&
      cert$theta$lo[j] == 0 && cert$theta$hi[j] == 0
  }, logical(1L))
  is_det <- seq_len(J) %in% deterministic & seq_len(J) %in% cert$structural_reference & identity
  s_lo <- as.double(cert$se[1L, ])
  # Refusal diagnostics (section 0b): the original fields first, then the target
  # labels and the already computed certificate state of this check. The
  # message and classes are unchanged.
  refuse <- function(msg, diag = list()) cdid_ipw_pipeline_abort(paste0("D20 scale-aware numerical acceptance refused: ", msg),
    c("cdid_ipw_pipeline_d20_error", "cdid_ipw_pipeline_accuracy_error"), c(list(level = level, lower_se_scale = stats::setNames(s_lo, ids)), diag,
      list(target_ids = ids, structural_zero_targets = ids[is_det], n_total = cert$n_total, n_clusters = cert$n_clusters,
        se_law = if (is.null(cert$se_law)) "reported SE law of the certificate (analytic cluster-sum unless bootstrap IQR)" else cert$se_law,
        critical_value = cert$critical_value, pointwise_constant = cert$pointwise_constant,
        refusal_meaning = "a refusal certifies only that a conservative error bound could not be shown to meet the level; it says nothing about estimate values, sampling uncertainty, population covariance rank or ATT identification")))
  err_ci <- if (!is.null(cert$ci) && !is.null(cert$ci_upper)) {
    zc <- cert$pointwise_constant
    rbind(lower = ipw_pipeline_interval_distance(as.double(estimates) - zc * as.double(se), cert$ci$lo, cert$ci$hi),
          upper = ipw_pipeline_interval_distance(as.double(estimates) + zc * as.double(se), cert$ci_upper$lo, cert$ci_upper$hi))
  } else matrix(0.0, 2L, J, dimnames = list(c("lower", "upper"), NULL))
  err_band_lo <- ipw_pipeline_interval_distance(band_lower, cert$band$lo, cert$band$hi)
  err_band_hi <- ipw_pipeline_interval_distance(band_upper, cert$band_upper$lo, cert$band_upper$hi)
  err_mult <- if (!is.null(cert$multiplier_error_by_target)) as.double(cert$multiplier_error_by_target) else rep(0.0, J)
  mult_scale <- if (!is.null(cert$raw_draw_scale_lower)) as.double(cert$raw_draw_scale_lower) else rep(NA_real_, J)
  per_target <- rbind(estimate = as.double(err_est), se = as.double(err_se), ci_lower = err_ci["lower", ], ci_upper = err_ci["upper", ],
                      band_lower = err_band_lo, band_upper = err_band_hi)
  dimnames(per_target) <- list(rownames(per_target), ids)   # labels only; the values are the same certified errors
  if (any(!is.finite(per_target)) || any(!is.finite(err_v)) || any(!is.finite(err_mult))) refuse("a certified error bound is not finite.")
  ratios <- matrix(NA_real_, nrow(per_target), J, dimnames = list(rownames(per_target), ids))
  ratio_mult <- rep(NA_real_, J); names(ratio_mult) <- ids
  contains_zero <- function(lo, hi) all(is.finite(c(lo, hi))) && all(lo <= 0) && all(hi >= 0)
  for (j in seq_len(J)) {
    if (is_det[j]) {
      # Exact structural zero: an identity, never a ratio. The estimate and its
      # IF column are exactly zero by construction, the SE interval is exactly
      # zero, and every certified box must contain zero (outward rounding of an
      # exact zero may leave a box of one subnormal around it; its width is not
      # an error of the identity and no scale is inverted).
      boxes_ok <- all(cert$se[, j] == 0) && estimates[j] == 0 && se[j] == 0 && err_mult[j] == 0 &&
        contains_zero(cert$theta$lo[j], cert$theta$hi[j]) && contains_zero(cert$band$lo[j], cert$band$hi[j]) &&
        contains_zero(cert$band_upper$lo[j], cert$band_upper$hi[j]) &&
        (is.null(cert$ci) || (contains_zero(cert$ci$lo[j], cert$ci$hi[j]) && contains_zero(cert$ci_upper$lo[j], cert$ci_upper$hi[j])))
      if (!boxes_ok)
        refuse(paste0("structural reference zero `", ids[j], "` is not an exact identity in the certificate (nonzero value, SE interval or a box excluding zero)."),
          list(target = ids[j], errors = per_target[, j]))
      ratios[, j] <- 0.0; ratio_mult[j] <- 0.0
      next
    }
    if (!is.finite(s_lo[j]) || s_lo[j] <= 0)
      refuse(paste0("target `", ids[j], "` has no certified positive lower SE scale (lower endpoint ", format(s_lo[j], digits = 6),
        "); the error cannot be normalized and the level is never relaxed."), list(target = ids[j]))
    for (q in seq_len(nrow(per_target))) ratios[q, j] <- gmm_iv_ratio_upper(per_target[q, j], s_lo[j])
    ratio_mult[j] <- if (err_mult[j] == 0) 0.0 else if (is.finite(mult_scale[j]) && mult_scale[j] > 0)
      gmm_iv_ratio_upper(err_mult[j], mult_scale[j]) else NA_real_

  }
  ratio_v <- matrix(NA_real_, J, J, dimnames = list(ids, ids))
  for (j in seq_len(J)) for (k in seq_len(J)) {
    if (is_det[j] || is_det[k]) {
      # Identity: the entry is exactly zero and its certified box contains zero.
      if (vcov[j, k] != 0 || !contains_zero(cert$vcov$lo[j, k], cert$vcov$hi[j, k]))
        refuse(paste0("covariance entry (", ids[j], ", ", ids[k], ") involving a structural zero is not an exact identity."))
      ratio_v[j, k] <- 0.0
    } else {
      den <- gmm_iv_mul(gmm_iv(s_lo[j]), gmm_iv(s_lo[k]))[1L]
      ratio_v[j, k] <- gmm_iv_ratio_upper(err_v[j, k], den)
    }
  }
  worst <- c(estimate = max(ratios["estimate", ]), se = max(ratios["se", ]), ci_lower = max(ratios["ci_lower", ]),
             ci_upper = max(ratios["ci_upper", ]), band_lower = max(ratios["band_lower", ]), band_upper = max(ratios["band_upper", ]),
             covariance = max(ratio_v))
  requested_worst <- worst
  if (any(!is.finite(requested_worst)) || any(requested_worst > level)) {
    refuse(paste0("the largest certified error ratio (", format(max(requested_worst), digits = 6), " for `", names(requested_worst)[which.max(requested_worst)],
      "`) exceeds the level ", format(level, digits = 6), "; the estimator, sample, weights and inference formulas are unchanged and no relaxation is applied."),
      list(worst_ratios = worst, ratios_by_target = ratios, covariance_ratios = ratio_v, multiplier_ratios = ratio_mult,
        # Already computed error bounds and scales behind the ratios (the same
        # objects the accepted certificate retains on success).
        errors_by_target = per_target, covariance_errors = structure(err_v, dimnames = list(ids, ids)),
        multiplier_errors_by_target = stats::setNames(err_mult, ids), raw_draw_scale_lower = stats::setNames(mult_scale, ids),
        worst_target_by_family = apply(ratios, 1L, function(z) if (all(is.na(z))) NA_character_ else ids[which.max(z)]),
        worst_covariance_entry = if (all(is.na(ratio_v))) c(NA_character_, NA_character_) else
          ids[as.integer(arrayInd(which.max(ratio_v), dim(ratio_v)))],
        theta_box = cert$theta, se_box = cert$se, analytic_se_box = cert$analytic_se))
  }
  cert$d20 <- list(decision = "D20", level = level, certified = TRUE, lower_se_scale = stats::setNames(s_lo, ids),
    structural_zero_targets = ids[is_det], worst_ratios = worst, acceptance_families = names(requested_worst),
    intermediate_diagnostics = list(multiplier_ratio_is_acceptance_cutoff = FALSE, multiplier_max_ratio = max(ratio_mult)), ratios_by_target = ratios, covariance_ratios = ratio_v,
    multiplier_ratios = ratio_mult, raw_draw_scale_lower = stats::setNames(mult_scale, ids),
    errors_by_target = per_target, covariance_errors = structure(err_v, dimnames = list(ids, ids)),
    multiplier_errors_by_target = stats::setNames(err_mult, ids),
    scale_rule = paste0("errors of estimates, SEs and interval endpoints are divided by the certified lower endpoint of the reported SE interval; ",
      "covariance errors by the product of the two lower scales; raw-draw errors and own-scale ratios are intermediate diagnostics propagated into requested SE/band boxes, not an additional cutoff; ",
      "exact structural zeros are identities with exactly zero error and no inversion"),
    refusal_meaning = "a refusal certifies only that a conservative error bound could not be shown to meet the level; it says nothing about estimate values, sampling uncertainty or ATT identification")
  cert
}

# ------------------------------------------------------------------------------
# 0b. Refusal diagnostic retention (diagnostics only; October 17, 2026)
# ------------------------------------------------------------------------------
# The three refusal stages demonstrated by the frozen public pilot (nuisance
# root certificate INCONCLUSIVE, D20 scale-aware refusal, Cholesky refusal in
# the GMM core) raised typed conditions that carried almost none of the already
# computed state. The helpers below attach compact packets of that state to the
# SAME condition objects (class, message and the original diagnostic fields are
# unchanged) inside the existing exception paths. No numerical operation,
# control value, acceptance rule, estimator, sample, weight, moment family or
# public option changes; nothing is written to disk, no global state exists,
# and no owner, environment, function or full IF array is retained. Packets
# that would exceed the already accepted retained-diagnostic allowance are
# reduced to the minimal valid packet, with the seam reported in the packet.

# Re-signal the identical condition with one added diagnostic field. `state` is
# evaluated only on a matching refusal; other errors pass through untouched.
ipw_pipeline_resignal_with_state <- function(expr, class, field, state) {
  tryCatch(expr, error = function(e) {
    if (inherits(e, class) && (is.null(e$diagnostics) || is.list(e$diagnostics)) && is.null(e$diagnostics[[field]])) {
      packet <- tryCatch(ipw_pipeline_scalar_state(state(e)), error = function(err)
        list(retention_failed = TRUE, reason = conditionMessage(err)))
      e$diagnostics[[field]] <- packet
    }
    stop(e)
  })
}

# Keep atomic vectors, matrices, data frames and plain lists; drop environments,
# functions, external pointers and sealed references wherever they appear.
ipw_pipeline_scalar_state <- function(x) {
  if (is.environment(x) || is.function(x) || typeof(x) %in% c("externalptr", "closure", "builtin", "special")) return(NULL)
  if (is.data.frame(x)) {
    keep <- vapply(x, function(col) is.atomic(col) || is.list(col), logical(1L))
    x <- x[, keep, drop = FALSE]
    for (nm in names(x)) {
      # Preserve list-column positions, including NULLs, to keep row alignment.
      x[[nm]] <- if (is.list(x[[nm]]) && !is.data.frame(x[[nm]])) lapply(x[[nm]], ipw_pipeline_scalar_state) else ipw_pipeline_scalar_state(x[[nm]])
    }
    attributes(x) <- list(names = names(x), row.names = attr(x, "row.names"), class = "data.frame")
    return(x)
  }
  if (is.list(x)) {
    out <- lapply(x, ipw_pipeline_scalar_state)
    out <- out[!vapply(out, is.null, logical(1L))]
    # Diagnostic lists need names only; arbitrary attributes can hide owners.
    attributes(out) <- list(names = names(out))
    return(out)
  }
  if (is.atomic(x)) {
    at <- attributes(x)
    allowed <- c("names", "dim", "dimnames", if (is.factor(x)) c("levels", "class"))
    keep <- at[names(at) %in% allowed]
    # Ordinary retained vectors/matrices keep their exact payload reference.
    # Only unsupported attributes are stripped; no matrix arithmetic is run.
    if (!identical(at, keep)) attributes(x) <- keep
    return(x)
  }
  NULL
}

# Nuisance-stage packet: the already computed root-certificate scalars, the
# fitted-view diagnostics and risk counts of the ONE refusing nuisance. The
# accepted retained-diagnostic contract reserves 24 entries per nuisance (at most
# one per link); at refusal only this nuisance is materialized, so the full
# packet must fit 24 x n_links entries. Otherwise the minimal valid packet (the
# radius/certificate record that fits 24 entries) is retained and the seam named.
ipw_pipeline_nuisance_refusal_state <- function(prob, fv, witness, link_id, n_certified_before, n_links) {
  cert <- prob$certificate
  num1 <- function(z) if (is.null(z) || !length(z)) NA_real_ else as.double(z)[1L]
  chr1 <- function(z) if (is.null(z) || !length(z)) NA_character_ else as.character(z)[1L]
  lgl1 <- function(z) if (is.null(z) || !length(z)) NA else isTRUE(z)
  radius_record <- list(
    nuisance_id = chr1(fv$nuisance_id), mode = chr1(cert$radius_selection$mode), branch = chr1(cert$radius_selection$branch),
    r = num1(prob$raw$r), delta = num1(cert$delta), t = num1(cert$t), a = num1(cert$a), x = num1(cert$x), q = num1(cert$q),
    selfmap = num1(cert$selfmap), rho = num1(cert$rho), inside_sufficient_region = lgl1(cert$radius_selection$inside_sufficient_region),
    certified = lgl1(prob$certified), inconclusive_type = chr1(prob$inconclusive_type),
    inequalities = list(a_gt_0 = lgl1(cert$inequalities$a_gt_0), delta_lt_1 = lgl1(cert$inequalities$delta_lt_1),
      x_lt_1 = lgl1(cert$inequalities$x_lt_1), q_lt_1 = lgl1(cert$inequalities$q_lt_1), selfmap_le_r = lgl1(cert$inequalities$selfmap_le_r)),
    arithmetic_refusal = list(occurred = lgl1(prob$arithmetic_refusal$occurred), reason = chr1(prob$arithmetic_refusal$reason)))
  n_radius <- 12L + 2L + 5L + 2L
  fit <- fv$diagnostics; counts <- fv$counts; dims <- fv$dimensions
  by_cohort_risk <- cert$dimensions$by_cohort_risk_count
  by_cohort_observed <- cert$dimensions$by_cohort_observed_count
  n_extended <- 11L + 9L + length(by_cohort_risk) + length(by_cohort_observed) + 4L + 3L
  allowance <- 24.0 * as.double(n_links)
  full_fits <- (n_radius + n_extended) <= allowance
  # Admit before constructing extended state; the minimal branch never forces
  # the omitted witness/fit/count payload. The conservative entry count covers
  # numeric/logical leaves including packet and resource headers (not heap/RSS).
  extended <- if (full_fits) list(
    fit_diagnostics = list(converged = lgl1(fit$converged), boundary = lgl1(fit$boundary), iter = num1(fit$iter), rank = num1(fit$rank),
      n_covariates = num1(fit$n_covariates), df_residual = num1(fit$df_residual), deviance = num1(fit$deviance),
      null_deviance = num1(fit$null_deviance), aic = num1(fit$aic),
      glm_control = list(epsilon = num1(fit$control$epsilon), maxit = num1(fit$control$maxit))),
    counts = list(n_a_units = num1(counts$n_a_units), n_risk_a_units = num1(counts$n_risk_a_units),
      n_outside_risk_a_units = num1(counts$n_outside_risk_a_units), n_observed_pair_units = num1(counts$n_observed_pair_units),
      n_observed_risk_units = num1(counts$n_observed_risk_units), n_nonobserved_risk_units = num1(counts$n_nonobserved_risk_units),
      n_treated_obs = num1(counts$n_treated_obs), n_control_obs = num1(counts$n_control_obs),
      m_cohorts = num1(dims$m_cohorts), by_cohort_risk_count = by_cohort_risk, by_cohort_observed_count = by_cohort_observed),
    witness = list(F_dim = dim(as.matrix(witness$F)),
      r_request = if (is.numeric(witness$r)) as.double(witness$r)[1L] else NA_real_,
      r_request_kind = if (is.numeric(witness$r)) "explicit_numeric_radius" else chr1(class(witness$r)[1L])),
    work = list(actual_interval_ops = num1(cert$work$actual_interval_ops),
      actual_reader_work_ops = num1(cert$work$actual_reader_work_ops), endpoints_read = num1(cert$work$endpoints_read))) else NULL
  packet <- c(list(stage = "nuisance_root_certificate", link_id = chr1(link_id), nuisance_index = as.integer(n_certified_before) + 1L,
      n_links = as.integer(n_links), radius_record = radius_record),
    if (full_fits) extended else list(),
    list(resource = list(retained_packet = if (full_fits) "full" else "minimal",
        packet_entries = as.double(n_radius + if (full_fits) n_extended else 0L),
        full_packet_entries = as.double(n_radius + n_extended), allowance_entries = allowance,
        contract = "accepted retained-diagnostic reservation: 24 entries per nuisance, at most one nuisance per link; only the refusing nuisance is materialized at refusal",
        seam = if (full_fits) "none" else paste0("the full packet (", n_radius + n_extended, " entries) exceeds the 24 x n_links = ",
          allowance, " reservation; fit diagnostics, counts, witness and work fields were not retained (no new budget)")),
      meaning = paste0("The sufficient root certificate failed on the stored coefficient at the selected radius; the reasons ",
        "are inequality failures, not evidence against finite-root existence, separation or identification. No radius ",
        "search, retry or fallback was performed.")))
  packet
}

# Theorem 12 attribution packet from the already computed transfer certificate:
# the certified scalars (epsilon_m, epsilon_V, l, l - epsilon_V, nu_I, delta_W),
# the three point and four row-IF error components, the per-link root/evaluation
# split of the moment error, the loading-column target labels and the per-nuisance
# radius record. Existing vectors are referenced; the per-link split adds up to
# four numeric scalars per link (including epsilon_m_vec) within the existing
# 24-per-link reservation. List/table metadata is not a heap/RSS guarantee.
ipw_pipeline_transfer_error_attribution <- function(num_transfer, gmm_fit = NULL) {
  if (is.null(num_transfer)) return(list(available = FALSE, reason = "numerical transfer was not requested"))
  gt <- num_transfer$gmm_transfer
  split <- NULL
  if (!is.null(num_transfer$link_enclosures)) {
    le <- num_transfer$link_enclosures
    pick <- function(e, nm) { z <- e$separated_errors$mu[[nm]]; if (is.null(z)) NA_real_ else as.double(z)[1L] }
    split <- data.frame(link_id = vapply(le, function(e) if (is.null(e$link_id)) NA_character_ else as.character(e$link_id), character(1L)),
      eval_error = vapply(le, pick, numeric(1L), nm = "eval_error"), root_error = vapply(le, pick, numeric(1L), nm = "root_error"),
      total_error_bound = vapply(le, pick, numeric(1L), nm = "total_error_bound"), stringsAsFactors = FALSE)
    if (length(gt$epsilon_m_vec) == nrow(split)) split$epsilon_m_vec <- as.double(gt$epsilon_m_vec)
  }
  loading_ids <- if (!is.null(gmm_fit) && is.data.frame(gmm_fit$targets)) gmm_fit$targets$target_id[gmm_fit$targets$estimable] else NULL
  list(available = TRUE, theorem_ref = gt$theorem_ref, weighting = num_transfer$weighting,
    epsilon_m = gt$epsilon_m, epsilon_V = gt$epsilon_V, nu_I = gt$nu_I, l_lower = gt$l_lower, u_upper = gt$u_upper,
    l_minus_eps_V = gt$l_minus_eps_V, delta_W = gt$delta_W,
    certified_lower_eigenvalue = gt$certified_lower_eigenvalue,
    loading_target_ids = loading_ids,
    point_error_components = gt$point_error_components, R_c_components = gt$R_c_components,
    total_target_error_cert = gt$total_target_error_cert, total_loading_error_cert = gt$total_loading_error_cert,
    R_c_cert = gt$R_c_cert, total_se_error_cert = gt$total_se_error_cert, total_cov_error_cert = gt$total_cov_error_cert,
    total_cov_diag_error_cert = gt$total_cov_diag_error_cert, projection_error_cert = gt$projection_error_cert,
    accumulation_se_error_cert = gt$accumulation_se_error_cert, accumulation_cov_error_cert = gt$accumulation_cov_error_cert,
    s_star_lower = gt$s_star_lower, s_star_upper = gt$s_star_upper,
    moment_error_split = split, nuisance_root_radii = num_transfer$nuisance_root_radii,
    scope = paste0("Upper endpoints of already certified Theorem 12 (thm:gmm) terms: attribution of conservative bounds, ",
      "not measured error. A D20 refusal downstream says nothing about estimate values, sampling uncertainty, ",
      "population covariance rank or ATT identification."))
}

# ------------------------------------------------------------------------------
# 1. Shared Nuisance & Point-Link Construction Over Qualified Catalog
# ------------------------------------------------------------------------------
# Iterates over retained comparison catalog rows `1..M` in `metadata_obj`:
#   - Maintains a registry of distinct `cdid_ipw_nuisance_view` objects indexed
#     by `operator_id` and actual pair-observation mask `O_l = R_l intersect S_l`.
#   - When two catalog rows share `operator_id`, tries the operator's views in
#     creation order. A differing outcome-free fingerprint (`fingerprint_ipw_view_mask`)
#     proves the masks differ, so that split is recorded without building the
#     link. Equal fingerprints prove nothing: the link is then constructed with
#     `verify_nuisance_mask = TRUE`, whose chunk-by-chunk zero-outcome
#     `verify_ipw_view_mask_compatibility` decides and resolves any collision.
#   - If the actual pair-observation masks match across all A-units, both links
#     share the exact same `nuisance_view` and the exact same `fitted_view`.
#   - If the actual masks differ (`cdid_ipw_shared_mask_error`), records the mask
#     split and constructs a separate `nuisance_view` for the new mask branch.
#   - Trust boundary: the catalog, frame owner and controls are validated at
#     entry with the arguments every per-row call previously repeated; the loop
#     runs no caller code; the owner constructors still validate the frame,
#     control and nuisance view on every call; each row keeps its preflights.
#     Fingerprints live only in this call and are never cached across calls.
#   - CP08 scaling action 2: every streamed scan (fingerprint pass, row
#     construction, exact collision resolution) is charged before it runs
#     against the caller's composed work, so adversarial fingerprint collisions
#     are counted as they occur. Split records and the sealed constructor-owned
#     mask-equivalence certificate are reserved at their worst case up front.
#     The certificate binds every row to its frame, operator, cohorts, ordered
#     controls, link view and nuisance/fit identity; within one operator its
#     distinct views have pairwise unequal masks (a fingerprint difference or an
#     exact streamed mismatch separated each pair). Structural certification
#     consumes it instead of re-streaming those masks.

build_pipeline_nuisance_and_link_views <- function(
  metadata_obj,
  frame_owner = NULL,
  owner_control,
  fitted_control,
  point_control,
  glm_control,
  retained_pipeline_entries = 0,
  fit_point_links = TRUE,
  external_work_ops = 0,
  progress = NULL
) {
  validate_owner_comparison_catalog(metadata_obj, frame_owner = frame_owner)
  own_ctrl <- validate_ipw_owner_control(owner_control)
  fit_ctrl <- validate_ipw_fitted_view_control(fitted_control)
  pt_ctrl <- validate_ipw_point_view_control(point_control)
  glm_ctrl <- validate_ipw_glm_control(glm_control)
  if (!gmm_scalar_logical(fit_point_links)) {
    cdid_ipw_pipeline_abort("`fit_point_links` must be a scalar logical.", "cdid_ipw_pipeline_input_error")
  }

  cat_rows <- metadata_obj$catalog$rows
  M <- nrow(cat_rows)
  if (M < 1L) {
    cdid_ipw_pipeline_abort(
      "Comparison catalog contains zero qualified links; cannot construct IPW views.",
      c("cdid_ipw_pipeline_empty_catalog_error", "cdid_gmm_target_estimability_error")
    )
  }

  f_owner <- if (!is.null(frame_owner)) {
    frame_owner
  } else if (!is.null(metadata_obj$frame_ref)) {
    get("frame_owner", envir = metadata_obj$frame_ref, inherits = FALSE)
  } else {
    cdid_ipw_pipeline_abort(
      "`frame_owner` must be supplied either on `metadata_obj` or explicitly.",
      "cdid_ipw_pipeline_input_error"
    )
  }
  validate_ipw_frame_owner(f_owner)
  # Checks each `construct_pipeline_catalog_row_views` call repeated with these
  # same immutable arguments, performed once for the whole trusted loop.
  validate_owner_comparison_catalog(metadata_obj, frame_owner = f_owner)
  validate_ipw_frame_owner(f_owner, panel_owner = get("panel_owner", envir = metadata_obj$panel_ref, inherits = FALSE))
  # Invocation-local fingerprint registry, charged before anything is stored:
  # per distinct view a list slot and four named doubles; for the current row
  # the same plus its pair descriptor and cohort index (at most every cohort).
  # Additional 128 slots per row cover its support-record wrapper/scalars,
  # fingerprint record and bounded scalar key, and at most two date facts/keys.
  # Cohort/support vectors are referenced, never copied. The previous 16 slots
  # per row still cover the original screening registry independently. The
  # extra reservation is composed with the constructor in split_cert_entries
  # below, preserving the established pre-construction refusal ordering.
  mask_registry_entries <- ipw_owner_safe_sum(c(ipw_owner_safe_product(16, M + 1), length(metadata_obj$catalog$spec$cohorts)),
    own_ctrl$max_retained_entries, "invocation-local mask fingerprint registry")
  # The returned split records (at most n_op (n_op - 1) / 2 per operator, eight
  # slots each, adversarial collisions included) and the certificate (three
  # identity slots per row; cohorts, ordered controls and two references per
  # view, at most one view per row) are retained with the registry; they join
  # every composed registry reservation below, before the first row is built.
  C_all <- length(metadata_obj$catalog$spec$cohorts)
  op_rows <- as.double(tabulate(match(cat_rows$operator_id, unique(cat_rows$operator_id))))
  split_bound <- sum(op_rows * (op_rows - 1) / 2)
  split_cert_entries <- ipw_owner_safe_sum(c(ipw_owner_safe_product(8, split_bound, label = "mask split records"),
    ipw_owner_safe_product(M, 11 + C_all, label = "mask-equivalence certificate"),
    ipw_owner_safe_product(128, M + 1, label = "construction support, fingerprint and date fact records")))
  retained_pipeline_entries <- ipw_owner_int_bound(retained_pipeline_entries, "retained_pipeline_entries", min_val = 0)
  external_work_ops <- ipw_owner_int_bound(external_work_ops, "external_work_ops", min_val = 0)
  # One streamed mask pass (fingerprint, support scan or exact comparison) is
  # charged with the accepted per-stream envelope before it starts.
  stream_n <- as.double(f_owner$n_a_units)
  stream_step <- max(1, min(stream_n, own_ctrl$max_chunk_rows))
  stream_cost <- function(m) 128 * stream_n * (1 + ceiling(log2(m + 1))) +
    128 * ceiling(stream_n / stream_step) * (1 + ceiling(log2(stream_n + 1)))
  scan_work <- 0
  charge_scans <- function(n_streams, m) {
    scan_work <<- ipw_owner_safe_sum(c(scan_work, ipw_owner_safe_product(n_streams, stream_cost(m))))
    ipw_owner_safe_sum(c(external_work_ops, scan_work), own_ctrl$max_work_ops,
      "pipeline work plus mask registry scans (including exact collision resolution)")
  }
  retained_with_masks <- ipw_owner_safe_sum(c(retained_pipeline_entries, mask_registry_entries),
    min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries),
    "pipeline envelope plus invocation-local mask fingerprint registry")
  # The per-row checks below compose the split/certificate reservation with
  # max(envelope, constructor), a bound on envelope + registry + records; the
  # fingerprint workspace is composed with the same live total.
  retained_with_masks <- ipw_owner_safe_sum(c(retained_with_masks, split_cert_entries))
  mismatch_reason <- ipw_owner_mask_mismatch_message()

  nuisance_buckets <- list()
  distinct_nuisance_views <- list()
  distinct_nuisance_keys <- character()
  distinct_operator_ids <- integer()
  distinct_fingerprints <- vector("list", M)
  link_views <- vector("list", M)
  row_specs <- vector("list", M)
  link_to_nuisance_idx <- integer(M)
  mask_splits <- list()
  n_compared <- n_fp_rejected <- n_exact <- n_exact_mismatch <- n_fp_passes <- 0
  # CP08 balanced performance: exact complete-history facts for this trusted
  # loop only (see ipw_owner_construction_context). A link reuses the support
  # scan of a nuisance view built here when their masks are equal by a posting
  # fact; a complete-pair fingerprint is streamed at most once per cohort index
  # while its bounded bucket entry is retained. Every
  # row keeps its constructors, validations, preflights and charges above; the
  # facts live in this call only and are dropped before fitting. Per distinct
  # view they hold one support record and at most one fingerprint record, which
  # reference already-built vectors, plus one logical per checked calendar date
  # (each row checks at most its two dates), in the explicit additional
  # reservation above. Completeness reads use chunks of the row's own step.
  # The context also records, after one more completed validation, that this
  # exact sealed frame owner was checked: the constructors of every row reuse
  # that check and the checked nuisance/support template of their operator
  # (recorded at the template's construction) instead of validating the same
  # sealed snapshots once per row. One reference slot, within the reservation.
  construction_ctx <- ipw_owner_construction_context(f_owner)
  ipw_owner_context_check_frame(construction_ctx, f_owner)
  # CP08 construction: the same descriptor representation as the CP04 registry
  # (pipeline_metadata_descriptor_source): the validated catalog's column
  # vectors read at each row, each referenced pool's K expanded once and
  # shared by its rows' retained descriptors. Dropped with the context before
  # fitting. The pure dimension preflights below are reused only for a row
  # whose (n_A, p, cohort count, control) tuple is identical to the previous
  # row's (one memo record, within the context's reserved block); every row
  # keeps its own two row preflights, its live-envelope admission and its
  # support, mask and owner checks in the original order.
  descriptor_src <- descriptor_charge <- dim_memo <- NULL

  for (m in seq_len(M)) {
    pipeline_metadata_row_spec_preflight(metadata_obj, m, own_ctrl)
    pipeline_metadata_row_spec_preflight(metadata_obj, m, metadata_obj$control$owner_control)
    if (is.null(descriptor_src)) {
      descriptor_src <- pipeline_metadata_descriptor_source(metadata_obj, own_ctrl,
        retained_with_masks, external_work_ops, M, memo_entries = 64)
      descriptor_charge <- descriptor_src$accounting
      split_cert_entries <- ipw_owner_safe_sum(c(split_cert_entries, descriptor_charge$retained_entries))
      retained_with_masks <- ipw_owner_safe_sum(c(retained_with_masks, descriptor_charge$retained_entries))
      external_work_ops <- ipw_owner_safe_sum(c(external_work_ops, descriptor_charge$work_ops_bound),
        own_ctrl$max_work_ops, "IPW registry plus descriptor preparation work")
    }
    # Existing constructor preflights include retained view descriptors and
    # their scratch. The caller's whole-pipeline envelope already covers those
    # buffers; standalone callers instead use the larger row constructor peak.
    # Add the new registry once, before extracting or constructing the row.
    row_cohorts <- 1 + metadata_obj$catalog$pools$size[cat_rows$pool_id[m]]
    dim_inputs <- list(f_owner$n_a_units, f_owner$n_covariates, row_cohorts, own_ctrl)
    if (is.null(dim_memo) || !identical(dim_memo$inputs, dim_inputs)) {
      row_nuisance_preflight <- preflight_ipw_nuisance_view(f_owner$n_a_units, f_owner$n_covariates, row_cohorts, own_ctrl)
      row_link_preflight <- preflight_ipw_link_view(f_owner$n_a_units, f_owner$n_covariates, row_cohorts, own_ctrl)
      dim_memo <- list(inputs = dim_inputs,
        peak = max(row_nuisance_preflight$concurrent_entries, row_link_preflight$concurrent_entries),
        retained = max(row_nuisance_preflight$total_retained_entries, row_link_preflight$total_retained_entries))
      row_nuisance_preflight <- row_link_preflight <- NULL
    }
    row_constructor_peak <- dim_memo$peak
    row_constructor_retained <- dim_memo$retained
    ipw_owner_safe_sum(c(mask_registry_entries, split_cert_entries, max(retained_pipeline_entries, row_constructor_retained)),
      own_ctrl$max_retained_entries, "mask fingerprint registry plus existing row constructor retention")
    ipw_owner_safe_sum(c(mask_registry_entries, split_cert_entries, max(retained_pipeline_entries, row_constructor_peak),
      descriptor_charge$scratch_entries), own_ctrl$max_matrix_entries,
      "mask fingerprint registry plus row constructor and descriptor scratch envelope")
    spec <- pipeline_metadata_row_spec_fields(metadata_obj, m, descriptor_src)
    op_id <- as.integer(cat_rows$operator_id[m])
    op_key <- paste0("op_", op_id)
    bucket <- nuisance_buckets[[op_key]]
    matched_idx <- NA_integer_
    row_res <- NULL
    row_fp <- NULL

    if (!is.null(bucket) && length(bucket) > 0L) {
      for (cand_idx in bucket) {
        cand_nv <- distinct_nuisance_views[[cand_idx]]
        n_compared <- n_compared + 1
        # Fingerprints decide only where the exact attempt could end in nothing
        # but its streamed mask comparison: same frame, identical g and ordered K.
        if (identical(cand_nv$frame_ref, f_owner$ref) &&
            identical(as.double(spec$g), cand_nv$cohorts$g) &&
            identical(as.double(spec$K), cand_nv$cohorts$K)) {
          if (is.null(row_fp) || !identical(row_fp$cohort_index, cand_nv$cohorts$owner_cohort_index)) {
            charge_scans(1, cand_nv$cohorts$m_cohorts)
            fp_row <- fingerprint_ipw_view_mask(cand_nv, own_ctrl, s = spec$s, t = spec$t,
              retained_entries = retained_with_masks, construction_context = construction_ctx)
            row_fp <- list(cohort_index = cand_nv$cohorts$owner_cohort_index, pair = fp_row$pair,
              fingerprint = fp_row$fingerprint)
            n_fp_passes <- n_fp_passes + 1
          }
          if (is.null(distinct_fingerprints[[cand_idx]])) {
            charge_scans(1, cand_nv$cohorts$m_cohorts)
            distinct_fingerprints[[cand_idx]] <- fingerprint_ipw_view_mask(cand_nv, own_ctrl,
              retained_entries = retained_with_masks, construction_context = construction_ctx)$fingerprint
            n_fp_passes <- n_fp_passes + 1
          }
          if (!identical(row_fp$fingerprint, distinct_fingerprints[[cand_idx]])) {
            n_fp_rejected <- n_fp_rejected + 1
            mask_splits[[length(mask_splits) + 1L]] <- list(
              row_id = m,
              operator_id = op_id,
              candidate_nuisance_id = cand_nv$nuisance_id,
              reason = mismatch_reason
            )
            next
          }
        }
        n_exact <- n_exact + 1
        # Link support scan plus the exact streamed comparison; equal
        # fingerprints (collisions included) never skip this charge.
        charge_scans(2, cand_nv$cohorts$m_cohorts)
        attempt <- tryCatch(
          pipeline_metadata_row_views_from_spec(spec, f_owner, own_ctrl,
            nuisance_view = cand_nv, verify_nuisance_mask = TRUE, construction_context = construction_ctx),
          cdid_ipw_shared_mask_error = function(e) e
        )
        if (!inherits(attempt, "error")) {
          matched_idx <- cand_idx
          row_res <- attempt
          break
        } else {
          n_exact_mismatch <- n_exact_mismatch + 1
          mask_splits[[length(mask_splits) + 1L]] <- list(
            row_id = m,
            operator_id = op_id,
            candidate_nuisance_id = cand_nv$nuisance_id,
            reason = conditionMessage(attempt)
          )
        }
      }
    }

    if (is.na(matched_idx)) {
      # Nuisance and link support scans plus the link's streamed verification.
      charge_scans(3, row_cohorts)
      row_res <- pipeline_metadata_row_views_from_spec(spec, f_owner, own_ctrl,
        nuisance_view = NULL, verify_nuisance_mask = TRUE, construction_context = construction_ctx)
      new_idx <- length(distinct_nuisance_views) + 1L
      distinct_nuisance_views[[new_idx]] <- row_res$nuisance_view
      distinct_operator_ids[new_idx] <- op_id
      branch_num <- if (is.null(bucket)) 1L else (length(bucket) + 1L)
      distinct_nuisance_keys[new_idx] <- paste0("op", op_id, ":mask", branch_num)
      nuisance_buckets[[op_key]] <- c(bucket, new_idx)
      matched_idx <- new_idx
      # The row fingerprint streamed this exact pair under these cohorts on this
      # frame; it is the new view's own fingerprint only when that literally holds.
      nv_new <- row_res$nuisance_view
      if (!is.null(row_fp) && identical(nv_new$frame_ref, f_owner$ref) &&
          identical(nv_new$cohorts$owner_cohort_index, row_fp$cohort_index) &&
          identical(nv_new$defining_pair, row_fp$pair)) {
        distinct_fingerprints[[new_idx]] <- row_fp$fingerprint
      }
    }

    row_specs[[m]] <- row_res$row_spec
    link_views[[m]] <- row_res$link_view
    link_to_nuisance_idx[m] <- matched_idx
    if (is.function(progress)) progress(m, M, length(distinct_nuisance_views))
  }

  # The screening cache has no role in fitting and must not persist into its
  # existing workspace envelope or the later GMM/IF phases.
  distinct_fingerprints <- row_fp <- fp_row <- NULL
  descriptor_src <- dim_memo <- dim_inputs <- NULL
  n_support_reused <- construction_ctx$support_reuses
  n_fp_reused <- construction_ctx$fingerprint_reuses
  construction_ctx <- NULL

  # Fit each distinct nuisance view once with bounded observer lifecycle
  fitted_views <- fit_ipw_nuisance_view_set(
    nuisance_views = distinct_nuisance_views,
    glm_control = glm_ctrl,
    control = fit_ctrl
  )
  names(fitted_views) <- distinct_nuisance_keys
  names(distinct_nuisance_views) <- distinct_nuisance_keys

  # Constructor-owned, sealed mask-equivalence certificate (reserved above).
  mask_certificate <- ipw_seal(list(
    kind = "cdid_ipw_registry_mask_equivalence",
    frame_ref = f_owner$ref, panel_ref = f_owner$panel_ref, metadata_ref = metadata_obj$ref,
    n_rows = as.double(M), row_operator = as.integer(cat_rows$operator_id),
    row_view = as.integer(link_to_nuisance_idx),
    row_link_refs = unname(lapply(link_views, `[[`, "ref")),
    view_operator = as.integer(distinct_operator_ids),
    view_g = unname(vapply(distinct_nuisance_views, function(v) as.double(v$cohorts$g), numeric(1L))),
    view_K = unname(lapply(distinct_nuisance_views, function(v) as.double(v$cohorts$K))),
    view_nuisance_refs = unname(lapply(distinct_nuisance_views, `[[`, "ref")),
    view_fit_refs = unname(lapply(fitted_views, `[[`, "ref")),
    separations = list(fingerprint_differences = n_fp_rejected, exact_streamed_mismatches = n_exact_mismatch),
    claim = paste0("Within one operator (target cohort and ordered control pool under the declared conventions) ",
      "the distinct views have pairwise unequal risk-set pair-observation masks: each pair was separated by a ",
      "fingerprint difference (a deterministic function of the exact mask) or by an exact streamed mismatch; ",
      "equal fingerprints were never accepted without the exact comparison (streamed, or decided by an exact ",
      "posting fact: identical pair or four complete dates). Each row's sealed link view carries a completed ",
      "verification against its own view."),
    scope = "invocation-local evidence of this registry only; no statistical equivalence, sample or estimate is inferred"),
    "mask_certificate")

  # Fit compact point link views reusing the shared fitted nuisance views. The
  # outcome-free reduction route (CP08 P3) fits only its certified retained rows.
  point_links <- NULL
  if (isTRUE(fit_point_links)) {
    point_links <- fit_pipeline_point_links(fitted_views, link_views, link_to_nuisance_idx, seq_len(M), pt_ctrl)
  }
  names(link_views) <- vapply(link_views, `[[`, character(1L), "link_id")

  list(
    nuisance_views = distinct_nuisance_views,
    fitted_views = fitted_views,
    nuisance_keys = distinct_nuisance_keys,
    operator_ids = distinct_operator_ids,
    row_specs = row_specs,
    link_views = link_views,
    point_links = point_links,
    link_to_nuisance_idx = link_to_nuisance_idx,
    link_nuisance_keys = distinct_nuisance_keys[link_to_nuisance_idx],
    mask_splits = mask_splits,
    mask_certificate = mask_certificate,
    construction_work = list(
      n_rows = M,
      n_distinct_nuisance_views = length(distinct_nuisance_views),
      candidate_comparisons = n_compared,
      fingerprint_rejections = n_fp_rejected,
      exact_candidate_attempts = n_exact,
      exact_mask_mismatches = n_exact_mismatch,
      fingerprint_passes = n_fp_passes,
      link_view_constructions = n_exact + length(distinct_nuisance_views),
      endpoints_read = 0L,
      endpoints_read_scope = "mask_screening_only",
      work_ops_charged = ipw_owner_safe_sum(c(scan_work, descriptor_charge$work_ops_bound)),
      registry_entries_reserved = mask_registry_entries + split_cert_entries,
      split_records_bound = split_bound,
      link_support_scans_reused = n_support_reused,
      complete_pair_fingerprints_reused = n_fp_reused,
      scope = paste0("Invocation-local counts. Every mask_splits entry is a fingerprint rejection (masks proved ",
        "different, no link built) or an exact streamed mismatch after equal fingerprints. ",
        "endpoints_read covers mask screening only and excludes nuisance and point fitting. ",
        "work_ops_charged includes descriptor preparation and the per-stream envelope charged before each streamed scan; the stream charge stays charged when ",
        "an exact posting fact (identical pair or complete dates) decides the scan without streaming it: ",
        "link_support_scans_reused links took their nuisance view's support scan, ",
        "complete_pair_fingerprints_reused fingerprints were the memoized complete-roster value, and mask ",
        "verifications of such pairs return without streaming.")
    )
  )
}

# Compact point links for the given catalog rows, in increasing row order, each
# reusing its row's single shared fitted nuisance view (no refit; reads only that
# link's outcome endpoints).
fit_pipeline_point_links <- function(fitted_views, link_views, link_to_nuisance_idx, rows, point_control) {
  out <- vector("list", length(rows))
  for (h in seq_along(rows)) {
    m <- rows[h]
    out[[h]] <- fit_ipw_point_link_view(
      fitted_view = fitted_views[[link_to_nuisance_idx[m]]],
      link_view = link_views[[m]],
      control = point_control
    )
  }
  names(out) <- vapply(out, `[[`, character(1L), "link_id")
  out
}

# CP08 stage 2: facts that select the trusted batch execution route (see
# Section 2b of R/chained_ipw_point_view.R) on the outcome-free reduction
# route, before any retained point fit. All retained rows share one fitted
# nuisance whose single target cohort has the never-treated pool as its only
# control cohort; the frame has one shared fixed-unit propensity design (the
# owner's n_A x p double matrix, p >= 1 columns including the intercept, already
# validated finite by the frame owner; `~x` has p = 2) and constant unit
# weights; every retained link view is bound to that nuisance's owners and
# snapshot with its cohort descriptor and one functional key; and the defining
# pair and every retained pair have complete dates (each a sealed posting slice
# exactly 1..n_A). NULL keeps the general route. Nothing here raises, reads an
# outcome, scans the design or fits; it reads only the existing chunked
# per-date posting facts, memoized for the batch fitter. Admission is
# structural only: never by coefficient, loading or output closeness.
ipw_pipeline_batch_eligibility <- function(registry, kept_idx, frame_owner, point_control) {
  if (!is.list(registry) || length(kept_idx) < 1L || !is.list(frame_owner)) return(NULL)
  nidx <- unique(registry$link_to_nuisance_idx[kept_idx])
  if (length(nidx) != 1L || is.na(nidx) || nidx < 1L || nidx > length(registry$fitted_views)) return(NULL)
  fv <- registry$fitted_views[[nidx]]
  nv <- registry$nuisance_views[[nidx]]
  if (!inherits(fv, "cdid_ipw_fitted_nuisance_view") || !inherits(nv, "cdid_ipw_nuisance_view") ||
      !is.environment(fv$nuisance_view_ref) || !identical(nv$ref, fv$nuisance_view_ref) ||
      !identical(fv$frame_ref, frame_owner$ref) || !identical(nv$frame_ref, frame_owner$ref)) return(NULL)
  coh <- nv$cohorts
  if (length(coh$g) != 1L || !identical(coh$K, Inf) || !identical(as.integer(coh$m_cohorts), 2L) ||
      !identical(fv$cohorts$cohorts, coh$cohorts)) return(NULL)
  Z <- frame_owner$Z_a
  w <- frame_owner$w_a
  p <- frame_owner$n_covariates
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 1 || p != floor(p) ||
      !is.matrix(Z) || !is.double(Z) ||
      !identical(dim(Z), c(as.integer(frame_owner$n_a_units), as.integer(p))) ||
      !is.double(w) || length(w) < 1L || !isTRUE(all(w == w[1L]))) return(NULL)
  po <- if (is.environment(frame_owner$panel_ref)) get0("panel_owner", envir = frame_owner$panel_ref, inherits = FALSE) else NULL
  if (!inherits(po, "cdid_ipw_panel_owner") || !identical(po$n_a_units, frame_owner$n_a_units)) return(NULL)
  step <- min(as.double(frame_owner$n_a_units), as.double(point_control$max_chunk_rows))
  memo <- list2env(list(panel_owner = po), parent = emptyenv())
  if (!isTRUE(ipw_owner_pair_complete(po, nv$defining_pair, step, memo))) return(NULL)
  key <- NULL
  for (m in kept_idx) {
    lv <- registry$link_views[[m]]
    if (!inherits(lv, "cdid_ipw_link_view") || !identical(lv$frame_ref, fv$frame_ref) ||
        !identical(lv$panel_ref, fv$panel_ref) || !identical(lv$nuisance_ref, fv$nuisance_view_ref) ||
        !isTRUE(lv$nuisance_mask_verified) || !identical(lv$cohorts, coh) ||
        !gmm_scalar_character(lv$functional_key) || (!is.null(key) && !identical(lv$functional_key, key)) ||
        !isTRUE(ipw_owner_pair_complete(po, lv$pair, step, memo))) return(NULL)
    key <- lv$functional_key
  }
  list(fitted_view = fv, nuisance_view = nv, date_memo = memo)
}

# ------------------------------------------------------------------------------
# 2. Compact Structural Class Certification & Redundancy Reduction
# ------------------------------------------------------------------------------
# Adapts `certify_ipw_structural_classes` in `R/chained_ipw_shared.R` to compact
# sealed `cdid_ipw_point_link_view` objects and returns a sealed
# `cdid_ipw_structural_certificate` directly compatible with `apply_ipw_reduction_map`:
#   - Refuses structural reduction when `weighting != "optimal"` with
#     `"cdid_ipw_reduction_weighting_error"`.
#   - Validates `design` and `design_link_index` via `validate_ipw_design_for_classes`
#     and `validate_ipw_design_incidence`.
#   - Verifies that any links sharing `(g, K, glm_control)` and an identical actual
#     pair-observation mask (checked chunk-by-chunk with zero outcome reads via
#     `verify_ipw_view_mask_compatibility`) also share the exact same sealed
#     `fitted_nuisance_ref` (refusing separate refits or mismatched coefficients
#     with `"cdid_ipw_shared_fit_mismatch_error"`).
#   - Within each `(nuisance_view_ref, functional_key)` structural class, selects
#     a spanning forest via `incidence_forest`, `ipw_forest_depth`, and
#     `ipw_forest_path`, and builds the sparse expansion map `E = list(i, j, x)`
#     so `apply_ipw_reduction_map(certificate, x_retained)` expands reduced moments
#     or IF rows back to all original links.
#   - Checks cross-class cycle-freeness among retained links sharing identical
#     `(g, K, glm_control, mask)` across operators (`"cdid_ipw_reduction_singular_error"`).

certify_ipw_structural_classes_view <- function(
  point_links,
  design,
  design_link_index = seq_along(point_links),
  weighting = "optimal",
  bounds,
  owner_control
) {
  if (!is.character(weighting) || length(weighting) != 1L || is.na(weighting) ||
      !weighting %in% c("optimal", "identity")) {
    cdid_abort("`weighting` must be exactly \"optimal\" or \"identity\".", "cdid_ipw_input_error")
  }
  if (!identical(weighting, "optimal")) {
    cdid_abort(
      "Structural reduction is licensed only for optimal weighting; identity weighting keeps every row.",
      c("cdid_ipw_reduction_weighting_error", "cdid_ipw_certificate_error")
    )
  }
  b <- validate_ipw_certificate_bounds(bounds)
  own_ctrl <- validate_ipw_owner_control(owner_control)

  if (!is.list(point_links) || is.object(point_links) || length(point_links) < 1L) {
    cdid_abort("`point_links` must be a nonempty plain list.", c("cdid_ipw_input_error", "cdid_ipw_certificate_error"))
  }
  n_links <- length(point_links)
  custom_links <- vapply(point_links, inherits, logical(1L), what = "cdid_custom_ipw_link")
  if (any(custom_links) && !all(custom_links)) {
    cdid_abort("Structural classes cannot mix custom and built-in links.", "cdid_ipw_certificate_error")
  }
  is_custom <- all(custom_links)
  if (n_links > b$max_links) {
    cdid_abort("Number of links exceeds `max_links`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }

  for (ell in seq_len(n_links)) {
    if (is_custom) validate_custom_ipw_link(point_links[[ell]]) else validate_ipw_point_link_view(point_links[[ell]])
    if (!isTRUE(point_links[[ell]]$nuisance_mask_verified)) {
      cdid_abort(
        paste0("Point link ", ell, " does not have `nuisance_mask_verified == TRUE`."),
        c("cdid_ipw_shared_mask_error", "cdid_ipw_certificate_error")
      )
    }
    if (ell > 1L) {
      if (!identical(point_links[[ell]]$frame_ref, point_links[[1L]]$frame_ref) ||
          !identical(point_links[[ell]]$panel_ref, point_links[[1L]]$panel_ref)) {
        cdid_abort(
          "Point links belong to different `frame_owner` or `panel_owner` snapshots.",
          c("cdid_ipw_context_mismatch_error", "cdid_ipw_certificate_error")
        )
      }
    }
  }

  ipw_view_structural_classes_core(point_links, design, design_link_index, b, own_ctrl, is_custom,
    if (is_custom) "custom_original_ipw" else "shared_view_fit")
}

# Binds a registry mask-equivalence certificate (build_pipeline_nuisance_and_
# link_views) to the links being certified: sealed snapshot, same frame, one
# row per link, and per link the same sealed link view, nuisance view, fitted
# nuisance, target cohort and ordered controls as recorded. Any mismatched
# reference refuses; nothing is re-derived from labels. Returns each link's view.
ipw_mask_certificate_bind <- function(cert, links, frame_owner) {
  bad <- function(msg) cdid_abort(msg, c("cdid_ipw_shared_mask_error", "cdid_ipw_certificate_error"))
  if (!is.list(cert) || !identical(cert$kind, "cdid_ipw_registry_mask_equivalence") ||
      !isTRUE(ipw_seal_ok(cert, "mask_certificate")))
    bad("Mask-equivalence certificate is malformed or differs from its sealed snapshot.")
  n <- length(links)
  rv <- cert$row_view; nv <- length(cert$view_nuisance_refs)
  if (!identical(cert$frame_ref, frame_owner$ref) || !identical(cert$n_rows, as.double(n)) ||
      !is.integer(rv) || length(rv) != n || anyNA(rv) || any(rv < 1L | rv > nv) ||
      length(cert$row_link_refs) != n || !is.integer(cert$row_operator) || length(cert$row_operator) != n ||
      length(cert$view_fit_refs) != nv || length(cert$view_K) != nv || length(cert$view_g) != nv ||
      !is.integer(cert$view_operator) || length(cert$view_operator) != nv ||
      any(cert$row_operator != cert$view_operator[rv]))
    bad("Mask-equivalence certificate does not describe these links (mismatched frame, row count or row/view map).")
  for (ell in seq_len(n)) {
    v <- rv[ell]; pl <- links[[ell]]
    if (!identical(pl$link_view_ref, cert$row_link_refs[[ell]]) ||
        !identical(pl$nuisance_view_ref, cert$view_nuisance_refs[[v]]) ||
        !identical(pl$fitted_nuisance_ref, cert$view_fit_refs[[v]]) ||
        !identical(as.double(pl$cohorts$g), cert$view_g[v]) ||
        !identical(as.double(pl$cohorts$K), cert$view_K[[v]])) {
      bad(paste0("Mask-equivalence certificate does not describe link ", ell,
        " (mismatched link view, nuisance/fit identity, cohort or ordered controls)."))
    }
  }
  rv
}

# Streamed comparisons a bound certificate leaves open: per link, the views
# with identical (g, ordered K) under a different operator (none within one
# catalog). NA when the certificate is absent or malformed (the core then
# refuses it before any stream, or streams every comparison without it).
ipw_mask_certificate_open_streams <- function(cert) {
  ok <- is.list(cert) && identical(cert$kind, "cdid_ipw_registry_mask_equivalence") &&
    isTRUE(tryCatch(ipw_seal_ok(cert, "mask_certificate"), error = function(e) FALSE)) &&
    is.integer(cert$row_view) && is.integer(cert$view_operator) && is.list(cert$view_K) &&
    length(cert$view_K) == length(cert$view_operator) && length(cert$view_g) == length(cert$view_operator) &&
    !anyNA(cert$row_view) && all(cert$row_view >= 1L & cert$row_view <= length(cert$view_operator))
  if (!ok) return(NA_real_)
  key <- vapply(seq_along(cert$view_K), function(v)
    paste(sprintf("%a", c(cert$view_g[v], cert$view_K[[v]]) + 0), collapse = ","), character(1L))
  kv <- match(key, unique(key))
  ko <- paste(kv, cert$view_operator)
  kov <- match(ko, unique(ko))
  per_link <- tabulate(kv)[kv[cert$row_view]] - tabulate(kov)[kov[cert$row_view]]
  sum(as.double(per_link))
}

# Shared certificate body. `links` are validated point links (above) or the
# outcome-free view descriptors of `certify_ipw_structural_classes_from_views`;
# both carry frame_ref, panel_ref, cohorts$g/K, pair$s/t, nuisance_view_ref,
# fitted_nuisance_ref, link_view_ref, glm_control and ref. Every mask, nuisance,
# vertex/operator, class, forest and E step below is the accepted certificate.
# CP08 scaling action 2: with a bound registry `mask_certificate`, the link's
# sealed completed verification replaces its re-streamed mask, and two
# distinct views of one operator are unequal by the certificate; every other
# comparison still streams. Classes, forests, E and the cycle check are unchanged.
# Shared design vertex/operator checks of the structural certificate (moved
# verbatim out of `ipw_view_structural_classes_core`; CP08 compact basis). Pure
# in the per-link descriptors (g, s, t, ordered K, functional key) and the
# design: every refusal, its class and its order are those of the view core.
# Returns the validated design row map and endpoint vertices for the links.
ipw_structural_design_checks <- function(link_g, link_s, link_t, link_K, link_fk, design, design_link_index, n_links) {
  # Validate intrinsic design vertex and link mapping using unchanged helpers in R/chained_ipw_shared.R
  dli <- validate_ipw_design_for_classes(design, n_links, design_link_index)
  DL <- design$links
  DV <- design$vertices
  cspec <- design$calendar_spec
  has_cal <- !is.null(cspec)
  is_pt_all <- FALSE
  vs <- DL$start_vertex_index[dli]
  vt <- DL$end_vertex_index[dli]
  fk <- DL$functional_key[dli]
  map_bad <- function(msg) cdid_abort(msg, c("cdid_ipw_class_vertex_error", "cdid_ipw_certificate_error"))
  if (has_cal) {
    if (!inherits(cspec, "cdid_calendar_spec") || !is.list(cspec)) {
      map_bad("Design `calendar_spec` is not a valid `cdid_calendar_spec` object.")
    }
    cspec_val <- tryCatch(
      validate_calendar_spec(
        cspec$calendar, cspec$cohorts, cspec$anticipation,
        pt_assumption = if (is.null(cspec$pt_assumption)) "post" else cspec$pt_assumption,
        pt_start = cspec$pt_start
      ),
      cdid_error = function(e) map_bad(conditionMessage(e))
    )
    if (!identical(cspec, cspec_val)) {
      map_bad("Design `calendar_spec` differs from its canonical validated representation.")
    }
    is_pt_all <- identical(cspec$pt_assumption, "all")
  }
  # One typed parse per distinct design pool string (CP08 scaling action 1).
  dl_pools <- if (has_cal && "control_pool" %in% names(DL)) design_pool_cache(DL$control_pool) else NULL

  for (ell in seq_len(n_links)) {
    j <- dli[ell]
    if (!identical(as.double(DL$g[j]), link_g[ell]) ||
        !identical(as.double(DL$s[j]), link_s[ell]) ||
        !identical(as.double(DL$t[j]), link_t[ell])) {
      map_bad("Design row cohort/dates differ from the link they are mapped to.")
    }
    if (has_cal) {
      delta <- cspec$anticipation
      K_ell <- link_K[[ell]]
      ref_clean <- tryCatch(
        all(vapply(K_ell, function(k) exact_cmp_sum(DL$ref[j], delta, k) < 0L, logical(1L))),
        cdid_error = function(e) map_bad(conditionMessage(e))
      )
      if (any(!K_ell %in% cspec$cohorts) || any(K_ell == link_g[ell]) || !ref_clean) {
        map_bad("Link control pool K is outside the design cohort universe or not clean at the cohort reference period.")
      }
      if ("control_pool" %in% names(DL)) {
        cp_vals <- dl_pools$values(j)
        if (is.null(cp_vals) || !setequal(as.double(cp_vals), K_ell)) {
          map_bad("Design row `control_pool` differs from the link's actual control pool K.")
        }
      }
    }
    for (side in 1:2) {
      v <- if (side == 1L) vs[ell] else vt[ell]
      u <- if (side == 1L) link_s[ell] else link_t[ell]
      qual_col <- if (side == 1L) "s_qualified" else "t_qualified"
      role_col <- if (side == 1L) "s_coordinate_role" else "t_coordinate_role"
      if (!has_cal) {
        qualified <- if (qual_col %in% names(DL)) DL[[qual_col]][j] else (u == DL$ref[j] || u >= link_g[ell])
        role <- if (u == DL$ref[j]) "reference" else if (u < link_g[ell]) "pre" else if (qualified) "post" else "post_nuisance"
        key <- if (role %in% c("pre", "post_nuisance")) fk[ell] else NA_character_
        expected_anchored <- (u == DL$ref[j])
      } else {
        cmp_ug <- tryCatch(exact_cmp_sum(u, delta, link_g[ell]), cdid_error = function(e) map_bad(conditionMessage(e)))
        qualified <- if (qual_col %in% names(DL)) DL[[qual_col]][j] else (u == DL$ref[j] || cmp_ug >= 0L)
        if (cmp_ug >= 0L && qualified) {
          u_clean <- tryCatch(
            exact_cmp_sum(u, delta, cspec$max_cohort) < 0L &&
              all(vapply(K_ell, function(k) exact_cmp_sum(u, delta, k) < 0L, logical(1L))),
            cdid_error = function(e) map_bad(conditionMessage(e))
          )
          if (!u_clean) {
            map_bad("Qualified causal endpoint is not clean under the link's actual control pool K and calendar anticipation.")
          }
        }
        role <- if (u == DL$ref[j]) {
          "reference"
        } else if (cmp_ug < 0L) {
          "pre"
        } else if (u < link_g[ell]) {
          if (qualified) "anticipation" else "anticipation_nuisance"
        } else {
          if (qualified) "post" else "post_nuisance"
        }
        key <- if (role %in% c("pre", "anticipation_nuisance", "post_nuisance")) fk[ell] else NA_character_
        expected_anchored <- (u == DL$ref[j]) || (is_pt_all && cmp_ug < 0L && u >= cspec$pt_start)
        if (role_col %in% names(DL)) {
          rc_val <- DL[[role_col]][j]
          ok_rc <- if (u == DL$ref[j]) {
            identical(rc_val, "reference")
          } else if (cmp_ug < 0L) {
            rc_val %in% c("pre", "pre_nuisance") || (is_pt_all && u >= cspec$pt_start && identical(rc_val, "pre_pt_all_zero"))
          } else if (u < link_g[ell]) {
            if (qualified) rc_val %in% c("anticipation", "anticipation_att") else identical(rc_val, "anticipation_nuisance")
          } else {
            if (qualified) rc_val %in% c("post", "post_att") else identical(rc_val, "post_nuisance")
          }
          if (!isTRUE(ok_rc)) {
            map_bad("Design row endpoint coordinate role disagrees with its calendar region or PT-all window.")
          }
        }
      }
      if (!identical(as.double(DV$time[v]), u) || !identical(as.double(DV$g[v]), link_g[ell]) ||
          !identical(DV$role[v], role) || !identical(DV$functional_key[v], key) ||
          !identical(is.na(DV$param_index[v]), expected_anchored)) {
        map_bad("Endpoint vertex is not the intrinsic qualified vertex of the link's cohort, date and operator.")
      }
    }
    if (vs[ell] == vt[ell]) map_bad("Loop link in design.")
    if (!is.null(link_fk[[ell]]) && !is.na(link_fk[[ell]]) && !identical(link_fk[[ell]], fk[ell])) {
      cdid_abort(
        "Link functional key differs from its design operator key.",
        c("cdid_ipw_class_operator_error", "cdid_ipw_certificate_error")
      )
    }
  }

  validate_ipw_design_incidence(design)
  list(dli = dli, vs = vs, vt = vt, fk = fk)
}

# Shared class/forest/expansion assembly of the structural certificate (moved
# verbatim out of `ipw_view_structural_classes_core`; CP08 compact basis). The
# per-link facts arrive established by the caller: dates and design vertices,
# the verified nuisance reference group `nu_gid` over `nu_refs`, each link's
# cohorts and glm control, and `same_mask_fn(ell, rep_ell)` deciding equality
# of link ell's mask with source representative rep_ell's nuisance mask from
# the caller's actual evidence. Classes, forests, E, the cross-class cycle
# check, the class records and the sealed certificate are unchanged.
ipw_structural_classes_assemble <- function(n_links, link_s, link_t, fk, vs, vt, dli, nu_gid, nu_refs,
                                            link_cohort_g, link_cohort_K, link_glm, same_mask_fn, b,
                                            link_kind, object_refs, context_ref) {
  map_bad <- function(msg) cdid_abort(msg, c("cdid_ipw_class_vertex_error", "cdid_ipw_certificate_error"))
  # Classes: (verified nuisance_view_ref, design operator functional_key)
  fk_levels <- sort(unique(fk), method = "radix")
  fk_id <- match(fk, fk_levels)
  ord <- order(nu_gid, fk_id, method = "radix")
  starts <- which(c(TRUE, diff(nu_gid[ord]) != 0L | diff(fk_id[ord]) != 0L))
  n_classes <- length(starts)
  if (n_classes > b$max_classes) {
    cdid_abort("Number of structural classes exceeds `max_classes`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  ends <- c(starts[-1L] - 1L, n_links)
  plans <- vector("list", n_classes)
  total_E <- 0

  for (ci in seq_len(n_classes)) {
    members <- ord[seq.int(starts[ci], ends[ci])]
    members <- members[order(link_s[members], link_t[members], method = "radix")]
    dates <- sort(unique(c(link_s[members], link_t[members])))
    if (length(dates) > b$max_class_dates) {
      cdid_abort("Class date count exceeds `max_class_dates`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
    }
    a_idx <- match(link_s[members], dates)
    b_idx <- match(link_t[members], dates)
    if (anyDuplicated(paste(a_idx, b_idx))) {
      cdid_abort(
        "Duplicate (s, t) link within one structural class.",
        c("cdid_ipw_class_duplicate_error", "cdid_ipw_certificate_error")
      )
    }
    vmap <- rep(NA_integer_, length(dates))
    for (side in 1:2) {
      di <- if (side == 1L) a_idx else b_idx
      vv <- if (side == 1L) vs[members] else vt[members]
      for (h in seq_along(di)) {
        if (is.na(vmap[di[h]])) vmap[di[h]] <- vv[h] else if (vmap[di[h]] != vv[h]) {
          map_bad("Class links map one date to different qualified vertices; no common vertex map.")
        }
      }
    }
    if (anyDuplicated(vmap)) map_bad("Class vertex map is not injective over dates.")
    forest <- incidence_forest(length(dates), a_idx, b_idx, integer(0L))
    tree <- sort(forest$parent_link[forest$parent_link > 0L])
    dropped <- setdiff(seq_along(members), tree)
    depth <- ipw_forest_depth(forest)
    n_e <- 0
    for (h in dropped) n_e <- n_e + ipw_forest_path(a_idx[h], b_idx[h], forest, depth, TRUE)
    total_E <- total_E + n_e
    plans[[ci]] <- list(
      members = members, dates = dates, a_idx = a_idx, b_idx = b_idx,
      vmap = vmap, forest = forest, depth = depth, tree = tree, dropped = dropped
    )
  }
  if (total_E > b$max_E_entries) {
    cdid_abort("Reduction map entries exceed `max_E_entries`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }

  # Cross-class date-cycle check among retained links sharing identical (g, K, glm_control, mask)
  # Two distinct nuisance_view_ref objects with the same (g, K, glm_control) are checked via
  # chunked zero-outcome `verify_ipw_view_mask_compatibility`.
  src_reps <- integer()
  src_id <- integer(n_links)
  for (ell in seq_len(n_links)) {
    hit_k <- 0L
    for (k in seq_along(src_reps)) {
      rep_ell <- src_reps[k]
      # Same verified nuisance reference group (identical `nuisance_view_ref`).
      if (nu_gid[ell] == nu_gid[rep_ell] &&
          identical(link_glm[[ell]], link_glm[[rep_ell]])) {
        hit_k <- k
        break
      }
      if (identical(link_cohort_g[[ell]], link_cohort_g[[rep_ell]]) &&
          identical(link_cohort_K[[ell]], link_cohort_K[[rep_ell]]) &&
          identical(link_glm[[ell]], link_glm[[rep_ell]])) {
        # Decided by the caller's actual mask evidence for (link ell, source
        # representative rep_ell): the view core streams or reads its bound
        # registry certificate; the compact route reads complete-date posting facts.
        same_mask <- same_mask_fn(ell, rep_ell)
        if (isTRUE(same_mask)) {
          hit_k <- k
          break
        }
      }
    }
    if (hit_k == 0L) {
      src_reps <- c(src_reps, ell)
      hit_k <- length(src_reps)
    }
    src_id[ell] <- hit_k
  }

  dropped_all <- unlist(lapply(plans, function(pl) pl$members[pl$dropped]), use.names = FALSE)
  retained_all <- setdiff(seq_len(n_links), dropped_all)
  for (k in seq_along(src_reps)) {
    rk <- retained_all[src_id[retained_all] == k]
    if (length(rk) < 2L) next
    dk <- sort(unique(c(link_s[rk], link_t[rk])))
    parent <- seq_along(dk)
    find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
    for (ell in rk) {
      x <- find(match(link_s[ell], dk))
      y <- find(match(link_t[ell], dk))
      if (x == y) {
        cdid_abort(
          "Retained links with identical nuisance inputs (same object across operators, duplicate objects or ordinary refits) form a date cycle; this singularity is not certified by any class.",
          c("cdid_ipw_reduction_singular_error", "cdid_ipw_certificate_error")
        )
      }
      parent[x] <- y
    }
  }

  E_i <- E_j <- E_x <- integer(total_E)
  pos <- 0L
  classes <- vector("list", n_classes)
  for (ci in seq_len(n_classes)) {
    pl <- plans[[ci]]
    for (h in pl$dropped) {
      path <- ipw_forest_path(pl$a_idx[h], pl$b_idx[h], pl$forest, pl$depth, FALSE)
      o <- order(pl$members[path$edges])
      idx <- pos + seq_along(o)
      E_i[idx] <- pl$members[h]
      E_j[idx] <- pl$members[path$edges[o]]
      E_x[idx] <- path$signs[o]
      pos <- pos + length(o)
    }
    first <- pl$members[1L]
    nv_first <- get("nuisance_view", envir = nu_refs[[nu_gid[first]]], inherits = FALSE)
    classes[[ci]] <- list(
      nuisance_id = nv_first$nuisance_id,
      functional_key = fk[first],
      g = nv_first$cohorts$g,
      K = nv_first$cohorts$K,
      dates = pl$dates,
      vertex_index = pl$vmap,
      link_index = pl$members,
      forest_link_index = pl$members[pl$tree],
      dropped_link_index = pl$members[pl$dropped],
      roots = pl$dates[pl$forest$roots],
      n_components = pl$forest$n_components,
      cycle_rank = length(pl$members) - (length(pl$dates) - pl$forest$n_components),
      n_risk_a_units = nv_first$support$n_risk_a_units,
      n_observers = nv_first$support$n_observed_risk_units,
      n_pair_nonobservers = nv_first$support$n_nonobserved_risk_units,
      n_outside_risk_a_units = nv_first$support$n_outside_risk_a_units
    )
  }

  ord_e <- order(E_i, E_j, method = "radix")
  dropped_index <- sort(as.integer(dropped_all))
  certificate <- structure(
    list(
      weighting = "optimal",
      n_links = n_links,
      link_kind = rep(link_kind, n_links),
      retained_index = setdiff(seq_len(n_links), dropped_index),
      dropped_index = dropped_index,
      E = list(i = E_i[ord_e], j = E_j[ord_e], x = E_x[ord_e]),
      classes = classes,
      design_link_index = dli,
      nuisance_refs = nu_refs,
      object_refs = object_refs,
      context_ref = context_ref,
      bounds = b,
      retained_covariance = "not_certified: PD/accuracy of the retained covariance must be validated by the caller (Cor. 6); otherwise refuse",
      scope = "exact real-algebra structural identities of Prop. 5 / Cor. 6 for optimal GMM; no numerical acceptance, sampling or catalog claim"
    ),
    class = "cdid_ipw_structural_certificate"
  )
  ipw_seal(certificate, "certificate")
}

ipw_view_structural_classes_core <- function(links, design, design_link_index, b, own_ctrl, is_custom, link_kind,
                                             mask_certificate = NULL) {
  n_links <- length(links)
  frame_owner <- get("frame_owner", envir = links[[1L]]$frame_ref, inherits = FALSE)
  n_a_units <- as.integer(frame_owner$n_a_units)
  mask_work <- as.double(n_links) * as.double(n_a_units)
  if (mask_work > b$max_mask_entries) {
    cdid_abort("Mask verification work exceeds `max_mask_entries`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }

  link_g <- link_s <- link_t <- numeric(n_links)
  link_fk <- vector("list", n_links)
  link_K <- vector("list", n_links)
  for (ell in seq_len(n_links)) {
    pl <- links[[ell]]
    lv <- get("link_view", envir = pl$link_view_ref, inherits = FALSE)
    link_g[ell] <- as.double(pl$cohorts$g)
    link_s[ell] <- as.double(pl$pair$s)
    link_t[ell] <- as.double(pl$pair$t)
    link_fk[ell] <- list(lv$functional_key)
    link_K[[ell]] <- as.double(pl$cohorts$K)
  }

  checks <- ipw_structural_design_checks(link_g, link_s, link_t, link_K, link_fk, design, design_link_index, n_links)
  dli <- checks$dli; vs <- checks$vs; vt <- checks$vt; fk <- checks$fk

  cert_view <- if (is.null(mask_certificate)) NULL else ipw_mask_certificate_bind(mask_certificate, links, frame_owner)

  # Group links by shared fitted nuisance reference and verify chunked mask compatibility.
  # The nuisance and fitted snapshots behind one exact immutable reference pair
  # are validated once per invocation: a reference seen before with the same
  # fitted reference is the identical locked snapshot already validated against
  # this frame owner, so re-running the same check on the same object is
  # removed. A known nuisance reference arriving with a different fitted
  # reference is a new pair: its fitted snapshot is validated first, as before,
  # and then refused as a conflict. Every link still runs its own seal, owner,
  # panel, cohort, pair, binding and completed-mask checks; the call-local
  # checked set lets the per-link validator skip only the recursive validation
  # of a nuisance snapshot this loop validated itself. Its scratch (two slots
  # per distinct nuisance reference plus three headers, alongside the two
  # reference-group slots per link) is charged against the explicit controls
  # before it is allocated and is dropped before the class computation.
  ipw_owner_safe_sum(c(ipw_owner_safe_product(4, n_links, label = "structural reference groups and checked nuisance set"), 4),
    min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries),
    "structural-class reference groups plus call-local checked nuisance set")
  nu_refs <- list()
  fit_refs <- list()
  nu_gid <- integer(n_links)
  nuisance_checked <- if (is.null(cert_view)) NULL else ipw_link_view_nuisance_checked_set(frame_owner)
  for (ell in seq_len(n_links)) {
    pl <- links[[ell]]
    r_nu <- pl$nuisance_view_ref
    r_fit <- pl$fitted_nuisance_ref
    nv_obj <- get("nuisance_view", envir = r_nu, inherits = FALSE)
    lv_obj <- get("link_view", envir = pl$link_view_ref, inherits = FALSE)
    hit <- 0L
    for (k in seq_along(nu_refs)) {
      if (identical(nu_refs[[k]], r_nu)) {
        hit <- k
        break
      }
    }
    nu_new <- hit == 0L
    pair_new <- nu_new || !identical(fit_refs[[hit]], r_fit)
    # Group membership alone is not a validation receipt: an unlocked or
    # foreign same-label container can return a different object on a later
    # read. Only the actual immutable aliased snapshot can reuse its check.
    nu_cached <- !nu_new && inherits(nv_obj, "cdid_ipw_nuisance_view") &&
      ipw_owner_immutable_view_snapshot(r_nu, nv_obj, "nuisance_view")
    if (!nu_cached) {
      validate_ipw_nuisance_view(nv_obj, frame_owner = frame_owner)
      ipw_link_view_nuisance_checked_add(nuisance_checked, r_nu, nv_obj, frame_owner$ref)
    }
    if (!is_custom) {
      fv_obj <- get("fitted_nuisance_view", envir = r_fit, inherits = FALSE)
      fit_cached <- !pair_new && nu_cached && inherits(fv_obj, "cdid_ipw_fitted_nuisance_view") &&
        ipw_owner_immutable_view_snapshot(r_fit, fv_obj, "fitted_nuisance_view")
      if (!fit_cached) validate_ipw_fitted_nuisance_view(fv_obj, nuisance_view = nv_obj, frame_owner = frame_owner)
    }
    # For the custom path this is the original IPW callback's semantic contract:
    # identical weighted design/mask gives the same nuisance, and normalized
    # outcome means and complete IFs are linear in endpoint outcomes. Shape
    # validation does not establish that contract; no built-in fit is substituted.
    if (is.null(cert_view)) {
      verify_ipw_view_mask_compatibility(nuisance_view = nv_obj, link_view = lv_obj, control = own_ctrl)
    } else {
      # The sealed link view records a completed streamed verification bound to
      # exactly this nuisance view snapshot; it is checked, not re-streamed.
      ipw_link_view_validate_body(lv_obj, frame_owner = frame_owner, nuisance_view = nv_obj,
        nuisance_checked = nuisance_checked)
    }

    if (!nu_new && pair_new) {
      cdid_abort(
        "Links sharing a `nuisance_view_ref` do not share the exact same `fitted_nuisance_ref`.",
        c("cdid_ipw_shared_fit_mismatch_error", "cdid_ipw_shared_fit_error", "cdid_ipw_certificate_error")
      )
    }
    if (nu_new) {
      nu_refs[[length(nu_refs) + 1L]] <- r_nu
      fit_refs[length(nu_refs)] <- list(r_fit)
      hit <- length(nu_refs)
    }
    nu_gid[ell] <- hit
  }
  nuisance_checked <- NULL

  # Equality of link ell's mask with the nuisance mask of source representative
  # rep_ell: a bound registry certificate proves two distinct views of one
  # operator unequal; every other comparison streams exactly as before.
  same_mask_fn <- function(ell, rep_ell) {
    pl <- links[[ell]]; pl_rep <- links[[rep_ell]]
    if (!is.null(cert_view) &&
        identical(mask_certificate$view_operator[cert_view[ell]], mask_certificate$view_operator[cert_view[rep_ell]])) {
      # Two distinct registry views of one operator: their masks were proved
      # unequal when the later view was created (fingerprint or exact stream).
      return(FALSE)
    }
    lv_ell <- get("link_view", envir = pl$link_view_ref, inherits = FALSE)
    nv_rep <- get("nuisance_view", envir = pl_rep$nuisance_view_ref, inherits = FALSE)
    tryCatch(
      {
        verify_ipw_view_mask_compatibility(nuisance_view = nv_rep, link_view = lv_ell, control = own_ctrl)
        TRUE
      },
      cdid_ipw_shared_mask_error = function(e) FALSE
    )
  }
  ipw_structural_classes_assemble(n_links, link_s, link_t, fk, vs, vt, dli, nu_gid, nu_refs,
    link_cohort_g = lapply(links, function(pl) pl$cohorts$g), link_cohort_K = lapply(links, function(pl) pl$cohorts$K),
    link_glm = lapply(links, `[[`, "glm_control"), same_mask_fn = same_mask_fn, b = b,
    link_kind = link_kind, object_refs = lapply(links, `[[`, "ref"), context_ref = links[[1L]]$frame_ref)
}

# Reserve the early certificate while the actual fitted nuisance/view registry
# and the accepted whole-pipeline envelope are live. This is a resource bound,
# never a substitute for the shared core's scientific/integrity checks. A path
# contains at most min(T - 1, M) edges, so M times that length bounds all E entries
# before any descriptor, forest, expansion or dropped-record allocation.
preflight_ipw_early_structural_reduction <- function(registry, design, bounds, owner_control,
                                                   retained_pipeline_entries = 0, external_work_ops = 0) {
  b <- validate_ipw_certificate_bounds(bounds)
  ctrl <- validate_ipw_owner_control(owner_control)
  retained_pipeline_entries <- ipw_owner_int_bound(retained_pipeline_entries, "retained_pipeline_entries", min_val = 0)
  external_work_ops <- ipw_owner_int_bound(external_work_ops, "external_work_ops", min_val = 0)
  M <- as.double(length(registry$link_views))
  K <- as.double(length(registry$fitted_views))
  if (M < 1 || K < 1 || M > b$max_links || !inherits(design, "cdid_chained_design") ||
      !is.data.frame(design$vertices)) {
    cdid_abort("Invalid early structural preflight dimensions.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  cap <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
  V <- as.double(nrow(design$vertices))
  # Reserve the scan/descriptor headers before examining registry contents.
  linear <- ipw_owner_safe_sum(c(256 * M, 64 * V, 256), cap, "early certificate linear plan")
  ipw_owner_safe_sum(c(retained_pipeline_entries, linear), cap, "pipeline plus early certificate headers")
  # Call-local checked set of the nuisance snapshots validated through their
  # fitted views (two slots per fitted view plus headers, charged before use);
  # each link then validates its own seal, owner, binding and completed mask
  # while its recursive nuisance validation runs once per sealed reference.
  # The wrapper's flags/set coexist with the core's reference groups/set.
  # Include their maximum simultaneous payload in the returned whole-pipeline
  # reservation, not only in an isolated local fit-to-cap check.
  validation_scratch <- ipw_owner_safe_sum(c(ipw_owner_safe_product(4, M, label = "structural reference groups and checked nuisance set"),
    ipw_owner_safe_product(3, K, label = "checked fitted views and nuisance set"), 12), cap,
    "simultaneous structural validation scratch")
  ipw_owner_safe_sum(c(retained_pipeline_entries, linear, validation_scratch),
    cap, "pipeline plus early certificate headers and checked nuisance set")
  fitted_entries <- view_entries <- 0
  nuisance_checked <- NULL
  for (fv in registry$fitted_views) {
    validate_ipw_fitted_nuisance_view(fv)
    fitted_entries <- ipw_owner_safe_sum(c(fitted_entries, fv$accounting$total_retained_entries), cap,
      "actual fitted nuisances already live before the early certificate")
    if (is.null(nuisance_checked)) {
      nuisance_checked <- ipw_link_view_nuisance_checked_set(get("frame_owner", envir = fv$frame_ref, inherits = FALSE))
    }
    nv_obj <- get("nuisance_view", envir = fv$nuisance_view_ref, inherits = FALSE)
    ipw_link_view_nuisance_checked_add(nuisance_checked, fv$nuisance_view_ref, nv_obj, nv_obj$frame_ref)
  }
  for (lv in registry$link_views) {
    ipw_link_view_validate_body(lv, nuisance_checked = nuisance_checked)
    view_entries <- ipw_owner_safe_sum(c(view_entries, lv$preflight$total_retained_entries), cap,
      "all original link views live before the early certificate")
  }
  nuisance_checked <- NULL
  fv <- registry$fitted_views[[1L]]
  frame <- get("frame_owner", envir = fv$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", envir = fv$panel_ref, inherits = FALSE)
  T <- as.double(panel$n_calendar_dates)
  C <- as.double(length(panel$cohorts))
  N <- as.double(frame$n_a_units)
  path_bound <- min(M, T - 1, b$max_class_dates - 1)
  E_bound <- min(b$max_E_entries, ipw_owner_safe_product(M, path_bound, label = "early E entry upper bound"))
  # P1's registry reservation is composed, rather than replacing it with a
  # certificate-only allowance. Fitted accounting includes its nuisance view
  # and character metadata. The original-row envelope remains charged as well.
  registry_entries <- 16 * (M + 1) + C
  # Live split records and the registry's mask-equivalence certificate.
  split_entries <- ipw_owner_safe_product(8, length(registry$mask_splits), label = "live mask split records")
  cert <- registry$mask_certificate
  cert_entries <- if (is.null(cert)) 0 else
    ipw_owner_safe_product(M, 11 + C, label = "live mask-equivalence certificate")
  live_registry <- ipw_owner_safe_sum(c(fitted_entries, view_entries, registry_entries, split_entries, cert_entries), cap,
    "actual fitted/view objects plus accepted P1 registry reservation")
  certificate_retained <- ipw_owner_safe_sum(c(5 * E_bound, linear, 4 * M * C), cap,
    "certificate E, two dropped-record expansion vectors, full-row maps and classes")
  ipw_owner_safe_sum(c(live_registry, certificate_retained), ctrl$max_output_entries,
    "early registry, certificate and reduced-link output bound")
  # 16 E slots cover construction/sort/copy buffers as well as the five
  # persistent E slots; linear includes all class plans and descriptors.
  mask_scratch <- 64 * min(N, ctrl$max_chunk_rows) + 32 * C + 256
  additional_entries <- ipw_owner_safe_sum(c(live_registry, linear, 4 * M * C, 16 * E_bound, mask_scratch, validation_scratch), cap,
    "early certificate construction, retained expansion records and live fitted registry")
  ipw_owner_safe_sum(c(retained_pipeline_entries, additional_entries), cap,
    "whole pipeline plus early certificate before allocation")
  # Without a registry certificate the core verifies M masks, then compares at
  # most M*K source representatives. With one it re-streams only the open
  # comparisons it cannot resolve (identical cohorts across operators) and
  # checks M sealed bindings plus the certificate's view keys. Reuse the
  # accepted streamed-mask work formula; include forest paths, design checks,
  # sorting and binary searches for the reduced-link expansion slices.
  mask_work <- 128 * N * (1 + ceiling(log2(C + 1))) +
    128 * ceiling(N / min(N, ctrl$max_chunk_rows)) * (1 + ceiling(log2(N + 1)))
  open_streams <- if (is.null(cert)) NA_real_ else ipw_mask_certificate_open_streams(cert)
  mask_streams <- if (is.na(open_streams)) M * (K + 1) else open_streams
  integrity_work <- if (is.na(open_streams)) 0 else 512 * M * (C + 2)
  work_ops_bound <- ipw_owner_safe_sum(c(mask_streams * mask_work, integrity_work,
    256 * M * (T + (K + 1) * (C + 1) + ceiling(log2(M + 1))), 64 * E_bound), ctrl$max_work_ops,
    "early structural mask/design/forest/expansion work")
  ipw_owner_safe_sum(c(external_work_ops, work_ops_bound), ctrl$max_work_ops,
    "whole pipeline plus early structural work before allocation")
  list(E_entries_upper = E_bound, fitted_nuisance_entries = fitted_entries,
    original_link_view_entries = view_entries, p1_registry_entries = registry_entries,
    mask_split_entries = split_entries, mask_certificate_entries = cert_entries,
    mask_streams_charged = mask_streams,
    mask_evidence = if (is.na(open_streams)) "streamed: every link mask and source-representative comparison" else
      "registry mask-equivalence certificate; only cross-operator comparisons with identical cohorts stream",
    certificate_retained_entries = certificate_retained, additional_entries = additional_entries,
    work_ops_bound = work_ops_bound)
}

# CP08 P3: the same certificate from outcome-free objects, before any point fit
# (cf. shared declarations in `certify_ipw_structural_classes`). Row m is its
# sealed link view, bound with a completed mask verification to the nuisance view
# whose single fitted nuisance it reuses (registry of
# `build_pipeline_nuisance_and_link_views(fit_point_links = FALSE)`). The
# point-link wrapper's checks that need no fitted outcome functional are made
# here; the shared core then recomputes masks, vertex/operator maps, classes,
# forests and E. Licensed only for built-in baseline IPW under optimal weighting.
certify_ipw_structural_classes_from_views <- function(registry, design, design_link_index, bounds, owner_control,
                                                      retained_pipeline_entries = 0, external_work_ops = 0) {
  b <- validate_ipw_certificate_bounds(bounds)
  own_ctrl <- validate_ipw_owner_control(owner_control)
  if (!is.list(registry)) registry <- list()
  lvs <- registry$link_views
  fvs <- registry$fitted_views
  nvs <- registry$nuisance_views
  idx <- registry$link_to_nuisance_idx
  if (!is.list(lvs) || is.object(lvs) || length(lvs) < 1L || !is.list(fvs) || !is.list(nvs) ||
      length(fvs) != length(nvs) || !is.numeric(idx) || length(idx) != length(lvs) || anyNA(idx) ||
      any(idx < 1 | idx > length(fvs) | idx != floor(idx)) || !is.null(registry$point_links)) {
    cdid_abort(
      "Outcome-free structural certification needs an unfitted view registry with one fitted nuisance per view.",
      c("cdid_ipw_input_error", "cdid_ipw_certificate_error")
    )
  }
  n_links <- length(lvs)
  if (n_links > b$max_links) {
    cdid_abort("Number of links exceeds `max_links`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  preflight_ipw_early_structural_reduction(registry, design, b, own_ctrl,
    retained_pipeline_entries, external_work_ops)
  frame_ref <- fvs[[idx[1L]]]$frame_ref
  if (!is.environment(frame_ref) || !exists("frame_owner", envir = frame_ref, inherits = FALSE)) {
    cdid_abort("Fitted nuisance view lacks a sealed frame reference.", c("cdid_ipw_context_mismatch_error", "cdid_ipw_certificate_error"))
  }
  frame_owner <- get("frame_owner", envir = frame_ref, inherits = FALSE)
  # Each fitted view and its nuisance are validated once; the nuisance sealed
  # reference is then recorded in a call-local checked set (two slots per view
  # plus headers, charged here next to the existing per-view flags) so each
  # link's validation keeps its own seal, owner, binding and completed-mask
  # checks without re-validating the same nuisance snapshot per link.
  ipw_owner_safe_sum(c(retained_pipeline_entries, ipw_owner_safe_product(3, length(fvs), label = "checked fitted views and nuisance set"), 4),
    min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries),
    "pipeline envelope plus structural checked-view scratch")
  checked <- logical(length(fvs))
  nuisance_checked <- ipw_link_view_nuisance_checked_set(frame_owner)
  links <- vector("list", n_links)
  for (ell in seq_len(n_links)) {
    k <- idx[ell]
    fv <- fvs[[k]]
    nv <- nvs[[k]]
    lv <- lvs[[ell]]
    if (!checked[k]) {
      validate_ipw_fitted_nuisance_view(fv, nuisance_view = nv, frame_owner = frame_owner)
      checked[k] <- TRUE
      if (is.environment(nv$ref)) ipw_link_view_nuisance_checked_add(nuisance_checked, nv$ref, nv, nv$frame_ref)
    }
    ipw_link_view_validate_body(lv, frame_owner = frame_owner, nuisance_view = nv, nuisance_checked = nuisance_checked)
    if (!isTRUE(lv$nuisance_mask_verified) || !identical(fv$frame_ref, frame_ref) ||
        !identical(lv$panel_ref, fv$panel_ref) ||
        !identical(as.double(lv$cohorts$g), as.double(fv$cohorts$g)) ||
        !identical(as.double(lv$cohorts$K), as.double(fv$cohorts$K))) {
      cdid_abort(
        paste0("Link view ", ell, " is not a verified outcome-free sharer of its fitted nuisance."),
        c("cdid_ipw_shared_mask_error", "cdid_ipw_certificate_error")
      )
    }
    links[[ell]] <- list(frame_ref = fv$frame_ref, panel_ref = fv$panel_ref, cohorts = fv$cohorts,
      pair = list(s = lv$pair$s, t = lv$pair$t), nuisance_view_ref = fv$nuisance_view_ref,
      fitted_nuisance_ref = fv$ref, link_view_ref = lv$ref, glm_control = fv$glm_control, ref = lv$ref)
  }
  # The registry's constructor-owned certificate (when present) is bound and
  # consumed by the core; a registry without one streams every comparison.
  ipw_view_structural_classes_core(links, design, design_link_index, b, own_ctrl, FALSE, "shared_view_declaration",
    mask_certificate = registry$mask_certificate)
}

# Outcome-free identity of a dropped original row on the early route: its link
# view, the fitted nuisance it shares and its certified expansion row. Not a
# point link: no estimate, influence function or fitted object exists for it.
ipw_pipeline_reduced_link_record <- function(registry, m, certificate) {
  lv <- registry$link_views[[m]]
  fv <- registry$fitted_views[[registry$link_to_nuisance_idx[m]]]
  # E is issued in (original row, retained row) order. Binary-search the
  # contiguous row, avoiding an E-length logical allocation/scan per dropped row.
  n_e <- length(certificate$E$i)
  first <- ipw_owner_lower_bound(certificate$E$i, 1, n_e, m)
  last <- ipw_owner_upper_bound(certificate$E$i, first, n_e, m)
  hit <- if (last >= first) seq.int(first, last) else integer(0L)
  structure(list(
    link_id = lv$link_id,
    catalog_row = as.integer(m),
    link_view_ref = lv$ref,
    nuisance_view_ref = fv$nuisance_view_ref,
    fitted_nuisance_ref = fv$ref,
    expansion = list(retained_rows = certificate$E$j[hit], signs = certificate$E$x[hit]),
    point_fit_performed = FALSE,
    outcomes_read = FALSE,
    status = "structurally_reduced_outcome_free",
    scope = "Original comparison recovered through the certified expansion E (Prop. 5/Cor. 6); never fitted, replayed or read."
  ), class = "cdid_ipw_reduced_link_record")
}

# ------------------------------------------------------------------------------
# 2b. CP08 compact canonical basis: the early structural certificate from
#     catalog descriptors on the complete-roster built-in baseline IPW route
# ------------------------------------------------------------------------------
# The accepted outcome-free reduction (P3) certifies classes, forests and E
# from sealed link views, one per original catalog row, built before the
# certificate reduces M rows to the L forest rows. On the exact route below the
# same certificate is determined earlier, from the validated catalog's canonical
# row descriptors (cohort, ordered controls, dates, operator key), the one
# fitted nuisance of the operator and actual complete-date posting facts of the
# sealed panel: when every catalog date posts exactly the whole A roster, every
# row's risk-set pair-observation mask equals the nuisance view's mask by the
# same posting fact the link constructor uses (`ipw_view_mask_compatibility_fact`,
# `ipw_owner_pair_masks_equal_fact`), so the registry's per-row streamed or
# fact-decided verification would return exactly this result for every row.
# Heavy link views are then constructed only for the L retained rows (and the
# auxiliary comparisons the reporting stage always builds); every dropped row
# keeps a sealed canonical descriptor and an outcome-free reduced-link record.
# The design/vertex/operator checks, class grouping, BFS forest, expansion E and
# cross-class cycle check are the shared functions above, unchanged; the
# certificate, kept indices, E, the original-M metric, the global nuisance IF
# and every downstream output are those of the accepted route.
#
# Route predicate (decided from actual sealed facts before any fit; every other
# request keeps its current route unchanged, and a refusal inside this route
# is never caught to retry another): built-in baseline IPW (no CP04 backend, no
# callback), optimal weighting with outcome-free reduction, PT-post without
# numerical transfer, exactly one catalog operator whose target is one finite
# cohort and whose control pool is the never-treated cohort alone, a frame whose
# supplied roster is exactly the A roster (no singleton, unobserved or
# frame-only unit), one fixed n_A x p double nuisance design with constant unit
# weights, and every catalog date complete on the sealed postings.

ipw_pipeline_compact_basis_facts <- function(metadata_obj, frame_owner, owner_control, pt_assumption,
                                             run_numerical_transfer, n_total) {
  if (!is.list(metadata_obj) || !inherits(metadata_obj, "cdid_pipeline_metadata") ||
      !is.list(frame_owner) || !inherits(frame_owner, "cdid_ipw_frame_owner") ||
      !identical(pt_assumption, "post") || isTRUE(run_numerical_transfer) ||
      !is.list(owner_control) || !gmm_scalar_numeric(owner_control$max_chunk_rows, positive = TRUE)) return(NULL)
  cat_obj <- metadata_obj$catalog
  rows <- cat_obj$rows
  ops <- cat_obj$operators
  pools <- cat_obj$pools
  if (!is.data.frame(rows) || nrow(rows) < 1L || !is.data.frame(ops) || nrow(ops) != 1L || !is.list(pools)) return(NULL)
  M <- nrow(rows)
  op_id <- ops$operator_id[1L]
  if (!is.numeric(op_id) || length(op_id) != 1L || is.na(op_id) || !isTRUE(all(rows$operator_id == op_id))) return(NULL)
  pid <- ops$pool_id[1L]
  if (!is.numeric(pid) || length(pid) != 1L || is.na(pid) || pid < 1 || pid > length(pools$size) ||
      !identical(as.integer(pools$size[pid]), 1L) || !isTRUE(pools$never_treated[pid])) return(NULL)
  coh <- cat_obj$spec$cohorts
  cal <- cat_obj$spec$calendar
  k_idx <- pools$entries[pools$start[pid]]
  if (length(k_idx) != 1L || is.na(k_idx)) return(NULL)
  K <- as.double(coh[k_idx])
  if (length(K) != 1L || !identical(K, Inf)) return(NULL)
  g <- as.double(coh[ops$target_index[1L]])
  if (length(g) != 1L || !is.finite(g) || !isTRUE(all(as.double(rows$g) == g))) return(NULL)
  # Owners: the metadata's sealed frame and panel; the supplied roster is the A roster.
  if (!is.environment(frame_owner$panel_ref) || !identical(frame_owner$ref, metadata_obj$frame_ref) ||
      !identical(frame_owner$panel_ref, metadata_obj$panel_ref)) return(NULL)
  po <- get0("panel_owner", envir = frame_owner$panel_ref, inherits = FALSE)
  if (!inherits(po, "cdid_ipw_panel_owner") || !is.list(po)) return(NULL)
  n_a <- frame_owner$n_a_units
  if (!gmm_scalar_numeric(n_a, positive = TRUE) || !identical(as.double(po$n_a_units), as.double(n_a)) ||
      !identical(as.double(n_total), as.double(n_a)) ||
      !identical(as.double(frame_owner$n_total_units), as.double(n_a)) ||
      !identical(as.double(po$n_total_units), as.double(n_a)) ||
      !identical(as.double(frame_owner$n_frame_only_units), 0) ||
      !identical(as.double(po$n_singleton_units), 0) || !identical(as.double(po$n_unobserved_units), 0)) return(NULL)
  Z <- frame_owner$Z_a
  w <- frame_owner$w_a
  p <- frame_owner$n_covariates
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 1 || p != floor(p) ||
      !is.matrix(Z) || !is.double(Z) || !identical(dim(Z), c(as.integer(n_a), as.integer(p))) ||
      !is.double(w) || length(w) != n_a || !isTRUE(w[1L] == 1) || !isTRUE(all(w == w[1L]))) return(NULL)
  # Catalog dates index the panel calendar; every one of them posts the whole A roster.
  if (!identical(as.double(cal), as.double(po$all_times))) return(NULL)
  dates <- sort(unique(c(as.integer(rows$s_index), as.integer(rows$t_index))))
  if (!length(dates) || anyNA(dates) || any(dates < 1L | dates > as.integer(po$n_calendar_dates))) return(NULL)
  step <- min(as.double(n_a), as.double(owner_control$max_chunk_rows))
  memo <- list2env(list(panel_owner = po), parent = emptyenv())
  for (d in dates) if (!isTRUE(ipw_owner_date_complete(po, d, step, memo))) return(NULL)
  list(
    panel_owner = po, frame_ref = frame_owner$ref, panel_ref = frame_owner$panel_ref, metadata_ref = metadata_obj$ref,
    operator_id = as.integer(op_id), first_row = which(rows$operator_id == op_id)[1L], g = g, K = K,
    p = as.integer(p), n_a_units = as.integer(n_a), step = step, date_indices = dates, n_rows = as.integer(M),
    scope = "structural facts of the sealed owners and validated catalog only; no estimate, closeness or label is read"
  )
}

# Simultaneous reservation of the compact route before any allocation: the
# canonical descriptor records, the P1 registry headers, the operator's
# nuisance view and fitted model, the retained basis link views and row
# descriptors (at most one class on this route, hence at most
# min(M, T - 1, max_class_dates - 1) forest rows), the certificate (five
# persistent slots per E entry, class plans, maps) with its construction peak
# (sixteen slots per E entry), the constructors' transient scratch, the
# date-fact read scratch and the validation scratch, all against the unchanged
# caps next to the live pipeline envelope. The returned additional entries
# persist through the downstream guards exactly like the view route's
# reservation; they are smaller because the M - L heavy views are never built.
ipw_pipeline_compact_basis_preflight <- function(facts, design, bounds, owner_control,
                                                 retained_pipeline_entries = 0, external_work_ops = 0) {
  b <- validate_ipw_certificate_bounds(bounds)
  ctrl <- validate_ipw_owner_control(owner_control)
  retained_pipeline_entries <- ipw_owner_int_bound(retained_pipeline_entries, "retained_pipeline_entries", min_val = 0)
  external_work_ops <- ipw_owner_int_bound(external_work_ops, "external_work_ops", min_val = 0)
  if (!is.list(facts) || !inherits(facts$panel_owner, "cdid_ipw_panel_owner") ||
      !gmm_scalar_numeric(facts$n_rows, positive = TRUE) || !gmm_scalar_numeric(facts$p, positive = TRUE) ||
      !gmm_scalar_numeric(facts$n_a_units, positive = TRUE) || !gmm_scalar_numeric(facts$step, positive = TRUE)) {
    cdid_abort("Compact basis preflight needs the route facts of a sealed panel owner.",
      c("cdid_ipw_input_error", "cdid_ipw_certificate_error"))
  }
  M <- as.double(facts$n_rows)
  K <- 1
  if (M < 1 || M > b$max_links || !inherits(design, "cdid_chained_design") || !is.data.frame(design$vertices)) {
    cdid_abort("Invalid compact basis preflight dimensions.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  cap <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
  V <- as.double(nrow(design$vertices))
  panel <- facts$panel_owner
  T <- as.double(panel$n_calendar_dates)
  C <- as.double(length(panel$cohorts))
  N <- as.double(facts$n_a_units)
  p <- as.double(facts$p)
  D <- as.double(length(facts$date_indices))
  step <- min(N, as.double(ctrl$max_chunk_rows))
  # Same linear plan and validation scratch as the view route's preflight.
  linear <- ipw_owner_safe_sum(c(256 * M, 64 * V, 256), cap, "compact certificate linear plan")
  ipw_owner_safe_sum(c(retained_pipeline_entries, linear), cap, "pipeline plus compact certificate headers")
  validation_scratch <- ipw_owner_safe_sum(c(ipw_owner_safe_product(4, M, label = "structural reference groups and checked nuisance set"),
    ipw_owner_safe_product(3, K, label = "checked fitted views and nuisance set"), 12), cap,
    "simultaneous structural validation scratch")
  # One sealed canonical descriptor per original row (fields, pair descriptor,
  # environment and binding: 32 slots) plus the accepted P1 registry headers.
  registry_entries <- ipw_owner_safe_sum(c(ipw_owner_safe_product(16, M + 1, label = "compact registry headers"), C), cap,
    "compact registry headers")
  descriptor_entries <- ipw_owner_safe_product(32, M + 1, cap, "canonical descriptor records")
  # Operator nuisance view and fitted model (p x p factors, p x m partials,
  # coefficients, cohort vectors and headers), from the standalone formulas.
  nuis_pre <- preflight_ipw_nuisance_view(N, p, 2L, ctrl)
  link_pre <- preflight_ipw_link_view(N, p, 2L, ctrl)
  fitted_allowance <- ipw_owner_safe_sum(c(ipw_owner_safe_product(4, p * p, label = "fitted nuisance factor matrices"),
    ipw_owner_safe_product(8, p + 2, label = "fitted nuisance coefficient and partial vectors"), 512), cap,
    "fitted nuisance model allowance")
  nuisance_entries <- ipw_owner_safe_sum(c(nuis_pre$total_retained_entries, fitted_allowance), cap,
    "operator nuisance view and fitted model")
  # Retained basis: one class, hence at most path_bound forest rows; each keeps
  # its sealed link view and its 30-field row descriptor.
  path_bound <- min(M, T - 1, b$max_class_dates - 1)
  if (path_bound < 1) {
    cdid_abort("Compact basis needs at least two calendar dates in one class.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  L_bound <- path_bound
  E_bound <- min(b$max_E_entries, ipw_owner_safe_product(M, path_bound, label = "compact E entry upper bound"))
  basis_entries <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(L_bound, link_pre$total_retained_entries, label = "retained basis link views"),
    ipw_owner_safe_product(32, L_bound, label = "retained basis row descriptors")), cap,
    "retained basis link views and row descriptors")
  certificate_retained <- ipw_owner_safe_sum(c(5 * E_bound, linear, 4 * M * C), cap,
    "certificate E, two dropped-record expansion vectors, full-row maps and classes")
  live_registry <- ipw_owner_safe_sum(c(registry_entries, descriptor_entries, nuisance_entries, basis_entries), cap,
    "compact registry: headers, descriptors, nuisance model and retained basis")
  ipw_owner_safe_sum(c(live_registry, certificate_retained), ctrl$max_output_entries,
    "compact registry, certificate and reduced-link output bound")
  # Transient construction scratch: the larger constructor's scratch (views are
  # built one at a time), the certificate's sixteen E slots, the date-fact and
  # mask scratch of one chunk, and the validation scratch.
  construction_peak <- max(nuis_pre$concurrent_entries - nuis_pre$total_retained_entries,
    link_pre$concurrent_entries - link_pre$total_retained_entries)
  mask_scratch <- 64 * min(N, ctrl$max_chunk_rows) + 32 * C + 256
  additional_entries <- ipw_owner_safe_sum(c(live_registry, linear, 4 * M * C, 16 * E_bound, mask_scratch,
    validation_scratch, construction_peak), cap,
    "compact certificate construction, retained basis, nuisance model and expansion records")
  ipw_owner_safe_sum(c(retained_pipeline_entries, additional_entries), cap,
    "whole pipeline plus compact basis before allocation")
  # Work: complete-date facts (one chunked read of each catalog date's posting
  # slice), two binary pair-descriptor searches per row, the design/forest/
  # expansion work of the view route (no streamed mask), the per-row descriptor
  # integrity checks, the nuisance construction, and at most L_bound link
  # constructions with their per-stream envelopes (support scan and
  # fact-decided verification, charged as on the registry loop).
  stream_cost <- 128 * N * (1 + ceiling(log2(2 + 1))) + 128 * ceiling(N / step) * (1 + ceiling(log2(N + 1)))
  date_fact_work <- ipw_owner_safe_product(D, 16 * N + 64 * ceiling(N / step) + 64, label = "complete-date posting facts")
  descriptor_work <- ipw_owner_safe_product(M, 64 * (1 + ceiling(log2(T + 1))) + 512 * (C + 2), label = "canonical descriptor checks")
  structure_work <- ipw_owner_safe_sum(c(256 * M * (T + (K + 1) * (C + 1) + ceiling(log2(M + 1))), 64 * E_bound),
    ctrl$max_work_ops, "compact design/forest/expansion work")
  construction_work <- ipw_owner_safe_sum(c(nuis_pre$work_ops_bound,
    ipw_owner_safe_product(L_bound, link_pre$work_ops_bound, label = "retained basis link constructions"),
    ipw_owner_safe_product(2 * L_bound + 1, stream_cost, label = "retained basis stream envelopes")),
    ctrl$max_work_ops, "compact basis construction work")
  work_ops_bound <- ipw_owner_safe_sum(c(date_fact_work, descriptor_work, structure_work, construction_work),
    ctrl$max_work_ops, "compact structural facts/design/forest/expansion/construction work")
  ipw_owner_safe_sum(c(external_work_ops, work_ops_bound), ctrl$max_work_ops,
    "whole pipeline plus compact basis work before allocation")
  list(route = "compact_canonical_basis",
    E_entries_upper = E_bound, fitted_nuisance_entries = fitted_allowance, nuisance_view_entries = nuis_pre$total_retained_entries,
    original_link_view_entries = 0, retained_basis_bound = L_bound, retained_basis_link_view_entries = basis_entries,
    canonical_descriptor_entries = descriptor_entries, p1_registry_entries = registry_entries,
    mask_split_entries = 0, mask_certificate_entries = 0, mask_streams_charged = 0,
    date_facts_checked = D,
    mask_evidence = "complete-date posting facts: every catalog date posts exactly the whole A roster, so every row mask equals the operator nuisance mask without a streamed comparison",
    certificate_retained_entries = certificate_retained, additional_entries = additional_entries,
    work_ops_bound = work_ops_bound)
}

# The same structural certificate from the validated catalog's canonical row
# descriptors and the operator's fitted nuisance, with no link view. Per row the
# actual facts are re-established here from the sealed objects (never read from
# a label or a stored flag): the fitted model is validated once with its
# nuisance view against this frame; the row's cohort and exact ordered controls
# must equal the nuisance view's; the row's pair descriptor is built by the
# owner's own constructor and its mask equality with the nuisance mask must be
# decided by the complete-date posting fact (`ipw_owner_pair_masks_equal_fact`
# on this panel's sealed postings), otherwise the certificate is refused (the
# caller never retries another route). Each row then receives a sealed
# canonical descriptor bound to the metadata, frame, panel, nuisance and fitted
# references; the shared design checks and class/forest/E assembly run on the
# descriptors exactly as on the view route.
certify_ipw_structural_classes_from_descriptors <- function(metadata_obj, frame_owner, nuisance_views, fitted_views,
                                                            link_to_nuisance_idx, design, design_link_index, bounds,
                                                            owner_control, retained_pipeline_entries = 0,
                                                            descriptor_src = NULL) {
  b <- validate_ipw_certificate_bounds(bounds)
  own_ctrl <- validate_ipw_owner_control(owner_control)
  retained_pipeline_entries <- ipw_owner_int_bound(retained_pipeline_entries, "retained_pipeline_entries", min_val = 0)
  validate_owner_comparison_catalog(metadata_obj, frame_owner = frame_owner)
  validate_ipw_frame_owner(frame_owner)
  if (!identical(frame_owner$ref, metadata_obj$frame_ref) || !identical(frame_owner$panel_ref, metadata_obj$panel_ref)) {
    cdid_abort("Catalog metadata and frame owner are different sealed snapshots.",
      c("cdid_ipw_context_mismatch_error", "cdid_ipw_certificate_error"))
  }
  panel <- get("panel_owner", envir = frame_owner$panel_ref, inherits = FALSE)
  validate_ipw_frame_owner(frame_owner, panel_owner = panel)
  rows <- metadata_obj$catalog$rows
  n_links <- nrow(rows)
  idx <- link_to_nuisance_idx
  if (!is.list(nuisance_views) || is.object(nuisance_views) || length(nuisance_views) < 1L || !is.list(fitted_views) ||
      length(fitted_views) != length(nuisance_views) || !is.numeric(idx) || length(idx) != n_links || anyNA(idx) ||
      any(idx < 1 | idx > length(fitted_views) | idx != floor(idx)) || n_links < 1L) {
    cdid_abort(
      "Compact structural certification needs the validated catalog, one fitted nuisance per nuisance view and a row-to-view map.",
      c("cdid_ipw_input_error", "cdid_ipw_certificate_error")
    )
  }
  if (n_links > b$max_links) {
    cdid_abort("Number of links exceeds `max_links`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  n_a_units <- as.integer(frame_owner$n_a_units)
  if (as.double(n_links) * as.double(n_a_units) > b$max_mask_entries) {
    cdid_abort("Mask verification work exceeds `max_mask_entries`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  K_views <- length(fitted_views)
  cap <- min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries)
  ipw_owner_safe_sum(c(retained_pipeline_entries, ipw_owner_safe_product(3, K_views, label = "checked fitted views and nuisance set"),
    ipw_owner_safe_product(4, n_links, label = "structural reference groups and canonical descriptors"), 4), cap,
    "pipeline envelope plus compact structural checked-view scratch")
  step <- min(as.double(n_a_units), as.double(own_ctrl$max_chunk_rows))
  # Fresh call-local date facts of this sealed panel only; no caller-populated memo.
  memo <- list2env(list(panel_owner = panel), parent = emptyenv())
  mask_bad <- function(msg) cdid_abort(msg, c("cdid_ipw_shared_mask_error", "cdid_ipw_certificate_error"))
  for (k in seq_len(K_views)) {
    fv <- fitted_views[[k]]
    nv <- nuisance_views[[k]]
    validate_ipw_fitted_nuisance_view(fv, nuisance_view = nv, frame_owner = frame_owner)
    if (!identical(nv$frame_ref, frame_owner$ref) || !identical(nv$panel_ref, frame_owner$panel_ref) ||
        !identical(fv$nuisance_view_ref, nv$ref) || !identical(fv$frame_ref, frame_owner$ref) ||
        !identical(fv$panel_ref, frame_owner$panel_ref)) {
      cdid_abort(paste0("Nuisance view ", k, " or its fitted model is not bound to this sealed frame and panel."),
        c("cdid_ipw_context_mismatch_error", "cdid_ipw_certificate_error"))
    }
    if (!isTRUE(ipw_owner_pair_complete(panel, nv$defining_pair, step, memo))) {
      mask_bad(paste0("Nuisance view ", k, " defining pair is not complete on the sealed postings; ",
        "the compact certificate is not licensed."))
    }
  }
  desc_src <- if (is.null(descriptor_src)) {
    pipeline_metadata_descriptor_source(metadata_obj, own_ctrl, retained_pipeline_entries, 0, n_links)
  } else {
    if (!is.environment(descriptor_src) || !identical(descriptor_src$rows, rows) ||
        !identical(descriptor_src$map, metadata_obj$row_mapping)) {
      cdid_abort("Descriptor source does not alias this validated catalog.", c("cdid_ipw_input_error", "cdid_ipw_certificate_error"))
    }
    descriptor_src
  }
  link_g <- link_s <- link_t <- numeric(n_links)
  link_fk <- vector("list", n_links)
  link_K <- vector("list", n_links)
  pairs <- vector("list", n_links)
  records <- vector("list", n_links)
  for (ell in seq_len(n_links)) {
    pipeline_metadata_row_spec_preflight(metadata_obj, ell, own_ctrl)
    pipeline_metadata_row_spec_preflight(metadata_obj, ell, metadata_obj$control$owner_control)
    spec <- pipeline_metadata_row_spec_fields(metadata_obj, ell, desc_src)
    k <- as.integer(idx[ell])
    nv <- nuisance_views[[k]]
    fv <- fitted_views[[k]]
    g <- as.double(spec$g)
    K_row <- as.double(spec$K)
    if (!identical(g, nv$cohorts$g) || !identical(K_row, nv$cohorts$K)) {
      mask_bad(paste0("Catalog row ", ell, " declares a cohort or ordered control pool differing from its nuisance view."))
    }
    pair <- ipw_owner_make_pair_descriptor(panel, spec$s, spec$t)
    if (!isTRUE(ipw_owner_pair_masks_equal_fact(panel, nv$defining_pair, pair, step, memo))) {
      mask_bad(paste0("Catalog row ", ell, " mask equality with its nuisance view is not decided by complete-date ",
        "posting facts; the compact certificate is not licensed."))
    }
    ipw_owner_check_interrupt()
    link_g[ell] <- g
    link_s[ell] <- pair$s
    link_t[ell] <- pair$t
    link_fk[ell] <- list(spec$functional_key)
    link_K[[ell]] <- K_row
    pairs[[ell]] <- pair
    records[[ell]] <- ipw_owner_seal(structure(list(
      link_id = paste0("link:r", spec$row_id),
      catalog_row = as.integer(ell),
      operator_id = as.integer(spec$operator_id),
      g = g, K = K_row, s = pair$s, t = pair$t, pair = pair,
      functional_key = spec$functional_key,
      frame_ref = frame_owner$ref, panel_ref = frame_owner$panel_ref, metadata_ref = metadata_obj$ref,
      nuisance_view_ref = nv$ref, fitted_nuisance_ref = fv$ref, glm_control = fv$glm_control,
      mask_evidence = "complete_dates",
      link_view_constructed = FALSE
    ), class = "cdid_ipw_canonical_link_descriptor"), "link_descriptor")
  }
  checks <- ipw_structural_design_checks(link_g, link_s, link_t, link_K, link_fk, design, design_link_index, n_links)
  nu_gid <- as.integer(idx)
  nu_refs <- lapply(nuisance_views, `[[`, "ref")
  # Equality of row ell's mask with source representative rep_ell's nuisance
  # mask, decided only by the same complete-date posting facts.
  same_mask_fn <- function(ell, rep_ell) {
    isTRUE(ipw_owner_pair_masks_equal_fact(panel, nuisance_views[[nu_gid[rep_ell]]]$defining_pair, pairs[[ell]], step, memo))
  }
  certificate <- ipw_structural_classes_assemble(n_links, link_s, link_t, checks$fk, checks$vs, checks$vt, checks$dli,
    nu_gid, nu_refs,
    link_cohort_g = lapply(seq_len(n_links), function(ell) fitted_views[[nu_gid[ell]]]$cohorts$g),
    link_cohort_K = lapply(seq_len(n_links), function(ell) fitted_views[[nu_gid[ell]]]$cohorts$K),
    link_glm = lapply(seq_len(n_links), function(ell) fitted_views[[nu_gid[ell]]]$glm_control),
    same_mask_fn = same_mask_fn, b = b, link_kind = "canonical_descriptor_declaration",
    object_refs = lapply(records, `[[`, "ref"), context_ref = frame_owner$ref)
  list(certificate = certificate, descriptors = records, pairs = pairs,
    link_ids = vapply(records, `[[`, character(1L), "link_id"),
    n_date_facts = sum(startsWith(ls(memo, all.names = TRUE), "complete:")))
}

# Compact registry of the route: one nuisance view and one fitted model for the
# operator (constructed from the operator's first catalog row exactly as the
# registry loop does, so the nuisance identity, defining pair and fit are the
# same), the certificate from canonical descriptors before any link view, then
# the sealed link views of the L retained rows only, through the unchanged
# constructor path of the registry loop (row preflights, live-envelope
# admission, descriptor source, construction context, support reuse and the
# fact-decided verification). Its result carries the registry fields the
# pipeline consumes; dropped rows hold no link view, only their sealed
# canonical descriptor (see `ipw_pipeline_compact_reduced_link_record`).
build_pipeline_compact_basis_registry <- function(
  metadata_obj,
  frame_owner,
  design,
  design_link_index,
  facts,
  owner_control,
  fitted_control,
  point_control,
  glm_control,
  class_bounds,
  retained_pipeline_entries = 0,
  external_work_ops = 0,
  progress = NULL
) {
  validate_owner_comparison_catalog(metadata_obj, frame_owner = frame_owner)
  own_ctrl <- validate_ipw_owner_control(owner_control)
  fit_ctrl <- validate_ipw_fitted_view_control(fitted_control)
  pt_ctrl <- validate_ipw_point_view_control(point_control)
  glm_ctrl <- validate_ipw_glm_control(glm_control)
  cat_rows <- metadata_obj$catalog$rows
  M <- nrow(cat_rows)
  if (M < 1L) {
    cdid_ipw_pipeline_abort(
      "Comparison catalog contains zero qualified links; cannot construct IPW views.",
      c("cdid_ipw_pipeline_empty_catalog_error", "cdid_gmm_target_estimability_error")
    )
  }
  f_owner <- frame_owner
  validate_ipw_frame_owner(f_owner)
  panel <- get("panel_owner", envir = metadata_obj$panel_ref, inherits = FALSE)
  validate_ipw_frame_owner(f_owner, panel_owner = panel)
  if (!is.list(facts) || !identical(facts$frame_ref, f_owner$ref) || !identical(facts$panel_ref, f_owner$panel_ref) ||
      !identical(facts$metadata_ref, metadata_obj$ref) || !identical(facts$panel_owner$ref, panel$ref) ||
      !identical(as.integer(facts$n_rows), as.integer(M)) || !gmm_scalar_numeric(facts$first_row, positive = TRUE)) {
    cdid_abort("Compact basis facts do not describe this validated catalog and these sealed owners.",
      c("cdid_ipw_input_error", "cdid_ipw_certificate_error"))
  }
  retained_pipeline_entries <- ipw_owner_int_bound(retained_pipeline_entries, "retained_pipeline_entries", min_val = 0)
  external_work_ops <- ipw_owner_int_bound(external_work_ops, "external_work_ops", min_val = 0)
  cap <- min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries)

  # Everything this route allocates is reserved here, simultaneously, before
  # the first descriptor, view, fit or certificate buffer exists.
  pre <- ipw_pipeline_compact_basis_preflight(facts, design, class_bounds, own_ctrl,
    retained_pipeline_entries, external_work_ops)
  live <- ipw_owner_safe_sum(c(retained_pipeline_entries, pre$additional_entries), cap,
    "pipeline envelope plus compact basis reservation")
  work <- external_work_ops
  stream_n <- as.double(f_owner$n_a_units)
  stream_step <- max(1, min(stream_n, own_ctrl$max_chunk_rows))
  stream_cost <- function(m) 128 * stream_n * (1 + ceiling(log2(m + 1))) +
    128 * ceiling(stream_n / stream_step) * (1 + ceiling(log2(stream_n + 1)))
  scan_work <- 0
  charge_scans <- function(n_streams, m) {
    scan_work <<- ipw_owner_safe_sum(c(scan_work, ipw_owner_safe_product(n_streams, stream_cost(m))))
    ipw_owner_safe_sum(c(external_work_ops, scan_work), own_ctrl$max_work_ops,
      "pipeline work plus compact registry support scans (fact-decided verifications stay charged)")
  }

  construction_ctx <- ipw_owner_construction_context(f_owner)
  ipw_owner_context_check_frame(construction_ctx, f_owner)
  descriptor_src <- pipeline_metadata_descriptor_source(metadata_obj, own_ctrl, live, external_work_ops, M, memo_entries = 64)
  descriptor_charge <- descriptor_src$accounting
  live <- ipw_owner_safe_sum(c(live, descriptor_charge$retained_entries), cap, "compact registry plus descriptor context")
  work <- ipw_owner_safe_sum(c(work, descriptor_charge$work_ops_bound), own_ctrl$max_work_ops,
    "compact registry plus descriptor preparation work")

  # Operator nuisance view: the first catalog row of the operator, as the registry loop.
  first_row <- as.integer(facts$first_row)
  pipeline_metadata_row_spec_preflight(metadata_obj, first_row, own_ctrl)
  pipeline_metadata_row_spec_preflight(metadata_obj, first_row, metadata_obj$control$owner_control)
  spec1 <- pipeline_metadata_row_spec_fields(metadata_obj, first_row, descriptor_src)
  op_id <- as.integer(spec1$operator_id)
  if (!identical(op_id, as.integer(facts$operator_id)) || !identical(as.double(spec1$g), facts$g) ||
      !identical(as.double(spec1$K), facts$K) || !isTRUE(all(cat_rows$operator_id == op_id))) {
    cdid_abort("Compact basis facts disagree with the catalog's operator rows.", c("cdid_ipw_input_error", "cdid_ipw_certificate_error"))
  }
  row_cohorts1 <- 1 + metadata_obj$catalog$pools$size[cat_rows$pool_id[first_row]]
  nuis_pre <- preflight_ipw_nuisance_view(f_owner$n_a_units, f_owner$n_covariates, row_cohorts1, own_ctrl)
  ipw_owner_safe_sum(c(pre$p1_registry_entries, pre$canonical_descriptor_entries, max(retained_pipeline_entries, nuis_pre$total_retained_entries)),
    own_ctrl$max_retained_entries, "compact registry plus nuisance constructor retention")
  ipw_owner_safe_sum(c(pre$p1_registry_entries, pre$canonical_descriptor_entries, max(retained_pipeline_entries, nuis_pre$concurrent_entries),
    descriptor_charge$scratch_entries), own_ctrl$max_matrix_entries, "compact registry plus nuisance constructor and descriptor scratch envelope")
  charge_scans(1, row_cohorts1)
  nv <- construct_ipw_nuisance_view(
    frame_owner = f_owner, g = spec1$g, K = spec1$K, s = spec1$s, t = spec1$t, control = own_ctrl,
    nuisance_id = paste0("nuisance:op", op_id, ":r", spec1$row_id), require_support = TRUE,
    construction_context = construction_ctx
  )
  # The actual fitted-model preflight composed with the live envelope before the fit allocates.
  fit_pre <- preflight_ipw_fitted_nuisance_view(nv, glm_ctrl, fit_ctrl, n_planned_fits = 1L)
  ipw_owner_safe_sum(c(live, fit_pre$peak_concurrent_entries), cap, "compact registry live envelope plus nuisance fit workspace")
  work <- ipw_owner_safe_sum(c(work, nv$preflight$work_ops_bound, fit_pre$total_work_ops_bound,
    fit_pre$character_metadata$work_ops_bound, 128), own_ctrl$max_work_ops,
    "compact registry plus nuisance construction and fitting work")
  fitted_views <- fit_ipw_nuisance_view_set(nuisance_views = list(nv), glm_control = glm_ctrl, control = fit_ctrl)
  key <- paste0("op", op_id, ":mask1")
  names(fitted_views) <- key
  nuisance_views <- list(nv)
  names(nuisance_views) <- key
  fit_actual <- fitted_views[[1L]]$accounting$total_retained_entries
  if (gmm_scalar_numeric(fit_actual, nonneg = TRUE) && fit_actual > pre$fitted_nuisance_entries) {
    live <- ipw_owner_safe_sum(c(live, fit_actual - pre$fitted_nuisance_entries), cap,
      "compact registry live envelope plus fitted nuisance beyond its allowance")
  }

  # The certificate from canonical descriptors: no link view exists yet.
  cert_res <- certify_ipw_structural_classes_from_descriptors(metadata_obj, f_owner, nuisance_views, fitted_views,
    rep(1L, M), design, design_link_index, class_bounds, own_ctrl,
    retained_pipeline_entries = live, descriptor_src = descriptor_src)
  certificate <- cert_res$certificate
  kept_idx <- certificate$retained_index
  dropped_idx <- certificate$dropped_index
  L <- length(kept_idx)
  if (L > pre$retained_basis_bound) {
    cdid_abort("Certified retained basis exceeds its preflight bound.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }

  # Heavy link views of the retained rows only, through the registry loop's path.
  link_views <- vector("list", M)
  row_specs <- vector("list", M)
  dim_memo <- NULL
  for (h in seq_along(kept_idx)) {
    m <- kept_idx[h]
    row_cohorts <- 1 + metadata_obj$catalog$pools$size[cat_rows$pool_id[m]]
    dim_inputs <- list(f_owner$n_a_units, f_owner$n_covariates, row_cohorts, own_ctrl)
    if (is.null(dim_memo) || !identical(dim_memo$inputs, dim_inputs)) {
      row_link_preflight <- preflight_ipw_link_view(f_owner$n_a_units, f_owner$n_covariates, row_cohorts, own_ctrl)
      dim_memo <- list(inputs = dim_inputs, peak = row_link_preflight$concurrent_entries,
        retained = row_link_preflight$total_retained_entries)
      row_link_preflight <- NULL
    }
    ipw_owner_safe_sum(c(pre$p1_registry_entries, pre$canonical_descriptor_entries, max(retained_pipeline_entries, dim_memo$retained)),
      own_ctrl$max_retained_entries, "compact registry plus existing row constructor retention")
    ipw_owner_safe_sum(c(pre$p1_registry_entries, pre$canonical_descriptor_entries, max(retained_pipeline_entries, dim_memo$peak),
      descriptor_charge$scratch_entries), own_ctrl$max_matrix_entries,
      "compact registry plus row constructor and descriptor scratch envelope")
    spec <- pipeline_metadata_row_spec_fields(metadata_obj, m, descriptor_src)
    # Link support scan plus its verification against the nuisance view; the
    # stream charge stays charged when the posting fact decides it.
    charge_scans(2, row_cohorts)
    row_res <- pipeline_metadata_row_views_from_spec(spec, f_owner, own_ctrl,
      nuisance_view = nv, verify_nuisance_mask = TRUE, construction_context = construction_ctx)
    lv <- row_res$link_view
    if (!inherits(lv, "cdid_ipw_link_view") || !isTRUE(lv$nuisance_mask_verified) || !identical(lv$nuisance_ref, nv$ref) ||
        !identical(lv$link_id, cert_res$link_ids[m])) {
      cdid_abort(paste0("Retained row ", m, " link view is not a verified sharer of the operator nuisance view."),
        c("cdid_ipw_shared_mask_error", "cdid_ipw_certificate_error"))
    }
    row_specs[[m]] <- row_res$row_spec
    link_views[[m]] <- lv
    if (is.function(progress)) progress(h, L, 1L)
  }
  names(link_views) <- cert_res$link_ids
  n_support_reused <- construction_ctx$support_reuses
  construction_ctx <- NULL
  descriptor_src <- dim_memo <- NULL

  mask_facts <- ipw_seal(list(
    kind = "cdid_ipw_compact_basis_mask_facts",
    frame_ref = f_owner$ref, panel_ref = f_owner$panel_ref, metadata_ref = metadata_obj$ref,
    n_rows = as.double(M), row_operator = as.integer(cat_rows$operator_id), row_view = rep(1L, M),
    view_operator = op_id, view_g = as.double(nv$cohorts$g), view_K = list(as.double(nv$cohorts$K)),
    view_nuisance_refs = list(nv$ref), view_fit_refs = list(fitted_views[[1L]]$ref),
    date_indices_complete = as.integer(facts$date_indices), n_date_facts_reread = cert_res$n_date_facts,
    claim = paste0("Every catalog date posts exactly the whole A roster on the sealed panel, read chunk by chunk from the ",
      "postings in this call; hence every original row's risk-set pair-observation mask equals the operator nuisance ",
      "view's mask by the owner's posting fact (`ipw_owner_pair_masks_equal_fact`), the same fact the link constructor ",
      "records as a completed verification. No mask was streamed and no label or flag was trusted."),
    scope = "invocation-local evidence of this compact registry only; no statistical equivalence, sample or estimate is inferred"),
    "mask_facts")

  registry <- list(
    nuisance_views = nuisance_views,
    fitted_views = fitted_views,
    nuisance_keys = key,
    operator_ids = op_id,
    row_specs = row_specs,
    link_views = link_views,
    point_links = NULL,
    link_to_nuisance_idx = rep(1L, M),
    link_nuisance_keys = rep(key, M),
    mask_splits = list(),
    mask_certificate = NULL,
    compact_basis = list(certificate = certificate, descriptors = cert_res$descriptors, mask_facts = mask_facts,
      retained_index = kept_idx, dropped_index = dropped_idx, route = "compact_canonical_basis"),
    construction_work = list(
      n_rows = M,
      n_distinct_nuisance_views = 1L,
      candidate_comparisons = 0,
      fingerprint_rejections = 0,
      exact_candidate_attempts = 0,
      exact_mask_mismatches = 0,
      fingerprint_passes = 0,
      link_view_constructions = L,
      heavy_link_views_constructed = L,
      canonical_descriptor_rows = M,
      dropped_rows_without_link_view = length(dropped_idx),
      endpoints_read = 0L,
      endpoints_read_scope = "mask_screening_only",
      work_ops_charged = ipw_owner_safe_sum(c(scan_work, descriptor_charge$work_ops_bound)),
      registry_entries_reserved = pre$p1_registry_entries + pre$canonical_descriptor_entries,
      split_records_bound = 0,
      link_support_scans_reused = n_support_reused,
      complete_pair_fingerprints_reused = 0,
      route = "compact_canonical_basis",
      scope = paste0("Invocation-local counts of the compact canonical basis route. The structural certificate was ",
        "determined from the validated catalog's canonical descriptors, the operator's single fitted nuisance and ",
        "complete-date posting facts before any link view existed; link views were constructed for the ",
        L, " retained basis rows only, each through the unchanged constructor path (support reuse and fact-decided ",
        "verification), and every dropped row keeps a sealed canonical descriptor. work_ops_charged includes descriptor ",
        "preparation and the per-stream envelope charged before each support scan/verification.")
    )
  )
  list(registry = registry, certificate = certificate, preflight = pre, live_entries = live, work_ops = work)
}

# Outcome-free identity of a dropped original row on the compact route: its
# sealed canonical descriptor (bound to the metadata, frame, panel, nuisance and
# fitted references and to the certificate's object reference for that row),
# the fitted nuisance it shares and its certified expansion row. No link view
# was constructed for it; no estimate, influence function or fitted object exists.
ipw_pipeline_compact_reduced_link_record <- function(registry, m, certificate) {
  desc <- registry$compact_basis$descriptors[[m]]
  fv <- registry$fitted_views[[registry$link_to_nuisance_idx[m]]]
  if (!is.list(desc) || !inherits(desc, "cdid_ipw_canonical_link_descriptor") ||
      !ipw_owner_seal_ok(desc, "link_descriptor") || !identical(desc$catalog_row, as.integer(m)) ||
      !identical(desc$fitted_nuisance_ref, fv$ref) || !identical(desc$nuisance_view_ref, fv$nuisance_view_ref) ||
      !identical(certificate$object_refs[[m]], desc$ref)) {
    cdid_abort(paste0("Dropped row ", m, " has no sealed canonical descriptor bound to this certificate and fitted nuisance."),
      c("cdid_ipw_owner_integrity_error", "cdid_ipw_certificate_error"))
  }
  n_e <- length(certificate$E$i)
  first <- ipw_owner_lower_bound(certificate$E$i, 1, n_e, m)
  last <- ipw_owner_upper_bound(certificate$E$i, first, n_e, m)
  hit <- if (last >= first) seq.int(first, last) else integer(0L)
  structure(list(
    link_id = desc$link_id,
    catalog_row = as.integer(m),
    link_view_ref = NULL,
    descriptor_ref = desc$ref,
    nuisance_view_ref = fv$nuisance_view_ref,
    fitted_nuisance_ref = fv$ref,
    g = desc$g, K = desc$K, s = desc$s, t = desc$t, functional_key = desc$functional_key,
    expansion = list(retained_rows = certificate$E$j[hit], signs = certificate$E$x[hit]),
    point_fit_performed = FALSE,
    link_view_constructed = FALSE,
    outcomes_read = FALSE,
    status = "structurally_reduced_outcome_free_compact",
    scope = paste0("Original comparison identified by its sealed canonical catalog descriptor and recovered through the ",
      "certified expansion E (Prop. 5/Cor. 6); no link view, fit, replay or outcome read.")
  ), class = "cdid_ipw_reduced_link_record")
}

# ------------------------------------------------------------------------------
# 3. Complete Numerical Transfer Composition (Root + Link + Covariance + GMM)
# ------------------------------------------------------------------------------
# Implements Theorem 12 (`thm:gmm`) in `docs/theory/cp02_ipw_gmm_numerical_propagation.tex`:
# Composes:
#   (a) Certified nuisance root enclosures (`bind_fitted_nuisance_to_problem_view`)
#   (b) Certified link enclosures (`bind_point_link_to_enclosure_view`)
#   (c) Certified moment covariance transfer (`certify_ipw_covariance_transfer_view`)
#   (d) Outward interval changing-covariance/loading GMM transfer bounds for
#       target estimates (`total_target_error_cert`), target loadings
#       (`total_loading_error_cert`), target covariance (`total_cov_error_cert`),
#       and target standard errors (`total_se_error_cert`).

certify_chained_ipw_numerical_transfer_view <- function(
  point_links,
  design,
  mhat,
  V,
  gmm_fit,
  weighting = c("optimal", "identity"),
  root_witness = NULL,
  certificate_control,
  problem_control,
  link_control,
  covariance_control,
  point_control,
  reader_control,
  gmm_accuracy,
  gmm_bounds,
  retained_pipeline_entries
) {
  wt_mode <- match.arg(weighting)
  cert_ctrl <- validate_logit_certificate_control(certificate_control)
  prob_ctrl <- validate_ipw_nuisance_problem_control(problem_control)
  lnk_ctrl <- validate_ipw_link_enclosure_control(link_control)
  cov_ctrl <- validate_ipw_covariance_transfer_control(covariance_control)
  pt_ctrl <- validate_ipw_point_view_control(point_control)
  rd_ctrl <- validate_ipw_owner_control(reader_control)
  acc <- validate_gmm_accuracy_control(gmm_accuracy, wt_mode)
  bnds <- validate_gmm_bounds_control(gmm_bounds)
  if (is.null(acc$max_loading_error_abs)) {
    cdid_ipw_pipeline_abort("Complete transfer requires an explicit loading-error budget.", "cdid_ipw_pipeline_control_error")
  }

  L <- length(point_links)
  if (L < 1L || design$incidence$n_links != L || length(mhat) != L || nrow(V) != L || ncol(V) != L) {
    cdid_ipw_pipeline_abort(
      "Dimension mismatch among `point_links`, `design`, `mhat`, and `V` in numerical transfer.",
      "cdid_ipw_pipeline_input_error"
    )
  }

  # 1. Certify each distinct fitted nuisance view once via Krawczyk contraction
  nuisance_problems <- list()
  nuisance_refs <- list()
  nuisance_ids <- character(0L)
  link_enclosures <- vector("list", L)

  for (l in seq_len(L)) {
    pl <- point_links[[l]]
    validate_ipw_point_link_view(pl)
    f_ref <- pl$fitted_nuisance_ref
    found_idx <- NA_integer_
    for (k in seq_along(nuisance_refs)) {
      if (identical(f_ref, nuisance_refs[[k]])) {
        found_idx <- k
        break
      }
    }
    if (is.na(found_idx)) {
      fv <- get("fitted_nuisance_view", envir = f_ref, inherits = FALSE)
      witness <- root_witness(fv)
      if (!is.list(witness) || !identical(sort(names(witness)), c("F", "r"))) {
        cdid_ipw_pipeline_abort("Each caller root witness must contain exactly F and r.", "cdid_ipw_pipeline_control_error")
      }
      prob_k <- bind_fitted_nuisance_to_problem_view(
        fitted_view = fv,
        F = witness$F,
        r = witness$r,
        certificate_control = cert_ctrl,
        problem_control = prob_ctrl,
        reader_control = rd_ctrl
      )
      if (!isTRUE(prob_k$certified)) {
        # Unchanged class and message; the already computed certificate scalars,
        # fit diagnostics and risk counts of this one nuisance are retained
        # (section 0b), within the accepted per-nuisance reservation.
        cdid_ipw_pipeline_abort(
          paste0("Nuisance root certificate is INCONCLUSIVE for `", fv$nuisance_id, "`."),
          c("cdid_ipw_pipeline_transfer_inconclusive", "cdid_logit_certificate_inconclusive"),
          diagnostics = list(nuisance_id = fv$nuisance_id, reasons = prob_k$inconclusive_reasons,
            refusal_state = ipw_pipeline_scalar_state(ipw_pipeline_nuisance_refusal_state(prob_k, fv, witness, pl$link_id,
              length(nuisance_problems), L)))
        )
      }
      found_idx <- length(nuisance_problems) + 1L
      nuisance_problems[[found_idx]] <- prob_k
      nuisance_refs[[found_idx]] <- f_ref
      nuisance_ids[found_idx] <- as.character(fv$nuisance_id)
    }

    enc_l <- bind_point_link_to_enclosure_view(
      problem = nuisance_problems[[found_idx]],
      point_link = pl,
      link_control = lnk_ctrl,
      reader_control = rd_ctrl,
      include_candidates = TRUE
    )
    if (!isTRUE(enc_l$certified)) {
      cdid_ipw_pipeline_abort(
        paste0("Link enclosure is INCONCLUSIVE for link `", pl$link_id, "`."),
        c("cdid_ipw_pipeline_transfer_inconclusive", "cdid_ipw_enclosure_inconclusive"),
        diagnostics = list(link_id = pl$link_id, reasons = enc_l$inconclusive_reasons)
      )
    }
    link_enclosures[[l]] <- enc_l
  }

  # 2. Stream covariance transfer certificate over all retained links
  frame_owner <- get("frame_owner", envir = point_links[[1L]]$frame_ref, inherits = FALSE)
  n_a_units <- as.integer(frame_owner$n_a_units)
  cand_spec <- make_ipw_point_candidate_rows_spec(
    point_links = point_links,
    control = pt_ctrl
  )
  cov_op <- ipw_covariance_operator_spec(
    purpose = "GMM",
    row_level = "unit",
    omega = rep(1.0, n_a_units),
    centered = TRUE,
    kappa = 1.0 / as.double(n_a_units)
  )
  frame_header <- ipw_cov_transfer_validate_frame_owner_header(frame_owner)
  links_header <- ipw_cov_transfer_validate_links_view_headers(link_enclosures,frame_header)
  cov_pre <- preflight_ipw_covariance_transfer_view(n=n_a_units,L=L,d=L,row_level="unit",
    n_clusters=n_a_units,p_vec=links_header$p_vec,m_vec=links_header$m_vec,
    certificate_controls=links_header$certificate_controls,control=cov_ctrl,reader_control=rd_ctrl,
    id_is_character=is.character(frame_owner$id),existing_4bd_retained_entries=links_header$existing_4bd_retained_entries)
  producer_pre <- preflight_ipw_point_if_view(point_links,pt_ctrl)
  combined_cov <- ipw_owner_safe_sum(c(retained_pipeline_entries,cov_pre$total_with_existing_4bd_entries,
    producer_pre$peak_concurrent_entries,gmm_fit$allocations$total_payload_entries,
    as.double(n_a_units)*ncol(gmm_fit$internal_gauge$Lambda)*32.0,as.double(L)*L*32.0),
    rd_ctrl$max_matrix_entries,"pipeline simultaneous point producer, interval reader and GMM consumer")
  cov_transfer <- certify_ipw_covariance_transfer_view(
    links = link_enclosures,
    candidate_rows = cand_spec,
    frame_owner = frame_owner,
    C = diag(1.0, nrow = L, ncol = L),
    op = cov_op,
    V_candidate = unname(V),
    control = cov_ctrl,
    reader_control = rd_ctrl
  )
  if (!(if (identical(wt_mode, "optimal")) isTRUE(cov_transfer$certified) else
        isTRUE(cov_transfer$enclosure_certified) && is.finite(cov_transfer$epsilon_V))) {
    cdid_ipw_pipeline_abort(
      "Covariance transfer certificate is INCONCLUSIVE.",
      c("cdid_ipw_pipeline_transfer_inconclusive", "cdid_ipw_covariance_transfer_inconclusive"),
      diagnostics = list(reasons = cov_transfer$inconclusive_reasons)
    )
  }

  # 3. Compose Theorem 12 (`thm:gmm`) outward interval bounds for changing-weight GMM
  eps_m_vec <- numeric(L)
  for (l in seq_len(L)) {
    sep_mu <- link_enclosures[[l]]$separated_errors$mu
    if (is.null(sep_mu) || !is.finite(sep_mu$total_error_bound)) {
      cdid_ipw_pipeline_abort(
        paste0("Link enclosure for `", point_links[[l]]$link_id, "` lacks finite `separated_errors$mu$total_error_bound`."),
        "cdid_ipw_pipeline_transfer_inconclusive"
      )
    }
    eps_m_vec[l] <- as.double(sep_mu$total_error_bound)
  }

  eps_m_sq_iv <- gmm_iv(0.0)
  mhat_sq_iv <- gmm_iv(0.0)
  for (l in seq_len(L)) {
    e_iv <- gmm_iv(eps_m_vec[l])
    eps_m_sq_iv <- gmm_iv_add(eps_m_sq_iv, gmm_iv_mul(e_iv, e_iv))
    m_iv <- gmm_iv(mhat[l])
    mhat_sq_iv <- gmm_iv_add(mhat_sq_iv, gmm_iv_mul(m_iv, m_iv))
  }
  eps_m_l2 <- gmm_iv_sqrt(eps_m_sq_iv)[2L]
  mhat_l2_upper <- gmm_iv_sqrt(mhat_sq_iv)[2L]

  eps_V <- as.double(cov_transfer$epsilon_V)
  nu_I <- as.double(cov_transfer$nu_total)
  hat_V_frob <- as.double(cov_transfer$v_candidate_norm_upper)
  omega_star_norm_upper <- gmm_iv_add(gmm_iv(hat_V_frob), gmm_iv(eps_V))[2L]
  sqrt_omega_star_iv <- gmm_iv_sqrt(gmm_iv(omega_star_norm_upper))

  grd <- gmm_fit$numerical_guard
  pb <- gmm_fit$propagation_bounds
  if (is.null(grd) || !isTRUE(grd$point_certified)) {
    cdid_ipw_pipeline_abort(
      "Fixed-input GMM numerical guard is missing or not point-certified.",
      "cdid_ipw_pipeline_transfer_inconclusive"
    )
  }

  Lambda_hat <- as.matrix(gmm_fit$internal_gauge$Lambda)
  J_est <- ncol(Lambda_hat)

  c_hat_norms <- as.double(pb$lambda_norm2_cert)
  e_r_norms <- as.double(pb$c_norm2_cert)
  delta_C_num <- as.double(pb$G_lambda)
  delta_theta_solve <- as.double(pb$G_c)

  if (identical(wt_mode, "optimal")) {
    l_lower <- as.double(pb$l_lambda_min_V_cert)
    u_upper <- as.double(pb$lambda_max_V_cert)
    l_minus_eps_V <- gmm_iv_sub(gmm_iv(l_lower), gmm_iv(eps_V))[1L]
    if (!is.finite(l_minus_eps_V) || l_minus_eps_V <= 0.0) {
      cdid_ipw_pipeline_abort(
        paste0(
          "Certified lower eigenvalue margin `l_lower - eps_V` (",
          format(l_minus_eps_V, digits = 6),
          ") is non-positive; optimal-weight GMM transfer is INCONCLUSIVE."
        ),
        c("cdid_ipw_pipeline_transfer_inconclusive", "cdid_gmm_pd_error"),
        diagnostics = list(l_lower = l_lower, eps_V = eps_V, l_minus_eps_V = l_minus_eps_V)
      )
    }
    sqrt_l_iv <- gmm_iv_sqrt(gmm_iv(l_lower))
    sqrt_l_m_eps_iv <- gmm_iv_sqrt(gmm_iv(l_minus_eps_V))
    sqrt_l_prod_iv <- gmm_iv_mul(sqrt_l_iv, sqrt_l_m_eps_iv)
    denom_loading_stat_iv <- gmm_iv_mul(sqrt_l_iv, gmm_iv(l_minus_eps_V))
    one_plus_eps_over_l_iv <- gmm_iv_add(gmm_iv(1.0), gmm_iv_div(gmm_iv(eps_V), gmm_iv(l_lower)))
    sqrt_one_plus_eps_over_l_iv <- gmm_iv_sqrt(one_plus_eps_over_l_iv)
    one_minus_eps_over_l_iv <- gmm_iv_div(gmm_iv(l_minus_eps_V), gmm_iv(l_lower))
    sqrt_one_minus_eps_over_l_iv <- gmm_iv_sqrt(one_minus_eps_over_l_iv)
    sqrt_Q_upper <- as.double(pb$sqrt_Q_cert_upper)
    s_tilde_upper <- as.double(pb$se_criterion_cert_upper)
    s_tilde_lower <- as.double(pb$se_criterion_cert_lower)
    delta_W_iv <- gmm_iv_div(gmm_iv(eps_V), gmm_iv_mul(gmm_iv(l_lower), gmm_iv(l_minus_eps_V)))
  } else {
    l_lower <- 1.0
    u_upper <- as.double(hat_V_frob)
    l_minus_eps_V <- 1.0
    sqrt_l_iv <- gmm_iv(1.0)
    sqrt_l_m_eps_iv <- gmm_iv(1.0)
    sqrt_l_prod_iv <- gmm_iv(1.0)
    denom_loading_stat_iv <- gmm_iv(1.0)
    sqrt_one_plus_eps_over_l_iv <- gmm_iv(1.0)
    sqrt_one_minus_eps_over_l_iv <- gmm_iv(1.0)
    sqrt_Q_upper <- 0.0
    s_tilde_upper <- c_hat_norms
    s_tilde_lower <- c_hat_norms
    delta_W_iv <- gmm_iv(0.0)
  }

  s_star_lower <- numeric(J_est)
  s_star_upper <- numeric(J_est)
  delta_C_stat <- numeric(J_est)
  total_loading_error_cert <- numeric(J_est)
  total_target_error_cert <- numeric(J_est)
  R_c_cert <- numeric(J_est)
  total_cov_diag_error_cert <- numeric(J_est)
  total_se_error_cert <- numeric(J_est)
  # Retained outward components of the Theorem 12 point and row-IF bounds
  # (diagnostic output only; the certified totals below are unchanged). Each
  # component is the upper endpoint of the corresponding already-formed interval
  # term; with J_est = 0 every vector is empty and no maximum is taken.
  point_err_solve <- numeric(J_est)
  point_err_moment <- numeric(J_est)
  point_err_weights <- numeric(J_est)
  rc_err_proj <- numeric(J_est)
  rc_err_rowif <- numeric(J_est)
  rc_err_loading <- numeric(J_est)
  rc_err_weights <- numeric(J_est)

  Sigma_sand <- as.matrix(gmm_fit$sandwich_cov)
  diag_Sigma <- diag(Sigma_sand)
  projection_sq <- rep(0.0, J_est)
  stored_rows <- gmm_fit$if_projection$materialized_target_if
  # D5: with zero estimable targets the covariance/PD/guard checks above still
  # apply; the target-specific projection/accumulation terms are vacuous (no
  # loading column exists), so nothing target-wise is replayed or certified.
  if (J_est > 0L) replay_ipw_point_if_view(point_links, control = pt_ctrl, callback_schema = TRUE,
    callback = function(batch) {
      idx <- batch$a_index
      if (is.null(idx)) idx <- match(batch$id, frame_owner$a_id)
      # Native directed dot products include all cancellation and underflow.
      residual <- gmm_iv_mat_residual(unname(batch$phi), unname(Lambda_hat),
        unname(stored_rows[idx,,drop=FALSE]), max_matrix_entries = bnds$max_matrix_entries)
      for (j in seq_len(J_est)) {
        norm_iv <- gmm_iv_norm(residual$lo[,j], residual$hi[,j],
                               max_matrix_entries = bnds$max_matrix_entries)
        projection_sq[j] <<- gmm_iv_add(gmm_iv(projection_sq[j]),
                                        gmm_iv_mul(gmm_iv(norm_iv[2L]), gmm_iv(norm_iv[2L])))[2L]
      }
    })
  projection_error <- vapply(projection_sq, function(z)
    gmm_iv_sqrt(gmm_iv_div(gmm_iv(z), gmm_iv(as.double(n_a_units))))[2L], numeric(1L))
  # Exact centered row covariance, evaluated with directed accumulation.
  stored_enclosure <- list(rows = list(lo = stored_rows, hi = stored_rows),
                            theta = list(lo = matrix(gmm_fit$targets$estimate[gmm_fit$targets$estimable],1L),
                                         hi = matrix(gmm_fit$targets$estimate[gmm_fit$targets$estimable],1L)))
  row_acc <- ipw_pipeline_inference_intervals(stored_enclosure, NULL)
  acc_sigma <- numeric(J_est)
  for (j in seq_len(J_est)) {
    si <- gmm_iv_sqrt(gmm_iv(max(0.0, Sigma_sand[j,j])))
    row_si <- gmm_iv_mul(row_acc$se[,j], gmm_iv_sqrt(gmm_iv(as.double(n_a_units))))
    acc_sigma[j] <- max(abs(gmm_iv_sub(si, row_si)))
  }
  row_cov_lo <- row_acc$vcov$lo * n_a_units
  row_cov_hi <- row_acc$vcov$hi * n_a_units
  for (k in seq_along(row_cov_lo)) {
    z <- gmm_iv_mul(c(row_acc$vcov$lo[k], row_acc$vcov$hi[k]), gmm_iv(as.double(n_a_units)))
    row_cov_lo[k] <- z[1L]; row_cov_hi[k] <- z[2L]
  }
  acc_cov_dist <- ipw_pipeline_interval_distance(Sigma_sand, row_cov_lo, row_cov_hi)
  acc_cov <- if (J_est > 0L) gmm_iv_norm(unname(acc_cov_dist), max_matrix_entries = bnds$max_matrix_entries)[2L] else 0.0

  for (r_idx in seq_len(J_est)) {
    s_t_hi_iv <- gmm_iv(s_tilde_upper[r_idx])
    s_t_lo_iv <- gmm_iv(max(0.0, s_tilde_lower[r_idx]))
    s_star_upper[r_idx] <- gmm_iv_mul(s_t_hi_iv, sqrt_one_plus_eps_over_l_iv)[2L]
    s_star_lower[r_idx] <- gmm_iv_mul(s_t_lo_iv, sqrt_one_minus_eps_over_l_iv)[1L]

    if (identical(wt_mode, "optimal")) {
      # Theorem 12.3: ||\tilde\lambda - \lambda^*||_2 <= G_\lambda + \epsilon \tilde s / (\sqrt{l} (l - \epsilon))
      d_stat_iv <- gmm_iv_div(gmm_iv_mul(gmm_iv(eps_V), s_t_hi_iv), denom_loading_stat_iv)
      delta_C_stat[r_idx] <- d_stat_iv[2L]
      tot_c_iv <- gmm_iv_add(gmm_iv(delta_C_num[r_idx]), d_stat_iv)
      total_loading_error_cert[r_idx] <- tot_c_iv[2L]

      # Theorem 12.2: |\tilde\tau - \tau^*| <= G_c + s^* [ \epsilon_m / \sqrt{l - \epsilon} + \epsilon \sqrt{\tilde Q} / \sqrt{l(l - \epsilon)} ]
      term_m_iv <- gmm_iv_div(gmm_iv(eps_m_l2), sqrt_l_m_eps_iv)
      term_q_iv <- gmm_iv_div(gmm_iv_mul(gmm_iv(eps_V), gmm_iv(sqrt_Q_upper)), sqrt_l_prod_iv)
      tot_theta_iv <- gmm_iv_add(
        gmm_iv(delta_theta_solve[r_idx]),
        gmm_iv_mul(gmm_iv(s_star_upper[r_idx]), gmm_iv_add(term_m_iv, term_q_iv))
      )
      total_target_error_cert[r_idx] <- tot_theta_iv[2L]
      point_err_solve[r_idx] <- delta_theta_solve[r_idx]
      point_err_moment[r_idx] <- gmm_iv_mul(gmm_iv(s_star_upper[r_idx]), term_m_iv)[2L]
      point_err_weights[r_idx] <- gmm_iv_mul(gmm_iv(s_star_upper[r_idx]), term_q_iv)[2L]

      # Theorem 12.4: R_c = eps_proj + \nu_I ||\tilde\lambda||_2 + (||\tilde V|| + \epsilon)^{1/2} G_\lambda + \epsilon \tilde s / \sqrt{l(l - \epsilon)}
      eps_proj_iv <- gmm_iv(projection_error[r_idx])
      nu_lam_iv <- gmm_iv_mul(gmm_iv(nu_I), gmm_iv(c_hat_norms[r_idx]))
      v_glam_iv <- gmm_iv_mul(sqrt_omega_star_iv, gmm_iv(delta_C_num[r_idx]))
      w_term_iv <- gmm_iv_div(gmm_iv_mul(gmm_iv(eps_V), s_t_hi_iv), sqrt_l_prod_iv)
      rc_iv <- gmm_iv_add(
        gmm_iv_add(eps_proj_iv, nu_lam_iv),
        gmm_iv_add(v_glam_iv, w_term_iv)
      )
      R_c_cert[r_idx] <- rc_iv[2L]
      rc_err_proj[r_idx] <- eps_proj_iv[2L]
      rc_err_rowif[r_idx] <- nu_lam_iv[2L]
      rc_err_loading[r_idx] <- v_glam_iv[2L]
      rc_err_weights[r_idx] <- w_term_iv[2L]
    } else {
      delta_C_stat[r_idx] <- 0.0
      total_loading_error_cert[r_idx] <- delta_C_num[r_idx]
      c_bar_iv <- gmm_iv_add(gmm_iv(c_hat_norms[r_idx]), gmm_iv(delta_C_num[r_idx]))
      tot_theta_iv <- gmm_iv_add(
        gmm_iv(delta_theta_solve[r_idx]),
        gmm_iv_mul(c_bar_iv, gmm_iv(eps_m_l2))
      )
      total_target_error_cert[r_idx] <- tot_theta_iv[2L]
      point_err_solve[r_idx] <- delta_theta_solve[r_idx]
      point_err_moment[r_idx] <- gmm_iv_mul(c_bar_iv, gmm_iv(eps_m_l2))[2L]
      point_err_weights[r_idx] <- 0.0  # identity weighting: no changing-weight term

      eps_proj_iv <- gmm_iv(projection_error[r_idx])
      nu_lam_iv <- gmm_iv_mul(gmm_iv(nu_I), gmm_iv(c_hat_norms[r_idx]))
      v_glam_iv <- gmm_iv_mul(sqrt_omega_star_iv, gmm_iv(delta_C_num[r_idx]))
      rc_iv <- gmm_iv_add(eps_proj_iv, gmm_iv_add(nu_lam_iv, v_glam_iv))
      R_c_cert[r_idx] <- rc_iv[2L]
      rc_err_proj[r_idx] <- eps_proj_iv[2L]
      rc_err_rowif[r_idx] <- nu_lam_iv[2L]
      rc_err_loading[r_idx] <- v_glam_iv[2L]
      rc_err_weights[r_idx] <- 0.0
    }

    # Theorem 12.4 SE and diagonal variance bounds:
    # |\tilde\sigma_r - \sigma^*_r| <= \epsilon_{acc,\sigma,r} + R_{c,r}
    v_hat_r <- max(0.0, as.double(diag_Sigma[r_idx]))
    sigma_hat_r_iv <- gmm_iv_sqrt(gmm_iv(v_hat_r))
    eps_acc_sigma_iv <- gmm_iv(acc_sigma[r_idx])
    se_err_iv <- gmm_iv_add(eps_acc_sigma_iv, gmm_iv(R_c_cert[r_idx]))
    total_se_error_cert[r_idx] <- se_err_iv[2L]

    # Diagonal variance error: 2 \tilde\sigma_r (\epsilon_{acc,\sigma,r} + R_{c,r}) + (\epsilon_{acc,\sigma,r} + R_{c,r})^2
    var_err_iv <- gmm_iv_add(
      gmm_iv_mul(gmm_iv(2.0), gmm_iv_mul(sigma_hat_r_iv, se_err_iv)),
      gmm_iv_mul(se_err_iv, se_err_iv)
    )
    total_cov_diag_error_cert[r_idx] <- var_err_iv[2L]
  }

  # Theorem 12.4 matrix covariance bound:
  # ||\tilde\Sigma - \Sigma^*||_2 <= \epsilon_{acc,\Sigma} + 2(||\tilde\Sigma|| + \epsilon_{acc,\Sigma})^{1/2} \bar R + \bar R^2
  R_bar_sq_iv <- gmm_iv(0.0)
  for (r_idx in seq_len(J_est)) {
    rc_iv <- gmm_iv(R_c_cert[r_idx])
    R_bar_sq_iv <- gmm_iv_add(R_bar_sq_iv, gmm_iv_mul(rc_iv, rc_iv))
  }
  R_bar_iv <- gmm_iv_sqrt(R_bar_sq_iv)
  Sigma_frob_iv <- gmm_iv(if (J_est > 0L) gmm_iv_norm(unname(Sigma_sand), max_matrix_entries = bnds$max_matrix_entries)[2L] else 0.0)
  eps_acc_Sigma_iv <- gmm_iv(acc_cov)
  sig_plus_acc_sqrt_iv <- gmm_iv_sqrt(gmm_iv_add(Sigma_frob_iv, eps_acc_Sigma_iv))
  total_cov_error_cert <- gmm_iv_add(
    eps_acc_Sigma_iv,
    gmm_iv_add(
      gmm_iv_mul(gmm_iv(2.0), gmm_iv_mul(sig_plus_acc_sqrt_iv, R_bar_iv)),
      R_bar_sq_iv
    )
  )[2L]

  # Enforce explicit accuracy budgets on the complete Theorem 12 transfer bounds
  if (J_est > 0L && acc$target_error_unit %in% c("abs", "both") && is.finite(acc$max_target_error_abs)) {
    if (any(!is.finite(total_target_error_cert)) || max(total_target_error_cert) > acc$max_target_error_abs) {
      cdid_ipw_pipeline_abort(
        paste0(
          "Complete Theorem 12 target error bound (",
          format(max(total_target_error_cert), digits = 6),
          ") exceeds `gmm_accuracy$max_target_error_abs` (",
          format(acc$max_target_error_abs, digits = 6), ")."
        ),
        c("cdid_ipw_pipeline_accuracy_error", "cdid_gmm_accuracy_error"),
        diagnostics = list(
          max_total_target_error_cert = max(total_target_error_cert),
          max_target_error_abs = acc$max_target_error_abs
        )
      )
    }
  }
  if (J_est > 0L && !is.null(acc$max_loading_error_abs) && is.finite(acc$max_loading_error_abs)) {
    if (any(!is.finite(total_loading_error_cert)) || max(total_loading_error_cert) > acc$max_loading_error_abs) {
      cdid_ipw_pipeline_abort(
        paste0(
          "Complete Theorem 12 loading error bound (",
          format(max(total_loading_error_cert), digits = 6),
          ") exceeds `gmm_accuracy$max_loading_error_abs` (",
          format(acc$max_loading_error_abs, digits = 6), ")."
        ),
        c("cdid_ipw_pipeline_accuracy_error", "cdid_gmm_accuracy_error"),
        diagnostics = list(
          max_total_loading_error_cert = max(total_loading_error_cert),
          max_loading_error_abs = acc$max_loading_error_abs
        )
      )
    }
  }

  if (acc$target_error_unit %in% c("model_se","both")) {
    # The accepted criterion-scale policy remains distinct from actual sandwich SE.
    ratios <- vapply(seq_len(J_est),function(j)
      gmm_iv_ratio_upper(total_target_error_cert[j],s_star_lower[j]),numeric(1L))
    if (any(!is.finite(ratios)) || any(ratios > acc$max_target_error_se)) {
      cdid_ipw_pipeline_abort("Complete target error exceeds the supplied criterion-SE budget.", "cdid_ipw_pipeline_accuracy_error")
    }
  }
  structure(
    list(
      certified = TRUE,
      weighting = wt_mode,
      nuisance_problems = nuisance_problems,
      link_enclosures = link_enclosures,
      covariance_transfer = cov_transfer,
      simultaneous_covariance_entries = combined_cov,
      gmm_transfer = list(
        certified = TRUE,
        theorem_ref = "thm:gmm",
        epsilon_m_vec = eps_m_vec,
        epsilon_m = eps_m_l2,
        epsilon_V = eps_V,
        nu_I = nu_I,
        l_lower = l_lower,
        u_upper = u_upper,
        l_minus_eps_V = l_minus_eps_V,
        s_star_lower = s_star_lower,
        s_star_upper = s_star_upper,
        delta_W = delta_W_iv[2L],
        delta_theta_solve = delta_theta_solve,
        delta_C_num = delta_C_num,
        delta_C_stat = delta_C_stat,
        total_loading_error_cert = total_loading_error_cert,
        total_target_error_cert = total_target_error_cert,
        R_c_cert = R_c_cert,
        total_cov_error_cert = total_cov_error_cert,
        total_cov_diag_error_cert = total_cov_diag_error_cert,
        total_se_error_cert = total_se_error_cert,
        raw_root_transfer_implemented = TRUE,
        covariance_transfer_implemented = TRUE,
        projection_error_cert = projection_error,
        accumulation_se_error_cert = acc_sigma,
        accumulation_cov_error_cert = acc_cov,
        # Diagnostic retention (no new bound, threshold or formula): the outward
        # components whose sums form `total_target_error_cert` and `R_c_cert`,
        # and the certified lower eigenvalue used by the optimal-weight terms.
        point_error_components = list(
          solve = point_err_solve,
          moment_change = point_err_moment,
          changing_weights = point_err_weights,
          scope = "Upper endpoints of the three Theorem 12.2 terms (G_c; s* eps_m/sqrt(l-eps); s* eps_V sqrt(Q)/sqrt(l(l-eps))); identity weighting has an exact-zero changing-weight term. Attribution of certified upper bounds, not measured error."
        ),
        R_c_components = list(
          projection = rc_err_proj,
          row_if_change = rc_err_rowif,
          numerical_loading = rc_err_loading,
          changing_weights = rc_err_weights,
          scope = "Upper endpoints of the four Theorem 12.4 R_c terms; identity weighting has an exact-zero changing-weight term."
        ),
        certified_lower_eigenvalue = list(
          value = if (identical(wt_mode, "optimal")) l_lower else NA_real_,
          applies = identical(wt_mode, "optimal"),
          source = if (identical(wt_mode, "optimal")) "gmm_fit$propagation_bounds$l_lambda_min_V_cert" else "identity weighting: no covariance inversion; l_lower placeholder 1.0 is not a certified eigenvalue"
        ),
        complete_inference_certified = FALSE
      ),
      # Per-nuisance stored radius record (already computed certificate scalars).
      nuisance_root_radii = data.frame(
        nuisance_id = nuisance_ids,
        mode = vapply(nuisance_problems, function(pr) as.character(pr$certificate$radius_selection$mode), character(1L)),
        branch = vapply(nuisance_problems, function(pr) { b <- pr$certificate$radius_selection$branch; if (is.null(b)) NA_character_ else as.character(b) }, character(1L)),
        r = vapply(nuisance_problems, function(pr) as.double(pr$raw$r), numeric(1L)),
        delta = vapply(nuisance_problems, function(pr) as.double(pr$certificate$delta), numeric(1L)),
        t = vapply(nuisance_problems, function(pr) as.double(pr$certificate$t), numeric(1L)),
        a = vapply(nuisance_problems, function(pr) as.double(pr$certificate$a), numeric(1L)),
        x = vapply(nuisance_problems, function(pr) as.double(pr$certificate$x), numeric(1L)),
        q = vapply(nuisance_problems, function(pr) as.double(pr$certificate$q), numeric(1L)),
        selfmap = vapply(nuisance_problems, function(pr) as.double(pr$certificate$selfmap), numeric(1L)),
        rho = vapply(nuisance_problems, function(pr) as.double(pr$certificate$rho), numeric(1L)),
        inside_sufficient_region = vapply(nuisance_problems, function(pr) { z <- pr$certificate$radius_selection$inside_sufficient_region; if (is.null(z)) NA else isTRUE(z) }, logical(1L)),
        stringsAsFactors = FALSE
      )
    ),
    class = "cdid_ipw_numerical_transfer_certificate"
  )
}

# ------------------------------------------------------------------------------
# 4. Full-Frame Influence-Function Transport (`lem:transport`)
# ------------------------------------------------------------------------------

transport_if_to_full_frame <- function(Phi_A, n_total, a_to_full) {
  if (!is.matrix(Phi_A) || !is.numeric(Phi_A) || any(!is.finite(Phi_A))) {
    cdid_ipw_pipeline_abort(
      "`Phi_A` must be a finite numeric matrix of admissible-frame influence functions.",
      "cdid_ipw_pipeline_input_error"
    )
  }
  n_a <- nrow(Phi_A)
  R_cols <- ncol(Phi_A)
  n_tot <- as.integer(n_total)
  idx_full <- as.integer(a_to_full)

  if (n_a < 1L || n_tot < n_a || length(idx_full) != n_a ||
      anyNA(idx_full) || any(idx_full < 1L | idx_full > n_tot) ||
      anyDuplicated(idx_full) > 0L) {
    cdid_ipw_pipeline_abort(
      "`a_to_full` must be a strictly valid 1-to-1 index map from `1..n_A` into `1..n_total`.",
      "cdid_ipw_pipeline_input_error"
    )
  }

  scale_fac <- as.double(n_tot) / as.double(n_a)
  Phi_full <- matrix(0.0, nrow = n_tot, ncol = R_cols, dimnames = list(NULL, colnames(Phi_A)))
  Phi_full[idx_full, ] <- scale_fac * Phi_A

  list(
    Phi_full = Phi_full,
    n_total = n_tot,
    n_a_units = n_a,
    scale_factor = scale_fac,
    a_to_full = idx_full
  )
}

# ------------------------------------------------------------------------------
# 5. Base-Period Reporting Map (`sec:maps` in `cp02_population_aggregation_if.tex`)
# ------------------------------------------------------------------------------

apply_base_period_reporting_map <- function(
  theta_hat, target_if_full, Sigma_theta_a, target_meta,
  base_period = c("varying", "universal"), anticipation = 0L,
  calendar, enclosure = NULL
) {
  bp <- match.arg(base_period)
  R <- nrow(target_meta)
  if (length(theta_hat) != R || ncol(target_if_full) != R ||
      nrow(Sigma_theta_a) != R || ncol(Sigma_theta_a) != R) {
    cdid_ipw_pipeline_abort("Dimension mismatch in reporting map.", "cdid_ipw_pipeline_input_error")
  }
  cal <- sort(unique(as.double(calendar)))
  maps <- list(); metas <- list(); bases <- numeric()
  append <- function(meta, row, base) {
    maps[[length(maps)+1L]] <<- row
    metas[[length(metas)+1L]] <<- meta
    bases[length(bases)+1L] <<- base
  }
  for (j in seq_len(R)) {
    m <- target_meta[j,,drop=FALSE]; t <- m$time; g <- m$g
    if (bp == "varying" && t == cal[1L]) next
    row <- numeric(R); row[j] <- 1.0
    base <- m$ref
    # Released did2.5.1 tests report < g, including anticipation dates.
    if (bp == "varying" && t < g) {
      base <- cal[match(t,cal)-1L]
      if (base != m$ref) {
        key <- m$functional_key
        hit <- which(target_meta$g == g & target_meta$time == base &
          (is.na(target_meta$functional_key) |
           (!is.na(key) & target_meta$functional_key == key)))
        if (length(hit) != 1L) {
          cdid_ipw_pipeline_abort("Varying-base contrast lacks one compatible previous-period coordinate.",
            "cdid_ipw_pipeline_reporting_support_error", list(target=m, base=base))
        }
        row[hit] <- -1.0
      }
    }
    append(m,row,base)
  }
  for (g in unique(target_meta$g)) {
    pos <- which(target_meta$g == g); ref <- unique(target_meta$ref[pos])
    if (length(ref) != 1L) cdid_ipw_pipeline_abort("Cohort has multiple reference dates.", "cdid_ipw_pipeline_reporting_support_error")
    m <- target_meta[pos[1L],,drop=FALSE]
    m$time <- ref; m$role <- "reference"; m$functional_key <- NA_character_
    m$target_id <- paste0("reference(g=",g,",t=",ref,")")
    m$param_id <- m$target_id; m$param_index <- NA_integer_
    if (bp == "universal") {
      append(m,numeric(R),ref)
    } else if (ref != cal[1L]) {
      base <- cal[match(ref,cal)-1L]
      # A distinct qualified early operator is a distinct pretrend coordinate.
      hit <- which(target_meta$g == g & target_meta$time == base)
      if (!length(hit)) cdid_ipw_pipeline_abort("Varying reference contrast is unsupported by the catalog.", "cdid_ipw_pipeline_reporting_support_error")
      for (h in hit) {
        row <- numeric(R); row[h] <- -1.0
        mh <- m; mh$role <- "pre"; mh$functional_key <- target_meta$functional_key[h]
        mh$target_id <- paste0(m$target_id,";fn=",target_meta$functional_key[h])
        append(mh,row,base)
      }
    }
  }
  L <- do.call(rbind,maps); meta <- do.call(rbind,metas)
  o <- order(meta$g,meta$time,meta$functional_key,na.last=TRUE)
  L <- L[o,,drop=FALSE]; meta <- meta[o,,drop=FALSE]; bases <- bases[o]
  rownames(meta) <- NULL
  meta$base_period_mode <- bp; meta$base_t <- bases; meta$is_post <- meta$time >= meta$g
  th <- as.double(L %*% theta_hat); names(th) <- meta$target_id
  rows <- target_if_full %*% t(L); colnames(rows) <- meta$target_id
  S <- L %*% Sigma_theta_a %*% t(L); dimnames(S) <- list(meta$target_id,meta$target_id)
  enc <- if (is.null(enclosure)) NULL else list(
    theta = ipw_pipeline_interval_product(enclosure$theta$lo,enclosure$theta$hi,t(L)),
    rows = ipw_pipeline_interval_product(enclosure$rows$lo,enclosure$rows$hi,t(L)))
  list(base_period=bp,L_map=L,target_meta=meta,theta_hat=th,if_full=rows,Sigma_a=S,enclosure=enc)
}

# ------------------------------------------------------------------------------
# 5a. D5 requested-cell ledger and post-contrast basis (catalog-relative support)
# ------------------------------------------------------------------------------
# The ledger is the declared request (target cohorts x requested post dates),
# normalized before any support-dependent step. It never derives cells from
# retained rows or successful fits. Statuses are catalog-relative: an
# unavailable cell is "not supported by the retained comparison catalog", never
# a population non-identification claim. Proof labels: prop:graph (component
# estimability), prop:def (unanchored coordinates are not ATT levels), D5.

ipw_pipeline_cell_id <- function(g, t) {
  paste0("ATT(g=", fmt_design_num(as.double(g)), ",t=", fmt_design_num(as.double(t)), ")")
}

ipw_pipeline_request_ledger <- function(target_cohorts, calendar, cohort_universe, anticipation,
                                        control_group, target_post_times = NULL, requested_post = NULL) {
  cal <- sort(unique(as.double(calendar)))
  coh <- sort(unique(as.double(cohort_universe)))
  tg <- sort(unique(as.double(target_cohorts)))
  delta <- as.double(anticipation)
  if (is.null(target_post_times)) {
    # did-style declared grid: every stored calendar date at or after g.
    parts <- lapply(tg, function(g) {
      tt <- cal[cal >= g]
      if (!length(tt)) NULL else data.frame(g = g, time = tt, requested_by = "calendar_grid", stringsAsFactors = FALSE)
    })
    parts <- parts[!vapply(parts, is.null, logical(1L))]
    rows <- if (length(parts)) do.call(rbind, parts) else data.frame(g = numeric(), time = numeric(), requested_by = character(), stringsAsFactors = FALSE)
  } else {
    tp <- target_post_times
    tc <- if (is.data.frame(tp) && "time" %in% names(tp)) "time" else "t"
    if (!is.data.frame(tp) || !all(c("g", tc) %in% names(tp)) || nrow(tp) < 1L ||
        !is.numeric(tp$g) || !is.numeric(tp[[tc]]) || anyNA(tp$g) || anyNA(tp[[tc]]) ||
        any(!is.finite(tp$g)) || any(!is.finite(tp[[tc]]))) {
      cdid_ipw_pipeline_abort("`target_post_times` must be a nonempty data.frame with finite numeric columns `g` and `t`/`time`.",
        c("cdid_ipw_pipeline_target_error", "cdid_ipw_pipeline_input_error"))
    }
    rows <- data.frame(g = as.double(tp$g), time = as.double(tp[[tc]]), requested_by = "target_post_times", stringsAsFactors = FALSE)
  }
  if (!is.null(requested_post)) {
    if (!is.data.frame(requested_post) || !all(c("g", "t") %in% names(requested_post)) ||
        !is.numeric(requested_post$g) || !is.numeric(requested_post$t) ||
        anyNA(requested_post$g) || anyNA(requested_post$t)) {
      cdid_ipw_pipeline_abort("`requested_post` must be NULL or a data.frame with finite numeric columns `g` and `t`.",
        c("cdid_ipw_pipeline_target_error", "cdid_ipw_pipeline_input_error"))
    }
    rows <- rbind(rows, data.frame(g = as.double(requested_post$g), time = as.double(requested_post$t),
      requested_by = "requested_post", stringsAsFactors = FALSE))
  }
  if (nrow(rows) < 1L) {
    cdid_ipw_pipeline_abort("No requested post cell: no stored calendar date at or after any target cohort's treatment date.",
      c("cdid_ipw_pipeline_target_error", "cdid_catalog_target_error"))
  }
  if (any(!rows$g %in% tg)) {
    cdid_ipw_pipeline_abort("Requested post cells must belong to the declared target cohorts.",
      c("cdid_ipw_pipeline_target_error", "cdid_ipw_pipeline_input_error"))
  }
  if (any(!rows$time %in% cal)) {
    cdid_ipw_pipeline_abort("Requested post dates must be stored dates of the full calendar.",
      c("cdid_ipw_pipeline_target_error", "cdid_design_target_error", "cdid_ipw_pipeline_input_error"))
  }
  if (any(rows$time < rows$g)) {
    cdid_ipw_pipeline_abort("Requested post cells require t >= g.",
      c("cdid_ipw_pipeline_target_error", "cdid_ipw_pipeline_input_error"))
  }
  key <- paste(match(rows$g, tg), match(rows$time, cal), sep = ":")
  first <- which(!duplicated(key))
  by <- vapply(first, function(i) paste(sort(unique(rows$requested_by[key == key[i]])), collapse = "+"), character(1L))
  rows <- rows[first, , drop = FALSE]
  rows$requested_by <- by
  o <- order(rows$g, rows$time, method = "radix")
  rows <- rows[o, , drop = FALSE]
  rownames(rows) <- NULL
  refs <- vapply(rows$g, function(g) tryCatch(calendar_reference_date(g, cal, delta), cdid_error = function(e) NA_real_), numeric(1L))
  pools <- lapply(seq_len(nrow(rows)), function(i) {
    e <- eligible_control_cohorts(rows$g[i], rows$time[i], coh, delta)
    if (identical(control_group, "nevertreated")) e[is.infinite(e)] else e
  })
  data.frame(
    ledger_index = seq_len(nrow(rows)),
    cell_id = ipw_pipeline_cell_id(rows$g, rows$time),
    g = rows$g, time = rows$time, ref = refs, event_time = rows$time - rows$g, role = "post",
    requested_by = rows$requested_by,
    has_clean_control = lengths(pools) > 0L,
    clean_control_pool = vapply(pools, function(p) if (length(p)) paste0("{", paste(fmt_design_num(p), collapse = ","), "}") else "{}", character(1L)),
    stringsAsFactors = FALSE
  )
}

# Join the ledger to the catalog and the (full) design by cohort/date identity.
ipw_pipeline_ledger_support <- function(ledger, metadata_obj = NULL, design = NULL, global_status = NULL) {
  R <- nrow(ledger)
  ledger$cohort_has_rows <- FALSE
  ledger$n_catalog_rows_touching <- 0L
  ledger$n_qualified_rows_touching <- 0L
  ledger$catalog_target_status <- NA_character_
  ledger$catalog_request_reason <- NA_character_
  ledger$design_param_index <- NA_integer_
  ledger$design_param_id <- NA_character_
  ledger$component_id <- NA_integer_
  ledger$is_isolated <- NA
  ledger$level_estimable <- FALSE
  if (!is.null(metadata_obj)) {
    rows <- metadata_obj$catalog$rows
    ledger$cohort_has_rows <- ledger$g %in% rows$g
    for (i in seq_len(R)) {
      same_g <- rows$g == ledger$g[i]
      ledger$n_catalog_rows_touching[i] <- sum(same_g & (rows$s == ledger$time[i] | rows$t == ledger$time[i]))
      ledger$n_qualified_rows_touching[i] <- sum(same_g & ((rows$s == ledger$time[i] & rows$s_role == "post_conditional") |
                                                           (rows$t == ledger$time[i] & rows$t_role == "post_conditional")))
    }
    ct <- metadata_obj$catalog$targets
    ledger$catalog_target_status <- ct$status[match(ledger$g, ct$g)]
    ru <- metadata_obj$catalog$requested_unsupported
    if (is.data.frame(ru) && nrow(ru) > 0L) {
      hit <- match(paste(ledger$g, ledger$time), paste(ru$g, ru$t))
      ledger$catalog_request_reason[!is.na(hit)] <- ru$reason[hit[!is.na(hit)]]
    }
  }
  if (!is.null(design)) {
    # Adapter use without a catalog object: the design's own link cohorts are the row-bearing cohorts.
    if (is.null(metadata_obj)) ledger$cohort_has_rows <- ledger$g %in% as.double(design$links$g)
    pt <- design$post_targets
    hit <- match(paste(ledger$g, ledger$time), paste(pt$g, pt$time))
    ok <- !is.na(hit)
    ledger$design_param_index[ok] <- as.integer(pt$param_index[hit[ok]])
    ledger$design_param_id[ok] <- as.character(pt$param_id[hit[ok]])
    ledger$component_id[ok] <- as.integer(pt$component_id[hit[ok]])
    ledger$is_isolated[ok] <- as.logical(pt$is_isolated[hit[ok]])
    ledger$level_estimable[ok] <- as.logical(pt$level_estimable[hit[ok]])
    design_status <- rep(NA_character_, R)
    design_status[ok] <- as.character(pt$status[hit[ok]])
  } else {
    design_status <- rep(NA_character_, R)
  }
  status <- character(R); reason <- character(R)
  for (i in seq_len(R)) {
    pool_txt <- paste0("clean control pool at t=", fmt_design_num(ledger$time[i]), " under the declared control group: ", ledger$clean_control_pool[i])
    if (!is.null(global_status)) {
      status[i] <- global_status$status
      reason[i] <- paste0(global_status$reason, "; ", pool_txt)
    } else if (!ledger$cohort_has_rows[i]) {
      status[i] <- "rowless_cohort"
      reason[i] <- paste0("no retained comparison catalog row exists for cohort g=", fmt_design_num(ledger$g[i]),
        if (!is.na(ledger$catalog_target_status[i])) paste0(" (catalog target status: ", ledger$catalog_target_status[i], ")") else "",
        "; ", pool_txt)
    } else if (is.na(design_status[i])) {
      status[i] <- "not_represented"
      reason[i] <- paste0("cell is not represented in the structural design; ", pool_txt)
    } else if (identical(design_status[i], "estimable")) {
      status[i] <- "estimable"
      reason[i] <- if (identical(design$calendar_spec$pt_assumption, "all"))
        "level anchored to the reference or a represented PT-all pre-date through retained comparison rows (prop:graph)" else
        "level anchored to the cohort reference through retained comparison rows (prop:graph)"
    } else if (identical(design_status[i], "isolated_unsupported")) {
      status[i] <- "isolated_unsupported"
      reason[i] <- paste0("no retained comparison row has (g=", fmt_design_num(ledger$g[i]), ",t=", fmt_design_num(ledger$time[i]),
        ") as a qualified post endpoint; catalog rows touching this date: ", ledger$n_catalog_rows_touching[i],
        ", qualified post rows touching it: ", ledger$n_qualified_rows_touching[i], "; ", pool_txt,
        if (!is.na(ledger$catalog_request_reason[i])) paste0("; catalog request reason: ", ledger$catalog_request_reason[i]) else "")
    } else {
      status[i] <- "disconnected_component"
      reason[i] <- paste0("cell lies in design component ", ledger$component_id[i],
        ", which no retained comparison row connects to the cohort reference; its level is a gauge coordinate, not an ATT estimate (prop:def); estimable contrasts inside the component are retained in the post-contrast basis; ", pool_txt)
    }
  }
  ledger$support_status <- status
  ledger$support_reason <- reason
  ledger$available <- status == "estimable"
  ledger$population_claim <- "none: catalog-relative support only; not a population non-identification statement"
  ledger
}

# Sparse estimable post-contrast basis: within every unanchored component with
# q >= 2 post vertices, the q - 1 differences from its lowest-index post vertex.
ipw_pipeline_post_basis <- function(design) {
  params <- design$parameters
  un <- as.integer(design$rank_diagnostics$unanchored_components)
  out <- list()
  for (u in un) {
    members <- sort(which(params$component_id == u & params$role == "post"))
    if (length(members) < 2L) next
    base <- members[1L]
    for (j in members[-1L]) {
      out[[length(out) + 1L]] <- data.frame(
        basis_id = paste0("contrast:", params$param_id[j], "-", params$param_id[base]),
        component_id = u, from_param_index = as.integer(base), to_param_index = as.integer(j),
        from_id = params$param_id[base], to_id = params$param_id[j],
        g = params$g[j], time_from = params$time[base], time_to = params$time[j],
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) {
    return(data.frame(basis_id = character(), component_id = integer(), from_param_index = integer(),
      to_param_index = integer(), from_id = character(), to_id = character(), g = numeric(),
      time_from = numeric(), time_to = numeric(), stringsAsFactors = FALSE))
  }
  b <- do.call(rbind, out)
  rownames(b) <- NULL
  b
}

ipw_pipeline_core_target_plan <- function(ledger, design_used) {
  pt_used <- design_used$post_targets
  used_hit <- match(paste(ledger$g, ledger$time), paste(pt_used$g, pt_used$time))
  if (!identical(is.na(used_hit), is.na(ledger$design_param_index)) ||
      any(pt_used$status[used_hit[!is.na(used_hit)]] != ledger$support_status[!is.na(used_hit)])) {
    cdid_ipw_pipeline_abort(
      "Structural reduction changed a requested cell's catalog-relative support status; the reduction is not lossless.",
      "cdid_ipw_pipeline_integrity_error"
    )
  }
  ledger$core_param_index <- NA_integer_
  ledger$core_param_index[!is.na(used_hit)] <- as.integer(pt_used$param_index[used_hit[!is.na(used_hit)]])
  basis <- ipw_pipeline_post_basis(design_used)
  n_basis <- nrow(basis)

  param_df <- design_used$parameters
  est_params <- param_df[param_df$level_estimable, , drop = FALSE]
  rownames(est_params) <- NULL
  by_coh <- design_used$rank_diagnostics$by_cohort
  est_params$ref <- by_coh$ref[match(est_params$g, by_coh$g)]
  est_params$target_id <- est_params$param_id
  unavailable_rows <- which(!ledger$available & !is.na(ledger$core_param_index))
  # Core target table: level-estimable parameters, represented unavailable post
  # levels (returned NA with catalog-relative reasons) and the contrast basis.
  core_targets <- data.frame(
    target_id = c(est_params$target_id, ledger$cell_id[unavailable_rows], basis$basis_id),
    kind = c(rep("level", nrow(est_params)), rep("level_unavailable", length(unavailable_rows)), rep("basis", n_basis)),
    ledger_index = c(match(est_params$param_id, ledger$design_param_id), unavailable_rows, rep(NA_integer_, n_basis)),
    stringsAsFactors = FALSE
  )
  # The immutable core requires a nonempty target specification even when its
  # estimable loading set is empty. If only rowless cells were requested and
  # no post basis/anchored coordinate exists, validate one actual unavailable
  # design coordinate. It remains an internal NA query, never a reported ATT.
  internal_unavailable <- if (!nrow(core_targets)) which(!param_df$level_estimable)[1L] else NA_integer_
  if (!is.na(internal_unavailable)) {
    core_targets <- rbind(core_targets, data.frame(
      target_id = param_df$param_id[internal_unavailable], kind = "internal_unavailable",
      ledger_index = NA_integer_, stringsAsFactors = FALSE
    ))
  }
  target_specs <- c(
    lapply(seq_len(nrow(est_params)), function(j) gmm_target_spec(as.character(est_params$target_id[j]), as.integer(est_params$param_index[j]), 1.0)),
    lapply(unavailable_rows, function(i) gmm_target_spec(ledger$cell_id[i], as.integer(ledger$core_param_index[i]), 1.0)),
    lapply(seq_len(n_basis), function(b) gmm_target_spec(basis$basis_id[b], c(basis$from_param_index[b], basis$to_param_index[b]), c(-1.0, 1.0)))
  )
  if (!is.na(internal_unavailable)) {
    target_specs <- c(target_specs, list(gmm_target_spec(
      param_df$param_id[internal_unavailable], as.integer(param_df$param_index[internal_unavailable]), 1.0
    )))
  }
  list(
    ledger = ledger,
    basis = basis,
    n_basis = n_basis,
    est_params = est_params,
    unavailable_rows = unavailable_rows,
    core_targets = core_targets,
    target_specs = target_specs,
    J_est_expected = nrow(est_params) + n_basis
  )
}

# ------------------------------------------------------------------------------
# 5b. Direct did-style auxiliary reporting over the full ledger template
# ------------------------------------------------------------------------------
# D4/D18 and moment-catalog contract Section 4: direct auxiliary pretrend cells
# use released cellwise controls. The union of both base-period modes is fitted
# independently of post GMM; a reporting selector changes neither its rows nor
# its weights. D5: every requested cohort owns a ledger template. A pair without
# clean controls or without declared pair support (`cdid_ipw_support_error`)
# is recorded with dates, controls and counts and skipped. Nuisance fitting,
# normalization, overlap, numerical, resource and inference failures propagate.
# `aux_spec` is the pipeline's one union plan (built before its reservations);
# NULL rebuilds it here for standalone callers. `compute_if = FALSE`
# (compute_inffunc = FALSE) fits and reports every auxiliary point estimate but
# replays no auxiliary IF and forms no reported IF matrix or covariance.
build_ipw_auxiliary_reporting <- function(ledger, calendar, cohorts, anticipation, control_group,
  registry, frame_owner, gmm, control, root_witness, numerical_transfer, transfer_requested,
  base_period, n_total, a_to_full, retained_pipeline_entries, custom_est_method = NULL, custom_dots = list(),
  pt_assumption = "post", pt_start = NULL, pt_catalog_links = NULL, aux_spec = NULL, compute_if = TRUE,
  retained_store = NULL) {
  if (is.null(aux_spec)) aux_spec <- ipw_auxiliary_spec(ledger, calendar, cohorts, anticipation, control_group)
  if (!isTRUE(compute_if) && isTRUE(transfer_requested))
    cdid_ipw_pipeline_abort("Numerical transfer certifies inference and needs the auxiliary influence functions.",
      "cdid_ipw_pipeline_control_error")
  descriptors <- aux_spec$descriptors; selections <- aux_spec$selections
  ant <- as.double(anticipation)
  Q <- length(descriptors)
  if (Q > control$catalog_bounds$max_rows) cdid_ipw_pipeline_abort("Auxiliary reporting union exceeds explicit row bound.", "cdid_ipw_pipeline_bound_error")
  n_a <- length(a_to_full)
  has_raw_aux <- !is.null(registry$raw_ipw) && length(registry$raw_aux_points) == Q
  nvs <- if (is.null(registry)) list() else if (has_raw_aux) registry$raw_aux_nuisance_views else registry$nuisance_views
  fits <- if (is.null(registry)) list() else if (has_raw_aux) registry$raw_aux_fitted_views else registry$fitted_views
  all_points <- if (is.null(registry)) list() else registry$point_links
  points <- vector("list", Q); problems <- list(); enclosures <- vector("list", Q)
  custom_aux_retained <- 0.0
  desc_status <- rep("supported", Q); desc_reason <- rep(NA_character_, Q); desc_counts <- vector("list", Q)
  for (j in seq_len(Q)) {
    d <- descriptors[[j]]; lv <- NULL; idx <- NA_integer_; support_err <- NULL
    if (!is.null(registry$cp04_reader)) {
      idx <- registry$cp04_aux_index[j]
      support_err <- registry$cp04_aux_errors[[j]]
      if (!is.na(idx)) { points[[j]] <- registry$cp04_points[[idx]]; next }
    } else if (has_raw_aux) {
      support_err <- registry$raw_aux_errors[[j]]
      if (is.null(support_err)) {
        if (!is.null(registry$raw_aux_fit_errors[[j]])) stop(registry$raw_aux_fit_errors[[j]])
        points[[j]] <- registry$raw_aux_points[[j]]
        next
      }
    } else if (is.function(custom_est_method)) {
      hit_pk <- integer(0L)
      for (k in seq_along(all_points)) {
        pk <- all_points[[k]]
        if (identical(as.double(pk$cohorts$g), as.double(d$g)) &&
            identical(as.double(pk$pair$s), as.double(d$s)) &&
            identical(as.double(pk$pair$t), as.double(d$t)) &&
            identical(as.double(pk$cohorts$K), as.double(d$K))) {
          hit_pk <- k
          break
        }
      }
      if (length(hit_pk) == 1L) {
        points[[j]] <- all_points[[hit_pk]]
        next
      }
      trial <- tryCatch(construct_ipw_link_view(frame_owner, d$s, d$t, g = d$g, K = d$K, nuisance_view = NULL,
        control = control$owner_control, link_id = paste0("reporting:", j),
        functional_key = paste0("reporting_operator:", gsub(" ", "", paste(d$g, paste(d$K, collapse = ","), sep = ":"))),
        require_support = TRUE, verify_nuisance_mask = FALSE),
        cdid_ipw_support_error = identity)
      if (inherits(trial, "error")) {
        support_err <- trial
      } else {
        lv <- trial
      }
    } else {
      for (k in seq_along(nvs)) {
        if (identical(nvs[[k]]$cohorts$g, as.double(d$g)) && identical(nvs[[k]]$cohorts$K, d$K)) {
          trial <- tryCatch(construct_ipw_link_view(frame_owner, d$s, d$t, g = d$g, K = d$K, nuisance_view = nvs[[k]],
            control = control$owner_control, link_id = paste0("reporting:", j),
            functional_key = paste0("reporting_operator:", gsub(" ", "", paste(d$g, paste(d$K, collapse = ","), sep = ":"))),
            require_support = TRUE, verify_nuisance_mask = TRUE),
            cdid_ipw_support_error = identity, cdid_ipw_shared_mask_error = identity)
          if (!inherits(trial, "error")) { lv <- trial; idx <- k; break }
          if (inherits(trial, "cdid_ipw_support_error")) { support_err <- trial; break }
        }
      }
      if (is.null(support_err) && is.na(idx)) {
        nv <- tryCatch(construct_ipw_nuisance_view(frame_owner, d$g, d$K, d$s, d$t, control$owner_control,
          nuisance_id = paste0("reporting_nuisance:", j)), cdid_ipw_support_error = identity)
        if (inherits(nv, "error")) {
          support_err <- nv
        } else {
          fv <- fit_ipw_nuisance_view(nv, control$glm_control, control$fitted_control)
          nvs[[length(nvs) + 1L]] <- nv; fits[[length(fits) + 1L]] <- fv; idx <- length(nvs)
          lv <- construct_ipw_link_view(frame_owner, d$s, d$t, g = d$g, K = d$K, nuisance_view = nv,
            control = control$owner_control, link_id = paste0("reporting:", j), functional_key = paste0("reporting_operator:", j),
            require_support = TRUE, verify_nuisance_mask = TRUE)
        }
      }
    }
    if (!is.null(support_err)) {
      # Only absent declared pair/cohort support reaches this branch. Record the
      # actual dates, controls and outcome-free counts; do not fit anything.
      desc_status[j] <- "unsupported_pair"
      desc_reason[j] <- paste0(conditionMessage(support_err), " [g=", fmt_design_num(d$g), ", s=", fmt_design_num(d$s),
        ", t=", fmt_design_num(d$t), ", K={", paste(fmt_design_num(d$K), collapse = ","), "}]")
      if (has_raw_aux && !is.null(registry$raw_aux_counts[[j]])) {
        desc_counts[[j]] <- registry$raw_aux_counts[[j]]
      } else {
        probe <- construct_ipw_nuisance_view(frame_owner, d$g, d$K, d$s, d$t, control$owner_control,
          nuisance_id = paste0("reporting_probe:", j), require_support = FALSE)
        desc_counts[[j]] <- list(counts_available = TRUE,
          n_a_units = probe$support$n_a_units, n_observed_pair_all = probe$support$n_observed_pair_all,
          n_observed_risk_units = probe$support$n_observed_risk_units,
          by_cohort_a_count = probe$support$by_cohort_a_count, by_cohort_pair_count = probe$support$by_cohort_pair_count)
      }
      next
    }
    if (is.function(custom_est_method)) {
      cb_pre <- preflight_custom_ipw_link_view(frame_owner, lv, control$point_control)
      ipw_owner_safe_sum(c(retained_pipeline_entries, custom_aux_retained, cb_pre$peak_concurrent_entries),
        min(control$owner_control$max_matrix_entries, control$point_control$max_matrix_entries),
        "simultaneous retained pipeline, auxiliary custom links and callback buffers")
      pl <- fit_custom_ipw_link_view(
        frame_owner = frame_owner,
        link_view = lv,
        control = control$point_control,
        custom_est_method = custom_est_method,
        custom_dots = custom_dots
      )
      points[[j]] <- pl
      custom_aux_retained <- custom_aux_retained + cb_pre$retained_entries_per_link
      next
    }
    pl <- fit_ipw_point_link_view(fits[[idx]], lv, control$point_control); points[[j]] <- pl
    if (isTRUE(transfer_requested)) {
      prob <- NULL
      if (!is.null(numerical_transfer)) for (k in seq_along(all_points)) if (identical(pl$fitted_nuisance_ref, all_points[[k]]$fitted_nuisance_ref)) {
        h <- which(vapply(numerical_transfer$link_enclosures, function(x) identical(x$link_ref, all_points[[k]]$link_view_ref), logical(1L)))
        if (length(h)) { prob <- get("nuisance_problem", numerical_transfer$link_enclosures[[h[1L]]]$problem_ref, inherits = FALSE); break }
      }
      if (is.null(prob)) {
        saved <- which(vapply(problems, function(x) identical(x$fit_ref, pl$fitted_nuisance_ref), logical(1L)))
        if (length(saved)) prob <- problems[[saved[1L]]]$problem else {
          witness <- root_witness(fits[[idx]])
          if (!is.list(witness) || !identical(sort(names(witness)), c("F", "r"))) cdid_ipw_pipeline_abort("Auxiliary root witness needs F and r.", "cdid_ipw_pipeline_control_error")
          prob <- bind_fitted_nuisance_to_problem_view(fits[[idx]], F = witness$F, r = witness$r,
            certificate_control = control$certificate_control, problem_control = control$problem_control, reader_control = control$owner_control)
          if (!isTRUE(prob$certified)) cdid_ipw_pipeline_abort("Auxiliary nuisance root enclosure is inconclusive.", "cdid_ipw_pipeline_transfer_inconclusive")
          problems[[length(problems) + 1L]] <- list(fit_ref = pl$fitted_nuisance_ref, problem = prob)
        }
      }
      enc <- bind_point_link_to_enclosure_view(prob, pl, link_control = control$link_control, reader_control = control$owner_control, include_candidates = TRUE)
      if (!isTRUE(enc$certified)) cdid_ipw_pipeline_abort("Auxiliary link enclosure is inconclusive.", "cdid_ipw_pipeline_transfer_inconclusive")
      enclosures[[j]] <- enc
    }
  }
  sup <- which(desc_status == "supported"); Qs <- length(sup)
  A <- NULL; lo <- hi <- if (transfer_requested) matrix(0.0, n_total, Qs) else NULL
  theta <- if (Qs) vapply(points[sup], `[[`, numeric(1L), "estimate") else numeric(0L); theta_lo <- theta_hi <- theta
  preflight <- NULL; backend_ids <- integer(0L)
  if (Qs && isTRUE(compute_if)) {
    sp <- points[sup]; se_ <- enclosures[sup]
    if (!is.null(registry$cp04) || has_raw_aux) {
      # Only requested reporting columns are materialized, directly into `out`
      # below. No dense F x union-of-auxiliary-moments IF is constructed.
      aux_backend <- if (!is.null(registry$cp04)) registry$cp04 else registry$raw_ipw
      requested <- selections[[base_period]]
      active_columns <- which(vapply(requested, function(d) !is.na(d$index) && d$index > 0L &&
        desc_status[d$index] == "supported", logical(1)))
      q_requested <- vapply(requested[active_columns], `[[`, integer(1), "index")
      backend_ids <- vapply(points[q_requested], `[[`, integer(1), "cp04_index")
      if (length(backend_ids)) {
        preflight <- cp04_preflight_links(aux_backend, unique(backend_ids),
          retained_pipeline_entries + 8 * n_total * (length(requested) + nrow(ledger)))
        ipw_owner_safe_sum(c(retained_pipeline_entries, preflight$peak_concurrent_entries,
          8 * n_total * (length(requested) + nrow(ledger))), aux_backend$control$max_matrix_entries,
          "CP04 auxiliary producer and requested reporting outputs")
      }
    } else if (is.function(custom_est_method)) {
      preflight <- preflight_custom_ipw_if_view(sp, control$point_control)
      ipw_owner_safe_sum(c(retained_pipeline_entries, preflight$peak_concurrent_entries, 8.0 * n_total * Qs),
        control$owner_control$max_matrix_entries, "simultaneous auxiliary producer and retained outputs (plus interval reader when requested)")
      A <- matrix(0.0, n_a, Qs)
      replay_custom_ipw_if_view(sp, control$point_control, callback_schema = TRUE, callback = function(batch) {
        ii <- batch$a_index; A[ii, ] <<- batch$phi
      })
    } else {
      preflight <- preflight_ipw_point_if_view(sp, control$point_control)
      interval_retained <- if (!isTRUE(transfer_requested)) 0.0 else sum(vapply(se_, function(x) x$preflight$total_retained_entries, numeric(1L))) +
        sum(vapply(problems, function(x) x$problem$preflight$total_retained_entries, numeric(1L)))
      ipw_owner_safe_sum(c(retained_pipeline_entries, preflight$peak_concurrent_entries, interval_retained,
        (if (transfer_requested) 32.0 else 8.0) * n_total * Qs), control$owner_control$max_matrix_entries, "simultaneous auxiliary producer and retained outputs (plus interval reader when requested)")
      A <- matrix(0.0, n_a, Qs)
      replay_ipw_point_if_view(sp, control$point_control, callback_schema = TRUE, callback = function(batch) {
        ii <- batch$a_index; A[ii, ] <<- batch$phi
        if (isTRUE(transfer_requested)) for (j in seq_len(Qs)) {
          pr <- get("nuisance_problem", se_[[j]]$problem_ref, inherits = FALSE)
          chunk_pre <- preflight_ipw_link_if_view_chunk(b_len = length(ii), p = pr$raw$p_covariates, m = length(pr$raw$cohorts),
            control = control$link_control, reader_control = control$owner_control,
            certificate_control = pr$controls$certificate_control, include_components = FALSE, include_shared_xi = FALSE,
            has_phi_candidate = TRUE, has_xi_beta = FALSE, read_outcomes = TRUE,
            panel_owner = get("panel_owner", frame_owner$panel_ref, inherits = FALSE), a_seq = ii,
            expected_a_id = batch$id, caller_phi_entries = length(ii), caller_xi_entries = 0)
          ipw_owner_safe_sum(c(retained_pipeline_entries, interval_retained, preflight$peak_concurrent_entries,
            chunk_pre$combined_concurrent_entries, 32.0 * n_total * Qs), control$owner_control$max_matrix_entries,
            "simultaneous auxiliary point producer, compact interval reader, nuisance cache and retained post/report outputs")
          ch <- enclose_ipw_link_if_view_chunk(se_[[j]], min(ii), max(ii), control$link_control,
            control$owner_control, phi_candidate = unname(batch$phi[, j]), include_components = FALSE, include_shared_xi = FALSE)
          if (!isTRUE(ch$certified)) cdid_ipw_pipeline_abort("Auxiliary IF enclosure is inconclusive.", "cdid_ipw_pipeline_transfer_inconclusive")
          fac <- gmm_iv_div(gmm_iv(as.double(n_total)), gmm_iv(as.double(n_a)))
          for (h in seq_along(ii)) {
            v <- gmm_iv_mul(ch$phi_star_intersected[h, ], fac)
            lo[a_to_full[ii[h]], j] <<- v[1L]; hi[a_to_full[ii[h]], j] <<- v[2L]
          }
        }
      })
      if (isTRUE(transfer_requested)) for (j in seq_len(Qs)) {
        theta_lo[j] <- se_[[j]]$enclosures_star$mu[1L]; theta_hi[j] <- se_[[j]]$enclosures_star$mu[2L]
      }
    }
  }
  full <- if (!is.null(registry$cp04) || has_raw_aux || !isTRUE(compute_if)) NULL else
    if (Qs) transport_if_to_full_frame(A, n_total, a_to_full)$Phi_full else matrix(0.0, n_total, 0L)
  # Reported rows: direct pretrend/reference cells of the selected base period
  # plus every declared post ledger cell. Post targets remain the original GMM
  # targets; unavailable cells are explicit NA columns, never zero.
  specs <- selections[[base_period]]; R_post <- nrow(ledger)
  R <- length(specs) + R_post
  out <- if (isTRUE(compute_if)) matrix(NA_real_, n_total, R) else NULL
  th <- rep(NA_real_, R); available <- logical(R)
  if (!is.null(registry$cp04) && Qs && length(backend_ids)) {
    ids <- unique(backend_ids)
    signs <- vapply(specs[active_columns], `[[`, numeric(1), "sign")
    # CP08 retained moment IF: an admitted auxiliary extension serves the exact
    # same-recipe columns from the retained raw store (frame scale on the
    # F = A, lambda = 1 frame, the existing sign, the original output slot) and
    # replays only the unmatched columns, in their original order, through the
    # ordinary producer under the planner's own partition with the store
    # charged. Without an extension the original replay runs unchanged.
    aux_retained <- retained_pipeline_entries + 8 * n_total * (length(requested) + nrow(ledger))
    service <- cp04_retained_if_auxiliary_service(retained_store, registry$cp04, backend_ids, aux_retained)
    if (is.null(service)) {
      cp04_replay_links(registry$cp04, ids, function(b) {
        for (j in seq_along(active_columns)) out[b$rows, active_columns[j]] <<-
          signs[j] * b$phi[, match(backend_ids[j], ids)]
      }, scale = "frame", chunk_rows = preflight$chunk_rows)
    } else {
      if (!isTRUE(preflight$chunk_rows == service$chunk_rows))
        cp04_abort("The auxiliary replay row partition differs from its admitted original plan.", "integrity")
      served <- which(!is.na(service$column)); streamed <- which(is.na(service$column))
      cp04_retained_if_auxiliary_slices(retained_store, registry$cp04, backend_ids, function(b) {
        for (j in served) out[b$rows, active_columns[j]] <<-
          signs[j] * b$phi[, match(backend_ids[j], b$links)]
      })
      if (length(streamed)) {
        ids_stream <- unique(backend_ids[streamed])
        pre_stream <- cp04_preflight_links(registry$cp04, ids_stream, aux_retained + retained_store$entries,
          chunk_rows = service$chunk_rows)
        ipw_owner_safe_sum(c(retained_pipeline_entries, pre_stream$peak_concurrent_entries,
          8 * n_total * (length(requested) + nrow(ledger)), retained_store$entries),
          control$owner_control$max_matrix_entries,
          "CP04 auxiliary producer, retained moment IF and requested reporting outputs")
        if (!identical(ids_stream, service$ids_stream) || !isTRUE(pre_stream$chunk_rows == service$chunk_rows_stream))
          cp04_abort("The unmatched auxiliary replay differs from its admitted plan.", "integrity")
        cp04_replay_links(registry$cp04, ids_stream, function(b) {
          for (j in streamed) out[b$rows, active_columns[j]] <<-
            signs[j] * b$phi[, match(backend_ids[j], ids_stream)]
        }, scale = "frame", chunk_rows = pre_stream$chunk_rows)
      }
    }
  } else if (has_raw_aux && Qs && length(backend_ids)) {
    ids <- unique(backend_ids)
    signs <- vapply(specs[active_columns], `[[`, numeric(1), "sign")
    aux_retained <- retained_pipeline_entries + 8 * n_total * (length(requested) + nrow(ledger))
    service <- cp04_retained_if_auxiliary_service(retained_store, registry$raw_ipw, backend_ids, aux_retained)
    A_req <- matrix(0.0, nrow = n_a, ncol = length(active_columns))
    if (is.null(service)) {
      cp04_replay_links(registry$raw_ipw, ids, function(b) {
        for (j in seq_along(active_columns)) A_req[b$rows, j] <<- b$phi[, match(backend_ids[j], ids)]
      }, scale = "frame", chunk_rows = preflight$chunk_rows)
    } else {
      if (!isTRUE(preflight$chunk_rows == service$chunk_rows))
        cp04_abort("The auxiliary replay row partition differs from its admitted original plan.", "integrity")
      served <- which(!is.na(service$column)); streamed <- which(is.na(service$column))
      cp04_retained_if_auxiliary_slices(retained_store, registry$raw_ipw, backend_ids, function(b) {
        for (j in served) A_req[b$rows, j] <<- b$phi[, match(backend_ids[j], b$links)]
      })
      if (length(streamed)) {
        ids_stream <- unique(backend_ids[streamed])
        pre_stream <- cp04_preflight_links(registry$raw_ipw, ids_stream, aux_retained + retained_store$entries,
          chunk_rows = service$chunk_rows)
        ipw_owner_safe_sum(c(retained_pipeline_entries, pre_stream$peak_concurrent_entries,
          8 * n_total * (length(requested) + nrow(ledger)), retained_store$entries),
          registry$raw_ipw$control$max_matrix_entries,
          "CP04 auxiliary producer, retained moment IF and requested reporting outputs")
        if (!identical(ids_stream, service$ids_stream) || !isTRUE(pre_stream$chunk_rows == service$chunk_rows_stream))
          cp04_abort("The unmatched auxiliary replay differs from its admitted plan.", "integrity")
        cp04_replay_links(registry$raw_ipw, ids_stream, function(b) {
          for (j in streamed) A_req[b$rows, j] <<- b$phi[, match(backend_ids[j], ids_stream)]
        }, scale = "frame", chunk_rows = pre_stream$chunk_rows)
      }
    }
    full_req <- transport_if_to_full_frame(A_req, n_total, a_to_full)$Phi_full
    for (j in seq_along(active_columns)) {
      out[, active_columns[j]] <- signs[j] * full_req[, j]
    }
  }
  el <- eh <- if (transfer_requested) matrix(NA_real_, 1L, R) else NULL
  rl <- rh <- if (transfer_requested) out else NULL; metas <- vector("list", R)
  meta_row <- function(g, time, role, target_id, ref, base_t, is_post, ledger_index, auxiliary_index, status, reason, loading_column, functional_key) {
    data.frame(g = g, time = time, role = role, target_id = target_id, param_id = target_id, param_index = NA_integer_,
      functional_key = functional_key, ref = ref, base_t = base_t, base_period_mode = base_period, is_post = is_post,
      ledger_index = ledger_index, auxiliary_index = auxiliary_index, support_status = status, support_reason = reason,
      loading_column = loading_column, stringsAsFactors = FALSE)
  }
  for (j in seq_along(specs)) {
    d <- specs[[j]]
    if (identical(d$status, "reference_zero")) {
      metas[[j]] <- meta_row(d$g, d$t, "reference", paste0("reference(g=", fmt_design_num(d$g), ",t=", fmt_design_num(d$t), ")"), d$ref, d$from, FALSE,
        NA_integer_, NA_integer_, "deterministic_reference_zero", "universal base-period reference cell: exact zero with zero influence function", NA_integer_, NA_character_)
      th[j] <- 0.0; if (isTRUE(compute_if)) out[, j] <- 0.0; available[j] <- TRUE
      if (isTRUE(transfer_requested)) { el[j] <- eh[j] <- 0.0; rl[, j] <- rh[, j] <- 0.0 }
      next
    }
    tid <- paste0("pre(g=", fmt_design_num(d$g), ",t=", fmt_design_num(d$t), ")")
    if (identical(d$status, "unavailable")) {
      metas[[j]] <- meta_row(d$g, d$t, "pre", tid, d$ref, d$from, FALSE, NA_integer_, NA_integer_, "unavailable_pretrend", d$reason, NA_integer_, NA_character_)
      next
    }
    q <- d$index
    if (!identical(desc_status[q], "supported")) {
      metas[[j]] <- meta_row(d$g, d$t, "pre", tid, d$ref, d$from, FALSE, NA_integer_, as.integer(q), "unavailable_pretrend",
        paste0("unsupported_pair: ", desc_reason[q]), NA_integer_, paste0("auxiliary:", q))
      next
    }
    qs <- match(q, sup)
    metas[[j]] <- meta_row(d$g, d$t, "pre", tid, d$ref, d$from, FALSE, NA_integer_, as.integer(q), "supported_direct_pretrend",
      paste0("direct auxiliary comparison ", fmt_design_num(d$from), "->", fmt_design_num(d$t), " with controls {", paste(fmt_design_num(d$K), collapse = ","), "}"),
      NA_integer_, paste0("auxiliary:", q))
    th[j] <- d$sign * theta[qs]
    if (is.null(registry$cp04) && !has_raw_aux && isTRUE(compute_if)) out[, j] <- d$sign * full[, qs]
    available[j] <- TRUE
    if (isTRUE(transfer_requested)) {
      z <- gmm_iv_mul(gmm_iv(d$sign), c(theta_lo[qs], theta_hi[qs])); el[j] <- z[1L]; eh[j] <- z[2L]
      for (i in seq_len(n_total)) { z <- gmm_iv_mul(gmm_iv(d$sign), c(lo[i, qs], hi[i, qs])); rl[i, j] <- z[1L]; rh[i, j] <- z[2L] }
    }
  }
  for (h in seq_len(R_post)) {
    j <- length(specs) + h
    col <- if (is.null(gmm)) NA_integer_ else ledger$loading_column[h]
    metas[[j]] <- meta_row(ledger$g[h], ledger$time[h], "post", ledger$cell_id[h], ledger$ref[h], ledger$ref[h], TRUE,
      as.integer(h), NA_integer_, ledger$support_status[h], ledger$support_reason[h], col, NA_character_)
    if (!is.na(col)) {
      th[j] <- gmm$estimates[col]; if (isTRUE(compute_if)) out[, j] <- gmm$if_full[, col]; available[j] <- TRUE
      if (isTRUE(transfer_requested)) { el[j] <- gmm$enclosure$theta$lo[col]; eh[j] <- gmm$enclosure$theta$hi[col]; rl[, j] <- gmm$enclosure$rows$lo[, col]; rh[, j] <- gmm$enclosure$rows$hi[, col] }
    }
  }
  meta <- do.call(rbind, metas); o <- order(meta$g, meta$time); meta <- meta[o, , drop = FALSE]; rownames(meta) <- NULL
  th <- th[o]; available <- available[o]; names(th) <- meta$target_id
  if (isTRUE(compute_if)) { out <- out[, o, drop = FALSE]; colnames(out) <- meta$target_id }
  meta <- chained_pt_reporting_tags(meta, pt_assumption, pt_start, ant, descriptors, pt_catalog_links)
  meta$available <- available
  meta$reported_column <- seq_len(R)
  av <- which(available)
  # Complete reported joint covariance on the supported block from the full common frame.
  S <- if (isTRUE(compute_if)) matrix(NA_real_, R, R, dimnames = list(meta$target_id, meta$target_id)) else NULL
  if (length(av) && isTRUE(compute_if)) {
    centered <- sweep(out[, av, drop = FALSE], 2L, colMeans(out[, av, drop = FALSE]), "-")
    S[av, av] <- crossprod(centered) * as.double(n_a) / (as.double(n_total)^2)
  }
  E <- if (!isTRUE(transfer_requested)) NULL else list(theta = list(lo = el[, o, drop = FALSE], hi = eh[, o, drop = FALSE]),
    rows = list(lo = rl[, o, drop = FALSE], hi = rh[, o, drop = FALSE]))
  SE <- NULL
  if (!is.null(E) && length(av)) {
    sub <- list(theta = list(lo = E$theta$lo[, av, drop = FALSE], hi = E$theta$hi[, av, drop = FALSE]),
      rows = list(lo = E$rows$lo[, av, drop = FALSE], hi = E$rows$hi[, av, drop = FALSE]))
    si <- ipw_pipeline_inference_intervals(sub, NULL)$vcov
    for (k in seq_along(si$lo)) { v <- gmm_iv_mul(c(si$lo[k], si$hi[k]), gmm_iv(as.double(n_a))); si$lo[k] <- v[1L]; si$hi[k] <- v[2L] }
    SE <- list(lo = S, hi = S); SE$lo[av, av] <- si$lo; SE$hi[av, av] <- si$hi
  }
  desc_table <- if (Q) data.frame(index = seq_len(Q), g = vapply(descriptors, `[[`, numeric(1L), "g"),
    s = vapply(descriptors, `[[`, numeric(1L), "s"), t = vapply(descriptors, `[[`, numeric(1L), "t"),
    K = vapply(descriptors, function(x) paste0("{", paste(fmt_design_num(x$K), collapse = ","), "}"), character(1L)),
    status = desc_status, reason = desc_reason, stringsAsFactors = FALSE) else
    data.frame(index = integer(), g = numeric(), s = numeric(), t = numeric(), K = character(), status = character(), reason = character(), stringsAsFactors = FALSE)
  list(base_period = base_period, target_meta = meta, theta_hat = th, if_full = out, available = available,
    Sigma_a = S, Sigma_a_enclosure = SE, enclosure = E, L_map = NULL,
    auxiliary = list(descriptors = descriptors, selections = selections, point_links = points, descriptor_status = desc_table,
      unsupported_counts = desc_counts, fitted_views = fits, nuisance_views = nvs, new_root_problems = problems,
      link_enclosures = enclosures, preflight = preflight, post_GMM_rows_added = 0L,
      n_descriptors = Q, n_supported = Qs, n_unsupported = Q - Qs,
      scope = "D4/D18 direct auxiliary reporting; no added post restriction or GMM weighting row; unsupported pairs are catalog-relative diagnostics, not population statements"))
}

# ------------------------------------------------------------------------------
# 6. Analytical & One-Way Cluster Inference + Multiplier Bootstrap
# ------------------------------------------------------------------------------

# Flag composition (public binding, October 16, 2026). The released did 2.5.1
# flags are distinct: `bstrap` selects multiplier IQR standard errors, `cband`
# selects the simultaneous sup-t critical value. `bstrap = NULL` keeps the
# existing private law for existing callers (bootstrap iff `cband`). An explicit
# `bstrap` separates the two without changing either meaning:
#   bstrap=FALSE,cband=FALSE  analytic SE, pointwise normal critical value, no draw;
#   bstrap=TRUE ,cband=FALSE  multiplier IQR SE, pointwise normal critical value;
#   bstrap=TRUE ,cband=TRUE   multiplier IQR SE, sup-t critical value;
#   bstrap=FALSE,cband=TRUE   bound by the released ENDPOINT (D6 familiar did
#                             behaviour; flag_reference/, att_gt.R:582-590,688-695,
#                             compute.aggte.R:365-370,509-514,642-647, getSE:811-835):
#       endpoint "cells"      (att_gt reported cells, post basis, simple
#                             aggregation, overall summaries): analytic SE, normal
#                             pointwise critical value, NO draw, no warning;
#       endpoint "aggte_bins" (group/dynamic/calendar bin families): analytic
#                             reported SE, multiplier sup-t critical value whose
#                             draws are standardized by their OWN IQR scale
#                             (released mboot convention), caller RNG advanced,
#                             released warning issued by the aggregation.
#   No upstream tiny-SE filtering, draw deletion or pointwise fallback is
#   imported (D5/D12): a refused multiplier result stays a refusal.
# RNG: generated signs either use the existing explicit-seed-and-restore
# contract (`seed`), or, with `caller_rng = TRUE`, the caller's current RNG
# stream, which is advanced honestly and never set or restored (the familiar
# `set.seed()` workflow). Fixed signs use no RNG. No silent seed exists.
ipw_pipeline_flag_law <- function(which = NULL) {
  laws <- c(private = "private_bootstrap_iff_cband",
            explicit = "explicit_bstrap_and_cband",
            cells_false_true = "released_bstrap_false_cband_true_cells_analytic_pointwise_no_draws",
            aggte_bins_false_true = "released_bstrap_false_cband_true_aggte_bins_analytic_se_multiplier_sup_t_band",
            empty = "empty_selection_no_draws")
  if (is.null(which)) laws else unname(laws[[which]])
}

ipw_pipeline_resolve_bootstrap_flags <- function(bstrap, cband, where = "inference", endpoint = c("cells", "aggte_bins", "simple")) {
  endpoint <- match.arg(endpoint)
  if (!gmm_scalar_logical(cband)) cdid_ipw_pipeline_abort(paste0("`cband` must be a scalar logical (", where, ")."), "cdid_ipw_pipeline_input_error")
  if (!is.null(bstrap) && !gmm_scalar_logical(bstrap)) cdid_ipw_pipeline_abort(paste0("`bstrap` must be NULL, TRUE or FALSE (", where, ")."), "cdid_ipw_pipeline_input_error")
  if (is.null(bstrap)) {
    return(list(bootstrap = isTRUE(cband), simultaneous = isTRUE(cband), draws = isTRUE(cband),
                law = ipw_pipeline_flag_law("private"), endpoint = endpoint, released_false_true = FALSE))
  }
  if (isTRUE(bstrap) || !isTRUE(cband)) {
    return(list(bootstrap = isTRUE(bstrap), simultaneous = isTRUE(cband) && !identical(endpoint, "simple"), draws = isTRUE(bstrap),
                law = ipw_pipeline_flag_law("explicit"), endpoint = endpoint, released_false_true = FALSE))
  }
  # Explicit bstrap = FALSE, cband = TRUE: the released endpoint distinction.
  if (endpoint %in% c("cells", "simple")) {
    list(bootstrap = FALSE, simultaneous = FALSE, draws = FALSE,
         law = ipw_pipeline_flag_law("cells_false_true"), endpoint = endpoint, released_false_true = TRUE)
  } else {
    list(bootstrap = FALSE, simultaneous = TRUE, draws = TRUE,
         law = ipw_pipeline_flag_law("aggte_bins_false_true"), endpoint = endpoint, released_false_true = TRUE)
  }
}

# Released compute.aggte warning text (l. 368/512/646), issued once per
# group/dynamic/calendar aggregation whose band draws actually ran.
ipw_pipeline_released_band_warning <- function() "Used bootstrap procedure to compute simultaneous confidence band"

# Authoritative universal-reference rows, independent of evaluated estimates/IFs.
ipw_pipeline_known_reference_rows <- function(meta) {
  fields <- c("role", "support_status", "time", "ref", "base_t")
  if (!is.data.frame(meta) || !all(fields %in% names(meta))) return(rep(FALSE, NROW(meta)))
  z <- meta$role == "reference" & meta$support_status == "deterministic_reference_zero" &
    meta$time == meta$ref & meta$time == meta$base_t
  !is.na(z) & z
}

compute_chained_ipw_inference <- function(
  estimates, if_full, cluster_full = NULL, alp = 0.05, cband = FALSE,
  biters = 1000L, btype = "rademacher", seed = NULL,
  multiplier_budget, inference_accuracy, enclosure = NULL,
  fixed_signs = NULL, eligible_full = NULL, bstrap = NULL, caller_rng = FALSE,
  endpoint = c("cells", "aggte_bins", "simple"), target_meta = NULL
) {
  est <- as.double(estimates); Psi <- as.matrix(if_full)
  n <- nrow(Psi); J <- ncol(Psi)
  mb <- validate_multiplier_budget(multiplier_budget)
  ia <- validate_ipw_pipeline_inference_accuracy(inference_accuracy)
  if (length(est) != J || n < 1L || any(!is.finite(Psi)) || any(!is.finite(est)) ||
      !gmm_scalar_numeric(alp, positive=TRUE) || alp >= 1.0) {
    cdid_ipw_pipeline_abort("Invalid inference inputs.", "cdid_ipw_pipeline_input_error")
  }
  if (!identical(btype,"rademacher")) {
    cdid_ipw_pipeline_abort("Only explicit Rademacher or caller-supplied fixed +/-1 signs are implemented.", "cdid_ipw_pipeline_bootstrap_error")
  }
  # D20's exception is a known reference OPERATOR, not a zero evaluation.
  # These fields come from build_ipw_auxiliary_reporting's reference_zero branch.
  # Private absolute inference retains its historical deterministic-zero handling.
  structural_reference <- integer(0L)
  if (ipw_pipeline_accuracy_is_relative(ia) && !is.null(target_meta)) {
    fields <- c("target_id", "role", "support_status", "time", "ref", "base_t")
    if (!is.data.frame(target_meta) || nrow(target_meta) != J || !all(fields %in% names(target_meta)) ||
        (!is.null(names(estimates)) && !identical(as.character(target_meta$target_id), names(estimates))))
      cdid_ipw_pipeline_abort("D20 reference metadata must align with the inference targets.", "cdid_ipw_pipeline_input_error")
    known <- ipw_pipeline_known_reference_rows(target_meta)
    # A nonempty weighted family of known reference identities remains exactly
    # zero, including the estimated-share derivative. Carry its source proof;
    # aggregate labels or a zero numerical evaluation never establish this role.
    if (all(c("reference_operator", "reference_sources") %in% names(target_meta))) {
      for (j in seq_len(J)) if (identical(target_meta$role[j], "reference_aggregate") &&
          identical(target_meta$support_status[j], "deterministic_reference_zero") &&
          identical(target_meta$reference_operator[j], "weighted_reference_identity")) {
        src <- target_meta$reference_sources[[j]]
        known[j] <- is.data.frame(src) && nrow(src) > 0L && all(ipw_pipeline_known_reference_rows(src))
      }
    }
    structural_reference <- which(known)
  }
  flags <- ipw_pipeline_resolve_bootstrap_flags(bstrap, cband, endpoint = endpoint)
  if (!gmm_scalar_logical(caller_rng)) cdid_ipw_pipeline_abort("`caller_rng` must be a scalar logical.", "cdid_ipw_pipeline_input_error")
  if (isTRUE(caller_rng) && !is.null(seed)) {
    cdid_ipw_pipeline_abort("`caller_rng = TRUE` uses the caller's RNG stream; an explicit `seed` cannot also be supplied.", "cdid_ipw_pipeline_bootstrap_error")
  }
  if (!is.null(cluster_full) && (length(cluster_full) != n || anyNA(cluster_full))) {
    cdid_ipw_pipeline_abort("Cluster IDs must cover the full unit frame.", "cdid_ipw_pipeline_cluster_error")
  }
  u_cl <- if (is.null(cluster_full)) seq_len(n) else unique(cluster_full)
  map <- if (is.null(cluster_full)) seq_len(n) else match(cluster_full,u_cl)
  G <- length(u_cl)
  if (!is.null(cluster_full) && G < 2L) cdid_ipw_pipeline_abort("At least two declared clusters are required.", "cdid_ipw_pipeline_cluster_error")
  if (is.null(eligible_full)) eligible_full <- rep(TRUE,n)
  if (!is.logical(eligible_full) || length(eligible_full)!=n || anyNA(eligible_full)) {
    cdid_ipw_pipeline_abort("Eligibility must be declared for every full-frame unit.", "cdid_ipw_pipeline_input_error")
  }
  # Check producer matrices before allocation. The accepted consumer does its
  # own workspace preflight as well; its allowance is supplied by the caller.
  infer_entries <- if (is.null(enclosure)) 4.0*n*J + 4.0*G*J + 4.0*J*J else 12.0*n*J + 12.0*G*J + 8.0*J*J
  if (n > mb$max_units || G > mb$max_clusters || J > mb$max_targets ||
      infer_entries > mb$max_workspace_cells || 8.0*infer_entries > mb$max_workspace_bytes) {
    cdid_ipw_pipeline_abort("Inference producer/interval/consumer storage exceeds the explicit multiplier budget.", "cdid_ipw_pipeline_bound_error")
  }
  S <- matrix(0.0,G,J)
  if (is.null(cluster_full) && typeof(Psi) == "double" && !is.object(Psi)) {
    S[,] <- 0.0 + Psi
  } else {
    for (i in seq_len(n)) S[map[i],] <- S[map[i],] + Psi[i,]
  }
  raw_cluster_sums <- S
  S <- sweep(S,2L,colMeans(S),"-")
  V <- crossprod(S)/(as.double(n)^2)
  dimnames(V) <- list(names(estimates),names(estimates))
  if (any(!is.finite(S)) || any(!is.finite(V)))
    cdid_ipw_pipeline_abort("Nonfinite inference covariance or cluster sums.", "cdid_ipw_pipeline_scale_error")
  analytical_se <- sqrt(pmax(0.0,diag(V)))
  se <- analytical_se
  deterministic <- which(colSums(abs(Psi)) == 0.0 & est == 0.0)
  # A reference-only selection has exact zero intervals and needs no quantile.
  zcrit <- if (length(deterministic) == J) 0.0 else stats::qnorm(1.0-alp/2.0)
  crit <- zcrit; mult <- NULL; rng <- NULL
  if (!is.null(enclosure) && is.null(ia))
    cdid_ipw_pipeline_abort("An explicit private accuracy budget is required with numerical enclosures.", "cdid_ipw_pipeline_control_error")
  cert <- if (is.null(enclosure)) NULL else ipw_pipeline_inference_intervals(enclosure,cluster_full)
  if (!is.null(cert)) {
    cert$multiplier_error <- 0.0
    # Retained for the D20 scale-aware check (no new bound): the analytic
    # cluster-sum SE interval (the raw multiplier draws' own dimension is
    # n/sqrt(G) times this scale), per-target multiplier errors, the exact
    # structural-zero set and the counts that define the raw-draw dimension.
    cert$analytic_se <- cert$se
    cert$multiplier_error_by_target <- rep(0.0, J)
    cert$raw_draw_scale_lower <- rep(NA_real_, J)
    cert$deterministic <- as.integer(deterministic)
    if (ipw_pipeline_accuracy_is_relative(ia)) cert$structural_reference <- structural_reference
    cert$n_total <- n; cert$n_clusters <- G
  }
  pos <- which(analytical_se > 0.0)
  if (any(!seq_len(J) %in% c(pos,deterministic))) {
    cdid_ipw_pipeline_abort("Nonconstant targets require positive marginal inference scales; covariance rank alone is not a refusal.", "cdid_ipw_pipeline_scale_error")
  }
  # Draws run when the law needs them: bootstrap SEs (bstrap = TRUE), the private
  # law (bootstrap iff cband), or the released aggte-bin band under
  # bstrap = FALSE, cband = TRUE. SE selection never follows draw presence.
  if (flags$draws && length(pos)) {
    B <- validate_positive_int_scalar(biters,"biters",matrix_dimension=TRUE)
    if (B > mb$max_draws) cdid_ipw_pipeline_abort("Draw count exceeds caller multiplier budget.", "cdid_ipw_pipeline_bound_error")
    br <- as.integer(min(mb$max_sign_batch_rows, floor(mb$max_sign_batch_cells/G), B))
    ir <- as.integer(min(mb$max_if_batch_rows, floor(mb$max_if_batch_cells/length(pos)), n))
    # Declared multiplier workers (inherited pl/cores) add one round of whole
    # worker and parent buffers; signs are still drawn here, in batch order.
    producer <- infer_entries + (if (is.null(cert)) 4.0 else 12.0)*B*length(pos) + 4.0*br*G +
                (if (is.null(fixed_signs)) 0.0 else length(fixed_signs)) +
                multiplier_parallel_cells(multiplier_budget_workers(mb), G, length(pos), br)
    if (br < 1L || ir < 1L || producer > mb$max_workspace_cells || 8.0*producer > mb$max_workspace_bytes) {
      cdid_ipw_pipeline_abort("Simultaneous IF/sign/interval/bootstrap producer and consumer payload exceeds caller limits.", "cdid_ipw_pipeline_bound_error")
    }
    if (!is.null(fixed_signs)) {
      if (!is.matrix(fixed_signs) || !is.numeric(fixed_signs) ||
          !identical(dim(fixed_signs),c(as.integer(B),as.integer(G))) || anyNA(fixed_signs) ||
          any(!fixed_signs %in% c(-1,1))) {
        cdid_ipw_pipeline_abort("fixed_signs must be a biters by declared-clusters matrix of +/-1.", "cdid_ipw_pipeline_bootstrap_error")
      }
      if ((!is.null(colnames(fixed_signs)) && !identical(colnames(fixed_signs),as.character(u_cl))) ||
          (!is.null(rownames(fixed_signs)) && !identical(rownames(fixed_signs),as.character(seq_len(B))))) {
        cdid_ipw_pipeline_abort("Fixed-sign dimnames must identify the declared cluster/draw order.", "cdid_ipw_pipeline_bootstrap_error")
      }
      rng <- list(mode="caller_fixed_signs",kind=NULL,seed=NULL,batch_rows=br)
    } else if (isTRUE(caller_rng)) {
      # Public-style request: draw from the caller's current stream and leave
      # it advanced. Nothing is seeded, saved or restored here.
      rng <- list(mode="caller_rng_stream_Rademacher",kind=RNGkind(),seed=NULL,batch_rows=br,
                  ordering="batches in increasing draw order; column-major within each batch",
                  caller_rng_restored=FALSE,
                  note="signs drawn from the caller's RNG stream (set.seed() before the call reproduces them); the stream is advanced, never reset")
    } else {
      if (!gmm_scalar_numeric(seed) || seed != floor(seed) || abs(seed) > .Machine$integer.max) {
        cdid_ipw_pipeline_abort("Generated Rademacher signs require an explicit integer seed (or `caller_rng = TRUE` for the caller's stream).", "cdid_ipw_pipeline_bootstrap_error")
      }
      old_exists <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
      if (old_exists) old_seed <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
      on.exit(if (old_exists) assign(".Random.seed",old_seed,.GlobalEnv) else
        if (exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
      rng <- list(mode="R_sample_Rademacher",kind=RNGkind(),seed=as.integer(seed),batch_rows=br,
                  ordering="batches in increasing draw order; column-major within each batch",caller_rng_restored=TRUE)
      set.seed(as.integer(seed))
    }
    rawlo <- rawhi <- if (is.null(cert)) NULL else matrix(0.0,B,length(pos))
    ids <- seq_len(n); cids <- seq_len(G); dids <- seq_len(B)
    tids <- if (is.null(names(estimates))) paste0("t",pos) else names(estimates)[pos]
    mult <- chained_cluster_multiplier(
      target_ids=tids,expected_unit_ids=ids,
      unit_cluster_map=list(unit_ids=ids,cluster_ids=as.integer(map)),cluster_ids=cids,
      n_denom=as.double(n),center=TRUE,
      if_batch_fn=function(k) {
        ii <- seq.int((k-1L)*ir+1L,min(k*ir,n))
        list(unit_ids=ids[ii],target_ids=tids,if_matrix=unname(Psi[ii,pos,drop=FALSE]))
      },n_if_batches=as.integer(ceiling(n/ir)),
      sign_batch_fn=function(k) {
        ii <- seq.int((k-1L)*br+1L,min(k*br,B))
        signs <- if (is.null(fixed_signs)) matrix(sample(c(-1.0,1.0),length(ii)*G,replace=TRUE),length(ii),G) else unname(fixed_signs[ii,,drop=FALSE])
        if (!is.null(cert)) {
          zz <- ipw_pipeline_interval_product(signs,signs,cert$cluster$lo[,pos,drop=FALSE],cert$cluster$hi[,pos,drop=FALSE])
          fac <- gmm_iv_div(gmm_iv(1.0),gmm_iv_sqrt(gmm_iv(as.double(G))))
          for (j in seq_along(pos)) for (h in seq_along(ii)) {
            z <- gmm_iv_mul(c(zz$lo[h,j],zz$hi[h,j]),fac)
            rawlo[ii[h],j] <<- z[1L]; rawhi[ii[h],j] <<- z[2L]
          }
        }
        list(draw_ids=dids[ii],cluster_ids=cids,sign_matrix=signs)
      },n_sign_batches=as.integer(ceiling(B/br)),expected_draw_ids=dids,n_draws=as.integer(B),
      alpha=as.double(alp),return_raw_draws=TRUE,budget=mb)
    if (!isTRUE(mult$accepted) || !is.finite(mult$crit_val)) {
      cdid_ipw_pipeline_abort("Supplied/generated multiplier inference is refused; no analytical-band fallback.", "cdid_ipw_pipeline_bootstrap_error",list(status=mult$target_status))
    }
    # did2.5.1 uses its normal-consistent bootstrap IQR scale when bootstrap SEs
    # are requested (bstrap = TRUE); under the released aggte-bin FALSE/TRUE law
    # the reported SEs stay analytic (getSE) while the draws serve the band only.
    if (flags$bootstrap) se[pos] <- mult$bootstrap_se
    # The sup-t quantile is the critical value only when a band is requested;
    # bootstrap SEs with cband = FALSE keep the pointwise normal constant. The
    # multiplier standardizes each draw by its own IQR scale (mboot.R:157-180),
    # never by the reported analytic SE.
    crit <- if (flags$simultaneous) mult$crit_val else zcrit
    if (!is.null(cert)) {
      mult_dist <- ipw_pipeline_interval_distance(mult$raw_draws,rawlo,rawhi)
      cert$multiplier_error <- max(mult_dist)
      cert$multiplier_error_by_target[pos] <- apply(as.matrix(mult_dist), 2L, max)
      # Raw-draw dimension: (1/sqrt(G)) sum_g sign_g S_gj has scale n se_j / sqrt(G)
      # with se_j the analytic cluster-sum SE; the lower endpoint is the D20 scale.
      raw_fac <- gmm_iv_div(gmm_iv(as.double(n)), gmm_iv_sqrt(gmm_iv(as.double(G))))
      for (j in seq_along(pos)) cert$raw_draw_scale_lower[pos[j]] <- gmm_iv_mul(raw_fac, cert$analytic_se[,pos[j]])[1L]
      scaleiv <- matrix(0.0,2L,length(pos))
      normalizer <- gmm_iv(stats::qnorm(.75)-stats::qnorm(.25))
      for (j in seq_along(pos)) {
        q25 <- c(sort(rawlo[,j])[ceiling(.25*B)],sort(rawhi[,j])[ceiling(.25*B)])
        q75 <- c(sort(rawlo[,j])[ceiling(.75*B)],sort(rawhi[,j])[ceiling(.75*B)])
        scaleiv[,j] <- gmm_iv_div(gmm_iv_sub(q75,q25),normalizer)
        if (scaleiv[1L,j]<=0) cdid_ipw_pipeline_abort("Positive bootstrap marginal scale is not numerically certified.", "cdid_ipw_pipeline_transfer_inconclusive")
        # Reported-SE certificate follows the reported SE's law: bootstrap IQR
        # scale when bstrap = TRUE, otherwise the analytic cluster-sum interval
        # already in `cert$se` (ipw_pipeline_inference_intervals).
        if (flags$bootstrap) cert$se[,pos[j]] <- gmm_iv_mul(scaleiv[,j],gmm_iv_div(gmm_iv_sqrt(gmm_iv(as.double(G))),gmm_iv(as.double(n))))
      }
      # Under the released aggte-bin FALSE/TRUE law the multiplier scales used for
      # the critical value are recorded apart from the analytic reported-SE
      # certificate (both certified; neither replaces the other). Existing laws
      # keep their certificate structure unchanged.
      if (!flags$bootstrap) {
        cert$multiplier_scale <- matrix(NA_real_,2L,J); cert$multiplier_scale[,pos] <- scaleiv
        cert$se_law <- "analytic_cluster_sum_reported_se; multiplier IQR scales only standardize the sup-t draws"
      }
      if (flags$simultaneous) {
        sup_lo <- sup_hi <- numeric(B)
        for (i in seq_len(B)) for (j in seq_along(pos)) {
          ri <- c(rawlo[i,j],rawhi[i,j])
          ab <- c(if(ri[1L]<=0 && ri[2L]>=0) 0.0 else min(abs(ri)),max(abs(ri)))
          zi <- gmm_iv_div(ab,scaleiv[,j])
          sup_lo[i] <- max(sup_lo[i],zi[1L]); sup_hi[i] <- max(sup_hi[i],zi[2L])
        }
        qi <- ceiling((1.0-alp)*B)
        cert$critical_value <- c(sort(sup_lo)[qi],sort(sup_hi)[qi])
      }
    }
  }
  band_family <- if (length(deterministic) == J) "deterministic_exact_zero_over_certified_ids" else
    if (flags$simultaneous && !is.null(mult)) "simultaneous_multiplier_sup_t_over_certified_ids" else
    "pointwise_normal_over_certified_ids"
  if (is.null(ia)) band_family <- sub("certified_ids", "supported_ids", band_family, fixed = TRUE)
  names(se) <- names(estimates)
  if (!is.null(cert)) {
    if (is.null(cert$critical_value)) cert$critical_value <- gmm_iv(zcrit)
    # Exact zero reference coordinates do not participate in scale inversions.
    if (length(deterministic)) cert$se[,deterministic] <- 0.0
    # Keep both endpoints of each bound, rather than a scalar radius prediction.
    blo <- bhi <- ulo <- uhi <- matrix(0.0,1L,J)
    for (j in seq_len(J)) {
      width <- gmm_iv_mul(cert$critical_value,cert$se[,j]); ti <- c(cert$theta$lo[j],cert$theta$hi[j])
      z <- gmm_iv_sub(ti,width); blo[j]<-z[1L];bhi[j]<-z[2L]
      z <- gmm_iv_add(ti,width); ulo[j]<-z[1L];uhi[j]<-z[2L]
    }
    cert$band <- list(lo=blo,hi=bhi);cert$band_upper <- list(lo=ulo,hi=uhi)
    # Pointwise interval endpoints est -/+ zcrit * se with the same accepted
    # primitives (zcrit is the stored normal constant, as for the band).
    clo <- chi <- culo <- cuhi <- matrix(0.0,1L,J)
    for (j in seq_len(J)) {
      width <- gmm_iv_mul(gmm_iv(zcrit),cert$se[,j]); ti <- c(cert$theta$lo[j],cert$theta$hi[j])
      z <- gmm_iv_sub(ti,width); clo[j]<-z[1L];chi[j]<-z[2L]
      z <- gmm_iv_add(ti,width); culo[j]<-z[1L];cuhi[j]<-z[2L]
    }
    cert$ci <- list(lo=clo,hi=chi); cert$ci_upper <- list(lo=culo,hi=cuhi); cert$pointwise_constant <- zcrit
  }
  bandlo <- est-crit*se; bandhi <- est+crit*se
  if (!is.null(cert)) cert <- ipw_pipeline_enforce_inference(cert,stats::setNames(est,names(estimates)),se,V,ia,bandlo,bandhi)
  z <- ifelse(se>0,est/se,NA_real_)
  list(estimates=stats::setNames(est,names(estimates)),se=se,analytical_se=analytical_se,
    z_stat=stats::setNames(z,names(estimates)),p_value=stats::setNames(2*stats::pnorm(-abs(z)),names(estimates)),
    ci_lower=stats::setNames(est-zcrit*se,names(estimates)),ci_upper=stats::setNames(est+zcrit*se,names(estimates)),
    cband_crit_val=crit,cband_lower=stats::setNames(bandlo,names(estimates)),cband_upper=stats::setNames(bandhi,names(estimates)),
    vcov=V,Omega=V*n,cluster_if=S,raw_cluster_sums=raw_cluster_sums,n_total=n,n_clusters=G,
    n_zero_eligible_clusters=sum(!seq_len(G) %in% unique(map[eligible_full])),
    clustered=!is.null(cluster_full),multiplier=mult,rng=rng,numerical_certificate=cert,
    # Effective inference record (separate from any requested parameter fields):
    # `bootstrap` is the reported-SE law, `simultaneous` the critical-value law;
    # whether draws ran is `!is.null(multiplier)`. They are not interchangeable:
    # under the released aggte-bin FALSE/TRUE law draws exist with analytic SEs.
    bootstrap=!is.null(mult) && flags$bootstrap,simultaneous=!is.null(mult) && flags$simultaneous,
    flag_law=flags$law,band_family=band_family,n_draws=if (is.null(mult)) NA_integer_ else as.integer(biters))
}

# D5 supported-subvector inference. Estimates, IF columns and enclosures are
# subset to the available columns BEFORE any transport, product, centering,
# scale inversion, multiplier or maximum; results are scattered into explicit
# NA-initialized positions. An empty supported selection returns a typed
# unavailable object with no quantile computation and no affirmative certificate.
ipw_pipeline_empty_inference <- function(ids, n, cluster_full, eligible_full, reason) {
  J <- length(ids); na_v <- stats::setNames(rep(NA_real_, J), ids)
  u_cl <- if (is.null(cluster_full)) seq_len(n) else unique(cluster_full)
  map <- if (is.null(cluster_full)) seq_len(n) else match(cluster_full, u_cl)
  if (is.null(eligible_full)) eligible_full <- rep(TRUE, n)
  list(status = "empty_supported_selection", reason = reason,
    estimates = na_v, se = na_v, analytical_se = na_v, z_stat = na_v, p_value = na_v, ci_lower = na_v, ci_upper = na_v,
    cband_crit_val = NA_real_, cband_lower = na_v, cband_upper = na_v,
    vcov = matrix(NA_real_, J, J, dimnames = list(ids, ids)), Omega = matrix(NA_real_, J, J, dimnames = list(ids, ids)),
    cluster_if = matrix(NA_real_, length(u_cl), J), raw_cluster_sums = matrix(NA_real_, length(u_cl), J),
    n_total = n, n_clusters = length(u_cl), n_zero_eligible_clusters = sum(!seq_along(u_cl) %in% unique(map[eligible_full])),
    clustered = !is.null(cluster_full), multiplier = NULL, rng = NULL,
    bootstrap = FALSE, simultaneous = FALSE, flag_law = ipw_pipeline_flag_law("empty"), band_family = "none", n_draws = NA_integer_,
    support = list(status = "empty_supported_selection", available_ids = character(0L), unavailable_ids = ids,
      n_available = 0L, n_requested = J, whole_request_complete = FALSE),
    numerical_certificate = list(certified = FALSE, status = "no_supported_targets", certified_target_ids = character(0L),
      unavailable_target_ids = ids, band_family = "none", whole_request_complete = FALSE,
      scope = "No supported target: no quantile, scale inversion or band was computed; nothing is certified."))
}

ipw_pipeline_supported_inference <- function(estimates, if_full, available, cluster_full = NULL, alp = 0.05,
  cband = FALSE, biters = 1000L, btype = "rademacher", seed = NULL, multiplier_budget, inference_accuracy,
  enclosure = NULL, fixed_signs = NULL, eligible_full = NULL, reason = "no requested target is supported",
  bstrap = NULL, caller_rng = FALSE, endpoint = c("cells", "aggte_bins", "simple"), target_meta = NULL) {
  endpoint <- match.arg(endpoint)
  est <- as.double(estimates); Psi <- as.matrix(if_full); J <- length(est); n <- nrow(Psi)
  ids <- names(estimates); if (is.null(ids)) ids <- paste0("t", seq_len(J))
  # The flag combination is validated before any selection so that an invalid
  # flag never depends on the support pattern.
  ipw_pipeline_resolve_bootstrap_flags(bstrap, cband, endpoint = endpoint)
  if (!is.logical(available) || length(available) != J || anyNA(available)) {
    cdid_ipw_pipeline_abort("`available` must be a complete logical selector over the requested targets.", "cdid_ipw_pipeline_input_error")
  }
  if (any(available & (!is.finite(est) | colSums(!is.finite(Psi)) > 0L))) {
    cdid_ipw_pipeline_abort("An available target must have a finite estimate and finite influence-function column.", "cdid_ipw_pipeline_integrity_error")
  }
  av <- which(available)
  if (!length(av)) return(ipw_pipeline_empty_inference(ids, n, cluster_full, eligible_full, reason))
  sub_enc <- if (is.null(enclosure)) NULL else list(
    theta = list(lo = enclosure$theta$lo[, av, drop = FALSE], hi = enclosure$theta$hi[, av, drop = FALSE]),
    rows = list(lo = enclosure$rows$lo[, av, drop = FALSE], hi = enclosure$rows$hi[, av, drop = FALSE]))
  if (!is.null(sub_enc) && (anyNA(sub_enc$theta$lo) || anyNA(sub_enc$theta$hi) || anyNA(sub_enc$rows$lo) || anyNA(sub_enc$rows$hi))) {
    cdid_ipw_pipeline_abort("Supported targets lack finite numerical enclosures.", "cdid_ipw_pipeline_integrity_error")
  }
  if (!is.null(target_meta) && (!is.data.frame(target_meta) || nrow(target_meta) != J))
    cdid_ipw_pipeline_abort("Reference metadata must cover the declared target rows before support selection.", "cdid_ipw_pipeline_input_error")
  sub_meta <- if (is.null(target_meta)) NULL else target_meta[av, , drop = FALSE]
  sub <- compute_chained_ipw_inference(stats::setNames(est[av], ids[av]), Psi[, av, drop = FALSE], cluster_full, alp, cband,
    biters, btype, seed, multiplier_budget, inference_accuracy, sub_enc, fixed_signs, eligible_full,
    bstrap = bstrap, caller_rng = caller_rng, endpoint = endpoint, target_meta = sub_meta)
  scatter_v <- function(v) { z <- stats::setNames(rep(NA_real_, J), ids); z[av] <- as.double(v); z }
  scatter_m <- function(m) { z <- matrix(NA_real_, J, J, dimnames = list(ids, ids)); z[av, av] <- m; z }
  scatter_cols <- function(m) { z <- matrix(NA_real_, nrow(m), J); z[, av] <- m; z }
  cert <- sub$numerical_certificate
  if (!is.null(cert)) {
    sc <- cert
    sc$theta <- list(lo = scatter_cols(cert$theta$lo), hi = scatter_cols(cert$theta$hi))
    sc$rows <- list(lo = scatter_cols(cert$rows$lo), hi = scatter_cols(cert$rows$hi))
    sc$cluster <- list(lo = scatter_cols(cert$cluster$lo), hi = scatter_cols(cert$cluster$hi))
    sc$vcov <- list(lo = scatter_m(cert$vcov$lo), hi = scatter_m(cert$vcov$hi))
    sc$se <- scatter_cols(cert$se)
    if (!is.null(cert$multiplier_scale)) sc$multiplier_scale <- scatter_cols(cert$multiplier_scale)
    sc$band <- list(lo = scatter_cols(cert$band$lo), hi = scatter_cols(cert$band$hi))
    sc$band_upper <- list(lo = scatter_cols(cert$band_upper$lo), hi = scatter_cols(cert$band_upper$hi))
    if (!is.null(cert$ci)) {
      sc$ci <- list(lo = scatter_cols(cert$ci$lo), hi = scatter_cols(cert$ci$hi))
      sc$ci_upper <- list(lo = scatter_cols(cert$ci_upper$lo), hi = scatter_cols(cert$ci_upper$hi))
    }
    if (!is.null(cert$analytic_se)) sc$analytic_se <- scatter_cols(cert$analytic_se)
    if (!is.null(cert$multiplier_error_by_target)) sc$multiplier_error_by_target <- scatter_v(cert$multiplier_error_by_target)
    if (!is.null(cert$raw_draw_scale_lower)) sc$raw_draw_scale_lower <- scatter_v(cert$raw_draw_scale_lower)
    if (!is.null(cert$deterministic)) sc$deterministic <- av[cert$deterministic]
    if (!is.null(cert$structural_reference)) sc$structural_reference <- av[cert$structural_reference]
    sc$certified_target_ids <- ids[av]
    sc$unavailable_target_ids <- ids[-av]
    sc$band_family <- sub$band_family
    sc$whole_request_complete <- length(av) == J
    sc$scope <- paste0(cert$scope, " Certified IDs: exactly the listed supported targets; unavailable IDs carry support status and no certificate; whole-request completeness: ", length(av) == J, ".")
    cert <- sc
  }
  list(status = if (length(av) == J) "all_requested_supported" else "supported_subvector",
    estimates = scatter_v(sub$estimates), se = scatter_v(sub$se), analytical_se = scatter_v(sub$analytical_se),
    z_stat = scatter_v(sub$z_stat), p_value = scatter_v(sub$p_value), ci_lower = scatter_v(sub$ci_lower), ci_upper = scatter_v(sub$ci_upper),
    cband_crit_val = sub$cband_crit_val, cband_lower = scatter_v(sub$cband_lower), cband_upper = scatter_v(sub$cband_upper),
    vcov = scatter_m(sub$vcov), Omega = scatter_m(sub$Omega), cluster_if = scatter_cols(sub$cluster_if),
    raw_cluster_sums = scatter_cols(sub$raw_cluster_sums), n_total = sub$n_total, n_clusters = sub$n_clusters,
    n_zero_eligible_clusters = sub$n_zero_eligible_clusters, clustered = sub$clustered, multiplier = sub$multiplier, rng = sub$rng,
    bootstrap = sub$bootstrap, simultaneous = sub$simultaneous, flag_law = sub$flag_law, band_family = sub$band_family, n_draws = sub$n_draws,
    support = list(status = if (length(av) == J) "all_requested_supported" else "supported_subvector", available_ids = ids[av],
      unavailable_ids = ids[-av], n_available = length(av), n_requested = J, whole_request_complete = length(av) == J),
    numerical_certificate = cert, compact = sub)
}

# ------------------------------------------------------------------------------
# 7. Population Cohort-Share Aggregation (`thm:agg` in `cp02_population_aggregation_if.tex`)
# ------------------------------------------------------------------------------
# D5: intended cell sets B are formed from the FULL declared ledger. With
# `na.rm = FALSE` an aggregate needing an unavailable cell is NA with the
# missing IDs and reasons; unrelated aggregates survive. With `na.rm = TRUE`
# the restricted set B* = B minus unavailable cells is used explicitly with
# renormalized weights (sec:maps ratio map, prop:cells(b): a changed target,
# displayed); empty B* stays unavailable. Outer summaries follow the released
# did 2.5.1 convention: under `na.rm = TRUE` empty bins are omitted and the
# surviving outer summary renormalized; under FALSE every intended bin is required.
# Arithmetic uses only supported finite columns (no zero-times-NA products).
# Unknown frame cohorts make share-dependent targets unavailable with a
# missing-share explanation (prop:cells(c)); no share is imputed.

ipw_pipeline_unavailable_aggregation <- function(typ, keys, labels, sets, target_meta, n, reason, na_rm,
                                                 cluster_full, eligible_full, shares_known) {
  M <- length(labels); ids <- labels
  na_v <- stats::setNames(rep(NA_real_, M), ids)
  inf <- ipw_pipeline_empty_inference(ids, n, cluster_full, eligible_full, reason)
  structure(list(type = typ, na_rm = na_rm, status = "unavailable", reason = reason, cell_keys = keys, cell_labels = labels,
    cell_target_ids = lapply(sets, function(rr) target_meta$target_id[rr]),
    cell_status = data.frame(label = labels, available = FALSE, restricted = FALSE, reason = reason, stringsAsFactors = FALSE),
    estimates = na_v, se = na_v, se_fixed_shares_only = na_v, ci_lower = na_v, ci_upper = na_v, cband_crit_val = NA_real_,
    cband_lower = na_v, cband_upper = na_v, if_full = matrix(NA_real_, n, M, dimnames = list(NULL, ids)),
    if_fixed_full = matrix(NA_real_, n, M), if_share_full = matrix(NA_real_, n, M), weights_theta = NULL, gradient_shares = NULL,
    cohort_shares = NULL, shares_known = shares_known, cohort_share_enclosure = NULL, if_shares_full = NULL, if_shares_enclosure = NULL,
    fixed_share_numerical_certificate = NULL, decomposition_enclosure = NULL, decomposition = NULL,
    inference = inf, numerical_certificate = inf$numerical_certificate,
    overall = list(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, if_full = rep(NA_real_, n),
      status = "unavailable", reason = reason, numerical_certificate = inf$numerical_certificate)),
    class = "cdid_ipw_aggregation_result")
}

# Released did 2.5.1 aggregation selectors (`balance_e`, `min_e`, `max_e`),
# bound to the chained ledger from the finite selector oracle
# (`plan/evidence/cp02_selector_reference`, compute.aggte.R l. 229-255, 279,
# 333-350, 448-467, 575-581). The selectors only choose cells among the
# reported ledger; they never change the estimation catalog. The intended
# cell/bin/cohort sets are formed from the FULL declared ledger before D5's
# missing-cell handling (`na.rm`).
#   simple/group : keep cells with g <= t and pos(t) - pos(g) <= max_e where pos
#                  is the position on the sorted union of calendar and cohorts
#                  (recoded grid; `min_e`/`balance_e` have no effect, no warning);
#   dynamic      : e = t - g in ORIGINAL date units; with balance_e, keep cohorts
#                  with max(t) - g >= balance_e and horizons in
#                  [balance_e - max(t) + min(grid), balance_e]; then [min_e, max_e];
#                  an empty window is refused like upstream;
#   calendar     : selectors are ignored with the released warning.
ipw_pipeline_aggregation_selection <- function(typ, g, tt, balance_e = NULL, min_e = -Inf, max_e = Inf,
                                               calendar = NULL, cohort_grid = NULL) {
  if (!is.numeric(min_e) || length(min_e) != 1L || is.na(min_e) || !is.numeric(max_e) || length(max_e) != 1L || is.na(max_e)) {
    cdid_ipw_pipeline_abort("`min_e` and `max_e` must be non-missing numeric scalars (possibly infinite).",
      c("cdid_ipw_pipeline_selector_error", "cdid_ipw_pipeline_input_error"))
  }
  if (!is.null(balance_e) && (!is.numeric(balance_e) || length(balance_e) != 1L || is.na(balance_e) || !is.finite(balance_e) ||
      balance_e < 0 || balance_e != floor(balance_e))) {
    cdid_ipw_pipeline_abort("`balance_e` must be NULL or a nonnegative whole number.",
      c("cdid_ipw_pipeline_selector_error", "cdid_ipw_pipeline_input_error"))
  }
  min_e <- as.double(min_e); max_e <- as.double(max_e)
  if (!is.null(balance_e)) balance_e <- as.double(balance_e)
  nondefault <- !is.null(balance_e) || is.finite(min_e) || is.finite(max_e)
  grid <- sort(unique(c(as.double(calendar), as.double(cohort_grid), g, tt)))
  grid <- grid[is.finite(grid)]
  pos <- function(x) match(as.double(x), grid)
  R <- length(g)
  post <- which(tt >= g)
  rec <- list(type = typ, balance_e = balance_e, min_e = min_e, max_e = max_e, nondefault = nondefault,
    applied = nondefault && typ != "calendar", ignored = nondefault && typ == "calendar", warning = NULL,
    recoded_grid = grid, convention = switch(typ,
      simple = "released did 2.5.1: cells with g <= t and recoded-grid horizon pos(t) - pos(g) <= max_e; min_e and balance_e have no effect",
      group = "released did 2.5.1: per cohort, cells with g <= t and recoded-grid horizon pos(t) - pos(g) <= max_e; min_e and balance_e have no effect",
      dynamic = "released did 2.5.1: e = t - g in original date units; balance_e keeps cohorts with max(t) - g >= balance_e and horizons in [balance_e - max(t) + min(grid), balance_e]; then [min_e, max_e]",
      calendar = "released did 2.5.1: min_e, max_e and balance_e are ignored for calendar aggregation (warning)"),
    scope = "selectors choose among already reported ledger cells; the estimation catalog, moments, weights and IFs are unchanged")
  if (typ == "calendar") {
    if (nondefault) rec$warning <- paste0("`min_e`, `max_e`, and `balance_e` are ignored for type = \"calendar\"; ",
                                          "returning the unrestricted calendar-time effects.")
    keys <- sort(unique(tt[post])); sets <- lapply(keys, function(h) post[tt[post] == h])
    rec$selected_cells <- post
    rec$eligible_cohorts <- sort(unique(g[post]))
  } else if (typ %in% c("simple", "group")) {
    horizon <- pos(tt) - pos(g)
    keep <- which(tt >= g & horizon <= max_e)
    rec$recoded_horizon <- horizon
    rec$selected_cells <- keep
    rec$eligible_cohorts <- sort(unique(g[post]))
    if (typ == "simple") { keys <- "simple"; sets <- list(keep) } else {
      keys <- as.character(sort(unique(g[post]))); sets <- lapply(sort(unique(g[post])), function(h) keep[g[keep] == h])
    }
  } else {
    e_all <- tt - g
    include <- rep(TRUE, R)
    max_t <- max(tt); min_grid <- min(grid)
    if (!is.null(balance_e)) {
      include <- (max_t - g >= balance_e)
      eseq <- sort(unique(e_all[include]))
      eseq <- eseq[eseq <= balance_e & eseq >= balance_e - max_t + min_grid]
    } else eseq <- sort(unique(e_all))
    eseq <- eseq[eseq >= min_e & eseq <= max_e]
    if (!length(eseq)) {
      cdid_ipw_pipeline_abort(paste0("No event times fall within the requested window. Adjust 'min_e'/'max_e' (and 'balance_e') ",
        "so at least one event time is included (released did 2.5.1 refusal; the reported cells are unchanged)."),
        c("cdid_ipw_pipeline_selector_error", "cdid_ipw_pipeline_input_error"),
        list(balance_e = balance_e, min_e = min_e, max_e = max_e, declared_event_times = sort(unique(e_all))))
    }
    if (nondefault && !any(eseq >= 0)) {
      cdid_ipw_pipeline_abort("No valid att_gt() estimates found for this aggregation: the explicit event window contains no nonnegative horizon (released did 2.5.1 selector refusal).",
        c("cdid_ipw_pipeline_selector_error", "cdid_ipw_pipeline_input_error"),
        list(balance_e = balance_e, min_e = min_e, max_e = max_e, selected_event_times = eseq))
    }
    keys <- eseq; sets <- lapply(eseq, function(e) which(e_all == e & include))
    rec$include_cell <- include
    rec$eligible_cohorts <- sort(unique(g[include]))
    rec$excluded_cohorts <- sort(unique(g[!include]))
    rec$event_window <- eseq
    rec$declared_event_times <- sort(unique(e_all))
    rec$selected_cells <- sort(unique(unlist(sets)))
  }
  list(keys = keys, sets = sets, record = rec)
}

aggregate_chained_ipw_targets <- function(
  theta_hat, target_if_full, target_meta, full_frame_g, full_frame_w,
  type = c("simple", "group", "dynamic", "calendar"), anticipation = 0L,
  cluster_full = NULL, alp = 0.05, cband = FALSE, biters = 1000L,
  btype = "rademacher", seed = NULL, multiplier_budget, inference_accuracy,
  enclosure = NULL, fixed_signs = NULL, eligible_full = NULL, na.rm = FALSE, available = NULL,
  bstrap = NULL, caller_rng = FALSE, balance_e = NULL, min_e = -Inf, max_e = Inf,
  calendar = NULL, cohort_grid = NULL, share_law = NULL, frame_ids = NULL
) {
  typ <- match.arg(type); th <- as.double(theta_hat); X <- as.matrix(target_if_full)
  n <- nrow(X); R <- ncol(X); g <- as.double(target_meta$g); tt <- as.double(target_meta$time)
  na_rm <- na.rm
  if (!gmm_scalar_logical(na_rm)) cdid_ipw_pipeline_abort("`na.rm` must be a scalar logical.", "cdid_ipw_pipeline_input_error")
  # Released endpoint of the bin family: group/dynamic/calendar bins may carry a
  # multiplier band under bstrap = FALSE, cband = TRUE (compute.aggte.R:365-370,
  # 509-514, 642-647); simple aggregation has no bin-band family (l. 306-315).
  bin_endpoint <- if (identical(typ, "simple")) "simple" else "aggte_bins"
  ipw_pipeline_resolve_bootstrap_flags(bstrap, cband, where = "aggregation", endpoint = bin_endpoint)
  # D61: under first-calendar masses V = w_t1 may be absent (NA) on units whose
  # mass no requested share reads; present values must still be positive.
  first_masses <- identical(share_law, "first_calendar_masses")
  w_bad <- if (first_masses) any(!is.na(full_frame_w) & (!is.finite(full_frame_w) | full_frame_w <= 0)) else
    any(!is.finite(full_frame_w)) || any(full_frame_w <= 0)
  if (length(th) != R || length(full_frame_g) != n || length(full_frame_w) != n ||
      w_bad || any(is.nan(full_frame_g) | full_frame_g == -Inf, na.rm = TRUE) ||
      (!is.null(frame_ids) && length(frame_ids) != n)) {
    cdid_ipw_pipeline_abort("Aggregation needs same-frame cohort labels (NA only for genuinely unknown frame-only cohorts) and positive fixed weights.", "cdid_ipw_pipeline_aggregation_error")
  }
  incomplete_masses <- first_masses && anyNA(full_frame_w)
  if (incomplete_masses && !is.null(enclosure))
    cdid_ipw_pipeline_abort("Interval enclosures are not available with incomplete first-calendar masses.", "cdid_ipw_pipeline_aggregation_error")
  if (is.null(available)) available <- is.finite(th) & colSums(!is.finite(X)) == 0L
  if (!is.logical(available) || length(available) != R || anyNA(available) ||
      any(available & (!is.finite(th) | colSums(!is.finite(X)) > 0L))) {
    cdid_ipw_pipeline_abort("`available` must be a complete logical selector; available cells need finite estimates and IF columns.", "cdid_ipw_pipeline_input_error")
  }
  ids_all <- target_meta$target_id
  post <- which(tt >= g); cohorts <- sort(unique(g[post])); K <- length(cohorts)
  shares_known <- !anyNA(full_frame_g)
  # Intended sets from the full declared ledger under the released selector
  # conventions, before any D5 missing-cell handling.
  sel <- ipw_pipeline_aggregation_selection(typ, g, tt, balance_e, min_e, max_e, calendar, cohort_grid)
  if (!is.null(sel$record$warning)) warning(sel$record$warning, call. = FALSE)
  keys <- sel$keys; sets <- sel$sets
  labels <- switch(typ, simple = "ATT_simple", group = paste0("group_", keys), calendar = paste0("calendar_", keys), paste0("event_", keys))
  selectors <- sel$record
  selectors$selected_cell_ids <- ids_all[selectors$selected_cells]
  if (!K) {
    out <- ipw_pipeline_unavailable_aggregation(typ, keys, labels, sets, target_meta, n, "no_declared_post_cells: the reporting ledger declares no post cell",
      na_rm, cluster_full, eligible_full, shares_known)
    out$selectors <- selectors
    return(out)
  }
  M <- length(sets)
  # Full-frame weighted cohort shares and their ratio IF (lem:transport, thm:agg).
  p <- rep(NA_real_, K); P <- matrix(0.0, n, K); piv <- if (is.null(enclosure)) NULL else matrix(NA_real_, 2L, K); Penc <- NULL
  names(p) <- colnames(P) <- as.character(cohorts)
  mass_complete <- stats::setNames(rep(TRUE, K), as.character(cohorts)); mass_missing <- character(K)
  if (shares_known && incomplete_masses) {
    # D61 (TeX (5)-(6)): masses m_g = H(I_g V) with IF u_g - m_g. Every share
    # map here is degree 0, so the unnormalized masses feed the same weights,
    # derivatives and cross-covariances; never-treated V is not demanded and a
    # cohort lacking V is refused only by the aggregates that read its mass.
    for (k in seq_len(K)) {
      ind <- full_frame_g == cohorts[k]
      gap <- which(ind & is.na(full_frame_w))
      if (length(gap)) {
        mass_complete[k] <- FALSE
        mass_missing[k] <- paste0("cohort ", fmt_design_num(cohorts[k]), " lacks w at the first calendar date for ",
          length(gap), " unit(s)", if (!is.null(frame_ids)) paste0(" (ids ", paste(utils::head(frame_ids[gap], 5L),
            collapse = ", "), if (length(gap) > 5L) ", ..." else "", ")"))
        next
      }
      u <- ifelse(ind, full_frame_w, 0.0)
      p[k] <- sum(u) / n
      P[, k] <- u - p[k]
    }
  } else if (shares_known) {
    if (!is.null(enclosure)) {
      totalw <- gmm_iv(0.0)
      for (i in seq_len(n)) totalw <- gmm_iv_add(totalw, gmm_iv(as.double(full_frame_w[i])))
    }
    meanw <- mean(full_frame_w)
    for (k in seq_len(K)) {
      ind <- as.double(full_frame_g == cohorts[k]); p[k] <- sum(full_frame_w * ind) / sum(full_frame_w)
      P[, k] <- full_frame_w * (ind - p[k]) / meanw
      if (!is.null(enclosure)) {
        z <- gmm_iv(0.0)
        for (i in seq_len(n)) if (ind[i] == 1.0) z <- gmm_iv_add(z, gmm_iv(as.double(full_frame_w[i])))
        piv[, k] <- gmm_iv_div(z, totalw)
      }
    }
    if (!is.null(enclosure)) {
      plo <- phi <- matrix(0.0, n, K)
      mw <- gmm_iv_div(totalw, gmm_iv(as.double(n)))
      for (i in seq_len(n)) for (k in seq_len(K)) {
        z <- gmm_iv_div(gmm_iv_mul(gmm_iv(as.double(full_frame_w[i])), gmm_iv_sub(gmm_iv(as.double(full_frame_g[i] == cohorts[k])), piv[, k])), mw)
        plo[i, k] <- z[1L]; phi[i, k] <- z[2L]
      }
      Penc <- list(lo = plo, hi = phi)
    }
  }
  n_unknown <- sum(is.na(full_frame_g))
  share_reason <- paste0("missing_frame_cohort_share: ", n_unknown, " frame unit(s) have unknown cohort, so population shares of cohorts {",
    paste(fmt_design_num(cohorts), collapse = ","), "} are not identified by the supplied frame; no share is inferred from eligible records")
  # Per-aggregate intended sets, restricted sets and local weights.
  W <- matrix(0.0, M, R); D <- matrix(0.0, M, K); estimates <- rep(NA_real_, M)
  Wlo <- Whi <- if (is.null(enclosure)) NULL else W
  Dlo <- Dhi <- if (is.null(enclosure)) NULL else D
  tlo <- thi <- if (is.null(enclosure)) NULL else matrix(NA_real_, 1L, M)
  cell_available <- logical(M); restricted <- logical(M); share_dep <- logical(M)
  share_denominator <- rep(NA_real_, M)
  reasons <- rep(NA_character_, M); used_sets <- vector("list", M); excluded <- vector("list", M)
  for (m in seq_len(M)) {
    rr <- sets[[m]]; miss <- rr[!available[rr]]
    excluded[[m]] <- ids_all[miss]
    if (!length(rr)) { reasons[m] <- "empty_intended_set"; next }
    if (!na_rm && length(miss)) {
      reasons[m] <- paste0("missing_targets (na.rm = FALSE): ", paste(paste0(ids_all[miss], " [", target_meta$support_status[miss], "]"), collapse = "; "))
      next
    }
    used <- rr[available[rr]]
    if (!length(used)) { reasons[m] <- paste0("empty_restricted_set (na.rm = TRUE): every intended cell is unavailable: ", paste(ids_all[miss], collapse = "; ")); next }
    restricted[m] <- length(miss) > 0L
    used_sets[[m]] <- used
    kk <- match(g[used], cohorts); n_coh <- length(unique(kk))
    share_dep[m] <- (typ != "group") && n_coh >= 2L
    if (share_dep[m] && !shares_known) { reasons[m] <- share_reason; next }
    if (share_dep[m] && !all(mass_complete[kk])) {
      reasons[m] <- paste0("missing_first_period_weight (D61): ", paste(unique(mass_missing[kk][!mass_complete[kk]]), collapse = "; "),
        "; first-calendar masses are never replaced by first-observed, mean or anchor weights")
      next
    }
    if (typ != "group" && shares_known && (!incomplete_masses || share_dep[m])) {
      share_denominator[m] <- sum(p[kk])
      if (!is.finite(share_denominator[m]) || share_denominator[m] <= 0) { reasons[m] <- "nonpositive_share_mass"; next }
    }
    if (typ == "group" || n_coh == 1L) {
      # Fixed uniform weights; the single-cohort ratio map is degree-0 with exact zero share derivative.
      w <- rep(1.0 / length(used), length(used)); den <- NA_real_
    } else {
      den <- share_denominator[m]
      w <- p[kk] / den
    }
    estimates[m] <- sum(w * th[used]); W[m, used] <- w; cell_available[m] <- TRUE
    if (share_dep[m]) for (j in seq_along(used)) D[m, kk[j]] <- D[m, kk[j]] + (th[used[j]] - estimates[m]) / den
    if (!is.null(enclosure)) {
      ti <- gmm_iv(0.0)
      if (share_dep[m]) {
        di <- gmm_iv(0.0)
        for (k in kk) di <- gmm_iv_add(di, piv[, k])
        if (di[1L] <= 0) cdid_ipw_pipeline_abort("Aggregate mass lacks a positive enclosure.", "cdid_ipw_pipeline_transfer_inconclusive")
      }
      for (j in seq_along(used)) {
        wi <- if (share_dep[m]) gmm_iv_div(piv[, kk[j]], di) else gmm_iv_div(gmm_iv(1.0), gmm_iv(as.double(length(used))))
        Wlo[m, used[j]] <- wi[1L]; Whi[m, used[j]] <- wi[2L]
        ti <- gmm_iv_add(ti, gmm_iv_mul(wi, c(enclosure$theta$lo[used[j]], enclosure$theta$hi[used[j]])))
      }
      tlo[m] <- ti[1L]; thi[m] <- ti[2L]
      if (share_dep[m]) for (j in seq_along(used)) {
        z <- gmm_iv_div(gmm_iv_sub(c(enclosure$theta$lo[used[j]], enclosure$theta$hi[used[j]]), ti), di)
        z <- gmm_iv_add(c(Dlo[m, kk[j]], Dhi[m, kk[j]]), z)
        Dlo[m, kk[j]] <- z[1L]; Dhi[m, kk[j]] <- z[2L]
      }
    }
  }
  names(estimates) <- labels; rownames(W) <- rownames(D) <- labels; colnames(W) <- ids_all; colnames(D) <- as.character(cohorts)
  av <- which(available); am <- which(cell_available)
  # Local products on supported columns only; unavailable aggregates are NA columns.
  fixed <- share <- total <- matrix(NA_real_, n, M); colnames(fixed) <- colnames(share) <- colnames(total) <- labels
  if (length(am)) {
    fixed[, am] <- X[, av, drop = FALSE] %*% t(W[am, av, drop = FALSE])
    share[, am] <- if (shares_known) P %*% t(D[am, , drop = FALSE]) else 0.0
    total[, am] <- fixed[, am, drop = FALSE] + share[, am, drop = FALSE]
  }
  make_enclosure <- function(wlo, whi, dlo, dhi, tl, tu, rows_sel) {
    # rows_sel: aggregate rows with finite enclosures; wlo/whi are rows x R (zero outside av).
    if (is.null(enclosure) || !length(rows_sel)) return(NULL)
    a <- ipw_pipeline_interval_product(enclosure$rows$lo[, av, drop = FALSE], enclosure$rows$hi[, av, drop = FALSE],
      t(wlo[rows_sel, av, drop = FALSE]), t(whi[rows_sel, av, drop = FALSE]))
    fixed_box <- a
    if (shares_known) {
      b <- ipw_pipeline_interval_product(Penc$lo, Penc$hi, t(dlo[rows_sel, , drop = FALSE]), t(dhi[rows_sel, , drop = FALSE]))
      for (k in seq_along(a$lo)) { z <- gmm_iv_add(c(a$lo[k], a$hi[k]), c(b$lo[k], b$hi[k])); a$lo[k] <- z[1L]; a$hi[k] <- z[2L] }
    } else b <- list(lo = matrix(0.0, n, length(rows_sel)), hi = matrix(0.0, n, length(rows_sel)))
    list(theta = list(lo = tl[, rows_sel, drop = FALSE], hi = tu[, rows_sel, drop = FALSE]), rows = a, fixed_rows = fixed_box, share_rows = b)
  }
  scatter_enc <- function(e, rows_sel) {
    if (is.null(e)) return(NULL)
    full <- function(z) { o <- matrix(NA_real_, n, M); o[, rows_sel] <- z; o }
    th_ <- function(z) { o <- matrix(NA_real_, 1L, M); o[, rows_sel] <- z; o }
    list(theta = list(lo = th_(e$theta$lo), hi = th_(e$theta$hi)), rows = list(lo = full(e$rows$lo), hi = full(e$rows$hi)),
      fixed_rows = list(lo = full(e$fixed_rows$lo), hi = full(e$fixed_rows$hi)), share_rows = list(lo = full(e$share_rows$lo), hi = full(e$share_rows$hi)))
  }
  enc_c <- make_enclosure(Wlo, Whi, Dlo, Dhi, tlo, thi, am)
  bin_meta <- NULL
  if (ipw_pipeline_accuracy_is_relative(inference_accuracy) && !is.null(enclosure)) {
    src_reference <- ipw_pipeline_known_reference_rows(target_meta) &
      as.double(enclosure$theta$lo) == 0 & as.double(enclosure$theta$hi) == 0 &
      colSums(enclosure$rows$lo != 0 | enclosure$rows$hi != 0, na.rm = FALSE) == 0L
    src_reference[is.na(src_reference)] <- FALSE
    ref_bins <- vapply(used_sets, function(used) length(used) > 0L && all(src_reference[used]), logical(1L)) & cell_available
    bin_meta <- data.frame(target_id = labels, role = ifelse(ref_bins, "reference_aggregate", "aggregate"),
      support_status = ifelse(ref_bins, "deterministic_reference_zero", "aggregate"),
      time = NA_real_, ref = NA_real_, base_t = NA_real_,
      reference_operator = ifelse(ref_bins, "weighted_reference_identity", "none"), stringsAsFactors = FALSE)
    bin_meta$reference_sources <- lapply(seq_len(M), function(m) if (ref_bins[m]) target_meta[used_sets[[m]], , drop = FALSE] else NULL)
    # This is an exact operator identity, not tightening from zero containment:
    # every source estimate and IF enclosure is the singleton zero, so both the
    # weighted effect and its share derivative vanish for all admissible shares.
    if (!is.null(enc_c)) for (h in which(ref_bins[am])) {
      enc_c$theta$lo[h] <- enc_c$theta$hi[h] <- 0.0
      for (family in c("rows", "fixed_rows", "share_rows")) {
        enc_c[[family]]$lo[, h] <- enc_c[[family]]$hi[, h] <- 0.0
      }
    }
  }
  enc <- scatter_enc(enc_c, am)
  cell_reason <- if (length(am) < M) paste0("unavailable aggregates: ", paste(labels[!cell_available], collapse = ", ")) else "all intended aggregates available"
  inf <- ipw_pipeline_supported_inference(estimates, total, cell_available, cluster_full, alp, cband, biters, btype, seed,
    multiplier_budget, inference_accuracy, enc, fixed_signs, eligible_full, reason = cell_reason,
    bstrap = bstrap, caller_rng = caller_rng, endpoint = bin_endpoint, target_meta = bin_meta)
  # Released warning (compute.aggte.R:368/512/646) exactly when the band draws
  # ran under bstrap = FALSE, cband = TRUE; exact-zero/empty selections draw
  # nothing and therefore warn about nothing.
  if (!is.null(inf$multiplier) && isTRUE(inf$simultaneous) && !isTRUE(inf$bootstrap) &&
      identical(inf$flag_law, ipw_pipeline_flag_law("aggte_bins_false_true"))) {
    warning(ipw_pipeline_released_band_warning(), call. = FALSE)
  }
  # Fixed-share-only SEs are an analytic diagnostic (no draws are invented for them).
  fixed_inf <- ipw_pipeline_supported_inference(estimates, fixed, cell_available, cluster_full, alp, FALSE, biters, btype, seed,
    multiplier_budget, inference_accuracy, if (is.null(enc)) NULL else list(theta = enc$theta, rows = enc$fixed_rows), NULL, eligible_full, reason = cell_reason,
    bstrap = FALSE, caller_rng = FALSE, target_meta = bin_meta)
  decomposition_enclosure <- NULL
  if (!is.null(enc_c) && length(am)) {
    vf <- fixed_inf$compact$numerical_certificate$vcov
    vs <- ipw_pipeline_inference_intervals(list(theta = enc_c$theta, rows = enc_c$share_rows), cluster_full)$vcov
    vt <- inf$compact$numerical_certificate$vcov
    vc <- list(lo = vt$lo, hi = vt$hi)
    for (k in seq_along(vc$lo)) { v <- gmm_iv_sub(gmm_iv_sub(c(vt$lo[k], vt$hi[k]), c(vf$lo[k], vf$hi[k])), c(vs$lo[k], vs$hi[k])); vc$lo[k] <- v[1L]; vc$hi[k] <- v[2L] }
    scale <- function(v) { for (k in seq_along(v$lo)) { z <- gmm_iv_mul(c(v$lo[k], v$hi[k]), gmm_iv(as.double(n))); v$lo[k] <- z[1L]; v$hi[k] <- z[2L] }; v }
    sc_m <- function(v) { o <- list(lo = matrix(NA_real_, M, M), hi = matrix(NA_real_, M, M)); o$lo[am, am] <- v$lo; o$hi[am, am] <- v$hi; o }
    decomposition_enclosure <- list(Omega_total = sc_m(scale(vt)), Omega_fixed = sc_m(scale(vf)), Omega_share = sc_m(scale(vs)), Omega_cross = sc_m(scale(vc)))
    if (shares_known) {
      joint <- ipw_pipeline_inference_intervals(list(theta = list(lo = matrix(c(th[av], p), 1L), hi = matrix(c(th[av], p), 1L)),
        rows = list(lo = cbind(enclosure$rows$lo[, av, drop = FALSE], Penc$lo), hi = cbind(enclosure$rows$hi[, av, drop = FALSE], Penc$hi))), cluster_full)$vcov
      na_ <- length(av)
      otp_c <- scale(list(lo = joint$lo[seq_len(na_), na_ + seq_len(K), drop = FALSE], hi = joint$hi[seq_len(na_), na_ + seq_len(K), drop = FALSE]))
      otp <- list(lo = matrix(NA_real_, R, K), hi = matrix(NA_real_, R, K))
      otp$lo[av, ] <- otp_c$lo; otp$hi[av, ] <- otp_c$hi
      decomposition_enclosure$Omega_theta_p <- otp
    }
  }
  # Both IF blocks use the declared inference operator, preserving all cross-record terms.
  center_sums <- function(z) {
    map <- if (is.null(cluster_full)) seq_len(n) else match(cluster_full, unique(cluster_full))
    out <- matrix(0.0, max(map), ncol(z))
    if (is.null(cluster_full) && typeof(z) == "double" && !is.object(z)) {
      out[, ] <- 0.0 + z
    } else {
      for (i in seq_len(n)) out[map[i], ] <- out[map[i], ] + z[i, ]
    }
    sweep(out, 2L, colMeans(out), "-")
  }
  decomposition <- NULL
  if (length(am)) {
    xc <- center_sums(X[, av, drop = FALSE]); Ott <- crossprod(xc) / n
    Wa <- W[am, av, drop = FALSE]
    Of <- Wa %*% Ott %*% t(Wa)
    full_m <- function(z) { o <- matrix(NA_real_, M, M, dimnames = list(labels, labels)); o[am, am] <- z; o }
    if (shares_known) {
      pc <- center_sums(P); Otp <- crossprod(xc, pc) / n; Opp <- crossprod(pc) / n
      Da <- D[am, , drop = FALSE]
      Os <- Da %*% Opp %*% t(Da); Oc <- Wa %*% Otp %*% t(Da) + Da %*% t(Otp) %*% t(Wa)
      Otp_full <- matrix(NA_real_, R, K, dimnames = list(ids_all, as.character(cohorts))); Otp_full[av, ] <- Otp
    } else { Os <- Oc <- matrix(0.0, length(am), length(am)); Otp_full <- NULL }
    decomposition <- list(Omega_total = full_m(Of + Os + Oc), Omega_fixed = full_m(Of), Omega_share = full_m(Os), Omega_cross = full_m(Oc), Omega_theta_p = Otp_full)
  }
  # Outer summary: intended bins (na.rm = FALSE) or surviving bins renormalized (na.rm = TRUE; did 2.5.1 convention).
  outer <- if (typ == "dynamic") which(as.double(keys) >= 0) else seq_len(M)
  outer_missing <- outer[!cell_available[outer]]
  surv <- if (na_rm) outer[cell_available[outer]] else if (length(outer_missing)) integer(0L) else outer
  overall_status <- "available"; overall_reason <- NA_character_
  if (!length(outer)) { overall_status <- "unavailable"; overall_reason <- "no_outer_bins: the declared grid has no nonnegative outer bin" }
  else if (!na_rm && length(outer_missing)) { overall_status <- "unavailable"; overall_reason <- paste0("missing_inner_aggregates (na.rm = FALSE): ", paste(labels[outer_missing], collapse = ", ")) }
  else if (!length(surv)) { overall_status <- "unavailable"; overall_reason <- paste0("empty_restricted_outer_set (na.rm = TRUE): every intended bin is unavailable: ", paste(labels[outer], collapse = ", ")) }
  else if (typ == "group" && length(surv) >= 2L && !shares_known) { overall_status <- "unavailable"; overall_reason <- share_reason }
  else if (typ == "group" && length(surv) >= 2L && !all(mass_complete[surv])) {
    overall_status <- "unavailable"
    overall_reason <- paste0("missing_first_period_weight (D61): ", paste(mass_missing[surv][!mass_complete[surv]], collapse = "; "))
  }
  # A one-bin group overall reads no mass (exact degree-0 identity), so an
  # absent first-calendar mass of that cohort leaves it available.
  reads_mass <- !incomplete_masses || length(surv) >= 2L
  overall_share_denominator <- if (typ == "group" && shares_known && length(surv) && reads_mass) sum(p[surv]) else NA_real_
  if (overall_status == "available" && typ == "group" && shares_known && reads_mass &&
      (!is.finite(overall_share_denominator) || overall_share_denominator <= 0)) {
    overall_status <- "unavailable"; overall_reason <- "nonpositive_share_mass"
  }
  ow <- rep(0.0, M); og <- numeric(K); overall <- NA_real_; overall_if <- rep(NA_real_, n); oe <- NULL
  overall_share_dep <- FALSE
  if (overall_status == "available") {
    if (typ == "group") {
      if (length(surv) >= 2L) {
        S_g <- sum(p[surv]); ow[surv] <- p[surv] / S_g; overall <- sum(ow[surv] * estimates[surv])
        og[surv] <- (estimates[surv] - overall) / S_g; overall_share_dep <- TRUE
      } else { ow[surv] <- 1.0; overall <- estimates[surv] }
    } else { ow[surv] <- 1.0 / length(surv); overall <- sum(ow[surv] * estimates[surv]) }
    overall_if <- as.double(total[, surv, drop = FALSE] %*% ow[surv]) + if (overall_share_dep) as.double(P %*% og) else 0.0
    if (!is.null(enc_c)) {
      srows <- match(surv, am)
      if (typ == "group" && length(surv) >= 2L) {
        den <- gmm_iv(0.0); for (k in surv) den <- gmm_iv_add(den, piv[, k])
        owlo <- owhi <- matrix(0.0, 1L, length(am)); ot <- gmm_iv(0.0)
        for (h in seq_along(surv)) { z <- gmm_iv_div(piv[, surv[h]], den); owlo[srows[h]] <- z[1L]; owhi[srows[h]] <- z[2L]; ot <- gmm_iv_add(ot, gmm_iv_mul(z, c(tlo[surv[h]], thi[surv[h]]))) }
        ww <- ipw_pipeline_interval_product(owlo, owhi, Wlo[am, , drop = FALSE], Whi[am, , drop = FALSE])
        dd <- ipw_pipeline_interval_product(owlo, owhi, Dlo[am, , drop = FALSE], Dhi[am, , drop = FALSE])
        for (h in seq_along(surv)) { k <- surv[h]; z <- gmm_iv_div(gmm_iv_sub(c(tlo[k], thi[k]), ot), den); z <- gmm_iv_add(c(dd$lo[k], dd$hi[k]), z); dd$lo[k] <- z[1L]; dd$hi[k] <- z[2L] }
      } else {
        oi <- matrix(0.0, 1L, length(am)); oi[srows] <- ow[surv]
        ww <- ipw_pipeline_interval_product(oi, oi, Wlo[am, , drop = FALSE], Whi[am, , drop = FALSE])
        dd <- ipw_pipeline_interval_product(oi, oi, Dlo[am, , drop = FALSE], Dhi[am, , drop = FALSE])
        z <- ipw_pipeline_interval_product(tlo[, am, drop = FALSE], thi[, am, drop = FALSE], t(oi)); ot <- c(z$lo[1L], z$hi[1L])
      }
      oe <- make_enclosure(ww$lo, ww$hi, dd$lo, dd$hi, matrix(ot[1L], 1L), matrix(ot[2L], 1L), 1L)
    }
  }
  # Overall summary: pointwise interval by upstream convention; its SE is the
  # multiplier IQR SE when `bstrap` is explicitly TRUE (did getSE), analytic
  # otherwise (and under the existing private law, which bootstraps only bands).
  reuse_simple_inference <- identical(typ, "simple") && !isTRUE(inf$simultaneous)
  oi_res <- if (reuse_simple_inference) inf else
    ipw_pipeline_supported_inference(c(overall = overall), matrix(overall_if, ncol = 1L), overall_status == "available", cluster_full, alp, FALSE,
      biters, btype, seed, multiplier_budget, inference_accuracy, oe, fixed_signs, eligible_full, reason = if (is.na(overall_reason)) "overall unavailable" else overall_reason,
      bstrap = if (is.null(bstrap)) FALSE else bstrap, caller_rng = caller_rng)
  cell_status <- data.frame(label = labels, available = cell_available, restricted = restricted, share_dependent = share_dep,
    n_intended = lengths(sets), n_used = lengths(used_sets),
    share_denominator = share_denominator,
    intended_ids = vapply(sets, function(rr) paste(ids_all[rr], collapse = "; "), character(1L)),
    used_ids = vapply(used_sets, function(rr) if (is.null(rr)) "" else paste(ids_all[rr], collapse = "; "), character(1L)),
    excluded_ids = vapply(excluded, function(x) paste(x, collapse = "; "), character(1L)),
    reason = ifelse(cell_available, if (na_rm) ifelse(restricted, "restricted_target (na.rm = TRUE): weights renormalized over the supported subset", "intended_target") else "intended_target", reasons),
    stringsAsFactors = FALSE)
  target_definition <- if (na_rm) "explicit restricted targets: each aggregate averages its SUPPORTED intended cells with renormalized weights; outer summaries omit empty bins and renormalize over surviving bins (did 2.5.1 na.rm = TRUE convention); the intended sets, exclusions and changed targets are displayed in `cell_status` and `overall`" else
    "intended targets: each aggregate requires every intended cell; an aggregate with an unavailable intended cell is NA with the missing IDs and reasons; unrelated aggregates survive"
  out <- structure(list(type = typ, na_rm = na_rm, status = if (length(am) == M && overall_status == "available") "available" else if (length(am)) "partially_available" else "unavailable",
    cell_keys = keys, cell_labels = labels, cell_target_ids = lapply(sets, function(rr) ids_all[rr]),
    cell_used_target_ids = lapply(used_sets, function(rr) if (is.null(rr)) character(0L) else ids_all[rr]),
    cell_status = cell_status, target_definition = target_definition, shares_known = shares_known,
    estimates = estimates, se = inf$se, se_fixed_shares_only = fixed_inf$se, ci_lower = inf$ci_lower, ci_upper = inf$ci_upper,
    cband_crit_val = inf$cband_crit_val, cband_lower = inf$cband_lower, cband_upper = inf$cband_upper,
    if_full = total, if_fixed_full = fixed, if_share_full = share, weights_theta = W, gradient_shares = D,
    cohort_shares = p, cohort_share_enclosure = if (is.null(enc) || !shares_known) NULL else list(lo = piv[1L, ], hi = piv[2L, ]),
    if_shares_full = if (shares_known) P else NULL, if_shares_enclosure = Penc, fixed_share_numerical_certificate = fixed_inf$numerical_certificate,
    # Retained aggregation boxes (already computed; diagnostic output only):
    # outward weight boxes W (M x R), share-derivative boxes D (M x K) and the
    # aggregate point boxes. Rows of unavailable aggregates and columns of
    # unavailable cells are exact zeros / NA exactly as in `weights_theta`.
    weight_enclosure = if (is.null(enclosure)) NULL else list(lo = structure(Wlo, dimnames = dimnames(W)), hi = structure(Whi, dimnames = dimnames(W))),
    share_gradient_enclosure = if (is.null(enclosure)) NULL else list(lo = structure(Dlo, dimnames = dimnames(D)), hi = structure(Dhi, dimnames = dimnames(D))),
    cell_estimate_enclosure = if (is.null(enclosure)) NULL else list(lo = stats::setNames(as.double(tlo), labels), hi = stats::setNames(as.double(thi), labels)),
    decomposition_enclosure = decomposition_enclosure, decomposition = decomposition,
    inference = inf, numerical_certificate = inf$numerical_certificate,
    selectors = selectors,
    overall = list(estimate = overall, se = as.double(oi_res$se), ci_lower = as.double(oi_res$ci_lower), ci_upper = as.double(oi_res$ci_upper),
      if_full = overall_if, status = overall_status, reason = overall_reason,
      inference_kind = if (isTRUE(oi_res$bootstrap)) "multiplier_bootstrap_iqr_se" else "analytic",
      n_draws = oi_res$n_draws, rng = oi_res$rng,
      draws_reused_from = if (reuse_simple_inference && !is.null(inf$multiplier)) "inference" else NULL,
      intended_bins = labels[outer], used_bins = labels[surv], excluded_bins = labels[setdiff(outer, surv)],
      weights = stats::setNames(ow, labels), share_gradient = stats::setNames(og, as.character(cohorts)), share_dependent = overall_share_dep,
      share_denominator = overall_share_denominator,
      numerical_certificate = oi_res$numerical_certificate)), class = "cdid_ipw_aggregation_result")
  if (first_masses) {
    # D61 record: cohort_shares then hold the masses H(I_g V) when V is
    # incomplete (degree-0 maps; never-treated V not demanded).
    out$share_law <- "first_calendar_masses"
    out$cohort_mass_complete <- mass_complete
    out$cohort_masses_unnormalized <- incomplete_masses
  }
  out
}

# ------------------------------------------------------------------------------
# 8. End-to-End Chained IPW-GMM Pipeline Entrypoint
# ------------------------------------------------------------------------------

# Typed support-only result: the full declared ledger with every cell
# unavailable, no accepted GMM fit, no quantile and no affirmative certificate.
ipw_pipeline_support_only_result <- function(status, reason, ledger, panel_owner, frame_owner, metadata_obj, n_total,
                                             reporting_map, cell_inference, agg_results, bp_mode, wt_mode, cg_mode,
                                             controls, accounting, nuisance_registry = NULL, selection = NULL) {
  structure(list(
    support_status = status, support_reason = reason,
    support = list(ledger = ledger, basis = NULL, n_requested = nrow(ledger), n_available = 0L,
      whole_request_complete = FALSE, unavailable = ledger[, c("cell_id", "g", "time", "support_status", "support_reason")],
      population_claim = "none: catalog-relative support only"),
    resource_accounting = accounting, panel_owner = panel_owner, frame_owner = frame_owner, metadata = metadata_obj,
    design_all = NULL, design_used = NULL, nuisance_registry = nuisance_registry, selection = selection, point_links_all = list(), point_links_used = list(),
    structural_certificate = NULL,
    reduction_summary = list(weighting = wt_mode, reduction_active = FALSE, n_links_all = 0L, n_links_used = 0L, kept_indices = integer(0L), dropped_indices = integer(0L), expanded_m = NULL, max_moment_expansion_diff = NA_real_),
    covariance_stream = NULL, gmm_fit = NULL, covariance_handling = NULL, numerical_transfer = NULL,
    auxiliary_reporting = if (is.null(reporting_map)) NULL else reporting_map$auxiliary,
    full_frame_transport = NULL,
    reporting = if (is.null(reporting_map)) list(base_period = bp_mode, target_meta = NULL, estimates = NULL, status = status) else list(
      base_period = bp_mode, target_meta = reporting_map$target_meta, L_map = NULL, estimates = cell_inference$estimates, se = cell_inference$se,
      analytical_se = cell_inference$analytical_se, numerical_certificate = cell_inference$numerical_certificate, z_stat = cell_inference$z_stat,
      p_value = cell_inference$p_value, ci_lower = cell_inference$ci_lower, ci_upper = cell_inference$ci_upper, cband_crit_val = cell_inference$cband_crit_val,
      cband_lower = cell_inference$cband_lower, cband_upper = cell_inference$cband_upper, vcov = cell_inference$vcov,
      Omega_a = reporting_map$Sigma_a, Omega_a_enclosure = reporting_map$Sigma_a_enclosure, if_full = reporting_map$if_full,
      available = reporting_map$available, support = cell_inference$support, inference = cell_inference,
      enclosure = reporting_map$enclosure),
    basis = NULL, aggregations = agg_results, controls = controls), class = "cdid_ipw_pipeline_result")
}

ipw_pipeline_use_raw_backend <- function(cp04, custom_est_method, run_numerical_transfer, cpt_mode, env) {
  if (!is.null(cp04) || is.function(custom_est_method) || isTRUE(run_numerical_transfer) || isTRUE(cpt_mode)) {
    return(FALSE)
  }
  if (is.environment(env) && !isNamespace(env) && !identical(env, .GlobalEnv)) {
    if (exists(".counts", envir = env, inherits = FALSE)) return(FALSE)
    legacy_spies <- c(
      "ipw_pipeline_batch_eligibility",
      "ipw_point_if_batch_plan",
      "ipw_point_if_general_plan",
      "ipw_point_if_retained_reserve",
      "ipw_point_if_general_state_reserve",
      "build_pipeline_nuisance_and_link_views",
      "build_pipeline_compact_basis_registry",
      "ipw_pipeline_compact_basis_facts",
      "fit_ipw_point_link_view_batch",
      "fit_pipeline_point_links"
    )
    p_env <- parent.env(env)
    for (sp in legacy_spies) {
      if (exists(sp, envir = env, inherits = FALSE)) {
        if (!exists("build_ipw_raw_backend_registry", envir = env, inherits = FALSE)) return(FALSE)
        if (isNamespace(p_env) && exists(sp, envir = p_env, inherits = FALSE) &&
            !identical(body(get(sp, envir = env, inherits = FALSE)),
                       body(get(sp, envir = p_env, inherits = FALSE)))) {
          return(FALSE)
        }
      }
    }
  }
  TRUE
}

run_chained_ipw_pipeline <- function(
  data,
  idname,
  tname,
  yname,
  gname,
  control,
  xformla = ~1,
  weightsname = NULL,
  clustername = NULL,
  cluster_vec = NULL,
  target_cohorts = NULL,
  anticipation = 0L,
  control_group = c("notyettreated", "nevertreated"),
  conventions = "D10_D13_admissible_frame",
  requested_post = NULL,
  declared_qualified = NULL,
  target_post_times = NULL,
  base_period = c("varying", "universal"),
  weighting = c("optimal", "identity", "two_step"),
  reduce_structural_redundancy = TRUE,
  run_numerical_transfer = TRUE,
  root_witness = NULL,
  aggregations = c("simple", "group", "dynamic", "calendar"),
  cband = FALSE,
  alp = 0.05,
  biters = 1000L,
  btype = "rademacher",
  seed = NULL,
  fixed_signs = NULL,
  frame_id = NULL,
  frame_Z = NULL,
  frame_weights = NULL,
  frame_g = NULL,
  allow_frame_only = FALSE,
  na.rm = FALSE,
  bstrap = NULL,
  caller_rng = FALSE,
  custom_est_method = NULL,
  custom_dots = list(),
  gmm_covariance = NULL,
  pt_assumption = "post", pt_start = NULL, cp04 = NULL,
  structural_reduction_stage = "outcome_free",
  compute_inffunc = TRUE,
  print_details = FALSE,
  conditional_pretest = FALSE,
  prepared = NULL
) {
  if (!gmm_scalar_logical(compute_inffunc) || !gmm_scalar_logical(print_details) ||
      !gmm_scalar_logical(conditional_pretest)) {
    cdid_ipw_pipeline_abort("`compute_inffunc`, `print_details` and `conditional_pretest` must be scalar logicals.",
      "cdid_ipw_pipeline_input_error")
  }
  # CP08 O14: the conditional pretrend test needs the fitted union registry
  # (post plus both auxiliary catalogs) and the auxiliary links only. It stops
  # before structural reduction, the GMM solve and every post-ATT inference or
  # identification step, and charges no GMM/inference/draw reservation.
  if (isTRUE(conditional_pretest) && (is.function(custom_est_method) || isTRUE(run_numerical_transfer) ||
      length(aggregations) || !identical(bstrap, FALSE) || !identical(cband, FALSE)))
    cdid_ipw_pipeline_abort("The conditional pretest state is built for built-in routes without transfer, draws or aggregations.",
      "cdid_ipw_pipeline_input_error")
  progress <- chained_progress_emitter(print_details)
  if (is.null(progress)) progress <- function(...) invisible(NULL) else progress("chained pipeline started")
  pt <- validate_chained_pt_request(pt_assumption, pt_start)
  pt_assumption <- pt$pt_assumption; pt_start <- pt$pt_start
  if (!is.null(cp04$or_method)) conventions <- paste0(conventions, "|OR-control-mean-", cp04$or_method)
  if (!is.null(cp04$dr_method)) conventions <- paste0(conventions, "|DR-", cp04$dr_method)
  # D60: dated baseline IPW and custom callbacks run on the full-roster reference
  # backend with recipe-indexed designs and anchored weights.
  if (!is.null(cp04$ipw_method)) conventions <- paste0(conventions, "|IPW-", cp04$ipw_method,
    if (is.function(custom_est_method)) "-custom-callback" else "", "-reference-recipes")
  # Pure imputation (D52) fits no observation-probability model; the augmented
  # correction (D52 seventh route) fits both the imputation and the
  # observation-probability arms on F.
  imputation <- identical(cp04$selection_est_method, "imputation")
  augmented <- identical(cp04$selection_est_method, "dr")
  if (imputation) {
    conventions <- paste0(conventions, "|IMP-", cp04$imputation_est_method, "-", cp04$imputation_pooling,
      "-F-actual-change-fixed-w")
  } else if (augmented) {
    conventions <- paste0(conventions, "|AUG-DR-", cp04$imputation_pooling, "-", cp04$selection_model, "-",
      cp04$selection_pooling, "-F-actual-change-fixed-w-Z=j+S(D-j)/q")
  } else if (!is.null(cp04) && !identical(cp04$or_method, "baseline") && !identical(cp04$dr_method, "baseline") &&
             is.null(cp04$ipw_method)) {
    sel_m <- if (is.null(cp04$selection_model)) "pair_logit" else cp04$selection_model
    sel_p <- if (is.null(cp04$selection_pooling)) "separate" else cp04$selection_pooling
    conventions <- if (identical(sel_m, "pair_logit") && identical(sel_p, "separate")) {
      paste0(conventions, "|CP04-", cp04$universe, "-weighted-cohort-pair-logit-aU/(rho*pU)")
    } else {
      paste0(conventions, "|CP04-", cp04$universe, "-", sel_m, "-", sel_p, "-aU/(rho*pU)")
    }
  }
  or_preprocess_entries <- if (is.null(cp04$or_preprocess)) 0 else cp04$or_preprocess$retained_entries
  or_preprocess_work <- if (is.null(cp04$or_preprocess)) 0 else cp04$or_preprocess$work_ops
  bp_mode <- match.arg(base_period)
  w_choice <- resolve_gmm_weighting_choice(
    weighting,
    arg = "weighting",
    error_fn = cdid_ipw_pipeline_abort,
    error_class = "cdid_ipw_pipeline_input_error"
  )
  wt_mode <- w_choice$canonical
  wt_requested <- w_choice$requested
  cg_mode <- match.arg(control_group)
  na_rm <- na.rm
  if (!is.null(custom_est_method) && !is.function(custom_est_method)) {
    cdid_ipw_pipeline_abort("`custom_est_method` must be NULL or a function.", "cdid_ipw_pipeline_input_error")
  }
  if (!is.list(custom_dots) || is.object(custom_dots)) {
    cdid_ipw_pipeline_abort("`custom_dots` must be a plain list.", "cdid_ipw_pipeline_input_error")
  }
  if (is.function(custom_est_method) && isTRUE(run_numerical_transfer)) {
    cdid_ipw_pipeline_abort(
      "Custom function-valued `est_method` does not support built-in interval numerical transfer (`run_numerical_transfer = TRUE`).",
      "cdid_ipw_pipeline_control_error"
    )
  }
  # O15 (D64): custom callbacks run on the baseline reference route and, with
  # their explicit log-weight gradient, on conditional_mtar IPW (selection
  # models); OR/DR/imputation/augmented callback interfaces stay unsupported.
  cp04_custom_ok <- identical(cp04$ipw_method, "baseline") ||
    (is.null(cp04$or_method) && is.null(cp04$dr_method) && is.null(cp04$selection_est_method) &&
      identical(cp04$canonical_assumption, "conditional_mtar"))
  if (!is.null(cp04) && ((is.function(custom_est_method) && !cp04_custom_ok) ||
      isTRUE(run_numerical_transfer)))
    cp04_abort("The full-roster backend supports built-in DR, IPW or OR without numerical transfer; corrected callbacks/known probabilities are unsupported.", "unsupported")
  # Flag composition is resolved once, before any owner is constructed:
  # `bstrap = NULL` keeps the private law (bootstrap iff `cband`); an explicit
  # `bstrap` separates SEs from bands; FALSE/TRUE follows the released endpoint
  # (reported cells/basis: pointwise, no draw; fit-time group/dynamic/calendar
  # aggregations: multiplier band with analytic SEs). Draw storage is reserved
  # whenever any endpoint of this call will draw.
  boot_flags <- ipw_pipeline_resolve_bootstrap_flags(bstrap, cband, where = "pipeline inference")
  if (!gmm_scalar_logical(caller_rng)) cdid_ipw_pipeline_abort("`caller_rng` must be a scalar logical.", "cdid_ipw_pipeline_input_error")
  if (isTRUE(caller_rng) && !is.null(seed)) {
    cdid_ipw_pipeline_abort("`caller_rng = TRUE` uses the caller's RNG stream; an explicit `seed` cannot also be supplied.", "cdid_ipw_pipeline_bootstrap_error")
  }
  ctrl <- validate_ipw_pipeline_control(control, weighting = wt_mode)
  if (!is.null(gmm_covariance)) {
    cov_method <- resolve_gmm_covariance_choice(
      gmm_covariance,
      arg = "gmm_covariance",
      error_fn = cdid_ipw_pipeline_abort,
      error_class = "cdid_ipw_pipeline_input_error"
    )
    if (!is.null(ctrl$gmm_accuracy$covariance_strategy)) {
      ctrl$gmm_accuracy$covariance_strategy$method <- cov_method
      ctrl$gmm_accuracy$covariance_strategy$weighting_requested <- wt_requested
      ctrl$gmm_accuracy$covariance_strategy$weighting_canonical <- wt_mode
    } else if (identical(ctrl$gmm_accuracy$method, "conventional_chol_qr")) {
      ctrl$gmm_accuracy$covariance_strategy <- gmm_covariance_strategy(
        method = cov_method,
        weighting_requested = wt_requested,
        weighting_canonical = wt_mode
      )
    }
  } else if (!is.null(ctrl$gmm_accuracy$covariance_strategy)) {
    ctrl$gmm_accuracy$covariance_strategy$weighting_requested <- wt_requested
    ctrl$gmm_accuracy$covariance_strategy$weighting_canonical <- wt_mode
  }
  if (isTRUE(run_numerical_transfer) && (is.null(ctrl$inference_accuracy) ||
      identical(ctrl$gmm_accuracy$method, "conventional_chol_qr")))
    cdid_ipw_pipeline_abort("Private numerical transfer requires explicit inference budgets and a strict GMM certificate method.", "cdid_ipw_pipeline_control_error")
  if (isTRUE(run_numerical_transfer) && (missing(root_witness) || !is.function(root_witness))) {
    cdid_ipw_pipeline_abort("Numerical transfer requires a caller-supplied root_witness(fitted_view) returning F and r.", "cdid_ipw_pipeline_control_error")
  }

  if (!gmm_scalar_logical(reduce_structural_redundancy) ||
      !gmm_scalar_logical(run_numerical_transfer) ||
      !gmm_scalar_logical(cband) ||
      !gmm_scalar_logical(allow_frame_only) ||
      !gmm_scalar_logical(na_rm)) {
    cdid_ipw_pipeline_abort(
      "Flags `reduce_structural_redundancy`, `run_numerical_transfer`, `cband`, `allow_frame_only` and `na.rm` must be scalar logicals.",
      "cdid_ipw_pipeline_input_error"
    )
  }
  # Private verification switch, not a public option: "post_point_fit" runs the
  # accepted path that fits every row before certifying, as an independent
  # reference for the licensed outcome-free route. Other routes ignore it.
  if (!gmm_scalar_character(structural_reduction_stage) ||
      !structural_reduction_stage %in% c("outcome_free", "post_point_fit")) {
    cdid_ipw_pipeline_abort("`structural_reduction_stage` must be \"outcome_free\" or \"post_point_fit\".",
      "cdid_ipw_pipeline_input_error")
  }

  if (xor(is.null(frame_id),is.null(frame_Z))) {
    cdid_ipw_pipeline_abort("Supply frame_id and frame_Z together; partial frame input is refused.", "cdid_ipw_pipeline_input_error")
  }

  # CP08 common boundary: a sealed prepared intake record of the request layer
  # is bound here by object identity to exactly the data and frame arrays this
  # call received (and, with a full-roster backend, its baseline reader arrays).
  # When bound, the panel/frame/reader constructors below reuse its validated
  # positions, orders and counts instead of rescanning and re-matching; every
  # other call (NULL, a different object, a mutated copy, data weights, cluster
  # column, frame-only units, custom callbacks, reference recipes) keeps the
  # standalone constructors unchanged. Nothing is retried after a failure.
  bound <- chained_request_prepared_bind(prepared, data, idname, tname, yname, gname, frame_id, frame_Z, frame_weights,
    frame_g, allow_frame_only, weightsname, clustername, cp04, custom_est_method)

  # 1. Construct sealed sparse panel owner
  panel_owner <- if (is.null(bound)) construct_ipw_panel_owner(
    data = data,
    idname = idname,
    tname = tname,
    yname = yname,
    gname = gname,
    control = ctrl$owner_control
  ) else construct_ipw_panel_owner_prepared(
    data = data,
    idname = idname,
    tname = tname,
    yname = yname,
    gname = gname,
    control = ctrl$owner_control,
    prepared = bound
  )

  # 2. Extract or validate unit-level frame (id, Z, weights, g, and optional cluster)
  if (is.null(frame_id) || is.null(frame_Z)) {
    id_col <- data[[idname]]
    first_pos <- which(!duplicated(id_col))
    unit_df <- data[first_pos, , drop = FALSE]
    rownames(unit_df) <- NULL

    f_id <- unit_df[[idname]]
    f_g <- unit_df[[gname]]
    f_Z <- stats::model.matrix(xformla, data = unit_df)
    if (nrow(f_Z) != nrow(unit_df)) {
      cdid_ipw_pipeline_abort(
        "Covariate formula `xformla` dropped rows (e.g. due to missing covariates); NA covariates are refused.",
        "cdid_ipw_pipeline_input_error"
      )
    }
    full_Z <- stats::model.matrix(xformla, data = data)
    if (nrow(full_Z) != nrow(data)) {
      cdid_ipw_pipeline_abort(
        "Covariate formula `xformla` contains missing values in panel rows.",
        "cdid_ipw_pipeline_input_error"
      )
    }
    u_match <- match(id_col, f_id)
    if (max(abs(full_Z - f_Z[u_match, , drop = FALSE])) > 0.0) {
      cdid_ipw_pipeline_abort(
        "Covariates in `xformla` vary within unit across time; unit-level time-invariant covariates are required.",
        "cdid_ipw_pipeline_covariate_error"
      )
    }

    if (!is.null(weightsname)) {
      if (!is.character(weightsname) || length(weightsname) != 1L || !(weightsname %in% names(data))) {
        cdid_ipw_pipeline_abort("`weightsname` must name a column in `data`.", "cdid_ipw_pipeline_input_error")
      }
      w_col <- as.double(data[[weightsname]])
      f_w <- as.double(unit_df[[weightsname]])
      if (anyNA(w_col) || max(abs(w_col - f_w[u_match])) > 0.0) {
        cdid_ipw_pipeline_abort(
          "Sampling weights in `weightsname` must be finite and time-invariant within unit.",
          "cdid_ipw_pipeline_weight_error"
        )
      }
    } else if (!is.null(frame_weights)) {
      f_w <- as.double(frame_weights)
    } else {
      f_w <- rep(1.0, length(f_id))
    }
  } else {
    f_id <- frame_id
    f_Z <- as.matrix(frame_Z)
    f_w <- if (!is.null(frame_weights)) as.double(frame_weights) else rep(1.0, length(f_id))
    f_g <- frame_g
  }

  frame_owner <- if (is.null(bound)) construct_ipw_frame_owner(
    panel_owner = panel_owner,
    id = f_id,
    Z = f_Z,
    weights = f_w,
    control = ctrl$owner_control,
    frame_g = f_g,
    allow_frame_only = allow_frame_only, design_on_A_only = !is.null(cp04)
  ) else construct_ipw_frame_owner_prepared(
    panel_owner = panel_owner,
    id = frame_id,
    Z = frame_Z,
    weights = frame_weights,
    control = ctrl$owner_control,
    frame_g = frame_g,
    prepared = bound, design_on_A_only = !is.null(cp04)
  )

  cp04_baseline <- !is.null(cp04) && (identical(cp04$or_method, "baseline") || identical(cp04$dr_method, "baseline") ||
    identical(cp04$ipw_method, "baseline"))
  # D60: on the reference route the reader keeps identities, cohorts, U and rho;
  # its mechanical unit weights and carrier design are blanked (NA) after the
  # identity checks, so every weight, design and rate is read at its recipe.
  cp04_reader <- if (is.null(cp04)) NULL else if (is.null(bound)) construct_cp04_frame_reader(frame_owner,
    f_id, cp04$Z, if (cp04_baseline) cp04$frame_g else f_g, f_w, ctrl$owner_control, U = cp04$U, rho = cp04$rho,
    universe = cp04$universe, justification = cp04$justification,
    allow_unknown_off_A = cp04_baseline, reference_carrier = !is.null(cp04$reference_source)) else
    construct_cp04_frame_reader_prepared(frame_owner,
      frame_id, cp04$Z, cp04$frame_g, frame_weights, ctrl$owner_control, prepared = bound, U = cp04$U, rho = cp04$rho,
      universe = cp04$universe, justification = cp04$justification,
      allow_unknown_off_A = TRUE, reference_carrier = FALSE)
  # Invocation-only maps have no later consumer. Retained replay closures must
  # not keep this sealed intake (and its source/ordering arrays) alive.
  prepared <- bound <- NULL
  replay_n <- if (is.null(cp04)) panel_owner$n_a_units else length(f_id)
  n_total <- as.integer(length(f_id))
  n_a_units <- as.integer(panel_owner$n_a_units)
  a_to_full <- match(panel_owner$a_id, f_id)
  if (anyNA(a_to_full) || length(a_to_full) != n_a_units) {
    cdid_ipw_pipeline_abort(
      "Failed to align `panel_owner$a_id` into the full unit frame `f_id`.",
      "cdid_ipw_pipeline_integrity_error"
    )
  }
  full_frame_g <- if (!is.null(f_g)) as.double(f_g) else {
    out_g <- rep(NA_real_, n_total)
    # The sealed owner retains cohort metadata for singleton and unobserved
    # panel units as well. Only genuinely unknown frame-only IDs remain NA.
    known <- match(f_id,panel_owner$units$id)
    out_g[!is.na(known)] <- as.double(panel_owner$units$g[known[!is.na(known)]])
    out_g
  }
  full_frame_w <- as.double(f_w)
  # D61: a reference-input fit aggregates with first-calendar masses V = w_t1 by
  # original id (NA where absent); the owner's carrier unit weights never enter.
  share_law <- NULL
  if (!is.null(cp04$reference_source)) {
    full_frame_w <- as.double(cp04$reference_source$V)
    share_law <- "first_calendar_masses"
  }

  # Optional one-way cluster assignment on the full frame
  cluster_full <- NULL
  if (!is.null(cluster_vec)) {
    if (length(cluster_vec) != n_total || anyNA(cluster_vec)) {
      cdid_ipw_pipeline_abort(
        "`cluster_vec` must be a non-NA vector of length `n_total` (matching the full unit frame).",
        "cdid_ipw_pipeline_cluster_error"
      )
    }
    cluster_full <- cluster_vec
  } else if (!is.null(clustername)) {
    if (!is.character(clustername) || length(clustername) != 1L || !(clustername %in% names(data))) {
      cdid_ipw_pipeline_abort("`clustername` must name a column in `data`.", "cdid_ipw_pipeline_cluster_error")
    }
    cl_panel <- data[[clustername]]
    if (anyNA(cl_panel)) {
      cdid_ipw_pipeline_abort("Cluster column contains NA values.", "cdid_ipw_pipeline_cluster_error")
    }
    id_col <- data[[idname]]
    u_pos <- match(f_id, id_col)
    if (anyNA(u_pos)) {
      cdid_ipw_pipeline_abort(
        "`clustername` cannot populate frame-only units absent from `data`; supply `cluster_vec` on the full frame.",
        "cdid_ipw_pipeline_cluster_error"
      )
    }
    cluster_full <- cl_panel[u_pos]
    row_to_frame <- match(id_col, f_id)
    if (!identical(as.character(cl_panel), as.character(cluster_full[row_to_frame]))) {
      cdid_ipw_pipeline_abort(
        "Cluster assignment varies within unit across panel rows; units must be nested within clusters.",
        "cdid_ipw_pipeline_cluster_error"
      )
    }
  }
  eligible_full <- seq_len(n_total) %in% a_to_full

  # 3. Declared target cohorts and the D5 requested-cell ledger (before any
  #    support-dependent step; independent of retained rows and fits).
  max_coh <- max(panel_owner$cohorts)
  if (is.null(target_cohorts)) {
    # All finite cohorts strictly below the latest clean-control cohort
    cand_g <- panel_owner$cohorts[is.finite(panel_owner$cohorts) & panel_owner$cohorts < max_coh]
    if (length(cand_g) < 1L) {
      cdid_ipw_pipeline_abort(
        "No treated target cohorts strictly precede the latest control cohort.",
        c("cdid_ipw_pipeline_target_error", "cdid_catalog_target_error")
      )
    }
    tgt_cohorts <- sort(unique(as.double(cand_g)))
  } else {
    tgt_cohorts <- sort(unique(as.double(target_cohorts)))
  }
  validate_chained_pt_request(pt_assumption, pt_start, panel_owner$all_times)
  ledger <- ipw_pipeline_request_ledger(tgt_cohorts, panel_owner$all_times, panel_owner$cohorts, anticipation,
    cg_mode, target_post_times = target_post_times, requested_post = requested_post)
  R_full <- nrow(ledger)
  controls_out <- list(
    anticipation = as.integer(anticipation), control_group = cg_mode, conventions = conventions, base_period = bp_mode,
    pt_assumption = pt_assumption, pt_start = pt_start,
    weighting = wt_mode, weighting_requested = wt_requested,
    gmm_covariance = if (!is.null(ctrl$gmm_accuracy$covariance_strategy)) ctrl$gmm_accuracy$covariance_strategy$method else NULL,
    covariance_strategy = ctrl$gmm_accuracy$covariance_strategy,
    reduce_structural_redundancy = reduce_structural_redundancy, run_numerical_transfer = run_numerical_transfer,
    root_witness_supplied = is.function(root_witness), custom_est_method_supplied = is.function(custom_est_method),
    cband = cband, alp = alp, biters = as.integer(biters), btype = btype,
    seed = if (is.null(seed)) NULL else as.integer(seed), fixed_signs_supplied = !is.null(fixed_signs), na_rm = na_rm,
    bstrap = bstrap, caller_rng = caller_rng, bootstrap_requested = boot_flags$bootstrap, flag_law = boot_flags$law,
    or_basis = cp04$or_basis, outcome_recipe = cp04$outcome_recipe, treatment_columns = cp04$treatment_columns,
    pipeline_control = ctrl)
  agg_types <- if (!is.null(aggregations) && length(aggregations) > 0L) unique(as.character(aggregations)) else character(0L)
  reserve_draws <- boot_flags$draws ||
    (isTRUE(boot_flags$released_false_true) && any(agg_types %in% c("group", "dynamic", "calendar")))
  controls_out$draws_reserved <- reserve_draws
  controls_out$compute_inffunc <- compute_inffunc
  controls_out$print_details <- print_details
  # compute_inffunc = FALSE (did 2.5.1): point estimates only. The GMM weights
  # still use the complete internal link IFs/covariance; nothing after the
  # solve replays an IF, so draws, aggregations and transfer cannot be served.
  if (!isTRUE(compute_inffunc) && (reserve_draws || length(agg_types) || isTRUE(run_numerical_transfer))) {
    cdid_ipw_pipeline_abort(paste0("compute_inffunc = FALSE computes point estimates only; multiplier draws, fit-time ",
      "aggregations and numerical transfer require influence functions."), "cdid_ipw_pipeline_input_error")
  }
  run_aggregations <- function(rm) {
    out <- list()
    for (agg_t in agg_types) {
      out[[agg_t]] <- aggregate_chained_ipw_targets(
        theta_hat = rm$theta_hat, target_if_full = rm$if_full, target_meta = rm$target_meta,
        full_frame_g = full_frame_g, full_frame_w = full_frame_w, type = agg_t, anticipation = anticipation,
        cluster_full = cluster_full, alp = alp, cband = cband, biters = biters, btype = btype, seed = seed,
        multiplier_budget = ctrl$multiplier_budget, inference_accuracy = ctrl$inference_accuracy,
        enclosure = rm$enclosure, fixed_signs = fixed_signs, eligible_full = eligible_full, na.rm = na_rm, available = rm$available,
        bstrap = bstrap, caller_rng = caller_rng, calendar = panel_owner$all_times, cohort_grid = panel_owner$cohorts,
        share_law = share_law, frame_ids = f_id)
    }
    out
  }
  run_cell_inference <- function(rm) {
    if (!isTRUE(compute_inffunc))
      return(ipw_pipeline_point_only_inference(rm$theta_hat, rm$available, n_total, cluster_full, eligible_full, alp))
    ipw_pipeline_supported_inference(rm$theta_hat, rm$if_full, rm$available, cluster_full, alp, cband, biters, btype, seed,
      ctrl$multiplier_budget, ctrl$inference_accuracy, rm$enclosure, fixed_signs, eligible_full,
      reason = "no reported cell is supported by the retained comparison catalog",
      bstrap = bstrap, caller_rng = caller_rng, target_meta = rm$target_meta)
  }

  conditional_state <- function(status, reason, registry = NULL, reporting_map = NULL, aux_plan = NULL,
                                metadata_obj = NULL, combined_initial = 0, combined_work = 0) {
    structure(list(status = status, reason = reason, panel_owner = panel_owner, frame_owner = frame_owner,
      cp04_reader = cp04_reader, cp04 = cp04, registry = registry, auxiliary = reporting_map$auxiliary,
      aux_plan = aux_plan, ledger = ledger, metadata = metadata_obj, base_period = bp_mode, control_group = cg_mode,
      anticipation = as.double(anticipation), n_total = n_total, n_a_units = n_a_units, a_to_full = a_to_full,
      frame_ids = f_id, cluster_full = cluster_full, eligible_full = eligible_full, control = ctrl,
      conventions = conventions,
      resource = list(retained_entries = combined_initial, work_ops = combined_work,
        registry_entries = if (is.null(registry$cp04)) 0 else registry$cp04$accounting$owned_entries,
        reader_entries = if (is.null(cp04_reader)) 0 else cp04_reader$accounting$owned_entries)),
      class = "cdid_conditional_pretest_state")
  }
  if (n_a_units < 2L && isTRUE(conditional_pretest)) {
    return(conditional_state("insufficient_pair_observing_units",
      paste0("Admissible pair-observing unit count `n_a_units` = ", n_a_units, " is less than 2; no auxiliary comparison ",
        "has treated and control support.")))
  }
  if (n_a_units < 2L) {
    ledger <- ipw_pipeline_ledger_support(ledger, NULL, NULL, global_status = list(status = "insufficient_pair_observing_units",
      reason = paste0("admissible pair-observing unit count n_a_units = ", n_a_units, " is below 2; no comparison pair can have treated and control support")))
    return(ipw_pipeline_support_only_result("insufficient_pair_observing_units",
      paste0("Admissible pair-observing unit count `n_a_units` = ", n_a_units, " is less than 2; no chained IPW-GMM comparison is supported."),
      ledger, panel_owner, frame_owner, NULL, n_total, NULL, NULL, list(), bp_mode, wt_mode, cg_mode, controls_out,
      list(reachable_numeric_payload_only = TRUE, peak_heap_or_RSS_claim = FALSE, d5_arrays = list(n_ledger_cells = R_full, n_basis_contrasts = 0L))))
  }

  # CP08 scaling action 3: the union of both auxiliary reporting plans is built
  # once, from the ledger alone, before any covariance/replay/report
  # reservation, and reused by every consumer. Its own payload and work are
  # charged with the live owners first; it persists through reporting.
  owners_live <- ipw_owner_safe_sum(c(or_preprocess_entries, panel_owner$accounting$actual_retained_entries,
    length(f_Z), length(f_w), length(f_id), length(f_g)))
  aux_env <- ipw_auxiliary_spec_envelope(length(unique(ledger$g)), length(panel_owner$all_times), length(panel_owner$cohorts))
  ipw_owner_safe_sum(c(owners_live, aux_env$retained_entries),
    min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_output_entries),
    "pipeline owners plus the auxiliary reporting plan")
  ipw_owner_safe_sum(c(owners_live, aux_env$retained_entries, aux_env$scratch_entries), ctrl$owner_control$max_matrix_entries,
    "pipeline owners plus auxiliary reporting plan construction")
  ipw_owner_safe_sum(c(or_preprocess_work, aux_env$work_ops), ctrl$owner_control$max_work_ops,
    "auxiliary reporting plan construction work")
  aux_plan <- ipw_auxiliary_spec(ledger, panel_owner$all_times, panel_owner$cohorts, anticipation, cg_mode)
  Q_aux <- as.double(length(aux_plan$descriptors))
  Q_rep <- as.double(length(aux_plan$selections[[bp_mode]]))
  aux_record <- list(n_descriptors_union = Q_aux,
    n_selections = c(varying = length(aux_plan$selections$varying), universal = length(aux_plan$selections$universal)),
    n_reported_pretrend = Q_rep, base_period = bp_mode,
    n_reference_zero = sum(vapply(aux_plan$selections[[bp_mode]], function(d) identical(d$status, "reference_zero"), logical(1L))),
    n_unavailable = sum(vapply(aux_plan$selections[[bp_mode]], function(d) identical(d$status, "unavailable"), logical(1L))),
    retained_entries_bound = aux_env$retained_entries, work_ops_bound = aux_env$work_ops,
    replaced_bound = 2 * length(tgt_cohorts) * length(panel_owner$all_times),
    meaning = paste("exact union of both base-period plans (fitted auxiliary comparisons) and the reported pretrend",
      "columns of the selected family; it replaces the former 2 J T reservation bound"))
  progress("auxiliary reporting plan: ", Q_aux, " distinct comparisons (", aux_record$n_selections[["varying"]],
    " varying, ", aux_record$n_selections[["universal"]], " universal reports; ", Q_rep, " reported for base_period = \"",
    bp_mode, "\")")

  # First build catalog (or if `declared_qualified` & `target_post_times` are NULL,
  # qualify all clean causal endpoints `(g, u)` present in the retained catalog rows).
  progress("building the comparison catalog for target cohorts ", paste(fmt_design_num(tgt_cohorts), collapse = ", "))
  meta_initial <- build_owner_comparison_catalog(
    panel_owner = panel_owner,
    targets = tgt_cohorts,
    anticipation = anticipation,
    control_group = cg_mode,
    conventions = conventions,
    requested_post = requested_post,
    control = ctrl$owner_control,
    catalog_bounds = ctrl$catalog_bounds,
    frame_owner = frame_owner,
    build_design = FALSE, pt_assumption = pt_assumption, pt_start = pt_start
  )
  if (nrow(meta_initial$catalog$rows) < 1L && isTRUE(conditional_pretest)) {
    # Auxiliary links only: owners, the plan and the auxiliary views/nuisances
    # with their capture metadata (no GMM, report, interval or draw buffer).
    ledger <- ipw_pipeline_ledger_support(ledger, meta_initial, NULL, global_status = list(status = "empty_catalog",
      reason = "the comparison catalog retained zero links for the requested target cohorts and control group"))
    cpt_initial <- ipw_owner_safe_sum(c(owners_live, aux_env$retained_entries,
      32.0 * Q_aux * (ncol(f_Z)^2 + ncol(f_Z) * length(panel_owner$cohorts) + 64.0),
      2048.0 * Q_aux + 2048.0 + 8.0 * min(n_a_units, ctrl$owner_control$max_chunk_rows)),
      min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_matrix_entries),
      "pipeline owners, auxiliary plan and auxiliary views/nuisances (conditional pretest, empty catalog)")
    cpt_work <- ipw_owner_safe_sum(c(or_preprocess_work, aux_env$work_ops, Q_aux * (128.0 * max(1, n_a_units) + 4096.0)),
      ctrl$owner_control$max_work_ops, "auxiliary fitting work (conditional pretest, empty catalog)")
    empty_registry <- if (is.null(cp04)) NULL else build_cp04_public_registry(NULL, frame_owner,
      cp04_reader, ledger, ctrl, cg_mode, anticipation = anticipation,
      retained_pipeline_entries = cpt_initial, external_work_ops = cpt_work, or_method = cp04$or_method,
      automatic_chunk_rows = isTRUE(cp04$automatic_chunk_rows), outcome_recipe = cp04$outcome_recipe,
      selection_model = if (is.null(cp04$selection_model)) "pair_logit" else cp04$selection_model,
      selection_pooling = if (is.null(cp04$selection_pooling)) "separate" else cp04$selection_pooling,
      dr_method = cp04$dr_method, treatment_columns = cp04$treatment_columns,
      imputation = if (imputation) cp04 else NULL, augmented = if (augmented) cp04 else NULL,
      reference_source = cp04$reference_source, ipw_method = cp04$ipw_method, aux_spec = aux_plan)
    reporting_map <- build_ipw_auxiliary_reporting(ledger, panel_owner$all_times, panel_owner$cohorts, anticipation, cg_mode,
      empty_registry, frame_owner, NULL, ctrl, NULL, NULL, FALSE, bp_mode, n_total, a_to_full, cpt_initial,
      pt_assumption = pt_assumption, pt_start = pt_start, aux_spec = aux_plan, compute_if = FALSE)
    return(conditional_state("available", NA_character_, empty_registry, reporting_map, aux_plan, meta_initial,
      cpt_initial, cpt_work))
  }
  if (nrow(meta_initial$catalog$rows) < 1L) {
    # Entirely empty catalog: support metadata, direct diagnostics where
    # supported, no accepted GMM fit/certificate.
    ledger <- ipw_pipeline_ledger_support(ledger, meta_initial, NULL, global_status = list(status = "empty_catalog",
      reason = "the comparison catalog retained zero links for the requested target cohorts and control group"))
    # Exact plan dimensions: Q_rep reported pretrend columns plus every ledger
    # cell; Q_aux fitted auxiliary comparisons (union of both families).
    J_empty <- R_full + Q_rep
    # Retained aggregation diagnostic boxes (W, D and point boxes per aggregation type).
    diag_entries_empty <- (if (run_numerical_transfer) 2.0 else 1.0) * length(agg_types) * J_empty * (J_empty + length(panel_owner$cohorts) + 1.0) +
      8.0 * length(agg_types) + 1.0
    combined_empty <- ipw_owner_safe_sum(c(owners_live, aux_env$retained_entries,
      (if (run_numerical_transfer) 64.0 else 16.0) * J_empty * J_empty,
      (if (run_numerical_transfer) 128.0 else 32.0) * n_total * J_empty,
      if (reserve_draws) 16.0 * as.double(biters) * J_empty * (length(agg_types) + 1.0) else 0.0,
      diag_entries_empty,
      # D27 compact capture outputs plus one bounded summary scratch buffer.
      if (is.function(custom_est_method)) 0.0 else
        2048.0 * Q_aux + 2048.0 + 8.0 * min(n_a_units, ctrl$owner_control$max_chunk_rows)),
      min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_matrix_entries),
      "pipeline owners, ledger and direct reporting arrays (empty catalog)")
    empty_work <- ipw_owner_safe_sum(c(or_preprocess_work, aux_env$work_ops,
      (if (run_numerical_transfer) 128.0 else 32.0) * n_total * (J_empty * J_empty + Q_aux * Q_aux + ncol(f_Z)^2),
      if (reserve_draws) (if (run_numerical_transfer) 64.0 else 8.0) * as.double(biters) * n_total * J_empty else 0.0,
      if (is.function(custom_est_method)) 0.0 else Q_aux * (128.0 * max(1, n_a_units) + 4096.0)),
      ctrl$owner_control$max_work_ops, "empty-catalog reporting/aggregation/interval/multiplier whole-path work envelope")
    empty_registry <- if (is.null(cp04)) NULL else build_cp04_public_registry(NULL, frame_owner,
      cp04_reader, ledger, ctrl, cg_mode, anticipation = anticipation,
      retained_pipeline_entries = combined_empty, external_work_ops = empty_work, or_method = cp04$or_method,
      automatic_chunk_rows = isTRUE(cp04$automatic_chunk_rows), outcome_recipe = cp04$outcome_recipe,
      selection_model = if (is.null(cp04$selection_model)) "pair_logit" else cp04$selection_model,
      selection_pooling = if (is.null(cp04$selection_pooling)) "separate" else cp04$selection_pooling,
      dr_method = cp04$dr_method, treatment_columns = cp04$treatment_columns,
      imputation = if (imputation) cp04 else NULL, augmented = if (augmented) cp04 else NULL,
      reference_source = cp04$reference_source, ipw_method = cp04$ipw_method,
      custom_est_method = custom_est_method, custom_dots = custom_dots, aux_spec = aux_plan)
    reporting_map <- build_ipw_auxiliary_reporting(ledger, panel_owner$all_times, panel_owner$cohorts, anticipation, cg_mode,
      empty_registry, frame_owner, NULL, ctrl, root_witness, NULL, isTRUE(run_numerical_transfer), bp_mode, n_total, a_to_full, combined_empty,
      custom_est_method = custom_est_method, custom_dots = custom_dots,
      pt_assumption = pt_assumption, pt_start = pt_start,
      pt_catalog_links = NULL, aux_spec = aux_plan, compute_if = compute_inffunc)
    cell_inference <- run_cell_inference(reporting_map)
    agg_results <- run_aggregations(reporting_map)
    return(ipw_pipeline_support_only_result("empty_catalog",
      "Comparison catalog retained zero links for the requested target cohorts and control group; every requested post cell is unavailable (catalog-relative).",
      ledger, panel_owner, frame_owner, meta_initial, n_total, reporting_map, cell_inference, agg_results, bp_mode, wt_mode, cg_mode, controls_out,
      list(initial_envelope_entries = combined_empty, downstream_work_ops = empty_work, auxiliary_plan = aux_record,
        reachable_numeric_payload_only = TRUE, peak_heap_or_RSS_claim = FALSE,
        d5_arrays = list(n_ledger_cells = R_full, n_basis_contrasts = 0L, n_reported_columns = nrow(reporting_map$target_meta),
          ledger_if_entries = as.double(n_total) * R_full, reported_if_entries = as.double(n_total) * nrow(reporting_map$target_meta),
          retained_diagnostic_entries = diag_entries_empty)), nuisance_registry = empty_registry,
      selection = cp04_selection_record(cp04, empty_registry$cp04, cp04_reader)))
  }

  # Determine `declared_qualified` if not explicitly provided
  if (is.null(declared_qualified)) {
    c_rows <- meta_initial$catalog$rows
    s_causal <- c_rows$s_role %in% c("anticipation_conditional", "post_conditional")
    t_causal <- c_rows$t_role %in% c("anticipation_conditional", "post_conditional")
    dq_g <- c(c_rows$g[s_causal], c_rows$g[t_causal])
    dq_t <- c(c_rows$s[s_causal], c_rows$t[t_causal])
    if (length(dq_g) > 0L) {
      dq_df <- unique(data.frame(g = as.double(dq_g), t = as.double(dq_t), stringsAsFactors = FALSE))
      dq_df <- dq_df[order(dq_df$g, dq_df$t), , drop = FALSE]
      rownames(dq_df) <- NULL
    } else {
      dq_df <- NULL
    }
  } else {
    dq_df <- declared_qualified
  }

  # Only row-bearing cohorts' requests enter the existing design bridge; rowless
  # cohorts keep their ledger cells as catalog-relative unsupported cells.
  row_cohorts <- unique(as.double(meta_initial$catalog$rows$g))
  in_design <- ledger$g %in% row_cohorts
  tpt_df <- data.frame(g = ledger$g[in_design], t = ledger$time[in_design], stringsAsFactors = FALSE)
  # A zero-row target table is valid: other declared cohorts may contribute
  # moments while every requested cell belongs to a rowless cohort. The design
  # still represents all catalog endpoints and the GMM checks still apply.

  metadata_obj <- build_owner_comparison_catalog(
    panel_owner = panel_owner,
    targets = tgt_cohorts,
    anticipation = anticipation,
    control_group = cg_mode,
    conventions = conventions,
    requested_post = requested_post,
    control = ctrl$owner_control,
    catalog_bounds = ctrl$catalog_bounds,
    occupancy = meta_initial$occupancy_object,
    frame_owner = frame_owner,
    declared_qualified = dq_df,
    target_post_times = tpt_df,
    build_design = TRUE, pt_assumption = pt_assumption, pt_start = pt_start,
    # CP08 compact basis: the validated first catalog, occupancy and CSR
    # mappings (sealed, same owners, identical arguments) are reused; only the
    # design step is constructed here.
    metadata_source = meta_initial
  )

  design_all <- metadata_obj$design
  ledger <- ipw_pipeline_ledger_support(ledger, metadata_obj, design_all)
  n_post_bound <- sum(!ledger$available & !is.na(ledger$design_param_index))

  # Reserve the whole integration's known arrays before fitting nuisance views.
  # This is a conservative numeric-payload envelope, using the existing owner
  # budget for the composition rather than granting every producer the full cap.
  # D5 arrays: the full ledger (R_full), the post-contrast basis (at most the
  # number of represented unavailable post vertices) and their IF/enclosure columns.
  # CP08 scaling action 3: exact plan dimensions replace the former 2 J T bound.
  # The Q_aux fitted auxiliary comparisons (union of both families) carry views,
  # nuisances and capture metadata with the M original rows; the reported
  # columns are the Q_rep pretrend cells of the selected family (including
  # reference zeros and unavailable cells) plus every ledger cell. Retained
  # covariance (GMM/replay on M, or on the retained basis L below) is reserved
  # separately from the auxiliary/report buffers: the reported joint covariance
  # (R_rep^2; three boxes under transfer) and one replay chunk of the union.
  M_rows <- as.double(nrow(metadata_obj$catalog$rows))
  M_bound <- M_rows + Q_aux
  R_rep <- Q_rep + R_full
  J_bound <- nrow(metadata_obj$design$parameters) + Q_rep + R_full + n_post_bound +
    length(unique(metadata_obj$design$parameters$g))
  p_bound <- ncol(f_Z)
  transfer_requested <- isTRUE(run_numerical_transfer)
  aux_report_entries <- ipw_owner_safe_sum(c((if (transfer_requested) 64.0 else 16.0)*R_rep*R_rep,
    32.0*ctrl$owner_control$max_chunk_rows*(Q_aux+p_bound+1.0)))
  # CP08 P3: the accepted Prop. 5/Cor. 6 class identities license certifying
  # built-in baseline IPW under optimal weighting from outcome-free objects and
  # never fitting, replaying or reading dropped rows. Identity weighting, custom
  # callbacks, every CP04 backend (OR/DR/selection/imputation/augmented) and
  # unreduced requests keep their accepted paths unchanged. On that route every
  # original row still carries its views, nuisance and D27 metadata (M_bound
  # terms); the covariance/GMM/replay terms are charged on the certified retained
  # basis plus auxiliary rows before any point fit, not on M_bound here.
  early_reduction <- is.null(cp04) && !is.function(custom_est_method) && identical(wt_mode, "optimal") &&
    isTRUE(reduce_structural_redundancy) && identical(structural_reduction_stage, "outcome_free")
  # CP08 O14: the conditional pretest state charges no GMM, reporting or draw buffer.
  cpt_mode <- isTRUE(conditional_pretest)
  # Keep complete full-frame panels on the existing scalar-leaf/batched replay
  # path. Eligibility is decided from sealed, outcome-free posting counts,
  # before any fit; the compact raw backend is reserved for incomplete panels.
  use_raw_ipw <- ipw_pipeline_use_raw_backend(
    cp04, custom_est_method, run_numerical_transfer, cpt_mode,
    environment(sys.function())
  ) && (n_a_units != n_total || any(panel_owner$postings$date_count != n_a_units))
  # CP08 compact canonical basis: on the outcome-free route, the exact
  # complete-roster one-operator/never-treated/PT-post/fixed-X/unit-weight
  # predicate is decided here from the sealed owners and validated catalog,
  # before any nuisance fit or view construction (NULL keeps the accepted
  # view registry route; nothing below ever retries the other route).
  compact_facts <- NULL
  compact_res <- NULL
  M_cov <- if (early_reduction || cpt_mode) 0.0 else M_rows
  combined_initial <- ipw_owner_safe_sum(c(
    owners_live, aux_env$retained_entries,
    4.0*as.double(nrow(data))*(length(data)+p_bound),
    32.0*M_bound*(p_bound*p_bound + p_bound*length(panel_owner$cohorts) + 64.0),
    (if (transfer_requested) 64.0 else 16.0)*M_cov*M_cov,
    if (cpt_mode) 0.0 else (if (transfer_requested) 128.0 else 32.0)*n_total*J_bound,
    if(reserve_draws)(if (transfer_requested) 16.0 else 8.0)*as.double(biters)*J_bound*(length(agg_types)+1.0) else 0.0,
    32.0*ctrl$owner_control$max_chunk_rows*(M_cov+p_bound+1.0),
    if (cpt_mode) 0.0 else aux_report_entries,
    # D27 capture metadata persists across fits; one summary chunk is live.
    if (is.function(custom_est_method)) 0.0 else
      2048.0*M_bound + 2048.0 + 8.0*min(n_a_units,ctrl$owner_control$max_chunk_rows)
  ),min(ctrl$owner_control$max_retained_entries,ctrl$owner_control$max_matrix_entries),
  if (early_reduction) "pipeline owners, original-row views/nuisances, ledger/basis, retained results and inference matrices before the structural basis" else
  "pipeline simultaneous owners, nuisances, ledger/basis, retained results and inference matrices")

  combined_work <- ipw_owner_safe_sum(c(
    or_preprocess_work, aux_env$work_ops,
    (if (transfer_requested) 128.0 else 32.0)*n_total*((if (cpt_mode) 0.0 else J_bound*J_bound) + M_cov*M_cov + p_bound*p_bound),
    if(reserve_draws)(if (transfer_requested) 64.0 else 8.0)*as.double(biters)*n_total*J_bound else 0.0,
    if (is.function(custom_est_method)) 0.0 else M_bound*(128.0*max(1,n_a_units)+4096.0)
  ),ctrl$owner_control$max_work_ops,
  if (early_reduction) "pipeline construction/reporting/aggregation/interval/multiplier work before the structural basis" else
  "pipeline reporting/aggregation/interval/multiplier whole-path work envelope")
  # The existing whole-pipeline envelope above covers the fact scan's N-weight
  # comparison and O(M + T + chunk_rows) date buffers inside its input-copy,
  # per-original-row and bounded-chunk allowances. Its original M-bound scan
  # work dominates these complete-date reads; the compact reservation below
  # also retains its explicit date-fact charge. Run only after both guards.
  compact_facts <- if (!use_raw_ipw && early_reduction && !cpt_mode)
    ipw_pipeline_compact_basis_facts(metadata_obj, frame_owner, ctrl$owner_control, pt_assumption,
      run_numerical_transfer, n_total) else NULL
  aux_record$report_buffer_entries <- aux_report_entries
  aux_record$n_reported_columns <- R_rep
  progress("comparison catalog: ", M_rows, " comparisons; design: ", nrow(design_all$parameters), " parameters, ",
    R_full, " requested post cells (", sum(ledger$available), " estimable)")

  # 4. Build shared nuisance views (with zero-outcome chunked mask verification) and point links
  if (!is.null(cp04)) {
    progress("fitting the full-roster backend on ", M_rows, " catalog comparisons and ", Q_aux, " auxiliary comparisons")
    nuis_res <- build_cp04_public_registry(metadata_obj, frame_owner, cp04_reader, ledger, ctrl, cg_mode, anticipation,
      retained_pipeline_entries = combined_initial, external_work_ops = combined_work, or_method = cp04$or_method,
      automatic_chunk_rows = isTRUE(cp04$automatic_chunk_rows), outcome_recipe = cp04$outcome_recipe,
      selection_model = if (is.null(cp04$selection_model)) "pair_logit" else cp04$selection_model,
      selection_pooling = if (is.null(cp04$selection_pooling)) "separate" else cp04$selection_pooling,
      dr_method = cp04$dr_method, treatment_columns = cp04$treatment_columns,
      imputation = if (imputation) cp04 else NULL, augmented = if (augmented) cp04 else NULL,
      reference_source = cp04$reference_source, ipw_method = cp04$ipw_method,
      custom_est_method = custom_est_method, custom_dots = custom_dots, aux_spec = aux_plan,
      # CP08 descriptor context: the registry decides the early canonical OR/DR
      # intake from these request facts; a standalone caller without them
      # keeps the link-view route.
      pt_assumption = pt_assumption, run_numerical_transfer = run_numerical_transfer)
    combined_initial <- ipw_owner_safe_sum(c(combined_initial,
      nuis_res$cp04$accounting$owned_entries, cp04_reader$accounting$owned_entries),
      min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_matrix_entries), "CP04 combined retained pipeline")
    # CP08/D78 work ownership: the registry composes this pipeline prefix exactly
    # once (cp04_registry_work_accounting). Baseline OR/DR report CUMULATIVE
    # bounds that already contain `combined_work` plus the registry's
    # descriptor/context preparation; adding the prefix again here charged the
    # whole upstream work twice. Every other backend keeps the unchanged legacy
    # composition inside the same helper.
    cp04_work_acc <- nuis_res$cp04_work_accounting
    if (is.null(cp04_work_acc) || !identical(as.double(cp04_work_acc$prefix), as.double(combined_work)))
      cdid_ipw_pipeline_abort("CP04 registry work accounting does not start from the pipeline prefix it was given.",
        c("cdid_ipw_pipeline_integrity_error", "cdid_ipw_owner_error"))
    combined_work <- ipw_owner_safe_sum(cp04_work_acc$composed_total, ctrl$owner_control$max_work_ops,
      "CP04 fitting plus downstream work")
  } else if (is.function(custom_est_method)) {
    progress("evaluating the custom IPW callback on ", M_rows, " catalog comparisons")
    nuis_res <- build_custom_ipw_catalog_links(
      metadata_obj = metadata_obj,
      frame_owner = frame_owner,
      owner_control = ctrl$owner_control,
      point_control = ctrl$point_control,
      custom_est_method = custom_est_method,
      custom_dots = custom_dots,
      retained_pipeline_entries = combined_initial
    )
    combined_initial <- ipw_owner_safe_sum(c(combined_initial, nuis_res$custom_retained_entries),
      min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_matrix_entries),
      "pipeline envelope plus all retained custom callback link arrays")
  } else if (use_raw_ipw) {
    progress("fitting the compact standalone raw-IPW backend on ", M_rows, " catalog comparisons and ", Q_aux, " auxiliary comparisons")
    nuis_res <- build_ipw_raw_backend_registry(
      metadata_obj = metadata_obj,
      frame_owner = frame_owner,
      design = design_all,
      design_link_index = match(seq_len(as.integer(M_rows)), design_all$links$catalog_row),
      owner_control = ctrl$owner_control,
      fitted_control = ctrl$fitted_control,
      point_control = ctrl$point_control,
      glm_control = ctrl$glm_control,
      class_bounds = ctrl$class_bounds,
      weighting = wt_mode,
      reduce_structural_redundancy = reduce_structural_redundancy,
      structural_reduction_stage = structural_reduction_stage,
      aux_spec = aux_plan,
      retained_pipeline_entries = combined_initial,
      external_work_ops = combined_work,
      progress = chained_progress_rows(if (isTRUE(print_details)) progress else NULL,
        "nuisance and link views")
    )
    combined_initial <- ipw_owner_safe_sum(
      c(combined_initial, nuis_res$raw_ipw$accounting$owned_entries),
      min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_matrix_entries),
      "pipeline envelope plus compact raw-IPW backend"
    )
    combined_work <- ipw_owner_safe_sum(
      nuis_res$raw_work_accounting$composed_total,
      ctrl$owner_control$max_work_ops,
      "compact raw-IPW backend construction plus pipeline work"
    )
  } else if (!is.null(compact_facts)) {
    # CP08 compact canonical basis: one nuisance view and fit for the operator,
    # the structural certificate from canonical catalog descriptors and
    # complete-date posting facts, then link views of the retained rows only.
    progress("compact canonical basis: certifying ", M_rows, " comparisons from catalog descriptors before any link view")
    compact_res <- build_pipeline_compact_basis_registry(
      metadata_obj = metadata_obj,
      frame_owner = frame_owner,
      design = design_all,
      design_link_index = match(seq_len(as.integer(M_rows)), design_all$links$catalog_row),
      facts = compact_facts,
      owner_control = ctrl$owner_control,
      fitted_control = ctrl$fitted_control,
      point_control = ctrl$point_control,
      glm_control = ctrl$glm_control,
      class_bounds = ctrl$class_bounds,
      retained_pipeline_entries = combined_initial,
      external_work_ops = combined_work,
      progress = chained_progress_rows(if (isTRUE(print_details)) progress else NULL,
        "retained basis link views")
    )
    nuis_res <- compact_res$registry
    # The registry's actually charged support scans join the composed work.
    combined_work <- ipw_owner_safe_sum(c(combined_work, nuis_res$construction_work$work_ops_charged),
      ctrl$owner_control$max_work_ops, "pipeline work plus charged compact registry scans")
  } else {
    nuis_res <- build_pipeline_nuisance_and_link_views(
      metadata_obj = metadata_obj,
      frame_owner = frame_owner,
      owner_control = ctrl$owner_control,
      fitted_control = ctrl$fitted_control,
      point_control = ctrl$point_control,
      glm_control = ctrl$glm_control,
      retained_pipeline_entries = combined_initial,
      fit_point_links = !early_reduction && !cpt_mode,
      external_work_ops = combined_work,
      progress = chained_progress_rows(if (isTRUE(print_details)) progress else NULL,
        "nuisance and link views")
    )
    # The registry's actually charged mask scans join the composed work.
    combined_work <- ipw_owner_safe_sum(c(combined_work, nuis_res$construction_work$work_ops_charged),
      ctrl$owner_control$max_work_ops, "pipeline work plus charged mask registry scans")
  }
  if (cpt_mode) {
    # CP08 O14: auxiliary links of the fitted union registry; no structural
    # reduction, GMM solve, post-ATT gate, reporting IF or draw.
    reporting_map <- build_ipw_auxiliary_reporting(ledger, panel_owner$all_times, panel_owner$cohorts, anticipation, cg_mode,
      nuis_res, frame_owner, NULL, ctrl, NULL, NULL, FALSE, bp_mode, n_total, a_to_full, combined_initial,
      pt_assumption = pt_assumption, pt_start = pt_start, aux_spec = aux_plan, compute_if = FALSE)
    return(conditional_state("available", NA_character_, nuis_res, reporting_map, aux_plan, metadata_obj,
      combined_initial, combined_work))
  }

  # 5. Certify structural equivalence classes and apply optional redundancy reduction
  #    (only under `weighting == "optimal"`; all rows are kept under `weighting == "identity"`).
  class_cert <- NULL
  batch_facts <- batch_fit_trace <- point_plan <- point_general <- point_custom <- NULL
  if (use_raw_ipw) {
    point_links_all <- nuis_res$point_links
    M_all <- length(point_links_all)
    class_cert <- nuis_res$class_cert
    kept_idx <- nuis_res$kept_idx
    dropped_idx <- nuis_res$dropped_idx
    if (early_reduction) {
      structural_preflight <- nuis_res$structural_preflight
      combined_initial <- combined_initial + structural_preflight$additional_entries
      combined_work <- combined_work + structural_preflight$work_ops_bound
      progress("structural certificate: ", length(kept_idx), " of ", M_all, " comparisons retained (",
        structural_preflight$mask_evidence, ")")
      L_cov <- as.double(length(kept_idx))
      combined_initial <- ipw_owner_safe_sum(c(combined_initial,
        (if (transfer_requested) 64.0 else 16.0)*L_cov*L_cov,
        32.0*ctrl$owner_control$max_chunk_rows*(L_cov+p_bound+1.0)),
        min(ctrl$owner_control$max_retained_entries,ctrl$owner_control$max_matrix_entries),
        "pipeline envelope plus retained structural-basis covariance/GMM/replay buffers, certificate and expansion records")
      combined_work <- ipw_owner_safe_sum(c(combined_work,
        (if (transfer_requested) 128.0 else 32.0)*n_total*L_cov*L_cov),
        ctrl$owner_control$max_work_ops, "pipeline work plus retained structural-basis covariance accumulation")
    }
  } else if (early_reduction) {
    # Outcome-free certificate on every original row (views, verified masks, the
    # one fitted nuisance per view, design vertices/operators), then the retained
    # basis payload, then point fits of retained rows only. Dropped rows keep an
    # outcome-free identity record and are recovered through E.
    M_all <- length(nuis_res$link_views)
    if (!is.null(compact_res)) {
      # Compact canonical basis: the certificate and its simultaneous
      # reservation were determined before any link view was constructed.
      structural_preflight <- compact_res$preflight
      class_cert <- compact_res$certificate
    } else {
      structural_preflight <- preflight_ipw_early_structural_reduction(nuis_res, design_all,
        ctrl$class_bounds, ctrl$owner_control, combined_initial, combined_work)
      class_cert <- certify_ipw_structural_classes_from_views(nuis_res, design_all,
        match(seq_len(M_all), design_all$links$catalog_row), ctrl$class_bounds, ctrl$owner_control,
        retained_pipeline_entries = combined_initial, external_work_ops = combined_work)
    }
    combined_initial <- combined_initial + structural_preflight$additional_entries
    combined_work <- combined_work + structural_preflight$work_ops_bound
    kept_idx <- class_cert$retained_index
    dropped_idx <- class_cert$dropped_index
    progress("structural certificate: ", length(kept_idx), " of ", M_all, " comparisons retained (",
      structural_preflight$mask_evidence, ")")
    # Retained covariance/GMM/replay on the L retained basis rows only; the
    # auxiliary/report buffers were reserved separately with their exact plan.
    L_cov <- as.double(length(kept_idx))
    combined_initial <- ipw_owner_safe_sum(c(combined_initial,
      (if (transfer_requested) 64.0 else 16.0)*L_cov*L_cov,
      32.0*ctrl$owner_control$max_chunk_rows*(L_cov+p_bound+1.0)),
      min(ctrl$owner_control$max_retained_entries,ctrl$owner_control$max_matrix_entries),
      "pipeline envelope plus retained structural-basis covariance/GMM/replay buffers, certificate and expansion records")
    combined_work <- ipw_owner_safe_sum(c(combined_work,
      (if (transfer_requested) 128.0 else 32.0)*n_total*L_cov*L_cov),
      ctrl$owner_control$max_work_ops, "pipeline work plus retained structural-basis covariance accumulation")
    point_links_all <- vector("list", M_all)
    # CP08 stage 2: the trusted batch route is chosen here from actual facts,
    # before any retained point fit; every other request fits link by link.
    # Neither route catches the other's failures.
    batch_facts <- if (identical(pt_assumption, "post") && !isTRUE(run_numerical_transfer))
      ipw_pipeline_batch_eligibility(nuis_res, kept_idx, frame_owner, ctrl$point_control) else NULL
    if (!is.null(batch_facts)) {
      batch_fit_trace <- new.env(parent = emptyenv())
      point_links_all[kept_idx] <- fit_ipw_point_link_view_batch(batch_facts$fitted_view, nuis_res$link_views[kept_idx],
        validate_ipw_point_view_control(ctrl$point_control), date_memo = batch_facts$date_memo,
        trace_env = batch_fit_trace)
    } else {
      point_links_all[kept_idx] <- fit_pipeline_point_links(nuis_res$fitted_views, nuis_res$link_views,
        nuis_res$link_to_nuisance_idx, kept_idx, validate_ipw_point_view_control(ctrl$point_control))
    }
    for (m in dropped_idx) point_links_all[[m]] <- if (!is.null(compact_res))
      ipw_pipeline_compact_reduced_link_record(nuis_res, m, class_cert) else
      ipw_pipeline_reduced_link_record(nuis_res, m, class_cert)
    names(point_links_all) <- names(nuis_res$link_views)
    nuis_res$point_links <- point_links_all
  } else {
    point_links_all <- nuis_res$point_links
    M_all <- length(point_links_all)
    # CP04 retains every actual catalog row. The legacy cross-date mask/forest
    # proof is A-only and does not certify equality of the new selection rows.
    if (is.null(cp04) && identical(wt_mode, "optimal")) {
      class_cert <- certify_ipw_structural_classes_view(
        point_links = point_links_all,
        design = design_all,
        design_link_index = match(seq_len(M_all), design_all$links$catalog_row),
        weighting = "optimal",
        bounds = ctrl$class_bounds,
        owner_control = ctrl$owner_control
      )
    } else if (!is.null(cp04) && identical(wt_mode, "optimal") && isTRUE(reduce_structural_redundancy) &&
               identical(pt_assumption, "post") && !isTRUE(run_numerical_transfer)) {
      class_cert <- cp04_complete_reduction(
        nuis_res$cp04,
        point_links_all,
        design_all,
        match(seq_len(M_all), design_all$links$catalog_row),
        frame_owner,
        ctrl$class_bounds,
        ctrl$owner_control,
        combined_initial,
        combined_work
      )
      if (!is.null(class_cert)) {
        combined_initial <- ipw_owner_safe_sum(
          c(combined_initial, class_cert$accounting$additional_entries),
          min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_matrix_entries),
          "pipeline envelope plus CP04 complete-support structural reduction map"
        )
        combined_work <- ipw_owner_safe_sum(
          c(combined_work, class_cert$accounting$work_ops_bound),
          ctrl$owner_control$max_work_ops,
          "pipeline work plus CP04 complete-support structural reduction"
        )
      }
    }

    if (!is.null(class_cert) && identical(wt_mode, "optimal") && isTRUE(reduce_structural_redundancy)) {
      kept_idx <- class_cert$retained_index
      dropped_idx <- class_cert$dropped_index
    } else {
      kept_idx <- seq_len(M_all)
      dropped_idx <- integer(0L)
    }
  }
  point_links_used <- point_links_all[kept_idx]
  L_used <- length(point_links_used)

  if (L_used == M_all) {
    design_used <- design_all
  } else {
    links_sub <- metadata_obj$design_links[kept_idx, , drop = FALSE]
    rownames(links_sub) <- NULL
    cspec_all <- design_all$calendar_spec
    design_used <- build_calendar_chained_design(
      links = links_sub,
      calendar = metadata_obj$catalog$spec$calendar,
      cohorts = metadata_obj$catalog$spec$cohorts,
      anticipation = metadata_obj$catalog$spec$anticipation,
      target_post_times = tpt_df,
      pt_assumption = if (!is.null(cspec_all) && !is.null(cspec_all$pt_assumption)) cspec_all$pt_assumption else "post",
      pt_start = if (!is.null(cspec_all)) cspec_all$pt_start else NULL
    )
  }

  design_order <- match(metadata_obj$catalog$rows$row_id[kept_idx],design_used$links$catalog_row)
  if (anyNA(design_order) || anyDuplicated(design_order)) {
    cdid_ipw_pipeline_abort("Catalog rows do not map bijectively to structural design rows.", "cdid_ipw_pipeline_integrity_error")
  }
  oo <- order(design_order)
  kept_idx <- kept_idx[oo]; point_links_used <- point_links_used[oo]

  # 6. Ledger alignment with the (possibly reduced) design and the post-contrast basis.
  #    The proved lossless reduction preserves component connectivity; any change
  #    of a cell's catalog-relative status is an integrity failure.
  target_plan <- ipw_pipeline_core_target_plan(ledger, design_used)
  ledger <- target_plan$ledger
  basis <- target_plan$basis
  n_basis <- target_plan$n_basis
  est_params <- target_plan$est_params
  unavailable_rows <- target_plan$unavailable_rows
  core_targets <- target_plan$core_targets
  target_specs <- target_plan$target_specs
  J_est_expected <- target_plan$J_est_expected
  active_backend <- if (!is.null(cp04)) nuis_res$cp04 else if (use_raw_ipw) nuis_res$raw_ipw else NULL
  cp04_used <- if (is.null(active_backend)) NULL else vapply(point_links_used, `[[`, integer(1), "cp04_index")
  cp04_auto_plan <- inherits(active_backend, c("cdid_or_corollary1", "cdid_dr_corollary1", "cdid_dr_baseline",
    "cdid_imputation_backend", "cdid_augmented_backend"))
  # The batch route's one normalization and general producer preflight; the
  # plan then serves both covariance passes and both GMM projection passes.
  if (!is.null(batch_facts)) point_plan <- ipw_point_if_batch_plan(point_links_used, ctrl$point_control,
    batch_facts$fitted_view)
  point_preflight <- if (cp04_auto_plan) {
    NULL # Plan after all simultaneous GMM and output reservations are known.
  } else if (!is.null(active_backend)) {
    cp04_preflight_links(active_backend, cp04_used)
  } else if (is.function(custom_est_method)) {
    if (transfer_requested) {
      preflight_custom_ipw_if_view(point_links_used, ctrl$point_control)
    } else {
      point_custom <- custom_ipw_if_plan(point_links_used, ctrl$point_control)
      point_custom$preflight
    }
  } else if (!is.null(point_plan)) {
    point_plan$preflight
  } else {
    # CP08 missing-panel IPW: the general route's one normalization and
    # whole-replay preflight, sealed with every validated object. The record
    # serves the two covariance passes and the two GMM projection passes below,
    # each re-verifying it at its boundaries; its preflight is the general
    # producer preflight unchanged.
    point_general <- ipw_point_if_general_plan(point_links_used, ctrl$point_control)
    point_general$preflight
  }
  if (!is.null(point_general)) combined_work <- ipw_owner_safe_sum(
    c(combined_work, point_general$general$additional_work_ops_bound),
    ctrl$owner_control$max_work_ops,
    "pipeline cumulative work plus general normalization record construction and verifications")
  if (!is.null(point_custom)) combined_work <- ipw_owner_safe_sum(
    c(combined_work, point_custom$custom$additional_work_ops_bound),
    ctrl$owner_control$max_work_ops,
    "pipeline cumulative work plus custom-IPW IF plan construction and verifications")
  n_dropped_triplets <- if (!is.null(class_cert) && isTRUE(L_used < M_all)) length(class_cert$E$i) else 0L
  gmm_preflight <- preflight_gmm_allocations(L_used, nrow(design_used$parameters), J_est_expected,
    wt_mode, ctrl$gmm_bounds, n_out_rows=replay_n,request_param_cov=TRUE,
    materialize_target_if=TRUE,n_input_rows=replay_n,n_parameters=nrow(design_used$parameters),
    certify_error = !identical(ctrl$gmm_accuracy$method, "conventional_chol_qr"),
    covariance_strategy = ctrl$gmm_accuracy$covariance_strategy,
    n_dropped_triplets = n_dropped_triplets)
  combined_work <- ipw_owner_safe_sum(c(combined_work, gmm_preflight$covariance_strategy_work_ops),
    ctrl$owner_control$max_work_ops, "pipeline composition plus D26 Gram/transform/eigen/factor/solve counted work")
  d5_arrays <- list(n_ledger_cells = R_full, n_represented_cells = sum(!is.na(ledger$core_param_index)),
    n_available_cells = sum(ledger$available), n_unavailable_level_targets = length(unavailable_rows),
    n_basis_contrasts = n_basis, n_core_targets = nrow(core_targets), n_estimable_loading_columns = J_est_expected,
    ledger_if_entries = as.double(n_total) * R_full, basis_if_entries = as.double(n_total) * n_basis,
    ledger_enclosure_entries = if (isTRUE(run_numerical_transfer)) 2.0 * as.double(n_total) * R_full else 0.0,
    basis_enclosure_entries = if (isTRUE(run_numerical_transfer)) 2.0 * as.double(n_total) * n_basis else 0.0,
    envelope_term_entries = (if (run_numerical_transfer) 128.0 else 32.0) * n_total * (R_full + n_post_bound),
    # Retained diagnostic outputs charged to the existing simultaneous envelope:
    # per aggregation type the W (M x R), D (M x K) and point boxes with
    # M, R <= J_bound; seven J_est component vectors; the per-nuisance radius
    # record plus new preflight headers (at most one nuisance per catalog/auxiliary
    # link, 24 entries each), matrix dimension headers (8 per aggregation),
    # and fixed diagnostic/accounting scalars (8 entries).
    retained_diagnostic_entries = if (isTRUE(run_numerical_transfer))
      2.0 * length(agg_types) * J_bound * (J_bound + length(panel_owner$cohorts) + 1.0) +
        7.0 * J_est_expected + 24.0 * M_bound + 8.0 * length(agg_types) + 8.0 else 0.0)
  # Final public diagnostic copying is checked against the recorded pipeline
  # peak. Reserve its existing one-output plus four-scratch charge before the
  # automatic OR/DR row planner consumes the remaining workspace.
  final_diagnostic_copy_entries <- if (cp04_auto_plan)
    5 * (4 * nuis_res$cp04$accounting$owned_entries + 4096) else 0
  # The automatic planner's retained argument (zero on every other CP04 route);
  # the retained moment IF admission below probes the same planner with it.
  cp04_planner_retained <- if (cp04_auto_plan || use_raw_ipw)
    final_diagnostic_copy_entries + combined_initial + gmm_preflight$total_payload_entries +
      d5_arrays$ledger_if_entries + d5_arrays$basis_if_entries +
      d5_arrays$ledger_enclosure_entries + d5_arrays$basis_enclosure_entries +
      d5_arrays$retained_diagnostic_entries else 0
  if (cp04_auto_plan) {
    point_preflight <- cp04_preflight_links(nuis_res$cp04, cp04_used, cp04_planner_retained)
  }
  combined_live <- ipw_owner_safe_sum(c(combined_initial,
    point_preflight$producer_workspace_entries,gmm_preflight$total_payload_entries,
    d5_arrays$ledger_if_entries, d5_arrays$basis_if_entries, d5_arrays$ledger_enclosure_entries, d5_arrays$basis_enclosure_entries,
    d5_arrays$retained_diagnostic_entries, if (is.null(point_plan)) 0 else point_plan$batch$additional_live_entries,
    if (is.null(point_general)) 0 else point_general$general$additional_live_entries,
    if (is.null(point_custom)) 0 else point_custom$custom$additional_live_entries),
    min(ctrl$owner_control$max_retained_entries,ctrl$owner_control$max_matrix_entries),
    "pipeline producer plus GMM consumer plus D5 ledger/basis simultaneous payload")

  # Carry the retained diagnostic reservation into every later simultaneous
  # covariance/auxiliary guard; `combined_live` alone protects only this stage.
  retained_with_diagnostics <- ipw_owner_safe_sum(
    c(combined_initial, d5_arrays$retained_diagnostic_entries),
    min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_matrix_entries),
    "pipeline persistent envelope plus retained diagnostics")

  cp04_retained <- NULL
  cp04_cov_retained <- NULL
  cp04_gmm_retained <- NULL
  cp04_cov_record <- NULL
  cp04_gmm_record <- NULL
  cp04_phase_eligible <- FALSE
  on.exit(if (!is.null(cp04_retained)) cp04_retained_if_release(cp04_retained), add = TRUE)
  on.exit(if (!is.null(cp04_cov_retained)) cp04_retained_if_release(cp04_cov_retained), add = TRUE)
  on.exit(if (!is.null(cp04_gmm_retained)) cp04_retained_if_release(cp04_gmm_retained), add = TRUE)
  if (!is.null(active_backend)) {
    cp04_work <- cp04_stream_plan(active_backend, cp04_used, retained_with_diagnostics)
    aux_ids <- if (!is.null(cp04)) {
      nuis_res$cp04_aux_index[!is.na(nuis_res$cp04_aux_index)]
    } else {
      nuis_res$raw_aux_index[!is.na(nuis_res$raw_aux_index)]
    }
    aux_sel_bp <- if (!is.null(cp04)) nuis_res$cp04_aux$selections[[bp_mode]] else aux_plan$selections[[bp_mode]]
    aux_work <- if (!length(aux_ids)) 0 else if (!cp04_auto_plan && !use_raw_ipw)
      cp04_preflight_links(active_backend, aux_ids)$work_ops_bound else {
      # Match the requested reporting replay's IDs, order and reservation.
      q <- vapply(aux_sel_bp, function(d)
        if (!is.na(d$index) && d$index > 0L) d$index else NA_integer_, integer(1))
      idx_map <- if (!is.null(cp04)) nuis_res$cp04_aux_index else nuis_res$raw_aux_index
      ids <- idx_map[q[!is.na(q)]]; ids <- unique(ids[!is.na(ids)])
      if (!length(ids)) 0 else cp04_preflight_links(active_backend, ids, retained_with_diagnostics +
        8 * n_total * (length(aux_sel_bp) + nrow(ledger)))$work_ops_bound
    }
    combined_work <- ipw_owner_safe_sum(c(combined_work, cp04_work$work,
      2 * point_preflight$work_ops_bound, aux_work), active_backend$control$max_work_ops,
      "CP04 fitting/covariance/GMM/auxiliary cumulative work")
    # CP08 retained moment IF: admitted only after every original charge above,
    # against the complete simultaneous envelope counted once more with the
    # final diagnostic copy and the auxiliary reporting output envelope, and
    # only when neither original row partition would change. Declined, the
    # ordinary stream below runs exactly as before; nothing is allocated here.
    cp04_retained <- cp04_retained_if_reserve(active_backend, cp04_used, point_preflight$chunk_rows, cp04_work,
      retained = retained_with_diagnostics, gmm_retained = cp04_planner_retained,
      # combined_live includes the GMM producer workspace. Covariance may
      # need a wider producer/consumer peak at the same B; charge that excess
      # as well. The store separately reserves all four full Gram buffers.
      live = ipw_owner_safe_sum(c(combined_live, final_diagnostic_copy_entries,
        max(0, ipw_owner_safe_sum(c(cp04_work$peak, active_backend$accounting$owned_entries,
          active_backend$accounting$borrowed_entries)) - point_preflight$producer_workspace_entries),
        8 * n_total * (length(aux_sel_bp) + nrow(ledger)))),
      work = combined_work, control = active_backend$control)
    if (!is.null(cp04_retained)) combined_work <- ipw_owner_safe_sum(c(combined_work, cp04_retained$work),
      active_backend$control$max_work_ops, "CP04 retained moment IF fill and slice work")
    # D92 covariance-only store: only when the whole-lifetime store declined
    # for a one-block, non-transfer, dated augmented or compact raw-IPW backend,
    # independently consider a covariance-only store at `cp04_work$chunk_rows`.
    cp04_phase_eligible <- is.null(cp04_retained) && !transfer_requested &&
      ((inherits(active_backend, "cdid_augmented_backend") &&
        !is.null(active_backend$reference_ref) && !is.null(active_backend$layout)) ||
       inherits(active_backend, "cdid_ipw_raw_backend")) &&
      length(cp04_work$blocks) == 1L
    if (cp04_phase_eligible) {
      cov_phase_inv <- cp04_retained_if_phase_inventory("covariance", active_backend, cp04_work, point_preflight,
        retained_with_diagnostics, cp04_planner_retained, d5_arrays, L_used, M_all)
      cp04_cov_retained <- cp04_retained_if_phase_reserve(active_backend, cp04_used, "covariance",
        cp04_work$chunk_rows, cp04_work, retained = cov_phase_inv$retained, live = cov_phase_inv$live,
        work = combined_work, control = active_backend$control, numerical_transfer = transfer_requested)
      cov_phase_inv <- NULL
      if (!is.null(cp04_cov_retained)) combined_work <- ipw_owner_safe_sum(c(combined_work, cp04_cov_retained$work),
        active_backend$control$max_work_ops, "CP04 covariance retained moment IF fill and slice work")
    }
    # CP08 auxiliary consumer: keep the admitted (still unallocated) store alive
    # through the auxiliary reporting replay only when exact same-backend
    # recipes map requested reporting columns onto retained columns, the frame
    # is the A roster with lambda = 1, no transfer certifies the auxiliary IF,
    # and the complete envelope through the transported target IF, the
    # auxiliary outputs and the diagnostic copies admits the longer lifetime
    # without changing either auxiliary planner's issued partition. Declined,
    # the store keeps its original release after the GMM solve.
    if (!is.null(cp04_retained) && isTRUE(compute_inffunc) && !transfer_requested) {
      q <- vapply(aux_sel_bp, function(d) if (!is.na(d$index) && d$index > 0L) d$index else NA_integer_, integer(1))
      idx_map <- if (!is.null(cp04)) nuis_res$cp04_aux_index else nuis_res$raw_aux_index
      aux_backend_ids <- idx_map[q[!is.na(q)]]
      aux_backend_ids <- aux_backend_ids[!is.na(aux_backend_ids)]
      if (length(aux_backend_ids)) {
        cp04_retained_if_auxiliary_admit(cp04_retained, active_backend, aux_backend_ids,
          retained = retained_with_diagnostics + 8 * n_total * (length(aux_sel_bp) + nrow(ledger)),
          live = cp04_retained$admission_total + 2 * as.double(n_total) * as.double(J_est_expected),
          work = combined_work, control = active_backend$control)
        if (!is.null(cp04_retained$auxiliary)) combined_work <- ipw_owner_safe_sum(
          c(combined_work, cp04_retained$auxiliary$work), active_backend$control$max_work_ops,
          "CP04 retained moment IF auxiliary slice work")
      }
    }
  }
  # CP08/D92 retained point IF (batch, general, and static custom IPW routes):
  # the plan's first covariance replay is kept once and sliced, in the plan's
  # sealed row partition, for the Gram pass and both GMM projection passes.
  # Admitted from the sealed plan's facts and the complete simultaneous
  # envelope (counted once more with the auxiliary reporting output envelope)
  # before any allocation; declined, the plan replays every pass as before
  # (or falls back to D90 sparse replay state on the general route).
  point_carrier <- NULL
  active_point_plan <- if (!is.null(point_plan)) point_plan else if (!is.null(point_custom)) point_custom else point_general
  if (!is.null(active_point_plan)) {
    point_carrier <- ipw_point_if_retained_reserve(active_point_plan, ctrl$point_control,
      live = combined_live + 8 * n_total * (Q_rep + nrow(ledger)), work = combined_work,
      envelope_control = ctrl$owner_control, n_total = n_total, clustered = !is.null(cluster_full),
      numerical_transfer = transfer_requested)
    if (!is.null(point_carrier)) {
      combined_live <- ipw_owner_safe_sum(
        c(combined_live, point_carrier$entries),
        min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_matrix_entries),
        "pipeline producer/GMM envelope plus retained point IF carrier"
      )
      combined_work <- ipw_owner_safe_sum(c(combined_work, point_carrier$work),
        ctrl$owner_control$max_work_ops, "CP08 retained point IF fill and slice work")
      on.exit(if (!is.null(point_carrier)) ipw_point_if_retained_release(point_carrier), add = TRUE)
    }
  }
  active_point_plan <- NULL
  # D90 general route: outcome-free sparse replay-state record bound to the
  # verified general normalization record when the raw IF carrier is declined.
  # Admitted against the complete simultaneous pipeline envelope (including
  # auxiliary reporting outputs) and the point producer/covariance/replay caps
  # before any state allocation or chunk query; typed capacity decline leaves
  # `point_general_state` NULL so the unchanged general replay serves all four passes.
  point_general_state <- NULL
  if (!is.null(point_general) && is.null(point_carrier)) {
    point_general_state_res <- ipw_point_if_general_state_reserve(
      point_general,
      ctrl$point_control,
      live = combined_live + 8 * n_total * (Q_rep + nrow(ledger)),
      work = combined_work,
      envelope_control = ctrl$owner_control,
      .verified_entry = TRUE
    )
    if (!is.null(point_general_state_res)) {
      combined_live <- ipw_owner_safe_sum(
        c(combined_live, point_general_state_res$entries),
        min(ctrl$owner_control$max_retained_entries, ctrl$owner_control$max_matrix_entries),
        "pipeline producer/GMM envelope plus D90 sparse replay state"
      )
      combined_work <- ipw_owner_safe_sum(
        c(combined_work, point_general_state_res$work),
        ctrl$owner_control$max_work_ops,
        "D90 general point-IF sparse replay state construction, slicing and check work"
      )
      point_general_state <- ipw_point_if_general_state_build(
        point_general_state_res,
        point_general,
        ctrl$point_control
      )
      on.exit(if (!is.null(point_general_state)) ipw_point_if_general_state_release(point_general_state), add = TRUE)
    }
    point_general_state_res <- NULL
  }
  # The one producer handed to the covariance stream and both GMM projection
  # passes: the sealed backend, or the backend bound to its retained store.
  cp04_producer <- if (is.null(active_backend)) NULL else if (!is.null(cp04_retained))
    cp04_retained_if_bind(active_backend, cp04_retained) else if (!is.null(cp04_cov_retained))
    cp04_retained_if_bind(active_backend, cp04_cov_retained) else active_backend
  # 7. Stream two-pass centered point moment covariance on `point_links_used`
  progress("streaming the ", L_used, "-moment covariance from complete unit link influence functions")
  cov_stream_used <- if (!is.null(active_backend)) {
    cov_out <- cp04_stream_covariance(cp04_producer, cp04_used, retained_with_diagnostics)
    # D92 covariance-only store: release the store and every reference to it
    # immediately after the covariance Gram pass, before allocating GMM consumers.
    if (!is.null(cp04_cov_retained)) {
      cp04_cov_record <- cp04_retained_if_release(cp04_cov_retained)
      cp04_cov_retained <- NULL
      cp04_producer <- NULL
    }
    cov_out
  } else if (!is.null(point_custom)) {
    stream_custom_ipw_covariance_plan(
      plan = point_custom,
      control = ctrl$point_control,
      carrier = point_carrier
    )
  } else if (is.function(custom_est_method)) {
    stream_custom_ipw_covariance_view(
      point_links = point_links_used,
      control = ctrl$point_control
    )
  } else if (!is.null(point_plan)) {
    stream_ipw_point_covariance_plan(
      plan = point_plan,
      control = ctrl$point_control,
      carrier = point_carrier
    )
  } else if (!is.null(point_general)) {
    # Both passes replay through the verified general record (one
    # normalization for this invocation) and, when admitted, its retained IF
    # carrier or outcome-free sparse replay state; values, order and refusals
    # are those of the standalone stream.
    stream_ipw_point_covariance_general(
      plan = point_general,
      control = ctrl$point_control,
      state = point_general_state,
      carrier = point_carrier
    )
  } else {
    stream_ipw_point_covariance_view(
      point_links = point_links_used,
      control = ctrl$point_control
    )
  }
  mhat_used <- stats::setNames(as.double(cov_stream_used$estimates), design_used$incidence$link_ids)
  V_used <- cov_stream_used$covariance
  dimnames(V_used) <- list(design_used$incidence$link_ids, design_used$incidence$link_ids)

  # Verify `apply_ipw_reduction_map` expansion back to all original links when `class_cert` is present
  expanded_m <- NULL
  max_moment_expansion_diff <- 0.0
  if (!is.null(class_cert)) {
    ret_est <- as.double(cov_stream_used$estimates)[ match(class_cert$retained_index, kept_idx) ]
    expanded_m <- if (inherits(class_cert, "cdid_ipw_structural_certificate")) {
      apply_ipw_reduction_map(class_cert, ret_est)
    } else {
      cp04_expand_reduction_map(class_cert, ret_est)
    }
    if (!early_reduction || !length(dropped_idx)) {
      all_est <- vapply(point_links_all, `[[`, numeric(1L), "estimate")
      max_moment_expansion_diff <- max(abs(all_est - as.double(expanded_m)))
    } else {
      # Dropped rows are never fitted on this route (Prop. 5/Cor. 6), so the
      # fitted-versus-expanded regression self-test has no value to report.
      max_moment_expansion_diff <- NA_real_
    }
  }

  # 8. Solve Chained GMM Core with two-pass IF replay over `point_links_used`.
  #    Every retained D10/D13 moment and its joint covariance enter; component
  #    metadata never licenses blockwise fits. The core still runs its PD,
  #    conditioning and accuracy checks when no requested level is available.
  scale_id_str <- if (is.null(cp04)) "scale:admissible_unit_centered_1_over_nA" else
    "scale:CP04_complete_F_psi_centered_1_over_nA"
  nuis_id_str <- paste0("nuisance_fit:", paste(nuis_res$nuisance_keys, collapse = "+"))

  prov_used <- gmm_covariance_provenance(
    row_level = "unit",
    n_units = replay_n,
    n_clusters = replay_n,
    n_positive_contributing_rows = replay_n,
    centered = TRUE,
    weighted = FALSE,
    weight_type = "unweighted",
    normalization_factor = 1.0 / as.double(n_a_units),
    scale_id = scale_id_str,
    nuisance_fit_id = nuis_id_str,
    reduction = list(
      mode = if (!is.null(cp04) && !is.null(class_cert)) "sparse_class_certificate" else "none",
      class_certificate = class_cert,
      kept_indices = kept_idx,
      n_links_all = M_all,
      reduction_active = (L_used < M_all)
    )
  )
  if (wt_mode == "optimal" && !is.null(ctrl$gmm_accuracy$covariance_strategy)) {
    # This stream was formed here from complete unit IFs. Bind its exact
    # covariance and nuisance identity using the existing snapshot seal;
    # caller metadata never grants the D27 roundoff-negative diagnosis.
    prov_used$complete_unit_if_gram <- ipw_seal(list(covariance = V_used,
      n_units = prov_used$n_units, nuisance_fit_id = prov_used$nuisance_fit_id), "complete_unit_if_gram")
  }
  # D92 GMM-only store: after the covariance store and its temporary buffers
  # have been released, independently consider a GMM-only store at the original
  # GMM partition (`point_preflight$chunk_rows`).
  if (cp04_phase_eligible) {
    gmm_phase_inv <- cp04_retained_if_phase_inventory("gmm", active_backend, cp04_work, point_preflight,
      retained_with_diagnostics, cp04_planner_retained, d5_arrays, L_used, M_all)
    cp04_gmm_retained <- cp04_retained_if_phase_reserve(active_backend, cp04_used, "gmm",
      point_preflight$chunk_rows, cp04_work, retained = gmm_phase_inv$retained, live = gmm_phase_inv$live,
      work = combined_work, control = active_backend$control, numerical_transfer = transfer_requested,
      active_store = cp04_cov_retained)
    gmm_phase_inv <- NULL
    if (!is.null(cp04_gmm_retained)) {
      combined_work <- ipw_owner_safe_sum(c(combined_work, cp04_gmm_retained$work),
        active_backend$control$max_work_ops, "CP04 GMM retained moment IF fill and slice work")
      cp04_producer <- cp04_retained_if_bind(active_backend, cp04_gmm_retained)
    } else {
      cp04_producer <- active_backend
    }
  }

  if_replay_spec <- gmm_if_replay_spec(
    callback = function(consume) {
      if (!is.null(active_backend)) {
        cp04_replay_links(cp04_producer, cp04_used, function(batch) {
          phi <- batch$phi; colnames(phi) <- design_used$incidence$link_ids
          consume(list(unit_id = batch$id, phi = phi, weight = rep(1, nrow(phi)), nuisance_fit_id = nuis_id_str))
        }, chunk_rows = point_preflight$chunk_rows)
        return(invisible(NULL))
      }
      emit <- function(batch) {
        phi_b <- batch$phi
        colnames(phi_b) <- design_used$incidence$link_ids
        consume(list(
          unit_id = batch$id,
          phi = phi_b,
          weight = rep(1.0, nrow(phi_b)),
          nuisance_fit_id = nuis_id_str
        ))
      }
      # Each projection pass regenerates every IF from raw endpoints through the
      # plan, which re-verifies its bindings at the pass boundary and after
      # every consumed batch; with an admitted carrier the pass serves the raw
      # slices of the first covariance replay in the same partition, with the
      # same plan verification at every boundary.
      if (!is.null(point_plan)) {
        if (!is.null(point_carrier)) return(replay_ipw_point_if_retained(point_carrier, point_plan,
          ctrl$point_control, emit))
        return(replay_ipw_point_if_plan(plan = point_plan, control = ctrl$point_control, callback = emit))
      }
      if (!is.null(point_custom)) {
        return(replay_custom_ipw_if_plan(plan = point_custom, control = ctrl$point_control,
          callback = emit, callback_schema = TRUE, carrier = point_carrier))
      }
      # General route: each projection pass replays through the verified
      # general record (same chunk partition, chunk evaluation, four IF terms
      # and callback order as the standalone replay), re-verifying the record
      # at the pass boundary and after every consumed batch.
      if (!is.null(point_general)) {
        return(replay_ipw_point_if_general(plan = point_general, control = ctrl$point_control,
          callback = emit, callback_schema = TRUE, state = point_general_state, carrier = point_carrier))
      }
      replay_fn <- if (is.function(custom_est_method)) replay_custom_ipw_if_view else replay_ipw_point_if_view
      replay_fn(
        point_links = point_links_used,
        control = ctrl$point_control,
        callback_schema = TRUE,
        callback = emit
      )
    },
    input_row_level = "unit",
    aggregate_to_clusters = FALSE,
    n_input_rows = replay_n,
    n_output_rows = replay_n,
    n_positive_contributing_rows = replay_n,
    centered = TRUE,
    weighted = FALSE,
    weight_type = "unweighted",
    normalization_factor = 1.0 / as.double(n_a_units),
    scale_id = scale_id_str,
    nuisance_fit_id = nuis_id_str,
    expected_unit_ids = if (is.null(cp04)) panel_owner$a_id else f_id
  )

  # Pipeline-level refusal context (section 0b): the catalog/class map, the
  # covariance-stream provenance, the design rank record and the core target
  # table that the core and transfer stages cannot see. Scalars, index vectors
  # and small tables already retained by this call; evaluated only on refusal.
  pipeline_refusal_context <- function(stage) {
    force(stage)
    function(e) {
      cat_rows <- metadata_obj$catalog$rows
      list(stage = stage, weighting = wt_mode, control_group = cg_mode, base_period = bp_mode,
        pt_assumption = pt_assumption, pt_start = pt_start,
        n_a_units = n_a_units, n_total = n_total,
        catalog_class_map = list(n_links_all = M_all, n_links_used = L_used, reduction_active = (L_used < M_all),
          reduce_structural_redundancy = reduce_structural_redundancy, kept_indices = kept_idx, dropped_indices = dropped_idx,
          catalog_rows_used = if (is.data.frame(cat_rows)) ipw_pipeline_scalar_state(cat_rows[kept_idx, , drop = FALSE]) else NULL,
          design_link_ids_used = design_used$incidence$link_ids,
          class_certificate = if (is.null(class_cert)) NULL else list(n_links = class_cert$n_links, link_kind = class_cert$link_kind,
            retained_index = class_cert$retained_index, dropped_index = class_cert$dropped_index,
            classes = ipw_pipeline_scalar_state(class_cert$classes))),
        covariance_stream = list(provenance = unclass(prov_used), nuisance_keys = nuis_res$nuisance_keys,
          n_moments = length(mhat_used), link_ids = names(mhat_used)),
        design_rank = design_used$rank_diagnostics,
        core_targets = core_targets,
        ledger_support = ledger[, intersect(c("cell_id", "g", "time", "available", "support_status", "support_reason"), names(ledger)), drop = FALSE],
        resource = list(simultaneous_entries = combined_live, retained_with_diagnostics = retained_with_diagnostics,
          downstream_work_ops = combined_work))
    }
  }

  progress("solving the chained GMM over ", L_used, " moments (weighting = \"", wt_requested, "\")")
  gmm_fit <- ipw_pipeline_resignal_with_state(solve_chained_gmm_core(
    design = design_used,
    mhat = mhat_used,
    V = V_used,
    provenance = prov_used,
    targets = target_specs,
    weighting = wt_requested,
    if_replay = if_replay_spec,
    accuracy = ctrl$gmm_accuracy,
    bounds = ctrl$gmm_bounds,
    gauge_rule = "lowest",
    request_param_cov = TRUE,
    materialize_target_if = TRUE
  ), "cdid_gmm_error", "pipeline_state", pipeline_refusal_context("gmm_core_solve"))
  # CP08 retained moment IF: release the raw storage and every reference to it
  # (store, bound producer, replay closure) now, before the auxiliary replay and
  # the final output, unless the admitted auxiliary extension holds the store
  # for the auxiliary reporting below (released right after it). The bound
  # producer and the replay closure are dropped here in every case.
  if (!is.null(active_backend)) {
    if (!is.null(cp04_gmm_retained)) {
      cp04_gmm_record <- cp04_retained_if_release(cp04_gmm_retained)
      cp04_gmm_retained <- NULL
    }
    if (!is.null(cp04_cov_record) || !is.null(cp04_gmm_record)) {
      point_preflight$retained_moment_if <- cp04_retained_if_phase_record(cp04_cov_record, cp04_gmm_record)
      cp04_cov_record <- NULL
      cp04_gmm_record <- NULL
    } else if (is.null(cp04_retained) || is.null(cp04_retained$auxiliary)) {
      point_preflight$retained_moment_if <- cp04_retained_if_release(cp04_retained)
      cp04_retained <- NULL
    }
    cp04_producer <- NULL
    if_replay_spec <- NULL
  }
  # CP08/D92 retained point IF: the GMM solve was its last consumer.
  point_if_record <- NULL
  if (!is.null(point_plan)) {
    point_if_record <- ipw_point_if_retained_release(point_carrier, point_plan)
    point_carrier <- NULL
  }
  point_custom_record <- if (is.null(point_custom)) NULL else {
    rec <- point_custom$custom
    rec$retained_if <- ipw_point_if_retained_release(point_carrier, point_custom)
    point_carrier <- NULL
    rec
  }
  point_custom <- NULL
  # CP08/D90/D92 general route: the GMM solve was the normalization record's,
  # retained IF carrier's and sparse replay state's last consumer; release
  # them, keep only scalar accounting and drop every reference they held.
  point_general_record <- if (is.null(point_general)) NULL else {
    rec <- point_general$general
    has_carrier <- !is.null(point_carrier)
    rec$retained_if <- ipw_point_if_retained_release(point_carrier, point_general)
    point_carrier <- NULL
    rec$replay_state <- ipw_point_if_general_state_release(
      point_general_state,
      point_general,
      superseded = has_carrier
    )
    rec
  }
  point_general_state <- NULL
  point_general <- NULL
  if_replay_spec <- NULL

  # Separate index maps: full ledger rows -> core target-table rows -> estimable
  # loading columns (the materialized IF / Lambda columns, in core order).
  est_idx <- which(gmm_fit$targets$estimable)
  if (!identical(gmm_fit$targets$target_id, core_targets$target_id) ||
      !identical(as.logical(gmm_fit$targets$estimable), core_targets$kind %in% c("level", "basis"))) {
    cdid_ipw_pipeline_abort("Core target estimability disagrees with the design-derived ledger/basis classification.",
      "cdid_ipw_pipeline_integrity_error", list(core = gmm_fit$targets[, c("target_id", "estimable", "status")], expected = core_targets))
  }
  core_targets$core_target_row <- seq_len(nrow(core_targets))
  core_targets$loading_column <- match(core_targets$core_target_row, est_idx)
  core_targets$status <- gmm_fit$targets$status
  core_targets$support_reason <- gmm_fit$targets$support_reason
  ledger$core_target_row <- NA_integer_; ledger$loading_column <- NA_integer_
  hit <- match(ledger$cell_id, core_targets$target_id)
  ledger$core_target_row[!is.na(hit)] <- core_targets$core_target_row[hit[!is.na(hit)]]
  ledger$loading_column[!is.na(hit)] <- core_targets$loading_column[hit[!is.na(hit)]]
  if (any(ledger$available != !is.na(ledger$loading_column))) {
    cdid_ipw_pipeline_abort("Ledger availability disagrees with core loading columns.", "cdid_ipw_pipeline_integrity_error")
  }
  basis$core_target_row <- match(basis$basis_id, core_targets$target_id)
  basis$loading_column <- core_targets$loading_column[basis$core_target_row]
  basis$available <- !is.na(basis$loading_column)
  basis$support_status <- if (n_basis) core_targets$status[basis$core_target_row] else character(0L)
  basis$support_reason <- if (n_basis) core_targets$support_reason[basis$core_target_row] else character(0L)
  if (n_basis && anyNA(basis$loading_column)) {
    cdid_ipw_pipeline_abort("A post-contrast basis element is not estimable in the core; the component basis construction is inconsistent.", "cdid_ipw_pipeline_integrity_error")
  }

  # 9. Complete Outward Interval Numerical Transfer (Root + Link + Covariance + Theorem 12 GMM)
  num_transfer <- NULL
  if (isTRUE(run_numerical_transfer)) {
    transfer_context <- pipeline_refusal_context("numerical_transfer")
    num_transfer <- ipw_pipeline_resignal_with_state(certify_chained_ipw_numerical_transfer_view(
      point_links = point_links_used,
      design = design_used,
      mhat = as.double(mhat_used),
      V = unname(V_used),
      gmm_fit = gmm_fit,
      weighting = wt_mode,
      root_witness = root_witness,
      certificate_control = ctrl$certificate_control,
      problem_control = ctrl$problem_control,
      link_control = ctrl$link_control,
      covariance_control = ctrl$covariance_control,
      point_control = ctrl$point_control,
      reader_control = ctrl$owner_control,
      gmm_accuracy = ctrl$gmm_accuracy,
      gmm_bounds = ctrl$gmm_bounds,
      retained_pipeline_entries = retained_with_diagnostics
    ), "cdid_ipw_pipeline_error", "pipeline_state", function(e) c(transfer_context(e),
      list(gmm_core = list(target_ids = gmm_fit$targets$target_id, estimable = gmm_fit$targets$estimable,
        status = gmm_fit$targets$status, n_estimable = length(est_idx),
        point_certified = isTRUE(gmm_fit$numerical_guard$point_certified)))))
  }

  # 10. Transport Admissible-Frame Target IFs (`n_a_units x J_est`) to the Full Frame (`n_total x J_est`)
  #     compute_inffunc = FALSE stops IF work here (post-GMM cutoff): the GMM
  #     solve above used the complete internal link IFs/covariance unchanged.
  transport_res <- NULL
  if (isTRUE(compute_inffunc)) {
    target_if_A <- as.matrix(gmm_fit$if_projection$materialized_target_if)
    if (ncol(target_if_A) != length(est_idx)) cdid_ipw_pipeline_abort("Materialized IF column count differs from the estimable target count.", "cdid_ipw_pipeline_integrity_error")
    colnames(target_if_A) <- gmm_fit$targets$target_id[est_idx]
    transport_res <- if (!is.null(cp04)) list(Phi_full = target_if_A * n_total / n_a_units,
      n_total = n_total, n_a_units = n_a_units, complete_frame_rows = TRUE,
      normalization = "core psi on all F rows converted by N/n_A; no zero padding") else transport_if_to_full_frame(
      Phi_A = target_if_A,
      n_total = n_total,
      a_to_full = a_to_full
    )
  } else {
    progress("compute_inffunc = FALSE: point estimates only; no target IF transport, auxiliary IF replay, covariance or inference")
  }

  target_enclosure <- if (is.null(num_transfer)) NULL else
    ipw_pipeline_target_enclosure(num_transfer,gmm_fit,n_total,a_to_full)
  gmm_handoff <- list(fit = gmm_fit, estimates = gmm_fit$targets$estimate[est_idx], if_full = transport_res$Phi_full,
    enclosure = target_enclosure)

  # 11. Direct did-style reporting over the full ledger template and supported-subvector inference
  progress("reporting ", Q_rep, " pretrend and ", R_full, " post cells; fitting ", Q_aux, " auxiliary comparisons")
  reporting_map <- build_ipw_auxiliary_reporting(ledger, metadata_obj$catalog$spec$calendar, metadata_obj$catalog$spec$cohorts,
    metadata_obj$catalog$spec$anticipation, cg_mode, nuis_res, frame_owner, gmm_handoff, ctrl, root_witness, num_transfer,
    isTRUE(run_numerical_transfer), bp_mode, n_total, a_to_full, retained_with_diagnostics,
    custom_est_method = custom_est_method, custom_dots = custom_dots,
    pt_assumption = pt_assumption, pt_start = pt_start, pt_catalog_links = metadata_obj$design_links,
    aux_spec = aux_plan, compute_if = compute_inffunc, retained_store = cp04_retained)
  # CP08 retained moment IF: the admitted auxiliary extension's last consumer
  # was the reporting replay above; release before inference and output.
  if (!is.null(active_backend) && !is.null(cp04_retained)) {
    point_preflight$retained_moment_if <- cp04_retained_if_release(cp04_retained)
    cp04_retained <- NULL
  }
  # Inference-stage refusals (D20/accuracy) carry the Theorem 12 attribution of
  # the already certified transfer (section 0b); evaluated only on refusal.
  attribution_context <- function(stage) {
    force(stage)
    function(e) c(list(stage = stage), ipw_pipeline_transfer_error_attribution(num_transfer, gmm_fit))
  }
  cell_inference <- ipw_pipeline_resignal_with_state(run_cell_inference(reporting_map),
    "cdid_ipw_pipeline_accuracy_error", "error_attribution", attribution_context("reported_cell_inference"))

  # 11b. Private post-contrast basis: estimates, complete IFs, covariance and
  #      certificates separately from level cells, plus cross-covariance with
  #      the supported reported cells (joint analytical covariance).
  basis_result <- NULL
  if (n_basis && !isTRUE(compute_inffunc)) {
    # Point estimates only: no basis IF, covariance or interval.
    b_avail <- !is.na(basis$loading_column)
    b_point <- ipw_pipeline_point_only_inference(stats::setNames(gmm_fit$targets$estimate[basis$core_target_row], basis$basis_id),
      b_avail, n_total, cluster_full, eligible_full, alp)
    basis_result <- list(meta = basis, estimates = b_point$estimates, se = b_point$se, analytical_se = b_point$analytical_se,
      ci_lower = b_point$ci_lower, ci_upper = b_point$ci_upper, cband_crit_val = b_point$cband_crit_val,
      cband_lower = b_point$cband_lower, cband_upper = b_point$cband_upper, vcov = NULL, if_full = NULL,
      cross_vcov_with_reported_cells = NULL, joint_vcov_reported_and_basis = NULL,
      witnesses = gmm_fit$witnesses[basis$basis_id], numerical_certificate = NULL, inference = b_point,
      scope = "Private estimable post-contrast basis (prop:graph(b)); point estimates only (compute_inffunc = FALSE).")
  } else if (n_basis) {
    b_est <- stats::setNames(gmm_fit$targets$estimate[basis$core_target_row], basis$basis_id)
    b_avail <- !is.na(basis$loading_column)
    b_if <- matrix(NA_real_, n_total, n_basis, dimnames = list(NULL, basis$basis_id))
    if (any(b_avail)) {
      b_if[, b_avail] <- transport_res$Phi_full[, basis$loading_column[b_avail], drop = FALSE]
    }
    b_enc <- if (is.null(target_enclosure)) NULL else {
      t_lo <- t_hi <- matrix(NA_real_, 1L, n_basis)
      r_lo <- r_hi <- matrix(NA_real_, n_total, n_basis)
      if (any(b_avail)) {
        cols_b <- basis$loading_column[b_avail]
        t_lo[, b_avail] <- target_enclosure$theta$lo[, cols_b, drop = FALSE]
        t_hi[, b_avail] <- target_enclosure$theta$hi[, cols_b, drop = FALSE]
        r_lo[, b_avail] <- target_enclosure$rows$lo[, cols_b, drop = FALSE]
        r_hi[, b_avail] <- target_enclosure$rows$hi[, cols_b, drop = FALSE]
      }
      list(theta = list(lo = t_lo, hi = t_hi), rows = list(lo = r_lo, hi = r_hi))
    }
    b_inf <- ipw_pipeline_resignal_with_state(ipw_pipeline_supported_inference(b_est, b_if, b_avail, cluster_full, alp, cband, biters, btype, seed,
      ctrl$multiplier_budget, ctrl$inference_accuracy, b_enc, fixed_signs, eligible_full,
      reason = "no post-contrast basis element is supported",
      bstrap = bstrap, caller_rng = caller_rng), "cdid_ipw_pipeline_accuracy_error", "error_attribution", attribution_context("post_basis_inference"))
    av_cells <- which(reporting_map$available)
    av_basis <- which(b_avail)
    cross <- matrix(NA_real_, nrow(reporting_map$target_meta), n_basis, dimnames = list(reporting_map$target_meta$target_id, basis$basis_id))
    joint_vcov <- NULL
    if (length(av_cells) || length(av_basis)) {
      joint_if <- cbind(reporting_map$if_full[, av_cells, drop = FALSE], b_if[, av_basis, drop = FALSE])
      joint_est <- stats::setNames(c(reporting_map$theta_hat[av_cells], b_est[av_basis]), colnames(joint_if))
      # Joint analytical covariance only: never a draw.
      joint <- ipw_pipeline_resignal_with_state(compute_chained_ipw_inference(joint_est, joint_if, cluster_full, alp, FALSE,
        biters, btype, seed, ctrl$multiplier_budget, ctrl$inference_accuracy, NULL, NULL, eligible_full, bstrap = FALSE, caller_rng = FALSE),
        "cdid_ipw_pipeline_accuracy_error", "error_attribution", attribution_context("joint_basis_covariance"))
      joint_vcov <- joint$vcov
      if (length(av_cells) && length(av_basis)) {
        cross[av_cells, av_basis] <- joint$vcov[seq_along(av_cells), length(av_cells) + seq_along(av_basis), drop = FALSE]
      }
    } else {
      joint_vcov <- matrix(NA_real_, 0L, 0L)
    }
    basis_result <- list(meta = basis, estimates = b_inf$estimates, se = b_inf$se, analytical_se = b_inf$analytical_se,
      ci_lower = b_inf$ci_lower, ci_upper = b_inf$ci_upper, cband_crit_val = b_inf$cband_crit_val, cband_lower = b_inf$cband_lower,
      cband_upper = b_inf$cband_upper, vcov = b_inf$vcov, if_full = b_if, cross_vcov_with_reported_cells = cross,
      joint_vcov_reported_and_basis = joint_vcov, witnesses = gmm_fit$witnesses[basis$basis_id],
      numerical_certificate = b_inf$numerical_certificate, inference = b_inf,
      scope = "Private estimable post-contrast basis (prop:graph(b)); internal gauge coordinates are never level estimates; any evaluated contrast composes its IF and bounds explicitly from these columns.")
  }

  # 12. Requested population cohort-share aggregations from the full declared ledger
  progress("inference: ", if (!isTRUE(compute_inffunc)) "not computed (compute_inffunc = FALSE)" else
    if (!is.null(cell_inference$multiplier)) paste0(cell_inference$n_draws, " multiplier draws",
      if (!is.null(ctrl$multiplier_budget$workers) && ctrl$multiplier_budget$workers > 1) paste0(" on ", ctrl$multiplier_budget$workers, " workers") else "")
    else "analytical", "; ", length(agg_types), " fit-time aggregation(s)")
  agg_results <- ipw_pipeline_resignal_with_state(run_aggregations(reporting_map),
    "cdid_ipw_pipeline_accuracy_error", "error_attribution", attribution_context("fit_time_aggregation"))

  n_available <- sum(ledger$available)
  whole_request_complete <- n_available == R_full && all(reporting_map$available)
  if (!is.null(num_transfer)) {
    num_transfer$reported_inference <- cell_inference$numerical_certificate
    num_transfer$basis_inference <- if (is.null(basis_result)) NULL else basis_result$numerical_certificate
    num_transfer$aggregation_inference <- lapply(agg_results,function(z)
      list(cells=z$numerical_certificate,overall=z$overall$numerical_certificate))
    agg_ok <- function(z) {
      # Every available aggregate output must be certified; unavailable outputs carry no certificate.
      cells_ok <- if (any(z$cell_status$available)) isTRUE(z$numerical_certificate$certified) else TRUE
      overall_ok <- if (identical(z$overall$status, "available")) isTRUE(z$overall$numerical_certificate$certified) else TRUE
      cells_ok && overall_ok
    }
    supported_ok <- (any(reporting_map$available) || n_basis > 0L) &&
      (if (any(reporting_map$available)) isTRUE(cell_inference$numerical_certificate$certified) else TRUE) &&
      (is.null(basis_result) || isTRUE(basis_result$numerical_certificate$certified)) &&
      all(vapply(agg_results, agg_ok, logical(1L)))
    agg_complete <- all(vapply(agg_results, function(z) all(z$cell_status$available) && identical(z$overall$status, "available"), logical(1L)))
    num_transfer$gmm_transfer$supported_subvector_certified <- supported_ok
    num_transfer$gmm_transfer$whole_request_complete <- whole_request_complete && agg_complete
    num_transfer$gmm_transfer$certified_target_ids <- c(cell_inference$numerical_certificate$certified_target_ids,
      if (is.null(basis_result)) character(0L) else basis_result$numerical_certificate$certified_target_ids)
    num_transfer$gmm_transfer$unavailable_target_ids <- cell_inference$numerical_certificate$unavailable_target_ids
    families <- list()
    if (any(reporting_map$available)) families$reported_cells <- list(
      certified_target_ids = cell_inference$numerical_certificate$certified_target_ids,
      band_family = cell_inference$numerical_certificate$band_family)
    if (!is.null(basis_result)) families$post_basis <- list(
      certified_target_ids = basis_result$numerical_certificate$certified_target_ids,
      band_family = basis_result$numerical_certificate$band_family)
    num_transfer$gmm_transfer$band_families <- families
    num_transfer$gmm_transfer$band_family <- if (!length(families)) "none" else if (length(families) == 1L) families[[1L]]$band_family else "separate_families_see_band_families"
    num_transfer$gmm_transfer$complete_inference_certified <- supported_ok && whole_request_complete && agg_complete
    num_transfer$gmm_transfer$certificate_scope <- "Certifies exactly `certified_target_ids` under the separate ID sets in `band_families`; separate sup-t families do not certify their union simultaneously. It never affirms unsupported cells or whole-request coverage when `whole_request_complete` is FALSE."
  }
  structure(
    list(
      support_status = if (whole_request_complete) "all_requested_cells_supported" else "partial_support",
      support = list(ledger = ledger, basis = basis, core_targets = core_targets, n_requested = R_full, n_available = n_available,
        n_unavailable = R_full - n_available, whole_request_complete = whole_request_complete,
        unavailable = ledger[!ledger$available, c("cell_id", "g", "time", "support_status", "support_reason"), drop = FALSE],
        population_claim = "none: catalog-relative support only; unavailable cells are not population non-identification claims",
        maps = list(ledger_to_core_target_row = ledger$core_target_row, ledger_to_loading_column = ledger$loading_column,
          basis_to_loading_column = basis$loading_column, core_target_to_loading_column = core_targets$loading_column)),
      resource_accounting = c(list(initial_envelope_entries=combined_initial,
        simultaneous_entries=combined_live,downstream_work_ops=combined_work,point_producer=point_preflight,
        final_diagnostic_copy_entries=final_diagnostic_copy_entries, auxiliary_plan=aux_record,
        mask_registry=if (is.null(nuis_res$construction_work)) NULL else nuis_res$construction_work[
          intersect(c("work_ops_charged", "registry_entries_reserved", "split_records_bound"), names(nuis_res$construction_work))],
        gmm_consumer=gmm_preflight,d5_arrays=d5_arrays,reachable_numeric_payload_only=TRUE,peak_heap_or_RSS_claim=FALSE,
        numerical_transfer_requested = run_numerical_transfer,
        lifetimes="Owners, registry, all original moments, the full ledger/basis arrays and accumulated outputs persist; producer chunks/cache and GMM consumer are reserved concurrently; interval readers/buffers only when requested."),
        # CP08/D78: ownership of the CP04 fit's work composition (prefix charged once).
        if (is.null(nuis_res$cp04_work_accounting)) list() else list(cp04_work_composition = nuis_res$cp04_work_accounting),
        if (early_reduction) list(early_structural = structural_preflight) else list(),
        # CP08 missing-panel IPW: the general route's one normalization record
        # (header and verification charges; the record itself is not retained).
        if (is.null(point_general_record)) list() else list(point_general = point_general_record),
        if (is.null(point_custom_record)) list() else list(point_custom = point_custom_record),
        if (is.null(point_plan)) list() else list(point_batch = list(
          point_fit = list(preflight = batch_fit_trace$preflight, links_fitted = batch_fit_trace$links_fitted,
            endpoints_read = batch_fit_trace$endpoints_read,
            propensity_solver_calls = batch_fit_trace$propensity_solver_calls,
            shared_state_released = batch_fit_trace$shared_state_released),
          if_plan = point_plan$batch,
          retained_if = point_if_record,
          scope = paste0("Batch route: one shared observer state for the retained point fits, released after them; ",
            "one normalization and one sealed plan serve the two covariance and two GMM projection passes, each ",
            "regenerating every IF from raw endpoints in one native call per chunk, or, with an admitted retained ",
            "point IF carrier, the first covariance replay's raw slices in the same partition for the other three.")))),
      panel_owner = panel_owner,
      frame_owner = frame_owner,
      metadata = metadata_obj,
      design_all = design_all,
      design_used = design_used,
      selection = cp04_selection_record(cp04, nuis_res$cp04, cp04_reader),
      nuisance_registry = nuis_res,
      point_links_all = point_links_all,
      point_links_used = point_links_used,
      structural_certificate = class_cert,
      reduction_summary = list(
        weighting = wt_mode,
        weighting_requested = wt_requested,
        reduce_structural_redundancy = reduce_structural_redundancy,
        reduction_active = (L_used < M_all),
        n_links_all = M_all,
        n_links_used = L_used,
        kept_indices = kept_idx,
        dropped_indices = dropped_idx,
        expanded_m = expanded_m,
        max_moment_expansion_diff = max_moment_expansion_diff,
        reduction_stage = if (early_reduction && !is.null(compact_res)) "outcome_free_compact_canonical_basis" else
          if (early_reduction) "outcome_free_before_point_fits" else if (!is.null(class_cert)) "after_point_fits" else "none",
        n_rows_fitted = if (early_reduction) length(kept_idx) else M_all,
        # Heavy link views actually constructed by the main registry (the
        # compact route builds them for the retained basis rows only).
        n_link_views_constructed = if (is.list(nuis_res$construction_work) &&
          !is.null(nuis_res$construction_work$link_view_constructions)) nuis_res$construction_work$link_view_constructions else NA_integer_,
        catalog_constructions = if (isTRUE(metadata_obj$diagnostics$catalog_reuse$reused)) 1L else 2L,
        expansion_self_test = if (is.null(class_cert)) "not_applicable" else if (is.na(max_moment_expansion_diff))
          "not_evaluated: dropped rows are never fitted (Prop. 5/Cor. 6)" else "evaluated_on_fitted_rows"
      ),
      covariance_stream = list(
        mhat = mhat_used,
        V = V_used,
        mean_phi = cov_stream_used$mean_phi,
        n_a_units = n_a_units
      ),
      gmm_fit = gmm_fit,
      covariance_handling = gmm_fit$covariance_handling,
      numerical_transfer = num_transfer,
      auxiliary_reporting = reporting_map$auxiliary,
      full_frame_transport = transport_res,
      reporting = list(
        base_period = bp_mode,
        target_meta = reporting_map$target_meta,
        L_map = reporting_map$L_map,
        estimates = cell_inference$estimates,
        se = cell_inference$se,
        analytical_se = cell_inference$analytical_se,
        numerical_certificate = cell_inference$numerical_certificate,
        z_stat = cell_inference$z_stat,
        p_value = cell_inference$p_value,
        ci_lower = cell_inference$ci_lower,
        ci_upper = cell_inference$ci_upper,
        cband_crit_val = cell_inference$cband_crit_val,
        cband_lower = cell_inference$cband_lower,
        cband_upper = cell_inference$cband_upper,
        vcov = cell_inference$vcov,
        Omega_a = reporting_map$Sigma_a,
        Omega_a_enclosure = reporting_map$Sigma_a_enclosure,
        if_full = reporting_map$if_full,
        available = reporting_map$available,
        support = cell_inference$support,
        inference = cell_inference,
        # Retained reported-cell enclosure (theta/rows boxes): the state a later
        # aggregation request needs to certify aggregates without a refit.
        enclosure = reporting_map$enclosure
      ),
      basis = basis_result,
      aggregations = agg_results,
      controls = controls_out
    ),
    class = "cdid_ipw_pipeline_result"
  )
}
