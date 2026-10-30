# Private CP02 raw weighted-logit arithmetic and local-root certificate component.
#
# Implements the conditional mathematical root certificate of:
#   * docs/theory/cp02_finite_logit_validation.tex (Section 1, Lemma var, Theorem root)
#   * docs/theory/cp02_ipw_gmm_numerical_propagation.tex (Section 1-2, Lemma incl,
#     Lemma exp, Lemma logis, Lemma dual)
#
# Scope and mathematical contract:
#   * Exact finite problem for stored raw A-unit arrays (Z, w, G, pair_observed)
#     and declared ordered cohorts (target first, controls).
#   * Normalizer n = nrow(Z) is the count of ALL eligible A-units (including
#     outside-risk units and pair non-observers), never pair count or weight mass.
#   * Risk set R = {i : G_i in cohorts}, observer set O = {i in R : pair_observed_i},
#     target indicator D_i = 1{G_i == cohorts[1]}.
#   * Exact cohort masses M_j = sum_{i in R, G_i == cohorts[j]} w_i > 0,
#     pair-observer masses Pi_j = sum_{i in O, G_i == cohorts[j]} w_i > 0,
#     marginal rates r_j = Pi_j / M_j, shares q_j = M_j / n, and fitting weights
#     c_i = w_i * M_{j(i)} / Pi_{j(i)}.
#   * All quantities (M_j, Pi_j, r_j, q_j, c_i, eta_i = z_i' beta0, o_i = exp(eta_i),
#     p_i, 1 - p_i, v_i, score g(beta0) = sum_O c_i (p_i - D_i) z_i / n, and
#     Hessian H(beta0) = sum_O c_i v_i z_i z_i' / n) are enclosed outward directly
#     from raw stored weights w_i, never from cached rounded masses or rates.
#   * Validated exponential (Lemma exp): dyadic range reduction y = 2^m u with
#     |u| <= 1/2, Taylor degree N >= 2, remainder bound 2 |u|^(N+1) / (N+1)!
#     evaluated via outward recurrence by exact represented integers k <= N + 1,
#     strictly positive lower endpoint, and m repeated outward squarings.
#     No system libm exp accuracy assumption is made.
#   * Metric is DEFINED by exact M = F' F for explicit nonsingular upper-triangular
#     stored F with positive diagonal, never by a rounded cross-product or inverse.
#   * Outward interval triangular solves enclose F^{-T} H F^{-1} - I, F^{-T} g, and
#     F^{-T} z_i (for i in O).
#   * Returns outward upper bounds delta, t, a, x = upper(a * r),
#     q = upper((delta + x) / (1 - x)), selfmap = upper(t + q * r), and
#     rho = upper(t / (1 - q)), accepting ONLY when a > 0, delta < 1, x < 1, q < 1,
#     and selfmap <= exact stored r.
#   * Every failed sufficient inequality or arithmetic domain/overflow/underflow/
#     zero-crossing refusal is INCONCLUSIVE, never evidence of separation or
#     non-identification. No optimizer, fit, clipping, row dropping, radius search,
#     or public policy is selected here.

cdid_logit_cert_abort <- function(message, class = character(), diagnostics = list()) {
  stop(structure(
    class = unique(c(
      class,
      "cdid_logit_certificate_refusal",
      "cdid_logit_certificate_error",
      "cdid_ipw_certificate_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

logit_cert_arith_refuse <- function(message, reason, diagnostics = list()) {
  stop(structure(
    class = unique(c(
      "cdid_logit_arithmetic_inconclusive",
      "cdid_logit_certificate_inconclusive",
      "cdid_logit_certificate_refusal",
      "cdid_logit_certificate_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, reason = reason, call = NULL, diagnostics = diagnostics)
  ))
}

logit_cert_require_gmm_core <- function() {
  ns_env <- environment(logit_cert_require_gmm_core)
  req_fns <- c(
    "gmm_next", "gmm_iv", "gmm_iv_check", "gmm_iv_add", "gmm_iv_neg",
    "gmm_iv_sub", "gmm_iv_mul", "gmm_iv_div", "gmm_iv_sqrt", "gmm_iv_abs",
    "gmm_iv_dot", "gmm_iv_norm", "gmm_scalar_logical", "gmm_scalar_numeric",
    "gmm_int_bound", "gmm_safe_product", "ldexp_safe"
  )
  for (fn in req_fns) {
    if (!exists(fn, envir = ns_env, mode = "function", inherits = TRUE)) {
      cdid_logit_cert_abort(
        paste0("Required private helper `", fn, "` from `R/chained_gmm_core.R` is not loaded."),
        "cdid_logit_certificate_dependency_error"
      )
    }
  }
  invisible(TRUE)
}

logit_cert_int_bound <- function(x, name, min_val = 1) {
  if (!gmm_scalar_numeric(x) || x != floor(x) || x < min_val || x > 2^53 - 1) {
    cdid_logit_cert_abort(
      paste0("`", name, "` must be an explicit finite integer in [", min_val, ", 2^53 - 1]."),
      c("cdid_logit_certificate_bound_error", "cdid_logit_certificate_control_error")
    )
  }
  as.double(x)
}

logit_cert_safe_product <- function(a, b, max_entries = NULL, label = "allocation") {
  if (!gmm_scalar_numeric(a, nonneg = TRUE) || a != floor(a) ||
      !gmm_scalar_numeric(b, nonneg = TRUE) || b != floor(b)) {
    cdid_logit_cert_abort(
      paste0("Invalid integer dimensions for ", label, "."),
      "cdid_logit_certificate_bound_error"
    )
  }
  da <- as.double(a)
  db <- as.double(b)
  limit <- 2^53 - 1
  if (da > 0 && db > limit / da) {
    cdid_logit_cert_abort(
      paste0("Count product for ", label, " overflows exact IEEE binary64 integer range."),
      "cdid_logit_certificate_bound_error",
      diagnostics = list(label = label, dim_a = da, dim_b = db)
    )
  }
  prod_val <- da * db
  if (!is.finite(prod_val) || prod_val > limit) {
    cdid_logit_cert_abort(
      paste0("Count product for ", label, " exceeds exact integer range."),
      "cdid_logit_certificate_bound_error"
    )
  }
  if (!is.null(max_entries) && prod_val > max_entries) {
    cdid_logit_cert_abort(
      paste0("Count for ", label, " (", prod_val, ") exceeds explicit bound (", max_entries, ")."),
      "cdid_logit_certificate_bound_error",
      diagnostics = list(label = label, requested = prod_val, bound = max_entries)
    )
  }
  prod_val
}

logit_cert_count_output <- function(x) {
  if (x <= .Machine$integer.max) as.integer(x) else as.double(x)
}

logit_cert_safe_sum <- function(vals, max_bound = NULL, label = "sum") {
  limit <- 2^53 - 1
  total <- 0
  for (v in vals) {
    dv <- as.double(v)
    if (!is.finite(dv) || dv < 0 || dv != floor(dv) || dv > limit - total) {
      cdid_logit_cert_abort(
        paste0("Cumulative count for ", label, " overflows exact IEEE binary64 integer range."),
        "cdid_logit_certificate_bound_error",
        diagnostics = list(label = label, total_so_far = total, term = dv)
      )
    }
    total <- total + dv
  }
  if (!is.null(max_bound) && total > max_bound) {
    cdid_logit_cert_abort(
      paste0("Total ", label, " (", total, ") exceeds explicit bound (", max_bound, ")."),
      "cdid_logit_certificate_bound_error",
      diagnostics = list(label = label, requested = total, bound = max_bound)
    )
  }
  total
}

# ---------------------------------------------------------------------------
# Explicit resource and precision controls (NO internal defaults)
# ---------------------------------------------------------------------------

logit_certificate_control <- function(
  taylor_degree,
  max_abs_eta,
  max_squarings,
  max_chunk_rows,
  max_matrix_entries,
  max_work_ops
) {
  logit_cert_require_gmm_core()
  if (missing(taylor_degree) || missing(max_abs_eta) || missing(max_squarings) ||
      missing(max_chunk_rows) || missing(max_matrix_entries) || missing(max_work_ops)) {
    cdid_logit_cert_abort(
      "All fields of `logit_certificate_control` (`taylor_degree`, `max_abs_eta`, `max_squarings`, `max_chunk_rows`, `max_matrix_entries`, `max_work_ops`) must be supplied explicitly by the caller (no internal defaults).",
      c("cdid_logit_certificate_control_error", "cdid_logit_certificate_bound_error")
    )
  }
  N_deg <- logit_cert_int_bound(taylor_degree, "taylor_degree", min_val = 2)
  if (N_deg > 2^53 - 2) {
    cdid_logit_cert_abort(
      "`taylor_degree + 1` must be an exact IEEE binary64 integer (<= 2^53 - 1).",
      c("cdid_logit_certificate_control_error", "cdid_logit_certificate_bound_error")
    )
  }
  if (!gmm_scalar_numeric(max_abs_eta, positive = TRUE)) {
    cdid_logit_cert_abort(
      "`max_abs_eta` must be an explicit finite positive numeric scalar.",
      c("cdid_logit_certificate_control_error", "cdid_logit_certificate_bound_error")
    )
  }
  max_sq <- logit_cert_int_bound(max_squarings, "max_squarings", min_val = 0)
  max_ch <- logit_cert_int_bound(max_chunk_rows, "max_chunk_rows", min_val = 1)
  max_me <- logit_cert_int_bound(max_matrix_entries, "max_matrix_entries", min_val = 1)
  max_wo <- logit_cert_int_bound(max_work_ops, "max_work_ops", min_val = 1)

  structure(
    list(
      taylor_degree = N_deg,
      max_abs_eta = as.double(max_abs_eta),
      max_squarings = max_sq,
      max_chunk_rows = max_ch,
      max_matrix_entries = max_me,
      max_work_ops = max_wo
    ),
    class = "cdid_logit_certificate_control"
  )
}

validate_logit_certificate_control <- function(control) {
  logit_cert_require_gmm_core()
  if (missing(control) || is.null(control) || !is.list(control) || length(control) != 6L ||
      is.null(names(control)) || anyDuplicated(names(control))) {
    cdid_logit_cert_abort(
      "Explicit `control` list of resource and precision bounds is required; no internal defaults exist.",
      c("cdid_logit_certificate_control_error", "cdid_logit_certificate_bound_error")
    )
  }
  req <- c(
    "taylor_degree", "max_abs_eta", "max_squarings",
    "max_chunk_rows", "max_matrix_entries", "max_work_ops"
  )
  if (!setequal(names(control), req)) {
    missing_f <- setdiff(req, names(control))
    extra_f <- setdiff(names(control), req)
    cdid_logit_cert_abort(
      paste0(
        "Invalid `control` fields.",
        if (length(missing_f)) paste0(" Missing: ", paste(missing_f, collapse = ", "), ".") else "",
        if (length(extra_f)) paste0(" Unexpected: ", paste(extra_f, collapse = ", "), ".") else ""
      ),
      c("cdid_logit_certificate_control_error", "cdid_logit_certificate_bound_error")
    )
  }
  logit_certificate_control(
    taylor_degree = control$taylor_degree,
    max_abs_eta = control$max_abs_eta,
    max_squarings = control$max_squarings,
    max_chunk_rows = control$max_chunk_rows,
    max_matrix_entries = control$max_matrix_entries,
    max_work_ops = control$max_work_ops
  )
}

# ---------------------------------------------------------------------------
# Preflight input and resource/work budget validation before allocations
# ---------------------------------------------------------------------------

validate_logit_certificate_cohorts <- function(cohorts) {
  if (missing(cohorts) || is.null(cohorts) || !is.atomic(cohorts) ||
      !is.null(dim(cohorts)) || is.object(cohorts) || !is.numeric(cohorts) ||
      length(cohorts) < 2L || anyNA(cohorts) || any(is.nan(cohorts) | cohorts == -Inf)) {
    cdid_logit_cert_abort(
      "`cohorts` must be an explicit ordered numeric vector `c(target, controls)` of length >= 2 with finite values or `+Inf` for controls.",
      c("cdid_logit_certificate_cohort_error", "cdid_logit_certificate_input_error")
    )
  }
  coh_dbl <- as.double(cohorts)
  if (!is.finite(coh_dbl[1L])) {
    cdid_logit_cert_abort(
      "First element of `cohorts` (target cohort) must be a finite numeric value.",
      c("cdid_logit_certificate_cohort_error", "cdid_logit_certificate_input_error")
    )
  }
  if (anyDuplicated(coh_dbl)) {
    cdid_logit_cert_abort(
      "Declared `cohorts` must be distinct (target first, followed by distinct controls).",
      c("cdid_logit_certificate_cohort_error", "cdid_logit_certificate_input_error")
    )
  }
  coh_dbl
}

validate_logit_certificate_metric_factor <- function(F_mat, p) {
  if (missing(F_mat) || is.null(F_mat) || !is.matrix(F_mat) || is.object(F_mat) ||
      !is.numeric(F_mat) || nrow(F_mat) != p || ncol(F_mat) != p ||
      anyNA(F_mat) || any(!is.finite(F_mat))) {
    cdid_logit_cert_abort(
      "`F` must be an explicit finite numeric p x p matrix.",
      c("cdid_logit_certificate_metric_error", "cdid_logit_certificate_input_error")
    )
  }
  if (!is.null( attr(F_mat, "is_pd") ) || !is.null( attr(F_mat, "certified") )) {
    cdid_logit_cert_abort(
      "Decorative metadata attributes on `F` are refused.",
      c("cdid_logit_certificate_metric_error", "cdid_logit_certificate_input_error")
    )
  }
  F_dbl <- matrix(as.double(F_mat), nrow = p, ncol = p)
  for (j in seq_len(p)) {
    if (F_dbl[j, j] <= 0) {
      cdid_logit_cert_abort(
        "Stored upper-triangular metric factor `F` must have strictly positive diagonal entries.",
        c("cdid_logit_certificate_metric_error", "cdid_logit_certificate_input_error"),
        diagnostics = list(col = j, diag_value = F_dbl[j, j])
      )
    }
    if (j < p) {
      for (i in seq.int(j + 1L, p)) {
        if (F_dbl[i, j] != 0) {
          cdid_logit_cert_abort(
            "Stored metric factor `F` must be strictly upper-triangular (all entries below diagonal must be exactly 0).",
            c("cdid_logit_certificate_metric_error", "cdid_logit_certificate_input_error"),
            diagnostics = list(row = i, col = j, value = F_dbl[i, j])
          )
        }
      }
    }
  }
  F_dbl
}

# ---------------------------------------------------------------------------
# Private automatic-radius request (compact-view path only)
# ---------------------------------------------------------------------------
# One explicit typed request object selects the accepted one-radius rule
#   r = upper(2 t / (1 - delta))            (certified t > 0)
#   r = lower((1 - delta) / (4 a))          (certified exact-zero score, t == 0)
# AFTER the radius-independent enclosure stage; the unchanged Theorem `root`
# inequalities then decide. The request is recognised ONLY by exact identity
# with the canonical object below. A bare `NA_real_`, `NaN`, `NA_integer_`,
# a logical NA, a decorated missing/nonfinite numeric, a modified copy of the
# request or any other invalid radius keeps the original typed refusal and never
# switches modes. Decorated positive numerics retain the existing manual path. The dense entry `certify_chained_logit_root` accepts no
# request (manual positive radius only). Nothing here is exported.
.logit_cert_automatic_radius_mode <- "automatic_one_radius"
.logit_cert_automatic_radius_ops <- 64
.logit_cert_automatic_radius_entries <- 16

logit_cert_automatic_radius_request <- function() {
  structure(
    list(
      request = .logit_cert_automatic_radius_mode,
      formula = "r = upper(2t/(1-delta)); exact-zero score branch r = lower((1-delta)/(4a)); original Theorem root inequalities decide"
    ),
    class = "cdid_logit_automatic_radius_request"
  )
}

logit_cert_is_automatic_radius_request <- function(r) {
  isTRUE(identical(r, logit_cert_automatic_radius_request()))
}

logit_cert_automatic_radius_flag <- function(automatic_radius) {
  if (!gmm_scalar_logical(automatic_radius)) {
    cdid_logit_cert_abort(
      "`automatic_radius` must be a scalar logical.",
      c("cdid_logit_certificate_radius_error", "cdid_logit_certificate_input_error")
    )
  }
  isTRUE(automatic_radius)
}

# One explicit outward radius from the certified upper endpoints delta, t, a.
# No search, refit, retry or alternate solver: a failure to represent a positive
# finite radius is an arithmetic refusal (caller returns INCONCLUSIVE).
# `inside_sufficient_region` is a scalar diagnostic of the proved sufficient
# region a t <= (1-delta)^2 / (2 (3+delta)); it never decides certification.
logit_cert_select_automatic_radius <- function(delta_iv, t_iv, a_iv, a_gt_0, delta_lt_1) {
  delta_upper <- delta_iv[2L]
  t_upper <- t_iv[2L]
  a_upper <- a_iv[2L]
  if (!isTRUE(a_gt_0)) logit_cert_arith_refuse("Automatic radius requires certified a > 0.", reason = "auto_radius_a_not_positive")
  if (!isTRUE(delta_lt_1)) logit_cert_arith_refuse("Automatic radius requires certified delta < 1.", reason = "auto_radius_delta_ge_1")
  if (!is.finite(t_upper) || t_upper < 0.0) logit_cert_arith_refuse("Automatic radius requires finite t >= 0.", reason = "auto_radius_t_invalid")
  d_iv <- logit_b64_sub(gmm_iv(1.0), c(max(0.0, delta_iv[1L]), delta_upper))
  if (d_iv[1L] <= 0.0) logit_cert_arith_refuse("Denominator 1 - delta nonpositive.", reason = "auto_radius_one_minus_delta_nonpositive")
  if (t_upper > 0.0) {
    num_iv <- logit_b64_mul(gmm_iv(2.0), c(max(0.0, t_iv[1L]), t_upper))
    r_sel <- logit_b64_div(num_iv, d_iv)[2L]
    branch <- "positive_score"
    formula <- "r = upper(2t/(1-delta))"
  } else {
    den_iv <- logit_b64_mul(gmm_iv(4.0), c(max(0.0, a_iv[1L]), a_upper))
    if (den_iv[1L] <= 0.0) logit_cert_arith_refuse("Denominator 4a nonpositive.", reason = "auto_radius_4a_nonpositive")
    r_sel <- logit_b64_div(d_iv, den_iv)[1L]
    branch <- "exact_zero_score"
    formula <- "r = lower((1-delta)/(4a))"
  }
  if (!is.finite(r_sel) || r_sel <= 0.0) logit_cert_arith_refuse("Automatic radius is not a positive finite binary64.", reason = "auto_radius_nonpositive")
  at_hi <- logit_b64_mul(gmm_iv(a_upper), gmm_iv(t_upper))[2L]
  bound_lo <- logit_b64_div(
    logit_b64_mul(d_iv, d_iv),
    logit_b64_mul(gmm_iv(2.0), logit_b64_add(gmm_iv(3.0), c(max(0.0, delta_iv[1L]), delta_upper)))
  )[1L]
  list(
    r = r_sel,
    branch = branch,
    formula = formula,
    a_times_t_upper = at_hi,
    sufficient_bound_lower = bound_lo,
    inside_sufficient_region = isTRUE(at_hi <= bound_lo)
  )
}

preflight_logit_certificate <- function(n, p, m, control, automatic_radius = FALSE) {
  auto_r <- logit_cert_automatic_radius_flag(automatic_radius)
  n_d <- logit_cert_int_bound(n, "nrow(Z)", min_val = 1)
  p_d <- logit_cert_int_bound(p, "ncol(Z)", min_val = 1)
  m_d <- logit_cert_int_bound(m, "length(cohorts)", min_val = 2)
  max_e <- control$max_matrix_entries
  max_w <- control$max_work_ops
  B_eff <- min(n_d, control$max_chunk_rows)
  N_deg <- control$taylor_degree
  S_max <- control$max_squarings

  # Conservative concurrent numeric/logical payload envelope, not peak RSS:
  # 4Bp+16B covers slices/coercions/masks/indices;48p^2 covers accumulators,
  # all7 interval matrix routes, F copy, norm temporaries and returned slices;
  # 32p+32m+128 covers vector/cohort/scalar scratch. Object headers are excluded.
  chunk_Z_entries <- logit_cert_safe_product(4, logit_cert_safe_product(B_eff, p_d, max_e, "streamed Z chunk"), max_e, "chunk slice/coercion/validation scratch")
  chunk_vec_entries <- logit_cert_safe_product(16, B_eff, max_e, "chunk vectors, indices, masks and coercion scratch")
  p2 <- logit_cert_safe_product(p_d, p_d, max_e, "p x p matrix block")
  hess_acc_entries <- logit_cert_safe_product(4, p2, max_e, "Hessian interval accumulators (2 x p x p x 2)")
  tri_scratch_entries <- logit_cert_safe_product(44, p2, max_e, "metric/triangular/norm/validation/returned-slice scratch")
  vec_acc_entries <- logit_cert_safe_product(32, p_d, max_e, "score, coefficient and triangular vector scratch")
  coh_acc_entries <- logit_cert_safe_product(32, m_d, max_e, "cohort accumulators, labels and validation scratch")

  # Automatic-radius request only: the selection record (radius, branch,
  # a*t and sufficient-region scalars, interval scratch) is a fixed 16-entry
  # payload and the one outward radius expression a fixed 64-operation charge.
  # Manual positive-radius calls reserve nothing extra (unchanged figures).
  auto_entries <- if (auto_r) .logit_cert_automatic_radius_entries else 0
  auto_ops <- if (auto_r) .logit_cert_automatic_radius_ops else 0
  concurrent_entries <- logit_cert_safe_sum(
    c(chunk_Z_entries, chunk_vec_entries, hess_acc_entries,
      tri_scratch_entries, vec_acc_entries, coh_acc_entries, 128, auto_entries),
    max_bound = max_e,
    label = "concurrent double/interval buffer entries"
  )

  # Conservative interval-operation charge:3n mass additions; per observer
  # 2p^2+9p+10N+2S+23 (two endpoints, each Horner3N+remainder2(N+1));
  # matrix/contraction stage6p^3+6p^2+5p+3m+13. Validation is a separate
  # traversal, with counters returned separately; no exponential work occurs there.
  pass1_ops <- logit_cert_safe_product(3, n_d, max_w, "Pass 1 total/cohort/pair mass interval operations")
  per_obs_ops <- logit_cert_safe_sum(
    c(
      logit_cert_safe_product(2, p2, max_w, "per-observer quadratic ops"),
      logit_cert_safe_product(9, p_d, max_w, "per-observer linear ops"),
      logit_cert_safe_product(10, N_deg, max_w, "two endpoints: Horner3N plus remainder2(N+1)"),
      logit_cert_safe_product(2, S_max, max_w, "per-observer squaring ops"),
      23
    ),
    max_bound = max_w,
    label = "per-observer interval operations"
  )
  pass2_ops <- logit_cert_safe_product(n_d, per_obs_ops, max_w, "Pass 2 worst-case observer interval operations")
  p3 <- logit_cert_safe_product(p2, p_d, max_w, "p^3 triangular solve operations")
  matrix_stage_ops <- logit_cert_safe_sum(
    c(
      logit_cert_safe_product(6, p3, max_w, "matrix triangular solve and F'F operations"),
      logit_cert_safe_product(6, p2, max_w, "matrix normalization and norm operations"),
      logit_cert_safe_product(5, p_d, max_w, "vector normalization and norm operations"),
      logit_cert_safe_product(3, m_d, max_w, "cohort ratio operations"),
      13,
      auto_ops
    ),
    max_bound = max_w,
    label = "post-loop matrix stage operations"
  )
  total_preflight_ops <- logit_cert_safe_sum(
    c(pass1_ops, pass2_ops, matrix_stage_ops),
    max_bound = max_w,
    label = "total preflight interval work operations"
  )

  list(
    n_a_units = n_d,
    p_covariates = p_d,
    m_cohorts = m_d,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(ceiling(n_d / B_eff)),
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    per_observer_ops_bound = per_obs_ops,
    matrix_stage_ops_bound = matrix_stage_ops,
    total_preflight_ops_bound = total_preflight_ops,
    automatic_radius_reserved = auto_r,
    automatic_radius_ops_bound = auto_ops,
    automatic_radius_entries = auto_entries
  )
}

# ---------------------------------------------------------------------------
# Validated outward interval helpers, Lemma exp, Lemma logis, and F^-T solves
# ---------------------------------------------------------------------------

logit_iv_intersect <- function(a, b, label = "interval intersection") {
  lo <- max(a[1L], b[1L])
  hi <- min(a[2L], b[2L])
  if (!is.finite(lo) || !is.finite(hi) || lo > hi) {
    logit_cert_arith_refuse(
      paste0("Empty or nonfinite outward ", label, "."),
      reason = "empty_interval_intersection",
      diagnostics = list(a = a, b = b, label = label)
    )
  }
  c(lo, hi)
}

# Outward binary64 operations that ALWAYS widen to adjacent binary64 neighbours
# (even when an operand is 0), preserving strict binary64 roundoff/subnormal accounting.
logit_b64_add <- function(a, b) {
  s_lo <- a[1L] + b[1L]
  s_hi <- a[2L] + b[2L]
  if (!is.finite(s_lo) || !is.finite(s_hi)) {
    logit_cert_arith_refuse("Interval addition overflowed.", reason = "interval_overflow")
  }
  c(gmm_next(s_lo, FALSE), gmm_next(s_hi, TRUE))
}

logit_b64_sub <- function(a, b) {
  d_lo <- a[1L] - b[2L]
  d_hi <- a[2L] - b[1L]
  if (!is.finite(d_lo) || !is.finite(d_hi)) {
    logit_cert_arith_refuse("Interval subtraction overflowed.", reason = "interval_overflow")
  }
  c(gmm_next(d_lo, FALSE), gmm_next(d_hi, TRUE))
}

logit_b64_mul <- function(a, b) {
  p1 <- a[1L] * b[1L]
  p2 <- a[1L] * b[2L]
  p3 <- a[2L] * b[1L]
  p4 <- a[2L] * b[2L]
  if (!is.finite(p1) || !is.finite(p2) || !is.finite(p3) || !is.finite(p4)) {
    logit_cert_arith_refuse("Interval multiplication overflowed.", reason = "interval_overflow")
  }
  z <- c(p1, p2, p3, p4)
  c(gmm_next(min(z), FALSE), gmm_next(max(z), TRUE))
}

logit_b64_div <- function(a, b) {
  if (b[1L] <= 0 && b[2L] >= 0) {
    logit_cert_arith_refuse("Interval division crosses zero.", reason = "interval_zero_crossing_division")
  }
  q1 <- a[1L] / b[1L]
  q2 <- a[1L] / b[2L]
  q3 <- a[2L] / b[1L]
  q4 <- a[2L] / b[2L]
  if (!is.finite(q1) || !is.finite(q2) || !is.finite(q3) || !is.finite(q4)) {
    logit_cert_arith_refuse("Interval division overflowed.", reason = "interval_overflow")
  }
  z <- c(q1, q2, q3, q4)
  c(gmm_next(min(z), FALSE), gmm_next(max(z), TRUE))
}

# Lemma exp (docs/theory/cp02_ipw_gmm_numerical_propagation.tex):
#   Given scalar binary64 y:
#   1. Require |y| <= max_abs_eta and |y| <= 1024.
#   2. Let m >= 0 be the least integer with |y| * 2^(-m) <= 1/2, and u = y * 2^(-m).
#      Verify exact dyadic scaling: |u| <= 1/2, ldexp_safe(u, m) == y, and
#      for m > 0, 1/4 < |u| <= 1/2 (normal range, zero mantissa loss).
#   3. Evaluate Taylor polynomial T_N(u) = sum_{k=0}^N u^k / k! for N >= 2 via
#      Horner's rule with outward binary64 operations using exact represented k.
#   4. Evaluate remainder epsilon_N = 2 |u|^(N+1) / (N+1)! via outward recurrence
#      R_0 = 2, R_k = (R_{k-1} * |u|) / k for k = 1, ..., N + 1, where each k is
#      verified to be an exact represented binary64 integer (never forming (N+1)!).
#   5. Form [T_N(u) - epsilon_N, T_N(u) + epsilon_N], require positive lower endpoint,
#      and square both endpoints m times with outward binary64 multiplication,
#      requiring finite upper endpoint and strictly positive lower endpoint at every step.
logit_iv_exp_scalar <- function(y, taylor_degree, max_abs_eta, max_squarings) {
  logit_cert_require_gmm_core()
  N <- logit_cert_int_bound(taylor_degree, "taylor_degree", min_val = 2)
  if (N > 2^53 - 2) cdid_logit_cert_abort("Taylor degree plus one exceeds exact integer range.", "cdid_logit_certificate_bound_error")
  if (length(y) != 1L || is.na(y) || !is.numeric(y) || !is.finite(y)) {
    logit_cert_arith_refuse(
      "Nonfinite scalar input to validated exponential.",
      reason = "exp_nonfinite_input",
      diagnostics = list(y = y)
    )
  }
  y_dbl <- as.double(y)
  if (y_dbl == 0.0) y_dbl <- 0.0
  ay <- abs(y_dbl)
  if (ay > max_abs_eta || ay > 1024.0) {
    logit_cert_arith_refuse(
      paste0("Exponential argument |y| = ", signif(ay, 6), " exceeds declared domain bound."),
      reason = "exp_domain_exceeded",
      diagnostics = list(y = y_dbl, abs_y = ay, max_abs_eta = max_abs_eta)
    )
  }

  m <- 0L
  while (ldexp_safe(ay, -m) > 0.5) {
    m <- m + 1L
    if (m > max_squarings || m > 1074L) {
      logit_cert_arith_refuse(
        paste0("Exponential dyadic range reduction required m = ", m, " > max_squarings = ", max_squarings, "."),
        reason = "exp_squarings_exceeded",
        diagnostics = list(y = y_dbl, m = m, max_squarings = max_squarings)
      )
    }
  }

  u <- if (m == 0L) y_dbl else ldexp_safe(y_dbl, -m)
  if (!is.finite(u) || abs(u) > 0.5 || ldexp_safe(u, m) != y_dbl ||
      (m > 0L && abs(u) <= 0.25)) {
    logit_cert_arith_refuse(
      "Exact dyadic range reduction y = 2^m u with |u| <= 1/2 failed verification.",
      reason = "exp_range_reduction_failed",
      diagnostics = list(y = y_dbl, m = m, u = u)
    )
  }

  u_iv <- gmm_iv(u)
  ops_used <- logit_cert_safe_sum(c(logit_cert_safe_product(5, N,
    label = "scalar Taylor operation count"), 4, m), label = "scalar exponential operations")
  one_iv <- gmm_iv(1.0)
  T_iv <- one_iv
  k <- N
  while (k >= 1) {
    k_dbl <- as.double(k)
    if (!is.finite(k_dbl) || k_dbl != floor(k_dbl) || k_dbl < 1 || k_dbl > 2^53 - 1) {
      logit_cert_arith_refuse("Non-exact integer divisor in Taylor evaluation.", reason = "exp_divisor_not_exact")
    }
    prod_iv <- logit_b64_mul(u_iv, T_iv)
    div_iv <- logit_b64_div(prod_iv, gmm_iv(k_dbl))
    T_iv <- logit_b64_add(one_iv, div_iv)
    k <- k - 1
  }

  au_iv <- gmm_iv(abs(u))
  rem_iv <- gmm_iv(2.0)
  k <- 1
  while (k <= N + 1) {
    k_dbl <- as.double(k)
    if (!is.finite(k_dbl) || k_dbl != floor(k_dbl) || k_dbl < 1 || k_dbl > 2^53 - 1) {
      logit_cert_arith_refuse("Non-exact integer divisor in remainder recurrence.", reason = "exp_divisor_not_exact")
    }
    prod_rem <- logit_b64_mul(rem_iv, au_iv)
    rem_iv <- logit_b64_div(prod_rem, gmm_iv(k_dbl))
    rem_iv[1L] <- max(0.0, rem_iv[1L])
    k <- k + 1
  }

  base_lo <- logit_b64_sub(T_iv, rem_iv)[1L]
  base_hi <- logit_b64_add(T_iv, rem_iv)[2L]
  if (!is.finite(base_lo) || !is.finite(base_hi) || base_lo <= 0.0 || base_lo > base_hi) {
    logit_cert_arith_refuse(
      "Validated exponential base interval on [-1/2, 1/2] failed positivity.",
      reason = "exp_nonpositive_base_lower",
      diagnostics = list(y = y_dbl, u = u, base_lo = base_lo, base_hi = base_hi)
    )
  }

  cur_iv <- c(base_lo, base_hi)
  if (m > 0L) {
    for (s in seq_len(m)) {
      cur_iv <- tryCatch(
        logit_b64_mul(cur_iv, cur_iv),
        cdid_logit_arithmetic_inconclusive = function(e) {
          logit_cert_arith_refuse(
            paste0("Exponential squaring overflowed at step ", s, " of ", m, "."),
            reason = "exp_overflow",
            diagnostics = list(y = y_dbl, m = m, step = s)
          )
        }
      )
      if (!is.finite(cur_iv[1L]) || !is.finite(cur_iv[2L]) || cur_iv[1L] <= 0.0 || cur_iv[1L] > cur_iv[2L]) {
        logit_cert_arith_refuse(
          paste0("Exponential lower endpoint underflowed to <= 0 at squaring step ", s, " of ", m, "."),
          reason = "exp_lower_underflow",
          diagnostics = list(y = y_dbl, m = m, step = s, lo = cur_iv[1L], hi = cur_iv[2L])
        )
      }
    }
  }

  list(iv = cur_iv, m = m, u = u, ops = as.double(ops_used))
}

logit_iv_exp <- function(eta_iv, taylor_degree, max_abs_eta, max_squarings) {
  logit_cert_require_gmm_core()
  if (!is.numeric(eta_iv) || length(eta_iv) < 1L || length(eta_iv) > 2L || anyNA(eta_iv)) {
    cdid_logit_cert_abort(
      "`eta_iv` must be a numeric scalar or length-2 interval `c(lo, hi)`.",
      "cdid_logit_certificate_input_error"
    )
  }
  if (length(eta_iv) == 1L) eta_iv <- c(as.double(eta_iv), as.double(eta_iv))
  if (!all(is.finite(eta_iv)) || eta_iv[1L] > eta_iv[2L]) {
    logit_cert_arith_refuse(
      "Invalid or nonfinite linear predictor interval `eta_iv`.",
      reason = "exp_nonfinite_input",
      diagnostics = list(eta_iv = eta_iv)
    )
  }
  res_lo <- logit_iv_exp_scalar(eta_iv[1L], taylor_degree, max_abs_eta, max_squarings)
  if (eta_iv[2L] == eta_iv[1L]) {
    res_hi <- res_lo
    eval_count <- 1L
    total_ops <- res_lo$ops
  } else {
    res_hi <- logit_iv_exp_scalar(eta_iv[2L], taylor_degree, max_abs_eta, max_squarings)
    eval_count <- 2L
    total_ops <- res_lo$ops + res_hi$ops
  }
  out_iv <- c(res_lo$iv[1L], res_hi$iv[2L])
  if (!all(is.finite(out_iv)) || out_iv[1L] <= 0.0 || out_iv[1L] > out_iv[2L]) {
    logit_cert_arith_refuse(
      "Validated exponential interval failed positivity or ordering.",
      reason = "exp_invalid_output",
      diagnostics = list(eta_iv = eta_iv, out_iv = out_iv)
    )
  }
  list(
    iv = out_iv,
    m_lo = res_lo$m,
    m_hi = res_hi$m,
    max_m = max(res_lo$m, res_hi$m),
    eval_count = eval_count,
    ops = total_ops
  )
}

# Lemma logis (docs/theory/cp02_ipw_gmm_numerical_propagation.tex):
#   Given o in [o_lo, o_hi] with o_lo > 0:
#     p in [o_lo / (1 + o_lo), o_hi / (1 + o_hi)]
#     1 - p = 1 / (1 + o) in [1 / (1 + o_hi), 1 / (1 + o_lo)]
#     v = p (1 - p) in [p_lo * (1 - p)_lo, p_hi * (1 - p)_hi]
logit_iv_logistic_from_odds <- function(o_iv) {
  if (!is.numeric(o_iv) || length(o_iv) != 2L || any(!is.finite(o_iv)) ||
      o_iv[1L] <= 0.0 || o_iv[1L] > o_iv[2L]) {
    logit_cert_arith_refuse(
      "Odds interval must be finite with strictly positive lower endpoint.",
      reason = "logistic_nonpositive_odds",
      diagnostics = list(o_iv = o_iv)
    )
  }
  one_iv <- gmm_iv(1.0)
  o_lo_iv <- gmm_iv(o_iv[1L])
  o_hi_iv <- gmm_iv(o_iv[2L])

  den_lo_iv <- gmm_iv_add(one_iv, o_lo_iv)
  den_hi_iv <- gmm_iv_add(one_iv, o_hi_iv)

  p_lo <- gmm_iv_div(o_lo_iv, den_lo_iv)[1L]
  p_hi <- gmm_iv_div(o_hi_iv, den_hi_iv)[2L]
  p_iv <- c(p_lo, p_hi)

  q_lo <- gmm_iv_div(one_iv, den_hi_iv)[1L]
  q_hi <- gmm_iv_div(one_iv, den_lo_iv)[2L]
  one_minus_p_iv <- c(q_lo, q_hi)

  if (p_iv[1L] <= 0.0 || one_minus_p_iv[1L] <= 0.0 || p_iv[1L] > p_iv[2L] ||
      one_minus_p_iv[1L] > one_minus_p_iv[2L]) {
    logit_cert_arith_refuse(
      "Logistic probability or complement lower bound underflowed to <= 0.",
      reason = "logistic_probability_underflow",
      diagnostics = list(o_iv = o_iv, p_iv = p_iv, one_minus_p_iv = one_minus_p_iv)
    )
  }

  v_iv <- gmm_iv_mul(p_iv, one_minus_p_iv)
  if (v_iv[1L] <= 0.0 || v_iv[1L] > v_iv[2L]) {
    logit_cert_arith_refuse(
      "Logistic variance lower bound underflowed to <= 0.",
      reason = "logistic_variance_underflow",
      diagnostics = list(o_iv = o_iv, v_iv = v_iv)
    )
  }

  list(odds = o_iv, prob = p_iv, one_minus_prob = one_minus_p_iv, var = v_iv, ops = 7)
}

# Solve F' y = b for upper-triangular nonsingular F (so F' is lower-triangular)
# using outward interval forward substitution on exact stored entries of F.
logit_iv_forwardsolve_ft <- function(F_mat, b_iv) {
  p <- nrow(F_mat)
  if (is.null(dim(b_iv))) {
    if (length(b_iv) != p) {
      cdid_logit_cert_abort("Dimension mismatch in `logit_iv_forwardsolve_ft`.", "cdid_logit_certificate_input_error")
    }
    b_iv <- cbind(as.double(b_iv), as.double(b_iv))
  }
  if (!is.matrix(b_iv) || nrow(b_iv) != p || ncol(b_iv) != 2L) {
    cdid_logit_cert_abort("`b_iv` must be a p x 2 interval matrix.", "cdid_logit_certificate_input_error")
  }
  y_iv <- matrix(0.0, nrow = p, ncol = 2L)
  ops <- 0
  for (i in seq_len(p)) {
    s_iv <- b_iv[i, ]
    if (i > 1L) {
      for (j in seq_len(i - 1L)) {
        f_ji <- F_mat[j, i]
        if (f_ji != 0.0) {
          s_iv <- gmm_iv_sub(s_iv, gmm_iv_mul(gmm_iv(f_ji), y_iv[j, ]))
          ops <- ops + 2
        }
      }
    }
    y_iv[i, ] <- gmm_iv_div(s_iv, gmm_iv(F_mat[i, i]))
    ops <- ops + 1
  }
  list(y = y_iv, ops = ops)
}

# Solve F x = b for upper-triangular nonsingular F using outward interval back-substitution.
logit_iv_backsolve_f <- function(F_mat, b_iv) {
  p <- nrow(F_mat)
  if (is.null(dim(b_iv))) {
    if (length(b_iv) != p) {
      cdid_logit_cert_abort("Dimension mismatch in `logit_iv_backsolve_f`.", "cdid_logit_certificate_input_error")
    }
    b_iv <- cbind(as.double(b_iv), as.double(b_iv))
  }
  if (!is.matrix(b_iv) || nrow(b_iv) != p || ncol(b_iv) != 2L) {
    cdid_logit_cert_abort("`b_iv` must be a p x 2 interval matrix.", "cdid_logit_certificate_input_error")
  }
  x_iv <- matrix(0.0, nrow = p, ncol = 2L)
  ops <- 0
  for (i in seq.int(p, 1L, by = -1L)) {
    s_iv <- b_iv[i, ]
    if (i < p) {
      for (j in seq.int(i + 1L, p)) {
        f_ij <- F_mat[i, j]
        if (f_ij != 0.0) {
          s_iv <- gmm_iv_sub(s_iv, gmm_iv_mul(gmm_iv(f_ij), x_iv[j, ]))
          ops <- ops + 2
        }
      }
    }
    x_iv[i, ] <- gmm_iv_div(s_iv, gmm_iv(F_mat[i, i]))
    ops <- ops + 1
  }
  list(x = x_iv, ops = ops)
}

# Outward enclosure of K = F^{-T} H F^{-1} - I = F^{-T} (H - F'F) F^{-1},
# where M = F'F is defined by the exact product of stored upper-triangular F.
logit_iv_whitened_hessian_discrepancy <- function(F_mat, H_iv) {
  p <- nrow(F_mat)
  M_iv <- array(0.0, dim = c(p, p, 2L))
  E_iv <- array(0.0, dim = c(p, p, 2L))
  ops <- 0
  solves <- 0L

  for (k in seq_len(p)) {
    for (l in seq.int(k, p)) {
      rmax <- min(k, l)
      m_kl <- gmm_iv_dot(F_mat[seq_len(rmax), k], F_mat[seq_len(rmax), l])
      ops <- ops + 2 * rmax
      if (k == l) m_kl[1L] <- max(0.0, m_kl[1L])
      M_iv[k, l, ] <- m_kl
      M_iv[l, k, ] <- m_kl

      e_kl <- gmm_iv_sub(H_iv[k, l, ], m_kl)
      ops <- ops + 1
      E_iv[k, l, ] <- e_kl
      E_iv[l, k, ] <- e_kl
    }
  }

  # Route 1: Direct interval solve of F^{-T} H F^{-1} - I
  U_A <- array(0.0, dim = c(p, p, 2L))
  for (l in seq_len(p)) {
    sv <- logit_iv_forwardsolve_ft(F_mat, matrix(H_iv[, l, ], nrow = p, ncol = 2L))
    U_A[, l, ] <- sv$y
    ops <- ops + sv$ops
    solves <- solves + 1L
  }
  K_A <- array(0.0, dim = c(p, p, 2L))
  for (k in seq_len(p)) {
    sv <- logit_iv_forwardsolve_ft(F_mat, matrix(U_A[k, , ], nrow = p, ncol = 2L))
    K_A[k, , ] <- sv$y
    ops <- ops + sv$ops
    solves <- solves + 1L
    K_A[k, k, ] <- gmm_iv_sub(K_A[k, k, ], gmm_iv(1.0))
    ops <- ops + 1
  }

  # Route 2: Residual interval solve of F^{-T} (H - F'F) F^{-1}
  U_B <- array(0.0, dim = c(p, p, 2L))
  for (l in seq_len(p)) {
    sv <- logit_iv_forwardsolve_ft(F_mat, matrix(E_iv[, l, ], nrow = p, ncol = 2L))
    U_B[, l, ] <- sv$y
    ops <- ops + sv$ops
    solves <- solves + 1L
  }
  K_B <- array(0.0, dim = c(p, p, 2L))
  for (k in seq_len(p)) {
    sv <- logit_iv_forwardsolve_ft(F_mat, matrix(U_B[k, , ], nrow = p, ncol = 2L))
    K_B[k, , ] <- sv$y
    ops <- ops + sv$ops
    solves <- solves + 1L
  }

  # Intersect both sound enclosures and enforce exact symmetry K_{kl} = K_{lk} (Lemma incl)
  K_iv <- array(0.0, dim = c(p, p, 2L))
  for (k in seq_len(p)) {
    for (l in seq.int(k, p)) {
      iv_kl <- logit_iv_intersect(K_A[k, l, ], K_B[k, l, ], "K = F^{-T} H F^{-1} - I")
      if (k != l) {
        iv_lk <- logit_iv_intersect(K_A[l, k, ], K_B[l, k, ], "K' = F^{-T} H F^{-1} - I")
        iv_kl <- logit_iv_intersect(iv_kl, iv_lk, "symmetric K_{kl}")
      }
      K_iv[k, l, ] <- iv_kl
      K_iv[l, k, ] <- iv_kl
    }
  }

  delta_iv <- gmm_iv_norm(K_iv[, , 1L], K_iv[, , 2L])
  ops <- ops + 3 * p * p + 2

  list(
    M_iv = M_iv,
    E_iv = E_iv,
    K_direct_iv = K_A,
    K_residual_iv = K_B,
    K_iv = K_iv,
    delta_iv = delta_iv,
    ops = ops,
    solves = solves
  )
}

# ---------------------------------------------------------------------------
# Internal sealed access context seam for dense arrays and compact owner views
# ---------------------------------------------------------------------------

.logit_access_seal_token <- local({
  e <- new.env(parent = emptyenv())
  e$id <- "cdid_logit_access_seal_v1"
  lockEnvironment(e, bindings = TRUE)
  e
})

logit_seal_access_context <- function(ctx, mode) {
  ref <- new.env(parent = emptyenv())
  ref$mode <- as.character(mode)
  ref$seal_token <- .logit_access_seal_token
  ctx$ref <- ref
  ctx <- structure(ctx, class = "cdid_logit_access_context")
  ref$context <- ctx
  lockEnvironment(ref, bindings = TRUE)
  ctx
}

logit_validate_access_context <- function(ctx) {
  if (missing(ctx) || is.null(ctx) || !is.list(ctx) ||
      !inherits(ctx, "cdid_logit_access_context")) {
    cdid_logit_cert_abort(
      "Internal logit access context is missing, unsealed, or malformed.",
      c("cdid_logit_certificate_access_context_error", "cdid_logit_certificate_integrity_error", "cdid_logit_certificate_input_error")
    )
  }
  ref <- ctx$ref
  if (!is.environment(ref) || !environmentIsLocked(ref) ||
      !exists("seal_token", envir = ref, inherits = FALSE) ||
      !exists("context", envir = ref, inherits = FALSE) ||
      !identical(get("seal_token", envir = ref, inherits = FALSE), .logit_access_seal_token)) {
    cdid_logit_cert_abort(
      "Foreign or forged logit access context is refused; only internally sealed contexts may enter shared workers.",
      c("cdid_logit_certificate_access_context_error", "cdid_logit_certificate_integrity_error", "cdid_logit_certificate_input_error")
    )
  }
  snap <- get("context", envir = ref, inherits = FALSE)
  if (!is.list(snap) || !identical(names(ctx), names(snap)) ||
      !all(vapply(names(ctx), function(nm) ((is.null(ctx[[nm]]) && is.null(snap[[nm]])) || gmm_sexp_aliases(ctx[[nm]], snap[[nm]])), logical(1)))) {
    cdid_logit_cert_abort(
      "Logit access context was modified after internal sealing; execution is refused.",
      c("cdid_logit_certificate_access_context_error", "cdid_logit_certificate_integrity_error", "cdid_logit_certificate_input_error")
    )
  }
  invisible(ctx)
}

preflight_logit_certificate_view <- function(
  nuisance_view,
  control,
  reader_control,
  expected_a_id = NULL,
  automatic_radius = FALSE
) {
  logit_cert_require_gmm_core()
  auto_r <- logit_cert_automatic_radius_flag(automatic_radius)
  ctrl <- validate_logit_certificate_control(control)
  if (is.null(expected_a_id) && !missing(reader_control) && !is.null(reader_control)) {
    expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  }
  r_ctrl <- validate_ipw_owner_control(reader_control)
  validate_ipw_nuisance_view(nuisance_view)

  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = nuisance_view$panel_ref, inherits = FALSE)

  n_d <- logit_cert_int_bound(frame_owner$n_a_units, "n_a_units", min_val = 1)
  p_d <- logit_cert_int_bound(frame_owner$n_covariates, "p_covariates", min_val = 1)
  m_d <- logit_cert_int_bound(nuisance_view$cohorts$m_cohorts, "m_cohorts", min_val = 2)

  B_eff <- min(n_d, ctrl$max_chunk_rows, r_ctrl$max_chunk_rows)
  n_chunks <- ceiling(n_d / B_eff)

  ctrl_eff <- ctrl
  ctrl_eff$max_chunk_rows <- B_eff
  arith_pre <- preflight_logit_certificate(n = n_d, p = p_d, m = m_d, control = ctrl_eff, automatic_radius = auto_r)

  r_pre_design_base <- preflight_ipw_range_query(
    b_len = B_eff,
    p_covariates = p_d,
    m_cohorts = m_d,
    control = r_ctrl,
    read_outcomes = FALSE,
    include_design = TRUE,
    enclose_differences = FALSE
  )
  r_pre_nodesign_base <- preflight_ipw_range_query(
    b_len = B_eff,
    p_covariates = p_d,
    m_cohorts = m_d,
    control = r_ctrl,
    read_outcomes = FALSE,
    include_design = FALSE,
    enclose_differences = FALSE
  )

  r_pre_design <- r_pre_design_base
  r_pre_nodesign <- r_pre_nodesign_base
  id_pass_work_ops <- 0
  id_max_chunk_scratch <- 0
  id_byte_summary <- NULL

  if (!is.null(expected_a_id)) {
    # No expected-ID content scan or slice before all base stage caps pass.
    preflight_logit_certificate_view(nuisance_view, ctrl, r_ctrl, expected_a_id = NULL, automatic_radius = auto_r)
    ipw_owner_refuse_decorative_attrs(expected_a_id, "expected_a_id")
    if (is.factor(expected_a_id) || is.object(expected_a_id) || !is.atomic(expected_a_id) ||
        !is.null(dim(expected_a_id)) || length(expected_a_id) != n_d) {
      cdid_ipw_owner_abort(
        "`expected_a_id` must be NULL or a plain atomic vector of length `n_a_units`.",
        c("cdid_ipw_owner_id_error", "cdid_ipw_owner_order_mismatch_error")
      )
    }
    if (identical(panel_owner$id_type, "character") && is.character(expected_a_id) && !anyNA(expected_a_id)) {
      max_chunk_b_id <- 0
      max_single_b_id <- 0
      total_pass_b_id <- 0
      c_start <- 1L
      while (c_start <= n_d) {
        c_end <- as.integer(min(as.double(n_d), as.double(c_start) + B_eff - 1))
        c_seq <- seq.int(c_start, c_end)
        eb <- ipw_owner_id_byte_metadata(expected_a_id[c_seq])
        ob <- ipw_owner_id_byte_metadata(panel_owner$a_id_key[c_seq])
        ch_tot <- ipw_owner_safe_sum(c(eb$total, ob$total), label = "chunk paired expected ID bytes")
        max_chunk_b_id <- max(max_chunk_b_id, ch_tot)
        max_single_b_id <- max(max_single_b_id, eb$max, ob$max)
        total_pass_b_id <- ipw_owner_safe_sum(c(total_pass_b_id, ch_tot), label = "pass paired expected ID bytes")
        c_start <- as.double(c_end) + 1
      }
      r_pre_design <- ipw_owner_preflight_id_bytes(
        r_pre_design_base, r_ctrl,
        list(total = max_chunk_b_id, max = max_single_b_id),
        owner_max_bytes = 0, lookup_count = 0, search_width = 1, simultaneous = TRUE
      )
      r_pre_nodesign <- ipw_owner_preflight_id_bytes(
        r_pre_nodesign_base, r_ctrl,
        list(total = max_chunk_b_id, max = max_single_b_id),
        owner_max_bytes = 0, lookup_count = 0, search_width = 1, simultaneous = TRUE
      )
      id_max_chunk_scratch <- r_pre_design$id_byte_preflight$additional_scratch_entries
      id_pass_work_ops <- ipw_owner_safe_product(
        128, total_pass_b_id, r_ctrl$max_work_ops,
        "whole-pass character ID comparison work"
      )
      id_byte_summary <- list(
        max_chunk_source_bytes = max_chunk_b_id,
        total_pass_source_bytes = total_pass_b_id,
        max_chunk_additional_scratch_entries = id_max_chunk_scratch,
        per_pass_additional_work_ops = id_pass_work_ops
      )
    }
  }

  # Three reader passes in certify_chained_logit_root_view:
  # Pass 1a (validation, include_design = TRUE),
  # Pass 1b (mass accumulation, include_design = FALSE),
  # Pass 2 (observer score/Hessian, include_design = TRUE).
  pass1a_reader_ops <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, r_pre_design_base$work_ops_bound, r_ctrl$max_work_ops, "Pass 1a reader chunk work"),
      id_pass_work_ops
    ),
    max_bound = r_ctrl$max_work_ops,
    label = "Pass 1a reader work operations"
  )
  pass1b_reader_ops <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, r_pre_nodesign_base$work_ops_bound, r_ctrl$max_work_ops, "Pass 1b reader chunk work"),
      id_pass_work_ops
    ),
    max_bound = r_ctrl$max_work_ops,
    label = "Pass 1b reader work operations"
  )
  pass2_reader_ops <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, r_pre_design_base$work_ops_bound, r_ctrl$max_work_ops, "Pass 2 reader chunk work"),
      id_pass_work_ops
    ),
    max_bound = r_ctrl$max_work_ops,
    label = "Pass 2 reader work operations"
  )
  total_reader_work_ops <- ipw_owner_safe_sum(
    c(pass1a_reader_ops, pass1b_reader_ops, pass2_reader_ops),
    max_bound = r_ctrl$max_work_ops,
    label = "total 3-pass reader work operations"
  )

  # Simultaneous stage payload sums across Stages 1a, 1b, 2, and Matrix stage:
  stage1a_concurrent <- logit_cert_safe_sum(
    c(arith_pre$concurrent_entries, r_pre_design$concurrent_entries),
    max_bound = ctrl$max_matrix_entries,
    label = "Pass 1a simultaneous arithmetic + reader concurrent payload entries"
  )
  stage1b_concurrent <- logit_cert_safe_sum(
    c(arith_pre$concurrent_entries, r_pre_nodesign$concurrent_entries),
    max_bound = ctrl$max_matrix_entries,
    label = "Pass 1b simultaneous arithmetic + reader concurrent payload entries"
  )
  stage2_concurrent <- logit_cert_safe_sum(
    c(arith_pre$concurrent_entries, r_pre_design$concurrent_entries),
    max_bound = ctrl$max_matrix_entries,
    label = "Pass 2 simultaneous arithmetic + reader concurrent payload entries"
  )
  stage_matrix_concurrent <- arith_pre$concurrent_entries
  combined_concurrent_entries <- max(
    stage1a_concurrent,
    stage1b_concurrent,
    stage2_concurrent,
    stage_matrix_concurrent
  )

  combined_preflight_ops <- logit_cert_safe_sum(
    c(arith_pre$total_preflight_ops_bound, total_reader_work_ops),
    max_bound = ctrl$max_work_ops,
    label = "combined arithmetic and reader preflight work operations"
  )

  list(
    n_a_units = n_d,
    p_covariates = p_d,
    m_cohorts = m_d,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(n_chunks),
    n_reader_passes = 3L,
    concurrent_entries = arith_pre$concurrent_entries,
    concurrent_bytes_8double = arith_pre$concurrent_bytes_8double,
    reader_concurrent_entries = r_pre_design$concurrent_entries,
    combined_concurrent_entries = combined_concurrent_entries,
    combined_concurrent_bytes_8double = 8 * combined_concurrent_entries,
    stage_concurrent_entries = list(
      pass1a = stage1a_concurrent,
      pass1b = stage1b_concurrent,
      pass2 = stage2_concurrent,
      matrix_stage = stage_matrix_concurrent
    ),
    per_observer_ops_bound = arith_pre$per_observer_ops_bound,
    matrix_stage_ops_bound = arith_pre$matrix_stage_ops_bound,
    total_preflight_ops_bound = arith_pre$total_preflight_ops_bound,
    automatic_radius_reserved = arith_pre$automatic_radius_reserved,
    automatic_radius_ops_bound = arith_pre$automatic_radius_ops_bound,
    automatic_radius_entries = arith_pre$automatic_radius_entries,
    reader_pass_work_ops_bounds = list(
      pass1a = pass1a_reader_ops,
      pass1b = pass1b_reader_ops,
      pass2 = pass2_reader_ops
    ),
    total_reader_work_ops_bound = total_reader_work_ops,
    reader_work_ops_bound = total_reader_work_ops,
    combined_preflight_ops_bound = combined_preflight_ops,
    combined_work_ops_bound = combined_preflight_ops,
    id_comparison_bytes_per_pass = if (!is.null(id_byte_summary)) id_byte_summary$total_pass_source_bytes else 0,
    id_byte_preflight = id_byte_summary,
    endpoints_read = 0L
  )
}

construct_dense_logit_access_context <- function(
  Z,
  w,
  G,
  pair_observed,
  coh_dbl,
  cohort_labels,
  beta0_dbl,
  F_dbl,
  r_dbl,
  ctrl,
  preflight
) {
  n <- nrow(Z)
  p <- ncol(Z)
  m <- length(coh_dbl)
  header <- list(
    mode = "dense",
    n = as.integer(n),
    p = as.integer(p),
    m = as.integer(m),
    z_colnames = colnames(Z),
    coh_dbl = coh_dbl,
    cohort_labels = cohort_labels,
    beta0_dbl = beta0_dbl,
    F_dbl = F_dbl,
    r_dbl = r_dbl,
    automatic_radius = FALSE,
    ctrl = ctrl,
    reader_ctrl = NULL,
    preflight = preflight,
    step = preflight$effective_chunk_rows
  )

  validate_chunk <- function(a_start, a_end) {
    idx <- seq.int(a_start, a_end)
    Z_ch <- Z[idx, , drop = FALSE]
    if (anyNA(Z_ch) || any(!is.finite(Z_ch))) {
      cdid_logit_cert_abort(
        "`Z` contains NA, NaN, or infinite values.",
        c("cdid_logit_certificate_covariate_error", "cdid_logit_certificate_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    w_ch <- w[idx]
    if (anyNA(w_ch) || any(!is.finite(w_ch) | w_ch <= 0)) {
      cdid_logit_cert_abort(
        "`w` must contain finite, strictly positive weights (`w > 0`).",
        c("cdid_logit_certificate_weight_error", "cdid_logit_certificate_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    G_ch <- G[idx]
    if (anyNA(G_ch) || any(is.nan(G_ch) | G_ch == -Inf)) {
      cdid_logit_cert_abort(
        "`G` must contain finite numeric cohort values or `+Inf` (never NA, NaN, or `-Inf`).",
        c("cdid_logit_certificate_cohort_error", "cdid_logit_certificate_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    po_ch <- pair_observed[idx]
    if (anyNA(po_ch)) {
      cdid_logit_cert_abort(
        "`pair_observed` cannot contain NA values.",
        c("cdid_logit_certificate_pair_error", "cdid_logit_certificate_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    j_ch <- match(G_ch, coh_dbl)
    list(
      j_ch = j_ch,
      po_ch = po_ch,
      outside_risk_pair_hit_count = NA_integer_,
      reader_ops = 0,
      n_endpoints_read = 0L
    )
  }

  read_mass_chunk <- function(a_start, a_end) {
    idx <- seq.int(a_start, a_end)
    b_len <- length(idx)
    w_ch <- as.double(w[idx])
    G_ch <- as.double(G[idx])
    po_ch <- as.logical(pair_observed[idx])
    j_ch <- match(G_ch, coh_dbl)
    list(
      b_len = b_len,
      w_ch = w_ch,
      j_ch = j_ch,
      po_ch = po_ch,
      reader_ops = 0,
      n_endpoints_read = 0L
    )
  }

  read_observer_chunk <- function(a_start, a_end) {
    idx <- seq.int(a_start, a_end)
    b_len <- length(idx)
    Z_ch <- matrix(as.double(Z[idx, , drop = FALSE]), nrow = b_len, ncol = p)
    w_ch <- as.double(w[idx])
    G_ch <- as.double(G[idx])
    po_ch <- as.logical(pair_observed[idx])
    j_ch <- match(G_ch, coh_dbl)
    list(
      idx = idx,
      b_len = b_len,
      Z_ch = Z_ch,
      w_ch = w_ch,
      j_ch = j_ch,
      po_ch = po_ch,
      reader_ops = 0,
      n_endpoints_read = 0L
    )
  }

  logit_seal_access_context(
    list(
      mode = "dense",
      header = header,
      validate_chunk = validate_chunk,
      read_mass_chunk = read_mass_chunk,
      read_observer_chunk = read_observer_chunk
    ),
    mode = "dense"
  )
}

construct_view_logit_access_context <- function(
  nuisance_view,
  coh_dbl,
  cohort_labels,
  beta0_dbl,
  F_dbl,
  r_dbl,
  ctrl,
  reader_ctrl,
  preflight,
  expected_a_id = NULL,
  trace_env = NULL,
  automatic_radius = FALSE
) {
  auto_r <- logit_cert_automatic_radius_flag(automatic_radius)
  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  n <- as.integer(frame_owner$n_a_units)
  p <- as.integer(frame_owner$n_covariates)
  m <- length(coh_dbl)
  if (!is.null(trace_env)) {
    if (!is.environment(trace_env)) {
      cdid_ipw_owner_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_owner_input_error")
    }
    trace_env$calls <- 0L
    trace_env$pos_s_index <- integer(0L)
    trace_env$pos_t_index <- integer(0L)
  }

  header <- list(
    mode = "compact_view",
    n = n,
    p = p,
    m = m,
    z_colnames = colnames(frame_owner$Z_a),
    coh_dbl = coh_dbl,
    cohort_labels = cohort_labels,
    beta0_dbl = beta0_dbl,
    F_dbl = F_dbl,
    r_dbl = r_dbl,
    automatic_radius = auto_r,
    ctrl = ctrl,
    reader_ctrl = reader_ctrl,
    preflight = preflight,
    step = preflight$effective_chunk_rows,
    nuisance_ref = nuisance_view$ref,
    frame_ref = nuisance_view$frame_ref,
    panel_ref = nuisance_view$panel_ref,
    defining_pair = nuisance_view$defining_pair
  )

  query_slice <- function(a_start, a_end, include_design) {
    exp_slice <- if (is.null(expected_a_id)) {
      NULL
    } else {
      expected_a_id[seq.int(a_start, a_end)]
    }
    ch_tr <- if (!is.null(trace_env)) new.env(parent = emptyenv()) else NULL
    ch <- query_ipw_view_chunk(
      view = nuisance_view,
      a_start = a_start,
      a_end = a_end,
      control = reader_ctrl,
      read_outcomes = FALSE,
      include_design = include_design,
      enclose_differences = FALSE,
      expected_a_id = exp_slice,
      trace_env = ch_tr
    )
    if (!is.null(trace_env) && isTRUE(ch_tr$calls > 0L)) {
      trace_env$calls <- trace_env$calls + ch_tr$calls
      trace_env$pos_s_index <- c(trace_env$pos_s_index, ch_tr$pos_s_index)
      trace_env$pos_t_index <- c(trace_env$pos_t_index, ch_tr$pos_t_index)
    }
    ch
  }

  validate_chunk <- function(a_start, a_end) {
    ch <- query_slice(a_start, a_end, include_design = TRUE)
    Z_ch <- ch$design$Z_chunk
    if (anyNA(Z_ch) || any(!is.finite(Z_ch))) {
      cdid_logit_cert_abort(
        "`Z` contains NA, NaN, or infinite values.",
        c("cdid_logit_certificate_covariate_error", "cdid_logit_certificate_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    w_ch <- ch$weight
    if (anyNA(w_ch) || any(!is.finite(w_ch) | w_ch <= 0)) {
      cdid_logit_cert_abort(
        "`w` must contain finite, strictly positive weights (`w > 0`).",
        c("cdid_logit_certificate_weight_error", "cdid_logit_certificate_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    G_ch <- ch$g
    if (anyNA(G_ch) || any(is.nan(G_ch) | G_ch == -Inf)) {
      cdid_logit_cert_abort(
        "`G` must contain finite numeric cohort values or `+Inf` (never NA, NaN, or `-Inf`).",
        c("cdid_logit_certificate_cohort_error", "cdid_logit_certificate_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    po_ch <- ch$in_observed
    if (anyNA(po_ch)) {
      cdid_logit_cert_abort(
        "`pair_observed` cannot contain NA values.",
        c("cdid_logit_certificate_pair_error", "cdid_logit_certificate_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    list(
      j_ch = ch$cohort_index,
      po_ch = po_ch,
      outside_risk_pair_hit_count = as.integer(ch$counts$outside_risk_pair_hit_count),
      reader_ops = ch$work$actual_work_ops,
      n_endpoints_read = ch$read_diagnostics$n_endpoints_read
    )
  }

  read_mass_chunk <- function(a_start, a_end) {
    ch <- query_slice(a_start, a_end, include_design = FALSE)
    list(
      b_len = ch$chunk_rows,
      w_ch = as.double(ch$weight),
      j_ch = ch$cohort_index,
      po_ch = as.logical(ch$in_observed),
      reader_ops = ch$work$actual_work_ops,
      n_endpoints_read = ch$read_diagnostics$n_endpoints_read
    )
  }

  read_observer_chunk <- function(a_start, a_end) {
    ch <- query_slice(a_start, a_end, include_design = TRUE)
    b_len <- ch$chunk_rows
    list(
      idx = ch$a_index,
      b_len = b_len,
      Z_ch = matrix(as.double(ch$design$Z_chunk), nrow = b_len, ncol = p),
      w_ch = as.double(ch$weight),
      j_ch = ch$cohort_index,
      po_ch = as.logical(ch$in_observed),
      reader_ops = ch$work$actual_work_ops,
      n_endpoints_read = ch$read_diagnostics$n_endpoints_read
    )
  }

  logit_seal_access_context(
    list(
      mode = "compact_view",
      header = header,
      validate_chunk = validate_chunk,
      read_mass_chunk = read_mass_chunk,
      read_observer_chunk = read_observer_chunk
    ),
    mode = "compact_view"
  )
}

# ---------------------------------------------------------------------------
# Shared weighted-logit root certification worker
# ---------------------------------------------------------------------------

logit_cert_run_with_context <- function(ctx) {
  logit_validate_access_context(ctx)
  hdr <- ctx$header
  is_view <- identical(ctx$mode, "compact_view")
  n <- hdr$n
  p <- hdr$p
  m <- hdr$m
  z_colnames <- hdr$z_colnames
  coh_dbl <- hdr$coh_dbl
  cohort_labels <- hdr$cohort_labels
  beta0_dbl <- hdr$beta0_dbl
  F_dbl <- hdr$F_dbl
  r_dbl <- hdr$r_dbl
  ctrl <- hdr$ctrl
  preflight <- hdr$preflight
  step <- hdr$step
  n_iv <- gmm_iv(as.double(n))
  # Explicit typed automatic-radius request (compact view only). The radius is
  # read once; before the selection block it is referenced only by the
  # provenance/inconclusive builders. The selection record is retained in every
  # returned certificate (CERTIFIED or INCONCLUSIVE).
  auto_radius <- isTRUE(hdr$automatic_radius)
  radius_selection <- list(
    mode = if (auto_radius) .logit_cert_automatic_radius_mode else "caller_supplied",
    selected = FALSE,
    branch = NA_character_,
    r = r_dbl
  )

  # 3. Pass 1: Streamed validation of raw chunks (Z, w, G, pair_observed) and outward
  #    cohort/pair mass enclosures M_j, Pi_j without materializing full-N vectors.
  M_cohort <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
  Pi_cohort <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
  count_R <- integer(m)
  count_O <- integer(m)
  names(count_R) <- names(count_O) <- cohort_labels
  n_outside_risk <- 0L
  n_observed_outside_risk <- 0L
  w_total_iv <- c(0.0, 0.0)

  actual_ops <- 0
  actual_reader_ops <- 0
  total_endpoints_read <- 0L
  chunks_processed <- 0
  passes_completed <- 0L
  validation_chunks_processed <- 0
  validation_passes_completed <- 0L
  exp_evals <- 0
  max_sq_used <- 0L
  tri_solves <- 0

  make_work_list <- function() {
    base_work <- list(
      chunk_rows_limit = ctrl$max_chunk_rows,
      effective_chunk_rows = step,
      n_chunks_per_pass = preflight$n_chunks_per_pass,
      n_passes_completed = passes_completed,
      n_validation_passes_completed = validation_passes_completed,
      n_validation_chunks_processed = logit_cert_count_output(validation_chunks_processed),
      n_chunks_processed = logit_cert_count_output(chunks_processed),
      exp_endpoint_evaluations = logit_cert_count_output(exp_evals),
      max_squarings_used = max_sq_used,
      taylor_degree = ctrl$taylor_degree,
      triangular_vector_solves = logit_cert_count_output(tri_solves),
      actual_interval_ops = actual_ops,
      interval_ops_definition = "Conservative interval-operation charges at completed helper boundaries, including worst-case norm charges; not CPU instructions. A helper interrupted by refusal is not fully charged.",
      preflight_max_interval_ops = preflight$total_preflight_ops_bound,
      preflight_concurrent_entries = preflight$concurrent_entries,
      operation_order_cost = "One validation traversal plus2 arithmetic traversals: Pass 1 O(n + m) raw cohort/pair mass enclosures; Pass 2 O(|O| * (p^2 + N + m_sq)) raw c_i, eta_i, Lemma-exp, g, H, and interval triangular solves for a_i; Matrix stage O(p^3) interval triangular solves for F^{-T} H F^{-1} - I and F^{-T} g; O(B*p + p^2 + p + m) concurrent extra space."
    )
    if (is_view) {
      base_work$reader_chunk_rows_limit <- hdr$reader_ctrl$max_chunk_rows
      base_work$reader_chunks_processed <- logit_cert_count_output(validation_chunks_processed + chunks_processed)
      base_work$actual_reader_work_ops <- actual_reader_ops
      base_work$preflight_reader_work_ops <- preflight$total_reader_work_ops_bound
      base_work$preflight_combined_work_ops <- preflight$combined_preflight_ops_bound
      base_work$preflight_reader_concurrent_entries <- preflight$reader_concurrent_entries
      base_work$preflight_combined_concurrent_entries <- preflight$combined_concurrent_entries
      base_work$endpoints_read <- total_endpoints_read
    }
    base_work
  }

  make_provenance_list <- function() {
    base_prov <- list(
      normalizer_n = as.double(n),
      normalizer_definition = "count of all eligible A-units (including outside-risk and pair non-observers)",
      cohorts = coh_dbl,
      cohort_labels = cohort_labels,
      target_cohort = coh_dbl[1L],
      control_cohorts = coh_dbl[-1L],
      beta0 = beta0_dbl,
      F_factor = F_dbl,
      radius_r = r_dbl,
      control = ctrl,
      weights_source = "raw_w_enclosures_via_M_j_over_Pi_j",
      metric_definition = "exact_F_transpose_F",
      exp_method = "Lemma_exp_dyadic_range_reduction_Taylor_outward_binary64",
      outcomes_used = FALSE,
      decorative_metadata_used = FALSE
    )
    if (is_view) {
      base_prov$storage_mode <- "compact_view"
      base_prov$reader_control <- hdr$reader_ctrl
      base_prov$defining_pair <- hdr$defining_pair
      base_prov$read_diagnostics <- list(
        read_outcomes = FALSE,
        endpoints_read = total_endpoints_read,
        n_endpoints_read = total_endpoints_read,
        risk_filtered_before_outcome_read = TRUE,
        off_risk_endpoints_read = FALSE,
        nonobserver_endpoints_read = FALSE,
        cached_rounded_dy_returned = FALSE
      )
    }
    base_prov
  }

  # Helper to build inconclusive return when an arithmetic refusal occurs
  make_arithmetic_inconclusive <- function(err, partial_enclosures, dims_list) {
    reason <- if (!is.null(err$reason) && nzchar(err$reason)) err$reason else "interval_arithmetic_refusal"
    prov_list <- make_provenance_list()
    out_list <- list(
      certified = FALSE,
      status = "INCONCLUSIVE",
      inconclusive = TRUE,
      inconclusive_type = "arithmetic",
      inconclusive_reasons = reason,
      separation_certified = FALSE,
      sampling_certified = FALSE,
      public_gate_adopted = FALSE,
      raw_score_hessian_enclosures_complete = FALSE,
      partial_enclosure_scope = "Accumulated prefixes/stages only; these diagnostics need not enclose the full requested sums.",
      delta = Inf,
      t = Inf,
      a = Inf,
      a_control = Inf,
      r = r_dbl,
      x = Inf,
      q = Inf,
      selfmap = Inf,
      rho = Inf,
      radius_selection = radius_selection,
      inequalities = list(
        a_gt_0 = FALSE,
        delta_lt_1 = FALSE,
        x_lt_1 = FALSE,
        q_lt_1 = FALSE,
        selfmap_le_r = FALSE
      ),
      arithmetic_refusal = list(
        occurred = TRUE,
        reason = reason,
        message = conditionMessage(err),
        condition_class = class(err),
        diagnostics = err$diagnostics
      ),
      enclosures = partial_enclosures,
      dimensions = dims_list,
      read_diagnostics = prov_list$read_diagnostics,
      work = make_work_list(),
      provenance = prov_list
    )
    if (!is_view) out_list$read_diagnostics <- NULL
    structure(out_list, class = "cdid_logit_certificate")
  }

  # Pass 1a: Streamed validation of raw array chunks (before cohort/pair support checks)
  a_start <- 1L
  while (a_start <= n) {
    a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
    val_ch <- ctx$validate_chunk(a_start, a_end)
    actual_reader_ops <- actual_reader_ops + val_ch$reader_ops
    total_endpoints_read <- total_endpoints_read + val_ch$n_endpoints_read

    j_ch <- val_ch$j_ch
    po_ch <- val_ch$po_ch
    for (u in seq_along(j_ch)) {
      ju <- j_ch[u]
      if (is.na(ju)) {
        n_outside_risk <- n_outside_risk + 1L
        if (is.na(val_ch$outside_risk_pair_hit_count) && po_ch[u]) {
          n_observed_outside_risk <- n_observed_outside_risk + 1L
        }
      } else {
        count_R[ju] <- count_R[ju] + 1L
        if (po_ch[u]) count_O[ju] <- count_O[ju] + 1L
      }
    }
    if (!is.na(val_ch$outside_risk_pair_hit_count)) {
      n_observed_outside_risk <- n_observed_outside_risk + val_ch$outside_risk_pair_hit_count
    }
    validation_chunks_processed <- validation_chunks_processed + 1
    a_start <- as.double(a_end) + 1
  }
  validation_passes_completed <- 1L

  # Pass 1b: Streamed support counts and outward interval accumulation of M_j and Pi_j
  pass1_status <- tryCatch({
    a_start <- 1L
    while (a_start <= n) {
      a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
      m_ch <- ctx$read_mass_chunk(a_start, a_end)
      actual_reader_ops <- actual_reader_ops + m_ch$reader_ops
      total_endpoints_read <- total_endpoints_read + m_ch$n_endpoints_read
      b_len <- m_ch$b_len
      w_ch <- m_ch$w_ch
      j_ch <- m_ch$j_ch
      po_ch <- m_ch$po_ch

      for (u in seq_len(b_len)) {
        wu_iv <- gmm_iv(w_ch[u])
        w_total_iv <- gmm_iv_add(w_total_iv, wu_iv)
        actual_ops <- actual_ops + 1
        ju <- j_ch[u]
        if (!is.na(ju)) {
          M_cohort[ju, ] <- gmm_iv_add(M_cohort[ju, ], wu_iv)
          actual_ops <- actual_ops + 1
          if (po_ch[u]) {
            Pi_cohort[ju, ] <- gmm_iv_add(Pi_cohort[ju, ], wu_iv)
            actual_ops <- actual_ops + 1
          }
        }
      }
      chunks_processed <- chunks_processed + 1L
      a_start <- as.double(a_end) + 1
    }
    passes_completed <- 1L
    NULL
  }, cdid_logit_arithmetic_inconclusive = function(e) e,
     cdid_gmm_accuracy_error = function(e) e)

  # Check cohort and pair support BEFORE returning any arithmetic status
  if (any(count_R < 1L)) {
    cdid_logit_cert_abort(
      "Every declared cohort in `cohorts` must have positive eligible risk-set unit count and mass.",
      c("cdid_logit_certificate_cohort_support_error", "cdid_logit_certificate_support_error"),
      diagnostics = list(by_cohort_risk_count = count_R, cohorts = coh_dbl)
    )
  }
  if (any(count_O < 1L)) {
    cdid_logit_cert_abort(
      "Every declared cohort in `cohorts` must have positive pair-observed unit count and mass.",
      c("cdid_logit_certificate_pair_support_error", "cdid_logit_certificate_support_error"),
      diagnostics = list(by_cohort_observed_count = count_O, cohorts = coh_dbl)
    )
  }
  dims_list <- list(
    n_a_units = as.integer(n),
    p_covariates = as.integer(p),
    m_cohorts = as.integer(m),
    n_risk_units = as.integer(sum(count_R)),
    n_outside_risk_units = as.integer(n_outside_risk),
    n_observed_risk_units = as.integer(sum(count_O)),
    n_nonobserved_risk_units = as.integer(sum(count_R) - sum(count_O)),
    n_observed_outside_risk_units = as.integer(n_observed_outside_risk),
    n_target_observed = as.integer(count_O[1L]),
    n_control_observed = as.integer(sum(count_O[-1L])),
    by_cohort_risk_count = count_R,
    by_cohort_observed_count = count_O
  )

  if (!is.null(pass1_status)) {
    return(make_arithmetic_inconclusive(
      pass1_status,
      list(M_cohort = M_cohort, Pi_cohort = Pi_cohort, w_a_total = w_total_iv),
      dims_list
    ))
  }

  # 4. Pass 2 + Matrix/Contraction Stage inside narrow arithmetic refusal handler
  r_cohort <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
  q_cohort <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
  M_over_Pi_cohort <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
  c_range_cohort <- cbind(
    lo = rep(Inf, m),
    hi = rep(-Inf, m)
  )
  rownames(c_range_cohort) <- cohort_labels

  g_sum_iv <- matrix(0.0, nrow = p, ncol = 2L)
  H_sum_iv <- array(0.0, dim = c(p, p, 2L))
  c_sum_iv <- c(0.0, 0.0)
  eta_range_iv <- c(Inf, -Inf)
  odds_range_iv <- c(Inf, -Inf)
  prob_range_iv <- c(Inf, -Inf)
  var_range_iv <- c(Inf, -Inf)
  a_iv <- c(0.0, 0.0)
  a_control_iv <- c(0.0, 0.0)
  any_nonzero_z <- FALSE

  arith_err <- tryCatch({
    for (j in seq_len(m)) {
      if (M_cohort[j, 1L] <= 0.0 || Pi_cohort[j, 1L] <= 0.0) {
        logit_cert_arith_refuse(
          paste0("Nonpositive lower mass bound for cohort index ", j, "."),
          reason = "mass_lower_underflow",
          diagnostics = list(cohort_index = j, M_j = M_cohort[j, ], Pi_j = Pi_cohort[j, ])
        )
      }
      r_cohort[j, ] <- gmm_iv_div(Pi_cohort[j, ], M_cohort[j, ])
      q_cohort[j, ] <- gmm_iv_div(M_cohort[j, ], n_iv)
      mop_try <- tryCatch(
        gmm_iv_div(M_cohort[j, ], Pi_cohort[j, ]),
        cdid_gmm_accuracy_error = function(e) c(NA_real_, NA_real_),
        cdid_logit_arithmetic_inconclusive = function(e) c(NA_real_, NA_real_)
      )
      M_over_Pi_cohort[j, ] <- mop_try
      actual_ops <- actual_ops + 3
      if (r_cohort[j, 1L] <= 0.0 || q_cohort[j, 1L] <= 0.0 ||
          (!anyNA(mop_try) && mop_try[1L] <= 0.0)) {
        logit_cert_arith_refuse(
          paste0("Rate or cohort share lower endpoint underflowed to <= 0 for cohort ", j, "."),
          reason = "rate_lower_underflow",
          diagnostics = list(cohort_index = j, r_j = r_cohort[j, ], q_j = q_cohort[j, ])
        )
      }
    }

    a_start <- 1L
    while (a_start <= n) {
      a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
      obs_ch <- ctx$read_observer_chunk(a_start, a_end)
      actual_reader_ops <- actual_reader_ops + obs_ch$reader_ops
      total_endpoints_read <- total_endpoints_read + obs_ch$n_endpoints_read
      idx <- obs_ch$idx
      b_len <- obs_ch$b_len
      Z_ch <- obs_ch$Z_ch
      w_ch <- obs_ch$w_ch
      j_ch <- obs_ch$j_ch
      po_ch <- obs_ch$po_ch

      for (u in seq_len(b_len)) {
        ju <- j_ch[u]
        if (is.na(ju) || !po_ch[u]) next

        wu_iv <- gmm_iv(w_ch[u])
        z_u <- Z_ch[u, ]
        if (!any_nonzero_z && any(z_u != 0.0)) any_nonzero_z <- TRUE
        D_u <- if (ju == 1L) 1L else 0L

        # Enclose c_i = w_i * M_j / Pi_j from raw w_i and outward cohort mass enclosures
        c_iv1 <- if (!anyNA(M_over_Pi_cohort[ju, ])) {
          tryCatch(
            gmm_iv_mul(wu_iv, M_over_Pi_cohort[ju, ]),
            cdid_gmm_accuracy_error = function(e) NULL,
            cdid_logit_arithmetic_inconclusive = function(e) NULL
          )
        } else {
          NULL
        }
        if (!is.null(c_iv1)) actual_ops <- actual_ops + 1

        wM_try <- tryCatch(
          gmm_iv_mul(wu_iv, M_cohort[ju, ]),
          cdid_gmm_accuracy_error = function(e) NULL,
          cdid_logit_arithmetic_inconclusive = function(e) NULL
        )
        c_iv2 <- if (!is.null(wM_try)) {
          tryCatch(
            gmm_iv_div(wM_try, Pi_cohort[ju, ]),
            cdid_gmm_accuracy_error = function(e) NULL,
            cdid_logit_arithmetic_inconclusive = function(e) NULL
          )
        } else {
          NULL
        }
        if (!is.null(c_iv2)) actual_ops <- actual_ops + 2

        if (!is.null(c_iv1) && !is.null(c_iv2)) {
          c_iv <- logit_iv_intersect(c_iv1, c_iv2, "c_i = w_i * M_j / Pi_j")
        } else if (!is.null(c_iv1)) {
          c_iv <- c_iv1
        } else if (!is.null(c_iv2)) {
          c_iv <- c_iv2
        } else {
          logit_cert_arith_refuse(
            "Fitting weight c_i = w_i * M_j / Pi_j overflowed binary64.",
            reason = "fitting_weight_overflow",
            diagnostics = list(unit_offset = idx[u], cohort_index = ju)
          )
        }
        if (c_iv[1L] <= 0.0) {
          logit_cert_arith_refuse(
            "Fitting weight lower endpoint c_i underflowed to <= 0.",
            reason = "fitting_weight_underflow",
            diagnostics = list(unit_offset = idx[u], c_iv = c_iv)
          )
        }
        c_sum_iv <- gmm_iv_add(c_sum_iv, c_iv)
        actual_ops <- actual_ops + 1
        c_range_cohort[ju, 1L] <- min(c_range_cohort[ju, 1L], c_iv[1L])
        c_range_cohort[ju, 2L] <- max(c_range_cohort[ju, 2L], c_iv[2L])

        # Outward dot product eta_i = z_i' beta0
        eta_iv <- gmm_iv_dot(z_u, beta0_dbl)
        actual_ops <- actual_ops + 2 * p
        eta_range_iv[1L] <- min(eta_range_iv[1L], eta_iv[1L])
        eta_range_iv[2L] <- max(eta_range_iv[2L], eta_iv[2L])

        # Validated exponential o_i = exp(eta_i) via Lemma exp
        exp_res <- logit_iv_exp(
          eta_iv,
          taylor_degree = ctrl$taylor_degree,
          max_abs_eta = ctrl$max_abs_eta,
          max_squarings = ctrl$max_squarings
        )
        o_iv <- exp_res$iv
        exp_evals <- exp_evals + exp_res$eval_count
        max_sq_used <- max(max_sq_used, exp_res$max_m)
        actual_ops <- actual_ops + exp_res$ops

        # Logistic quantities p_i, 1 - p_i, v_i via Lemma logis
        logis_res <- logit_iv_logistic_from_odds(o_iv)
        actual_ops <- actual_ops + logis_res$ops
        p_iv <- logis_res$prob
        omp_iv <- logis_res$one_minus_prob
        v_iv <- logis_res$var

        odds_range_iv[1L] <- min(odds_range_iv[1L], o_iv[1L])
        odds_range_iv[2L] <- max(odds_range_iv[2L], o_iv[2L])
        prob_range_iv[1L] <- min(prob_range_iv[1L], p_iv[1L])
        prob_range_iv[2L] <- max(prob_range_iv[2L], p_iv[2L])
        var_range_iv[1L] <- min(var_range_iv[1L], v_iv[1L])
        var_range_iv[2L] <- max(var_range_iv[2L], v_iv[2L])

        # Score contribution c_i * (p_i - D_i) * z_i (using -(1 - p_i) when D_i == 1 to avoid cancellation)
        diff_iv <- if (D_u == 1L) gmm_iv_neg(omp_iv) else p_iv
        c_diff_iv <- gmm_iv_mul(c_iv, diff_iv)
        actual_ops <- actual_ops + 1
        for (k in seq_len(p)) {
          term_gk <- gmm_iv_mul(c_diff_iv, gmm_iv(z_u[k]))
          g_sum_iv[k, ] <- gmm_iv_add(g_sum_iv[k, ], term_gk)
          actual_ops <- actual_ops + 2
        }

        # Hessian contribution c_i * v_i * z_i * z_i' (upper triangle)
        c_v_iv <- gmm_iv_mul(c_iv, v_iv)
        actual_ops <- actual_ops + 1
        for (k in seq_len(p)) {
          cv_zk_iv <- gmm_iv_mul(c_v_iv, gmm_iv(z_u[k]))
          actual_ops <- actual_ops + 1
          for (l in seq.int(k, p)) {
            term_hkl <- gmm_iv_mul(cv_zk_iv, gmm_iv(z_u[l]))
            if (k == l) term_hkl[1L] <- max(0.0, term_hkl[1L])
            H_sum_iv[k, l, ] <- gmm_iv_add(H_sum_iv[k, l, ], term_hkl)
            actual_ops <- actual_ops + 2
          }
        }

        # Dual covariate norm ||F^{-T} z_i||_2 via interval forward substitution
        if (all(z_u == 0.0)) {
          ai_iv <- c(0.0, 0.0)
        } else {
          sv_z <- logit_iv_forwardsolve_ft(F_mat = F_dbl, b_iv = z_u)
          tri_solves <- tri_solves + 1L
          actual_ops <- actual_ops + sv_z$ops
          ai_iv <- gmm_iv_norm(sv_z$y[, 1L], sv_z$y[, 2L])
          actual_ops <- actual_ops + 3 * p + 2
        }
        a_iv[1L] <- max(a_iv[1L], ai_iv[1L])
        a_iv[2L] <- max(a_iv[2L], ai_iv[2L])
        if (D_u == 0L) {
          a_control_iv[1L] <- max(a_control_iv[1L], ai_iv[1L])
          a_control_iv[2L] <- max(a_control_iv[2L], ai_iv[2L])
        }
      }

      chunks_processed <- chunks_processed + 1L
      a_start <- as.double(a_end) + 1
    }
    passes_completed <- 2L
    NULL
  }, cdid_logit_arithmetic_inconclusive = function(e) e,
     cdid_gmm_accuracy_error = function(e) e)

  if (!is.null(arith_err)) {
    return(make_arithmetic_inconclusive(
      arith_err,
      list(
        M_cohort = M_cohort,
        Pi_cohort = Pi_cohort,
        r_cohort = r_cohort,
        q_cohort = q_cohort,
        M_over_Pi_cohort = M_over_Pi_cohort,
        c_fitting_weight_range_by_cohort = c_range_cohort,
        c_fitting_weight_sum = c_sum_iv,
        eta_range = eta_range_iv,
        odds_range = odds_range_iv,
        prob_range = prob_range_iv,
        var_range = var_range_iv
      ),
      dims_list
    ))
  }

  # 5. Normalize g(beta0) and H(beta0) by eligible-unit count n, then compute
  #    whitened quantities and Theorem root enclosures.
  g_iv <- matrix(0.0, nrow = p, ncol = 2L, dimnames = list(z_colnames, c("lo", "hi")))
  H_iv <- array(0.0, dim = c(p, p, 2L), dimnames = list(z_colnames, z_colnames, c("lo", "hi")))

  post_err <- tryCatch({
    for (k in seq_len(p)) {
      g_iv[k, ] <- gmm_iv_div(g_sum_iv[k, ], n_iv)
      actual_ops <- actual_ops + 1
      for (l in seq.int(k, p)) {
        h_kl <- gmm_iv_div(H_sum_iv[k, l, ], n_iv)
        if (k == l) h_kl[1L] <- max(0.0, h_kl[1L])
        H_iv[k, l, ] <- h_kl
        H_iv[l, k, ] <- h_kl
        actual_ops <- actual_ops + 1
      }
    }

    wh_H <- logit_iv_whitened_hessian_discrepancy(F_dbl, H_iv)
    actual_ops <- actual_ops + wh_H$ops
    tri_solves <- tri_solves + wh_H$solves

    sv_g <- logit_iv_forwardsolve_ft(F_dbl, g_iv)
    actual_ops <- actual_ops + sv_g$ops
    tri_solves <- tri_solves + 1L
    t_iv <- gmm_iv_norm(sv_g$y[, 1L], sv_g$y[, 2L])
    actual_ops <- actual_ops + 3 * p + 2

    list(wh_H = wh_H, whitened_g = sv_g$y, t_iv = t_iv)
  }, cdid_logit_arithmetic_inconclusive = function(e) e,
     cdid_gmm_accuracy_error = function(e) e)

  if (inherits(post_err, "error")) {
    return(make_arithmetic_inconclusive(
      post_err,
      list(
        M_cohort = M_cohort,
        Pi_cohort = Pi_cohort,
        r_cohort = r_cohort,
        q_cohort = q_cohort,
        M_over_Pi_cohort = M_over_Pi_cohort,
        c_fitting_weight_range_by_cohort = c_range_cohort,
        c_fitting_weight_sum = c_sum_iv,
        eta_range = eta_range_iv,
        odds_range = odds_range_iv,
        prob_range = prob_range_iv,
        var_range = var_range_iv,
        g = g_iv,
        H = H_iv
      ),
      dims_list
    ))
  }

  wh_H <- post_err$wh_H
  whitened_g_iv <- post_err$whitened_g
  delta_iv <- wh_H$delta_iv
  t_iv <- post_err$t_iv

  delta_upper <- delta_iv[2L]
  t_upper <- t_iv[2L]
  a_upper <- a_iv[2L]
  a_control_upper <- a_control_iv[2L]

  # Theorem root contraction quantities using upper endpoints (monotone outward enclosures):
  #   x = upper(a * r)
  #   q = upper((delta + x) / (1 - x))  [valid ONLY when x < 1]
  #   selfmap = upper(t + q * r)        [valid ONLY when x < 1]
  #   rho = upper(t / (1 - q))          [valid ONLY when x < 1 and q < 1]
  a_gt_0 <- isTRUE(any_nonzero_z && is.finite(a_upper) && a_upper > 0.0)
  delta_lt_1 <- isTRUE(is.finite(delta_upper) && delta_upper < 1.0)

  # Automatic radius (explicit request only). The radius-independent enclosures
  # above are complete; select ONE explicit outward radius here and let the
  # unchanged inequalities below decide. No grid, search, refit or retry. A
  # selection refusal is an arithmetic INCONCLUSIVE, never nonexistence.
  if (auto_radius) {
    sel <- tryCatch(
      logit_cert_select_automatic_radius(delta_iv, t_iv, a_iv, a_gt_0, delta_lt_1),
      cdid_logit_arithmetic_inconclusive = function(e) e,
      cdid_gmm_accuracy_error = function(e) e
    )
    actual_ops <- actual_ops + .logit_cert_automatic_radius_ops
    if (inherits(sel, "error")) {
      radius_selection$selected <- FALSE
      radius_selection$refusal_reason <- if (!is.null(sel$reason)) sel$reason else "interval_arithmetic_refusal"
      return(make_arithmetic_inconclusive(
        sel,
        list(
          M_cohort = M_cohort,
          Pi_cohort = Pi_cohort,
          r_cohort = r_cohort,
          q_cohort = q_cohort,
          M_over_Pi_cohort = M_over_Pi_cohort,
          c_fitting_weight_range_by_cohort = c_range_cohort,
          c_fitting_weight_sum = c_sum_iv,
          eta_range = eta_range_iv,
          odds_range = odds_range_iv,
          prob_range = prob_range_iv,
          var_range = var_range_iv,
          g = g_iv,
          H = H_iv,
          delta_iv = delta_iv,
          t_iv = t_iv,
          a_iv = a_iv,
          a_control_iv = a_control_iv
        ),
        dims_list
      ))
    }
    r_dbl <- sel$r
    radius_selection <- c(
      list(mode = .logit_cert_automatic_radius_mode, selected = TRUE),
      sel,
      list(metric = "exact_F_transpose_F from the caller-supplied upper-triangular F (current fitted H_chol in the pipeline)")
    )
  }

  x_iv <- c(NA_real_, Inf)
  x_upper <- Inf
  if (a_gt_0) {
    x_try <- tryCatch({
      actual_ops <- actual_ops + 2
      x_lo <- max(0.0, logit_b64_mul(gmm_iv(a_iv[1L]), gmm_iv(r_dbl))[1L])
      x_hi <- logit_b64_mul(gmm_iv(a_upper), gmm_iv(r_dbl))[2L]
      c(x_lo, x_hi)
    }, cdid_gmm_accuracy_error = function(e) NULL,
       cdid_logit_arithmetic_inconclusive = function(e) NULL)
    if (!is.null(x_try)) {
      x_iv <- x_try
      x_upper <- x_try[2L]
    }
  }
  x_lt_1 <- isTRUE(a_gt_0 && is.finite(x_upper) && x_upper < 1.0)

  q_iv <- c(NA_real_, Inf)
  q_upper <- Inf
  if (x_lt_1 && is.finite(delta_upper)) {
    q_try <- tryCatch({
      num_q <- logit_b64_add(c(delta_iv[1L], delta_upper), c(x_iv[1L], x_upper))
      den_q <- logit_b64_sub(gmm_iv(1.0), c(x_iv[1L], x_upper))
      actual_ops <- actual_ops + 3
      if (den_q[1L] <= 0.0) {
        logit_cert_arith_refuse("Denominator 1 - x nonpositive.", reason = "one_minus_x_nonpositive")
      }
      q_res <- logit_b64_div(num_q, den_q)
      c(max(0.0, q_res[1L]), q_res[2L])
    }, cdid_gmm_accuracy_error = function(e) NULL,
       cdid_logit_arithmetic_inconclusive = function(e) NULL)
    if (!is.null(q_try)) {
      q_iv <- q_try
      q_upper <- q_try[2L]
    }
  }
  q_lt_1 <- isTRUE(x_lt_1 && is.finite(q_upper) && q_upper < 1.0)

  selfmap_iv <- c(NA_real_, Inf)
  selfmap_upper <- Inf
  if (x_lt_1 && is.finite(q_upper) && is.finite(t_upper)) {
    sm_try <- tryCatch({
      qr_iv <- logit_b64_mul(c(q_iv[1L], q_upper), gmm_iv(r_dbl))
      actual_ops <- actual_ops + 2
      sm_res <- logit_b64_add(c(t_iv[1L], t_upper), qr_iv)
      c(max(0.0, sm_res[1L]), sm_res[2L])
    }, cdid_gmm_accuracy_error = function(e) NULL,
       cdid_logit_arithmetic_inconclusive = function(e) NULL)
    if (!is.null(sm_try)) {
      selfmap_iv <- sm_try
      selfmap_upper <- sm_try[2L]
    }
  }
  selfmap_le_r <- isTRUE(delta_lt_1 && x_lt_1 && q_lt_1 && is.finite(selfmap_upper) && selfmap_upper <= r_dbl)

  rho_iv <- c(NA_real_, Inf)
  rho_upper <- Inf
  if (q_lt_1 && is.finite(t_upper)) {
    rho_try <- tryCatch({
      den_rho <- logit_b64_sub(gmm_iv(1.0), c(q_iv[1L], q_upper))
      actual_ops <- actual_ops + 2
      if (den_rho[1L] <= 0.0) {
        logit_cert_arith_refuse("Denominator 1 - q nonpositive.", reason = "one_minus_q_nonpositive")
      }
      rho_res <- logit_b64_div(c(t_iv[1L], t_upper), den_rho)
      c(max(0.0, rho_res[1L]), rho_res[2L])
    }, cdid_gmm_accuracy_error = function(e) NULL,
       cdid_logit_arithmetic_inconclusive = function(e) NULL)
    if (!is.null(rho_try)) {
      rho_iv <- rho_try
      rho_upper <- rho_try[2L]
    }
  }

  reasons <- character(0L)
  if (!a_gt_0) reasons <- c(reasons, "a_not_positive")
  if (!delta_lt_1) reasons <- c(reasons, "delta_ge_1")
  if (!x_lt_1) reasons <- c(reasons, "x_ge_1")
  if (!q_lt_1) reasons <- c(reasons, "q_ge_1")
  if (!selfmap_le_r) reasons <- c(reasons, "selfmap_gt_r")
  if (!is.finite(rho_upper)) reasons <- c(reasons, "rho_nonfinite")

  certified <- (length(reasons) == 0L)

  final_cert_list <- list(
    certified = certified,
      status = if (certified) "CERTIFIED" else "INCONCLUSIVE",
      inconclusive = !certified,
      inconclusive_type = if (certified) "none" else "inequality",
      inconclusive_reasons = reasons,
      separation_certified = FALSE,
      sampling_certified = FALSE,
      public_gate_adopted = FALSE,
      raw_score_hessian_enclosures_complete = TRUE,
      delta = delta_upper,
      t = t_upper,
      a = a_upper,
      a_control = a_control_upper,
      r = r_dbl,
      x = x_upper,
      q = q_upper,
      selfmap = selfmap_upper,
      rho = rho_upper,
      radius_selection = radius_selection,
      inequalities = list(
        a_gt_0 = a_gt_0,
        delta_lt_1 = delta_lt_1,
        x_lt_1 = x_lt_1,
        q_lt_1 = q_lt_1,
        selfmap_le_r = selfmap_le_r
      ),
      arithmetic_refusal = list(
        occurred = FALSE,
        reason = "none",
        message = NULL,
        condition_class = NULL,
        diagnostics = NULL
      ),
      enclosures = list(
        M_cohort = M_cohort,
        Pi_cohort = Pi_cohort,
        r_cohort = r_cohort,
        q_cohort = q_cohort,
        M_over_Pi_cohort = M_over_Pi_cohort,
        c_fitting_weight_range_by_cohort = c_range_cohort,
        c_fitting_weight_sum = c_sum_iv,
        w_a_total = w_total_iv,
        eta_range = eta_range_iv,
        odds_range = odds_range_iv,
        prob_range = prob_range_iv,
        var_range = var_range_iv,
        g = g_iv,
        H = H_iv,
        H_lo = H_iv[, , 1L],
        H_hi = H_iv[, , 2L],
        M = wh_H$M_iv,
        M_lo = wh_H$M_iv[, , 1L],
        M_hi = wh_H$M_iv[, , 2L],
        E_hessian_minus_M = wh_H$E_iv,
        whitened_g = whitened_g_iv,
        whitened_H_minus_I_direct = wh_H$K_direct_iv,
        whitened_H_minus_I_residual = wh_H$K_residual_iv,
        whitened_H_minus_I = wh_H$K_iv,
        K_lo = wh_H$K_iv[, , 1L],
        K_hi = wh_H$K_iv[, , 2L],
        delta_iv = delta_iv,
        t_iv = t_iv,
        a_iv = a_iv,
        a_control_iv = a_control_iv,
        x_iv = x_iv,
        q_iv = q_iv,
        selfmap_iv = selfmap_iv,
        rho_iv = rho_iv
      ),
      dimensions = dims_list,
      read_diagnostics = if (is_view) make_provenance_list()$read_diagnostics else NULL,
      work = make_work_list(),
      provenance = make_provenance_list()
    )
  if (!is_view) final_cert_list$read_diagnostics <- NULL
  structure(final_cert_list, class = "cdid_logit_certificate")
}

# ---------------------------------------------------------------------------
# Main private raw weighted-logit root certificate entrypoints
# ---------------------------------------------------------------------------

certify_chained_logit_root <- function(
  Z,
  w,
  G,
  pair_observed,
  cohorts,
  beta0,
  F,
  r,
  control
) {
  logit_cert_require_gmm_core()
  if (missing(Z) || missing(w) || missing(G) || missing(pair_observed) ||
      missing(cohorts) || missing(beta0) || missing(F) || missing(r) ||
      missing(control)) {
    cdid_logit_cert_abort(
      "All arguments (`Z`, `w`, `G`, `pair_observed`, `cohorts`, `beta0`, `F`, `r`, `control`) must be supplied explicitly.",
      "cdid_logit_certificate_input_error"
    )
  }

  # 1. Validate control and lightweight metadata BEFORE any full-array coercion or allocation
  ctrl <- validate_logit_certificate_control(control)
  # Header-only cohort checks precede any content masks/coercion/labels.
  if (is.null(cohorts) || !is.atomic(cohorts) || !is.null(dim(cohorts)) ||
      is.object(cohorts) || !is.numeric(cohorts) || length(cohorts) < 2)
    cdid_logit_cert_abort("An explicit ordered numeric cohort vector is required.",
      c("cdid_logit_certificate_cohort_error", "cdid_logit_certificate_input_error"))
  m <- length(cohorts)

  if (is.null(Z) || !is.matrix(Z) || is.object(Z) || !is.numeric(Z)) {
    cdid_logit_cert_abort(
      "`Z` must be a plain numeric matrix (one row per eligible A-unit).",
      c("cdid_logit_certificate_covariate_error", "cdid_logit_certificate_input_error")
    )
  }
  n <- nrow(Z)
  p <- ncol(Z)
  if (is.na(n) || is.na(p) || n < 1L || p < 1L) {
    cdid_logit_cert_abort(
      "`Z` must have `nrow(Z) >= 1` and `ncol(Z) >= 1`.",
      c("cdid_logit_certificate_covariate_error", "cdid_logit_certificate_input_error")
    )
  }

  if (is.null(w) || !is.atomic(w) || !is.null(dim(w)) || is.object(w) ||
      !is.numeric(w) || length(w) != n) {
    cdid_logit_cert_abort(
      "`w` must be an atomic numeric vector of length `nrow(Z)`.",
      c("cdid_logit_certificate_weight_error", "cdid_logit_certificate_input_error")
    )
  }
  if (is.null(G) || !is.atomic(G) || !is.null(dim(G)) || is.object(G) ||
      !is.numeric(G) || length(G) != n) {
    cdid_logit_cert_abort(
      "`G` must be an atomic numeric vector of length `nrow(Z)`.",
      c("cdid_logit_certificate_cohort_error", "cdid_logit_certificate_input_error")
    )
  }
  if (is.null(pair_observed) || !is.atomic(pair_observed) || !is.null(dim(pair_observed)) ||
      is.object(pair_observed) || !is.logical(pair_observed) || length(pair_observed) != n) {
    cdid_logit_cert_abort(
      "`pair_observed` must be an atomic logical vector of length `nrow(Z)`.",
      c("cdid_logit_certificate_pair_error", "cdid_logit_certificate_input_error")
    )
  }
  if (is.null(beta0) || !is.atomic(beta0) || !is.null(dim(beta0)) || is.object(beta0) ||
      !is.numeric(beta0) || length(beta0) != p) {
    cdid_logit_cert_abort(
      "`beta0` must be a finite numeric vector of length `ncol(Z)`.",
      c("cdid_logit_certificate_beta_error", "cdid_logit_certificate_input_error")
    )
  }
  if (!is.null(attr(beta0, "converged")) || !is.null(attr(beta0, "certified"))) {
    cdid_logit_cert_abort(
      "Decorative convergence/certification attributes on `beta0` are refused.",
      c("cdid_logit_certificate_beta_error", "cdid_logit_certificate_input_error")
    )
  }
  if (!is.matrix(F) || is.object(F) || !is.numeric(F) || nrow(F) != p || ncol(F) != p)
    cdid_logit_cert_abort("`F` must be an explicit numeric p x p matrix.",
      c("cdid_logit_certificate_metric_error", "cdid_logit_certificate_input_error"))
  if (!gmm_scalar_numeric(r, positive = TRUE)) {
    cdid_logit_cert_abort(
      "`r` must be an explicit finite positive numeric scalar radius.",
      c("cdid_logit_certificate_radius_error", "cdid_logit_certificate_input_error")
    )
  }
  r_dbl <- as.double(r)

  # 2. Check preflight chunk, matrix-entry, and work-operation bounds BEFORE materialization
  preflight <- preflight_logit_certificate(n = n, p = p, m = m, control = ctrl)
  if (m > n) cdid_logit_cert_abort("Declared cohorts cannot all have positive mass when their count exceeds n_A.",
    c("cdid_logit_certificate_cohort_support_error", "cdid_logit_certificate_support_error"))
  coh_dbl <- validate_logit_certificate_cohorts(cohorts)
  cohort_labels <- c("g", paste0("k_", seq_len(m - 1L)))
  if (anyNA(beta0) || any(!is.finite(beta0)))
    cdid_logit_cert_abort("`beta0` must be finite.", c("cdid_logit_certificate_beta_error", "cdid_logit_certificate_input_error"))
  beta0_dbl <- as.double(beta0)
  F_dbl <- validate_logit_certificate_metric_factor(F, p)

  ctx <- construct_dense_logit_access_context(
    Z = Z,
    w = w,
    G = G,
    pair_observed = pair_observed,
    coh_dbl = coh_dbl,
    cohort_labels = cohort_labels,
    beta0_dbl = beta0_dbl,
    F_dbl = F_dbl,
    r_dbl = r_dbl,
    ctrl = ctrl,
    preflight = preflight
  )
  logit_cert_run_with_context(ctx)
}

certify_chained_logit_root_view <- function(
  nuisance_view,
  beta0,
  F,
  r,
  control,
  reader_control
) {
  logit_cert_require_gmm_core()
  if (missing(nuisance_view) || missing(beta0) || missing(F) ||
      missing(r) || missing(control) || missing(reader_control)) {
    cdid_logit_cert_abort(
      "All arguments (`nuisance_view`, `beta0`, `F`, `r`, `control`, `reader_control`) must be supplied explicitly.",
      "cdid_logit_certificate_input_error"
    )
  }

  ctrl <- validate_logit_certificate_control(control)
  expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  trace_env <- attr(reader_control, "trace_env", exact = TRUE)
  r_ctrl <- validate_ipw_owner_control(reader_control)
  validate_ipw_nuisance_view(nuisance_view)

  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  n <- as.integer(frame_owner$n_a_units)
  p <- as.integer(frame_owner$n_covariates)
  cohorts <- nuisance_view$cohorts$cohorts
  m <- length(cohorts)

  if (is.null(beta0) || !is.atomic(beta0) || !is.null(dim(beta0)) || is.object(beta0) ||
      !is.numeric(beta0) || length(beta0) != p) {
    cdid_logit_cert_abort(
      "`beta0` must be a finite numeric vector of length `ncol(Z)`.",
      c("cdid_logit_certificate_beta_error", "cdid_logit_certificate_input_error")
    )
  }
  if (!is.null(attr(beta0, "converged")) || !is.null(attr(beta0, "certified"))) {
    cdid_logit_cert_abort(
      "Decorative convergence/certification attributes on `beta0` are refused.",
      c("cdid_logit_certificate_beta_error", "cdid_logit_certificate_input_error")
    )
  }
  if (!is.matrix(F) || is.object(F) || !is.numeric(F) || nrow(F) != p || ncol(F) != p) {
    cdid_logit_cert_abort(
      "`F` must be an explicit numeric p x p matrix.",
      c("cdid_logit_certificate_metric_error", "cdid_logit_certificate_input_error")
    )
  }
  # Explicit typed automatic-radius request, or an explicit positive radius.
  # Any invalid radius (including NA_real_, NaN, decorated missing/nonfinite
  # or multi-valued numerics) keeps the original typed refusal. Decorated
  # positive numerics keep the existing manual path at this view entry.
  r_auto <- logit_cert_is_automatic_radius_request(r)
  if (!r_auto && !gmm_scalar_numeric(r, positive = TRUE)) {
    cdid_logit_cert_abort(
      "`r` must be an explicit finite positive numeric scalar radius.",
      c("cdid_logit_certificate_radius_error", "cdid_logit_certificate_input_error")
    )
  }
  r_dbl <- if (r_auto) NA_real_ else as.double(r)

  preflight <- preflight_logit_certificate_view(
    nuisance_view = nuisance_view,
    control = ctrl,
    reader_control = r_ctrl,
    expected_a_id = expected_a_id,
    automatic_radius = r_auto
  )
  if (m > n) {
    cdid_logit_cert_abort(
      "Declared cohorts cannot all have positive mass when their count exceeds n_A.",
      c("cdid_logit_certificate_cohort_support_error", "cdid_logit_certificate_support_error")
    )
  }
  coh_dbl <- validate_logit_certificate_cohorts(cohorts)
  cohort_labels <- c("g", paste0("k_", seq_len(m - 1L)))
  if (anyNA(beta0) || any(!is.finite(beta0))) {
    cdid_logit_cert_abort(
      "`beta0` must be finite.",
      c("cdid_logit_certificate_beta_error", "cdid_logit_certificate_input_error")
    )
  }
  beta0_dbl <- as.double(beta0)
  beta0_dbl[beta0_dbl == 0.0] <- 0.0
  F_dbl <- validate_logit_certificate_metric_factor(F, p)
  F_dbl[F_dbl == 0.0] <- 0.0
  if (!r_auto) r_dbl[r_dbl == 0.0] <- 0.0

  if (!is.null(nuisance_view$candidates)) {
    vc <- nuisance_view$candidates
    if (!is.null(vc$beta0)) {
      cand_b <- as.double(vc$beta0)
      cand_b[cand_b == 0.0] <- 0.0
      if (!identical(beta0_dbl, cand_b)) {
        cdid_logit_cert_abort(
          "Explicit `beta0` disagrees with optional candidate stored on `nuisance_view`.",
          c("cdid_logit_certificate_candidate_mismatch_error", "cdid_logit_certificate_beta_error", "cdid_logit_certificate_input_error")
        )
      }
    }
    if (!is.null(vc$F)) {
      cand_F <- matrix(as.double(vc$F), nrow = p, ncol = p)
      cand_F[cand_F == 0.0] <- 0.0
      if (!identical(F_dbl, cand_F)) {
        cdid_logit_cert_abort(
          "Explicit `F` disagrees with optional candidate stored on `nuisance_view`.",
          c("cdid_logit_certificate_candidate_mismatch_error", "cdid_logit_certificate_metric_error", "cdid_logit_certificate_input_error")
        )
      }
    }
    if (!is.null(vc$r)) {
      if (r_auto) {
        cdid_logit_cert_abort(
          "An automatic-radius request conflicts with an explicit candidate radius stored on `nuisance_view`.",
          c("cdid_logit_certificate_candidate_mismatch_error", "cdid_logit_certificate_radius_error", "cdid_logit_certificate_input_error")
        )
      }
      cand_r <- as.double(vc$r)
      cand_r[cand_r == 0.0] <- 0.0
      if (!identical(r_dbl, cand_r)) {
        cdid_logit_cert_abort(
          "Explicit `r` disagrees with optional candidate stored on `nuisance_view`.",
          c("cdid_logit_certificate_candidate_mismatch_error", "cdid_logit_certificate_radius_error", "cdid_logit_certificate_input_error")
        )
      }
    }
  }

  ctx <- construct_view_logit_access_context(
    nuisance_view = nuisance_view,
    coh_dbl = coh_dbl,
    cohort_labels = cohort_labels,
    beta0_dbl = beta0_dbl,
    F_dbl = F_dbl,
    r_dbl = r_dbl,
    ctrl = ctrl,
    reader_ctrl = r_ctrl,
    preflight = preflight,
    expected_a_id = expected_a_id,
    trace_env = trace_env,
    automatic_radius = r_auto
  )
  logit_cert_run_with_context(ctx)
}

chained_logit_certificate <- certify_chained_logit_root
certify_weighted_logit_root <- certify_chained_logit_root
chained_logit_certificate_view <- certify_chained_logit_root_view
certify_weighted_logit_root_view <- certify_chained_logit_root_view

