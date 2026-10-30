# Private CP02 raw-endpoint IPW link functional and complete IF enclosure layer.
#
# Implements the conditional numerical propagation primitives of:
#   * docs/theory/cp02_ipw_gmm_numerical_propagation.tex
#     (Section 1 Definition if & Lemma struct; Section 2 Lemmas incl, exp, logis, dual;
#      Section 3 Lemma mult; Section 4 Theorem A; Section 5 Theorem B & Corollary sep)
#   * docs/theory/cp02_finite_logit_validation.tex (Theorem root, Corollary point)
#   * docs/theory/cp02_weighted_pair_if.tex (complete 4-term case-weight derivative)
#
# Scope and boundaries:
#   * Reusable, explicitly bounded certified nuisance-problem constructor from raw
#     all-A unit IDs/order, Z, w, G, pair_observed (S), ordered cohorts (target first),
#     beta0, upper-triangular F, radius r, raw logit certificate control, and explicit
#     problem retention/work controls. Calls certify_chained_logit_root() on the SAME
#     retained raw inputs; never trusts foreign certificates or decorative labels.
#   * Seals nuisance problem and link enclosure objects with locked environments and
#     structural/copy-on-modify identity verification (ipw_seal pattern); no global cache.
#   * For a single explicit link, takes raw stored left/right endpoint vectors (y_left,
#     y_right in all-A order) and stored loading/action candidates plus explicit resource
#     controls. Differences Delta_i = y_{it} - y_{is} are enclosed outward directly from
#     stored endpoints, never from cached rounded dy. Never reads unavailable endpoints
#     or off-risk/non-observer outcome indices.
#   * Encloses T, C, mu, d^+, d^-, v (-G_beta), gamma (G_r), L (A_beta_r for ALL cohorts
#     including target), lambda_beta = -H^{-1} v, and lambda_r = gamma + L' lambda_beta
#     at beta* (Theorem A + Lemma dual(iii)) and at beta0 (Route B), plus separated
#     root-perturbation bounds (Theorem B) and evaluation-vs-root error bounds (Cor. sep).
#   * Emits complete 4-term IF row enclosures (rate, score, treated, control) and shared
#     xi_rate / xi_beta actions by bounded row chunks, storing outcome sufficient
#     statistics/loadings and two compact observer endpoint copies per link.
#     Retaining many links costs O(sum_l |O_l| + n*links) including position maps;
#     this is not a shared-panel/period-endpoint storage implementation.
#   * Does NOT implement covariance accumulation, changing-weight GMM propagation,
#     public adapters, catalog selection (D10), cluster weights (D11), or acceptance
#     policy (D12). Scalar/interval domain or positivity failures return explicit
#     INCONCLUSIVE, never separation or non-identification claims.

cdid_ipw_enclosure_abort <- function(message, class = character(), diagnostics = list()) {
  stop(structure(
    class = unique(c(
      class,
      "cdid_ipw_enclosure_refusal",
      "cdid_ipw_enclosure_error",
      "cdid_ipw_certificate_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

ipw_enclosure_arith_refuse <- function(message, reason, diagnostics = list()) {
  stop(structure(
    class = unique(c(
      "cdid_ipw_enclosure_arithmetic_inconclusive",
      "cdid_logit_arithmetic_inconclusive",
      "cdid_ipw_enclosure_inconclusive",
      "cdid_ipw_enclosure_refusal",
      "cdid_ipw_enclosure_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, reason = reason, call = NULL, diagnostics = diagnostics)
  ))
}

ipw_enclosure_require_deps <- function() {
  ns_env <- environment(ipw_enclosure_require_deps)
  req_fns <- c(
    "gmm_next", "gmm_iv", "gmm_iv_check", "gmm_iv_add", "gmm_iv_neg",
    "gmm_iv_sub", "gmm_iv_mul", "gmm_iv_div", "gmm_iv_sqrt", "gmm_iv_abs",
    "gmm_iv_dot", "gmm_iv_norm", "gmm_scalar_logical", "gmm_scalar_numeric",
    "gmm_scalar_character", "ldexp_safe",
    "logit_cert_int_bound", "logit_cert_safe_product", "logit_cert_safe_sum",
    "logit_cert_count_output", "logit_b64_add", "logit_b64_sub", "logit_b64_mul",
    "logit_b64_div", "logit_iv_intersect", "logit_iv_exp",
    "logit_iv_logistic_from_odds", "logit_iv_forwardsolve_ft",
    "logit_iv_backsolve_f", "validate_logit_certificate_control",
    "validate_logit_certificate_cohorts", "validate_logit_certificate_metric_factor",
    "preflight_logit_certificate", "certify_chained_logit_root"
  )
  for (fn in req_fns) {
    if (!exists(fn, envir = ns_env, mode = "function", inherits = TRUE)) {
      cdid_ipw_enclosure_abort(
        paste0("Required private helper `", fn, "` is not loaded in the namespace."),
        "cdid_ipw_enclosure_dependency_error"
      )
    }
  }
  invisible(TRUE)
}

ipw_enclosure_check_interrupt <- function() {
  Sys.sleep(0)
}

# Typed arithmetic failures provide no certificate and no partial numerical output.
ipw_enclosure_capture_arithmetic <- function(expr) {
  tryCatch(expr,
    cdid_ipw_enclosure_arithmetic_inconclusive = function(e) e,
    cdid_logit_arithmetic_inconclusive = function(e) e,
    cdid_gmm_accuracy_error = function(e) e)
}

ipw_enclosure_inconclusive_result <- function(err, result_class, fields, binding = NULL) {
  reason <- if (!is.null(err$reason) && nzchar(err$reason)) err$reason else "arithmetic_refusal"
  obj <- structure(c(list(
    certified = FALSE, status = "INCONCLUSIVE", inconclusive = TRUE,
    inconclusive_type = "arithmetic", inconclusive_reasons = reason,
    separation_certified = FALSE, sampling_certified = FALSE, public_gate_adopted = FALSE,
    arithmetic_refusal = list(occurred = TRUE, reason = reason,
      message = conditionMessage(err), condition_class = class(err), diagnostics = err$diagnostics)
  ), fields), class = result_class)
  if (is.null(binding)) obj else ipw_enclosure_seal(obj, binding)
}

# Upper charges for the frozen interval helpers, including BOTH exp endpoints.
# Forward substitution <= p^2; explicit F^-T product <= 2p^2; norm <= 3p+2.
# Residual solve adds 2p^2+3p+1; ball dot adds 2p+3.
ipw_enclosure_helper_work <- function(p, certificate_control) {
  cc <- validate_logit_certificate_control(certificate_control)
  p2 <- ipw_enclosure_safe_product(p, p, label = "helper p^2")
  dual <- ipw_enclosure_safe_sum(c(3*p2, 3*p, 2), label = "dual helper work")
  list(dual = dual,
    solve = ipw_enclosure_safe_sum(c(2*p2, dual, 3*p, 1), label = "inverse-action work"),
    ball = ipw_enclosure_safe_sum(c(dual, 2*p, 3), label = "ball-dot work"),
    exp = ipw_enclosure_safe_sum(c(10*cc$taylor_degree, 2*cc$max_squarings, 8), label = "two-endpoint exp work"))
}

ipw_enclosure_int_bound <- function(x, name, min_val = 1) {
  if (!gmm_scalar_numeric(x) || x != floor(x) || x < min_val || x > 2^53 - 1) {
    cdid_ipw_enclosure_abort(
      paste0("`", name, "` must be an explicit finite integer in [", min_val, ", 2^53 - 1]."),
      c("cdid_ipw_enclosure_bound_error", "cdid_ipw_enclosure_control_error")
    )
  }
  as.double(x)
}

ipw_enclosure_safe_product <- function(a, b, max_entries = NULL, label = "allocation") {
  if (!gmm_scalar_numeric(a, nonneg = TRUE) || a != floor(a) ||
      !gmm_scalar_numeric(b, nonneg = TRUE) || b != floor(b)) {
    cdid_ipw_enclosure_abort(
      paste0("Invalid integer dimensions for ", label, "."),
      "cdid_ipw_enclosure_bound_error"
    )
  }
  da <- as.double(a)
  db <- as.double(b)
  limit <- 2^53 - 1
  if (da > 0 && db > limit / da) {
    cdid_ipw_enclosure_abort(
      paste0("Count product for ", label, " overflows exact IEEE binary64 integer range."),
      "cdid_ipw_enclosure_bound_error",
      diagnostics = list(label = label, dim_a = da, dim_b = db)
    )
  }
  prod_val <- da * db
  if (!is.finite(prod_val) || prod_val > limit) {
    cdid_ipw_enclosure_abort(
      paste0("Count product for ", label, " exceeds exact integer range."),
      "cdid_ipw_enclosure_bound_error"
    )
  }
  if (!is.null(max_entries) && prod_val > max_entries) {
    cdid_ipw_enclosure_abort(
      paste0("Count for ", label, " (", prod_val, ") exceeds explicit bound (", max_entries, ")."),
      "cdid_ipw_enclosure_bound_error",
      diagnostics = list(label = label, requested = prod_val, bound = max_entries)
    )
  }
  prod_val
}

ipw_enclosure_safe_sum <- function(vals, max_bound = NULL, label = "sum") {
  limit <- 2^53 - 1
  total <- 0
  for (v in vals) {
    dv <- as.double(v)
    if (!is.finite(dv) || dv < 0 || dv != floor(dv) || dv > limit - total) {
      cdid_ipw_enclosure_abort(
        paste0("Cumulative count for ", label, " overflows exact IEEE binary64 integer range."),
        "cdid_ipw_enclosure_bound_error",
        diagnostics = list(label = label, total_so_far = total, term = dv)
      )
    }
    total <- total + dv
  }
  if (!is.null(max_bound) && total > max_bound) {
    cdid_ipw_enclosure_abort(
      paste0("Total ", label, " (", total, ") exceeds explicit bound (", max_bound, ")."),
      "cdid_ipw_enclosure_bound_error",
      diagnostics = list(label = label, requested = total, bound = max_bound)
    )
  }
  total
}

ipw_enclosure_seal <- function(obj, binding) {
  ref <- new.env(parent = emptyenv())
  obj$ref <- ref
  assign(binding, obj, envir = ref)
  lockEnvironment(ref, bindings = TRUE)
  obj
}

ipw_enclosure_seal_ok <- function(obj, binding) {
  ref <- if (is.list(obj)) obj$ref else NULL
  is.environment(ref) && identical(parent.env(ref), emptyenv()) &&
    environmentIsLocked(ref) &&
    exists(binding, envir = ref, inherits = FALSE) &&
    bindingIsLocked(binding, ref) &&
    identical(obj, get(binding, envir = ref, inherits = FALSE))
}

# ---------------------------------------------------------------------------
# Explicit resource controls (NO internal defaults)
# ---------------------------------------------------------------------------

ipw_nuisance_problem_control <- function(
  max_chunk_rows,
  max_retained_entries,
  max_matrix_entries,
  max_work_ops
) {
  ipw_enclosure_require_deps()
  if (missing(max_chunk_rows) || missing(max_retained_entries) ||
      missing(max_matrix_entries) || missing(max_work_ops)) {
    cdid_ipw_enclosure_abort(
      "All fields of `ipw_nuisance_problem_control` (`max_chunk_rows`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`) must be supplied explicitly by the caller (no internal defaults).",
      c("cdid_ipw_enclosure_control_error", "cdid_ipw_enclosure_bound_error")
    )
  }
  max_ch <- ipw_enclosure_int_bound(max_chunk_rows, "max_chunk_rows", min_val = 1)
  max_ret <- ipw_enclosure_int_bound(max_retained_entries, "max_retained_entries", min_val = 1)
  max_me <- ipw_enclosure_int_bound(max_matrix_entries, "max_matrix_entries", min_val = 1)
  max_wo <- ipw_enclosure_int_bound(max_work_ops, "max_work_ops", min_val = 1)

  structure(
    list(
      max_chunk_rows = max_ch,
      max_retained_entries = max_ret,
      max_matrix_entries = max_me,
      max_work_ops = max_wo
    ),
    class = "cdid_ipw_nuisance_problem_control"
  )
}

validate_ipw_nuisance_problem_control <- function(control) {
  ipw_enclosure_require_deps()
  if (missing(control) || is.null(control) || !is.list(control) || length(control) != 4L ||
      is.null(names(control)) || anyDuplicated(names(control))) {
    cdid_ipw_enclosure_abort(
      "Explicit `problem_control` list (`max_chunk_rows`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`) is required; no internal defaults exist.",
      c("cdid_ipw_enclosure_control_error", "cdid_ipw_enclosure_bound_error")
    )
  }
  req <- c("max_chunk_rows", "max_retained_entries", "max_matrix_entries", "max_work_ops")
  if (!setequal(names(control), req)) {
    missing_f <- setdiff(req, names(control))
    extra_f <- setdiff(names(control), req)
    cdid_ipw_enclosure_abort(
      paste0(
        "Invalid `problem_control` fields.",
        if (length(missing_f)) paste0(" Missing: ", paste(missing_f, collapse = ", "), ".") else "",
        if (length(extra_f)) paste0(" Unexpected: ", paste(extra_f, collapse = ", "), ".") else ""
      ),
      c("cdid_ipw_enclosure_control_error", "cdid_ipw_enclosure_bound_error")
    )
  }
  ipw_nuisance_problem_control(
    max_chunk_rows = control$max_chunk_rows,
    max_retained_entries = control$max_retained_entries,
    max_matrix_entries = control$max_matrix_entries,
    max_work_ops = control$max_work_ops
  )
}

ipw_link_enclosure_control <- function(
  max_chunk_rows,
  max_retained_entries,
  max_matrix_entries,
  max_work_ops,
  max_output_entries
) {
  ipw_enclosure_require_deps()
  if (missing(max_chunk_rows) || missing(max_retained_entries) ||
      missing(max_matrix_entries) || missing(max_work_ops) ||
      missing(max_output_entries)) {
    cdid_ipw_enclosure_abort(
      "All fields of `ipw_link_enclosure_control` (`max_chunk_rows`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`) must be supplied explicitly by the caller (no internal defaults).",
      c("cdid_ipw_enclosure_control_error", "cdid_ipw_enclosure_bound_error")
    )
  }
  max_ch <- ipw_enclosure_int_bound(max_chunk_rows, "max_chunk_rows", min_val = 1)
  max_ret <- ipw_enclosure_int_bound(max_retained_entries, "max_retained_entries", min_val = 1)
  max_me <- ipw_enclosure_int_bound(max_matrix_entries, "max_matrix_entries", min_val = 1)
  max_wo <- ipw_enclosure_int_bound(max_work_ops, "max_work_ops", min_val = 1)
  max_out <- ipw_enclosure_int_bound(max_output_entries, "max_output_entries", min_val = 1)

  structure(
    list(
      max_chunk_rows = max_ch,
      max_retained_entries = max_ret,
      max_matrix_entries = max_me,
      max_work_ops = max_wo,
      max_output_entries = max_out
    ),
    class = "cdid_ipw_link_enclosure_control"
  )
}

validate_ipw_link_enclosure_control <- function(control) {
  ipw_enclosure_require_deps()
  if (missing(control) || is.null(control) || !is.list(control) || length(control) != 5L ||
      is.null(names(control)) || anyDuplicated(names(control))) {
    cdid_ipw_enclosure_abort(
      "Explicit `control` list (`max_chunk_rows`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`) is required; no internal defaults exist.",
      c("cdid_ipw_enclosure_control_error", "cdid_ipw_enclosure_bound_error")
    )
  }
  req <- c("max_chunk_rows", "max_retained_entries", "max_matrix_entries", "max_work_ops", "max_output_entries")
  if (!setequal(names(control), req)) {
    missing_f <- setdiff(req, names(control))
    extra_f <- setdiff(names(control), req)
    cdid_ipw_enclosure_abort(
      paste0(
        "Invalid link enclosure `control` fields.",
        if (length(missing_f)) paste0(" Missing: ", paste(missing_f, collapse = ", "), ".") else "",
        if (length(extra_f)) paste0(" Unexpected: ", paste(extra_f, collapse = ", "), ".") else ""
      ),
      c("cdid_ipw_enclosure_control_error", "cdid_ipw_enclosure_bound_error")
    )
  }
  ipw_link_enclosure_control(
    max_chunk_rows = control$max_chunk_rows,
    max_retained_entries = control$max_retained_entries,
    max_matrix_entries = control$max_matrix_entries,
    max_work_ops = control$max_work_ops,
    max_output_entries = control$max_output_entries
  )
}

# ---------------------------------------------------------------------------
# Preflight accounting (executed BEFORE content scans or copy allocations)
# ---------------------------------------------------------------------------

preflight_ipw_nuisance_problem <- function(n, p, m, certificate_control, problem_control) {
  ipw_enclosure_require_deps()
  cert_ctrl <- validate_logit_certificate_control(certificate_control)
  prob_ctrl <- validate_ipw_nuisance_problem_control(problem_control)

  n_d <- ipw_enclosure_int_bound(n, "nrow(Z)", min_val = 1)
  p_d <- ipw_enclosure_int_bound(p, "ncol(Z)", min_val = 1)
  m_d <- ipw_enclosure_int_bound(m, "length(cohorts)", min_val = 2)

  # Raw certificate preflight covers only certify_chained_logit_root's own work
  raw_pre <- preflight_logit_certificate(n = n_d, p = p_d, m = m_d, control = cert_ctrl)

  max_ret <- prob_ctrl$max_retained_entries
  max_me <- prob_ctrl$max_matrix_entries
  max_wo <- prob_ctrl$max_work_ops
  B_eff <- min(n_d, prob_ctrl$max_chunk_rows)
  N_deg <- cert_ctrl$taylor_degree
  S_max <- cert_ctrl$max_squarings

  np <- ipw_enclosure_safe_product(n_d, p_d, max_ret, "retained Z (n x p)")
  p2 <- ipw_enclosure_safe_product(p_d, p_d, max_ret, "p x p matrix block")
  pm <- ipw_enclosure_safe_product(p_d, m_d, max_ret, "p x cohorts Jacobian block")

  # Retained raw inputs + derived risk/observer masks + retained nuisance enclosures
  retained_raw_entries <- ipw_enclosure_safe_sum(
    c(
      np,
      ipw_enclosure_safe_product(5, n_d, max_ret, "retained id, id_key, w, G, pair_observed"),
      m_d,
      p_d,
      p2,
      1
    ),
    max_bound = max_ret,
    label = "retained raw nuisance problem inputs"
  )
  retained_mask_entries <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(5, n_d, max_ret, "retained cohort_index, in_risk, in_observed, risk/observed index maps"),
      ipw_enclosure_safe_product(2, m_d, max_ret, "retained cohort support counts")
    ),
    max_bound = max_ret,
    label = "retained risk and observer mask entries"
  )
  retained_enclosure_entries <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(22, p2, max_ret, "certificate H/M/E/K variants and slices (18p^2), H_star/F_inv_T (4p^2)"),
      ipw_enclosure_safe_product(4, pm, max_ret, "retained L_0 and L_star (p x cohorts x 2)"),
      ipw_enclosure_safe_product(8, p_d, max_ret, "retained g, whitened_g, coord_dual_iv, beta_star_box"),
      ipw_enclosure_safe_product(24, m_d, max_ret, "retained cohort mass/rate/share/pi/e_L/ell enclosures"),
      64
    ),
    max_bound = max_ret,
    label = "retained nuisance and certificate enclosures"
  )
  total_retained_entries <- ipw_enclosure_safe_sum(
    c(retained_raw_entries, retained_mask_entries, retained_enclosure_entries),
    max_bound = max_ret,
    label = "total retained nuisance problem entries"
  )

  # Chunk scratch during problem validation and outcome-free Route-A/B pass
  chunk_B_by_p <- ipw_enclosure_safe_product(4, ipw_enclosure_safe_product(B_eff, p_d, max_me, "B x p chunk"), max_me, "chunk B x p scratch")
  chunk_B_by_m <- ipw_enclosure_safe_product(2, ipw_enclosure_safe_product(B_eff, m_d, max_me, "B x cohorts chunk"), max_me, "chunk B x cohorts scratch")
  chunk_vec_scratch <- ipw_enclosure_safe_product(16, B_eff, max_me, "chunk vector validation/enclosure scratch")
  matrix_scratch <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(8, p2, max_me, "p^2 triangular solve scratch"),
      ipw_enclosure_safe_product(4, pm, max_me, "p x cohorts accumulator scratch"),
      ipw_enclosure_safe_product(16, p_d, max_me, "p vector scratch"),
      ipw_enclosure_safe_product(16, m_d, max_me, "cohort vector scratch"),
      64
    ),
    max_bound = max_me,
    label = "problem matrix/vector scratch"
  )

  concurrent_entries <- ipw_enclosure_safe_sum(
    c(total_retained_entries, chunk_B_by_p, chunk_B_by_m, chunk_vec_scratch, matrix_scratch),
    max_bound = max_me,
    label = "concurrent nuisance problem double/logical payload entries"
  )

  # During raw certification its frozen, separately bounded workspace coexists
  # with retained constructor inputs. ID normalization/hash and mask validation
  # can also use full-n temporary vectors; account before either traversal.
  raw_stage_entries <- ipw_enclosure_safe_sum(c(retained_raw_entries,
    raw_pre$concurrent_entries,
    ipw_enclosure_safe_product(4, n_d, max_me, "full-n validation/hash temporary payload")),
    max_me, "raw certificate stage plus retained inputs")
  input_copy_stage_entries <- ipw_enclosure_safe_sum(c(retained_raw_entries,
    np, ipw_enclosure_safe_product(4, n_d, max_me, "input coercion/ID validation temporary payload")),
    max_me, "retained raw inputs plus Z coercion/matrix temporary")
  concurrent_entries <- max(concurrent_entries, raw_stage_entries, input_copy_stage_entries)

  # Added work operations in the nuisance problem layer (beyond raw certificate)
  per_obs_added_ops <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(4, p2, max_wo, "per-observer H_star and dual norm F^-T ops"),
      ipw_enclosure_safe_product(12, p_d, max_wo, "per-observer L_0, L_star, eta, z ops"),
      ipw_enclosure_safe_product(10, N_deg, max_wo, "per-observer Lemma exp Taylor/remainder ops"),
      ipw_enclosure_safe_product(2, S_max, max_wo, "per-observer squaring ops"),
      40
    ),
    max_bound = max_wo,
    label = "per-observer added nuisance interval operations"
  )
  pass_added_ops <- ipw_enclosure_safe_product(n_d, per_obs_added_ops, max_wo, "streamed observer added interval operations")
  p3 <- ipw_enclosure_safe_product(p2, p_d, max_wo, "p^3 F^-T basis operations")
  post_added_ops <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(4, p3, max_wo, "F^-T basis and dual norm solves"),
      ipw_enclosure_safe_product(6, ipw_enclosure_safe_product(m_d, p2, max_wo, "m * p^2 L dual norms"), max_wo, "L_0 and L_star F^-T solves"),
      ipw_enclosure_safe_product(8, pm, max_wo, "L normalization and norm ops"),
      ipw_enclosure_safe_product(16, m_d, max_wo, "cohort rate/pi/bound ops"),
      ipw_enclosure_safe_product(8, p_d, max_wo, "coordinate box projection ops"),
      64
    ),
    max_bound = max_wo,
    label = "post-loop nuisance problem interval operations"
  )
  added_work_ops_bound <- ipw_enclosure_safe_sum(
    c(ipw_enclosure_safe_product(2, n_d, max_wo, "mask/weight pass ops"), pass_added_ops, post_added_ops),
    max_bound = max_wo,
    label = "total added nuisance problem interval operations"
  )
  total_work_ops_bound <- ipw_enclosure_safe_sum(
    c(raw_pre$total_preflight_ops_bound, added_work_ops_bound), max_wo,
    "raw certificate plus nuisance transfer interval work")

  list(
    n_a_units = n_d,
    p_covariates = p_d,
    m_cohorts = m_d,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(ceiling(n_d / B_eff)),
    retained_raw_entries = retained_raw_entries,
    retained_mask_entries = retained_mask_entries,
    retained_enclosure_entries = retained_enclosure_entries,
    total_retained_entries = total_retained_entries,
    chunk_B_by_p_entries = chunk_B_by_p,
    chunk_B_by_cohorts_entries = chunk_B_by_m,
    p_squared_entries = p2,
    p_by_cohorts_entries = pm,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    added_work_ops_bound = added_work_ops_bound,
    total_work_ops_bound = total_work_ops_bound,
    raw_stage_entries = raw_stage_entries,
    input_copy_stage_entries = input_copy_stage_entries,
    raw_certificate_preflight = raw_pre,
    peak_rss_promised = FALSE,
    accounting_scope = "8-byte numeric/logical payload entries and conservative helper-boundary interval operation charges; excludes R SEXP object headers and process peak RSS."
  )
}

preflight_ipw_link_enclosure <- function(n, p, m, control, certificate_control, has_candidates = FALSE) {
  ipw_enclosure_require_deps()
  ctrl <- validate_ipw_link_enclosure_control(control)

  n_d <- ipw_enclosure_int_bound(n, "n_a_units", min_val = 1)
  p_d <- ipw_enclosure_int_bound(p, "p_covariates", min_val = 1)
  m_d <- ipw_enclosure_int_bound(m, "m_cohorts", min_val = 2)

  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  B_eff <- min(n_d, ctrl$max_chunk_rows)

  p2 <- ipw_enclosure_safe_product(p_d, p_d, max_me, "p x p matrix block")
  pm <- ipw_enclosure_safe_product(p_d, m_d, max_me, "p x cohorts block")

  # Retained per-link state: observed endpoint values (at most n each for y_left, y_right)
  # plus outcome sufficient statistics, loadings, and Theorem B bounds (no n x links IF cache).
  retained_endpoint_entries <- ipw_enclosure_safe_product(3, n_d, max_ret, "two observed endpoint vectors and all-A observer-position map")
  retained_state_entries <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(4, p2, max_ret, "link H_0, H_star references"),
      ipw_enclosure_safe_product(4, pm, max_ret, "link L_0, L_star references"),
      ipw_enclosure_safe_product(32, p_d, max_ret, "link v, G_beta, residual/coordinate balls, candidates and separated errors"),
      ipw_enclosure_safe_product(48, m_d, max_ret, "link gamma, lambda_r, rate_coef, Theorem B and separated cohort vectors"),
      if (isTRUE(has_candidates)) ipw_enclosure_safe_sum(c(ipw_enclosure_safe_product(4, pm, max_ret, "candidate L and two entrywise error matrices"), 8 * p_d, 8 * m_d, 32), max_ret, "candidate storage") else 0,
      64
    ),
    max_bound = max_ret,
    label = "retained per-link sufficient statistics and loading enclosures"
  )
  total_retained_entries <- ipw_enclosure_safe_sum(
    c(retained_endpoint_entries, retained_state_entries),
    max_bound = max_ret,
    label = "total per-link retained entries"
  )
  if (retained_state_entries > max_out) {
    cdid_ipw_enclosure_abort(
      paste0("Link summary output entries (", retained_state_entries, ") exceed `max_output_entries` (", max_out, ")."),
      "cdid_ipw_enclosure_bound_error",
      diagnostics = list(requested = retained_state_entries, bound = max_out)
    )
  }

  chunk_B_by_p <- ipw_enclosure_safe_product(4, ipw_enclosure_safe_product(B_eff, p_d, max_me, "B x p chunk"), max_me, "link chunk B x p scratch")
  chunk_B_by_m <- ipw_enclosure_safe_product(2, ipw_enclosure_safe_product(B_eff, m_d, max_me, "B x cohorts chunk"), max_me, "link chunk B x cohorts scratch")
  chunk_vec_scratch <- ipw_enclosure_safe_product(16, B_eff, max_me, "link chunk endpoint/weight scratch")

  concurrent_entries <- ipw_enclosure_safe_sum(
    c(total_retained_entries, chunk_B_by_p, chunk_B_by_m, chunk_vec_scratch,
      ipw_enclosure_safe_product(4, n_d, max_me, "all-A ID normalization/hash and optional-mask temporary payload"),
      4 * p2, 4 * pm, 64),
    max_bound = max_me,
    label = "concurrent link enclosure payload entries"
  )

  # Two worst-case observer passes, each including c_i, dual norm, exp,
  # endpoint subtraction and all scalar/vector accumulation. Post-loop charges
  # include both loading solves, all cohort ball dots and every Corollary sep field.
  hw <- ipw_enclosure_helper_work(p_d, certificate_control)
  per_obs_link_ops <- ipw_enclosure_safe_sum(
    c(2*hw$dual, 2*hw$exp, 16*p_d, 80), max_wo, "two-pass per-row link work")
  stream_ops <- ipw_enclosure_safe_product(n_d, per_obs_link_ops, max_wo, "two-pass link work")
  loading_stage_ops <- ipw_enclosure_safe_sum(
    c(2*hw$dual, 2*hw$solve,
      ipw_enclosure_safe_product(m_d, 2*hw$ball + 60, max_wo, "all cohort loading work"),
      4*pm, 40*p_d, 128), max_wo, "loading and separated error work")
  total_work_ops_bound <- ipw_enclosure_safe_sum(
    c(stream_ops, loading_stage_ops), max_wo, "total link interval work")

  list(
    n_a_units = n_d,
    p_covariates = p_d,
    m_cohorts = m_d,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(ceiling(n_d / B_eff)),
    retained_endpoint_entries = retained_endpoint_entries,
    retained_state_entries = retained_state_entries,
    total_retained_entries = total_retained_entries,
    output_entries = retained_state_entries,
    chunk_B_by_p_entries = chunk_B_by_p,
    chunk_B_by_cohorts_entries = chunk_B_by_m,
    p_squared_entries = p2,
    p_by_cohorts_entries = pm,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    total_work_ops_bound = total_work_ops_bound,
    dense_n_by_links_matrix_allocated = FALSE,
    retained_many_links_cost = "Two observer endpoint copies plus an n-position map per link; aggregate O(sum_l |O_l| + n*links), not shared period storage.",
    n_by_targets_cache_allocated = FALSE,
    peak_rss_promised = FALSE
  )
}

preflight_ipw_link_if_chunk <- function(
  b_len,
  p,
  m,
  control,
  certificate_control,
  include_components = TRUE,
  include_shared_xi = FALSE,
  has_phi_candidate = FALSE,
  has_xi_beta = FALSE
) {
  ipw_enclosure_require_deps()
  ctrl <- validate_ipw_link_enclosure_control(control)

  B_d <- ipw_enclosure_int_bound(b_len, "chunk_rows", min_val = 1)
  p_d <- ipw_enclosure_int_bound(p, "p_covariates", min_val = 1)
  m_d <- ipw_enclosure_int_bound(m, "m_cohorts", min_val = 2)

  if (B_d > ctrl$max_chunk_rows) {
    cdid_ipw_enclosure_abort(
      paste0("Requested chunk length (", B_d, ") exceeds `control$max_chunk_rows` (", ctrl$max_chunk_rows, ")."),
      "cdid_ipw_enclosure_bound_error",
      diagnostics = list(requested_chunk_rows = B_d, max_chunk_rows = ctrl$max_chunk_rows)
    )
  }

  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries

  Bp <- ipw_enclosure_safe_product(B_d, p_d, max_me, "chunk B x p")
  Bm <- ipw_enclosure_safe_product(B_d, m_d, max_me, "chunk B x cohorts")
  p2 <- ipw_enclosure_safe_product(p_d, p_d, max_me, "p x p block")
  pm <- ipw_enclosure_safe_product(p_d, m_d, max_me, "p x cohorts block")

  # Output IF entries for this chunk:
  # Base: phi_star (2B), phi_0 (2B), phi_star_via_root_bound (2B), phi_star_intersected (2B), E_root (B), metadata (8B) = 17B
  # Components: 4 terms * 2 routes * 2 endpoints (16B) + 4 E_root terms (4B) = 20B
  # Candidate sep: tilde_phi, eval_error, total_error_bound, direct_star_error = 4B
  # Shared xi: xi_rate_0 (2Bm), xi_rate_star (2Bm), u_0 (2Bp), u_star (2Bp), rhs_0 (2Bp), rhs_star (2Bp),
  #            gamma_dot_xi_r (4B), lam_dot_rhs (4B), and if has_xi_beta: xi_beta_0 (2Bp), xi_beta_star (2Bp),
  #            xi_beta_rad (2B), neg_v_dot_xi_beta (4B), phi_shared (4B).
  if_cols_per_row <- 17 + (if (isTRUE(include_components)) 20 else 0) + (if (isTRUE(has_phi_candidate)) 4 else 0)
  chunk_B_by_if <- ipw_enclosure_safe_product(if_cols_per_row, B_d, max_out, "chunk B x IF output entries")

  shared_out_entries <- 0
  if (isTRUE(include_shared_xi)) {
    shared_out_entries <- ipw_enclosure_safe_sum(
      c(
        ipw_enclosure_safe_product(4, Bm, max_out, "shared xi_rate chunk output (B x cohorts x 2 x 2)"),
        ipw_enclosure_safe_product(8, Bp, max_out, "shared u and rhs chunk output (B x p x 2 x 2)"),
        ipw_enclosure_safe_product(16, B_d, max_out, "shared chunk metadata and rate/score action output"),
        if (isTRUE(has_xi_beta)) {
          ipw_enclosure_safe_sum(
            c(
              ipw_enclosure_safe_product(5, Bp, max_out, "shared xi_beta interval output and stored candidate centers"),
              ipw_enclosure_safe_product(10, B_d, max_out, "shared xi_beta radius and link action output")
            ),
            max_bound = max_out,
            label = "shared xi_beta chunk output"
          )
        } else 0
      ),
      max_bound = max_out,
      label = "shared xi chunk output entries"
    )
  }

  output_entries <- ipw_enclosure_safe_sum(
    c(chunk_B_by_if, shared_out_entries),
    max_bound = max_out,
    label = "total chunk output entries"
  )

  chunk_scratch <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(6, Bp, max_me, "outer/inner Z slices and candidate slice/coercion/matrix temporaries"),
      ipw_enclosure_safe_product(2, Bm, max_me, "chunk cohort scratch"),
      ipw_enclosure_safe_product(8, B_d, max_me, "chunk vector scratch"),
      ipw_enclosure_safe_product(4, p2, max_me, "chunk p^2 scratch"),
      ipw_enclosure_safe_product(4, pm, max_me, "chunk p x cohorts scratch"),
      32
    ),
    max_bound = max_me,
    label = "chunk workspace scratch"
  )

  concurrent_entries <- ipw_enclosure_safe_sum(
    c(output_entries, chunk_scratch),
    max_bound = max_me,
    label = "concurrent chunk IF enclosure entries"
  )

  hw <- ipw_enclosure_helper_work(p_d, certificate_control)
  # Base four-term IF and root-error transfer, plus optional outcome-free
  # recomputation, inverse actions, link actions and separated candidate errors.
  per_row_ops <- ipw_enclosure_safe_sum(c(
    hw$dual, hw$exp, 16*p_d, 128,
    if (isTRUE(include_shared_xi)) ipw_enclosure_safe_sum(c(
      hw$dual, hw$exp, 12*p_d, 64, 2*hw$ball,
      if (isTRUE(has_xi_beta)) 2*hw$solve + 4*p_d + 16 else 0
    ), max_wo, "optional shared actions") else 0,
    if (isTRUE(has_phi_candidate)) 5 else 0
  ), max_wo, "per-row IF and shared action work")
  work_ops_bound <- ipw_enclosure_safe_product(B_d, per_row_ops, max_wo, "chunk IF interval work")

  list(
    chunk_rows = logit_cert_count_output(B_d),
    p_covariates = p_d,
    m_cohorts = m_d,
    chunk_B_by_p_entries = Bp,
    chunk_B_by_cohorts_entries = Bm,
    chunk_B_by_if_entries = chunk_B_by_if,
    shared_xi_output_entries = shared_out_entries,
    p_squared_entries = p2,
    p_by_cohorts_entries = pm,
    output_entries = output_entries,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    work_ops_bound = work_ops_bound,
    peak_rss_promised = FALSE
  )
}

preflight_ipw_link_if_chunk_cached <- function(
  b_len, p, m, control, certificate_control,
  include_components = TRUE, include_shared_xi = FALSE,
  has_phi_candidate = FALSE, has_xi_beta = FALSE
) {
  pre <- preflight_ipw_link_if_chunk(b_len, p, m, control, certificate_control,
    include_components, include_shared_xi, has_phi_candidate, has_xi_beta)
  ctrl <- validate_ipw_link_enclosure_control(control)
  B <- as.double(b_len)
  # At most44B primitive/flag scalars, <=16B+Bp cached metadata/design,
  # and safe simultaneous replacements; headers/list pointers are excluded.
  cache_entries <- ipw_enclosure_safe_sum(c(80 * B, 2 * ipw_enclosure_safe_product(B, p)),
    ctrl$max_matrix_entries, "explicit primitive cache payload plus replacement scratch")
  pre$concurrent_entries <- ipw_enclosure_safe_sum(c(pre$concurrent_entries, cache_entries),
    ctrl$max_matrix_entries, "IF workspace plus explicit primitive cache")
  if (isTRUE(include_shared_xi)) pre$output_entries <- ipw_enclosure_safe_sum(c(pre$output_entries, cache_entries),
    ctrl$max_output_entries, "shared output including explicit primitive cache")
  copy_work <- ipw_enclosure_safe_product(4, ipw_enclosure_safe_product(B, B),
    ctrl$max_work_ops, "bounded primitive-cache list/vector copy work")
  pre$work_ops_bound <- ipw_enclosure_safe_sum(c(pre$work_ops_bound, copy_work), ctrl$max_work_ops,
    "IF arithmetic plus explicit cache copy work")
  pre$concurrent_bytes_8double <- 8 * pre$concurrent_entries
  pre$cache_payload_entries <- cache_entries
  pre$cache_copy_work_ops_bound <- copy_work
  pre
}

ipw_if_scratch_work <- function(scratch) {
  if (is.null(scratch)) 0 else scratch$copy_work_ops
}

# ---------------------------------------------------------------------------
# Internal sealed access context seam and compact-view combined preflights
# ---------------------------------------------------------------------------

.ipw_access_seal_token <- local({
  e <- new.env(parent = emptyenv())
  e$id <- "cdid_ipw_access_seal_v1"
  lockEnvironment(e, bindings = TRUE)
  e
})

ipw_seal_access_context <- function(ctx, mode = ctx$mode, phase) {
  ref <- new.env(parent = emptyenv())
  ref$mode <- as.character(mode)
  ref$phase <- as.character(phase)
  ref$seal_token <- .ipw_access_seal_token
  ctx$ref <- ref
  ctx <- structure(ctx, class = "cdid_ipw_access_context")
  ref$context <- ctx
  lockEnvironment(ref, bindings = TRUE)
  ctx
}

ipw_validate_access_context <- function(ctx, expected_phase = NULL) {
  if (missing(ctx) || is.null(ctx) || !is.list(ctx) ||
      !inherits(ctx, "cdid_ipw_access_context")) {
    cdid_ipw_enclosure_abort(
      "Internal IPW access context is missing, unsealed, or malformed.",
      c("cdid_ipw_enclosure_access_context_error", "cdid_ipw_enclosure_integrity_error", "cdid_ipw_enclosure_input_error")
    )
  }
  ref <- ctx$ref
  if (!is.environment(ref) || !environmentIsLocked(ref) ||
      !exists("seal_token", envir = ref, inherits = FALSE) ||
      !exists("context", envir = ref, inherits = FALSE) ||
      !identical(get("seal_token", envir = ref, inherits = FALSE), .ipw_access_seal_token)) {
    cdid_ipw_enclosure_abort(
      "Foreign or forged IPW access context is refused; only internally sealed contexts may enter shared workers.",
      c("cdid_ipw_enclosure_access_context_error", "cdid_ipw_enclosure_integrity_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.null(expected_phase) && !identical(get("phase", envir = ref, inherits = FALSE), expected_phase)) {
    cdid_ipw_enclosure_abort(
      "IPW access context phase mismatch.",
      c("cdid_ipw_enclosure_access_context_error", "cdid_ipw_enclosure_integrity_error", "cdid_ipw_enclosure_input_error")
    )
  }
  snap <- get("context", envir = ref, inherits = FALSE)
  if (!is.list(snap) || !identical(names(ctx), names(snap)) ||
      !all(vapply(names(ctx), function(nm) ((is.null(ctx[[nm]]) && is.null(snap[[nm]])) || gmm_sexp_aliases(ctx[[nm]], snap[[nm]])), logical(1)))) {
    cdid_ipw_enclosure_abort(
      "IPW access context was modified after internal sealing; execution is refused.",
      c("cdid_ipw_enclosure_access_context_error", "cdid_ipw_enclosure_integrity_error", "cdid_ipw_enclosure_input_error")
    )
  }
  invisible(ctx)
}

ipw_enclosure_id_bytes_for_pass <- function(panel_owner, expected_a_id, n_d, B_eff, r_ctrl) {
  if (is.null(expected_a_id)) {
    return(list(
      active = FALSE,
      max_chunk_source_bytes = 0,
      max_single_id_bytes = 0,
      total_pass_source_bytes = 0,
      max_chunk_additional_scratch_entries = 0,
      per_pass_additional_work_ops = 0
    ))
  }
  ipw_owner_refuse_decorative_attrs(expected_a_id, "expected_a_id")
  if (is.factor(expected_a_id) || is.object(expected_a_id) || !is.atomic(expected_a_id) ||
      !is.null(dim(expected_a_id)) || length(expected_a_id) != n_d) {
    cdid_ipw_owner_abort(
      "`expected_a_id` must be NULL or a plain atomic vector of length `n_a_units`.",
      c("cdid_ipw_owner_id_error", "cdid_ipw_owner_order_mismatch_error")
    )
  }
  if (!identical(panel_owner$id_type, "character") || !is.character(expected_a_id) || anyNA(expected_a_id)) {
    return(list(
      active = TRUE,
      max_chunk_source_bytes = 0,
      max_single_id_bytes = 0,
      total_pass_source_bytes = 0,
      max_chunk_additional_scratch_entries = 0,
      per_pass_additional_work_ops = 0
    ))
  }
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
  add_scratch <- ipw_owner_safe_product(
    32, max_chunk_b_id, r_ctrl$max_matrix_entries,
    "UTF8/raw/integer ID temporary slots"
  )
  add_work <- ipw_owner_safe_product(
    128, total_pass_b_id, r_ctrl$max_work_ops,
    "whole-pass character ID comparison work"
  )
  list(
    active = TRUE,
    max_chunk_source_bytes = max_chunk_b_id,
    max_single_id_bytes = max_single_b_id,
    total_pass_source_bytes = total_pass_b_id,
    max_chunk_additional_scratch_entries = add_scratch,
    per_pass_additional_work_ops = add_work
  )
}

preflight_ipw_nuisance_problem_view <- function(
  nuisance_view,
  certificate_control,
  problem_control,
  reader_control,
  expected_a_id = NULL,
  automatic_radius = FALSE
) {
  ipw_enclosure_require_deps()
  auto_r <- logit_cert_automatic_radius_flag(automatic_radius)
  cert_ctrl <- validate_logit_certificate_control(certificate_control)
  prob_ctrl <- validate_ipw_nuisance_problem_control(problem_control)
  if (is.null(expected_a_id) && !missing(reader_control) && !is.null(reader_control)) {
    expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  }
  r_ctrl <- validate_ipw_owner_control(reader_control)
  validate_ipw_nuisance_view(nuisance_view)

  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = nuisance_view$panel_ref, inherits = FALSE)

  n_d <- ipw_enclosure_int_bound(frame_owner$n_a_units, "n_a_units", min_val = 1)
  p_d <- ipw_enclosure_int_bound(frame_owner$n_covariates, "p_covariates", min_val = 1)
  m_d <- ipw_enclosure_int_bound(nuisance_view$cohorts$m_cohorts, "m_cohorts", min_val = 2)

  if (!is.null(expected_a_id)) preflight_ipw_nuisance_problem_view(nuisance_view, cert_ctrl, prob_ctrl, r_ctrl, expected_a_id = NULL, automatic_radius = auto_r)

  # 1. Raw certificate view preflight (3 reader passes + raw arithmetic; the
  #    automatic-radius selection charge is reserved only when requested)
  raw_view_pre <- preflight_logit_certificate_view(
    nuisance_view = nuisance_view,
    control = cert_ctrl,
    reader_control = r_ctrl,
    expected_a_id = expected_a_id,
    automatic_radius = auto_r
  )

  max_ret <- prob_ctrl$max_retained_entries
  max_me <- prob_ctrl$max_matrix_entries
  max_wo <- prob_ctrl$max_work_ops
  B_eff <- min(n_d, prob_ctrl$max_chunk_rows, r_ctrl$max_chunk_rows)
  n_chunks_prob <- ceiling(n_d / B_eff)
  N_deg <- cert_ctrl$taylor_degree
  S_max <- cert_ctrl$max_squarings

  p2 <- ipw_enclosure_safe_product(p_d, p_d, max_ret, "p x p matrix block")
  pm <- ipw_enclosure_safe_product(p_d, m_d, max_ret, "p x cohorts Jacobian block")

  # Compact retained state: ZERO Z copy, ZERO id/w/G/pair_observed copies, ZERO full-N masks!
  retained_raw_entries <- ipw_enclosure_safe_sum(
    c(m_d, p_d, p2, 1),
    max_bound = max_ret,
    label = "retained compact nuisance problem header"
  )
  retained_mask_entries <- ipw_enclosure_safe_product(
    2, m_d, max_ret,
    "retained cohort support counts"
  )
  retained_enclosure_entries <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(22, p2, max_ret, "certificate H/M/E/K variants and slices (18p^2), H_star/F_inv_T (4p^2)"),
      ipw_enclosure_safe_product(4, pm, max_ret, "retained L_0 and L_star (p x cohorts x 2)"),
      ipw_enclosure_safe_product(8, p_d, max_ret, "retained g, whitened_g, coord_dual_iv, beta_star_box"),
      ipw_enclosure_safe_product(24, m_d, max_ret, "retained cohort mass/rate/share/pi/e_L/ell enclosures"),
      64
    ),
    max_bound = max_ret,
    label = "retained nuisance and certificate enclosures"
  )
  total_retained_entries <- ipw_enclosure_safe_sum(
    c(retained_raw_entries, retained_mask_entries, retained_enclosure_entries),
    max_bound = max_ret,
    label = "total retained compact nuisance problem entries"
  )

  # Nuisance transfer chunk scratch (size B_eff) and matrix scratch
  chunk_B_by_p <- ipw_enclosure_safe_product(4, ipw_enclosure_safe_product(B_eff, p_d, max_me, "B x p chunk"), max_me, "chunk B x p scratch")
  chunk_B_by_m <- ipw_enclosure_safe_product(2, ipw_enclosure_safe_product(B_eff, m_d, max_me, "B x cohorts chunk"), max_me, "chunk B x cohorts scratch")
  chunk_vec_scratch <- ipw_enclosure_safe_product(16, B_eff, max_me, "chunk vector validation/enclosure scratch")
  matrix_scratch <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(8, p2, max_me, "p^2 triangular solve scratch"),
      ipw_enclosure_safe_product(4, pm, max_me, "p x cohorts accumulator scratch"),
      ipw_enclosure_safe_product(16, p_d, max_me, "p vector scratch"),
      ipw_enclosure_safe_product(16, m_d, max_me, "cohort vector scratch"),
      64
    ),
    max_bound = max_me,
    label = "problem matrix/vector scratch"
  )

  # Reader preflight for the 4th pass (nuisance transfer pass at B_eff)
  r_pre_prob_base <- preflight_ipw_range_query(
    b_len = B_eff,
    p_covariates = p_d,
    m_cohorts = m_d,
    control = r_ctrl,
    read_outcomes = FALSE,
    include_design = TRUE,
    enclose_differences = FALSE
  )
  id_prob_info <- ipw_enclosure_id_bytes_for_pass(panel_owner, expected_a_id, n_d, B_eff, r_ctrl)
  r_pre_prob <- if (id_prob_info$max_chunk_source_bytes > 0) {
    ipw_owner_preflight_id_bytes(
      r_pre_prob_base, r_ctrl,
      list(total = id_prob_info$max_chunk_source_bytes, max = id_prob_info$max_single_id_bytes),
      owner_max_bytes = 0, lookup_count = 0, search_width = 1, simultaneous = TRUE
    )
  } else {
    r_pre_prob_base
  }

  prob_reader_pass_ops <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks_prob, r_pre_prob_base$work_ops_bound, r_ctrl$max_work_ops, "nuisance transfer reader chunk work"),
      id_prob_info$per_pass_additional_work_ops
    ),
    max_bound = r_ctrl$max_work_ops,
    label = "nuisance transfer pass reader work operations"
  )
  total_reader_work_ops <- ipw_owner_safe_sum(
    c(raw_view_pre$total_reader_work_ops_bound, prob_reader_pass_ops),
    max_bound = r_ctrl$max_work_ops,
    label = "total 4-pass nuisance problem reader work operations"
  )

  # Simultaneous stage payload sums checked against prob_ctrl$max_matrix_entries:
  raw_stage_entries <- ipw_enclosure_safe_sum(
    c(retained_raw_entries, raw_view_pre$combined_concurrent_entries),
    max_bound = max_me,
    label = "raw certificate stage plus compact header"
  )
  transfer_arith_concurrent <- ipw_enclosure_safe_sum(
    c(total_retained_entries, chunk_B_by_p, chunk_B_by_m, chunk_vec_scratch, matrix_scratch),
    max_bound = max_me,
    label = "compact nuisance transfer arithmetic payload entries"
  )
  transfer_stage_entries <- ipw_enclosure_safe_sum(
    c(transfer_arith_concurrent, r_pre_prob$concurrent_entries),
    max_bound = max_me,
    label = "simultaneous compact nuisance transfer arithmetic and reader payload entries"
  )
  matrix_stage_entries <- ipw_enclosure_safe_sum(
    c(total_retained_entries, matrix_scratch),
    max_bound = max_me,
    label = "compact nuisance matrix stage payload entries"
  )
  combined_concurrent_entries <- max(raw_stage_entries, transfer_stage_entries, matrix_stage_entries)

  # Added arithmetic work operations in nuisance transfer
  per_obs_added_ops <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(4, p2, max_wo, "per-observer H_star and dual norm F^-T ops"),
      ipw_enclosure_safe_product(12, p_d, max_wo, "per-observer L_0, L_star, eta, z ops"),
      ipw_enclosure_safe_product(10, N_deg, max_wo, "per-observer Lemma exp Taylor/remainder ops"),
      ipw_enclosure_safe_product(2, S_max, max_wo, "per-observer squaring ops"),
      40
    ),
    max_bound = max_wo,
    label = "per-observer added nuisance interval operations"
  )
  pass_added_ops <- ipw_enclosure_safe_product(n_d, per_obs_added_ops, max_wo, "streamed observer added interval operations")
  p3 <- ipw_enclosure_safe_product(p2, p_d, max_wo, "p^3 F^-T basis operations")
  post_added_ops <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(4, p3, max_wo, "F^-T basis and dual norm solves"),
      ipw_enclosure_safe_product(6, ipw_enclosure_safe_product(m_d, p2, max_wo, "m * p^2 L dual norms"), max_wo, "L_0 and L_star F^-T solves"),
      ipw_enclosure_safe_product(8, pm, max_wo, "L normalization and norm ops"),
      ipw_enclosure_safe_product(16, m_d, max_wo, "cohort rate/pi/bound ops"),
      ipw_enclosure_safe_product(8, p_d, max_wo, "coordinate box projection ops"),
      64
    ),
    max_bound = max_wo,
    label = "post-loop nuisance problem interval operations"
  )
  added_work_ops_bound <- ipw_enclosure_safe_sum(
    c(pass_added_ops, post_added_ops),
    max_bound = max_wo,
    label = "total added compact nuisance problem interval operations"
  )
  total_arith_work_ops_bound <- ipw_enclosure_safe_sum(
    c(raw_view_pre$total_preflight_ops_bound, added_work_ops_bound),
    max_bound = max_wo,
    label = "raw certificate plus nuisance transfer interval work"
  )
  combined_work_ops_bound <- ipw_enclosure_safe_sum(
    c(total_arith_work_ops_bound, total_reader_work_ops),
    max_bound = max_wo,
    label = "combined raw certificate, nuisance transfer, and 4-pass reader work operations"
  )

  list(
    storage_mode = "compact_view",
    n_a_units = n_d,
    p_covariates = p_d,
    m_cohorts = m_d,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(n_chunks_prob),
    n_reader_passes = 4L,
    retained_raw_entries = retained_raw_entries,
    retained_mask_entries = retained_mask_entries,
    retained_enclosure_entries = retained_enclosure_entries,
    total_retained_entries = total_retained_entries,
    chunk_B_by_p_entries = chunk_B_by_p,
    chunk_B_by_cohorts_entries = chunk_B_by_m,
    p_squared_entries = p2,
    p_by_cohorts_entries = pm,
    concurrent_entries = transfer_arith_concurrent,
    concurrent_bytes_8double = 8 * transfer_arith_concurrent,
    reader_concurrent_entries = r_pre_prob$concurrent_entries,
    combined_concurrent_entries = combined_concurrent_entries,
    combined_concurrent_bytes_8double = 8 * combined_concurrent_entries,
    stage_concurrent_entries = list(
      raw_certificate_stage = raw_stage_entries,
      transfer_chunk_stage = transfer_stage_entries,
      matrix_stage = matrix_stage_entries
    ),
    automatic_radius_reserved = raw_view_pre$automatic_radius_reserved,
    automatic_radius_ops_bound = raw_view_pre$automatic_radius_ops_bound,
    automatic_radius_entries = raw_view_pre$automatic_radius_entries,
    added_work_ops_bound = added_work_ops_bound,
    total_work_ops_bound = total_arith_work_ops_bound,
    transfer_reader_work_ops_bound = prob_reader_pass_ops,
    total_reader_work_ops_bound = total_reader_work_ops,
    combined_work_ops_bound = combined_work_ops_bound,
    raw_stage_entries = raw_stage_entries,
    input_copy_stage_entries = 0,
    raw_certificate_preflight = raw_view_pre,
    transfer_reader_preflight = r_pre_prob,
    id_byte_preflight = id_prob_info,
    full_n_masks_retained = FALSE,
    copied_Z_retained = FALSE,
    endpoints_read = 0L,
    peak_rss_promised = FALSE,
    accounting_scope = "Compact O(p^2 + p*m + m) retained nuisance state plus combined arithmetic and 4-pass bounded reader concurrent/work charges; zero Z copy or full-N masks retained."
  )
}

preflight_ipw_link_enclosure_view <- function(
  problem,
  link_view,
  control,
  reader_control,
  has_candidates = FALSE,
  expected_a_id = NULL
) {
  ipw_enclosure_require_deps()
  ctrl <- validate_ipw_link_enclosure_control(control)
  if (is.null(expected_a_id) && !missing(reader_control) && !is.null(reader_control)) {
    expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  }
  r_ctrl <- validate_ipw_owner_control(reader_control)
  ipw_validate_nuisance_problem_mode(problem, require_certified = FALSE, expected_storage_mode = "compact_view")
  frame_owner <- get("frame_owner", envir = problem$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = problem$panel_ref, inherits = FALSE)
  nuisance_view <- get("nuisance_view", envir = problem$nuisance_ref, inherits = FALSE)
  validate_ipw_link_view(link_view, frame_owner = frame_owner, nuisance_view = nuisance_view)

  n_d <- ipw_enclosure_int_bound(problem$raw$n_a_units, "n_a_units", min_val = 1)
  p_d <- ipw_enclosure_int_bound(problem$raw$p_covariates, "p_covariates", min_val = 1)
  m_d <- ipw_enclosure_int_bound(length(problem$raw$cohorts), "m_cohorts", min_val = 2)

  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  B_eff <- min(n_d, ctrl$max_chunk_rows, r_ctrl$max_chunk_rows)
  n_chunks <- ceiling(n_d / B_eff)

  p2 <- ipw_enclosure_safe_product(p_d, p_d, max_me, "p x p matrix block")
  pm <- ipw_enclosure_safe_product(p_d, m_d, max_me, "p x cohorts block")

  # Compact per-link state: ZERO retained endpoint vectors, ZERO retained N-position map!
  retained_endpoint_entries <- 0
  retained_state_entries <- ipw_enclosure_safe_sum(
    c(
      ipw_enclosure_safe_product(4, p2, max_ret, "link H_0, H_star references"),
      ipw_enclosure_safe_product(4, pm, max_ret, "link L_0, L_star references"),
      ipw_enclosure_safe_product(32, p_d, max_ret, "link v, G_beta, residual/coordinate balls, candidates and separated errors"),
      ipw_enclosure_safe_product(48, m_d, max_ret, "link gamma, lambda_r, rate_coef, Theorem B and separated cohort vectors"),
      if (isTRUE(has_candidates)) ipw_enclosure_safe_sum(c(ipw_enclosure_safe_product(4, pm, max_ret, "candidate L and two entrywise error matrices"), 8 * p_d, 8 * m_d, 32), max_ret, "candidate storage") else 0,
      64
    ),
    max_bound = max_ret,
    label = "retained per-link sufficient statistics and loading enclosures"
  )
  total_retained_entries <- retained_state_entries
  if (retained_state_entries > max_out) {
    cdid_ipw_enclosure_abort(
      paste0("Link summary output entries (", retained_state_entries, ") exceed `max_output_entries` (", max_out, ")."),
      "cdid_ipw_enclosure_bound_error",
      diagnostics = list(requested = retained_state_entries, bound = max_out)
    )
  }

  chunk_B_by_p <- ipw_enclosure_safe_product(4, ipw_enclosure_safe_product(B_eff, p_d, max_me, "B x p chunk"), max_me, "link chunk B x p scratch")
  chunk_B_by_m <- ipw_enclosure_safe_product(2, ipw_enclosure_safe_product(B_eff, m_d, max_me, "B x cohorts chunk"), max_me, "link chunk B x cohorts scratch")
  chunk_vec_scratch <- ipw_enclosure_safe_product(16, B_eff, max_me, "link chunk endpoint/weight scratch")

  arith_concurrent_entries <- ipw_enclosure_safe_sum(
    c(total_retained_entries, chunk_B_by_p, chunk_B_by_m, chunk_vec_scratch, 4 * p2, 4 * pm, 64),
    max_bound = max_me,
    label = "compact link arithmetic concurrent payload entries"
  )

  r_pre_base <- preflight_ipw_range_query(
    b_len = B_eff,
    p_covariates = p_d,
    m_cohorts = m_d,
    control = r_ctrl,
    read_outcomes = TRUE,
    include_design = TRUE,
    enclose_differences = FALSE
  )
  if (!is.null(expected_a_id)) preflight_ipw_link_enclosure_view(problem, link_view, ctrl, r_ctrl, has_candidates, expected_a_id = NULL)
  id_info <- ipw_enclosure_id_bytes_for_pass(panel_owner, expected_a_id, n_d, B_eff, r_ctrl)
  r_pre_chunk <- if (id_info$max_chunk_source_bytes > 0) {
    ipw_owner_preflight_id_bytes(
      r_pre_base, r_ctrl,
      list(total = id_info$max_chunk_source_bytes, max = id_info$max_single_id_bytes),
      owner_max_bytes = 0, lookup_count = 0, search_width = 1, simultaneous = TRUE
    )
  } else {
    r_pre_base
  }

  pass_reader_ops <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, r_pre_base$work_ops_bound, r_ctrl$max_work_ops, "link pass reader chunk work"),
      id_info$per_pass_additional_work_ops
    ),
    max_bound = r_ctrl$max_work_ops,
    label = "single-pass link reader work operations"
  )
  total_reader_work_ops <- ipw_owner_safe_sum(
    c(pass_reader_ops, pass_reader_ops),
    max_bound = r_ctrl$max_work_ops,
    label = "total 2-pass link reader work operations"
  )

  stage_pass_concurrent <- ipw_enclosure_safe_sum(
    c(arith_concurrent_entries, r_pre_chunk$concurrent_entries),
    max_bound = max_me,
    label = "simultaneous compact link arithmetic and reader concurrent payload entries"
  )
  stage_loading_concurrent <- ipw_enclosure_safe_sum(
    c(total_retained_entries, 4 * p2, 4 * pm, 64),
    max_bound = max_me,
    label = "compact link loading stage concurrent payload entries"
  )
  combined_concurrent_entries <- max(stage_pass_concurrent, stage_loading_concurrent)

  certificate_control <- problem$controls$certificate_control
  hw <- ipw_enclosure_helper_work(p_d, certificate_control)
  per_obs_link_ops <- ipw_enclosure_safe_sum(
    c(2 * hw$dual, 2 * hw$exp, 16 * p_d, 80), max_wo, "two-pass per-row link work"
  )
  stream_ops <- ipw_enclosure_safe_product(n_d, per_obs_link_ops, max_wo, "two-pass link work")
  loading_stage_ops <- ipw_enclosure_safe_sum(
    c(
      2 * hw$dual, 2 * hw$solve,
      ipw_enclosure_safe_product(m_d, 2 * hw$ball + 60, max_wo, "all cohort loading work"),
      4 * pm, 40 * p_d, 128
    ),
    max_wo, "loading and separated error work"
  )
  total_work_ops_bound <- ipw_enclosure_safe_sum(
    c(stream_ops, loading_stage_ops), max_wo, "total link interval work"
  )
  combined_work_ops_bound <- ipw_enclosure_safe_sum(
    c(total_work_ops_bound, total_reader_work_ops),
    max_bound = max_wo,
    label = "combined 2-pass link interval and reader work"
  )

  list(
    storage_mode = "compact_view",
    n_a_units = n_d,
    p_covariates = p_d,
    m_cohorts = m_d,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(n_chunks),
    n_reader_passes = 2L,
    retained_endpoint_entries = retained_endpoint_entries,
    retained_position_map_entries = 0,
    retained_state_entries = retained_state_entries,
    total_retained_entries = total_retained_entries,
    output_entries = retained_state_entries,
    chunk_B_by_p_entries = chunk_B_by_p,
    chunk_B_by_cohorts_entries = chunk_B_by_m,
    p_squared_entries = p2,
    p_by_cohorts_entries = pm,
    concurrent_entries = arith_concurrent_entries,
    concurrent_bytes_8double = 8 * arith_concurrent_entries,
    reader_concurrent_entries = r_pre_chunk$concurrent_entries,
    combined_concurrent_entries = combined_concurrent_entries,
    combined_concurrent_bytes_8double = 8 * combined_concurrent_entries,
    stage_concurrent_entries = list(
      pass1 = stage_pass_concurrent,
      pass2 = stage_pass_concurrent,
      loading_stage = stage_loading_concurrent
    ),
    total_work_ops_bound = total_work_ops_bound,
    reader_pass_work_ops_bound = pass_reader_ops,
    total_reader_work_ops_bound = total_reader_work_ops,
    combined_work_ops_bound = combined_work_ops_bound,
    reader_chunk_preflight = r_pre_chunk,
    id_byte_preflight = id_info,
    dense_n_by_links_matrix_allocated = FALSE,
    per_link_endpoint_copies_retained = FALSE,
    full_n_position_map_retained = FALSE,
    retained_many_links_cost = "O(p^2 + p*m + m) sufficient statistics and loading enclosures per link; zero endpoint copies or N-position maps retained.",
    n_by_targets_cache_allocated = FALSE,
    peak_rss_promised = FALSE
  )
}

preflight_ipw_link_if_view_chunk <- function(
  b_len,
  p,
  m,
  control,
  reader_control,
  certificate_control,
  include_components = TRUE,
  include_shared_xi = FALSE,
  has_phi_candidate = FALSE,
  has_xi_beta = FALSE,
  read_outcomes = TRUE,
  panel_owner = NULL,
  a_seq = NULL,
  expected_a_id = NULL,
  caller_phi_entries = 0,
  caller_xi_entries = 0
) {
  ipw_enclosure_require_deps()
  ctrl <- validate_ipw_link_enclosure_control(control)
  r_ctrl <- validate_ipw_owner_control(reader_control)

  B_d <- ipw_enclosure_int_bound(b_len, "chunk_rows", min_val = 1)
  if (B_d > ctrl$max_chunk_rows) {
    cdid_ipw_enclosure_abort(
      paste0("Requested chunk length (", B_d, ") exceeds `control$max_chunk_rows` (", ctrl$max_chunk_rows, ")."),
      "cdid_ipw_enclosure_bound_error",
      diagnostics = list(requested_chunk_rows = B_d, max_chunk_rows = ctrl$max_chunk_rows)
    )
  }
  if (B_d > r_ctrl$max_chunk_rows) {
    cdid_ipw_owner_abort(
      paste0("Requested contiguous A-range length (", B_d, ") exceeds `reader_control$max_chunk_rows` (", r_ctrl$max_chunk_rows, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error"),
      diagnostics = list(requested_chunk_rows = B_d, max_chunk_rows = r_ctrl$max_chunk_rows)
    )
  }

  arith_pre <- preflight_ipw_link_if_chunk_cached(
    b_len = B_d,
    p = p,
    m = m,
    control = ctrl,
    certificate_control = certificate_control,
    include_components = include_components,
    include_shared_xi = include_shared_xi,
    has_phi_candidate = has_phi_candidate,
    has_xi_beta = has_xi_beta
  )

  r_pre_base <- preflight_ipw_range_query(
    b_len = B_d,
    p_covariates = p,
    m_cohorts = m,
    control = r_ctrl,
    read_outcomes = read_outcomes,
    include_design = TRUE,
    enclose_differences = FALSE
  )

  ipw_enclosure_safe_sum(c(arith_pre$concurrent_entries, r_pre_base$concurrent_entries), ctrl$max_matrix_entries,
    "base simultaneous arithmetic and reader payload before optional ID metadata")
  ipw_enclosure_safe_sum(c(arith_pre$work_ops_bound, r_pre_base$work_ops_bound), ctrl$max_work_ops,
    "base combined work before optional ID metadata")
  r_pre <- r_pre_base
  if (!is.null(expected_a_id) && !is.null(panel_owner) && !is.null(a_seq)) {
    exp_sub <- if (length(expected_a_id) == B_d) expected_a_id else expected_a_id[a_seq]
    if (identical(panel_owner$id_type, "character") && is.character(exp_sub) && !anyNA(exp_sub)) {
      exp_bytes <- ipw_owner_id_byte_metadata(exp_sub)
      own_bytes <- ipw_owner_id_byte_metadata(panel_owner$a_id_key[a_seq])
      merged <- list(
        total = ipw_owner_safe_sum(c(exp_bytes$total, own_bytes$total), label = "paired expected ID bytes"),
        max = max(exp_bytes$max, own_bytes$max)
      )
      r_pre <- ipw_owner_preflight_id_bytes(
        r_pre_base, r_ctrl, merged,
        owner_max_bytes = 0, lookup_count = 0, search_width = 1, simultaneous = TRUE
      )
    }
  }

  disclosed_caller_entries <- ipw_enclosure_safe_sum(
    c(caller_phi_entries, caller_xi_entries),
    label = "disclosed caller-owned candidate entries"
  )

  combined_concurrent_entries <- ipw_enclosure_safe_sum(
    c(arith_pre$concurrent_entries, r_pre$concurrent_entries),
    max_bound = ctrl$max_matrix_entries,
    label = "combined chunk IF arithmetic and reader concurrent payload entries"
  )
  combined_work_ops_bound <- ipw_enclosure_safe_sum(
    c(arith_pre$work_ops_bound, r_pre$work_ops_bound),
    max_bound = ctrl$max_work_ops,
    label = "combined chunk IF arithmetic and reader work operations"
  )

  c(
    arith_pre,
    list(
      storage_mode = "compact_view",
      reader_base_work_ops_bound = r_pre_base$work_ops_bound,
      reader_work_ops_bound = r_pre$work_ops_bound,
      reader_concurrent_entries = r_pre$concurrent_entries,
      combined_concurrent_entries = combined_concurrent_entries,
      combined_concurrent_bytes_8double = 8 * combined_concurrent_entries,
      arith_work_ops_bound = arith_pre$work_ops_bound,
      combined_work_ops_bound = combined_work_ops_bound,
      caller_supplied_candidate_entries = disclosed_caller_entries,
      caller_candidate_retained_by_result = FALSE,
      reader_preflight = r_pre
    )
  )
}

# ---------------------------------------------------------------------------
# Internal validated interval algebra helpers
# ---------------------------------------------------------------------------

ipw_enclosure_validate_id <- function(id, n, step, label = "id") {
  if (missing(id) || is.null(id) || is.factor(id) || is.object(id) ||
      !is.atomic(id) || !is.null(dim(id)) || length(id) != n) {
    cdid_ipw_enclosure_abort(
      paste0("`", label, "` must be a plain non-object atomic vector (numeric or character) of length ", n, "."),
      c("cdid_ipw_enclosure_id_error", "cdid_ipw_enclosure_input_error")
    )
  }
  is_num <- is.numeric(id)
  is_chr <- is.character(id)
  if (!is_num && !is_chr) {
    cdid_ipw_enclosure_abort(
      paste0("`", label, "` must be finite numeric or nonempty character identifiers."),
      c("cdid_ipw_enclosure_id_error", "cdid_ipw_enclosure_input_error")
    )
  }
  a_start <- 1L
  while (a_start <= n) {
    ipw_enclosure_check_interrupt()
    a_end <- min(as.double(n), as.double(a_start) + step - 1)
    idx <- seq.int(a_start, a_end)
    ch <- id[idx]
    if (is_num) {
      if (anyNA(ch) || any(!is.finite(ch))) {
        cdid_ipw_enclosure_abort(
          paste0("`", label, "` contains NA, NaN, or infinite numeric identifiers."),
          c("cdid_ipw_enclosure_id_error", "cdid_ipw_enclosure_input_error"),
          diagnostics = list(chunk_start = a_start, chunk_end = a_end)
        )
      }
    } else {
      if (anyNA(ch) || any(!nzchar(ch))) {
        cdid_ipw_enclosure_abort(
          paste0("`", label, "` contains NA or empty character identifiers."),
          c("cdid_ipw_enclosure_id_error", "cdid_ipw_enclosure_input_error"),
          diagnostics = list(chunk_start = a_start, chunk_end = a_end)
        )
      }
    }
    a_start <- as.double(a_end) + 1
  }
  id_key <- if (is_num) {
    d <- as.double(id)
    d[d == 0.0] <- 0.0
    d
  } else {
    as.character(id)
  }
  if (anyDuplicated(id_key)) {
    cdid_ipw_enclosure_abort(
      paste0("Duplicate unit identifiers in `", label, "` are forbidden."),
      c("cdid_ipw_enclosure_id_error", "cdid_ipw_enclosure_input_error")
    )
  }
  list(
    id_copy = if (is_num) (if (is.integer(id)) as.integer(id) else id_key) else id_key,
    id_key = id_key,
    id_type = if (is_num) "numeric" else "character",
    id_storage = typeof(id)
  )
}

ipw_enclosure_check_optional_masks <- function(
  n,
  in_risk,
  in_observed,
  cohort_index,
  pair_observed_ref,
  risk_mask = NULL,
  observed_mask = NULL,
  cohort_index_arg = NULL,
  pair_observed_arg = NULL
) {
  if (!is.null(risk_mask)) {
    if (!is.atomic(risk_mask) || !is.null(dim(risk_mask)) || is.object(risk_mask) ||
        !is.logical(risk_mask) || length(risk_mask) != n || anyNA(risk_mask) ||
        !identical(as.logical(risk_mask), in_risk)) {
      cdid_ipw_enclosure_abort(
        "Supplied `risk_mask` disagrees with derived risk set R = {G in cohorts}.",
        c("cdid_ipw_enclosure_mask_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (!is.null(observed_mask)) {
    if (!is.atomic(observed_mask) || !is.null(dim(observed_mask)) || is.object(observed_mask) ||
        !is.logical(observed_mask) || length(observed_mask) != n || anyNA(observed_mask) ||
        !identical(as.logical(observed_mask), in_observed)) {
      cdid_ipw_enclosure_abort(
        "Supplied `observed_mask` disagrees with derived observer set O = R intersect S.",
        c("cdid_ipw_enclosure_mask_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (!is.null(cohort_index_arg)) {
    if (!is.atomic(cohort_index_arg) || !is.null(dim(cohort_index_arg)) || is.object(cohort_index_arg) ||
        !is.numeric(cohort_index_arg) || length(cohort_index_arg) != n ||
        !identical(as.integer(cohort_index_arg), cohort_index)) {
      cdid_ipw_enclosure_abort(
        "Supplied `cohort_index` disagrees with `match(G, cohorts)`.",
        c("cdid_ipw_enclosure_mask_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (!is.null(pair_observed_arg)) {
    if (!is.atomic(pair_observed_arg) || !is.null(dim(pair_observed_arg)) || is.object(pair_observed_arg) ||
        !is.logical(pair_observed_arg) || length(pair_observed_arg) != n || anyNA(pair_observed_arg)) {
      cdid_ipw_enclosure_abort(
        "Supplied `pair_observed` must be a non-NA atomic logical vector of length n.",
        c("cdid_ipw_enclosure_mask_error", "cdid_ipw_enclosure_input_error")
      )
    }
    if (!identical(as.logical(pair_observed_arg[in_risk]), pair_observed_ref[in_risk])) {
      cdid_ipw_enclosure_abort(
        "Link `pair_observed` mask on the risk set R = {G in cohorts} disagrees with the certified nuisance problem.",
        c("cdid_ipw_enclosure_mask_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  invisible(TRUE)
}

# Enclose c_i = w_i * M_j / Pi_j outward from raw w_i and cohort mass intervals
ipw_enclosure_c_weight <- function(wu_iv, M_j_iv, Pi_j_iv, M_over_Pi_j_iv, unit_offset, cohort_index) {
  ops <- 3 # Conservative charge for both attempted algebraic routes.
  c_iv1 <- if (!anyNA(M_over_Pi_j_iv)) {
    tryCatch(
      gmm_iv_mul(wu_iv, M_over_Pi_j_iv),
      cdid_gmm_accuracy_error = function(e) NULL,
      cdid_logit_arithmetic_inconclusive = function(e) NULL
    )
  } else NULL

  wM_try <- tryCatch(
    gmm_iv_mul(wu_iv, M_j_iv),
    cdid_gmm_accuracy_error = function(e) NULL,
    cdid_logit_arithmetic_inconclusive = function(e) NULL
  )
  c_iv2 <- if (!is.null(wM_try)) {
    tryCatch(
      gmm_iv_div(wM_try, Pi_j_iv),
      cdid_gmm_accuracy_error = function(e) NULL,
      cdid_logit_arithmetic_inconclusive = function(e) NULL
    )
  } else NULL

  if (!is.null(c_iv1) && !is.null(c_iv2)) {
    c_iv <- logit_iv_intersect(c_iv1, c_iv2, "c_i = w_i * M_j / Pi_j")
  } else if (!is.null(c_iv1)) {
    c_iv <- c_iv1
  } else if (!is.null(c_iv2)) {
    c_iv <- c_iv2
  } else {
    ipw_enclosure_arith_refuse(
      "Fitting weight c_i = w_i * M_j / Pi_j overflowed binary64.",
      reason = "fitting_weight_overflow",
      diagnostics = list(unit_offset = unit_offset, cohort_index = cohort_index)
    )
  }
  if (c_iv[1L] <= 0.0) {
    ipw_enclosure_arith_refuse(
      "Fitting weight lower endpoint c_i underflowed to <= 0.",
      reason = "fitting_weight_underflow",
      diagnostics = list(unit_offset = unit_offset, c_iv = c_iv)
    )
  }
  list(c_iv = c_iv, ops = ops)
}

# Outward enclosure of F^{-T} and coordinate dual norms ||F^{-T} e_k||_2 = sqrt((M^{-1})_{kk})
ipw_enclosure_ft_inv_basis <- function(F_dbl) {
  p <- nrow(F_dbl)
  F_inv_T <- array(0.0, dim = c(p, p, 2L))
  coord_dual_iv <- matrix(0.0, nrow = p, ncol = 2L)
  ops <- 0
  solves <- 0L
  for (k in seq_len(p)) {
    ek <- numeric(p)
    ek[k] <- 1.0
    sv <- logit_iv_forwardsolve_ft(F_dbl, ek)
    F_inv_T[, k, ] <- sv$y
    ops <- ops + sv$ops
    solves <- solves + 1L
    nrm <- gmm_iv_norm(sv$y[, 1L], sv$y[, 2L])
    coord_dual_iv[k, ] <- nrm
    ops <- ops + 3 * p + 2
  }
  list(F_inv_T = F_inv_T, coord_dual_iv = coord_dual_iv, ops = ops, solves = solves)
}

# Outward upper bound on sup_{v in v_iv} ||F^{-T} v||_2 via interval forward substitution
# intersected with explicit F^{-T} interval product (Lemma dual & Lemma incl)
ipw_enclosure_dual_norm_iv <- function(F_dbl, F_inv_T, v_iv) {
  p <- nrow(F_dbl)
  if (is.null(dim(v_iv))) {
    v_iv <- cbind(as.double(v_iv), as.double(v_iv))
  }
  sv <- logit_iv_forwardsolve_ft(F_dbl, v_iv)
  y_iv <- sv$y
  ops <- sv$ops

  # Second route via pre-enclosed F^{-T} rows
  if (!is.null(F_inv_T)) {
    for (a in seq_len(p)) {
      s_a <- c(0.0, 0.0)
      for (l in seq_len(p)) {
        f_al <- F_inv_T[a, l, ]
        if (f_al[1L] != 0.0 || f_al[2L] != 0.0) {
          s_a <- gmm_iv_add(s_a, gmm_iv_mul(f_al, v_iv[l, ]))
          ops <- ops + 2
        }
      }
      y_iv[a, ] <- logit_iv_intersect(y_iv[a, ], s_a, "F^{-T} v")
    }
  }
  nrm <- gmm_iv_norm(y_iv[, 1L], y_iv[, 2L])
  ops <- ops + 3 * p + 2
  list(y_iv = y_iv, norm_iv = nrm, upper = nrm[2L], ops = ops, solves = 1L)
}

# Widen an interval [lo, hi] by +/- rad (rad >= 0) with outward binary64 rounding
ipw_enclosure_widen_iv <- function(iv, rad) {
  if (!is.finite(rad) || rad < 0.0) {
    ipw_enclosure_arith_refuse("Nonfinite or negative radius in interval widening.", reason = "invalid_widen_radius")
  }
  if (rad == 0.0) return(iv)
  r_iv <- gmm_iv(as.double(rad))
  lo <- logit_b64_sub(gmm_iv(iv[1L]), r_iv)[1L]
  hi <- logit_b64_add(gmm_iv(iv[2L]), r_iv)[2L]
  c(lo, hi)
}

# Distance sup_{y in iv} |x - y| for stored scalar x and interval iv = [lo, hi]
ipw_enclosure_iv_dist <- function(iv, x) {
  if (!gmm_scalar_numeric(x)) {
    cdid_ipw_enclosure_abort("Candidate scalar for interval distance must be finite numeric.", "cdid_ipw_enclosure_candidate_error")
  }
  d_iv <- gmm_iv_abs(logit_b64_sub(gmm_iv(as.double(x)), iv))
  d_iv[2L]
}

# Lemma dual(iii): given interval matrix H_iv (p x p x 2), interval vector rhs_iv (p x 2)
# where target is x = -H^{-1} rhs (or H^{-1} (-rhs)), stored candidate center_dbl (length p),
# and certified lower eigenvalue bound kappa_lo > 0 of M^{-1/2} H M^{-1/2}:
#   s = H * center + rhs,  rad = sup_{s in [s]} ||F^{-T} s||_2 / kappa_lo.
ipw_enclosure_lemma_dual_solve <- function(F_dbl, F_inv_T, coord_dual_iv, H_iv, rhs_iv, center_dbl, kappa_lo) {
  if (!is.finite(kappa_lo) || kappa_lo <= 0.0) {
    ipw_enclosure_arith_refuse(
      "Lemma dual(iii) lower spectral constant kappa is nonpositive.",
      reason = "lemma_dual_kappa_nonpositive",
      diagnostics = list(kappa_lo = kappa_lo)
    )
  }
  p <- nrow(F_dbl)
  s_iv <- matrix(0.0, nrow = p, ncol = 2L)
  ops <- 0
  for (k in seq_len(p)) {
    acc <- rhs_iv[k, ]
    for (l in seq_len(p)) {
      if (center_dbl[l] != 0.0) {
        acc <- gmm_iv_add(acc, gmm_iv_mul(H_iv[k, l, ], gmm_iv(center_dbl[l])))
        ops <- ops + 2
      }
    }
    s_iv[k, ] <- acc
  }
  dn <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, s_iv)
  ops <- ops + dn$ops
  rad_iv <- logit_b64_div(gmm_iv(dn$upper), gmm_iv(kappa_lo))
  ops <- ops + 1
  rad_hi <- rad_iv[2L]
  if (!is.finite(rad_hi) || rad_hi < 0.0) {
    ipw_enclosure_arith_refuse("Lemma dual(iii) radius overflowed.", reason = "lemma_dual_radius_overflow")
  }

  coord_iv <- matrix(0.0, nrow = p, ncol = 2L)
  for (k in seq_len(p)) {
    rk <- logit_b64_mul(gmm_iv(coord_dual_iv[k, 2L]), gmm_iv(rad_hi))[2L]
    ops <- ops + 1
    coord_iv[k, ] <- ipw_enclosure_widen_iv(gmm_iv(center_dbl[k]), rk)
    ops <- ops + 2
  }

  list(
    center = center_dbl,
    residual_s_iv = s_iv,
    whitened_residual_norm_upper = dn$upper,
    kappa_lower = kappa_lo,
    radius_m_norm = rad_hi,
    coord_iv = coord_iv,
    ops = ops,
    solves = dn$solves
  )
}

# Enclose dot product vec_iv' lambda where lambda in {center + d : ||d||_M <= rad}
ipw_enclosure_ball_dot <- function(F_dbl, F_inv_T, vec_iv, ball_sol, vec_dual_upper = NULL) {
  p <- nrow(F_dbl)
  is_pt <- is.null(dim(vec_iv))
  center_dot <- if (is_pt) {
    gmm_iv_dot(as.double(vec_iv), ball_sol$center)
  } else {
    acc <- c(0.0, 0.0)
    for (k in seq_len(p)) {
      if (ball_sol$center[k] != 0.0) {
        acc <- gmm_iv_add(acc, gmm_iv_mul(vec_iv[k, ], gmm_iv(ball_sol$center[k])))
      }
    }
    acc
  }
  ops <- 2 * p
  if (is.null(vec_dual_upper)) {
    dn <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, vec_iv)
    vec_dual_upper <- dn$upper
    ops <- ops + dn$ops
  }
  rad_dot <- logit_b64_mul(gmm_iv(vec_dual_upper), gmm_iv(ball_sol$radius_m_norm))[2L]
  ops <- ops + 1
  iv <- ipw_enclosure_widen_iv(center_dot, rad_dot)
  ops <- ops + 2
  list(iv = iv, dual_upper = vec_dual_upper, ops = ops)
}

# ---------------------------------------------------------------------------
# Component 1: Reusable, explicitly bounded certified nuisance-problem constructor
# ---------------------------------------------------------------------------

construct_certified_ipw_nuisance_problem <- function(
  id,
  Z,
  w,
  G,
  pair_observed,
  cohorts,
  beta0,
  F,
  r,
  certificate_control,
  problem_control,
  risk_mask = NULL,
  observed_mask = NULL,
  cohort_index = NULL
) {
  ipw_enclosure_require_deps()
  if (missing(id) || missing(Z) || missing(w) || missing(G) ||
      missing(pair_observed) || missing(cohorts) || missing(beta0) ||
      missing(F) || missing(r) || missing(certificate_control) ||
      missing(problem_control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`id`, `Z`, `w`, `G`, `pair_observed`, `cohorts`, `beta0`, `F`, `r`, `certificate_control`, `problem_control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }

  # 1. Validate controls and lightweight headers BEFORE any content scan or allocation
  cert_ctrl <- validate_logit_certificate_control(certificate_control)
  prob_ctrl <- validate_ipw_nuisance_problem_control(problem_control)

  # Refuse decorative/foreign certificate attributes on inputs
  for (obj_item in list(
    list(x = id, nm = "id"), list(x = Z, nm = "Z"), list(x = w, nm = "w"),
    list(x = G, nm = "G"), list(x = pair_observed, nm = "pair_observed"),
    list(x = cohorts, nm = "cohorts"), list(x = beta0, nm = "beta0"),
    list(x = F, nm = "F"), list(x = r, nm = "r")
  )) {
    ax <- attributes(obj_item$x)
    if (!is.null(ax)) {
      bad_attr <- intersect(names(ax), c("certified", "certificate", "converged", "is_pd", "status", "rho"))
      if (length(bad_attr) > 0L) {
        cdid_ipw_enclosure_abort(
          paste0("Decorative or foreign certificate attribute (`", bad_attr[1L], "`) on `", obj_item$nm, "` is refused."),
          c("cdid_ipw_enclosure_foreign_certificate_error", "cdid_ipw_enclosure_input_error")
        )
      }
    }
  }

  if (is.null(cohorts) || !is.atomic(cohorts) || !is.null(dim(cohorts)) ||
      is.object(cohorts) || !is.numeric(cohorts) || length(cohorts) < 2L) {
    cdid_ipw_enclosure_abort(
      "`cohorts` must be an explicit ordered numeric vector `c(target, controls)` of length >= 2.",
      c("cdid_ipw_enclosure_cohort_error", "cdid_ipw_enclosure_input_error")
    )
  }
  m <- length(cohorts)

  if (is.null(Z) || !is.matrix(Z) || is.object(Z) || !is.numeric(Z)) {
    cdid_ipw_enclosure_abort(
      "`Z` must be a plain numeric matrix (one row per eligible A-unit).",
      c("cdid_ipw_enclosure_covariate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  n <- nrow(Z)
  p <- ncol(Z)
  if (is.na(n) || is.na(p) || n < 1L || p < 1L) {
    cdid_ipw_enclosure_abort(
      "`Z` must have `nrow(Z) >= 1` and `ncol(Z) >= 1`.",
      c("cdid_ipw_enclosure_covariate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (is.null(id) || is.factor(id) || is.object(id) || !is.atomic(id) ||
      !is.null(dim(id)) || !(is.numeric(id) || is.character(id)) || length(id) != n) {
    cdid_ipw_enclosure_abort(
      "`id` must be a plain atomic numeric or character vector of length `nrow(Z)`.",
      c("cdid_ipw_enclosure_id_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (is.null(w) || !is.atomic(w) || !is.null(dim(w)) || is.object(w) ||
      !is.numeric(w) || length(w) != n) {
    cdid_ipw_enclosure_abort(
      "`w` must be an atomic numeric vector of length `nrow(Z)`.",
      c("cdid_ipw_enclosure_weight_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (is.null(G) || !is.atomic(G) || !is.null(dim(G)) || is.object(G) ||
      !is.numeric(G) || length(G) != n) {
    cdid_ipw_enclosure_abort(
      "`G` must be an atomic numeric vector of length `nrow(Z)`.",
      c("cdid_ipw_enclosure_cohort_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (is.null(pair_observed) || !is.atomic(pair_observed) || !is.null(dim(pair_observed)) ||
      is.object(pair_observed) || !is.logical(pair_observed) || length(pair_observed) != n) {
    cdid_ipw_enclosure_abort(
      "`pair_observed` must be an atomic logical vector of length `nrow(Z)`.",
      c("cdid_ipw_enclosure_pair_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (is.null(beta0) || !is.atomic(beta0) || !is.null(dim(beta0)) || is.object(beta0) ||
      !is.numeric(beta0) || length(beta0) != p) {
    cdid_ipw_enclosure_abort(
      "`beta0` must be a finite numeric vector of length `ncol(Z)`.",
      c("cdid_ipw_enclosure_beta_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.matrix(F) || is.object(F) || !is.numeric(F) || nrow(F) != p || ncol(F) != p) {
    cdid_ipw_enclosure_abort(
      "`F` must be an explicit numeric p x p matrix.",
      c("cdid_ipw_enclosure_metric_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!gmm_scalar_numeric(r, positive = TRUE)) {
    cdid_ipw_enclosure_abort(
      "`r` must be an explicit finite positive numeric scalar radius.",
      c("cdid_ipw_enclosure_radius_error", "cdid_ipw_enclosure_input_error")
    )
  }
  r_dbl <- as.double(r)

  # 2. Preflight BEFORE content scans or copy allocations
  preflight <- preflight_ipw_nuisance_problem(
    n = n, p = p, m = m,
    certificate_control = cert_ctrl,
    problem_control = prob_ctrl
  )
  if (m > n) {
    cdid_ipw_enclosure_abort(
      "Declared cohorts cannot all have positive mass when their count exceeds n_A.",
      c("cdid_ipw_enclosure_cohort_support_error", "cdid_ipw_enclosure_support_error")
    )
  }

  coh_dbl <- validate_logit_certificate_cohorts(cohorts)
  cohort_labels <- c("g", paste0("k_", seq_len(m - 1L)))
  if (anyNA(beta0) || any(!is.finite(beta0))) {
    cdid_ipw_enclosure_abort(
      "`beta0` must be finite.",
      c("cdid_ipw_enclosure_beta_error", "cdid_ipw_enclosure_input_error")
    )
  }
  beta0_dbl <- as.double(beta0)
  F_dbl <- validate_logit_certificate_metric_factor(F, p)
  step <- preflight$effective_chunk_rows

  # 3. Validate IDs and create owned copies of raw inputs
  id_info <- ipw_enclosure_validate_id(id, n = n, step = step, label = "id")

  # Call certify_chained_logit_root on the raw inputs (validates Z, w, G, pair_observed in chunks)
  Z_copy <- matrix(as.double(Z), nrow = n, ncol = p, dimnames = list(NULL, colnames(Z)))
  w_copy <- as.double(w)
  G_copy <- as.double(G)
  po_copy <- as.logical(pair_observed)

  cert <- certify_chained_logit_root(
    Z = Z_copy,
    w = w_copy,
    G = G_copy,
    pair_observed = po_copy,
    cohorts = coh_dbl,
    beta0 = beta0_dbl,
    F = F_dbl,
    r = r_dbl,
    control = cert_ctrl
  )

  # Derive canonical masks R = {G in cohorts}, O = R intersect S
  cohort_idx <- as.integer(match(G_copy, coh_dbl))
  in_risk <- !is.na(cohort_idx)
  in_observed <- in_risk & po_copy
  risk_indices <- which(in_risk)
  observed_indices <- which(in_observed)

  ipw_enclosure_check_optional_masks(
    n = n,
    in_risk = in_risk,
    in_observed = in_observed,
    cohort_index = cohort_idx,
    pair_observed_ref = po_copy,
    risk_mask = risk_mask,
    observed_mask = observed_mask,
    cohort_index_arg = cohort_index
  )

  # 4. Run shared outcome-free nuisance transfer worker over sealed dense context
  ctx <- construct_dense_ipw_nuisance_access_context(
    Z_copy = Z_copy,
    w_copy = w_copy,
    po_copy = po_copy,
    cohort_idx = cohort_idx,
    coh_dbl = coh_dbl,
    cohort_labels = cohort_labels,
    beta0_dbl = beta0_dbl,
    F_dbl = F_dbl,
    r_dbl = r_dbl,
    cert_ctrl = cert_ctrl,
    prob_ctrl = prob_ctrl,
    preflight = preflight
  )
  tr_res <- ipw_enclosure_nuisance_transfer_with_context(ctx, cert)
  nuisance_enclosures <- tr_res$nuisance_enclosures
  arith_refusal <- tr_res$arith_refusal
  extra_reasons <- tr_res$extra_reasons
  added_ops <- tr_res$added_ops
  added_chunks <- tr_res$added_chunks
  added_passes <- tr_res$added_passes
  added_solves <- tr_res$added_solves

  stopifnot(added_ops <= preflight$added_work_ops_bound,
    cert$work$actual_interval_ops + added_ops <= preflight$total_work_ops_bound)

  certified <- isTRUE(cert$certified) && !is.null(nuisance_enclosures)
  inconclusive_reasons <- unique(c(cert$inconclusive_reasons, extra_reasons))
  inconclusive_type <- if (certified) {
    "none"
  } else if (!is.null(arith_refusal) || identical(cert$inconclusive_type, "arithmetic")) {
    "arithmetic"
  } else {
    cert$inconclusive_type
  }

  cert_sealed <- ipw_enclosure_seal(cert, "certificate")

  problem_obj <- structure(
    list(
      certified = certified,
      status = if (certified) "CERTIFIED" else "INCONCLUSIVE",
      inconclusive = !certified,
      inconclusive_type = inconclusive_type,
      inconclusive_reasons = inconclusive_reasons,
      separation_certified = FALSE,
      sampling_certified = FALSE,
      public_gate_adopted = FALSE,
      raw = list(
        id = id_info$id_copy,
        id_key = id_info$id_key,
        id_type = id_info$id_type,
        id_storage = id_info$id_storage,
        Z = Z_copy,
        w = w_copy,
        G = G_copy,
        pair_observed = po_copy,
        cohorts = coh_dbl,
        cohort_labels = cohort_labels,
        beta0 = beta0_dbl,
        F = F_dbl,
        r = r_dbl
      ),
      masks = list(
        cohort_index = cohort_idx,
        in_risk = in_risk,
        in_observed = in_observed,
        risk_indices = as.integer(risk_indices),
        observed_indices = as.integer(observed_indices)
      ),
      dimensions = cert$dimensions,
      certificate = cert,
      certificate_sealed = cert_sealed,
      nuisance_enclosures = nuisance_enclosures,
      arithmetic_refusal = if (!is.null(arith_refusal)) {
        list(
          occurred = TRUE,
          reason = if (!is.null(arith_refusal$reason)) arith_refusal$reason else "arithmetic_refusal",
          message = conditionMessage(arith_refusal),
          condition_class = class(arith_refusal),
          diagnostics = arith_refusal$diagnostics
        )
      } else {
        cert$arithmetic_refusal
      },
      controls = list(
        certificate_control = cert_ctrl,
        problem_control = prob_ctrl
      ),
      preflight = preflight,
      work = list(
        raw_certificate_work = cert$work,
        added_passes_completed = added_passes,
        added_chunks_processed = logit_cert_count_output(added_chunks),
        added_triangular_solves = logit_cert_count_output(added_solves),
        added_interval_ops = added_ops,
        total_interval_ops = cert$work$actual_interval_ops + added_ops,
        interval_ops_definition = "Conservative interval-operation charges at completed helper boundaries; excludes CPU instructions, R object headers, and process peak RSS."
      ),
      provenance = list(
        normalizer_n = as.double(n),
        normalizer_definition = "count of all eligible A-units (including outside-risk and pair non-observers)",
        cohorts = coh_dbl,
        cohort_labels = cohort_labels,
        target_cohort = coh_dbl[1L],
        control_cohorts = coh_dbl[-1L],
        immutability_contract = "Owned copy-on-modify R arrays sealed in a locked environment; any mutation of problem fields or certificate is refused. Ordinary R in-memory integrity check, not a portable cryptographic identity or sampling-eligibility proof."
      )
    ),
    class = "cdid_ipw_nuisance_problem"
  )

  ipw_enclosure_seal(problem_obj, "nuisance_problem")
}

construct_dense_ipw_nuisance_access_context <- function(
  Z_copy,
  w_copy,
  po_copy,
  cohort_idx,
  coh_dbl,
  cohort_labels,
  beta0_dbl,
  F_dbl,
  r_dbl,
  cert_ctrl,
  prob_ctrl,
  preflight
) {
  n <- nrow(Z_copy)
  p <- ncol(Z_copy)
  m <- length(coh_dbl)
  header <- list(
    mode = "dense",
    n = as.integer(n),
    p = as.integer(p),
    m = as.integer(m),
    z_colnames = colnames(Z_copy),
    coh_dbl = coh_dbl,
    cohort_labels = cohort_labels,
    beta0_dbl = beta0_dbl,
    F_dbl = F_dbl,
    r_dbl = r_dbl,
    cert_ctrl = cert_ctrl,
    prob_ctrl = prob_ctrl,
    reader_ctrl = NULL,
    preflight = preflight,
    step = preflight$effective_chunk_rows
  )
  read_chunk <- function(a_start, a_end) {
    idx <- seq.int(a_start, a_end)
    list(
      idx = idx,
      b_len = length(idx),
      Z_ch = Z_copy[idx, , drop = FALSE],
      w_ch = w_copy[idx],
      po_ch = po_copy[idx],
      j_ch = cohort_idx[idx],
      reader_ops = 0,
      n_endpoints_read = 0L
    )
  }
  ipw_seal_access_context(
    list(mode = "dense", header = header, read_chunk = read_chunk),
    mode = "dense",
    phase = "nuisance_transfer"
  )
}

construct_view_ipw_nuisance_access_context <- function(
  nuisance_view,
  coh_dbl,
  cohort_labels,
  beta0_dbl,
  F_dbl,
  r_dbl,
  cert_ctrl,
  prob_ctrl,
  reader_ctrl,
  preflight,
  expected_a_id = NULL,
  trace_env = NULL
) {
  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  n <- as.integer(frame_owner$n_a_units)
  p <- as.integer(frame_owner$n_covariates)
  m <- length(coh_dbl)
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
    cert_ctrl = cert_ctrl,
    prob_ctrl = prob_ctrl,
    reader_ctrl = reader_ctrl,
    preflight = preflight,
    step = preflight$effective_chunk_rows
  )
  read_chunk <- function(a_start, a_end) {
    exp_slice <- if (is.null(expected_a_id)) NULL else expected_a_id[seq.int(a_start, a_end)]
    ch_tr <- if (!is.null(trace_env)) new.env(parent = emptyenv()) else NULL
    ch <- query_ipw_view_chunk(
      view = nuisance_view,
      a_start = a_start,
      a_end = a_end,
      control = reader_ctrl,
      read_outcomes = FALSE,
      include_design = TRUE,
      enclose_differences = FALSE,
      expected_a_id = exp_slice,
      trace_env = ch_tr
    )
    if (!is.null(trace_env) && isTRUE(ch_tr$calls > 0L)) {
      trace_env$calls <- trace_env$calls + ch_tr$calls
      trace_env$pos_s_index <- c(trace_env$pos_s_index, ch_tr$pos_s_index)
      trace_env$pos_t_index <- c(trace_env$pos_t_index, ch_tr$pos_t_index)
    }
    list(
      idx = ch$a_index,
      b_len = ch$chunk_rows,
      Z_ch = ch$design$Z_chunk,
      w_ch = ch$weight,
      po_ch = ch$in_observed,
      j_ch = ch$cohort_index,
      reader_ops = ch$work$actual_work_ops,
      n_endpoints_read = ch$read_diagnostics$n_endpoints_read
    )
  }
  ipw_seal_access_context(
    list(mode = "compact_view", header = header, read_chunk = read_chunk),
    mode = "compact_view",
    phase = "nuisance_transfer"
  )
}

ipw_enclosure_nuisance_transfer_with_context <- function(ctx, cert) {
  ipw_validate_access_context(ctx, expected_phase = "nuisance_transfer")
  hdr <- ctx$header
  n <- hdr$n
  p <- hdr$p
  m <- hdr$m
  z_colnames <- hdr$z_colnames
  cohort_labels <- hdr$cohort_labels
  beta0_dbl <- hdr$beta0_dbl
  F_dbl <- hdr$F_dbl
  cert_ctrl <- hdr$cert_ctrl
  step <- hdr$step

  added_ops <- 0
  added_reader_ops <- 0
  added_endpoints_read <- 0L
  added_chunks <- 0
  added_passes <- 0L
  added_solves <- 0L
  n_iv <- gmm_iv(as.double(n))

  nuisance_enclosures <- NULL
  arith_refusal <- NULL
  extra_reasons <- character(0L)

  if (isTRUE(cert$certified)) {
    nuis_try <- tryCatch({
      M_cohort <- cert$enclosures$M_cohort
      Pi_cohort <- cert$enclosures$Pi_cohort
      r_cohort <- cert$enclosures$r_cohort
      q_cohort <- cert$enclosures$q_cohort
      M_over_Pi_cohort <- cert$enclosures$M_over_Pi_cohort
      H_0_iv <- cert$enclosures$H

      delta_hi <- cert$delta
      a_hi <- cert$a
      a_c_hi <- cert$a_control
      rho_hi <- cert$rho

      x_bar_iv <- logit_b64_mul(gmm_iv(a_hi), gmm_iv(rho_hi))
      x_bar_hi <- x_bar_iv[2L]
      x_bar_c_iv <- logit_b64_mul(gmm_iv(a_c_hi), gmm_iv(rho_hi))
      x_bar_c_hi <- min(x_bar_hi, x_bar_c_iv[2L])
      added_ops <- added_ops + 2

      if (!is.finite(x_bar_hi) || x_bar_hi >= 1.0) {
        ipw_enclosure_arith_refuse(
          "Root-ball covariate factor a * rho >= 1 in Theorem A transfer.",
          reason = "x_bar_ge_1",
          diagnostics = list(a = a_hi, rho = rho_hi, x_bar = x_bar_hi)
        )
      }

      one_minus_x_bar_iv <- logit_b64_sub(gmm_iv(1.0), gmm_iv(x_bar_hi))
      one_minus_x_bar_c_iv <- logit_b64_sub(gmm_iv(1.0), gmm_iv(x_bar_c_hi))
      one_minus_delta_iv <- logit_b64_sub(gmm_iv(1.0), gmm_iv(delta_hi))
      one_plus_delta_iv <- logit_b64_add(gmm_iv(1.0), gmm_iv(delta_hi))
      added_ops <- added_ops + 4

      if (one_minus_x_bar_iv[1L] <= 0.0 || one_minus_x_bar_c_iv[1L] <= 0.0 || one_minus_delta_iv[1L] <= 0.0) {
        ipw_enclosure_arith_refuse(
          "Denominator (1 - a*rho) or (1 - delta) nonpositive.",
          reason = "nonpositive_transfer_denominator",
          diagnostics = list(
            one_minus_x_bar_lo = one_minus_x_bar_iv[1L],
            one_minus_delta_lo = one_minus_delta_iv[1L]
          )
        )
      }

      kappa_0_lo <- one_minus_delta_iv[1L]
      kappa_star_iv <- logit_b64_mul(gmm_iv(one_minus_x_bar_iv[1L]), gmm_iv(one_minus_delta_iv[1L]))
      kappa_star_lo <- kappa_star_iv[1L]
      added_ops <- added_ops + 1
      if (!is.finite(kappa_star_lo) || kappa_star_lo <= 0.0) {
        ipw_enclosure_arith_refuse(
          "Lemma dual(iii) lower constant (1 - a*rho)(1 - delta) <= 0.",
          reason = "kappa_star_nonpositive",
          diagnostics = list(kappa_star_lo = kappa_star_lo)
        )
      }

      # u1(x_bar) = x_bar / (1 - x_bar) and u2(x_bar_c) = 1/(1 - x_bar_c)^2 - 1 = x_bar_c*(2 - x_bar_c)/(1 - x_bar_c)^2
      u1_xbar_hi <- logit_b64_div(gmm_iv(x_bar_hi), gmm_iv(one_minus_x_bar_iv[1L]))[2L]
      one_m_xbc_sq_lo <- logit_b64_mul(gmm_iv(one_minus_x_bar_c_iv[1L]), gmm_iv(one_minus_x_bar_c_iv[1L]))[1L]
      if (!is.finite(one_m_xbc_sq_lo) || one_m_xbc_sq_lo <= 0.0) {
        ipw_enclosure_arith_refuse("Denominator (1 - a_c*rho)^2 nonpositive.", reason = "one_minus_xbc_sq_nonpositive")
      }
      u2_num_iv <- logit_b64_mul(gmm_iv(x_bar_c_hi), logit_b64_sub(gmm_iv(2.0), gmm_iv(x_bar_c_hi)))
      u2_xbar_c_hi <- logit_b64_div(u2_num_iv, gmm_iv(one_m_xbc_sq_lo))[2L]
      added_ops <- added_ops + 5

      # Enclose F^{-T} and coordinate dual norms ||F^{-T} e_k||_2 = sqrt((M^{-1})_{kk})
      ft_basis <- ipw_enclosure_ft_inv_basis(F_dbl)
      F_inv_T <- ft_basis$F_inv_T
      coord_dual_iv <- ft_basis$coord_dual_iv
      added_ops <- added_ops + ft_basis$ops
      added_solves <- added_solves + ft_basis$solves

      # Coordinate projections beta*_k in beta_{0k} +/- rho * ||F^{-T} e_k||_2 (Lemma dual(ii))
      beta_star_box <- matrix(0.0, nrow = p, ncol = 2L, dimnames = list(z_colnames, c("lo", "hi")))
      for (k in seq_len(p)) {
        rad_bk <- logit_b64_mul(gmm_iv(coord_dual_iv[k, 2L]), gmm_iv(rho_hi))[2L]
        beta_star_box[k, ] <- ipw_enclosure_widen_iv(gmm_iv(beta0_dbl[k]), rad_bk)
        added_ops <- added_ops + 3
      }

      # Streamed pass over chunks to enclose Sigma_plus, d^+, Sigma_minus_0, d^{-,0},
      # Sigma_minus_star, d^{-,*}, H^*, L^0, L^* (all cohorts including target), and b_j^0 sums
      Sigma_plus_iv <- c(0.0, 0.0)
      Sigma_minus_0_iv <- c(0.0, 0.0)
      Sigma_minus_star_iv <- c(0.0, 0.0)
      b_cohort_0_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
      H_star_sum_iv <- array(0.0, dim = c(p, p, 2L))
      L_0_sum_iv <- array(0.0, dim = c(p, m, 2L))
      L_star_sum_iv <- array(0.0, dim = c(p, m, 2L))

      a_start <- 1L
      while (a_start <= n) {
        ipw_enclosure_check_interrupt()
        a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
        ch <- ctx$read_chunk(a_start, a_end)
        added_reader_ops <- added_reader_ops + ch$reader_ops
        added_endpoints_read <- added_endpoints_read + ch$n_endpoints_read
        idx <- ch$idx
        b_len <- ch$b_len
        Z_ch <- ch$Z_ch
        w_ch <- ch$w_ch
        po_ch <- ch$po_ch
        j_ch <- ch$j_ch

        for (u in seq_len(b_len)) {
          ju <- j_ch[u]
          if (is.na(ju) || !po_ch[u]) next

          wu_iv <- gmm_iv(w_ch[u])
          z_u <- Z_ch[u, ]
          D_u <- if (ju == 1L) 1L else 0L

          cw <- ipw_enclosure_c_weight(
            wu_iv, M_cohort[ju, ], Pi_cohort[ju, ], M_over_Pi_cohort[ju, ],
            unit_offset = idx[u], cohort_index = ju
          )
          c_iv <- cw$c_iv
          added_ops <- added_ops + cw$ops

          # Dual norm a_i = min(a, ||F^{-T} z_i||_2) ( clipped by a_c for controls, Lemma mult )
          if (all(z_u == 0.0)) {
            ai_hi <- 0.0
          } else {
            dn_z <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, z_u)
            added_ops <- added_ops + dn_z$ops
            added_solves <- added_solves + dn_z$solves
            ai_hi <- min(dn_z$upper, if (D_u == 0L) min(a_c_hi, a_hi) else a_hi)
          }
          xi_hi <- min(x_bar_hi, logit_b64_mul(gmm_iv(ai_hi), gmm_iv(rho_hi))[2L])
          one_m_xi_lo <- logit_b64_sub(gmm_iv(1.0), gmm_iv(xi_hi))[1L]
          added_ops <- added_ops + 2
          if (!is.finite(one_m_xi_lo) || one_m_xi_lo <= 0.0) {
            ipw_enclosure_arith_refuse(
              "Per-row denominator 1 - x_i nonpositive.",
              reason = "one_minus_xi_nonpositive",
              diagnostics = list(unit_offset = idx[u], xi_hi = xi_hi)
            )
          }

          # Logistic primitives at beta0 (Route B)
          eta_iv <- gmm_iv_dot(z_u, beta0_dbl)
          added_ops <- added_ops + 2 * p
          exp_res <- logit_iv_exp(
            eta_iv,
            taylor_degree = cert_ctrl$taylor_degree,
            max_abs_eta = cert_ctrl$max_abs_eta,
            max_squarings = cert_ctrl$max_squarings
          )
          o0_iv <- exp_res$iv
          added_ops <- added_ops + exp_res$ops
          logis0 <- logit_iv_logistic_from_odds(o0_iv)
          p0_iv <- logis0$prob
          omp0_iv <- logis0$one_minus_prob
          added_ops <- added_ops + logis0$ops

          # Theorem A step 1: [o_i^*] = [o_i^0_lo * (1 - x_i), o_i^0_hi / (1 - x_i)]
          o_star_lo <- logit_b64_mul(gmm_iv(o0_iv[1L]), gmm_iv(one_m_xi_lo))[1L]
          o_star_hi <- logit_b64_div(gmm_iv(o0_iv[2L]), gmm_iv(one_m_xi_lo))[2L]
          added_ops <- added_ops + 2
          if (!is.finite(o_star_lo) || !is.finite(o_star_hi) || o_star_lo <= 0.0 || o_star_lo > o_star_hi) {
            ipw_enclosure_arith_refuse(
              "Root-ball transferred odds interval [o_i^*] failed positivity.",
              reason = "odds_star_nonpositive",
              diagnostics = list(unit_offset = idx[u], o_star_lo = o_star_lo, o_star_hi = o_star_hi)
            )
          }
          o_star_iv <- c(o_star_lo, o_star_hi)
          logis_star <- logit_iv_logistic_from_odds(o_star_iv)
          p_star_iv <- logis_star$prob
          omp_star_iv <- logis_star$one_minus_prob
          v_star_iv <- logis_star$var
          added_ops <- added_ops + logis_star$ops

          if (D_u == 1L) {
            Sigma_plus_iv <- gmm_iv_add(Sigma_plus_iv, c_iv)
            added_ops <- added_ops + 1
          } else {
            b0_u <- gmm_iv_mul(c_iv, o0_iv)
            b_star_u <- gmm_iv_mul(c_iv, o_star_iv)
            Sigma_minus_0_iv <- gmm_iv_add(Sigma_minus_0_iv, b0_u)
            Sigma_minus_star_iv <- gmm_iv_add(Sigma_minus_star_iv, b_star_u)
            b_cohort_0_iv[ju, ] <- gmm_iv_add(b_cohort_0_iv[ju, ], b0_u)
            added_ops <- added_ops + 5
          }

          # Accumulate L_j^0 and L_j^*: sum_{O, j(i)=j} c_i * (D_i - p_i) * z_i
          dmp0_iv <- if (D_u == 1L) omp0_iv else gmm_iv_neg(p0_iv)
          dmp_star_iv <- if (D_u == 1L) omp_star_iv else gmm_iv_neg(p_star_iv)
          c_dmp0_iv <- gmm_iv_mul(c_iv, dmp0_iv)
          c_dmp_star_iv <- gmm_iv_mul(c_iv, dmp_star_iv)
          added_ops <- added_ops + 2

          for (k in seq_len(p)) {
            zk_iv <- gmm_iv(z_u[k])
            L_0_sum_iv[k, ju, ] <- gmm_iv_add(L_0_sum_iv[k, ju, ], gmm_iv_mul(c_dmp0_iv, zk_iv))
            L_star_sum_iv[k, ju, ] <- gmm_iv_add(L_star_sum_iv[k, ju, ], gmm_iv_mul(c_dmp_star_iv, zk_iv))
            added_ops <- added_ops + 4
          }

          # Accumulate H^*: sum_O c_i * v_i^* * z_i z_i'
          c_v_star_iv <- gmm_iv_mul(c_iv, v_star_iv)
          added_ops <- added_ops + 1
          for (k in seq_len(p)) {
            cv_zk <- gmm_iv_mul(c_v_star_iv, gmm_iv(z_u[k]))
            added_ops <- added_ops + 1
            for (l in seq.int(k, p)) {
              term_h <- gmm_iv_mul(cv_zk, gmm_iv(z_u[l]))
              if (k == l) term_h[1L] <- max(0.0, term_h[1L])
              H_star_sum_iv[k, l, ] <- gmm_iv_add(H_star_sum_iv[k, l, ], term_h)
              added_ops <- added_ops + 2
            }
          }
        }
        added_chunks <- added_chunks + 1
        a_start <- as.double(a_end) + 1
      }
      added_passes <- 1L

      if (Sigma_plus_iv[1L] <= 0.0 || Sigma_minus_0_iv[1L] <= 0.0 || Sigma_minus_star_iv[1L] <= 0.0) {
        ipw_enclosure_arith_refuse(
          "Treated or control normalization denominator Sigma_+ or Sigma_- underflowed to <= 0.",
          reason = "denominator_sigma_nonpositive",
          diagnostics = list(
            Sigma_plus = Sigma_plus_iv,
            Sigma_minus_0 = Sigma_minus_0_iv,
            Sigma_minus_star = Sigma_minus_star_iv
          )
        )
      }

      d_plus_iv <- gmm_iv_div(Sigma_plus_iv, n_iv)
      d_minus_0_iv <- gmm_iv_div(Sigma_minus_0_iv, n_iv)
      d_minus_star_iv <- gmm_iv_div(Sigma_minus_star_iv, n_iv)
      added_ops <- added_ops + 3
      if (d_plus_iv[1L] <= 0.0 || d_minus_0_iv[1L] <= 0.0 || d_minus_star_iv[1L] <= 0.0) {
        ipw_enclosure_arith_refuse(
          "Normalized denominator d^+ or d^- underflowed to <= 0.",
          reason = "denominator_d_nonpositive",
          diagnostics = list(d_plus = d_plus_iv, d_minus_0 = d_minus_0_iv, d_minus_star = d_minus_star_iv)
        )
      }

      H_star_iv <- array(0.0, dim = c(p, p, 2L), dimnames = list(z_colnames, z_colnames, c("lo", "hi")))
      for (k in seq_len(p)) {
        for (l in seq.int(k, p)) {
          h_kl <- gmm_iv_div(H_star_sum_iv[k, l, ], n_iv)
          if (k == l) h_kl[1L] <- max(0.0, h_kl[1L])
          H_star_iv[k, l, ] <- h_kl
          H_star_iv[l, k, ] <- h_kl
          added_ops <- added_ops + 1
        }
      }

      L_0_iv <- array(0.0, dim = c(p, m, 2L), dimnames = list(z_colnames, cohort_labels, c("lo", "hi")))
      L_star_iv <- array(0.0, dim = c(p, m, 2L), dimnames = list(z_colnames, cohort_labels, c("lo", "hi")))
      ell_0_upper <- numeric(m)
      ell_star_upper <- numeric(m)
      e_L_upper <- numeric(m)
      names(ell_0_upper) <- names(ell_star_upper) <- names(e_L_upper) <- cohort_labels

      for (j in seq_len(m)) {
        nr_j_iv <- gmm_iv_mul(n_iv, r_cohort[j, ])
        added_ops <- added_ops + 1
        for (k in seq_len(p)) {
          L_0_iv[k, j, ] <- gmm_iv_neg(gmm_iv_div(L_0_sum_iv[k, j, ], nr_j_iv))
          L_star_iv[k, j, ] <- gmm_iv_neg(gmm_iv_div(L_star_sum_iv[k, j, ], nr_j_iv))
          added_ops <- added_ops + 2
        }
        dn_L0 <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, matrix(L_0_iv[, j, ], nrow = p, ncol = 2L))
        dn_Ls <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, matrix(L_star_iv[, j, ], nrow = p, ncol = 2L))
        ell_0_upper[j] <- dn_L0$upper
        ell_star_upper[j] <- dn_Ls$upper
        added_ops <- added_ops + dn_L0$ops + dn_Ls$ops
        added_solves <- added_solves + dn_L0$solves + dn_Ls$solves

        # Theorem B: e_{L,j} = (1 + delta) * rho / ((1 - x_bar) * r_j_lo)
        num_eL <- logit_b64_mul(gmm_iv(one_plus_delta_iv[2L]), gmm_iv(rho_hi))
        den_eL <- logit_b64_mul(gmm_iv(one_minus_x_bar_iv[1L]), gmm_iv(r_cohort[j, 1L]))
        e_L_upper[j] <- logit_b64_div(num_eL, den_eL)[2L]
        added_ops <- added_ops + 3
      }

      # Control cohort weight share upper bounds pi_{(j)}^0 and bar{pi}_j for Theorem B
      pi_cohort_0_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
      pi_bar_cohort <- numeric(m)
      names(pi_bar_cohort) <- cohort_labels
      for (j in seq.int(2L, m)) {
        pij_iv <- gmm_iv_div(b_cohort_0_iv[j, ], Sigma_minus_0_iv)
        pij_iv[1L] <- max(0.0, pij_iv[1L])
        pij_iv[2L] <- min(1.0, pij_iv[2L])
        pi_cohort_0_iv[j, ] <- pij_iv
        pbar_hi <- logit_b64_div(gmm_iv(pij_iv[2L]), gmm_iv(one_m_xbc_sq_lo))[2L]
        pi_bar_cohort[j] <- min(1.0, pbar_hi)
        added_ops <- added_ops + 2
      }

      list(
        Sigma_plus = Sigma_plus_iv,
        d_plus = d_plus_iv,
        Sigma_minus_0 = Sigma_minus_0_iv,
        d_minus_0 = d_minus_0_iv,
        Sigma_minus_star = Sigma_minus_star_iv,
        d_minus_star = d_minus_star_iv,
        b_cohort_0 = b_cohort_0_iv,
        pi_cohort_0 = pi_cohort_0_iv,
        pi_bar_cohort = pi_bar_cohort,
        H_0 = H_0_iv,
        H_star = H_star_iv,
        L_0 = L_0_iv,
        L_star = L_star_iv,
        ell_0_upper = ell_0_upper,
        ell_star_upper = ell_star_upper,
        e_L_upper = e_L_upper,
        F_inv_T = F_inv_T,
        coord_dual_iv = coord_dual_iv,
        beta_star_box = beta_star_box,
        x_bar = x_bar_hi,
        x_bar_c = x_bar_c_hi,
        u1_x_bar = u1_xbar_hi,
        u2_x_bar_c = u2_xbar_c_hi,
        one_minus_x_bar_lo = one_minus_x_bar_iv[1L],
        one_minus_x_bar_c_lo = one_minus_x_bar_c_iv[1L],
        one_minus_x_bar_c_sq_lo = one_m_xbc_sq_lo,
        one_minus_delta_lo = one_minus_delta_iv[1L],
        one_plus_delta_hi = one_plus_delta_iv[2L],
        kappa_0_lower = kappa_0_lo,
        kappa_star_lower = kappa_star_lo
      )
    }, cdid_ipw_enclosure_arithmetic_inconclusive = function(e) e,
       cdid_logit_arithmetic_inconclusive = function(e) e,
       cdid_gmm_accuracy_error = function(e) e)

    if (inherits(nuis_try, "error")) {
      arith_refusal <- nuis_try
      reason <- if (!is.null(nuis_try$reason) && nzchar(nuis_try$reason)) nuis_try$reason else "nuisance_transfer_arithmetic_refusal"
      extra_reasons <- c(extra_reasons, reason)
    } else {
      nuisance_enclosures <- nuis_try
    }
  }

  list(
    nuisance_enclosures = nuisance_enclosures,
    arith_refusal = arith_refusal,
    extra_reasons = extra_reasons,
    added_ops = added_ops,
    added_reader_ops = added_reader_ops,
    added_endpoints_read = added_endpoints_read,
    added_chunks = added_chunks,
    added_passes = added_passes,
    added_solves = added_solves
  )
}

construct_certified_ipw_nuisance_problem_view <- function(
  nuisance_view,
  beta0,
  F,
  r,
  certificate_control,
  problem_control,
  reader_control
) {
  ipw_enclosure_require_deps()
  if (missing(nuisance_view) || missing(beta0) || missing(F) || missing(r) ||
      missing(certificate_control) || missing(problem_control) || missing(reader_control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`nuisance_view`, `beta0`, `F`, `r`, `certificate_control`, `problem_control`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }

  cert_ctrl <- validate_logit_certificate_control(certificate_control)
  prob_ctrl <- validate_ipw_nuisance_problem_control(problem_control)
  expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  trace_env <- attr(reader_control, "trace_env", exact = TRUE)
  r_ctrl <- validate_ipw_owner_control(reader_control)
  validate_ipw_nuisance_view(nuisance_view)

  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = nuisance_view$panel_ref, inherits = FALSE)
  n <- as.integer(frame_owner$n_a_units)
  p <- as.integer(frame_owner$n_covariates)
  z_colnames <- colnames(frame_owner$Z_a)
  cohorts <- nuisance_view$cohorts$cohorts
  m <- length(cohorts)

  if (is.null(beta0) || !is.atomic(beta0) || !is.null(dim(beta0)) || is.object(beta0) ||
      !is.numeric(beta0) || length(beta0) != p) {
    cdid_ipw_enclosure_abort(
      "`beta0` must be a finite numeric vector of length `ncol(Z)`.",
      c("cdid_ipw_enclosure_beta_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.null(attr(beta0, "converged")) || !is.null(attr(beta0, "certified"))) {
    cdid_ipw_enclosure_abort(
      "Decorative convergence/certification attributes on `beta0` are refused.",
      c("cdid_ipw_enclosure_beta_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.matrix(F) || is.object(F) || !is.numeric(F) || nrow(F) != p || ncol(F) != p) {
    cdid_ipw_enclosure_abort(
      "`F` must be an explicit numeric p x p matrix.",
      c("cdid_ipw_enclosure_metric_error", "cdid_ipw_enclosure_input_error")
    )
  }
  # Explicit typed automatic-radius request (compact view only) or an explicit
  # positive radius; every other value keeps the original typed refusal.
  r_auto <- logit_cert_is_automatic_radius_request(r)
  if (!r_auto && !gmm_scalar_numeric(r, positive = TRUE)) {
    cdid_ipw_enclosure_abort(
      "`r` must be an explicit finite positive numeric scalar radius.",
      c("cdid_ipw_enclosure_radius_error", "cdid_ipw_enclosure_input_error")
    )
  }
  r_dbl <- if (r_auto) NA_real_ else as.double(r)

  # Preflight combined 4-pass reader + arithmetic work and simultaneous concurrent payload BEFORE any chunk reads
  preflight <- preflight_ipw_nuisance_problem_view(
    nuisance_view = nuisance_view,
    certificate_control = cert_ctrl,
    problem_control = prob_ctrl,
    reader_control = r_ctrl,
    expected_a_id = expected_a_id,
    automatic_radius = r_auto
  )
  if (m > n) {
    cdid_ipw_enclosure_abort(
      "Declared cohorts cannot all have positive mass when their count exceeds n_A.",
      c("cdid_ipw_enclosure_cohort_support_error", "cdid_ipw_enclosure_support_error")
    )
  }

  coh_dbl <- validate_logit_certificate_cohorts(cohorts)
  cohort_labels <- c("g", paste0("k_", seq_len(m - 1L)))
  if (anyNA(beta0) || any(!is.finite(beta0))) {
    cdid_ipw_enclosure_abort(
      "`beta0` must be finite.",
      c("cdid_ipw_enclosure_beta_error", "cdid_ipw_enclosure_input_error")
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
        cdid_ipw_enclosure_abort(
          "Explicit `beta0` disagrees with optional candidate stored on `nuisance_view`.",
          c("cdid_ipw_enclosure_candidate_mismatch_error", "cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_beta_error", "cdid_ipw_enclosure_input_error")
        )
      }
    }
    if (!is.null(vc$F)) {
      cand_F <- matrix(as.double(vc$F), nrow = p, ncol = p)
      cand_F[cand_F == 0.0] <- 0.0
      if (!identical(F_dbl, cand_F)) {
        cdid_ipw_enclosure_abort(
          "Explicit `F` disagrees with optional candidate stored on `nuisance_view`.",
          c("cdid_ipw_enclosure_candidate_mismatch_error", "cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_metric_error", "cdid_ipw_enclosure_input_error")
        )
      }
    }
    if (!is.null(vc$r)) {
      if (r_auto) {
        cdid_ipw_enclosure_abort(
          "An automatic-radius request conflicts with an explicit candidate radius stored on `nuisance_view`.",
          c("cdid_ipw_enclosure_candidate_mismatch_error", "cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_radius_error", "cdid_ipw_enclosure_input_error")
        )
      }
      cand_r <- as.double(vc$r)
      cand_r[cand_r == 0.0] <- 0.0
      if (!identical(r_dbl, cand_r)) {
        cdid_ipw_enclosure_abort(
          "Explicit `r` disagrees with optional candidate stored on `nuisance_view`.",
          c("cdid_ipw_enclosure_candidate_mismatch_error", "cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_radius_error", "cdid_ipw_enclosure_input_error")
        )
      }
    }
  }

  r_ctrl_with_attrs <- r_ctrl
  if (!is.null(expected_a_id)) attr(r_ctrl_with_attrs, "expected_a_id") <- expected_a_id
  if (!is.null(trace_env)) attr(r_ctrl_with_attrs, "trace_env") <- trace_env

  # The typed request itself (never a numeric placeholder) is forwarded to the
  # certificate; the certificate selects the one radius. The selected binary64
  # scalar `cert$r` is then adopted as the problem radius so that `raw$r`,
  # `certificate$r` and `provenance$radius_r` are one identical value for every
  # consumer identity check. An INCONCLUSIVE certificate keeps `r = NA_real_`
  # inside the certificate object only; the problem is not certified and the
  # pipeline refuses it before any consumer sees a radius.
  cert <- certify_chained_logit_root_view(
    nuisance_view = nuisance_view,
    beta0 = beta0_dbl,
    F = F_dbl,
    r = if (r_auto) r else r_dbl,
    control = cert_ctrl,
    reader_control = r_ctrl_with_attrs
  )
  if (r_auto) {
    if (!identical(cert$radius_selection$mode, .logit_cert_automatic_radius_mode) ||
        !identical(cert$r, cert$provenance$radius_r)) {
      cdid_ipw_enclosure_abort(
        "Automatic-radius certificate did not record one selected radius consistently.",
        c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error")
      )
    }
    r_dbl <- as.double(cert$r)
  }

  ctx <- construct_view_ipw_nuisance_access_context(
    nuisance_view = nuisance_view,
    coh_dbl = coh_dbl,
    cohort_labels = cohort_labels,
    beta0_dbl = beta0_dbl,
    F_dbl = F_dbl,
    r_dbl = r_dbl,
    cert_ctrl = cert_ctrl,
    prob_ctrl = prob_ctrl,
    reader_ctrl = r_ctrl,
    preflight = preflight,
    expected_a_id = expected_a_id,
    trace_env = trace_env
  )
  tr_res <- ipw_enclosure_nuisance_transfer_with_context(ctx, cert)
  nuisance_enclosures <- tr_res$nuisance_enclosures
  arith_refusal <- tr_res$arith_refusal
  extra_reasons <- tr_res$extra_reasons
  added_ops <- tr_res$added_ops
  added_reader_ops <- tr_res$added_reader_ops
  added_endpoints_read <- tr_res$added_endpoints_read
  added_chunks <- tr_res$added_chunks
  added_passes <- tr_res$added_passes
  added_solves <- tr_res$added_solves

  total_interval_ops <- cert$work$actual_interval_ops + added_ops
  total_reader_ops <- cert$work$actual_reader_work_ops + added_reader_ops
  total_endpoints_read <- cert$work$endpoints_read + added_endpoints_read

  stopifnot(
    added_ops <= preflight$added_work_ops_bound,
    total_interval_ops <= preflight$total_work_ops_bound,
    total_reader_ops <= preflight$total_reader_work_ops_bound,
    total_interval_ops + total_reader_ops <= preflight$combined_work_ops_bound,
    total_endpoints_read == 0L
  )

  certified <- isTRUE(cert$certified) && !is.null(nuisance_enclosures)
  inconclusive_reasons <- unique(c(cert$inconclusive_reasons, extra_reasons))
  inconclusive_type <- if (certified) {
    "none"
  } else if (!is.null(arith_refusal) || identical(cert$inconclusive_type, "arithmetic")) {
    "arithmetic"
  } else {
    cert$inconclusive_type
  }

  cert_sealed <- ipw_enclosure_seal(cert, "certificate")

  problem_obj <- structure(
    list(
      storage_mode = "compact_view",
      certified = certified,
      status = if (certified) "CERTIFIED" else "INCONCLUSIVE",
      inconclusive = !certified,
      inconclusive_type = inconclusive_type,
      inconclusive_reasons = inconclusive_reasons,
      separation_certified = FALSE,
      sampling_certified = FALSE,
      public_gate_adopted = FALSE,
      nuisance_ref = nuisance_view$ref,
      frame_ref = nuisance_view$frame_ref,
      panel_ref = nuisance_view$panel_ref,
      defining_pair = nuisance_view$defining_pair,
      raw = list(
        n_a_units = as.integer(n),
        p_covariates = as.integer(p),
        z_colnames = z_colnames,
        id = NULL,
        id_key = NULL,
        id_type = panel_owner$id_type,
        id_storage = panel_owner$id_storage,
        Z = NULL,
        w = NULL,
        G = NULL,
        pair_observed = NULL,
        cohorts = coh_dbl,
        cohort_labels = cohort_labels,
        # Preserve supplied candidate bits; arithmetic above uses beta0_dbl.
        beta0 = as.double(beta0),
        F = F_dbl,
        r = r_dbl
      ),
      masks = list(
        cohort_index = NULL,
        in_risk = NULL,
        in_observed = NULL,
        risk_indices = NULL,
        observed_indices = NULL,
        by_cohort_risk_count = cert$dimensions$by_cohort_risk_count,
        by_cohort_observed_count = cert$dimensions$by_cohort_observed_count
      ),
      dimensions = cert$dimensions,
      certificate = cert,
      certificate_sealed = cert_sealed,
      nuisance_enclosures = nuisance_enclosures,
      arithmetic_refusal = if (!is.null(arith_refusal)) {
        list(
          occurred = TRUE,
          reason = if (!is.null(arith_refusal$reason)) arith_refusal$reason else "arithmetic_refusal",
          message = conditionMessage(arith_refusal),
          condition_class = class(arith_refusal),
          diagnostics = arith_refusal$diagnostics
        )
      } else {
        cert$arithmetic_refusal
      },
      controls = list(
        certificate_control = cert_ctrl,
        problem_control = prob_ctrl,
        reader_control = r_ctrl
      ),
      preflight = preflight,
      read_diagnostics = list(
        read_outcomes = FALSE,
        endpoints_read = total_endpoints_read,
        n_endpoints_read = total_endpoints_read,
        risk_filtered_before_outcome_read = TRUE,
        off_risk_endpoints_read = FALSE,
        nonobserver_endpoints_read = FALSE,
        cached_rounded_dy_returned = FALSE
      ),
      work = list(
        raw_certificate_work = cert$work,
        added_passes_completed = added_passes,
        added_chunks_processed = logit_cert_count_output(added_chunks),
        added_triangular_solves = logit_cert_count_output(added_solves),
        added_interval_ops = added_ops,
        total_interval_ops = total_interval_ops,
        added_reader_work_ops = added_reader_ops,
        total_reader_work_ops = total_reader_ops,
        total_combined_work_ops = total_interval_ops + total_reader_ops,
        endpoints_read = total_endpoints_read,
        interval_ops_definition = "Conservative interval-operation charges at completed helper boundaries; excludes CPU instructions, R object headers, and process peak RSS."
      ),
      provenance = list(
        storage_mode = "compact_view",
        normalizer_n = as.double(n),
        normalizer_definition = "count of all eligible A-units (including outside-risk and pair non-observers)",
        cohorts = coh_dbl,
        cohort_labels = cohort_labels,
        target_cohort = coh_dbl[1L],
        control_cohorts = coh_dbl[-1L],
        defining_pair = nuisance_view$defining_pair,
        immutability_contract = "Compact O(p^2 + p*m + m) nuisance state bound to sealed nuisance/frame/panel snapshots; zero owned Z copy or full-N masks retained."
      )
    ),
    class = "cdid_ipw_nuisance_problem"
  )

  ipw_enclosure_seal(problem_obj, "nuisance_problem")
}

certify_ipw_nuisance_problem <- construct_certified_ipw_nuisance_problem
construct_ipw_enclosure_problem <- construct_certified_ipw_nuisance_problem
certify_ipw_nuisance_problem_view <- construct_certified_ipw_nuisance_problem_view
construct_ipw_enclosure_problem_view <- construct_certified_ipw_nuisance_problem_view

ipw_validate_nuisance_problem_mode <- function(
  problem,
  require_certified = FALSE,
  expected_storage_mode = NULL
) {
  ipw_enclosure_require_deps()
  if (missing(problem) || is.null(problem) || !inherits(problem, "cdid_ipw_nuisance_problem") ||
      !is.list(problem) || !ipw_enclosure_seal_ok(problem, "nuisance_problem")) {
    cdid_ipw_enclosure_abort(
      "Nuisance problem is malformed, foreign, or differs from its sealed snapshot; reuse is refused.",
      c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error")
    )
  }
  actual_mode <- if (identical(problem$storage_mode, "compact_view")) "compact_view" else "dense"
  if (!is.null(expected_storage_mode) && !identical(actual_mode, expected_storage_mode)) {
    cdid_ipw_enclosure_abort(
      paste0("Expected nuisance problem with storage mode `", expected_storage_mode, "`, but got `", actual_mode, "`."),
      c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error", "cdid_ipw_enclosure_input_error")
    )
  }

  if (!is.list(problem$certificate_sealed) ||
      !ipw_enclosure_seal_ok(problem$certificate_sealed, "certificate") ||
      !inherits(problem$certificate, "cdid_logit_certificate")) {
    cdid_ipw_enclosure_abort(
      "Nuisance problem certificate seal is missing or altered.",
      c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error")
    )
  }
  cert_unsealed <- problem$certificate_sealed
  cert_unsealed$ref <- NULL
  if (!identical(problem$certificate, cert_unsealed)) {
    cdid_ipw_enclosure_abort(
      "Nuisance problem certificate differs from its sealed snapshot.",
      c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error")
    )
  }

  raw <- problem$raw
  cert <- problem$certificate
  prov <- cert$provenance
  dims <- cert$dimensions

  if (identical(actual_mode, "compact_view")) {
    if (!is.list(raw) || !is.list(prov) || !is.list(dims) ||
        !is.null(raw$Z) || !is.null(raw$id) || !is.null(raw$id_key) ||
        !is.null(raw$w) || !is.null(raw$G) || !is.null(raw$pair_observed) ||
        !is.null(problem$masks$cohort_index) || !is.null(problem$masks$in_risk) ||
        !is.null(problem$masks$in_observed) || !is.null(problem$masks$risk_indices) ||
        !is.null(problem$masks$observed_indices)) {
      cdid_ipw_enclosure_abort(
        "Compact nuisance problem must not retain dense `Z`, `id`, `w`, `G`, `pair_observed`, or full-N mask arrays.",
        c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error")
      )
    }
    nref <- problem$nuisance_ref
    fref <- problem$frame_ref
    pref <- problem$panel_ref
    if (!is.environment(nref) || !environmentIsLocked(nref) ||
        !exists("nuisance_view", envir = nref, inherits = FALSE) ||
        !is.environment(fref) || !environmentIsLocked(fref) ||
        !exists("frame_owner", envir = fref, inherits = FALSE) ||
        !is.environment(pref) || !environmentIsLocked(pref) ||
        !exists("panel_owner", envir = pref, inherits = FALSE)) {
      cdid_ipw_enclosure_abort(
        "Compact nuisance problem lacks valid sealed `nuisance_ref`, `frame_ref`, or `panel_ref`.",
        c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error")
      )
    }
    nuisance_view <- get("nuisance_view", envir = nref, inherits = FALSE)
    frame_owner <- get("frame_owner", envir = fref, inherits = FALSE)
    validate_ipw_nuisance_view(nuisance_view, frame_owner = frame_owner)
    if (!identical(pref, frame_owner$panel_ref) ||
        !identical(as.double(raw$n_a_units), prov$normalizer_n) ||
        !identical(as.integer(raw$n_a_units), as.integer(frame_owner$n_a_units)) ||
        !identical(as.integer(raw$p_covariates), as.integer(frame_owner$n_covariates)) ||
        !identical(raw$z_colnames, colnames(frame_owner$Z_a)) ||
        !identical(dims$n_a_units, as.integer(raw$n_a_units)) ||
        !identical(raw$cohorts, prov$cohorts) ||
        !identical(raw$cohorts, nuisance_view$cohorts$cohorts) ||
        !identical(raw$cohort_labels, prov$cohort_labels) ||
        !identical(raw$beta0, prov$beta0) ||
        !identical(raw$F, prov$F_factor) ||
        !identical(raw$r, prov$radius_r) ||
        !identical(problem$controls$certificate_control, prov$control)) {
      cdid_ipw_enclosure_abort(
        "Compact nuisance problem header or bound owner references disagree with the retained root certificate provenance.",
        c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error")
      )
    }
  } else {
    if (!is.list(raw) || !is.list(prov) || !is.list(dims) ||
        !identical(as.double(nrow(raw$Z)), prov$normalizer_n) ||
        !identical(raw$cohorts, prov$cohorts) ||
        !identical(raw$cohort_labels, prov$cohort_labels) ||
        !identical(raw$beta0, prov$beta0) ||
        !identical(raw$F, prov$F_factor) ||
        !identical(raw$r, prov$radius_r) ||
        !identical(problem$controls$certificate_control, prov$control)) {
      cdid_ipw_enclosure_abort(
        "Nuisance problem raw inputs disagree with the retained root certificate provenance.",
        c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error")
      )
    }

    # Constructor verified all raw arrays and derived masks before sealing.
    # Under ordinary R copy-on-modify, an unchanged locked snapshot retains those
    # identities. Do not allocate/recompute full-n match/which masks per row chunk.
    # This is an integrity check for ordinary objects, not an adversarial security boundary.
    if (!identical(dims$n_a_units, as.integer(nrow(raw$Z))) ||
        length(problem$masks$cohort_index) != nrow(raw$Z) ||
        length(problem$masks$in_risk) != nrow(raw$Z) ||
        length(problem$masks$in_observed) != nrow(raw$Z)) {
      cdid_ipw_enclosure_abort("Sealed nuisance dimensions disagree with raw array headers.",
        c("cdid_ipw_enclosure_problem_error", "cdid_ipw_enclosure_integrity_error"))
    }
  }

  if (isTRUE(require_certified) && !isTRUE(problem$certified)) {
    cdid_ipw_enclosure_abort(
      "Certified nuisance problem is required, but problem status is INCONCLUSIVE.",
      c("cdid_ipw_enclosure_inconclusive", "cdid_ipw_enclosure_problem_error"),
      diagnostics = list(inconclusive_reasons = problem$inconclusive_reasons)
    )
  }
  invisible(problem)
}

validate_certified_ipw_nuisance_problem <- function(problem, require_certified = FALSE) {
  ipw_validate_nuisance_problem_mode(problem, require_certified = require_certified)
}
validate_ipw_nuisance_problem <- validate_certified_ipw_nuisance_problem

# ---------------------------------------------------------------------------
# Component 2: Single-link outcome functional, loading, and Theorem B / Cor. sep enclosure
# ---------------------------------------------------------------------------

ipw_enclosure_validate_candidates <- function(candidates, p, m) {
  if (is.null(candidates)) return(NULL)
  if (!is.list(candidates) || is.object(candidates)) {
    cdid_ipw_enclosure_abort(
      "`candidates` must be NULL or a plain named list of stored candidate quantities.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  allowed <- c(
    "mu", "treated_mean", "control_mean",
    "G_beta", "v", "G_r", "gamma", "A_beta_r", "L",
    "loading_beta", "loading_rate", "H_chol"
  )
  if (is.null(names(candidates)) || any(!nzchar(names(candidates))) ||
      anyDuplicated(names(candidates)) || any(!names(candidates) %in% allowed)) {
    cdid_ipw_enclosure_abort(
      paste0("Invalid `candidates` fields. Allowed fields: ", paste(allowed, collapse = ", "), "."),
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  out <- list()
  for (sc in c("mu", "treated_mean", "control_mean")) {
    if (!is.null(candidates[[sc]])) {
      if (!gmm_scalar_numeric(candidates[[sc]])) {
        cdid_ipw_enclosure_abort(
          paste0("`candidates$", sc, "` must be a finite numeric scalar."),
          c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
        )
      }
      out[[sc]] <- as.double(candidates[[sc]])
    }
  }
  if (!is.null(out$mu) && !is.null(out$treated_mean) && !is.null(out$control_mean)) {
    # Check candidate internal consistency up to binary64 subtraction
    if (abs((out$treated_mean - out$control_mean) - out$mu) > 16 * .Machine$double.eps * max(1, abs(out$treated_mean), abs(out$control_mean))) {
      cdid_ipw_enclosure_abort(
        "`candidates$mu` is inconsistent with `candidates$treated_mean - candidates$control_mean`.",
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  for (vp in c("G_beta", "v", "loading_beta")) {
    if (!is.null(candidates[[vp]])) {
      x <- candidates[[vp]]
      if (!is.atomic(x) || !is.null(dim(x)) || is.object(x) || !is.numeric(x) ||
          length(x) != p || anyNA(x) || any(!is.finite(x))) {
        cdid_ipw_enclosure_abort(
          paste0("`candidates$", vp, "` must be a finite numeric vector of length p = ", p, "."),
          c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
        )
      }
      out[[vp]] <- as.double(x)
    }
  }
  if (!is.null(out$G_beta) && !is.null(out$v) && !identical(out$G_beta, -out$v)) {
    cdid_ipw_enclosure_abort(
      "`candidates$G_beta` and `candidates$v` disagree (`G_beta` must equal `-v`).",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  for (vm in c("G_r", "gamma", "loading_rate")) {
    if (!is.null(candidates[[vm]])) {
      x <- candidates[[vm]]
      if (!is.atomic(x) || !is.null(dim(x)) || is.object(x) || !is.numeric(x) ||
          length(x) != m || anyNA(x) || any(!is.finite(x))) {
        cdid_ipw_enclosure_abort(
          paste0("`candidates$", vm, "` must be a finite numeric vector of length m = ", m, "."),
          c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
        )
      }
      if (vm %in% c("G_r", "gamma") && as.double(x[1L]) != 0.0) {
        cdid_ipw_enclosure_abort(
          paste0("Target cohort entry `candidates$", vm, "[1]` must be exactly 0."),
          c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
        )
      }
      out[[vm]] <- as.double(x)
    }
  }
  if (!is.null(out$G_r) && !is.null(out$gamma) && !identical(out$G_r, out$gamma)) {
    cdid_ipw_enclosure_abort(
      "`candidates$G_r` and `candidates$gamma` disagree.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  for (ml in c("A_beta_r", "L")) {
    if (!is.null(candidates[[ml]])) {
      M_cand <- candidates[[ml]]
      if (!is.matrix(M_cand) || is.object(M_cand) || !is.numeric(M_cand) ||
          nrow(M_cand) != p || ncol(M_cand) != m || anyNA(M_cand) || any(!is.finite(M_cand))) {
        cdid_ipw_enclosure_abort(
          paste0("`candidates$", ml, "` must be a finite numeric p x m matrix."),
          c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
        )
      }
      out[[ml]] <- matrix(as.double(M_cand), nrow = p, ncol = m)
    }
  }
  if (!is.null(out$A_beta_r) && !is.null(out$L) && !identical(out$A_beta_r, out$L)) {
    cdid_ipw_enclosure_abort(
      "`candidates$A_beta_r` and `candidates$L` disagree.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.null(candidates$H_chol)) {
    out$H_chol <- validate_logit_certificate_metric_factor(candidates$H_chol, p)
  }
  out
}

construct_dense_ipw_link_access_context <- function(
  problem,
  ctrl,
  preflight,
  step,
  tilde_lambda,
  cand_norm,
  link_id,
  functional_key,
  obs_idx,
  a_to_obs_pos,
  y_left_obs,
  y_right_obs,
  validation_chunks
) {
  raw <- problem$raw
  n <- nrow(raw$Z)
  p <- ncol(raw$Z)
  m <- length(raw$cohorts)
  cohort_idx <- problem$masks$cohort_index

  read_chunk <- function(a_start, a_end) {
    a_seq <- seq.int(as.integer(a_start), as.integer(a_end))
    obs_pos_ch <- a_to_obs_pos[a_seq]
    obs_local <- which(obs_pos_ch > 0L)
    obs_global <- obs_pos_ch[obs_local]
    chunk_to_obs_pos <- integer(length(a_seq))
    if (length(obs_local) > 0L) chunk_to_obs_pos[obs_local] <- seq_along(obs_local)
    list(
      a_index = a_seq,
      Z_chunk = raw$Z[a_seq, , drop = FALSE],
      weight = raw$w[a_seq],
      cohort_index = cohort_idx[a_seq],
      chunk_to_obs_pos = chunk_to_obs_pos,
      y_left_obs = if (length(obs_global) > 0L) y_left_obs[obs_global] else numeric(0L),
      y_right_obs = if (length(obs_global) > 0L) y_right_obs[obs_global] else numeric(0L),
      n_endpoints_read = 2L * length(obs_global),
      reader_work_ops = 0
    )
  }

  ipw_seal_access_context(
    list(
      mode = "dense",
      problem = problem,
      header = list(n = n, p = p, m = m, z_colnames = colnames(raw$Z), cohort_labels = raw$cohort_labels),
      ctrl = ctrl,
      reader_ctrl = NULL,
      preflight = preflight,
      step = step,
      tilde_lambda = tilde_lambda,
      cand_norm = cand_norm,
      link_id = link_id,
      functional_key = functional_key,
      link_ref = NULL,
      defining_pair = NULL,
      endpoints_list = list(
        obs_a_index = obs_idx,
        a_to_obs_pos = a_to_obs_pos,
        y_left_obs = y_left_obs,
        y_right_obs = y_right_obs,
        off_observed_endpoints_read = FALSE
      ),
      validation_chunks = validation_chunks,
      read_chunk = read_chunk
    ),
    phase = "link_enclosure"
  )
}

construct_view_ipw_link_access_context <- function(
  problem,
  link_view,
  ctrl,
  reader_ctrl,
  preflight,
  step,
  tilde_lambda,
  cand_norm,
  expected_a_id = NULL,
  trace_env = NULL
) {
  raw <- problem$raw
  n <- as.integer(raw$n_a_units)
  p <- as.integer(raw$p_covariates)
  m <- length(raw$cohorts)

  read_chunk <- function(a_start, a_end) {
    ch <- query_ipw_view_chunk(
      view = link_view,
      a_start = a_start,
      a_end = a_end,
      control = reader_ctrl,
      read_outcomes = isTRUE(problem$certified),
      include_design = isTRUE(problem$certified),
      enclose_differences = FALSE,
      expected_a_id = expected_a_id,
      trace_env = trace_env
    )
    if (!is.null(trace_env) && isTRUE(problem$certified)) {
      trace_env$obs_s_read <- ch$read_diagnostics$read_obs_s_index
      trace_env$obs_t_read <- ch$read_diagnostics$read_obs_t_index
    }
    list(
      a_index = ch$a_index,
      Z_chunk = ch$design$Z_chunk,
      weight = ch$weight,
      cohort_index = ch$cohort_index,
      chunk_to_obs_pos = ch$chunk_to_obs_pos,
      y_left_obs = ch$observed_risk$y_left_obs,
      y_right_obs = ch$observed_risk$y_right_obs,
      n_endpoints_read = ch$read_diagnostics$n_endpoints_read,
      reader_work_ops = ch$work$actual_work_ops
    )
  }

  ipw_seal_access_context(
    list(
      mode = "compact_view",
      problem = problem,
      header = list(n = n, p = p, m = m, z_colnames = raw$z_colnames, cohort_labels = raw$cohort_labels),
      ctrl = ctrl,
      reader_ctrl = reader_ctrl,
      preflight = preflight,
      step = step,
      tilde_lambda = tilde_lambda,
      cand_norm = cand_norm,
      link_id = link_view$link_id,
      functional_key = link_view$functional_key,
      link_ref = link_view$ref,
      defining_pair = link_view$pair,
      endpoints_list = list(
        obs_a_index = NULL,
        a_to_obs_pos = NULL,
        y_left_obs = NULL,
        y_right_obs = NULL,
        n_observed_risk_units = problem$dimensions$n_observed_risk_units,
        retained_endpoint_entries = 0L,
        retained_position_map_entries = 0L,
        off_observed_endpoints_read = FALSE
      ),
      validation_chunks = 0L,
      read_chunk = read_chunk
    ),
    phase = "link_enclosure"
  )
}

ipw_enclosure_link_run_with_context <- function(ctx) {
  ipw_validate_access_context(ctx, expected_phase = "link_enclosure")

  problem <- ctx$problem
  raw <- problem$raw
  n <- ctx$header$n
  p <- ctx$header$p
  m <- ctx$header$m
  z_colnames <- ctx$header$z_colnames
  cohort_labels <- ctx$header$cohort_labels
  ctrl <- ctx$ctrl
  preflight <- ctx$preflight
  step <- ctx$step
  tilde_lambda <- ctx$tilde_lambda
  cand_norm <- ctx$cand_norm
  link_id <- ctx$link_id
  functional_key <- ctx$functional_key
  endpoints_list <- ctx$endpoints_list
  validation_chunks <- ctx$validation_chunks
  is_compact <- identical(ctx$mode, "compact_view")

  total_reader_ops <- 0
  total_endpoints_read <- 0
  reader_passes_completed <- 0L

  # If the nuisance problem itself is INCONCLUSIVE, in compact_view mode we still
  # run 1 outcome-free pass of ctx$read_chunk so that any structural expected_a_id
  # refusal takes precedence over INCONCLUSIVE, with 0 outcome endpoints read.
  if (!isTRUE(problem$certified)) {
    if (is_compact) {
      a_start <- 1L
      while (a_start <= n) {
        ipw_enclosure_check_interrupt()
        a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
        ch <- ctx$read_chunk(a_start, a_end)
        total_reader_ops <- total_reader_ops + ch$reader_work_ops
        total_endpoints_read <- total_endpoints_read + ch$n_endpoints_read
        validation_chunks <- validation_chunks + 1
        a_start <- as.double(a_end) + 1
      }
      reader_passes_completed <- 1L
    }
    inc_fields <- list(
      certified = FALSE,
      status = "INCONCLUSIVE",
      inconclusive = TRUE,
      inconclusive_type = problem$inconclusive_type,
      inconclusive_reasons = unique(c("nuisance_problem_inconclusive", problem$inconclusive_reasons)),
      separation_certified = FALSE,
      sampling_certified = FALSE,
      public_gate_adopted = FALSE,
      storage_mode = if (is_compact) "compact_view" else NULL,
      link_id = link_id,
      functional_key = functional_key,
      link_ref = if (is_compact) ctx$link_ref else NULL,
      nuisance_ref = if (is_compact) problem$nuisance_ref else NULL,
      frame_ref = if (is_compact) problem$frame_ref else NULL,
      panel_ref = if (is_compact) problem$panel_ref else NULL,
      problem_ref = problem$ref,
      endpoints = endpoints_list,
      loading_beta_candidate = tilde_lambda,
      candidates = cand_norm,
      enclosures_star = NULL,
      enclosures_0 = NULL,
      root_perturbation_bounds = NULL,
      separated_error_budget = NULL,
      arithmetic_refusal = problem$arithmetic_refusal,
      control = ctrl,
      reader_control = if (is_compact) ctx$reader_ctrl else NULL,
      preflight = preflight,
      read_diagnostics = if (is_compact) list(
        endpoints_read = as.integer(total_endpoints_read),
        n_endpoints_read = as.integer(total_endpoints_read),
        reader_passes_completed = reader_passes_completed,
        actual_reader_work_ops = total_reader_ops
      ) else NULL,
      work = if (is_compact) {
        list(
          validation_chunks_processed = logit_cert_count_output(validation_chunks),
          arithmetic_passes_completed = 0L,
          arithmetic_chunks_processed = 0L,
          reader_passes_completed = reader_passes_completed,
          triangular_solves = 0L,
          actual_interval_ops = 0,
          actual_reader_work_ops = total_reader_ops,
          combined_work_ops = total_reader_ops,
          endpoints_read = logit_cert_count_output(total_endpoints_read)
        )
      } else {
        list(
          validation_chunks_processed = logit_cert_count_output(validation_chunks),
          arithmetic_passes_completed = 0L,
          arithmetic_chunks_processed = 0L,
          triangular_solves = 0L,
          actual_interval_ops = 0
        )
      }
    )
    if (!is_compact) {
      inc_fields$storage_mode <- NULL
      inc_fields$link_ref <- NULL
      inc_fields$nuisance_ref <- NULL
      inc_fields$frame_ref <- NULL
      inc_fields$panel_ref <- NULL
      inc_fields$reader_control <- NULL
      inc_fields$read_diagnostics <- NULL
    }
    inc_obj <- structure(inc_fields, class = "cdid_ipw_link_enclosure")
    return(ipw_enclosure_seal(inc_obj, "link_enclosure"))
  }

  # 5. Streamed Pass 1 (T, C^0, C^*, mu^0, mu^*, omega) & Pass 2 (v^0, v^*, gamma^0, gamma^*)
  #    plus Lemma dual(iii) loadings and Theorem B / Cor. sep bounds
  nuis <- problem$nuisance_enclosures
  cert <- problem$certificate
  M_cohort <- cert$enclosures$M_cohort
  Pi_cohort <- cert$enclosures$Pi_cohort
  r_cohort <- cert$enclosures$r_cohort
  q_cohort <- cert$enclosures$q_cohort
  M_over_Pi_cohort <- cert$enclosures$M_over_Pi_cohort
  F_dbl <- raw$F
  F_inv_T <- nuis$F_inv_T
  coord_dual_iv <- nuis$coord_dual_iv
  beta0_dbl <- raw$beta0
  cert_ctrl <- problem$controls$certificate_control
  n_iv <- gmm_iv(as.double(n))

  actual_ops <- 0
  arith_chunks <- 0
  arith_passes <- 0L
  tri_solves <- 0L
  exp_evals <- 0L

  link_try <- tryCatch({
    T_num_iv <- c(0.0, 0.0)
    C0_num_iv <- c(0.0, 0.0)
    Cstar_num_iv <- c(0.0, 0.0)
    min_Delta_c <- c(Inf, Inf)
    max_Delta_c <- c(-Inf, -Inf)

    a_hi <- cert$a
    a_c_hi <- cert$a_control
    rho_hi <- cert$rho
    delta_hi <- cert$delta
    x_bar_hi <- nuis$x_bar
    x_bar_c_hi <- nuis$x_bar_c

    # Pass 1: Streamed outcome means T, C^0, C^* and control outcome range omega.
    # In compact_view mode, if an arithmetic refusal occurs on an earlier chunk,
    # we finish reading remaining Pass 1 chunks so any structural endpoint/ID error
    # on a later chunk still takes precedence over arithmetic INCONCLUSIVE.
    first_pass1_arith_err <- NULL
    a_start <- 1L
    while (a_start <= n) {
      ipw_enclosure_check_interrupt()
      a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
      ch <- ctx$read_chunk(a_start, a_end)
      total_reader_ops <- total_reader_ops + ch$reader_work_ops
      total_endpoints_read <- total_endpoints_read + ch$n_endpoints_read
      if (is_compact) validation_chunks <- validation_chunks + 1

      if (is.null(first_pass1_arith_err)) {
        ch_err <- tryCatch({
          idx <- ch$a_index
          b_len <- length(idx)
          Z_ch <- ch$Z_chunk
          w_ch <- ch$weight
          j_ch <- ch$cohort_index
          obs_pos_ch <- ch$chunk_to_obs_pos
          yl_obs <- ch$y_left_obs
          yr_obs <- ch$y_right_obs

          for (u in seq_len(b_len)) {
            op <- obs_pos_ch[u]
            if (op == 0L) next
            ju <- j_ch[u]
            D_u <- if (ju == 1L) 1L else 0L
            wu_iv <- gmm_iv(w_ch[u])
            z_u <- Z_ch[u, ]

            # Real endpoint difference Delta_i = y_{it} - y_{is} enclosed outward
            Delta_iv <- gmm_iv_sub(gmm_iv(yr_obs[op]), gmm_iv(yl_obs[op]))
            actual_ops <- actual_ops + 1

            cw <- ipw_enclosure_c_weight(
              wu_iv, M_cohort[ju, ], Pi_cohort[ju, ], M_over_Pi_cohort[ju, ],
              unit_offset = idx[u], cohort_index = ju
            )
            c_iv <- cw$c_iv
            actual_ops <- actual_ops + cw$ops

            if (D_u == 1L) {
              T_num_iv <- gmm_iv_add(T_num_iv, gmm_iv_mul(c_iv, Delta_iv))
              actual_ops <- actual_ops + 2
            } else {
              if (Delta_iv[1L] < min_Delta_c[1L]) min_Delta_c[1L] <- Delta_iv[1L]
              if (Delta_iv[2L] < min_Delta_c[2L]) min_Delta_c[2L] <- Delta_iv[2L]
              if (Delta_iv[1L] > max_Delta_c[1L]) max_Delta_c[1L] <- Delta_iv[1L]
              if (Delta_iv[2L] > max_Delta_c[2L]) max_Delta_c[2L] <- Delta_iv[2L]

              if (all(z_u == 0.0)) {
                ai_hi <- 0.0
              } else {
                dn_z <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, z_u)
                actual_ops <- actual_ops + dn_z$ops
                tri_solves <- tri_solves + dn_z$solves
                ai_hi <- min(dn_z$upper, a_c_hi, a_hi)
              }
              xi_hi <- min(x_bar_hi, logit_b64_mul(gmm_iv(ai_hi), gmm_iv(rho_hi))[2L])
              one_m_xi_lo <- logit_b64_sub(gmm_iv(1.0), gmm_iv(xi_hi))[1L]
              actual_ops <- actual_ops + 2

              eta_iv <- gmm_iv_dot(z_u, beta0_dbl)
              actual_ops <- actual_ops + 2 * p
              exp_res <- logit_iv_exp(
                eta_iv,
                taylor_degree = cert_ctrl$taylor_degree,
                max_abs_eta = cert_ctrl$max_abs_eta,
                max_squarings = cert_ctrl$max_squarings
              )
              o0_iv <- exp_res$iv
              actual_ops <- actual_ops + exp_res$ops
              exp_evals <- exp_evals + 2L

              o_star_lo <- logit_b64_mul(gmm_iv(o0_iv[1L]), gmm_iv(one_m_xi_lo))[1L]
              o_star_hi <- logit_b64_div(gmm_iv(o0_iv[2L]), gmm_iv(one_m_xi_lo))[2L]
              o_star_iv <- c(o_star_lo, o_star_hi)
              actual_ops <- actual_ops + 2

              b0_u <- gmm_iv_mul(c_iv, o0_iv)
              b_star_u <- gmm_iv_mul(c_iv, o_star_iv)
              C0_num_iv <- gmm_iv_add(C0_num_iv, gmm_iv_mul(b0_u, Delta_iv))
              Cstar_num_iv <- gmm_iv_add(Cstar_num_iv, gmm_iv_mul(b_star_u, Delta_iv))
              actual_ops <- actual_ops + 6
            }
          }
          arith_chunks <- arith_chunks + 1
          NULL
        }, cdid_ipw_enclosure_arithmetic_inconclusive = function(e) e,
           cdid_logit_arithmetic_inconclusive = function(e) e,
           cdid_gmm_accuracy_error = function(e) e)

        if (inherits(ch_err, "error")) {
          first_pass1_arith_err <- ch_err
          if (!is_compact) stop(ch_err)
        }
      }
      a_start <- as.double(a_end) + 1
    }
    reader_passes_completed <- 1L
    if (!is.null(first_pass1_arith_err)) stop(first_pass1_arith_err)
    arith_passes <- 1L

    Sigma_plus_iv <- nuis$Sigma_plus
    Sigma_minus_0_iv <- nuis$Sigma_minus_0
    Sigma_minus_star_iv <- nuis$Sigma_minus_star

    T_iv <- gmm_iv_div(T_num_iv, Sigma_plus_iv)
    C_0_iv <- gmm_iv_div(C0_num_iv, Sigma_minus_0_iv)
    C_star_direct_iv <- gmm_iv_div(Cstar_num_iv, Sigma_minus_star_iv)
    actual_ops <- actual_ops + 3

    omega_iv <- gmm_iv_sub(max_Delta_c, min_Delta_c)
    omega_iv[1L] <- max(0.0, omega_iv[1L])
    omega_hi <- omega_iv[2L]
    actual_ops <- actual_ops + 1

    # Point root perturbation bound e_mu^root = upper(omega * a_c * rho) (Cor. point / Lemma mult.5)
    point_root_bound <- logit_b64_mul(
      logit_b64_mul(gmm_iv(omega_hi), gmm_iv(a_c_hi)),
      gmm_iv(rho_hi)
    )[2L]
    actual_ops <- actual_ops + 2

    # Theorem A step 2: intersect direct interval evaluation of C^* with [C^0] +/- omega * a_c * rho
    C_star_ball_iv <- ipw_enclosure_widen_iv(C_0_iv, point_root_bound)
    actual_ops <- actual_ops + 2
    C_star_iv <- logit_iv_intersect(C_star_direct_iv, C_star_ball_iv, "C_star = C([o^*]) intersect ([C^0] +/- omega*a_c*rho)")

    mu_0_iv <- gmm_iv_sub(T_iv, C_0_iv)
    mu_star_iv <- gmm_iv_sub(T_iv, C_star_iv)
    actual_ops <- actual_ops + 2

    # Pass 2: Streamed outcome partials v^0, v^*, gamma^0, gamma^*
    v_0_sum_iv <- matrix(0.0, nrow = p, ncol = 2L)
    v_star_sum_iv <- matrix(0.0, nrow = p, ncol = 2L)
    gamma_0_sum_iv <- matrix(0.0, nrow = m, ncol = 2L)
    gamma_star_sum_iv <- matrix(0.0, nrow = m, ncol = 2L)

    a_start <- 1L
    while (a_start <= n) {
      ipw_enclosure_check_interrupt()
      a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
      ch <- ctx$read_chunk(a_start, a_end)
      total_reader_ops <- total_reader_ops + ch$reader_work_ops
      total_endpoints_read <- total_endpoints_read + ch$n_endpoints_read

      idx <- ch$a_index
      b_len <- length(idx)
      Z_ch <- ch$Z_chunk
      w_ch <- ch$weight
      j_ch <- ch$cohort_index
      obs_pos_ch <- ch$chunk_to_obs_pos
      yl_obs <- ch$y_left_obs
      yr_obs <- ch$y_right_obs

      for (u in seq_len(b_len)) {
        op <- obs_pos_ch[u]
        if (op == 0L) next
        ju <- j_ch[u]
        if (ju == 1L) next # v and gamma sum over O_c (controls only)

        wu_iv <- gmm_iv(w_ch[u])
        z_u <- Z_ch[u, ]
        Delta_iv <- gmm_iv_sub(gmm_iv(yr_obs[op]), gmm_iv(yl_obs[op]))
        actual_ops <- actual_ops + 1

        cw <- ipw_enclosure_c_weight(
          wu_iv, M_cohort[ju, ], Pi_cohort[ju, ], M_over_Pi_cohort[ju, ],
          unit_offset = idx[u], cohort_index = ju
        )
        c_iv <- cw$c_iv
        actual_ops <- actual_ops + cw$ops

        if (all(z_u == 0.0)) {
          ai_hi <- 0.0
        } else {
          dn_z <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, z_u)
          actual_ops <- actual_ops + dn_z$ops
          tri_solves <- tri_solves + dn_z$solves
          ai_hi <- min(dn_z$upper, a_c_hi, a_hi)
        }
        xi_hi <- min(x_bar_hi, logit_b64_mul(gmm_iv(ai_hi), gmm_iv(rho_hi))[2L])
        one_m_xi_lo <- logit_b64_sub(gmm_iv(1.0), gmm_iv(xi_hi))[1L]
        actual_ops <- actual_ops + 2

        eta_iv <- gmm_iv_dot(z_u, beta0_dbl)
        actual_ops <- actual_ops + 2 * p
        exp_res <- logit_iv_exp(
          eta_iv,
          taylor_degree = cert_ctrl$taylor_degree,
          max_abs_eta = cert_ctrl$max_abs_eta,
          max_squarings = cert_ctrl$max_squarings
        )
        o0_iv <- exp_res$iv
        actual_ops <- actual_ops + exp_res$ops
        exp_evals <- exp_evals + 2L
        o_star_iv <- c(
          logit_b64_mul(gmm_iv(o0_iv[1L]), gmm_iv(one_m_xi_lo))[1L],
          logit_b64_div(gmm_iv(o0_iv[2L]), gmm_iv(one_m_xi_lo))[2L]
        )
        actual_ops <- actual_ops + 2

        b0_u <- gmm_iv_mul(c_iv, o0_iv)
        b_star_u <- gmm_iv_mul(c_iv, o_star_iv)
        res_0_u <- gmm_iv_sub(Delta_iv, C_0_iv)
        res_star_u <- gmm_iv_sub(Delta_iv, C_star_iv)
        b_res_0 <- gmm_iv_mul(b0_u, res_0_u)
        b_res_star <- gmm_iv_mul(b_star_u, res_star_u)
        actual_ops <- actual_ops + 6

        gamma_0_sum_iv[ju, ] <- gmm_iv_add(gamma_0_sum_iv[ju, ], b_res_0)
        gamma_star_sum_iv[ju, ] <- gmm_iv_add(gamma_star_sum_iv[ju, ], b_res_star)
        actual_ops <- actual_ops + 2

        for (k in seq_len(p)) {
          zk_iv <- gmm_iv(z_u[k])
          v_0_sum_iv[k, ] <- gmm_iv_add(v_0_sum_iv[k, ], gmm_iv_mul(b_res_0, zk_iv))
          v_star_sum_iv[k, ] <- gmm_iv_add(v_star_sum_iv[k, ], gmm_iv_mul(b_res_star, zk_iv))
          actual_ops <- actual_ops + 4
        }
      }
      arith_chunks <- arith_chunks + 1
      a_start <- as.double(a_end) + 1
    }
    reader_passes_completed <- 2L
    arith_passes <- 2L

    v_0_iv <- matrix(0.0, nrow = p, ncol = 2L, dimnames = list(z_colnames, c("lo", "hi")))
    v_star_iv <- matrix(0.0, nrow = p, ncol = 2L, dimnames = list(z_colnames, c("lo", "hi")))
    G_beta_0_iv <- matrix(0.0, nrow = p, ncol = 2L, dimnames = list(z_colnames, c("lo", "hi")))
    G_beta_star_iv <- matrix(0.0, nrow = p, ncol = 2L, dimnames = list(z_colnames, c("lo", "hi")))

    for (k in seq_len(p)) {
      v_0_iv[k, ] <- gmm_iv_div(v_0_sum_iv[k, ], Sigma_minus_0_iv)
      v_star_iv[k, ] <- gmm_iv_div(v_star_sum_iv[k, ], Sigma_minus_star_iv)
      G_beta_0_iv[k, ] <- gmm_iv_neg(v_0_iv[k, ])
      G_beta_star_iv[k, ] <- gmm_iv_neg(v_star_iv[k, ])
      actual_ops <- actual_ops + 2
    }

    gamma_0_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
    gamma_star_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
    for (j in seq.int(2L, m)) {
      den_g0 <- gmm_iv_mul(r_cohort[j, ], Sigma_minus_0_iv)
      den_gs <- gmm_iv_mul(r_cohort[j, ], Sigma_minus_star_iv)
      gamma_0_iv[j, ] <- gmm_iv_div(gamma_0_sum_iv[j, ], den_g0)
      gamma_star_iv[j, ] <- gmm_iv_div(gamma_star_sum_iv[j, ], den_gs)
      actual_ops <- actual_ops + 4
    }

    # Dual norms of v^0 and v^* (used for shared xi_beta action bounds)
    dn_v0 <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, v_0_iv)
    dn_vs <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, v_star_iv)
    actual_ops <- actual_ops + dn_v0$ops + dn_vs$ops
    tri_solves <- tri_solves + dn_v0$solves + dn_vs$solves

    # Lemma dual(iii) inverse-Hessian action for lambda_beta = -H^{-1} v
    # At beta0: kappa_0 = 1 - delta; at beta*: kappa_* = (1 - x_bar)(1 - delta)
    lam_sol_0 <- ipw_enclosure_lemma_dual_solve(
      F_dbl = F_dbl, F_inv_T = F_inv_T, coord_dual_iv = coord_dual_iv,
      H_iv = nuis$H_0, rhs_iv = v_0_iv, center_dbl = tilde_lambda,
      kappa_lo = nuis$kappa_0_lower
    )
    lam_sol_star <- ipw_enclosure_lemma_dual_solve(
      F_dbl = F_dbl, F_inv_T = F_inv_T, coord_dual_iv = coord_dual_iv,
      H_iv = nuis$H_star, rhs_iv = v_star_iv, center_dbl = tilde_lambda,
      kappa_lo = nuis$kappa_star_lower
    )
    actual_ops <- actual_ops + lam_sol_0$ops + lam_sol_star$ops
    tri_solves <- tri_solves + lam_sol_0$solves + lam_sol_star$solves

    rownames(lam_sol_0$coord_iv) <- rownames(lam_sol_star$coord_iv) <- z_colnames
    colnames(lam_sol_0$coord_iv) <- colnames(lam_sol_star$coord_iv) <- c("lo", "hi")

    # Enclose lambda_{r, j} = gamma_j + L_j' lambda_beta for ALL cohorts j = 1..m (target included!)
    lambda_r_0_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
    lambda_r_star_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
    rate_coef_s0_0_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
    rate_coef_s1_0_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
    rate_coef_s0_star_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))
    rate_coef_s1_star_iv <- matrix(0.0, nrow = m, ncol = 2L, dimnames = list(cohort_labels, c("lo", "hi")))

    for (j in seq_len(m)) {
      dot_L0 <- ipw_enclosure_ball_dot(
        F_dbl, F_inv_T,
        vec_iv = matrix(nuis$L_0[, j, ], nrow = p, ncol = 2L),
        ball_sol = lam_sol_0,
        vec_dual_upper = nuis$ell_0_upper[j]
      )
      dot_Ls <- ipw_enclosure_ball_dot(
        F_dbl, F_inv_T,
        vec_iv = matrix(nuis$L_star[, j, ], nrow = p, ncol = 2L),
        ball_sol = lam_sol_star,
        vec_dual_upper = nuis$ell_star_upper[j]
      )
      lambda_r_0_iv[j, ] <- gmm_iv_add(gamma_0_iv[j, ], dot_L0$iv)
      lambda_r_star_iv[j, ] <- gmm_iv_add(gamma_star_iv[j, ], dot_Ls$iv)
      actual_ops <- actual_ops + dot_L0$ops + dot_Ls$ops + 2

      # rate_coef_s0 = -(r_j / q_j) * lambda_{r,j}, rate_coef_s1 = ((1 - r_j) / q_j) * lambda_{r,j}
      r_over_q <- gmm_iv_div(r_cohort[j, ], q_cohort[j, ])
      one_m_r_over_q <- gmm_iv_div(gmm_iv_sub(gmm_iv(1.0), r_cohort[j, ]), q_cohort[j, ])
      rate_coef_s0_0_iv[j, ] <- gmm_iv_neg(gmm_iv_mul(r_over_q, lambda_r_0_iv[j, ]))
      rate_coef_s1_0_iv[j, ] <- gmm_iv_mul(one_m_r_over_q, lambda_r_0_iv[j, ])
      rate_coef_s0_star_iv[j, ] <- gmm_iv_neg(gmm_iv_mul(r_over_q, lambda_r_star_iv[j, ]))
      rate_coef_s1_star_iv[j, ] <- gmm_iv_mul(one_m_r_over_q, lambda_r_star_iv[j, ])
      actual_ops <- actual_ops + 7
    }

    # Theorem B root-perturbation bounds (computed outward from beta0 endpoints)
    # 1. ||F tilde_lambda||_2
    F_lam_iv <- matrix(0.0, nrow = p, ncol = 2L)
    for (k in seq_len(p)) {
      F_lam_iv[k, ] <- gmm_iv_dot(F_dbl[k, ], tilde_lambda)
      actual_ops <- actual_ops + 2 * p
    }
    norm_F_lam_hi <- gmm_iv_norm(F_lam_iv[, 1L], F_lam_iv[, 2L])[2L]
    actual_ops <- actual_ops + 3 * p + 2

    Lam0_cand1 <- logit_b64_add(gmm_iv(norm_F_lam_hi), gmm_iv(lam_sol_0$radius_m_norm))[2L]
    omega_ac_hi <- logit_b64_mul(gmm_iv(omega_hi), gmm_iv(a_c_hi))[2L]
    Lam0_cand2 <- logit_b64_div(gmm_iv(omega_ac_hi), gmm_iv(nuis$one_minus_delta_lo))[2L]
    Lambda_0 <- min(Lam0_cand1, Lam0_cand2)
    actual_ops <- actual_ops + 3

    # e_v = omega * a_c^2 * rho
    e_v <- logit_b64_mul(gmm_iv(point_root_bound), gmm_iv(a_c_hi))[2L]
    # e_lambda = (e_v + (1 + delta) * u1(x_bar) * Lambda_0) / ((1 - x_bar)(1 - delta))
    term_elam <- logit_b64_mul(
      logit_b64_mul(gmm_iv(nuis$one_plus_delta_hi), gmm_iv(nuis$u1_x_bar)),
      gmm_iv(Lambda_0)
    )[2L]
    num_elam <- logit_b64_add(gmm_iv(e_v), gmm_iv(term_elam))[2L]
    e_lambda <- logit_b64_div(gmm_iv(num_elam), gmm_iv(nuis$kappa_star_lower))[2L]
    actual_ops <- actual_ops + 5

    # Lambda_* = min(Lambda_0 + e_lambda, omega * a_c / ((1 - x_bar)(1 - delta)))
    Lams_cand1 <- logit_b64_add(gmm_iv(Lambda_0), gmm_iv(e_lambda))[2L]
    Lams_cand2 <- logit_b64_div(gmm_iv(omega_ac_hi), gmm_iv(nuis$kappa_star_lower))[2L]
    Lambda_star <- min(Lams_cand1, Lams_cand2)
    actual_ops <- actual_ops + 2

    # e_{gamma, j} = 3 * omega * a_c * rho * bar{pi}_j / r_j (j >= 2), e_{gamma, 1} = 0
    e_gamma <- numeric(m)
    e_r <- numeric(m)
    e_lambda_beta_coord <- numeric(p)
    names(e_gamma) <- names(e_r) <- cohort_labels
    names(e_lambda_beta_coord) <- z_colnames

    three_pt_bound <- logit_b64_mul(gmm_iv(3.0), gmm_iv(point_root_bound))[2L]
    actual_ops <- actual_ops + 1
    for (j in seq_len(m)) {
      if (j == 1L) {
        e_gamma[j] <- 0.0
      } else {
        num_eg <- logit_b64_mul(gmm_iv(three_pt_bound), gmm_iv(nuis$pi_bar_cohort[j]))[2L]
        e_gamma[j] <- logit_b64_div(gmm_iv(num_eg), gmm_iv(r_cohort[j, 1L]))[2L]
        actual_ops <- actual_ops + 2
      }
      t1_er <- logit_b64_mul(gmm_iv(nuis$e_L_upper[j]), gmm_iv(Lambda_star))[2L]
      t2_er <- logit_b64_mul(gmm_iv(nuis$ell_0_upper[j]), gmm_iv(e_lambda))[2L]
      e_r[j] <- logit_b64_add(logit_b64_add(gmm_iv(e_gamma[j]), gmm_iv(t1_er)), gmm_iv(t2_er))[2L]
      actual_ops <- actual_ops + 4
    }
    for (k in seq_len(p)) {
      e_lambda_beta_coord[k] <- logit_b64_mul(gmm_iv(coord_dual_iv[k, 2L]), gmm_iv(e_lambda))[2L]
      actual_ops <- actual_ops + 1
    }

    list(
      T_iv = T_iv,
      C_0_iv = C_0_iv,
      C_star_direct_iv = C_star_direct_iv,
      C_star_ball_iv = C_star_ball_iv,
      C_star_iv = C_star_iv,
      mu_0_iv = mu_0_iv,
      mu_star_iv = mu_star_iv,
      omega_iv = omega_iv,
      v_0_iv = v_0_iv,
      v_star_iv = v_star_iv,
      v_0_dual_upper = dn_v0$upper,
      v_star_dual_upper = dn_vs$upper,
      G_beta_0_iv = G_beta_0_iv,
      G_beta_star_iv = G_beta_star_iv,
      gamma_0_iv = gamma_0_iv,
      gamma_star_iv = gamma_star_iv,
      lam_sol_0 = lam_sol_0,
      lam_sol_star = lam_sol_star,
      lambda_r_0_iv = lambda_r_0_iv,
      lambda_r_star_iv = lambda_r_star_iv,
      rate_coef_s0_0_iv = rate_coef_s0_0_iv,
      rate_coef_s1_0_iv = rate_coef_s1_0_iv,
      rate_coef_s0_star_iv = rate_coef_s0_star_iv,
      rate_coef_s1_star_iv = rate_coef_s1_star_iv,
      root_bounds = list(
        omega = omega_hi,
        a = a_hi,
        a_c = a_c_hi,
        rho = rho_hi,
        x_bar = x_bar_hi,
        x_bar_c = x_bar_c_hi,
        u1_x_bar = nuis$u1_x_bar,
        u2_x_bar_c = nuis$u2_x_bar_c,
        point_mu = point_root_bound,
        e_C = point_root_bound,
        e_T = 0.0,
        e_v = e_v,
        Lambda_0 = Lambda_0,
        e_lambda = e_lambda,
        e_lambda_beta_coord = e_lambda_beta_coord,
        Lambda_star = Lambda_star,
        e_gamma = e_gamma,
        e_L = nuis$e_L_upper,
        ell_0 = nuis$ell_0_upper,
        ell_star = nuis$ell_star_upper,
        e_r = e_r
      )
    )
  }, cdid_ipw_enclosure_arithmetic_inconclusive = function(e) e,
     cdid_logit_arithmetic_inconclusive = function(e) e,
     cdid_gmm_accuracy_error = function(e) e)

  if (inherits(link_try, "error")) {
    reason <- if (!is.null(link_try$reason) && nzchar(link_try$reason)) link_try$reason else "link_arithmetic_refusal"
    inc_fields <- list(
      certified = FALSE,
      status = "INCONCLUSIVE",
      inconclusive = TRUE,
      inconclusive_type = "arithmetic",
      inconclusive_reasons = reason,
      separation_certified = FALSE,
      sampling_certified = FALSE,
      public_gate_adopted = FALSE,
      storage_mode = if (is_compact) "compact_view" else NULL,
      link_id = link_id,
      functional_key = functional_key,
      link_ref = if (is_compact) ctx$link_ref else NULL,
      nuisance_ref = if (is_compact) problem$nuisance_ref else NULL,
      frame_ref = if (is_compact) problem$frame_ref else NULL,
      panel_ref = if (is_compact) problem$panel_ref else NULL,
      problem_ref = problem$ref,
      endpoints = endpoints_list,
      loading_beta_candidate = tilde_lambda,
      candidates = cand_norm,
      enclosures_star = NULL,
      enclosures_0 = NULL,
      root_perturbation_bounds = NULL,
      separated_error_budget = NULL,
      arithmetic_refusal = list(
        occurred = TRUE,
        reason = reason,
        message = conditionMessage(link_try),
        condition_class = class(link_try),
        diagnostics = link_try$diagnostics
      ),
      control = ctrl,
      reader_control = if (is_compact) ctx$reader_ctrl else NULL,
      preflight = preflight,
      read_diagnostics = if (is_compact) list(
        read_outcomes = TRUE,
        endpoints_read = logit_cert_count_output(total_endpoints_read),
        n_endpoints_read = logit_cert_count_output(total_endpoints_read),
        reader_passes_completed = reader_passes_completed,
        actual_reader_work_ops = total_reader_ops,
        risk_filtered_before_outcome_read = TRUE,
        off_risk_endpoints_read = FALSE,
        nonobserver_endpoints_read = FALSE,
        cached_rounded_dy_returned = FALSE
      ) else NULL,
      work = if (is_compact) {
        list(
          validation_chunks_processed = logit_cert_count_output(validation_chunks),
          arithmetic_passes_completed = arith_passes,
          arithmetic_chunks_processed = logit_cert_count_output(arith_chunks),
          reader_passes_completed = reader_passes_completed,
          triangular_solves = logit_cert_count_output(tri_solves),
          exp_endpoint_evaluations = logit_cert_count_output(exp_evals),
          actual_interval_ops = actual_ops,
          actual_reader_work_ops = total_reader_ops,
          combined_work_ops = actual_ops + total_reader_ops,
          endpoints_read = logit_cert_count_output(total_endpoints_read)
        )
      } else {
        list(
          validation_chunks_processed = logit_cert_count_output(validation_chunks),
          arithmetic_passes_completed = arith_passes,
          arithmetic_chunks_processed = logit_cert_count_output(arith_chunks),
          triangular_solves = logit_cert_count_output(tri_solves),
          actual_interval_ops = actual_ops
        )
      }
    )
    if (!is_compact) {
      inc_fields$storage_mode <- NULL
      inc_fields$link_ref <- NULL
      inc_fields$nuisance_ref <- NULL
      inc_fields$frame_ref <- NULL
      inc_fields$panel_ref <- NULL
      inc_fields$reader_control <- NULL
      inc_fields$read_diagnostics <- NULL
    }
    inc_obj <- structure(inc_fields, class = "cdid_ipw_link_enclosure")
    return(ipw_enclosure_seal(inc_obj, "link_enclosure"))
  }

  # 6. Corollary sep separated evaluation-error and root-error budget for stored candidates
  rb <- link_try$root_bounds
  sep_try <- ipw_enclosure_capture_arithmetic({
    sep_dist <- function(iv, x) { actual_ops <<- actual_ops + 2; ipw_enclosure_iv_dist(iv, x) }
    sep_add <- function(a, b) { actual_ops <<- actual_ops + 1; logit_b64_add(a, b) }
    sep_budget <- list(
      loading_beta = list(
        candidate = tilde_lambda,
        eval_error_coord = vapply(seq_len(p), function(k) sep_dist(link_try$lam_sol_0$coord_iv[k, ], tilde_lambda[k]), double(1L)),
        root_error_coord = rb$e_lambda_beta_coord,
        total_error_coord = vapply(seq_len(p), function(k) {
          ev <- sep_dist(link_try$lam_sol_0$coord_iv[k, ], tilde_lambda[k])
          sep_add(gmm_iv(ev), gmm_iv(rb$e_lambda_beta_coord[k]))[2L]
        }, double(1L)),
        eval_error_m_norm = link_try$lam_sol_0$radius_m_norm,
        root_error_m_norm = rb$e_lambda,
        total_error_m_norm = min(
          link_try$lam_sol_star$radius_m_norm,
          sep_add(gmm_iv(link_try$lam_sol_0$radius_m_norm), gmm_iv(rb$e_lambda))[2L]
        ),
        direct_star_error_coord = vapply(seq_len(p), function(k) sep_dist(link_try$lam_sol_star$coord_iv[k, ], tilde_lambda[k]), double(1L))
      )
    )
    if (!is.null(cand_norm)) {
      if (!is.null(cand_norm$mu)) {
        ev_mu <- sep_dist(link_try$mu_0_iv, cand_norm$mu)
        tot_mu <- sep_add(gmm_iv(ev_mu), gmm_iv(rb$point_mu))[2L]
        sep_budget$mu <- list(
          candidate = cand_norm$mu,
          eval_error = ev_mu,
          root_error = rb$point_mu,
          total_error_bound = tot_mu,
          direct_star_error = sep_dist(link_try$mu_star_iv, cand_norm$mu)
        )
      }
      if (!is.null(cand_norm$treated_mean)) {
        ev_T <- sep_dist(link_try$T_iv, cand_norm$treated_mean)
        sep_budget$treated_mean <- list(
          candidate = cand_norm$treated_mean,
          eval_error = ev_T,
          root_error = 0.0,
          total_error_bound = ev_T
        )
      }
      if (!is.null(cand_norm$control_mean)) {
        ev_C <- sep_dist(link_try$C_0_iv, cand_norm$control_mean)
        tot_C <- sep_add(gmm_iv(ev_C), gmm_iv(rb$e_C))[2L]
        sep_budget$control_mean <- list(
          candidate = cand_norm$control_mean,
          eval_error = ev_C,
          root_error = rb$e_C,
          total_error_bound = tot_C,
          direct_star_error = sep_dist(link_try$C_star_iv, cand_norm$control_mean)
        )
      }
      v_cand <- if (!is.null(cand_norm$v)) cand_norm$v else if (!is.null(cand_norm$G_beta)) -cand_norm$G_beta else NULL
      if (!is.null(v_cand)) {
        ev_v_coord <- vapply(seq_len(p), function(k) sep_dist(link_try$v_0_iv[k, ], v_cand[k]), double(1L))
        sep_budget$v <- list(
          candidate = v_cand,
          G_beta_candidate = -v_cand,
          eval_error_coord = ev_v_coord,
          root_error_dual_norm = rb$e_v,
          direct_star_error_coord = vapply(seq_len(p), function(k) sep_dist(link_try$v_star_iv[k, ], v_cand[k]), double(1L))
        )
      }
      gam_cand <- if (!is.null(cand_norm$gamma)) cand_norm$gamma else cand_norm$G_r
      if (!is.null(gam_cand)) {
        ev_g <- vapply(seq_len(m), function(j) sep_dist(link_try$gamma_0_iv[j, ], gam_cand[j]), double(1L))
        tot_g <- vapply(seq_len(m), function(j) sep_add(gmm_iv(ev_g[j]), gmm_iv(rb$e_gamma[j]))[2L], double(1L))
        sep_budget$gamma <- list(
          candidate = gam_cand,
          eval_error = ev_g,
          root_error = rb$e_gamma,
          total_error_bound = tot_g,
          direct_star_error = vapply(seq_len(m), function(j) sep_dist(link_try$gamma_star_iv[j, ], gam_cand[j]), double(1L))
        )
      }
      if (!is.null(cand_norm$loading_rate)) {
        lr_cand <- cand_norm$loading_rate
        ev_lr <- vapply(seq_len(m), function(j) sep_dist(link_try$lambda_r_0_iv[j, ], lr_cand[j]), double(1L))
        tot_lr <- vapply(seq_len(m), function(j) sep_add(gmm_iv(ev_lr[j]), gmm_iv(rb$e_r[j]))[2L], double(1L))
        sep_budget$loading_rate <- list(
          candidate = lr_cand,
          eval_error = ev_lr,
          root_error = rb$e_r,
          total_error_bound = tot_lr,
          direct_star_error = vapply(seq_len(m), function(j) sep_dist(link_try$lambda_r_star_iv[j, ], lr_cand[j]), double(1L))
        )
      }
      L_cand <- if (!is.null(cand_norm$L)) cand_norm$L else cand_norm$A_beta_r
      if (!is.null(L_cand)) {
        ev_L <- matrix(0.0, nrow = p, ncol = m, dimnames = list(z_colnames, cohort_labels))
        st_L <- matrix(0.0, nrow = p, ncol = m, dimnames = list(z_colnames, cohort_labels))
        for (j in seq_len(m)) {
          for (k in seq_len(p)) {
            ev_L[k, j] <- sep_dist(nuis$L_0[k, j, ], L_cand[k, j])
            st_L[k, j] <- sep_dist(nuis$L_star[k, j, ], L_cand[k, j])
          }
        }
        sep_budget$L <- list(
          candidate = L_cand,
          eval_error_entrywise = ev_L,
          root_error_dual_norm_by_cohort = rb$e_L,
          direct_star_error_entrywise = st_L
        )
      }
    }

    sep_budget
  })
  if (inherits(sep_try, "error")) {
    return(ipw_enclosure_inconclusive_result(
      sep_try, "cdid_ipw_link_enclosure",
      list(
        storage_mode = if (is_compact) "compact_view" else NULL,
        link_ref = if (is_compact) ctx$link_ref else NULL,
        nuisance_ref = if (is_compact) problem$nuisance_ref else NULL,
        frame_ref = if (is_compact) problem$frame_ref else NULL,
        panel_ref = if (is_compact) problem$panel_ref else NULL,
        problem_ref = problem$ref,
        preflight = preflight,
        control = ctrl,
        work = list(actual_interval_ops = actual_ops)
      ),
      "link_enclosure"
    ))
  }
  sep_budget <- sep_try
  stopifnot(actual_ops <= preflight$total_work_ops_bound)
  if (is_compact) {
    stopifnot(total_reader_ops <= preflight$total_reader_work_ops_bound)
    stopifnot((actual_ops + total_reader_ops) <= preflight$combined_work_ops_bound)
  }

  out_fields <- list(
    certified = TRUE,
    status = "CERTIFIED",
    inconclusive = FALSE,
    inconclusive_type = "none",
    inconclusive_reasons = character(0L),
    separation_certified = FALSE,
    sampling_certified = FALSE,
    public_gate_adopted = FALSE,
    storage_mode = if (is_compact) "compact_view" else NULL,
    link_id = link_id,
    functional_key = functional_key,
    link_ref = if (is_compact) ctx$link_ref else NULL,
    nuisance_ref = if (is_compact) problem$nuisance_ref else NULL,
    frame_ref = if (is_compact) problem$frame_ref else NULL,
    panel_ref = if (is_compact) problem$panel_ref else NULL,
    problem_ref = problem$ref,
    endpoints = endpoints_list,
    loading_beta_candidate = tilde_lambda,
    candidates = cand_norm,
    enclosures_star = list(
      T = link_try$T_iv,
      C = link_try$C_star_iv,
      C_direct = link_try$C_star_direct_iv,
      C_ball = link_try$C_star_ball_iv,
      mu = link_try$mu_star_iv,
      Sigma_plus = nuis$Sigma_plus,
      Sigma_minus = nuis$Sigma_minus_star,
      d_plus = nuis$d_plus,
      d_minus = nuis$d_minus_star,
      omega = link_try$omega_iv,
      v = link_try$v_star_iv,
      v_dual_upper = link_try$v_star_dual_upper,
      G_beta = link_try$G_beta_star_iv,
      gamma = link_try$gamma_star_iv,
      G_r = link_try$gamma_star_iv,
      H = nuis$H_star,
      L = nuis$L_star,
      A_beta_r = nuis$L_star,
      lambda_beta = link_try$lam_sol_star$coord_iv,
      loading_beta = link_try$lam_sol_star$coord_iv,
      lambda_beta_ball = link_try$lam_sol_star,
      lambda_r = link_try$lambda_r_star_iv,
      loading_rate = link_try$lambda_r_star_iv,
      rate_coef_s0 = link_try$rate_coef_s0_star_iv,
      rate_coef_s1 = link_try$rate_coef_s1_star_iv
    ),
    enclosures_0 = list(
      T = link_try$T_iv,
      C = link_try$C_0_iv,
      mu = link_try$mu_0_iv,
      Sigma_plus = nuis$Sigma_plus,
      Sigma_minus = nuis$Sigma_minus_0,
      d_plus = nuis$d_plus,
      d_minus = nuis$d_minus_0,
      omega = link_try$omega_iv,
      v = link_try$v_0_iv,
      v_dual_upper = link_try$v_0_dual_upper,
      G_beta = link_try$G_beta_0_iv,
      gamma = link_try$gamma_0_iv,
      G_r = link_try$gamma_0_iv,
      H = nuis$H_0,
      L = nuis$L_0,
      A_beta_r = nuis$L_0,
      lambda_beta = link_try$lam_sol_0$coord_iv,
      loading_beta = link_try$lam_sol_0$coord_iv,
      lambda_beta_ball = link_try$lam_sol_0,
      lambda_r = link_try$lambda_r_0_iv,
      loading_rate = link_try$lambda_r_0_iv,
      rate_coef_s0 = link_try$rate_coef_s0_0_iv,
      rate_coef_s1 = link_try$rate_coef_s1_0_iv
    ),
    root_perturbation_bounds = rb,
    separated_error_budget = sep_budget,
    separated_errors = sep_budget,
    arithmetic_refusal = list(
      occurred = FALSE,
      reason = "none",
      message = NULL,
      condition_class = NULL,
      diagnostics = NULL
    ),
    control = ctrl,
    reader_control = if (is_compact) ctx$reader_ctrl else NULL,
    preflight = preflight,
    read_diagnostics = if (is_compact) list(
      read_outcomes = TRUE,
      endpoints_read = logit_cert_count_output(total_endpoints_read),
      n_endpoints_read = logit_cert_count_output(total_endpoints_read),
      reader_passes_completed = reader_passes_completed,
      actual_reader_work_ops = total_reader_ops,
      risk_filtered_before_outcome_read = TRUE,
      off_risk_endpoints_read = FALSE,
      nonobserver_endpoints_read = FALSE,
      cached_rounded_dy_returned = FALSE
    ) else NULL,
    work = if (is_compact) {
      list(
        validation_chunks_processed = logit_cert_count_output(validation_chunks),
        arithmetic_passes_completed = arith_passes,
        arithmetic_chunks_processed = logit_cert_count_output(arith_chunks),
        reader_passes_completed = reader_passes_completed,
        triangular_solves = logit_cert_count_output(tri_solves),
        exp_endpoint_evaluations = logit_cert_count_output(exp_evals),
        actual_interval_ops = actual_ops,
        actual_reader_work_ops = total_reader_ops,
        combined_work_ops = actual_ops + total_reader_ops,
        endpoints_read = logit_cert_count_output(total_endpoints_read),
        interval_ops_definition = "Conservative interval-operation charges at completed helper boundaries; excludes CPU instructions, R object headers, and process peak RSS."
      )
    } else {
      list(
        validation_chunks_processed = logit_cert_count_output(validation_chunks),
        arithmetic_passes_completed = arith_passes,
        arithmetic_chunks_processed = logit_cert_count_output(arith_chunks),
        triangular_solves = logit_cert_count_output(tri_solves),
        actual_interval_ops = actual_ops,
        interval_ops_definition = "Conservative interval-operation charges at completed helper boundaries; excludes CPU instructions, R object headers, and process peak RSS."
      )
    },
    provenance = if (is_compact) {
      list(
        storage_mode = "compact_view",
        normalizer_n = as.double(n),
        cohorts = raw$cohorts,
        cohort_labels = cohort_labels,
        defining_pair = ctx$defining_pair,
        differences_source = "exact_real_subtraction_of_stored_y_right_minus_y_left_enclosed_outward",
        off_observed_endpoints_read = FALSE,
        per_link_endpoint_copies_retained = FALSE,
        full_n_position_map_retained = FALSE,
        shared_nuisance_reuse_supported = TRUE,
        coherent_shared_error_rule = "Links sharing one certified nuisance problem share beta* - beta0; linear combinations bound root error by sum_l |kappa_l| * e_root_l (never independent root-sum-square).",
        covariance_or_gmm_computed = FALSE
      )
    } else {
      list(
        normalizer_n = as.double(n),
        cohorts = raw$cohorts,
        cohort_labels = cohort_labels,
        differences_source = "exact_real_subtraction_of_stored_y_right_minus_y_left_enclosed_outward",
        off_observed_endpoints_read = FALSE,
        shared_nuisance_reuse_supported = TRUE,
        coherent_shared_error_rule = "Links sharing one certified nuisance problem share beta* - beta0; linear combinations bound root error by sum_l |kappa_l| * e_root_l (never independent root-sum-square).",
        covariance_or_gmm_computed = FALSE
      )
    }
  )
  if (!is_compact) {
    out_fields$storage_mode <- NULL
    out_fields$link_ref <- NULL
    out_fields$nuisance_ref <- NULL
    out_fields$frame_ref <- NULL
    out_fields$panel_ref <- NULL
    out_fields$reader_control <- NULL
    out_fields$read_diagnostics <- NULL
  }

  link_obj <- structure(out_fields, class = "cdid_ipw_link_enclosure")
  ipw_enclosure_seal(link_obj, "link_enclosure")
}

enclose_ipw_link <- function(
  problem,
  id,
  y_left,
  y_right,
  loading_beta_candidate,
  control,
  pair_observed = NULL,
  risk_mask = NULL,
  observed_mask = NULL,
  cohorts = NULL,
  candidates = NULL,
  link_id = NULL,
  functional_key = NULL
) {
  ipw_enclosure_require_deps()
  if (missing(problem) || missing(id) || missing(y_left) || missing(y_right) ||
      missing(loading_beta_candidate) || missing(control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`problem`, `id`, `y_left`, `y_right`, `loading_beta_candidate`, `control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }

  # 1. Validate sealed problem, control, and lightweight headers BEFORE scanning endpoints
  ipw_validate_nuisance_problem_mode(problem, require_certified = FALSE, expected_storage_mode = "dense")
  ctrl <- validate_ipw_link_enclosure_control(control)

  raw <- problem$raw
  n <- nrow(raw$Z)
  p <- ncol(raw$Z)
  m <- length(raw$cohorts)

  if (!is.null(link_id) && !gmm_scalar_character(link_id)) {
    cdid_ipw_enclosure_abort("`link_id` must be NULL or a nonempty character scalar.", "cdid_ipw_enclosure_input_error")
  }
  if (!is.null(functional_key) && !gmm_scalar_character(functional_key)) {
    cdid_ipw_enclosure_abort("`functional_key` must be NULL or a nonempty character scalar.", "cdid_ipw_enclosure_input_error")
  }
  if (!is.null(cohorts)) {
    if (!is.atomic(cohorts) || !is.null(dim(cohorts)) || is.object(cohorts) ||
        !is.numeric(cohorts) || length(cohorts) != m) {
      cdid_ipw_enclosure_abort(
        "Link `cohorts` header disagrees with the certified nuisance problem's ordered cohorts.",
        c("cdid_ipw_enclosure_cohort_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }

  # Header-only checks on id, y_left, y_right, loading_beta_candidate, candidates
  if (is.null(id) || is.factor(id) || is.object(id) || !is.atomic(id) ||
      !is.null(dim(id)) || length(id) != n) {
    cdid_ipw_enclosure_abort(
      "`id` on link call must be a plain atomic vector of length `n_a_units`.",
      c("cdid_ipw_enclosure_id_error", "cdid_ipw_enclosure_input_error")
    )
  }
  for (ep_item in list(list(v = y_left, nm = "y_left"), list(v = y_right, nm = "y_right"))) {
    v <- ep_item$v; nm <- ep_item$nm
    if (is.null(v) || is.object(v) || !is.atomic(v) || !is.null(dim(v)) ||
        !is.numeric(v) || length(v) != n) {
      cdid_ipw_enclosure_abort(
        paste0("`", nm, "` must be a plain atomic numeric vector of length `n_a_units` (", n, ") in all-A order."),
        c("cdid_ipw_enclosure_endpoint_error", "cdid_ipw_enclosure_input_error")
      )
    }
    ax <- attributes(v)
    if (!is.null(ax) && any(names(ax) %in% c("dy", "pre_subtracted", "certified"))) {
      cdid_ipw_enclosure_abort(
        paste0("Decorative or cached-difference attributes on `", nm, "` are refused."),
        c("cdid_ipw_enclosure_endpoint_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (is.null(loading_beta_candidate) || is.object(loading_beta_candidate) ||
      !is.atomic(loading_beta_candidate) || !is.null(dim(loading_beta_candidate)) ||
      !is.numeric(loading_beta_candidate) || length(loading_beta_candidate) != p) {
    cdid_ipw_enclosure_abort(
      paste0("`loading_beta_candidate` must be a finite atomic numeric vector of length p = ", p, "."),
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.null(candidates) && (!is.list(candidates) || is.object(candidates))) {
    cdid_ipw_enclosure_abort(
      "`candidates` must be NULL or a plain named list of candidate link scalars/loadings.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }

  # 2. Preflight BEFORE scanning content or allocating retained endpoint buffers
  preflight <- preflight_ipw_link_enclosure(
    n = n, p = p, m = m,
    control = ctrl,
    certificate_control = problem$controls$certificate_control,
    has_candidates = !is.null(candidates)
  )
  step <- preflight$effective_chunk_rows

  if (!is.null(cohorts)) {
    coh_arg <- validate_logit_certificate_cohorts(cohorts)
    if (!identical(coh_arg, raw$cohorts)) {
      cdid_ipw_enclosure_abort(
        "Link `cohorts` disagrees with the certified nuisance problem's ordered cohorts.",
        c("cdid_ipw_enclosure_cohort_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (anyNA(loading_beta_candidate) || any(!is.finite(loading_beta_candidate))) {
    cdid_ipw_enclosure_abort(
      paste0("`loading_beta_candidate` must be a finite atomic numeric vector of length p = ", p, "."),
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  tilde_lambda <- as.double(loading_beta_candidate)

  cand_norm <- ipw_enclosure_validate_candidates(candidates, p = p, m = m)
  if (!is.null(cand_norm$loading_beta) && !identical(cand_norm$loading_beta, tilde_lambda)) {
    cdid_ipw_enclosure_abort(
      "`candidates$loading_beta` disagrees with `loading_beta_candidate`.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }

  # 3. Validate IDs and optional masks in chunks
  id_info <- ipw_enclosure_validate_id(id, n = n, step = step, label = "id")
  if (!identical(id_info$id_type, raw$id_type) || !identical(id_info$id_key, raw$id_key)) {
    cdid_ipw_enclosure_abort(
      "Link `id` vector differs in value, type, or unit order from the certified nuisance problem.",
      c("cdid_ipw_enclosure_id_error", "cdid_ipw_enclosure_input_error")
    )
  }
  ipw_enclosure_check_optional_masks(
    n = n,
    in_risk = problem$masks$in_risk,
    in_observed = problem$masks$in_observed,
    cohort_index = problem$masks$cohort_index,
    pair_observed_ref = raw$pair_observed,
    risk_mask = risk_mask,
    observed_mask = observed_mask,
    pair_observed_arg = pair_observed
  )

  # 4. Validate finite endpoints ONLY on O = R intersect S (never read y_left/y_right outside O!)
  obs_idx <- problem$masks$observed_indices
  n_obs <- length(obs_idx)
  if (n_obs < 2L || problem$dimensions$n_target_observed < 1L || problem$dimensions$n_control_observed < 1L) {
    cdid_ipw_enclosure_abort(
      "Link requires nonempty observed target and control risk-set support.",
      c("cdid_ipw_enclosure_support_error", "cdid_ipw_enclosure_endpoint_error")
    )
  }

  y_left_obs <- numeric(n_obs)
  y_right_obs <- numeric(n_obs)
  validation_chunks <- 0
  a_start <- 1L
  while (a_start <= n) {
    ipw_enclosure_check_interrupt()
    a_end <- min(as.double(n), as.double(a_start) + step - 1)
    # Find positions in obs_idx that fall in [a_start, a_end]
    in_ch <- which(obs_idx >= a_start & obs_idx <= a_end)
    if (length(in_ch) > 0L) {
      u_idx <- obs_idx[in_ch]
      yl_sub <- y_left[u_idx]
      yr_sub <- y_right[u_idx]
      if (anyNA(yl_sub) || any(!is.finite(yl_sub)) || anyNA(yr_sub) || any(!is.finite(yr_sub))) {
        cdid_ipw_enclosure_abort(
          "Required observed risk-set endpoints (`y_left` or `y_right` on O = R intersect S) contain NA, NaN, or Inf.",
          c("cdid_ipw_enclosure_endpoint_error", "cdid_ipw_enclosure_input_error"),
          diagnostics = list(chunk_start = a_start, chunk_end = a_end)
        )
      }
      y_left_obs[in_ch] <- as.double(yl_sub)
      y_right_obs[in_ch] <- as.double(yr_sub)
    }
    validation_chunks <- validation_chunks + 1
    a_start <- as.double(a_end) + 1
  }

  # Map from all-A index 1..n to position in 1..n_obs (0L if not in O)
  a_to_obs_pos <- integer(n)
  a_to_obs_pos[obs_idx] <- seq_len(n_obs)

  ctx <- construct_dense_ipw_link_access_context(
    problem = problem,
    ctrl = ctrl,
    preflight = preflight,
    step = step,
    tilde_lambda = tilde_lambda,
    cand_norm = cand_norm,
    link_id = link_id,
    functional_key = functional_key,
    obs_idx = obs_idx,
    a_to_obs_pos = a_to_obs_pos,
    y_left_obs = y_left_obs,
    y_right_obs = y_right_obs,
    validation_chunks = validation_chunks
  )
  ipw_enclosure_link_run_with_context(ctx)
}

enclose_ipw_link_view <- function(
  problem,
  link_view,
  loading_beta_candidate,
  control,
  reader_control,
  candidates = NULL
) {
  ipw_enclosure_require_deps()
  if (missing(problem) || missing(link_view) || missing(loading_beta_candidate) ||
      missing(control) || missing(reader_control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`problem`, `link_view`, `loading_beta_candidate`, `control`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }

  # 1. Validate sealed compact problem, link_view, and controls BEFORE any chunk scan
  ipw_validate_nuisance_problem_mode(problem, require_certified = FALSE, expected_storage_mode = "compact_view")
  ctrl <- validate_ipw_link_enclosure_control(control)
  expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  trace_env <- attr(reader_control, "trace_env", exact = TRUE)
  reader_ctrl <- validate_ipw_owner_control(reader_control)

  frame_owner <- get("frame_owner", envir = problem$frame_ref, inherits = FALSE)
  nuisance_view <- get("nuisance_view", envir = problem$nuisance_ref, inherits = FALSE)
  validate_ipw_link_view(link_view, frame_owner = frame_owner, nuisance_view = nuisance_view)

  raw <- problem$raw
  n <- as.integer(raw$n_a_units)
  p <- as.integer(raw$p_covariates)
  m <- length(raw$cohorts)

  if (!identical(link_view$panel_ref, problem$panel_ref) ||
      !identical(link_view$cohorts$cohorts, raw$cohorts) ||
      !identical(link_view$cohorts$g, nuisance_view$cohorts$g) ||
      !identical(link_view$cohorts$K, nuisance_view$cohorts$K) ||
      !isTRUE(link_view$nuisance_mask_verified)) {
    cdid_ipw_enclosure_abort(
      "Link view disagrees with the certified nuisance problem's frame, ordered cohorts, or verified observer mask.",
      c("cdid_ipw_enclosure_link_error", "cdid_ipw_enclosure_integrity_error", "cdid_ipw_shared_incompatible_error")
    )
  }

  ipw_owner_refuse_decorative_attrs(loading_beta_candidate, "loading_beta_candidate")
  ipw_owner_refuse_decorative_attrs(candidates, "candidates")
  if (is.null(loading_beta_candidate) || is.object(loading_beta_candidate) ||
      !is.atomic(loading_beta_candidate) || !is.null(dim(loading_beta_candidate)) ||
      !is.numeric(loading_beta_candidate) || length(loading_beta_candidate) != p) {
    cdid_ipw_enclosure_abort(
      paste0("`loading_beta_candidate` must be a finite atomic numeric vector of length p = ", p, "."),
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.null(candidates) && (!is.list(candidates) || is.object(candidates))) {
    cdid_ipw_enclosure_abort(
      "`candidates` must be NULL or a plain named list of candidate link scalars/loadings.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }

  # 2. Combined arithmetic + 2-pass reader preflight BEFORE scanning candidate contents or chunks
  preflight <- preflight_ipw_link_enclosure_view(
    problem = problem,
    link_view = link_view,
    control = ctrl,
    reader_control = reader_ctrl,
    has_candidates = !is.null(candidates),
    expected_a_id = expected_a_id
  )
  step <- preflight$effective_chunk_rows

  if (anyNA(loading_beta_candidate) || any(!is.finite(loading_beta_candidate))) {
    cdid_ipw_enclosure_abort(
      paste0("`loading_beta_candidate` must be a finite atomic numeric vector of length p = ", p, "."),
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  tilde_lambda <- as.double(loading_beta_candidate)

  cand_norm <- ipw_enclosure_validate_candidates(candidates, p = p, m = m)
  if (!is.null(cand_norm$loading_beta) && !identical(cand_norm$loading_beta, tilde_lambda)) {
    cdid_ipw_enclosure_abort(
      "`candidates$loading_beta` disagrees with `loading_beta_candidate`.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }

  # Optional candidates already stored on `link_view` must agree with explicit inputs;
  # none may supply a hidden default.
  if (!is.null(link_view$loading_beta_candidate) &&
      !identical(link_view$loading_beta_candidate, tilde_lambda)) {
    cdid_ipw_enclosure_abort(
      "Explicit `loading_beta_candidate` disagrees with the candidate stored on `link_view`.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.null(link_view$candidates) &&
      (is.null(cand_norm) || !identical(link_view$candidates, cand_norm))) {
    cdid_ipw_enclosure_abort(
      "Explicit `candidates` argument is missing or disagrees with the candidates stored on `link_view`.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }

  # Require actual support prerequisites even if `link_view` was constructed with `require_support = FALSE`
  sup <- link_view$support
  if (problem$dimensions$n_observed_risk_units < 2L ||
      problem$dimensions$n_target_observed < 1L ||
      problem$dimensions$n_control_observed < 1L ||
      sup$n_observed_risk_units < 2L ||
      sup$by_cohort_pair_count[1L] < 1L ||
      sum(sup$by_cohort_pair_count[-1L]) < 1L ||
      any(sup$by_cohort_a_count < 1L) ||
      any(sup$by_cohort_pair_count < 1L) ||
      any(sup$by_cohort_a_mass <= 0.0) ||
      any(sup$by_cohort_pair_mass <= 0.0)) {
    cdid_ipw_enclosure_abort(
      "Link requires nonempty observed target and control risk-set support across all cohorts.",
      c("cdid_ipw_enclosure_support_error", "cdid_ipw_enclosure_endpoint_error")
    )
  }

  ctx <- construct_view_ipw_link_access_context(
    problem = problem,
    link_view = link_view,
    ctrl = ctrl,
    reader_ctrl = reader_ctrl,
    preflight = preflight,
    step = step,
    tilde_lambda = tilde_lambda,
    cand_norm = cand_norm,
    expected_a_id = expected_a_id,
    trace_env = trace_env
  )
  ipw_enclosure_link_run_with_context(ctx)
}

construct_ipw_link_enclosure <- enclose_ipw_link
certify_ipw_link_enclosure <- enclose_ipw_link
construct_ipw_link_enclosure_view <- enclose_ipw_link_view
certify_ipw_link_enclosure_view <- enclose_ipw_link_view

ipw_validate_link_enclosure_mode <- function(
  link_enclosure,
  problem = NULL,
  require_certified = FALSE,
  expected_storage_mode = NULL
) {
  ipw_enclosure_require_deps()
  if (missing(link_enclosure) || is.null(link_enclosure) ||
      !inherits(link_enclosure, "cdid_ipw_link_enclosure") ||
      !is.list(link_enclosure) || !ipw_enclosure_seal_ok(link_enclosure, "link_enclosure")) {
    cdid_ipw_enclosure_abort(
      "Link enclosure object is malformed, foreign, or differs from its sealed snapshot.",
      c("cdid_ipw_enclosure_link_error", "cdid_ipw_enclosure_integrity_error")
    )
  }
  actual_mode <- if (identical(link_enclosure$storage_mode, "compact_view")) "compact_view" else "dense"
  if (!is.null(expected_storage_mode) && !identical(actual_mode, expected_storage_mode)) {
    cdid_ipw_enclosure_abort(
      paste0("Expected link enclosure with storage mode `", expected_storage_mode, "`, but got `", actual_mode, "`."),
      c("cdid_ipw_enclosure_link_error", "cdid_ipw_enclosure_integrity_error", "cdid_ipw_enclosure_input_error")
    )
  }
  pref <- link_enclosure$problem_ref
  if (!is.environment(pref) || !environmentIsLocked(pref) ||
      !exists("nuisance_problem", envir = pref, inherits = FALSE)) {
    cdid_ipw_enclosure_abort(
      "Link enclosure lacks a valid sealed nuisance problem reference.",
      c("cdid_ipw_enclosure_link_error", "cdid_ipw_enclosure_integrity_error")
    )
  }
  prob <- get("nuisance_problem", envir = pref, inherits = FALSE)
  ipw_validate_nuisance_problem_mode(prob, require_certified = require_certified, expected_storage_mode = actual_mode)
  if (!is.null(problem)) {
    ipw_validate_nuisance_problem_mode(problem, require_certified = require_certified, expected_storage_mode = actual_mode)
    if (!identical(problem$ref, pref)) {
      cdid_ipw_enclosure_abort(
        "Link enclosure belongs to a different certified nuisance problem.",
        c("cdid_ipw_enclosure_problem_mismatch_error", "cdid_ipw_enclosure_integrity_error")
      )
    }
  }
  if (identical(actual_mode, "compact_view")) {
    ep <- link_enclosure$endpoints
    if (!is.list(ep) || !is.null(ep$obs_a_index) || !is.null(ep$a_to_obs_pos) ||
        !is.null(ep$y_left_obs) || !is.null(ep$y_right_obs)) {
      cdid_ipw_enclosure_abort(
        "Compact link enclosure must not retain `obs_a_index`, `a_to_obs_pos`, `y_left_obs`, or `y_right_obs` arrays.",
        c("cdid_ipw_enclosure_link_error", "cdid_ipw_enclosure_integrity_error")
      )
    }
    lref <- link_enclosure$link_ref
    if (!is.environment(lref) || !environmentIsLocked(lref) ||
        !exists("link_view", envir = lref, inherits = FALSE) ||
        !identical(link_enclosure$frame_ref, prob$frame_ref) ||
        !identical(link_enclosure$panel_ref, prob$panel_ref) ||
        !identical(link_enclosure$nuisance_ref, prob$nuisance_ref)) {
      cdid_ipw_enclosure_abort(
        "Compact link enclosure lacks valid sealed `link_ref`, `frame_ref`, `panel_ref`, or `nuisance_ref`.",
        c("cdid_ipw_enclosure_link_error", "cdid_ipw_enclosure_integrity_error")
      )
    }
    link_view <- get("link_view", envir = lref, inherits = FALSE)
    frame_owner <- get("frame_owner", envir = prob$frame_ref, inherits = FALSE)
    nuisance_view <- get("nuisance_view", envir = prob$nuisance_ref, inherits = FALSE)
    validate_ipw_link_view(link_view, frame_owner = frame_owner, nuisance_view = nuisance_view)
    if (!identical(link_view$cohorts$cohorts, prob$raw$cohorts) ||
        !isTRUE(link_view$nuisance_mask_verified)) {
      cdid_ipw_enclosure_abort(
        "Compact link enclosure's bound `link_view` disagrees with problem cohorts or lacks verified observer mask.",
        c("cdid_ipw_enclosure_link_error", "cdid_ipw_enclosure_integrity_error")
      )
    }
  }
  if (isTRUE(require_certified) && !isTRUE(link_enclosure$certified)) {
    cdid_ipw_enclosure_abort(
      "Certified link enclosure is required, but link status is INCONCLUSIVE.",
      c("cdid_ipw_enclosure_inconclusive", "cdid_ipw_enclosure_link_error"),
      diagnostics = list(inconclusive_reasons = link_enclosure$inconclusive_reasons)
    )
  }
  invisible(prob)
}

# ---------------------------------------------------------------------------
# Component 3: Outcome-free shared (xi_rate, xi_beta) chunk primitive (Lemma struct(b))
#              and bounded row-chunk 4-term IF enclosure emitter
# ---------------------------------------------------------------------------

new_ipw_if_chunk_scratch <- function(problem, a_start, a_end, control) {
  ipw_enclosure_require_deps()
  ipw_validate_nuisance_problem_mode(problem, require_certified = FALSE)
  n <- if (identical(problem$storage_mode, "compact_view")) {
    as.integer(problem$raw$n_a_units)
  } else {
    nrow(problem$raw$Z)
  }
  s_idx <- as.integer(ipw_enclosure_int_bound(a_start, "a_start", min_val = 1))
  e_idx <- as.integer(ipw_enclosure_int_bound(a_end, "a_end", min_val = 1))
  if (s_idx > e_idx || e_idx > n) {
    cdid_ipw_enclosure_abort(
      paste0("Chunk range `[a_start, a_end]` must satisfy `1 <= a_start <= a_end <= ", n, "`."),
      "cdid_ipw_enclosure_bound_error"
    )
  }
  b_len <- e_idx - s_idx + 1L
  if (missing(control)) cdid_ipw_enclosure_abort("Scratch requires an explicit resource control.", "cdid_ipw_enclosure_input_error")
  ctrl <- validate_ipw_link_enclosure_control(control)
  if (b_len > ctrl$max_chunk_rows) cdid_ipw_enclosure_abort("Scratch range exceeds explicit chunk cap.", "cdid_ipw_enclosure_bound_error")
  p <- if (identical(problem$storage_mode, "compact_view")) problem$raw$p_covariates else ncol(problem$raw$Z)
  entries <- ipw_enclosure_safe_sum(c(80 * b_len, 2 * ipw_enclosure_safe_product(b_len, p)),
    min(ctrl$max_matrix_entries, ctrl$max_retained_entries), "sealed scratch payload")
  ipw_enclosure_safe_product(4, b_len, ctrl$max_work_ops, "scratch initialization work")
  env <- new.env(parent = emptyenv())
  env$seal_token <- .ipw_access_seal_token
  env$problem_ref <- problem$ref
  env$storage_mode <- if (identical(problem$storage_mode, "compact_view")) "compact_view" else "dense"
  env$a_start <- s_idx
  env$a_end <- e_idx
  env$b_len <- b_len
  # Own cache values once; public fields only read this private state.
  state <- new.env(parent = emptyenv())
  state$rate_done <- logical(b_len)
  state$rate <- vector("list", b_len)
  state$obs_done <- logical(b_len)
  state$obs <- vector("list", b_len)
  state$outcome_free_chunk <- NULL
  state$shared_xi_chunk <- NULL
  state$copy_work_ops <- 0
  lockEnvironment(state, bindings = FALSE)
  for (name in ls(state, all.names = TRUE)) {
    reader_ctx <- new.env(parent = baseenv())
    reader_ctx$name <- name
    reader_ctx$state <- state
    lockEnvironment(reader_ctx, bindings = TRUE)
    reader <- function(value) {
      if (!missing(value)) stop("Read-only scratch cache field.", call. = FALSE)
      get(name, envir = state, inherits = FALSE)
    }
    # Never retain the constructor frame (problem, control, or temporaries).
    environment(reader) <- reader_ctx
    makeActiveBinding(name, reader, env)
  }
  env <- structure(env, class = "cdid_ipw_if_chunk_scratch")
  attr(state, "scratch_ref") <- env
  lockEnvironment(env, bindings = TRUE)
  env
}

ipw_if_scratch_set <- function(scratch, name, value) {
  if (!is.environment(scratch) || !inherits(scratch, "cdid_ipw_if_chunk_scratch") ||
      !environmentIsLocked(scratch) || !gmm_scalar_character(name) ||
      !name %in% c("rate_done", "rate", "obs_done", "obs", "outcome_free_chunk", "shared_xi_chunk", "copy_work_ops") ||
      !exists(name, envir = scratch, inherits = FALSE) ||
      !bindingIsLocked(name, scratch) || !bindingIsActive(name, scratch) ||
      !exists("seal_token", envir = scratch, inherits = FALSE) ||
      !bindingIsLocked("seal_token", scratch) ||
      !identical(scratch$seal_token, .ipw_access_seal_token)) {
    cdid_ipw_enclosure_abort("Invalid scratch cache update.",
      c("cdid_ipw_enclosure_scratch_error", "cdid_ipw_enclosure_integrity_error"))
  }
  reader_ctx <- environment(activeBindingFunction(name, scratch))
  if (!is.environment(reader_ctx) || !environmentIsLocked(reader_ctx) ||
      !identical(parent.env(reader_ctx), baseenv()) ||
      !identical(ls(reader_ctx, all.names = TRUE), c("name", "state")) ||
      !bindingIsLocked("name", reader_ctx) || !bindingIsLocked("state", reader_ctx) ||
      !identical(reader_ctx$name, name)) {
    cdid_ipw_enclosure_abort("Invalid scratch cache reader.",
      c("cdid_ipw_enclosure_scratch_error", "cdid_ipw_enclosure_integrity_error"))
  }
  state <- reader_ctx$state
  if (!is.environment(state) || !environmentIsLocked(state) ||
      !identical(parent.env(state), emptyenv()) ||
      !identical(attr(state, "scratch_ref", exact = TRUE), scratch) ||
      !exists(name, envir = state, inherits = FALSE) ||
      bindingIsLocked(name, state) || bindingIsActive(name, state)) {
    cdid_ipw_enclosure_abort("Invalid scratch cache owner.",
      c("cdid_ipw_enclosure_scratch_error", "cdid_ipw_enclosure_integrity_error"))
  }
  assign(name, value, envir = state)
  invisible(NULL)
}

ipw_validate_if_chunk_scratch <- function(scratch, problem, s_idx, e_idx) {
  if (!is.environment(scratch) || !inherits(scratch, "cdid_ipw_if_chunk_scratch") ||
      !environmentIsLocked(scratch) ||
      !identical(sort(ls(scratch, all.names = TRUE)), sort(c("seal_token", "problem_ref", "storage_mode", "a_start", "a_end", "b_len", "rate_done", "rate", "obs_done", "obs", "outcome_free_chunk", "shared_xi_chunk", "copy_work_ops"))) ||
      !all(vapply(ls(scratch, all.names = TRUE), function(nm) bindingIsLocked(nm, scratch), logical(1))) ||
      !exists("seal_token", envir = scratch, inherits = FALSE) ||
      !identical(get("seal_token", envir = scratch, inherits = FALSE), .ipw_access_seal_token)) {
    cdid_ipw_enclosure_abort(
      "Malformed or foreign `chunk_scratch` object; reuse is refused.",
      c("cdid_ipw_enclosure_scratch_error", "cdid_ipw_enclosure_integrity_error")
    )
  }
  cache_names <- c("rate_done", "rate", "obs_done", "obs", "outcome_free_chunk", "shared_xi_chunk", "copy_work_ops")
  state <- NULL
  for (nm in cache_names) {
    if (!bindingIsActive(nm, scratch)) {
      cdid_ipw_enclosure_abort("Missing scratch cache reader.",
        c("cdid_ipw_enclosure_scratch_error", "cdid_ipw_enclosure_integrity_error"))
    }
    reader_ctx <- environment(activeBindingFunction(nm, scratch))
    if (!is.environment(reader_ctx) || !environmentIsLocked(reader_ctx) ||
        !identical(parent.env(reader_ctx), baseenv()) ||
        !identical(ls(reader_ctx, all.names = TRUE), c("name", "state")) ||
        !bindingIsLocked("name", reader_ctx) || !bindingIsLocked("state", reader_ctx) ||
        !identical(reader_ctx$name, nm) ||
        !is.environment(reader_ctx$state) || !environmentIsLocked(reader_ctx$state) ||
        !identical(parent.env(reader_ctx$state), emptyenv()) ||
        !identical(attr(reader_ctx$state, "scratch_ref", exact = TRUE), scratch) ||
        !identical(sort(ls(reader_ctx$state, all.names = TRUE)), sort(cache_names)) ||
        bindingIsLocked(nm, reader_ctx$state) || bindingIsActive(nm, reader_ctx$state) ||
        (!is.null(state) && !identical(reader_ctx$state, state))) {
      cdid_ipw_enclosure_abort("Malformed or foreign scratch cache state.",
        c("cdid_ipw_enclosure_scratch_error", "cdid_ipw_enclosure_integrity_error"))
    }
    state <- reader_ctx$state
  }
  if (!identical(scratch$problem_ref, problem$ref) ||
      !identical(scratch$a_start, as.integer(s_idx)) ||
      !identical(scratch$a_end, as.integer(e_idx))) {
    cdid_ipw_enclosure_abort(
      "`chunk_scratch` belongs to a different nuisance problem or chunk range `[a_start, a_end]`; reuse is refused.",
      c("cdid_ipw_enclosure_scratch_error", "cdid_ipw_enclosure_integrity_error", "cdid_ipw_enclosure_problem_mismatch_error")
    )
  }
  invisible(scratch)
}

ipw_if_scratch_get_rate <- function(
  scratch,
  u,
  wu_iv,
  su,
  ju,
  r_cohort,
  q_cohort,
  M_cohort,
  n_iv
) {
  if (!is.null(scratch) && scratch$rate_done[u]) {
    return(list(val = scratch$rate[[u]], reused = TRUE))
  }
  s_minus_r <- if (su) gmm_iv_sub(gmm_iv(1.0), r_cohort[ju, ]) else gmm_iv_neg(r_cohort[ju, ])
  xir_1 <- gmm_iv_div(gmm_iv_mul(wu_iv, s_minus_r), q_cohort[ju, ])
  xir_2 <- gmm_iv_mul(gmm_iv_div(gmm_iv_mul(n_iv, wu_iv), M_cohort[ju, ]), s_minus_r)
  xir_iv <- logit_iv_intersect(xir_1, xir_2, "xi_rate_ij")
  val <- list(
    s_minus_r = s_minus_r,
    xir_iv = xir_iv,
    ops = 6
  )
  if (!is.null(scratch)) {
    values <- scratch$rate; values[[u]] <- val
    done <- scratch$rate_done; done[u] <- TRUE
    ipw_if_scratch_set(scratch, "rate", values)
    ipw_if_scratch_set(scratch, "rate_done", done)
    ipw_if_scratch_set(scratch, "copy_work_ops", scratch$copy_work_ops + 2 * scratch$b_len)
  }
  list(val = val, reused = FALSE)
}

ipw_if_scratch_get_obs <- function(
  scratch,
  u,
  a_u,
  ju,
  wu_iv,
  z_u,
  M_cohort,
  Pi_cohort,
  M_over_Pi_cohort,
  F_dbl,
  F_inv_T,
  a_hi,
  a_c_hi,
  rho_hi,
  x_bar_hi,
  beta0_dbl,
  cert_ctrl,
  p
) {
  if (!is.null(scratch) && scratch$obs_done[u]) {
    return(list(val = scratch$obs[[u]], reused = TRUE))
  }
  D_u <- if (ju == 1L) 1L else 0L
  cw <- ipw_enclosure_c_weight(
    wu_iv, M_cohort[ju, ], Pi_cohort[ju, ], M_over_Pi_cohort[ju, ],
    unit_offset = a_u, cohort_index = ju
  )
  c_iv <- cw$c_iv
  dn_ops <- 0
  dn_solves <- 0L
  if (all(z_u == 0.0)) {
    ai_hi <- 0.0
  } else {
    dn_z <- ipw_enclosure_dual_norm_iv(F_dbl, F_inv_T, z_u)
    dn_ops <- dn_z$ops
    dn_solves <- dn_z$solves
    ai_hi <- min(dn_z$upper, if (D_u == 0L) min(a_c_hi, a_hi) else a_hi)
  }
  xi_hi <- min(x_bar_hi, logit_b64_mul(gmm_iv(ai_hi), gmm_iv(rho_hi))[2L])
  one_m_xi_lo <- logit_b64_sub(gmm_iv(1.0), gmm_iv(xi_hi))[1L]

  eta_iv <- gmm_iv_dot(z_u, beta0_dbl)
  exp_res <- logit_iv_exp(
    eta_iv,
    taylor_degree = cert_ctrl$taylor_degree,
    max_abs_eta = cert_ctrl$max_abs_eta,
    max_squarings = cert_ctrl$max_squarings
  )
  o0_iv <- exp_res$iv
  logis0 <- logit_iv_logistic_from_odds(o0_iv)
  o_star_iv <- c(
    logit_b64_mul(gmm_iv(o0_iv[1L]), gmm_iv(one_m_xi_lo))[1L],
    logit_b64_div(gmm_iv(o0_iv[2L]), gmm_iv(one_m_xi_lo))[2L]
  )
  logis_star <- logit_iv_logistic_from_odds(o_star_iv)

  dmp0_iv <- if (D_u == 1L) logis0$one_minus_prob else gmm_iv_neg(logis0$prob)
  dmp_star_iv <- if (D_u == 1L) logis_star$one_minus_prob else gmm_iv_neg(logis_star$prob)

  base_ops <- cw$ops + dn_ops + 2 + 2 * p + exp_res$ops + logis0$ops + logis_star$ops + 2
  val <- list(
    D_u = D_u,
    c_iv = c_iv,
    ai_hi = ai_hi,
    xi_hi = xi_hi,
    one_m_xi_lo = one_m_xi_lo,
    eta_iv = eta_iv,
    o0_iv = o0_iv,
    o_star_iv = o_star_iv,
    logis0 = logis0,
    logis_star = logis_star,
    dmp0_iv = dmp0_iv,
    dmp_star_iv = dmp_star_iv,
    base_ops = base_ops,
    solves = dn_solves,
    exp_evals = 2L
  )
  if (!is.null(scratch)) {
    values <- scratch$obs; values[[u]] <- val
    done <- scratch$obs_done; done[u] <- TRUE
    ipw_if_scratch_set(scratch, "obs", values)
    ipw_if_scratch_set(scratch, "obs_done", done)
    ipw_if_scratch_set(scratch, "copy_work_ops", scratch$copy_work_ops + 2 * scratch$b_len)
  }
  list(val = val, reused = FALSE)
}

construct_dense_ipw_shared_xi_access_context <- function(
  problem,
  s_idx,
  e_idx,
  ctrl,
  pre,
  has_xi_cand,
  R_H_cand,
  xi_cand_mat,
  scratch = NULL
) {
  raw <- problem$raw
  n <- nrow(raw$Z)
  p <- ncol(raw$Z)
  m <- length(raw$cohorts)
  a_seq <- seq.int(as.integer(s_idx), as.integer(e_idx))

  read_chunk <- function() {
    list(
      a_index = a_seq,
      id = raw$id[a_seq],
      g = raw$G[a_seq],
      weight = raw$w[a_seq],
      pair_observed = raw$pair_observed[a_seq],
      cohort_index = problem$masks$cohort_index[a_seq],
      in_risk = !is.na(problem$masks$cohort_index[a_seq]),
      in_observed = !is.na(problem$masks$cohort_index[a_seq]) & raw$pair_observed[a_seq],
      Z_chunk = raw$Z[a_seq, , drop = FALSE],
      n_endpoints_read = 0L,
      reader_work_ops = 0
    )
  }

  ipw_seal_access_context(
    list(
      mode = "dense",
      problem = problem,
      header = list(
        n = n, p = p, m = m,
        s_idx = as.integer(s_idx), e_idx = as.integer(e_idx),
        b_len = as.integer(e_idx - s_idx + 1L),
        z_colnames = colnames(raw$Z),
        cohort_labels = raw$cohort_labels
      ),
      ctrl = ctrl,
      reader_ctrl = NULL,
      pre = pre,
      has_xi_cand = has_xi_cand,
      R_H_cand = R_H_cand,
      xi_cand_mat = xi_cand_mat,
      scratch = scratch,
      read_chunk = read_chunk
    ),
    phase = "shared_xi_chunk"
  )
}

construct_view_ipw_shared_xi_access_context <- function(
  problem,
  s_idx,
  e_idx,
  ctrl,
  reader_ctrl,
  pre,
  has_xi_cand,
  R_H_cand,
  xi_cand_mat,
  expected_a_id = NULL,
  trace_env = NULL,
  scratch = NULL
) {
  raw <- problem$raw
  n <- as.integer(raw$n_a_units)
  p <- as.integer(raw$p_covariates)
  m <- length(raw$cohorts)
  nuisance_view <- get("nuisance_view", envir = problem$nuisance_ref, inherits = FALSE)

  read_chunk <- function() {
    if (!is.null(scratch) && !is.null(scratch$outcome_free_chunk)) {
      cached <- scratch$outcome_free_chunk
      cached$reader_work_ops <- 0
      return(cached)
    }
    ch <- query_ipw_view_chunk(
      view = nuisance_view,
      a_start = s_idx,
      a_end = e_idx,
      control = reader_ctrl,
      read_outcomes = FALSE,
      include_design = TRUE,
      enclose_differences = FALSE,
      expected_a_id = expected_a_id,
      trace_env = trace_env
    )
    res <- list(
      a_index = ch$a_index,
      id = ch$id,
      g = ch$g,
      weight = ch$weight,
      pair_observed = ch$in_observed,
      cohort_index = ch$cohort_index,
      in_risk = ch$in_risk,
      in_observed = ch$in_observed,
      Z_chunk = ch$design$Z_chunk,
      n_endpoints_read = ch$read_diagnostics$n_endpoints_read,
      reader_work_ops = ch$work$actual_work_ops
    )
    if (!is.null(scratch)) ipw_if_scratch_set(scratch, "outcome_free_chunk", res)
    res
  }

  ipw_seal_access_context(
    list(
      mode = "compact_view",
      problem = problem,
      header = list(
        n = n, p = p, m = m,
        s_idx = as.integer(s_idx), e_idx = as.integer(e_idx),
        b_len = as.integer(e_idx - s_idx + 1L),
        z_colnames = raw$z_colnames,
        cohort_labels = raw$cohort_labels
      ),
      ctrl = ctrl,
      reader_ctrl = reader_ctrl,
      pre = pre,
      has_xi_cand = has_xi_cand,
      R_H_cand = R_H_cand,
      xi_cand_mat = xi_cand_mat,
      scratch = scratch,
      read_chunk = read_chunk
    ),
    phase = "shared_xi_chunk"
  )
}

ipw_enclosure_shared_xi_run_with_context <- function(ctx) {
  ipw_validate_access_context(ctx, expected_phase = "shared_xi_chunk")

  problem <- ctx$problem
  raw <- problem$raw
  n <- ctx$header$n
  p <- ctx$header$p
  m <- ctx$header$m
  s_idx <- ctx$header$s_idx
  e_idx <- ctx$header$e_idx
  b_len <- ctx$header$b_len
  z_colnames <- ctx$header$z_colnames
  cohort_labels <- ctx$header$cohort_labels
  pre <- ctx$pre
  has_xi_cand <- ctx$has_xi_cand
  R_H_cand <- ctx$R_H_cand
  xi_cand_mat <- ctx$xi_cand_mat
  is_compact <- identical(ctx$mode, "compact_view")

  scratch <- ctx$scratch
  if (is.null(scratch) && is_compact) {
    scratch <- new_ipw_if_chunk_scratch(problem, s_idx, e_idx, ctx$ctrl)
  } else if (!is.null(scratch)) {
    ipw_validate_if_chunk_scratch(scratch, problem, s_idx, e_idx)
  }
  cache_work_at_entry <- ipw_if_scratch_work(scratch)

  ch <- ctx$read_chunk()
  a_seq <- ch$a_index

  if (!isTRUE(problem$certified)) {
    inc_fields <- list(
      certified = FALSE,
      status = "INCONCLUSIVE",
      inconclusive = TRUE,
      inconclusive_reasons = problem$inconclusive_reasons,
      storage_mode = if (is_compact) "compact_view" else NULL,
      a_index = a_seq,
      id = ch$id,
      problem_ref = problem$ref,
      preflight = pre
    )
    if (!is_compact) inc_fields$storage_mode <- NULL
    return(structure(inc_fields, class = "cdid_ipw_shared_xi_chunk"))
  }

  nuis <- problem$nuisance_enclosures
  cert <- problem$certificate
  M_cohort <- cert$enclosures$M_cohort
  Pi_cohort <- cert$enclosures$Pi_cohort
  r_cohort <- cert$enclosures$r_cohort
  q_cohort <- cert$enclosures$q_cohort
  M_over_Pi_cohort <- cert$enclosures$M_over_Pi_cohort
  F_dbl <- raw$F
  F_inv_T <- nuis$F_inv_T
  coord_dual_iv <- nuis$coord_dual_iv
  beta0_dbl <- raw$beta0
  cert_ctrl <- problem$controls$certificate_control

  xi_rate_iv <- array(0.0, dim = c(b_len, m, 2L), dimnames = list(NULL, cohort_labels, c("lo", "hi")))
  u_0_iv <- array(0.0, dim = c(b_len, p, 2L), dimnames = list(NULL, z_colnames, c("lo", "hi")))
  u_star_iv <- array(0.0, dim = c(b_len, p, 2L), dimnames = list(NULL, z_colnames, c("lo", "hi")))
  rhs_0_iv <- array(0.0, dim = c(b_len, p, 2L), dimnames = list(NULL, z_colnames, c("lo", "hi")))
  rhs_star_iv <- array(0.0, dim = c(b_len, p, 2L), dimnames = list(NULL, z_colnames, c("lo", "hi")))

  xi_beta_0_iv <- if (has_xi_cand) array(0.0, dim = c(b_len, p, 2L), dimnames = list(NULL, z_colnames, c("lo", "hi"))) else NULL
  xi_beta_star_iv <- if (has_xi_cand) array(0.0, dim = c(b_len, p, 2L), dimnames = list(NULL, z_colnames, c("lo", "hi"))) else NULL
  xi_beta_0_rad <- if (has_xi_cand) numeric(b_len) else NULL
  xi_beta_star_rad <- if (has_xi_cand) numeric(b_len) else NULL
  xi_beta_center <- if (has_xi_cand) matrix(0.0, nrow = b_len, ncol = p, dimnames = list(NULL, z_colnames)) else NULL

  Z_ch <- ch$Z_chunk
  w_ch <- ch$weight
  po_ch <- ch$pair_observed
  j_ch <- ch$cohort_index
  n_iv <- gmm_iv(as.double(n))

  a_hi <- cert$a
  a_c_hi <- cert$a_control
  rho_hi <- cert$rho
  x_bar_hi <- nuis$x_bar
  actual_ops <- 0
  physical_ops <- 0
  reused_ops <- 0
  tri_solves <- 0L
  physical_tri_solves <- 0L
  reused_tri_solves <- 0L
  exp_evals <- 0L
  physical_exp_evals <- 0L
  reused_exp_evals <- 0L

  shared_err <- ipw_enclosure_capture_arithmetic({
    for (u in seq_len(b_len)) {
      if (u > 1L && (u %% 256L) == 0L) ipw_enclosure_check_interrupt()
      ju <- j_ch[u]
      if (is.na(ju)) next # Outside R: xi_rate = 0, u = 0, xi_beta = 0

      wu_iv <- gmm_iv(w_ch[u])
      su <- po_ch[u]

      # xi_{r, ij} = w_i * (S_i - r_j) / q_j (nonzero for risk non-observers S_i = FALSE!)
      rate_entry <- ipw_if_scratch_get_rate(
        scratch = scratch, u = u, wu_iv = wu_iv, su = su, ju = ju,
        r_cohort = r_cohort, q_cohort = q_cohort, M_cohort = M_cohort, n_iv = n_iv
      )
      xir_iv <- rate_entry$val$xir_iv
      xi_rate_iv[u, ju, ] <- xir_iv
      actual_ops <- actual_ops + rate_entry$val$ops
      if (rate_entry$reused) {
        reused_ops <- reused_ops + rate_entry$val$ops
      } else {
        physical_ops <- physical_ops + rate_entry$val$ops
      }

      if (su) {
        z_u <- Z_ch[u, ]
        obs_entry <- ipw_if_scratch_get_obs(
          scratch = scratch, u = u, a_u = a_seq[u], ju = ju, wu_iv = wu_iv, z_u = z_u,
          M_cohort = M_cohort, Pi_cohort = Pi_cohort, M_over_Pi_cohort = M_over_Pi_cohort,
          F_dbl = F_dbl, F_inv_T = F_inv_T, a_hi = a_hi, a_c_hi = a_c_hi,
          rho_hi = rho_hi, x_bar_hi = x_bar_hi, beta0_dbl = beta0_dbl,
          cert_ctrl = cert_ctrl, p = p
        )
        obs_u <- obs_entry$val
        actual_ops <- actual_ops + obs_u$base_ops
        tri_solves <- tri_solves + obs_u$solves
        exp_evals <- exp_evals + obs_u$exp_evals
        if (obs_entry$reused) {
          reused_ops <- reused_ops + obs_u$base_ops
          reused_tri_solves <- reused_tri_solves + obs_u$solves
          reused_exp_evals <- reused_exp_evals + obs_u$exp_evals
        } else {
          physical_ops <- physical_ops + obs_u$base_ops
          physical_tri_solves <- physical_tri_solves + obs_u$solves
          physical_exp_evals <- physical_exp_evals + obs_u$exp_evals
        }

        c_dmp0 <- gmm_iv_mul(obs_u$c_iv, obs_u$dmp0_iv)
        c_dmp_star <- gmm_iv_mul(obs_u$c_iv, obs_u$dmp_star_iv)
        actual_ops <- actual_ops + 2
        physical_ops <- physical_ops + 2

        for (k in seq_len(p)) {
          zk_iv <- gmm_iv(z_u[k])
          u_0_iv[u, k, ] <- gmm_iv_mul(c_dmp0, zk_iv)
          u_star_iv[u, k, ] <- gmm_iv_mul(c_dmp_star, zk_iv)
          actual_ops <- actual_ops + 2
          physical_ops <- physical_ops + 2
        }
      }

      # rhs_i = u_i + L_{j(i)} * xi_{r, i, j(i)}
      for (k in seq_len(p)) {
        l0_xir <- gmm_iv_mul(nuis$L_0[k, ju, ], xir_iv)
        ls_xir <- gmm_iv_mul(nuis$L_star[k, ju, ], xir_iv)
        rhs_0_iv[u, k, ] <- gmm_iv_add(u_0_iv[u, k, ], l0_xir)
        rhs_star_iv[u, k, ] <- gmm_iv_add(u_star_iv[u, k, ], ls_xir)
        actual_ops <- actual_ops + 4
        physical_ops <- physical_ops + 4
      }

      if (has_xi_cand) {
        # Determine candidate center tilde_xi_beta_i (either supplied directly or via R_H_cand)
        cand_u <- if (!is.null(xi_cand_mat)) {
          xi_cand_mat[u, ]
        } else {
          rhs_mid <- 0.5 * (rhs_0_iv[u, , 1L] + rhs_0_iv[u, , 2L])
          drop(backsolve(R_H_cand, forwardsolve(t(R_H_cand), rhs_mid)))
        }
        if (any(!is.finite(cand_u))) {
          ipw_enclosure_arith_refuse("Generated xi_beta candidate overflowed.",
            reason = "generated_xi_beta_candidate_overflow")
        }
        xi_beta_center[u, ] <- cand_u

        # Solve H * xi_beta = rhs_i <=> H * xi_beta + (-rhs_i) = 0 via Lemma dual(iii)
        neg_rhs_0 <- cbind(-rhs_0_iv[u, , 2L], -rhs_0_iv[u, , 1L])
        neg_rhs_star <- cbind(-rhs_star_iv[u, , 2L], -rhs_star_iv[u, , 1L])

        sol0 <- ipw_enclosure_lemma_dual_solve(
          F_dbl, F_inv_T, coord_dual_iv, nuis$H_0, neg_rhs_0, cand_u, nuis$kappa_0_lower
        )
        sols <- ipw_enclosure_lemma_dual_solve(
          F_dbl, F_inv_T, coord_dual_iv, nuis$H_star, neg_rhs_star, cand_u, nuis$kappa_star_lower
        )
        xi_beta_0_iv[u, , ] <- sol0$coord_iv
        xi_beta_star_iv[u, , ] <- sols$coord_iv
        xi_beta_0_rad[u] <- sol0$radius_m_norm
        xi_beta_star_rad[u] <- sols$radius_m_norm
        actual_ops <- actual_ops + sol0$ops + sols$ops
        physical_ops <- physical_ops + sol0$ops + sols$ops
        tri_solves <- tri_solves + sol0$solves + sols$solves
        physical_tri_solves <- physical_tri_solves + sol0$solves + sols$solves
      }
    }

    NULL
  })
  if (inherits(shared_err, "error")) {
    return(ipw_enclosure_inconclusive_result(
      shared_err, "cdid_ipw_shared_xi_chunk",
      list(
        storage_mode = if (is_compact) "compact_view" else NULL,
        a_index = a_seq,
        id = ch$id,
        problem_ref = problem$ref,
        preflight = pre
      ),
      "shared_xi_chunk"
    ))
  }
  stopifnot(actual_ops + ipw_if_scratch_work(scratch) - cache_work_at_entry <= pre$work_ops_bound)
  if (is_compact) {
    stopifnot((actual_ops + ch$reader_work_ops + ipw_if_scratch_work(scratch) - cache_work_at_entry) <= pre$combined_work_ops_bound)
  }

  out_fields <- list(
    certified = TRUE,
    status = "CERTIFIED",
    storage_mode = if (is_compact) "compact_view" else NULL,
    a_index = a_seq,
    id = ch$id,
    g = ch$g,
    weight = w_ch,
    pair_observed = po_ch,
    cohort_index = j_ch,
    in_risk = ch$in_risk,
    in_observed = ch$in_observed,
    xi_rate = xi_rate_iv,
    u_0 = u_0_iv,
    u_star = u_star_iv,
    rhs_0 = rhs_0_iv,
    rhs_star = rhs_star_iv,
    has_xi_beta = has_xi_cand,
    xi_beta_center = xi_beta_center,
    xi_beta_0 = xi_beta_0_iv,
    xi_beta_star = xi_beta_star_iv,
    xi_beta_0_radius_m_norm = xi_beta_0_rad,
    xi_beta_star_radius_m_norm = xi_beta_star_rad,
    chunk_scratch = scratch,
    problem_ref = problem$ref,
    preflight = pre,
    read_diagnostics = if (is_compact) list(
      endpoints_read = 0L,
      n_endpoints_read = 0L,
      reader_work_ops = ch$reader_work_ops
    ) else NULL,
    work = list(
      chunk_rows = b_len,
      triangular_solves = logit_cert_count_output(tri_solves),
      physical_triangular_solves = logit_cert_count_output(physical_tri_solves),
      reused_triangular_solves = logit_cert_count_output(reused_tri_solves),
      exp_endpoint_evaluations = logit_cert_count_output(exp_evals),
      physical_exp_endpoint_evaluations = logit_cert_count_output(physical_exp_evals),
      reused_exp_endpoint_evaluations = logit_cert_count_output(reused_exp_evals),
      actual_interval_ops = actual_ops,
      logical_interval_ops = actual_ops,
      physical_interval_ops = physical_ops,
      reused_shared_interval_ops = reused_ops,
      actual_reader_work_ops = ch$reader_work_ops,
      cache_copy_work_ops_charged = ipw_if_scratch_work(scratch) - cache_work_at_entry,
      combined_work_ops = actual_ops + ch$reader_work_ops + ipw_if_scratch_work(scratch) - cache_work_at_entry,
      endpoints_read = 0L
    )
  )
  if (!is_compact) {
    out_fields$storage_mode <- NULL
    out_fields$read_diagnostics <- NULL
  }
  out_chunk <- structure(out_fields, class = "cdid_ipw_shared_xi_chunk")
  sealed_chunk <- ipw_enclosure_seal(out_chunk, "shared_xi_chunk")
  if (!is.null(scratch)) ipw_if_scratch_set(scratch, "shared_xi_chunk", sealed_chunk)
  sealed_chunk
}

enclose_ipw_shared_xi_chunk <- function(
  problem,
  a_start,
  a_end,
  control,
  xi_beta_candidate = NULL,
  H_chol_candidate = NULL,
  require_xi_beta = TRUE
) {
  ipw_enclosure_require_deps()
  if (missing(problem) || missing(a_start) || missing(a_end) || missing(control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`problem`, `a_start`, `a_end`, `control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }
  ipw_validate_nuisance_problem_mode(problem, require_certified = FALSE, expected_storage_mode = "dense")
  scratch <- attr(control, "chunk_scratch", exact = TRUE)
  ctrl <- validate_ipw_link_enclosure_control(control)
  if (!gmm_scalar_logical(require_xi_beta)) {
    cdid_ipw_enclosure_abort("`require_xi_beta` must be a scalar logical.", "cdid_ipw_enclosure_input_error")
  }

  raw <- problem$raw
  n <- nrow(raw$Z)
  p <- ncol(raw$Z)
  m <- length(raw$cohorts)

  s_idx <- ipw_enclosure_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- ipw_enclosure_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > n) {
    cdid_ipw_enclosure_abort(
      paste0("Chunk range `[a_start, a_end]` must satisfy `1 <= a_start <= a_end <= ", n, "`."),
      "cdid_ipw_enclosure_bound_error"
    )
  }
  b_len <- e_idx - s_idx + 1

  has_xi_cand <- !is.null(xi_beta_candidate) || !is.null(H_chol_candidate)
  if (isTRUE(require_xi_beta) && !has_xi_cand) {
    cdid_ipw_enclosure_abort(
      "Enclosing `xi_beta = H^{-1}(u + L * xi_rate)` via Lemma dual(iii) requires an explicit stored candidate (`xi_beta_candidate` or `H_chol_candidate`); general interval matrix inversion is not supported.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.null(H_chol_candidate)) {
    if (!is.matrix(H_chol_candidate) || is.object(H_chol_candidate) || !is.numeric(H_chol_candidate) ||
        nrow(H_chol_candidate) != p || ncol(H_chol_candidate) != p) {
      cdid_ipw_enclosure_abort(
        paste0("`H_chol_candidate` must be a numeric p x p matrix (p = ", p, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (!is.null(xi_beta_candidate)) {
    if (!is.matrix(xi_beta_candidate) || is.object(xi_beta_candidate) || !is.numeric(xi_beta_candidate) ||
        ncol(xi_beta_candidate) != p || !(nrow(xi_beta_candidate) %in% c(b_len, n))) {
      cdid_ipw_enclosure_abort(
        paste0("`xi_beta_candidate` must be a finite numeric matrix with ncol = p (", p, ") and nrow = chunk length (", b_len, ") or n (", n, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }

  pre <- (if (is.null(scratch)) preflight_ipw_link_if_chunk else preflight_ipw_link_if_chunk_cached)(
    b_len = b_len, p = p, m = m, control = ctrl,
    certificate_control = problem$controls$certificate_control,
    include_components = FALSE, include_shared_xi = TRUE,
    has_phi_candidate = FALSE, has_xi_beta = has_xi_cand
  )

  a_seq <- seq.int(as.integer(s_idx), as.integer(e_idx))
  R_H_cand <- NULL
  if (!is.null(H_chol_candidate)) {
    R_H_cand <- validate_logit_certificate_metric_factor(H_chol_candidate, p)
  }
  xi_cand_mat <- NULL
  if (!is.null(xi_beta_candidate)) {
    xi_cand_mat <- if (nrow(xi_beta_candidate) == n) {
      matrix(as.double(xi_beta_candidate[a_seq, , drop = FALSE]), nrow = b_len, ncol = p)
    } else matrix(as.double(xi_beta_candidate), nrow = b_len, ncol = p)
    if (anyNA(xi_cand_mat) || any(!is.finite(xi_cand_mat))) {
      cdid_ipw_enclosure_abort("Queried `xi_beta_candidate` rows must be finite.",
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error"))
    }
  }

  ctx <- construct_dense_ipw_shared_xi_access_context(
    problem = problem,
    s_idx = s_idx,
    e_idx = e_idx,
    ctrl = ctrl,
    pre = pre,
    has_xi_cand = has_xi_cand,
    R_H_cand = R_H_cand,
    xi_cand_mat = xi_cand_mat,
    scratch = scratch
  )
  ipw_enclosure_shared_xi_run_with_context(ctx)
}

enclose_ipw_shared_xi_view_chunk <- function(
  problem,
  a_start,
  a_end,
  control,
  reader_control,
  xi_beta_candidate = NULL,
  H_chol_candidate = NULL,
  require_xi_beta = TRUE
) {
  ipw_enclosure_require_deps()
  if (missing(problem) || missing(a_start) || missing(a_end) ||
      missing(control) || missing(reader_control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`problem`, `a_start`, `a_end`, `control`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }
  ipw_validate_nuisance_problem_mode(problem, require_certified = FALSE, expected_storage_mode = "compact_view")
  scratch <- attr(control, "chunk_scratch", exact = TRUE)
  ctrl <- validate_ipw_link_enclosure_control(control)
  expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  trace_env <- attr(reader_control, "trace_env", exact = TRUE)
  reader_ctrl <- validate_ipw_owner_control(reader_control)
  if (!gmm_scalar_logical(require_xi_beta)) {
    cdid_ipw_enclosure_abort("`require_xi_beta` must be a scalar logical.", "cdid_ipw_enclosure_input_error")
  }

  raw <- problem$raw
  n <- as.integer(raw$n_a_units)
  p <- as.integer(raw$p_covariates)
  m <- length(raw$cohorts)
  panel_owner <- get("panel_owner", envir = problem$panel_ref, inherits = FALSE)

  s_idx <- ipw_enclosure_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- ipw_enclosure_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > n) {
    cdid_ipw_enclosure_abort(
      paste0("Chunk range `[a_start, a_end]` must satisfy `1 <= a_start <= a_end <= ", n, "`."),
      "cdid_ipw_enclosure_bound_error"
    )
  }
  b_len <- e_idx - s_idx + 1

  ipw_owner_refuse_decorative_attrs(xi_beta_candidate, "xi_beta_candidate")
  ipw_owner_refuse_decorative_attrs(H_chol_candidate, "H_chol_candidate")
  has_xi_cand <- !is.null(xi_beta_candidate) || !is.null(H_chol_candidate)
  if (isTRUE(require_xi_beta) && !has_xi_cand) {
    cdid_ipw_enclosure_abort(
      "Enclosing `xi_beta = H^{-1}(u + L * xi_rate)` via Lemma dual(iii) requires an explicit stored candidate (`xi_beta_candidate` or `H_chol_candidate`); general interval matrix inversion is not supported.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }
  if (!is.null(H_chol_candidate)) {
    if (!is.matrix(H_chol_candidate) || is.object(H_chol_candidate) || !is.numeric(H_chol_candidate) ||
        nrow(H_chol_candidate) != p || ncol(H_chol_candidate) != p) {
      cdid_ipw_enclosure_abort(
        paste0("`H_chol_candidate` must be a numeric p x p matrix (p = ", p, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  caller_xi_entries <- 0
  if (!is.null(xi_beta_candidate)) {
    if (!is.matrix(xi_beta_candidate) || is.object(xi_beta_candidate) || !is.numeric(xi_beta_candidate) ||
        ncol(xi_beta_candidate) != p || !(nrow(xi_beta_candidate) %in% c(b_len, n))) {
      cdid_ipw_enclosure_abort(
        paste0("`xi_beta_candidate` must be a finite numeric matrix with ncol = p (", p, ") and nrow = chunk length (", b_len, ") or n (", n, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
    caller_xi_entries <- nrow(xi_beta_candidate) * p
  }
  if (!is.null(expected_a_id)) {
    ipw_owner_refuse_decorative_attrs(expected_a_id, "expected_a_id")
    if (is.factor(expected_a_id) || is.object(expected_a_id) || !is.atomic(expected_a_id) ||
        !is.null(dim(expected_a_id)) || !(length(expected_a_id) %in% c(b_len, n))) {
      cdid_ipw_owner_abort(
        "`expected_a_id` must be NULL or a plain atomic vector of length `b_len` or `n_a_units`.",
        c("cdid_ipw_owner_id_error", "cdid_ipw_owner_order_mismatch_error")
      )
    }
  }

  pre <- preflight_ipw_link_if_view_chunk(
    b_len = b_len,
    p = p,
    m = m,
    control = ctrl,
    reader_control = reader_ctrl,
    certificate_control = problem$controls$certificate_control,
    include_components = FALSE,
    include_shared_xi = TRUE,
    has_phi_candidate = FALSE,
    has_xi_beta = has_xi_cand,
    read_outcomes = FALSE,
    panel_owner = panel_owner,
    a_seq = seq.int(as.integer(s_idx), as.integer(e_idx)),
    expected_a_id = expected_a_id,
    caller_phi_entries = 0,
    caller_xi_entries = caller_xi_entries
  )

  a_seq <- seq.int(as.integer(s_idx), as.integer(e_idx))

  R_H_cand <- NULL
  if (!is.null(H_chol_candidate)) {
    R_H_cand <- validate_logit_certificate_metric_factor(H_chol_candidate, p)
  }
  xi_cand_mat <- NULL
  if (!is.null(xi_beta_candidate)) {
    xi_cand_mat <- if (nrow(xi_beta_candidate) == n) {
      matrix(as.double(xi_beta_candidate[a_seq, , drop = FALSE]), nrow = b_len, ncol = p)
    } else matrix(as.double(xi_beta_candidate), nrow = b_len, ncol = p)
    if (anyNA(xi_cand_mat) || any(!is.finite(xi_cand_mat))) {
      cdid_ipw_enclosure_abort("Queried `xi_beta_candidate` rows must be finite.",
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error"))
    }
  }

  ctx <- construct_view_ipw_shared_xi_access_context(
    problem = problem,
    s_idx = s_idx,
    e_idx = e_idx,
    ctrl = ctrl,
    reader_ctrl = reader_ctrl,
    pre = pre,
    has_xi_cand = has_xi_cand,
    R_H_cand = R_H_cand,
    xi_cand_mat = xi_cand_mat,
    expected_a_id = expected_a_id,
    trace_env = trace_env,
    scratch = scratch
  )
  ipw_enclosure_shared_xi_run_with_context(ctx)
}

construct_dense_ipw_if_chunk_access_context <- function(
  link_enclosure,
  prob,
  s_idx,
  e_idx,
  ctrl,
  pre,
  include_components,
  include_shared_xi,
  has_phi_cand,
  phi_cand_sub,
  xi_beta_candidate,
  R_H_cand,
  scratch = NULL,
  shared_xi_chunk = NULL
) {
  raw <- prob$raw
  n <- nrow(raw$Z)
  p <- ncol(raw$Z)
  m <- length(raw$cohorts)
  a_seq <- seq.int(as.integer(s_idx), as.integer(e_idx))
  ep <- link_enclosure$endpoints

  read_chunk <- function() {
    obs_pos_ch <- ep$a_to_obs_pos[a_seq]
    obs_local <- which(obs_pos_ch > 0L)
    obs_global <- obs_pos_ch[obs_local]
    chunk_to_obs_pos <- integer(length(a_seq))
    if (length(obs_local) > 0L) chunk_to_obs_pos[obs_local] <- seq_along(obs_local)
    list(
      a_index = a_seq,
      id = raw$id[a_seq],
      g = raw$G[a_seq],
      weight = raw$w[a_seq],
      pair_observed = raw$pair_observed[a_seq],
      cohort_index = prob$masks$cohort_index[a_seq],
      in_risk = prob$masks$in_risk[a_seq],
      in_observed = prob$masks$in_observed[a_seq],
      Z_chunk = raw$Z[a_seq, , drop = FALSE],
      chunk_to_obs_pos = chunk_to_obs_pos,
      y_left_obs = if (length(obs_global) > 0L) ep$y_left_obs[obs_global] else numeric(0L),
      y_right_obs = if (length(obs_global) > 0L) ep$y_right_obs[obs_global] else numeric(0L),
      n_endpoints_read = 2L * length(obs_global),
      reader_work_ops = 0
    )
  }

  ipw_seal_access_context(
    list(
      mode = "dense",
      link_enclosure = link_enclosure,
      problem = prob,
      header = list(
        n = n, p = p, m = m,
        s_idx = as.integer(s_idx), e_idx = as.integer(e_idx),
        b_len = as.integer(e_idx - s_idx + 1L),
        z_colnames = colnames(raw$Z),
        cohort_labels = raw$cohort_labels
      ),
      ctrl = ctrl,
      reader_ctrl = NULL,
      pre = pre,
      include_components = include_components,
      include_shared_xi = include_shared_xi,
      has_phi_cand = has_phi_cand,
      phi_cand_sub = phi_cand_sub,
      xi_beta_candidate = xi_beta_candidate,
      R_H_cand = R_H_cand,
      scratch = scratch,
      shared_xi_chunk = shared_xi_chunk,
      read_chunk = read_chunk
    ),
    phase = "link_if_chunk"
  )
}

construct_view_ipw_if_chunk_access_context <- function(
  link_enclosure,
  prob,
  s_idx,
  e_idx,
  ctrl,
  reader_ctrl,
  pre,
  include_components,
  include_shared_xi,
  has_phi_cand,
  phi_cand_sub,
  xi_beta_candidate,
  R_H_cand,
  expected_a_id = NULL,
  trace_env = NULL,
  scratch = NULL,
  shared_xi_chunk = NULL
) {
  raw <- prob$raw
  n <- as.integer(raw$n_a_units)
  p <- as.integer(raw$p_covariates)
  m <- length(raw$cohorts)
  link_view <- get("link_view", envir = link_enclosure$link_ref, inherits = FALSE)

  read_chunk <- function() {
    ch <- query_ipw_view_chunk(
      view = link_view,
      a_start = s_idx,
      a_end = e_idx,
      control = reader_ctrl,
      read_outcomes = isTRUE(link_enclosure$certified),
      include_design = isTRUE(link_enclosure$certified),
      enclose_differences = FALSE,
      expected_a_id = expected_a_id,
      trace_env = trace_env
    )
    if (!is.null(trace_env) && isTRUE(link_enclosure$certified)) {
      trace_env$obs_s_read <- ch$read_diagnostics$read_obs_s_index
      trace_env$obs_t_read <- ch$read_diagnostics$read_obs_t_index
    }
    if (!is.null(scratch) && is.null(scratch$outcome_free_chunk)) {
      ipw_if_scratch_set(scratch, "outcome_free_chunk", list(
        a_index = ch$a_index,
        id = ch$id,
        g = ch$g,
        weight = ch$weight,
        pair_observed = ch$in_observed,
        cohort_index = ch$cohort_index,
        in_risk = ch$in_risk,
        in_observed = ch$in_observed,
        Z_chunk = ch$design$Z_chunk,
        n_endpoints_read = 0L,
        reader_work_ops = 0
      ))
    }
    list(
      a_index = ch$a_index,
      id = ch$id,
      g = ch$g,
      weight = ch$weight,
      pair_observed = ch$in_observed,
      cohort_index = ch$cohort_index,
      in_risk = ch$in_risk,
      in_observed = ch$in_observed,
      Z_chunk = ch$design$Z_chunk,
      chunk_to_obs_pos = ch$chunk_to_obs_pos,
      y_left_obs = ch$observed_risk$y_left_obs,
      y_right_obs = ch$observed_risk$y_right_obs,
      n_endpoints_read = ch$read_diagnostics$n_endpoints_read,
      reader_work_ops = ch$work$actual_work_ops
    )
  }

  ipw_seal_access_context(
    list(
      mode = "compact_view",
      link_enclosure = link_enclosure,
      problem = prob,
      header = list(
        n = n, p = p, m = m,
        s_idx = as.integer(s_idx), e_idx = as.integer(e_idx),
        b_len = as.integer(e_idx - s_idx + 1L),
        z_colnames = raw$z_colnames,
        cohort_labels = raw$cohort_labels
      ),
      ctrl = ctrl,
      reader_ctrl = reader_ctrl,
      pre = pre,
      include_components = include_components,
      include_shared_xi = include_shared_xi,
      has_phi_cand = has_phi_cand,
      phi_cand_sub = phi_cand_sub,
      xi_beta_candidate = xi_beta_candidate,
      R_H_cand = R_H_cand,
      expected_a_id = expected_a_id,
      trace_env = trace_env,
      scratch = scratch,
      shared_xi_chunk = shared_xi_chunk,
      read_chunk = read_chunk
    ),
    phase = "link_if_chunk"
  )
}

ipw_enclosure_if_chunk_run_with_context <- function(ctx) {
  ipw_validate_access_context(ctx, expected_phase = "link_if_chunk")

  link_enclosure <- ctx$link_enclosure
  prob <- ctx$problem
  raw <- prob$raw
  n <- ctx$header$n
  p <- ctx$header$p
  m <- ctx$header$m
  s_idx <- ctx$header$s_idx
  e_idx <- ctx$header$e_idx
  b_len <- ctx$header$b_len
  ctrl <- ctx$ctrl
  pre <- ctx$pre
  include_components <- ctx$include_components
  include_shared_xi <- ctx$include_shared_xi
  has_phi_cand <- ctx$has_phi_cand
  phi_cand_sub <- ctx$phi_cand_sub
  xi_beta_candidate <- ctx$xi_beta_candidate
  R_H_cand <- ctx$R_H_cand
  is_compact <- identical(ctx$mode, "compact_view")

  scratch <- ctx$scratch
  if (is.null(scratch) && is_compact) {
    scratch <- new_ipw_if_chunk_scratch(prob, s_idx, e_idx, ctrl)
  } else if (!is.null(scratch)) {
    ipw_validate_if_chunk_scratch(scratch, prob, s_idx, e_idx)
  }
  cache_work_at_entry <- ipw_if_scratch_work(scratch)

  if (!is.null(ctx$shared_xi_chunk)) {
    xi <- ctx$shared_xi_chunk
    if (!inherits(xi, "cdid_ipw_shared_xi_chunk") || !is.list(xi) ||
        !ipw_enclosure_seal_ok(xi, "shared_xi_chunk") ||
        !identical(xi$problem_ref, prob$ref) || !isTRUE(xi$certified) ||
        !identical(xi$a_index, seq.int(as.integer(s_idx), as.integer(e_idx))) ||
        !identical(identical(xi$storage_mode, "compact_view"), is_compact)) {
      cdid_ipw_enclosure_abort("Shared xi reuse requires a sealed chunk from the same problem, storage mode and exact range.",
        c("cdid_ipw_enclosure_shared_xi_error", "cdid_ipw_enclosure_integrity_error"))
    }
  }
  ch <- ctx$read_chunk()
  a_seq <- ch$a_index
  in_r_sub <- ch$in_risk

  if (has_phi_cand && any(phi_cand_sub[!in_r_sub] != 0.0)) {
    cdid_ipw_enclosure_abort(
      "`phi_candidate` has nonzero entries for units outside the risk set R = {G in cohorts}.",
      c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
    )
  }

  if (!isTRUE(link_enclosure$certified)) {
    inc_fields <- list(
      certified = FALSE,
      status = "INCONCLUSIVE",
      inconclusive = TRUE,
      inconclusive_reasons = link_enclosure$inconclusive_reasons,
      storage_mode = if (is_compact) "compact_view" else NULL,
      a_index = a_seq,
      id = ch$id,
      g = ch$g,
      weight = ch$weight,
      preflight = pre
    )
    if (!is_compact) inc_fields$storage_mode <- NULL
    return(structure(inc_fields, class = "cdid_ipw_link_if_chunk"))
  }

  nuis <- prob$nuisance_enclosures
  cert <- prob$certificate
  M_cohort <- cert$enclosures$M_cohort
  Pi_cohort <- cert$enclosures$Pi_cohort
  r_cohort <- cert$enclosures$r_cohort
  q_cohort <- cert$enclosures$q_cohort
  M_over_Pi_cohort <- cert$enclosures$M_over_Pi_cohort
  F_dbl <- raw$F
  F_inv_T <- nuis$F_inv_T
  coord_dual_iv <- nuis$coord_dual_iv
  beta0_dbl <- raw$beta0
  cert_ctrl <- prob$controls$certificate_control

  enc_star <- link_enclosure$enclosures_star
  enc_0 <- link_enclosure$enclosures_0
  rb <- link_enclosure$root_perturbation_bounds

  Z_ch <- ch$Z_chunk
  w_ch <- ch$weight
  po_ch <- ch$pair_observed
  j_ch <- ch$cohort_index
  obs_pos_ch <- ch$chunk_to_obs_pos
  yl_obs <- ch$y_left_obs
  yr_obs <- ch$y_right_obs
  n_iv <- gmm_iv(as.double(n))

  phi_star <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
  phi_0 <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
  phi_star_via_root <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
  phi_star_intersected <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
  E_root <- numeric(b_len)

  rate_star <- score_star <- treated_star <- control_star <- NULL
  rate_0 <- score_0 <- treated_0 <- control_0 <- NULL
  E_rate <- E_score <- E_treated <- E_control <- NULL
  if (include_components) {
    rate_star <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    score_star <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    treated_star <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    control_star <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    rate_0 <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    score_0 <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    treated_0 <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    control_0 <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    E_rate <- numeric(b_len)
    E_score <- numeric(b_len)
    E_treated <- numeric(b_len)
    E_control <- numeric(b_len)
  }

  a_hi <- cert$a
  a_c_hi <- cert$a_control
  rho_hi <- cert$rho
  x_bar_hi <- nuis$x_bar
  actual_ops <- 0
  physical_ops <- 0
  reused_ops <- 0
  tri_solves <- 0L
  physical_tri_solves <- 0L
  reused_tri_solves <- 0L
  exp_evals <- 0L
  physical_exp_evals <- 0L
  reused_exp_evals <- 0L

  chunk_err <- tryCatch({
    for (u in seq_len(b_len)) {
      if (u > 1L && (u %% 256L) == 0L) ipw_enclosure_check_interrupt()
      ju <- j_ch[u]
      if (is.na(ju)) next # Outside R: all 4 terms and E_root are identically 0; endpoints not read

      wu_iv <- gmm_iv(w_ch[u])
      su <- po_ch[u]
      op <- obs_pos_ch[u]

      # 1. Rate term: 1_R * w_i * (S_i - r_j) / q_j * lambda_{r, j}
      rate_entry <- ipw_if_scratch_get_rate(
        scratch = scratch, u = u, wu_iv = wu_iv, su = su, ju = ju,
        r_cohort = r_cohort, q_cohort = q_cohort, M_cohort = M_cohort, n_iv = n_iv
      )
      s_minus_r <- rate_entry$val$s_minus_r
      xir_iv <- rate_entry$val$xir_iv
      if (rate_entry$reused) {
        reused_ops <- reused_ops + rate_entry$val$ops
      } else {
        physical_ops <- physical_ops + rate_entry$val$ops
      }

      r_0_u <- gmm_iv_mul(xir_iv, enc_0$lambda_r[ju, ])
      r_star_u <- gmm_iv_mul(xir_iv, enc_star$lambda_r[ju, ])
      actual_ops <- actual_ops + 8
      physical_ops <- physical_ops + 2

      # Theorem B rate root error: w_i * |S_i - r_j|_hi / q_j_lo * e_{r, j}
      abs_smr_hi <- gmm_iv_abs(s_minus_r)[2L]
      e_rate_u <- logit_b64_mul(
        logit_b64_div(logit_b64_mul(wu_iv, gmm_iv(abs_smr_hi)), gmm_iv(q_cohort[ju, 1L])),
        gmm_iv(rb$e_r[ju])
      )[2L]
      actual_ops <- actual_ops + 3
      physical_ops <- physical_ops + 3

      sc_0_u <- c(0.0, 0.0)
      sc_star_u <- c(0.0, 0.0)
      tr_0_u <- c(0.0, 0.0)
      tr_star_u <- c(0.0, 0.0)
      cn_0_u <- c(0.0, 0.0)
      cn_star_u <- c(0.0, 0.0)
      e_sc_u <- 0.0
      e_cn_u <- 0.0

      if (su) {
        z_u <- Z_ch[u, ]
        D_u <- if (ju == 1L) 1L else 0L
        # CRITICAL REFUSAL PRECEDENCE: compute Delta_iv FIRST before requesting obs_u from scratch!
        Delta_iv <- gmm_iv_sub(gmm_iv(yr_obs[op]), gmm_iv(yl_obs[op]))
        actual_ops <- actual_ops + 1
        physical_ops <- physical_ops + 1

        obs_entry <- ipw_if_scratch_get_obs(
          scratch = scratch, u = u, a_u = a_seq[u], ju = ju, wu_iv = wu_iv, z_u = z_u,
          M_cohort = M_cohort, Pi_cohort = Pi_cohort, M_over_Pi_cohort = M_over_Pi_cohort,
          F_dbl = F_dbl, F_inv_T = F_inv_T, a_hi = a_hi, a_c_hi = a_c_hi,
          rho_hi = rho_hi, x_bar_hi = x_bar_hi, beta0_dbl = beta0_dbl,
          cert_ctrl = cert_ctrl, p = p
        )
        obs_u <- obs_entry$val
        c_iv <- obs_u$c_iv
        ai_hi <- obs_u$ai_hi
        xi_hi <- obs_u$xi_hi
        one_m_xi_lo <- obs_u$one_m_xi_lo
        o0_iv <- obs_u$o0_iv
        o_star_iv <- obs_u$o_star_iv
        logis0 <- obs_u$logis0
        logis_star <- obs_u$logis_star
        dmp0_iv <- obs_u$dmp0_iv
        dmp_star_iv <- obs_u$dmp_star_iv

        actual_ops <- actual_ops + obs_u$base_ops
        tri_solves <- tri_solves + obs_u$solves
        exp_evals <- exp_evals + obs_u$exp_evals
        if (obs_entry$reused) {
          reused_ops <- reused_ops + obs_u$base_ops
          reused_tri_solves <- reused_tri_solves + obs_u$solves
          reused_exp_evals <- reused_exp_evals + obs_u$exp_evals
        } else {
          physical_ops <- physical_ops + obs_u$base_ops
          physical_tri_solves <- physical_tri_solves + obs_u$solves
          physical_exp_evals <- physical_exp_evals + obs_u$exp_evals
        }

        # 2. Score term: c_i * (D_i - p_i) * (z_i' lambda_beta)
        z_lam_0 <- ipw_enclosure_ball_dot(F_dbl, F_inv_T, z_u, enc_0$lambda_beta_ball, vec_dual_upper = ai_hi)
        z_lam_star <- ipw_enclosure_ball_dot(F_dbl, F_inv_T, z_u, enc_star$lambda_beta_ball, vec_dual_upper = ai_hi)
        sc_0_u <- gmm_iv_mul(gmm_iv_mul(c_iv, dmp0_iv), z_lam_0$iv)
        sc_star_u <- gmm_iv_mul(gmm_iv_mul(c_iv, dmp_star_iv), z_lam_star$iv)
        actual_ops <- actual_ops + z_lam_0$ops + z_lam_star$ops + 4
        physical_ops <- physical_ops + z_lam_0$ops + z_lam_star$ops + 4

        # Theorem B score root error: c_i_hi * a_i * [ (v_i^0_hi * x_i / (1 - x_i)_lo) * Lambda_* + |D_i - p_i^0|_hi * e_lambda ]
        v0_hi <- logis0$var[2L]
        dmp0_abs_hi <- if (D_u == 1L) logis0$one_minus_prob[2L] else logis0$prob[2L]
        t1_sc <- logit_b64_mul(
          logit_b64_div(logit_b64_mul(gmm_iv(v0_hi), gmm_iv(xi_hi)), gmm_iv(one_m_xi_lo)),
          gmm_iv(rb$Lambda_star)
        )[2L]
        t2_sc <- logit_b64_mul(gmm_iv(dmp0_abs_hi), gmm_iv(rb$e_lambda))[2L]
        e_sc_u <- logit_b64_mul(
          logit_b64_mul(gmm_iv(c_iv[2L]), gmm_iv(ai_hi)),
          logit_b64_add(gmm_iv(t1_sc), gmm_iv(t2_sc))
        )[2L]
        actual_ops <- actual_ops + 7
        physical_ops <- physical_ops + 7

        if (D_u == 1L) {
          # 3. Treated term: c_i * (Delta_i - T) / d^+ (no root error)
          tr_u <- gmm_iv_div(gmm_iv_mul(c_iv, gmm_iv_sub(Delta_iv, enc_star$T)), enc_star$d_plus)
          tr_0_u <- tr_u
          tr_star_u <- tr_u
          actual_ops <- actual_ops + 3
          physical_ops <- physical_ops + 3
        } else {
          # 4. Control term: -b_i * (Delta_i - C) / d^- = -n * pi_i * (Delta_i - C)
          b0_u <- gmm_iv_mul(c_iv, o0_iv)
          b_star_u <- gmm_iv_mul(c_iv, o_star_iv)
          res0_u <- gmm_iv_sub(Delta_iv, enc_0$C)
          ress_u <- gmm_iv_sub(Delta_iv, enc_star$C)
          cn_0_u <- gmm_iv_neg(gmm_iv_div(gmm_iv_mul(b0_u, res0_u), enc_0$d_minus))
          cn_star_u <- gmm_iv_neg(gmm_iv_div(gmm_iv_mul(b_star_u, ress_u), enc_star$d_minus))
          actual_ops <- actual_ops + 8
          physical_ops <- physical_ops + 8

          # Theorem B control root error:
          # n * pi_i^0_hi * [ u2(x_bar_c) * min(omega, |Delta_i - C^0|_hi + omega*a_c*rho) + omega*a_c*rho ]
          pi0_hi <- gmm_iv_div(b0_u, enc_0$Sigma_minus)[2L]
          dC_hi <- gmm_iv_abs(res0_u)[2L]
          dC_plus_pt <- logit_b64_add(gmm_iv(dC_hi), gmm_iv(rb$point_mu))[2L]
          min_dc_om <- min(rb$omega, dC_plus_pt)
          inner_cn <- logit_b64_add(
            logit_b64_mul(gmm_iv(rb$u2_x_bar_c), gmm_iv(min_dc_om)),
            gmm_iv(rb$point_mu)
          )[2L]
          e_cn_u <- logit_b64_mul(
            logit_b64_mul(n_iv, gmm_iv(pi0_hi)),
            gmm_iv(inner_cn)
          )[2L]
          actual_ops <- actual_ops + 6
          physical_ops <- physical_ops + 6
        }
      }

      row_0 <- gmm_iv_add(gmm_iv_add(tr_0_u, cn_0_u), gmm_iv_add(sc_0_u, r_0_u))
      row_star <- gmm_iv_add(gmm_iv_add(tr_star_u, cn_star_u), gmm_iv_add(sc_star_u, r_star_u))
      e_tot_u <- logit_b64_add(logit_b64_add(gmm_iv(e_rate_u), gmm_iv(e_sc_u)), gmm_iv(e_cn_u))[2L]
      row_via_root <- ipw_enclosure_widen_iv(row_0, e_tot_u)
      row_inter <- logit_iv_intersect(row_star, row_via_root, "phi_star intersect ([phi_0] +/- E_root)")
      actual_ops <- actual_ops + 10
      physical_ops <- physical_ops + 10

      phi_0[u, ] <- row_0
      phi_star[u, ] <- row_star
      E_root[u] <- e_tot_u
      phi_star_via_root[u, ] <- row_via_root
      phi_star_intersected[u, ] <- row_inter

      if (include_components) {
        rate_0[u, ] <- r_0_u
        score_0[u, ] <- sc_0_u
        treated_0[u, ] <- tr_0_u
        control_0[u, ] <- cn_0_u
        rate_star[u, ] <- r_star_u
        score_star[u, ] <- sc_star_u
        treated_star[u, ] <- tr_star_u
        control_star[u, ] <- cn_star_u
        E_rate[u] <- e_rate_u
        E_score[u] <- e_sc_u
        E_treated[u] <- 0.0
        E_control[u] <- e_cn_u
      }
    }
    NULL
  }, cdid_ipw_enclosure_arithmetic_inconclusive = function(e) e,
     cdid_logit_arithmetic_inconclusive = function(e) e,
     cdid_gmm_accuracy_error = function(e) e)

  if (inherits(chunk_err, "error")) {
    reason <- if (!is.null(chunk_err$reason) && nzchar(chunk_err$reason)) chunk_err$reason else "chunk_if_arithmetic_refusal"
    inc_fields <- list(
      certified = FALSE,
      status = "INCONCLUSIVE",
      inconclusive = TRUE,
      inconclusive_type = "arithmetic",
      inconclusive_reasons = reason,
      storage_mode = if (is_compact) "compact_view" else NULL,
      a_index = a_seq,
      id = ch$id,
      g = ch$g,
      weight = ch$weight,
      arithmetic_refusal = list(
        occurred = TRUE,
        reason = reason,
        message = conditionMessage(chunk_err),
        condition_class = class(chunk_err),
        diagnostics = chunk_err$diagnostics
      ),
      preflight = pre
    )
    if (!is_compact) inc_fields$storage_mode <- NULL
    return(structure(inc_fields, class = "cdid_ipw_link_if_chunk"))
  }

  optional_err <- ipw_enclosure_capture_arithmetic({
    # Optional shared xi decomposition on this chunk (Lemma struct(b))
    shared_xi_out <- NULL
    if (include_shared_xi) {
      xi_chunk <- ctx$shared_xi_chunk
      if (!is.null(xi_chunk)) {
        # Genuinely reuse precomputed shared_xi_chunk without recomputing H^-1 or row primitives!
        actual_ops <- actual_ops + xi_chunk$work$actual_interval_ops
        reused_ops <- reused_ops + xi_chunk$work$actual_interval_ops
        tri_solves <- tri_solves + xi_chunk$work$triangular_solves
        reused_tri_solves <- reused_tri_solves + xi_chunk$work$triangular_solves
        if (!is.null(xi_chunk$work$exp_endpoint_evaluations)) {
          exp_evals <- exp_evals + xi_chunk$work$exp_endpoint_evaluations
          reused_exp_evals <- reused_exp_evals + xi_chunk$work$exp_endpoint_evaluations
        }
      } else {
        ctrl_with_scratch <- ctrl
        attr(ctrl_with_scratch, "chunk_scratch") <- scratch
        xi_chunk <- if (is_compact) {
          enclose_ipw_shared_xi_view_chunk(
            problem = prob,
            a_start = s_idx,
            a_end = e_idx,
            control = ctrl_with_scratch,
            reader_control = ctx$reader_ctrl,
            xi_beta_candidate = xi_beta_candidate,
            H_chol_candidate = R_H_cand,
            require_xi_beta = FALSE
          )
        } else {
          enclose_ipw_shared_xi_chunk(
            problem = prob,
            a_start = s_idx,
            a_end = e_idx,
            control = ctrl_with_scratch,
            xi_beta_candidate = xi_beta_candidate,
            H_chol_candidate = R_H_cand,
            require_xi_beta = FALSE
          )
        }
        if (!isTRUE(xi_chunk$certified)) {
          ipw_enclosure_arith_refuse(
            "Shared xi chunk is INCONCLUSIVE.", reason = "shared_xi_inconclusive",
            diagnostics = list(reasons = xi_chunk$inconclusive_reasons)
          )
        }
        actual_ops <- actual_ops + xi_chunk$work$actual_interval_ops
        physical_ops <- physical_ops + xi_chunk$work$physical_interval_ops
        reused_ops <- reused_ops + xi_chunk$work$reused_shared_interval_ops
        tri_solves <- tri_solves + xi_chunk$work$triangular_solves
        physical_tri_solves <- physical_tri_solves + xi_chunk$work$physical_triangular_solves
        reused_tri_solves <- reused_tri_solves + xi_chunk$work$reused_triangular_solves
        exp_evals <- exp_evals + xi_chunk$work$exp_endpoint_evaluations
        physical_exp_evals <- physical_exp_evals + xi_chunk$work$physical_exp_endpoint_evaluations
        reused_exp_evals <- reused_exp_evals + xi_chunk$work$reused_exp_endpoint_evaluations
      }

      gamma_dot_xi_r_0 <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
      gamma_dot_xi_r_star <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
      lam_dot_rhs_0 <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
      lam_dot_rhs_star <- matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
      neg_v_dot_xi_beta_0 <- if (xi_chunk$has_xi_beta) matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi"))) else NULL
      neg_v_dot_xi_beta_star <- if (xi_chunk$has_xi_beta) matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi"))) else NULL
      phi_shared_0 <- if (xi_chunk$has_xi_beta) matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi"))) else NULL
      phi_shared_star <- if (xi_chunk$has_xi_beta) matrix(0.0, nrow = b_len, ncol = 2L, dimnames = list(NULL, c("lo", "hi"))) else NULL

      for (u in seq_len(b_len)) {
        ju <- j_ch[u]
        if (is.na(ju)) next
        xir_u <- xi_chunk$xi_rate[u, ju, ]
        gamma_dot_xi_r_0[u, ] <- gmm_iv_mul(enc_0$gamma[ju, ], xir_u)
        gamma_dot_xi_r_star[u, ] <- gmm_iv_mul(enc_star$gamma[ju, ], xir_u)

        bd0 <- ipw_enclosure_ball_dot(F_dbl, F_inv_T, matrix(xi_chunk$rhs_0[u, , ], nrow = p, ncol = 2L), enc_0$lambda_beta_ball)
        bds <- ipw_enclosure_ball_dot(F_dbl, F_inv_T, matrix(xi_chunk$rhs_star[u, , ], nrow = p, ncol = 2L), enc_star$lambda_beta_ball)
        lam_dot_rhs_0[u, ] <- bd0$iv
        lam_dot_rhs_star[u, ] <- bds$iv
        actual_ops <- actual_ops + bd0$ops + bds$ops + 2
        physical_ops <- physical_ops + bd0$ops + bds$ops + 2

        if (xi_chunk$has_xi_beta) {
          cand_xi <- xi_chunk$xi_beta_center[u, ]
          # -v' xi_beta in (-v)' tilde_xi_beta +/- ||F^{-T} v||_2 * rad_xi_beta
          c0_dot <- c(0.0, 0.0)
          cs_dot <- c(0.0, 0.0)
          for (k in seq_len(p)) {
            if (cand_xi[k] != 0.0) {
              c0_dot <- gmm_iv_add(c0_dot, gmm_iv_mul(enc_0$G_beta[k, ], gmm_iv(cand_xi[k])))
              cs_dot <- gmm_iv_add(cs_dot, gmm_iv_mul(enc_star$G_beta[k, ], gmm_iv(cand_xi[k])))
            }
          }
          r0_act <- logit_b64_mul(gmm_iv(enc_0$v_dual_upper), gmm_iv(xi_chunk$xi_beta_0_radius_m_norm[u]))[2L]
          rs_act <- logit_b64_mul(gmm_iv(enc_star$v_dual_upper), gmm_iv(xi_chunk$xi_beta_star_radius_m_norm[u]))[2L]
          act0 <- ipw_enclosure_widen_iv(c0_dot, r0_act)
          acts <- ipw_enclosure_widen_iv(cs_dot, rs_act)
          neg_v_dot_xi_beta_0[u, ] <- act0
          neg_v_dot_xi_beta_star[u, ] <- acts
          actual_ops <- actual_ops + 4 * p + 6
          physical_ops <- physical_ops + 4 * p + 6

          t_0_u <- if (include_components) gmm_iv_add(treated_0[u, ], control_0[u, ]) else gmm_iv_sub(phi_0[u, ], gmm_iv_add(gamma_dot_xi_r_0[u, ], lam_dot_rhs_0[u, ]))
          t_s_u <- if (include_components) gmm_iv_add(treated_star[u, ], control_star[u, ]) else gmm_iv_sub(phi_star[u, ], gmm_iv_add(gamma_dot_xi_r_star[u, ], lam_dot_rhs_star[u, ]))
          phi_shared_0[u, ] <- gmm_iv_add(t_0_u, gmm_iv_add(gamma_dot_xi_r_0[u, ], act0))
          phi_shared_star[u, ] <- gmm_iv_add(t_s_u, gmm_iv_add(gamma_dot_xi_r_star[u, ], acts))
          actual_ops <- actual_ops + 6
          physical_ops <- physical_ops + 6
        }
      }

      shared_xi_out <- list(
        xi_chunk = xi_chunk,
        gamma_dot_xi_rate_0 = gamma_dot_xi_r_0,
        gamma_dot_xi_rate_star = gamma_dot_xi_r_star,
        lambda_beta_dot_rhs_0 = lam_dot_rhs_0,
        lambda_beta_dot_rhs_star = lam_dot_rhs_star,
        neg_v_dot_xi_beta_0 = neg_v_dot_xi_beta_0,
        neg_v_dot_xi_beta_star = neg_v_dot_xi_beta_star,
        phi_shared_0 = phi_shared_0,
        phi_shared_star = phi_shared_star
      )
    }

    # Optional Corollary sep per-row evaluation and total error budget for stored phi_candidate
    sep_rows <- NULL
    if (has_phi_cand) {
      eval_err <- numeric(b_len)
      total_err <- numeric(b_len)
      direct_star_err <- numeric(b_len)
      for (u in seq_len(b_len)) {
        if (!in_r_sub[u]) {
          eval_err[u] <- 0.0
          total_err[u] <- 0.0
          direct_star_err[u] <- 0.0
        } else {
          ev_u <- ipw_enclosure_iv_dist(phi_0[u, ], phi_cand_sub[u])
          eval_err[u] <- ev_u
          total_err[u] <- logit_b64_add(gmm_iv(ev_u), gmm_iv(E_root[u]))[2L]
          direct_star_err[u] <- ipw_enclosure_iv_dist(phi_star[u, ], phi_cand_sub[u])
          actual_ops <- actual_ops + 5
          physical_ops <- physical_ops + 5
        }
      }
      sep_rows <- list(
        phi_candidate = phi_cand_sub,
        eval_error = eval_err,
        root_error = E_root,
        total_error_bound = total_err,
        direct_star_error = direct_star_err,
        max_eval_error = max(eval_err),
        max_root_error = max(E_root),
        max_total_error_bound = max(total_err)
      )
    }

    NULL
  })
  if (inherits(optional_err, "error")) {
    return(ipw_enclosure_inconclusive_result(
      optional_err, "cdid_ipw_link_if_chunk",
      list(
        storage_mode = if (is_compact) "compact_view" else NULL,
        a_index = a_seq,
        id = ch$id,
        problem_ref = prob$ref,
        preflight = pre
      )
    ))
  }
  stopifnot(actual_ops + ipw_if_scratch_work(scratch) - cache_work_at_entry <= pre$work_ops_bound)
  if (is_compact) {
    stopifnot((actual_ops + ch$reader_work_ops + ipw_if_scratch_work(scratch) - cache_work_at_entry) <= pre$combined_work_ops_bound)
  }

  out_fields <- list(
    certified = TRUE,
    status = "CERTIFIED",
    inconclusive = FALSE,
    storage_mode = if (is_compact) "compact_view" else NULL,
    a_index = a_seq,
    id = ch$id,
    g = ch$g,
    weight = w_ch,
    pair_observed = po_ch,
    cohort_index = j_ch,
    in_risk = in_r_sub,
    in_observed = ch$in_observed,
    phi_star = phi_star,
    phi_0 = phi_0,
    E_root = E_root,
    phi_star_via_root_bound = phi_star_via_root,
    phi_star_intersected = phi_star_intersected,
    components_star = if (include_components) {
      list(
        rate = rate_star,
        score = score_star,
        treated = treated_star,
        control = control_star
      )
    } else NULL,
    components_0 = if (include_components) {
      list(
        rate = rate_0,
        score = score_0,
        treated = treated_0,
        control = control_0
      )
    } else NULL,
    E_root_components = if (include_components) {
      list(
        rate = E_rate,
        score = E_score,
        treated = E_treated,
        control = E_control
      )
    } else NULL,
    separated_row_errors = sep_rows,
    shared_xi = shared_xi_out,
    preflight = pre,
    work = list(
      chunk_rows = b_len,
      triangular_solves = logit_cert_count_output(tri_solves),
      physical_triangular_solves = logit_cert_count_output(physical_tri_solves),
      reused_triangular_solves = logit_cert_count_output(reused_tri_solves),
      exp_endpoint_evaluations = logit_cert_count_output(exp_evals),
      physical_exp_endpoint_evaluations = logit_cert_count_output(physical_exp_evals),
      reused_exp_endpoint_evaluations = logit_cert_count_output(reused_exp_evals),
      actual_interval_ops = actual_ops,
      logical_interval_ops = actual_ops,
      physical_interval_ops = physical_ops,
      reused_shared_interval_ops = reused_ops,
      actual_reader_work_ops = ch$reader_work_ops,
      cache_copy_work_ops_charged = ipw_if_scratch_work(scratch) - cache_work_at_entry,
      combined_work_ops = actual_ops + ch$reader_work_ops + ipw_if_scratch_work(scratch) - cache_work_at_entry,
      endpoints_read = logit_cert_count_output(ch$n_endpoints_read)
    )
  )
  if (!is_compact) out_fields$storage_mode <- NULL
  structure(out_fields, class = "cdid_ipw_link_if_chunk")
}

# Apply a common outcome-free shared_xi_chunk to a link on the SAME nuisance problem
# via Lemma struct(b): phi_{il} = t_{il} + gamma_l' xi_{r,i} - v_l' xi_{beta,i}
apply_ipw_shared_xi_to_link <- function(
  link_enclosure,
  shared_xi_chunk,
  control,
  phi_candidate = NULL
) {
  ipw_enclosure_require_deps()
  if (missing(link_enclosure) || missing(shared_xi_chunk) || missing(control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`link_enclosure`, `shared_xi_chunk`, `control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }
  prob <- ipw_validate_link_enclosure_mode(link_enclosure, require_certified = TRUE, expected_storage_mode = "dense")
  ctrl <- validate_ipw_link_enclosure_control(control)
  if (!inherits(shared_xi_chunk, "cdid_ipw_shared_xi_chunk") || !is.list(shared_xi_chunk) ||
      !ipw_enclosure_seal_ok(shared_xi_chunk, "shared_xi_chunk") ||
      identical(shared_xi_chunk$storage_mode, "compact_view") ||
      !identical(shared_xi_chunk$problem_ref, prob$ref) ||
      !isTRUE(shared_xi_chunk$certified) || !isTRUE(shared_xi_chunk$has_xi_beta)) {
    cdid_ipw_enclosure_abort(
      "`shared_xi_chunk` must be a valid sealed dense `cdid_ipw_shared_xi_chunk` from the same nuisance problem with `has_xi_beta == TRUE`.",
      c("cdid_ipw_enclosure_shared_xi_error", "cdid_ipw_enclosure_integrity_error")
    )
  }

  a_seq <- shared_xi_chunk$a_index
  ctrl_with_reuse <- ctrl
  attr(ctrl_with_reuse, "chunk_scratch") <- shared_xi_chunk$chunk_scratch
  attr(ctrl_with_reuse, "shared_xi_chunk") <- shared_xi_chunk
  enclose_ipw_link_if_chunk(
    link_enclosure = link_enclosure,
    a_start = a_seq[1L],
    a_end = a_seq[length(a_seq)],
    control = ctrl_with_reuse,
    phi_candidate = phi_candidate,
    xi_beta_candidate = shared_xi_chunk$xi_beta_center,
    include_components = TRUE,
    include_shared_xi = TRUE
  )
}

apply_ipw_shared_xi_to_link_view <- function(
  link_enclosure,
  shared_xi_chunk,
  control,
  reader_control,
  phi_candidate = NULL
) {
  ipw_enclosure_require_deps()
  if (missing(link_enclosure) || missing(shared_xi_chunk) ||
      missing(control) || missing(reader_control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`link_enclosure`, `shared_xi_chunk`, `control`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }
  prob <- ipw_validate_link_enclosure_mode(link_enclosure, require_certified = TRUE, expected_storage_mode = "compact_view")
  ctrl <- validate_ipw_link_enclosure_control(control)
  validate_ipw_owner_control(reader_control)
  if (!inherits(shared_xi_chunk, "cdid_ipw_shared_xi_chunk") || !is.list(shared_xi_chunk) ||
      !ipw_enclosure_seal_ok(shared_xi_chunk, "shared_xi_chunk") ||
      !identical(shared_xi_chunk$storage_mode, "compact_view") ||
      !identical(shared_xi_chunk$problem_ref, prob$ref) ||
      !isTRUE(shared_xi_chunk$certified) || !isTRUE(shared_xi_chunk$has_xi_beta)) {
    cdid_ipw_enclosure_abort(
      "`shared_xi_chunk` must be a valid sealed compact `cdid_ipw_shared_xi_chunk` from the same nuisance problem with `has_xi_beta == TRUE`.",
      c("cdid_ipw_enclosure_shared_xi_error", "cdid_ipw_enclosure_integrity_error")
    )
  }

  a_seq <- shared_xi_chunk$a_index
  ctrl_with_reuse <- ctrl
  attr(ctrl_with_reuse, "chunk_scratch") <- shared_xi_chunk$chunk_scratch
  attr(ctrl_with_reuse, "shared_xi_chunk") <- shared_xi_chunk
  enclose_ipw_link_if_view_chunk(
    link_enclosure = link_enclosure,
    a_start = a_seq[1L],
    a_end = a_seq[length(a_seq)],
    control = ctrl_with_reuse,
    reader_control = reader_control,
    phi_candidate = phi_candidate,
    xi_beta_candidate = shared_xi_chunk$xi_beta_center,
    include_components = TRUE,
    include_shared_xi = TRUE
  )
}

enclose_ipw_link_if_chunk <- function(
  link_enclosure,
  a_start,
  a_end,
  control,
  phi_candidate = NULL,
  xi_beta_candidate = NULL,
  H_chol_candidate = NULL,
  include_components = TRUE,
  include_shared_xi = FALSE
) {
  ipw_enclosure_require_deps()
  if (missing(link_enclosure) || missing(a_start) || missing(a_end) || missing(control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`link_enclosure`, `a_start`, `a_end`, `control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }

  prob <- ipw_validate_link_enclosure_mode(link_enclosure, require_certified = FALSE, expected_storage_mode = "dense")
  scratch <- attr(control, "chunk_scratch", exact = TRUE)
  shared_xi_chunk <- attr(control, "shared_xi_chunk", exact = TRUE)
  ctrl <- validate_ipw_link_enclosure_control(control)
  if (!gmm_scalar_logical(include_components) || !gmm_scalar_logical(include_shared_xi)) {
    cdid_ipw_enclosure_abort(
      "`include_components` and `include_shared_xi` must be scalar logicals.",
      "cdid_ipw_enclosure_input_error"
    )
  }

  raw <- prob$raw
  n <- nrow(raw$Z)
  p <- ncol(raw$Z)
  m <- length(raw$cohorts)

  s_idx <- ipw_enclosure_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- ipw_enclosure_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > n) {
    cdid_ipw_enclosure_abort(
      paste0("Chunk range `[a_start, a_end]` must satisfy `1 <= a_start <= a_end <= ", n, "`."),
      "cdid_ipw_enclosure_bound_error"
    )
  }
  b_len <- e_idx - s_idx + 1

  has_phi_cand <- !is.null(phi_candidate)
  if (has_phi_cand) {
    if (is.object(phi_candidate) || !is.atomic(phi_candidate) || !is.null(dim(phi_candidate)) ||
        !is.numeric(phi_candidate) || !(length(phi_candidate) %in% c(b_len, n))) {
      cdid_ipw_enclosure_abort(
        paste0("`phi_candidate` must be a finite atomic numeric vector of length `b_len` (", b_len, ") or `n` (", n, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (!is.null(H_chol_candidate)) {
    if (!is.matrix(H_chol_candidate) || is.object(H_chol_candidate) || !is.numeric(H_chol_candidate) ||
        nrow(H_chol_candidate) != p || ncol(H_chol_candidate) != p) {
      cdid_ipw_enclosure_abort(
        paste0("`H_chol_candidate` must be a numeric p x p matrix (p = ", p, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (!is.null(xi_beta_candidate)) {
    if (!is.matrix(xi_beta_candidate) || is.object(xi_beta_candidate) || !is.numeric(xi_beta_candidate) ||
        ncol(xi_beta_candidate) != p || !(nrow(xi_beta_candidate) %in% c(b_len, n))) {
      cdid_ipw_enclosure_abort(
        paste0("`xi_beta_candidate` must be a finite numeric matrix with ncol = p (", p, ") and nrow = chunk length (", b_len, ") or n (", n, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  has_xi_beta <- isTRUE(include_shared_xi) && (!is.null(xi_beta_candidate) || !is.null(H_chol_candidate) || !is.null(link_enclosure$candidates$H_chol))

  # Preflight BEFORE allocating chunk output matrices or scanning candidate content
  pre <- (if (is.null(scratch)) preflight_ipw_link_if_chunk else preflight_ipw_link_if_chunk_cached)(
    b_len = b_len, p = p, m = m, control = ctrl,
    certificate_control = prob$controls$certificate_control,
    include_components = include_components,
    include_shared_xi = include_shared_xi,
    has_phi_candidate = has_phi_cand,
    has_xi_beta = has_xi_beta
  )

  a_seq <- seq.int(as.integer(s_idx), as.integer(e_idx))
  phi_cand_sub <- NULL
  if (has_phi_cand) {
    phi_cand_sub <- if (length(phi_candidate) == n) as.double(phi_candidate[a_seq]) else as.double(phi_candidate)
    if (anyNA(phi_cand_sub) || any(!is.finite(phi_cand_sub))) {
      cdid_ipw_enclosure_abort("Queried `phi_candidate` rows must be finite.",
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error"))
    }
  }

  R_H_cand <- if (!is.null(H_chol_candidate)) {
    validate_logit_certificate_metric_factor(H_chol_candidate, p)
  } else if (!is.null(link_enclosure$candidates$H_chol)) {
    link_enclosure$candidates$H_chol
  } else NULL

  ctx <- construct_dense_ipw_if_chunk_access_context(
    link_enclosure = link_enclosure,
    prob = prob,
    s_idx = s_idx,
    e_idx = e_idx,
    ctrl = ctrl,
    pre = pre,
    include_components = include_components,
    include_shared_xi = include_shared_xi,
    has_phi_cand = has_phi_cand,
    phi_cand_sub = phi_cand_sub,
    xi_beta_candidate = xi_beta_candidate,
    R_H_cand = R_H_cand,
    scratch = scratch,
    shared_xi_chunk = shared_xi_chunk
  )
  ipw_enclosure_if_chunk_run_with_context(ctx)
}

enclose_ipw_link_if_view_chunk <- function(
  link_enclosure,
  a_start,
  a_end,
  control,
  reader_control,
  phi_candidate = NULL,
  xi_beta_candidate = NULL,
  H_chol_candidate = NULL,
  include_components = TRUE,
  include_shared_xi = FALSE
) {
  ipw_enclosure_require_deps()
  if (missing(link_enclosure) || missing(a_start) || missing(a_end) ||
      missing(control) || missing(reader_control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`link_enclosure`, `a_start`, `a_end`, `control`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }

  prob <- ipw_validate_link_enclosure_mode(link_enclosure, require_certified = FALSE, expected_storage_mode = "compact_view")
  scratch <- attr(control, "chunk_scratch", exact = TRUE)
  shared_xi_chunk <- attr(control, "shared_xi_chunk", exact = TRUE)
  ctrl <- validate_ipw_link_enclosure_control(control)
  expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  trace_env <- attr(reader_control, "trace_env", exact = TRUE)
  reader_ctrl <- validate_ipw_owner_control(reader_control)
  if (!gmm_scalar_logical(include_components) || !gmm_scalar_logical(include_shared_xi)) {
    cdid_ipw_enclosure_abort(
      "`include_components` and `include_shared_xi` must be scalar logicals.",
      "cdid_ipw_enclosure_input_error"
    )
  }

  raw <- prob$raw
  n <- as.integer(raw$n_a_units)
  p <- as.integer(raw$p_covariates)
  m <- length(raw$cohorts)
  panel_owner <- get("panel_owner", envir = prob$panel_ref, inherits = FALSE)

  s_idx <- ipw_enclosure_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- ipw_enclosure_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > n) {
    cdid_ipw_enclosure_abort(
      paste0("Chunk range `[a_start, a_end]` must satisfy `1 <= a_start <= a_end <= ", n, "`."),
      "cdid_ipw_enclosure_bound_error"
    )
  }
  b_len <- e_idx - s_idx + 1

  ipw_owner_refuse_decorative_attrs(phi_candidate, "phi_candidate")
  ipw_owner_refuse_decorative_attrs(xi_beta_candidate, "xi_beta_candidate")
  ipw_owner_refuse_decorative_attrs(H_chol_candidate, "H_chol_candidate")
  has_phi_cand <- !is.null(phi_candidate)
  caller_phi_entries <- 0
  if (has_phi_cand) {
    if (is.object(phi_candidate) || !is.atomic(phi_candidate) || !is.null(dim(phi_candidate)) ||
        !is.numeric(phi_candidate) || !(length(phi_candidate) %in% c(b_len, n))) {
      cdid_ipw_enclosure_abort(
        paste0("`phi_candidate` must be a finite atomic numeric vector of length `b_len` (", b_len, ") or `n` (", n, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
    caller_phi_entries <- length(phi_candidate)
  }
  if (!is.null(H_chol_candidate)) {
    if (!is.matrix(H_chol_candidate) || is.object(H_chol_candidate) || !is.numeric(H_chol_candidate) ||
        nrow(H_chol_candidate) != p || ncol(H_chol_candidate) != p) {
      cdid_ipw_enclosure_abort(
        paste0("`H_chol_candidate` must be a numeric p x p matrix (p = ", p, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  caller_xi_entries <- 0
  if (!is.null(xi_beta_candidate)) {
    if (!is.matrix(xi_beta_candidate) || is.object(xi_beta_candidate) || !is.numeric(xi_beta_candidate) ||
        ncol(xi_beta_candidate) != p || !(nrow(xi_beta_candidate) %in% c(b_len, n))) {
      cdid_ipw_enclosure_abort(
        paste0("`xi_beta_candidate` must be a finite numeric matrix with ncol = p (", p, ") and nrow = chunk length (", b_len, ") or n (", n, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
    caller_xi_entries <- nrow(xi_beta_candidate) * p
  }
  if (!is.null(expected_a_id)) {
    ipw_owner_refuse_decorative_attrs(expected_a_id, "expected_a_id")
    if (is.factor(expected_a_id) || is.object(expected_a_id) || !is.atomic(expected_a_id) ||
        !is.null(dim(expected_a_id)) || !(length(expected_a_id) %in% c(b_len, n))) {
      cdid_ipw_owner_abort(
        "`expected_a_id` must be NULL or a plain atomic vector of length `b_len` or `n_a_units`.",
        c("cdid_ipw_owner_id_error", "cdid_ipw_owner_order_mismatch_error")
      )
    }
  }

  has_xi_beta <- isTRUE(include_shared_xi) && (!is.null(xi_beta_candidate) || !is.null(H_chol_candidate) || !is.null(link_enclosure$candidates$H_chol))

  # Combined arithmetic + reader preflight BEFORE allocating chunk output matrices or scanning candidate content
  pre <- preflight_ipw_link_if_view_chunk(
    b_len = b_len,
    p = p,
    m = m,
    control = ctrl,
    reader_control = reader_ctrl,
    certificate_control = prob$controls$certificate_control,
    include_components = include_components,
    include_shared_xi = include_shared_xi,
    has_phi_candidate = has_phi_cand,
    has_xi_beta = has_xi_beta,
    read_outcomes = TRUE,
    panel_owner = panel_owner,
    a_seq = seq.int(as.integer(s_idx), as.integer(e_idx)),
    expected_a_id = expected_a_id,
    caller_phi_entries = caller_phi_entries,
    caller_xi_entries = caller_xi_entries
  )

  a_seq <- seq.int(as.integer(s_idx), as.integer(e_idx))

  phi_cand_sub <- NULL
  if (has_phi_cand) {
    phi_cand_sub <- if (length(phi_candidate) == n) as.double(phi_candidate[a_seq]) else as.double(phi_candidate)
    if (anyNA(phi_cand_sub) || any(!is.finite(phi_cand_sub))) {
      cdid_ipw_enclosure_abort("Queried `phi_candidate` rows must be finite.",
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error"))
    }
  }

  R_H_cand <- if (!is.null(H_chol_candidate)) {
    validate_logit_certificate_metric_factor(H_chol_candidate, p)
  } else if (!is.null(link_enclosure$candidates$H_chol)) {
    link_enclosure$candidates$H_chol
  } else NULL

  ctx <- construct_view_ipw_if_chunk_access_context(
    link_enclosure = link_enclosure,
    prob = prob,
    s_idx = s_idx,
    e_idx = e_idx,
    ctrl = ctrl,
    reader_ctrl = reader_ctrl,
    pre = pre,
    include_components = include_components,
    include_shared_xi = include_shared_xi,
    has_phi_cand = has_phi_cand,
    phi_cand_sub = phi_cand_sub,
    xi_beta_candidate = xi_beta_candidate,
    R_H_cand = R_H_cand,
    expected_a_id = expected_a_id,
    trace_env = trace_env,
    scratch = scratch,
    shared_xi_chunk = shared_xi_chunk
  )
  ipw_enclosure_if_chunk_run_with_context(ctx)
}

enclose_ipw_link_if_chunk_with_scratch <- function(
  link_enclosure,
  a_start,
  a_end,
  control,
  scratch,
  phi_candidate = NULL,
  xi_beta_candidate = NULL,
  H_chol_candidate = NULL,
  include_components = TRUE,
  include_shared_xi = FALSE
) {
  ctrl_with_scratch <- validate_ipw_link_enclosure_control(control)
  attr(ctrl_with_scratch, "chunk_scratch") <- scratch
  enclose_ipw_link_if_chunk(
    link_enclosure = link_enclosure,
    a_start = a_start,
    a_end = a_end,
    control = ctrl_with_scratch,
    phi_candidate = phi_candidate,
    xi_beta_candidate = xi_beta_candidate,
    H_chol_candidate = H_chol_candidate,
    include_components = include_components,
    include_shared_xi = include_shared_xi
  )
}

enclose_ipw_link_if_view_chunk_with_scratch <- function(
  link_enclosure,
  a_start,
  a_end,
  control,
  reader_control,
  scratch,
  phi_candidate = NULL,
  xi_beta_candidate = NULL,
  H_chol_candidate = NULL,
  include_components = TRUE,
  include_shared_xi = FALSE
) {
  ctrl_with_scratch <- validate_ipw_link_enclosure_control(control)
  attr(ctrl_with_scratch, "chunk_scratch") <- scratch
  enclose_ipw_link_if_view_chunk(
    link_enclosure = link_enclosure,
    a_start = a_start,
    a_end = a_end,
    control = ctrl_with_scratch,
    reader_control = reader_control,
    phi_candidate = phi_candidate,
    xi_beta_candidate = xi_beta_candidate,
    H_chol_candidate = H_chol_candidate,
    include_components = include_components,
    include_shared_xi = include_shared_xi
  )
}

stream_ipw_link_if_enclosure <- function(
  link_enclosure,
  control,
  callback,
  phi_candidate = NULL,
  xi_beta_candidate = NULL,
  H_chol_candidate = NULL,
  include_components = TRUE,
  include_shared_xi = FALSE
) {
  ipw_enclosure_require_deps()
  if (missing(link_enclosure) || missing(control) || missing(callback)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`link_enclosure`, `control`, `callback`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }
  prob <- ipw_validate_link_enclosure_mode(link_enclosure, require_certified = TRUE, expected_storage_mode = "dense")
  ctrl <- validate_ipw_link_enclosure_control(control)
  if (!is.function(callback)) {
    cdid_ipw_enclosure_abort("`callback` must be a function.", "cdid_ipw_enclosure_input_error")
  }
  n <- nrow(prob$raw$Z)
  step <- min(as.double(n), ctrl$max_chunk_rows)

  if (!is.null(phi_candidate) && length(phi_candidate) != n) {
    cdid_ipw_enclosure_abort(
      paste0("`phi_candidate` in `stream_ipw_link_if_enclosure` must have length `n_a_units` (", n, ")."),
      "cdid_ipw_enclosure_candidate_error"
    )
  }
  if (!is.null(xi_beta_candidate) && (!is.matrix(xi_beta_candidate) || nrow(xi_beta_candidate) != n)) {
    cdid_ipw_enclosure_abort(
      paste0("`xi_beta_candidate` in `stream_ipw_link_if_enclosure` must be a matrix with `nrow == n_a_units` (", n, ")."),
      "cdid_ipw_enclosure_candidate_error"
    )
  }

  # Whole traversal work is bounded before the first candidate slice/callback.
  # Chunk output is released after each callback; caller-retained callback results
  # belong to the caller and are not a process-RSS guarantee.
  stream_pre <- preflight_ipw_link_if_chunk(step, ncol(prob$raw$Z),
    length(prob$raw$cohorts), ctrl, prob$controls$certificate_control,
    include_components = include_components, include_shared_xi = include_shared_xi,
    has_phi_candidate = !is.null(phi_candidate),
    has_xi_beta = isTRUE(include_shared_xi) && (!is.null(xi_beta_candidate) ||
      !is.null(H_chol_candidate) || !is.null(link_enclosure$candidates$H_chol)))
  stream_work_bound <- ipw_enclosure_safe_product(n,
    stream_pre$work_ops_bound / step + 2, ctrl$max_work_ops, "whole IF traversal work")

  emitted_units <- 0L
  batches <- 0L
  max_batch_units <- 0L
  total_ops <- 0
  phi_star_sum_iv <- c(0.0, 0.0)
  phi_0_sum_iv <- c(0.0, 0.0)
  max_E_root <- 0.0
  max_eval_err <- if (!is.null(phi_candidate)) 0.0 else NA_real_
  max_total_err <- if (!is.null(phi_candidate)) 0.0 else NA_real_

  a_start <- 1L
  while (a_start <= n) {
    ipw_enclosure_check_interrupt()
    a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
    b_len <- a_end - a_start + 1L

    ch <- enclose_ipw_link_if_chunk(
      link_enclosure = link_enclosure,
      a_start = a_start,
      a_end = a_end,
      control = ctrl,
      phi_candidate = phi_candidate,
      xi_beta_candidate = xi_beta_candidate,
      H_chol_candidate = H_chol_candidate,
      include_components = include_components,
      include_shared_xi = include_shared_xi
    )
    if (!isTRUE(ch$certified)) {
      cdid_ipw_enclosure_abort(
        "Streamed IF chunk became INCONCLUSIVE during traversal.",
        "cdid_ipw_enclosure_inconclusive",
        diagnostics = list(batch_start = a_start, batch_end = a_end, reasons = ch$inconclusive_reasons)
      )
    }

    sum_err <- ipw_enclosure_capture_arithmetic({
      for (u in seq_len(b_len)) {
        phi_star_sum_iv <- gmm_iv_add(phi_star_sum_iv, ch$phi_star[u, ])
        phi_0_sum_iv <- gmm_iv_add(phi_0_sum_iv, ch$phi_0[u, ])
        total_ops <- total_ops + 2
      }
      NULL
    })
    if (inherits(sum_err, "error")) cdid_ipw_enclosure_abort(
      "Stream aggregate arithmetic is INCONCLUSIVE.", "cdid_ipw_enclosure_inconclusive",
      diagnostics = list(batch_start = a_start, batch_end = a_end, arithmetic = conditionMessage(sum_err)))
    max_E_root <- max(max_E_root, max(ch$E_root))
    if (!is.null(ch$separated_row_errors)) {
      max_eval_err <- max(max_eval_err, ch$separated_row_errors$max_eval_error)
      max_total_err <- max(max_total_err, ch$separated_row_errors$max_total_error_bound)
    }

    batches <- batches + 1L
    ch$batch_index <- batches
    callback(ch)

    emitted_units <- emitted_units + b_len
    max_batch_units <- max(max_batch_units, b_len)
    total_ops <- total_ops + ch$work$actual_interval_ops
    rm(ch)
    a_start <- as.double(a_end) + 1
  }

  if (emitted_units != n) {
    cdid_ipw_enclosure_abort("Streamed IF enclosure unit count mismatch.", "cdid_ipw_enclosure_error")
  }

  invisible(list(
    n_a_units = logit_cert_count_output(emitted_units),
    n_batches = logit_cert_count_output(batches),
    max_batch_units = logit_cert_count_output(max_batch_units),
    phi_star_sum_iv = phi_star_sum_iv,
    phi_star_sum_contains_zero = (phi_star_sum_iv[1L] <= 0.0 && phi_star_sum_iv[2L] >= 0.0),
    phi_0_sum_iv = phi_0_sum_iv,
    max_E_root = max_E_root,
    max_eval_error = max_eval_err,
    max_total_error_bound = max_total_err,
    total_chunk_interval_ops = total_ops,
    whole_stream_work_bound = stream_work_bound
  ))
}

stream_ipw_link_if_view_enclosure <- function(
  link_enclosure,
  control,
  callback,
  reader_control,
  phi_candidate = NULL,
  xi_beta_candidate = NULL,
  H_chol_candidate = NULL,
  include_components = TRUE,
  include_shared_xi = FALSE
) {
  ipw_enclosure_require_deps()
  if (missing(link_enclosure) || missing(control) || missing(callback) || missing(reader_control)) {
    cdid_ipw_enclosure_abort(
      "All required arguments (`link_enclosure`, `control`, `callback`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_enclosure_input_error"
    )
  }
  prob <- ipw_validate_link_enclosure_mode(link_enclosure, require_certified = TRUE, expected_storage_mode = "compact_view")
  ctrl <- validate_ipw_link_enclosure_control(control)
  expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  trace_env <- attr(reader_control, "trace_env", exact = TRUE)
  reader_ctrl <- validate_ipw_owner_control(reader_control)
  if (!is.function(callback)) {
    cdid_ipw_enclosure_abort("`callback` must be a function.", "cdid_ipw_enclosure_input_error")
  }
  if (!gmm_scalar_logical(include_components) || !gmm_scalar_logical(include_shared_xi)) {
    cdid_ipw_enclosure_abort(
      "`include_components` and `include_shared_xi` must be scalar logicals.",
      "cdid_ipw_enclosure_input_error"
    )
  }

  raw <- prob$raw
  n <- as.integer(raw$n_a_units)
  p <- as.integer(raw$p_covariates)
  m <- length(raw$cohorts)
  panel_owner <- get("panel_owner", envir = prob$panel_ref, inherits = FALSE)
  step <- min(as.double(n), ctrl$max_chunk_rows, reader_ctrl$max_chunk_rows)
  n_chunks <- ceiling(n / step)

  ipw_owner_refuse_decorative_attrs(phi_candidate, "phi_candidate")
  ipw_owner_refuse_decorative_attrs(xi_beta_candidate, "xi_beta_candidate")
  ipw_owner_refuse_decorative_attrs(H_chol_candidate, "H_chol_candidate")
  if (!is.null(phi_candidate)) {
    if (is.object(phi_candidate) || !is.atomic(phi_candidate) || !is.null(dim(phi_candidate)) ||
        !is.numeric(phi_candidate) || length(phi_candidate) != n) {
      cdid_ipw_enclosure_abort(
        paste0("`phi_candidate` in `stream_ipw_link_if_view_enclosure` must have length `n_a_units` (", n, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (!is.null(xi_beta_candidate)) {
    if (!is.matrix(xi_beta_candidate) || is.object(xi_beta_candidate) || !is.numeric(xi_beta_candidate) ||
        nrow(xi_beta_candidate) != n || ncol(xi_beta_candidate) != p) {
      cdid_ipw_enclosure_abort(
        paste0("`xi_beta_candidate` in `stream_ipw_link_if_view_enclosure` must be a numeric matrix with `nrow == n_a_units` (", n, ") and `ncol == p` (", p, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }
  if (!is.null(H_chol_candidate)) {
    if (!is.matrix(H_chol_candidate) || is.object(H_chol_candidate) || !is.numeric(H_chol_candidate) ||
        nrow(H_chol_candidate) != p || ncol(H_chol_candidate) != p) {
      cdid_ipw_enclosure_abort(
        paste0("`H_chol_candidate` must be a numeric p x p matrix (p = ", p, ")."),
        c("cdid_ipw_enclosure_candidate_error", "cdid_ipw_enclosure_input_error")
      )
    }
  }

  # Whole traversal work and concurrent payload are preflighted before the first slice or callback.
  has_xi_beta <- isTRUE(include_shared_xi) && (!is.null(xi_beta_candidate) ||
    !is.null(H_chol_candidate) || !is.null(link_enclosure$candidates$H_chol))

  stream_pre <- preflight_ipw_link_if_view_chunk(
    b_len = step,
    p = p,
    m = m,
    control = ctrl,
    reader_control = reader_ctrl,
    certificate_control = prob$controls$certificate_control,
    include_components = include_components,
    include_shared_xi = include_shared_xi,
    has_phi_candidate = !is.null(phi_candidate),
    has_xi_beta = has_xi_beta,
    read_outcomes = TRUE,
    panel_owner = NULL,
    a_seq = NULL,
    expected_a_id = NULL,
    caller_phi_entries = if (!is.null(phi_candidate)) n else 0,
    caller_xi_entries = if (!is.null(xi_beta_candidate)) n * p else 0
  )
  base_reader_work <- ipw_owner_safe_product(n_chunks, stream_pre$reader_base_work_ops_bound, reader_ctrl$max_work_ops, "base whole-stream reader work")
  base_arith_work <- ipw_enclosure_safe_product(n, stream_pre$work_ops_bound / step + 2, ctrl$max_work_ops, "base whole-stream arithmetic/cache work")
  ipw_enclosure_safe_sum(c(base_reader_work, base_arith_work), ctrl$max_work_ops, "base whole-stream combined work before ID metadata")
  id_info <- ipw_enclosure_id_bytes_for_pass(panel_owner, expected_a_id, n, step, reader_ctrl)
  if (id_info$max_chunk_source_bytes > 0) {
    r_pre_id <- ipw_owner_preflight_id_bytes(
      stream_pre$reader_preflight, reader_ctrl,
      list(total = id_info$max_chunk_source_bytes, max = id_info$max_single_id_bytes),
      owner_max_bytes = 0, lookup_count = 0, search_width = 1, simultaneous = TRUE
    )
    ipw_enclosure_safe_sum(
      c(stream_pre$concurrent_entries, r_pre_id$concurrent_entries),
      max_bound = ctrl$max_matrix_entries,
      label = "whole-stream simultaneous chunk IF arithmetic and character-ID reader payload entries"
    )
  }

  stream_reader_work_bound <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, stream_pre$reader_base_work_ops_bound, reader_ctrl$max_work_ops, "whole IF traversal reader chunk work"),
      id_info$per_pass_additional_work_ops
    ),
    max_bound = reader_ctrl$max_work_ops,
    label = "whole IF traversal reader work operations"
  )
  stream_arith_work_bound <- ipw_enclosure_safe_product(
    n,
    stream_pre$work_ops_bound / step + 2,
    ctrl$max_work_ops,
    "whole IF traversal interval work"
  )
  stream_combined_work_bound <- ipw_enclosure_safe_sum(
    c(stream_arith_work_bound, stream_reader_work_bound),
    max_bound = ctrl$max_work_ops,
    label = "combined whole IF traversal interval and reader work"
  )

  emitted_units <- 0L
  batches <- 0L
  max_batch_units <- 0L
  total_ops <- 0
  total_physical_ops <- 0
  total_reused_ops <- 0
  total_reader_ops <- 0
  total_cache_work_ops <- 0
  total_endpoints_read <- 0
  phi_star_sum_iv <- c(0.0, 0.0)
  phi_0_sum_iv <- c(0.0, 0.0)
  max_E_root <- 0.0
  max_eval_err <- if (!is.null(phi_candidate)) 0.0 else NA_real_
  max_total_err <- if (!is.null(phi_candidate)) 0.0 else NA_real_

  a_start <- 1L
  while (a_start <= n) {
    ipw_enclosure_check_interrupt()
    a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
    b_len <- a_end - a_start + 1L

    r_ctrl_call <- reader_ctrl
    if (!is.null(expected_a_id)) attr(r_ctrl_call, "expected_a_id") <- expected_a_id
    if (!is.null(trace_env)) attr(r_ctrl_call, "trace_env") <- trace_env

    ch <- enclose_ipw_link_if_view_chunk(
      link_enclosure = link_enclosure,
      a_start = a_start,
      a_end = a_end,
      control = ctrl,
      reader_control = r_ctrl_call,
      phi_candidate = phi_candidate,
      xi_beta_candidate = xi_beta_candidate,
      H_chol_candidate = H_chol_candidate,
      include_components = include_components,
      include_shared_xi = include_shared_xi
    )
    if (!isTRUE(ch$certified)) {
      cdid_ipw_enclosure_abort(
        "Streamed IF view chunk became INCONCLUSIVE during traversal.",
        "cdid_ipw_enclosure_inconclusive",
        diagnostics = list(batch_start = a_start, batch_end = a_end, reasons = ch$inconclusive_reasons)
      )
    }

    sum_err <- ipw_enclosure_capture_arithmetic({
      for (u in seq_len(b_len)) {
        phi_star_sum_iv <- gmm_iv_add(phi_star_sum_iv, ch$phi_star[u, ])
        phi_0_sum_iv <- gmm_iv_add(phi_0_sum_iv, ch$phi_0[u, ])
        total_ops <- total_ops + 2
        total_physical_ops <- total_physical_ops + 2
      }
      NULL
    })
    if (inherits(sum_err, "error")) cdid_ipw_enclosure_abort(
      "Stream aggregate arithmetic is INCONCLUSIVE.", "cdid_ipw_enclosure_inconclusive",
      diagnostics = list(batch_start = a_start, batch_end = a_end, arithmetic = conditionMessage(sum_err)))
    max_E_root <- max(max_E_root, max(ch$E_root))
    if (!is.null(ch$separated_row_errors)) {
      max_eval_err <- max(max_eval_err, ch$separated_row_errors$max_eval_error)
      max_total_err <- max(max_total_err, ch$separated_row_errors$max_total_error_bound)
    }

    batches <- batches + 1L
    ch$batch_index <- batches
    callback(ch)

    emitted_units <- emitted_units + b_len
    max_batch_units <- max(max_batch_units, b_len)
    total_ops <- total_ops + ch$work$actual_interval_ops
    total_physical_ops <- total_physical_ops + ch$work$physical_interval_ops
    total_reused_ops <- total_reused_ops + ch$work$reused_shared_interval_ops
    total_reader_ops <- total_reader_ops + ch$work$actual_reader_work_ops
    total_cache_work_ops <- total_cache_work_ops + ch$work$cache_copy_work_ops_charged
    total_endpoints_read <- total_endpoints_read + ch$work$endpoints_read
    rm(ch)
    a_start <- as.double(a_end) + 1
  }

  if (emitted_units != n) {
    cdid_ipw_enclosure_abort("Streamed IF enclosure unit count mismatch.", "cdid_ipw_enclosure_error")
  }

  invisible(list(
    storage_mode = "compact_view",
    n_a_units = logit_cert_count_output(emitted_units),
    n_batches = logit_cert_count_output(batches),
    max_batch_units = logit_cert_count_output(max_batch_units),
    phi_star_sum_iv = phi_star_sum_iv,
    phi_star_sum_contains_zero = (phi_star_sum_iv[1L] <= 0.0 && phi_star_sum_iv[2L] >= 0.0),
    phi_0_sum_iv = phi_0_sum_iv,
    max_E_root = max_E_root,
    max_eval_error = max_eval_err,
    max_total_error_bound = max_total_err,
    total_chunk_interval_ops = total_ops,
    total_physical_interval_ops = total_physical_ops,
    total_reused_shared_interval_ops = total_reused_ops,
    total_reader_work_ops = total_reader_ops,
    total_cache_copy_work_ops_charged = total_cache_work_ops,
    total_combined_work_ops = total_ops + total_reader_ops + total_cache_work_ops,
    endpoints_read = logit_cert_count_output(total_endpoints_read),
    read_diagnostics = list(
      read_outcomes = TRUE,
      endpoints_read = logit_cert_count_output(total_endpoints_read),
      n_endpoints_read = logit_cert_count_output(total_endpoints_read),
      risk_filtered_before_outcome_read = TRUE,
      off_risk_endpoints_read = FALSE,
      nonobserver_endpoints_read = FALSE
    ),
    whole_stream_work_bound = stream_arith_work_bound,
    whole_stream_reader_work_bound = stream_reader_work_bound,
    whole_stream_combined_work_bound = stream_combined_work_bound
  ))
}

traverse_ipw_link_if_enclosure <- stream_ipw_link_if_enclosure
traverse_ipw_link_if_view_enclosure <- stream_ipw_link_if_view_enclosure

validate_ipw_link_enclosure <- function(link_enclosure, problem = NULL, require_certified = FALSE) {
  ipw_validate_link_enclosure_mode(link_enclosure, problem, require_certified)
}
