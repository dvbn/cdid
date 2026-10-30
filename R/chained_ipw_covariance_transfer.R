# Private CP02 complete-IF covariance-accumulation and root-perturbation PD-transfer layer.
#
# Implements the conditional numerical covariance propagation and streaming Gram
# enclosure primitives of:
#   * docs/theory/cp02_ipw_gmm_numerical_propagation.tex
#     (Section 6 Lemmas psd, quad & Proposition cov)
#   * docs/theory/cp02_cluster_inference.tex
#     (Section 1 Lemma samp; Section 3 Proposition ker; Section 4 Theorem cov;
#      Section 5 Theorem gmm & Corollary d11; Section 8 Proposition num)
#
# Scope and boundaries:
#   * Unexported entry point:
#       certify_ipw_covariance_transfer(links, candidate_rows, frame, C, op, V_candidate, control)
#     All seven arguments are required; no internal defaults, public API choices,
#     moment-catalog selection (D10), cluster-vs-unit GMM policy adoption (D11),
#     or public acceptance budget levels (D12) are chosen here.
#   * First scope: every sealed 4bd link enclosure (`cdid_ipw_link_enclosure`) and
#     its sealed nuisance problem (`cdid_ipw_nuisance_problem`) declare the SAME
#     complete A-unit source set, typed IDs/order, case weights `w`, cohort vector `G`,
#     and count normalizer `n_A`. Distinct links may reference either one shared
#     sealed nuisance problem or distinct sealed nuisance problems on that identical
#     complete A-unit frame.
#   * `frame` declares an explicit `"identity"` embedding into the identical complete
#     source A-frame (`normalizer_n == n_A`, identical typed IDs in unit order).
#     No absent rows are silently padded with zeros, no A-normalized IF rows are
#     silently rescaled by `n_F / n_A` or weight-mass ratios, and case weights `w`
#     already inside `Phi` are never multiplied again. Larger-frame transport,
#     frame-only share IF, and estimated aggregation-weight IF remain required
#     future pipeline obligations.
#   * `candidate_rows` supplies bounded complete-IF point blocks (`phi` of dimension
#     `b x L`), unit IDs (`id`), and ordered link IDs (`link_ids`) for requested
#     contiguous source unit ranges `[a_start, a_end]`. Numerical provenance comes
#     from passing those exact stored candidate columns into `enclose_ipw_link_if_chunk()`
#     on each sealed link to obtain separated evaluation discrepancy (`sep_rows$eval_error`
#     from `phi(beta0)`) and root perturbation (`sep_rows$root_error` from `beta0` to
#     `beta*`). `direct_star_error` is never relabeled as a root-only radius.
#   * `C` (`d x L`) and `op` are fixed along nuisance root balls. Under D11, `"GMM"`
#     purpose requires unit aggregation (`row_level = "unit"`, `A = I_n`, `omega = 1`);
#     `"inference"` purpose may declare either unit identity or complete cluster
#     sums (`row_level = "cluster"`, keeping every declared cluster in `cluster_order`
#     including zero-IF clusters).
#   * Implements the reviewed streaming Gram route using only existing native and
#     scalar outward interval primitives (`gmm_iv_mat_residual`, `gmm_iv_dot`,
#     `gmm_iv_norm`, `gmm_iv_add`, `gmm_iv_sub`, `gmm_iv_mul`, `gmm_iv_div`,
#     `gmm_iv_sqrt`, `certify_symmetric_eigen_bounds`).
#   * Keeps `epsilon_acc` distinct from `nu_eval` and `nu_root`, bounds the baseline
#     covariance norm with `v_upper = outward_norm(V_candidate) + epsilon_acc`,
#     retains the quadratic evaluation/root cross term `2 * nu_eval * nu_root` and
#     quadratic terms `nu_eval^2 + 2 * nu_eval * nu_root + nu_root^2`, and transfers
#     positive definiteness only when `epsilon_V < lambda_min_cert(V_candidate)`.
#     Failed sufficient conditions or uncertified/arithmetic failures return explicit
#     `INCONCLUSIVE`, never singularity/non-identification claims or covariance substitution.

cdid_ipw_covariance_transfer_abort <- function(message, class = character(), diagnostics = list()) {
  stop(structure(
    class = unique(c(
      class,
      "cdid_ipw_covariance_transfer_refusal",
      "cdid_ipw_covariance_transfer_error",
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

ipw_cov_transfer_arith_refuse <- function(message, reason, diagnostics = list()) {
  stop(structure(
    class = unique(c(
      "cdid_ipw_covariance_transfer_arithmetic_inconclusive",
      "cdid_ipw_enclosure_arithmetic_inconclusive",
      "cdid_logit_arithmetic_inconclusive",
      "cdid_ipw_covariance_transfer_inconclusive",
      "cdid_ipw_enclosure_inconclusive",
      "cdid_ipw_covariance_transfer_refusal",
      "cdid_ipw_covariance_transfer_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, reason = reason, call = NULL, diagnostics = diagnostics)
  ))
}

ipw_cov_transfer_require_deps <- function() {
  ns_env <- environment(ipw_cov_transfer_require_deps)
  req_fns <- c(
    "gmm_next", "gmm_iv", "gmm_iv_check", "gmm_iv_add", "gmm_iv_neg",
    "gmm_iv_sub", "gmm_iv_mul", "gmm_iv_div", "gmm_iv_sqrt", "gmm_iv_abs",
    "gmm_iv_dot", "gmm_iv_norm", "gmm_iv_mat_residual",
    "certify_symmetric_eigen_bounds", "gmm_sexp_aliases",
    "gmm_scalar_logical", "gmm_scalar_numeric", "gmm_scalar_character",
    "logit_cert_count_output", "logit_b64_add", "logit_b64_sub",
    "logit_b64_mul", "logit_b64_div",
    "ipw_enclosure_int_bound", "ipw_enclosure_safe_product",
    "ipw_enclosure_safe_sum", "ipw_enclosure_seal", "ipw_enclosure_seal_ok",
    "ipw_enclosure_iv_dist", "validate_ipw_link_enclosure_control",
    "preflight_ipw_link_if_chunk", "validate_ipw_link_enclosure",
    "enclose_ipw_link_if_chunk"
  )
  for (fn in req_fns) {
    if (!exists(fn, envir = ns_env, mode = "function", inherits = TRUE)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Required private helper `", fn, "` is not loaded in the namespace."),
        "cdid_ipw_covariance_transfer_dependency_error"
      )
    }
  }
  invisible(TRUE)
}

ipw_cov_transfer_check_interrupt <- function() {
  Sys.sleep(0)
}

ipw_cov_transfer_capture_arithmetic <- function(expr) {
  tryCatch(
    expr,
    cdid_ipw_covariance_transfer_arithmetic_inconclusive = function(e) e,
    cdid_ipw_enclosure_arithmetic_inconclusive = function(e) e,
    cdid_logit_arithmetic_inconclusive = function(e) e,
    cdid_gmm_accuracy_error = function(e) e
  )
}

ipw_cov_transfer_int_bound <- function(x, name, min_val = 1) {
  if (!gmm_scalar_numeric(x) || x != floor(x) || x < min_val || x > 2^53 - 1) {
    cdid_ipw_covariance_transfer_abort(
      paste0("`", name, "` must be an explicit finite integer in [", min_val, ", 2^53 - 1]."),
      c("cdid_ipw_covariance_transfer_bound_error", "cdid_ipw_covariance_transfer_control_error")
    )
  }
  as.double(x)
}

ipw_cov_transfer_safe_product <- function(a, b, max_entries = NULL, label = "allocation") {
  if (!gmm_scalar_numeric(a, nonneg = TRUE) || a != floor(a) ||
      !gmm_scalar_numeric(b, nonneg = TRUE) || b != floor(b)) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Invalid integer dimensions for ", label, "."),
      "cdid_ipw_covariance_transfer_bound_error"
    )
  }
  da <- as.double(a)
  db <- as.double(b)
  limit <- 2^53 - 1
  if (da > 0 && db > limit / da) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Count product for ", label, " overflows exact IEEE binary64 integer range."),
      "cdid_ipw_covariance_transfer_bound_error",
      diagnostics = list(label = label, dim_a = da, dim_b = db)
    )
  }
  prod_val <- da * db
  if (!is.finite(prod_val) || prod_val > limit) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Count product for ", label, " exceeds exact integer range."),
      "cdid_ipw_covariance_transfer_bound_error"
    )
  }
  if (!is.null(max_entries) && prod_val > max_entries) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Count for ", label, " (", prod_val, ") exceeds explicit bound (", max_entries, ")."),
      "cdid_ipw_covariance_transfer_bound_error",
      diagnostics = list(label = label, requested = prod_val, bound = max_entries)
    )
  }
  prod_val
}

ipw_cov_transfer_safe_sum <- function(vals, max_bound = NULL, label = "sum") {
  limit <- 2^53 - 1
  total <- 0
  for (v in vals) {
    dv <- as.double(v)
    if (!is.finite(dv) || dv < 0 || dv != floor(dv) || dv > limit - total) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Cumulative count for ", label, " overflows exact IEEE binary64 integer range."),
        "cdid_ipw_covariance_transfer_bound_error",
        diagnostics = list(label = label, total_so_far = total, term = dv)
      )
    }
    total <- total + dv
  }
  if (!is.null(max_bound) && total > max_bound) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Total ", label, " (", total, ") exceeds explicit bound (", max_bound, ")."),
      "cdid_ipw_covariance_transfer_bound_error",
      diagnostics = list(label = label, requested = total, bound = max_bound)
    )
  }
  total
}

ipw_cov_transfer_refuse_decorative_attrs <- function(x, name, allowed_attrs = character()) {
  ax <- attributes(x)
  if (is.null(ax)) return(invisible(TRUE))
  forbidden_names <- c(
    "certified", "is_pd", "certified_pd", "bypass_pd_guard", "bypass_guard",
    "trust_metadata", "trust_metadata_pd", "same_frame", "same_complete_a_frame",
    "verified", "allow_pad_zeros", "allow_rescale", "pad_outside_a",
    "is_exact", "pre_validated"
  )
  bad <- intersect(names(ax), forbidden_names)
  if (length(bad) > 0L) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Decorative or shortcut metadata attributes (`", paste(bad, collapse = "`, `"), "`) on `", name, "` are refused."),
      c("cdid_ipw_covariance_transfer_metadata_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (length(allowed_attrs) > 0L) {
    extra <- setdiff(names(ax), allowed_attrs)
    if (length(extra) > 0L) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Unexpected attributes (`", paste(extra, collapse = "`, `"), "`) on `", name, "` are refused."),
        c("cdid_ipw_covariance_transfer_metadata_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  invisible(TRUE)
}

# Outward distance sup_{y in iv} |y - x| for finite point x and interval iv = [lo, hi].
# When an endpoint equals x in binary64, their real difference is identically 0;
# otherwise outward binary64 subtraction encloses |endpoint - x|.
ipw_cov_transfer_iv_point_dist <- function(iv, x) {
  x_dbl <- as.double(x)
  lo <- as.double(iv[1L])
  hi <- as.double(iv[2L])
  if (!is.finite(x_dbl) || !is.finite(lo) || !is.finite(hi) || lo > hi) {
    ipw_cov_transfer_arith_refuse(
      "Nonfinite or inverted interval/point in outward distance evaluation.",
      reason = "invalid_interval_point_distance",
      diagnostics = list(iv = iv, x = x)
    )
  }
  d_lo <- if (lo == x_dbl) 0.0 else gmm_iv_abs(logit_b64_sub(gmm_iv(lo), gmm_iv(x_dbl)))[2L]
  d_hi <- if (hi == x_dbl) 0.0 else gmm_iv_abs(logit_b64_sub(gmm_iv(hi), gmm_iv(x_dbl)))[2L]
  max(d_lo, d_hi)
}

# Chunked character-byte temporary cost accounting before string copies/comparisons.
ipw_cov_transfer_account_char_bytes <- function(
  vec,
  n_len,
  step,
  label,
  max_matrix_entries,
  base_concurrent_entries,
  char_bytes_so_far = 0
) {
  total_bytes <- as.double(char_bytes_so_far)
  a_start <- 1L
  while (a_start <= n_len) {
    ipw_cov_transfer_check_interrupt()
    a_end <- as.integer(min(as.double(n_len), as.double(a_start) + step - 1))
    ch <- vec[seq.int(a_start, a_end)]
    if (anyNA(ch)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`", label, "` contains NA character identifiers."),
        c("cdid_ipw_covariance_transfer_id_error", "cdid_ipw_covariance_transfer_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    nb <- nchar(ch, type = "bytes", allowNA = FALSE)
    if (any(nb < 1L)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`", label, "` contains empty character identifiers."),
        c("cdid_ipw_covariance_transfer_id_error", "cdid_ipw_covariance_transfer_input_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    chunk_bytes <- sum(as.double(nb))
    # Account for both character byte payload and comparison/hash temporary overhead (2x bytes)
    total_bytes <- ipw_cov_transfer_safe_sum(
      c(total_bytes, 2 * chunk_bytes),
      label = paste0("character temporary bytes for ", label)
    )
    char_entries_8byte <- ceiling(total_bytes / 8)
    ipw_cov_transfer_safe_sum(
      c(base_concurrent_entries, char_entries_8byte),
      max_bound = max_matrix_entries,
      label = paste0("concurrent entries plus character temporary cost for ", label)
    )
    a_start <- as.double(a_end) + 1
  }
  list(
    total_char_temp_bytes = total_bytes,
    total_char_temp_entries_8byte = ceiling(total_bytes / 8)
  )
}

# ---------------------------------------------------------------------------
# Explicit caller control, frame, operator, and candidate-row constructors
# (NO internal defaults)
# ---------------------------------------------------------------------------

ipw_covariance_transfer_control <- function(
  max_chunk_rows,
  max_operator_block_rows,
  max_cluster_rows,
  max_matrix_entries,
  max_work_ops,
  max_output_entries,
  link_control
) {
  ipw_cov_transfer_require_deps()
  if (missing(max_chunk_rows) || missing(max_operator_block_rows) ||
      missing(max_cluster_rows) || missing(max_matrix_entries) ||
      missing(max_work_ops) || missing(max_output_entries) ||
      missing(link_control)) {
    cdid_ipw_covariance_transfer_abort(
      "All fields of `ipw_covariance_transfer_control` (`max_chunk_rows`, `max_operator_block_rows`, `max_cluster_rows`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`, `link_control`) must be supplied explicitly by the caller.",
      c("cdid_ipw_covariance_transfer_control_error", "cdid_ipw_covariance_transfer_bound_error")
    )
  }

  max_ch <- ipw_cov_transfer_int_bound(max_chunk_rows, "max_chunk_rows", min_val = 1)
  max_op <- ipw_cov_transfer_int_bound(max_operator_block_rows, "max_operator_block_rows", min_val = 1)
  max_cl <- ipw_cov_transfer_int_bound(max_cluster_rows, "max_cluster_rows", min_val = 0)
  max_me <- ipw_cov_transfer_int_bound(max_matrix_entries, "max_matrix_entries", min_val = 1)
  max_wo <- ipw_cov_transfer_int_bound(max_work_ops, "max_work_ops", min_val = 1)
  max_out <- ipw_cov_transfer_int_bound(max_output_entries, "max_output_entries", min_val = 1)
  lk_ctrl <- validate_ipw_link_enclosure_control(link_control)

  structure(
    list(
      max_chunk_rows = max_ch,
      max_operator_block_rows = max_op,
      max_cluster_rows = max_cl,
      max_matrix_entries = max_me,
      max_work_ops = max_wo,
      max_output_entries = max_out,
      link_control = lk_ctrl
    ),
    class = "cdid_ipw_covariance_transfer_control"
  )
}

validate_ipw_covariance_transfer_control <- function(control) {
  ipw_cov_transfer_require_deps()
  if (missing(control) || is.null(control) || !is.list(control)) {
    cdid_ipw_covariance_transfer_abort(
      "Explicit `control` list (`max_chunk_rows`, `max_operator_block_rows`, `max_cluster_rows`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`, `link_control`) is required; no internal defaults exist.",
      c("cdid_ipw_covariance_transfer_control_error", "cdid_ipw_covariance_transfer_bound_error")
    )
  }
  ipw_cov_transfer_refuse_decorative_attrs(control, "control", allowed_attrs = c("names", "class"))
  req <- c(
    "max_chunk_rows", "max_operator_block_rows", "max_cluster_rows",
    "max_matrix_entries", "max_work_ops", "max_output_entries", "link_control"
  )
  if (length(control) != length(req) || is.null(names(control)) ||
      anyDuplicated(names(control)) || !setequal(names(control), req)) {
    missing_f <- setdiff(req, names(control))
    extra_f <- setdiff(names(control), req)
    cdid_ipw_covariance_transfer_abort(
      paste0(
        "Invalid covariance-transfer `control` fields.",
        if (length(missing_f)) paste0(" Missing: ", paste(missing_f, collapse = ", "), ".") else "",
        if (length(extra_f)) paste0(" Unexpected: ", paste(extra_f, collapse = ", "), ".") else ""
      ),
      c("cdid_ipw_covariance_transfer_control_error", "cdid_ipw_covariance_transfer_bound_error")
    )
  }
  ipw_covariance_transfer_control(
    max_chunk_rows = control$max_chunk_rows,
    max_operator_block_rows = control$max_operator_block_rows,
    max_cluster_rows = control$max_cluster_rows,
    max_matrix_entries = control$max_matrix_entries,
    max_work_ops = control$max_work_ops,
    max_output_entries = control$max_output_entries,
    link_control = control$link_control
  )
}

ipw_covariance_frame_spec <- function(
  id,
  normalizer_n,
  embedding,
  w = NULL,
  G = NULL
) {
  ipw_cov_transfer_require_deps()
  if (missing(id) || missing(normalizer_n) || missing(embedding)) {
    cdid_ipw_covariance_transfer_abort(
      "All required fields of `ipw_covariance_frame_spec` (`id`, `normalizer_n`, `embedding`) must be supplied explicitly.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  ipw_cov_transfer_refuse_decorative_attrs(id, "frame$id")
  norm_n <- ipw_cov_transfer_int_bound(normalizer_n, "frame$normalizer_n", min_val = 1)
  if (!gmm_scalar_character(embedding)) {
    cdid_ipw_covariance_transfer_abort(
      "`frame$embedding` must be a non-empty character scalar.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (!identical(embedding, "identity")) {
    cdid_ipw_covariance_transfer_abort(
      paste0(
        "Unsupported covariance frame embedding `", embedding, "`. ",
        "This component requires an explicit `'identity'` embedding into the identical complete source A-frame; ",
        "larger-frame transport, silent zero-row padding, and A-row rescaling are refused."
      ),
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error"),
      diagnostics = list(embedding = embedding)
    )
  }
  if (is.null(id) || is.factor(id) || is.object(id) || !is.atomic(id) ||
      !is.null(dim(id)) || !(is.numeric(id) || is.character(id)) ||
      as.double(length(id)) != norm_n) {
    cdid_ipw_covariance_transfer_abort(
      paste0("`frame$id` must be a plain atomic numeric or character vector of length `normalizer_n` (", norm_n, ")."),
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_id_error")
    )
  }
  if (!is.null(w)) {
    ipw_cov_transfer_refuse_decorative_attrs(w, "frame$w")
    if (is.object(w) || !is.atomic(w) || !is.null(dim(w)) || !is.numeric(w) ||
        as.double(length(w)) != norm_n) {
      cdid_ipw_covariance_transfer_abort(
        "`frame$w` must be NULL or a plain numeric vector of length `normalizer_n`.",
        c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  if (!is.null(G)) {
    ipw_cov_transfer_refuse_decorative_attrs(G, "frame$G")
    if (is.object(G) || !is.atomic(G) || !is.null(dim(G)) || !is.numeric(G) ||
        as.double(length(G)) != norm_n) {
      cdid_ipw_covariance_transfer_abort(
        "`frame$G` must be NULL or a plain numeric vector of length `normalizer_n`.",
        c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  structure(
    list(
      id = id,
      normalizer_n = norm_n,
      embedding = "identity",
      w = w,
      G = G
    ),
    class = "cdid_ipw_covariance_frame_spec"
  )
}

ipw_covariance_operator_spec <- function(
  purpose,
  row_level,
  omega,
  centered,
  kappa,
  cluster_id = NULL,
  cluster_order = NULL
) {
  ipw_cov_transfer_require_deps()
  if (missing(purpose) || missing(row_level) || missing(omega) ||
      missing(centered) || missing(kappa)) {
    cdid_ipw_covariance_transfer_abort(
      "All required fields of `ipw_covariance_operator_spec` (`purpose`, `row_level`, `omega`, `centered`, `kappa`) must be supplied explicitly.",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (!gmm_scalar_character(purpose) || !purpose %in% c("GMM", "inference")) {
    cdid_ipw_covariance_transfer_abort(
      "`op$purpose` must be explicitly `'GMM'` or `'inference'`.",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (!gmm_scalar_character(row_level) || !row_level %in% c("unit", "cluster")) {
    cdid_ipw_covariance_transfer_abort(
      "`op$row_level` must be explicitly `'unit'` or `'cluster'`.",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (identical(purpose, "GMM") && !identical(row_level, "unit")) {
    cdid_ipw_covariance_transfer_abort(
      "Under D11, `op$purpose == 'GMM'` requires unit aggregation (`row_level = 'unit'`); cluster covariance operators are permitted only for `purpose = 'inference'`.",
      c("cdid_ipw_covariance_transfer_d11_error", "cdid_ipw_covariance_transfer_operator_error")
    )
  }
  if (!gmm_scalar_logical(centered)) {
    cdid_ipw_covariance_transfer_abort(
      "`op$centered` must be a non-missing scalar logical.",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (!gmm_scalar_numeric(kappa, positive = TRUE)) {
    cdid_ipw_covariance_transfer_abort(
      "`op$kappa` must be a finite positive numeric scalar.",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  ipw_cov_transfer_refuse_decorative_attrs(omega, "op$omega")
  if (is.null(omega) || is.object(omega) || !is.atomic(omega) ||
      !is.null(dim(omega)) || !is.numeric(omega) || length(omega) < 1L) {
    cdid_ipw_covariance_transfer_abort(
      "`op$omega` must be a non-empty plain numeric vector of outer operator weights.",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (identical(row_level, "cluster")) {
    if (is.null(cluster_id) || is.null(cluster_order)) {
      cdid_ipw_covariance_transfer_abort(
        "Cluster covariance operator (`row_level == 'cluster'`) requires explicit `cluster_id` and `cluster_order` vectors.",
        c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
    ipw_cov_transfer_refuse_decorative_attrs(cluster_id, "op$cluster_id")
    ipw_cov_transfer_refuse_decorative_attrs(cluster_order, "op$cluster_order")
    if (is.factor(cluster_id) || is.object(cluster_id) || !is.atomic(cluster_id) ||
        !is.null(dim(cluster_id)) || !(is.numeric(cluster_id) || is.character(cluster_id)) ||
        length(cluster_id) < 1L) {
      cdid_ipw_covariance_transfer_abort(
        "`op$cluster_id` must be a plain non-object atomic numeric or character vector.",
        c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
    if (is.factor(cluster_order) || is.object(cluster_order) || !is.atomic(cluster_order) ||
        !is.null(dim(cluster_order)) || !(is.numeric(cluster_order) || is.character(cluster_order)) ||
        length(cluster_order) < 1L || is.numeric(cluster_id) != is.numeric(cluster_order)) {
      cdid_ipw_covariance_transfer_abort(
        "`op$cluster_order` must be a plain non-object atomic vector of the same type as `op$cluster_id`.",
        c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
    if (length(omega) != length(cluster_order)) {
      cdid_ipw_covariance_transfer_abort(
        "`op$omega` length must equal `length(op$cluster_order)` for cluster aggregation.",
        c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  } else {
    if (!is.null(cluster_id) || !is.null(cluster_order)) {
      if (identical(purpose, "GMM")) {
        cdid_ipw_covariance_transfer_abort(
          "Under D11, `op$purpose == 'GMM'` uses unit identity `A = I_n` directly and forbids `cluster_id`/`cluster_order`.",
          c("cdid_ipw_covariance_transfer_d11_error", "cdid_ipw_covariance_transfer_operator_error")
        )
      }
      if (is.null(cluster_id) || is.null(cluster_order)) {
        cdid_ipw_covariance_transfer_abort(
          "If supplied on a unit-level operator, both `cluster_id` and `cluster_order` must be present and represent the identity order.",
          c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
      ipw_cov_transfer_refuse_decorative_attrs(cluster_id, "op$cluster_id")
      ipw_cov_transfer_refuse_decorative_attrs(cluster_order, "op$cluster_order")
      if (is.factor(cluster_id) || is.object(cluster_id) || !is.atomic(cluster_id) || !is.null(dim(cluster_id)) ||
          is.factor(cluster_order) || is.object(cluster_order) || !is.atomic(cluster_order) || !is.null(dim(cluster_order)) ||
          !(is.numeric(cluster_id) || is.character(cluster_id)) ||
          is.numeric(cluster_id) != is.numeric(cluster_order) ||
          length(cluster_id) != length(omega) || length(cluster_order) != length(omega)) {
        cdid_ipw_covariance_transfer_abort(
          "Unit-level `cluster_id` and `cluster_order` headers must match `omega` length and type.",
          c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
    }
  }

  structure(
    list(
      purpose = purpose,
      row_level = row_level,
      omega = omega,
      centered = as.logical(centered),
      kappa = as.double(kappa),
      cluster_id = cluster_id,
      cluster_order = cluster_order
    ),
    class = "cdid_ipw_covariance_operator_spec"
  )
}

ipw_candidate_rows_spec <- function(callback, link_ids = NULL) {
  ipw_cov_transfer_require_deps()
  if (missing(callback) || !is.function(callback)) {
    cdid_ipw_covariance_transfer_abort(
      "`candidate_rows` requires an explicit function `callback(a_start, a_end)`.",
      c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  ipw_cov_transfer_refuse_decorative_attrs(callback, "candidate_rows$callback")
  if (!is.null(link_ids)) {
    ipw_cov_transfer_refuse_decorative_attrs(link_ids, "candidate_rows$link_ids")
    if (is.object(link_ids) || !is.atomic(link_ids) || !is.null(dim(link_ids)) ||
        !is.character(link_ids) || length(link_ids) < 1L || anyNA(link_ids) ||
        any(!nzchar(link_ids)) || anyDuplicated(link_ids)) {
      cdid_ipw_covariance_transfer_abort(
        "`candidate_rows$link_ids` must be NULL or a vector of unique non-empty character link IDs.",
        c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  structure(
    list(
      callback = callback,
      link_ids = link_ids
    ),
    class = "cdid_ipw_candidate_rows_spec"
  )
}

# ---------------------------------------------------------------------------
# Preflight accounting (executed BEFORE content scans, callbacks, or copies)
# ---------------------------------------------------------------------------

preflight_ipw_covariance_transfer <- function(
  n,
  L,
  d,
  row_level,
  n_clusters,
  p_vec,
  m_vec,
  certificate_controls,
  control,
  id_is_character = FALSE,
  cluster_is_character = FALSE,
  existing_4bd_retained_entries = 0
) {
  ipw_cov_transfer_require_deps()
  ctrl <- validate_ipw_covariance_transfer_control(control)

  n_d <- ipw_cov_transfer_int_bound(n, "n_a_units", min_val = 1)
  L_d <- ipw_cov_transfer_int_bound(L, "n_links", min_val = 1)
  d_d <- ipw_cov_transfer_int_bound(d, "n_moments", min_val = 1)
  if (!gmm_scalar_character(row_level) || !row_level %in% c("unit", "cluster")) {
    cdid_ipw_covariance_transfer_abort(
      "`row_level` in preflight must be `'unit'` or `'cluster'`.",
      "cdid_ipw_covariance_transfer_operator_error"
    )
  }
  N_c <- ipw_cov_transfer_int_bound(n_clusters, "n_clusters", min_val = 1)
  if (identical(row_level, "unit") && N_c != n_d) {
    cdid_ipw_covariance_transfer_abort(
      "Unit operator requires `n_clusters == n_a_units`.",
      "cdid_ipw_covariance_transfer_operator_error"
    )
  }
  if (identical(row_level, "cluster") && N_c > ctrl$max_cluster_rows) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Declared cluster count (", N_c, ") exceeds `control$max_cluster_rows` (", ctrl$max_cluster_rows, ")."),
      "cdid_ipw_covariance_transfer_bound_error",
      diagnostics = list(requested_clusters = N_c, max_cluster_rows = ctrl$max_cluster_rows)
    )
  }

  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  b_in <- min(n_d, ctrl$max_chunk_rows)
  b_op <- if (identical(row_level, "unit")) {
    min(b_in, ctrl$max_operator_block_rows)
  } else {
    min(N_c, ctrl$max_operator_block_rows)
  }
  n_input_chunks <- ceiling(n_d / b_in)
  n_operator_blocks <- if (identical(row_level, "unit")) {
    n_input_chunks * ceiling(b_in / b_op)
  } else {
    ceiling(N_c / b_op)
  }

  # 1. Simultaneous 4bd helper buffers and whole-stream 4bd chunk work
  if (length(p_vec) != L_d || length(m_vec) != L_d || length(certificate_controls) != L_d) {
    cdid_ipw_covariance_transfer_abort(
      "Link dimension/control vectors in preflight must have length `L`.",
      "cdid_ipw_covariance_transfer_input_error"
    )
  }
  max_4bd_concurrent_entries <- 0
  max_4bd_output_entries <- 0
  total_4bd_stream_work_ops <- 0
  link_chunk_preflights <- vector("list", as.integer(L_d))
  for (l in seq_len(as.integer(L_d))) {
    pre_l <- preflight_ipw_link_if_chunk(
      b_len = b_in,
      p = p_vec[l],
      m = m_vec[l],
      control = ctrl$link_control,
      certificate_control = certificate_controls[[l]],
      include_components = FALSE,
      include_shared_xi = FALSE,
      has_phi_candidate = TRUE,
      has_xi_beta = FALSE
    )
    link_chunk_preflights[[l]] <- pre_l
    if (pre_l$concurrent_entries > max_4bd_concurrent_entries) {
      max_4bd_concurrent_entries <- pre_l$concurrent_entries
    }
    if (pre_l$output_entries > max_4bd_output_entries) {
      max_4bd_output_entries <- pre_l$output_entries
    }
    per_row_4bd_l <- ceiling(pre_l$work_ops_bound / b_in) + 2
    stream_4bd_l <- ipw_cov_transfer_safe_product(
      n_d, per_row_4bd_l, ctrl$link_control$max_work_ops,
      paste0("whole-stream 4bd IF chunk work for link ", l)
    )
    total_4bd_stream_work_ops <- ipw_cov_transfer_safe_sum(
      c(total_4bd_stream_work_ops, stream_4bd_l),
      max_bound = max_wo,
      label = "cumulative 4bd IF chunk work across all links"
    )
  }

  # 2. Component concurrent payload entries (8-byte double equivalents)
  bL <- ipw_cov_transfer_safe_product(b_in, L_d, max_me, "b x L candidate/error chunk block")
  b_by_L_entries <- ipw_cov_transfer_safe_sum(
    c(ipw_cov_transfer_safe_product(3, bL, max_me, "b x L Phi_block, E_eval_block, E_root_block"), b_in),
    max_bound = max_me,
    label = "bounded b x L IF and separated error buffers"
  )

  bd_in <- ipw_cov_transfer_safe_product(b_in, d_d, max_me, "b x d input transform block")
  bd_op <- ipw_cov_transfer_safe_product(b_op, d_d, max_me, "b_op x d operator Gram block")
  b_by_d_entries <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(14, bd_in, max_me, "b x d Uhat, [U], residuals, and mapped error blocks"),
      ipw_cov_transfer_safe_product(12, bd_op, max_me, "five operator input slices, Yhat, |Yhat|, E_mat, and native transpose temporaries")
    ),
    max_bound = max_me,
    label = "bounded b x d transform and operator buffers"
  )

  dL <- ipw_cov_transfer_safe_product(d_d, L_d, max_me, "d x L moment map C")
  d_by_L_entries <- ipw_cov_transfer_safe_product(4, dL, max_me, "caller C, copied C, t(C), and |t(C)| matrices")

  d2 <- ipw_cov_transfer_safe_product(d_d, d_d, max_me, "d x d covariance matrix block")
  d_by_d_entries <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(24, d2, max_me, "V_candidate, S, [L,U], dist, point Gram, D1..D3, and spectral workspace"),
      ipw_cov_transfer_safe_product(4, d_d, max_me, "global t accumulator and moment vector scratch"),
      256
    ),
    max_bound = max_me,
    label = "d x d covariance accumulators and spectral scratch"
  )

  id_membership_entries <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(6, n_d, max_me, "full-frame IDs, masks, and unit validation vectors"),
      if (identical(row_level, "cluster")) {
        ipw_cov_transfer_safe_sum(
          c(
            ipw_cov_transfer_safe_product(2, n_d, max_me, "cluster_id key and unit-to-cluster index map"),
            ipw_cov_transfer_safe_product(2, N_c, max_me, "cluster_order key and cluster omega vector")
          ),
          max_bound = max_me,
          label = "cluster membership and weight vectors"
        )
      } else {
        ipw_cov_transfer_safe_product(1, n_d, max_me, "unit omega vector")
      },
      ipw_cov_transfer_safe_product(66, L_d, max_me, "ordered link IDs, link headers, and retained per-link preflight scalar payload")
    ),
    max_bound = max_me,
    label = "full-frame IDs, membership, and weight vectors"
  )

  cluster_buffer_entries <- if (identical(row_level, "cluster")) {
    Ncd <- ipw_cov_transfer_safe_product(N_c, d_d, max_me, "N_c x d cluster accumulator matrix")
    ipw_cov_transfer_safe_product(
      5, Ncd, max_me,
      "explicit N_c x d cluster accumulators (z_hat, z_lo, z_hi, ae_eval, ae_root)"
    )
  } else {
    0
  }

  min_char_temp_bytes <- ipw_cov_transfer_safe_sum(
    c(
      if (isTRUE(id_is_character)) 2 * (n_d + b_in) else 0,
      if (isTRUE(cluster_is_character)) 2 * (n_d + N_c) else 0,
      2 * L_d
    ),
    label = "minimum character comparison temporary bytes"
  )
  min_char_temp_entries <- ceiling(min_char_temp_bytes / 8)

  base_component_concurrent_entries <- ipw_cov_transfer_safe_sum(
    c(
      b_by_L_entries,
      b_by_d_entries,
      d_by_L_entries,
      d_by_d_entries,
      id_membership_entries,
      cluster_buffer_entries,
      max_4bd_concurrent_entries
    ),
    max_bound = max_me,
    label = "simultaneous component and 4bd chunk buffers"
  )

  concurrent_entries <- ipw_cov_transfer_safe_sum(
    c(base_component_concurrent_entries, min_char_temp_entries),
    max_bound = max_me,
    label = "total concurrent payload entries (including 4bd chunk buffers and character temporaries)"
  )

  # 3. Output payload entry bound
  output_entries <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(6, d2, max_out, "enclosure [L,U], distance, V_candidate, and both retained S endpoints"),
      dL,
      ipw_cov_transfer_safe_product(4, d_d, max_out, "global t interval and moment summaries"),
      if (identical(row_level, "cluster")) {
        ipw_cov_transfer_safe_product(3, Ncd, max_out, "retained cluster point and both interval endpoint arrays")
      } else 0,
      ipw_cov_transfer_safe_product(64, L_d, max_out, "retained per-link preflight scalar payload"),
      256
    ),
    max_bound = max_out,
    label = "covariance transfer output entries"
  )

  # 4. Whole-stream component interval operation bound
  per_input_row_transform_ops <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(6, dL, max_wo, "per-row Uhat residual and |E_eval| |C|', |E_root| |C|' dot ops"),
      ipw_cov_transfer_safe_product(32, d_d, max_wo, "per-row transform residual subtraction, native returned norms, interval addition, and error extraction ops"),
      16,
      if (identical(row_level, "cluster")) ipw_cov_transfer_safe_product(3, d_d, max_wo, "per-row outward cluster accumulation ops") else 0
    ),
    max_bound = max_wo,
    label = "per-input-row transform and cluster accumulation ops"
  )
  input_transform_ops <- ipw_cov_transfer_safe_product(
    n_d, per_input_row_transform_ops, max_wo, "whole-stream input transform and cluster sum ops"
  )

  per_operator_row_ops <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(8, d2, max_wo, "per-operator-row point Gram residual and D1, D2, D3 dot ops"),
      ipw_cov_transfer_safe_product(16, d_d, max_wo, "per-operator-row sqrt(omega), y_iv, E_mat, t, and nu_kind^2 ops"),
      8
    ),
    max_bound = max_wo,
    label = "per-operator-row Gram and error accumulation ops"
  )
  operator_row_ops <- ipw_cov_transfer_safe_product(
    N_c, per_operator_row_ops, max_wo, "whole-stream operator-row Gram ops"
  )
  operator_block_merge_ops <- ipw_cov_transfer_safe_product(
    n_operator_blocks,
    ipw_cov_transfer_safe_sum(
      c(ipw_cov_transfer_safe_product(32, d2, max_wo, "per-operator-block residual subtraction, native returned norms, D sum, widening, and S addition ops"), 32),
      max_bound = max_wo, label = "operator-block native norm constants"
    ),
    max_wo,
    "whole-stream operator-block merge ops"
  )

  d3 <- ipw_cov_transfer_safe_product(d2, d_d, max_wo, "d^3 spectral residual ops")
  post_gram_ops <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(6, d3, max_wo, "spectral certificate Q'Q and S*Q-Q*Lambda residual ops"),
      ipw_cov_transfer_safe_product(64, d2, max_wo, "global centering, kappa scaling, distance, returned norms, and eigen transfer ops"),
      ipw_cov_transfer_safe_product(64, d_d, max_wo, "post-Gram vector ops"),
      256
    ),
    max_bound = max_wo,
    label = "post-Gram centering, norm, and eigen transfer ops"
  )

  component_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(input_transform_ops, operator_row_ops, operator_block_merge_ops, post_gram_ops),
    max_bound = max_wo,
    label = "component covariance transfer interval work bound"
  )

  total_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(total_4bd_stream_work_ops, component_work_ops_bound),
    max_bound = max_wo,
    label = "whole-stream 4bd plus covariance transfer interval work bound"
  )

  list(
    n_a_units = n_d,
    n_links = L_d,
    n_moments = d_d,
    row_level = row_level,
    n_clusters = N_c,
    effective_chunk_rows = logit_cert_count_output(b_in),
    effective_operator_block_rows = logit_cert_count_output(b_op),
    n_input_chunks = logit_cert_count_output(n_input_chunks),
    n_operator_blocks = logit_cert_count_output(n_operator_blocks),
    b_by_L_entries = b_by_L_entries,
    b_by_d_entries = b_by_d_entries,
    d_by_L_entries = d_by_L_entries,
    d_by_d_entries = d_by_d_entries,
    id_and_membership_entries = id_membership_entries,
    cluster_buffer_entries = cluster_buffer_entries,
    max_4bd_chunk_concurrent_entries = max_4bd_concurrent_entries,
    max_4bd_chunk_output_entries = max_4bd_output_entries,
    base_component_concurrent_entries = base_component_concurrent_entries,
    character_temporary_bytes = min_char_temp_bytes,
    character_temporary_entries_8byte = min_char_temp_entries,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    output_entries = output_entries,
    total_4bd_stream_work_ops_bound = total_4bd_stream_work_ops,
    component_work_ops_bound = component_work_ops_bound,
    total_work_ops_bound = total_work_ops_bound,
    existing_4bd_retained_entries = as.double(existing_4bd_retained_entries),
    total_with_existing_4bd_entries = ipw_cov_transfer_safe_sum(
      c(existing_4bd_retained_entries, concurrent_entries),
      label = "pre-existing 4bd retained plus component concurrent payload"
    ),
    total_with_existing_4bd_bytes_8double = 8 * ipw_cov_transfer_safe_sum(
      c(existing_4bd_retained_entries, concurrent_entries),
      label = "pre-existing 4bd retained plus component concurrent payload bytes"
    ),
    link_chunk_preflights = link_chunk_preflights,
    n_by_L_cache_allocated = FALSE,
    dense_A_incidence_allocated = FALSE,
    per_block_full_n_validation = FALSE,
    peak_rss_promised = FALSE,
    storage_and_work_disclosure = paste0(
      "Preflight bounds simultaneous 4bd chunk helper buffers, bounded b x L IF/error arrays, ",
      "b x d and b_op x d transforms, d x L map C, d x d accumulators/scratch, full-frame IDs/membership, ",
      "character comparison temporaries, and explicit N_c x d cluster buffers (when clustered), ",
      "plus whole-stream 4bd and covariance interval operations. Pre-existing 4bd retained problem/link ",
      "storage and opaque caller callback internal work/storage are disclosed separately; ",
      "payload entry counts are 8-byte double equivalents, not process peak RSS."
    )
  )
}

# ---------------------------------------------------------------------------
# Internal validators for links, common A-frame, operator, and matrices
# ---------------------------------------------------------------------------

ipw_cov_transfer_validate_links_headers <- function(links) {
  if (missing(links) || is.null(links) || !is.list(links) || is.object(links) || length(links) < 1L) {
    cdid_ipw_covariance_transfer_abort(
      "`links` must be a non-empty plain list of sealed `cdid_ipw_link_enclosure` objects.",
      c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  ipw_cov_transfer_refuse_decorative_attrs(links, "links", allowed_attrs = "names")
  L <- length(links)
  problems <- vector("list", L)
  p_vec <- numeric(L)
  m_vec <- numeric(L)
  cert_ctrls <- vector("list", L)
  link_ids <- character(L)
  list_names <- names(links)
  if (!is.null(list_names) && (anyNA(list_names) || any(!nzchar(list_names)) || anyDuplicated(list_names))) {
    cdid_ipw_covariance_transfer_abort(
      "If `links` is a named list, its names must be unique non-empty character link IDs.",
      c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }

  for (l in seq_len(L)) {
    lk <- links[[l]]
    prob_l <- validate_ipw_link_enclosure(lk, require_certified = FALSE)
    problems[[l]] <- prob_l
    p_vec[l] <- as.double(ncol(prob_l$raw$Z))
    m_vec[l] <- as.double(length(prob_l$raw$cohorts))
    cert_ctrls[[l]] <- prob_l$controls$certificate_control

    lid <- lk$link_id
    if (!is.null(lid)) {
      if (!gmm_scalar_character(lid)) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`links[[", l, "]]$link_id` must be a non-empty character scalar."),
          c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
      if (!is.null(list_names) && !identical(list_names[l], lid)) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`names(links)[", l, "]` ('", list_names[l], "') disagrees with `links[[", l, "]]$link_id` ('", lid, "')."),
          c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
      link_ids[l] <- lid
    } else if (!is.null(list_names)) {
      link_ids[l] <- list_names[l]
    } else {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` must declare an explicit `link_id` (or `names(links)` must supply explicit link IDs)."),
        c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }

  if (anyDuplicated(link_ids)) {
    cdid_ipw_covariance_transfer_abort(
      "Declared `link_id` values across `links` must be unique.",
      c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error"),
      diagnostics = list(link_ids = link_ids)
    )
  }

  prob_1 <- problems[[1L]]
  n_1 <- nrow(prob_1$raw$Z)
  norm_1 <- prob_1$provenance$normalizer_n
  id_type_1 <- prob_1$raw$id_type
  if (!identical(as.double(n_1), norm_1)) {
    cdid_ipw_covariance_transfer_abort(
      "Link 1 nuisance problem normalizer disagrees with `nrow(raw$Z)`.",
      c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_integrity_error")
    )
  }

  # Compute existing 4bd retained storage across unique nuisance problems and links
  seen_prob_refs <- list(prob_1$ref)
  existing_retained <- as.double(prob_1$preflight$total_retained_entries) +
    as.double(links[[1L]]$preflight$total_retained_entries)

  for (l in seq_len(L)[-1L]) {
    prob_l <- problems[[l]]
    if (nrow(prob_l$raw$Z) != n_1 || !identical(prob_l$provenance$normalizer_n, norm_1) ||
        !identical(prob_l$raw$id_type, id_type_1) ||
        length(prob_l$raw$id) != n_1 || length(prob_l$raw$w) != n_1 || length(prob_l$raw$G) != n_1) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` does not share the same A-unit count, normalizer, or ID type as `links[[1]]`."),
        c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_link_error")
      )
    }
    already_seen <- any(vapply(seen_prob_refs, identical, logical(1L), y = prob_l$ref))
    if (!already_seen) {
      seen_prob_refs[[length(seen_prob_refs) + 1L]] <- prob_l$ref
      existing_retained <- existing_retained + as.double(prob_l$preflight$total_retained_entries)
    }
    existing_retained <- existing_retained + as.double(links[[l]]$preflight$total_retained_entries)
  }

  list(
    L = L,
    n = n_1,
    normalizer_n = norm_1,
    id_type = id_type_1,
    link_ids = link_ids,
    problems = problems,
    n_unique_problems = length(seen_prob_refs),
    p_vec = p_vec,
    m_vec = m_vec,
    certificate_controls = cert_ctrls,
    existing_4bd_retained_entries = existing_retained
  )
}

# Verify masks, endpoint mappings, and common A-unit arrays across all links ONCE before streaming.
ipw_cov_transfer_validate_links_content <- function(links, link_hdr) {
  L <- link_hdr$L
  n <- link_hdr$n
  probs <- link_hdr$problems
  raw_1 <- probs[[1L]]$raw
  ref_1 <- probs[[1L]]$ref

  seen_prob_refs <- list()
  for (l in seq_len(L)) {
    ipw_cov_transfer_check_interrupt()
    lk <- links[[l]]
    prob_l <- probs[[l]]
    raw_l <- prob_l$raw
    masks_l <- prob_l$masks
    ep_l <- lk$endpoints

    # Verify problem masks once per unique nuisance problem
    already_verified_prob <- any(vapply(seen_prob_refs, identical, logical(1L), y = prob_l$ref))
    if (!already_verified_prob) {
      cohort_idx_calc <- match(raw_l$G, raw_l$cohorts)
      in_risk_calc <- !is.na(cohort_idx_calc)
      in_obs_calc <- in_risk_calc & raw_l$pair_observed
      if (!identical(masks_l$cohort_index, as.integer(cohort_idx_calc)) ||
          !identical(masks_l$in_risk, as.logical(in_risk_calc)) ||
          !identical(masks_l$in_observed, as.logical(in_obs_calc)) ||
          !identical(masks_l$risk_indices, as.integer(which(in_risk_calc))) ||
          !identical(masks_l$observed_indices, as.integer(which(in_obs_calc)))) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`links[[", l, "]]` nuisance problem masks disagree with retained raw `G`, `cohorts`, or `pair_observed`."),
          c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_integrity_error")
        )
      }
      if (!identical(prob_l$ref, ref_1)) {
        if (!identical(raw_l$id_type, raw_1$id_type) ||
            !identical(raw_l$id_key, raw_1$id_key) ||
            !identical(raw_l$w, raw_1$w) ||
            !identical(raw_l$G, raw_1$G) ||
            !identical(prob_l$provenance$normalizer_n, probs[[1L]]$provenance$normalizer_n)) {
          cdid_ipw_covariance_transfer_abort(
            paste0("`links[[", l, "]]` nuisance problem differs in A-unit `id`, `id_type`, `w`, `G`, or `normalizer_n` from `links[[1]]`."),
            c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_link_error")
          )
        }
      }
      seen_prob_refs[[length(seen_prob_refs) + 1L]] <- prob_l$ref
    }

    # Verify link provenance normalizer if present
    if (!is.null(lk$provenance) && !is.null(lk$provenance$normalizer_n) &&
        !identical(lk$provenance$normalizer_n, link_hdr$normalizer_n)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]$provenance$normalizer_n` disagrees with nuisance problem normalizer."),
        c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_link_error")
      )
    }

    # Verify link endpoint mapping against problem masks (required for certified links or when endpoints are present)
    if (isTRUE(lk$certified) || !is.null(ep_l)) {
      obs_idx <- masks_l$observed_indices
      n_obs <- length(obs_idx)
      if (!is.list(ep_l) ||
          !identical(ep_l$obs_a_index, obs_idx) ||
          !is.integer(ep_l$a_to_obs_pos) || length(ep_l$a_to_obs_pos) != n ||
          !identical(ep_l$a_to_obs_pos[obs_idx], seq_len(n_obs)) ||
          any(ep_l$a_to_obs_pos[!masks_l$in_observed] != 0L) ||
          !is.numeric(ep_l$y_left_obs) || length(ep_l$y_left_obs) != n_obs ||
          !is.numeric(ep_l$y_right_obs) || length(ep_l$y_right_obs) != n_obs ||
          anyNA(ep_l$y_left_obs) || any(!is.finite(ep_l$y_left_obs)) ||
          anyNA(ep_l$y_right_obs) || any(!is.finite(ep_l$y_right_obs)) ||
          !identical(ep_l$off_observed_endpoints_read, FALSE)) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`links[[", l, "]]$endpoints` mapping or retained observer array is malformed or inconsistent with `masks$observed_indices`."),
          c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_integrity_error")
        )
      }
    }
  }
  invisible(TRUE)
}

ipw_cov_transfer_validate_frame_header <- function(frame, n, id_type_ref) {
  if (missing(frame) || is.null(frame) || !is.list(frame)) {
    cdid_ipw_covariance_transfer_abort(
      "Explicit `frame` specification is required.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  ipw_cov_transfer_refuse_decorative_attrs(frame, "frame", allowed_attrs = c("names", "class"))
  if (is.null(names(frame)) || any(!nzchar(names(frame))) || anyDuplicated(names(frame))) {
    cdid_ipw_covariance_transfer_abort(
      "`frame` must be a uniquely named list.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  forbidden_fields <- c(
    "same_frame", "same_complete_a_frame", "certified", "verified",
    "allow_pad_zeros", "allow_rescale", "pad_outside_a", "trust_metadata",
    "rescale_factor", "pad_zeros"
  )
  bad_f <- intersect(names(frame), forbidden_fields)
  if (length(bad_f) > 0L) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Shortcut or unsupported larger-frame fields (`", paste(bad_f, collapse = "`, `"), "`) in `frame` are refused."),
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_metadata_error")
    )
  }
  allowed_fields <- c("id", "normalizer_n", "embedding", "w", "G", "n_a_units", "n_frame_units")
  extra_f <- setdiff(names(frame), allowed_fields)
  if (length(extra_f) > 0L) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Unexpected fields (`", paste(extra_f, collapse = "`, `"), "`) in `frame`."),
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (!all(c("id", "normalizer_n", "embedding") %in% names(frame))) {
    cdid_ipw_covariance_transfer_abort(
      "`frame` must explicitly supply `id`, `normalizer_n`, and `embedding`.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (!gmm_scalar_character(frame$embedding) || !identical(frame$embedding, "identity")) {
    cdid_ipw_covariance_transfer_abort(
      "Only `frame$embedding == 'identity'` into the identical complete source A-frame is supported; larger-frame transport or rescaling is refused.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  norm_f <- ipw_cov_transfer_int_bound(frame$normalizer_n, "frame$normalizer_n", min_val = 1)
  if (norm_f != as.double(n)) {
    cdid_ipw_covariance_transfer_abort(
      paste0("`frame$normalizer_n` (", norm_f, ") disagrees with source A-frame normalizer (", n, "); silent rescaling or zero padding is forbidden."),
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error"),
      diagnostics = list(frame_normalizer_n = norm_f, source_normalizer_n = n)
    )
  }
  for (opt_n in c("n_a_units", "n_frame_units")) {
    if (!is.null(frame[[opt_n]])) {
      val_n <- ipw_cov_transfer_int_bound(frame[[opt_n]], paste0("frame$", opt_n), min_val = 1)
      if (val_n != as.double(n)) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`frame$", opt_n, "` (", val_n, ") disagrees with source A-frame count (", n, ")."),
          c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
    }
  }
  f_spec <- ipw_covariance_frame_spec(
    id = frame$id,
    normalizer_n = norm_f,
    embedding = frame$embedding,
    w = frame$w,
    G = frame$G
  )
  f_id_type <- if (is.numeric(f_spec$id)) "numeric" else "character"
  if (!identical(f_id_type, id_type_ref)) {
    cdid_ipw_covariance_transfer_abort(
      paste0("`frame$id` type ('", f_id_type, "') disagrees with source `id_type` ('", id_type_ref, "')."),
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_id_error")
    )
  }
  f_spec
}

ipw_cov_transfer_validate_frame_content <- function(
  frame_spec,
  raw_ref,
  n,
  step,
  max_matrix_entries,
  base_concurrent_entries,
  char_bytes_so_far = 0
) {
  fid <- frame_spec$id
  is_num <- is.numeric(fid)
  char_info <- list(
    total_char_temp_bytes = as.double(char_bytes_so_far),
    total_char_temp_entries_8byte = ceiling(as.double(char_bytes_so_far) / 8)
  )

  # Even if `fid` aliases `raw_ref$id`, validate integrity and character byte budget
  if (!is_num) {
    char_info <- ipw_cov_transfer_account_char_bytes(
      vec = fid,
      n_len = n,
      step = step,
      label = "frame$id",
      max_matrix_entries = max_matrix_entries,
      base_concurrent_entries = base_concurrent_entries,
      char_bytes_so_far = char_bytes_so_far
    )
  } else {
    a_start <- 1L
    while (a_start <= n) {
      ipw_cov_transfer_check_interrupt()
      a_end <- as.integer(min(as.double(n), as.double(a_start) + step - 1))
      ch <- fid[seq.int(a_start, a_end)]
      if (anyNA(ch) || any(!is.finite(ch))) {
        cdid_ipw_covariance_transfer_abort(
          "`frame$id` contains NA, NaN, or infinite numeric identifiers.",
          c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_id_error")
        )
      }
      a_start <- as.double(a_end) + 1
    }
  }

  fid_key <- if (is_num) {
    if (gmm_sexp_aliases(fid, raw_ref$id) && typeof(fid) == "double") {
      raw_ref$id_key
    } else {
      d <- as.double(fid)
      d[d == 0.0] <- 0.0
      d
    }
  } else {
    if (gmm_sexp_aliases(fid, raw_ref$id)) raw_ref$id_key else as.character(fid)
  }

  if (anyDuplicated(fid_key)) {
    cdid_ipw_covariance_transfer_abort(
      "Duplicate unit identifiers in `frame$id` are forbidden.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_id_error")
    )
  }
  if (!identical(fid_key, raw_ref$id_key)) {
    cdid_ipw_covariance_transfer_abort(
      "`frame$id` differs in value or unit order from the complete source A-frame.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_id_error")
    )
  }

  if (!is.null(frame_spec$w)) {
    fw <- frame_spec$w
    if (anyNA(fw) || any(!is.finite(fw)) || !identical(as.double(fw), raw_ref$w)) {
      cdid_ipw_covariance_transfer_abort(
        "`frame$w` disagrees with the source A-frame case weights `w`.",
        c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  if (!is.null(frame_spec$G)) {
    fG <- frame_spec$G
    if (anyNA(fG) || !identical(as.double(fG), raw_ref$G)) {
      cdid_ipw_covariance_transfer_abort(
        "`frame$G` disagrees with the source A-frame cohort vector `G`.",
        c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }

  list(
    frame_id_key = fid_key,
    id_aliased_to_source = gmm_sexp_aliases(fid, raw_ref$id),
    char_info = char_info
  )
}

ipw_cov_transfer_validate_operator_header <- function(op, n) {
  if (missing(op) || is.null(op) || !is.list(op)) {
    cdid_ipw_covariance_transfer_abort(
      "Explicit covariance operator `op` is required.",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  ipw_cov_transfer_refuse_decorative_attrs(op, "op", allowed_attrs = c("names", "class"))
  if (is.null(names(op)) || any(!nzchar(names(op))) || anyDuplicated(names(op))) {
    cdid_ipw_covariance_transfer_abort(
      "`op` must be a uniquely named list.",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  forbidden_op <- c(
    "apply_case_weights", "multiply_case_weights", "drop_zero_clusters",
    "center_per_batch", "certified", "is_pd", "certified_pd", "trust_metadata"
  )
  bad_op <- intersect(names(op), forbidden_op)
  if (length(bad_op) > 0L) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Forbidden or shortcut fields (`", paste(bad_op, collapse = "`, `"), "`) in `op` are refused."),
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_metadata_error")
    )
  }
  allowed_op <- c("purpose", "row_level", "omega", "centered", "kappa", "cluster_id", "cluster_order")
  extra_op <- setdiff(names(op), allowed_op)
  if (length(extra_op) > 0L) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Unexpected fields (`", paste(extra_op, collapse = "`, `"), "`) in `op`."),
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  op_spec <- ipw_covariance_operator_spec(
    purpose = op$purpose,
    row_level = op$row_level,
    omega = op$omega,
    centered = op$centered,
    kappa = op$kappa,
    cluster_id = op$cluster_id,
    cluster_order = op$cluster_order
  )
  if (identical(op_spec$row_level, "unit")) {
    if (length(op_spec$omega) != n) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Unit-level `op$omega` must have length `n_a_units` (", n, ")."),
        c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  } else {
    if (length(op_spec$cluster_id) != n) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Cluster-level `op$cluster_id` must have length `n_a_units` (", n, ")."),
        c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  op_spec
}

ipw_cov_transfer_validate_operator_content <- function(
  op_spec,
  n,
  step,
  max_matrix_entries,
  base_concurrent_entries,
  char_bytes_so_far = 0
) {
  char_info <- list(
    total_char_temp_bytes = as.double(char_bytes_so_far),
    total_char_temp_entries_8byte = ceiling(as.double(char_bytes_so_far) / 8)
  )

  om <- op_spec$omega
  if (anyNA(om) || any(!is.finite(om)) || any(om < 0.0)) {
    cdid_ipw_covariance_transfer_abort(
      "`op$omega` must contain finite nonnegative outer weights (`omega >= 0`).",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  om_dbl <- as.double(om)
  om_dbl[om_dbl == 0.0] <- 0.0
  if (!any(om_dbl > 0.0)) {
    cdid_ipw_covariance_transfer_abort(
      "`op$omega` must have at least one strictly positive weight (`sum(omega) > 0`).",
      c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (identical(op_spec$purpose, "GMM") && any(om_dbl != 1.0)) {
    cdid_ipw_covariance_transfer_abort(
      "Under D11, `op$purpose == 'GMM'` fixes unit point-GMM weights (`omega == 1` for all units).",
      c("cdid_ipw_covariance_transfer_d11_error", "cdid_ipw_covariance_transfer_operator_error")
    )
  }

  unit_to_cluster <- NULL
  cluster_order_key <- NULL
  cluster_id_key <- NULL
  n_zero_unit_clusters <- 0L

  if (identical(op_spec$row_level, "cluster") || !is.null(op_spec$cluster_id)) {
    cid <- op_spec$cluster_id
    cord <- op_spec$cluster_order
    N_c <- length(cord)
    is_num_c <- is.numeric(cid)

    if (!is_num_c) {
      char_info <- ipw_cov_transfer_account_char_bytes(
        vec = cid,
        n_len = n,
        step = step,
        label = "op$cluster_id",
        max_matrix_entries = max_matrix_entries,
        base_concurrent_entries = base_concurrent_entries,
        char_bytes_so_far = char_info$total_char_temp_bytes
      )
      char_info <- ipw_cov_transfer_account_char_bytes(
        vec = cord,
        n_len = N_c,
        step = step,
        label = "op$cluster_order",
        max_matrix_entries = max_matrix_entries,
        base_concurrent_entries = base_concurrent_entries,
        char_bytes_so_far = char_info$total_char_temp_bytes
      )
      cluster_id_key <- as.character(cid)
      cluster_order_key <- as.character(cord)
    } else {
      if (anyNA(cid) || any(!is.finite(cid)) || anyNA(cord) || any(!is.finite(cord))) {
        cdid_ipw_covariance_transfer_abort(
          "`op$cluster_id` and `op$cluster_order` must be finite numeric identifiers without NA/NaN/Inf.",
          c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
      cluster_id_key <- as.double(cid)
      cluster_id_key[cluster_id_key == 0.0] <- 0.0
      cluster_order_key <- as.double(cord)
      cluster_order_key[cluster_order_key == 0.0] <- 0.0
    }

    if (anyDuplicated(cluster_order_key)) {
      cdid_ipw_covariance_transfer_abort(
        "`op$cluster_order` must not contain duplicate cluster identifiers.",
        c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }

    if (identical(op_spec$row_level, "unit")) {
      if (!identical(cluster_id_key, cluster_order_key)) {
        cdid_ipw_covariance_transfer_abort(
          "For `row_level == 'unit'`, any supplied `cluster_id` and `cluster_order` must represent the exact identity map `A = I_n`.",
          c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
    } else {
      unit_to_cluster <- match(cluster_id_key, cluster_order_key)
      if (anyNA(unit_to_cluster)) {
        cdid_ipw_covariance_transfer_abort(
          "Every unit's `op$cluster_id` must belong to declared `op$cluster_order`; dropping units or undeclared clusters is forbidden.",
          c("cdid_ipw_covariance_transfer_operator_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
      unit_counts <- tabulate(unit_to_cluster, nbins = N_c)
      n_zero_unit_clusters <- sum(unit_counts == 0L)
    }
  }

  list(
    purpose = op_spec$purpose,
    row_level = op_spec$row_level,
    omega = om_dbl,
    centered = op_spec$centered,
    kappa = op_spec$kappa,
    n_clusters = if (identical(op_spec$row_level, "cluster")) length(op_spec$cluster_order) else n,
    n_positive_weight_rows = sum(om_dbl > 0.0),
    n_zero_unit_declared_clusters = n_zero_unit_clusters,
    unit_to_cluster = unit_to_cluster,
    cluster_id_key = cluster_id_key,
    cluster_order_key = cluster_order_key,
    char_info = char_info
  )
}

# ---------------------------------------------------------------------------
# Reviewed Streaming Gram & Covariance-Transfer Kernel
# ---------------------------------------------------------------------------

ipw_cov_transfer_stream_kernel <- function(
  n,
  L,
  d,
  b_in,
  b_op,
  fetch_chunk,
  C_dbl,
  op_norm,
  V_cand_dbl,
  max_matrix_entries,
  preflight = NULL
) {
  tC <- t(C_dbl)
  abs_tC <- abs(tC)
  is_cluster <- identical(op_norm$row_level, "cluster")
  N_c <- op_norm$n_clusters
  omega <- op_norm$omega
  centered <- op_norm$centered
  kappa <- op_norm$kappa
  kappa_iv <- gmm_iv(kappa)
  zero_d_d <- matrix(0.0, nrow = d, ncol = d)

  actual_component_ops <- 0
  actual_4bd_ops <- 0
  n_input_chunks_seen <- 0L
  n_operator_blocks_seen <- 0L
  emitted_units <- 0L
  weight_rounding_occurred <- FALSE
  max_eval_err_seen <- 0.0
  max_root_err_seen <- 0.0
  max_direct_star_dist_seen <- 0.0

  # Explicit N_c x d cluster accumulators ONLY when op_norm$row_level == "cluster"
  cluster_z_hat <- NULL
  cluster_z_lo <- NULL
  cluster_z_hi <- NULL
  cluster_ae_eval <- NULL
  cluster_ae_root <- NULL
  if (is_cluster) {
    cluster_z_hat <- matrix(0.0, nrow = N_c, ncol = d)
    cluster_z_lo <- matrix(0.0, nrow = N_c, ncol = d)
    cluster_z_hi <- matrix(0.0, nrow = N_c, ncol = d)
    cluster_ae_eval <- matrix(0.0, nrow = N_c, ncol = d)
    cluster_ae_root <- matrix(0.0, nrow = N_c, ncol = d)
  }

  S_iv <- array(0.0, dim = c(d, d, 2L))
  t_iv <- matrix(0.0, nrow = d, ncol = 2L)
  W_iv <- c(0.0, 0.0)
  sum_eval_sq_iv <- c(0.0, 0.0)
  sum_root_sq_iv <- c(0.0, 0.0)

  # Helper to process one operator block of rows (either directly from a unit chunk or from completed clusters)
  process_operator_block <- function(z_hat_blk, z_lo_blk, z_hi_blk, ae_eval_blk, ae_root_blk, omega_blk) {
    m_blk <- nrow(z_hat_blk)
    Yhat <- matrix(0.0, nrow = m_blk, ncol = d)
    E_mat <- matrix(0.0, nrow = m_blk, ncol = d)

    for (r in seq_len(m_blk)) {
      w_r <- omega_blk[r]
      w_iv <- gmm_iv(w_r)
      sw_iv <- gmm_iv_sqrt(w_iv)
      sw_pt <- sqrt(w_r)
      if (!is.finite(sw_pt)) {
        ipw_cov_transfer_arith_refuse("Nonfinite sqrt(omega) in operator block.", "sqrt_omega_nonfinite")
      }
      if (sw_iv[1L] != sw_iv[2L]) {
        weight_rounding_occurred <<- TRUE
      }
      W_iv <<- gmm_iv_add(W_iv, w_iv)
      actual_component_ops <<- actual_component_ops + 2

      for (j in seq_len(d)) {
        z_rj_iv <- c(z_lo_blk[r, j], z_hi_blk[r, j])
        y_rj_iv <- gmm_iv_mul(sw_iv, z_rj_iv)
        yh_rj <- sw_pt * z_hat_blk[r, j]
        if (!is.finite(yh_rj)) {
          ipw_cov_transfer_arith_refuse(
            "Point weighted operator row Yhat = sqrt(omega) * z_hat overflowed.",
            reason = "yhat_overflow",
            diagnostics = list(row = r, col = j, omega = w_r, z_hat = z_hat_blk[r, j])
          )
        }
        Yhat[r, j] <- yh_rj
        E_mat[r, j] <- ipw_cov_transfer_iv_point_dist(y_rj_iv, yh_rj)

        t_iv[j, ] <<- gmm_iv_add(t_iv[j, ], gmm_iv_mul(w_iv, z_rj_iv))

        ae_ev <- ae_eval_blk[r, j]
        if (ae_ev != 0.0) {
          ev_sq <- gmm_iv_mul(gmm_iv(ae_ev), gmm_iv(ae_ev))
          sum_eval_sq_iv <<- gmm_iv_add(sum_eval_sq_iv, gmm_iv_mul(w_iv, ev_sq))
        }
        ae_rt <- ae_root_blk[r, j]
        if (ae_rt != 0.0) {
          rt_sq <- gmm_iv_mul(gmm_iv(ae_rt), gmm_iv(ae_rt))
          sum_root_sq_iv <<- gmm_iv_add(sum_root_sq_iv, gmm_iv_mul(w_iv, rt_sq))
        }
        actual_component_ops <<- actual_component_ops + 11
      }
    }

    abs_Yhat <- abs(Yhat)
    gram_pt <- gmm_iv_mat_residual(t(Yhat), Yhat, zero_d_d, max_matrix_entries = max_matrix_entries)
    D1 <- gmm_iv_mat_residual(t(abs_Yhat), E_mat, zero_d_d, max_matrix_entries = max_matrix_entries)$hi
    D2 <- gmm_iv_mat_residual(t(E_mat), abs_Yhat, zero_d_d, max_matrix_entries = max_matrix_entries)$hi
    D3 <- gmm_iv_mat_residual(t(E_mat), E_mat, zero_d_d, max_matrix_entries = max_matrix_entries)$hi
    # Charge worst-case interval primitive work, including the norms returned by
    # all four native residual calls (even though only Gram endpoints are used).
    actual_component_ops <<- actual_component_ops + 8 * m_blk * d * d + 16 * d * d + 8

    for (j in seq_len(d)) {
      for (k in seq_len(d)) {
        d1_jk <- max(0.0, D1[j, k])
        d2_jk <- max(0.0, D2[j, k])
        d3_jk <- max(0.0, D3[j, k])
        D_jk <- gmm_iv_add(gmm_iv_add(gmm_iv(d1_jk), gmm_iv(d2_jk)), gmm_iv(d3_jk))[2L]
        widened_jk <- gmm_iv_add(c(gram_pt$lo[j, k], gram_pt$hi[j, k]), c(-D_jk, D_jk))
        S_iv[j, k, ] <<- gmm_iv_add(S_iv[j, k, ], widened_jk)
        actual_component_ops <<- actual_component_ops + 4
      }
    }
    n_operator_blocks_seen <<- n_operator_blocks_seen + 1L
  }

  # Stream over contiguous source A-unit blocks
  a_start <- 1L
  while (a_start <= n) {
    ipw_cov_transfer_check_interrupt()
    a_end <- as.integer(min(as.double(n), as.double(a_start) + b_in - 1))
    b_len <- a_end - a_start + 1L

    ch_data <- fetch_chunk(a_start, a_end)
    if (!isTRUE(ch_data$certified)) {
      return(list(
        certified = FALSE,
        enclosure_certified = FALSE,
        pd_transferred = FALSE,
        status = "INCONCLUSIVE",
        inconclusive = TRUE,
        inconclusive_type = ch_data$inconclusive_type,
        inconclusive_reasons = ch_data$inconclusive_reasons,
        chunk_diagnostics = ch_data$diagnostics
      ))
    }

    Phi_block <- ch_data$Phi_block
    E_eval_block <- ch_data$E_eval_block
    E_root_block <- ch_data$E_root_block
    actual_4bd_ops <- actual_4bd_ops + ch_data$chunk_4bd_ops
    max_eval_err_seen <- max(max_eval_err_seen, max(E_eval_block))
    max_root_err_seen <- max(max_root_err_seen, max(E_root_block))
    if (!is.null(ch_data$max_direct_star_dist)) {
      max_direct_star_dist_seen <- max(max_direct_star_dist_seen, ch_data$max_direct_star_dist)
    }

    arith_step <- ipw_cov_transfer_capture_arithmetic({
      Uhat <- Phi_block %*% tC
      if (anyNA(Uhat) || any(!is.finite(Uhat))) {
        ipw_cov_transfer_arith_refuse(
          "Point moment transform Uhat = Phi_block %*% t(C) produced nonfinite values.",
          reason = "uhat_nonfinite",
          diagnostics = list(chunk_start = a_start, chunk_end = a_end)
        )
      }
      res_U <- gmm_iv_mat_residual(Phi_block, tC, Uhat, max_matrix_entries = max_matrix_entries)
      U_lo <- matrix(0.0, nrow = b_len, ncol = d)
      U_hi <- matrix(0.0, nrow = b_len, ncol = d)
      for (u in seq_len(b_len)) {
        for (j in seq_len(d)) {
          u_iv <- gmm_iv_add(gmm_iv(Uhat[u, j]), c(res_U$lo[u, j], res_U$hi[u, j]))
          U_lo[u, j] <- u_iv[1L]
          U_hi[u, j] <- u_iv[2L]
        }
      }

      zero_b_d <- matrix(0.0, nrow = b_len, ncol = d)
      res_E_eval <- gmm_iv_mat_residual(E_eval_block, abs_tC, zero_b_d, max_matrix_entries = max_matrix_entries)
      res_E_root <- gmm_iv_mat_residual(E_root_block, abs_tC, zero_b_d, max_matrix_entries = max_matrix_entries)
      ae_eval_blk <- pmax(0.0, res_E_eval$hi)
      dim(ae_eval_blk) <- c(b_len, d)
      ae_root_blk <- pmax(0.0, res_E_root$hi)
      dim(ae_root_blk) <- c(b_len, d)
      actual_component_ops <- actual_component_ops + b_len * (6 * d * L + 16 * d) + 6

      if (is_cluster) {
        c_idx_chunk <- op_norm$unit_to_cluster[seq.int(a_start, a_end)]
        for (u in seq_len(b_len)) {
          cu <- c_idx_chunk[u]
          for (j in seq_len(d)) {
            zh_new <- cluster_z_hat[cu, j] + Uhat[u, j]
            if (!is.finite(zh_new)) {
              ipw_cov_transfer_arith_refuse(
                "Cluster point sum overflowed binary64.",
                reason = "cluster_point_sum_overflow",
                diagnostics = list(cluster_index = cu, moment = j)
              )
            }
            cluster_z_hat[cu, j] <- zh_new
            z_sum_iv <- gmm_iv_add(c(cluster_z_lo[cu, j], cluster_z_hi[cu, j]), c(U_lo[u, j], U_hi[u, j]))
            cluster_z_lo[cu, j] <- z_sum_iv[1L]
            cluster_z_hi[cu, j] <- z_sum_iv[2L]
            if (ae_eval_blk[u, j] != 0.0) {
              cluster_ae_eval[cu, j] <- gmm_iv_add(gmm_iv(cluster_ae_eval[cu, j]), gmm_iv(ae_eval_blk[u, j]))[2L]
            }
            if (ae_root_blk[u, j] != 0.0) {
              cluster_ae_root[cu, j] <- gmm_iv_add(gmm_iv(cluster_ae_root[cu, j]), gmm_iv(ae_root_blk[u, j]))[2L]
            }
            actual_component_ops <- actual_component_ops + 3
          }
        }
      } else {
        # Unit-level operator: stream directly in sub-blocks of size <= b_op
        r_start <- 1L
        while (r_start <= b_len) {
          r_end <- as.integer(min(as.double(b_len), as.double(r_start) + b_op - 1))
          r_seq <- seq.int(r_start, r_end)
          u_seq <- seq.int(a_start + r_start - 1L, a_start + r_end - 1L)
          process_operator_block(
            z_hat_blk = Uhat[r_seq, , drop = FALSE],
            z_lo_blk = U_lo[r_seq, , drop = FALSE],
            z_hi_blk = U_hi[r_seq, , drop = FALSE],
            ae_eval_blk = ae_eval_blk[r_seq, , drop = FALSE],
            ae_root_blk = ae_root_blk[r_seq, , drop = FALSE],
            omega_blk = omega[u_seq]
          )
          r_start <- as.double(r_end) + 1
        }
      }
      NULL
    })

    if (inherits(arith_step, "error")) {
      reason <- if (!is.null(arith_step$reason) && nzchar(arith_step$reason)) arith_step$reason else "stream_gram_arithmetic_refusal"
      return(list(
        certified = FALSE,
        enclosure_certified = FALSE,
        pd_transferred = FALSE,
        status = "INCONCLUSIVE",
        inconclusive = TRUE,
        inconclusive_type = "arithmetic",
        inconclusive_reasons = reason,
        arithmetic_refusal = list(
          occurred = TRUE,
          reason = reason,
          message = conditionMessage(arith_step),
          condition_class = class(arith_step),
          diagnostics = arith_step$diagnostics
        )
      ))
    }

    # End the input-chunk lifetime before the next fetch allocates its buffers.
    # Release aliases in ch_data as well as every derived chunk-sized array.
    rm(ch_data, Phi_block, E_eval_block, E_root_block, Uhat, res_U,
       U_lo, U_hi, zero_b_d, res_E_eval, res_E_root, ae_eval_blk, ae_root_blk)
    if (is_cluster) rm(c_idx_chunk) else rm(r_seq, u_seq)

    n_input_chunks_seen <- n_input_chunks_seen + 1L
    emitted_units <- emitted_units + b_len
    a_start <- as.double(a_end) + 1
  }

  if (emitted_units != n) {
    cdid_ipw_covariance_transfer_abort(
      "Streamed input unit count mismatch during covariance accumulation.",
      "cdid_ipw_covariance_transfer_error"
    )
  }

  # If clustered, traverse completed clusters in declared `cluster_order` in blocks of <= b_op
  if (is_cluster) {
    cluster_arith <- ipw_cov_transfer_capture_arithmetic({
      c_start <- 1L
      while (c_start <= N_c) {
        ipw_cov_transfer_check_interrupt()
        c_end <- as.integer(min(as.double(N_c), as.double(c_start) + b_op - 1))
        c_seq <- seq.int(c_start, c_end)
        process_operator_block(
          z_hat_blk = cluster_z_hat[c_seq, , drop = FALSE],
          z_lo_blk = cluster_z_lo[c_seq, , drop = FALSE],
          z_hi_blk = cluster_z_hi[c_seq, , drop = FALSE],
          ae_eval_blk = cluster_ae_eval[c_seq, , drop = FALSE],
          ae_root_blk = cluster_ae_root[c_seq, , drop = FALSE],
          omega_blk = omega[c_seq]
        )
        c_start <- as.double(c_end) + 1
      }
      NULL
    })
    if (inherits(cluster_arith, "error")) {
      reason <- if (!is.null(cluster_arith$reason) && nzchar(cluster_arith$reason)) cluster_arith$reason else "cluster_gram_arithmetic_refusal"
      return(list(
        certified = FALSE,
        enclosure_certified = FALSE,
        pd_transferred = FALSE,
        status = "INCONCLUSIVE",
        inconclusive = TRUE,
        inconclusive_type = "arithmetic",
        inconclusive_reasons = reason,
        arithmetic_refusal = list(
          occurred = TRUE,
          reason = reason,
          message = conditionMessage(cluster_arith),
          condition_class = class(cluster_arith),
          diagnostics = cluster_arith$diagnostics
        )
      ))
    }
  }

  # The separated cluster-error accumulators are dead after the operator pass.
  if (is_cluster) rm(cluster_ae_eval, cluster_ae_root)

  # Global centering, normalization by kappa, error budget, and PD transfer
  post_arith <- ipw_cov_transfer_capture_arithmetic({
    if (!is.finite(W_iv[1L]) || W_iv[1L] <= 0.0) {
      ipw_cov_transfer_arith_refuse(
        "Total operator weight enclosure W has a nonpositive lower endpoint.",
        reason = "weight_sum_lower_nonpositive",
        diagnostics = list(W_iv = W_iv)
      )
    }

    L_cov <- matrix(0.0, nrow = d, ncol = d, dimnames = dimnames(V_cand_dbl))
    U_cov <- matrix(0.0, nrow = d, ncol = d, dimnames = dimnames(V_cand_dbl))
    dist_mat <- matrix(0.0, nrow = d, ncol = d, dimnames = dimnames(V_cand_dbl))

    for (j in seq_len(d)) {
      for (k in seq_len(d)) {
        s_jk <- S_iv[j, k, ]
        if (centered) {
          tt_jk <- gmm_iv_mul(t_iv[j, ], t_iv[k, ])
          cent_jk <- gmm_iv_div(tt_jk, W_iv)
          s_jk <- gmm_iv_sub(s_jk, cent_jk)
          actual_component_ops <- actual_component_ops + 3
        }
        v_jk <- gmm_iv_mul(kappa_iv, s_jk)
        actual_component_ops <- actual_component_ops + 1
        L_cov[j, k] <- v_jk[1L]
        U_cov[j, k] <- v_jk[2L]
        dist_mat[j, k] <- ipw_cov_transfer_iv_point_dist(v_jk, V_cand_dbl[j, k])
        actual_component_ops <- actual_component_ops + 4
      }
    }
    # Symmetrize interval enclosure conservatively (since exact V is symmetric, intersection is valid;
    # refuse explicitly if floating interval intersection is empty)
    for (j in seq_len(d)) {
      for (k in seq.int(j, d)) {
        lo_sym <- max(L_cov[j, k], L_cov[k, j])
        hi_sym <- min(U_cov[j, k], U_cov[k, j])
        if (lo_sym > hi_sym) {
          ipw_cov_transfer_arith_refuse(
            "Symmetric covariance interval entries have empty intersection.",
            reason = "symmetric_interval_intersection_empty",
            diagnostics = list(j = j, k = k, iv_jk = c(L_cov[j, k], U_cov[j, k]), iv_kj = c(L_cov[k, j], U_cov[k, j]))
          )
        }
        L_cov[j, k] <- L_cov[k, j] <- lo_sym
        U_cov[j, k] <- U_cov[k, j] <- hi_sym
        d_sym <- ipw_cov_transfer_iv_point_dist(c(lo_sym, hi_sym), V_cand_dbl[j, k])
        dist_mat[j, k] <- dist_mat[k, j] <- d_sym
      }
    }

    eps_acc_iv <- gmm_iv_norm(dist_mat, max_matrix_entries = max_matrix_entries)
    epsilon_acc <- eps_acc_iv[2L]
    v_cand_norm_iv <- gmm_iv_norm(V_cand_dbl, max_matrix_entries = max_matrix_entries)
    v_cand_norm_upper <- v_cand_norm_iv[2L]
    v_upper <- gmm_iv_add(gmm_iv(v_cand_norm_upper), gmm_iv(epsilon_acc))[2L]
    sqrt_v_upper <- gmm_iv_sqrt(c(max(0.0, v_upper), max(0.0, v_upper)))[2L]
    two_sqrt_v_upper <- gmm_iv_mul(gmm_iv(2.0), gmm_iv(sqrt_v_upper))[2L]
    actual_component_ops <- actual_component_ops + 10 * d * d + 16

    nu_eval_sq_iv <- gmm_iv_mul(kappa_iv, sum_eval_sq_iv)
    nu_root_sq_iv <- gmm_iv_mul(kappa_iv, sum_root_sq_iv)
    nu_eval <- gmm_iv_sqrt(c(max(0.0, nu_eval_sq_iv[1L]), max(0.0, nu_eval_sq_iv[2L])))[2L]
    nu_root <- gmm_iv_sqrt(c(max(0.0, nu_root_sq_iv[1L]), max(0.0, nu_root_sq_iv[2L])))[2L]
    nu_total <- gmm_iv_add(gmm_iv(nu_eval), gmm_iv(nu_root))[2L]

    linear_eval_term <- gmm_iv_mul(gmm_iv(two_sqrt_v_upper), gmm_iv(nu_eval))[2L]
    linear_root_term <- gmm_iv_mul(gmm_iv(two_sqrt_v_upper), gmm_iv(nu_root))[2L]
    linear_total_term <- max(
      gmm_iv_mul(gmm_iv(two_sqrt_v_upper), gmm_iv(nu_total))[2L],
      gmm_iv_add(gmm_iv(linear_eval_term), gmm_iv(linear_root_term))[2L]
    )

    quad_eval_term <- gmm_iv_mul(gmm_iv(nu_eval), gmm_iv(nu_eval))[2L]
    quad_cross_term <- gmm_iv_mul(gmm_iv(2.0), gmm_iv_mul(gmm_iv(nu_eval), gmm_iv(nu_root)))[2L]
    quad_root_term <- gmm_iv_mul(gmm_iv(nu_root), gmm_iv(nu_root))[2L]
    quad_total_term <- max(
      gmm_iv_mul(gmm_iv(nu_total), gmm_iv(nu_total))[2L],
      gmm_iv_add(
        gmm_iv_add(gmm_iv(quad_eval_term), gmm_iv(quad_cross_term)),
        gmm_iv(quad_root_term)
      )[2L]
    )

    perturbation_total <- gmm_iv_add(gmm_iv(linear_total_term), gmm_iv(quad_total_term))[2L]
    epsilon_V <- gmm_iv_add(gmm_iv(epsilon_acc), gmm_iv(perturbation_total))[2L]
    actual_component_ops <- actual_component_ops + 20

    list(
      L_cov = L_cov,
      U_cov = U_cov,
      dist_mat = dist_mat,
      epsilon_acc = epsilon_acc,
      v_candidate_norm_upper = v_cand_norm_upper,
      v_upper = v_upper,
      sqrt_v_upper = sqrt_v_upper,
      nu_eval = nu_eval,
      nu_root = nu_root,
      nu_total = nu_total,
      linear_eval_term = linear_eval_term,
      linear_root_term = linear_root_term,
      linear_total_term = linear_total_term,
      quad_eval_term = quad_eval_term,
      quad_cross_term = quad_cross_term,
      quad_root_term = quad_root_term,
      quad_total_term = quad_total_term,
      perturbation_total = perturbation_total,
      epsilon_V = epsilon_V
    )
  })

  if (inherits(post_arith, "error")) {
    reason <- if (!is.null(post_arith$reason) && nzchar(post_arith$reason)) post_arith$reason else "post_gram_arithmetic_refusal"
    return(list(
      certified = FALSE,
      enclosure_certified = FALSE,
      pd_transferred = FALSE,
      status = "INCONCLUSIVE",
      inconclusive = TRUE,
      inconclusive_type = "arithmetic",
      inconclusive_reasons = reason,
      arithmetic_refusal = list(
        occurred = TRUE,
        reason = reason,
        message = conditionMessage(post_arith),
        condition_class = class(post_arith),
        diagnostics = post_arith$diagnostics
      )
    ))
  }

  # Spectral certification on V_candidate using the accepted primitive
  eig_try <- tryCatch(
    certify_symmetric_eigen_bounds(V_cand_dbl, label = "V_candidate", max_matrix_entries = max_matrix_entries),
    cdid_gmm_pd_error = function(e) e,
    cdid_gmm_accuracy_error = function(e) e
  )
  actual_component_ops <- actual_component_ops + 6 * d * d * d + 16 * d * d + 64
  if (!is.null(preflight)) {
    stopifnot(
      actual_4bd_ops <= preflight$total_4bd_stream_work_ops_bound,
      actual_component_ops <= preflight$component_work_ops_bound,
      actual_4bd_ops + actual_component_ops <= preflight$total_work_ops_bound
    )
  }

  if (inherits(eig_try, "error") || !isTRUE(eig_try$valid) ||
      !is.finite(eig_try$lambda_min_cert) || !is.finite(eig_try$lambda_max_cert)) {
    reason <- if (inherits(eig_try, "error")) "candidate_eigen_certificate_failed" else "candidate_eigen_orthogonality_inconclusive"
    return(list(
      certified = FALSE,
      enclosure_certified = TRUE,
      pd_transferred = FALSE,
      status = "INCONCLUSIVE",
      inconclusive = TRUE,
      inconclusive_type = "eigen_certification_inconclusive",
      inconclusive_reasons = reason,
      singularity_proved = FALSE,
      nonidentification_proved = FALSE,
      rank_deficiency_proved = FALSE,
      enclosure = list(lo = post_arith$L_cov, hi = post_arith$U_cov, distance_from_candidate = post_arith$dist_mat),
      S_unnormalized_iv = S_iv,
      t_iv = t_iv,
      W_iv = W_iv,
      epsilon_acc = post_arith$epsilon_acc,
      v_candidate_norm_upper = post_arith$v_candidate_norm_upper,
      v_upper = post_arith$v_upper,
      nu_eval = post_arith$nu_eval,
      nu_root = post_arith$nu_root,
      nu_total = post_arith$nu_total,
      epsilon_V = post_arith$epsilon_V,
      error_budget = post_arith[c(
        "epsilon_acc", "v_candidate_norm_upper", "v_upper", "sqrt_v_upper",
        "nu_eval", "nu_root", "nu_total", "linear_eval_term", "linear_root_term",
        "linear_total_term", "quad_eval_term", "quad_cross_term", "quad_root_term",
        "quad_total_term", "perturbation_total", "epsilon_V"
      )],
      eigen_bounds = if (inherits(eig_try, "error")) NULL else eig_try,
      work = list(
        n_input_chunks = logit_cert_count_output(n_input_chunks_seen),
        n_operator_blocks = logit_cert_count_output(n_operator_blocks_seen),
        actual_4bd_interval_ops = actual_4bd_ops,
        actual_component_interval_ops = actual_component_ops,
        total_actual_interval_ops = actual_4bd_ops + actual_component_ops,
        weight_rounding_occurred = weight_rounding_occurred
      )
    ))
  }

  l_min_cert <- eig_try$lambda_min_cert
  l_max_cert <- eig_try$lambda_max_cert
  eps_V <- post_arith$epsilon_V
  eps_acc <- post_arith$epsilon_acc

  lambda_min_star_lower <- tryCatch(
    logit_b64_sub(gmm_iv(l_min_cert), gmm_iv(eps_V))[1L],
    error = function(e) -Inf
  )
  lambda_max_star_upper <- tryCatch(
    logit_b64_add(gmm_iv(l_max_cert), gmm_iv(eps_V))[2L],
    error = function(e) Inf
  )
  lambda_min_stored_lower <- tryCatch(
    logit_b64_sub(gmm_iv(l_min_cert), gmm_iv(eps_acc))[1L],
    error = function(e) -Inf
  )
  lambda_max_stored_upper <- tryCatch(
    logit_b64_add(gmm_iv(l_max_cert), gmm_iv(eps_acc))[2L],
    error = function(e) Inf
  )

  pd_ok <- is.finite(l_min_cert) && (l_min_cert > 0.0) &&
    is.finite(eps_V) && (eps_V < l_min_cert) &&
    is.finite(lambda_min_star_lower) && (lambda_min_star_lower > 0.0)

  inc_reasons <- character(0L)
  inc_type <- "none"
  if (!pd_ok) {
    inc_type <- "pd_transfer_inconclusive"
    inc_reasons <- if (l_min_cert <= 0.0) {
      "candidate_lambda_min_cert_nonpositive"
    } else {
      "epsilon_V_not_less_than_lambda_min_cert"
    }
  }

  # For finite endpoints, a row has zero sum of absolute endpoint magnitudes
  # exactly when every endpoint is zero (including either signed zero).
  # Scalar traversal avoids full N_c x d diagnostic temporaries.
  n_zero_if_clusters <- 0L
  if (is_cluster) {
    cu <- 1.0
    while (cu <= N_c) {
      if ((cu - 1) %% b_op == 0) ipw_cov_transfer_check_interrupt()
      row_is_zero <- TRUE
      j <- 1.0
      while (j <= d) {
        if (cluster_z_lo[cu, j] != 0.0 || cluster_z_hi[cu, j] != 0.0) {
          row_is_zero <- FALSE
          break
        }
        j <- j + 1.0
      }
      if (row_is_zero) n_zero_if_clusters <- n_zero_if_clusters + 1L
      cu <- cu + 1.0
    }
  }

  list(
    certified = pd_ok,
    enclosure_certified = TRUE,
    pd_transferred = pd_ok,
    status = if (pd_ok) "CERTIFIED" else "INCONCLUSIVE",
    inconclusive = !pd_ok,
    inconclusive_type = inc_type,
    inconclusive_reasons = inc_reasons,
    singularity_proved = FALSE,
    nonidentification_proved = FALSE,
    rank_deficiency_proved = FALSE,
    enclosure = list(
      lo = post_arith$L_cov,
      hi = post_arith$U_cov,
      distance_from_candidate = post_arith$dist_mat
    ),
    S_unnormalized_iv = S_iv,
    t_iv = t_iv,
    W_iv = W_iv,
    epsilon_acc = eps_acc,
    v_candidate_norm_upper = post_arith$v_candidate_norm_upper,
    v_upper = post_arith$v_upper,
    nu_eval = post_arith$nu_eval,
    nu_root = post_arith$nu_root,
    nu_total = post_arith$nu_total,
    epsilon_V = eps_V,
    error_budget = post_arith[c(
      "epsilon_acc", "v_candidate_norm_upper", "v_upper", "sqrt_v_upper",
      "nu_eval", "nu_root", "nu_total", "linear_eval_term", "linear_root_term",
      "linear_total_term", "quad_eval_term", "quad_cross_term", "quad_root_term",
      "quad_total_term", "perturbation_total", "epsilon_V"
    )],
    eigen_bounds = list(
      candidate = eig_try,
      lambda_min_candidate_cert = l_min_cert,
      lambda_max_candidate_cert = l_max_cert,
      lambda_min_stored_lower = lambda_min_stored_lower,
      lambda_max_stored_upper = lambda_max_stored_upper,
      lambda_min_star_lower = lambda_min_star_lower,
      lambda_max_star_upper = lambda_max_star_upper
    ),
    row_error_summary = list(
      max_eval_error = max_eval_err_seen,
      max_root_error = max_root_err_seen,
      max_direct_star_interval_distance_not_root_radius = max_direct_star_dist_seen
    ),
    cluster_summary = if (is_cluster) {
      list(
        n_declared_clusters = N_c,
        n_zero_unit_declared_clusters = op_norm$n_zero_unit_declared_clusters,
        n_zero_if_clusters = n_zero_if_clusters,
        cluster_point = cluster_z_hat,
        cluster_interval_lo = cluster_z_lo,
        cluster_interval_hi = cluster_z_hi
      )
    } else NULL,
    work = list(
      n_input_chunks = logit_cert_count_output(n_input_chunks_seen),
      n_operator_blocks = logit_cert_count_output(n_operator_blocks_seen),
      actual_4bd_interval_ops = actual_4bd_ops,
      actual_component_interval_ops = actual_component_ops,
      total_actual_interval_ops = actual_4bd_ops + actual_component_ops,
      weight_rounding_occurred = weight_rounding_occurred
    )
  )
}

# ---------------------------------------------------------------------------
# Required Unexported Entry Point:
#   certify_ipw_covariance_transfer(links, candidate_rows, frame, C, op, V_candidate, control)
# ---------------------------------------------------------------------------

certify_ipw_covariance_transfer <- function(
  links,
  candidate_rows,
  frame,
  C,
  op,
  V_candidate,
  control
) {
  ipw_cov_transfer_require_deps()
  if (missing(links) || missing(candidate_rows) || missing(frame) ||
      missing(C) || missing(op) || missing(V_candidate) || missing(control)) {
    cdid_ipw_covariance_transfer_abort(
      "All seven required arguments (`links`, `candidate_rows`, `frame`, `C`, `op`, `V_candidate`, `control`) must be supplied explicitly.",
      "cdid_ipw_covariance_transfer_input_error"
    )
  }

  # 1. Validate control and lightweight headers BEFORE scanning content or calling callbacks
  ctrl <- validate_ipw_covariance_transfer_control(control)
  # Bound the header vectors before allocating/scanning one header per link.
  # The fuller preflight below also accounts for retained per-link summaries.
  ipw_cov_transfer_safe_product(66, length(links), ctrl$max_matrix_entries, "link header vectors and per-link summaries before header validation")
  link_hdr <- ipw_cov_transfer_validate_links_headers(links)
  n <- link_hdr$n
  L <- link_hdr$L
  expected_link_ids <- link_hdr$link_ids

  ipw_cov_transfer_refuse_decorative_attrs(C, "C", allowed_attrs = c("dim", "dimnames"))
  if (is.null(C) || is.object(C) || !is.atomic(C) || !is.numeric(C) ||
      !is.matrix(C) || nrow(C) < 1L || ncol(C) != L) {
    cdid_ipw_covariance_transfer_abort(
      paste0("`C` must be a plain numeric d x L matrix with L = ", L, " columns and d >= 1 rows."),
      c("cdid_ipw_covariance_transfer_map_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  d <- nrow(C)

  ipw_cov_transfer_refuse_decorative_attrs(V_candidate, "V_candidate", allowed_attrs = c("dim", "dimnames"))
  if (is.null(V_candidate) || is.object(V_candidate) || !is.atomic(V_candidate) ||
      !is.numeric(V_candidate) || !is.matrix(V_candidate) ||
      nrow(V_candidate) != d || ncol(V_candidate) != d) {
    cdid_ipw_covariance_transfer_abort(
      paste0("`V_candidate` must be a plain numeric d x d matrix (d = ", d, ")."),
      c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }

  frame_spec <- ipw_cov_transfer_validate_frame_header(frame, n = n, id_type_ref = link_hdr$id_type)
  op_spec <- ipw_cov_transfer_validate_operator_header(op, n = n)

  ipw_cov_transfer_refuse_decorative_attrs(
    candidate_rows,
    "candidate_rows",
    allowed_attrs = if (is.function(candidate_rows)) "srcref" else c("names", "class")
  )
  cand_spec <- if (is.function(candidate_rows)) {
    ipw_candidate_rows_spec(callback = candidate_rows, link_ids = NULL)
  } else if (is.list(candidate_rows)) {
    if (is.null(names(candidate_rows)) || any(!nzchar(names(candidate_rows))) || anyDuplicated(names(candidate_rows))) {
      cdid_ipw_covariance_transfer_abort(
        "`candidate_rows` list must have unique non-empty field names.",
        c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
    forbidden_cand <- c("certified", "trusted", "is_exact", "bypass_enclosure", "trust_metadata")
    if (any(names(candidate_rows) %in% forbidden_cand)) {
      cdid_ipw_covariance_transfer_abort(
        "Shortcut or certificate flags on `candidate_rows` are refused.",
        c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_metadata_error")
      )
    }
    extra_cand <- setdiff(names(candidate_rows), c("callback", "link_ids"))
    if (length(extra_cand) > 0L) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Unexpected fields (`", paste(extra_cand, collapse = "`, `"), "`) in `candidate_rows`."),
        c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
    ipw_candidate_rows_spec(callback = candidate_rows$callback, link_ids = candidate_rows$link_ids)
  } else {
    cdid_ipw_covariance_transfer_abort(
      "`candidate_rows` must be a callback function `function(a_start, a_end)` or an `ipw_candidate_rows_spec` list.",
      c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (!is.null(cand_spec$link_ids) && !identical(as.character(cand_spec$link_ids), expected_link_ids)) {
    cdid_ipw_covariance_transfer_abort(
      "`candidate_rows$link_ids` disagrees with the ordered `link_id` sequence of `links`.",
      c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }

  # 2. Preflight BEFORE content scans, character byte scans, copies, or callbacks
  n_clusters_hdr <- if (identical(op_spec$row_level, "cluster")) length(op_spec$cluster_order) else n
  preflight <- preflight_ipw_covariance_transfer(
    n = n,
    L = L,
    d = d,
    row_level = op_spec$row_level,
    n_clusters = n_clusters_hdr,
    p_vec = link_hdr$p_vec,
    m_vec = link_hdr$m_vec,
    certificate_controls = link_hdr$certificate_controls,
    control = ctrl,
    id_is_character = identical(link_hdr$id_type, "character"),
    cluster_is_character = !is.null(op_spec$cluster_id) && is.character(op_spec$cluster_id),
    existing_4bd_retained_entries = link_hdr$existing_4bd_retained_entries
  )
  b_in <- preflight$effective_chunk_rows
  b_op <- preflight$effective_operator_block_rows
  max_me <- ctrl$max_matrix_entries

  # 3. Validate C, V_candidate, links content, frame content, and operator content
  if (anyNA(C) || any(!is.finite(C))) {
    cdid_ipw_covariance_transfer_abort(
      "`C` must be a finite numeric d x L matrix without NA/NaN/Inf.",
      c("cdid_ipw_covariance_transfer_map_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (!is.null(colnames(C)) && !identical(colnames(C), expected_link_ids)) {
    cdid_ipw_covariance_transfer_abort(
      "`colnames(C)` disagrees with the ordered `link_id` sequence of `links`.",
      c("cdid_ipw_covariance_transfer_map_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  C_dbl <- matrix(as.double(C), nrow = d, ncol = L, dimnames = list(rownames(C), expected_link_ids))

  if (anyNA(V_candidate) || any(!is.finite(V_candidate))) {
    cdid_ipw_covariance_transfer_abort(
      "`V_candidate` must be a finite numeric d x d matrix without NA/NaN/Inf.",
      c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  for (j in seq_len(d)) {
    if (!identical(as.double(V_candidate[, j]), as.double(V_candidate[j, ]))) {
      cdid_ipw_covariance_transfer_abort(
        "`V_candidate` must be an exactly symmetric numeric d x d matrix.",
        c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  if (!is.null(rownames(C))) {
    if ((!is.null(rownames(V_candidate)) && !identical(rownames(V_candidate), rownames(C))) ||
        (!is.null(colnames(V_candidate)) && !identical(colnames(V_candidate), rownames(C)))) {
      cdid_ipw_covariance_transfer_abort(
        "`V_candidate` row/column names disagree with `rownames(C)`.",
        c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  V_cand_dbl <- matrix(as.double(V_candidate), nrow = d, ncol = d, dimnames = dimnames(V_candidate))

  # Account for link_ids character temporary bytes
  link_char_info <- ipw_cov_transfer_account_char_bytes(
    vec = expected_link_ids,
    n_len = L,
    step = max(1L, b_in),
    label = "links$link_id",
    max_matrix_entries = max_me,
    base_concurrent_entries = preflight$base_component_concurrent_entries,
    char_bytes_so_far = 0
  )

  ipw_cov_transfer_validate_links_content(links, link_hdr)

  raw_1 <- link_hdr$problems[[1L]]$raw
  frame_val <- ipw_cov_transfer_validate_frame_content(
    frame_spec = frame_spec,
    raw_ref = raw_1,
    n = n,
    step = b_in,
    max_matrix_entries = max_me,
    base_concurrent_entries = preflight$base_component_concurrent_entries,
    char_bytes_so_far = link_char_info$total_char_temp_bytes
  )

  op_norm <- ipw_cov_transfer_validate_operator_content(
    op_spec = op_spec,
    n = n,
    step = b_in,
    max_matrix_entries = max_me,
    base_concurrent_entries = preflight$base_component_concurrent_entries,
    char_bytes_so_far = frame_val$char_info$total_char_temp_bytes
  )

  # Update preflight character temporary accounting with actual validated static string bytes
  static_char_bytes <- op_norm$char_info$total_char_temp_bytes
  preflight$character_temporary_bytes <- static_char_bytes
  preflight$character_temporary_entries_8byte <- ceiling(static_char_bytes / 8)
  preflight$concurrent_entries <- ipw_cov_transfer_safe_sum(
    c(preflight$base_component_concurrent_entries, preflight$character_temporary_entries_8byte),
    max_bound = max_me,
    label = "concurrent entries after frame/operator character byte accounting"
  )
  preflight$concurrent_bytes_8double <- 8 * preflight$concurrent_entries
  preflight$total_with_existing_4bd_entries <- ipw_cov_transfer_safe_sum(
    c(preflight$existing_4bd_retained_entries, preflight$concurrent_entries),
    label = "validated existing 4bd plus component concurrent payload"
  )
  preflight$total_with_existing_4bd_bytes_8double <- 8 * preflight$total_with_existing_4bd_entries

  # 4. Check whether any sealed source nuisance problem or link is INCONCLUSIVE
  uncert_links <- which(!vapply(links, function(lk) isTRUE(lk$certified), logical(1L)))
  uncert_probs <- which(!vapply(link_hdr$problems, function(pr) isTRUE(pr$certified), logical(1L)))
  if (length(uncert_links) > 0L || length(uncert_probs) > 0L) {
    reasons <- unique(unlist(lapply(seq_len(L), function(l) {
      c(
        if (!isTRUE(link_hdr$problems[[l]]$certified)) link_hdr$problems[[l]]$inconclusive_reasons else character(0L),
        if (!isTRUE(links[[l]]$certified)) links[[l]]$inconclusive_reasons else character(0L)
      )
    })))
    if (length(reasons) == 0L) reasons <- "source_link_or_problem_inconclusive"
    inc_obj <- structure(
      list(
        certified = FALSE,
        enclosure_certified = FALSE,
        pd_transferred = FALSE,
        status = "INCONCLUSIVE",
        inconclusive = TRUE,
        inconclusive_type = "source_link_inconclusive",
        inconclusive_reasons = reasons,
        separation_certified = FALSE,
        sampling_certified = FALSE,
        identification_certified = FALSE,
        raw_panel_eligibility_certified = FALSE,
        public_gate_adopted = FALSE,
        singularity_proved = FALSE,
        nonidentification_proved = FALSE,
        rank_deficiency_proved = FALSE,
        link_ids = expected_link_ids,
        n_a_units = as.double(n),
        n_links = as.double(L),
        n_moments = as.double(d),
        frame = list(
          normalizer_n = as.double(n),
          embedding = "identity",
          id_type = link_hdr$id_type,
          id_aliased_to_source = frame_val$id_aliased_to_source
        ),
        operator = list(
          purpose = op_norm$purpose,
          row_level = op_norm$row_level,
          n_clusters = as.double(op_norm$n_clusters),
          n_positive_weight_rows = as.double(op_norm$n_positive_weight_rows),
          centered = op_norm$centered,
          kappa = op_norm$kappa
        ),
        C = C_dbl,
        V_candidate = V_cand_dbl,
        enclosure = NULL,
        epsilon_acc = NA_real_,
        v_upper = NA_real_,
        nu_eval = NA_real_,
        nu_root = NA_real_,
        nu_total = NA_real_,
        epsilon_V = NA_real_,
        error_budget = NULL,
        eigen_bounds = NULL,
        control = ctrl,
        preflight = preflight
      ),
      class = "cdid_ipw_covariance_transfer"
    )
    return(ipw_enclosure_seal(inc_obj, "covariance_transfer"))
  }

  # 5. Callback chunk fetcher that validates batch headers BEFORE copying and calls enclose_ipw_link_if_chunk
  cb_fn <- cand_spec$callback
  is_num_id <- identical(link_hdr$id_type, "numeric")
  max_chunk_char_bytes <- 0

  fetch_chunk_from_links <- function(a_start, a_end) {
    b_len <- a_end - a_start + 1L
    a_seq <- seq.int(a_start, a_end)

    batch <- tryCatch(
      withCallingHandlers(
        cb_fn(a_start, a_end),
        warning = function(w) {
          cdid_ipw_covariance_transfer_abort(
            paste0("`candidate_rows` callback emitted a warning: ", conditionMessage(w)),
            "cdid_ipw_covariance_transfer_callback_error",
            diagnostics = list(a_start = a_start, a_end = a_end)
          )
        }
      ),
      error = function(e) {
        if (inherits(e, "cdid_ipw_covariance_transfer_refusal") ||
            inherits(e, "cdid_ipw_enclosure_refusal")) {
          stop(e)
        }
        cdid_ipw_covariance_transfer_abort(
          paste0("`candidate_rows` callback failed: ", conditionMessage(e)),
          "cdid_ipw_covariance_transfer_callback_error",
          diagnostics = list(a_start = a_start, a_end = a_end)
        )
      }
    )

    # Refuse oversized or malformed callback outputs BEFORE any component copies/materialization
    if (is.null(batch) || !is.list(batch) || is.object(batch)) {
      cdid_ipw_covariance_transfer_abort(
        "`candidate_rows` callback must return a plain named list (`phi`, `id`, `link_ids`).",
        "cdid_ipw_covariance_transfer_callback_error",
        diagnostics = list(a_start = a_start, a_end = a_end)
      )
    }
    ipw_cov_transfer_refuse_decorative_attrs(batch, "candidate_rows batch", allowed_attrs = "names")
    if (is.null(names(batch)) || any(!nzchar(names(batch))) || anyDuplicated(names(batch))) {
      cdid_ipw_covariance_transfer_abort(
        "`candidate_rows` callback batch must have unique non-empty names.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    forbidden_batch <- c("certified", "trusted", "is_exact", "bypass_enclosure", "trust_metadata")
    if (any(names(batch) %in% forbidden_batch)) {
      cdid_ipw_covariance_transfer_abort(
        "Decorative certificate flags on `candidate_rows` callback output are refused.",
        c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_metadata_error")
      )
    }
    if (!all(c("phi", "id", "link_ids") %in% names(batch))) {
      cdid_ipw_covariance_transfer_abort(
        "`candidate_rows` callback batch must contain `phi`, `id`, and `link_ids`.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    extra_b <- setdiff(names(batch), c("phi", "id", "link_ids", "a_index"))
    if (length(extra_b) > 0L) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Unexpected fields (`", paste(extra_b, collapse = "`, `"), "`) in `candidate_rows` callback batch."),
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }

    phi_raw <- batch$phi
    id_raw <- batch$id
    lids_raw <- batch$link_ids
    ipw_cov_transfer_refuse_decorative_attrs(phi_raw, "batch$phi", allowed_attrs = c("dim", "dimnames"))
    ipw_cov_transfer_refuse_decorative_attrs(id_raw, "batch$id")
    ipw_cov_transfer_refuse_decorative_attrs(lids_raw, "batch$link_ids")

    if (is.null(phi_raw) || is.object(phi_raw) || !is.atomic(phi_raw) ||
        !is.numeric(phi_raw) || !is.matrix(phi_raw)) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$phi` from `candidate_rows` must be a plain numeric matrix.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    nr_p <- nrow(phi_raw)
    nc_p <- ncol(phi_raw)
    if (is.na(nr_p) || is.na(nc_p) || nr_p > ctrl$max_chunk_rows ||
        as.double(nr_p) * as.double(nc_p) > max_me) {
      cdid_ipw_covariance_transfer_abort(
        "Oversized `batch$phi` returned by `candidate_rows` callback exceeds explicit chunk/matrix bounds (refused before copy).",
        c("cdid_ipw_covariance_transfer_bound_error", "cdid_ipw_covariance_transfer_callback_error"),
        diagnostics = list(nrow = nr_p, ncol = nc_p, expected_rows = b_len, expected_cols = L)
      )
    }
    if (nr_p != b_len || nc_p != L) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`batch$phi` dimensions (", nr_p, " x ", nc_p, ") disagree with requested chunk (", b_len, " x ", L, ")."),
        "cdid_ipw_covariance_transfer_callback_error",
        diagnostics = list(nrow = nr_p, ncol = nc_p, expected_rows = b_len, expected_cols = L)
      )
    }

    if (is.null(id_raw) || is.factor(id_raw) || is.object(id_raw) || !is.atomic(id_raw) ||
        !is.null(dim(id_raw)) || length(id_raw) != b_len || is.numeric(id_raw) != is_num_id ||
        !(is.numeric(id_raw) || is.character(id_raw))) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$id` from `candidate_rows` must be a plain vector of length `b_len` matching the source ID type.",
        c("cdid_ipw_covariance_transfer_id_error", "cdid_ipw_covariance_transfer_callback_error")
      )
    }
    if (is.null(lids_raw) || is.object(lids_raw) || !is.atomic(lids_raw) ||
        !is.null(dim(lids_raw)) || !is.character(lids_raw) || length(lids_raw) != L) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$link_ids` from `candidate_rows` must be a plain character vector of length `L`.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    if (!is.null(batch$a_index)) {
      if (is.object(batch$a_index) || !is.atomic(batch$a_index) || !is.null(dim(batch$a_index)) ||
          !is.numeric(batch$a_index) || length(batch$a_index) != b_len ||
          !identical(as.integer(batch$a_index), a_seq)) {
        cdid_ipw_covariance_transfer_abort(
          "`batch$a_index` from `candidate_rows` disagrees with requested `[a_start, a_end]`.",
          "cdid_ipw_covariance_transfer_callback_error"
        )
      }
    }

    # Account for callback character temporary bytes before string comparisons
    chunk_char_bytes <- 0
    if (!is_num_id) {
      ch_id_info <- ipw_cov_transfer_account_char_bytes(
        vec = id_raw,
        n_len = b_len,
        step = b_len,
        label = "batch$id",
        max_matrix_entries = max_me,
        base_concurrent_entries = preflight$base_component_concurrent_entries,
        char_bytes_so_far = static_char_bytes
      )
      chunk_char_bytes <- ch_id_info$total_char_temp_bytes - static_char_bytes
    }
    ch_lid_info <- ipw_cov_transfer_account_char_bytes(
      vec = lids_raw,
      n_len = L,
      step = L,
      label = "batch$link_ids",
      max_matrix_entries = max_me,
      base_concurrent_entries = preflight$base_component_concurrent_entries,
      char_bytes_so_far = static_char_bytes + chunk_char_bytes
    )
    chunk_char_bytes <- ch_lid_info$total_char_temp_bytes - static_char_bytes
    if (chunk_char_bytes > max_chunk_char_bytes) {
      max_chunk_char_bytes <<- chunk_char_bytes
    }

    if (!identical(as.character(lids_raw), expected_link_ids)) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$link_ids` from `candidate_rows` disagrees with the ordered `link_id` sequence of `links`.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    if (!is.null(colnames(phi_raw)) && !identical(colnames(phi_raw), expected_link_ids)) {
      cdid_ipw_covariance_transfer_abort(
        "`colnames(batch$phi)` disagrees with `batch$link_ids`.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }

    batch_id_key <- if (is_num_id) {
      if (anyNA(id_raw) || any(!is.finite(id_raw))) {
        cdid_ipw_covariance_transfer_abort(
          "`batch$id` contains NA, NaN, or infinite numeric identifiers.",
          c("cdid_ipw_covariance_transfer_id_error", "cdid_ipw_covariance_transfer_callback_error")
        )
      }
      d_id <- as.double(id_raw)
      d_id[d_id == 0.0] <- 0.0
      d_id
    } else {
      as.character(id_raw)
    }
    if (!identical(batch_id_key, raw_1$id_key[a_seq])) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$id` from `candidate_rows` disagrees with source A-frame unit IDs on `[a_start, a_end]`.",
        c("cdid_ipw_covariance_transfer_id_error", "cdid_ipw_covariance_transfer_callback_error"),
        diagnostics = list(a_start = a_start, a_end = a_end)
      )
    }

    if (anyNA(phi_raw) || any(!is.finite(phi_raw))) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$phi` from `candidate_rows` contains NA, NaN, or Inf.",
        c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_callback_error"),
        diagnostics = list(a_start = a_start, a_end = a_end)
      )
    }
    Phi_block <- if (typeof(phi_raw) == "double") phi_raw else matrix(as.double(phi_raw), nrow = b_len, ncol = L)

    E_eval_block <- matrix(0.0, nrow = b_len, ncol = L)
    E_root_block <- matrix(0.0, nrow = b_len, ncol = L)
    chunk_4bd_ops <- 0
    max_direct_star_dist <- 0.0

    for (l in seq_len(L)) {
      phi_col_l <- as.double(Phi_block[, l])
      ch_l <- enclose_ipw_link_if_chunk(
        link_enclosure = links[[l]],
        a_start = a_start,
        a_end = a_end,
        control = ctrl$link_control,
        phi_candidate = phi_col_l,
        include_components = FALSE,
        include_shared_xi = FALSE
      )
      if (!isTRUE(ch_l$certified)) {
        return(list(
          certified = FALSE,
          inconclusive_type = "source_chunk_inconclusive",
          inconclusive_reasons = unique(c("link_if_chunk_inconclusive", ch_l$inconclusive_reasons)),
          diagnostics = list(link_index = l, link_id = expected_link_ids[l], a_start = a_start, a_end = a_end)
        ))
      }
      sep_rows <- ch_l$separated_row_errors
      if (is.null(sep_rows) ||
          !identical(ch_l$a_index, a_seq) ||
          !identical(sep_rows$phi_candidate, phi_col_l) ||
          length(sep_rows$eval_error) != b_len ||
          length(sep_rows$root_error) != b_len) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`enclose_ipw_link_if_chunk` output for link ", l, " lacks verified `separated_row_errors`."),
          "cdid_ipw_covariance_transfer_integrity_error"
        )
      }
      # Never relabel direct_star_error as root_error:
      # eval_error bounds discrepancy from phi(beta0); root_error bounds change from beta0 to beta*.
      E_eval_block[, l] <- sep_rows$eval_error
      E_root_block[, l] <- sep_rows$root_error
      if (!is.null(sep_rows$direct_star_error)) {
        max_direct_star_dist <- max(max_direct_star_dist, max(sep_rows$direct_star_error))
      }
      chunk_4bd_ops <- chunk_4bd_ops + ch_l$work$actual_interval_ops
      # No previous link result/aliases may survive into the next 4bd call.
      rm(ch_l, sep_rows, phi_col_l)
    }

    list(
      certified = TRUE,
      Phi_block = Phi_block,
      E_eval_block = E_eval_block,
      E_root_block = E_root_block,
      chunk_4bd_ops = chunk_4bd_ops,
      max_direct_star_dist = max_direct_star_dist
    )
  }

  # 6. Run the streaming Gram and covariance-transfer kernel
  kern <- ipw_cov_transfer_stream_kernel(
    n = n,
    L = L,
    d = d,
    b_in = b_in,
    b_op = b_op,
    fetch_chunk = fetch_chunk_from_links,
    C_dbl = C_dbl,
    op_norm = op_norm,
    V_cand_dbl = V_cand_dbl,
    max_matrix_entries = max_me,
    preflight = preflight
  )

  total_char_bytes <- static_char_bytes + max_chunk_char_bytes
  preflight$character_temporary_bytes <- total_char_bytes
  preflight$character_temporary_entries_8byte <- ceiling(total_char_bytes / 8)
  preflight$concurrent_entries <- ipw_cov_transfer_safe_sum(
    c(preflight$base_component_concurrent_entries, preflight$character_temporary_entries_8byte),
    max_bound = max_me,
    label = "final concurrent entries including callback character temporary bytes"
  )
  preflight$concurrent_bytes_8double <- 8 * preflight$concurrent_entries
  preflight$total_with_existing_4bd_entries <- ipw_cov_transfer_safe_sum(
    c(preflight$existing_4bd_retained_entries, preflight$concurrent_entries),
    label = "validated existing 4bd plus component concurrent payload"
  )
  preflight$total_with_existing_4bd_bytes_8double <- 8 * preflight$total_with_existing_4bd_entries

  res_obj <- structure(
    list(
      certified = isTRUE(kern$certified),
      enclosure_certified = isTRUE(kern$enclosure_certified),
      pd_transferred = isTRUE(kern$pd_transferred),
      status = kern$status,
      inconclusive = isTRUE(kern$inconclusive),
      inconclusive_type = kern$inconclusive_type,
      inconclusive_reasons = kern$inconclusive_reasons,
      separation_certified = FALSE,
      sampling_certified = FALSE,
      identification_certified = FALSE,
      raw_panel_eligibility_certified = FALSE,
      public_gate_adopted = FALSE,
      singularity_proved = FALSE,
      nonidentification_proved = FALSE,
      rank_deficiency_proved = FALSE,
      link_ids = expected_link_ids,
      n_a_units = as.double(n),
      n_links = as.double(L),
      n_moments = as.double(d),
      n_unique_nuisance_problems = as.integer(link_hdr$n_unique_problems),
      frame = list(
        normalizer_n = as.double(n),
        embedding = "identity",
        id_type = link_hdr$id_type,
        id_aliased_to_source = frame_val$id_aliased_to_source,
        larger_frame_transport_applied = FALSE,
        a_rows_rescaled = FALSE,
        zero_rows_invented = FALSE
      ),
      operator = list(
        purpose = op_norm$purpose,
        row_level = op_norm$row_level,
        n_clusters = as.double(op_norm$n_clusters),
        n_positive_weight_rows = as.double(op_norm$n_positive_weight_rows),
        n_zero_unit_declared_clusters = as.integer(op_norm$n_zero_unit_declared_clusters),
        centered = op_norm$centered,
        kappa = op_norm$kappa,
        case_weights_reapplied = FALSE
      ),
      C = C_dbl,
      V_candidate = V_cand_dbl,
      enclosure = kern$enclosure,
      S_unnormalized_iv = kern$S_unnormalized_iv,
      t_iv = kern$t_iv,
      W_iv = kern$W_iv,
      epsilon_acc = if (!is.null(kern$epsilon_acc)) kern$epsilon_acc else NA_real_,
      v_candidate_norm_upper = if (!is.null(kern$v_candidate_norm_upper)) kern$v_candidate_norm_upper else NA_real_,
      v_upper = if (!is.null(kern$v_upper)) kern$v_upper else NA_real_,
      nu_eval = if (!is.null(kern$nu_eval)) kern$nu_eval else NA_real_,
      nu_root = if (!is.null(kern$nu_root)) kern$nu_root else NA_real_,
      nu_total = if (!is.null(kern$nu_total)) kern$nu_total else NA_real_,
      epsilon_V = if (!is.null(kern$epsilon_V)) kern$epsilon_V else NA_real_,
      error_budget = kern$error_budget,
      eigen_bounds = kern$eigen_bounds,
      row_error_summary = kern$row_error_summary,
      cluster_summary = kern$cluster_summary,
      arithmetic_refusal = if (!is.null(kern$arithmetic_refusal)) {
        kern$arithmetic_refusal
      } else {
        list(occurred = FALSE, reason = "none", message = NULL, condition_class = NULL, diagnostics = NULL)
      },
      chunk_diagnostics = kern$chunk_diagnostics,
      control = ctrl,
      preflight = preflight,
      work = kern$work,
      provenance = list(
        source_link_refs = lapply(links, `[[`, "ref"),
        source_problem_refs = lapply(link_hdr$problems, `[[`, "ref"),
        sep_rows_semantics = "eval_error bounds discrepancy from phi(beta0); root_error is the Theorem B E_root bound from beta0 to beta*; direct_star_error is never relabeled as root_error.",
        conditional_scope = "Certifies the finite covariance functional and sufficient Weyl PD transfer conditional on supplied sealed 4bd link enclosures in the identical complete A-frame. Does not validate sampling, raw panel eligibility, missing-trends/PT/overlap, population identification, larger-frame transport, estimated aggregation-weight IF, GMM/multiplier transfer, or coverage."
      )
    ),
    class = "cdid_ipw_covariance_transfer"
  )

  ipw_enclosure_seal(res_obj, "covariance_transfer")
}

validate_ipw_covariance_transfer <- function(transfer, require_certified = FALSE) {
  ipw_cov_transfer_require_deps()
  if (missing(transfer) || is.null(transfer) ||
      !inherits(transfer, "cdid_ipw_covariance_transfer") ||
      !is.list(transfer) || !ipw_enclosure_seal_ok(transfer, "covariance_transfer")) {
    cdid_ipw_covariance_transfer_abort(
      "Covariance transfer certificate is malformed, foreign, or differs from its sealed snapshot.",
      c("cdid_ipw_covariance_transfer_integrity_error", "cdid_ipw_covariance_transfer_error")
    )
  }
  if (isTRUE(require_certified) && !isTRUE(transfer$certified)) {
    cdid_ipw_covariance_transfer_abort(
      "Certified covariance transfer is required, but transfer status is INCONCLUSIVE.",
      c("cdid_ipw_covariance_transfer_inconclusive", "cdid_ipw_covariance_transfer_error"),
      diagnostics = list(inconclusive_reasons = transfer$inconclusive_reasons)
    )
  }
  invisible(transfer)
}
