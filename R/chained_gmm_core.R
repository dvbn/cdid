# Private CP02 exact-GMM computation core.
#
# Implements the patched conditional private contract (docs/theory/cp02_gmm_computation.tex
# and docs/theory/cp02_gmm_computation_contract.md):
#   * deterministic unanchored component gauge (lowest parameter index per
#     unanchored component, including isolated vertices) without dropping rows;
#   * exact IEEE 754 binary64 dyadic zero-sum check for target estimability;
#   * full joint solve across all cohorts and components via unpivoted Cholesky
#     whitening and LAPACK column-pivoted QR (no rank truncation, no inverse,
#     no normal-equations solve, no fallback);
#   * explicit identity criterion when requested, with streamed sandwich SEs
#     and no structural reduction;
#   * conservative sufficient numerical guard for positive definiteness,
#     factorization residuals, conditioning, and certified a posteriori target
#     error bounds (including the approximate-adjoint residual term
#     ||s||_2 ||K^{-1}||_2 ||r||_2 and IEEE binary64 roundoff margins);
#   * bounded two-pass streamed IF projection with cluster accumulation across
#     batch boundaries and weighted centering;
#   * explicit caller controls/bounds/provenance with no internal defaults.
# The D26 conventional strategy adds reconstruction-aware covariance ridge
# after proved exact reductions (R/chained_gmm_weighting.R). The ordinary
# unregularized calculations and private interval route remain explicit.

cdid_gmm_abort <- function(message, class, diagnostics = list()) {
  stop(structure(
    class = unique(c(class, "cdid_gmm_refusal", "cdid_gmm_error", "cdid_error", "error", "condition")),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

gmm_scalar_logical <- function(x) {
  is.logical(x) && is.null(dim(x)) && !is.object(x) && length(x) == 1L && !is.na(x)
}

gmm_scalar_character <- function(x) {
  is.character(x) && is.null(dim(x)) && !is.object(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

gmm_scalar_numeric <- function(x, positive = FALSE, nonneg = FALSE) {
  is.numeric(x) && is.null(dim(x)) && !is.object(x) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && (!positive || x > 0) && (!nonneg || x >= 0)
}

gmm_int_bound <- function(x, name, allow_zero = FALSE) {
  min_val <- if (allow_zero) 0 else 1
  if (!gmm_scalar_numeric(x) || x != floor(x) || x < min_val || x > 2^53 - 1) {
    cdid_gmm_abort(
      paste0("`", name, "` must be an explicit finite integer bound in [", min_val, ", 2^53 - 1]."),
      c("cdid_gmm_bound_error", "cdid_gmm_control_error")
    )
  }
  as.double(x)
}

gmm_safe_product <- function(a, b, max_entries = NULL, label = "allocation") {
  if (!gmm_scalar_numeric(a, nonneg = TRUE) || a != floor(a) ||
      !gmm_scalar_numeric(b, nonneg = TRUE) || b != floor(b)) {
    cdid_gmm_abort(
      paste0("Invalid dimensions for ", label, "."),
      "cdid_gmm_bound_error"
    )
  }
  da <- as.double(a)
  db <- as.double(b)
  limit <- 2^53 - 1
  if (da > 0 && db > limit / da) {
    cdid_gmm_abort(
      paste0("Entry count for ", label, " overflows exact IEEE binary64 integer range."),
      "cdid_gmm_bound_error",
      diagnostics = list(label = label, dim_a = da, dim_b = db)
    )
  }
  prod_val <- da * db
  if (!is.finite(prod_val) || prod_val > limit) cdid_gmm_abort("Allocation product exceeds exact integer range.", "cdid_gmm_bound_error")
  if (!is.null(max_entries) && prod_val > max_entries) {
    cdid_gmm_abort(
      paste0("Allocation for ", label, " (", prod_val, " entries) exceeds explicit `max_matrix_entries` (", max_entries, ")."),
      "cdid_gmm_bound_error",
      diagnostics = list(label = label, requested_entries = prod_val, max_matrix_entries = max_entries)
    )
  }
  prod_val
}

gmm_gamma <- function(n) {
  # Upper endpoint for the real gamma_n; not used to certify BLAS products.
  nd <- gmm_int_bound(n, "gamma dimension", allow_zero = TRUE)
  if (nd == 0) return(0)
  nu <- gmm_iv_mul(gmm_iv(nd), gmm_iv(2^-53))
  if (nu[2L] >= 0.5) cdid_gmm_abort("Dimension exceeds gamma bound domain.", "cdid_gmm_accuracy_error")
  gmm_iv_div(nu, gmm_iv_sub(gmm_iv(1), nu))[2L]
}

# ---------------------------------------------------------------------------
# Exact IEEE 754 binary64 dyadic summation (no tolerance, no extra dependency)
# ---------------------------------------------------------------------------
# Every finite IEEE 754 binary64 value is an exact integer multiple of 2^-1074
# whose bits lie in positions [-1074, 1023] (2098 bit positions = 70 bins of
# width 30 bits). Accumulating the exact 53-bit significands into 71 radix-2^30
# bins with canonical nonnegative remainders decides sum(x) == 0 and sign(sum(x))
# with 100% mathematical exactness across the entire binary64 exponent range.

exact_binary64_dyadic_sum <- function(x) {
  if (!is.atomic(x) || !is.null(dim(x)) || is.object(x) || !is.numeric(x) ||
      anyNA(x) || any(!is.finite(x))) {
    cdid_gmm_abort(
      "Target coefficients must be finite atomic numeric vectors without NA/NaN/Inf.",
      "cdid_gmm_target_error"
    )
  }
  n <- length(x)
  if (n == 0L || all(x == 0)) {
    return(list(is_zero = TRUE, sign = 0L, fast_path = TRUE, approx_double = 0))
  }
  xd <- as.double(x)
  # Fast path for sparse {0, -1, +1} or small exact integers whose absolute sum
  # cannot exceed 2^53 - 1.
  if (all(xd == 0 | xd == 1 | xd == -1)) {
    s <- sum(xd)
    return(list(is_zero = (s == 0), sign = as.integer(sign(s)), fast_path = TRUE, approx_double = s))
  }
  if (n <= 1048576L && all(xd == floor(xd) & abs(xd) <= 1073741824)) {
    s <- sum(xd)
    return(list(is_zero = (s == 0), sign = as.integer(sign(s)), fast_path = TRUE, approx_double = s))
  }

  nz <- xd[xd != 0]
  m_nz <- length(nz)
  base_radix <- 1073741824 # 2^30
  n_bins <- 71L
  bins <- numeric(n_bins)
  top_carry <- 0

  normalize_bins <- function() {
    carry <- 0
    for (b in seq_len(n_bins)) {
      val <- bins[b] + carry
      rem <- val %% base_radix
      if (rem < 0) rem <- rem + base_radix
      carry <- (val - rem) / base_radix
      bins[b] <<- rem
    }
    top_carry <<- top_carry + carry
  }

  pow2_26 <- 2^(0:25)
  pow2_11 <- 2L^(0:10)
  chunk_len <- 1024L
  n_chunks <- as.integer(ceiling(m_nz / chunk_len))

  for (ch in seq_len(n_chunks)) {
    idx <- seq.int((ch - 1L) * chunk_len + 1L, min(ch * chunk_len, m_nz))
    sub_x <- nz[idx]
    bits <- matrix(as.integer(numToBits(sub_x)), nrow = 64L)
    low26 <- as.double(colSums(bits[1:26, , drop = FALSE] * pow2_26))
    high26 <- as.double(colSums(bits[27:52, , drop = FALSE] * pow2_26))
    exp_biased <- as.integer(colSums(bits[53:63, , drop = FALSE] * pow2_11))
    sgn <- ifelse(bits[64L, ] == 1L, -1, 1)

    high27 <- high26 + ifelse(exp_biased > 0L, 67108864, 0)
    shift_lo <- ifelse(exp_biased > 0L, exp_biased - 1L, 0L)
    shift_hi <- shift_lo + 26L

    for (part in list(list(v = low26, sh = shift_lo), list(v = high27, sh = shift_hi))) {
      v <- part$v
      sh <- part$sh
      b_idx <- (sh %/% 30L) + 1L
      rem_bit <- sh %% 30L
      split_mod <- 2^(30L - rem_bit)
      lo_val <- (v %% split_mod) * (2^rem_bit)
      hi_val <- floor(v / split_mod)
      for (i in seq_along(v)) {
        if (lo_val[i] != 0) bins[b_idx[i]] <- bins[b_idx[i]] + sgn[i] * lo_val[i]
        if (hi_val[i] != 0) bins[b_idx[i] + 1L] <- bins[b_idx[i] + 1L] + sgn[i] * hi_val[i]
      }
    }
    normalize_bins()
  }

  is_zero <- (top_carry == 0) && all(bins == 0)
  if (is_zero) {
    return(list(is_zero = TRUE, sign = 0L, fast_path = FALSE, approx_double = 0, bins = bins, top_carry = top_carry))
  }
  out_sign <- if (top_carry < 0) -1L else 1L

  # Reconstruct an approximate double from the highest active bins for diagnostics
  mag_bins <- bins
  if (out_sign < 0) {
    # Convert two's-complement radix-2^30 representation to positive magnitude
    borrow <- 0
    for (b in seq_len(n_bins)) {
      val <- -mag_bins[b] - borrow
      rem <- val %% base_radix
      if (rem < 0) rem <- rem + base_radix
      borrow <- if (val < 0) 1 else 0
      mag_bins[b] <- rem
    }
  }
  nz_bins <- which(mag_bins > 0)
  top_b <- max(nz_bins)
  lead <- mag_bins[top_b]
  if (top_b >= 2L) lead <- lead + mag_bins[top_b - 1L] / base_radix
  if (top_b >= 3L) lead <- lead + mag_bins[top_b - 2L] / (base_radix^2)
  bit_exp <- (top_b - 1L) * 30L - 1074L
  approx_double <- out_sign * ldexp_safe(lead, bit_exp)

  list(
    is_zero = FALSE,
    sign = out_sign,
    fast_path = FALSE,
    approx_double = approx_double,
    bins = bins,
    top_carry = top_carry
  )
}

ldexp_safe <- function(m, e) {
  # Multiply m by 2^e in steps to avoid intermediate 2^e overflow/underflow
  res <- as.double(m)
  rem <- as.integer(e)
  while (rem > 0L) {
    step <- min(rem, 512L)
    res <- res * (2^step)
    rem <- rem - step
  }
  while (rem < 0L) {
    step <- max(rem, -512L)
    res <- res * (2^step)
    rem <- rem - step
  }
  res
}

# ---------------------------------------------------------------------------
# Explicit caller control constructors and validators (NO internal defaults)
# ---------------------------------------------------------------------------

gmm_bounds_control <- function(
  max_matrix_entries,
  max_batch_rows,
  max_targets,
  max_open_clusters,
  allow_materialize_if,
  max_materialized_if_entries
) {
  if (missing(max_matrix_entries) || missing(max_batch_rows) ||
      missing(max_targets) || missing(max_open_clusters) ||
      missing(allow_materialize_if) || missing(max_materialized_if_entries)) {
    cdid_gmm_abort(
      "All fields of `gmm_bounds_control` must be supplied explicitly by the caller.",
      c("cdid_gmm_bound_error", "cdid_gmm_control_error")
    )
  }
  out <- list(
    max_matrix_entries = gmm_int_bound(max_matrix_entries, "max_matrix_entries"),
    max_batch_rows = gmm_int_bound(max_batch_rows, "max_batch_rows"),
    max_targets = gmm_int_bound(max_targets, "max_targets"),
    max_open_clusters = gmm_int_bound(max_open_clusters, "max_open_clusters"),
    allow_materialize_if = if (gmm_scalar_logical(allow_materialize_if)) allow_materialize_if else
      cdid_gmm_abort("`allow_materialize_if` must be a non-missing scalar logical.", "cdid_gmm_control_error"),
    max_materialized_if_entries = gmm_int_bound(max_materialized_if_entries, "max_materialized_if_entries", allow_zero = TRUE)
  )
  class(out) <- "cdid_gmm_bounds_control"
  out
}

validate_gmm_bounds_control <- function(bounds) {
  if (missing(bounds) || is.null(bounds) || !is.list(bounds) || anyDuplicated(names(bounds))) {
    cdid_gmm_abort("Explicit `bounds` control list is required.", c("cdid_gmm_bound_error", "cdid_gmm_control_error"))
  }
  req <- c("max_matrix_entries", "max_batch_rows", "max_targets",
           "max_open_clusters", "allow_materialize_if", "max_materialized_if_entries")
  if (!all(req %in% names(bounds))) {
    cdid_gmm_abort(
      paste0("Missing required `bounds` control fields: ", paste(setdiff(req, names(bounds)), collapse = ", "), "."),
      c("cdid_gmm_bound_error", "cdid_gmm_control_error")
    )
  }
  gmm_bounds_control(
    max_matrix_entries = bounds$max_matrix_entries,
    max_batch_rows = bounds$max_batch_rows,
    max_targets = bounds$max_targets,
    max_open_clusters = bounds$max_open_clusters,
    allow_materialize_if = bounds$allow_materialize_if,
    max_materialized_if_entries = bounds$max_materialized_if_entries
  )
}

gmm_accuracy_control <- function(
  method,
  target_error_unit,
  max_target_error_abs,
  max_target_error_se,
  max_orthogonality_defect,
  min_certified_lambda_rel,
  max_factor_residual_rel,
  min_rcond_R,
  min_rcond_U,
  min_sigma_rel_A,
  max_kkt_residual_rel,
  max_loading_residual_rel,
  require_matching_scale,
  max_loading_error_abs = NULL,
  covariance_strategy = NULL
) {
  if (missing(method) || missing(target_error_unit) ||
      missing(max_target_error_abs) || missing(max_target_error_se) ||
      missing(max_orthogonality_defect) || missing(min_certified_lambda_rel) ||
      missing(max_factor_residual_rel) || missing(min_rcond_R) ||
      missing(min_rcond_U) || missing(min_sigma_rel_A) ||
      missing(max_kkt_residual_rel) || missing(max_loading_residual_rel) ||
      missing(require_matching_scale)) {
    cdid_gmm_abort(
      "All fields of `gmm_accuracy_control` must be supplied explicitly by the caller (no internal defaults).",
      "cdid_gmm_control_error"
    )
  }
  out <- list(
    method = method,
    target_error_unit = target_error_unit,
    max_target_error_abs = max_target_error_abs,
    max_target_error_se = max_target_error_se,
    max_orthogonality_defect = max_orthogonality_defect,
    min_certified_lambda_rel = min_certified_lambda_rel,
    max_factor_residual_rel = max_factor_residual_rel,
    min_rcond_R = min_rcond_R,
    min_rcond_U = min_rcond_U,
    min_sigma_rel_A = min_sigma_rel_A,
    max_kkt_residual_rel = max_kkt_residual_rel,
    max_loading_residual_rel = max_loading_residual_rel,
    require_matching_scale = require_matching_scale,
    max_loading_error_abs = max_loading_error_abs,
    covariance_strategy = covariance_strategy
  )
  class(out) <- "cdid_gmm_accuracy_control"
  out
}

validate_gmm_accuracy_control <- function(accuracy, weighting) {
  if (missing(accuracy) || is.null(accuracy) || !is.list(accuracy) || anyDuplicated(names(accuracy))) {
    cdid_gmm_abort(
      "Caller-explicit `accuracy` specification is required; no internal numerical defaults exist.",
      "cdid_gmm_control_error"
    )
  }
  req <- c(
    "method", "target_error_unit", "max_target_error_abs", "max_target_error_se",
    "max_orthogonality_defect", "min_certified_lambda_rel", "max_factor_residual_rel",
    "min_rcond_R", "min_rcond_U", "min_sigma_rel_A",
    "max_kkt_residual_rel", "max_loading_residual_rel", "require_matching_scale"
  )
  if (!all(req %in% names(accuracy))) {
    cdid_gmm_abort(
      paste0("Missing required `accuracy` fields: ", paste(setdiff(req, names(accuracy)), collapse = ", "), "."),
      "cdid_gmm_control_error"
    )
  }
  if (!gmm_scalar_character(accuracy$method)) {
    cdid_gmm_abort("`accuracy$method` must be a non-empty character scalar.", "cdid_gmm_control_error")
  }
  allowed_methods <- if (weighting == "optimal") {
    c("certified_eigen_chol_kkt", "conventional_chol_qr")
  } else {
    c("certified_eigen_chol_kkt", "certified_qr_kkt", "conventional_chol_qr")
  }
  if (!accuracy$method %in% allowed_methods) {
    cdid_gmm_abort(
      paste0("Unsupported or uncertified validation method `", accuracy$method, "`. ",
             "Cholesky success, rcond alone, metadata flags, or ordinary computed indicators cannot certify PD/accuracy."),
      c("cdid_gmm_accuracy_error", "cdid_gmm_control_error"),
      diagnostics = list(method = accuracy$method, weighting = weighting, allowed_methods = allowed_methods)
    )
  }
  if (!gmm_scalar_character(accuracy$target_error_unit) ||
      !accuracy$target_error_unit %in% c("abs", "model_se", "both")) {
    cdid_gmm_abort("`accuracy$target_error_unit` must be explicitly 'abs', 'model_se', or 'both'.", "cdid_gmm_control_error")
  }
  if (weighting == "identity" && accuracy$target_error_unit %in% c("model_se", "both")) {
    cdid_gmm_abort(
      "Identity weighting has no model-variance SE interpretation for target-error budgeting; use `target_error_unit = 'abs'`.",
      c("cdid_gmm_scale_error", "cdid_gmm_control_error")
    )
  }
  if (accuracy$target_error_unit %in% c("abs", "both") &&
      !gmm_scalar_numeric(accuracy$max_target_error_abs, positive = TRUE)) {
    cdid_gmm_abort("`max_target_error_abs` must be a finite positive numeric scalar.", "cdid_gmm_control_error")
  }
  if (accuracy$target_error_unit == "model_se" &&
      !(is.na(accuracy$max_target_error_abs) || gmm_scalar_numeric(accuracy$max_target_error_abs, positive = TRUE))) {
    cdid_gmm_abort("`max_target_error_abs` must be finite positive or NA when `target_error_unit == 'model_se'`.", "cdid_gmm_control_error")
  }
  if (accuracy$target_error_unit %in% c("model_se", "both") &&
      !gmm_scalar_numeric(accuracy$max_target_error_se, positive = TRUE)) {
    cdid_gmm_abort("`max_target_error_se` must be a finite positive numeric scalar.", "cdid_gmm_control_error")
  }
  if (accuracy$target_error_unit == "abs" &&
      !(is.na(accuracy$max_target_error_se) || gmm_scalar_numeric(accuracy$max_target_error_se, positive = TRUE))) {
    cdid_gmm_abort("`max_target_error_se` must be finite positive or NA when `target_error_unit == 'abs'`.", "cdid_gmm_control_error")
  }
  if (!gmm_scalar_numeric(accuracy$max_orthogonality_defect, positive = TRUE) ||
      accuracy$max_orthogonality_defect >= 0.5) {
    cdid_gmm_abort("`max_orthogonality_defect` must be a positive scalar in (0, 0.5).", "cdid_gmm_control_error")
  }
  for (nm in c("max_factor_residual_rel", "min_rcond_U", "min_sigma_rel_A",
               "max_kkt_residual_rel", "max_loading_residual_rel")) {
    if (!gmm_scalar_numeric(accuracy[[nm]], positive = TRUE) || accuracy[[nm]] >= 1) {
      cdid_gmm_abort(paste0("`accuracy$", nm, "` must be a finite positive scalar in (0, 1)."), "cdid_gmm_control_error")
    }
  }
  if (weighting == "optimal") {
    for (nm in c("min_certified_lambda_rel", "min_rcond_R")) {
      if (!gmm_scalar_numeric(accuracy[[nm]], positive = TRUE) || accuracy[[nm]] >= 1) {
        cdid_gmm_abort(paste0("`accuracy$", nm, "` must be a finite positive scalar in (0, 1) for optimal weighting."), "cdid_gmm_control_error")
      }
    }
  } else {
    for (nm in c("min_certified_lambda_rel", "min_rcond_R")) {
      val <- accuracy[[nm]]
      if (!(is.na(val) || (gmm_scalar_numeric(val, positive = TRUE) && val < 1))) {
        cdid_gmm_abort(paste0("`accuracy$", nm, "` must be finite positive in (0, 1) or NA for identity weighting."), "cdid_gmm_control_error")
      }
    }
  }
  if (!gmm_scalar_logical(accuracy$require_matching_scale)) {
    cdid_gmm_abort("`accuracy$require_matching_scale` must be a non-missing scalar logical.", "cdid_gmm_control_error")
  }
  # Optional private NULL on `max_loading_error_abs` is an explicit absence state
  # ("budget not supplied"), never a numerical default or repurposed relative ratio.
  if (!is.null(accuracy[["max_loading_error_abs"]]) &&
      !gmm_scalar_numeric(accuracy[["max_loading_error_abs"]], positive = TRUE)) {
    cdid_gmm_abort(
      "`accuracy$max_loading_error_abs` must be NULL (absence state: not supplied) or a finite positive numeric scalar.",
      "cdid_gmm_control_error"
    )
  }
  accuracy$max_loading_error_abs <- accuracy[["max_loading_error_abs"]]
  accuracy$covariance_strategy <- validate_gmm_covariance_strategy(accuracy[["covariance_strategy"]], weighting)
  accuracy
}

gmm_covariance_provenance <- function(
  row_level,
  n_units,
  n_clusters,
  n_positive_contributing_rows,
  centered,
  weighted,
  weight_type,
  normalization_factor,
  scale_id,
  nuisance_fit_id,
  reduction = list(mode = "none")
) {
  if (missing(row_level) || missing(n_units) || missing(n_clusters) ||
      missing(n_positive_contributing_rows) || missing(centered) ||
      missing(weighted) || missing(weight_type) ||
      missing(normalization_factor) || missing(scale_id) ||
      missing(nuisance_fit_id)) {
    cdid_gmm_abort(
      "All fields of `gmm_covariance_provenance` must be supplied explicitly by the caller.",
      "cdid_gmm_control_error"
    )
  }
  out <- list(
    row_level = row_level,
    n_units = n_units,
    n_clusters = n_clusters,
    n_positive_contributing_rows = n_positive_contributing_rows,
    centered = centered,
    weighted = weighted,
    weight_type = weight_type,
    normalization_factor = normalization_factor,
    scale_id = scale_id,
    nuisance_fit_id = nuisance_fit_id,
    reduction = reduction
  )
  class(out) <- "cdid_gmm_covariance_provenance"
  out
}

validate_gmm_covariance_provenance <- function(provenance, weighting) {
  if (missing(provenance) || is.null(provenance) || !is.list(provenance) || anyDuplicated(names(provenance))) {
    cdid_gmm_abort("Explicit covariance `provenance` is required.", "cdid_gmm_control_error")
  }
  if (isTRUE(provenance$certified_pd) || isTRUE(provenance$is_pd) || isTRUE(provenance$trust_metadata_pd)) {
    cdid_gmm_abort(
      "Untrusted metadata PD flags (`is_pd`/`certified_pd`) cannot certify positive definiteness.",
      c("cdid_gmm_covariance_error", "cdid_gmm_control_error")
    )
  }
  req <- c("row_level", "n_units", "n_clusters", "n_positive_contributing_rows",
           "centered", "weighted", "weight_type", "normalization_factor",
           "scale_id", "nuisance_fit_id", "reduction")
  if (!all(req %in% names(provenance))) {
    cdid_gmm_abort(
      paste0("Missing required `provenance` fields: ", paste(setdiff(req, names(provenance)), collapse = ", "), "."),
      "cdid_gmm_control_error"
    )
  }
  if (!gmm_scalar_character(provenance$row_level) || !provenance$row_level %in% c("unit", "cluster")) {
    cdid_gmm_abort("`provenance$row_level` must be 'unit' or 'cluster'.", "cdid_gmm_control_error")
  }
  n_units <- gmm_int_bound(provenance$n_units, "provenance$n_units")
  n_pos <- gmm_int_bound(provenance$n_positive_contributing_rows, "provenance$n_positive_contributing_rows")
  if (provenance$row_level == "cluster") {
    n_clusters <- gmm_int_bound(provenance$n_clusters, "provenance$n_clusters")
    if (n_pos > n_clusters || n_clusters > n_units) {
      cdid_gmm_abort("Invalid cluster/contributing row counts in `provenance`.", "cdid_gmm_control_error")
    }
  } else {
    if (!is.na(provenance$n_clusters)) {
      n_clusters <- gmm_int_bound(provenance$n_clusters, "provenance$n_clusters")
      if (n_clusters > n_units) cdid_gmm_abort("`n_clusters` cannot exceed `n_units`.", "cdid_gmm_control_error")
    }
    if (n_pos > n_units) {
      cdid_gmm_abort("`n_positive_contributing_rows` cannot exceed `n_units`.", "cdid_gmm_control_error")
    }
  }
  if (!gmm_scalar_logical(provenance$centered) || !gmm_scalar_logical(provenance$weighted)) {
    cdid_gmm_abort("`provenance$centered` and `provenance$weighted` must be scalar logicals.", "cdid_gmm_control_error")
  }
  if (!gmm_scalar_character(provenance$weight_type) ||
      !gmm_scalar_character(provenance$scale_id) ||
      !gmm_scalar_character(provenance$nuisance_fit_id)) {
    cdid_gmm_abort("`weight_type`, `scale_id`, and `nuisance_fit_id` must be non-empty character scalars.", "cdid_gmm_control_error")
  }
  if (!gmm_scalar_numeric(provenance$normalization_factor, positive = TRUE)) {
    cdid_gmm_abort("`provenance$normalization_factor` must be a finite positive scalar.", "cdid_gmm_control_error")
  }
  validate_gmm_reduction_spec(provenance$reduction, weighting)
  provenance
}

validate_gmm_reduction_spec <- function(reduction, weighting) {
  if (is.null(reduction) || !is.list(reduction) || !gmm_scalar_character(reduction$mode)) {
    cdid_gmm_abort("`provenance$reduction` must be a list with an explicit `mode` field.", "cdid_gmm_reduction_error")
  }
  mode <- reduction$mode
  if (mode == "none") {
    if (weighting == "identity" && isTRUE(reduction$reduction_active)) {
      cdid_gmm_abort(
        "Structural reduction is forbidden under identity weighting because dropping rows changes the identity criterion.",
        "cdid_gmm_reduction_error",
        diagnostics = list(mode = mode, weighting = weighting)
      )
    }
    return(invisible(reduction))
  }
  if (weighting == "identity") {
    cdid_gmm_abort(
      "Structural reduction is forbidden under identity weighting because dropping rows changes the identity criterion.",
      "cdid_gmm_reduction_error",
      diagnostics = list(mode = mode, weighting = weighting)
    )
  }
  if (mode == "sparse_class_certificate") {
    if (is.null(reduction$class_certificate) && is.null(reduction$E) && is.null(reduction$G)) {
      cdid_gmm_abort(
        "`reduction$mode == 'sparse_class_certificate'` requires `class_certificate`, sparse `E`, or `G`.",
        "cdid_gmm_reduction_error"
      )
    }
    return(invisible(reduction))
  }
  if (mode %in% c("auto_detect_ipw", "numeric_agreement", "class_label_only", "sample_rank")) {
    cdid_gmm_abort(
      paste0("Reduction mode `", mode, "` is refused: the IPW kernel does not yet expose a shared nuisance-object API, ",
             "and numeric fit agreement, sample rank, or bare class labels never certify structural reduction."),
      "cdid_gmm_reduction_error",
      diagnostics = list(mode = mode)
    )
  }
  if (mode != "pre_reduced") {
    cdid_gmm_abort(
      paste0("Unknown reduction mode `", mode, "`. Only `none`, `sparse_class_certificate`, or already-proved `pre_reduced` provenance is accepted."),
      "cdid_gmm_reduction_error"
    )
  }
  cdid_gmm_abort("Pre-reduced membership cannot be verified from string IDs/certified_by. Supply the unreduced problem or await the actual shared-object/mask/map verifier.", "cdid_gmm_reduction_error")
}

# ---------------------------------------------------------------------------
# Spanning forest across all connected components of a certified class
# (Private algebra helper for proved structural classes; never auto-certifies)
# ---------------------------------------------------------------------------

# Conditional graph algebra only: class membership must ALREADY be proved externally.
# Neither string equality nor TRUE metadata verifies an IPW class.
build_certified_class_forest <- function(class_catalog, weighting, max_matrix_entries) {
  if (missing(max_matrix_entries)) cdid_gmm_abort("Conditional forest requires an explicit allocation bound.", "cdid_gmm_bound_error")
  forest_bound <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
  if (missing(weighting) || !identical(weighting, "optimal")) {
    cdid_gmm_abort(
      "Certified structural reduction requires `weighting = 'optimal'`; never reduce identity weighting.",
      "cdid_gmm_reduction_error"
    )
  }
  if (!is.data.frame(class_catalog) || nrow(class_catalog) < 1L) {
    cdid_gmm_abort("`class_catalog` must be a non-empty data.frame.", "cdid_gmm_reduction_error")
  }
  req <- c("s", "t", "context_id", "g", "operator_id", "vertex_map_id",
           "pair_s_id", "nuisance_object_id")
  if (!all(req %in% names(class_catalog))) {
    cdid_gmm_abort("Missing structural key columns in `class_catalog`.", "cdid_gmm_reduction_error")
  }
  m <- nrow(class_catalog)
  m_sq_forest <- gmm_safe_product(m, m, forest_bound, "conditional forest link-pair storage")
  m_lin_forest <- gmm_safe_product(10, m, forest_bound, "conditional forest linear link storage")
  gmm_safe_product(2, m_sq_forest + m_lin_forest, forest_bound, "conditional forest E/path/graph storage")
  gmm_safe_product(4, m, forest_bound, "conditional forest graph/key buffers")
  if (!is.numeric(class_catalog$g) || anyNA(class_catalog$g) || any(!is.finite(class_catalog$g)))
    cdid_gmm_abort("Invalid conditional cohort values.", "cdid_gmm_reduction_error")
  if (!all(is.finite(class_catalog$s)) || !all(is.finite(class_catalog$t)) ||
      any(class_catalog$s >= class_catalog$t)) {
    cdid_gmm_abort("Class links require finite endpoints with s < t.", "cdid_gmm_reduction_error")
  }
  for (nm in c("context_id", "operator_id", "vertex_map_id", "pair_s_id", "nuisance_object_id")) {
    col <- class_catalog[[nm]]
    if (!is.character(col) || anyNA(col) || any(!nzchar(col))) {
      cdid_gmm_abort(paste0("Invalid structural key column `", nm, "`."), "cdid_gmm_reduction_error")
    }
  }

  # Collision-free equality of complete tuples; delimiter concatenation is not identity.
  keys <- c("context_id", "g", "operator_id", "vertex_map_id", "pair_s_id", "nuisance_object_id")
  class_key <- integer(m)
  next_key <- 0L
  for (i in seq_len(m)) {
    if (class_key[i] != 0L) next
    next_key <- next_key + 1L
    same <- rep(TRUE, m)
    for (nm in keys) same <- same & class_catalog[[nm]] == class_catalog[[nm]][i]
    class_key[same] <- next_key
  }
  u_keys <- unique(class_key)
  retained <- integer(0)
  e_paths <- vector("list", m)

  for (ck in u_keys) {
    ls <- which(class_key == ck)
    if (anyDuplicated(data.frame(s=class_catalog$s[ls],t=class_catalog$t[ls])))
      cdid_gmm_abort("Duplicate conditional forest links.", "cdid_gmm_reduction_error")
    ls <- ls[order(class_catalog$s[ls], class_catalog$t[ls], ls)]
    dates <- sort(unique(c(class_catalog$s[ls], class_catalog$t[ls])))
    seen <- numeric(0)
    par <- list()
    tree_edges <- integer(0)

    # Restart BFS at every unvisited smallest date so disconnected date components
    # within the same class are all traversed.
    for (root in dates) {
      if (root %in% seen) next
      seen <- c(seen, root)
      queue <- root
      while (length(queue) > 0L) {
        u <- queue[1L]
        queue <- queue[-1L]
        for (l in ls) {
          s_l <- class_catalog$s[l]
          t_l <- class_catalog$t[l]
          if (s_l == u && !(t_l %in% seen)) {
            seen <- c(seen, t_l)
            par[[as.character(match(t_l, dates))]] <- list(parent = u, edge = l, sign = 1L)
            tree_edges <- c(tree_edges, l)
            queue <- c(queue, t_l)
          } else if (t_l == u && !(s_l %in% seen)) {
            seen <- c(seen, s_l)
            par[[as.character(match(s_l, dates))]] <- list(parent = u, edge = l, sign = -1L)
            tree_edges <- c(tree_edges, l)
            queue <- c(queue, s_l)
          }
        }
      }
    }
    tree_edges <- sort(unique(tree_edges))
    retained <- c(retained, tree_edges)
    for (l in ls) {
      if (l %in% tree_edges) {
        vec <- integer(m)
        vec[l] <- 1L
        e_paths[[l]] <- vec
      } else {
        vec <- integer(m)
        for (endpoint in list(list(d = class_catalog$t[l], sgn = 1L), list(d = class_catalog$s[l], sgn = -1L))) {
          cur <- endpoint$d
          repeat {
            info <- par[[as.character(match(cur, dates))]]
            if (is.null(info)) break
            vec[info$edge] <- vec[info$edge] + endpoint$sgn * info$sign
            cur <- info$parent
          }
        }
        e_paths[[l]] <- vec
      }
    }
  }

  retained <- sort(unique(retained))
  E <- matrix(0L, nrow = m, ncol = length(retained))
  for (l in seq_len(m)) {
    E[l, ] <- e_paths[[l]][retained]
  }
  list(retained_links = retained, E = E, class_key = class_key,
       structural_proof_verified = FALSE,
       scope = "Conditional forest algebra on externally proved classes; labels and TRUE metadata supply no IPW certificate.")
}

# ---------------------------------------------------------------------------
# Deterministic unanchored component gauge (without dropping rows)
# ---------------------------------------------------------------------------

build_gmm_gauge <- function(design, gauge_rule, custom_roots = NULL) {
  if (!inherits(design, "cdid_chained_design")) {
    cdid_gmm_abort("`design` must be a `cdid_chained_design` object.", "cdid_gmm_design_error")
  }
  if (missing(gauge_rule) || !gmm_scalar_character(gauge_rule) ||
      !gauge_rule %in% c("lowest", "highest", "custom")) {
    cdid_gmm_abort(
      "`gauge_rule` must be supplied explicitly as 'lowest', 'highest', or 'custom'.",
      "cdid_gmm_control_error"
    )
  }
  inc <- design$incidence
  params <- design$parameters
  rk <- design$rank_diagnostics
  if (!gmm_scalar_numeric(inc$n_links, positive=TRUE) || !gmm_scalar_numeric(inc$n_parameters,positive=TRUE) ||
      inc$n_links > .Machine$integer.max || inc$n_parameters > .Machine$integer.max ||
      inc$n_links != floor(inc$n_links) || inc$n_parameters != floor(inc$n_parameters))
    cdid_gmm_abort("Invalid R incidence dimensions.", "cdid_gmm_design_error")
  m <- as.integer(inc$n_links)
  p <- as.integer(inc$n_parameters)
  if (m < 1L || p < 1L || length(inc$start) != m || length(inc$end) != m || nrow(params) != p) {
    cdid_gmm_abort("Malformed `design` incidence or parameter table.", "cdid_gmm_design_error")
  }
  if (anyNA(inc$start) || anyNA(inc$end) || any(inc$start != floor(inc$start)) || any(inc$end != floor(inc$end))) {
    cdid_gmm_abort("Malformed incidence indices/IDs.", "cdid_gmm_design_error")
  }
  zero_inc <- inc$start == 0L & inc$end == 0L
  if (any(zero_inc)) {
    L <- design$links
    V <- design$vertices
    cspec <- design$calendar_spec
    if (!is.data.frame(L) || !is.data.frame(V) || nrow(L) != m ||
        !all(c("start_vertex_index", "end_vertex_index", "g", "s", "t", "ref", "functional_key") %in% names(L)) ||
        !all(c("vertex_index", "param_index", "g", "time", "role", "functional_key", "anchored_value") %in% names(V)) ||
        is.null(cspec) || !identical(cspec$pt_assumption, "all") ||
        !gmm_scalar_numeric(cspec$pt_start, positive = FALSE)) {
      cdid_gmm_abort("`design` contains invalid or self-loop incidence endpoints.", "cdid_gmm_design_error")
    }
    zi <- which(zero_inc)
    sv <- L$start_vertex_index[zi]
    ev <- L$end_vertex_index[zi]
    if (anyNA(sv) || anyNA(ev) || any(sv < 1L | sv > nrow(V) | ev < 1L | ev > nrow(V) | sv == ev)) {
      cdid_gmm_abort("`design` contains invalid or self-loop incidence endpoints.", "cdid_gmm_design_error")
    }
    is_valid_anchor <- function(vi, l_idx) {
      is.na(V$param_index[vi]) & !is.na(V$anchored_value[vi]) & V$anchored_value[vi] == 0 &
        ((V$role[vi] == "reference" & is.na(V$functional_key[vi]) & V$time[vi] == L$ref[l_idx]) |
         (V$role[vi] == "pre" & !is.na(V$functional_key[vi]) &
          V$functional_key[vi] == L$functional_key[l_idx] &
          V$time[vi] >= cspec$pt_start & V$time[vi] < L$ref[l_idx]))
    }
    if (any(!is_valid_anchor(sv, zi) | !is_valid_anchor(ev, zi) |
            V$g[sv] != L$g[zi] | V$g[ev] != L$g[zi] |
            V$time[sv] != L$s[zi] | V$time[ev] != L$t[zi] |
            V$time[sv] >= V$time[ev])) {
      cdid_gmm_abort("`design` contains invalid or self-loop incidence endpoints.", "cdid_gmm_design_error")
    }
  }
  if (any(inc$start < 0L | inc$start > p | inc$end < 0L | inc$end > p | (!zero_inc & inc$start == inc$end))) {
    cdid_gmm_abort("`design` contains invalid or self-loop incidence endpoints.", "cdid_gmm_design_error")
  }

  if (!identical(as.integer(params$param_index), seq_len(p)) || anyNA(params$component_id) ||
      anyDuplicated(inc$link_ids) || anyDuplicated(inc$param_ids) || length(inc$link_ids)!=m || length(inc$param_ids)!=p ||
      anyNA(inc$link_ids) || anyNA(inc$param_ids) || any(!nzchar(inc$link_ids)) || any(!nzchar(inc$param_ids)))
    cdid_gmm_abort("Malformed incidence indices/IDs.", "cdid_gmm_design_error")
  # Independently recompute combinatorial components, treating ground edges as anchors.
  adj <- vector("list", p); touched <- grounded <- logical(p)
  for (ii in seq_len(m)) {
    a <- inc$start[ii]; b <- inc$end[ii]
    if (a == 0L && b == 0L) next
    if (a>0) touched[a] <- TRUE
    if (b>0) touched[b] <- TRUE
    if (a==0) grounded[b] <- TRUE else if (b==0) grounded[a] <- TRUE else {
      adj[[a]] <- c(adj[[a]],b); adj[[b]] <- c(adj[[b]],a)
    }
  }
  comp_actual <- integer(p); level_actual <- logical(p); nc <- 0L
  for (v in seq_len(p)) {
    if (comp_actual[v] != 0L) next
    nc <- nc+1L; queue <- integer(p); queue[1L] <- v; head <- tail <- 1L; comp_actual[v] <- nc
    while (head <= tail) {
      u <- queue[head]; head <- head+1L
      for (w in adj[[u]]) if (comp_actual[w]==0L) { tail <- tail+1L; queue[tail] <- w; comp_actual[w] <- nc }
    }
    members <- queue[seq_len(tail)]
    level_actual[members] <- any(grounded[members])
    if (!level_actual[v]) {
      label <- params$component_id[v]
      if (!setequal(which(params$component_id==label),members) || !label %in% rk$unanchored_components)
        cdid_gmm_abort("Stored unanchored component partition disagrees with graph.", "cdid_gmm_design_error")
    }
  }
  actual_unanch <- unique(params$component_id[!level_actual])
  if (!setequal(actual_unanch,rk$unanchored_components) || anyNA(params$level_estimable) ||
      !identical(as.logical(params$level_estimable),level_actual) || anyNA(params$is_isolated) ||
      !identical(as.logical(params$is_isolated),!touched) || rk$rank != p-length(actual_unanch))
    cdid_gmm_abort("Stored level/rank/isolation claims disagree with graph.", "cdid_gmm_design_error")

  unanch <- as.integer(rk$unanchored_components)
  nu <- length(unanch)
  p_comp <- as.integer(params$component_id)

  root_params <- integer(nu)
  if (nu > 0L) {
    if (gauge_rule == "custom") {
      if (!is.numeric(custom_roots) || length(custom_roots) != nu || anyNA(custom_roots) ||
          any(custom_roots != floor(custom_roots))) {
        cdid_gmm_abort("`custom_roots` must supply one integer parameter index per unanchored component.", "cdid_gmm_control_error")
      }
      for (r in seq_len(nu)) {
        u <- unanch[r]
        members <- which(p_comp == u)
        cand <- as.integer(custom_roots[r])
        if (!cand %in% members) {
          cdid_gmm_abort("Each custom gauge root must belong to its unanchored component.", "cdid_gmm_control_error")
        }
        root_params[r] <- cand
      }
    } else {
      for (r in seq_len(nu)) {
        u <- unanch[r]
        members <- which(p_comp == u)
        if (!length(members)) {
          cdid_gmm_abort("Unanchored component has no parameter vertices.", "cdid_gmm_design_error")
        }
        root_params[r] <- if (gauge_rule == "lowest") min(members) else max(members)
      }
    }
  }

  keep_params <- setdiff(seq_len(p), root_params)
  k <- length(keep_params)
  if (k != as.integer(rk$rank)) {
    cdid_gmm_abort(
      "Gauge column count does not equal structural design rank p - nu.",
      "cdid_gmm_design_error",
      diagnostics = list(k = k, p = p, nu = nu, design_rank = rk$rank)
    )
  }

  # Map original parameter indices 1..p to gauge columns 1..k (0L for dropped roots and ground)
  param_to_col <- integer(p)
  param_to_col[keep_params] <- seq_len(k)
  start_col <- integer(m)
  pos_s <- inc$start > 0L
  start_col[pos_s] <- param_to_col[inc$start[pos_s]]
  end_col <- integer(m)
  pos_e <- inc$end > 0L
  end_col[pos_e] <- param_to_col[inc$end[pos_e]]

  # Proposition 1(a): no originally active row of W_G may collapse to zero or self-loop, and no row is dropped.
  if (any(!zero_inc & ((start_col == 0L & end_col == 0L) | (start_col == end_col)))) {
    cdid_gmm_abort(
      "Gauge incidence matrix W_G has a zero row; check grounded graph components.",
      "cdid_gmm_design_error"
    )
  }

  structure(
    list(
      n_links = m,
      n_parameters = p,
      rank_k = k,
      nullity_nu = nu,
      gauge_rule = gauge_rule,
      unanchored_components = unanch,
      root_params = root_params,
      keep_params = keep_params,
      param_to_col = param_to_col,
      start_col = start_col,
      end_col = end_col,
      param_component_id = p_comp,
      level_estimable = as.logical(params$level_estimable),
      is_isolated = as.logical(params$is_isolated),
      link_ids = inc$link_ids,
      param_ids = inc$param_ids,
      keep_param_ids = inc$param_ids[keep_params]
    ),
    class = "cdid_gmm_gauge"
  )
}

gmm_gauge_multiply <- function(gauge, x, transpose = FALSE) {
  if (!inherits(gauge, "cdid_gmm_gauge")) {
    cdid_gmm_abort("Expected `cdid_gmm_gauge` object.", "cdid_gmm_design_error")
  }
  n_in <- if (transpose) gauge$n_links else gauge$rank_k
  vec_in <- is.null(dim(x))
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
      (vec_in && (length(x) != n_in || is.object(x))) ||
      (!vec_in && (!is.matrix(x) || nrow(x) != n_in))) {
    cdid_gmm_abort("Invalid input to `gmm_gauge_multiply`.", "cdid_gmm_design_error")
  }
  if (vec_in) x <- matrix(as.double(x), ncol = 1L)
  if (!transpose) {
    grounded <- rbind(rep(0, ncol(x)), x)
    out <- grounded[gauge$end_col + 1L, , drop = FALSE] - grounded[gauge$start_col + 1L, , drop = FALSE]
    rownames(out) <- gauge$link_ids
  } else {
    out <- matrix(0, nrow = gauge$rank_k, ncol = ncol(x))
    for (side in c("start_col", "end_col")) {
      idx <- gauge[[side]]
      keep <- idx > 0L
      if (any(keep)) {
        sums <- rowsum(x[keep, , drop = FALSE], idx[keep], reorder = FALSE)
        dest <- as.integer(rownames(sums))
        sgn <- if (side == "start_col") -1 else 1
        out[dest, ] <- out[dest, , drop = FALSE] + sgn * sums
      }
    }
    rownames(out) <- gauge$keep_param_ids
  }
  if (any(!is.finite(out))) {
    cdid_gmm_abort("Nonfinite value in gauge incidence product.", "cdid_gmm_design_error")
  }
  if (vec_in) drop(out) else out
}

materialize_gmm_gauge_matrix <- function(gauge, max_matrix_entries) {
  if (!inherits(gauge, "cdid_gmm_gauge")) {
    cdid_gmm_abort("Expected `cdid_gmm_gauge` object.", "cdid_gmm_design_error")
  }
  bound <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
  m <- gauge$n_links
  k <- gauge$rank_k
  gmm_safe_product(m, k, bound, "W_G dense matrix")
  W_G <- matrix(0, nrow = m, ncol = k, dimnames = list(gauge$link_ids, gauge$keep_param_ids))
  s_pos <- which(gauge$start_col > 0L)
  if (length(s_pos)) W_G[cbind(s_pos, gauge$start_col[s_pos])] <- -1
  e_pos <- which(gauge$end_col > 0L)
  if (length(e_pos)) W_G[cbind(e_pos, gauge$end_col[e_pos])] <- 1
  W_G
}

# Transform embedded internal coordinates from gauge G to centered minimum-norm
# gauge or another gauge G' using Proposition 1(c),(e), for internal verification only.
transform_gmm_gauge_coordinates <- function(gauge, theta_keep, target_gauge = c("centered", "roots"), new_roots = NULL) {
  target_gauge <- match.arg(target_gauge)
  p <- gauge$n_parameters
  full <- numeric(p)
  full[gauge$keep_params] <- as.double(theta_keep)
  for (r in seq_along(gauge$unanchored_components)) {
    u <- gauge$unanchored_components[r]
    members <- which(gauge$param_component_id == u)
    if (target_gauge == "centered") {
      shift <- sum(full[members]) / length(members)
    } else {
      shift <- full[new_roots[r]]
    }
    full[members] <- full[members] - shift
  }
  full
}

# ---------------------------------------------------------------------------
# Sparse exact target specification, estimability, and BFS path witnesses
# ---------------------------------------------------------------------------

gmm_target_spec <- function(target_id, param_index, coeff) {
  if (!gmm_scalar_character(target_id)) {
    cdid_gmm_abort("`target_id` must be a non-empty character scalar.", "cdid_gmm_target_error")
  }
  if (!is.numeric(param_index) || !is.null(dim(param_index)) || is.object(param_index) ||
      anyNA(param_index) || any(!is.finite(param_index) | param_index != floor(param_index) | param_index < 1 | param_index > .Machine$integer.max) || length(param_index) < 1L) {
    cdid_gmm_abort("`param_index` in target spec must be finite integer indices.", "cdid_gmm_target_error")
  }
  if (!is.numeric(coeff) || !is.null(dim(coeff)) || is.object(coeff) ||
      length(coeff) != length(param_index) || anyNA(coeff) || any(!is.finite(coeff))) {
    cdid_gmm_abort("`coeff` in target spec must be finite numeric values aligned with `param_index`.", "cdid_gmm_target_error")
  }
  list(target_id = target_id, param_index = as.integer(param_index), coeff = as.double(coeff))
}

validate_gmm_targets <- function(design, gauge, targets, bounds) {
  if (missing(targets) || is.null(targets)) {
    cdid_gmm_abort("Explicit `targets` specification is required.", "cdid_gmm_target_error")
  }
  p <- gauge$n_parameters
  m <- gauge$n_links
  k <- gauge$rank_k

  raw_entries <- if (is.matrix(targets)) length(targets) else if (is.data.frame(targets)) nrow(targets) else if (is.list(targets)) sum(vapply(targets, function(x) length(x$coeff), double(1L))) else 0
  J_bound <- if (is.matrix(targets)) ncol(targets) else if (is.data.frame(targets)) min(nrow(targets), bounds$max_targets) else length(targets)
  gmm_safe_product(4, raw_entries, bounds$max_matrix_entries, "raw/canonical target conversion")
  gmm_safe_product(160, min(max(raw_entries,p),1024), bounds$max_matrix_entries, "exact dyadic bit/radix workspace")
  gmm_safe_product(p, J_bound, bounds$max_matrix_entries, "all target coefficient vectors")
  gmm_safe_product(m, J_bound, bounds$max_matrix_entries, "all target witnesses")
  gmm_safe_product(p, p, bounds$max_matrix_entries, "exact path term workspace")

  # Convert supported representations (list of sparse specs, data.frame, or p x J matrix)
  # into a canonical list of sparse target specs.
  spec_list <- list()
  if (is.matrix(targets)) {
    gmm_safe_product(nrow(targets),ncol(targets),bounds$max_matrix_entries,"target validation masks")
    if (!is.numeric(targets) || anyNA(targets) || any(!is.finite(targets)) || nrow(targets) != p || ncol(targets) < 1L) {
      cdid_gmm_abort("Matrix `targets` must be a finite numeric p x J matrix with J >= 1.", "cdid_gmm_target_error")
    }
    J <- ncol(targets)
    if (J > bounds$max_targets) {
      cdid_gmm_abort("Number of targets exceeds `bounds$max_targets`.", "cdid_gmm_bound_error")
    }
    gmm_safe_product(p, J, bounds$max_matrix_entries, "target matrix")
    t_names <- colnames(targets)
    if (is.null(t_names)) t_names <- paste0("target_", seq_len(J))
    if (anyNA(t_names) || any(!nzchar(t_names)) || anyDuplicated(t_names))
      cdid_gmm_abort("Duplicate/empty target matrix names.", "cdid_gmm_target_error")
    for (j in seq_len(J)) {
      nz <- which(targets[, j] != 0)
      spec_list[[j]] <- gmm_target_spec(t_names[j], nz, targets[nz, j])
    }
  } else if (is.data.frame(targets)) {
    if (!all(c("target_id", "param_index", "coeff") %in% names(targets)) || nrow(targets) < 1L) {
      cdid_gmm_abort("Data-frame `targets` must contain `target_id`, `param_index`, and `coeff`.", "cdid_gmm_target_error")
    }
    if (!is.character(targets$target_id) || anyNA(targets$target_id) || any(!nzchar(targets$target_id))) {
      cdid_gmm_abort("`targets$target_id` must be non-empty character strings.", "cdid_gmm_target_error")
    }
    u_ids <- unique(targets$target_id)
    if (length(u_ids) > bounds$max_targets) {
      cdid_gmm_abort("Number of targets exceeds `bounds$max_targets`.", "cdid_gmm_bound_error")
    }
    for (j in seq_along(u_ids)) {
      sub <- targets[targets$target_id == u_ids[j], , drop = FALSE]
      spec_list[[j]] <- gmm_target_spec(u_ids[j], sub$param_index, sub$coeff)
    }
  } else if (is.list(targets) && length(targets) >= 1L) {
    if (length(targets) > bounds$max_targets) {
      cdid_gmm_abort("Number of targets exceeds `bounds$max_targets`.", "cdid_gmm_bound_error")
    }
    for (j in seq_along(targets)) {
      item <- targets[[j]]
      if (!is.list(item) || !all(c("target_id", "param_index", "coeff") %in% names(item))) {
        cdid_gmm_abort("Each element of `targets` list must have `target_id`, `param_index`, and `coeff`.", "cdid_gmm_target_error")
      }
      spec_list[[j]] <- gmm_target_spec(item$target_id, item$param_index, item$coeff)
    }
  } else {
    cdid_gmm_abort("Unsupported `targets` format.", "cdid_gmm_target_error")
  }

  J <- length(spec_list)
  target_ids <- vapply(spec_list, `[[`, character(1L), "target_id")
  if (anyDuplicated(target_ids)) {
    cdid_gmm_abort("Duplicate `target_id` in `targets`.", "cdid_gmm_target_error")
  }

  estimable <- logical(J)
  status <- character(J)
  support_reason <- character(J)
  c_vectors <- vector("list", J)
  witnesses <- vector("list", J)
  unanch <- gauge$unanchored_components
  p_comp <- gauge$param_component_id

  for (j in seq_len(J)) {
    sp <- spec_list[[j]]
    if (any(sp$param_index < 1L | sp$param_index > p)) {
      cdid_gmm_abort(paste0("Target `", sp$target_id, "` references out-of-range parameter index."), "cdid_gmm_target_error")
    }
    # Check exact zero sum on every unanchored component U_r
    failed_iso <- integer(0)
    failed_disc <- integer(0)
    for (u in unanch) {
      members <- which(p_comp == u)
      # Include all raw occurrences for members of U_r so even intra-component
      # expressions like (2^100, 0.1, -2^100) are summed in exact dyadic bins!
      raw_in_u <- sp$coeff[sp$param_index %in% members]
      if (length(raw_in_u) > 0L) {
        ds <- exact_binary64_dyadic_sum(raw_in_u)
        if (!ds$is_zero) {
          if (length(members) == 1L && gauge$is_isolated[members[1L]]) {
            failed_iso <- c(failed_iso, u)
          } else {
            failed_disc <- c(failed_disc, u)
          }
        }
      }
    }

    if (length(failed_iso) == 0L && length(failed_disc) == 0L) {
      # Consolidate repeated coefficients only if their exact sum is representable.
    u_p <- sort(unique(sp$param_index))
    c_vec <- numeric(p)
    for (pi in u_p) {
      vals <- sp$coeff[sp$param_index == pi]
      if (length(vals) == 1L) {
        c_vec[pi] <- vals
      } else {
        ds <- exact_binary64_dyadic_sum(vals)
        if (ds$is_zero) {
          c_vec[pi] <- 0
        } else {
          s_fp <- ds$approx_double
          if (!is.finite(s_fp) || s_fp == 0 || !exact_binary64_dyadic_sum(c(vals, -s_fp))$is_zero) {
            cdid_gmm_abort(
              paste0("Target `", sp$target_id, "` duplicate entries on parameter ", pi, " cannot be represented in one binary64 float without underflow/cancellation."),
              "cdid_gmm_target_error"
            )
          }
          c_vec[pi] <- s_fp
        }
      }
    }
    c_vectors[[j]] <- c_vec

      estimable[j] <- TRUE
      status[j] <- "estimable"
      support_reason[j] <- "exact_zero_sum_on_all_unanchored_components"
      # Build sparse BFS path witness w = sum_v c_v pi(v) satisfying W' w = c
      nz_v <- which(c_vec != 0)
      if (length(nz_v) == 0L) {
        witnesses[[j]] <- data.frame(link_index = integer(0), sign_or_weight = numeric(0))
      } else {
        w_full <- numeric(m)
        w_terms <- vector("list", m)
        for (v_idx in nz_v) {
          vert <- design$parameters$vertex_index[v_idx]
          cert <- root_path_certificate(design, vert)
          if (nrow(cert) > 0L) {
            for (ii in seq_len(nrow(cert))) {
              li <- cert$link_index[ii]
              w_terms[[li]] <- c(w_terms[[li]], c_vec[v_idx] * cert$sign[ii])
            }
          }
        }
        for (li in seq_len(m)) if (length(w_terms[[li]])) {
          ds <- exact_binary64_dyadic_sum(w_terms[[li]])
          v <- ds$approx_double
          if (!is.finite(v) || !exact_binary64_dyadic_sum(c(w_terms[[li]], -v))$is_zero)
            cdid_gmm_abort("Exact path witness is not representable in binary64.", "cdid_gmm_target_error")
          w_full[li] <- v
        }
        if (any(!is.finite(w_full))) {
          cdid_gmm_abort(paste0("Path witness overflow for target `", sp$target_id, "`."), "cdid_gmm_target_error")
        }
        for (pi in seq_len(p)) {
          vals <- c(w_full[design$incidence$end==pi], -w_full[design$incidence$start==pi], -c_vec[pi])
          if (!exact_binary64_dyadic_sum(vals)$is_zero) cdid_gmm_abort("BFS witness does not satisfy the exact original incidence identity.", "cdid_gmm_target_error")
        }
        nz_l <- which(w_full != 0)
        witnesses[[j]] <- data.frame(link_index = nz_l, weight = w_full[nz_l])
      }
    } else {
      c_vectors[j] <- list(NULL)
      estimable[j] <- FALSE
      if (length(failed_iso) > 0L && length(failed_disc) == 0L) {
        status[j] <- "isolated_unsupported"
        support_reason[j] <- paste0("nonzero_sum_on_isolated_component:", paste(failed_iso, collapse = ","))
      } else if (length(failed_disc) > 0L && length(failed_iso) == 0L) {
        status[j] <- "disconnected_component"
        support_reason[j] <- paste0("nonzero_sum_on_unanchored_component:", paste(failed_disc, collapse = ","))
      } else {
        status[j] <- "disconnected_component"
        support_reason[j] <- paste0(
          "nonzero_sum_on_isolated(", paste(failed_iso, collapse = ","),
          ")_and_unanchored(", paste(failed_disc, collapse = ","), ")"
        )
      }
      witnesses[j] <- list(NULL)
    }
  }

  est_idx <- which(estimable)
  J_est <- length(est_idx)
  gmm_safe_product(k, max(J_est, 1L), bounds$max_matrix_entries, "C_gauge target matrix")
  C_gauge <- matrix(0, nrow = k, ncol = J_est,
                    dimnames = list(gauge$keep_param_ids, target_ids[est_idx]))
  for (jj in seq_len(J_est)) {
    j <- est_idx[jj]
    C_gauge[, jj] <- c_vectors[[j]][gauge$keep_params]
  }

  list(
    n_targets = J,
    n_estimable = J_est,
    target_ids = target_ids,
    estimable = estimable,
    status = status,
    support_reason = support_reason,
    c_vectors = c_vectors,
    C_gauge = C_gauge,
    estimable_indices = est_idx,
    witnesses = witnesses
  )
}

# ---------------------------------------------------------------------------
# Preflight allocation and payload/scratch budget accounting
# ---------------------------------------------------------------------------

preflight_gmm_allocations <- function(
  m_R, k, J_est, weighting, bounds,
  n_out_rows = 0L, request_param_cov = FALSE, materialize_target_if = FALSE,
  n_input_rows = n_out_rows, n_parameters = k, certify_error = TRUE,
  covariance_strategy = NULL, n_dropped_triplets = 0L
) {
  max_e <- bounds$max_matrix_entries
  B <- bounds$max_batch_rows

  m2 <- if (weighting == "optimal") gmm_safe_product(m_R, m_R, max_e, "V_R / R (m_R x m_R)") else 0
  mk <- gmm_safe_product(m_R, k, max_e, "W_G / A (m_R x k)")
  qr_comp <- mk + k
  if (qr_comp > max_e) {
    cdid_gmm_abort("Compact QR payload (m_R*k + k) exceeds `max_matrix_entries`.", "cdid_gmm_bound_error")
  }
  kJ <- if (J_est > 0L) gmm_safe_product(k, J_est, max_e, "Z1 / Mc (k x J)") else 0
  mJ <- if (J_est > 0L) gmm_safe_product(m_R, J_est, max_e, "Lambda (m_R x J)") else 0
  Bm <- gmm_safe_product(B, m_R, max_e, "IF batch Phi_b (B x m_R)")
  BJ <- if (J_est > 0L) gmm_safe_product(B, J_est, max_e, "Projected IF batch Psi_b (B x J)") else 0
  J2 <- if (J_est > 0L) gmm_safe_product(J_est, J_est, max_e, "Target covariance accumulator (J x J)") else 0
  k2 <- if (isTRUE(request_param_cov)) gmm_safe_product(k, k, max_e, "Full parameter covariance (k x k)") else 0

  mat_if <- 0
  if (isTRUE(materialize_target_if)) {
    if (!isTRUE(bounds$allow_materialize_if)) {
      cdid_gmm_abort(
        "Target IF materialization requested but `bounds$allow_materialize_if` is FALSE.",
        "cdid_gmm_bound_error"
      )
    }
    mat_if <- gmm_safe_product(n_out_rows, J_est, bounds$max_materialized_if_entries, "Materialized target IF")
    if (mat_if > max_e) {
      cdid_gmm_abort("Materialized target IF exceeds `bounds$max_matrix_entries`.", "cdid_gmm_bound_error")
    }
  }

  payload <- data.frame(
    buffer = c(
      "V_R", "R_chol", "W_G_dense", "A_whitened", "qr_compact",
      "b_whitened", "y_multiplier", "r1_residual",
      "Z1", "Mc", "Lambda",
      "if_batch_phi", "if_batch_psi", "target_cov_accumulator",
      "open_cluster_state", "full_param_cov", "materialized_target_if"
    ),
    rows = c(
      if (weighting == "optimal") m_R else 0,
      if (weighting == "optimal") m_R else 0,
      m_R, m_R, m_R,
      m_R, m_R, m_R,
      k, k, m_R,
      B, B, J_est,
      1, if (isTRUE(request_param_cov)) k else 0,
      if (isTRUE(materialize_target_if)) n_out_rows else 0
    ),
    cols = c(
      if (weighting == "optimal") m_R else 0,
      if (weighting == "optimal") m_R else 0,
      k, k, k,
      1, 1, 1,
      J_est, J_est, J_est,
      m_R, J_est, J_est,
      J_est, if (isTRUE(request_param_cov)) k else 0,
      if (isTRUE(materialize_target_if)) J_est else 0
    ),
    entries = c(
      m2, m2, mk, mk, qr_comp,
      m_R, m_R, m_R,
      kJ, kJ, mJ,
      Bm, BJ, J2,
      J_est, k2, mat_if
    ),
    stringsAsFactors = FALSE
  )
  payload$bytes_8double <- 8 * payload$entries

  kk <- gmm_safe_product(k, k, max_e, "Gram matrix A'A (k x k)")
  scratch <- data.frame(
    buffer = c(
      "Q1_Z1_workspace",
      "guard_eigen_V_vectors",
      "guard_orth_and_residual_V",
      "guard_gram_A_and_eig",
      "lapack_qr_chol_copy_allowance"
    ),
    entries = c(
      mJ,
      m2,
      2 * m2,
      2 * kk + mk,
      m2 + mk + 4 * k
    ),
    stringsAsFactors = FALSE
  )
  if (!certify_error) scratch$entries[2:4] <- 0 # no interval eigen/Gram buffers; computed diagnostics reserved below
  scratch$bytes_8double <- 8 * scratch$entries

  m2_guard <- gmm_safe_product(m_R, m_R, max_e, "guard m_R x m_R matrices")
  guard_base <- gmm_safe_product(1, m2_guard + kk + mk, max_e, "guard base matrices")
  extra <- c(
    target_coefficients = gmm_safe_product(n_parameters, J_est, max_e, "target coefficient cache"),
    target_witnesses = gmm_safe_product(m_R, J_est, max_e, "target witness cache"),
    path_terms = gmm_safe_product(n_parameters, n_parameters, max_e, "path terms"),
    projected_replay_cache = gmm_safe_product(n_input_rows, J_est, max_e, "projected replay cache"),
    replay_identifiers_weights = gmm_safe_product(8, n_input_rows + n_out_rows, max_e, "replay ID/weight buffers"),
    gauge_validation = gmm_safe_product(20, m_R + n_parameters, max_e, "gauge validation buffers"),
    validated_guard_buffers = if (certify_error) gmm_safe_product(12, guard_base, max_e, "interval guard matrices (including native lo/hi buffers)") else 0,
    conventional_diagnostic_buffers = if (certify_error) 0 else gmm_safe_product(4, guard_base, max_e, "computed factor/solve diagnostic matrices"),
    optional_covariance_workspace = if (request_param_cov) gmm_safe_product(3, k, max_e, "covariance streamed solves") else 0)
  if (certify_error) extra <- extra[names(extra) != "conventional_diagnostic_buffers"]
  strategy_work <- 0
  if (!is.null(covariance_strategy) && identical(weighting, "optimal")) {
    # Additional simultaneously reachable D26 payload. Reserve copies and
    # expression temporaries as well as retained matrices; the base diagnostic
    # reservation above is deliberately not borrowed for these new buffers.
    # LAPACK eigen workspace includes its copied input/output and linear work.
    nt <- gmm_int_bound(n_dropped_triplets, "n_dropped_triplets", allow_zero = TRUE)
    P <- function(a, b, label) gmm_safe_product(a, b, max_e, label)
    extra <- c(extra,
      D26_reconstruction_G_C = P(2, m2, "D26 G and C"),
      D26_reconstruction_scratch = P(4, m2, "D26 Gram submatrix/product/copy scratch"),
      D26_sparse_triplets_indices = P(12, nt, "D26 triplet/order/group buffers"),
      D26_transformed_covariance_products = P(4, m2, "D26 B, transform products and copies"),
      D26_eigensystem_matrix_workspace = P(8, m2, "D26 eigenvectors, copies, workspace and checked QR/Rayleigh refinement"),
      D26_eigenvalues_linear_workspace = P(40, m_R, "D26 eigenvalues and LAPACK linear workspace"),
      D26_eigen_residual_norm_buffers = P(6, m2, "D26 reconstructed eigen product, residual and scaled norms"),
      D26_ridge_factor_residual_buffers = P(6, m2, "D26 shifted covariance, factor, product/residual/copies"),
      D26_design_solve_buffers = P(6, mk, "D26 transformed design and solve buffers"),
      D26_loading_solve_buffers = P(8, mJ, "D26 loading transforms, triangular solves and residuals"),
      D26_coefficient_solve_buffers = P(6, kJ, "D26 coefficient solves and residuals"),
      D26_design_gram_workspace = P(4, kk, "D26 design Gram and eigensystem"),
      D26_vector_workspace = P(16, m_R + k + J_est, "D26 solve/vector/index scratch"),
      D26_gram_provenance_snapshot = P(1, m2, "D26 sealed covariance snapshot allowance"))
    # Conservative counted algebra work envelope, not CPU instructions/time
    # or a bound on LAPACK iteration counts. The pipeline composes and admits
    # this charge under its existing finite work limit before covariance/GMM.
    W <- function(a, b, label) gmm_safe_product(a, b, NULL, label)
    strategy_work <- W(1, sum(c(
      W(160, W(m2, m_R, "D26 cubic matrix work"), "D26 transforms/eigen/refinement/factors work"),
      W(64, W(m2, k + J_est + 1, "D26 matrix solve work"), "D26 solve/residual work"),
      W(64, W(m_R, kk + J2 + 1, "D26 QR/loading work"), "D26 QR/Gram work"),
      W(16, W(nt, nt + m_R + 1, "D26 sparse assembly work"), "D26 triplet/group work")
    )), "D26 counted covariance-strategy work envelope")
  }
  total_payload_entries <- sum(c(payload$entries, extra, scratch$entries))
  if (!is.finite(total_payload_entries) || total_payload_entries > 2^53-1) cdid_gmm_abort("Aggregate allocation count overflow.", "cdid_gmm_bound_error")
  if (total_payload_entries > max_e) {
    cdid_gmm_abort(
      paste0("Total simultaneous double payload entries (", total_payload_entries,
             ") exceed `bounds$max_matrix_entries` (", max_e, ")."),
      "cdid_gmm_bound_error",
      diagnostics = list(payload = payload, scratch = scratch, extra_payload_entries = extra,
        total_payload_entries = total_payload_entries, max_matrix_entries = max_e)
    )
  }

  list(
    payload = payload,
    scratch = scratch,
    extra_payload_entries = extra,
    covariance_strategy_work_ops = strategy_work,
    total_payload_entries = total_payload_entries,
    total_payload_bytes = 8 * total_payload_entries,
    total_scratch_entries = sum(scratch$entries),
    total_scratch_bytes = 8 * sum(scratch$entries),
    forbidden_allocations = list(
      V_inverse_allocated = FALSE,
      full_n_by_m_if_allocated = FALSE,
      dense_null_basis_allocated = FALSE
    ),
    peak_memory_promised = FALSE,
    note = "Figures are payload and explicit scratch/workspace budgets (8-byte doubles), not process peak RSS."
  )
}

# ---------------------------------------------------------------------------
# Conservative sufficient numerical guard:
#   1. Symmetric eigen-residual + orthogonality Weyl/Rayleigh bounds
#   2. Cholesky factor residual bound (rejects indefinite/singular chol passes)
#   3. Whitened design singular-value lower bound & KKT inverse 2-norm bound
#   4. Primal & adjoint KKT residuals with IEEE binary64 roundoff margins
#   5. Rigorous target forward-error bound (including ||s||_2 ||K^{-1}||_2 ||r||_2)
#      explicitly distinguished from ordinary computed indicator
# ---------------------------------------------------------------------------

gmm_native_dll_cache <- new.env(parent = emptyenv())

gmm_require_native_dll <- function() {
  owner <- environment(sys.function())
  dlls <- if (isNamespace(owner)) getNamespaceInfo(owner, "DLLs") else base::getLoadedDLLs()
  dll <- dlls[["cdid"]]
  if (is.null(dll)) {
    cdid_gmm_abort(
      "Production native C backend (`cdid` DLL) is not loaded; automatic pure-R fallback is forbidden.",
      c("cdid_gmm_native_error", "cdid_gmm_accuracy_error")
    )
  }
  info <- if (inherits(dll, "DLLInfo")) dll[["info"]] else NULL
  cache <- get0("gmm_native_dll_cache", envir = owner, mode = "environment", inherits = FALSE)
  if (!is.null(info) && is.environment(cache) && identical(cache$owner, owner) &&
      identical(cache$info, info) && is.list(cache$syms)) {
    syms <- cache$syms
    # C_gmm_next already checks live arithmetic in older valid backends too;
    # their identity-only C_gmm_sexp_aliases cannot serve as the warm guard.
    ok <- tryCatch(.Call(syms$C_gmm_next, 1.0, TRUE), error = function(e) {
      if (inherits(e, "cdid_gmm_refusal")) stop(e)
      NULL
    })
    if (identical(ok, 1.0 + .Machine$double.eps)) return(syms)
    cache$owner <- NULL
    cache$info <- NULL
    cache$syms <- NULL
  }
    reg <- tryCatch(base::getDLLRegisteredRoutines(dll, addNames = TRUE)[['.Call']], error = function(e) NULL)
    req_syms <- c(
      C_gmm_native_status = 0L,
      C_gmm_test_arithmetic_probe = 5L,
      C_gmm_sexp_aliases = 2L,
      C_gmm_next = 2L,
      C_gmm_iv_dot = 3L,
      C_gmm_iv_norm = 3L,
      C_gmm_iv_mat_residual = 4L,
      C_gmm_iv_eigen_residual = 4L,
      C_gmm_iv_mat_product = 5L
    )
    if (is.null(reg) || !all(names(req_syms) %in% names(reg))) {
      cdid_gmm_abort(
        "Loaded `cdid` DLL is missing required registered `.Call` routines.",
        c("cdid_gmm_native_error", "cdid_gmm_accuracy_error")
      )
    }
    syms <- vector("list", length(req_syms))
    names(syms) <- names(req_syms)
    for (nm in names(req_syms)) {
      si <- tryCatch(base::getNativeSymbolInfo(nm, PACKAGE = dll), error = function(e) NULL)
      if (is.null(si) || !inherits(si, "NativeSymbolInfo") || is.null(si$address) ||
          !identical(as.integer(si$numParameters), req_syms[[nm]])) {
        cdid_gmm_abort(
          paste0("Native symbol registration/forcing verification failed for `", nm, "`."),
          c("cdid_gmm_native_error", "cdid_gmm_accuracy_error")
        )
      }
      syms[[nm]] <- si$address
    }
    st <- .Call(syms$C_gmm_native_status)
    if (!is.list(st) || !identical(st$backend, "native_c_ieee754_binary64") ||
        !isTRUE(st$c_registration) || !isTRUE(st$force_symbols) ||
        !identical(st$fallback_permitted, FALSE) ||
        !identical(st$flt_radix, 2L) || !identical(st$dbl_mant_dig, 53L) ||
        !identical(st$rounding_mode, "FE_TONEAREST") ||
        !isTRUE(st$gradual_underflow_verified) || !isTRUE(st$volatile_fp_stores)) {
      cdid_gmm_abort(
        "Native C backend arithmetic/registration status check failed.",
        c("cdid_gmm_native_error", "cdid_gmm_accuracy_error")
      )
    }
  if (!is.null(info) && is.environment(cache)) {
    cache$owner <- owner
    cache$info <- info
    cache$syms <- syms
  }
  syms
}

gmm_native_backend_info <- function() {
  syms <- gmm_require_native_dll()
  st <- .Call(syms$C_gmm_native_status)
  owner <- environment(sys.function())
  dll <- if (isNamespace(owner)) getNamespaceInfo(owner, "DLLs")[["cdid"]] else base::getLoadedDLLs()[["cdid"]]
  st$dll_name <- "cdid"
  st$dll_loaded <- TRUE
  st$dll_path <- if (!is.null(dll[["path"]])) as.character(dll[["path"]]) else NA_character_
  st$outer_derived_controls <- "R_scalar_interval"
  st$pure_r_reference_role <- "frozen_independent_test_reference_only"
  st
}

gmm_sexp_aliases <- function(a, b) {
  syms <- gmm_require_native_dll()
  isTRUE(.Call(syms$C_gmm_sexp_aliases, a, b))
}

# Scalar outward rounding. Assumptions: IEEE binary64, round-to-nearest ties-to-even,
# correctly rounded scalar +,-,*,/,sqrt, gradual underflow, no FTZ/DAZ/reassociation.
# BLAS/LAPACK candidates are untrusted; certificates use no BLAS reductions.
gmm_next <- function(x, up = TRUE) {
  if (length(x) != 1L || is.na(x)) cdid_gmm_abort("Invalid interval endpoint.", "cdid_gmm_accuracy_error")
  if (x == 0) return(if (up) 2^-1074 else -2^-1074)
  if (is.infinite(x)) {
    if ((x > 0 && up) || (x < 0 && !up)) return(x)
    return(if (x > 0) .Machine$double.xmax else -.Machine$double.xmax)
  }
  b <- as.integer(writeBin(as.double(x), raw(), size = 8L, endian = "little"))
  inc <- (x > 0) == up
  for (i in seq_len(8L)) {
    if (inc) {
      if (b[i] < 255L) { b[i] <- b[i] + 1L; break }
      b[i] <- 0L
    } else {
      if (b[i] > 0L) { b[i] <- b[i] - 1L; break }
      b[i] <- 255L
    }
  }
  readBin(as.raw(b), "double", n = 1L, size = 8L, endian = "little")
}
gmm_iv <- function(x) c(x, x)
gmm_iv_check <- function(z) {
  if (length(z) != 2L || any(!is.finite(z)) || z[1L] > z[2L])
    cdid_gmm_abort("Validated scalar arithmetic overflowed or became inconclusive.", "cdid_gmm_accuracy_error")
  z
}
gmm_iv_add <- function(a, b) {
  if (identical(unname(a), c(0, 0))) return(b)
  if (identical(unname(b), c(0, 0))) return(a)
  gmm_iv_check(c(gmm_next(a[1L] + b[1L], FALSE), gmm_next(a[2L] + b[2L])))
}
gmm_iv_neg <- function(a) -rev(a)
gmm_iv_sub <- function(a, b) gmm_iv_add(a, gmm_iv_neg(b))
gmm_iv_mul <- function(a, b) {
  if (all(a == 0) || all(b == 0)) return(c(0, 0))
  z <- c(a[1L]*b[1L], a[1L]*b[2L], a[2L]*b[1L], a[2L]*b[2L])
  gmm_iv_check(c(gmm_next(min(z), FALSE), gmm_next(max(z))))
}
gmm_iv_div <- function(a, b) {
  if (b[1L] <= 0 && b[2L] >= 0) cdid_gmm_abort("Interval division crosses zero.", "cdid_gmm_accuracy_error")
  z <- c(a[1L]/b[1L], a[1L]/b[2L], a[2L]/b[1L], a[2L]/b[2L])
  gmm_iv_check(c(gmm_next(min(z), FALSE), gmm_next(max(z))))
}
gmm_iv_sqrt <- function(a) {
  if (a[1L] < 0) cdid_gmm_abort("Negative interval square root.", "cdid_gmm_accuracy_error")
  gmm_iv_check(c(if (a[1L] == 0) 0 else max(0, gmm_next(sqrt(a[1L]), FALSE)),
                 if (a[2L] == 0) 0 else gmm_next(sqrt(a[2L]))))
}
gmm_iv_abs <- function(a) c(if (a[1L] <= 0 && a[2L] >= 0) 0 else min(abs(a)), max(abs(a)))
gmm_iv_dot <- function(x, y, max_matrix_entries = NULL) {
  syms <- gmm_require_native_dll()
  if (!is.null(max_matrix_entries)) {
    max_matrix_entries <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
  }
  if (is.null(x) || !is.atomic(x) || is.object(x) || !is.numeric(x) || !is.null(dim(x)) ||
      is.null(y) || !is.atomic(y) || is.object(y) || !is.numeric(y) || !is.null(dim(y))) {
    cdid_gmm_abort("`x` and `y` in `gmm_iv_dot` must be non-object atomic numeric vectors.", "cdid_gmm_accuracy_error")
  }
  nx <- length(x); ny <- length(y)
  if (nx < 1L || ny < 1L || nx != ny) {
    cdid_gmm_abort("`x` and `y` in `gmm_iv_dot` must have identical nonempty lengths.", "cdid_gmm_bound_error")
  }
  gmm_safe_product(nx, 1, max_matrix_entries, "gmm_iv_dot(x)")
  gmm_safe_product(2, nx, max_matrix_entries, "gmm_iv_dot(x + y)")
  if (anyNA(x) || any(!is.finite(x)) || anyNA(y) || any(!is.finite(y))) {
    cdid_gmm_abort("Nonfinite element in `gmm_iv_dot`.", "cdid_gmm_accuracy_error")
  }
  if (typeof(x) != "double") x <- as.double(x)
  if (typeof(y) != "double") y <- as.double(y)
  .Call(syms$C_gmm_iv_dot, x, y, max_matrix_entries)
}
# Scaled validated norm, avoiding squaring huge/small values before scaling.
gmm_iv_norm <- function(lo, hi = lo, max_matrix_entries = NULL) {
  syms <- gmm_require_native_dll()
  has_hi <- !missing(hi) && !is.null(hi)
  if (!is.null(max_matrix_entries)) {
    max_matrix_entries <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
  }
  if (is.null(lo) || !is.atomic(lo) || is.object(lo) || !is.numeric(lo)) {
    cdid_gmm_abort("`lo` in `gmm_iv_norm` must be a non-object atomic numeric vector or matrix.", "cdid_gmm_accuracy_error")
  }
  n <- length(lo)
  if (n < 1L) {
    cdid_gmm_abort("Invalid empty input length in `gmm_iv_norm`.", "cdid_gmm_bound_error")
  }
  d_lo <- dim(lo)
  if (!is.null(d_lo)) {
    if (length(d_lo) != 2L || anyNA(d_lo) || any(d_lo < 1L)) {
      cdid_gmm_abort("Invalid matrix dimensions on `lo` in `gmm_iv_norm`.", "cdid_gmm_bound_error")
    }
    prod_lo <- gmm_safe_product(d_lo[1L], d_lo[2L], max_matrix_entries, "gmm_iv_norm(lo)")
    if (prod_lo != n) cdid_gmm_abort("Dimension product on `lo` disagrees with length.", "cdid_gmm_bound_error")
  } else {
    gmm_safe_product(n, 1, max_matrix_entries, "gmm_iv_norm(lo)")
  }
  if (has_hi) {
    if (!is.atomic(hi) || is.object(hi) || !is.numeric(hi)) {
      cdid_gmm_abort("`hi` in `gmm_iv_norm` must be a non-object atomic numeric vector or matrix.", "cdid_gmm_accuracy_error")
    }
    if (length(hi) != n || !identical(dim(lo), dim(hi))) {
      cdid_gmm_abort("`lo` and `hi` in `gmm_iv_norm` must have identical lengths and dimensions.", "cdid_gmm_bound_error")
    }
    gmm_safe_product(2, n, max_matrix_entries, "gmm_iv_norm(lo + hi)")
    }
  if (anyNA(lo) || any(!is.finite(lo))) {
    cdid_gmm_abort("Nonfinite element in `lo` for `gmm_iv_norm`.", "cdid_gmm_accuracy_error")
  }
  if (typeof(lo) != "double") storage.mode(lo) <- "double"
  if (has_hi) {
    if (anyNA(hi) || any(!is.finite(hi))) {
      cdid_gmm_abort("Nonfinite element in `hi` for `gmm_iv_norm`.", "cdid_gmm_accuracy_error")
    }
    if (any(lo > hi)) {
      cdid_gmm_abort("Inverted interval endpoints (`lo > hi`) in `gmm_iv_norm`.", "cdid_gmm_accuracy_error")
    }
    if (typeof(hi) != "double") storage.mode(hi) <- "double"
  }
  .Call(syms$C_gmm_iv_norm, lo, if (has_hi) hi else NULL, max_matrix_entries)
}
gmm_iv_mat_residual <- function(X, Y, D, max_matrix_entries = NULL) {
  syms <- gmm_require_native_dll()
  if (!is.null(max_matrix_entries)) {
    max_matrix_entries <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
  }
  for (item in list(list(M = X, nm = "X"), list(M = Y, nm = "Y"), list(M = D, nm = "D"))) {
    M <- item$M; nm <- item$nm
    if (is.null(M) || !is.atomic(M) || is.object(M) || !is.numeric(M) || !is.matrix(M)) {
      cdid_gmm_abort(paste0("`", nm, "` in `gmm_iv_mat_residual` must be a non-object numeric matrix."), "cdid_gmm_accuracy_error")
    }
    nr <- nrow(M); nc <- ncol(M)
    if (is.na(nr) || is.na(nc) || nr < 1L || nc < 1L) {
      cdid_gmm_abort(paste0("`", nm, "` in `gmm_iv_mat_residual` has invalid matrix dimensions."), "cdid_gmm_bound_error")
    }
    prod_m <- gmm_safe_product(nr, nc, max_matrix_entries, paste0("gmm_iv_mat_residual(", nm, ")"))
    if (prod_m != length(M)) {
      cdid_gmm_abort(paste0("`", nm, "` matrix dimension product disagrees with length."), "cdid_gmm_bound_error")
    }
  }
  if (ncol(X) != nrow(Y) || nrow(D) != nrow(X) || ncol(D) != ncol(Y)) {
    cdid_gmm_abort("Conformability error in `gmm_iv_mat_residual`: require ncol(X)==nrow(Y), dim(D)==c(nrow(X),ncol(Y)).", "cdid_gmm_bound_error")
  }
  out_two <- gmm_safe_product(2, length(D), max_matrix_entries, "gmm_iv_mat_residual interval output (lo + hi)")
  if (anyNA(X) || any(!is.finite(X)) || anyNA(Y) || any(!is.finite(Y)) || anyNA(D) || any(!is.finite(D))) {
    cdid_gmm_abort("Nonfinite element in `gmm_iv_mat_residual`.", "cdid_gmm_accuracy_error")
  }
  if (typeof(X) != "double") storage.mode(X) <- "double"
  if (typeof(Y) != "double") storage.mode(Y) <- "double"
  if (typeof(D) != "double") storage.mode(D) <- "double"
  .Call(syms$C_gmm_iv_mat_residual, X, Y, D, max_matrix_entries)
}
# Directed interval matrix product [xlo, xhi] (nr x nk) %*% [ylo, yhi] (nk x nc)
# on the native backend (`C_gmm_iv_mat_product`): entry (i, j) is the same-order
# recurrence z_k = gmm_iv_add(z_{k-1}, gmm_iv_mul(X[i,k], Y[k,j])) from z_0 = [0, 0]
# with the frozen R reference's zero shortcuts (an all-zero row of X or column
# of Y leaves [0, 0]). The frozen pure-R recurrence is a TEST reference only
# (tests/testthat/reference/cp02_interval_product_R_reference.R); there is no R
# fallback. Engineering note: docs/theory/cp02_interval_product_native.tex.
gmm_iv_mat_product <- function(xlo, xhi, ylo, yhi = ylo, max_matrix_entries = NULL) {
  syms <- gmm_require_native_dll()
  if (!is.null(max_matrix_entries)) {
    max_matrix_entries <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
  }
  for (item in list(list(M = xlo, nm = "xlo"), list(M = xhi, nm = "xhi"), list(M = ylo, nm = "ylo"), list(M = yhi, nm = "yhi"))) {
    M <- item$M; nm <- item$nm
    if (is.null(M) || !is.atomic(M) || is.object(M) || !is.numeric(M) || !is.matrix(M)) {
      cdid_gmm_abort(paste0("`", nm, "` in `gmm_iv_mat_product` must be a non-object numeric matrix."), "cdid_gmm_accuracy_error")
    }
    nr <- nrow(M); nc <- ncol(M)
    if (is.na(nr) || is.na(nc) || nr < 1L || nc < 1L) {
      cdid_gmm_abort(paste0("`", nm, "` in `gmm_iv_mat_product` has invalid matrix dimensions."), "cdid_gmm_bound_error")
    }
    prod_m <- gmm_safe_product(nr, nc, max_matrix_entries, paste0("gmm_iv_mat_product(", nm, ")"))
    if (prod_m != length(M)) {
      cdid_gmm_abort(paste0("`", nm, "` matrix dimension product disagrees with length."), "cdid_gmm_bound_error")
    }
  }
  if (!identical(dim(xlo), dim(xhi))) {
    cdid_gmm_abort("`xlo` and `xhi` in `gmm_iv_mat_product` must have identical dimensions.", "cdid_gmm_bound_error")
  }
  if (!identical(dim(ylo), dim(yhi))) {
    cdid_gmm_abort("`ylo` and `yhi` in `gmm_iv_mat_product` must have identical dimensions.", "cdid_gmm_bound_error")
  }
  if (ncol(xlo) != nrow(ylo)) {
    cdid_gmm_abort("Conformability error in `gmm_iv_mat_product`: require ncol(x) == nrow(y).", "cdid_gmm_bound_error")
  }
  out_len <- gmm_safe_product(nrow(xlo), ncol(ylo), max_matrix_entries, "gmm_iv_mat_product output (nrow(x) x ncol(y))")
  gmm_safe_product(2, out_len, max_matrix_entries, "gmm_iv_mat_product interval output (lo + hi)")
  if (anyNA(xlo) || any(!is.finite(xlo)) || anyNA(xhi) || any(!is.finite(xhi)) ||
      anyNA(ylo) || any(!is.finite(ylo)) || anyNA(yhi) || any(!is.finite(yhi))) {
    cdid_gmm_abort("Nonfinite element in `gmm_iv_mat_product`.", "cdid_gmm_accuracy_error")
  }
  if (any(xlo > xhi)) {
    cdid_gmm_abort("Inverted interval endpoints (`lo > hi`) in `gmm_iv_mat_product` operand `x`.", "cdid_gmm_accuracy_error")
  }
  if (any(ylo > yhi)) {
    cdid_gmm_abort("Inverted interval endpoints (`lo > hi`) in `gmm_iv_mat_product` operand `y`.", "cdid_gmm_accuracy_error")
  }
  if (typeof(xlo) != "double") storage.mode(xlo) <- "double"
  if (typeof(xhi) != "double") storage.mode(xhi) <- "double"
  if (typeof(ylo) != "double") storage.mode(ylo) <- "double"
  if (typeof(yhi) != "double") storage.mode(yhi) <- "double"
  .Call(syms$C_gmm_iv_mat_product, xlo, xhi, ylo, yhi, max_matrix_entries)
}
gmm_iv_ratio_upper <- function(num, den) {
  if (num == 0) return(0)
  if (den <= 0) return(Inf)
  gmm_iv_div(gmm_iv(num), gmm_iv(den))[2L]
}
gmm_iv_residual <- function(X, y, D, max_matrix_entries = NULL) {
  syms <- gmm_require_native_dll()
  if (is.null(y) || !is.atomic(y) || is.object(y) || !is.numeric(y) || !is.null(dim(y)) ||
      is.null(D) || !is.atomic(D) || is.object(D) || !is.numeric(D) || !is.null(dim(D))) {
    cdid_gmm_abort("`y` and `D` in `gmm_iv_residual` must be non-object numeric vectors.", "cdid_gmm_accuracy_error")
  }
  if (!is.matrix(X) || length(y) != ncol(X) || length(D) != nrow(X))
    cdid_gmm_abort("Nonconformable residual vectors.", "cdid_gmm_bound_error")
  gmm_safe_product(length(y), 1, max_matrix_entries, "residual vector y")
  gmm_safe_product(length(D), 1, max_matrix_entries, "residual vector D")
  gmm_safe_product(length(X), 1, max_matrix_entries, "residual matrix X")
  gmm_safe_product(2, length(D), max_matrix_entries, "residual interval output")
  z <- gmm_iv_mat_residual(X, matrix(as.double(y), ncol = 1L), matrix(as.double(D), ncol = 1L),
                           max_matrix_entries = max_matrix_entries)
  list(lo = drop(z$lo), hi = drop(z$hi), abs = pmax(abs(z$lo), abs(z$hi)), norm = z$norm)
}
certify_symmetric_eigen_bounds <- function(S, label = "matrix", max_matrix_entries = NULL) {
  syms <- gmm_require_native_dll()
  if (!is.null(max_matrix_entries)) {
    max_matrix_entries <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
  }
  if (!is.matrix(S) || is.object(S) || !is.numeric(S) || nrow(S) < 1L || nrow(S) != ncol(S))
    cdid_gmm_abort("Invalid spectral certificate input.", "cdid_gmm_pd_error")
  m <- nrow(S)
  m2_s <- gmm_safe_product(m, m, max_matrix_entries, paste0("spectral matrix ", label))
  gmm_safe_product(6, m2_s, max_matrix_entries, paste0("spectral certificate workspace for ", label))
  if (anyNA(S) || any(!is.finite(S)))
    cdid_gmm_abort("Invalid spectral certificate input.", "cdid_gmm_pd_error")
  for (j in seq_len(ncol(S))) if (!identical(as.double(S[, j]), as.double(S[j, ])))
    cdid_gmm_abort("Spectral input is not exactly symmetric.", "cdid_gmm_pd_error")
  if (typeof(S) != "double") storage.mode(S) <- "double"
  eg <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eg) || any(!is.finite(eg$values)) || any(!is.finite(eg$vectors)))
    cdid_gmm_abort(paste0("Eigen candidate failed: ", label), "cdid_gmm_pd_error")
  Q <- eg$vectors; lam <- as.double(eg$values)
  if (typeof(Q) != "double") storage.mode(Q) <- "double"
  G <- gmm_iv_mat_residual(t(Q), Q, diag(m), max_matrix_entries = max_matrix_entries)
  E_eig <- .Call(syms$C_gmm_iv_eigen_residual, S, Q, lam, max_matrix_entries)
  lo <- E_eig$lo; hi <- E_eig$hi
  delta <- G$norm[2L]; eps <- E_eig$norm[2L]
  ans <- list(valid = FALSE, lambda_min_est = min(lam), lambda_max_est = max(lam),
              lambda_min_cert = -Inf, lambda_max_cert = Inf,
              delta_Q_est = norm(crossprod(Q) - diag(m), "F"), delta_Q_cert = delta,
              eps_eig_est = norm(S %*% Q - sweep(Q, 2L, lam, `*`), "F"), eps_eig_cert = eps,
              backend = "native_c_ieee754_binary64")
  if (delta >= 1) return(ans)
  den <- c(gmm_iv_sub(gmm_iv(1), gmm_iv(delta))[1L], gmm_iv_add(gmm_iv(1), gmm_iv(delta))[2L])
  pert <- gmm_iv_add(gmm_iv_mul(gmm_iv_sqrt(gmm_iv_add(gmm_iv(1), gmm_iv(delta))), gmm_iv(eps)),
                     gmm_iv_mul(gmm_iv(delta), gmm_iv(max(abs(lam)))))[2L]
  ans$lambda_min_cert <- gmm_iv_div(gmm_iv_sub(gmm_iv(min(lam)), gmm_iv(pert)), den)[1L]
  ans$lambda_max_cert <- gmm_iv_div(gmm_iv_add(gmm_iv(max(lam)), gmm_iv(pert)), den)[2L]
  ans$pert_cert <- pert; ans$valid <- TRUE
  ans
}

certify_gmm_numerical_guard <- function(
  V_R, R_chol, W_G, A_whitened, U_qr, piv, mhat_R,
  theta, y, C_gauge, Lambda, Mc, Z1, weighting, accuracy,
  max_matrix_entries = NULL
) {
  conventional <- identical(accuracy$method, "conventional_chol_qr")
  nb_info <- if (conventional) list(backend = "base_R_BLAS_LAPACK_computed_diagnostics") else gmm_native_backend_info()
  if (!is.null(max_matrix_entries)) {
    max_matrix_entries <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
  }
  if (missing(weighting) || !gmm_scalar_character(weighting) || !weighting %in% c("optimal", "identity")) {
    cdid_gmm_abort("`weighting` must be explicitly 'optimal' or 'identity'.", "cdid_gmm_control_error")
  }
  accuracy <- validate_gmm_accuracy_control(accuracy, weighting)
  if (!is.matrix(W_G) || is.object(W_G) || !is.numeric(W_G) || nrow(W_G) < 1L || ncol(W_G) < 1L ||
      !is.matrix(A_whitened) || is.object(A_whitened) || !is.numeric(A_whitened) ||
      !identical(dim(W_G), dim(A_whitened)) ||
      !is.matrix(U_qr) || is.object(U_qr) || !is.numeric(U_qr) ||
      nrow(U_qr) != ncol(W_G) || ncol(U_qr) != ncol(W_G) ||
      !is.matrix(C_gauge) || is.object(C_gauge) || !is.numeric(C_gauge) || nrow(C_gauge) != ncol(W_G) ||
      !is.matrix(Lambda) || is.object(Lambda) || !is.numeric(Lambda) ||
      nrow(Lambda) != nrow(W_G) || ncol(Lambda) != ncol(C_gauge) ||
      !is.matrix(Mc) || is.object(Mc) || !is.numeric(Mc) || !identical(dim(Mc), dim(C_gauge)) ||
      !is.matrix(Z1) || is.object(Z1) || !is.numeric(Z1) || !identical(dim(Z1), dim(C_gauge))) {
    cdid_gmm_abort("Malformed matrix dimensions or types in `certify_gmm_numerical_guard`.", "cdid_gmm_accuracy_error")
  }
  m <- nrow(W_G); k <- ncol(W_G); J_est <- ncol(C_gauge)
  if (!is.atomic(mhat_R) || is.object(mhat_R) || !is.numeric(mhat_R) || !is.null(dim(mhat_R)) || length(mhat_R) != m ||
      !is.atomic(theta) || is.object(theta) || !is.numeric(theta) || !is.null(dim(theta)) || length(theta) != k ||
      !is.atomic(y) || is.object(y) || !is.numeric(y) || !is.null(dim(y)) || length(y) != m) {
    cdid_gmm_abort("Malformed vector dimensions or types in `certify_gmm_numerical_guard`.", "cdid_gmm_accuracy_error")
  }
  if (conventional) {
    # D23: the same stored Cholesky/QR candidate, using the existing condition
    # and factor-residual thresholds on ordinary computed diagnostics. No
    # forward-error certificate or relative-residual rejection is requested.
    gmm_safe_product(4, m*m + k*k + m*k, max_matrix_entries, "computed factor/solve diagnostic matrices")
    candidates <- list(W_G, A_whitened, U_qr, mhat_R, theta, y, C_gauge, Lambda, Mc, Z1)
    if (weighting == "optimal") candidates <- c(candidates, list(V_R, R_chol))
    if (any(vapply(candidates, function(x) !is.numeric(x) || any(!is.finite(x)), logical(1L))))
      cdid_gmm_abort("Nonfinite numerical candidate.", "cdid_gmm_accuracy_error")
    if (!identical(sort(as.integer(piv)), seq_len(k))) cdid_gmm_abort("Invalid QR pivot permutation.", "cdid_gmm_qr_error")
    if (weighting == "optimal") {
      if (!is.matrix(V_R) || !identical(dim(V_R), c(m,m)) ||
          !is.matrix(R_chol) || !identical(dim(R_chol), c(m,m)) || any(diag(R_chol) <= 0))
        cdid_gmm_abort("Invalid Cholesky factor/covariance.", "cdid_gmm_pd_error")
      ev <- eigen(V_R, symmetric = TRUE, only.values = TRUE)$values
      rc_R <- rcond(R_chol, triangular = TRUE)
      factor_rel <- norm(crossprod(R_chol) - V_R, "F") / norm(V_R, "F")
      whitening_rel <- norm(crossprod(R_chol, A_whitened) - W_G, "F") / norm(W_G, "F")
      if (any(!is.finite(ev)) || min(ev) <= 0 || min(ev)/max(ev) < accuracy$min_certified_lambda_rel ||
          !is.finite(rc_R) || rc_R < accuracy$min_rcond_R ||
          !is.finite(factor_rel) || factor_rel > accuracy$max_factor_residual_rel)
        cdid_gmm_abort("Weight covariance / Cholesky conditioning or factor residual refused.", "cdid_gmm_pd_error")
      if (!is.finite(whitening_rel) || whitening_rel > accuracy$max_factor_residual_rel)
        cdid_gmm_abort("Whitened design residual exceeds caller bound.", "cdid_gmm_qr_error")
      Vcrit <- V_R
    } else {
      if (!identical(unname(A_whitened), unname(W_G))) cdid_gmm_abort("Identity criterion requires A = W.", "cdid_gmm_qr_error")
      ev <- NULL; rc_R <- factor_rel <- whitening_rel <- NULL
      Vcrit <- diag(m)
    }
    rc_U <- rcond(U_qr, triangular = TRUE)
    gram <- crossprod(A_whitened)
    if (any(!is.finite(gram))) cdid_gmm_abort("Gram candidate overflowed.", "cdid_gmm_qr_error")
    ea <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
    sigma_rel <- sqrt(max(0, min(ea))/max(ea))
    if (any(!is.finite(ea)) || min(ea) <= 0 || !is.finite(sigma_rel) || sigma_rel < accuracy$min_sigma_rel_A ||
        !is.finite(rc_U) || rc_U < accuracy$min_rcond_U)
      cdid_gmm_abort("Whitened design / QR conditioning refused.", "cdid_gmm_qr_error")
    r1 <- mhat_R - drop(Vcrit %*% y) - drop(W_G %*% theta)
    r2 <- -drop(crossprod(W_G, y))
    scale1 <- norm(Vcrit, "F")*sqrt(sum(y*y)) + norm(W_G, "F")*sqrt(sum(theta*theta)) + sqrt(sum(mhat_R*mhat_R))
    scale2 <- norm(W_G, "F")*sqrt(sum(y*y))
    ratio <- function(r, scale) if (scale > 0) sqrt(sum(r*r))/scale else if (all(r == 0)) 0 else Inf
    s1 <- W_G %*% Mc - Vcrit %*% Lambda
    s2 <- C_gauge - crossprod(W_G, Lambda)
    tau <- drop(crossprod(C_gauge, theta))
    if (any(!is.finite(c(r1, r2, s1, s2, tau)))) cdid_gmm_abort("Nonfinite solve/loading diagnostics.", "cdid_gmm_accuracy_error")
    return(list(certified = FALSE, point_certified = FALSE, loading_certified = FALSE,
      loading_budget_supplied = FALSE, loading_validation_status = "not_requested_D23",
      point_validation_status = "not_requested_D23", complete_inference_certified = FALSE,
      full_inference_status = "complete_per_fit_error_certification_not_requested_or_performed",
      backend = nb_info$backend, native_backend = nb_info, method = accuracy$method, weighting = weighting,
      computed_target_values = tau, propagation_bounds = NULL,
      conditioning = list(lambda_V = ev, rcond_R_1norm = rc_R, rcond_U_1norm = rc_U, sigma_rel_A_est_2norm = sigma_rel),
      factorization = list(cholesky_relative_residual = factor_rel, whitening_relative_residual = whitening_rel),
      primal_kkt = list(r1_est = r1, r2_est = r2, solve_rel_est = ratio(r1, scale1), stat_rel_est = ratio(r2, scale2),
        residual_ratio_role = "diagnostic_only_under_D23"),
      targets = list(target_error_est = rep(NA_real_, J_est), target_error_cert = rep(NA_real_, J_est),
        loading_error_cert = rep(NA_real_, J_est), se_criterion_est = if (weighting == "optimal") sqrt(colSums(Z1*Z1)) else rep(NA_real_, J_est),
        loading_primal_residual = s1, loading_calibration_residual = s2)))
  }
  m2_g <- gmm_safe_product(m, m, max_matrix_entries, "guard m x m")
  kk_g <- gmm_safe_product(k, k, max_matrix_entries, "guard k x k")
  mk_g <- gmm_safe_product(m, k, max_matrix_entries, "guard m x k")
  gmm_safe_product(12, gmm_safe_product(1, m2_g + kk_g + mk_g, max_matrix_entries, "guard base"), max_matrix_entries, "guard interval workspace")
  if (J_est > 0L) {
    gmm_safe_product(m, J_est, max_matrix_entries, "guard Lambda")
    gmm_safe_product(k, J_est, max_matrix_entries, "guard C_gauge/Mc/Z1")
  }
  if (weighting == "optimal") {
    if (!is.matrix(V_R) || is.object(V_R) || !is.numeric(V_R) || nrow(V_R) != m || ncol(V_R) != m ||
        !is.matrix(R_chol) || is.object(R_chol) || !is.numeric(R_chol) || nrow(R_chol) != m || ncol(R_chol) != m) {
      cdid_gmm_abort("Malformed `V_R` or `R_chol` in `certify_gmm_numerical_guard`.", "cdid_gmm_pd_error")
    }
    if (gmm_sexp_aliases(V_R, R_chol)) {
      cdid_gmm_abort("`V_R` and `R_chol` must not alias in `certify_gmm_numerical_guard`.", "cdid_gmm_bound_error")
    }
  }
  # Preconditions are stored finite matrices; no numerical rank decision is made.
  candidates <- list(W_G, A_whitened, U_qr, mhat_R, theta, y, C_gauge, Lambda, Mc, Z1)
  if (any(vapply(candidates, function(x) !is.numeric(x) || any(!is.finite(x)), logical(1L))))
    cdid_gmm_abort("Nonfinite numerical candidate.", "cdid_gmm_accuracy_error")
  if (!identical(sort(as.integer(piv)), seq_len(k))) cdid_gmm_abort("Invalid QR pivot permutation.", "cdid_gmm_qr_error")
  norm_W <- gmm_iv_norm(W_G, max_matrix_entries = max_matrix_entries); norm_W_F <- norm(W_G, "F")
  if (weighting == "optimal") {
    norm_V <- gmm_iv_norm(V_R, max_matrix_entries = max_matrix_entries); norm_V_F <- norm(V_R, "F")
    eig_V <- certify_symmetric_eigen_bounds(V_R, "V_R", max_matrix_entries = max_matrix_entries)
    chol_res <- gmm_iv_mat_residual(t(R_chol), R_chol, V_R, max_matrix_entries = max_matrix_entries)
    eps_chol_est <- norm(crossprod(R_chol) - V_R, "F")
    eps_chol_cert <- chol_res$norm[2L]
    chol_rel_cert <- gmm_iv_ratio_upper(eps_chol_cert, norm_V[1L])
    rcond_R <- rcond(R_chol, triangular = TRUE)
    sigma_min_R_sq_cert <- gmm_iv_sub(gmm_iv(eig_V$lambda_min_cert), gmm_iv(eps_chol_cert))[1L]
    eta_R_cert <- if (sigma_min_R_sq_cert > 0) gmm_iv_ratio_upper(eps_chol_cert, sigma_min_R_sq_cert) else Inf
    lambda_rel_cert <- if (eig_V$lambda_max_cert > 0) gmm_iv_div(gmm_iv(eig_V$lambda_min_cert), gmm_iv(eig_V$lambda_max_cert))[1L] else -Inf
    pd_reasons <- character(0)
    if (!eig_V$valid || eig_V$delta_Q_cert > accuracy$max_orthogonality_defect) pd_reasons <- c(pd_reasons, "eigenvector_orthogonality_defect_exceeds_bound")
    if (eig_V$lambda_min_cert <= 0) pd_reasons <- c(pd_reasons, "certified_min_eigenvalue_not_positive")
    # eta < 1 is an additional condition needed for the upper A* bound.
    if (sigma_min_R_sq_cert <= 0 || !is.finite(eta_R_cert) || eta_R_cert >= 1) pd_reasons <- c(pd_reasons, "cholesky_factor_residual_exceeds_certified_min_eigenvalue")
    if (lambda_rel_cert < accuracy$min_certified_lambda_rel) pd_reasons <- c(pd_reasons, "certified_relative_eigenvalue_below_caller_threshold")
    if (chol_rel_cert > accuracy$max_factor_residual_rel) pd_reasons <- c(pd_reasons, "cholesky_relative_factor_residual_exceeds_bound")
    if (!is.finite(rcond_R) || rcond_R < accuracy$min_rcond_R) pd_reasons <- c(pd_reasons, "rcond_R_below_caller_threshold")
    if (length(pd_reasons)) cdid_gmm_abort(paste0("Positive-definiteness / Cholesky certification refused: ", paste(pd_reasons, collapse = "; ")),
      "cdid_gmm_pd_error", list(reasons = pd_reasons, lambda_min_cert = eig_V$lambda_min_cert,
        lambda_max_cert = eig_V$lambda_max_cert, eps_chol_est = eps_chol_est, eps_chol_cert = eps_chol_cert,
        sigma_min_R_sq_cert = sigma_min_R_sq_cert, eta_R_cert = eta_R_cert))
    sigma_min_R_cert <- gmm_iv_sqrt(gmm_iv(sigma_min_R_sq_cert))[1L]
    E_A <- gmm_iv_mat_residual(t(R_chol), A_whitened, W_G, max_matrix_entries = max_matrix_entries)
    eps_A_est <- norm(crossprod(R_chol, A_whitened) - W_G, "F"); eps_A_cert <- E_A$norm[2L]
    whitening_rel_cert <- gmm_iv_ratio_upper(eps_A_cert, norm_W[1L])
    if (whitening_rel_cert > accuracy$max_factor_residual_rel) cdid_gmm_abort("Whitened design residual exceeds caller bound.", "cdid_gmm_qr_error")
  } else {
    norm_V <- gmm_iv_sqrt(gmm_iv(m)); norm_V_F <- sqrt(m)
    eig_V <- list(valid = TRUE, lambda_min_est = 1, lambda_max_est = 1, lambda_min_cert = 1, lambda_max_cert = 1,
      delta_Q_est = 0, delta_Q_cert = 0, eps_eig_est = 0, eps_eig_cert = 0, backend = nb_info$backend)
    eps_chol_est <- eps_chol_cert <- chol_rel_cert <- eta_R_cert <- eps_A_est <- eps_A_cert <- whitening_rel_cert <- 0
    rcond_R <- sigma_min_R_cert <- 1
    if (!identical(unname(A_whitened), unname(W_G))) cdid_gmm_abort("Identity criterion requires A = W.", "cdid_gmm_qr_error")
  }
  rcond_U <- rcond(U_qr, triangular = TRUE)
  if (!is.finite(rcond_U) || rcond_U < accuracy$min_rcond_U) cdid_gmm_abort("rcond(U) below caller threshold.", "cdid_gmm_qr_error")
  # Gram is used only to certify singular values, never to solve normal equations.
  AtA_sym <- crossprod(A_whitened)
  if (any(!is.finite(AtA_sym))) cdid_gmm_abort("Gram candidate overflowed.", "cdid_gmm_qr_error")
  for (j in seq_len(k)) for (i in seq_len(j)) AtA_sym[j,i] <- AtA_sym[i,j]
  gram_res <- gmm_iv_mat_residual(t(A_whitened), A_whitened, AtA_sym, max_matrix_entries = max_matrix_entries)
  eps_gram_cert <- gram_res$norm[2L]
  eig_AtA <- certify_symmetric_eigen_bounds(AtA_sym, "A'A", max_matrix_entries = max_matrix_entries)
  if (!eig_AtA$valid || eig_AtA$delta_Q_cert > accuracy$max_orthogonality_defect) cdid_gmm_abort("Gram orthogonality certificate refused.", "cdid_gmm_qr_error")
  amin2 <- gmm_iv_sub(gmm_iv(eig_AtA$lambda_min_cert), gmm_iv(eps_gram_cert))[1L]
  amax2 <- gmm_iv_add(gmm_iv(eig_AtA$lambda_max_cert), gmm_iv(eps_gram_cert))[2L]
  if (amin2 <= 0) cdid_gmm_abort("Certified minimum singular value of A not positive.", "cdid_gmm_qr_error")
  amin <- gmm_iv_sqrt(gmm_iv(amin2))[1L]; amax <- gmm_iv_sqrt(gmm_iv(amax2))[2L]
  if (weighting == "optimal") {
    e <- gmm_iv_div(gmm_iv(eps_A_cert), gmm_iv(sigma_min_R_cert))[2L]
    sigma_min_A_star_cert <- gmm_iv_div(gmm_iv_sub(gmm_iv(amin), gmm_iv(e)),
      gmm_iv_sqrt(gmm_iv_add(gmm_iv(1), gmm_iv(eta_R_cert))))[1L]
    sigma_max_A_star_cert <- gmm_iv_div(gmm_iv_add(gmm_iv(amax), gmm_iv(e)),
      gmm_iv_sqrt(gmm_iv_sub(gmm_iv(1), gmm_iv(eta_R_cert))))[2L]
  } else {
    sigma_min_A_star_cert <- amin; sigma_max_A_star_cert <- amax
  }
  sigma_rel_A_cert <- gmm_iv_div(gmm_iv(sigma_min_A_star_cert), gmm_iv(sigma_max_A_star_cert))[1L]
  sigma_rel_A_est <- sqrt(max(0, eig_AtA$lambda_min_est) / eig_AtA$lambda_max_est)
  if (sigma_min_A_star_cert <= 0 || sigma_rel_A_cert < accuracy$min_sigma_rel_A) cdid_gmm_abort("Exact whitened design certificate refused.", "cdid_gmm_qr_error")
  # Absolute Rayleigh bound: a||x||^2+2sqrt(ad)||x||||y||+d||y||^2 <= (a+d)||z||^2.
  inv_a <- gmm_iv_div(gmm_iv(1), gmm_iv(eig_V$lambda_min_cert))
  inv_s <- gmm_iv_div(gmm_iv(1), gmm_iv(sigma_min_A_star_cert))
  d_inv_sigma_min_A_star_sq_cert <- gmm_iv_mul(inv_s, inv_s)[2L]
  kkt_inv_2norm_cert <- gmm_iv_add(inv_a, gmm_iv_mul(inv_s, inv_s))[2L]
  Vcrit <- if (weighting == "optimal") V_R else diag(m)
  r1_lo <- numeric(m); r1_hi <- numeric(m)
  for (i in seq_len(m)) {
    z <- gmm_iv_sub(gmm_iv_sub(gmm_iv(mhat_R[i]), gmm_iv_dot(Vcrit[i, ], y)), gmm_iv_dot(W_G[i, ], theta))
    r1_lo[i] <- z[1L]; r1_hi[i] <- z[2L]
  }
  r2 <- gmm_iv_residual(t(W_G), y, rep(0,k), max_matrix_entries = max_matrix_entries)
  r2_lo <- -r2$hi; r2_hi <- -r2$lo
  r1_abs_cert <- pmax(abs(r1_lo), abs(r1_hi)); r2_abs_cert <- pmax(abs(r2_lo), abs(r2_hi))
  r1_est <- as.double(mhat_R) - drop(Vcrit %*% y) - drop(W_G %*% theta)
  r2_est <- -drop(crossprod(W_G, y))
  r_norm2_est <- norm(matrix(c(r1_est,r2_est)), "F")
  r_norm2_cert <- gmm_iv_norm(c(r1_lo,r2_lo),c(r1_hi,r2_hi))[2L]
  fit_res <- gmm_iv_residual(W_G, theta, mhat_R, max_matrix_entries = max_matrix_entries)
  fit_residual_norm2_cert <- fit_res$norm[2L]
  sqrt_Q_cert_upper <- if (weighting == "optimal") {
    gmm_iv_div(gmm_iv(fit_residual_norm2_cert), gmm_iv_sqrt(gmm_iv(eig_V$lambda_min_cert)))[2L]
  } else {
    NA_real_
  }
  scale1_iv <- gmm_iv_add(gmm_iv_add(gmm_iv_mul(norm_V, gmm_iv_norm(y)), gmm_iv_mul(norm_W, gmm_iv_norm(theta))), gmm_iv_norm(mhat_R))
  scale2_iv <- gmm_iv_mul(norm_W, gmm_iv_norm(y))
  solve_rel_cert <- gmm_iv_ratio_upper(gmm_iv_norm(r1_lo,r1_hi)[2L], scale1_iv[1L])
  stat_rel_cert <- gmm_iv_ratio_upper(gmm_iv_norm(r2_lo,r2_hi)[2L], scale2_iv[1L])
  solve_rel_est <- gmm_iv_ratio_upper(norm(matrix(r1_est),"F"), scale1_iv[2L])
  stat_rel_est <- gmm_iv_ratio_upper(norm(matrix(r2_est),"F"), scale2_iv[2L])
  kkt_rel_cert <- max(solve_rel_cert, stat_rel_cert); kkt_rel_est <- max(solve_rel_est, stat_rel_est)
  # D12: primal/stationarity relative residual ratios are preserved as diagnostics
  # along with their scales and whether the former diagnostic threshold is exceeded.
  solve_rel_exceeds_threshold <- (!is.finite(solve_rel_cert) || solve_rel_cert > accuracy$max_kkt_residual_rel)
  stat_rel_exceeds_threshold <- (!is.finite(stat_rel_cert) || stat_rel_cert > accuracy$max_kkt_residual_rel)
  kkt_rel_exceeds_threshold <- (!is.finite(kkt_rel_cert) || kkt_rel_cert > accuracy$max_kkt_residual_rel)
  loading_kkt_rel_est <- loading_kkt_rel_cert <- loading_cal_abs_est <- loading_cal_rel_est <- loading_cal_rel_cert <-
    s_norm2_est <- s_norm2_cert <- s1_norm2_cert <- s2_norm2_cert <- G_lambda_cert <-
    c_norm2_cert <- lambda_norm2_cert <- scale_L1_lo <- scale_L1_hi <- scale_L2_vec <-
    t1_cert <- t2_cert <- tau_round_cert <-
    adjoint_missing_term_cert <- target_error_est <- target_error_cert <-
    se_criterion_est <- se_criterion_cert_lower <- se_criterion_cert_upper <-
    target_error_over_se_est <- target_error_over_se_cert <- numeric(J_est)
  loading_kkt_rel_exceeds_threshold <- loading_cal_rel_exceeds_threshold <- loading_rel_exceeds_threshold <- logical(J_est)
  tau_computed <- drop(crossprod(C_gauge, theta))
  for (j in seq_len(J_est)) {
    lam_j <- Lambda[,j]; mc_j <- Mc[,j]; c_j <- C_gauge[,j]
    s1_lo <- numeric(m); s1_hi <- numeric(m)
    for (i in seq_len(m)) {
      z <- gmm_iv_sub(gmm_iv_dot(W_G[i, ], mc_j), gmm_iv_dot(Vcrit[i, ], lam_j))
      s1_lo[i] <- z[1L]; s1_hi[i] <- z[2L]
    }
    s2 <- gmm_iv_residual(t(W_G), lam_j, c_j, max_matrix_entries = max_matrix_entries)
    s2_lo <- -s2$hi; s2_hi <- -s2$lo
    s1_est_j <- drop(W_G %*% mc_j - Vcrit %*% lam_j)
    s2_est_j <- c_j - drop(crossprod(W_G,lam_j))
    s1_norm2_cert[j] <- gmm_iv_norm(s1_lo,s1_hi)[2L]
    s2_norm2_cert[j] <- gmm_iv_norm(s2_lo,s2_hi)[2L]
    s_norm2_cert[j] <- gmm_iv_norm(c(s1_lo,s2_lo),c(s1_hi,s2_hi))[2L]
    s_norm2_est[j] <- norm(matrix(c(s1_est_j,s2_est_j)), "F")
    G_lambda_cert[j] <- gmm_iv_mul(gmm_iv(kkt_inv_2norm_cert), gmm_iv(s_norm2_cert[j]))[2L]
    c_norm2_cert[j] <- gmm_iv_norm(c_j)[2L]
    lambda_norm2_cert[j] <- gmm_iv_norm(lam_j)[2L]
    scale_L1 <- gmm_iv_add(gmm_iv_mul(norm_V, gmm_iv_norm(lam_j)), gmm_iv_mul(norm_W, gmm_iv_norm(mc_j)))
    scale_L2 <- max(gmm_iv_norm(c_j)[1L], gmm_iv_mul(norm_W,gmm_iv_norm(lam_j))[1L], 1)
    scale_L1_lo[j] <- scale_L1[1L]; scale_L1_hi[j] <- scale_L1[2L]; scale_L2_vec[j] <- scale_L2
    loading_kkt_rel_cert[j] <- gmm_iv_ratio_upper(s1_norm2_cert[j],scale_L1[1L])
    loading_cal_rel_cert[j] <- gmm_iv_ratio_upper(s2_norm2_cert[j],scale_L2)
    loading_kkt_rel_est[j] <- gmm_iv_ratio_upper(norm(matrix(s1_est_j),"F"),scale_L1[2L])
    loading_cal_abs_est[j] <- max(abs(s2_est_j)); loading_cal_rel_est[j] <- norm(matrix(s2_est_j),"F") / scale_L2
    loading_kkt_rel_exceeds_threshold[j] <- (!is.finite(loading_kkt_rel_cert[j]) || loading_kkt_rel_cert[j] > accuracy$max_loading_residual_rel)
    loading_cal_rel_exceeds_threshold[j] <- (!is.finite(loading_cal_rel_cert[j]) || loading_cal_rel_cert[j] > accuracy$max_loading_residual_rel)
    loading_rel_exceeds_threshold[j] <- (loading_kkt_rel_exceeds_threshold[j] || loading_cal_rel_exceeds_threshold[j])
    target_error_est[j] <- abs(sum(lam_j*r1_est)) + abs(sum(mc_j*r2_est))
    adjoint_missing_term_cert[j] <- gmm_iv_mul(gmm_iv_mul(gmm_iv(s_norm2_cert[j]),gmm_iv(kkt_inv_2norm_cert)),gmm_iv(r_norm2_cert))[2L]
    t1 <- gmm_iv_dot(abs(lam_j),r1_abs_cert)[2L]; t2 <- gmm_iv_dot(abs(mc_j),r2_abs_cert)[2L]
    # Include rounding in the exported target c'theta, not just the exact linear functional.
    tau_round <- gmm_iv_abs(gmm_iv_sub(gmm_iv(tau_computed[j]),gmm_iv_dot(c_j,theta)))[2L]
    t1_cert[j] <- t1; t2_cert[j] <- t2; tau_round_cert[j] <- tau_round
    target_error_cert[j] <- gmm_iv_add(gmm_iv_add(gmm_iv_add(gmm_iv(t1),gmm_iv(t2)),gmm_iv(adjoint_missing_term_cert[j])),gmm_iv(tau_round))[2L]
    if (weighting == "optimal") {
      se_criterion_est[j] <- norm(matrix(Z1[,j]),"F")
      vpen <- gmm_iv_mul(gmm_iv_mul(gmm_iv(c_norm2_cert[j]),gmm_iv(kkt_inv_2norm_cert)),gmm_iv(s_norm2_cert[j]))
      cmc_iv <- gmm_iv_dot(c_j,mc_j)
      vlower <- gmm_iv_sub(cmc_iv,vpen)[1L]
      se_criterion_cert_lower[j] <- if (vlower > 0) gmm_iv_sqrt(gmm_iv(vlower))[1L] else 0
      vupper_adj <- max(0, gmm_iv_add(cmc_iv, vpen)[2L])
      vupper_d <- gmm_iv_mul(gmm_iv(d_inv_sigma_min_A_star_sq_cert), gmm_iv_mul(gmm_iv(c_norm2_cert[j]), gmm_iv(c_norm2_cert[j])))[2L]
      se_criterion_cert_upper[j] <- gmm_iv_sqrt(gmm_iv(min(vupper_adj, vupper_d)))[2L]
      target_error_over_se_est[j] <- if (se_criterion_est[j] > 0) target_error_est[j] / se_criterion_est[j] else Inf
      target_error_over_se_cert[j] <- gmm_iv_ratio_upper(target_error_cert[j],se_criterion_cert_lower[j])
    } else {
      se_criterion_est[j] <- se_criterion_cert_lower[j] <- se_criterion_cert_upper[j] <- target_error_over_se_est[j] <- target_error_over_se_cert[j] <- NA_real_
    }
  }
  loading_budget_supplied <- !is.null(accuracy$max_loading_error_abs)
  if (J_est > 0) {
    if (any(!is.finite(target_error_cert)) || any(!is.finite(G_lambda_cert)))
      cdid_gmm_abort("Undefined or nonfinite target or loading forward-error bound.", "cdid_gmm_accuracy_error",
        list(reasons = "nonfinite_forward_error_bound", target_error_cert = target_error_cert, G_lambda = G_lambda_cert))
    if (accuracy$target_error_unit %in% c("abs","both") && any(target_error_cert > accuracy$max_target_error_abs))
      cdid_gmm_abort("Certified target forward-error bound exceeds caller budget.", "cdid_gmm_accuracy_error",
        list(reasons = "certified_target_error_abs_exceeded", target_error_est = target_error_est,
             target_error_cert = target_error_cert, adjoint_missing_term_cert = adjoint_missing_term_cert,
             tau_round_cert = tau_round_cert, G_lambda = G_lambda_cert))
    if (accuracy$target_error_unit %in% c("model_se","both") && any(!is.finite(target_error_over_se_cert) | target_error_over_se_cert > accuracy$max_target_error_se))
      cdid_gmm_abort("Certified relative target forward-error bound exceeds caller budget.", "cdid_gmm_accuracy_error",
        list(reasons = "certified_target_error_se_exceeded", target_error_over_se_cert = target_error_over_se_cert,
             se_criterion_cert_lower = se_criterion_cert_lower))
    if (loading_budget_supplied && any(G_lambda_cert > accuracy$max_loading_error_abs))
      cdid_gmm_abort("Certified target loading forward-error bound exceeds caller budget.", "cdid_gmm_accuracy_error",
        list(reasons = "certified_loading_error_abs_exceeded", G_lambda = G_lambda_cert,
             loading_error_cert = G_lambda_cert, s_norm2_cert = s_norm2_cert,
             s1_norm2_cert = s1_norm2_cert, s2_norm2_cert = s2_norm2_cert,
             kkt_inv_2norm_cert = kkt_inv_2norm_cert, max_loading_error_abs = accuracy$max_loading_error_abs))
  }
  point_certified_vec <- rep(TRUE, J_est)
  loading_certified_vec <- rep(loading_budget_supplied, J_est)
  loading_certified <- loading_budget_supplied
  loading_validation_status <- if (loading_budget_supplied) {
    "certified_conditional_stored_array"
  } else {
    "incomplete_missing_caller_budget"
  }
  list(
    certified = TRUE,
    point_certified = TRUE,
    point_validation_status = "certified_conditional_stored_array",
    loading_budget_supplied = loading_budget_supplied,
    loading_certified = loading_certified,
    loading_validation_status = loading_validation_status,
    complete_inference_certified = FALSE,
    full_inference_status = "pending_unvalidated_raw_root_covariance_projection_multiplier",
    backend = nb_info$backend,
    native_backend = nb_info,
    arithmetic_assumptions = "IEEE binary64 round-to-nearest ties-to-even scalar operations and sqrt; gradual underflow; no FTZ/DAZ/reassociation. BLAS/LAPACK output is an untrusted candidate.",
    sampling_provenance_certified = FALSE,
    computed_target_values = tau_computed,
    method = accuracy$method,
    weighting = weighting,
    eigen_V = eig_V,
    cholesky = list(
      eps_chol_est = eps_chol_est,
      eps_chol_cert = eps_chol_cert,
      chol_rel_cert = chol_rel_cert,
      sigma_min_R_cert = sigma_min_R_cert,
      eta_R_cert = eta_R_cert,
      rcond_R_1norm = rcond_R
    ),
    whitening_and_qr = list(
      eps_A_est = eps_A_est,
      eps_A_cert = eps_A_cert,
      whitening_rel_cert = whitening_rel_cert,
      rcond_U_1norm = rcond_U,
      rcond_U_norm = "1-norm reciprocal condition estimate (LAPACK dtrcon)",
      sigma_rel_A_est_2norm = sigma_rel_A_est,
      sigma_min_A_star_cert = sigma_min_A_star_cert,
      sigma_max_A_star_cert = sigma_max_A_star_cert,
      sigma_rel_A_cert_2norm = sigma_rel_A_cert,
      d_inv_sigma_min_A_star_sq_cert = d_inv_sigma_min_A_star_sq_cert,
      kkt_inv_2norm_cert = kkt_inv_2norm_cert
    ),
    primal_kkt = list(
      r1_est = r1_est,
      r1_abs_cert = r1_abs_cert,
      r2_est = r2_est,
      r2_abs_cert = r2_abs_cert,
      r_norm2_est = r_norm2_est,
      r_norm2_cert = r_norm2_cert,
      scale1_iv = scale1_iv,
      scale2_iv = scale2_iv,
      solve_rel_est = solve_rel_est,
      solve_rel_cert = solve_rel_cert,
      stat_rel_est = stat_rel_est,
      stat_rel_cert = stat_rel_cert,
      kkt_rel_est = kkt_rel_est,
      kkt_rel_cert = kkt_rel_cert,
      max_kkt_residual_rel_threshold = accuracy$max_kkt_residual_rel,
      solve_rel_exceeds_threshold = solve_rel_exceeds_threshold,
      stat_rel_exceeds_threshold = stat_rel_exceeds_threshold,
      kkt_rel_exceeds_threshold = kkt_rel_exceeds_threshold,
      residual_ratio_role = "diagnostic_only_under_d12"
    ),
    targets = list(
      criterion_scale_only = TRUE,
      loading_kkt_rel_est = loading_kkt_rel_est,
      loading_kkt_rel_cert = loading_kkt_rel_cert,
      loading_cal_abs_est = loading_cal_abs_est,
      loading_cal_rel_est = loading_cal_rel_est,
      loading_cal_rel_cert = loading_cal_rel_cert,
      s1_norm2_cert = s1_norm2_cert,
      s2_norm2_cert = s2_norm2_cert,
      s_norm2_est = s_norm2_est,
      s_norm2_cert = s_norm2_cert,
      scale_L1_lo = scale_L1_lo,
      scale_L1_hi = scale_L1_hi,
      scale_L2 = scale_L2_vec,
      max_loading_residual_rel_threshold = accuracy$max_loading_residual_rel,
      loading_kkt_rel_exceeds_threshold = loading_kkt_rel_exceeds_threshold,
      loading_cal_rel_exceeds_threshold = loading_cal_rel_exceeds_threshold,
      loading_rel_exceeds_threshold = loading_rel_exceeds_threshold,
      G_lambda = G_lambda_cert,
      loading_error_cert = G_lambda_cert,
      max_loading_error_abs = accuracy$max_loading_error_abs,
      loading_budget_supplied = loading_budget_supplied,
      point_certified = point_certified_vec,
      loading_certified = loading_certified_vec,
      complete_inference_certified = rep(FALSE, J_est),
      direct_primal_term_cert = t1_cert,
      direct_stationarity_term_cert = t2_cert,
      tau_round_cert = tau_round_cert,
      target_round_cert = tau_round_cert,
      adjoint_missing_term_cert = adjoint_missing_term_cert,
      target_error_est = target_error_est,
      target_error_cert = target_error_cert,
      se_criterion_est = se_criterion_est,
      se_criterion_cert_lower = se_criterion_cert_lower,
      se_criterion_cert_upper = se_criterion_cert_upper,
      target_error_over_se_est = target_error_over_se_est,
      target_error_over_se_cert = target_error_over_se_cert,
      bounds_distinction = paste0(
        "`target_error_est` is an ordinary computed indicator (|Lambda'r1| + |Mc'r2|); ",
        "`target_error_cert` is a conservative sufficient forward-error bound including ",
        "outward-rounded scalar interval residual bounds, exported target rounding, and ||s||_2 ||K^{-1}||_2 ||r||_2; ",
        "`G_lambda` is the outward loading forward-error bound upper(||K^{-1}||_2) upper(||s||_2) including both adjoint blocks. ",
        "Criterion SE is not actual sandwich/inference SE."
      )
    ),
    propagation_bounds = list(
      l_lambda_min_V_cert = eig_V$lambda_min_cert,
      lambda_max_V_cert = eig_V$lambda_max_cert,
      a_inv_lambda_min_V_cert = inv_a[2L],
      d_inv_sigma_min_A_star_sq_cert = d_inv_sigma_min_A_star_sq_cert,
      kkt_inv_2norm_cert = kkt_inv_2norm_cert,
      r_norm2_cert = r_norm2_cert,
      fit_residual_norm2_cert = fit_residual_norm2_cert,
      sqrt_Q_cert_upper = sqrt_Q_cert_upper,
      G_c = target_error_cert,
      G_lambda = G_lambda_cert,
      s_norm2_cert = s_norm2_cert,
      s1_norm2_cert = s1_norm2_cert,
      s2_norm2_cert = s2_norm2_cert,
      c_norm2_cert = c_norm2_cert,
      lambda_norm2_cert = lambda_norm2_cert,
      se_criterion_cert_lower = se_criterion_cert_lower,
      se_criterion_cert_upper = se_criterion_cert_upper,
      raw_root_transfer_implemented = FALSE,
      covariance_transfer_implemented = FALSE,
      projection_evaluation_bound_implemented = FALSE,
      sandwich_accumulation_bound_implemented = FALSE,
      multiplier_transfer_implemented = FALSE,
      complete_inference_certified = FALSE,
      status = "pending_raw_root_covariance_and_full_inference_transfer"
    ),
    restriction_classification = list(
      theorem_prerequisites = c(
        "delta_Q_cert < 1",
        "lambda_min_cert > 0",
        "sigma_min_R_sq_cert = lambda_min_cert - eps_chol_cert > 0",
        "eta_R_cert = eps_chol_cert / sigma_min_R_sq_cert < 1",
        "sigma_min_A_star_cert > 0",
        "finite nonzero triangular pivots in R_chol and U_qr",
        "finite outward interval enclosures"
      ),
      additional_caller_restrictions = c(
        "max_orthogonality_defect < 0.5",
        "min_certified_lambda_rel",
        "max_factor_residual_rel",
        "min_rcond_R",
        "min_rcond_U",
        "min_sigma_rel_A"
      ),
      diagnostic_only_ratios = c(
        "solve_rel_cert",
        "stat_rel_cert",
        "kkt_rel_cert",
        "loading_kkt_rel_cert",
        "loading_cal_rel_cert"
      ),
      conditional_stored_array_budgets = c(
        "max_target_error_abs",
        "max_target_error_se (conditional criterion scale only, not sandwich/inference SE)",
        "max_loading_error_abs (optional; NULL means not supplied / incomplete loading certification)"
      )
    )
  )
}

# ---------------------------------------------------------------------------
# Streamed IF replay specification and two-pass projected cluster accumulation
# ---------------------------------------------------------------------------

gmm_if_replay_spec <- function(
  callback,
  input_row_level,
  aggregate_to_clusters,
  n_input_rows,
  n_output_rows,
  n_positive_contributing_rows,
  centered,
  weighted,
  weight_type,
  normalization_factor,
  scale_id,
  nuisance_fit_id,
  expected_unit_ids = NULL,
  expected_cluster_ids = NULL
) {
  if (missing(callback) || missing(input_row_level) || missing(aggregate_to_clusters) ||
      missing(n_input_rows) || missing(n_output_rows) || missing(n_positive_contributing_rows) ||
      missing(centered) || missing(weighted) || missing(weight_type) ||
      missing(normalization_factor) || missing(scale_id) || missing(nuisance_fit_id)) {
    cdid_gmm_abort(
      "All fields of `gmm_if_replay_spec` must be supplied explicitly by the caller.",
      c("cdid_gmm_replay_error", "cdid_gmm_control_error")
    )
  }
  if (!is.function(callback)) {
    cdid_gmm_abort("`if_replay$callback` must be a function `function(consume_batch)`.", "cdid_gmm_replay_error")
  }
  if (!gmm_scalar_character(input_row_level) || !input_row_level %in% c("unit", "cluster")) {
    cdid_gmm_abort("`input_row_level` must be 'unit' or 'cluster'.", "cdid_gmm_replay_error")
  }
  if (!gmm_scalar_logical(aggregate_to_clusters)) {
    cdid_gmm_abort("`aggregate_to_clusters` must be a scalar logical.", "cdid_gmm_replay_error")
  }
  if (input_row_level == "cluster" && isTRUE(aggregate_to_clusters)) {
    cdid_gmm_abort(
      "IF replay already supplies complete cluster sums (`input_row_level = 'cluster'`); never sum clusters twice.",
      "cdid_gmm_replay_error"
    )
  }
  n_in <- gmm_int_bound(n_input_rows, "n_input_rows")
  n_out <- gmm_int_bound(n_output_rows, "n_output_rows")
  n_pos <- gmm_int_bound(n_positive_contributing_rows, "n_positive_contributing_rows")
  if (n_out > n_in || n_pos > n_out) {
    cdid_gmm_abort("Invalid `n_input_rows`, `n_output_rows`, or `n_positive_contributing_rows` in `if_replay`.", "cdid_gmm_replay_error")
  }
  if (!isTRUE(aggregate_to_clusters) && n_out != n_in) {
    cdid_gmm_abort("When `aggregate_to_clusters` is FALSE, `n_output_rows` must equal `n_input_rows`.", "cdid_gmm_replay_error")
  }
  if (!gmm_scalar_logical(centered) || !gmm_scalar_logical(weighted)) {
    cdid_gmm_abort("`centered` and `weighted` must be scalar logicals in `if_replay`.", "cdid_gmm_replay_error")
  }
  if (!gmm_scalar_character(weight_type) || !gmm_scalar_character(scale_id) || !gmm_scalar_character(nuisance_fit_id)) {
    cdid_gmm_abort("`weight_type`, `scale_id`, and `nuisance_fit_id` must be non-empty character scalars in `if_replay`.", "cdid_gmm_replay_error")
  }
  if (!gmm_scalar_numeric(normalization_factor, positive = TRUE)) {
    cdid_gmm_abort("`normalization_factor` must be a finite positive numeric scalar in `if_replay`.", "cdid_gmm_replay_error")
  }

  structure(
    list(
      callback = callback,
      input_row_level = input_row_level,
      aggregate_to_clusters = aggregate_to_clusters,
      output_row_level = if (input_row_level == "cluster" || isTRUE(aggregate_to_clusters)) "cluster" else "unit",
      n_input_rows = n_in,
      n_output_rows = n_out,
      n_positive_contributing_rows = n_pos,
      centered = centered,
      weighted = weighted,
      weight_type = weight_type,
      normalization_factor = normalization_factor,
      scale_id = scale_id,
      nuisance_fit_id = nuisance_fit_id,
      expected_unit_ids = expected_unit_ids,
      expected_cluster_ids = expected_cluster_ids
    ),
    class = "cdid_gmm_if_replay_spec"
  )
}

validate_gmm_if_replay_spec <- function(if_replay) {
  if (missing(if_replay) || is.null(if_replay) || !inherits(if_replay, "cdid_gmm_if_replay_spec")) {
    cdid_gmm_abort("`if_replay` must be constructed explicitly via `gmm_if_replay_spec(...)`.", "cdid_gmm_replay_error")
  }
  gmm_if_replay_spec(
    callback = if_replay$callback,
    input_row_level = if_replay$input_row_level,
    aggregate_to_clusters = if_replay$aggregate_to_clusters,
    n_input_rows = if_replay$n_input_rows,
    n_output_rows = if_replay$n_output_rows,
    n_positive_contributing_rows = if_replay$n_positive_contributing_rows,
    centered = if_replay$centered,
    weighted = if_replay$weighted,
    weight_type = if_replay$weight_type,
    normalization_factor = if_replay$normalization_factor,
    scale_id = if_replay$scale_id,
    nuisance_fit_id = if_replay$nuisance_fit_id,
    expected_unit_ids = if_replay$expected_unit_ids,
    expected_cluster_ids = if_replay$expected_cluster_ids
  )
}

gmm_if_unit_batch_eligible <- function(if_replay, J = 1L) {
  inherits(if_replay, "cdid_gmm_if_replay_spec") &&
    identical(if_replay$input_row_level, "unit") &&
    !isTRUE(if_replay$aggregate_to_clusters) &&
    identical(if_replay$output_row_level, "unit") &&
    !isTRUE(if_replay$weighted) &&
    is.numeric(J) && length(J) == 1L && !is.na(J) && J >= 1L
}

gmm_require_unit_batch_symbol <- function(syms = NULL) {
  # Preserve the old nine-symbol generic resolver/cache. This extra kernel is
  # authenticated only on eligible unit replays, against that resolver's DLL.
  checked_syms <- gmm_require_native_dll()
  owner <- environment(gmm_require_native_dll)
  dlls <- if (isNamespace(owner)) getNamespaceInfo(owner, "DLLs") else base::getLoadedDLLs()
  dll <- dlls[["cdid"]]
  si <- tryCatch(base::getNativeSymbolInfo("C_gmm_if_unit_batch", PACKAGE = dll), error = function(e) NULL)
  sym <- if (is.null(si)) NULL else si$address
  refuse <- function() cdid_gmm_abort(
    "Native symbol registration/forcing verification failed for `C_gmm_if_unit_batch`.",
    c("cdid_gmm_native_error", "cdid_gmm_accuracy_error")
  )
  if (is.null(si) || !inherits(si, "NativeSymbolInfo") || is.null(sym) ||
      typeof(sym) != "externalptr" || !inherits(sym, "NativeSymbol") ||
      !identical(as.integer(si$numParameters), 10L) ||
      any(vapply(checked_syms, identical, logical(1L), sym))) refuse()
  if (!is.null(syms)) {
    supplied <- if (identical(syms, checked_syms)) sym else if (
      typeof(syms) == "externalptr" && inherits(syms, "NativeSymbol")) syms else if (
      is.list(syms)) syms[["C_gmm_if_unit_batch"]] else NULL
    other <- if (is.list(syms)) syms[names(syms) != "C_gmm_if_unit_batch"] else list()
    if (is.null(supplied) || typeof(supplied) != "externalptr" ||
        !inherits(supplied, "NativeSymbol") || !identical(supplied, sym) ||
        any(vapply(other, identical, logical(1L), sym))) refuse()
  }
  sym
}

gmm_if_unit_batch <- function(
  pass, psi, w, flags, counts,
  weighted_sum = NULL, raw_sum = NULL, psi_mean = NULL, S = NULL,
  max_matrix_entries = NULL, syms = NULL
) {
  sym <- gmm_require_unit_batch_symbol(syms)
  max_e <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
  if (!is.matrix(psi) || is.object(psi) || typeof(psi) != "double") {
    cdid_gmm_abort("`psi` in `C_gmm_if_unit_batch` must be a non-object double matrix.", "cdid_gmm_accuracy_error")
  }
  B <- nrow(psi); J <- ncol(psi)
  if (is.na(B) || is.na(J) || B < 1L || J < 1L) {
    cdid_gmm_abort("`psi` in `C_gmm_if_unit_batch` must have positive finite matrix dimensions.", "cdid_gmm_bound_error")
  }
  bj <- gmm_safe_product(B, J, max_e, "unit batch psi (B x J)")
  b2 <- gmm_safe_product(2, B, max_e, "unit batch weights/flags (2 x B)")
  if (identical(pass, 1L)) {
    j4 <- gmm_safe_product(4, J, max_e, "unit batch pass-1 sums (4 x J)")
    gmm_safe_product(1, bj + b2 + j4 + 8, max_e, "unit batch pass-1 simultaneous buffers")
  } else if (identical(pass, 2L)) {
    j2 <- gmm_safe_product(2, J, max_e, "unit batch pass-2 centering (2 x J)")
    jj <- gmm_safe_product(J, J, max_e, "unit batch covariance (J x J)")
    jj2 <- gmm_safe_product(2, jj, max_e, "unit batch pass-2 covariance buffers (2 x J x J)")
    gmm_safe_product(1, bj + b2 + j2 + jj2 + 8, max_e, "unit batch pass-2 simultaneous buffers")
  }
  .Call(sym, pass, psi, w, flags, counts, weighted_sum, raw_sum, psi_mean, S, max_e)
}

project_gmm_if_stream <- function(if_replay, Lambda, m_R, bounds, provenance, materialize_target_if = FALSE) {
  if_replay <- validate_gmm_if_replay_spec(if_replay); bounds <- validate_gmm_bounds_control(bounds)
  J <- ncol(Lambda); ni <- if_replay$n_input_rows; no <- if_replay$n_output_rows
  if (!is.matrix(Lambda) || nrow(Lambda) != m_R || any(!is.finite(Lambda))) cdid_gmm_abort("Invalid projection loadings.", "cdid_gmm_replay_error")
  gmm_safe_product(ni,J,bounds$max_matrix_entries,"projected replay cache")
  gmm_safe_product(8,ni+no,bounds$max_matrix_entries,"replay identifiers/weights")
  gmm_safe_product(bounds$max_batch_rows,m_R,bounds$max_matrix_entries,"replay batch")
  gmm_safe_product(bounds$max_batch_rows,J,bounds$max_matrix_entries,"projected batch")
  gmm_safe_product(J,J,bounds$max_matrix_entries,"sandwich accumulator")
  if (materialize_target_if) {
    if (!bounds$allow_materialize_if) cdid_gmm_abort("IF materialization is not authorized by bounds.", "cdid_gmm_bound_error")
    gmm_safe_product(no,J,bounds$max_materialized_if_entries,"materialized IF")
    gmm_safe_product(no,J,bounds$max_matrix_entries,"materialized IF")
  }
  if (!identical(if_replay$nuisance_fit_id,provenance$nuisance_fit_id)) cdid_gmm_abort("Replay nuisance label differs.", "cdid_gmm_replay_error")
  keys <- function(x, n, label) {
    if (!is.null(dim(x)) || is.object(x) || length(x) != n || anyNA(x) ||
        !(is.numeric(x) || is.character(x)) || (is.numeric(x) && any(!is.finite(x))) ||
        (is.character(x) && any(!nzchar(x)))) cdid_gmm_abort(paste("Invalid",label), "cdid_gmm_replay_error")
    if (is.numeric(x)) sprintf("%.17g",as.double(x)) else x
  }
  if (ni > .Machine$integer.max || no > .Machine$integer.max) cdid_gmm_abort("Replay dimensions exceed supported R integer dimensions.", "cdid_gmm_bound_error")
  unit_batch <- gmm_if_unit_batch_eligible(if_replay, J)
  unit_batch_sym <- NULL
  if (unit_batch) {
    B_max <- min(ni, bounds$max_batch_rows)
    bj_max <- gmm_safe_product(B_max, J, bounds$max_matrix_entries, "unit batch projected rows")
    b2_max <- gmm_safe_product(2, B_max, bounds$max_matrix_entries, "unit batch weights/flags")
    j4_max <- gmm_safe_product(4, J, bounds$max_matrix_entries, "unit batch pass-1 sums")
    j2_max <- gmm_safe_product(2, J, bounds$max_matrix_entries, "unit batch pass-2 centering")
    jj2_max <- gmm_safe_product(2, gmm_safe_product(J, J, bounds$max_matrix_entries, "sandwich accumulator"),
                                bounds$max_matrix_entries, "unit batch pass-2 covariance buffers")
    gmm_safe_product(1, bj_max + b2_max + j4_max + 8, bounds$max_matrix_entries, "unit batch pass-1 simultaneous buffers")
    gmm_safe_product(1, bj_max + b2_max + j2_max + jj2_max + 8, bounds$max_matrix_entries, "unit batch pass-2 simultaneous buffers")
    unit_batch_sym <- gmm_require_unit_batch_symbol(gmm_require_native_dll())
  }
  unit_cache <- cluster_cache <- character(ni); weight_cache <- numeric(ni)
  psi_cache <- matrix(0,ni,J)
  output_keys <- character(no); output_weights <- numeric(no)
  batch_sizes <- numeric(ni); sum_w <- 0; weighted_sum <- raw_sum <- numeric(J)
  npos <- 0L; nb1 <- 0L; max_B <- 0L
  closed <- character(no); nclosed <- 0L
  run <- function(pass) {
    nin <- nout <- nb <- 0L; open_id <- NULL; open_w <- 0; open_vec <- numeric(J)
    # `closed_before`, when supplied, is the caller's exact answer to
    # `id %in% closed[seq_len(nclosed)]` for this flush.
    flush <- function(id,w,vec,closed_before=NULL) {
      nout <<- nout+1L
      if (nout > no) cdid_gmm_abort("Replay exceeds declared output count.", "cdid_gmm_replay_error")
      if (any(!is.finite(vec))) cdid_gmm_abort("Cluster sum overflowed.", "cdid_gmm_replay_error")
      if (pass == 1L) {
        if (if (is.null(closed_before)) id %in% closed[seq_len(nclosed)] else closed_before)
          cdid_gmm_abort("Closed cluster reappeared.", "cdid_gmm_replay_error")
        nclosed <<- nclosed+1L; closed[nclosed] <<- id
        output_keys[nout] <<- id; output_weights[nout] <<- w
        if (w > 0) npos <<- npos+1L
        sum_w <<- sum_w+w; weighted_sum <<- weighted_sum+w*vec; raw_sum <<- raw_sum+vec
        if (!is.finite(sum_w) || any(!is.finite(weighted_sum)) || any(!is.finite(raw_sum)))
          cdid_gmm_abort("Replay mean/weight accumulation overflowed.", "cdid_gmm_replay_error")
      } else {
        if (!identical(id,output_keys[nout]) || !identical(w,output_weights[nout])) cdid_gmm_abort("Replay output metadata changed.", "cdid_gmm_replay_error")
        z <- vec-psi_mean
        S <<- S+w*tcrossprod(z)
        if (any(!is.finite(S))) cdid_gmm_abort("Replay covariance overflowed.", "cdid_gmm_replay_error")
        if (!is.null(materialized)) materialized[nout,] <<- vec
      }
    }
    consume <- function(batch) {
      if (!is.list(batch) || anyDuplicated(names(batch)) || !all(c("unit_id","phi") %in% names(batch))) cdid_gmm_abort("Invalid IF batch.", "cdid_gmm_replay_error")
      phi <- batch$phi
      if (!is.matrix(phi) || !is.numeric(phi) || ncol(phi) != m_R || nrow(phi) < 1L) cdid_gmm_abort("Invalid IF batch shape.", "cdid_gmm_replay_error")
      B <- nrow(phi)
      if (B > bounds$max_batch_rows) cdid_gmm_abort("Replay batch exceeds row bound.", c("cdid_gmm_bound_error","cdid_gmm_replay_error"))
      if (nin+B > ni) cdid_gmm_abort("Replay exceeds declared input count.", "cdid_gmm_replay_error")
      gmm_safe_product(B,m_R,bounds$max_matrix_entries,"IF finite/projection validation")
      if (any(!is.finite(phi))) cdid_gmm_abort("Nonfinite IF batch.", "cdid_gmm_replay_error")
      if (!is.null(batch$nuisance_fit_id) && !identical(batch$nuisance_fit_id,if_replay$nuisance_fit_id)) cdid_gmm_abort("Changed batch nuisance label.", "cdid_gmm_replay_error")
      u <- keys(batch$unit_id,B,"unit IDs")
      cid <- if (if_replay$aggregate_to_clusters) keys(batch$cluster_id,B,"cluster IDs") else u
      if (if_replay$input_row_level == "cluster" && !is.null(batch$cluster_id)) {
        if (!identical(keys(batch$cluster_id,B,"cluster IDs"),u)) cdid_gmm_abort("Complete cluster IDs disagree with row IDs.", "cdid_gmm_replay_error")
      }
      if (if_replay$weighted) {
        w <- batch$weight
        if (!is.numeric(w) || !is.null(dim(w)) || is.object(w) || length(w)!=B || any(!is.finite(w)) || any(w<0)) cdid_gmm_abort("Invalid replay weights.", "cdid_gmm_replay_error")
        w <- as.double(w)
      } else {
        if (!is.null(batch$weight) && (!is.numeric(batch$weight) || length(batch$weight)!=B || anyNA(batch$weight) || any(batch$weight!=1))) cdid_gmm_abort("Invalid unweighted replay weights.", "cdid_gmm_replay_error")
        w <- rep(1,B)
      }
      psi <- phi %*% Lambda
      if (any(!is.finite(psi))) cdid_gmm_abort("Projected IF overflowed.", "cdid_gmm_replay_error")
      idx <- seq.int(nin+1L,nin+B); nb <<- nb+1L
      if (pass==1L) {
        unit_cache[idx] <<- u; cluster_cache[idx] <<- cid; weight_cache[idx] <<- w; psi_cache[idx,] <<- psi
        batch_sizes[nb] <<- B; max_B <<- max(max_B,B)
      } else {
        if (nb > nb1 || B != batch_sizes[nb] || !identical(u,unit_cache[idx]) || !identical(cid,cluster_cache[idx]) ||
            !identical(w,weight_cache[idx]) || !identical(as.double(psi),as.double(psi_cache[idx,,drop=FALSE])))
          cdid_gmm_abort("IF replay changed projected rows, weights, IDs or batch structure.", "cdid_gmm_replay_error")
      }
      nin <<- nin+B
      # Unit rows are flushed, and closed, in this order, so row i's membership
      # among the keys closed before it is its duplicated() flag in that same
      # sequence: one hashed pass per batch replaces a scan of every closed key
      # per row. This is O(N+B) per batch and O(N^2/B) overall, not linear.
      # In this unit branch no=ni and nclosed+B<=ni. Retained ID/weight arrays
      # use 7*ni slots; combined keys, hash/logical scratch and batch keys/weights
      # fit the existing 8*(ni+no)=16*ni identifier reservation. The hash/input
      # are released before the indexed flag result and before row accumulation.
      # Same check, order and refusal; floating-point additions are unchanged.
      closed_flags <- if (!if_replay$aggregate_to_clusters && pass == 1L)
        duplicated(c(closed[seq_len(nclosed)], u))[nclosed + seq_len(B)] else NULL
      if (unit_batch) {
        counts_vec <- c(as.double(nout), as.double(no), as.double(npos), sum_w)
        if (pass == 1L) {
          res_b <- gmm_if_unit_batch(
            1L, psi, w, closed_flags, counts_vec,
            weighted_sum = weighted_sum, raw_sum = raw_sum,
            max_matrix_entries = bounds$max_matrix_entries, syms = unit_batch_sym
          )
          out_idx <- seq.int(nout + 1L, nout + B)
          closed[seq.int(nclosed + 1L, nclosed + B)] <<- u
          nclosed <<- nclosed + B
          output_keys[out_idx] <<- u
          output_weights[out_idx] <<- w
          nout <<- as.integer(res_b$counts[1L])
          npos <<- as.integer(res_b$counts[3L])
          sum_w <<- res_b$counts[4L]
          weighted_sum <<- res_b$weighted_sum
          raw_sum <<- res_b$raw_sum
        } else {
          meta_flags <- if (nout + B <= no) {
            out_idx <- seq.int(nout + 1L, nout + B)
            (u != output_keys[out_idx]) | (w != output_weights[out_idx])
          } else {
            logical(B)
          }
          res_b <- gmm_if_unit_batch(
            2L, psi, w, meta_flags, counts_vec,
            psi_mean = psi_mean, S = S,
            max_matrix_entries = bounds$max_matrix_entries, syms = unit_batch_sym
          )
          if (!is.null(materialized)) materialized[seq.int(nout + 1L, nout + B), ] <<- psi
          nout <<- as.integer(res_b$counts[1L])
          S <<- res_b$S
        }
        return(invisible(NULL))
      }
      for (i in seq_len(B)) {
        if (!if_replay$aggregate_to_clusters) { flush(u[i],w[i],psi[i,],closed_flags[i]); next }
        if (is.null(open_id)) {
          open_id <<- cid[i]; open_w <<- w[i]; open_vec <<- psi[i,]
        } else if (cid[i]==open_id) {
          if (w[i]!=open_w) cdid_gmm_abort("Within-cluster weight changed.", "cdid_gmm_replay_error")
          open_vec <<- open_vec+psi[i,]
        } else {
          flush(open_id,open_w,open_vec)
          open_id <<- cid[i]; open_w <<- w[i]; open_vec <<- psi[i,]
        }
      }
    }
    tryCatch(withCallingHandlers(if_replay$callback(consume),
      warning = function(w) cdid_gmm_abort(paste0("IF replay warning: ", conditionMessage(w)), "cdid_gmm_replay_error")),
      error = function(e) {
        if (inherits(e, "cdid_gmm_refusal")) stop(e)
        cdid_gmm_abort(paste0("IF replay callback failed: ", conditionMessage(e)), "cdid_gmm_replay_error")
      })
    if (!is.null(open_id)) flush(open_id,open_w,open_vec)
    if (nin!=ni || nout!=no || nb<1L || (pass==2L && nb!=nb1)) cdid_gmm_abort("Replay counts/batches differ.", "cdid_gmm_replay_error")
    nb
  }
  nb1 <- run(1L)
  if (anyDuplicated(unit_cache)) cdid_gmm_abort("Duplicate replay unit IDs.", "cdid_gmm_replay_error")
  if (npos!=if_replay$n_positive_contributing_rows) cdid_gmm_abort("Positive-weight count differs.", "cdid_gmm_replay_error")
  if (!is.null(if_replay$expected_unit_ids) && !identical(unit_cache,keys(if_replay$expected_unit_ids,ni,"expected unit IDs"))) cdid_gmm_abort("Expected unit sequence differs.", "cdid_gmm_replay_error")
  if (!is.null(if_replay$expected_cluster_ids) && !identical(output_keys,keys(if_replay$expected_cluster_ids,no,"expected cluster IDs"))) cdid_gmm_abort("Expected cluster sequence differs.", "cdid_gmm_replay_error")
  if (sum_w <= 0) cdid_gmm_abort("No positive replay weight.", "cdid_gmm_replay_error")
  psi_mean <- if (if_replay$centered) weighted_sum/sum_w else numeric(J)
  if (any(!is.finite(psi_mean))) cdid_gmm_abort("Replay mean overflowed.", "cdid_gmm_replay_error")
  S <- matrix(0,J,J,dimnames=list(colnames(Lambda),colnames(Lambda)))
  materialized <- if (materialize_target_if) matrix(0,no,J,dimnames=list(output_keys,colnames(Lambda))) else NULL
  run(2L)
  sandwich <- if_replay$normalization_factor*S
  if (any(!is.finite(sandwich)) || any(diag(sandwich)<0)) cdid_gmm_abort("Nonfinite/negative sandwich variance.", "cdid_gmm_replay_error")
  list(sandwich_cov=sandwich,se_sandwich=sqrt(diag(sandwich)),psi_weighted_mean=psi_mean,
    psi_raw_sum=raw_sum,total_weight=sum_w,n_batches=nb1,max_batch_rows_seen=max_B,
    n_input_rows=ni,n_output_rows=no,n_positive_contributing_rows=npos,materialized_target_if=materialized,
    replay_verified_scope="Exact stored projected rows and complete IDs/weights/counts, across two measured passes. Raw IF/nuisance immutability and sampling independence remain caller assumptions.",
    sampling_provenance_certified=FALSE)
}

# ---------------------------------------------------------------------------
# Cholesky-stage refusal state (diagnostics only; October 17, 2026)
# ---------------------------------------------------------------------------
# Compact record attached to the unchanged `cdid_gmm_pd_error` raised when the
# unpivoted Cholesky factorization of the weight covariance `V` fails. Every
# field is already computed state of this call or a direct inspection of it:
# no eigenvalue, condition number, pivoting, regularization, rank threshold or
# retry is evaluated. The V/W_G matrix payloads reuse objects already
# materialized under this call's preflight (`V validation buffers`, `gauge`
# allocations); small named-vector/list metadata may allocate. No large matrix
# or IF payload is duplicated. The record distinguishes the
# parameter dimension `p` of the retained design from its certified graph rank
# `k = p - nullity` (rank of `W_G`); neither is the rank of the sample moment
# covariance, which this refusal does not determine.
gmm_cholesky_refusal_state <- function(V, mhat_R, W_G, gauge, design, provenance, target_info,
                                       chol_result, chol_message) {
  rk <- design$rank_diagnostics
  scalar_rank <- rk[setdiff(names(rk), "by_cohort")]
  outcome <- if (is.null(chol_result)) "factorization_error" else
    if (anyNA(chol_result) || any(!is.finite(chol_result))) "nonfinite_factor" else "nonpositive_pivot"
  dV <- diag(V)
  list(
    stage = "optimal_weight_cholesky_whitening",
    factorization = list(routine = "base::chol (unpivoted, no fallback)", outcome = outcome,
      error_message = if (is.null(chol_message)) NA_character_ else as.character(chol_message)),
    dimensions = list(m_R = gauge$n_links, p = gauge$n_parameters, k = gauge$rank_k, nullity = gauge$nullity_nu,
      n_targets = length(target_info$estimable), J_est = target_info$n_estimable,
      M_positive_rows = provenance$n_positive_contributing_rows,
      meaning = paste0("m_R retained D10/D13 moments; p = parameter (level) dimension of the retained design; ",
        "k = p - nullity = certified structural design rank (rank of W_G); J_est = estimable requested loading columns. ",
        "None of these is the rank of the m_R x m_R sample covariance V.")),
    design_rank = c(scalar_rank, list(by_cohort = rk$by_cohort)),
    covariance = list(provenance = unclass(provenance), exactly_symmetric = TRUE,
      diag_min = min(dV), diag_max = max(dV), link_ids = gauge$link_ids),
    moment_state = list(mhat = stats::setNames(as.double(mhat_R), gauge$link_ids), V = V, W_G = W_G,
      param_ids = gauge$param_ids, keep_param_ids = gauge$keep_param_ids, root_params = gauge$root_params,
      level_estimable = gauge$level_estimable, is_isolated = gauge$is_isolated,
      target_ids = target_info$target_ids, target_estimable = target_info$estimable),
    retained_payload = list(entries = as.double(gauge$n_links)^2 + as.double(gauge$n_links) * gauge$rank_k + as.double(gauge$n_links),
      admitted_by = "this call's `V validation buffers` (m_R x m_R) and gauge allocation preflights; V/W_G matrix payloads reused, small vector/list metadata may allocate"),
    meaning = paste0("Factorization refusal before the accepted PD/transfer certificate. It establishes neither exact ",
      "singularity of V, nor the population covariance rank, nor ATT nonidentification; no fallback, pivoting, ",
      "regularization or identity weighting is substituted.")
  )
}

# ---------------------------------------------------------------------------
# Main private exact-GMM computation core entrypoint
# ---------------------------------------------------------------------------

solve_chained_gmm_core <- function(
  design,
  mhat,
  V,
  provenance,
  targets,
  weighting,
  if_replay,
  accuracy,
  bounds,
  gauge_rule,
  custom_gauge_roots = NULL,
  request_param_cov = FALSE,
  materialize_target_if = FALSE,
  gmm_covariance = NULL,
  covariance_tuning = NULL,
  reduction_spec = NULL
) {
  if (missing(design) || missing(mhat) || missing(V) || missing(provenance) ||
      missing(targets) || missing(weighting) || missing(if_replay) ||
      missing(accuracy) || missing(bounds) || missing(gauge_rule)) {
    cdid_gmm_abort(
      "All core arguments (`design`, `mhat`, `V`, `provenance`, `targets`, `weighting`, `if_replay`, `accuracy`, `bounds`, `gauge_rule`) must be supplied explicitly.",
      "cdid_gmm_control_error"
    )
  }
  w_choice <- resolve_gmm_weighting_choice(weighting, arg = "weighting")
  weighting_requested <- w_choice$requested
  weighting <- w_choice$canonical
  if (!gmm_scalar_logical(request_param_cov) || !gmm_scalar_logical(materialize_target_if)) {
    cdid_gmm_abort("`request_param_cov` and `materialize_target_if` must be scalar logicals.", "cdid_gmm_control_error")
  }

  nb_info <- if (identical(accuracy$method, "conventional_chol_qr"))
    list(backend = "base_R_BLAS_LAPACK_computed_diagnostics") else gmm_native_backend_info()
  bounds <- validate_gmm_bounds_control(bounds)
  accuracy <- validate_gmm_accuracy_control(accuracy, weighting)
  cov_strategy <- accuracy$covariance_strategy
  if (!is.null(gmm_covariance) || !is.null(covariance_tuning)) {
    cov_method <- if (!is.null(gmm_covariance)) {
      resolve_gmm_covariance_choice(gmm_covariance)
    } else if (!is.null(cov_strategy)) {
      cov_strategy$method
    } else {
      "tikhonov"
    }
    tuning <- resolve_gmm_covariance_tuning(covariance_tuning)
    if (!is.null(cov_strategy) && is.null(covariance_tuning)) {
      tuning <- list(
        condition_trigger = cov_strategy$condition_trigger,
        ridge_relative = cov_strategy$ridge_relative,
        overrides = cov_strategy$overrides
      )
    }
    cov_strategy <- gmm_covariance_strategy(
      method = cov_method,
      condition_trigger = tuning$condition_trigger,
      ridge_relative = tuning$ridge_relative,
      overrides = tuning$overrides,
      weighting_requested = weighting_requested,
      weighting_canonical = weighting
    )
  } else if (!is.null(cov_strategy)) {
    cov_strategy$weighting_requested <- weighting_requested
    cov_strategy$weighting_canonical <- weighting
    cov_strategy$alias_used <- w_choice$alias_used
  }
  provenance <- validate_gmm_covariance_provenance(provenance, weighting)
  if (weighting == "optimal" && !is.null(cov_strategy) && provenance$row_level != "unit") {
    cdid_gmm_abort("Automatic covariance handling requires complete unit-IF covariance; clusters change inference only.",
      "cdid_gmm_covariance_error")
  }
  if (!is.null(reduction_spec)) {
    validate_gmm_reduction_spec(
      if (is.list(reduction_spec) && !is.null(reduction_spec$mode)) reduction_spec else list(mode = "sparse_class_certificate", class_certificate = reduction_spec),
      weighting
    )
  }
  if_replay <- validate_gmm_if_replay_spec(if_replay)

  if (!inherits(design, "cdid_chained_design")) cdid_gmm_abort("Invalid design.", "cdid_gmm_design_error")
  m_pre <- gmm_int_bound(design$incidence$n_links, "design links")
  p_pre <- gmm_int_bound(design$incidence$n_parameters, "design parameters")
  k_pre <- gmm_int_bound(design$rank_diagnostics$rank, "design rank")
  gmm_safe_product(20, m_pre + p_pre, bounds$max_matrix_entries, "gauge/design validation vectors")
  gmm_safe_product(m_pre, m_pre, bounds$max_matrix_entries, "V validation buffers")
  gmm_safe_product(8, if_replay$n_input_rows + if_replay$n_output_rows, bounds$max_matrix_entries, "replay IDs/weights/count buffers")
  J_pre <- if (is.matrix(targets)) ncol(targets) else if (is.data.frame(targets)) min(nrow(targets),bounds$max_targets) else length(targets)
  if (J_pre < 1 || J_pre > bounds$max_targets) cdid_gmm_abort("Invalid target count/budget.", "cdid_gmm_bound_error")
  eff_red <- if (!is.null(reduction_spec)) reduction_spec else provenance$reduction
  n_trip_pre <- if (!is.null(eff_red$class_certificate$E$i)) length(eff_red$class_certificate$E$i) else if (!is.null(eff_red$E$i)) length(eff_red$E$i) else 0L
  preflight_gmm_allocations(m_pre, k_pre, J_pre, weighting, bounds,
    n_out_rows = if_replay$n_output_rows, n_input_rows = if_replay$n_input_rows,
    n_parameters = p_pre, request_param_cov = request_param_cov, materialize_target_if = materialize_target_if,
    certify_error = !identical(accuracy$method, "conventional_chol_qr"),
    covariance_strategy = cov_strategy, n_dropped_triplets = n_trip_pre)

  # 1. Gauge construction and exact target estimability
  gauge <- build_gmm_gauge(design, gauge_rule = gauge_rule, custom_roots = custom_gauge_roots)
  m_R <- gauge$n_links
  p <- gauge$n_parameters
  k <- gauge$rank_k

  if (m_R < k) {
    cdid_gmm_abort(
      paste0("Underidentified retained system: m_R (", m_R, ") < rank k (", k, ")."),
      "cdid_gmm_design_error",
      diagnostics = list(m_R = m_R, k = k, p = p)
    )
  }

  # Validate mhat and V
  if (!is.atomic(mhat) || !is.null(dim(mhat)) || is.object(mhat) || !is.numeric(mhat) ||
      length(mhat) != m_R || anyNA(mhat) || any(!is.finite(mhat))) {
    cdid_gmm_abort("`mhat` must be a finite numeric vector of length m_R.", "cdid_gmm_design_error")
  }
  if (!is.null(names(mhat)) && !identical(names(mhat), gauge$link_ids)) {
    cdid_gmm_abort("Named `mhat` does not match `design$incidence$link_ids` order.", "cdid_gmm_design_error")
  }
  mhat_R <- as.double(mhat)

  if (!is.matrix(V) || !is.numeric(V) || is.object(V) || nrow(V) != m_R || ncol(V) != m_R ||
      anyNA(V) || any(!is.finite(V))) {
    cdid_gmm_abort("`V` must be a finite numeric m_R x m_R matrix.", "cdid_gmm_covariance_error")
  }
  if (isTRUE(attr(V, "is_pd")) || isTRUE(attr(V, "certified_pd")) || isTRUE(attr(V, "bypass_pd_guard"))) {
    cdid_gmm_abort("Untrusted PD metadata attributes on `V` are refused.", "cdid_gmm_covariance_error")
  }
  for (jj in seq_len(m_R)) if (!identical(as.double(V[,jj]), as.double(V[jj,])))
    cdid_gmm_abort("Weight covariance `V` must be exactly symmetric.", "cdid_gmm_covariance_error")
  if ((!is.null(rownames(V)) && !identical(rownames(V),gauge$link_ids)) ||
      (!is.null(colnames(V)) && !identical(colnames(V),gauge$link_ids)))
    cdid_gmm_abort("Covariance row/column labels disagree with incidence order.", "cdid_gmm_covariance_error")

  # Rank bound check for optimal weighting:
  # Under legacy unregularized private controls (`is.null(cov_strategy)`), m_R > max_rank_allowed
  # is refused immediately. Under D26 (`!is.null(cov_strategy)`), sample-rank deficiency triggers
  # the selected covariance handler rather than blanket refusal.
  M_pos <- provenance$n_positive_contributing_rows
  max_rank_allowed <- if (provenance$centered) M_pos - 1 else M_pos
  if (is.null(cov_strategy) && weighting == "optimal" && m_R > max_rank_allowed) {
    cdid_gmm_abort(
      paste0("Sample covariance rank bound refused: m_R = ", m_R, " exceeds ",
             if (provenance$centered) "M - 1 = " else "M = ", max_rank_allowed,
             " for M = ", M_pos, " positive contributing ", provenance$row_level, " rows."),
      "cdid_gmm_rank_bound_error",
      diagnostics = list(
        m_R = m_R,
        k = k,
        M_positive_rows = M_pos,
        centered = provenance$centered,
        row_level = provenance$row_level,
        max_rank_allowed = max_rank_allowed
      )
    )
  }

  if (if_replay$input_row_level == "unit" && provenance$n_units != if_replay$n_input_rows)
    cdid_gmm_abort("Full unit counts disagree between covariance and replay.", "cdid_gmm_replay_error")
  if (provenance$row_level == "cluster" && if_replay$output_row_level == "cluster" && provenance$n_clusters != if_replay$n_output_rows)
    cdid_gmm_abort("Full cluster counts disagree between covariance and replay.", "cdid_gmm_replay_error")

  # Scale & convention compatibility check between weight covariance V and inference replay
  convention_matches <- (
    weighting == "optimal" &&
    provenance$row_level == if_replay$output_row_level &&
    provenance$centered == if_replay$centered &&
    provenance$weighted == if_replay$weighted &&
    provenance$weight_type == if_replay$weight_type &&
    provenance$n_positive_contributing_rows == if_replay$n_positive_contributing_rows &&
    provenance$scale_id == if_replay$scale_id &&
    provenance$nuisance_fit_id == if_replay$nuisance_fit_id
  )
  scale_matches_inference <- convention_matches &&
    identical(provenance$normalization_factor, if_replay$normalization_factor)

  if (isTRUE(accuracy$require_matching_scale) && !scale_matches_inference) {
    cdid_gmm_abort(
      "Weight covariance convention/scale does not match inference replay while `require_matching_scale` is TRUE.",
      "cdid_gmm_scale_error",
      diagnostics = list(
        weight_provenance = provenance,
        replay_provenance = list(
          output_row_level = if_replay$output_row_level,
          centered = if_replay$centered,
          weighted = if_replay$weighted,
          weight_type = if_replay$weight_type,
          normalization_factor = if_replay$normalization_factor,
          scale_id = if_replay$scale_id
        )
      )
    )
  }

  target_info <- validate_gmm_targets(design, gauge, targets, bounds)
  J_est <- target_info$n_estimable

  # 2. Preflight allocation check before constructing dense W_G, R, A, Lambda
  alloc_diag <- preflight_gmm_allocations(
    m_R = m_R,
    k = k,
    J_est = J_est,
    weighting = weighting,
    bounds = bounds,
    n_out_rows = if_replay$n_output_rows,
    n_input_rows = if_replay$n_input_rows, n_parameters = p,
    request_param_cov = request_param_cov,
    materialize_target_if = materialize_target_if,
    certify_error = !identical(accuracy$method, "conventional_chol_qr"),
    covariance_strategy = cov_strategy,
    n_dropped_triplets = n_trip_pre
  )

  # 3. Materialize W_G and solve GMM system
  W_G <- materialize_gmm_gauge_matrix(gauge, bounds$max_matrix_entries)
  gram_info <- if (weighting == "optimal" && !is.null(cov_strategy)) assemble_gmm_gram_matrix(
    m_R = m_R,
    weighting = weighting,
    reduction = eff_red,
    max_matrix_entries = bounds$max_matrix_entries
  ) else list(G = NULL, C = NULL, M_full = m_R, L_retained = m_R,
    reduction_active = FALSE, is_identity_gram = TRUE, is_diagonal_gram = TRUE,
    n_triplets = 0L, kept_indices = seq_len(m_R), gram_provenance = "identity_unreduced")
  diag_B <- NULL
  handler_used <- "none"
  trigger_reason <- if (weighting == "identity") "not_applicable_identity" else "legacy_unregularized"
  lambda_used <- NA_real_
  weighted_design_rank <- k

  if (weighting == "optimal" && !is.null(cov_strategy)) {
    diag_B <- diagnose_gmm_transformed_covariance(
      V = V,
      gram_info = gram_info,
      provenance = provenance,
      strategy = cov_strategy,
      max_matrix_entries = bounds$max_matrix_entries
    )
    if (m_R > max_rank_allowed && !diag_B$triggered) {
      diag_B$triggered <- TRUE
      diag_B$trigger_reason <- "known_unit_covariance_rank_bound"
    }
    trigger_reason <- diag_B$trigger_reason
    if (!diag_B$triggered) {
      # Attempt ordinary unregularized Cholesky solve on well-conditioned PD covariance
      R_chol_try <- if (any(diag(V) <= 0.0)) NULL else tryCatch(base::chol(V), error = function(e) NULL)
      if (is.null(R_chol_try) || anyNA(R_chol_try) || any(!is.finite(R_chol_try)) || any(diag(R_chol_try) <= 0.0)) {
        diag_B$triggered <- TRUE
        diag_B$trigger_reason <- "cholesky_failure_on_valid_psd"
        trigger_reason <- diag_B$trigger_reason
      } else {
        R_chol <- R_chol_try
      }
    }
    if (diag_B$triggered) {
      h_res <- solve_gmm_covariance_handler(
        V = V,
        mhat_R = mhat_R,
        W_G = W_G,
        gauge = gauge,
        target_info = target_info,
        gram_info = gram_info,
        diag_B = diag_B,
        strategy = cov_strategy,
        accuracy = accuracy,
        max_matrix_entries = bounds$max_matrix_entries
      )
      handler_used <- h_res$handler_used
      lambda_used <- h_res$lambda
      weighted_design_rank <- h_res$weighted_design_rank
      qa <- h_res$qa
      U_qr <- h_res$U_qr
      piv <- h_res$piv
      R_chol <- h_res$R_chol
      A_whitened <- h_res$A_whitened
      b_whitened <- h_res$b_whitened
      theta_keep <- h_res$theta_keep
      Z1 <- h_res$Z1
      Lambda <- h_res$Lambda
      Mc <- h_res$Mc
      y_mult <- h_res$y_mult
      guard <- h_res$guard
      raw_Z1tZ1 <- h_res$raw_Z1tZ1
      raw_quad <- h_res$raw_quad
      target_info <- h_res$target_info
      J_est <- target_info$n_estimable
    }
  }

  if (identical(handler_used, "none")) {
    if (weighting == "optimal") {
      if (any(diag(V) <= 0)) {
        cdid_gmm_abort(
          "Weight covariance `V` has nonpositive diagonal entries.",
          "cdid_gmm_pd_error",
          diagnostics = list(min_diag_V = min(diag(V)))
        )
      }
      if (is.null(cov_strategy)) {
        chol_message <- NULL
        R_chol <- tryCatch(base::chol(V), error = function(e) { chol_message <<- conditionMessage(e); NULL })
        if (is.null(R_chol) || anyNA(R_chol) || any(!is.finite(R_chol)) || any(diag(R_chol) <= 0)) {
          cdid_gmm_abort(
            "Unpivoted Cholesky factorization failed on `V`; no fallback or regularization is permitted.",
            "cdid_gmm_pd_error",
            diagnostics = list(m_R = m_R, k = k,
              refusal_state = gmm_cholesky_refusal_state(V, mhat_R, W_G, gauge, design, provenance, target_info,
                                                         R_chol, chol_message))
          )
        }
      }
      A_whitened <- backsolve(R_chol, W_G, transpose = TRUE)
      b_whitened <- backsolve(R_chol, mhat_R, transpose = TRUE)
      if (any(!is.finite(A_whitened)) || any(!is.finite(b_whitened))) {
        cdid_gmm_abort("Nonfinite values in Cholesky-whitened system.", "cdid_gmm_pd_error")
      }
    } else {
      R_chol <- NULL
      A_whitened <- W_G
      b_whitened <- mhat_R
    }

    qa <- base::qr(A_whitened, LAPACK = TRUE)
    U_qr <- base::qr.R(qa)
    piv <- qa$pivot
    if (anyNA(U_qr) || any(!is.finite(U_qr)) || any(diag(U_qr) == 0)) {
      cdid_gmm_abort(
        "LAPACK QR factor U has a zero or nonfinite pivot.",
        "cdid_gmm_qr_error",
        diagnostics = list(min_abs_diag_U = min(abs(diag(U_qr))))
      )
    }

    qtb <- base::qr.qty(qa, b_whitened)
    theta_keep <- numeric(k)
    theta_keep[piv] <- backsolve(U_qr, qtb[seq_len(k)])
    if (any(!is.finite(theta_keep))) {
      cdid_gmm_abort("Nonfinite parameter solution `theta` in triangular solve.", "cdid_gmm_qr_error")
    }
    names(theta_keep) <- gauge$keep_param_ids

    if (J_est > 0L) {
      C_piv <- target_info$C_gauge[piv, , drop = FALSE]
      Z1 <- backsolve(U_qr, C_piv, transpose = TRUE)
      Z2 <- base::qr.qy(qa, rbind(Z1, matrix(0, nrow = m_R - k, ncol = J_est)))
      Lambda <- if (weighting == "optimal") backsolve(R_chol, Z2) else Z2
      Mc <- matrix(0, nrow = k, ncol = J_est, dimnames = list(gauge$keep_param_ids, colnames(target_info$C_gauge)))
      Mc[piv, ] <- backsolve(U_qr, Z1)
      rownames(Lambda) <- gauge$link_ids
      colnames(Lambda) <- colnames(target_info$C_gauge)
      colnames(Z1) <- colnames(target_info$C_gauge)
      if (any(!is.finite(Z1)) || any(!is.finite(Lambda)) || any(!is.finite(Mc))) {
        cdid_gmm_abort("Nonfinite target loadings `Lambda` or `Mc`.", "cdid_gmm_qr_error")
      }
    } else {
      Z1 <- matrix(0, nrow = k, ncol = 0L)
      Lambda <- matrix(0, nrow = m_R, ncol = 0L, dimnames = list(gauge$link_ids, character(0)))
      Mc <- matrix(0, nrow = k, ncol = 0L, dimnames = list(gauge$keep_param_ids, character(0)))
    }

    e_resid <- mhat_R - gmm_gauge_multiply(gauge, theta_keep, transpose = FALSE)
    y_mult <- if (weighting == "optimal") {
      backsolve(R_chol, backsolve(R_chol, e_resid, transpose = TRUE))
    } else {
      e_resid
    }
    if (any(!is.finite(y_mult))) {
      cdid_gmm_abort("Nonfinite KKT multiplier `y`.", "cdid_gmm_qr_error")
    }

    # 4. Complete numerical guard (PD, conditioning, KKT & loading residuals, certified target error)
    acc_guard <- accuracy
    if (!is.null(cov_strategy) && weighting == "optimal") {
      # rho(B) was already checked on B = C V C' (B = V for identity G).
      # A values-only recomputation must not contradict the selected ordinary route.
      ev_V <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
      if (all(is.finite(ev_V)) && min(ev_V) > 0.0) {
        acc_guard$min_certified_lambda_rel <- min(acc_guard$min_certified_lambda_rel, (min(ev_V) / max(ev_V)) * 0.5)
      }
    }
    guard <- certify_gmm_numerical_guard(
      V_R = V,
      R_chol = R_chol,
      W_G = W_G,
      A_whitened = A_whitened,
      U_qr = U_qr,
      piv = piv,
      mhat_R = mhat_R,
      theta = theta_keep,
      y = y_mult,
      C_gauge = target_info$C_gauge,
      Lambda = Lambda,
      Mc = Mc,
      Z1 = Z1,
      weighting = weighting,
      accuracy = acc_guard,
      max_matrix_entries = bounds$max_matrix_entries
    )

    raw_Z1tZ1 <- if (J_est > 0L) crossprod(Z1) else matrix(0, 0L, 0L)
    raw_quad <- if (m_R > k) sum(qtb[-seq_len(k)]^2) else 0
    if (any(!is.finite(raw_Z1tZ1)) || !is.finite(raw_quad))
      cdid_gmm_abort("Criterion covariance or quadratic residual overflowed.", "cdid_gmm_accuracy_error")
  }

  cov_handling <- build_gmm_covariance_handling_record(
    weighting = weighting,
    strategy = cov_strategy,
    gram_info = gram_info,
    diag_B = diag_B,
    handler_used = handler_used,
    trigger_reason = trigger_reason,
    lambda = lambda_used,
    weighted_design_rank = weighted_design_rank,
    gauge = gauge,
    target_info = target_info,
    provenance = provenance
  )

  # 5. Streamed IF projection and cluster sandwich accumulation
  stream_res <- project_gmm_if_stream(
    if_replay = if_replay,
    Lambda = Lambda,
    m_R = m_R,
    bounds = bounds,
    provenance = provenance,
    materialize_target_if = materialize_target_if
  )

  # 6. Assemble target and parameter outputs, respecting NA rules for non-estimable levels
  J_total <- target_info$n_targets
  est_vals <- rep(NA_real_, J_total)
  se_sand_vals <- rep(NA_real_, J_total)
  se_model_vals <- rep(NA_real_, J_total)
  se_criterion_vals <- rep(NA_real_, J_total)
  err_est_vals <- rep(NA_real_, J_total)
  err_cert_vals <- rep(NA_real_, J_total)
  loading_err_cert_vals <- rep(NA_real_, J_total)
  point_cert_vals <- rep(NA, J_total)
  loading_cert_vals <- rep(NA, J_total)
  complete_inf_cert_vals <- rep(FALSE, J_total)

  if (J_est > 0L) {
    tau_est <- guard$computed_target_values
    est_idx <- target_info$estimable_indices
    est_vals[est_idx] <- tau_est
    se_sand_vals[est_idx] <- stream_res$se_sandwich
    err_est_vals[est_idx] <- guard$targets$target_error_est
    err_cert_vals[est_idx] <- guard$targets$target_error_cert
    loading_err_cert_vals[est_idx] <- guard$targets$loading_error_cert
    point_cert_vals[est_idx] <- rep(isTRUE(guard$point_certified), J_est)
    loading_cert_vals[est_idx] <- rep(isTRUE(guard$loading_certified), J_est)
    complete_inf_cert_vals[est_idx] <- rep(isTRUE(guard$complete_inference_certified), J_est)
    if (weighting == "optimal") se_criterion_vals[est_idx] <- guard$targets$se_criterion_est
  }

  target_table <- data.frame(
    target_id = target_info$target_ids,
    estimable = target_info$estimable,
    status = target_info$status,
    support_reason = target_info$support_reason,
    estimate = est_vals,
    se_sandwich = se_sand_vals,
    se_model = se_model_vals,
    se_criterion_conditional = se_criterion_vals,
    target_error_est = err_est_vals,
    target_error_cert = err_cert_vals,
    loading_error_cert = loading_err_cert_vals,
    G_lambda = loading_err_cert_vals,
    point_certified = point_cert_vals,
    loading_certified = loading_cert_vals,
    complete_inference_certified = complete_inf_cert_vals,
    stringsAsFactors = FALSE
  )

  # Parameter table: only connected_to_ref (level_estimable) parameters have level estimates;
  # unanchored roots, non-root unanchored vertices, and isolated vertices remain strictly NA!
  param_table <- design$parameters[, c("param_index", "param_id", "g", "time", "role",
                                       "component_id", "is_isolated", "level_estimable"), drop = FALSE]
  param_table$status <- ifelse(
    param_table$level_estimable,
    "estimable",
    ifelse(param_table$is_isolated, "isolated_unsupported", "disconnected_component")
  )
  param_table$is_gauge_root <- param_table$param_index %in% gauge$root_params
  param_table$level_estimate <- NA_real_
  anchored_keep <- which(gauge$level_estimable[gauge$keep_params])
  if (length(anchored_keep) > 0L) {
    param_table$level_estimate[gauge$keep_params[anchored_keep]] <- theta_keep[anchored_keep]
  }

  # Declared metadata cannot verify the actual construction or sampling meaning of V.
  # Raw criterion quantities remain available; qualified model covariance/SE/J do not.
  model_cov <- NULL
  variance_gap_abs <- if (scale_matches_inference && J_est > 0L) {
    max(abs(stream_res$sandwich_cov - raw_Z1tZ1))
  } else {
    NA_real_
  }
  variance_gap_rel <- if (scale_matches_inference && J_est > 0L && max(abs(raw_Z1tZ1)) > 0) {
    variance_gap_abs / max(abs(raw_Z1tZ1))
  } else {
    NA_real_
  }

  J_stat_optimal <- NA_real_

  param_cov_gauge <- NULL
  if (isTRUE(request_param_cov)) {
      P_cov <- matrix(0, nrow = k, ncol = k, dimnames = list(gauge$keep_param_ids, gauge$keep_param_ids))
      for (ii in seq_len(k)) {
        ei <- numeric(k); ei[ii] <- 1
        zi <- backsolve(U_qr, ei[piv], transpose = TRUE)
        for (jj in seq_len(ii)) {
          ej <- numeric(k); ej[jj] <- 1
          zj <- backsolve(U_qr, ej[piv], transpose = TRUE)
          P_cov[ii,jj] <- P_cov[jj,ii] <- sum(zi*zj)
        }
      }
      if (any(!is.finite(P_cov))) cdid_gmm_abort("Parameter criterion covariance overflowed.", "cdid_gmm_accuracy_error")
      param_cov_gauge <- P_cov
  }

  alloc_diag$measured_retained_bytes <- list(
    W_G = as.numeric(utils::object.size(W_G)),
    A_whitened = as.numeric(utils::object.size(A_whitened)),
    qr_object = if (!is.null(qa)) as.numeric(utils::object.size(qa)) else 0,
    R_chol = if (!is.null(R_chol)) as.numeric(utils::object.size(R_chol)) else 0,
    Lambda = as.numeric(utils::object.size(Lambda)),
    sandwich_cov = as.numeric(utils::object.size(stream_res$sandwich_cov))
  )

  structure(
    list(
      targets = target_table,
      parameters = param_table,
      witnesses = stats::setNames(target_info$witnesses, target_info$target_ids),
      sandwich_cov = stream_res$sandwich_cov,
      model_cov = model_cov,
      raw_Z1tZ1 = raw_Z1tZ1,
      J_stat_optimal = J_stat_optimal,
      raw_quadratic_residual = raw_quad,
      overid_df = as.integer(m_R - weighted_design_rank),
      reference_distribution_certified = FALSE,
      model_covariance_provenance_certified = FALSE,
      point_certified = guard$point_certified,
      loading_certified = guard$loading_certified,
      loading_budget_supplied = guard$loading_budget_supplied,
      loading_validation_status = guard$loading_validation_status,
      complete_inference_certified = FALSE,
      full_inference_status = guard$full_inference_status,
      propagation_bounds = guard$propagation_bounds,
      covariance_handling = cov_handling,
      provenance_scope = "Caller assumptions only: labels do not verify independent sampling rows, actual nuisance identity, or V construction. Model covariance/J are conditional algebra.",
      scale_matches_inference = FALSE,
      declared_scale_matches_inference = scale_matches_inference,
      criterion_interpretation = "raw_Z1tZ1, se_criterion_conditional, internal param_cov_gauge and raw_quadratic_residual describe the supplied criterion only. A sampling/model/J interpretation additionally requires actual covariance construction and sampling assumptions, which metadata do not verify.",
      variance_gap = list(abs = variance_gap_abs, rel = variance_gap_rel),
      internal_gauge = list(
        gauge = gauge,
        theta_gauge_coords = theta_keep,
        y_multiplier = y_mult,
        Lambda = Lambda,
        Mc = Mc,
        Z1 = Z1,
        U_qr = U_qr,
        pivot = piv,
        R_chol = R_chol,
        param_cov_gauge = param_cov_gauge,
        note = "Internal gauge coordinates on unanchored components are root contrasts, never ATT levels."
      ),
      if_projection = stream_res,
      numerical_guard = guard,
      allocations = alloc_diag,
      provenance = list(
        weighting = weighting,
        weighting_requested = weighting_requested,
        native_backend = nb_info,
        weight_covariance = provenance,
        inference_replay = list(
          input_row_level = if_replay$input_row_level,
          output_row_level = if_replay$output_row_level,
          aggregate_to_clusters = if_replay$aggregate_to_clusters,
          n_input_rows = if_replay$n_input_rows,
          n_output_rows = if_replay$n_output_rows,
          n_positive_contributing_rows = if_replay$n_positive_contributing_rows,
          centered = if_replay$centered,
          weighted = if_replay$weighted,
          weight_type = if_replay$weight_type,
          normalization_factor = if_replay$normalization_factor,
          scale_id = if_replay$scale_id,
          nuisance_fit_id = if_replay$nuisance_fit_id
        ),
        reduction = provenance$reduction
      )
    ),
    class = "cdid_gmm_core_result"
  )
}

chained_gmm_core <- solve_chained_gmm_core

