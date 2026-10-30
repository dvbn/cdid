# D26/D27: feasible two-step GMM (optimal; two_step alias), explicit identity,
# and automatic covariance-ridge Tikhonov after proved exact reductions.
# G = E'E = C'C, B = C V_r C', s = tr(B)/M, lambda = ridge_relative*s.
# No ATT penalty, spectral clipping, automatic identity or rank repair.

gmm_weighting_public_choices <- c("optimal", "identity", "two_step")
gmm_weighting_canonical_choices <- c("optimal", "identity")
gmm_covariance_choices <- "tikhonov"
gmm_covariance_tuning_names <- c("condition_trigger", "ridge_relative")

gmm_covariance_default_tuning <- function() {
  list(condition_trigger = 1e-10, ridge_relative = 1e-4)
}

resolve_gmm_weighting_choice <- function(
  weighting,
  arg = "gmm_weighting",
  error_fn = cdid_gmm_abort,
  error_class = "cdid_gmm_control_error"
) {
  if (missing(weighting) || is.null(weighting) || !is.character(weighting) ||
      is.object(weighting) || !is.null(dim(weighting)) || length(weighting) < 1L ||
      anyNA(weighting) || any(!nzchar(weighting))) {
    error_fn(
      paste0("`", arg, "` must be one of 'optimal', 'identity', or 'two_step'."),
      error_class
    )
  }
  if (length(weighting) > 1L) {
    if (identical(weighting, c("optimal", "identity")) ||
        identical(weighting, c("optimal", "identity", "two_step"))) {
      weighting <- weighting[1L]
    } else {
      error_fn(
        paste0("`", arg, "` must be a single string ('optimal', 'identity', or 'two_step')."),
        error_class
      )
    }
  }
  if (!weighting %in% gmm_weighting_public_choices) {
    error_fn(
      paste0("`", arg, "` must be 'optimal', 'identity', or 'two_step'."),
      error_class
    )
  }
  weighting <- unname(weighting)
  canonical <- if (identical(weighting, "two_step")) "optimal" else weighting
  list(
    requested = weighting,
    canonical = canonical,
    alias_used = identical(weighting, "two_step")
  )
}

resolve_gmm_covariance_choice <- function(
  gmm_covariance = "tikhonov", arg = "gmm_covariance",
  error_fn = cdid_gmm_abort, error_class = "cdid_gmm_control_error"
) {
  if (identical(gmm_covariance, "moore_penrose")) {
    error_fn("Moore-Penrose covariance handling is deferred; `gmm_covariance` currently accepts only 'tikhonov'.", error_class)
  }
  if (!gmm_scalar_character(gmm_covariance) || !identical(gmm_covariance, "tikhonov")) {
    error_fn(paste0("`", arg, "` must be the single string 'tikhonov'."), error_class)
  }
  gmm_covariance
}

resolve_gmm_covariance_tuning <- function(
  control_list = NULL, error_fn = cdid_gmm_abort,
  error_class = "cdid_gmm_control_error"
) {
  out <- gmm_covariance_default_tuning()
  if (!is.null(control_list) && (!is.list(control_list) || is.object(control_list) ||
      (length(control_list) && (is.null(names(control_list)) || anyNA(names(control_list)) ||
       any(!nzchar(names(control_list))) || anyDuplicated(names(control_list)))))) {
    error_fn("Covariance tuning must be NULL or a plain uniquely named list.", error_class)
  }
  # Other advanced public fields are validated by the public resolver.
  if (any(c("mp_cutoff", "moore_penrose_cutoff") %in% names(control_list))) {
    error_fn("Moore-Penrose tuning is deferred and is not accepted.", error_class)
  }
  if (any(!names(control_list) %in% names(out))) {
    error_fn("Unknown covariance tuning entry.", error_class)
  }
  overrides <- as.character(intersect(names(out), names(control_list)))
  for (nm in overrides) {
    val <- control_list[[nm]]
    if (!gmm_scalar_numeric(val, positive = TRUE) || (nm == "condition_trigger" && val >= 1)) {
      error_fn(paste0("`control$", nm, "` must be a finite positive numeric scalar",
        if (nm == "condition_trigger") " strictly below 1." else "."), error_class)
    }
    out[[nm]] <- as.double(val)
  }
  c(out, list(overrides = overrides))
}

gmm_covariance_strategy <- function(
  method = "tikhonov", condition_trigger = 1e-10, ridge_relative = 1e-4,
  overrides = character(0L), weighting_requested = "optimal",
  weighting_canonical = "optimal"
) {
  method <- resolve_gmm_covariance_choice(method)
  w <- resolve_gmm_weighting_choice(weighting_requested)
  if (!identical(weighting_canonical, w$canonical)) {
    cdid_gmm_abort("Inconsistent requested/canonical weighting in covariance strategy.", "cdid_gmm_control_error")
  }
  tuning <- resolve_gmm_covariance_tuning(list(condition_trigger = condition_trigger, ridge_relative = ridge_relative))
  if (!is.character(overrides) || is.object(overrides) || !is.null(dim(overrides)) || anyNA(overrides) ||
      anyDuplicated(overrides) || any(!overrides %in% gmm_covariance_tuning_names)) {
    cdid_gmm_abort("Invalid covariance tuning override record.", "cdid_gmm_control_error")
  }
  structure(list(method = method, condition_trigger = tuning$condition_trigger,
    ridge_relative = tuning$ridge_relative, overrides = overrides,
    weighting_requested = w$requested, weighting_canonical = w$canonical,
    alias_used = w$alias_used), class = "cdid_gmm_covariance_strategy")
}

validate_gmm_covariance_strategy <- function(strategy, weighting = "optimal") {
  if (is.null(strategy)) return(NULL)
  if (is.character(strategy)) return(gmm_covariance_strategy(method = strategy,
    weighting_requested = weighting, weighting_canonical = resolve_gmm_weighting_choice(weighting)$canonical))
  required <- c("method", "condition_trigger", "ridge_relative", "overrides",
    "weighting_requested", "weighting_canonical", "alias_used")
  if (!inherits(strategy, "cdid_gmm_covariance_strategy") || !is.list(strategy) ||
      anyDuplicated(names(strategy)) || !setequal(names(strategy), required)) {
    cdid_gmm_abort("Invalid or obsolete covariance strategy; only Tikhonov is implemented.", "cdid_gmm_control_error")
  }
  gmm_covariance_strategy(method = strategy$method, condition_trigger = strategy$condition_trigger,
    ridge_relative = strategy$ridge_relative, overrides = strategy$overrides,
    weighting_requested = strategy$weighting_requested, weighting_canonical = strategy$weighting_canonical)
}

gmm_scale_safe_frobenius <- function(X) {
  if (!length(X)) return(0.0)
  if (anyNA(X) || any(!is.finite(X))) return(NaN)
  m_max <- max(abs(X))
  if (m_max == 0.0) return(0.0)
  Z <- as.double(X) / m_max
  m_max * sqrt(sum(Z * Z))
}

# ------------------------------------------------------------------------------
# Sparse Gram reconstruction G = E'E = I_L + E_drop' E_drop in `kept_idx` order
# without allocating the M x L dense matrix E or M x M unreduced covariance.
# ------------------------------------------------------------------------------

assemble_gmm_gram_matrix <- function(m_R, weighting, reduction = NULL, max_matrix_entries = NULL) {
  L <- gmm_int_bound(m_R, "m_R")
  if (L > .Machine$integer.max) cdid_gmm_abort("Gram dimension exceeds R's matrix dimension range.", "cdid_gmm_bound_error")
  L <- as.integer(L)
  if (!is.null(max_matrix_entries)) {
    max_matrix_entries <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
    gmm_safe_product(6, gmm_safe_product(L, L, max_matrix_entries, "Gram matrix G"), max_matrix_entries, "Gram G/C and assembly/factor scratch")
  }
  if (is.null(reduction)) {
    I_L <- diag(1.0, L)
    return(list(
      G = I_L,
      C = I_L,
      M_full = L,
      L_retained = L,
      reduction_active = FALSE,
      is_identity_gram = TRUE,
      is_diagonal_gram = TRUE,
      n_triplets = 0L,
      kept_indices = seq_len(L),
      gram_provenance = "identity_unreduced"
    ))
  }
  if (!is.list(reduction)) cdid_gmm_abort("Reconstruction specification must be a list.", "cdid_gmm_reduction_error")

  # Allow caller to supply an explicit Gram matrix `G` with `M_full` for low-level algebra checks,
  # or a `cdid_ipw_structural_certificate` / sparse triplet spec (`E`, `retained_index`, `dropped_index`, `n_links`).
  if (!is.null(reduction$G)) {
    G <- reduction$G
    if (!is.matrix(G) || !is.numeric(G) || is.object(G) || nrow(G) != L || ncol(G) != L ||
        anyNA(G) || any(!is.finite(G))) {
      cdid_gmm_abort("Explicit reconstruction Gram matrix `G` must be a finite L x L numeric matrix.", "cdid_gmm_reduction_error")
    }
    for (j in seq_len(L)) if (!identical(as.double(G[, j]), as.double(G[j, ]))) {
      cdid_gmm_abort("Explicit reconstruction Gram matrix `G` must be symmetric.", "cdid_gmm_reduction_error")
    }
    M_full <- if (!is.null(reduction$n_links_all)) gmm_int_bound(reduction$n_links_all, "n_links_all") else L
    if (M_full > .Machine$integer.max) cdid_gmm_abort("Full catalog dimension exceeds R's index range.", "cdid_gmm_bound_error")
    M_full <- as.integer(M_full)
    if (M_full < L) cdid_gmm_abort("`n_links_all` cannot be smaller than `m_R`.", "cdid_gmm_reduction_error")
    C <- tryCatch(base::chol(G), error = function(e) NULL)
    if (is.null(C) || anyNA(C) || any(!is.finite(C)) || any(diag(C) <= 0.0)) {
      cdid_gmm_abort("Reconstruction Gram matrix `G` is not positive definite.", "cdid_gmm_reduction_error")
    }
    is_id <- identical(unname(G), diag(1.0, L))
    is_diag <- all(vapply(seq_len(L), function(j) all(G[-j, j] == 0), logical(1L)))
    return(list(
      G = G,
      C = C,
      M_full = M_full,
      L_retained = L,
      reduction_active = (M_full > L || !is_id),
      is_identity_gram = is_id,
      is_diagonal_gram = is_diag,
      n_triplets = 0L,
      kept_indices = if (!is.null(reduction$kept_indices)) as.integer(reduction$kept_indices) else seq_len(L),
      gram_provenance = if (!is.null(reduction$gram_provenance)) as.character(reduction$gram_provenance) else "explicit_gram_matrix"
    ))
  }

  cert <- if (inherits(reduction, "cdid_ipw_structural_certificate")) {
    reduction
  } else if (!is.null(reduction$class_certificate)) {
    reduction$class_certificate
  } else if (!is.null(reduction$E) && !is.null(reduction$retained_index)) {
    reduction
  } else {
    NULL
  }

  if (is.null(cert)) {
    I_L <- diag(1.0, L)
    return(list(
      G = I_L,
      C = I_L,
      M_full = L,
      L_retained = L,
      reduction_active = FALSE,
      is_identity_gram = TRUE,
      is_diagonal_gram = TRUE,
      n_triplets = 0L,
      kept_indices = seq_len(L),
      gram_provenance = "identity_unreduced"
    ))
  }

  # If `reduce_structural_redundancy` was explicitly FALSE, all M_all links were retained.
  if (isFALSE(reduction$reduction_active) && identical(as.double(L), gmm_int_bound(cert$n_links, "certificate$n_links"))) {
    I_L <- diag(1.0, L)
    return(list(
      G = I_L,
      C = I_L,
      M_full = L,
      L_retained = L,
      reduction_active = FALSE,
      is_identity_gram = TRUE,
      is_diagonal_gram = TRUE,
      n_triplets = 0L,
      kept_indices = if (!is.null(reduction$kept_indices)) as.integer(reduction$kept_indices) else seq_len(L),
      gram_provenance = "identity_unreduced"
    ))
  }

  if (inherits(cert, "cdid_ipw_structural_certificate")) {
    if (!exists("ipw_seal_ok", mode = "function") || !ipw_seal_ok(cert, "certificate")) {
      cdid_gmm_abort("Structural class certificate failed seal verification.", "cdid_gmm_reduction_error")
    }
  }

  M_full <- gmm_int_bound(cert$n_links, "certificate$n_links")
  if (M_full > .Machine$integer.max) cdid_gmm_abort("Full catalog dimension exceeds R's index range.", "cdid_gmm_bound_error")
  M_full <- as.integer(M_full)
  valid_index <- function(x) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    !anyNA(x) && all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= M_full)
  if (!valid_index(cert$retained_index) || !valid_index(cert$dropped_index) ||
      (!is.null(reduction$kept_indices) && !valid_index(reduction$kept_indices))) {
    cdid_gmm_abort("Reconstruction indices must be finite catalog integers.", "cdid_gmm_reduction_error")
  }
  ret_idx <- as.integer(cert$retained_index)
  drop_idx <- as.integer(cert$dropped_index)
  if (length(ret_idx) != L || M_full < L || length(drop_idx) != (M_full - L) ||
      anyDuplicated(ret_idx) || anyDuplicated(drop_idx) ||
      length(intersect(ret_idx, drop_idx)) > 0L ||
      !setequal(c(ret_idx, drop_idx), seq_len(M_full))) {
    cdid_gmm_abort("Malformed `retained_index` / `dropped_index` partition in structural reduction certificate.", "cdid_gmm_reduction_error")
  }

  kept_idx <- if (!is.null(reduction$kept_indices)) as.integer(reduction$kept_indices) else ret_idx
  if (length(kept_idx) != L || anyNA(kept_idx) || anyDuplicated(kept_idx) || !setequal(kept_idx, ret_idx)) {
    cdid_gmm_abort("`kept_indices` must be a permutation of `certificate$retained_index`.", "cdid_gmm_reduction_error")
  }

  E_trip <- cert$E
  if (!is.list(E_trip) || !all(c("i", "j", "x") %in% names(E_trip))) {
    cdid_gmm_abort("`certificate$E` must be a sparse triplet list with fields `i`, `j`, `x`.", "cdid_gmm_reduction_error")
  }
  if (!valid_index(E_trip$i) || !valid_index(E_trip$j) || !is.numeric(E_trip$x) ||
      is.object(E_trip$x) || !is.null(dim(E_trip$x))) {
    cdid_gmm_abort("Malformed sparse reconstruction triplet types or indices.", "cdid_gmm_reduction_error")
  }
  if (!is.null(max_matrix_entries)) {
    gmm_safe_product(1, 6 * as.double(L)^2 + 12 * length(E_trip$i) + 40 * L,
      max_matrix_entries, "Gram matrices plus sparse reconstruction buffers")
  }
  E_i <- as.integer(E_trip$i)
  E_j <- as.integer(E_trip$j)
  E_x <- as.double(E_trip$x)
  n_trip <- length(E_i)
  if (length(E_j) != n_trip || length(E_x) != n_trip || anyNA(E_i) || anyNA(E_j) || anyNA(E_x) || any(!is.finite(E_x))) {
    cdid_gmm_abort("Sparse triplets `certificate$E` must have matching lengths and finite values.", "cdid_gmm_reduction_error")
  }
  if (!is.null(max_matrix_entries)) {
    gmm_safe_product(3, n_trip, max_matrix_entries, "sparse reduction triplets E")
  }
  if (n_trip > 0L && (!all(E_i %in% drop_idx) || !all(E_j %in% ret_idx))) {
    cdid_gmm_abort("Sparse triplets `certificate$E` must map dropped rows (`i`) to retained columns (`j`).", "cdid_gmm_reduction_error")
  }

  G <- diag(1.0, L)
  if (n_trip > 0L) {
    col_pos <- match(E_j, kept_idx)
    ord <- order(E_i, col_pos)
    ei_s <- E_i[ord]
    cj_s <- col_pos[ord]
    ex_s <- E_x[ord]
    starts <- which(c(TRUE, diff(ei_s) != 0L))
    ends <- c(starts[-1L] - 1L, n_trip)
    for (r_idx in seq_along(starts)) {
      rng <- starts[r_idx]:ends[r_idx]
      cols <- cj_s[rng]
      vals <- ex_s[rng]
      if (anyDuplicated(cols)) {
        u_cols <- unique(cols)
        vals <- vapply(u_cols, function(cc) sum(vals[cols == cc]), numeric(1L))
        cols <- u_cols
      }
      G[cols, cols] <- G[cols, cols] + tcrossprod(vals)
    }
  }
  for (j in seq_len(L)) for (i in seq_len(j)) G[j, i] <- G[i, j]

  C <- tryCatch(base::chol(G), error = function(e) NULL)
  if (is.null(C) || anyNA(C) || any(!is.finite(C)) || any(diag(C) <= 0.0)) {
    cdid_gmm_abort("Assembled reconstruction Gram matrix `G = E'E` is not positive definite.", "cdid_gmm_reduction_error")
  }
  is_id <- (M_full == L && n_trip == 0L)
  is_diag <- all(vapply(seq_len(L), function(j) all(G[-j, j] == 0), logical(1L)))
  list(
    G = G,
    C = C,
    M_full = M_full,
    L_retained = L,
    reduction_active = (M_full > L),
    is_identity_gram = is_id,
    is_diagonal_gram = is_diag,
    n_triplets = n_trip,
    kept_indices = kept_idx,
    gram_provenance = if (M_full > L) "sparse_structural_reduction_E_t_E" else "sparse_structural_certificate_identity"
  )
}

# ------------------------------------------------------------------------------
# Spectral condition trigger and roundoff-negative diagnosis on B = C V_r C'
# ------------------------------------------------------------------------------

diagnose_gmm_transformed_covariance <- function(V, gram_info, provenance, strategy, max_matrix_entries = NULL) {
  if (!is.matrix(V) || !is.numeric(V) || is.object(V) || nrow(V) < 1L || nrow(V) != ncol(V) ||
      anyNA(V) || any(!is.finite(V)) || !identical(unname(V), unname(t(V)))) {
    cdid_gmm_abort("Weight covariance must be a finite exactly symmetric numeric matrix.", "cdid_gmm_covariance_error")
  }
  L <- nrow(V)
  if (!is.null(max_matrix_entries)) {
    max_matrix_entries <- gmm_int_bound(max_matrix_entries, "max_matrix_entries")
    gmm_safe_product(1, 22 * gmm_safe_product(L, L, max_matrix_entries, "transformed covariance B") + 40 * L,
      max_matrix_entries, "Gram/covariance/transform/eigen/residual simultaneous buffers")
  }
  # Only the internally formed complete-unit-IF Gram, bound to these exact
  # covariance bytes and construction metadata using the existing seal, admits
  # the D27 numerical-negative diagnosis. Caller booleans/attributes do not.
  src <- provenance$complete_unit_if_gram
  gram_prov_flag <- identical(provenance$row_level, "unit") && is.list(src) &&
    exists("ipw_seal_ok", mode = "function") && ipw_seal_ok(src, "complete_unit_if_gram") &&
    identical(src$covariance, V) && identical(src$n_units, provenance$n_units) &&
    identical(src$nuisance_fit_id, provenance$nuisance_fit_id)
  if (!is.null(src) && !gram_prov_flag) {
    cdid_gmm_abort("Complete-unit-IF Gram construction provenance does not match this covariance.",
      c("cdid_gmm_pd_error", "cdid_gmm_covariance_error"))
  }

  B <- if (gram_info$is_identity_gram) {
    V
  } else {
    gram_info$C %*% V %*% t(gram_info$C)
  }
  if (anyNA(B) || any(!is.finite(B))) {
    cdid_gmm_abort(
      "Transformed covariance `B = C V C'` contains nonfinite values.",
      c("cdid_gmm_pd_error", "cdid_gmm_covariance_error")
    )
  }
  for (j in seq_len(L)) for (i in seq_len(j)) B[j, i] <- B[i, j]

  is_zero_cov <- all(V == 0.0) && all(B == 0.0)
  eig <- tryCatch(eigen(B, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig) || anyNA(eig$values) || any(!is.finite(eig$values)) ||
      anyNA(eig$vectors) || any(!is.finite(eig$vectors))) {
    cdid_gmm_abort(
      "Symmetric eigendecomposition of transformed covariance `B = C V C'` failed.",
      c("cdid_gmm_pd_error", "cdid_gmm_covariance_error")
    )
  }
  d <- as.double(eig$values)
  U <- eig$vectors
  if (typeof(U) != "double") storage.mode(U) <- "double"
  d_max <- max(d)
  d_min <- min(d)

  norm_B_F <- gmm_scale_safe_frobenius(B)
  b_bound <- 100.0 * as.double(L) * .Machine$double.eps * norm_B_F
  if (is_zero_cov) {
    r_eig <- 0.0
  } else {
    U_D_Ut <- U %*% sweep(t(U), 1L, d, `*`)
    r_eig <- gmm_scale_safe_frobenius(B - U_D_Ut)
  }
  eigensystem_refinement <- list(method = "none", original_residual = r_eig)

  # Some symmetric LAPACK paths lose a little orthogonality even for a
  # well-conditioned B. Reorthogonalize once and recompute Rayleigh values;
  # accept only if the original, unchanged reconstruction bound is met.
  # B itself, the ridge policy and all covariance definiteness guards stay intact.
  if (!is_zero_cov && is.finite(b_bound) && is.finite(r_eig) && r_eig > b_bound) {
    U_D_Ut <- NULL
    qU <- tryCatch(qr(U), error = function(e) NULL)
    if (!is.null(qU) && qU$rank == L) {
      U_refined <- qr.Q(qU)
      qU <- NULL
      d_refined <- colSums(U_refined * (B %*% U_refined))
      ord <- order(d_refined, decreasing = TRUE)
      d_refined <- d_refined[ord]
      U_refined <- U_refined[, ord, drop = FALSE]
      r_refined <- gmm_scale_safe_frobenius(
        B - U_refined %*% sweep(t(U_refined), 1L, d_refined, `*`))
      if (all(is.finite(d_refined)) && is.finite(r_refined) && r_refined <= b_bound) {
        eigensystem_refinement$method <- "QR_Rayleigh_reorthogonalization"
        eigensystem_refinement$refined_residual <- r_refined
        U <- U_refined
        d <- d_refined
        d_min <- min(d)
        d_max <- max(d)
        r_eig <- r_refined
      }
    }
  }

  if (!is_zero_cov) {
    if (!is.finite(b_bound) || !is.finite(r_eig) || r_eig > b_bound) {
      cdid_gmm_abort(
        "Transformed covariance eigendecomposition residual `r = ||B - U D U'||_F` exceeds bound `b = 100 * L * eps * ||B||_F`.",
        c("cdid_gmm_pd_error", "cdid_gmm_covariance_error"),
        diagnostics = list(r_eigen_residual = r_eig, b_roundoff_bound = b_bound, d_min = d_min, d_max = d_max)
      )
    }
    if (d_min < -(b_bound + r_eig)) {
      cdid_gmm_abort(
        "Weight covariance `V` is materially indefinite: minimum eigenvalue `d_min < -(b + r)`.",
        c("cdid_gmm_pd_error", "cdid_gmm_covariance_error"),
        diagnostics = list(d_min = d_min, d_max = d_max, b_roundoff_bound = b_bound, r_eigen_residual = r_eig)
      )
    }
    if (any(diag(B) < -(b_bound + r_eig))) {
      cdid_gmm_abort(
        "Transformed covariance `B` has a materially negative diagonal entry.",
        c("cdid_gmm_pd_error", "cdid_gmm_covariance_error"),
        diagnostics = list(min_diag_B = min(diag(B)), b_roundoff_bound = b_bound, r_eigen_residual = r_eig)
      )
    }
    if (d_min < 0.0) {
      if (!isTRUE(gram_prov_flag)) {
        cdid_gmm_abort(
          "Negative eigenvalue in `B` is refused because certified unit-IF Gram provenance is absent.",
          c("cdid_gmm_pd_error", "cdid_gmm_covariance_error"),
          diagnostics = list(d_min = d_min, b_roundoff_bound = b_bound, r_eigen_residual = r_eig, gram_provenance_certified = FALSE)
        )
      }
      neg_diag <- "consistent_with_roundoff"
      diag_num_neg <- TRUE
    } else {
      neg_diag <- "nonnegative_spectrum"
      diag_num_neg <- FALSE
    }
    trace_B <- sum(diag(B))
    scale_s <- trace_B / as.double(gram_info$M_full)
    if (!is.finite(trace_B) || trace_B <= 0.0 || !is.finite(scale_s) || scale_s <= 0.0 || d_max <= 0.0) {
      cdid_gmm_abort(
        "Transformed covariance `B` has nonpositive trace or maximum eigenvalue outside the exact-zero covariance branch.",
        c("cdid_gmm_pd_error", "cdid_gmm_covariance_error"),
        diagnostics = list(trace_B = trace_B, scale_s = scale_s, d_max = d_max)
      )
    }
  } else {
    trace_B <- 0.0
    scale_s <- 0.0
    neg_diag <- "exact_zero_matrix"
    diag_num_neg <- FALSE
  }

  rho <- if (is_zero_cov || d_max <= 0.0) 0.0 else max(d_min, 0.0) / d_max
  rho_raw <- if (is_zero_cov || d_max == 0.0) 0.0 else d_min / d_max
  cond_trig <- strategy$condition_trigger
  triggered <- is_zero_cov || d_min <= 0.0 || rho < cond_trig
  trigger_reason <- if (is_zero_cov) {
    "zero_covariance"
  } else if (d_min < 0.0) {
    "roundoff_negative_psd"
  } else if (d_min == 0.0) {
    "exact_zero_eigenvalue"
  } else if (rho < cond_trig) {
    "spectral_condition_below_trigger"
  } else {
    "well_conditioned_pd"
  }

  list(
    B = B,
    eig = list(values = d, vectors = U),
    is_zero_cov = is_zero_cov,
    d_min = d_min,
    d_max = d_max,
    rho = rho,
    rho_raw = rho_raw,
    norm_B_F = norm_B_F,
    b_bound = b_bound,
    r_eig = r_eig,
    eigensystem_refinement = eigensystem_refinement,
    trace_B = trace_B,
    scale_s = scale_s,
    gram_prov_flag = gram_prov_flag,
    negative_eigenvalue_diagnosis = neg_diag,
    diagnosed_numerical_negative = diag_num_neg,
    nonpositive_eigenvalue_count = as.integer(sum(d <= 0.0)),
    triggered = triggered,
    trigger_reason = trigger_reason,
    covariance_rank = as.integer(sum(d > b_bound + r_eig)),
    covariance_rank_rule = "computed eigenvalues above b + r; numerical diagnostic only"
  )
}

# ------------------------------------------------------------------------------
# Handler execution: covariance-ridge Tikhonov
# ------------------------------------------------------------------------------

solve_gmm_covariance_handler <- function(
  V, mhat_R, W_G, gauge, target_info, gram_info, diag_B, strategy, accuracy,
  max_matrix_entries = NULL
) {
  m_R <- gauge$n_links
  k <- gauge$rank_k
  design_est_idx <- target_info$estimable_indices
  J_design_est <- length(design_est_idx)
  C_gauge_design <- target_info$C_gauge
  nb_info <- list(backend = "base_R_BLAS_LAPACK_computed_diagnostics")
  if (!is.null(max_matrix_entries)) {
    gmm_safe_product(1, 30 * as.double(m_R)^2 + 6 * m_R * k +
      8 * m_R * J_design_est + 6 * k * J_design_est + 4 * as.double(k)^2 +
      40 * m_R + 16 * (m_R + k + J_design_est), max_matrix_entries,
      "covariance handler retained matrices, factor/solve/loading/residual buffers")
  }

  if (!identical(accuracy$method, "conventional_chol_qr")) {
    cdid_gmm_abort(
      "Private interval numerical certification (`certified_eigen_chol_kkt`) requires an unregularized strictly positive-definite covariance; regularized certificates are not fabricated.",
      "cdid_gmm_pd_error"
    )
  }

  if (identical(strategy$method, "tikhonov")) {
    if (diag_B$is_zero_cov) {
      lambda <- 0.0
      R_B_lam <- gram_info$C
      rc_R <- rcond(R_B_lam, triangular = TRUE)
      factor_rel <- if (gram_info$is_identity_gram) 0.0 else gmm_scale_safe_frobenius(crossprod(gram_info$C) - gram_info$G) / gmm_scale_safe_frobenius(gram_info$G)
      A_whitened <- if (gram_info$is_identity_gram) W_G else gram_info$C %*% W_G
      b_whitened <- if (gram_info$is_identity_gram) mhat_R else as.double(gram_info$C %*% mhat_R)
      whitening_rel <- 0.0
    } else {
      lambda <- strategy$ridge_relative * diag_B$scale_s
      if (!is.finite(lambda) || lambda <= 0.0) {
        cdid_gmm_abort("Computed covariance-ridge penalty `lambda` is nonpositive or nonfinite.", "cdid_gmm_pd_error")
      }
      B_lam <- diag_B$B + diag(lambda, m_R)
      R_B_lam <- tryCatch(base::chol(B_lam), error = function(e) NULL)
      if (is.null(R_B_lam) || anyNA(R_B_lam) || any(!is.finite(R_B_lam)) || any(diag(R_B_lam) <= 0.0)) {
        cdid_gmm_abort(
          "Cholesky factorization of regularized covariance `B + lambda * I` failed; automatic lambda increase or method fallback is forbidden.",
          "cdid_gmm_pd_error",
          diagnostics = list(lambda = lambda, scale = diag_B$scale_s, d_min = diag_B$d_min, d_max = diag_B$d_max)
        )
      }
      rc_R <- rcond(R_B_lam, triangular = TRUE)
      factor_rel <- gmm_scale_safe_frobenius(crossprod(R_B_lam) - B_lam) / gmm_scale_safe_frobenius(B_lam)
      if (!is.finite(rc_R) || rc_R < min(accuracy$min_rcond_R, sqrt(strategy$condition_trigger)) ||
          !is.finite(factor_rel) || factor_rel > accuracy$max_factor_residual_rel) {
        cdid_gmm_abort(
          "Regularized covariance Cholesky conditioning or factor residual refused.",
          "cdid_gmm_pd_error",
          diagnostics = list(rcond_R = rc_R, factor_rel = factor_rel, lambda = lambda)
        )
      }
      CW <- if (gram_info$is_identity_gram) W_G else gram_info$C %*% W_G
      Cm <- if (gram_info$is_identity_gram) mhat_R else as.double(gram_info$C %*% mhat_R)
      A_whitened <- backsolve(R_B_lam, CW, transpose = TRUE)
      b_whitened <- backsolve(R_B_lam, Cm, transpose = TRUE)
      if (any(!is.finite(A_whitened)) || any(!is.finite(b_whitened))) {
        cdid_gmm_abort("Nonfinite values in covariance-ridge whitened system.", "cdid_gmm_pd_error")
      }
      whitening_rel <- gmm_scale_safe_frobenius(crossprod(R_B_lam, A_whitened) - CW) / gmm_scale_safe_frobenius(CW)
      if (!is.finite(whitening_rel) || whitening_rel > accuracy$max_factor_residual_rel) {
        cdid_gmm_abort("Covariance-ridge whitened design residual exceeds caller bound.", "cdid_gmm_qr_error")
      }
    }

    qa <- base::qr(A_whitened, LAPACK = TRUE)
    U_qr <- base::qr.R(qa)
    piv <- qa$pivot
    if (anyNA(U_qr) || any(!is.finite(U_qr)) || any(diag(U_qr) == 0.0)) {
      cdid_gmm_abort("LAPACK QR factor U has a zero or nonfinite pivot under covariance ridge.", "cdid_gmm_qr_error")
    }
    rc_U <- rcond(U_qr, triangular = TRUE)
    gram_A <- crossprod(A_whitened)
    if (any(!is.finite(gram_A))) cdid_gmm_abort("Gram candidate overflowed under covariance ridge.", "cdid_gmm_qr_error")
    ea <- eigen(gram_A, symmetric = TRUE, only.values = TRUE)$values
    sigma_rel <- sqrt(max(0.0, min(ea)) / max(ea))
    if (any(!is.finite(ea)) || min(ea) <= 0.0 || !is.finite(sigma_rel) || sigma_rel < accuracy$min_sigma_rel_A ||
        !is.finite(rc_U) || rc_U < accuracy$min_rcond_U) {
      cdid_gmm_abort("Whitened design / QR conditioning refused under covariance ridge.", "cdid_gmm_qr_error")
    }

    qtb <- base::qr.qty(qa, b_whitened)
    theta_keep <- numeric(k)
    theta_keep[piv] <- backsolve(U_qr, qtb[seq_len(k)])
    if (any(!is.finite(theta_keep))) {
      cdid_gmm_abort("Nonfinite parameter solution `theta` in covariance-ridge triangular solve.", "cdid_gmm_qr_error")
    }
    names(theta_keep) <- gauge$keep_param_ids

    if (J_design_est > 0L) {
      C_piv <- C_gauge_design[piv, , drop = FALSE]
      Z1 <- backsolve(U_qr, C_piv, transpose = TRUE)
      Z2 <- base::qr.qy(qa, rbind(Z1, matrix(0.0, nrow = m_R - k, ncol = J_design_est)))
      Z3 <- if (diag_B$is_zero_cov) Z2 else backsolve(R_B_lam, Z2)
      Lambda <- if (gram_info$is_identity_gram) Z3 else crossprod(gram_info$C, Z3)
      Mc <- matrix(0.0, nrow = k, ncol = J_design_est, dimnames = list(gauge$keep_param_ids, colnames(C_gauge_design)))
      Mc[piv, ] <- backsolve(U_qr, Z1)
      rownames(Lambda) <- gauge$link_ids
      colnames(Lambda) <- colnames(C_gauge_design)
      colnames(Z1) <- colnames(C_gauge_design)
      if (any(!is.finite(Z1)) || any(!is.finite(Lambda)) || any(!is.finite(Mc))) {
        cdid_gmm_abort("Nonfinite target loadings `Lambda` or `Mc` under covariance ridge.", "cdid_gmm_qr_error")
      }
    } else {
      Z1 <- matrix(0.0, nrow = k, ncol = 0L)
      Lambda <- matrix(0.0, nrow = m_R, ncol = 0L, dimnames = list(gauge$link_ids, character(0)))
      Mc <- matrix(0.0, nrow = k, ncol = 0L, dimnames = list(gauge$keep_param_ids, character(0)))
    }

    e_resid <- mhat_R - as.double(W_G %*% theta_keep)
    if (diag_B$is_zero_cov) {
      y_mult <- if (gram_info$is_identity_gram) e_resid else as.double(crossprod(gram_info$C, gram_info$C %*% e_resid))
    } else {
      y_trans <- backsolve(R_B_lam, backsolve(R_B_lam, if (gram_info$is_identity_gram) e_resid else as.double(gram_info$C %*% e_resid), transpose = TRUE))
      y_mult <- if (gram_info$is_identity_gram) as.double(y_trans) else as.double(crossprod(gram_info$C, y_trans))
    }
    if (any(!is.finite(y_mult))) {
      cdid_gmm_abort("Nonfinite KKT multiplier `y` under covariance ridge.", "cdid_gmm_qr_error")
    }

    Vcrit_apply <- function(X) {
      GX <- if (gram_info$is_identity_gram) X else backsolve(gram_info$C, backsolve(gram_info$C, X, transpose = TRUE))
      if (diag_B$is_zero_cov) GX else (V %*% X + lambda * GX)
    }
    Vcrit_y <- drop(Vcrit_apply(y_mult))
    r1 <- mhat_R - Vcrit_y - drop(W_G %*% theta_keep)
    r2 <- -drop(crossprod(W_G, y_mult))
    s1 <- if (J_design_est > 0L) W_G %*% Mc - Vcrit_apply(Lambda) else matrix(0.0, m_R, 0L)
    s2 <- if (J_design_est > 0L) C_gauge_design - crossprod(W_G, Lambda) else matrix(0.0, k, 0L)
    tau <- if (J_design_est > 0L) drop(crossprod(C_gauge_design, theta_keep)) else numeric(0L)
    scale1 <- gmm_scale_safe_frobenius(Vcrit_y) + gmm_scale_safe_frobenius(W_G) * gmm_scale_safe_frobenius(theta_keep) + gmm_scale_safe_frobenius(mhat_R)
    scale2 <- gmm_scale_safe_frobenius(W_G) * gmm_scale_safe_frobenius(y_mult)
    ratio <- function(r, sc) if (sc > 0.0) gmm_scale_safe_frobenius(r) / sc else if (all(r == 0.0)) 0.0 else Inf
    if (any(!is.finite(c(r1, r2, s1, s2, tau)))) {
      cdid_gmm_abort("Nonfinite solve/loading diagnostics under covariance ridge.", "cdid_gmm_accuracy_error")
    }

    guard <- list(
      certified = FALSE, point_certified = FALSE, loading_certified = FALSE,
      loading_budget_supplied = FALSE, loading_validation_status = "not_requested_D23",
      point_validation_status = "not_requested_D23", complete_inference_certified = FALSE,
      full_inference_status = "complete_per_fit_error_certification_not_requested_or_performed",
      backend = nb_info$backend, native_backend = nb_info, method = accuracy$method, weighting = "optimal",
      covariance_handler = "tikhonov",
      computed_target_values = tau, propagation_bounds = NULL,
      conditioning = list(eigenvalues_B = diag_B$eig$values,
        eigenvalues_B_regularized = if (diag_B$is_zero_cov) NULL else diag_B$eig$values + lambda,
        factor_role = if (diag_B$is_zero_cov) "reconstruction_Gram_C" else "shifted_B_Cholesky",
        rcond_R_1norm = rc_R, rcond_U_1norm = rc_U, sigma_rel_A_est_2norm = sigma_rel),
      factorization = list(cholesky_relative_residual = factor_rel, whitening_relative_residual = whitening_rel),
      primal_kkt = list(r1_est = r1, r2_est = r2, solve_rel_est = ratio(r1, scale1), stat_rel_est = ratio(r2, scale2),
        residual_ratio_role = "diagnostic_only_under_D23"),
      targets = list(target_error_est = rep(NA_real_, J_design_est), target_error_cert = rep(NA_real_, J_design_est),
        loading_error_cert = rep(NA_real_, J_design_est), se_criterion_est = if (J_design_est > 0L) sqrt(colSums(Z1 * Z1)) else numeric(0L),
        loading_primal_residual = s1, loading_calibration_residual = s2)
    )
    raw_quad <- if (m_R > k) sum(qtb[-seq_len(k)]^2) else 0.0
    raw_Z1tZ1 <- if (J_design_est > 0L) crossprod(Z1) else matrix(0.0, 0L, 0L)
    if (!is.finite(raw_quad) || any(!is.finite(raw_Z1tZ1)) ||
        any(!is.finite(guard$targets$se_criterion_est))) {
      cdid_gmm_abort("Covariance-ridge criterion covariance or quadratic residual overflowed.", "cdid_gmm_accuracy_error")
    }

    return(list(
      handler_used = "tikhonov",
      lambda = lambda,
      weighted_design_rank = k,
      qa = qa,
      U_qr = U_qr,
      piv = piv,
      R_chol = R_B_lam,
      A_whitened = A_whitened,
      b_whitened = b_whitened,
      theta_keep = theta_keep,
      Z1 = Z1,
      Lambda = Lambda,
      Mc = Mc,
      y_mult = y_mult,
      guard = guard,
      raw_Z1tZ1 = raw_Z1tZ1,
      raw_quad = raw_quad,
      target_info = target_info
    ))
  }

}

build_gmm_covariance_handling_record <- function(
  weighting, strategy, gram_info, diag_B, handler_used, trigger_reason,
  lambda, weighted_design_rank, gauge, target_info, provenance
) {
  M_pos <- provenance$n_positive_contributing_rows
  max_rank_allowed <- if (provenance$centered) M_pos - 1L else M_pos
  w_req <- if (!is.null(strategy)) strategy$weighting_requested else weighting
  w_can <- if (!is.null(strategy)) strategy$weighting_canonical else weighting
  h_sel <- if (identical(w_can, "identity") || is.null(strategy)) "none" else strategy$method
  h_inv <- !identical(handler_used, "none")

  list(
    weighting_requested = w_req,
    weighting_canonical = w_can,
    alias_used = identical(w_req, "two_step"),
    handler_selected = h_sel,
    handler_invoked = h_inv,
    handler_used = handler_used,
    trigger_reason = trigger_reason,
    dimensions = list(
      M_full = as.integer(gram_info$M_full),
      L_retained = as.integer(gram_info$L_retained),
      p_parameters = as.integer(gauge$n_parameters),
      k_design_rank = as.integer(gauge$rank_k),
      n_targets = as.integer(target_info$n_targets),
      n_estimable = as.integer(sum(target_info$estimable)),
      M_positive_rows = as.integer(M_pos),
      sample_rank_upper_bound = as.integer(max_rank_allowed),
      exceeds_sample_rank_bound = (identical(w_can, "optimal") && gauge$n_links > max_rank_allowed)
    ),
    gram_provenance = gram_info$gram_provenance,
    gram_is_identity = gram_info$is_identity_gram,
    gram_is_diagonal = gram_info$is_diagonal_gram,
    gram_n_triplets = as.integer(gram_info$n_triplets),
    kept_indices = as.integer(gram_info$kept_indices),
    condition_trigger = if (!is.null(strategy)) strategy$condition_trigger else NA_real_,
    ridge_relative = if (!is.null(strategy)) strategy$ridge_relative else NA_real_,
    scale = if (!is.null(diag_B)) diag_B$scale_s else NA_real_,
    trace_B = if (!is.null(diag_B)) diag_B$trace_B else NA_real_,
    lambda = lambda,
    covariance_rank = if (!is.null(diag_B)) diag_B$covariance_rank else NA_integer_,
    covariance_rank_rule = if (!is.null(diag_B)) diag_B$covariance_rank_rule else "not_computed",
    weighted_design_rank = as.integer(weighted_design_rank),
    conditioning = if (!is.null(diag_B)) list(
      rho = diag_B$rho,
      rho_raw = diag_B$rho_raw,
      d_min = diag_B$d_min,
      d_max = diag_B$d_max,
      eigenvalues_B = diag_B$eig$values,
      norm_B_frobenius = diag_B$norm_B_F,
      b_roundoff_bound = diag_B$b_bound,
      r_eigen_residual = diag_B$r_eig,
      eigensystem_refinement = diag_B$eigensystem_refinement,
      diagnosed_numerical_negative = diag_B$diagnosed_numerical_negative,
      negative_eigenvalue_diagnosis = diag_B$negative_eigenvalue_diagnosis,
      nonpositive_eigenvalue_count = diag_B$nonpositive_eigenvalue_count,
      gram_provenance_certified = diag_B$gram_prov_flag
    ) else NULL,
    zero_covariance = if (!is.null(diag_B)) diag_B$is_zero_cov else FALSE,
    zero_covariance_branch = if (!is.null(diag_B) && diag_B$is_zero_cov) {
      "scale_free_gram_limit"
    } else "not_applicable",
    inference_claim = if (!is.null(diag_B) && diag_B$is_zero_cov) "none_zero_sample_covariance" else "complete_unit_if_sandwich_conditional_on_stable_loading",
    target_availability = data.frame(
      target_id = target_info$target_ids,
      estimable = target_info$estimable,
      status = target_info$status,
      stringsAsFactors = FALSE
    ),
    tuning_overrides = if (!is.null(strategy)) strategy$overrides else character(0L),
    scope_note = paste0(
      "Regularized criterion covariance is used only for the GMM weighting matrix; ",
      "inference replays the complete unit target IF and forms the actual sandwich covariance. ",
      "Feasible two-step weights estimate the optimal weighting matrix; no finite-sample exact-unbiasedness ",
      "or unconditional inference claim is made."
    )
  )
}
