# D21 full declared nonreference pretrend null. D23 public computation uses
# ordinary Cholesky/solve with the existing factor conditioning threshold.
# Required missing cells, nonfinite inputs or unusable full covariance make
# the joint test unavailable; effects remain available. No restriction dropping,
# tiny-SE filtering or pseudoinverse. The explicit private enclosure path below
# remains research material and does not supply a public per-fit guarantee.

chained_pretest_abort <- function(message, class = character()) {
  cdid_abort(message, unique(c(class, "cdid_chained_pretest_error", "cdid_chained_request_error")))
}

chained_pretest_not_conditional <- "this is the unconditional full-set pretrend Wald test (did-style W/Wpval); it is never the conditional pretrend test, which is the separate cdid::conditional_did_pretest()"

# Declared nonreference pretrend set of a reporting target_meta table.
chained_pretest_declared_set <- function(target_meta) {
  if (!is.data.frame(target_meta) || !all(c("role", "target_id", "available") %in% names(target_meta)))
    chained_pretest_abort("`target_meta` must carry `role`, `target_id` and `available`.", "cdid_chained_pretest_input_error")
  pre <- which(target_meta$role == "pre")
  ref <- which(target_meta$role == "reference")
  list(pre = pre, reference = ref, ids = as.character(target_meta$target_id[pre]),
       available = as.logical(target_meta$available[pre]))
}

# Core computation on explicit inputs: point vector/covariance and their
# certified enclosures (matrices lo/hi of the same shape).
chained_pretest_wald <- function(theta, V, theta_lo = NULL, theta_hi = NULL, V_lo = NULL, V_hi = NULL, ids, relative_accuracy = NULL, max_matrix_entries = NULL) {
  if (!is.null(relative_accuracy) && (!is.numeric(relative_accuracy) || length(relative_accuracy) != 1L || !is.finite(relative_accuracy) ||
      relative_accuracy <= 0 || relative_accuracy >= 1))
    chained_pretest_abort("`relative_accuracy` must be a positive finite level below 1 (D20).", "cdid_chained_pretest_input_error")
  J <- length(theta)
  unavailable <- function(reason, extra = list()) c(list(status = "unavailable", reason = reason, df = J, cells = ids,
    W = NULL, Wpval = NULL, not_conditional_pretest = chained_pretest_not_conditional,
    identification_claim = "none: covariance or numerical failure of this joint test is not a statement that any ATT is unidentified"), extra)
  if (J < 1L) return(list(status = "not_applicable", reason = "no declared nonreference pretrend comparison exists in the selected base-period reporting",
    df = 0L, cells = character(0L), W = NULL, Wpval = NULL, not_conditional_pretest = chained_pretest_not_conditional))
  a <- as.double(theta); A <- unname(as.matrix(V))
  if (!is.matrix(A) || !identical(dim(A), c(J, J)) || !is.numeric(A) || anyNA(a) || anyNA(A) || any(!is.finite(a)) || any(!is.finite(A)))
    return(unavailable("the declared pretrend estimates or their joint covariance are not finite"))
  if (is.null(relative_accuracy)) {
    # Conventional full-set Wald. The old enclosure and its width are not a
    # public acceptance predicate under D23. Reuse the Cholesky conditioning
    # policy; do not choose another rank tolerance.
    gmm_safe_product(4, J*J, max_matrix_entries, "pretrend covariance/solve diagnostics")
    if (!all(A == t(A))) return(unavailable("the joint pretrend covariance is not symmetric"))
    F <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(F) || any(!is.finite(F)) || any(diag(F) <= 0))
      return(unavailable("the full joint pretrend covariance has no usable Cholesky factor; no pseudoinverse or restriction removal is applied"))
    rc <- rcond(F, triangular = TRUE)
    if (!is.finite(rc) || rc < chained_public_carried_settings()$gmm_guard$min_rcond_R)
      return(unavailable("the full joint pretrend covariance fails the existing Cholesky conditioning policy; no restrictions are removed"))
    x <- tryCatch(as.double(backsolve(F, backsolve(F, a, transpose = TRUE))), error = function(e) NULL)
    if (is.null(x) || any(!is.finite(x))) return(unavailable("the full joint pretrend covariance solve is unusable"))
    W <- sum(a*x); residual <- drop(A %*% x) - a
    p <- stats::pchisq(W, df = J, lower.tail = FALSE)
    if (!is.finite(W) || W < 0 || !is.finite(p) || any(!is.finite(residual)))
      return(unavailable("the full-set Wald computation is not finite and nonnegative"))
    return(list(status = "computed", W = W, Wpval = p, df = J, cells = ids, certified = FALSE,
      certification_status = "complete_per_fit_error_certification_not_requested_or_performed",
      covariance = list(rcond_cholesky = rc), solve_residual = residual,
      description = "full declared nonreference pretrend Wald test under conventional computation; structural reference zeros excluded as identities",
      not_conditional_pretest = chained_pretest_not_conditional))
  }
  tlo <- as.double(theta_lo); thi <- as.double(theta_hi); Vlo <- unname(as.matrix(V_lo)); Vhi <- unname(as.matrix(V_hi))
  if (length(tlo) != J || length(thi) != J || !identical(dim(Vlo), c(J, J)) || !identical(dim(Vhi), c(J, J)) ||
      anyNA(tlo) || anyNA(thi) || anyNA(Vlo) || anyNA(Vhi) || any(!is.finite(c(tlo, thi, Vlo, Vhi))) || any(tlo > a) || any(thi < a) || any(Vlo > A) || any(Vhi < A) || any(Vlo > Vhi))
    return(unavailable("the certified enclosures of the declared pretrend estimates/covariance are missing or do not contain the reported values"))
  mme <- max_matrix_entries
  compute <- function() {
    # Symmetric candidate. The exact covariance is symmetric and lies in the
    # box, so the box distance of the symmetrized candidate bounds its 2-norm
    # error (Frobenius >= spectral norm).
    A_sym <- if (all(A == t(A))) A else 0.5 * (A + t(A))
    dev_V <- ipw_pipeline_interval_distance(A_sym, Vlo, Vhi)
    eps_V <- gmm_iv_norm(unname(dev_V), max_matrix_entries = mme)[2L]
    dev_t <- ipw_pipeline_interval_distance(a, tlo, thi)
    eps_t <- gmm_iv_norm(unname(dev_t), max_matrix_entries = mme)[2L]
    eg <- certify_symmetric_eigen_bounds(A_sym, "declared pretrend covariance", max_matrix_entries = mme)
    l <- eg$lambda_min_cert
    if (!isTRUE(eg$valid) || !is.finite(l) || l <= 0) {
      return(unavailable(paste0("the joint covariance of the declared pretrend comparisons is not certified positive definite (eigen-residual certificate: ",
        "estimated lambda_min ", format(eg$lambda_min_est, digits = 6), ", certified lower bound ", format(l, digits = 6),
        "); no pseudoinverse, restriction removal or tiny-SE filter is applied"),
        list(lambda_min_est = eg$lambda_min_est, lambda_min_cert = l, eps_V = eps_V)))
    }
    l_star <- gmm_iv_sub(gmm_iv(l), gmm_iv(eps_V))[1L]
    if (!is.finite(l_star) || l_star <= 0) {
      return(unavailable(paste0("the certified lower eigenvalue of the declared pretrend covariance (", format(l, digits = 6),
        ") does not exceed its certified enclosure width (", format(eps_V, digits = 6), "); the exact covariance is not certified nonsingular"),
        list(lambda_min_est = eg$lambda_min_est, lambda_min_cert = l, eps_V = eps_V)))
    }
    # Candidate solve and its certified residual: ||x* - x|| <= ||r|| / l.
    x <- tryCatch(as.double(solve(A_sym, a)), error = function(e)
      chained_pretest_abort(paste0("the candidate Wald solve failed: ", conditionMessage(e)), "cdid_chained_pretest_solve_error"))
    res <- gmm_iv_mat_residual(A_sym, matrix(x, ncol = 1L), matrix(a, ncol = 1L), max_matrix_entries = mme)
    r_norm <- res$norm[2L]
    a_norm <- gmm_iv_norm(a, max_matrix_entries = mme)
    W_iv <- gmm_iv_dot(a, x, max_matrix_entries = mme)
    # (i) candidate quadratic form: |a'x* - a'x| <= ||a|| ||r|| / l
    e_solve <- gmm_iv_div(gmm_iv_mul(gmm_iv(a_norm[2L]), gmm_iv(r_norm)), gmm_iv(l))[2L]
    # (ii) exact inputs within their boxes (Proposition in the note):
    #      eps_t (2||a|| + eps_t) / l_star + ||a||^2 eps_V / (l l_star)
    term1 <- gmm_iv_div(gmm_iv_mul(gmm_iv(eps_t), gmm_iv_add(gmm_iv_mul(gmm_iv(2.0), gmm_iv(a_norm[2L])), gmm_iv(eps_t))), gmm_iv(l_star))[2L]
    term2 <- gmm_iv_div(gmm_iv_mul(gmm_iv_mul(gmm_iv(a_norm[2L]), gmm_iv(a_norm[2L])), gmm_iv(eps_V)),
                        gmm_iv_mul(gmm_iv(l), gmm_iv(l_star)))[2L]
    E <- gmm_iv_add(gmm_iv(e_solve), gmm_iv_add(gmm_iv(term1), gmm_iv(term2)))[2L]
    W_lo <- max(0.0, gmm_iv_sub(gmm_iv(W_iv[1L]), gmm_iv(E))[1L])
    W_hi <- gmm_iv_add(gmm_iv(W_iv[2L]), gmm_iv(E))[2L]
    W_hat <- sum(a * x)
    if (!is.finite(W_hat) || W_hat < 0) return(unavailable("the computed Wald statistic is not a finite nonnegative number"))
    # Finite arithmetic is required, but D21 adopts no additional Wald scale
    # or p-value accuracy cutoff. Keep the valid enclosure as a diagnostic.
    width <- gmm_iv_sub(gmm_iv(W_hi), gmm_iv(W_lo))[2L]
    if (any(!is.finite(c(W_lo, W_hi, E, width))) || W_lo > W_hi)
      return(unavailable("the Wald arithmetic has no finite ordered enclosure"))
    p_hat <- stats::pchisq(W_hat, df = J, lower.tail = FALSE)
    p_int <- c(stats::pchisq(W_hi, df = J, lower.tail = FALSE), stats::pchisq(W_lo, df = J, lower.tail = FALSE))
    list(status = "computed", W = W_hat, Wpval = p_hat, df = J, cells = ids,
      W_interval = c(W_lo, W_hi), p_interval_evaluated = p_int, certified = TRUE,
      certification_scope = "finite Wald arithmetic enclosure only; pchisq values are ordinary numerical evaluations, not certified CDF bounds; no additional Wald/p-value accuracy cutoff",
      arithmetic = list(width = width, error_bound = E),
      covariance = list(lambda_min_est = eg$lambda_min_est, lambda_min_cert = l, lambda_min_exact_lower = l_star, eps_V = eps_V,
        eigen_orthogonality_defect_cert = eg$delta_Q_cert, symmetrized = !all(A == t(A))),
      enclosure = list(eps_theta = eps_t, solve_residual_norm = r_norm, error_solve = e_solve, error_inputs = gmm_iv_add(gmm_iv(term1), gmm_iv(term2))[2L]),
      description = paste0("full declared nonreference pretrend Wald test: W = theta' V^-1 theta over ", J, " declared comparison",
        if (J == 1L) "" else "s", ", joint covariance certified positive definite, arithmetic enclosed with outward primitives; p-value evaluated by ordinary stats::pchisq; ",
        "exact structural reference zeros excluded as identities; ", chained_pretest_not_conditional),
      not_conditional_pretest = chained_pretest_not_conditional)
  }
  tryCatch(compute(), cdid_error = function(e) unavailable(paste0("the Wald arithmetic could not be certified: ", conditionMessage(e)),
    list(error_class = class(e)[1L])))
}

# From the pipeline reporting list (target_meta, estimates, inference$vcov,
# numerical_certificate$theta/$vcov) to the D21 record carried on the cdid_MP.
chained_pretest_record <- function(reporting, relative_accuracy = NULL, max_matrix_entries = NULL) {
  tm <- reporting$target_meta
  set <- chained_pretest_declared_set(tm)
  base <- list(decision = "D21", base_period = reporting$base_period, declared_cells = set$ids,
    structural_reference_cells = as.character(tm$target_id[set$reference]),
    n_declared_pretrend = length(set$pre), n_available_pretrend = sum(set$available),
    null = "every declared nonreference pretrend comparison of the selected base-period reporting equals zero; exact structural reference zeros are identities, not restrictions")
  if (!length(set$pre)) {
    rec <- chained_pretest_wald(numeric(0L), matrix(0, 0L, 0L), numeric(0L), numeric(0L), matrix(0, 0L, 0L), matrix(0, 0L, 0L), character(0L), relative_accuracy)
    return(c(base, rec))
  }
  if (!all(set$available)) {
    miss <- set$ids[!set$available]
    return(c(base, list(status = "unavailable", df = length(set$pre), cells = set$ids, W = NULL, Wpval = NULL,
      reason = paste0("declared pretrend comparison", if (length(miss) == 1L) " " else "s ", paste(miss, collapse = ", "),
        if (length(miss) == 1L) " is" else " are", " unavailable (", paste(unique(tm$support_status[match(miss, tm$target_id)]), collapse = ", "),
        "); the full declared null cannot be tested and no available-only subset is substituted; available effects and their inference are retained"),
      unavailable_cells = miss, not_conditional_pretest = chained_pretest_not_conditional,
      identification_claim = "none: an unavailable pretrend comparison is catalog-relative support information, not a statement about ATT identification")))
  }
  cert <- reporting$numerical_certificate
  # Centered unit/cluster influence sums span at most M-1 dimensions. A
  # rounded Gram can acquire tiny positive pivots, so Cholesky alone must
  # not authorize a full-null test that violates this exact rank bound.
  inf <- reporting$inference
  n_rows <- if (isTRUE(inf$clustered)) inf$n_clusters else inf$n_total
  if (is.numeric(n_rows) && length(n_rows) == 1L && is.finite(n_rows) &&
      n_rows >= 1 && n_rows == floor(n_rows) && length(set$pre) > n_rows - 1) {
    return(c(base, list(status = "unavailable", df = length(set$pre), cells = set$ids,
      W = NULL, Wpval = NULL, covariance_rank_upper = n_rows - 1,
      n_independent_covariance_rows = n_rows,
      reason = paste0("the ", length(set$pre), " declared pretrend restrictions exceed the known centered covariance rank bound ",
        n_rows - 1, " from ", n_rows, if (isTRUE(inf$clustered)) " cluster sums" else " unit influence rows",
        "; no pseudoinverse or restriction removal is applied"),
      not_conditional_pretest = chained_pretest_not_conditional)))
  }
  V <- reporting$inference$vcov
  if (is.null(V)) V <- reporting$vcov
  pre <- set$pre
  if (is.null(relative_accuracy)) {
    A <- tryCatch(as.matrix(V), error = function(e) NULL)
    if (is.null(A) || !is.numeric(A) || !identical(dim(A), c(nrow(tm), nrow(tm))))
      return(c(base, list(status = "unavailable", df = length(pre), cells = set$ids, W = NULL, Wpval = NULL,
        reason = "the full reported covariance is missing or malformed; the declared null is unavailable without restriction removal",
        not_conditional_pretest = chained_pretest_not_conditional)))
    return(c(base, chained_pretest_wald(theta = as.double(reporting$estimates[pre]),
      V = unname(A[pre, pre, drop = FALSE]), ids = set$ids, max_matrix_entries = max_matrix_entries)))
  }
  if (is.null(cert) || is.null(cert$theta) || is.null(cert$vcov) || !isTRUE(cert$certified)) {
    return(c(base, list(status = "unavailable", df = length(pre), cells = set$ids, W = NULL, Wpval = NULL,
      reason = "no certified numerical enclosure of the reported cells is available, so the Wald arithmetic cannot be certified",
      not_conditional_pretest = chained_pretest_not_conditional)))
  }
  rec <- chained_pretest_wald(theta = as.double(reporting$estimates[pre]), V = unname(as.matrix(V)[pre, pre, drop = FALSE]),
    theta_lo = cert$theta$lo[, pre], theta_hi = cert$theta$hi[, pre], V_lo = cert$vcov$lo[pre, pre, drop = FALSE], V_hi = cert$vcov$hi[pre, pre, drop = FALSE],
    ids = set$ids, relative_accuracy = relative_accuracy, max_matrix_entries = max_matrix_entries)
  c(base, rec)
}
