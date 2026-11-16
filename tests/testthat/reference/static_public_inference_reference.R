# ==============================================================================
# Narrow Independent Reference for Public Static Inference & Switch-Share Aggregate
# ==============================================================================
#
# Pure base-R/stats implementation of:
#   1. Centered cluster-sum analytical covariance and common-sign Rademacher
#      multiplier bootstrap.
#   2. Switcher-frequency aggregate point estimate and random switch-share
#      influence function.
#
# Uses ZERO production `cdid` functions.
# ==============================================================================

ref_static_public_cluster_inference <- function(phi_mat,
                                                cluster_vec = NULL,
                                                signs_mat = NULL,
                                                alp = 0.05) {
  Phi <- as.matrix(phi_mat)
  storage.mode(Phi) <- "double"
  n <- nrow(Phi)
  J <- ncol(Phi)
  cnames <- colnames(Phi)
  if (is.null(cnames)) cnames <- paste0("t", seq_len(J))

  if (is.null(cluster_vec)) {
    G <- n
    S_raw <- Phi
  } else {
    cl_chr <- as.character(cluster_vec)
    u_cl <- unique(cl_chr)
    G <- length(u_cl)
    S_raw <- matrix(0.0, nrow = G, ncol = J, dimnames = list(u_cl, cnames))
    for (g in seq_len(G)) {
      idx <- which(cl_chr == u_cl[g])
      S_raw[g, ] <- colSums(Phi[idx, , drop = FALSE])
    }
  }

  # Centered cluster sums: \tilde{S}_g = \hat{S}_g - G^{-1} \sum_{h=1}^G \hat{S}_h
  S_tilde <- sweep(S_raw, 2L, colMeans(S_raw), FUN = "-")

  # Analytical covariance: \hat{V} = n^{-2} \sum_{g=1}^G \tilde{S}_g \tilde{S}_g'
  V_an <- crossprod(S_tilde) / (as.double(n)^2)
  dimnames(V_an) <- list(cnames, cnames)
  an_var <- diag(V_an)
  an_se <- ifelse(is.finite(an_var) & an_var > 0.0, sqrt(an_var), NA_real_)
  names(an_se) <- cnames

  if (is.null(signs_mat)) {
    return(list(
      n = n,
      G = G,
      S_raw = S_raw,
      S_tilde = S_tilde,
      V_analytical = V_an,
      analytical_se = an_se,
      raw_draws = NULL,
      iqr_scale = stats::setNames(rep(NA_real_, J), cnames),
      bootstrap_se = stats::setNames(rep(NA_real_, J), cnames),
      crit_val = NA_real_,
      accepted = all(is.finite(an_se) & an_se > 0.0)
    ))
  }

  B <- nrow(signs_mat)
  stopifnot(ncol(signs_mat) == G)

  # Multiplier draws: R^*_b = G^{-1/2} \sum_{g=1}^G \eta_{b,g} \tilde{S}_g
  R_draws <- (signs_mat %*% S_tilde) / sqrt(as.double(G))
  colnames(R_draws) <- cnames

  # Type-1 empirical quartiles: k = ceiling(p * B)
  k25 <- max(1L, as.integer(ceiling(0.25 * B)))
  k75 <- max(1L, as.integer(ceiling(0.75 * B)))
  k_crit <- max(1L, as.integer(ceiling((1.0 - alp) * B)))
  iqr_norm <- stats::qnorm(0.75) - stats::qnorm(0.25)
  se_factor <- sqrt(as.double(G)) / as.double(n)

  iqr_scale <- stats::setNames(rep(NA_real_, J), cnames)
  boot_se <- stats::setNames(rep(NA_real_, J), cnames)
  target_ok <- logical(J)

  for (j in seq_len(J)) {
    ord_j <- sort(R_draws[, j], method = "quick")
    iqr_j <- ord_j[k75] - ord_j[k25]
    sc_j <- iqr_j / iqr_norm
    se_j <- sc_j * se_factor
    if (is.finite(an_se[j]) && an_se[j] > 0.0 && is.finite(sc_j) && sc_j > 0.0) {
      iqr_scale[j] <- sc_j
      boot_se[j] <- se_j
      target_ok[j] <- TRUE
    } else {
      target_ok[j] <- FALSE
    }
  }

  all_ok <- all(target_ok) && all(is.finite(V_an))
  crit_val <- NA_real_
  max_t <- rep(NA_real_, B)
  if (all_ok) {
    t_mat <- sweep(abs(R_draws), 2L, iqr_scale, FUN = "/")
    max_t <- apply(t_mat, 1L, max)
    ord_max <- sort(max_t, method = "quick")
    crit_val <- unname(ord_max[k_crit])
  }

  list(
    n = n,
    G = G,
    S_raw = S_raw,
    S_tilde = S_tilde,
    V_analytical = V_an,
    analytical_se = an_se,
    raw_draws = R_draws,
    iqr_scale = iqr_scale,
    bootstrap_se = boot_se,
    max_t = max_t,
    crit_val = crit_val,
    accepted = all_ok
  )
}

ref_static_public_aggregate_if <- function(theta_hat, phi_theta, G_mat, signs) {
  Phi <- as.matrix(phi_theta)
  G_m <- as.matrix(G_mat)
  n <- nrow(Phi)
  J <- ncol(Phi)
  cids <- colnames(Phi)
  th <- as.double(theta_hat)
  s_vec <- as.double(signs)

  pi_hat <- colMeans(G_m)
  Pi_hat <- sum(pi_hat)
  omega_hat <- pi_hat / Pi_hat
  signed_theta <- s_vec * th
  delta_hat <- sum(omega_hat * signed_theta)

  cell_term <- drop(Phi %*% (omega_hat * s_vec))
  G_cent <- sweep(G_m, 2L, pi_hat, FUN = "-")
  share_term <- drop(G_cent %*% (signed_theta - delta_hat)) / Pi_hat
  phi_agg <- cell_term + share_term
  names(phi_agg) <- rownames(Phi)

  list(
    delta_hat = delta_hat,
    omega_hat = stats::setNames(omega_hat, cids),
    pi_hat = stats::setNames(pi_hat, cids),
    Pi_hat = Pi_hat,
    phi = phi_agg,
    cell_term = cell_term,
    share_term = share_term
  )
}
