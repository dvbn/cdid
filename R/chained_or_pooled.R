# Private D47 shared one-period outcome-change OR helpers. A pure layer over the
# accepted or_wls_* / or_pooled_wls_* kernels: the explicit outcome recipe,
# bounded support accumulation, ONE joint blockwise fit with the augmented-model
# alias convention, and complete original-unit IF rows. Owner scans, weights,
# forcing, budgets and seals stay in the route backends. Designs and outcome
# changes reach these helpers only on actual observed pair rows.

or_pooled_abort <- function(message, class = "input", route = "pooled",
                            also = character(), diagnostics = list()) {
  stop(structure(list(message = message, call = NULL, diagnostics = diagnostics),
    class = unique(c(paste0("cdid_or_", route, "_", class, "_error"),
      paste0("cdid_or_pooled_", class, "_error"), also,
      "cdid_ipw_fit_error", "cdid_error", "error", "condition"))))
}

# Explicit private column roles of the common frame Z_a. An intercept recipe has
# exactly one intercept column (verified == 1 on every scanned observed row; the
# kernel supplies its own 1-column); a removed intercept is never re-added.
or_pooled_outcome_recipe <- function(columns, intercept, slope_columns) {
  names_ok <- function(x) is.character(x) && !is.object(x) && is.null(dim(x)) &&
    !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  if (missing(columns) || missing(intercept) || missing(slope_columns) ||
      !names_ok(columns) || !length(columns) || !identical(intercept, isTRUE(intercept)) ||
      !names_ok(slope_columns) || !all(slope_columns %in% columns))
    or_pooled_abort(paste("Supply the explicit outcome recipe: unique Z_a column names,",
      "a TRUE/FALSE formula intercept flag and the slope columns."))
  slope_index <- which(columns %in% slope_columns)
  other <- setdiff(seq_along(columns), slope_index)
  if (intercept && length(other) != 1L)
    or_pooled_abort("An intercept recipe maps exactly one Z_a column to the formula intercept.")
  if (!intercept && length(other))
    or_pooled_abort("A no-intercept recipe maps every Z_a column to a slope; none is added.")
  if (!intercept && !length(slope_index))
    or_pooled_abort("Empty outcome formulas remain refused: no intercept and no slope column.")
  structure(list(columns = columns, intercept = intercept,
    intercept_index = if (intercept) other else 0L, slope_index = slope_index,
    slope_columns = columns[slope_index], p = length(slope_index),
    duration = "L_e = t_calendar_index - s_calendar_index on the stored modeled calendar",
    design = "d_e = L_e * z(X) with fixed-unit z", key = paste0("fixed-X-summed|intercept:",
      intercept, "|z:", paste(columns[slope_index], collapse = ","))),
    class = "cdid_or_pooled_recipe")
}

or_pooled_check_recipe <- function(recipe, columns, route = "pooled") {
  same <- inherits(recipe, "cdid_or_pooled_recipe") && is.list(recipe) &&
    identical(tryCatch(or_pooled_outcome_recipe(recipe$columns, recipe$intercept,
      recipe$slope_columns), error = function(e) NULL), recipe) &&
    identical(recipe$columns, columns)
  if (!isTRUE(same))
    or_pooled_abort("The outcome recipe must be the private recipe of exactly the frame Z_a columns.",
      "input", route)
  invisible(recipe)
}

# Modeled one-period transitions between the actual pair dates; raw numeric
# label differences (pair$time_gap) are never used.
or_pooled_transitions <- function(pair, route = "pooled") {
  ds <- pair$s_calendar_index; dt <- pair$t_calendar_index
  if (!is.numeric(ds) || !is.numeric(dt) || length(ds) != 1L || length(dt) != 1L ||
      !is.finite(ds) || !is.finite(dt) || ds != floor(ds) || dt != floor(dt) ||
      ds < 1 || dt <= ds)
    or_pooled_abort("Actual pair dates must lie on the stored modeled calendar with s before t.",
      "support", route)
  as.double(dt - ds)
}

# Raw-weight design summary of ONE (pair, cohort) support: [1, L z, dy] chunks
# with exact original-row constancy, or uncentered [L z, dy] without intercept.
or_pooled_support_start <- function(recipe, L) {
  p <- recipe$p
  D <- matrix(numeric(), 0L, p, dimnames = if (p) list(NULL, recipe$slope_columns))
  list(L = L, summary = if (recipe$intercept) or_pooled_wls_block_summary(D, numeric(), numeric())
    else or_wls_summary(D, numeric(), numeric()), weight_mass = 0, n = 0,
    sumZ = numeric(p), sumY = 0, colmax = numeric(p))
}

or_pooled_support_add <- function(acc, recipe, Z, dy, weights, route = "pooled") {
  n <- length(dy)
  if (!n) return(acc)
  if (!is.matrix(Z) || nrow(Z) != n || ncol(Z) != length(recipe$columns) || length(weights) != n)
    or_pooled_abort("Observed pair rows disagree with their design rows.", "integrity", route)
  if (recipe$intercept && any(Z[, recipe$intercept_index] != 1))
    or_pooled_abort("The declared intercept column must equal one on every observed pair row.",
      "input", route)
  z <- Z[, recipe$slope_index, drop = FALSE]
  D <- acc$L * z
  dimnames(D) <- if (recipe$p) list(NULL, recipe$slope_columns)
  acc$summary <- if (recipe$intercept)
    or_pooled_wls_merge_chunks(list(acc$summary, or_pooled_wls_block_summary(D, dy, weights))) else
    or_wls_merge(list(acc$summary, or_wls_summary(D, dy, weights)))
  acc$weight_mass <- acc$weight_mass + sum(weights); acc$n <- acc$n + n
  acc$sumZ <- acc$sumZ + unname(colSums(weights * z)); acc$sumY <- acc$sumY + sum(weights * dy)
  if (recipe$p) acc$colmax <- pmax(acc$colmax, unname(apply(abs(z), 2L, max)))
  acc
}

# Sum_i w x_i r_i and Sum_i w r_i^2 over ONE summary's rows, x = [1, d_keep] or
# d_keep, r = dy - alpha - d_keep' beta. Factor crossproducts only (F'F = M'M).
or_pooled_factor_moments <- function(recipe, summary, keep, alpha, beta) {
  F <- if (recipe$intercept) summary$raw$factor else summary$factor
  xi <- if (recipe$intercept) c(1L, 1L + keep) else keep
  if (!nrow(F)) return(list(g = numeric(length(xi)), e2 = 0))
  rho <- F[, recipe$p + 1L + recipe$intercept] +
    as.vector(F[, xi, drop = FALSE] %*% c(if (recipe$intercept) -alpha, -beta))
  list(g = as.vector(crossprod(F[, xi, drop = FALSE], rho)), e2 = sum(rho^2))
}

# Residual ratio of one link on its target support's raw-weight moments; a
# scalar support weight v cancels from both terms of this ONE ratio.
or_pooled_link_moments <- function(recipe, keep, mass, sumZ, sumY, L, alpha, beta) {
  fitted <- (if (recipe$intercept) alpha * mass else 0) + L * sum(sumZ[keep] * beta)
  mu <- (sumY - fitted) / mass
  list(mu = mu, c = L * sumZ[keep] / mass, residual = sumY - fitted - mu * mass)
}

# Strict reconstruction of proposed slope columns on actual original rows:
# z_P = c_b + z_R' A on every control row of block b and every treated row of its
# consuming links (c_b = 0 without intercept). Same scaled roundoff convention
# and threshold as or_public_basis.
or_pooled_alias_check <- function(recipe, retained, proposed, map, constants, roles,
                                  colmax, visit_rows) {
  ctl <- trt <- numeric(length(proposed)); rows <- 0
  thr <- 64 * .Machine$double.eps * max(1L, length(recipe$columns))
  for (sid in seq_along(roles)) {
    r <- roles[[sid]]
    if (!length(r$control) && !length(r$treated)) next
    visit_rows(sid, function(Z) {
      if (!nrow(Z)) return(invisible(NULL))
      z <- Z[, recipe$slope_index, drop = FALSE]; n <- nrow(z)
      u <- z[, proposed, drop = FALSE] - z[, retained, drop = FALSE] %*% map
      s0 <- abs(z[, proposed, drop = FALSE]) + abs(z[, retained, drop = FALSE]) %*% abs(map)
      for (b in union(r$control, r$treated)) {
        cc <- constants[b, ]
        ratio <- apply(abs(u - rep(cc, each = n)) / pmax(s0 + rep(abs(cc), each = n),
          rep(colmax, each = n), .Machine$double.xmin), 2L, max)
        ratio[!is.finite(ratio)] <- Inf
        if (b %in% r$control) ctl <<- pmax(ctl, ratio)
        if (b %in% r$treated) trt <<- pmax(trt, ratio)
      }
      rows <<- rows + n
    })
  }
  list(control = ctl, treated = trt, threshold = thr, rows_checked = rows)
}

# Reduced accepted summary on the kept slope columns: selecting factor columns
# preserves F'F for that sub-design (and blockwise exact zero support).
or_pooled_select <- function(recipe, pooled, keep) {
  F <- if (recipe$intercept) pooled$residual$factor else pooled$factor
  fac <- or_wls_compress(F[, c(keep, recipe$p + 1L), drop = FALSE])
  cols <- if (length(keep)) recipe$slope_columns[keep]
  if (!recipe$intercept) {
    ans <- structure(list(factor = fac, n = pooled$n, columns = cols), class = "cdid_or_wls_summary")
    or_wls_check_summary(ans)
    return(ans)
  }
  pooled$residual <- if (length(keep)) structure(list(factor = fac, n = pooled$n,
    columns = cols), class = "cdid_or_wls_summary") else list(factor = fac, n = pooled$n, columns = NULL)
  pooled$p <- length(keep); pooled$columns <- cols
  pooled$blocks <- lapply(pooled$blocks, function(b) {
    b$r1d <- b$r1d[keep]; b$mean_d <- b$mean_d[keep]; b
  })
  pooled$accounting$p <- length(keep); pooled$accounting$residual_dim <- dim(fac)
  or_pooled_wls_check_summary(pooled)
  pooled
}

# ONE joint fit over unique blocks. supports[[s]] = list(summary = raw-weight
# support summary, scale = fixed weight multiplier v_s, colmax); blocks[[b]] =
# list(name, support_ids, L); links[[l]] = list(regression, support_id).
# visit_rows(sid, visit) replays bounded chunks of a support's observed Z_a rows.
or_pooled_fit_blocks <- function(recipe, supports, blocks, links, rank_tol, N,
                                 visit_rows, route = "pooled") {
  p <- recipe$p; icpt <- recipe$intercept
  scaled <- function(s) {
    x <- supports[[s]]$summary
    if (icpt) x$raw$factor <- x$raw$factor * sqrt(supports[[s]]$scale) else
      x$factor <- x$factor * sqrt(supports[[s]]$scale)
    x
  }
  parts <- lapply(blocks, function(b) {
    merged <- if (icpt) or_pooled_wls_merge_chunks(lapply(b$support_ids, scaled)) else
      or_wls_merge(lapply(b$support_ids, scaled))
    if (icpt) or_pooled_wls_partial_block(merged) else merged
  })
  names(parts) <- vapply(blocks, `[[`, character(1), "name")
  pooled <- if (icpt) or_pooled_wls_merge(parts) else or_wls_merge(unname(parts))
  keep <- seq_len(p); retained <- keep; proposed <- integer()
  check <- NULL; aliased <- candidates <- integer()
  map <- matrix(0, p, 0L); constants <- matrix(0, length(blocks), 0L)
  if (p) {
    F <- if (icpt) pooled$residual$factor else pooled$factor
    X <- F[, seq_len(p), drop = FALSE]
    proposal <- qr(X, tol = rank_tol, LAPACK = FALSE)
    retained <- sort(proposal$pivot[seq_len(proposal$rank)])
    proposed <- setdiff(seq_len(p), retained)
  }
  if (length(proposed)) {
    map <- matrix(0, length(retained), length(proposed))
    if (length(retained)) {
      kept_qr <- qr(X[, retained, drop = FALSE], tol = rank_tol, LAPACK = FALSE)
      map <- qr.coef(kept_qr, X[, proposed, drop = FALSE])
      if (!anyNA(map)) map <- map + qr.coef(kept_qr,
        X[, proposed, drop = FALSE] - X[, retained, drop = FALSE] %*% map)
    }
    map <- matrix(map, length(retained), length(proposed))
    constants <- matrix(0, length(blocks), length(proposed))
    if (icpt) for (b in seq_along(blocks)) {
      zbar <- parts[[b]]$blocks[[1L]]$mean_d / blocks[[b]]$L
      constants[b, ] <- zbar[proposed] - as.vector(zbar[retained] %*% map)
    }
    roles <- lapply(seq_along(supports), function(s) list(
      control = which(vapply(blocks, function(b) s %in% b$support_ids, logical(1))),
      treated = sort(unique(unlist(lapply(links, function(l)
        if (l$support_id == s) l$regression else integer()))))))
    used <- vapply(roles, function(r) length(r$control) + length(r$treated) > 0L, logical(1))
    colmax <- do.call(pmax, lapply(supports[used], function(u) u$colmax))[proposed]
    check <- or_pooled_alias_check(recipe, retained, proposed, map, constants, roles,
      colmax, visit_rows)
    exact <- check$control <= check$threshold
    if (any(exact & check$treated > check$threshold))
      or_pooled_abort(paste0("Pooled control design has an exact null direction (",
        paste(recipe$slope_columns[proposed[exact & check$treated > check$threshold]],
          collapse = ", "), ") that does not reconstruct the required treated predictions; ",
        "no ridge or singular-contrast solver is used; model unchanged."), "support", route,
        "cdid_or_wls_rank_error", list(control_scaled = check$control,
          treated_scaled = check$treated, threshold = check$threshold))
    alias <- exact & check$treated <= check$threshold
    aliased <- proposed[alias]; candidates <- proposed[!alias]
    keep <- sort(c(retained, candidates))
    if (!length(keep) && !icpt)
      or_pooled_abort(paste("Outcome design has no retained columns on actual control and treated",
        "rows; no reduced control regression can be computed."), "support", route,
        "cdid_or_wls_rank_error")
    if (length(keep) < p) pooled <- or_pooled_select(recipe, pooled, keep)
    map <- map[, alias, drop = FALSE]; constants <- constants[, alias, drop = FALSE]
  }
  fit <- if (icpt) or_pooled_wls_fit(pooled, rank_tol) else or_wls_fit(pooled, rank_tol)
  beta <- unname(fit$coefficients)
  dimnames(map) <- list(recipe$slope_columns[retained], recipe$slope_columns[aliased])
  dimnames(constants) <- list(names(parts), recipe$slope_columns[aliased])
  regressions <- lapply(seq_along(blocks), function(b) {
    rec <- if (icpt) fit$blocks[[b]]
    list(name = blocks[[b]]$name, L = blocks[[b]]$L, n = parts[[b]]$n,
      alpha = if (icpt) unname(fit$intercepts[b]) else 0,
      weight_sum = if (icpt) rec$weight_sum else NA_real_,
      mass = if (icpt) rec$weight_sum / N else NA_real_,
      mean_d = if (icpt) unname(rec$mean_d) else numeric(length(keep)),
      mean_y = if (icpt) rec$mean_y else NA_real_,
      residual_rows = if (icpt) parts[[b]]$blocks[[1L]]$residual_rows else nrow(parts[[b]]$factor))
  })
  basis <- list(slope_columns = recipe$slope_columns, kept_columns = recipe$slope_columns[keep],
    kept_index = keep, proposed_columns = recipe$slope_columns[proposed],
    aliased_columns = recipe$slope_columns[aliased],
    nontrivial_candidates_retained = recipe$slope_columns[candidates],
    alias_map = map, block_constants = constants,
    reconstruction_control_scaled = if (is.null(check)) numeric() else check$control,
    reconstruction_treated_scaled = if (is.null(check)) numeric() else check$treated,
    reconstruction_roundoff_threshold = if (is.null(check)) NA_real_ else check$threshold,
    rows_checked = if (is.null(check)) 0 else check$rows_checked,
    proposal = "R stats::lm.fit / dqrdc2 pivoted QR at rank_tol on the pooled (partialled) factor",
    rank_check = "unchanged relative-diagonal rule of the accepted kernel on the reduced factor",
    frozen_during_differentiation = TRUE, numerical_certificate = FALSE)
  list(fit = fit, summary = pooled, keep = keep, beta = beta, regressions = regressions,
    basis = basis, rss = fit$rss, objective = fit$rss / N, verification_scan = !is.null(check))
}

# Common-design collapsed OR fit on the eligible complete-support route:
# - For p = 0: direct normalized means alpha_b = H(q * rawDelta_b) / H(q) with
#   no QR decomposition.
# - For p >= 1: ONE Householder QR summary of sqrt(q) [1, z, ystar] with
#   D = sum_b v_b L_b^2 and ystar = sum_b (v_b L_b / D) rawDelta_b, factor
#   scaling sqrt(D), existing alias proposal/verification, triangular solve for
#   beta and alpha_b = mean_q(Delta_b) - L_b mean_q(z)' beta BEFORE Pass 2, and
#   exact objective reconstruction Q(alpha, beta) = K + D * Q_star(beta) from
#   Pass 2's sum-of-squares K.
or_pooled_fit_common <- function(recipe, supports, blocks, links, rank_tol, N,
                                 visit_rows, common_part = NULL, pass2_fn,
                                 route = "baseline") {
  p <- recipe$p; G <- length(blocks)
  block_names <- vapply(blocks, `[[`, character(1), "name")
  L_vec <- vapply(blocks, `[[`, numeric(1), "L")
  s_vec <- vapply(blocks, function(b) b$support_ids[1L], integer(1))
  v_vec <- vapply(s_vec, function(s) supports[[s]]$scale, numeric(1))
  D <- sum(v_vec * L_vec^2)
  keep <- seq_len(p); retained <- keep; proposed <- integer()
  check <- NULL; aliased <- candidates <- integer()
  map <- matrix(0, p, 0L); constants <- matrix(0, G, 0L)
  F_scaled <- NULL
  if (p > 0L) {
    F_scaled <- common_part$residual$factor * sqrt(D)
    X <- F_scaled[, seq_len(p), drop = FALSE]
    proposal <- qr(X, tol = rank_tol, LAPACK = FALSE)
    retained <- sort(proposal$pivot[seq_len(proposal$rank)])
    proposed <- setdiff(seq_len(p), retained)
  }
  if (length(proposed)) {
    map <- matrix(0, length(retained), length(proposed))
    if (length(retained)) {
      kept_qr <- qr(X[, retained, drop = FALSE], tol = rank_tol, LAPACK = FALSE)
      map <- qr.coef(kept_qr, X[, proposed, drop = FALSE])
      if (!anyNA(map)) map <- map + qr.coef(kept_qr,
        X[, proposed, drop = FALSE] - X[, retained, drop = FALSE] %*% map)
    }
    map <- matrix(map, length(retained), length(proposed))
    constants <- matrix(0, G, length(proposed))
    zbar <- common_part$blocks[[1L]]$mean_d
    c_vec <- zbar[proposed] - as.vector(zbar[retained] %*% map)
    for (b in seq_len(G)) constants[b, ] <- c_vec
    roles <- lapply(seq_along(supports), function(s) list(
      control = which(vapply(blocks, function(b) s %in% b$support_ids, logical(1))),
      treated = sort(unique(unlist(lapply(links, function(l)
        if (l$support_id == s) l$regression else integer()))))))
    used <- vapply(roles, function(r) length(r$control) + length(r$treated) > 0L, logical(1))
    colmax <- do.call(pmax, lapply(supports[used], function(u) u$colmax))[proposed]
    check <- or_pooled_alias_check(recipe, retained, proposed, map, constants, roles,
      colmax, visit_rows)
    exact <- check$control <= check$threshold
    if (any(exact & check$treated > check$threshold))
      or_pooled_abort(paste0("Pooled control design has an exact null direction (",
        paste(recipe$slope_columns[proposed[exact & check$treated > check$threshold]],
          collapse = ", "), ") that does not reconstruct the required treated predictions; ",
        "no ridge or singular-contrast solver is used; model unchanged."), "support", route,
        "cdid_or_wls_rank_error", list(control_scaled = check$control,
          treated_scaled = check$treated, threshold = check$threshold))
    alias <- exact & check$treated <= check$threshold
    aliased <- proposed[alias]; candidates <- proposed[!alias]
    keep <- sort(c(retained, candidates))
    map <- map[, alias, drop = FALSE]; constants <- constants[, alias, drop = FALSE]
  }
  p_keep <- length(keep)
  cols_keep <- if (p_keep) recipe$slope_columns[keep] else NULL
  blocks_rec <- lapply(seq_len(G), function(b) {
    u <- supports[[s_vec[b]]]; v_b <- v_vec[b]; L_b <- L_vec[b]
    m_b <- u$weight_mass * v_b
    sumY_b <- u$sumY * v_b
    mean_y_b <- sumY_b / m_b
    r11_b <- sqrt(m_b)
    mean_d_b <- if (p_keep) L_b * ((u$sumZ[keep] * v_b) / m_b) else numeric(0)
    names(mean_d_b) <- cols_keep
    r1d_b <- if (p_keep) setNames(r11_b * mean_d_b, cols_keep) else numeric(0)
    r1y_b <- r11_b * mean_y_b
    aug_rows <- if (p > 0L) common_part$blocks[[1L]]$augmented_rows else 2L
    res_rows <- if (p > 0L) common_part$blocks[[1L]]$residual_rows else 1L
    list(r11 = r11_b, r1d = r1d_b, r1y = r1y_b, n = as.integer(u$n), weight_sum = m_b,
      mean_d = mean_d_b, mean_y = mean_y_b, augmented_rows = aug_rows,
      residual_rows = res_rows)
  })
  names(blocks_rec) <- block_names
  n_total <- sum(vapply(blocks_rec, `[[`, integer(1), "n"))
  F_keep <- NULL; slope_fit_pre <- NULL
  if (p_keep > 0L) {
    F_keep <- if (p_keep < p) or_wls_compress(F_scaled[, c(keep, p + 1L), drop = FALSE]) else F_scaled
    slope_fit_pre <- or_wls_fit(structure(list(factor = F_keep, n = n_total, columns = cols_keep),
      class = "cdid_or_wls_summary"), rank_tol)
    beta <- unname(slope_fit_pre$coefficients)
  } else {
    beta <- numeric(0)
  }
  alpha <- vapply(blocks_rec, function(rec) rec$mean_y - sum(rec$mean_d * beta), numeric(1))
  names(alpha) <- block_names
  mean_y_vec <- vapply(blocks_rec, `[[`, numeric(1), "mean_y")
  y_star_bar <- if (p_keep > 0L) sum((v_vec * L_vec / D) * mean_y_vec) else 0
  e_bar <- if (p_keep > 0L) mean_y_vec - L_vec * y_star_bar else numeric(G)

  pass2 <- pass2_fn(keep, alpha, beta, e_bar)

  if (p_keep > 0L) {
    K_total <- sum(pass2$moments[3L, ])
    q0 <- qr(F_keep[, seq_len(p_keep), drop = FALSE], LAPACK = TRUE)
    R_unpiv <- qr.R(q0)[seq_len(p_keep), order(q0$pivot), drop = FALSE]
    qty0 <- qr.qty(q0, F_keep[, p_keep + 1L])
    r_tail <- sqrt(sum(qty0[-seq_len(p_keep)]^2) + K_total)
    fac_final <- unname(rbind(cbind(R_unpiv, qty0[seq_len(p_keep)]), c(rep.int(0, p_keep), r_tail)))
    if (any(!is.finite(fac_final)))
      cdid_abort("Nonfinite compressed WLS factor.", "cdid_or_wls_numeric_error")
    res_summary <- structure(list(factor = fac_final, n = n_total, columns = cols_keep),
      class = "cdid_or_wls_summary")
    rss_final <- slope_fit_pre$rss + K_total
    slope_fit <- slope_fit_pre; slope_fit$rss <- rss_final
  } else {
    rss_final <- sum(v_vec * pass2$moments[2L, ])
    fac_final <- matrix(sqrt(rss_final), 1L, 1L)
    if (any(!is.finite(fac_final)))
      cdid_abort("Nonfinite compressed WLS factor.", "cdid_or_wls_numeric_error")
    res_summary <- list(factor = fac_final, n = n_total, columns = NULL)
    slope_fit <- NULL
  }
  pooled_acct <- list(
    n = n_total,
    p = p_keep,
    n_blocks = G,
    residual_dim = dim(fac_final),
    max_augmented_rows = if (p > 0L) common_part$accounting$max_augmented_rows else 2L,
    max_residual_rows = if (p > 0L) common_part$accounting$max_residual_rows else 1L,
    max_chunk_merge_rows = if (p > 0L) common_part$accounting$max_chunk_merge_rows else 1L,
    max_residual_merge_rows = nrow(fac_final),
    n_compressions = if (p > 0L) common_part$accounting$n_compressions + (if (p_keep < p) 1L else 0L) else 0L
  )
  pooled <- structure(list(residual = res_summary, blocks = blocks_rec, n = n_total,
    p = p_keep, columns = cols_keep, accounting = pooled_acct), class = "cdid_or_pooled_wls_summary")
  or_pooled_wls_check_summary(pooled)
  if (any(!is.finite(c(alpha, beta, rss_final))))
    cdid_abort("Nonfinite final pooled WLS solution.", "cdid_or_wls_numeric_error")
  fit <- structure(list(
    coefficients = if (p_keep > 0L) slope_fit$coefficients else setNames(numeric(0), cols_keep),
    intercepts = alpha,
    rss = rss_final,
    R = if (p_keep > 0L) slope_fit$R else matrix(numeric(0), 0L, 0L),
    pivot = if (p_keep > 0L) slope_fit$pivot else integer(0),
    n = n_total,
    p = p_keep,
    n_blocks = G,
    rank = if (p_keep > 0L) slope_fit$rank else 0L,
    rank_tol = rank_tol,
    cutoff = if (p_keep > 0L) slope_fit$cutoff else 0,
    blocks = blocks_rec,
    slope_fit = slope_fit,
    accounting = pooled_acct
  ), class = "cdid_or_pooled_wls_fit")

  dimnames(map) <- list(recipe$slope_columns[retained], recipe$slope_columns[aliased])
  dimnames(constants) <- list(block_names, recipe$slope_columns[aliased])
  regressions <- lapply(seq_len(G), function(b) {
    rec <- fit$blocks[[b]]
    list(name = blocks[[b]]$name, L = blocks[[b]]$L, n = rec$n,
      alpha = unname(fit$intercepts[b]),
      weight_sum = rec$weight_sum,
      mass = rec$weight_sum / N,
      mean_d = unname(rec$mean_d),
      mean_y = rec$mean_y,
      residual_rows = rec$residual_rows)
  })
  basis <- list(slope_columns = recipe$slope_columns, kept_columns = recipe$slope_columns[keep],
    kept_index = keep, proposed_columns = recipe$slope_columns[proposed],
    aliased_columns = recipe$slope_columns[aliased],
    nontrivial_candidates_retained = recipe$slope_columns[candidates],
    alias_map = map, block_constants = constants,
    reconstruction_control_scaled = if (is.null(check)) numeric() else check$control,
    reconstruction_treated_scaled = if (is.null(check)) numeric() else check$treated,
    reconstruction_roundoff_threshold = if (is.null(check)) NA_real_ else check$threshold,
    rows_checked = if (is.null(check)) 0 else check$rows_checked,
    proposal = "R stats::lm.fit / dqrdc2 pivoted QR at rank_tol on the pooled (partialled) factor",
    rank_check = "unchanged relative-diagonal rule of the accepted kernel on the reduced factor",
    frozen_during_differentiation = TRUE, numerical_certificate = FALSE)
  list(fit = fit, summary = pooled, keep = keep, beta = beta, regressions = regressions,
    basis = basis, rss = fit$rss, objective = fit$rss / N, verification_scan = !is.null(check),
    pass2 = pass2)
}

# Rows needed by one link batch: its local blocks, and EVERY block and control
# support whenever a shared slope exists (xi_beta is global) or details are on.
or_pooled_if_sets <- function(model, links, details = FALSE) {
  gids <- sort(unique(vapply(model$links[links], `[[`, integer(1), "regression")))
  beta_gids <- if (model$shared$p > 0L || details) seq_along(model$regression_specs) else integer()
  sids <- sort(unique(c(vapply(model$links[links], `[[`, integer(1), "support_id"),
    unlist(lapply(model$regression_specs[union(gids, beta_gids)], `[[`, "support_ids")))))
  list(gids = gids, beta_gids = beta_gids, sids = sids)
}

# Upstream loadings H(T t'), H(T r t') and H(T d t') of one link's treated ratio
# terms on its support's upstream coordinates. A baseline rate is one coordinate
# whose t = -1/r is constant on the support rows (log_weight_derivative).
or_pooled_link_loadings <- function(l, p) {
  if (!is.null(l$weight_loading)) return(list(weight = l$weight_loading,
    numerator = l$numerator_loading, design = l$design_loading))
  t <- l$log_weight_derivative
  list(weight = l$d * t, numerator = l$numerator * t, design = matrix(l$c * l$d * t, p, 1L))
}

# Complete original-unit rows for one bounded row chunk. hits[[j]] holds the
# observed rows (chunk positions), dy and kept slope z of support pre$sids[j],
# plus their fitted factors v when these vary by row (else the support's v),
# and their raw recipe weights w when the support reads a reference recipe
# (D60; else the chunk's common weights).
# Without pre$uids, up[, j] is that support's single upstream IF column; with
# pre$uids, up holds those coordinates and support s reads its
# model$supports[[s]]$channels, a coordinate shared by supports being ONE column.
# Unit scores of every block are summed before any outer product, retaining
# reused units. D90: beta_only = TRUE (with details = FALSE only) returns
# list(beta = xi_beta) right after the unchanged global shared-beta score and
# solve, for a caller that consumes no block-alpha, direct or standalone-OR
# link phi; the default returns every output unchanged.
or_pooled_if_core <- function(model, N, w, hits, up, pre, links, details = FALSE, beta_only = FALSE) {
  if (!identical(beta_only, FALSE) && (!identical(beta_only, TRUE) || !identical(details, FALSE)))
    or_pooled_abort("A private beta-only OR core request requires details = FALSE.")
  B <- length(w); p <- model$shared$p; icpt <- model$outcome_recipe$intercept
  beta <- model$shared$beta
  col <- function(s) match(s, pre$sids)
  uid <- function(s) if (is.null(pre$uids)) s else
    sort(unique(unlist(lapply(model$supports[s], `[[`, "channels"))))
  ucol <- function(s) if (is.null(pre$uids)) col(s) else match(uid(s), pre$uids)
  cache <- list()
  resid <- function(b, s) {
    key <- paste(b, s)
    if (!is.null(cache[[key]])) return(cache[[key]])
    h <- hits[[col(s)]]; blk <- model$regressions[[b]]
    d <- blk$L * h$z
    e <- list(ix = h$obs, d = d, q = (if (is.null(h$w)) w[h$obs] else h$w) *
        (if (is.null(h$v)) model$supports[[s]]$v else h$v),
      r = h$dy - blk$alpha - as.vector(d %*% beta))
    cache[[key]] <<- e
    e
  }
  h_beta <- xi_beta <- matrix(0, B, p)
  block_scores <- if (details) lapply(pre$beta_gids, function(b) matrix(0, B, p)) else list()
  if (p) {
    for (k in seq_along(pre$beta_gids)) {
      b <- pre$beta_gids[k]; hb <- matrix(0, B, p)
      for (s in model$regression_specs[[b]]$support_ids) {
        e <- resid(b, s)
        if (length(e$ix)) hb[e$ix, ] <- hb[e$ix, ] + (e$q * e$r) *
          (if (icpt) sweep(e$d, 2L, model$regressions[[b]]$mean_d) else e$d)
      }
      h_beta <- h_beta + hb
      if (details) block_scores[[k]] <- hb
    }
    ctrl <- model$shared$control_support_ids
    rhs <- h_beta + up[, ucol(ctrl), drop = FALSE] %*%
      t(model$shared$forcing_beta[, uid(ctrl), drop = FALSE])
    xi_beta <- t(if (icpt) or_pooled_wls_gram_solve(model$shared$fit, t(rhs), normalizer = N) else
      or_wls_gram_solve(model$shared$fit, t(rhs), normalizer = N))
  }
  if (beta_only) return(list(beta = xi_beta))
  h_alpha <- alpha <- matrix(0, B, length(pre$gids))
  if (icpt) for (j in seq_along(pre$gids)) {
    b <- pre$gids[j]; blk <- model$regressions[[b]]; sids <- model$regression_specs[[b]]$support_ids
    for (s in sids) {
      e <- resid(b, s)
      if (length(e$ix)) h_alpha[e$ix, j] <- h_alpha[e$ix, j] + e$q * e$r
    }
    alpha[, j] <- (h_alpha[, j] + as.vector(up[, ucol(sids), drop = FALSE] %*% blk$forcing_alpha)) /
      blk$mass - as.vector(xi_beta %*% blk$mean_d)
  }
  phi <- direct <- matrix(0, B, length(links))
  if (details) {
    den <- num <- matrix(0, B, length(links)); design <- vector("list", length(links))
  }
  for (j in seq_along(links)) {
    l <- model$links[[links[j]]]; b <- l$regression; s <- l$support_id
    e <- resid(b, s); ix <- e$ix
    if (length(ix)) direct[ix, j] <- e$q * (e$r - l$mu) / l$d
    us <- up[, ucol(s), drop = FALSE]
    phi[, j] <- direct[, j] - (if (icpt) alpha[, match(b, pre$gids)] else 0) -
      as.vector(xi_beta %*% l$c) + as.vector(us %*% l$rate_loading)
    if (details) {
      # Treated denominator, numerator and design-sum rows of this one ratio.
      tw <- numeric(B); tr <- numeric(B); tw[ix] <- e$q; tr[ix] <- e$q * e$r
      hd <- l$c * l$d; num_l <- l$numerator; lo <- or_pooled_link_loadings(l, p)
      den[, j] <- tw - l$d + as.vector(us %*% lo$weight)
      num[, j] <- tr - num_l + as.vector(us %*% lo$numerator) -
        (if (icpt) l$d * alpha[, match(b, pre$gids)] else 0) - as.vector(xi_beta %*% hd)
      dm <- matrix(0, B, p)
      if (length(ix) && p) dm[ix, ] <- e$q * e$d
      design[[j]] <- dm - matrix(hd, B, p, byrow = TRUE) + us %*% t(lo$design)
    }
  }
  ans <- list(phi = phi, direct = direct, beta = xi_beta, alpha = alpha)
  if (details) {
    energy <- numeric(B)
    for (b in pre$beta_gids) for (s in model$regression_specs[[b]]$support_ids) {
      e <- resid(b, s)
      if (length(e$ix)) energy[e$ix] <- energy[e$ix] + e$q * e$r^2
    }
    ctrl <- model$shared$control_support_ids
    ans <- c(ans, list(beta_score = h_beta, beta_block_scores = block_scores,
      alpha_score = h_alpha, denominator = den, numerator = num, design = design,
      objective = energy - model$shared$objective +
        as.vector(up[, ucol(ctrl), drop = FALSE] %*% model$shared$forcing_energy[uid(ctrl)])))
  }
  ans
}
