# Frozen ordinary CP04 evaluator, exact body from PR11 4a8db4807a811fb85b87f09f22380682955e8e4c.
# Function binding renamed only; no numerical expression changed.
cp04_v56_oracle <- function(backend, rows, links, scale = "frame", details = FALSE) {
  validate_cp04_backend(backend)
  if (isTRUE(backend$uses_selection_models))
    return(evaluate_cp04_selection_if_chunk(backend, rows, links, scale = scale, details = details))
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state; ctrl <- backend$control
  if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || !length(rows) || anyNA(rows) || any(rows != floor(rows)) ||
      any(rows < 1 | rows > st$N) || anyDuplicated(rows) || length(rows) > ctrl$max_chunk_rows ||
      !is.numeric(links) || is.object(links) || !is.null(dim(links)) || !length(links) || anyNA(links) || any(links != floor(links)) ||
      any(links < 1 | links > length(backend$links)) || anyDuplicated(links) ||
      !identical(details, isTRUE(details)) || !is.character(scale) || length(scale) != 1L || is.na(scale) ||
      !scale %in% c("frame", "core")) cp04_abort("Invalid explicit IF row/link batch, scale or details flag.")
  pre <- cp04_if_budget(backend, length(rows), links)
  tids <- pre$tids; sids <- pre$sids
  B <- length(rows); d <- ncol(st$Z)
  # The original32-fold row scratch leaves24-fold for all existing
  # rows/solves/results; only the other8-fold can fund shared position reads.
  pair_op <- NULL
  if (!details && is.null(backend$reference_ref) && is.null(backend$layout) && is.null(backend$pair_outcomes)) {
    frame_p <- get("frame_owner", reader$frame_ref, inherits = FALSE)$n_covariates
    widths <- length(sids) * (d + 2) + sum(vapply(backend$treatment_specs[tids],
      function(s) 2 * frame_p + 10 + length(s$cohorts) * (d + 1), numeric(1))) +
      3 * length(st$cohorts) + length(links) + 16
    pair_op <- cp04_pair_read_operation(reader, rows, 8 * B * (widths + frame_p + d + 4))
  }
  on.exit(if (!is.null(pair_op)) pair_op$abort(), add = TRUE)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE); p <- frame$n_covariates
  # Scores evaluated once for each actual model needed by this requested batch.
  ci <- st$cohort_index[rows]; w <- st$w[rows]; A <- as.double(st$a_index[rows] > 0L)
  a_if <- matrix(0, B, length(st$cohorts))
  for (k in unique(ci)) {
    rr <- which(ci == k & st$U[rows])
    a_if[rr, k] <- w[rr] * (A[rr] - st$a[k]) / st$mass_U[k]
  }
  alpha <- probabilities <- responses <- vector("list", length(backend$selections))
  for (sj in sids) {
    s <- backend$selection_specs[[sj]]; f <- backend$selections[[sj]]
    rr <- which(ci == s$cohort & st$U[rows])
    xi <- matrix(0, B, d); pr <- S <- numeric(B)
    if (length(rr)) {
      pr[rr] <- cp04_predict_selection(reader, rows[rr], s, f)
      S[rr] <- (if (is.null(pair_op)) cp04_pair_rows(reader, rows[rr], s$pair) else
        pair_op$pair(rows[rr], s$pair))$S
      score <- w[rr] * (S[rr] - pr[rr]) * st$Z[rows[rr], , drop = FALSE]
      xi[rr, ] <- t(cp04_solve(f$H_chol, t(score)))
    }
    alpha[[sj]] <- xi; probabilities[[sj]] <- pr; responses[[sj]] <- S
  }
  models <- vector("list", length(backend$treatments))
  for (tj in tids) {
    spec <- backend$treatment_specs[[tj]]; f <- backend$treatments[[tj]]
    u <- length(spec$cohorts) * (d + 1L)
    upstream <- matrix(0, B, u)
    S <- pp <- numeric(B)
    for (k in seq_along(spec$cohorts)) {
      sj <- spec$selection_ids[k]; rr <- which(ci == spec$cohorts[k])
      ix <- (k - 1L) * (d + 1L) + seq_len(d + 1L)
      upstream[, ix] <- cbind(a_if[, spec$cohorts[k]], alpha[[sj]])
      S[rr] <- responses[[sj]][rr]; pp[rr] <- probabilities[[sj]][rr]
    }
    rr <- which(S == 1)
    beta_score <- matrix(0, B, p); direct <- matrix(0, B, 4L)
    direct[, 3L] <- -f$means["dp"]; direct[, 4L] <- -f$means["dm"]
    if (length(rr)) {
      X <- frame$Z_a[st$a_index[rows[rr]], , drop = FALSE]
      hit <- if (is.null(pair_op)) cp04_pair_rows(reader, rows[rr], spec$pair, TRUE) else
        pair_op$pair(rows[rr], spec$pair, TRUE)
      lv <- ipw_logit_values(X, f$coefficients)
      cp04_finite(c(lv$prob, lv$eta, lv$odds), "treatment replay")
      D <- as.double(ci[rr] == spec$cohorts[1L])
      fw <- w[rr] * st$a[ci[rr]] / (st$rho[rows[rr]] * ipw_clip_floor(pp[rr]))
      wp <- fw * D; wm <- fw * (1 - D) * lv$odds
      beta_score[rr, ] <- fw * (D - lv$prob) * X
      direct[rr, 1L] <- wp * (hit$dy - f$means["T"]) / f$means["dp"]
      direct[rr, 2L] <- wm * (hit$dy - f$means["M"]) / f$means["dm"]
      direct[rr, 3L] <- direct[rr, 3L] + wp
      direct[rr, 4L] <- direct[rr, 4L] + wm
    }
    xb <- t(cp04_solve(f$H_chol, t(beta_score + upstream %*% t(f$cross))))
    coord <- direct + upstream %*% t(f$partials$upstream) + xb %*% t(f$partials$beta)
    colnames(coord) <- c("T", "M", "dp", "dm")
    cp04_finite(c(xb, coord), "complete frame IF")
    models[[tj]] <- if (details) list(beta = xb, coordinates = coord, direct = direct,
      beta_score = beta_score, upstream = upstream,
      normalizers_A = coord[, 3:4, drop = FALSE] / st$lambda -
        outer(A - st$lambda, f$means[c("dp", "dm")] / st$lambda^2)) else
      list(phi = coord[, 1L] - coord[, 2L])
  }
  phi <- matrix(0, B, length(links), dimnames = list(NULL,
    vapply(backend$links[links], `[[`, character(1), "link_id")))
  for (j in seq_along(links)) {
    m <- models[[backend$links[[links[j]]]$model]]
    phi[, j] <- if (details) m$coordinates[, 1L] - m$coordinates[, 2L] else m$phi
  }
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi,
    scale = scale, normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_selection_score_evaluations = length(sids),
      n_treatment_score_evaluations = length(tids), output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (details) {
    # All diagnostic coordinates are on the frame scale, even with core phi.
    tau <- eta <- matrix(0, B, length(st$cohorts))
    for (k in seq_along(st$cohorts)) {
      tau[, k] <- w * (as.double(ci == k) - st$tau[k]) / st$mean_w
      eta[, k] <- w * A * (as.double(ci == k) - st$eta[k]) / st$mean_wA
    }
    out$details <- list(a = a_if, alpha = alpha, models = models,
      tau = tau, eta = eta, lambda = A - st$lambda)
  }
  actual_output <- unname(cp04_payload(out)["entries"])
  if (actual_output > pre$output) cp04_abort("IF output bound underestimated.", "bound")
  if (!is.null(pair_op)) pair_op$finish()
  out
}

