# Private computational summaries for ONE declared control-WLS recipe.
# `weights` are the complete fixed positive row weights (e.g. w * v), applied
# once. Extraction, pool keys, upstream fitting, links and IF replay live later.

or_wls_compress <- function(augmented) {
  if (any(!is.finite(augmented))) {
    cdid_abort("Nonfinite weighted WLS rows or merged factor.", "cdid_or_wls_numeric_error")
  }
  if (nrow(augmented) == 0L) return(augmented)
  # LAPACK's column-pivoted QR retains ALL columns, including dependent ones.
  # A[, pivot] = Q R: restore the ORIGINAL order before any subsequent merge.
  q <- qr(augmented, LAPACK = TRUE)
  factor <- qr.R(q)[, order(q$pivot), drop = FALSE]
  if (any(!is.finite(factor))) {
    cdid_abort("Nonfinite compressed WLS factor.", "cdid_or_wls_numeric_error")
  }
  unname(factor)
}

or_wls_summary <- function(Z, delta, weights) {
  if (!is.matrix(Z) || !is.numeric(Z) || is.object(Z) || ncol(Z) < 1L ||
      any(!is.finite(Z)) || !is.numeric(delta) || is.object(delta) ||
      !is.null(dim(delta)) || length(delta) != nrow(Z) ||
      any(!is.finite(delta)) || !is.numeric(weights) || is.object(weights) ||
      !is.null(dim(weights)) || length(weights) != nrow(Z) ||
      any(!is.finite(weights) | weights <= 0)) {
    cdid_abort(
      "Require finite numeric Z, delta and strictly positive row weights with matching sizes.",
      "cdid_or_wls_input_error"
    )
  }
  # Keeping the response column preserves ||sqrt(W) (delta - Z beta)||^2
  # for EVERY beta, including the residual constant, not only its minimizer.
  structure(list(factor = or_wls_compress(cbind(Z, delta) * sqrt(weights)),
                 n = nrow(Z), columns = colnames(Z)),
            class = "cdid_or_wls_summary")
}

or_wls_check_summary <- function(summary) {
  if (!inherits(summary, "cdid_or_wls_summary") || !is.list(summary) ||
      !is.matrix(summary$factor) || !is.numeric(summary$factor) ||
      ncol(summary$factor) < 2L ||
      nrow(summary$factor) > ncol(summary$factor) ||
      any(!is.finite(summary$factor)) ||
      !is.numeric(summary$n) || length(summary$n) != 1L ||
      !is.finite(summary$n) || summary$n < nrow(summary$factor) ||
      summary$n != floor(summary$n)) {
    cdid_abort("Expected a compact WLS summary.", "cdid_or_wls_input_error")
  }
  invisible(summary)
}

or_wls_merge <- function(summaries) {
  if (!is.list(summaries) || length(summaries) == 0L) {
    cdid_abort("Supply a nonempty list of WLS summaries.", "cdid_or_wls_input_error")
  }
  acc <- summaries[[1L]]
  or_wls_check_summary(acc)
  # Each QR sees at most 2(p+1) rows. No all-cohort stacked factor is built.
  for (i in seq_along(summaries)[-1L]) {
    next_summary <- summaries[[i]]
    or_wls_check_summary(next_summary)
    if (ncol(acc$factor) != ncol(next_summary$factor) ||
        !identical(acc$columns, next_summary$columns)) {
      cdid_abort("WLS summaries must have the same predictor order and labels.",
                 "cdid_or_wls_input_error")
    }
    acc$factor <- or_wls_compress(rbind(acc$factor, next_summary$factor))
    acc$n <- acc$n + next_summary$n
  }
  acc
}

or_wls_fit <- function(summary, rank_tol) {
  or_wls_check_summary(summary)
  if (missing(rank_tol) || !is.numeric(rank_tol) || is.object(rank_tol) ||
      !is.null(dim(rank_tol)) || length(rank_tol) != 1L ||
      !is.finite(rank_tol) || rank_tol <= 0 || rank_tol >= 1) {
    cdid_abort("Supply explicit finite private rank_tol in (0, 1).",
               "cdid_or_wls_input_error")
  }
  p <- ncol(summary$factor) - 1L
  m <- nrow(summary$factor)
  if (m == 0L) {
    cdid_abort(sprintf("Final WLS rank 0 < %d (empty pool); model unchanged.", p),
               "cdid_or_wls_rank_error")
  }
  q <- qr(summary$factor[, seq_len(p), drop = FALSE], LAPACK = TRUE)
  R <- qr.R(q)
  diagonal <- abs(diag(R))
  cutoff <- rank_tol * max(diagonal)
  rank <- sum(diagonal > cutoff)
  # LAPACK's `rank` slot does not diagnose deficiency. This explicit private
  # relative-diagonal criterion is applied only to the final pooled design.
  if (rank < p) {
    cdid_abort(sprintf(
      "Final WLS rank %d < %d (n=%g, rank_tol=%g, cutoff=%g); model unchanged.",
      rank, p, summary$n, rank_tol, cutoff), "cdid_or_wls_rank_error")
  }
  qty <- qr.qty(q, summary$factor[, p + 1L])
  coefficients <- numeric(p)
  coefficients[q$pivot] <- backsolve(R, qty[seq_len(p)])
  names(coefficients) <- summary$columns
  rss <- sum(qty[-seq_len(p)]^2)
  if (any(!is.finite(c(coefficients, rss, R)))) {
    cdid_abort("Nonfinite final WLS solution.", "cdid_or_wls_numeric_error")
  }
  structure(list(coefficients = coefficients, rss = rss, R = unname(R),
                 pivot = q$pivot, n = summary$n,
                 rank = rank, rank_tol = rank_tol, cutoff = cutoff),
            class = "cdid_or_wls_fit")
}

or_wls_gram_solve <- function(fit, b, normalizer) {
  if (!inherits(fit, "cdid_or_wls_fit")) {
    cdid_abort("Expected a regular full-rank WLS fit.", "cdid_or_wls_input_error")
  }
  p <- ncol(fit$R)
  vector_b <- is.null(dim(b))
  if (!is.numeric(b) || is.object(b) ||
      !(vector_b || is.matrix(b)) || any(!is.finite(b)) ||
      (vector_b && length(b) != p) || (!vector_b && nrow(b) != p) ||
      missing(normalizer) || !is.numeric(normalizer) || is.object(normalizer) ||
      !is.null(dim(normalizer)) || length(normalizer) != 1L ||
      !is.finite(normalizer) || normalizer <= 0) {
    cdid_abort("Require finite p-row b and an explicit positive normalizer.",
               "cdid_or_wls_input_error")
  }
  rhs <- if (vector_b) matrix(b, ncol = 1L) else b
  # Z' W Z in pivot order is R' R. Reuse the final factor; never invert it.
  sol <- backsolve(fit$R, forwardsolve(t(fit$R),
                   rhs[fit$pivot, , drop = FALSE])) * normalizer
  ans <- sol
  ans[fit$pivot, ] <- sol
  dimnames(ans) <- dimnames(rhs)
  if (any(!is.finite(ans))) {
    cdid_abort("Nonfinite normalized WLS Gram loading.", "cdid_or_wls_numeric_error")
  }
  if (vector_b) {
    ans <- as.vector(ans)
    names(ans) <- names(b)
  }
  ans
}
