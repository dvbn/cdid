# Private model-agnostic shared-slope, block-intercept WLS kernel.
# Minimizes sum_b sum_i w_bi (y_bi - alpha_b - d_bi' beta)^2 over local block
# intercepts alpha_b and one common slope vector beta using small augmented QR
# factors, per-block intercept partialling, and pairwise residual factor merges.

or_pooled_wls_block_summary <- function(D, y, weights) {
  if (!is.matrix(D) || !is.numeric(D) || is.object(D) ||
      any(!is.finite(D)) || !is.numeric(y) || is.object(y) ||
      !is.null(dim(y)) || length(y) != nrow(D) ||
      any(!is.finite(y)) || !is.numeric(weights) || is.object(weights) ||
      !is.null(dim(weights)) || length(weights) != nrow(D) ||
      any(!is.finite(weights) | weights <= 0)) {
    cdid_abort(
      "Require finite numeric D, y and strictly positive row weights with matching sizes.",
      "cdid_or_wls_input_error"
    )
  }
  p <- ncol(D)
  raw <- or_wls_summary(cbind(rep.int(1, nrow(D)), unname(D)), y, weights)
  # Exact support facts come from original rows, not projection roundoff.
  # Empty chunks are neutral; only O(p) facts persist beside the QR factor.
  anchor <- if (nrow(D) > 0L) as.numeric(D[1L, ]) else numeric(0)
  constant <- if (nrow(D) > 0L) {
    vapply(seq_len(p), function(j) all(D[, j] == anchor[j]), logical(1))
  } else {
    rep.int(TRUE, p)
  }
  structure(
    list(
      raw = raw,
      n = nrow(D),
      p = p,
      columns = colnames(D),
      anchor = anchor,
      constant = constant,
      accounting = list(
        input_rows = nrow(D),
        augmented_cols = p + 2L,
        augmented_rows = nrow(raw$factor),
        max_chunk_merge_rows = nrow(raw$factor),
        n_compressions = if (nrow(D) > 0L) 1L else 0L
      )
    ),
    class = "cdid_or_pooled_wls_chunk"
  )
}

or_pooled_wls_check_chunk <- function(chunk) {
  if (!inherits(chunk, "cdid_or_pooled_wls_chunk") || !is.list(chunk) ||
      !is.numeric(chunk$p) || length(chunk$p) != 1L ||
      !is.finite(chunk$p) || chunk$p < 0L || chunk$p != floor(chunk$p)) {
    cdid_abort("Expected a compact augmented block WLS chunk.",
               "cdid_or_wls_input_error")
  }
  or_wls_check_summary(chunk$raw)
  p <- as.integer(chunk$p)
  if (ncol(chunk$raw$factor) != p + 2L ||
      !identical(chunk$n, chunk$raw$n) ||
      !is.logical(chunk$constant) || length(chunk$constant) != p ||
      anyNA(chunk$constant) ||
      !is.numeric(chunk$anchor) || is.object(chunk$anchor) ||
      !is.null(dim(chunk$anchor)) ||
      length(chunk$anchor) != (if (chunk$n > 0L) p else 0L) ||
      any(!is.finite(chunk$anchor)) ||
      (!is.null(chunk$columns) &&
       (!is.character(chunk$columns) || length(chunk$columns) != p))) {
    cdid_abort("Expected a compact augmented block WLS chunk.",
               "cdid_or_wls_input_error")
  }
  invisible(chunk)
}

or_pooled_wls_merge_chunks <- function(chunks) {
  if (!is.list(chunks) || length(chunks) == 0L) {
    cdid_abort("Supply a nonempty list of augmented block WLS chunks.",
               "cdid_or_wls_input_error")
  }
  acc <- chunks[[1L]]
  or_pooled_wls_check_chunk(acc)
  max_merge_rows <- acc$accounting$max_chunk_merge_rows
  n_compressions <- acc$accounting$n_compressions
  for (i in seq_along(chunks)[-1L]) {
    next_chunk <- chunks[[i]]
    or_pooled_wls_check_chunk(next_chunk)
    if (acc$p != next_chunk$p || !identical(acc$columns, next_chunk$columns)) {
      cdid_abort("Block WLS chunks must have the same predictor order and labels.",
                 "cdid_or_wls_input_error")
    }
    stacked_rows <- nrow(acc$raw$factor) + nrow(next_chunk$raw$factor)
    if (stacked_rows > max_merge_rows) {
      max_merge_rows <- stacked_rows
    }
    if (next_chunk$n > 0L) {
      if (acc$n == 0L) {
        acc$anchor <- next_chunk$anchor
        acc$constant <- next_chunk$constant
      } else {
        acc$constant <- acc$constant & next_chunk$constant &
          acc$anchor == next_chunk$anchor
      }
    }
    acc$raw <- or_wls_merge(list(acc$raw, next_chunk$raw))
    acc$n <- acc$raw$n
    n_compressions <- n_compressions + next_chunk$accounting$n_compressions +
      if (stacked_rows > 0L) 1L else 0L
  }
  acc$accounting <- list(
    input_rows = acc$n,
    augmented_cols = acc$p + 2L,
    augmented_rows = nrow(acc$raw$factor),
    max_chunk_merge_rows = max_merge_rows,
    n_compressions = n_compressions
  )
  acc
}

or_pooled_wls_partial_block <- function(chunk) {
  or_pooled_wls_check_chunk(chunk)
  p <- as.integer(chunk$p)
  factor <- chunk$raw$factor
  if (chunk$n == 0L || nrow(factor) == 0L) {
    cdid_abort("Empty intercept block (rank 0 < 1); model unchanged.",
               "cdid_or_wls_rank_error")
  }
  # Partial out the block intercept column on the small (<= p + 2)-row factor.
  q_int <- qr(factor[, 1L, drop = FALSE], LAPACK = TRUE)
  r11 <- qr.R(q_int)[1L, 1L]
  if (!is.finite(r11)) {
    cdid_abort("Nonfinite block intercept QR factor.",
               "cdid_or_wls_numeric_error")
  }
  if (abs(r11) <= 0) {
    cdid_abort("Empty intercept block (rank 0 < 1); model unchanged.",
               "cdid_or_wls_rank_error")
  }
  qty <- qr.qty(q_int, factor[, -1L, drop = FALSE])
  if (any(!is.finite(qty))) {
    cdid_abort("Nonfinite partialled block WLS factor.",
               "cdid_or_wls_numeric_error")
  }
  r1d <- if (p > 0L) qty[1L, seq_len(p)] else numeric(0)
  names(r1d) <- chunk$columns
  r1y <- qty[1L, p + 1L]
  residual_raw <- qty[-1L, , drop = FALSE]
  if (p > 0L && any(chunk$constant)) {
    # A truly block-constant predictor has exactly zero centered support.
    # Keep the column: another block can identify its common coefficient.
    residual_raw[, which(chunk$constant)] <- 0
  }
  residual_factor <- or_wls_compress(residual_raw)
  residual <- if (p > 0L) {
    structure(
      list(factor = residual_factor, n = chunk$n, columns = chunk$columns),
      class = "cdid_or_wls_summary"
    )
  } else {
    list(factor = residual_factor, n = chunk$n, columns = chunk$columns)
  }
  block_record <- list(
    r11 = r11,
    r1d = r1d,
    r1y = r1y,
    n = chunk$n,
    weight_sum = r11^2,
    mean_d = r1d / r11,
    mean_y = r1y / r11,
    augmented_rows = nrow(factor),
    residual_rows = nrow(residual_factor)
  )
  structure(
    list(
      residual = residual,
      blocks = list(block_record),
      n = chunk$n,
      p = p,
      columns = chunk$columns,
      accounting = list(
        n = chunk$n,
        p = p,
        n_blocks = 1L,
        residual_dim = dim(residual_factor),
        max_augmented_rows = nrow(factor),
        max_residual_rows = nrow(residual_factor),
        max_chunk_merge_rows = chunk$accounting$max_chunk_merge_rows,
        max_residual_merge_rows = nrow(residual_factor),
        n_compressions = chunk$accounting$n_compressions +
          if (nrow(residual_raw) > 0L) 1L else 0L
      )
    ),
    class = "cdid_or_pooled_wls_summary"
  )
}

or_pooled_wls_block <- function(D, y, weights) {
  or_pooled_wls_partial_block(or_pooled_wls_block_summary(D, y, weights))
}

or_pooled_wls_check_summary <- function(summary) {
  if (!inherits(summary, "cdid_or_pooled_wls_summary") || !is.list(summary) ||
      !is.numeric(summary$p) || length(summary$p) != 1L ||
      !is.finite(summary$p) || summary$p < 0L || summary$p != floor(summary$p) ||
      !is.numeric(summary$n) || length(summary$n) != 1L ||
      !is.finite(summary$n) || summary$n < 1L || summary$n != floor(summary$n) ||
      !is.list(summary$blocks) || length(summary$blocks) < 1L) {
    cdid_abort("Expected a compact pooled WLS summary.",
               "cdid_or_wls_input_error")
  }
  p <- as.integer(summary$p)
  if (!is.null(summary$columns) &&
      (!is.character(summary$columns) || length(summary$columns) != p)) {
    cdid_abort("Expected a compact pooled WLS summary.",
               "cdid_or_wls_input_error")
  }
  if (p > 0L) {
    or_wls_check_summary(summary$residual)
    if (ncol(summary$residual$factor) != p + 1L ||
        !identical(summary$residual$columns, summary$columns) ||
        !identical(summary$residual$n, summary$n)) {
      cdid_abort("Expected a compact pooled WLS summary.",
                 "cdid_or_wls_input_error")
    }
  } else {
    if (!is.list(summary$residual) || !is.matrix(summary$residual$factor) ||
        !is.numeric(summary$residual$factor) ||
        ncol(summary$residual$factor) != 1L ||
        nrow(summary$residual$factor) > 1L ||
        any(!is.finite(summary$residual$factor)) ||
        !identical(summary$residual$n, summary$n)) {
      cdid_abort("Expected a compact pooled WLS summary.",
                 "cdid_or_wls_input_error")
    }
  }
  total_n <- 0L
  for (b in summary$blocks) {
    if (!is.list(b) || !is.numeric(b$r11) || length(b$r11) != 1L ||
        !is.finite(b$r11) || abs(b$r11) <= 0 ||
        !is.numeric(b$r1d) || length(b$r1d) != p || any(!is.finite(b$r1d)) ||
        !is.numeric(b$r1y) || length(b$r1y) != 1L || !is.finite(b$r1y) ||
        !is.numeric(b$n) || length(b$n) != 1L || !is.finite(b$n) ||
        b$n < 1L || b$n != floor(b$n)) {
      cdid_abort("Expected valid partialled block intercept summaries.",
                 "cdid_or_wls_input_error")
    }
    total_n <- total_n + as.integer(b$n)
  }
  if (total_n != summary$n) {
    cdid_abort("Block row counts do not match pooled summary total.",
               "cdid_or_wls_input_error")
  }
  invisible(summary)
}

or_pooled_wls_name_blocks <- function(summaries) {
  outer_names <- names(summaries)
  if (!is.null(outer_names) && any(is.na(outer_names))) {
    cdid_abort("Block names must not be NA.", "cdid_or_wls_input_error")
  }
  blocks <- list()
  for (i in seq_along(summaries)) {
    s_blocks <- summaries[[i]]$blocks
    if (!is.null(outer_names) && nzchar(outer_names[i])) {
      if (length(s_blocks) == 1L && is.null(names(s_blocks))) {
        names(s_blocks) <- outer_names[i]
      }
    }
    blocks <- c(blocks, s_blocks)
  }
  block_names <- names(blocks)
  if (!is.null(block_names)) {
    if (any(is.na(block_names) | !nzchar(block_names)) ||
        anyDuplicated(block_names) > 0L) {
      cdid_abort("Supplied block names must be non-empty and unique.",
                 "cdid_or_wls_input_error")
    }
  }
  blocks
}

or_pooled_wls_merge <- function(summaries) {
  if (!is.list(summaries) || length(summaries) == 0L) {
    cdid_abort("Supply a nonempty list of pooled WLS summaries.",
               "cdid_or_wls_input_error")
  }
  for (s in summaries) {
    or_pooled_wls_check_summary(s)
  }
  acc <- summaries[[1L]]
  p <- as.integer(acc$p)
  max_aug <- acc$accounting$max_augmented_rows
  max_res <- acc$accounting$max_residual_rows
  max_chunk_merge <- acc$accounting$max_chunk_merge_rows
  max_res_merge <- acc$accounting$max_residual_merge_rows
  n_compressions <- acc$accounting$n_compressions
  for (i in seq_along(summaries)[-1L]) {
    next_summary <- summaries[[i]]
    if (p != next_summary$p || !identical(acc$columns, next_summary$columns)) {
      cdid_abort("Pooled WLS summaries must have the same predictor order and labels.",
                 "cdid_or_wls_input_error")
    }
    stacked_rows <- nrow(acc$residual$factor) + nrow(next_summary$residual$factor)
    if (stacked_rows > max_res_merge) {
      max_res_merge <- stacked_rows
    }
    if (p > 0L) {
      acc$residual <- or_wls_merge(list(acc$residual, next_summary$residual))
    } else {
      acc$residual$factor <- or_wls_compress(
        rbind(acc$residual$factor, next_summary$residual$factor)
      )
      acc$residual$n <- acc$residual$n + next_summary$residual$n
    }
    acc$n <- acc$n + next_summary$n
    max_aug <- max(max_aug, next_summary$accounting$max_augmented_rows)
    max_res <- max(max_res, next_summary$accounting$max_residual_rows)
    max_chunk_merge <- max(max_chunk_merge,
                           next_summary$accounting$max_chunk_merge_rows)
    max_res_merge <- max(max_res_merge,
                         next_summary$accounting$max_residual_merge_rows)
    n_compressions <- n_compressions + next_summary$accounting$n_compressions +
      if (stacked_rows > 0L) 1L else 0L
  }
  acc$blocks <- or_pooled_wls_name_blocks(summaries)
  acc$accounting <- list(
    n = acc$n,
    p = p,
    n_blocks = length(acc$blocks),
    residual_dim = dim(acc$residual$factor),
    max_augmented_rows = max_aug,
    max_residual_rows = max_res,
    max_chunk_merge_rows = max_chunk_merge,
    max_residual_merge_rows = max_res_merge,
    n_compressions = n_compressions
  )
  acc
}

or_pooled_wls_fit <- function(summary, rank_tol) {
  or_pooled_wls_check_summary(summary)
  if (missing(rank_tol) || !is.numeric(rank_tol) || is.object(rank_tol) ||
      !is.null(dim(rank_tol)) || length(rank_tol) != 1L ||
      !is.finite(rank_tol) || rank_tol <= 0 || rank_tol >= 1) {
    cdid_abort("Supply explicit finite private rank_tol in (0, 1).",
               "cdid_or_wls_input_error")
  }
  p <- as.integer(summary$p)
  n_blocks <- length(summary$blocks)
  if (p > 0L) {
    m <- nrow(summary$residual$factor)
    if (m == 0L) {
      cdid_abort(sprintf(
        "Final pooled WLS slope rank 0 < %d (n=%g, blocks=%d); model unchanged.",
        p, summary$n, n_blocks), "cdid_or_wls_rank_error")
    }
    q <- qr(summary$residual$factor[, seq_len(p), drop = FALSE], LAPACK = TRUE)
    R <- qr.R(q)
    k_diag <- min(nrow(R), ncol(R))
    diagonal <- abs(R[cbind(seq_len(k_diag), seq_len(k_diag))])
    cutoff <- rank_tol * max(diagonal)
    rank <- if (cutoff > 0) sum(diagonal > cutoff) else 0L
    if (m < p || rank < p) {
      cdid_abort(sprintf(
        "Final pooled WLS slope rank %d < %d (n=%g, blocks=%d, rank_tol=%g, cutoff=%g); model unchanged.",
        rank, p, summary$n, n_blocks, rank_tol, cutoff),
        "cdid_or_wls_rank_error")
    }
    slope_fit <- or_wls_fit(summary$residual, rank_tol = rank_tol)
    slope_fit$cutoff <- cutoff
    coefficients <- slope_fit$coefficients
    rss <- slope_fit$rss
    R_out <- slope_fit$R
    pivot <- slope_fit$pivot
  } else {
    slope_fit <- NULL
    coefficients <- numeric(0)
    if (!is.null(summary$columns)) {
      names(coefficients) <- summary$columns
    }
    rss <- sum(summary$residual$factor[, 1L]^2)
    R_out <- matrix(numeric(0), 0L, 0L)
    pivot <- integer(0)
    rank <- 0L
    cutoff <- 0
  }
  # Recover local block intercepts via 1x1 triangular backsolve on each block's
  # intercept QR row: r11 * alpha_b + r1d' * beta = r1y.
  intercepts <- vapply(summary$blocks, function(b) {
    rhs <- b$r1y - sum(b$r1d * coefficients)
    as.numeric(backsolve(matrix(b$r11, 1L, 1L), rhs))
  }, numeric(1))
  names(intercepts) <- names(summary$blocks)
  if (any(!is.finite(c(intercepts, coefficients, rss, R_out)))) {
    cdid_abort("Nonfinite final pooled WLS solution.",
               "cdid_or_wls_numeric_error")
  }
  structure(
    list(
      coefficients = coefficients,
      intercepts = intercepts,
      rss = rss,
      R = R_out,
      pivot = pivot,
      n = summary$n,
      p = p,
      n_blocks = n_blocks,
      rank = rank,
      rank_tol = rank_tol,
      cutoff = cutoff,
      blocks = summary$blocks,
      slope_fit = slope_fit,
      accounting = summary$accounting
    ),
    class = "cdid_or_pooled_wls_fit"
  )
}

or_pooled_wls_gram_solve <- function(fit, b, normalizer) {
  if (!inherits(fit, "cdid_or_pooled_wls_fit")) {
    cdid_abort("Expected a regular full-rank pooled WLS fit.",
               "cdid_or_wls_input_error")
  }
  p <- as.integer(fit$p)
  if (p > 0L) {
    return(or_wls_gram_solve(fit$slope_fit, b, normalizer = normalizer))
  }
  vector_b <- is.null(dim(b))
  if (!is.numeric(b) || is.object(b) ||
      !(vector_b || is.matrix(b)) || any(!is.finite(b)) ||
      (vector_b && length(b) != 0L) || (!vector_b && nrow(b) != 0L) ||
      missing(normalizer) || !is.numeric(normalizer) || is.object(normalizer) ||
      !is.null(dim(normalizer)) || length(normalizer) != 1L ||
      !is.finite(normalizer) || normalizer <= 0) {
    cdid_abort("Require finite p-row b and an explicit positive normalizer.",
               "cdid_or_wls_input_error")
  }
  b * normalizer
}

or_pooled_wls_objective <- function(summary, beta, intercepts = NULL) {
  or_pooled_wls_check_summary(summary)
  p <- as.integer(summary$p)
  if (!is.numeric(beta) || is.object(beta) || !is.null(dim(beta)) ||
      length(beta) != p || any(!is.finite(beta))) {
    cdid_abort("Require a finite numeric slope vector matching summary$p.",
               "cdid_or_wls_input_error")
  }
  factor <- summary$residual$factor
  profiled <- if (p > 0L) {
    sum((factor[, p + 1L] - factor[, seq_len(p), drop = FALSE] %*% beta)^2)
  } else {
    sum(factor[, 1L]^2)
  }
  if (is.null(intercepts)) {
    return(profiled)
  }
  n_blocks <- length(summary$blocks)
  if (!is.numeric(intercepts) || is.object(intercepts) ||
      !is.null(dim(intercepts)) || length(intercepts) != n_blocks ||
      any(!is.finite(intercepts))) {
    cdid_abort("Require finite numeric intercepts matching the number of blocks.",
               "cdid_or_wls_input_error")
  }
  intercept_penalty <- sum(vapply(seq_len(n_blocks), function(j) {
    b <- summary$blocks[[j]]
    (b$r1y - b$r11 * intercepts[j] - sum(b$r1d * beta))^2
  }, numeric(1)))
  profiled + intercept_penalty
}

or_pooled_wls_predict <- function(fit, D_new, block = 1L) {
  if (!inherits(fit, "cdid_or_pooled_wls_fit")) {
    cdid_abort("Expected a regular full-rank pooled WLS fit.",
               "cdid_or_wls_input_error")
  }
  p <- as.integer(fit$p)
  if (!is.matrix(D_new) || !is.numeric(D_new) || is.object(D_new) ||
      ncol(D_new) != p || any(!is.finite(D_new))) {
    cdid_abort("Require a finite numeric D_new matrix with p columns.",
               "cdid_or_wls_input_error")
  }
  n_new <- nrow(D_new)
  if (is.character(block)) {
    if (is.null(names(fit$intercepts)) || any(is.na(block)) ||
        !all(block %in% names(fit$intercepts))) {
      cdid_abort("Unknown block name in pooled WLS prediction.",
                 "cdid_or_wls_input_error")
    }
    alpha <- unname(fit$intercepts[block])
  } else if (is.numeric(block) && !is.object(block) && is.null(dim(block)) &&
             all(is.finite(block)) && all(block == floor(block)) &&
             all(block >= 1L & block <= fit$n_blocks)) {
    alpha <- unname(fit$intercepts[as.integer(block)])
  } else {
    cdid_abort("Invalid block index or name in pooled WLS prediction.",
               "cdid_or_wls_input_error")
  }
  if (!(length(alpha) == 1L || length(alpha) == n_new)) {
    cdid_abort("Block specifier must have length 1 or nrow(D_new).",
               "cdid_or_wls_input_error")
  }
  slope_part <- if (p > 0L) drop(D_new %*% fit$coefficients) else rep(0, n_new)
  as.numeric(alpha + slope_part)
}
