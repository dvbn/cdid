# Optional automatic baseline separate incomplete-panel fit optimization.
# Owner reads and scalar moment additions retain their original chunk order;
# only unaggregated QR input rows receive a separately admitted larger buffer.

or_baseline_qr_buffer_plan <- function(control, base, bound, scratch, p, n_A) {
  cap <- control$max_matrix_entries
  used <- ipw_owner_safe_sum(c(base, bound, scratch))
  fixed <- 128 * (p + 2)^2 + 4096
  per_row <- 64 * (p + 4)
  rows <- min(n_A, chained_public_carried_settings()$chunk_rows_cap,
    floor((cap - used - fixed) / per_row))
  if (!is.finite(rows) || rows <= control$max_chunk_rows) return(NULL)
  extra <- ipw_owner_safe_sum(c(per_row * rows, fixed))
  if (ipw_owner_safe_sum(c(used, extra)) > cap) return(NULL)
  list(rows = as.integer(rows), extra = extra)
}

# Same original per-owner-chunk validation and scalar statements as
# or_pooled_support_add. QR validation is performed before storing raw rows,
# including the augmented weighted-row overflow check of or_wls_compress.
or_baseline_qr_buffer_moments <- function(acc, recipe, Z, dy, weights) {
  n <- length(dy)
  if (!n) return(list(acc = acc, D = NULL))
  if (!is.matrix(Z) || nrow(Z) != n || ncol(Z) != length(recipe$columns) || length(weights) != n)
    or_pooled_abort("Observed pair rows disagree with their design rows.", "integrity", "baseline")
  if (recipe$intercept && any(Z[, recipe$intercept_index] != 1))
    or_pooled_abort("The declared intercept column must equal one on every observed pair row.",
      "input", "baseline")
  z <- Z[, recipe$slope_index, drop = FALSE]
  D <- acc$L * z
  dimnames(D) <- if (recipe$p) list(NULL, recipe$slope_columns)
  if (recipe$intercept) {
    if (!is.matrix(D) || !is.numeric(D) || is.object(D) || any(!is.finite(D)) ||
        !is.numeric(dy) || is.object(dy) || !is.null(dim(dy)) || length(dy) != nrow(D) ||
        any(!is.finite(dy)) || !is.numeric(weights) || is.object(weights) || !is.null(dim(weights)) ||
        length(weights) != nrow(D) || any(!is.finite(weights) | weights <= 0))
      cdid_abort("Require finite numeric D, y and strictly positive row weights with matching sizes.",
        "cdid_or_wls_input_error")
  } else {
    if (!is.matrix(D) || !is.numeric(D) || is.object(D) || ncol(D) < 1L || any(!is.finite(D)) ||
        !is.numeric(dy) || is.object(dy) || !is.null(dim(dy)) || length(dy) != nrow(D) ||
        any(!is.finite(dy)) || !is.numeric(weights) || is.object(weights) || !is.null(dim(weights)) ||
        length(weights) != nrow(D) || any(!is.finite(weights) | weights <= 0))
      cdid_abort("Require finite numeric Z, delta and strictly positive row weights with matching sizes.",
        "cdid_or_wls_input_error")
  }
  augmented <- if (recipe$intercept) cbind(1, unname(D), dy) else cbind(D, dy)
  if (any(!is.finite(augmented * sqrt(weights))))
    cdid_abort("Nonfinite weighted WLS rows or merged factor.", "cdid_or_wls_numeric_error")
  acc$weight_mass <- acc$weight_mass + sum(weights); acc$n <- acc$n + n
  acc$sumZ <- acc$sumZ + unname(colSums(weights * z)); acc$sumY <- acc$sumY + sum(weights * dy)
  if (recipe$p) acc$colmax <- pmax(acc$colmax, unname(apply(abs(z), 2L, max)))
  list(acc = acc, D = D)
}

or_baseline_qr_buffered_support <- function(reader, spec, recipe, control, buffer_rows, or_frame_index = NULL) {
  st <- reader$state
  acc <- or_pooled_support_start(recipe, spec$transitions)
  repeat_mass <- 0; n_repeat <- 0
  D_buffer <- matrix(0, buffer_rows, recipe$p,
    dimnames = if (recipe$p) list(NULL, recipe$slope_columns))
  dy_buffer <- weight_buffer <- numeric(buffer_rows)
  filled <- 0L
  flush <- function() {
    if (!filled) return(invisible(NULL))
    ix <- seq_len(filled)
    D <- D_buffer[ix, , drop = FALSE]; dy <- dy_buffer[ix]; weights <- weight_buffer[ix]
    acc$summary <<- if (recipe$intercept)
      or_pooled_wls_merge_chunks(list(acc$summary, or_pooled_wls_block_summary(D, dy, weights))) else
      or_wls_merge(list(acc$summary, or_wls_summary(D, dy, weights)))
    filled <<- 0L
    invisible(NULL)
  }
  visit <- function(rows, S, obs, dy, Z) {
    A <- st$a_index[rows] > 0L
    w <- st$w[rows]
    repeat_mass <<- repeat_mass + sum(w[A]); n_repeat <<- n_repeat + sum(A)
    added <- or_baseline_qr_buffer_moments(acc, recipe, Z, dy, w[obs])
    acc <<- added$acc
    n <- length(dy); first <- 1L
    while (first <= n) {
      take <- min(n - first + 1L, buffer_rows - filled)
      source <- seq.int(first, first + take - 1L); dest <- seq.int(filled + 1L, filled + take)
      if (recipe$p) D_buffer[dest, ] <<- added$D[source, , drop = FALSE]
      dy_buffer[dest] <<- dy[source]; weight_buffer[dest] <<- w[obs[source]]
      filled <<- filled + take; first <- first + take
      if (filled == buffer_rows) flush()
    }
    invisible(NULL)
  }
  if (is.null(or_frame_index)) or_baseline_visit(reader, spec, control, visit) else
    dr_support_visit(reader, spec, control, or_frame_index, visit)
  flush()
  or_baseline_support_rate(acc, repeat_mass, n_repeat, st$N)
}
