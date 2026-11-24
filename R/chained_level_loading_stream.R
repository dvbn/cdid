# Streamed signed level loadings for `gmm_weighting = "within"`.
#
# The general cohort producers record every signed loading contribution in
# program order instead of writing dense per-date `n_a x L` matrices. A record
# adds its values at date `tp` and, for a pair difference, subtracts the same
# values at date `tn`. Values are regenerated from the factors the producer
# already formed, with the same scalar operations: a vector, `outer(a, b)`, or
# `left %*% right`, each with an optional source-row selection and row
# multiplier. Replaying the records of one date on one block of A rows, in
# program order, reproduces every entry of the dense accumulation exactly.

# `limit` bounds the retained log entries; a record that would exceed it is
# refused with the established workspace error before it is retained.
level_loading_log <- function(n_cal, n_a, L, limit = Inf) {
  e <- new.env(parent = emptyenv())
  e$limit <- as.double(limit)
  e$n_cal <- as.integer(n_cal)
  e$n_a <- as.integer(n_a)
  e$L <- as.integer(L)
  e$recs <- vector("list", 256L)
  e$n <- 0L
  e$entries <- 0
  e$by_tau <- NULL
  e
}

loading_log_push <- function(log, rec) {
  if (!length(rec$rows) || !length(rec$cols)) return(invisible(NULL))
  n <- log$n + 1L
  if (n > length(log$recs)) length(log$recs) <- 2L * length(log$recs)
  log$recs[[n]] <- rec
  log$n <- n
  log$by_tau <- NULL
  sz <- length(rec$rows) + length(rec$cols) + length(rec$val) + length(rec$a) +
    length(rec$b) + length(rec$src) + length(rec$pre) + length(rec$Lf) + length(rec$R)
  if (is.finite(log$limit)) ipw_owner_safe_sum(c(log$entries, sz), log$limit, "CP04 within loading log entries")
  log$entries <- log$entries + sz
  invisible(NULL)
}

# Pair difference of a single link column: `+val` at `tp`, `-val` at `tn`.
loading_log_pair_vec <- function(log, rows, tp, tn, col, val) {
  loading_log_push(log, list(kind = "v", tp = as.integer(tp), tn = as.integer(tn),
    rows = as.integer(rows), cols = as.integer(col), val = as.double(val)))
}

# Pair difference of `pre * outer(a, b)[src, ]` on link columns `cols`.
loading_log_pair_outer <- function(log, rows, tp, tn, cols, a, b, src = NULL, pre = NULL) {
  loading_log_push(log, list(kind = "o", tp = as.integer(tp), tn = as.integer(tn),
    rows = as.integer(rows), cols = as.integer(cols), a = a, b = b, src = src, pre = pre))
}

# Pair difference (`tn > 0`) or single-date addition (`tn = 0`) of
# `pre * (left %*% right)[src, ]` on link columns `cols`.
loading_log_lowrank <- function(log, rows, tp, tn, cols, left, right, src = NULL, pre = NULL) {
  loading_log_push(log, list(kind = "lr", tp = as.integer(tp), tn = as.integer(tn),
    rows = as.integer(rows), cols = as.integer(cols), Lf = left, R = right, src = src, pre = pre))
}

loading_log_index <- function(log) {
  if (!is.null(log$by_tau)) return(log$by_tau)
  n <- log$n
  tp <- integer(n)
  tn <- integer(n)
  for (k in seq_len(n)) {
    tp[k] <- log$recs[[k]]$tp
    tn[k] <- log$recs[[k]]$tn
  }
  # Record order within a date is program order; the two dates of a pair
  # difference are distinct matrices, so per-date order is all that matters.
  by_tau <- vector("list", log$n_cal)
  for (tau in seq_len(log$n_cal)) {
    k_p <- which(tp == tau)
    k_n <- which(tn == tau)
    k_all <- c(k_p, k_n)
    sg <- c(rep.int(1L, length(k_p)), rep.int(-1L, length(k_n)))
    o <- order(k_all)
    by_tau[[tau]] <- list(k = k_all[o], sgn = sg[o])
  }
  log$by_tau <- by_tau
  by_tau
}

loading_rec_values <- function(rec, sel, keep_c) {
  if (identical(rec$kind, "v")) return(rec$val[sel])
  if (identical(rec$kind, "o")) {
    a <- if (is.null(rec$src)) rec$a[sel] else rec$a[rec$src[sel]]
    M <- outer(a, rec$b[keep_c])
  } else {
    Lf <- if (is.null(rec$src)) rec$Lf[sel, , drop = FALSE] else rec$Lf[rec$src[sel], , drop = FALSE]
    M <- Lf %*% rec$R[, keep_c, drop = FALSE]
  }
  if (!is.null(rec$pre)) M <- rec$pre[sel] * M
  M
}

# Dense loadings of date `tau` on A rows `r0:r1` for catalog columns mapped by
# `col_pos` (0 = not requested) into an `(r1 - r0 + 1) x n_col` matrix.
loading_log_replay <- function(log, tau, r0, r1, col_pos, n_col) {
  b <- r1 - r0 + 1L
  out <- matrix(0.0, b, n_col)
  if (b < 1L || n_col < 1L) return(out)
  idx <- loading_log_index(log)[[tau]]
  if (is.null(idx) || !length(idx$k)) return(out)
  for (ii in seq_along(idx$k)) {
    rec <- log$recs[[idx$k[ii]]]
    oc <- col_pos[rec$cols]
    keep_c <- which(oc > 0L)
    if (!length(keep_c)) next
    rws <- rec$rows
    sel <- which(rws >= r0 & rws <= r1)
    if (!length(sel)) next
    loc <- rws[sel] - (r0 - 1L)
    val <- loading_rec_values(rec, sel, keep_c)
    oc_k <- oc[keep_c]
    if (identical(rec$kind, "v")) {
      if (idx$sgn[ii] > 0L) out[loc, oc_k] <- out[loc, oc_k] + val else out[loc, oc_k] <- out[loc, oc_k] - val
    } else if (idx$sgn[ii] > 0L) {
      out[loc, oc_k] <- out[loc, oc_k, drop = FALSE] + val
    } else {
      out[loc, oc_k] <- out[loc, oc_k, drop = FALSE] - val
    }
  }
  out
}

# Dense per-date list equal to the former accumulation (reference use only).
loading_log_dense <- function(log, cols = seq_len(log$L)) {
  col_pos <- integer(log$L)
  col_pos[cols] <- seq_along(cols)
  lapply(seq_len(log$n_cal), function(tau)
    loading_log_replay(log, tau, 1L, log$n_a, col_pos, length(cols)))
}

# Smallest replay block, in rows. Every block rescans the replay records, so
# very small blocks multiply the replay cost without bounding memory usefully.
stream_min_block_rows <- 64L

# Rows per replay block. `live_entries` counts arrays that stay allocated while
# blocks are formed; `per_row` counts the doubles one block row occupies across
# all simultaneous block buffers; the row selection of one replay record adds
# at most four vectors of length `n_rows`. Per block row, a replay holds its
# output, one record's values and the selected entries with their update
# (4 columns per replayed column); the within Gram holds the replayed block
# and a cohort's column subset, its nonzero mask or row subset and scaled
# rows (4 L); the witness holds the replayed block, a group's tree and
# candidate columns, the predicted and observed candidates and the
# difference, its absolute value and the transposed copy formed by `apply()`
# (7 M); the optimal Gram holds a row subset of the centered rows and its
# scaled copy, charged as 3 L. Buffers of the previous block are released before the
# next one is formed. The callers of the streamed Grams count as live, beside
# the log and the influence-function rows, 5 L^2: the Gram accumulators, a
# block cross-product and their sum during the loop, and at the end the
# accumulators, their transpose, sum and halved copy and the result. The block never exceeds `cap` or `n_rows` rows and is
# never refused. It holds the rows that fit in the workspace left
# after the live arrays, but at least `stream_min_block_rows` rows (or all
# rows, when fewer remain). The block buffers can therefore exceed the
# remaining counted workspace (zero when the live arrays already exceed it) by
# at most `stream_min_block_rows * per_row` entries.
stream_block_rows <- function(n_rows, per_row, max_entries, live_entries, cap = n_rows) {
  n_rows <- as.integer(n_rows)
  if (n_rows < 1L) return(1L)
  cap <- as.integer(max(1, cap))
  per_row <- max(1.0, as.double(per_row))
  spare <- as.double(max_entries) - as.double(live_entries) - 4.0 * n_rows
  b <- if (is.finite(spare)) floor(spare / per_row) else n_rows
  as.integer(max(min(n_rows, cap, stream_min_block_rows), min(n_rows, cap, b)))
}

# Streamed loading witness for the post-fit forest reduction. `checks` holds,
# per group, the tree columns, surviving candidate columns and their tree-path
# coefficients. Returns, per group, the logical vector of candidates whose
# loadings fail to telescope at some date. Per date and per candidate column
# the maxima equal those of the dense check (maxima are exact and rows without
# loadings contribute zero), and thresholds are applied per date as before.
loading_log_witness <- function(log, checks, tol, block_rows) {
  bad <- lapply(checks, function(ch) logical(length(ch$cand)))
  if (!length(checks)) return(bad)
  need <- sort(unique(unlist(lapply(checks, function(ch) c(ch$tree, ch$cand)), use.names = FALSE)))
  col_pos <- integer(log$L)
  col_pos[need] <- seq_along(need)
  n_a <- log$n_a
  starts <- seq.int(1L, max(1L, n_a), by = max(1L, block_rows))
  for (tau in seq_len(log$n_cal)) {
    acc <- lapply(checks, function(ch) {
      d <- length(ch$cand)
      list(has = logical(d), mx_ca = numeric(d), mx_pr = numeric(d), mx_err = numeric(d), any_nz = FALSE)
    })
    for (r0 in starts) {
      r1 <- min(n_a, r0 + block_rows - 1L)
      if (r1 < r0) next
      # Release the previous block's buffers before the next replay.
      U_b <- U_tr <- U_ca <- U_pred <- U_ca_nz <- NULL
      U_b <- loading_log_replay(log, tau, r0, r1, col_pos, length(need))
      for (gi in seq_along(checks)) {
        ch <- checks[[gi]]
        U_tr <- U_b[, col_pos[ch$tree], drop = FALSE]
        U_ca <- U_b[, col_pos[ch$cand], drop = FALSE]
        nz_u <- which(rowSums(U_tr != 0) > 0L | rowSums(U_ca != 0) > 0L)
        if (!length(nz_u)) next
        U_pred <- U_tr[nz_u, , drop = FALSE] %*% ch$E
        U_ca_nz <- U_ca[nz_u, , drop = FALSE]
        a <- acc[[gi]]
        a$any_nz <- TRUE
        a$has <- a$has | (colSums(U_ca_nz != 0) > 0L) | (colSums(U_pred != 0) > 0L)
        a$mx_ca <- pmax(a$mx_ca, apply(abs(U_ca_nz), 2L, max))
        a$mx_pr <- pmax(a$mx_pr, apply(abs(U_pred), 2L, max))
        a$mx_err <- pmax(a$mx_err, apply(abs(U_ca_nz - U_pred), 2L, max))
        acc[[gi]] <- a
      }
    }
    for (gi in seq_along(checks)) {
      a <- acc[[gi]]
      if (!a$any_nz || !any(a$has)) next
      u_scale <- pmax(1.0, a$mx_ca, a$mx_pr)
      bad_u <- a$has & (!is.finite(a$mx_err) | a$mx_err > tol * u_scale)
      bad[[gi]] <- bad[[gi]] | bad_u
    }
  }
  bad
}

# Streamed within covariance on retained catalog columns `cols` (cohort labels
# `g_vec`). Equals `within_cov_from_level_loadings()` applied to the dense
# loadings of `cols` with the same `scale_by_g`, up to the order of the row
# sums in each cross-product.
within_cov_stream <- function(log, cols, g_vec, scale_by_g = NULL, block_rows) {
  L <- length(g_vec)
  if (L == 0L) return(matrix(0.0, 0L, 0L))
  g_num <- as.double(g_vec)
  u_g <- unique(g_num)
  cols_by_g <- split(seq_len(L), match(g_num, u_g))
  sc_list <- lapply(seq_along(u_g), function(gi) {
    if (is.null(scale_by_g)) NULL else if (is.list(scale_by_g)) scale_by_g[[as.character(u_g[gi])]] else scale_by_g
  })
  S <- lapply(cols_by_g, function(cg) matrix(0.0, length(cg), length(cg)))
  col_pos <- integer(log$L)
  col_pos[cols] <- seq_along(cols)
  n_a <- log$n_a
  starts <- seq.int(1L, max(1L, n_a), by = max(1L, block_rows))
  for (tau in seq_len(log$n_cal)) {
    for (r0 in starts) {
      r1 <- min(n_a, r0 + block_rows - 1L)
      if (r1 < r0) next
      # Release the previous block's buffers before the next replay.
      U_b <- U_sub <- U_nz <- NULL
      U_b <- loading_log_replay(log, tau, r0, r1, col_pos, L)
      for (gi in seq_along(u_g)) {
        U_sub <- U_b[, cols_by_g[[gi]], drop = FALSE]
        nz_r <- which(rowSums(U_sub != 0) > 0L)
        if (!length(nz_r)) next
        sc_g <- sc_list[[gi]]
        U_nz <- if (is.null(sc_g)) U_sub[nz_r, , drop = FALSE] else
          U_sub[nz_r, , drop = FALSE] * sc_g[r0 - 1L + nz_r]
        S[[gi]] <- S[[gi]] + crossprod(U_nz)
      }
    }
  }
  V <- matrix(0.0, L, L)
  for (gi in seq_along(u_g)) {
    cg <- cols_by_g[[gi]]
    V[cg, cg] <- 0.5 * (S[[gi]] + t(S[[gi]]))
  }
  V
}
