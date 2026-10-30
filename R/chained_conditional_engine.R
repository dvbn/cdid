# CP08 O14 (D63/D66): streaming engine of the chained conditional pretrend
# CvM test. Every declared cell c is one direct auxiliary link
#   mu = sum_r s_r N_r / D_r,  N_r = H_N[nu_r rho_r],  D_r = H_N[nu_r] > 0,
# of the accepted union-catalog fit; its process inserts the frozen grid
# indicator h_u(z_c) into every numerator integrand (no nuisance refit). The
# backend extractors (R/chained_conditional_extract.R) return, per bounded row
# chunk of the full original frame F, the term weights nu_r and responses
# rho_r, the complete nuisance IF rows xi (N-normalized, upstream forcing
# included) and the per-row loading integrands
#   P_i = sum_r (s_r / D_r) nu_r,i (rho_r,i lambda_r,i + pi_r,i),
#   blam_r = sum_i nu_r,i lambda_r,i,
# with lambda_r = d log nu_r / d eta and pi_r = d rho_r / d eta. By the
# accepted factorization (S4)
#   phi_i(u) = h_u(z_i) A_i - sum_r R_r(u) kappa_r,i + xi_i' L(u),
#   A_i = sum_r s_r nu_r,i rho_r,i / D_r,  kappa_r,i = s_r nu_r,i / D_r,
#   R_r(u) = H_N[nu_r h_u rho_r] / D_r,
#   L(u) = H_N[h_u P] - sum_r (s_r / D_r) R_r(u) H_N[nu_r lambda_r].
# The grid is one row per original unit (M = N, multiplicities kept); all
# u-dependence is through indicator sums over the DISTINCT grid points, so no
# original-unit by grid IF array is formed. Multiplier draws use one B x C
# Rademacher sign matrix shared by every cell and grid chunk:
#   J*_b(u) = N^-1 [eta_b' Agg(h_u A) - sum_r R_r(u) eta_b' Agg kappa_r
#                   + (eta_b' Agg Xi) L(u)],
#   CvM = (N/M) sum_c sum_m J_c(u_m)^2,  CvM_b = (N/M) sum_c sum_m J*_bc(u_m)^2.

cpt_abort <- function(message, class = character(), diagnostics = NULL) {
  stop(structure(
    class = unique(c(class, "cdid_conditional_pretest_error", "cdid_error", "error", "condition")),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

# Distinct grid points of one recipe design (exact binary64 identity, -0 = 0),
# in first-occurrence order of the original units, with each unit's map and the
# multiplicities. M is the number of original units (duplicates keep mass).
cpt_distinct_grid <- function(Z) {
  Z <- unname(as.matrix(Z))
  storage.mode(Z) <- "double"
  if (!nrow(Z) || any(!is.finite(Z))) cpt_abort("The conditional grid design must be finite on every original unit.",
    "cdid_conditional_pretest_integrity_error")
  key <- if (ncol(Z)) do.call(paste, c(lapply(seq_len(ncol(Z)), function(j) sprintf("%a", Z[, j] + 0)), sep = "|")) else
    rep("", nrow(Z))
  first <- which(!duplicated(key))
  map <- match(key, key[first])
  list(z = Z[first, , drop = FALSE], first = first, map = map, mult = tabulate(map, length(first)),
    M = nrow(Z), n_distinct = length(first))
}

# H[i, m] = 1{zrow_i <= zcol_m componentwise}, the did orthant indicator.
cpt_indicator <- function(zrow, zcol) {
  H <- matrix(TRUE, nrow(zrow), nrow(zcol))
  for (j in seq_len(ncol(zrow))) H <- H & outer(zrow[, j], zcol[, j], "<=")
  storage.mode(H) <- "double"
  H
}

cpt_chunks <- function(n, size) {
  size <- max(1L, as.integer(size))
  if (n < 1L) return(list())
  split(seq_len(n), ceiling(seq_len(n) / size))
}

# Empty accumulator of one cell on its recipe grid. `keep_rows` additionally
# retains the original-unit rows (A, kappa, xi) for small diagnostic checks.
cpt_cell_acc <- function(spec, grid, n_clusters, n_total, keep_rows = FALSE) {
  R <- length(spec$s); Dl <- length(spec$coords); Md <- grid$n_distinct
  acc <- list(Vd = matrix(0, Md, R), Pd = matrix(0, Md, Dl), Bl = matrix(0, Dl, R),
    AggX = matrix(0, n_clusters, Dl), AggKu = matrix(0, n_clusters, R),
    A_key = numeric(0L), A_val = numeric(0L), sum_nu = numeric(R), n_rows_seen = 0)
  if (isTRUE(keep_rows)) {
    acc$rows_A <- numeric(n_total); acc$rows_nu <- matrix(0, n_total, R); acc$rows_xi <- matrix(0, n_total, Dl)
  }
  acc
}

cpt_rowsum_into <- function(target, x, group) {
  if (!length(group) || !ncol(x)) return(target)
  s <- rowsum(x, group, reorder = FALSE)
  idx <- as.integer(rownames(s))
  target[idx, ] <- target[idx, , drop = FALSE] + s
  target
}

# One extractor chunk of one cell. Rows are original F positions; nu/rho are
# rows x terms; P/xi rows x coordinates; blam coordinates x terms (chunk sums).
cpt_cell_add <- function(acc, spec, piece, grid, cmap) {
  rows <- piece$rows
  if (!length(rows)) return(acc)
  R <- length(spec$s); Dl <- length(spec$coords)
  if (!is.matrix(piece$nu) || !identical(dim(piece$nu), c(length(rows), R)) ||
      !identical(dim(piece$rho), c(length(rows), R)) || !identical(dim(piece$P), c(length(rows), Dl)) ||
      !identical(dim(piece$xi), c(length(rows), Dl)) || !identical(dim(piece$blam), c(Dl, R)) ||
      any(!is.finite(piece$nu)) || any(!is.finite(piece$rho)) || any(!is.finite(piece$P)) ||
      any(!is.finite(piece$xi)) || any(!is.finite(piece$blam)))
    cpt_abort(paste0("Process extractor returned malformed or nonfinite rows for cell ", spec$id, "."),
      "cdid_conditional_pretest_integrity_error")
  d <- grid$map[rows]; k <- cmap[rows]
  nr <- piece$nu * piece$rho
  acc$Vd <- cpt_rowsum_into(acc$Vd, nr, d)
  acc$Pd <- cpt_rowsum_into(acc$Pd, piece$P, d)
  acc$Bl <- acc$Bl + piece$blam
  acc$AggX <- cpt_rowsum_into(acc$AggX, piece$xi, k)
  acc$AggKu <- cpt_rowsum_into(acc$AggKu, piece$nu, k)
  acc$sum_nu <- acc$sum_nu + colSums(piece$nu)
  A <- as.vector(nr %*% (spec$s / spec$D))
  nz <- which(A != 0)
  if (length(nz)) {
    acc$A_key <- c(acc$A_key, (as.double(k[nz]) - 1) * grid$n_distinct + d[nz])
    acc$A_val <- c(acc$A_val, A[nz])
  }
  acc$n_rows_seen <- acc$n_rows_seen + length(rows)
  if (!is.null(acc$rows_A)) {
    # Layout rows of one original unit are summed (duplicate row indices).
    acc$rows_A <- as.vector(cpt_rowsum_into(matrix(acc$rows_A), matrix(A), rows))
    acc$rows_nu <- cpt_rowsum_into(acc$rows_nu, piece$nu, rows)
    acc$rows_xi <- cpt_rowsum_into(acc$rows_xi, piece$xi, rows)
  }
  acc
}

# Process values at every distinct grid point (grid chunks of `gchunk`
# columns), the scalar link (h = 1) and the compact multiplier state.
cpt_cell_finalize <- function(acc, spec, grid, N, gchunk) {
  R <- length(spec$s); Dl <- length(spec$coords); Md <- grid$n_distinct
  sD <- spec$s / spec$D
  Vn <- sweep(acc$Vd, 2L, spec$D, "/")
  Rmat <- matrix(0, Md, R); Lmat <- matrix(0, Dl, Md)
  for (cols in cpt_chunks(Md, gchunk)) {
    H <- cpt_indicator(grid$z, grid$z[cols, , drop = FALSE])
    Rc <- crossprod(H, Vn) / N
    Rmat[cols, ] <- Rc
    if (Dl) Lmat[, cols] <- (crossprod(acc$Pd, H) - acc$Bl %*% (t(Rc) * sD)) / N
  }
  J <- spec$sigma * as.vector(Rmat %*% spec$s)
  R1 <- colSums(Vn) / N
  L1 <- if (Dl) as.vector(colSums(acc$Pd) - acc$Bl %*% (R1 * sD)) / N else numeric(0L)
  A_key <- acc$A_key; A_val <- acc$A_val
  if (length(A_key)) {
    s <- rowsum(A_val, A_key, reorder = TRUE)
    A_key <- as.double(rownames(s)); A_val <- as.vector(s)
    keep <- A_val != 0
    A_key <- A_key[keep]; A_val <- A_val[keep]
  }
  A_k <- as.integer(floor((A_key - 1) / Md) + 1)
  A_d <- as.integer(A_key - (A_k - 1) * Md)
  cell <- list(id = spec$id, sigma = spec$sigma, J = J, mult = as.double(grid$mult), M = grid$M,
    zd = grid$z, Rmat = Rmat, Lmat = Lmat, AggX = acc$AggX, AggK = sweep(acc$AggKu, 2L, sD, "*"),
    A_k = A_k, A_d = A_d, A_v = A_val, mu1 = sum(spec$s * R1), L1 = L1, R1 = R1,
    denominator_check = list(supplied = spec$D, accumulated = acc$sum_nu / N),
    coords = spec$coords, n_coordinates = Dl, n_terms = R)
  if (!is.null(acc$rows_A)) cell$rows <- list(A = acc$rows_A,
    kappa = sweep(acc$rows_nu, 2L, sD, "*"), xi = acc$rows_xi, map = grid$map)
  cell
}

# CvM = (N/M) sum_c sum_m J_c(u_m)^2 with M = N grid rows per cell.
cpt_statistic <- function(cells, N) {
  sum(vapply(cells, function(cl) (N / cl$M) * sum(cl$mult * cl$J^2), numeric(1L)))
}

# CvM contribution of one sign batch E (b x C) for every cell. Base functions
# only (the same expression runs serially or on declared worker processes).
cpt_batch_cvmb <- function(E, cells, N, gchunk, echunk) {
  b <- nrow(E)
  out <- numeric(b)
  for (cl in cells) {
    Md <- nrow(cl$zd)
    EX <- E %*% cl$AggX
    EK <- E %*% cl$AggK
    Q <- matrix(0, b, Md)
    ne <- length(cl$A_v)
    if (ne) for (st in seq.int(1L, ne, by = echunk)) {
      ii <- seq.int(st, min(ne, st + echunk - 1L))
      Em <- E[, cl$A_k[ii], drop = FALSE] * rep(cl$A_v[ii], each = b)
      s <- rowsum(t(Em), cl$A_d[ii], reorder = FALSE)
      dd <- as.integer(rownames(s))
      Q[, dd] <- Q[, dd, drop = FALSE] + t(s)
    }
    qd <- which(colSums(Q != 0) > 0)
    for (cs in seq.int(1L, Md, by = gchunk)) {
      cols <- seq.int(cs, min(Md, cs + gchunk - 1L))
      T1 <- if (length(qd)) {
        H <- matrix(TRUE, length(qd), length(cols))
        for (j in seq_len(ncol(cl$zd))) H <- H & outer(cl$zd[qd, j], cl$zd[cols, j], "<=")
        storage.mode(H) <- "double"
        Q[, qd, drop = FALSE] %*% H
      } else matrix(0, b, length(cols))
      Js <- T1 - EK %*% t(cl$Rmat[cols, , drop = FALSE])
      if (ncol(cl$AggX)) Js <- Js + EX %*% cl$Lmat[, cols, drop = FALSE]
      Js <- Js / N
      out <- out + (N / cl$M) * as.vector((Js * Js) %*% cl$mult[cols])
    }
  }
  out
}

# Actual cluster process rows, streamed over grid and event chunks. A finite
# Rademacher sample may annihilate a nonzero process; realized draws therefore
# cannot decide whether the process is globally degenerate. No RNG is used.
cpt_process_nonzero <- function(cells, n_clusters, gchunk, echunk) {
  for (cl in cells) for (cols in cpt_chunks(nrow(cl$zd), gchunk)) {
    Q <- matrix(0, n_clusters, length(cols))
    ne <- length(cl$A_v)
    if (ne) for (st in seq.int(1L, ne, by = echunk)) {
      ii <- seq.int(st, min(ne, st + echunk - 1L))
      H <- cpt_indicator(cl$zd[cl$A_d[ii], , drop = FALSE], cl$zd[cols, , drop = FALSE])
      Q <- cpt_rowsum_into(Q, H * cl$A_v[ii], cl$A_k[ii])
    }
    Q <- Q - cl$AggK %*% t(cl$Rmat[cols, , drop = FALSE])
    if (ncol(cl$AggX)) Q <- Q + cl$AggX %*% cl$Lmat[, cols, drop = FALSE]
    if (any(Q != 0)) return(TRUE)
  }
  FALSE
}

# One common B x C sign stream (caller RNG, batches in draw order, column-major
# within a batch, exactly as the pipeline's multiplier) or caller-fixed signs
# (private arithmetic checks only). Declared workers multiply rounds of drawn
# batches; every draw, statistic and quantile equals the serial result.
cpt_multiplier <- function(cells, N, n_clusters, B, br, gchunk, echunk, fixed_signs = NULL, workers = 1) {
  B <- as.integer(B); br <- as.integer(br)
  if (!is.null(fixed_signs) && (!is.matrix(fixed_signs) || !identical(dim(fixed_signs), c(B, as.integer(n_clusters))) ||
      anyNA(fixed_signs) || any(!fixed_signs %in% c(-1, 1))))
    cpt_abort("fixed_signs must be a biters by declared-clusters matrix of +/-1.", "cdid_conditional_pretest_input_error")
  batches <- cpt_chunks(B, br)
  out <- numeric(B)
  f <- cpt_batch_cvmb
  route <- if (workers > 1) chained_parallel_route() else "serial"
  cl_psock <- NULL
  if (identical(route, "psock")) {
    cl_psock <- parallel::makePSOCKcluster(as.integer(workers))
    on.exit(parallel::stopCluster(cl_psock), add = TRUE)
  }
  pending <- list()
  run_round <- function() {
    smats <- lapply(pending, `[[`, "E")
    res <- if (identical(route, "fork")) {
      parallel::mclapply(smats, f, cells = cells, N = N, gchunk = gchunk, echunk = echunk,
        mc.cores = as.integer(workers), mc.preschedule = TRUE, mc.set.seed = FALSE, mc.silent = TRUE)
    } else {
      fw <- f
      environment(fw) <- baseenv()
      tryCatch(parallel::parLapply(cl_psock, smats, fw, cells = cells, N = N, gchunk = gchunk, echunk = echunk),
        error = function(e) cpt_abort(paste0("A PSOCK conditional-pretest multiplier worker failed (", conditionMessage(e),
          "); no serial fallback or sign regeneration is attempted."), "cdid_conditional_pretest_parallel_error"))
    }
    if (!is.list(res) || length(res) != length(smats))
      cpt_abort("Parallel conditional-pretest workers did not return one result per sign batch; no serial fallback is attempted.",
        "cdid_conditional_pretest_parallel_error")
    for (j in seq_along(res)) {
      r <- res[[j]]
      if (!is.numeric(r) || length(r) != length(pending[[j]]$ii) || any(!is.finite(r)))
        cpt_abort("A conditional-pretest multiplier worker failed or returned a malformed draw vector; no serial fallback is attempted.",
          "cdid_conditional_pretest_parallel_error")
      out[pending[[j]]$ii] <<- r
    }
    pending <<- list()
  }
  for (k in seq_along(batches)) {
    ii <- batches[[k]]
    E <- if (is.null(fixed_signs)) matrix(sample(c(-1.0, 1.0), length(ii) * n_clusters, replace = TRUE), length(ii), n_clusters) else
      unname(fixed_signs[ii, , drop = FALSE]) * 1.0
    if (identical(route, "serial")) {
      out[ii] <- f(E, cells, N, gchunk, echunk)
    } else {
      pending[[length(pending) + 1L]] <- list(ii = ii, E = E)
      if (length(pending) >= workers || k == length(batches)) run_round()
    }
  }
  list(CvMb = out, batch_rows = br, n_batches = length(batches), route = route)
}

# Complete original-unit process rows of one cell (diagnostic, small inputs
# only): N x n_distinct, phi_i(u_m) at every distinct grid point.
cpt_cell_process_rows <- function(cell) {
  if (is.null(cell$rows)) cpt_abort("Process rows were not retained for this cell.", "cdid_conditional_pretest_input_error")
  r <- cell$rows
  H <- cpt_indicator(cell$zd[r$map, , drop = FALSE], cell$zd)
  phi <- H * r$A - r$kappa %*% t(cell$Rmat)
  if (ncol(r$xi)) phi <- phi + r$xi %*% cell$Lmat
  cell$sigma * phi
}

# Scalar link row (h = 1), oriented by the cell's sign.
cpt_cell_scalar_rows <- function(cell) {
  if (is.null(cell$rows)) cpt_abort("Process rows were not retained for this cell.", "cdid_conditional_pretest_input_error")
  r <- cell$rows
  phi <- r$A - as.vector(r$kappa %*% cell$R1)
  if (ncol(r$xi)) phi <- phi + as.vector(r$xi %*% cell$L1)
  cell$sigma * phi
}
