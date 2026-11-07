# Independent dense reference for test-static-nuisance-models.R (test-only).
# Builds full stacked designs [cell dummies, slopes] per declared block and
# uses stats::glm.fit / stats::lm.wfit; does NOT call any package helper.

ref_stack <- function(cells, keys, kind) {
  parts <- lapply(seq_along(keys), function(b) {
    cl <- cells[[keys[b]]]
    rows <- if (kind == "logit") which(cl$mask) else which(cl$mask & cl$observed %in% TRUE)
    X <- if (is.null(cl$design)) matrix(numeric(0), length(rows), 0L) else
      unname(cl$design[rows, , drop = FALSE])
    y <- if (kind == "logit") as.double(cl$response[rows]) else as.double(cl$dy[rows])
    list(unit = rows, cell = rep.int(b, length(rows)), y = y, X = X)
  })
  list(unit = unlist(lapply(parts, `[[`, "unit")), cell = unlist(lapply(parts, `[[`, "cell")),
       y = unlist(lapply(parts, `[[`, "y")), X = do.call(rbind, lapply(parts, `[[`, "X")),
       B = length(keys))
}

ref_full_design <- function(st) cbind(outer(st$cell, seq_len(st$B), "==") * 1, st$X)

ref_block_coef <- function(cells, keys, kind, w_unit) {
  st <- ref_stack(cells, keys, kind)
  Xf <- ref_full_design(st)
  ww <- w_unit[st$unit]
  if (kind == "logit") {
    f <- suppressWarnings(stats::glm.fit(Xf, st$y, weights = ww, family = stats::binomial(),
                                         control = stats::glm.control(epsilon = 1e-15, maxit = 200)))
    stopifnot(f$converged)
    unname(f$coefficients)
  } else {
    unname(stats::lm.wfit(Xf, st$y, ww)$coefficients)
  }
}

ref_fit <- function(cells, partition, kind, w_unit) {
  unlist(lapply(partition, function(keys) ref_block_coef(cells, keys, kind, w_unit)), use.names = FALSE)
}

# Analytic dense IF: N * H_sum^{-1} (unit-summed stacked score), all N units.
ref_if_dense <- function(cells, partition, kind, coef, N) {
  out <- NULL
  off <- 0L
  for (keys in partition) {
    st <- ref_stack(cells, keys, kind)
    Xf <- ref_full_design(st)
    th <- coef[off + seq_len(ncol(Xf))]
    off <- off + ncol(Xf)
    eta <- drop(Xf %*% th)
    if (kind == "logit") {
      p <- stats::plogis(eta)
      H <- crossprod(Xf, (p * (1 - p)) * Xf)
      S <- (st$y - p) * Xf
    } else {
      H <- crossprod(Xf)
      S <- (st$y - eta) * Xf
    }
    Su <- matrix(0, N, ncol(Xf))
    rs <- rowsum(S, st$unit)
    Su[as.integer(rownames(rs)), ] <- rs
    out <- cbind(out, N * t(solve(H, t(Su))))
  }
  out
}

# Empirical original-unit weight perturbation with independent refits.
ref_if_fd <- function(cells, partition, kind, N, units, h = 1e-4) {
  t(vapply(units, function(i) {
    wp <- wm <- rep.int(1, N)
    wp[i] <- 1 + h
    wm[i] <- 1 - h
    N * (ref_fit(cells, partition, kind, wp) - ref_fit(cells, partition, kind, wm)) / (2 * h)
  }, numeric(length(ref_fit(cells, partition, kind, rep.int(1, N))))))
}

ref_expand <- function(pred, P) {
  D <- matrix(0, nrow(pred$deriv), P)
  D[, pred$param_index] <- pred$deriv
  D
}
