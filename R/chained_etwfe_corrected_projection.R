# The catalog retains every pair with observed change support for each cohort.
# Its reconstruction is evaluated on all units of that cohort. A missing pair
# is an absent equation and is reported in the fitted support record.
etwfe_corrected_catalog <- function(index, cohort, periods, entries) {
  n <- length(cohort)
  groups <- sort(unique(cohort))
  cost <- 4 * as.double(n) * periods + 4 * as.double(length(groups)) * periods^2
  if (cost > entries) etwfe_tm_workspace_error()
  Y <- matrix(NA_real_, n, periods)
  Y[cbind(index$unit, index$period)] <- index$y
  units <- split(seq_len(n), factor(cohort, levels = groups))
  counts <- lapply(units, function(rows) crossprod(!is.na(Y[rows, , drop = FALSE])))
  cells <- do.call(rbind, lapply(seq_along(groups), function(k) {
    count <- counts[[k]]
    pairs <- which(upper.tri(count) & count > 0L, arr.ind = TRUE)
    data.frame(cohort = rep.int(groups[k], nrow(pairs)), s = pairs[, 1L], t = pairs[, 2L],
               observed = count[pairs])
  }))
  cells$group <- match(cells$cohort, groups)
  observed <- vapply(counts, function(count) diag(count) > 0L, logical(periods))
  list(Y = Y, cells = cells, units = units, groups = groups, entries = cost, observed = observed,
       omitted_pairs = length(groups) * choose(periods, 2) - nrow(cells))
}

etwfe_corrected_features <- function(Z) {
  centered <- sweep(Z, 2L, colMeans(Z), "-")
  scale <- vapply(seq_len(ncol(centered)), function(k) max(abs(centered[, k])), 0)
  scale[scale == 0] <- 1
  if (ncol(centered)) centered <- sweep(centered, 2L, scale, "/")
  cbind(1, centered)
}

etwfe_corrected_projection_geometry <- function(catalog, slot, use, requested_slots, entries) {
  cells <- catalog$cells
  A <- max(slot)
  if (3 * as.double(A)^2 + catalog$entries > entries) etwfe_tm_workspace_error()
  left <- slot[cbind(cells$s[use], cells$group[use])]
  right <- slot[cbind(cells$t[use], cells$group[use])]
  link <- matrix(FALSE, A, A)
  link[cbind(left, right)] <- TRUE
  link[cbind(right, left)] <- TRUE
  connected <- etwfe_core_components(link)
  estimable <- connected[requested_slots[, 1L]] == connected[requested_slots[, 2L]]
  used_components <- unique(connected[requested_slots[estimable, 1L]])
  active <- which(connected[left] %in% used_components & left != right)
  list(use = use[active], left = left[active], right = right[active],
       keep_slots = which(duplicated(connected) & connected %in% used_components),
       estimable = estimable, n_slots = A)
}

# Within a cohort the full-frame feature Gram is identical on every admitted
# edge. A symmetric linked-slot Gram and its RHS are sufficient.
etwfe_corrected_projection <- function(catalog, features, geometry, obs,
                                      compute_scores, chunk_rows, entries, extra = 0) {
  n <- nrow(features)
  q <- ncol(features)
  cells <- catalog$cells
  A <- geometry$n_slots
  estimable <- geometry$estimable
  use <- geometry$use
  left <- geometry$left
  right <- geometry$right
  keep_slots <- geometry$keep_slots
  p <- as.double(A) * q
  if (p > .Machine$integer.max) etwfe_tm_workspace_error()
  coords <- matrix(seq_len(p), A, q)
  keep <- as.vector(coords[keep_slots, , drop = FALSE])
  beta <- matrix(0, A, q)
  if (!length(keep)) {
    if (extra + catalog$entries + p > entries) etwfe_tm_workspace_error()
    return(list(beta = beta, estimable = estimable, live_entries = length(beta),
      score = if (compute_scores) function(rows, nuisance) matrix(0, length(rows), p)))
  }

  nu <- obs$n_parameters
  per_unit <- 4 * p + 2 * nu + nrow(cells) + 3 * q + 8
  retained_obs <- if (is.null(obs$diagnostics$budget$retained_entries)) 0 else obs$diagnostics$budget$retained_entries
  retained <- as.double(retained_obs) + catalog$entries + 6 * as.double(p)^2 + 3 * as.double(A)^2 +
    compute_scores * 3 * as.double(p) * nu
  if (retained + extra + per_unit > entries) etwfe_tm_workspace_error()
  block_size <- min(chunk_rows, floor((entries - retained - extra) / per_unit))
  blocks <- function(ids) split(ids, (seq_along(ids) - 1L) %/% block_size)

  gram <- matrix(0, p, p)
  rhs <- numeric(p)
  loading <- if (compute_scores) matrix(0, p, nu)
  feature_grams <- lapply(catalog$units, function(rows) crossprod(features[rows, , drop = FALSE]))
  for (i in seq_along(use)) {
    cell <- use[i]
    group <- cells$group[cell]
    s <- coords[left[i], ]
    t <- coords[right[i], ]
    G <- feature_grams[[group]]
    gram[t, t] <- gram[t, t, drop = FALSE] + G
    gram[s, s] <- gram[s, s, drop = FALSE] + G
    gram[t, s] <- gram[t, s, drop = FALSE] - G
    gram[s, t] <- gram[s, t, drop = FALSE] - G
    for (rows in blocks(catalog$units[[group]])) {
      out <- obs$read(cell, rows, gradient = compute_scores)
      w <- features[rows, , drop = FALSE]
      a <- drop(crossprod(w, out$value))
      rhs[t] <- rhs[t] + a
      rhs[s] <- rhs[s] - a
      if (compute_scores && length(out$coordinates)) {
        derivative <- crossprod(w, out$gradient)
        loading[t, out$coordinates] <- loading[t, out$coordinates, drop = FALSE] + derivative
        loading[s, out$coordinates] <- loading[s, out$coordinates, drop = FALSE] - derivative
      }
    }
  }
  G <- gram[keep, keep, drop = FALSE]
  D <- sqrt(diag(G))
  R <- if (all(is.finite(D) & D > 0)) {
    tryCatch(chol(t(t(G / D) / D)), error = function(e) NULL)
  }
  rc <- if (is.null(R)) NA_real_ else rcond(R, triangular = TRUE)
  if (!is.finite(rc) || rc^2 <= .Machine$double.eps) {
    cdid_abort("The observation-corrected cohort-time projection is rank deficient after fixing calendar gauges.",
               "cdid_error_etwfe_numeric")
  }
  solve_scaled <- function(z) cp04_solve(R, z / D) / D
  beta[keep] <- solve_scaled(rhs[keep])
  if (!compute_scores) return(list(beta = beta, estimable = estimable, score = NULL,
                                   live_entries = length(beta)))

  # Fold the nuisance loading into the same slot factor used for the point.
  # Each row's outcome part is evaluated once on its active cohort edges.
  adjusted <- if (nu) solve_scaled(loading[keep, , drop = FALSE]) else matrix(0, length(keep), 0L)
  group_of_row <- integer(n)
  group_of_row[unlist(catalog$units, use.names = FALSE)] <- rep(seq_along(catalog$units), lengths(catalog$units))
  active_in_group <- split(seq_along(use), factor(cells$group[use], levels = seq_along(catalog$groups)))
  score <- function(rows, nuisance) {
    values <- matrix(0, length(rows), p)
    for (group in unique(group_of_row[rows])) {
      selected <- which(group_of_row[rows] == group)
      ix <- rows[selected]
      w <- features[ix, , drop = FALSE]
      for (i in active_in_group[[group]]) {
        s <- left[i]
        t <- right[i]
        dy <- obs$read(use[i], ix)$value
        e <- n * (dy - drop(w %*% (beta[t, ] - beta[s, ])))
        dx <- w * e
        values[selected, coords[t, ]] <- values[selected, coords[t, ], drop = FALSE] + dx
        values[selected, coords[s, ]] <- values[selected, coords[s, ], drop = FALSE] - dx
      }
    }
    scores <- matrix(0, length(rows), p)
    scores[, keep] <- t(solve_scaled(t(values[, keep, drop = FALSE]))) + tcrossprod(nuisance, adjusted)
    scores
  }
  # The projection closure only needs the factor, coefficients and loadings.
  rm(gram, rhs, G, loading, feature_grams)
  list(beta = beta, estimable = estimable, score = score, block_rows = block_size,
       live_entries = length(R) + length(adjusted) + 3 * p + n + 6 * length(use))
}

etwfe_corrected_regression_geometry <- function(catalog, control_group, entries) {
  groups <- catalog$groups
  periods <- ncol(catalog$Y)
  dates <- seq_len(periods)
  treated <- which(groups > 0L)
  free <- if (identical(control_group, "nevertreated")) outer(dates, groups[treated] - 1L, "!=") else
    outer(dates, groups[treated], ">=")
  slot <- matrix(dates, periods, length(groups))
  slot[, treated] <- ifelse(free, outer(dates, (treated - 1L) * periods, "+"), dates)
  target <- which(free, arr.ind = TRUE)
  target_date <- target[, 1L]
  target_group <- treated[target[, 2L]]
  nodes <- cbind(slot[cbind(target_date, target_group)], target_date)
  geometry <- etwfe_corrected_projection_geometry(
    catalog, slot, seq_len(nrow(catalog$cells)), nodes, entries
  )
  geometry$nodes <- nodes
  geometry$target_group <- target_group
  geometry$target_date <- target_date
  geometry
}

etwfe_corrected_regression <- function(catalog, features, control_group, obs,
                                       compute_scores, chunk_rows, entries,
                                       geometry = etwfe_corrected_regression_geometry(catalog, control_group, entries),
                                       report = NULL, Z = NULL) {
  n <- nrow(catalog$Y)
  groups <- catalog$groups
  target_group <- geometry$target_group
  target_date <- geometry$target_date
  nodes <- geometry$nodes
  J <- length(target_date)
  reserve <- if (compute_scores) as.double(n) * J + as.double(geometry$n_slots) * ncol(features) * J else 0
  fit <- etwfe_corrected_projection(catalog, features, geometry, obs,
                                    compute_scores, chunk_rows, entries, extra = reserve)
  target_x <- t(matrix(vapply(catalog$units, function(i) colMeans(features[i, , drop = FALSE]),
                              numeric(ncol(features))), nrow = ncol(features)))
  ids <- which(fit$estimable)
  diff <- fit$beta[nodes[ids, 1L], , drop = FALSE] - fit$beta[nodes[ids, 2L], , drop = FALSE]
  att <- rep(NA_real_, J)
  att[ids] <- rowSums(diff * target_x[target_group[ids], , drop = FALSE])

  observed <- catalog$observed[cbind(target_date, target_group)]
  status <- ifelse(fit$estimable, "identified", ifelse(observed, "unidentified", "no_support"))
  targets <- data.frame(cohort = groups[target_group], period = target_date,
                        is_pre = target_date < groups[target_group], att = att, status = status)
  if (reserve > entries) etwfe_tm_workspace_error()
  influence <- if (compute_scores) matrix(NA_real_, n, J)
  if (compute_scores && length(ids)) {
    p <- length(fit$beta)
    coord <- matrix(seq_len(p), nrow(fit$beta))
    mix <- target_x[target_group[ids], , drop = FALSE]
    C <- matrix(0, p, length(ids))
    for (a in seq_len(ncol(features))) {
      C[cbind(coord[nodes[ids, 1L], a], seq_along(ids))] <- mix[, a]
      C[cbind(coord[nodes[ids, 2L], a], seq_along(ids))] <- -mix[, a]
    }
    for (first in seq.int(1L, n, by = fit$block_rows)) {
      rows <- seq.int(first, min(n, first + fit$block_rows - 1L))
      V <- fit$score(rows, obs$score_rows(rows)) %*% C
      for (g in unique(target_group[ids])) {
        u <- which(rows %in% catalog$units[[g]])
        if (!length(u)) next
        v <- which(target_group[ids] == g)
        W <- sweep(features[rows[u], , drop = FALSE], 2L, target_x[g, ], "-")
        V[u, v] <- V[u, v, drop = FALSE] + n / length(catalog$units[[g]]) *
          tcrossprod(W, diff[v, , drop = FALSE])
      }
      influence[rows, ids] <- V
    }
  }
  result <- list(targets = targets, inffunc = influence,
    diagnostics = list(n_moments = nrow(catalog$cells), omitted_pairs = catalog$omitted_pairs))
  if (!is.null(report)) {
    cells <- catalog$cells[geometry$use, , drop = FALSE]
    comparisons <- etwfe_nuisance_comparisons(cells$cohort, cells$s, cells$t,
      rep(NA_integer_, nrow(cells)), rep("calendar", nrow(cells)), report)
    result$nuisance <- list(outcome = etwfe_nuisance_projection(catalog, geometry, fit, obs,
      Z, report, entries - reserve, "corrected_levels", comparisons))
  }
  result
}
