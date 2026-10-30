# D52 independent dense-small pure actual-change imputation reference.
# No cdid/did imports. Fixed positive case weights; original F identities retained.
# Dense implicit derivatives are N*d/dc_i at c=1 for w_i -> c_i*w_i.
# Only unchanged identity-GMM/reporting maps are supplied by the public helpers.
# This is a supported-design finite oracle, not a package or rank-policy emulator.

imputation_public_helpers <- function(helper_dir = getOption(
    "imputation.public.helper_dir", file.path("input", "reference"))) {
  env <- environment(imputation_public_helpers)
  needed <- c("cp06_pooled_public_catalog", "sp_ref_public_eval", "cp07_public_dr_fit_logit")
  files <- c("cp06_pooled_public_or_reference.R", "selection_public_reference.R",
             "cp07_public_dr_reference.R")
  for (i in seq_along(files)) {
    if (!exists(needed[i], envir = env, mode = "function", inherits = TRUE)) {
      path <- file.path(helper_dir, files[i])
      if (!file.exists(path)) stop("Missing supplied independent helper: ", path)
      sys.source(path, envir = env)
    }
  }
  invisible(NULL)
}

imputation_public_raw <- function(raw) {
  # Public 0 and fixture NT are both labels for the internal Inf cohort.
  g <- as.character(raw$cohort)
  g[g %in% c("NT", "0")] <- "Inf"
  raw$cohort <- as.numeric(g)
  stopifnot(!anyNA(raw$cohort), !anyDuplicated(raw$id),
            all(is.finite(raw$w) & raw$w > 0))
  raw
}

imputation_public_fixture <- function(raw, dates = 1:4) {
  raw <- imputation_public_raw(raw)
  N <- nrow(raw)
  Y <- as.matrix(raw[paste0("Y", seq_along(dates))])
  ids <- paste0("unit-", raw$id)
  g <- ifelse(is.infinite(raw$cohort), 0, raw$cohort)
  panel <- data.frame(id = rep(ids, each = length(dates)), time = rep(dates, N),
    y = as.vector(t(Y)), g = rep(g, each = length(dates)), stringsAsFactors = FALSE)
  zero <- which(rowSums(!is.na(Y)) == 0L)
  if (length(zero) > 1L) panel <- panel[!panel$id %in% ids[zero[-1L]], , drop = FALSE]
  frame <- data.frame(id = ids, g = g, x = raw$x, w = raw$w,
    cl = paste0("cl-", (raw$id - 1L) %/% 3L), stringsAsFactors = FALSE)
  frame[["(Intercept)"]] <- raw$x + 1.5
  list(raw = raw, frame = frame, data = panel,
    args = list(yname = "y", tname = "time", idname = "id", gname = "g",
      data = panel, frame = frame, xformla = ~x, weightsname = "w",
      selection = "conditional_mtar", selection_est_method = "imputation",
      imputation_pooling = "separate", est_method = "dr", control_group = "nevertreated",
      gmm_weighting = "identity", bstrap = FALSE, cband = FALSE))
}

imputation_public_blueprint <- function(raw, imputation_pooling = "separate",
    intercept = TRUE, X_slope = cbind(x = raw$x),
    V_treat = if (intercept) cbind(`(Intercept)` = rep(1, nrow(raw)), X_slope) else X_slope,
    dates = 1:4, control_group = "nevertreated", anticipation = 0L,
    pt_assumption = "post", pt_start = NULL, base_period = "varying") {
  imputation_public_helpers()
  raw <- imputation_public_raw(raw)
  imputation_pooling <- match.arg(imputation_pooling,
    c("separate", "within_cohort", "across_cohorts"))
  stopifnot(length(dates) >= 2L, identical(as.numeric(dates), sort(unique(as.numeric(dates)))))
  N <- nrow(raw)
  Z <- if (is.null(X_slope)) matrix(0, N, 0L) else as.matrix(X_slope)
  V <- if (is.null(V_treat)) matrix(0, N, 0L) else as.matrix(V_treat)
  stopifnot(nrow(Z) == N, nrow(V) == N, all(is.finite(Z)), all(is.finite(V)))
  if (!intercept && !ncol(Z)) stop("Empty outcome design has no supplied public convention")
  Y <- as.matrix(raw[paste0("Y", seq_along(dates))])
  R <- !is.na(Y)
  if (all(paste0("R", seq_along(dates)) %in% names(raw)))
    stopifnot(all(R == as.matrix(raw[paste0("R", seq_along(dates))])))
  A <- rowSums(R) >= 2L
  if ("A" %in% names(raw)) stopifnot(all(A == raw$A))
  if ("B" %in% names(raw)) stopifnot(all((rowSums(R) >= 1L) == raw$B))
  ci <- cp06_pooled_public_catalog(raw, dates, control_group, anticipation,
                                 pt_assumption, pt_start, base_period)
  links <- ci$links_df
  if (!nrow(links)) stop("No supported public comparison catalog")
  cohorts <- ci$cohorts_all
  if (!all(raw$cohort %in% cohorts)) stop("Full-F cohort absent from supplied public catalog")
  controls <- lapply(links$controls, function(s) as.numeric(strsplit(s, "|", fixed = TRUE)[[1L]]))
  cells <- do.call(rbind, lapply(seq_len(nrow(links)), function(l)
    data.frame(k = c(links$g[l], controls[[l]]), is = links$is[l], it = links$it[l])))
  cells <- unique(cells)
  cells <- cells[order(match(cells$k, cohorts), cells$is, cells$it), , drop = FALSE]
  rownames(cells) <- NULL
  cells$L <- cells$it - cells$is
  key <- function(k, is, it) paste(k, is, it, sep = ":")
  cell_id <- function(k, is, it) match(key(k, is, it), key(cells$k, cells$is, cells$it))
  groups <- switch(imputation_pooling,
    separate = as.list(seq_len(nrow(cells))),
    within_cohort = lapply(unique(cells$k), function(k) which(cells$k == k)),
    across_cohorts = list(seq_len(nrow(cells))))
  designs <- vector("list", nrow(cells))
  train <- vector("list", length(groups))
  for (r in seq_along(groups)) {
    cs <- groups[[r]]
    rows <- list(); yy <- list(); xx <- list()
    for (b in seq_along(cs)) {
      c <- cs[b]
      E <- matrix(0, N, if (intercept) length(cs) else 0L)
      if (intercept) E[, b] <- 1
      D <- cbind(E, cells$L[c] * Z)
      colnames(D) <- c(if (intercept) paste0("alpha.cell", cs) else character(),
                      if (ncol(Z)) paste0("beta.", seq_len(ncol(Z))) else character())
      designs[[c]] <- D
      ii <- which(raw$cohort == cells$k[c] & R[, cells$is[c]] & R[, cells$it[c]])
      if (!length(ii)) stop("No observed training pair for actual cell ", c)
      rows[[b]] <- ii
      yy[[b]] <- Y[ii, cells$it[c]] - Y[ii, cells$is[c]]
      xx[[b]] <- D[ii, , drop = FALSE]
    }
    train[[r]] <- list(rows = unlist(rows), y = unlist(yy), X = do.call(rbind, xx))
  }
  # Same pair/block order and deduplication as the supplied CP07 public helper.
  pairs <- unique(links[c("s", "t", "is", "it")])
  pairs <- pairs[order(pairs$is, pairs$it), , drop = FALSE]
  rownames(pairs) <- NULL
  pairs$pair_id <- seq_len(nrow(pairs)); pairs$L <- pairs$it - pairs$is
  link_pair <- match(paste(links$s, links$t), paste(pairs$s, pairs$t))
  blocks <- unique(data.frame(pair_id = link_pair, controls = links$controls))
  pool_sort <- vapply(blocks$controls, function(s) paste(sprintf("%012d",
    match(as.numeric(strsplit(s, "|", fixed = TRUE)[[1L]]), cohorts)), collapse = "|"), character(1))
  blocks <- blocks[order(sprintf("%012d|%s", blocks$pair_id, pool_sort), method = "radix"), , drop = FALSE]
  rownames(blocks) <- NULL
  link_block <- match(paste(link_pair, links$controls), paste(blocks$pair_id, blocks$controls))
  or_design <- lapply(seq_len(nrow(blocks)), function(b) {
    E <- matrix(0, N, if (intercept) nrow(blocks) else 0L)
    if (intercept) E[, b] <- 1
    cbind(E, pairs$L[blocks$pair_id[b]] * Z)
  })
  mapping_args <- list(anticipation = anticipation, pt_assumption = pt_assumption,
                       pt_start = pt_start, base_period = base_period)
  map <- function(x) do.call(cp06_pooled_public_cells_from_mu,
    c(list(cat_info = ci, mu_all = x), mapping_args))
  zero <- map(rep(0, nrow(links)))
  loading <- vapply(seq_len(nrow(links)), function(l) {
    e <- numeric(nrow(links)); e[l] <- 1; map(e)$att
  }, numeric(nrow(zero)))
  list(raw = raw, N = N, A = A, R = R, Y = Y, Z = Z, V = V, cells = cells,
    groups = groups, designs = designs, train = train, cohorts = cohorts, controls = controls,
    links = links, blocks = blocks, pairs = pairs, link_pair = link_pair,
    link_block = link_block, or_design = or_design, cell_id = cell_id, cat_info = ci,
    mapping_args = mapping_args, loading = loading, pooling = imputation_pooling)
}

# Verbatim aggregation arithmetic from CP07, factored into a callable map.
imputation_public_aggregate <- function(cells_out, tau) {
  post <- which(cells_out$time >= cells_out$g)
  agg_overall <- c(simple = NA_real_, group = NA_real_, calendar = NA_real_, dynamic = NA_real_)
  if (length(post) && all(!is.na(cells_out$att[post]))) {
    g_post <- as.character(cells_out$g[post])
    agg_overall["simple"] <- sum(tau[g_post] * cells_out$att[post]) / sum(tau[g_post])
    ug <- sort(unique(cells_out$g[post]))
    est_g <- vapply(ug, function(gv) mean(cells_out$att[post[cells_out$g[post] == gv]]), numeric(1))
    agg_overall["group"] <- sum(tau[as.character(ug)] * est_g) / sum(tau[as.character(ug)])
    ut <- sort(unique(cells_out$time[post]))
    est_t <- vapply(ut, function(tv) {
      rt <- post[cells_out$time[post] == tv]
      gt <- as.character(cells_out$g[rt])
      sum(tau[gt] * cells_out$att[rt]) / sum(tau[gt])
    }, numeric(1))
    agg_overall["calendar"] <- mean(est_t)
    ue <- sort(unique(cells_out$time[post] - cells_out$g[post]))
    est_e <- vapply(ue, function(ev) {
      re <- post[(cells_out$time[post] - cells_out$g[post]) == ev]
      ge <- as.character(cells_out$g[re])
      sum(tau[ge] * cells_out$att[re]) / sum(tau[ge])
    }, numeric(1))
    agg_overall["dynamic"] <- mean(est_e)
  }
  agg_overall
}

# Dense differential of precisely the four unchanged rational reporting maps.
# This differentiates their ATT/share inputs, not an estimated GMM weight.
imputation_public_aggregate_rows <- function(cells, tau, phi, dtau) {
  N <- nrow(phi)
  out <- matrix(NA_real_, N, 4L,
                dimnames = list(NULL, c("simple", "group", "calendar", "dynamic")))
  post <- which(cells$time >= cells$g)
  if (!length(post) || anyNA(cells$att[post])) return(out)
  weighted <- function(att, rows, cohorts) {
    idx <- match(as.character(cohorts), names(tau))
    weights <- tau[idx]
    value <- sum(weights * att) / sum(weights)
    drop(rows %*% weights + t(dtau[idx, , drop = FALSE]) %*% (att - value)) / sum(weights)
  }
  out[, "simple"] <- weighted(cells$att[post], phi[, post, drop = FALSE], cells$g[post])
  ug <- sort(unique(cells$g[post]))
  att_g <- vapply(ug, function(g) mean(cells$att[post[cells$g[post] == g]]), numeric(1))
  rows_g <- vapply(ug, function(g) rowMeans(phi[, post[cells$g[post] == g], drop = FALSE]), numeric(N))
  out[, "group"] <- weighted(att_g, rows_g, ug)
  ut <- sort(unique(cells$time[post]))
  rows_t <- vapply(ut, function(t) {
    rr <- post[cells$time[post] == t]
    weighted(cells$att[rr], phi[, rr, drop = FALSE], cells$g[rr])
  }, numeric(N))
  out[, "calendar"] <- rowMeans(rows_t)
  ue <- sort(unique(cells$time[post] - cells$g[post]))
  rows_e <- vapply(ue, function(e) {
    rr <- post[cells$time[post] - cells$g[post] == e]
    weighted(cells$att[rr], phi[, rr, drop = FALSE], cells$g[rr])
  }, numeric(N))
  out[, "dynamic"] <- rowMeans(rows_e)
  out
}

imputation_public_dense_eval <- function(bp, w = bp$raw$w,
    est_method = c("ipw", "reg", "dr"), derivatives = TRUE) {
  methods <- match.arg(est_method, c("ipw", "reg", "dr"), several.ok = TRUE)
  N <- bp$N; g <- bp$raw$cohort; K <- bp$cohorts; C <- nrow(bp$cells)
  stopifnot(length(w) == N, all(is.finite(w) & w > 0))
  use_or <- any(methods %in% c("reg", "dr"))
  use_tr <- any(methods %in% c("ipw", "dr"))
  pk <- vapply(K, function(k) sum(w[g == k]), numeric(1))
  qk <- vapply(K, function(k) sum(w[g == k & bp$A]), numeric(1))
  if (any(pk <= 0 | qk <= 0)) stop("Nonpositive cohort/repeat mass")
  a <- qk / pk; tau <- pk / sum(pk); eta <- qk / sum(qk)
  names(a) <- names(tau) <- names(eta) <- as.character(K)
  da <- dtau <- deta <- matrix(0, length(K), N)
  rownames(da) <- rownames(dtau) <- rownames(deta) <- as.character(K)
  if (derivatives) for (k in seq_along(K)) {
    da[k, ] <- N * w * (g == K[k]) * (bp$A - a[k]) / pk[k]
    dtau[k, ] <- N * w * ((g == K[k]) - tau[k]) / sum(pk)
    deta[k, ] <- N * w * bp$A * ((g == K[k]) - eta[k]) / sum(qk)
  }
  snapshot <- numeric(); rows <- matrix(0, 0L, N)
  add <- function(prefix, value, phi) {
    if (!length(value)) return(invisible(NULL))
    nm <- if (length(value) == 1L) prefix else paste0(prefix, ".", seq_along(value))
    snapshot <<- c(snapshot, setNames(as.numeric(value), nm))
    if (derivatives) rows <<- rbind(rows, matrix(phi, nrow = length(value), ncol = N,
                                               dimnames = list(nm, NULL)))
  }
  add("a", a, da); add("tau", tau, dtau); add("eta", eta, deta)
  roots <- numeric(); ranks <- list()
  j <- vector("list", C); dj <- vector("list", C); actual <- list()
  for (r in seq_along(bp$groups)) {
    tr <- bp$train[[r]]; X <- tr$X; ww <- w[tr$rows]
    fit <- stats::lm.wfit(X, tr$y, ww, tol = 1e-12)
    if (fit$rank != ncol(X) || any(!is.finite(fit$coefficients)))
      stop("Actual-imputation design rank failure in group ", r)
    coef <- unname(fit$coefficients)
    G <- crossprod(X, ww * X)
    roots <- c(roots, max(abs(crossprod(X, ww * (tr$y - drop(X %*% coef))))))
    phi <- matrix(0, ncol(X), N)
    if (derivatives) {
      score <- rowsum(X * (ww * (tr$y - drop(X %*% coef))), tr$rows, reorder = FALSE)
      rhs <- matrix(0, ncol(X), N)
      rhs[, as.integer(rownames(score))] <- N * t(score)
      phi <- qr.solve(G, rhs, tol = 1e-12)
    }
    add(paste0("imputation.group", r), coef, phi)
    actual[[r]] <- list(coef = coef, phi = phi, rank = fit$rank, design = colnames(X))
    for (c in bp$groups[[r]]) {
      j[[c]] <- drop(bp$designs[[c]] %*% coef)
      dj[[c]] <- if (derivatives) bp$designs[[c]] %*% phi else NULL
    }
  }
  # Synthetic F mass and its complete case-weight derivative, one cohort at a time.
  mass <- lapply(seq_along(K), function(k) w * (g == K[k]) * a[k])
  dmass <- lapply(seq_along(K), function(k) if (derivatives)
    diag(N * mass[[k]], N) + outer(w * (g == K[k]), da[k, ]) else NULL)
  theta <- numeric(); dtheta <- matrix(0, 0L, N)
  if (use_or) {
    d <- ncol(bp$or_design[[1L]])
    G <- matrix(0, d, d); rhs <- numeric(d)
    for (b in seq_len(nrow(bp$blocks))) {
      D <- bp$or_design[[b]]; p <- bp$pairs[bp$blocks$pair_id[b], ]
      for (k in as.numeric(strsplit(bp$blocks$controls[b], "|", fixed = TRUE)[[1L]])) {
        ki <- match(k, K); c <- bp$cell_id(k, p$is, p$it)
        G <- G + crossprod(D, mass[[ki]] * D)
        rhs <- rhs + drop(crossprod(D, mass[[ki]] * j[[c]]))
      }
    }
    theta <- drop(qr.solve(G, rhs, tol = 1e-12))
    roots <- c(roots, max(abs(G %*% theta - rhs)))
    dtheta <- matrix(0, d, N)
    if (derivatives) {
      forcing <- matrix(0, d, N)
      for (b in seq_len(nrow(bp$blocks))) {
        D <- bp$or_design[[b]]; p <- bp$pairs[bp$blocks$pair_id[b], ]
        for (k in as.numeric(strsplit(bp$blocks$controls[b], "|", fixed = TRUE)[[1L]])) {
          ki <- match(k, K); c <- bp$cell_id(k, p$is, p$it)
          res <- j[[c]] - drop(D %*% theta)
          forcing <- forcing + crossprod(D, res * dmass[[ki]] + mass[[ki]] * dj[[c]])
        }
      }
      dtheta <- qr.solve(G, forcing, tol = 1e-12)
    }
    add("OR.theta", theta, dtheta)
  }
  # H_N raw moments: dq includes direct empirical mass and every a path;
  # dh includes fitted j, treatment odds and/or shared-OR paths as appropriate.
  moment <- function(ki, h, dh = NULL) {
    value <- sum(mass[[ki]] * h) / N
    phi <- numeric(N)
    if (derivatives) {
      phi <- colSums(h * dmass[[ki]]) / N
      if (!is.null(dh)) phi <- phi + drop(crossprod(mass[[ki]], dh)) / N
    }
    list(value = value, phi = phi)
  }
  combine <- function(xs) list(value = sum(vapply(xs, `[[`, numeric(1), "value")),
    phi = Reduce(`+`, lapply(xs, `[[`, "phi")))
  ratio <- function(num, den) {
    if (!is.finite(den$value) || den$value <= 0) stop("Nonpositive explicit ratio normalizer")
    value <- num$value / den$value
    list(value = value, phi = (num$phi - value * den$phi) / den$value)
  }
  L <- nrow(bp$links)
  mu <- matrix(0, L, length(methods), dimnames = list(NULL, methods))
  phi_mu <- setNames(lapply(methods, function(m) matrix(0, N, L)), methods)
  treatment <- list(); denominators <- matrix(NA_real_, L, 2L,
                                            dimnames = list(NULL, c("treated", "control")))
  for (l in seq_len(L)) {
    link <- bp$links[l, ]; ks <- c(link$g, bp$controls[[l]])
    kis <- match(ks, K); cs <- bp$cell_id(ks, link$is, link$it)
    plus <- moment(kis[1L], rep(1, N))
    numT <- moment(kis[1L], j[[cs[1L]]], dj[[cs[1L]]])
    T <- ratio(numT, plus)
    denominators[l, 1L] <- plus$value
    values <- list(denplus = plus, numT = numT, T = T)
    if (use_or) {
      D <- bp$or_design[[bp$link_block[l]]]
      pred <- drop(D %*% theta)
      dpred <- if (derivatives) D %*% dtheta else NULL
      numBg <- moment(kis[1L], pred, dpred)
      Bg <- ratio(numBg, plus)
      values$numBg <- numBg; values$Bg <- Bg
    }
    if (use_tr) {
      q <- Reduce(`+`, mass[kis]); V <- bp$V; treated <- as.numeric(g == link$g)
      ii <- which(q > 0)
      gamma <- cp07_public_dr_fit_logit(V[ii, , drop = FALSE], treated[ii], q[ii])
      odds <- if (ncol(V)) exp(drop(V %*% gamma)) else rep(1, N)
      if (any(!is.finite(odds))) stop("Nonfinite treatment odds")
      pr <- odds / (1 + odds)
      roots <- c(roots, if (ncol(V)) max(abs(crossprod(V, q * (treated - pr)))) else 0)
      dgamma <- matrix(0, ncol(V), N)
      dodds <- NULL
      if (derivatives && ncol(V)) {
        info <- crossprod(V, (q * pr * (1 - pr)) * V)
        dgamma <- qr.solve(info, crossprod(V, (treated - pr) * Reduce(`+`, dmass[kis])), tol = 1e-12)
        dodds <- odds * (V %*% dgamma)
      }
      add(paste0("treatment.link", l), gamma, dgamma)
      treatment[[l]] <- gamma
      minus <- combine(lapply(kis[-1L], function(ki) moment(ki, odds, dodds)))
      numC <- combine(lapply(seq_along(kis)[-1L], function(z) {
        c <- cs[z]
        dh <- if (derivatives) odds * dj[[c]] + (if (is.null(dodds)) 0 else j[[c]] * dodds) else NULL
        moment(kis[z], odds * j[[c]], dh)
      }))
      Cmean <- ratio(numC, minus)
      values$denminus <- minus; values$numC <- numC; values$C <- Cmean
      denominators[l, 2L] <- minus$value
      if (use_or) {
        dh <- if (derivatives) odds * dpred + (if (is.null(dodds)) 0 else pred * dodds) else NULL
        numBo <- combine(lapply(kis[-1L], function(ki) moment(ki, odds * pred, dh)))
        Bo <- ratio(numBo, minus)
        values$numBo <- numBo; values$Bo <- Bo
      }
    }
    for (nm in names(values)) add(paste0(nm, ".link", l), values[[nm]]$value, values[[nm]]$phi)
    for (m in methods) {
      vv <- switch(m, ipw = T$value - Cmean$value, reg = T$value - Bg$value,
                   dr = T$value - Bg$value - Cmean$value + Bo$value)
      pp <- switch(m, ipw = T$phi - Cmean$phi, reg = T$phi - Bg$phi,
                   dr = T$phi - Bg$phi - Cmean$phi + Bo$phi)
      mu[l, m] <- vv; phi_mu[[m]][, l] <- pp
      add(paste0(m, ".mu.link", l), vv, pp)
    }
  }
  # Map links with the unchanged identity-GMM helper, including auxiliary rows.
  # Basis evaluation applies its fixed linear map; no new GMM derivative is fit.
  map <- function(x) do.call(cp06_pooled_public_cells_from_mu,
                             c(list(cat_info = bp$cat_info, mu_all = x), bp$mapping_args))
  loading <- bp$loading
  results <- list()
  for (m in methods) {
    cells <- map(mu[, m])
    phi <- phi_mu[[m]] %*% t(loading)
    agg <- imputation_public_aggregate(cells, tau)
    pagg <- matrix(NA_real_, N, length(agg), dimnames = list(NULL, names(agg)))
    if (derivatives) pagg <- imputation_public_aggregate_rows(cells, tau, phi, dtau)
    for (r in seq_along(cells$att)) add(paste0(m, ".att", r), cells$att[r], phi[, r])
    for (r in seq_along(agg)) add(paste0(m, ".agg.", names(agg)[r]), agg[r], pagg[, r])
    results[[m]] <- list(mu = mu[, m], phi_mu = phi_mu[[m]], cells = cells,
      att = cells$att, inffunc = phi, agg_overall = agg, phi_agg_overall = pagg)
  }
  out <- list(results = results, a = a, a_U = a, tau = tau, eta = eta,
    phi_a_U = t(da), phi_tau = t(dtau), phi_eta = t(deta), actual = actual,
    j = j, phi_j = dj, theta = theta, phi_theta = t(dtheta), treatment = treatment,
    denominators = denominators, snapshot = snapshot, phi_snapshot = t(rows),
    max_root_residual = max(roots), n_imputation_fits = length(bp$groups),
    n_actual_cells = C, n_treatment_fits = length(treatment), n_OR_fits = as.integer(use_or),
    n_observation_logits = 0L, blueprint = bp)
  if (!derivatives) {
    out[c("phi_a_U", "phi_tau", "phi_eta", "phi_j", "phi_theta", "phi_snapshot")] <- NULL
    out$actual <- lapply(out$actual, function(x) { x$phi <- NULL; x })
    out$results <- lapply(out$results, function(x) {
      x[c("phi_mu", "inffunc", "phi_agg_overall")] <- NULL
      x
    })
  }
  out
}

imputation_public_covariance <- function(phi, cluster = NULL) {
  N <- nrow(phi)
  if (!is.null(cluster)) {
    stopifnot(length(cluster) == N, !anyNA(cluster))
    phi <- rowsum(phi, cluster, reorder = FALSE)
  }
  # Same centered original-unit-first / cluster-second mapping as public CP07.
  centered <- sweep(phi, 2L, colMeans(phi), "-")
  crossprod(centered) / N
}

imputation_public_reference <- function(raw, est_method = "dr", imputation_pooling = "separate",
    intercept = TRUE, X_slope = cbind(x = raw$x),
    V_treat = if (intercept) cbind(`(Intercept)` = rep(1, nrow(raw)), X_slope) else X_slope,
    dates = 1:4, control_group = "nevertreated", anticipation = 0L,
    pt_assumption = "post", pt_start = NULL, base_period = "varying",
    derivatives = TRUE, cluster = NULL, blueprint = NULL) {
  bp <- if (is.null(blueprint)) imputation_public_blueprint(raw, imputation_pooling,
    intercept, X_slope, V_treat, dates, control_group, anticipation,
    pt_assumption, pt_start, base_period) else blueprint
  ev <- imputation_public_dense_eval(bp, est_method = est_method, derivatives = derivatives)
  N <- bp$N; nA <- sum(bp$A)
  ids <- paste0("unit-", bp$raw$id)
  if (derivatives) for (m in names(ev$results)) {
    res <- ev$results[[m]]
    rownames(res$inffunc) <- rownames(res$phi_mu) <- rownames(res$phi_agg_overall) <- ids
    res$V_analytical <- imputation_public_covariance(res$inffunc, cluster)
    res$se <- sqrt(diag(res$V_analytical) / N)
    res$V_agg <- imputation_public_covariance(res$phi_agg_overall, cluster)
    res$se_agg <- sqrt(diag(res$V_agg) / N)
    res$core_inffunc <- (nA / N) * res$inffunc
    res$n <- N; res$nA <- nA
    ev$results[[m]] <- res
  }
  ev$metadata <- list(selection = "conditional_mtar", selection_est_method = "imputation",
    imputation_pooling = bp$pooling, gmm_weighting = "identity", original_F_rows = N,
    nA = nA, core_scale = nA / N,
    derivative = if (derivatives) "dense analytic N*d/dc_i" else "not computed",
    inactive_observation_probability_arguments = c("selection_model", "selection_pooling", "selection_xformla"),
    retained_basis = "explicit supplied full-rank dense design; no alias policy invented")
  if (length(ev$results) == 1L) ev <- c(ev, ev$results[[1L]])
  ev
}
