# Independent D60/D61 reference for every remaining reference-input route.
# Base R only (stats::model.matrix, solve/qr.solve and a plain Newton logit); no
# cdid function, catalog, fit, replay, IF, aggregation or GMM helper is called.
# Source AFTER reference/d60_or_reference.R into the same environment: the
# approved catalog, auxiliary descriptors and recipes of d60_ref_setup() are
# reused; designs, fits and functionals below are written from the accepted
# contracts (CP02 IPW and custom composition, CP04/CP06/CP07 baseline and
# conditional selection routes, pure actual-change imputation and the augmented
# observation correction) with every context reading
#   * X at its recipe date (post links of g: b_g; auxiliary: the earlier pair
#     endpoint), frozen ONCE per formula on the canonical union of the route's
#     required (date, unit) rows (pair rows; selection roster U; full roster F);
#   * W_g = w at t_1 (first_period) or at b_g (base_period) for EVERY unit of
#     the context, including controls and bridges; rates, repeat rates, masses,
#     selection cells, treatment logits and imputation cells keyed by recipe;
#   * ONE global shared-OR slope over all (pair, pool, recipe) blocks; selection
#     and imputation pooling only inside one recipe.
# The point functional is the static accepted functional with these reads; the
# influence functions are central differences of original-unit case weights
# c_i scaling ALL of unit i's dated weights (and its D61 mass):
#   phi_i = N (psi(c + eps e_i) - psi(c - eps e_i)) / (2 eps).
# Fixed-fixture oracle: identified designs, interior probabilities and full-rank
# fits are assumed (no rank/alias/support validator is replicated).
#
# route = list(kind = "baseline" | "corollary1" | "imputation" | "augmented",
#              est = "ipw" | "reg" | "dr" | "custom", custom = callback,
#              sel_model = "pair_logit" | "independent_period_logit",
#              sel_pooling = "separate" | "within_cohort" | "across_cohorts",
#              universe = "F" | "B" | "A", rho = per-unit known ratio (unit order),
#              imp_pooling = "separate" | "within_cohort" | "across_cohorts")
# `variant` negative controls: "static_x" (X at each unit's first date for every
# context), "own_anchor" (control cohorts read w at their own b_k),
# "pool_across_recipes" (selection / imputation pooling ignores recipes).

d60r_setup <- function(d, route, formula = ~x, sformula = NULL, fix_weights = "base_period",
                       control_group = "nevertreated", anticipation = 0, pt_assumption = "post",
                       pt_start = NULL, variant = "none") {
  if (!exists("d60_ref_setup", mode = "function")) stop("source reference/d60_or_reference.R first")
  S <- d60_ref_setup(d, formula = formula, weightsname = "w", fix_weights = fix_weights,
    control_group = control_group, anticipation = anticipation, pt_assumption = pt_assumption,
    pt_start = pt_start)
  if (is.null(route$sel_model)) route$sel_model <- "pair_logit"
  if (is.null(route$sel_pooling)) route$sel_pooling <- "separate"
  if (is.null(route$universe)) route$universe <- "F"
  if (is.null(route$imp_pooling)) route$imp_pooling <- "separate"
  if (is.null(sformula)) sformula <- formula
  S$route <- route; S$variant <- variant
  N <- S$N; T_len <- length(S$cal)
  ids <- S$ids; first_row <- match(ids, d$id)
  wide <- function(v) { m <- matrix(NA, N, T_len); m[cbind(match(d$id, ids), match(as.numeric(d$time), S$cal))] <- v; m }
  changes <- function(v) {
    m <- if (v == "time") matrix(rep(S$cal, each = N), N) else wide(d[[v]])
    any(apply(m, 1L, function(r) { r <- r[!is.na(r)]; length(r) > 1L && any(r != r[1L]) }))
  }
  dated_f <- function(f) length(all.vars(f)) > 0L && any(vapply(all.vars(f), changes, logical(1)))
  uses_s <- route$kind %in% c("corollary1", "augmented")
  any_dated <- dated_f(formula) || (uses_s && dated_f(sformula))
  # Recipe design index: b_g (post) or the earlier endpoint (aux) when ANY
  # formula of the route changes; anchors are those of d60_ref_setup().
  for (l in seq_along(S$links)) {
    L <- S$links[[l]]; b <- match(S$ref_of(L$g), S$cal)
    S$links[[l]]$d <- if (!any_dated) 0L else if (L$role == "post") b else L$is
    S$links[[l]]$rkey <- paste(S$links[[l]]$d, S$links[[l]]$a)
  }
  S$B <- rowSums(S$Rm) >= 1L
  S$U <- switch(route$universe, F = rep(TRUE, N), B = S$B, A = S$A)
  S$rho <- if (is.null(route$rho)) rep(1, N) else as.numeric(route$rho)
  full <- route$kind %in% c("imputation", "augmented")
  # Required (date, unit) rows of each role.
  need <- function(role) {
    out <- list()
    for (L in S$links) for (k in c(L$g, L$K)) {
      rows <- if (role == "x") {
        if (full) which(S$G == k) else which(S$A & S$G == k & S$Rm[, L$is] & S$Rm[, L$it])
      } else if (route$kind == "augmented") which(S$G == k) else which(S$U & S$G == k)
      if (length(rows)) out[[length(out) + 1L]] <- cbind(L$d, rows)
    }
    unique(do.call(rbind, out))
  }
  build <- function(f, req, dated) {
    vars <- all.vars(f)
    if (!dated) req[, 1L] <- 0L
    req <- unique(req); req <- req[order(req[, 1L], req[, 2L]), , drop = FALSE]
    df <- as.data.frame(lapply(stats::setNames(vars, vars), function(v) {
      if (!dated || identical(variant, "static_x")) return(d[[v]][first_row[req[, 2L]]])
      wide(d[[v]])[cbind(req[, 2L], req[, 1L])]
    }), stringsAsFactors = FALSE)
    if (!length(vars)) df <- data.frame(row.names = seq_len(nrow(req)))
    if (anyNA(df)) stop("reference: missing required design value")
    tf <- stats::terms(f)
    Z <- stats::model.matrix(tf, stats::model.frame(tf, df))
    idx <- matrix(NA_integer_, T_len + 1L, N)
    idx[cbind(req[, 1L] + 1L, req[, 2L])] <- seq_len(nrow(req))
    list(req = req, Z = Z, idx = idx, dated = dated, intercept = attr(tf, "intercept") == 1L,
         slope = colnames(Z) != "(Intercept)")
  }
  rx <- need("x")
  same <- identical(paste(deparse(formula), collapse = ""), paste(deparse(sformula), collapse = ""))
  if (uses_s) {
    rs <- need("s")
    if (same) {
      S$dx <- build(formula, rbind(rx, rs), dated_f(formula)); S$ds <- S$dx
    } else {
      S$dx <- build(formula, rx, dated_f(formula)); S$ds <- build(sformula, rs, dated_f(sformula))
    }
  } else S$dx <- build(formula, rx, dated_f(formula))
  S$dated_any <- any_dated
  S
}

# Design rows (all columns) of `des` at recipe date index dd for units `rows`.
d60r_X <- function(des, dd, rows) {
  hit <- des$idx[(if (des$dated) dd else 0L) + 1L, rows]
  if (anyNA(hit)) stop("reference: design row outside the required union")
  des$Z[hit, , drop = FALSE]
}
d60r_z <- function(des, dd, rows) d60r_X(des, dd, rows)[, des$slope, drop = FALSE]

d60r_logit <- function(X, y, w) {
  b <- numeric(ncol(X))
  for (it in 1:200) {
    p <- stats::plogis(drop(X %*% b))
    step <- solve(crossprod(X, X * (w * p * (1 - p))), crossprod(X, w * (y - p)))
    b <- b + drop(step)
    if (max(abs(step)) < 1e-14 * max(1, max(abs(b)))) break
  }
  p <- stats::plogis(drop(X %*% b))
  if (max(abs(crossprod(X, w * (y - p)))) > 1e-10 * max(1, sum(w))) stop("reference: logit did not converge")
  b
}

d60r_wls <- function(X, y, w) drop(solve(crossprod(X, X * w), crossprod(X, w * y)))

d60r_eval <- function(S, cw = rep(1, S$N)) {
  rt <- S$route; N <- S$N; Rm <- S$Rm; G <- S$G; A <- S$A; Y <- S$Y; links <- S$links
  W <- S$W0 * cw; Wu <- S$w_unit * cw
  wa <- function(a) if (a == 0L) Wu else W[, a]
  anchor <- function(L, k) {
    if (identical(S$variant, "own_anchor") && is.finite(k) && k != L$g && L$a > 0L) match(S$ref_of(k), S$cal) else L$a
  }
  obs <- function(k, is, it) which(A & G == k & Rm[, is] & Rm[, it])
  dy <- function(rows, is, it) Y[rows, it] - Y[rows, is]
  icpt <- S$dx$intercept
  # ---- repeat rates (selection / imputation) ---------------------------------
  a_rate <- function(L, k, U) {
    w <- wa(anchor(L, k)); u <- which(U & G == k)
    sum(w[u] * A[u]) / sum(w[u])
  }
  # ---- selection probabilities (conditional selection and augmented) ---------
  sel <- NULL
  if (rt$kind %in% c("corollary1", "augmented")) {
    Usel <- if (rt$kind == "augmented") rep(TRUE, N) else S$U
    period <- identical(rt$sel_model, "independent_period_logit")
    cells <- list(); ckeys <- character()
    for (L in links) for (k in c(L$g, L$K)) {
      dates <- if (period) c(L$is, L$it) else list(c(L$is, L$it))
      for (tt in dates) {
        key <- paste(k, paste(tt, collapse = "-"), L$rkey)
        if (key %in% ckeys) next
        ckeys <- c(ckeys, key)
        fam <- if (identical(S$variant, "pool_across_recipes")) "all" else L$rkey
        blk <- switch(rt$sel_pooling, separate = key, within_cohort = paste(k, fam), across_cohorts = fam)
        cells[[length(cells) + 1L]] <- list(key = key, k = k, tt = tt, d = L$d, a = L$a, blk = blk)
      }
    }
    fits <- list()
    for (bk in unique(vapply(cells, `[[`, character(1), "blk"))) {
      cb <- which(vapply(cells, `[[`, character(1), "blk") == bk)
      XX <- NULL; yy <- numeric(); ww <- numeric()
      for (j in seq_along(cb)) {
        cl <- cells[[cb[j]]]; u <- which(Usel & G == cl$k)
        Zs <- d60r_z(S$ds, cl$d, u)
        E <- if (S$ds$intercept) { m <- matrix(0, length(u), length(cb)); m[, j] <- 1; m } else matrix(0, length(u), 0L)
        resp <- if (length(cl$tt) == 1L) Rm[u, cl$tt] else Rm[u, cl$tt[1L]] & Rm[u, cl$tt[2L]]
        XX <- rbind(XX, cbind(E, Zs)); yy <- c(yy, as.numeric(resp)); ww <- c(ww, wa(cl$a)[u])
      }
      b <- d60r_logit(XX, yy, ww)
      ni <- if (S$ds$intercept) length(cb) else 0L
      for (j in seq_along(cb)) fits[[cells[[cb[j]]]$key]] <- list(alpha = if (ni) b[j] else 0, beta = b[ni + seq_len(length(b) - ni)])
    }
    sel <- function(L, k, rows) {
      pr <- function(tt) {
        f <- fits[[paste(k, paste(tt, collapse = "-"), L$rkey)]]
        stats::plogis(f$alpha + drop(d60r_z(S$ds, L$d, rows) %*% f$beta))
      }
      if (period) pr(L$is) * pr(L$it) else pr(c(L$is, L$it))
    }
  }
  # ---- per-row observation weights v on observed pair rows (non-imputation) --
  vw <- function(L, k, rows) {
    w <- wa(anchor(L, k))
    if (rt$kind == "baseline") {
      u <- which(A & G == k)
      return(w[rows] / (sum(w[u] * (Rm[u, L$is] & Rm[u, L$it])) / sum(w[u])))
    }
    w[rows] * a_rate(L, k, S$U) / (S$rho[rows] * sel(L, k, rows))
  }
  # ---- imputation predictions j_c (pure imputation and augmented) ------------
  jhat <- NULL; ahat <- NULL
  if (rt$kind %in% c("imputation", "augmented")) {
    cells <- list(); ckeys <- character()
    for (L in links) for (k in c(L$g, L$K)) {
      key <- paste(k, L$is, L$it, L$rkey)
      if (key %in% ckeys) next
      ckeys <- c(ckeys, key)
      fam <- if (identical(S$variant, "pool_across_recipes")) "all" else L$rkey
      grp <- switch(rt$imp_pooling, separate = key, within_cohort = paste(k, fam), across_cohorts = fam)
      cells[[length(cells) + 1L]] <- list(key = key, k = k, L = L, grp = grp)
    }
    coef <- list()
    for (gk in unique(vapply(cells, `[[`, character(1), "grp"))) {
      cg <- which(vapply(cells, `[[`, character(1), "grp") == gk)
      XX <- NULL; yy <- numeric(); ww <- numeric()
      for (j in seq_along(cg)) {
        cl <- cells[[cg[j]]]; L <- cl$L; o <- obs(cl$k, L$is, L$it)
        E <- if (icpt) { m <- matrix(0, length(o), length(cg)); m[, j] <- 1; m } else matrix(0, length(o), 0L)
        XX <- rbind(XX, cbind(E, (L$it - L$is) * d60r_z(S$dx, L$d, o)))
        yy <- c(yy, dy(o, L$is, L$it)); ww <- c(ww, wa(L$a)[o])
      }
      b <- d60r_wls(XX, yy, ww); ni <- if (icpt) length(cg) else 0L
      for (j in seq_along(cg)) coef[[cells[[cg[j]]]$key]] <- list(alpha = if (ni) b[j] else 0, beta = b[ni + seq_len(length(b) - ni)])
    }
    jhat <- function(L, k, rows) {
      f <- coef[[paste(k, L$is, L$it, L$rkey)]]
      f$alpha + (L$it - L$is) * drop(d60r_z(S$dx, L$d, rows) %*% f$beta)
    }
    ahat <- function(L, k) a_rate(L, k, rep(TRUE, N))
  }
  # Augmented change Z = j + S/q (Delta - j) on all F rows of cohort k.
  zaug <- function(L, k, rows) {
    j <- jhat(L, k, rows)
    if (rt$kind != "augmented") return(j)
    s <- Rm[rows, L$is] & Rm[rows, L$it]
    out <- j
    if (any(s)) {
      q <- sel(L, k, rows[s])
      out[s] <- j[s] + (dy(rows[s], L$is, L$it) - j[s]) / q
    }
    out
  }
  full <- rt$kind %in% c("imputation", "augmented")
  # ---- ONE global shared OR over every (pair, pool, recipe) block ------------
  or_fit <- NULL
  if (rt$est %in% c("reg", "dr")) {
    blocks <- unique(vapply(links, function(L) paste(L$is, L$it, paste(L$K, collapse = "|"), L$rkey), character(1)))
    XX <- NULL; YY <- numeric(); WW <- numeric(); done <- character()
    for (L in links) {
      bk <- paste(L$is, L$it, paste(L$K, collapse = "|"), L$rkey)
      if (bk %in% done) next
      done <- c(done, bk); bi <- match(bk, blocks)
      for (k in L$K) {
        if (full) {
          rows <- which(G == k); yv <- zaug(L, k, rows); q <- wa(anchor(L, k))[rows] * ahat(L, k)
        } else {
          rows <- obs(k, L$is, L$it); yv <- dy(rows, L$is, L$it); q <- vw(L, k, rows)
        }
        E <- if (icpt) { m <- matrix(0, length(rows), length(blocks)); m[, bi] <- 1; m } else matrix(0, length(rows), 0L)
        XX <- rbind(XX, cbind(E, (L$it - L$is) * d60r_z(S$dx, L$d, rows))); YY <- c(YY, yv); WW <- c(WW, q)
      }
    }
    th <- d60r_wls(XX, YY, WW); nb <- if (icpt) length(blocks) else 0L
    or_fit <- list(blocks = blocks, alpha = if (nb) th[seq_len(nb)] else numeric(length(blocks)),
      beta = th[nb + seq_len(length(th) - nb)])
  }
  or_pred <- function(L, rows) {
    bi <- match(paste(L$is, L$it, paste(L$K, collapse = "|"), L$rkey), or_fit$blocks)
    or_fit$alpha[bi] + (L$it - L$is) * drop(d60r_z(S$dx, L$d, rows) %*% or_fit$beta)
  }
  # ---- link functionals ------------------------------------------------------
  mu <- vapply(links, function(L) {
    ks <- c(L$g, L$K)
    if (!full) {
      rows <- unlist(lapply(ks, obs, is = L$is, it = L$it)); rows <- sort(rows)
      kk <- G[rows]; D <- as.numeric(kk == L$g); fw <- numeric(length(rows))
      for (k in ks) fw[kk == k] <- vw(L, k, rows[kk == k])
      dv <- dy(rows, L$is, L$it)
      if (rt$est == "custom") {
        X <- d60r_X(S$dx, L$d, rows)
        res <- rt$custom(y1 = Y[rows, L$it], y0 = Y[rows, L$is], D = D, covariates = X, i.weights = fw, inffunc = TRUE)
        return(as.numeric(res$ATT))
      }
      if (rt$est == "reg") {
        r <- dv - or_pred(L, rows)
        return(sum(fw * D * r) / sum(fw * D))
      }
      X <- d60r_X(S$dx, L$d, rows)
      odds <- exp(drop(X %*% d60r_logit(X, D, fw)))
      wp <- fw * D; wm <- fw * (1 - D) * odds
      r <- if (rt$est == "dr") dv - or_pred(L, rows) else dv
      return(sum(wp * r) / sum(wp) - sum(wm * r) / sum(wm))
    }
    # Pure imputation / augmented: synthetic roster masses w a_k on all F rows.
    rows <- sort(which(G %in% ks)); kk <- G[rows]; D <- as.numeric(kk == L$g)
    ws <- numeric(length(rows)); Zc <- numeric(length(rows))
    for (k in ks) {
      ws[kk == k] <- wa(anchor(L, k))[rows[kk == k]] * ahat(L, k)
      Zc[kk == k] <- zaug(L, k, rows[kk == k])
    }
    Tg <- sum(ws * D * Zc) / sum(ws * D)
    bpred <- if (rt$est %in% c("reg", "dr")) or_pred(L, rows) else NULL
    if (rt$est == "reg") return(Tg - sum(ws * D * bpred) / sum(ws * D))
    X <- d60r_X(S$dx, L$d, rows)
    odds <- exp(drop(X %*% d60r_logit(X, D, ws)))
    wm <- ws * (1 - D) * odds
    Co <- sum(wm * Zc) / sum(wm)
    if (rt$est == "ipw") return(Tg - Co)
    (Tg - sum(ws * D * bpred) / sum(ws * D)) - (Co - sum(wm * bpred) / sum(wm))
  }, numeric(1))
  d60r_report(S, mu, cw)
}

# Identity GMM on post links, reported cells and D61 aggregates: the
# d60_or_reference.R construction applied to the route's link estimates.
d60r_report <- function(S, mu, cw) {
  links <- S$links
  post <- which(vapply(links, `[[`, character(1), "role") == "post")
  cname <- function(L, u) {
    ref <- S$ref_of(L$g); ctl <- paste(L$K, collapse = "|")
    if (u == ref) return("ZERO")
    if (identical(S$pt_assumption, "all") && u >= S$pt_start && u < ref && min(L$K) > u + S$anticipation) return("ZERO")
    if (u > ref && min(L$K) > u + S$anticipation) return(paste0("shared:", L$g, ":", u))
    paste0("priv:", L$g, ":", ctl, ":", u)
  }
  sn <- vapply(post, function(l) cname(links[[l]], S$cal[links[[l]]$is]), character(1))
  tn <- vapply(post, function(l) cname(links[[l]], S$cal[links[[l]]$it]), character(1))
  coords <- setdiff(unique(c(sn, tn)), "ZERO")
  Dm <- matrix(0, length(post), length(coords), dimnames = list(NULL, coords))
  for (r in seq_along(post)) {
    if (sn[r] != "ZERO") Dm[r, sn[r]] <- Dm[r, sn[r]] - 1
    if (tn[r] != "ZERO") Dm[r, tn[r]] <- Dm[r, tn[r]] + 1
  }
  psi <- stats::setNames(as.vector(qr.solve(crossprod(Dm), crossprod(Dm, mu[post]))), coords)
  cells <- list()
  for (g in S$targets) for (u in S$cal[S$cal >= g]) {
    cn <- paste0("shared:", g, ":", u)
    cells[[length(cells) + 1L]] <- data.frame(g = g, t = u, att = if (cn %in% coords) unname(psi[cn]) else NA_real_)
  }
  cells <- do.call(rbind, cells)
  V <- S$V0 * cw
  mass <- vapply(S$targets, function(g) { u <- which(S$G == g); sum(V[u]) / S$N }, numeric(1))
  mg <- mass[match(cells$g, S$targets)]
  agg <- c(simple = sum(mg * cells$att) / sum(mg))
  th_g <- vapply(S$targets, function(g) mean(cells$att[cells$g == g]), numeric(1))
  agg["group"] <- sum(mass * th_g) / sum(mass)
  ev <- cells$t - cells$g; ue <- sort(unique(ev))
  th_e <- vapply(ue, function(e) { r <- which(ev == e); sum(mg[r] * cells$att[r]) / sum(mg[r]) }, numeric(1))
  agg["dynamic"] <- mean(th_e[ue >= 0])
  ut <- sort(unique(cells$t))
  th_t <- vapply(ut, function(t) { r <- which(cells$t == t); sum(mg[r] * cells$att[r]) / sum(mg[r]) }, numeric(1))
  agg["calendar"] <- mean(th_t)
  list(mu = mu, psi = psi, cells = cells, mass = mass, agg = agg)
}

# `units`: original units whose IF rows are differentiated (default all).
d60r_if <- function(S, eps = 1e-4, units = seq_len(S$N)) {
  base <- d60r_eval(S)
  vec <- function(e) c(e$mu, e$cells$att, e$agg)
  v0 <- vec(base)
  phi <- matrix(NA_real_, S$N, length(v0))
  for (i in units) {
    cp <- rep(1, S$N); cp[i] <- 1 + eps
    cm <- rep(1, S$N); cm[i] <- 1 - eps
    phi[i, ] <- S$N * (vec(d60r_eval(S, cp)) - vec(d60r_eval(S, cm))) / (2 * eps)
  }
  nm <- length(base$mu); nc <- nrow(base$cells)
  list(base = base, phi_mu = phi[, seq_len(nm), drop = FALSE], phi_cells = phi[, nm + seq_len(nc), drop = FALSE],
       phi_agg = phi[, nm + nc + seq_len(4L), drop = FALSE])
}
