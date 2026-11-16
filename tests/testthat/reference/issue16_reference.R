# Independent reference for issue #16 (missing_reference = "exclude") on the
# dated baseline IPW route with never-treated controls. Base R only (stats::
# model.matrix and a plain Newton logit); no cdid function, catalog, reader,
# fit, replay, IF, GMM or aggregation helper is called.
#
# Contract written out directly:
#   * b_g = last stored date d with d + anticipation < g;
#   * E_g = units whose raw covariates of `formula` are all supplied at b_g (an
#     absent (unit, b_g) row or an NA value is unknown; the outcome plays no
#     role); E_g never reads another date;
#   * every main link (s, t) of target g reads X at b_g and W_g = w at b_g
#     (base_period) or at the first stored date (first_period) for every unit;
#     its observers, rate denominators and risk set are restricted to E_g:
#       r_k = H(W_g I_k E_g S_st) / H(W_g I_k E_g A),
#       mu = sum wp dy / sum wp - sum wm dy / sum wm,
#       wp = W_g D / r_G, wm = W_g (1 - D) / r_G * exp(x' beta),
#     beta the W_g / r_G-weighted logit of D on X over the restricted risk set;
#   * a link exists when target g and the never-treated cohort both have an
#     eligible unit observed at s and t (per-context occupancy);
#   * identity GMM per target: coordinates ZERO at b_g, shared:g:u for u > b_g
#     and private:g:u for u < b_g; ATT(g, t) = shared:g:t (t >= g);
#   * D61 masses m_g = H(V I_g E_g), V = w at the first stored date; simple,
#     group, dynamic and calendar aggregates as in did.
# Influence functions are central differences of original-unit case weights c_i
# scaling ALL of unit i's weights (and its mass):
#   phi_i = N (psi(c + eps e_i) - psi(c - eps e_i)) / (2 eps).
# Fixed-fixture oracle: identified designs and interior probabilities assumed.

i16_setup <- function(d, formula = ~x, fix_weights = "base_period", anticipation = 0, exclude = TRUE) {
  ids <- unique(d$id); N <- length(ids); cal <- sort(unique(as.numeric(d$time))); Tn <- length(cal)
  ui <- match(d$id, ids); ti <- match(as.numeric(d$time), cal)
  wide <- function(v) { m <- matrix(NA, N, Tn); m[cbind(ui, ti)] <- v; m }
  G <- as.numeric(d$G[match(ids, d$id)]); G[G == 0] <- Inf
  Y <- wide(as.numeric(d$y)); R <- !is.na(Y); A <- rowSums(R) >= 2L
  W <- if (is.null(d$w)) matrix(1, N, Tn) else wide(as.numeric(d$w))
  vars <- all.vars(formula)
  Xraw <- lapply(stats::setNames(vars, vars), function(v) wide(d[[v]]))
  targets <- sort(unique(G[is.finite(G) & G <= max(cal)]))
  ref_of <- function(g) max(cal[cal + anticipation < g])
  E <- lapply(targets, function(g) {
    bi <- match(ref_of(g), cal)
    ok <- rep(TRUE, N)
    for (v in vars) ok <- ok & !is.na(Xraw[[v]][, bi])
    if (!exclude && !all(ok)) stop("reference: strict policy with a missing reference value")
    ok
  })
  names(E) <- as.character(targets)
  links <- list()
  for (g in targets) {
    e <- E[[as.character(g)]]
    for (s in seq_len(Tn - 1L)) for (t in seq.int(s + 1L, Tn)) {
      obs <- A & e & R[, s] & R[, t]
      if (any(obs & G == g) && any(obs & is.infinite(G)))
        links[[length(links) + 1L]] <- list(g = g, is = s, it = t, b = match(ref_of(g), cal))
    }
  }
  list(d = d, ids = ids, N = N, cal = cal, G = G, Y = Y, R = R, A = A, W = W, V = W[, 1L], Xraw = Xraw,
    formula = formula, fix_weights = fix_weights, targets = targets, E = E, links = links, ref_of = ref_of)
}

i16_logit <- function(X, y, w) {
  b <- numeric(ncol(X))
  for (it in 1:200) {
    p <- stats::plogis(drop(X %*% b))
    step <- solve(crossprod(X, X * (w * p * (1 - p))), crossprod(X, w * (y - p)))
    b <- b + drop(step)
    if (max(abs(step)) < 1e-14 * max(1, max(abs(b)))) break
  }
  b
}

# Link points, cells and aggregates under original-unit case weights `cw`.
i16_eval <- function(S, cw = rep(1, S$N)) {
  G <- S$G; R <- S$R; A <- S$A; Y <- S$Y
  mu <- vapply(S$links, function(L) {
    e <- S$E[[as.character(L$g)]]
    a <- if (identical(S$fix_weights, "first_period")) 1L else L$b
    Wg <- S$W[, a] * cw
    obs <- A & e & R[, L$is] & R[, L$it]
    rate <- function(k) { u <- which(A & e & G == k); sum(Wg[u] * obs[u]) / sum(Wg[u]) }
    rows <- which(obs & (G == L$g | is.infinite(G)))
    D <- as.numeric(G[rows] == L$g)
    fw <- Wg[rows] / ifelse(D == 1, rate(L$g), rate(Inf))
    df <- as.data.frame(lapply(S$Xraw, function(m) m[rows, L$b]))
    X <- stats::model.matrix(S$formula, stats::model.frame(S$formula, df))
    odds <- exp(drop(X %*% i16_logit(X, D, fw)))
    dy <- Y[rows, L$it] - Y[rows, L$is]
    wp <- fw * D; wm <- fw * (1 - D) * odds
    sum(wp * dy) / sum(wp) - sum(wm * dy) / sum(wm)
  }, numeric(1))
  # Identity GMM, one decoupled least-squares system per target.
  cells <- list()
  for (g in S$targets) {
    li <- which(vapply(S$links, `[[`, numeric(1), "g") == g)
    if (!length(li)) next
    b <- S$links[[li[1L]]]$b
    us <- sort(unique(unlist(lapply(S$links[li], function(L) c(L$is, L$it)))))
    us <- setdiff(us, b)
    Dm <- matrix(0, length(li), length(us))
    for (r in seq_along(li)) {
      L <- S$links[[li[r]]]
      if (L$is != b) Dm[r, match(L$is, us)] <- -1
      if (L$it != b) Dm[r, match(L$it, us)] <- 1
    }
    psi <- as.vector(qr.solve(crossprod(Dm), crossprod(Dm, mu[li])))
    for (u in us[S$cal[us] >= g]) cells[[length(cells) + 1L]] <- data.frame(g = g, t = S$cal[u], att = psi[match(u, us)])
  }
  cells <- do.call(rbind, cells)
  V <- S$V * cw
  mass <- vapply(S$targets, function(g) sum(V[S$G == g & S$E[[as.character(g)]]]) / S$N, numeric(1))
  mg <- mass[match(cells$g, S$targets)]
  agg <- c(simple = sum(mg * cells$att) / sum(mg))
  tg <- unique(cells$g)
  th_g <- vapply(tg, function(g) mean(cells$att[cells$g == g]), numeric(1))
  agg["group"] <- sum(mass[match(tg, S$targets)] * th_g) / sum(mass[match(tg, S$targets)])
  ev <- cells$t - cells$g; ue <- sort(unique(ev))
  th_e <- vapply(ue, function(x) { r <- which(ev == x); sum(mg[r] * cells$att[r]) / sum(mg[r]) }, numeric(1))
  agg["dynamic"] <- mean(th_e[ue >= 0])
  ut <- sort(unique(cells$t))
  th_t <- vapply(ut, function(x) { r <- which(cells$t == x); sum(mg[r] * cells$att[r]) / sum(mg[r]) }, numeric(1))
  agg["calendar"] <- mean(th_t)
  list(mu = mu, cells = cells, mass = mass, agg = agg)
}

i16_if <- function(S, eps = 1e-4) {
  base <- i16_eval(S)
  vec <- function(e) c(e$mu, e$cells$att, e$agg)
  phi <- matrix(NA_real_, S$N, length(vec(base)))
  for (i in seq_len(S$N)) {
    cp <- rep(1, S$N); cp[i] <- 1 + eps
    cm <- rep(1, S$N); cm[i] <- 1 - eps
    phi[i, ] <- S$N * (vec(i16_eval(S, cp)) - vec(i16_eval(S, cm))) / (2 * eps)
  }
  nm <- length(base$mu); nc <- nrow(base$cells)
  list(base = base, phi_mu = phi[, seq_len(nm), drop = FALSE], phi_cells = phi[, nm + seq_len(nc), drop = FALSE],
    phi_agg = phi[, nm + nc + seq_len(4L), drop = FALSE])
}
