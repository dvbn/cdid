# Independent D60/D61 reference for the shared one-period outcome-change OR
# baseline route (est_method = "reg", selection = "none") with reference-date
# covariates and fixed-anchor changing sampling weights. Base R only
# (stats::model.matrix, stats::lm.wfit, qr.solve); no cdid function, catalog,
# fit, replay, IF, aggregation or GMM helper is called.
#   * catalog: approved comparison catalog (CP02 Definition 1) and the varying /
#     universal auxiliary descriptors, from outcome observation patterns only;
#   * recipes: every post link of target g (bridges and controls included)
#     extracts X at b_g = max{d : d + delta < g}; an auxiliary report extracts at
#     the earlier endpoint of its ordered pair; W_g = w at t_1 (first_period) or
#     at b_g (base_period) for EVERY unit of target g's context;
#   * design: the formula is evaluated ONCE on the canonical union of required
#     (date, id) rows = observed pair rows of every support a link or block uses;
#   * one stacked WLS with a dummy per unique (pair, pool, recipe) block and the
#     regressor L_e * z(C^r) (L_e = calendar-index steps), weights W / r with
#     r = sum(W S) / sum(W A) per (pair, cohort, recipe);
#   * identity GMM by incidence least squares; D61 masses V = w_t1 by unit;
#   * influence functions: central differences of multiplicative original-unit
#     case weights c_i scaling ALL of unit i's dated weights,
#       phi_i = N (psi(c + eps e_i) - psi(c - eps e_i)) / (2 eps).
# Fixed-fixture oracle: supported positive weights and identified designs are
# assumed (no rank validator); simple aliasing is not replicated.

# `variant` builds deliberately WRONG functionals for negative controls:
#   "own_anchor"      each control cohort k reads w at its own b_k (not the target's);
#   "summed_x"        regressor sum_{u in (s,t]} x_u (evolving X summed; linear ~x only);
#   "merge_recipes"   blocks pooled across recipes (one intercept per pair/pool);
#   "duplicate_aux"   auxiliary-recipe blocks counted twice (multiplier two);
#   "drop_no_t1"      units without a first-calendar weight silently dropped;
#   "frozen_beta"     slope frozen at the full-sample root (omits xi_beta);
#   "frozen_shares"   D61 masses frozen at the full-sample values (omits shares).
# `frame_extra` (id, G) appends frame-only roster units after the panel units.
d60_ref_setup <- function(d, formula = ~x, weightsname = "w", fix_weights = "base_period",
                          control_group = "nevertreated", anticipation = 0, tname = "time",
                          pt_assumption = "post", pt_start = NULL, variant = "none", frame_extra = NULL) {
  if (identical(variant, "drop_no_t1")) {
    t1 <- min(d$time)
    keep_id <- d$id[d$time == t1 & !is.na(d[[weightsname]])]
    d <- d[d$id %in% keep_id, , drop = FALSE]
  }
  ids <- c(unique(d$id), frame_extra$id)
  cal <- sort(unique(as.numeric(d$time)))
  N <- length(ids); T_len <- length(cal)
  ii <- match(d$id, ids); jj <- match(as.numeric(d$time), cal)
  wide <- function(v) { m <- matrix(NA, N, T_len); m[cbind(ii, jj)] <- v; m }
  Y <- wide(d$y); Rm <- !is.na(Y)
  G <- as.numeric(d$G[match(ids, d$id)])
  if (!is.null(frame_extra)) G[match(frame_extra$id, ids)] <- as.numeric(frame_extra$G)
  G[G == 0] <- Inf
  A <- rowSums(Rm) >= 2L
  W0 <- if (is.null(weightsname)) matrix(1, N, T_len) else wide(d[[weightsname]])
  w_unit <- apply(W0, 1L, function(r) r[!is.na(r)][1L])
  cohorts <- sort(unique(G))
  ref_of <- function(g) { pre <- cal[cal + anticipation < g]; if (length(pre)) max(pre) else NA_real_ }
  occ <- function(k, is, it) any(A & G == k & Rm[, is] & Rm[, it])
  targets <- sort(unique(G[is.finite(G)]))
  targets <- targets[vapply(targets, function(g) any(cal >= g) && !is.na(ref_of(g)), logical(1))]
  # CP02 Definition 1 catalog.
  rows <- list()
  for (g in targets) for (is in seq_len(T_len - 1L)) for (it in (is + 1L):T_len) {
    if (!occ(g, is, it)) next
    Ce <- cohorts[vapply(cohorts, occ, logical(1), is = is, it = it)]
    Ie <- Ce[is.infinite(Ce)]; Fe <- Ce[is.finite(Ce) & Ce > g]
    pools <- if (control_group == "nevertreated") list(Ie) else {
      H <- cal[cal >= g & cal >= cal[it]]
      c(list(c(Ie, Fe[cal[it] + anticipation < Fe])), lapply(H, function(h) c(Ie, Fe[h + anticipation < Fe])))
    }
    pools <- unique(lapply(pools, sort))
    for (K in pools) if (length(K)) rows[[length(rows) + 1L]] <- list(role = "post", g = g, is = is, it = it, K = K)
  }
  # Auxiliary descriptors (union of varying and universal), dedup by (g, s, t, K).
  aux <- list(); seen <- character()
  for (bp in c("varying", "universal")) for (g in targets) {
    ref <- ref_of(g)
    times <- if (bp == "varying") cal[-1L] else cal
    for (u in times[times < g]) {
      from <- if (bp == "varying") cal[match(u, cal) - 1L] else ref
      if (bp == "universal" && u == ref) next
      s <- min(u, from); t <- max(u, from)
      K <- if (control_group == "nevertreated") cohorts[is.infinite(cohorts)] else cohorts[cohorts != g & cohorts > t + anticipation]
      key <- paste(g, s, t, paste(K, collapse = ","))
      if (key %in% seen || !length(K)) next
      seen <- c(seen, key)
      is <- match(s, cal); it <- match(t, cal)
      if (!occ(g, is, it) || !all(vapply(K, occ, logical(1), is = is, it = it))) next
      aux[[length(aux) + 1L]] <- list(role = "aux", g = g, is = is, it = it, K = sort(K))
    }
  }
  links <- c(rows, aux)
  # Recipes: (design date index or 0 when no formula variable changes, weight
  # anchor index or 0 when weights do not change within any unit).
  vars <- all.vars(formula)
  changes <- function(v) {
    m <- if (v == tname) matrix(rep(cal, each = N), N) else wide(d[[v]])
    any(apply(m, 1L, function(r) { r <- r[!is.na(r)]; length(r) > 1L && any(r != r[1L]) }))
  }
  dated_x <- length(vars) > 0L && any(vapply(vars, changes, logical(1)))
  dated_w <- !is.null(weightsname) && any(apply(W0, 1L, function(r) { r <- r[!is.na(r)]; length(r) > 1L && any(r != r[1L]) }))
  for (l in seq_along(links)) {
    L <- links[[l]]; b <- match(ref_of(L$g), cal)
    links[[l]]$d <- if (!dated_x) 0L else if (L$role == "post") b else L$is
    links[[l]]$a <- if (!dated_w) 0L else if (fix_weights == "first_period") 1L else b
    links[[l]]$rkey <- paste(links[[l]]$d, links[[l]]$a)
    links[[l]]$block <- paste(L$is, L$it, paste(L$K, collapse = "|"),
      if (identical(variant, "merge_recipes")) "merged" else links[[l]]$rkey)
  }
  blocks <- unique(vapply(links, `[[`, character(1), "block"))
  # Canonical required design rows: observed pair rows of every used support.
  req <- list()
  for (L in links) for (k in c(L$g, L$K)) {
    obs <- which(A & G == k & Rm[, L$is] & Rm[, L$it])
    req[[length(req) + 1L]] <- cbind(L$d, obs)
  }
  req <- unique(do.call(rbind, req))
  req <- req[order(req[, 1L], req[, 2L]), , drop = FALSE]
  frame_row <- match(ids, d$id)
  union_df <- as.data.frame(lapply(stats::setNames(vars, vars), function(v) {
    if (v == tname) return(cal[pmax(1L, req[, 1L])])
    if (!dated_x) return(d[[v]][frame_row[req[, 2L]]])
    m <- wide(d[[v]])
    if (any(req[, 1L] == 0L)) stop("reference: static rows in a dated design")
    m[cbind(req[, 2L], req[, 1L])]
  }), stringsAsFactors = FALSE)
  if (!length(vars)) union_df <- data.frame(row.names = seq_len(nrow(req)))
  if (anyNA(union_df)) stop("reference: missing required design value")
  tf <- stats::terms(formula)
  Zu <- stats::model.matrix(tf, stats::model.frame(tf, union_df))
  icpt <- attr(tf, "intercept") == 1L
  slope <- Zu[, colnames(Zu) != "(Intercept)", drop = FALSE]
  if (identical(variant, "summed_x")) xw <- wide(d$x)
  list(d = d, ids = ids, cal = cal, N = N, Y = Y, Rm = Rm, G = G, A = A, W0 = W0, cohorts = cohorts,
       targets = targets, links = links, blocks = blocks, req = req, slope = slope, icpt = icpt,
       variant = variant, xw = if (identical(variant, "summed_x")) xw else NULL,
       pt_assumption = pt_assumption, pt_start = if (is.null(pt_start)) min(cal) else pt_start,
       anticipation = anticipation, ref_of = ref_of, dated_x = dated_x, dated_w = dated_w,
       fix_weights = fix_weights, w_unit = w_unit, V0 = if (dated_w) W0[, 1L] else w_unit)
}

d60_ref_z <- function(S, dd, i) {
  hit <- match(paste(dd, i), paste(S$req[, 1L], S$req[, 2L]))
  if (anyNA(hit)) stop("reference: design row outside the required union")
  S$slope[hit, , drop = FALSE]
}

d60_ref_eval <- function(S, cw = rep(1, S$N), frozen = NULL) {
  W <- S$W0 * cw; Wu <- S$w_unit * cw
  links <- S$links; nb <- length(S$blocks); p <- ncol(S$slope)
  support <- function(L, k) {
    is <- L$is; it <- L$it
    Sv <- S$Rm[, is] & S$Rm[, it]
    a <- L$a
    if (identical(S$variant, "own_anchor") && is.finite(k) && k != L$g) a <- match(S$ref_of(k), S$cal)
    wa <- if (a == 0L) Wu else W[, a]
    inA <- which(S$A & S$G == k)
    rate <- sum(wa[inA] * Sv[inA]) / sum(wa[inA])
    obs <- which(S$A & S$G == k & Sv)
    z <- if (identical(S$variant, "summed_x")) {
      cbind(x = rowSums(S$xw[obs, (is + 1L):it, drop = FALSE]) / (it - is))
    } else d60_ref_z(S, L$d, obs)
    list(rows = obs, q = wa[obs] / rate, dy = S$Y[obs, it] - S$Y[obs, is], z = z, L = it - is)
  }
  XX <- NULL; YY <- numeric(); WW <- numeric(); done <- character()
  for (L in links) {
    bi <- match(L$block, S$blocks)
    if (L$block %in% done) next
    done <- c(done, L$block)
    for (k in L$K) {
      s <- support(L, k)
      E <- if (S$icpt) { m <- matrix(0, length(s$rows), nb); m[, bi] <- 1; m } else matrix(0, length(s$rows), 0L)
      times <- if (identical(S$variant, "duplicate_aux") && L$role == "aux") 2L else 1L
      for (r in seq_len(times)) { XX <- rbind(XX, cbind(E, s$L * s$z)); YY <- c(YY, s$dy); WW <- c(WW, s$q) }
    }
  }
  fit <- stats::lm.wfit(XX, YY, WW, tol = 1e-12)
  coef <- unname(fit$coefficients)
  alpha <- if (S$icpt) coef[seq_len(nb)] else numeric(nb)
  beta <- coef[(if (S$icpt) nb else 0L) + seq_len(p)]
  if (!is.null(frozen$beta)) {
    # Frozen-slope negative control: block intercepts refit given the frozen beta.
    beta <- frozen$beta
    if (S$icpt) alpha <- unname(stats::lm.wfit(XX[, seq_len(nb), drop = FALSE],
      YY - as.vector(XX[, nb + seq_len(p), drop = FALSE] %*% beta), WW, tol = 1e-12)$coefficients)
  }
  mu <- vapply(links, function(L) {
    s <- support(L, L$g); bi <- match(L$block, S$blocks)
    sum(s$q * (s$dy - alpha[bi] - s$L * as.vector(s$z %*% beta))) / sum(s$q)
  }, numeric(1))
  # Identity GMM on post links (CP02 coordinate naming; no PT-all zeros).
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
  mass <- vapply(S$targets, function(g) { u <- which(S$G == g); if (anyNA(V[u])) NA_real_ else sum(V[u]) / S$N }, numeric(1))
  names(mass) <- S$targets
  if (!is.null(frozen$mass)) mass <- frozen$mass
  mg <- mass[as.character(cells$g)]
  agg <- c(simple = sum(mg * cells$att) / sum(mg))
  ug <- S$targets
  th_g <- vapply(ug, function(g) mean(cells$att[cells$g == g]), numeric(1))
  agg["group"] <- sum(mass * th_g) / sum(mass)
  ev <- cells$t - cells$g; ue <- sort(unique(ev))
  th_e <- vapply(ue, function(e) { r <- which(ev == e); sum(mg[r] * cells$att[r]) / sum(mg[r]) }, numeric(1))
  agg["dynamic"] <- mean(th_e[ue >= 0])
  ut <- sort(unique(cells$t))
  th_t <- vapply(ut, function(t) { r <- which(cells$t == t); sum(mg[r] * cells$att[r]) / sum(mg[r]) }, numeric(1))
  agg["calendar"] <- mean(th_t)
  list(mu = mu, alpha = alpha, beta = beta, psi = psi, cells = cells, mass = mass, agg = agg,
       dynamic = stats::setNames(th_e, ue), group = stats::setNames(th_g, ug))
}

d60_ref_if <- function(S, eps = 1e-4) {
  base <- d60_ref_eval(S)
  frozen <- list(beta = if (identical(S$variant, "frozen_beta")) base$beta,
                 mass = if (identical(S$variant, "frozen_shares")) base$mass)
  vec <- function(e) c(e$mu, e$cells$att, e$agg)
  v0 <- vec(base)
  phi <- matrix(0, S$N, length(v0))
  for (i in seq_len(S$N)) {
    cp <- rep(1, S$N); cp[i] <- 1 + eps
    cm <- rep(1, S$N); cm[i] <- 1 - eps
    phi[i, ] <- S$N * (vec(d60_ref_eval(S, cp, frozen)) - vec(d60_ref_eval(S, cm, frozen))) / (2 * eps)
  }
  nm <- length(base$mu); nc <- nrow(base$cells)
  list(base = base, phi_mu = phi[, seq_len(nm), drop = FALSE], phi_cells = phi[, nm + seq_len(nc), drop = FALSE],
       phi_agg = phi[, nm + nc + seq_len(4L), drop = FALSE])
}
