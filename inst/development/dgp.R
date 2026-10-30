# Development simulation support for cdid (not package API).
# Requires seed.R (with_seed).
#
# Staggered-adoption population with known truth (CP01: balanced designs only;
# no attrition, no observation selection).
#
# Population, for each unit i:
#   X_i ~ Bernoulli(p_x)                                binary covariate
#   P(G_i = g | X_i = x) = cohort_prob[x + 1, g]        g in {0, cohorts}; 0 = never treated
#   Z_i ~ N(0, 1), U_i ~ N(0, sd_unit^2)                independent of (X, G)
#   C_i uniform on 1..n_clusters                        cluster label, independent of all else
#   eps_it ~ N(0, sd_eps^2), independent over i and t
#
# Untreated potential outcome, periods t in `periods` (t1 = first period):
#   Y_it(0) = cohort_level[G_i] + time_effect[t] + x_trend * X_i * (t - t1)
#             + z_coef * Z_i + U_i + eps_it
# Treated potential outcome for cohort g (no anticipation):
#   Y_it(g) = Y_it(0) + 1{t >= g} * (theta[g] + delta * (t - g) + kappa * X_i)
# Observed outcome: Y_it = Y_it(G_i), with Y_it = Y_it(0) for never-treated units.
#
# Parallel trends hold conditional on X: untreated trends depend on X only.
# If x_trend is nonzero, they fail unconditionally between cohorts with
# different E[X | G = g]. The default design therefore requires xformla = ~ x.
#
# Population cohort ATT (event time e = t - g, in period-label units):
#   ATT(g, t) = 1{t >= g} * (theta[g] + delta * (t - g) + kappa * E[X | G = g])
#   E[X | G = g] = p_x P(g | 1) / pi_g,  pi_g = (1 - p_x) P(g | 0) + p_x P(g | 1)
# Pre-treatment cells have truth 0 under both base_period conventions.
#
# Aggregation truth: population analogues of did's aggte() estimands. They use
# the population cohort probabilities pi_g over treated cohorts and the same
# cell sets did uses with default arguments (no balance_e, min_e or max_e):
#   simple   : sum over post cells of pi_g ATT(g,t) / sum over post cells of pi_g
#   group    : theta(g) = mean_{t >= g} ATT(g,t); overall = sum pi_g theta(g) / sum pi_g
#   dynamic  : theta(e) = sum_g pi_g ATT(g,g+e) / sum_g pi_g over cohorts with g+e observed;
#              overall = mean of theta(e) over e >= 0
#   calendar : theta(t) = sum_{g <= t} pi_g ATT(g,t) / sum_{g <= t} pi_g for t >= min(cohorts);
#              overall = mean of theta(t)

sim_spec <- function(periods = 1:5,
                     cohorts = c(3, 4, 5),
                     p_x = 0.5,
                     cohort_prob = rbind(c(0.40, 0.20, 0.20, 0.20),
                                         c(0.10, 0.30, 0.30, 0.30)),
                     theta = c(1, 2, 1.5),
                     delta = 0.5,
                     kappa = 1,
                     cohort_level = c(0, 0.5, -0.5, 1),
                     time_effect = 0.25 * seq_along(periods),
                     x_trend = 0.5,
                     z_coef = 1,
                     sd_unit = 1,
                     sd_eps = 1,
                     n_clusters = 20L) {
  periods <- as.numeric(periods)
  cohorts <- as.numeric(cohorts)
  K <- length(cohorts)
  stopifnot(
    length(periods) >= 2L, all(is.finite(periods)), all(periods > 0),
    !is.unsorted(periods, strictly = TRUE),
    K >= 1L, !anyDuplicated(cohorts), all(cohorts %in% periods),
    all(cohorts > min(periods)),
    length(p_x) == 1L, p_x > 0, p_x < 1,
    is.matrix(cohort_prob), dim(cohort_prob) == c(2L, K + 1L),
    all(cohort_prob >= 0),
    all(abs(rowSums(cohort_prob) - 1) < 1e-12),
    length(theta) == K, length(cohort_level) == K + 1L,
    length(time_effect) == length(periods),
    sd_unit >= 0, sd_eps >= 0, n_clusters >= 1
  )
  codes <- c(0, cohorts)
  dimnames(cohort_prob) <- list(c("x0", "x1"), codes)
  pi_g <- (1 - p_x) * cohort_prob[1, ] + p_x * cohort_prob[2, ]
  if (any(pi_g <= 0)) {
    stop("Every cohort, including never-treated (0), needs positive ",
         "population probability.", call. = FALSE)
  }
  list(periods = periods, cohorts = cohorts, p_x = p_x,
       cohort_prob = cohort_prob,
       theta = stats::setNames(theta, cohorts), delta = delta, kappa = kappa,
       cohort_level = stats::setNames(cohort_level, codes),
       time_effect = stats::setNames(time_effect, periods),
       x_trend = x_trend, z_coef = z_coef, sd_unit = sd_unit, sd_eps = sd_eps,
       n_clusters = as.integer(n_clusters))
}

sim_truth <- function(spec) {
  codes <- as.character(c(0, spec$cohorts))
  joint <- rbind((1 - spec$p_x) * spec$cohort_prob[1, ],
                 spec$p_x * spec$cohort_prob[2, ])
  pi_g <- colSums(joint)
  cohort <- data.frame(group = c(0, spec$cohorts), pi = unname(pi_g),
                       mean_x = unname(joint[2, ] / pi_g))

  cells <- expand.grid(t = spec$periods, group = spec$cohorts)[, 2:1]
  gk <- as.character(cells$group)
  post <- cells$t >= cells$group
  mean_x <- stats::setNames(cohort$mean_x, codes)
  cells$att <- ifelse(post, spec$theta[gk] + spec$delta * (cells$t - cells$group) +
                        spec$kappa * mean_x[gk], 0)
  cells$post <- post
  cells$e <- cells$t - cells$group
  w <- unname(pi_g[gk])

  simple <- sum(w[post] * cells$att[post]) / sum(w[post])

  g_att <- vapply(spec$cohorts, function(g) {
    mean(cells$att[cells$group == g & post])
  }, numeric(1))
  wg <- unname(pi_g[as.character(spec$cohorts)])
  group <- list(by = data.frame(group = spec$cohorts, att = g_att),
                overall = sum(wg * g_att) / sum(wg))

  es <- sort(unique(cells$e))
  e_att <- vapply(es, function(e) {
    k <- cells$e == e
    sum(w[k] * cells$att[k]) / sum(w[k])
  }, numeric(1))
  dynamic <- list(by = data.frame(e = es, att = e_att),
                  overall = mean(e_att[es >= 0]))

  ts <- spec$periods[spec$periods >= min(spec$cohorts)]
  t_att <- vapply(ts, function(tt) {
    k <- cells$t == tt & post
    sum(w[k] * cells$att[k]) / sum(w[k])
  }, numeric(1))
  calendar <- list(by = data.frame(t = ts, att = t_att),
                   overall = mean(t_att))

  rownames(cells) <- NULL
  list(cohort = cohort, att_gt = cells, simple = simple, group = group,
       dynamic = dynamic, calendar = calendar)
}

# Draw a sample. panel = TRUE: n units, each observed in every period
# (balanced). panel = FALSE: repeated cross sections with n independent new
# units per period and one row per unit (ids unique across rows).
# Columns: id, period, G, x, z, cluster, y (observed), y0 (untreated potential
# outcome, for checks only; never pass it to an estimator).
sim_draw <- function(spec, n, panel = TRUE, seed) {
  stopifnot(length(n) == 1L, n >= 1, n == round(n), is.logical(panel),
            length(panel) == 1L, !is.na(panel))
  periods <- spec$periods
  T <- length(periods)
  codes <- c(0, spec$cohorts)
  with_seed(seed, {
    N <- if (panel) n else n * T
    x <- stats::rbinom(N, 1L, spec$p_x)
    cum <- t(apply(spec$cohort_prob, 1L, cumsum))
    u <- stats::runif(N)
    idx <- 1L + rowSums(u > cum[x + 1L, -ncol(cum), drop = FALSE])
    G <- codes[idx]
    z <- stats::rnorm(N)
    a <- stats::rnorm(N, 0, spec$sd_unit)
    cl <- sample.int(spec$n_clusters, N, replace = TRUE)
    if (panel) {
      unit <- rep(seq_len(N), each = T)
      period <- rep(periods, times = N)
    } else {
      unit <- seq_len(N)
      period <- rep(periods, each = n)
    }
    eps <- stats::rnorm(length(unit), 0, spec$sd_eps)
    out <- data.frame(id = unit, period = period, G = G[unit], x = x[unit],
                      z = z[unit], cluster = cl[unit])
    gk <- as.character(out$G)
    out$y0 <- unname(spec$cohort_level[gk]) +
      unname(spec$time_effect[as.character(out$period)]) +
      spec$x_trend * out$x * (out$period - periods[1L]) +
      spec$z_coef * out$z + a[unit] + eps
    treated <- out$G > 0 & out$period >= out$G
    effect <- numeric(nrow(out))
    effect[treated] <- unname(spec$theta[gk[treated]]) +
      spec$delta * (out$period[treated] - out$G[treated]) +
      spec$kappa * out$x[treated]
    out$y <- out$y0 + effect
    out[, c("id", "period", "G", "x", "z", "cluster", "y", "y0")]
  })
}
