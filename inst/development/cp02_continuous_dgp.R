# Development simulation support for cdid (not package API).
# Requires seed.R (with_seed). Standalone: does not use dgp.R.
#
# CP02 continuous-covariate, heavy-tail family exactly as proposed in
# plan/evidence/cp02_checks/continuous_heavytail_development_proposal.md.
# The design is fixed; there are no tuning arguments. Definitions and proofs:
# docs/theory/cp02_continuous_heavytail_dgp.tex.
#
# Population, for each unit i and period t in 1:5:
#   G_i in {3, 4, 5, 0} (0 = never treated), population mass 1/4 each
#   X_i | G = 3  : log(1/2 + 3V/2) / log(2), V ~ U(0, 1);
#                  density (2 log 2 / 3) 2^x on [-1, 1], mean 5/3 - 1/log 2
#   X_i | G != 3 : U(-1, 1)
#   Z1, Z2, Z3 iid U(-1, 1), U_i ~ N(0, 1), all independent of (G, X)
#   eps_it = (1 + X/4) sqrt(3/5) T_it, T_it iid Student t5 (built as
#            N / sqrt(chi2_5 / 5) from standard normals)
#   Y_it(0) = alpha_G + t/4 + 0.4 sin(pi X)(t - 1) + 0.2 X^2 (t^2 - 1)
#             + 0.1 Z2^2 (t - 1) + Z1 + U_i + eps_it + 1{t = 1} s_G(X)
#   alpha = 0, 1/2, -1/2, 1 for G = 0, 3, 4, 5
#   s_0 = 0, s_3 = 2 + 0.7 X, s_4 = -1 + 0.2 X, s_5 = 1 - 0.4 X
#   Y_it(G) = Y_it(0) + 1{G > 0, t >= G} (theta_G + 0.5 (t - G) + X),
#   theta = 1, 2, 3/2 for G = 3, 4, 5 (no anticipation)
# Observation pattern, independent of (G, X, Z, U, T): each of the ten
# two-date subsets of 1:5 w.p. 0.08, singleton {3} w.p. 0.10, none w.p. 0.10.
# Y_it is observed iff t is in the unit's pattern. Nothing is repaired,
# redrawn, trimmed or imputed, and no cell count is forced.
#
# Truth (cp02_continuous_truth) is analytic and uses population masses only;
# realized or selected counts never define it. Pre-treatment causal effects
# are zero, but the period-1 shifts make pre-period DiD reporting estimands
# that involve period 1 nonzero; no pre-period reporting truth is declared.

cp02_continuous_patterns <- function() {
  pairs <- utils::combn(5L, 2L)
  data.frame(
    pattern = c(paste(pairs[1L, ], pairs[2L, ], sep = "-"), "3", "none"),
    type = c(rep("pair", 10L), "singleton", "none"),
    date1 = c(pairs[1L, ], 3L, NA_integer_),
    date2 = c(pairs[2L, ], NA_integer_, NA_integer_),
    slots = c(rep(4L, 10L), 5L, 5L),
    prob = c(rep(0.08, 10L), 0.10, 0.10),
    stringsAsFactors = FALSE
  )
}

cp02_continuous_settings <- function() {
  list(
    family = "cp02_continuous_heavytail",
    status = paste("Development validation DGP proposal; no estimator,",
                   "coverage, calibration or certification role."),
    periods = 1:5,
    cohorts = c(3, 4, 5),
    never_treated_code = 0,
    cohort_mass = c("3" = 0.25, "4" = 0.25, "5" = 0.25, "0" = 0.25),
    alpha = c("0" = 0, "3" = 0.5, "4" = -0.5, "5" = 1),
    theta = c("3" = 1, "4" = 2, "5" = 1.5),
    delta = 0.5,
    shift = data.frame(G = c(0, 3, 4, 5), intercept = c(0, 2, -1, 1),
                       slope_x = c(0, 0.7, 0.2, -0.4)),
    noise = list(df = 5, scale = sqrt(3 / 5), sd_factor = "1 + x/4"),
    patterns = cp02_continuous_patterns(),
    covariates = list(
      x_only = "~ x",
      x_plus_z = "~ x + z1 + z2 + z3",
      propensity_correct = paste("Both specifications are correct logits for",
                                 "the declared future-control operators",
                                 "(truth$propensity); Z coefficients are 0."),
      outcome_regression_linear_correct = FALSE,
      fitted = FALSE
    ),
    rng = list(
      kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
      sample.kind = "Rejection",
      draw_order = c("runif(n): cohort", "runif(n): X", "runif(n, -1, 1): z1",
                     "runif(n, -1, 1): z2", "runif(n, -1, 1): z3",
                     "rnorm(n): unit effect", "rnorm(5n): t5 numerators",
                     "5 x rnorm(5n): t5 chi-square terms",
                     "runif(n): observation pattern"),
      note = paste("Bit-identical only for the same seed, n and R RNG",
                   "implementation; every block depends on n. No",
                   "nonoverlap claim across seeds or streams.")
    )
  )
}

# Inverse-CDF maps from one U(0,1) draw (exposed for checks).
cp02_continuous_cohort <- function(v) {
  c(3, 4, 5, 0)[pmin(floor(4 * v), 3) + 1]
}

cp02_continuous_x_quantile <- function(v, g3) {
  ifelse(rep_len(g3, length(v)), log(0.5 + 1.5 * v) / log(2), 2 * v - 1)
}

cp02_continuous_pattern_index <- function(v) {
  slot <- pmin(floor(50 * v), 49)
  ifelse(slot < 40, slot %/% 4 + 1, ifelse(slot < 45, 11, 12))
}

cp02_continuous_truth <- function() {
  s <- cp02_continuous_settings()
  mass <- s$cohort_mass
  mean_x <- c("3" = 5 / 3 - 1 / log(2), "4" = 0, "5" = 0, "0" = 0)
  cohort <- data.frame(G = c(3, 4, 5, 0), mass = unname(mass),
                       mean_x = unname(mean_x))

  # Post cells; ATT = const + mx_coef * E[X | G = 3].
  cells <- expand.grid(t = s$periods, group = s$cohorts)[, 2:1]
  cells <- cells[cells$t >= cells$group, ]
  gk <- as.character(cells$group)
  cells$e <- cells$t - cells$group
  cells$att_const <- unname(s$theta[gk]) + s$delta * cells$e
  cells$att_mx_coef <- as.numeric(cells$group == 3)
  cells$att <- cells$att_const + unname(mean_x[gk])
  rownames(cells) <- NULL

  pre <- expand.grid(t = s$periods, group = s$cohorts)[, 2:1]
  pre <- pre[pre$t < pre$group, ]
  pre$causal_att <- 0
  pre$reporting_did_truth <- NA_real_
  rownames(pre) <- NULL

  # Aggregate weights over post cells, population masses only.
  w <- unname(mass[gk])
  wts <- list()
  add <- function(agg, comp, k, weight) {
    wts[[length(wts) + 1L]] <<- data.frame(
      aggregate = agg, component = comp, group = cells$group[k],
      t = cells$t[k], weight = weight)
  }
  add("simple", "overall", seq_len(nrow(cells)), w / sum(w))
  wg <- mass[as.character(s$cohorts)] / sum(mass[as.character(s$cohorts)])
  for (g in s$cohorts) {
    k <- which(cells$group == g)
    add("group", paste0("g=", g), k, rep(1 / length(k), length(k)))
    add("group", "overall", k, unname(wg[as.character(g)]) / length(k))
  }
  es <- sort(unique(cells$e))
  for (e in es) {
    k <- which(cells$e == e)
    add("dynamic", paste0("e=", e), k, w[k] / sum(w[k]))
    add("dynamic", "overall", k, w[k] / sum(w[k]) / length(es))
  }
  ts <- sort(unique(cells$t))
  for (tt in ts) {
    k <- which(cells$t == tt)
    add("calendar", paste0("t=", tt), k, w[k] / sum(w[k]))
    add("calendar", "overall", k, w[k] / sum(w[k]) / length(ts))
  }
  wts <- do.call(rbind, wts)
  key <- match(paste(wts$group, wts$t), paste(cells$group, cells$t))
  comb <- function(v) tapply(wts$weight * v[key],
                             paste(wts$aggregate, wts$component), sum)
  vals <- comb(cells$att)
  aggregates <- data.frame(
    aggregate = sub(" .*", "", names(vals)),
    component = sub("^[^ ]* ", "", names(vals)),
    att_const = unname(comb(cells$att_const)),
    att_mx_coef = unname(comb(cells$att_mx_coef)),
    att = unname(vals), stringsAsFactors = FALSE)

  # Declared future-control operators for post cells: NT = {0};
  # NYT = {c : c > t} U {0}. Equal masses cancel in the propensity.
  ops <- list()
  for (i in seq_len(nrow(cells))) {
    g <- cells$group[i]
    tt <- cells$t[i]
    for (ctl in c("NT", "NYT")) {
      pool <- if (ctl == "NT") 0 else c(s$cohorts[s$cohorts > tt], 0)
      k <- length(pool)
      ops[[length(ops) + 1L]] <- data.frame(
        group = g, t = tt, control = ctl, pool = paste(pool, collapse = ","),
        pool_size = k,
        intercept = if (g == 3) log(4 * log(2) / (3 * k)) else -log(k),
        slope_x = if (g == 3) log(2) else 0, slope_z = 0,
        correct_x_only = TRUE, correct_x_plus_z = TRUE,
        stringsAsFactors = FALSE)
    }
  }
  propensity <- do.call(rbind, ops)

  pat <- s$patterns
  list(
    cohort = cohort,
    att_gt = cells,
    pre = pre,
    aggregates = aggregates,
    aggregate_weights = wts,
    propensity = propensity,
    propensity_note = paste("Correctness is proved only for the declared",
                            "future-control operators. Pools containing G=3",
                            "for targets 4 or 5 (earlier-treated reporting",
                            "comparisons) are not covered."),
    density_ratio = list(formula = "(4 log 2 / 3) 2^x",
                         lower = 2 * log(2) / 3, upper = 8 * log(2) / 3),
    selection = list(pattern_prob = stats::setNames(pat$prob, pat$pattern),
                     p_pair = sum(pat$prob[pat$type == "pair"]),
                     q_pair_given_g = 0.08 / sum(pat$prob[pat$type == "pair"]),
                     expected_units_per_cohort_pair = "n / 50"),
    noise = list(var_standardized_t5 = 1, fourth_standardized_t5 = 9,
                 conditional_var = "(1 + x/4)^2",
                 conditional_fourth = "9 (1 + x/4)^4",
                 note = "Absolute moments of order >= 5 are infinite."),
    pre_note = paste("Causal pre-treatment effects are 0. Pre-period DiD",
                     "reporting estimands involving period 1 need not be 0",
                     "because of cohort-specific period-1 shifts; none is",
                     "declared as truth.")
  )
}

cp02_check_n <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n != round(n) || n < 1) {
    stop("`n` must be a single whole number >= 1.", call. = FALSE)
  }
  if (5 * n > .Machine$integer.max) {
    stop("`n` too large: 5 * n rows must not exceed .Machine$integer.max.",
         call. = FALSE)
  }
  as.integer(n)
}

# Draw n units. Returns list(data, oracle, settings): `data` holds only what
# an analyst sees (full long calendar with y = NA where unobserved, unit frame,
# realized G x pattern counts); `oracle` holds potential outcomes, noise and
# the analytic truth; never pass oracle columns to an estimator.
cp02_continuous_draw <- function(n, seed) {
  n <- cp02_check_n(n)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed) || abs(seed) > .Machine$integer.max) {
    stop("`seed` must be a single whole number in integer range.",
         call. = FALSE)
  }
  s <- cp02_continuous_settings()
  pat <- s$patterns
  m <- 5L * n
  with_seed(seed, {
    G <- cp02_continuous_cohort(stats::runif(n))
    x <- cp02_continuous_x_quantile(stats::runif(n), G == 3)
    z1 <- stats::runif(n, -1, 1)
    z2 <- stats::runif(n, -1, 1)
    z3 <- stats::runif(n, -1, 1)
    unit_effect <- stats::rnorm(n)
    t5 <- stats::rnorm(m)
    chi2 <- 0
    for (k in 1:5) chi2 <- chi2 + stats::rnorm(m)^2
    t5 <- t5 / sqrt(chi2 / 5)
    pidx <- cp02_continuous_pattern_index(stats::runif(n))
  })

  id <- rep(seq_len(n), each = 5L)
  period <- rep(s$periods, times = n)
  Gi <- G[id]
  xi <- x[id]
  gk <- as.character(Gi)
  shift <- ifelse(period == 1,
                  s$shift$intercept[match(Gi, s$shift$G)] +
                    s$shift$slope_x[match(Gi, s$shift$G)] * xi, 0)
  eps <- (1 + xi / 4) * s$noise$scale * t5
  y0 <- unname(s$alpha[gk]) + period / 4 + 0.4 * sin(pi * xi) * (period - 1) +
    0.2 * xi^2 * (period^2 - 1) + 0.1 * z2[id]^2 * (period - 1) + z1[id] +
    unit_effect[id] + eps + shift
  treated <- Gi > 0 & period >= Gi
  effect <- numeric(m)
  effect[treated] <- unname(s$theta[gk[treated]]) +
    s$delta * (period[treated] - Gi[treated]) + xi[treated]
  y_g <- y0 + effect

  seen <- matrix(FALSE, nrow(pat), 5L)
  for (j in seq_len(nrow(pat))) {
    seen[j, stats::na.omit(c(pat$date1[j], pat$date2[j]))] <- TRUE
  }
  observed <- seen[cbind(pidx[id], period)]

  units <- data.frame(id = seq_len(n), G = G, x = x, z1 = z1, z2 = z2,
                      z3 = z3, pattern = pat$pattern[pidx],
                      pattern_type = pat$type[pidx],
                      n_obs = rowSums(seen)[pidx], stringsAsFactors = FALSE)
  panel <- data.frame(id = id, period = period, G = Gi, x = xi, z1 = z1[id],
                      z2 = z2[id], z3 = z3[id], pattern = units$pattern[id],
                      pattern_type = units$pattern_type[id],
                      observed = observed, y = ifelse(observed, y_g, NA_real_),
                      stringsAsFactors = FALSE)
  counts <- as.data.frame(table(
    G = factor(G, levels = c(0, 3, 4, 5)),
    pattern = factor(units$pattern, levels = pat$pattern)),
    responseName = "units", stringsAsFactors = FALSE)
  counts$G <- as.numeric(counts$G)

  s$n <- n
  s$seed <- seed
  s$rng$r_version <- R.version.string
  list(
    data = list(panel = panel, units = units, counts = counts),
    oracle = list(
      panel = data.frame(id = id, period = period, y0 = y0, y_g = y_g,
                         effect = effect, eps = eps, t5 = t5, shift = shift),
      units = data.frame(id = seq_len(n), unit_effect = unit_effect),
      truth = cp02_continuous_truth()
    ),
    settings = s
  )
}
