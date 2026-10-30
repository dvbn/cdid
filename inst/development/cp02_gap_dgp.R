# D27-adopted CP02 gap laws: development only, not package API or acceptance.
# Source seed.R, dgp.R and cp02_baseline_variants_dgp.R first. No fitting,
# sim_spec mutation, selection repair, attrition correction or new dependency.
# Literal dates are retained. Oracle columns must never be estimator inputs.

cp02_gap_catalog <- function() {
  data.frame(
    law = c("CAL_B", "CAL_S", "MON_B", "MON_S", "NONADJ_B", "NONADJ_S",
            "NO_NT_B", "NO_NT_S"),
    family = rep(c("CAL", "MON", "NONADJ", "NO_NT"), each = 2L),
    innovation_family = rep(c("gaussian", "centered_exponential"), 4L),
    anticipation = c(1, 2, rep(0, 6L)),
    first_date_shift = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
    has_never_treated = c(rep(TRUE, 6L), FALSE, FALSE),
    status = "D27 adopted for development; independent review/integration pending",
    stringsAsFactors = FALSE)
}

cp02_gap_settings <- function(law) {
  cat <- cp02_gap_catalog()
  if (!is.character(law) || length(law) != 1L || is.na(law) ||
      !law %in% cat$law) stop("`law` must name one of the eight D27 laws.", call. = FALSE)
  row <- cat[match(law, cat$law), ]
  b <- sim_spec()                    # read canonical components, never mutate b
  cal <- if (row$family == "CAL") c(0, 3, 7, 9, 12) else as.numeric(1:5)
  gs <- if (row$family == "CAL") c(6, 10, 12) else b$cohorts
  codes <- c(0, gs)
  joint <- rbind((1 - b$p_x) * b$cohort_prob[1L, ], b$p_x * b$cohort_prob[2L, ])
  colnames(joint) <- as.character(codes)
  if (!row$has_never_treated) {
    joint <- joint[, -1L, drop = FALSE] / sum(joint[, -1L])
    codes <- gs
  }
  # Reuse accepted pattern and shift definitions where identical.
  if (row$family %in% c("CAL", "NO_NT")) {
    pat <- cp02_baseline_patterns(if (row$family == "CAL") "general_pairs" else "balanced")
    mask <- matrix(FALSE, nrow(pat), 5L)
    for (i in seq_len(nrow(pat))) {
      j <- if (pat$type[i] == "full") 1:5 else
        as.integer(stats::na.omit(c(pat$date1[i], pat$date2[i])))
      mask[i, j] <- TRUE
    }
  } else {
    dates <- if (row$family == "MON") lapply(2:5, seq_len) else
      list(c(1, 3), c(1, 4), c(1, 5), c(2, 4), c(2, 5), c(3, 5))
    mask <- t(vapply(dates, function(v) 1:5 %in% v, logical(5L)))
    pat <- data.frame(pattern = vapply(dates, paste, "", collapse = "-"),
      type = if (row$family == "MON") "prefix" else "pair",
      date1 = vapply(dates, min, 0), date2 = vapply(dates, max, 0),
      n_obs = lengths(dates), eligible_ge2 = TRUE, slots = 1L,
      prob = 1 / length(dates), stringsAsFactors = FALSE)
  }
  # The stored mask and labels use literal dates, not index relabeling in data.
  pat$date1 <- cal[pat$date1]
  pat$date2 <- cal[pat$date2]
  pat$pattern <- vapply(seq_len(nrow(pat)), function(i) {
    if (any(mask[i, ])) paste(cal[mask[i, ]], collapse = "-") else "none"
  }, "")
  dimnames(mask) <- list(pat$pattern, as.character(cal))
  shift <- cp02_default_period1_shift(gs, active = row$first_date_shift)
  shift <- shift[shift$G %in% codes, ]
  refs <- stats::setNames(vapply(gs, function(g) max(cal[cal + row$anticipation < g]), 0), gs)
  # Direct categorical sampling: row-major joint (X0 cohorts, X1 cohorts).
  pop <- data.frame(x = rep(0:1, each = length(codes)), G = rep(codes, 2L),
                    prob = as.vector(t(joint)))
  slots <- if (row$has_never_treated) c(4, 2, 2, 2, 1, 3, 3, 3) else c(2, 2, 2, 3, 3, 3)
  pop$slots <- slots
  list(law = law, family = row$family, status = row$status,
    periods = cal, cohorts = gs, cohort_codes = codes, never_treated_code = 0,
    has_never_treated = row$has_never_treated, anticipation = row$anticipation,
    reference = refs, joint = joint, population = pop,
    cohort_level = stats::setNames(unname(b$cohort_level)[match(codes, c(0, gs))], codes),
    theta = stats::setNames(unname(b$theta), gs), delta = b$delta, kappa = b$kappa,
    x_trend = b$x_trend, z_coef = b$z_coef,
    time_effect = stats::setNames(unname(b$time_effect), cal),
    innovation_family = row$innovation_family,
    time_scale = if (row$innovation_family == "gaussian") rep(1, 5L) else (2 + 1:5) / 4,
    first_date_shift = shift, patterns = pat, observation_mask = mask,
    sampling_scheme = if (row$has_never_treated) "iid_target" else "direct_conditioned_population",
    rng = list(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection",
      draw_order = c("runif(n): direct finite (X,G) categorical",
        "rnorm(n): Z", "rnorm(n): U",
        if (row$innovation_family == "gaussian") "rnorm(5*n): eps_raw" else
          "runif(5*n): accepted centered-exponential quantile",
        "runif(n): independent observation pattern"),
      note = "Replay requires the same R RNG implementation; no stream-nonoverlap claim."),
    scope = paste("Population truth; D10 originals plus D13 fixed horizons only.",
      "Conditional untreated increments from the actual references are common across cohorts.",
      "No pre-PT-all, estimator/inference, final campaign or released-path compatibility claim."))
}

cp02_gap_slot_index <- function(v, slots) {
  stopifnot(is.numeric(v), all(is.finite(v)), all(v >= 0), all(v <= 1))
  total <- sum(slots)
  slot <- pmin(floor(v * total), total - 1)
  1L + rowSums(outer(slot, utils::head(cumsum(slots), -1L), ">="))
}

cp02_gap_population <- function(v, law) {
  p <- cp02_gap_settings(law)$population
  p[cp02_gap_slot_index(v, p$slots), c("x", "G"), drop = FALSE]
}

cp02_gap_pattern_index <- function(v, law) {
  s <- cp02_gap_settings(law)
  if (s$family == "CAL") return(cp02_baseline_pattern_index(v, "general_pairs"))
  if (s$family == "NO_NT") return(cp02_baseline_pattern_index(v, "balanced"))
  cp02_gap_slot_index(v, s$patterns$slots)
}

# Both potential paths at assigned G, conditional on explicit innovations.
# This pure finite function also exposes formula checks without random draws.
cp02_gap_potential <- function(law, G, period, x, z = 0, unit_effect = 0, eps_raw = 0) {
  s <- cp02_gap_settings(law)
  args <- list(G, period, x, z, unit_effect, eps_raw)
  n <- max(lengths(args))
  if (!all(lengths(args) %in% c(1L, n)) || n < 1L ||
      !all(vapply(args, function(a) is.numeric(a) && all(is.finite(a)), TRUE)))
    stop("Potential inputs must be finite numeric scalars or equal-length vectors.", call. = FALSE)
  G <- rep_len(G, n); period <- rep_len(period, n); x <- rep_len(x, n)
  if (any(!G %in% s$cohort_codes) || any(!period %in% s$periods) || any(!x %in% 0:1))
    stop("Potential inputs outside the literal finite support.", call. = FALSE)
  j <- match(period, s$periods); gk <- as.character(G)
  sh <- match(G, s$first_date_shift$G)
  shift <- (period == s$periods[1L]) *
    (s$first_date_shift$intercept[sh] + s$first_date_shift$slope_x[sh] * x)
  eps_scale <- s$time_scale[j]
  eps <- eps_scale * rep_len(eps_raw, n)
  y0 <- unname(s$cohort_level[gk]) + unname(s$time_effect[as.character(period)]) +
    s$x_trend * x * (period - s$periods[1L]) + s$z_coef * rep_len(z, n) +
    rep_len(unit_effect, n) + eps + shift
  effect <- numeric(n)
  post <- G > 0 & period >= G
  ant <- G > 0 & period < G & period >= G - s$anticipation
  effect[ant] <- 1 / 4 + x[ant] / 2
  effect[post] <- unname(s$theta[gk[post]]) + s$delta * (period[post] - G[post]) + s$kappa * x[post]
  data.frame(y0 = y0, y_g = y0 + effect, effect = effect,
    eps = eps, eps_raw = rep_len(eps_raw, n), eps_scale = eps_scale, shift = shift)
}

# q_k f_k(x) is proportional to the declared joint probability because masks
# are independent, with the same P(A) for every population cell.
cp02_gap_operator_truth <- function(law, g, pool, u, reference = NULL) {
  s <- cp02_gap_settings(law)
  if (length(g) != 1L || !g %in% s$cohorts || !length(pool) ||
      anyDuplicated(pool) || any(!pool %in% s$cohort_codes) || g %in% pool)
    stop("Invalid target or fixed control pool.", call. = FALSE)
  if (is.null(reference)) reference <- unname(s$reference[as.character(g)])
  if (length(u) != 1L || length(reference) != 1L ||
      !u %in% s$periods || !reference %in% s$periods)
    stop("Operator endpoints must be literal dates.", call. = FALSE)
  fg <- s$joint[, as.character(g)] / sum(s$joint[, as.character(g)])
  total <- 0
  for (x in 0:1) {
    weights <- s$joint[x + 1L, as.character(pool)]
    weights <- weights / sum(weights)
    inc <- function(k) cp02_gap_potential(law, k, u, x)$y_g -
      cp02_gap_potential(law, k, reference, x)$y_g
    total <- total + fg[x + 1L] * (inc(g) - sum(weights * vapply(pool, inc, 0)))
  }
  unname(total)
}

# Population pair rates, not realized sample occupancies.
cp02_gap_pair_rates <- function(law) {
  s <- cp02_gap_settings(law); p <- utils::combn(s$periods, 2L)
  data.frame(s = p[1L, ], t = p[2L, ], prob = vapply(seq_len(ncol(p)), function(i) {
    sum(s$patterns$prob[s$observation_mask[, match(p[1L, i], s$periods)] &
      s$observation_mask[, match(p[2L, i], s$periods)]])
  }, 0))
}

# Fixed five-date population truth catalog; not an estimator catalog builder.
# Enumerates precisely the adopted D10/D13 family for these eight laws. Every
# private endpoint retains its (g,pool,u) identity. No covariance/fit is used.
cp02_gap_comparison_truth <- function(law, control = c("NT", "NYT")) {
  control <- match.arg(control); s <- cp02_gap_settings(law)
  pairs <- cp02_gap_pair_rates(law); pairs <- pairs[pairs$prob > 0, ]
  rows <- list(); vertices <- list()
  put <- function(g, u, pool, force_post = FALSE) {
    ref <- s$reference[as.character(g)]
    if (!force_post && u == ref) return(NA_character_)
    qualified <- force_post || (u + s$anticipation >= g &&
      all(pool == 0 | pool > u + s$anticipation))
    role <- if (qualified) if (u >= g) "post" else "anticipation" else
      if (u < g - s$anticipation) "pre_private" else
        if (u < g) "anticipation_private" else "post_private"
    pool_key <- paste(sort(pool), collapse = ",")
    key <- if (qualified) paste("ATT", g, u, sep = "|") else
      paste("private", g, pool_key, u, sep = "|")
    if (is.null(vertices[[key]])) {
      beta <- if (qualified) {
        px <- s$joint[, as.character(g)] / sum(s$joint[, as.character(g)])
        sum(px * cp02_gap_potential(law, g, u, 0:1)$effect)
      } else cp02_gap_operator_truth(law, g, pool, u)
      vertices[[key]] <<- data.frame(coordinate = key, group = g, t = u,
        role = role, pool = if (qualified) "" else pool_key, beta = beta)
    }
    key
  }
  for (g in s$cohorts) for (u in s$periods[s$periods >= g]) put(g, u, numeric(), TRUE)
  for (g in s$cohorts) for (i in seq_len(nrow(pairs))) {
    a <- pairs$s[i]; b <- pairs$t[i]
    nt <- s$cohort_codes[s$cohort_codes == 0]
    cutoffs <- if (control == "NT") b else unique(c(b, s$periods[s$periods >= g & s$periods >= b]))
    pools <- list(); original <- logical(); horizons <- list()
    for (h in cutoffs) {
      K <- if (control == "NT") nt else c(nt, s$cohorts[s$cohorts > g & s$cohorts > h + s$anticipation])
      if (!length(K)) next
      key <- paste(sort(K), collapse = ",")
      if (is.null(pools[[key]])) {
        pools[[key]] <- K; original[key] <- FALSE; horizons[[key]] <- numeric()
      }
      if (h == b) original[key] <- TRUE
      # h=b may be both original and a D13 horizon.
      if (control == "NYT" && h >= g) horizons[[key]] <- c(horizons[[key]], h)
    }
    for (key in names(pools)) {
      K <- pools[[key]]
      from <- put(g, a, K); to <- put(g, b, K)
      bs <- cp02_gap_operator_truth(law, g, K, a)
      bt <- cp02_gap_operator_truth(law, g, K, b)
      r <- data.frame(group = g, s = a, t = b, reference = unname(s$reference[as.character(g)]),
        pool = key, original = unname(original[key]), from = from, to = to,
        beta_s = bs, beta_t = bt, mu = bt - bs, pair_probability = pairs$prob[i])
      r$horizons <- list(horizons[[key]])
      rows[[length(rows) + 1L]] <- r
    }
  }
  v <- do.call(rbind, unname(vertices)); rownames(v) <- NULL
  r <- if (length(rows)) do.call(rbind, rows) else data.frame(
    group = numeric(), s = numeric(), t = numeric(), reference = numeric(), pool = character(),
    original = logical(), from = character(), to = character(), beta_s = numeric(),
    beta_t = numeric(), mu = numeric(), pair_probability = numeric(), horizons = I(list()))
  W <- matrix(0, nrow(r), nrow(v), dimnames = list(NULL, v$coordinate))
  for (i in seq_len(nrow(r))) {
    if (!is.na(r$from[i])) W[i, r$from[i]] <- -1
    if (!is.na(r$to[i])) W[i, r$to[i]] <- 1
  }
  rk <- if (nrow(W)) qr(W)$rank else 0L
  available <- vapply(seq_len(ncol(W)), function(j) {
    e <- numeric(ncol(W)); e[j] <- 1
    (if (nrow(W)) qr(rbind(W, e))$rank else 1L) == rk
  }, TRUE)
  v$available <- available
  list(control = control, rows = r, coordinates = v, incidence = W, rank = rk)
}

# Auxiliary reporting operators from the approved copied-development rule.
# A known population functional can have no direct observed pair. This flag
# does not override post-level availability in the chained population catalog.
cp02_gap_reporting_truth <- function(law) {
  s <- cp02_gap_settings(law); pairs <- cp02_gap_pair_rates(law); out <- list()
  for (base in c("varying", "universal")) for (control in c("NT", "NYT"))
    for (g in s$cohorts) for (j in if (base == "varying") 2:5 else 1:5) {
      u <- s$periods[j]; ref <- unname(s$reference[as.character(g)])
      a <- if (base == "varying" && u < g) s$periods[j - 1L] else ref
      cutoff <- max(a, u) + s$anticipation
      K <- s$cohort_codes[s$cohort_codes == 0 | (s$cohort_codes != g & s$cohort_codes > cutoff)]
      if (control == "NT") K <- K[K == 0]
      normalized <- base == "universal" && u == ref
      defined <- length(K) > 0L
      val <- if (normalized) 0 else if (defined) cp02_gap_operator_truth(law, g, K, u, a) else NA_real_
      k <- which(pairs$s == min(a, u) & pairs$t == max(a, u))
      prob <- if (a == u) 1 else pairs$prob[k]
      out[[length(out) + 1L]] <- data.frame(base_period = base, control = control,
        group = g, t = u, base = a, cutoff = cutoff, pool = paste(K, collapse = ","),
        normalized = normalized, operator_defined = defined,
        direct_pair_supported = normalized || (defined && prob > 0), truth = val)
    }
  do.call(rbind, out)
}

# Exact population maps over the explicitly supplied post-cell selector B.
cp02_gap_aggregate_truth <- function(cells, cohort, target_set) {
  weights <- list(); w <- cohort$pi[match(cells$group, cohort$group)]
  add <- function(agg, component, k, a) {
    weights[[length(weights) + 1L]] <<- data.frame(target_set = target_set,
      aggregate = agg, component = component, group = cells$group[k], t = cells$t[k], weight = a)
  }
  add("simple", "overall", seq_len(nrow(cells)), w / sum(w))
  gs <- unique(cells$group); wg <- cohort$pi[match(gs, cohort$group)]
  for (i in seq_along(gs)) {
    k <- which(cells$group == gs[i])
    add("group", paste0("g=", gs[i]), k, rep(1 / length(k), length(k)))
    add("group", "overall", k, rep(wg[i] / sum(wg) / length(k), length(k)))
  }
  for (agg in c("dynamic", "calendar")) {
    index <- if (agg == "dynamic") cells$e else cells$t
    levels <- sort(unique(index))
    for (z in levels) {
      k <- which(index == z); a <- w[k] / sum(w[k])
      add(agg, paste0(if (agg == "dynamic") "e=" else "t=", z), k, a)
      add(agg, "overall", k, a / length(levels))
    }
  }
  weights <- do.call(rbind, weights)
  hit <- match(paste(weights$group, weights$t), paste(cells$group, cells$t))
  keys <- unique(paste(weights$aggregate, weights$component))
  values <- do.call(rbind, lapply(keys, function(key) {
    k <- paste(weights$aggregate, weights$component) == key
    data.frame(target_set = target_set, aggregate = weights$aggregate[which(k)[1L]],
      component = weights$component[which(k)[1L]], att = sum(weights$weight[k] * cells$att[hit[k]]),
      available = all(cells$available[hit[k]]))
  }))
  list(values = values, weights = weights)
}

cp02_gap_truth <- function(law) {
  s <- cp02_gap_settings(law); mass <- colSums(s$joint)
  cohort <- data.frame(group = s$cohort_codes, pi = unname(mass), mean_x = unname(s$joint[2L, ] / mass))
  causal <- expand.grid(t = s$periods, group = s$cohorts)[, 2:1]
  causal$e <- causal$t - causal$group
  causal$att <- vapply(seq_len(nrow(causal)), function(i) {
    g <- causal$group[i]; px <- s$joint[, as.character(g)] / mass[as.character(g)]
    sum(px * cp02_gap_potential(law, g, causal$t[i], 0:1)$effect)
  }, 0)
  causal$role <- ifelse(causal$t >= causal$group, "post", ifelse(
    causal$t >= causal$group - s$anticipation, "causal_anticipation", "causal_pre"))
  comparison <- list(NT = cp02_gap_comparison_truth(law, "NT"), NYT = cp02_gap_comparison_truth(law, "NYT"))
  post <- causal[causal$role == "post", ]; rownames(post) <- NULL
  for (ctl in c("NT", "NYT")) {
    v <- comparison[[ctl]]$coordinates
    post[[paste0("available_", ctl)]] <- v$available[match(paste("ATT", post$group, post$t, sep = "|"), v$coordinate)]
  }
  post$available <- post$available_NYT
  ag <- cp02_gap_aggregate_truth(post, cohort, "full_post")
  restricted <- if (!s$has_never_treated)
    cp02_gap_aggregate_truth(post[post$available, ], cohort, "explicit_restricted_available_post") else NULL
  list(cohort = cohort, att_gt = post, causal = causal,
    anticipation = causal[causal$role == "causal_anticipation", ],
    comparison = comparison, reporting = cp02_gap_reporting_truth(law),
    aggregates = ag$values, aggregate_weights = ag$weights, restricted = restricted,
    selection = list(pair_rates = cp02_gap_pair_rates(law),
      p_eligible_ge2 = sum(s$patterns$prob[s$patterns$eligible_ge2]), independent_of_all_outcomes = TRUE),
    noise = list(mean_eps = rep(0, 5L), var_eps = s$time_scale^2,
      third_eps = if (s$innovation_family == "gaussian") rep(0, 5L) else 2 * s$time_scale^3,
      fourth_eps = if (s$innovation_family == "gaussian") rep(3, 5L) else 9 * s$time_scale^4),
    note = paste("Analytic population truth, independent of observed counts.",
      "Availability is within the adopted fixed-operator class, not fitted validity.",
      "Causal anticipation is separate from varying/universal reporting operators.",
      "NO_NT retains all six truths; full and explicitly restricted maps are distinct."))
}

cp02_gap_validate_draw_inputs <- function(law, n, seed) {
  cp02_gap_settings(law)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n != round(n) ||
      n > .Machine$integer.max %/% 5L) stop("`n` must be whole, >=1, with 5*n in integer range.", call. = FALSE)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != round(seed) ||
      abs(seed) > .Machine$integer.max) stop("`seed` must be whole and in integer range.", call. = FALSE)
  invisible(TRUE)
}

# Sole top-level random constructor. Full long calendar, observed rows and
# frame-only units coexist; only oracle holds the potential outcomes/noise.
cp02_gap_draw <- function(law, n, seed) {
  cp02_gap_validate_draw_inputs(law, n, seed)
  n <- as.integer(n); s <- cp02_gap_settings(law)
  with_seed(seed, {
    pop <- cp02_gap_population(stats::runif(n), law)
    z <- stats::rnorm(n); unit_effect <- stats::rnorm(n)
    eps_raw <- if (s$innovation_family == "gaussian") stats::rnorm(5L * n) else
      cp02_centered_exp_quantile(stats::runif(5L * n))
    pidx <- cp02_gap_pattern_index(stats::runif(n), law)
  })
  id <- rep(seq_len(n), each = 5L); period <- rep(s$periods, n)
  pat <- s$patterns
  potential <- cp02_gap_potential(law, pop$G[id], period, pop$x[id], z[id], unit_effect[id], eps_raw)
  mask <- s$observation_mask[pidx, , drop = FALSE]
  observed <- as.vector(t(mask))
  units <- data.frame(id = seq_len(n), G = pop$G, x = pop$x, z = z,
    cluster = seq_len(n), sample_weight = 1, norm_weight = 1,
    pattern = pat$pattern[pidx], pattern_type = pat$type[pidx],
    n_obs = as.integer(rowSums(mask)), eligible_ge2 = rowSums(mask) >= 2L)
  panel <- units[id, ]; rownames(panel) <- NULL
  panel$period <- period; panel$observed <- observed
  panel$y <- ifelse(observed, potential$y_g, NA_real_)
  counts <- as.data.frame(table(G = factor(pop$G, s$cohort_codes),
    pattern = factor(units$pattern, pat$pattern)), responseName = "units", stringsAsFactors = FALSE)
  counts$G <- as.numeric(counts$G)
  s$n <- n; s$seed <- seed; s$rng$r_version <- R.version.string
  list(data = list(panel = panel, observed = panel[observed, , drop = FALSE],
      units = units, observation_mask = mask, counts = counts),
    oracle = list(panel = cbind(data.frame(id = id, period = period), potential),
      units = data.frame(id = seq_len(n), unit_effect = unit_effect), truth = cp02_gap_truth(law)),
    settings = s)
}
