# Transparent synthetic panel generator for notebooks/cdid_practical_guide.qmd.
# Base R only; no package wrapper calls or hidden estimation helpers.

default_dgp_config <- function(...) {
  cfg <- list(
    n_units           = 2000L,                  # Study population size N in `units` (frame)
    n_periods         = 6L,                     # Calendar dates t = 1, ..., T
    seed              = 20261029L,              # RNG seed for drawing the synthetic panel
    cohorts           = c(3L, 4L, 5L),          # Positive treatment cohorts (0 = never treated)
    cohort_intercepts = c(-0.35, -0.45, -0.40), # Baseline log-odds vs G = 0 when x = 0
    treat_x_slope     = c(0.55, 0.35, -0.25),   # Log-odds shift per cohort when binary x = 1
    prob_x            = 0.50,                   # P(x = 1) for binary unit covariate x
    n_clusters        = 100L,                   # Number of regional clusters in `units`
    cluster_trend_sd  = 0.00,                   # SD of mean-zero cluster-level trend slope (0 = IID default)
    unit_sd           = 1.00,                   # SD of time-invariant unit fixed effect alpha_i
    unit_x_shift      = 0.40,                   # Shift in alpha_i when x = 1 (cancels in changes)
    time_trend        = 0.20,                   # Common per-period untreated trend
    trend_x           = 0.30,                   # Per-period untreated trend slope on binary x
    trend_z           = 0.15,                   # Per-period untreated trend slope on normal z
    effect_level      = 1.00,                   # Treatment effect level at onset e = t - g = 0
    effect_slope      = 0.25,                   # Treatment effect slope per event-time step e
    effect_x          = 0.50,                   # Treatment effect heterogeneity on binary x
    noise_sd          = 0.50,                   # SD of independent normal shock eps_{i,t}
    obs_intercept     = 0.25,                   # Baseline date-observation log-odds at t = 1, x = 0
    obs_time_slope    = 0.08,                   # Per-period shift in date-observation log-odds
    obs_x_slope       = 0.55,                   # Shift in date-observation log-odds when x = 1
    obs_z_slope       = 0.00,                   # Shift in date-observation log-odds on normal z
    obs_treat_shift   = 0.00                    # Shift in date-observation log-odds for G > 0
  )
  overrides <- list(...)
  if (length(overrides) > 0L) {
    unknown <- setdiff(names(overrides), names(cfg))
    if (length(unknown) > 0L) stop("Unknown DGP config field(s): ", paste(unknown, collapse = ", "))
    cfg[names(overrides)] <- overrides
  }
  cfg
}

simulate_panel <- function(cfg = default_dgp_config()) {
  stopifnot(
    is.list(cfg),
    length(cfg$n_units) == 1L, !is.na(cfg$n_units), cfg$n_units >= 1L, cfg$n_units == as.integer(cfg$n_units),
    length(cfg$n_periods) == 1L, !is.na(cfg$n_periods), cfg$n_periods >= 2L, cfg$n_periods == as.integer(cfg$n_periods),
    length(cfg$n_clusters) == 1L, !is.na(cfg$n_clusters), cfg$n_clusters >= 1L, cfg$n_clusters == as.integer(cfg$n_clusters),
    length(cfg$cohorts) >= 1L, all(!is.na(cfg$cohorts)), all(cfg$cohorts == as.integer(cfg$cohorts)),
    all(cfg$cohorts >= 1L & cfg$cohorts <= cfg$n_periods), !anyDuplicated(cfg$cohorts),
    cfg$prob_x > 0, cfg$prob_x < 1, cfg$noise_sd >= 0, cfg$unit_sd >= 0, cfg$cluster_trend_sd >= 0,
    length(cfg$cohort_intercepts) == length(cfg$cohorts),
    length(cfg$treat_x_slope) == length(cfg$cohorts)
  )
  set.seed(as.integer(cfg$seed))
  N <- as.integer(cfg$n_units); T_dates <- as.integer(cfg$n_periods)
  pos_g <- as.integer(cfg$cohorts); all_g <- c(0L, pos_g)

  # Exact conditional cohort probabilities P(G = g | X = x) for x in {0, 1}
  odds_x0 <- c(1, exp(cfg$cohort_intercepts))
  odds_x1 <- c(1, exp(cfg$cohort_intercepts + cfg$treat_x_slope))
  prob_g_x0 <- stats::setNames(odds_x0 / sum(odds_x0), as.character(all_g))
  prob_g_x1 <- stats::setNames(odds_x1 / sum(odds_x1), as.character(all_g))
  pop_prob_g <- (1 - cfg$prob_x) * prob_g_x0 + cfg$prob_x * prob_g_x1
  pop_ex_given_g <- (cfg$prob_x * prob_g_x1) / pop_prob_g

  # Draw unit table `units` (full study population F)
  x <- stats::rbinom(N, size = 1L, prob = cfg$prob_x)
  z <- stats::rnorm(N, mean = 0, sd = 1)
  cluster <- rep(seq_len(as.integer(cfg$n_clusters)), length.out = N)
  u_draw <- stats::runif(N)
  cdf_x0 <- cumsum(prob_g_x0); cdf_x1 <- cumsum(prob_g_x1)
  g <- vapply(seq_len(N), function(i) {
    all_g[which(u_draw[i] <= (if (x[i] == 1L) cdf_x1 else cdf_x0))[1L]]
  }, integer(1L))
  units <- data.frame(id = seq_len(N), g = g, x = x, z = z, cluster = cluster)

  # Latent complete outcomes on the N x T panel grid
  alpha_i <- cfg$unit_x_shift * x + stats::rnorm(N, mean = 0, sd = cfg$unit_sd)
  cluster_slope <- stats::rnorm(as.integer(cfg$n_clusters), mean = 0, sd = cfg$cluster_trend_sd)
  id_rep <- rep(units$id, each = T_dates); time_rep <- rep(seq_len(T_dates), times = N)
  base_cols <- data.frame(
    id = id_rep, time = time_rep, g = units$g[id_rep],
    x = units$x[id_rep], z = units$z[id_rep], cluster = units$cluster[id_rep]
  )
  eps_it <- stats::rnorm(N * T_dates, mean = 0, sd = cfg$noise_sd)
  y0 <- alpha_i[id_rep] +
    (cfg$time_trend + cfg$trend_x * base_cols$x + cfg$trend_z * base_cols$z + cluster_slope[base_cols$cluster]) * time_rep +
    eps_it
  treated_now <- base_cols$g > 0L & time_rep >= base_cols$g
  tau <- ifelse(treated_now, cfg$effect_level + cfg$effect_slope * (time_rep - base_cols$g) + cfg$effect_x * base_cols$x, 0)
  y_complete <- y0 + tau

  balanced_panel <- base_cols; balanced_panel$y <- y_complete

  # Mask 1: independent-date Bernoulli observation conditional on (x, z, G, t)
  eta_obs <- cfg$obs_intercept + cfg$obs_time_slope * (time_rep - 1L) +
    cfg$obs_x_slope * base_cols$x + cfg$obs_z_slope * base_cols$z + cfg$obs_treat_shift * (base_cols$g > 0L)
  observed_indep <- stats::runif(N * T_dates) < stats::plogis(eta_obs)
  missing_panel <- base_cols; missing_panel$y <- ifelse(observed_indep, y_complete, NA_real_)

  # Mask 2: rotating consecutive two-period window {w, w + 1} on the same latent outcomes
  window_start <- sample.int(T_dates - 1L, size = N, replace = TRUE)
  in_window <- (time_rep == window_start[id_rep]) | (time_rep == window_start[id_rep] + 1L)
  rotating_panel <- base_cols; rotating_panel$y <- ifelse(in_window, y_complete, NA_real_)

  # Analytic superpopulation ATT truth and finite-sample latent mean effects
  roster_x_mean <- stats::setNames(vapply(pos_g, function(gv) mean(units$x[units$g == gv]), numeric(1L)), as.character(pos_g))
  roster_g_share <- stats::setNames(vapply(pos_g, function(gv) mean(units$g == gv), numeric(1L)), as.character(pos_g))
  pop_g_share <- pop_prob_g[as.character(pos_g)]

  cells_truth <- do.call(rbind, lapply(pos_g, function(gv) {
    t_post <- seq.int(gv, T_dates); ch <- as.character(gv)
    data.frame(
      target.id          = sprintf("ATT(g=%d,t=%d)", gv, t_post),
      group = gv, time = t_post, event_time = t_post - gv,
      population_att     = cfg$effect_level + cfg$effect_slope * (t_post - gv) + cfg$effect_x * pop_ex_given_g[ch],
      roster_latent_mean = cfg$effect_level + cfg$effect_slope * (t_post - gv) + cfg$effect_x * roster_x_mean[ch],
      pop_cohort_share   = unname(pop_g_share[ch]),
      roster_share       = unname(roster_g_share[ch]),
      stringsAsFactors   = FALSE
    )
  }))
  rownames(cells_truth) <- NULL

  w_avg <- function(vals, w) sum(vals * w) / sum(w)
  dyn_e <- sort(unique(cells_truth$event_time))
  dynamic_truth <- data.frame(
    event_time         = dyn_e,
    population_att     = vapply(dyn_e, function(e) w_avg(cells_truth$population_att[cells_truth$event_time == e], cells_truth$pop_cohort_share[cells_truth$event_time == e]), numeric(1L)),
    roster_latent_mean = vapply(dyn_e, function(e) w_avg(cells_truth$roster_latent_mean[cells_truth$event_time == e], cells_truth$roster_share[cells_truth$event_time == e]), numeric(1L))
  )
  group_truth <- data.frame(
    group              = pos_g,
    population_att     = vapply(pos_g, function(gv) mean(cells_truth$population_att[cells_truth$group == gv]), numeric(1L)),
    roster_latent_mean = vapply(pos_g, function(gv) mean(cells_truth$roster_latent_mean[cells_truth$group == gv]), numeric(1L))
  )
  cal_t <- sort(unique(cells_truth$time))
  calendar_truth <- data.frame(
    time               = cal_t,
    population_att     = vapply(cal_t, function(tt) w_avg(cells_truth$population_att[cells_truth$time == tt], cells_truth$pop_cohort_share[cells_truth$time == tt]), numeric(1L)),
    roster_latent_mean = vapply(cal_t, function(tt) w_avg(cells_truth$roster_latent_mean[cells_truth$time == tt], cells_truth$roster_share[cells_truth$time == tt]), numeric(1L))
  )
  overall_truth <- data.frame(
    type = c("simple", "group", "dynamic", "calendar"),
    population_att = c(w_avg(cells_truth$population_att, cells_truth$pop_cohort_share), w_avg(group_truth$population_att, pop_g_share), mean(dynamic_truth$population_att), mean(calendar_truth$population_att)),
    roster_latent_mean = c(w_avg(cells_truth$roster_latent_mean, cells_truth$roster_share), w_avg(group_truth$roster_latent_mean, roster_g_share), mean(dynamic_truth$roster_latent_mean), mean(calendar_truth$roster_latent_mean)),
    stringsAsFactors = FALSE
  )

  list(
    config = cfg, units = units, balanced_panel = balanced_panel,
    missing_panel = missing_panel, rotating_panel = rotating_panel,
    latent = list(y_complete = y_complete, y0 = y0, tau = tau, q_indep = stats::plogis(eta_obs), observed_indep = observed_indep, window_start = window_start),
    truth = list(cells = cells_truth, dynamic = dynamic_truth, group = group_truth, calendar = calendar_truth, overall = overall_truth,
                 pop_prob_g = pop_prob_g, pop_ex_given_g = pop_ex_given_g, roster_g_share = roster_g_share, roster_x_mean = roster_x_mean)
  )
}
