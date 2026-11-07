# ==============================================================================
# Runnable Demo: Static Binary Entry, Exit, and Re-Entry (`att_gt(dname = ...)`)
# ==============================================================================
#
# Run with cdid 1.3.0.0 installed in your R library (before merge to main,
# install from the development pull request branch with
# remotes::install_github("dvbn/cdid", ref = "codex/static-entry-exit-1.3.0.0");
# the submitted 1.2.0.1 cohort package on main does not include dname):
#   Rscript --vanilla demos/static_entry_exit.R
#
# User-adjustable parameters are grouped at the top of this script.
# ==============================================================================

# --- 1. User-Adjustable Parameters --------------------------------------------
DEMO_SEED           <- 20261104L
N_UNITS_PER_PATH    <- 12L              # 8 fixed treatment paths x 12 = 96 units
PANEL_DATES         <- c(1, 2, 3, 4)    # Exactly 4 ordered dates for this fixed path catalog
TRUE_TAU_EFFECT     <- 1.25             # True static treatment effect Y_t(1) - Y_t(0)
BOOTSTRAP_ITERS     <- 200L             # Multiplier bootstrap draws for demo speed
CONFIDENCE_ALPHA    <- 0.05             # Significance level (95% intervals / bands)
RUN_OPTIONAL_PLOT   <- TRUE             # Display optional ggdid() plot object if ggplot2 is available

# --- 2. Load the Installed Package --------------------------------------------
library(cdid)
cat("Loaded cdid", as.character(packageVersion("cdid")), "from", find.package("cdid"), "\n")
if (!("dname" %in% names(formals(cdid::att_gt)))) {
  stop("This demo requires an installed cdid >= 1.3.0.0 build with att_gt(dname = ...).")
}
stopifnot(is.numeric(PANEL_DATES), length(PANEL_DATES) == 4L,
          all(is.finite(PANEL_DATES)), all(diff(PANEL_DATES) > 0),
          length(N_UNITS_PER_PATH) == 1L, is.finite(N_UNITS_PER_PATH),
          N_UNITS_PER_PATH >= 12L, N_UNITS_PER_PATH == as.integer(N_UNITS_PER_PATH))
# The fixed four-date observation catalog includes all pair types plus zero-
# and single-outcome units at the default 12 repeats per path. Smaller samples
# are outside this demo's supported parameter range.

cat("==============================================================================\n")
cat("Static Binary Entry, Exit, and Re-Entry Demo (`att_gt(dname = ...)`)\n")
cat("==============================================================================\n\n")

cat("Conceptual Overview:\n")
cat("  1. Source-Static No-Carryover Benchmark:\n")
cat("     Potential outcomes Y_{i,t}(0) and Y_{i,t}(1) depend on current treatment\n")
cat("     status d in {0, 1} without dynamic carryover from earlier treatment history.\n")
cat("     Each transition cell c = (t, d) compares units who switch at date t from\n")
cat("     prior status D_{i,t-1} = d to D_{i,t} = 1 - d against stayers who remain\n")
cat("     at D_{i,t-1} = D_{i,t} = d, conditioning on reference covariates X_{i,t-1}.\n\n")
cat("  2. Sign Orientation:\n")
cat("     With sigma(d) = 1 - 2*d:\n")
cat("       - Entry (d = 0, 0 -> 1, sigma = +1): tau_{t,0} = +theta_{t,0}\n")
cat("       - Exit  (d = 1, 1 -> 0, sigma = -1): tau_{t,1} = -theta_{t,1}\n")
cat("     Both entry and exit effects are reported on the treated-vs-untreated scale\n")
cat("     E[Y_t(1) - Y_t(0) | D_{t-1} = d, D_t = 1 - d].\n\n")
cat("  3. Known Binary Treatment Requirement:\n")
cat("     Binary treatment D_{i,t} in {0, 1} must be observed across the study panel\n")
cat("     for every unit even at dates where outcome Y_{i,t} is missing, because\n")
cat("     switcher (G_{i,c}) and stayer (C_{i,c}) memberships are defined by\n")
cat("     (D_{i,t-1}, D_{i,t}). Units with zero or only one observed outcome are\n")
cat("     retained in the full n-unit study population.\n\n")

# --- 3. Simulate Matched Complete and Missing-Outcome Static Panel ------------
build_demo_static_data <- function(seed, n_rep, dates, tau_true) {
  set.seed(seed)
  # 8 binary treatment histories covering never-treated, always-treated,
  # staggered entry, exit from initial treatment, and entry-exit-reentry (0->1->0->1)
  path_types <- rbind(
    c(0L, 0L, 0L, 0L),  # stayer at 0 across all transitions
    c(0L, 1L, 1L, 1L),  # entry at t = 2, then stayer at 1
    c(0L, 0L, 1L, 1L),  # entry at t = 3, then stayer at 1
    c(0L, 0L, 0L, 1L),  # entry at t = 4
    c(1L, 1L, 1L, 1L),  # stayer at 1 across all transitions
    c(1L, 0L, 0L, 0L),  # exit at t = 2, then stayer at 0
    c(1L, 1L, 0L, 0L),  # exit at t = 3, then stayer at 0
    c(0L, 1L, 0L, 1L)   # entry at t = 2, exit at t = 3, re-entry at t = 4
  )
  # Deterministic observation templates ensuring positive pair support across
  # all 6 date pairs plus explicit single-outcome and zero-outcome units
  obs_templates <- rbind(
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 0L),
    c(1L, 1L, 0L, 1L),
    c(1L, 0L, 1L, 1L),
    c(0L, 1L, 1L, 1L),
    c(1L, 0L, 1L, 0L),
    c(1L, 0L, 0L, 0L),  # single observed outcome
    c(0L, 0L, 0L, 0L)   # zero observed outcomes
  )

  D_mat <- path_types[rep(seq_len(nrow(path_types)), each = n_rep), , drop = FALSE]
  n_units <- nrow(D_mat)
  rows_complete <- vector("list", n_units)
  rows_missing  <- vector("list", n_units)

  for (i in seq_len(n_units)) {
    tpl_idx <- ((i - 1L) %% nrow(obs_templates)) + 1L
    o_vec <- obs_templates[tpl_idx, ]
    x1 <- rnorm(length(dates), sd = 0.35)
    x2 <- runif(length(dates), -0.4, 0.4)
    d_vec <- D_mat[i, ]
    unit_fe <- rnorm(1L, sd = 0.3)
    y_vec <- unit_fe + 0.45 * dates + 0.40 * x1 - 0.25 * x2 + tau_true * d_vec +
      rnorm(length(dates), sd = 0.20)

    cl_id <- sprintf("cluster_%02d", ((i - 1L) %% 12L) + 1L)
    uid <- sprintf("unit_%03d", i)

    rows_complete[[i]] <- data.frame(
      id = uid, time = dates, y = y_vec, d = d_vec,
      x1 = x1, x2 = x2, cluster = cl_id, stringsAsFactors = FALSE
    )
    rows_missing[[i]] <- data.frame(
      id = uid, time = dates, y = ifelse(o_vec == 1L, y_vec, NA_real_), d = d_vec,
      x1 = x1, x2 = x2, cluster = cl_id, stringsAsFactors = FALSE
    )
  }

  list(
    complete = do.call(rbind, rows_complete),
    missing  = do.call(rbind, rows_missing)
  )
}

sim <- build_demo_static_data(
  seed     = DEMO_SEED,
  n_rep    = N_UNITS_PER_PATH,
  dates    = PANEL_DATES,
  tau_true = TRUE_TAU_EFFECT
)

cat(sprintf(
  "Simulated panel: N = %d units across T = %d dates (true static tau = %.2f).\n\n",
  length(unique(sim$complete$id)), length(PANEL_DATES), TRUE_TAU_EFFECT
))

# --- 4. Complete-Data Static Estimation ---------------------------------------
cat("------------------------------------------------------------------------------\n")
cat("Example 1: Complete-Data Static DR Estimation (`selection = \"none\"`)\n")
cat("------------------------------------------------------------------------------\n")
fit_complete_dr <- att_gt(
  yname      = "y",
  tname      = "time",
  idname     = "id",
  dname      = "d",
  xformla    = ~ x1 + x2,
  data       = sim$complete,
  est_method = "dr",
  outcome_pooling = "separate", # Same separate default as the cohort route in 1.2.0.0
  bstrap     = FALSE,
  cband      = FALSE
)
print(summary(fit_complete_dr))
cat("\nComplete-data kernel equivalence is checked in the package tests.\n\n")

# --- 5. Missing-Outcome Static Estimation Across Observation Routes -----------
cat("------------------------------------------------------------------------------\n")
cat("Example 2: Missing-Outcome Doubly Robust Estimation with Clustered Bootstrap\n")
cat("           (`selection = \"conditional_mtar\"`, `selection_est_method = \"dr\"`)\n")
cat("------------------------------------------------------------------------------\n")
fit_missing_aug <- att_gt(
  yname                = "y",
  tname                = "time",
  idname               = "id",
  dname                = "d",
  xformla              = ~ x1 + x2,
  selection_xformla    = ~ x1,
  data                 = sim$missing,
  est_method           = "dr",
  selection            = "conditional_mtar",
  selection_est_method = "dr",
  selection_model      = "pair_logit",
  selection_pooling    = "within_group",
  imputation_pooling   = "within_group",
  imputation_scale     = "per_period",
  outcome_pooling      = stats::setNames(
    c("entry", "entry", "entry", "exit", "exit"),
    c(paste0("t", PANEL_DATES[2:4], "_entry"),
      paste0("t", PANEL_DATES[2:3], "_exit"))),
  clustervars          = "cluster",
  bstrap               = TRUE,
  cband                = TRUE,
  biters               = BOOTSTRAP_ITERS,
  alp                  = CONFIDENCE_ALPHA
)
print(summary(fit_missing_aug))

cat("\nTidy Table (`tidy(fit_missing_aug)`):\n")
print(tidy(fit_missing_aug))

cat("\nGlance Summary (`glance(fit_missing_aug)`):\n")
print(glance(fit_missing_aug))

cat("\nComparison of Aggregate Estimates Across Missing-Outcome Routes:\n")
fit_missing_base <- att_gt(
  yname = "y", tname = "time", idname = "id", dname = "d",
  xformla = ~ x1 + x2, data = sim$missing,
  est_method = "dr", selection = "none",
  bstrap = FALSE, cband = FALSE
)
fit_missing_wt <- att_gt(
  yname = "y", tname = "time", idname = "id", dname = "d",
  xformla = ~ x1 + x2, selection_xformla = ~ x1, data = sim$missing,
  est_method = "ipw", selection = "conditional_mtar", selection_est_method = "ipw",
  selection_pooling = "within_group",
  bstrap = FALSE, cband = FALSE
)
fit_missing_imp <- att_gt(
  yname = "y", tname = "time", idname = "id", dname = "d",
  xformla = ~ x1 + x2, data = sim$missing,
  est_method = "reg", selection = "conditional_mtar", selection_est_method = "imputation",
  imputation_pooling = "within_group", imputation_scale = "per_change",
  outcome_pooling = "pooled",
  bstrap = FALSE, cband = FALSE
)
fit_missing_aug_ipw <- att_gt(
  yname = "y", tname = "time", idname = "id", dname = "d",
  xformla = ~ x1 + x2, selection_xformla = ~ x1, data = sim$missing,
  est_method = "ipw", selection = "conditional_mtar", selection_est_method = "dr",
  selection_model = "pair_logit", selection_pooling = "within_group",
  imputation_pooling = "within_group", imputation_scale = "per_period",
  bstrap = FALSE, cband = FALSE
)
fit_missing_aug_reg <- att_gt(
  yname = "y", tname = "time", idname = "id", dname = "d",
  xformla = ~ x1 + x2, selection_xformla = ~ x1, data = sim$missing,
  est_method = "reg", selection = "conditional_mtar", selection_est_method = "dr",
  selection_model = "pair_logit", selection_pooling = "within_group",
  imputation_pooling = "within_group", imputation_scale = "per_period",
  outcome_pooling = "pooled",
  bstrap = FALSE, cband = FALSE
)

route_comp <- data.frame(
  route      = c(
    "baseline (none / dr)",
    "weighting (ipw / ipw)",
    "imputation (imputation / reg)",
    "augmented (dr / ipw)",
    "augmented (dr / reg)",
    "augmented (dr / dr)"
  ),
  agg_est    = c(
    fit_missing_base$aggregate$estimate,
    fit_missing_wt$aggregate$estimate,
    fit_missing_imp$aggregate$estimate,
    fit_missing_aug_ipw$aggregate$estimate,
    fit_missing_aug_reg$aggregate$estimate,
    fit_missing_aug$aggregate$estimate
  ),
  agg_an_se  = c(
    fit_missing_base$aggregate$analytical_se,
    fit_missing_wt$aggregate$analytical_se,
    fit_missing_imp$aggregate$analytical_se,
    fit_missing_aug_ipw$aggregate$analytical_se,
    fit_missing_aug_reg$aggregate$analytical_se,
    fit_missing_aug$aggregate$analytical_se
  ),
  stringsAsFactors = FALSE
)
print(route_comp, row.names = FALSE)
cat("\n")

# --- 6. Direction / Date Subsets and Explicit Empty-Cell Diagnostics ----------
cat("------------------------------------------------------------------------------\n")
cat("Example 3: Explicit Direction/Date Request Including an Empty Exit Cell\n")
cat("           (`switchers = \"exit\"`, `switch_times = c(3, 4)`)\n")
cat("------------------------------------------------------------------------------\n")
fit_exit_subset <- att_gt(
  yname        = "y",
  tname        = "time",
  idname       = "id",
  dname        = "d",
  xformla      = ~ x1,
  data         = sim$missing,
  switchers    = "exit",
  switch_times = PANEL_DATES[3:4],
  bstrap       = FALSE,
  cband        = FALSE
)
print(summary(fit_exit_subset))

# --- 7. Optional Visualization via ggdid() ------------------------------------
if (isTRUE(RUN_OPTIONAL_PLOT) && requireNamespace("ggplot2", quietly = TRUE)) {
  p_static <- ggdid(fit_missing_aug)
  print(p_static)
  cat(sprintf(
    "\nConstructed optional `ggdid()` plot object (`%s`) distinguished by transition direction.\n",
    paste(class(p_static), collapse = "/")
  ))
}

cat("\nDemo completed successfully.\n")
