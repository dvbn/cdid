# Real-owner and base-R CSV fixture loader for CP07 Chained DR (dr_finite_v1).
# Uses only base-R read.csv and locked installed cdid namespace constructors.
# Never calls any DR/OR/IPW fit or replay method and never modifies or injects
# bindings into the installed package namespace.

dr_fixture_expected_hashes <- c(
  actual_y_values.csv = "fde063b70d13df454b0581bc87885e23e3a986a731a5bf714704957217414d97",
  atoms.csv = "f6869b76895e3c0c91b8d7152bf38f6dbff65ab5400f033c07730c378a03d8d1",
  calendar.csv = "77f699998e69567d62fb55d91e3813f53b39f3cd4ae1fbe9724f66012c2e0039",
  catalog_links.csv = "7c2ebe64591ba45cb471160833fa76477025dbd736874d83f8c58bd4ae00ab4e",
  cohorts.csv = "48b5310372a41bcbe3823488f14a3f52f1f8b3ceade74d7a97ea81482bb16129",
  covariances.csv = "a1a41cd389642acbe26555666be0d2111f0bd779b5f2b5f27038ce798d6f394d",
  fixture_specs.csv = "12e9bf775b9d3c775b03c3d6febdc95f4c82de7878c56895c0b36175fd565708",
  link_summaries.csv = "9e9c17f581148a2f89cbc72e707f43c4ec5a7ab14073843cf2be25a398824048",
  nuisance_if_rows.csv = "a607bf85c4603ef3ee0d276e4cb2e50c41bb0d5edc788b72de450ff44a68ad9b",
  or_blocks.csv = "c6b23b0ec1488ab9af4694dd73e047eb9afd85031b6bda51dde48b4ed75a02ae",
  points.csv = "96163dc9a4c8444d067191bad9f1cfcff473ebfa948c180d28a256b2688e9d0f",
  population_examples.csv = "0378ee9ff9c5dc85794a4680eb272d3705c862779eb31ba9092c91728b27f770",
  qualifications.csv = "99d3d573399124b8c6fcd77342c078f7281560c33afd0469708d5a26d858782d",
  unit_if_rows.csv = "aa5618592cd9bfa345c32cbd1b8af73a286c98e14bcd783778f592ead2356048"
)

dr_fixture_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, envir = ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}

dr_fixture_verify_hashes <- function(directory) {
  if (!requireNamespace("digest", quietly = TRUE)) return(invisible(TRUE))
  for (nm in names(dr_fixture_expected_hashes)) {
    path <- file.path(directory, nm)
    if (!file.exists(path)) stop("Missing packaged DR fixture CSV: ", nm, call. = FALSE)
    got <- digest::digest(file = path, algo = "sha256")
    if (!identical(got, dr_fixture_expected_hashes[[nm]])) {
      stop("SHA-256 mismatch for packaged DR fixture CSV: ", nm, call. = FALSE)
    }
  }
  invisible(TRUE)
}

dr_fixture_csv <- function(directory, name) {
  path <- file.path(directory, name)
  if (!file.exists(path)) stop("Missing packaged DR fixture CSV: ", name, call. = FALSE)
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "NaN", "nan")
  )
}

dr_split_pipe <- function(x) {
  if (is.na(x) || !nzchar(x)) character(0) else strsplit(x, "|", fixed = TRUE)[[1L]]
}

# Construct observed-only outcome matrix (Y_t = NA when R_t == 0) and observed-only
# pair changes (NA when S_e == 0) from atoms.csv, keeping A over all 4 stored dates.
dr_prepare_observed_atoms <- function(atoms, calendar) {
  N <- nrow(atoms)
  M <- nrow(calendar)
  r_cols <- paste0("R_", calendar$time)
  y_cols <- paste0("Y_", calendar$time)
  R_mat <- as.matrix(atoms[r_cols])
  storage.mode(R_mat) <- "double"
  Y_latent <- as.matrix(atoms[y_cols])
  storage.mode(Y_latent) <- "double"

  # Mask latent outcomes where R == 0 so unobserved outcomes are strictly NA_real_
  Y_obs <- Y_latent
  Y_obs[R_mat == 0] <- NA_real_

  # Verify A is defined over ALL 4 stored calendar dates (including units 10, 22, 34)
  n_obs_dates <- rowSums(!is.na(Y_obs))
  A_from_obs <- as.integer(n_obs_dates >= 2L)
  B_from_obs <- as.integer(n_obs_dates >= 1L)
  stopifnot(
    identical(A_from_obs, as.integer(atoms$A)),
    identical(B_from_obs, as.integer(atoms$B)),
    sum(A_from_obs) == 24L,
    all(A_from_obs[atoms$unit_id %in% c(10L, 22L, 34L)] == 1L),
    all(atoms$S_e12[atoms$unit_id %in% c(10L, 22L, 34L)] == 0L),
    all(atoms$S_e23[atoms$unit_id %in% c(10L, 22L, 34L)] == 0L),
    all(atoms$S_e13[atoms$unit_id %in% c(10L, 22L, 34L)] == 0L)
  )

  # Observed-only pair changes: strictly NA_real_ when S_e == 0
  pair_outcomes_delta <- list(
    "2000:2002" = ifelse(atoms$S_e12 == 1L, as.double(atoms$Delta_e12), NA_real_),
    "2002:2005" = ifelse(atoms$S_e23 == 1L, as.double(atoms$Delta_e23), NA_real_),
    "2000:2005" = ifelse(atoms$S_e13 == 1L, as.double(atoms$Delta_e13), NA_real_)
  )
  pair_outcomes_y <- list(
    "2000:2002" = ifelse(atoms$S_e12 == 1L, Y_obs[, 2L] - Y_obs[, 1L], NA_real_),
    "2002:2005" = ifelse(atoms$S_e23 == 1L, Y_obs[, 3L] - Y_obs[, 2L], NA_real_),
    "2000:2005" = ifelse(atoms$S_e13 == 1L, Y_obs[, 3L] - Y_obs[, 1L], NA_real_)
  )

  list(
    R = R_mat,
    Y_obs = Y_obs,
    A = A_from_obs,
    B = B_from_obs,
    pair_outcomes_delta = pair_outcomes_delta,
    pair_outcomes_y = pair_outcomes_y
  )
}

# Parse expectations for one fixture_id aligned to target unit_ids
dr_load_expectations <- function(directory, fixture_id, unit_ids, link_ids) {
  all_points <- dr_fixture_csv(directory, "points.csv")
  pts <- all_points[all_points$fixture_id == fixture_id, , drop = FALSE]
  point_group <- function(grp) {
    sub <- pts[pts$group == grp, , drop = FALSE]
    if (!nrow(sub)) return( stats::setNames(numeric(0), character(0)) )
    ord <- order(as.integer(sub$index))
    stats::setNames(as.double(sub$value[ord]), sub$name[ord])
  }

  all_cov <- dr_fixture_csv(directory, "covariances.csv")
  cov_sub <- all_cov[all_cov$fixture_id == fixture_id, , drop = FALSE]
  cov_matrix <- function(mtype) {
    sub <- cov_sub[cov_sub$matrix_type == mtype, , drop = FALSE]
    L <- length(link_ids)
    mat <- matrix(NA_real_, L, L, dimnames = list(link_ids, link_ids))
    for (r in seq_len(nrow(sub))) {
      i <- match(sub$row_link[r], link_ids)
      j <- match(sub$col_link[r], link_ids)
      stopifnot(!anyNA(c(i, j)))
      mat[i, j] <- as.double(sub$value[r])
    }
    stopifnot(!anyNA(mat))
    mat
  }

  all_uif <- dr_fixture_csv(directory, "unit_if_rows.csv")
  uif <- all_uif[all_uif$fixture_id == fixture_id, , drop = FALSE]
  N <- length(unit_ids)
  L <- length(link_ids)
  uif_matrix <- function(col) {
    mat <- matrix(NA_real_, N, L, dimnames = list(NULL, link_ids))
    for (j in seq_len(L)) {
      lid <- link_ids[j]
      sub <- uif[uif$link_id == lid, , drop = FALSE]
      idx <- match(unit_ids, sub$unit_id)
      stopifnot(!anyNA(idx), identical(as.integer(sub$unit_id[idx]), as.integer(unit_ids)))
      mat[, j] <- as.double(sub[[col]][idx])
    }
    mat
  }

  all_nif <- dr_fixture_csv(directory, "nuisance_if_rows.csv")
  nif <- all_nif[all_nif$fixture_id == fixture_id, , drop = FALSE]
  nif_matrix <- function(ntype, pnames) {
    if (!length(pnames)) return(matrix(0, N, 0L))
    mat <- matrix(NA_real_, N, length(pnames), dimnames = list(NULL, pnames))
    sub_t <- nif[nif$nuisance_type == ntype, , drop = FALSE]
    for (j in seq_along(pnames)) {
      sub <- sub_t[sub_t$param_name == pnames[j], , drop = FALSE]
      idx <- match(unit_ids, sub$unit_id)
      stopifnot(!anyNA(idx), identical(as.integer(sub$unit_id[idx]), as.integer(unit_ids)))
      mat[, j] <- as.double(sub$analytic_if[idx])
    }
    mat
  }

  all_lsum <- dr_fixture_csv(directory, "link_summaries.csv")
  lsum <- all_lsum[all_lsum$fixture_id == fixture_id, , drop = FALSE]
  lsum <- lsum[match(link_ids, lsum$link_id), , drop = FALSE]

  list(
    points = pts,
    beta_hat = point_group("beta_hat"),
    alpha_hat = point_group("alpha_hat"),
    mu_dr = point_group("mu_dr"),
    mu_or = point_group("mu_or"),
    mu_ipw = point_group("mu_ipw"),
    theta_gmm_fixed_W = point_group("theta_gmm_fixed_W"),
    V_unit = cov_matrix("V_unit"),
    V_cluster = cov_matrix("V_cluster"),
    Omega_cluster = cov_matrix("Omega_cluster"),
    psi_dr = uif_matrix("psi_dr"),
    psi_dr_core = uif_matrix("psi_dr_core"),
    psi_dir = uif_matrix("psi_dir"),
    psi_or = uif_matrix("psi_or"),
    psi_ps = uif_matrix("psi_ps"),
    psi_sel = uif_matrix("psi_sel"),
    psi_ipw = uif_matrix("psi_ipw"),
    fd_1e6_quad = uif_matrix("fd_1e6_quad"),
    fd_5e5_quad = uif_matrix("fd_5e5_quad"),
    fd_1e5_float64 = uif_matrix("fd_1e5_float64"),
    abs_err_1e6 = uif_matrix("abs_err_1e6"),
    tau_if = nif_matrix("tau", c("tau_0", "tau_2005", "tau_2010")),
    eta_if = nif_matrix("eta", c("eta_0", "eta_2005", "eta_2010")),
    beta_if_fn = function(or_slopes) {
      bnames <- if (length(or_slopes)) paste0("beta_", or_slopes) else character(0)
      nif_matrix("beta", bnames)
    },
    link_summaries = lsum
  )
}

# Build real locked installed cdid owners (panel, frame, reader, link_views)
# and aligned reference expectations for any of the 8 fixtures F1..F8.
load_dr_finite_reference <- function(directory,
                                     fixture_id = "F1_baseline_none_intercept_slope",
                                     order = NULL,
                                     max_chunk_rows = 13L) {
  dr_fixture_verify_hashes(directory)
  atoms <- dr_fixture_csv(directory, "atoms.csv")
  atoms$cohort <- as.double(atoms$cohort)
  N <- nrow(atoms)
  if (is.null(order)) order <- seq_len(N)
  stopifnot(
    is.numeric(order),
    length(order) == N,
    identical(sort(as.integer(order)), seq_len(N))
  )

  specs <- dr_fixture_csv(directory, "fixture_specs.csv")
  meta <- specs[specs$fixture_id == fixture_id, , drop = FALSE]
  stopifnot(nrow(meta) == 1L)

  calendar <- dr_fixture_csv(directory, "calendar.csv")
  cohorts <- dr_fixture_csv(directory, "cohorts.csv")
  links <- dr_fixture_csv(directory, "catalog_links.csv")
  blocks <- dr_fixture_csv(directory, "or_blocks.csv")
  obs_prep <- dr_prepare_observed_atoms(atoms, calendar)

  control <- dr_fixture_ns("ipw_owner_control")(
    max_chunk_rows = as.integer(max_chunk_rows),
    max_retained_entries = 1e7,
    max_matrix_entries = 2e7,
    max_work_ops = 1e10,
    max_output_entries = 1e7
  )

  # Construct long panel with Y = NA_real_ whenever R == 0.
  # Keep one explicit all-NA zero-observer unit in the panel; the other 5 zero-observer
  # units enter through the linked frame with allow_frame_only = TRUE.
  Y_obs <- obs_prep$Y_obs
  long <- data.frame(
    id = rep(as.double(atoms$unit_id), each = nrow(calendar)),
    time = rep(as.double(calendar$time), N),
    y = as.vector(t(Y_obs)),
    g = rep(as.double(atoms$cohort), each = nrow(calendar))
  )
  zero <- which(rowSums(!is.na(Y_obs)) == 0L)
  if (length(zero) > 1L) {
    long <- long[!long$id %in% atoms$unit_id[zero[-1L]], , drop = FALSE]
  }
  panel <- dr_fixture_ns("construct_ipw_panel_owner")(
    long, "id", "time", "y", "g", control
  )
  stopifnot(
    panel$n_a_units == 24L,
    all(c(10, 22, 34) %in% panel$a_id)
  )

  or_slopes <- dr_split_pipe(meta$or_slopes)
  or_intercept <- as.integer(meta$or_intercept) == 1L
  trt_raw <- dr_split_pipe(meta$trt_cols)
  trt_cols <- ifelse(trt_raw == "1", "(Intercept)", trt_raw)
  or_cols <- c(if (or_intercept) "(Intercept)", or_slopes)
  frame_cols <- unique(c(if (or_intercept || "(Intercept)" %in% trt_cols) "(Intercept)", or_slopes, setdiff(trt_cols, "(Intercept)")))

  all_X <- cbind(`(Intercept)` = rep(1, N), x1 = as.double(atoms$x1), x2 = as.double(atoms$x2))
  Z_frame <- all_X[, frame_cols, drop = FALSE]
  # Set Z_frame to NA off A so design_on_A_only = TRUE is strictly enforced
  Z_frame[obs_prep$A == 0L, ] <- NA_real_

  frame <- dr_fixture_ns("construct_ipw_frame_owner")(
    panel,
    as.double(atoms$unit_id[order]),
    Z_frame[order, , drop = FALSE],
    weights = as.double(atoms$w[order]),
    control = control,
    frame_g = as.double(atoms$cohort[order]),
    allow_frame_only = TRUE,
    design_on_A_only = TRUE
  )

  # Selection design Q and universe U / transport rho
  sel_slopes <- dr_split_pipe(meta$sel_slopes)
  sel_intercept <- as.integer(meta$sel_intercept) == 1L
  if (meta$route == "none") {
    Q <- cbind(`(Intercept)` = rep(1, N))
    has_icpt_decl <- NULL
  } else {
    sel_cols <- c(if (sel_intercept) "(Intercept)", sel_slopes)
    Q <- all_X[, sel_cols, drop = FALSE]
    has_icpt_decl <- sel_intercept
  }
  U_vec <- as.logical(atoms[[meta$U]] == 1L)
  rho_vec <- as.double(atoms[[paste0("rho_", meta$U)]])
  Q[!U_vec, ] <- NA_real_
  rho_vec[!U_vec] <- NA_real_
  just <- if (meta$U %in% c("B", "A")) {
    "known transport under the frozen independent fixed atom law"
  } else {
    NULL
  }

  reader <- dr_fixture_ns("construct_cp04_frame_reader")(
    frame,
    as.double(atoms$unit_id[order]),
    Q[order, , drop = FALSE],
    as.double(atoms$cohort[order]),
    as.double(atoms$w[order]),
    control,
    U = U_vec[order],
    rho = rho_vec[order],
    universe = meta$U,
    justification = just,
    has_intercept = has_icpt_decl
  )

  views <- lapply(seq_len(nrow(links)), function(i) {
    lk <- links[i, , drop = FALSE]
    K <- as.double(dr_split_pipe(lk$K))
    dr_fixture_ns("construct_ipw_link_view")(
      frame,
      as.double(lk$s),
      as.double(lk$t),
      g = as.double(lk$g),
      K = K,
      control = control,
      link_id = lk$link_id,
      functional_key = lk$functional_key
    )
  })
  names(views) <- links$link_id

  outcome_recipe <- dr_fixture_ns("or_pooled_outcome_recipe")(
    or_cols, or_intercept, or_slopes
  )

  ids <- as.integer(reader$state$id)
  ord_idx <- match(ids, atoms$unit_id)
  exp <- dr_load_expectations(directory, fixture_id, ids, links$link_id)
  exp$beta_if <- exp$beta_if_fn(or_slopes)

  list(
    fixture_id = fixture_id,
    spec = meta,
    panel = panel,
    frame = frame,
    reader = reader,
    views = views,
    control = control,
    glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE),
    rank_tol = 1e-10,
    outcome_recipe = outcome_recipe,
    treatment_columns = trt_cols,
    selection = if (meta$route == "none") "none" else "conditional_mtar",
    selection_model = meta$selection_model,
    selection_pooling = meta$selection_pooling,
    has_intercept = has_icpt_decl,
    raw = atoms[ord_idx, , drop = FALSE],
    R = obs_prep$R[ord_idx, , drop = FALSE],
    Y_obs = obs_prep$Y_obs[ord_idx, , drop = FALSE],
    pair_outcomes_delta = lapply(obs_prep$pair_outcomes_delta, `[`, ord_idx),
    pair_outcomes_y = lapply(obs_prep$pair_outcomes_y, `[`, ord_idx),
    Q = Q[ord_idx, , drop = FALSE],
    X_or = all_X[ord_idx, or_cols, drop = FALSE],
    V_trt = all_X[ord_idx, trt_cols, drop = FALSE],
    U = U_vec[ord_idx],
    rho = rho_vec[ord_idx],
    catalog = list(
      calendar = calendar,
      cohorts = cohorts,
      links = links,
      blocks = blocks,
      link_list = lapply(seq_len(nrow(links)), function(i) {
        list(
          link_id = links$link_id[i],
          g = as.double(links$g[i]),
          K = as.double(dr_split_pipe(links$K[i])),
          s = as.double(links$s[i]),
          t = as.double(links$t[i]),
          or_block = links$or_block[i]
        )
      })
    ),
    expectations = exp
  )
}
