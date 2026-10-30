# Real-owner and base-R CSV fixture loader for D52 Chained Pure Imputation
# (imputation_finite_v1).
# Uses only base-R read.csv and locked installed cdid namespace constructors.
# Masks all Y_t by R_t BEFORE constructing the panel so unobserved outcome levels
# are strictly NA_real_ and never used.

imp_fixture_expected_hashes <- c(
  atoms.csv = "bd67374134a87db27fb229d20528268c4146b5eda340001aa3bf9a94b1cc9b30",
  source_observed_panel_atoms.csv = "f6869b76895e3c0c91b8d7152bf38f6dbff65ab5400f033c07730c378a03d8d1",
  points.csv = "0e39b01a20ba3a2c82f8635049e5ef264fd83ae12bfb13b6c7083ab49e2ff930",
  nuisance.csv = "3676de0c1666952999aeed5bcc241cb912ba29ac54b1be92def62bed63e9b32c",
  unit_if.csv = "b26936adaea024d1a1b3209c069624dc86800d2adc6f2f12a2d1100c7b3ec655",
  components.csv = "2ed9ef42254870010eb86ffd09b8a4627d16c157cc8bae69cdfd32d434dd4302",
  covariance.csv = "5cc745cf9eafae60d772bc5be6ad7424b7c5b048009b8b6df67fb7829024c2b5",
  predictions.csv = "de44fc5e355e47d8d2101324f2d9c35afcf6edb9e10581d47200d44fcba6eb93",
  projection.csv = "35f06e36f459380610c763036f2c58e65184b413cbdc94cd16000bb11cc36369",
  projection_if.csv = "d97bed4fd7249e730826d053fdc9fdd3940bff46a9329acfa99f260006bc2bea",
  finite_contract.json = "2e5c74d7accedacef81472ec1f56fa46c00e3e749f5a60aaed62117d116b5955",
  check_results.json = "271d770687a13ea6ec8c261d6546c2c9c81a63bbd2133fbdc66f9af5421822cc"
)

imp_fixture_specs <- data.frame(
  fixture_id = c(
    "F1_sep_int_x1x2",
    "F2_sep_noint_x1x2",
    "F3_sep_int_only",
    "F4_within_int_x1x2",
    "F5_within_noint_x1x2",
    "F6_within_int_only",
    "F7_across_int_x1x2",
    "F8_across_noint_x1x2",
    "F9_across_int_only"
  ),
  imputation_pooling = c(
    rep("separate", 3L),
    rep("within_cohort", 3L),
    rep("across_cohorts", 3L)
  ),
  recipe_key = rep(c("intercept_x1_x2", "no_intercept", "intercept_only"), 3L),
  or_intercept = rep(c(1L, 0L, 1L), 3L),
  or_slopes = rep(c("x1|x2", "x1|x2", ""), 3L),
  trt_cols = rep(c("(Intercept)|x1|x2", "x1|x2", "(Intercept)"), 3L),
  stringsAsFactors = FALSE
)

imp_catalog_links <- data.frame(
  link_id = c("l1", "l2", "l3", "l4"),
  g = c(2005, 2010, 2005, 2005),
  K = c("0", "0", "0|2010", "0|2010"),
  s = c(2000, 2000, 2002, 2000),
  t = c(2002, 2002, 2005, 2005),
  L = c(1L, 1L, 1L, 2L),
  pair_name = c("e12", "e12", "e23", "e13"),
  or_block = c("0", "0", "1", "2"),
  functional_key = c(
    "cmpcat:IMP:l1:g2005:K0:2000->2002",
    "cmpcat:IMP:l2:g2010:K0:2000->2002",
    "cmpcat:IMP:l3:g2005:K0,2010:2002->2005",
    "cmpcat:IMP:l4:g2005:K0,2010:2000->2005"
  ),
  stringsAsFactors = FALSE
)

imp_calendar <- data.frame(
  time_index = 1:4,
  time = c(2000, 2002, 2005, 2010),
  stringsAsFactors = FALSE
)

imp_fixture_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, envir = ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}

imp_fixture_resolve_path <- function(directory, name) {
  path <- file.path(directory, name)
  if (file.exists(path)) return(path)
  if (identical(name, "source_observed_panel_atoms.csv")) {
    alt1 <- file.path(dirname(directory), name)
    if (file.exists(alt1)) return(alt1)
    alt2 <- file.path(dirname(directory), "dr_finite_v1", "atoms.csv")
    if (file.exists(alt2)) return(alt2)
  }
  stop("Missing packaged imputation fixture file: ", name, call. = FALSE)
}

imp_fixture_verify_hashes <- function(directory) {
  if (!requireNamespace("digest", quietly = TRUE)) return(invisible(TRUE))
  for (nm in names(imp_fixture_expected_hashes)) {
    path <- imp_fixture_resolve_path(directory, nm)
    got <- digest::digest(file = path, algo = "sha256")
    if (!identical(got, imp_fixture_expected_hashes[[nm]])) {
      stop("SHA-256 mismatch for packaged imputation fixture file: ", nm, call. = FALSE)
    }
  }
  invisible(TRUE)
}

imp_fixture_csv <- function(directory, name) {
  path <- imp_fixture_resolve_path(directory, name)
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "NaN", "nan")
  )
}

imp_split_pipe <- function(x) {
  if (is.na(x) || !nzchar(x)) character(0) else strsplit(x, "|", fixed = TRUE)[[1L]]
}

# Construct observed-only outcome matrix (Y_t = NA when R_t == 0) BEFORE panel
# construction, and observed-only pair changes (NA when S_e == 0) from atoms.csv.
imp_prepare_observed_atoms <- function(atoms_ref, atoms_panel, calendar = imp_calendar) {
  N <- nrow(atoms_ref)
  stopifnot(
    nrow(atoms_panel) == N,
    identical(as.integer(atoms_ref$unit_id), as.integer(atoms_panel$unit_id)),
    identical(as.integer(atoms_ref$cluster_id), as.integer(atoms_panel$cluster_id)),
    identical(as.integer(atoms_ref$cohort), as.integer(atoms_panel$cohort))
  )
  r_cols <- paste0("R_", calendar$time)
  y_cols <- paste0("Y_", calendar$time)
  R_mat <- as.matrix(atoms_ref[r_cols])
  storage.mode(R_mat) <- "double"
  R_panel <- as.matrix(atoms_panel[r_cols])
  storage.mode(R_panel) <- "double"
  stopifnot(identical(R_mat, R_panel))

  # Mask all Y_t by R_t immediately so unobserved outcomes are strictly NA_real_
  Y_raw <- as.matrix(atoms_panel[y_cols])
  storage.mode(Y_raw) <- "double"
  Y_obs <- ifelse(R_mat == 1, Y_raw, NA_real_)
  storage.mode(Y_obs) <- "double"
  colnames(Y_obs) <- y_cols
  stopifnot(all(is.na(Y_obs[R_mat == 0])))

  # Verify A is defined over ALL 4 stored calendar dates (including units 10, 22, 34)
  n_obs_dates <- rowSums(!is.na(Y_obs))
  A_from_obs <- as.integer(n_obs_dates >= 2L)
  B_from_obs <- as.integer(n_obs_dates >= 1L)
  stopifnot(
    identical(A_from_obs, as.integer(atoms_panel$A)),
    identical(B_from_obs, as.integer(atoms_panel$B)),
    sum(A_from_obs) == 24L,
    all(A_from_obs[atoms_ref$unit_id %in% c(10L, 22L, 34L)] == 1L),
    all(atoms_ref$S_e12[atoms_ref$unit_id %in% c(10L, 22L, 34L)] == 0L),
    all(atoms_ref$S_e23[atoms_ref$unit_id %in% c(10L, 22L, 34L)] == 0L),
    all(atoms_ref$S_e13[atoms_ref$unit_id %in% c(10L, 22L, 34L)] == 0L),
    all(is.na(atoms_ref$Delta_e12[atoms_ref$S_e12 == 0L])),
    all(is.na(atoms_ref$Delta_e23[atoms_ref$S_e23 == 0L])),
    all(is.na(atoms_ref$Delta_e13[atoms_ref$S_e13 == 0L]))
  )

  pair_outcomes_delta <- list(
    "2000:2002" = ifelse(atoms_ref$S_e12 == 1L, as.double(atoms_ref$Delta_e12), NA_real_),
    "2002:2005" = ifelse(atoms_ref$S_e23 == 1L, as.double(atoms_ref$Delta_e23), NA_real_),
    "2000:2005" = ifelse(atoms_ref$S_e13 == 1L, as.double(atoms_ref$Delta_e13), NA_real_)
  )
  pair_outcomes_y <- list(
    "2000:2002" = ifelse(atoms_ref$S_e12 == 1L, Y_obs[, 2L] - Y_obs[, 1L], NA_real_),
    "2002:2005" = ifelse(atoms_ref$S_e23 == 1L, Y_obs[, 3L] - Y_obs[, 2L], NA_real_),
    "2000:2005" = ifelse(atoms_ref$S_e13 == 1L, Y_obs[, 3L] - Y_obs[, 1L], NA_real_)
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

# Load reference tables for one fixture_id aligned to target unit_ids
imp_load_expectations <- function(directory, fixture_id, unit_ids, link_ids = imp_catalog_links$link_id) {
  N <- length(unit_ids)
  L <- length(link_ids)
  meta <- imp_fixture_specs[imp_fixture_specs$fixture_id == fixture_id, , drop = FALSE]
  reference_id <- paste0(meta$imputation_pooling, "__", switch(meta$recipe_key,
    intercept_x1_x2 = "intercept_slopes", no_intercept = "no_intercept",
    intercept_only = "intercept_only"))
  read_fixture <- function(name) {
    tab <- imp_fixture_csv(directory, name)
    ans <- tab[tab$fixture == reference_id, , drop = FALSE]
    stopifnot(nrow(ans) > 0L)
    ans
  }
  pts <- read_fixture("points.csv")
  nuis <- read_fixture("nuisance.csv")
  nuis_vec <- stats::setNames(as.double(nuis$value), nuis$coordinate)
  uif <- read_fixture("unit_if.csv")
  uif_vec <- function(nm, col = "psi") {
    sub <- uif[uif$coordinate == nm, , drop = FALSE]
    idx <- match(unit_ids, sub$unit_id)
    stopifnot(nrow(sub) == N, !anyNA(idx), !anyDuplicated(sub$unit_id))
    as.double(sub[[col]][idx])
  }
  uif_method_mat <- function(method, col = "analytic") {
    field <- switch(col, analytic = "psi", fd_central_1e6 = "symmetric_FD", col)
    vapply(paste0(method, ".", link_ids), uif_vec, numeric(N), col = field)
  }
  comp <- read_fixture("components.csv")
  channels <- c("direct", "a_direct", "j_direct", "treatment_direct", "treatment_a", "OR_direct", "OR_a", "OR_j")
  comp_method_link <- function(method, lid) {
    sub <- comp[comp$coordinate == paste0(method, ".", lid), , drop = FALSE]
    idx <- match(unit_ids, sub$unit_id)
    stopifnot(nrow(sub) == N, !anyNA(idx))
    mat <- as.matrix(sub[idx, channels, drop = FALSE])
    storage.mode(mat) <- "double"
    mat
  }
  cov_sub <- read_fixture("covariance.csv")
  joint_coords <- unique(cov_sub$row_coordinate)
  cov_matrix <- function(mtype, method = NULL) {
    sub <- cov_sub[cov_sub$scale == mtype, , drop = FALSE]
    coords <- if (is.null(method)) joint_coords else paste0(method, ".", link_ids)
    mat <- matrix(NA_real_, length(coords), length(coords), dimnames = list(coords, coords))
    sub <- sub[sub$row_coordinate %in% coords & sub$col_coordinate %in% coords, , drop = FALSE]
    stopifnot(nrow(sub) == length(mat))
    mat[cbind(match(sub$row_coordinate, coords), match(sub$col_coordinate, coords))] <- as.double(sub$value)
    stopifnot(!anyNA(mat))
    mat
  }
  pred <- read_fixture("predictions.csv")
  pred$cell <- paste0("k", pred$actual_cohort, ".", pred$pair)
  pred$j_hat <- pred$jhat
  proj <- read_fixture("projection.csv")
  proj_if <- read_fixture("projection_if.csv")
  point_method <- function(method, field = "estimate") {
    if (field == "estimate") {
      sub <- pts[pts$method == method, , drop = FALSE]
      idx <- match(link_ids, sub$link)
      stopifnot(!anyNA(idx))
      return(stats::setNames(as.double(sub$value[idx]), link_ids))
    }
    prefix <- switch(field, d_plus = "denplus", d_minus = "denminus", field)
    ans <- nuis_vec[paste0(prefix, ".", link_ids)]
    stopifnot(!anyNA(ans))
    stats::setNames(as.double(ans), link_ids)
  }
  list(points = pts, point_method = point_method, nuisance_df = nuis,
    nuisance = nuis_vec, unit_if_df = uif, uif_vec = uif_vec,
    uif_method_mat = uif_method_mat, components_df = comp,
    comp_method_link = comp_method_link, covariance_df = cov_sub,
    cov_matrix = cov_matrix, joint_coords = joint_coords, predictions_df = pred,
    projection_df = proj, projection_if_df = proj_if)

}

# Build real locked installed cdid owners (panel, frame, reader, design_owner, link_views)
# and aligned reference expectations for any of the 9 fixtures F1..F9.
load_imputation_finite_reference <- function(directory,
                                             fixture_id = "F1_sep_int_x1x2",
                                             order = NULL,
                                             max_chunk_rows = 13L) {
  imp_fixture_verify_hashes(directory)
  atoms_ref <- imp_fixture_csv(directory, "atoms.csv")
  atoms_panel <- imp_fixture_csv(directory, "source_observed_panel_atoms.csv")
  atoms_ref$cohort <- as.double(atoms_ref$cohort)
  atoms_panel$cohort <- as.double(atoms_panel$cohort)
  N <- nrow(atoms_ref)
  if (is.null(order)) order <- seq_len(N)
  stopifnot(
    is.numeric(order),
    length(order) == N,
    identical(sort(as.integer(order)), seq_len(N))
  )

  meta <- imp_fixture_specs[imp_fixture_specs$fixture_id == fixture_id, , drop = FALSE]
  stopifnot(nrow(meta) == 1L)

  obs_prep <- imp_prepare_observed_atoms(atoms_ref, atoms_panel, imp_calendar)

  control <- imp_fixture_ns("ipw_owner_control")(
    max_chunk_rows = as.integer(max_chunk_rows),
    max_retained_entries = 1e7,
    max_matrix_entries = 2e7,
    max_work_ops = 1e10,
    max_output_entries = 1e7
  )

  # Construct long panel from Y_obs where unobserved outcomes (R == 0) are strictly NA_real_.
  # Keep one explicit all-NA zero-observer unit in the panel; the other 5 zero-observer
  # units enter through the linked frame with allow_frame_only = TRUE.
  Y_obs <- obs_prep$Y_obs
  long <- data.frame(
    id = rep(as.double(atoms_ref$unit_id), each = nrow(imp_calendar)),
    time = rep(as.double(imp_calendar$time), N),
    y = as.vector(t(Y_obs)),
    g = rep(as.double(atoms_ref$cohort), each = nrow(imp_calendar))
  )
  zero <- which(rowSums(!is.na(Y_obs)) == 0L)
  if (length(zero) > 1L) {
    long <- long[!long$id %in% atoms_ref$unit_id[zero[-1L]], , drop = FALSE]
  }
  panel <- imp_fixture_ns("construct_ipw_panel_owner")(
    long, "id", "time", "y", "g", control
  )
  stopifnot(
    panel$n_a_units == 24L,
    all(c(10, 22, 34) %in% panel$a_id)
  )

  or_slopes <- imp_split_pipe(meta$or_slopes)
  or_intercept <- as.integer(meta$or_intercept) == 1L
  trt_cols <- imp_split_pipe(meta$trt_cols)
  or_cols <- c(if (or_intercept) "(Intercept)", or_slopes)
  frame_cols <- unique(c(if (or_intercept || "(Intercept)" %in% trt_cols) "(Intercept)", or_slopes, setdiff(trt_cols, "(Intercept)")))

  all_X <- cbind(`(Intercept)` = rep(1, N), x1 = as.double(atoms_ref$x1), x2 = as.double(atoms_ref$x2))
  X_full <- all_X[, frame_cols, drop = FALSE]
  Z_frame <- X_full
  # Set Z_frame to NA off A so design_on_A_only = TRUE is strictly enforced on frame_owner
  Z_frame[obs_prep$A == 0L, ] <- NA_real_

  frame <- imp_fixture_ns("construct_ipw_frame_owner")(
    panel,
    as.double(atoms_ref$unit_id[order]),
    Z_frame[order, , drop = FALSE],
    weights = as.double(atoms_ref$w[order]),
    control = control,
    frame_g = as.double(atoms_ref$cohort[order]),
    allow_frame_only = TRUE,
    design_on_A_only = TRUE
  )

  # Reader holds an intercept-only probability recipe placeholder on Z to verify
  # that pure imputation uses the explicit full-F design_owner rather than reader$state$Z.
  Q_prob <- cbind(`(Intercept)` = rep(1, N))
  reader <- imp_fixture_ns("construct_cp04_frame_reader")(
    frame,
    as.double(atoms_ref$unit_id[order]),
    Q_prob[order, , drop = FALSE],
    as.double(atoms_ref$cohort[order]),
    as.double(atoms_ref$w[order]),
    control,
    U = rep(TRUE, N),
    rho = rep(1, N),
    universe = "F",
    justification = NULL,
    has_intercept = NULL
  )

  ids <- as.integer(reader$state$id)
  ord_idx <- match(ids, atoms_ref$unit_id)

  # Sealed full-F outcome/treatment design owner aligned to reader$state$id
  design_owner <- imp_fixture_ns("construct_imputation_design_owner")(
    frame_owner = frame,
    reader = reader,
    id = reader$state$id,
    X = X_full[ord_idx, , drop = FALSE],
    control = control
  )

  views <- lapply(seq_len(nrow(imp_catalog_links)), function(i) {
    lk <- imp_catalog_links[i, , drop = FALSE]
    K <- as.double(imp_split_pipe(lk$K))
    imp_fixture_ns("construct_ipw_link_view")(
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
  names(views) <- imp_catalog_links$link_id

  outcome_recipe <- imp_fixture_ns("or_pooled_outcome_recipe")(
    or_cols, or_intercept, or_slopes
  )

  exp <- imp_load_expectations(directory, fixture_id, ids, imp_catalog_links$link_id)

  list(
    fixture_id = fixture_id,
    spec = meta,
    imputation_pooling = meta$imputation_pooling,
    panel = panel,
    frame = frame,
    reader = reader,
    design_owner = design_owner,
    views = views,
    control = control,
    glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE),
    rank_tol = 1e-10,
    outcome_recipe = outcome_recipe,
    treatment_columns = trt_cols,
    raw = atoms_ref[ord_idx, , drop = FALSE],
    raw_panel = atoms_panel[ord_idx, , drop = FALSE],
    R = obs_prep$R[ord_idx, , drop = FALSE],
    Y_obs = obs_prep$Y_obs[ord_idx, , drop = FALSE],
    A = obs_prep$A[ord_idx],
    B = obs_prep$B[ord_idx],
    pair_outcomes_delta = lapply(obs_prep$pair_outcomes_delta, `[`, ord_idx),
    pair_outcomes_y = lapply(obs_prep$pair_outcomes_y, `[`, ord_idx),
    X_full = X_full[ord_idx, , drop = FALSE],
    X_or = all_X[ord_idx, or_cols, drop = FALSE],
    V_trt = all_X[ord_idx, trt_cols, drop = FALSE],
    catalog = list(
      calendar = imp_calendar,
      links = imp_catalog_links,
      link_list = lapply(seq_len(nrow(imp_catalog_links)), function(i) {
        list(
          link_id = imp_catalog_links$link_id[i],
          g = as.double(imp_catalog_links$g[i]),
          K = as.double(imp_split_pipe(imp_catalog_links$K[i])),
          s = as.double(imp_catalog_links$s[i]),
          t = as.double(imp_catalog_links$t[i]),
          or_block = imp_catalog_links$or_block[i]
        )
      })
    ),
    expectations = exp
  )
}
