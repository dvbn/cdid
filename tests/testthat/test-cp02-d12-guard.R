# Focused packaged CP02 D12 numerical-guard tests:
# - Hard-refuses absent or foreign namespaces and globalenv/dyn.load fallbacks;
#   uses only actual private functions from `asNamespace("cdid")`
# - Reproduces original zero-dual exact-fit refusal on frozen pre-D12 reference
#   (`tests/testthat/reference/chained_gmm_core_original_pre_d12_reference.R`)
#   in an isolated comparison environment
# - Verifies conditional point acceptance and caller-explicit loading-error budgeting on unchanged candidate
# - Verifies outward G_lambda = upper(||K^-1||_2) * upper(||s||_2) including both s_1 and s_2
# - Discriminates omitted s_1, omitted s_2, omitted adjoint remainder, omitted target rounding,
#   corrupt loadings despite accurate points, too-tight point/loading budgets, and nonfinite/prerequisite refusals
# - Preserves D11 separation of unit-based criterion SE from clustered sandwich SE
# - Exercises native C and pure-R interval arithmetic on identical stored candidates and scale transformations

find_d12_ref_file <- function(filename) {
  candidates <- c(
    file.path("reference", filename),
    file.path("tests", "testthat", "reference", filename)
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Cannot locate packaged D12 reference file: ", filename)
  }
  hit[1L]
}

validate_cdid_namespace_contract <- function(
  pkg = "cdid",
  eval_env = parent.frame(),
  loaded_ns = loadedNamespaces(),
  loaded_dlls = base::getLoadedDLLs(),
  global_env = globalenv()
) {
  if (!identical(pkg, "cdid") || !(pkg %in% loaded_ns)) {
    stop("Absent `cdid` namespace: packaged tests require the loaded `cdid` namespace and refuse source/dyn.load fallback.")
  }
  ns <- asNamespace(pkg)
  if (!isNamespace(ns) || !identical(unname(getNamespaceName(ns)), "cdid")) {
    stop("Foreign or invalid `cdid` namespace object.")
  }
  ns_dlls <- getNamespaceInfo(ns, "DLLs")
  if (!("cdid" %in% names(loaded_dlls)) || is.null(ns_dlls[["cdid"]])) {
    stop("Absent or unregistered `cdid` native DLL in namespace; dyn.load fallback is forbidden.")
  }
  req_fns <- c(
    "solve_chained_gmm_core",
    "certify_gmm_numerical_guard",
    "gmm_accuracy_control",
    "validate_gmm_accuracy_control",
    "gmm_bounds_control",
    "gmm_covariance_provenance",
    "gmm_if_replay_spec",
    "gmm_target_spec",
    "construct_chained_design",
    "materialize_gmm_gauge_matrix",
    "validate_gmm_targets",
    "gmm_native_backend_info",
    "gmm_require_native_dll"
  )
  for (fn_nm in req_fns) {
    if (exists(fn_nm, envir = global_env, inherits = FALSE)) {
      stop("Forbidden globalenv shadowing of private function `", fn_nm, "`; source fallback is refused.")
    }
    if (!exists(fn_nm, envir = ns, inherits = FALSE)) {
      stop("Incomplete or foreign `cdid` namespace: missing private function `", fn_nm, "`.")
    }
    fn_obj <- get(fn_nm, envir = ns, inherits = FALSE)
    if (!is.function(fn_obj) || !identical(environment(fn_obj), ns)) {
      stop("Foreign closure environment for `", fn_nm, "`; function must belong to `asNamespace('cdid')`.")
    }
    if (exists(fn_nm, envir = eval_env, inherits = TRUE)) {
      resolved <- get(fn_nm, envir = eval_env, inherits = TRUE)
      if (!identical(resolved, fn_obj)) {
        stop("Foreign binding resolved for `", fn_nm, "` outside `asNamespace('cdid')`.")
      }
    }
  }
  ns
}

read_tsv_stored_objects <- function(tsv_path) {
  tab <- utils::read.delim(tsv_path, stringsAsFactors = FALSE)
  get_vec <- function(nm) as.numeric(tab$hex[tab$name == nm])
  get_mat <- function(nm) {
    sub <- tab[tab$name == nm, , drop = FALSE]
    dm <- strsplit(sub$dim[1L], "x", fixed = TRUE)[[1L]]
    matrix(as.numeric(sub$hex), nrow = as.integer(dm[1L]), ncol = as.integer(dm[2L]))
  }
  list(tab = tab, get_vec = get_vec, get_mat = get_mat)
}

make_d12_accuracy <- function(
  ns,
  weighting = "optimal",
  target_error_unit = "abs",
  max_target_error_abs = 1e-10,
  max_target_error_se = NA_real_,
  max_loading_error_abs = NULL,
  max_kkt_residual_rel = 1e-12,
  max_loading_residual_rel = 1e-12,
  require_matching_scale = FALSE
) {
  ns$gmm_accuracy_control(
    method = if (weighting == "optimal") "certified_eigen_chol_kkt" else "certified_qr_kkt",
    target_error_unit = target_error_unit,
    max_target_error_abs = max_target_error_abs,
    max_target_error_se = max_target_error_se,
    max_orthogonality_defect = 1e-11,
    min_certified_lambda_rel = if (weighting == "optimal") 1e-8 else NA_real_,
    max_factor_residual_rel = 1e-12,
    min_rcond_R = if (weighting == "optimal") 1e-6 else NA_real_,
    min_rcond_U = 1e-12,
    min_sigma_rel_A = 1e-6,
    max_kkt_residual_rel = max_kkt_residual_rel,
    max_loading_residual_rel = max_loading_residual_rel,
    require_matching_scale = require_matching_scale,
    max_loading_error_abs = max_loading_error_abs
  )
}

test_that("D12.1 Preserved exact-fit refusal on original source and conditional acceptance on updated core", {
  ns <- validate_cdid_namespace_contract()

  # Verify hard refusal of absent namespace, foreign namespace, missing DLL, and globalenv shadowing
  expect_error(
    validate_cdid_namespace_contract(loaded_ns = character(0)),
    "Absent `cdid` namespace"
  )
  expect_error(
    validate_cdid_namespace_contract(pkg = "stats"),
    "Absent `cdid` namespace"
  )
  expect_error(
    validate_cdid_namespace_contract(loaded_dlls = list()),
    "Absent or unregistered `cdid` native DLL"
  )
  fake_global <- new.env(parent = emptyenv())
  fake_global$solve_chained_gmm_core <- function(...) NULL
  expect_error(
    validate_cdid_namespace_contract(global_env = fake_global),
    "Forbidden globalenv shadowing"
  )
  fake_eval <- new.env(parent = environment())
  fake_eval$certify_gmm_numerical_guard <- function(...) NULL
  expect_error(
    validate_cdid_namespace_contract(eval_env = fake_eval),
    "Foreign binding resolved"
  )

  # Separate process probe: confirm absent namespace is refused without source/dyn.load fallback
  probe_ns <- tempfile(fileext = ".R")
  on.exit(unlink(probe_ns), add = TRUE)
  writeLines(c(
    paste("validate_cdid_namespace_contract <-", paste(deparse(validate_cdid_namespace_contract), collapse = "\n")),
    "e <- tryCatch(validate_cdid_namespace_contract(), error = identity)",
    "stopifnot(inherits(e, 'error'), grepl('Absent `cdid` namespace', conditionMessage(e), fixed = TRUE))",
    "stopifnot(!('cdid' %in% loadedNamespaces()), !('cdid' %in% names(getLoadedDLLs())))",
    "cat('absent_namespace_refused_cleanly\\n')"
  ), probe_ns)
  out_ns <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(probe_ns)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(out_ns, "status"))
  expect_match(paste(out_ns, collapse = " "), "absent_namespace_refused_cleanly")

  tsv_path <- find_d12_ref_file("zero_dual_stored_binary64.tsv")
  expect_identical(as.integer(file.info(tsv_path)$size), 3556L)
  expect_identical(unname(tools::md5sum(tsv_path)), "784ef5f57589c7128d568f7da2985835")
  zd <- read_tsv_stored_objects(tsv_path)

  orig_ref_path <- find_d12_ref_file("chained_gmm_core_original_pre_d12_reference.R")
  expect_identical(as.integer(file.info(orig_ref_path)$size), 118998L)
  expect_identical(unname(tools::md5sum(orig_ref_path)), "51773bf8a53099c49f3459667a7ad5f5")
  orig_env <- new.env(parent = ns)
  sys.source(orig_ref_path, envir = orig_env)
  expect_identical(environment(orig_env$certify_gmm_numerical_guard), orig_env)
  expect_identical(environment(ns$certify_gmm_numerical_guard), ns)

  acc_orig <- orig_env$gmm_accuracy_control(
    method = "certified_eigen_chol_kkt",
    target_error_unit = "abs",
    max_target_error_abs = 1e-10,
    max_target_error_se = NA_real_,
    max_orthogonality_defect = 1e-11,
    min_certified_lambda_rel = 1e-8,
    max_factor_residual_rel = 1e-12,
    min_rcond_R = 1e-6,
    min_rcond_U = 1e-12,
    min_sigma_rel_A = 1e-6,
    max_kkt_residual_rel = 1e-12,
    max_loading_residual_rel = 1e-12,
    require_matching_scale = FALSE
  )

  # 1. Original pre-D12 reference refuses the exact-fit zero-dual candidate when kkt_rel_cert > max_kkt_residual_rel
  err_orig <- expect_error(
    orig_env$certify_gmm_numerical_guard(
      V_R = zd$get_mat("V_R"), R_chol = zd$get_mat("R_chol"), W_G = zd$get_mat("W_G"),
      A_whitened = zd$get_mat("A_whitened"), U_qr = zd$get_mat("U_qr"),
      piv = as.integer(zd$get_vec("piv")), mhat_R = zd$get_vec("mhat_R"),
      theta = zd$get_vec("theta"), y = zd$get_vec("y"), C_gauge = zd$get_mat("C_gauge"),
      Lambda = zd$get_mat("Lambda"), Mc = zd$get_mat("Mc"), Z1 = zd$get_mat("Z1"),
      weighting = "optimal", accuracy = acc_orig, max_matrix_entries = 50000L
    ),
    class = "cdid_gmm_accuracy_error"
  )
  expect_identical(
    conditionMessage(err_orig),
    "Certified primal/stationarity KKT residual exceeds caller threshold."
  )
  expect_identical(sprintf("%a", err_orig$diagnostics$kkt_rel_cert), "0x1.43d136248491ep-1")

  # 2. Updated core with NULL loading budget (not supplied) accepts point effects conditionally,
  # preserves kkt_rel_cert and scale intervals bit-for-bit, records threshold exceedance,
  # and leaves loading and complete inference explicitly uncertified.
  acc_pt <- make_d12_accuracy(ns, "optimal", max_target_error_abs = 1e-12, max_loading_error_abs = NULL)
  g_pt <- ns$certify_gmm_numerical_guard(
    V_R = zd$get_mat("V_R"), R_chol = zd$get_mat("R_chol"), W_G = zd$get_mat("W_G"),
    A_whitened = zd$get_mat("A_whitened"), U_qr = zd$get_mat("U_qr"),
    piv = as.integer(zd$get_vec("piv")), mhat_R = zd$get_vec("mhat_R"),
    theta = zd$get_vec("theta"), y = zd$get_vec("y"), C_gauge = zd$get_mat("C_gauge"),
    Lambda = zd$get_mat("Lambda"), Mc = zd$get_mat("Mc"), Z1 = zd$get_mat("Z1"),
    weighting = "optimal", accuracy = acc_pt, max_matrix_entries = 50000L
  )
  expect_true(g_pt$point_certified)
  expect_false(g_pt$loading_certified)
  expect_false(g_pt$loading_budget_supplied)
  expect_identical(g_pt$loading_validation_status, "incomplete_missing_caller_budget")
  expect_false(g_pt$complete_inference_certified)
  expect_identical(
    g_pt$full_inference_status,
    "pending_unvalidated_raw_root_covariance_projection_multiplier"
  )
  expect_identical(sprintf("%a", g_pt$primal_kkt$kkt_rel_cert), "0x1.43d136248491ep-1")
  expect_identical(sprintf("%a", g_pt$primal_kkt$scale1_iv), sprintf("%a", zd$get_vec("scale1_iv")))
  expect_identical(sprintf("%a", g_pt$primal_kkt$scale2_iv), sprintf("%a", zd$get_vec("scale2_iv")))
  expect_identical(sprintf("%a", g_pt$whitening_and_qr$kkt_inv_2norm_cert), sprintf("%a", zd$get_vec("kkt_inv_2norm_cert")))
  expect_true(g_pt$primal_kkt$kkt_rel_exceeds_threshold)
  expect_identical(g_pt$computed_target_values, c(-2 * (1 - 2^-52), 2, 5))
  expect_true(all(g_pt$targets$target_error_cert <= 1e-13))
  expect_true(all(g_pt$targets$G_lambda <= 1e-13))
  expect_identical(g_pt$targets$G_lambda, g_pt$targets$loading_error_cert)

  # 3. Updated core with caller-explicit max_loading_error_abs = 1e-12 certifies both point and loading
  # bounds on the stored candidate while keeping complete_inference_certified == FALSE.
  acc_both <- make_d12_accuracy(ns, "optimal", max_target_error_abs = 1e-12, max_loading_error_abs = 1e-12)
  g_both <- ns$certify_gmm_numerical_guard(
    V_R = zd$get_mat("V_R"), R_chol = zd$get_mat("R_chol"), W_G = zd$get_mat("W_G"),
    A_whitened = zd$get_mat("A_whitened"), U_qr = zd$get_mat("U_qr"),
    piv = as.integer(zd$get_vec("piv")), mhat_R = zd$get_vec("mhat_R"),
    theta = zd$get_vec("theta"), y = zd$get_vec("y"), C_gauge = zd$get_mat("C_gauge"),
    Lambda = zd$get_mat("Lambda"), Mc = zd$get_mat("Mc"), Z1 = zd$get_mat("Z1"),
    weighting = "optimal", accuracy = acc_both, max_matrix_entries = 50000L
  )
  expect_true(g_both$point_certified)
  expect_true(g_both$loading_certified)
  expect_true(g_both$loading_budget_supplied)
  expect_identical(g_both$loading_validation_status, "certified_conditional_stored_array")
  expect_false(g_both$complete_inference_certified)
})

test_that("D12.2 Outward G_lambda = upper(||K^-1||_2) * upper(||s||_2) discriminates omitted s1, omitted s2, and corrupt loadings", {
  ns <- validate_cdid_namespace_contract()
  zd <- read_tsv_stored_objects(find_d12_ref_file("zero_dual_stored_binary64.tsv"))

  W_G <- zd$get_mat("W_G")
  V_R <- zd$get_mat("V_R")
  R_chol <- zd$get_mat("R_chol")
  A_whitened <- zd$get_mat("A_whitened")
  U_qr <- zd$get_mat("U_qr")
  piv <- as.integer(zd$get_vec("piv"))
  C_gauge <- zd$get_mat("C_gauge")
  Z1 <- zd$get_mat("Z1")

  # Exact primal fit (theta* = (-2, 2, 5), mhat_R = W_G %*% theta* = (4, 2, 3), y* = (0, 0, 0)):
  # Then r_1 = 0 and r_2 = 0, so point error is ~0 for any (Lambda, Mc)!
  theta_exact <- c(-2.0, 2.0, 5.0)
  mhat_exact <- c(4.0, 2.0, 3.0)
  y_exact <- c(0.0, 0.0, 0.0)

  # True exact loadings for C_gauge = I_3 when W_G is invertible 3x3:
  Lambda_star <- solve(t(W_G))
  Mc_star <- solve(W_G, V_R %*% Lambda_star)

  # Adversary 1 (Omitted s_2 / calibration block):
  # Set Lambda_zero = 0, Mc_zero = 0. Then s_1 = W_G Mc_zero - V_R Lambda_zero = 0,
  # while s_2 = c_j - W_G' Lambda_zero = c_j != 0, and actual loading error is ||Lambda_*[, j]||_2 >= 1!
  Lambda_zero <- matrix(0.0, 3L, 3L)
  Mc_zero <- matrix(0.0, 3L, 3L)
  acc_pt_only <- make_d12_accuracy(ns, "optimal", max_target_error_abs = 1e-12, max_loading_error_abs = NULL)
  g_omitted_s2 <- ns$certify_gmm_numerical_guard(
    V_R = V_R, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_qr,
    piv = piv, mhat_R = mhat_exact, theta = theta_exact, y = y_exact,
    C_gauge = C_gauge, Lambda = Lambda_zero, Mc = Mc_zero, Z1 = Z1,
    weighting = "optimal", accuracy = acc_pt_only
  )
  # s_1 is identically 0, s_2 has norm >= 1, and G_lambda >= ||Lambda_zero - Lambda_*||_2
  expect_identical(g_omitted_s2$targets$s1_norm2_cert, c(0.0, 0.0, 0.0))
  expect_true(all(g_omitted_s2$targets$s2_norm2_cert >= 1.0))
  true_err_s2 <- sqrt(colSums((Lambda_zero - Lambda_star)^2))
  expect_true(all(g_omitted_s2$targets$G_lambda >= true_err_s2))
  # Point-only is accepted when max_loading_error_abs is NULL, but loading is NOT certified
  expect_true(g_omitted_s2$point_certified)
  expect_false(g_omitted_s2$loading_certified)
  expect_true(all(g_omitted_s2$targets$loading_cal_rel_cert >= 1.0))
  expect_true(all(g_omitted_s2$targets$loading_rel_exceeds_threshold))

  # When caller supplies a loading budget (e.g. 1e-3), the s_2-only corruption is REFUSED:
  acc_load_strict <- make_d12_accuracy(ns, "optimal", max_target_error_abs = 1e-12, max_loading_error_abs = 1e-3)
  err_s2 <- expect_error(
    ns$certify_gmm_numerical_guard(
      V_R = V_R, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_qr,
      piv = piv, mhat_R = mhat_exact, theta = theta_exact, y = y_exact,
      C_gauge = C_gauge, Lambda = Lambda_zero, Mc = Mc_zero, Z1 = Z1,
      weighting = "optimal", accuracy = acc_load_strict
    ),
    class = "cdid_gmm_accuracy_error"
  )
  expect_true("certified_loading_error_abs_exceeded" %in% err_s2$diagnostics$reasons)

  # Adversary 2 (Omitted s_1 / adjoint KKT block):
  # Keep Lambda = Lambda_* (so exact s_2 = C_gauge - W_G' Lambda_* = 0),
  # but perturb Mc_bad = Mc_* + 0.05 * I_3 so s_1 = W_G Mc_bad - V_R Lambda_* != 0!
  Mc_bad <- Mc_star + diag(0.05, 3L)
  g_omitted_s1 <- ns$certify_gmm_numerical_guard(
    V_R = V_R, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_qr,
    piv = piv, mhat_R = mhat_exact, theta = theta_exact, y = y_exact,
    C_gauge = C_gauge, Lambda = Lambda_star, Mc = Mc_bad, Z1 = Z1,
    weighting = "optimal", accuracy = acc_pt_only
  )
  expect_identical(g_omitted_s1$targets$loading_cal_abs_est, c(0.0, 0.0, 0.0))
  expect_true(all(g_omitted_s1$targets$s2_norm2_cert < 1e-15))
  expect_true(all(g_omitted_s1$targets$s1_norm2_cert >= 0.05))
  true_mc_err <- sqrt(colSums((Mc_bad - Mc_star)^2))
  expect_true(all(g_omitted_s1$targets$G_lambda >= true_mc_err))
  err_s1 <- expect_error(
    ns$certify_gmm_numerical_guard(
      V_R = V_R, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_qr,
      piv = piv, mhat_R = mhat_exact, theta = theta_exact, y = y_exact,
      C_gauge = C_gauge, Lambda = Lambda_star, Mc = Mc_bad, Z1 = Z1,
      weighting = "optimal", accuracy = acc_load_strict
    ),
    class = "cdid_gmm_accuracy_error"
  )
  expect_true("certified_loading_error_abs_exceeded" %in% err_s1$diagnostics$reasons)

  # Adversary 3 (Exact 2x1 counterexample from reference/d12/exact_reference.py):
  # W = [1; 0], V = I_2, m = (2, 0), theta = 2, y = (0, 0), c = 1, lam = (3/2, 0), mc = 3/2:
  # Here r_1 = 0, r_2 = 0, s_1 = (0, 0), s_2 = -1/2, ||lam - lam_*||_2 = 1/2 <= G_lambda!
  W2 <- matrix(c(1.0, 0.0), nrow = 2L, ncol = 1L)
  U1 <- matrix(1.0, 1L, 1L)
  C1 <- matrix(1.0, 1L, 1L)
  Lam_bad2 <- matrix(c(1.5, 0.0), nrow = 2L, ncol = 1L)
  Mc_bad2 <- matrix(1.5, 1L, 1L)
  g_ex2 <- ns$certify_gmm_numerical_guard(
    V_R = diag(2L), R_chol = diag(2L), W_G = W2, A_whitened = W2, U_qr = U1,
    piv = 1L, mhat_R = c(2.0, 0.0), theta = 2.0, y = c(0.0, 0.0),
    C_gauge = C1, Lambda = Lam_bad2, Mc = Mc_bad2, Z1 = C1,
    weighting = "optimal", accuracy = acc_pt_only
  )
  expect_identical(drop(W2 %*% Mc_bad2 - diag(2L) %*% Lam_bad2), c(0.0, 0.0))
  expect_lt(g_ex2$targets$s1_norm2_cert, 1e-15)
  expect_gte(g_ex2$targets$s2_norm2_cert, 0.5)
  expect_gte(g_ex2$targets$G_lambda, 0.5)
})

test_that("D12.3 Overidentified Cartesian adversary verifies G_c and G_lambda enclosures and budget refusals", {
  ns <- validate_cdid_namespace_contract()

  # Overidentified 4-link, 3-parameter design where W_G' lambda = c has a 1D nullspace:
  links <- data.frame(
    g = c(2, 2, 2, 2),
    s = c(0, 0, 1, 0),
    t = c(1, 2, 2, 3),
    ref = c(1, 1, 1, 1),
    functional_key = rep("k1", 4L),
    stringsAsFactors = FALSE
  )
  des <- ns$construct_chained_design(links)
  V <- diag(c(1.0, 2.0, 4.0, 1.0))
  mhat <- c(1.0, 2.5, 1.0, 3.0)
  targets <- list(
    ns$gmm_target_spec("t1", 1L, 1.0),
    ns$gmm_target_spec("t2", 2L, 1.0),
    ns$gmm_target_spec("t3", 3L, 1.0),
    ns$gmm_target_spec("t_sum", c(1L, 2L, 3L), c(1.0, 1.0, 1.0))
  )
  Phi <- matrix(c(
     1,  0,  0,  1,
    -1,  0,  0, -1,
     0,  2,  0,  0,
     0, -2,  0,  0,
     0,  0,  2,  1,
     0,  0, -2, -1,
     1, -1,  2,  0,
    -1,  1, -2,  0
  ), nrow = 8L, byrow = TRUE)
  uid <- sprintf("u%d", 1:8)
  rep_spec <- ns$gmm_if_replay_spec(
    callback = function(consume) consume(list(unit_id = uid, phi = Phi, weight = rep(1, 8L))),
    input_row_level = "unit", aggregate_to_clusters = FALSE,
    n_input_rows = 8L, n_output_rows = 8L, n_positive_contributing_rows = 8L,
    centered = TRUE, weighted = TRUE, weight_type = "unit_uniform",
    normalization_factor = 1 / 7, scale_id = "s_unit", nuisance_fit_id = "nf1",
    expected_unit_ids = uid
  )
  prov <- ns$gmm_covariance_provenance(
    row_level = "unit", n_units = 8L, n_clusters = NA_real_, n_positive_contributing_rows = 8L,
    centered = TRUE, weighted = TRUE, weight_type = "unit_uniform",
    normalization_factor = 1 / 7, scale_id = "s_unit", nuisance_fit_id = "nf1"
  )
  bnds <- ns$gmm_bounds_control(50000L, 8L, 16L, 1L, FALSE, 0L)

  for (wt in c("optimal", "identity")) {
    acc <- make_d12_accuracy(ns, wt, max_target_error_abs = 1e-10, max_loading_error_abs = 1e-10)
    res <- ns$solve_chained_gmm_core(des, mhat, V, prov, targets, wt, rep_spec, acc, bnds, "lowest")
    expect_true(res$point_certified)
    expect_true(res$loading_certified)
    expect_false(res$complete_inference_certified)

    ig <- res$internal_gauge
    W_G <- ns$materialize_gmm_gauge_matrix(ig$gauge, 50000L)
    A_ok <- if (wt == "optimal") backsolve(ig$R_chol, W_G, transpose = TRUE) else W_G
    t_info <- ns$validate_gmm_targets(des, ig$gauge, targets, bnds)

    # Exact KKT matrix K = [W_G, -V_used; 0, W_G'] where V_used = V (optimal) or I (identity):
    # Both primal (theta, y) and adjoint (Mc, Lambda) solve K * [u; v] = rhs!
    V_used <- if (wt == "optimal") V else diag(4L)
    K_mat <- rbind(cbind(W_G, -V_used), cbind(matrix(0, 3L, 3L), t(W_G)))
    sol_primal <- solve(K_mat, c(mhat, rep(0, 3L)))
    theta_star <- sol_primal[1:3]
    y_star <- sol_primal[4:7]
    sol_adj <- solve(K_mat, rbind(matrix(0, 4L, 4L), t_info$C_gauge))
    Mc_star <- sol_adj[1:3, , drop = FALSE]
    Lambda_star <- sol_adj[4:7, , drop = FALSE]
    tau_star <- drop(crossprod(t_info$C_gauge, theta_star))

    # Cartesian perturbation grid over primal and adjoint directions
    for (eps_th in c(0, 1e-6, -1e-6)) {
      for (eps_lam in c(0, 1e-5, -1e-5)) {
        th_p <- ig$theta_gauge_coords + eps_th * c(1, -0.5, 0.25)
        # Nullspace vector of W_G' on triangle {1,2,3} is (1, -1, 1, 0):
        # Shifting Lambda along (1, -1, 1, 0) keeps s_2 = c - W_G' Lambda == 0 while changing s_1!
        Lam_p <- ig$Lambda + eps_lam * matrix(c(1, -1, 1, 0), 4L, 4L)
        acc_open <- make_d12_accuracy(ns, wt, max_target_error_abs = 1.0, max_loading_error_abs = 1.0)
        g_p <- ns$certify_gmm_numerical_guard(
          V_R = V, R_chol = ig$R_chol, W_G = W_G, A_whitened = A_ok,
          U_qr = ig$U_qr, piv = ig$pivot, mhat_R = mhat,
          theta = th_p, y = ig$y_multiplier, C_gauge = t_info$C_gauge,
          Lambda = Lam_p, Mc = ig$Mc, Z1 = ig$Z1,
          weighting = wt, accuracy = acc_open
        )
        actual_pt_err <- abs(g_p$computed_target_values - tau_star)
        actual_lam_err <- sqrt(colSums((Lam_p - Lambda_star)^2) + colSums((ig$Mc - Mc_star)^2))
        expect_true(all(actual_pt_err <= g_p$targets$target_error_cert + 1e-15))
        expect_true(all(actual_lam_err <= g_p$targets$G_lambda + 1e-15))
      }
    }

    # Too-tight point budget (abs and se) and too-tight loading budget refuse cleanly
    acc_tight_pt <- make_d12_accuracy(ns, wt, max_target_error_abs = 1e-30, max_loading_error_abs = 1e-8)
    err_pt <- expect_error(
      ns$solve_chained_gmm_core(des, mhat, V, prov, targets, wt, rep_spec, acc_tight_pt, bnds, "lowest"),
      class = "cdid_gmm_accuracy_error"
    )
    expect_true("certified_target_error_abs_exceeded" %in% err_pt$diagnostics$reasons)

    if (wt == "optimal") {
      acc_tight_se <- make_d12_accuracy(
        ns, "optimal", target_error_unit = "both",
        max_target_error_abs = 1e-8, max_target_error_se = 1e-30,
        max_loading_error_abs = 1e-8, require_matching_scale = TRUE
      )
      err_se <- expect_error(
        ns$solve_chained_gmm_core(des, mhat, V, prov, targets, "optimal", rep_spec, acc_tight_se, bnds, "lowest"),
        class = "cdid_gmm_accuracy_error"
      )
      expect_true("certified_target_error_se_exceeded" %in% err_se$diagnostics$reasons)
    }

    acc_tight_ld <- make_d12_accuracy(ns, wt, max_target_error_abs = 1e-8, max_loading_error_abs = 1e-30)
    err_ld <- expect_error(
      ns$solve_chained_gmm_core(des, mhat, V, prov, targets, wt, rep_spec, acc_tight_ld, bnds, "lowest"),
      class = "cdid_gmm_accuracy_error"
    )
    expect_true("certified_loading_error_abs_exceeded" %in% err_ld$diagnostics$reasons)
  }
})

test_that("D12.4 Target rounding term in G_c, adjoint remainder witness, and D11 separation of unit criterion SE from clustered sandwich SE", {
  ns <- validate_cdid_namespace_contract()

  # 1. Verify target rounding term |tau_round - c' theta| is included in G_c even when r == 0 and s == 0:
  # Construct 3x3 exact system where theta = (1, 2^53, -2^53), c = (1, 1, 1):
  # Exact c' theta = 1 + 2^53 - 2^53 = 1, while left-to-right fl((1 + 2^53) - 2^53) = 0 (error = 1)!
  W_G <- diag(3L)
  V_R <- diag(3L)
  R_chol <- diag(3L)
  A_whitened <- diag(3L)
  U_qr <- diag(3L)
  piv <- 1:3
  theta <- c(1.0, 2^53, -2^53)
  mhat_R <- theta
  y <- c(0.0, 0.0, 0.0)
  C_gauge <- matrix(c(1.0, 1.0, 1.0), nrow = 3L, ncol = 1L)
  Lambda <- C_gauge
  Mc <- C_gauge
  Z1 <- C_gauge

  acc_open <- make_d12_accuracy(ns, "optimal", max_target_error_abs = 20.0, max_loading_error_abs = 1e-10)
  g_rnd <- ns$certify_gmm_numerical_guard(
    V_R = V_R, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_qr,
    piv = piv, mhat_R = mhat_R, theta = theta, y = y,
    C_gauge = C_gauge, Lambda = Lambda, Mc = Mc, Z1 = Z1,
    weighting = "optimal", accuracy = acc_open
  )
  expect_identical(g_rnd$computed_target_values, 0.0)
  expect_gte(g_rnd$targets$target_round_cert, 1.0)
  expect_gte(g_rnd$targets$target_error_cert, 1.0)

  # Tight point budget (1e-6) refuses because of the target rounding error!
  acc_tight <- make_d12_accuracy(ns, "optimal", max_target_error_abs = 1e-6, max_loading_error_abs = 1e-10)
  expect_error(
    ns$certify_gmm_numerical_guard(
      V_R = V_R, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_qr,
      piv = piv, mhat_R = mhat_R, theta = theta, y = y,
      C_gauge = C_gauge, Lambda = Lambda, Mc = Mc, Z1 = Z1,
      weighting = "optimal", accuracy = acc_tight
    ),
    class = "cdid_gmm_accuracy_error"
  )

  # 2. Verify approximate-adjoint remainder term when direct products |Lambda'r1| + |Mc'r2| == 0
  # (matching the 2x1 witness in reference/d12/exact_reference.py: W = [1; 1], V = I_2, m = (1, 1),
  # theta = 0, y = (0, 0), lam = (0, 0), mc = 0, c = 1 -> direct indicator = 0, true error = 1):
  W_rem <- matrix(c(1.0, 1.0), nrow = 2L, ncol = 1L)
  U_rem <- matrix(sqrt(2.0), nrow = 1L, ncol = 1L)
  C_rem <- matrix(1.0, nrow = 1L, ncol = 1L)
  g_rem <- ns$certify_gmm_numerical_guard(
    V_R = diag(2L), R_chol = diag(2L), W_G = W_rem, A_whitened = W_rem, U_qr = U_rem,
    piv = 1L, mhat_R = c(1.0, 1.0), theta = 0.0, y = c(0.0, 0.0),
    C_gauge = C_rem, Lambda = matrix(0.0, 2L, 1L), Mc = matrix(0.0, 1L, 1L), Z1 = C_rem,
    weighting = "optimal", accuracy = make_d12_accuracy(ns, "optimal", max_target_error_abs = 10.0, max_loading_error_abs = NULL)
  )
  expect_identical(g_rem$targets$target_error_est, 0.0)
  expect_identical(g_rem$targets$direct_primal_term_cert, 0.0)
  expect_identical(g_rem$targets$direct_stationarity_term_cert, 0.0)
  expect_gte(g_rem$targets$adjoint_missing_term_cert, 1.0)
  expect_gte(g_rem$targets$target_error_cert, 1.0)

  # 3. D11 separation (matching reference/covariance/reference.py):
  # unit-based GMM point weighting V_unit vs cluster sandwich V_cluster
  Psi_u <- matrix(c(
     2,  1,  1,
     2,  1, -1,
    -2, -1,  1,
    -2, -1, -1,
     1,  2,  1,
     1,  2, -1,
    -1, -2,  1,
    -1, -2, -1
  ), nrow = 8L, byrow = TRUE)
  V_unit <- crossprod(Psi_u)
  links_cov <- data.frame(
    g = c(2, 2, 2), s = c(0, 0, 1), t = c(1, 2, 2),
    ref = c(0, 0, 0), functional_key = rep("k1", 3L), stringsAsFactors = FALSE
  )
  des_cov <- ns$construct_chained_design(links_cov)
  mhat_cov <- c(1.0, 2.0, 1.5)
  targets_cov <- list(ns$gmm_target_spec("tau1", 1L, 1.0), ns$gmm_target_spec("tau2", 2L, 1.0))

  uid_cov <- sprintf("u%02d", 1:8)
  cid_cov <- rep(sprintf("c%d", 1:4), each = 2L)
  rep_cluster <- ns$gmm_if_replay_spec(
    callback = function(consume) {
      consume(list(unit_id = uid_cov, cluster_id = cid_cov, phi = Psi_u, weight = rep(1.0, 8L), nuisance_fit_id = "nf_cov"))
    },
    input_row_level = "unit", aggregate_to_clusters = TRUE,
    n_input_rows = 8L, n_output_rows = 4L, n_positive_contributing_rows = 4L,
    centered = TRUE, weighted = TRUE, weight_type = "cluster_multiplier",
    normalization_factor = 1.0, scale_id = "cluster_sum_scale", nuisance_fit_id = "nf_cov",
    expected_unit_ids = uid_cov, expected_cluster_ids = sprintf("c%d", 1:4)
  )
  prov_unit <- ns$gmm_covariance_provenance(
    row_level = "unit", n_units = 8L, n_clusters = NA_real_, n_positive_contributing_rows = 8L,
    centered = TRUE, weighted = TRUE, weight_type = "unit_uniform",
    normalization_factor = 1.0, scale_id = "unit_sum_scale", nuisance_fit_id = "nf_cov"
  )
  bnds_cov <- ns$gmm_bounds_control(50000L, 8L, 16L, 1L, FALSE, 0L)
  acc_cov <- make_d12_accuracy(ns, "optimal", max_target_error_abs = 1e-10, max_loading_error_abs = 1e-10)

  res_d11 <- ns$solve_chained_gmm_core(
    design = des_cov, mhat = mhat_cov, V = V_unit, provenance = prov_unit,
    targets = targets_cov, weighting = "optimal", if_replay = rep_cluster,
    accuracy = acc_cov, bounds = bnds_cov, gauge_rule = "lowest"
  )
  expect_equal(res_d11$targets$estimate, c(7 / 8, 17 / 8), tolerance = 1e-14)
  expect_equal(res_d11$targets$se_criterion_conditional, rep(sqrt(19), 2L), tolerance = 1e-14)
  expect_equal(res_d11$targets$se_sandwich, rep(sqrt(37), 2L), tolerance = 1e-14)
  expect_gt(min(res_d11$targets$se_sandwich / res_d11$targets$se_criterion_conditional), 1.39)
  expect_false(res_d11$scale_matches_inference)
  expect_false(res_d11$declared_scale_matches_inference)
  expect_false(res_d11$complete_inference_certified)

  # 4. Nonfinite candidates, overflowed interval bounds, and violated PD/QR prerequisites still refuse under D12:
  theta_inf <- c(Inf, 0.0, 0.0)
  expect_error(
    ns$certify_gmm_numerical_guard(
      V_R = V_R, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_qr,
      piv = piv, mhat_R = mhat_R, theta = theta_inf, y = y,
      C_gauge = C_gauge, Lambda = Lambda, Mc = Mc, Z1 = Z1,
      weighting = "optimal", accuracy = acc_open
    ),
    class = "cdid_gmm_accuracy_error"
  )
  theta_huge <- rep(.Machine$double.xmax, 3L)
  expect_error(
    ns$certify_gmm_numerical_guard(
      V_R = V_R, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_qr,
      piv = piv, mhat_R = mhat_R, theta = theta_huge, y = y,
      C_gauge = C_gauge, Lambda = Lambda, Mc = Mc, Z1 = Z1,
      weighting = "optimal", accuracy = acc_open
    ),
    class = "cdid_gmm_accuracy_error"
  )
  V_nonpd <- diag(c(-1.0, 1.0, 1.0))
  expect_error(
    ns$certify_gmm_numerical_guard(
      V_R = V_nonpd, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_qr,
      piv = piv, mhat_R = mhat_R, theta = theta, y = y,
      C_gauge = C_gauge, Lambda = Lambda, Mc = Mc, Z1 = Z1,
      weighting = "optimal", accuracy = acc_open
    ),
    class = "cdid_gmm_pd_error"
  )
  U_sing <- diag(c(1.0, 0.0, 1.0))
  expect_error(
    ns$certify_gmm_numerical_guard(
      V_R = V_R, R_chol = R_chol, W_G = W_G, A_whitened = A_whitened, U_qr = U_sing,
      piv = piv, mhat_R = mhat_R, theta = theta, y = y,
      C_gauge = C_gauge, Lambda = Lambda, Mc = Mc, Z1 = Z1,
      weighting = "optimal", accuracy = acc_open
    ),
    class = "cdid_gmm_qr_error"
  )
})

test_that("D12.5 Native C and pure-R interval primitives agree bit-for-bit on G_lambda and scale transformations", {
  ns <- validate_cdid_namespace_contract()
  pure_r_path <- find_d12_ref_file("chained_gmm_core_pure_r_reference.R")
  expect_identical(as.integer(file.info(pure_r_path)$size), 106617L)
  expect_identical(unname(tools::md5sum(pure_r_path)), "21307628a8060a6c223dcabfc0c74600")
  pure_r <- new.env(parent = ns)
  sys.source(pure_r_path, envir = pure_r)
  raw_bits <- function(x) writeBin(as.double(x), raw(), size = 8L, endian = "little")

  zd <- read_tsv_stored_objects(find_d12_ref_file("zero_dual_stored_binary64.tsv"))
  W_G <- zd$get_mat("W_G")
  V_R <- zd$get_mat("V_R")
  R_chol <- zd$get_mat("R_chol")
  A_whitened <- zd$get_mat("A_whitened")
  U_qr <- zd$get_mat("U_qr")
  piv <- as.integer(zd$get_vec("piv"))
  mhat_R <- zd$get_vec("mhat_R")
  theta <- zd$get_vec("theta")
  y <- zd$get_vec("y")
  C_gauge <- zd$get_mat("C_gauge")
  Lambda <- zd$get_mat("Lambda")
  Mc <- zd$get_mat("Mc")
  Z1 <- zd$get_mat("Z1")

  for (sc in c(2^-6, 1.0, 2^6)) {
    V_sc <- V_R * (sc^2)
    R_sc <- R_chol * sc
    A_sc <- A_whitened / sc
    U_sc <- U_qr / sc
    m_sc <- mhat_R * sc
    th_sc <- theta * sc
    y_sc <- y / sc
    Mc_sc <- Mc * (sc^2)
    Z1_sc <- Z1 * sc

    acc <- make_d12_accuracy(ns, "optimal", max_target_error_abs = 1e-10 * sc, max_loading_error_abs = 1e-10 * (1 + sc^2))
    g_c <- ns$certify_gmm_numerical_guard(
      V_R = V_sc, R_chol = R_sc, W_G = W_G, A_whitened = A_sc, U_qr = U_sc,
      piv = piv, mhat_R = m_sc, theta = th_sc, y = y_sc,
      C_gauge = C_gauge, Lambda = Lambda, Mc = Mc_sc, Z1 = Z1_sc,
      weighting = "optimal", accuracy = acc
    )
    expect_true(g_c$point_certified)
    expect_true(g_c$loading_certified)

    # Pure-R computation of s_norm2_cert and G_lambda on the identical candidate matches native C bit-for-bit
    kkt_inv_c <- g_c$whitening_and_qr$kkt_inv_2norm_cert
    for (j in 1:3) {
      c_j <- C_gauge[, j]; lam_j <- Lambda[, j]; mc_j <- Mc_sc[, j]
      s1_lo_r <- s1_hi_r <- numeric(3L)
      for (i in 1:3) {
        w_mc_r <- pure_r$gmm_iv_dot(W_G[i, ], mc_j)
        v_lam_r <- pure_r$gmm_iv_dot(V_sc[i, ], lam_j)
        d_r <- pure_r$gmm_iv_sub(w_mc_r, v_lam_r)
        s1_lo_r[i] <- d_r[1L]; s1_hi_r[i] <- d_r[2L]
      }
      s2_r <- pure_r$gmm_iv_residual(t(W_G), lam_j, c_j)
      s2_lo_r <- -s2_r$hi; s2_hi_r <- -s2_r$lo
      s_norm_r <- pure_r$gmm_iv_norm(c(s1_lo_r, s2_lo_r), c(s1_hi_r, s2_hi_r))[2L]
      G_lam_r <- pure_r$gmm_iv_mul(pure_r$gmm_iv(kkt_inv_c), pure_r$gmm_iv(s_norm_r))[2L]
      expect_identical(raw_bits(g_c$targets$s_norm2_cert[j]), raw_bits(s_norm_r))
      expect_identical(raw_bits(g_c$targets$G_lambda[j]), raw_bits(G_lam_r))
    }
  }
})
