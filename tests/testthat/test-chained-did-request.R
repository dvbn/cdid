# Focused integration tests for the CP02 development request path
# (R/chained_did_request.R) and its narrow pipeline/adapter composition
# (bstrap/cband separation including the released bstrap = FALSE, cband = TRUE
# endpoint laws, caller-RNG draws, released selectors, no-refit aggregation).
# Toy RNG use is explicit: seeds 3/5/7/11/31 draw Rademacher signs with
# biters 24/32 (software verification only, no coverage claim). Sections:
#   A. pure-R composition checks (no did, no native DLL): selector oracle,
#      flag/RNG laws, request preparation, option resolution, frozen N200
#      no-refit aggregation parity;
#   B. installed-host checks (did 2.5.1 + cdid DLL): bounded end-to-end requests
#      on the reviewed 48-unit development panel with the reviewed controls.
# Controls below are explicit development inputs copied from
# inst/scripts/cp02_pipeline_dev.R; they are not public defaults.

req_locate_fixture <- function(filename) {
  candidates <- c(
    tryCatch(testthat::test_path("fixtures", filename), error = function(e) character(0L)),
    file.path("fixtures", filename), file.path("tests", "testthat", "fixtures", filename),
    file.path("..", "..", "tests", "testthat", "fixtures", filename))
  for (p in candidates) if (nzchar(p) && file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  stop("Fixture file not found: ", filename)
}

req_budget <- function() list(max_units = 2000, max_clusters = 2000, max_targets = 100, max_draws = 2000,
  max_if_batch_rows = 64, max_if_batch_cells = 6400, max_sign_batch_rows = 32, max_sign_batch_cells = 64000,
  max_workspace_cells = 1000000, max_workspace_bytes = 8000000)
req_accuracy <- function() list(max_estimate_error_abs = 1e-3, max_se_error_abs = 1e-2, max_vcov_error_abs = 1e-2,
  max_multiplier_error_abs = 1e-1, max_band_error_abs = 1e-1)

make_request_test_control <- function(require_matching_scale = TRUE) {
  lnk_ctrl <- ipw_link_enclosure_control(max_chunk_rows = 64L, max_retained_entries = 2000000L,
    max_matrix_entries = 2000000L, max_work_ops = 500000000L, max_output_entries = 500000L)
  ipw_pipeline_control(
    owner_control = ipw_owner_control(max_chunk_rows = 64L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = 500000000L, max_output_entries = 500000L),
    catalog_bounds = list(max_candidates = 500L, max_rows = 500L, max_pool_entries = 5000L, max_scratch = 100000L, max_work = 500000000L),
    class_bounds = list(max_links = 500L, max_classes = 50L, max_class_dates = 50L, max_E_entries = 5000L, max_mask_entries = 100000L),
    fitted_control = ipw_fitted_view_control(max_chunk_rows = 64L, max_observer_units = 2000L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = 500000000L, max_output_entries = 500000L),
    point_control = ipw_point_view_control(max_chunk_rows = 64L, max_observer_units = 2000L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = 500000000L, max_output_entries = 500000L),
    glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    certificate_control = logit_certificate_control(taylor_degree = 12L, max_abs_eta = 50.0, max_squarings = 16L,
      max_chunk_rows = 64L, max_matrix_entries = 2000000L, max_work_ops = 500000000L),
    problem_control = ipw_nuisance_problem_control(max_chunk_rows = 64L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = 500000000L),
    link_control = lnk_ctrl,
    covariance_control = ipw_covariance_transfer_control(max_chunk_rows = 64L, max_operator_block_rows = 64L,
      max_cluster_rows = 2000L, max_matrix_entries = 2000000L, max_work_ops = 500000000L, max_output_entries = 1000000L,
      link_control = lnk_ctrl),
    gmm_accuracy = gmm_accuracy_control(method = "certified_eigen_chol_kkt", target_error_unit = "abs",
      max_target_error_abs = 1e-4, max_target_error_se = 1e-2, max_orthogonality_defect = 1e-8,
      min_certified_lambda_rel = 1e-10, max_factor_residual_rel = 1e-10, min_rcond_R = 1e-10, min_rcond_U = 1e-10,
      min_sigma_rel_A = 1e-10, max_kkt_residual_rel = 1e-10, max_loading_residual_rel = 1e-10,
      require_matching_scale = require_matching_scale, max_loading_error_abs = 1e-4),
    gmm_bounds = gmm_bounds_control(max_matrix_entries = 2000000L, max_batch_rows = 64L, max_targets = 100L,
      max_open_clusters = 2000L, allow_materialize_if = TRUE, max_materialized_if_entries = 1000000L),
    multiplier_budget = req_budget(), inference_accuracy = req_accuracy())
}

req_auto_witness <- function(fv) list(F = unname(fv$factors$H_chol), r = logit_cert_automatic_radius_request())

# Reviewed deterministic 48-unit development panel (inst/scripts/cp02_pipeline_dev.R),
# in PUBLIC coding: never treated = 0, plus a 52-unit frame with four frame-only
# units and a zero-eligible cluster 13.
req_dev_inputs <- function() {
  ids <- 1:48; gs <- rep(c(3, 4, 0), each = 16); xs <- rep(c(-1, 1, -.5, .5), 12); w <- rep(c(1, 1.5, .8, 1.2), 12)
  panel <- expand.grid(id = ids, t = 1:4)
  panel$g <- gs[panel$id]; panel$x <- xs[panel$id]; panel$weight <- w[panel$id]
  panel$y <- with(panel, .25 * id + .5 * t + .3 * x * t +
    ifelse(g == 3 & t >= 3, 2 + .5 * (t - 3), ifelse(g == 4 & t >= 4, 1.5, 0)) + .2 * sin(id * t + t * t))
  panel$cluster <- rep(1:12, each = 4)[panel$id]
  frame <- data.frame(id = 1:52, g = c(gs, 3, 4, 0, 0), x = c(xs, -1, 1, -.5, .5), weight = c(w, 1, 1, 1, 1),
    cluster = c(rep(1:12, each = 4), rep(13L, 4)))
  list(panel = panel, frame = frame)
}

capture_warnings <- function(expr) {
  w <- character()
  val <- withCallingHandlers(expr, warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  list(value = val, warnings = w)
}

# Frozen released-did 2.5.1 selector oracle (plan/evidence/cp02_selector_reference;
# inputs.rds sha256 dde74ccf92a70d5380a99ce5de1fd596620eb27fed61a907ce161ac6cd02ecca,
# results.rds sha256 8007cae41be145f14ef8407de6666b87167c82d50b69020a548f4b7f6331d5fa).
# Generated by verification/generate_oracle_constants.R; numbers are exact copies.
oracle_inputs <- function() {
  att <- c(0, 6.75, -1.5, 8.25, -3, 7.5, -0.25, 9, -1.75, 10.5, 1, 9.75, -0.5, 11.25, -2)
  names(att) <- c("g3_t2", "g3_t3", "g3_t4", "g3_t5", "g3_t6", "g4_t2", "g4_t3", "g4_t4", "g4_t5", "g4_t6", "g5_t2", "g5_t3", "g5_t4", "g5_t5", "g5_t6")
  inffunc <- matrix(c(-1.0625, -0.3125, -0.5625, -0.1875, -0.4375, 0.3125, -0.3125, 0.4375, 0.1875, 0.5625, 0.3125, 1.0625, -1.59375, -0.46875, -0.84375, -0.28125, -0.65625, 0.46875, -0.46875, 0.65625, 0.28125, 0.84375, 0.46875, 1.59375, -1.75, -1, -0.75, -0.75, -0.5, 0.25, -0.25, 0.5, 0.75, 0.75, 1, 1.75, -2.03125, -1.15625, -1.28125, -0.59375, -0.71875, 0.15625, -0.15625, 0.71875, 0.59375, 1.28125, 1.15625, 2.03125, -2.5625, -1.3125, -1.5625, -0.6875, -0.9375, 0.3125, -0.3125, 0.9375, 0.6875, 1.5625, 1.3125, 2.5625, -2.71875, -1.84375, -1.46875, -1.15625, -0.78125, 0.09375, -0.09375, 0.78125, 1.15625, 1.46875, 1.84375, 2.71875, -3.25, -2, -1.75, -1.25, -1, 0.25, -0.25, 1, 1.25, 1.75, 2, 3.25, -3.53125, -2.15625, -2.28125, -1.09375, -1.21875, 0.15625, -0.15625, 1.21875, 1.09375, 2.28125, 2.15625, 3.53125, -3.6875, -2.6875, -2.1875, -1.5625, -1.0625, -0.0625, 0.0625, 1.0625, 1.5625, 2.1875, 2.6875, 3.6875, -4.21875, -2.84375, -2.46875, -1.65625, -1.28125, 0.09375, -0.09375, 1.28125, 1.65625, 2.46875, 2.84375, 4.21875, -4.75, -3, -2.75, -1.75, -1.5, 0.25, -0.25, 1.5, 1.75, 2.75, 3, 4.75, -4.65625, -3.53125, -2.90625, -1.96875, -1.34375, -0.21875, 0.21875, 1.34375, 1.96875, 2.90625, 3.53125, 4.65625, -5.1875, -3.6875, -3.1875, -2.0625, -1.5625, -0.0625, 0.0625, 1.5625, 2.0625, 3.1875, 3.6875, 5.1875, -5.71875, -3.84375, -3.46875, -2.15625, -1.78125, 0.09375, -0.09375, 1.78125, 2.15625, 3.46875, 3.84375, 5.71875, -5.875, -4.375, -3.375, -2.625, -1.625, -0.125, 0.125, 1.625, 2.625, 3.375, 4.375, 5.875), 12, 15)
  dimnames(inffunc) <- list(c("101", "104", "107", "110", "113", "116", "119", "122", "125", "128", "131", "134"), c("g3_t2", "g3_t3", "g3_t4", "g3_t5", "g3_t6", "g4_t2", "g4_t3", "g4_t4", "g4_t5", "g4_t6", "g5_t2", "g5_t3", "g5_t4", "g5_t5", "g5_t6"))
  list(att = att, inffunc = inffunc, frame_id = c(101, 104, 107, 110, 113, 116, 119, 122, 125, 128, 131, 134), frame_cohort = c(0, 0, 0, 3, 3, 3, 4, 4, 5, 5, 5, 5), frame_weight = c(1, 2, 3, 1, 4, 3, 2, 4, 1, 2, 4, 5),
       cell_group = c(3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5), cell_t = c(2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6), cell_id = c("g3_t2", "g3_t3", "g3_t4", "g3_t5", "g3_t6", "g4_t2", "g4_t3", "g4_t4", "g4_t5", "g4_t6", "g5_t2", "g5_t3", "g5_t4", "g5_t5", "g5_t6"),
       times_regular = c(1, 2, 3, 4, 5, 6), times_irregular = c(2000, 2002, 2005, 2009, 2010, 2014),
       irregular_cohort = c(0, 0, 0, 2005, 2005, 2005, 2009, 2009, 2010, 2010, 2010, 2010), irregular_cell_group = c(2005, 2005, 2005, 2005, 2005, 2009, 2009, 2009, 2009, 2009, 2010, 2010, 2010, 2010, 2010), irregular_cell_t = c(2002, 2005, 2009, 2010, 2014, 2002, 2005, 2009, 2010, 2014, 2002, 2005, 2009, 2010, 2014), irregular_cell_id = c("g2005_t2002", "g2005_t2005", "g2005_t2009", "g2005_t2010", "g2005_t2014", "g2009_t2002", "g2009_t2005", "g2009_t2009", "g2009_t2010", "g2009_t2014", "g2010_t2002", "g2010_t2005", "g2010_t2009", "g2010_t2010", "g2010_t2014"))
}
oracle_expectations <- function() list(
  regular_simple = list(input = "regular", type = "simple", balance_e = NULL, min_e = -Inf, max_e = Inf, status = "ok", warnings = c(), bin_ids = list(c("g3_t3", "g3_t4", "g3_t5", "g3_t6", "g4_t4", "g4_t5", "g4_t6", "g5_t5", "g5_t6")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t3", "g3_t4", "g3_t5", "g3_t6", "g4_t4", "g4_t5", "g4_t6", "g5_t5", "g5_t6"), overall_att = 4.0743243243243246, overall_se = 0.93000136250812182, egt = NULL, att_egt = NULL, se_egt = NULL, crit = NULL, overall_if = c(-3.6655405405405408, -2.3817567567567566, -2.1520270270270272, -2.3151022644265886, -4.9056793279766246, -2.6817932797662523, 1.6540357925493065, 4.7303688823959087, 1.5535975164353542, 2.5092220598977359, 3.0961468224981736, 4.5585281227173127), bin_if = NULL, result_min_e = NULL, result_max_e = NULL, result_balance_e = NULL),
  regular_simple_max1_other_selectors = list(input = "regular", type = "simple", balance_e = 2, min_e = 1, max_e = 1, status = "ok", warnings = c(), bin_ids = list(c("g3_t3", "g3_t4", "g4_t4", "g4_t5", "g5_t5", "g5_t6")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t3", "g3_t4", "g4_t4", "g4_t5", "g5_t5", "g5_t6"), overall_att = 3.7788461538461537, overall_se = 0.89359806418597587, egt = NULL, att_egt = NULL, se_egt = NULL, crit = NULL, overall_if = c(-4.0228365384615392, -2.6814903846153846, -2.3401442307692308, -2.1010539940828408, -3.3573409763313613, -1.483450443786982, -0.25619452662721942, 0.9431397928994093, 1.9590421597633134, 3.1212093195266273, 4.243620562130177, 5.975499260355031), bin_if = NULL, result_min_e = NULL, result_max_e = NULL, result_balance_e = NULL),
  regular_group = list(input = "regular", type = "group", balance_e = NULL, min_e = -Inf, max_e = Inf, status = "ok", warnings = c(), bin_ids = list(c("g3_t3", "g3_t4", "g3_t5", "g3_t6"), c("g4_t4", "g4_t5", "g4_t6"), c("g5_t5", "g5_t6")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t3", "g3_t4", "g3_t5", "g3_t6", "g4_t4", "g4_t5", "g4_t6", "g5_t5", "g5_t6"), overall_att = 4.3076923076923075, overall_se = 0.93934520095057983, egt = c(3, 4, 5), att_egt = c(2.625, 5.916666666666667, 4.625), se_egt = c(0.31300956718494333, 0.64474033856015767, 0.99174312008334042), crit = 1.9599639845400534, overall_if = c(-4.165865384615385, -2.7908653846153846, -2.4543269230769234, -2.3896079881656811, -4.3829511834319526, -2.231323964497042, 1.3866494082840235, 4.2468565088757408, 1.759430473372781, 2.7472263313609471, 3.3766642011834311, 4.8981139053254443), bin_if = matrix(c(-1.984375, -0.984375, -1.109375, -0.578125, -0.703125, 0.296875, -0.296875, 0.703125, 0.578125, 1.109375, 0.984375, 1.984375, -3.8125, -2.5625, -2.3125, -1.4375, -1.1875, 0.062499999999999993, -0.062499999999999993, 1.1875, 1.4375, 2.3125, 2.5625, 3.8125, -5.796875, -4.109375, -3.421875, -2.390625, -1.703125, -0.015625, 0.015625, 1.703125, 2.390625, 3.421875, 4.109375, 5.796875), 12, 3), result_min_e = NULL, result_max_e = NULL, result_balance_e = NULL),
  regular_group_max1_other_selectors = list(input = "regular", type = "group", balance_e = 2, min_e = 1, max_e = 1, status = "ok", warnings = c(), bin_ids = list(c("g3_t3", "g3_t4"), c("g4_t4", "g4_t5"), c("g5_t5", "g5_t6")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t3", "g3_t4", "g4_t4", "g4_t5", "g5_t5", "g5_t6"), overall_att = 3.7788461538461537, overall_se = 0.89359806418597587, egt = c(3, 4, 5), att_egt = c(2.625, 3.625, 4.625), se_egt = c(0.25546065332905099, 0.6118350894230451, 0.99174312008334042), crit = 1.9599639845400534, overall_if = c(-4.0228365384615383, -2.6814903846153846, -2.3401442307692308, -2.1010539940828408, -3.3573409763313613, -1.4834504437869827, -0.2561945266272192, 0.94313979289940875, 1.9590421597633134, 3.1212093195266273, 4.243620562130177, 5.9754992603550292), bin_if = matrix(c(-1.671875, -0.734375, -0.796875, -0.515625, -0.578125, 0.359375, -0.359375, 0.578125, 0.515625, 0.796875, 0.734375, 1.671875, -3.609375, -2.421875, -2.234375, -1.328125, -1.140625, 0.046875, -0.046875, 1.140625, 1.328125, 2.234375, 2.421875, 3.609375, -5.796875, -4.109375, -3.421875, -2.390625, -1.703125, -0.015625, 0.015625, 1.703125, 2.390625, 3.421875, 4.109375, 5.796875), 12, 3), result_min_e = NULL, result_max_e = NULL, result_balance_e = NULL),
  regular_dynamic = list(input = "regular", type = "dynamic", balance_e = NULL, min_e = -Inf, max_e = Inf, status = "ok", warnings = c(), bin_ids = list(c("g5_t2"), c("g4_t2", "g5_t3"), c("g3_t2", "g4_t3", "g5_t4"), c("g3_t3", "g4_t4", "g5_t5"), c("g3_t4", "g4_t5", "g5_t6"), c("g3_t5", "g4_t6"), c("g3_t6")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t2", "g3_t3", "g3_t4", "g3_t5", "g3_t6", "g4_t2", "g4_t3", "g4_t4", "g4_t5", "g4_t6", "g5_t2", "g5_t3", "g5_t4", "g5_t5", "g5_t6"), overall_att = 3.4429945054945059, overall_se = 0.71250157280232096, egt = c(-3, -2, -1, 0, 1, 2, 3), att_egt = c(1, 9, -0.28846153846153849, 9.3461538461538467, -1.7884615384615385, 9.2142857142857153, -3), se_egt = c(0.78616509433805026, 0.95128592052278194, 0.53943109681143719, 1.2063781190787524, 0.65870945618837951, 0.80173648357935978, 0.41129347963332574), crit = 1.9599639845400534, overall_if = c(-3.3942307692307696, -2.1387362637362637, -2.0082417582417582, -1.6913114358169306, -2.9795314575534357, -1.2511321096485928, 0.31243207342108475, 2.0479410699190925, 1.4136728655959425, 2.3987743026204567, 2.9198013524936597, 4.3705621301775155), bin_if = matrix(c(-4.75, -3, -2.75, -1.75, -1.5, 0.25, -0.25, 1.5, 1.75, 2.75, 3, 4.75, -4.0104166666666661, -2.96875, -2.427083333333333, -1.6979166666666665, -1.15625, -0.11458333333333333, -1.8854166666666667, -2.8437499999999982, 2.1979166666666665, 3.427083333333333, 4.96875, 6.5104166666666679, -3.4711538461538463, -2.2596153846153846, -2.0480769230769234, -1.1649408284023668, -0.55399408284023677, 0.52440828402366868, -0.089497041420118328, 1.1575443786982249, 1.2004437869822486, 1.8528106508875741, 1.8690828402366864, 2.9829881656804735, -3.9447115384615388, -2.4158653846153846, -2.3870192307692308, -2.5323594674556222, -6.0981878698224854, -3.3711168639053253, -0.54308431952662761, 0.66623520710059292, 2.2128328402366861, 4.1444156804733732, 5.9306582840236679, 8.3382026627218959, -4.1009615384615392, -2.9471153846153846, -2.2932692307692308, -1.6697485207100591, -0.61649408284023677, 0.40421597633136108, 0.030695266272189423, 1.2200443786982247, 1.7052514792899409, 2.0980029585798814, 2.5565828402366861, 3.6127958579881665, -2.9687499999999996, -1.8794642857142856, -1.7901785714285714, -1.8756377551020407, -4.2659438775510203, -2.3501275510204067, 2.0746173469387772, 5.3679846938775526, 1.0491071428571428, 1.7901785714285714, 1.8794642857142856, 2.9687499999999996, -2.5625, -1.3125, -1.5625, -0.6875, -0.9375, 0.3125, -0.3125, 0.9375, 0.6875, 1.5625, 1.3125, 2.5625), 12, 7), result_min_e = -Inf, result_max_e = Inf, result_balance_e = NULL),
  regular_dynamic_window = list(input = "regular", type = "dynamic", balance_e = NULL, min_e = -2, max_e = 1, status = "ok", warnings = c(), bin_ids = list(c("g4_t2", "g5_t3"), c("g3_t2", "g4_t3", "g5_t4"), c("g3_t3", "g4_t4", "g5_t5"), c("g3_t4", "g4_t5", "g5_t6")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t2", "g3_t3", "g3_t4", "g4_t2", "g4_t3", "g4_t4", "g4_t5", "g5_t3", "g5_t4", "g5_t5", "g5_t6"), overall_att = 3.7788461538461542, overall_se = 0.89359806418597587, egt = c(-2, -1, 0, 1), att_egt = c(9, -0.28846153846153849, 9.3461538461538467, -1.7884615384615385), se_egt = c(0.95128592052278194, 0.53943109681143719, 1.2063781190787524, 0.65870945618837951), crit = 1.9599639845400534, overall_if = c(-4.0228365384615392, -2.6814903846153846, -2.3401442307692308, -2.1010539940828408, -3.3573409763313613, -1.4834504437869822, -0.25619452662721909, 0.94313979289940875, 1.9590421597633134, 3.1212093195266273, 4.243620562130177, 5.975499260355031), bin_if = matrix(c(-4.0104166666666661, -2.96875, -2.427083333333333, -1.6979166666666665, -1.15625, -0.11458333333333333, -1.8854166666666667, -2.8437499999999982, 2.1979166666666665, 3.427083333333333, 4.96875, 6.5104166666666679, -3.4711538461538463, -2.2596153846153846, -2.0480769230769234, -1.1649408284023668, -0.55399408284023677, 0.52440828402366868, -0.089497041420118328, 1.1575443786982249, 1.2004437869822486, 1.8528106508875741, 1.8690828402366864, 2.9829881656804735, -3.9447115384615388, -2.4158653846153846, -2.3870192307692308, -2.5323594674556222, -6.0981878698224854, -3.3711168639053253, -0.54308431952662761, 0.66623520710059292, 2.2128328402366861, 4.1444156804733732, 5.9306582840236679, 8.3382026627218959, -4.1009615384615392, -2.9471153846153846, -2.2932692307692308, -1.6697485207100591, -0.61649408284023677, 0.40421597633136108, 0.030695266272189423, 1.2200443786982247, 1.7052514792899409, 2.0980029585798814, 2.5565828402366861, 3.6127958579881665), 12, 4), result_min_e = -2, result_max_e = 1, result_balance_e = NULL),
  regular_dynamic_balance2 = list(input = "regular", type = "dynamic", balance_e = 2, min_e = -Inf, max_e = Inf, status = "ok", warnings = c(), bin_ids = list(c("g4_t2"), c("g3_t2", "g4_t3"), c("g3_t3", "g4_t4"), c("g3_t4", "g4_t5"), c("g3_t5", "g4_t6")), eligible_cohorts = c(3, 4), selected_cell_ids = c("g3_t2", "g3_t3", "g3_t4", "g3_t5", "g4_t2", "g4_t3", "g4_t4", "g4_t5", "g4_t6"), overall_att = 5.1071428571428577, overall_se = 0.60188215317043259, egt = c(-2, -1, 0, 1, 2), att_egt = c(7.5, -0.10714285714285714, 7.7142857142857135, -1.6071428571428572, 9.2142857142857153), se_egt = c(0.45497705795335319, 0.31054649953338292, 0.72305166518774322, 0.42240250401275664, 0.80173648357935978), crit = 1.9599639845400534, overall_if = c(-2.6577380952380949, -1.5982142857142856, -1.5386904761904761, -1.4460034013605438, -2.947704081632653, -1.3677721088435366, 1.194302721088436, 3.6415816326530615, 0.92559523809523803, 1.5386904761904761, 1.5982142857142856, 2.6577380952380949), bin_if = matrix(c(-2.71875, -1.84375, -1.46875, -1.15625, -0.78125, 0.09375, -0.09375, 0.78125, 1.15625, 1.46875, 1.84375, 2.71875, -1.9999999999999998, -1.0357142857142856, -1.0714285714285714, -0.55102040816326525, -0.31122448979591838, 0.56122448979591832, -0.53061224489795922, 0.18877551020408168, 0.64285714285714279, 1.0714285714285714, 1.0357142857142856, 1.9999999999999998, -2.4241071428571428, -1.1919642857142856, -1.4598214285714284, -1.4559948979591835, -4.2034438775510203, -2.1447704081632639, 1.8692602040816333, 5.3054846938775508, 0.6294642857142857, 1.4598214285714284, 1.1919642857142856, 2.4241071428571428, -2.5803571428571428, -1.7232142857142856, -1.3660714285714286, -1.006377551020408, -0.37372448979591838, 0.39158163265306073, -0.36096938775510234, 0.25127551020408134, 1.0982142857142856, 1.3660714285714286, 1.7232142857142856, 2.5803571428571428, -2.9687499999999996, -1.8794642857142856, -1.7901785714285714, -1.8756377551020407, -4.2659438775510203, -2.3501275510204067, 2.0746173469387772, 5.3679846938775526, 1.0491071428571428, 1.7901785714285714, 1.8794642857142856, 2.9687499999999996), 12, 5), result_min_e = -Inf, result_max_e = Inf, result_balance_e = 2),
  regular_dynamic_positive = list(input = "regular", type = "dynamic", balance_e = NULL, min_e = 1, max_e = 2, status = "ok", warnings = c(), bin_ids = list(c("g3_t4", "g4_t5", "g5_t6"), c("g3_t5", "g4_t6")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t4", "g3_t5", "g4_t5", "g4_t6", "g5_t6"), overall_att = 3.7129120879120885, overall_se = 0.67659270675095962, egt = c(1, 2), att_egt = c(-1.7884615384615385, 9.2142857142857153), se_egt = c(0.65870945618837951, 0.80173648357935978), crit = 1.9599639845400534, overall_if = c(-3.5348557692307692, -2.4132898351648349, -2.041723901098901, -1.77269313790605, -2.4412189801956288, -0.97295578734452282, 1.0526563066054833, 3.2940145362878885, 1.377179311073542, 1.9440907650042263, 2.2180235629754859, 3.290772928994083), bin_if = matrix(c(-4.1009615384615392, -2.9471153846153846, -2.2932692307692308, -1.6697485207100591, -0.61649408284023677, 0.40421597633136108, 0.030695266272189423, 1.2200443786982247, 1.7052514792899409, 2.0980029585798814, 2.5565828402366861, 3.6127958579881665, -2.9687499999999996, -1.8794642857142856, -1.7901785714285714, -1.8756377551020407, -4.2659438775510203, -2.3501275510204067, 2.0746173469387772, 5.3679846938775526, 1.0491071428571428, 1.7901785714285714, 1.8794642857142856, 2.9687499999999996), 12, 2), result_min_e = 1, result_max_e = 2, result_balance_e = NULL),
  regular_dynamic_negative_only = list(input = "regular", type = "dynamic", balance_e = NULL, min_e = -3, max_e = -1, status = "refused", warnings = c(), bin_ids = list(c("g5_t2"), c("g4_t2", "g5_t3"), c("g3_t2", "g4_t3", "g5_t4")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t2", "g4_t2", "g4_t3", "g5_t2", "g5_t3", "g5_t4"), error = "No valid att_gt() estimates found for this aggregation. This may happen if all estimates for a particular group or time period are NA."),
  regular_dynamic_empty = list(input = "regular", type = "dynamic", balance_e = NULL, min_e = 20, max_e = 30, status = "refused", warnings = c(), bin_ids = list(), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c(), error = "No event times fall within the requested window. Adjust 'min_e'/'max_e' (and 'balance_e') so at least one event time is included."),
  regular_calendar = list(input = "regular", type = "calendar", balance_e = NULL, min_e = -Inf, max_e = Inf, status = "ok", warnings = c(), bin_ids = list(c("g3_t3"), c("g3_t4", "g4_t4"), c("g3_t5", "g4_t5", "g5_t5"), c("g3_t6", "g4_t6", "g5_t6")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t3", "g3_t4", "g3_t5", "g3_t6", "g4_t4", "g4_t5", "g4_t6", "g5_t5", "g5_t6"), overall_att = 4.4134615384615383, overall_se = 1.0250853260708908, egt = c(3, 4, 5, 6), att_egt = c(6.75, 2.9999999999999996, 7.3269230769230775, 0.57692307692307709), se_egt = c(0.24155663802719229, 2.6296801854057446, 2.1517641764497801, 2.0576336937951307), crit = 1.9599639845400534, overall_if = c(-3.174021291208792, -1.9484031593406594, -1.8395432692307692, -2.4010790891800511, -6.1037153951817409, -3.6075985312764161, 2.5627971787827559, 6.5551101542688084, 1.2859058009298392, 2.1501941568047336, 2.5697049344885881, 3.9506485101437026), bin_if = matrix(c(-1.59375, -0.46875, -0.84375, -0.28125, -0.65625, 0.46875, -0.46875, 0.65625, 0.28125, 0.84375, 0.46875, 1.59375, -2.5133928571428572, -1.4955357142857142, -1.40625, -4.7544642857142856, -16.236607142857142, -11.361607142857142, 10.075892857142858, 21.379464285714285, 0.89732142857142849, 1.40625, 1.4955357142857142, 2.5133928571428572, -4.1153846153846159, -2.75, -2.5, -1.1124260355029594, 0.41568047337278058, 1.3550295857988164, -8.4556213017751478, -15.468934911242604, 3.3491124260355027, 6.1213017751479288, 9.9926035502958577, 13.168639053254438, -4.4735576923076925, -3.0793269230769234, -2.6081730769230771, -3.4561760355029589, -7.937684911242604, -4.8925665680473367, 9.0996671597633139, 19.653661242603555, 0.61593934911242609, 0.22947485207100593, -1.678069526627219, -1.4731878698224854), 12, 4), result_min_e = NULL, result_max_e = NULL, result_balance_e = NULL),
  regular_calendar_ignored_selectors = list(input = "regular", type = "calendar", balance_e = 2, min_e = 1, max_e = 1, status = "ok", warnings = c("`min_e`, `max_e`, and `balance_e` are ignored for type = \"calendar\"; returning the unrestricted calendar-time effects."), bin_ids = list(c("g3_t3"), c("g3_t4", "g4_t4"), c("g3_t5", "g4_t5", "g5_t5"), c("g3_t6", "g4_t6", "g5_t6")), eligible_cohorts = c(3, 4, 5), selected_cell_ids = c("g3_t3", "g3_t4", "g3_t5", "g3_t6", "g4_t4", "g4_t5", "g4_t6", "g5_t5", "g5_t6"), overall_att = 4.4134615384615383, overall_se = 1.0250853260708908, egt = c(3, 4, 5, 6), att_egt = c(6.75, 2.9999999999999996, 7.3269230769230775, 0.57692307692307709), se_egt = c(0.24155663802719229, 2.6296801854057446, 2.1517641764497801, 2.0576336937951307), crit = 1.9599639845400534, overall_if = c(-3.174021291208792, -1.9484031593406594, -1.8395432692307692, -2.4010790891800511, -6.1037153951817409, -3.6075985312764161, 2.5627971787827559, 6.5551101542688084, 1.2859058009298392, 2.1501941568047336, 2.5697049344885881, 3.9506485101437026), bin_if = matrix(c(-1.59375, -0.46875, -0.84375, -0.28125, -0.65625, 0.46875, -0.46875, 0.65625, 0.28125, 0.84375, 0.46875, 1.59375, -2.5133928571428572, -1.4955357142857142, -1.40625, -4.7544642857142856, -16.236607142857142, -11.361607142857142, 10.075892857142858, 21.379464285714285, 0.89732142857142849, 1.40625, 1.4955357142857142, 2.5133928571428572, -4.1153846153846159, -2.75, -2.5, -1.1124260355029594, 0.41568047337278058, 1.3550295857988164, -8.4556213017751478, -15.468934911242604, 3.3491124260355027, 6.1213017751479288, 9.9926035502958577, 13.168639053254438, -4.4735576923076925, -3.0793269230769234, -2.6081730769230771, -3.4561760355029589, -7.937684911242604, -4.8925665680473367, 9.0996671597633139, 19.653661242603555, 0.61593934911242609, 0.22947485207100593, -1.678069526627219, -1.4731878698224854), 12, 4), result_min_e = NULL, result_max_e = NULL, result_balance_e = NULL),
  irregular_simple_max1 = list(input = "irregular", type = "simple", balance_e = NULL, min_e = -Inf, max_e = 1, status = "ok", warnings = c(), bin_ids = list(c("g2005_t2005", "g2005_t2009", "g2009_t2009", "g2009_t2010", "g2010_t2010", "g2010_t2014")), eligible_cohorts = c(2005, 2009, 2010), selected_cell_ids = c("g2005_t2005", "g2005_t2009", "g2009_t2009", "g2009_t2010", "g2010_t2010", "g2010_t2014"), overall_att = 3.7788461538461537, overall_se = 0.89359806418597587, egt = NULL, att_egt = NULL, se_egt = NULL, crit = NULL, overall_if = c(-4.0228365384615392, -2.6814903846153846, -2.3401442307692308, -2.1010539940828408, -3.3573409763313613, -1.483450443786982, -0.25619452662721942, 0.9431397928994093, 1.9590421597633134, 3.1212093195266273, 4.243620562130177, 5.975499260355031), bin_if = NULL, result_min_e = NULL, result_max_e = NULL, result_balance_e = NULL),
  irregular_group_max1 = list(input = "irregular", type = "group", balance_e = NULL, min_e = -Inf, max_e = 1, status = "ok", warnings = c(), bin_ids = list(c("g2005_t2005", "g2005_t2009"), c("g2009_t2009", "g2009_t2010"), c("g2010_t2010", "g2010_t2014")), eligible_cohorts = c(2005, 2009, 2010), selected_cell_ids = c("g2005_t2005", "g2005_t2009", "g2009_t2009", "g2009_t2010", "g2010_t2010", "g2010_t2014"), overall_att = 3.7788461538461537, overall_se = 0.89359806418597587, egt = c(2005, 2009, 2010), att_egt = c(2.625, 3.625, 4.625), se_egt = c(0.25546065332905099, 0.6118350894230451, 0.99174312008334042), crit = 1.9599639845400534, overall_if = c(-4.0228365384615383, -2.6814903846153846, -2.3401442307692308, -2.1010539940828408, -3.3573409763313613, -1.4834504437869827, -0.2561945266272192, 0.94313979289940875, 1.9590421597633134, 3.1212093195266273, 4.243620562130177, 5.9754992603550292), bin_if = matrix(c(-1.671875, -0.734375, -0.796875, -0.515625, -0.578125, 0.359375, -0.359375, 0.578125, 0.515625, 0.796875, 0.734375, 1.671875, -3.609375, -2.421875, -2.234375, -1.328125, -1.140625, 0.046875, -0.046875, 1.140625, 1.328125, 2.234375, 2.421875, 3.609375, -5.796875, -4.109375, -3.421875, -2.390625, -1.703125, -0.015625, 0.015625, 1.703125, 2.390625, 3.421875, 4.109375, 5.796875), 12, 3), result_min_e = NULL, result_max_e = NULL, result_balance_e = NULL),
  irregular_dynamic = list(input = "irregular", type = "dynamic", balance_e = NULL, min_e = -Inf, max_e = Inf, status = "ok", warnings = c(), bin_ids = list(c("g2010_t2002"), c("g2009_t2002"), c("g2010_t2005"), c("g2009_t2005"), c("g2005_t2002"), c("g2010_t2009"), c("g2005_t2005", "g2009_t2009", "g2010_t2010"), c("g2009_t2010"), c("g2005_t2009", "g2010_t2014"), c("g2005_t2010", "g2009_t2014"), c("g2005_t2014")), eligible_cohorts = c(2005, 2009, 2010), selected_cell_ids = c("g2005_t2002", "g2005_t2005", "g2005_t2009", "g2005_t2010", "g2005_t2014", "g2009_t2002", "g2009_t2005", "g2009_t2009", "g2009_t2010", "g2009_t2014", "g2010_t2002", "g2010_t2005", "g2010_t2009", "g2010_t2010", "g2010_t2014"), overall_att = 2.4020879120879122, overall_se = 0.69247041423248001, egt = c(-8, -7, -5, -4, -3, -1, 0, 1, 4, 5, 9), att_egt = c(1, 7.5, 9.75, -0.25, 0, -0.5, 9.3461538461538467, -1.75, -1.8, 9.2142857142857153, -3), se_egt = c(0.78616509433805026, 0.45497705795335319, 0.81926603224912375, 0.5303300858899106, 0.16103775868479486, 0.89261777545350041, 1.2063781190787524, 0.63661089960823003, 0.66787546660203723, 0.80173648357935978, 0.41129347963332574), crit = 1.9599639845400534, overall_if = c(-3.4776923076923079, -2.2640659340659344, -2.0504395604395604, -1.6705994445115326, -2.5638263494747013, -0.98124888298514645, 0.25130660548242995, 1.8418439801956292, 1.4533879966187659, 2.3539188503803889, 2.871024513947591, 4.2363905325443794), bin_if = matrix(c(-4.75, -3, -2.75, -1.75, -1.5, 0.25, -0.25, 1.5, 1.75, 2.75, 3, 4.75, -2.71875, -1.84375, -1.46875, -1.15625, -0.78125, 0.09375, -0.09375, 0.78125, 1.15625, 1.46875, 1.84375, 2.71875, -4.65625, -3.53125, -2.90625, -1.96875, -1.34375, -0.21875, 0.21875, 1.34375, 1.96875, 2.90625, 3.53125, 4.65625, -3.25, -2, -1.75, -1.25, -1, 0.25, -0.25, 1, 1.25, 1.75, 2, 3.25, -1.0625, -0.3125, -0.5625, -0.1875, -0.4375, 0.3125, -0.3125, 0.4375, 0.1875, 0.5625, 0.3125, 1.0625, -5.1875, -3.6875, -3.1875, -2.0625, -1.5625, -0.0625, 0.0625, 1.5625, 2.0625, 3.1875, 3.6875, 5.1875, -3.9447115384615388, -2.4158653846153846, -2.3870192307692308, -2.5323594674556222, -6.0981878698224854, -3.3711168639053253, -0.54308431952662761, 0.66623520710059292, 2.2128328402366861, 4.1444156804733732, 5.9306582840236679, 8.3382026627218959, -3.6875, -2.6875, -2.1875, -1.5625, -1.0625, -0.0625, 0.0625, 1.0625, 1.5625, 2.1875, 2.6875, 3.6875, -4.2249999999999996, -3.0249999999999999, -2.3250000000000002, -1.6950000000000001, -0.45500000000000029, 0.5650000000000005, -0.025000000000000008, 1.175, 1.7550000000000001, 2.0850000000000004, 2.5449999999999999, 3.625, -2.9687499999999996, -1.8794642857142856, -1.7901785714285714, -1.8756377551020407, -4.2659438775510203, -2.3501275510204067, 2.0746173469387772, 5.3679846938775526, 1.0491071428571428, 1.7901785714285714, 1.8794642857142856, 2.9687499999999996, -2.5625, -1.3125, -1.5625, -0.6875, -0.9375, 0.3125, -0.3125, 0.9375, 0.6875, 1.5625, 1.3125, 2.5625), 12, 11), result_min_e = -Inf, result_max_e = Inf, result_balance_e = NULL),
  irregular_dynamic_balance5_window = list(input = "irregular", type = "dynamic", balance_e = 5, min_e = -4, max_e = 5, status = "ok", warnings = c(), bin_ids = list(c("g2009_t2005"), c("g2005_t2002"), c("g2005_t2005", "g2009_t2009"), c("g2009_t2010"), c("g2005_t2009"), c("g2005_t2010", "g2009_t2014")), eligible_cohorts = c(2005, 2009), selected_cell_ids = c("g2005_t2002", "g2005_t2005", "g2005_t2009", "g2005_t2010", "g2009_t2005", "g2009_t2009", "g2009_t2010", "g2009_t2014"), overall_att = 3.4196428571428568, overall_se = 0.56417053939077577, egt = c(-4, -3, 0, 1, 4, 5), att_egt = c(-0.25, 0, 7.7142857142857135, -1.75, -1.5, 9.2142857142857153), se_egt = c(0.5303300858899106, 0.16103775868479486, 0.72305166518774322, 0.63661089960823003, 0.27638539919628335, 0.80173648357935978), crit = 1.9599639845400534, overall_if = c(-2.7075892857142856, -1.6897321428571428, -1.546875, -1.4110331632653059, -2.5079719387755102, -1.0768494897959178, 0.93909438775510257, 3.0589923469387759, 0.9977678571428571, 1.546875, 1.6897321428571428, 2.7075892857142856), bin_if = matrix(c(-3.25, -2, -1.75, -1.25, -1, 0.25, -0.25, 1, 1.25, 1.75, 2, 3.25, -1.0625, -0.3125, -0.5625, -0.1875, -0.4375, 0.3125, -0.3125, 0.4375, 0.1875, 0.5625, 0.3125, 1.0625, -2.4241071428571428, -1.1919642857142856, -1.4598214285714284, -1.4559948979591835, -4.2034438775510203, -2.1447704081632639, 1.8692602040816333, 5.3054846938775508, 0.6294642857142857, 1.4598214285714284, 1.1919642857142856, 2.4241071428571428, -3.6875, -2.6875, -2.1875, -1.5625, -1.0625, -0.0625, 0.0625, 1.0625, 1.5625, 2.1875, 2.6875, 3.6875, -1.75, -1, -0.75, -0.75, -0.5, 0.25, -0.25, 0.5, 0.75, 0.75, 1, 1.75, -2.9687499999999996, -1.8794642857142856, -1.7901785714285714, -1.8756377551020407, -4.2659438775510203, -2.3501275510204067, 2.0746173469387772, 5.3679846938775526, 1.0491071428571428, 1.7901785714285714, 1.8794642857142856, 2.9687499999999996), 12, 6), result_min_e = -4, result_max_e = 5, result_balance_e = 5)
)

# ==============================================================================
# A. Pure-R composition checks
# ==============================================================================

test_that("released selectors reproduce the executed did 2.5.1 oracle: sets, values and complete IFs", {
  oi <- oracle_inputs(); oe <- oracle_expectations()
  tm_reg <- data.frame(g = oi$cell_group, time = oi$cell_t, target_id = oi$cell_id,
    role = ifelse(oi$cell_t >= oi$cell_group, "post", "pre"), support_status = "estimable", stringsAsFactors = FALSE)
  tm_irr <- data.frame(g = oi$irregular_cell_group, time = oi$irregular_cell_t, target_id = oi$irregular_cell_id,
    role = ifelse(oi$irregular_cell_t >= oi$irregular_cell_group, "post", "pre"), support_status = "estimable", stringsAsFactors = FALSE)
  for (nm in names(oe)) {
    e <- oe[[nm]]
    irregular <- identical(e$input, "irregular")
    tm <- if (irregular) tm_irr else tm_reg
    coh <- if (irregular) oi$irregular_cohort else oi$frame_cohort
    fg <- coh; fg[fg == 0] <- Inf                                  # public 0 -> internal Inf
    att <- oi$att; names(att) <- tm$target_id
    IF <- oi$inffunc; dimnames(IF) <- NULL
    cal <- if (irregular) oi$times_irregular else oi$times_regular
    run <- function() aggregate_chained_ipw_targets(att, IF, tm, fg, oi$frame_weight, type = e$type,
      multiplier_budget = req_budget(), inference_accuracy = req_accuracy(), balance_e = e$balance_e,
      min_e = e$min_e, max_e = e$max_e, calendar = cal, cohort_grid = unique(coh[coh != 0]),
      bstrap = FALSE, caller_rng = TRUE)
    if (identical(nm, "regular_dynamic_empty")) {
      # Upstream refuses an empty event window; so does the chained binding (typed).
      expect_error(run(), "No event times", class = "cdid_ipw_pipeline_selector_error")
      next
    }
    if (identical(nm, "regular_dynamic_negative_only")) {
      # Explicit selector choices preserve the released whole-call refusal;
      # D5 support preservation does not redefine an intentionally empty overall.
      expect_error(run(), "No valid att_gt", class = "cdid_ipw_pipeline_selector_error")
      next
    }
    cw <- capture_warnings(run()); a <- cw$value
    expect_identical(cw$warnings, as.character(e$warnings), info = nm)     # calendar warning only
    expect_identical(sort(a$selectors$selected_cell_ids), sort(e$selected_cell_ids), info = nm)
    expect_identical(as.double(a$selectors$eligible_cohorts), e$eligible_cohorts, info = nm)
    expect_identical(unname(lapply(a$cell_target_ids, sort)), unname(lapply(e$bin_ids, sort)), info = nm)
    nondefault <- !is.null(e$balance_e) || is.finite(e$min_e) || is.finite(e$max_e)
    expect_identical(a$selectors$applied, nondefault && !identical(e$type, "calendar"), info = nm)
    expect_identical(a$selectors$ignored, nondefault && identical(e$type, "calendar"), info = nm)
    expect_identical(e$status, "ok")
    expect_equal(a$overall$estimate, e$overall_att, tolerance = 1e-12, info = nm)
    expect_equal(a$overall$se, e$overall_se, tolerance = 1e-12, info = nm)
    expect_equal(unname(a$overall$if_full), e$overall_if, tolerance = 1e-12, info = nm)
    if (!identical(e$type, "simple")) {
      expect_identical(as.double(a$cell_keys), e$egt, info = nm)
      expect_equal(unname(a$estimates), e$att_egt, tolerance = 1e-12, info = nm)
      expect_equal(unname(a$se), e$se_egt, tolerance = 1e-12, info = nm)
      expect_equal(unname(a$if_full), unname(e$bin_if), tolerance = 1e-12, info = nm)
      expect_equal(a$cband_crit_val, e$crit, tolerance = 1e-12, info = nm)
    }
    expect_identical(a$inference$bootstrap, FALSE)
  }
  # Analytic aggregation never touches the caller RNG.
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  invisible(aggregate_chained_ipw_targets(stats::setNames(oi$att, tm_reg$target_id), unname(oi$inffunc), tm_reg,
    replace(oi$frame_cohort, oi$frame_cohort == 0, Inf), oi$frame_weight, type = "dynamic",
    multiplier_budget = req_budget(), inference_accuracy = req_accuracy(), bstrap = FALSE, caller_rng = TRUE))
  expect_false(exists(".Random.seed", envir = globalenv()))
  # Selector validation mirrors upstream's: nonnumeric or negative balance_e refused.
  expect_error(aggregate_chained_ipw_targets(stats::setNames(oi$att, tm_reg$target_id), unname(oi$inffunc), tm_reg,
    replace(oi$frame_cohort, oi$frame_cohort == 0, Inf), oi$frame_weight, type = "dynamic", balance_e = -1,
    multiplier_budget = req_budget(), inference_accuracy = req_accuracy()), class = "cdid_ipw_pipeline_selector_error")
})

# Ten units in five clusters (consecutive pairs) with distinct nonzero cluster
# sums, so every multiplier marginal has a positive empirical IQR scale.
req_fixed_if <- function() {
  X <- rbind(c(-2, 4), c(1, 1), c(-1, 2), c(1, -3), c(1, -2), c(2, 1), c(2, -4), c(-1, 2), c(0.5, 1), c(-3.5, -3))
  colnames(X) <- c("a", "b"); X
}

test_that("bstrap and cband are separate flags with the released meanings; RNG laws are honest", {
  X <- req_fixed_if(); th <- c(a = 0.5, b = -0.2); cl <- rep(1:5, each = 2)
  run <- function(...) compute_chained_ipw_inference(th, X, cl, 0.05, biters = 32L, multiplier_budget = req_budget(),
    inference_accuracy = req_accuracy(), ...)
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  # (FALSE, FALSE): analytic SE, pointwise constant, no draw, RNG untouched
  ff <- run(cband = FALSE, bstrap = FALSE, caller_rng = TRUE)
  expect_false(ff$bootstrap); expect_null(ff$multiplier); expect_equal(ff$cband_crit_val, stats::qnorm(0.975))
  expect_identical(unname(ff$se), unname(ff$analytical_se)); expect_identical(ff$band_family, "pointwise_normal_over_certified_ids")
  expect_false(exists(".Random.seed", envir = globalenv()))
  # (TRUE, FALSE): multiplier IQR SE, pointwise normal constant, caller RNG advanced
  set.seed(5); s0 <- .Random.seed
  tf <- run(cband = FALSE, bstrap = TRUE, caller_rng = TRUE); s1 <- .Random.seed
  expect_true(tf$bootstrap); expect_false(tf$simultaneous); expect_false(identical(s0, s1))
  expect_equal(tf$cband_crit_val, stats::qnorm(0.975)); expect_false(isTRUE(all.equal(unname(tf$se), unname(tf$analytical_se))))
  expect_identical(tf$band_family, "pointwise_normal_over_certified_ids"); expect_identical(tf$rng$mode, "caller_rng_stream_Rademacher")
  expect_false(isTRUE(tf$rng$caller_rng_restored)); expect_identical(tf$n_draws, 32L)
  # (TRUE, TRUE): same draws from the same caller seed -> same SEs, sup-t constant
  set.seed(5); tt <- run(cband = TRUE, bstrap = TRUE, caller_rng = TRUE)
  expect_true(tt$simultaneous); expect_equal(tt$se, tf$se); expect_false(isTRUE(all.equal(tt$cband_crit_val, stats::qnorm(0.975))))
  expect_identical(tt$band_family, "simultaneous_multiplier_sup_t_over_certified_ids")
  # the same caller seed reproduces the public-style draws
  set.seed(5); tt2 <- run(cband = TRUE, bstrap = TRUE, caller_rng = TRUE)
  expect_identical(tt2$se, tt$se); expect_identical(tt2$cband_crit_val, tt$cband_crit_val)
  # (FALSE, TRUE) at the att_gt-style endpoint (released att_gt.R:582-590,688-695):
  # analytic SEs, normal pointwise constant, NO draw, caller RNG untouched.
  set.seed(5); s0 <- .Random.seed
  ft <- run(cband = TRUE, bstrap = FALSE, caller_rng = TRUE)
  expect_identical(.Random.seed, s0)
  expect_false(ft$bootstrap); expect_false(ft$simultaneous); expect_null(ft$multiplier); expect_null(ft$rng)
  expect_identical(ft$se, ff$se); expect_identical(ft$analytical_se, ff$analytical_se)
  expect_identical(ft$cband_crit_val, ff$cband_crit_val); expect_identical(ft$cband_lower, ff$cband_lower)
  expect_identical(ft$vcov, ff$vcov); expect_identical(ft$band_family, "pointwise_normal_over_certified_ids")
  expect_identical(ft$flag_law, "released_bstrap_false_cband_true_cells_analytic_pointwise_no_draws")
  expect_true(is.na(ft$n_draws))
  # (FALSE, TRUE) at the aggte-bin endpoint (released compute.aggte.R:365-370 etc.):
  # analytic reported SEs (getSE), multiplier sup-t critical value from the
  # caller's stream (same draws as TRUE/TRUE from the same seed), RNG advanced.
  set.seed(5); fb <- run(cband = TRUE, bstrap = FALSE, caller_rng = TRUE, endpoint = "aggte_bins"); s_fb <- .Random.seed
  set.seed(5); tt3 <- run(cband = TRUE, bstrap = TRUE, caller_rng = TRUE); s_tt <- .Random.seed
  expect_identical(s_fb, s_tt); expect_false(identical(s_fb, s0))             # same stream consumption as TRUE/TRUE
  expect_false(fb$bootstrap); expect_true(fb$simultaneous); expect_false(is.null(fb$multiplier))
  expect_identical(fb$se, ff$se)                                               # analytic, not the IQR SEs
  expect_false(isTRUE(all.equal(fb$se, tt3$se)))
  expect_identical(fb$cband_crit_val, tt3$cband_crit_val)                      # the multiplier's own sup-t quantile
  expect_identical(fb$cband_crit_val, fb$multiplier$crit_val)
  expect_identical(unname(fb$cband_lower), unname(ff$estimates - fb$cband_crit_val * ff$analytical_se))
  expect_identical(unname(fb$cband_upper), unname(ff$estimates + fb$cband_crit_val * ff$analytical_se))
  expect_identical(fb$ci_lower, ff$ci_lower)                                   # pointwise interval stays analytic
  expect_identical(fb$band_family, "simultaneous_multiplier_sup_t_over_certified_ids")
  expect_identical(fb$flag_law, "released_bstrap_false_cband_true_aggte_bins_analytic_se_multiplier_sup_t_band")
  expect_identical(fb$n_draws, 32L); expect_identical(fb$rng$mode, "caller_rng_stream_Rademacher")
  # released mboot convention: draws standardized by their own IQR scale (type-1 quantiles), never by the reported SE
  expect_identical(fb$multiplier$metadata$quantile_type, 1L)
  expect_equal(fb$multiplier$metadata$iqr_normalizer, stats::qnorm(.75) - stats::qnorm(.25))
  expect_identical(fb$multiplier$bootstrap_iqr_scale, tt3$multiplier$bootstrap_iqr_scale)
  expect_false(isTRUE(all.equal(unname(fb$multiplier$bootstrap_iqr_scale), unname(fb$se))))
  # existing private contract (bstrap = NULL): bootstrap iff cband, explicit seed, caller RNG restored
  set.seed(7); s0 <- .Random.seed
  leg <- run(cband = TRUE, seed = 99L)
  expect_identical(.Random.seed, s0); expect_identical(leg$rng$mode, "R_sample_Rademacher"); expect_identical(leg$flag_law, "private_bootstrap_iff_cband")
  expect_true(leg$bootstrap); expect_true(leg$simultaneous)
  leg0 <- run(cband = FALSE, seed = 99L); expect_false(leg0$bootstrap)
  # explicit flags with the same explicit seed are the identical computation
  exp_ <- run(cband = TRUE, bstrap = TRUE, seed = 99L)
  expect_identical(exp_$se, leg$se); expect_identical(exp_$cband_crit_val, leg$cband_crit_val)
  expect_identical(exp_$cband_lower, leg$cband_lower)
  expect_identical(exp_$flag_law, "explicit_bstrap_and_cband")
  # no silent seed: generated signs need a seed or the caller stream; both together refused
  expect_error(run(cband = TRUE, bstrap = TRUE), class = "cdid_ipw_pipeline_bootstrap_error")
  expect_error(run(cband = TRUE, bstrap = TRUE, seed = 1L, caller_rng = TRUE), class = "cdid_ipw_pipeline_bootstrap_error")
  # exact-reference-only selection: no draw under any bootstrap request
  Z <- matrix(0, 10, 1, dimnames = list(NULL, "ref")); rm(".Random.seed", envir = globalenv())
  det <- compute_chained_ipw_inference(c(ref = 0), Z, cl, 0.05, cband = TRUE, bstrap = TRUE, caller_rng = TRUE, biters = 32L,
    multiplier_budget = req_budget(), inference_accuracy = req_accuracy())
  expect_false(det$bootstrap); expect_identical(det$cband_crit_val, 0); expect_false(exists(".Random.seed", envir = globalenv()))
  expect_identical(det$band_family, "deterministic_exact_zero_over_certified_ids")
})

test_that("fixed signs give matched original-vs-revised inference across flag forms and use no RNG", {
  X <- req_fixed_if(); th <- c(a = 0.5, b = -0.2); cl <- rep(1:5, each = 2)
  signs <- unname(as.matrix(expand.grid(rep(list(c(-1.0, 1.0)), 5L)))); storage.mode(signs) <- "double"
  enc <- list(theta = list(lo = matrix(th, 1L), hi = matrix(th, 1L)), rows = list(lo = X, hi = X))
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  args <- list(th, X, cl, 0.05, biters = 32L, multiplier_budget = req_budget(), inference_accuracy = req_accuracy(),
    fixed_signs = signs, enclosure = enc)
  original <- do.call(compute_chained_ipw_inference, c(args, list(cband = TRUE)))
  revised <- do.call(compute_chained_ipw_inference, c(args, list(cband = TRUE, bstrap = TRUE)))
  pointwise <- do.call(compute_chained_ipw_inference, c(args, list(cband = FALSE, bstrap = TRUE)))
  expect_false(exists(".Random.seed", envir = globalenv()))
  expect_identical(revised$se, original$se); expect_identical(revised$cband_crit_val, original$cband_crit_val)
  expect_identical(revised$numerical_certificate$errors, original$numerical_certificate$errors)
  expect_identical(pointwise$se, original$se)                       # same IQR SEs
  expect_equal(pointwise$cband_crit_val, stats::qnorm(0.975))       # pointwise constant
  expect_true(pointwise$numerical_certificate$certified)
  expect_identical(pointwise$numerical_certificate$critical_value, c(stats::qnorm(0.975), stats::qnorm(0.975)))
  expect_true(all(pointwise$numerical_certificate$se == original$numerical_certificate$se))
  expect_identical(pointwise$rng$mode, "caller_fixed_signs")
  # Released aggte-bin FALSE/TRUE law with the same fixed signs, checked against
  # the supplied upstream formulas (mboot.R:157-188, getSE 811-835): analytic
  # cluster-sum SEs; crit = type-1 (1 - alp) order statistic of max_j |b_j| / s_j
  # with s_j = (q75_j - q25_j) / (qnorm(.75) - qnorm(.25)) from the draws
  # themselves; band = est +/- crit * analytic SE.
  analytic <- do.call(compute_chained_ipw_inference, c(args, list(cband = FALSE, bstrap = FALSE)))
  band_only <- do.call(compute_chained_ipw_inference, c(args, list(cband = TRUE, bstrap = FALSE, endpoint = "aggte_bins")))
  expect_false(exists(".Random.seed", envir = globalenv()))
  S <- rowsum(X, cl); Sc <- sweep(S, 2L, colMeans(S), "-")
  raw <- signs %*% Sc / sqrt(5)
  scale <- apply(raw, 2L, function(b) diff(stats::quantile(b, c(.25, .75), type = 1, names = FALSE))) / diff(stats::qnorm(c(.25, .75)))
  crit <- stats::quantile(apply(abs(sweep(raw, 2L, scale, "/")), 1L, max), .95, type = 1, names = FALSE)
  expect_equal(unname(band_only$se), unname(sqrt(colSums(Sc^2)) / 10), tolerance = 1e-15)
  expect_identical(band_only$se, analytic$se); expect_identical(band_only$analytical_se, analytic$analytical_se)
  expect_equal(band_only$cband_crit_val, crit, tolerance = 1e-15)
  expect_identical(band_only$cband_crit_val, original$cband_crit_val)        # same sup-t statistic as the band under TRUE/TRUE
  expect_false(isTRUE(all.equal(band_only$se, original$se)))                 # but not the IQR SEs
  expect_equal(unname(band_only$cband_lower), unname(th - crit * band_only$se), tolerance = 1e-15)
  expect_equal(unname(band_only$cband_upper), unname(th + crit * band_only$se), tolerance = 1e-15)
  expect_false(band_only$bootstrap); expect_true(band_only$simultaneous); expect_identical(band_only$rng$mode, "caller_fixed_signs")
  # Mixed law certificate: reported SEs certified against their analytic law,
  # the multiplier scales kept apart, the sup-t constant and the band certified.
  cb <- band_only$numerical_certificate
  expect_true(cb$certified)
  expect_identical(cb$se, analytic$numerical_certificate$se)
  expect_identical(cb$critical_value, original$numerical_certificate$critical_value)
  # the retained multiplier scale intervals enclose mboot's IQR scales and are the
  # same intervals that TRUE/TRUE turns into its SE certificate (times sqrt(G)/n)
  expect_true(all(cb$multiplier_scale[1L, ] <= scale & scale <= cb$multiplier_scale[2L, ]))
  # Interval scaling must itself use outward-rounded operations. Ordinary
  # multiplication cannot invert or reproduce that operation bit for bit.
  pure <- new.env(parent = baseenv())
  sys.source(testthat::test_path("reference", "chained_gmm_core_pure_r_reference.R"), envir = pure)
  factor <- pure$gmm_iv_div(pure$gmm_iv_sqrt(pure$gmm_iv(5.0)), pure$gmm_iv(10.0))
  expected_se_interval <- vapply(seq_len(ncol(cb$multiplier_scale)), function(j)
    pure$gmm_iv_mul(cb$multiplier_scale[, j], factor), numeric(2L))
  expect_identical(expected_se_interval, original$numerical_certificate$se)
  expect_true(all(cb$band$lo <= band_only$cband_lower & band_only$cband_lower <= cb$band$hi))
  expect_true(all(cb$band_upper$lo <= band_only$cband_upper & band_only$cband_upper <= cb$band_upper$hi))
  expect_identical(names(cb$errors), names(analytic$numerical_certificate$errors))
  expect_identical(cb$errors[["max_multiplier_error_abs"]], original$numerical_certificate$errors[["max_multiplier_error_abs"]])
  expect_identical(cb$errors[["max_se_error_abs"]], analytic$numerical_certificate$errors[["max_se_error_abs"]])
  expect_null(analytic$numerical_certificate$multiplier_scale); expect_null(original$numerical_certificate$multiplier_scale)
})

test_that("aggregation threads the flags: bootstrap SEs for bins and overall, analytic fixed-share diagnostic", {
  oi <- oracle_inputs()
  tm <- data.frame(g = oi$cell_group, time = oi$cell_t, target_id = oi$cell_id,
    role = ifelse(oi$cell_t >= oi$cell_group, "post", "pre"), support_status = "estimable", stringsAsFactors = FALSE)
  fg <- replace(oi$frame_cohort, oi$frame_cohort == 0, Inf); att <- stats::setNames(oi$att, tm$target_id)
  run <- function(...) aggregate_chained_ipw_targets(att, unname(oi$inffunc), tm, fg, oi$frame_weight, type = "group",
    multiplier_budget = req_budget(), inference_accuracy = req_accuracy(), biters = 24L, ...)
  set.seed(3); a_tf <- run(bstrap = TRUE, cband = FALSE, caller_rng = TRUE)
  set.seed(3); a_tt <- run(bstrap = TRUE, cband = TRUE, caller_rng = TRUE)
  a_ff <- run(bstrap = FALSE, cband = FALSE, caller_rng = TRUE)
  expect_true(a_tf$inference$bootstrap); expect_false(a_tf$inference$simultaneous)
  expect_equal(a_tf$cband_crit_val, stats::qnorm(0.975)); expect_false(isTRUE(all.equal(a_tf$se, a_ff$se)))
  expect_identical(a_tf$se, a_tt$se)                                        # same caller draws
  expect_false(isTRUE(all.equal(a_tt$cband_crit_val, stats::qnorm(0.975))))
  expect_identical(a_tf$overall$inference_kind, "multiplier_bootstrap_iqr_se")
  expect_identical(a_ff$overall$inference_kind, "analytic")
  expect_false(isTRUE(all.equal(a_tf$overall$se, a_ff$overall$se)))
  expect_identical(a_tf$se_fixed_shares_only, a_ff$se_fixed_shares_only)      # diagnostic stays analytic
  expect_identical(a_tf$estimates, a_ff$estimates)
  # (FALSE, TRUE) group bins: released compute.aggte.R:365-370 warning, analytic
  # reported SEs and overall (getSE), band from the caller's draws (same seed 3 ->
  # same sup-t constant as TRUE/TRUE), RNG advanced for the bins only.
  set.seed(3); s0 <- .Random.seed
  cw <- capture_warnings(run(bstrap = FALSE, cband = TRUE, caller_rng = TRUE)); a_ft <- cw$value; s_ft <- .Random.seed
  expect_identical(cw$warnings, "Used bootstrap procedure to compute simultaneous confidence band")
  expect_false(identical(s_ft, s0))
  expect_false(a_ft$inference$bootstrap); expect_true(a_ft$inference$simultaneous)
  expect_identical(a_ft$se, a_ff$se); expect_identical(a_ft$overall$se, a_ff$overall$se)
  expect_identical(a_ft$overall$inference_kind, "analytic"); expect_identical(a_ft$overall$ci_lower, a_ff$overall$ci_lower)
  expect_identical(a_ft$cband_crit_val, a_tt$cband_crit_val)
  expect_identical(unname(a_ft$cband_lower), unname(a_ff$estimates - a_tt$cband_crit_val * a_ff$se))
  expect_identical(a_ft$se_fixed_shares_only, a_ff$se_fixed_shares_only); expect_identical(a_ft$estimates, a_ff$estimates)
  expect_identical(a_ft$inference$flag_law, "released_bstrap_false_cband_true_aggte_bins_analytic_se_multiplier_sup_t_band")
  expect_null(a_ft$numerical_certificate)   # no enclosure in this oracle run; certificates are checked with enclosures above
  # the bins consume exactly the draws TRUE/TRUE used for its bins (its overall draws a second set)
  set.seed(3); s0 <- .Random.seed
  bins_tt <- compute_chained_ipw_inference(a_tt$estimates, a_tt$if_full, NULL, 0.05, cband = TRUE, bstrap = TRUE, caller_rng = TRUE,
    biters = 24L, multiplier_budget = req_budget(), inference_accuracy = req_accuracy())
  expect_identical(.Random.seed, s_ft); expect_identical(bins_tt$cband_crit_val, a_ft$cband_crit_val)
  # flag_reference parity (actual host did 2.5.1, group FALSE/TRUE, same regular MP):
  # reported se.egt are the analytical oracle values; the band constant is a
  # finite sup-t quantile above the normal constant on this fixture.
  expect_equal(unname(a_ft$se), oracle_expectations()$regular_group$se_egt, tolerance = 1e-12)
  expect_true(is.finite(a_ft$cband_crit_val) && a_ft$cband_crit_val > 0)
  # simple aggregation has no bin-band family: FALSE/TRUE is analytic, no warning, no draw
  rm(".Random.seed", envir = globalenv())
  cw_s <- capture_warnings(aggregate_chained_ipw_targets(att, unname(oi$inffunc), tm, fg, oi$frame_weight, type = "simple",
    multiplier_budget = req_budget(), inference_accuracy = req_accuracy(), biters = 24L, bstrap = FALSE, cband = TRUE, caller_rng = TRUE))
  a_s <- cw_s$value
  expect_identical(cw_s$warnings, character()); expect_false(exists(".Random.seed", envir = globalenv()))
  expect_false(a_s$inference$bootstrap); expect_false(a_s$inference$simultaneous); expect_null(a_s$inference$multiplier)
  expect_identical(a_s$inference$flag_law, "released_bstrap_false_cband_true_cells_analytic_pointwise_no_draws")
  expect_equal(a_s$overall$se, oracle_expectations()$regular_simple$overall_se, tolerance = 1e-12)
  expect_equal(a_s$cband_crit_val, stats::qnorm(0.975))
  # dynamic and calendar bins warn and draw under the same law; explicit seed law too
  for (typ in c("dynamic", "calendar")) {
    set.seed(5); s0 <- .Random.seed
    cw_t <- capture_warnings(aggregate_chained_ipw_targets(att, unname(oi$inffunc), tm, fg, oi$frame_weight, type = typ,
      multiplier_budget = req_budget(), inference_accuracy = req_accuracy(), biters = 24L, bstrap = FALSE, cband = TRUE, caller_rng = TRUE))
    expect_identical(cw_t$warnings, "Used bootstrap procedure to compute simultaneous confidence band", info = typ)
    expect_false(identical(.Random.seed, s0), info = typ)
    expect_false(cw_t$value$inference$bootstrap, info = typ); expect_true(cw_t$value$inference$simultaneous, info = typ)
    expect_equal(unname(cw_t$value$se), oracle_expectations()[[paste0("regular_", typ)]]$se_egt, tolerance = 1e-12, info = typ)
    expect_identical(cw_t$value$overall$inference_kind, "analytic", info = typ)
  }
  set.seed(7); s0 <- .Random.seed
  cw_e <- capture_warnings(run(bstrap = FALSE, cband = TRUE, seed = 11L))
  expect_identical(.Random.seed, s0); expect_identical(cw_e$value$inference$rng$mode, "R_sample_Rademacher")
  expect_identical(cw_e$warnings, "Used bootstrap procedure to compute simultaneous confidence band")
  # legacy private law (bstrap = NULL): overall stays analytic, bins bootstrap only under cband
  leg <- run(cband = TRUE, seed = 11L)
  expect_true(leg$inference$bootstrap); expect_identical(leg$overall$inference_kind, "analytic")
  expect_null(leg$selectors$warning); expect_false(leg$selectors$applied)
  expect_identical(leg$inference$flag_law, "private_bootstrap_iff_cband")
})

test_that("request preparation binds 0 to Inf, aligns by identity in any order and refuses unsupported inputs", {
  d <- expand.grid(id = 1:6, period = 1:4)
  d$G <- c(0, 0, 3, 3, 4, 4)[d$id]; d$x <- c(0.5, -1, 1, 0, 2, -2)[d$id]; d$w <- c(1, 2, 1, 1.5, 1, 1)[d$id]
  d$cl <- c(1, 1, 2, 2, 3, 3)[d$id]; d$y <- d$id + d$period
  d <- rbind(d, data.frame(id = 7L, period = 3:4, G = 0, x = 1, w = 1, cl = 1, y = c(1, 2)),   # first appears in period 3
             data.frame(id = 8L, period = 2, G = 3, x = 0, w = 1, cl = 2, y = 5))              # singleton
  d$id <- as.integer(d$id)
  perm <- d[c(seq(nrow(d), 1L, by = -3L), setdiff(seq_len(nrow(d)), seq(nrow(d), 1L, by = -3L))), ]
  p <- chained_request_prepare_inputs(perm, NULL, "y", "period", "id", "G", ~x, "w", c("id", "cl"))
  expect_true(all(is.infinite(p$data_internal$G[p$data_internal$id %in% c(1, 2, 7)])))
  expect_false(any(p$data_internal$G == 0))
  expect_identical(names(p$data_internal), c("id", "period", "y", "G"))
  expect_identical(p$frame_ids, unique(perm$id))                      # first appearance in the supplied order
  expect_true(all(c(7L, 8L) %in% p$frame_ids))
  i <- match(1:8, p$frame_ids)
  expect_identical(p$frame_w[i], c(1, 2, 1, 1.5, 1, 1, 1, 1))
  expect_identical(unname(p$frame_Z[i, "x"]), c(0.5, -1, 1, 0, 2, -2, 1, 0))
  expect_identical(as.double(p$cluster_full[i]), c(1, 1, 2, 2, 3, 3, 1, 2)); expect_identical(p$clustervar, "cl")
  expect_identical(p$recoding$n_units_recoded, 3L); expect_identical(p$target_cohorts, c(3, 4))
  expect_identical(p$calendar, c(1, 2, 3, 4)); expect_identical(p$units$n_obs[match(c(7L, 8L), p$units$id)], c(2L, 1L))
  expect_identical(p$frame_g_arg[i], c(Inf, Inf, 3, 3, 4, 4, Inf, 3))
  # explicit permuted frame with a frame-only unit (its own cohort code 0, weight and cluster)
  fr <- data.frame(id = c(9L, 8:1), G = c(0, 3, 0, 4, 4, 3, 3, 0, 0), x = c(0.25, 0, 1, -2, 2, 0, 1, -1, 0.5),
    w = c(3, 1, 1, 1, 1, 1.5, 1, 2, 1), cl = c(4, 2, 1, 3, 3, 2, 2, 1, 1))
  p2 <- chained_request_prepare_inputs(d, fr, "y", "period", "id", "G", ~x, "w", "cl")
  expect_identical(p2$frame_ids, fr$id); expect_identical(p2$n_frame_only, 1L)
  expect_true(is.infinite(p2$frame_g[1L])); expect_identical(p2$frame_w[1L], 3); expect_identical(p2$cluster_full[1L], 4)
  expect_true(p2$units$frame_only[1L]); expect_identical(p2$recoding$n_units_recoded, 4L)
  expect_null(chained_request_prepare_inputs(d, fr[, setdiff(names(fr), "G")], "y", "period", "id", "G", ~x, "w", "cl")$frame_g_arg)
  # D1/D2/D5 keep early-treated units, frame weights and targets. The catalog
  # classifies their absent reference; public composition never drops them.
  d3 <- d; d3$G[d3$id == 3] <- 1
  cw <- capture_warnings(chained_request_prepare_inputs(d3, NULL, "y", "period", "id", "G"))
  expect_identical(cw$warnings, character())
  expect_identical(cw$value$frame_ids, unique(d3$id))
  expect_identical(nrow(cw$value$data_internal), nrow(d3))
  expect_identical(cw$value$target_cohorts, c(1, 3, 4))
  expect_identical(cw$value$dropped_units$n, 0L)
  expect_true(3L %in% cw$value$frame_ids)
  # Anticipation changes reference eligibility, never frame membership.
  pa <- chained_request_prepare_inputs(d3, NULL, "y", "period", "id", "G", anticipation = 2)
  expect_identical(pa$frame_ids, unique(d3$id)); expect_identical(pa$target_cohorts, c(1, 3, 4))
  # explicit refusals (never silent): changing covariates/weights, clusters, frame coverage, cohort codes.
  # D60: a genuinely changing covariate becomes a dated reference recipe of the
  # route (no anchor is needed while weights are constant); changing weights
  # need an explicit anchor.
  dv <- d; dv$x[dv$id == 1 & dv$period == 2] <- 9
  pv <- chained_request_prepare_inputs(dv, NULL, "y", "period", "id", "G", ~x)
  expect_identical(pv$weight_reduction, "reference_recipes")
  expect_identical(pv$reference_route, "baseline")
  expect_identical(pv$cp04$ipw_method, "baseline")
  expect_identical(pv$reference_source$prov$kind[pv$reference_source$prov$variable == "x"], "dated")
  dw <- d; dw$w[dw$id == 2 & dw$period == 3] <- 7
  expect_error(chained_request_prepare_inputs(dw, NULL, "y", "period", "id", "G", NULL, "w"), "first_period", class = "cdid_chained_request_weight_error")
  expect_error(chained_request_prepare_inputs(d, NULL, "y", "period", "id", "G", NULL, NULL, c("cl", "x")), class = "cdid_chained_request_cluster_error")
  expect_error(chained_request_prepare_inputs(d, fr[, setdiff(names(fr), "cl")], "y", "period", "id", "G", ~x, "w", "cl"), class = "cdid_chained_request_cluster_error")
  expect_error(chained_request_prepare_inputs(d, fr[, setdiff(names(fr), "w")], "y", "period", "id", "G", ~x, "w"), class = "cdid_chained_request_weight_error")
  expect_error(chained_request_prepare_inputs(d, fr[fr$id != 2, ], "y", "period", "id", "G", ~x), class = "cdid_chained_request_frame_error")
  f2 <- fr; f2$G[f2$id == 4] <- 0
  expect_error(chained_request_prepare_inputs(d, f2, "y", "period", "id", "G", ~x), class = "cdid_chained_request_cohort_error")
  di <- d; di$G[di$id == 1] <- Inf
  expect_error(chained_request_prepare_inputs(di, NULL, "y", "period", "id", "G"), class = "cdid_chained_request_cohort_error")
  expect_error(chained_request_prepare_inputs(rbind(d, d[1, ]), NULL, "y", "period", "id", "G"), class = "cdid_chained_request_input_error")
  expect_error(chained_request_prepare_inputs(d, NULL, "y", "period", NULL, "G"), class = "cdid_chained_request_input_error")
})

test_that("the executor refuses out-of-scope requests before touching data", {
  ctrl_stub <- structure(list(), class = "cdid_ipw_pipeline_control")
  base <- list(yname = "y", tname = "period", idname = "id", gname = "G", data = quote(stop("data must not be evaluated")),
    est_method = "ipw", control = ctrl_stub, root_witness = function(fv) NULL)
  refuse <- function(cls, ...) expect_error(do.call(chained_att_gt_request, utils::modifyList(base, list(...))), class = cls)
  refuse("cdid_chained_request_mode_error", mode = "ordinary")
  refuse("cdid_chained_panel", panel = FALSE, mode = "chained")
  refuse("cdid_chained_balanced", allow_unbalanced_panel = FALSE)
  refuse("cdid_chained_request_est_method", est_method = "stdipw")    # unsupported estimator strings are never substituted
  # Implemented OR and DR pass the estimator gate and reach the data promise.
  for (em in c("reg", "dr")) {
    expect_error(do.call(chained_att_gt_request, utils::modifyList(base, list(est_method = em))),
      "data must not be evaluated")
  }
  refuse("cdid_chained_request_unsupported_selection", selection = "ipw")
  refuse("cdid_chained_request_unsupported_selection", selection_xformla = ~x)
  refuse("cdid_chained_request_pt_start", pt_assumption = "all")
  refuse("cdid_chained_request_unsupported_weights", fix_weights = "varying")
  # CP08: point estimates only and one declared worker pass; more workers than the
  # explicit budget declares refuse before the data are read.
  for (ok in list(list(compute_inffunc = FALSE, run_numerical_transfer = FALSE), list(pl = TRUE), list(cores = 2)))
    expect_error(do.call(chained_att_gt_request, utils::modifyList(base, ok)), "data must not be evaluated")
  refuse("cdid_chained_request_parallel_error", pl = TRUE, cores = 2)
  refuse("cdid_chained_request_input_error", compute_inffunc = NA)
  # Private interval transfer certifies inference, which point-only fits omit.
  refuse("cdid_chained_request_input_error", compute_inffunc = FALSE)
  # bstrap = FALSE, cband = TRUE is a bound released endpoint, not a refusal: the
  # gate passes and the request proceeds to its data (which this stub refuses).
  expect_error(do.call(chained_att_gt_request, utils::modifyList(base, list(bstrap = FALSE, cband = TRUE))), "data must not be evaluated")
  refuse("cdid_chained_request_input_error", foo = 1)
  expect_error(do.call(chained_att_gt_request, base[setdiff(names(base), "control")]), class = "cdid_chained_request_control_error")
  expect_error(do.call(chained_att_gt_request, base[setdiff(names(base), "root_witness")]), class = "cdid_chained_request_control_error")
})

test_that("aggte options inherit from the MP request unless explicitly overridden; cluster switches are refused", {
  dp <- list(bstrap = TRUE, cband = TRUE, biters = 1000, alp = 0.05, clustervars = "cl")
  r1 <- chained_request_resolve_aggte_options(dp, idname = "id", mp_clustered = TRUE)
  expect_true(r1$bstrap); expect_true(r1$cband); expect_identical(r1$biters, 1000); expect_identical(r1$alp, 0.05)
  expect_identical(r1$clustervars, "cl"); expect_true(all(r1$sources == "inherited"))
  r2 <- chained_request_resolve_aggte_options(dp, cband = FALSE, biters = 50, alp = 0.1, idname = "id", mp_clustered = TRUE)
  expect_false(r2$cband); expect_identical(r2$biters, 50); expect_identical(r2$alp, 0.1)
  expect_identical(unname(r2$sources[c("bstrap", "cband", "biters", "alp")]), c("inherited", "explicit", "explicit", "explicit"))
  # explicit bstrap = FALSE with inherited cband = TRUE resolves to the released aggte law (no refusal)
  r3 <- chained_request_resolve_aggte_options(dp, bstrap = FALSE, idname = "id", mp_clustered = TRUE)
  expect_false(r3$bstrap); expect_true(r3$cband)
  expect_identical(unname(r3$sources[c("bstrap", "cband")]), c("explicit", "inherited"))
  r4 <- chained_request_resolve_aggte_options(list(bstrap = FALSE, cband = TRUE, biters = 10, alp = .05, clustervars = NULL), idname = "id")
  expect_false(r4$bstrap); expect_true(r4$cband); expect_true(all(r4$sources[c("bstrap", "cband")] == "inherited"))
  expect_error(chained_request_resolve_aggte_options(list(bstrap = FALSE, cband = FALSE, biters = 10, alp = .05, clustervars = NULL),
    clustervars = "cl", idname = "id"), class = "cdid_chained_request_cluster_error")
  expect_error(chained_request_resolve_aggte_options(dp, clustervars = "other", idname = "id", mp_clustered = TRUE), class = "cdid_chained_request_cluster_error")
  expect_identical(chained_request_resolve_aggte_options(dp, clustervars = c("id", "cl"), idname = "id", mp_clustered = TRUE)$clustervars, "cl")
})

test_that("no-refit aggregation on the frozen N200 result fields reproduces the fit-time aggregations", {
  z <- readRDS(req_locate_fixture("cp02_result_adapter_N200.rds"))
  rep_ <- z$reporting; fr <- z$frame; G <- fr$G; G[G == 0] <- Inf
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  expect_if_mixed <- function(x, y) {
    schema_ok <- typeof(x) == "double" && typeof(y) == "double" &&
      length(x) > 0L && identical(length(x), length(y)) &&
      identical(attributes(x), attributes(y), num.eq = FALSE)
    expect_true(schema_ok)
    if (!schema_ok) return(invisible(NULL))
    finite_ok <- all(is.finite(x)) && all(is.finite(y))
    expect_true(finite_ok)
    if (!finite_ok) return(invisible(NULL))
    expect_identical(x == 0, y == 0)
    expect_identical(sign(x), sign(y))
    expect_identical(1 / x[x == 0], 1 / y[y == 0])
    expect_true(all(abs(x - y) <= 1e-12 * pmax(1, abs(x), abs(y))))
  }
  for (typ in c("simple", "group", "dynamic", "calendar")) {
    a <- aggregate_chained_ipw_targets(rep_$estimates, rep_$if_full, rep_$target_meta, G, fr$sample_weight, type = typ,
      multiplier_budget = req_budget(), inference_accuracy = req_accuracy(), available = rep_$available,
      bstrap = FALSE, caller_rng = TRUE, calendar = z$calendar, cohort_grid = sort(unique(G)))
    o <- z$aggregations[[typ]]
    expect_identical(a$cell_target_ids, o$cell_target_ids); expect_identical(a$status, o$status)
    expect_identical(a$cell_status$available, o$cell_status$available); expect_identical(a$overall$used_bins, o$overall$used_bins)
    expect_identical(unname(a$estimates), unname(o$estimates)); expect_identical(unname(a$cohort_shares), unname(o$cohort_shares))
    expect_if_mixed(a$if_full, o$if_full)
    expect_if_mixed(a$overall$if_full, o$overall$if_full)
    # SEs are sqrt of cluster-sum cross products; BLAS summation order may differ by platform (one ulp).
    expect_equal(unname(a$se), unname(o$se), tolerance = 1e-14); expect_equal(a$overall$se, o$overall$se, tolerance = 1e-14)
    expect_equal(unname(a$se_fixed_shares_only), unname(o$se_fixed_shares_only), tolerance = 1e-14)
  }
  expect_false(exists(".Random.seed", envir = globalenv()))
})

# ==============================================================================
# B. Installed-host checks (did 2.5.1 + cdid DLL): bounded end-to-end requests
# ==============================================================================

req_fit <- function(bstrap, cband, biters = 200L, ..., inputs = req_dev_inputs(), control = make_request_test_control()) {
  chained_att_gt_request(yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~x, data = inputs$panel,
    control_group = "notyettreated", weightsname = "weight", bstrap = bstrap, cband = cband, biters = biters,
    clustervars = "cluster", est_method = "ipw", base_period = "varying", frame = inputs$frame,
    control = control, root_witness = req_auto_witness, ...)
}

test_that("analytic pointwise request: complete did-compatible MP, retained state, no draw, selectors on aggte", {
  inputs <- req_dev_inputs()
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  mp <- req_fit(bstrap = FALSE, cband = FALSE)
  expect_false(exists(".Random.seed", envir = globalenv()))                 # analytical call invents no draw
  expect_s3_class(mp, c("cdid_MP", "MP"), exact = TRUE)
  tm <- mp$chained$target_meta
  expect_identical(mp$group, tm$g); expect_identical(mp$t, tm$time); expect_true(all(tm$g %in% c(3, 4)))
  expect_identical(length(mp$att), nrow(tm)); expect_identical(length(mp$se), nrow(tm))
  expect_identical(dim(mp$inffunc), c(52L, nrow(tm))); expect_identical(rownames(mp$inffunc), as.character(1:52))
  expect_identical(colnames(mp$inffunc), tm$target_id); expect_identical(mp$n, 52L)
  expect_identical(dim(mp$V_analytical), c(nrow(tm), nrow(tm)))
  av <- mp$chained$available
  expect_equal(mp$se[av], unname(sqrt(diag(mp$V_analytical)[av] / mp$n)))
  expect_equal(mp$c, stats::qnorm(0.975)); expect_identical(mp$alp, 0.05)
  expect_null(mp$W); expect_null(mp$Wpval); expect_null(mp$aggte)
  dp <- mp$DIDparams
  expect_identical(dp$est_method, "ipw"); expect_identical(dp$base_period, "varying"); expect_identical(dp$control_group, "notyettreated")
  expect_false(dp$bstrap); expect_false(dp$cband); expect_identical(dp$biters, 200L); expect_identical(dp$clustervars, "cluster")
  expect_identical(dp$tlist, as.double(1:4)); expect_identical(dp$glist, c(3, 4)); expect_identical(dp$n, 52L)
  expect_identical(dp$cluster_vector, inputs$frame$cluster); expect_null(dp$data); expect_identical(dp$cdid$mode, "chained")
  # public never-treated code bound explicitly; frame-only units and the zero-eligible cluster are retained
  fr <- mp$chained$frame
  expect_identical(fr$id, 1:52); expect_identical(fr$cohort[c(1, 17, 33, 51)], c(3, 4, Inf, Inf))
  expect_identical(sum(fr$frame_only), 4L); expect_identical(sum(!fr$eligible), 4L)
  expect_identical(mp$chained$inference$n_clusters, 13L); expect_identical(mp$chained$inference$n_zero_eligible_clusters, 1L)
  expect_identical(mp$chained$recoding$n_units_recoded, 18L); expect_identical(mp$chained$inference$kind, "analytic")
  # pretrend comparisons survive; the unconditional Wald test is honestly not computed, never a conditional pretest
  pre <- mp$chained$pretest
  expect_identical(pre$status, "not_computed"); expect_true(pre$n_available_pretrend >= 1L)
  expect_true(all(pre$candidate_cells %in% tm$target_id[tm$role == "pre"]))
  expect_match(pre$policy_entry_point, "chained_request_bind_mp"); expect_match(pre$not_conditional_pretest, "cdid::conditional_did_pretest()", fixed = TRUE)
  expect_identical(names(pre$unresolved_policies), c("null_set", "covariance_rule", "acceptance_threshold"))
  # retained aggregation state, no owners/fits
  st <- mp$chained$aggregation_state
  expect_identical(names(st$enclosure), c("theta", "rows")); expect_identical(dim(st$enclosure$rows$lo), c(52L, nrow(tm)))
  expect_identical(st$calendar, as.double(1:4)); expect_identical(st$cohort_grid, c(3, 4, Inf))
  expect_true(all(c("multiplier_budget", "inference_accuracy") %in% names(st)))
  expect_false(any(c("panel_owner", "frame_owner", "gmm_fit", "nuisance_registry", "point_links_all") %in% names(mp$chained)))
  expect_true(isTRUE(mp$chained$numerical_transfer$complete_inference_certified))
  expect_true(mp$chained$numerical_transfer$supported_subvector_certified)
  expect_identical(mp$chained$reduction_summary$n_links_used < mp$chained$reduction_summary$n_links_all, TRUE)
  # methods
  expect_silent(invisible(utils::capture.output(summary(mp))))
  td <- tidy(mp); expect_identical(td$target.id, tm$target_id); expect_true(all(td$available == av))
  expect_identical(glance(mp)$inference, "analytic")
  # aggte requests: inherited analytic flags, selectors with the released conventions, no refit, no draw
  g <- chained_aggte_request(mp, type = "group")
  expect_s3_class(g, c("cdid_AGGTEobj", "AGGTEobj"), exact = TRUE)
  expect_identical(g$egt, c(3, 4)); expect_false(g$DIDparams$bstrap); expect_identical(g$chained$request$resolved$sources[["bstrap"]], "inherited")
  expect_false(g$chained$selectors$applied); expect_false(g$chained$request$refit)
  expect_identical(rownames(g$inf.function$selective.inf.func.g), as.character(1:52))
  g1 <- chained_aggte_request(mp, type = "group", max_e = 0)
  expect_true(g1$chained$selectors$applied)
  expect_identical(g1$chained$cell$intended_ids, c("ATT(g=3,t=3)", "ATT(g=4,t=4)"))      # recoded-grid horizon 0
  dy <- chained_aggte_request(mp, type = "dynamic", balance_e = 1)
  expect_identical(dy$min_e, -Inf); expect_identical(dy$max_e, Inf); expect_identical(dy$balance_e, 1)
  expect_identical(dy$chained$selectors$eligible_cohorts, c(3)); expect_true(all(dy$egt <= 1))
  dw <- chained_aggte_request(mp, type = "dynamic", min_e = 0, max_e = 1)
  expect_identical(dw$egt, c(0, 1)); expect_identical(dw$min_e, 0); expect_identical(dw$max_e, 1)
  expect_error(chained_aggte_request(mp, type = "dynamic", min_e = 10), class = "cdid_ipw_pipeline_selector_error")
  cw <- capture_warnings(chained_aggte_request(mp, type = "calendar", max_e = 1))
  expect_match(cw$warnings, "ignored for type = \"calendar\"", fixed = TRUE)
  expect_true(cw$value$chained$selectors$ignored); expect_identical(cw$value$egt, c(3, 4))
  s <- chained_aggte_request(mp, type = "simple", na.rm = TRUE)
  expect_true(s$chained$na_rm); expect_null(s$egt); expect_true(is.finite(s$overall.att))
  expect_false(exists(".Random.seed", envir = globalenv()))
  # explicit overrides with the familiar meaning
  g2 <- chained_aggte_request(mp, type = "group", alp = 0.1)
  expect_identical(g2$DIDparams$alp, 0.1); expect_equal(g2$crit.val.egt, stats::qnorm(0.95))
  expect_identical(g2$chained$request$resolved$sources[["alp"]], "explicit")
  # explicit cband = TRUE on the analytic MP (inherited bstrap = FALSE): released
  # aggte law, warning, analytic SEs, multiplier band from the caller's stream
  set.seed(31); s0 <- .Random.seed
  cwb <- capture_warnings(chained_aggte_request(mp, type = "group", cband = TRUE)); gb <- cwb$value
  expect_identical(cwb$warnings, "Used bootstrap procedure to compute simultaneous confidence band")
  expect_false(identical(.Random.seed, s0))
  expect_false(gb$DIDparams$bstrap); expect_true(gb$DIDparams$cband)
  expect_identical(gb$se.egt, g$se.egt); expect_identical(gb$overall.se, g$overall.se)
  expect_identical(gb$chained$inference$kind, "analytic"); expect_true(gb$chained$inference$cband)
  expect_true(is.finite(gb$crit.val.egt) && gb$crit.val.egt != stats::qnorm(0.975))
  expect_identical(gb$chained$band_family, "simultaneous_multiplier_sup_t_over_certified_ids")
  expect_true(gb$chained$inference$released_false_true)
  txt <- paste(utils::capture.output(summary(gb)), collapse = "\n")
  expect_match(txt, "Simult.", fixed = TRUE); expect_match(txt, "analytic standard errors", fixed = TRUE)
  expect_match(txt, "standardized by their own IQR scale", fixed = TRUE)
  expect_error(chained_aggte_request(mp, type = "group", clustervars = "x"), class = "cdid_chained_request_cluster_error")
  expect_silent(invisible(utils::capture.output(summary(g1))))
  expect_identical(nrow(tidy(dy)), length(dy$egt))
})

test_that("bstrap = FALSE, cband = TRUE request: released att_gt endpoint at the fit, released aggte endpoint on inheritance", {
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  mp_ft <- req_fit(bstrap = FALSE, cband = TRUE)
  expect_false(exists(".Random.seed", envir = globalenv()))                 # no draw at the att_gt endpoint
  expect_s3_class(mp_ft, c("cdid_MP", "MP"), exact = TRUE)
  expect_false(mp_ft$DIDparams$bstrap); expect_true(mp_ft$DIDparams$cband)   # requested flags retained
  expect_identical(mp_ft$DIDparams$cdid$requested_cband, TRUE); expect_identical(mp_ft$DIDparams$cdid$effective_cband, FALSE)
  expect_identical(mp_ft$chained$inference$kind, "analytic"); expect_false(mp_ft$chained$inference$cband)
  expect_identical(mp_ft$chained$band_family, "pointwise_normal_over_certified_ids")
  expect_identical(mp_ft$chained$inference$flag_law, "released_bstrap_false_cband_true_cells_analytic_pointwise_no_draws")
  expect_equal(mp_ft$c, stats::qnorm(0.975))
  av <- mp_ft$chained$available
  expect_identical(mp_ft$se[av], unname(mp_ft$chained$analytical_se[av]))
  expect_true(isTRUE(mp_ft$chained$numerical_transfer$complete_inference_certified))
  mp_ff <- req_fit(bstrap = FALSE, cband = FALSE)
  expect_identical(mp_ft$att, mp_ff$att); expect_identical(mp_ft$se, mp_ff$se); expect_identical(mp_ft$inffunc, mp_ff$inffunc)
  expect_identical(mp_ft$V_analytical, mp_ff$V_analytical); expect_identical(mp_ft$chained$cband_lower, mp_ff$chained$cband_lower)
  txt <- paste(utils::capture.output(summary(mp_ft)), collapse = "\n")
  expect_match(txt, "Pointwise", fixed = TRUE); expect_no_match(txt, "Simult.", fixed = TRUE)
  expect_match(txt, "released did 2.5.1 att_gt", fixed = TRUE)
  expect_identical(glance(mp_ft)$inference, "analytic")
  # inherited FALSE/TRUE on aggte: group/dynamic/calendar warn and bootstrap the band from the caller stream
  set.seed(11); s0 <- .Random.seed
  cw <- capture_warnings(chained_aggte_request(mp_ft, type = "dynamic")); ag <- cw$value; s1 <- .Random.seed
  expect_identical(cw$warnings, "Used bootstrap procedure to compute simultaneous confidence band")
  expect_false(identical(s0, s1))
  expect_false(ag$DIDparams$bstrap); expect_true(ag$DIDparams$cband)
  expect_identical(ag$chained$request$resolved$sources[["cband"]], "inherited")
  expect_identical(ag$chained$inference$kind, "analytic"); expect_true(ag$chained$inference$cband)
  expect_identical(ag$chained$overall$inference_kind, "analytic")
  ag_ff <- chained_aggte_request(mp_ff, type = "dynamic")
  expect_identical(ag$se.egt, ag_ff$se.egt); expect_identical(ag$overall.se, ag_ff$overall.se); expect_identical(ag$att.egt, ag_ff$att.egt)
  expect_true(is.finite(ag$crit.val.egt) && ag$crit.val.egt != stats::qnorm(0.975))
  set.seed(11); cw2 <- capture_warnings(chained_aggte_request(mp_ft, type = "dynamic"))
  expect_identical(cw2$value$crit.val.egt, ag$crit.val.egt)                 # reproducible from the caller seed
  set.seed(11); ag_tt <- chained_aggte_request(mp_ft, type = "dynamic", bstrap = TRUE)
  expect_identical(ag_tt$crit.val.egt, ag$crit.val.egt)                     # same bin draws, bootstrap SEs instead
  expect_false(isTRUE(all.equal(ag_tt$se.egt, ag$se.egt)))
  txt <- paste(utils::capture.output(summary(ag)), collapse = "\n")
  expect_match(txt, "Simult.", fixed = TRUE); expect_match(txt, "analytic standard errors", fixed = TRUE)
  # simple inherits FALSE/TRUE: no bin-band family, no warning, no draw, analytic
  s2 <- .Random.seed
  cws <- capture_warnings(chained_aggte_request(mp_ft, type = "simple")); sm <- cws$value
  expect_identical(cws$warnings, character()); expect_identical(.Random.seed, s2)
  expect_true(sm$DIDparams$cband); expect_false(sm$chained$inference$cband); expect_identical(sm$chained$inference$kind, "analytic")
  expect_identical(sm$overall.se, chained_aggte_request(mp_ff, type = "simple")$overall.se)
  # explicit overrides still work on this MP
  ag_a <- chained_aggte_request(mp_ft, type = "group", cband = FALSE)
  expect_identical(.Random.seed, s2); expect_false(ag_a$chained$inference$cband); expect_equal(ag_a$crit.val.egt, stats::qnorm(0.975))
})

test_that("bootstrap requests use the caller's RNG honestly for both pointwise and simultaneous inference", {
  set.seed(2026); s0 <- .Random.seed
  mp_tf <- req_fit(bstrap = TRUE, cband = FALSE, biters = 120L); s1 <- .Random.seed
  expect_false(identical(s0, s1))                                          # advanced, never reset
  set.seed(2026); mp_tt <- req_fit(bstrap = TRUE, cband = TRUE, biters = 120L); s2 <- .Random.seed
  expect_identical(s1, s2)                                                 # same stream consumption
  set.seed(2026); mp_tf2 <- req_fit(bstrap = TRUE, cband = FALSE, biters = 120L)
  expect_identical(mp_tf2$se, mp_tf$se)                                    # reproducible from the caller seed
  expect_identical(mp_tt$se, mp_tf$se)                                     # same multiplier IQR SEs
  expect_equal(mp_tf$c, stats::qnorm(0.975)); expect_true(mp_tt$c > stats::qnorm(0.975) || mp_tt$c != mp_tf$c)
  expect_true(mp_tf$DIDparams$bstrap); expect_false(mp_tf$DIDparams$cband); expect_true(mp_tt$DIDparams$cband)
  expect_identical(mp_tf$chained$inference$kind, "multiplier_bootstrap_iqr_se"); expect_identical(mp_tf$chained$inference$biters, 120L)
  expect_identical(mp_tf$chained$band_family, "pointwise_normal_over_certified_ids")
  expect_identical(mp_tt$chained$band_family, "simultaneous_multiplier_sup_t_over_certified_ids")
  expect_identical(mp_tf$chained$inference$rng$mode, "caller_rng_stream_Rademacher")
  expect_true(isTRUE(mp_tt$chained$numerical_transfer$complete_inference_certified))
  # analytical SEs retained separately; requested fields stay separate from effective fields
  expect_false(isTRUE(all.equal(mp_tf$se, unname(mp_tf$chained$analytical_se))))
  expect_identical(mp_tf$DIDparams$cdid$requested_bstrap, TRUE); expect_identical(mp_tf$DIDparams$cdid$actual_draws, 120L)
  txt <- paste(utils::capture.output(summary(mp_tt)), collapse = "\n"); expect_match(txt, "Simult.", fixed = TRUE)
  txt <- paste(utils::capture.output(summary(mp_tf)), collapse = "\n"); expect_match(txt, "Pointwise", fixed = TRUE)
  # aggte inherits bstrap = TRUE: draws come from the caller stream (advanced); explicit bstrap = FALSE draws nothing
  set.seed(9); a0 <- .Random.seed; ag <- chained_aggte_request(mp_tf, type = "dynamic"); a1 <- .Random.seed
  expect_false(identical(a0, a1)); expect_true(ag$DIDparams$bstrap); expect_identical(ag$chained$inference$kind, "multiplier_bootstrap_iqr_se")
  expect_identical(ag$chained$overall$inference_kind, "multiplier_bootstrap_iqr_se"); expect_equal(ag$crit.val.egt, stats::qnorm(0.975))
  set.seed(9); ag2 <- chained_aggte_request(mp_tf, type = "dynamic"); expect_identical(ag2$se.egt, ag$se.egt)
  b0 <- .Random.seed; ag_a <- chained_aggte_request(mp_tf, type = "dynamic", bstrap = FALSE); expect_identical(.Random.seed, b0)
  expect_false(ag_a$DIDparams$bstrap); expect_identical(ag_a$chained$inference$kind, "analytic")
  expect_identical(ag_a$chained$request$resolved$sources[["bstrap"]], "explicit")
  ag_c <- chained_aggte_request(mp_tt, type = "group")
  expect_true(ag_c$DIDparams$cband); expect_identical(ag_c$chained$band_family, "simultaneous_multiplier_sup_t_over_certified_ids")
  expect_true(ag_c$crit.val.egt != stats::qnorm(0.975))
})

test_that("public wrapper activates default DR and explicit IPW", {
  expect_error(cdid::att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = stop("data must not be evaluated")),
    "data must not be evaluated")
  expect_error(cdid::att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = stop("public IPW reaches data"), est_method = "ipw", mode = "chained"),
    "public IPW reaches data")
  expect_identical(formals(cdid::aggte), formals(did::aggte))
  expect_false("chained_att_gt_request" %in% getNamespaceExports("cdid"))
  expect_true(exists("chained_att_gt_request", envir = asNamespace("cdid"), inherits = FALSE))
})
