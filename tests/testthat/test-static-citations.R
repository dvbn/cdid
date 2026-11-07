# ==============================================================================
# Narrow Reference and Citation Tests for Static (`dname`) vs Cohort (`gname`)
# ==============================================================================

local_find_repo_file <- function(rel_path) {
  candidates <- c(
    rel_path,
    file.path("..", "..", rel_path),
    file.path("..", rel_path)
  )
  for (p in candidates) {
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  if ("cdid" %in% loadedNamespaces() || nzchar(system.file(package = "cdid"))) {
    pkg_file <- system.file(sub("^inst/", "", rel_path), package = "cdid")
    if (nzchar(pkg_file) && file.exists(pkg_file)) {
      return(pkg_file)
    }
  }
  stop("Could not locate file: ", rel_path)
}

local_load_results_helpers <- function() {
  if ("cdid" %in% loadedNamespaces()) {
    ns <- asNamespace("cdid")
    if (exists("new_cdid_static_att_gt", envir = ns, inherits = FALSE) &&
        exists("chained_result_references", envir = ns, inherits = FALSE)) {
      return(ns)
    }
  }
  env <- new.env(parent = globalenv())
  source(local_find_repo_file("R/chained_did_results.R"), local = env)
  source(local_find_repo_file("R/static_results.R"), local = env)
  # Register S3 methods in globalenv if running standalone outside namespace
  for (fn in c(
    "tidy.cdid_static_att_gt",
    "glance.cdid_static_att_gt",
    "summary.cdid_static_att_gt",
    "print.summary.cdid_static_att_gt",
    "print.cdid_static_att_gt"
  )) {
    assign(fn, get(fn, envir = env, inherits = FALSE), envir = globalenv())
  }
  env
}

local_make_minimal_static_result <- function(helpers) {
  new_fn <- get("new_cdid_static_att_gt", envir = helpers, inherits = FALSE)
  tid <- c("t2_entry", "t2_exit")
  req <- list(
    call = quote(cdid::att_gt(yname = "y", tname = "time", idname = "id", dname = "d", data = df)),
    yname = "y",
    tname = "time",
    idname = "id",
    dname = "d",
    xformla = ~ x,
    est_method_requested = "dr",
    est_method = "dr",
    route = "augmented",
    selection_requested = "conditional_mtar",
    selection = "conditional_mtar",
    selection_est_method_requested = "dr",
    selection_est_method = "dr",
    selection_model = "pair_logit",
    selection_pooling = "within_group",
    imputation_pooling = "within_group",
    imputation_scale = "per_period",
    outcome_pooling_mode = "separate",
    gmm_weighting = list(canonical = "optimal", requested = "optimal"),
    clustervars = NULL,
    faster_mode = TRUE,
    print_details = FALSE,
    pl = FALSE,
    cores = 1L,
    compute_inffunc = TRUE
  )
  inf <- list(
    se_reported = "analytical",
    band_reported = "pointwise_normal",
    clustered = FALSE,
    clustername = NULL,
    n_clusters = 4L,
    n_units = 4L,
    qualifications = "Asymptotic validity requires compatible link moments and regular interior roots."
  )
  diag_list <- list(
    n_units = 4L,
    n_dates = 3L,
    n_repeat_outcome_units = 3L,
    n_single_outcome_units = 1L,
    n_zero_outcome_units = 0L,
    n_declared_targets = 2L,
    n_available_targets = 2L,
    n_unavailable_targets = 0L,
    n_total_links = 6L
  )
  agg <- list(
    available = TRUE,
    status = "ok",
    reason = NA_character_,
    estimate = 1.5,
    se = 0.1,
    ci_lower = 1.304,
    ci_upper = 1.696,
    total_switchers = 2L,
    scope = "static switcher-frequency weighted average across the declared transition set A"
  )
  new_fn(
    target_id = tid,
    t = c(2, 2),
    t_prev = c(1, 1),
    d = c(0L, 1L),
    direction = c("entry", "exit"),
    sigma = c(1, -1),
    available = c(TRUE, TRUE),
    support_status = c("ok", "ok"),
    support_reason = c(NA_character_, NA_character_),
    n_risk = c(2L, 2L),
    n_switchers = c(1L, 1L),
    n_stayers = c(1L, 1L),
    att = stats::setNames(c(1.5, 1.5), tid),
    theta_hat = stats::setNames(c(1.5, -1.5), tid),
    se = stats::setNames(c(0.12, 0.14), tid),
    analytical_se = stats::setNames(c(0.12, 0.14), tid),
    bootstrap_se = stats::setNames(c(NA_real_, NA_real_), tid),
    V_analytical = diag(c(0.12^2, 0.14^2)),
    c = 1.959964,
    pointwise_crit_val = 1.959964,
    cband_crit_val = NA_real_,
    ci_lower = stats::setNames(c(1.5 - 1.959964 * 0.12, 1.5 - 1.959964 * 0.14), tid),
    ci_upper = stats::setNames(c(1.5 + 1.959964 * 0.12, 1.5 + 1.959964 * 0.14), tid),
    cband_lower = stats::setNames(c(NA_real_, NA_real_), tid),
    cband_upper = stats::setNames(c(NA_real_, NA_real_), tid),
    alp = 0.05,
    bstrap = FALSE,
    cband = FALSE,
    biters = 1000L,
    n = 4L,
    unit_id = 1:4,
    dates = 1:3,
    dates_num = c(1, 2, 3),
    inffunc = matrix(0, nrow = 4L, ncol = 2L, dimnames = list(as.character(1:4), tid)),
    aggregate = agg,
    graphs = list(rule = "all_positive_observed_count_pairs_per_role_no_pruning"),
    request = req,
    inference = inf,
    diagnostics = diag_list,
    fit = NULL
  )
}

testthat::test_that("inst/CITATION lists BBD (2025), CS (2021), and de Chaisemartin et al. (2026) with route guidance", {
  cit_path <- local_find_repo_file("inst/CITATION")
  cit <- utils::readCitationFile(cit_path)
  testthat::expect_length(cit, 3L)

  titles <- vapply(cit, function(e) as.character(e$title), character(1L))
  testthat::expect_identical(
    titles,
    c(
      "The Chained Difference-in-Differences",
      "Difference-in-Differences with Multiple Time Periods",
      "Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period"
    )
  )

  # Entry 1: Bellego, Benatia, Dortet-Bernadet (2025)
  aut1 <- format(cit[[1L]]$author)
  testthat::expect_true(any(grepl("Bellego", aut1, fixed = TRUE)))
  testthat::expect_true(any(grepl("Benatia", aut1, fixed = TRUE)))
  testthat::expect_true(any(grepl("Dortet-Bernadet", aut1, fixed = TRUE)))
  testthat::expect_identical(as.character(cit[[1L]]$doi), "10.1016/j.jeconom.2024.105783")

  # Entry 2: Callaway and Sant'Anna (2021) for cohort gname and ordinary mode
  aut2 <- format(cit[[2L]]$author)
  testthat::expect_true(any(grepl("Callaway", aut2, fixed = TRUE)))
  testthat::expect_true(any(grepl("Sant'Anna", aut2, fixed = TRUE)))
  testthat::expect_identical(as.character(cit[[2L]]$doi), "10.1016/j.jeconom.2020.12.001")
  testthat::expect_match(cit[[2L]]$header, "gname", fixed = TRUE)
  testthat::expect_match(cit[[2L]]$header, "ordinary", fixed = TRUE)

  # Entry 3: de Chaisemartin, D'Haultfoeuille, Pasquier, Sow, Vazquez-Bare (2026) for static dname route
  aut3 <- format(cit[[3L]]$author)
  for (surname in c("de Chaisemartin", "D'Haultfoeuille", "Pasquier", "Sow", "Vazquez-Bare")) {
    testthat::expect_true(any(grepl(surname, aut3, fixed = TRUE)), info = surname)
  }
  testthat::expect_match(cit[[3L]]$note, "arXiv:2201.06898v7", fixed = TRUE)
  testthat::expect_match(cit[[3L]]$note, "2026-07-22", fixed = TRUE)
  testthat::expect_match(cit[[3L]]$note, "static binary specialization", fixed = TRUE)
  testthat::expect_match(cit[[3L]]$header, "dname", fixed = TRUE)
  testthat::expect_false(grepl("Intertemporal", cit[[3L]]$title, fixed = TRUE))
})

testthat::test_that("static result metadata and summary cite de Chaisemartin et al. alongside BBD (2025) and never cite CS (2021)", {
  helpers <- local_load_results_helpers()
  obj <- local_make_minimal_static_result(helpers)
  s_obj <- summary(obj)

  for (ref_slot in list(obj$references, obj$request$references, obj$diagnostics$references, obj$DIDparams$references, s_obj$references)) {
    testthat::expect_true(is.list(ref_slot))
    testthat::expect_match(ref_slot$cdid, "Bellego", fixed = TRUE)
    testthat::expect_match(ref_slot$cdid, "Benatia", fixed = TRUE)
    testthat::expect_match(ref_slot$cdid, "Dortet-Bernadet", fixed = TRUE)
    testthat::expect_match(ref_slot$cdid, "The Chained Difference-in-Differences", fixed = TRUE)
    testthat::expect_match(ref_slot$cdid, "10.1016/j.jeconom.2024.105783", fixed = TRUE)

    testthat::expect_match(ref_slot$static, "de Chaisemartin", fixed = TRUE)
    testthat::expect_match(ref_slot$static, "D'Haultfoeuille", fixed = TRUE)
    testthat::expect_match(ref_slot$static, "Pasquier", fixed = TRUE)
    testthat::expect_match(ref_slot$static, "Sow", fixed = TRUE)
    testthat::expect_match(ref_slot$static, "Vazquez-Bare", fixed = TRUE)
    testthat::expect_match(
      ref_slot$static,
      "Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period",
      fixed = TRUE
    )
    testthat::expect_match(ref_slot$static, "arXiv:2201.06898v7", fixed = TRUE)
    testthat::expect_match(ref_slot$static, "2026-07-22", fixed = TRUE)

    joined_refs <- paste(unlist(ref_slot), collapse = " ")
    testthat::expect_false(grepl("Callaway", joined_refs, fixed = TRUE))
    testthat::expect_false(grepl("Sant'Anna", joined_refs, fixed = TRUE))
    testthat::expect_false(grepl("2020.12.001", joined_refs, fixed = TRUE))
    testthat::expect_false(grepl("Intertemporal", joined_refs, fixed = TRUE))
  }

  out_sum <- paste(utils::capture.output(print(s_obj)), collapse = "\n")
  out_prt <- paste(utils::capture.output(print(obj)), collapse = "\n")
  for (out_txt in c(out_sum, out_prt)) {
    testthat::expect_match(out_txt, "Bellego C., Benatia D. and V. Dortet-Bernadet (2025)", fixed = TRUE)
    testthat::expect_match(out_txt, "The Chained Difference-in-Differences", fixed = TRUE)
    testthat::expect_match(
      out_txt,
      "de Chaisemartin C., D'Haultfoeuille X., Pasquier F., Sow D. and G. Vazquez-Bare (2026)",
      fixed = TRUE
    )
    testthat::expect_match(
      out_txt,
      "Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period",
      fixed = TRUE
    )
    testthat::expect_match(out_txt, "arXiv:2201.06898v7", fixed = TRUE)
    testthat::expect_false(grepl("Callaway", out_txt, fixed = TRUE))
    testthat::expect_false(grepl("Sant'Anna", out_txt, fixed = TRUE))
    testthat::expect_false(grepl("Difference-in-Differences with Multiple Time Periods", out_txt, fixed = TRUE))
    testthat::expect_false(grepl("Intertemporal", out_txt, fixed = TRUE))
  }
})

testthat::test_that("existing cohort (gname) reference helper retains BBD (2025) and Callaway and Sant'Anna (2021) unchanged", {
  helpers <- local_load_results_helpers()
  ch_refs_fn <- get("chained_result_references", envir = helpers, inherits = FALSE)

  for (m in c("ipw", "reg", "dr")) {
    refs <- ch_refs_fn(est_method = m)
    testthat::expect_identical(names(refs), c("estimator", "cdid", "did"))
    testthat::expect_match(refs$cdid, "Bellego C., Benatia D. and V. Dortet-Bernadet (2025)", fixed = TRUE)
    testthat::expect_match(refs$cdid, "The Chained Difference-in-Differences", fixed = TRUE)
    testthat::expect_match(refs$cdid, "10.1016/j.jeconom.2024.105783", fixed = TRUE)
    testthat::expect_match(refs$did, "Callaway B. and P. H. C. Sant'Anna (2021)", fixed = TRUE)
    testthat::expect_match(refs$did, "Difference-in-Differences with Multiple Time Periods", fixed = TRUE)
    testthat::expect_match(refs$did, "10.1016/j.jeconom.2020.12.001", fixed = TRUE)
    testthat::expect_false(grepl("de Chaisemartin", paste(unlist(refs), collapse = " "), fixed = TRUE))
  }
})
