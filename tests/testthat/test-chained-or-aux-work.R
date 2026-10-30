# R1 accounting-unit regression with real installed owners and sealed backends.
# The installed local accounting expression is evaluated in an ordinary test
# scope containing an actual registry. No namespace, owner, reader or backend
# is synthesized or mocked; the installed function body is never changed.
aux_work_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  f <- get(name, ns, inherits = FALSE)
  stopifnot(is.function(f), identical(environment(f), ns)); f
}
aux_work_registry <- function(method) {
  raw <- read.csv(testthat::test_path("fixtures", "cp06_or_corollary1", "raw.csv"))
  raw$cohort <- as.numeric(replace(raw$cohort, raw$cohort == "NT", "Inf"))
  N <- nrow(raw); Y <- as.matrix(raw[paste0("Y", 1:4)])
  long <- data.frame(id = rep(raw$id, each = 4L), time = rep(1:4, N),
    y = as.vector(t(Y)), g = rep(raw$cohort, each = 4L))
  zero <- which(rowSums(!is.na(Y)) == 0L)
  long <- long[!long$id %in% raw$id[zero[-1L]], ]
  ctrl <- aux_work_ns("ipw_owner_control")(23L, 134217728, 134217728, 1e11, 134217728)
  panel <- aux_work_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  X <- cbind(intercept = 1, x = raw$x); X[raw$A == 0, ] <- NA_real_
  frame <- aux_work_ns("construct_ipw_frame_owner")(panel, raw$id, X,
    weights = raw$w, control = ctrl, frame_g = raw$cohort,
    allow_frame_only = TRUE, design_on_A_only = TRUE)
  reader <- aux_work_ns("construct_cp04_frame_reader")(frame, raw$id,
    cbind(one = 1, x = raw$x), raw$cohort, raw$w, ctrl, U = rep(TRUE, N),
    rho = rep(1, N), universe = "F", justification = "accepted deterministic finite fixture")
  ledger <- data.frame(g = c(3, 4))
  registry <- aux_work_ns("build_cp04_public_registry")(NULL, frame, reader, ledger,
    list(owner_control = ctrl, glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE)),
    "notyettreated", 0L, or_method = method, automatic_chunk_rows = TRUE)
  list(registry = registry, ledger = ledger, N = N, control = ctrl)
}
aux_work_assignment <- function() {
  found <- list()
  visit <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    if (identical(e[[1L]], as.name("<-")) && identical(e[[2L]], as.name("aux_work")))
      found[[length(found) + 1L]] <<- e
    if (is.symbol(e[[1L]]) && as.character(e[[1L]]) %in% c("{", "if", "for", "while", "repeat"))
      for (child in as.list(e)[-1L]) if (is.call(child)) visit(child)
    invisible(NULL)
  }
  visit(body(aux_work_ns("run_chained_ipw_pipeline")))
  stopifnot(length(found) == 1L); found[[1L]]
}
aux_work_account <- function(fx, bp, retained) {
  r <- fx$registry
  scope <- list2env(list(nuis_res = r, cp04 = TRUE, active_backend = r$cp04,
    cp04_auto_plan = inherits(r$cp04, c("cdid_or_corollary1", "cdid_dr_corollary1", "cdid_dr_baseline")),
    use_raw_ipw = FALSE, aux_sel_bp = r$cp04_aux$selections[[bp]],
    aux_ids = r$cp04_aux_index[!is.na(r$cp04_aux_index)], bp_mode = bp,
    retained_with_diagnostics = retained, n_total = fx$N, ledger = fx$ledger),
    parent = asNamespace("cdid"))
  value <- eval(aux_work_assignment(), scope)
  list(work = value, ids = if (exists("ids", scope, inherits = FALSE)) scope$ids else NULL)
}
aux_work_reporting_ids <- function(fx, bp) {
  # Independently follow the reporting path: supported descriptor -> point ->
  # backend cp04_index. Keep selected order and remove duplicates only at replay.
  r <- fx$registry; requested <- r$cp04_aux$selections[[bp]]
  active <- vapply(requested, function(d) !is.na(d$index) && d$index > 0L &&
    !is.na(r$cp04_aux_index[d$index]), logical(1))
  descriptor_ids <- vapply(requested[active], `[[`, integer(1), "index")
  points <- lapply(descriptor_ids, function(j) r$cp04_points[[r$cp04_aux_index[j]]])
  unique(vapply(points, `[[`, integer(1), "cp04_index"))
}

testthat::test_that("C1 OR auxiliary accounting uses the exact active reporting plan", {
  fx <- aux_work_registry("corollary1"); r <- fx$registry; b <- r$cp04
  S <- r$cp04_aux_index[!is.na(r$cp04_aux_index)]
  own <- b$accounting$owned_entries + b$accounting$borrowed_entries
  envelope <- aux_work_ns("or_corollary1_if_envelope")(b, 12L, S)
  reserved <- fx$control$max_matrix_entries - own - envelope$output - envelope$scratch
  for (bp in c("varying", "universal")) {
    T <- aux_work_reporting_ids(fx, bp)
    testthat::expect_lt(length(T), length(S))
    retained <- reserved - 8 * fx$N * (length(r$cp04_aux$selections[[bp]]) + nrow(fx$ledger))
    testthat::expect_gte(retained, 0)
    replay <- aux_work_ns("cp04_preflight_links")(b, T, reserved)
    superset <- aux_work_ns("cp04_preflight_links")(b, S, reserved)
    testthat::expect_equal(superset$chunk_rows, 12, tolerance = 0)
    testthat::expect_gt(replay$chunk_rows, superset$chunk_rows)
    actual <- aux_work_account(fx, bp, retained)
    testthat::expect_identical(actual$ids, T)
    testthat::expect_equal(actual$work, replay$work_ops_bound, tolerance = 0)
    actual_plan <- aux_work_ns("cp04_preflight_links")(b, actual$ids, reserved)
    testthat::expect_equal(actual_plan$chunk_rows, replay$chunk_rows, tolerance = 0)
    testthat::expect_equal(actual_plan$work_ops_bound, replay$work_ops_bound, tolerance = 0)
    testthat::expect_false(identical(superset$work_ops_bound, replay$work_ops_bound))
  }
  testthat::expect_identical(aux_work_reporting_ids(fx, "universal"), c(1L, 4L, 3L))
})

testthat::test_that("baseline OR and CP04 IPW retain their original unreserved accounting preflight", {
  for (method in list("baseline", NULL)) {
    fx <- aux_work_registry(method); r <- fx$registry; b <- r$cp04; bp <- "universal"
    S <- r$cp04_aux_index[!is.na(r$cp04_aux_index)]; T <- aux_work_reporting_ids(fx, bp)
    testthat::expect_lt(length(T), length(S))
    own <- b$accounting$owned_entries + b$accounting$borrowed_entries
    eS <- aux_work_ns("chained_full_if_budget")(b, 23L, S)
    eT <- aux_work_ns("chained_full_if_budget")(b, 23L, T)
    testthat::expect_gt(eS$output + eS$scratch, eT$output + eT$scratch)
    reserved <- fx$control$max_matrix_entries - own - eS$output - eS$scratch + 1
    retained <- reserved - 8 * fx$N * (length(r$cp04_aux$selections[[bp]]) + nrow(fx$ledger))
    testthat::expect_gte(retained, 0)
    testthat::expect_error(aux_work_ns("cp04_preflight_links")(b, S, reserved),
      class = "cdid_ipw_owner_bound_error")
    reporting <- aux_work_ns("cp04_preflight_links")(b, T, reserved)
    testthat::expect_equal(reporting$chunk_rows, 23, tolerance = 0)
    original <- aux_work_ns("cp04_preflight_links")(b, S)
    actual <- aux_work_account(fx, bp, retained)
    testthat::expect_equal(actual$work, original$work_ops_bound, tolerance = 0)
  }
})
