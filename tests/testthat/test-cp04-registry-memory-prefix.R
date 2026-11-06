# Source-owned prefix accounting witnesses; no estimator/solver reference is changed.
mpx_ns <- function(name) getFromNamespace(name, "cdid")
mpx_case <- function() {
  N <- 12L; id <- seq_len(N); g <- rep(c(3, Inf), each = 6L)
  x <- rep(c(-2, -1, 0, 1, 2, 3), 2L); w <- 1 + (id %% 3L) / 10
  Y <- vapply(1:3, function(d) d + d * x / 3 + sin(id + d) / 10 +
    as.double(g == 3) * d, numeric(N))
  ctrl <- mpx_ns("ipw_owner_control")(3L, 1e8, 1e8, 1e13, 1e8)
  data <- data.frame(id = rep(id, each = 3L), time = rep(1:3, N),
    y = as.vector(t(Y)), g = rep(g, each = 3L))
  panel <- mpx_ns("construct_ipw_panel_owner")(data, "id", "time", "y", "g", ctrl)
  X <- cbind(`(Intercept)` = 1, x = x)
  frame <- mpx_ns("construct_ipw_frame_owner")(panel, id, X, w, ctrl)
  reader <- mpx_ns("construct_cp04_frame_reader")(frame, id,
    matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept")), g, w, ctrl)
  views <- lapply(1:2, function(j) mpx_ns("construct_ipw_link_view")(frame,
    j, 3, 3, Inf, control = ctrl, link_id = paste0("memory-", j), functional_key = paste0("cmpcat:memory-", j)))
  list(reader = reader, views = views, control = ctrl,
    recipe = mpx_ns("or_pooled_outcome_recipe")(colnames(X), TRUE, "x"))
}
mpx_fit <- function(a, method, prefix, control = a$control) {
  external <- max(0, prefix - a$reader$accounting$borrowed_owner_entries)
  args <- list(reader = a$reader, link_views = a$views, control = control,
    rank_tol = 1e-10, outcome_recipe = a$recipe, retained_external_entries = external,
    outcome_pooling = "separate")
  if (method == "dr") args <- c(args, list(glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    treatment_columns = "(Intercept)"))
  do.call(mpx_ns(if (method == "reg") "fit_or_baseline" else "fit_dr_baseline"), args)
}
mpx_registry <- function(a, fit, prefix) list(cp04 = fit, cp04_reader = a$reader,
  cp04_memory_accounting = mpx_ns("cp04_registry_memory_accounting")(prefix, a$reader, fit))

testthat::test_that("same registry subtracts exactly its embedded external prefix once", {
  a <- mpx_case(); E0 <- 1e6
  for (method in c("reg", "dr")) {
    high <- mpx_fit(a, method, E0); registry <- mpx_registry(a, high, E0)
    X <- max(0, E0 - a$reader$accounting$borrowed_owner_entries)
    R <- a$reader$accounting$owned_entries; B <- high$accounting$owned_entries
    retained <- E0 + B + R + 5000
    corrected <- mpx_ns("cp04_registry_retained_entries")(registry, retained)
    testthat::expect_equal(retained - corrected, X)
    testthat::expect_equal(high$accounting$borrowed_entries,
      R + a$reader$accounting$borrowed_owner_entries + X)
    pre <- mpx_ns("cp04_preflight_links")(high, 1:2, corrected)
    required <- corrected + pre$peak_concurrent_entries
    # A real newly fitted sealed backend under this unchanged finite cap.
    cap <- ceiling(max(high$accounting$retained_bound + high$accounting$fit_scratch_bound,
      required)) + 64
    testthat::expect_lt(cap, required + X)
    lower <- a$control; lower$max_matrix_entries <- cap
    fit <- mpx_fit(a, method, E0, lower); reg <- mpx_registry(a, fit, E0)
    actual_retained <- E0 + fit$accounting$owned_entries + R + 5000
    adjusted <- mpx_ns("cp04_registry_retained_entries")(reg, actual_retained)
    admitted <- mpx_ns("cp04_preflight_links")(fit, 1:2, adjusted)
    testthat::expect_lte(adjusted + admitted$peak_concurrent_entries, cap)
    testthat::expect_error(mpx_ns("cp04_preflight_links")(fit, 1:2, actual_retained),
      class = "cdid_ipw_owner_bound_error")
    testthat::expect_error(mpx_ns("cp04_preflight_links")(fit, 1:2, adjusted + cap),
      class = "cdid_ipw_owner_bound_error")
    # The standalone API and a registry without source lineage keep the full
    # caller reservation; equal numeric amounts do not authorize subtraction.
    plain <- list(cp04 = fit, cp04_reader = a$reader)
    testthat::expect_equal(mpx_ns("cp04_registry_retained_entries")(plain, actual_retained), actual_retained)
    testthat::expect_error(mpx_ns("cp04_preflight_links")(plain$cp04, 1:2, actual_retained),
      class = "cdid_ipw_owner_bound_error")
    testthat::expect_equal(fit$control$max_matrix_entries, cap)
    fn <- mpx_ns(if (method == "reg") "evaluate_or_baseline_if_chunk" else "evaluate_dr_if_chunk")
    testthat::expect_equal(fn(fit, 3:1, 2:1)$phi, fn(high, 3:1, 2:1)$phi, tolerance = 1e-12)
    wrong <- reg; wrong$cp04 <- high
    testthat::expect_error(mpx_ns("cp04_registry_retained_entries")(wrong, actual_retained),
      class = "cdid_cp04_integrity_error")
  }
})
