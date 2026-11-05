# Static legacy CP04 lineage only. Numerical fitting and replay are unchanged.
lmp_ns <- function(name) getFromNamespace(name, "cdid")
lmp_case <- function() {
  patterns <- rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(1, 0, 1, 0),
    c(0, 1, 0, 1), c(0, 0, 1, 1), c(1, 0, 0, 1),
    c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 0, 0))
  id <- seq_len(54L); x <- rep(rep(-1:1, each = 9L), 2L)
  h <- rep(1:9, 6L); g <- rep(c(3, Inf), each = 27L)
  w <- 1 + (id %% 3L) / 10; obs <- patterns[h, , drop = FALSE]
  Y <- vapply(1:4, function(d) d * (x + h / 20) + d^2 * x^2 / 10 +
    as.double(g == 3) * d * (1 + x / 4), numeric(length(id)))
  Y[obs == 0] <- NA_real_
  data <- data.frame(id = rep(id, each = 4L), time = rep(1:4, length(id)),
    y = as.vector(t(Y)), g = rep(g, each = 4L))
  zero <- which(rowSums(obs) == 0L)
  data <- data[!data$id %in% zero[-1L], ]
  control <- lmp_ns("ipw_owner_control")(7L, 1e8, 1e8, 1e13, 1e8)
  panel <- lmp_ns("construct_ipw_panel_owner")(data, "id", "time", "y", "g", control)
  Z <- cbind(`(Intercept)` = 1, x = x)
  frame <- lmp_ns("construct_ipw_frame_owner")(panel, id, Z, w, control,
    frame_g = g, allow_frame_only = TRUE)
  reader <- lmp_ns("construct_cp04_frame_reader")(frame, id, Z, g, w, control)
  views <- lapply(list(c(1, 2), c(2, 4)), function(pair)
    lmp_ns("construct_ipw_link_view")(frame, pair[1L], pair[2L], 3, Inf,
      control = control, link_id = paste0("legacy-memory-", pair[2L]),
      functional_key = paste0("cmpcat:legacy-memory-", pair[2L])))
  list(reader = reader, views = views, control = control, frame = frame, panel = panel,
    zero = zero, single = which(rowSums(obs) == 1L))
}
lmp_fit <- function(a, prefix = 0, control = a$control, pooling = "separate")
  lmp_ns("fit_cp04_corollary1")(a$reader, a$views, control,
    list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    retained_external_entries = max(0, prefix - a$reader$accounting$borrowed_owner_entries),
    selection_pooling = pooling)
lmp_registry <- function(a, b, prefix) list(cp04 = b, cp04_reader = a$reader,
  cp04_memory_accounting = lmp_ns("cp04_registry_memory_accounting")(prefix, a$reader, b))
lmp_collect <- function(b) {
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  phi <- matrix(NA_real_, reader$state$N, length(b$links))
  info <- lmp_ns("replay_cp04_if")(b, 7L, 2L, function(chunk)
    phi[chunk$rows, chunk$links] <<- chunk$phi, scale = "frame")
  list(phi = phi, info = info,
    point = vapply(seq_along(b$links), function(j) lmp_ns("chained_full_link_estimate")(b, j), numeric(1)),
    covariance = crossprod(sweep(phi, 2L, colMeans(phi))) / reader$state$N^2)
}

testthat::test_that("legacy CP04 subtracts only its exact embedded registry prefix", {
  a <- lmp_case(); E0 <- 1e6; b <- lmp_fit(a, E0)
  X <- max(0, E0 - a$reader$accounting$borrowed_owner_entries)
  R <- a$reader$accounting$owned_entries
  reg <- lmp_registry(a, b, E0); retained <- E0 + b$accounting$owned_entries + R + 5000
  testthat::expect_false(isTRUE(b$uses_selection_models))
  testthat::expect_identical(reg$cp04_memory_accounting$kind, "cp04_registry_field_memory_prefix")
  testthat::expect_equal(b$accounting$retained_external_entries, X)
  testthat::expect_equal(b$accounting$borrowed_entries,
    R + a$reader$accounting$borrowed_owner_entries + X)
  testthat::expect_equal(b$accounting$owned_entries, unname(lmp_ns("cp04_payload")(b)["entries"]))
  without_scalar <- b; without_scalar$accounting$retained_external_entries <- NULL
  testthat::expect_lt(unname(lmp_ns("cp04_payload")(without_scalar)["entries"]), b$accounting$owned_entries)
  adjusted <- lmp_ns("cp04_registry_retained_entries")(reg, retained)
  testthat::expect_equal(adjusted, retained - X)
  pre <- lmp_ns("cp04_preflight_links")(b, 1:2, adjusted)
  cap <- ceiling(max(b$accounting$retained_bound + b$accounting$fit_scratch_bound,
    adjusted + pre$peak_concurrent_entries)) + 64
  testthat::expect_lt(cap, retained + pre$peak_concurrent_entries)
  lower <- a$control; lower$max_matrix_entries <- cap
  fit <- lmp_fit(a, E0, lower); registry <- lmp_registry(a, fit, E0)
  corrected <- lmp_ns("cp04_registry_retained_entries")(registry, retained)
  admitted <- lmp_ns("cp04_preflight_links")(fit, 1:2, corrected)
  testthat::expect_lte(corrected + admitted$peak_concurrent_entries, cap)
  testthat::expect_identical(admitted$chunk_rows, pre$chunk_rows)
  testthat::expect_identical(admitted$work_ops_bound, pre$work_ops_bound)
  testthat::expect_error(lmp_ns("cp04_preflight_links")(fit, 1:2, retained),
    class = "cdid_ipw_owner_bound_error")
  testthat::expect_error(lmp_ns("cp04_preflight_links")(fit, 1:2, corrected + cap),
    class = "cdid_ipw_owner_bound_error")
  plain <- registry; plain$cp04_memory_accounting <- NULL
  testthat::expect_identical(lmp_ns("cp04_registry_retained_entries")(plain, retained), retained)
  testthat::expect_equal(lmp_ns("cp04_registry_retained_entries")(registry, E0), E0 - X)
  testthat::expect_error(lmp_ns("cp04_registry_retained_entries")(registry, E0 - 1),
    class = "cdid_cp04_integrity_error")
  zero <- lmp_fit(a, 0); zero_reg <- lmp_registry(a, zero, 0)
  testthat::expect_null(zero_reg$cp04_memory_accounting)
  testthat::expect_equal(zero$accounting$retained_external_entries, 0)
  testthat::expect_identical(lmp_ns("cp04_registry_retained_entries")(zero_reg, retained), retained)
})

testthat::test_that("legacy prefix lineage rejects foreign, tampered and resealed components", {
  a <- lmp_case(); E0 <- 1e6; b <- lmp_fit(a, E0)
  reg <- lmp_registry(a, b, E0); retained <- 2 * E0
  wrong <- reg; wrong$cp04 <- lmp_fit(a, E0)
  testthat::expect_error(lmp_ns("cp04_registry_retained_entries")(wrong, retained), class = "cdid_cp04_integrity_error")
  wrong <- reg; wrong$cp04_reader <- lmp_case()$reader
  testthat::expect_error(lmp_ns("cp04_registry_retained_entries")(wrong, retained), class = "cdid_cp04_integrity_error")
  for (component in c("record", "backend", "reader")) {
    bad <- reg
    if (component == "record") bad$cp04_memory_accounting$embedded_external_entries <- 0 else
      if (component == "backend") bad$cp04$accounting$borrowed_entries <- 0 else
        bad$cp04_reader$accounting$owned_entries <- 0
    testthat::expect_error(lmp_ns("cp04_registry_retained_entries")(bad, retained), class = "cdid_cp04_integrity_error")
  }
  # A valid replacement seal does not establish the component equation.
  bad <- b; bad$accounting$borrowed_entries <- bad$accounting$borrowed_entries - 1
  bad <- lmp_ns("ipw_owner_seal")(bad, "cp04_backend")
  testthat::expect_error(lmp_registry(a, bad, E0), class = "cdid_cp04_integrity_error")
  forged <- reg; forged$cp04 <- bad
  forged$cp04_memory_accounting$backend_ref <- bad$ref
  forged$cp04_memory_accounting <- lmp_ns("ipw_owner_seal")(
    forged$cp04_memory_accounting, "cp04_registry_field_memory_prefix")
  testthat::expect_error(lmp_ns("cp04_registry_retained_entries")(forged, retained), class = "cdid_cp04_integrity_error")
  bad <- reg; bad$cp04_memory_accounting$pipeline_prefix_entries <- E0 + 1
  bad$cp04_memory_accounting <- lmp_ns("ipw_owner_seal")(
    bad$cp04_memory_accounting, "cp04_registry_field_memory_prefix")
  testthat::expect_error(lmp_ns("cp04_registry_retained_entries")(bad, retained), class = "cdid_cp04_integrity_error")
})

testthat::test_that("legacy lineage leaves numerical F rows and excluded routes unchanged", {
  a <- lmp_case(); E0 <- 1e6; b <- lmp_fit(a, E0); standalone <- lmp_fit(a, 0)
  testthat::expect_gt(a$frame$n_frame_only_units, 0)
  testthat::expect_gt(a$panel$n_singleton_units, 0)
  testthat::expect_gt(a$panel$n_unobserved_units, 0)
  got <- lmp_collect(b); before <- lmp_collect(standalone)
  testthat::expect_identical(got, before)
  testthat::expect_equal(nrow(got$phi), a$reader$state$N)
  testthat::expect_true(all(is.finite(got$phi[c(a$zero, a$single), , drop = FALSE])))
  testthat::expect_gt(max(abs(got$phi[a$zero, , drop = FALSE])), 0)
  testthat::expect_gt(max(abs(got$phi[a$single, , drop = FALSE])), 0)
  diagnostic <- function(fit) lmp_ns("cp04_public_diagnostics")(
    list(nuisance_registry = list(cp04 = fit), resource_accounting = list(),
      frame_owner = a$frame, reporting = list(inference = list(n_total = a$reader$state$N)),
      auxiliary_reporting = list()), list(owner_control = a$control))
  current <- diagnostic(b); old <- diagnostic(standalone)
  testthat::expect_identical(current$records, old$records)
  testthat::expect_identical(current$status, old$status)
  testthat::expect_identical(current$links, old$links)
  general <- lmp_fit(a, E0, pooling = "within_cohort")
  testthat::expect_true(general$uses_selection_models)
  testthat::expect_null(lmp_ns("cp04_registry_memory_accounting")(E0, a$reader, general))
  for (field in c("reference_ref", "layout")) {
    excluded <- b; excluded[[field]] <- new.env(parent = emptyenv())
    testthat::expect_null(lmp_ns("cp04_registry_memory_accounting")(E0, a$reader, excluded))
  }
  for (cls in c("cdid_cp04_custom_backend", "cdid_reference_ipw_backend", "cdid_ipw_raw_backend")) {
    excluded <- structure(list(), class = cls)
    testthat::expect_null(lmp_ns("cp04_registry_memory_accounting")(E0, a$reader, excluded))
  }
})
