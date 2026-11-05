# Source-owned accounting witnesses. These use the unchanged fit/replay paths;
# numerical reference fixtures and solver tolerances are not rewritten.
fma_ns <- function(n) getFromNamespace(n, "cdid")
fma_case <- function(long_labels = FALSE) {
  N <- 24L; id <- seq_len(N); g <- rep(c(3, Inf), each = N / 2L)
  x <- rep(rep(-2:3, 2L), 2L); w <- 1 + (id %% 3L) / 10
  Y <- vapply(1:3, function(d) d + d * x / 3 + sin(id + d) / 10 +
    as.double(g == 3) * d, numeric(N))
  ctrl <- fma_ns("ipw_owner_control")(3L, 1e9, 1e9, 1e14, 1e9)
  data <- data.frame(id = rep(id, each = 3L), time = rep(1:3, N),
    y = as.vector(t(Y)), g = rep(g, each = 3L))
  data <- data[!((data$id %% 4L) == (data$time - 1L)), ]
  panel <- fma_ns("construct_ipw_panel_owner")(data, "id", "time", "y", "g", ctrl)
  label <- if (long_labels) paste(rep("long_covariate_label", 80L), collapse = "") else "x"
  X <- cbind(1, x); colnames(X) <- c("(Intercept)", label)
  frame <- fma_ns("construct_ipw_frame_owner")(panel, id, X, w, ctrl)
  reader <- fma_ns("construct_cp04_frame_reader")(frame, id,
    matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept")), g, w, ctrl)
  views <- lapply(1:2, function(j) fma_ns("construct_ipw_link_view")(frame,
    j, 3, 3, Inf, control = ctrl, link_id = paste0("field-memory-", j, if (long_labels) label else ""),
    functional_key = paste0("cmpcat:field-memory-", j)))
  list(reader = reader, views = views, control = ctrl, X = X,
    recipe = fma_ns("or_pooled_outcome_recipe")(colnames(X), TRUE, label))
}
fma_fit <- function(a, route, method = "dr", outcome = "separate", model = "pair_logit",
                    pooling = "separate", design = "X_frame", prefix = 1e6, control = a$control) {
  args <- list(reader = a$reader, link_views = a$views, control = control,
    rank_tol = 1e-10, outcome_recipe = a$recipe,
    retained_external_entries = max(0, prefix - a$reader$accounting$borrowed_owner_entries),
    outcome_pooling = outcome)
  gc <- list(epsilon = 1e-12, maxit = 50L, trace = FALSE)
  if (route %in% c("imputation", "augmented")) {
    args <- c(args, list(est_method = method, glm_control = gc,
      treatment_columns = "(Intercept)", imputation_pooling = pooling))
    if (design == "X_frame") args$X_frame <- a$X else if (design == "design_owner")
      args$design_owner <- fma_ns("construct_imputation_design_owner")(
        get("frame_owner", a$reader$frame_ref), a$reader, a$reader$state$id, a$X, control)
    if (route == "augmented") args <- c(args, list(selection_model = model, selection_pooling = pooling))
    return(do.call(fma_ns(if (route == "imputation") "fit_imputation_backend" else "fit_augmented_backend"), args))
  }
  args <- c(args, list(glm_control = gc, selection_model = model, selection_pooling = pooling))
  if (method == "dr") args$treatment_columns <- "(Intercept)"
  do.call(fma_ns(if (method == "dr") "fit_dr_corollary1" else "fit_or_corollary1"), args)
}
fma_registry <- function(a, b, prefix = 1e6) list(cp04 = b, cp04_reader = a$reader,
  cp04_memory_accounting = fma_ns("cp04_registry_memory_accounting")(prefix, a$reader, b))
fma_res <- function(a, b) list(nuisance_registry = fma_registry(a, b),
  resource_accounting = list(initial_envelope_entries = 0, simultaneous_entries = 0),
  controls = list(or_basis = NULL), selection = NULL,
  frame_owner = get("frame_owner", a$reader$frame_ref),
  reporting = list(inference = list(n_total = a$reader$state$N)),
  auxiliary_reporting = list(descriptor_status = NULL))

# Finite branch coverage, not a Cartesian campaign: pure IPW/OR/DR; separate
# and pooled outcomes; both observation models and all selection poolings.
fma_specs <- list(
  list("imputation", "ipw", "separate", "pair_logit", "separate"),
  list("imputation", "reg", "pooled", "pair_logit", "within_cohort"),
  list("imputation", "dr", "separate", "pair_logit", "across_cohorts"),
  list("augmented", "dr", "separate", "pair_logit", "within_cohort"),
  list("augmented", "dr", "pooled", "independent_period_logit", "across_cohorts"),
  list("selection", "reg", "separate", "pair_logit", "separate"),
  list("selection", "dr", "pooled", "pair_logit", "within_cohort"),
  list("selection", "dr", "separate", "independent_period_logit", "across_cohorts"))

testthat::test_that("static registry removes only its exact embedded external reservation", {
  a <- fma_case(); prefix <- 1e6
  for (spec in fma_specs) {
    b <- do.call(fma_fit, c(list(a = a), stats::setNames(spec, c("route", "method", "outcome", "model", "pooling"))))
    reg <- fma_registry(a, b, prefix); rec <- reg$cp04_memory_accounting
    testthat::expect_identical(rec$kind, "cp04_registry_field_memory_prefix")
    X <- max(0, prefix - a$reader$accounting$borrowed_owner_entries)
    retained <- prefix + b$accounting$owned_entries + a$reader$accounting$owned_entries + 5000
    adjusted <- fma_ns("cp04_registry_retained_entries")(reg, retained)
    testthat::expect_equal(retained - adjusted, X)
    testthat::expect_equal(b$accounting$borrowed_entries, a$reader$accounting$owned_entries +
      a$reader$accounting$borrowed_owner_entries + X +
      sum(unlist(b$accounting[c("design_owned_entries", "design_borrowed_entries")])))
    pre <- fma_ns("cp04_preflight_links")(b, 1:2, adjusted)
    standalone <- fma_ns("cp04_preflight_links")(b, 1:2, retained)
    testthat::expect_equal(pre$peak_concurrent_entries, standalone$peak_concurrent_entries)
    plain <- reg; plain$cp04_memory_accounting <- NULL
    testthat::expect_identical(fma_ns("cp04_registry_retained_entries")(plain, retained), retained)
    testthat::expect_error(fma_ns("cp04_registry_retained_entries")(reg, prefix - 1), class = "cdid_cp04_integrity_error")
    wrong <- reg; wrong$cp04_reader <- fma_case()$reader
    testthat::expect_error(fma_ns("cp04_registry_retained_entries")(wrong, retained), class = "cdid_cp04_integrity_error")
    bad <- reg; bad$cp04_memory_accounting$embedded_external_entries <- X - 1
    testthat::expect_error(fma_ns("cp04_registry_retained_entries")(bad, retained), class = "cdid_cp04_integrity_error")
    bad <- reg; bad$cp04$accounting$borrowed_entries <- bad$cp04$accounting$borrowed_entries - 1
    testthat::expect_error(fma_ns("cp04_registry_retained_entries")(bad, retained))
  }
})

testthat::test_that("field envelopes cover complete diagnostic outputs and partial public results", {
  a <- fma_case(TRUE)
  for (spec in fma_specs) {
    b <- do.call(fma_fit, c(list(a = a), stats::setNames(spec, c("route", "method", "outcome", "model", "pooling"))))
    res <- fma_res(a, b); public <- list(partial = a$X, alias = a$X)
    plan <- fma_ns("cp04_public_diagnostic_plan")(b, res$selection, res$controls$or_basis)
    out <- fma_ns("cp04_public_diagnostics")(res, list(owner_control = a$control), public)
    testthat::expect_lte(unname(fma_ns("cp04_payload")(out)["entries"]), plan$output)
    testthat::expect_equal(plan$peak, plan$output + plan$scratch)
    testthat::expect_equal(plan$scratch, 4 * plan$output)
    # OR/DR retain complete links/supports; reconstruction's catalog has only
    # selected summaries and does not return its link adjoints.
    if (spec[[1L]] == "selection") {
      testthat::expect_identical(out$links, b$links)
      testthat::expect_length(out$supports, length(b$supports))
    } else testthat::expect_null(out$links)
    design <- if (is.null(b$design_ref)) 0 else
      get("imputation_design_owner", b$design_ref)$accounting$owned_entries
    retained <- unname(fma_ns("cp04_payload")(res)["entries"]) + design +
      unname(fma_ns("cp04_payload")(public)["entries"])
    cap <- retained + plan$peak
    exact <- a$control; exact$max_matrix_entries <- cap
    testthat::expect_equal(fma_ns("cp04_public_diagnostic_admit")(res,
      list(owner_control = exact), public, "field boundary"), plan$output)
    exact$max_matrix_entries <- cap - 1
    testthat::expect_error(fma_ns("cp04_public_diagnostic_admit")(res,
      list(owner_control = exact), public, "field boundary"), class = "cdid_ipw_owner_bound_error")
    bad <- b; bad$links[[1L]]$link_id <- "mutated"
    testthat::expect_error(fma_ns("cp04_public_diagnostic_plan")(bad))
  }
})

testthat::test_that("external design-owner prefix retains its separate source reservation", {
  a <- fma_case()
  b <- fma_fit(a, "imputation", design = "design_owner")
  reg <- fma_registry(a, b)
  d <- get("imputation_design_owner", b$design_ref)
  testthat::expect_equal(reg$cp04_memory_accounting$design_owned_entries, 0)
  testthat::expect_equal(reg$cp04_memory_accounting$design_borrowed_entries, d$accounting$owned_entries)
  wrong <- reg; wrong$cp04 <- fma_fit(a, "imputation", design = "design_owner")
  testthat::expect_error(fma_ns("cp04_registry_retained_entries")(wrong, 2e6), class = "cdid_cp04_integrity_error")
})

testthat::test_that("complete auxiliary tables include unsupported rows and unbounded source messages", {
  a <- fma_case(); Q <- 257L
  descriptors <- lapply(seq_len(Q), function(j) list(g = 3, s = 1, t = 2,
    K = c(Inf, j + 4)))
  index <- rep(NA_integer_, Q); index[1L] <- 1L
  message <- paste(rep("long independent support reason", 100L), collapse = " ")
  errors <- lapply(seq_len(Q), function(j) if (j == 1L) NULL else
    structure(list(message = paste0(message, j), call = NULL),
      class = c("cdid_ipw_support_error", "error", "condition")))
  # Independent reproduction of the existing reporting table, without a fit
  # for any unsupported auxiliary row or any truncation of its support reason.
  reasons <- vapply(seq_len(Q), function(j) if (j == 1L) NA_character_ else
    sprintf("%s [g=3, s=1, t=2, K={Inf,%d}]", errors[[j]]$message, j + 4L), character(1))
  table <- data.frame(index = seq_len(Q), g = rep(3, Q), s = rep(1, Q), t = rep(2, Q),
    K = vapply(seq_len(Q), function(j) sprintf("{Inf,%d}", j + 4L), character(1)),
    status = c("supported", rep("unsupported_pair", Q - 1L)), reason = reasons, stringsAsFactors = FALSE)
  # All four diagnostics return this field, including Q much larger than the
  # fitted backend's link count. Its arbitrary message bytes need source funds.
  for (route in c("imputation", "augmented", "or", "dr")) {
    b <- fma_fit(a, if (route %in% c("or", "dr")) "selection" else route,
      method = if (route == "or") "reg" else "dr")
    res <- fma_res(a, b)
    res$nuisance_registry$cp04_aux <- list(descriptors = descriptors)
    res$nuisance_registry$cp04_aux_index <- index
    res$nuisance_registry$cp04_aux_errors <- errors
    res$auxiliary_reporting$descriptor_status <- table
    aux <- fma_ns("cp04_public_auxiliary_diagnostic_bound")(res$nuisance_registry)
    testthat::expect_lte(unname(fma_ns("cp04_payload")(table)["entries"]), aux)
    omitted <- fma_ns("cp04_public_diagnostic_plan")(b)
    plan <- fma_ns("cp04_public_diagnostic_plan")(b, auxiliary_entries = aux)
    testthat::expect_equal(plan$output - omitted$output, aux)
    out <- fma_ns("cp04_public_diagnostics")(res, list(owner_control = a$control))
    testthat::expect_identical(out$auxiliary_descriptor_status, table)
    testthat::expect_gt(unname(fma_ns("cp04_payload")(out)["entries"]), omitted$output)
    testthat::expect_lte(unname(fma_ns("cp04_payload")(out)["entries"]), plan$output)
    design <- if (is.null(b$design_ref)) 0 else
      get("imputation_design_owner", b$design_ref)$accounting$owned_entries
    exact <- a$control
    exact$max_matrix_entries <- unname(fma_ns("cp04_payload")(res)["entries"]) + design + plan$peak
    testthat::expect_equal(fma_ns("cp04_public_diagnostic_admit")(res,
      list(owner_control = exact), NULL, "full auxiliary boundary"), plan$output)
    exact$max_matrix_entries <- exact$max_matrix_entries - 1
    testthat::expect_error(fma_ns("cp04_public_diagnostic_admit")(res,
      list(owner_control = exact), NULL, "full auxiliary boundary"), class = "cdid_ipw_owner_bound_error")
    bad <- res; bad$auxiliary_reporting$descriptor_status$reason[2L] <- paste(rep(message, 300L), collapse = " ")
    testthat::expect_error(fma_ns("cp04_public_diagnostic_admit")(bad,
      list(owner_control = a$control), NULL, "mutated auxiliary table"), class = "cdid_cp04_bound_error")
  }
})
