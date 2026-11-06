# Admission and payload checks use real sealed owners and unchanged fits.
res_ns <- function(name) getFromNamespace(name, "cdid")
res_gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)

res_case <- function(universe = "F", intercept = TRUE, slopes = TRUE, dated = FALSE, route_kind = "corollary1", extra_reference = FALSE, or_only = FALSE) {
  grid <- expand.grid(R3 = 0:1, R2 = 0:1, R1 = 0:1, z = c(-1, 1),
    u = c(-1, 1), k = 1:3)
  N <- nrow(grid); id <- seq_len(N); calendar <- c(2, 5, 11)
  obs <- as.matrix(grid[c("R1", "R2", "R3")]) != 0
  A <- rowSums(obs) >= 2; U <- switch(universe, F = rep(TRUE, N),
    B = rowSums(obs) >= 1, A = A)
  z <- grid$z; k <- grid$k; g <- c(11, 23, Inf)[k]
  x <- grid$u + .25 * (grid$R1 - grid$R3) + .1 * k * z
  prob <- vapply(1:3, function(d) stats::plogis(-.6 + .3 * d + .2 * k + .35 * z), numeric(N))
  w <- apply(prob^obs * (1 - prob)^(!obs), 1L, prod) * (1 + .15 * grid$u)
  pu <- if (universe == "F") rep(1, N) else if (universe == "B")
    1 - apply(1 - prob, 1L, prod) else
      prob[, 1] * prob[, 2] + prob[, 1] * prob[, 3] + prob[, 2] * prob[, 3] -
        2 * apply(prob, 1L, prod)
  rho <- pu
  for (kk in 1:3) rho[k == kk] <- pu[k == kk] /
    (sum(w[k == kk] * U[k == kk]) / sum(w[k == kk]))
  Y <- vapply(1:3, function(d) 1 + d + k * d / 3 +
    (.3 * k * d + .09 * d^2) * x + .13 * d * x^2 + .07 * d * z * grid$u, numeric(N))
  Y[!obs] <- NA_real_
  X <- matrix(numeric(), N, 0)
  if (intercept) X <- cbind(`(Intercept)` = rep(1, N))
  if (slopes) X <- cbind(X, x = x)
  Q <- cbind(`(Intercept)` = 1, z = z)
  profile <- res_ns("chained_public_profile_defaults")()
  entries <- floor(profile$numeric_workspace_bytes / 8)
  ctrl <- res_ns("ipw_owner_control")(11L, entries, entries, res_ns("chained_public_work_arithmetic_ceiling")(), entries)
  long <- data.frame(id = rep(id, each = 3), time = rep(calendar, N),
    y = as.vector(t(Y)), g = rep(g, each = 3))
  panel <- res_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  treatment_v <- grid$u + .3 * z
  V <- cbind(`(Intercept)` = 1, treatment_v = treatment_v)
  Xa <- if (intercept) cbind(X, treatment_v = treatment_v) else cbind(`(Intercept)` = 1, X, treatment_v = treatment_v)
  if (or_only) Xa <- X
  Xa[!A, ] <- NA_real_; Qa <- Q; Qa[!U, ] <- NA_real_
  frame <- res_ns("construct_ipw_frame_owner")(panel, id, Xa, w, ctrl,
    frame_g = g, design_on_A_only = TRUE)
  reader <- res_ns("construct_cp04_frame_reader")(frame, id, Qa, g, w, ctrl,
    U = U, rho = rho, universe = universe,
    justification = "known finite-law observation-subset transport ratios")
  defs <- list(list(s = 2, t = 11, g = 11, K = c(23, Inf)),
    list(s = 5, t = 11, g = 11, K = Inf), list(s = 2, t = 11, g = 23, K = Inf))
  views <- lapply(seq_along(defs), function(i) {
    a <- defs[[i]]
    res_ns("construct_ipw_link_view")(frame, a$s, a$t, a$g, a$K, control = ctrl,
      link_id = paste0("selection-pooling-", i), functional_key = paste0("cmpcat:selection-pooling-", i))
  })
  recipe <- res_ns("or_pooled_outcome_recipe")(colnames(X), intercept, setdiff(colnames(X), "(Intercept)"))
  ref <- NULL; ref_X <- ref_Q <- ref_W <- ref_V <- NULL
  if (dated) {
    stopifnot(universe == "F", intercept, slopes)
    di <- rep(1:3, N)
    long$x <- rep(x, each = 3) + di * rep(z * grid$u / 7, each = 3)
    long$z <- rep(z, each = 3) * (1 + di / 20)
    long$treatment_v <- rep(treatment_v, each = 3) + di * rep(z / 30, each = 3)
    long$w <- rep(w, each = 3) * (1 + di * rep((grid$u + 1) / 10, each = 3))
    key <- match(long$time, calendar) * (N + 1) + long$id; wo <- order(key)
    route <- res_ns("chained_reference_route")(route_kind, "F"); route$or_basis <- TRUE
    formulas <- list(x = list(formula = ~x + treatment_v))
    if (!route_kind %in% c("baseline", "imputation")) formulas$s <- list(formula = ~z)
    reference_vars <- c("x", "z", "treatment_v")
    if (extra_reference) {
      long$unused_a <- rep(x^2 + .2 * grid$R2, each = 3) + di / 13
      long$unused_b <- rep(z * grid$u + .3 * grid$R1, each = 3) + di * rep(grid$R3 / 17, each = 3)
      formulas$x$formula <- ~x + treatment_v + unused_a + unused_b
      reference_vars <- c(reference_vars, "unused_a", "unused_b")
    }
    src <- res_ns("chained_reference_source")(formulas,
      data.frame(variable = reference_vars, kind = "dated", source = "data"), long,
      data.frame(id = id), FALSE, long$id, long$time, "time", id, calendar, 0,
      list(kind = "varying", name = "w", key = key[wo], value = long$w[wo], V = w,
        example = "finite-law dated reference weights"), "base_period", route)
    ref <- res_ns("build_chained_reference_reader")(src, reader, views, length(views), ctrl)
    recipe <- res_ns("or_pooled_outcome_recipe")(c("(Intercept)", "x"), TRUE, "x")
    ref_V <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, treatment_v = treatment_v + d * z / 30))
    ref_X <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, x = x + d * z * grid$u / 7))
    ref_Q <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, z = z * (1 + d / 20)))
    ref_W <- lapply(ref$recipes$weight_index, function(d) w * (1 + d * (grid$u + 1) / 10))
  }
  list(reader = reader, views = views, control = ctrl, recipe = recipe, reference = ref,
    X = X, V = V, X_full = if (intercept) cbind(X, treatment_v = treatment_v) else cbind(`(Intercept)` = 1, X, treatment_v = treatment_v), Q = Q, Y = Y, obs = obs, g = g, w = w, U = U, A = A,
    rho = rho, calendar = calendar, ref_X = ref_X, ref_Q = ref_Q, ref_W = ref_W, ref_V = ref_V)
}

res_fit <- function(a, kind, control = a$control) {
  base <- list(reader = a$reader, link_views = a$views, control = control,
    rank_tol = 1e-10, outcome_recipe = a$recipe, outcome_pooling = "separate")
  if (kind != "or") base$glm_control <- res_gc
  if (kind %in% c("dr", "selection_dr", "imputation", "augmented"))
    base$treatment_columns <- c("(Intercept)", "treatment_v")
  if (kind %in% c("imputation", "augmented")) {
    base$X_frame <- if (is.null(a$reference)) a$X_full else NULL
    base$reference <- a$reference
  }
  if (kind == "selection_or") base$selection_pooling <- "within_cohort"
  if (kind == "selection_dr") base$selection_pooling <- "within_cohort"
  fn <- switch(kind, or = "fit_or_baseline", selection_or = "fit_or_corollary1",
    dr = "fit_dr_baseline", selection_dr = "fit_dr_corollary1",
    imputation = "fit_imputation_backend", augmented = "fit_augmented_backend")
  do.call(res_ns(fn), base)
}
res_eval <- function(fit, kind, rows, links, details) {
  fn <- switch(kind, or = "evaluate_or_baseline_if_chunk", selection_or = "evaluate_or_corollary1_if_chunk",
    dr = "evaluate_dr_if_chunk", selection_dr = "evaluate_dr_if_chunk",
    imputation = "evaluate_imputation_if_chunk", augmented = "evaluate_augmented_if_chunk")
  res_ns(fn)(fit, rows, links, details = details)
}
res_budget <- function(fit, kind, B, links, details) {
  fn <- switch(kind, or = "or_baseline_if_budget", selection_or = "or_corollary1_if_budget",
    dr = "dr_if_budget", selection_dr = "dr_if_budget", imputation = "imputation_if_budget", augmented = "augmented_if_budget")
  res_ns(fn)(fit, B, links, details)
}

testthat::test_that("separate local envelopes admit normal-profile batches and dominate actual payload", {
  for (kind in c("or", "selection_or", "dr", "selection_dr", "imputation", "augmented")) {
    a <- res_case(or_only = kind %in% c("or", "selection_or")); fit <- res_fit(a, kind)
    testthat::expect_lte(fit$accounting$owned_entries, fit$accounting$retained_bound)
    ordinary <- res_budget(fit, kind, 3L, c(3L, 1L), FALSE)
    detailed <- res_budget(fit, kind, 3L, c(3L, 1L), TRUE)
    testthat::expect_lt(ordinary$output, detailed$output)
    for (details in c(FALSE, TRUE)) {
      pre <- if (details) detailed else ordinary
      out <- res_eval(fit, kind, c(10L, 3L, 1L), c(3L, 1L), details)
      testthat::expect_lte(unname(res_ns("cp04_payload")(out)["entries"]), pre$output)
      testthat::expect_lt(fit$accounting$owned_entries + fit$accounting$borrowed_entries + pre$scratch,
        fit$control$max_matrix_entries)
    }
    # Explicit limits still refuse before an IF buffer is evaluated. Build a
    # new properly sealed fit with the caller's lower output allowance.
    lower <- a$control; lower$max_output_entries <- ceiling(ordinary$output) - 1
    # Fit admission can itself exceed this limit, which is also a valid early refusal.
    testthat::expect_error({
      bounded <- res_fit(a, kind, lower)
      res_budget(bounded, kind, 3L, c(3L, 1L), FALSE)
    }, class = "cdid_ipw_owner_bound_error")
    broken <- fit; broken$law <- paste0(broken$law, ":changed")
    testthat::expect_error(res_eval(broken, kind, 1L, 1L, FALSE))
  }
})

testthat::test_that("partial DR admission has no global treatment by upstream row product", {
  a <- res_case(); fit <- res_fit(a, "dr")
  # Arithmetic-only catalog dimensions for the documented G=T=190,S=D=380
  # example. This record is never fitted, evaluated or presented as a backend.
  # Its reader is real; only metadata counts enter the pure envelope helper.
  metadata <- list(reader_ref = fit$reader_ref, control = fit$control,
    accounting = fit$accounting, outcome_recipe = fit$outcome_recipe,
    treatment_recipe = fit$treatment_recipe, uses_selection_models = FALSE,
    regression_specs = rep(list(list(support_ids = 2L)), 190L),
    supports = rep(list(list()), 380L), treatments = rep(list(list()), 190L),
    treatment_specs = rep(list(list(support_ids = c(1L, 2L))), 190L),
    links = list(list(link_id = "partial", regression = 1L, model = 1L, support_id = 1L)),
    outcome_groups = list(groups = rep(list(list(keep = 1L)), 190L)))
  pre <- res_ns("dr_separate_if_envelope")(metadata, 1L, 1L, FALSE)
  profile <- res_ns("chained_public_profile_defaults")()
  testthat::expect_lt(pre$scratch * 8, profile$numeric_workspace_bytes / 8)
  testthat::expect_identical(pre$gids, 1L)
  testthat::expect_identical(pre$tids, 1L)
  testthat::expect_identical(pre$sids, c(1L, 2L))
})

testthat::test_that("long link IDs and unused dated source columns retain sufficient admission", {
  for (kind in c("dr", "imputation", "augmented")) {
    a <- res_case(); frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
    long_id <- paste0("long-", strrep("q", 4096L))
    extra <- res_ns("construct_ipw_link_view")(frame, 2, 11, 11, c(23, Inf), control = a$control,
      link_id = long_id, functional_key = paste0("cmpcat:", long_id))
    a$views <- c(a$views, list(extra)); fit <- res_fit(a, kind)
    for (details in c(FALSE, TRUE)) {
      pre <- res_budget(fit, kind, 3L, c(4L, 2L), details)
      out <- res_eval(fit, kind, c(10L, 1L, 3L), c(4L, 2L), details)
      testthat::expect_lte(unname(res_ns("cp04_payload")(out)["entries"]), pre$output)
    }
  }
  narrow <- res_case(dated = TRUE, route_kind = "imputation")
  wide <- res_case(dated = TRUE, route_kind = "imputation", extra_reference = TRUE)
  f1 <- res_fit(narrow, "imputation"); f2 <- res_fit(wide, "imputation")
  p1 <- res_budget(f1, "imputation", 3L, c(3L, 1L), FALSE)
  p2 <- res_budget(f2, "imputation", 3L, c(3L, 1L), FALSE)
  testthat::expect_gt(p2$scratch, p1$scratch)
  o1 <- res_eval(f1, "imputation", c(10L, 1L, 3L), c(3L, 1L), FALSE)
  o2 <- res_eval(f2, "imputation", c(10L, 1L, 3L), c(3L, 1L), FALSE)
  testthat::expect_equal(o1$phi, o2$phi, tolerance = 1e-10)
})


testthat::test_that("metadata admission counts repeated treatment-support incidences before concatenation", {
  a <- res_case(); fit <- res_fit(a, "dr")
  # Metadata-only envelope witness, never fitted or evaluated as a backend:
  # 500 models share 500 controls. The temporary concatenation has 250,500
  # entries although the catalog has only 1,000 distinct supports.
  metadata <- list(reader_ref = fit$reader_ref, control = fit$control,
    accounting = list(owned_entries = 0, borrowed_entries = 0, work_ops_bound = 0),
    outcome_recipe = fit$outcome_recipe, treatment_recipe = fit$treatment_recipe,
    uses_selection_models = FALSE, regression_specs = list(list(support_ids = 1:500)),
    supports = rep(list(list()), 1000L), treatments = rep(list(list()), 500L),
    treatment_specs = lapply(501:1000, function(s) list(support_ids = c(s, 1:500))),
    links = lapply(seq_len(500L), function(j) list(link_id = paste0("partial-", j),
      regression = 1L, model = j, support_id = j + 500L)),
    outcome_groups = list(groups = list(list(keep = 1L))))
  metadata$control$max_matrix_entries <- 200000
  testthat::expect_error(res_ns("dr_separate_if_envelope")(
    metadata, 1L, seq_len(500L), FALSE), class = "cdid_ipw_owner_bound_error")
})
