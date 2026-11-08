# Source-owned ordinary buffer inventory; the established raw finite-frequency
# references and every existing numerical tolerance remain unchanged.
rob_ns <- function(name) getFromNamespace(name, "cdid")
rob_helpers <- function(file, wanted) {
  e <- new.env(parent = asNamespace("cdid"))
  for (x in parse(file = testthat::test_path(file))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  stopifnot(all(vapply(wanted, exists, logical(1), envir = e, inherits = FALSE)))
  e
}
rob_pure <- rob_helpers("test-reconstructed-outcome-pooling.R",
  c("rcp_ns", "rcp_gc", "rcp_case", "rcp_fit", "rcp_reference", "rcp_point"))
rob_aug <- rob_helpers("test-augmented-outcome-pooling.R",
  c("rap_ns", "rap_gc", "rap_case", "rap_fit", "rap_reference", "rap_point"))
rob_old <- local({
  e <- new.env(parent = asNamespace("cdid"))
  for (name in c("imputation_if_envelope", "imputation_if_budget", "evaluate_imputation_if_chunk",
      "augmented_if_envelope", "augmented_if_budget", "evaluate_augmented_if_chunk")) {
    f <- rob_ns(name); environment(f) <- e; assign(name, f, e)
  }
  e$reconstruction_ordinary_if_envelope <- function(backend, B, links, pre, width, metadata, selection_pre = NULL) pre
  e
})
rob_fun <- function(kind, suffix) paste0(if (kind == "augmented") "augmented" else "imputation", suffix)
rob_inventory <- function(fit, kind, B, links) {
  pre <- rob_old[[rob_fun(kind, "_if_envelope")]](fit, B, links, FALSE)
  sel <- if (kind == "augmented") rob_ns("selection_if_budget")(
    fit$selection_backend, B, pre$q_ids, details = FALSE) else NULL
  rob_ns("reconstruction_ordinary_inventory")(fit, B, links, pre, sel)
}

# Observe core and nested construction frames. Immutable source aliases are
# excluded, including source fields behind V40 state/zero-template snapshots;
# cached row records, maps, RHS/solve arrays and their repeated views stay in.
rob_payload <- function(fit, kind, rows, links) {
  e <- new.env(parent = asNamespace("cdid")); outer <- NULL; peak <- 0
  model_shapes <- list(); saw_previous <- FALSE; selection_shapes <- list()
  borrowed <- c("st", "X_F", "core_fit", "backend", "reader", "frame", "recipe", "trt_recipe",
    "cf", "gf", "cf_g", "gf_g", "cf_k", "gf_k", "lf", "tf", "tif", "sh", "bf", "fit", "cells", "groups")
  strip <- function(x) {
    if (is.list(x)) {
      if (any(startsWith(class(x), "cdid_"))) return(NULL)
      x[intersect(names(x), c("core", "cells", "groups", "group"))] <- NULL
      return(lapply(x, strip))
    }
    x
  }
  own <- function(f) {
    x <- as.list(f, all.names = TRUE); x[intersect(names(x), borrowed)] <- NULL
    strip(x)
  }
  remember <- function(local = NULL, extra = list()) {
    if (is.null(outer)) return(invisible(NULL))
    frames <- list(core = own(outer), local = if (!is.null(local) && !identical(local, outer)) own(local) else NULL)
    peak <<- max(peak, unname(rob_ns("cp04_payload")(c(frames, strip(extra)))["entries"]))
    if (exists("trt_if", outer, inherits = FALSE)) {
      models <- get("trt_if", outer, inherits = FALSE)
      model_shapes <<- lapply(models, function(m) if (is.null(m)) NULL else lapply(m, dim))
    }
    if (exists("comp", outer, inherits = FALSE) || exists("comp_ipw", outer, inherits = FALSE)) saw_previous <<- TRUE
  }
  for (name in c("matrix", "numeric", "rep", "cbind", "t", "sweep")) local({
    nm <- name; original <- get(nm, baseenv())
    assign(nm, function(...) { value <- original(...); remember(parent.frame(), list(new_value = value)); value }, e)
  })
  e[["%*%"]] <- function(x, y) {
    value <- base::`%*%`(x, y); remember(parent.frame(), list(x = x, y = y, product = value)); value
  }
  for (name in c("cp04_solve", "or_pooled_wls_gram_solve", "or_wls_gram_solve", "ipw_logit_values")) local({
    nm <- name; original <- rob_ns(nm)
    assign(nm, function(...) { value <- original(...); remember(parent.frame(), list(result = value)); value }, e)
  })
  factory <- rob_ns("reconstruction_cached_cell_reader")
  e$reconstruction_cached_cell_reader <- function(read, plan, n_cells) {
    outer <<- environment(read); factory(read, plan, n_cells)
  }
  upstream <- rob_ns("evaluate_selection_upstream_chunk")
  e$evaluate_selection_upstream_chunk <- function(...) {
    outer <<- parent.frame()
    value <- upstream(...); selection_shapes <<- lapply(value$block_ifs, function(b) if (is.null(b)) NULL else lapply(b, dim))
    remember(parent.frame(), list(up_res = value)); value
  }
  for (name in c("reconstruction_score_cache_view", "reconstruction_complete_cache_witness", "reconstruction_sparse_adjoint_ids")) local({
    nm <- name; original <- rob_ns(nm)
    assign(nm, function(...) { value <- original(...); remember(parent.frame(), list(cache_or_closure = value)); value }, e)
  })
  core <- rob_ns(if (kind == "augmented") "augmented_evaluate_if_chunk_core" else "imputation_evaluate_if_chunk_core")
  environment(core) <- e
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  pre <- rob_ns(rob_fun(kind, "_if_envelope"))(fit, length(rows), links, FALSE)
  pair_reader <- function(r, pair, read_outcomes = TRUE) rob_ns("cp04_pair_rows")(reader, r, pair, read_outcomes)
  args <- list(st = reader$state, X_F = if (kind == "augmented") rob_ns("augmented_backend_X_F")(fit) else
    rob_ns("imputation_backend_X_F")(fit, reader), pair_reader = pair_reader,
    rows = rows, links = links, scale = "core", details = FALSE, sets = pre)
  args <- if (kind == "augmented") c(args, list(backend = fit)) else c(args,
    list(core_fit = fit, recipe = fit$outcome_recipe, or_frame_index = fit$or_frame_index,
      trt_recipe = fit$treatment_recipe, est_method = fit$est_method))
  out <- do.call(core, args)
  testthat::expect_lte(peak, pre$scratch)
  for (model in model_shapes) if (!is.null(model)) {
    testthat::expect_length(model, 5L)
    for (shape in model) testthat::expect_equal(shape, c(length(rows), fit$treatment_recipe$p))
  }
  list(out = out, peak = peak, previous_link_locals = saw_previous, selection_shapes = selection_shapes)
}

testthat::test_that("generic ordinary reconstruction preserves full F numerical and affine contracts", {
  specs <- list(list(kind = "imputation", method = "ipw", imp = "separate", p0 = FALSE, intercept = TRUE),
    list(kind = "imputation", method = "reg", imp = "within_cohort", p0 = FALSE, intercept = TRUE),
    list(kind = "imputation", method = "dr", imp = "across_cohorts", p0 = FALSE, intercept = TRUE),
    list(kind = "imputation", method = "dr", imp = "separate", p0 = TRUE, intercept = TRUE),
    list(kind = "imputation", method = "dr", imp = "within_cohort", p0 = FALSE, intercept = FALSE),
    list(kind = "augmented", model = "pair_logit", imp = "within_cohort", p0 = FALSE, intercept = TRUE),
    list(kind = "augmented", model = "independent_period_logit", imp = "across_cohorts", p0 = FALSE, intercept = TRUE))
  for (s in specs) {
    h <- if (s$kind == "augmented") rob_aug else rob_pure
    a <- h[[if (s$kind == "augmented") "rap_case" else "rcp_case"]](intercept = s$intercept, slopes = !s$p0,
      route_kind = s$kind)
    fit <- if (s$kind == "augmented") h$rap_fit(a, imp = s$imp, model = s$model, mode = "pooled") else
      h$rcp_fit(a, method = s$method, pooling = s$imp, mode = "pooled")
    envelope <- rob_ns(rob_fun(s$kind, "_if_envelope")); old_envelope <- rob_old[[rob_fun(s$kind, "_if_envelope")]]
    evaluate <- rob_ns(paste0("evaluate_", s$kind, "_if_chunk"))
    old_evaluate <- rob_old[[paste0("evaluate_", s$kind, "_if_chunk")]]
    links <- rev(seq_along(fit$links)); N <- nrow(a$Y)
    one <- envelope(fit, 1L, links); two <- envelope(fit, 2L, links)
    for (B in c(1L, 2L, 3L, 7L, 11L)) {
      pre <- envelope(fit, B, links); old <- old_envelope(fit, B, links)
      inv <- rob_inventory(fit, s$kind, B, links)
      testthat::expect_identical(pre$work, old$work)
      testthat::expect_identical(pre$score_cache, old$score_cache)
      testthat::expect_identical(envelope(fit, B, links, TRUE), old_envelope(fit, B, links, TRUE))
      for (field in c("output", "scratch")) testthat::expect_equal(pre[[field]],
        one[[field]] + (B - 1) * (two[[field]] - one[[field]]))
      expected_cache <- 1024 + sum(vapply(pre$cids, function(cid) {
        g <- fit$groups[[fit$cells[[cid]]$group_id]]
        4 * B * (g$p + 2) + 256 + 8 * sum(nchar(g$kept_columns, type = "bytes"))
      }, numeric(1)))
      testthat::expect_equal(inv$cache_full_needed_entries, expected_cache)
      testthat::expect_gte(inv$cache_full_needed_entries, pre$score_cache$entries)
      testthat::expect_equal(inv$treatment_entries - inv$coefficients$treatment[["constant"]],
        5 * B * length(pre$tids) * fit$treatment_recipe$p)
      for (phase in inv$phases) testthat::expect_true(all(inv$coefficients$phase >= phase))
      testthat::expect_lte(pre$output, old$output); testthat::expect_lte(pre$scratch, old$scratch)
    }
    pre <- envelope(fit, 11L, links)
    measured <- rob_payload(fit, s$kind, 11:1, links)
    testthat::expect_true(measured$previous_link_locals)
    if (s$kind == "augmented") {
      for (m in seq_along(measured$selection_shapes)) if (!is.null(measured$selection_shapes[[m]])) {
        bl <- fit$selection_backend$blocks[[m]]
        testthat::expect_equal(measured$selection_shapes[[m]]$xi_alpha,
          c(11L, if (bl$has_intercept) bl$B_m else 0L))
        testthat::expect_equal(measured$selection_shapes[[m]]$xi_beta, c(11L, bl$d_s))
        testthat::expect_equal(measured$selection_shapes[[m]]$xi_theta, c(11L, length(fit$up_map$theta_coords[[m]])))
      }
    }
    actual <- old_phi <- detailed <- matrix(NA_real_, N, length(links))
    for (first in seq.int(1L, N, by = 11L)) {
      rows <- rev(seq.int(first, min(N, first + 10L)))
      value <- evaluate(fit, rows, links, "frame", FALSE)
      testthat::expect_lte(unname(rob_ns("cp04_payload")(value)["entries"]), envelope(fit, length(rows), links)$output)
      actual[rows, ] <- value$phi; old_phi[rows, ] <- old_evaluate(fit, rows, links, "frame", FALSE)$phi
      detailed[rows, ] <- evaluate(fit, rows, links, "frame", TRUE)$phi
    }
    testthat::expect_equal(actual, old_phi, tolerance = 1e-10)
    testthat::expect_equal(actual, detailed, tolerance = 1e-10)
    testthat::expect_equal(crossprod(actual), crossprod(old_phi), tolerance = 1e-10)
    point <- h[[if (s$kind == "augmented") "rap_point" else "rcp_point"]]
    reference <- h[[if (s$kind == "augmented") "rap_reference" else "rcp_reference"]]
    testthat::expect_equal(point(fit), reference(a, fit), tolerance = 2e-8)
    zero <- which(rowSums(a$obs) == 0L)[1L]; single <- which(rowSums(a$obs) == 1L)[1L]
    for (i in c(zero, single)) {
      delta <- rep(-1, N); delta[i] <- N - 1
      # The fourth-order stencil reduces cancellation noise near zero while
      # retaining the independent empirical-frequency refit and tolerance.
      h <- 5e-5
      derivative <- tail((-reference(a, fit, 1 + 2 * h * delta) +
        8 * reference(a, fit, 1 + h * delta) -
        8 * reference(a, fit, 1 - h * delta) +
        reference(a, fit, 1 - 2 * h * delta)) / (12 * h), length(fit$links))
      testthat::expect_equal(actual[i, ], derivative[links], tolerance = 1e-5)
    }
    owned <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
    reseal <- function(control) {
      b <- fit; b$control <- control
      rob_ns("ipw_owner_seal")(b, if (s$kind == "augmented") "augmented_backend" else "imputation_backend")
    }
    ctrl <- fit$control; ctrl$max_matrix_entries <- owned + pre$scratch
    testthat::expect_equal(evaluate(reseal(ctrl), 1:11, links)$phi, old_evaluate(fit, 1:11, links)$phi, tolerance = 1e-12)
    ctrl$max_matrix_entries <- ctrl$max_matrix_entries - 1
    testthat::expect_error(evaluate(reseal(ctrl), 1:11, links), class = "cdid_ipw_owner_bound_error")
    ctrl <- fit$control; ctrl$max_output_entries <- pre$output
    testthat::expect_silent(evaluate(reseal(ctrl), 1:11, links))
    ctrl$max_output_entries <- ctrl$max_output_entries - 1
    testthat::expect_error(evaluate(reseal(ctrl), 1:11, links), class = "cdid_ipw_owner_bound_error")
    if (length(fit$treatments)) {
      broken <- fit; broken$treatments[[1L]]$coefficients[] <- Inf
      broken <- rob_ns("ipw_owner_seal")(broken,
        if (s$kind == "augmented") "augmented_backend" else "imputation_backend")
      current_error <- tryCatch(evaluate(broken, 1:11, links), error = identity)
      old_error <- tryCatch(old_evaluate(broken, 1:11, links), error = identity)
      testthat::expect_s3_class(current_error, "error")
      testthat::expect_identical(class(current_error), class(old_error))
      testthat::expect_identical(conditionMessage(current_error), conditionMessage(old_error))
    }
  }
})

testthat::test_that("cache selection and affine admission stay independent, including long consumers", {
  a <- rob_pure$rcp_case(route_kind = "imputation")
  frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
  views <- a$views
  a$views <- lapply(seq_along(c(1L, 2L, 3L, 1L, 1L)), function(j) {
    v <- views[[c(1L, 2L, 3L, 1L, 1L)[j]]]
    id <- paste0("ordinary-budget-", j, strrep("\u03bb\u6f22", 200L))
    rob_ns("construct_ipw_link_view")(frame, v$pair$s, v$pair$t, v$cohorts$g, v$cohorts$K,
      control = a$control, link_id = id, functional_key = paste0("cmpcat:", id))
  })
  fit <- rob_pure$rcp_fit(a, method = "dr", pooling = "across_cohorts", mode = "pooled")
  links <- seq_along(fit$links); attr(links, "source_notes") <- strrep("\u03bb\u6f22", 500L)
  envelope <- rob_ns("imputation_if_envelope")
  genuine <- rob_ns("reconstruction_score_cache_plan")
  inventory <- rob_ns("reconstruction_ordinary_inventory")
  for (cache in c("partial", "declined")) {
    scope <- function() {
      testthat::local_mocked_bindings(reconstruction_score_cache_plan = function(...) {
        plan <- genuine(...)
        if (cache == "declined") list(enabled = FALSE, cids = integer(0L), entries = 0) else {
          plan$cids <- head(plan$cids, 1L); plan
        }
      }, .package = "cdid")
      one <- envelope(fit, 1L, links); two <- envelope(fit, 2L, links)
      for (B in c(1L, 2L, 3L, 7L, 11L)) {
        pre <- envelope(fit, B, links); inv <- inventory(fit, B, links, pre)
        for (field in c("output", "scratch")) testthat::expect_equal(pre[[field]],
          one[[field]] + (B - 1) * (two[[field]] - one[[field]]))
        testthat::expect_equal(inv$cache_full_needed_entries, rob_inventory(fit, "imputation", B, links)$cache_full_needed_entries)
      }
      actual <- rob_ns("evaluate_imputation_if_chunk")(fit, 11:1, links)$phi
      testthat::expect_equal(actual, rob_old$evaluate_imputation_if_chunk(fit, 11:1, links)$phi, tolerance = 1e-10)
    }
    scope()
  }
  # If an optional replacement is larger at B=1 or has a larger slope, keep
  # the old line for every B, rather than choosing a B-dependent minimum.
  scope <- function() {
    testthat::local_mocked_bindings(reconstruction_ordinary_inventory = function(...) {
      inv <- inventory(...)
      inv$coefficients$output[] <- 1e12; inv$coefficients$scratch[] <- 1e12
      inv
    }, .package = "cdid")
    for (B in c(1L, 2L, 11L)) testthat::expect_identical(envelope(fit, B, links),
      rob_old$imputation_if_envelope(fit, B, links))
  }
  scope()
})

testthat::test_that("positive-row affine dominance admits larger constants without crossing-line choices", {
  a <- rob_pure$rcp_case(route_kind = "imputation")
  fit <- rob_pure$rcp_fit(a, method = "ipw", mode = "pooled")
  links <- seq_along(fit$links)
  first <- rob_old$imputation_if_envelope(fit, 1L, links)
  second <- rob_old$imputation_if_envelope(fit, 2L, links)
  coefficients <- lapply(c("output", "scratch"), function(field)
    c(constant = 2 * first[[field]] - second[[field]], row = second[[field]] - first[[field]]))
  names(coefficients) <- c("output", "scratch")
  testthat::expect_true(all(vapply(coefficients, function(x) x[["constant"]] > 3 && x[["row"]] > 3, logical(1))))
  deltas <- list(
    larger_constant = list(delta = c(constant = 2, row = -3), selected = TRUE),
    equal_at_one = list(delta = c(constant = 2, row = -2), selected = TRUE),
    crosses_after_one = list(delta = c(constant = 3, row = -2), selected = FALSE),
    larger_slope = list(delta = c(constant = -3, row = 2), selected = FALSE))
  frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
  width <- coefficients$output[["row"]] - 16 - frame$id_byte_metadata$max
  choose <- rob_ns("reconstruction_ordinary_if_envelope")
  for (case in deltas) {
    candidate <- lapply(coefficients, function(x) x + case$delta)
    scope <- function() {
      testthat::local_mocked_bindings(reconstruction_ordinary_inventory = function(backend, B, links, pre, selection_pre = NULL)
        list(output = candidate$output[["constant"]] + B * candidate$output[["row"]],
          scratch = candidate$scratch[["constant"]] + B * candidate$scratch[["row"]], coefficients = candidate),
        .package = "cdid")
      seen <- list()
      for (B in c(1L, 2L, 3L, 7L, 11L)) {
        old <- rob_old$imputation_if_envelope(fit, B, links)
        metadata <- coefficients$output[["constant"]] + B * frame$id_byte_metadata$max
        got <- choose(fit, B, links, old, width, metadata)
        seen[[as.character(B)]] <- got
        for (field in c("output", "scratch")) {
          expected <- if (case$selected) candidate[[field]][["constant"]] + B * candidate[[field]][["row"]] else old[[field]]
          testthat::expect_equal(got[[field]], expected)
          if (case$selected) testthat::expect_lte(got[[field]], old[[field]])
        }
        unchanged <- setdiff(names(old), c("output", "scratch"))
        testthat::expect_identical(got[unchanged], old[unchanged])
      }
      for (B in c(3L, 7L, 11L)) for (field in c("output", "scratch"))
        testthat::expect_equal(seen[[as.character(B)]][[field]],
          seen[["1"]][[field]] + (B - 1) * (seen[["2"]][[field]] - seen[["1"]][[field]]))
    }
    scope()
    for (field in c("output", "scratch")) {
      # This identity is the proof used by the production choice, including
      # equality at one. Neither a lower value at a later B nor at one alone
      # licenses a line whose positive-row domain is not wholly dominated.
      gap_one <- sum(coefficients[[field]]) - sum(candidate[[field]])
      gap_slope <- coefficients[[field]][["row"]] - candidate[[field]][["row"]]
      testthat::expect_identical(gap_one >= 0 && gap_slope >= 0, case$selected)
      for (B in c(1L, 2L, 11L, 100000L)) testthat::expect_equal(
        coefficients[[field]][["constant"]] + B * coefficients[[field]][["row"]] -
          candidate[[field]][["constant"]] - B * candidate[[field]][["row"]], gap_one + (B - 1) * gap_slope)
    }
  }
})

testthat::test_that("separate-outcome and reference-layout reconstruction budgets remain unchanged", {
  for (kind in c("imputation", "augmented")) {
    h <- if (kind == "augmented") rob_aug else rob_pure
    make <- h[[if (kind == "augmented") "rap_case" else "rcp_case"]]
    for (dated in c(FALSE, TRUE)) {
      a <- make(dated = dated, route_kind = kind, extra_reference = dated)
      fit <- if (kind == "augmented") h$rap_fit(a, mode = if (dated) "pooled" else "separate") else
        h$rcp_fit(a, mode = if (dated) "pooled" else "separate")
      current <- rob_ns(rob_fun(kind, "_if_envelope")); old <- rob_old[[rob_fun(kind, "_if_envelope")]]
      for (B in c(1L, 2L, 11L)) for (details in c(FALSE, TRUE))
        testthat::expect_identical(current(fit, B, 1:3, details), old(fit, B, 1:3, details))
    }
  }
})
