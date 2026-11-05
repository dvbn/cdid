# V38 reuses the accepted raw finite-law nuisance refits without copying a new
# fixture or running another test file's expressions. Only named helpers load.
d38_ns <- function(name) getFromNamespace(name, "cdid")
d38_fixture <- local({
  e <- new.env(parent = environment())
  wanted <- c("drp_ns", "drp_gc", "drp_case", "drp_reference", "drp_fit", "drp_point")
  for (x in parse(file = testthat::test_path("test-dr-outcome-pooling.R"))) {
    if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
        is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  }
  stopifnot(all(vapply(wanted, exists, logical(1L), envir = e, inherits = FALSE)))
  e
})
# A test-local exact old envelope; numerical execution is the same function.
d38_old <- local({
  e <- new.env(parent = asNamespace("cdid"))
  for (name in c("dr_separate_if_envelope", "evaluate_dr_separate_if_chunk")) {
    f <- d38_ns(name); environment(f) <- e; assign(name, f, e)
  }
  e$dr_conditional_separate_ordinary_if_envelope <- function(backend, B, links, pre, selection_pre, maps) pre
  e
})
d38_maps <- function(fit, links) {
  incidence <- 0
  for (r in fit$regression_specs) incidence <- incidence + length(r$support_ids)
  for (t in fit$treatment_specs) incidence <- incidence + length(t$support_ids)
  for (s in fit$supports) incidence <- incidence + length(s$channels)
  for (t in fit$treatments) incidence <- incidence + length(t$upstream_ids)
  4096 + 32 * (length(fit$regression_specs) + length(fit$supports) + fit$upstream$D +
    length(fit$treatments) + length(links) + incidence)
}
d38_inventory <- function(fit, B, links) {
  pre <- d38_old$dr_separate_if_envelope(fit, B, links, FALSE)
  sel <- d38_ns("selection_if_budget")(fit$selection_backend, B,
    unique(fit$q_of_s[pre$sids]), details = FALSE)
  d38_ns("dr_conditional_separate_ordinary_inventory")(fit, B, links, pre, sel, d38_maps(fit, links))
}

# Observe actual live caller + local-core arrays at returns/products/solves.
# Excluded fields are original sealed owners/factors or input aliases already
# present in the caller; only maxima/scalars survive this test invocation.
d38_payload <- function(fit, rows, links) {
  e <- new.env(parent = asNamespace("cdid")); outer <- NULL
  peak <- 0; core_calls <- solves <- 0L; saw_previous_model <- FALSE
  own <- function(env, exclude) {
    x <- as.list(env, all.names = TRUE); x[intersect(names(x), exclude)] <- NULL; x
  }
  borrowed <- c("backend", "reader", "st", "frame", "ref", "spec", "f", "of", "or_model")
  remember <- function(extra = list()) {
    if (!is.null(outer)) peak <<- max(peak,
      unname(d38_ns("cp04_payload")(c(own(outer, borrowed), extra))["entries"]))
  }
  evaluate <- d38_ns("evaluate_dr_separate_if_chunk"); environment(evaluate) <- e
  envelope <- d38_ns("dr_separate_if_envelope")
  e$dr_separate_if_envelope <- function(...) {
    outer <<- parent.frame()
    envelope(...)
  }
  e$matrix <- function(...) { value <- base::matrix(...); remember(list(new_matrix = value)); value }
  e[["%*%"]] <- function(x, y) {
    value <- base::`%*%`(x, y)
    remember(list(product_x = x, product_y = y, product_result = value)); value
  }
  e$cp04_solve <- function(R, b) {
    value <- d38_ns("cp04_solve")(R, b); solves <<- solves + 1L
    remember(list(solve_rhs = b, solve_result = value)); value
  }
  e$ipw_logit_values <- function(...) {
    value <- d38_ns("ipw_logit_values")(...); remember(list(logit_result = value)); value
  }
  e$cp04_finite <- function(x, label, positive = FALSE) {
    remember(list(finite_input = x)); d38_ns("cp04_finite")(x, label, positive)
  }
  core <- d38_ns("or_selection_separate_if_core")
  ce <- new.env(parent = environment(core)); environment(core) <- ce
  core_own <- function(env) own(env, c("model", "N", "w", "hits", "up", "pre", "links",
    "details", "block_fits", "f", "rec", "sids", "h"))
  ce$matrix <- function(...) {
    value <- base::matrix(...)
    remember(list(core_live = core_own(parent.frame()), new_core_matrix = value)); value
  }
  ce[["%*%"]] <- function(x, y) {
    value <- base::`%*%`(x, y)
    remember(list(core_live = core_own(parent.frame()), product_x = x,
      product_y = y, product_result = value)); value
  }
  for (name in c("or_pooled_wls_gram_solve", "or_wls_gram_solve")) local({
    original <- d38_ns(name)
    assign(name, function(...) {
      value <- original(...); args <- list(...)
      remember(list(core_live = core_own(parent.frame()), gram_rhs = args[[2L]], gram_result = value)); value
    }, ce)
  })
  e$or_selection_separate_if_core <- function(...) {
    core_calls <<- core_calls + 1L
    if (exists("gamma", outer, inherits = FALSE)) saw_previous_model <<- TRUE
    remember()
    value <- core(...); remember(list(returned_core = value)); value
  }
  out <- evaluate(fit, rows, links, "frame", FALSE)
  pre <- d38_ns("dr_separate_if_envelope")(fit, length(rows), links, FALSE)
  testthat::expect_lte(peak, pre$scratch)
  testthat::expect_lte(unname(d38_ns("cp04_payload")(out)["entries"]), pre$output)
  testthat::expect_identical(core_calls, length(pre$gids))
  testthat::expect_identical(solves, length(pre$tids))
  if (length(pre$gids) > 1L) testthat::expect_true(saw_previous_model)
  list(out = out, peak = peak, pre = pre)
}

testthat::test_that("conditional separate ordinary DR accounts local phases and full original-unit results", {
  cases <- list(list(u = "F", model = "pair_logit", pool = "separate", dated = FALSE, p0 = FALSE),
    list(u = "B", model = "pair_logit", pool = "within_cohort", dated = FALSE, p0 = FALSE),
    list(u = "A", model = "pair_logit", pool = "across_cohorts", dated = FALSE, p0 = FALSE),
    list(u = "F", model = "independent_period_logit", pool = "across_cohorts", dated = TRUE, p0 = FALSE),
    list(u = "F", model = "pair_logit", pool = "separate", dated = FALSE, p0 = TRUE))
  for (spec in cases) {
    a <- d38_fixture$drp_case(spec$u, slopes = !spec$p0, dated = spec$dated, extra_reference = spec$dated)
    if (spec$u == "F" && !spec$dated && !spec$p0) {
      frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
      # Repeated comparison consumers deliberately make local link counts
      # uneven; a global nl/ng average is not a valid block peak.
      a$views <- lapply(c(1L, 2L, 3L, 1L, 1L), function(j) a$views[[j]])
      a$views <- lapply(seq_along(a$views), function(j) {
        v <- a$views[[j]]; id <- paste0("budget-", j, strrep("\u03bb\u6f22", 120L))
        d38_ns("construct_ipw_link_view")(frame, v$pair$s, v$pair$t, v$cohorts$g, v$cohorts$K,
          control = a$control, link_id = id, functional_key = paste0("cmpcat:", id))
      })
    }
    fit <- d38_fixture$drp_fit(a, "conditional_mtar", spec$model, spec$pool)
    links <- rev(seq_along(fit$links)); N <- nrow(a$Y)
    if (spec$u == "F" && !spec$dated && !spec$p0)
      attr(links, "source_notes") <- strrep("\u03bb\u6f22", 500L)
    pre <- d38_ns("dr_separate_if_envelope")(fit, 11L, links, FALSE)
    old <- d38_old$dr_separate_if_envelope(fit, 11L, links, FALSE)
    inv <- d38_inventory(fit, 11L, links)
    testthat::expect_identical(pre$work, old$work)
    testthat::expect_identical(d38_ns("dr_separate_if_envelope")(fit, 11L, links, TRUE),
      d38_old$dr_separate_if_envelope(fit, 11L, links, TRUE))
    testthat::expect_equal(inv$output, pre$output); testthat::expect_equal(inv$scratch, pre$scratch)
    testthat::expect_identical(inv$context_metadata, 2048)
    active <- which(vapply(fit$up_map$theta_coords, function(tc) any(tc %in% pre$uids), logical(1)))
    W <- sum(vapply(active, function(m) {
      b <- fit$selection_backend$blocks[[m]]
      (if (b$has_intercept) b$B_m else 0) + b$d_s + length(fit$up_map$theta_coords[[m]])
    }, numeric(1)))
    testthat::expect_equal(inv$selection_columns, W)
    testthat::expect_equal(inv$upstream, 11 * (length(pre$uids) + W))
    testthat::expect_equal(inv$hits, 11 * length(pre$sids) * (fit$outcome_recipe$p + 4))
    testthat::expect_equal(inv$max_block_links,
      max(tabulate(vapply(fit$links[links], `[[`, integer(1), "regression"))))
    testthat::expect_gte(inv$phase_peak, inv$previous_model_state + inv$core_peak)
    testthat::expect_gte(inv$phase_peak, inv$core_return + inv$current_model_peak)
    selection <- d38_ns("selection_if_budget")(fit$selection_backend, 11L, unique(fit$q_of_s[pre$sids]), details = FALSE)
    testthat::expect_identical(inv$selection_peak, selection$output + selection$scratch)
    if (!is.null(attr(links, "source_notes"))) testthat::expect_equal(inv$link_attribute_entries,
      unname(d38_ns("cp04_payload")(attributes(links))["entries"]))
    if (spec$dated) {
      testthat::expect_equal(inv$reference_x_columns, length(a$reference$x_columns))
      testthat::expect_gt(inv$reference_x_columns, length(a$recipe$columns))
      testthat::expect_equal(inv$reference_s_columns,
        length(a$reference$designs[[a$reference$role_slot["s"]]]$columns))
    }
    subset <- d38_inventory(fit, 11L, links[1L])
    testthat::expect_lte(subset$max_block_links, inv$max_block_links)
    testthat::expect_lte(subset$upstream, inv$upstream)
    testthat::expect_equal(d38_fixture$drp_point(fit), d38_fixture$drp_reference(a, fit), tolerance = 2e-8)
    actual <- raw <- detailed <- matrix(0, N, length(links))
    for (first in seq.int(1L, N, by = 11L)) {
      rows <- rev(seq.int(first, min(N, first + 10L)))
      value <- d38_payload(fit, rows, links)$out
      reference <- d38_old$evaluate_dr_separate_if_chunk(fit, rows, links, "frame", FALSE)
      rich <- d38_ns("evaluate_dr_separate_if_chunk")(fit, rows, links, "frame", TRUE)
      actual[rows, ] <- value$phi; raw[rows, ] <- reference$phi; detailed[rows, ] <- rich$phi
    }
    testthat::expect_equal(actual, raw, tolerance = 1e-10)
    testthat::expect_equal(actual, detailed, tolerance = 1e-10)
    testthat::expect_equal(crossprod(actual), crossprod(raw), tolerance = 1e-10)
    offset <- length(d38_fixture$drp_point(fit)) - 2 * length(fit$links) - 4 * length(fit$treatments)
    for (i in c(1L, N)) {
      delta <- rep(-1, N); delta[i] <- N - 1
      derivative <- (d38_fixture$drp_reference(a, fit, 1 + 1e-6 * delta) -
        d38_fixture$drp_reference(a, fit, 1 - 1e-6 * delta)) / 2e-6
      testthat::expect_equal(actual[i, ], derivative[offset + links], tolerance = 1e-5)
    }
    one <- d38_ns("evaluate_dr_separate_if_chunk")(fit, 11:1, links[1L])
    testthat::expect_equal(unname(one$phi), actual[11:1, 1L, drop = FALSE], tolerance = 1e-10)
    small <- d38_ns("dr_separate_if_envelope")(fit, 1L, links, FALSE)
    frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
    testthat::expect_equal(pre$output - small$output,
      10 * (length(links) + 4 + frame$id_byte_metadata$max))
  }
})

testthat::test_that("fresh conditional separate accounting preserves baseline budgets and exact cap and numeric refusals", {
  a <- d38_fixture$drp_case(); fit <- d38_fixture$drp_fit(a, "conditional_mtar"); links <- 1:3
  baseline <- d38_fixture$drp_fit(a)
  for (details in c(FALSE, TRUE)) testthat::expect_identical(
    d38_ns("dr_separate_if_envelope")(baseline, 11L, links, details),
    d38_old$dr_separate_if_envelope(baseline, 11L, links, details))
  pre <- d38_ns("dr_separate_if_envelope")(fit, 11L, links, FALSE)
  old <- d38_old$dr_separate_if_envelope(fit, 11L, links, FALSE)
  testthat::expect_lt(pre$output + pre$scratch, old$output + old$scratch)
  owned <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
  reseal <- function(control = fit$control, treatments = fit$treatments) {
    x <- fit; x["ref"] <- list(NULL); x$control <- control; x$treatments <- treatments
    d38_ns("ipw_owner_seal")(x, "dr_corollary1_backend")
  }
  ctrl <- fit$control; ctrl$max_matrix_entries <- owned + pre$scratch
  testthat::expect_equal(d38_ns("evaluate_dr_separate_if_chunk")(reseal(ctrl), 1:11, links)$phi,
    d38_old$evaluate_dr_separate_if_chunk(fit, 1:11, links)$phi, tolerance = 1e-12)
  ctrl$max_matrix_entries <- ctrl$max_matrix_entries - 1
  testthat::expect_error(d38_ns("evaluate_dr_separate_if_chunk")(reseal(ctrl), 1:11, links),
    class = "cdid_ipw_owner_bound_error")
  ctrl <- fit$control; ctrl$max_output_entries <- pre$output
  testthat::expect_silent(d38_ns("evaluate_dr_separate_if_chunk")(reseal(ctrl), 1:11, links))
  ctrl$max_output_entries <- ctrl$max_output_entries - 1
  testthat::expect_error(d38_ns("evaluate_dr_separate_if_chunk")(reseal(ctrl), 1:11, links),
    class = "cdid_ipw_owner_bound_error")
  for (kind in c("logit", "loading")) {
    tr <- fit$treatments
    if (kind == "logit") tr[[1L]]$coefficients[] <- Inf else tr[[1L]]$partials$beta[2L, 1L] <- Inf
    b <- reseal(treatments = tr)
    actual <- tryCatch(d38_ns("evaluate_dr_separate_if_chunk")(b, 1:11, links), error = identity)
    reference <- tryCatch(d38_old$evaluate_dr_separate_if_chunk(b, 1:11, links), error = identity)
    testthat::expect_s3_class(actual, "cdid_cp04_numeric_error")
    testthat::expect_identical(class(actual), class(reference))
    testthat::expect_identical(conditionMessage(actual), conditionMessage(reference))
  }
  broken <- fit; broken$links[[1L]]$mu <- broken$links[[1L]]$mu + 1
  testthat::expect_error(d38_ns("evaluate_dr_separate_if_chunk")(broken, 1:11, links),
    class = "cdid_dr_corollary1_integrity_error")
})
