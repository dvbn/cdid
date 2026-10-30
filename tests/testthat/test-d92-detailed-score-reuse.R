# test_file() does not source sibling test files. Import only the required
# helper assignments, without sourcing their fixture setup or test bodies.
local({
  target <- parent.env(environment())
  helpers <- new.env(parent = target)
  helpers$test_path <- testthat::test_path
  providers <- list(
    "test-chained-augmented-backend.R" = c("aug_ns", "aug_has_installed_cdid"),
    "test-chained-conditional-pretest.R" = c("cpt_path", "cpt_cache", "cpt_rds",
      "cpt_csv", "cpt_signs", "cpt_cases", "cpt_defaults", "cpt_fixed",
      "cpt_case_data", "cpt_state")
  )
  for (file in names(providers)) {
    definitions <- parse(testthat::test_path(file), keep.source = FALSE)
    for (name in providers[[file]]) {
      hit <- which(vapply(definitions, function(x) is.call(x) &&
        identical(x[[1L]], quote(`<-`)) && identical(x[[2L]], as.name(name)), logical(1L)))
      stopifnot(length(hit) == 1L)
      eval(definitions[[hit]], envir = helpers)
      assign(name, get(name, envir = helpers, inherits = FALSE), envir = target)
    }
  }
})

# Installed-namespace regression checks use those existing aug_* and cpt_* helpers.
# Only copies of installed functions are edited. The old path restores both the
# unconditional helper read and the original alpha-then-score detailed loop.
d92_core <- function(kind, old = FALSE, reads = NULL) {
  fn <- removeSource(aug_ns(paste0(kind, "_evaluate_if_chunk_core")))
  b <- as.list(body(fn))
  at <- function(name) which(vapply(b, function(x) is.call(x) &&
    identical(x[[1L]], quote(`<-`)) && identical(x[[2L]], as.name(name)), logical(1L)))
  if (old) {
    a <- b[[at("alpha_if")]][[3L]]
    stopifnot(identical(a[[2L]], as.pairlist(alist(cid = , obs_rec = NULL))))
    a[[2L]] <- as.pairlist(alist(cid = ))
    ab <- as.list(a[[3L]])
    k <- which(vapply(ab, identical, logical(1L),
      quote(if (missing(obs_rec)) obs_rec <- read_cell_score(cid))))
    stopifnot(length(k) == 1L)
    ab[[k]] <- quote(obs_rec <- read_cell_score(cid))
    a[[3L]] <- as.call(ab)
    b[[at("alpha_if")]][[3L]] <- a
    db <- b[[length(b) - 1L]][[3L]]
    stopifnot(identical(db[[4L]][[4L]], utils::removeSource(quote({
      obs_rec <- read_cell_score(cid)
      xi_alpha_imp[, cid] <- alpha_if(cid, obs_rec)
      if (!is.null(obs_rec)) cell_scores[obs_rec$rr, cid] <- obs_rec$s_rr
    }))))
    db[[4L]][[4L]] <- quote({
      xi_alpha_imp[, cid] <- alpha_if(cid)
      obs_rec <- read_cell_score(cid)
      if (!is.null(obs_rec)) cell_scores[obs_rec$rr, cid] <- obs_rec$s_rr
    })
    b[[length(b) - 1L]][[3L]] <- db
  }
  if (!is.null(reads)) {
    r <- b[[at("read_cell_score")]][[3L]]
    rb <- as.list(r[[3L]])
    last <- rb[[length(rb)]]
    rb[[length(rb)]] <- substitute({
      record <- LAST
      d92_reads$record_entries <- max(d92_reads$record_entries,
        unname(cp04_payload(record)["entries"]))
      record
    }, list(LAST = last))
    r[[3L]] <- as.call(c(rb[1L], list(quote({
      d92_reads$all <- c(d92_reads$all, cid)
      if (exists("xi_alpha_imp", envir = parent.env(environment()), inherits = FALSE))
        d92_reads$detailed <- c(d92_reads$detailed, cid)
    })), rb[-1L]))
    b[[at("read_cell_score")]][[3L]] <- r
    environment(fn) <- list2env(list(d92_reads = reads), parent = environment(fn))
  }
  body(fn) <- as.call(b)
  fn
}

d92_run <- function(fit, kind, rows, old = FALSE, details = TRUE, scale = "frame") {
  reads <- list2env(list(all = integer(), detailed = integer(), record_entries = 0))
  core <- d92_core(kind, old, reads)
  evaluate <- aug_ns(paste0("evaluate_", kind, "_if_chunk"))
  environment(evaluate) <- list2env(stats::setNames(list(core),
    paste0(kind, "_evaluate_if_chunk_core")), parent = environment(evaluate))
  value <- evaluate(fit, rows, seq_along(fit$links), scale, details)
  list(value = value, reads = reads)
}

d92_fit <- function(case, formula = ~x + z) {
  cd <- cpt_case_data(case)
  ov <- utils::modifyList(cd$cs$args, list(est_method = "dr", xformla = formula))
  state <- cpt_state
  environment(state) <- list2env(list(chained_conditional_did_pretest =
    aug_ns("chained_conditional_did_pretest")), parent = environment(state))
  state(cd$data, ov)$registry$cp04
}

testthat::test_that("D92 preserves complete phi and every detailed coordinate", {
  testthat::skip_if_not(aug_has_installed_cdid(), "Installed cdid namespace required.")
  for (case in c("B2", "B3")) for (formula in list(~x + z, ~1, ~0 + x + z)) {
    kind <- if (case == "B2") "imputation" else "augmented"
    fit <- d92_fit(case, formula)
    st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
    testthat::expect_null(fit$layout)
    testthat::expect_identical(fit$outcome_recipe$intercept, attr(stats::terms(formula), "intercept") == 1L)
    if (!length(all.vars(formula)))
      testthat::expect_true(all(vapply(fit$groups, `[[`, integer(1L), "p") == 0L))
    chunk <- aug_ns("cpt_details_rows")(fit, kind, seq_along(fit$links), 0)
    for (start in seq.int(1L, st$N, by = chunk)) {
      rows <- seq.int(start, min(st$N, start + chunk - 1L))
      fast <- d92_run(fit, kind, rows)
      old <- d92_run(fit, kind, rows, old = TRUE)
      # Whole objects include all nested OR/global-beta, selection, treatment,
      # rate, link-component and nuisance channels, names, ids and normalization.
      testthat::expect_equal(fast$value, old$value, tolerance = 1e-12)
      testthat::expect_identical(fast$value$id, st$id[rows])
      pre <- aug_ns(paste0(kind, "_if_envelope"))(fit, length(rows), details = TRUE)
      testthat::expect_lte(unname(aug_ns("cp04_payload")(fast$value)["entries"]), pre$output)
      # Include the previous outer record while its replacement is evaluated,
      # and alpha's output/subtraction temporaries; payload includes metadata.
      # Reserve 64 further logical entries for helper/argument bookkeeping.
      local_live <- 2 * fast$reads$record_entries +
        4 * unname(aug_ns("cp04_payload")(numeric(length(rows)))["entries"]) + 64
      testthat::expect_lte(local_live, pre$scratch)
      core <- d92_run(fit, kind, rows, scale = "core")$value
      testthat::expect_equal(core$phi, st$lambda * fast$value$phi, tolerance = 1e-12)
      testthat::expect_identical(core$normalization_count, st$n_A)
    }
  }
})

testthat::test_that("D92 counts real score reads including supplied NULL and branch exits", {
  testthat::skip_if_not(aug_has_installed_cdid(), "Installed cdid namespace required.")
  for (case in c("B2", "B3")) for (formula in list(~x + z, ~1, ~0 + x + z)) {
    kind <- if (case == "B2") "imputation" else "augmented"
    fit <- d92_fit(case, formula)
    st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
    # A singleton chunk forces NULL for every other cohort. Where present,
    # off-A rows also force the no-observer return within their own cohort.
    off <- which(st$a_index == 0L)
    null_row <- if (length(off)) off[1L] else 1L
    testthat::expect_true(any(vapply(fit$cells, `[[`, integer(1L), "cohort_index") != st$cohort_index[null_row]))
    chunk <- aug_ns("cpt_details_rows")(fit, kind, seq_along(fit$links), 0)
    for (rows in list(seq_len(chunk), null_row)) {
      fast <- d92_run(fit, kind, rows)
      old <- d92_run(fit, kind, rows, old = TRUE)
      nc <- length(fit$cells)
      intercept <- fit$outcome_recipe$intercept
      testthat::expect_identical(tabulate(fast$reads$detailed, nc), rep(1L, nc))
      testthat::expect_identical(tabulate(old$reads$detailed, nc), rep(if (intercept) 2L else 1L, nc))
      testthat::expect_identical(length(old$reads$all) - length(fast$reads$all), if (intercept) nc else 0L)
      testthat::expect_equal(fast$value, old$value, tolerance = 1e-12)
      if (length(off) && identical(rows, off[1L]))
        testthat::expect_true(all(fast$value$details$imputation_cell_scores == 0))
      a <- d92_run(fit, kind, rows, details = FALSE)
      b <- d92_run(fit, kind, rows, old = TRUE, details = FALSE)
      testthat::expect_identical(a$reads$all, b$reads$all)
      testthat::expect_equal(a$value, b$value, tolerance = 1e-12)
    }
  }
})

testthat::test_that("D92 retains conditional processes, fixed signs and existing cap errors", {
  testthat::skip_if_not(aug_has_installed_cdid(), "Installed cdid namespace required.")
  fixed <- cpt_fixed
  environment(fixed) <- list2env(list(chained_conditional_did_pretest =
    aug_ns("chained_conditional_did_pretest")), parent = environment(fixed))
  for (case in c("B2", "B3")) {
    cd <- cpt_case_data(case)
    signs <- cpt_signs(cd$fx$unit_rademacher_signs)[seq_len(9L), , drop = FALSE]
    ov <- utils::modifyList(cd$cs$args, list(est_method = "dr"))
    seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    fast <- fixed(cd$data, ov, signs)
    old <- local({
      testthat::local_mocked_bindings(
        imputation_evaluate_if_chunk_core = d92_core("imputation", old = TRUE),
        augmented_evaluate_if_chunk_core = d92_core("augmented", old = TRUE), .package = "cdid")
      fixed(cd$data, ov, signs)
    })
    testthat::expect_identical(get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE), seed)
    testthat::expect_identical(fast$chained$status, "available")
    for (nm in c("CvM", "CvMb", "CvMpval", "CvMcval"))
      testthat::expect_equal(fast[[nm]], old[[nm]], tolerance = 1e-12)
    a <- attr(fast, "cells"); b <- attr(old, "cells")
    testthat::expect_gt(length(a), 0L)
    testthat::expect_equal(a, b, tolerance = 1e-12)
    for (j in seq_along(a)) {
      testthat::expect_equal(aug_ns("cpt_cell_process_rows")(a[[j]]),
        aug_ns("cpt_cell_process_rows")(b[[j]]), tolerance = 1e-12)
      testthat::expect_identical(a[[j]]$mult, b[[j]]$mult)
    }
    fit <- attr(fast, "state")$registry$cp04
    kind <- if (case == "B2") "imputation" else "augmented"
    links <- seq_along(fit$links)
    pre <- aug_ns(paste0(kind, "_if_envelope"))(fit, 1L, links, TRUE)
    own <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
    retained <- fit$control$max_matrix_entries - own - pre$output - pre$scratch
    testthat::expect_gte(retained, 0)
    testthat::expect_identical(aug_ns("cpt_details_rows")(fit, kind, links, retained), 1L)
    testthat::expect_error(aug_ns("cpt_details_rows")(fit, kind, links, retained + 1),
      class = "cdid_ipw_owner_bound_error")
    evaluate <- aug_ns(paste0("evaluate_", kind, "_if_chunk"))
    testthat::expect_error(evaluate(fit, c(1L, 1L), details = TRUE), class = paste0("cdid_", kind, "_input_error"))
    testthat::expect_error(evaluate(list(), 1L, details = TRUE), class = paste0("cdid_", kind, "_integrity_error"))
    testthat::expect_error(evaluate(fit, seq_len(fit$control$max_chunk_rows + 1L), details = TRUE),
      class = "cdid_ipw_owner_bound_error")
  }
})
