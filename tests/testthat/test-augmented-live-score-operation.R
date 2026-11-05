v58_ns <- function(nm) getFromNamespace(nm, "cdid")
v58_e <- new.env(parent = asNamespace("cdid"))
for (x in parse(file = testthat::test_path("test-augmented-adjoint-replay.R"))) if (is.call(x) &&
    identical(x[[1L]], as.name("<-"))) eval(x, v58_e)

# Author policy: private namespace/body/lexical-frame mutation authenticity
# is DEVELOPMENT-only. The old frozen V58 test retains that historical gate;
# production tests retain actual issued data, shapes, source-reader and math.
v58_probe <- function(current, reader, adj, ids, core) {
  reference <- current(reader, adj, ids)
  testthat::expect_type(reference, "double")
  e <- environment(current)
  view <- get("issued_view", e, inherits = FALSE)
  state <- get("issued_state", e, inherits = FALSE)
  contract <- get("contract", e, inherits = FALSE)
  testthat::expect_equal(contract(view, state, adj, v58_ns("reconstruction_sparse_adjoint_ids")(
    get("complete", e, inherits = FALSE), ids)), reference, tolerance = 1e-12)
  cache <- view$cache
  k <- which(!vapply(cache, is.null, logical(1L)))[1L]
  testthat::expect_false(is.na(k))
  changed <- view; changed$cache[[k]]$s_rr <- changed$cache[[k]]$s_rr + 1
  testthat::expect_null(contract(changed, state, adj, ids))
  changed <- state; cid <- state$cids[1L]; changed$cells[[cid]]$mean_d <- changed$cells[[cid]]$mean_d + 1
  testthat::expect_null(contract(view, changed, adj, ids))
  changed <- state; changed$B <- changed$B + 1L
  testthat::expect_null(contract(view, changed, adj, ids))
  testthat::expect_null(current(reader, adj, as.integer(length(state$cells) + 1L)))
  # A different supplied callback still takes its checked old path.
  bad <- reader; body(bad) <- quote(NULL)
  testthat::expect_null(current(bad, adj, ids))
  bad <- reader; attr(bad, "private") <- TRUE
  testthat::expect_null(current(bad, adj, ids))
  testthat::expect_equal(current(reader, adj, ids), reference, tolerance = 1e-12)
  invisible(TRUE)
}

v58_run <- function(fit, links, fast = TRUE, probe = FALSE) {
  plan <- v58_ns("chained_full_if_plan")(fit, links, fit$control$max_chunk_rows)
  on.exit(v58_ns("augmented_adjoint_replay_release")(plan), add = TRUE)
  testthat::expect_type(plan, "list")
  access <- v58_ns("augmented_adjoint_replay_access")
  attempts <- hits <- probes <- 0L
  testthat::local_mocked_bindings(augmented_adjoint_replay_access = function(...) {
    op <- access(...); if (is.null(op)) return(NULL)
    old <- op$get; tapped <- FALSE
    op$get <- function(...) {
      out <- old(...)
      if (!tapped) for (core in rev(sys.frames())) {
        if (!exists("eval_saved_adjoint", core, inherits = FALSE) ||
            !exists("coarse_operation", core, inherits = FALSE)) next
        native <- get("coarse_operation", core, inherits = FALSE)
        if (is.null(native)) next
        current <- native$contract_current
        native$contract_current <- function(reader, adj, ids) {
          attempts <<- attempts + 1L
          val <- if (fast) current(reader, adj, ids) else NULL
          if (!is.null(val)) {
            hits <<- hits + 1L
            if (probe && !probes) { v58_probe(current, reader, adj, ids, core); probes <<- probes + 1L }
          }
          val
        }
        assign("coarse_operation", native, core); tapped <<- TRUE
        break
      }
      out
    }
    op
  }, .package = "cdid")
  phi <- v58_e$v57_collect(fit, links, plan)
  list(phi = phi, attempts = attempts, hits = hits, probes = probes)
}

testthat::test_that("owned production operation preserves full/tail scores and actual data guards", {
  for (model in c("pair_logit", "independent_period_logit")) {
    a <- v58_e$v57_wide_case(); fit <- v58_e$v57_aug$rap_fit(a, imp = "separate", model = model, mode = "pooled")
    links <- rev(seq_along(fit$links))
    old <- v58_run(fit, links, FALSE)
    new <- v58_run(fit, links, TRUE, probe = TRUE)
    testthat::expect_gt(new$attempts, 0L)
    # Direct existing C6 has no R4.5 private-binding API dependency.
    testthat::expect_gt(new$hits, 0L); testthat::expect_identical(new$probes, 1L)
    testthat::expect_equal(new$phi, old$phi, tolerance = 1e-12)
    testthat::expect_identical(dim(new$phi), dim(old$phi))
    testthat::expect_identical(new$phi == 0, old$phi == 0)
    testthat::expect_equal(crossprod(new$phi), crossprod(old$phi), tolerance = 1e-12)
    testthat::expect_equal(sqrt(diag(crossprod(new$phi))), sqrt(diag(crossprod(old$phi))), tolerance = 1e-12)
    testthat::expect_gt(max(abs(new$phi[!a$A, , drop = FALSE])), 1e-4)
    N <- nrow(a$Y); offset <- length(v58_e$v57_aug$rap_point(fit)) - length(fit$links)
    for (i in c(1L, N)) {
      d <- rep(-1, N); d[i] <- N - 1
      ref <- (v58_e$v57_aug$rap_reference(a, fit, 1 + 1e-6 * d) -
        v58_e$v57_aug$rap_reference(a, fit, 1 - 1e-6 * d)) / 2e-6
      testthat::expect_equal(new$phi[i, ], ref[offset + links], tolerance = 1e-5)
    }
  }
})


testthat::test_that("production score calls do not authenticate private code per arm", {
  a <- v58_e$v57_wide_case()
  fit <- v58_e$v57_aug$rap_fit(a, imp = "separate", model = "pair_logit", mode = "pooled")
  links <- rev(seq_along(fit$links)); reference <- v58_run(fit, links, FALSE)
  sparse <- v58_ns("reconstruction_sparse_adjoint_ids"); calls <- 0L
  testthat::local_mocked_bindings(
    reconstruction_coarse_source = function(...) stop("development-only call-stack authentication ran"),
    reconstruction_sparse_adjoint_ids = function(...) { calls <<- calls + 1L; sparse(...) },
    .package = "cdid")
  actual <- v58_run(fit, links, TRUE)
  testthat::expect_gt(actual$hits, 0L)
  testthat::expect_gte(calls, actual$hits)
  testthat::expect_equal(actual$phi, reference$phi, tolerance = 1e-12)
  testthat::expect_identical(actual$phi == 0, reference$phi == 0)
})
