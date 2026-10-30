# Explicit private enclosure tests retained unchanged under D23. Public
# conventional full-null tests are in test-cp02-d23-public.R.
# D21 full declared unconditional pretrend Wald test (R/chained_pretest.R).
# Deterministic helper tests on hand matrices, a formula mutant, refusal cases
# and the already computed frozen N200 reporting fields (no fit, DGP, RNG or
# bootstrap). Requires the copied native DLL for the certified eigen/residual
# primitives; did is not needed except for the final MP-binding test.

pt_locate_fixture <- function(filename) {
  candidates <- c(tryCatch(testthat::test_path("fixtures", filename), error = function(e) character(0L)),
    file.path("fixtures", filename), file.path("tests", "testthat", "fixtures", filename), file.path("..", "..", "tests", "testthat", "fixtures", filename))
  for (p in candidates) if (nzchar(p) && file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  stop("Fixture file not found: ", filename)
}

pt_meta <- function(roles, available = rep(TRUE, length(roles)), status = NULL) {
  data.frame(g = 3, time = seq_along(roles), role = roles, target_id = paste0(roles, "(g=3,t=", seq_along(roles), ")"),
    available = available, support_status = if (is.null(status)) ifelse(available, "supported", "unavailable_pretrend") else status,
    support_reason = "fixture", stringsAsFactors = FALSE)
}

# Hand case: theta = (0.3, -0.2), V = [[0.04, 0.01], [0.01, 0.09]].
# det V = 0.0036 - 0.0001 = 0.0035; theta' adj(V) theta = 0.09*0.09 + 0.04*0.04 + 2*0.01*0.06 = 0.0081 + 0.0016 + 0.0012 = 0.0109
# W = 0.0109 / 0.0035 = 3.1142857142857...; df = 2; p = exp(-W/2) = 0.21075...
pt_hand <- list(theta = c(0.3, -0.2), V = matrix(c(0.04, 0.01, 0.01, 0.09), 2L), W = 0.0109 / 0.0035)

test_that("the declared nonreference set excludes structural reference zeros and keeps unavailable pre cells declared", {
  tm <- pt_meta(c("pre", "reference", "pre", "post", "post"), available = c(TRUE, TRUE, FALSE, TRUE, TRUE))
  s <- chained_pretest_declared_set(tm)
  expect_identical(s$pre, c(1L, 3L)); expect_identical(s$reference, 2L)
  expect_identical(s$ids, tm$target_id[c(1, 3)]); expect_identical(s$available, c(TRUE, FALSE))
  expect_error(chained_pretest_declared_set(data.frame(g = 1)), class = "cdid_chained_pretest_input_error")
})

test_that("the full-set Wald statistic matches hand arithmetic and the did-style n theta' Omega^-1 theta form", {
  h <- pt_hand
  r <- chained_pretest_wald(h$theta, h$V, h$theta, h$theta, h$V, h$V, c("pre(g=3,t=1)", "pre(g=3,t=2)"), 1e-6)
  expect_identical(r$status, "computed"); expect_identical(r$df, 2L)
  expect_equal(r$W, h$W, tolerance = 1e-13)
  expect_equal(r$Wpval, exp(-h$W / 2), tolerance = 1e-13)                 # chi-square(2) survival function
  expect_true(r$W_interval[1L] <= r$W && r$W <= r$W_interval[2L])
  expect_true(diff(r$W_interval) < 1e-12)
  expect_true(r$p_interval_evaluated[1L] <= r$Wpval && r$Wpval <= r$p_interval_evaluated[2L])
  expect_true(r$certified); expect_match(r$certification_scope, "not certified CDF bounds")
  expect_true(is.finite(r$arithmetic$width)); expect_null(r$normalization)
  expect_match(r$description, "full declared nonreference pretrend Wald test", fixed = TRUE); expect_match(r$not_conditional_pretest, "cdid::conditional_did_pretest()", fixed = TRUE)
  # did 2.5.1 writes W = n * t(pre) %*% solve(n * V) %*% pre: the same number for any n
  n <- 200
  expect_equal(r$W, as.double(n * t(h$theta) %*% solve(n * h$V) %*% h$theta), tolerance = 1e-13)
  # a one-dimensional set reduces to z^2
  r1 <- chained_pretest_wald(0.3, matrix(0.04), 0.3, 0.3, matrix(0.04), matrix(0.04), "pre(g=3,t=1)", 1e-6)
  expect_equal(r1$W, 0.3^2 / 0.04, tolerance = 1e-14); expect_identical(r1$df, 1L)
  expect_equal(r1$Wpval, 2 * stats::pnorm(-abs(0.3 / 0.2)), tolerance = 1e-13)
})

test_that("the hand expectation discriminates a mutated formula (A x instead of A^-1 a)", {
  h <- pt_hand
  mutant <- chained_pretest_wald
  b <- body(mutant)
  mutated <- do.call(substitute, list(b, list(solve = as.name("crossprod"))))   # solve(A_sym, a) -> crossprod(A_sym, a) = A_sym a
  expect_false(identical(mutated, b))
  body(mutant) <- mutated
  environment(mutant) <- environment(chained_pretest_wald)
  rm_ <- mutant(h$theta, h$V, h$theta, h$theta, h$V, h$V, c("a", "b"), 1e-6)
  # the mutant still runs (its residual certificate reports an inconsistent solve) or returns another number;
  # either way it cannot reproduce the hand value, so the expectation above is meaningful
  if (identical(rm_$status, "computed")) {
    expect_false(isTRUE(all.equal(rm_$W, h$W, tolerance = 1e-8)))
    expect_equal(rm_$W, as.double(t(h$theta) %*% h$V %*% h$theta), tolerance = 1e-12)   # the mutant computes theta' V theta
  } else {
    expect_identical(rm_$status, "unavailable")                      # residual-certified refusal of the wrong solve
  }
  # the correct implementation and the hand value agree
  expect_equal(chained_pretest_wald(h$theta, h$V, h$theta, h$theta, h$V, h$V, c("a", "b"), 1e-6)$W, h$W, tolerance = 1e-13)
})

test_that("unavailable and not-applicable statuses carry exact reasons, keep W/Wpval NULL and never filter restrictions", {
  h <- pt_hand
  # singular covariance: no pseudoinverse, explicit non-identification disclaimer
  S <- matrix(c(1, 1, 1, 1), 2L)
  rs <- chained_pretest_wald(c(0.1, 0.2), S, c(0.1, 0.2), c(0.1, 0.2), S, S, c("a", "b"), 1e-6)
  expect_identical(rs$status, "unavailable"); expect_null(rs$W); expect_null(rs$Wpval); expect_identical(rs$df, 2L)
  expect_match(rs$reason, "not certified positive definite"); expect_match(rs$reason, "no pseudoinverse")
  expect_match(rs$identification_claim, "not a statement that any ATT is unidentified")
  # certified eigenvalue below the covariance enclosure width: exact covariance not certified nonsingular
  Vw <- h$V; rw <- chained_pretest_wald(h$theta, Vw, h$theta, h$theta, Vw - 0.05, Vw + 0.05, c("a", "b"), 1e-6)
  expect_identical(rw$status, "unavailable"); expect_match(rw$reason, "does not exceed its certified enclosure width")
  # A wide finite PD-compatible enclosure remains diagnostic; no unapproved Wald cutoff.
  rt <- chained_pretest_wald(h$theta, h$V, h$theta - 1e-4, h$theta + 1e-4, h$V, h$V, c("a", "b"), 1e-6)
  expect_identical(rt$status, "computed"); expect_true(rt$arithmetic$width > 1e-6)
  expect_true(rt$W_interval[1L] <= h$W && h$W <= rt$W_interval[2L])
  # The D20 level is not a new rank or Wald/p-value policy.
  expect_identical(chained_pretest_wald(h$theta, h$V, h$theta - 1e-4, h$theta + 1e-4, h$V, h$V, c("a", "b"), 1e-2)$status, "computed")
  # enclosures that do not contain the reported values are refused as inconsistent
  expect_identical(chained_pretest_wald(h$theta, h$V, h$theta + 1, h$theta + 2, h$V, h$V, c("a", "b"), 1e-6)$status, "unavailable")
  # Covariance boxes must contain the reported covariance, not merely be ordered.
  bad <- chained_pretest_wald(h$theta, h$V, h$theta, h$theta, h$V + 1, h$V + 2, c("a", "b"), 1e-6)
  expect_identical(bad$status, "unavailable"); expect_match(bad$reason, "do not contain")
  # unsupported levels refuse (typed), nonfinite inputs are unavailable
  expect_error(chained_pretest_wald(h$theta, h$V, h$theta, h$theta, h$V, h$V, c("a", "b"), 0), class = "cdid_chained_pretest_input_error")
  expect_error(chained_pretest_wald(h$theta, h$V, h$theta, h$theta, h$V, h$V, c("a", "b"), 1), class = "cdid_chained_pretest_input_error")
  expect_identical(chained_pretest_wald(c(0.3, NA), h$V, h$theta, h$theta, h$V, h$V, c("a", "b"), 1e-6)$status, "unavailable")
  # empty declared set
  r0 <- chained_pretest_wald(numeric(0L), matrix(0, 0L, 0L), numeric(0L), numeric(0L), matrix(0, 0L, 0L), matrix(0, 0L, 0L), character(0L), 1e-6)
  expect_identical(r0$status, "not_applicable"); expect_identical(r0$df, 0L)
})

test_that("the reporting-level record tests the FULL declared set: a missing pretrend cell makes the test unavailable, effects untouched", {
  tm <- pt_meta(c("pre", "pre", "post"), available = c(TRUE, FALSE, TRUE))
  est <- c(0.1, NA, 2); names(est) <- tm$target_id
  rep_ <- list(base_period = "varying", target_meta = tm, estimates = est, inference = list(vcov = diag(c(0.01, NA, 0.04))),
    numerical_certificate = list(certified = TRUE, theta = list(lo = matrix(est, 1L), hi = matrix(est, 1L)),
      vcov = list(lo = diag(c(0.01, NA, 0.04)), hi = diag(c(0.01, NA, 0.04)))))
  rec <- chained_pretest_record(rep_, 1e-6)
  expect_identical(rec$status, "unavailable"); expect_null(rec$W); expect_null(rec$Wpval)
  expect_identical(rec$declared_cells, tm$target_id[1:2]); expect_identical(rec$unavailable_cells, tm$target_id[2L])
  expect_match(rec$reason, "no available-only subset is substituted"); expect_identical(rec$n_declared_pretrend, 2L); expect_identical(rec$n_available_pretrend, 1L)
  expect_identical(rec$decision, "D21")
  # structural reference zeros are identities: excluded from the declared set, never tested
  tm2 <- pt_meta(c("reference", "pre", "post"), status = c("deterministic_reference_zero", "supported_direct_pretrend", "estimable"))
  est2 <- c(0, 0.3, 2); names(est2) <- tm2$target_id; V2 <- diag(c(0, 0.04, 0.01))
  rep2 <- list(base_period = "universal", target_meta = tm2, estimates = est2, inference = list(vcov = V2),
    numerical_certificate = list(certified = TRUE, theta = list(lo = matrix(est2, 1L), hi = matrix(est2, 1L)), vcov = list(lo = V2, hi = V2)))
  rec2 <- chained_pretest_record(rep2, 1e-6)
  expect_identical(rec2$status, "computed"); expect_identical(rec2$df, 1L); expect_identical(rec2$cells, tm2$target_id[2L])
  expect_identical(rec2$structural_reference_cells, tm2$target_id[1L]); expect_equal(rec2$W, 0.3^2 / 0.04, tolerance = 1e-14)
  # no declared comparison at all
  tm3 <- pt_meta(c("reference", "post"), status = c("deterministic_reference_zero", "estimable"))
  rec3 <- chained_pretest_record(list(base_period = "universal", target_meta = tm3, estimates = c(0, 1), inference = list(vcov = diag(c(0, 1))),
    numerical_certificate = list(certified = TRUE, theta = list(lo = matrix(c(0, 1), 1L), hi = matrix(c(0, 1), 1L)), vcov = list(lo = diag(c(0, 1)), hi = diag(c(0, 1))))), 1e-6)
  expect_identical(rec3$status, "not_applicable")
  # no certified enclosure: unavailable with the reason, nothing invented
  rec4 <- chained_pretest_record(list(base_period = "varying", target_meta = pt_meta("pre"), estimates = c(`pre(g=3,t=1)` = 0.1),
    inference = list(vcov = matrix(0.01)), numerical_certificate = NULL), 1e-6)
  expect_identical(rec4$status, "unavailable"); expect_match(rec4$reason, "no certified numerical enclosure")
})

test_that("the frozen N200 reporting fields give a computed test equal to the did-style statistic", {
  z <- readRDS(pt_locate_fixture("cp02_result_adapter_N200.rds"))
  rep_ <- z$reporting
  rec <- chained_pretest_record(rep_, 1e-6)
  expect_identical(rec$status, "computed"); expect_identical(rec$df, 1L); expect_identical(rec$cells, "pre(g=3,t=2)")
  pre <- which(rep_$target_meta$role == "pre")
  n <- rep_$inference$n_total
  did_style <- as.double(n * t(rep_$estimates[pre]) %*% solve(rep_$inference$Omega[pre, pre, drop = FALSE]) %*% rep_$estimates[pre])
  expect_equal(rec$W, did_style, tolerance = 1e-12)
  expect_equal(rec$W, (rep_$estimates[[pre]] / rep_$se[[pre]])^2, tolerance = 1e-12)
  expect_equal(rec$Wpval, stats::pchisq(rec$W, 1L, lower.tail = FALSE), tolerance = 1e-14)
  expect_true(rec$W_interval[1L] <= rec$W && rec$W <= rec$W_interval[2L])
  expect_true(is.finite(rec$arithmetic$error_bound)); expect_null(rec$normalization)
  expect_true(rec$covariance$lambda_min_cert > 0)
  expect_false(rec$covariance$symmetrized)
})

test_that("a computed pretrend record binds to the cdid_MP W/Wpval fields and the summary reports the test", {
  skip_if_not_installed("did")
  z <- readRDS(pt_locate_fixture("cp02_result_adapter_N200.rds"))
  rep_ <- z$reporting; fr <- z$frame
  rec <- chained_pretest_record(rep_, 1e-6)
  units <- data.frame(id = fr$id, cohort = replace(fr$G, fr$G == 0, Inf), weight = fr$sample_weight, eligible = fr$eligible_ge2)
  params <- list(yname = "y", tname = "period", idname = "id", gname = "G", xformla = ~x, control_group = "nevertreated", anticipation = 0,
    alp = 0.05, tlist = z$calendar, glist = c(3, 4, 5), est_method = "ipw", base_period = "varying", bstrap = FALSE, cband = FALSE, biters = 1000)
  mp <- chained_did_MP(rep_, fr$id, quote(att_gt()), params, units = units, pretest = rec)
  expect_identical(mp$W, rec$W); expect_identical(mp$Wpval, rec$Wpval); expect_identical(mp$chained$pretest$status, "computed")
  txt <- paste(utils::capture.output(summary(mp)), collapse = "\n")
  expect_match(txt, paste0("Parallel-trends pre-test: W = ", format(rec$W, digits = 6), ", df = 1, p = ",
                           format(rec$Wpval, digits = 4), " (Wald test of the pretrend cell)."), fixed = TRUE)
  expect_no_match(txt, "Tested cells", fixed = TRUE)
  # an unavailable record is never bound as a statistic
  expect_error(chained_did_MP(rep_, fr$id, quote(att_gt()), params, units = units, pretest = list(status = "unavailable", W = 1, Wpval = 0.5,
    cells = "pre(g=3,t=2)", df = 1, description = "x")), class = "cdid_result_adapter_error")
  # post cells are never tested restrictions
  expect_error(chained_did_MP(rep_, fr$id, quote(att_gt()), params, units = units, pretest = list(W = 1, Wpval = 0.5,
    cells = "ATT(g=3,t=3)", df = 1, description = "x")), class = "cdid_result_adapter_error")
})
