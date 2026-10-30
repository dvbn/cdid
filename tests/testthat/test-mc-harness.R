# Hand-constructed replications with known outcomes. Seeds select the
# scripted behaviour, so each count below can be verified by reading the table.
#
# seed  estimate  se    stage outcome                     covered (truth 1, 95%)
# 101   1.0       0.1   complete                          yes
# 102   1.5       0.1   complete                          no
# 103   -         -     estimate() error                  -
# 104   1.1       -     infer() error                     -
# 105   1.2       NA    invalid SE                        -
# 106   0.8       0     invalid SE (zero is kept, not used) -
# 107   0.9       0.2   complete, with a warning          yes
# 108   Inf       0.1   invalid estimate                  -
# 109   -         -     generate() error                  -
scripted <- list(
  `101` = list(est = 1.0, se = 0.1), `102` = list(est = 1.5, se = 0.1),
  `103` = list(est_error = TRUE), `104` = list(est = 1.1, inf_error = TRUE),
  `105` = list(est = 1.2, se = NA), `106` = list(est = 0.8, se = 0),
  `107` = list(est = 0.9, se = 0.2, warn = TRUE),
  `108` = list(est = Inf, se = 0.1), `109` = list(gen_error = TRUE)
)
scripted_generate <- function(seed) {
  s <- scripted[[as.character(seed)]]
  if (isTRUE(s$gen_error)) stop("generation broke")
  s
}
scripted_estimate <- function(s) {
  if (isTRUE(s$est_error)) stop("estimator broke")
  if (isTRUE(s$warn)) warning("fitted with a warning")
  list(estimate = c(att = s$est), s = s)
}
scripted_infer <- function(fit) {
  if (isTRUE(fit$s$inf_error)) stop("variance broke")
  list(se = c(att = fit$s$se))
}

run_scripted <- function() {
  mc_run(seeds = 101:109, generate = scripted_generate,
         estimate = scripted_estimate, infer = scripted_infer,
         truth = c(att = 1))
}

test_that("every attempted replication is retained with its seed and stage", {
  rec <- run_scripted()
  expect_identical(nrow(rec), 9L)
  expect_identical(rec$seed, 101:109)
  expect_identical(rec$status, c(
    "complete", "complete", "estimation_failed", "inference_failed",
    "inference_failed", "inference_failed", "complete", "estimation_failed",
    "generation_failed"))
  expect_identical(rec$message[3], "estimator broke")
  expect_identical(rec$message[4], "variance broke")
  expect_identical(rec$message[9], "generation broke")
  expect_identical(rec$covered, c(TRUE, FALSE, NA, NA, NA, NA, TRUE, NA, NA))
  # invalid SEs are recorded as returned, never replaced
  expect_identical(rec$se[5], NA_real_)
  expect_identical(rec$se[6], 0)
  expect_identical(rec$estimate[4], 1.1)
  expect_identical(rec$n_warnings, c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L))
  expect_identical(rec$warnings[7], "estimate: fitted with a warning")
  z <- qnorm(0.975)
  expect_equal(rec$lower[1:2], c(1.0, 1.5) - z * 0.1)
  expect_equal(rec$upper[7], 0.9 + z * 0.2)
})

test_that("summary keeps failures in unconditional denominators", {
  s <- mc_summary(run_scripted())
  expect_identical(s$n_attempted, 9L)
  expect_identical(s$n_generation_failed, 1L)
  expect_identical(s$n_estimation_failed, 2L)
  expect_identical(s$n_inference_failed, 3L)
  expect_identical(s$n_complete, 3L)
  expect_identical(s$n_with_warnings, 1L)
  expect_identical(s$n_covered, 2L)
  expect_equal(s$coverage_unconditional, 2 / 9)
  expect_equal(s$coverage_conditional, 2 / 3)
  expect_equal(s$mcse_coverage_unconditional, sqrt((2 / 9) * (7 / 9) / 9))
  # conditional on a valid estimate: seeds 101, 102, 104, 105, 106, 107
  est <- c(1.0, 1.5, 1.1, 1.2, 0.8, 0.9)
  expect_identical(s$n_estimate, 6L)
  expect_equal(s$bias, mean(est) - 1)
  expect_equal(s$emp_sd, sd(est))
  expect_equal(s$rmse, sqrt(mean((est - 1)^2)))
  # complete replications only: seeds 101, 102, 107
  expect_equal(s$mean_se, mean(c(0.1, 0.1, 0.2)))
  expect_equal(s$se_ratio, mean(c(0.1, 0.1, 0.2)) / sd(c(1.0, 1.5, 0.9)))
  expect_match(attr(s, "run_type"), "not a D7 certification")
})

test_that("targets fail independently within a replication", {
  rec <- mc_run(seeds = 1:2,
                generate = function(seed) seed,
                estimate = function(d) list(estimate = c(a = 1, b = if (d == 1) NA else 2)),
                infer = function(fit) list(se = c(a = 0.5, b = 0.5), crit = 2),
                truth = c(a = 0, b = 2))
  expect_identical(rec$target, c("a", "b", "a", "b"))
  expect_identical(rec$status, c("complete", "estimation_failed",
                                 "complete", "complete"))
  expect_identical(rec$covered, c(TRUE, NA, TRUE, TRUE))
  s <- mc_summary(rec)
  expect_identical(s$n_attempted, c(2L, 2L))
  expect_equal(s$coverage_unconditional, c(1, 0.5))
  expect_equal(s$coverage_conditional, c(1, 1))
})

test_that("replications replay deterministically and keep caller RNG state", {
  sp <- sim_spec()
  truth <- c(diff = 0)
  gen <- function(seed) sim_draw(sp, n = 60, panel = TRUE, seed = seed)
  # an estimator that itself draws random numbers, like a bootstrap
  est <- function(d) list(estimate = c(diff = mean(d$y[d$period == 2]) -
                                         mean(d$y[d$period == 1]) + rnorm(1, 0, 1e-3)))
  inf <- function(fit) list(se = c(diff = 0.1 + runif(1, 0, 1e-3)))
  set.seed(99)
  before <- .Random.seed
  seeds <- mc_seeds(4, first = 500)
  expect_identical(seeds, 500:503)
  r1 <- mc_run(seeds, gen, est, inf, truth)
  r2 <- mc_run(seeds, gen, est, inf, truth)
  expect_identical(.Random.seed, before)
  drop_time <- function(r) {
    r <- r[, setdiff(names(r), "elapsed")]
    rownames(r) <- NULL
    r
  }
  expect_identical(drop_time(r1), drop_time(r2))
  r3 <- mc_replay(r1, seeds = 502, gen, est, inf)
  expect_identical(drop_time(r3), drop_time(r1[r1$seed == 502, ]))
  expect_identical(r3$rep, 3L)
})

test_that("data and randomized estimation use distinct reproducible streams", {
  gen <- function(seed) with_seed(seed, runif(1L))
  est <- function(d) list(estimate = c(a = runif(1L) - d))
  inf <- function(fit) list(se = c(a = 1))
  out <- mc_run(101:110, gen, est, inf, c(a = 0))
  # The original nested with_seed implementation reused the exact draw in
  # every replication. An independent next-stream reference fixes its value.
  reference <- vapply(101:110, function(seed) with_seed(seed, {
    next_state <- parallel::nextRNGStream(.Random.seed)
    data_draw <- runif(1L)
    assign(".Random.seed", next_state, globalenv())
    runif(1L) - data_draw
  }), numeric(1))
  expect_equal(out$estimate, reference, tolerance = 0)
  expect_true(all(out$estimate != 0))
})

test_that("malformed returns and invalid critical values remain failed attempts", {
  gen <- function(seed) seed
  est <- function(d) list(estimate = c(a = 0), seed = d)
  inf <- function(fit) {
    if (fit$seed == 1) return(1)
    if (fit$seed == 2) return(list(se = c(a = 1), crit = -1))
    if (fit$seed == 3) return(list(se = c(a = 1), crit = 0))
    if (fit$seed == 4) return(list(se = c(a = "bad")))
    list(se = c(a = 1), crit = 2)
  }
  out <- mc_run(1:5, gen, est, inf, c(a = 0))
  expect_identical(out$status, c(rep("inference_failed",4), "complete"))
  expect_identical(out$critical_value[2:3], c(-1,0))
  expect_identical(out$covered, c(NA,NA,NA,NA,TRUE))
  expect_equal(mc_summary(out)$coverage_unconditional, 1/5)
  broken <- mc_run(1:2, gen, function(d) 1, inf, c(a=0))
  expect_true(all(broken$status == "estimation_failed"))
  expect_error(mc_run(1:2,gen,est,inf,c(a=NA_real_)))
  expect_error(mc_run(c(1,NA),gen,est,inf,c(a=0)))
  expect_error(mc_seeds(2, first=.Machine$integer.max))
})
