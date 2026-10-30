# Independent checks of the development DGP. Expected values are computed by
# hand below (or directly from the potential-outcome definition), not with
# sim_truth() or any estimator.

hand_spec <- function() {
  sim_spec(periods = 1:4, cohorts = c(3, 4), p_x = 0.4,
           cohort_prob = rbind(c(0.5, 0.3, 0.2),    # P(G = 0, 3, 4 | X = 0)
                               c(0.2, 0.2, 0.6)),   # P(G = 0, 3, 4 | X = 1)
           theta = c(1, 2), delta = 0.5, kappa = 3,
           cohort_level = c(0, 0.5, -0.5), n_clusters = 5L)
}

test_that("population cohort weights and E[X | G] match hand calculation", {
  tr <- sim_truth(hand_spec())
  # joint P(X, G): X=0 row 0.6 * (0.5, 0.3, 0.2) = (0.30, 0.18, 0.12)
  #                X=1 row 0.4 * (0.2, 0.2, 0.6) = (0.08, 0.08, 0.24)
  expect_equal(tr$cohort$group, c(0, 3, 4))
  expect_equal(tr$cohort$pi, c(0.38, 0.26, 0.36))
  expect_equal(tr$cohort$mean_x, c(0.08 / 0.38, 0.08 / 0.26, 0.24 / 0.36))
})

test_that("ATT(g,t) truth matches hand calculation, zero before treatment", {
  tr <- sim_truth(hand_spec())$att_gt
  att <- function(g, t) tr$att[tr$group == g & tr$t == t]
  expect_equal(nrow(tr), 8L)
  for (t in 1:2) expect_identical(att(3, t), 0)
  for (t in 1:3) expect_identical(att(4, t), 0)
  expect_equal(att(3, 3), 1 + 3 * 4 / 13)          # theta_3 + kappa E[X|3]
  expect_equal(att(3, 4), 1 + 0.5 + 3 * 4 / 13)    # + delta * 1
  expect_equal(att(4, 4), 2 + 3 * 2 / 3)
})

test_that("aggregation truth matches hand calculation with unequal cohort weights", {
  tr <- sim_truth(hand_spec())
  a33 <- 25 / 13; a34 <- 63 / 26; a44 <- 4
  expect_equal(tr$simple, (0.26 * a33 + 0.26 * a34 + 0.36 * a44) / 0.88)
  th3 <- (a33 + a34) / 2
  expect_equal(tr$group$by$att, c(th3, a44))
  expect_equal(tr$group$overall, (0.26 * th3 + 0.36 * a44) / 0.62)
  e0 <- (0.26 * a33 + 0.36 * a44) / 0.62
  dyn <- tr$dynamic$by
  expect_equal(dyn$e, -3:1)
  expect_equal(dyn$att, c(0, 0, 0, e0, a34))
  expect_equal(tr$dynamic$overall, (e0 + a34) / 2)
  expect_equal(tr$calendar$by$t, c(3, 4))
  expect_equal(tr$calendar$by$att, c(a33, (0.26 * a34 + 0.36 * a44) / 0.62))
  expect_equal(tr$calendar$overall, mean(tr$calendar$by$att))
})

test_that("balanced panel: support, cohort timing and potential outcomes", {
  sp <- hand_spec()
  d <- sim_draw(sp, n = 400, panel = TRUE, seed = 7)
  expect_identical(names(d), c("id", "period", "G", "x", "z", "cluster", "y", "y0"))
  expect_identical(nrow(d), 1600L)
  expect_true(all(table(d$id) == 4L))
  expect_true(all(tapply(d$period, d$id, function(p) identical(as.numeric(p), c(1, 2, 3, 4)))))
  for (v in c("G", "x", "z", "cluster")) {
    expect_true(all(tapply(d[[v]], d$id, function(u) length(unique(u)) == 1L)))
  }
  expect_true(all(d$G %in% c(0, 3, 4)))
  expect_true(all(d$cluster %in% 1:5))
  # y - y0 is the unit-level effect from the stated potential-outcome equation
  treated <- d$G > 0 & d$period >= d$G
  theta <- c(`3` = 1, `4` = 2)
  effect <- ifelse(treated,
                   theta[as.character(d$G)] + 0.5 * (d$period - d$G) + 3 * d$x, 0)
  expect_equal(d$y - d$y0, unname(effect))
  expect_true(all(d$y[!treated] == d$y0[!treated]))
})

test_that("untreated outcome follows the stated equation up to unit and idiosyncratic noise", {
  sp <- hand_spec()
  sp$sd_unit <- 0
  sp$sd_eps <- 0
  d <- sim_draw(sp, n = 50, panel = TRUE, seed = 3)
  lev <- c(`0` = 0, `3` = 0.5, `4` = -0.5)
  expected <- lev[as.character(d$G)] + 0.25 * d$period +
    0.5 * d$x * (d$period - 1) + 1 * d$z
  expect_equal(d$y0, unname(expected))
})

test_that("repeated cross sections: one row per unit, n units per period", {
  d <- sim_draw(hand_spec(), n = 300, panel = FALSE, seed = 9)
  expect_identical(nrow(d), 1200L)
  expect_identical(anyDuplicated(d$id), 0L)
  expect_true(all(table(d$period) == 300L))
  treated <- d$G > 0 & d$period >= d$G
  effect <- ifelse(treated, c(`3` = 1, `4` = 2)[as.character(d$G)] +
                     0.5 * (d$period - d$G) + 3 * d$x, 0)
  expect_equal(d$y - d$y0, unname(effect))
})

test_that("large draw reproduces population cohort shares and cohort ATT", {
  sp <- hand_spec()
  n <- 100000
  d <- sim_draw(sp, n = n, panel = TRUE, seed = 2026)
  u <- d[d$period == 1, ]
  pi <- c(0.38, 0.26, 0.36)
  share <- as.numeric(table(factor(u$G, levels = c(0, 3, 4)))) / n
  expect_true(all(abs(share - pi) < 5 * sqrt(pi * (1 - pi) / n)))
  # sample average of Y(g) - Y(0) for cohort 4 at t = 4 versus 2 + 3 E[X | G = 4]
  c4 <- d[d$G == 4 & d$period == 4, ]
  diff <- c4$y - c4$y0
  expect_true(abs(mean(diff) - (2 + 3 * 2 / 3)) < 5 * sd(diff) / sqrt(nrow(c4)))
})

test_that("seed contract: reproducible draws and caller RNG state preserved", {
  sp <- hand_spec()
  set.seed(123)
  before <- .Random.seed
  kind <- RNGkind()
  a <- sim_draw(sp, 30, TRUE, seed = 5)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  b <- sim_draw(sp, 30, TRUE, seed = 5)
  expect_identical(a, b)
  expect_false(identical(a, sim_draw(sp, 30, TRUE, seed = 6)))
  expect_error(sim_draw(sp, 30, TRUE, seed = 1.5), "whole number")
  expect_error(sim_draw(sp, 30, TRUE), "seed")

  # no .Random.seed before the call: none afterwards either
  genv <- globalenv()
  saved <- get(".Random.seed", envir = genv)
  rm(".Random.seed", envir = genv)
  c1 <- sim_draw(sp, 30, TRUE, seed = 5)
  expect_false(exists(".Random.seed", envir = genv, inherits = FALSE))
  assign(".Random.seed", saved, envir = genv)
  expect_identical(c1, a)
})

test_that("specification errors are reported", {
  expect_error(sim_spec(cohort_prob = rbind(c(0.5, 0.2, 0.2, 0.2),
                                            c(0.1, 0.3, 0.3, 0.3))))
  expect_error(sim_spec(periods = 1:4, cohorts = c(1, 3), theta = c(1, 1),
                        cohort_prob = rbind(c(0.4, 0.3, 0.3), c(0.4, 0.3, 0.3)),
                        cohort_level = c(0, 0, 0)))
  expect_error(sim_spec(cohort_prob = rbind(c(0, 0.4, 0.3, 0.3),
                                            c(0, 0.4, 0.3, 0.3))),
               "positive")
})
