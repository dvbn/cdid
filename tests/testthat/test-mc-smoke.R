# The harness driving ordinary cdid on the development DGP. This checks the
# plumbing only: five replications are a smoke run, not statistical evidence.

test_that("harness runs ordinary att_gt end to end and labels the run", {
  sp <- parity_spec()
  tr <- sim_truth(sp)
  cells <- tr$att_gt[tr$att_gt$t >= 2, ]
  truth <- c(stats::setNames(cells$att, paste0("ATT(", cells$group, ",", cells$t, ")")),
             simple = tr$simple)
  gen <- function(seed) sim_draw(sp, n = 300, panel = TRUE, seed = seed)
  est <- function(data) {
    mp <- cdid::att_gt(yname = "y", tname = "period", idname = "id",
                       gname = "G", xformla = ~x, data = data, bstrap = FALSE,
                       mode = "ordinary")
    agg <- cdid::aggte(mp, type = "simple", bstrap = FALSE)
    list(estimate = c(stats::setNames(mp$att, paste0("ATT(", mp$group, ",", mp$t, ")")),
                      simple = agg$overall.att),
         mp = mp, agg = agg)
  }
  inf <- function(fit) {
    list(se = c(stats::setNames(fit$mp$se, paste0("ATT(", fit$mp$group, ",", fit$mp$t, ")")),
                simple = fit$agg$overall.se))
  }
  rec <- mc_run(mc_seeds(5, first = 1000), gen, est, inf, truth)
  expect_identical(nrow(rec), 5L * length(truth))
  # This prespecified moderate-size valid fixture must actually estimate
  # every target. Retention under deliberate failure is checked separately.
  expect_true(all(rec$status == "complete"))
  s <- mc_summary(rec)
  expect_identical(s$target, names(truth))
  expect_true(all(s$n_attempted == 5L))
  expect_match(attr(s, "run_type"), "not a D7 certification")
})
