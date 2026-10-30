# Repeated cross sections: mode = NULL resolves to ordinary.

r <- parity_rc()

test_that("RC default mode: DR, bootstrap with uniform band, both faster_mode", {
  for (fm in c(TRUE, FALSE)) {
    expect_att_gt_parity(att_gt(yname = "y", tname = "period", gname = "G",
                                xformla = ~x, data = r, panel = FALSE,
                                est_method = "dr", bstrap = TRUE, biters = 99,
                                faster_mode = fm),
                         mode = NULL)
  }
})

test_that("RC explicit ordinary mode: IPW, not-yet-treated, analytical SE", {
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", gname = "G",
                              xformla = ~x, data = r, panel = FALSE,
                              est_method = "ipw",
                              control_group = "notyettreated", bstrap = FALSE,
                              faster_mode = FALSE))
})

test_that("RC: reg, universal base, weights, idname, cluster bootstrap", {
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                              gname = "G", xformla = ~x, data = r,
                              panel = FALSE, est_method = "reg",
                              base_period = "universal", weightsname = "w",
                              clustervars = "cluster", bstrap = TRUE,
                              biters = 99),
                       mode = NULL)
})

test_that("RC: anticipation, analytical cluster SE, custom estimator with dots", {
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", gname = "G",
                              data = r, panel = FALSE, anticipation = 1,
                              clustervars = "cluster", bstrap = FALSE),
                       mode = NULL)
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", gname = "G",
                              data = r, panel = FALSE, est_method = custom_rc,
                              shift = 0.5, bstrap = FALSE),
                       mode = NULL)
})

test_that("RC: compute_inffunc = FALSE", {
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", gname = "G",
                              xformla = ~x, data = r, panel = FALSE,
                              compute_inffunc = FALSE),
                       mode = NULL)
})
