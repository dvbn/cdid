# Ordinary mode versus did::att_gt on panel data. Each case changes several
# options at once so that a small set of calls still exercises every
# forwarded argument; see helper-parity.R for the comparison.

d <- parity_panel()
u <- parity_unbalanced()

test_that("DR, covariates, never-treated, bootstrap with uniform band, both faster_mode", {
  for (fm in c(TRUE, FALSE)) {
    expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                                gname = "G", xformla = ~x, data = d,
                                est_method = "dr", bstrap = TRUE, cband = TRUE,
                                biters = 99, faster_mode = fm))
  }
})

test_that("IPW, not-yet-treated, universal base period, analytical SE, both faster_mode", {
  for (fm in c(TRUE, FALSE)) {
    expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                                gname = "G", xformla = ~x, data = d,
                                est_method = "ipw",
                                control_group = "notyettreated",
                                base_period = "universal", bstrap = FALSE,
                                faster_mode = fm))
  }
})

test_that("reg, anticipation, sampling weights, unit and cluster bootstrap", {
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                              gname = "G", xformla = ~x, data = d,
                              est_method = "reg", anticipation = 1,
                              weightsname = "w", clustervars = c("id", "cluster"),
                              bstrap = TRUE, cband = TRUE, biters = 99))
})

test_that("analytical cluster-robust SE, alp and fix_weights, both faster_mode", {
  for (fm in c(TRUE, FALSE)) {
    expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                                gname = "G", xformla = ~x + z, data = d,
                                clustervars = "cluster", bstrap = FALSE,
                                alp = 0.1, weightsname = "w",
                                fix_weights = "base_period", faster_mode = fm))
  }
})

test_that("no covariates, bootstrap without uniform band", {
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                              gname = "G", data = d, bstrap = TRUE,
                              cband = FALSE, biters = 99))
})

test_that("compute_inffunc = FALSE: point estimates only, both faster_mode", {
  for (fm in c(TRUE, FALSE)) {
    p <- expect_att_gt_parity(att_gt(yname = "y", tname = "period",
                                     idname = "id", gname = "G", xformla = ~x,
                                     data = d, compute_inffunc = FALSE,
                                     faster_mode = fm))
    full <- did::att_gt(yname = "y", tname = "period", idname = "id",
                        gname = "G", xformla = ~x, data = d, bstrap = FALSE,
                        faster_mode = fm)
    expect_equal(p$ours$att, full$att)
    e_ours <- tryCatch(cdid::aggte(p$ours), error = function(e) e)
    e_ref <- tryCatch(did::aggte(p$ref), error = function(e) e)
    expect_s3_class(e_ref, "error")
    expect_identical(conditionMessage(e_ours), conditionMessage(e_ref))
  }
})

test_that("unbalanced panel with allow_unbalanced_panel = TRUE, both faster_mode", {
  for (fm in c(TRUE, FALSE)) {
    expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                                gname = "G", xformla = ~x, data = u,
                                allow_unbalanced_panel = TRUE,
                                clustervars = "cluster", bstrap = FALSE,
                                faster_mode = fm))
  }
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                              gname = "G", data = u,
                              allow_unbalanced_panel = TRUE,
                              est_method = custom_rc, shift = 1, bstrap = FALSE))
})

test_that("unbalanced panel with allow_unbalanced_panel = FALSE balances like did", {
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                              gname = "G", xformla = ~x, data = u,
                              allow_unbalanced_panel = FALSE, bstrap = FALSE))
})

test_that("did's own input errors are unchanged", {
  expect_same_error(att_gt(yname = "y", tname = "period", idname = "id",
                           gname = "G", data = d, est_method = "nope"))
  expect_same_error(att_gt(yname = "y", tname = "period", gname = "G",
                           data = d))
  expect_same_error(att_gt(yname = "y", tname = "period", idname = "id",
                           gname = "G", data = d, control_group = "everyone"))
  expect_same_error(att_gt(yname = "y", tname = "period", idname = "id",
                           gname = "G", data = d,
                           clustervars = c("cluster", "x")))
})
