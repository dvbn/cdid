# aggte() and output methods on ordinary results.

d <- parity_panel()
pair_boot <- local({
  set.seed(11)
  ours <- cdid::att_gt(yname = "y", tname = "period", idname = "id",
                       gname = "G", xformla = ~x, data = d, biters = 99,
                       clustervars = "cluster", mode = "ordinary")
  set.seed(11)
  ref <- did::att_gt(yname = "y", tname = "period", idname = "id",
                     gname = "G", xformla = ~x, data = d, biters = 99,
                     clustervars = "cluster")
  list(ours = ours, ref = ref)
})
pair_analytic <- local({
  ours <- cdid::att_gt(yname = "y", tname = "period", idname = "id",
                       gname = "G", xformla = ~x, data = d, bstrap = FALSE,
                       control_group = "notyettreated", mode = "ordinary")
  ref <- did::att_gt(yname = "y", tname = "period", idname = "id",
                     gname = "G", xformla = ~x, data = d, bstrap = FALSE,
                     control_group = "notyettreated")
  list(ours = ours, ref = ref)
})

test_that("all aggregation types match did, bootstrap and analytical", {
  for (type in c("simple", "group", "dynamic", "calendar")) {
    expect_aggte_parity(pair_boot, type = type)
    expect_aggte_parity(pair_analytic, type = type)
  }
})

test_that("aggregation options are forwarded", {
  expect_aggte_parity(pair_boot, type = "dynamic", min_e = -2, max_e = 1,
                      balance_e = 1)
  expect_aggte_parity(pair_boot, type = "dynamic", na.rm = TRUE)
  expect_aggte_parity(pair_boot, type = "group", bstrap = FALSE)
  expect_aggte_parity(pair_boot, type = "calendar", cband = FALSE, alp = 0.1,
                      biters = 49)
  expect_aggte_parity(pair_boot, type = "simple", clustervars = "cluster")
  # an override did cannot honour: same warning, same result
  expect_aggte_parity(pair_analytic, type = "dynamic", clustervars = "cluster")
})

test_that("RC aggregation matches did", {
  r <- parity_rc()
  set.seed(12)
  ours <- cdid::att_gt(yname = "y", tname = "period", gname = "G",
                       xformla = ~x, data = r, panel = FALSE, biters = 99)
  set.seed(12)
  ref <- did::att_gt(yname = "y", tname = "period", gname = "G",
                     xformla = ~x, data = r, panel = FALSE, biters = 99)
  expect_aggte_parity(list(ours = ours, ref = ref), type = "dynamic")
})

test_that("did's output methods give the same results on delegated objects", {
  ours <- pair_boot$ours
  ref <- pair_boot$ref
  expect_identical(cdid::tidy(ours), did::tidy(ref))
  expect_identical(cdid::glance(ours), did::glance(ref))
  expect_identical(nobs(ours), nobs(ref))

  # printed summaries differ only in the recorded call
  ours_call_text <- ours
  ours_call_text$DIDparams$call <- ref$DIDparams$call
  expect_identical(capture.output(summary(ours_call_text)),
                   capture.output(summary(ref)))
  expect_identical(capture.output(print(ours_call_text)),
                   capture.output(print(ref)))

  set.seed(13)
  a_ours <- cdid::aggte(ours, type = "dynamic")
  set.seed(13)
  a_ref <- did::aggte(ref, type = "dynamic")
  expect_identical(cdid::tidy(a_ours), did::tidy(a_ref))
  expect_identical(cdid::glance(a_ours), did::glance(a_ref))
  a_ours$call <- a_ref$call
  a_ours$DIDparams$call <- a_ref$DIDparams$call
  expect_identical(capture.output(summary(a_ours)),
                   capture.output(summary(a_ref)))

  p_ours <- cdid::ggdid(a_ours)
  p_ref <- did::ggdid(a_ref)
  expect_identical(class(p_ours), class(p_ref))
  if (inherits(p_ref, "ggplot")) expect_equal(p_ours$data, p_ref$data)
  expect_identical(class(cdid::ggdid(ours)), class(did::ggdid(ref)))
})
