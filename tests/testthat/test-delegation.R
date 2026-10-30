# How ordinary mode forwards arguments: once, lazily, in the caller's
# environment, with missing arguments left missing.

test_that("each argument expression is evaluated exactly once", {
  d <- parity_panel()
  n_data <- 0
  n_panel <- 0
  res <- cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                      data = { n_data <- n_data + 1; d },
                      panel = { n_panel <- n_panel + 1; TRUE },
                      bstrap = FALSE, mode = "ordinary")
  expect_s3_class(res, "MP")
  expect_identical(n_data, 1)
  expect_identical(n_panel, 1)

  # mode = NULL evaluates `panel` to choose the mode, still only once
  r <- parity_rc()
  n_panel <- 0
  res <- cdid::att_gt(yname = "y", tname = "period", gname = "G", data = r,
                      panel = { n_panel <- n_panel + 1; FALSE }, bstrap = FALSE)
  expect_s3_class(res, "MP")
  expect_identical(n_panel, 1)
})

test_that("explicit mode = 'ordinary' does not evaluate panel before did does", {
  d <- parity_panel()
  order <- character()
  res <- cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                      data = d, bstrap = FALSE,
                      control_group = { order <- c(order, "control_group"); "nevertreated" },
                      panel = { order <- c(order, "panel"); TRUE },
                      mode = "ordinary")
  ref_order <- character()
  ref <- did::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                     data = d, bstrap = FALSE,
                     control_group = { ref_order <- c(ref_order, "control_group"); "nevertreated" },
                     panel = { ref_order <- c(ref_order, "panel"); TRUE })
  expect_identical(order, ref_order)
})

test_that("arguments are found in the caller's environment through wrappers", {
  d <- parity_panel()
  # A wrapper whose own environment holds a decoy `ctrl`: forwarding must use
  # the caller's `ctrl`, never re-evaluate the expression in the wrapper.
  wenv <- new.env(parent = baseenv())
  wenv$ctrl <- "decoy-control-group"
  my_fit <- function(dat, outcome, ...) {
    cdid::att_gt(yname = outcome, tname = "period", idname = "id", gname = "G",
                 data = dat, mode = "ordinary", ...)
  }
  environment(my_fit) <- wenv
  fx <- ~x  # one formula object (did stores it, with its environment)
  caller <- function(dat, form) {
    ctrl <- "notyettreated"
    my_fit(dat, "y", control_group = ctrl, xformla = form, bstrap = FALSE)
  }
  ours <- with_conditions(caller(d, fx))
  ref <- with_conditions(did::att_gt(yname = "y", tname = "period", idname = "id",
                                     gname = "G", xformla = fx, data = d,
                                     control_group = "notyettreated",
                                     bstrap = FALSE))
  expect_same_result(ours, ref, mp_fields)
  expect_identical(ours$value$DIDparams$control_group, "notyettreated")
})

test_that("do.call() forwarding matches did", {
  d <- parity_panel()
  args <- list(yname = "y", tname = "period", idname = "id", gname = "G",
               xformla = ~x, data = d, est_method = "ipw", bstrap = TRUE,
               biters = 99)
  set.seed(5)
  ours <- with_conditions(do.call(cdid::att_gt, c(args, mode = "ordinary")))
  set.seed(5)
  ref <- with_conditions(do.call(did::att_gt, args))
  expect_same_result(ours, ref, mp_fields)

  # a custom estimator and its extra argument, both through do.call()
  args2 <- c(args[c("yname", "tname", "idname", "gname", "data")],
             list(est_method = custom_panel, shift = 2, bstrap = FALSE))
  ours2 <- with_conditions(do.call(cdid::att_gt, c(args2, mode = "ordinary")))
  ref2 <- with_conditions(do.call(did::att_gt, args2))
  expect_same_result(ours2, ref2, mp_fields)
})

test_that("custom est_method receives arguments passed through ...", {
  d <- parity_panel()
  p0 <- expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                                    gname = "G", data = d, bstrap = FALSE,
                                    est_method = custom_panel))
  p2 <- expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                                    gname = "G", data = d, bstrap = FALSE,
                                    est_method = custom_panel, shift = 2))
  expect_equal(p2$ours$att, p0$ours$att + 2)
  expect_equal(p2$ours$se, p0$ours$se)

  # an extra argument with a built-in estimator gets did's own warning
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                              gname = "G", data = d, bstrap = FALSE,
                              est_method = "reg", shift = 2))
})

test_that("missing arguments stay missing; explicit NULL sentinel means FALSE", {
  u <- parity_unbalanced()
  # omitted: did's own default
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                              gname = "G", data = u, bstrap = FALSE))
  # explicit NULL: identical to did with FALSE (including did's warning)
  set.seed(3)
  ours <- with_conditions(cdid::att_gt(yname = "y", tname = "period",
                                       idname = "id", gname = "G", data = u,
                                       bstrap = FALSE, mode = "ordinary",
                                       allow_unbalanced_panel = NULL))
  set.seed(3)
  ref <- with_conditions(did::att_gt(yname = "y", tname = "period",
                                     idname = "id", gname = "G", data = u,
                                     bstrap = FALSE,
                                     allow_unbalanced_panel = FALSE))
  expect_same_result(ours, ref, mp_fields)
  expect_true(length(ref$warnings) > 0L)
})

test_that("did records the delegated call; that is the only wrapper text", {
  d <- parity_panel()
  res <- cdid::att_gt(yname = "y", tname = "period", idname = "id",
                      gname = "G", data = d, bstrap = FALSE, mode = "ordinary")
  cl <- res$DIDparams$call
  expect_identical(cl[[1L]], quote(did::att_gt))
  expect_identical(names(as.list(cl))[-1L],
                   c("yname", "tname", "idname", "gname", "data", "bstrap"))
  agg <- cdid::aggte(res, type = "simple")
  expect_identical(agg$call[[1L]], quote(did::aggte))
})

test_that("formula normalization is limited to did's generated intercept", {
  d <- parity_panel()
  generated <- cdid::att_gt(yname = "y", tname = "period", idname = "id",
    gname = "G", data = d, bstrap = FALSE, mode = "ordinary")
  expect_true(did_default_formula(generated))
  normalized <- normalize_did(generated)
  expect_identical(model.matrix(generated$DIDparams$xformla, d),
                   model.matrix(normalized$DIDparams$xformla, d))
  for (fx in list(local({ marker <- 19; ~1 }), local({ marker <- 23; ~x }))) {
    supplied <- cdid::att_gt(yname = "y", tname = "period", idname = "id",
      gname = "G", xformla = fx, data = d, bstrap = FALSE, mode = "ordinary")
    expect_false(did_default_formula(supplied))
    expect_identical(environment(normalize_did(supplied)$DIDparams$xformla),
                     environment(fx))
  }
})
