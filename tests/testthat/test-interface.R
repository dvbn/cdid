# Public interface tests (installed host with did 2.5.1).
#
# CP02 candidate (October 16, 2026): the chained branch of att_gt() now runs
# for explicit est_method = "ipw" or an original-IPW callback. The CP01 expectation that EVERY chained
# request stops with "cdid_chained_unavailable" is superseded for these IPW requests;
# the retained evidence is this file's history (input pin
# 0b3564... in input_manifest.json of the CP02 public-binding run). The final
# "dr" default still refuses with the familiar
# "cdid_chained_unavailable" class (now also "cdid_chained_request_est_method")
# before the data promise is forced.

test_that("att_gt keeps did's formals, order and defaults, plus the extensions", {
  ours <- as.list(formals(cdid::att_gt))
  ref <- as.list(formals(did::att_gt))
  n_up <- length(ref)
  expect_identical(names(ours)[seq_len(n_up)], names(ref))
  expect_identical(names(ref)[n_up], "...")
  same <- setdiff(names(ref), "allow_unbalanced_panel")
  expect_identical(ours[same], ref[same])
  expect_identical(ref[["allow_unbalanced_panel"]], FALSE)
  expect_true("allow_unbalanced_panel" %in% names(ours))
  expect_null(ours[["allow_unbalanced_panel"]])
  expect_identical(ours[["est_method"]], "dr")
  expect_identical(ours[["base_period"]], "varying")
  expect_identical(
    ours[-seq_len(n_up)],
    list(mode = NULL, frame = NULL, selection = "none",
         selection_model = "pair_logit", selection_pooling = "separate",
         selection_xformla = NULL, selection_universe = "F", selection_rho = NULL, selection_justification = NULL,
         selection_est_method = "ipw", imputation_pooling = "separate",
         gmm_weighting = quote(c("optimal", "identity")),
         gmm_covariance = quote("tikhonov"),
         pt_assumption = quote(c("post", "all")),
         pt_start = NULL,
         control = NULL)
  )
  # the forwarded formal list tracks the loaded did
  expect_identical(cdid:::upstream_att_gt_formals, setdiff(names(ref), "..."))
})

test_that("aggte keeps did's signature", {
  expect_identical(formals(cdid::aggte), formals(did::aggte))
  expect_identical(cdid:::upstream_aggte_formals, names(formals(did::aggte)))
})

test_that("workflow functions are did's own objects", {
  expect_identical(cdid::ggdid, did::ggdid)
  expect_identical(cdid::tidy, did::tidy)
  expect_identical(cdid::glance, did::glance)
  # conditional_did_pretest is cdid's own wrapper (O14/D66): did's formals first, in order.
  up <- names(formals(did::conditional_did_pretest))
  expect_identical(names(formals(cdid::conditional_did_pretest))[seq_along(up)], up)
  expect_identical(sort(getNamespaceExports("cdid")),
                   sort(c("att_gt", "aggte", "ggdid", "tidy", "glance",
                          "conditional_did_pretest")))
})

test_that("panel default resolves to chained DR; unsupported est_method strings are refused before touching data", {
  d <- parity_panel()
  e <- expect_error(
    cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                 data = stop("data must not be evaluated"), est_method = "stdipw"),
    class = "cdid_chained_unavailable"
  )
  expect_s3_class(e, "cdid_chained_request_est_method")
  expect_match(conditionMessage(e), "not implemented")
  expect_match(conditionMessage(e), "est_method = \"dr\"", fixed = TRUE)
  expect_match(conditionMessage(e), "est_method = \"ipw\"", fixed = TRUE)
  expect_match(conditionMessage(e), "mode = \"ordinary\"", fixed = TRUE)
  # Unsupported string chained requests are refused; IPW is never substituted
  for (m in list("stdipw", "dr_imp")) {
    expect_error(
      cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                   data = d, est_method = m),
      class = "cdid_chained_unavailable"
    )
  }
  # D48 CP07: default est_method = "dr" and explicit "dr", "ipw", "reg", or custom IPW callback pass the gate and reach the data
  expect_error(
    cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                 data = stop("data must not be evaluated")),
    "data must not be evaluated"
  )
  for (m in list("dr", "ipw", "reg", custom_panel)) {
    expect_error(
      cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                   data = stop("data must not be evaluated"), est_method = m),
      "data must not be evaluated"
    )
  }
  mp_dr <- cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                        data = d, mode = "chained", bstrap = FALSE, cband = FALSE)
  expect_s3_class(mp_dr, c("cdid_MP", "MP"))
  expect_identical(mp_dr$chained$estimator, "chained doubly robust GMM")
})

test_that("explicit allow_unbalanced_panel = FALSE in chained mode points to ordinary mode", {
  d <- parity_panel()
  for (m in list(NULL, "chained")) {
    e <- expect_error(
      cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                   data = d, allow_unbalanced_panel = FALSE, mode = m),
      class = "cdid_chained_balanced"
    )
    expect_match(conditionMessage(e), "mode = \"ordinary\"", fixed = TRUE)
  }
  # TRUE and NULL are not this error: with an unsupported est_method they reach the est_method refusal,
  # and with default est_method = "dr" they pass the gate and reach data.
  for (aup in list(TRUE, NULL)) {
    expect_error(
      cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                   data = stop("data must not be evaluated"),
                   allow_unbalanced_panel = aup, est_method = "stdipw"),
      class = "cdid_chained_unavailable"
    )
    expect_error(
      cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                   data = stop("data must not be evaluated"),
                   allow_unbalanced_panel = aup),
      "data must not be evaluated"
    )
  }
})

test_that("chained mode with repeated cross sections and invalid modes are refused", {
  d <- parity_rc()
  expect_error(
    cdid::att_gt(yname = "y", tname = "period", gname = "G", data = d,
                 panel = FALSE, mode = "chained"),
    class = "cdid_chained_panel"
  )
  for (m in list("Ordinary", "", NA_character_, c("ordinary", "chained"), 1)) {
    expect_error(
      cdid::att_gt(yname = "y", tname = "period", gname = "G", data = d,
                   panel = FALSE, mode = m),
      class = "cdid_mode_error"
    )
  }
  expect_error(
    cdid::att_gt(yname = "y", tname = "period", gname = "G", data = d,
                 panel = NA),
    class = "cdid_mode_error"
  )
})

test_that("ordinary mode delegates without validating chained-only arguments", {
  d <- parity_panel()
  bad <- list(frame = list(frame = d), selection = list(selection = "ipw"),
              selection_model = list(selection_model = quote(stop("selection model promise was forced"))),
              selection_pooling = list(selection_pooling = quote(stop("selection pooling promise was forced"))),
              selection_xformla = list(selection_xformla = ~x),
              selection_universe = list(selection_universe = "invalid"),
              selection_rho = list(selection_rho = list(NA)),
              selection_justification = list(selection_justification = 123),
              gmm_weighting = list(gmm_weighting = "identity"),
              gmm_covariance = list(gmm_covariance = "invalid_handler"),
              pt_assumption = list(pt_assumption = "all"),
              pt_start = list(pt_start = list(NA_real_, Inf)),
              control = list(control = list(bogus = 1)))
  base <- list(yname = "y", tname = "period", idname = "id", gname = "G",
               data = d, bstrap = FALSE, mode = "ordinary")
  for (nm in names(bad)) {
    ours <- with_conditions(do.call(cdid::att_gt, c(base, bad[[nm]])))
    ref <- with_conditions(did::att_gt(yname = "y", tname = "period",
                                       idname = "id", gname = "G", data = d,
                                       bstrap = FALSE))
    expect_same_result(ours, ref, mp_fields)
  }
})

test_that("ordinary dispatch never forces chained-only promises", {
  d <- parity_panel()
  ours <- with_conditions(cdid::att_gt(yname = "y", tname = "period",
    idname = "id", gname = "G", data = d, bstrap = FALSE, mode = "ordinary",
    frame = stop("forced frame"), selection = stop("forced selection"),
    selection_xformla = stop("forced selection formula"),
    selection_universe = stop("forced fitting universe"), selection_rho = stop("forced rho"),
    selection_justification = stop("forced transport justification"),
    gmm_weighting = stop("forced weighting"),
    gmm_covariance = stop("forced gmm_covariance"),
    pt_assumption = stop("forced PT"), pt_start = stop("forced pt_start"),
    control = stop("forced control")))
  ref <- with_conditions(did::att_gt(yname = "y", tname = "period",
    idname = "id", gname = "G", data = d, bstrap = FALSE))
  expect_same_result(ours, ref, mp_fields)
})

test_that("neutral extension values are accepted in ordinary mode", {
  d <- parity_panel()
  expect_att_gt_parity(att_gt(yname = "y", tname = "period", idname = "id",
                              gname = "G", data = d, bstrap = FALSE),
                       mode = "ordinary")
  set.seed(1)
  a <- cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                    data = d, bstrap = FALSE, mode = "ordinary", frame = NULL,
                    selection = "none", selection_xformla = NULL,
                    gmm_weighting = "optimal", gmm_covariance = "tikhonov",
                    pt_assumption = "post",
                    control = NULL)
  set.seed(1)
  b <- did::att_gt(yname = "y", tname = "period", idname = "id", gname = "G",
                   data = d, bstrap = FALSE)
  expect_equal(normalize_did(a), normalize_did(b), tolerance = 0)
})
