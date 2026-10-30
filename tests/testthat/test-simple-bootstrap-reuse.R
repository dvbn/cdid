# Literal aggregation arithmetic with an isolated inference spy: no draws or fits.
simple_reuse_spy <- function() {
  ns <- asNamespace("cdid")
  e <- new.env(parent = ns)
  for (nm in c("aggregate_chained_ipw_targets", "ipw_pipeline_supported_inference")) {
    f <- get(nm, envir = ns, inherits = FALSE)
    environment(f) <- e
    assign(nm, f, envir = e)
  }
  e$calls <- list()
  e$compute_chained_ipw_inference <- function(
      estimates, if_full, cluster_full = NULL, alp = 0.05, cband = FALSE,
      biters = 1000L, btype = "rademacher", seed = NULL, multiplier_budget,
      inference_accuracy, enclosure = NULL, fixed_signs = NULL,
      eligible_full = NULL, bstrap = NULL, caller_rng = FALSE,
      endpoint = c("cells", "aggte_bins", "simple"), target_meta = NULL) {
    flags <- get("ipw_pipeline_resolve_bootstrap_flags", envir = ns)(
      bstrap, cband, endpoint = endpoint)
    e$calls[[length(e$calls) + 1L]] <- list(flags = flags, IF = if_full)
    n <- nrow(if_full)
    centered <- sweep(if_full, 2L, colMeans(if_full), "-")
    analytical <- sqrt(colSums(centered^2)) / n
    # Distinct fixed sentinels expose an accidental analytic-SE substitution.
    se <- analytical * if (flags$bootstrap) 1.5 else 1
    crit <- if (flags$simultaneous) 9 else 2
    V <- crossprod(centered) / n^2
    list(estimates = estimates, se = se, analytical_se = analytical,
      z_stat = estimates / se, p_value = rep(0.25, length(se)),
      ci_lower = estimates - 2 * se, ci_upper = estimates + 2 * se,
      cband_crit_val = crit, cband_lower = estimates - crit * se,
      cband_upper = estimates + crit * se, vcov = V, Omega = V * n,
      cluster_if = if_full, raw_cluster_sums = if_full,
      n_total = n, n_clusters = n, n_zero_eligible_clusters = 0L,
      clustered = FALSE,
      multiplier = if (flags$draws) list(metadata = list(n_draws = biters)) else NULL,
      rng = NULL, bootstrap = flags$bootstrap, simultaneous = flags$simultaneous,
      flag_law = flags$law,
      band_family = if (flags$simultaneous)
        "simultaneous_multiplier_sup_t_over_certified_ids" else
        "pointwise_normal_over_certified_ids",
      n_draws = if (flags$draws) biters else NA_integer_,
      numerical_certificate = NULL)
  }
  e
}

test_that("simple aggregation reuses one pointwise bootstrap and complete IF", {
  e <- simple_reuse_spy()
  meta <- data.frame(g = c(3, 3, 6, 6), time = c(3, 4, 6, 7),
    target_id = c("a", "b", "c", "d"), support_status = "available")
  X <- matrix(c(1, 2, 3, -2, -2, 1, 2, 3, 3, -1, -2, 1,
    -1, 4, 1, -4, 2, -2, -3, 4, -3, -4, -1, -2,
    1, -1, 2, -2, -1, 1, -2, 2), 8L, 4L, byrow = TRUE)
  frame_g <- c(3, 3, 6, 6, Inf, Inf, Inf, Inf)
  frame_w <- c(1, 3, 1, 3, 1, 3, 1, 3)
  mp <- structure(list(n = 8L, DIDparams = list(cdid = list()),
    chained = list(frame = data.frame(id = as.character(1:8)),
      target_id = meta$target_id, versions = list(),
      numerical_policy = list(decision = "D23"))), class = "cdid_MP")
  # Binary-exact reference: masses 1/4, denominator 1, estimate 3.
  # Derivatives are (1-3)+(2-3)=-3 and (4-3)+(5-3)=3.
  share <- (frame_w / 2) * 3 *
    (as.double(frame_g == 6) - as.double(frame_g == 3))
  expected_if <- rowMeans(X) + share
  flags <- list(FF = c(FALSE, FALSE), TF = c(TRUE, FALSE),
    TT = c(TRUE, TRUE), FT = c(FALSE, TRUE))
  for (flag in names(flags)) {
    fl <- flags[[flag]]
    e$calls <- list()
    a <- e$aggregate_chained_ipw_targets(c(1, 2, 4, 5), X, meta,
      frame_g, frame_w, type = "simple", bstrap = fl[1L], cband = fl[2L],
      multiplier_budget = list(), inference_accuracy = NULL)
    requested <- list(bstrap = fl[1L], cband = fl[2L], biters = 1000L,
      alp = 0.05, raw = list(bstrap = fl[1L], cband = fl[2L]),
      sources = c(bstrap = "explicit", cband = "explicit"))
    obj <- get("chained_did_AGGTEobj", envir = asNamespace("cdid"))(
      a, mp, quote(literal_simple_request()), requested = requested)
    inf <- obj$chained$inference
    expect_length(e$calls, 2L) # scalar + analytic fixed-share diagnostic
    expect_equal(sum(vapply(e$calls, function(z) z$flags$draws, logical(1L))),
      as.integer(fl[1L]), info = flag)
    expect_false(inf$cband, info = flag)
    expect_identical(inf$bstrap, fl[1L], info = flag)
    expect_identical(inf$requested$cband, fl[2L], info = flag)
    expect_identical(obj$DIDparams$cband, fl[2L], info = flag)
    expect_identical(obj$chained$requested, requested, info = flag)
    expect_equal(a$if_fixed_full[, 1L], rowMeans(X), tolerance = 0, info = flag)
    expect_equal(a$if_share_full[, 1L], share, tolerance = 0, info = flag)
    expect_equal(unname(obj$inf.function$simple.att), expected_if,
      tolerance = 0, info = flag)
    expect_equal(a$overall$if_full, expected_if, tolerance = 0, info = flag)
    expected_se <- sqrt(sum((expected_if - mean(expected_if))^2)) / 8 *
      if (fl[1L]) 1.5 else 1
    expect_equal(obj$overall.se, expected_se, tolerance = 0, info = flag)
    expect_equal(a$overall$se, expected_se, tolerance = 0, info = flag)
    expect_equal(unname(a$ci_lower), 3 - 2 * expected_se, tolerance = 0, info = flag)
    if (fl[1L]) {
      expect_identical(obj$chained$overall$draws_reused_from, "inference", info = flag)
      expect_identical(obj$chained$overall$inference_kind,
        "multiplier_bootstrap_iqr_se", info = flag)
      expect_identical(obj$chained$overall$n_draws, 1000L, info = flag)
    } else {
      expect_null(obj$chained$overall$draws_reused_from, info = flag)
      expect_identical(obj$chained$overall$inference_kind, "analytic", info = flag)
    }
  }
})
