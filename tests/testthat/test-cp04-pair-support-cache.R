# Expected results come from the unchanged no-context owner scanner. These
# witnesses cover CP04's shared constructor seam independently of fit method.
psc <- function(name) getFromNamespace(name, "cdid")
psc_control <- function(chunk = 5L, cap = 1e8, work = 1e12) psc("ipw_owner_control")(
  max_chunk_rows = chunk, max_retained_entries = cap, max_matrix_entries = cap,
  max_work_ops = work, max_output_entries = 1e8)
psc_case <- function(zero_cohort = FALSE, changed_weights = FALSE, char_ids = FALSE) {
  N <- 48L; id <- seq_len(N); g <- rep(c(3, 4, 5, Inf), length.out = N)
  units <- data.frame(id = if (char_ids) sprintf("unit-%03d", id) else id, g = g,
    w = 1 + (id %% 7L) / 10 + 1 / 3, x = (id %% 5L) / 4)
  if (changed_weights) units$w <- units$w + (id %% 3L) / 20
  k <- rep(id, each = 4L); tt <- rep(1:4, N)
  data <- data.frame(id = units$id[k], t = tt, g = g[k], y = 0.3 * k + tt + (k * tt %% 7L) / 10)
  data$y[k == 13L & tt == 4L] <- NA_real_
  drop <- (k == 2L & tt == 2L) | (k == 9L & tt == 3L)
  if (zero_cohort) drop <- drop | (g[k] == 5 & tt == 1L)
  data <- data[!drop, ]
  data <- data[order((seq_len(nrow(data)) * 37L) %% nrow(data)), ]
  ctrl <- psc_control()
  panel <- psc("construct_ipw_panel_owner")(data, "id", "t", "y", "g", ctrl)
  frame <- psc("construct_ipw_frame_owner")(panel, units$id,
    cbind(`(Intercept)` = 1, x = units$x), weights = units$w, control = ctrl)
  list(frame = frame, panel = panel, units = units, data = data)
}
psc_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.list(x)) { a <- attributes(x); x <- lapply(x, psc_strip); attributes(x) <- a }
  x
}
psc_spec <- function(g, K, s = 1, t = 2) list(g = g, K = K, s = s, t = t,
  functional_key = paste0("cmpcat:pair-support:", g, ":", paste(K, collapse = ","), ":", s, ":", t))
psc_plain <- function(frame, specs, control) lapply(seq_along(specs), function(j) {
  s <- specs[[j]]
  psc("construct_ipw_link_view")(frame, s$s, s$t, s$g, s$K, control = control,
    link_id = paste0("cp04:catalog:", j), functional_key = s$functional_key)
})
psc_context <- function(frame, control, n = 20) {
  ctx <- psc("ipw_owner_link_support_context")(frame, n, 0, control, max_cohorts = 4)
  psc("cp04_pair_support_context")(ctx, frame, n, 0, control, 4)
}
psc_view <- function(a, ctx, g = 3, K = Inf, s = 1, t = 2, control = psc_control(), rule = TRUE) {
  psc("construct_ipw_link_view")(a$frame, s, t, g, K, control = control,
    construction_context = ctx, require_support = rule)
}

testthat::test_that("incomplete raw pair supports reuse across g/K with bitwise original summaries", {
  specs <- list(psc_spec(3, Inf), psc_spec(4, Inf), psc_spec(3, c(4, Inf)),
    psc_spec(4, c(Inf, 3)), psc_spec(5, c(Inf, 3)), psc_spec(4, c(3, Inf)),
    psc_spec(3, Inf, 2, 3), psc_spec(5, c(Inf, 4), 2, 3))
  for (char_ids in c(FALSE, TRUE)) for (chunk in c(1L, 5L, 11L)) {
    a <- psc_case(char_ids = char_ids); ctrl <- psc_control(chunk)
    got <- psc("cp04_public_link_views")(a$frame, specs, list(descriptors = list()), ctrl,
      reuse_pair_support = TRUE)
    want <- psc_plain(a$frame, specs, ctrl)
    testthat::expect_true(identical(psc_strip(got$views), psc_strip(want), num.eq = FALSE))
    testthat::expect_equal(got$construction_work$pair_support$scans, 2)
    testthat::expect_equal(got$construction_work$pair_support$reuses, 6)
    testthat::expect_equal(got$construction_work$pair_support$records, 2)
    testthat::expect_true(all(vapply(got$views, function(v) identical(v$accounting$endpoints_read, 0L), logical(1))))
    # The default private helper preserves its historical path/diagnostics.
    old <- psc("cp04_public_link_views")(a$frame, specs, list(descriptors = list()), ctrl)
    testthat::expect_null(old$construction_work$pair_support)
    testthat::expect_true(identical(psc_strip(old$views), psc_strip(want), num.eq = FALSE))
  }
})

testthat::test_that("zero support is checked on the requested subset in original order", {
  a <- psc_case(zero_cohort = TRUE); ctrl <- psc_control()
  specs <- list(psc_spec(3, Inf), psc_spec(4, c(Inf, 3)))
  aux <- list(descriptors = list(psc_spec(5, Inf), psc_spec(4, Inf)))
  got <- psc("cp04_public_link_views")(a$frame, specs, aux, ctrl, reuse_pair_support = TRUE)
  testthat::expect_true(identical(psc_strip(got$views[1:2]), psc_strip(psc_plain(a$frame, specs, ctrl)), num.eq = FALSE))
  original_error <- tryCatch(psc_view(a, NULL, g = 5), error = identity)
  testthat::expect_s3_class(original_error, "cdid_ipw_pair_support_error")
  testthat::expect_identical(got$aux_errors[[1L]], original_error)
  testthat::expect_identical(got$aux_index, c(NA_integer_, 3L))
  testthat::expect_equal(got$construction_work$pair_support$scans, 1)
  testthat::expect_equal(got$construction_work$pair_support$reuses, 3)
  # Unchecked standalone semantics retain absent/zero cohorts and raw rates.
  ctx <- psc_context(a$frame, ctrl)
  for (K in list(Inf, c(5, Inf), c(99, Inf))) {
    got_v <- psc_view(a, ctx, K = K, rule = FALSE)
    ref_v <- psc_view(a, NULL, K = K, rule = FALSE)
    testthat::expect_true(identical(psc_strip(got_v), psc_strip(ref_v), num.eq = FALSE))
  }
  # A failed auxiliary support can populate raw counts without licensing its
  # failure for a different valid subset.
  ctx <- psc_context(a$frame, ctrl)
  testthat::expect_identical(tryCatch(psc_view(a, ctx, g = 5), error = identity), original_error)
  testthat::expect_identical(psc_strip(psc_view(a, ctx, g = 4)), psc_strip(psc_view(a, NULL, g = 4)))
  testthat::expect_equal(ctx$pair_support$scans, 1)
})

testthat::test_that("unadmitted pair reservation and exhausted capacity keep original scans and refusals", {
  a <- psc_case(); ctrl <- psc_control(); n <- 2
  old_ctx <- psc("ipw_owner_link_support_context")(a$frame, n, 0, ctrl, max_cohorts = 4)
  pre <- psc("preflight_ipw_link_view")(a$frame$n_a_units, a$frame$n_covariates, 4, ctrl)
  added <- (256 + 64 * length(a$frame$a_cohorts)) * (n + 1)
  cap_at <- old_ctx$link_entries + added + pre$concurrent_entries
  for (offset in c(-1, 0)) {
    tight <- psc_control(cap = cap_at + offset)
    ctx <- psc("ipw_owner_link_support_context")(a$frame, n, 0, tight, max_cohorts = 4)
    ctx <- psc("cp04_pair_support_context")(ctx, a$frame, n, 0, tight, 4)
    testthat::expect_identical(is.environment(ctx$pair_support), offset == 0)
    testthat::expect_true(identical(psc_strip(psc_view(a, ctx, control = tight)),
      psc_strip(psc_view(a, NULL, control = tight)), num.eq = FALSE))
  }
  ctx <- psc_context(a$frame, ctrl, n)
  for (pair in list(c(1, 2), c(2, 3), c(3, 4))) {
    testthat::expect_identical(psc_strip(psc_view(a, ctx, s = pair[1], t = pair[2])),
      psc_strip(psc_view(a, NULL, s = pair[1], t = pair[2])))
  }
  testthat::expect_equal(ctx$pair_support$n_records, n)
  testthat::expect_equal(ctx$pair_support$scans, n)
  low <- psc_control(work = 1)
  e0 <- tryCatch(psc_view(a, NULL, control = low), error = identity)
  e1 <- tryCatch(psc_view(a, ctx, control = low), error = identity)
  testthat::expect_s3_class(e0, "cdid_ipw_owner_bound_error")
  testthat::expect_identical(e1, e0)
})

testthat::test_that("pair records require intact seals, exact owner/weight identity and chunk step", {
  a <- psc_case(); ctrl <- psc_control(); ctx <- psc_context(a$frame, ctrl)
  psc_view(a, ctx)
  testthat::expect_equal(ctx$pair_support$scans, 1)
  # A different step gets its own original-order raw scan.
  other_step <- psc_control(7L)
  testthat::expect_identical(psc_strip(psc_view(a, ctx, g = 4, control = other_step)),
    psc_strip(psc_view(a, NULL, g = 4, control = other_step)))
  testthat::expect_equal(ctx$pair_support$scans, 2)
  for (other in list(psc_case(), psc_case(changed_weights = TRUE))) {
    testthat::expect_identical(psc_strip(psc_view(other, ctx, g = 4)), psc_strip(psc_view(other, NULL, g = 4)))
  }
  testthat::expect_equal(ctx$pair_support$scans, 2)
  bad_frame <- a; bad_frame$frame$w_a[1L] <- bad_frame$frame$w_a[1L] + 1
  e0 <- tryCatch(psc_view(bad_frame, NULL), error = identity)
  testthat::expect_s3_class(e0, "cdid_ipw_owner_integrity_error")
  testthat::expect_identical(tryCatch(psc_view(bad_frame, ctx), error = identity), e0)
  key <- "cp04-pair:5:1:2"
  raw <- get(key, envir = ctx$pair_support$records)
  changed <- raw; changed$pair_mass[1L] <- changed$pair_mass[1L] + 1
  assign(key, changed, envir = ctx$pair_support$records)
  testthat::expect_false(psc("ipw_owner_seal_ok")(changed, "cp04_pair_support"))
  before <- ctx$pair_support$reuses
  testthat::expect_identical(psc_strip(psc_view(a, ctx, g = 5, K = c(Inf, 3))),
    psc_strip(psc_view(a, NULL, g = 5, K = c(Inf, 3))))
  testthat::expect_equal(ctx$pair_support$reuses, before)
  testthat::expect_identical(get(key, envir = ctx$pair_support$records), changed)
})
