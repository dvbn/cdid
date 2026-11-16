# Focused tests for the private nuisance-model output proposal, initial slice:
# fitted observation-selection logits (R/chained_nuisance_models.R).
#
# Pure-helper and independent-reference tests run against the installed cdid
# namespace, or, when `nuisance_models_helpers` (an environment holding the
# package R sources) exists in the calling environment, against those sources.
# They drive the production pass/table code through the package's own
# pure-helper selection driver (evaluate_selection_table_if: the same
# partition_selection_cells -> fit_selection_block_slices ->
# evaluate_selection_block_if chain as the sealed backend) on one tiny
# deterministic fixture, and compare with an independent weighted glm() plus a
# hand-written unit/cluster sandwich. The installed-only tests run att_gt()
# end to end and are skipped when the installed cdid lacks nuisance_models().

nm_pure_mode <- function() exists("nuisance_models_helpers", inherits = TRUE)

nm_has_installed_cdid <- function() {
  if (nm_pure_mode()) return(FALSE)
  (isNamespaceLoaded("cdid") || requireNamespace("cdid", quietly = TRUE)) &&
    exists("nuisance_models", envir = asNamespace("cdid"), inherits = FALSE)
}

nm_ns <- function(name) {
  env <- if (nm_pure_mode()) get("nuisance_models_helpers", inherits = TRUE) else asNamespace("cdid")
  stopifnot(exists(name, envir = env, inherits = FALSE))
  get(name, envir = env, inherits = FALSE)
}

# ---- one tiny deterministic fixture (no RNG) ----------------------------------
# 60 units, cohorts 3 and never treated (0), dates 1:3, unequal integer weights,
# a factor covariate, 15 nested clusters, and covariate-dependent missing
# outcomes from a fixed hash (every unit keeps at least one observed date).
# `N` scales the same generator (the augmented DR+DR fixture uses N = 240).
nm_fixture <- function(N = 60L) {
  id <- seq_len(N)
  g <- ifelse(id <= N %/% 2L, 3, 0)
  x <- ((id * 37L) %% 23L - 11) / 5
  region <- ifelse(id %% 3L == 0L, "B", "A")
  w <- as.double(1L + id %% 3L)
  cl <- ceiling(id / 4)
  calendar <- c(1, 2, 3)
  u <- outer(id, calendar, function(i, t) ((i * 7919 + t * 104729) %% 1000 + 0.5) / 1000)
  pr <- stats::plogis(outer(0.9 + 0.7 * x + 0.5 * (region == "B"), calendar, function(a, t) a - 0.35 * t))
  obs <- (u < pr) * 1
  obs[rowSums(obs) == 0, 1L] <- 1
  frame <- data.frame(id = id, g = g, x = x, region = region, w = w, cl = cl, stringsAsFactors = FALSE)
  long <- do.call(rbind, lapply(calendar, function(t) {
    keep <- obs[, t] == 1
    data.frame(id = id[keep], t = t, y = 0.5 * x[keep] + 0.2 * t + (g[keep] == 3 & t == 3) +
      ((id[keep] * 13 + t * 7) %% 10) / 10, g = g[keep], x = x[keep], region = region[keep], w = w[keep],
      cl = cl[keep], stringsAsFactors = FALSE)
  }))
  long <- long[order(long$id, long$t), ]
  rownames(long) <- NULL
  list(N = N, frame = frame, obs = obs, calendar = calendar, data = long)
}

nm_fx <- nm_fixture()

# Independent reference for one fitted selection block: weighted glm on the
# stacked (unit x cell) rows and a hand-written sandwich with scores summed by
# unit (or cluster). `cells`: data frame with type, cohort (public), s, t.
nm_reference_block <- function(fx, cells, intercept = TRUE, slopes = c("x", "regionB")) {
  fr <- fx$frame
  rows <- do.call(rbind, lapply(seq_len(nrow(cells)), function(j) {
    ii <- which(fr$g == cells$cohort[j])
    y <- if (identical(cells$type[j], "date")) fx$obs[ii, match(cells$t[j], fx$calendar)] else
      fx$obs[ii, match(cells$s[j], fx$calendar)] * fx$obs[ii, match(cells$t[j], fx$calendar)]
    data.frame(unit = ii, cell = j, y = y, w = fr$w[ii], x = fr$x[ii],
      regionB = as.double(fr$region[ii] == "B"), cl = fr$cl[ii])
  }))
  rows$cell <- factor(rows$cell, levels = seq_len(nrow(cells)))
  rhs <- c(if (intercept && nrow(cells) > 1L) "0 + cell" else if (intercept) "1" else "0", slopes)
  f <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))
  fit <- stats::glm(f, family = stats::binomial(), data = rows, weights = w,
    control = stats::glm.control(epsilon = 1e-15, maxit = 100L))
  X <- stats::model.matrix(fit); p <- stats::fitted(fit)
  Ainv <- solve(crossprod(X, X * (rows$w * p * (1 - p))))
  S <- X * (rows$w * (rows$y - p))
  sw <- function(G) sqrt(diag(Ainv %*% crossprod(rowsum(S, G)) %*% Ainv))
  list(coef = unname(stats::coef(fit)), se_unit = unname(sw(rows$unit)), se_cluster = unname(sw(rows$cl)),
    se_row = unname(sqrt(diag(Ainv %*% crossprod(S) %*% Ainv))), se_model = unname(sqrt(diag(Ainv))))
}

# Production pass/table code on the pure-helper selection driver.
nm_table_capture <- function(fx, selection_model, selection_pooling, B = 7L, clustered = FALSE,
                             formula = ~ x + region) {
  N <- fx$N; fr <- fx$frame
  Z <- stats::model.matrix(formula, fr)
  attr(Z, "assign") <- NULL; attr(Z, "contrasts") <- NULL
  attr(Z, "has_intercept") <- identical(attr(stats::terms(formula), "intercept"), 1L)
  g <- ifelse(fr$g == 0, Inf, fr$g)
  cohorts <- sort(unique(g)); ci <- match(g, cohorts)
  pairs <- list(c(1, 2), c(2, 3), c(1, 3))
  mk_pair <- function(pq) list(s = pq[1], t = pq[2], s_calendar_index = pq[1], t_calendar_index = pq[2])
  pair_specs <- unlist(lapply(pairs, function(pq) lapply(1:2, function(k) list(cohort = k, pair = mk_pair(pq)))), recursive = FALSE)
  pair_resp <- lapply(pairs, function(pq) fx$obs[, pq[1]] * fx$obs[, pq[2]])
  names(pair_resp) <- vapply(pairs, function(pq) paste0(pq[1], ":", pq[2]), "")
  tab <- nm_ns("evaluate_selection_table_if")(Z, g, fr$w, as.double(rowSums(fx$obs) >= 2), cohorts = cohorts,
    all_times = fx$calendar, date_responses = fx$obs, pair_responses = pair_resp, pair_specs = pair_specs,
    selection_model = selection_model, selection_pooling = selection_pooling,
    glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE), details = TRUE)
  class(tab) <- c("cdid_selection_models_table", class(tab))
  up_map <- nm_ns("build_selection_upstream_map")(tab)
  up_all <- nm_ns("selection_table_upstream_if")(tab, up_map)
  tab$links <- lapply(seq_along(pairs), function(j) list(link_id = paste0("cp04:catalog:", j), cohorts = c(1L, 2L),
    pair = mk_pair(pairs[[j]]), pair_cohort_ids = vapply(1:2, function(k)
      nm_ns("resolve_pair_cohort_id")(tab, list(cohort = k, pair = mk_pair(pairs[[j]]))), integer(1))))
  blk <- nm_ns("chained_nuisance_block_geometry")(tab, up_map)
  cnt <- tabulate(ci, length(cohorts)); ends <- cumsum(as.double(cnt))
  st <- list(cohort_start = ends - cnt + 1, cohort_end = ends, cohort_order = order(ci, method = "radix"), U = rep(TRUE, N))
  ds <- tab$design_spec
  fns <- list(
    influence = function(rows, m) up_all[rows, blk[[m]]$uids, drop = FALSE],
    cell = function(rows, cid) {
      cell <- tab$cells[[cid]]; bl <- tab$blocks[[cell$block_id]]
      eta <- (if (bl$has_intercept) bl$alpha[cell$local_b] else 0) + drop(Z[rows, ds$slope_cols, drop = FALSE] %*% bl$beta)
      y <- if (cell$type == "date") fx$obs[rows, cell$calendar_index] else pair_resp[[paste0(cell$pair$s, ":", cell$pair$t)]][rows]
      list(p = nm_ns("ipw_clip_logit")(eta)$prob, y = y)
    },
    pair = function(rows, q) nm_ns("predict_table_pair_selection")(tab, Z, fr$w, rep(1, N), rows, q)$prob)
  group <- if (clustered) match(fr$cl, unique(fr$cl)) else NULL
  ng <- if (clustered) max(group) else N
  stats_ <- nm_ns("chained_nuisance_selection_passes")(tab, blk, st, B, N, fns, group, ng)
  meta <- list(cohorts = cohorts, mass_U = as.double(tapply(fr$w, ci, sum)) / N, N = N, design_columns = colnames(Z),
    universe = "F", formula = paste(deparse(formula), collapse = " "), population = "fixture frame",
    weights = list(text = "w_i = w", column = "w", estimated = FALSE), transport = "rho = 1",
    glm_control = list(epsilon = 1e-13, maxit = 60L), route = "pure-helper table driver",
    covariance = nm_ns("chained_nuisance_covariance_record")(TRUE, NULL, NULL, clustered, if (clustered) "cl", ng, N),
    accounting = list())
  list(tab = tab, comp = nm_ns("chained_nuisance_selection_tables")(tab, blk, stats_, meta))
}

# Legacy separate pair-logit records as fit_cp04_corollary1() holds them: one
# (pair, cohort) spec per fit, coefficients over every design column and the
# N-normalized factor chol(H/N). The records come from a hand-written Newton
# logit (not glm), the oracle (nm_reference_block) refits glm independently.
nm_newton_logit <- function(X, y, w) {
  b <- numeric(ncol(X))
  for (it in seq_len(100L)) {
    p <- stats::plogis(drop(X %*% b))
    step <- drop(solve(crossprod(X, X * (w * p * (1 - p))), crossprod(X, w * (y - p))))
    b <- b + step
    if (max(abs(step)) < 1e-13) break
  }
  list(b = b, p = stats::plogis(drop(X %*% b)), iter = it)
}

nm_legacy_capture <- function(fx, formula = ~ x + region, perm = NULL, B = 7L, clustered = FALSE) {
  N <- fx$N; fr <- fx$frame
  Z <- stats::model.matrix(formula, fr)
  attr(Z, "assign") <- NULL; attr(Z, "contrasts") <- NULL
  if (!is.null(perm)) Z <- Z[, perm, drop = FALSE]
  has_int <- identical(attr(stats::terms(formula), "intercept"), 1L)
  g <- ifelse(fr$g == 0, Inf, fr$g)
  cohorts <- sort(unique(g)); ci <- match(g, cohorts)
  pairs <- list(c(1, 2), c(2, 3), c(1, 3))
  mk_pair <- function(pq) list(s = pq[1], t = pq[2], s_calendar_index = pq[1], t_calendar_index = pq[2])
  specs <- unlist(lapply(pairs, function(pq) lapply(1:2, function(k) list(pair = mk_pair(pq), cohort = k))),
    recursive = FALSE)
  resp <- function(rows, j) fx$obs[rows, specs[[j]]$pair$s] * fx$obs[rows, specs[[j]]$pair$t]
  mass_U <- as.double(tapply(fr$w, ci, sum)) / N
  fits <- lapply(seq_along(specs), function(j) {
    rows <- which(ci == specs[[j]]$cohort); X <- Z[rows, , drop = FALSE]; y <- resp(rows, j); w <- fr$w[rows]
    nf <- nm_newton_logit(X, y, w)
    list(coefficients = nf$b, H_chol = chol(crossprod(X, X * (w * nf$p * (1 - nf$p))) / N),
      diagnostics = list(converged = TRUE, boundary = FALSE, rank = ncol(X), iter = nf$iter),
      pi = sum(w * y) / (mass_U[specs[[j]]$cohort] * N))
  })
  # One comparison per pair, target g = 3 against the never-treated, interned
  # as fit_cp04_corollary1() does: cohorts (target, controls) and selection_ids.
  b <- structure(list(selection_specs = specs, selections = fits,
    treatment_specs = lapply(seq_along(pairs), function(q) list(pair = mk_pair(pairs[[q]]), cohorts = c(1L, 2L),
      selection_ids = c(2L * q - 1L, 2L * q), n_observers = 1L)),
    links = lapply(seq_along(pairs), function(q) list(model = q, link_id = paste0("cp04:catalog:", q)))),
    class = "cdid_cp04_corollary1")
  links <- nm_ns("chained_nuisance_legacy_links")(b)
  view <- nm_ns("chained_nuisance_legacy_view")(specs, fits, cohorts, colnames(Z), has_int, links, "F")
  blk <- nm_ns("chained_nuisance_block_geometry")(view$sel, view$up_map)
  n_predict <- 0L
  predict <- function(rows, j) {
    n_predict <<- n_predict + 1L
    stats::plogis(drop(Z[rows, , drop = FALSE] %*% fits[[j]]$coefficients))
  }
  fns <- nm_ns("chained_nuisance_legacy_fns")(Z, fr$w, lapply(fits, `[[`, "H_chol"), view$coord_order, predict, resp)
  cnt <- tabulate(ci, length(cohorts)); ends <- cumsum(as.double(cnt))
  st <- list(cohort_start = ends - cnt + 1, cohort_end = ends, cohort_order = order(ci, method = "radix"), U = rep(TRUE, N))
  group <- if (clustered) match(fr$cl, unique(fr$cl)) else NULL
  ng <- if (clustered) max(group) else N
  stats_ <- nm_ns("chained_nuisance_selection_passes")(view$sel, blk, st, B, N, fns, group, ng)
  meta <- list(cohorts = cohorts, mass_U = mass_U, N = N, design_columns = colnames(Z),
    universe = "F", formula = paste(deparse(formula), collapse = " "), population = "fixture frame",
    weights = list(text = "w_i = w", column = "w", estimated = FALSE), transport = "rho = 1",
    glm_control = list(), route = "legacy records on the fixture",
    covariance = nm_ns("chained_nuisance_covariance_record")(TRUE, NULL, NULL, clustered, if (clustered) "cl", ng, N),
    accounting = list())
  list(b = b, view = view, fits = fits, n_predict = n_predict,
    comp = nm_ns("chained_nuisance_selection_tables")(view$sel, blk, stats_, meta))
}

nm_fake_fit <- function(comp, est_label = "ipw") {
  cap <- nm_ns("chained_nuisance_capture_assemble")(comp,
    list(selection_est_method = "ipw", selection = "conditional_mtar"), est_label)
  structure(list(chained = list(nuisance_models = cap)), class = c("cdid_MP", "MP"))
}

nm_print <- function(x) {
  cls <- class(x)[1L]
  utils::capture.output(nm_ns(paste0("print.", cls))(x))
}

# ---------------------------------------------------------------------------------
test_that("control$keep_nuisance is consumed before the strict normalizer; other entries are unchanged", {
  ctl <- nm_ns("chained_nuisance_control")
  resolve <- nm_ns("chained_public_resolve_control")
  expect_identical(ctl(NULL), list(keep = FALSE, control = NULL, supplied = FALSE))
  plain <- list(workers = 2, chunk_rows = 64)
  expect_identical(ctl(plain)$control, plain)
  expect_false(ctl(plain)$keep)
  k <- ctl(list(workers = 2, keep_nuisance = TRUE, chunk_rows = 64))
  expect_true(k$keep)
  expect_identical(k$control, plain)
  expect_identical(resolve(k$control), resolve(plain))
  # Only the opt-in entry: the remainder is exactly the user's empty list.
  expect_identical(resolve(ctl(list(keep_nuisance = TRUE))$control), resolve(list()))
  expect_false(ctl(list(keep_nuisance = FALSE))$keep)
  # The normalizer itself stays strict (no new accepted name).
  expect_error(resolve(list(keep_nuisance = TRUE)), class = "cdid_chained_public_control_error")
  expect_error(resolve(ctl(list(keep_nuisance = TRUE, unknown_entry = 1))$control),
    class = "cdid_chained_public_control_error")
  for (bad in list(NA, "yes", c(TRUE, TRUE), 1)) {
    expect_error(ctl(list(keep_nuisance = bad)), class = "cdid_chained_public_control_error")
  }
  expect_error(ctl(list(keep_nuisance = TRUE, keep_nuisance = FALSE)), class = "cdid_chained_public_control_error")
})

test_that("selection SEs reuse the block IF and equal an independent unit and cluster sandwich", {
  for (cfg in list(c("pair_logit", "within_cohort"), c("independent_period_logit", "within_cohort"),
                   c("pair_logit", "separate"), c("pair_logit", "across_cohorts"),
                   c("independent_period_logit", "separate"))) {
    r <- nm_table_capture(nm_fx, cfg[1], cfg[2])
    r_cl <- nm_table_capture(nm_fx, cfg[1], cfg[2], B = 5L, clustered = TRUE)
    r_one <- nm_table_capture(nm_fx, cfg[1], cfg[2], B = nm_fx$N)
    r_tiny <- nm_table_capture(nm_fx, cfg[1], cfg[2], B = 1L)
    co <- r$comp$coefficients
    expect_true(all(is.finite(co$std_error)))
    # Chunking never changes the accumulated diagonal.
    expect_equal(r_one$comp$coefficients$std_error, co$std_error, tolerance = 1e-12)
    expect_equal(r_tiny$comp$coefficients$std_error, co$std_error, tolerance = 1e-12)
    expect_identical(r_cl$comp$coefficients$estimate, co$estimate)
    for (m in seq_along(r$tab$blocks)) {
      id <- paste0("S", m)
      cells <- r$comp$cells[r$comp$cells$block == id, , drop = FALSE]
      ref <- nm_reference_block(nm_fx, cells)
      cm <- co[co$block == id, , drop = FALSE]
      cc <- r_cl$comp$coefficients[r_cl$comp$coefficients$block == id, , drop = FALSE]
      expect_equal(cm$estimate, ref$coef, tolerance = 1e-8)
      expect_equal(cm$std_error, ref$se_unit, tolerance = 1e-8)
      expect_equal(cc$std_error, ref$se_cluster, tolerance = 1e-8)
      expect_equal(cm$z_value, cm$estimate / cm$std_error, tolerance = 1e-14)
      # The reference discriminates: a model-based inverse Hessian is not the reported SE.
      expect_true(max(abs(ref$se_model / ref$se_unit - 1)) > 1e-2)
      # Pooled blocks reuse a unit across cells: row-level scores would be wrong.
      if (nrow(cells) > 1L) expect_true(max(abs(ref$se_row / ref$se_unit - 1)) > 1e-2)
    }
  }
})

test_that("equations and tables carry the actual response, cells, terms, sample and weights", {
  r <- nm_table_capture(nm_fx, "pair_logit", "within_cohort")
  comp <- r$comp
  expect_identical(comp$template, "logit Pr(S_i(s,t) = 1 | X_i, G_i = g, U_i = 1) = a[g,(s,t)] + X_i'b[g]")
  expect_identical(comp$response, "S_i(s,t) = 1{Y_is and Y_it both observed}")
  expect_identical(comp$slope_terms, c("x", "regionB"))
  expect_identical(comp$blocks$label, c("within_cohort[g=3]", "within_cohort[g=0]"))
  expect_identical(comp$blocks$equation[1L], "logit Pr(S_i(s,t) = 1 | X_i, G_i = 3, U_i = 1) = a[(s,t)] + b1*x + b2*regionB")
  expect_identical(comp$coefficients$term[comp$coefficients$block == "S1"],
    c("(Intercept)[g=3,(1,2)]", "(Intercept)[g=3,(1,3)]", "(Intercept)[g=3,(2,3)]", "x", "regionB"))
  expect_identical(comp$coefficients$role[comp$coefficients$block == "S1"],
    c(rep("local_intercept", 3L), rep("shared_slope", 2L)))
  expect_equal(comp$blocks$n_units, c(30, 30))
  expect_equal(comp$blocks$n_fitting_rows, c(90, 90))
  expect_equal(comp$blocks$n_comparisons, c(3, 3))
  expect_equal(comp$cells$n_observed[1L], sum(nm_fx$obs[nm_fx$frame$g == 3, 1] * nm_fx$obs[nm_fx$frame$g == 3, 2]))
  expect_null(comp$derived)

  d <- nm_table_capture(nm_fx, "independent_period_logit", "across_cohorts")$comp
  expect_identical(d$template, "logit Pr(R_it = 1 | X_i, G_i = g) = a[g,t] + X_i'b")
  expect_identical(d$blocks$label, "across_cohorts[g in {0,3}]")
  expect_match(d$derived, "transformation of the date logits, not a second fit", fixed = TRUE)
  expect_true(all(grepl("derived product", d$pair_cohorts$probability, fixed = TRUE)))
  expect_identical(d$pair_cohorts$blocks, rep("S1", 6L))
  # The derived product never exceeds either date probability on the fitting sample.
  expect_true(all(d$pair_cohorts$p_min <= min(d$cells$p_max)))
  sep <- nm_table_capture(nm_fx, "independent_period_logit", "separate")$comp
  expect_identical(sep$blocks$equation[2L], "logit Pr(R_i[t=2] = 1 | X_i, G_i = 3) = a + b1*x + b2*regionB")
})

test_that("intercept-only, no-intercept and zero-column designs are reported as fitted", {
  for (cfg in list(c("pair_logit", "within_cohort"), c("independent_period_logit", "across_cohorts"))) {
    io <- nm_table_capture(nm_fx, cfg[1], cfg[2], formula = ~ 1)$comp
    expect_identical(io$design_branch, "intercept_only")
    expect_false(any(grepl("X_i'b", io$template, fixed = TRUE)))
    expect_true(all(io$coefficients$role == "local_intercept"))
    ni <- nm_table_capture(nm_fx, cfg[1], cfg[2], formula = ~ 0 + x)$comp
    expect_identical(ni$design_branch, "no_intercept")
    expect_identical(unique(ni$coefficients$term), "x")
    expect_true(all(endsWith(ni$blocks$equation, "= b1*x")))
    for (m in seq_len(nrow(io$blocks))) {
      id <- io$blocks$block[m]
      ref_io <- nm_reference_block(nm_fx, io$cells[io$cells$block == id, , drop = FALSE], slopes = character(0))
      expect_equal(io$coefficients$estimate[io$coefficients$block == id], ref_io$coef, tolerance = 1e-8)
      expect_equal(io$coefficients$std_error[io$coefficients$block == id], ref_io$se_unit, tolerance = 1e-8)
      ref_ni <- nm_reference_block(nm_fx, ni$cells[ni$cells$block == id, , drop = FALSE], intercept = FALSE, slopes = "x")
      expect_equal(ni$coefficients$estimate[ni$coefficients$block == id], ref_ni$coef, tolerance = 1e-8)
      expect_equal(ni$coefficients$std_error[ni$coefficients$block == id], ref_ni$se_unit, tolerance = 1e-8)
    }
    zc <- nm_table_capture(nm_fx, cfg[1], cfg[2], formula = ~ 0)$comp
    expect_identical(zc$design_branch, "zero_col")
    expect_identical(nrow(zc$coefficients), 0L)
    expect_true(all(grepl("= 0  (fixed: Pr = 1/2", zc$blocks$equation, fixed = TRUE)))
    expect_true(all(zc$cells$p_min == 0.5 & zc$cells$p_max == 0.5))
    out <- nm_print(nm_ns("nuisance_models")(nm_fake_fit(zc), block = 1L))
    expect_true(any(grepl("No fitted coefficients", out, fixed = TRUE)))
  }
})

test_that("nuisance_models() gives an index, component and block tables, and typed statuses", {
  nm <- nm_ns("nuisance_models")
  comp <- nm_table_capture(nm_fx, "pair_logit", "within_cohort")$comp
  fit <- nm_fake_fit(comp)
  idx <- nm(fit)
  expect_s3_class(idx, "cdid_nuisance_index")
  expect_identical(idx$components$status, c("available", "not_applicable", "not_captured", "not_applicable"))
  expect_identical(idx$blocks$block, c("S1", "S2"))
  out <- nm_print(idx)
  expect_true(any(grepl("a[g,(s,t)] + X_i'b[g]", out, fixed = TRUE)))
  expect_true(any(grepl("treatment: not_captured", out, fixed = TRUE)))

  b1 <- nm(fit, "selection", block = "S1")
  expect_s3_class(b1, "cdid_nuisance_block")
  expect_identical(colnames(b1$coefficients), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_identical(rownames(b1$coefficients), comp$coefficients$term[comp$coefficients$block == "S1"])
  expect_identical(nm(fit, block = 1L)$coefficients, b1$coefficients)
  expect_identical(nm(fit, "selection", block = "within_cohort[g=3]")$coefficients, b1$coefficients)
  expect_identical(b1$intercept_cells, c("g=3,(1,2)", "g=3,(1,3)", "g=3,(2,3)"))
  bo <- nm_print(b1)
  expect_true(any(grepl("logit Pr(S_i(s,t) = 1 | X_i, G_i = 3, U_i = 1) = a[(s,t)] + b1*x + b2*regionB", bo, fixed = TRUE)))
  expect_true(any(grepl("one score contribution per original unit", bo, fixed = TRUE)))
  expect_s3_class(nm(fit, "selection"), "cdid_nuisance_component")
  expect_true(length(nm_print(nm(fit, "selection"))) > length(bo))

  st <- nm(fit, "outcome")
  expect_s3_class(st, "cdid_nuisance_status")
  expect_identical(st$status, "not_applicable")
  expect_error(nm(fit, "selection", block = "S9"), class = "cdid_nuisance_input_error")
  expect_error(nm(fit, "propensity"), class = "cdid_nuisance_input_error")

  # Fits without a capture, ordinary delegation, aggregations and other objects.
  expect_identical(nm(structure(list(chained = list()), class = c("cdid_MP", "MP")))$code, "not_requested")
  expect_identical(nm(structure(list(), class = "MP"))$code, "ordinary_delegation")
  expect_identical(nm(structure(list(), class = c("cdid_AGGTEobj", "AGGTEobj")))$code, "aggregation_object")
  expect_error(nm(list()), class = "cdid_nuisance_input_error")
  expect_true(any(grepl("keep_nuisance = TRUE", nm_print(nm(structure(list(chained = list()),
    class = c("cdid_MP", "MP")))), fixed = TRUE)))
})

test_that("routes outside the initial slice return typed statuses without touching the backend", {
  scope <- nm_ns("chained_nuisance_selection_scope")
  gate <- list(selection = "conditional_mtar", selection_est_method = "ipw", is_custom_est = FALSE)
  shared <- structure(list(uses_selection_models = TRUE, selection_backend = list(blocks = list(list(), list())),
    up_map = list(theta_coords = list(1:2, 3:4))), class = "cdid_dr_selection_backend")
  expect_null(scope(shared, list(), gate))
  expect_identical(scope(structure(list(uses_selection_models = TRUE, selection_backend = list(blocks = list(list()))),
    class = "cdid_or_selection_backend"), list(), gate)$code, "unrecognized_backend")
  expect_identical(scope(NULL, list(), utils::modifyList(gate, list(selection = "none")))$code, "no_selection_model")
  expect_identical(scope(shared, list(), utils::modifyList(gate, list(is_custom_est = TRUE)))$code, "custom_callback_route")
  expect_identical(scope(NULL, list(), utils::modifyList(gate, list(selection_est_method = "imputation")))$code, "pure_imputation")
  expect_identical(scope(NULL, list(), utils::modifyList(gate, list(selection_est_method = "dr")))$code, "no_full_roster_backend")
  expect_identical(scope(shared, list(reference_source = list()), gate)$code, "reference_recipe_route")
  # Legacy consumers are routed to their own predicate (tested below); a legacy
  # class without its fit records is unrecognized, a non-legacy class keeps its code.
  legacy <- structure(list(selections = list()), class = "cdid_cp04_corollary1")
  expect_identical(scope(legacy, list(), gate)$code, "unrecognized_backend")
  expect_identical(scope(structure(list(), class = "cdid_dr_corollary1"), list(), gate)$code,
    "legacy_separate_pair_consumer")
  expect_identical(scope(structure(list(uses_selection_models = TRUE), class = "cdid_or_selection_backend"),
    list(), gate)$code, "unrecognized_backend")
  other <- nm_ns("chained_nuisance_other_components")
  expect_identical(vapply(other(gate, "dr"), `[[`, "", "status"),
    c(outcome = "not_captured", treatment = "not_captured", imputation = "not_applicable"))
  expect_identical(other(gate, "reg")$treatment$status, "not_applicable")
})

test_that("augmented observation logits: only the ordinary static shared backend is eligible", {
  scope <- nm_ns("chained_nuisance_selection_scope")
  gate <- list(selection = "conditional_mtar", selection_est_method = "dr", is_custom_est = FALSE)
  env <- new.env(parent = emptyenv())
  assign("cp04_reader", list(state = list(universe = "F", id = 1:4, U = rep(TRUE, 4L), rho = rep(1, 4L))), env)
  aug <- structure(list(uses_selection_models = TRUE, universe = "F", reader_ref = env,
    selection_backend = structure(list(blocks = list(list(), list()), reader_ref = env), class = "cdid_selection_models"),
    up_map = list(theta_coords = list(1:3, 4:6))), class = c("cdid_augmented_backend", "cdid_chained_augmented"))
  prep <- list(frame_ids = 1:4)
  code <- function(b, p = prep, g = gate) scope(b, p, g)$code
  mod <- function(...) utils::modifyList(aug, list(...))
  expect_null(scope(aug, prep, gate))
  # Recipe/reference layout: different weights W_g and layout rows.
  expect_identical(code(mod(layout = list(N_layout = 8))), "reference_recipe_route")
  expect_identical(code(mod(reference_ref = new.env())), "reference_recipe_route")
  expect_identical(code(aug, c(prep, list(reference_source = list()))), "reference_recipe_route")
  expect_identical(code(aug, g = utils::modifyList(gate, list(is_custom_est = TRUE))), "custom_callback_route")
  expect_identical(code(mod(pair_outcomes = list(`1:2` = c(0, 1, NA, 2)))), "augmented_route")
  expect_identical(code(structure(unclass(aug), class = "cdid_imputation_backend")), "unrecognized_backend")
  short_map <- aug
  short_map$up_map$theta_coords <- list(1:3)
  for (bad in list(mod(uses_selection_models = FALSE), mod(selection_backend = list(reader_ref = new.env())),
                   mod(selection_backend = list(reference_ref = new.env())), short_map))
    expect_identical(code(bad), "unrecognized_backend")
  roster <- scope(aug, list(frame_ids = c(1:3, 5L)), gate)
  expect_identical(roster$code, "unrecognized_backend")
  expect_match(roster$reason, "U = F with rho = 1", fixed = TRUE)

  other <- nm_ns("chained_nuisance_other_components")
  dr <- other(gate, "dr")
  expect_identical(vapply(dr, `[[`, "", "status"),
    c(outcome = "not_captured", treatment = "not_captured", imputation = "not_captured"))
  expect_match(dr$outcome$reason, "reconstructed response Z_c", fixed = TRUE)
  expect_match(dr$imputation$reason, "predictor in Z_c", fixed = TRUE)
  expect_identical(other(gate, "ipw")$outcome$status, "not_applicable")
  expect_identical(other(gate, "reg")$treatment$status, "not_applicable")
})

test_that("legacy separate pair logits: stored-factor SEs equal an independent unit and cluster sandwich", {
  for (perm in list(NULL, c(2L, 1L, 3L))) {
    r <- nm_legacy_capture(nm_fx, perm = perm)
    r_cl <- nm_legacy_capture(nm_fx, perm = perm, B = 5L, clustered = TRUE)
    r_one <- nm_legacy_capture(nm_fx, perm = perm, B = nm_fx$N)
    co <- r$comp$coefficients
    expect_true(all(is.finite(co$std_error)))
    expect_equal(r_one$comp$coefficients$std_error, co$std_error, tolerance = 1e-12)
    expect_identical(r_cl$comp$coefficients$estimate, co$estimate)
    # Each fitting chunk is predicted once and shared by the IF and the cell summary.
    expect_identical(r$n_predict, 6L * as.integer(ceiling(30 / 7)))
    expect_identical(nrow(r$comp$blocks), 6L)
    for (m in seq_len(6L)) {
      id <- paste0("S", m)
      cells <- r$comp$cells[r$comp$cells$block == id, , drop = FALSE]
      expect_identical(nrow(cells), 1L)
      ref <- nm_reference_block(nm_fx, cells)
      cm <- co[co$block == id, , drop = FALSE]
      cc <- r_cl$comp$coefficients[r_cl$comp$coefficients$block == id, , drop = FALSE]
      # Exactly the stored record, shown as (intercept, slopes) whatever the column order.
      zc <- r$view$coord_order
      expect_identical(cm$estimate, as.double(r$fits[[m]]$coefficients[zc]))
      expect_identical(cm$term, c(paste0("(Intercept)[", cells$key, "]"), "x", "regionB"))
      expect_identical(cm$role, c("local_intercept", "slope", "slope"))
      expect_equal(cm$estimate, ref$coef, tolerance = 1e-8)
      expect_equal(cm$std_error, ref$se_unit, tolerance = 1e-8)
      expect_equal(cc$std_error, ref$se_cluster, tolerance = 1e-8)
      expect_true(max(abs(ref$se_model / ref$se_unit - 1)) > 1e-2)
      expect_true(max(abs(ref$se_cluster / ref$se_unit - 1)) > 1e-2)
    }
  }
  # No-intercept legacy design.
  ni <- nm_legacy_capture(nm_fx, formula = ~ 0 + x)$comp
  expect_identical(ni$design_branch, "no_intercept")
  ref_ni <- nm_reference_block(nm_fx, ni$cells[1L, , drop = FALSE], intercept = FALSE, slopes = "x")
  expect_equal(ni$coefficients$std_error[ni$coefficients$block == "S1"], ref_ni$se_unit, tolerance = 1e-8)
})

test_that("legacy separate pair logits: equations, comparison map and scope typing", {
  r <- nm_legacy_capture(nm_fx)
  comp <- r$comp
  expect_identical(comp$template, "logit Pr(S_i(s,t) = 1 | X_i, G_i = g, U_i = 1) = a[g,(s,t)] + X_i'b[g,(s,t)]")
  expect_identical(comp$blocks$label[1:2], c("pair[g=3,(1,2)]", "pair[g=0,(1,2)]"))
  expect_identical(comp$blocks$equation[1L], "logit Pr(S_i(1,2) = 1 | X_i, G_i = 3, U_i = 1) = a + b1*x + b2*regionB")
  expect_equal(comp$blocks$n_units, rep(30, 6L))
  expect_equal(comp$blocks$n_comparisons, rep(1, 6L))
  expect_identical(comp$comparisons$controls, rep("0", 3L))
  expect_identical(comp$comparisons$pair_cohorts, c("1,2", "3,4", "5,6"))
  expect_equal(comp$cells$n_observed[2L], sum(nm_fx$obs[nm_fx$frame$g == 0, 1] * nm_fx$obs[nm_fx$frame$g == 0, 2]))

  # OR: target support plus control-pool supports give the same map.
  links <- nm_ns("chained_nuisance_legacy_links")
  b <- r$b
  or_b <- structure(list(selection_specs = b$selection_specs, selections = b$selections,
    support_specs = lapply(b$selection_specs, function(s) c(s, list(transitions = 1L))),
    regression_specs = lapply(1:3, function(q) list(support_ids = 2L * q)),
    links = lapply(1:3, function(q) list(link_id = paste0("cp04:catalog:", q), cohort = 1L,
      support_id = 2L * q - 1L, regression = q))), class = "cdid_or_corollary1")
  expect_identical(links(or_b), links(b))
  bad <- or_b; bad$regression_specs[[1L]]$support_ids <- 4L # control fit of another pair
  expect_null(links(bad))
  bad <- b; bad$treatment_specs[[1L]]$selection_ids <- c(1L, 7L)
  expect_null(links(bad))

  scope <- nm_ns("chained_nuisance_selection_scope")
  gate <- list(selection = "conditional_mtar", selection_est_method = "ipw", is_custom_est = FALSE)
  env <- new.env(parent = emptyenv())
  assign("cp04_reader", list(state = list(Z = matrix(1, 2L, 3L))), env)
  for (x in list(b, or_b)) {
    x$reader_ref <- env
    expect_null(scope(x, list(), gate))
    expect_identical(scope(utils::modifyList(x, list(reference_ref = new.env())), list(), gate)$code,
      "reference_recipe_route")
    short <- x; short$selections[[2L]]$coefficients <- 1
    expect_identical(scope(short, list(), gate)$code, "unrecognized_backend")
    expect_identical(scope(x, list(reference_source = list()), gate)$code, "reference_recipe_route")
    expect_identical(scope(x, list(), utils::modifyList(gate, list(is_custom_est = TRUE)))$code,
      "custom_callback_route")
  }
})

test_that("the capture is compact and serializable (no environments, functions or owners)", {
  comp <- nm_table_capture(nm_fx, "independent_period_logit", "within_cohort", clustered = TRUE)$comp
  cap <- nm_fake_fit(comp)$chained$nuisance_models
  has_ref <- function(x) is.environment(x) || is.function(x) || is.externalptr(x) ||
    (is.list(x) && any(vapply(x, has_ref, logical(1))))
  is.externalptr <- function(x) identical(typeof(x), "externalptr")
  expect_false(has_ref(unclass(cap)))
  raw <- serialize(cap, NULL)
  expect_identical(unserialize(raw), cap)
  expect_true(length(raw) < 65536)
})

# ---------------------------------------------------------------------------------
# Installed package: att_gt() end to end.
# ---------------------------------------------------------------------------------
nm_att_args <- function(fx, est_method, selection_model, selection_pooling, ...) {
  utils::modifyList(list(yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x + region, data = fx$data,
    weightsname = "w", est_method = est_method, selection = "conditional_mtar",
    selection_model = selection_model, selection_pooling = selection_pooling,
    bstrap = TRUE, cband = TRUE, biters = 40L), list(...))
}

nm_routes <- list(
  list(em = "ipw", sm = "independent_period_logit", sp = "within_cohort"),
  list(em = "reg", sm = "pair_logit", sp = "within_cohort"),
  list(em = "dr", sm = "pair_logit", sp = "separate"),
  list(em = "ipw", sm = "pair_logit", sp = "across_cohorts"))

# The built-in default selection_model/pooling (pair_logit, separate) with a
# nonempty design: the legacy per-(pair, cohort) consumers.
nm_legacy_routes <- list(
  list(em = "ipw", sm = "pair_logit", sp = "separate", class = "cdid_cp04_corollary1"),
  list(em = "reg", sm = "pair_logit", sp = "separate", class = "cdid_or_corollary1"))

test_that("installed att_gt: opt-in capture leaves estimates, IF, inference and RNG unchanged", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  att <- getExportedValue("cdid", "att_gt")
  agg <- getExportedValue("cdid", "aggte")
  for (rt in c(nm_routes, nm_legacy_routes)) {
    args <- nm_att_args(nm_fx, rt$em, rt$sm, rt$sp)
    set.seed(20261115L)
    f0 <- do.call(att, c(args, list(control = list())))
    a0 <- agg(f0, type = "dynamic"); s0 <- .Random.seed
    set.seed(20261115L)
    f1 <- do.call(att, c(args, list(control = list(keep_nuisance = TRUE))))
    a1 <- agg(f1, type = "dynamic"); s1 <- .Random.seed
    expect_identical(s1, s0)
    for (nm_f in c("att", "se", "inffunc", "V_analytical", "c", "n", "group", "t", "W", "Wpval"))
      expect_identical(f1[[nm_f]], f0[[nm_f]])
    expect_identical(a1$overall.att, a0$overall.att)
    expect_identical(a1$overall.se, a0$overall.se)
    expect_identical(a1$att.egt, a0$att.egt)
    expect_identical(f1$chained$selection, f0$chained$selection)
    expect_identical(f1$chained$nuisance_diagnostics, f0$chained$nuisance_diagnostics)
    expect_identical(f1$chained$public, f0$chained$public)
    expect_null(f0$chained$nuisance_models)
    expect_identical(f1$chained$nuisance_models$components$selection$status, "available")
    expect_identical(getExportedValue("cdid", "nuisance_models")(f0)$code, "not_requested")
  }
})

test_that("installed att_gt: captured coefficients are the estimator's own and match independent glm/sandwich", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  att <- getExportedValue("cdid", "att_gt")
  for (rt in nm_routes) {
    for (clustered in c(FALSE, TRUE)) {
      args <- nm_att_args(nm_fx, rt$em, rt$sm, rt$sp, control = list(keep_nuisance = TRUE),
        clustervars = if (clustered) "cl" else NULL, bstrap = FALSE)
      f <- do.call(att, args)
      comp <- f$chained$nuisance_models$components$selection
      expect_identical(comp$status, "available")
      expect_identical(comp$covariance$status, "available")
      expect_identical(comp$covariance$type, if (clustered) "cluster" else "unit")
      expect_identical(comp$selection_model, rt$sm)
      expect_identical(comp$selection_pooling, rt$sp)
      expect_identical(nrow(comp$cells), as.integer(f$chained$selection$n_active_cells))
      expect_identical(nrow(comp$blocks), as.integer(f$chained$selection$n_selection_blocks))
      recs <- f$chained$nuisance_diagnostics$records$selection
      for (m in seq_len(nrow(comp$blocks))) {
        id <- comp$blocks$block[m]
        cm <- comp$coefficients[comp$coefficients$block == id, , drop = FALSE]
        # Exactly the backend's fitted block (the record the default fit already reports).
        if (length(recs) == nrow(comp$blocks) && !is.null(recs[[m]]$coefficients))
          expect_identical(cm$estimate, unname(as.double(recs[[m]]$coefficients)))
        ref <- nm_reference_block(nm_fx, comp$cells[comp$cells$block == id, , drop = FALSE])
        expect_equal(cm$estimate, ref$coef, tolerance = 1e-6)
        expect_equal(cm$std_error, if (clustered) ref$se_cluster else ref$se_unit, tolerance = 1e-6)
      }
    }
  }
})

test_that("installed att_gt: legacy separate pair logits are the actual fit records with independent SEs", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  att <- getExportedValue("cdid", "att_gt")
  nm <- getExportedValue("cdid", "nuisance_models")
  for (rt in nm_legacy_routes) {
    for (clustered in c(FALSE, TRUE)) {
      f <- do.call(att, nm_att_args(nm_fx, rt$em, rt$sm, rt$sp, control = list(keep_nuisance = TRUE),
        clustervars = if (clustered) "cl" else NULL, bstrap = FALSE))
      comp <- nm(f, "selection")
      expect_identical(comp$status, "available")
      expect_match(comp$route, paste0(rt$class, " (legacy per-(pair, cohort) separate pair-logit consumer"), fixed = TRUE)
      expect_identical(c(comp$selection_model, comp$selection_pooling), c(rt$sm, rt$sp))
      expect_identical(comp$covariance$status, "available")
      expect_identical(comp$covariance$type, if (clustered) "cluster" else "unit")
      expect_match(comp$covariance$source, "stored factor H_chol", fixed = TRUE)
      expect_true(all(comp$blocks$n_cells == 1L))
      expect_true(all(comp$blocks$n_comparisons >= 1L))
      expect_identical(nrow(comp$comparisons), length(f$chained$nuisance_diagnostics$links))
      # One block per fitted (pair, cohort) record of the default diagnostics, same order and values.
      recs <- f$chained$nuisance_diagnostics$records$selection
      expect_identical(length(recs), nrow(comp$blocks))
      for (m in seq_len(nrow(comp$blocks))) {
        id <- comp$blocks$block[m]
        cm <- comp$coefficients[comp$coefficients$block == id, , drop = FALSE]
        cell <- comp$cells[comp$cells$block == id, , drop = FALSE]
        expect_identical(c(cell$s, cell$t), as.double(c(recs[[m]]$spec$pair$s, recs[[m]]$spec$pair$t)))
        expect_identical(cm$estimate, unname(as.double(recs[[m]]$coefficients)))
        ref <- nm_reference_block(nm_fx, cell)
        expect_equal(cm$estimate, ref$coef, tolerance = 1e-6)
        expect_equal(cm$std_error, if (clustered) ref$se_cluster else ref$se_unit, tolerance = 1e-6)
      }
    }
  }
})

test_that("legacy report charges the model view while building the tables", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  budget <- nm_ns("cp04_budget")
  att <- getExportedValue("cdid", "att_gt")
  nm <- getExportedValue("cdid", "nuisance_models")
  args <- nm_att_args(nm_fx, "reg", "pair_logit", "separate", bstrap = FALSE)
  base <- do.call(att, args)
  view <- NULL
  testthat::local_mocked_bindings(.package = "cdid",
    cp04_budget = function(ctrl, retained, scratch, work, output, label) {
      if (identical(label, "nuisance selection view")) view <<- list(retained = retained, work = work)
      if (identical(label, "nuisance selection tables")) {
        expect_false(is.null(view))
        expect_identical(retained[seq_along(view$retained)], view$retained)
        expect_identical(work[seq_along(view$work)], view$work)
        ctrl$max_retained_entries <- sum(retained) - retained[2L]
      }
      budget(ctrl, retained, scratch, work, output, label)
    })
  args$control <- list(keep_nuisance = TRUE)
  fit <- do.call(att, args)
  expect_identical(fit$att, base$att)
  expect_identical(nm(fit, "selection")$code, "budget")
})

test_that("installed att_gt: legacy, no-selection and ordinary routes report typed statuses", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  att <- getExportedValue("cdid", "att_gt")
  nm <- getExportedValue("cdid", "nuisance_models")
  legacy <- do.call(att, nm_att_args(nm_fx, "ipw", "pair_logit", "separate", control = list(keep_nuisance = TRUE),
    bstrap = FALSE))
  expect_identical(nm(legacy, "selection")$status, "available")
  expect_identical(nm(legacy)$components$status, c("available", "not_applicable", "not_captured", "not_applicable"))
  none <- do.call(att, list(yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~ x + region,
    data = nm_fx$data, est_method = "ipw", bstrap = FALSE, control = list(keep_nuisance = TRUE)))
  expect_identical(nm(none, "selection")$status, "not_applicable")
  expect_error(do.call(att, nm_att_args(nm_fx, "ipw", "pair_logit", "within_cohort",
    control = list(keep_nuisance = "yes"))), class = "cdid_chained_public_control_error")
  ord <- tryCatch(suppressWarnings(do.call(att, list(yname = "y", tname = "t", idname = "id", gname = "g",
    data = nm_fx$data, mode = "ordinary", allow_unbalanced_panel = TRUE, bstrap = FALSE))), error = function(e) NULL)
  skip_if_not(!is.null(ord), "did::att_gt() could not fit the tiny fixture in ordinary mode")
  expect_identical(nm(ord)$code, "ordinary_delegation")
})


test_that("one-cluster coefficient inference is unavailable without affecting point fits", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  att <- getExportedValue("cdid", "att_gt")
  nm <- getExportedValue("cdid", "nuisance_models")
  fx <- nm_fx
  fx$data$one_cluster <- 1L
  fx$frame$one_cluster <- 1L
  for (cl in c("one_cluster", "g")) {
    args <- nm_att_args(fx, "ipw", "independent_period_logit", "within_cohort",
                        clustervars = cl, bstrap = FALSE)
    args$compute_inffunc <- FALSE
    base <- do.call(att, args)
    args$control <- list(keep_nuisance = TRUE)
    fit <- do.call(att, args)
    expect_identical(fit$att, base$att)
    sel <- nm(fit, "selection")
    expect_identical(sel$status, "available")
    expect_identical(sel$covariance$code, "insufficient_clusters")
    expect_true(all(is.na(sel$coefficients$std_error)))
    if (cl == "g") expect_true(all(sel$blocks$n_fitting_clusters == 1L))
  }
})

test_that("optional workspace refusals leave the estimate and fitted cluster metadata truthful", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  budget <- nm_ns("cp04_budget")
  att <- getExportedValue("cdid", "att_gt")
  nm <- getExportedValue("cdid", "nuisance_models")
  args <- nm_att_args(nm_fx, "ipw", "independent_period_logit", "within_cohort",
                      clustervars = "cl", bstrap = FALSE)
  base <- do.call(att, args)
  refused <- "nuisance selection tables"
  testthat::local_mocked_bindings(.package = "cdid",
    cp04_budget = function(ctrl, retained, scratch, work, output, label) {
      if (identical(label, refused)) {
        stop(structure(list(message = "Diagnostic workspace refused.", call = NULL),
                       class = c("cdid_ipw_owner_bound_error", "error", "condition")))
      }
      budget(ctrl, retained, scratch, work, output, label)
    })
  args$control <- list(keep_nuisance = TRUE)
  fit <- do.call(att, args)
  expect_identical(fit$att, base$att)
  status <- nm(fit, "selection")
  expect_identical(status$status, "unavailable")
  expect_identical(status$code, "budget")
  expect_null(status$coefficients)

  refused <- "nuisance selection capture"
  fit <- do.call(att, args)
  expect_identical(fit$att, base$att)
  component <- nm(fit, "selection")
  expect_identical(component$status, "available")
  expect_identical(component$covariance$code, "budget")
  expect_true(is.na(component$covariance$n_groups))
  expect_true(all(is.na(component$coefficients$std_error)))
})

# ---------------------------------------------------------------------------------
# Augmented observation correction (selection_est_method = "dr"), DR+DR.
# ---------------------------------------------------------------------------------
nm_fx_aug <- nm_fixture(240L)

test_that("installed att_gt: DR+DR captures the augmented observation logits without changing the fit", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  att <- getExportedValue("cdid", "att_gt")
  agg <- getExportedValue("cdid", "aggte")
  nm <- getExportedValue("cdid", "nuisance_models")
  for (cfg in list(c("pair_logit", "separate"), c("independent_period_logit", "within_cohort"))) {
    args <- nm_att_args(nm_fx_aug, "dr", cfg[1], cfg[2], selection_est_method = "dr")
    set.seed(20261115L)
    f0 <- do.call(att, c(args, list(control = list())))
    a0 <- agg(f0, type = "dynamic"); s0 <- .Random.seed
    set.seed(20261115L)
    f1 <- do.call(att, c(args, list(control = list(keep_nuisance = TRUE))))
    a1 <- agg(f1, type = "dynamic"); s1 <- .Random.seed
    expect_identical(s1, s0)
    for (nm_f in c("att", "se", "inffunc", "V_analytical", "c", "n", "group", "t", "W", "Wpval"))
      expect_identical(f1[[nm_f]], f0[[nm_f]])
    expect_identical(a1$overall.att, a0$overall.att)
    expect_identical(a1$overall.se, a0$overall.se)
    expect_identical(a1$att.egt, a0$att.egt)
    expect_identical(f1$chained$selection, f0$chained$selection)
    expect_identical(f1$chained$nuisance_diagnostics, f0$chained$nuisance_diagnostics)
    expect_identical(f1$chained$public, f0$chained$public)
    expect_null(f0$chained$nuisance_models)
    expect_identical(nm(f0)$code, "not_requested")

    # DR+DR on the static route: every role is reported (tested below).
    idx <- nm(f1)
    expect_identical(idx$selection_est_method, "dr")
    expect_identical(idx$components$status, rep("available", 4L))
    expect_true(any(grepl("Z_c = j_c + S_e / max(q_c, eps) * (Delta_e - j_c)", nm_print(idx), fixed = TRUE)))
    comp <- nm(f1, "selection")
    expect_identical(c(comp$selection_model, comp$selection_pooling), cfg)
    expect_match(comp$route, "cdid_augmented_backend$selection_backend (shared cdid_selection_models backend",
      fixed = TRUE)
    expect_match(comp$role, "augmented observation correction", fixed = TRUE)
    expect_match(comp$weights$note, "1 / max(q_c, eps)", fixed = TRUE)
    expect_identical(comp$universe, "F")
    expect_identical(comp$slope_terms, c("x", "regionB"))
    expect_identical(comp$template, if (cfg[1] == "pair_logit")
      "logit Pr(S_i(s,t) = 1 | X_i, G_i = g, U_i = 1) = a[g,(s,t)] + X_i'b[g,(s,t)]" else
      "logit Pr(R_it = 1 | X_i, G_i = g) = a[g,t] + X_i'b[g]")
    expect_identical(is.null(comp$derived), cfg[1] == "pair_logit")
    expect_true(any(grepl("use: q_c", nm_print(nm(f1, block = 1L)), fixed = TRUE)))

    # Same pair/date logit estimators as before: the independent dense glm and
    # unit/cluster sandwich on the actual active cells of every block.
    fc <- do.call(att, utils::modifyList(args, list(control = list(keep_nuisance = TRUE), clustervars = "cl",
      bstrap = FALSE)))
    for (cc in list(comp, nm(fc, "selection"))) {
      clustered <- identical(cc$covariance$type, "cluster")
      expect_identical(cc$covariance$status, "available")
      for (m in seq_len(nrow(cc$blocks))) {
        id <- cc$blocks$block[m]
        cm <- cc$coefficients[cc$coefficients$block == id, , drop = FALSE]
        ref <- nm_reference_block(nm_fx_aug, cc$cells[cc$cells$block == id, , drop = FALSE])
        expect_equal(cm$estimate, ref$coef, tolerance = 1e-6)
        expect_equal(cm$std_error, if (clustered) ref$se_cluster else ref$se_unit, tolerance = 1e-6)
      }
    }
    expect_identical(nm(fc, "selection")$covariance$type, "cluster")
  }
})

# ---------------------------------------------------------------------------------
# Augmented DR+DR: outcome regression on Z, treatment propensity and imputation.
# ---------------------------------------------------------------------------------
nm_aug_roles <- c("outcome", "treatment", "imputation")

nm_aug_slopes <- function(cm) cm$estimate[cm$role %in% c("shared_slope", "slope")]

test_that("installed att_gt: DR+DR roles are the backend's fits with IDs, finite SEs and unchanged inference", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  att <- getExportedValue("cdid", "att_gt")
  agg <- getExportedValue("cdid", "aggte")
  nm <- getExportedValue("cdid", "nuisance_models")
  for (op in c("pooled", "separate")) for (ip in c("separate", "within_cohort")) {
    args <- nm_att_args(nm_fx_aug, "dr", "pair_logit", "within_cohort", selection_est_method = "dr",
      outcome_pooling = op, imputation_pooling = ip)
    set.seed(20261116L)
    f0 <- do.call(att, c(args, list(control = list())))
    a0 <- agg(f0, type = "group"); s0 <- .Random.seed
    set.seed(20261116L)
    f1 <- do.call(att, c(args, list(control = list(keep_nuisance = TRUE))))
    a1 <- agg(f1, type = "group"); s1 <- .Random.seed
    expect_identical(s1, s0)
    for (nm_f in c("att", "se", "inffunc", "V_analytical", "c", "W", "Wpval"))
      expect_identical(f1[[nm_f]], f0[[nm_f]])
    expect_identical(a1$overall.att, a0$overall.att)
    expect_identical(a1$overall.se, a0$overall.se)
    expect_identical(f1$chained$nuisance_diagnostics, f0$chained$nuisance_diagnostics)

    dg <- f1$chained$nuisance_diagnostics
    expect_false(is.null(dg$coefficients))
    comps <- lapply(stats::setNames(nm_aug_roles, nm_aug_roles), function(r) nm(f1, r))
    for (r in nm_aug_roles) {
      cc <- comps[[r]]
      expect_identical(cc$status, "available")
      expect_identical(cc$covariance$status, "available")
      expect_identical(cc$covariance$type, "unit")
      # A nontrivial deterministic example: finite positive SEs in every fit.
      expect_true(all(is.finite(cc$coefficients$std_error) & cc$coefficients$std_error > 0))
      expect_true(all(cc$blocks$n_comparisons >= 1L))
      expect_identical(sort(unique(cc$coefficients$block)), sort(cc$blocks$block))
      expect_true(all(nzchar(cc$blocks$equation)))
      # The equation is printed immediately above each coefficient table.
      out <- nm_print(nm(f1, r, block = 1L))
      e <- which(out == cc$blocks$equation[1L])
      expect_length(e, 1L)
      expect_match(out[e + 1L], "Estimate", fixed = TRUE)
    }
    # Imputation: one fit per pooling group; exactly the backend's group slopes
    # and cell intercepts (the default diagnostic records).
    imp <- comps$imputation
    expect_identical(nrow(imp$blocks), length(dg$groups))
    expect_identical(imp$cells$n_observed_pairs, as.double(dg$cells$n_observed_pairs))
    for (m in seq_len(nrow(imp$blocks))) {
      cm <- imp$coefficients[imp$coefficients$block == imp$blocks$block[m], , drop = FALSE]
      expect_identical(nm_aug_slopes(cm), unname(as.double(dg$coefficients$imputation_slopes[[imp$blocks$backend_id[m]]])))
      ii <- cm$role == "local_intercept"
      expect_identical(cm$estimate[ii], unname(as.double(dg$coefficients$imputation_local_intercepts[cm$cell[ii]])))
    }
    if (ip == "separate") expect_identical(nrow(imp$blocks), nrow(imp$cells))
    # Outcome: one pooled fit (common slope once, block intercepts listed) or one per control block.
    out_c <- comps$outcome
    if (op == "pooled") {
      expect_identical(nrow(out_c$blocks), 1L)
      expect_identical(nm_aug_slopes(out_c$coefficients), unname(as.double(dg$coefficients$shared_slopes)))
      expect_identical(sum(out_c$coefficients$role == "shared_slope"), length(dg$coefficients$shared_slopes))
    } else {
      expect_identical(nrow(out_c$blocks), nrow(out_c$control_blocks))
      for (m in seq_len(nrow(out_c$blocks)))
        expect_identical(nm_aug_slopes(out_c$coefficients[out_c$coefficients$block == out_c$blocks$block[m], ]),
          unname(as.double(dg$coefficients$block_slopes[[out_c$blocks$backend_id[m]]])))
    }
    ii <- out_c$coefficients$role == "local_intercept"
    expect_identical(out_c$coefficients$estimate[ii], unname(as.double(dg$coefficients$local_intercepts[out_c$coefficients$cell[ii]])))
    # Treatment: one fit per distinct (target, risk set); every reusing spec has the same coefficients.
    trt <- comps$treatment
    expect_identical(nrow(trt$models), length(dg$coefficients$treatment))
    for (s in seq_len(nrow(trt$models))) {
      cm <- trt$coefficients[trt$coefficients$block == trt$models$block[s], , drop = FALSE]
      expect_identical(cm$estimate, unname(as.double(dg$coefficients$treatment[[s]])))
    }
    expect_identical(sum(!trt$models$reuses_fit), nrow(trt$blocks))
  }
})

test_that("installed att_gt: DR+DR role SEs are clustered, refuse one cluster and degrade under budget refusals", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  budget <- nm_ns("cp04_budget")
  att <- getExportedValue("cdid", "att_gt")
  nm <- getExportedValue("cdid", "nuisance_models")
  args <- nm_att_args(nm_fx_aug, "dr", "pair_logit", "within_cohort", selection_est_method = "dr",
    outcome_pooling = "pooled", bstrap = FALSE, control = list(keep_nuisance = TRUE))
  fu <- do.call(att, args)
  fc <- do.call(att, utils::modifyList(args, list(clustervars = "cl")))
  for (r in nm_aug_roles) {
    u <- nm(fu, r); cl <- nm(fc, r)
    expect_identical(cl$covariance$type, "cluster")
    expect_identical(cl$coefficients$estimate, u$coefficients$estimate)
    expect_true(all(is.finite(cl$coefficients$std_error)))
    expect_true(max(abs(cl$coefficients$std_error / u$coefficients$std_error - 1)) > 1e-3)
    expect_true(all(cl$blocks$n_contributing_clusters >= 2L))
    expect_match(paste(nm_print(nm(fc, r, block = 1L)), collapse = "\n"),
      paste0(cl$blocks$n_contributing_clusters[1L], " with nonzero unit contributions"), fixed = TRUE)
  }
  fx1 <- nm_fx_aug; fx1$data$one_cluster <- 1L
  point_only <- utils::modifyList(args, list(data = fx1$data, clustervars = "one_cluster", compute_inffunc = FALSE))
  one <- do.call(att, point_only)
  expect_identical(one$att, do.call(att, utils::modifyList(point_only, list(control = list())))$att)
  for (r in nm_aug_roles) {
    x <- nm(one, r)
    expect_identical(x$status, "available")
    expect_identical(x$covariance$code, "insufficient_clusters")
    expect_true(all(is.na(x$coefficients$std_error)))
  }
  base <- do.call(att, utils::modifyList(args, list(control = list())))
  refused <- "nuisance augmented coefficient replay"
  observation_work <- NULL
  testthat::local_mocked_bindings(.package = "cdid",
    cp04_budget = function(ctrl, retained, scratch, work, output, label) {
      if (identical(label, "nuisance selection capture")) observation_work <<- sum(work[-1L])
      if (identical(label, "nuisance augmented role tables")) {
        expect_identical(work[2L], observation_work)
        expect_gt(output[1L], 0)
      }
      if (identical(label, refused))
        stop(structure(list(message = "Diagnostic workspace refused.", call = NULL),
          class = c("cdid_ipw_owner_bound_error", "error", "condition")))
      budget(ctrl, retained, scratch, work, output, label)
    })
  fit <- do.call(att, args)
  expect_identical(fit$att, base$att)
  expect_identical(fit$inffunc, base$inffunc)
  expect_identical(nm(fit, "selection")$covariance$status, "available")
  for (r in nm_aug_roles) {
    x <- nm(fit, r)
    expect_identical(x$status, "available")
    expect_identical(x$covariance$code, "budget")
    expect_identical(x$coefficients$estimate, nm(fu, r)$coefficients$estimate)
    expect_true(all(is.na(x$coefficients$std_error)))
  }
  refused <- "nuisance augmented role tables"
  fit <- do.call(att, args)
  expect_identical(fit$att, base$att)
  for (r in nm_aug_roles) expect_identical(nm(fit, r)$code, "budget")
})

test_that("a named predictor is still a slope in a no-intercept DR+DR report", {
  skip_if_not(nm_has_installed_cdid(), "installed cdid with nuisance_models() required")
  att <- getExportedValue("cdid", "att_gt")
  nm <- getExportedValue("cdid", "nuisance_models")
  fx <- nm_fx_aug
  fx$data[["(Intercept)"]] <- fx$data$x
  fit <- do.call(att, nm_att_args(fx, "dr", "pair_logit", "separate",
    selection_est_method = "dr", xformla = ~ 0 + `(Intercept)`, selection_xformla = ~ x,
    bstrap = FALSE, control = list(keep_nuisance = TRUE)))
  treatment <- nm(fit, "treatment")
  expect_identical(treatment$status, "available")
  expect_identical(treatment$slope_terms, "(Intercept)")
  expect_null(treatment$intercept_term)
  block <- nm(fit, "treatment", block = 1L)
  expect_identical(block$table$role, "slope")
  expect_match(block$equation, "c0*(Intercept)", fixed = TRUE)
  expect_true(all(is.finite(block$coefficients[, "Std. Error"])))
})

test_that("augmented roles outside DR+DR keep typed statuses", {
  other <- nm_ns("chained_nuisance_other_components")
  gate <- list(selection = "conditional_mtar", selection_est_method = "dr", is_custom_est = FALSE)
  expect_identical(other(gate, "ipw")$outcome$status, "not_applicable")
  expect_identical(other(gate, "reg")$treatment$status, "not_applicable")
  expect_identical(other(gate, "reg")$imputation$code, "outside_initial_slice")
  roles <- nm_ns("chained_nuisance_augmented_roles")
  st <- nm_ns("chained_nuisance_status")("unavailable", "reference_recipe_route", "layout")
  r <- roles(list(), NULL, list(), NULL, st)
  expect_identical(vapply(r, `[[`, "", "code"), c(outcome = "reference_recipe_route", treatment = "reference_recipe_route",
    imputation = "reference_recipe_route"))
  bad <- roles(list(), structure(list(est_method = "reg"), class = "cdid_augmented_backend"), list(), NULL, NULL)
  expect_identical(bad$outcome$code, "unrecognized_backend")
})
