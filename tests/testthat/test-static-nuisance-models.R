# Private static observation-nuisance kernels (R/static_nuisance_models.R).
# Functions come from the test/namespace environment; only the independent
# dense reference is sourced via testthat::test_path().
source(testthat::test_path("reference", "static_nuisance_reference.R"), local = TRUE)

gc_sn <- list(epsilon = 1e-12, maxit = 100, trace = FALSE)

# Deterministic fixture: three dates, overlapping role populations G, C, C2;
# date-observation O_is; reference-date covariates x_s differ by date, so the
# per-cell designs differ across cells for the same unit.
sn_fixture <- function(N = 150, seed = 20261104L, common_x = FALSE) {
  set.seed(seed)
  id <- sample(sprintf("u%04d", 1000:9999), N)
  xb <- stats::rnorm(N)
  wv <- stats::runif(N)
  x <- sapply(1:3, function(s) xb + if (common_x) 0 else 0.5 * stats::rnorm(N))
  mG <- stats::runif(N) < 0.45
  mC <- !mG
  mC2 <- stats::runif(N) < 0.55            # overlaps both G and C
  O <- sapply(1:3, function(s) stats::rbinom(N, 1, stats::plogis(0.7 + 0.8 * x[, s] - 0.6 * wv + 0.2 * s)))
  O[1:6, ] <- 0                             # zero-outcome units (kept everywhere)
  O[7:12, 1] <- 1                           # single-outcome units
  O[7:12, 2:3] <- 0
  mC2[1:12] <- TRUE
  outside <- 13:15                          # units in no mask of the C2-only fits
  mC2[outside] <- FALSE
  des <- function(s, mask) {
    d <- cbind(x = x[, s], w = wv)
    d[!mask, ] <- NA                        # never read outside the needed rows
    d
  }
  resp <- function(v, mask) { v <- as.double(v); v[!mask] <- NA; v }
  date_cells <- list()
  for (r in list(list("G", mG), list("C", mC), list("C2", mC2))) for (s in 1:3) {
    date_cells[[paste0(r[[1]], ":d", s)]] <- list(mask = r[[2]], response = resp(O[, s], r[[2]]),
                                                 design = des(s, r[[2]]))
  }
  pc <- function(mask, s, u) list(mask = mask, response = resp(O[, s] * O[, u], mask), design = des(s, mask))
  pair_cells <- list(`G:12` = pc(mG, 1, 2), `C:12` = pc(mC, 1, 2), `C:23` = pc(mC, 2, 3),
                     `C2:13` = pc(mC2, 1, 3))
  pair_cells$`Cdup:12` <- pair_cells$`C:12` # different key, identical mask/values
  y <- sapply(1:3, function(s) 1 + 0.5 * s + 0.9 * x[, s] + 0.4 * wv + stats::rnorm(N))
  ic <- function(mask, s, u) {
    obs <- mask & O[, s] == 1 & O[, u] == 1
    dy <- y[, u] - y[, s]
    dy[!obs] <- NA                          # dy exists only on A & O
    list(mask = mask, observed = obs, dy = dy, design = des(s, mask))
  }
  imp_cells <- list(`G:12` = ic(mG, 1, 2), `C:12` = ic(mC, 1, 2), `C:23` = ic(mC, 2, 3),
                    `C2:13` = ic(mC2, 1, 3))
  list(N = N, id = id, O = O, outside = outside, date_cells = date_cells,
       pair_cells = pair_cells, imp_cells = imp_cells)
}

fx <- sn_fixture()
pair_parts <- list(
  separate = list(b1 = "G:12", b2 = "C:12", b3 = "C:23", b4 = "C2:13", b5 = "Cdup:12"),
  intermediate = list(g = "G:12", c = c("C:12", "C:23", "Cdup:12"), c2 = "C2:13"),
  shared = list(all = c("G:12", "C:12", "C:23", "C2:13", "Cdup:12")))
fd_units <- c(1L, 7L, 13L, 20L, 33L, 101L)   # zero-, single-outcome, outside C2, overlapping

strip_design <- function(cells) lapply(cells, function(cl) { cl$design <- NULL; cl })

test_that("pair logits match full stacked glm for separate/intermediate/shared blocks", {
  for (nm in names(pair_parts)) {
    part <- pair_parts[[nm]]
    f <- static_nuisance_logit_fit(fx$id, fx$pair_cells, part, gc_sn)
    expect_equal(unname(f$coefficients), ref_fit(fx$pair_cells, part, "logit", rep(1, fx$N)),
                 tolerance = 1e-8, info = nm)
    nfit <- sum(vapply(unlist(part), function(k) sum(fx$pair_cells[[k]]$mask), numeric(1)))
    expect_equal(sum(vapply(f$block_fits, `[[`, numeric(1), "n_fit")), nfit)
  }
  # Distinct keys are never deduplicated even with identical masks/values.
  f <- static_nuisance_logit_fit(fx$id, fx$pair_cells, pair_parts$intermediate, gc_sn)
  expect_true(all(c("c|alpha|C:12", "c|alpha|Cdup:12") %in% names(f$coefficients)))
  expect_equal(f$coefficients[["c|alpha|C:12"]], f$coefficients[["c|alpha|Cdup:12"]], tolerance = 1e-8)
  # Intercept-only: local intercepts are the cell logits of full-mask rates.
  f0 <- static_nuisance_logit_fit(fx$id, strip_design(fx$pair_cells), pair_parts$shared, gc_sn)
  rate <- vapply(pair_parts$shared$all, function(k) {
    cl <- fx$pair_cells[[k]]; mean(cl$response[cl$mask]) }, numeric(1))
  expect_equal(unname(f0$coefficients), unname(stats::qlogis(rate)), tolerance = 1e-9)
  # Chunked provider gives the same fit.
  fc <- static_nuisance_logit_fit(fx$id, fx$pair_cells, pair_parts$shared, gc_sn, max_chunk_rows = 7)
  f1 <- static_nuisance_logit_fit(fx$id, fx$pair_cells, pair_parts$shared, gc_sn)
  expect_equal(fc$coefficients, f1$coefficients, tolerance = 1e-10)
})

test_that("logit IF: dense reference, weight-perturbation refits, chunks and retained units", {
  ids_all <- fx$id
  for (nm in names(pair_parts)) {
    part <- pair_parts[[nm]]
    f <- static_nuisance_logit_fit(fx$id, fx$pair_cells, part, gc_sn)
    IF <- static_nuisance_if(f, ids_all)
    expect_equal(unname(IF), ref_if_dense(fx$pair_cells, part, "logit", unname(f$coefficients), fx$N),
                 tolerance = 1e-8, info = nm)
    expect_lt(max(abs(colMeans(IF))), 1e-8)
    fd <- ref_if_fd(fx$pair_cells, part, "logit", fx$N, fd_units)
    expect_equal(unname(IF[fd_units, ]), fd, tolerance = 1e-5, info = nm)
  }
  f <- static_nuisance_logit_fit(fx$id, fx$pair_cells, pair_parts$shared, gc_sn)
  IF <- static_nuisance_if(f, fx$id)
  # zero-outcome units are failures with nonzero scores; units in no mask give 0
  expect_true(all(rowSums(abs(IF[1:6, ])) > 0))
  fC2 <- static_nuisance_logit_fit(fx$id, fx$pair_cells["C2:13"], list(b = "C2:13"), gc_sn)
  expect_equal(unname(static_nuisance_if(fC2, fx$id[fx$outside])), matrix(0, 3, 3))
  # chunked, reordered requested rows reproduce the full evaluation
  perm <- rev(seq_len(fx$N))
  ch <- split(perm, ceiling(seq_along(perm) / 17))
  IFc <- do.call(rbind, lapply(ch, function(r) static_nuisance_if(f, fx$id[r])))
  expect_equal(IFc, IF[perm, ], tolerance = 0)
  expect_identical(rownames(IFc), fx$id[perm])
})

test_that("date-cell logits, raw pair q and date products with keyed derivatives", {
  dc <- fx$date_cells
  parts <- list(same = list(G = c("G:d1", "G:d2", "G:d3"), C = c("C:d1", "C:d2", "C:d3"),
                            C2 = c("C2:d1", "C2:d2", "C2:d3")),
                diff = stats::setNames(as.list(names(dc)), paste0("b_", seq_along(dc))))
  rows <- fx$id[which(dc$`G:d1`$mask)][c(5, 1, 9, 3)]
  pos <- match(rows, fx$id)
  for (nm in names(parts)) {
    f <- static_nuisance_logit_fit(fx$id, dc, parts[[nm]], gc_sn)
    th <- unname(f$coefficients)
    expect_equal(th, ref_fit(dc, parts[[nm]], "logit", rep(1, fx$N)), tolerance = 1e-8)
    pr <- static_nuisance_q_predict(f, rows, dates = c("G:d1", "G:d3"))
    lin <- function(k) { bl <- f$layout$blocks[[f$layout$cell_block[[k]]]]
      th[bl$alpha_index[f$layout$cell_local[[k]]]] +
        drop(dc[[k]]$design[pos, ] %*% th[bl$beta_index]) }
    expect_equal(unname(pr$value), stats::plogis(lin("G:d1")) * stats::plogis(lin("G:d3")), tolerance = 1e-12)
    expect_identical(pr$floored, FALSE)
    expect_match(pr$assumption, "not tested")
    nslope <- sum(grepl("|beta|", pr$param_keys, fixed = TRUE))
    expect_equal(nslope, if (nm == "same") 2L else 4L)  # shared slopes appear once
    # central differences of the actual predictor in the full parameter vector
    P <- length(th)
    num <- sapply(seq_len(P), function(j) { e <- replace(numeric(P), j, 1e-6)
      (static_nuisance_q_predict(f, rows, dates = c("G:d1", "G:d3"), coefficients = th + e)$value -
         static_nuisance_q_predict(f, rows, dates = c("G:d1", "G:d3"), coefficients = th - e)$value) / 2e-6 })
    expect_equal(ref_expand(pr, P), unname(num), tolerance = 1e-7)
    p1 <- static_nuisance_q_predict(f, rows, pair = "G:d2")
    expect_equal(unname(p1$value), stats::plogis(lin("G:d2")), tolerance = 1e-12)
    expect_equal(unname(p1$deriv), unname(p1$value * (1 - p1$value) * cbind(1, dc$`G:d2`$design[pos, ])),
                 tolerance = 1e-12)
  }
  # No floor/cap on dates or product: tiny raw products pass through.
  f <- static_nuisance_logit_fit(fx$id, dc, parts$same, gc_sn)
  th <- f$coefficients
  th[c("G|alpha|G:d1", "G|alpha|G:d3")] <- -12
  low <- static_nuisance_q_predict(f, rows, dates = c("G:d1", "G:d3"), coefficients = th)
  expect_true(all(low$value < 1e-6 & low$value > 0))
  expect_true(all(low$q_dates < 1e-4))
  # Single-outcome units contribute successes/failures in date fits, zero-outcome units failures.
  expect_true(all(dc$`C2:d1`$mask[1:12]))
  expect_true(all(rowSums(abs(static_nuisance_if(f, fx$id[1:12]))) > 0))
})

test_that("imputation WLS: stacked references, NA isolation, predictions and IF", {
  ic <- fx$imp_cells
  parts <- list(separate = list(a = "G:12", b = "C:12", c = "C:23", d = "C2:13"),
                shared = list(all = c("G:12", "C:12", "C:23", "C2:13")))
  for (nm in names(parts)) {
    for (cells in list(ic, strip_design(ic))) {
      f <- static_nuisance_wls_fit(fx$id, cells, parts[[nm]], rank_tol = 1e-10)
      expect_equal(unname(f$coefficients), ref_fit(cells, parts[[nm]], "wls", rep(1, fx$N)),
                   tolerance = 1e-10)
      IF <- static_nuisance_if(f, fx$id)
      expect_equal(unname(IF), ref_if_dense(cells, parts[[nm]], "wls", unname(f$coefficients), fx$N),
                   tolerance = 1e-9)
      fd <- ref_if_fd(cells, parts[[nm]], "wls", fx$N, fd_units)
      expect_equal(unname(IF[fd_units, ]), fd, tolerance = 1e-6)
      # units without any observed pair have no imputation score
      expect_equal(unname(IF[1:12, ]), matrix(0, 12, ncol(IF)))
    }
  }
  f <- static_nuisance_wls_fit(fx$id, ic, parts$shared, rank_tol = 1e-10, max_chunk_rows = 5)
  f1 <- static_nuisance_wls_fit(fx$id, ic, parts$shared, rank_tol = 1e-10)
  expect_equal(f$coefficients, f1$coefficients, tolerance = 1e-10)
  # predictions over the full role population A (including unobserved rows)
  A_ids <- fx$id[which(ic$`C:23`$mask)]
  pj <- static_nuisance_j_predict(f, rev(A_ids), "C:23")
  pos <- match(rev(A_ids), fx$id)
  th <- f$coefficients
  expect_equal(unname(pj$value), unname(th[["all|alpha|C:23"]] +
                 drop(ic$`C:23`$design[pos, ] %*% th[c("all|beta|x", "all|beta|w")])), tolerance = 1e-12)
  expect_equal(unname(pj$deriv), unname(cbind(1, ic$`C:23`$design[pos, ])), tolerance = 0)
  expect_identical(pj$param_keys, c("all|alpha|C:23", "all|beta|x", "all|beta|w"))
  expect_true(any(!ic$`C:23`$observed[pos]))
})

test_that("mutations: common-Z shortcut and n_A normalization fail the refit check", {
  part <- list(all = c("C:12", "C:23", "C2:13"))
  ev <- function(fxx) {
    f <- static_nuisance_logit_fit(fxx$id, fxx$pair_cells[part$all], part, gc_sn)
    bf <- f$block_fits$all
    cs <- sapply(part$all, function(k) { cl <- fxx$pair_cells[[k]]; s <- numeric(fxx$N)
      r <- which(cl$mask); s[r] <- cl$response[r] - static_nuisance_q_predict(f, fxx$id[r], pair = k)$value; s })
    Z1 <- fxx$pair_cells$`C:12`$design            # one common Z per unit (shortcut)
    Z1[is.na(Z1[, 1]), ] <- fxx$pair_cells$`C2:13`$design[is.na(Z1[, 1]), ]
    Z1[is.na(Z1[, 1]), ] <- 0
    short <- evaluate_selection_block_if(bf, Z1, cs, fxx$N)$xi_theta
    list(f = f, IF = static_nuisance_if(f, fxx$id), short = short,
         fd = ref_if_fd(fxx$pair_cells, part, "logit", fxx$N, fd_units))
  }
  r <- ev(fx)
  expect_equal(unname(r$IF[fd_units, ]), r$fd, tolerance = 1e-5)
  expect_gt(max(abs(r$short[fd_units, ] - r$fd)), 1e-2)      # wrong score detected
  nA <- sum(Reduce(`|`, lapply(part$all, function(k) fx$pair_cells[[k]]$mask)))
  expect_lt(nA, fx$N)
  expect_gt(max(abs(r$IF[fd_units, ] * nA / fx$N - r$fd)), 1e-2) # wrong normalization detected
  # With one common design per unit across cells (time-invariant x), the existing
  # common-Z helper and the cell-centred kernel coincide (reuse equivalence).
  fxc <- sn_fixture(common_x = TRUE)
  rc <- ev(fxc)
  expect_equal(unname(rc$IF), unname(rc$short), tolerance = 1e-10)
})

test_that("genuine refusals and layout checks", {
  dc <- fx$date_cells
  bad <- dc["G:d1"]; bad$`G:d1`$response[bad$`G:d1`$mask] <- 1
  expect_error(static_nuisance_logit_fit(fx$id, bad, list(b = "G:d1"), gc_sn),
               class = "cdid_cp04_separation_error")
  bad$`G:d1`$response[bad$`G:d1`$mask] <- 0
  expect_error(static_nuisance_logit_fit(fx$id, bad, list(b = "G:d1"), gc_sn),
               class = "cdid_cp04_separation_error")
  two <- dc[c("G:d1", "G:d2")]; two$`G:d2`$response[two$`G:d2`$mask] <- 1
  expect_error(static_nuisance_logit_fit(fx$id, two, list(b = c("G:d1", "G:d2")), gc_sn),
               class = "cdid_cp04_separation_error")    # pooled block keeps local intercepts
  one <- dc["G:d1"]; one$`G:d1`$design <- cbind(one$`G:d1`$design, k = 1)
  expect_error(static_nuisance_logit_fit(fx$id, one, list(b = "G:d1"), gc_sn), class = "cdid_cp04_rank_error")
  nb <- dc["G:d1"]; nb$`G:d1`$response[which(nb$`G:d1`$mask)[1]] <- NA
  expect_error(static_nuisance_logit_fit(fx$id, nb, list(b = "G:d1"), gc_sn), class = "cdid_ipw_fit_error")
  ic <- fx$imp_cells["G:12"]; ic$`G:12`$observed[] <- FALSE
  expect_error(static_nuisance_wls_fit(fx$id, ic, list(b = "G:12"), 1e-10), class = "cdid_or_wls_rank_error")
  ie <- "cdid_static_nuisance_input_error"
  expect_error(static_nuisance_logit_fit(fx$id, dc[1:2], list(b = "G:d1"), gc_sn), class = ie)
  expect_error(static_nuisance_logit_fit(fx$id, dc[1:2], list(b = names(dc)[1:2], c = names(dc)[1]), gc_sn), class = ie)
  mis <- dc[1:2]; colnames(mis[[2]]$design) <- c("x", "v")
  expect_error(static_nuisance_logit_fit(fx$id, mis, list(b = names(dc)[1:2]), gc_sn), class = ie)
  mis[[2]]$design <- mis[[2]]$design[, 1, drop = FALSE]
  expect_error(static_nuisance_logit_fit(fx$id, mis, list(b = names(dc)[1:2]), gc_sn), class = ie)
  expect_error(static_nuisance_logit_fit(fx$id, dc[1:2], list(b = names(dc)[1:2]), list(epsilon = 1e-8)),
               class = "cdid_ipw_control_error")
  expect_error(static_nuisance_wls_fit(fx$id, fx$imp_cells[1], list(b = names(fx$imp_cells)[1])), "rank_tol")
  f <- static_nuisance_logit_fit(fx$id, dc, list(G = c("G:d1", "G:d2", "G:d3"),
                                                 C = c("C:d1", "C:d2", "C:d3"),
                                                 C2 = c("C2:d1", "C2:d2", "C2:d3")), gc_sn)
  outside <- fx$id[which(!dc$`G:d1`$mask)[1]]
  expect_error(static_nuisance_q_predict(f, outside, pair = "G:d1"), class = ie)
  expect_error(static_nuisance_q_predict(f, fx$id[which(dc$`G:d1`$mask)[1]], dates = c("G:d1", "G:d1")), class = ie)
  expect_error(static_nuisance_q_predict(f, fx$id[1:2], pair = "G:d1", dates = c("G:d1", "G:d2")), class = ie)
  expect_error(static_nuisance_if(f, c(fx$id[1], fx$id[1])), class = ie)
  expect_error(static_nuisance_if(f, "nope"), class = ie)
  numeric_fit <- f
  numeric_fit$unit_id <- seq_along(f$unit_id)
  expect_error(static_nuisance_if(numeric_fit, TRUE), class = ie)
})
