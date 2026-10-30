# CP02 continuous heavy-tail development DGP (inst/development; not package
# API). Loads the actual inst files. Expected values are derived here
# (closed forms, quadrature, root finding, dt-based moments) or copied from
# the independent Python reference
# plan/evidence/cp02_checks/cp02_continuous_dgp_reference.py; none are taken
# from generator output. Fixed seeds are small generator checks only.

cp02_dev_dir <- function() {
  dir <- Sys.getenv("CDID_CP02_DEV_DIR", "")
  if (nzchar(dir)) return(dir)
  cands <- c(test_path("..", "..", "inst", "development"),
             system.file("development", package = "cdid"))
  for (d in cands) {
    if (nzchar(d) && file.exists(file.path(d, "cp02_continuous_dgp.R"))) {
      return(d)
    }
  }
  stop("inst/development/cp02_continuous_dgp.R not found.", call. = FALSE)
}

E <- local({
  dir <- cp02_dev_dir()
  env <- new.env(parent = baseenv())
  sys.source(file.path(dir, "seed.R"), envir = env)
  sys.source(file.path(dir, "cp02_continuous_dgp.R"), envir = env)
  env
})

# Independent definitions -------------------------------------------------
L2 <- log(2)
MX <- 5 / 3 - 1 / L2
MX_PY <- 0.22397162577770325931       # Python reference, 20 digits
f3 <- function(x) ifelse(abs(x) <= 1, 2 * L2 / 3 * 2^x, 0)
fk <- function(x) ifelse(abs(x) <= 1, 1 / 2, 0)
PATS <- list("1-2" = 1:2, "1-3" = c(1, 3), "1-4" = c(1, 4), "1-5" = c(1, 5),
             "2-3" = 2:3, "2-4" = c(2, 4), "2-5" = c(2, 5), "3-4" = 3:4,
             "3-5" = c(3, 5), "4-5" = 4:5, "3" = 3, "none" = integer(0))
PROB <- c(rep(0.08, 10), 0.10, 0.10)
DATA_COLS <- c("id", "period", "G", "x", "z1", "z2", "z3", "pattern",
               "pattern_type", "observed", "y")

# Returns a character vector of violated identities (empty = all hold).
check_draw <- function(d, tol = 1e-10) {
  bad <- character(0)
  flag <- function(ok, msg) if (!isTRUE(ok)) bad <<- c(bad, msg)
  flag(identical(names(d), c("data", "oracle", "settings")), "parts")
  p <- d$data$panel
  u <- d$data$units
  o <- d$oracle$panel
  n <- d$settings$n
  flag(identical(names(p), DATA_COLS), "data columns (oracle leak?)")
  flag(!any(c("y0", "y_g", "effect", "eps", "t5", "shift", "unit_effect") %in%
              c(names(p), names(u), names(d$data$counts))), "oracle columns")
  flag(identical(sort(u$id), seq_len(n)), "unit ids")
  flag(nrow(p) == 5 * n && nrow(o) == 5 * n, "panel rows")
  key_p <- paste(p$id, p$period)
  flag(!anyDuplicated(key_p) && setequal(key_p, paste(rep(seq_len(n), each = 5),
                                                     1:5)), "full calendar")
  ui <- match(p$id, u$id)
  for (v in c("G", "x", "z1", "z2", "z3", "pattern", "pattern_type")) {
    flag(identical(p[[v]], u[[v]][ui]), paste("unit constant", v))
  }
  flag(all(u$G %in% c(0, 3, 4, 5)), "cohort codes")
  flag(all(abs(c(u$x, u$z1, u$z2, u$z3)) < 1), "covariate support")
  flag(all(u$pattern %in% names(PATS)), "pattern labels")
  nd <- unname(lengths(PATS[u$pattern]))
  typ <- ifelse(nd == 2, "pair", ifelse(nd == 1, "singleton", "none"))
  flag(identical(u$pattern_type, typ), "pattern type")
  flag(isTRUE(all.equal(u$n_obs, as.numeric(lengths(PATS[u$pattern])))),
       "n_obs")
  want_obs <- mapply(function(pt, t) t %in% PATS[[pt]], p$pattern, p$period,
                     USE.NAMES = FALSE)
  flag(identical(p$observed, want_obs), "observed vs pattern")
  flag(identical(is.na(p$y), !p$observed), "mask")
  flag(max(tapply(p$observed, p$id, sum)) <= 2, "fully observed unit")
  # Outcome identities, joined by key.
  oi <- match(key_p, paste(o$id, o$period))
  t <- p$period
  x <- p$x
  G <- p$G
  alpha <- c(0, 0.5, -0.5, 1)[match(G, c(0, 3, 4, 5))]
  sh <- (t == 1) * ifelse(G == 3, 2 + 0.7 * x, ifelse(G == 4, -1 + 0.2 * x,
                          ifelse(G == 5, 1 - 0.4 * x, 0)))
  eps <- (1 + x / 4) * sqrt(0.6) * o$t5[oi]
  ue <- d$oracle$units$unit_effect[match(p$id, d$oracle$units$id)]
  y0 <- alpha + t / 4 + 0.4 * sin(pi * x) * (t - 1) +
    0.2 * x^2 * (t^2 - 1) + 0.1 * p$z2^2 * (t - 1) + p$z1 + ue + eps + sh
  th <- c(1, 2, 1.5)[match(G, c(3, 4, 5))]
  eff <- ifelse(G > 0 & t >= G, th + 0.5 * (t - G) + x, 0)
  flag(max(abs(o$eps[oi] - eps)) < tol, "eps identity")
  flag(max(abs(o$shift[oi] - sh)) < tol, "shift identity")
  flag(max(abs(o$y0[oi] - y0)) < tol, "y0 identity")
  flag(max(abs(o$effect[oi] - eff)) < tol, "effect identity")
  flag(max(abs(o$y_g[oi] - y0 - eff)) < tol, "y_g identity")
  ob <- which(p$observed)
  flag(max(abs(p$y[ob] - o$y_g[oi][ob])) < tol, "observed y = oracle Y(G)")
  cnt <- table(factor(u$G, c(0, 3, 4, 5)), factor(u$pattern, names(PATS)))
  ci <- cbind(match(d$data$counts$G, c(0, 3, 4, 5)),
              match(d$data$counts$pattern, names(PATS)))
  flag(nrow(d$data$counts) == 48 && sum(d$data$counts$units) == n &&
         all(d$data$counts$units == cnt[ci]), "counts")
  bad
}

check_truth <- function(tr, tol = 1e-12) {
  bad <- character(0)
  flag <- function(ok, msg) if (!isTRUE(ok)) bad <<- c(bad, msg)
  flag(abs(tr$cohort$mean_x[tr$cohort$G == 3] - MX_PY) < tol, "mean x")
  flag(all(tr$cohort$mean_x[tr$cohort$G != 3] == 0), "mean x others")
  flag(all(tr$cohort$mass == 0.25), "masses")
  a <- tr$att_gt
  want <- data.frame(group = c(3, 3, 3, 4, 4, 5), t = c(3, 4, 5, 4, 5, 5),
                     att = c(1 + MX, 1.5 + MX, 2 + MX, 2, 2.5, 1.5))
  flag(nrow(a) == 6 && identical(paste(a$group, a$t), paste(want$group, want$t)) &&
         max(abs(a$att - want$att)) < tol, "post ATT")
  flag(nrow(tr$pre) == 9 && all(tr$pre$causal_att == 0) &&
         all(is.na(tr$pre$reporting_did_truth)), "pre cells")
  # Hand-derived aggregate rationals: (const, coefficient on mean x).
  ref <- rbind(
    c("simple", "overall", 7 / 4, 1 / 2), c("group", "g=3", 3 / 2, 1),
    c("group", "g=4", 9 / 4, 0), c("group", "g=5", 3 / 2, 0),
    c("group", "overall", 7 / 4, 1 / 3), c("dynamic", "e=0", 3 / 2, 1 / 3),
    c("dynamic", "e=1", 2, 1 / 2), c("dynamic", "e=2", 2, 1),
    c("dynamic", "overall", 11 / 6, 11 / 18), c("calendar", "t=3", 1, 1),
    c("calendar", "t=4", 7 / 4, 1 / 2), c("calendar", "t=5", 2, 1 / 3),
    c("calendar", "overall", 19 / 12, 11 / 18))
  ag <- tr$aggregates
  k <- match(paste(ref[, 1], ref[, 2]), paste(ag$aggregate, ag$component))
  flag(nrow(ag) == 13 && !anyNA(k) &&
         max(abs(ag$att[k] - (as.numeric(ref[, 3]) +
                                as.numeric(ref[, 4]) * MX))) < tol,
       "aggregates")
  w <- tr$aggregate_weights
  flag(all(w$weight > 0) && all(w$t >= w$group) &&
         max(abs(tapply(w$weight, paste(w$aggregate, w$component), sum) - 1)) <
           tol, "aggregate weights")
  attw <- a$att[match(paste(w$group, w$t), paste(a$group, a$t))]
  re <- tapply(w$weight * attw, paste(w$aggregate, w$component), sum)
  flag(max(abs(re[paste(ag$aggregate, ag$component)] - ag$att)) < tol,
       "aggregate = weights x ATT")
  bad
}

d_small <- E$cp02_continuous_draw(60, 11)
d_mid <- E$cp02_continuous_draw(2400, 20261015)
d_big <- E$cp02_continuous_draw(24000, 7)
truth <- E$cp02_continuous_truth()

# Analytic pieces -----------------------------------------------------------

test_that("G3 density, CDF, inverse and mean are independently confirmed", {
  expect_equal(integrate(f3, -1, 1, rel.tol = 1e-12)$value, 1, tolerance = 1e-10)
  m <- integrate(function(x) x * f3(x), -1, 1, rel.tol = 1e-12)$value
  expect_equal(m, MX_PY, tolerance = 1e-10)
  expect_equal(MX, MX_PY, tolerance = 1e-14)
  v <- c(1e-6, 0.01, 0.2, 0.5, 0.77, 0.999999)
  xq <- E$cp02_continuous_x_quantile(v, TRUE)
  root <- vapply(v, function(p) uniroot(function(x) integrate(
    f3, -1, x, rel.tol = 1e-12)$value - p, c(-1, 1), tol = 1e-13)$root,
    numeric(1))
  expect_equal(xq, root, tolerance = 1e-8)
  expect_equal(E$cp02_continuous_x_quantile(0.5, TRUE), 0.3219280948873623,
               tolerance = 1e-14)  # Python reference median
  expect_equal(E$cp02_continuous_x_quantile(v, FALSE), 2 * v - 1)
  ratio <- f3(c(-1, 0, 1)) / fk(c(-1, 0, 1))
  expect_equal(ratio, 4 * L2 / 3 * 2^c(-1, 0, 1))
  expect_equal(c(truth$density_ratio$lower, truth$density_ratio$upper),
               c(0.4620981203732969, 1.8483924814931875), tolerance = 1e-14)
})

test_that("cohort and pattern maps give exact slot probabilities", {
  mid4 <- (0:3 + 0.5) / 4
  expect_identical(E$cp02_continuous_cohort(mid4), c(3, 4, 5, 0))
  mid <- (0:49 + 0.5) / 50
  idx <- E$cp02_continuous_pattern_index(mid)
  expect_identical(as.vector(table(factor(idx, 1:12))) / 50, PROB)
  edge <- c(0, 4, 8, 36, 40, 45) / 50 + 1e-12
  expect_identical(E$cp02_continuous_pattern_index(edge), c(1, 2, 3, 10, 11, 12))
  s <- E$cp02_continuous_settings()$patterns
  expect_identical(s$pattern, names(PATS))
  expect_equal(s$prob, PROB)
  expect_equal(s$slots / 50, PROB)
})

test_that("truth: post ATT, aggregates, weights match independent derivations", {
  expect_identical(check_truth(truth), character(0))
  py <- c(simple = 1.86198581288885162965, group = 1.82465720859256775310,
          dynamic = 1.97020488241970754735, calendar = 1.72020488241970754735)
  ov <- truth$aggregates[truth$aggregates$component == "overall", ]
  expect_equal(ov$att[match(names(py), ov$aggregate)], unname(py),
               tolerance = 1e-14)
  # Truth does not depend on any realized sample.
  expect_identical(d_small$oracle$truth, truth)
  expect_identical(d_big$oracle$truth, truth)
})

test_that("pre-period: causal effect 0, period-1 DiD reporting not zero", {
  # E[Y2(0) - Y1(0) | G = g] - same for NT = -E[s_g(X) | G = g] (X ~ U(-1,1)).
  s4 <- integrate(function(x) (-1 + 0.2 * x) * fk(x), -1, 1)$value
  s5 <- integrate(function(x) (1 - 0.4 * x) * fk(x), -1, 1)$value
  expect_equal(c(-s4, -s5), c(1, -1))
  expect_true(all(truth$pre$causal_att == 0))
  expect_true(all(is.na(truth$pre$reporting_did_truth)))
})

test_that("propensity for declared future-control operators is the stated logit", {
  pr <- truth$propensity
  expect_identical(nrow(pr), 12L)
  for (i in seq_len(nrow(pr))) {
    g <- pr$group[i]
    tt <- pr$t[i]
    pool <- if (pr$control[i] == "NT") 0 else c(c(3, 4, 5)[c(3, 4, 5) > tt], 0)
    expect_identical(pr$pool[i], paste(pool, collapse = ","))
    expect_false(3 %in% pool)
    xs <- seq(-1, 1, by = 0.05)
    num <- if (g == 3) f3(xs) else fk(xs)     # equal masses cancel
    den <- num + length(pool) * fk(xs)
    lp <- qlogis(num / den)
    expect_equal(lp, pr$intercept[i] + pr$slope_x[i] * xs, tolerance = 1e-12)
  }
  expect_true(all(pr$slope_z == 0 & pr$correct_x_only & pr$correct_x_plus_z))
  py <- c(-0.07883084812988339957, -0.77197802868982870899,
          -1.17744313679799309097)
  g3 <- pr[pr$group == 3, ]
  expect_equal(g3$intercept, py[g3$pool_size], tolerance = 1e-14)
  # Not claimed: target 4 vs earlier-treated pool {3, 0} is not affine.
  lg <- function(x) qlogis(fk(x) / (fk(x) + f3(x) + fk(x)))
  expect_equal(lg(0.5) - 2 * lg(0) + lg(-0.5), -0.02983353592308929588,
               tolerance = 1e-12)
  expect_match(truth$propensity_note, "not covered")
})

test_that("selection: pair MTAR, P(A) = .8, Q(pair|G) = .1, N/50", {
  sel <- truth$selection
  expect_equal(unname(sel$pattern_prob), PROB)
  expect_identical(names(sel$pattern_prob), names(PATS))
  expect_equal(sel$p_pair, 0.8)
  expect_equal(sel$q_pair_given_g, 0.1)
  expect_equal(c(2400, 9600, 38400) / 4 * 0.08, c(48, 192, 768))
  expect_false(any(lengths(PATS) > 2))
})

test_that("standardized t5 noise has variance 1 and fourth moment 9", {
  m2 <- 2 * integrate(function(x) x^2 * dt(x, 5), 0, Inf, rel.tol = 1e-12)$value
  m4 <- 2 * integrate(function(x) x^4 * dt(x, 5), 0, Inf, rel.tol = 1e-12)$value
  expect_equal(c(m2, m4), c(5 / 3, 25), tolerance = 1e-8)
  expect_equal(c(m2, m4) * c(0.6, 0.36), c(1, 9), tolerance = 1e-8)
  expect_equal(c(truth$noise$var_standardized_t5,
                 truth$noise$fourth_standardized_t5), c(1, 9))
  expect_equal(E$cp02_continuous_settings()$noise$scale^2, 0.6)
})

# Generator ----------------------------------------------------------------

test_that("draws satisfy every outcome, mask and metadata identity", {
  for (d in list(d_small, d_mid, d_big)) {
    expect_identical(check_draw(d), character(0))
  }
  # Key-based checks are unaffected by permuting rows.
  d <- d_mid
  pp <- order((seq_len(nrow(d$data$panel)) * 7919) %% 12007)
  d$data$panel <- d$data$panel[pp, ]
  d$oracle$panel <- d$oracle$panel[rev(seq_len(nrow(d$oracle$panel))), ]
  d$data$units <- d$data$units[rev(seq_len(nrow(d$data$units))), ]
  expect_identical(check_draw(d), character(0))
})

test_that("unit frame keeps singleton/none units; nothing fully observed", {
  u <- d_mid$data$units
  expect_identical(u$id, seq_len(2400L))
  expect_true(all(c("pair", "singleton", "none") %in% u$pattern_type))
  p <- d_mid$data$panel
  none_ids <- u$id[u$pattern == "none"]
  expect_true(all(is.na(p$y[p$id %in% none_ids])))
  expect_identical(sum(p$id %in% none_ids), 5L * length(none_ids))
  expect_true(all(p$observed[p$id %in% u$id[u$pattern == "3"]] == (p$period[
    p$id %in% u$id[u$pattern == "3"]] == 3)))
})

test_that("small n leaves empty cells unforced; counts are realized only", {
  cn <- d_small$data$counts
  expect_identical(sum(cn$units), 60L)
  expect_true(any(cn$units == 0L))
  expect_identical(nrow(cn), 48L)
  expect_identical(d_small$settings$n, 60L)
})

test_that("fixed-seed distribution checks (generator only)", {
  u <- d_big$data$units
  F3 <- function(x) (2^(x + 1) - 1) / 3      # integral of f3
  expect_gt(ks.test(u$x[u$G == 3], F3)$p.value, 1e-3)
  expect_gt(ks.test(u$x[u$G != 3], "punif", -1, 1)$p.value, 1e-3)
  for (z in c("z1", "z2", "z3")) {
    expect_gt(ks.test(u[[z]][u$G == 3], "punif", -1, 1)$p.value, 1e-3)
  }
  expect_gt(chisq.test(table(u$G))$p.value, 1e-3)
  expect_gt(chisq.test(table(factor(u$pattern, names(PATS))), p = PROB)$p.value,
            1e-3)
  expect_gt(chisq.test(table(u$G, u$pattern))$p.value, 1e-3)
  expect_gt(t.test(u$x[u$pattern_type == "pair"],
                   u$x[u$pattern_type != "pair"])$p.value, 1e-3)
  o <- d_big$oracle$panel
  expect_gt(ks.test(o$t5, "pt", 5)$p.value, 1e-3)
})

test_that("inputs are validated before allocation", {
  for (bad_n in list(0, -1, 2.5, NA, NA_real_, Inf, c(10, 20), "10", TRUE,
                     numeric(0))) {
    expect_error(E$cp02_continuous_draw(bad_n, 1), "`n`")
  }
  expect_error(E$cp02_continuous_draw(1e12, 1), "too large")
  expect_error(E$cp02_continuous_draw(.Machine$integer.max %/% 5 + 1, 1),
               "too large")
  for (bad_s in list(1.5, NA, Inf, "1", c(1, 2), 2^31, numeric(0))) {
    expect_error(E$cp02_continuous_draw(10, bad_s), "`seed`")
  }
  expect_error(E$cp02_continuous_draw(1e12, 1.5), "`n`")
  expect_identical(nrow(E$cp02_continuous_draw(1, 3)$data$panel), 5L)
})

test_that("seed reproducibility and caller RNG state restoration", {
  expect_identical(E$cp02_continuous_draw(500, 42), E$cp02_continuous_draw(500, 42))
  expect_false(identical(E$cp02_continuous_draw(500, 42)$data,
                         E$cp02_continuous_draw(500, 43)$data))
  genv <- globalenv()
  run <- function() {
    old_kind <- RNGkind()
    had <- exists(".Random.seed", envir = genv, inherits = FALSE)
    if (had) old <- get(".Random.seed", envir = genv)
    on.exit({
      rng_formals <- names(formals(RNGkind))
      n_kind <- min(length(old_kind), length(rng_formals))
      do.call(RNGkind, stats::setNames(as.list(old_kind[seq_len(n_kind)]), rng_formals[seq_len(n_kind)]))
      if (had) {
        assign(".Random.seed", old, envir = genv)
      } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
        rm(".Random.seed", envir = genv)
      }
    })
    RNGkind("Mersenne-Twister", "Box-Muller", "Rejection")
    expected_kind <- RNGkind()
    expect_identical(expected_kind[1:3], c("Mersenne-Twister", "Box-Muller", "Rejection"))
    set.seed(99)
    before <- get(".Random.seed", envir = genv)
    E$cp02_continuous_draw(50, 5)
    expect_identical(get(".Random.seed", envir = genv), before)
    expect_identical(RNGkind(), expected_kind)
    rm(".Random.seed", envir = genv)
    E$cp02_continuous_draw(50, 5)
    expect_false(exists(".Random.seed", envir = genv, inherits = FALSE))
    expect_identical(RNGkind(), expected_kind)
  }
  run()
})

test_that("documented draw order replays exactly (same R RNG only)", {
  n <- 2400
  replay <- local({
    genv <- globalenv()
    old_kind <- RNGkind()
    had <- exists(".Random.seed", envir = genv, inherits = FALSE)
    if (had) old <- get(".Random.seed", envir = genv)
    on.exit({
      RNGkind(old_kind[1], old_kind[2], old_kind[3])
      if (had) {
        assign(".Random.seed", old, envir = genv)
      } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
        rm(".Random.seed", envir = genv)
      }
    })
    set.seed(20261015, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
             sample.kind = "Rejection")
    r <- list(vg = runif(n), vx = runif(n), z1 = runif(n, -1, 1),
              z2 = runif(n, -1, 1), z3 = runif(n, -1, 1), ue = rnorm(n),
              num = rnorm(5 * n))
    r$chi <- rowSums(sapply(1:5, function(k) rnorm(5 * n)^2))
    r$vp <- runif(n)
    r
  })
  u <- d_mid$data$units
  G <- c(3, 4, 5, 0)[findInterval(replay$vg, c(0.25, 0.5, 0.75)) + 1]
  expect_identical(u$G, G)
  x <- ifelse(G == 3, log2((3 * replay$vx + 1) / 2), 2 * replay$vx - 1)
  expect_equal(u$x, x, tolerance = 1e-14)
  expect_identical(list(u$z1, u$z2, u$z3), list(replay$z1, replay$z2, replay$z3))
  expect_identical(d_mid$oracle$units$unit_effect, replay$ue)
  expect_equal(d_mid$oracle$panel$t5, replay$num / sqrt(replay$chi / 5),
               tolerance = 1e-13)
  cum <- cumsum(PROB)[-12]
  pt <- names(PATS)[findInterval(replay$vp, cum) + 1]
  expect_identical(u$pattern, pt)
  expect_match(paste(E$cp02_continuous_settings()$rng$note), "No")
})

# Mutation controls: the checkers must reject corrupted draws/truth ----------

test_that("checkers reject corrupted missingness, outcomes and truth", {
  d0 <- d_mid
  p <- d0$data$panel
  hid <- which(!p$observed)[1]
  mut <- list(
    leak_value = function(d) {
      d$data$panel$y[hid] <- d$oracle$panel$y_g[hid]; d },
    leak_column = function(d) {
      d$data$panel$y0 <- d$oracle$panel$y0; d },
    full_unit = function(d) {
      k <- d$data$panel$id == 1
      d$data$panel$observed[k] <- TRUE
      d$data$panel$y[k] <- d$oracle$panel$y_g[k]; d },
    relabel_pattern = function(d) {
      i <- which(d$data$units$pattern == "1-2")[1]
      d$data$units$pattern[i] <- "4-5"
      d$data$panel$pattern[d$data$panel$id == i] <- "4-5"; d },
    drop_row = function(d) {
      d$data$panel <- d$data$panel[-hid, ]; d },
    drop_none_unit = function(d) {
      i <- d$data$units$id[d$data$units$pattern == "none"][1]
      d$data$units <- d$data$units[d$data$units$id != i, ]; d },
    homoskedastic = function(d) {
      d$oracle$panel$eps <- sqrt(0.6) * d$oracle$panel$t5; d },
    wrong_shift = function(d) {
      k <- d$oracle$panel$period == 1
      d$oracle$panel$y0[k] <- d$oracle$panel$y0[k] + 0.1; d },
    wrong_effect = function(d) {
      d$oracle$panel$effect <- d$oracle$panel$effect * 1.01; d },
    stale_counts = function(d) {
      d$data$counts$units[1:2] <- d$data$counts$units[2:1] + c(1L, -1L); d }
  )
  for (nm in names(mut)) {
    expect_gt(length(check_draw(mut[[nm]](d0))), 0, label = nm)
  }
  tmut <- list(
    mean_x = function(t) { t$cohort$mean_x[1] <- 5 / 3 - 1 / log(10); t },
    pretrend_zero = function(t) { t$pre$reporting_did_truth <- 0; t },
    realized_mass = function(t) {
      w <- t$aggregate_weights
      k <- w$aggregate == "simple"
      w$weight[k] <- c(0.2, 0.2, 0.2, 0.15, 0.15, 0.1); t$aggregate_weights <- w; t },
    att = function(t) { t$att_gt$att[4] <- 2.1; t },
    aggregate = function(t) { t$aggregates$att[1] <- t$aggregates$att[1] + 1e-6; t }
  )
  for (nm in names(tmut)) {
    expect_gt(length(check_truth(tmut[[nm]](truth))), 0, label = nm)
  }
})
