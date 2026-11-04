# An empirical-weight derivative refits every nuisance model. This reference
# is independent of the analytic influence-function implementation and keeps
# all probabilities away from clipping/trimming boundaries.
test_that("ordinary DR dependency differentiates pre-period outcome fitting", {
  n <- 120L
  x <- seq(-1.5, 1.5, length.out = n)
  X <- cbind(1, x)
  D <- rep(c(0, 1, 1, 0), length.out = n)
  post <- rep(c(0, 0, 1, 1, 1), length.out = n)
  y <- .2 + .7 * x + .5 * post + D * post * (1 + .3 * x) + .2 * sin(seq_len(n))
  fit <- function(w) DRDID::drdid_rc(y, post, D, X, i.weights = w, inffunc = TRUE)
  ref <- fit(rep(1, n))
  cells <- split(seq_len(n), interaction(D, post))
  selected <- unname(unlist(lapply(cells, function(i) i[c(1L, length(i))])))
  for (h in c(.001, .0005)) {
    derivative <- vapply(selected, function(i) {
      plus <- minus <- rep(1, n)
      plus[i] <- plus[i] + h
      minus[i] <- minus[i] - h
      (fit(plus)$ATT - fit(minus)$ATT) / (h / (n + h) + h / (n - h))
    }, numeric(1))
    # Central finite differences have O(h^2) error; both step sizes must pass.
    expect_lt(max(abs(derivative - as.vector(ref$att.inf.func)[selected])), 1e-6)
  }
  expect_equal(as.numeric(ref$se)^2, sum(ref$att.inf.func^2) / n^2,
               tolerance = 1e-12)
})
