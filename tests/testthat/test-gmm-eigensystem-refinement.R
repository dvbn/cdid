test_that("a checked QR/Rayleigh refinement preserves the original residual bound", {
  V <- diag(c(2, 1))
  gram <- assemble_gmm_gram_matrix(2, "optimal")
  strategy <- gmm_covariance_strategy()
  f <- diagnose_gmm_transformed_covariance
  seam <- new.env(parent = environment(f))
  seam$eigen <- function(...) {
    z <- base::eigen(...)
    z$vectors[, 1] <- z$vectors[, 1] * (1 + 1e-8)
    z
  }
  environment(f) <- seam
  result <- f(V, gram, list(), strategy)
  expect_identical(result$eigensystem_refinement$method, "QR_Rayleigh_reorthogonalization")
  expect_gt(result$eigensystem_refinement$original_residual, result$b_bound)
  expect_lte(result$r_eig, result$b_bound)
  expect_identical(result$B, V)
  expect_equal(result$eig$values, c(2, 1), tolerance = 1e-14)
  expect_equal(crossprod(result$eig$vectors), diag(2), tolerance = 1e-14)
  expect_false(result$triggered)
  expect_error(f(diag(c(2, -1)), gram, list(), strategy),
               "indefinite", class = "cdid_gmm_pd_error")
  seam$eigen <- function(...) {
    z <- base::eigen(...); z$vectors[, 1] <- 0; z
  }
  expect_error(f(V, gram, list(), strategy), "residual", class = "cdid_gmm_pd_error")
  seam$eigen <- function(...) {
    z <- base::eigen(...); z$vectors[2, 1] <- 0.1; z
  }
  expect_error(f(V, gram, list(), strategy), "residual", class = "cdid_gmm_pd_error")
})

test_that("the well-conditioned campaign covariance passes without changing B", {
  B <- readRDS(test_path("fixtures", "gmm_orthogonality_ADJ8000_376.rds"))
  raw <- eigen(B, symmetric = TRUE)
  residual <- norm(B - raw$vectors %*% diag(raw$values) %*% t(raw$vectors), "F")
  result <- diagnose_gmm_transformed_covariance(B,
    assemble_gmm_gram_matrix(nrow(B), "optimal"), list(), gmm_covariance_strategy())
  expect_identical(result$B, B)
  expect_lte(result$r_eig, result$b_bound)
  expect_gt(result$d_min, 7)
  expect_lt(result$d_max, 41)
  expect_false(result$triggered)
  expect_equal(result$b_bound, 100 * nrow(B) * .Machine$double.eps * norm(B, "F"), tolerance = 1e-14)
  if (residual > result$b_bound) {
    expect_identical(result$eigensystem_refinement$method, "QR_Rayleigh_reorthogonalization")
    expect_lt(result$r_eig, residual / 10)
  }
})
