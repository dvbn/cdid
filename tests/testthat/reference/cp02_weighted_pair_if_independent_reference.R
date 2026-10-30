# Independent numerical derivative of the FULL estimating-equation stack.
# Run from the copied workspace with:
# R_LIBS_USER=/opt/tools/r-libs Rscript --vanilla \
#   plan/evidence/cp02_checks/cp02_weighted_pair_if_independent_reference.R
# All tolerances/steps/Newton controls below are diagnostic fixture choices ONLY.
# No estimator package, population identification, MC, CLT or inference claim.
options(digits = 12, warn = 2)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
out_dir <- dirname(sub("^--file=", "", script_arg))
root_path <- file.path(out_dir, "weighted_pair_ipw_derivative_reference.R")

# The supplied seed/data construction is copied solely to align the comparisons.
# No supplied fit, IF, analytic Jacobian or perturbation is used by the reference
# equations, numerical Jacobian, Newton solve or their case-weight refits.
set.seed(141026)
n <- 240L
x <- rnorm(n)
Z <- cbind(intercept = 1, x = x)
G <- sample(c(3, 4, 5, Inf), n, replace = TRUE,
            prob = c(.32, .26, .15, .27))
K <- c(4, Inf)
cohorts <- c(3, K)
risk <- G %in% cohorts
D <- as.numeric(G == 3)
C <- as.numeric(G %in% K)
patterns <- list(c(1, 2, 3), c(1, 2), c(1, 3), c(2, 3),
                 c(1, 4), c(3, 4))
patt <- sample(seq_along(patterns), n, replace = TRUE)
observed <- matrix(FALSE, n, 4)
for (i in seq_len(n)) observed[i, patterns[[patt[i]]]] <- TRUE
alpha <- rnorm(n)
pre <- ifelse(G == 3, 2, ifelse(G == 4, -1, ifelse(G == 5, .5, 1)))
Y <- cbind(alpha, alpha + pre + .3*x + rnorm(n),
           alpha + pre + 1 + .8*x + rnorm(n),
           alpha + pre + 2 + .5*x + rnorm(n))
Y[G == 3, 3] <- Y[G == 3, 3] + 1 + .4*x[G == 3]
stopifnot(all(rowSums(observed) >= 2))

I <- vapply(cohorts, function(k) as.numeric(G == k), numeric(n))
colnames(I) <- c("g", "k4", "kinf")
m <- ncol(I)
p <- ncol(Z)
ir <- seq_len(m)
ib <- m + seq_len(p)
it <- m + p + 1L
ic <- m + p + 2L
param_names <- c("r_g", "r_4", "r_inf", "beta_0", "beta_x", "T", "C")
q <- colMeans(I)
gi <- match(G, cohorts)
pairs <- list(c(1, 2), c(1, 3))
names(pairs) <- c("12", "13")
fixture_tol <- 1e-6
solver_tol <- 2e-12
jac_steps <- c(1e-5, 5e-6)
case_steps <- c(1e-4, 5e-5)

pair_data <- function(ell) {
  S <- as.numeric(observed[, ell[1]] & observed[, ell[2]])
  delta <- numeric(n)
  delta[S == 1] <- Y[S == 1, ell[2]] - Y[S == 1, ell[1]]
  list(S = S, delta = delta)
}

# Rows: I_k(S-r_k), inverse-rate weighted logit score, and the TWO
# unnormalized residual equations for normalized means. Averaging uses n,
# the entire Q input universe, not an observed-pair or risk-set sample size.
unit_equations <- function(theta, dat) {
  stopifnot(length(theta) == length(param_names), all(is.finite(theta)),
            all(theta[ir] > 0), all(theta[ir] < 1))
  rr <- rep(1, n)
  rr[risk] <- theta[ir][gi[risk]]
  eta <- drop(Z %*% theta[ib])
  prob <- plogis(eta)
  w <- as.numeric(risk) * dat$S / rr
  a <- D * dat$S / rr
  b <- C * dat$S * exp(eta) / rr
  ans <- cbind(I * (dat$S - rep(theta[ir], each = n)),
               Z * as.numeric(w * (D - prob)),
               a * (dat$delta - theta[it]),
               b * (dat$delta - theta[ic]))
  colnames(ans) <- param_names
  stopifnot(all(is.finite(ans)))
  ans
}
equations <- function(theta, dat, case = rep(1, n)) {
  colSums(unit_equations(theta, dat) * case) / n
}

# Black-box coordinate differences of the complete moment vector.
# This function contains NO analytic derivative blocks.
numerical_jacobian <- function(theta, dat, case = rep(1, n), step = jac_steps[1]) {
  ans <- vapply(seq_along(theta), function(j) {
    h <- step * max(1, abs(theta[j]))
    plus <- minus <- theta
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    (equations(plus, dat, case) - equations(minus, dat, case)) / (2*h)
  }, numeric(length(theta)))
  dimnames(ans) <- list(param_names, param_names)
  ans
}

# Local fixture-only full-stack Newton solution, including rate and mean rows.
solve_stack <- function(dat, case = rep(1, n), start = NULL) {
  if (is.null(start)) start <- c(rep(.5, m), rep(0, p), 0, 0)
  theta <- start
  for (iter in 0:30) {
    val <- equations(theta, dat, case)
    if (max(abs(val)) <= solver_tol)
      return(list(theta = theta, iterations = iter, residual = max(abs(val))))
    if (iter == 30) stop("Fixture Newton iteration limit reached")
    change <- solve(numerical_jacobian(theta, dat, case), -val)
    accepted <- FALSE
    for (bt in 0:24) {
      trial <- theta + change * 2^(-bt)
      trial_val <- tryCatch(equations(trial, dat, case),
                            error = function(e) rep(Inf, length(theta)))
      if (all(is.finite(trial_val)) && max(abs(trial_val)) < max(abs(val))) {
        theta <- trial
        accepted <- TRUE
        break
      }
    }
    if (!accepted) stop("Fixture Newton backtracking failed")
  }
}
solve_if <- function(J, psi) t(solve(J, -t(psi)))
link <- function(full_if) full_if[, it] - full_if[, ic]

# A separate analytic comparator; NEVER called by the numerical reference.
analytic_comparator <- function(theta, dat) {
  r <- theta[ir]
  rr <- rep(1, n)
  rr[risk] <- r[gi[risk]]
  prob <- plogis(drop(Z %*% theta[ib]))
  odds <- exp(drop(Z %*% theta[ib]))
  w <- as.numeric(risk)*dat$S/rr
  a <- D*dat$S/rr
  b <- C*dat$S*odds/rr
  dp <- mean(a)
  dm <- mean(b)
  score <- Z * as.numeric(w*(D-prob))
  H <- crossprod(Z, Z * as.numeric(w*prob*(1-prob))) / n
  L <- vapply(ir, function(j)
    -colMeans(Z * as.numeric(I[, j]*dat$S*(D-prob)/r[j]^2)), numeric(p))
  rate_if <- I * (dat$S - rep(r, each = n)) / rep(q, each = n)
  beta_rate_if <- t(solve(H, t(rate_if %*% t(L))))
  beta_score_if <- t(solve(H, t(score)))
  beta_if <- beta_score_if + beta_rate_if
  gamma <- c(0, vapply(ir[-1], function(j)
    mean(I[, j]*dat$S*odds*(dat$delta-theta[ic]))/(r[j]^2*dm), numeric(1)))
  kappa <- -colMeans(Z * as.numeric(b*(dat$delta-theta[ic]))) / dm
  treated_if <- a*(dat$delta-theta[it])/dp
  control_direct <- -b*(dat$delta-theta[ic])/dm
  rate_direct <- drop(rate_if %*% gamma)
  beta_total <- drop(beta_if %*% kappa)
  components <- cbind(rate = rate_direct + drop(beta_rate_if %*% kappa),
                      beta_score = drop(beta_score_if %*% kappa),
                      treated = treated_if, control = control_direct)
  phi <- rowSums(components)
  control_if <- -control_direct - rate_direct - beta_total
  theta_if <- cbind(rate_if, beta_if, treated_if, control_if)
  colnames(theta_if) <- param_names
  J <- matrix(0, length(theta), length(theta),
              dimnames = list(param_names, param_names))
  J[ir, ir] <- -diag(q)
  J[ib, ir] <- L
  J[ib, ib] <- -H
  J[it, ir[1]] <- -mean(D*dat$S*(dat$delta-theta[it]))/r[1]^2
  J[it, it] <- -dp
  J[ic, ir] <- -dm*gamma
  J[ic, ib] <- -dm*kappa
  J[ic, ic] <- -dm
  list(phi = phi, theta_if = theta_if, J = J, components = components,
       L = L, gamma = gamma, kappa = kappa, dp = dp, dm = dm)
}

# Full parameter refits, also checking that rate reconstruction already
# carries empirical cohort composition. Normalize by sum(case) for the q law.
finite_refits <- function(theta, dat, h) {
  refit <- matrix(NA_real_, n, length(theta),
                  dimnames = list(NULL, param_names))
  q_derivative_error <- mass_identity_error <- max_residual <- 0
  masses <- function(fit, case) {
    weighted <- vapply(ir, function(j)
      sum(case*I[, j]*dat$S/fit[ir[j]])/sum(case), numeric(1))
    raw <- colSums(I*case)/sum(case)
    mass_identity_error <<- max(mass_identity_error, abs(weighted-raw))
    weighted
  }
  for (i in seq_len(n)) {
    plus <- minus <- rep(1, n)
    plus[i] <- 1+h
    minus[i] <- 1-h
    tp <- solve_stack(dat, plus, theta)
    tm <- solve_stack(dat, minus, theta)
    max_residual <- max(max_residual, tp$residual, tm$residual)
    refit[i, ] <- n*(tp$theta-tm$theta)/(2*h)
    q_if <- n*(masses(tp$theta, plus)-masses(tm$theta, minus))/(2*h)
    q_derivative_error <- max(q_derivative_error, abs(q_if-(I[i, ]-q)))
  }
  list(full_if = refit, mass_error = mass_identity_error,
       q_error = q_derivative_error, residual = max_residual)
}

cat("Independent weighted-pair IPW stacked-equation audit\n")
cat("R version:", R.version.string, "\n")
cat("n =", n, "; cohort Q counts:", colSums(I),
    "; outside-risk count:", sum(!risk), "\n")
cat("theta = (r_g, r_4, r_inf, beta_0, beta_x, T, C)\n")
cat("Black-box central Jacobian steps:", jac_steps,
    "; case steps:", case_steps, "\n")
cat("Fixture-only Newton residual target:", solver_tol,
    "; comparison threshold:", fixture_tol, "\n")
results <- list()
summary_rows <- list()
unit_rows <- list()
for (pair_name in names(pairs)) {
  dat <- pair_data(pairs[[pair_name]])
  base <- solve_stack(dat)
  theta <- base$theta
  names(theta) <- param_names
  J1 <- numerical_jacobian(theta, dat, step = jac_steps[1])
  J2 <- numerical_jacobian(theta, dat, step = jac_steps[2])
  psi <- unit_equations(theta, dat)
  numeric_if <- solve_if(J2, psi)
  phi <- link(numeric_if)
  analytic <- analytic_comparator(theta, dat)
  f1 <- finite_refits(theta, dat, case_steps[1])
  f2 <- finite_refits(theta, dat, case_steps[2])
  sources <- list(rate = ir, beta_score = ib, treated = it, control = ic)
  numeric_components <- vapply(sources, function(index) {
    forcing <- psi * 0
    forcing[, index] <- psi[, index]
    link(solve_if(J2, forcing))
  }, numeric(n))

  # Substantive mutations applied to the independent NUMERICAL Jacobian.
  J_omit <- J2
  J_omit[ib, ir] <- 0
  omitted_phi <- link(solve_if(J_omit, psi))
  J_omit_target <- J2
  J_omit_target[ib, ir[1]] <- 0
  omitted_target_phi <- link(solve_if(J_omit_target, psi))
  J_reverse_beta <- J2
  J_reverse_beta[ic, ib] <- -J_reverse_beta[ic, ib]
  reversed_phi <- link(solve_if(J_reverse_beta, psi))

  errors <- c(J_vs_analytic = max(abs(J2-analytic$J)),
    J_step_resolution = max(abs(J2-J1)),
    full_IF_vs_analytic = max(abs(numeric_if-analytic$theta_if)),
    link_IF_vs_analytic = max(abs(phi-analytic$phi)),
    IF_source_components = max(abs(numeric_components-analytic$components)),
    full_IF_vs_stack_refits = max(abs(numeric_if-f2$full_if)),
    link_IF_vs_stack_refits = max(abs(phi-link(f2$full_if))),
    refit_step_resolution = max(abs(f1$full_if-f2$full_if)),
    normalized_mass_identity = max(f1$mass_error, f2$mass_error),
    empirical_q_IF_identity = max(f1$q_error, f2$q_error))
  mutation_errors <- c(omit_all_score_rate_blocks = max(abs(omitted_phi-phi)),
    omit_target_score_rate_block = max(abs(omitted_target_phi-phi)),
    reverse_control_beta_block = max(abs(reversed_phi-phi)))
  stopifnot(all(errors < fixture_tol), all(mutation_errors > 1e-4),
            max(abs(phi[!risk])) == 0,
            max(abs(numeric_if[!risk, ])) == 0,
            max(abs(colMeans(numeric_if))) < 1e-9,
            any(abs(phi[risk & dat$S == 0]) > 1e-4),
            abs(J2[it, ir[1]]) < 1e-9,
            max(abs(J2[ib, ir[1]])) > 1e-3)

  cat("\nPAIR", pair_name, "\n")
  print(theta)
  cat("Base iterations:", base$iterations, "; residual:", base$residual,
      "; maximum refit residual:", max(f1$residual, f2$residual), "\n")
  cat("Numerical FULL Jacobian (moment rows, parameter columns):\n")
  print(J2)
  cat("Comparison errors (maximum absolute across all units/parameters):\n")
  print(errors)
  cat("Full parameter IF errors vs analytic:\n")
  print(apply(abs(numeric_if-analytic$theta_if), 2, max))
  cat("REJECTED mutations: link IF discrepancies vs independent full stack:\n")
  print(mutation_errors)
  cat("Treated outcome rate row:", J2[it, ir[1]],
      "; treated rate inside score:", J2[ib, ir[1]], "\n")
  cat("Risk pair non-observers:", sum(risk & dat$S == 0),
      "; max non-observer link IF:", max(abs(phi[risk & dat$S == 0])),
      "; outside-risk full IF:", max(abs(numeric_if[!risk, ])),
      "; mean link IF:", mean(phi), "\n")
  write.csv(J2, file.path(out_dir, paste0("cp02_weighted_pair_jacobian_", pair_name, ".csv")))
  results[[pair_name]] <- list(theta = theta, numeric_if = numeric_if, phi = phi,
                              finite_if = f2$full_if, omitted_phi = omitted_phi)
  summary_rows[[pair_name]] <- data.frame(pair = pair_name,
    check = c(names(errors), names(mutation_errors)),
    max_abs_error = c(errors, mutation_errors),
    status = c(rep("PASS", length(errors)), rep("EXPECTED_REJECTION", length(mutation_errors))))
  unit_rows[[pair_name]] <- data.frame(pair = pair_name, unit = seq_len(n),
    cohort = G, S = dat$S, risk = risk, phi_numeric = phi,
    phi_analytic = analytic$phi, phi_stack_refit = link(f2$full_if),
    phi_omitted_score_rates = omitted_phi, numeric_components)
}

# The supplied reference is loaded ONLY AFTER all independent calculations.
# It is an external comparator, with its output captured verbatim below.
root <- new.env(parent = globalenv())
root_output <- capture.output(sys.source(root_path, envir = root))
stopifnot(identical(G, root$G), identical(Z, root$Z), identical(Y, root$Y),
          identical(observed, root$observed))
numeric_phi <- do.call(cbind, lapply(results, `[[`, "phi"))
numeric_finite <- do.call(cbind, lapply(results, function(z) link(z$finite_if)))
numeric_omitted <- do.call(cbind, lapply(results, `[[`, "omitted_phi"))
root_wrong <- do.call(cbind, lapply(root$fits, `[[`, "wrong"))
root_errors <- c(numeric_IF_vs_root_analytic = max(abs(numeric_phi-root$phi)),
                numeric_IF_vs_root_refits = max(abs(numeric_phi-root$d2)),
                stack_refits_vs_root_refits = max(abs(numeric_finite-root$d2)),
                omitted_block_vs_root_mutation = max(abs(numeric_omitted-root_wrong)),
                bridge_difference_vs_root_refits = max(abs(
                  (numeric_phi[, 2]-numeric_phi[, 1])-(root$d2[, 2]-root$d2[, 1]))))
stopifnot(all(root_errors < fixture_tol),
          max(abs(root_wrong-numeric_phi)) > 1e-4)
cat("\nSUPPLIED ROOT COMPARISON (identical seed/data verified)\n")
print(root_errors)
cat("Root omitted-rate mutation error vs independent IF:",
    max(abs(root_wrong-numeric_phi)), "\n")
cat("Captured supplied reference output:\n", paste(root_output, collapse = "\n"), "\n", sep = "")
write.csv(do.call(rbind, summary_rows),
          file.path(out_dir, "cp02_weighted_pair_independent_comparisons.csv"), row.names = FALSE)
write.csv(do.call(rbind, unit_rows),
          file.path(out_dir, "cp02_weighted_pair_independent_unit_checks.csv"), row.names = FALSE)
cat("\nACCEPT: full stacked derivative, analytic terms, finite refits, n-scaling,\n")
cat("pair non-observer contributions and outside-risk zeros; mutations rejected.\n")
cat("Finite derivative evidence only. No population-assumption, consistency, MC,\n")
cat("CLT, inference, weighted/cluster/anticipation or production-solver certification.\n")
