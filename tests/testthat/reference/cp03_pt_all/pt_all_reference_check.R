# Base-R deterministic cross-check for spec/pt_all_reference.py
# Runs with: R_LIBS_USER=/opt/tools/r-libs Rscript --vanilla spec/pt_all_reference_check.R
# Uses only base R matrix algebra (no installed package or production source calls).

# Fixture 1: Restricted GLS (theta = H * eta) with zero row in WH
W1 <- matrix(c(-1,  0, 0,
                1, -1, 0,
                0,  0, 1), nrow = 3, ncol = 3)
H1 <- matrix(c(1, 0, 0,
               0, 0, 1), nrow = 3, ncol = 2)
R1 <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
WH1 <- W1 %*% H1
V1 <- matrix(c( 4, -1,  0,
               -1,  3, -1,
                0, -1,  2), nrow = 3, ncol = 3)
m1 <- c(1.5, 0.5, 2.5)

V1_inv <- solve(V1)
cov_eta1 <- solve(t(WH1) %*% V1_inv %*% WH1)
Lambda_eta1 <- V1_inv %*% WH1 %*% cov_eta1
eta_hat1 <- as.numeric(t(Lambda_eta1) %*% m1)
theta_sub1 <- as.numeric(H1 %*% eta_hat1)

cov_unr1 <- solve(t(W1) %*% V1_inv %*% W1)
theta_unr1 <- as.numeric(cov_unr1 %*% t(W1) %*% V1_inv %*% m1)
gain1 <- cov_unr1 %*% t(R1) %*% solve(R1 %*% cov_unr1 %*% t(R1))
theta_lag1 <- as.numeric(theta_unr1 - gain1 %*% (R1 %*% theta_unr1))

stopifnot(max(abs(eta_hat1 - c(-5 / 3, 8 / 3))) < 1e-15)
stopifnot(max(abs(theta_sub1 - theta_lag1)) < 1e-15)
stopifnot(max(abs(Lambda_eta1[2, ] - c(-1 / 3, 1 / 3))) < 1e-15)

# Fixture 2: Multi-baseline single-ATT star vs chain tree with two zero rows in WH
V_star <- matrix(c(3, 1, 1,
                   1, 2, 1,
                   1, 1, 2), nrow = 3, ncol = 3)
d_hat <- c(1.4, 1.1, 1.0)
ones3 <- matrix(1, nrow = 3, ncol = 1)
V_star_inv <- solve(V_star)
w_star <- as.numeric(V_star_inv %*% ones3) / as.numeric(t(ones3) %*% V_star_inv %*% ones3)
tau_star <- sum(w_star * d_hat)
var_star <- 1 / as.numeric(t(ones3) %*% V_star_inv %*% ones3)

A_chain <- matrix(c(1,  0,  0,
                   -1,  1,  0,
                    0, -1,  1), nrow = 3, ncol = 3)
WH_chain <- A_chain %*% ones3
m_chain <- as.numeric(A_chain %*% d_hat)
V_chain <- A_chain %*% V_star %*% t(A_chain)
V_chain_inv <- solve(V_chain)
var_chain <- 1 / as.numeric(t(WH_chain) %*% V_chain_inv %*% WH_chain)
Lambda_chain <- as.numeric(V_chain_inv %*% WH_chain) * var_chain
tau_chain <- sum(Lambda_chain * m_chain)

stopifnot(abs(tau_star - 28 / 25) < 1e-15)
stopifnot(abs(tau_chain - tau_star) < 1e-15)
stopifnot(abs(var_chain - var_star) < 1e-15)
stopifnot(max(abs(w_star - c(0.2, 0.4, 0.4))) < 1e-15)
stopifnot(max(abs(Lambda_chain - c(0.2, 0.6, 1.0))) < 1e-15)

cat("Base-R cross-check passed: all exact identities verified to < 1e-15.\n")
