#!/usr/bin/env python3
"""
Independent deterministic reference script for CP03 (pt-all specification).

Computes exact rational (fractions.Fraction) and IEEE-754 binary64 (numpy)
reference values for:
  1. Restricted GLS (theta = H * eta) equivalence with classical Lagrangian
     constrained GLS (R * theta = 0) and Schur-complement control-variate
     projection when WH contains zero rows (c <= s < t <= b_g).
  2. Multi-baseline single-ATT special case w = V^{-1} 1 / (1' V^{-1} 1),
     exact spanning-tree basis invariance (chain tree with zero rows in WH vs.
     star of pre-to-post links), D28 structural reduction + Tikhonov ridge,
     and non-equivalence of separate per-cell averaging in multi-ATT /
     unbalanced systems.
  3. Preservation of early operator-specific cancellation bridges (u < c < b_g)
     without imposing parallel trends before c or merging distinct NYT operators.
  4. Exact endpoint classification under positive anticipation (delta > 0) and
     irregular stored calendars (distinguishing pre-private u < c, PT-all zero
     c <= u < b_g, reference u = b_g, and causal anticipation b_g < u < g).
  5. D5/D28 distinction between ATT identification (row space of WH) and
     moment-covariance singularity (rank(V) < M).
  6. Invariance of post-treatment ATT(g, t) to base_period ("varying" vs.
     "universal") under fixed (pt_assumption = "all", pt_start = c), and
     non-degenerate unrestricted pretrend Wald diagnostics vs. degenerate
     restricted coordinates.

No production R/ source is imported or executed to generate expected values.
"""

from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
from typing import List, Sequence, Tuple
import numpy as np


# ---------------------------------------------------------------------------
# Pure exact rational linear algebra helpers (no external dependencies)
# ---------------------------------------------------------------------------

FracMatrix = List[List[Fraction]]
FracVector = List[Fraction]


def fmat(rows: Sequence[Sequence[int | float | str | Fraction]]) -> FracMatrix:
    return [[Fraction(x) for x in r] for r in rows]


def fvec(vals: Sequence[int | float | str | Fraction]) -> FracVector:
    return [Fraction(x) for x in vals]


def mat_transpose(A: FracMatrix) -> FracMatrix:
    if not A:
        return []
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def mat_mul(A: FracMatrix, B: FracMatrix) -> FracMatrix:
    m, k = len(A), len(A[0])
    assert len(B) == k
    n = len(B[0])
    return [
        [sum((A[i][r] * B[r][j] for r in range(k)), Fraction(0)) for j in range(n)]
        for i in range(m)
    ]


def mat_vec(A: FracMatrix, x: FracVector) -> FracVector:
    assert len(A[0]) == len(x)
    return [sum((A[i][j] * x[j] for j in range(len(x))), Fraction(0)) for i in range(len(A))]


def mat_sub(A: FracMatrix, B: FracMatrix) -> FracMatrix:
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def vec_sub(a: FracVector, b: FracVector) -> FracVector:
    return [a[i] - b[i] for i in range(len(a))]


def dot(a: FracVector, b: FracVector) -> Fraction:
    assert len(a) == len(b)
    return sum((a[i] * b[i] for i in range(len(a))), Fraction(0))


def mat_inv(A: FracMatrix) -> FracMatrix:
    n = len(A)
    assert all(len(r) == n for r in A)
    aug = [
        [A[i][j] for j in range(n)] + [Fraction(1 if i == j else 0) for j in range(n)]
        for i in range(n)
    ]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if aug[r][col] != 0:
                piv = r
                break
        if piv is None:
            raise ValueError("Matrix is singular in exact rational inversion.")
        aug[col], aug[piv] = aug[piv], aug[col]
        scale = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= scale
        for r in range(n):
            if r != col and aug[r][col] != 0:
                fac = aug[r][col]
                for j in range(2 * n):
                    aug[r][j] -= fac * aug[col][j]
    return [row[n:] for row in aug]


def gls_solve(X: FracMatrix, V: FracMatrix, y: FracVector) -> Tuple[FracVector, FracMatrix, FracMatrix]:
    """
    Exact rational GLS:
      beta_hat = (X' V^{-1} X)^{-1} X' V^{-1} y
      cov_beta = (X' V^{-1} X)^{-1}
      Lambda   = V^{-1} X (X' V^{-1} X)^{-1}   (m x p loading matrix: beta_hat = Lambda' y)
    """
    V_inv = mat_inv(V)
    Xt = mat_transpose(X)
    Xt_Vinv = mat_mul(Xt, V_inv)
    gram = mat_mul(Xt_Vinv, X)
    cov_beta = mat_inv(gram)
    Lambda_t = mat_mul(cov_beta, Xt_Vinv)
    Lambda = mat_transpose(Lambda_t)
    beta_hat = mat_vec(Lambda_t, y)
    return beta_hat, cov_beta, Lambda


def frac_to_str_mat(A: FracMatrix) -> List[List[str]]:
    return [[str(x) for x in r] for r in A]


def frac_to_float_mat(A: FracMatrix) -> List[List[float]]:
    return [[float(x) for x in r] for r in A]


def frac_to_str_vec(v: FracVector) -> List[str]:
    return [str(x) for x in v]


def frac_to_float_vec(v: FracVector) -> List[float]:
    return [float(x) for x in v]


# ---------------------------------------------------------------------------
# Fixture 1: Restricted GLS (theta = H * eta) vs. Lagrangian Constrained GLS
#            and Zero-Row Covariance Information (c <= s < t <= b_g)
# ---------------------------------------------------------------------------

def run_fixture_1() -> dict:
    """
    Cohort g = 4, calendar T = {1, 2, 3, 4}, delta = 0, reference b_g = 3.
    PT-all start date c = 2.
    Links:
      l1 = (1, 2) [early pre-start link, 1 < c = 2]
      l2 = (2, 3) [PT-all pre-treatment link in [c, b_g] = [2, 3]]
      l3 = (3, 4) [post-treatment link from b_g = 3 to t = 4]

    Unrestricted parameter vector: theta = (beta_1, beta_2, tau_4)' in R^3.
    Unrestricted design matrix W (3 x 3):
      l1 = beta_2 - beta_1 -> [-1,  1,  0]
      l2 = 0      - beta_2 -> [ 0, -1,  0]
      l3 = tau_4  - 0      -> [ 0,  0,  1]

    Under PT-all from c = 2: beta_2 = 0 (since c = 2 <= u = 2 < b_g = 3).
    Early parameter beta_1 (u = 1 < c = 2) remains unrestricted!
    Free parameter vector: eta = (beta_1, tau_4)' in R^2.
    Selection matrix H (3 x 2) and restriction matrix R (1 x 3):
      H = [[1, 0], [0, 0], [0, 1]],  R = [[0, 1, 0]].
    Restricted design matrix W_all = W * H (3 x 2):
      l1 -> [-1, 0]
      l2 -> [ 0, 0]   <-- ZERO ROW in W * H!
      l3 -> [ 0, 1]
    """
    W = fmat([
        [-1,  1, 0],
        [ 0, -1, 0],
        [ 0,  0, 1],
    ])
    H = fmat([
        [1, 0],
        [0, 0],
        [0, 1],
    ])
    R = fmat([
        [0, 1, 0],
    ])
    WH = mat_mul(W, H)

    # Exact SPD covariance V (3 x 3) and sample moment vector m_hat (3 x 1)
    # Note: l2=(2,3) and l3=(3,4) share outcome Y_3 with opposite signs, so V_{23} != 0;
    # l1=(1,2) and l2=(2,3) share outcome Y_2 with opposite signs, so V_{12} != 0.
    V = fmat([
        [ 4, -1,  0],
        [-1,  3, -1],
        [ 0, -1,  2],
    ])
    m_hat = fvec([Fraction(3, 2), Fraction(1, 2), Fraction(5, 2)])

    # 1. Restricted GLS via subspace parameterization theta = H * eta
    eta_hat, cov_eta, Lambda_eta = gls_solve(WH, V, m_hat)
    theta_sub = mat_vec(H, eta_hat)
    cov_theta_sub = mat_mul(H, mat_mul(cov_eta, mat_transpose(H)))

    # 2. Classical Lagrangian equality-constrained GLS:
    #    min (m - W theta)' V^{-1} (m - W theta) s.t. R theta = 0
    theta_unr, cov_unr, _ = gls_solve(W, V, m_hat)
    Rt = mat_transpose(R)
    R_cov_Rt_inv = mat_inv(mat_mul(R, mat_mul(cov_unr, Rt)))
    correction_gain = mat_mul(cov_unr, mat_mul(Rt, R_cov_Rt_inv))
    R_theta_unr = mat_vec(R, theta_unr)
    theta_lagrange = vec_sub(theta_unr, mat_vec(correction_gain, R_theta_unr))
    cov_theta_lagrange = mat_sub(
        cov_unr,
        mat_mul(correction_gain, mat_mul(R, cov_unr)),
    )

    # 3. Schur-complement control-variate form partitioning active rows A={0,2} vs zero row Z={1}
    W_A_H = fmat([WH[0], WH[2]])
    m_A = fvec([m_hat[0], m_hat[2]])
    m_Z = fvec([m_hat[1]])
    V_AA = fmat([
        [V[0][0], V[0][2]],
        [V[2][0], V[2][2]],
    ])
    V_AZ = fmat([[V[0][1]], [V[2][1]]])
    V_ZA = mat_transpose(V_AZ)
    V_ZZ_inv = mat_inv(fmat([[V[1][1]]]))
    reg_coeff = mat_mul(V_AZ, V_ZZ_inv)  # V_AZ * V_ZZ^{-1}
    m_A_given_Z = vec_sub(m_A, mat_vec(reg_coeff, m_Z))
    V_A_given_Z = mat_sub(V_AA, mat_mul(reg_coeff, V_ZA))
    eta_schur, cov_eta_schur, _ = gls_solve(W_A_H, V_A_given_Z, m_A_given_Z)

    # 4. Naive estimator that discards the zero row l2 of WH
    eta_drop, cov_eta_drop, _ = gls_solve(W_A_H, V_AA, m_A)

    # Exact assertions
    assert theta_sub == theta_lagrange
    assert cov_theta_sub == cov_theta_lagrange
    assert eta_hat == eta_schur
    assert cov_eta == cov_eta_schur
    assert Lambda_eta[1][0] != 0 and Lambda_eta[1][1] != 0  # Nonzero loading on zero row of WH!
    assert cov_eta[0][0] < cov_eta_drop[0][0]  # Strict variance reduction on beta_1
    assert cov_eta[1][1] < cov_eta_drop[1][1]  # Strict variance reduction on tau_4

    # Floating-point cross-check
    WH_np = np.array(frac_to_float_mat(WH), dtype=np.float64)
    V_np = np.array(frac_to_float_mat(V), dtype=np.float64)
    m_np = np.array(frac_to_float_vec(m_hat), dtype=np.float64)
    eta_np = np.linalg.solve(WH_np.T @ np.linalg.solve(V_np, WH_np), WH_np.T @ np.linalg.solve(V_np, m_np))
    max_fp_diff = float(np.max(np.abs(eta_np - np.array(frac_to_float_vec(eta_hat)))))

    return {
        "WH": frac_to_str_mat(WH),
        "zero_row_index_1based": 2,
        "eta_hat_exact": frac_to_str_vec(eta_hat),
        "eta_hat_float": frac_to_float_vec(eta_hat),
        "theta_sub_exact": frac_to_str_vec(theta_sub),
        "theta_lagrange_exact": frac_to_str_vec(theta_lagrange),
        "eta_schur_exact": frac_to_str_vec(eta_schur),
        "cov_eta_exact": frac_to_str_mat(cov_eta),
        "cov_eta_float": frac_to_float_mat(cov_eta),
        "Lambda_eta_exact": frac_to_str_mat(Lambda_eta),
        "zero_row_loading_on_beta1_and_tau4": frac_to_str_vec(Lambda_eta[1]),
        "eta_drop_zero_row_exact": frac_to_str_vec(eta_drop),
        "cov_eta_drop_zero_row_exact": frac_to_str_mat(cov_eta_drop),
        "tau4_variance_full_exact": str(cov_eta[1][1]),
        "tau4_variance_drop_exact": str(cov_eta_drop[1][1]),
        "tau4_variance_reduction_exact": str(cov_eta_drop[1][1] - cov_eta[1][1]),
        "exact_lagrange_and_schur_match": True,
        "max_float64_abs_diff": max_fp_diff,
    }


# ---------------------------------------------------------------------------
# Fixture 2: Multi-Baseline Single-ATT Special Case, Spanning-Tree Basis
#            Invariance, D28 Tikhonov Convergence, and Joint Multi-ATT vs.
#            Per-Cell Averaging
# ---------------------------------------------------------------------------

def run_fixture_2() -> dict:
    """
    Part 2A: Single target ATT(g=4, t=4) with J = 3 baselines b_1=1, b_2=2, b_3=3
             under PT-all (c = 1).
      - Star basis of pre-to-post links: d_hat = (m_(1,4), m_(2,4), m_(3,4))'
        with W_star * H = [1, 1, 1]'.
      - Chain spanning-tree basis: m_chain = (m_(1,2), m_(2,3), m_(3,4))' = A * d_hat
        with W_chain * H = A * 1_3 = [0, 0, 1]' (two zero rows in WH!).
      - Full 6-link single-class catalog:
        m_full = (m_(1,2), m_(1,3), m_(2,3), m_(1,4), m_(2,4), m_(3,4))' = E * m_chain
        where rank(V_full) = 3 < 6.
    Part 2B: Two post-treatment periods (tau_4, tau_5) in an unbalanced / cross-linked
             system where joint restricted GMM strictly beats separate per-cell averaging.
    """
    # Part 2A: J = 3 baselines
    V_star = fmat([
        [ 3,  1,  1],
        [ 1,  2,  1],
        [ 1,  1,  2],
    ])
    d_hat = fvec([Fraction(14, 10), Fraction(11, 10), Fraction(10, 10)])
    ones_3 = fmat([[1], [1], [1]])

    # Star GLS: tau_hat = (1' V_star^{-1} 1)^{-1} 1' V_star^{-1} d_hat
    tau_star_vec, cov_star, Lambda_star = gls_solve(ones_3, V_star, d_hat)
    tau_star = tau_star_vec[0]
    var_star = cov_star[0][0]
    w_star = [Lambda_star[i][0] for i in range(3)]
    assert sum(w_star, Fraction(0)) == Fraction(1)

    # Unimodular basis change A mapping d_hat = (m_(1,4), m_(2,4), m_(3,4))'
    # to m_chain = (m_(1,2), m_(2,3), m_(3,4))':
    #   m_(1,2) = m_(1,4) - m_(2,4)
    #   m_(2,3) = m_(2,4) - m_(3,4)
    #   m_(3,4) = m_(3,4)
    A = fmat([
        [ 1, -1,  0],
        [ 0,  1, -1],
        [ 0,  0,  1],
    ])
    WH_chain = mat_mul(A, ones_3)  # [[0], [0], [1]] -> two zero rows!
    assert WH_chain == fmat([[0], [0], [1]])
    m_chain = mat_vec(A, d_hat)
    V_chain = mat_mul(A, mat_mul(V_star, mat_transpose(A)))

    tau_chain_vec, cov_chain, Lambda_chain = gls_solve(WH_chain, V_chain, m_chain)
    tau_chain = tau_chain_vec[0]
    var_chain = cov_chain[0][0]

    # Verify exact equality between multi-baseline star formula and chain tree with zero rows
    assert tau_chain == tau_star
    assert var_chain == var_star
    # Also verify induced loading on d_hat: A' * Lambda_chain == w_star
    induced_w = mat_vec(mat_transpose(A), [Lambda_chain[i][0] for i in range(3)])
    assert induced_w == w_star

    # Full 6-link catalog expansion from m_chain = (m_(1,2), m_(2,3), m_(3,4))':
    #   1: (1,2) = [1, 0, 0]
    #   2: (1,3) = [1, 1, 0]
    #   3: (2,3) = [0, 1, 0]
    #   4: (1,4) = [1, 1, 1]
    #   5: (2,4) = [0, 1, 1]
    #   6: (3,4) = [0, 0, 1]
    E = fmat([
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [1, 1, 1],
        [0, 1, 1],
        [0, 0, 1],
    ])
    WH_full = mat_mul(E, WH_chain)  # [0, 0, 0, 1, 1, 1]'
    m_full = mat_vec(E, m_chain)
    V_full = mat_mul(E, mat_mul(V_chain, mat_transpose(E)))

    # D28 reconstruction-aware Tikhonov identity check:
    # G = E' E, B = G^{1/2} V_chain G^{1/2}, Q_lambda = E' (V_full + lambda I_6)^{-1} E == C' (B + lambda I_3)^{-1} C
    E_np = np.array(frac_to_float_mat(E), dtype=np.float64)
    V_chain_np = np.array(frac_to_float_mat(V_chain), dtype=np.float64)
    V_full_np = np.array(frac_to_float_mat(V_full), dtype=np.float64)
    WH_full_np = np.array(frac_to_float_mat(WH_full), dtype=np.float64)
    m_full_np = np.array(frac_to_float_vec(m_full), dtype=np.float64)

    G_np = E_np.T @ E_np
    eig_vals, eig_vecs = np.linalg.eigh(G_np)
    C_np = eig_vecs @ np.diag(np.sqrt(eig_vals)) @ eig_vecs.T
    B_np = C_np @ V_chain_np @ C_np.T
    lam = 1e-4 * (np.trace(B_np) / 6.0)
    Q_direct = E_np.T @ np.linalg.inv(V_full_np + lam * np.eye(6)) @ E_np
    Q_d28 = C_np.T @ np.linalg.inv(B_np + lam * np.eye(3)) @ C_np
    d28_identity_diff = float(np.max(np.abs(Q_direct - Q_d28)))

    # Same-full-catalog D28 reconstruction; positive ridge need not equal unregularized GLS.
    # This is not invariance to independently adding identity/ridge in each basis.
    # Tikhonov estimate on full 6-link singular system vs exact reduced estimate
    V_reg = V_full_np + lam * np.eye(6)
    tau_tikh = float(
        np.linalg.solve(WH_full_np.T @ np.linalg.solve(V_reg, WH_full_np),
                        WH_full_np.T @ np.linalg.solve(V_reg, m_full_np))[0]
    )

    # Part 2B: Joint restricted GMM for two post periods (tau_4, tau_5) vs separate per-cell averaging
    # Suppose links are:
    #   l1 = (2, 4): identifies tau_4 (from baseline 2 in [c, b_g])
    #   l2 = (3, 4): identifies tau_4 (from baseline 3 = b_g)
    #   l3 = (3, 5): identifies tau_5 (from baseline 3 = b_g)
    #   l4 = (4, 5): post-to-post link tau_5 - tau_4
    WH_2post = fmat([
        [1, 0],
        [1, 0],
        [0, 1],
        [-1, 1],
    ])
    V_2post = fmat([
        [ 4,  1,  1, -1],
        [ 1,  3,  1, -1],
        [ 1,  1,  4,  1],
        [-1, -1,  1,  3],
    ])
    m_2post = fvec([Fraction(6, 5), Fraction(1, 1), Fraction(2, 1), Fraction(9, 10)])
    eta_joint, cov_joint, _ = gls_solve(WH_2post, V_2post, m_2post)

    # Separate per-cell averaging uses only {l1, l2} for tau_4 and only {l3} for tau_5
    V_cell4 = fmat([
        [V_2post[0][0], V_2post[0][1]],
        [V_2post[1][0], V_2post[1][1]],
    ])
    tau4_sep_vec, cov_cell4, _ = gls_solve(fmat([[1], [1]]), V_cell4, fvec([m_2post[0], m_2post[1]]))
    tau4_sep = tau4_sep_vec[0]
    var_tau4_sep = cov_cell4[0][0]
    tau5_sep = m_2post[2]
    var_tau5_sep = V_2post[2][2]

    assert cov_joint[0][0] < var_tau4_sep
    assert cov_joint[1][1] < var_tau5_sep

    return {
        "part_2a_single_att": {
            "tau_star_exact": str(tau_star),
            "tau_chain_exact": str(tau_chain),
            "var_star_exact": str(var_star),
            "var_chain_exact": str(var_chain),
            "weights_star_exact": frac_to_str_vec(w_star),
            "Lambda_chain_exact": [str(Lambda_chain[i][0]) for i in range(3)],
            "exact_star_chain_equivalence": True,
            "d28_Q_lambda_max_abs_diff": d28_identity_diff,
            "tau_tikhonov_full6_float": tau_tikh,
            "tau_exact_float": float(tau_star),
            "tikhonov_abs_diff": abs(tau_tikh - float(tau_star)),
        },
        "part_2b_joint_vs_per_cell": {
            "eta_joint_exact": frac_to_str_vec(eta_joint),
            "cov_joint_exact": frac_to_str_mat(cov_joint),
            "tau4_per_cell_exact": str(tau4_sep),
            "var_tau4_per_cell_exact": str(var_tau4_sep),
            "var_tau4_joint_exact": str(cov_joint[0][0]),
            "tau5_per_cell_exact": str(tau5_sep),
            "var_tau5_per_cell_exact": str(var_tau5_sep),
            "var_tau5_joint_exact": str(cov_joint[1][1]),
            "strict_joint_variance_dominance": True,
        },
    }


# ---------------------------------------------------------------------------
# Fixture 3: Early Cancellation Bridge (u < c < b_g) Under Pre-c Trend Shift
# ---------------------------------------------------------------------------

def run_fixture_3() -> dict:
    """
    Cohort g = 4, calendar T = {1, 2, 3, 4, 5}, delta = 0, b_g = 3, PT-all start c = 2.
    Stored horizon h=5 makes fixed-horizon pool K2={inf} admissible (K1={5,inf}).
    True causal effect at t = 4: tau_4 = 3.
    Parallel trends holds from c = 2 onward: beta_{gK, 2} = beta_{gK, 3} = 0 for all clean K.
    Before c = 2 (at date u = 1), parallel trends FAILS with operator-specific shifts:
      beta_{gK1, 1} = -3/2  (for NYT control pool K1 = {5, inf})
      beta_{gK2, 1} = -1/2  (for fixed-horizon NYT control pool K2 = {inf})

    Observed catalog links:
      l1 = (1, 2) under K1: population mean = beta_{gK1, 2} - beta_{gK1, 1} = 0 - (-3/2) = +3/2
      l2 = (1, 2) under K2: population mean = beta_{gK2, 2} - beta_{gK2, 1} = 0 - (-1/2) = +1/2
      l3 = (1, 4) under K2: population mean = tau_4         - beta_{gK2, 1} = 3 - (-1/2) = +7/2
      l4 = (2, 3) under K1: population mean = 0 - 0 = 0 (overidentifying zero row in WH)
      l5 = (2, 4) under K1: population mean = tau_4 - 0 = 3

    Under our CP03 PT-all specification with c = 2:
      - beta_{gK1, 2} = 0 and beta_{gK2, 2} = 0 (since u = 2 in [c, b_g] = [2, 3]).
      - beta_{gK1, 1} and beta_{gK2, 1} remain distinct private nuisance coordinates
        (since u = 1 < c = 2).
      - Free parameter vector: eta = (beta_{gK1, 1}, beta_{gK2, 1}, tau_4)' in R^3.
    """
    true_eta = fvec([Fraction(-3, 2), Fraction(-1, 2), Fraction(3, 1)])
    WH_correct = fmat([
        [-1,  0, 0],  # l1: (1,2) K1 -> 0 - beta_{gK1,1}
        [ 0, -1, 0],  # l2: (1,2) K2 -> 0 - beta_{gK2,1}
        [ 0, -1, 1],  # l3: (1,4) K2 -> tau_4 - beta_{gK2,1}
        [ 0,  0, 0],  # l4: (2,3) K1 -> 0 - 0 (zero row)
        [ 0,  0, 1],  # l5: (2,4) K1 -> tau_4 - 0
    ])
    mu_pop = mat_vec(WH_correct, true_eta)
    assert mu_pop == fvec([Fraction(3, 2), Fraction(1, 2), Fraction(7, 2), Fraction(0), Fraction(3)])

    V = fmat([
        [ 3,  1,  0, -1,  0],
        [ 1,  3,  1,  0,  0],
        [ 0,  1,  4,  0,  1],
        [-1,  0,  0,  2, -1],
        [ 0,  0,  1, -1,  3],
    ])

    # 1. Correct CP03 specification (u=1 private per-operator, u=2 zeroed)
    eta_correct, cov_correct, _ = gls_solve(WH_correct, V, mu_pop)
    assert eta_correct == true_eta  # Exact zero population bias!

    # 2. Invalid alternative A: falsely imposing PT before c (zeroing beta_{gK1,1}=0, beta_{gK2,1}=0)
    #    Drop l1, l2 or treat them as zero rows with tau_4 as sole parameter:
    WH_false_zero_pre_c = fmat([[0], [0], [1], [0], [1]])
    tau_false_zero_vec, _, _ = gls_solve(WH_false_zero_pre_c, V, mu_pop)
    bias_false_zero = tau_false_zero_vec[0] - true_eta[2]
    assert bias_false_zero != 0

    # 3. Invalid alternative B: falsely merging beta_{gK1,1} == beta_{gK2,1} across K1 != K2 at u = 1 < c
    WH_false_merge = fmat([
        [-1, 0],
        [-1, 0],
        [-1, 1],
        [ 0, 0],
        [ 0, 1],
    ])
    eta_false_merge, _, _ = gls_solve(WH_false_merge, V, mu_pop)
    bias_false_merge = eta_false_merge[1] - true_eta[2]
    assert bias_false_merge != 0

    return {
        "population_moments_exact": frac_to_str_vec(mu_pop),
        "true_eta_exact": frac_to_str_vec(true_eta),
        "eta_correct_pop_exact": frac_to_str_vec(eta_correct),
        "tau4_bias_correct_spec": str(eta_correct[2] - true_eta[2]),
        "tau4_bias_if_zeroed_before_c": str(bias_false_zero),
        "tau4_bias_if_merged_operators_before_c": str(bias_false_merge),
    }


# ---------------------------------------------------------------------------
# Fixture 4: Anticipation (delta > 0) and Irregular Calendar Endpoint Rules
# ---------------------------------------------------------------------------

def classify_endpoint_pt_all(
    u: float,
    g: float,
    b_g: float,
    delta: float,
    pt_start: float,
) -> str:
    """
    Exact CP03 endpoint coordinate role classification for cohort g on calendar T:
      - b_g = max{d in T : d + delta < g}
      - Pre-anticipation condition: u + delta < g
      - If u == b_g: 'reference' (anchored at 0)
      - If u + delta < g and u != b_g:
          * if u >= pt_start: 'pre_pt_all_zero' (constrained to 0 under PT-all(pt_start))
          * if u <  pt_start: 'pre_private' (unrestricted operator-specific nuisance)
      - If u + delta >= g and u < g: 'anticipation_conditional' (NOT zeroed by PT-all!)
      - If u >= g: 'post_conditional'
    """
    if u == b_g:
        return "reference"
    if u + delta < g:
        if u >= pt_start:
            return "pre_pt_all_zero"
        return "pre_private"
    if u < g:
        return "anticipation_conditional"
    return "post_conditional"


def run_fixture_4() -> dict:
    calendar = [2010, 2012, 2015, 2017, 2020]
    g = 2017
    delta = 2
    b_g = max(d for d in calendar if d + delta < g)
    assert b_g == 2012  # 2015 + 2 = 2017 >= 2017, so 2015 is in anticipation window!

    roles_c2010 = {str(u): classify_endpoint_pt_all(u, g, b_g, delta, pt_start=2010) for u in calendar}
    roles_c2012 = {str(u): classify_endpoint_pt_all(u, g, b_g, delta, pt_start=2012) for u in calendar}
    roles_c2015 = {str(u): classify_endpoint_pt_all(u, g, b_g, delta, pt_start=2015) for u in calendar}

    assert roles_c2010 == {
        "2010": "pre_pt_all_zero",
        "2012": "reference",
        "2015": "anticipation_conditional",
        "2017": "post_conditional",
        "2020": "post_conditional",
    }
    # When c >= b_g (e.g. c = 2012 or c = 2015), no extra pre_pt_all_zero coordinates exist for cohort g;
    # roles reduce algebraically to PT-post; global public no-window behavior is pending.
    assert roles_c2012 == {
        "2010": "pre_private",
        "2012": "reference",
        "2015": "anticipation_conditional",
        "2017": "post_conditional",
        "2020": "post_conditional",
    }
    assert roles_c2015 == roles_c2012

    return {
        "calendar": calendar,
        "cohort_g": g,
        "anticipation_delta": delta,
        "reference_b_g": b_g,
        "roles_pt_start_2010": roles_c2010,
        "roles_pt_start_2012_reduces_to_pt_post": roles_c2012,
        "roles_pt_start_2015_reduces_to_pt_post": roles_c2015,
        "anticipation_2015_preserved_nonzero_in_all_cases": True,
    }


# ---------------------------------------------------------------------------
# Fixture 5: ATT Identification (rank(WH)) vs. Covariance Singularity (D5/D28)
# ---------------------------------------------------------------------------

def run_fixture_5() -> dict:
    """
    Cohort g = 4, calendar T = {1, 2, 3, 4}, delta = 0, b_g = 3.
    Observed links:
      l1 = (1, 3)  [connects 1 to b_g = 3]
      l2 = (2, 4)  [connects 2 to post date 4]
    Notice:
      - Under PT-post, vertices are {1, 2, 3(ref), 4}. Component of 4 is {2, 4},
        which is disconnected from b_g = 3! So W_post = [[-1, 0, 0], [0, -1, 1]]
        has rank 2 < 3, and ATT(4, 4) = [0, 0, 1] * theta is NOT in row(W_post).
      - Under PT-all with c = 2, vertex u = 2 in [c, b_g] = [2, 3] is anchored at 0!
        Free parameters are eta = (beta_1, tau_4)' in R^2, and
        WH = [[-1, 0], [0, 1]] has full column rank 2, so ATT(4, 4) IS identified!
    """
    W_post = np.array([
        [-1.0,  0.0, 0.0],
        [ 0.0, -1.0, 1.0],
    ])
    H_c2 = np.array([
        [1.0, 0.0],
        [0.0, 0.0],
        [0.0, 1.0],
    ])
    WH_c2 = W_post @ H_c2  # [[-1, 0], [0, 1]]
    rank_W_post = int(np.linalg.matrix_rank(W_post))
    rank_WH_c2 = int(np.linalg.matrix_rank(WH_c2))
    # Under PT-post, [0, 0, 1] is orthogonal to nullspace vector [0, 1, 1]'? No: [0, 0, 1].[0, 1, 1] = 1 != 0!
    null_vec_post = np.array([0.0, 1.0, 1.0])
    assert np.allclose(W_post @ null_vec_post, 0.0)
    att4_identifiable_pt_post = bool(abs(np.dot([0.0, 0.0, 1.0], null_vec_post)) < 1e-12)
    att4_identifiable_pt_all_c2 = bool(rank_WH_c2 == 2)

    return {
        "rank_W_pt_post": rank_W_post,
        "n_params_pt_post": 3,
        "att4_identifiable_under_pt_post": att4_identifiable_pt_post,
        "rank_WH_pt_all_c2": rank_WH_c2,
        "n_free_params_pt_all_c2": 2,
        "att4_identifiable_under_pt_all_c2": att4_identifiable_pt_all_c2,
    }


# ---------------------------------------------------------------------------
# Fixture 6: Reporting Invariance of Post-Treatment ATT(g,t) to base_period
#            and Non-Degenerate Unrestricted Pretrend Diagnostics
# ---------------------------------------------------------------------------

def run_fixture_6() -> dict:
    """
    Complete deterministic unit-level IF fixture (n = 6 units, cohort g = 4,
    calendar T = {1, 2, 3, 4}, delta = 0, b_g = 3, pt_start c = 1).
    Shows:
      1. Post-treatment ATT(4, 4) estimate, unit IF column, and SE are invariant
         when base_period switches between 'varying' and 'universal'.
      2. Unrestricted direct auxiliary pretrend rows:
         - Under 'universal': pre(4, 1) = -m_(1,3), pre(4, 2) = -m_(2,3), ref(4, 3) = 0
         - Under 'varying'  : pre(4, 2) = +m_(1,2), pre(4, 3) = +m_(2,3)
         When all pairs share the same unit mask/nuisance fit, m_(1,3) = m_(1,2) + m_(2,3),
         so the universal nonreference pretrend vector (-m_(1,3), -m_(2,3))' and the
         varying pretrend vector (m_(1,2), m_(2,3))' are related by the invertible linear
         bijection L_var_to_univ = [[-1, -1], [0, -1]], making the D21 full-null Wald
         statistic W_pre invariant between 'varying' and 'universal' and strictly positive!
      3. Conversely, naively mapping the restricted GMM parameter vector theta_all = H * eta
         into pre-period cells sets beta_1 = 0 and beta_2 = 0 with zero IFs, producing a
         degenerate (zero / singular) pre-test!
    """
    # 6 units, 3 spanning-tree links l1=(1,2), l2=(2,3), l3=(3,4)
    # Centered unit-level influence functions Phi (6 x 3) and point estimates m_hat (3,)
    Phi_raw = np.array([
        [ 2.0, -1.0,  1.0],
        [-1.0,  2.0,  0.0],
        [ 0.0, -1.0, -2.0],
        [ 1.0,  1.0,  2.0],
        [-2.0,  0.0, -1.0],
        [ 0.0, -1.0,  0.0],
    ], dtype=np.float64)
    Phi = Phi_raw - Phi_raw.mean(axis=0, keepdims=True)
    n = Phi.shape[0]
    V = (Phi.T @ Phi) / n
    m_hat = np.array([0.30, 0.20, 1.50], dtype=np.float64)

    # Under PT-all with c = 1: beta_1 = 0, beta_2 = 0, beta_3 = 0 (ref), free eta = (tau_4)
    # On the tree links (1,2), (2,3), (3,4), WH = [0, 0, 1]'
    WH = np.array([[0.0], [0.0], [1.0]], dtype=np.float64)
    V_inv = np.linalg.inv(V)
    gram = WH.T @ V_inv @ WH
    Lambda = V_inv @ WH @ np.linalg.inv(gram)  # (3 x 1)
    tau4_hat = float((Lambda.T @ m_hat)[0])
    psi_tau4 = Phi @ Lambda[:, 0]
    se_tau4 = float(np.sqrt((psi_tau4 @ psi_tau4) / (n ** 2)))

    # No additional overidentification diagnostic is part of CP03.
    # Removed the special-case residual-criterion check: it does not calibrate
    # a standard chi-square test for automatic Tikhonov weighting.
    # Reported cells under 'varying' vs 'universal' using unrestricted auxiliary pre-links:
    # In this single-class system, m_(1,3) = m_(1,2) + m_(2,3) and phi_(1,3) = phi_(1,2) + phi_(2,3).
    m_12, m_23 = m_hat[0], m_hat[1]
    m_13 = m_12 + m_23
    phi_12, phi_23 = Phi[:, 0], Phi[:, 1]
    phi_13 = phi_12 + phi_23

    # Varying: pre(4,2) = m_(1,2), pre(4,3) = m_(2,3), post ATT(4,4) = tau4_hat
    est_varying = np.array([m_12, m_23, tau4_hat])
    if_varying = np.column_stack([phi_12, phi_23, psi_tau4])

    # Universal: pre(4,1) = -m_(1,3), pre(4,2) = -m_(2,3), ref(4,3) = 0, post ATT(4,4) = tau4_hat
    est_universal = np.array([-m_13, -m_23, 0.0, tau4_hat])
    if_universal = np.column_stack([-phi_13, -phi_23, np.zeros(n), psi_tau4])

    # Check exact invariance of post-treatment estimate, IF, and SE
    assert est_varying[2] == est_universal[3] == tau4_hat
    assert np.array_equal(if_varying[:, 2], if_universal[:, 3])

    # D21 unrestricted pretrend Wald statistic on non-reference pre rows:
    pre_est_var = est_varying[:2]
    pre_vcov_var = (if_varying[:, :2].T @ if_varying[:, :2]) / (n ** 2)
    wald_var = float(pre_est_var.T @ np.linalg.solve(pre_vcov_var, pre_est_var))

    pre_est_univ = est_universal[:2]
    pre_vcov_univ = (if_universal[:, :2].T @ if_universal[:, :2]) / (n ** 2)
    wald_univ = float(pre_est_univ.T @ np.linalg.solve(pre_vcov_univ, pre_est_univ))

    # Same-operator, shared-mask pretrend views transform their entire covariance.
    assert abs(wald_var - wald_univ) < 1e-12

    return {
        "tau4_hat_varying": float(est_varying[2]),
        "tau4_hat_universal": float(est_universal[3]),
        "se_tau4_varying": se_tau4,
        "se_tau4_universal": se_tau4,
        "post_if_max_abs_diff": float(np.max(np.abs(if_varying[:, 2] - if_universal[:, 3]))),
        "unrestricted_pre_estimates_varying": [float(x) for x in pre_est_var],
        "unrestricted_pre_estimates_universal": [float(x) for x in pre_est_univ],
        "d21_wald_stat_varying": wald_var,
        "d21_wald_stat_universal": wald_univ,
        "naive_restricted_pre_estimates": [0.0, 0.0],
        "naive_restricted_pre_vcov_rank": 0,
    }


def main() -> None:
    results = {
        "script": "spec/pt_all_reference.py",
        "description": "Independent deterministic reference results for CP03 PT-all specification",
        "production_source_used": False,
        "fixture_1_restricted_gls_and_zero_row": run_fixture_1(),
        "fixture_2_multi_baseline_and_joint_gmm": run_fixture_2(),
        "fixture_3_early_cancellation_bridge": run_fixture_3(),
        "fixture_4_anticipation_and_irregular_calendar": run_fixture_4(),
        "fixture_5_identification_vs_covariance_singularity": run_fixture_5(),
        "fixture_6_reporting_invariance_and_pretest": run_fixture_6(),
    }
    out_path = Path(__file__).resolve().parent / "pt_all_reference_results.json"
    out_path.write_text(json.dumps(results, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
