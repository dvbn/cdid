"""Exact rational check of the Prop. 5 class-potential algebra (CP02).

Independent of the R package. At ONE fixed shared nuisance (rates from the
common mask, fixed positive odds, H, L) the accepted weighted-pair link IF
(cp02_weighted_pair_if.tex) is evaluated in exact rationals for links
(1,2), (2,3), (1,3). With an identical risk-set mask the identities
mu_13 = mu_12 + mu_23 and phi_13 = phi_12 + phi_23 hold EXACTLY for every unit,
including pair non-observers and an outside-risk unit. When one risk unit
misses date 3, per-pair masks differ and the identity fails exactly.
Fixture values are arbitrary test data, not a sampling or solver claim.
Run: python3 -I cp02_shared_ipw_exact_reference.py
"""
from fractions import Fraction as F
import sys

G_T, G_C1, G_C2, G_OUT = 3, 5, 99, 7  # target, controls K=(5, 99), outside
COHORTS = [G_T, G_C1, G_C2]


def unit(g, w, x, odds, obs):
    return {"g": g, "w": F(w), "z": (F(1), F(x)), "o": F(odds), "obs": obs}


def fixture(drop_unit_date3=False):
    # obs: date -> outcome. Risk units observe all of {1,2,3} or at most one.
    U = [
        unit(G_T, "3/2", "1/3", "2", {1: F(1), 2: F(3), 3: F(7), 4: F(2)}),
        unit(G_T, "1", "-1/2", "1/2", {1: F(-2), 2: F(1), 3: F(4)}),
        unit(G_T, "2", "2", "3", {1: F(5), 2: F(5, 2), 3: F(9)}),
        unit(G_T, "1/2", "0", "5/4", {2: F(1), 4: F(3)}),          # non-observer
        unit(G_C1, "1", "1", "1/3", {1: F(0), 2: F(2), 3: F(1)}),
        unit(G_C1, "3", "-1", "2/5", {1: F(4), 2: F(1), 3: F(-1), 4: F(1)}),
        unit(G_C1, "2/3", "1/4", "3/7", {3: F(2), 4: F(6)}),        # non-observer
        unit(G_C2, "1", "-2", "1/4", {1: F(1), 2: F(1, 2), 3: F(3)}),
        unit(G_C2, "5/4", "3/2", "2/3", {1: F(-1), 2: F(2), 3: F(2)}),
        unit(G_C2, "1", "1/2", "1", {1: F(6), 4: F(1)}),            # non-observer
        unit(G_C2, "2", "0", "1/2", {1: F(2), 2: F(0), 3: F(5)}),
        unit(G_C1, "1/2", "3", "1/5", {1: F(1), 2: F(4), 3: F(2)}),
        unit(G_T, "1", "-1", "7/5", {1: F(0), 2: F(2), 3: F(5)}),
        unit(G_OUT, "4", "1", "9", {1: F(10), 2: F(-3), 3: F(8)}),  # outside risk
    ]
    if drop_unit_date3:
        del U[11]["obs"][3]
        U[11]["obs"][4] = F(0)
    return U


def solve2(A, b):
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    if det == 0:
        raise SystemExit("singular H")
    return ((A[1][1] * b[0] - A[0][1] * b[1]) / det, (-A[1][0] * b[0] + A[0][0] * b[1]) / det)


def link(U, s, t, S=None):
    n = len(U)
    if S is None:
        S = [1 if (s in u["obs"] and t in u["obs"]) else 0 for u in U]
    delta = [(u["obs"][t] - u["obs"][s]) if S[i] else F(0) for i, u in enumerate(U)]
    I = [[1 if u["g"] == k else 0 for k in COHORTS] for u in U]
    m = len(COHORTS)
    q = [sum(U[i]["w"] * I[i][j] for i in range(n)) / n for j in range(m)]
    r = [sum(U[i]["w"] * I[i][j] * S[i] for i in range(n)) / sum(U[i]["w"] * I[i][j] for i in range(n))
         for j in range(m)]
    R = [sum(I[i]) for i in range(n)]
    D = [I[i][0] for i in range(n)]
    Cc = [R[i] - D[i] for i in range(n)]
    rG = [sum(I[i][j] * r[j] for j in range(m)) if R[i] else F(1) for i in range(n)]
    p = [u["o"] / (1 + u["o"]) for u in U]
    H = [[sum(U[i]["w"] * R[i] * S[i] * p[i] * (1 - p[i]) * U[i]["z"][a] * U[i]["z"][b] / rG[i]
              for i in range(n)) / n for b in range(2)] for a in range(2)]
    L = [[-sum(U[i]["w"] * I[i][j] * S[i] * U[i]["z"][a] * (D[i] - p[i]) for i in range(n)) / n / r[j] ** 2
          for j in range(m)] for a in range(2)]
    xi_r = [[U[i]["w"] * I[i][j] * (S[i] - r[j]) / q[j] for j in range(m)] for i in range(n)]
    xi_b = []
    for i in range(n):
        u_i = [U[i]["w"] * R[i] * S[i] * U[i]["z"][a] * (D[i] - p[i]) / rG[i] for a in range(2)]
        rhs = [u_i[a] + sum(L[a][j] * xi_r[i][j] for j in range(m)) for a in range(2)]
        xi_b.append(solve2(H, rhs))
    a_p = [U[i]["w"] * D[i] * S[i] / r[0] for i in range(n)]
    a_m = [U[i]["w"] * Cc[i] * S[i] * U[i]["o"] / rG[i] for i in range(n)]
    dp, dm = sum(a_p) / n, sum(a_m) / n
    T = sum(a_p[i] * delta[i] for i in range(n)) / sum(a_p)
    M = sum(a_m[i] * delta[i] for i in range(n)) / sum(a_m)
    v = [sum(a_m[i] * U[i]["z"][a] * (delta[i] - M) for i in range(n)) / n / dm for a in range(2)]
    gam = [F(0)] + [sum(U[i]["w"] * I[i][j] * S[i] * U[i]["o"] * (delta[i] - M) for i in range(n)) / n
                    / (r[j] ** 2 * dm) for j in range(1, m)]
    phi = [a_p[i] / dp * (delta[i] - T) - a_m[i] / dm * (delta[i] - M)
           + sum(gam[j] * xi_r[i][j] for j in range(m)) - sum(v[a] * xi_b[i][a] for a in range(2))
           for i in range(n)]
    return {"S": S, "mu": T - M, "phi": phi, "R": R}


def main():
    ok = True
    U = fixture()
    l12, l23, l13 = link(U, 1, 2), link(U, 2, 3), link(U, 1, 3)
    risk = [i for i in range(len(U)) if l12["R"][i]]
    same_mask = all(l12["S"][i] == l23["S"][i] == l13["S"][i] for i in risk)
    mu_gap = l13["mu"] - l12["mu"] - l23["mu"]
    phi_gap = [l13["phi"][i] - l12["phi"][i] - l23["phi"][i] for i in range(len(U))]
    nonobs = [i for i in risk if l12["S"][i] == 0]
    outside = [i for i in range(len(U)) if not l12["R"][i]]
    print("balanced class: identical risk mask:", same_mask, "; non-observers:", len(nonobs),
          "; outside-risk units:", len(outside))
    print("exact mu_13 - mu_12 - mu_23 =", mu_gap)
    print("exact max |phi_13 - phi_12 - phi_23| =", max(abs(x) for x in phi_gap))
    print("mu_12 =", l12["mu"])
    ok &= same_mask and mu_gap == 0 and all(x == 0 for x in phi_gap)
    ok &= all(l12["phi"][i] == 0 and l13["phi"][i] == 0 for i in outside)
    ok &= any(l12["phi"][i] != 0 for i in nonobs)
    # The fixed odds are an approximate shared fit, not a logit root: the
    # algebra must hold without score-root (or IF-mean-zero) properties.
    print("outside-risk IF exactly zero; non-observer IF nonzero:", ok)
    # Near-redundant mutation: one risk unit misses date 3 -> masks differ.
    V = fixture(drop_unit_date3=True)
    m12, m23, m13 = link(V, 1, 2), link(V, 2, 3), link(V, 1, 3)
    diff_mask = any(m12["S"][i] != m13["S"][i] for i in range(len(V)))
    gap = max(abs(m13["phi"][i] - m12["phi"][i] - m23["phi"][i]) for i in range(len(V)))
    mgap = m13["mu"] - m12["mu"] - m23["mu"]
    print("mutated (unit 12 misses date 3): masks differ:", diff_mask,
          "; exact mu gap:", mgap, "; max phi gap:", gap)
    ok &= diff_mask and gap != 0 and mgap != 0
    if not ok:
        print("FAIL")
        sys.exit(1)
    print("EXACT PASS: Prop. 5 telescoping identities hold exactly at one shared nuisance;"
          " the split-mask mutation breaks them.")


if __name__ == "__main__":
    main()
