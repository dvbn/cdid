"""Independent fixed-atom date-logit reference. No RNG, package or candidate code."""
import csv
import itertools
import json
import platform
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOL, ROOT_TOL = 3e-7, 2e-11
STEPS = (1e-4, 5e-5)
PAIRS = ((0, 1), (0, 2), (1, 2))
ATOMS = np.array([(k, x, *r) for k in range(2) for x in (-1., 1.)
                  for r in itertools.product((0., 1.), repeat=3)])
K, X, R = ATOMS[:, 0].astype(int), ATOMS[:, 1], ATOMS[:, 2:]
N = len(K)
A = (R.sum(axis=1) >= 2).astype(float)
PROBES = tuple(itertools.product(range(2), (-1., 1.)))


def expit(v):
    return 1 / (1 + np.exp(-v))


def design(mode, k, x):
    size = 6 + (2 if mode == 'within_cohort' else 1 if mode == 'across_cohorts' else 0)
    z = np.zeros((len(k), 3, size))
    for i in range(len(k)):
        for t in range(3):
            z[i, t, 3 * k[i] + t] = 1.
            if mode != 'date_only':
                z[i, t, 6 + (k[i] if mode == 'within_cohort' else 0)] = x[i]
    return z


def exact_law(mode):
    theta = np.array([-.8, .2, .7, .4, -.5, 1.1] +
                     ([.55, -.35] if mode == 'within_cohort' else [.4] if mode == 'across_cohorts' else []))
    z = design(mode, K, X)
    p = expit(z @ theta)
    cell_mass = np.array([[.21, .29], [.31, .19]])
    atom_mass = cell_mass[K, (X > 0).astype(int)] * np.prod(p**R * (1-p)**(1-R), axis=1)
    assert abs(atom_mass.sum() - 1.) < ROOT_TOL
    # H is uniform on the enumerated atoms; fixed w turns H into the declared P.
    w = N * atom_mass
    return theta, z, w


def fit(z, w, cw):
    theta = np.zeros(z.shape[2])
    for iteration in range(40):
        p = expit(z @ theta)
        unit_score = w[:, None] * np.einsum('ntd,nt->nd', z, R-p)
        score = cw @ unit_score / cw.sum()
        info = np.einsum('n,nt,ntd,nte->de', cw*w, p*(1-p), z, z) / cw.sum()
        step = np.linalg.solve(info, score)
        theta += step
        if np.max(np.abs(step)) < 1e-13:
            break
    else:
        raise AssertionError('Dense Newton did not reach its finite score root')
    assert np.max(np.abs(np.einsum('n,nt,ntd->d', cw*w, R-expit(z@theta), z))) < ROOT_TOL
    return theta


def finite(mode, z, w, cw=None):
    cw = np.ones(N) if cw is None else cw
    theta = fit(z, w, cw)
    a = np.array([np.dot(cw*w*(K==k), A) / np.sum(cw*w*(K==k)) for k in range(2)])
    zp = design(mode, np.array([k for k, x in PROBES]), np.array([x for k, x in PROBES]))
    pp = expit(zp @ theta)
    pairp = np.array([pp[:, s] * pp[:, t] for s, t in PAIRS]).T.reshape(-1)
    repeat = np.repeat(a[[k for k, x in PROBES]], 3)
    mult = repeat / pairp
    return np.r_[theta, a, pairp, mult]


def check(mode):
    truth, z, w = exact_law(mode)
    v = finite(mode, z, w); d = len(truth)
    theta, a = v[:d], v[d:d+2]
    assert np.max(np.abs(theta-truth)) < ROOT_TOL
    p = expit(z @ theta)
    dated_scores = w[:, None, None] * z * (R-p)[:, :, None]
    unit_score = dated_scores.sum(axis=1)
    info = np.einsum('n,nt,ntd,nte->de', w, p*(1-p), z, z) / N
    x_theta = np.linalg.solve(info, unit_score.T).T
    x_date = np.stack([np.linalg.solve(info, dated_scores[:, t].T).T for t in range(3)], axis=1)
    x_a = np.column_stack([w*(K==k)*(A-a[k]) / np.mean(w*(K==k)) for k in range(2)])
    zp = design(mode, np.array([k for k, x in PROBES]), np.array([x for k, x in PROBES]))
    pp = expit(zp @ theta)
    rows, omitted_endpoint, omitted_other_date = [], [], []
    point_pair = []
    for j, (k, x) in enumerate(PROBES):
        for s, t in PAIRS:
            prob = pp[j, s] * pp[j, t]
            grad = prob * ((1-pp[j, s])*zp[j, s] + (1-pp[j, t])*zp[j, t])
            rows.append(x_theta @ grad)
            omitted_endpoint.append(x_theta @ (prob*(1-pp[j, s])*zp[j, s]))
            other = next(q for q in range(3) if q not in (s, t))
            omitted_other_date.append((x_theta-x_date[:, other]) @ grad)
            point_pair.append(prob)
            # P is explicitly constructed to be mutually independent conditional on k,X.
            cell = (K==k) & (X==x)
            joint = np.sum(w[cell]*R[cell, s]*R[cell, t]) / np.sum(w[cell])
            assert abs(joint-prob) < ROOT_TOL
    x_pair = np.array(rows).T
    repeat_if = np.column_stack([x_a[:, k] for k, x in PROBES for _ in PAIRS])
    repeat = np.repeat(a[[k for k, x in PROBES]], 3)
    point_pair = np.array(point_pair)
    x_mult = repeat_if / point_pair - repeat*x_pair / point_pair**2
    analytic = np.column_stack([x_theta, x_a, x_pair, x_mult])
    assert np.max(np.abs(analytic.mean(axis=0))) < ROOT_TOL
    fd_results = []
    for h in STEPS:
        fd = np.zeros_like(analytic)
        for i in range(N):
            cp = np.ones(N); cm = cp.copy(); cp[i] += h; cm[i] -= h
            fd[i] = N*(finite(mode, z, w, cp)-finite(mode, z, w, cm))/(2*h)
        err = float(np.max(np.abs(fd-analytic)))
        assert err < TOL, (mode, h, err)
        fd_results.append(dict(step=h, stack_max_error=err,
            coefficient_max_error=float(np.max(np.abs(fd[:, :d]-x_theta))),
            repeat_max_error=float(np.max(np.abs(fd[:, d:d+2]-x_a))),
            pair_max_error=float(np.max(np.abs(fd[:, d+2:d+14]-x_pair))),
            multiplier_max_error=float(np.max(np.abs(fd[:, d+14:]-x_mult)))))
    cov = x_theta.T @ x_theta / N**2
    copied_cov = sum(x_date[:, t].T @ x_date[:, t] for t in range(3)) / N**2
    cross_cov_error = float(np.max(np.abs(cov-copied_cov)))
    endpoint_error = float(np.max(np.abs(x_pair-np.array(omitted_endpoint).T)))
    copied_score_error = float(np.max(np.abs(x_theta-x_date[:, 0])))
    other_date_error = float(np.max(np.abs(x_pair-np.array(omitted_other_date).T)))
    assert min(cross_cov_error, endpoint_error, copied_score_error) > TOL
    if mode != 'date_only':
        assert other_date_error > TOL
    else:
        assert other_date_error < ROOT_TOL
    off_a = A==0
    assert np.max(np.abs(x_theta[off_a])) > TOL and np.max(np.abs(x_a[off_a])) > TOL
    assert np.max(np.abs(x_mult[off_a])) > TOL
    pair_cov = x_pair.T @ x_pair / N**2
    assert abs(pair_cov[0, 1]) > TOL
    names = ([f'theta{j}' for j in range(d)]+['a0','a1']+
             [f'pair_k{k}_x{x:g}_{s+1}{t+1}' for k,x in PROBES for s,t in PAIRS]+
             [f'mult_k{k}_x{x:g}_{s+1}{t+1}' for k,x in PROBES for s,t in PAIRS])
    np.savetxt(HERE/f'{mode}_if.csv', analytic, delimiter=',', header=','.join(names), comments='', fmt='%.17g')
    return w, dict(coefficients=theta.tolist(), repeat_rates=a.tolist(),
        pair_probabilities=point_pair.tolist(), combined_multipliers=(repeat/point_pair).tolist(),
        information=info.tolist(), finite_derivatives=fd_results,
        max_center_error=float(np.max(np.abs(analytic.mean(axis=0)))),
        mutants=dict(omitted_endpoint=endpoint_error, independent_date_copy_score=copied_score_error,
                     independent_date_covariance=cross_cov_error, unused_endpoint_date_score=other_date_error),
        covariance=cov.tolist(),
        shared_pair_covariance_12_13=float(pair_cov[0,1]),
        off_A_max=dict(coefficient=float(np.max(np.abs(x_theta[off_a]))),
                       empirical_repeat=float(np.max(np.abs(x_a[off_a]))),
                       combined_multiplier=float(np.max(np.abs(x_mult[off_a])))),
        dimension=len(v), records=N)


def main():
    results, weights = {}, []
    for mode in ('within_cohort','across_cohorts','date_only'):
        w, result = check(mode); weights.append(w); results[mode] = result
    with (HERE/'atoms.csv').open('w',newline='') as out:
        writer=csv.writer(out)
        writer.writerow(['id','cohort','x','R1','R2','R3','A','w_within','w_across','w_date_only'])
        writer.writerows([[i+1,*ATOMS[i].tolist(),A[i],*[w[i] for w in weights]] for i in range(N)])
    h=np.array(list(itertools.product((0,1),repeat=3)))
    aa=h.sum(axis=1)>=2
    marginal=float(h[aa,0].mean()); joint=float((h[aa,0]*h[aa,1]).mean())
    assert marginal==.75 and joint==.5 and marginal**2==.5625
    payload=dict(status='PASS',tolerance=TOL,root_tolerance=ROOT_TOL,cases=results,
        repeat_conditioning_half_witness=dict(F_marginal=.5,F_joint=.25,
            A_marginal=marginal,A_joint=joint,product_A_marginals=marginal**2),
        no_rng_or_candidate_code=True,runtime=platform.python_version(),numpy=np.__version__)
    (HERE/'results.json').write_text(json.dumps(payload,indent=2)+'\n')
    print('PASS: 3 fixed-atom models, 32 original units each, 2 finite-difference steps; no RNG.')
    print('Maximum derivative error:',max(r['stack_max_error'] for c in results.values() for r in c['finite_derivatives']))
    print('All endpoint, date-copy score/covariance mutants fail; shared-nuisance off-pair date effect is correctly zero in date-only case.')
    print('Conditioning witness: A marginals=3/4, joint=1/2, product=9/16.')


if __name__=='__main__':
    main()
