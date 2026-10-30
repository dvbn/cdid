"""Finite direct-pair pooling reference. No RNG, production or candidate code."""
import sys
sys.dont_write_bytecode = True
import csv
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
DATE_REF = HERE.parent/'independent_reference'
pin = json.loads((DATE_REF/'receipt.json').read_text())['output_hashes']['reference.py']
assert hashlib.sha256((DATE_REF/'reference.py').read_bytes()).hexdigest() == pin
module_spec = importlib.util.spec_from_file_location('frozen_dense_logit',DATE_REF/'reference.py')
dense = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(dense)
TOL, ROOT_TOL, STEPS = dense.TOL, dense.ROOT_TOL, dense.STEPS
N, K, X, R = dense.N, dense.K, dense.X, dense.R
A, B = (R.sum(axis=1)>=2).astype(float), (R.sum(axis=1)>=1).astype(float)
S = np.column_stack([R[:,s]*R[:,t] for s,t in dense.PAIRS])
dense.R = S  # In-memory response argument adaptation of the frozen dense solver.
PROBES = dense.PROBES
cell_mass = np.array([[.21,.29],[.31,.19]])
population_pair = dense.expit(np.array([-2.2,-1.9])[K]+.3*X)
# P111=.5p and each pair-only atom=.5p; repeat mass=2p, observed mass=2.5p.
history_mass = np.zeros(N)
for i, history in enumerate(R):
    ones = int(history.sum()); p = population_pair[i]
    history_mass[i] = .5*p if ones>=2 else (1-2.5*p if ones==0 else
        p*np.array([.2,.15,.15])[int(np.argmax(history))])
W = N*cell_mass[K,(X>0).astype(int)]*history_mass
assert abs(W.mean()-1)<ROOT_TOL and np.all(W>0)
MASKS = dict(F=np.ones(N), B=B, A=A)
RHO = {}
for u, mask in MASKS.items():
    up = np.ones(N) if u=='F' else (2.5 if u=='B' else 2.)*population_pair
    rho = np.ones(N)
    for k in range(2):
        ik = K==k
        ubar = np.sum(W*ik*mask)/np.sum(W*ik)
        rho[ik] = up[ik]/ubar
        assert abs(np.sum(W*ik*rho)/np.sum(W*ik)-1)<ROOT_TOL
    RHO[u]=rho


def design(mode,k,x):
    slopes = 6 if mode=='separate' else 2 if mode=='within_cohort' else 1
    z=np.zeros((len(k),3,6+slopes))
    for i in range(len(k)):
        for e in range(3):
            cell=3*k[i]+e
            z[i,e,cell]=1
            slope=cell if mode=='separate' else k[i] if mode=='within_cohort' else 0
            z[i,e,6+slope]=x[i]
    return z


def separate_closed(w,cw):
    theta=np.zeros(12)
    for k in range(2):
        for e in range(3):
            logodds=[]
            for x in (-1.,1.):
                cell=(K==k)&(X==x)
                pr=np.dot(cw*w*cell,S[:,e])/np.sum(cw*w*cell)
                logodds.append(np.log(pr/(1-pr)))
            b=3*k+e
            theta[b]=(logodds[0]+logodds[1])/2
            theta[6+b]=(logodds[1]-logodds[0])/2
    return theta


def finite(mode,u,cw=None):
    cw=np.ones(N) if cw is None else cw
    mask=MASKS[u]; z=design(mode,K,X)
    theta=dense.fit(z,W*mask,cw)
    if mode=='separate':
        assert np.max(np.abs(theta-separate_closed(W*mask,cw)))<ROOT_TOL
    a=np.array([np.dot(cw*W*(K==k)*mask,A)/np.sum(cw*W*(K==k)*mask) for k in range(2)])
    k=np.array([k for k,x in PROBES]); x=np.array([x for k,x in PROBES])
    pp=dense.expit(design(mode,k,x)@theta).reshape(-1)
    rho=np.repeat([RHO[u][np.flatnonzero((K==kk)&(X==xx))[0]] for kk,xx in PROBES],3)
    v=np.repeat(a[k],3)/(rho*pp)
    return np.r_[theta,a,pp,v]


def check(mode,u):
    value=finite(mode,u); z=design(mode,K,X); d=z.shape[2]
    theta,a=value[:d],value[d:d+2]
    p=dense.expit(z@theta); mask=MASKS[u]
    expected=np.zeros(d)
    if u=='F':
        expected[:6]=np.repeat([-2.2,-1.9],3); expected[6:]=.3
    elif u=='B':
        expected[:6]=np.log(.4/.6)
    assert np.max(np.abs(theta-expected))<ROOT_TOL
    by_pair=W[:,None,None]*mask[:,None,None]*z*(S-p)[:,:,None]
    score=by_pair.sum(axis=1)
    J=np.einsum('n,ne,ned,nef->df',W*mask,p*(1-p),z,z)/N
    xt=np.linalg.solve(J,score.T).T
    parts=np.stack([np.linalg.solve(J,by_pair[:,e].T).T for e in range(3)],axis=1)
    xa=np.column_stack([W*(K==k)*mask*(A-a[k])/np.mean(W*(K==k)*mask) for k in range(2)])
    kp=np.array([k for k,x in PROBES]); xp=np.array([x for k,x in PROBES])
    zp=design(mode,kp,xp)
    pp=dense.expit(zp@theta).reshape(-1)
    grad=(pp*(1-pp))[:,None]*zp.reshape(-1,d)
    xprob=xt@grad.T
    rho=np.repeat([RHO[u][np.flatnonzero((K==kk)&(X==xx))[0]] for kk,xx in PROBES],3)
    repeat=np.repeat(a[kp],3)
    xrepeat=np.column_stack([xa[:,k] for k,x in PROBES for _ in range(3)])
    xv=xrepeat/(rho*pp)-repeat*xprob/(rho*pp**2)
    analytic=np.c_[xt,xa,xprob,xv]
    assert np.max(np.abs(analytic.mean(axis=0)))<ROOT_TOL
    errors=[]
    for h in STEPS:
        fd=np.zeros_like(analytic)
        for i in range(N):
            cp=np.ones(N); cm=cp.copy(); cp[i]+=h; cm[i]-=h
            fd[i]=N*(finite(mode,u,cp)-finite(mode,u,cm))/(2*h)
        err=float(np.max(np.abs(fd-analytic)))
        assert err<TOL,(mode,u,h,err)
        errors.append(dict(step=h,max_error=err))
    # Pair 12 predictions must also use the other pairs' scores when slopes share.
    retained_pair=np.linalg.solve(J,by_pair[:,0].T).T@grad[[0,3]].T
    omission=float(np.max(np.abs(xprob[:,[0,3]]-retained_pair)))
    assert omission<TOL if mode=='separate' else omission>TOL
    cross_cohort=float(np.max(np.abs(xprob[K==1,0])))
    assert cross_cohort>TOL if mode=='across_cohorts' else cross_cohort<ROOT_TOL
    covariance=xt.T@xt/N**2
    copied=sum(parts[:,e].T@parts[:,e] for e in range(3))/N**2
    covariance_loss=float(np.max(np.abs(covariance-copied)))
    if u=='A':
        # Under this constructed A law, the three S_e are pairwise uncorrelated:
        # each marginal is1/2 and each pairwise joint is1/4. Retain truthful zero.
        assert covariance_loss<ROOT_TOL,(mode,u,covariance_loss)
    else:
        assert covariance_loss>TOL,(mode,u,covariance_loss)
    if u=='A':
        assert np.max(np.abs(xa))<ROOT_TOL and np.max(np.abs(analytic[A==0]))<ROOT_TOL
    else:
        assert np.max(np.abs(analytic[A==0]))>TOL
    separate_error=None
    if mode=='separate':
        separate_error=0.
        for b in range(6):
            cols=[b,6+b]
            local=np.linalg.solve(J[np.ix_(cols,cols)],score[:,cols].T).T
            separate_error=max(separate_error,float(np.max(np.abs(local-xt[:,cols]))))
        assert separate_error<ROOT_TOL
    np.savetxt(HERE/f'{mode}_{u}_if.csv',analytic,delimiter=',',fmt='%.17g')
    return dict(coefficients=theta.tolist(),a=a.tolist(),pair_probabilities=pp.tolist(),multipliers=value[d+14:].tolist(),
        finite_derivatives=errors,center_max=float(np.max(np.abs(analytic.mean(axis=0)))),
        dropped_other_pair_score_error=omission,cross_cohort_pair_prediction_IF=cross_cohort,
        independent_pair_copy_covariance_error=covariance_loss,separate_CP04_score_error=separate_error,
        off_A_row_max=float(np.max(np.abs(analytic[A==0]))),dimension=len(value))


def main():
    cases={f'{mode}_{u}':check(mode,u) for mode in ('separate','within_cohort','across_cohorts') for u in ('F','B','A')}
    for u in ('F','B','A'):
        base=np.array(cases[f'separate_{u}']['multipliers'])
        assert max(np.max(np.abs(base-cases[f'{m}_{u}']['multipliers'])) for m in ('within_cohort','across_cohorts'))<ROOT_TOL
    for m in ('separate','within_cohort','across_cohorts'):
        base=np.array(cases[f'{m}_F']['multipliers'])
        assert max(np.max(np.abs(base-cases[f'{m}_{u}']['multipliers'])) for u in ('B','A'))<ROOT_TOL
    with (HERE/'atoms.csv').open('w',newline='') as out:
        writer=csv.writer(out)
        writer.writerow(['id','cohort','x','R1','R2','R3','S12','S13','S23','A','B','w','rho_B','rho_A'])
        writer.writerows([[i+1,K[i],X[i],*R[i],*S[i],A[i],B[i],W[i],RHO['B'][i],RHO['A'][i]] for i in range(N)])
    result=dict(status='PASS',records=N,cases=cases,derivative_tolerance=TOL,root_tolerance=ROOT_TOL,
                frozen_solver_sha256=pin,no_rng_or_production_code=True)
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: 32 fixed original units, 3 pair-pooling modes x F/B/A, 2 finite-difference steps.')
    print('Max full-stack derivative error:',max(r['max_error'] for c in cases.values() for r in c['finite_derivatives']))
    print('Separate modes reduce to CP04; dropped cross-pair scores fail in shared fits. Copy covariance fails on F/B; the constructed A law has truthful zero cross-pair covariance.')
    assert hashlib.sha256((DATE_REF/'reference.py').read_bytes()).hexdigest()==pin


if __name__=='__main__':
    main()
