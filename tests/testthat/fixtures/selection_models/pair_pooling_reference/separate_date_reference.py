"""Only the missing separate-date branch, preserving the frozen date reference."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import numpy as np

HERE=Path(__file__).resolve().parent
REF=HERE.parent/'independent_reference'
pins=json.loads((REF/'receipt.json').read_text())['output_hashes']
assert hashlib.sha256((REF/'reference.py').read_bytes()).hexdigest()==pins['reference.py']
spec=importlib.util.spec_from_file_location('frozen_date_solver',REF/'reference.py')
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
N,K,X,R,A=base.N,base.K,base.X,base.R,base.A
theta0=np.array([-.8,.2,.7,.4,-.5,1.1,.55,-.2,.35,-.35,.45,.15])


def design(k,x):
    z=np.zeros((len(k),3,12))
    for i in range(len(k)):
        for t in range(3):
            b=3*k[i]+t
            z[i,t,b]=1; z[i,t,6+b]=x[i]
    return z


Z=design(K,X)
P=base.expit(Z@theta0)
W=N*np.array([[.21,.29],[.31,.19]])[K,(X>0).astype(int)]*np.prod(P**R*(1-P)**(1-R),axis=1)
KP=np.array([k for k,x in base.PROBES]); XP=np.array([x for k,x in base.PROBES])
ZP=design(KP,XP)


def finite(cw=None):
    cw=np.ones(N) if cw is None else cw
    theta=base.fit(Z,W,cw)
    a=np.array([np.dot(cw*W*(K==k),A)/np.sum(cw*W*(K==k)) for k in range(2)])
    p=base.expit(ZP@theta)
    pair=np.array([p[:,s]*p[:,t] for s,t in base.PAIRS]).T.reshape(-1)
    return np.r_[theta,a,pair,np.repeat(a[KP],3)/pair]


def main():
    value=finite(); theta=value[:12]; a=value[12:14]
    assert np.max(np.abs(theta-theta0))<base.ROOT_TOL
    pp=base.expit(Z@theta)
    score_date=W[:,None,None]*Z*(R-pp)[:,:,None]
    J=np.einsum('n,nt,ntd,nte->de',W,pp*(1-pp),Z,Z)/N
    xt=np.linalg.solve(J,score_date.sum(axis=1).T).T
    parts=np.stack([np.linalg.solve(J,score_date[:,t].T).T for t in range(3)],axis=1)
    xa=np.column_stack([W*(K==k)*(A-a[k])/np.mean(W*(K==k)) for k in range(2)])
    p=base.expit(ZP@theta); gradients=[]
    for j in range(4):
        for s,t in base.PAIRS:
            gradients.append(p[j,s]*p[j,t]*((1-p[j,s])*ZP[j,s]+(1-p[j,t])*ZP[j,t]))
    gradients=np.array(gradients)
    xpair=xt@gradients.T; pair=value[14:26]
    xrepeat=np.column_stack([xa[:,k] for k,x in base.PROBES for _ in range(3)])
    xv=xrepeat/pair-np.repeat(a[KP],3)*xpair/pair**2
    analytic=np.c_[xt,xa,xpair,xv]
    assert np.max(np.abs(analytic.mean(axis=0)))<base.ROOT_TOL
    errors=[]
    for h in base.STEPS:
        fd=np.zeros_like(analytic)
        for i in range(N):
            cp=np.ones(N); cm=cp.copy(); cp[i]+=h; cm[i]-=h
            fd[i]=N*(finite(cp)-finite(cm))/(2*h)
        error=float(np.max(np.abs(fd-analytic)))
        assert error<base.TOL
        errors.append(dict(step=h,max_error=error))
    cov=xt.T@xt/N**2
    copy=sum(parts[:,t].T@parts[:,t] for t in range(3))/N**2
    omitted_cov=float(np.max(np.abs(cov-copy)))
    assert omitted_cov>base.TOL
    # A nonendpoint has no predictor loading with fully separate parameters.
    assert np.max(np.abs(parts[:,2]@gradients[0]))<base.ROOT_TOL
    pair_cov=xpair.T@xpair/N**2
    assert abs(pair_cov[0,1])>base.TOL
    np.savetxt(HERE/'separate_date_F_if.csv',analytic,delimiter=',',fmt='%.17g')
    with (HERE/'separate_date_atoms.csv').open('w',newline='') as output:
        writer=csv.writer(output); writer.writerow(['id','cohort','x','R1','R2','R3','A','w'])
        writer.writerows([[i+1,K[i],X[i],*R[i],A[i],W[i]] for i in range(N)])
    result=dict(status='PASS',records=N,mode='separate_date_F',distinct_date_slopes=True,
        coefficients=theta.tolist(),a=a.tolist(),finite_derivatives=errors,
        center_max=float(np.max(np.abs(analytic.mean(axis=0)))),
        independent_date_copy_covariance_error=omitted_cov,
        pair12_pair13_covariance=float(pair_cov[0,1]),nonendpoint_predictor_loading=0.,
        derivative_tolerance=base.TOL,root_tolerance=base.ROOT_TOL,
        frozen_solver_sha256=pins['reference.py'],dimension=len(value),no_rng_or_production_code=True)
    (HERE/'separate_date_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: only missing separate-date F branch;32 atoms, distinct date/cohort slopes; max derivative error',max(r['max_error'] for r in errors))


if __name__=='__main__':
    main()
