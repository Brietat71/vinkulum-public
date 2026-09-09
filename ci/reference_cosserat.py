"""Référence statique continue : équations d'équilibre de Cosserat par tir.

Indépendante des éléments finis et des deux moteurs. Charges terminales
conservatives en force ; pas de charge distribuée. Matériau linéaire en
déformations, cinématique tridimensionnelle géométriquement exacte.
Equilibre : n'=0, m'=-r'×n, r'=R(e1+Cn^-1 R^T n), R'=R[Cm^-1 R^T m]×.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import root


def cross(v):
    x,y,z=v
    return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])


def reference(tolerance):
    length=.508
    cn=np.array([2.84191e6,6.40131e5,9.03881e5])
    cm=np.array([3.10338,36.2794,2.42873])
    force=np.array([0.,8.896/np.sqrt(2),8.896/np.sqrt(2)])
    root_moment=np.zeros(3)
    stages=[]
    start=time.perf_counter()
    for stage in range(1,17):
        applied=force*stage/16
        def integrate(m0,dense=False):
            def rhs(s,y):
                r=y[3:12].reshape(3,3)
                dr=r@(np.array([1.,0.,0.])+(r.T@applied)/cn)
                curvature=(r.T@y[12:])/cm
                return np.r_[dr,(r@cross(curvature)).ravel(),-np.cross(dr,applied)]
            y0=np.r_[np.zeros(3),np.eye(3).ravel(),m0]
            sol=solve_ivp(rhs,[0.,length],y0,method='DOP853',rtol=tolerance,atol=tolerance*.01,dense_output=dense)
            if not sol.success: raise ValueError(sol.message)
            return sol
        solved=root(lambda m0:integrate(m0).y[12:,-1],root_moment,tol=1e-10)
        sol=integrate(solved.x,dense=stage==16)
        residual=float(np.max(np.abs(sol.y[12:,-1])))
        if residual>1e-10: raise ValueError(f'tir non convergé au palier {stage}: {residual}')
        root_moment=solved.x
        stages.append(dict(palier=stage,residu_moment_nm=residual,evaluations_tir=solved.nfev))
    grid=np.linspace(0,length,961)
    y=sol.sol(grid)
    rotations=y[3:12].T.reshape(-1,3,3)
    orthogonality=float(np.max(np.abs(np.einsum('nji,njk->nik',rotations,rotations)-np.eye(3))))
    # n=F constant ; conservation du moment spatial m+r×F.
    balance=y[12:].T+np.cross(y[:3].T,force)-root_moment
    return dict(tolerance_ode=tolerance,temps_s=time.perf_counter()-start,
                coordonnee=grid.tolist(),positions=y[:3].T.tolist(),rotations=rotations.tolist(),
                moments=y[12:].T.tolist(),moment_racine=root_moment.tolist(),paliers=stages,
                defaut_orthogonalite=orthogonality,defaut_bilan_moment_nm=float(np.max(np.abs(balance))))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('sortie',type=Path)
    args=p.parse_args()
    runs=[reference(tol) for tol in (1e-10,2e-12)]
    change=float(np.max(np.abs(np.asarray(runs[0]['positions'])-runs[1]['positions'])))
    if change>1e-9: raise ValueError(f'référence spatiale insuffisamment convergée: {change}')
    for r in runs:
        if r['defaut_orthogonalite']>1e-8 or r['defaut_bilan_moment_nm']>1e-9:
            raise ValueError('invariant du continuum non vérifié')
    report=dict(definition=__doc__,programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                ecart_raffinement_m=change,references=runs)
    args.sortie.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(ecart_raffinement_m=change,position=runs[-1]['positions'][-1],
                         defaut_bilan_moment_nm=runs[-1]['defaut_bilan_moment_nm'])))


if __name__=='__main__': main()
