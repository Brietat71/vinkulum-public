"""Équilibres et bilans indépendants pour la globalisation de Newton."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq
from scipy.spatial.transform import Rotation
from vinkulum import Noyau


def poutre(spec):
    nb = spec['n']
    length_unit, force_unit = spec.get('unites', [1., 1.])
    q = Rotation.from_rotvec([.4, -.3, .2]).as_matrix() if spec.get('tourne') else np.eye(3)
    offset = q @ [.12, -.21, .34] if spec.get('tourne') else np.zeros(3)
    length = .508
    n = Noyau([0., 0., 0.])
    for i in range(nb+1):
        pos = length_unit*(q @ [length*i/nb, 0., 0.]+offset)
        n.corps(str(i), 1., np.eye(3).ravel().tolist(), pos.tolist(), rot=q.ravel().tolist())
    n.liaison('enc', None, 0, pa=(length_unit*offset).tolist())
    for i in range(nb):
        n.poutre(str(i), i, i+1, 2.84191e6*force_unit, 6.40131e5*force_unit,
                  3.10338*force_unit*length_unit**2, 36.2794*force_unit*length_unit**2,
                  ei3=2.42873*force_unit*length_unit**2, ga3=9.03881e5*force_unit,
                  formulation=spec['formulation'])
    if spec.get('moment'):
        force = np.zeros(3)
        torque = np.array([0., 0., .7*2.42873/length])
    else:
        angle = np.radians(spec.get('angle', 45.))
        force = 8.896*np.array([0., np.sin(angle), np.cos(angle)])
        torque = np.zeros(3)
    n.effort(nb, [0., 0., 0.], [0., 0., 0.])
    loads = .5*(1-np.cos(np.pi*np.arange(1, 51)/50)) if spec.get('rampe') else [1.]
    tol = 1e-8*min(1., length_unit, force_unit)

    def charge(fraction):
        n.pose_effort(0, (fraction*force_unit*(q @ force)).tolist(),
                       (fraction*force_unit*length_unit*(q @ torque)).tolist())

    def canonique():
        s = n.etat_precis()
        p = (np.asarray(s[1], dtype=np.longdouble)+np.asarray(s[6], dtype=np.longdouble))/length_unit
        p = (p-offset) @ q
        r = np.asarray(s[2]).reshape(-1, 3, 3)
        r = np.einsum('ij,bjk->bik', q.T, r)
        return np.asarray(p, dtype=float), r

    def controle():
        p, r = canonique()
        reaction = np.asarray(n.reactions()[0][1])
        rf = q.T @ reaction[:3]/force_unit
        rm = q.T @ reaction[3:]/(force_unit*length_unit)
        result = {'equilibre_force_N': float(np.max(np.abs(rf-force))),
                  'equilibre_moment_Nm': float(np.max(np.abs(rm-(np.cross(p[-1], force)+torque))))}
        if spec.get('moment') and spec['formulation'] == 'integree':
            th = np.linspace(0., .7, nb+1)
            reference = length/.7*np.column_stack((np.sin(th), 1-np.cos(th), np.zeros(nb+1)))
            rotation_ref = Rotation.from_rotvec(np.column_stack((np.zeros(nb+1),np.zeros(nb+1),th))).as_matrix()
            result.update(position_arc_m=float(np.max(np.abs(p-reference))),
                          rotation_arc=float(np.max(np.abs(r-rotation_ref))))
        return result
    return n, loads, charge, canonique, controle, dict(tol=tol, iters=100, paliers_max=1, strict=True)


def cable(spec):
    nb = spec['n']; length, span, mass = 10., 8., 30./nb
    ell = length/nb
    p0 = np.array([[span*k/nb, 0., -3*(1-abs(2*k/nb-1))] for k in range(nb+1)])
    n = Noyau([0., 0., -9.81])
    for k in range(1, nb):
        n.corps(str(k), mass, (1e-9*np.eye(3)).ravel().tolist(), p0[k].tolist())
        n.liaison(f'plan{k}', None, k-1, bloque_t=[1], bloque_r=[0, 1, 2])
    for k in range(nb):
        n.distance(f'barre{k}', None if k==0 else k-1, None if k==nb-1 else k,
                   p0[0].tolist() if k==0 else [0.,0.,0.],
                   p0[-1].tolist() if k==nb-1 else [0.,0.,0.], ell)

    # Polygone funiculaire discret EXACT : efforts verticaux connus par
    # symétrie et poids nodaux ; une seule inconnue, la tension horizontale H.
    vertical = mass*9.81*(np.arange(nb)-(nb-1)/2)
    horizontal = brentq(lambda h: np.sum(ell*h/np.hypot(h,vertical))-span, 1e-6, 1e6, xtol=1e-12)
    segments = ell*np.column_stack((np.full(nb,horizontal), np.zeros(nb), vertical))/np.hypot(horizontal,vertical)[:,None]
    reference = np.vstack((p0[0],p0[0]+np.cumsum(segments,axis=0)))

    def canonique():
        s = n.etat_precis()
        p = np.asarray(s[1])+np.asarray(s[6])
        return p,np.asarray(s[2]).reshape(-1,3,3)

    def controle():
        p,_ = canonique()
        pos = np.vstack((p0[0],p,p0[-1]));diff=np.diff(pos,axis=0)
        u = diff/np.linalg.norm(diff,axis=1)[:,None]
        values=dict(n.reactions()); tensions=np.array([values[f'barre{k}'][0] for k in range(nb)])
        return {'position_reference_m':float(np.max(np.abs(p-reference[1:-1]))),
                'longueur_m':float(np.max(np.abs(np.linalg.norm(diff,axis=1)-ell))),
                'tension_horizontale_N':float(np.max(np.abs(tensions*u[:,0]-horizontal))),
                'equilibre_noeuds_N':float(np.max(np.abs(np.diff(tensions[:,None]*u,axis=0)-[0.,0.,mass*9.81])))}
    return n,[1.],lambda _:None,canonique,controle,dict(tol=1e-10,iters=60,paliers_max=2,strict=True)


def modele(spec):
    if spec['famille']=='poutre':return poutre(spec)
    if spec['famille']=='cable':return cable(spec)
    if spec['famille']=='cascade':
        from diagnostic_cascade import modele as m, controle as c
        q = Rotation.from_rotvec([.4,-.3,.2]).as_matrix() if spec.get('tourne') else np.eye(3)
        n,joints = m(spec['n'],q)
        control=lambda:c(n.etat_precis(),n.reactions(),spec['n'],joints,q)
    elif spec['famille']=='contact':
        n=Noyau([0.,0.,-9.81]);radius=.1;stiffness=1e5
        for i in range(spec['n']):
            n.corps(str(i),1.,np.eye(3).ravel().tolist(),[i*.3,0.,radius-.005])
            n.liaison(f'l{i}',None,i,bloque_t=[0,1],bloque_r=[0,1,2])
            n.contact(f'c{i}',i,[0.,0.,0.],radius,k=stiffness)
        control=lambda:{'penetration_m':float(np.max(np.abs(np.array(n.etat()[1])[:,2]-(radius-(9.81/stiffness)**(2/3)))))}
    elif spec['famille']=='sans_equilibre':
        n=Noyau([0.,0.,0.]);n.corps('libre',1.,np.eye(3).ravel().tolist(),[0.,0.,0.])
        n.effort(0,[1.,0.,0.],[0.,0.,0.]);control=lambda:{}
    else:raise ValueError(spec)
    def canonique():
        s=n.etat_precis()
        return np.asarray(s[1])+np.asarray(s[6]),np.asarray(s[2]).reshape(-1,3,3)
    return n,[1.],lambda _:None,canonique,control,dict(tol=1e-8,iters=100,paliers_max=1,strict=True)


def run(spec):
    n,loads,charge,canonical,control,options=modele(spec)
    before=n.etat_precis();records=[];failure=None;restored=None
    start=time.perf_counter()
    for load in loads:
        charge(float(load));checkpoint=n.etat_precis()
        try:n.statique(**options)
        except Exception as e:
            failure=f'{type(e).__name__}: {e}';restored=n.etat_precis()==checkpoint
        records.append(n.statique_info())
        if failure:break
    elapsed=time.perf_counter()-start
    pos,rot=canonical()
    return dict(spec=spec,temps_s=elapsed,erreur=failure,restaure=restored,
                rapports=records,controles=control() if failure is None else {},
                position=pos.tolist(),rotation=rot.tolist(),
                residu_statique_max=float(np.max(np.abs(n.residu_statique()))),
                etat_initial=before, etat_final=n.etat_precis())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec');p.add_argument('--echauffement',action='store_true')
    a=p.parse_args();spec=json.loads(a.spec)
    if a.echauffement:run(spec)
    print(json.dumps(run(spec),ensure_ascii=False,allow_nan=False,separators=(',',':')))
