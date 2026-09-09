"""Masses consistantes : processus frais, contrôle Loewner et LU corrigée.

Une chauffe puis trois mesures, sans reprise. Paramètres déclarés avant
calcul ; tous les refus restent dans le bilan. Aucune comparaison Exudyn.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time
import warnings

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import LinearOperator,eigsh,splu
from comparaison_masse import compiler_comparaison,BibliothequeComparaison
from inertie_binaire_compile import BibliothequeInertie,certifier_inertie_compilee
from ordre_separateurs import ordre_separateurs
from condensation_energie import CondensationEnergie
from krylov_contraint import KrylovContraint
from controle_masse_comparee import ControleMasseComparee
from confronte_ports_exudyn import lire,normes
from juge_masse_couplee import juger_couple
from experience_controle_facteurs import (sha,ecrire,environnement,controle_oracle,
    qualifier,qualifier_controle,SOURCES as ANCIENS,THREADS,FILS_SUPPLEMENTAIRES,BUDGETS_INERTIE)

CI=Path(__file__).resolve().parent
SOURCES=tuple(dict.fromkeys(ANCIENS+('experience_masse_comparee.py','comparaison_masse.py',
    'comparaison_masse_native.cpp','inertie_binaire_native.cpp','inertie_binaire_compile.py',
    'matrices_certificat.py','ordre_separateurs.py','controle_masse_comparee.py',
    'juge_masse_couplee.py','racine_masse_juge.py','modeles_masse_consistante.py',
    'oracle_champs_couples.py')))
VARIANTES=('masse_comparee','lu_corrigee')
PARAMETRES=dict(alpha=.49,beta=1.51,blocs=4,max_directions=128,profondeur=8,
    seuil_hz=80.,eigsh=dict(k=1,ncv=20,tol=1e-11,which='LA'),budget_inertie=BUDGETS_INERTIE)


def selection(qr,mi,n):
    root=np.sqrt(mi.diagonal());sc=diags(1/root);a=(sc@mi@sc).tocsr();calls=[0,0]
    def mult(v):
        calls[0]+=1
        return (qr.produit(np.asarray(v).reshape(-1,1)/root[:,None])/root[:,None]).ravel()
    def inv(v):
        calls[1]+=1
        return (root[:,None]*qr.solve(root[:,None]*np.asarray(v).reshape(-1,1))).ravel()
    k=LinearOperator(mi.shape,matvec=mult,dtype=float)
    ki=LinearOperator(mi.shape,matvec=inv,dtype=float)
    val,vec=eigsh(a,M=k,Minv=ki,v0=np.random.default_rng(107+n).normal(size=len(root)),**PARAMETRES['eigsh'])
    phi=vec/root[:,None];b=mi@phi
    if val[0]<=0 or not np.all(np.isfinite(b)):raise ArithmeticError('direction proposée invalide')
    return phi,b,dict(valeur_inverse=val.tolist(),appels=calls,graine=107+n,
        probleme='L^-1/2 M L^-1/2 v = mu L^-1/2 K_Q L^-1/2 v ; Phi=L^-1/2 v')


def worker(inputs,out,n,variant,lib):
    out.mkdir(parents=True,exist_ok=False);os.sched_setaffinity(0,{8})
    bib=BibliothequeComparaison(lib);inertia=BibliothequeInertie(lib)
    entry=inputs/f'n{n}-f40.npz';refpath=inputs/f'n{n}-f40.ref.npy'
    manifest=json.loads((inputs/'manifest.json').read_text())
    for path in (entry,refpath):
        if sha(path)!=manifest[path.name]:raise ValueError('entrée/référence altérée')
    d,m,metric,forces,omega=lire(entry);size=d.shape[1]
    result=dict(n=n,variante=variant,statut='refus_ou_echec',parametres=PARAMETRES,
        environnement=environnement('facteurs_controle' if variant=='masse_comparee' else variant),
        sources_sha256={x:sha(CI/x) for x in SOURCES},entree_sha256=sha(entry),reference_sha256=sha(refpath),
        bibliotheque=bib.identite,nombre_frequences=len(omega),charges=6,certification_machine_reponses=False)
    fields=np.full((len(omega),size,6),np.nan);returns=[];phases={};control=None;reduction=None
    start=last=time.perf_counter();phase='preparation'
    def mark(name):
        nonlocal last
        now=time.perf_counter();phases[name]=now-last;last=now
    try:
        with warnings.catch_warnings(record=True) as observed:
            warnings.simplefilter('always')
            if variant=='masse_comparee':
                qr=CondensationEnergie(d,np.arange(size-6),np.arange(size-6,size),metric);mark('condensation')
                phase='selection';mi=m[qr.i][:,qr.i].tocsr();phi,b,info=selection(qr,mi,n);mark(phase)
                result['selection']=info;phase='inertie'
                perm=ordre_separateurs(d[:,qr.i],mi,taille_feuille=8)
                cert=certifier_inertie_compilee(d[:,qr.i],mi,b,float((2*np.pi*80)**2),
                    bibliotheque=inertia,permutation=perm,**BUDGETS_INERTIE)
                result['certificat_complement']=cert;mark(phase);phase='krylov'
                reduction=KrylovContraint(qr,m,b,phi,cert['lambda_min'],omega[-1],blocs=4,max_directions=128)
                mark(phase);phase='controle'
                control=ControleMasseComparee(reduction,bibliotheque=bib,alpha=.49,beta=1.51,profondeur=8)
                mark(phase);result['comparaison_masse']=control.comparaison_masse
                result['directions']=reduction.taille_complement_reduit
            else:
                k=(d.T@d).tocsc();eq=1/np.sqrt(k.diagonal());s=diags(eq)
                ke,me=(s@k@s).tocsc(),(s@m@s).tocsc()
                f=np.zeros((size,6));f[-6:]=forces;fe=eq[:,None]*f;mark('preparation_lu')
            prepared=last;phase='reponses'
            for j,w in enumerate(omega):
                if control is not None:
                    rep=control.reponses(w,forces);fields[j]=rep['champ']
                    returns.append({k:v for k,v in rep.items() if k not in ('champ','coordonnees')})
                else:
                    factor=splu(ke-w*w*me,permc_spec='MMD_AT_PLUS_A');x=eq[:,None]*factor.solve(fe)
                    for _ in range(2):
                        residual=f-d.T@(d@x)+w*w*(m@x)
                        x+=eq[:,None]*factor.solve(eq[:,None]*residual)
                    fields[j]=x;normes(d,m,x)
            mark(phase)
            result.update(statut='termine',preparation_s=prepared-start,reponses_s=last-prepared,total_s=last-start,
                avertissements=sorted(set(str(w.message) for w in observed)))
        if not np.all(np.isfinite(fields)):raise ArithmeticError('champ non fini')
    except Exception as exc:
        result.update(statut='refus_ou_echec',phase=phase,erreur=dict(type=type(exc).__name__,message=str(exc),
            diagnostic=getattr(exc,'diagnostic',None)))
    result['phases_s']=phases
    np.save(out/'champs.npy',fields,allow_pickle=False);result['champs_sha256']=sha(out/'champs.npy')
    result['retours_par_frequence']=returns
    if reduction is not None:
        np.savez_compressed(out/'bases.npz',B=b,Phi=phi,base=reduction.base)
    if result['statut']=='termine':
        try:
            ref=np.load(refpath,mmap_mode='r',allow_pickle=False)
            result['juge']=juger_couple(d,m,metric,fields,ref,bibliotheque=bib)
            if control is not None:result['controle_oracle']=controle_oracle(d,m,fields,ref,returns)
        except Exception as exc:
            result.update(statut='refus_ou_echec',phase='juge',erreur=dict(type=type(exc).__name__,message=str(exc)))
    result['qualifie_champs']=qualifier(result)
    if control is not None:result['qualifie_controle']=qualifier_controle(result)
    ecrire(out/'resultat.json',result)


def campagne(inputs,out,passages):
    out.mkdir(parents=True,exist_ok=False);(out/'sources').mkdir()
    for name in SOURCES:shutil.copy2(CI/name,out/'sources'/name)
    library=compiler_comparaison(out/'build');results=[]
    ecrire(out/'protocole.json',dict(parametres=PARAMETRES,passages=passages,variantes=VARIANTES,
        maillages=[32,128,512],python=sys.executable,entrees=str(inputs),cpu=8,
        threads=THREADS,certification_machine_reponses=False,portee='nouvelle masse P1 ; aucune comparaison inter-modèle'))
    env=dict(os.environ,**THREADS,**FILS_SUPPLEMENTAIRES)
    for passage in range(passages):
        for n in (32,128,512):
            for variant in VARIANTES:
                name=f'n{n}-{variant}-p{passage}';folder=out/name
                command=[sys.executable,str(CI/'experience_masse_comparee.py'),'worker',str(inputs),str(folder),
                         '--n',str(n),'--variante',variant,'--bibliotheque',str(library)]
                with (out/f'{name}.log').open('x') as log:
                    run=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
                ecrire(out/f'{name}.terminal.json',dict(code=run.returncode,commande=command))
                result=json.loads((folder/'resultat.json').read_text()) if (folder/'resultat.json').exists() else dict(statut='crash')
                result.update(passage=passage,role='chauffe' if passage==0 else 'mesure',dossier=name)
                results.append(result);ecrire(out/'bilan.json',results,remplacer=(out/'bilan.json').exists())
                print(name,result['statut'],result.get('qualifie_champs'),result.get('qualifie_controle'),
                    result.get('total_s'),result.get('erreur'),flush=True)
    ecrire(out/'manifest.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['worker','campagne'])
    p.add_argument('entrees',type=Path);p.add_argument('sortie',type=Path);p.add_argument('--passages',type=int,default=4)
    p.add_argument('--n',type=int);p.add_argument('--variante',choices=VARIANTES);p.add_argument('--bibliotheque',type=Path)
    a=p.parse_args()
    if a.mode=='worker':worker(a.entrees,a.sortie,a.n,a.variante,a.bibliotheque)
    else:campagne(a.entrees,a.sortie,a.passages)
