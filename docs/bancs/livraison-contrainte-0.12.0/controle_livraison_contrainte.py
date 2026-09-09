"""Qualification de la roue : six processus frais, 257 fréquences et six charges.

Validation de livraison, sans classement de vitesse ni reprise d'un essai.
Les références Decimal 70/90 déjà produites sont identifiées et relues.
Le candidat emploie uniquement l'API installée ; le juge reste indépendant.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

CI=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ecrire(path,value):
    def encoder(x):
        if hasattr(x,'tolist'):return x.tolist()
        raise TypeError(type(x).__name__)
    with Path(path).open('x') as f:json.dump(value,f,ensure_ascii=False,indent=2,default=encoder,allow_nan=False)


def worker(args):
    import numpy as np
    import scipy
    import vinkulum
    from vinkulum import _vinkulum
    from vinkulum.reduction_contrainte import ReductionContrainte,AssemblageContraint
    from comparaison_masse import BibliothequeComparaison
    from confronte_ports_exudyn import lire
    from modeles_assemblage_complements import charger_cas,espace_juge,empiler
    from juge_masse_couplee import juger_couple
    from experience_controle_facteurs import controle_oracle,qualifier,qualifier_controle
    os.sched_setaffinity(0,{8})
    args.sortie.mkdir(parents=True,exist_ok=False)
    result=dict(n=args.n,mode=args.mode,statut='refus_ou_echec',nombre_frequences=257,charges=6,
                certification_machine_reponses=False,executable=sys.executable,
                version=vinkulum.__version__,numpy=np.__version__,scipy=scipy.__version__,
                plateforme=platform.platform(),cpu=sorted(os.sched_getaffinity(0)),
                module=str(Path(vinkulum.__file__).resolve()),extension_sha256=sha(_vinkulum.__file__),
                parametres=dict(modes_retenus=1,blocs=4,max_directions=128,profondeur=8,
                                alpha_masse=.49,beta_masse=1.51,gamma_hz=80.))
    if vinkulum.__version__!='0.12.0' or 'site-packages' not in Path(vinkulum.__file__).parts:
        raise ValueError('roue 0.12.0 non éditable requise')
    dist=importlib.metadata.distribution('vinkulum')
    result['paquet_sha256']={str(p):sha(dist.locate_file(p)) for p in dist.files
                            if str(p).endswith(('.py','.so','/METADATA','/RECORD'))}
    try:
        source_manifest=json.loads((args.entrees/'manifest.json').read_text())
        if args.mode=='simple':
            entry=args.entrees/f'n{args.n}-f40.npz';reference=args.entrees/f'n{args.n}-f40.ref.npy'
            d,m,h,forces,omega=lire(entry);pieces=[(d,m,h)]
            result['entrees_sha256']={entry.name:sha(entry)}
            if sha(entry)!=source_manifest[entry.name]:raise ValueError('entrée altérée')
        else:
            case=args.entrees/f'n{args.n}';pieces,j=charger_cas(case)
            reference=case/'reference90.npy';d,m=espace_juge(pieces,j)
            h,forces,omega=j['metrique'],j['forces'],j['omega']
            result['entrees_sha256']=json.loads((case/'modele.json').read_text())['fichiers_sha256']
        result['reference_sha256']=sha(reference)
        if sha(reference)!=source_manifest[str(reference.relative_to(args.entrees))]:raise ValueError('référence altérée')
        reference=np.load(reference,mmap_mode='r',allow_pickle=False)
        # Les métriques historiques étaient traitées par Cholesky inférieur.
        # Leur triangle supérieur peut différer par arrondi. Rendre explicite
        # la même métrique symétrique, sans modifier D, M ni les références.
        result['metriques_normalisation']=[]
        normalisees=[]
        for dd,mm,metric in pieces:
            symmetric=np.tril(metric)+np.tril(metric,-1).T
            result['metriques_normalisation'].append(dict(
                convention='recopie_triangle_inferieur_comme_cholesky_historique',
                ecart_max=float(np.max(abs(symmetric-metric))),
                source_sha256=hashlib.sha256(np.asarray(metric,dtype='<f8').tobytes()).hexdigest(),
                symetrique_sha256=hashlib.sha256(np.asarray(symmetric,dtype='<f8').tobytes()).hexdigest()))
            normalisees.append((dd,mm,symmetric))
        start=time.perf_counter();reductions=[]
        for dd,mm,metric in normalisees:
            n=dd.shape[1]
            reductions.append(ReductionContrainte(dd,mm,np.arange(n-6),np.arange(n-6,n),metric,
                float(omega[-1]),gamma=float((2*np.pi*80)**2),alpha_masse=.49,beta_masse=1.51))
        service=reductions[0] if args.mode=='simple' else AssemblageContraint(
            reductions,j['applications'],h,facteur_externe=j['facteur_externe'],masse_externe=j['masse_externe'])
        prepared=time.perf_counter();fields=[];returns=[]
        for w in omega:
            rep=service.reponses(w,forces)
            fields.append(rep['champ'] if args.mode=='simple' else empiler(rep))
            returns.append({k:rep[k] for k in ('marge','bornes','certification_machine','statut_controle')})
        end=time.perf_counter()
        fields=np.array(fields);np.save(args.sortie/'champs.npy',fields,allow_pickle=False)
        result.update(preparation_s=prepared-start,reponses_s=end-prepared,total_s=end-start,
                      champs_sha256=sha(args.sortie/'champs.npy'),retours_par_frequence=returns,
                      certificats=[r.certificats for r in reductions],diagnostics=[r.diagnostic for r in reductions])
        result['juge']=juger_couple(d,m,h,fields,reference,bibliotheque=BibliothequeComparaison(args.bibliotheque))
        result['controle_oracle']=controle_oracle(d,m,fields,reference,returns)
        result['statut']='termine'
        result['qualifie_champs']=qualifier(result)
        result['qualifie_controle']=qualifier_controle(result)
    except Exception as exc:
        result.update(statut='refus_ou_echec',erreur=dict(type=type(exc).__name__,message=str(exc)))
        traceback.print_exc()
    ecrire(args.sortie/'resultat.json',result)
    return 0 if result.get('qualifie_champs') and all(result.get('qualifie_controle',(False,))) else 1


def campagne(args):
    from comparaison_masse import compiler_comparaison
    args.sortie.mkdir(parents=True,exist_ok=False)
    lib=compiler_comparaison(args.sortie/'juge-build')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',RAYON_NUM_THREADS='1')
    ecrire(args.sortie/'protocole.json',dict(cas=[(mode,n) for mode in ('simple','assemblage') for n in (32,128,512)],
        repetitions=1,portee='qualification de livraison ; aucune comparaison de vitesse',
        python=str(args.python),roue_sha256=sha(args.roue),script_sha256=sha(__file__),
        metriques='recopie explicite du triangle inférieur utilisé par le Cholesky historique',
        sources_juge_sha256={p.name:sha(p) for p in CI.glob('*.py')},
        entrees_simple=str(args.simple),entrees_assemblage=str(args.assemblage)))
    rows=[]
    for mode in ('simple','assemblage'):
        for n in (32,128,512):
            name=f'{mode}-n{n}';out=args.sortie/name
            command=[str(args.python),str(Path(__file__).resolve()),'worker','--mode',mode,'--n',str(n),
                     '--entrees',str(args.simple if mode=='simple' else args.assemblage),
                     '--sortie',str(out),'--bibliotheque',str(lib)]
            with (args.sortie/(name+'.log')).open('x') as log:
                run=subprocess.run(command,cwd=args.sortie,env=env,stdout=log,stderr=subprocess.STDOUT)
            ecrire(args.sortie/(name+'.terminal.json'),dict(code=run.returncode,commande=command))
            result=json.loads((out/'resultat.json').read_text()) if (out/'resultat.json').exists() else dict(statut='crash')
            rows.append(dict(dossier=name,code=run.returncode,statut=result['statut'],
                             qualifie_champs=result.get('qualifie_champs'),qualifie_controle=result.get('qualifie_controle')))
            print(rows[-1],flush=True)
    ecrire(args.sortie/'bilan.json',rows)
    ecrire(args.sortie/'manifest.json',{str(p.relative_to(args.sortie)):sha(p) for p in args.sortie.rglob('*') if p.is_file()})
    return int(any(r['code'] for r in rows))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    w=sub.add_parser('worker');w.add_argument('--mode',choices=['simple','assemblage'],required=True);w.add_argument('--n',type=int,required=True)
    for name in ('entrees','sortie','bibliotheque'):w.add_argument('--'+name,type=Path,required=True)
    c=sub.add_parser('campagne')
    for name in ('python','roue','simple','assemblage','sortie'):c.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();raise SystemExit(worker(args) if args.action=='worker' else campagne(args))
