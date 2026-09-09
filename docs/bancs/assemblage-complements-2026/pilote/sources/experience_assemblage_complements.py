"""24 processus frais : trois branches assemblées, contrôle global et LU corrigée."""
import argparse,json,os,shutil,subprocess,sys,time,warnings
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix,diags,vstack
from scipy.sparse.linalg import splu
from assemblage_complements import AssemblageComplements
from modeles_assemblage_complements import charger_cas,espace_juge,empiler
from comparaison_masse import compiler_comparaison,BibliothequeComparaison
from inertie_binaire_compile import BibliothequeInertie,certifier_inertie_compilee
from ordre_separateurs import ordre_separateurs
from condensation_energie import CondensationEnergie
from krylov_contraint import KrylovContraint
from controle_masse_comparee import ControleMasseComparee
from experience_masse_comparee import selection,SOURCES as PRECEDENTS,PARAMETRES
from experience_controle_facteurs import (sha,ecrire,environnement,qualifier,qualifier_controle,
    controle_oracle,THREADS,FILS_SUPPLEMENTAIRES,BUDGETS_INERTIE)
from confronte_ports_exudyn import normes
from juge_masse_couplee import juger_couple

CI=Path(__file__).resolve().parent
VARIANTES=('assemblage_controle','lu_corrigee')
SOURCES=tuple(dict.fromkeys(PRECEDENTS+('assemblage_complements.py','modeles_assemblage_complements.py',
    'oracle_assemblage_complements.py','references_assemblage_complements.py','experience_assemblage_complements.py')))


def applications_completes(pieces,apps):
    total=6+sum(d.shape[1]-6 for d,m,h in pieces);offset=6;maps=[]
    for (d,m,h),a in zip(pieces,apps,strict=True):
        ni=d.shape[1]-6;rows=list(range(ni));cols=list(range(offset,offset+ni));values=[1.]*ni;offset+=ni
        for i in range(6):
            for j in range(6):
                if a[i,j]!=0:rows.append(ni+i);cols.append(j);values.append(a[i,j])
        maps.append(csr_matrix((values,(rows,cols)),shape=(ni+6,total)))
    return maps,total


def worker(inputs,out,n,variant,library):
    out.mkdir(parents=True,exist_ok=False);os.sched_setaffinity(0,{8})
    case=inputs/f'n{n}';pieces,j=charger_cas(case);refpath=case/'reference90.npy'
    manifest=json.loads((inputs/'manifest.json').read_text())
    if sha(refpath)!=manifest[f'n{n}/reference90.npy']:raise ValueError('référence altérée')
    bib=BibliothequeComparaison(library);inertia=BibliothequeInertie(library)
    info=json.loads((case/'modele.json').read_text());size=sum(d.shape[1] for d,m,h in pieces)+6
    result=dict(n=n,variante=variant,statut='refus_ou_echec',parametres=PARAMETRES,
        nombre_frequences=len(j['omega']),charges=6,certification_machine_reponses=False,
        environnement=environnement('facteurs_controle' if variant=='assemblage_controle' else variant),
        sources_sha256={name:sha(CI/name) for name in SOURCES},modele=info,reference_sha256=sha(refpath),
        bibliotheque=bib.identite,portee='trois branches avec jonction rigide linéarisée ; sans comparaison externe')
    fields=np.full((len(j['omega']),size,6),np.nan);returns=[];controls=[];certs=[];phases={};start=last=time.perf_counter();phase='preparation'
    def mark(name):
        nonlocal last
        now=time.perf_counter();phases[name]=now-last;last=now
    try:
        with warnings.catch_warnings(record=True) as observed:
            warnings.simplefilter('always')
            if variant=='assemblage_controle':
                for index,(d,m,h) in enumerate(pieces):
                    phase=f'piece{index}';total=d.shape[1]
                    qr=CondensationEnergie(d,np.arange(total-6),np.arange(total-6,total),h)
                    mi=m[:-6,:-6].tocsr();phi,b,sel=selection(qr,mi,n+100*index)
                    perm=ordre_separateurs(d[:,:-6],mi,taille_feuille=8)
                    cert=certifier_inertie_compilee(d[:,:-6],mi,b,float((2*np.pi*80)**2),bibliotheque=inertia,permutation=perm,**BUDGETS_INERTIE)
                    reduction=KrylovContraint(qr,m,b,phi,cert['lambda_min'],j['omega'][-1],blocs=4,max_directions=128)
                    control=ControleMasseComparee(reduction,bibliotheque=bib,alpha=.49,beta=1.51,profondeur=8)
                    controls.append(control);certs.append(dict(inertie=cert,comparaison=control.comparaison_masse,selection=sel));mark(phase)
                phase='assemblage';assembly=AssemblageComplements(controls,j['applications'],j['metrique'],
                    facteur_externe=j['facteur_externe'],masse_externe=j['masse_externe']);mark(phase)
                result['taille_conservee']=assembly.taille_conservee
            else:
                maps,total=applications_completes(pieces,j['applications']);dp=[d@p for (d,m,h),p in zip(pieces,maps,strict=True)]
                k=sum(x.T@x for x in dp).tocsc();mass=sum(p.T@m@p for (d,m,h),p in zip(pieces,maps,strict=True)).tocsc()
                pe=csr_matrix((np.ones(6),(np.arange(6),np.arange(6))),shape=(6,total))
                k+=pe.T@csr_matrix(j['facteur_externe'].T@j['facteur_externe'])@pe;mass+=pe.T@csr_matrix(j['masse_externe'])@pe
                eq=1/np.sqrt(k.diagonal());scale=diags(eq);ke,me=(scale@k@scale).tocsc(),(scale@mass@scale).tocsc()
                f=np.zeros((total,6));f[:6]=j['forces'];fe=eq[:,None]*f;mark('assemblage_lu')
            prepared=last;phase='reponses'
            for index,w in enumerate(j['omega']):
                if variant=='assemblage_controle':
                    rep=assembly.reponses(w,j['forces']);fields[index]=empiler(rep)
                    returns.append({k:v for k,v in rep.items() if k not in ('champs','champ_ports','coordonnees','schur','enveloppe')})
                else:
                    factor=splu(ke-w*w*me,permc_spec='MMD_AT_PLUS_A');u=eq[:,None]*factor.solve(fe)
                    for _ in range(2):
                        residual=f.copy()
                        for (d,m,h),p in zip(pieces,maps,strict=True):
                            x=p@u;residual+=p.T@(-d.T@(d@x)+w*w*(m@x))
                        residual[:6]+=-j['facteur_externe'].T@(j['facteur_externe']@u[:6])+w*w*(j['masse_externe']@u[:6])
                        u+=eq[:,None]*factor.solve(eq[:,None]*residual)
                    xs=[p@u for p in maps];fields[index]=np.vstack(xs+[u[:6]])
                    energy=np.zeros((2,6))
                    for (d,m,h),x in zip(pieces,xs,strict=True):
                        nm,nd=normes(d,m,x);energy+=np.array([nm*nm,nd*nd])
                    energy[0]+=np.sum(u[:6]*(j['masse_externe']@u[:6]),axis=0)
                    energy[1]+=np.sum((j['facteur_externe']@u[:6])**2,axis=0);np.sqrt(energy)
            mark(phase);result.update(statut='termine',preparation_s=prepared-start,reponses_s=last-prepared,total_s=last-start,
                avertissements=sorted(set(str(x.message) for x in observed)))
        if not np.all(np.isfinite(fields)):raise ArithmeticError('champs non finis')
    except Exception as exc:
        result.update(statut='refus_ou_echec',phase=phase,erreur=dict(type=type(exc).__name__,message=str(exc),diagnostic=getattr(exc,'diagnostic',None)))
    result.update(phases_s=phases,certificats=certs,retours_par_frequence=returns)
    np.save(out/'champs.npy',fields,allow_pickle=False);result['champs_sha256']=sha(out/'champs.npy')
    for index,c in enumerate(controls):np.savez_compressed(out/f'base{index}.npz',B=c.reduction.b,Phi=c.reduction.phi,base=c.reduction.base)
    if result['statut']=='termine':
        try:
            d,m=espace_juge(pieces,j);ref=np.load(refpath,mmap_mode='r',allow_pickle=False)
            result['juge']=juger_couple(d,m,j['metrique'],fields,ref,bibliotheque=bib)
            if variant=='assemblage_controle':result['controle_oracle']=controle_oracle(d,m,fields,ref,returns)
        except Exception as exc:result.update(statut='refus_ou_echec',phase='juge',erreur=dict(type=type(exc).__name__,message=str(exc)))
    result['qualifie_champs']=qualifier(result)
    if variant=='assemblage_controle':result['qualifie_controle']=qualifier_controle(result)
    ecrire(out/'resultat.json',result)


def campagne(inputs,out,passages):
    out.mkdir(parents=True,exist_ok=False);(out/'sources').mkdir()
    for name in SOURCES:shutil.copy2(CI/name,out/'sources'/name)
    library=compiler_comparaison(out/'build');results=[]
    ecrire(out/'protocole.json',dict(parametres=PARAMETRES,passages=passages,variantes=VARIANTES,
        maillages=[32,128,512],branches=3,python=sys.executable,entrees=str(inputs),cpu=8,threads=THREADS))
    env=dict(os.environ,**THREADS,**FILS_SUPPLEMENTAIRES)
    for passage in range(passages):
        for n in (32,128,512):
            for variant in VARIANTES:
                name=f'n{n}-{variant}-p{passage}';folder=out/name
                command=[sys.executable,str(CI/'experience_assemblage_complements.py'),'worker',str(inputs),str(folder),
                    '--n',str(n),'--variante',variant,'--bibliotheque',str(library)]
                with (out/f'{name}.log').open('x') as log:run=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
                ecrire(out/f'{name}.terminal.json',dict(code=run.returncode,commande=command))
                result=json.loads((folder/'resultat.json').read_text()) if (folder/'resultat.json').exists() else dict(statut='crash')
                result.update(passage=passage,role='chauffe' if passage==0 else 'mesure',dossier=name)
                results.append(result);ecrire(out/'bilan.json',results,remplacer=(out/'bilan.json').exists())
                print(name,result['statut'],result.get('qualifie_champs'),result.get('qualifie_controle'),result.get('total_s'),result.get('erreur'),flush=True)
    ecrire(out/'manifest.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['worker','campagne']);p.add_argument('entrees',type=Path);p.add_argument('sortie',type=Path)
    p.add_argument('--passages',type=int,default=4);p.add_argument('--n',type=int);p.add_argument('--variante',choices=VARIANTES);p.add_argument('--bibliotheque',type=Path);a=p.parse_args()
    if a.mode=='worker':worker(a.entrees,a.sortie,a.n,a.variante,a.bibliotheque)
    else:campagne(a.entrees,a.sortie,a.passages)
