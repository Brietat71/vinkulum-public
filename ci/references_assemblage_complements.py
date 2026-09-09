"""Génération indépendante 70/90 chiffres des trois assemblages physiques."""
import argparse,json,os,shutil,time
from pathlib import Path
import numpy as np
from modeles_assemblage_complements import creer_cas,charger_cas,espace_juge,empiler
from oracle_assemblage_complements import OracleAssemblage
from comparaison_masse import compiler_comparaison,BibliothequeComparaison
from juge_masse_couplee import juger_couple
from experience_controle_facteurs import sha,ecrire,environnement

SOURCES=('references_assemblage_complements.py','modeles_assemblage_complements.py',
    'oracle_assemblage_complements.py','oracle_champs_couples.py','oracle_champs_ports.py',
    'reference_ports_precision.py','comparaison_masse.py','comparaison_masse_native.cpp',
    'inertie_binaire_native.cpp','inertie_binaire_compile.py','inertie_complement_dirigee.py',
    'matrices_certificat.py','racine_masse_juge.py','juge_masse_couplee.py','confronte_ports_exudyn.py',
    'experience_controle_facteurs.py')

def generer(source,out):
    out.mkdir(parents=True,exist_ok=False);(out/'sources').mkdir();os.sched_setaffinity(0,{8})
    for name in SOURCES:shutil.copy2(Path(__file__).with_name(name),out/'sources'/name)
    bib=BibliothequeComparaison(compiler_comparaison(out/'build'));ecrire(out/'environnement.json',environnement('facteurs_controle'))
    for n in (32,128,512):
        case=out/f'n{n}';creer_cas(source,case,n);pieces,j=charger_cas(case);d,m=espace_juge(pieces,j)
        oracles=[OracleAssemblage([(d,m) for d,m,h in pieces],j['applications'],j['facteur_externe'],j['masse_externe'],dps=p) for p in (70,90)]
        fields=[[],[]];start=time.perf_counter()
        for index,w in enumerate(j['omega']):
            for k,o in enumerate(oracles):fields[k].append(empiler(o.reponses(w,j['forces'])))
            if index%64==0:print(n,index,'/257',flush=True)
        arrays=list(map(np.asarray,fields))
        for precision,a in zip((70,90),arrays,strict=True):np.save(case/f'reference{precision}.npy',a,allow_pickle=False)
        judge=juger_couple(d,m,j['metrique'],*arrays,bibliotheque=bib)
        error=max(*judge['maxima'].values(),*judge['maxima_operateurs'].values())
        info=dict(precisions=[70,90],maxima=judge['maxima'],maxima_operateurs=judge['maxima_operateurs'],
            ecart_relatif_max=error,reference_sha256=sha(case/'reference90.npy'),audit_masse=judge['audit_masse'],
            certification_machine=False,temps_hors_mesure_s=time.perf_counter()-start)
        ecrire(case/'reference.json',info)
        if error>1e-8:raise ArithmeticError('références non convergées')
        print('terminé',n,error,info['temps_hors_mesure_s'],flush=True)
    ecrire(out/'manifest.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path);p.add_argument('sortie',type=Path);a=p.parse_args();generer(a.source,a.sortie)
