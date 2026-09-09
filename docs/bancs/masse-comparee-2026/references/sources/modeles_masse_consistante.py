"""Masses à interpolation linéaire des vitesses de section, assemblage complet.

La rigidité matérielle native est conservée. La masse est une nouvelle
forme cinétique de poutre, avec couplages entre sections voisines. Ce
modèle exige ses propres références ; celles de la masse concentrée ne
sont jamais réutilisées pour juger les champs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.sparse import diags,kron

from confronte_ports_exudyn import lire,sparse_pack
from juge_masse_couplee import juger_couple
from comparaison_masse import compiler_comparaison,BibliothequeComparaison
from experience_controle_facteurs import DONNEES_FIGEES,sha,environnement,jsonable
from oracle_champs_couples import OracleChampsCouples


def masse_consistante(n,modele):
    length,width,height,rho=(float(modele[k]) for k in ('longueur_m','largeur_y_m','hauteur_z_m','rho_kg_m3'))
    step=length/n;element_mass=rho*width*height*step
    inertia=np.array([(width**2+height**2)/12,height**2/12,width**2/12])
    weights=element_mass*np.r_[np.ones(3),inertia]
    # Intégrales de N_i N_j sur chaque segment : h/6 [[2,1],[1,2]].
    # Les six coordonnées encastrées sont supprimées après assemblage.
    diagonal=np.full(n,2/3);diagonal[-1]=1/3
    chain=diags([np.full(n-1,1/6),diagonal,np.full(n-1,1/6)],[-1,0,1],shape=(n,n),format='csr')
    mass=kron(chain,diags(weights),format='csr');mass.eliminate_zeros();mass.sort_indices()
    return mass,dict(type='vitesses_sections_P1',masse_segment_kg=element_mass,
        moments_quadratiques_section_m2=inertia.tolist(),coefficients_segment='[[2,1],[1,2]]/6',
        couplages_conserves=True,inertie_longitudinale_segment_exclue=True)


def generer(source,destination):
    source,destination=Path(source),Path(destination);destination.mkdir(parents=True,exist_ok=False)
    captures=('modeles_masse_consistante.py','oracle_champs_couples.py','oracle_champs_ports.py',
              'reference_ports_precision.py','confronte_ports_exudyn.py','experience_controle_facteurs.py',
              'juge_masse_couplee.py','racine_masse_juge.py','comparaison_masse.py','comparaison_masse_native.cpp',
              'inertie_binaire_native.cpp','inertie_binaire_compile.py','inertie_complement_dirigee.py','matrices_certificat.py')
    (destination/'sources').mkdir()
    for name in captures:shutil.copy2(Path(__file__).with_name(name),destination/'sources'/name)
    bibliotheque=BibliothequeComparaison(compiler_comparaison(destination/'build'))
    env=environnement('facteurs_controle')
    (destination/'environnement.json').write_text(json.dumps(jsonable(env),ensure_ascii=False,indent=2))
    results=[]
    for n in (32,128,512):
        original=source/f'n{n}-f40.npz'
        if sha(original)!=DONNEES_FIGEES[n][0]:raise ValueError('rigidité native de provenance différente')
        metadata=json.loads(original.with_name(original.name+'.json').read_text())
        d,_,metric,forces,omega=lire(original);m,kinetic=masse_consistante(n,metadata['modele'])
        target=destination/f'n{n}-f40.npz'
        np.savez_compressed(target,**sparse_pack('d',d),**sparse_pack('m',m),metrique=metric,forces=forces,omega=omega)
        metadata.update(masse=kinetic,source_rigidite=dict(fichier=str(original),sha256=sha(original)),
                        entree_sha256=sha(target),type_modele='poutre_masse_consistante_P1')
        metadata['modele']['masse_discrete']='Masse consistante P1 des vitesses de section ; voir masse et formule cinétique.'
        for obsolete in ('masse_noeud_interieur_kg','inerties_noeud_interieur_kg_m2'):
            metadata['modele'].pop(obsolete,None)
        (destination/f'n{n}-f40.npz.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2))
        start=time.perf_counter();refs=[OracleChampsCouples(d,m,p=6,dps=p) for p in (70,90)]
        fields=[[],[]]
        for i,w in enumerate(omega):
            for k,ref in enumerate(refs):fields[k].append(ref.reponse(w,forces))
            if i%32==0:print('oracle',n,i,'/257',flush=True)
        arrays=list(map(np.asarray,fields))
        # Préserver les deux calculs, y compris si le juge refuse ensuite.
        reference=destination/f'n{n}-f40.ref.npy';np.save(reference,arrays[1],allow_pickle=False)
        np.save(destination/f'n{n}-f40.ref70.npy',arrays[0],allow_pickle=False)
        judgement=juger_couple(d,m,metric,*arrays,bibliotheque=bibliotheque)
        error=max(*judgement['maxima'].values(),*judgement['maxima_operateurs'].values())
        if error>1e-8:raise ArithmeticError('oracles 70/90 non convergés')
        info=dict(precision=[70,90],ecart_relatif_max=error,maxima=judgement['maxima'],
            maxima_operateurs=judgement['maxima_operateurs'],fichier_sha256=sha(reference),
            temps_hors_mesures_s=time.perf_counter()-start,oracle='OracleChampsCouples',
            controle='colonnes et opérateurs masse/déformation/port après conversion binary64',
            certification_machine=False,compteurs=[x.compteurs for x in refs],audit_masse=judgement['audit_masse'])
        (destination/f'n{n}-f40.ref.npy.json').write_text(json.dumps(info,ensure_ascii=False,indent=2))
        results.append(dict(n=n,entree_sha256=sha(target),reference=info))
        print('oracle terminé',n,error,info['temps_hors_mesures_s'],flush=True)
        (destination/'bilan.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
    manifest={str(p.relative_to(destination)):sha(p) for p in destination.rglob('*') if p.is_file()}
    (destination/'manifest.json').write_text(json.dumps(manifest,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('destination');a=p.parse_args()
    os.sched_setaffinity(0,{8});generer(a.source,a.destination)
