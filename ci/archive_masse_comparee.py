"""Archive compacte, requalification des 24 essais et rejeu des preuves.

Les NPY sont exclus et identifiés ; les champs ne sont pas rejugés. Les
sources et bibliothèques capturées ne sont jamais exécutées par le rejeu.
"""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tempfile

import numpy as np
from archive_inertie_binaire import inventaire,octets,lire_json,preuve_stable
from comparaison_masse import compiler_comparaison,BibliothequeComparaison,certifier_comparaison_masse
from inertie_binaire_compile import BibliothequeInertie,certifier_inertie_compilee
from experience_masse_comparee import PARAMETRES,VARIANTES,SOURCES,CI
from experience_controle_facteurs import qualifier,qualifier_controle,BUDGETS_INERTIE,THREADS,FILS_SUPPLEMENTAIRES
from confronte_ports_exudyn import lire

ARCHIVE=CI.parent/'docs/bancs/masse-comparee-2026'

def digest(data):return hashlib.sha256(data).hexdigest()


def archiver(destination,groupes):
    destination.mkdir(parents=True,exist_ok=False);files={};excluded={}
    for prefix,source in groupes.items():
        source=Path(source)
        paths=sorted(p for p in source.rglob('*') if p.is_file()) if source.is_dir() else [source]
        for p in paths:
            name=prefix+'/'+str(p.relative_to(source)) if source.is_dir() else prefix
            data=p.read_bytes();meta=dict(source_sha256=digest(data),source_octets=len(data))
            if p.suffix=='.npy':excluded[name]=meta;continue
            compress=p.suffix in ('.json','.log','.so') and len(data)>4096
            target=destination/(name+'.gz' if compress else name);target.parent.mkdir(parents=True,exist_ok=True)
            packed=gzip.compress(data,mtime=0) if compress else data;target.write_bytes(packed)
            files[str(target.relative_to(destination))]=dict(**meta,compression='gzip' if compress else None,
                sha256=digest(packed),octets=len(packed))
    (destination/'manifest.json').write_text(json.dumps(dict(schema=1,fichiers=files,exclus_npy=excluded),indent=2)+'\n')


def verifier_rapport(rows,passages=4):
    expected=[(p,n,v) for p in range(passages) for n in (32,128,512) for v in VARIANTES]
    if [(r['passage'],r['n'],r['variante']) for r in rows]!=expected:raise ValueError('plan incomplet ou dupliqué')
    admitted=controls=0;envs=set();sources=None;binary=None
    for r in rows:
        if r['parametres']!=PARAMETRES:raise ValueError('paramètres modifiés')
        if r['role']!=('chauffe' if r['passage']==0 else 'mesure'):raise ValueError('rôle incorrect')
        if sources is None:sources=r['sources_sha256'];binary=r['bibliotheque']['sha256']
        if r['sources_sha256']!=sources or r['bibliotheque']['sha256']!=binary:raise ValueError('sources ou binaire instables')
        e=r['environnement']
        if e['cpu']!=[8] or e['fils']!=THREADS or e['fils_supplementaires']!=FILS_SUPPLEMENTAIRES:
            raise ValueError('budget de calcul différent')
        envs.add(json.dumps({k:e[k] for k in ('python','numpy','scipy','executable')},sort_keys=True))
        if r['certification_machine_reponses'] is not False:raise ValueError('certification non justifiée')
        q=qualifier(r)
        if q!=r['qualifie_champs']:raise ValueError('qualification de champ incohérente')
        admitted+=q
        if r['variante']=='masse_comparee':
            qc=qualifier_controle(r)
            if list(qc)!=r.get('qualifie_controle'):raise ValueError('qualification de contrôle incohérente')
            controls+=qc[0]
    if len(envs)!=1:raise ValueError('environnements numériques différents')
    return dict(essais=len(rows),champs_admis=admitted,controles_admis=controls)


def verifier_archive(racine=ARCHIVE,*,rejouer=True):
    racine=Path(racine);manifest=inventaire(racine)
    for group in ('campagne','references','pilote'):
        original=lire_json(racine,group+'/manifest.json')
        for name,expected in original.items():
            full=group+'/'+name
            actual=manifest['exclus_npy'][full]['source_sha256'] if full in manifest['exclus_npy'] else digest(octets(racine,full))
            if actual!=expected:raise ValueError('objet différent du manifeste original')
    rows=lire_json(racine,'campagne/bilan.json');result=verifier_rapport(rows)
    verifier_rapport(lire_json(racine,'pilote/bilan.json'),passages=1)
    cache={}
    for r in rows:
        prefix='campagne/'+r['dossier'];single=lire_json(racine,prefix+'/resultat.json')
        if single!={k:v for k,v in r.items() if k not in ('passage','role','dossier')}:raise ValueError('résultat individuel différent')
        if lire_json(racine,prefix+'.terminal.json')['code']!=0:raise ValueError('processus en échec')
        if manifest['exclus_npy'][prefix+'/champs.npy']['source_sha256']!=r['champs_sha256']:raise ValueError('champ exclu incohérent')
        for name,expected in r['sources_sha256'].items():
            if digest(octets(racine,'campagne/sources/'+name))!=expected:raise ValueError('source capturée différente')
            if digest((CI/name).read_bytes())!=expected:raise ValueError('source courante différente de la mesure')
        for suffix,key in (('.npz','entree_sha256'),('.ref.npy','reference_sha256')):
            name=f'references/n{r["n"]}-f40{suffix}'
            value=manifest['exclus_npy'][name]['source_sha256'] if name in manifest['exclus_npy'] else digest(octets(racine,name))
            if value!=r[key]:raise ValueError('modèle ou référence différents')
        identity=lire_json(racine,'campagne/build/identite.json')
        if r['bibliotheque']['sha256']!=identity['bibliotheque_sha256'] or digest(octets(racine,'campagne/build/comparaison_masse.so'))!=identity['bibliotheque_sha256']:
            raise ValueError('identité de compilation différente')
    if rejouer:
        with tempfile.TemporaryDirectory(prefix='vinkulum-rejeu-masse-') as tmp:
            lib=compiler_comparaison(Path(tmp)/'build');bib=BibliothequeComparaison(lib);inertia=BibliothequeInertie(lib)
            for r in rows:
                if r['variante']!='masse_comparee' or r['statut']!='termine':continue
                n=r['n'];raw=octets(racine,'campagne/'+r['dossier']+'/bases.npz')
                with np.load(io.BytesIO(raw),allow_pickle=False) as z:b=z['B']
                key=(n,digest(b.tobytes()))
                if key not in cache:
                    d,m,*_=lire(racine/f'references/n{n}-f40.npz');mi=m[:-6,:-6].tocsr()
                    comp=certifier_comparaison_masse(mi,.49,1.51,bibliotheque=bib)
                    cert=r['certificat_complement']
                    actual=certifier_inertie_compilee(d[:,:-6],mi,b,float((2*np.pi*80)**2),bibliotheque=inertia,
                        permutation=cert['permutation_physique'],**BUDGETS_INERTIE)
                    cache[key]=(preuve_stable(comp),preuve_stable(actual))
                if cache[key]!=(preuve_stable(r['comparaison_masse']),preuve_stable(r['certificat_complement'])):
                    raise ValueError('preuve différente du rejeu')
    return dict(**result,paires_certificats_rejouees=len(cache))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verifier',type=Path,default=ARCHIVE)
    print(verifier_archive(p.parse_args().verifier))
