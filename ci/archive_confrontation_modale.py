"""Relecture des spectres et du classement, sans exécutable concurrent.

Les temps et RSS sont des observations conservées, pas des valeurs que ce
vérificateur pourrait démontrer. Il contrôle leur cohérence et interdit de
classer des sorties qui ne satisfont pas le critère numérique commun.
"""
import argparse
import gzip
import json
import math
from pathlib import Path
import statistics

from confronte_modes import (SCHEMA, TAILLES, MOTEURS, SEUIL, REPETITIONS,
                             DELAI, MEMOIRE, ROOT, juger, oracle, sha, ecrire, modele_mbdyn)


def proche(a, b):
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=8e-15)


def verifier_document(d, references=None):
    contrat = dict(schema=SCHEMA, seuil=SEUIL, tailles=list(TAILLES), moteurs=list(MOTEURS),
                   repetitions=REPETITIONS, delai_s=DELAI, memoire_octets=MEMOIRE)
    if any(d.get(k) != v for k, v in contrat.items()):
        raise ValueError('contrat de campagne incorrect')
    refs = references if references is not None else {str(n): oracle(n) for n in TAILLES}
    if set(d['references']) != set(refs):
        raise ValueError('références incomplètes')
    for key, ref in refs.items():
        old = d['references'][key]
        for field in ('frequences', 'bande_hz'):
            if len(old[field]) != len(ref[field]) or not all(proche(x,y) for x,y in zip(old[field],ref[field])):
                raise ValueError('référence indépendante différente')
    expected = {(m,n,r) for m in MOTEURS for n in TAILLES for r in range(REPETITIONS+1)}
    seen = set()
    by_case = {}
    for row in d['essais']:
        key = row['moteur'], row['taille'], row['repetition']
        if key not in expected or key in seen:
            raise ValueError('essai inattendu ou dupliqué')
        seen.add(key)
        if row['echauffement'] != (row['repetition'] == 0):
            raise ValueError('échauffement mal identifié')
        elapsed = row['secondes_processus']
        rss = row['rss_kio']
        if not math.isfinite(elapsed) or elapsed <= 0 or (rss is not None and (type(rss) is not int or rss <= 0)):
            raise ValueError('mesure invalide')
        if 'resultat' in row:
            if row['code_retour'] != 0:
                raise ValueError('résultat associé à un échec de processus')
            if not row['moteur'].startswith('mbdyn'):
                package='vinkulum' if row['moteur']=='vinkulum' else 'exudyn'
                if row['resultat']['version'] != d['identites'][package]['version']:
                    raise ValueError('version de paquet différente')
            key_ref=str(row['taille'])
            # Reproduire les métriques avec leur référence archivée, puis
            # contrôler séparément la décision avec la référence recalculée.
            # Le parallélisme BLAS peut changer ses derniers bits ; ce n'est
            # pas une altération du résultat conservé. Le seuil reste SEUIL.
            result = juger(row['moteur'],row['resultat'],d['references'][key_ref])
            independent = juger(row['moteur'],row['resultat'],refs[key_ref])
            if independent['decision'] != result['decision']:
                raise ValueError('décision dépendant des derniers bits de la référence')
            if result['decision'] != row['decision']:
                raise ValueError('décision spectrale altérée')
            for field in ('erreur_relative_frequence','partie_reelle_relative'):
                if field in result and not proche(result[field],row[field]):
                    raise ValueError('métrique spectrale altérée')
        elif row['decision'] not in ('delai_depasse','echec_execution','sortie_invalide'):
            raise ValueError('décision sans résultat')
        if row['decision']=='delai_depasse' and (elapsed < DELAI*.99 or row['code_retour'] == 0):
            raise ValueError('dépassement de délai incohérent')
        if row['decision']=='echec_execution' and row['code_retour']==0:
            raise ValueError('échec sans code de retour')
        by_case.setdefault(key[:2],[]).append(row)
    if seen != expected:
        raise ValueError('campagne incomplète')
    summary=[]
    for nb in TAILLES:
        for moteur in MOTEURS:
            rows=by_case[moteur,nb]
            accepted=all(r['decision']=='qualifie' for r in rows)
            measured=[r for r in rows if not r['echauffement']]
            record=dict(moteur=moteur,taille=nb,classe=accepted,
                        decisions=[r['decision'] for r in sorted(rows,key=lambda r:r['repetition'])])
            if accepted:
                times=[r['secondes_processus'] for r in measured]
                rss=[r['rss_kio'] for r in measured]
                if any(r is None for r in rss):
                    raise ValueError('RSS manquant pour un essai classé')
                record.update(mediane_s=statistics.median(times),minimum_s=min(times),maximum_s=max(times),
                              rss_max_kio=max(rss),erreur_relative_max=max(r['erreur_relative_frequence'] for r in rows))
            summary.append(record)
    return summary


def archiver(source, destination):
    source=Path(source)
    data=json.loads((source/'campagne.json').read_text())
    summary=verifier_document(data)
    target=Path(destination)
    target.mkdir(parents=True,exist_ok=True)
    if (target/'manifest.json').exists():
        raise ValueError('archive existante')
    # Capturer les sources réellement utilisées : refuser un changement
    # intervenu pendant la campagne, avant de créer les instantanés.
    for name,digest in data['sources'].items():
        p=ROOT/'ci'/name
        if p.name != name or sha(p)!=digest:
            raise ValueError('producteur modifié depuis la campagne')
        (target/name).write_bytes(p.read_bytes())
    (target/'protocole.md').write_bytes((ROOT/'docs/PROTOCOLE_CONFRONTATION_MODALE.md').read_bytes())
    (target/'campagne.json.gz').write_bytes(gzip.compress(json.dumps(data,allow_nan=False).encode(),mtime=0))
    ecrire(target/'bilan.json',summary)
    inputs=target/'modeles';inputs.mkdir()
    for nb in TAILLES:
        for moteur in MOTEURS:
            if moteur.startswith('mbdyn'):
                p=source/f'{moteur}-{nb}-0'/'modele.mbd'
                for repetition in range(1,REPETITIONS+1):
                    if sha(source/f'{moteur}-{nb}-{repetition}'/'modele.mbd')!=sha(p):
                        raise ValueError('modèle changé entre répétitions')
                (inputs/f'{moteur}-{nb}.mbd').write_bytes(p.read_bytes())
    files={str(p.relative_to(target)):sha(p) for p in sorted(target.rglob('*')) if p.is_file()}
    ecrire(target/'manifest.json',dict(schema=SCHEMA,fichiers=files))
    return summary


def verifier_archive(dossier):
    dossier=Path(dossier)
    manifest=json.loads((dossier/'manifest.json').read_text())
    if manifest['schema']!=SCHEMA:
        raise ValueError('schéma incorrect')
    expected={'campagne.json.gz','bilan.json','confronte_modes.py','modeles_modes.py','protocole.md'}
    expected.update(f'modeles/{m}-{n}.mbd' for m in MOTEURS if m.startswith('mbdyn') for n in TAILLES)
    if set(manifest['fichiers'])!=expected:
        raise ValueError('manifeste incomplet')
    for name,digest in manifest['fichiers'].items():
        p=Path(name)
        if p.is_absolute() or '..' in p.parts or sha(dossier/p)!=digest:
            raise ValueError('archive altérée')
    data=json.loads(gzip.decompress((dossier/'campagne.json.gz').read_bytes()))
    for name,digest in data['sources'].items():
        if name not in ('confronte_modes.py','modeles_modes.py') or sha(dossier/name)!=digest:
            raise ValueError('instantané de source différent')
    for nb in TAILLES:
        for moteur in MOTEURS:
            if moteur.startswith('mbdyn'):
                text=(dossier/'modeles'/f'{moteur}-{nb}.mbd').read_text()
                if text!=modele_mbdyn(nb,moteur,data['references'][str(nb)]['bande_hz']):
                    raise ValueError('modèle archivé différent du protocole')
    summary=verifier_document(data)
    old=json.loads((dossier/'bilan.json').read_text())
    if len(old)!=len(summary):
        raise ValueError('bilan incomplet')
    for a,b in zip(old,summary):
        if a.keys()!=b.keys():
            raise ValueError('bilan altéré')
        for key in a:
            if isinstance(a[key],float):
                if not proche(a[key],b[key]):
                    raise ValueError('bilan numérique altéré')
            elif a[key]!=b[key]:
                raise ValueError('classement altéré')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source');p.add_argument('--sortie');p.add_argument('--verifier')
    a=p.parse_args()
    if a.verifier:
        print(json.dumps(verifier_archive(a.verifier),indent=2))
    elif a.source and a.sortie:
        print(json.dumps(archiver(a.source,a.sortie),indent=2))
    else:
        p.error('--verifier ou --source et --sortie requis')
