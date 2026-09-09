"""Vérifie l'intégrité et les décisions du dossier 0.14.1, bibliothèque standard.

Ce contrôle d'archive n'est pas une nouvelle exécution des modèles et n'est
pas une preuve d'erreur des trajectoires. Les limites restent dans le dossier.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path


def exige(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def controle_physique(p):
    familles = p['familles']
    exige(len(familles) == 10 and len({f['famille'] for f in familles}) == 10, 'familles incomplètes')
    complet = all(f['code_sortie'] == 0 and f['erreur'] is None for f in familles)
    exige(p['execution_complete'] == complet, 'statut d’exécution erroné')
    lignes = [r for f in familles for r in f['lignes']]
    comptes = {s: sum(r['statut'] == s for r in lignes) for s in ('ok', 'ecart', 'erreur_execution')}
    exige(comptes == p['comptes'], 'comptes de validation erronés')
    exige(sum(comptes.values()) == len(lignes), 'statut physique inconnu')
    index = {(r['cas'], r['source'], r['grandeur']): r for r in lignes}
    exige(len(index) == len(lignes), 'identité physique dupliquée')
    for r in lignes:
        if r['statut'] == 'erreur_execution' or r['nature'] == 'banc_assertions':
            continue
        exige(r['tolerance_relative'] is not None, 'tolérance absente')
        exige(all(math.isfinite(r[k]) for k in ('reference', 'valeur', 'ecart_relatif', 'tolerance_relative')),
              'mesure non finie')
        exige(r['tolerance_relative'] >= 0, 'tolérance négative')
        e = (r['valeur']-r['reference'])/abs(r['reference'])
        exige(e == r['ecart_relatif'], 'écart physique incorrect')
        exige((abs(e) <= r['tolerance_relative']) == (r['statut'] == 'ok'), 'verdict physique incorrect')
    exige(len(p['suivi_historique']) == 62, 'inventaire historique incomplet')
    for s in p['suivi_historique']:
        h = s['historique']
        r = index.get((h['cas'], h['source'], h['grandeur']))
        exige(s['actuel'] == r and s['suivi'] == ('absent' if r is None else r['statut']),
              'rapprochement historique incorrect')
    exige(p['historique_entier_rejoue'] == all(s['actuel'] is not None for s in p['suivi_historique']),
          'couverture historique erronée')
    exige(p['certification_mathematique'] is False, 'qualification physique présentée comme preuve')
    return index


def violations_lisses(d):
    c, t = d['criteres'], d['toupie']['cas']
    exige(len(t) == 32 and len({(r['echelle'], r['cadre'], r['pas']) for r in t}) == 32, 'toupies manquantes')
    n = 0
    for r in t:
        n += r['moment_relatif'] > c['toupie_moment_relatif_max']
        n += r['energie_relative'] > c['toupie_energie_relative_max']
        n += r['pas'] == .001 and r['erreur_omega'] > c['toupie_erreur_omega_max']
    for a, b in zip(t[::2], t[1::2], strict=True):
        n += math.log2(a['erreur_omega']/b['erreur_omega']) < c['toupie_ordre_min']
    exige(len(d['chute']) == 12, 'chutes manquantes')
    n += sum(r['erreur'] > c['chute_erreur_max'] for r in d['chute'])
    exige(len(d['double_pendule']) == len(d['quatre_barres']) == 3, 'raffinement manquant')
    for a, b in zip(d['double_pendule'][:-1], d['double_pendule'][1:]):
        n += math.log2(a['erreur']/b['erreur']) < c['double_pendule_ordre_min']
    n += sum(r['erreur'] > c['quatre_barres_ecart_max'] for r in d['quatre_barres'])
    a = d['assemblage']
    exige(len(a['projections']) == 2304 and len(a['pendules']) == 6, 'assemblage incomplet')
    n += not a['admis']
    exige(len(d['echecs']) == n and d['admis'] == (n == 0), 'décision du corpus lisse erronée')
    exige(not d['certification_trajectoire'], 'comparaison lisse présentée comme certificat')
    return n


def verifie(manifeste):
    m = json.loads(manifeste.read_text())
    data = (manifeste.parent/m['archive']['nom']).read_bytes()
    exige(sha(data) == m['archive']['sha256'], 'empreinte archive incorrecte')
    files = json.loads(gzip.decompress(data))['fichiers']
    exige(set(files) == set(m['archive']['fichiers_sha256']), 'inventaire d’archive incorrect')
    for k, v in files.items():
        exige(sha(v.encode()) == m['archive']['fichiers_sha256'][k], f'empreinte incorrecte : {k}')
    read = lambda n: json.loads(files[n])
    base, courant = read('physique-reference.json'), read('physique-candidat.json')
    ib, ic = controle_physique(base), controle_physique(courant)
    exige(base['execution_complete'] and courant['execution_complete'], 'campagne incomplète')
    exige(base['historique_entier_rejoue'] and courant['historique_entier_rejoue'], 'histoire non rejouée')
    exige(ib.keys() == ic.keys(), 'changement de corpus')
    for k in ib:
        for champ in ('reference', 'tolerance_relative', 'nature'):
            exige(ib[k][champ] == ic[k][champ], f'critère changé : {k}, {champ}')
    l0, l1 = read('lisse-reference.json'), read('lisse-candidat.json')
    exige(l0['criteres'] == l1['criteres'], 'critères lisses modifiés')
    exige(violations_lisses(l0) == 32 and violations_lisses(l1) == 0, 'contre-épreuve mécanique non reproduite')
    for role, rapport, lisse in [('reference', base, l0), ('candidat', courant, l1)]:
        roue = m['roues'][role]
        exige(roue['sha256'] == rapport['environnement']['roue_sha256'], 'mauvaise roue de campagne')
        exige(roue['fichiers_sha256'] == rapport['environnement']['fichiers_sha256'], 'mauvais paquet de campagne')
        exige(lisse['extension_sha256'] == rapport['environnement']['extension_sha256'], 'mauvais noyau lisse')
    exige(m['roues']['reference']['fichiers_sha256'] == m['roues']['reconstruction']['fichiers_sha256'],
          'reconstruction de référence différente')
    exige(m['roues']['candidat']['fichiers_sha256'] == m['roues']['livraison']['fichiers_sha256'],
          'paquet livré différent du paquet qualifié')
    exige(not m['noyau_entier_certifie'] and not m['preuve_formelle_implementation'], 'portée excessive')
    for name, total in [('ci-reference.log', 175), ('ci-candidat.log', 181)]:
        log = files[name]
        for phrase in ('87 passed; 0 failed', '8 passed; 0 failed', f'Ran {total} tests',
                       '41/41', '46/46', '9/9', '╚ CI locale OK'):
            exige(phrase in log, f'CI incomplète : {name}, {phrase}')
    exige('Ran 181 tests' in files['tests-roue-finale.log'] and '\nOK\n' in files['tests-roue-finale.log'],
          'tests roue finale incomplets')
    for k in ('reference', 'candidat'):
        r = read(f'reperes-{k}.json')
        exige(r['anomalie_long_terme_close'] is False and not r['amplification_arrondis_certifiee'],
              'diagnostic de repères présenté comme preuve')
    print('Dossier vérifié : deux campagnes physiques, corpus lisse avant/après, identités des paquets et CI.')
    return m


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifeste', type=Path)
    verifie(p.parse_args().manifeste)
