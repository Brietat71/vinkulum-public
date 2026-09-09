"""Épreuves indépendantes de l'assemblage : covariance et pendule physique.

Les seuils sont fixés dans CRITERES, avant exécution. Le programme ne lit
aucune matrice ni aucun résidu interne pour construire ses références.
Usage : python ci/audit_assemblage_physique.py --sortie bilan.json
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

import vinkulum
from vinkulum import Noyau


CRITERES = {
    'projection_erreur_absolue_max': 1e-10,
    'longueurs_projection': [2.**i for i in (-20, -10, 0, 10, 16, 26)],
    'longueurs_exploration': [2.**40, 2.**48],
    'pendule_erreur_adimensionnelle_max_256_pas': 2e-4,
    'pendule_derive_energie_relative_max_256_pas': 2e-4,
    'pendule_rapport_erreurs_128_sur_256_min': 3.,
    'pendule_rapport_erreurs_128_sur_256_max': 5.,
    'reference_ecart_entre_tolerances_max': 1e-10,
}


def reperes():
    """Les 24 rotations propres à coefficients 0, ±1, exacts en binary64."""
    for permutation in itertools.permutations(range(3)):
        for signes in itertools.product((-1., 1.), repeat=3):
            q = np.eye(3)[:, permutation] @ np.diag(signes)
            if np.linalg.det(q) > 0:
                yield q


def projection(q, longueur, origine, copies, inverse, rotation_bloquee=True):
    n = Noyau([0., 0., 0.])
    n.corps('mobile', 1., np.eye(3).ravel().tolist(), origine.tolist(),
            rot=q.ravel().tolist(), v=(q @ [2., 1., 3.]).tolist(),
            w=(q @ [.1, .2, 0.]).tolist())
    types = ([True, False] if inverse else [False, True]) if rotation_bloquee else [False]
    for copie in range(copies):
        for rotation in types:
            n.liaison(f'{copie}_{rotation}', None, 0,
                      pa=(origine + q @ [longueur, 0., 0.]).tolist(),
                      ra=q.ravel().tolist(), bloque_t=[] if rotation else [1],
                      bloque_r=[2] if rotation else [])
    avant = n.etat_precis()
    try:
        n.assemble()
    except ValueError as e:
        return {'statut': 'refus', 'raison': str(e), 'etat_restitue': n.etat_precis() == avant}
    v, w = q.T @ n.etat()[3][0], q.T @ n.etat()[4][0]
    # Solution fermée de vy + L*wz = 0, avec ou sans wz = 0.
    vy = 0. if rotation_bloquee else longueur**2/(1.+longueur**2)
    wz = 0. if rotation_bloquee else -longueur/(1.+longueur**2)
    erreur = float(max(np.max(np.abs(v - [2., vy, 3.])),
                       np.max(np.abs(w - [.1, .2, wz]))))
    return {'statut': 'succes', 'erreur_reference': erreur,
            'vy': float(v[1]), 'wz': float(w[2])}


def pendule(longueur, nb):
    # Solide de masse 1 kg, Jzz = .2*m*L², articulé à distance L.
    # Son équation exacte réduite est theta'' + sin(theta) = 0,
    # avec tau = t*sqrt(g/(1.2*L)). La référence ne fait pas appel à Vinkulum.
    g, alpha = 9.81, .2
    pulsation = math.sqrt(g/((1.+alpha)*longueur))
    n = Noyau([0., -g, 0.])
    n.corps('pendule', 1., (alpha*longueur**2*np.eye(3)).ravel().tolist(),
            [0., -longueur, 0.], v=[.2*longueur*pulsation, 0., 0.],
            w=[0., 0., .3*pulsation])
    n.liaison('pivot', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1])
    c, s = math.cos(.4), math.sin(.4)
    rot = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    etat = n.etat_precis()
    etat[1][0] = (rot @ [0., -longueur, 0.] + [.01*longueur, 0., 0.]).tolist()
    etat[2][0] = rot.ravel().tolist()
    n.pose_etat(*etat)
    n.assemble()
    r0 = np.array(n.etat()[2][0]).reshape(3, 3)
    initial = [math.atan2(r0[1, 0], r0[0, 0]), n.etat()[4][0][2]/pulsation]
    solutions = [solve_ivp(lambda t, y: [y[1], -math.sin(y[0])], (0., 4.), initial,
                           method='DOP853', rtol=tol, atol=tol/100)
                 for tol in (1e-11, 1e-13)]
    if not all(s.success and s.t[-1] == 4. for s in solutions):
        raise RuntimeError('la référence indépendante ne termine pas son intervalle')
    refs = [s.y[:, -1] for s in solutions]
    ecart_reference = float(np.max(np.abs(refs[0]-refs[1])))
    trajectoire = n.simule(4./pulsation, 4./(nb*pulsation), rho=.9, tol=1e-11)
    if n.t() != 4./pulsation or len(trajectoire) != nb:
        raise RuntimeError('la trajectoire ne suit pas la grille de convergence annoncée')
    r = np.array(n.etat()[2][0]).reshape(3, 3)
    final = np.array([math.atan2(r[1, 0], r[0, 0]), n.etat()[4][0][2]/pulsation])
    energie = lambda y: .5*y[1]**2 + 1.-math.cos(y[0])
    derive = abs(energie(final)-energie(initial))/energie(initial)
    return {'longueur': longueur, 'pas': nb,
            'pas_acceptes': len(trajectoire), 'temps_final': n.t(),
            'initial_adimensionnel': initial, 'final_adimensionnel': final.tolist(),
            'reference': refs[1].tolist(), 'ecart_references': ecart_reference,
            'erreur_reference': float(np.max(np.abs(final-refs[1]))),
            'derive_energie_relative': float(derive)}


def campagne():
    sorties, explorations, trajectoires, echecs = [], [], [], []
    for ir, q in enumerate(reperes()):
        for longueur, decalage, copies, inverse, bloque in itertools.product(
                CRITERES['longueurs_projection'], (0., 2.**20), (1, 3), (False, True), (False, True)):
            cas = {'repere': ir, 'longueur': longueur, 'origine': decalage,
                   'copies': copies, 'ordre_inverse': inverse, 'rotation_bloquee': bloque}
            cas.update(projection(q, longueur, np.full(3, decalage), copies, inverse, bloque))
            sorties.append(cas)
            if cas['statut'] != 'succes' or cas['erreur_reference'] > CRITERES['projection_erreur_absolue_max']:
                echecs.append({'projection': len(sorties)-1})
    for longueur in CRITERES['longueurs_exploration']:
        cas = {'longueur': longueur}
        cas.update(projection(np.eye(3), longueur, np.zeros(3), 1, False))
        explorations.append(cas)
        # Hors du domaine principal, un refus est admis s'il est atomique.
        if ((cas['statut'] == 'refus' and not cas['etat_restitue']) or
                (cas['statut'] == 'succes' and cas['erreur_reference'] > CRITERES['projection_erreur_absolue_max'])):
            echecs.append({'exploration': len(explorations)-1})
    for longueur in (.001, 1., 1000.):
        cas = []
        for nb in (128, 256):
            try:
                cas.append(pendule(longueur, nb))
            except (ValueError, RuntimeError) as e:
                cas.append({'longueur': longueur, 'pas': nb, 'erreur': str(e)})
        trajectoires.extend(cas)
        if any('erreur' in c for c in cas):
            echecs.append({'pendule': longueur, 'raison': 'refus'})
            continue
        ratio = cas[0]['erreur_reference']/max(cas[1]['erreur_reference'], 1e-30)
        cas[1]['rapport_erreurs'] = ratio
        c = CRITERES
        if not (cas[1]['erreur_reference'] <= c['pendule_erreur_adimensionnelle_max_256_pas']
                and cas[1]['derive_energie_relative'] <= c['pendule_derive_energie_relative_max_256_pas']
                and c['pendule_rapport_erreurs_128_sur_256_min'] <= ratio <= c['pendule_rapport_erreurs_128_sur_256_max']
                and max(x['ecart_references'] for x in cas) <= c['reference_ecart_entre_tolerances_max']):
            echecs.append({'pendule': longueur, 'raison': 'critere physique'})
    return {'schema': 1, 'version': vinkulum.__version__, 'criteres': CRITERES,
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'projections': sorties, 'explorations': explorations,
            'pendules': trajectoires, 'echecs': echecs, 'admis': not echecs}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    result = campagne()
    args.sortie.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'version': result['version'], 'projections': len(result['projections']),
                      'explorations': result['explorations'], 'pendules': result['pendules'],
                      'echecs': len(result['echecs']), 'admis': result['admis']}, ensure_ascii=False, indent=2))
    raise SystemExit(0 if result['admis'] else 1)
