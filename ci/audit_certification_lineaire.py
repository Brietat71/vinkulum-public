"""Campagne de certificats linéaires et d'initialisation à état figé."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

import numpy as np
from vinkulum import Noyau
from vinkulum.certification import certifier_initialisation, certifier_systeme, verifier_certificat
from vinkulum.certification import CertificationImpossible


CRITERES = {'borne_initialisation_max': 1e-9, 'borne_systemes_denses_max': 1e-12}


def pendule(longueur, angle):
    n = Noyau([0., -9.81, 0.])
    c, s = math.cos(angle), math.sin(angle)
    rot = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    pos = rot @ [0., -longueur, 0.]
    vitesse = np.cross([0., 0., .3], pos)
    n.corps('pendule', 1., (.2*longueur**2*np.eye(3)).ravel().tolist(), pos.tolist(),
            rot=rot.ravel().tolist(), v=vitesse.tolist(), w=[0., 0., .3])
    n.liaison('pivot', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1])
    n.assemble()
    return n


def campagne():
    cas = []
    def conserve(nom, produit):
        debut = time.perf_counter()
        c = produit()
        generation = time.perf_counter()-debut
        debut = time.perf_counter()
        q = verifier_certificat(c)
        cas.append({'nom': nom, 'certificat': c, 'dimension': q['dimension'],
                    'borne_erreur_inf': float(q['borne_erreur_inf']),
                    'eta_exact': str(q['eta_exact']),
                    'secondes_generation': generation,
                    'secondes_verification_independante': time.perf_counter()-debut})
    for longueur in (.001, 1., 1000.):
        for angle in (0., .4):
            n = pendule(longueur, angle)
            avant = n.etat_precis(), n.stats(), n.schema()
            conserve(f'pendule_L{longueur}_angle{angle}',
                     lambda: certifier_initialisation(n, erreur_max=CRITERES['borne_initialisation_max']))
            if (n.etat_precis(), n.stats(), n.schema()) != avant:
                raise AssertionError('la certification a modifié le modèle')
    rng = np.random.default_rng(13002026)
    for n in (2, 6, 12, 24):
        a = rng.integers(-3, 4, (n, n)).astype(float)+4*n*np.eye(n)
        b = rng.integers(-10, 11, n).astype(float)
        conserve(f'dense_{n}', lambda: certifier_systeme(
            a, b, erreur_max=CRITERES['borne_systemes_denses_max']))
    conserve('residu_petit_erreur_un', lambda: certifier_systeme(
        np.diag([1., 2.**-60]), [1., 2.**-60], [1., 0.]))
    faux = cas[-1]['certificat']
    if F(float.fromhex(faux['borne_erreur_inf'])) != 1:
        raise AssertionError('borne unitaire du contre-exemple absente')
    try:
        certifier_systeme(np.diag([1., 2.**-60]), [1., 2.**-60], [1., 0.], erreur_max=1e-12)
    except CertificationImpossible:
        refus_precision = True
    else:
        raise AssertionError('erreur unitaire admise à 1e-12')
    petit = math.ulp(0.)
    conserve('contraction_1_moins_sous_normal', lambda: certifier_systeme(
        [[petit]], [petit], [1.], inverse_approche=[[1.]], erreur_max=0.))
    return {'criteres': CRITERES, 'cas': cas, 'certificats_verifies': len(cas),
            'refus_precision_insuffisante': refus_precision,
            'portee': 'Systèmes linéaires à coefficients stockés ; géométrie et trajectoires non certifiées.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    resultat = campagne()
    args.sortie.write_text(json.dumps(resultat, ensure_ascii=False, indent=2)+'\n')
    print(resultat['certificats_verifies'], 'certificats vérifiés indépendamment.')
