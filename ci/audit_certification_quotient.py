"""Campagne de quotients exacts ; critères déclarés avant l'exécution."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from vinkulum import Noyau
from vinkulum.certification import CertificationImpossible, certifier_initialisation, verifier_certificat
from audit_certification_lineaire import pendule


CRITERES = {'borne_representant_max': 1e-9, 'rang_pendule': 5,
            'contraintes_pendule': 10, 'faux_rang_doit_etre_refuse': True}


def campagne():
    cas = []

    def conserve(nom, n, rang, contraintes):
        avant = n.etat_precis(), n.stats(), n.schema()
        debut = time.perf_counter()
        d = certifier_initialisation(n, redondances=True,
                                    erreur_max=CRITERES['borne_representant_max'])
        generation = time.perf_counter()-debut
        debut = time.perf_counter()
        q = verifier_certificat(d)
        verification = time.perf_counter()-debut
        assert (q['rang_contraintes'], q['nombre_contraintes']) == (rang, contraintes)
        assert (n.etat_precis(), n.stats(), n.schema()) == avant
        cas.append({'nom': nom, 'certificat': d, 'dimension': q['dimension'],
                    'rang': rang, 'contraintes': contraintes,
                    'borne_erreur_inf': float(q['borne_erreur_inf']),
                    'bornes_reaction_generalisee': list(map(float, q['bornes_reaction_generalisee'])),
                    'secondes_generation': generation, 'secondes_verification_independante': verification})

    for longueur in (.001, 1., 1000.):
        for angle in (0., .4):
            n = pendule(longueur, angle)
            n.liaison('copie', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1])
            n.assemble()
            conserve(f'pendule_double_L{longueur}_angle{angle}', n,
                     CRITERES['rang_pendule'], CRITERES['contraintes_pendule'])

    n = Noyau([0., 0., -9.81])
    n.corps('libre', 1., np.eye(3).ravel().tolist(), [0., 0., 0.])
    conserve('sans_contrainte', n, 0, 0)

    n = Noyau([0., 0., -9.81])
    n.corps('non_holonome', 1., np.eye(3).ravel().tolist(), [0., 0., 0.], v=[1., 0., 0.])
    for nom in ('commande', 'copie'):
        n.liaison(nom, None, 0, bloque_t=[0], bloque_r=[], nh=True,
                  cible_t=([1., 0., 0.], ('lineaire', [0., 1.])))
    n.assemble()
    conserve('non_holonome_double', n, 1, 2)

    n = Noyau([0., 0., 0.])
    n.corps('mobile', 1., np.eye(3).ravel().tolist(), [0., 0., 0.])
    n.liaison('a', None, 0, bloque_t=[0], bloque_r=[])
    n.liaison('b', None, 0, pa=[0., 2.**-80, 0.], bloque_t=[0], bloque_r=[])
    n.assemble()
    avant = n.etat_precis(), n.stats(), n.schema()
    try:
        certifier_initialisation(n, redondances=True)
    except CertificationImpossible as e:
        assert 'indépendante' in str(e)
        refus = {'raison': str(e), 'etat_restitue': (n.etat_precis(), n.stats(), n.schema()) == avant,
                 'bras_de_levier_hex': (2.**-80).hex()}
        assert refus['etat_restitue']
    else:
        raise AssertionError('une redondance seulement numérique a été certifiée')
    return {'criteres': CRITERES, 'certificats_verifies': len(cas), 'cas': cas,
            'faux_rang_refuse': refus,
            'portee': 'Dépendances des coefficients stockés, compatibilité, accélération et réaction uniques ; jauge native des multiplicateurs.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    bilan = campagne()
    args.sortie.write_text(json.dumps(bilan, ensure_ascii=False, indent=2)+'\n')
    print(bilan['certificats_verifies'], 'quotients vérifiés indépendamment ; faux rang refusé.')
