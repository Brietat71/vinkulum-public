"""Coût de l'assemblage déjà convergé, avec ou sans garde entier exact.

Exécuter ce même fichier avec les interpréteurs des roues à comparer.
Médianes de 9 lots ; le modèle et son premier assemblage sont hors mesure.
Ce banc isole le coût ajouté au contrôle, sans extrapoler à la simulation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import time

import vinkulum
from vinkulum import Noyau, _vinkulum


def mesure(nombre, repetitions, mobile):
    n = Noyau([0., 0., 0.])
    j = [1., 0., 0., 0., 1., 0., 0., 0., 1.]
    for i in range(nombre):
        n.corps(str(i), 1., j, [0., 0., 0.],
                v=[0., -.6 if mobile else 2., 3.], w=[0., .1, .2])
        # Le cas mobile force des produits entiers non nuls qui se compensent,
        # pour ne pas mesurer uniquement le raccourci des facteurs nuls.
        n.liaison(str(i), None, i, pa=[3. if mobile else 0., 0., 0.],
                  bloque_t=[1 if mobile else 0], bloque_r=[])
    n.assemble()
    lots = []
    for _ in range(9):
        debut = time.perf_counter_ns()
        for _ in range(repetitions):
            n.assemble()
        lots.append((time.perf_counter_ns()-debut)/repetitions)
    return {'famille': 'levier_mobile' if mobile else 'translation_immobile',
            'corps': nombre, 'contraintes': nombre, 'colonnes': 6*nombre,
            'repetitions_par_lot': repetitions, 'lots_ns_par_appel': lots,
            'mediane_ns_par_appel': statistics.median(lots)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    r = {'version': vinkulum.__version__,
         'extension_sha256': hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
         'portee': 'Assemblage répété d’un état déjà convergé ; pas un temps de simulation.',
         'cas': [mesure(n, reps, mobile) for mobile in (False, True)
                 for n, reps in ((1, 2000), (10, 500), (50, 80), (100, 30))]}
    args.sortie.write_text(json.dumps(r, ensure_ascii=False, indent=2)+'\n')
    for c in r['cas']:
        print(c['famille'], c['corps'], c['mediane_ns_par_appel']/1000, 'µs/appel')


if __name__ == '__main__':
    main()
