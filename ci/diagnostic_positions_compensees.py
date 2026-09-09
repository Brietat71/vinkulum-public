"""Compare l'équilibre strict et archive les positions avec leurs petits termes.

Même programme utilisable avec l'extension antérieure : l'absence d'état
précis est enregistrée. Les échecs restent dans le rapport.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
from vinkulum import Noyau, _vinkulum
from confronte_mbdyn import princeton
from diagnostic_statut_statique import appelle
from diagnostic_precision_statique import rotations_moyennes, residu_translation

ROOT = Path(__file__).resolve().parents[1]


def etat(n):
    return n.etat_precis() if hasattr(n, 'etat_precis') else n.etat()


def barre(element):
    n = Noyau([0., 0., 0.])
    for i in range(2):
        n.corps(str(i), 1., np.eye(3).ravel().tolist(), [float(i), 0., 0.])
    n.liaison('enc', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    if element == 'super':
        k = np.eye(12)
        k[6, 6] = 2.**30
        n.superelement('barre', [0, 1], k.ravel().tolist())
    else:
        n.poutre('barre', 0, 1, 2.**30, 1., 1., 1., formulation=element)
    n.effort(1, [2.**-30, 0., 0.], [0., 0., 0.])
    before = etat(n)
    report = appelle(n, True, iters=40, tol=1e-14)
    report.update(element=element, etat=etat(n),
                  residu_courant=n.residu_statique(),
                  restauration_precise_verifiee=(etat(n) == before if report['erreur'] else None),
                  allongement_reference_m=2.**-60)
    return report


def parcours(ne, formulation):
    n, ids = princeton(ne, formulation)
    initial = n.etat()
    loads = [.5 * (1 - math.cos(math.pi * k / 50)) for k in range(1, 51)]
    out = dict(intervalles=ne, formulation=formulation, appels=[])
    previous = 0.
    for index, load in enumerate(loads, 1):
        force = 8.896 * (load - previous) / math.sqrt(2)
        n.effort(ids[-1], [0., force, force], [0., 0., 0.])
        before = etat(n)
        report = appelle(n, True)
        report.update(palier_externe=index, charge=load,
                      restauration_precise_verifiee=(etat(n) == before if report['erreur'] else None))
        out['appels'].append(report)
        if index == 1:
            out['etat_premier_palier'] = etat(n)
        if report['erreur'] is not None:
            break
        previous = load
    out.update(charge_finale_acceptee=previous, etat_final=etat(n),
               residu_final=n.residu_statique(),
               tolerance_tous_paliers=(previous == 1. and all(
                   x['rapport']['statut'] == 'tolerance' for x in out['appels'])))
    # Contrôle des forces au premier palier aux poses réellement conservées,
    # sans reconstruire l'équilibre. Calcul limité aux translations de milieu.
    if (formulation == 'milieu' and len(out['etat_premier_palier']) == 7
            and out['appels'][0]['erreur'] is None
            and np.finfo(np.longdouble).nmant > 52):
        state = out['etat_premier_palier']
        positions = np.array(state[1], dtype=np.longdouble) + np.array(state[6], dtype=np.longdouble)
        means = rotations_moyennes(state[2], np.longdouble)
        force = 8.896 * loads[0] / math.sqrt(2)
        r = residu_translation(positions, means, np.diff(np.array(initial[1])[:, 0]), [0., force, force])
        out['residu_translations_premier_palier_calcul_etendu_N'] = float(np.max(np.abs(r)))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sortie', required=True, type=Path)
    parser.add_argument('--barres-seules', action='store_true')
    args = parser.parse_args()
    files = ['src/lib.rs', 'src/tangent.rs', 'src/numerique.rs', 'src/robustesse.rs',
             'src/superelement.rs', 'ci/diagnostic_positions_compensees.py',
             'ci/confronte_mbdyn.py', 'ci/diagnostic_statut_statique.py',
             'ci/diagnostic_precision_statique.py']
    out = dict(python=platform.python_version(), numpy=np.__version__,
               extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
               # Ces empreintes décrivent le worktree du lanceur. Le binaire
               # antérieur est identifié par son empreinte séparée ci-dessus.
               worktree_lanceur_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
               etat_precis_disponible=hasattr(Noyau, 'etat_precis'), barres=[], princeton=[])
    def save():
        args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n')
    for element in ('milieu', 'integree', 'super'):
        result = barre(element)
        out['barres'].append(result)
        save()
        print(element, result['rapport']['statut'], flush=True)
    if not args.barres_seules:
        for formulation in ('milieu', 'integree'):
            for ne in (10, 20, 40, 60):
                result = parcours(ne, formulation)
                out['princeton'].append(result)
                save()
                print(formulation, ne, result['tolerance_tous_paliers'],
                      result['charge_finale_acceptee'], flush=True)


if __name__ == '__main__':
    main()
