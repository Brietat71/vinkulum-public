"""Coût des corrections de jacobien, à pas et sorties identiques."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import time

from vinkulum import _vinkulum
from vinkulum.test_noyau import _pale_inflow_spatial, _roulement_charge


def modele(nom):
    if nom == 'roulement':
        return _roulement_charge()
    n = _pale_inflow_spatial('harmoniques' if nom == 'uniforme' else nom)
    if nom == 'uniforme':
        n.pitt_peters(0, [0., 0., 0.], [1., 0., 0.], 10., actif=False)
    return n


def run(nom):
    n = modele(nom)
    t = time.perf_counter()
    n.simule(.25, .001, tous=1000000000)
    elapsed = time.perf_counter()-t
    return dict(secondes=elapsed, stats=n.stats(), etat=n.etat_precis())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    out = dict(cpu=cpu, pas_s=.001, duree_s=.25, cas={},
               extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
               programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    for nom in ('uniforme', 'harmoniques', 'profil', 'carte', 'roulement'):
        run(nom)
        samples = [run(nom) for _ in range(5)]
        out['cas'][nom] = dict(echantillons=samples, mediane_s=statistics.median(s['secondes'] for s in samples))
        print(nom, out['cas'][nom]['mediane_s'], samples[0]['stats'], flush=True)
        args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n')


if __name__ == '__main__':
    main()
