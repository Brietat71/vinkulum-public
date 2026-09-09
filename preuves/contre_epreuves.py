"""Confrontation Lean compilé / Rust installé / rationnels CPython indépendants."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'ci'))
from certifie_noyau import echantillons, UNITE, MARGE
from vinkulum import _vinkulum


def mot(x):
    return str(int.from_bytes(struct.pack('>d', x), 'big'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', default=str(Path(__file__).parent / '.lake/build/bin/garde'))
    args = parser.parse_args()
    cas = list(echantillons())
    lignes = [' '.join(map(mot, [tol, dt] + [x for pair in zip(g, u, strict=True) for x in pair]))
              for _, g, u, dt, tol in cas]
    invalides = ['0 0 1', '0 0 18446744073709551616 0', '0 0 9218868437227405312 0',
                 '0 0 0 9221120237041090560', '0 18442240474082181120',
                 '9218868437227405312 0', '13830554455654793216 0']
    resultat = subprocess.run([args.executable], input='\n'.join(lignes + invalides) + '\n',
                              text=True, capture_output=True, check=True, timeout=180)
    sorties = resultat.stdout.splitlines()
    if len(sorties) != len(lignes) + len(invalides):
        raise AssertionError('nombre de réponses Lean incorrect')
    for (nom, g, u, dt, tol), sortie in zip(cas, sorties):
        ok, r, e, t = sortie.split()
        if ok not in ('true', 'false'):
            raise AssertionError(f'{nom}: décision mal formée')
        lean = (ok == 'true', int(r), int(e), int(t))
        okr, rr, er, tr = _vinkulum._certificat_ligne_vitesse(g, u, dt, tol)
        rust = (okr, int(rr, 16), int(er, 16), int(tr, 16))
        produits = [F(a) * F(b) for a, b in zip(g, u, strict=True)]
        residu = F(dt) + sum(produits, F(0))
        echelle = abs(F(dt)) + sum(map(abs, produits), F(0))
        rationnel = (abs(residu) <= F(tol) + MARGE * echelle,
                     residu / UNITE, echelle / UNITE, F(tol) / UNITE)
        if lean != rust or lean != rationnel:
            raise AssertionError(f'{nom}: divergence Lean/Rust/Fraction')
    if any(not s.startswith('ERREUR ') for s in sorties[len(cas):]):
        raise AssertionError('entrée hors contrat acceptée par Lean')
    print(json.dumps({'cas_valides': len(cas), 'invalides_lean_refuses': len(invalides),
                      'divergences': 0, 'entrees_sha256': hashlib.sha256(
                          ('\n'.join(lignes + invalides) + '\n').encode()).hexdigest()}, indent=2))


if __name__ == '__main__':
    main()
