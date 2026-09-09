"""Propositions SuperLU et matrices mécaniques natives pour le prototype."""
from fractions import Fraction as F
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import subprocess
import sys

import vinkulum
import numpy as np
from scipy.sparse import bmat, csc_matrix
from scipy.sparse.linalg import splu
from vinkulum import Noyau

from prototype_lineaire_creux import certifier
import archive_lineaire_creux as archive


def lignes(a):
    a = a.tocsr(copy=True)
    a.sum_duplicates()
    return [{int(a.indices[k]): F(float(a.data[k])) for k in range(a.indptr[i], a.indptr[i+1])
             if a.data[k] != 0} for i in range(a.shape[0])]


def proposition(a, b):
    debut = time.perf_counter()
    lu = splu(a.tocsc())
    x = lu.solve(np.asarray(b, dtype=float))
    pr, pc = np.argsort(lu.perm_r), np.argsort(lu.perm_c)
    # Les deux permutations sont explicites. La vérification du défaut
    # LU-A porte sur le système ainsi permuté, pas sur une convention crue.
    entree = {'a': lignes(a[pr, :][:, pc]), 'b': [F(float(b[i])) for i in pr],
              'x': [F(float(x[i])) for i in pc], 'l': lignes(lu.L), 'u': lignes(lu.U)}
    return entree, pr.tolist(), pc.tolist(), time.perf_counter()-debut


def chaine(nb):
    n = Noyau([0., 0., -9.81])
    for i in range(nb):
        n.corps(str(i), 1., np.eye(3).ravel().tolist(), [0., 0., -i-.5])
        n.liaison(str(i), None if i == 0 else i-1, i,
                  pa=[0., 0., 0. if i == 0 else -.5], bloque_r=[0, 2])
    etat = n.etat_precis()
    export = n.k_c_m_g_creux()
    if n.etat_precis() != etat:
        raise AssertionError('export non immuable')
    def csc(e):
        nr, nc, ptr, idx, val = e
        return csc_matrix((val, idx, ptr), shape=(nr, nc))
    m, g = csc(export[2]), csc(export[3])
    a = bmat([[m, g.T], [g, None]], format='csc')
    b = np.zeros(a.shape[0])
    b[2:6*nb:6] = -9.81
    return a, b


def experience(nb):
    a, b = chaine(nb)
    e, pr, pc, preparation = proposition(a, b)
    debut = time.perf_counter()
    try:
        preuve = certifier(**e)
        resultat = {'decision': 'certifie', 'contraction': str(preuve['contraction']),
                    'erreur_max': str(max(preuve['bornes_composantes'])),
                    'produits_facteurs': preuve['produits_facteurs'],
                    'termes_defaut': preuve['termes_defaut']}
    except ValueError as exc:
        resultat = {'decision': 'refus_de_preuve', 'raison': str(exc)}
    resultat.update(corps=nb, dimension=a.shape[0], nnz=a.nnz,
                    preparation_s=preparation, verification_s=time.perf_counter()-debut)
    return resultat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sortie', type=Path)
    parser.add_argument('--roue', type=Path)
    args = parser.parse_args()
    if args.sortie is None:
        print(json.dumps([experience(n) for n in (1, 4, 16, 64, 256)], indent=2))
        return
    if args.roue is None:
        parser.error('--roue requis pour archiver une qualification')
    args.sortie.mkdir(parents=True, exist_ok=False)
    rapport = {'version': vinkulum.__version__, 'python': sys.version,
               'paquet': str(Path(vinkulum.__file__).resolve()),
               'roue_sha256': hashlib.sha256(args.roue.read_bytes()).hexdigest(), 'cas': {},
               'modules_sha256': {nom: hashlib.sha256(Path(__file__).with_name(nom).read_bytes()).hexdigest()
                  for nom in ('prototype_lineaire_creux.py', 'archive_lineaire_creux.py', 'qualifie_lineaire_creux.py')}}
    for nb in (1, 4, 16, 64, 256):
        debut = time.perf_counter()
        a, b = chaine(nb)
        modele_export = time.perf_counter()-debut
        e, pr, pc, preparation = proposition(a, b)
        debut = time.perf_counter()
        doc = archive.document(e, pr, pc)
        # Contrôle de la reconstruction des indices depuis le système permuté.
        if archive.lire_matrice(doc['a'], a.shape[0]) != lignes(a) or list(map(F, doc['b'])) != list(map(F, b)):
            raise AssertionError('système original non conservé par les permutations')
        construction = time.perf_counter()-debut
        chemin = args.sortie/('chaine-'+str(nb)+'.json.gz')
        contenu = json.dumps(doc, separators=(',', ':'), allow_nan=False).encode()
        chemin.write_bytes(gzip.compress(contenu, mtime=0))
        debut = time.perf_counter()
        result = subprocess.run([sys.executable, '-S', archive.__file__, str(chemin)],
                                capture_output=True, text=True, check=True)
        temps_verification = time.perf_counter()-debut
        rapport['cas'][str(nb)] = {'sha256': hashlib.sha256(chemin.read_bytes()).hexdigest(),
            'verification_autonome': json.loads(result.stdout), 'modele_export_s': modele_export,
            'proposition_s': preparation, 'construction_document_s': construction,
            'relecture_processus_s': temps_verification, 'octets_json': len(contenu),
            'octets_gzip': chemin.stat().st_size}
    (args.sortie/'qualification.json').write_text(json.dumps(rapport, indent=2)+'\n')
    print(json.dumps(rapport, indent=2))


if __name__ == '__main__': main()
