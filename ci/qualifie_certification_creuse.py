"""Qualification des API installées : aucun import des preuves source de ci."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.sparse import bmat, csc_matrix, coo_matrix, eye
import vinkulum
from vinkulum import Noyau, _verification_lineaire
from vinkulum.certification import certifier_systeme_creux, certifier_quotient_structurel, verifier_certificat


def chaine(nb):
    n = Noyau([0., 0., -9.81])
    for i in range(nb):
        n.corps(str(i), 1., np.eye(3).ravel().tolist(), [0., 0., -i-.5])
        n.liaison(str(i), None if i == 0 else i-1, i, pa=[0., 0., 0. if i == 0 else -.5], bloque_r=[0, 2])
    avant = n.etat_precis()
    exports = n.k_c_m_g_creux()
    if n.etat_precis() != avant:
        raise AssertionError('export non immuable')
    def csc(e):
        nr, nc, ptr, idx, data = e
        return csc_matrix((data, idx, ptr), shape=(nr, nc))
    m, g = csc(exports[2]), csc(exports[3])
    a = bmat([[m, g.T], [g, None]], format='csc')
    b = np.zeros(a.shape[0]); b[2:6*nb:6] = -9.81
    avant = (a.data.copy(), a.indices.copy(), a.indptr.copy(), b.copy())
    doc = certifier_systeme_creux(a, b)
    for x, y in zip(avant, (a.data, a.indices, a.indptr, b), strict=True):
        np.testing.assert_array_equal(x, y)
    return doc


def cas():
    for nb in (1, 4, 16, 64, 256):
        yield 'chaine-'+str(nb), chaine(nb)
    dup = coo_matrix(([1e16, 1., -1e16], ([0, 0, 0], [0, 0, 0])), shape=(1, 1))
    yield 'contributions', certifier_systeme_creux(dup, [1.])
    yield 'structure-rang-1', certifier_quotient_structurel(
        eye(2, format='csc'), csc_matrix([[1., 0.]]), csc_matrix([[1.], [2.], [-3.]]),
        [0], [0., 1.], [0.], delta_m=eye(2, format='csc')*.001,
        delta_c=csc_matrix([[.01, .01]]), delta_t=csc_matrix([[0.], [10.], [10.]]),
        delta_force=[.002, .002], delta_d=[.001])
    yield 'structure-rang-2', certifier_quotient_structurel(
        eye(3, format='csc'), csc_matrix([[1., 0., .2], [0., 1., -.3]]),
        csc_matrix([[0., 1.], [1., -1.], [1., 0.]]), [2, 0], [1., 2., 3.], [0., 0.],
        delta_c=csc_matrix(np.full((2, 3), .001)),
        delta_t=csc_matrix([[0., 0.], [.2, .3], [0., 0.]]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sortie', type=Path, required=True)
    parser.add_argument('--roue', type=Path, required=True)
    args = parser.parse_args()
    args.sortie.mkdir(parents=True, exist_ok=False)
    package = Path(vinkulum.__file__).resolve().parent
    rapport = {'version': vinkulum.__version__, 'python': sys.version, 'paquet': str(package),
               'roue_sha256': hashlib.sha256(args.roue.read_bytes()).hexdigest(), 'cas': {},
               'modules_sha256': {nom: hashlib.sha256((package/nom).read_bytes()).hexdigest()
                   for nom in ('_api_certification_creuse.py', '_certification_creuse.py',
                               '_verification_creuse.py', '_verification_structurelle.py',
                               '_verification_lineaire.py', 'certification.py')}}
    for nom, doc in cas():
        preuve = verifier_certificat(doc)
        contenu = json.dumps(doc, separators=(',', ':'), allow_nan=False).encode()
        p = args.sortie/(nom+'.json.gz')
        p.write_bytes(gzip.compress(contenu, mtime=0))
        with tempfile.TemporaryDirectory() as tmp:
            fichier = Path(tmp)/'preuve.json'; fichier.write_bytes(contenu)
            r = subprocess.run([sys.executable, '-S', _verification_lineaire.__file__, str(fichier)],
                               capture_output=True, text=True, check=True)
        rapport['cas'][nom] = {'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
            'schema': doc['schema'], 'dimension': preuve['dimension'],
            'borne_erreur_inf': str(preuve['borne_erreur_inf']), 'verification_autonome': r.stdout.strip()}
    (args.sortie/'qualification.json').write_text(json.dumps(rapport, indent=2)+'\n')
    print(json.dumps(rapport, indent=2))


if __name__ == '__main__': main()
