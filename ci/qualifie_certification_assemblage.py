"""Qualification de l'API installée ; aucun import des modules source de ci."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import vinkulum
from vinkulum import Noyau
from vinkulum.certification import certifier_assemblage, verifier_certificat
import vinkulum._verification_lineaire as verification


J = [1., 0., 0., 0., 1., 0., 0., 0., 1.]


def cas():
    for nom, longueur in [('micro', 1e-6), ('metre', 1.), ('grand', 1e6)]:
        n = Noyau([0., 0., -9.81])
        n.corps('pendule', 1., J, [0., 0., -longueur])
        n.liaison('pivot', None, 0, bloque_r=[0, 2])
        yield nom, n, [(5, 0)], [F(longueur)*F(1, 10**8)]*3+[F(1, 10**8)]*4
    import numpy as np
    from scipy.spatial.transform import Rotation
    repere = Rotation.from_rotvec([.2, -.3, .4]).as_matrix()
    n = Noyau([0., 0., -9.81])
    n.corps('tourne', 1., J, (repere @ np.array([0., 0., -1.])).tolist(), rot=repere.ravel().tolist())
    n.liaison('pivot', None, 0, ra=repere.ravel().tolist(), bloque_r=[0, 2])
    qy = float(Rotation.from_matrix(repere).as_quat()[1])
    yield 'repere-tourne', n, [(5, qy)], [F(1,10**6)]*7
    n = Noyau([0., 0., -9.81])
    n.corps('a', 1., J, [0., 0., -.5])
    n.corps('b', 1., J, [0., 0., -1.5])
    n.liaison('sol', None, 0, bloque_r=[0, 2])
    n.liaison('coude', 0, 1, pa=[0., 0., -.5], bloque_r=[0, 2])
    yield 'double', n, [(5, 0), (12, 0)], [F(1, 10**6)]*14


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sortie', type=Path, required=True)
    parser.add_argument('--roue', type=Path, required=True)
    args = parser.parse_args()
    args.sortie.mkdir(parents=True, exist_ok=False)
    package = Path(vinkulum.__file__).resolve().parent
    rapport = {'version': vinkulum.__version__, 'python': sys.version, 'paquet': str(package),
               'roue_sha256': hashlib.sha256(args.roue.read_bytes()).hexdigest(),
               'cas': {}, 'modules_sha256': {}}
    for nom, n, jauges, rayons in cas():
        avant = n.etat_precis()
        doc = certifier_assemblage(n, jauges=jauges, rayons=rayons)
        resultat = verifier_certificat(doc)
        if n.etat_precis() != avant:
            raise AssertionError('certification non immuable')
        p = args.sortie/(nom+'.json')
        p.write_text(json.dumps(doc, indent=2)+'\n')
        r = subprocess.run([sys.executable, '-S', verification.__file__, str(p)],
                           capture_output=True, text=True, check=True)
        if 'Certificat géométrique vérifié' not in r.stdout:
            raise AssertionError('vérification autonome absente')
        rapport['cas'][nom] = {'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
            'dimension': resultat['dimension'], 'etat_inchange': True,
            'verification_autonome': r.stdout.strip()}
    for nom in ('_verification_lineaire.py', '_verification_assemblage.py', '_oracle_assemblage.py',
                '_geometrie_certifiee.py', 'certification.py'):
        rapport['modules_sha256'][nom] = hashlib.sha256((package/nom).read_bytes()).hexdigest()
    (args.sortie/'qualification.json').write_text(json.dumps(rapport, indent=2)+'\n')
    print(json.dumps(rapport, indent=2))


if __name__ == '__main__':
    main()
