"""Relecture publique des documents livrés et des prototypes historiques."""
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import unittest

from vinkulum.certification import verifier_certificat

RACINE = Path(__file__).resolve().parents[1]/'docs/bancs'


class Archives(unittest.TestCase):
    def test_api_de_la_roue(self):
        dossier = RACINE/'certification-creuse-0.16.0'
        rapport = json.loads((dossier/'qualification.json').read_text())
        self.assertEqual(rapport['version'], '0.16.0')
        self.assertEqual(set(rapport['cas']), {'chaine-1', 'chaine-4', 'chaine-16', 'chaine-64',
                                              'chaine-256', 'contributions', 'structure-rang-1', 'structure-rang-2'})
        for nom, attendu in rapport['cas'].items():
            p = dossier/(nom+'.json.gz')
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), attendu['sha256'])
            r = verifier_certificat(json.loads(gzip.decompress(p.read_bytes())))
            self.assertEqual(r['dimension'], attendu['dimension'])
            self.assertEqual(r['borne_erreur_inf'], F(attendu['borne_erreur_inf']))

    def test_prototypes_historiques(self):
        dossier = RACINE/'lineaire-creux-prototype'
        rapport = json.loads((dossier/'qualification.json').read_text())
        for nb, attendu in rapport['cas'].items():
            p = dossier/('chaine-'+nb+'.json.gz')
            r = verifier_certificat(json.loads(gzip.decompress(p.read_bytes())))
            self.assertEqual(r['borne_erreur_inf'], F(attendu['verification_autonome']['erreur_max']))


if __name__ == '__main__': unittest.main()
