"""Relecture des certificats produits par la roue isolée 0.15.0."""
import hashlib
import json
from pathlib import Path
import unittest

from vinkulum.certification import verifier_certificat


RACINE = Path(__file__).resolve().parents[1]/'docs/bancs/assemblage-certifie-0.15.0'


class ArchivePublique(unittest.TestCase):
    def test_documents_de_la_roue_isolee(self):
        rapport = json.loads((RACINE/'qualification.json').read_text())
        self.assertEqual(rapport['version'], '0.15.0')
        self.assertEqual(set(rapport['cas']), {'micro', 'metre', 'grand', 'double', 'repere-tourne'})
        for nom, preuve in rapport['cas'].items():
            with self.subTest(cas=nom):
                p = RACINE/(nom+'.json')
                self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), preuve['sha256'])
                doc = json.loads(p.read_text())
                resultat = verifier_certificat(doc)
                self.assertTrue(resultat['existence_unique_locale'])
                self.assertEqual(resultat['dimension'], preuve['dimension'])
                self.assertEqual(len(resultat['bornes_pose']), doc['modele']['corps'])


if __name__ == '__main__':
    unittest.main()
