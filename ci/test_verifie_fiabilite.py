"""Altérer une décision ou supprimer un résultat doit invalider le dossier."""
import gzip
import json
from pathlib import Path
import unittest

from verifie_fiabilite import controle_physique, violations_lisses

ARCHIVE = Path(__file__).resolve().parents[1]/'docs/bancs/fiabilite-0.14.1-preuves.json.gz'


class ContreEpreuves(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.files = json.loads(gzip.decompress(ARCHIVE.read_bytes()))['fichiers']

    def physique(self):
        return json.loads(self.files['physique-candidat.json'])

    def test_un_ecart_ne_peut_pas_etre_reclasse_vert(self):
        p = self.physique()
        ligne = next(r for f in p['familles'] for r in f['lignes'] if r['statut'] == 'ecart')
        ligne['statut'] = 'ok'
        p['comptes']['ok'] += 1
        p['comptes']['ecart'] -= 1
        with self.assertRaisesRegex(ValueError, 'verdict physique'):
            controle_physique(p)

    def test_une_reference_absente_ne_peut_pas_rester_complete(self):
        p = self.physique()
        p['familles'][0]['lignes'].pop()
        p['comptes']['ok'] -= 1
        with self.assertRaises(ValueError):
            controle_physique(p)

    def test_un_echec_mecanique_ne_peut_pas_etre_declare_admis(self):
        p = json.loads(self.files['lisse-reference.json'])
        p['admis'], p['echecs'] = True, []
        with self.assertRaisesRegex(ValueError, 'décision du corpus lisse'):
            violations_lisses(p)

    def test_corpus_candidat_et_reference(self):
        self.assertEqual(violations_lisses(json.loads(self.files['lisse-reference.json'])), 32)
        self.assertEqual(violations_lisses(json.loads(self.files['lisse-candidat.json'])), 0)
        controle_physique(self.physique())


if __name__ == '__main__':
    unittest.main()
