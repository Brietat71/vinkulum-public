"""Archives, refus de faux témoins et vérification sans le module natif."""
import copy
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import random
import struct
import subprocess
import sys
import tempfile
import unittest

import archive_trajectoire_certifiee as a


RACINE = Path(__file__).resolve().parents[1]/'docs/bancs/trajectoire-certifiee-pendule'


class Documents(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = a.lire_document(RACINE/'pendule-64.json')

    def test_archives_et_identites(self):
        rapport = json.loads((RACINE/'qualification.json').read_text())
        self.assertEqual(rapport['version'], '0.15.0')
        self.assertEqual(set(rapport['cas']), {'64', '128', '256'})
        for nom, empreinte in rapport['modules_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(__file__).with_name(nom).read_bytes()).hexdigest(), empreinte)
        for den, attendu in rapport['cas'].items():
            chemin = RACINE/('pendule-'+den+'.json')
            self.assertEqual(hashlib.sha256(chemin.read_bytes()).hexdigest(), attendu['sha256'])
            self.assertEqual(a.verifier_document(a.lire_document(chemin)), attendu['verification_autonome'])

    def test_majorant_sous_estime_et_trace_falsifiee(self):
        doc = copy.deepcopy(self.doc)
        doc['majorants']['position'] = '0'
        with self.assertRaisesRegex(ValueError, 'sous-estimé'): a.verifier_document(doc)
        doc = copy.deepcopy(self.doc)
        doc['noeuds'][-1]['trace']['position'][0] = '1'
        with self.assertRaisesRegex(ValueError, 'sous-estimé'): a.verifier_document(doc)

    def test_temoins_et_contrat_falsifies(self):
        mutations = [
            lambda d: d.update(ordre=2),
            lambda d: d.update(bits=64),
            lambda d: d.update(pas='1/32'),
            lambda d: d.update(modele='autre'),
            lambda d: d.update(coefficient=['0', '0']),
            lambda d: d.update(initial=[['0', '0'], ['0', '0']]),
            lambda d: d['noeuds'].pop(0),
            lambda d: d['noeuds'][1].update(t='0'),
            lambda d: d['noeuds'][1].update(domaine=[['0', '0'], ['0', '0']]),
            lambda d: d['noeuds'][1].update(boite=[['0', '0'], ['0', '0']]),
            lambda d: d['noeuds'][0]['trace'].update(rotation=['0']*8),
            lambda d: d['noeuds'][0]['trace']['position'].__setitem__(0, '1/3'),
            lambda d: d.update(champ_inconnu=True),
        ]
        for mutation in mutations:
            doc = copy.deepcopy(self.doc)
            mutation(doc)
            with self.assertRaises(ValueError): a.verifier_document(doc)

    def test_autonome_sans_site_et_refus_json(self):
        p = subprocess.run([sys.executable, '-S', a.__file__, str(RACINE/'pendule-64.json')],
                           capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(p.stdout)['noeuds'], 65)
        with tempfile.TemporaryDirectory() as dossier:
            fichier = Path(dossier)/'faux.json'
            fichier.write_text('{"schema": 1, "schema": 2}')
            with self.assertRaisesRegex(ValueError, 'dupliquée'): a.lire_document(fichier)
            p = subprocess.run([sys.executable, '-S', a.__file__, str(fichier)],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 1)

    def test_rationnels_stricts(self):
        for x in (True, .5, None, '0.5', '2/4', 'NaN', '1/0', '1'*1301, '1e999999999'):
            with self.assertRaises(ValueError): a.rationnel(x)
        with self.assertRaises(ValueError): a.lire_intervalle(['1/3', '1'])

    def test_binary64_contre_decodage_binaire(self):
        r = random.Random(72461)
        for _ in range(4096):
            x = struct.unpack('>d', r.getrandbits(64).to_bytes(8, 'big'))[0]
            if math.isfinite(x):
                self.assertTrue(a.est_binary64(F(x)))
                self.assertEqual(a.rationnel(str(F(x))), F(x))
        for q in (F(1, 3), F(1, 2**1075), F(2**1024), 1+F(1, 2**53)):
            self.assertFalse(a.est_binary64(q))


if __name__ == '__main__':
    unittest.main()
