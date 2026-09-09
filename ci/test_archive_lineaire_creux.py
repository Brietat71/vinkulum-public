"""Archives natives, permutations et falsifications du certificat creux."""
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import unittest

import archive_lineaire_creux as a

RACINE = Path(__file__).resolve().parents[1]/'docs/bancs/lineaire-creux-prototype'


class Documents(unittest.TestCase):
    def test_cinq_archives_et_empreintes(self):
        rapport = json.loads((RACINE/'qualification.json').read_text())
        self.assertEqual(set(rapport['cas']), {'1', '4', '16', '64', '256'})
        for nom, empreinte in rapport['modules_sha256'].items():
            self.assertEqual(hashlib.sha256(Path(__file__).with_name(nom).read_bytes()).hexdigest(), empreinte)
        for nb, cas in rapport['cas'].items():
            p = RACINE/('chaine-'+nb+'.json.gz')
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), cas['sha256'])
            r = a.verifier_document(a.lire_document(p))
            self.assertEqual(r['dimension'], 11*int(nb))
            r['erreur_max'] = str(max(map(F, r.pop('bornes_composantes'))))
            self.assertEqual(r, cas['verification_autonome'])

    def test_falsifications(self):
        original = a.lire_document(RACINE/'chaine-4.json.gz')
        mutations = [lambda d: d.update(contraction='0'),
                     lambda d: d.update(bornes_composantes=['0']*44),
                     lambda d: d['pr'].__setitem__(0, d['pr'][1]),
                     lambda d: d['x'].__setitem__(0, '1000'),
                     lambda d: d['l'][0].__setitem__(0, [0, '0']),
                     lambda d: d.update(champ_inconnu=0)]
        for mutation in mutations:
            d = copy.deepcopy(original)
            mutation(d)
            with self.assertRaises(ValueError): a.verifier_document(d)

    def test_enveloppe_archivee(self):
        d = {'schema': a.SCHEMA, 'dimension': 1, 'a': [[[0, '2']]], 'b': ['1'], 'x': ['1/2'],
             'l': [[[0, '1']]], 'u': [[[0, '2']]], 'pr': [0], 'pc': [0], 'poids': ['1'],
             'delta_a': [[[0, '1/4']]], 'delta_b': ['1/10'], 'contraction': '1/8',
             'bornes_composantes': ['13/100']}
        self.assertLessEqual(F(a.verifier_document(d)['bornes_composantes'][0]), F(13, 100))
        d['bornes_composantes'] = ['1/8']
        with self.assertRaisesRegex(ValueError, 'sous-estimée'): a.verifier_document(d)

    def test_sans_module_natif(self):
        r = subprocess.run([sys.executable, '-S', a.__file__, str(RACINE/'chaine-256.json.gz')],
                           capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(r.stdout)['dimension'], 2816)


if __name__ == '__main__': unittest.main()
