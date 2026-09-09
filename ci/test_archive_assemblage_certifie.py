import copy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from archive_assemblage_certifie import fabriquer, verifier_document
from prototype_assemblage_certifie import proposer_inverse
from test_prototype_assemblage_certifie import pendule, CENTRE


class ArchiveGeometrique(unittest.TestCase):
    def setUp(self):
        m = pendule()
        self.document = fabriquer(m, CENTRE, [F(1, 10**6)]*7, proposer_inverse(m, CENTRE))

    def test_verification_sans_site_packages(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'preuve.json'
            p.write_text(json.dumps(self.document))
            r = subprocess.run([sys.executable, '-S', str(Path(__file__).with_name(
                'archive_assemblage_certifie.py')), str(p)], capture_output=True, text=True, check=True)
            self.assertTrue(json.loads(r.stdout)['inclusion_stricte'])

    def test_documents_mecaniques_conserves(self):
        racine = Path(__file__).resolve().parents[1]/'docs/bancs/assemblage-certifie-prototype'
        noms = {'pendule-micro.json', 'pendule.json', 'pendule-grand.json', 'double-pendule.json'}
        self.assertEqual({p.name for p in racine.glob('*.json')}, noms)
        for nom in sorted(noms):
            with self.subTest(document=nom):
                self.assertTrue(verifier_document(json.loads((racine/nom).read_text())))

    def test_alterations_refusees(self):
        alterations = [
            lambda d: d['centre'].__setitem__(0, '1'),
            lambda d: d['rayons'].__setitem__(0, '0'),
            lambda d: d['inverse_approche'][0].__setitem__(0, '0'),
            lambda d: d['modele']['contraintes'][0]['bt'].pop(),
            lambda d: d['modele']['jauges'].clear(),
            lambda d: d['centre'].__setitem__(0, '0.0'),
            lambda d: d.__setitem__('inclusion_stricte', True),
        ]
        for change in alterations:
            d = copy.deepcopy(self.document)
            change(d)
            with self.assertRaises(ValueError):
                verifier_document(d)


if __name__ == '__main__':
    unittest.main()
