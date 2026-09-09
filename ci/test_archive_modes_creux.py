"""Relecture indépendante et contre-épreuves des archives modales."""
import copy
import gzip
import json
from pathlib import Path
import unittest

from qualifie_modes_creux import verifier, verifier_cas

DOSSIER = Path(__file__).resolve().parents[1]/'docs/bancs/modes-creux-0.17.0'


class ArchivesModales(unittest.TestCase):
    def test_grille_complete(self):
        self.assertEqual(len(verifier(DOSSIER)), 18)

    def test_alterations_physiques(self):
        cas = json.loads(gzip.decompress((DOSSIER/'articulee-4.json.gz').read_bytes()))
        for field in ('forme', 'reactions', 'valeur_propre'):
            with self.subTest(field=field):
                faux = copy.deepcopy(cas)
                mode = faux['resultat']['modes'][0]
                if field == 'valeur_propre':
                    mode[field] *= 1.01
                else:
                    mode[field][0] += 1.
                with self.assertRaises(ValueError):
                    verifier_cas(faux)

    def test_non_fini_et_fausse_certification(self):
        cas = json.loads(gzip.decompress((DOSSIER/'articulee-4.json.gz').read_bytes()))
        faux = copy.deepcopy(cas)
        faux['resultat']['modes'][0]['forme'][0] = float('nan')
        with self.assertRaises(ValueError):
            verifier_cas(faux)
        cas['resultat']['rang_certifie'] = True
        with self.assertRaises(ValueError):
            verifier_cas(cas)


if __name__ == '__main__':
    unittest.main()
