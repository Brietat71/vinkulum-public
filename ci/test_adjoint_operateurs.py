"""La comparaison d'adjoints exige des observables et des roues identiques."""
from copy import deepcopy
import gzip
import json
from pathlib import Path
import unittest

from mesure_adjoint_operateurs import analyse

ARCHIVE=Path(__file__).resolve().parents[1]/'docs/bancs/adjoint-operateurs-0.9.0-livraison-essais.json.gz'


class ComparaisonAdjoint(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original=json.loads(gzip.decompress(ARCHIVE.read_bytes()))

    def setUp(self):
        self.report=deepcopy(self.original)

    def candidate(self,rep=0):
        return next(r for r in self.report['runs'] if
                    (r['family'],r['elements'],r['variant'],r['rep'])==('pas',120,'produits',rep))

    def test_campagne_complete(self):
        self.assertTrue(analyse(self.report)['all_verified'])

    def test_grand_covecteur_ne_masque_pas_gradient_errone(self):
        self.candidate()['result']['values'][-1]+=1.
        self.assertFalse(analyse(self.report)['all_verified'])

    def test_echauffement_errone_refuse(self):
        self.candidate(-1)['result']['values'][-1]+=1.
        self.assertFalse(analyse(self.report)['all_verified'])

    def test_roue_differente_refusee(self):
        self.candidate()['result']['extension_sha256']='0'*64
        with self.assertRaises(ValueError):analyse(self.report)

    def test_resolution_adjointe_imprecise_refusee(self):
        self.candidate()['result']['adjoint_backward_error']=1e-7
        with self.assertRaises(ValueError):analyse(self.report)

    def test_essai_manquant_refuse(self):
        self.report['runs'].remove(self.candidate())
        with self.assertRaises(ValueError):analyse(self.report)

    def test_essai_echoue_conserve_sans_valider(self):
        r=self.candidate();r['ok']=False;r['error']='échec injecté'
        out=analyse(self.report)
        self.assertFalse(out['all_verified']);self.assertEqual(len(out['failures']),1)


if __name__=='__main__':unittest.main()
