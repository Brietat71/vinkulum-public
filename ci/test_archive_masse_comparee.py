"""Contre-épreuves de la qualification et rejeu des données conservées."""
import copy
import unittest
from archive_masse_comparee import ARCHIVE,verifier_archive,verifier_rapport,lire_json

class ArchiveMasse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.rows=lire_json(ARCHIVE,'campagne/bilan.json')
    def test_rejeu_et_inventaire(self):
        self.assertEqual(verifier_archive(),dict(essais=24,champs_admis=24,controles_admis=12,paires_certificats_rejouees=3))
    def test_couverture(self):
        with self.assertRaisesRegex(ValueError,'plan'):verifier_rapport(self.rows[:-1])
        rows=self.rows.copy();rows[-1]=rows[0]
        with self.assertRaisesRegex(ValueError,'plan'):verifier_rapport(rows)
    def test_parametres_budget_et_certification(self):
        for change in ('alpha','cpu','certificat'):
            rows=copy.deepcopy(self.rows)
            if change=='alpha':rows[0]['parametres']['alpha']=.5
            elif change=='cpu':rows[0]['environnement']['cpu']=[9]
            else:rows[0]['certification_machine_reponses']=True
            with self.assertRaises(ValueError):verifier_rapport(rows)
    def test_precision_non_masquee(self):
        rows=copy.deepcopy(self.rows);rows[0]['juge']['maxima']['masse']=1e-2
        with self.assertRaises(ValueError):verifier_rapport(rows)
    def test_majorants_non_masques(self):
        rows=copy.deepcopy(self.rows);rows[0]['retours_par_frequence'][0]['bornes']['masse']['relatives'][0]=1e-2
        with self.assertRaises(ValueError):verifier_rapport(rows)

if __name__=='__main__':unittest.main()
