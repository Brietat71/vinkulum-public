"""Couverture des expériences, altérations de verdict et rejeu des preuves."""
import copy,unittest
from archive_assemblage_complements import ARCHIVE,lire_json,verifier_rapport,verifier_archive

class ArchiveAssemblage(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.rows=lire_json(ARCHIVE,'campagne/bilan.json')
    def test_archive_et_neuf_paires_rejouees(self):
        self.assertEqual(verifier_archive(),dict(essais=24,champs_admis=24,controles_admis=12,paires_certificats_rejouees=9))
    def test_plan_et_pieces(self):
        with self.assertRaises(ValueError):verifier_rapport(self.rows[:-1])
        rows=copy.deepcopy(self.rows);rows[0]['taille_conservee']=7
        with self.assertRaises(ValueError):verifier_rapport(rows)
    def test_sources_budget_et_certification(self):
        for change in ('source','cpu','machine'):
            rows=copy.deepcopy(self.rows)
            if change=='source':rows[0]['sources_sha256']['assemblage_complements.py']='altéré'
            elif change=='cpu':rows[0]['environnement']['cpu']=[9]
            else:rows[0]['certification_machine_reponses']=True
            with self.assertRaises(ValueError):verifier_rapport(rows)
    def test_majorant_et_precision(self):
        for change in ('borne','precision'):
            rows=copy.deepcopy(self.rows)
            if change=='borne':rows[0]['retours_par_frequence'][0]['bornes']['masse']['relatives'][0]=.1
            else:rows[0]['juge']['maxima']['deformation']=.1
            with self.assertRaises(ValueError):verifier_rapport(rows)

if __name__=='__main__':unittest.main()
