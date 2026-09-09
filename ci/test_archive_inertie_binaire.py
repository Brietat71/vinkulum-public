"""Contre-épreuves sémantiques : les empreintes seules ne prouvent rien."""
import copy
from pathlib import Path
import tempfile
import unittest

from archive_inertie_binaire import (ARCHIVE,SONDES,lire_json,inventaire,
    verifier_archive,verifier_rapport,verifier_certificats)
from inertie_binaire_compile import compiler,BibliothequeInertie


class ArchiveInertieBinaire(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rapport=lire_json(ARCHIVE,'rapport.json')
        cls.temp=tempfile.TemporaryDirectory(prefix='vinkulum-archive-binaire-test-')
        cls.bib=BibliothequeInertie(compiler(Path(cls.temp.name)/'build'))

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def copie_premier(self):
        r=dict(self.rapport);r['essais']=list(r['essais'])
        r['essais'][0]=copy.deepcopy(r['essais'][0])
        return r,r['essais'][0]['resultat']

    def test_archive_complete_et_six_preuves_rejouees(self):
        self.assertEqual(verifier_archive(rejouer=False),dict(essais=36,admis=36,controles_admis=24,certificats_rejoues=0))
        self.assertEqual(verifier_certificats(self.rapport,ARCHIVE,self.bib),6)
        self.assertEqual(len(inventaire(SONDES)['fichiers']),53)

    def test_couverture_dupliquee_malgre_resultats_valides(self):
        r=dict(self.rapport);r['essais']=list(r['essais']);r['essais'][-1]=r['essais'][0]
        with self.assertRaisesRegex(ValueError,'dupliqué|imprévu'):verifier_rapport(r)
        r['essais']=r['essais'][:-1]
        with self.assertRaisesRegex(ValueError,'36 essais'):verifier_rapport(r)

    def test_pivot_falsifie_avec_signature_globale_conservee(self):
        r,e=self.copie_premier()
        e['certificat_complement']['preuve_kkt']['pivots'][0]['diagonal']=['0x1p+100','0x1p+101']
        with self.assertRaisesRegex(ValueError,'preuve différente'):verifier_certificats(r,ARCHIVE,self.bib)

    def test_permutation_valide_mais_differente_du_protocole(self):
        r,e=self.copie_premier();p=e['certificat_complement']['permutation_physique'];p[0],p[1]=p[1],p[0]
        with self.assertRaisesRegex(ValueError,'preuve différente'):verifier_certificats(r,ARCHIVE,self.bib)

    def test_binaire_different_malgre_preuve_compatible(self):
        r,e=self.copie_premier();e['bibliotheque_apres_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'bibliothèque ou arithmétique'):verifier_rapport(r)

    def test_fausse_certification_de_reponse_et_champ_change(self):
        r,e=self.copie_premier();e['retours_par_frequence'][0]['certification_machine']=True
        with self.assertRaisesRegex(ValueError,'certification machine'):verifier_rapport(r)
        r,e=self.copie_premier();e['champs']['sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'champs différents'):verifier_rapport(r)

    def test_minorant_et_normes_modifies(self):
        r,e=self.copie_premier();e['certificat_complement']['lambda_min']*=2
        with self.assertRaisesRegex(ValueError,'seuil ou inertie'):verifier_rapport(r)
        r,e=self.copie_premier();e['normes_candidates']['masse'][0][0]*=2
        with self.assertRaisesRegex(ValueError,'normes ou majorants différents'):verifier_rapport(r)


if __name__=='__main__':unittest.main()
