"""Contrat public, relecture autonome et conservation de l'état natif."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from vinkulum import Noyau
from vinkulum.certification import certifier_assemblage, verifier_certificat, CertificationImpossible


J = [1., 0., 0., 0., 1., 0., 0., 0., 1.]


def pendule():
    n = Noyau([0., 0., -9.81])
    n.corps('mobile', 1., J, [0., 0., -1.])
    n.liaison('pivot', None, 0, bloque_r=[0, 2])
    return n


class CertificationAssemblage(unittest.TestCase):
    def setUp(self):
        self.n = pendule()
        self.rayons = [F(1, 10**6)]*7

    def certifie(self, **kw):
        return certifier_assemblage(self.n, jauges=kw.get('jauges', [(5, 0)]),
                                   rayons=kw.get('rayons', self.rayons))

    def test_etat_et_parametres_inchanges(self):
        avant = self.n.etat_precis()
        doc = self.certifie()
        self.assertEqual(self.n.etat_precis(), avant)
        retour = verifier_certificat(json.loads(json.dumps(doc)))
        self.assertTrue(retour['existence_unique_locale'])
        self.assertEqual(retour['bornes_pose'][0]['position'], self.rayons[:3])
        self.assertGreater(retour['bornes_pose'][0]['rotation_frobenius_carree'], 0)
        self.assertEqual(self.rayons, [F(1, 10**6)]*7)

    def test_bornes_annoncees_et_pose_initiale_alterees_refusees(self):
        document = self.certifie()
        for champ in ('position', 'rotation'):
            d = copy.deepcopy(document)
            if champ == 'position':
                d['bornes_pose'][0]['position'][0] = '0'
            else:
                d['bornes_pose'][0]['rotation_frobenius_carree'] = '0'
            with self.assertRaisesRegex(ValueError, 'distance'):
                verifier_certificat(d)
        d = copy.deepcopy(document)
        d['pose_initiale'][0]['rotation'][0][0] = '2'
        with self.assertRaisesRegex(ValueError, 'distance'):
            verifier_certificat(d)

    def test_entrees_hors_contrat_refusees_sans_mutation(self):
        avant = self.n.etat_precis()
        for rayons in ([0]*7, [-1]*7, [float('inf')]*7, [float('nan')]*7,
                       [True]*7, ['1']*7, [1]*6, [1+0j]*7):
            with self.subTest(rayons=rayons):
                with self.assertRaises(CertificationImpossible):
                    self.certifie(rayons=rayons)
        for jauges in ([], [(5, 0), (5, 0)], [(True, 0)], [(99, 0)], [(5, '0')]):
            with self.subTest(jauges=jauges):
                with self.assertRaises(CertificationImpossible):
                    self.certifie(jauges=jauges)
        self.assertEqual(self.n.etat_precis(), avant)

    def test_refus_de_preuve_ne_repare_pas_la_pose(self):
        etat = self.n.etat_precis()
        etat[1][0][0] += .01
        self.n.pose_etat(*etat)
        avant = self.n.etat_precis()
        with self.assertRaises(CertificationImpossible):
            self.certifie()
        self.assertEqual(self.n.etat_precis(), avant)

    def test_verificateur_public_sans_import_du_paquet(self):
        import vinkulum._verification_lineaire as module
        with tempfile.TemporaryDirectory() as tmp:
            document = Path(tmp)/'geometrie.json'
            document.write_text(json.dumps(self.certifie()))
            r = subprocess.run([sys.executable, '-S', module.__file__, str(document)],
                               capture_output=True, text=True, check=True)
            self.assertIn('Certificat géométrique vérifié', r.stdout)

    def test_reperes_tournes_et_export_non_modifie(self):
        import numpy as np
        from scipy.spatial.transform import Rotation
        r = Rotation.from_rotvec([.2, -.3, .4]).as_matrix()
        n = Noyau([0., 0., -9.81])
        n.corps('tourne', 1., J, (r @ np.array([0., 0., -1.])).tolist(), rot=r.ravel().tolist())
        n.liaison('pivot', None, 0, ra=r.ravel().tolist(), bloque_r=[0, 2])
        avant = n.etat_precis()
        qy = float(Rotation.from_matrix(r).as_quat()[1])
        d = certifier_assemblage(n, jauges=[(5, qy)], rayons=self.rayons)
        self.assertEqual(n.etat_precis(), avant)
        self.assertTrue(verifier_certificat(d)['existence_unique_locale'])
        for i in range(3):
            for j in range(3):
                self.assertEqual(F(d['modele']['contraintes'][0]['ra'][i][j]), F(float(r[i,j])))
        mauvais = copy.deepcopy(d)
        rb = mauvais['modele']['contraintes'][0]['rb']
        rb[1][2] = str(F(rb[1][2])+F(1,10))
        with self.assertRaises(ValueError):
            verifier_certificat(mauvais)
        del d['modele']['contraintes'][0]['ra']
        with self.assertRaises(ValueError):
            verifier_certificat(d)

    def test_schemas_anterieurs_restent_disponibles(self):
        from vinkulum.certification import certifier_systeme
        d = certifier_systeme([[2.]], [4.])
        self.assertEqual(verifier_certificat(d)['dimension'], 1)


if __name__ == '__main__':
    unittest.main()
