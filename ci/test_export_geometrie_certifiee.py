"""Contre-épreuves du prototype sur la géométrie réellement stockée par Rust."""
from fractions import Fraction as F
import math
import unittest

from vinkulum import Noyau
from prototype_assemblage_certifie import depuis_noyau, inclusion, proposer_inverse, bornes_pose, evaluation


J = [1., 0., 0., 0., 1., 0., 0., 0., 1.]


def pendule(longueur=1.):
    n = Noyau([0., 0., -9.81])
    n.corps('mobile', 1., J, [0., 0., -longueur])
    n.liaison('pivot', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 2])
    return n


class ExportGeometrique(unittest.TestCase):
    def test_inclusion_et_etat_inchange_avec_translation_compensee(self):
        n = pendule()
        etat = n.etat_precis()
        etat[6][0][0] = 2.**-70
        n.pose_etat(*etat)
        avant = n.etat_precis()
        m, x, photo = depuis_noyau(n, [(5, 0)])
        self.assertEqual(x[0], F(1, 2**70))
        self.assertEqual(photo['lignes_originales'], 5)
        preuve = inclusion(m, x, [F(1, 10**6)]*7, proposer_inverse(m, x))
        self.assertLess(preuve['contraction'], 1)
        self.assertEqual(n.etat_precis(), avant)

    def test_distance_rotation_inclut_conversion_de_la_matrice(self):
        n = pendule()
        m, x, photo = depuis_noyau(n, [(5, 0)])
        rayons = [F(1, 10**6)]*7
        base = bornes_pose(photo, x, rayons)[0]
        self.assertTrue(all(v == F(1, 10**6) for v in base['position']))
        # Une matrice d'origine différente doit contribuer à la distance,
        # même si le centre quaternion et la boîte restent identiques.
        photo['poses'][0]['rotation'][0][0] += 1e-3
        autre = bornes_pose(photo, x, rayons)[0]
        self.assertGreater(autre['rotation_frobenius_carree'], F(1, 10**6))
        self.assertGreater(autre['rotation_frobenius_carree'], base['rotation_frobenius_carree'])

    def test_echelles_de_longueur(self):
        for longueur in (1e-6, 1., 1e6):
            with self.subTest(longueur=longueur):
                n = pendule(longueur)
                m, x, _ = depuis_noyau(n, [(5, 0)])
                rayons = [F(longueur)*F(1, 10**8)]*3 + [F(1, 10**8)]*4
                preuve = inclusion(m, x, rayons, proposer_inverse(m, x))
                self.assertLess(preuve['contraction'], 1)

    def test_double_pendule(self):
        n = Noyau([0., 0., -9.81])
        n.corps('a', 1., J, [0., 0., -.5])
        n.corps('b', 1., J, [0., 0., -1.5])
        n.liaison('sol', None, 0, bloque_r=[0, 2])
        n.liaison('coude', 0, 1, pa=[0., 0., -.5], bloque_r=[0, 2])
        m, x, photo = depuis_noyau(n, [(5, 0), (12, 0)])
        self.assertEqual(photo['lignes_originales'], 10)
        preuve = inclusion(m, x, [F(1, 10**6)]*14, proposer_inverse(m, x))
        self.assertLess(preuve['contraction'], 1)

    def test_trilateration_distance(self):
        n = Noyau([0., 0., 0.])
        n.corps('point', 1., J, [0., 0., 1.])
        for k, p in enumerate(([0., 0., 0.], [1., 0., 0.], [0., 1., 0.])):
            n.distance(str(k), None, 0, p, [0., 0., 0.], l=1. if k == 0 else math.sqrt(2.))
        m, x, _ = depuis_noyau(n, [(4, 0), (5, 0), (6, 0)])
        preuve = inclusion(m, x, [F(1, 10**6)]*7, proposer_inverse(m, x))
        # Le x exact vient de la différence des deux premières sphères.
        exact_x = (2-F(math.sqrt(2.))**2)/2
        self.assertLess(abs(exact_x-x[0]), preuve['rayons'][0])

    def test_residus_et_derives_contre_le_noyau_hors_equilibre(self):
        import numpy as np
        from scipy.spatial.transform import Rotation
        rng = np.random.default_rng(0xCE4705)
        n = Noyau([0., 0., -9.81])
        n.corps('mobile', 1., J, [0., 0., -1.])
        repere = Rotation.from_rotvec([.2, -.3, .4]).as_matrix()
        n.liaison('pivot', None, 0, ra=repere.ravel().tolist(), bloque_r=[0, 2])
        for _ in range(32):
            etat = n.etat_precis()
            etat[1][0] = rng.normal(size=3).tolist()
            etat[2][0] = Rotation.from_rotvec(rng.normal(size=3)).as_matrix().ravel().tolist()
            etat[3][0] = rng.normal(size=3).tolist()
            etat[4][0] = rng.normal(size=3).tolist()
            n.pose_etat(*etat)
            m, x, _ = depuis_noyau(n, [(5, 0)])
            f, j, _ = evaluation(m, x)
            w, a, b, c = x[3:]
            vx, vy, vz = map(F, etat[4][0])
            qdot = [-(vx*a+vy*b+vz*c)/2,
                    (w*vx+vy*c-vz*b)/2,
                    (w*vy+vz*a-vx*c)/2,
                    (w*vz+vx*b-vy*a)/2]
            dx = list(map(F, etat[3][0])) + qdot
            derivees = [sum(v.bas*d for v, d in zip(row, dx, strict=True)) for row in j[:5]]
            # Tolérance fixée pour la confrontation de deux évaluations
            # arrondies ; l'inclusion du certificat reste strictement exacte.
            np.testing.assert_allclose([float(v.bas) for v in f[:5]], n.phi(), rtol=0., atol=1e-12)
            np.testing.assert_allclose(list(map(float, derivees)), n.phi_dot(), rtol=0., atol=1e-12)

    def test_contrainte_redondante_conservee_et_refusee(self):
        n = pendule()
        n.liaison('double', None, 0, bloque_t=[0], bloque_r=[])
        m, x, photo = depuis_noyau(n, [(5, 0)])
        self.assertEqual(photo['lignes_originales'], 6)
        self.assertEqual(len(m['contraintes']), 2)
        avant = n.etat_precis()
        with self.assertRaisesRegex(ValueError, 'système carré'):
            proposer_inverse(m, x)
        self.assertEqual(n.etat_precis(), avant)

    def test_loi_imposee_refusee_sans_mutation(self):
        n = pendule()
        n.liaison('loi', None, 0, bloque_t=[0], bloque_r=[],
                  cible_t=([1., 0., 0.], ('lineaire', [0., 1.])))
        avant = n.etat_precis()
        with self.assertRaisesRegex(ValueError, 'domaine'):
            depuis_noyau(n, [(5, 0)])
        self.assertEqual(n.etat_precis(), avant)


if __name__ == '__main__':
    unittest.main()
