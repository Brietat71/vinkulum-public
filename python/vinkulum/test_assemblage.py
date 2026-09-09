"""Assemblage : références cinématiques fermées et atomicité des refus."""
import unittest

import numpy as np

from vinkulum import Noyau


J = np.eye(3).ravel().tolist()


def mobile(v=(0., 0., 0.), w=(0., 0., 0.)):
    n = Noyau([0., 0., 0.])
    n.corps('mobile', 1., J, [0., 0., 0.], v=list(v), w=list(w))
    return n


def commande(n, nom, taux):
    n.liaison(nom, None, 0, bloque_t=[0], bloque_r=[],
              cible_t=([1., 0., 0.], ('lineaire', [0., taux])))


class Assemblage(unittest.TestCase):
    def test_frontiere_tolerance_verifiee_exactement(self):
        # L'addition flottante arrondit le budget vers u, mais le budget
        # rationnel exact lui est strictement inférieur : il faut corriger.
        u = float(np.nextafter(1., 0.))
        tol = u * (1. - 64.*np.finfo(float).eps)
        self.assertLessEqual(u, tol + 64.*np.finfo(float).eps*u)
        n = mobile(v=(u, 0., 0.))
        commande(n, 'immobile', 0.)
        n.assemble(tol=tol)
        self.assertEqual(n.etat()[3][0][0], 0.)

    def test_derniere_correction_comptee(self):
        n = mobile()
        commande(n, 'translation', 1.)
        residu, iterations = n.assemble(t=1., iters=1)
        self.assertEqual(iterations, 1)
        self.assertLessEqual(residu, 1e-12)
        np.testing.assert_allclose(n.etat()[1][0], [1., 0., 0.], atol=1e-12)

    def test_commandes_incompatibles_restaurent_tout_etat(self):
        n = mobile(v=(3., 4., 5.), w=(.1, .2, .3))
        commande(n, 'lent', 1.)
        commande(n, 'rapide', 2.)
        s = n.etat_precis()
        s[1][0][0] = .25
        s[6][0][0] = 2.**-60
        n.pose_etat(*s)
        avant = n.etat_precis()
        for _ in range(2):
            with self.assertRaisesRegex(ValueError, 'vitesses.*incompatibles'):
                n.assemble()
            self.assertEqual(n.etat_precis(), avant)
        # Le mode explicitement limité aux positions conserve les vitesses.
        n.assemble(vitesses=False)
        self.assertLess(abs(n.phi()[0]), 1e-12)
        self.assertEqual(n.etat()[3:5], avant[3:5])

    def test_redondance_compatible_et_composantes_libres(self):
        for repetitions in (1, 2, 10):
            n = mobile(v=(3., 4., 5.), w=(.1, .2, .3))
            for i in range(repetitions):
                commande(n, str(i), 1.)
            n.assemble()
            np.testing.assert_allclose(n.etat()[3][0], [1., 4., 5.], atol=1e-12)
            np.testing.assert_allclose(n.etat()[4][0], [.1, .2, .3], atol=1e-12)

    def test_refus_preserve_horloge_schema_et_reactions(self):
        n = mobile()
        for i, debut in enumerate((1., 1.5)):
            n.liaison(str(i), None, 0, bloque_t=[0], bloque_r=[],
                      cible_t=([1., 0., 0.], ('table', [0., 0., debut, 0., 2., .25])))
        n.simule(.02, .01)
        avant = (n.etat_precis(), n.schema(), n.reactions(), n.stats())
        with self.assertRaisesRegex(ValueError, 'vitesses.*incompatibles'):
            n.assemble(t=2.)
        self.assertEqual((n.etat_precis(), n.schema(), n.reactions(), n.stats()), avant)

    def test_grand_bras_et_rotation_bloquee(self):
        # vy + L*wz = 0 et wz = 0 imposent vy = wz = 0, quel que soit L.
        for bras in (1., 1e5, 1e8):
            with self.subTest(bras=bras):
                n = mobile(v=(2., 1., 3.), w=(.1, .2, 0.))
                n.liaison('bras', None, 0, pa=[bras, 0., 0.],
                          bloque_t=[1], bloque_r=[2])
                n.assemble()
                np.testing.assert_allclose(n.etat()[3][0], [2., 0., 3.], atol=1e-12)
                np.testing.assert_allclose(n.etat()[4][0], [.1, .2, 0.], atol=1e-12)

    def test_correction_de_norme_minimale(self):
        # Projection orthogonale sur la droite vy + L*wz = 0.
        for bras in (0., 3., 1e5):
            n = mobile(v=(0., 2., 0.), w=(0., 0., 3.))
            n.liaison('bras', None, 0, pa=[bras, 0., 0.],
                      bloque_t=[1], bloque_r=[])
            n.assemble()
            vy = (2.*bras**2 - 3.*bras)/(1. + bras**2)
            wz = (3. - 2.*bras)/(1. + bras**2)
            np.testing.assert_allclose([n.etat()[3][0][1], n.etat()[4][0][2]],
                                       [vy, wz], rtol=1e-12, atol=1e-12)

    def test_position_avec_grand_bras(self):
        for bras in (1., 1e5, 1e8):
            n = mobile()
            n.liaison('bras', None, 0, pa=[bras, 0., 0.],
                      bloque_t=[1], bloque_r=[2])
            s = n.etat_precis()
            s[1][0][1] = .25
            n.pose_etat(*s)
            n.assemble(iters=3)
            self.assertLess(abs(n.etat()[1][0][1]), 1e-12)
            np.testing.assert_allclose(n.etat()[2][0], J, atol=1e-12)

    def test_derivation_temporelle_analytique(self):
        n = mobile()
        n.liaison('translation', None, 0, bloque_t=[0], bloque_r=[],
                  cible_t=([1., 0., 0.], ('lineaire', [-1e300, 1.])))
        n.assemble(t=1e300)
        self.assertEqual(n.etat()[3][0][0], 1.)
        np.testing.assert_allclose(n.phi_dot(t=1e300), [0.], atol=1e-12)
        n = mobile()
        n.liaison('rotation', None, 0, bloque_t=[], bloque_r=[0, 1, 2],
                  cible_r=([0., 0., 1.], ('lineaire', [0., 1e8])))
        n.assemble()
        np.testing.assert_allclose(n.etat()[4][0], [0., 0., 1e8], rtol=1e-14)
        np.testing.assert_allclose(n.phi_dot(), [0., 0., 0.], atol=1e-12)

    def test_pente_table_sans_moyenner_les_segments(self):
        # Convention de Loi::derivee : segment gauche au nœud intérieur.
        for date, pente in ((0., 1.), (1., 1.), (1.+1e-10, 3.), (2., 3.), (3., 0.)):
            n = mobile()
            n.liaison('table', None, 0, bloque_t=[0], bloque_r=[],
                      cible_t=([1., 0., 0.], ('table', [0., 0., 1., 1., 2., 4.])))
            n.assemble(t=date)
            self.assertAlmostEqual(n.etat()[3][0][0], pente, delta=1e-12)
            np.testing.assert_allclose(n.phi_dot(t=date), [0.], atol=1e-12)

    def test_non_holonome_ne_contraint_pas_la_position(self):
        n = mobile(v=(3., 2., 1.))
        n.liaison('vitesse_x', None, 0, bloque_t=[0], bloque_r=[], nh=True,
                  cible_t=([1., 0., 0.], ('lineaire', [0., 1.])))
        n.liaison('position_y', None, 0, bloque_t=[1], bloque_r=[])
        s = n.etat_precis()
        s[1][0] = [5., .25, 2.]
        n.pose_etat(*s)
        n.assemble(t=1.)
        np.testing.assert_allclose(n.etat()[1][0], [5., 0., 2.], atol=1e-12)
        # Au contact situé 2 m sous le centre, vx - 2*wy = 1.
        # La projection minimale de (vx, wy) = (3, 0) vaut (2.6, .8).
        np.testing.assert_allclose(n.etat()[3][0], [2.6, 0., 1.], atol=1e-12)
        np.testing.assert_allclose(n.etat()[4][0], [0., .8, 0.], atol=1e-12)

    def test_debordement_de_commande_restaure_la_pose(self):
        n = mobile()
        n.liaison('debordement', None, 0, bloque_t=[0], bloque_r=[],
                  cible_t=([1e150, 0., 0.], ('lineaire', [0., 1e200])))
        s = n.etat_precis()
        s[1][0][0] = .25
        n.pose_etat(*s)
        avant = n.etat_precis()
        with self.assertRaisesRegex(ValueError, 'non finie'):
            n.assemble()
        self.assertEqual(n.etat_precis(), avant)

    def test_diagnostic_non_holonome_sans_fausse_position(self):
        n = mobile(v=(1., 0., 0.))
        n.liaison('vitesse_x', None, 0, bloque_t=[0], bloque_r=[], nh=True,
                  cible_t=([1., 0., 0.], ('lineaire', [0., 1.])))
        # vx = 1 est satisfait ; x = t ne fait pas partie de ce modèle.
        self.assertEqual(n.controle(t=1.), [])

    def test_echec_position_restaure_la_pose(self):
        n = mobile(v=(3., 4., 5.))
        commande(n, 'lent', 1.)
        commande(n, 'rapide', 2.)
        avant = n.etat_precis()
        with self.assertRaisesRegex(ValueError, 'pas convergé'):
            n.assemble(t=1., iters=3)
        self.assertEqual(n.etat_precis(), avant)


if __name__ == '__main__':
    unittest.main()
