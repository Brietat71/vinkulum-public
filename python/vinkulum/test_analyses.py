"""Régressions indépendantes des analyses et diagnostics livrés en 0.6.0."""
import math
import unittest
from unittest.mock import patch

import numpy as np

from vinkulum import Noyau, convergence


def oscillateur(k=4., c=0., deplacement=0.):
    """Une masse unitaire en translation : x'' + c x' + k x = 0."""
    n = Noyau(g=[0, 0, 0])
    for i in range(2):
        n.corps(str(i), 1., np.eye(3).ravel().tolist(), [i, 0, 0])
    n.liaison('fixe', None, 0)
    n.liaison('glissiere', None, 1, pa=[1, 0, 0], bloque_t=[1, 2])
    K = np.zeros((12, 12))
    K[0, 0] = K[6, 6] = k
    K[0, 6] = K[6, 0] = -k
    n.superelement('ressort', [0, 1], K.ravel().tolist(), beta=c/k if c else 0.)
    if deplacement:
        etat = list(n.etat())
        etat[1][1][0] += deplacement
        n.pose_etat(*etat)
    return n


class AnalysesFiables(unittest.TestCase):
    def test_spectre_oscillateurs_analytiques(self):
        cas = [
            (4., 1., [-.5-1j*math.sqrt(3.75), -.5+1j*math.sqrt(3.75)], 'decroissance_spectrale'),
            (4., 5., [-4., -1.], 'decroissance_spectrale'),
            (4., 4., [-2., -2.], 'decroissance_spectrale'),
            (-4., 0., [-2., 2.], 'croissance_detectee'),
            (4., 0., [-2j, 2j], 'marginal_ou_indetermine'),
        ]
        for k, c, attendu, statut in cas:
            with self.subTest(k=k, c=c):
                n = oscillateur(k, c)
                etat = n.etat()
                bilan = n.bilan_stabilite()
                racines = n.spectre()
                z = np.array([complex(*v) for v in racines])
                # Appariement indépendant, y compris multiplicité critique.
                from scipy.optimize import linear_sum_assignment
                cout = abs(z[:, None] - np.array(attendu)[None, :])
                i, j = linear_sum_assignment(cout)
                self.assertEqual(len(z), 2)
                # À la racine double critique, l'erreur sur la racine varie
                # comme la racine carrée de celle de la tangente numérique.
                self.assertLess(cout[i, j].max(), 5e-6 if c == 4 else 1e-8)
                self.assertLess(np.max(abs(z*z + c*z + k)), 1e-8)
                self.assertEqual(racines, sorted(racines, key=lambda v: (-v[0], v[1])))
                self.assertEqual(bilan['racines'], racines)
                self.assertEqual(bilan['statut'], statut)
                self.assertAlmostEqual(bilan['abscisse_spectrale'], max(z.real))
                self.assertTrue(bilan['limites'])
                self.assertEqual(n.etat(), etat)
                anciens = n.modes_complexes()
                if c > 4 or k < 0:
                    self.assertEqual(anciens, [])
                elif c < 4 and k > 0:
                    self.assertEqual(len(anciens), 1)
                    f, zeta, sigma = anciens[0]
                    self.assertAlmostEqual(f, math.sqrt(k-c*c/4)/(2*math.pi), places=8)
                    self.assertAlmostEqual(zeta, c/(2*math.sqrt(k)), places=8)
                    self.assertAlmostEqual(sigma, -c/2, places=8)

    def test_spectre_rigide_vide_et_tol(self):
        n = Noyau(g=[0, 0, 0])
        for complet in (False, True):
            if complet:
                n.corps('libre', 1, np.eye(3).ravel().tolist(), [0, 0, 0])
                self.assertEqual(len(n.spectre()), 12)
                np.testing.assert_allclose(n.spectre(), 0, atol=1e-12)
                self.assertEqual(n.bilan_stabilite()['statut'], 'marginal_ou_indetermine')
                n.liaison('encastrement', None, 0)
            self.assertEqual(n.spectre(), [])
            b = n.bilan_stabilite()
            self.assertEqual(b['statut'], 'sans_ddl')
            self.assertIsNone(b['abscisse_spectrale'])
        n = oscillateur(-4)
        self.assertEqual(n.bilan_stabilite(tol=3)['statut'], 'marginal_ou_indetermine')
        for tol in (-1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                n.bilan_stabilite(tol=tol)

    def test_spectre_toupie_tournante(self):
        jt, ja, m, l, w = .001, .002, 1., .1, 300.
        n = Noyau(g=[0, 0, -9.81])
        n.corps('t', m, np.diag([jt, jt, ja]).ravel().tolist(), [0, 0, l], w=[0, 0, w])
        n.liaison('rotule', None, 0, bloque_r=[])
        racines = n.spectre()
        self.assertEqual(len(racines), 6)
        attendu = np.sort(abs(np.roots([jt+m*l*l, -ja*w, m*9.81*l])))
        calcule = sorted(im for re, im in racines if im > 1e-7)
        np.testing.assert_allclose(calcule, attendu, rtol=1e-6)
        self.assertLess(max(abs(re) for re, im in racines), 1e-7)

    def test_spectre_hors_domaine_et_dates(self):
        for option in ('contact', 'sphere', 'appariement', 'inflow'):
            n = oscillateur()
            if option == 'contact':
                n.contact('sol', 1, [0, 0, 0], .01)
            elif option == 'sphere':
                n.sphere(1, [0, 0, 0], .01)
            elif option == 'appariement':
                n.appariement()
            else:
                n.inflow([0, 0, 1], 1., .1)
            with self.assertRaisesRegex(RuntimeError, 'hors domaine'):
                n.spectre()
            bilan = n.bilan_stabilite()
            self.assertEqual(bilan['statut'], 'hors_domaine')
            self.assertIsNone(bilan['abscisse_spectrale'])
            self.assertEqual(bilan['racines'], [])
            self.assertGreater(len(bilan['limites']), 1)
            for f in (n.spectre, n.bilan_stabilite, n.controle, n.phi_dot):
                with self.assertRaises(RuntimeError):
                    f(t=math.nan)

    def test_diagnostic_cible_temporelle_sans_mutation(self):
        n = Noyau(g=[0, 0, 0])
        n.corps('mobile', 1, np.eye(3).ravel().tolist(), [0, 0, 0], v=[.5, 0, 0])
        n.liaison('guide', None, 0, cible_t=([1, 0, 0], ('lineaire', [0, .5])))
        # Obtenir un état de schéma et des multiplicateurs effectivement calculés.
        n.simule(.2, .01)
        avant = (n.etat(), n.schema(), n.reactions(), n.stats())
        self.assertEqual(n.controle(), [])
        np.testing.assert_allclose(n.phi_dot(), 0, atol=1e-9)
        self.assertTrue(any('contraintes non satisfaites' in x for x in n.controle(t=0)))
        self.assertEqual(n.controle(t=.2), n.controle())
        self.assertEqual(avant, (n.etat(), n.schema(), n.reactions(), n.stats()))
        e = list(n.etat())
        e[3][0][0] += .2
        n.pose_etat(*e)
        self.assertAlmostEqual(n.phi_dot()[0], .2, places=8)
        self.assertTrue(any('vitesses incompatibles' in x for x in n.controle()))

    def test_diagnostic_date_explicite_et_tables(self):
        n = Noyau(g=[0, 0, 0])
        n.corps('mobile', 1, np.eye(3).ravel().tolist(), [0, 0, 0])
        n.liaison('guide', None, 0, cible_t=([1, 0, 0], ('table', [0, 0, 1, .5, 2, 2])))
        e = list(n.etat())
        e[1][0][0] = 1.25
        e[3][0][0] = 1.5
        n.pose_etat(*e)  # horloge laissée à zéro
        self.assertEqual(n.controle(t=1.5), [])
        np.testing.assert_allclose(n.phi_dot(t=1.5), 0, atol=1e-9)
        self.assertAlmostEqual(n.phi_dot(t=.5)[0], 1., places=8)
        self.assertEqual(n.t(), 0.)
        # Une table est constante hors domaine, même à la dernière date
        # finie. La dérivée analytique ne requiert plus de date t + epsilon.
        extreme = float.fromhex('0x1.fffffffffffffp+1023')
        np.testing.assert_allclose(n.phi_dot(t=extreme), [1.5, 0., 0., 0., 0., 0.])
        self.assertTrue(any('vitesses incompatibles' in x for x in n.controle(t=extreme)))

    def test_diagnostic_corps_interactions(self):
        for interaction in ('super', 'effort', 'couple', 'sphere', 'contact', 'aucune'):
            n = Noyau(g=[0, 0, 0])
            for i in range(2):
                n.corps(str(i), 1, np.eye(3).ravel().tolist(), [i, 0, 0])
            if interaction == 'super':
                n.superelement('s', [0, 1], np.zeros((12, 12)).ravel().tolist())
            elif interaction == 'effort':
                for i in range(2):
                    n.effort(i, [1, 0, 0], [0, 0, 0])
            elif interaction == 'couple':
                n.couple('c', 0, 1, [0, 0, 1], ('constant', [1]))
            elif interaction == 'sphere':
                for i in range(2):
                    n.sphere(i, [0, 0, 0], .01)
            elif interaction == 'contact':
                n.contact('c', 0, [0, 0, 0], .01, b=1, rayon_b=.01)
            orphelins = [s for s in n.controle() if 'aucune liaison' in s]
            self.assertEqual(len(orphelins), 2 if interaction == 'aucune' else 0)

    def test_convergence_reference_analytique(self):
        monte = lambda: oscillateur(deplacement=.1)
        observable = lambda n: [n.pose(1)[0][0]-1]
        table = convergence.etude(monte, 1., .04, grandeur=observable)
        self.assertTrue(1.8 < table[0][2] < 2.2, table)
        tol = 1e-5
        self.assertFalse(convergence.verdict(monte, 1., .1, tol, grandeur=observable)[0])
        self.assertTrue(convergence.verdict(monte, 1., .002, tol, grandeur=observable)[0])
        n = monte()
        n.simule(1., .0005)
        self.assertLess(abs(observable(n)[0] - .1*math.cos(2)), 1e-7)
        irregulier = convergence.etude(monte, 1., .04, facteurs=(1, 2, 5), grandeur=observable)
        self.assertTrue(all(o is None for _, _, o in irregulier))
        self.assertGreater(irregulier[0][1], 0)

    def test_convergence_entrees_invalides(self):
        for facteurs in ((), (1,), (1, 1), (2, 1), (0, 1), (1, math.nan), (1, math.inf), ((1, 2), (3, 4))):
            with self.subTest(facteurs=facteurs), self.assertRaises(ValueError):
                convergence.etude(oscillateur, 1., .01, facteurs=facteurs)
        for h in (0, -1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                convergence.etude(oscillateur, 1., h)
        for fin in (0, -1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                convergence.etude(oscillateur, fin, .01)
        for tol in (-1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                convergence.verdict(oscillateur, 1., .01, tol)
        with self.assertRaisesRegex(ValueError, 'raffinés'):
            convergence.etude(oscillateur, 1., 5e-324)
        with self.assertRaisesRegex(ValueError, 'raffinés'):
            convergence.etude(oscillateur, 1., 1e-323, facteurs=(1, 1.1))

    def test_convergence_observables_et_modeles(self):
        n = oscillateur()
        with self.assertRaisesRegex(ValueError, 'indépendant'):
            convergence.etude(lambda: n, .02, .01)
        compteur = 0
        def decale():
            nonlocal compteur
            n = oscillateur()
            e = list(n.etat())
            e[0] = .001 * compteur
            n.pose_etat(*e)
            compteur += 1
            return n
        with self.assertRaisesRegex(ValueError, 'même instant'):
            convergence.etude(decale, .02, .01)
        for observable in (lambda n: [], lambda n: [math.nan], lambda n: [math.inf]):
            with self.assertRaisesRegex(ValueError, 'observable'):
                convergence.etude(oscillateur, .02, .01, grandeur=observable)
        dimensions = iter(([1], [1, 2]))
        with self.assertRaisesRegex(ValueError, 'dimensions'):
            convergence.etude(oscillateur, .02, .01, grandeur=lambda n: next(dimensions))
        extreme = iter(([1e308], [-1e308], [0]))
        with self.assertRaisesRegex(ValueError, 'écart'):
            convergence.etude(oscillateur, .02, .01, grandeur=lambda n: next(extreme))
        tampon = np.zeros(1)
        def partage(n):
            tampon[0] += 1
            return tampon
        table = convergence.etude(oscillateur, .02, .01, grandeur=partage)
        self.assertEqual([ligne[1] for ligne in table], [2., 1.])

    def test_convergence_simulation_incomplete(self):
        # Double de la frontière native pour exercer la garde de progression,
        # le noyau réel disposant déjà de ses propres refus et restaurations.
        class Incomplet:
            def t(self):
                return 0.
            def simule(self, *args, **kwargs):
                pass
        with patch.object(convergence, 'Noyau', Incomplet):
            with self.assertRaisesRegex(RuntimeError, 'date finale'):
                convergence.etude(Incomplet, 1., .01)


if __name__ == '__main__':
    unittest.main()
