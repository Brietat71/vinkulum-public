"""Oracles de Schur Decimal : physique, coefficients flottants et convergence."""
from decimal import Context, Decimal, ROUND_UP, getcontext, localcontext
import unittest

import numpy as np
from scipy.linalg import solve
from scipy.sparse import coo_matrix, csc_matrix, diags

from energie_ports import energie_chaine, energie_console
from modeles_ports import borne_lambda_trace_console, console, metrique_ports
from reference_ports_precision import schur_console, schur_facteur


def _norme_metric(delta, metric):
    w = np.linalg.solve(np.linalg.cholesky(metric).T, np.eye(len(metric)))
    return float(np.linalg.norm(w.T@delta@w, 2))


def _dense(k, m, ii, ss, omega):
    a = k-omega**2*m
    return (a[np.ix_(ss, ss)]-a[np.ix_(ss, ii)]
            @solve(a[np.ix_(ii, ii)], a[np.ix_(ii, ss)], assume_a="sym"))


class ReferencePortsPrecision(unittest.TestCase):
    def test_console_statique_et_signes(self):
        for n in (2, 8, 32):
            with self.subTest(n=n):
                model = console(n)
                reference = schur_console(model.metadata, 0.)
                metric = metrique_ports(model)
                self.assertLess(_norme_metric(reference-metric, metric), 3e-14)
                np.testing.assert_array_equal(reference, reference.T)
                self.assertLess(reference[1, 5], 0.)
                self.assertGreater(reference[2, 4], 0.)
                self.assertEqual(reference[1, 1], reference[2, 2])
                self.assertEqual(reference[4, 4], reference[5, 5])

    def test_frequences_petites_consoles_contre_dense(self):
        for n in (2, 4):
            model = console(n)
            energie = energie_console(model.metadata)
            metric = metrique_ports(model)
            k, m = model.k.toarray(), model.m.toarray()
            minimum = borne_lambda_trace_console(model.metadata)
            for rho in (0., .2, .6):
                with self.subTest(n=n, rho=rho):
                    omega = np.sqrt(rho*minimum)
                    reference = schur_console(model.metadata, omega)
                    fourni = schur_facteur(energie.d, energie.m, energie.interieur,
                                           energie.interface, omega)
                    dense = _dense(k, m, model.interieur, model.interface, omega)
                    self.assertLess(_norme_metric(reference-fourni, metric), 3e-13)
                    self.assertLess(_norme_metric(reference-dense, metric), 3e-11)

    def test_convergence_precision_128_elements_sans_matrice_dense(self):
        model = console(128)
        energie = energie_console(model.metadata)
        omega = np.sqrt(.5*borne_lambda_trace_console(model.metadata))
        a = schur_facteur(energie.d, energie.m, energie.interieur, energie.interface,
                          omega, dps=70, retour_decimal=True)
        b = schur_facteur(energie.d, energie.m, energie.interieur, energie.interface,
                          omega, dps=90, retour_decimal=True)
        with localcontext(Context(prec=100)):
            echelle = max(abs(x) for ligne in b for x in ligne)
            erreur = max(abs(a[i][j]-b[i][j]) for i in range(6) for j in range(6))
            self.assertLess(erreur/echelle, Decimal("1e-50"))
        self.assertIsInstance(a[0][0], Decimal)
        np.testing.assert_array_equal(np.asarray(a, float), np.asarray(b, float))
        physique = schur_console(model.metadata, omega)
        self.assertLess(_norme_metric(np.asarray(a, float)-physique, metrique_ports(model)), 3e-12)

    def test_coefficients_flottants_sont_des_entrees_exactes(self):
        n = 11
        energie = energie_chaine(n, 2., .017)
        reference = schur_facteur(energie.d, energie.m, energie.interieur, energie.interface,
                                  0., retour_decimal=True)[0][0]
        with localcontext(Context(prec=100)):
            racine_binaire = Decimal.from_float(float(np.sqrt(2.)))
            attendu = racine_binaire*racine_binaire/n
            # sqrt(2) arrondie au double puis mise au carré n'est pas le nombre 2.
            self.assertLess(abs(reference-attendu), Decimal("1e-65"))
            self.assertGreater(abs(reference-Decimal(2)/n), Decimal("1e-20"))

    def test_masse_demi_port_et_permutation(self):
        model = console(3)
        energie = energie_console(model.metadata)
        omega = 31.7
        reference = schur_facteur(energie.d, energie.m, energie.interieur, energie.interface, omega)
        masses = energie.m.diagonal()
        doubles = masses.copy()
        doubles[energie.interface] *= 2
        plein = schur_facteur(energie.d, diags(doubles, format="csc"), energie.interieur,
                              energie.interface, omega)
        delta = -omega**2*np.diag(masses[energie.interface])
        self.assertLess(_norme_metric((plein-reference)-delta, metrique_ports(model)), 2e-13)
        ordre = np.array([2, 0, 5, 1, 4, 3])
        permute = schur_facteur(energie.d, energie.m, energie.interieur[::-1],
                                energie.interface[ordre], omega)
        np.testing.assert_array_equal(permute, reference[np.ix_(ordre, ordre)])

    def test_blocs_couples_generiques(self):
        rng = np.random.default_rng(7743)
        n, p = 4, 2
        d = np.zeros((n*p, n*p))
        for j in range(n):
            d[j*p:(j+1)*p, j*p:(j+1)*p] = np.array([[2., .3], [-.1, 1.7]])
            if j:
                d[j*p:(j+1)*p, (j-1)*p:j*p] = .2*rng.normal(size=(p, p))
        m = np.diag(np.linspace(.2, .5, n*p))
        ii, ss = np.arange((n-1)*p), np.arange((n-1)*p, n*p)
        reference = schur_facteur(csc_matrix(d), csc_matrix(m), ii, ss, .8)
        attendu = _dense(d.T@d, m, ii, ss, .8)
        np.testing.assert_allclose(reference, attendu, rtol=2e-13, atol=2e-14)

    def test_contexte_decimal_local_et_rejets(self):
        model = console(2)
        with localcontext() as contexte:
            contexte.prec, contexte.rounding = 17, ROUND_UP
            a = schur_console(model.metadata, 0., retour_decimal=True)
            self.assertEqual(getcontext().prec, 17)
            self.assertEqual(getcontext().rounding, ROUND_UP)
        b = schur_console(model.metadata, 0., retour_decimal=True)
        self.assertEqual(a, b)
        for dps in (10, 30.5, True):
            with self.subTest(dps=dps), self.assertRaises(ValueError):
                schur_console(model.metadata, 0., dps=dps)
        with self.assertRaises(ValueError):
            schur_console(dict(model.metadata, formulation_poutre="milieu"), 0.)
        with self.assertRaises(ValueError):
            schur_console(model.metadata, -1.)
        d = np.eye(4)
        d[0, 2] = .1
        with self.assertRaisesRegex(ValueError, "non voisins"):
            schur_facteur(d, np.eye(4), [0, 1, 2], [3], 0.)
        with self.assertRaisesRegex(ValueError, "diagonale"):
            schur_facteur(np.eye(4), np.eye(4)+.1*np.ones((4, 4)), [0, 1, 2], [3], 0.)
        with self.assertRaisesRegex(ValueError, "partition"):
            schur_facteur(np.eye(4), np.eye(4), [1, 2, 3], [0], 0.)
        doublons = coo_matrix(([.1, .2, 1., 1., 1.],
                               ([0, 0, 1, 2, 3], [0, 0, 1, 2, 3])), shape=(4, 4))
        with self.assertRaisesRegex(ValueError, "canonique"):
            schur_facteur(doublons, np.eye(4), [0, 1, 2], [3], 0.)


if __name__ == "__main__":
    unittest.main()
