"""Contre-épreuves rationnelles du certificat contraint réellement calculé."""
from decimal import Decimal, ROUND_HALF_DOWN, localcontext
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, diags

from trace_complement_dirigee import certifier_complement
from vinkulum._ports.certificat_spectral import CertificationSpectraleImpossible
from vinkulum._ports.inverse_selectionnee import BudgetInverseSelectionnee
from test_trace_complement import (
    matrice, transpose, produit, difference, inverse_contrainte, restreindre,
    echelle, gram, identite, determinant)


def rationnel(a):
    return [[F(float(v)) for v in row] for row in a]


class TraceComplementDirigee(unittest.TestCase):
    def donnees(self):
        r = np.array([[2., 1., 0., 1.], [0., 3., 1., 0.],
                      [0., 0., 2., 1.], [0., 0., 0., 4.]])
        e = np.array([.5, 2., .25, 4.])
        d = np.vstack((r/e, [1e-8, -2e-8, 1e-8, 0.]))
        masse = np.array(gram(matrice([[1, 2, 0, 0], [0, 2, 1, 0],
                                      [0, 0, 3, 1], [0, 0, 0, 2]])), dtype=float)
        b = np.array([[1., 0.], [0., 1.], [2., -1.], [1., 3.]])
        z = matrice([[-2, -1], [1, -3], [1, 0], [0, 1]])
        return d, r, e, masse, b, z

    def contient(self, inf, sup, exact):
        self.assertLessEqual(F(Decimal(inf)), exact)
        self.assertGreaterEqual(F(Decimal(sup)), exact)

    def verifier_borne_exacte(self, d, m, z, resultat):
        q = difference(restreindre(gram(rationnel(d)), z),
                       echelle(restreindre(rationnel(m), z), F(resultat["lambda_min"])))
        # Les mineurs principaux dominants prouvent SPD par Sylvester.
        for j in range(1, len(q)+1):
            self.assertGreater(determinant([row[:j] for row in q[:j]]), 0)

    def test_contraintes_obliques_masse_couplee_et_facteur_perturbe(self):
        d, r, e, m, b, z = self.donnees()
        resultat = certifier_complement(d, r, e, m, b)
        _, tc, tau, chi = inverse_contrainte(gram(rationnel(r/e)), rationnel(m), rationnel(b))
        base = resultat["certificat_total"]
        self.contient(base["trace_inferieure_decimal"], base["trace_superieure_decimal"], tau)
        self.contient(resultat["correction_inferieure_decimal"], resultat["correction_superieure_decimal"], chi)
        self.contient(resultat["trace_complement_inferieure_decimal"], resultat["trace_complement_superieure_decimal"], tc)
        self.verifier_borne_exacte(d, m, z, resultat)
        self.assertEqual(resultat["rang_contraintes"], 2)
        self.assertTrue(resultat["certification_machine"])

    def test_conversion_finale_dirigee_et_base_du_meme_noyau(self):
        d, r, e, m, b, z = self.donnees()
        a = certifier_complement(d, r, e, m, b)
        bb = b@np.array([[2., 1.], [1., 1.]])
        autre = certifier_complement(d, r, e, m, bb)
        self.assertNotEqual(a["contraintes_sha256"], autre["contraintes_sha256"])
        for c in (a, autre):
            self.assertLessEqual(Decimal.from_float(c["lambda_min"]),
                                 Decimal(c["lambda_inferieur_decimal"]))
            self.verifier_borne_exacte(d, m, z, c)
        self.assertAlmostEqual(a["lambda_min"], autre["lambda_min"], places=13)

    def test_cancellation_et_precision_insuffisante_restent_sures(self):
        d = np.diag([2.**-50, 1., 2., 3.])
        m, b = np.eye(4), np.eye(4)[:, :1]
        exact = F(49, 36)
        bas = certifier_complement(d, d, np.ones(4), m, b, precision=16)
        haut = certifier_complement(d, d, np.ones(4), m, b, precision=80)
        for c in (bas, haut):
            self.contient(c["trace_complement_inferieure_decimal"], c["trace_complement_superieure_decimal"], exact)
            self.assertLessEqual(F(c["lambda_min"]), 1/exact)
            self.assertGreater(c["lambda_min"], 0.)
        self.assertLess(bas["lambda_min"], .01*haut["lambda_min"])
        self.assertAlmostEqual(haut["lambda_min"], float(1/exact), places=14)
        tau = F(2**100)+exact
        self.assertEqual(float(tau)-float(F(2**100)), 0.)

    def test_quasi_redondance_separe_rang_borne_et_bande(self):
        d, r, e, m, _, _ = self.donnees()
        b = np.array([[1., 1.], [0., 2.**-52], [0., 0.], [0., 0.]])
        z = matrice([[0, 0], [0, 0], [1, 0], [0, 1]])
        _, tc, _, chi = inverse_contrainte(
            gram(rationnel(r/e)), rationnel(m), rationnel(b))
        # Le rang exact vaut deux, mais 16 chiffres ne le démontrent pas.
        with self.assertRaisesRegex(CertificationSpectraleImpossible, "rang"):
            certifier_complement(d, r, e, m, b, precision=16)
        faible = certifier_complement(d, r, e, m, b, precision=40)
        precis = certifier_complement(d, r, e, m, b, precision=80)
        for c in (faible, precis):
            self.contient(c["trace_complement_inferieure_decimal"],
                          c["trace_complement_superieure_decimal"], tc)
            self.contient(c["correction_inferieure_decimal"],
                          c["correction_superieure_decimal"], chi)
            self.assertEqual(c["rang_contraintes"], 2)
            self.assertGreater(c["lambda_min"], 0.)
            self.assertLessEqual(F(c["lambda_min"]), 1/tc)
            self.verifier_borne_exacte(d, m, z, c)
        # Une cible de pulsation carrée choisie en fractions : obtenir un
        # certificat positif ne suffit pas à démontrer cette même bande.
        omega_cible_carre = 1/(2*tc)
        self.assertLess(F(faible["lambda_min"]), omega_cible_carre/100)
        self.assertGreater(F(precis["lambda_min"]), omega_cible_carre)
        self.assertLessEqual(F(Decimal(precis["trace_complement_superieure_decimal"])),
                             tc*(1+F(1, 10**10)))

    def test_reechelonnement_extreme_des_contraintes_preserve_la_borne(self):
        d, r, e, m, b, z = self.donnees()
        _, tc, _, chi = inverse_contrainte(
            gram(rationnel(r/e)), rationnel(m), rationnel(b))
        reference = certifier_complement(d, r, e, m, b, precision=80)
        for exposants in ((-400, 400), (400, -400)):
            with self.subTest(exposants=exposants):
                # Ces puissances de deux préservent exactement le noyau,
                # sans présumer la stabilité d'une normalisation flottante.
                bb = b*np.array([2.**p for p in exposants])
                c = certifier_complement(d, r, e, m, bb, precision=80)
                self.assertNotEqual(c["contraintes_sha256"], reference["contraintes_sha256"])
                self.contient(c["trace_complement_inferieure_decimal"],
                              c["trace_complement_superieure_decimal"], tc)
                self.contient(c["correction_inferieure_decimal"],
                              c["correction_superieure_decimal"], chi)
                self.assertLessEqual(F(c["lambda_min"]), 1/tc)
                self.assertGreater(F(c["lambda_min"]), F(9, 10)/tc)
                self.assertLessEqual(F(Decimal(c["trace_complement_superieure_decimal"])),
                                     tc*(1+F(1, 10**10)))
                self.verifier_borne_exacte(d, m, z, c)

    def test_contexte_decimal_ambiant_ne_change_pas_la_preuve(self):
        d, r, e, m, b, _ = self.donnees()
        a = certifier_complement(d, r, e, m, b)
        with localcontext() as ctx:
            ctx.prec, ctx.rounding = 6, ROUND_HALF_DOWN
            bb = certifier_complement(d, r, e, m, b)
        for cle in ("lambda_min", "lambda_inferieur_decimal", "trace_complement_superieure_decimal",
                    "correction_inferieure_decimal", "pivots_h_inferieurs_decimal"):
            self.assertEqual(a[cle], bb[cle], cle)

    def test_rang_insuffisant_et_mauvais_facteur_refusent(self):
        d, r, e, m, b, _ = self.donnees()
        for mauvais in (np.zeros((4, 1)), np.column_stack((b[:, 0], 2*b[:, 0]))):
            with self.assertRaisesRegex(CertificationSpectraleImpossible, "rang"):
                certifier_complement(d, r, e, m, mauvais)
        with self.assertRaisesRegex(CertificationSpectraleImpossible, "écart relatif"):
            certifier_complement(.01*np.eye(4), np.eye(4), np.ones(4), m, b)

    def test_donnees_invalides_et_budgets_explicites(self):
        d, r, e, m, b, _ = self.donnees()
        for mauvais in (b.astype(complex), np.full_like(b, np.nan), np.ones((3, 2)), np.eye(4), np.empty((4, 0))):
            with self.assertRaises(ValueError):
                certifier_complement(d, r, e, m, mauvais)
        with self.assertRaises(BudgetInverseSelectionnee):
            certifier_complement(d, r, e, m, b, budget_rectangulaire=2)
        a = certifier_complement(d, r, e, m, b)
        with self.assertRaises(BudgetInverseSelectionnee):
            certifier_complement(d, r, e, m, b, budget_operations=a["operations_decimal"]-1)

    def test_grande_chaine_sans_matrice_dense_avec_trace_rationnelle(self):
        n = 256
        r = diags([np.ones(n), .25*np.ones(n-1)], [0, 1], format="csr")
        b = np.zeros((n, 2)); b[0, 0] = b[1, 1] = 1.
        # Le complément donne un K tridiagonal Toeplitz : trace(K^-1)
        # = (d/da det(K))/det(K), avec deux récurrences scalaires exactes.
        a, carre = F(17, 16), F(1, 16)
        dm2, dm1, pm2, pm1 = F(1), a, F(0), F(1)
        for _ in range(2, n-1):
            dd = a*dm1-carre*dm2
            pp = dm1+a*pm1-carre*pm2
            dm2, dm1, pm2, pm1 = dm1, dd, pm1, pp
        exact = pm1/dm1
        def surveiller(original):
            def dense(matrice, *args, **kwargs):
                if min(matrice.shape) > 6:
                    raise AssertionError("grande matrice dense interdite")
                return original(matrice, *args, **kwargs)
            return dense
        with patch.object(csr_matrix, "toarray", surveiller(csr_matrix.toarray)), \
             patch.object(csc_matrix, "toarray", surveiller(csc_matrix.toarray)):
            c = certifier_complement(r, r, np.ones(n), diags(np.ones(n), format="csr"), b)
        self.contient(c["trace_complement_inferieure_decimal"], c["trace_complement_superieure_decimal"], exact)
        self.assertEqual(c["dimension_complement"], n-2)
        self.assertLessEqual(F(c["lambda_min"]), 1/exact)


if __name__ == "__main__":
    unittest.main()
