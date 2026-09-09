"""Témoins algébriques et exacts de la minoration spectrale sélectionnée."""
from decimal import Decimal, ROUND_HALF_DOWN, localcontext
from fractions import Fraction
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import cholesky, eigvalsh
from scipy.sparse import csc_matrix, coo_matrix, csr_matrix, diags, lil_matrix

from inverse_selectionnee import InverseSelectionnee, BudgetInverseSelectionnee
from certificat_spectral import (CertificationSpectraleImpossible, _Intervalles,
                                 certifier_spectral, certifier_spectre)


def _console_native(n):
    """D réellement exporté du noyau, racine supprimée, bout conservé."""
    from vinkulum import Noyau
    from condensation_energie import CondensationEnergie
    rayon, young, rho, poisson = .005, 70e9, 2700., .3
    aire, moment = np.pi*rayon**2, np.pi*rayon**4/4
    polaire, pas = 2*moment, 1./n
    cisaillement = young/(2*(1+poisson))
    noyau = Noyau([0., 0., 0.])
    ids = []
    for j in range(n+1):
        poids = .5 if j in (0, n) else 1.
        ids.append(noyau.corps(
            f"section{j}", rho*aire*pas*poids,
            (rho*pas*poids*np.diag([polaire, moment, moment])).ravel().tolist(),
            [j*pas, 0., 0.]))
    noyau.liaison("encastrement", None, ids[0])
    for j in range(n):
        noyau.poutre(f"poutre{j}", ids[j], ids[j+1], young*aire,
                     cisaillement*aire, cisaillement*polaire, young*moment,
                     formulation="integree")
    export = noyau.facteurs_materiels_poutres()
    def csc(t):
        nr, nc, ptr, indices, valeurs = t
        return csc_matrix((valeurs, indices, ptr), shape=(nr, nc))
    d, masse = csc(export["d"])[:, 6:], csc(export["masse"])[6:, 6:]
    interieur, interface = np.arange(6*(n-1)), np.arange(6*(n-1), 6*n)
    qr = CondensationEnergie(d, interieur, interface, np.eye(6))
    return qr, masse[interieur][:, interieur]


def _fractions(a):
    return [[Fraction(float(x)) for x in row] for row in a]


def _inverse_exacte(a):
    n = len(a)
    x = [list(row)+[Fraction(i == j) for j in range(n)] for i, row in enumerate(a)]
    for i in range(n):
        pivot = next(j for j in range(i, n) if x[j][i])
        x[i], x[pivot] = x[pivot], x[i]
        valeur = x[i][i]
        x[i] = [z/valeur for z in x[i]]
        for j in range(n):
            if i != j:
                facteur = x[j][i]
                x[j] = [a-facteur*b for a, b in zip(x[j], x[i])]
    return [row[n:] for row in x]


def _gram_exact(a):
    return [[sum(row[i]*row[j] for row in a) for j in range(len(a[0]))]
            for i in range(len(a[0]))]


def _pivots_ldl_exacts(a):
    a = [row.copy() for row in a]
    pivots = []
    for k in range(len(a)):
        pivot = a[k][k]
        pivots.append(pivot)
        if pivot <= 0:
            return pivots
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                a[i][j] -= a[i][k]*a[k][j]/pivot
    return pivots


class TestInverseSelectionnee(unittest.TestCase):
    def test_spd_rectangulaire_masse_blocs_echelles(self):
        rng = np.random.default_rng(501)
        d = rng.normal(size=(21, 9))
        k = d.T@d
        e = np.geomspace(.03, 20., 9)
        masse = np.diag(np.linspace(1., 3., 9))
        ids = [0, 4, 8]
        bloc = np.array([[2., .2, -.1], [.2, 1.5, .3], [-.1, .3, 2.5]])
        masse[np.ix_(ids, ids)] = bloc
        r = cholesky(e[:, None]*k*e[None, :], lower=False)
        obj = InverseSelectionnee(r, masse, e)
        s = np.linalg.inv(r.T@r)
        for i, row in enumerate(obj.valeurs):
            for j, v in row.items():
                self.assertAlmostEqual(v/s[i, j], 1., places=10)
        trace = np.trace(masse@np.linalg.inv(k))
        self.assertAlmostEqual(obj.trace_masse/trace, 1., places=12)
        self.assertLessEqual(obj.lambda_trace, eigvalsh(k, masse)[0])
        self.assertFalse(obj.bilan()["certification_machine"])
        self.assertEqual(obj.bilan()["taille_bloc_masse_max"], 3)

    def test_fermeture_non_chordale_et_masse_hors_motif(self):
        n = 9
        r = 2*np.eye(n)
        r[0, 1:] = np.linspace(-.2, .3, n-1)
        obj = InverseSelectionnee(r, np.ones(n))
        self.assertGreater(obj.bilan()["coefficients_fermeture"], 0)
        exact = np.linalg.inv(r.T@r)
        for i in range(n):
            for j in range(i, n):
                self.assertAlmostEqual(obj.entree(i, j), exact[i, j], places=14)
        r = diags((np.full(n, 2.), np.full(n-1, -.3)), (0, 1)).toarray()
        masse = np.eye(n)
        masse[0, -1] = masse[-1, 0] = .2
        obj = InverseSelectionnee(r, masse)
        self.assertGreater(obj.bilan()["coefficients_fermeture"], 0)
        self.assertAlmostEqual(obj.trace_masse, np.trace(masse@np.linalg.inv(r.T@r)), places=13)

    def test_refus_budgets_et_entree_non_selectionnee(self):
        r = np.eye(20)
        r[0, 1:] = .1
        with self.assertRaises(BudgetInverseSelectionnee) as cm:
            InverseSelectionnee(r, np.ones(20), budget_coefficients=50)
        self.assertEqual(cm.exception.bilan["phase"], "fermeture")
        with self.assertRaises(BudgetInverseSelectionnee):
            InverseSelectionnee(r, np.ones(20), budget_operations=3)
        diagonal = InverseSelectionnee(np.eye(3), np.ones(3))
        with self.assertRaises(KeyError):
            diagonal.entree(0, 2)

    def test_grande_bande_sans_inverse_dense_ni_resolutions(self):
        n = 4096
        r = diags((np.full(n, 2.), np.full(n-1, -.3)), (0, 1), format="csr")
        masse = np.geomspace(.7, 1.7, n)
        attendu = np.empty(n)
        attendu[-1] = .25
        for i in range(n-2, -1, -1):
            attendu[i] = .25+.0225*attendu[i+1]
        with patch("numpy.linalg.inv", side_effect=AssertionError("inverse dense")), \
             patch("numpy.linalg.solve", side_effect=AssertionError("résolution dense")), \
             patch.object(csr_matrix, "toarray", side_effect=AssertionError("densification")):
            obj = InverseSelectionnee(r, masse)
        self.assertEqual(obj.coefficients, 2*n-1)
        self.assertLess(obj.bilan()["operations_total"], 30*n)
        self.assertAlmostEqual(obj.trace_masse/float(masse@attendu), 1., places=14)

    def test_permutation_physique_conserve_trace(self):
        rng = np.random.default_rng(940)
        d = rng.normal(size=(18, 10))
        k, masse = d.T@d, np.linspace(.5, 2., 10)
        avant = InverseSelectionnee(cholesky(k), masse)
        p = rng.permutation(10)
        apres = InverseSelectionnee(cholesky(k[np.ix_(p, p)]), masse[p])
        self.assertAlmostEqual(avant.trace_masse/apres.trace_masse, 1., places=13)

    def test_masses_et_facteurs_hors_contrat(self):
        for masse, echelles in ((np.array([1.+2.j, 1.]), [1., 1.]),
                                ([1., 1.], np.array([1.+2.j, 1.]))):
            with self.assertRaises(ValueError):
                certifier_spectre(np.eye(2), np.eye(2), echelles, masse)
        with self.assertRaises(ValueError):
            InverseSelectionnee(np.eye(2), [1., -1.])
        with self.assertRaises(ValueError):
            InverseSelectionnee(np.eye(2), [[1., 2.], [2., 1.]])
        with self.assertRaises(ValueError):
            InverseSelectionnee([[1., 0.], [1., 1.]], [1., 1.])
        m = np.eye(7)+.1*np.ones((7, 7))
        with self.assertRaises(BudgetInverseSelectionnee):
            InverseSelectionnee(np.eye(7), m)


class TestCertificatSpectral(unittest.TestCase):
    def test_fraction_trace_et_minorant_original_perturbe(self):
        r = np.array([[2., .3, -.2], [0., 1.3, .4], [0., 0., .9]])
        e = np.array([.7, 1.2, 2.1])
        masse = np.array([[1.2, .1, 0.], [.1, 2., -.2], [0., -.2, 1.7]])
        d = r/e[None, :]
        d[0, 1] += 3e-5
        inv = InverseSelectionnee(r, masse, e)
        cert = certifier_spectral(d, inv)
        rf, ef, mf = _fractions(r), [Fraction(float(x)) for x in e], _fractions(masse)
        qr_exact = _gram_exact([[row[j]/ef[j] for j in range(3)] for row in rf])
        sqr = _inverse_exacte(qr_exact)
        tau = sum(mf[i][j]*sqr[j][i] for i in range(3) for j in range(3))
        self.assertLessEqual(Fraction(Decimal(cert["trace_inferieure_decimal"])), tau)
        self.assertGreaterEqual(Fraction(Decimal(cert["trace_superieure_decimal"])), tau)
        lam = Fraction(cert["lambda_min"])
        kd = _gram_exact(_fractions(d))
        difference = [[kd[i][j]-lam*mf[i][j] for j in range(3)] for i in range(3)]
        self.assertTrue(all(x > 0 for x in _pivots_ldl_exacts(difference)))
        self.assertLessEqual(lam, Fraction(Decimal(cert["lambda_inferieur_decimal"])))
        self.assertTrue(cert["certification_machine"])
        self.assertTrue(cert["minoration_D_original_etablie"])

    def test_contexte_ambiant_hostile_et_racine_dirigee(self):
        r = np.array([[1.3, -.7], [0., .9]])
        inv = InverseSelectionnee(r, [1.1, 2.3], [.8, 1.4])
        d = r/inv.echelles[None, :]
        avant = certifier_spectral(d, inv, precision=80)
        with localcontext() as contexte:
            contexte.prec = 6
            contexte.rounding = ROUND_HALF_DOWN
            contexte.Emax = 9
            contexte.Emin = -9
            apres = certifier_spectral(d, inv, precision=80)
            self.assertEqual(contexte.prec, 6)
            self.assertEqual(contexte.rounding, ROUND_HALF_DOWN)
            for texte in ("2", "3", "1e-700", "7e700"):
                x = Decimal(texte)
                racine = _Intervalles(16, 100).racine_superieure(x)
                self.assertGreaterEqual(Fraction(racine)**2, Fraction(x))
        for cle in ("lambda_min", "lambda_inferieur_decimal", "eta_superieur_decimal",
                    "trace_superieure_decimal", "operations_decimal"):
            self.assertEqual(avant[cle], apres[cle])

    def test_mauvais_facteur_refuse_et_budgets(self):
        inv = InverseSelectionnee(np.eye(3), np.ones(3))
        with self.assertRaises(CertificationSpectraleImpossible) as cm:
            certifier_spectral(.1*np.eye(3), inv)
        self.assertGreaterEqual(Decimal(cm.exception.bilan["eta_superieur_decimal"]), 1)
        with self.assertRaises(BudgetInverseSelectionnee):
            certifier_spectral(np.eye(3), inv, budget_operations=10)
        with self.assertRaises(BudgetInverseSelectionnee):
            certifier_spectral(np.ones((2, 3)), inv, budget_coefficients=2)

    def test_certificat_reconstruit_les_caches_mutables(self):
        inv = InverseSelectionnee([[1.]], [1.])
        inv.r.data[:] = .5
        cert = certifier_spectral([[.5]], inv)
        self.assertLessEqual(cert["lambda_min"], .25)
        inv.diagonale[:] = 100.
        inv.voisins = [[(0, 123.)]]
        inv.valeurs = [{0: 1e-30}]
        apres = certifier_spectral([[.5]], inv)
        self.assertEqual(cert["lambda_min"], apres["lambda_min"])
        inv.masse.data[:] = 2.
        self.assertLessEqual(certifier_spectral([[.5]], inv)["lambda_min"], .125)
        inv.masse.data[:] = -1.
        with self.assertRaises(ValueError):
            certifier_spectral([[.5]], inv)

    def test_doublons_refuses_avant_toute_somme_flottante(self):
        doublons = csr_matrix(([1., 1., 1.], [0, 0, 1], [0, 2, 3]), shape=(2, 2))
        self.assertFalse(doublons.has_canonical_format)
        for r, m, d in ((doublons, np.eye(2), np.eye(2)),
                        (np.eye(2), doublons, np.eye(2)),
                        (np.eye(2), np.eye(2), doublons)):
            with self.assertRaises(ValueError):
                certifier_spectre(d, r, np.ones(2), m)
        # LIL n'expose pas has_canonical_format ; la conversion CSR doit
        # être contrôlée elle aussi avant de sommer les produits de Gram.
        double_lil = lil_matrix((1, 1))
        double_lil.rows[0] = [0, 0]
        double_lil.data[0] = [1., -.1]
        for r, m, d in (([[np.sqrt(.91)]], [[1.]], double_lil),
                        (double_lil, [[1.]], [[1.]]),
                        ([[1.]], double_lil, [[1.]])):
            with self.assertRaises(ValueError):
                certifier_spectre(d, r, [1.], m)
        double_coo = coo_matrix(([1., -.1], ([0, 0], [0, 0])), shape=(1, 1))
        for r, m, d in (([[np.sqrt(.91)]], [[1.]], double_coo),
                        (double_coo, [[1.]], [[1.]]),
                        ([[1.]], double_coo, [[1.]])):
            with self.assertRaises(ValueError):
                certifier_spectre(d, r, [1.], m)
        # Un utilisateur peut rendre le drapeau périmé par une mutation
        # directe des indices ; seule l'inspection des données est fiable.
        cache_perime = csr_matrix([[1., -.1], [0., 1.]])
        self.assertTrue(cache_perime.has_canonical_format)
        cache_perime.indices[1] = 0
        self.assertTrue(cache_perime.has_canonical_format)
        with self.assertRaises(ValueError):
            certifier_spectre(cache_perime, np.eye(2), [1., 1.], np.eye(2))

    def test_console_native_512_sans_densification_globale(self):
        from vinkulum import Noyau
        if not hasattr(Noyau, "facteurs_materiels_poutres"):
            self.skipTest("export matériel natif absent de ce noyau installé")
        original = csr_matrix.toarray
        def petit_seulement(a, *args, **kwargs):
            if min(a.shape) > 64:
                raise AssertionError("densification d'une grande matrice")
            return original(a, *args, **kwargs)
        with patch.object(csr_matrix, "toarray", petit_seulement), \
             patch("numpy.linalg.inv", side_effect=AssertionError("inverse globale")):
            qr, masse = _console_native(512)
            inverse = InverseSelectionnee(qr.r, masse, qr.echelles)
            cert = certifier_spectral(qr.di, inverse)
        self.assertLess(inverse.coefficients, 6*inverse.n)
        self.assertLess(cert["coefficients_gram_ecart"], 6*inverse.n)
        self.assertLess(cert["eta_superieur"], 1e-3)
        self.assertGreater(cert["lambda_min"], 0)
        self.assertTrue(cert["certification_machine"])
        # Le minorant corrigé ne peut dépasser la trace exacte du facteur.
        borne_qr = Decimal(1)/Decimal(cert["trace_superieure_decimal"])
        self.assertLessEqual(Fraction(cert["lambda_min"]), Fraction(borne_qr))


if __name__ == "__main__":
    unittest.main()
