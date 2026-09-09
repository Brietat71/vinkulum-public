"""Tests de champ de l'oracle Decimal, indépendants de la réduction."""
from decimal import Decimal, ROUND_DOWN, localcontext
from fractions import Fraction
import unittest
from unittest.mock import patch

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import spsolve

from oracle_champs_ports import OracleChamps


def _resout_fraction(d, masses, omega, forces, p):
    d = [[Fraction(float(x)) for x in row] for row in d]
    masses = [Fraction(float(x)) for x in masses]
    z = Fraction(float(omega))**2
    n, q = len(masses), forces.shape[1]
    a = [[sum(row[i]*row[j] for row in d)-(z*masses[i] if i == j else 0)
          for j in range(n)] for i in range(n)]
    rhs = [[Fraction(0) for _ in range(q)] for _ in range(n-p)]
    rhs += [[Fraction(float(x)) for x in row] for row in forces]
    aug = [a[i]+rhs[i] for i in range(n)]
    for i in range(n):
        pivot = next(j for j in range(i, n) if aug[j][i])
        aug[i], aug[pivot] = aug[pivot], aug[i]
        valeur = aug[i][i]
        aug[i] = [x/valeur for x in aug[i]]
        for j in range(n):
            if j != i:
                valeur = aug[j][i]
                aug[j] = [a-valeur*b for a, b in zip(aug[j], aug[i])]
    return [row[n:] for row in aug]


def _facteur_blocs(n, p, seed):
    rng = np.random.default_rng(seed)
    d = np.zeros((n*p, n*p))
    for j in range(n):
        bloc = np.triu(rng.integers(-3, 4, size=(p, p))/8.)+2*np.eye(p)
        d[j*p:(j+1)*p, j*p:(j+1)*p] = bloc
        if j+1 < n:
            d[j*p:(j+1)*p, (j+1)*p:(j+2)*p] = rng.integers(-3, 4, size=(p, p))/8.
    return d


class TestOracleChamps(unittest.TestCase):
    def test_champs_fraction_et_nombre_arbitraire_de_forces(self):
        p, n, omega = 2, 3, .375
        d = _facteur_blocs(n, p, 88)
        masses = np.linspace(1., 2., n*p)
        oracle = OracleChamps(csr_matrix(d), masses, p=p, dps=70)
        for q in (1, 3, 6):
            forces = np.arange(1, p*q+1).reshape(p, q)/4.
            reference = _resout_fraction(d, masses, omega, forces, p)
            calcule = oracle.reponse(omega, forces, retour_decimal=True)
            for i in range(n*p):
                for j in range(q):
                    self.assertLess(abs(Fraction(calcule[i][j])-reference[i][j]), Fraction(1, 10**60))
        vecteur = oracle.reponse(omega, forces[:, 0])
        np.testing.assert_array_equal(vecteur, oracle.reponse(omega, forces)[:, 0])
        self.assertEqual(oracle.compteurs["assemblages_raideur"], 1)
        self.assertEqual(oracle.compteurs["condensations_frequence"], 1)
        self.assertEqual(oracle.compteurs["resolutions_blocs"], n)
        self.assertFalse(oracle.certification_machine)

    def test_poutre_native_statique_timoshenko(self):
        from vinkulum import Noyau
        if not hasattr(Noyau, "facteurs_materiels_poutres"):
            self.skipTest("export matériel absent du noyau installé")
        from experience_reduction_native import modele, parametres, reference_statique
        p = parametres(8)
        noyau, libres, _, _ = modele(p)
        donnees = noyau.facteurs_materiels_poutres()
        def lire(t):
            nr, nc, ptr, ids, valeurs = t
            return csc_matrix((valeurs, ids, ptr), shape=(nr, nc))
        d = lire(donnees["d"])[:, libres]
        m = lire(donnees["masse"])[libres][:, libres]
        oracle = OracleChamps(d, m, p=6)
        champ = oracle.reponse(0., p["force_port_physique"])
        np.testing.assert_allclose(champ, reference_statique(p), rtol=2e-12, atol=1e-16)
        charges = np.diag([.1, .1, .1, .01, .01, .01])
        six = oracle.reponse(0., charges)
        self.assertEqual(six.shape, (48, 6))
        self.assertAlmostEqual(six[-6, 0]/(.1*p["longueur_m"]/p["ea_n"]), 1., places=12)
        self.assertAlmostEqual(six[-3, 3]/(.01*p["longueur_m"]/p["gj_nm2"]), 1., places=12)

    def test_deux_precisions_scipy_et_contexte_prive(self):
        p, n = 3, 4
        dense = _facteur_blocs(n, p, 458)
        d = csc_matrix(dense)
        masses = np.linspace(.7, 1.2, n*p)
        m = diags(masses, format="csr")
        o70, o90 = OracleChamps(d, m, p=p, dps=70), OracleChamps(d, m, p=p, dps=90)
        charges = np.arange(1, 16).reshape(3, 5)/10.
        for omega in (0., .2, .6):
            rhs = np.zeros((n*p, charges.shape[1]))
            rhs[-p:] = charges
            scipy = spsolve((d.T@d-omega**2*m).tocsc(), rhs)
            a, b = o70.reponse(omega, charges), o90.reponse(omega, charges)
            np.testing.assert_array_equal(a, b)
            np.testing.assert_allclose(a, scipy, rtol=2e-13, atol=2e-15)
            with localcontext() as contexte:
                contexte.prec = 5
                contexte.rounding = ROUND_DOWN
                np.testing.assert_array_equal(OracleChamps(d, m, p=p).reponse(omega, charges), a)
                self.assertEqual(contexte.prec, 5)
        self.assertEqual(o70.compteurs["assemblages_raideur"], 1)
        self.assertEqual(o70.compteurs["condensations_frequence"], 3)

    def test_refus_non_voisins_et_entrees_hors_domaine(self):
        with self.assertRaisesRegex(ValueError, "non voisins"):
            OracleChamps(csr_matrix([[1., 0., 1.]]), np.ones(3), p=1)
        with self.assertRaisesRegex(ValueError, "diagonale"):
            OracleChamps(csr_matrix(np.eye(2)), csr_matrix([[1., .1], [.1, 1.]]), p=1)
        doublons = csr_matrix(([1., -.1, 1.], [0, 0, 1], [0, 2, 3]), shape=(2, 2))
        doublons.has_canonical_format = True
        with self.assertRaises(ValueError):
            OracleChamps(doublons, np.ones(2), p=1)
        with self.assertRaises(ValueError):
            OracleChamps(lil_matrix(np.eye(2)), np.ones(2), p=1)
        oracle = OracleChamps(csr_matrix(np.eye(2)), np.ones(2), p=1)
        for omega in (-1., np.inf, np.nan, 1.+1.j):
            with self.assertRaises((ValueError, TypeError)):
                oracle.reponse(omega, [1.])
        for force in ([np.inf], np.ones((2, 1)), np.ones((1, 0)), np.array([1.+1.j])):
            with self.assertRaises(ValueError):
                oracle.reponse(0., force)

    def test_poles_et_pivot_interieur_singulier_refuses(self):
        oracle = OracleChamps(csr_matrix(np.eye(2)), np.ones(2), p=1)
        with self.assertRaisesRegex(ValueError, "pivot intérieur"):
            oracle.reponse(1., [0.])
        terminal = OracleChamps(csr_matrix([[1.]]), np.ones(1), p=1)
        with self.assertRaisesRegex(ValueError, "Schur terminal singulier"):
            terminal.reponse(1., [0.])
        # A=[[0,1],[1,1]] est inversible : le refus vient bien du schéma
        # de blocs sans permutation globale, pas d'une preuve de singularité.
        bloc = OracleChamps(csr_matrix([[1., 1.], [0., 1.]]), np.ones(2), p=1)
        with self.assertRaisesRegex(ValueError, "pivot intérieur"):
            bloc.reponse(1., [1.])
        self.assertTrue(np.all(np.isfinite(terminal.reponse(1.-1e-12, [1.]))))

    def test_chaine1024_sans_matrice_globale_dense(self):
        n = 1024
        d = diags((-np.ones(n-1), np.ones(n)), (-1, 0), format="csr")
        with patch.object(csr_matrix, "toarray", side_effect=AssertionError("densification globale")), \
             patch("numpy.linalg.inv", side_effect=AssertionError("inverse globale")):
            oracle = OracleChamps(d, np.ones(n), p=1)
            champs = oracle.reponse(0., np.array([[1., 2., -.5]]))
        attendu = np.arange(1., n+1)[:, None]*np.array([[1., 2., -.5]])
        np.testing.assert_allclose(champs, attendu, rtol=1e-14, atol=0.)
        self.assertEqual(oracle.compteurs["assemblages_raideur"], 1)
        self.assertEqual(oracle.compteurs["condensations_frequence"], 1)
        self.assertEqual(oracle.compteurs["resolutions_blocs"], n)


if __name__ == "__main__":
    unittest.main()
