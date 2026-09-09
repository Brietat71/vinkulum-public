"""Déformations élémentaires et relèvements analytiques contre le noyau figé."""
import unittest

import numpy as np

from energie_ports import (energie_chaine, energie_console, lifting_chaine,
                           lifting_console)
from modeles_ports import chaine, console, metrique_ports


def _ecart_raideur(a, b):
    difference = (a-b).tocoo()
    echelle = np.sqrt(b.diagonal())
    return float(np.max(np.abs(difference.data)/echelle[difference.row]
                        /echelle[difference.col], initial=0.))


class EnergieElementairePorts(unittest.TestCase):
    def test_energie_chaine_et_masses_physiques(self):
        for n in (2, 8, 32):
            with self.subTest(n=n):
                model = chaine(n, 23., .017)
                energie = energie_chaine(n, 23., .017)
                self.assertEqual(energie.d.shape, (n, n))
                self.assertEqual(energie.d.nnz, 2*n-1)
                self.assertLess(_ecart_raideur(energie.d.T@energie.d, model.k), 5e-15)
                np.testing.assert_array_equal(energie.m.diagonal(), model.m.diagonal())
                np.testing.assert_array_equal(energie.interieur, model.interieur)
                np.testing.assert_array_equal(energie.interface, model.interface)
                u = np.random.default_rng(n).normal(size=n)
                differences = np.diff(np.r_[0., u])
                np.testing.assert_allclose(np.linalg.norm(energie.d@u)**2,
                                           23*np.dot(differences, differences), rtol=2e-15)

    def test_energie_console_native_integree(self):
        for n in (2, 8, 32):
            with self.subTest(n=n):
                model = console(n)
                energie = energie_console(model.metadata)
                self.assertEqual(energie.d.shape, (6*n, 6*n))
                self.assertEqual(energie.d.nnz, 16*n-8)
                self.assertLess(_ecart_raideur(energie.d.T@energie.d, model.k), 5e-14)
                np.testing.assert_array_equal(energie.m.diagonal(), model.m.diagonal())
                np.testing.assert_array_equal(energie.interieur, model.interieur)
                np.testing.assert_array_equal(energie.interface, model.interface)
                self.assertNotIn("representation_energie", model.metadata)
                # Les deux plans sont couplés aux rotations avec des signes opposés.
                premier = energie.d[:6, :6].toarray()
                self.assertLess(premier[1, 5], 0.)
                self.assertGreater(premier[2, 4], 0.)

    def test_lifting_chaine_equilibre_et_energie_du_port(self):
        for n in (2, 8, 32):
            with self.subTest(n=n):
                model = chaine(n, 23., .017)
                energie = energie_chaine(n, 23., .017)
                psi = lifting_chaine(n)
                self.assertEqual(psi.shape, (n-1, 1))
                u = np.vstack((psi, [[1.]]))
                residu = model.k@u
                np.testing.assert_allclose(residu[:-1], 0., atol=2e-14)
                np.testing.assert_allclose(residu[-1], [23/n], rtol=2e-14)
                d = energie.d@u
                np.testing.assert_allclose(d.T@d, metrique_ports(model), rtol=2e-14)
                # Valeurs nodales et incréments connus, sans résolution globale.
                np.testing.assert_allclose(np.diff(np.r_[0., u[:, 0]]), 1/n, rtol=2e-14)

    def test_lifting_console_equilibre_et_energie_du_port(self):
        for n in (2, 8, 32):
            with self.subTest(n=n):
                model = console(n)
                energie = energie_console(model.metadata)
                psi = lifting_console(model.metadata)
                self.assertEqual(psi.shape, (6*(n-1), 6))
                u = np.vstack((psi, np.eye(6)))
                metric = metrique_ports(model)
                normalisation = np.sqrt(np.diag(metric))
                residu = (model.k@u)[model.interieur]
                rows = np.sqrt(model.k.diagonal()[model.interieur])
                erreur = np.max(np.abs(residu)/rows[:, None]/normalisation[None, :])
                self.assertLess(erreur, 3e-12)
                du = energie.d@u
                np.testing.assert_allclose((du.T@du)/normalisation[:, None]/normalisation[None, :],
                                           metric/normalisation[:, None]/normalisation[None, :],
                                           rtol=2e-13, atol=2e-13)
                # La raideur condensée est évaluée comme une somme d'énergies.
                for j in range(6):
                    np.testing.assert_allclose(np.dot(du[:, j], du[:, j]), metric[j, j], rtol=2e-13)

    def test_polynomes_contre_flexibilite_root_free(self):
        for longueur, rayon in ((.4, .04), (2., .0001)):
            with self.subTest(longueur=longueur, rayon=rayon):
                n = 8
                model = console(n, longueur=longueur, rayon=rayon)
                metadata = model.metadata
                ga, ei = metadata["GA_N"], metadata["EI_N_m2"]
                fll = np.array([[longueur/ga+longueur**3/(3*ei), longueur**2/(2*ei)],
                                [longueur**2/(2*ei), longueur/ei]])
                psi = lifting_console(metadata).reshape(n-1, 6, 6)
                for j in range(1, n):
                    x = longueur*j/n
                    fxl = np.array([[x/ga+x*x*(3*longueur-x)/(6*ei), x*x/(2*ei)],
                                    [(longueur*x-x*x/2)/ei, x/ei]])
                    # Résolution locale 2x2, indépendante des polynômes factorisés.
                    attendu = np.linalg.solve(fll.T, fxl.T).T
                    np.testing.assert_allclose(psi[j-1][np.ix_([1, 5], [1, 5])], attendu,
                                               rtol=3e-13, atol=3e-14)
                    signes = np.array([1., -1.])
                    np.testing.assert_allclose(psi[j-1][np.ix_([2, 4], [2, 4])],
                                               signes[:, None]*attendu*signes[None, :],
                                               rtol=3e-13, atol=3e-14)
                energie = energie_console(metadata)
                self.assertLess(_ecart_raideur(energie.d.T@energie.d, model.k), 5e-14)

    def test_rejets_hors_fixture(self):
        for n in (1, 2.5, True):
            with self.subTest(n=n), self.assertRaises(ValueError):
                energie_chaine(n)
            with self.subTest(n=n), self.assertRaises(ValueError):
                lifting_chaine(n)
        for raideur, masse in ((0., 1.), (1., -1.), (np.inf, 1.)):
            with self.subTest(raideur=raideur, masse=masse), self.assertRaises(ValueError):
                energie_chaine(2, raideur, masse)
        metadata = console(2).metadata
        for changement in ({"formulation_poutre": "milieu"}, {"poids_extremites_masse_et_inertie": 1.},
                            {"ports": [0, 1, 2, 3, 4, 5]}, {"EA_N": -1.}):
            invalide = dict(metadata, **changement)
            with self.subTest(changement=changement), self.assertRaises(ValueError):
                energie_console(invalide)
            with self.subTest(changement=changement), self.assertRaises(ValueError):
                lifting_console(invalide)


if __name__ == "__main__":
    unittest.main()
