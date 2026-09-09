"""Contre-épreuves physiques des fixtures de ports et de la référence HCB."""
from dataclasses import replace
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import eigh, solve
from scipy.sparse import csc_matrix

from modeles_ports import (borne_lambda_trace_console, chaine, console,
                           flexibilites_conditionnees_console, metrique_ports,
                           reference_chaine)
import reference_ports_hcb as hcb


def _couple_spd():
    """Intérieur couplé et masse non diagonale, sans couplage massique au port."""
    rng = np.random.default_rng(8742)
    ni, p = 9, 3
    a, b = rng.normal(size=(ni, ni)), rng.normal(size=(ni, p))
    c = rng.normal(size=(ni, ni))
    kii, mii = a.T@a+np.eye(ni), c.T@c+2*np.eye(ni)
    metrique = np.array([[3., .2, -.1], [.2, 1.7, .3], [-.1, .3, 2.1]])
    kss = metrique+b.T@solve(kii, b, assume_a="pos")
    k = np.block([[kii, b], [b.T, kss]])
    m = np.block([[mii, np.zeros((ni, p))],
                  [np.zeros((p, ni)), np.diag([.3, .8, 1.4])]])
    # Borne inférieure en arithmétique exacte, sans prendre un Ritz supérieur.
    lambda_min = 1/np.trace(solve(kii, mii, assume_a="pos"))
    return k, m, np.arange(ni), np.arange(ni, ni+p), metrique, lambda_min


def _schur_dense(k, m, ii, ss, omega, w):
    a = k-omega**2*m
    return w.T@(a[np.ix_(ss, ss)]-a[np.ix_(ss, ii)]
                @solve(a[np.ix_(ii, ii)], a[np.ix_(ii, ss)], assume_a="sym"))@w


class ReferencePortsHCB(unittest.TestCase):
    def test_chaine_statique_et_tous_modes(self):
        n, ressort, masse = 7, 23., .017
        model = chaine(n, ressort, masse)
        metric = metrique_ports(model)
        guyan = hcb.reduire(model.k, model.m, model.interieur, model.interface, 0, metric)
        np.testing.assert_allclose(guyan.schur(0), [[1.]], rtol=2e-14)
        # Somme exacte des masses physiques le long de la déformée i/n.
        omega = .31*np.sqrt(model.metadata["mu_min_interieur_s_moins_2"])
        attendu = 1-omega**2*masse*(n+1)*(2*n+1)/(6*ressort)
        np.testing.assert_allclose(guyan.schur(omega), [[attendu]], rtol=2e-14)
        complet = hcb.reduire(model.k, model.m, model.interieur, model.interface, n-1, metric)
        for rapport in (0., .2, .7):
            with self.subTest(rapport=rapport):
                omega = rapport*np.sqrt(model.metadata["mu_min_interieur_s_moins_2"])
                oracle = reference_chaine(n, omega, ressort, masse)
                np.testing.assert_allclose(complet.schur(omega), [[oracle["schur"]/metric[0, 0]]],
                                           rtol=3e-13, atol=2e-14)
                port = .61
                champ = complet.reconstruit(omega, [port])
                np.testing.assert_allclose(champ, oracle["rapports_deplacement"]*complet.w[0, 0]*port,
                                           rtol=3e-13, atol=2e-14)

    def test_champ_et_schur_conjugues(self):
        k, m, ii, ss, metric, lambda_min = _couple_spd()
        reduction = hcb.reduire(k, m, ii, ss, 3, metric)
        omega = np.sqrt(.4*lambda_min)
        a = k-omega**2*m
        u = reduction.reconstruit(omega, np.eye(len(ss)))
        s = reduction.schur(omega)
        np.testing.assert_allclose(u[ss], reduction.w, rtol=1e-14, atol=1e-14)
        # Travail conjugué au port, avec énergie du champ reconstruit.
        np.testing.assert_allclose(u.T@a@u, s, rtol=2e-12, atol=2e-13)
        residu = a[np.ix_(ii, ii)]@u[ii]+a[np.ix_(ii, ss)]@u[ss]
        np.testing.assert_allclose(reduction.phi.T@residu, 0., atol=3e-13)
        delta = s-_schur_dense(k, m, ii, ss, omega, reduction.w)
        energie = residu.T@solve(a[np.ix_(ii, ii)], residu, assume_a="pos")
        np.testing.assert_allclose(delta, energie, rtol=3e-9, atol=5e-13)
        self.assertEqual(reduction.compteurs["factorisations_M_II_controle"], 1)

    def test_borne_residuelle_et_gram_reel(self):
        k, m, ii, ss, metric, lambda_min = _couple_spd()
        for r in (0, 3):
            reduction = hcb.reduire(k, m, ii, ss, r, metric)
            omega_max = np.sqrt(.7*lambda_min)
            if r:
                triangulaire = hcb.solve_triangular

                def coordonnees_non_orthonormales(a, b, **kwargs):
                    # Changement inversible de coordonnées dans le même espace.
                    # Il force le calcul de la borne à traiter le Gram réel.
                    return np.array([.7, 1.3, 2.1])[:, None]*triangulaire(a, b, **kwargs)

                with patch.object(hcb, "solve_triangular", side_effect=coordonnees_non_orthonormales):
                    borne = reduction.borne_uniforme(lambda_min, omega_max)
                self.assertGreater(reduction.bilan_borne["defaut_gram"], 1.)
            else:
                borne = reduction.borne_uniforme(lambda_min, omega_max)
            for omega in np.linspace(0., omega_max, 9):
                with self.subTest(r=r, omega=omega):
                    a = k-omega**2*m
                    u = reduction.reconstruit(omega, np.eye(len(ss)))
                    residu = a[np.ix_(ii, ii)]@u[ii]+a[np.ix_(ii, ss)]@u[ss]
                    gram = residu.T@solve(a[np.ix_(ii, ii)], residu, assume_a="pos")
                    erreur = eigh((gram+gram.T)*.5, eigvals_only=True)[-1]
                    self.assertLessEqual(erreur, borne*(1+1e-10)+1e-25)
            self.assertEqual(reduction.compteurs["resolutions_borne_K_II"], 1)
            self.assertEqual(reduction.bilan_borne["statut"], "majorant_conditionnel_flottant_non_certifie")

    def test_modes_creux_et_compteurs(self):
        n, r = 80, 4
        model = chaine(n, 23., .017)
        reduction = hcb.reduire(model.k, model.m, model.interieur, model.interface, r,
                                metrique_ports(model))
        j = np.arange(1, r+1)
        exactes = 4*23/.017*np.sin(j*np.pi/(2*n))**2
        np.testing.assert_allclose(reduction.valeurs_modales, exactes, rtol=2e-11)
        c = reduction.compteurs
        self.assertEqual(c["methode_modes"], "eigsh_shift_inverse_zero")
        self.assertEqual(c["eigensolves_denses_interieurs"], 0)
        self.assertEqual(c["taille_dense_interieure"], 0)
        self.assertEqual(c["factorisations_K_II"], 1)
        self.assertEqual(c["resolutions_statiques_K_II"], 1)
        self.assertGreater(c["resolutions_modes_K_II"], 0)
        reduction.borne_uniforme(exactes[0], np.sqrt(.5*exactes[0]))
        self.assertEqual(c["resolutions_K_II"], c["resolutions_modes_K_II"]+2)
        self.assertEqual(c["colonnes_rhs_K_II"], c["resolutions_modes_K_II"]+1+(2+r))
        solves = c["resolutions_K_II"]
        for omega in (0., .1, .2):
            reduction.schur(omega)
        self.assertEqual(c["resolutions_K_II"], solves)
        self.assertEqual(c["factorisations_frequence"], 0)
        self.assertEqual(c["appels_schur"], 3)
        self.assertEqual(c["divisions_modales_rhs"], 3*r)

    def test_flexibilite_console_et_trace(self):
        for n in (2, 3, 8):
            with self.subTest(n=n):
                model = console(n)
                ii = model.interieur
                kii = model.k[ii, :][:, ii].toarray()
                mii = model.m[ii, :][:, ii].toarray()
                inverse = solve(kii, np.eye(len(ii)), assume_a="pos")
                analytique = flexibilites_conditionnees_console(model.metadata).ravel()
                np.testing.assert_allclose(analytique, np.diag(inverse), rtol=2e-11)
                borne = borne_lambda_trace_console(model.metadata)
                np.testing.assert_allclose(borne, 1/np.trace(mii@inverse), rtol=2e-11)
                minimum = eigh(kii, mii, eigvals_only=True)[0]
                self.assertGreater(borne, 0.)
                self.assertLessEqual(borne, minimum*(1+1e-11))
                self.assertEqual(model.metadata["formulation_poutre"], "integree")
                self.assertAlmostEqual(model.m.diagonal()[::6].sum(), model.metadata["masse_mobile_kg"])

    def test_metrique_ports_et_permutation(self):
        model = console(5)
        metric = metrique_ports(model)
        reference = _schur_dense(model.k.toarray(), model.m.toarray(), model.interieur,
                                 model.interface, 0., np.eye(6))
        d = np.sqrt(np.diag(metric))
        np.testing.assert_allclose(reference/d[:, None]/d[None, :],
                                   metric/d[:, None]/d[None, :], rtol=3e-12, atol=3e-12)
        self.assertTrue(np.all(eigh(metric, eigvals_only=True) > 0))
        ordre = np.array([2, 0, 5, 1, 4, 3])
        permute = replace(model, interface=model.interface[ordre])
        np.testing.assert_array_equal(metrique_ports(permute), metric[np.ix_(ordre, ordre)])
        np.testing.assert_array_equal(metrique_ports(chaine(500, 23., .017)), [[23/500]])

    def test_oracle_chaine_bande_et_poles(self):
        for n in (2, 7, 19):
            model = chaine(n, 23., .017)
            k, m = model.k.toarray(), model.m.toarray()
            for rho in (0., .1, 1.9, 2., 2.1, 10.):
                with self.subTest(n=n, rho=rho):
                    omega = rho*np.sqrt(23/.017)
                    oracle = reference_chaine(n, omega, 23., .017, force=1.7)
                    f = np.zeros(n)
                    f[-1] = 1.7
                    champ = solve(k-omega**2*m, f, assume_a="sym")
                    np.testing.assert_allclose(oracle["deplacement"], champ, rtol=2e-11, atol=2e-13)
                    direct = _schur_dense(k, m, model.interieur, model.interface, omega, np.eye(1))[0, 0]
                    self.assertLess(abs(oracle["schur"]-direct), 2e-11*max(abs(direct), 23/n))
        for pole in (2*np.sin(np.pi/14), 2*np.sin(np.pi/30)):
            with self.assertRaisesRegex(ValueError, "pôle|résonance"):
                reference_chaine(7, pole)

    def test_rejets(self):
        for k, m, ii, ss in ((np.diag([-1., 1., 1.]), np.eye(3), [0, 1], [2]),
                             (np.eye(3), np.ones((3, 3))+np.eye(3), [0, 1], [2]),
                             (np.eye(3), np.eye(3), [0, 1], [1]),
                             (np.eye(3), np.eye(3), [0., 1.], [2])):
            with self.subTest(ii=ii, ss=ss), self.assertRaises(ValueError):
                hcb.reduire(csc_matrix(k), csc_matrix(m), ii, ss, 1, np.eye(len(ss)))
        model = chaine(8)
        for r in (-1, 8, 2.5, True):
            with self.subTest(r=r), self.assertRaises(ValueError):
                hcb.reduire(model.k, model.m, model.interieur, model.interface, r, metrique_ports(model))
        reduction = hcb.reduire(model.k, model.m, model.interieur, model.interface, 2,
                                metrique_ports(model))
        with self.assertRaisesRegex(ValueError, "quotient de Rayleigh"):
            reduction.borne_uniforme(2*model.metadata["mu_min_interieur_s_moins_2"], 0.)
        with self.assertRaisesRegex(ValueError, "omega_max"):
            reduction.borne_uniforme(1., 1.)


if __name__ == "__main__":
    unittest.main()
