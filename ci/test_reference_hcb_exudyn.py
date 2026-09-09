"""Réponses réduites : tests algébriques obligatoires, API Exudyn optionnelle."""
import importlib.util
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import eigh, solve
from scipy.sparse import csr_matrix

from reference_hcb_exudyn import _preparer_reponse_reduite, prepare, reponse


def _donnees():
    rng = np.random.default_rng(7402)
    d = rng.normal(size=(39, 18))
    b = rng.normal(size=(31, 18))
    m = b.T@b+np.eye(18)
    return d, d.T@d, m


def _preparation_algebrique():
    d, k, m = _donnees()
    statiques = -solve(k[:-6, :-6], k[:-6, -6:], assume_a="pos")
    _, modes = eigh(k[:-6, :-6], m[:-6, :-6])
    t = np.block([[statiques, modes[:, :4]], [np.eye(6), np.zeros((6, 4))]])
    kr, mr = t.T@k@t, t.T@m@t
    kr, mr = (kr+kr.T)/2, (mr+mr.T)/2
    p = dict(Tphys=t, Kred=kr, Mred=mr, **_preparer_reponse_reduite(kr, mr, t))
    return p, d, k, m


def _dense(p, omega, forces):
    """Oracle de Galerkin par résolution directe, sans caches spectraux."""
    t = p["Tphys"]
    f = np.zeros((t.shape[0],)+forces.shape[1:])
    f[-6:] = forces
    return t@np.linalg.solve(p["Kred"]-omega**2*p["Mred"], t.T@f)


class ReponseAlgebrique(unittest.TestCase):
    def setUp(self):
        self.p, self.d, self.k, self.m = _preparation_algebrique()
        self.w = .4*np.sqrt(self.p["valeurs_reduites"][0])
        self.forces = np.diag([.1, .2, .3, .01, .02, .03])

    def test_masse_couplee_reponse_independante(self):
        self.assertGreater(np.linalg.norm(self.m[:-6, -6:]), 1.)
        x = reponse(self.p, self.w, self.forces, resolution="spectrale_corrigee")
        assert_allclose(x, _dense(self.p, self.w, self.forces), rtol=3e-12, atol=2e-16)
        # Ce témoin doit détecter une projection supprimant M_IS.
        m_faux = self.m.copy()
        m_faux[:-6, -6:] = 0
        m_faux[-6:, :-6] = 0
        faux = dict(self.p, Mred=self.p["Tphys"].T@m_faux@self.p["Tphys"])
        self.assertGreater(np.linalg.norm(x-_dense(faux, self.w, self.forces))/np.linalg.norm(x), 1e-3)

    def test_vecteur_et_six_colonnes(self):
        x = reponse(self.p, self.w, self.forces, resolution="spectrale_corrigee")
        for col in range(6):
            f = self.forces[:, col]
            v = reponse(self.p, self.w, f, resolution="spectrale_corrigee")
            self.assertEqual(v.shape, (18,))
            assert_allclose(v, x[:, col], rtol=2e-12, atol=2e-16)
            assert_allclose(v, _dense(self.p, self.w, f), rtol=3e-12, atol=2e-16)

    def test_melanges_de_charges_et_pas_mutation(self):
        f = np.random.default_rng(841).normal(size=(6, 9))
        sauvegardes = {k: v.copy() for k, v in self.p.items() if isinstance(v, np.ndarray)}
        f_avant = f.copy()
        x = reponse(self.p, self.w, f, resolution="spectrale_corrigee")
        assert_allclose(x, _dense(self.p, self.w, f), rtol=3e-12, atol=5e-15)
        np.testing.assert_array_equal(f, f_avant)
        for nom, avant in sauvegardes.items():
            np.testing.assert_array_equal(self.p[nom], avant, err_msg=nom)

    def test_proximite_resonance_controlee(self):
        # Ecart relatif 10^-4 en omega² : amplification sensible, sans pôle.
        w = np.sqrt(self.p["valeurs_reduites"][0]*(1-1e-4))
        ref = _dense(self.p, w, self.forces)
        x = reponse(self.p, w, self.forces, resolution="spectrale_corrigee")
        self.assertLess(np.linalg.norm(x-ref)/np.linalg.norm(ref), 2e-10)

    def test_correction_utilise_paire_initiale_pas_spectre_arrondi(self):
        # Simule une erreur du préconditionneur modal. Un résidu calculé
        # dans la paire reconstruite depuis ce spectre erroné serait nul.
        p = dict(self.p, valeurs_reduites=self.p["valeurs_reduites"]*(1+2e-4))
        ref = _dense(p, self.w, self.forces)
        brut = reponse(p, self.w, self.forces, resolution="spectrale")
        corrige = reponse(p, self.w, self.forces, resolution="spectrale_corrigee")
        eb = np.linalg.norm(brut-ref)/np.linalg.norm(ref)
        ec = np.linalg.norm(corrige-ref)/np.linalg.norm(ref)
        self.assertGreater(eb, 1e-5)
        self.assertLess(ec, 1e-9)
        self.assertLess(ec, eb*1e-6)

    def test_changement_echelles_base_preserve_reponse(self):
        facteurs = np.logspace(-2, 2, self.p["Kred"].shape[0])
        t = self.p["Tphys"]*facteurs[None, :]
        kr = facteurs[:, None]*self.p["Kred"]*facteurs[None, :]
        mr = facteurs[:, None]*self.p["Mred"]*facteurs[None, :]
        p = dict(Tphys=t, Kred=kr, Mred=mr, **_preparer_reponse_reduite(kr, mr, t))
        assert_allclose(reponse(p, self.w, self.forces, resolution="spectrale_corrigee"),
                        _dense(self.p, self.w, self.forces), rtol=3e-12, atol=2e-16)

    def test_spectrale_et_directe_restent_auditables(self):
        for resolution in ("spectrale", "directe", "spectrale_corrigee"):
            with self.subTest(resolution=resolution):
                assert_allclose(reponse(self.p, 0., self.forces, resolution=resolution),
                                _dense(self.p, 0., self.forces), rtol=3e-12, atol=2e-16)

    def test_pole_exact_refuse_meme_charge_nulle(self):
        valeurs = self.p["valeurs_reduites"].copy()
        valeurs[0] = 4.
        p = dict(self.p, valeurs_reduites=valeurs)
        for resolution in ("spectrale", "spectrale_corrigee"):
            with self.assertRaisesRegex(ArithmeticError, "dénominateur modal nul"):
                reponse(p, 2., np.zeros(6), resolution=resolution)

    def test_entrees_non_reelles_ou_non_finies_refusees(self):
        for w, f in ((self.w+0j, self.forces),
                     (self.w, self.forces.astype(complex)),
                     (float("nan"), self.forces),
                     (self.w, np.full(6, np.inf)),
                     (self.w, np.empty((6, 0)))):
            with self.subTest(w=w, forme=f.shape):
                with self.assertRaises(ValueError):
                    reponse(self.p, w, f, resolution="spectrale_corrigee")


@unittest.skipUnless(importlib.util.find_spec("exudyn") is not None,
                     "intégration officielle optionnelle : Exudyn absent")
class IntegrationExudyn(unittest.TestCase):
    def test_audit_optionnel_preserve_base_et_reponse(self):
        d, _, m = _donnees()
        preparations = []
        for complet in (True, False):
            np.random.seed(42123)
            preparations.append(prepare(csr_matrix(d), csr_matrix(m), 4,
                projection="energie", diagnostics_complets=complet))
        a, b = preparations
        self.assertTrue(a["diagnostics"]["complet"])
        self.assertFalse(b["diagnostics"]["complet"])
        for nom in ("residu_statique_equilibre_relatif", "residu_modal_equilibre_relatif",
                    "ecart_projection_K_vers_D_relatif", "ecart_projection_statique_relatif",
                    "asymetrie_projection_K", "asymetrie_projection_M"):
            self.assertIsNone(b["diagnostics"][nom])
            self.assertIsNotNone(a["diagnostics"][nom])
        # Le solveur propre peut retourner des signes différents : seule
        # l'identité de l'espace et de la réponse physique est pertinente.
        qa, _ = np.linalg.qr(a["Tphys"])
        qb, _ = np.linalg.qr(b["Tphys"])
        assert_allclose(qa@qa.T, qb@qb.T, rtol=3e-12, atol=2e-13)
        assert_allclose(reponse(a, .1, np.eye(6), resolution="spectrale_corrigee"),
                        reponse(b, .1, np.eye(6), resolution="spectrale_corrigee"),
                        rtol=3e-11, atol=2e-14)

    def test_paire_officielle_masse_couplee_et_deux_projections(self):
        d, _, m = _donnees()
        self.assertGreater(np.linalg.norm(m[:-6, -6:]), 1.)
        f = np.random.default_rng(68).normal(size=(6, 7))
        for projection in ("standard", "energie"):
            with self.subTest(projection=projection):
                np.random.seed(42123)
                p = prepare(csr_matrix(d), csr_matrix(m), 4, projection=projection)
                w = .3*np.sqrt(p["valeurs_reduites"][0])
                assert_allclose(p["Mred"], p["Tphys"].T@m@p["Tphys"], rtol=3e-12, atol=2e-14)
                assert_allclose(reponse(p, w, f, resolution="spectrale_corrigee"),
                                _dense(p, w, f), rtol=4e-12, atol=5e-15)
                self.assertEqual(p["diagnostics"]["dimensions"], [18, 10])
                self.assertEqual(p["diagnostics"]["reponse_spectrale_corrigee"]["nombre_corrections"], 2)

    def test_geometrie_auxiliaire_ne_change_pas_probleme_algebrique(self):
        d, k, m = _donnees()
        geometries = (np.zeros((6, 3)), np.random.default_rng(18).normal(size=(6, 3)))
        reponses = []
        for positions in geometries:
            np.random.seed(42123)
            p = prepare(csr_matrix(d), csr_matrix(m), 4,
                        projection="energie", positions=positions)
            reponses.append(reponse(p, .1, np.eye(6), resolution="spectrale_corrigee"))
            # Trace physique S_boundary [I,0], sans liaison rigide RBE.
            attendu = np.zeros((6, 10))
            attendu[:, :6] = np.diag(1/np.sqrt(np.diag(k)[-6:]))
            assert_allclose(p["Tphys"][-6:], attendu, rtol=1e-13, atol=1e-15)
        assert_allclose(reponses[0], reponses[1], rtol=3e-11, atol=2e-14)


if __name__ == "__main__":
    unittest.main()
