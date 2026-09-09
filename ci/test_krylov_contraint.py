"""Réponses physiques et majorants du Krylov contraint face à Fraction.

Les oracles inversent le modèle D.T D-z M exact de petite taille, sans
réutiliser de base projetée, de KKT ou de factorisation numérique candidate.
"""
from fractions import Fraction as F
import unittest

import numpy as np
from scipy.sparse import csr_matrix

from vinkulum._ports.condensation_energie import CondensationEnergie
from krylov_contraint import KrylovContraint
from test_trace_complement import gram, inverse, produit, difference, echelle


def rationnel(a):
    return [[F(float(v)) for v in row] for row in a]


def flottant(a):
    return np.array([[float(v) for v in row] for row in a])


def oracle(d, m, ports, omega, forces):
    raideur = gram(rationnel(d))
    dynamique = difference(raideur, echelle(rationnel(m), F(float(omega))**2))
    charges = [[F(0) for _ in range(forces.shape[1])] for _ in range(d.shape[1])]
    for i, coordonnee in enumerate(ports):
        charges[coordonnee] = [F(float(v)) for v in forces[i]]
    return produit(inverse(dynamique), charges)


def energie(v, a):
    return sum((v[i]*a[i][j]*v[j] for i in range(len(v)) for j in range(len(v))), F(0))


def donnees(ni=4):
    n = ni+2
    d = np.diag(np.arange(4., n+4))
    for i in range(n-1):
        d[i, i+1] = .25
    for i in range(n-2):
        d[i, i+2] = -.125
    # Blocs massiques couplés, dont le dernier joint intérieur et ports.
    # Chaque bloc a au plus quatre DDL, compatible avec la racine locale.
    groupes = [list(range(i, min(i+3, n))) for i in range(0, n, 3)]
    if len(groupes[-1]) == 1 and len(groupes) > 1:
        groupes[-2].extend(groupes.pop())
    l = np.eye(n)
    for groupe in groupes:
        for i, j in zip(groupe[:-1], groupe[1:]):
            l[i, j] = .125 if i % 2 else -.25
    m = l.T@l
    if ni == 4:
        b = np.array([[1.], [2.], [-1.], [1.]])
        phi = np.array([[-1.], [1.], [2.], [3.]])
    else:
        b = np.zeros((ni, 1)); b[0, 0] = 1.
        phi = b.copy()
    i, s = np.arange(ni), np.arange(ni, n)
    qr = CondensationEnergie(csr_matrix(d), i, s, np.diag([4., .25]))
    # ||D-diag(D)||_2<=3/8 et ||L||_2<=5/4 : lambda>=4 est conservateur
    # même sur l'espace complet ; aucune valeur propre numérique présumée.
    return d, m, i, s, b, phi, qr


class KrylovContraintPhysique(unittest.TestCase):
    def forces(self):
        epsilon = 2.**-30
        return np.array([[1., 0., 1., epsilon, 0., 1.],
                         [0., 1., 1., -epsilon, 0., -1.]])

    def comparer(self, resultat, exact, tolerance=3e-11):
        attendu = flottant(exact)
        erreur = np.linalg.norm(resultat["champ"]-attendu, axis=0)
        normes = np.linalg.norm(attendu, axis=0)
        np.testing.assert_array_less(erreur, tolerance*normes+1e-30)
        self.assertFalse(resultat["certification_machine"])

    def test_modele_couple_oblique_non_modal_et_charges_multiples(self):
        d, m, i, s, b, phi, qr = donnees()
        reduction = KrylovContraint(qr, m, b, phi, 4., .5, blocs=3, max_directions=3)
        self.assertEqual(reduction.taille_complement_reduit, 3)
        # Ces couplages sont physiques et ne doivent pas être supprimés.
        self.assertGreater(np.linalg.norm(m[np.ix_(i, s)]), .01)
        ki = d[:, i].T@d[:, i]
        modal = np.column_stack(((ki@phi)[:, 0], (m[np.ix_(i, i)]@phi)[:, 0]))
        self.assertEqual(np.linalg.matrix_rank(modal), 2)
        forces = self.forces()
        for omega in (0., .125, .5):
            with self.subTest(omega=omega):
                rep = reduction.reponses(omega, forces)
                self.comparer(rep, oracle(d, m, s, omega, forces))
                self.assertGreater(np.linalg.norm(rep["champ"][:, 3]), 0.)
                np.testing.assert_array_equal(rep["champ"][:, 4], np.zeros(d.shape[1]))

    def test_pole_interieur_traversable_et_vrai_pole_global(self):
        from controle_complement import ControleComplement
        for global_singulier in (False, True):
            with self.subTest(global_singulier=global_singulier):
                d = np.array([[1., 0., 0. if global_singulier else 1.],
                              [0., 2., 1.], [0., 0., 1.]])
                m, b = np.eye(3), np.array([[1.], [0.]])
                qr = CondensationEnergie(csr_matrix(d), np.arange(2), np.array([2]), np.eye(1))
                reduction = KrylovContraint(qr, m, b, b, 3., 1.25, blocs=2, max_directions=1)
                controle = ControleComplement(reduction)
                forces = np.array([[1., -.5, 2.**-30]])
                # K_ii-I est singulier dans les deux cas ; seul le premier
                # système complet est inversible grâce au couplage retenu.
                self.assertEqual(F(d[0, 0])**2-F(1), 0)
                if global_singulier:
                    with self.assertRaises(ValueError):
                        oracle(d, m, [2], 1., forces)
                    for solveur in (reduction, controle):
                        with self.assertRaises(np.linalg.LinAlgError):
                            solveur.reponses(1., forces)
                else:
                    exact = oracle(d, m, [2], 1., forces)
                    for solveur in (reduction, controle):
                        self.comparer(solveur.reponses(1., forces), exact)

    def test_masse_complete_psd_et_port_sans_masse(self):
        from controle_complement import ControleComplement
        d = np.array([[1., 0., 1.], [0., 2., 1.], [0., 0., 1.]])
        m = np.diag([1., 1., 0.])
        b = np.array([[1.], [0.]])
        qr = CondensationEnergie(csr_matrix(d), np.arange(2), np.array([2]), np.eye(1))
        reduction = KrylovContraint(qr, m, b, b, 3., 1.25, blocs=2, max_directions=1)
        controle = ControleComplement(reduction)
        self.assertEqual(m[2, 2], 0.)
        forces = np.array([[1., -.5, 2.**-30, 0.]])
        for omega in (0., .5, 1.):
            with self.subTest(omega=omega):
                exact = oracle(d, m, [2], omega, forces)
                for solveur in (reduction, controle):
                    self.comparer(solveur.reponses(omega, forces), exact)
                rep = controle.reponses(omega, forces)
                self.assertGreater(rep["marge"], 0.)
                for j in range(forces.shape[1]):
                    x = [F(float(v)) for v in rep["champ"][:, j]]
                    norme2 = x[0]**2+x[1]**2  # Le port ne contribue aucune masse.
                    calculee2 = F(float(rep["bornes"]["masse"]["normes"][j]))**2
                    self.assertLessEqual(abs(calculee2-norme2), F(1, 10**12)*norme2)

    def test_complement_krylov_vide_dans_un_modele_decouple(self):
        from controle_complement import ControleComplement
        d, m = np.diag([1., 2., 3.]), np.diag([1., 1., 0.])
        b = np.array([[1.], [0.]])
        qr = CondensationEnergie(csr_matrix(d), np.arange(2), np.array([2]), np.eye(1))
        reduction = KrylovContraint(qr, m, b, b, 4., .5, blocs=2, max_directions=1)
        controle = ControleComplement(reduction)
        self.assertEqual(reduction.v.shape, (2, 0))
        self.assertEqual(reduction.taille_complement_reduit, 0)
        self.assertEqual(reduction.taille_reduite, 2)
        self.assertEqual(controle.delta, 0.)
        forces = np.array([[1., -.5, 2.**-30, 0.]])
        for omega in (0., .25, .5):
            with self.subTest(omega=omega):
                exact = oracle(d, m, [2], omega, forces)
                for solveur in (reduction, controle):
                    rep = solveur.reponses(omega, forces)
                    self.comparer(rep, exact)
                    self.assertEqual(rep["taille_complement_reduit"], 0)
                    np.testing.assert_array_equal(rep["champ"][:2], np.zeros((2, 4)))
                    np.testing.assert_allclose(rep["champ"][2], forces[0]/9, rtol=2e-15, atol=0.)
                    self.assertFalse(rep["certification_machine"])
                rep = controle.reponses(omega, forces)
                self.assertEqual(rep["taille_conservee"], 2)
                np.testing.assert_array_equal(rep["bornes"]["masse"]["normes"], np.zeros(4))
        # delta=0 décrit l'absence d'erreur de réduction en exact. Il ne
        # certifie pas l'arrondi binary64 du déplacement non représentable 1/9.

    def test_reparation_finale_apres_normalisation_preserve_les_covecteurs(self):
        d, m, i, _, _, _, qr = donnees(8)
        b = np.array([[1.], [2.], [-1.], [3.], [-2.], [1.], [4.], [-3.]])
        reduction = KrylovContraint(qr, m, b, np.ones_like(b), 4., .5,
                                     blocs=1, max_directions=7)
        taille = reduction.v.shape[1]
        # Une direction presque déjà présente : l'orthogonalisation laisse
        # une petite composante qui doit être normalisée, sans amplifier un
        # défaut résiduel de B.T en violation significative des contraintes.
        candidat = reduction.v[:, [0]]+2.**-30*np.sin(np.arange(8)+.25)[:, None]
        nouveau = reduction._ajouter(candidat)
        self.assertEqual(nouveau.shape, (8, 1))
        self.assertEqual(reduction.v.shape[1], taille+1)
        energie_nouvelle = float(np.linalg.norm(d[:, i]@nouveau))
        self.assertGreater(energie_nouvelle, .99)
        self.assertLess(energie_nouvelle, 1.01)
        seuil = 128*np.finfo(float).eps*np.linalg.norm(b)*np.linalg.norm(nouveau)
        self.assertLess(float(np.linalg.norm(b.T@nouveau)), seuil)

    def test_enrichissement_reutilise_le_facteur_et_perime_le_controle(self):
        from controle_complement import ControleComplement
        d, m, i, s, b, phi, qr = donnees(8)
        reduction = KrylovContraint(qr, m, b, phi, 4., .5, blocs=1, max_directions=7)
        controle = ControleComplement(reduction)
        forces = self.forces()
        exact = flottant(oracle(d, m, s, .5, forces))
        avant = np.linalg.norm(d@(reduction.reponses(.5, forces)["champ"]-exact))
        taille, compte = reduction.taille_complement_reduit, reduction.resolutions_statiques
        facteur = reduction.inverse
        self.assertTrue(reduction.enrichir())
        self.assertIs(reduction.inverse, facteur)
        self.assertIs(reduction.qr, qr)
        self.assertGreater(reduction.resolutions_statiques, compte)
        self.assertGreater(reduction.taille_complement_reduit, taille)
        with self.assertRaisesRegex(ValueError, "périmé"):
            controle.reponses(.5, forces)
        apres = np.linalg.norm(d@(reduction.reponses(.5, forces)["champ"]-exact))
        self.assertLess(apres, avant/10)
        while reduction.enrichir():
            self.assertLessEqual(reduction.taille_complement_reduit, 7)
        self.comparer(reduction.reponses(.5, forces), oracle(d, m, s, .5, forces))

    def test_majorants_complement_incomplet_face_aux_erreurs_exactes(self):
        from controle_complement import ControleComplement
        d, m, i, s, b, phi, qr = donnees(8)
        reduction = KrylovContraint(qr, m, b, phi, 4., .5, blocs=1, max_directions=7)
        self.assertLess(reduction.taille_complement_reduit, len(i)-b.shape[1])
        controle = ControleComplement(reduction)
        forces = self.forces()
        matrices = dict(masse=rationnel(m), deformation=gram(rationnel(d)))
        erreur_non_nulle = False
        for omega in (0., .25, .5):
            rep = controle.reponses(omega, forces)
            self.assertGreater(rep["marge"], 0)
            exact = oracle(d, m, s, omega, forces)
            candidat = rationnel(rep["champ"])
            for j in range(forces.shape[1]):
                x = [ligne[j] for ligne in candidat]
                erreur = [ligne[j]-ref[j] for ligne, ref in zip(candidat, exact)]
                for nom, matrice in matrices.items():
                    erreur2, norme2 = energie(erreur, matrice), energie(x, matrice)
                    borne = rep["bornes"][nom]["absolues"][j]
                    self.assertIsNotNone(borne)
                    # Comparaison de l'erreur exacte du champ binary64 renvoyé,
                    # sans calculer les petites différences en doubles.
                    self.assertLessEqual(erreur2, F(float(borne))**2,
                                         (omega, j, nom, float(erreur2), borne))
                    norme_calculee2 = F(float(rep["bornes"][nom]["normes"][j]))**2
                    self.assertLessEqual(abs(norme_calculee2-norme2),
                                         F(1, 10**12)*norme2,
                                         (omega, j, nom))
                    erreur_non_nulle |= erreur2 > F(1, 10**24)
            self.assertFalse(rep["certification_machine"])
        self.assertTrue(erreur_non_nulle)

    def test_budgets_et_forces_invalides_refusent(self):
        _, m, _, _, b, phi, qr = donnees(8)
        reduction = KrylovContraint(qr, m, b, phi, 4., .5, blocs=4, max_directions=2)
        self.assertEqual(reduction.taille_complement_reduit, 2)
        self.assertFalse(reduction.enrichir())
        self.assertEqual(reduction.taille_complement_reduit, 2)
        for options in (dict(blocs=0), dict(max_directions=0), dict(blocs=True)):
            with self.assertRaises(ValueError):
                KrylovContraint(qr, m, b, phi, 4., .5, **options)
        for forces in (np.ones(2), np.empty((2, 0)), np.ones((1, 2)),
                       np.ones((2, 1), dtype=complex), np.full((2, 1), np.inf)):
            with self.assertRaises(ValueError):
                reduction.reponses(.25, forces)
        for omega in (-1., .6, np.nan, 1j):
            with self.assertRaises(ValueError):
                reduction.reponses(omega, self.forces())


if __name__ == "__main__":
    unittest.main()
