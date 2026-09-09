"""Contrôleur à facteurs : champs physiques et contre-épreuves rationnelles."""
from fractions import Fraction as F
import unittest

import numpy as np
from scipy.sparse import csr_matrix

from controle_complement import ControleComplement
from controle_facteurs import ControleFacteurs, ImageFaibleRang, normes_majorees
from krylov_contraint import KrylovContraint
from vinkulum._ports.condensation_energie import CondensationEnergie
from test_krylov_contraint import donnees, oracle, rationnel, flottant, energie
from test_inertie_complement import inertie_fraction
from test_trace_complement import gram, produit


MODES = ('gram', 'qr')


class ControleFacteursPhysique(unittest.TestCase):
    def forces(self):
        e = 2.**-40
        return np.array([[1., 0., 1., e, 0., 1., 1.],
                         [0., 1., 1., -e, 0., -1., -1.+2.**-40]])

    def identiques(self, ancien, nouveau):
        for cle in ('champ', 'coordonnees'):
            np.testing.assert_array_equal(nouveau[cle], ancien[cle])
        for nom in ('masse', 'deformation'):
            np.testing.assert_array_equal(nouveau['bornes'][nom]['normes'],
                                          ancien['bornes'][nom]['normes'])
        self.assertEqual(nouveau['residu_coefficients'], ancien['residu_coefficients'])
        self.assertEqual(nouveau['erreur_coefficients'], ancien['erreur_coefficients'])
        self.assertFalse(nouveau['certification_machine'])

    def proche_oracle(self, rep, exact, tolerance=3e-11):
        attendu = flottant(exact)
        erreurs = np.linalg.norm(rep['champ']-attendu, axis=0)
        normes = np.linalg.norm(attendu, axis=0)
        np.testing.assert_array_less(erreurs, tolerance*normes+1e-32)

    def images_multiples(self):
        d, m, i, s, _, _, qr = donnees()
        b = np.array([[1., 0.], [0., 1.], [2., -1.], [1., 3.]])
        phi = b+np.array([[0., .125], [.25, 0.], [0., -.25], [.125, 0.]])
        reduction = KrylovContraint(qr, m, b, phi, 4., .5, blocs=3, max_directions=2)
        return d, m, i, s, b, phi, reduction

    def test_masse_couplee_oblique_non_modale_deux_contraintes(self):
        d, m, i, s, b, phi, reduction = self.images_multiples()
        self.assertGreater(np.linalg.norm(m[np.ix_(i, s)]), .01)
        self.assertEqual(reduction.taille_complement_reduit, len(i)-b.shape[1])
        ancien = ControleComplement(reduction)
        forces = self.forces()
        for extension in MODES:
            nouveau = ControleFacteurs(reduction, extension=extension)
            self.assertFalse(nouveau.facteurs['certification_machine'])
            for nom in ('masse', 'deformation'):
                self.assertEqual(nouveau.facteurs['reparation'][nom]['rang_facteur'], 2)
            for omega in (0., .125, .5):
                with self.subTest(extension=extension, omega=omega):
                    rep = nouveau.reponses(omega, forces)
                    self.identiques(ancien.reponses(omega, forces), rep)
                    self.proche_oracle(rep, oracle(d, m, s, omega, forces))
                    np.testing.assert_array_equal(rep['champ'][:, 4], np.zeros(d.shape[1]))
                    self.assertGreater(np.linalg.norm(rep['champ'][:, 3]), 0.)

    def test_pole_interieur_traversable_et_pole_global_refuse(self):
        for singulier in (False, True):
            d = np.array([[1., 0., 0. if singulier else 1.],
                          [0., 2., 1.], [0., 0., 1.]])
            m, b = np.eye(3), np.array([[1.], [0.]])
            qr = CondensationEnergie(csr_matrix(d), np.arange(2), np.array([2]), np.eye(1))
            reduction = KrylovContraint(qr, m, b, b, 3., 1.25, blocs=2, max_directions=1)
            ancien = ControleComplement(reduction)
            forces = np.array([[1., -.5, 2.**-40]])
            for extension in MODES:
                with self.subTest(singulier=singulier, extension=extension):
                    nouveau = ControleFacteurs(reduction, extension=extension)
                    if singulier:
                        with self.assertRaises(ValueError):
                            oracle(d, m, [2], 1., forces)
                        for controle in (ancien, nouveau):
                            with self.assertRaises(np.linalg.LinAlgError):
                                controle.reponses(1., forces)
                    else:
                        rep = nouveau.reponses(1., forces)
                        self.identiques(ancien.reponses(1., forces), rep)
                        self.proche_oracle(rep, oracle(d, m, [2], 1., forces))

    def test_masse_psd_port_sans_masse_et_complement_vide(self):
        for decouple in (False, True):
            d = (np.diag([1., 2., 3.]) if decouple else
                 np.array([[1., 0., 1.], [0., 2., 1.], [0., 0., 1.]]))
            m, b = np.diag([1., 1., 0.]), np.array([[1.], [0.]])
            qr = CondensationEnergie(csr_matrix(d), np.arange(2), np.array([2]), np.eye(1))
            reduction = KrylovContraint(qr, m, b, b, 3., .5, blocs=2, max_directions=1)
            ancien = ControleComplement(reduction)
            forces = np.array([[1., -.5, 2.**-40, 0.]])
            self.assertEqual(reduction.taille_complement_reduit, 0 if decouple else 1)
            for extension in MODES:
                nouveau = ControleFacteurs(reduction, extension=extension)
                for omega in (0., .25, .5):
                    with self.subTest(decouple=decouple, extension=extension, omega=omega):
                        rep = nouveau.reponses(omega, forces)
                        self.identiques(ancien.reponses(omega, forces), rep)
                        self.proche_oracle(rep, oracle(d, m, [2], omega, forces))
                        for j in range(forces.shape[1]):
                            x = [F(float(v)) for v in rep['champ'][:, j]]
                            attendu = x[0]**2+x[1]**2
                            observe = F(float(rep['bornes']['masse']['normes'][j]))**2
                            self.assertLessEqual(abs(attendu-observe), F(1, 10**12)*attendu)
                        if decouple:
                            np.testing.assert_array_equal(rep['champ'][:2], np.zeros((2, 4)))
                            np.testing.assert_array_equal(rep['bornes']['masse']['normes'], np.zeros(4))
                            self.assertEqual(rep['taille_complement_reduit'], 0)
        # Un majorant de réduction nul ne certifie pas l'arrondi de 1/9.

    def test_majorants_complement_incomplet_face_a_fraction(self):
        d, m, i, s, b, phi, qr = donnees(8)
        reduction = KrylovContraint(qr, m, b, phi, 4., .5, blocs=1, max_directions=7)
        self.assertLess(reduction.taille_complement_reduit, len(i)-b.shape[1])
        ancien = ControleComplement(reduction)
        forces = self.forces()
        matrices = dict(masse=rationnel(m), deformation=gram(rationnel(d)))
        non_nul = False
        for extension in MODES:
            controle = ControleFacteurs(reduction, extension=extension)
            for omega in (0., .25, .5):
                with self.subTest(extension=extension, omega=omega):
                    rep = controle.reponses(omega, forces)
                    self.identiques(ancien.reponses(omega, forces), rep)
                    self.assertGreater(rep['marge'], 0.)
                    exact = oracle(d, m, s, omega, forces)
                    candidat = rationnel(rep['champ'])
                    for j in range(forces.shape[1]):
                        erreur = [ligne[j]-ref[j] for ligne, ref in zip(candidat, exact)]
                        for nom, matrice in matrices.items():
                            erreur2 = energie(erreur, matrice)
                            borne = rep['bornes'][nom]['absolues'][j]
                            self.assertIsNotNone(borne)
                            self.assertLessEqual(erreur2, F(float(borne))**2,
                                                 (extension, omega, j, nom))
                            non_nul |= erreur2 > F(1, 10**24)
        self.assertTrue(non_nul)

    def test_controle_perime_et_entrees_invalides(self):
        _, m, _, _, b, phi, qr = donnees(8)
        reduction = KrylovContraint(qr, m, b, phi, 4., .5, blocs=1, max_directions=7)
        controles = [ControleFacteurs(reduction, extension=mode) for mode in MODES]
        for controle in controles:
            for forces in (np.ones(2), np.empty((2, 0)), np.ones((1, 2)),
                           np.ones((2, 1), dtype=complex), np.full((2, 1), np.inf)):
                with self.assertRaises(ValueError):
                    controle.reponses(.25, forces)
            for omega in (-1., .6, np.nan, 1j):
                with self.assertRaises(ValueError):
                    controle.reponses(omega, self.forces())
        self.assertTrue(reduction.enrichir())
        for controle in controles:
            with self.assertRaisesRegex(ValueError, 'périmé'):
                controle.reponses(.25, self.forces())
        with self.assertRaisesRegex(ValueError, 'extension'):
            ControleFacteurs(reduction, extension='inconnue')

    def test_minorant_trop_conservateur_ne_fabrique_pas_de_majorants(self):
        d, m, _, s, b, phi, qr = donnees(8)
        # Le modèle admet déjà lambda>=4 ; ce minorant valide presque au
        # bord de la bande rend le contrôle incomplet trop pessimiste.
        reduction = KrylovContraint(qr, m, b, phi, .25000001, .5,
                                     blocs=1, max_directions=1)
        forces = np.eye(2)
        self.assertTrue(np.all(np.isfinite(flottant(oracle(d, m, s, .5, forces)))))
        reference = ControleComplement(reduction).reponses(.5, forces)
        self.assertLess(reference['marge'], 0.)
        for extension in MODES:
            with self.subTest(extension=extension):
                rep = ControleFacteurs(reduction, extension=extension).reponses(.5, forces)
                self.identiques(reference, rep)
                self.assertLess(rep['marge'], 0.)
                self.assertIsNone(rep['borne_coordonnees'])
                for nom in ('masse', 'deformation'):
                    self.assertIsNone(rep['bornes'][nom]['absolues'])
                    self.assertEqual(rep['bornes'][nom]['relatives'], [None, None])

    def test_majoration_operateur_dx_avec_gram_fraction(self):
        _, m, _, _, b, phi, qr = donnees(8)
        reduction = KrylovContraint(qr, m, b, phi, 4., .5, blocs=1, max_directions=7)
        for extension in MODES:
            controle = ControleFacteurs(reduction, extension=extension)
            mu = .5**2/reduction.lambda_complement
            y = np.linalg.solve(np.eye(len(controle.theta))-mu*controle.theta,
                                controle.a0+mu*controle.a1)
            essais = [y, np.zeros_like(y), y*2.**30,
                      y+2.**-30*np.sin(np.arange(y.size)).reshape(y.shape)]
            for j, yy in enumerate(essais):
                with self.subTest(extension=extension, essai=j):
                    x = reduction.t.copy()
                    x[reduction.i] -= controle.w@yy
                    dx = reduction.d@x
                    borne = controle._norme_extension_d(yy, dx, dx.T@dx)
                    self.assertTrue(np.isfinite(borne))
                    exact = gram(rationnel(dx))
                    # b² I − DXᵀ DX PSD équivaut à ||DX||₂≤b. L'oracle
                    # d'inertie rationnel n'emploie aucune SVD candidate.
                    reste = [[(F(float(borne))**2 if i == k else F(0))-exact[i][k]
                              for k in range(len(exact))] for i in range(len(exact))]
                    self.assertEqual(inertie_fraction(reste)[1], 0)
                    self.assertGreater(borne, 0.)
                    self.assertLess(borne, 2*np.linalg.norm(dx, 2))

    def test_image_faible_rang_defaut_reel_et_compensations(self):
        gauche = np.array([[1.], [2.], [-1.], [.5]])
        droite = np.array([[1., -1., 2.]])
        # La petite perturbation ne provient pas des facteurs proposés.
        # Sans son défaut, la première combinaison serait faussement nulle.
        image = gauche@droite
        image[2, 1] += 2.**-30
        e = 2.**-40
        v = np.array([[1., 1., e, 0., 1.], [1., -1., e, 0., 1.],
                      [0., 0., 0., 0., e]])
        for puissance in (0, -300, 300):
            with self.subTest(reechelonnement=puissance):
                facteur = ImageFaibleRang(image, gauche*2.**puissance, droite*2.**-puissance)
                bornes = facteur.normes(v)
                exact = produit(rationnel(image), rationnel(v))
                flottantes = rationnel(image@v)
                for j, borne in enumerate(bornes):
                    for reference in (exact, flottantes):
                        norme2 = sum((ligne[j]**2 for ligne in reference), F(0))
                        self.assertLessEqual(norme2, F(float(borne))**2)
                self.assertGreater(bornes[0], 0.)
                self.assertLess(bornes[0], 1e-6)
                self.assertEqual(bornes[3], 0.)
                self.assertEqual(facteur.diagnostic['rang_facteur'], 1)
                self.assertGreater(max(facteur.diagnostic['defaut_colonnes']), 0.)
                self.assertFalse(facteur.diagnostic['certification_machine'])

    def test_norme_audit_stridee_ne_perd_pas_les_petits_carres(self):
        # Tous les produits sont normaux : ce cas ne dépend pas d'un
        # sous-flux. La réduction stridée sur128 lignes perdait127 carrés
        # 2^-54 ajoutés à1 ; le seul gamma du nombre de colonnes ne suffit pas.
        image = np.full((128, 2), 2.**-27)
        image[0] = 1.
        facteur = ImageFaibleRang(image, np.zeros((128, 1)), np.zeros((1, 2)))
        borne = float(facteur.normes(np.array([[1.], [0.]]))[0])
        exact2 = F(1)+127*F(2)**-54
        self.assertGreater(exact2, F(float.fromhex('0x1.0000000000002p+0'))**2)
        self.assertTrue(np.isfinite(borne))
        self.assertLessEqual(exact2, F(borne)**2)
        self.assertLess(borne, 1.+1e-10)

    def test_normes_majorees_echelles_extremes_oracle_fraction(self):
        for n in (1, 2, 3, 16, 128, 3072):
            for exposant in (-1074, -1000, -600, -300, 0, 300, 600, 1000):
                with self.subTest(n=n, exposant=exposant):
                    x = np.ldexp(np.ones(n), exposant)
                    x[1:] *= .125
                    borne = float(normes_majorees(x))
                    self.assertTrue(np.isfinite(borne))
                    # L'oracle porte sur les composantes effectivement
                    # stockées, y compris leurs éventuels sous-flux à zéro.
                    exact2 = sum((F(float(v))**2 for v in x), F(0))
                    self.assertGreaterEqual(F(borne)**2, exact2)

    def test_gram_nul_ou_subnormal_ne_donne_pas_zero_pour_dx_non_nul(self):
        _, m, _, _, b, phi, qr = donnees()
        reduction = KrylovContraint(qr, m, b, phi, 4., .5, blocs=1, max_directions=3)
        controle = ControleFacteurs(reduction, extension='gram')
        lignes, colonnes = reduction.d.shape[0], reduction.t.shape[1]
        y = np.zeros((len(controle.theta), colonnes))
        for valeur in (2.**-600, 1.2*2.**-537):
            with self.subTest(valeur=valeur.hex()):
                dx = np.zeros((lignes, colonnes))
                dx[np.arange(colonnes), np.arange(colonnes)] = valeur
                g = dx.T@dx
                self.assertTrue(np.all(np.diag(g) < np.finfo(float).tiny))
                borne = controle._norme_extension_d(y, dx, g)
                self.assertTrue(np.isfinite(borne))
                self.assertGreater(borne, 0.)
                # Les colonnes sont orthogonales : ||DX||2 = valeur en exact.
                self.assertGreaterEqual(F(float(borne)), F(valeur))
                self.assertLess(borne, 2*valeur)


if __name__ == '__main__':
    unittest.main()
