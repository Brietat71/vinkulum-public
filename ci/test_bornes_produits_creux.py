"""Contre-épreuve rationnelle d'une borne de produit matriciel publiée.

Ce fichier ne teste ni INTLAB ni son implémentation. Il confronte le lemme
2.3, (2.6), du manuscrit TUHH Part I à un calcul binary64 explicite. Une
borne indépendante, volontairement conservative, sert de témoin positif.
Voir docs/VERIFICATION_CREUSE_2026.md pour les hypothèses et la provenance.
"""
from fractions import Fraction as F
import math
import sys
import unittest


U = F(1, 2**53)
MIN_NORMAL = F.from_float(sys.float_info.min)
MAX_FINITE = F.from_float(sys.float_info.max)


def _normal_ou_zero(x):
    if x and not MIN_NORMAL <= abs(x) <= MAX_FINITE:
        raise ValueError("sous-flux ou débordement hors des hypothèses")


def _rationnels(a):
    if not a or not a[0] or any(len(row) != len(a[0]) for row in a):
        raise ValueError("matrice rectangulaire non vide requise")
    resultat = []
    for row in a:
        ligne = []
        for x in row:
            if not isinstance(x, float) or not math.isfinite(x):
                raise ValueError("données binary64 finies requises")
            exact = F.from_float(x)
            _normal_ou_zero(exact)
            ligne.append(exact)
        resultat.append(ligne)
    return resultat


def _dimensions(a, b):
    if len(a[0]) != len(b):
        raise ValueError("dimensions incompatibles")


def produit_exact(a, b):
    """Interprète exactement les coefficients binary64 d'entrée."""
    aa, bb = _rationnels(a), _rationnels(b)
    _dimensions(aa, bb)
    return [[sum((x*y for x, y in zip(row, col)), F(0))
             for col in zip(*bb)] for row in aa]


def produit_binary64(a, b, *, inverse=False):
    """Produits puis sommes séquentiels ; aucun FMA ni appel à un BLAS.

    Les opérations exactes intermédiaires doivent rester normales ou
    nulles. Ces gardes conservatrices rendent explicites les hypothèses
    du modèle relatif d'arrondi utilisé par le témoin ci-dessous.
    """
    aa, bb = _rationnels(a), _rationnels(b)
    _dimensions(aa, bb)
    indices = list(range(len(b)))
    if inverse:
        indices.reverse()
    resultat = []
    for i in range(len(a)):
        ligne = []
        for j in range(len(b[0])):
            somme = 0.0
            for k in indices:
                if not aa[i][k] or not bb[k][j]:
                    continue
                _normal_ou_zero(aa[i][k]*bb[k][j])
                terme = a[i][k]*b[k][j]
                _normal_ou_zero(F.from_float(somme)+F.from_float(terme))
                somme = somme+terme
            ligne.append(somme)
        resultat.append(ligne)
    return resultat


def gamma(k):
    """γ_k = ku/(1-ku), rationnel exact ; ku < 1 obligatoire."""
    if isinstance(k, bool) or not isinstance(k, int) or k < 0 or k*U >= 1:
        raise ValueError("entier k >= 0 avec k*u < 1 requis")
    return k*U/(1-k*U)


def majorant_independant(a, b):
    """W_ij=γ_(2k) Σ|a_il b_lj| ; k compte les produits non nuls.

    La validité suppose le calcul décrit par produit_binary64, en
    arrondi au plus proche et sans sous-flux/débordement. Chaque chemin
    de calcul contient au plus 2k opérations arrondies. Le dénominateur
    de γ n'est pas omis au profit d'une approximation par k*u.

    Renvoie W et min(||W||_F², ||W||_1 ||W||_inf), majorant de ||E||_2².
    Tous les nombres renvoyés sont des Fraction ; aucune racine approchée.
    """
    aa, bb = _rationnels(a), _rationnels(b)
    _dimensions(aa, bb)
    w = []
    for row in aa:
        ligne = []
        for col in zip(*bb):
            termes = [abs(x*y) for x, y in zip(row, col) if x and y]
            ligne.append(gamma(2*len(termes))*sum(termes, F(0)))
        w.append(ligne)
    frobenius2 = sum((x*x for row in w for x in row), F(0))
    norme_inf = max(sum(row, F(0)) for row in w)
    norme_1 = max(sum(col, F(0)) for col in zip(*w))
    return w, min(frobenius2, norme_1*norme_inf)


class BornesProduitsCreux(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if (sys.float_info.radix, sys.float_info.mant_dig,
                sys.float_info.rounds) != (2, 53, 1):
            raise RuntimeError("binary64 en arrondi au plus proche requis")

    @staticmethod
    def donnees(d=2.0**-10):
        q = 2.0**-27
        return [[1.0, q], [0.0, d]], [[d, 1.0], [0.0, q]]

    def erreur(self, a, b, **kwargs):
        exact = produit_exact(a, b)
        calcule = produit_binary64(a, b, **kwargs)
        return [[F.from_float(x)-y for x, y in zip(rx, ry)]
                for rx, ry in zip(calcule, exact)]

    def test_facteurs_inversibles_et_arrondi_exact(self):
        a, b = self.donnees()
        q, d = F(1, 2**27), F(1, 2**10)
        aa, bb = _rationnels(a), _rationnels(b)
        self.assertEqual(aa[0][0]*aa[1][1]-aa[0][1]*aa[1][0], d)
        self.assertEqual(bb[0][0]*bb[1][1]-bb[0][1]*bb[1][0], d*q)
        self.assertEqual(produit_exact(a, b), [[d, 1+q*q], [F(0), d*q]])
        # 1+2^-54 est plus près de 1 que de son successeur 1+2^-52.
        self.assertLess(q*q, U)
        self.assertEqual(F.from_float(math.nextafter(1.0, math.inf))-1, 2*U)
        for inverse in (False, True):
            with self.subTest(inverse=inverse):
                self.assertEqual(self.erreur(a, b, inverse=inverse),
                                 [[F(0), -q*q], [F(0), F(0)]])

    def test_inegalite_publiee_refutee_sans_racine_flottante(self):
        a, b = self.donnees()
        aa, bb = _rationnels(a), _rationnels(b)
        mu = [sum(x != 0 for x in row) for row in aa]
        nu = [sum(x != 0 for x in col) for col in zip(*bb)]
        rho2 = [sum((x*x for x in row), F(0)) for row in aa]
        sigma2 = [sum((x*x for x in col), F(0)) for col in zip(*bb)]
        q, d = F(1, 2**27), F(1, 2**10)
        self.assertEqual((mu, nu), ([2, 1], [1, 2]))
        self.assertEqual((rho2, sigma2), ([1+q*q, d*d], [d*d, 1+q*q]))
        # Le membre droit publié est 2*u*d*sqrt(1+q²).
        # 1+q² < 4 démontre qu'il est STRICTEMENT inférieur à 4*u*d.
        self.assertLess(1+q*q, 4)
        majorant_strict_du_membre_droit = 4*U*d
        self.assertLess(majorant_strict_du_membre_droit, q*q)
        # L'erreur n'a qu'une entrée non nulle : ||E||_2 = q² exactement.
        self.assertEqual(self.erreur(a, b), [[F(0), -q*q], [F(0), F(0)]])

    def test_variante_membre_droit_nul(self):
        q = 2.0**-27
        a, b = [[1.0, q], [0.0, 0.0]], [[0.0, 1.0], [0.0, q]]
        mu = [sum(x != 0 for x in row) for row in a]
        nu = [sum(x != 0 for x in col) for col in zip(*b)]
        self.assertEqual([min(x, y) for x, y in zip(mu, nu)], [0, 0])
        self.assertEqual(self.erreur(a, b),
                         [[F(0), -F(1, 2**54)], [F(0), F(0)]])

    def test_borne_correcte_contient_les_erreurs(self):
        cas = [self.donnees(),
               ([[1.0, 2.0**-27], [0.0, 0.0]],
                [[0.0, 1.0], [0.0, 2.0**-27]]),
               ([[1.0, 2.0**-54, -1.0], [2.0**-27, -2.0, 3.0]],
                [[1.0, 2.0**-27], [1.0, -3.0], [1.0, 0.5]])]
        for indice, (a, b) in enumerate(cas):
            w, borne2 = majorant_independant(a, b)
            for inverse in (False, True):
                with self.subTest(cas=indice, inverse=inverse):
                    e = self.erreur(a, b, inverse=inverse)
                    for re, rw in zip(e, w):
                        for erreur, borne in zip(re, rw):
                            self.assertLessEqual(abs(erreur), borne)
                    # Pour ces matrices à deux colonnes, vérifier exactement
                    # borne2*I-E.T*E PSD par ses trois mineurs principaux.
                    gram = [[sum((row[i]*row[j] for row in e), F(0))
                             for j in range(2)] for i in range(2)]
                    p, r = borne2-gram[0][0], borne2-gram[1][1]
                    self.assertGreaterEqual(p, 0)
                    self.assertGreaterEqual(r, 0)
                    self.assertGreaterEqual(p*r-gram[0][1]**2, 0)

    def test_zero_et_produit_rectangulaire(self):
        a, b = [[0.0, 0.0, 0.0]], [[1.0], [-2.0], [3.0]]
        self.assertEqual(produit_binary64(a, b), [[0.0]])
        self.assertEqual(majorant_independant(a, b), ([[F(0)]], F(0)))

    def test_gamma_et_hypothese_sur_le_nombre_operations(self):
        self.assertEqual(gamma(0), 0)
        self.assertEqual(gamma(4), 4*U/(1-4*U))
        self.assertGreater(gamma(4), 4*U)
        for k in (-1, True, 2.0, 2**53, 2**53+1):
            with self.subTest(k=k), self.assertRaises(ValueError):
                gamma(k)

    def test_operations_hors_modele_refusees(self):
        cas = [([[sys.float_info.min]], [[0.5]]),
               ([[sys.float_info.max]], [[2.0]]),
               ([[sys.float_info.max, sys.float_info.max]], [[1.0], [1.0]]),
               ([[sys.float_info.min, -math.nextafter(sys.float_info.min, math.inf)]],
                [[1.0], [1.0]])]
        for a, b in cas:
            with self.subTest(a=a, b=b), self.assertRaisesRegex(ValueError, "hors des hypothèses"):
                produit_binary64(a, b)


if __name__ == "__main__":
    unittest.main()
