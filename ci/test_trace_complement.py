"""Contre-épreuves exactes de la trace sur un complément contraint.

Petites matrices Fraction uniquement : aucune promesse de solveur certifié
en binary64. La voie Z explicite sert de référence à l'identité sans Z.
"""
from fractions import Fraction as F
from itertools import combinations
import unittest


def matrice(a):
    return [[F(v) for v in ligne] for ligne in a]


def transpose(a):
    return [list(ligne) for ligne in zip(*a)]


def produit(a, b):
    return [[sum((x*y for x, y in zip(ligne, col)), F(0))
             for col in zip(*b)] for ligne in a]


def difference(a, b):
    return [[x-y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def echelle(a, c):
    return [[c*v for v in ligne] for ligne in a]


def identite(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def inverse(a):
    n = len(a)
    aug = [ligne.copy()+un for ligne, un in zip(a, identite(n))]
    for i in range(n):
        pivot = next((j for j in range(i, n) if aug[j][i]), None)
        if pivot is None:
            raise ValueError("matrice singulière : rang insuffisant")
        aug[i], aug[pivot] = aug[pivot], aug[i]
        div = aug[i][i]
        aug[i] = [v/div for v in aug[i]]
        for j in range(n):
            if j != i:
                c = aug[j][i]
                aug[j] = [x-c*y for x, y in zip(aug[j], aug[i])]
    return [ligne[n:] for ligne in aug]


def determinant(a):
    if len(a) == 1:
        return a[0][0]
    return sum(((-1)**j*a[0][j]*determinant(
        [ligne[:j]+ligne[j+1:] for ligne in a[1:]])
        for j in range(len(a))), F(0))


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


def gram(a):
    return produit(transpose(a), a)


def restreindre(a, z):
    return produit(transpose(z), produit(a, z))


def inverse_contrainte(k, m, b):
    ki = inverse(k)
    x = produit(ki, b)
    hi = inverse(produit(transpose(b), x))
    correction = produit(produit(x, hi), transpose(x))
    s = difference(ki, correction)
    tau = trace(produit(m, ki))
    chi = trace(produit(hi, restreindre(m, x)))
    return s, tau-chi, tau, chi


class TraceComplement(unittest.TestCase):
    def donnees(self):
        k = gram(matrice([[2, 1, 0, 1], [0, 3, 1, 0],
                          [0, 0, 2, 1], [0, 0, 0, 4]]))
        m = gram(matrice([[1, 2, 0, 0], [0, 2, 1, 0],
                          [0, 0, 3, 1], [0, 0, 0, 2]]))
        b = matrice([[1, 0], [0, 1], [2, -1], [1, 3]])
        z = matrice([[-2, -1], [1, -3], [1, 0], [0, 1]])
        return k, m, b, z

    def verifier_psd(self, a):
        self.assertEqual(a, transpose(a))
        # Pour une matrice réelle symétrique, tous les mineurs principaux
        # non négatifs sont un critère exact de semi-définie positivité.
        for taille in range(1, len(a)+1):
            for ids in combinations(range(len(a)), taille):
                self.assertGreaterEqual(determinant([[a[i][j] for j in ids] for i in ids]), 0)

    def test_identite_oblique_et_masse_couplee(self):
        k, m, b, z = self.donnees()
        self.assertEqual(produit(transpose(b), z), matrice([[0, 0], [0, 0]]))
        s, tc, _, _ = inverse_contrainte(k, m, b)
        reference = produit(produit(z, inverse(restreindre(k, z))), transpose(z))
        self.assertEqual(s, reference)
        self.assertEqual(tc, trace(produit(m, reference)))
        self.assertGreater(tc, 0)
        self.assertNotEqual(m[0][1], 0)
        self.verifier_psd(s)

    def test_minorant_sur_le_complement(self):
        k, m, b, z = self.donnees()
        _, tc, _, _ = inverse_contrainte(k, m, b)
        self.verifier_psd(difference(restreindre(k, z),
                                    echelle(restreindre(m, z), 1/tc)))
        # Une autre base du même espace de contraintes donne le même S.
        change = matrice([[2, 1], [1, 1]])
        sb, tb, _, _ = inverse_contrainte(k, m, produit(b, change))
        sa, ta, _, _ = inverse_contrainte(k, m, b)
        self.assertEqual((sb, tb), (sa, ta))

    def test_restriction_apres_perturbation_energetique(self):
        kq, m, b, z = self.donnees()
        eta = F(1, 10)
        ajout = gram(matrice([[1, -2, 1, 0]]))
        kd = difference(echelle(kq, 1-eta), echelle(ajout, -1))
        self.verifier_psd(difference(kd, echelle(kq, 1-eta)))
        _, tc, _, _ = inverse_contrainte(kq, m, b)
        borne = (1-eta)/tc
        self.verifier_psd(difference(restreindre(kd, z),
                                    echelle(restreindre(m, z), borne)))

    def test_cancellation_interdit_soustraction_binary64_naive(self):
        k = matrice([[F(1, 10**30), 0, 0, 0], [0, 1, 0, 0],
                     [0, 0, 2, 0], [0, 0, 0, 3]])
        b = matrice([[1], [0], [0], [0]])
        _, tc, tau, chi = inverse_contrainte(k, identite(4), b)
        self.assertEqual(tc, F(11, 6))
        self.assertEqual(1/tc, F(6, 11))
        self.assertEqual(float(tau)-float(chi), 0.)

    def test_complement_unidimensionnel_borne_exacte(self):
        k = matrice([[2, 1], [1, 3]])
        m = matrice([[4, 1], [1, 2]])
        b, z = matrice([[1], [2]]), matrice([[-2], [1]])
        _, tc, _, _ = inverse_contrainte(k, m, b)
        self.assertEqual(1/tc, restreindre(k, z)[0][0]/restreindre(m, z)[0][0])

    def test_contrainte_redondante_refusee(self):
        k, m, _, _ = self.donnees()
        with self.assertRaisesRegex(ValueError, "rang insuffisant"):
            inverse_contrainte(k, m, matrice([[1, 2], [2, 4], [3, 6], [4, 8]]))


if __name__ == "__main__":
    unittest.main()
