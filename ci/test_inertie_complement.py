"""Identité d'inertie du complément, éprouvée en fractions exactes.

L'inertie est calculée par congruences et pivots 1×1 ou 2×2, sans valeurs
propres flottantes et sans supposer le complément inversible. Ce test
algébrique n'est ni un certificat machine à grande échelle ni un benchmark.
"""
from fractions import Fraction as F
import unittest

from test_trace_complement import (
    matrice, transpose, produit, difference, echelle, identite,
    inverse, determinant, gram, restreindre,
)


def inertie_fraction(a, pivots=None):
    """Renvoie (positives, négatives, nulles) par élimination congruente.

    Un pivot diagonal non nul donne un Schur exact. Si toute la diagonale
    restante est nulle, un coefficient hors diagonale b != 0 donne le pivot
    [[0,b],[b,0]], d'inertie (1,1,0). Une matrice restante nulle fournit sa
    nullité entière. Aucun seuil numérique ni extraction de racine.
    """
    a = matrice(a)
    n = len(a)
    if any(len(ligne) != n for ligne in a) or a != transpose(a):
        raise ValueError("matrice carrée symétrique requise")
    positifs = negatifs = 0
    while a:
        n = len(a)
        diagonal = next((i for i in range(n) if a[i][i] != 0), None)
        if diagonal is not None:
            ordre = [diagonal]+[i for i in range(n) if i != diagonal]
            a = [[a[i][j] for j in ordre] for i in ordre]
            pivot = a[0][0]
            positifs += int(pivot > 0)
            negatifs += int(pivot < 0)
            if pivots is not None:
                pivots.append(1)
            a = [[a[i][j]-a[i][0]*a[0][j]/pivot
                  for j in range(1, n)] for i in range(1, n)]
            continue
        paire = next(((i, j) for i in range(n) for j in range(i+1, n)
                      if a[i][j] != 0), None)
        if paire is None:
            return positifs, negatifs, n
        i, j = paire
        ordre = [i, j]+[k for k in range(n) if k not in paire]
        a = [[a[k][l] for l in ordre] for k in ordre]
        b = a[0][1]
        # P^-1=[[0,1/b],[1/b,0]] : soustraction exacte du Schur.
        a = [[a[k][l]-(a[k][0]*a[1][l]+a[k][1]*a[0][l])/b
              for l in range(2, n)] for k in range(2, n)]
        positifs += 1
        negatifs += 1
        if pivots is not None:
            pivots.append(2)
    return positifs, negatifs, 0


def selle(a, b):
    n, s = len(a), len(b[0])
    return [a[i]+b[i] for i in range(n)] + [
        [b[i][j] for i in range(n)]+[F(0)]*s for j in range(s)]


def permutation(a, ordre):
    return [[a[i][j] for j in ordre] for i in ordre]


def inertie_predite(a, b, z):
    """La formule simple nécessite réellement le rang colonne plein de B."""
    s = len(b[0])
    if determinant(gram(b)) == 0:
        raise ValueError("B redondant : rang colonne plein requis")
    if len(z[0]) != len(a)-s or determinant(gram(z)) == 0:
        raise ValueError("Z doit être une base complète du complément")
    if any(v != 0 for ligne in produit(transpose(b), z) for v in ligne):
        raise ValueError("Z doit appartenir à ker(B.T)")
    ip, im, iz = inertie_fraction(restreindre(a, z))
    return ip+s, im+s, iz


class InertieComplement(unittest.TestCase):
    def donnees_diagonales(self):
        k = matrice([[F(1, 10**30), 0, 0, 0], [0, 1, 0, 0],
                     [0, 0, 2, 0], [0, 0, 0, 3]])
        b = matrice([[1], [0], [0], [0]])
        z = matrice([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        return k, identite(4), b, z

    def donnees_obliques(self):
        k = gram(matrice([[2, 1, 0, 1], [0, 3, 1, 0],
                          [0, 0, 2, 1], [0, 0, 0, 4]]))
        m = gram(matrice([[1, 2, 0, 0], [0, 2, 1, 0],
                          [0, 0, 3, 1], [0, 0, 0, 2]]))
        b = matrice([[1, 0], [0, 1], [2, -1], [1, 3]])
        z = matrice([[-2, -1], [1, -3], [1, 0], [0, 1]])
        return k, m, b, z

    def test_elimination_independante_pivots_un_et_deux(self):
        pivots = []
        a = matrice([[0, 2, 0], [2, 0, 0], [0, 0, 0]])
        self.assertEqual(inertie_fraction(a, pivots), (1, 1, 1))
        self.assertEqual(pivots, [2])
        a = matrice([[2, 0, 1], [0, -3, 1], [1, 1, F(1, 6)]])
        # Le Schur des deux pivots diagonaux est exactement nul.
        pivots = []
        self.assertEqual(inertie_fraction(a, pivots), (1, 1, 1))
        self.assertEqual(pivots, [1, 1])
        self.assertEqual(inertie_fraction(matrice([[0, 0], [0, 0]])), (0, 0, 2))
        with self.assertRaisesRegex(ValueError, "symétrique"):
            inertie_fraction([[1, 2], [0, 1]])

    def test_gamma_un_demi_coercivite_acceptee(self):
        k, m, b, z = self.donnees_diagonales()
        a = difference(k, echelle(m, F(1, 2)))
        self.assertEqual(inertie_fraction(a), (3, 1, 0))
        self.assertEqual(inertie_fraction(restreindre(a, z)), (3, 0, 0))
        self.assertEqual(inertie_fraction(selle(a, b)), (4, 1, 0))
        self.assertEqual(inertie_fraction(selle(a, b)), inertie_predite(a, b, z))

    def test_gamma_trois_demis_coercivite_refusee(self):
        k, m, b, z = self.donnees_diagonales()
        a = difference(k, echelle(m, F(3, 2)))
        self.assertEqual(inertie_fraction(restreindre(a, z)), (2, 1, 0))
        self.assertEqual(inertie_fraction(selle(a, b)), (3, 2, 0))
        self.assertNotEqual(inertie_fraction(selle(a, b)), (4, 1, 0))
        self.assertEqual(inertie_fraction(selle(a, b)), inertie_predite(a, b, z))

    def test_gamma_un_frontiere_singuliere(self):
        k, m, b, z = self.donnees_diagonales()
        a = difference(k, m)
        self.assertEqual(determinant(restreindre(a, z)), 0)
        self.assertEqual(inertie_fraction(restreindre(a, z)), (2, 0, 1))
        self.assertEqual(inertie_fraction(selle(a, b)), (3, 1, 1))
        self.assertEqual(inertie_fraction(selle(a, b)), inertie_predite(a, b, z))

    def test_masse_couplee_contraintes_obliques_et_permutations(self):
        k, m, b, z = self.donnees_obliques()
        self.assertNotEqual(m[0][1], 0)
        for gamma in (F(0), F(1, 10), F(1), F(5), F(100)):
            with self.subTest(gamma=gamma):
                a = difference(k, echelle(m, gamma))
                attendu = inertie_predite(a, b, z)
                self.assertEqual(inertie_fraction(selle(a, b)), attendu)
                # Permutation physique cohérente de A, B et Z.
                ordre = [2, 0, 3, 1]
                ap = permutation(a, ordre)
                bp, zp = [b[i] for i in ordre], [z[i] for i in ordre]
                self.assertEqual(inertie_predite(ap, bp, zp), attendu)
                self.assertEqual(inertie_fraction(selle(ap, bp)), attendu)
                # Un ordre mêlant physique et multiplicateurs force une
                # autre élimination, sans changer la signature.
                cp = permutation(selle(a, b), [4, 2, 0, 5, 3, 1])
                self.assertEqual(inertie_fraction(cp), attendu)

    def test_congruence_explicite_sans_inverse_du_complement(self):
        k, m, b, z = self.donnees_obliques()
        y = produit(b, inverse(gram(b)))  # B.T Y = I.
        self.assertEqual(produit(transpose(b), y), identite(2))
        # Deux compléments réguliers et un complément identiquement nul.
        for a in (k, difference(k, echelle(m, F(5))), difference(m, m)):
            with self.subTest(inertie=inertie_fraction(restreindre(a, z))):
                yz = produit(transpose(y), produit(a, z))
                yy = restreindre(a, y)
                p = [z[i]+y[i]+[F(0), F(0)] for i in range(4)] + [
                    [-v for v in yz[j]]+[-v/2 for v in yy[j]]+identite(2)[j]
                    for j in range(2)]
                self.assertNotEqual(determinant(p), 0)
                h = restreindre(a, z)
                cible = [ligne+[F(0)]*4 for ligne in h] + [
                    [F(0)]*4+ligne for ligne in identite(2)] + [
                    [F(0)]*2+ligne+[F(0)]*2 for ligne in identite(2)]
                self.assertEqual(restreindre(selle(a, b), p), cible)
                self.assertEqual(inertie_fraction(selle(a, b)), inertie_predite(a, b, z))
        # La nullité du complément reste visible dans l'inertie entière.
        self.assertEqual(inertie_fraction(selle(difference(m, m), b)), (2, 2, 2))

    def test_changement_de_base_des_contraintes(self):
        k, m, b, z = self.donnees_obliques()
        a = difference(k, echelle(m, F(3, 2)))
        change = matrice([[2, 1], [1, 1]])
        bc = produit(b, change)
        self.assertEqual(inertie_fraction(selle(a, b)), inertie_fraction(selle(a, bc)))
        self.assertEqual(inertie_predite(a, b, z), inertie_predite(a, bc, z))

    def test_b_redondant_invalide_identite_simple(self):
        k, m, b, z = self.donnees_diagonales()
        a = difference(k, echelle(m, F(1, 2)))
        redondant = [ligne+[2*ligne[0]] for ligne in b]
        self.assertEqual(determinant(gram(redondant)), 0)
        with self.assertRaisesRegex(ValueError, "B redondant"):
            inertie_predite(a, redondant, z)
        observe = inertie_fraction(selle(a, redondant))
        self.assertEqual(observe, (4, 1, 1))
        # Rang réel 1 : une direction de multiplicateur est nulle.
        self.assertEqual(observe, (3+1, 0+1, 0+1))
        # Ajouter (s,s,0) avec s=2 donnerait même une mauvaise dimension.
        self.assertNotEqual(observe, (3+2, 0+2, 0))

    def test_tres_petit_mode_n_exige_aucune_soustraction_de_traces(self):
        k, m, b, z = self.donnees_diagonales()
        tau = sum((F(1)/k[i][i] for i in range(4)), F(0))
        retire = F(1)/k[0][0]
        self.assertEqual(tau-retire, F(11, 6))
        # Cette comparaison illustre la cancellation, pas un calcul
        # d'inertie flottant : l'inertie ci-dessous utilise Fraction seule.
        self.assertEqual(float(tau)-float(retire), 0.)
        a = difference(k, echelle(m, F(1, 2)))
        self.assertEqual(inertie_fraction(selle(a, b)), (4, 1, 0))
        self.assertEqual(inertie_predite(a, b, z), (4, 1, 0))

    def test_decision_possible_au_dela_du_minorant_par_trace(self):
        k, m, b, z = self.donnees_diagonales()
        trace_complement = F(1)+F(1, 2)+F(1, 3)
        gamma = F(3, 4)
        self.assertGreater(gamma, 1/trace_complement)
        a = difference(k, echelle(m, gamma))
        self.assertEqual(inertie_fraction(selle(a, b)), (4, 1, 0))
        self.assertEqual(inertie_predite(a, b, z), (4, 1, 0))


if __name__ == "__main__":
    unittest.main()
