"""Preuves/contre-preuves Fraction du chemin LDL à intervalles.

Petit témoin dense indépendant : aucun import du prototype creux, aucune
valeur propre flottante et aucun chronométrage. Les opérations exactes
d'intervalles sont arrondies vers l'extérieur à une précision Decimal
fixée pour éprouver l'invariant d'inclusion sur de petites matrices.
"""
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction as F
import math
import unittest

from test_inertie_complement import inertie_fraction, selle, permutation
from test_trace_complement import (matrice, produit, transpose, gram,
                                  difference, echelle, identite, inverse)


class RefusPreuve(ValueError):
    pass


class IntervallesTemoin:
    """Arithmétique d'essai ; Fraction reste l'oracle des extrémités."""
    def __init__(self, precision=16):
        self.bas = Context(prec=precision, rounding=ROUND_FLOOR)
        self.haut = Context(prec=precision, rounding=ROUND_CEILING)

    def arrondir(self, lo, hi):
        return (F(self.bas.divide(Decimal(lo.numerator), Decimal(lo.denominator))),
                F(self.haut.divide(Decimal(hi.numerator), Decimal(hi.denominator))))

    def point(self, x):
        return F(x), F(x)

    def add(self, x, y):
        return self.arrondir(x[0]+y[0], x[1]+y[1])

    def sub(self, x, y):
        return self.arrondir(x[0]-y[1], x[1]-y[0])

    def mul(self, x, y):
        candidats = [a*b for a in x for b in y]
        return self.arrondir(min(candidats), max(candidats))

    def square(self, x):
        lo = F(0) if x[0] <= 0 <= x[1] else min(x[0]**2, x[1]**2)
        return self.arrondir(lo, max(x[0]**2, x[1]**2))

    def div(self, x, y):
        if y[0] <= 0 <= y[1]:
            raise RefusPreuve("division par intervalle contenant zéro")
        candidats = [a/b for a in x for b in y]
        return self.arrondir(min(candidats), max(candidats))


def signature_pivot(p, ar):
    if len(p) == 1:
        lo, hi = p[0][0]
        if lo > 0:
            return (1, 0, 0), [[ar.div(ar.point(1), p[0][0])]]
        if hi < 0:
            return (0, 1, 0), [[ar.div(ar.point(1), p[0][0])]]
        raise RefusPreuve("pivot scalaire non séparé")
    a, b, c = p[0][0], p[0][1], p[1][1]
    det = ar.sub(ar.mul(a, c), ar.square(b))
    trace = ar.add(a, c)
    if det[1] < 0:
        sig = (1, 1, 0)
    elif det[0] > 0 and trace[0] > 0:
        sig = (2, 0, 0)
    elif det[0] > 0 and trace[1] < 0:
        sig = (0, 2, 0)
    else:
        raise RefusPreuve("pivot double non séparé")
    moins_b = ar.sub(ar.point(0), b)
    return sig, [[ar.div(c, det), ar.div(moins_b, det)],
                 [ar.div(moins_b, det), ar.div(a, det)]]


def etape(w, indices, ar):
    reste = [j for j in range(len(w)) if j not in indices]
    p = [[w[i][j] for j in indices] for i in indices]
    sig, pinv = signature_pivot(p, ar)
    s = [[ar.point(0) for _ in reste] for _ in reste]
    for ii, i in enumerate(reste):
        for jj in range(ii, len(reste)):
            j = reste[jj]
            v = w[i][j]
            for k, ik in enumerate(indices):
                for l, il in enumerate(indices):
                    v = ar.sub(v, ar.mul(ar.mul(w[i][ik], pinv[k][l]), w[il][j]))
            s[ii][jj] = s[jj][ii] = v
    return sig, s


def inertie_temoin(w, ar):
    signature = [0, 0, 0]
    while w:
        choix = None
        candidats = [[i] for i in range(len(w))] + [
            [i, j] for i in range(len(w)) for j in range(i+1, len(w))]
        for indices in candidats:
            try:
                choix = etape(w, indices, ar)
                break
            except RefusPreuve:
                continue
        if choix is None:
            raise RefusPreuve("aucun pivot séparé ; inertie inconnue")
        sig, w = choix
        signature = [x+y for x, y in zip(signature, sig)]
    return tuple(signature)


def assembler(d, m, b, gamma, ar):
    n, s = len(m), len(b[0])
    w = [[ar.point(0) for _ in range(n+s)] for _ in range(n+s)]
    for i in range(n):
        for j in range(i, n):
            v = ar.point(0)
            for ligne in d:
                v = ar.add(v, ar.mul(ar.point(ligne[i]), ar.point(ligne[j])))
            v = ar.sub(v, ar.mul(ar.point(gamma), ar.point(m[i][j])))
            w[i][j] = w[j][i] = v
        for j in range(s):
            w[i][n+j] = w[n+j][i] = ar.point(b[i][j])
    return w


class InertieDirigeePreuves(unittest.TestCase):
    def setUp(self):
        self.ar = IntervallesTemoin(16)

    def points(self, a):
        return [[self.ar.point(x) for x in ligne] for ligne in a]

    def contient(self, intervalles, exacte):
        self.assertEqual(len(intervalles), len(exacte))
        for ligne, valeurs in zip(intervalles, exacte):
            for (lo, hi), valeur in zip(ligne, valeurs):
                self.assertLessEqual(lo, valeur)
                self.assertLessEqual(valeur, hi)

    def test_signatures_pivots_un_et_deux(self):
        cas = [([[2]], (1, 0, 0)), ([[-2]], (0, 1, 0)),
               ([[0, 2], [2, 0]], (1, 1, 0)),
               ([[2, 1], [1, 2]], (2, 0, 0)),
               ([[-2, -1], [-1, -2]], (0, 2, 0)),
               ([[1, 2], [2, 1]], (1, 1, 0))]
        for p, attendu in cas:
            with self.subTest(p=p):
                sig, pinv = signature_pivot(self.points(p), self.ar)
                self.assertEqual(sig, attendu)
                self.assertEqual(sig, inertie_fraction(p))
                self.contient(pinv, inverse(matrice(p)))
        for p in ([[0]], [[1, 1], [1, 1]]):
            with self.assertRaises(RefusPreuve):
                signature_pivot(self.points(p), self.ar)

    def test_determinant_negatif_sans_diagonales_separees(self):
        p = [[(F(-1), F(1)), self.ar.point(2)],
             [self.ar.point(2), self.ar.point(0)]]
        self.assertEqual(signature_pivot(p, self.ar)[0], (1, 1, 0))
        for a in (F(-1), F(0), F(1)):
            self.assertEqual(inertie_fraction([[a, 2], [2, 0]]), (1, 1, 0))

    def test_trace_ou_determinant_seul_insuffisant(self):
        self.assertEqual(inertie_fraction([[2, 0], [0, -1]]), (1, 1, 0))
        self.assertEqual(inertie_fraction([[1, 0], [0, 1]]), (2, 0, 0))
        self.assertEqual(inertie_fraction([[-1, 0], [0, -1]]), (0, 2, 0))

    def test_diagonale_positive_ne_certifie_pas_une_masse_couplee(self):
        m = matrice([[1, 2], [2, 1]])
        self.assertTrue(all(m[i][i] > 0 for i in range(2)))
        self.assertEqual(inertie_fraction(m), (1, 1, 0))
        # Couper implicitement le couplage donne un autre problème SPD.
        self.assertEqual(inertie_fraction(identite(2)), (2, 0, 0))

    def test_carre_traversant_zero(self):
        x = (F(-2), F(3))
        self.assertEqual(self.ar.square(x), (F(0), F(9)))
        # min((-2)**2,3**2)=4 exclurait le carré de zéro.
        self.assertGreater(min(x[0]**2, x[1]**2), F(0))

    def test_intervalle_pivot_ne_se_reduit_pas_a_son_milieu(self):
        b = (F(0), F(3, 2))
        p = [[self.ar.point(1), b], [b, self.ar.point(1)]]
        with self.assertRaises(RefusPreuve):
            signature_pivot(p, self.ar)
        self.assertEqual(inertie_fraction([[1, F(3, 4)], [F(3, 4), 1]]), (2, 0, 0))
        self.assertEqual(inertie_fraction([[1, F(3, 2)], [F(3, 2), 1]]), (1, 1, 0))

    def test_dependance_perdue_elargit_sans_fausse_preuve(self):
        t = (F(1), F(2))
        p = [[t, t], [t, self.ar.add(t, self.ar.point(1))]]
        with self.assertRaises(RefusPreuve):
            signature_pivot(p, self.ar)
        # La famille corrélée est pourtant SPD : det([[t,t],[t,t+1]])=t.
        for t in (F(1), F(3, 2), F(2)):
            self.assertEqual(inertie_fraction([[t, t], [t, t+1]]), (2, 0, 0))

    def test_schur_double_enclos_et_congruence_exacte(self):
        a = matrice([[0, 2, F(1, 3), 2], [2, 3, 1, F(-2, 7)],
                     [F(1, 3), 1, 5, 1], [2, F(-2, 7), 1, 7]])
        p = [ligne[:2] for ligne in a[:2]]
        f = [ligne[:2] for ligne in a[2:]]
        q = [ligne[2:] for ligne in a[2:]]
        exact = difference(q, produit(produit(f, inverse(p)), transpose(f)))
        sig, schur = etape(self.points(a), [0, 1], self.ar)
        self.contient(schur, exact)
        self.assertEqual(inertie_fraction(a), tuple(x+y for x, y in zip(sig, inertie_fraction(exact))))

    def test_inverse_du_milieu_peut_fausser_la_signature(self):
        p = (F(1, 2), F(2))
        schur = self.ar.sub(self.ar.point(1), self.ar.div(self.ar.point(1), p))
        self.assertLess(schur[0], 0)
        self.assertGreater(schur[1], 0)
        faux_schur = 1-1/((p[0]+p[1])/2)
        self.assertGreater(faux_schur, 0)
        self.assertEqual(inertie_fraction([[F(1, 2), 1], [1, 1]]), (1, 1, 0))

    def test_assemblage_gram_arrondi_change_le_probleme(self):
        eps = F(1, 2**30)
        d = matrice([[1, 1, 0], [eps, 0, 0], [0, 0, 1]])
        b = matrice([[0], [0], [1]])
        exact = gram(d)
        faux = [[F(float(v)) for v in ligne] for ligne in exact]
        self.assertNotEqual(exact, faux)
        self.assertEqual(inertie_fraction(selle(exact, b)), (3, 1, 0))
        self.assertEqual(inertie_fraction(selle(faux, b)), (2, 1, 1))
        w = assembler(d, identite(3), b, F(0), self.ar)
        self.contient(w, selle(exact, b))
        # À cette précision la preuve peut refuser : ce n'est pas un
        # permis de remplacer le Gram exact par son arrondi ponctuel.
        with self.assertRaises(RefusPreuve):
            inertie_temoin(w, self.ar)
        ar = IntervallesTemoin(40)
        self.assertEqual(inertie_temoin(assembler(d, identite(3), b, F(0), ar), ar), (3, 1, 0))

    def test_binary64_exact_ne_signifie_pas_decimal_court(self):
        x, gamma = F(.1), F(.01)
        self.assertGreater(x*x-gamma, 0)
        self.assertEqual(Decimal(str(.1))**2-Decimal(str(.01)), Decimal(0))
        d = matrice([[x, 0], [0, 1]])
        b = matrice([[0], [1]])
        ar = IntervallesTemoin(40)
        self.assertEqual(inertie_temoin(assembler(d, identite(2), b, gamma, ar), ar), (2, 1, 0))

    def test_a_singuliere_kkt_regulier(self):
        d, m, b = matrice([[1, 0], [0, 2]]), identite(2), matrice([[1], [0]])
        a = difference(gram(d), m)
        self.assertEqual(inertie_fraction(a), (1, 0, 1))
        w = assembler(d, m, b, F(1), self.ar)
        sig, schur = etape(w, [0, 2], self.ar)
        self.assertEqual(sig, (1, 1, 0))
        self.assertEqual(inertie_temoin(schur, self.ar), (1, 0, 0))
        self.assertEqual(inertie_temoin(w, self.ar), (2, 1, 0))

    def test_acceptation_refus_et_frontiere(self):
        d = matrice([[F(1, 2**50), 0, 0, 0], [0, 1, 0, 0],
                     [0, 0, 2, 0], [0, 0, 0, 3]])
        b, m = matrice([[1], [0], [0], [0]]), identite(4)
        for gamma, attendu in ((F(1, 2), (4, 1, 0)), (F(3, 2), (3, 2, 0))):
            a = difference(gram(d), echelle(m, gamma))
            self.assertEqual(inertie_fraction(selle(a, b)), attendu)
            self.assertEqual(inertie_temoin(assembler(d, m, b, gamma, self.ar), self.ar), attendu)
        with self.assertRaises(RefusPreuve):
            inertie_temoin(assembler(d, m, b, F(1), self.ar), self.ar)
        self.assertEqual(inertie_fraction(selle(difference(gram(d), m), b)), (3, 1, 1))

    def test_masse_couplee_oblique_permutation(self):
        d = matrice([[2, 1, 0, 1], [0, 3, 1, 0], [0, 0, 2, 1], [0, 0, 0, 4]])
        m = gram(matrice([[1, 2, 0, 0], [0, 2, 1, 0], [0, 0, 3, 1], [0, 0, 0, 2]]))
        b = matrice([[1, 0], [0, 1], [2, -1], [1, 3]])
        for gamma in (F(1, 8), F(1), F(5)):
            exact = selle(difference(gram(d), echelle(m, gamma)), b)
            w = assembler(d, m, b, gamma, self.ar)
            self.contient(w, exact)
            attendu = inertie_fraction(exact)
            self.assertEqual(inertie_temoin(w, self.ar), attendu)
            ordre = [4, 2, 0, 5, 3, 1]
            self.assertEqual(inertie_temoin(permutation(w, ordre), self.ar), attendu)

    def test_rang_b_ne_decoule_pas_des_premiers_pivots(self):
        b = matrice([[1, 2], [0, 0], [0, 0]])
        c = selle(identite(3), b)
        self.assertEqual(inertie_fraction(c), (3, 1, 1))
        w = self.points(c)
        sig, reste = etape(w, [0], self.ar)
        self.assertEqual(sig, (1, 0, 0))
        with self.assertRaises(RefusPreuve):
            inertie_temoin(reste, self.ar)

    def test_supprimer_une_arete_de_remplissage_fausse_linertie(self):
        a = matrice([[1, 1, 1], [1, 2, 0], [1, 0, 2]])
        self.assertEqual(inertie_fraction(a), (2, 0, 1))
        _, s = etape(self.points(a), [0], self.ar)
        self.assertEqual(s, self.points([[1, -1], [-1, 1]]))
        # Omettre les nouvelles entrées hors diagonale ferait annoncer
        # trois pivots positifs à partir d'une matrice singulière.
        faux = [[s[0][0], self.ar.point(0)], [self.ar.point(0), s[1][1]]]
        self.assertEqual(inertie_temoin(faux, self.ar), (2, 0, 0))

    def test_pivot_multiplicateur_precoce_cree_une_clique(self):
        c = selle(identite(5), matrice([[1]]*5))
        _, s = etape(self.points(c), [0, 5], self.ar)
        self.assertEqual(s, self.points([[2 if i == j else 1 for j in range(4)] for i in range(4)]))
        self.assertEqual(sum(s[i][j] != self.ar.point(0) for i in range(4) for j in range(i+1, 4)), 6)

    def test_carre_binary64_ne_prouve_pas_la_bande_exacte(self):
        omega = math.nextafter(1., math.inf)
        gamma = omega*omega
        self.assertLess(F(gamma), F(omega)**2)
        self.assertEqual(F(omega)**2-F(gamma), F(1, 2**104))
        self.assertGreater(F(math.nextafter(gamma, math.inf)), F(omega)**2)


if __name__ == "__main__":
    unittest.main()
