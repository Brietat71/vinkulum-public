"""Contre-épreuves Fraction des transferts de normes, sans solveur."""
from fractions import Fraction as F
import math
import unittest

from test_trace_complement import (
    matrice, produit, transpose, difference, echelle, identite, gram,
    determinant,
)
from test_inertie_complement import inertie_fraction

U = F(1, 2**53)


def gamma(n):
    return n*U/(1-n*U)


def absolu(a):
    return [[abs(x) for x in ligne] for ligne in a]


def somme(a, b):
    return [[x+y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def produit_fl(a, b):
    return [[F(sum(float(x)*float(y) for x, y in zip(ligne, col)))
             for col in zip(*b)] for ligne in a]


def difference_fl(a, b):
    return [[F(float(x)-float(y)) for x, y in zip(ra, rb)]
            for ra, rb in zip(a, b)]


def arrondi(a):
    return [[F(float(x)) for x in ligne] for ligne in a]


def norme_carree(a):
    return sum((x*x for ligne in a for x in ligne), F(0))


class ControleFacteursPreuves(unittest.TestCase):
    def exemple_gram(self):
        eps = F(1, 2**30)
        z = matrice([[1, 1], [0, eps]])
        return eps, z, gram(z), arrondi(gram(z))

    def test_gram_annule_une_petite_charge(self):
        eps, z, g, gh = self.exemple_gram()
        q = matrice([[1], [-1]])
        self.assertEqual(produit(z, q), matrice([[0], [-eps]]))
        self.assertEqual(produit(transpose(q), produit(gh, q))[0][0], 0)
        self.assertEqual(produit(transpose(q), produit(g, q))[0][0], eps**2)

    def test_gershgorin_sans_budget_sous_estime_meme_norme_operateur(self):
        eps, _, g, gh = self.exemple_gram()
        brut = max(sum(abs(x) for x in ligne) for ligne in gh)
        q = matrice([[1], [1]])
        rayleigh = produit(transpose(q), produit(g, q))[0][0]/2
        self.assertEqual(brut, 2)
        self.assertEqual(rayleigh, 2+eps**2/2)
        self.assertGreater(rayleigh, brut)
        budget = gamma(2)*sum(gh[i][i] for i in range(2))/(1-gamma(2))
        self.assertEqual(inertie_fraction(difference(echelle(identite(2), brut+budget), g)), (2, 0, 0))

    def test_gershgorin_pondere_tout_vecteur_strictement_positif(self):
        h = matrice([[1, 2], [2, 4]])
        bornes = []
        for valeurs in ((1, 1), (1, 2), (3, 7)):
            w = matrice([[v] for v in valeurs])
            hw = produit(h, w)
            beta = max(hw[i][0]/w[i][0] for i in range(2))
            bornes.append(beta)
            self.assertEqual(inertie_fraction(difference(echelle(identite(2), beta), h))[1], 0)
        self.assertEqual(bornes[:2], [6, 5])

    def test_poids_signe_et_rayleigh_ne_sont_pas_des_majorants(self):
        h = matrice([[1, 2], [2, 4]])
        w = matrice([[1], [-1]])
        hw = produit(h, w)
        faux = max(hw[i][0]/w[i][0] for i in range(2))
        self.assertEqual(faux, 2)
        self.assertEqual(inertie_fraction(difference(echelle(identite(2), faux), h)), (1, 1, 0))
        self.assertLess(F(1+4, 2), 4)  # Rayleigh de diag(1,4) en (1,1).

    def test_domination_symetrique_depuis_gram_asymetrique(self):
        eta = F(1, 100)
        g = matrice([[1, 1], [1, 1]])
        gh = matrice([[1, 1+eta], [1-eta, 1]])
        # d rationnel majore sqrt(1/(1-eta)).
        d = 1/(1-eta)
        h = [[min(abs(gh[i][j]), abs(gh[j][i]))+eta*d*d
              for j in range(2)] for i in range(2)]
        self.assertTrue(all(h[i][j] >= abs(g[i][j]) for i in range(2) for j in range(2)))
        beta = max(map(sum, h))
        self.assertEqual(inertie_fraction(difference(echelle(identite(2), beta), g)), (2, 0, 0))

    def test_norme_operateur_ne_donne_pas_un_minorant_par_charge(self):
        a, q = matrice([[1, 2]]), matrice([[2], [-1]])
        self.assertEqual(produit(a, q), matrice([[0]]))
        self.assertGreater(5*norme_carree(q), 0)

    def test_arrondi_du_produit_exterieur_augmente_son_rang(self):
        eps = F(1, 2**27)
        c, j = matrice([[1], [1+eps]]), matrice([[1, 1-eps]])
        cj, image = produit(c, j), produit_fl(c, j)
        y = matrice([[1-eps], [-1]])
        self.assertEqual(determinant(cj), 0)
        self.assertEqual(determinant(image), eps**2)
        self.assertEqual(produit(j, y), matrice([[0]]))
        self.assertEqual(produit(image, y), matrice([[0], [-eps**2]]))
        self.assertEqual(difference(image, cj), matrice([[0, 0], [0, eps**2]]))

    def test_residu_direct_de_reparation_sur_charges_signees(self):
        q = matrice([[F(3, 5)], [F(4, 5)]])
        l = matrice([[2, -1, F(1, 3)]])
        erreur = matrice([[F(1, 100), 0, -F(1, 50)], [0, F(1, 200), 0]])
        image = somme(produit(q, l), erreur)
        e = [sum(abs(erreur[i][j]) for i in range(2)) for j in range(3)]
        self.assertEqual(gram(q), identite(1))
        for valeurs in ((1, 2, 3), (1, -2, -3), (0, 1, 3), (-1, -1, -1)):
            y = matrice([[v] for v in valeurs])
            borne = abs(produit(l, y)[0][0])+sum(e[j]*abs(y[j][0]) for j in range(3))
            self.assertLessEqual(norme_carree(produit(image, y)), borne**2)

    def test_residu_qr_observe_nul_mais_exact_non_nul(self):
        eps = F(1, 2**27)
        q, r = 1+eps, 1-eps
        self.assertEqual(float(q)*float(r), 1.)
        self.assertEqual(1-q*r, eps**2)
        self.assertLessEqual(eps**2, gamma(1)*abs(q*r))

    def test_defaut_orthogonalite_observe_nul(self):
        eps = F(1, 2**30)
        q = matrice([[1], [eps]])
        self.assertEqual(produit_fl(transpose(q), q), identite(1))
        self.assertEqual(gram(q)[0][0]-1, eps**2)
        self.assertLessEqual(eps**2, gamma(2)*(1+eps**2))

    def test_facteur_kappa_ne_peut_pas_etre_omis(self):
        eps = F(1, 100)
        q = matrice([[1+eps, 0], [0, 1]])
        y = matrice([[1], [0]])
        zeta = 2*eps+eps**2
        self.assertGreater(norme_carree(produit(q, y)), norme_carree(y))
        self.assertEqual(norme_carree(produit(q, y)), (1+zeta)*norme_carree(y))

    def test_formation_du_petit_facteur_peut_annuler(self):
        eps = F(1, 2**27)
        rw, y = 1+eps, 1-eps
        self.assertEqual(1.-float(rw)*float(y), 0.)
        self.assertEqual(1-rw*y, eps**2)
        self.assertLessEqual(eps**2, gamma(1)*abs(rw*y))

    def test_borne_de_formation_preparee_par_colonnes(self):
        eps = F(1, 2**27)
        d = matrice([[1+eps, -1, 2], [2, 1-eps, -1], [1, 0, 1]])
        t = matrice([[1, 1-eps], [-1, 2], [eps, 1]])
        w = matrice([[1+eps, 1], [-1, 1-eps]])
        di = [ligne[:2] for ligne in d]
        ah, bh = produit_fl(d, t), produit_fl(di, w)
        gd, gi, gk = gamma(3), gamma(2), gamma(2)
        ck = gk+U*(1+gk)
        ap = produit(absolu(d), absolu(t))
        aip = produit(absolu(di), absolu(t[:2]))
        bp = produit(absolu(di), absolu(w))
        h0 = somme(echelle(ap, 2*gd), echelle(aip, U*(1+gd)))
        h1 = echelle(bp, gd*(1+ck)+ck+gi)
        for y in (matrice([[1-eps, 2], [1, -1]]), matrice([[0, -3], [2, 0]])):
            x = difference_fl(t[:2], produit_fl(w, y))+[t[2]]
            zh = produit_fl(d, x)
            ideal = difference(ah, produit(bh, y))
            erreur = absolu(difference(zh, ideal))
            borne = somme(h0, produit(h1, absolu(y)))
            self.assertTrue(all(erreur[i][j] <= borne[i][j] for i in range(3) for j in range(2)))

    def test_ymax_exact_ne_borne_pas_y_calcule(self):
        exact, calcule, erreur = F(1), F(101, 100), F(1, 100)
        self.assertGreater(abs(calcule), abs(exact))
        self.assertLessEqual(abs(calcule), abs(exact)+erreur)

    def test_norme_colonne_sous_estimee_par_sommation(self):
        petite = F(1, 2**27)
        valeurs = [F(1)]+[petite]*127
        accumule = 0.
        for x in valeurs:
            accumule += float(x)*float(x)
        observe = F(accumule)
        exact = sum(x*x for x in valeurs)
        self.assertEqual(observe, 1)
        self.assertEqual(exact, 1+127*petite**2)
        self.assertGreater(exact, observe)
        self.assertLessEqual(exact, observe/(1-gamma(len(valeurs))))

    def test_sous_flux_invalide_budget_purement_relatif(self):
        x = math.ldexp(1., -600)
        self.assertEqual(x*x, 0.)
        exact = F(x)**2
        self.assertGreater(exact, 0)
        self.assertEqual(gamma(1)*F(x*x), 0)

    def test_budget_normalise_mixte_sur_echelles_extremes(self):
        tau = math.ldexp(1., -1074)
        colonnes = ([1.]+[math.ldexp(1., -27)]*127,
                    [math.ldexp(1., 500), math.ldexp(1., 499)],
                    [math.ldexp(1., -600), math.ldexp(1., -601)],
                    [tau, tau, tau], [1., 3*tau, math.ldexp(1., -500)])
        for valeurs in colonnes:
            m = max(abs(x) for x in valeurs)
            accumule = 0.
            for x in valeurs:
                normalise = x/m
                accumule += normalise*normalise
            n = len(valeurs)
            borne = F(m)**2*(F(accumule)+(3*n+2)*F(tau))/(1-gamma(n+3))
            self.assertLessEqual(sum(F(x)**2 for x in valeurs), borne)


if __name__ == '__main__':
    unittest.main()
