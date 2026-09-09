"""Contre-épreuves rationnelles d'une enveloppe Krylov contrainte uniforme.

Les calculs restent de petite dimension. La preuve uniforme utilise
résolvante et convexité de Bernstein ; les points de contrôle complètent
ces identités et ne remplacent pas la preuve de couverture de la bande.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, isqrt
import unittest

from test_trace_complement import (
    matrice, transpose, produit, difference, echelle, identite,
    inverse, determinant, gram, restreindre, inverse_contrainte,
)


def somme(a, b):
    return difference(a, echelle(b, -1))


def zeros(n, p):
    return [[F(0) for _ in range(p)] for _ in range(n)]


def racine_superieure(q, bits=32):
    """Racine majorée par un rationnel dyadique, sans flottants."""
    if q < 0:
        raise ValueError("racine d'un rationnel négatif")
    echelle2 = 1 << (2*bits)
    k = isqrt(q.numerator*echelle2//q.denominator)
    if k*k*q.denominator < q.numerator*echelle2:
        k += 1
    return F(k, 1 << bits)


def norme_superieure(a, metrique):
    """La norme de Frobenius pondérée majore la norme d'opérateur."""
    g = restreindre(metrique, a)
    return racine_superieure(sum((g[i][i] for i in range(len(g))), F(0)))


def polynome(coefficients, x):
    valeur = zeros(len(coefficients[0]), len(coefficients[0][0]))
    for a in reversed(coefficients):
        valeur = somme(echelle(valeur, x), a)
    return valeur


def bernstein(coefficients, rho):
    """Contrôles sur [0,rho], pour un polynôme matriciel en puissances."""
    degre = len(coefficients)-1
    controles = []
    for k in range(degre+1):
        b = zeros(len(coefficients[0]), len(coefficients[0][0]))
        for j in range(k+1):
            b = somme(b, echelle(coefficients[j], F(comb(k, j), comb(degre, j))*rho**j))
        controles.append(b)
    return controles


def evaluer_bernstein(controles, x):
    degre = len(controles)-1
    valeur = zeros(len(controles[0]), len(controles[0][0]))
    for k, b in enumerate(controles):
        valeur = somme(valeur, echelle(b, comb(degre, k)*x**k*(1-x)**(degre-k)))
    return valeur


class KrylovContraintPreuves(unittest.TestCase):
    def assert_psd(self, a):
        self.assertEqual(a, transpose(a))
        for n in range(1, len(a)+1):
            for ids in combinations(range(len(a)), n):
                self.assertGreaterEqual(determinant([[a[i][j] for j in ids] for i in ids]), 0)

    def donnees(self):
        d = matrice([[1, 0, 0, 1], [0, 2, 1, 1],
                     [0, 0, 3, 1], [0, 0, 0, 2]])
        mf = gram(matrice([[1, 1, 0, 0], [0, 1, 0, 1],
                           [0, 0, 1, 1], [0, 0, 0, 1]]))
        kf = gram(d)
        k, m = [r[:3] for r in kf[:3]], [r[:3] for r in mf[:3]]
        b = matrice([[1], [0], [0]])
        z = matrice([[0, 0], [1, 0], [0, 1]])
        # Relèvement arbitraire et direction retenue non modale.
        t = matrice([[F(1, 2), 1], [F(-1, 4), 1], [F(1, 3), 0], [2, 0]])
        v = matrice([[0], [F(1, 2)], [0]])
        self.assertEqual(restreindre(k, v), identite(1))
        lam, rho = F(1), F(1, 2)
        self.assert_psd(difference(restreindre(k, z), echelle(restreindre(m, z), lam)))
        g0, g1 = produit(kf[:3], t), produit(mf[:3], t)
        theta = echelle(restreindre(m, v), lam)
        a0 = produit(transpose(v), g0)
        a1 = echelle(produit(transpose(v), g1), -lam)
        return dict(d=d, kf=kf, mf=mf, k=k, m=m, b=b, z=z, t=t, v=v,
                    lam=lam, rho=rho, g0=g0, g1=g1, theta=theta, a0=a0, a1=a1)

    def residus(self, a, v=None):
        v = a["v"] if v is None else v
        p, q = produit(a["k"], v), echelle(produit(a["m"], v), a["lam"])
        e0 = difference(a["g0"], produit(p, a["a0"]))
        e1 = difference(echelle(a["g1"], -a["lam"]), produit(p, a["a1"]))
        f = difference(q, produit(p, a["theta"]))
        c1 = somme(e1, produit(f, a["a0"]))
        dd = somme(produit(a["theta"], a["a0"]), a["a1"])
        return e0, e1, f, c1, dd

    def coefficients(self, a, mu):
        res = inverse(difference(identite(len(a["theta"])), echelle(a["theta"], mu)))
        return produit(res, somme(a["a0"], echelle(a["a1"], mu)))

    def residu_direct(self, a, mu, v=None):
        v = a["v"] if v is None else v
        az = difference(a["k"], echelle(a["m"], a["lam"]*mu))
        g = difference(a["g0"], echelle(a["g1"], a["lam"]*mu))
        return difference(g, produit(az, produit(v, self.coefficients(a, mu))))

    def test_identite_residuelle_signes_et_couplages(self):
        a = self.donnees()
        e0, e1, f, c1, dd = self.residus(a)
        self.assertNotEqual(a["mf"][1][3], 0)
        self.assertNotEqual(produit(a["m"], matrice([[1], [1], [0]])), a["b"])
        for mu in (F(0), F(1, 8), F(1, 4), F(1, 2)):
            res = inverse(difference(identite(1), echelle(a["theta"], mu)))
            forme1 = somme(somme(e0, echelle(e1, mu)),
                           echelle(produit(f, self.coefficients(a, mu)), mu))
            forme2 = somme(somme(e0, echelle(c1, mu)),
                           echelle(produit(f, produit(res, dd)), mu**2))
            self.assertEqual(self.residu_direct(a, mu), forme1)
            self.assertEqual(forme1, forme2)

    def test_enveloppe_uniforme_polynome_et_queue(self):
        a = self.donnees()
        sm, _, _, _ = inverse_contrainte(a["m"], a["m"], a["b"])
        e0, _, f, c1, dd = self.residus(a)
        rho, t = a["rho"], a["theta"][0][0]
        self.assertLess(rho*t, 1)
        for q in (2, 4, 6):
            with self.subTest(q=q):
                poly = [e0, c1]
                puissance = identite(1)
                for j in range(2, q):
                    poly.append(produit(f, produit(puissance, dd)))
                    puissance = produit(puissance, a["theta"])
                controles = bernstein(poly, rho)
                borne_poly = max(norme_superieure(c, sm) for c in controles)
                queue = (rho**q*norme_superieure(produit(f, puissance), sm)
                         *norme_superieure(dd, identite(1))/(1-rho*t))
                delta = borne_poly+queue
                # Chaque contrôle satisfait la majoration matricielle
                # requise à la preuve uniforme par convexité.
                for c in controles:
                    self.assert_psd(difference(echelle(identite(2), borne_poly**2),
                                               restreindre(sm, c)))
                for j in range(17):
                    mu = rho*F(j, 16)
                    res = inverse(difference(identite(1), echelle(a["theta"], mu)))
                    reste = echelle(produit(f, produit(puissance, produit(res, dd))), mu**q)
                    self.assertEqual(self.residu_direct(a, mu), somme(polynome(poly, mu), reste))
                    self.assertEqual(polynome(poly, mu), evaluer_bernstein(controles, mu/rho))
                    self.assert_psd(difference(echelle(identite(2), delta**2),
                                               restreindre(sm, self.residu_direct(a, mu))))

    def test_certificat_polynomial_sur_toute_la_bande(self):
        a = self.donnees()
        sm, _, _, _ = inverse_contrainte(a["m"], a["m"], a["b"])
        e0, e1, f, c1, _ = self.residus(a)
        h, rho = a["theta"][0][0], a["rho"]
        # Ici la résolvante est scalaire. N=(1-mu*h)R est polynomial :
        numerateur = [e0, difference(c1, echelle(e0, h)),
                      difference(produit(f, a["a1"]), echelle(e1, h))]
        gram_poly = [zeros(2, 2) for _ in range(5)]
        for i, ni in enumerate(numerateur):
            for j, nj in enumerate(numerateur):
                gram_poly[i+j] = somme(gram_poly[i+j], produit(transpose(ni), produit(sm, nj)))
        denom = [identite(2), echelle(identite(2), -2*h), echelle(identite(2), h*h),
                 zeros(2, 2), zeros(2, 2)]
        gd, gg = bernstein(denom, rho), bernstein(gram_poly, rho)
        # Chaque coefficient de denom est un multiple positif de I.
        # Frobenius borne tous les signes des coefficients de N.T S_M N.
        delta2 = max(norme_superieure(g, identite(2))/d[0][0] for d, g in zip(gd, gg))
        for d, g in zip(gd, gg):
            self.assertGreater(d[0][0], 0)
            self.assert_psd(difference(echelle(d, delta2), g))
        # Par convexité, delta²*(1-mu*h)² I - N.T S_M N est PSD
        # sur TOUT [0,rho], sans échantillonnage pour ce certificat.
        for mu in (F(0), rho/2, rho):
            self.assertEqual(polynome(numerateur, mu),
                             echelle(self.residu_direct(a, mu), 1-mu*h))
            self.assert_psd(difference(echelle(identite(2), delta2),
                                       restreindre(sm, self.residu_direct(a, mu))))

    def test_bernstein_conserve_compensations_polynomiales(self):
        poly = [matrice([[F(1, 4)]]), matrice([[-1]]), matrice([[1]])]
        c = bernstein(poly, F(1))
        self.assertEqual(c, [matrice([[F(1, 4)]]), matrice([[F(-1, 4)]]), matrice([[F(1, 4)]])])
        triangle = sum(abs(a[0][0]) for a in poly)
        borne = max(abs(a[0][0]) for a in c)
        self.assertEqual(triangle, F(9, 4))
        self.assertEqual(borne, F(1, 4))

    def test_racine_duale_au_lieu_dune_difference_de_gram(self):
        eps = F(1, 2**30)
        r, b = matrice([[1], [eps]]), matrice([[1], [0]])
        sm, _, _, _ = inverse_contrainte(identite(2), identite(2), b)
        reference = restreindre(sm, r)[0][0]
        self.assertEqual(reference, eps**2)
        # Orthogonalisation/Householder : les coordonnées restantes sont
        # ici directement la seconde ligne, sans différence de grandes normes.
        restant = matrice([[eps]])
        self.assertEqual(gram(restant)[0][0], reference)
        total, reaction = gram(r)[0][0], F(1)
        self.assertEqual(float(total)-float(reaction), 0.)
        self.assertGreater(float(gram(restant)[0][0]), 0.)

    def test_reparation_contrainte_et_defaut_energetique(self):
        a = self.donnees()
        mu = F(1, 3)
        delta = F(1, 100)
        c0 = matrice([[2], [0], [0]])
        c = produit(c0, inverse(produit(transpose(a["b"]), c0)))
        self.assertEqual(produit(transpose(a["b"]), c), identite(1))
        brut = somme(a["v"], echelle(c, delta))
        h = produit(c, produit(transpose(a["b"]), brut))
        repare = difference(brut, h)
        self.assertEqual(repare, a["v"])
        self.assertEqual(produit(transpose(a["b"]), repare), zeros(1, 1))
        y = self.coefficients(a, mu)
        xb = difference(a["t"], produit(brut+[ [F(0)] ], y))
        xr = difference(a["t"], produit(repare+[ [F(0)] ], y))
        correction = produit(h+[ [F(0)] ], y)
        self.assertEqual(somme(xb, correction), xr)
        az = difference(a["k"], echelle(a["m"], mu*a["lam"]))
        rr, rb = self.residu_direct(a, mu, repare), self.residu_direct(a, mu, brut)
        self.assertEqual(rr, somme(rb, produit(az, produit(h, y))))
        # L'expansion garde ses signes même si les coefficients n'ont
        # pas été recalculés par Galerkin sur la base réparée.
        e0, _, f, c1, dd = self.residus(a, repare)
        res = inverse(difference(identite(1), echelle(a["theta"], mu)))
        self.assertEqual(rr, somme(somme(e0, echelle(c1, mu)),
                                   echelle(produit(f, produit(res, dd)), mu**2)))
        rz = produit(produit(a["z"], inverse(restreindre(az, a["z"]))), transpose(a["z"]))
        g = difference(a["g0"], echelle(a["g1"], mu*a["lam"]))
        xe = difference(a["t"], produit(rz, g)+[ [F(0), F(0)] ])
        af = difference(a["kf"], echelle(a["mf"], mu*a["lam"]))
        exact = restreindre(af, xe)
        self.assertEqual(difference(restreindre(af, xr), exact), restreindre(rz, rr))
        self.assertNotEqual(difference(restreindre(af, xb), exact), restreindre(rz, rb))
        sm, _, _, _ = inverse_contrainte(a["m"], a["m"], a["b"])
        eta = norme_superieure(rr, sm)
        alpha = a["lam"]*(1-mu)
        for metrique, facteur in ((a["mf"], F(1)), (a["kf"], a["lam"])):
            bborne = racine_superieure(facteur)*eta/alpha+norme_superieure(correction, metrique)
            self.assert_psd(difference(echelle(identite(2), bborne**2),
                                       restreindre(metrique, difference(xb, xe))))

    def test_modele_q_nefface_pas_defaut_d(self):
        kq = matrice([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        kd = matrice([[1, 0, 0], [0, 3, 0], [0, 0, 3]])
        b, charge = matrice([[1], [0], [0]]), matrice([[0], [1], [0]])
        candidat = produit(inverse(kq), charge)
        rq = difference(charge, produit(kq, candidat))
        rd = difference(charge, produit(kd, candidat))
        self.assertEqual(rq, zeros(3, 1))
        self.assertNotEqual(rd, rq)
        sm, _, _, _ = inverse_contrainte(identite(3), identite(3), b)
        erreur = difference(candidat, produit(inverse(kd), charge))
        self.assertEqual(gram(erreur)[0][0], F(1, 36))
        self.assertEqual(restreindre(sm, rd)[0][0]/F(3)**2, F(1, 36))

    def test_marge_du_bloc_mode_et_ports_entier(self):
        # Le port seul ne détecte pas la résonance de la coordonnée retenue.
        eps = F(1, 100)
        shat = matrice([[1, 0], [0, eps]])
        sexact = matrice([[1, 0], [0, 0]])
        self.assertEqual(determinant(sexact), 0)
        erreur_schur = eps
        self.assertGreater(shat[0][0]-erreur_schur, 0)
        sigma_bloc = min(shat[0][0], shat[1][1])
        self.assertEqual(sigma_bloc-erreur_schur, 0)

    def test_action_anisotrope_du_gap_ne_vaut_pas_eta_carre(self):
        eps = F(1, 10**6)
        r = matrice([[eps, 1]])
        gap = gram(r)  # Complément scalaire A_N=1, donc alpha=1.
        q = matrice([[1], [0]])
        eta2 = gram(produit(r, q))[0][0]
        delta_op2 = 1+eps**2
        action2 = gram(produit(gap, q))[0][0]
        self.assertEqual(eta2, eps**2)
        self.assertEqual(action2, eps**2*(1+eps**2))
        # ||gap q||=eps*sqrt(1+eps²), exactement égal à delta_op*eta.
        self.assertEqual(action2, delta_op2*eta2)
        self.assertEqual(restreindre(gap, q)[0][0], eta2)
        # L'énergie q.T gap q=eta² n'est pas la norme de l'action.
        self.assertGreater(action2, eta2**2)
        # Le majorant isotrope delta_op²||q|| est ici bien plus grand.
        isotrope2 = delta_op2**2*gram(q)[0][0]
        self.assertGreater(isotrope2, 10**12*action2)

    def test_majorants_absolus_ne_donnent_pas_loewner_rang_un(self):
        r = matrice([[1, -1]])
        gap = gram(r)
        delta = matrice([[1], [1]])
        majorant_absolu = produit(delta, transpose(delta))
        for i in range(2):
            for j in range(2):
                self.assertLessEqual(abs(gap[i][j]), majorant_absolu[i][j])
        q = matrice([[1], [-1]])
        self.assertEqual(restreindre(gap, q)[0][0], 4)
        self.assertEqual(restreindre(majorant_absolu, q)[0][0], 0)
        ecart = difference(majorant_absolu, gap)
        self.assertLess(restreindre(ecart, q)[0][0], 0)
        self.assertLess(determinant(ecart), 0)
        # À l'inverse, les seules normes de colonnes perdent cette
        # compensation exacte entre charges : R*(1,1)=0, mais eta_col=2.
        compensation = matrice([[1], [1]])
        self.assertEqual(produit(r, compensation), matrice([[0]]))
        self.assertEqual(sum(abs(v[0]) for v in compensation), 2)

    def test_action_anisotrope_sur_charges_signees(self):
        r = matrice([[3, 0, -3], [4, 5, 4]])
        alpha = F(3, 2)
        an = matrice([[alpha, 0], [0, 2*alpha]])
        ainv = inverse(an)
        gap = produit(transpose(r), produit(ainv, r))
        gr = gram(r)
        delta_op2 = F(81)
        delta_colonnes = [F(5)]*3
        # Ces prémisses exactes établissent la coercivité et les normes
        # requises ; les charges ne servent pas à les estimer.
        self.assert_psd(difference(echelle(identite(2), 1/alpha), ainv))
        self.assert_psd(difference(echelle(identite(3), delta_op2), gr))
        for j in range(3):
            self.assertEqual(gr[j][j], delta_colonnes[j]**2)
        for valeurs in ((1, 0, 0), (-1, 0, 0), (1, -1, 2),
                        (-2, 3, -1), (0, 0, 0), (1, F(-8, 5), 1)):
            with self.subTest(charge=valeurs):
                q = matrice([[v] for v in valeurs])
                colonne = sum((d*abs(v) for d, v in zip(delta_colonnes, valeurs)), F(0))
                # eta=min(delta_op*||q||, sum(delta_j*|q_j|)).
                # Son carré est rationnel : aucune racine n'est évaluée.
                eta2 = min(delta_op2*gram(q)[0][0], colonne**2)
                residu2 = gram(produit(r, q))[0][0]
                action2 = gram(produit(gap, q))[0][0]
                self.assertLessEqual(residu2, eta2)
                self.assertLessEqual(restreindre(gap, q)[0][0], eta2/alpha)
                self.assertLessEqual(action2, delta_op2*eta2/alpha**2)

    def test_racine_rationnelle_est_majoree(self):
        for q in (F(0), F(1, 3), F(2), F(17, 19), F(1, 2**80)):
            b = racine_superieure(q)
            self.assertGreaterEqual(b*b, q)
            if b:
                self.assertLess((b-F(1, 2**32))**2, q)


if __name__ == "__main__":
    unittest.main()
