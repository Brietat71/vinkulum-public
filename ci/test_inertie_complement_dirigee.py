"""Contre-épreuves Fraction du certificat creux, sans modèle natif ni banc."""
from decimal import Decimal, Inexact, ROUND_HALF_DOWN, Rounded, localcontext
from fractions import Fraction as F
import unittest

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix

from inertie_complement_dirigee import (
    InertieImpossible, IntervallesSignes, MatriceIntervalles,
    certifier_inertie_complement, eliminer,
)
from vinkulum._ports.inverse_selectionnee import BudgetInverseSelectionnee
from test_inertie_complement import inertie_fraction, selle
from test_trace_complement import (
    difference, echelle, gram, inverse, produit, transpose,
)


def rationnel(a):
    return [[F(float(v)) for v in ligne] for ligne in np.asarray(a)]


def kkt_exact(d, m, b, gamma):
    return selle(difference(gram(rationnel(d)), echelle(rationnel(m), F(float(gamma)))),
                 rationnel(b))


class InertieComplementDirigee(unittest.TestCase):
    def contient(self, intervalle, exact):
        self.assertLessEqual(F(Decimal(intervalle[0])), exact)
        self.assertGreaterEqual(F(Decimal(intervalle[1])), exact)

    def preuve_exacte(self, a, preuve):
        """Rejoue les congruences rapportées avec une inverse Fraction générale."""
        self.assertEqual(tuple(preuve['signature']), inertie_fraction(a))
        actifs = list(range(len(a)))
        observe = [0, 0, 0]
        for p in preuve['pivots']:
            ids = p['indices']
            self.assertTrue(set(ids).issubset(actifs))
            bloc = [[a[i][j] for j in ids] for i in ids]
            sig = inertie_fraction(bloc)
            self.assertEqual(p['signature'], list(sig[:2]))
            self.assertEqual(sig[2], 0)
            observe = [u+v for u, v in zip(observe, sig)]
            if len(ids) == 1:
                self.contient(p['diagonal'], bloc[0][0])
            else:
                self.assertEqual(len(ids), 2)
                for cle, v in [('a', bloc[0][0]), ('b', bloc[0][1]), ('c', bloc[1][1])]:
                    self.contient(p[cle], v)
                self.contient(p['determinant'], bloc[0][0]*bloc[1][1]-bloc[0][1]**2)
            restant = [i for i in actifs if i not in ids]
            croise = [[a[i][j] for j in ids] for i in restant]
            correction = produit(produit(croise, inverse(bloc)), transpose(croise)) if restant else []
            nouveau = [[a[i][j]-correction[ii][jj] for jj, j in enumerate(restant)]
                       for ii, i in enumerate(restant)]
            for ii, i in enumerate(restant):
                for jj, j in enumerate(restant):
                    a[i][j] = nouveau[ii][jj]
            actifs = restant
        self.assertFalse(actifs)
        self.assertEqual(observe, preuve['signature'])

    def verifier_certificat(self, d, m, b, gamma, **options):
        resultat = certifier_inertie_complement(d, m, b, gamma, **options)
        self.assertTrue(resultat['certification_machine'])
        self.assertEqual(resultat['lambda_min'], gamma)
        self.assertEqual(resultat['dimension_complement'], len(m)-b.shape[1])
        self.assertEqual(resultat['preuve_masse']['signature'], [len(m), 0, 0])
        self.assertIn('aucune réponse certifiée', resultat['portee'])
        self.preuve_exacte(kkt_exact(d, m, b, gamma), resultat['preuve_kkt'])
        return resultat

    def donnees_obliques(self):
        d = np.array([[2., 1., 0., 1.], [0., 3., 1., 0.],
                      [0., 0., 2., 1.], [0., 0., 0., 4.],
                      [.125, -.25, .5, .125]])
        l = np.array([[1., 2., 0., 0.], [0., 2., 1., 0.],
                      [0., 0., 3., 1.], [0., 0., 0., 2.]])
        return d, l.T@l, np.array([[1., 0.], [0., 1.], [2., -1.], [1., 3.]])

    def test_produits_tous_quadrants_carres_et_contexte_ambiant(self):
        intervalles = [tuple(map(Decimal, x)) for x in [
            ('0', '0'), ('1.234567890123456789', '3.456789012345678901'),
            ('-3.456789012345678901', '-1.234567890123456789'),
            ('-2.345678901234567890', '3.456789012345678901'),
            ('0', '3.456789012345678901'), ('-3.456789012345678901', '0')]]
        with localcontext() as ctx:
            ctx.prec, ctx.rounding = 3, ROUND_HALF_DOWN
            ctx.traps[Inexact] = ctx.traps[Rounded] = True
            ar = IntervallesSignes(16, 100_000)
            for a in intervalles:
                aa = tuple(map(F, a))
                carre = ar.carre(a)
                minimum = F(0) if aa[0] <= 0 <= aa[1] else min(x*x for x in aa)
                maximum = max(x*x for x in aa)
                self.contient(carre, minimum)
                self.contient(carre, maximum)
                for b in intervalles:
                    exacts = [x*F(y) for x in aa for y in b]
                    observe = ar.mul(a, b)
                    self.contient(observe, min(exacts))
                    self.contient(observe, max(exacts))
            self.assertEqual(ctx.prec, 3)
            self.assertEqual(ctx.rounding, ROUND_HALF_DOWN)

    def test_signatures_et_gamma_face_aux_fractions(self):
        d, m, b = self.donnees_obliques()
        acceptes = refuses = 0
        for gamma in (.125, 1., 5., 100.):
            with self.subTest(gamma=gamma):
                attendu = inertie_fraction(kkt_exact(d, m, b, gamma))
                if attendu == (4, 2, 0):
                    acceptes += 1
                    resultat = self.verifier_certificat(d, m, b, gamma, precision=60)
                    self.preuve_exacte(rationnel(m), resultat['preuve_masse'])
                else:
                    refuses += 1
                    with self.assertRaisesRegex(InertieImpossible, 'non coercif') as cm:
                        certifier_inertie_complement(d, m, b, gamma, precision=60)
                    self.assertEqual(tuple(cm.exception.bilan['signature']), attendu)
                    self.assertEqual(cm.exception.bilan['signature_attendue'], [4, 2, 0])
        self.assertGreater(acceptes, 0)
        self.assertGreater(refuses, 0)

    def test_pivot_deux_indispensable_a_singulier_kkt_regulier(self):
        d, m, b = np.diag([1., 2., 3.]), np.eye(3), np.eye(3)[:, :1]
        self.assertEqual(inertie_fraction(difference(gram(rationnel(d)), rationnel(m))), (2, 0, 1))
        resultat = self.verifier_certificat(d, m, b, 1.)
        blocs = [p for p in resultat['preuve_kkt']['pivots'] if len(p['indices']) == 2]
        self.assertEqual(len(blocs), 1)
        self.assertEqual(blocs[0]['signature'], [1, 1])
        self.assertLess(F(Decimal(blocs[0]['determinant'][1])), 0)

    def test_masse_couplee_generale_de_taille_sept(self):
        n = 7
        l = np.eye(n)+np.diag(np.full(n-1, .25), 1)
        m = l.T@l
        d = np.diag(np.arange(3., n+3.))+np.diag(np.full(n-1, -.125), 1)
        b = np.array([[1., 0.], [0., 1.], [1., -1.], [2., 1.],
                      [-1., 2.], [3., -2.], [1., 1.]])
        resultat = self.verifier_certificat(d, m, b, .5, precision=60)
        self.assertEqual(resultat['preuve_masse']['methode'], 'congruences_dirigees')
        self.preuve_exacte(rationnel(m), resultat['preuve_masse'])

    def test_permutation_physique_et_reechelonnement_des_contraintes(self):
        d, m, b = self.donnees_obliques()
        reference = self.verifier_certificat(d, m, b, .125, precision=80)
        for exposants in ((-400, 400), (400, -400)):
            with self.subTest(exposants=exposants):
                bb = b*np.array([2.**e for e in exposants])
                resultat = self.verifier_certificat(d, m, bb, .125, precision=80,
                                                    permutation=np.array([2, 0, 3, 1]))
                self.assertEqual(resultat['permutation_physique'], [2, 0, 3, 1])
                self.assertNotEqual(resultat['contraintes_sha256'], reference['contraintes_sha256'])
                self.assertEqual(resultat['d_sha256'], reference['d_sha256'])
                self.assertEqual(resultat['masse_sha256'], reference['masse_sha256'])

    def test_petit_b_stocke_ne_peut_pas_etre_supprime(self):
        d, m = np.diag([2.**-50, 1., 2.]), np.eye(3)
        b = np.array([[2.**-1000], [0.], [0.]])
        resultat = self.verifier_certificat(d, m, b, .5, precision=80)
        self.assertEqual(resultat['rang_contraintes'], 1)
        with self.assertRaises(InertieImpossible) as cm:
            certifier_inertie_complement(d, m, np.zeros_like(b), .5, precision=80)
        self.assertIsInstance(cm.exception.bilan, dict)
        self.assertTrue(cm.exception.bilan)

    def test_rang_deficient_et_frontiere_singuliere_refuses(self):
        d, m = np.diag([1., 2., 3., 4.]), np.eye(4)
        redondant = np.array([[1., 2.], [0., 0.], [0., 0.], [0., 0.]])
        self.assertEqual(inertie_fraction(kkt_exact(d, m, redondant, .5)), (4, 1, 1))
        with self.assertRaisesRegex(InertieImpossible, 'zéro') as cm:
            certifier_inertie_complement(d, m, redondant, .5)
        self.assertGreaterEqual(cm.exception.bilan['pivots_effectues'], 1)
        b = np.eye(4)[:, :1]
        self.assertEqual(inertie_fraction(kkt_exact(d, m, b, 4.)), (3, 1, 1))
        with self.assertRaisesRegex(InertieImpossible, 'zéro') as cm:
            certifier_inertie_complement(d, m, b, 4.)
        self.assertIn('signature_partielle', cm.exception.bilan)

    def test_pivot_intervalle_ambigu_refuse_sans_signe_du_centre(self):
        ar = IntervallesSignes(32, 1000)
        matrice = MatriceIntervalles(2, 10)
        matrice.set(0, 0, (Decimal('-1e-30'), Decimal('2e-30')))
        matrice.set(1, 1, ar.un)
        with self.assertRaisesRegex(InertieImpossible, 'zéro') as cm:
            eliminer(matrice, ar)
        self.assertEqual(cm.exception.bilan['pivot'], 0)
        self.assertEqual(cm.exception.bilan['signature_partielle'], [0, 0, 0])
        self.assertEqual(cm.exception.bilan['pivots_effectues'], 0)

    def test_budget_coefficients_operations_et_bilans(self):
        d, m, b = self.donnees_obliques()
        resultat = self.verifier_certificat(d, m, b, .125)
        with self.assertRaisesRegex(InertieImpossible, 'coefficients') as cm:
            certifier_inertie_complement(d, m, b, .125, budget_coefficients=6)
        self.assertGreater(cm.exception.bilan['coefficients'], cm.exception.bilan['budget'])
        with self.assertRaisesRegex(BudgetInverseSelectionnee, 'opérations') as cm:
            certifier_inertie_complement(d, m, b, .125,
                                        budget_operations=resultat['operations_decimal']-1)
        self.assertGreater(cm.exception.bilan['operations_demandees'], cm.exception.bilan['budget_operations'])
        with self.assertRaisesRegex(InertieImpossible, 'dimension') as cm:
            certifier_inertie_complement(d, m, b, .125, budget_coefficients=1)
        self.assertGreater(cm.exception.bilan['dimension'], cm.exception.bilan['budget_coefficients'])
        with self.assertRaisesRegex(InertieImpossible, 'contraintes') as cm:
            certifier_inertie_complement(d, m, b, .125, budget_rectangulaire=1)
        self.assertGreater(cm.exception.bilan['coefficients_contraintes'],
                           cm.exception.bilan['budget_rectangulaire'])

    def test_original_d_refuse_le_faux_positif_du_gram_arrondi(self):
        a = 1.1
        d = np.array([[a, a+3*np.spacing(a)]])
        m, b, gamma = np.eye(2), np.ones((2, 1)), 2.**-54
        exact = kkt_exact(d, m, b, gamma)
        gram_arrondi = rationnel(d.T@d)
        faux = selle(difference(gram_arrondi, echelle(rationnel(m), F(gamma))), rationnel(b))
        self.assertEqual((F(d[0, 0])-F(d[0, 1]))**2, 9*F(2)**-104)
        self.assertEqual(gram_arrondi[0][0]+gram_arrondi[1][1]-2*gram_arrondi[0][1], F(2)**-52)
        self.assertEqual(inertie_fraction(exact), (1, 2, 0))
        self.assertEqual(inertie_fraction(faux), (2, 1, 0))
        with self.assertRaisesRegex(InertieImpossible, 'non coercif') as cm:
            certifier_inertie_complement(d, m, b, gamma, precision=80)
        self.assertEqual(cm.exception.bilan['signature'], [1, 2, 0])
        self.preuve_exacte(exact, cm.exception.bilan['preuve'])

    def test_certificat_independant_du_contexte_decimal_ambiant(self):
        d, m, b = self.donnees_obliques()
        reference = certifier_inertie_complement(d, m, b, .125, precision=40)
        with localcontext() as ctx:
            ctx.prec, ctx.rounding = 3, ROUND_HALF_DOWN
            ctx.traps[Inexact] = ctx.traps[Rounded] = True
            observe = certifier_inertie_complement(d, m, b, .125, precision=40)
        for cle in ('preuve_kkt', 'preuve_masse', 'operations_decimal', 'lambda_min'):
            self.assertEqual(observe[cle], reference[cle])

    def test_masses_invalides_sparse_copie_et_permutation_invalide(self):
        d, m, b = self.donnees_obliques()
        ds, ms = csr_matrix(d), csc_matrix(m)
        avant = [(x.data.copy(), x.indices.copy(), x.indptr.copy()) for x in (ds, ms)]
        resultat = certifier_inertie_complement(ds, ms, b, .125)
        self.assertTrue(resultat['certification_machine'])
        for x, (data, indices, indptr) in zip((ds, ms), avant):
            np.testing.assert_array_equal(x.data, data)
            np.testing.assert_array_equal(x.indices, indices)
            np.testing.assert_array_equal(x.indptr, indptr)
        for invalide in (np.diag([1., 1., 1., 0.]),
                         np.array([[1., 2., 0., 0.], [2., 1., 0., 0.],
                                   [0., 0., 1., 0.], [0., 0., 0., 1.]])):
            with self.subTest(masse=invalide.tolist()):
                with self.assertRaises(InertieImpossible) as cm:
                    certifier_inertie_complement(d, invalide, b, .125)
                self.assertTrue(cm.exception.bilan)
        for perm in ([0, 0, 2, 3], [0., 1., 2., 3.], [0, 1, 2]):
            with self.assertRaisesRegex(ValueError, 'permutation'):
                certifier_inertie_complement(d, m, b, .125, permutation=perm)


if __name__ == '__main__':
    unittest.main()
