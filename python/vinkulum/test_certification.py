"""Contre-épreuves des certificats : oracle rationnel et chemin mécanique."""
import copy
from fractions import Fraction as F
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from vinkulum import Noyau
from vinkulum import _verification_lineaire
from vinkulum.certification import (
    CertificationImpossible, certifier_initialisation, certifier_systeme, verifier_certificat,
)


def gauss_exact(a, b):
    """Élimination rationnelle indépendante de la construction par contraction."""
    n = len(b)
    m = [[F(float(v)) for v in row]+[F(float(bi))] for row, bi in zip(a, b, strict=True)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if m[i][j])
        m[j], m[pivot] = m[pivot], m[j]
        d = m[j][j]
        m[j] = [x/d for x in m[j]]
        for i in range(n):
            if i != j:
                d = m[i][j]
                m[i] = [u-d*v for u, v in zip(m[i], m[j], strict=True)]
    return [row[-1] for row in m]


def pendule():
    n = Noyau([0., -9.81, 0.])
    n.corps('pendule', 1., (.2*np.eye(3)).ravel().tolist(), [0., -1., 0.],
            v=[.3, 0., 0.], w=[0., 0., .3])
    n.liaison('pivot', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1])
    n.assemble()
    return n


def etat(n):
    try:
        reactions = n.reactions()
    except RuntimeError as e:
        reactions = ('non_initialisees', str(e))
    return n.etat_precis(), n.schema(), reactions, n.stats()


class CertificationLineaire(unittest.TestCase):
    def contient(self, c):
        q = verifier_certificat(c)
        a = [[float.fromhex(v) for v in l] for l in c['matrice']]
        b = list(map(float.fromhex, c['second_membre']))
        x = list(map(float.fromhex, c['solution']))
        reference = gauss_exact(a, b)
        erreurs = [abs(r-F(v)) for r, v in zip(reference, x, strict=True)]
        self.assertLessEqual(max(erreurs), q['borne_erreur_inf'])
        for e, borne in zip(erreurs, q['bornes_composantes'], strict=True):
            self.assertLessEqual(e, borne)
        return q

    def test_residu_minuscule_erreur_unitaire(self):
        a = np.diag([1., 2.**-60])
        b, x = [1., 2.**-60], [1., 0.]
        self.assertLess(np.max(np.abs(a@x-b)), 1e-12)
        c = certifier_systeme(a, b, x)
        self.assertEqual(self.contient(c)['borne_erreur_inf'], 1)
        with self.assertRaisesRegex(CertificationImpossible, 'supérieure'):
            certifier_systeme(a, b, x, erreur_max=1e-12)

    def test_systemes_denses_et_permutations(self):
        rng = np.random.default_rng(13002026)
        for n in (2, 4, 7):
            for _ in range(6):
                a = rng.integers(-3, 4, size=(n, n)).astype(float)+5*n*np.eye(n)
                a = a[rng.permutation(n)][:, rng.permutation(n)]
                b = rng.integers(-5, 6, size=n).astype(float)
                x = np.linalg.solve(a, b)
                x[0] = math.nextafter(x[0], math.inf)
                self.contient(certifier_systeme(a, b, x, erreur_max=1e-12))

    def test_matrice_indefinie_et_changements_echelle(self):
        a = np.array([[2., .3, 1.], [.3, 3., -1.], [1., -1., 0.]])
        for exposant in (-500, 0, 500):
            for p in ([0, 1, 2], [2, 0, 1]):
                aa = a[p][:, p]*2.**exposant
                b = np.array([1., 2., 0.])[p]*2.**exposant
                self.contient(certifier_systeme(aa, b, erreur_max=1e-12))

    def test_tres_mal_conditionne_peut_etre_certifie(self):
        for k in (10, 30, 50):
            a = [[1., 1.], [1., 1.+2.**-k]]
            self.contient(certifier_systeme(a, [1., 2.], x=[0., 0.]))

    def test_contraction_juste_sous_unite(self):
        petit = math.ulp(0.)
        # eta = 1 - 2^-1074 serait arrondi à 1 en binary64. L'inverse
        # exacte déborde, mais R=1 suffit ici à prouver unicité et erreur nulle.
        c = certifier_systeme([[petit]], [petit], [1.], inverse_approche=[[1.]])
        q = self.contient(c)
        self.assertEqual(q['eta_exact'], 1-F(petit))
        self.assertEqual(q['borne_erreur_inf'], 0)

    def test_singulier_et_fausse_inverse_refuses(self):
        with self.assertRaises(CertificationImpossible):
            certifier_systeme([[1., 1.], [1., 1.]], [0., 0.])
        for a, r in (([[1.]], [[0.]]), ([[0.]], [[1.]]), ([[1.]], [[2.]])):
            with self.assertRaisesRegex(CertificationImpossible, 'contraction'):
                certifier_systeme(a, [0.], [0.], inverse_approche=r)

    def test_alterations_et_frontiere_un_ulp(self):
        c = certifier_systeme([[3.]], [1.], [0.])
        self.contient(c)
        # La borne majorante de 1/3 doit dépasser le flottant arrondi au plus près.
        self.assertGreater(F(float.fromhex(c['borne_erreur_inf'])), F(1, 3))
        for cle, valeur in (
            ('borne_erreur_inf', math.nextafter(float.fromhex(c['borne_erreur_inf']), 0.).hex()),
            ('bornes_composantes', [0..hex()]),
            ('solution', [1..hex()]),
            ('second_membre', [2..hex()]),
            ('inverse_approche', [[0..hex()]]),
            ('matrice', [[0..hex()]]),
        ):
            modifie = copy.deepcopy(c)
            modifie[cle] = valeur
            with self.subTest(cle=cle), self.assertRaises(ValueError):
                verifier_certificat(modifie)

    def test_budget_dimensions_et_finitude(self):
        for kw in ({'dimension_max': 0}, {'dimension_max': 129}, {'dimension_max': 1.5},
                   {'erreur_max': -1}, {'erreur_max': math.inf}, {'erreur_max': math.nan}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                certifier_systeme([[1.]], [1.], **kw)
        for a, b, x in (([[math.nan]], [0.], [0.]), ([[1.]], [math.inf], [0.]),
                        ([[1.]], [0.], [math.nan]), ([[1., 0.]], [0.], [0.]),
                        ([[1.]], [0., 0.], [0.]), ([], [], [])):
            with self.assertRaises(ValueError):
                certifier_systeme(a, b, x)
        self.assertEqual(self.contient(certifier_systeme([[1.]], [0.], erreur_max=0))['borne_erreur_inf'], 0)

    def test_sous_normal_et_arrondi_de_borne(self):
        petit = math.ulp(0.)
        c = certifier_systeme([[2.]], [petit], [0.], inverse_approche=[[.5]])
        self.assertEqual(float.fromhex(c['borne_erreur_inf']), petit)
        self.contient(c)

    def test_verificateur_sans_site_packages(self):
        c = certifier_systeme([[3.]], [1.], [0.])
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'certificat.json'
            p.write_text(json.dumps(c))
            cmd = [sys.executable, '-I', '-S', _verification_lineaire.__file__, str(p)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            c['borne_erreur_inf'] = 0..hex()
            p.write_text(json.dumps(c))
            r = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(r.returncode, 1, r.stdout)


class CertificationInitialisation(unittest.TestCase):
    def test_solution_native_et_acceleration_du_premier_pas(self):
        n = pendule()
        n.enregistre_schema()
        avant = etat(n)
        with patch('numpy.linalg.solve', side_effect=AssertionError('solution externe interdite')):
            c = certifier_initialisation(n, erreur_max=1e-11)
        self.assertEqual(etat(n), avant)
        CertificationLineaire.contient(self, c)
        self.assertEqual(c['origine']['composantes_dynamiques'], 6)
        n.simule(.001, .001)
        accelerations, _, reactions = n.schema()[0]
        self.assertEqual(c['solution'], [v.hex() for v in accelerations+reactions])

    def test_chute_libre_exacte(self):
        n = Noyau([0., 0., -9.81])
        n.corps('libre', 1., np.eye(3).ravel().tolist(), [0., 0., 0.])
        c = certifier_initialisation(n, erreur_max=0.)
        self.assertEqual(list(map(float.fromhex, c['solution'])), [0., 0., -9.81, 0., 0., 0.])
        self.assertEqual(verifier_certificat(c)['borne_erreur_inf'], 0)

    def test_refus_sans_mutation(self):
        n = pendule()
        avant = etat(n)
        for kw in ({'erreur_max': 0.}, {'dimension_max': 5}, {'t': math.nan}):
            with self.assertRaises(ValueError):
                certifier_initialisation(n, **kw)
            self.assertEqual(etat(n), avant)
        n.liaison('double', None, 0, bloque_t=[0], bloque_r=[])
        avant = etat(n)
        with self.assertRaisesRegex(CertificationImpossible, 'redondantes'):
            certifier_initialisation(n)
        self.assertEqual(etat(n), avant)

    def test_contribution_modifiee_refusee(self):
        c = certifier_initialisation(pendule())
        c['contributions'][0][2] = 999..hex()
        with self.assertRaisesRegex(ValueError, 'contributions'):
            verifier_certificat(c)

    def test_contact_non_lisse_non_certifie(self):
        n = Noyau()
        n.corps('bille', 1., np.eye(3).ravel().tolist(), [0., 0., .5])
        n.contact('sol', 0, [0., 0., 0.], .5, nonlisse=True)
        avant = etat(n)
        with self.assertRaisesRegex(CertificationImpossible, 'contact non lisse'):
            certifier_initialisation(n)
        self.assertEqual(etat(n), avant)


if __name__ == '__main__':
    unittest.main()
