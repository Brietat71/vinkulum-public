"""Quotients exacts : rang, compatibilité, jauge et réaction généralisée."""
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

from vinkulum import Noyau, _verification_lineaire
from vinkulum.certification import (
    CertificationImpossible, _certifie_quotient, certifier_initialisation,
    certifier_systeme, verifier_certificat,
)
from vinkulum.test_certification import etat, gauss_exact, pendule


def construit(g, c, actives, masse=None, force=None, perturbe=False):
    g = np.array(g, dtype=float)
    m, n = g.shape
    masse = np.eye(n) if masse is None else np.array(masse, dtype=float)
    force = np.arange(1., n+1.) if force is None else np.array(force, dtype=float)
    a = np.zeros((n+m, n+m))
    a[:n, :n] = masse
    b = np.r_[force, np.where(actives, c, 0.)]
    for i, actif in enumerate(actives):
        if actif:
            a[n+i, :n] = g[i]
            a[:n, n+i] = g[i]
        else:
            a[n+i, n+i] = 1.
    x = np.linalg.solve(a, b)
    if perturbe:
        x = np.array([math.nextafter(v, math.inf) for v in x])
    s = certifier_systeme(a, b, x)
    trip = [(i, j, float(g[i, j])) for i in range(m) for j in range(n)]
    return _certifie_quotient(s, g.tolist(), list(map(float, c)), actives, trip)


def exemple(perturbe=False):
    return construit([[3., 0., 0.], [0., 2., 0.], [1., 0., 0.], [1., 2., 0.]],
                     [0.]*4, [True, True, False, False],
                     masse=[[2., .25, 0.], [.25, 1., .5], [0., .5, 3.]],
                     perturbe=perturbe)


class QuotientExact(unittest.TestCase):
    def test_reference_rationnelle_et_reaction(self):
        d = exemple(perturbe=True)
        q = verifier_certificat(d)
        self.assertEqual(q['rang_contraintes'], 2)
        self.assertFalse(q['multiplicateurs_uniques'])
        self.assertTrue(q['acceleration_unique'])
        s = d['systeme']
        a = [[float.fromhex(v) for v in l] for l in s['matrice']]
        b = list(map(float.fromhex, s['second_membre']))
        x = list(map(float.fromhex, s['solution']))
        ref = gauss_exact(a, b)
        for i, h in enumerate(q['bornes_composantes']):
            self.assertLessEqual(abs(ref[i]-F(x[i])), h)
        g = [[F(float.fromhex(v)) for v in l] for l in d['contraintes']]
        for i, h in enumerate(q['bornes_reaction_generalisee']):
            erreur = abs(sum((g[j][i]*(ref[3+j]-F(x[3+j])) for j in range(4)), F(0)))
            self.assertLessEqual(erreur, h)
        # Ajouter un vecteur non nul de ker(G.T) change les multiplicateurs,
        # sans changer l'accélération ni la réaction mécanique résultante.
        delta = [F(-1, 3), F(0), F(1), F(0)]
        self.assertTrue(any(delta))
        self.assertEqual([sum(g[j][i]*delta[j] for j in range(4)) for i in range(3)], [0]*3)

    def test_dependance_non_dyadique(self):
        d = exemple()
        self.assertEqual(d['dependances'][2], ['1/3', '0'])
        self.assertEqual(d['dependances'][3], ['1/3', '1'])
        verifier_certificat(d)

    def test_permutations_et_echelles(self):
        for scale in (2.**-100, 1., 2.**100):
            for p in ([0, 1, 2], [2, 0, 1]):
                g = np.array([[3., 0.], [0., 2.], [1., 2.]])[p]*scale
                actifs = [i != 2 for i in p]
                d = construit(g, [0.]*3, actifs)
                self.assertEqual(verifier_certificat(d)['rang_contraintes'], 2)

    def test_ligne_presque_dependante_refusee(self):
        with self.assertRaisesRegex(CertificationImpossible, 'indépendante'):
            construit([[1., 0.], [1., 2.**-80]], [0., 0.], [True, False])

    def test_second_membre_incompatible_meme_sous_normal(self):
        with self.assertRaisesRegex(CertificationImpossible, 'incompatibles'):
            construit([[1.], [1.]], [0., math.ulp(0.)], [True, False])

    def test_rang_zero(self):
        d = construit([[0., 0.], [0., 0.]], [0., 0.], [False, False])
        q = verifier_certificat(d)
        self.assertEqual(q['rang_contraintes'], 0)
        self.assertEqual(q['bornes_reaction_generalisee'], [0., 0.])

    def test_temoins_et_structure_alteres(self):
        d = exemple(perturbe=True)
        changements = (
            ('dependances', lambda x: x['dependances'][2].__setitem__(0, str(F(1, 3)+F(math.ulp(0.))))),
            ('c', lambda x: x['second_membre_contraintes'].__setitem__(2, math.ulp(0.).hex())),
            ('g', lambda x: x['contraintes'][2].__setitem__(2, math.ulp(0.).hex())),
            ('masque', lambda x: x['actives'].__setitem__(2, True)),
            ('n', lambda x: x.__setitem__('composantes_dynamiques', 2)),
            ('contributions', lambda x: x['contributions_contraintes'][0].__setitem__(2, (4.).hex())),
            ('borne', lambda x: x.__setitem__('bornes_reaction_generalisee', [(0.).hex()]*3)),
            ('fraction', lambda x: x['dependances'][2].__setitem__(0, '1/0')),
        )
        for nom, change in changements:
            faux = copy.deepcopy(d)
            change(faux)
            with self.subTest(nom=nom), self.assertRaises(ValueError):
                verifier_certificat(faux)

    def test_systeme_valide_mais_autre_jauge_refuse(self):
        d = exemple()
        s = d['systeme']
        a = [[float.fromhex(v) for v in l] for l in s['matrice']]
        b = list(map(float.fromhex, s['second_membre']))
        a[-1][-1] = 2.
        # Le sous-certificat est valide, mais son bloc inactif n'est plus celui déclaré.
        d['systeme'] = certifier_systeme(a, b)
        with self.assertRaisesRegex(ValueError, 'jauge'):
            verifier_certificat(d)

    def test_verificateur_autonome(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'quotient.json'
            p.write_text(json.dumps(exemple()))
            r = subprocess.run([sys.executable, '-I', '-S', _verification_lineaire.__file__, str(p)],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)


class QuotientNatif(unittest.TestCase):
    def test_fausse_redondance_du_qr_refusee(self):
        n = Noyau([0., 0., 0.])
        n.corps('mobile', 1., np.eye(3).ravel().tolist(), [0., 0., 0.])
        n.liaison('a', None, 0, bloque_t=[0], bloque_r=[])
        n.liaison('b', None, 0, pa=[0., 2.**-80, 0.], bloque_t=[0], bloque_r=[])
        n.assemble()
        avant = etat(n)
        with self.assertRaisesRegex(CertificationImpossible, 'indépendante'):
            certifier_initialisation(n, redondances=True)
        self.assertEqual(etat(n), avant)

    def test_contact_non_lisse_reste_exclu(self):
        n = Noyau()
        n.corps('bille', 1., np.eye(3).ravel().tolist(), [0., 0., .5])
        n.contact('sol', 0, [0., 0., 0.], .5, nonlisse=True)
        avant = etat(n)
        with self.assertRaisesRegex(CertificationImpossible, 'contact non lisse'):
            certifier_initialisation(n, redondances=True)
        self.assertEqual(etat(n), avant)

    def test_pendule_double_et_solution_native(self):
        n = pendule()
        n.liaison('copie', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.assemble()
        n.enregistre_schema()
        avant = etat(n)
        with patch('numpy.linalg.solve', side_effect=AssertionError('solution externe interdite')):
            d = certifier_initialisation(n, redondances=True, erreur_max=1e-11)
        self.assertEqual(etat(n), avant)
        q = verifier_certificat(d)
        self.assertEqual(q['rang_contraintes'], 5)
        self.assertEqual(q['nombre_contraintes'], 10)
        self.assertFalse(q['multiplicateurs_uniques'])
        n.simule(.001, .001)
        acc, _, lam = n.schema()[0]
        self.assertEqual(d['systeme']['solution'], [v.hex() for v in acc+lam])

    def test_ancien_contrat_conserve(self):
        n = pendule()
        n.liaison('copie', None, 0, bloque_t=[0], bloque_r=[])
        with self.assertRaisesRegex(CertificationImpossible, 'redondantes'):
            certifier_initialisation(n)
        self.assertEqual(verifier_certificat(certifier_initialisation(n, redondances=True))['rang_contraintes'], 5)

    def test_aucune_contrainte(self):
        n = Noyau([0., 0., -9.81])
        n.corps('libre', 1., np.eye(3).ravel().tolist(), [0., 0., 0.])
        q = verifier_certificat(certifier_initialisation(n, redondances=True, erreur_max=0.))
        self.assertEqual(q['rang_contraintes'], 0)
        self.assertTrue(q['multiplicateurs_uniques'])

    def test_refus_precision_et_budget_sans_mutation(self):
        n = pendule()
        n.liaison('copie', None, 0, bloque_t=[0], bloque_r=[])
        avant = etat(n)
        for kw in ({'erreur_max': 0.}, {'dimension_max': 4}, {'redondances': 1}):
            with self.assertRaises(ValueError):
                certifier_initialisation(n, **({'redondances': True} | kw))
            self.assertEqual(etat(n), avant)


if __name__ == '__main__':
    unittest.main()
