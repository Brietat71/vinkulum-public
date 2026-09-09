"""API publique creuse : contributions exactes, enveloppes, ordre et autonomie."""
import copy
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix, eye

from vinkulum.certification import certifier_systeme_creux, certifier_quotient_structurel, verifier_certificat, CertificationImpossible
from vinkulum import _verification_lineaire


class CertificationCreuse(unittest.TestCase):
    def structure(self):
        return certifier_quotient_structurel(eye(2, format='csc'), csc_matrix([[1., 0.]]),
                    csc_matrix([[1.], [2.], [-3.]]), [0], [0., 1.], [0.], [0., 1., 0.],
                    delta_c=csc_matrix([[.01, .01]]), delta_t=csc_matrix([[0.], [10.], [10.]]))

    def test_structure_contre_solution_analytique(self):
        d = self.structure()
        r = verifier_certificat(d)
        self.assertTrue(r['acceleration_unique'])
        self.assertFalse(r['multiplicateurs_uniques'])
        for a in (F(99, 100), F(101, 100)):
            for b in (F(-1, 100), F(1, 100)):
                den = a*a+b*b
                exact = [-a*b/den, a*a/den, b/den]
                for i, cible in enumerate((0, 1, 0)):
                    self.assertLessEqual(abs(exact[i]-cible), r['bornes_composantes'][i])
                for i, v in enumerate((a, b)):
                    self.assertLessEqual(abs(v*exact[2]), r['bornes_reaction_generalisee'][i])

    def test_structure_falsifiee_et_refus(self):
        d = self.structure()
        for mutation in (lambda x: x.update(base=[1]),
                         lambda x: x.update(bornes_reaction=['0', '0']),
                         lambda x: x['delta_t'].__setitem__(0, [[0, '1/100']]),
                         lambda x: x['c'].__setitem__(0, [[0, '2']]),
                         lambda x: x['systeme'].pop('a')):
            faux = copy.deepcopy(d); mutation(faux)
            with self.assertRaises(ValueError): verifier_certificat(faux)
        with self.assertRaises(CertificationImpossible):
            certifier_quotient_structurel(eye(2, format='csc'), csc_matrix([[1., 0.]]),
                csc_matrix([[1.], [2.]]), [0], [0., 1.], [0.], delta_c=csc_matrix([[1., 0.]]))

    def test_doublons_conserves_exactement_et_nonmutation(self):
        a = coo_matrix(([1e16, 1., -1e16], ([0, 0, 0], [0, 0, 0])), shape=(1, 1))
        avant = (a.data.copy(), a.row.copy(), a.col.copy())
        d = certifier_systeme_creux(a, [1.])
        self.assertEqual(d['a'], [[[0, '1']]])
        self.assertEqual(d['x'], ['1'])
        self.assertEqual(verifier_certificat(d)['borne_erreur_inf'], 0)
        for x, y in zip(avant, (a.data, a.row, a.col), strict=True): np.testing.assert_array_equal(x, y)

    def test_enveloppes_et_facteurs_fournis(self):
        d = certifier_systeme_creux(csc_matrix([[2.]]), [1.], [.5],
                  facteurs=(eye(1, format='csc'), csc_matrix([[2.]])),
                  delta_a=csc_matrix([[.25]]), delta_b=[F(1, 10)])
        r = verifier_certificat(d)
        self.assertTrue(r['unicite_uniforme'])
        self.assertGreaterEqual(r['borne_erreur_inf'], F(9, 70))
        self.assertLess(r['borne_erreur_inf'], F(13, 100))
        d['bornes_composantes'] = ['1/8']
        with self.assertRaises(ValueError): verifier_certificat(d)

    def test_permutations_et_solution_originale(self):
        a = csc_matrix([[0., 2., 0.], [1., 0., 0.], [0., 0., 4.]])
        d = certifier_systeme_creux(a, [4., 1., 12.])
        self.assertEqual(list(map(F, d['x'])), [1, 2, 3])
        self.assertEqual(verifier_certificat(d)['bornes_composantes'], [0, 0, 0])
        self.assertNotEqual((d['pr'], d['pc']), ([0, 1, 2], [0, 1, 2]))

    def test_refus_de_domaine_et_de_preuve(self):
        a = eye(2, format='csc')
        for args, kwargs in [((np.eye(2), [0, 0]), {}), ((a, [0]), {}),
                             ((a.astype(complex), [0, 0]), {}), ((a, [True, 0]), {}),
                             ((a, [0, 0]), {'poids': [0, 1]}),
                             ((a, [0, 0]), {'delta_a': -a}),
                             ((a, [0, 0]), {'delta_a': a}),
                             ((a, [0, 0]), {'permutations': ([0, 1], [0, 1])}),
                             ((a, [0, 0]), {'facteurs': (a, a)})]:
            with self.assertRaises(CertificationImpossible): certifier_systeme_creux(*args, **kwargs)

    def test_verification_autonome(self):
        d = certifier_systeme_creux(eye(160, format='csc'), [1.]*160)
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'preuve.json'
            p.write_text(json.dumps(d))
            r = subprocess.run([sys.executable, '-S', _verification_lineaire.__file__, str(p)],
                               capture_output=True, text=True, check=True)
            self.assertIn('dimension 160', r.stdout)
            p.write_text(json.dumps(self.structure()))
            r = subprocess.run([sys.executable, '-S', _verification_lineaire.__file__, str(p)],
                               capture_output=True, text=True, check=True)
            self.assertIn('rang 1/3', r.stdout)


if __name__ == '__main__': unittest.main()
