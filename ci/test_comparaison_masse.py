"""Preuves rationnelles de la comparaison de normes et de son certificat."""
from fractions import Fraction as F
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.sparse import csr_matrix,diags

from comparaison_masse import compiler_comparaison,BibliothequeComparaison,certifier_comparaison_masse
from inertie_complement_dirigee import InertieImpossible
import test_inertie_complement_dirigee as temoin
from test_inertie_complement import inertie_fraction
from test_trace_complement import produit,transpose,inverse,echelle,difference


def rationnel(a):return [[F(float(v)) for v in row] for row in np.asarray(a)]
def dual(m,z):return produit(produit(z,inverse(produit(produit(transpose(z),m),z))),transpose(z))


class ComparaisonMasse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='vinkulum-comparaison-masse-test-')
        cls.bib=BibliothequeComparaison(compiler_comparaison(Path(cls.temp.name)/'build'))
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()
    def contient(self,a,x):
        self.assertLessEqual(F(float.fromhex(a[0])),x);self.assertGreaterEqual(F(float.fromhex(a[1])),x)
    def preuve_exacte(self,a,p):temoin.InertieComplementDirigee.preuve_exacte(self,a,p)

    def test_certificat_masses_connectees_et_echelles_physiques(self):
        for n in (3,7,12):
            for reverse in (False,True):
                root=np.eye(n)+np.diag(np.full(n-1,.25),1)
                scale=np.ldexp(np.ones(n),np.arange(n)-n//2)
                m=(root.T@root)*scale[:,None]*scale[None,:]
                p=np.arange(n)[::-1] if reverse else np.arange(n)
                result=certifier_comparaison_masse(csr_matrix(m),.25,2.,bibliotheque=self.bib,permutation=p)
                l=rationnel(np.diag(m.diagonal()));mm=rationnel(m)
                self.preuve_exacte(difference(mm,echelle(l,F(1,4))),result['preuve_inferieure'])
                self.preuve_exacte(difference(echelle(l,F(2)),mm),result['preuve_superieure'])
                self.assertTrue(result['certification_machine'])

    def test_ordre_loewner_inverse_sur_noyau_oblique(self):
        n=5
        root=np.eye(n)+np.diag(np.full(n-1,.25),1)
        m=rationnel(root.T@root);l=[[m[i][i] if i==j else F(0) for j in range(n)] for i in range(n)]
        z=rationnel([[1,2,-1],[-1,1,2],[1,0,0],[0,1,0],[0,0,1]])
        b=rationnel([[1,0],[0,1],[-1,1],[-2,-1],[1,-2]])
        self.assertEqual(produit(transpose(b),z),[[F(0)]*3 for _ in range(2)])
        sm,sl=dual(m,z),dual(l,z)
        for gap in (difference(m,echelle(l,F(1,4))),difference(echelle(l,F(2)),m)):
            self.assertEqual(inertie_fraction(gap),(5,0,0))
        for gap in (difference(sl,echelle(sm,F(1,4))),difference(echelle(sm,F(2)),sl)):
            self.assertEqual(inertie_fraction(gap),(3,0,2))
        # Une réaction B nu disparaît dans la norme duale contrainte.
        self.assertEqual(produit(sm,b),[[F(0)]*2 for _ in range(n)])
        self.assertEqual(produit(sl,b),[[F(0)]*2 for _ in range(n)])

    def test_comparaison_ne_remplace_pas_la_masse_physique(self):
        m=rationnel([[2,1],[1,2]]);l=rationnel([[2,0],[0,2]])
        k=rationnel([[4,0],[0,9]]);f=rationnel([[1],[0]])
        physical=produit(inverse(difference(k,m)),f)
        diagonal=produit(inverse(difference(k,l)),f)
        self.assertNotEqual(physical,diagonal)
        self.assertEqual(physical,[[F(7,13)],[F(1,13)]])

    def test_refus_bornes_fausses_singularites_et_budgets(self):
        m=csr_matrix([[2.,1.],[1.,2.]])
        for alpha,beta in ((.75,2.),(.25,1.25),(.5,2.),(.25,1.5)):
            with self.assertRaises(InertieImpossible):certifier_comparaison_masse(m,alpha,beta,bibliotheque=self.bib)
        for kw in (dict(budget_operations=1),dict(budget_coefficients=2)):
            with self.assertRaises(InertieImpossible):certifier_comparaison_masse(m,.25,2.,bibliotheque=self.bib,**kw)
        for alpha,beta in ((0,2),(.5,float('inf')),(1j,2),(2,1)):
            with self.assertRaises(ValueError):certifier_comparaison_masse(m,alpha,beta,bibliotheque=self.bib)
        for bad in (csr_matrix([[2.,1.],[0.,2.]]),csr_matrix([[0.,0.],[0.,2.]])):
            with self.assertRaises(ValueError):certifier_comparaison_masse(bad,.25,2.,bibliotheque=self.bib)
        with self.assertRaises(ValueError):certifier_comparaison_masse(m,.25,2.,bibliotheque=self.bib,permutation=[0,0])

    def test_diagonale_auxiliaire_arbitraire_et_identite(self):
        m=csr_matrix([[2.,.5],[.5,3.]])
        result=certifier_comparaison_masse(m,.25,4.,diagonale=np.array([1.,2.]),bibliotheque=self.bib)
        self.assertEqual(result['diagonale'],[1.,2.])
        self.preuve_exacte(rationnel([[1.75,.5],[.5,2.5]]),result['preuve_inferieure'])
        self.preuve_exacte(rationnel([[2.,-.5],[-.5,5.]]),result['preuve_superieure'])
        for ell in ([1.],[-1.,2.],[np.nan,2.]):
            with self.assertRaises(ValueError):certifier_comparaison_masse(m,.25,4.,diagonale=ell,bibliotheque=self.bib)


if __name__=='__main__':unittest.main()
