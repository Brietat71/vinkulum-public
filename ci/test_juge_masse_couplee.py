"""Défaut du facteur du juge vérifié sur les matrices rationnelles exactes."""
from fractions import Fraction as F
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy.sparse import csr_matrix

from comparaison_masse import compiler_comparaison,BibliothequeComparaison
from racine_masse_juge import RacineMasseJuge
from juge_masse_couplee import juger_couple
from test_trace_complement import gram,difference,echelle
from test_krylov_contraint import rationnel
from test_inertie_complement import inertie_fraction


class JugeMasseCouplee(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='vinkulum-juge-couple-test-')
        cls.bib=BibliothequeComparaison(compiler_comparaison(Path(cls.temp.name)/'build'))
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def test_defaut_rationnel_et_correction_relative(self):
        for exponent,shift in ((0,0.),(20,0.),(0,.125),(20,.125)):
            n=8;root=np.eye(n)+np.diag(np.full(n-1,.125),1)
            scale=np.ldexp(np.ones(n),np.arange(n)*exponent)
            m=(root.T@root+shift*np.eye(n))*scale[:,None]*scale[None,:]
            r=RacineMasseJuge(csr_matrix(m),bibliotheque=self.bib)
            a=r.audit;rho=F(int(a['rho_numerateur']),int(a['rho_denominateur']));alpha=F(a['alpha'])
            if shift:self.assertGreater(rho,0)
            else:self.assertEqual(rho,0)
            self.assertGreaterEqual(F(r.correction)**2,(alpha+rho)/(alpha-rho))
            rr=gram(rationnel(r.r.toarray()));mm=rationnel(m);epsilon=rho/alpha
            for gap in (difference(echelle(mm,1+epsilon),rr),difference(rr,echelle(mm,1-epsilon))):
                self.assertEqual(inertie_fraction(gap)[1],0)

    def test_jugement_champs_et_operateurs_couples(self):
        n=12;root=np.eye(n)+np.diag(np.full(n-1,.125),1);m=csr_matrix(root.T@root)
        d=csr_matrix(np.diag(np.arange(1.,n+1)));reference=np.eye(n)[:,-6:][None,:,:]
        exact=juger_couple(d,m,np.eye(6),reference,reference,bibliotheque=self.bib)
        self.assertTrue(exact['accepte']);self.assertEqual(max(exact['maxima'].values()),0.)
        self.assertGreaterEqual(exact['audit_masse']['correction_relative'],1.)
        for amplitude,accepted in ((1e-8,True),(1e-4,False)):
            candidate=reference.copy();candidate[0,0,0]=amplitude
            result=juger_couple(d,m,np.eye(6),candidate,reference,bibliotheque=self.bib)
            self.assertEqual(result['accepte'],accepted)
            self.assertGreater(result['maxima']['masse'],0.)
            self.assertGreater(result['maxima_operateurs']['masse'],0.)

    def test_budgets_et_rang_refuses(self):
        m=csr_matrix(np.eye(8)+np.diag(np.full(7,.125),1)+np.diag(np.full(7,.125),-1))
        for kw in (dict(budget_coefficients=1),dict(budget_produits=1)):
            with self.assertRaisesRegex(ValueError,'budget'):RacineMasseJuge(m,bibliotheque=self.bib,**kw)
        reference=np.ones((1,8,6))
        with self.assertRaisesRegex(ArithmeticError,'rang inférieur'):
            juger_couple(csr_matrix(np.eye(8)),m,np.eye(6),reference,reference,bibliotheque=self.bib)


if __name__=='__main__':unittest.main()
