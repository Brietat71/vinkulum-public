"""Forme cinétique P1 et sélection généralisée avec M complet."""
import unittest
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from modeles_masse_consistante import masse_consistante
from experience_masse_comparee import selection

class MasseP1(unittest.TestCase):
    def test_integrale_vitesses_et_encastrement(self):
        n=7;model=dict(longueur_m=1.,largeur_y_m=.04,hauteur_z_m=.006,rho_kg_m3=2700.)
        m,info=masse_consistante(n,model)
        v=np.random.default_rng(812).normal(size=(n+1,6));v[0]=0
        weights=info['masse_segment_kg']*np.r_[np.ones(3),info['moments_quadratiques_section_m2']]
        expected=sum(np.dot(weights,a*a+a*b+b*b)/3 for a,b in zip(v[:-1],v[1:],strict=True))
        x=v[1:].ravel();self.assertAlmostEqual(x@(m@x),expected,places=14)
        diagonal=m.diagonal();eig=eigh(m.toarray()/np.sqrt(diagonal[:,None]*diagonal[None,:]),eigvals_only=True)
        self.assertGreaterEqual(eig[0],.5-1e-14);self.assertLessEqual(eig[-1],1.5+1e-14)
        self.assertGreater(m.nnz,len(diagonal))

    def test_selection_utilise_masse_complete(self):
        n=30;rng=np.random.default_rng(918)
        a=rng.normal(size=(n,n));k=a.T@a+np.eye(n)
        b=rng.normal(size=(n,n));m=b.T@b+np.eye(n)
        class Q:
            def produit(self,v):return k@v
            def solve(self,v):return np.linalg.solve(k,v)
        phi,con,info=selection(Q(),csr_matrix(m),32)
        exact=eigh(k,m,eigvals_only=True)[0]
        self.assertAlmostEqual(1/info['valeur_inverse'][0],exact,places=10)
        np.testing.assert_allclose(con,m@phi,rtol=1e-14,atol=1e-14)
        self.assertLess(np.linalg.norm(k@phi-exact*m@phi)/np.linalg.norm(k@phi),1e-8)

if __name__=='__main__':unittest.main()
