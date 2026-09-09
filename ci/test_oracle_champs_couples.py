"""Oracle bloc couplé face à l'inverse Fraction générale."""
from fractions import Fraction as F
import unittest
import numpy as np
from scipy.sparse import csr_matrix
from oracle_champs_couples import OracleChampsCouples
from test_krylov_contraint import oracle,flottant


class OracleCouple(unittest.TestCase):
    def test_masse_couplee_et_charges_multiples(self):
        n,p=6,2
        d=np.diag(np.arange(2.,n+2.))+np.diag(np.full(n-1,.125),1)
        c=np.eye(n)+np.diag(np.full(n-1,.25),1);m=c.T@c
        f=np.array([[1.,0.,2.**-40,0.],[0.,1.,-2.**-40,0.]])
        for precision in (50,70):
            reference=OracleChampsCouples(csr_matrix(d),csr_matrix(m),p=p,dps=precision)
            for w in (0.,.25,1.):
                actual=reference.reponse(w,f,retour_decimal=True);expected=oracle(d,m,[4,5],w,f)
                for row,target in zip(actual,expected):
                    for x,y in zip(row,target):self.assertLessEqual(abs(F(x)-y),F(1,10**(precision-8)))
                np.testing.assert_array_equal(reference.reponse(w,f)[:,3],np.zeros(n))
        diagonal=OracleChampsCouples(csr_matrix(d),csr_matrix(np.diag(m.diagonal())),p=p,dps=70)
        self.assertGreater(np.linalg.norm(diagonal.reponse(1.,f)-reference.reponse(1.,f)),1e-5)

    def test_blocs_scalaires_et_zero_masse_permis(self):
        d=np.array([[2.,.25,0.],[0.,3.,.125],[0.,0.,4.]])
        m=np.array([[1.,.125,0.],[.125,2.,0.],[0.,0.,0.]])
        f=np.array([[1.,-1.]])
        ref=OracleChampsCouples(csr_matrix(d),csr_matrix(m),p=1,dps=70)
        np.testing.assert_allclose(ref.reponse(.5,f),flottant(oracle(d,m,[2],.5,f)),rtol=1e-15,atol=0)

    def test_couplage_hors_bande_asymetrie_et_singularite_refuses(self):
        d=csr_matrix(np.eye(4));m=np.eye(4);m[0,3]=m[3,0]=.125
        with self.assertRaisesRegex(ValueError,'non voisins'):OracleChampsCouples(d,csr_matrix(m),p=1)
        m[0,3]=0
        with self.assertRaisesRegex(ValueError,'symétrique'):OracleChampsCouples(d,csr_matrix(m),p=1)
        ref=OracleChampsCouples(d,csr_matrix(np.eye(4)),p=1)
        with self.assertRaisesRegex(ValueError,'pivot intérieur'):ref.reponse(1.,np.ones((1,1)))


if __name__=='__main__':unittest.main()
