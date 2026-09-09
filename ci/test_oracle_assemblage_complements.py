from fractions import Fraction as F
import unittest
from types import SimpleNamespace
import numpy as np
from scipy.sparse import csr_matrix
from oracle_assemblage_complements import OracleAssemblage
from test_assemblage_complements import physique
from test_krylov_contraint import flottant

class OracleAssemblageExact(unittest.TestCase):
    def test_assemblage_oblique_vs_fraction(self):
        pieces=[];controls=[]
        for scale in (1.,1.5):
            d=scale*(np.diag(np.arange(2.,8))+np.diag(np.full(5,.25),1))
            r=np.eye(6)+np.diag(np.full(5,.125),1);m=r.T@r
            pieces.append((csr_matrix(d),csr_matrix(m)))
            controls.append(SimpleNamespace(reduction=SimpleNamespace(d=csr_matrix(d),m=csr_matrix(m),i=np.arange(4),s=np.arange(4,6))))
        apps=[np.eye(2),np.array([[1.,.25],[-.5,1.]])];de=np.diag([.5,1.]);me=np.diag([.125,.25]);f=np.array([[1.,0.,.5],[0.,1.,-.5]])
        for precision in (50,70):
            oracle=OracleAssemblage(pieces,apps,de,me,dps=precision)
            for w in (0.,.25,1.):
                exact,fields,*_=physique(controls,apps,de,me,w,f);rep=oracle.reponses(w,f,retour_decimal=True)
                for a,b in [(rep['champ_ports'],exact[:2])]+list(zip(rep['champs'],fields,strict=True)):
                    for row,ref in zip(a,b,strict=True):
                        for v,e in zip(row,ref,strict=True):self.assertLessEqual(abs(F(v)-e),F(1,10**(precision-10)))
    def test_pole_local_admis_pole_global_refuse(self):
        for port,me,accepted in ((1.,0.,True),(2.,3.,False)):
            d=csr_matrix(np.diag([2.,3.,port]));m=csr_matrix(np.eye(3));de=np.array([[2. if accepted else 0.]])
            oracle=OracleAssemblage([(d,m)],[np.ones((1,1))],de,np.array([[me]]))
            if accepted:np.testing.assert_array_equal(oracle.reponses(1.,np.ones((1,1)))['champ_ports'],[[.25]])
            else:
                with self.assertRaises(ValueError):oracle.reponses(1.,np.ones((1,1)))

if __name__=='__main__':unittest.main()
