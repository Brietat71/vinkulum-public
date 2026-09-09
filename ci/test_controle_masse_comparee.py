"""Champs et majorants sur masses connectées, contre oracles Fraction."""
from fractions import Fraction as F
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from comparaison_masse import compiler_comparaison,BibliothequeComparaison
from controle_masse_comparee import ControleMasseComparee
from controle_facteurs import ControleFacteurs
from krylov_contraint import KrylovContraint
from vinkulum._ports.condensation_energie import CondensationEnergie
from test_krylov_contraint import donnees,oracle,rationnel,flottant,energie
from test_trace_complement import gram


class ControleMasseCompareePhysique(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='vinkulum-controle-masse-test-')
        cls.bib=BibliothequeComparaison(compiler_comparaison(Path(cls.temp.name)/'build'))
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()
    def forces(self):return np.array([[1.,0.,1.,2.**-30,0.,1.],[0.,1.,1.,-2.**-30,0.,-1.+2.**-30]])

    def verifier_bornes(self,r,d,m,ports,w,f):
        rep=r.reponses(w,f);exact=oracle(d,m,ports,w,f);x=rationnel(rep['champ'])
        self.assertGreater(rep['marge'],0.);self.assertFalse(rep['certification_machine'])
        matrices=dict(masse=rationnel(m),deformation=gram(rationnel(d)))
        for j in range(f.shape[1]):
            erreur=[v[j]-e[j] for v,e in zip(x,exact)]
            for name,a in matrices.items():
                b=rep['bornes'][name]['absolues'][j]
                self.assertIsNotNone(b);self.assertLessEqual(energie(erreur,a),F(float(b))**2)
                norm=rep['bornes'][name]['normes'][j]
                physical=energie([v[j] for v in x],a)
                self.assertLessEqual(abs(F(float(norm))**2-physical),F(1,10**11)*physical)
        return rep

    def test_masse_connectee_taille_huit_complement_incomplet(self):
        d,_,i,s,b,phi,qr=donnees(8);n=d.shape[1]
        factor=np.eye(n)+np.diag(np.full(n-1,.125),1);m=factor.T@factor
        self.assertEqual(connected_components(csr_matrix(m[np.ix_(i,i)]),directed=False)[0],1)
        reduction=KrylovContraint(qr,m,b,phi,4.,.5,blocs=1,max_directions=7)
        self.assertLess(reduction.taille_complement_reduit,7)
        with self.assertRaisesRegex(ValueError,'taille de bloc'):ControleFacteurs(reduction)
        control=ControleMasseComparee(reduction,bibliotheque=self.bib)
        self.assertGreater(np.linalg.norm(m[np.ix_(i,s)]),.01)
        for w in (0.,.25,.5):self.verifier_bornes(control,d,m,s,w,self.forces())
        np.testing.assert_array_equal(reduction.m.toarray(),m)

    def test_champs_identiques_ancien_sur_domaine_commun(self):
        d,m,i,s,b,phi,qr=donnees()
        reduction=KrylovContraint(qr,m,b,phi,4.,.5,blocs=1,max_directions=2)
        a=ControleFacteurs(reduction);b=ControleMasseComparee(reduction,bibliotheque=self.bib)
        for w in (0.,.125,.5):
            x=a.reponses(w,self.forces());y=b.reponses(w,self.forces())
            if y['marge']>0:self.verifier_bornes(b,d,m,s,w,self.forces())
            else:
                for name in ('masse','deformation'):
                    self.assertIsNone(y['bornes'][name]['absolues'])
                    self.assertTrue(all(v is None for v in y['bornes'][name]['relatives']))
            for name in ('champ','coordonnees'):np.testing.assert_array_equal(x[name],y[name])
            for name in ('masse','deformation'):np.testing.assert_array_equal(x['bornes'][name]['normes'],y['bornes'][name]['normes'])

    def test_deux_contraintes_obliques_et_changement_unites(self):
        d,_,i,s,_,_,_=donnees(8);n=d.shape[1]
        factor=np.eye(n)+np.diag(np.full(n-1,.125),1);m=factor.T@factor
        b=np.array([[1.,0.],[0.,1.],[1.,-1.],[2.,1.],[-1.,2.],[1.,1.],[0.,1.],[1.,0.]])
        phi=b.copy();f=self.forces()
        for scale in (np.ones(n),np.ldexp(np.ones(n),np.arange(n)-5)):
            dd=d*scale;mm=m*scale[:,None]*scale[None,:]
            qr=CondensationEnergie(csr_matrix(dd),i,s,np.diag([4.,.25])*scale[s,None]*scale[None,s])
            reduction=KrylovContraint(qr,mm,scale[i,None]*b,phi/scale[i,None],4.,.5,blocs=1,max_directions=4)
            control=ControleMasseComparee(reduction,bibliotheque=self.bib)
            self.verifier_bornes(control,dd,mm,s,.5,scale[s,None]*f)

    def test_masse_psd_ports_sans_masse_et_poles(self):
        d=np.array([[1.,0.,1.],[0.,2.,1.],[0.,0.,1.]])
        for zero_port in (False,True):
            m=np.diag([1.,1.,0. if zero_port else 1.]);b=np.array([[1.],[0.]])
            qr=CondensationEnergie(csr_matrix(d),np.arange(2),np.array([2]),np.eye(1))
            reduction=KrylovContraint(qr,m,b,b,3.,1.25,blocs=2,max_directions=1)
            control=ControleMasseComparee(reduction,bibliotheque=self.bib)
            self.verifier_bornes(control,d,m,[2],1.,np.array([[1.,-.5,2.**-30,0.]]))
        d=np.diag([1.,2.,1.]);m=np.eye(3)
        qr=CondensationEnergie(csr_matrix(d),np.arange(2),np.array([2]),np.eye(1))
        reduction=KrylovContraint(qr,m,b,b,3.,1.25,blocs=2,max_directions=1)
        control=ControleMasseComparee(reduction,bibliotheque=self.bib)
        with self.assertRaises(np.linalg.LinAlgError):control.reponses(1.,np.ones((1,1)))

    def test_controle_perime_et_comparaison_trop_large(self):
        d,m,i,s,b,phi,qr=donnees(8)
        reduction=KrylovContraint(qr,m,b,phi,4.,.5,blocs=1,max_directions=7)
        control=ControleMasseComparee(reduction,bibliotheque=self.bib,alpha=2.**-20,beta=2.)
        rep=control.reponses(.5,self.forces())
        self.assertTrue(any(v is None or v>1e-6 for v in rep['bornes']['masse']['relatives']))
        self.assertTrue(reduction.enrichir())
        with self.assertRaisesRegex(ValueError,'périmé'):control.reponses(.5,self.forces())


if __name__=='__main__':unittest.main()
