"""Confrontation à l'assemblage physique rationnel, sans base réduite oracle."""
from fractions import Fraction as F
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy.sparse import csr_matrix
from comparaison_masse import compiler_comparaison,BibliothequeComparaison
from controle_masse_comparee import ControleMasseComparee
from krylov_contraint import KrylovContraint
from condensation_energie import CondensationEnergie
from assemblage_complements import AssemblageComplements
from test_krylov_contraint import donnees,rationnel,flottant,energie
from test_trace_complement import gram,produit,inverse,difference,echelle,restreindre


def physique(controls,maps,de,me,omega,forces):
    g=maps[0].shape[1];n=g+sum(len(c.reduction.i) for c in controls)
    k=[[F(0)]*n for _ in range(n)];m=[[F(0)]*n for _ in range(n)];emb=[];offset=g
    for c,a in zip(controls,maps,strict=True):
        r=c.reduction;p=np.zeros((r.d.shape[1],n));p[r.s,:g]=a;p[r.i,offset:offset+len(r.i)]=np.eye(len(r.i));offset+=len(r.i)
        p=rationnel(p);emb.append(p)
        kj=restreindre(gram(rationnel(r.d.toarray())),p);mj=restreindre(rationnel(r.m.toarray()),p)
        for i in range(n):
            for j in range(n):k[i][j]+=kj[i][j];m[i][j]+=mj[i][j]
    ke=gram(rationnel(de));mm=rationnel(me)
    for i in range(g):
        for j in range(g):k[i][j]+=ke[i][j];m[i][j]+=mm[i][j]
    rhs=rationnel(np.vstack((forces,np.zeros((n-g,forces.shape[1])))))
    u=produit(inverse(difference(k,echelle(m,F(float(omega))**2))),rhs)
    return u,[produit(p,u) for p in emb],k,m


class AssemblagePhysique(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory(prefix='vinkulum-assemblage-test-')
        cls.bib=BibliothequeComparaison(compiler_comparaison(Path(cls.tmp.name)/'build'))
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def composant(self,factor=1.):
        d,_,i,s,b,phi,_=donnees(8);d*=factor;n=d.shape[1]
        root=np.eye(n)+np.diag(np.full(n-1,.125),1);m=root.T@root
        qr=CondensationEnergie(csr_matrix(d),i,s,np.diag([4.,.25]))
        r=KrylovContraint(qr,m,b,phi,4.,.5,blocs=1,max_directions=7)
        return ControleMasseComparee(r,bibliotheque=self.bib)
    def verifier(self,a,w,f):
        rep=a.reponses(w,f);exact,fields,k,m=physique(a.controles,a.applications,a.de,a.me,w,f)
        self.assertGreater(rep['marge'],0.);self.assertFalse(rep['certification_machine'])
        for column in range(f.shape[1]):
            err=[F(float(x))-v[column] for x,v in zip(rep['champ_ports'][:,column],exact[:a.g],strict=True)]
            b=rep['borne_coordonnees'][column]
            self.assertLessEqual(energie(err,rationnel(a.metrique)),F(float(b))**2)
            energy={name:F(0) for name in ('masse','deformation')}
            for c,candidate,ref in zip(a.controles,rep['champs'],fields,strict=True):
                e=[F(float(x))-v[column] for x,v in zip(candidate[:,column],ref,strict=True)]
                energy['masse']+=energie(e,rationnel(c.reduction.m.toarray()))
                energy['deformation']+=energie(e,gram(rationnel(c.reduction.d.toarray())))
            energy['masse']+=energie(err,rationnel(a.me));energy['deformation']+=energie(err,gram(rationnel(a.de)))
            for name in energy:
                self.assertLessEqual(energy[name],F(float(rep['bornes'][name]['absolues'][column]))**2)
        return rep
    def test_deux_pieces_obliques_masses_couplees_et_energies_externes(self):
        cs=[self.composant(),self.composant(1.5)]
        maps=[np.array([[1.,.25],[-.5,1.]]),np.array([[.5,-1.],[1.,.5]])]
        de=np.array([[.5,.25],[0.,1.]]);me=np.diag([.125,.25])
        a=AssemblageComplements(cs,maps,np.diag([2.,.5]),facteur_externe=de,masse_externe=me)
        self.assertEqual(a.taille_conservee,4)
        f=np.array([[1.,0.,1.,2.**-30,0.],[0.,1.,1.,-2.**-30,0.]])
        for w in (0.,.25,.5):
            rep=self.verifier(a,w,f)
            for field in rep['champs']:np.testing.assert_array_equal(field[:,-1],0.)
    def test_un_composant_retrouve_service_isole(self):
        c=self.composant();a=AssemblageComplements([c],[np.eye(2)],np.diag([4.,.25]))
        f=np.array([[1.,.5],[.25,1.]])
        old=c.reponses(.5,f);new=a.reponses(.5,f)
        np.testing.assert_array_equal(old['champ'],new['champs'][0])
        np.testing.assert_array_equal(old['coordonnees'],new['coordonnees'])
        for name in ('masse','deformation'):
            np.testing.assert_allclose(old['bornes'][name]['absolues'],new['bornes'][name]['absolues'],rtol=1e-14)
    def test_reseau_ports_partiellement_partages(self):
        cs=[self.composant(),self.composant(1.5)]
        maps=[np.array([[1.,0.,0.],[0.,1.,0.]]),np.array([[0.,1.,0.],[0.,0.,1.]])]
        a=AssemblageComplements(cs,maps,np.diag([4.,4.25,.25]),facteur_externe=np.eye(3)*.125)
        self.assertEqual(a.taille_conservee,5)
        self.verifier(a,.5,np.eye(3))
    def simple(self,port):
        d=np.diag([2.,3.,port]);qr=CondensationEnergie(csr_matrix(d),np.arange(2),np.array([2]),np.eye(1))
        b=np.array([[1.],[0.]])
        return ControleMasseComparee(KrylovContraint(qr,np.eye(3),b,b,8.,1.25,blocs=1,max_directions=1),bibliotheque=self.bib)
    def test_pole_local_traverse_sans_inversion_locale(self):
        c=self.simple(1.)
        with self.assertRaises(np.linalg.LinAlgError):c.reponses(1.,np.ones((1,1)))
        a=AssemblageComplements([c],[np.ones((1,1))],np.eye(1),facteur_externe=np.array([[2.]]))
        rep=self.verifier(a,1.,np.array([[1.,0.,-.5]]))
        np.testing.assert_array_equal(rep['champ_ports'],[[.25,0.,-.125]])
    def test_pole_global_refuse_alors_que_piece_reguliere(self):
        c=self.simple(2.);c.reponses(1.,np.ones((1,1)))
        a=AssemblageComplements([c],[np.ones((1,1))],np.eye(1),masse_externe=np.array([[3.]]))
        with self.assertRaises(np.linalg.LinAlgError):a.reponses(1.,np.ones((1,1)))
    def test_changement_unites_globales_et_modes_prives(self):
        c=self.composant();f=np.array([[1.,.5],[.25,1.]])
        maps=[np.eye(2),np.array([[1.,.25],[-.25,1.]])];metric=np.diag([2.,.5]);scale=np.diag([8.,.125])
        a=AssemblageComplements([c,c],maps,metric);b=AssemblageComplements([c,c],[p@scale for p in maps],scale@metric@scale)
        x=a.reponses(.5,f);y=b.reponses(.5,scale@f)
        self.assertEqual(a.taille_conservee,4)
        for u,v in zip(x['champs'],y['champs'],strict=True):np.testing.assert_allclose(u,v,rtol=1e-13,atol=1e-15)
        np.testing.assert_allclose(x['champ_ports'],scale@y['champ_ports'],rtol=1e-13,atol=1e-15)
    def test_refus_entrees_et_peremption(self):
        c=self.composant()
        for options in (dict(budget_conserve=1),dict(masse_externe=-np.eye(2))):
            with self.assertRaises(ValueError):AssemblageComplements([c],[np.eye(2)],np.eye(2),**options)
        a=AssemblageComplements([c],[np.eye(2)],np.eye(2))
        for w,f in ((.6,np.eye(2)),(.5,np.eye(2,dtype=complex)),(.5,np.ones((3,1)))):
            with self.assertRaises(ValueError):a.reponses(w,f)
        self.assertTrue(c.reduction.enrichir())
        with self.assertRaisesRegex(ValueError,'périmé'):a.reponses(.5,np.eye(2))
    def test_marge_non_positive_refuse_les_bornes_sans_changer_le_champ(self):
        c=self.composant()
        weak=ControleMasseComparee(c.reduction,bibliotheque=self.bib,alpha=2.**-40)
        good=AssemblageComplements([c],[np.eye(2)],np.eye(2)).reponses(.5,np.eye(2))
        bad=AssemblageComplements([weak],[np.eye(2)],np.eye(2)).reponses(.5,np.eye(2))
        self.assertGreater(good['marge'],0.)
        self.assertLessEqual(bad['marge'],0.)
        np.testing.assert_array_equal(good['champs'][0],bad['champs'][0])
        self.assertIsNone(bad['borne_coordonnees'])
        for name in ('masse','deformation'):
            self.assertIsNone(bad['bornes'][name]['absolues'])
            self.assertTrue(all(v is None for v in bad['bornes'][name]['relatives']))

if __name__=='__main__':unittest.main()
