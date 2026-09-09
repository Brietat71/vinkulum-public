"""Matrices creuses et produits adjoints : contre-épreuves mécaniques."""
import unittest
from unittest.mock import patch

import numpy as np
from scipy.sparse import csc_matrix

from vinkulum import Noyau
from vinkulum.adjoint_temps import PontNoyau, _pendule_spherique, _pendule_flexible
from vinkulum.test_analyses import oscillateur


def matrices(n, t=None):
    return [csc_matrix((data, indices, pointers), shape=(nr,nc))
            for nr,nc,pointers,indices,data in n.k_c_m_g_creux(t)]


class OperateursMecaniques(unittest.TestCase):
    def test_etats_non_portes_par_adjoint_refuses(self):
        for option in ('sphere','appariement','inflow','non_holonome'):
            def model(p):
                n=oscillateur()
                if option=='sphere':n.sphere(1,[0,0,0],.1)
                elif option=='appariement':n.appariement()
                elif option=='non_holonome':n.liaison('roulement',None,1,bloque_r=[],nh=True)
                else:n.inflow([0,0,1],1.,.1)
                return n
            with self.assertRaisesRegex(ValueError,'hors domaine'):
                PontNoyau(model,lambda N,p:np.zeros((12,1))).aller([1.],.01,.01)

    def test_csc_egale_dense_et_ne_modifie_pas_etat(self):
        models = [Noyau([0,0,0]), _pendule_spherique([9.81]),
                  _pendule_flexible(.75,4), oscillateur(4.,.3,.1)]
        libre = Noyau([0,0,0])
        libre.corps('libre',2.,[1.,.2,0,.2,2.,0,0,0,3.],[0,0,0],w=[2.,3.,4.])
        models.append(libre)
        contact = Noyau([0,0,-9.81])
        contact.corps('contact',1.,np.eye(3).ravel().tolist(),[0,0,.08])
        contact.contact('sol',0,[0,0,0],.1)
        models.append(contact)
        for i,n in enumerate(models):
            with self.subTest(modele=i):
                etat=n.etat_precis()
                ks,cs,ms,gs=matrices(n)
                kd,cd,md,_,gd=n.k_c_m_z()
                for sparse,dense in zip((ks,cs,ms,gs),(kd,cd,md,gd)):
                    self.assertTrue(sparse.has_canonical_format)
                    self.assertTrue(sparse.has_sorted_indices)
                    np.testing.assert_array_equal(sparse.toarray(),np.asarray(dense).reshape(sparse.shape))
                self.assertEqual(n.etat_precis(),etat)
                with self.assertRaises(RuntimeError):
                    n.k_c_m_g_creux(float('nan'))
        # Précontrainte et orientation modifiées par une trajectoire réelle.
        n=models[1]; n.simule(.032,.004)
        kd,cd,md,_,gd=n.k_c_m_z()
        for sparse,dense in zip(matrices(n),(kd,cd,md,gd)):
            np.testing.assert_array_equal(sparse.toarray(),np.asarray(dense).reshape(sparse.shape))

    def test_produit_adjoint_contre_jacobien_complet(self):
        rng=np.random.default_rng(417)
        for model,p,h,steps in [(_pendule_spherique,[9.81],.02,8),
                                (lambda p:_pendule_flexible(p[0],8),[.75],.001,5)]:
            n=6*len(model(p).etat()[1])
            rp=rng.normal(size=(n,1))
            reference=PontNoyau(model,lambda N,p:rp,creux=False)
            N,fr,sch=reference.aller(p,h*steps,h)
            args=(N,fr[-2],fr[-1],np.asarray(sch[-2][0]),np.asarray(sch[-2][1]),
                  np.asarray(sch[-1][0]),np.asarray(sch[-1][1]),p)
            aa,bb=reference._pas(*args)
            for creux in (False,True):
                pont=PontNoyau(model,lambda N,p:rp,creux=creux)
                for _ in range(3):
                    mu=rng.normal(size=4*n)
                    previous,grad=pont._recul(*args,mu)
                    np.testing.assert_allclose(previous,aa.T@mu,rtol=3e-9,atol=2e-7)
                    np.testing.assert_allclose(grad,bb.T@mu,rtol=3e-9,atol=2e-7)
                self.assertLessEqual(pont.residu_adjoint_max,1e-11)
                self.assertGreaterEqual(pont.n_resolutions_adjoint,3)

    def test_gradient_de_trajectoire_contre_differences_finies(self):
        def model(p):
            n=oscillateur(p[0],.4,.1)
            n.effort(1,[p[1],0,0],[0,0,0])
            return n
        def rp(N,p):
            x=N.etat()[1][1][0]-1.
            out=np.zeros((12,2))
            out[0,0],out[6,0],out[6,1]=-x,x,-1.
            return out
        def observable(N):
            frame=N.etat()
            return frame[1][1][0]+.2*frame[3][1][0]
        def dj(N,frame):
            out=np.zeros(48);out[6]=1.;out[18]=.2
            return out
        p=np.array([40.,2.]); t,h,rho=.083,.004,.7
        df=[]
        for j in range(2):
            vals=[];eps=1e-4*max(1.,abs(p[j]))
            for sign in (1.,-1.):
                pp=p.copy();pp[j]+=sign*eps
                n=model(pp);n.simule(t,h,rho=rho)
                vals.append(observable(n))
            df.append((vals[0]-vals[1])/(2*eps))
        for creux in (False,True,None):
            pont=PontNoyau(model,rp,rho=rho,creux=creux)
            with patch.object(pont,'_pas',side_effect=AssertionError('matrice de transition interdite')):
                grad=pont.gradient(p,t,h,dj)
            np.testing.assert_allclose(grad,df,rtol=2e-6,atol=2e-9)
        def initial(p,mu):
            total=mu[30]+mu[42]
            return np.array([-.1*total,total])
        pont=PontNoyau(model,rho=rho,creux=True,
                      d_residu_transpose=lambda N,p,ud,psi:rp(N,p).T@psi[:12],
                      d_initial_transpose=initial)
        grad=pont.gradient(p,t,h,dj)
        np.testing.assert_allclose(grad,df,rtol=2e-6,atol=2e-9)

    def test_sensibilites_poutres_produit_transpose(self):
        rng=np.random.default_rng(542)
        for formulation in ('milieu','integree'):
            n=Noyau([0,0,0])
            for i in range(4):
                n.corps(str(i),1.,np.eye(3).ravel().tolist(),[.3*i,0,0])
            for i in range(3):
                n.poutre(str(i),i,i+1,4000.,800.,30.,60.,ei3=20.,ga3=700.,formulation=formulation)
            frame=list(n.etat())
            for i in range(4):
                frame[1][i][1]=.02*i*i
            n.pose_etat(*frame)
            poids=rng.normal(size=24)
            for quoi in ('ea','ga','gj','ei','ei_e2','ei_e3'):
                ref=np.array([poids@np.asarray(n.d_residu_poutre(i,quoi)) for i in range(3)])
                actual=n.d_residu_poutres_transpose(poids.tolist(),quoi)
                np.testing.assert_allclose(actual,ref,rtol=3e-13,atol=1e-15)
            with self.assertRaises(ValueError):n.d_residu_poutres_transpose([1.])
            with self.assertRaises(ValueError):n.d_residu_poutres_transpose([float('nan')]*24)
            with self.assertRaises(ValueError):n.d_residu_poutres_transpose(poids.tolist(),'inconnue')

    def test_gradient_flexible_creux_contre_difference_finie(self):
        for ne in (8,30,120):
            def model(p):
                return _pendule_flexible(p[0],ne)
            initial=model([.75])
            # Aucun paramètre de masse/liaison ; poutres sans déformation.
            for i in range(ne):
                np.testing.assert_array_equal(initial.d_residu_poutre(i,'ei'),np.zeros(6*(ne+1)))
            def rtp(n,p,ud,psi):
                return np.array([sum(n.d_residu_poutres_transpose(psi[:len(ud)].tolist()))])
            def dj(n,frame):
                out=np.zeros(24*(ne+1));out[6*ne+2]=1.
                return out
            pont=PontNoyau(model,creux=True,d_residu_transpose=rtp,
                          d_initial_transpose=lambda p,mu:np.zeros(1))
            grad=pont.gradient([.75],.02,.001,dj)[0]
            for eps in (.001,.0003,.0001):
                vals=[]
                for sign in (1.,-1.):
                    n=model([.75+sign*eps]);n.simule(.02,.001)
                    vals.append(n.etat()[1][-1][2])
                df=(vals[0]-vals[1])/(2*eps)
                np.testing.assert_allclose(grad,df,rtol=2e-5,atol=1e-10)

    def test_equilibrage_adjoint_et_echelles(self):
        a=np.array([[2.,-1.,.2],[1.,3.,-1.],[.1,.2,4.]])
        rows=np.array([1e-12,1.,1e12]);cols=np.array([1e4,1e-4,1.])
        a=rows[:,None]*a*cols[None,:]
        expected=np.array([1.,-2.,3.])/rows
        rhs=a.T@expected
        p=PontNoyau(None)
        for creux in (False,True):
            actual=p._resout_adjoint(csc_matrix(a) if creux else a,rhs,creux)
            np.testing.assert_allclose(actual,expected,rtol=2e-14,atol=0.)

    def test_resolution_adjointe_singuliere_et_non_finie_refusee(self):
        p=PontNoyau(None,None)
        for creux in (False,True):
            a=np.array([[1.,1.],[1.,1.]])
            with self.assertRaisesRegex(RuntimeError,'non inversible'):
                p._resout_adjoint(csc_matrix(a) if creux else a,np.ones(2),creux)
            a[0,0]=float('nan')
            with self.assertRaises(ValueError):
                p._resout_adjoint(csc_matrix(a) if creux else a,np.ones(2),creux)

    def test_produit_sans_liaison(self):
        def model(p):
            n=Noyau([0,0,0])
            n.corps('libre',1.,np.diag([1.,2.,3.]).ravel().tolist(),[0,0,0],w=[.1,.2,.3])
            return n
        pont=PontNoyau(model,lambda N,p:np.zeros((6,1)),creux=True)
        N,fr,sch=pont.aller([1.],.02,.01)
        args=(N,fr[0],fr[1],np.asarray(sch[0][0]),np.asarray(sch[0][1]),np.asarray(sch[1][0]),np.asarray(sch[1][1]),[1.])
        aa,bb=pont._pas(*args)
        mu=np.arange(24,dtype=float)
        previous,grad=pont._recul(*args,mu)
        np.testing.assert_allclose(previous,aa.T@mu,rtol=1e-12,atol=1e-12)
        np.testing.assert_array_equal(grad,bb.T@mu)


if __name__=='__main__':
    unittest.main()
