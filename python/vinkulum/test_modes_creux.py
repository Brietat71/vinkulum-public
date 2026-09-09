"""Contre-épreuves des modes signés et des contraintes en métrique de masse."""
import math
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import eigh,null_space,subspace_angles
from scipy.sparse import csc_matrix
from scipy.spatial.transform import Rotation

from vinkulum import Noyau
from vinkulum.modes_creux import analyse,_Contraintes,_ritz_physique
from vinkulum.test_analyses import oscillateur
from vinkulum.test_noyau import _cascade_initiale


def console(ne,unites=(1.,1.,1.),rot=None,libre=False):
    length,mass,temps=unites
    rot=np.eye(3) if rot is None else rot
    n=Noyau([0,0,0]);dl=1/ne
    for i in range(ne+1):
        j=np.diag([.002,.001,.001])*mass*length**2
        n.corps(str(i),dl*mass,j.ravel().tolist(),(rot@[i*dl*length,0,0]).tolist(),rot=rot.ravel().tolist())
    if not libre:n.liaison('root',None,0)
    for i in range(ne):
        n.poutre(str(i),i,i+1,2e5*mass*length/temps**2,8e4*mass*length/temps**2,
                 50.*mass*length**3/temps**2,100.*mass*length**3/temps**2)
    return n


class ModesCreux(unittest.TestCase):
    def test_ritz_physique_sur_petite_valeur_analytique(self):
        # K = H[[1,-1],[-1,1]] + delta I : (1,1) a exactement
        # la valeur propre delta. Les grands termes s'annulent.
        h,delta=2.**20,2.**-10
        k=csc_matrix([[h+delta,-h],[-h,h+delta]])
        x=np.full((2,1),1/np.sqrt(2.))
        v,_,bits=_ritz_physique(k,csc_matrix(np.eye(2)),x)
        limite=1e-8 if bits>53 else 1e-6
        self.assertLess(abs(v[0]/delta-1),limite)

    def test_masse_couplee_refusee(self):
        # Export synthétique : empêche une future extension de M de rendre
        # silencieusement fausse la normalisation par blocs des corps.
        def export(a):
            a=csc_matrix(a)
            return (*a.shape,a.indptr.tolist(),a.indices.tolist(),a.data.tolist())
        class Modele:
            def _verifie_domaine_adjoint(self):pass
            def k_c_m_g_creux(self,t):
                m=np.eye(12);m[0,6]=m[6,0]=.5
                return (export(np.zeros((12,12))),export(np.zeros((12,12))),export(m),export(np.zeros((0,12))))
        with self.assertRaisesRegex(ValueError,'masse couplée'):analyse(Modele(),12)

    def test_rang_numerique_distinct_du_rang_exact(self):
        eps=1e-18
        b=csc_matrix([[1.,0.,0.],[1.,eps,0.]])
        p=_Contraintes(b)
        self.assertEqual(p.rank,1)
        self.assertGreater(p.defaut_reduction,0.)
        # Le mineur des deux premières colonnes vaut eps != 0 exactement.
        # Ce test conserve la limitation : le rang numérique n'est pas
        # une preuve de rang du tableau binary64 interprété exactement.
        n=_cascade_initiale(2,np.eye(3));r=analyse(n,2)
        self.assertFalse(r['rang_certifie'])
        self.assertTrue(r['reduction_contraintes_numerique'])
        self.assertFalse(r['certification_spectrale'])

    def test_refus_decalage_singulier_et_etat_conserve(self):
        n=oscillateur(4.);avant=n.etat_precis()
        with self.assertRaisesRegex(RuntimeError,'singulier'):n.modes_creux(1,decalage=4.)
        self.assertEqual(n.etat_precis(),avant)
        r=n.modes_creux(0)
        self.assertIsNone(r['rang_contraintes'])
        self.assertIsNone(r['ddl_admissibles'])

    def check_original(self,n,r):
        k,_,m,_,g=[np.asarray(x) for x in n.k_c_m_z()]
        g=g.reshape(-1,m.shape[0])
        x=np.array([m['forme'] for m in r['modes']]).T
        if x.size:
            np.testing.assert_allclose(x.T@m@x,np.eye(x.shape[1]),atol=1e-10)
        for mode in r['modes']:
            v=np.array(mode['forme']);eta=np.array(mode['reactions']);lam=mode['valeur_propre']
            scale=np.linalg.norm(k)*np.linalg.norm(v)+np.linalg.norm(g)*np.linalg.norm(eta)+abs(lam)*np.linalg.norm(m)*np.linalg.norm(v)
            err=np.linalg.norm((k+k.T)@v/2+g.T@eta-lam*m@v)
            self.assertLessEqual(err,1e-8*scale+1e-12)
            self.assertLessEqual(np.linalg.norm(g@v),1e-9*max(np.linalg.norm(g)*np.linalg.norm(v),1e-30))

    def test_valeurs_signees_analytiques(self):
        for stiffness in (4.,-4.,0.):
            n=oscillateur(stiffness)
            r=analyse(n,6)
            self.assertEqual(len(r['modes']),1)
            self.assertAlmostEqual(r['modes'][0]['valeur_propre'],stiffness,places=12)
            self.assertEqual(r['modes'][0]['statut'],'positif' if stiffness>0 else 'negatif' if stiffness<0 else 'signe_indetermine')
            self.check_original(n,r)

    def test_probleme_vide_ou_entierement_bloque(self):
        n=Noyau([0,0,0]);self.assertEqual(analyse(n)['modes'],[])
        n.corps('b',1.,np.eye(3).ravel().tolist(),[0,0,0])
        r=analyse(n);self.assertEqual(len(r['modes']),6)
        np.testing.assert_array_equal([m['valeur_propre'] for m in r['modes']],np.zeros(6))
        n.liaison('fixe',None,0)
        self.assertEqual(analyse(n)['modes'],[])

    def test_console_et_multiplicites(self):
        n=console(30);before=n.etat_precis();r=analyse(n,8)
        self.assertEqual(before,n.etat_precis())
        dense=n.modes(8)
        np.testing.assert_allclose([m['frequence_hz'] for m in r['modes']],[f for f,v in dense],rtol=3e-8)
        self.check_original(n,r)
        # Comparer le sous-espace, pas les signes ou les bases des paires doubles.
        x=np.array([v for f,v in dense]).T;y=np.array([m['forme'] for m in r['modes']]).T
        self.assertLess(np.max(subspace_angles(x,y)),1e-7)

    def test_modes_rigides_et_flexibles_libres(self):
        n=console(8,libre=True);r=analyse(n,10)
        vals=np.sort([m['valeur_propre'] for m in r['modes']])
        self.assertLess(np.max(abs(vals[:6])),1e-6)
        self.assertGreater(vals[6],100.)
        self.check_original(n,r)
        rigid=np.zeros((54,6))
        for i in range(9):
            rigid[6*i:6*i+3,:3]=np.eye(3)
            rigid[6*i+3:6*i+6,3:]=np.eye(3)
            for j in range(3):rigid[6*i:6*i+3,3+j]=np.cross(np.eye(3)[j],[i/8,0,0])
        x=np.array([m['forme'] for m in r['modes'][:6]]).T
        self.assertLess(np.max(subspace_angles(x,rigid)),1e-8)

    def test_contraintes_redondantes_generales(self):
        n=_cascade_initiale(4,np.eye(3));before=n.etat_precis()
        k,_,m,_,g=[np.asarray(x) for x in n.k_c_m_z()]
        z=null_space(g);ref=eigh(z.T@(k+k.T)@z/2,z.T@m@z,eigvals_only=True)
        r=analyse(n,10)
        self.assertEqual(n.etat_precis(),before)
        self.assertEqual(r['ddl_admissibles'],4)
        np.testing.assert_allclose(sorted(m['valeur_propre'] for m in r['modes']),ref,rtol=1e-9,atol=1e-8)
        self.check_original(n,r)

    def test_projection_grande_famille_et_reactions(self):
        # Chaîne de différences : référence explicite du noyau, puis une
        # contrainte dépendante qui ne se réduit pas à une ligne dupliquée.
        m,n=80,81
        b=np.zeros((m,n))
        for i in range(m):b[i,i:i+2]=[1.,-1.]
        rng=np.random.default_rng(671);v=rng.normal(size=(n,3))
        for extra in (False,True):
            bb=np.vstack((b,b[0]+b[1])) if extra else b
            p=_Contraintes(csc_matrix(bb))
            self.assertEqual(p.rank,m)
            np.testing.assert_allclose(p.projette(v),np.broadcast_to(v.mean(axis=0),v.shape),atol=5e-13)
            q=p.original.T@p.reactions(v)
            np.testing.assert_allclose(q,v-v.mean(axis=0),atol=5e-13)

    def test_changements_unites(self):
        ref=analyse(console(8),6)
        reference=np.array([m['valeur_propre'] for m in ref['modes']])
        for units in ((1e-3,1e5,.1),(1e3,1e-5,10.)):
            r=analyse(console(8,units),6)
            np.testing.assert_allclose([m['valeur_propre']*units[2]**2 for m in r['modes']],reference,rtol=1e-9)

    def test_decalage_interieur(self):
        n=console(8);dense=n.modes(20)
        vals=np.array([(2*np.pi*f)**2 for f,v in dense]);shift=50000.
        r=analyse(n,4,decalage=shift)
        expected=vals[np.argsort(abs(vals-shift))[:4]]
        np.testing.assert_allclose([m['valeur_propre'] for m in r['modes']],expected,rtol=1e-9)

    def test_repere_tourne_et_interface_native(self):
        q=Rotation.from_rotvec([.3,-.7,.4]).as_matrix()
        ref=analyse(console(8),8)
        n=console(8,rot=q);r=n.modes_creux(8)
        np.testing.assert_allclose([m['valeur_propre'] for m in r['modes']],
                                   [m['valeur_propre'] for m in ref['modes']],rtol=1e-9)
        self.check_original(n,r)

    def test_domaine_et_parametres_invalides(self):
        n=console(2)
        for kwargs in (dict(combien=-1),dict(combien=1.5),dict(tol=0),dict(tol=float('nan')),dict(t=float('nan')),dict(decalage=float('inf')),dict(maxiter=0)):
            with self.assertRaises(ValueError):analyse(n,**kwargs)
        n.sphere(0,[0,0,0],.1)
        with self.assertRaisesRegex(ValueError,'hors domaine'):analyse(n)


if __name__=='__main__':unittest.main()
