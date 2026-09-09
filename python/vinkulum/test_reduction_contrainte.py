"""Contre-épreuves de l'API installable, sans import des scripts de recherche."""
from fractions import Fraction as F
import unittest
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.sparse import csr_matrix
from vinkulum import _vinkulum
from vinkulum.reduction_contrainte import ReductionContrainte, AssemblageContraint, CertificatIndisponible


def resoudre_exact(d, m, omega, rhs):
    d = [[F(float(x)) for x in row] for row in d]
    n = len(m)
    a = [[sum((row[i]*row[j] for row in d), F(0))-F(float(omega))**2*F(float(m[i,j]))
          for j in range(n)] for i in range(n)]
    a = [row+[F(float(x)) for x in force] for row, force in zip(a, rhs, strict=True)]
    for i in range(n):
        pivot = next(j for j in range(i, n) if a[j][i])
        a[i], a[pivot] = a[pivot], a[i]
        div = a[i][i]; a[i] = [v/div for v in a[i]]
        for j in range(n):
            if j != i:
                coeff = a[j][i]; a[j] = [x-coeff*y for x,y in zip(a[j], a[i], strict=True)]
    return np.array([[float(v) for v in row[n:]] for row in a])


def donnees():
    n=10
    d=np.diag(np.arange(4., n+4))+np.diag(np.full(n-1,.25),1)+np.diag(np.full(n-2,-.125),2)
    l=np.eye(n)+np.diag(np.full(n-1,.125),1)
    return d,l.T@l,np.arange(n-2),np.arange(n-2,n),np.diag([4.,.25])


def reduction(**options):
    d,m,i,s,h=donnees()
    return ReductionContrainte(csr_matrix(d),csr_matrix(m),i,s,h,.5,gamma=4.,**options)


class ApiContrainte(unittest.TestCase):
    def test_champ_physique_face_au_systeme_rationnel(self):
        d,m,i,s,h=donnees(); r=reduction(blocs=1)
        f=np.array([[1.,0.,1.,2.**-30,0.],[0.,1.,1.,-2.**-30,0.]])
        rhs=np.zeros((len(m),len(f.T)));rhs[s]=f
        for w in (0.,.25,.5):
            result=r.reponses(w,f);exact=resoudre_exact(d,m,w,rhs)
            self.assertEqual(result['statut_controle'],'borne_numerique_disponible')
            self.assertFalse(result['certification_machine'])
            e=result['champ']-exact
            errors={'masse':np.sqrt(np.sum(e*(m@e),axis=0)),'deformation':np.linalg.norm(d@e,axis=0)}
            for name,err in errors.items():
                self.assertTrue(np.all(err<=result['bornes'][name]['absolues']))
        self.assertTrue(r.certificats['complement']['certification_machine'])
        self.assertTrue(r.certificats['masse_complete']['certification_machine'])
        self.assertFalse(r.diagnostic['selection']['certification_machine'])

    def test_basses_directions_traversent_une_resonance_interieure(self):
        # Pôle Dirichlet à omega=1, mais le système complet reste inversible.
        d=np.array([[1.,0.,.5],[0.,4.,.25],[0.,0.,2.]])
        m=np.eye(3);f=np.array([[1.]])
        r=ReductionContrainte(d,m,[0,1],[2],[[1.]],1.25,gamma=4.,modes_retenus=1)
        ans=r.reponses(1.,f)
        expected=resoudre_exact(d,m,1.,np.array([[0.],[0.],[1.]]))
        np.testing.assert_allclose(ans['champ'],expected,atol=2e-13)
        self.assertGreater(ans['marge'],0)

    def test_assemblage_physique_independant(self):
        d,m,i,s,h=donnees();rs=[reduction(),reduction()]
        maps=[np.array([[1.,.25],[-.5,1.]]),np.array([[.5,-1.],[1.,.5]])]
        de=np.array([[.5,.25],[0.,1.]]);me=np.diag([.125,.25])
        a=AssemblageContraint(rs,maps,np.diag([2.,.5]),facteur_externe=de,masse_externe=me)
        # Assemblage direct des champs physiques, sans base de réduction.
        p=[]
        for j,map_ in enumerate(maps):
            e=np.zeros((len(m),2+2*len(i)));e[s,:2]=map_;e[i,2+j*len(i):2+(j+1)*len(i)]=np.eye(len(i));p.append(e)
        dd=np.vstack([d@e for e in p]+[np.pad(de,((0,0),(0,2*len(i))))])
        mm=sum(e.T@m@e for e in p);mm[:2,:2]+=me
        f=np.array([[1.,0.,1.],[0.,1.,-1.]])
        rhs=np.zeros((len(mm),3));rhs[:2]=f
        for w in (0.,.5):
            ans=a.reponses(w,f);exact=resoudre_exact(dd,mm,w,rhs)
            np.testing.assert_allclose(ans['champ_ports'],exact[:2],atol=1e-10,rtol=1e-8)
            for field,e in zip(ans['champs'],p,strict=True):
                np.testing.assert_allclose(field,e@exact,atol=1e-10,rtol=1e-8)
            self.assertGreater(ans['marge'],0)
            self.assertFalse(ans['certification_machine'])
        self.assertTrue(a.certificat_masse_externe['certification_machine'])

    def test_pole_local_et_global(self):
        r=ReductionContrainte(np.diag([2.,3.,1.]),np.eye(3),[0,1],[2],[[1.]],1.1,gamma=2.,directions_retenues=[[1.],[0.]])
        with self.assertRaises(np.linalg.LinAlgError):r.reponses(1.,[[1.]])
        a=AssemblageContraint([r],[[[1.]]],[[1.]],facteur_externe=[[2.]])
        np.testing.assert_allclose(a.reponses(1.,[[1.]])['champ_ports'],[[.25]],atol=1e-15)
        a=AssemblageContraint([r],[[[1.]]],[[1.]])
        with self.assertRaises(np.linalg.LinAlgError):a.reponses(1.,[[1.]])

    def test_refus_masse_complete_et_ports_sans_masse(self):
        d,m,i,s,h=donnees()
        for bad in (-1.,0.):
            mm=m.copy();mm[-1,-1]=bad
            with self.assertRaises(ValueError):ReductionContrainte(d,mm,i,s,h,.5)
        mm=m.copy();mm[-1]=0.;mm[:,-1]=0.
        r=ReductionContrainte(d,mm,i,s,h,.5)
        self.assertEqual(r.certificats['masse_complete']['indices_sans_masse'],[len(m)-1])
        # Mii est SPD, mais le couplage aux ports rend la masse complète indéfinie.
        mm=m.copy();mm[0,-1]=mm[-1,0]=10.
        with self.assertRaises(CertificatIndisponible):ReductionContrainte(d,mm,i,s,h,.5)
        # Noyau oblique PSD hors du domaine de cette API.
        mm=np.eye(3);mm[:2,:2]=1.
        with self.assertRaises(CertificatIndisponible):ReductionContrainte(np.eye(3),mm,[0,1],[2],[[1.]],.1)

    def test_copies_et_entrees_invalides(self):
        d,m,i,s,h=donnees();r=ReductionContrainte(d,m,i,s,h,.5)
        ref=r.reponses(.25,np.eye(2))['champ'];d[:]=0;m[:]=0;h[:]=0;i[:]=0;s[:]=0
        cert=r.certificats;cert['complement']['lambda_min']=-1
        np.testing.assert_array_equal(ref,r.reponses(.25,np.eye(2))['champ'])
        self.assertGreater(r.certificats['complement']['lambda_min'],0)
        for kw in (dict(modes_retenus=0),dict(modes_retenus=True),dict(modes_retenus=8),
                   dict(budget_operations=0),dict(profondeur=1),dict(budget_qr=False),
                   dict(directions_retenues=np.zeros((8,1)))):
            with self.assertRaises((ValueError,CertificatIndisponible)):reduction(**kw)
        for w,f in ((-.1,np.eye(2)),(.6,np.eye(2)),(.2,[[1j],[0]]),(.2,[[np.nan],[1.]]),(.2,[1.,0.])):
            with self.assertRaises(ValueError):r.reponses(w,f)
        bad=csr_matrix(np.eye(3));bad.indices[0]=9;bad.has_canonical_format=True
        with self.assertRaises(ValueError):ReductionContrainte(bad,np.eye(3),[0,1],[2],[[1.]],.1)

    def test_marge_non_positive_ne_devient_pas_une_preuve(self):
        good=reduction(blocs=1)
        weak=reduction(blocs=1,alpha_masse=2.**-40)
        a=good.reponses(.5,np.eye(2));b=weak.reponses(.5,np.eye(2))
        np.testing.assert_array_equal(a['champ'],b['champ'])
        self.assertGreater(a['marge'],0)
        self.assertLessEqual(b['marge'],0)
        self.assertEqual(b['statut_controle'],'marge_non_positive')
        self.assertIsNone(b['borne_coordonnees'])
        for value in b['bornes'].values():self.assertIsNone(value['absolues'])

    def test_permutation_physique_et_deux_directions_retenues(self):
        d,m,i,s,h=donnees();f=np.array([[1.,.25],[.5,1.]])
        r=ReductionContrainte(d,m,i[::-1],s[::-1],h[::-1,::-1],.5,gamma=4.,modes_retenus=2)
        rhs=np.zeros((len(m),2));rhs[s]=f
        expected=resoudre_exact(d,m,.25,rhs)
        actual=r.reponses(.25,f[::-1])
        np.testing.assert_allclose(actual['champ'],expected,atol=1e-12,rtol=1e-9)
        self.assertEqual(r.diagnostic['directions_retenues'],2)
        self.assertGreater(actual['marge'],0.)

    def test_extraction_noyau_et_copie(self):
        from vinkulum.test_reduction_ports import _console,_donnees
        n=_console();d,m,free,ports,i,s=_donnees(n)
        r=ReductionContrainte.depuis_noyau(n,free,ports,np.eye(len(s)),.05)
        expected=resoudre_exact(d,m,.025,np.vstack((np.zeros((len(i),len(s))),np.eye(len(s)))))
        np.testing.assert_allclose(r.reponses(.025,np.eye(len(s)))['champ'],expected,rtol=1e-9,atol=1e-11)
        self.assertEqual(r.origine['modele'],'partie_materielle_des_poutres')
        np.testing.assert_array_equal(r.coordonnees_physiques,free)
        with self.assertRaisesRegex(ValueError,'non admissibles'):
            ReductionContrainte.depuis_noyau(n,np.arange(d.shape[1]+6),ports,np.eye(len(s)),.05)

    def test_refus_budget_et_seuil_non_demontre(self):
        with self.assertRaises(CertificatIndisponible):reduction(budget_operations=1)
        with self.assertRaises(CertificatIndisponible):reduction(budget_rectangulaire=1)
        d,m,i,s,h=donnees()
        with self.assertRaises(CertificatIndisponible):
            ReductionContrainte(d,m,i,s,h,.5,gamma=10000.)


class FrontiereNative(unittest.TestCase):
    def arguments(self):
        return dict(n=3,rang=1,d=([0,1,2,3],[0,1,2],[1.,2.,3.]),
                    m=([0,1,2,3],[0,1,2],[1.,1.,1.]),b=[1.,0.,0.],gamma=.5,
                    permutation=[0,1,2],budget_operations=100000,budget_coefficients=1000,budget_rectangulaire=10)

    def test_csr_corrompus_et_dimensions(self):
        bad=[([],[],[]),([0,1,2,9],[0,1,2],[1.,2.,3.]),
             ([0,-1,2,3],[0,1,2],[1.,2.,3.]),([0,2,1,3],[0,1,2],[1.,2.,3.]),
             ([0,1,2,3],[0,1,9],[1.,2.,3.]),([0,1,2,3],[0,1,-1],[1.,2.,3.]),
             ([0,2,2,3],[0,0,2],[1.,2.,3.]),([0,2,2,3],[1,0,2],[1.,2.,3.]),
             ([0,1,2,3],[0,1,2],[1.,np.nan,3.]),([0,1,2,3],[0,1,2],[1.,2.])]
        for name in ('d','m'):
            for value in bad:
                with self.subTest(name=name,value=value):
                    args=self.arguments();args[name]=value
                    with self.assertRaises(ValueError):_vinkulum._certificat_inertie_binaire(**args)
        for update in (dict(n=0),dict(rang=0),dict(rang=3),dict(rang=2**63),dict(b=[]),
                       dict(b=[1.,np.inf,0.]),dict(permutation=[0,0,2]),dict(permutation=[]),
                       dict(budget_operations=0),dict(budget_coefficients=3),dict(budget_rectangulaire=2),
                       dict(gamma=np.nan),dict(gamma=-1.),dict(n=2**64)):
            args=self.arguments();args.update(update)
            with self.assertRaises((ValueError,OverflowError)):_vinkulum._certificat_inertie_binaire(**args)
        args=self.arguments();args['m']=([0,2,3,4],[0,1,1,2],[1.,.5,1.,1.])
        with self.assertRaises(ValueError):_vinkulum._certificat_inertie_binaire(**args)

    def test_comparaison_entrees_invalides(self):
        args=dict(n=3,m=self.arguments()['m'],ell=[1.,1.,1.],alpha=.25,beta=2.,
                  permutation=[0,1,2],budget_operations=10000,budget_coefficients=1000)
        for update in (dict(ell=[1.]),dict(ell=[1.,0.,1.]),dict(ell=[1.,np.inf,1.]),
                       dict(alpha=0.),dict(alpha=2.),dict(beta=np.inf),dict(permutation=[0,0,2]),
                       dict(m=([0,1,2,4],[0,1,2],[1.,1.,1.]))):
            with self.assertRaises(ValueError):_vinkulum._certificat_comparaison_masse(**(args|update))

    @unittest.skipUnless(sys.platform == 'linux', 'constantes fenv Linux')
    def test_mode_arrondi_refuse_sans_modification_du_mode(self):
        # Le mode du processus de tests principal reste inchangé.
        script = """
import ctypes
from vinkulum import _vinkulum
lib = ctypes.CDLL(None)
initial = lib.fegetround()
try:
    for mode in (0x400, 0x800, 0xc00):
        assert lib.fesetround(mode) == 0
        status, _ = _vinkulum._intervalle_binaire(0, 1., 1., 2., 2.)
        assert status == 8, status
        assert lib.fegetround() == mode
finally:
    assert lib.fesetround(initial) == 0
"""
        result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_arrondis_face_aux_rationnels_et_parallele(self):
        rng=np.random.default_rng(8173)
        for _ in range(120):
            a,b,c,e=rng.uniform(-10,10,4);alo,ahi=sorted((a,b));blo,bhi=sorted((c,e))
            for op in range(4):
                if op==3 and blo<=0<=bhi:continue
                status,out=_vinkulum._intervalle_binaire(op,alo,ahi,blo,bhi)
                self.assertEqual(status,0)
                fa,fb=map(F,(alo,ahi));fc,fd=map(F,(blo,bhi))
                if op==0:exact=[fa+fc,fb+fd]
                elif op==1:exact=[fa-fd,fb-fc]
                elif op==2:exact=[x*y for x in (fa,fb) for y in (fc,fd)]
                else:exact=[x/y for x in (fa,fb) for y in (fc,fd)]
                self.assertLessEqual(F(out[0]),min(exact));self.assertGreaterEqual(F(out[1]),max(exact))
        with ThreadPoolExecutor(max_workers=4) as pool:
            results=list(pool.map(lambda _: _vinkulum._certificat_inertie_binaire(**self.arguments()),range(24)))
        self.assertTrue(all(r[0]==0 for r in results))
        for r in results:self.assertEqual(r[1:6],results[0][1:6])


if __name__=='__main__':unittest.main()
