"""Oracles rationnels des opérations compilées et des congruences complètes."""
import ctypes as ct
from fractions import Fraction as F
import math
import subprocess
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.sparse import csr_matrix

from inertie_binaire_compile import BibliothequeInertie, compiler, certifier_inertie_compilee
from inertie_complement_binaire import certifier_inertie_binaire
from inertie_complement_dirigee import InertieImpossible
from intervalles_binary64 import IntervallesBinaires
from ordre_separateurs import ordre_separateurs
from matrices_certificat import matrice_entree
import test_inertie_complement_dirigee as temoin
from test_inertie_complement_dirigee import kkt_exact, rationnel
from test_inertie_complement import inertie_fraction


class InertieBinaireCompilee(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='vinkulum-test-inertie-binaire-')
        cls.bib = BibliothequeInertie(compiler(Path(cls.temp.name)/'build'))
        cls.op = cls.bib.lib.vinkulum_interval_binary64
        cls.op.argtypes = [ct.c_int]+[ct.c_double]*4+[ct.POINTER(ct.c_double)]
        cls.op.restype = ct.c_int

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def contient(self, intervalle, exact):
        self.assertLessEqual(F(float.fromhex(intervalle[0])),exact)
        self.assertGreaterEqual(F(float.fromhex(intervalle[1])),exact)

    def preuve_exacte(self,a,preuve):
        temoin.InertieComplementDirigee.preuve_exacte(self,a,preuve)

    def certificat(self,d,m,b,gamma,**kw):
        r=certifier_inertie_compilee(csr_matrix(d),csr_matrix(m),b,gamma,bibliotheque=self.bib,**kw)
        self.assertTrue(r['certification_machine'])
        self.assertEqual(r['encodage_intervalles'],'hex_binary64')
        self.preuve_exacte(kkt_exact(d,m,b,gamma),r['preuve_kkt'])
        if r['preuve_masse']['methode'] != 'diagonale_binary64_positive':
            self.preuve_exacte(rationnel(m),r['preuve_masse'])
        return r

    def test_operations_quadrants_compensees_et_subnormales(self):
        points=[0.,1.,-1.,2.**-1074,-2.**-1074,2.**-1022,2.**-537,2.**-300,2.**300,
                np.nextafter(1.,2.),np.nextafter(1.,0.),1.2345678901234567,-3.456789012345678]
        intervals=[(float(x),float(x)) for x in points]+[(-3.,2.),(-3.,-1.),(1.,3.),(0.,2.),(-2.,0.),(-2.**-1074,2.**-1074)]
        ar=IntervallesBinaires(1000000);operations=(ar.add,ar.sub,ar.mul,ar.div,lambda a,b:ar.carre(a))
        accepted=0
        for a in intervals:
            for b in intervals:
                for op,method in enumerate(operations):
                    out=(ct.c_double*2)();status=self.op(op,*a,*b,out)
                    if op==3 and b[0]<=0<=b[1]:
                        self.assertEqual(status,4);continue
                    aa,bb=tuple(map(F,a)),tuple(map(F,b))
                    if op==0: exacts=[x+y for x in aa for y in bb]
                    elif op==1:exacts=[x-y for x in aa for y in bb]
                    elif op==2:exacts=[x*y for x in aa for y in bb]
                    elif op==3:exacts=[x/y for x in aa for y in bb]
                    else:exacts=[x*x for x in aa]+([F(0)] if aa[0]<=0<=aa[1] else [])
                    if status:
                        self.assertEqual(status,3)
                        with self.assertRaises(ArithmeticError):method(a,b)
                    else:
                        accepted+=1
                        self.assertLessEqual(F(out[0]),min(exacts))
                        self.assertGreaterEqual(F(out[1]),max(exacts))
                        self.assertEqual(tuple(out),method(a,b))
        self.assertGreater(accepted,1000)

    def test_pivot_double_a_singulier(self):
        d,m,b=np.diag([1.,2.,3.]),np.eye(3),np.eye(3)[:,:1]
        r=self.certificat(d,m,b,1.)
        self.assertTrue(any(len(p['indices'])==2 for p in r['preuve_kkt']['pivots']))

    def test_masses_couplees_obliques_et_permutations(self):
        d,m,b=temoin.InertieComplementDirigee.donnees_obliques(self)
        for p in (None,np.array([2,0,3,1])):
            r=self.certificat(d,m,b,.125,permutation=p)
            python=certifier_inertie_binaire(csr_matrix(d),csr_matrix(m),b,.125,permutation=p)
            self.assertEqual(r['operations_binary64'],python['operations_binary64'])
            self.assertEqual(r['preuve_kkt']['signature'],python['preuve_kkt']['signature'])
            for actual,expected in zip(r['preuve_kkt']['pivots'],python['preuve_kkt']['pivots'],strict=True):
                for key in ('indices','signature'):self.assertEqual(actual[key],expected[key])
                for key in set(actual)-{'indices','signature'}:
                    self.assertEqual([float.fromhex(x) for x in actual[key]],[float(x) for x in expected[key]])
        n=7;l=np.eye(n)+np.diag(np.full(n-1,.25),1);m=l.T@l
        d=np.diag(np.arange(3.,n+3.))+np.diag(np.full(n-1,-.125),1)
        b=np.array([[1.,0.],[0.,1.],[1.,-1.],[2.,1.],[-1.,2.],[3.,-2.],[1.,1.]])
        self.certificat(d,m,b,.5)

    def test_echantillon_rationnel_signatures_et_refutations(self):
        rng=np.random.default_rng(651932);acceptes=refutes=0
        for n in (3,4,6):
            for repetition in range(8):
                d=rng.integers(-2,3,size=(n+2,n)).astype(float)+np.vstack((4*np.eye(n),np.zeros((2,n))))
                l=np.eye(n)+np.diag(rng.integers(-2,3,size=n-1)/8,1);m=l.T@l
                b=rng.integers(-2,3,size=(n,1)).astype(float);b[0,0]=1.
                for gamma in (.125,8.,64.):
                    exact=kkt_exact(d,m,b,gamma);expected=inertie_fraction(exact)
                    if expected==(n,1,0):
                        self.certificat(d,m,b,gamma);acceptes+=1
                    else:
                        with self.assertRaises(InertieImpossible) as cm:
                            certifier_inertie_compilee(csr_matrix(d),csr_matrix(m),b,gamma,bibliotheque=self.bib)
                        self.assertEqual(cm.exception.bilan['code'],7)
                        self.preuve_exacte(exact,cm.exception.bilan['preuve_kkt']);refutes+=1
        self.assertGreater(acceptes,24);self.assertGreater(refutes,24)

    def test_rang_deficient_frontiere_masse_indefinie(self):
        d=np.diag([1.,2.,3.,4.]);m=np.eye(4)
        for b,gamma in ((np.zeros((4,1)),.5),(np.array([[1.,2.],[0.,0.],[0.,0.],[0.,0.]]),.5),(np.eye(4)[:,:1],4.)):
            with self.assertRaises(InertieImpossible) as cm:
                certifier_inertie_compilee(csr_matrix(d),csr_matrix(m),b,gamma,bibliotheque=self.bib)
            self.assertEqual(cm.exception.bilan['code'],6)
        m[:2,:2]=[[1.,2.],[2.,1.]]
        with self.assertRaises(InertieImpossible) as cm:
            certifier_inertie_compilee(csr_matrix(d),csr_matrix(m),np.eye(4)[:,:1],.5,bibliotheque=self.bib)
        self.assertEqual(cm.exception.bilan['code'],5)

    def test_budgets_et_entrees_invalides(self):
        d,m,b=temoin.InertieComplementDirigee.donnees_obliques(self)
        for kw in (dict(budget_operations=1),dict(budget_coefficients=6),dict(budget_rectangulaire=1)):
            with self.assertRaises(InertieImpossible):
                certifier_inertie_compilee(csr_matrix(d),csr_matrix(m),b,.125,bibliotheque=self.bib,**kw)
        for kw in (dict(permutation=[0,0,2,3]),dict(permutation=[0.,1.,2.,3.]),dict(budget_operations=2**64)):
            with self.assertRaises(ValueError):
                certifier_inertie_compilee(csr_matrix(d),csr_matrix(m),b,.125,bibliotheque=self.bib,**kw)
        for value in (float('nan'),float('inf'),0.,-1.,1j):
            with self.assertRaises(ValueError):
                certifier_inertie_compilee(csr_matrix(d),csr_matrix(m),b,value,bibliotheque=self.bib)

    def test_reechelonnement_et_refus_extremes(self):
        d,m,b=temoin.InertieComplementDirigee.donnees_obliques(self)
        for exponent in (-100,100):self.certificat(d,m,b*2.**exponent,.125)
        # Une congruence exacte ne garantit pas la réussite d'une arithmétique
        # de largeur finie : ici un carré intermédiaire sous-flue ou déborde.
        for exponent in (-300,300):
            with self.assertRaises(InertieImpossible) as cm:
                certifier_inertie_compilee(csr_matrix(d),csr_matrix(m),b*2.**exponent,.125,bibliotheque=self.bib)
            self.assertIn(cm.exception.bilan['code'],(3,6))
        # Une donnée non nulle ne peut être supprimée parce que son carré sous-flue.
        tiny=np.array([[2.**-1000],[0.],[0.]])
        with self.assertRaises(InertieImpossible) as cm:
            certifier_inertie_compilee(csr_matrix(np.diag([2.**-50,1.,2.])),csr_matrix(np.eye(3)),tiny,.5,bibliotheque=self.bib)
        self.assertIn(cm.exception.bilan['code'],(3,6))
        with self.assertRaises(InertieImpossible) as cm:
            certifier_inertie_compilee(csr_matrix(np.diag([2.**600,1.,2.])),csr_matrix(np.eye(3)),np.eye(3)[:,:1],.5,bibliotheque=self.bib)
        self.assertEqual(cm.exception.bilan['code'],3)

    def test_ordre_generique_motifs_deconnectes_et_masses(self):
        for n in (0,1,9,80,2000):
            d=csr_matrix(np.ones((1,n))) if n<100 else csr_matrix((0,n))
            m=csr_matrix((n,n));p=ordre_separateurs(d,m)
            np.testing.assert_array_equal(np.sort(p),np.arange(n))
            np.testing.assert_array_equal(p,ordre_separateurs(d,m))
        n=31;d=np.eye(n)+np.diag(np.full(n-1,-1.),1);m=np.eye(n)
        p=ordre_separateurs(csr_matrix(d),csr_matrix(m));self.assertFalse(np.array_equal(p,np.arange(n)))
        self.certificat(3*d,m,np.eye(n)[:,:1],.0001,permutation=p)

    def test_csr_frontieres_vides_et_drapeaux_mensongers(self):
        x=csr_matrix(([1.,2.,3.,0.],[0,2,1,3],[0,0,2,2,4,4]),shape=(5,4))
        original=x.copy();actual=matrice_entree(x,'D')
        np.testing.assert_array_equal(actual.toarray(),original.toarray())
        self.assertEqual(x.nnz,4);self.assertEqual(actual.nnz,3)
        for ids in ([2,0,1,3],[0,0,1,3]):
            bad=x.copy();bad.indices[:]=ids
            bad.has_sorted_indices=bad.has_canonical_format=True
            with self.assertRaisesRegex(ValueError,'doublons ou indices'):
                matrice_entree(bad,'D')
        for ptr in ([0,0,4,2,4,4],[0,0,2,2,4,5]):
            bad=x.copy();bad.indptr[:]=ptr
            with self.assertRaisesRegex(ValueError,'stockage CSR'):
                matrice_entree(bad,'D')
        for n in (0,1,3):
            self.assertEqual(matrice_entree(csr_matrix((n,4)),'vide').nnz,0)
        with self.assertRaises(ValueError):matrice_entree(x.astype(complex),'complexe')
        bad=x.copy();bad.data[0]=np.inf
        with self.assertRaises(ValueError):matrice_entree(bad,'infini')

    def test_gram_arrondi_ne_remplace_pas_d_exact(self):
        a=1.1;d=np.array([[a,a+3*np.spacing(a)]])
        m,b,gamma=np.eye(2),np.ones((2,1)),2.**-54
        self.assertEqual(inertie_fraction(kkt_exact(d,m,b,gamma)),(1,2,0))
        with self.assertRaises(InertieImpossible) as cm:
            certifier_inertie_compilee(csr_matrix(d),csr_matrix(m),b,gamma,bibliotheque=self.bib)
        self.assertIn(cm.exception.bilan['code'],(6,7))

    def test_modes_arrondis_et_ftz_daz_refuses_sans_mutation(self):
        # Processus séparé : aucun changement du mode de NumPy dans le testeur.
        directory=Path(self.temp.name);source=directory/'environnement.cpp'
        source.write_text(r'''
#include <cfenv>
#include <cstdio>
#ifdef __SSE__
#include <xmmintrin.h>
#endif
extern "C" int vinkulum_interval_binary64(int,double,double,double,double,double*);
int main() {
    double out[2];int initial=std::fegetround();
    for(int mode:{FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}) {
        if(std::fesetround(mode))return 1;
        if(vinkulum_interval_binary64(0,1,1,1,1,out)!=8||std::fegetround()!=mode)return 2;
    }
    if(std::fesetround(initial))return 3;
#ifdef __SSE__
    unsigned original=_mm_getcsr();
    for(unsigned bits:{64u,32768u,32832u}) {
        unsigned changed=original|bits;_mm_setcsr(changed);
        int status=vinkulum_interval_binary64(0,1,1,1,1,out);
        unsigned after=_mm_getcsr();_mm_setcsr(original);
        // Les drapeaux d'exceptions peuvent changer, pas les bits de contrôle.
        if(status!=8||(after&~63u)!=(changed&~63u))return 4;
    }
#endif
    return vinkulum_interval_binary64(0,1,1,1,1,out);
}
'''.replace('#include <cstdio>','#include <cstdio>\n#include <initializer_list>'))
        exe=directory/'environnement'
        command=['c++','-std=c++17',str(source),self.bib.identite['chemin'],'-o',str(exe)]
        result=subprocess.run(command,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        result=subprocess.run([str(exe)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)


if __name__=='__main__':unittest.main()
