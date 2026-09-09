"""Enveloppes par congruence et conservation du travail aux jonctions."""
from fractions import Fraction as F
import unittest
import numpy as np
from test_trace_complement import matrice,restreindre,difference,echelle,gram,identite
from test_inertie_complement import inertie_fraction
from modeles_assemblage_complements import geometrie
from test_krylov_contraint import rationnel,energie

class PreuvesAssemblage(unittest.TestCase):
    def test_enveloppe_globale_pour_defauts_locaux_indefinis(self):
        # Les erreurs locales n'ont pas besoin d'être positives. Le transport
        # ±bI survit à des applications obliques et aux modes privés distincts.
        maps=[matrice([[1,2,0,0],[0,0,1,0]]),matrice([[2,-1,0,0],[0,0,0,1]])]
        defects=[matrice([[1,1],[1,-1]]),matrice([[-1,0],[0,1]])]
        constants=[F(2),F(1)];global_error=matrice([[0]*4 for _ in range(4)]);envelope=matrice([[0]*4 for _ in range(4)])
        for p,d,b in zip(maps,defects,constants,strict=True):
            for gap in (difference(echelle(identite(2),b),d),difference(echelle(identite(2),b),echelle(d,-1))):
                self.assertEqual(inertie_fraction(gap)[1],0)
            addition=restreindre(d,p);bound=echelle(gram(p),b)
            global_error=difference(global_error,echelle(addition,-1));envelope=difference(envelope,echelle(bound,-1))
        for sign in (-1,1):self.assertEqual(inertie_fraction(difference(envelope,echelle(global_error,sign)))[1],0)
        # Somme de bornes scalaires sans congruence : 3I ne couvre pas
        # cette application, qui amplifie la direction du premier port.
        self.assertGreater(inertie_fraction(difference(echelle(identite(4),3),global_error))[1],0)
    def test_metrique_assemblee_evite_le_facteur_nombre_pieces(self):
        # Trois contributions scalaires au port global et un mode privé chacune.
        # H=sum(a_j² H_j)=9, W=1/3 est rationnel : identité exacte de repère.
        maps=[matrice([[F(1,3),0,0,0],[0,1,0,0]]),
              matrice([[F(2,3),0,0,0],[0,0,1,0]]),
              matrice([[F(2,3),0,0,0],[0,0,0,1]])]
        total=matrice([[0]*4 for _ in range(4)])
        for p in maps:total=difference(total,echelle(gram(p),-1))
        self.assertEqual(total,identite(4))
        bs=[F(1,4),F(1,8),F(1,16)];envelope=matrice([[0]*4 for _ in range(4)])
        for p,b in zip(maps,bs,strict=True):envelope=difference(envelope,echelle(gram(p),-b))
        self.assertEqual(inertie_fraction(difference(echelle(identite(4),max(bs)),envelope))[1],0)
    def test_travail_et_bras_de_levier(self):
        q=np.array([.25,-.5,1.,.125,.25,-.125]);f=np.array([2.,-1.,.5,.25,1.,-.5])
        offsets=[np.array([0.,0.,0.]),np.array([.125,0.,.0625]),np.array([0.,-.125,.0625])]
        for a,offset in zip(geometrie(),offsets,strict=True):
            rotation=a[:3,:3];local=np.r_[rotation@(q[:3]+np.cross(q[3:],offset)),rotation@q[3:]]
            np.testing.assert_array_equal(a@q,local)
            world_force=rotation.T@f[:3];world_moment=rotation.T@f[3:]+np.cross(offset,world_force)
            np.testing.assert_array_equal(a.T@f,np.r_[world_force,world_moment])
            self.assertEqual(sum(F(float(x))*F(float(y)) for x,y in zip(f,local)),sum(F(float(x))*F(float(y)) for x,y in zip(a.T@f,q)))

if __name__=='__main__':unittest.main()
