"""Identités exactes et refus de faux certificats sur trajectoires synthétiques."""
import copy
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import tempfile
import unittest

from algebre_moments_em import eye,inv,mm,ma,ms,mv,skew
from defauts_moments_em import cayley,valider_trace
from attribue_moments_em import attribuer
from attribue_rotation_em import RotationVitesse,transpose,flat


class AttributionRotation(unittest.TestCase):
    def packet(self,shift=0.):
        identity=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
        r=[[float(v) for v in row] for row in cayley([F(0),F(.001)/2,F(0)])]
        p=[shift,1.,0.]
        w=[float(v) for v in mv([[F(v) for v in row] for row in r],[F(shift),F(1,2),F(0)])]
        return dict(cadre='initial',duree=.001,pas_demande=.001,monde=identity,matiere=identity,
            trace=dict(schema='vinkulum.trace.moments_em.1',inerties=[[1.,0.,0.,0.,2.,0.,0.,0.,3.]],
            echantillons=[[0.,[[0.,1.,0.]],[flat(identity)],[[0.,.5,0.]]],
                          [.001,[p],[flat(r)],[w]]]))

    def run_pair(self,b):
        b=copy.deepcopy(b);b['cadre']='materiel'
        with tempfile.TemporaryDirectory() as td:
            paths=[Path(td)/'a.gz',Path(td)/'b.gz']
            for path,data in zip(paths,(self.packet(),b)):
                path.write_bytes(gzip.compress(json.dumps(data).encode()))
            return attribuer(*paths,observateur=RotationVitesse())

    def test_resolvante_et_similarite_non_orthogonale(self):
        x=skew([F(1,3),F(-2,7),F(3,5)])
        c=[[F(2),F(1,7),F(0)],[F(0),F(3),F(1,11)],[F(0),F(0),F(4)]]
        ct=transpose(c);cit=inv(ct)
        y0=skew([F(-2,3),F(4,9),F(1,5)])
        ca=cayley([F(1,3),F(-2,7),F(3,5)])
        self.assertEqual(ca,mm(inv(ma(eye(),ms(F(-1,2),x))),ma(eye(),ms(F(1,2),x))))
        cb=mm(mm(cit,cayley([F(-2,3),F(4,9),F(1,5)])),ct)
        y=mm(mm(cit,y0),ct)
        self.assertEqual(ma(cb,ms(-1,ca)),mm(mm(ms(F(1,2),ma(cb,eye())),ma(y,ms(-1,x))),ms(F(1,2),ma(ca,eye()))))

    def test_petite_perturbation_encadree(self):
        result=self.run_pair(self.packet(2.**-60))
        self.assertTrue(result['attribution_rotation_vitesse_verifiee'])
        self.assertFalse(result['rotation_vitesse']['borne_erreur_ode_continue'])

    def test_grands_defauts_refuses(self):
        for component in (2,3):
            with self.subTest(component=component):
                b=self.packet();b['trace']['echantillons'][1][component][0][0]+=1e-6
                result=self.run_pair(b)
                self.assertTrue(result['attribution_moment_verifiee'])
                self.assertFalse(result['attribution_rotation_vitesse_verifiee'])

    def test_reflexion_refusee_malgre_orthogonalite(self):
        b=self.packet();b['monde'][0][0]=-1.
        # Le constructeur synthétique partage la matrice monde/matière :
        # les deux réflexions restent orthogonales mais hors SO(3).
        result=self.run_pair(b)
        self.assertTrue(result['attribution_moment_verifiee'])
        self.assertFalse(result['rotation_vitesse']['determinants_positifs'])
        self.assertFalse(result['attribution_rotation_vitesse_verifiee'])

    def test_traces_malformees_refusees(self):
        cases=[]
        b=self.packet();b['trace']['echantillons'].pop();cases.append(b)
        b=self.packet();b['trace']['echantillons'][1][1][0].append(0.);cases.append(b)
        b=self.packet();b['trace']['echantillons'][1][3][0][0]=float('nan');cases.append(b)
        b=self.packet();b['trace']['inerties'][0][0]=-1.;cases.append(b)
        b=self.packet();b['pas_demande']=.0001;cases.append(b)
        b=self.packet();b['trace']['echantillons'][1][1][0][0]=True;cases.append(b)
        for i,b in enumerate(cases):
            with self.subTest(case=i),self.assertRaises(ValueError):valider_trace(b)


if __name__=='__main__':unittest.main()
