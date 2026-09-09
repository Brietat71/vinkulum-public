"""Enclosures dirigées et séparation des causes sur une trace synthétique."""
from fractions import Fraction as F
import gzip
import itertools
import json
from pathlib import Path
import tempfile
import unittest

from algebre_moments_em import mv
from attribue_moments_em import imv,plancher,plafond,attribuer


class AttributionMoments(unittest.TestCase):
    def test_arrondis_et_coins(self):
        for v in (F(1,3),F(-1,7),F(1,2**200),F(-1,2**200),F(0)):
            self.assertLessEqual(plancher(v),v)
            self.assertGreaterEqual(plafond(v),v)
        a=[[F(1,3),F(-2,7),F(5)],[F(-3),F(4),F(-1,11)],[F(0),F(1),F(2)]]
        box=[(F(-1,3),F(2,7)),(F(-2),F(3)),(F(1,13),F(1,7))]
        bounds=imv(a,box)
        for vertex in itertools.product(*box):
            for (l,h),value in zip(bounds,mv(a,vertex)):
                self.assertLessEqual(l,value);self.assertGreaterEqual(h,value)

    def cas(self,perturbation):
        identity=[1.,0.,0.,0.,1.,0.,0.,0.,1.]
        def packet(cadre,shift):
            return dict(cadre=cadre,duree=.001,pas_demande=.001,monde=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],matiere=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],
                trace=dict(schema='vinkulum.trace.moments_em.1',inerties=[[1.,0.,0.,0.,2.,0.,0.,0.,3.]],
                echantillons=[[0.,[[0.,1.,0.]],[identity],[[0.,.5,0.]]],
                              [.001,[[shift,1.,0.]],[identity],[[0.,.5,0.]]]]))
        with tempfile.TemporaryDirectory() as td:
            a,b=Path(td)/'a.gz',Path(td)/'b.gz'
            a.write_bytes(gzip.compress(json.dumps(packet('initial',0.)).encode()))
            b.write_bytes(gzip.compress(json.dumps(packet('spatial',perturbation)).encode()))
            return attribuer(a,b)

    def test_defaut_seul_et_recomposition(self):
        perturbation=2.**-60
        d=self.cas(perturbation)
        self.assertTrue(d['attribution_moment_verifiee'])
        self.assertEqual(d['maxima_contributions'][:2],[F(0),F(0)])
        for (lo,hi),x in zip(d['contributions_finales'][2],[F(perturbation),F(0),F(0)]):
            self.assertLessEqual(lo,x);self.assertGreaterEqual(hi,x)

    def test_identite_ne_suffit_pas_si_defaut_trop_grand(self):
        d=self.cas(1e-6)
        self.assertTrue(d['complet'])
        self.assertFalse(d['attribution_moment_verifiee'])
        self.assertGreater(max(d['defauts_locaux_max']),d['budget_local'])


if __name__=='__main__':unittest.main()
