"""Contre-épreuves du juge spectral commun, sans solveur concurrent."""
import copy
import math
import unittest
from confronte_modes import juger, frequences_mbdyn, modele_mbdyn, commande_mbdyn


def cayley(freqs, real=0.):
    c=.01
    z=[(1+c*complex(real,2*math.pi*f))/(1-c*complex(real,2*math.pi*f)) for f in freqs]
    return dict(dCoef=c,alpha=[[v.real,v.imag,1.] for v in z])


class JugeModal(unittest.TestCase):
    def setUp(self):
        self.ref=dict(frequences=[1.,2.,3.],bande_hz=[.5,3.5])

    def test_transformation_homogene(self):
        r=cayley([1.,2.,3.])
        for factor in (1.,-1.,1e200,1e-200):
            scaled=dict(dCoef=r['dCoef'],alpha=[[x*factor for x in row] for row in r['alpha']])
            self.assertEqual(juger('mbdyn_lapack',scaled,self.ref)['decision'],'qualifie')

    def test_modes_manquants_et_repetes(self):
        for frequencies in ([1.,3.],[1.,2.,2.],[1.,2.,3.,3.]):
            self.assertNotEqual(juger('mbdyn_lapack',cayley(frequencies),self.ref)['decision'],'qualifie')

    def test_imprecision_et_partie_reelle(self):
        self.assertEqual(juger('mbdyn_lapack',cayley([1.,2.,3.000001]),self.ref)['decision'],'imprecis')
        self.assertEqual(juger('mbdyn_lapack',cayley([1.,2.,3.],real=.001),self.ref)['decision'],'imprecis')

    def test_contraintes_et_bande(self):
        r=cayley([1.,2.,3.,100.])
        r['alpha'] += [[0.,0.,0.],[-1.,0.,1.],[1.,0.,1.]]
        d=juger('mbdyn_lapack',r,self.ref)
        self.assertEqual(d['decision'],'qualifie')
        self.assertEqual(d['hors_bande'],1)

    def test_spectre_symetrique_invalide(self):
        for x in (-1.,0.,float('nan'),float('inf')):
            self.assertEqual(juger('vinkulum',dict(valeurs_propres=[x]),self.ref)['decision'],'spectre_invalide')
        r=dict(valeurs_propres=[(2*math.pi*f)**2 for f in self.ref['frequences']])
        self.assertEqual(juger('exudyn_arbre_creux',r,self.ref)['decision'],'qualifie')

    def test_nom_spectral_avec_extension(self):
        cmd=commande_mbdyn('/tmp/une.version/mbdyn')
        self.assertEqual(cmd[-2:],['-o','natif.out'])

    def test_modele_et_configuration_fixes(self):
        s=modele_mbdyn(4,'mbdyn_arpack_large',[.5,3.5])
        self.assertIn('joints: 9;',s)
        self.assertIn('use arpack, 82, 86, 1e-12',s)
        self.assertIn('parameter, 0.02',s)
        self.assertEqual(s.count('linear elastic generic, diag, 0., 40, 0.'),4)


if __name__=='__main__':
    unittest.main()
