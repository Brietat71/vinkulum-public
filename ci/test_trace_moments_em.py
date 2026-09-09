"""La trace lit le moment porté sans modifier le modèle source."""
import unittest
import numpy as np
from vinkulum import Noyau


class TraceMoments(unittest.TestCase):
    def modele(self):
        n=Noyau([0.,0.,0.])
        n.corps('toupie',.7,[.001,0,0,0,.002,0,0,0,.003],[0,0,0],w=[.05,20,.03])
        return n

    def test_initial_moment_et_sorties_identiques(self):
        n=self.modele();initial=n.etat_precis()
        d=n._trace_em(.02,.001)
        self.assertEqual(initial,n.etat_precis())
        self.assertFalse(d['propagation_certifiee'])
        self.assertEqual(d['echantillons'][0][1][0],[.001*.05,.002*20,.003*.03])
        rows=n.simule_em(.02,.001)
        self.assertEqual(len(d['echantillons']),len(rows)+1)
        for trace,row in zip(d['echantillons'][1:],rows,strict=True):
            self.assertEqual(trace[0],row[0])
            np.testing.assert_array_equal(np.asarray(trace[2]).view(np.uint64),np.asarray(row[2]).view(np.uint64))
            np.testing.assert_array_equal(np.asarray(trace[3]).view(np.uint64),np.asarray(row[4]).view(np.uint64))

    def test_refus_sans_mutation(self):
        for t,h in ((1.,0.),(float('nan'),.001),(.1,float('inf'))):
            n=self.modele();initial=n.etat_precis()
            with self.assertRaises(ValueError):n._trace_em(t,h)
            self.assertEqual(initial,n.etat_precis())
        n=self.modele();n.liaison('enc',None,0);initial=n.etat_precis()
        with self.assertRaises(RuntimeError):n._trace_em(.1,.001)
        self.assertEqual(initial,n.etat_precis())

    def test_etat_initial_seul(self):
        n=self.modele();d=n._trace_em(0.,.001)
        self.assertEqual(len(d['echantillons']),1)

    def test_refus_apres_modification_de_la_copie(self):
        n=self.modele()
        # Le premier corps est avancé avant l'échec arithmétique du second.
        n.corps('debordement',1.,[1.,0.,0.,0.,2.,0.,0.,0.,3.],
                [0.,0.,0.],w=[1e200,2e200,3e200])
        initial=n.etat_precis()
        with self.assertRaises(RuntimeError):n._trace_em(.001,.001)
        self.assertEqual(initial,n.etat_precis())

    def test_trace_depuis_un_etat_deja_avance(self):
        n=self.modele();n.simule_em(.01,.001)
        initial=n.etat_precis();d=n._trace_em(.02,.001)
        self.assertEqual(initial,n.etat_precis())
        self.assertEqual(d['echantillons'][0][0],.01)
        rows=n.simule_em(.02,.001)
        for a,b in zip(d['echantillons'][1:],rows,strict=True):
            self.assertEqual(a[0],b[0])
            np.testing.assert_array_equal(np.asarray(a[2]).view(np.uint64),np.asarray(b[2]).view(np.uint64))
            np.testing.assert_array_equal(np.asarray(a[3]).view(np.uint64),np.asarray(b[4]).view(np.uint64))


if __name__=='__main__':unittest.main()
