"""Références physiques et transformations qui doivent commuter avec le calcul."""
import unittest

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

from vinkulum import Noyau


J = np.diag([.001, .002, .003])
W = np.array([.05, 20., .03])
ROT = Rotation.from_rotvec([.43, -1.11, 2.07]).as_matrix()
CADRES = (np.eye(3), Rotation.from_rotvec([1.82, -.47, .15]).as_matrix())
ECHELLES = (1e-20, 1e-16, 1e-12, 1e-8, 1., 1e8, 1e16, 1e20)


def toupie(echelle, cadre, pas, w=W):
    n = Noyau([0., 0., 0.])
    n.corps('toupie', .7, (echelle*cadre.T@J@cadre).ravel().tolist(), [0., 0., 0.],
            rot=(ROT@cadre).ravel().tolist(), w=(ROT@w).tolist())
    n.simule_em(.1, pas, tous=10**9)
    etat = n.etat()
    r = np.array(etat[2][0]).reshape(3, 3)@cadre.T
    omega = np.array(etat[4][0])
    return r, omega


class FiabiliteLisse(unittest.TestCase):
    def test_toupie_changement_echelle_inertie_et_repere(self):
        # Multiplier J par s ne change pas l'équation d'Euler libre.
        # Référence indépendante : aucune matrice exportée par le noyau.
        rhs = lambda t, w: np.linalg.solve(J, np.cross(J@w, w))
        refs = [solve_ivp(rhs, (0., .1), W, method='DOP853', rtol=tol, atol=tol/100)
                for tol in (1e-11, 1e-13)]
        self.assertTrue(all(s.success and s.t[-1] == .1 for s in refs))
        np.testing.assert_allclose(refs[0].y[:, -1], refs[1].y[:, -1], rtol=0., atol=1e-11)
        reference = refs[1].y[:, -1]
        r0, w0 = toupie(1., np.eye(3), .001)
        energie0 = .5*W@J@W
        moment0 = ROT@J@W
        for scale in ECHELLES:
            for cadre in CADRES:
                with self.subTest(echelle=scale, cadre=cadre.tolist()):
                    erreurs = []
                    for h in (.002, .001):
                        r, w = toupie(scale, cadre, h)
                        wb = r.T@w
                        erreurs.append(np.linalg.norm(wb-reference))
                        self.assertLess(abs(.5*wb@J@wb-energie0)/energie0, 1e-11)
                        self.assertLess(np.linalg.norm(r@J@wb-moment0)/np.linalg.norm(moment0), 1e-11)
                        if h == .001:
                            np.testing.assert_allclose(r, r0, rtol=0., atol=1e-11)
                            np.testing.assert_allclose(w, w0, rtol=0., atol=1e-10)
                    self.assertLess(erreurs[-1], 1e-6)
                    self.assertGreater(erreurs[0]/erreurs[1], 3.8)
                    self.assertLess(erreurs[0]/erreurs[1], 4.2)

    def test_moment_nul_ne_divise_pas_par_zero(self):
        for scale in (1e-20, 1., 1e20):
            r, w = toupie(scale, np.eye(3), .001, w=np.zeros(3))
            np.testing.assert_allclose(r, ROT, rtol=0., atol=1e-14)
            np.testing.assert_array_equal(w, np.zeros(3))


if __name__ == '__main__':
    unittest.main()
