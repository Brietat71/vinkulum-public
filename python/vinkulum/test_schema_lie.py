"""Références physiques indépendantes et domaine public de la cinématique sigma."""
import unittest

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

from vinkulum import Noyau


def _libre():
    n = Noyau([0., 0., 0.])
    n.corps('sphere', 1., np.eye(3).ravel().tolist(), [0., 0., 0.], w=[1.1, -.4, .7])
    n.effort(0, [0., 0., 0.], [-.3, .8, .2])
    return n


class SchemaLie(unittest.TestCase):
    def test_rotation_affine_contre_edo(self):
        # Inertie sphérique : omega(t) = omega0 + t*M est EXACTE.
        # L'EDO Rdot=[omega(t)]R ne contient ni Newmark ni Jacobien de Lie.
        def rhs(t, y):
            x, v, z = np.array([1.1, -.4, .7])+t*np.array([-.3, .8, .2])
            s = np.array([[0., -z, v], [z, 0., -x], [-v, x, 0.]])
            return (s @ y.reshape(3, 3)).ravel()
        sol = solve_ivp(rhs, (0., 1.), np.eye(3).ravel(), method='DOP853',
                        rtol=2e-13, atol=2e-15)
        self.assertTrue(sol.success)
        reference = sol.y[:, -1].reshape(3, 3)
        for rho in (0., .5, .9):
            gamma = .5+(1-rho)/(1+rho)
            beta = .25*(gamma+.5)**2
            errors = []
            for h in (.08, .04, .02):
                n = _libre()
                n.simule(1., h, rho=rho, sigma_lie=gamma/(3*beta))
                r = np.array(n.etat()[2][0]).reshape(3, 3)
                errors.append(Rotation.from_matrix(r @ reference.T).magnitude())
                np.testing.assert_allclose(n.etat()[4][0], [.8, .4, .9], atol=2e-12, rtol=0.)
            self.assertLess(errors[-1], 6e-10)
            self.assertGreater(errors[1]/errors[2], 12.)
        n = _libre()
        n.simule(1., .02)
        r = np.array(n.etat()[2][0]).reshape(3, 3)
        self.assertGreater(Rotation.from_matrix(r @ reference.T).magnitude(), 3e-5)

    def test_options_invalides_et_audit_restaure(self):
        n = _libre()
        before = n.etat_precis()
        with self.assertRaisesRegex(ValueError, 'estimateur adaptatif'):
            n.simule(.1, .01, sigma_lie=1., adaptatif=1e-4)
        self.assertEqual(before, n.etat_precis())
        for sigma in (-.1, 1.1, float('nan'), float('inf')):
            with self.assertRaisesRegex(ValueError, 'sigma_lie'):
                n.simule(.1, .01, sigma_lie=sigma)
            with self.assertRaisesRegex(ValueError, 'sigma_lie'):
                n.audit_jacobien(.01, sigma_lie=sigma)
            self.assertEqual(before, n.etat_precis())
        for sigma in (.6, 1.):
            blocks, _ = n.audit_jacobien(.1, sigma_lie=sigma)
            self.assertLess(max(e for _, e, _ in blocks), 1e-7)
            self.assertEqual(before, n.etat_precis())
        a, b = _libre(), _libre()
        a.simule(.2, .01)
        b.simule(.2, .01, sigma_lie=0.)
        self.assertEqual(a.etat_precis(), b.etat_precis())
        # Le réglage d'un appel ne doit pas contaminer l'appel suivant.
        n.simule(.1, .01, sigma_lie=1.)
        a.pose_etat(*n.etat_precis())
        n.simule(.2, .01)
        a.simule(.2, .01, sigma_lie=0.)
        self.assertEqual(n.etat_precis(), a.etat_precis())

    def test_plan_et_ggl(self):
        results = []
        for sigma in (0., .665, 1.):
            n = Noyau()
            n.corps('pendule', .2, (.001*np.eye(3)).ravel().tolist(), [.5, 0., -.8])
            n.liaison('pivot', None, 0, bloque_r=[])
            n.simule(.2, .01, sigma_lie=sigma)
            results.append(n.etat_precis())
        self.assertEqual(results[0], results[1])
        self.assertEqual(results[0], results[2])
        # Rotule spatiale : les colonnes zeta doivent différer des colonnes
        # d'accélération, et le solveur doit tenir G*u à la fin du pas.
        n = Noyau()
        n.corps('solide', 1., np.diag([.1, .2, .3]).ravel().tolist(), [.2, .1, .4],
                w=[.3, -.2, .5], v=[-.13, -.02, .07])
        n.liaison('rotule', None, 0, bloque_r=[])
        for sigma in (.6, 1.):
            n.simule(n.etat()[0]+.1, .01, sigma_lie=sigma, ggl=True)
            self.assertLess(np.linalg.norm(n.phi()), 1e-10)
            self.assertLess(np.linalg.norm(n.phi_dot()), 1e-9)
            blocks, _ = n.audit_jacobien(.03, sigma_lie=sigma, ggl=True)
            self.assertLess(max(e for _, e, magnitude in blocks if magnitude > 1e-7), 1e-6)

    def test_pas_hors_carte_rejoue_et_echec_restaure(self):
        # Rotation d'axe fixe exacte ; h*omega > pi force le rejeu, sans
        # transformer un résidu NaN en convergence ni garder un état d'essai.
        n = Noyau([0., 0., 0.])
        n.corps('rotor', 1., np.eye(3).ravel().tolist(), [0., 0., 0.], w=[0., 0., 10.])
        n.simule(.4, .4, sigma_lie=1.)
        self.assertGreater(n.adapt_stats()[0], 0)
        target = Rotation.from_rotvec([0., 0., 4.]).as_matrix()
        np.testing.assert_allclose(np.array(n.etat()[2][0]).reshape(3, 3), target, atol=2e-14)
        n = Noyau([0., 0., 0.])
        n.corps('trop_rapide', 1., np.eye(3).ravel().tolist(), [0., 0., 0.], w=[0., 0., 1e8])
        before = n.etat_precis()
        with self.assertRaises(RuntimeError):
            n.simule(1., 1., sigma_lie=1.)
        self.assertEqual(before, n.etat_precis())


if __name__ == '__main__':
    unittest.main()
