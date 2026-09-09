"""Régressions de l'audit du noyau : références physiques et sorties d'erreur."""
import math
import subprocess
import sys
import unittest

import numpy as np
from scipy.linalg import expm

from vinkulum import Noyau

J = np.eye(3).ravel().tolist()


def masses(beta=0.0):
    n = Noyau(g=[0, 0, 0])
    for i in range(2):
        n.corps(str(i), 1, J, [i, 0, 0])
    k = np.zeros((12, 12))
    k[0, 0] = k[6, 6] = 100
    k[0, 6] = k[6, 0] = -100
    n.superelement("ressort", [0, 1], k.ravel().tolist(), beta=beta)
    t, r, rot, v, w, vi = n.etat()
    r[1][0] += 0.1
    n.pose_etat(t, r, rot, v, w, vi)
    return n


def _pale_inflow_spatial(kind):
    n = Noyau([0., 0., 0.])
    n.corps('pale', 2., np.diag([.2, .3, .4]).ravel().tolist(), [.35, .22, .1],
            v=[-5., .4, -.1], w=[.1, .2, 1.])
    i = n.inflow([0., 0., 1.], np.pi, .1)
    n.pose_inflow(i, .3, impose=True)
    kw = dict(centre=[0., 0., 0.], e1=[1., 0., 0.], vtip=10.)
    if kind == 'harmoniques':
        n.pose_inflow_harmoniques(i, .3, 1.2, -.8, **kw)
    elif kind == 'profil':
        n.pose_inflow_profil(i, [0., 1.], [.1, 1.4], **kw)
    else:
        n.pose_inflow_carte(i, [0., 1.], [0., np.pi/2, np.pi, 3*np.pi/2],
                           [0., .2, .3, .1, 1.5, 1., .2, .3], **kw)
    a = [-90., -45., 0., 45., 90.]
    n.pale('p', 0, [0., 0., 0.], [1., 0., 0.], [0., 1., 0.], .25, .1,
           (a, [-1., -.7, 0., .7, 1.], [.2, .05, .01, .05, .2], [0.]*5), inflow=i)
    return n


def _roulement_charge():
    n = Noyau([0., 0., -9.81])
    n.corps('roue', 2., np.diag([.2, .3, .4]).ravel().tolist(), [0., 0., 1.],
            v=[1., 0., 0.], w=[0., 1., 0.])
    n.liaison('roulement', None, 0, pa=[0., 0., 0.], bloque_t=[0, 2], bloque_r=[], nh=True)
    return n


def _couple_precharge():
    n = Noyau([0., 0., 0.])
    n.corps('arbre', 1., J, [0., 0., 0.])
    n.couple('ressort', None, 0, [0., 0., 1.], ('ressort', [1e12, 1., 1.]))
    return n


def _tangentes_champ_fd(n):
    """Contrôle global indépendant des tangentes locales, sans contraintes."""
    import copy
    state = n.etat_precis()
    dim = 6 * len(state[1])
    results = []
    for velocity, eps in ((False, 1e-5), (True, 1e-6)):
        matrix = np.zeros((dim, dim))
        for col in range(dim):
            body, component = divmod(col, 6)
            values = []
            for sign in (1., -1.):
                s = copy.deepcopy(state)
                if velocity:
                    s[3 + component//3][body][component % 3] += sign*eps
                elif component < 3:
                    s[1][body][component] += sign*eps
                else:
                    d = np.zeros(3); d[component-3] = sign*eps
                    skew = np.array([[0., -d[2], d[1]], [d[2], 0., -d[0]],
                                     [-d[1], d[0], 0.]])
                    s[2][body] = (expm(skew) @ np.asarray(s[2][body]).reshape(3, 3)).ravel().tolist()
                n.pose_etat(*s)
                values.append(np.asarray(n.residu_statique()))
            matrix[:, col] = (values[0]-values[1])/(2*eps)
        results.append(matrix)
    n.pose_etat(*state)
    return results


def _cascade_initiale(nb, q):
    alpha = math.pi/3
    length, width, mass, k = 1./nb, .25, 1./(3*nb), 100./nb
    arm = length*np.array([math.cos(alpha), math.sin(alpha), 0.])/2
    n, base = Noyau((q @ [0., -9.81, 0.]).tolist()), np.zeros(3)
    for i in range(nb):
        parent = 3*i-2 if i else None
        for j, p in enumerate([base+arm, base+2*arm+[width/2, 0., 0.],
                                base+arm+[width, 0., 0.]]):
            n.corps(f'b{i}_{j}', mass, (mass*np.eye(3)/12).ravel().tolist(),
                    (q @ p).tolist(), rot=q.ravel().tolist())
        specs = [(parent, 3*i, [-width/2, 0., 0.] if i else [0., 0., 0.]),
                 (3*i, 3*i+1, arm), (3*i+1, 3*i+2, [width/2, 0., 0.]),
                 (parent, 3*i+2, [width/2, 0., 0.] if i else [width, 0., 0.])]
        for j, (a, b, pa) in enumerate(specs):
            n.liaison(f'p{i}_{j}', a, b, pa=list(q @ pa if a is None else pa),
                      ra=(q if a is None else np.eye(3)).ravel().tolist(),
                      bloque_t=[0, 1, 2], bloque_r=[0, 1])
        axe = q @ [0., 0., 1.] if parent is None else [0., 0., 1.]
        n.couple(f'k{i}', parent, 3*i, list(axe), ('ressort', [k, 0., 0.]))
        base += 2*arm
    return n


class AuditNoyau(unittest.TestCase):
    def test_acceleration_initiale_conservee_en_chute_libre(self):
        g = np.array([1., -2., -9.81])
        r0, v0 = np.array([2., 3., 4.]), np.array([.2, -.3, .4])
        for h in (1., .1, .01):
            for par_effort in (False, True):
                n = Noyau([0., 0., 0.] if par_effort else g.tolist())
                n.corps('libre', 2., J, r0.tolist(), v=v0.tolist())
                if par_effort:
                    n.effort(0, (2.*g).tolist(), [0., 0., 0.])
                n.simule(1., h)
                np.testing.assert_allclose(n.etat()[1][0], r0+v0+.5*g, atol=1e-11, rtol=0.)
                np.testing.assert_allclose(n.etat()[3][0], v0+g, atol=1e-11, rtol=0.)

    def test_acceleration_angulaire_initiale_conservee(self):
        for h in (1., .1, .01):
            n = Noyau([0., 0., 0.])
            n.corps('rotor', 1., J, [0., 0., 0.])
            n.effort(0, [0., 0., 0.], [0., 0., 3.])
            n.simule(1., h)
            np.testing.assert_allclose(n.etat()[4][0], [0., 0., 3.], atol=1e-11, rtol=0.)
            expected = expm(np.array([[0., -1.5, 0.], [1.5, 0., 0.], [0., 0., 0.]]))
            np.testing.assert_allclose(np.array(n.etat()[2][0]).reshape(3, 3), expected,
                                       atol=1e-11, rtol=0.)

    def test_acceleration_initiale_bielle_sans_pas_de_perturbation(self):
        for length, speed in [(1., 1.), (1e-6, 10.), (1., 1e8), (1e6, 1e8)]:
            n = Noyau([0., 0., 0.])
            n.corps('masse', 1., J, [length, 0., 0.], v=[0., speed, 0.])
            n.distance('bielle', None, 0, [0., 0., 0.], [0., 0., 0.])
            n.enregistre_schema()
            before = n.etat_precis()
            n.simule(0., .001)
            self.assertEqual(n.etat_precis(), before)
            a = n.schema()[0][0]
            np.testing.assert_allclose(a, [-speed*speed/length, 0., 0., 0., 0., 0.], rtol=2e-14, atol=0.)

    def test_initialisation_pendule_vertical_urdf(self):
        from pathlib import Path
        from vinkulum import urdf
        n, _ = urdf.charge(str(Path(urdf.__file__).parent/'donnees/urdfs/cartpole.urdf'))
        n.enregistre_schema()
        n.simule(0., .001)
        np.testing.assert_allclose(n.schema()[0][0], 0., atol=2e-14)

    def test_acceleration_commandes_sur_deux_axes(self):
        # R_b = R_x(2t) R_z(3t) : alpha_b(0) = 2 e_x × 3 e_z = -6 e_y.
        # Les poses et les vitesses satisfont les commandes dès t = 0.
        n = Noyau([0., 0., 0.])
        n.corps('porteur', 1., J, [0., 0., 0.], w=[2., 0., 0.])
        n.corps('arbre', 1., J, [0., 0., 0.], w=[2., 0., 3.])
        n.liaison('x', None, 0, cible_r=([1., 0., 0.], ('lineaire', [0., 2.])))
        n.liaison('z', 0, 1, cible_r=([0., 0., 1.], ('lineaire', [0., 3.])))
        np.testing.assert_allclose(n.phi_dot(), 0., atol=1e-12)
        n.enregistre_schema()
        n.simule(0., .001)
        a = np.array(n.schema()[0][0]).reshape(2, 6)
        np.testing.assert_allclose(a, [[0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., -6., 0.]], atol=1e-13)

    def test_contact_nonlisse_decollement_sans_attraction(self):
        n = Noyau([0., 0., 0.])
        n.corps('a', 1., J, [0., 0., 0.])
        n.corps('b', 1., J, [1., 0., 0.], v=[0., 1., 0.])
        n.contact('tangence', 0, [0., 0., 0.], .5, b=1, pb=[0., 0., 0.], rayon_b=.5, nonlisse=True)
        n.enregistre_schema()
        n.simule(.0001, .0001)
        np.testing.assert_allclose(n.schema()[-1][0], 0., atol=1e-14)
        n.simule(.001, .0001)
        np.testing.assert_allclose(n.etat()[3], [[0., 0., 0.], [0., 1., 0.]], atol=1e-12)

        # Rebond : dès que la bille s'écarte du plan, elle retrouve -g.
        n = Noyau([0., 0., -9.81])
        n.corps('bille', 1., J, [0., 0., .5], v=[0., 0., -1.])
        n.contact('sol', 0, [0., 0., 0.], .5, nonlisse=True, restitution=.8)
        n.enregistre_schema()
        n.simule(.0001, .0001)
        self.assertGreater(n.etat()[3][0][2], 0.)
        np.testing.assert_allclose(n.schema()[-1][0], [0., 0., -9.81, 0., 0., 0.], atol=1e-12)

    def test_analyse_et_initialisation_cascade_redondante(self):
        from scipy.spatial.transform import Rotation
        nb, alpha = 8, math.pi/3
        length, mass, k = 1./nb, 1./(3*nb), 100./nb
        for q in (np.eye(3), Rotation.from_rotvec([.31, -.47, .22]).as_matrix()):
            n = _cascade_initiale(nb, q)
            before = n.etat_precis()
            K, C, M, Z, G = map(np.asarray, n.k_c_m_z())
            self.assertEqual(n.etat_precis(), before)
            # Paramétrage indépendant par les nb angles absolus des bras.
            S = np.zeros((18*nb, nb))
            d = q @ (length*np.array([-math.sin(alpha), math.cos(alpha), 0.])/2)
            for i in range(nb):
                for j in range(i+1):
                    for b in range(3):
                        S[18*i+6*b:18*i+6*b+3, j] = d*(1. if j == i and b != 1 else 2.)
                    if j == i:
                        S[18*i+3:18*i+6, j] = q[:, 2]
                        S[18*i+15:18*i+18, j] = q[:, 2]
            masses_ref = np.tile([mass]*3+[mass/12]*3, 3*nb)
            np.testing.assert_allclose(M, np.diag(masses_ref), atol=1e-15)
            np.testing.assert_allclose(G @ S, 0., atol=2e-14)
            self.assertEqual(Z.shape, (18*nb, nb))
            np.testing.assert_allclose(G @ Z, 0., atol=2e-14)
            np.testing.assert_allclose(Z.T @ Z, np.eye(nb), atol=2e-14)
            # Ici les accélérations sont tangentes et la courbure des poses
            # radiale : SᵀKS égale la Hessienne du potentiel réduit.
            kr = k-mass*9.81*length*(3*(nb-np.arange(nb))-1)*math.sin(alpha)
            np.testing.assert_allclose(S.T @ K @ S, np.diag(kr), atol=2e-11)
            np.testing.assert_allclose(C, 0., atol=1e-15)
            force = np.zeros(18*nb)
            for i in range(3*nb):
                force[6*i:6*i+3] = mass*(q @ [0., -9.81, 0.])
            acceleration = S @ np.linalg.solve(S.T @ (masses_ref[:, None]*S), S.T @ force)
            n.enregistre_schema()
            n.simule(0., .001)
            np.testing.assert_allclose(n.schema()[0][0], acceleration, rtol=0., atol=2e-11)
            self.assertEqual(n.etat_precis(), before)

    def test_statique_cascade_de_parallelogrammes(self):
        self._statique_cascade(np.eye(3))

    def test_statique_cascade_tournee_et_restauration(self):
        for v in ([.31, -.47, .22], [-1.1, .6, 2.2]):
            x, y, z = v
            q = expm(np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]]))
            with self.subTest(rotation=v):
                self._statique_cascade(q)

    def _statique_cascade(self, q):
        # Graphe connexe : 8 mobilités et 24 dépendances qui ne sont pas
        # des copies de lignes. Référence obtenue par l'énergie mécanique.
        nb, alpha = 8, math.pi/3
        length, width, mass, k = 1./nb, .25, 1./(3*nb), 100./nb
        arm = length*np.array([math.cos(alpha), math.sin(alpha), 0.])/2
        n = _cascade_initiale(nb, q)
        initial = n.etat_precis()
        with self.assertRaises(RuntimeError):
            n.statique(strict=True, tol=1e-8, iters=1, paliers_max=1)
        self.assertEqual(n.etat_precis(), initial)
        n.statique(strict=True, tol=1e-8, iters=30, paliers_max=1)
        self.assertEqual(n.statique_info()['statut'], 'tolerance')
        state = n.etat_precis()
        pos = np.asarray(state[1], dtype=np.longdouble)+np.asarray(state[6], dtype=np.longdouble)
        pos = pos @ q.astype(np.longdouble)
        rot = np.einsum('ab,nbc->nac', q.T, np.asarray(state[2]).reshape(3*nb, 3, 3))
        base = np.zeros(3, dtype=np.longdouble)
        for i in range(nb):
            # V_i = mgL(3(N-i)-1) sin(alpha+delta) + k delta²/2.
            # V_i'' > 0 sur [-.4, 0] : racine unique, encadrée par bissection.
            lo, hi = -.4, 0.
            for _ in range(60):
                delta = (lo+hi)/2
                if k*delta+mass*9.81*length*(3*(nb-i)-1)*math.cos(alpha+delta) < 0:
                    lo = delta
                else:
                    hi = delta
            delta = (lo+hi)/2
            arm = length*np.array([math.cos(alpha+delta), math.sin(alpha+delta), 0.])/2
            expected = [base+arm, base+2*arm+[width/2, 0., 0.], base+arm+[width, 0., 0.]]
            c, s = math.cos(delta), math.sin(delta)
            r = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
            np.testing.assert_allclose(pos[3*i:3*i+3], expected, rtol=0., atol=1e-9)
            np.testing.assert_allclose(rot[3*i:3*i+3], [r, np.eye(3), r], rtol=0., atol=1e-9)
            base += 2*arm

    def test_statique_mobile_liaisons_equivalentes(self):
        # Même chaîne mobile, puis multiplicité 1/2/3 avec repères opposés.
        # Les réactions individuelles doivent se partager la charge ; les
        # ressorts doivent équilibrer le moment des poids de toute la suite.
        nb = 18
        solutions = []
        for repetitions in ([1]*nb, [1+i % 3 for i in range(nb)]):
            n = Noyau([0., -9.81, 0.])
            for i, count in enumerate(repetitions):
                n.corps(str(i), 1./nb, (np.eye(3)/(12*nb**3)).ravel().tolist(),
                        [(i+.5)/nb, 0., 0.])
                for j in range(count):
                    s = -1. if j % 2 else 1.
                    n.liaison(f'j{i}_{j}', i-1 if i else None, i,
                              pa=[.5/nb, 0., 0.] if i else [0., 0., 0.],
                              ra=np.diag([s, s, 1.]).ravel().tolist(),
                              bloque_t=[0, 1, 2], bloque_r=[0, 1])
                n.couple(f'k{i}', i-1 if i else None, i, [0., 0., 1.],
                         ('ressort', [10.*nb, 0., 0.]))
            n.statique(strict=True, tol=1e-8, iters=30, paliers_max=1)
            self.assertEqual(n.statique_info()['statut'], 'tolerance')
            self.assertGreater(n.statique_info()['evaluations'], 1)
            state = n.etat_precis()
            pos = np.asarray(state[1], dtype=np.longdouble) + np.asarray(state[6], dtype=np.longdouble)
            rot = np.asarray(state[2]).reshape(nb, 3, 3)
            angles = np.arctan2(rot[:, 1, 0], rot[:, 0, 0])
            self.assertLess(pos[-1, 1], -.1)
            points = np.zeros((nb, 3), dtype=np.longdouble)
            points[1:] = pos[:-1] + rot[:-1, :, 0]/(2*nb)
            np.testing.assert_allclose(pos-rot[:, :, 0]/(2*nb), points, rtol=0., atol=1e-12)
            torque = np.array([9.81/nb*np.sum(pos[i:, 0]-points[i, 0]) for i in range(nb)])
            np.testing.assert_allclose(-10.*nb*(angles-np.r_[0., angles[:-1]]), torque,
                                       rtol=0., atol=1e-9)
            reactions = dict(n.reactions())
            for i, count in enumerate(repetitions):
                r = np.eye(3) if i == 0 else rot[i-1]
                force = r.T @ [0., -9.81*(nb-i)/nb, 0.]
                for j in range(count):
                    s = -1. if j % 2 else 1.
                    expected = np.r_[np.array([s, s, 1.])*force/count, 0., 0.]
                    np.testing.assert_allclose(reactions[f'j{i}_{j}'], expected, rtol=0., atol=1e-9)
            solutions.append((pos, rot))
        for original, redundant in zip(*solutions):
            np.testing.assert_allclose(original, redundant, rtol=0., atol=1e-11)

    def test_statique_equivalentes_conserve_contradiction_et_restaure(self):
        n = Noyau([0., 0., 0.])
        n.corps('corps', 1., J, [0., 0., 0.])
        for i, target in enumerate([0., 0., .1]):
            n.liaison(str(i), None, 0, bloque_t=[0], bloque_r=[],
                      cible_t=([1., 0., 0.], ('constante', [target])))
        n.liaison('reste', None, 0, bloque_t=[1, 2], bloque_r=[0, 1, 2])
        initial = n.etat_precis()
        with self.assertRaises(RuntimeError):
            n.statique(strict=True, tol=1e-8, iters=4, paliers_max=1)
        self.assertEqual(n.etat_precis(), initial)
        self.assertEqual(n.statique_info()['statut'], 'echec')
        self.assertGreater(n.statique_info()['contraintes'], .049)

    def test_statique_princeton_encastrement_double(self):
        # Le système singulier désactivait le mérite naturel : dupliquer
        # l'encastrement faisait échouer cette même poutre à charge pleine.
        ne = 10
        force = np.array([0., 8.896/math.sqrt(2), 8.896/math.sqrt(2)])
        for formulation in ('milieu', 'integree'):
            states = []
            for copies in (1, 2):
                n = Noyau([0., 0., 0.])
                for i in range(ne+1):
                    n.corps(str(i), 1., J, [.508*i/ne, 0., 0.])
                for j in range(copies):
                    n.liaison(f'enc{j}', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
                for i in range(ne):
                    n.poutre(str(i), i, i+1, 2.84191e6, 6.40131e5, 3.10338, 36.2794,
                             ei3=2.42873, ga3=9.03881e5, formulation=formulation)
                n.effort(ne, force.tolist(), [0., 0., 0.])
                n.statique(strict=True, tol=1e-8, iters=100, paliers_max=1)
                self.assertEqual(n.statique_info()['statut'], 'tolerance')
                state = n.etat_precis()
                pos = np.asarray(state[1], dtype=np.longdouble)+np.asarray(state[6], dtype=np.longdouble)
                reaction = np.r_[force, np.cross(pos[-1], force)]/copies
                np.testing.assert_allclose([v for _, v in n.reactions()],
                                           np.tile(reaction, (copies, 1)), rtol=0., atol=1e-8)
                states.append((pos, np.asarray(state[2])))
            for original, doubled in zip(*states):
                np.testing.assert_allclose(original, doubled, rtol=0., atol=1e-10)

    def test_statique_boucle_reactions_norme_minimale(self):
        # Chaîne soudée fermée au sol à ses deux extrémités. La pesanteur
        # ne sollicite que Fz et My ; leurs deux auto-contraintes donnent
        # une référence de norme minimale par un calcul de taille 2.
        nb, length, weight = 24, .25, 9.81
        n = Noyau([0., 0., -weight])
        for i in range(nb):
            n.corps(str(i), 1., J, [length*i, 0., 0.])
            n.liaison('j'+str(i), i-1 if i else None, i,
                      pa=[length, 0., 0.] if i else [0., 0., 0.],
                      bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
        n.liaison('fermeture', None, nb-1, pa=[length*(nb-1), 0., 0.],
                  bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
        initial = n.etat_precis()
        n.statique(strict=True, tol=1e-8, iters=2, paliers_max=1)
        self.assertEqual(n.etat_precis(), initial)
        self.assertEqual(n.statique_info()['statut'], 'tolerance')
        tail = np.arange(nb, 0, -1)
        f = -weight*tail
        moment = weight*length*tail*(tail-1)/2
        arm = length*(tail-1)
        coefficients = [[nb+1+np.dot(arm, arm), -np.sum(arm)], [-np.sum(arm), nb+1]]
        c, d = np.linalg.solve(coefficients, [-np.sum(-f+arm*moment), np.sum(moment)])
        reference = np.zeros((nb+1, 6))
        reference[:-1, 2], reference[:-1, 4] = f-c, moment+arm*c-d
        reference[-1, 2], reference[-1, 4] = c, d
        np.testing.assert_allclose([v for _, v in n.reactions()], reference, rtol=2e-13, atol=2e-11)

    def test_statique_directions_presque_dependantes(self):
        # Deux guides bloquent x et une direction presque parallèle.
        # Le corps est immobile : ses réactions valent 2-cot(a), 1/sin(a).
        # Former GGᵀ perdait la seconde direction dès delta=1e-6 sur ce témoin.
        import math
        for delta in (1e-3, 1e-6, 1e-9, 1e-12, 0.):
            n = Noyau([0., 0., 0.])
            n.corps('support', 1., J, [0., 0., 0.])
            angle = math.atan(delta)
            for name, a in (('x', 0.), ('proche', angle)):
                c, s = math.cos(a), math.sin(a)
                n.liaison(name, None, 0, bloque_t=[0], bloque_r=[],
                          ra=[c, -s, 0., s, c, 0., 0., 0., 1.])
            n.liaison('reste', None, 0, bloque_t=[2], bloque_r=[0, 1, 2])
            n.effort(0, [2., 1., 3.], [0., 0., 0.])
            initial = n.etat_precis()
            if delta == 0.:
                # À coïncidence exacte, y reste libre et aucun équilibre n'existe.
                with self.assertRaises(RuntimeError):
                    n.statique(strict=True, tol=1e-8, iters=2, paliers_max=1)
                self.assertEqual(n.statique_info()['statut'], 'echec')
            else:
                residual, iterations = n.statique(strict=True, tol=1e-8, iters=10, paliers_max=1)
                self.assertEqual(iterations, 0)
                self.assertLess(residual, 1e-8)
                reactions = dict(n.reactions())
                np.testing.assert_allclose([reactions['x'][0], reactions['proche'][0]],
                    [2.-1./math.tan(angle), 1./math.sin(angle)], rtol=2e-15, atol=1e-12)
            self.assertEqual(n.etat_precis(), initial)

    def test_tangente_roulement_charge(self):
        n = _roulement_charge()
        before = n.etat_precis()
        for h in (.01, .05):
            blocs, _ = n.audit_jacobien(h)
            self.assertLess(max(e for _, e, _ in blocs), 1e-9, (h, blocs))
            self.assertEqual(n.etat_precis(), before)

    def test_amortissement_sous_precharge(self):
        # τ=10^12(1-θ)-ω : la pente en vitesse vaut -1, même quand
        # une différence de deux grands couples f64 perd cette variation.
        n = _couple_precharge()
        self.assertAlmostEqual(n.k_c_m_z()[1][5][5], 1., places=12)

    def test_tangentes_inflow_spatial(self):
        # Trois distributions, à l'intérieur de leurs cellules d'interpolation.
        # Le contrôle porte aussi sur la dynamique : le précédent jacobien
        # omettait les dérivées de translation, écart 4 à 6e-6 à h=.05.
        for kind in ('harmoniques', 'profil', 'carte'):
            n = _pale_inflow_spatial(kind)
            k, c, *_ = n.k_c_m_z()
            kfd, cfd = _tangentes_champ_fd(n)
            self.assertGreater(np.linalg.norm(kfd[:, :3]), .001)
            np.testing.assert_allclose(k, kfd, rtol=2e-7, atol=2e-8)
            np.testing.assert_allclose(c, cfd, rtol=2e-7, atol=2e-8)
            before = n.etat_precis()
            for h in (.01, .05):
                blocs, _ = n.audit_jacobien(h)
                self.assertLess(max(e for _, e, _ in blocs), 1e-9, (kind, h, blocs))
                self.assertEqual(n.etat_precis(), before)

    def test_tangentes_locales_super_et_contacts(self):
        # Superélément à trois nœuds non contigus, contacts croisés et bâti,
        # amortissement, vitesses de rotation et géométrie déformée.
        n = Noyau([0., 0., -9.81])
        for i, p in enumerate(([.7, .06, .46], [0., 0., .46], [.35, .04, .46])):
            n.corps(str(i), 1., np.diag([.2, .3, .4]).ravel().tolist(), p,
                    v=[.2, -.3, .1], w=[.1, .2, .3])
        k0 = np.diag(np.arange(1., 19.))
        n.superelement('s', [1, 2, 0], k0.ravel().tolist(), beta=.07)
        n.contact('paire', 0, [.03, .01, 0.], .3, b=2, pb=[0., .02, 0.],
                  rayon_b=.2, k=170., c=.3, mu=.2, v_eps=.05)
        n.contact('sol', 1, [.02, .03, 0.], .5, k=110., c=.2, mu=.3, v_eps=.05)
        s = n.etat_precis()
        s[1][2][1] += .03
        s[2][0] = expm(np.array([[0., -.05, .03], [.05, 0., -.02],
                                [-.03, .02, 0.]])).ravel().tolist()
        n.pose_etat(*s)
        k, c, *_ = n.k_c_m_z()
        kfd, cfd = _tangentes_champ_fd(n)
        np.testing.assert_allclose(k, kfd, rtol=3e-6, atol=2e-6)
        np.testing.assert_allclose(c, cfd, rtol=3e-6, atol=2e-6)

    def test_contraintes_consument_partie_basse(self):
        for element in ('liaison', 'distance', 'vis'):
            n = Noyau([0., 0., 0.])
            n.corps('mobile', 1., J, [1., 0., 0.])
            if element == 'liaison':
                n.liaison('fixe', None, 0, pa=[1., 0., 0.], bloque_t=[0], bloque_r=[])
            elif element == 'distance':
                n.distance('bielle', None, 0, [0., 0., 0.], [0., 0., 0.], l=1.)
            else:
                n.vis('vis', None, 0, [1., 0., 0.], 1.)
            state = n.etat_precis()
            state[6][0][0] = 2.**-60
            n.pose_etat(*state)
            self.assertEqual(n.phi()[0], 2.**-60)
            n.assemble(tol=1e-22)
            self.assertLess(abs(n.phi()[0]), 1e-22)

    def test_contact_translations_compensees(self):
        for mobile in (False, True):
            results = []
            for origin in (0., 2.**40):
                n = Noyau([0., 0., 0.])
                n.corps('sphere', 1., J, [origin + 1., 0., 0.])
                options = dict(origine=[origin, 0., 0.])
                if mobile:
                    n.corps('support', 1., J, [origin, 0., 0.])
                    options = dict(b=1)
                n.contact('plan', 0, [0., 0., 0.], 1., normale=[1., 0., 0.],
                          k=1e6, expo=1.5, **options)
                state = n.etat_precis()
                state[6][0][0] = -2.**-20
                n.pose_etat(*state)
                results.append(np.asarray(n.residu_statique()))
            self.assertGreater(abs(results[0][0]), 0.)
            np.testing.assert_array_equal(results[0], results[1])

    def test_positions_compensees_statique_et_restauration(self):
        # Allongement physique 2^-60, invisible dans la position f64 1+u.
        # Poutre et superélément doivent garder la même solution analytique.
        for element in ('poutre', 'super'):
            n = Noyau([0., 0., 0.])
            for i in range(2):
                n.corps(str(i), 1., J, [float(i), 0., 0.])
            n.liaison('enc', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
            if element == 'poutre':
                n.poutre('barre', 0, 1, 2.**30, 1., 1., 1.)
            else:
                k = np.eye(12)
                k[6, 6] = 2.**30
                n.superelement('barre', [0, 1], k.ravel().tolist())
            force = 2.**-30
            n.effort(1, [force, 0., 0.], [0., 0., 0.])
            n.statique(tol=1e-14, strict=True)
            state = n.etat_precis()
            residual = np.array(n.residu_statique())
            self.assertLess(abs(residual[6]), 1e-20)
            self.assertEqual(state[1][1][0], 1.)
            self.assertAlmostEqual(state[6][1][0] / 2.**-60, 1., places=12)
            n.pose_etat(*n.etat())
            self.assertEqual(n.etat_precis()[6][1][0], 0.)
            self.assertAlmostEqual(abs(n.residu_statique()[6]) / force, 1., places=12)
            n.pose_etat(*state)
            self.assertEqual(n.etat_precis(), state)
            np.testing.assert_array_equal(n.residu_statique(), residual)
            # Une tentative modifie effectivement les petits termes, puis
            # le refus pour budget insuffisant doit tous les restaurer.
            n.effort(1, [force, 0., 0.], [0., 0., 0.])
            with self.assertRaises(RuntimeError):
                n.statique(tol=1e-14, iters=1, paliers_max=1, strict=True)
            self.assertEqual(n.etat_precis(), state)

    def test_etat_precis_validation_atomique(self):
        n = masses()
        state = n.etat_precis()
        for problem in ('longueur', 'nan', 'rotation', 'debordement'):
            import copy
            bad = copy.deepcopy(list(state))
            if problem == 'longueur':
                bad[6].pop()
            elif problem == 'nan':
                bad[6][-1][0] = float('nan')
            elif problem == 'rotation':
                bad[2][-1][0] = 2.
            else:
                bad[1][-1][0] = bad[6][-1][0] = 1.7e308
            with self.assertRaises(ValueError):
                n.pose_etat(*bad)
            self.assertEqual(n.etat_precis(), state)

    def test_positions_compensees_integrateurs(self):
        for scheme in ('alpha', 'energie_moment', 'multirythme'):
            n = Noyau([0., 0., 0.])
            n.corps('libre', 1., J, [1., 0., 0.], v=[2.**-60, 0., 0.], w=[.1, .2, .3])
            if scheme == 'alpha':
                n.simule(1., .125)
            elif scheme == 'energie_moment':
                n.simule_em(1., .125)
            else:
                n.simule_multirythme(1., .125, rapides=[0], k=2)
            state = n.etat_precis()
            self.assertEqual(state[1][0][0], 1.)
            self.assertAlmostEqual(state[6][0][0] / 2.**-60, 1., places=12)
            n.pose_etat(*n.etat())
            n.pose_etat(*state)
            self.assertEqual(n.etat_precis(), state)

    def test_raideur_locale_contre_champ_global(self):
        # Deux poutres partagent un nœud ; les indices ne suivent pas leur
        # ordre géométrique. Déformation, couple et inertie tournante doivent
        # tous apparaître dans la tangente, sans double comptage.
        n = Noyau([0., 0., -9.81])
        for i, p in enumerate(([.8, 0., 0.], [0., 0., 0.], [.4, 0., 0.])):
            n.corps(str(i), 1., [1.,0.,0.,0.,2.,0.,0.,0.,2.], p,
                    w=[.1, .2, .3])
        n.poutre('gauche', 1, 2, 2000., 700., 30., 50., ei3=20.)
        n.poutre('droite', 2, 0, 3000., 800., 40., 70., ei3=25.)
        n.couple('rappel', None, 0, [0.,0.,1.], ('ressort', [5.,.2,.1]))
        t, p, rot, v, w, vi = n.etat()
        p[2][1] += .01
        rot[2] = expm(np.array([[0., -.08, .03], [.08, 0., -.02],
                              [-.03, .02, 0.]])).ravel().tolist()
        n.pose_etat(t, p, rot, v, w, vi)
        k = np.asarray(n.k_c_m_z()[0])
        reference = np.zeros_like(k)
        eps = 1e-5
        for col in range(18):
            body, component = divmod(col, 6)
            values = []
            for sign in (1., -1.):
                pp, rr = np.array(p), np.array(rot).reshape(3,3,3).copy()
                if component < 3:
                    pp[body, component] += sign*eps
                else:
                    d = np.zeros(3); d[component-3] = sign*eps
                    skew = np.array([[0.,-d[2],d[1]], [d[2],0.,-d[0]],
                                     [-d[1],d[0],0.]])
                    rr[body] = expm(skew) @ rr[body]
                n.pose_etat(t, pp.tolist(), rr.reshape(3,9).tolist(), v, w, vi)
                values.append(np.array(n.residu_statique()))
            reference[:,col] = (values[0]-values[1])/(2*eps)
        n.pose_etat(t, p, rot, v, w, vi)
        np.testing.assert_allclose(k, reference, rtol=1e-7, atol=2e-7)
        damping = np.asarray(n.k_c_m_z()[1])
        reference = np.zeros_like(damping)
        eps = 1e-6
        for col in range(18):
            body, component = divmod(col, 6)
            values = []
            for sign in (1., -1.):
                vv, ww = np.array(v), np.array(w)
                if component < 3:
                    vv[body, component] += sign*eps
                else:
                    ww[body, component-3] += sign*eps
                n.pose_etat(t, p, rot, vv.tolist(), ww.tolist(), vi)
                values.append(np.array(n.residu_statique()))
            reference[:,col] = (values[0]-values[1])/(2*eps)
        n.pose_etat(t, p, rot, v, w, vi)
        np.testing.assert_allclose(damping, reference, rtol=1e-6, atol=2e-7)

    def test_poutre_raideur_timoshenko(self):
        # Matrice nodale analytique, dans les deux plans anisotropes.
        n = Noyau([0., 0., 0.])
        length, ea, gay, gaz, gj, eiy, eiz = .9, 10000., 700., 900., 30., 50., 20.
        for i in range(2):
            n.corps(str(i), 1., J, [i*length, 0., 0.])
        n.poutre('p', 0, 1, ea, gay, gj, eiy, ei3=eiz, ga3=gaz, formulation="integree")
        expected = np.zeros((12, 12))
        for indices, stiff in (([0, 6], ea), ([3, 9], gj)):
            expected[np.ix_(indices, indices)] = stiff/length*np.array([[1., -1.], [-1., 1.]])
        for indices, ei, ga, signs in (([1, 5, 7, 11], eiz, gay, [1, 1, 1, 1]),
                                       ([2, 4, 8, 10], eiy, gaz, [1, -1, 1, -1])):
            l, phi = length, 12*ei/(ga*length**2)
            block = ei/(l**3*(1+phi))*np.array([
                [12, 6*l, -12, 6*l],
                [6*l, (4+phi)*l*l, -6*l, (2-phi)*l*l],
                [-12, -6*l, 12, -6*l],
                [6*l, (2-phi)*l*l, -6*l, (4+phi)*l*l]])
            expected[np.ix_(indices, indices)] = block*np.outer(signs, signs)
        np.testing.assert_allclose(n.k_c_m_z()[0], expected, rtol=1e-6, atol=1e-6)

    def test_poutre_arc_moment_pur(self):
        # Un moment constant donne un arc circulaire. La même solution doit
        # tenir après rotation et translation du modèle dans le monde.
        skew = np.array([[0., -.3, .2], [.3, 0., -.1], [-.2, .1, 0.]])
        origin = np.array([-.4, .2, .1])
        for frame in (np.eye(3), expm(skew)):
            for theta in (.01, 1., 2.5):
                n = Noyau([0., 0., 0.])
                for i in range(2):
                    n.corps(str(i), 1., J, (origin+frame@np.array([float(i), 0., 0.])).tolist(),
                            rot=frame.ravel().tolist())
                n.liaison('enc', None, 0, pa=origin.tolist(), bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
                n.poutre('p', 0, 1, 1e6, 1e5, 50., 100., formulation='integree')
                n.effort(1, [0., 0., 0.], (frame@np.array([0., 0., 100.*theta])).tolist())
                n.statique(tol=1e-9, iters=100)
                exact = origin+frame@np.array([math.sin(theta)/theta, (1-math.cos(theta))/theta, 0.])
                np.testing.assert_allclose(n.pose(1)[0], exact, atol=2e-8, rtol=0)
                self.assertLess(np.max(np.abs(n.residu_statique()[6:])), 1e-5)
                self.assertLess(np.max(np.abs(n.phi())), 1e-9)

    def test_poutre_sensibilites_condensation(self):
        # La condensation dépend de GA ET EI : dériver seulement les termes
        # de courbure donne une fausse sensibilité, même si l'équilibre passe.
        base = np.array([2000., 700., 30., 50., 20., 900.])
        def modele(params):
            ea, ga, gj, ei, ei3, ga3 = params
            n = Noyau([0., 0., 0.])
            for i in range(2):
                n.corps(str(i), 1., J, [.9*i, 0., 0.])
            n.poutre('p', 0, 1, ea, ga, gj, ei, ei3=ei3, ga3=ga3, formulation="integree")
            state = list(n.etat())
            state[1][1] = [.91, .03, -.02]
            state[2][1] = expm(np.array([[0., -.2, .1], [.2, 0., -.3], [-.1, .3, 0.]])).ravel().tolist()
            n.pose_etat(*state)
            return n
        n = modele(base)
        for key, indices in {'ea':[0], 'ga':[1, 5], 'gj':[2], 'ei':[3, 4],
                             'ei_e2':[3], 'ei_e3':[4]}.items():
            delta = min(base[indices])*1e-4
            pp, pm = base.copy(), base.copy()
            pp[indices] += delta
            pm[indices] -= delta
            plus, minus = modele(pp), modele(pm)
            force = (np.array(plus.residu_statique())-minus.residu_statique())/(2*delta)
            stiffness = (np.array(plus.k_c_m_z()[0])-minus.k_c_m_z()[0])/(2*delta)
            np.testing.assert_allclose(n.d_residu_poutre(0, key), force, rtol=1e-6, atol=1e-8)
            np.testing.assert_allclose(n.d_raideur_poutre(0, key), stiffness, rtol=2e-5, atol=3e-5)

    def test_statique_poutre_anisotrope_amortie(self):
        # Charge hors des axes principaux : l'ancien amortissement stagne,
        # même par 512 paliers. Référence MBDyn raffinée à 160 intervalles.
        ref = np.array([0.4936103410461536, 0.1090430831170913,
                        0.00888888595701552])
        for formulation in ("milieu", "integree"):
            erreurs = []
            for ne in (10, 20, 40):
                n = Noyau([0., 0., 0.])
                ids = [n.corps(str(i), 1., J, [.508*i/ne, 0., 0.])
                       for i in range(ne+1)]
                n.liaison('enc', None, ids[0], bloque_t=[0, 1, 2],
                          bloque_r=[0, 1, 2])
                for i in range(ne):
                    n.poutre(str(i), ids[i], ids[i+1], 2.84191e6, 6.40131e5,
                             3.10338, 36.2794, ei3=2.42873, ga3=9.03881e5, formulation=formulation)
                n.effort(ids[-1], [0., 8.896/math.sqrt(2), 8.896/math.sqrt(2)],
                         [0., 0., 0.])
                _, iterations = n.statique(tol=1e-8, iters=60, paliers_max=1)
                self.assertLess(iterations, 60)
                self.assertLess(max(map(abs, n.phi())), 1e-10)
                self.assertLess(np.max(np.abs(n.residu_statique()[6:])), 1e-7)
                erreurs.append(np.max(np.abs(np.array(n.pose(ids[-1])[0])-ref)))
            self.assertLess(erreurs[1], 7e-5)
            self.assertLess(erreurs[1], erreurs[0]/3)
            self.assertLess(erreurs[2], 2e-5)
            self.assertLess(erreurs[2], erreurs[1]/3)
            if formulation == "integree":
                for error, limit in zip(erreurs, (5.5e-6, 1.4e-6, 4e-7)):
                    self.assertLess(error, limit)

    def test_statique_petite_charge_anisotrope_sans_reprise(self):
        # Une charge faible ne doit pas provoquer une longue stagnation
        # dans les directions raides avant de changer de mérite.
        for formulation in ('milieu', 'integree'):
            ne = 8
            n = Noyau([0., 0., 0.])
            for i in range(ne+1):
                n.corps(str(i), 1., J, [.508*i/ne, 0., 0.])
            n.liaison('enc', None, 0)
            for i in range(ne):
                n.poutre(str(i), i, i+1, 2.84191e6, 6.40131e5,
                         3.10338, 36.2794, ei3=2.42873, ga3=9.03881e5,
                         formulation=formulation)
            force = .5*(1-math.cos(math.pi/50))*8.896/math.sqrt(2)
            n.effort(ne, [0., force, force], [0., 0., 0.])
            n.statique(tol=1e-8, iters=8, paliers_max=1, strict=True)
            self.assertEqual(n.statique_info()['tentatives'], 1)
            self.assertLessEqual(n.statique_info()['evaluations'], 5)
            self.assertLess(np.max(np.abs(n.residu_statique()[6:])), 1e-8)
            self.assertLess(np.max(np.abs(n.phi())), 1e-8)
            reaction = np.asarray(n.reactions()[0][1])
            np.testing.assert_allclose(reaction[:3], [0., force, force], atol=1e-8, rtol=0)
            np.testing.assert_allclose(reaction[3:], np.cross(n.pose(ne)[0], [0., force, force]),
                                       atol=1e-8, rtol=0)

    def test_statique_funiculaire_discret_et_liaisons(self):
        from scipy.optimize import brentq
        # Un mérite de correction seul pouvait accepter l'ouverture des
        # barres, puis diverger. Référence discrète indépendante, sans erreur
        # de discrétisation d'une caténaire continue.
        for nb in (80, 120):
            ell, mass = 10./nb, 30./nb
            initial = np.array([[8.*k/nb, 0., -3*(1-abs(2*k/nb-1))]
                                for k in range(nb+1)])
            n = Noyau([0., 0., -9.81])
            for k in range(1, nb):
                n.corps(str(k), mass, (1e-9*np.eye(3)).ravel().tolist(), initial[k].tolist())
                n.liaison(f'plan{k}', None, k-1, bloque_t=[1], bloque_r=[0, 1, 2])
            for k in range(nb):
                n.distance(f'barre{k}', None if k == 0 else k-1, None if k == nb-1 else k,
                           initial[0].tolist() if k == 0 else [0., 0., 0.],
                           initial[-1].tolist() if k == nb-1 else [0., 0., 0.], ell)
            vertical = mass*9.81*(np.arange(nb)-(nb-1)/2)
            horizontal = brentq(lambda h: np.sum(ell*h/np.hypot(h, vertical))-8.,
                                1e-6, 1e6, xtol=1e-12)
            segments = ell*np.column_stack((np.full(nb, horizontal), np.zeros(nb), vertical))
            segments /= np.hypot(horizontal, vertical)[:, None]
            reference = np.cumsum(segments, axis=0)[:-1]
            n.statique(tol=1e-10, iters=60, paliers_max=2, strict=True)
            p = np.asarray(n.etat()[1])
            np.testing.assert_allclose(p, reference, atol=1e-9, rtol=0)
            directions = np.diff(np.vstack((initial[0], p, initial[-1])), axis=0)
            np.testing.assert_allclose(np.linalg.norm(directions, axis=1), ell, atol=1e-10, rtol=0)
            directions /= np.linalg.norm(directions, axis=1)[:, None]
            values = dict(n.reactions())
            tensions = np.array([values[f'barre{k}'][0] for k in range(nb)])
            np.testing.assert_allclose(tensions*directions[:, 0], horizontal, atol=1e-8, rtol=0)
            np.testing.assert_allclose(np.diff(tensions[:, None]*directions, axis=0),
                                       np.tile([0., 0., mass*9.81], (nb-1, 1)), atol=1e-8, rtol=0)

    def test_statut_statique_et_mode_strict(self):
        import warnings
        n = Noyau([0., 0., 0.])
        n.corps('libre', 1., J, [0., 0., 0.])
        self.assertIsNone(n.statique_info())
        n.statique(strict=True, paliers_max=1)
        info = n.statique_info()
        self.assertEqual(info['statut'], 'tolerance')
        self.assertTrue(info['tolerance_finale_atteinte'])
        self.assertEqual(info['evaluations'], 1)
        self.assertEqual(info['tentatives'], 1)
        info['statut'] = 'modifie'
        self.assertEqual(n.statique_info()['statut'], 'tolerance')

        # Un corps libre sous force constante n'a AUCUN équilibre. Le seuil
        # de stagnation historique l'acceptait car la force est petite.
        n.effort(0, [5e-7, 0., 0.], [0., 0., 0.])
        before = n.etat()
        with self.assertWarnsRegex(RuntimeWarning, 'stagnation'):
            residual, iterations = n.statique(tol=1e-12, iters=12, paliers_max=1)
        info = n.statique_info()
        self.assertEqual(info['statut'], 'stagnation')
        self.assertFalse(info['tolerance_finale_atteinte'])
        self.assertEqual(info['paliers_stagnation'], 1)
        self.assertEqual(info['echelle_force'], 1.)
        self.assertAlmostEqual(info['residu_libre'], 5e-7)
        self.assertEqual(info['residu_libre'], residual)
        self.assertEqual(info['evaluations'], iterations+1)
        self.assertFalse(info['etat_restaure'])

        with self.assertRaises(RuntimeError):
            n.statique(tol=1e-12, iters=12, paliers_max=1, strict=True)
        info = n.statique_info()
        self.assertEqual(info['statut'], 'echec')
        self.assertTrue(info['strict'])
        self.assertTrue(info['etat_restaure'])
        self.assertFalse(info['tolerance_finale_atteinte'])
        self.assertGreaterEqual(info['tentatives'], 2)
        self.assertEqual(info['paliers_stagnation'], 0)
        self.assertEqual(n.etat(), before)

        # Transformer l'avertissement en exception doit aussi restaurer le
        # modèle, comme tout autre appel statique qui échoue.
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            with self.assertRaises(RuntimeWarning):
                n.statique(tol=1e-12, iters=12, paliers_max=1)
        self.assertTrue(n.statique_info()['etat_restaure'])
        self.assertEqual(n.statique_info()['statut'], 'echec')
        self.assertEqual(n.etat(), before)
        with self.assertRaises(ValueError):
            n.statique(tol=-1.)
        self.assertIsNone(n.statique_info())

    def test_statique_warning_restaure_pose_modifiee(self):
        import warnings
        n = masses()
        n.effort(1, [0., 5e-7, 0.], [0., 0., 0.])
        before = n.etat()
        # Le ressort se détend effectivement avant que la direction neutre
        # chargée fasse stagner le solveur : la restauration n'est pas vide.
        with self.assertWarns(RuntimeWarning):
            n.statique(tol=1e-12, iters=15, paliers_max=1)
        self.assertNotEqual(n.etat(), before)
        self.assertAlmostEqual(n.pose(0)[0][0], .05, places=10)
        n.pose_etat(*before)
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            with self.assertRaises(RuntimeWarning):
                n.statique(tol=1e-12, iters=15, paliers_max=1)
        self.assertEqual(n.etat(), before)
        self.assertTrue(n.statique_info()['etat_restaure'])
        self.assertEqual(n.statique_info()['statut'], 'echec')

    def test_statique_reaction_ne_masque_pas_equilibre(self):
        # Une charge entièrement reprise par le guide ne change pas la
        # flèche libre F/k et ne doit pas changer son seuil de convergence.
        for reaction in (0., 1e12):
            for initial in (0., .09):
                n = Noyau([0., 0., 0.])
                for i in range(2):
                    n.corps(str(i), 1., J, [0., 0., 0.])
                n.liaison('base', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
                n.liaison('guide', None, 1, bloque_t=[0, 1], bloque_r=[0, 1, 2])
                k = np.zeros((12, 12))
                k[np.ix_([2, 8], [2, 8])] = 10*np.array([[1., -1.], [-1., 1.]])
                n.superelement('rappel', [0, 1], k.ravel().tolist())
                n.effort(1, [reaction, 0., 1.], [0., 0., 0.])
                state = list(n.etat())
                state[1][1][2] = initial
                n.pose_etat(*state)
                before = n.etat()
                # Budget insuffisant : refuser le faux succès et restaurer
                # exactement la pose avant de permettre le vrai Newton.
                with self.assertRaises(RuntimeError):
                    n.statique(tol=1e-10, iters=1, paliers_max=1)
                self.assertEqual(n.etat(), before)
                residual, iterations = n.statique(tol=1e-10, iters=20, paliers_max=1, strict=True)
                self.assertGreater(iterations, 0)
                self.assertLess(residual, 1e-10)
                self.assertAlmostEqual(n.pose(1)[0][2], .1, places=12)
                self.assertLess(abs(n.residu_statique()[8]), 1e-10)
                self.assertLess(max(map(abs, n.phi())), 1e-10)

    def test_statique_singuliere_grande_dimension(self):
        # Au-delà du seuil creux : chaque corps garde cinq directions neutres.
        # Le repli doit résoudre l'équilibre sans déplacer ces directions.
        n = Noyau([0., 0., 0.])
        for i in range(30):
            c = n.corps(str(i), 1., J, [.1*i, 0., 0.])
            n.couple(str(i), None, c, [0.,0.,1.], ('ressort', [100.,0.,.1]))
        n.statique(tol=1e-10, iters=30, paliers_max=1)
        self.assertLess(np.max(np.abs(n.residu_statique())), 1e-8)
        for i in range(30):
            p, r = n.pose(i)
            np.testing.assert_allclose(p, [.1*i, 0., 0.], atol=1e-12)
            self.assertAlmostEqual(math.atan2(r[3], r[0]), .1, places=10)

    def test_engrenages_fractionnaires_plusieurs_tours(self):
        from concurrent.futures import ThreadPoolExecutor
        def tourne(ratio):
            n = Noyau(g=[0, 0, 0])
            for i, w in enumerate((ratio * 2, 2)):
                n.corps(str(i), 1, J, [i, 0, 0], w=[0, 0, w])
                n.liaison('p'+str(i), None, i, pa=[i, 0, 0], bloque_r=[0, 1])
            n.engrenage('g', 0, 1, [0, 0, 1], [0, 0, 1], ratio)
            n.simule(20., .02, tous=100)
            for i, vitesse in enumerate((ratio * 2, 2)):
                R = np.array(n.etat()[2][i]).reshape(3, 3)
                angle = vitesse * n.t()
                ref = np.array([[math.cos(angle), -math.sin(angle), 0],
                                [math.sin(angle), math.cos(angle), 0], [0, 0, 1]])
                np.testing.assert_allclose(R, ref, atol=1e-9)
                self.assertAlmostEqual(n.etat()[4][i][2], vitesse, places=10)
            self.assertLess(np.linalg.norm(n.phi()), 1e-10)
            return n.etat()
        rapports = (1.5, -2.3, math.sqrt(2), 2.)
        serie = [tourne(r) for r in rapports]
        with ThreadPoolExecutor(max_workers=4) as pool:
            parallele = list(pool.map(tourne, rapports))
        self.assertEqual(serie, parallele)

    def test_temps_petits_et_reprise(self):
        for methode in ('simule', 'simule_em', 'simule_multirythme'):
            for t0 in (0., 1e-6):
                n = Noyau(g=[0, 0, -9.81])
                for i in range(2):
                    n.corps(str(i), 1, J, [0, 0, 0], v=[1, 0, 0])
                etat = list(n.etat())
                etat[0] = t0
                n.pose_etat(*etat)
                fin = t0 + 5e-13
                kw = dict(rapides=[1], k=2) if methode == 'simule_multirythme' else {}
                getattr(n, methode)(fin, 1e-13, **kw)
                self.assertEqual(n.t(), fin, methode)
                for r, v in zip(n.etat()[1], n.etat()[3]):
                    self.assertAlmostEqual(r[0] / (fin-t0), 1., places=12)
                    self.assertAlmostEqual(v[2] / (fin-t0), -9.81, places=10)

    def test_temps_non_representable_restaure(self):
        for methode in ('simule', 'simule_em', 'simule_multirythme'):
            n = Noyau(g=[0, 0, 0])
            n.corps('a', 1, J, [0, 0, 0])
            etat = list(n.etat())
            etat[0] = 1e16
            n.pose_etat(*etat)
            before = n.etat()
            kw = dict(rapides=[0], k=2) if methode == 'simule_multirythme' else {}
            with self.assertRaisesRegex(RuntimeError, 'progresser le temps'):
                getattr(n, methode)(1e16+4, .1, **kw)
            self.assertEqual(n.etat(), before)

    def test_modes_zero_et_temps_invalides(self):
        for n in (Noyau(), masses()):
            for f in (n.modes, n.modes_complexes):
                self.assertEqual(f(0), [])
                with self.assertRaises(RuntimeError):
                    f(0, t=math.nan)
        self.assertEqual(Noyau().modes(), [])
        self.assertEqual(Noyau().modes_complexes(), [])

    def test_debordement_modal_exception_bornee(self):
        code = '''
from vinkulum import Noyau
import numpy as np
n = Noyau(g=[0, 0, 0])
for i in range(2):
    n.corps(str(i), 1e-300, np.eye(3).ravel().tolist(), [i, 0, 0])
k = np.zeros((12, 12))
k[0, 0] = k[6, 6] = 1e100
k[0, 6] = k[6, 0] = -1e100
n.superelement('s', [0, 1], k.ravel().tolist())
for f in (n.modes, n.modes_complexes, n.spectre, n.bilan_stabilite):
    try:
        f()
    except RuntimeError:
        pass
    else:
        raise AssertionError('débordement accepté')
'''
        p = subprocess.run([sys.executable, '-c', code], capture_output=True,
                           text=True, timeout=10)
        self.assertEqual(p.returncode, 0, p.stderr)

    def test_statique_geometrie_independante_charge(self):
        for masse in (1., 1e10):
            n = Noyau(g=[0, 0, -9.81])
            n.corps('a', masse, J, [0, 0, 0])
            n.liaison('fixe', None, 0)
            etat = list(n.etat())
            etat[1][0][0] = .1
            n.pose_etat(*etat)
            n.statique(paliers_max=1)
            self.assertLess(np.linalg.norm(n.phi()), 1e-10)
            self.assertAlmostEqual(abs(n.reactions()[0][1][2]) / masse, 9.81, places=10)

    def test_statique_et_tangentes_date_servo(self):
        n = Noyau(g=[0, 0, 0])
        n.corps('a', 1, J, [0, 0, 0])
        n.liaison('pivot', None, 0, bloque_r=[0, 1])
        n.couple('servo', None, 0, [0, 0, 1], ('pd', [10, 1, 2, 100]),
                 cible=('lineaire', [0, .5]))
        # Tangentes analytiques de la saturation tanh et de la limite de vitesse.
        K, C, *_ = n.k_c_m_z(t=1.)
        self.assertAlmostEqual(K[5][5], 10 / math.cosh(2.5)**2, places=8)
        self.assertAlmostEqual(C[5][5], 1 / math.cosh(2.5)**2 + math.tanh(2.5)/100, places=8)
        self.assertGreater(abs(n.k_c_m_z(t=0.)[0][5][5]), 9.)
        # Une limite de couple large évite une tangente nulle au Newton initial.
        n2 = Noyau(g=[0, 0, 0])
        n2.corps('a', 1, J, [0, 0, 0])
        n2.liaison('pivot', None, 0, bloque_r=[0, 1])
        n2.couple('servo', None, 0, [0, 0, 1], ('pd', [10, 1, 100, 100]),
                  cible=('lineaire', [0, .5]))
        n2.statique(t=1.)
        R = np.array(n2.etat()[2][0]).reshape(3, 3)
        self.assertAlmostEqual(math.atan2(R[1, 0], R[0, 0]), .5, places=9)
        self.assertEqual(n2.t(), 0.)  # date d'analyse, sans avancer l'horloge

    def test_base_admissible_bras_et_redondance(self):
        from scipy.linalg import null_space
        for L in (1e-5, 1., 1e4, 1e5, 1e8):
            for double in (False, True):
                n = Noyau(g=[0, 0, 0])
                n.corps('a', 1, J, [L, 0, 0])
                n.liaison('rotule', None, 0, bloque_r=[])
                if double:
                    n.liaison('doublon', None, 0, bloque_r=[])
                _, _, _, Z, G = n.k_c_m_z()
                Z, G = np.array(Z), np.array(G)
                self.assertEqual(Z.shape, (6, 3))
                G /= np.linalg.norm(G, axis=1)[:, None]
                reference = null_space(G)
                np.testing.assert_allclose(Z.T @ Z, np.eye(3), atol=2e-14)
                self.assertLess(np.linalg.norm(G @ Z), 1e-13)
                np.testing.assert_allclose(Z @ Z.T, reference @ reference.T, atol=1e-13)

    def test_sensibilites_poutre_isolees(self):
        def modele(ea, parasite):
            n = Noyau(g=[0, 0, 0])
            for i in range(3):
                n.corps(str(i), 1, J, [i, 0, 0])
            n.poutre('p', 0, 1, ea, 20, 30, 40)
            if parasite:
                k = np.kron([[1., -1.], [-1., 1.]], np.eye(6)) * 100
                n.superelement('s', [1, 2], k.ravel().tolist())
                n.contact('sol', 2, [0, 0, 0], 0, normale=[0, 0, 1],
                          origine=[0, 0, .1], k=1000)
            etat = list(n.etat())
            etat[1][1][0] += .1
            etat[1][1][1] += .02
            n.pose_etat(*etat)
            return n
        n = modele(10, True)
        seul = modele(10, False)
        for methode in ('d_residu_poutre', 'd_raideur_poutre'):
            np.testing.assert_allclose(getattr(n, methode)(0, 'ea'),
                                       getattr(seul, methode)(0, 'ea'), atol=1e-13)
        delta = .01
        plus, moins = modele(10+delta, True), modele(10-delta, True)
        ref = (np.array(plus.residu_statique()) - moins.residu_statique()) / (2*delta)
        np.testing.assert_allclose(n.d_residu_poutre(0, 'ea'), ref, atol=1e-10)
        ref_k = (np.array(plus.k_m_z()[0]) - moins.k_m_z()[0]) / (2*delta)
        np.testing.assert_allclose(n.d_raideur_poutre(0, 'ea'), ref_k, atol=1e-6)

    def test_liaison_axes_repetes_et_axes_extremes(self):
        n = Noyau()
        n.corps('a', 1, J, [0, 0, 0])
        for kw in (dict(bloque_r=[0, 1, 2, 0]), dict(bloque_t=[0, 0]),
                   dict(cible_r=([1e308, 0, 0], ('constante', [1])))):
            with self.assertRaises(ValueError):
                n.liaison('invalide', None, 0, **kw)
        self.assertEqual(n.phi(), [])
        with self.assertRaises(ValueError):
            n.couple('invalide', None, 0, [1e308, 0, 0], ('constant', [1]))
        self.assertEqual(n.couples(), [])

    def test_contact_maillage_meme_corps_refuse(self):
        n = Noyau()
        n.corps('a', 1, J, [0, 0, 0])
        mi = n.maillage('m', 0, [0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 1, 2])
        with self.assertRaises(ValueError):
            n.contact('auto', 0, [0, 0, 0], .1, maille=mi)
        self.assertEqual(n.contacts(), [])

    def test_parametres_non_finis_et_inertie_relative(self):
        n = Noyau(g=[0, 0, 0])
        mauvais = np.eye(3) * 1e-15
        mauvais[0, 1] = 1e-16
        with self.assertRaises(ValueError):
            n.corps('asymetrique', 1, mauvais.ravel().tolist(), [0, 0, 0])
        self.assertEqual(len(n.etat()[1]), 0)
        n.corps('valide', 1, (np.eye(3)*1e-15).ravel().tolist(), [0, 0, 0])
        for params in ([1, 1, 0, 1], [1, 1, 1, 0]):
            with self.assertRaises(ValueError):
                n.couple('pd', None, 0, [0, 0, 1], ('pd', params),
                         cible=('constante', [0]))
        for kw in (dict(v_eps=math.inf), dict(d_hat=math.nan),
                   dict(normale=[1e308, 0, 0])):
            with self.assertRaises(ValueError):
                n.contact('invalide', 0, [0, 0, 0], .1, **kw)
        with self.assertRaises(ValueError):
            n.inflow([1e308, 0, 0], 1, 1)
        self.assertEqual(n.contacts(), [])
        self.assertEqual(n.couples(), [])

    def test_oscillateur_super_convergence(self):
        for beta in (0.0, 0.02):
            ref = expm(np.array([[0, 1], [-200, -200 * beta]]) * 0.1) @ [0.1, 0]
            errors = []
            for h in (0.002, 0.001, 0.0005):
                n = masses(beta)
                n.simule(0.1, h)
                _, r, _, v, _, _ = n.etat()
                got = np.array([r[1][0] - r[0][0] - 1, v[1][0] - v[0][0]])
                errors.append(np.linalg.norm(got - ref))
                self.assertLess(abs(r[0][0] + r[1][0] - 1.1), 1e-12)
            self.assertLess(errors[-1], 2e-5)
            self.assertTrue(all(a / b > 3.7 for a, b in zip(errors, errors[1:])), errors)

    def test_super_jacobien(self):
        n = masses(0.02)
        n.liaison("support", None, 0)
        for ggl in (False, True):
            before = n.etat()
            blocs, _ = n.audit_jacobien(0.01, ggl=ggl)
            self.assertEqual(before, n.etat())
            self.assertLess(max(b[1] for b in blocs), 1e-6)

    def test_super_branche_pi(self):
        n = masses()
        # Un second superélément porte une raideur de torsion.
        k = np.zeros((12, 12))
        k[3, 3] = k[9, 9] = 1
        k[3, 9] = k[9, 3] = -1
        n.superelement("torsion", [0, 1], k.ravel().tolist())
        t, r, rot, v, w, vi = n.etat()
        rot[1] = [1, 0, 0, 0, -1, 0, 0, 0, -1]
        n.pose_etat(t, r, rot, v, w, vi)
        self.assertAlmostEqual(abs(n.residu_statique()[9]), math.pi, places=12)
        before = n.etat()
        with self.assertRaises(RuntimeError):
            n.audit_jacobien(.001)
        self.assertEqual(before, n.etat())

    def test_statique_debordement_restaure(self):
        n = Noyau(g=[0, 0, 0])
        n.corps("a", 1, J, [0, 0, 0])
        n.effort(0, [1e308, 1e308, 0], [0, 0, 0])
        before = n.etat()
        with self.assertRaises(RuntimeError):
            n.statique(paliers_max=1, iters=3)
        self.assertEqual(before, n.etat())

    def test_poutre_pi_modes_refuses_sans_panique(self):
        n = Noyau(g=[0, 0, 0])
        n.corps("a", 1, J, [0, 0, 0])
        n.corps("b", 1, J, [1, 0, 0])
        n.poutre("p", 0, 1, 1, 1, 1, 1)
        t, r, rot, v, w, vi = n.etat()
        rot[1] = [1, 0, 0, 0, -1, 0, 0, 0, -1]
        n.pose_etat(t, r, rot, v, w, vi)
        with self.assertRaisesRegex(RuntimeError, "tangente non différentiable"):
            n.modes()

    def test_super_rotation_rigide_amortissement(self):
        n = Noyau(g=[0, 0, 0])
        n.corps("a", 1, J, [0, 0, 0], w=[0, 0, 1])
        n.corps("b", 1, J, [1, 0, 0], v=[0, 1, 0], w=[0, 0, 1])
        k = np.kron([[1., -1.], [-1., 1.]], np.eye(6))
        n.superelement("s", [0, 1], k.ravel().tolist(), beta=1)
        self.assertLess(np.linalg.norm(n.residu_statique()), 1e-13)

    def test_super_domaine_energie(self):
        for action in (lambda n: n.simule_em(.1, .01), lambda n: n.simule(.1, .01, energie=True)):
            n = masses()
            before = n.etat()
            with self.assertRaises((ValueError, RuntimeError)):
                action(n)
            self.assertEqual(before, n.etat())

    def test_energie_moment_corps_distincts_et_reprise(self):
        def chapeau(v):
            x, y, z = v
            return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])

        g = np.array([.4, -.7, -9.81])
        r0 = expm(chapeau([.43, -1.11, 2.07]))
        p = expm(chapeau([-.83, .21, -.62]))
        erreurs = []
        for h in (.02, .01, .005):
            n = Noyau(g.tolist())
            references = []
            for i in range(3):
                j = np.diag(np.array([1., 2., 3.])*(i+1))
                wb = np.eye(3)[i]*(2+i)
                x0, v0 = np.array([i, 0., .2]), np.array([.1*i, -.2, .3])
                n.corps(str(i), i+1., (p.T@j@p).ravel().tolist(), x0.tolist(),
                        rot=(r0@p).ravel().tolist(), v=v0.tolist(), w=(r0@wb).tolist())
                references.append((x0+.4*v0+.5*.4**2*g, v0+.4*g,
                                   r0@expm(.4*chapeau(wb))@p, r0@wb))
            initial = n.etat_precis()
            n.simule_em(0., h)
            self.assertEqual(initial, n.etat_precis())
            n.simule_em(.2, h)
            n.simule_em(.4, h)
            t, xs, rs, vs, ws, _, _ = n.etat_precis()
            self.assertEqual(t, .4)
            errors = []
            for i, (x, v, r, w) in enumerate(references):
                self.assertLess(np.linalg.norm(np.asarray(xs[i])-x), 2e-14)
                self.assertLess(np.linalg.norm(np.asarray(vs[i])-v), 2e-14)
                self.assertLess(np.linalg.norm(np.asarray(ws[i])-w), 2e-13)
                errors.append(np.linalg.norm(np.asarray(rs[i]).reshape(3, 3)-r))
            erreurs.append(errors)
        ordres = np.log2(np.asarray(erreurs[:-1])/erreurs[1:])
        self.assertTrue(np.all((1.99 < ordres) & (ordres < 2.01)), ordres)

    def test_energie_moment_divergence(self):
        for h in (.1, 1., 3., 10., 100.):
            n = Noyau(g=[0, 0, 0])
            n.corps("a", 1, np.diag([1, 2, 3]).ravel().tolist(), [0, 0, 0], w=[1, 2, 3])
            before = n.etat()
            l = np.array(n.moment())
            e = n.energie()
            try:
                n.simule_em(h, h)
            except RuntimeError:
                self.assertEqual(before, n.etat())
            else:
                self.assertLess(abs(n.energie() / e - 1), 1e-12)
                self.assertLess(np.linalg.norm(np.array(n.moment()) - l), 1e-11)

    def test_reprise_multi_echec(self):
        n = Noyau(g=[0, 0, -9.81])
        for i in range(2):
            n.corps(str(i), 1, J, [i, 0, 0])
        # Une partition lente avec liaison non linéaire exige >1 correction.
        n.distance("d", None, 1, [0, 0, 0], [0, 0, 0], l=1)
        before = n.etat()
        with self.assertRaises(RuntimeError):
            n.simule_multirythme(.1, .01, [0], 2, newton_max=1, tol=1e-30)
        self.assertEqual(before, n.etat())
        n.simule(.01, .001)
        self.assertAlmostEqual(n.etat()[1][0][2], -.0004905, places=12)

    def test_multi_partitions_et_contraintes(self):
        def monte():
            n = Noyau(g=[0, 0, -9.81])
            for i in range(4):
                n.corps(str(i), 1, J, [i, 0, -1], v=[0, .2, 0])
                n.distance('d'+str(i), None, i, [i, 0, 0], [0, 0, 0], l=1)
            # Redondance dans une partition, sans liaison traversante.
            n.distance('double', None, 2, [2, 0, 0], [0, 0, 0], l=1)
            return n
        reference = monte()
        reference.simule(.0125, .001)
        for rapides in ([], [0, 1, 2, 3], [0, 2]):
            n = monte()
            n.simule_multirythme(.0125, .001, rapides, 1)
            self.assertEqual(n.t(), .0125)
            for champ in (1, 2, 3, 4):
                np.testing.assert_allclose(n.etat()[champ], reference.etat()[champ], atol=1e-9, rtol=1e-9)
            np.testing.assert_allclose([r[1] for r in n.reactions()],
                                       [r[1] for r in reference.reactions()], atol=1e-8, rtol=1e-8)
            self.assertLess(np.linalg.norm(n.phi()), 1e-10)
            n.simule_multirythme(.014, .001, [1, 3], 2)
            self.assertEqual(n.t(), .014)
            self.assertLess(np.linalg.norm(n.phi()), 1e-10)

    def test_reprise_multi_echec_rapide(self):
        n = Noyau(g=[0, 0, -9.81])
        for i in range(2):
            n.corps(str(i), 1, J, [i, 0, 0])
        n.distance('d', None, 1, [0, 0, 0], [0, 0, 0], l=1)
        before = n.etat()
        with self.assertRaises(RuntimeError):
            n.simule_multirythme(.1, .01, [1], 2, newton_max=1, tol=1e-30)
        self.assertEqual(before, n.etat())
        n.simule(.01, .001)
        self.assertAlmostEqual(n.etat()[1][0][2], -.0004905, places=12)

    def test_multi_echec_apres_macro_pas(self):
        # Un couple déborde seulement après deux macro-pas acceptés.
        def monte(rapide):
            n = Noyau(g=[0, 0, -9.81])
            for i in range(2):
                n.corps(str(i), 1, J, [i, 0, 0])
            n.couple('retarde', None, int(rapide), [0, 0, 1],
                     ('pd', [1e200, 0, 1e200, 1e200]),
                     cible=('table', [0, 0, .002, 0, .003, 1]))
            return n
        for rapide in (False, True):
            n, ref = monte(rapide), monte(rapide)
            ref.simule_multirythme(.002, .001, [1], 2)
            with self.assertRaises(RuntimeError):
                n.simule_multirythme(.004, .001, [1], 2)
            self.assertEqual(n.t(), .002)
            self.assertEqual(n.etat(), ref.etat())

    def test_appariement_modification_loi(self):
        n = Noyau(g=[0, 0, 0])
        for i in range(2):
            n.corps(str(i), 1, J, [i * 1.9, 0, 0])
            n.sphere(i, [0, 0, 0], 1)
        n.appariement(k=1, expo=1)
        n.simule(1e-5, 1e-5)
        n.appariement(k=100, expo=1)
        n.simule(2e-5, 1e-5)
        _, d, f = n.contacts()[0]
        self.assertAlmostEqual(f, 100 * d, places=12)

    def test_gardes_atomiques(self):
        actions = [
            lambda n: n.effort(0, [math.nan, 0, 0], [0, 0, 0]),
            lambda n: n.superelement("bad", [0, 1], [math.nan] * 144),
            lambda n: n.superelement("bad", [0, 1], [0.] * 144, alpha=1),
            lambda n: n.simule(.1, .01, pas_contact=(0, 0, .1)),
            lambda n: n.simule(.1, .01, tol=math.nan),
            lambda n: n.simule(.1, .01, bornes=(2, 1)),
            lambda n: n.simule(.1, .01, newton_max=0),
            lambda n: n.statique(tol=math.nan),
            lambda n: n.couple("vide", None, 0, [0, 0, 1], ("ressort", [])),
        ]
        for method in ("simule", "simule_em", "simule_multirythme"):
            for bad in (math.nan, math.inf, -.1):
                def action(n, method=method, bad=bad):
                    kw = dict(rapides=[0], k=2) if method == "simule_multirythme" else {}
                    getattr(n, method)(bad, .01, **kw)
                actions.append(action)
        for action in actions:
            n = masses()
            before = n.etat()
            with self.assertRaises(ValueError):
                action(n)
            self.assertEqual(before, n.etat())

    def test_pas_contact_termine(self):
        code = '''from vinkulum import Noyau
n=Noyau(g=[0,0,-9.81]);n.corps("a",1,[1,0,0,0,1,0,0,0,1],[1,0,0])
n.distance("d",None,0,[0,0,0],[0,0,0],l=1)
try: n.simule(.1,.01,newton_max=1,tol=1e-30,pas_contact=(.01,.01,.01))
except RuntimeError: pass
print("termine")
'''
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=5)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout.strip(), "termine")


if __name__ == "__main__":
    unittest.main()
