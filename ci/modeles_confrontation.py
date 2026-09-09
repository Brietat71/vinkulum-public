"""Modèles de confrontation : paramètres repris des bancs Vinkulum existants.
Aucune lecture du code du solveur MBDyn.
"""
import numpy as np
from vinkulum import Noyau

def six_barres():
    lg, mm, om = 1.0, 1.0, 1.0
    j1, j2 = mm * lg * lg / 12.0, 1e-9
    ref = np.array([-0.049323, -0.835610])
    pos = {0: [0, 0, lg / 2], 1: [lg, 0, lg / 2], 2: [2 * lg, 0, lg / 2],
           101: [lg / 2, 0, lg], 102: [3 * lg / 2, 0, lg]}
    vit = {k: ([lg / 2 * om, 0, 0] if k < 3 else [lg * om, 0, 0]) for k in pos}
    omg = {k: ([0, om, 0] if k < 3 else [0, 0, 0]) for k in pos}
    ine = {0: (j1, j1, j2), 1: (j1, j1, j2), 2: (j1, j1, j2),
           101: (j2, j1, j1), 102: (j2, j1, j1)}
    n = Noyau([0.0, 0.0, -9.81])
    idx = {}
    for k in (0, 1, 2, 101, 102):
        a, b, c = ine[k]
        idx[k] = n.corps(f"c{k}", mm, [a, 0, 0, 0, b, 0, 0, 0, c], pos[k],
                         v=vit[k], w=omg[k])
    for nom, a, b, pt, bt, br in [
            ("j0", None, 0, [0, 0, 0], [0, 1, 2], [0, 2]),
            ("j1", None, 1, [lg, 0, 0], [0, 2], []),
            ("j101", 0, 101, [0, 0, lg], [0, 1, 2], [0, 2]),
            ("j1101", 1, 101, [lg, 0, lg], [0, 1, 2], [0, 2]),
            ("j2", None, 2, [2 * lg, 0, 0], [0, 2], []),
            ("j102", 1, 102, [lg, 0, lg], [0, 1, 2], [0, 2]),
            ("j1102", 2, 102, [2 * lg, 0, lg], [0, 1, 2], [0, 2])]:
        pa = list(pt) if a is None else list(np.asarray(pt) - np.asarray(pos[a]))
        n.liaison(nom, None if a is None else idx[a], idx[b],
                  pa=pa, bloque_t=bt, bloque_r=br)
    phi0 = max(abs(x) for x in n.phi())
    phid0 = max(abs(x) for x in n.phi_dot())
    return n, [idx[k] for k in (0, 1, 2, 101, 102)]

def spatial():
    L_C, L_R, OM = 0.08, 0.30, 6.0
    A = np.array([0.0, 0.1, 0.12])
    B = A + np.array([0.0, 0.0, L_C])
    C = np.array([np.sqrt(L_R ** 2 - 0.1 ** 2 - (0.12 + L_C) ** 2), 0.0, 0.0])
    REF = {"rod": np.array([0.1273241, 0.0142809, 0.0780041]),
           "block": np.array([0.2546481, 0.0, 0.0])}
    REF_W = np.array([1.92, -0.96, 0.48])          # omega de la bielle, .mov de MBDyn a t = 0
    def rc(v1, v2):                # convention MBDyn : 1er axe exact, 2e orthogonalise
        e1 = np.asarray(v1, float) / np.linalg.norm(v1)
        e2 = np.asarray(v2, float) - np.dot(v2, e1) * e1
        return e1, e2 / np.linalg.norm(e2)
    e1, e2 = rc([1, 0, 0], [0, 0, 1])
    r_cr = np.column_stack([e1, e2, np.cross(e1, e2)])
    f2, f3 = rc(B - C, [1, 0, 0])                  # axe 2 = C -> B, le long de la bielle
    r_ro = np.column_stack([np.cross(f2, f3), f2, f3])
    cm_r = C + L_R / 2 * f2
    w_cr = np.array([OM, 0.0, 0.0])
    d, ua, ub = B - C, np.array([1.0, 0.0, 0.0]), r_ro[:, 0]
    m, y = np.zeros((4, 4)), np.zeros(4)
    m[:3, 0] = [1, 0, 0]
    for k in range(3):
        ek = np.zeros(3)
        ek[k] = 1.0
        m[:3, 1 + k] = np.cross(ek, d)             # v_C + omega x (B - C) = v_B
    m[3, 1:] = np.cross(ub, ua)                    # d/dt (u_a . u_b) = 0
    y[:3] = np.cross(w_cr, B - A)
    s = np.linalg.solve(m, y)
    v_bl = np.array([s[0], 0.0, 0.0])
    w_ro = s[1:]
    v_ro = v_bl + np.cross(w_ro, cm_r - C)
    def monte(cardan=True):
        n = Noyau([0.0, 0.0, -9.81])
        cr = n.corps("crank", 0.12, [1e-4, 0, 0, 0, 1e-5, 0, 0, 0, 1e-4], list(A),
                     rot=list(r_cr.flatten()), w=list(w_cr))
        ro = n.corps("rod", 0.5, [4e-3, 0, 0, 0, 4e-4, 0, 0, 0, 4e-3], list(cm_r),
                     rot=list(r_ro.flatten()), v=list(v_ro), w=list(w_ro))
        bl = n.corps("block", 2.0, [1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4], list(C),
                     v=list(v_bl))
        n.liaison("A", None, cr, pa=list(A), bloque_t=[0, 1, 2], bloque_r=[1, 2])
        n.liaison("B", cr, ro, pa=list(r_cr.T @ (B - A)), bloque_t=[0, 1, 2], bloque_r=[])
        n.liaison("C", bl, ro, pa=[0, 0, 0], bloque_t=[0, 1, 2], bloque_r=[])
        if cardan:
            n.cardan("croix", bl, ro, [1, 0, 0], [1, 0, 0])
        n.liaison("gliss", None, bl, pa=list(C), bloque_t=[1, 2], bloque_r=[0, 1, 2])
        return n, (ro, bl)
    n, _ = monte()
    return n, [0, 1, 2]
