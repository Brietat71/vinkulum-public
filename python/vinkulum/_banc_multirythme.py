"""Modèle commun au banc physique et à la mesure isolée du multi-rythme."""
from math import pi
import numpy as np
from vinkulum import Noyau

G = 9.81
I3 = np.eye(3)


def monte(n_ch=8, enjambe=False):
    L = 0.25                              # chaîne : 8 maillons de 0,25 m
    m_c, m_b, l1, lb = 0.05, 0.01, 0.05, 0.05
    f_ax, f_mont = 300.0, 5.0
    ea_r = m_b * lb * (2 * pi * f_ax) ** 2         # poutre raide C–B : mode axial à f_ax
    ei_r = (2 * pi * f_ax) ** 2 * m_b * lb ** 3 / 3.0
    ea_s = (m_c + m_b) * l1 * (2 * pi * f_mont) ** 2   # montage souple A–C à f_mont
    ei_s = (2 * pi * f_mont) ** 2 * (m_c + m_b) * l1 ** 3 / 3.0
    dv, w0 = 0.5, 0.5                              # excitation axiale de B ; rotation d'ensemble

    N = Noyau([0.0, -G, 0.0])
    ch = []
    for i in range(n_ch):
        x = (i + 0.5) * L
        ch.append(N.corps(f"m{i}", 0.2, list((0.2 * L * L / 12.0 * I3).ravel()), [x, 0.0, 0.0],
                          v=[0.0, w0 * x, 0.0], w=[0.0, 0.0, w0]))
    N.liaison("piv", None, ch[0], pa=[0.0, 0.0, 0.0], bloque_r=[0, 1])
    for i in range(n_ch - 1):
        N.liaison(f"l{i}", ch[i], ch[i + 1], pa=[0.5 * L, 0.0, 0.0], bloque_r=[0, 1])
    xa, xc, xb = n_ch * L, n_ch * L + l1, n_ch * L + l1 + lb
    c = N.corps("C", m_c, list((1e-6 * I3).ravel()), [xc, 0.0, 0.0], v=[0.0, w0 * xc, 0.0], w=[0, 0, w0])
    b = N.corps("B", m_b, list((1e-7 * I3).ravel()), [xb, 0.0, 0.0], v=[dv, w0 * xb, 0.0], w=[0, 0, w0])
    N.poutre("souple", ch[-1], c, ea_s, 1e2 * ea_s, 1e2 * ei_s, ei_s)
    N.poutre("raide", c, b, ea_r, 1e2 * ea_r, 1e2 * ei_r, ei_r)
    if enjambe:
        N.liaison("faux", ch[-1], c, pa=[0.5 * L + l1, 0.0, 0.0], bloque_t=[1, 2], bloque_r=[])
    return N, ch[-1], c, b

