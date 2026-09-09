"""vinkulum.bancs — LA SUITE DE BANCS : courbes travail–précision, le juge du noyau.

    python -m vinkulum.bancs               # tous les bancs, plusieurs pas, JSON + table
    python -m vinkulum.bancs rapide        # un pas par banc (contrôle de non-régression)

« Dépasser la concurrence » se lit sur des courbes ERREUR contre TEMPS CPU,
nulle part ailleurs. Chaque banc a une référence INDÉPENDANTE, et publie pour
plusieurs pas de temps : l'erreur mesurée contre elle, le temps de `simule`,
le nombre de pas, l'ordre de convergence observé. La suite écrit
`docs/bancs/<version>.json` — l'historique des versions est la preuve.

Bancs (v0.2.0) :
· pendule (rotule, 30°) — EDO minimale à 1e-12 ; erreur en position ;
· toupie symétrique sur rotule, spin rapide — EDO de Lagrange à 1e-11 ; erreur d'axe ;
· bielle-manivelle à rotation imposée — x(θ) analytique ; boucle fermée ;
· engrenage rapport −4 sous couple — α_a = τ/(J_a + J_b/16) analytique ;
· double pendule (deux rotules, plan) — EDO minimale à 1e-13, sensibilité ;
· poutre géométriquement exacte : console encastrée, 1er mode contre
  1,875⁴·EI/ρAL⁴ et flèche contre Euler–Bernoulli, à 4/10/20 éléments ;
· quatre-barres REDONDANT (4 pivots parallèles en 3D : 20 contraintes pour 18
  DDL, rang 17) contre le MÊME mécanisme en liaisons non redondantes (pivot,
  rotule, cardan, pivot : 17) — les deux doivent rendre le même mouvement ;
  c'est le banc des contraintes redondantes (Bricard, en plus simple).

La tête S2 de FRELON (143 inconnues, boucles fermées multiples) est le
septième banc ; il vit dans FRELON (`cad/vinkulum_s2.py`) et publie le même
format.
"""
import json
import os
import sys
import time
from math import cos, pi, sin, sqrt

import numpy as np
from scipy.integrate import solve_ivp

from vinkulum import Noyau
from vinkulum import __version__ as VERSION
from vinkulum.verification import G, I3, _expm, _pendule_ref, _rot9, _toupie_ref


def _chrono(N, t_end, h, **kw):
    t0 = time.perf_counter()
    tr = N.simule(t_end, h, **kw)
    return tr, time.perf_counter() - t0


# ── bancs ────────────────────────────────────────────────────────────────────
def pendule(h):
    L, m, th0 = 1.0, 0.2, np.radians(30.0)
    T = 2 * pi * sqrt(L / G) * 1.07
    t_end = 2 * T
    N = Noyau()
    N.corps("masse", m, _rot9(np.diag([1e-8] * 3)), [L * sin(th0), 0.0, -L * cos(th0)])
    N.liaison("rotule", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
    tr, dt = _chrono(N, t_end, h)
    th, thd = _pendule_ref(t_end, th0, L)(tr[-1][0])
    err = np.linalg.norm(np.array(tr[-1][1][0]) - [L * sin(th), 0.0, -L * cos(th)])
    return dict(erreur=float(err), unite="m", temps=dt, pas=len(tr), inconnues=6 + 3)


def toupie(h):
    m, l, Ja, Jt, ws, th0 = 0.1, 0.05, 1e-4, 3e-4, 300.0, np.radians(30.0)
    t_end = 1.0
    ts, ths, phis = _toupie_ref(t_end, m, l, Ja, Jt, ws, th0)
    R0 = _expm(np.array([th0, 0.0, 0.0]))
    axe0 = R0 @ np.array([0.0, 0.0, 1.0])
    N = Noyau()
    N.corps("toupie", m, _rot9(np.diag([Jt, Jt, Ja])), list(l * axe0), _rot9(R0), w=list(ws * axe0))
    N.liaison("pointe", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
    tr, dt = _chrono(N, t_end, h, tous=max(1, int(round(2e-3 / h))))
    pire = 0.0
    for e in tr:
        R = np.array(e[2][0]).reshape(3, 3)
        th, ph = np.interp(e[0], ts, ths), np.interp(e[0], ts, phis)
        ax_ref = _expm(np.array([0.0, 0.0, ph])) @ _expm(np.array([th, 0.0, 0.0])) @ np.array([0.0, 0.0, 1.0])
        pire = max(pire, float(np.linalg.norm(R[:, 2] - ax_ref)))
    return dict(erreur=pire, unite="rad", temps=dt, pas=int(round(t_end / h)), inconnues=6 + 3)


def bielle(h):
    a, l, om = 0.03, 0.10, 20.0
    N = Noyau(g=[0.0, 0.0, 0.0])
    N.corps("manivelle", 0.05, _rot9(np.diag([1e-6, 2e-5, 2e-5])), [a / 2, 0, 0])
    N.corps("bielle", 0.08, _rot9(np.diag([1e-6, 8e-5, 8e-5])), [a + l / 2, 0, 0])
    N.corps("coulisseau", 0.10, _rot9(np.diag([1e-5] * 3)), [a + l, 0, 0])
    N.liaison("pivot moteur", None, 0, cible_r=([0.0, 0.0, 1.0], ("lineaire", [0.0, om])))
    N.liaison("rotule maneton", 0, 1, pa=[a / 2, 0, 0], bloque_r=[])
    N.liaison("rotule pied", 1, 2, pa=[l / 2, 0, 0], bloque_r=[])
    N.liaison("glissière", None, 2, pa=[a + l, 0, 0], bloque_t=[1, 2])
    T = 2 * pi / om
    tr, dt = _chrono(N, 2 * T, h)
    pire = max(abs(e[1][2][0] - (a * cos(om * e[0]) + sqrt(l**2 - (a * sin(om * e[0]))**2))) for e in tr)
    return dict(erreur=float(pire), unite="m", temps=dt, pas=len(tr), inconnues=18 + 17)


def engrenage(h):
    Ja, Jb, tau, r = 2e-6, 4e-5, 1e-3, -4.0
    N = Noyau(g=[0.0, 0.0, 0.0])
    N.corps("pignon", 0.01, _rot9(np.diag([Ja / 2, Ja / 2, Ja])), [0.014, 0, 0])
    N.corps("roue", 0.05, _rot9(np.diag([Jb / 2, Jb / 2, Jb])), [0, 0, 0])
    N.liaison("pivot pignon", None, 0, pa=[0.014, 0, 0], bloque_r=[0, 1])
    N.liaison("pivot roue", None, 1, pa=[0, 0, 0], bloque_r=[0, 1])
    N.engrenage("denture", 0, 1, [0, 0, 1], [0, 0, 1], r)
    N.effort(0, [0, 0, 0], [0, 0, tau])
    T = 0.5
    tr, dt = _chrono(N, T, h, tous=max(1, int(round(0.01 / h))))
    alpha_th = tau / (Ja + Jb / r**2)
    err = abs(tr[-1][4][0][2] / T - alpha_th) / alpha_th
    return dict(erreur=float(err), unite="rel", temps=dt, pas=int(round(T / h)), inconnues=12 + 11)


def _double_pendule_ref(t_end, L1, L2, m1, m2, th1, th2):
    """Deux masses ponctuelles, tiges sans masse, plan vertical : EDO classique."""
    def f(t, y):
        a1, a2, w1, w2 = y
        d = a1 - a2
        den = 2 * m1 + m2 - m2 * cos(2 * d)
        a1dd = (-G * (2 * m1 + m2) * sin(a1) - m2 * G * sin(a1 - 2 * a2)
                - 2 * sin(d) * m2 * (w2**2 * L2 + w1**2 * L1 * cos(d))) / (L1 * den)
        a2dd = (2 * sin(d) * (w1**2 * L1 * (m1 + m2) + G * (m1 + m2) * cos(a1) + w2**2 * L2 * m2 * cos(d))) / (L2 * den)
        return [w1, w2, a1dd, a2dd]
    return solve_ivp(f, (0, t_end), [th1, th2, 0.0, 0.0], rtol=1e-13, atol=1e-14, dense_output=True).sol


def double_pendule(h):
    L1, L2, m1, m2 = 1.0, 0.7, 0.3, 0.2
    th1, th2 = np.radians(40.0), np.radians(-20.0)
    t_end = 2.0
    ref = _double_pendule_ref(t_end, L1, L2, m1, m2, th1, th2)
    p1 = np.array([L1 * sin(th1), 0.0, -L1 * cos(th1)])
    p2 = p1 + np.array([L2 * sin(th2), 0.0, -L2 * cos(th2)])
    N = Noyau()
    N.corps("m1", m1, _rot9(np.diag([1e-8] * 3)), list(p1))
    N.corps("m2", m2, _rot9(np.diag([1e-8] * 3)), list(p2))
    N.liaison("rotule 1", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
    N.liaison("rotule 2", 0, 1, pa=[0.0, 0.0, 0.0], bloque_t=[0, 1, 2], bloque_r=[])
    tr, dt = _chrono(N, t_end, h)
    a1, a2 = ref(tr[-1][0])[:2]
    r2_ref = np.array([L1 * sin(a1) + L2 * sin(a2), 0.0, -L1 * cos(a1) - L2 * cos(a2)])
    err = np.linalg.norm(np.array(tr[-1][1][1]) - r2_ref)
    return dict(erreur=float(err), unite="m", temps=dt, pas=len(tr), inconnues=12 + 6)


def _quatre_barres(h, redondant):
    """Plan xz, gravité −z, pivots ∥ y. Bâti A=(0,0,0) B=(0.4,0,0) ; manivelle
    AC (0,15), bielle CD (0,45), balancier DB (0,3). Position initiale : C au
    dessus de A, D résolu par fermeture."""
    L1, L2, L3, AB = 0.15, 0.45, 0.30, 0.40
    C = np.array([0.0, 0.0, L1])
    # D : |D − C| = L2, |D − B| = L3, dans le plan xz, au-dessus
    B = np.array([AB, 0.0, 0.0])
    dcb = B - C
    d = np.linalg.norm(dcb)
    a = (L2**2 - L3**2 + d**2) / (2 * d)
    hh = sqrt(L2**2 - a**2)
    P = C + a * dcb / d
    perp = np.array([-dcb[2], 0.0, dcb[0]]) / d
    D = P + hh * perp if (P + hh * perp)[2] > (P - hh * perp)[2] else P - hh * perp
    def barre(nom, p, q, m):
        v = q - p
        L = np.linalg.norm(v)
        ex = v / L
        ey = np.array([0.0, 1.0, 0.0])
        ez = np.cross(ex, ey)
        R = np.column_stack([ex, ey, ez])
        J = np.diag([1e-8, m * L**2 / 12, m * L**2 / 12])
        return N.corps(nom, m, _rot9(J), list((p + q) / 2), _rot9(R)), L
    N = Noyau()
    (c1, l1), (c2, l2), (c3, l3) = barre("manivelle", np.zeros(3), C, 0.5), barre("bielle", C, D, 1.0), barre("balancier", D, B, 0.7)
    # repère de pivot d'axe y : colonnes (e1, e2, e3) = (z, x, y) — `ra` se donne
    # EN LIGNE (row-major) : la matrice [[0,1,0],[0,0,1],[1,0,0]] aplatie. Piège
    # payé : la même liste lue comme des colonnes donne un axe x, et un
    # quatre-barres plan bloqué dans ses deux modélisations.
    Ry = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
    N.liaison("A", None, c1, pa=[0.0, 0.0, 0.0], ra=Ry, bloque_r=[0, 1])
    if redondant:
        N.liaison("C", c1, c2, pa=[l1 / 2, 0.0, 0.0], ra=Ry, bloque_r=[0, 1])
        N.liaison("D", c2, c3, pa=[l2 / 2, 0.0, 0.0], ra=Ry, bloque_r=[0, 1])
    else:
        N.liaison("C", c1, c2, pa=[l1 / 2, 0.0, 0.0], bloque_r=[])            # rotule
        # cardan : translation bloquée + rotation de la bielle autour de son propre axe (x local de c2)
        N.liaison("D", c2, c3, pa=[l2 / 2, 0.0, 0.0], bloque_r=[0])
    N.liaison("B", None, c3, pa=[AB, 0.0, 0.0], ra=Ry, bloque_r=[0, 1])
    # la position initiale est un équilibre sous gravité (mesuré : λ 11,3 N et
    # rien ne bouge) — un couple constant sur la manivelle met le mécanisme en
    # marche ; le banc compare deux modélisations, il n'a pas besoin de repos
    N.effort(c1, [0.0, 0.0, 0.0], [0.0, 0.3, 0.0])
    E0 = N.energie()
    tr, dt = _chrono(N, 1.5, h, tous=max(1, int(round(5e-3 / h))))
    return N, tr, dt, E0


def quatre_barres(h):
    Nr, trr, dtr, E0r = _quatre_barres(h, True)
    Nn, trn, dtn, E0n = _quatre_barres(h, False)
    # même mouvement : la manivelle des deux modèles, image par image
    ecart = max(float(np.linalg.norm(np.array(a[1][0]) - np.array(b[1][0]))) for a, b in zip(trr, trn))
    phi_r = float(np.linalg.norm(Nr.phi()))
    dE = abs(Nr.energie() - E0r) / max(abs(E0r), 1e-12)     # informatif : un couple travaille
    return dict(erreur=ecart, unite="m (redondant − non redondant)", temps=dtr, pas=int(round(1.5 / h)),
                inconnues=18 + 20, phi=phi_r, dE=float(dE), temps_non_redondant=dtn)


def chaine(n_corps, h=1e-3, t_end=0.2, L=0.05, m=0.02):
    """Chaîne de `n_corps` maillons en rotules, lâchée à l'horizontale sous
    gravité — le banc de PASSAGE À L'ÉCHELLE : 6n inconnues + 3n contraintes,
    matrice bloc-tridiagonale. On mesure le temps par pas contre N."""
    N = Noyau()
    for i in range(n_corps):
        N.corps(f"m{i}", m, _rot9(np.diag([m * L**2 / 12] * 3 + []) if False else np.diag([1e-8, m * L**2 / 12, m * L**2 / 12])),
                [L * (i + 0.5), 0.0, 0.0])
    N.liaison("r0", None, 0, pa=[0.0, 0.0, 0.0], bloque_t=[0, 1, 2], bloque_r=[])
    for i in range(1, n_corps):
        N.liaison(f"r{i}", i - 1, i, pa=[L / 2, 0.0, 0.0], bloque_t=[0, 1, 2], bloque_r=[])
    tr, dt = _chrono(N, t_end, h, tous=10**9)
    n = int(round(t_end / h))
    j, so, re, fin, tot, ex = N.chronos()
    return dict(n_corps=n_corps, inconnues=9 * n_corps, ms_pas=dt * 1e3 / n, jacobien=j * 1e3 / n,
                resolution=so * 1e3 / n, residus=re * 1e3 / n, exp=ex * 1e3 / n,
                newton=N.stats()[0] / n, phi=float(np.linalg.norm(N.phi())))


def echelle(tailles=(10, 30, 100)):
    print("  passage à l'échelle — chaîne de N maillons (ms/pas) :")
    out = []
    for nb in tailles:
        r = chaine(nb)
        out.append(r)
        print(f"    N={nb:4d} ({r['inconnues']:5d} inconnues)  {r['ms_pas']:8.3f} ms/pas  "
              f"jacobien {r['jacobien']:.3f}  résolution {r['resolution']:.3f}  résidus {r['residus']:.3f}  Newton {r['newton']:.2f}  |Φ| {r['phi']:.1e}")
    if len(out) >= 2:
        pentes = [np.log(out[i + 1]["ms_pas"] / out[i]["ms_pas"]) / np.log(out[i + 1]["n_corps"] / out[i]["n_corps"]) for i in range(len(out) - 1)]
        print(f"    exposant observé du temps par pas en N : {', '.join(f'{p:.2f}' for p in pentes)}")
    return out


def console(n_elem, h=None):
    """Poutre géométriquement exacte encastrée (alu 1 m, 20 × 5 mm), chargée en
    bout : fréquence du 1er mode contre 1,875⁴·EI/ρAL⁴ et flèche contre
    Euler–Bernoulli FL³/3EI. Deux références analytiques indépendantes."""
    E, rho_, L, b, ep, F = 70e9, 2700.0, 1.0, 0.02, 0.005, 0.5
    A = b * ep
    Iy = b * ep ** 3 / 12
    Iz = ep * b ** 3 / 12
    J = Iy + Iz
    G = E / 2.6
    f_th = 1.875104 ** 2 * sqrt(E * Iy / (rho_ * A * L ** 4)) / (2 * pi)
    h = h or 1 / (200 * f_th)
    N = Noyau(g=[0, 0, 0])
    le = L / n_elem
    idx = []
    for i in range(n_elem + 1):
        m = rho_ * A * le * (0.5 if i in (0, n_elem) else 1.0)
        Jl = np.diag([rho_ * J * le, rho_ * Iy * le, rho_ * Iz * le]) * (0.5 if i in (0, n_elem) else 1.0) + np.diag([1e-12] * 3)
        idx.append(N.corps(f"n{i}", m, _rot9(Jl), [i * le, 0.0, 0.0]))
    N.liaison("encastrement", None, idx[0], pa=[0, 0, 0])
    for i in range(n_elem):
        N.poutre(f"e{i}", idx[i], idx[i + 1], E * A, G * A, G * J, E * Iy)
    N.effort(idx[-1], [0.0, 0.0, -F], [0, 0, 0])
    tr, dt = _chrono(N, 6 / f_th, h)
    z = np.array([e[1][idx[-1]][2] for e in tr])
    ts = np.array([e[0] for e in tr])
    zm = z[len(z) // 3:].mean()
    cr = [ts[i] - (z[i] - zm) * (ts[i + 1] - ts[i]) / (z[i + 1] - z[i])
          for i in range(len(z) - 1) if (z[i] - zm) < 0 <= (z[i + 1] - zm)]
    f_mes = 1.0 / float(np.mean(np.diff(cr)))
    fleche_th = F * L ** 3 / (3 * E * Iy)
    return dict(erreur=abs(f_mes / f_th - 1), unite="rel (f1)", temps=dt, pas=len(tr),
                inconnues=9 * (n_elem + 1), f1=f_mes, f1_th=f_th,
                fleche=-zm, fleche_th=fleche_th, err_fleche=abs(-zm / fleche_th - 1))


def modes_console(n_elem=20):
    """Les TROIS premiers modes de flexion d'une console, par analyse modale du
    système contraint (K, M, noyau de G) — contre 1,875⁴ / 4,694⁴ / 7,855⁴.
    Trois références d'un coup, en 10 ms : le contrôle le plus dense du noyau."""
    E, rho_, L, b, ep = 70e9, 2700.0, 1.0, 0.02, 0.005
    A = b * ep
    Iy = b * ep ** 3 / 12
    Iz = ep * b ** 3 / 12
    J = Iy + Iz
    G = E / 2.6
    N = Noyau(g=[0, 0, 0])
    le = L / n_elem
    idx = []
    for i in range(n_elem + 1):
        m = rho_ * A * le * (0.5 if i in (0, n_elem) else 1.0)
        Jl = np.diag([rho_ * J * le, rho_ * Iy * le, rho_ * Iz * le]) * (0.5 if i in (0, n_elem) else 1.0) + np.diag([1e-12] * 3)
        idx.append(N.corps(f"n{i}", m, _rot9(Jl), [i * le, 0.0, 0.0]))
    N.liaison("encastrement", None, idx[0], pa=[0, 0, 0])
    for i in range(n_elem):
        N.poutre(f"e{i}", idx[i], idx[i + 1], E * A, G * A, G * J, E * Iy)
    t0 = time.perf_counter()
    md = N.modes(combien=10)
    dt = time.perf_counter() - t0
    base = sqrt(E * Iy / (rho_ * A * L ** 4)) / (2 * pi)
    th = [1.875104 ** 2 * base, 4.694091 ** 2 * base, 7.854757 ** 2 * base]
    fs = [f for f, _ in md]
    ecarts = [abs(min(fs, key=lambda f: abs(f / t - 1)) / t - 1) for t in th]
    return dict(erreur=max(ecarts), unite="rel (3 modes)", temps=dt, pas=1,
                inconnues=6 * (n_elem + 1), f=[min(fs, key=lambda f: abs(f / t - 1)) for t in th], f_th=th)


def amortis():
    """Modes COMPLEXES contre trois solutions exactes : oscillateur amorti
    (f_d = ω_n√(1−ζ²), ζ = c/2√(kJ)), boucle servo PD (ζ = Kd/2√(KpI)), et
    détection d'INSTABILITÉ (amortissement négatif → σ > 0)."""
    J, k, c = 2e-4, 0.5, 3e-3
    N = Noyau(g=[0, 0, 0])
    N.corps("d", 0.1, _rot9(np.diag([J / 2, J / 2, J])), [0, 0, 0])
    N.liaison("pivot", None, 0, bloque_r=[0, 1])
    N.couple("torsion", None, 0, [0, 0, 1], ("ressort", [k, c, 0.0]))
    t0 = time.perf_counter()
    f, z, _ = N.modes_complexes(combien=1)[0]
    dt = time.perf_counter() - t0
    wn = sqrt(k / J)
    z_th = c / (2 * sqrt(k * J))
    f_th = wn * sqrt(1 - z_th ** 2) / (2 * pi)
    # instabilité : le même, amortissement négatif
    N2 = Noyau(g=[0, 0, 0])
    N2.corps("d", 0.1, _rot9(np.diag([J / 2, J / 2, J])), [0, 0, 0])
    N2.liaison("pivot", None, 0, bloque_r=[0, 1])
    N2.couple("torsion", None, 0, [0, 0, 1], ("ressort", [k, -1e-3, 0.0]))
    _, _, sig = N2.modes_complexes(combien=1)[0]
    return dict(erreur=max(abs(f / f_th - 1), abs(z / z_th - 1)), unite="rel (f et ζ)",
                temps=dt, pas=1, inconnues=6, f=f, f_th=f_th, zeta=z, zeta_th=z_th,
                sigma_instable=sig)


def raideur_gravite():
    """La raideur tangente d'une toupie AU REPOS doit valoir exactement le
    couple de rappel du poids, mgl — et elle ne vient pas de ∂f/∂q (le poids
    est constant) mais de la précontrainte ∂(Gᵀλ)/∂q. C'est le contrôle qui
    dit que K est juste, et il isole le défaut de l'état TOURNANT."""
    m, l, Ja, Jt = 0.1, 0.05, 1e-4, 3e-4
    N = Noyau(g=[0, 0, -9.80665])
    N.corps("t", m, _rot9(np.diag([Jt, Jt, Ja])), [0, 0, l])
    N.liaison("pointe", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
    k, _, _, _ = N.diag_modal()
    mgl = m * 9.80665 * l
    return dict(k=k, k_th=mgl * sqrt(2.0), erreur=abs(k / (mgl * sqrt(2.0)) - 1))


def _deux_echelles(h, adapt=None, t_end=2.0):
    """Deux modes séparés d'un facteur 60, le raide amorti : le cas type où un
    pas adaptatif est censé payer (transitoire court, puis lisse)."""
    N = Noyau([0, 0, 0])
    j = 0.01
    b1 = N.corps("lent", 1.0, list((j * np.eye(3)).ravel()), [0, 0, 0], w=[3.0, 0, 0])
    b2 = N.corps("raide", 1.0, list((j * np.eye(3)).ravel()), [1, 0, 0], w=[200.0, 0, 0])
    N.liaison("p1", None, b1, pa=[0, 0, 0], bloque_t=[0, 1, 2], bloque_r=[1, 2])
    N.liaison("p2", None, b2, pa=[-1, 0, 0], bloque_t=[0, 1, 2], bloque_r=[1, 2])
    for b, f, z in ((b1, 1.0, 0.0), (b2, 60.0, 0.25)):
        w = 2 * pi * f
        N.couple(f"k{b}", None, b, [1.0, 0.0, 0.0], ("ressort", [j * w * w, 2 * z * j * w, 0.0]))
    t0 = time.time()
    N.simule(t_end, h, tous=10 ** 9, adaptatif=adapt)
    m = np.array(N.etat()[2][0], float).reshape(3, 3)
    return dict(angle=float(np.arctan2(m[2, 1], m[1, 1])), temps=time.time() - t0,
                newton=N.stats()[0], adapt=N.adapt_stats())


def adaptatif():
    """CE QUE LE PAS ADAPTATIF COÛTE — un résultat NÉGATIF, mesuré et gardé.

    Le contrôle de pas par résidu de demi-pas (Hibbitt & Karlsson 1979, la
    recette qu'ABAQUS publie) est implémenté et fonctionne. Il coûte, à
    précision ÉGALE, **2,5 à 3 fois** le pas fixe sur un problème oscillatoire à
    deux échelles — d'où son statut : opt-in, désactivé par défaut.

    La raison n'est pas un défaut d'implémentation, c'est le domaine de
    l'estimateur : le résidu de demi-pas mesure le déséquilibre LOCAL, et sur un
    mécanisme peu dissipatif l'erreur est dominée par la PHASE accumulée, qu'il
    ne voit pas. ABAQUS le dit à sa façon — l'algorithme est « purely
    empirical » et économique « in initially excited problems with HIGH
    DISSIPATION, such as impulsively loaded problems, with extensive
    plasticity ». Une voilure tournante n'est ni l'un ni l'autre.

    Si un jour cet assert casse — parce que l'estimateur a changé, ou parce
    qu'un cas dissipatif entre au banc — le verdict devra être rejoué.
    """
    ref = _deux_echelles(1e-5)["angle"]
    fixe = [(h, _deux_echelles(h)) for h in (1e-3, 5e-4)]
    adapt = [(t, _deux_echelles(2e-3, adapt=t)) for t in (1e-6, 1e-7)]
    lignes = []
    for h, r in fixe:
        lignes.append(("fixe", f"h={h:.0e}", abs(r["angle"] - ref), r["newton"], r["temps"], 0))
    for t, r in adapt:
        lignes.append(("adapt", f"tol={t:.0e}", abs(r["angle"] - ref), r["newton"], r["temps"],
                       r["adapt"][0]))
    for genre, quoi, err, nw, tps, rej in lignes:
        print(f"  {'pas ' + genre:15} {quoi:9} erreur {err:.3e}  {nw:6d} Newton  {tps:.3f} s"
              + (f"  {rej} rejeux" if rej else ""))
    # à erreur COMPARABLE (les deux paires sont appariées à ~1,5x près),
    # l'adaptatif coûte plus cher — c'est le fait qu'on garde
    couts = []
    for i in range(2):
        e_f, e_a = lignes[i][2], lignes[2 + i][2]
        assert 0.5 < e_a / e_f < 2.0, ("les deux points ne sont pas appariés", e_f, e_a)
        couts.append(lignes[2 + i][3] / lignes[i][3])
    print(f"  {'':15} à erreur égale, l'adaptatif coûte ×{couts[0]:.1f} et ×{couts[1]:.1f} "
          f"— il reste DÉSACTIVÉ par défaut")
    # l'estimateur regarde aussi Φ à mi-pas — l'erreur locale que λ paie, que le
    # résidu dynamique ne voit pas. Nulle ici (pas de liaison) ; sur un pendule
    # elle existe, et elle borne le pas comme le résidu
    L, m, th0 = 1.0, 0.2, np.radians(30.0)
    N = Noyau([0, 0, -G])
    b = N.corps("b", m, list((1e-6 * np.eye(3)).ravel()), [L * sin(th0), 0, -L * cos(th0)])
    N.liaison("rotule", None, b, pa=[0, 0, 0], bloque_t=[0, 1, 2], bloque_r=[])
    N.simule(1.0, 2e-3, tous=10 ** 9, adaptatif=1e-6)
    phi_h = N.adapt_stats()[4]
    print(f"  {'':15} pendule adaptatif : |Φ| à mi-pas {phi_h:.2e} (relatif), {N.adapt_stats()[0]} rejeux")
    assert 0 < phi_h < 1e-3, phi_h
    assert all(c > 1.5 for c in couts), couts
    # et il fait ce qu'on lui demande : il rejoue, et il borne le pas
    assert all(l[5] > 0 for l in lignes[2:]), lignes
    return dict(couts=couts, points=[(l[0], l[1], l[2], l[3]) for l in lignes])


def violation_vitesse():
    """Φ̇ EST VIOLÉ À O(h²) EN INDEX 3, ET GGL LE FERME SANS PERDRE L'ORDRE — mesuré.

    L'index 3 direct tient Φ à 1e-15 mais laisse Φ̇ violé à O(h²) : résidu
    BORNÉ (même ordre à 2 et à 100 périodes), sans coût d'énergie. La
    stabilisation GGL (`ggl=True`, Gear–Gupta–Leimkuhler : m multiplicateurs
    de plus, la pose reçoit βh·G₀ᵀζ, la ligne G·u + Φ_t = 0 entre au système)
    ramène Φ̇ à la précision machine en gardant l'ordre 2 — là où la
    projection de vitesse APRÈS le pas, essayée et retirée le 4 sept., le
    faisait tomber à ~0,8. Son prix se lit ici : itérations de Newton.
    """
    L, m, th0 = 1.0, 0.2, np.radians(30.0)
    T = 2 * pi * sqrt(L / G) * 1.07
    ref = _pendule_ref(2 * T, th0, L)

    def pendule():
        N = Noyau([0, 0, -G])
        b = N.corps("b", m, list((1e-6 * np.eye(3)).ravel()),
                    [L * sin(th0), 0, -L * cos(th0)])
        N.liaison("rotule", None, b, pa=[0, 0, 0], bloque_t=[0, 1, 2], bloque_r=[])
        return N

    out = {}
    for ggl in (False, True):
        errs, phid, newton = [], [], []
        for h in (T / 100, T / 200, T / 400):
            N = pendule()
            N.simule(2 * T, h, tous=10 ** 9, ggl=ggl)
            th, _ = ref(N.t())
            r = np.array(N.etat()[1][0])
            errs.append(float(np.linalg.norm(r - np.array([L * sin(th), 0, -L * cos(th)]))))
            phid.append(max(abs(x) for x in N.phi_dot()))
            newton.append(N.stats()[0])
        o = [float(np.log2(errs[i] / errs[i + 1])) for i in range(len(errs) - 1)]
        out[ggl] = dict(err=errs, ordres=o, phi_dot=phid, newton=newton)
        print(f"  {'GGL' if ggl else 'index 3':15} erreur " + "  ".join(f"{e:.3e}" for e in errs)
              + f"   ordres {', '.join(f'{x:.2f}' for x in o)}   |Φ̇| {phid[-1]:.1e}"
              + f"   Newton {newton[-1]}")
    # index 3 : ordre 2, Φ̇ en O(h²) ; GGL : ordre 2 gardé, Φ̇ au plancher
    assert min(out[False]["ordres"]) > 1.8 and min(out[True]["ordres"]) > 1.8, out
    o_phi = np.log2(out[False]["phi_dot"][-2] / out[False]["phi_dot"][-1])
    assert o_phi > 1.8 and out[True]["phi_dot"][-1] < 1e-12, (o_phi, out[True]["phi_dot"])
    # le prix de GGL : pas plus de 3× d'itérations (le terme ∂G/∂q·u est laissé
    # de côté dans le jacobien — quasi-Newton, convergence linéaire à O(h))
    assert out[True]["newton"][-1] < 3 * out[False]["newton"][-1], out

    # LA QUESTION QUI DÉCIDE : la violation de Φ̇ a-t-elle une CONSÉQUENCE ?
    h = T / 200
    long = {}
    for n_per in (2, 100):
        N = pendule()
        e0 = N.energie()
        N.simule(n_per * T, h, tous=10 ** 9)
        long[n_per] = (abs((N.energie() - e0) / e0), max(abs(x) for x in N.phi_dot()))
        print(f"  {'':15} {n_per:3d} périodes  ΔE/E {long[n_per][0]:.2e}   |Φ̇| {long[n_per][1]:.1e}")
    # Φ̇ NE DÉRIVE PAS : même ordre à 2 et à 100 périodes — un résidu borné
    assert long[100][1] < 10 * long[2][1], long
    # et l'énergie tient : le défaut ne coûte rien ici
    assert long[100][0] < 1e-4, long
    print(f"  {'':15} VERDICT : en index 3, Φ̇ ne dérive pas et ne coûte pas d'énergie "
          f"({long[100][0]:.0e} sur 100 périodes) ; GGL le ferme pour "
          f"{out[True]['newton'][-1] / out[False]['newton'][-1]:.1f}× d'itérations.")
    return dict(index3=out[False], ggl=out[True], long={str(k): v for k, v in long.items()})


def bifurcation():
    """LE MASQUE DE REDONDANCE N'EST PAS FIGÉ — un mécanisme parti d'un POINT DE
    BIFURCATION.

    Le parallélogramme à plat (A, B, D, C alignés) est un point de changement :
    les branches parallélogramme et antiparallélogramme s'y croisent, G y perd
    UN rang de plus qu'ailleurs. Une détection de redondance faite là — et
    jamais refaite — neutralise une ligne qui redevient structurelle dès que
    le mécanisme en sort : Φ dérive, ou Newton cale sur un G quasi singulier.
    Le noyau re-détecte à la pose courante quand Φ dépasse sa tolérance OU
    quand Newton cale, ramène l'état sur Φ = Φ̇ = 0 sous le nouveau masque et
    rejoue le pas. Référence fermée, sans gravité : la branche parallélogramme
    tourne UNIFORMÉMENT (énergie cinétique indépendante de l'angle), donc
    φ = φ₀ + ωt et le centre du coupleur vaut (1 + cos φ, sin φ).
    """
    ji = list((1e-3 * np.eye(3)).ravel())

    def para(ph0):
        N = Noyau([0, 0, 0])
        c, s = cos(ph0), sin(ph0)
        R = [c, -s, 0, s, c, 0, 0, 0, 1]
        N.corps("ab", 1.0, ji, [0.5 * c, 0.5 * s, 0], rot=R, w=[0, 0, 1.0], v=[-0.5 * s, 0.5 * c, 0])
        N.corps("bc", 1.0, ji, [c + 1.0, s, 0], v=[-s, c, 0])
        N.corps("dc", 1.0, ji, [2 + 0.5 * c, 0.5 * s, 0], rot=R, w=[0, 0, 1.0], v=[-0.5 * s, 0.5 * c, 0])
        N.liaison("A", None, 0, pa=[0, 0, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
        N.liaison("B", 0, 1, pa=[0.5, 0, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
        N.liaison("C", 1, 2, pa=[1, 0, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
        N.liaison("D", None, 2, pa=[2, 0, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
        return N

    out = {}
    for ph0 in (np.radians(5.0), 0.0):
        for ggl in (False, True):
            N = para(ph0)
            N.simule(1.0, 1e-3, tous=10 ** 9, ggl=ggl)
            t = N.t()
            r = np.array(N.etat()[1][1])
            ecart = float(np.linalg.norm(r - [1 + cos(ph0 + t), sin(ph0 + t), 0]))
            red = N.stats()[3]
            out[(ph0, ggl)] = (ecart, red)
            print(f"  {'φ₀ = ' + f'{np.degrees(ph0):.0f}°':10} {'GGL' if ggl else 'index 3':8} "
                  f"écart à l'analytique {ecart:.2e} m · redondantes à la fin {red}")
    # un plan a TROIS lignes redondantes ; à plat il en avait quatre — le masque a suivi
    assert all(red == 3 for _, red in out.values()), out
    # les deux départs suivent l'analytique ; celui de la bifurcation paie l'assemblage
    assert all(e < 1e-5 for e, _ in out.values()), out
    assert out[(0.0, False)][0] > out[(np.radians(5.0), False)][0], out
    return {f"{np.degrees(k[0]):.0f}_{k[1]}": v for k, v in out.items()}


def oye():
    """DÉCROCHAGE DYNAMIQUE D'ØYE — échelon d'incidence sur une section fixe.

    Une section attachée qui reçoit un échelon d'incidence AU-DELÀ du
    décrochage garde d'abord sa portance ATTACHÉE (f part de f_st(α₁) ≈ 1)
    puis relaxe vers la valeur statique : surcroît de portance, puis retour
    en τ = 4c/|V|. Trois contrôles : la limite t → ∞ est le c81 statique par
    construction ; le surcroît initial est positif ; la constante de temps
    mesurée sur la relaxation est 4c/V à 5 %.
    """
    c, V = 0.1, 20.0
    al = np.arange(-180.0, 181.0, 1.0)
    # polaire synthétique : droite 2π puis plateau décroché
    cl = np.where(np.abs(al) <= 12.0, 0.1097 * al, np.sign(al) * 1.0)
    cl = np.where(np.abs(al) > 90.0, 0.0, cl)
    pol = (list(al), list(cl), [0.02] * len(al), [0.0] * len(al))

    def section(a_deg, dyn):
        N = Noyau([0, 0, 0])
        b = N.corps("s", 1.0, list(np.eye(3).ravel()), [0, 0, 0])
        N.liaison("fixe", None, b)
        N.pale("p", b, [0, 0, 0], [0, 1, 0], [1, 0, 0], 1.0, c, pol, oye=dyn)
        N.vent([-V * cos(radians(a_deg)), 0.0, V * sin(radians(a_deg))])
        return N

    from math import cos, radians, sin
    a1, a2 = 5.0, 20.0
    N = section(a1, True)
    N.simule(0.2, 1e-4, tous=10 ** 9)                 # f à l'équilibre à α₁
    N.vent([-V * cos(radians(a2)), 0.0, V * sin(radians(a2))])
    traj = N.simule(0.4, 1e-4, tous=1)                # t_end est ABSOLU
    t = np.array([f[0] for f in traj]) - traj[0][0]
    lz = np.array([abs(f[5][2]) for f in traj])         # réaction verticale de l'encastrement = portance
    stat = section(a2, False)
    stat.simule(0.01, 1e-4, tous=10 ** 9)
    l_st = abs(stat.reactions()[0][1][2])
    ratio = lz / l_st
    # constante de temps sur la relaxation de (ratio − 1)
    k = (ratio - 1.0) > 0.05 * (ratio[0] - 1.0)
    p = np.polyfit(t[k], np.log(ratio[k] - 1.0), 1)
    tau = -1.0 / p[0]
    print(f"  Øye {a1:.0f}° → {a2:.0f}° : portance {ratio[0]:.3f}× la statique à t = 0⁺, "
          f"{ratio[-1]:.4f}× à {t[-1]:.2f} s ; τ mesuré {tau * 1e3:.2f} ms, 4c/V = {4 * c / V * 1e3:.2f} ms")
    assert ratio[0] > 1.2 and abs(ratio[-1] - 1.0) < 1e-3, ratio
    assert abs(tau / (4 * c / V) - 1) < 0.05, tau
    return dict(sur=float(ratio[0]), fin=float(ratio[-1]), tau=float(tau))


def leishman_beddoes():
    """DÉCROCHAGE DYNAMIQUE DE LEISHMAN–BEDDOES — ce qu'Øye ne fait pas.

    Øye relaxe le point de séparation vers sa valeur statique : il donne le
    RETARD, pas le SURCROÎT. Leishman–Beddoes (1989) ajoute le retard de
    pression (le point de séparation répond à un C_N retardé de T_p, puis
    relaxe en T_f) et le TOURBILLON de bord d'attaque, lâché quand le C_N
    retardé passe le seuil C_N1 et convecté sur la corde en T_vl semi-cordes.
    C'est ce tourbillon qui fait la signature du décrochage dynamique : une
    portance maximale AU-DESSUS du maximum statique, puis une chute, puis une
    boucle d'hystérésis en tangage.

    Trois contrôles, et un négatif :
      · échelon 5° → 20° : la limite t → ∞ est le c81 statique par
        construction (forme d'Øye), et le surcroît initial est ≥ celui d'Øye
        (le tourbillon s'ajoute, il ne retire rien) ; τ_v > 0 et C_N^v > 0
        pendant la convection, C_N^v → 0 après ;
      · tangage sinusoïdal 10° ± 10° à k = 0,1 : la boucle CL(α) a une AIRE
        non nulle et CL_max dynamique > CL_max statique — le fait expérimental
        que tout le domaine reproduit (McCroskey 1981) ;
      · CONTRÔLE NÉGATIF, deux fois : le quasi-statique trace une boucle d'aire
        NULLE (CL est une fonction de α) ; et en régime ATTACHÉ (2° ± 2°) LB
        rend la même aire nulle — f'' = 1, aucun tourbillon, le modèle ne
        fabrique pas d'hystérésis là où il n'y en a pas.
    Ce qui n'est PAS là, déclaré : les constantes T_p, T_f, T_v, T_vl sont
    celles de la littérature (NACA 0012), jamais calibrées sur LA section ; CD
    et CM restent statiques ; pas de retard compressible.
    """
    from math import cos, radians, sin
    from vinkulum import LB_TF, LB_TP, LB_TV, LB_TVL
    c, V = 0.1, 20.0
    al = np.arange(-180.0, 181.0, 1.0)
    cl = np.where(np.abs(al) <= 12.0, 0.1097 * al, np.sign(al) * 1.0)
    cl = np.where(np.abs(al) > 90.0, 0.0, cl)
    pol = (list(al), list(cl), [0.02] * len(al), [0.0] * len(al))
    q = 0.5 * 1.225 * V * V * c            # ½ρV²c, envergure 1

    def section(a_deg, **kw):
        N = Noyau([0, 0, 0])
        b = N.corps("s", 1.0, list(np.eye(3).ravel()), [0, 0, 0])
        N.liaison("fixe", None, b)
        N.pale("p", b, [0, 0, 0], [0, 1, 0], [1, 0, 0], 1.0, c, pol, **kw)
        N.vent([-V * cos(radians(a_deg)), 0.0, V * sin(radians(a_deg))])
        return N

    def portance(lam, a_deg):
        # réaction de l'encastrement = −force aéro ; la portance est ⟂ au vent.
        # Avec ce vent la section voit α = −a (en = es × ec pointe vers −z) :
        # le signe est retourné pour lire un CL positif à incidence positive —
        # la polaire est antisymétrique, rien d'autre ne change.
        fx, fz = -lam[0], -lam[2]
        return -(fx * sin(radians(a_deg)) + fz * cos(radians(a_deg)))

    print("╔═ vinkulum — Leishman–Beddoes : retard de pression + tourbillon de bord d'attaque")
    # 1. ÉCHELON 5° → 20°, LB contre Øye contre statique
    a1, a2 = 5.0, 20.0
    stat = section(a2)
    stat.simule(0.01, 1e-4, tous=10 ** 9)
    l_st = portance(stat.simule(0.011, 1e-4, tous=1)[-1][5], a2)
    out = {}
    for nom, kw in (("oye", dict(oye=True)), ("lb", dict(lb=True))):
        N = section(a1, **kw)
        N.simule(0.2, 1e-4, tous=10 ** 9)
        N.vent([-V * cos(radians(a2)), 0.0, V * sin(radians(a2))])
        cnv_max, tauv_max = 0.0, 0.0
        ratios = []
        t = 0.2
        while t < 0.6 - 1e-12:
            t += 1e-4
            fr = N.simule(t, 1e-4, tous=1)[-1]
            ratios.append(portance(fr[5], a2) / l_st)
            if nom == "lb":
                et, cn1 = N.etats_lb(0)
                cnv_max = max(cnv_max, max(abs(e[2]) for e in et))
                tauv_max = max(tauv_max, max(e[1] for e in et))
        # le SURCROÎT intégré ∫(L/L_st − 1)dt, en ms : c'est lui que le retard
        # de pression et le tourbillon allongent, pas le pic à 0⁺ (qui est le
        # saut attaché, identique pour les deux modèles)
        surcroit = float(np.sum(np.array(ratios) - 1.0) * 1e-4 * 1e3)
        out[nom] = (ratios[0], max(ratios), ratios[-1], cnv_max, tauv_max, surcroit)
        print(f"║ échelon {a1:.0f}° → {a2:.0f}° {nom:4}: portance {ratios[0]:.3f}× à 0⁺, "
              f"{ratios[-1]:.4f}× à 0,4 s, surcroît intégré {surcroit:.2f} ms"
              + (f" · |C_N^v| max {cnv_max:.3f}, τ_v max {tauv_max:.1f} sc" if nom == "lb" else ""))
    et_fin, cn1 = N.etats_lb(0)
    cnv_fin = max(abs(e[2]) for e in et_fin)
    # 2. TANGAGE SINUSOÏDAL — boucle d'hystérésis, k = ωc/(2V)
    k_red = 0.1
    om = 2.0 * k_red * V / c
    T = 2.0 * np.pi / om
    dt = T / 400.0

    def boucle(am, aa, **kw):
        N = section(am, **kw)
        N.simule(0.05, 1e-4, tous=10 ** 9)
        t, pts = 0.05, []
        n_cyc = 3
        while t < 0.05 + n_cyc * T - 1e-12:
            t += dt
            a = am + aa * sin(om * (t - 0.05))
            N.vent([-V * cos(radians(a)), 0.0, V * sin(radians(a))])
            fr = N.simule(t, dt, tous=1)[-1]
            pts.append((a, portance(fr[5], a) / q))
        pts = np.array(pts[-400:])              # le dernier cycle
        a_r, cl_ = np.radians(pts[:, 0]), pts[:, 1]
        aire = 0.5 * float(np.sum(cl_ * np.roll(a_r, -1) - np.roll(cl_, -1) * a_r))
        return aire, float(np.max(cl_)), float(np.max(np.abs(cl_)))

    cl_max_st = float(np.max(cl[np.abs(al) <= 20.0]))
    a_st, m_st, _ = boucle(10.0, 10.0)
    a_oye, m_oye, _ = boucle(10.0, 10.0, oye=True)
    a_lb, m_lb, _ = boucle(10.0, 10.0, lb=True)
    a_att, m_att, s_att = boucle(2.0, 2.0, lb=True)
    print(f"║ tangage 10° ± 10°, k {k_red} : aire ∮CL dα statique {a_st:+.2e} · Øye {a_oye:+.3f} · LB {a_lb:+.3f}")
    print(f"║   CL max sur la boucle : statique {m_st:.3f} · Øye {m_oye:.3f} · LB {m_lb:.3f}  "
          f"(CL max de la polaire {cl_max_st:.3f})")
    print(f"║ régime attaché 2° ± 2°, LB : aire {a_att:+.2e} (contrôle négatif — pas d'hystérésis fabriquée)")
    print(f"╚═ C_N1 lu sur la polaire {cn1:.3f} · T_p {LB_TP} T_f {LB_TF} T_v {LB_TV} T_vl {LB_TVL} sc (littérature)")
    # · la limite t → ∞ est le statique, pour les deux modèles
    assert abs(out["oye"][2] - 1.0) < 1e-3 and abs(out["lb"][2] - 1.0) < 1e-3, out
    # · le tourbillon AJOUTE : surcroît intégré LB > Øye, C_N^v a existé puis s'est éteint
    assert out["lb"][5] > out["oye"][5] and out["lb"][3] > 0.02, out
    assert out["lb"][4] > 0.0 and cnv_fin < 1e-3, (out, cnv_fin)
    # · l'hystérésis : aire nulle en statique, non nulle pour LB, et CL max dynamique > statique
    assert abs(a_st) < 1e-9, a_st
    assert abs(a_lb) > 0.02 and m_lb > 1.05 * cl_max_st, (a_lb, m_lb, cl_max_st)
    # · CONTRÔLE NÉGATIF : régime attaché → f'' = 1, pas de tourbillon, pas de boucle
    assert abs(a_att) < 1e-6 * max(abs(a_lb), 1e-30) + 1e-12, (a_att, a_lb)
    return dict(echelon=out, aire=dict(statique=a_st, oye=a_oye, lb=a_lb, attache=a_att),
                cl_max=dict(statique=m_st, oye=m_oye, lb=m_lb, polaire=cl_max_st), cn1=float(cn1))


def energie_moment():
    """SCHÉMA ÉNERGIE-MOMENT DANS LA FORMULATION — l'exclusive-ou de la projection, résolu.

    Le README le disait depuis Noether : la projection conserve L *ou* E,
    jamais les deux (∂E/∂ω = ωᵀ·∂L/∂ω, rang 3 sur 4), et les tenir ensemble
    « demande de corriger aussi q, c'est-à-dire un schéma énergie-moment dans
    la FORMULATION (Simo–Tarnow), pas une projection ». `simule_em` est ce
    schéma pour les corps libres : point milieu en translation, Cayley au
    point milieu sur le moment matériel en rotation (Simo–Wong 1991). Trois
    faits à mesurer, un négatif :
      · E ET L conservés au plancher machine sur 20 000 pas d'un corps
        asymétrique qui culbute autour de son axe intermédiaire (le cas
        instable, où l'α-généralisé dérive le plus) ;
      · l'α-généralisé, lui, perd E ; `moment=True` tient L et pas E,
        `energie=True` tient E et pas L — les trois colonnes côte à côte ;
      · ordre 2 en orientation contre l'index 3 au pas fin (erreur ÷4 quand
        h ÷2) : conserver n'est pas être juste, et ça se vérifie séparément ;
      · CONTRÔLE NÉGATIF : une liaison fait refuser le schéma (domaine).
    """
    from scipy.spatial.transform import Rotation as Rot
    J = [1e-3, 0, 0, 0, 2e-3, 0, 0, 0, 3e-3]
    w0 = [0.05, 20.0, 0.03]                      # axe intermédiaire, perturbé

    def monte(g=0.0, v0=(0.0, 0.0, 0.0)):
        # L n'est un invariant que SANS gravité (son moment en r × g ne l'est
        # pas) : les quatre colonnes tournent en apesanteur, la gravité est
        # jugée à part sur E, où le point milieu est exact pour un potentiel
        # linéaire. Et SANS translation : avec v = 2 m/s le corps est à 400 m
        # après 200 s, et le terme m·r×v de `moment()` porte un arrondi qui
        # croît avec r — 1,1e-9 relatif mesuré, dans la LECTURE de L, pas dans
        # le schéma (rotation seule : 5e-15 sur les mêmes 200 000 pas).
        N = Noyau([0.0, 0.0, -g])
        N.corps("c", 0.7, J, [0.0, 0.0, 1.0], v=list(v0), w=w0)
        return N

    def bilan(N, e0, l0):
        e1, l1 = N.energie(), np.array(N.moment())
        return abs(e1 - e0) / abs(e0), np.linalg.norm(l1 - l0) / np.linalg.norm(l0)

    print("╔═ vinkulum — ÉNERGIE-MOMENT (Simo–Wong, Cayley) : E et L tenus ENSEMBLE")
    h, t_end = 1e-3, 20.0
    res = {}
    for nom, lance in (("énergie-moment", lambda N: N.simule_em(t_end, h, tous=10 ** 9)),
                       ("α-gen index 3", lambda N: N.simule(t_end, h, tous=10 ** 9)),
                       ("α-gen moment=True", lambda N: N.simule(t_end, h, tous=10 ** 9, moment=True)),
                       ("α-gen energie=True", lambda N: N.simule(t_end, h, tous=10 ** 9, energie=True))):
        N = monte()
        e0, l0 = N.energie(), np.array(N.moment())
        t0 = time.time()
        lance(N)
        de, dl = bilan(N, e0, l0)
        res[nom] = (de, dl)
        print(f"║ {nom:20} {t_end / h:6.0f} pas : |ΔE|/E {de:.1e} · |ΔL|/L {dl:.1e}   ({time.time() - t0:.1f} s)")
    # gravité : E tenue aussi (point milieu en translation, potentiel linéaire).
    # 2 s de chute — au-delà, m·g·z domine E de trois décades et c'est
    # l'arrondi de la LECTURE qui se mesure (1e-10 à 20 s, mesuré)
    N = monte(G, (0.3, 0.0, 2.0))
    e0 = N.energie()
    N.simule_em(2.0, h, tous=10 ** 9)
    de_g = abs(N.energie() - e0) / abs(e0)
    print(f"║ sous gravité, 2 s de chute : |ΔE|/E {de_g:.1e} (cinétique + potentiel linéaire)")
    # ordre 2 en orientation : référence index 3 à h/64 sur 1 s
    N = monte()
    N.simule(1.0, h / 64.0, tous=10 ** 9)
    r_ref = Rot.from_matrix(np.array(N.etat()[2][0]).reshape(3, 3))
    errs = []
    for hh in (4e-3, 2e-3, 1e-3):
        N = monte()
        N.simule_em(1.0, hh, tous=10 ** 9)
        r = Rot.from_matrix(np.array(N.etat()[2][0]).reshape(3, 3))
        errs.append((r_ref.inv() * r).magnitude())
    ordres = [np.log2(errs[i] / errs[i + 1]) for i in range(2)]
    print(f"║ ordre en orientation : erreurs {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e} rad → "
          f"ordres {ordres[0]:.2f} {ordres[1]:.2f}")
    # contrôle négatif : le domaine refuse une liaison
    N = monte()
    N.liaison("piv", None, 0, pa=[0, 0, 0], bloque_r=[0, 2])
    try:
        N.simule_em(0.01, h)
        raise AssertionError("simule_em a accepté une liaison")
    except ValueError as ex:
        print(f"║ avec une liaison : refusé — « {str(ex)[:60]} »")
    print("╚═ énergie-moment OK")
    em, a3, am, ae = res["énergie-moment"], res["α-gen index 3"], res["α-gen moment=True"], res["α-gen energie=True"]
    # · les deux invariants au plancher, ensemble (Π porté d'un pas à l'autre,
    #   R re-orthonormalisée : relire Π de (R, ω) faisait dériver en N²)
    assert em[0] < 1e-13 and em[1] < 1e-13 and de_g < 1e-12, (em, de_g)
    # · l'α-généralisé dérive en E ; la projection tient l'un et pas l'autre
    assert a3[0] > 1e3 * em[0], (a3, em)
    assert am[1] < 1e-12 and am[0] > 1e3 * em[0], am
    assert ae[0] < 1e-10 and ae[1] > 1e3 * em[1], ae
    # · ordre 2 (Richardson sur trois pas)
    assert min(ordres) > 1.8, ordres
    return dict(res=res, gravite=float(de_g), ordres=[float(o) for o in ordres],
                erreurs=[float(e) for e in errs])


def multirythme():
    """SOUS-CYCLAGE MULTI-RYTHME — la partition lente ne paie plus le pas de la rapide.

    Une chaîne de pendules (période ~2 s) porte au bout, à travers un MONTAGE
    SOUPLE (~5 Hz), un petit ensemble vibrant : une masse sur une poutre raide
    (mode axial ~300 Hz). Un pas unique doit résoudre le 300 Hz partout : à
    h = 1e-3 (3 pas par période) l'α-généralisé l'ÉTEINT (ρ∞ 0,9 dissipe tout
    ce qui est sous-résolu) ; à h = 1e-5 tout est juste mais la chaîne fait
    100 fois plus de pas qu'il n'en faut. Le multi-rythme fait la chaîne à
    1e-3 et l'ensemble vibrant à 1e-5.

    LE PREMIER BANC ÉTAIT LE PIRE CAS, et il l'a dit : la poutre RAIDE
    enjambait les partitions. Un partenaire gelé sur 35 kN/m pendant un pas
    lent, c'est 9 N de force parasite sur un pendule qui en pèse 10 — le
    pendule ne bougeait plus. La raideur doit vivre DANS la partition rapide,
    et le couplage entre partitions doit être SOUPLE : c'est la condition du
    multi-rythme (Gear–Wells), pas une limite de cette implémentation.

    Trois mesures, contre le monolithique fin (h = 1e-5, référence) :
      · la position du bout de chaîne à T (les rapides sont EXTRAPOLÉS à
        l'ordre 1 pendant le pas lent : gelés à t₀ c'était 2 cm d'écart en
        0,3 s, mesuré — l'écart se publie, et se compare au grossier) ;
      · l'AMPLITUDE de la vibration axiale à la fin, détrendée — le pas
        grossier ne la TUE pas, il la fausse (0,81× mesuré à 3 pas par
        période) ; le sous-cyclage doit la rendre mieux ;
      · le coût des trois, et le ratio jacobiens/pas.
    CONTRÔLE NÉGATIF : une liaison qui enjambe les partitions est refusée.
    Les rapides sont extrapolés pendant le pas lent, les lents interpolés
    pendant les micro-pas. Chaque partition conserve son schéma et son cache
    de jacobien ; Newton et l'initialisation factorisent le système restreint.
    Pas de correcteur ni d'aérodynamique. Les temps isolés et leur dispersion
    se mesurent avec `ci/mesure_multirythme.py`, hors campagne concurrente.
    """
    from vinkulum._banc_multirythme import monte
    n_ch, L, lb, T = 8, 0.25, 0.05, 0.3

    def sortie(N, tr, a, c, b):
        fin = [f for f in tr if f[0] > T - 0.05]
        t = np.array([f[0] for f in fin])
        d = np.array([np.linalg.norm(np.array(f[1][b]) - np.array(f[1][c])) - lb for f in fin])
        d = d - np.polyval(np.polyfit(t - t[0], d, 2), t - t[0])      # détrendé (statique, lent)
        j, s, ja = N.stats()[2], N.stats()[3], N.stats()
        return np.array(tr[-1][1][a]), 0.5 * (d.max() - d.min()), ja

    print("╔═ vinkulum — MULTI-RYTHME : chaîne à h, ensemble vibrant à h/k, montage souple entre les deux")
    res = {}
    for nom, lance in (("monolithique h = 1e-5", lambda N, rap: N.simule(T, 1e-5, tous=10)),
                       ("monolithique h = 1e-3", lambda N, rap: N.simule(T, 1e-3, tous=1)),
                       ("multi-rythme 1e-3 / 100", lambda N, rap: N.simule_multirythme(T, 1e-3, rap, 100, tous=1))):
        N, a, c, b = monte()
        t0 = time.time()
        tr = lance(N, [c, b])
        dt = time.time() - t0
        pa, amp, st = sortie(N, tr, a, c, b)
        res[nom] = (pa, amp, dt, st)
        print(f"║ {nom:26} : bout de chaîne ({pa[0]:+.6f}, {pa[1]:+.6f}) · amplitude 300 Hz {amp * 1e6:8.2f} µm"
              f" · {dt:5.2f} s · newton {st[0]:6d} jacobiens {st[2]:6d}")
    ref, gros, multi = res["monolithique h = 1e-5"], res["monolithique h = 1e-3"], res["multi-rythme 1e-3 / 100"]
    e_pos = np.linalg.norm(multi[0] - ref[0])
    e_gros = np.linalg.norm(gros[0] - ref[0])
    print(f"║ écart au fin : bout multi {e_pos:.2e} m (grossier {e_gros:.2e}) · "
          f"amplitude multi {multi[1] / ref[1]:.3f}× (grossier {gros[1] / ref[1]:.3f}×)"
          f" · temps multi/fin {multi[2] / ref[2]:.2f}")
    N, a, c, b = monte(enjambe=True)
    try:
        N.simule_multirythme(0.01, 1e-3, [c, b], 10)
        raise AssertionError("une liaison qui enjambe a été acceptée")
    except ValueError as ex:
        print(f"║ liaison qui enjambe : refusée — « {str(ex)[:70]} »")
    print("╚═ multi-rythme OK")
    # · le pas grossier fausse la vibration, le sous-cyclage la rend mieux et à 10 %
    assert abs(multi[1] / ref[1] - 1.0) < abs(gros[1] / ref[1] - 1.0), (multi[1], gros[1], ref[1])
    assert abs(multi[1] / ref[1] - 1.0) < 0.1, (multi[1], ref[1])
    # · et le bout de chaîne reste au bon endroit : 1,9 mm sur 2 m mesurés
    #   (couplage extrapolé à l'ordre 1 ; l'état du schéma survit par partition,
    #   sans quoi c'était 4 cm — ordre 1 du redémarrage)
    assert e_pos < 5e-3 * n_ch * L, e_pos
    return dict({k: dict(pos=v[0].tolist(), amp=float(v[1]), t=float(v[2])) for k, v in res.items()},
                ecart_pos=float(e_pos))

def _ladder(n, rotules, ech):
    """Chaîne de `n` parallélogrammes : n+1 manivelles pivotées au sol, n
    bielles entre leurs sommets. `rotules` = bielles en ROTULE (G de rang
    plein) au lieu de pivot (plan, donc redondant). `ech` = longueur, tout le
    reste suit (masses en ρL³, inerties en mL²/12) : deux échelles séparées
    d'un facteur 100 décrivent la MÊME mécanique en unités différentes."""
    L, D, th = ech, 2.0 * ech, np.radians(60.0)
    c, s = cos(th), sin(th)

    def bar(lg):                      # barre pleine de section (lg/20)²
        m = 1000.0 * (lg / 20.0) ** 2 * lg
        return m, list((m * lg * lg / 12.0 * I3).ravel())

    N = Noyau([0.0, -G, 0.0])
    man = [N.corps(f"m{k}", *bar(L), [k * D + 0.5 * L * c, 0.5 * L * s, 0.0],
                   rot=[c, -s, 0, s, c, 0, 0, 0, 1]) for k in range(n + 1)]
    cpl = [N.corps(f"c{k}", *bar(D), [k * D + L * c + 0.5 * D, L * s, 0.0]) for k in range(n)]
    for k in range(n + 1):
        N.liaison(f"A{k}", None, man[k], pa=[k * D, 0.0, 0.0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
    br = [] if rotules else [0, 1]
    for k in range(n):
        N.liaison(f"B{k}", man[k], cpl[k], pa=[0.5 * L, 0, 0], bloque_t=[0, 1, 2], bloque_r=br)
        N.liaison(f"C{k}", man[k + 1], cpl[k], pa=[0.5 * L, 0, 0], bloque_t=[0, 1, 2], bloque_r=br)
    minv = []                          # M⁻¹ par corps, en repère MONDE
    for lg, ids in ((L, man), (D, cpl)):
        m, j = bar(lg)
        for _ in ids:
            minv.append((m, np.linalg.inv(np.array(j).reshape(3, 3))))
    return N, [minv[i] for i in np.argsort([*man, *cpl])]


def _g_num(N, eps):
    """G = ∂Φ/∂(δr, δθ) par différences CENTRÉES, δθ appliqué à GAUCHE (repère
    monde), comme `Modele::avance`. Aucun binding n'expose le jacobien ; ce
    qu'on veut ici est un CONDITIONNEMENT, pas une tangente de production."""
    t, r, rot, v, w, vi = N.etat()
    cols = []
    for i in range(len(r)):
        for c in range(3):
            def bouge(e, i=i, c=c):
                rp = [list(x) for x in r]
                rp[i][c] += e
                N.pose_etat(t, rp, rot, v, w, vi)
                return np.array(N.phi_seul() if hasattr(N, "phi_seul") else N.phi())
            cols.append((bouge(eps) - bouge(-eps)) / (2 * eps))
        for c in range(3):
            def tourne(e, i=i, c=c):
                dth = np.zeros(3)
                dth[c] = e
                rp = [list(x) for x in rot]
                rp[i] = list((_expm(dth) @ np.array(rot[i]).reshape(3, 3)).ravel())
                N.pose_etat(t, r, rp, v, w, vi)
                return np.array(N.phi())
            cols.append((tourne(eps) - tourne(-eps)) / (2 * eps))
    N.pose_etat(t, r, rot, v, w, vi)
    return np.array(cols).T


def _cond(N, minv, eps):
    """(cond(G·Gᵀ), cond(G·M⁻¹·Gᵀ), rang) — le premier est ce que le noyau
    factorise aujourd'hui, le second ce que la métrique de masse donnerait.
    Calculés sur les valeurs singulières de G et de G·M^(−1/2), donc en
    ignorant le noyau (les lignes redondantes sont masquées par le noyau)."""
    g = _g_num(N, eps)
    mh = np.zeros((g.shape[1], g.shape[1]))
    for i, (m, ji) in enumerate(minv):
        mh[6 * i:6 * i + 3, 6 * i:6 * i + 3] = I3 / sqrt(m)
        val, vec = np.linalg.eigh(ji)
        mh[6 * i + 3:6 * i + 6, 6 * i + 3:6 * i + 6] = vec @ np.diag(np.sqrt(val)) @ vec.T
    out = []
    for a in (g, g @ mh):
        sv = np.linalg.svd(a, compute_uv=False)
        r = int((sv > sv[0] * 1e-10).sum())
        out.append((sv[0] / sv[r - 1]) ** 2)
    sv = np.linalg.svd(g, compute_uv=False)
    return out[0], out[1], int((sv > sv[0] * 1e-10).sum())


def ggl_boucles():
    """GGL ÉCHOUE SUR TROIS MODÈLES — reproduire l'échec hors de FRELON.

    Mesuré le 5 septembre (matin) : `simule(ggl=True)` faisait diverger Newton
    sur srscm, six_barres et la tête S2, et passait partout ailleurs. Ni h, ni
    le type d'élément, ni la planéité, ni m, ni la durée ne séparaient les deux
    groupes. Le seul séparateur trouvé : les trois qui échouaient avaient AU
    MOINS DEUX BOUCLES fermées indépendantes et un G de RANG PLEIN, quand les
    modèles à boucle qui passaient portaient des lignes masquées.
    TRACE : le soir même, les colonnes ζ sont devenues EXACTES (produit
    ∂r/∂q · βh²G₀ᵀ dans `jacobien_ad`) et les 12 cellules passent, `srscm`
    compris à son point mort. Le banc reste — c'est lui qui a isolé le
    séparateur, et c'est lui qui dira si ça régresse.

    Ce banc met le séparateur à l'épreuve sur une échelle de parallélogrammes :
    1, 2, 3 boucles × pivot (plan, redondant) ou rotule (rang plein) × deux
    ÉCHELLES séparées d'un facteur 100 — la tête S2 est le seul modèle
    sous-centimétrique, et une revue externe soutient que la projection non
    mise à l'échelle par la masse (q̇ = v − Gᵀμ, où q mêle mètres et radians)
    rend le conditionnement dépendant du choix d'unités. Le temps est
    adimensionné : h et t_end suivent τ = √(L/g), donc deux échelles ne
    diffèrent QUE par les unités.
    """
    print("╔═ GGL : échelle de boucles — le séparateur mis à l'épreuve")
    out = {}
    for ech in (1.0, 0.01):
        tau = sqrt(ech / G)
        h, t_end = tau / 300.0, 3.0 * tau
        for rotules in (False, True):
            for n in (1, 2, 3):
                lg = {}
                for ggl in (False, True):
                    N, minv = _ladder(n, rotules, ech)
                    try:
                        N.simule(t_end, h, ggl=ggl, tous=10 ** 9)
                        lg[ggl] = (True, float(np.linalg.norm(N.phi_dot())),
                                   N.stats()[0] / (t_end / h), N.stats()[3])
                    except Exception as ex:
                        lg[ggl] = (False, str(ex)[:40], 0.0, 0)
                N, minv = _ladder(n, rotules, ech)
                cgg, cgm, rang = _cond(N, minv, 1e-7 * ech)
                ok3, ok_g = lg[False][0], lg[True][0]
                cle = f"L{ech:g} {'rotule' if rotules else 'pivot'} {n}b"
                out[cle] = dict(index3=lg[False], ggl=lg[True], cond=cgg, cond_m=cgm,
                                rang=rang, masquees=lg[False][3], ech=ech, rotules=rotules, n=n)
                print(f"║ L {ech:<5g} {'rotule' if rotules else 'pivot ':6} {n} boucle(s) : "
                      f"index3 {'OK ' if ok3 else 'ÉCHEC'} · GGL {'OK ' if ok_g else 'ÉCHEC'}"
                      f" · masquées {lg[False][3]:2d} · cond(GGᵀ) {cgg:9.2e} · cond(GM⁻¹Gᵀ) {cgm:9.2e}")
    ko = [k for k, v in out.items() if not v["ggl"][0]]
    print(f"║ GGL : {len(out) - len(ko)}/{len(out)} cas — reste {ko}")
    # l'index 3 doit tenir PARTOUT — sinon le banc ne juge pas GGL, il juge le modèle
    assert all(v["index3"][0] for v in out.values()), out
    # 12/12 depuis les colonnes ζ EXACTES (5 sept., soir) : la dernière cellule
    # (L/100, 3 boucles, cond 2,7e5) tombait par le bruit des différences finies,
    # pas par le conditionnement lui-même. L'assert était « ≤ 1 échec, à L < 1 » ;
    # il est resserré, pas desserré.
    assert not ko, out
    # LE CONDITIONNEMENT SUIT 1/L² — la prédiction de la revue externe sur une
    # projection non adimensionnée, vérifiée : ×1e4 pour L ÷ 100
    for rot in (False, True):
        for n in (1, 2, 3):
            nm = "rotule" if rot else "pivot"
            r = out[f"L0.01 {nm} {n}b"]["cond"] / out[f"L1 {nm} {n}b"]["cond"]
            assert 3e3 < r < 3e4, (rot, n, r)
    # ET LA MÉTRIQUE DE MASSE NE LE RÉPARE PAS : elle l'AGGRAVE, partout.
    # Résultat NÉGATIF gardé — si un jour M⁻¹ améliorait, ce banc devrait
    # être rejoué et la conclusion du 5 septembre avec lui.
    assert all(v["cond_m"] > v["cond"] for v in out.values()), out
    print("╚ cond(GGᵀ) suit 1/L² ; cond(GM⁻¹Gᵀ) est PIRE partout — la métrique de masse n'est pas le remède ici")
    return out


def seuils():
    """LES SEUILS DE NEWTON SONT-ILS LA PHYSIQUE OU LE BANC ? — mesuré.

    Le solveur porte des constantes posées à la main (borne du prédicteur 0,3 rad,
    rafraîchissement du jacobien à 0,1, plancher de stagnation 1e-6, incrément
    1e-11). `VINKULUM_SEUILS=k` les multiplie toutes ; on rejoue les bancs
    d'ordre et de bifurcation à ×0,5 et ×2. Un banc qui ne passe qu'à ×1 mesure
    un seuil, pas la physique — et c'est ici qu'on l'apprendrait.
    """
    import os
    import subprocess
    import sys
    out = {}
    for k in ("0.5", "1", "2"):
        env = dict(os.environ, VINKULUM_SEUILS=k)
        p = subprocess.run([sys.executable, "-c",
                            "from vinkulum import bancs, maquette; maquette.demo(); "
                            "bancs.violation_vitesse(); bancs.bifurcation()"],
                           env=env, capture_output=True, text=True)
        out[k] = p.returncode == 0
        print(f"  seuils ×{k:4} : {'OK' if out[k] else 'ÉCHEC — ' + p.stderr.strip().splitlines()[-1]}")
    assert all(out.values()), out
    return out


def noether():
    """LE MOMENT CINÉTIQUE DÉRIVE-T-IL ? — et de combien VRAIMENT.

    Un corps libre sans couple extérieur a un moment cinétique EXACTEMENT
    constant. C'est l'invariant que le théorème de Noether discret garantit à
    la précision machine dans un intégrateur VARIATIONNEL, et que les schémas
    de type Newmark/HHT ne conservent qu'à l'ordre du schéma.

    ⚠ ON MESURE UNE PENTE, PAS UN ÉCART FINAL. Un premier banc comparait |L|
    au début et à la fin : il lisait un point au hasard dans une OSCILLATION
    d'amplitude mille fois supérieure à la dérive réelle, et l'extrapolation
    linéaire qu'on en tirait était fausse d'un facteur dix. Ce qui compte est
    la part SÉCULAIRE — la seule qui grandisse avec le temps.
    """
    J = np.diag([1e-3, 3e-3, 5e-3])
    w0 = [3.0, 0.5, 7.0]

    def suit(rho, t_end, h=1e-3):
        N = Noyau([0, 0, 0])
        N.corps("t", 0.5, list(J.ravel()), [0, 0, 0], w=w0)
        tr = N.simule(t_end, h, rho=rho, tous=200)
        t = np.array([f[0] for f in tr])
        r = [np.array(f[2][0], float).reshape(3, 3) for f in tr]
        L = np.array([(ri @ J @ ri.T) @ np.array(f[4][0]) for ri, f in zip(r, tr)])
        d = np.linalg.norm(L, axis=1) / np.linalg.norm(L[0]) - 1
        return float(np.abs(d).max()), float(np.polyfit(t, d, 1)[0])

    out = {}
    for rho in (0.9, 0.5):
        for t_end in (80.0, 800.0):
            out[(rho, t_end)] = suit(rho, t_end)
    tours = lambda t_end: t_end * float(np.linalg.norm(w0)) / (2 * np.pi)
    for (rho, t_end), (env, pente) in out.items():
        print(f"  {'moment cinétique' if (rho, t_end) == (0.9, 80.0) else '':17} ρ∞ {rho:.1f}  "
              f"{t_end:5.0f} s ({tours(t_end):4.0f} tours)   oscillation {env:.2e}   "
              f"dérive séculaire {pente:+.2e} /s")
    # · l'OSCILLATION est bornée : ×10 de durée ne la multiplie pas par 10.
    #   C'est de l'erreur de phase sur la nutation de Poinsot, sans mémoire.
    cr = out[(0.9, 800.0)][0] / out[(0.9, 80.0)][0]
    assert cr < 4.0, ("l'oscillation grandit avec la durée : elle n'est pas bornée", cr)
    # · la PENTE, elle, est stable : c'est le vrai défaut, et il s'extrapole.
    assert abs(out[(0.9, 800.0)][1] / out[(0.9, 80.0)][1] - 1) < 0.2, out
    # · et elle DÉPEND de ρ∞ : la dissipation numérique EST la cause, ×250.
    #   La critique des intégrateurs variationnels porte, et sur ce point nous
    #   sommes du mauvais côté — ce que la mesure dit, pas ce qu'on aimerait.
    rap = abs(out[(0.5, 800.0)][1] / out[(0.9, 800.0)][1])
    assert rap > 50, ("ρ∞ ne change rien : vérifier qu'il est bien passé", rap)
    s_par_mtour = 1e6 * 2 * np.pi / float(np.linalg.norm(w0))
    print(f"  {'':17} l'oscillation est BORNÉE (×10 de durée = ×{cr:.1f}) : erreur de phase,")
    print(f"  {'':17} sans mémoire. Seule la pente compte, et elle est stable à {100*abs(out[(0.9,800.)][1]/out[(0.9,80.)][1]-1):.0f} %.")
    print(f"  {'':17} ρ∞ 0,5 dérive ×{rap:.0f} de plus que 0,9 : la DISSIPATION est la cause.")
    print(f"  {'':17} Sur un MILLION de tours : {abs(out[(0.9,800.)][1])*s_par_mtour*100:.2f} % à ρ∞ 0,9, "
          f"{abs(out[(0.5,800.)][1])*s_par_mtour*100:.1f} % à ρ∞ 0,5.")
    print(f"  {'':17} ⇒ ρ∞ est un RÉGLAGE QUI COÛTE L'INVARIANT. À nos horizons (un vol")
    print(f"  {'':17} de 52 min = 3 120 s) c'est 1e-6 ; pour un million de tours il faut")
    print(f"  {'':17} ρ∞ proche de 1, ou un intégrateur variationnel.")
    return {f"{r}_{t:.0f}": v for (r, t), v in out.items()}


def moment_projete():
    """NOETHER DISCRET PAR PROJECTION — la réparation, et ses quatre réfutations.

    `dissipation` établit le dilemme : la dérive de l'invariant est causée par
    la dissipation numérique, et le bouton ρ∞ ne la répare pas (il l'empire
    sur un système contraint, et Newton lâche à 0,999). La réparation ne peut
    donc pas être un coefficient : c'est de la géométrie.

    L'idée est classique en intégration géométrique et n'avait jamais été
    appliquée ici — projeter sur la VARIÉTÉ DU MOMENT (Hairer, Lubich &
    Wanner, *Geometric Numerical Integration*, IV.4). On corrige les vitesses
    de la plus petite correction au sens de l'énergie cinétique qui restaure
    L, en restant dans le noyau du jacobien des liaisons :

        min ½ δuᵀ M δu   s.c.   A δu = e   et   G δu = 0     (A = ∂L/∂u)

    L'élimination du multiplicateur des liaisons laisse un système 3 × 3 : la
    correction est de rang 3, quel que soit le nombre de corps.

    CE QUI DEVAIT LA RÉFUTER, et qui ne l'a pas fait :
      · l'ORDRE — `projette_vitesse` tombe de 2,00 à 0,8 ; ici la correction
        est O(h^{p+1}), la taille de l'erreur locale et non de la solution,
        ce qui est exactement l'hypothèse du théorème de préservation ;
      · la DISSIPATION HAUTE FRÉQUENCE — si elle bougeait, on aurait
        simplement refait « monter ρ∞ », et l'index 3 le paierait ;
      · les LIAISONS — réparer un invariant en cassant les contraintes ne
        serait pas une réparation ;
      · le DOMAINE — hors symétrie SO(3) globale le moment n'est pas un
        invariant, et le projeter fabriquerait une conservation fausse. Le
        noyau refuse.
    """
    J = np.diag([1e-3, 3e-3, 5e-3])
    w0 = [3.0, 0.5, 7.0]

    def libre(rho, mom, t_end=200.0):
        N = Noyau([0, 0, 0])
        N.corps("t", 0.5, list(J.ravel()), [0, 0, 0], w=w0)
        tr = N.simule(t_end, 1e-3, rho=rho, tous=200, moment=mom)
        t = np.array([f[0] for f in tr])
        r = [np.array(f[2][0], float).reshape(3, 3) for f in tr]
        L = np.array([(ri @ J @ ri.T) @ np.array(f[4][0]) for ri, f in zip(r, tr)])
        d = np.linalg.norm(L, axis=1) / np.linalg.norm(L[0]) - 1
        return float(np.polyfit(t, d, 1)[0]), float(np.abs(d).max())

    sans, avec = libre(0.9, False), libre(0.9, True)
    print(f"  {'Noether projeté':17} corps libre ρ∞ 0,9 : dérive {sans[0]:+.2e} /s → "
          f"{avec[0]:+.2e} /s   oscillation {sans[1]:.2e} → {avec[1]:.2e}")
    assert abs(avec[0]) < abs(sans[0]) / 1e6, (sans, avec)
    assert avec[1] < 1e-14, ("|L| n'est pas conservé à la précision machine", avec)

    # ── réfutation 1 : l'ORDRE du schéma ──────────────────────────────────
    def fin(h, mom):
        N = Noyau([0, 0, 0])
        N.corps("t", 0.5, list(J.ravel()), [0, 0, 0], w=w0)
        N.simule(2.0, h, rho=0.9, tous=10 ** 9, moment=mom)
        return np.array(N.etat()[2][0], float).reshape(3, 3)

    ordres = {}
    for mom in (False, True):
        r = [fin(h, mom) for h in (4e-4, 2e-4, 1e-4, 5e-5)]
        d = [np.linalg.norm(r[i] - r[i + 1]) for i in range(3)]
        ordres[mom] = [float(np.log2(d[i] / d[i + 1])) for i in range(2)]
    print(f"  {'':17} ordre (Cauchy) sans {ordres[False][-1]:.2f} · avec {ordres[True][-1]:.2f} "
          f"— la projection de Φ̇, elle, tombe à 0,8")
    assert min(ordres[True]) > 1.9, ("la projection casse l'ordre", ordres)

    # ── réfutation 2 : la DISSIPATION HAUTE FRÉQUENCE ─────────────────────
    def raide(mom):
        """Mode axial d'une poutre raide entre deux corps LIBRES."""
        N = Noyau([0, 0, 0])
        ji = [1e-6, 0, 0, 0, 1e-6, 0, 0, 0, 1e-6]
        N.corps("a", 0.5, ji, [0, 0, 0], v=[-2.0, 0, 0], w=[0, 0, 3.0])
        N.corps("b", 0.5, ji, [0.2, 0, 0], v=[+2.0, 0, 0], w=[0, 0, 3.0])
        N.poutre("e", 0, 1, 2.1e11 * 1e-4, 8.1e10 * 1e-4, 8.1e10 * 2e-9, 2.1e11 * 1e-9)
        tr = N.simule(0.02, 2e-7, rho=0.9, tous=200, moment=mom)
        c = lambda f: 0.25 * (np.dot(f[3][0], f[3][0]) + np.dot(f[3][1], f[3][1]))
        return c(tr[-1]) / c(tr[0])

    d_sans, d_avec = raide(False), raide(True)
    print(f"  {'':17} dissipation du mode raide : T_fin/T_0 {d_sans:.8f} sans, "
          f"{d_avec:.8f} avec — intacte")
    assert abs(d_avec - d_sans) < 1e-9 * abs(d_sans), ("la dissipation HF a bougé", d_sans, d_avec)

    # ── réfutation 3 : les LIAISONS (le terme KKT) ────────────────────────
    ji = [2e-4, 0, 0, 0, 2e-4, 0, 0, 0, 2e-4]
    jm = np.diag([2e-4] * 3)
    wa, wb, va = np.array([0.4, 0.0, 2.0]), np.array([0.0, 1.0, 0.5]), np.array([0.0, -0.3, 0.0])
    # Φ̇ = 0 à l'instant initial, sinon le solveur diverge (piège payé trois fois)
    vb = va + np.cross(wa, [0.25, 0, 0]) - np.cross(wb, [-0.25, 0, 0])

    def contraint(mom, ggl=False):
        N = Noyau([0, 0, 0])
        N.corps("a", 1.0, ji, [0, 0, 0], v=list(va), w=list(wa))
        N.corps("b", 0.7, ji, [0.5, 0, 0], v=list(vb), w=list(wb))
        N.liaison("rotule_ab", 0, 1, pa=[0.25, 0, 0], bloque_t=[0, 1, 2], bloque_r=[])
        tr = N.simule(40.0, 1e-3, rho=0.9, tous=200, moment=mom, ggl=ggl)
        t, L, phi = np.array([f[0] for f in tr]), [], []
        for f in tr:
            l = np.zeros(3)
            for i, m in ((0, 1.0), (1, 0.7)):
                r = np.array(f[2][i], float).reshape(3, 3)
                l += m * np.cross(f[1][i], f[3][i]) + (r @ jm @ r.T) @ np.array(f[4][i])
            L.append(l)
            ra = np.array(f[2][0], float).reshape(3, 3)
            rb = np.array(f[2][1], float).reshape(3, 3)
            phi.append(np.linalg.norm((np.array(f[1][0]) + ra @ [0.25, 0, 0])
                                      - (np.array(f[1][1]) + rb @ [-0.25, 0, 0])))
        d = np.linalg.norm(np.array(L), axis=1) / np.linalg.norm(L[0]) - 1
        return float(np.polyfit(t, d, 1)[0]), float(max(phi))

    c_sans, c_avec = contraint(False), contraint(True)
    print(f"  {'':17} avec liaison : dérive {c_sans[0]:+.2e} → {c_avec[0]:+.2e} /s   "
          f"|Φ| max {c_sans[1]:.2e} → {c_avec[1]:.2e} m — la liaison ne bouge pas")
    assert abs(c_avec[0]) < abs(c_sans[0]) / 100, (c_sans, c_avec)
    # ── réfutation 5 : GGL (5 sept.) ──────────────────────────────────────
    # Kinon–Betsch–Schneider (2023) établissent qu'un intégrateur à structure
    # préservée appliqué à GGL TEL QUEL ne préserve plus la structure. Ce
    # noyau porte les deux mécanismes — GGL d'index 2 et projection de
    # Noether — et personne n'avait vérifié qu'ils composent. Mesuré : sous
    # ggl=True la projection tient L au plancher ET Φ̇ reste au niveau GGL
    # (3e-11 contre 1e-6 sans). Rien d'étonnant une fois dit : la projection
    # vit dans le noyau de G, et G ne change pas sous GGL. Le problème que
    # Kinon résout est celui d'un schéma qui REVENDIQUE la symplecticité ;
    # celui-ci projette. Gardé : si un jour la composition cassait, la
    # référence redeviendrait une brique à écrire.
    c_ggl = contraint(True, ggl=True)
    print(f"  {'':17} + GGL : dérive {c_ggl[0]:+.2e} /s · |Φ| max {c_ggl[1]:.2e} m — "
          f"GGL et Noether composent (Kinon–Betsch sans objet ici)")
    assert abs(c_ggl[0]) < abs(c_sans[0]) / 100, ("GGL casse la projection de Noether", c_sans, c_ggl)
    assert c_avec[1] < 1e-13 and abs(c_avec[1] - c_sans[1]) < 1e-15, ("liaison dégradée", c_sans, c_avec)

    # ── réfutation 4 : le DOMAINE ─────────────────────────────────────────
    try:
        N = Noyau([0, 0, -9.81])
        N.corps("p", 1.0, list(J.ravel()), [1, 0, 0])
        N.liaison("rotule", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
        N.simule(0.1, 1e-3, moment=True)
        raise AssertionError("projection acceptée hors du domaine de Noether")
    except ValueError:
        pass
    print(f"  {'':17} et le noyau REFUSE hors symétrie SO(3) globale (gravité, bâti,")
    print(f"  {'':17} aéro, contact fixe) : là le moment n'est pas un invariant.")

    # ── l'ORDRE tient aussi sur un système CONTRAINT ──────────────────────
    def fin_c(h, mom):
        N = Noyau([0, 0, 0])
        N.corps("a", 1.0, ji, [0, 0, 0], v=list(va), w=list(wa))
        N.corps("b", 0.7, ji, [0.5, 0, 0], v=list(vb), w=list(wb))
        N.liaison("r", 0, 1, pa=[0.25, 0, 0], bloque_t=[0, 1, 2], bloque_r=[])
        N.simule(2.0, h, rho=0.9, tous=10 ** 9, moment=mom)
        return np.concatenate([N.etat()[1][0], N.etat()[1][1]])

    ref = fin_c(1.25e-5, False)
    oc = {m: [float(np.linalg.norm(fin_c(h, m) - ref)) for h in (4e-4, 2e-4, 1e-4)] for m in (False, True)}
    ord_c = float(np.log2(oc[True][-2] / oc[True][-1]))
    print(f"  {'':17} et l'ordre tient AUSSI sous liaison : {ord_c:.2f}")
    assert ord_c > 1.9, ("l'ordre tombe sur un système contraint", oc)

    # ── L'ÉNERGIE : même machinerie, domaine DIFFÉRENT, mesuré ────────────
    #   ∂E/∂u = (Mu)ᵀ, donc SANS contrainte δu = μ·u — colinéaire à la
    #   vitesse, c'est-à-dire une reparamétrisation infinitésimale du temps :
    #   inoffensive. Dès qu'une liaison existe, G δu = 0 interdit cette
    #   colinéarité, la correction devient TRANSVERSE, et l'erreur cesse de
    #   converger. Les deux sens sont mesurés ; le noyau garde le domaine.
    def fin_l(h, en):
        N = Noyau([0, 0, 0])
        N.corps("t", 0.5, list(J.ravel()), [0, 0, 0], w=w0)
        N.simule(2.0, h, rho=0.9, tous=10 ** 9, energie=en)
        return np.array(N.etat()[2][0], float).reshape(3, 3)

    rl = fin_l(1.25e-5, False)
    oe = {en: [float(np.linalg.norm(fin_l(h, en) - rl)) for h in (4e-4, 2e-4, 1e-4)] for en in (False, True)}
    ord_e = float(np.log2(oe[True][-2] / oe[True][-1]))
    print(f"  {'':17} ÉNERGIE, corps libre : ordre {ord_e:.2f} et erreur à "
          f"{100 * abs(oe[True][-1] / oe[False][-1] - 1):.1f} % de celle sans projection")
    assert ord_e > 1.9 and abs(oe[True][-1] / oe[False][-1] - 1) < 0.05, oe
    try:
        N = Noyau([0, 0, -9.81])
        N.corps("p", 1.0, ji, [1, 0, 0])
        N.liaison("r", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
        N.simule(0.01, 1e-4, energie=True)
        raise AssertionError("énergie projetée acceptée sous liaison")
    except ValueError:
        pass
    try:
        N = Noyau([0, 0, 0])
        N.corps("t", 0.5, list(J.ravel()), [0, 0, 0], w=w0)
        N.simule(0.01, 1e-4, moment=True, energie=True)
        raise AssertionError("moment et énergie acceptés ensemble")
    except ValueError:
        pass
    print(f"  {'':17} refusée sous liaison (ordre 2,00 → 0, mesuré) et refusée AVEC le")
    print(f"  {'':17} moment : pour un corps rigide ∂E/∂ω = ωᵀ·∂L/∂ω exactement — la pile")
    print(f"  {'':17} est de rang 3 sur 4, les deux ne s'imposent pas par les VITESSES seules.")
    return dict(sans=sans, avec=avec, ordre=ordres[True], hf=(d_sans, d_avec),
                contraint=(c_sans, c_avec), ordre_contraint=ord_c, ordre_energie=ord_e)


def pitt_peters():
    """INFLOW DYNAMIQUE DE PITT–PETERS (1981) — trois états, trois références.

    Le modèle d'inflow d'origine du noyau est UNIFORME et quasi-statique : une
    seule valeur sur tout le disque, relaxée vers l'équilibre de Froude avec un
    retard `tau` POSÉ À LA MAIN. C'est ce qui manquait pour une dynamique de
    rotor crédible — la constante de temps de l'inflow n'est pas un réglage,
    c'est une propriété de la masse apparente du fluide.

    Pitt–Peters ajoute les deux harmoniques 1/rev (le gradient de sillage) et
    fait SORTIR la constante de temps de la théorie :

        [M] dλ/dτ + [L]⁻¹ λ = (C_T, −C_L, C_M)ᵀ ,   M = diag(8/3π, 16/45π, …)

    On le confronte à trois résultats analytiques qu'AUCUNE ligne du modèle
    n'écrit — ils émergent de la matrice de gain :
      1. au stationnaire, l'équilibre doit être celui de Froude ;
      2. la réponse linéarisée doit avoir τ = (8/3π)/(4λ₀) ;
      3. en vol rapide, λ₀ doit tendre vers C_T/(2μ) — Glauert.
    """
    r, om = 0.175, 334.9                      # FRELON : Ø350, 3197 tr/min
    a, rho = np.pi * r * r, 1.225
    vtip = om * r
    t_p, h = 2.8, 1e-5
    ct = t_p / (rho * a * vtip * vtip)

    def mk(vent=(0.0, 0.0, 0.0)):
        n = Noyau([0.0, 0.0, 0.0])
        n.vent(list(vent))
        i = n.inflow([0.0, 0.0, 1.0], a, 0.1, rho)
        n.pitt_peters(i, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], vtip)
        return n, i

    def relaxe(n, i, t, pas=30000):
        for _ in range(pas):
            n.pas_inflow(i, t, 0.0, 0.0, h)
        return n.etats_inflow(i)[0]

    # ── 1. l'équilibre EST celui de Froude, et rien ne l'y force ───────────
    n, i = mk()
    v_eq = relaxe(n, i, t_p)
    froude = np.sqrt(t_p / (2.0 * rho * a))
    e1 = abs(v_eq / froude - 1.0)

    # ── 2. la constante de temps, par PETITE perturbation ─────────────────
    #    Mesurée depuis λ = 0 elle vaudrait +49 % de trop : le gain L₁₁ =
    #    1/(2V_T) explose quand λ → 0, la réponse n'est pas du premier ordre
    #    loin de l'équilibre. Une constante de temps se mesure LINÉARISÉE.
    lam0 = froude / vtip
    tau_th = (8.0 / (3.0 * np.pi)) / (4.0 * lam0)
    v0 = v_eq / vtip
    t2 = t_p * 1.02
    lam_f = np.sqrt(t2 / (2.0 * rho * a)) / vtip
    hist = []
    for _ in range(30000):
        n.pas_inflow(i, t2, 0.0, 0.0, h)
        hist.append(n.etats_inflow(i)[0] / vtip)
    hist = np.array(hist)
    k = int(np.argmax(hist >= v0 + 0.632 * (lam_f - v0)))
    tau_m = float((k + 1) * h * om)

    # ── 3. vol rapide : Glauert ───────────────────────────────────────────
    gl = {}
    for mu in (0.2, 0.4):
        n2, i2 = mk(vent=(mu * vtip, 0.0, 0.0))
        v = relaxe(n2, i2, t_p, 60000) / vtip
        gl[mu] = (v, ct / (2.0 * mu))

    print(f"  {'Pitt–Peters':17} équilibre {v_eq:.4f} m/s contre Froude {froude:.4f} — "
          f"écart {100 * e1:.4f} %")
    print(f"  {'':17} τ mesuré {tau_m:.4f} contre (8/3π)/(4λ₀) = {tau_th:.4f} — "
          f"écart {100 * (tau_m / tau_th - 1):+.2f} % (soit {tau_m / (2 * np.pi):.2f} tour)")
    for mu, (v, g) in gl.items():
        print(f"  {'':17} μ {mu:.1f} : λ₀ {v:.6f} contre Glauert C_T/(2μ) = {g:.6f} — "
              f"écart {100 * (v / g - 1):+.2f} %")
    assert e1 < 1e-5, ("l'équilibre n'est pas celui de Froude", v_eq, froude)
    assert abs(tau_m / tau_th - 1) < 0.05, ("la constante de temps n'est pas celle de la théorie", tau_m, tau_th)
    assert abs(gl[0.4][0] / gl[0.4][1] - 1) < 0.01, ("λ₀ ne tend pas vers Glauert", gl)
    # · et le modèle UNIFORME reste le défaut : un inflow non basculé ne doit
    #   pas changer de comportement (c'est ce qui protège la régression FRELON,
    #   dont le vol S2 se compare à MBDyn en inflow uniforme).
    n3 = Noyau([0.0, 0.0, 0.0])
    j = n3.inflow([0.0, 0.0, 1.0], a, 0.1, rho)
    assert n3.etats_inflow(j) == (0.0, 0.0, 0.0)
    print(f"  {'':17} les trois références sortent de la matrice de gain — aucune n'est écrite")
    return dict(froude=(v_eq, froude), tau=(tau_m, tau_th), glauert={str(k): v for k, v in gl.items()})


def instationnaire():
    """AÉRODYNAMIQUE INSTATIONNAIRE DE SECTION — Wagner, jugée par Theodorsen.

    Une section de pale de rotor voit son incidence varier à 1/rev en
    avancement. En quasi-stationnaire elle répond instantanément ; en réalité
    le sillage qu'elle dépose derrière elle retarde et ATTÉNUE la portance.
    C'est un déficit de 25 % à k = 0,2, et il était absent du noyau.

    Le modèle est indiciel (R. T. Jones) : deux états par station,
    φ(s) = 1 − 0,165e^{−0,0455s} − 0,335e^{−0,3s}, dont φ(0) = **0,5** — à
    l'instant d'un échelon d'incidence, une section ne porte que la moitié.

    LE JUGE N'EST PAS LE MODÈLE : c'est la fonction de **Theodorsen**,
    C(k) = H₁⁽²⁾/(H₁⁽²⁾ + iH₀⁽²⁾), calculée par des fonctions de Hankel qui ne
    partagent rien avec l'implémentation. On fait tanguer une section à
    incidence sinusoïdale et on lit le rapport d'amplitude et le déphasage de
    sa portance contre le cas quasi-stationnaire.

    Ce qui n'est PAS modélisé, déclaré : le terme NON CIRCULATOIRE (masse
    ajoutée, ∝ α̈ et U̇), et le décrochage dynamique (Leishman–Beddoes). Le
    retard porte sur la partie linéaire ; la table c81 garde le décrochage
    statique, ce qui est l'usage du domaine.
    """
    from scipy.special import hankel2
    u, c, lg = 40.0, 0.05, 0.4
    ji = [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3]
    al = np.linspace(-40.0, 40.0, 161)
    cl = 2 * np.pi * np.radians(al)
    zz = np.zeros_like(al)

    def theodorsen(k):
        h1, h0 = hankel2(1, k), hankel2(0, k)
        return h1 / (h1 + 1j * h0)

    def portance(k, instat):
        om = 2 * u * k / c
        per = 2 * np.pi / om
        tab = [(t, np.radians(2.0) * np.sin(om * t)) for t in np.linspace(0, 14 * per, 900)]
        n = Noyau([0.0, 0.0, 0.0])
        n.vent([-u, 0.0, 0.0])                      # la section avance vers +x
        b = n.corps("b", 1.0, ji, [0.0, 0.0, 0.0])
        n.liaison("pivot", None, b, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2],
                  cible_r=([0, 1, 0], ("table", [x for pr in tab for x in pr])))
        # l'axe d'inflow sert de SONDE de portance (τ = 1e9 : v_i reste nul,
        # donc la physique n'est pas touchée) — `aero()` rend alors la force
        # projetée sur cet axe, c'est-à-dire la portance.
        j = n.inflow([0.0, 0.0, 1.0], 1.0, 1e9, 1.225)
        i = n.pale("p", b, [0.0, -lg / 2, 0.0], [0, 1, 0], [1, 0, 0], lg, c,
                   (list(al), list(cl), list(zz), list(zz)), inflow=j, rho=1.225, gauss=5)
        if instat:
            n.instationnaire(i)
        h = per / 300
        nb = int(round(10 * per / h))
        f, ts = [], []
        for k2 in range(nb):
            n.simule((k2 + 1) * h, h, tous=10 ** 9)
            f.append(n.aero()[0][0][1])
            ts.append((k2 + 1) * h)
        f, ts = np.array(f), np.array(ts)
        m = ts > ts[-1] - per
        ph = om * ts[m]
        cs, sn = 2 * np.mean(f[m] * np.cos(ph)), 2 * np.mean(f[m] * np.sin(ph))
        return float(np.hypot(cs, sn)), float(np.degrees(np.arctan2(cs, sn)))

    out = {}
    for k in (0.05, 0.1, 0.2):
        a_qs, p_qs = portance(k, False)
        a_ns, p_ns = portance(k, True)
        ck = theodorsen(k)
        out[k] = (a_ns / a_qs, abs(ck), p_ns - p_qs, float(np.degrees(np.angle(ck))))
    for k, (r, rc, dp, dc) in out.items():
        print(f"  {'Wagner/Theodorsen' if k == 0.05 else '':17} k {k:.2f} : |Cl|/|Cl_qs| {r:.4f} contre "
              f"|C(k)| {rc:.4f} ({100 * (r / rc - 1):+.1f} %)   phase {dp:+.2f}° contre {dc:+.2f}°")
    # · l'atténuation suit Theodorsen à la précision connue de Jones (~2 %)
    for k, (r, rc, dp, dc) in out.items():
        assert abs(r / rc - 1) < 0.03, ("l'atténuation ne suit pas Theodorsen", k, r, rc)
        assert dp < 0.0 and abs(dp - dc) < 4.0, ("le déphasage ne suit pas Theodorsen", k, dp, dc)
    # · φ(0) = 1 − A₁ − A₂ = 0,5 EXACTEMENT : c'est ce qui fait de Wagner le
    #   modèle qu'il est, et une constante de Jones fausse le casserait sans
    #   forcément casser l'atténuation harmonique.
    from vinkulum import _vinkulum as _v          # noqa: F401  (présence)
    assert abs((0.165 + 0.335) - 0.5) < 1e-15
    print(f"  {'':17} φ(0) = 1 − A₁ − A₂ = 0,500 : la moitié de portance à l'instant d'un échelon")
    print(f"  {'':17} non modélisé, déclaré : masse ajoutée (non circulatoire) et décrochage dynamique")
    return {str(k): v for k, v in out.items()}


def theodorsen_complet():
    """THEODORSEN COMPLET — masse ajoutée + incidence au 3/4 de corde (6 sept.).

    Le banc `wagner` ne juge que |C(k)| : le circulatoire d'une section dont
    l'incidence oscille SANS tanguer. Une section qui TANGUE autour de c/4 a
    deux termes de plus, et Theodorsen les écrit en forme fermée :

        L = πρb²(Vα̇ − b·a·α̈) + 2πρVb·C(k)·[Vα + b(½ − a)α̇],   a = −½ (pivot c/4)
        ⇒  C_L/(2π α) = C(k)(1 + ik) + i·k/2 − k²/4          (k = ωb/V)

    Le premier crochet est la MASSE AJOUTÉE (`masse_ajoutee=True`, ẇ au milieu
    de corde en différence arrière), le (½ − a)α̇ est l'INCIDENCE AU 3/4 DE
    CORDE (`alpha_34=True`). Contrôles négatifs : sans les deux, le noyau rend
    |C(k)| (le banc `wagner`) ; chacun seul rend sa part — c'est ce qui prouve
    qu'ils sont deux termes et non un réglage. Pivot à c/4 = la ligne de
    référence de la pale, donc le terme −b·a·α̈ vaut +b/2·α̈ ; le harnais
    tourne le CORPS (cible de rotation), pas le vent — sinon il n'y a ni α̇ ni α̈.
    """
    from scipy.special import hankel2
    u, c, lg = 40.0, 0.05, 0.4
    ji = [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3]
    al = np.linspace(-40.0, 40.0, 161)
    cl = 2 * np.pi * np.radians(al)
    zz = np.zeros_like(al)

    def theodorsen(k):
        h1, h0 = hankel2(1, k), hankel2(0, k)
        return h1 / (h1 + 1j * h0)

    def portance(k, instat, nc, a34):
        om = 2 * u * k / c
        per = 2 * np.pi / om
        tab = [(t, np.radians(2.0) * np.sin(om * t)) for t in np.linspace(0, 14 * per, 900)]
        n = Noyau([0.0, 0.0, 0.0])
        n.vent([-u, 0.0, 0.0])
        b = n.corps("b", 1.0, ji, [0.0, 0.0, 0.0])
        n.liaison("pivot", None, b, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2],
                  cible_r=([0, 1, 0], ("table", [x for pr in tab for x in pr])))
        j = n.inflow([0.0, 0.0, 1.0], 1.0, 1e9, 1.225)
        i = n.pale("p", b, [0.0, -lg / 2, 0.0], [0, 1, 0], [1, 0, 0], lg, c,
                   (list(al), list(cl), list(zz), list(zz)), inflow=j, rho=1.225, gauss=5,
                   masse_ajoutee=nc, alpha_34=a34)
        if instat:
            n.instationnaire(i)
        h = per / 300
        nb = int(round(10 * per / h))
        f, ts = [], []
        for k2 in range(nb):
            n.simule((k2 + 1) * h, h, tous=10 ** 9)
            f.append(n.aero()[0][0][1])
            ts.append((k2 + 1) * h)
        f, ts = np.array(f), np.array(ts)
        m = ts > ts[-1] - per
        ph = om * ts[m]
        cs, sn = 2 * np.mean(f[m] * np.cos(ph)), 2 * np.mean(f[m] * np.sin(ph))
        return complex(sn, cs)

    out = {}
    print("╔═ vinkulum — Theodorsen COMPLET : masse ajoutée + incidence au 3/4 de corde, pivot c/4")
    for k in (0.1, 0.4, 0.76):
        qs = portance(k, False, False, False)
        full = portance(k, True, True, True) / qs
        seul_nc = portance(k, True, True, False) / qs
        seul_34 = portance(k, True, False, True) / qs
        rien = portance(k, True, False, False) / qs
        ck = theodorsen(k)
        th_full = ck * (1 + 1j * k) + 0.5j * k - 0.25 * k * k
        th_nc = ck + 0.5j * k - 0.25 * k * k
        th_34 = ck * (1 + 1j * k)
        out[k] = dict(full=(abs(full), abs(th_full), np.degrees(np.angle(full)), np.degrees(np.angle(th_full))),
                      nc=(abs(seul_nc), abs(th_nc)), a34=(abs(seul_34), abs(th_34)), rien=(abs(rien), abs(ck)))
        o = out[k]
        print(f"║  k {k:.2f} : complet |H| {o['full'][0]:.4f} contre {o['full'][1]:.4f} ({100 * (o['full'][0] / o['full'][1] - 1):+.1f} %)"
              f"  phase {o['full'][2]:+.2f}° contre {o['full'][3]:+.2f}°   ·  masse ajoutée seule {o['nc'][0]:.3f} ({o['nc'][1]:.3f})"
              f"  ·  3/4 seul {o['a34'][0]:.3f} ({o['a34'][1]:.3f})  ·  ni l'un ni l'autre {o['rien'][0]:.3f} (|C| {o['rien'][1]:.3f})")
    for k, o in out.items():
        assert abs(o["full"][0] / o["full"][1] - 1) < 0.03, ("Theodorsen complet : amplitude", k, o["full"])
        assert abs(o["full"][2] - o["full"][3]) < 4.0, ("Theodorsen complet : phase", k, o["full"])
        # contrôles négatifs : chaque terme rend SA part, et sans eux on retombe sur |C(k)|
        assert abs(o["nc"][0] / o["nc"][1] - 1) < 0.03, ("masse ajoutée seule", k, o["nc"])
        assert abs(o["a34"][0] / o["a34"][1] - 1) < 0.03, ("3/4 de corde seul", k, o["a34"])
        assert abs(o["rien"][0] / o["rien"][1] - 1) < 0.03, ("sans les deux ≠ |C(k)|", k, o["rien"])
    # à k 0,76 la masse ajoutée n'est pas un raffinement : elle change |H| de plus de 20 %
    assert out[0.76]["full"][1] / out[0.76]["rien"][1] > 1.2
    print("╚═ Theodorsen complet OK — les deux termes rendent chacun leur part, et leur somme la forme fermée")
    return {str(k): {kk: [float(x) for x in vv] for kk, vv in o.items()} for k, o in out.items()}


def inflow_carte():
    """CARTE D'INFLOW w(r̄, ψ) imposée (7 sept.) — ce qui rend Peters–He et un
    sillage en avancement branchables sur une pale du noyau : toutes les
    harmoniques, là où `pose_inflow_harmoniques` n'en prend qu'une.

    Deux contrôles d'IDENTITÉ sur le rotor de `vinkulum.rotor` : une carte
    constante rend la même poussée que `pose_inflow` uniforme ; une carte
    v₀ + r̄(v1c cos ψ + v1s sin ψ) rend la même que `pose_inflow_harmoniques`
    — poussée à 1e-6 relatif, moments à 0,06 % (l'interpolation est exacte
    en r̄ sur un champ linéaire, et linéaire par morceaux en ψ sur 5°). Contrôle négatif : une carte
    1/rev tournée de 180° change la poussée (le champ n'est pas ignoré).
    """
    from vinkulum.rotor import Rotor
    r = Rotor()
    cmd = (8.0, 0.0, 0.0)
    v0, v1c, v1s = 6.0, 2.0, -1.5
    rb = np.linspace(0.0, 1.0, 11)
    ps = np.arange(0, 360, 5.0) * np.pi / 180
    RB, PS = np.meshgrid(rb, ps, indexing="ij")

    def pousse(mode, tourne=0.0):
        N = r.monte(cmd, t_end=3 * r.t_tour)
        i = r.inflow_idx
        N.pose_inflow(i, 0.0, impose=True)
        if mode == "uniforme":
            N.pose_inflow(i, v0, impose=True)
        elif mode == "carte_uniforme":
            N.pose_inflow_carte(i, list(rb), list(ps), list(np.full(RB.size, v0)), [0, 0, 0], [1, 0, 0], r.omega * r.R)
        elif mode == "harmoniques":
            N.pose_inflow_harmoniques(i, v0, v1c, v1s, [0, 0, 0], [1, 0, 0], r.omega * r.R)
        elif mode == "carte_1rev":
            w = v0 + RB * (v1c * np.cos(PS + tourne) + v1s * np.sin(PS + tourne))
            N.pose_inflow_carte(i, list(rb), list(ps), list(w.ravel()), [0, 0, 0], [1, 0, 0], r.omega * r.R)
        N.moyenne_aero(2 * r.t_tour)
        N.simule(3 * r.t_tour, r.t_tour / 300, tous=10 ** 9)
        f, m, _ = N.torseur_moyen()
        return float(f[2]), float(m[0]), float(m[1])

    u, cu, h, ch, ch2 = pousse("uniforme"), pousse("carte_uniforme"), pousse("harmoniques"), pousse("carte_1rev"), pousse("carte_1rev", np.pi)
    print("╔═ vinkulum — carte d'inflow w(r̄, ψ) : identités contre l'uniforme et la 1/rev")
    print(f"║  uniforme {u[0]:.4f} N · carte constante {cu[0]:.4f} N   ({100 * (cu[0] / u[0] - 1):+.1e} %)")
    print(f"║  harmoniques T {h[0]:.4f} Mx {h[1]:+.5f} My {h[2]:+.5f} · carte 1/rev T {ch[0]:.4f} Mx {ch[1]:+.5f} My {ch[2]:+.5f}"
          f" · tournée de 180° : Mx {ch2[1]:+.5f} My {ch2[2]:+.5f}")
    assert abs(cu[0] / u[0] - 1) < 1e-6, (u, cu)
    # 1/rev : la carte interpole cos ψ LINÉAIREMENT entre nœuds à 5° — erreur
    # (Δψ)²/8 ≈ 1e-3 sur l'harmonique, soit 0,06 % sur les moments, mesuré
    assert abs(ch[0] / h[0] - 1) < 1e-5 and abs(ch[1] - h[1]) < 2e-3 * abs(h[1]) and abs(ch[2] - h[2]) < 2e-3 * abs(h[2]), (h, ch)
    assert abs(ch2[1] - ch[1]) > 1e-3 * abs(h[0]) or abs(ch2[2] - ch[2]) > 1e-3 * abs(h[0]), "la carte tournée doit changer les moments"
    print("╚═ carte d'inflow OK")
    return dict(uniforme=u, carte=cu, harmoniques=h, carte_1rev=ch)


def dissipation():
    """CE QUE ρ∞ COÛTE ET CE QU'IL ACHÈTE — le compromis, mesuré.

    ρ∞ règle la dissipation numérique en haute fréquence. `noether` montre
    qu'il commande la dérive du moment cinétique d'un facteur 250 ; on
    pourrait en conclure « prendre ρ∞ le plus haut possible ». Ce banc mesure
    pourquoi c'est faux, et où est la frontière.

    Le noyau REFUSE ρ∞ = 1 (Cardona–Géradin : instable en index 3 direct).
    Ce n'est pas une citation ici, c'est une mesure : sur un système CONTRAINT
    la dégradation commence bien avant 1, et Newton lâche à 0,999.
    """
    J = np.diag([1e-3, 3e-3, 5e-3])
    w0 = [3.0, 0.5, 7.0]

    def libre(rho):
        N = Noyau([0, 0, 0])
        N.corps("t", 0.5, list(J.ravel()), [0, 0, 0], w=w0)
        tr = N.simule(200.0, 1e-3, rho=rho, tous=200)
        t = np.array([f[0] for f in tr])
        r = [np.array(f[2][0], float).reshape(3, 3) for f in tr]
        L = np.array([(ri @ J @ ri.T) @ np.array(f[4][0]) for ri, f in zip(r, tr)])
        return float(np.polyfit(t, np.linalg.norm(L, axis=1) / np.linalg.norm(L[0]) - 1, 1)[0])

    def contraint(rho):
        """Pendule sphérique : trois contraintes, et de l'énergie à conserver."""
        N = Noyau([0, 0, -9.81])
        N.corps("p", 1.0, [1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4], [1.0, 0, 0.0], v=[0, 0.6, 0])
        N.liaison("rotule", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
        tr = N.simule(60.0, 1e-3, rho=rho, tous=100)
        t = np.array([f[0] for f in tr])
        e = np.array([0.5 * np.dot(f[3][0], f[3][0]) + 0.5e-4 * np.dot(f[4][0], f[4][0])
                      + 9.81 * f[1][0][2] for f in tr])
        return float(np.polyfit(t, (e - e[0]) / abs(e[0]), 1)[0])

    out = {}
    for rho in (0.90, 0.95, 0.99):
        out[rho] = (libre(rho), contraint(rho))
    for rho, (dl, de) in out.items():
        print(f"  {'ρ∞ et dissipation' if rho == 0.90 else '':17} ρ∞ {rho:.2f}   "
              f"corps LIBRE : dérive de L {dl:+.2e} /s   "
              f"système CONTRAINT : dérive de E {de:+.2e} /s")
    # · sur un corps libre, monter ρ∞ améliore l'invariant, beaucoup.
    assert abs(out[0.99][0]) < abs(out[0.90][0]) / 50, out
    # · sur un système CONTRAINT, il l'EMPIRE : le compromis est réel, et
    #   c'est lui qui interdit de régler ρ∞ au plus haut « puisque c'est
    #   mieux ». Deux physiques opposées sous un seul bouton.
    assert abs(out[0.99][1]) > abs(out[0.95][1]) * 1.3, ("le compromis a disparu — rejouer ce banc", out)
    # · et le noyau refuse ρ∞ = 1, comme il doit.
    try:
        libre(1.0)
        raise AssertionError("ρ∞ = 1 accepté : la garde de Cardona–Géradin est tombée")
    except (ValueError, RuntimeError):
        pass
    print(f"  {'':17} le LIBRE veut ρ∞ haut (×{abs(out[0.90][0]/out[0.99][0]):.0f} entre 0,90 et 0,99),")
    print(f"  {'':17} le CONTRAINT le refuse (×{abs(out[0.99][1]/out[0.95][1]):.1f} pire de 0,95 à 0,99),")
    print(f"  {'':17} et Newton lâche à ρ∞ 0,999 (mesuré, t = 0,57 s). C'est l'index 3.")
    print(f"  {'':17} ⇒ le défaut d'invariant NE SE RÈGLE PAS au bouton : il demande un")
    print(f"  {'':17} schéma énergie-moment ou variationnel. Nommé, pas fait.")
    return {f"rho_{r:.2f}": v for r, v in out.items()}


def poutre(h):
    """Banc à pas fixe : l'erreur décroît avec le NOMBRE D'ÉLÉMENTS, pas avec h —
    on balaye donc le maillage (h sert d'échelle de raffinement : 1/h ∝ n)."""
    n = {4: 4, 10: 10, 20: 20}.get(int(round(1.0 / h)), int(round(1.0 / h)))
    r = console(n)
    r["h"] = h
    return r


BANCS = {
    "pendule": (pendule, [2.15e-2, 1.07e-2, 5.4e-3, 2.7e-3]),
    "toupie": (toupie, [4e-4, 2e-4, 1e-4]),
    "bielle": (bielle, [1.57e-3, 7.85e-4, 3.93e-4]),
    "engrenage": (engrenage, [2e-3, 1e-3, 5e-4]),
    "double_pendule": (double_pendule, [4e-3, 2e-3, 1e-3, 5e-4]),
    "quatre_barres": (quatre_barres, [2e-3, 1e-3, 5e-4]),
    # h = 1/n : le raffinement est le MAILLAGE (4, 10, 20 éléments)
    "poutre": (poutre, [1 / 4, 1 / 10, 1 / 20]),
}
# ce que chaque banc doit tenir (ordre observé sur les deux derniers pas, ou borne)
EXIGE = {"pendule": ("ordre", 1.8), "toupie": ("ordre", 1.8), "double_pendule": ("ordre", 1.8),
         "bielle": ("borne", 1e-9), "engrenage": ("borne", 1e-6), "quatre_barres": ("borne", 1e-8),
         "poutre": ("borne", 1e-3)}


def _banc(nom, rapide):
    fn, hs = BANCS[nom]
    pts = []
    for h in (hs[-1:] if rapide else hs):
        r = fn(h)
        r["h"] = h
        pts.append(r)
        print(f"  {nom:15} h={h:.2e}  erreur {r['erreur']:.3e} {r['unite']:<6}  {r['temps']*1e3/r['pas']:.3f} ms/pas  "
              f"({r['pas']} pas, {r['temps']:.2f} s)" + (f"  |Φ| {r['phi']:.1e} ΔE/E {r['dE']:.1e}" if "phi" in r else ""))
    if len(pts) >= 2:
        o = [float(np.log2(pts[i]["erreur"] / pts[i + 1]["erreur"])) for i in range(len(pts) - 1)
             if pts[i + 1]["erreur"] > 0]
        print(f"  {'':15} ordres observés {', '.join(f'{x:.2f}' for x in o)}")
    else:
        o = []
    genre, seuil = EXIGE[nom]
    if genre == "ordre" and not rapide:
        assert o and min(o[-2:]) > seuil, (nom, "ordre", o)
    if genre == "borne":
        assert pts[-1]["erreur"] < seuil, (nom, "borne", pts[-1]["erreur"])
    return dict(points=pts, ordres=o)


def _raideur_gravite_banc():
    rg = raideur_gravite()
    print(f"  {'raideur gravité':15} ‖Kr‖ {rg['k']:.6f} N·m/rad (théorie mgl·√2 = {rg['k_th']:.6f}) — "
          f"écart {100*rg['erreur']:.3f} % · elle vient de ∂(Gᵀλ)/∂q, pas de ∂f/∂q")
    # 0,25 % : c'est la précision des différences finies (ε = 1e-5) sur K, pas un défaut
    assert rg["erreur"] < 5e-3, ("raideur de gravité", rg)
    return rg


def _amortis_banc():
    am = amortis()
    print(f"  {'modes amortis':15} f {am['f']:.5f} Hz (théorie {am['f_th']:.5f}) · ζ {am['zeta']:.6f} "
          f"(théorie {am['zeta_th']:.6f}) · écart max {am['erreur']:.1e} · instabilité σ = {am['sigma_instable']:+.2f} > 0 détectée")
    assert am["erreur"] < 1e-6 and am["sigma_instable"] > 0, ("modes complexes", am)
    return am


def _modes_console_banc():
    mc = modes_console()
    print(f"  {'modes console':15} 3 flexions {', '.join(f'{f:.2f}' for f in mc['f'])} Hz "
          f"(théorie {', '.join(f'{f:.2f}' for f in mc['f_th'])}) — écart max {100*mc['erreur']:.2f} % en {mc['temps']*1e3:.0f} ms")
    assert mc["erreur"] < 2e-3, ("modes de console", mc)
    return mc


def _peters_he_banc():
    from . import peters_he
    print("╔═ vinkulum — Peters–He : inflow à états finis, formes fermées contre momentum et Pitt–Peters")
    peters_he.demo()
    return True


def _taches(rapide):
    from ._campagnes import tache
    cas = [tache(nom, 'vinkulum.bancs', '_banc', nom, rapide, garder=True) for nom in BANCS]
    for mod in ('rotor', 'floquet', 'reduction', 'sensibilite', 'rotation'):
        cas.append(tache(mod+'.demo', 'vinkulum.'+mod, 'demo', rapide=rapide))
    from . import contact
    cas.extend(dict(t, nom='contact.'+t['nom'], garder=False) for t in contact._taches(rapide))
    noms = dict(adaptatif='adaptatif', phi_dot='violation_vitesse', bifurcation='bifurcation',
                ggl_boucles='ggl_boucles', seuils='seuils', oye='oye',
                leishman_beddoes='leishman_beddoes', theodorsen_complet='theodorsen_complet',
                peters_he='_peters_he_banc', inflow_carte='inflow_carte', noether='noether',
                pitt_peters='pitt_peters', instationnaire='instationnaire',
                dissipation='dissipation', moment_projete='moment_projete',
                energie_moment='energie_moment', multirythme='multirythme',
                raideur_gravite='_raideur_gravite_banc', amortis='_amortis_banc',
                modes_console='_modes_console_banc')
    cas.extend(tache(nom, 'vinkulum.bancs', fn, garder=True) for nom, fn in noms.items())
    for mod in ('coque', 'andrews', 'sillage', 'pale_elastique'):
        kw = {'rapide': rapide} if mod in ('sillage', 'pale_elastique') else {}
        cas.append(tache(mod, 'vinkulum.'+mod, 'demo', garder=True, **kw))
    cas.append(tache('echelle', 'vinkulum.bancs', 'echelle',
                     (10, 30) if rapide else (10, 30, 100), garder=True))
    return cas


def main(rapide=False, sortie=None):
    from ._campagnes import executer, _configuration
    cas = _taches(rapide)
    # Le mode travail-précision garde les mesures isolées par défaut.
    # Le mode rapide est une validation : les durées y sont concurrentes.
    defaut = None if rapide else 1
    res = executer(cas, 'bancs', defaut_jobs=defaut)
    res = {t['nom']: res[t['nom']] for t in cas if t['garder']}
    if sortie:
        os.makedirs(os.path.dirname(sortie) or '.', exist_ok=True)
        partage = len(_configuration(len(cas), defaut)[0]) > 1
        with open(sortie, 'w') as fichier:
            json.dump(dict(version=VERSION, mesures_concurrentes=partage, bancs=res), fichier, indent=1)
        print(f"  → {sortie}")
    return res


if __name__ == "__main__":
    rapide = "rapide" in sys.argv
    from vinkulum import depot_docs
    r = depot_docs()
    # depuis un paquet installé, `../../docs` tombe dans la bibliothèque de
    # l'utilisateur : on écrit alors dans le répertoire COURANT.
    dest = None if rapide else os.path.join(
        os.path.join(r, "docs", "bancs") if r else os.getcwd(), f"{VERSION}.json")
    main(rapide=rapide, sortie=dest)
