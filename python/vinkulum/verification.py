"""vinkulum.verification — la vérification APPROFONDIE du jalon 1 (hors gate, minutes).

    python -m vinkulum.verification

La démo (`maquette.py`) juge trois cas sur des grandeurs intégrales (période,
précession moyenne, x(θ)). Ici on juge ce que la théorie PROMET, chiffre par
chiffre, contre des références indépendantes :

1. ORDRE DE CONVERGENCE — l'α-généralisé est annoncé d'ordre 2 en positions et
   vitesses (Arnold–Brüls 2007). Pendule à 30° (non linéaire), référence :
   l'EDO minimale θ̈ = −(g/L)·sin θ intégrée à 1e-12. Erreur en h, h/2, h/4,
   h/8 → les rapports doivent tendre vers 4. Idem pour le multiplicateur
   (tension du fil = m·(g·cos θ + L·θ̇²)).
2. CONTRAINTE EN VITESSE — l'index 3 direct tient Φ(q) = 0 ; il ne tient PAS
   Φ̇ = G·u = 0 (c'est l'argument des formulations GGL / index 2). On MESURE
   la violation r·v sur le pendule, et son ordre en h.
3. TOUPIE CONTRE EDO — trajectoire complète de l'axe (pas la précession
   moyenne) contre les équations de Lagrange de la toupie symétrique en angles
   d'Euler, intégrées à 1e-11. Erreur à deux pas, rapport.
4. DÉTERMINISME — un fil ou tous : résultat identique au bit près (chaque
   colonne du jacobien est indépendante, la réduction ne mélange rien).
5. AMORTISSEMENT NUMÉRIQUE — ρ∞ = 1 est REFUSÉ (trapèze : instable en
   index 3, mesuré — Cardona–Géradin 1989) ; à 0,95 / 0,9 / 0,5 la
   dissipation d'énergie doit croître quand ρ∞ baisse.

Ce que ce fichier NE vérifie pas, déclaré : la stabilité sur un système
raide (pas de corps flexible au jalon 1) ; le comportement près d'une rotation
relative de π (le résidu de rotation vee(skew(E)) s'annule aussi à E = rotation
de π — solution parasite hors bassin de Newton, documentée) ; les contraintes
redondantes (moindres carrés testés seulement sur un cas non redondant).
"""
import re
import os
import sys
import time
from math import cos, degrees, pi, radians, sin, sqrt

import numpy as np
from scipy.integrate import solve_ivp

from vinkulum import Noyau

G = 9.80665
I3 = np.eye(3)


def _rot9(R):
    return [float(x) for x in np.asarray(R).reshape(9)]


def _expm(w):
    t = float(np.linalg.norm(w))
    K = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]], float)
    if t < 1e-12:
        return I3 + K
    return I3 + sin(t) / t * K + (1 - cos(t)) / t**2 * K @ K


# ── 1 & 2 & 5 : pendule ──────────────────────────────────────────────────────
def _pendule(h, t_end, rho=0.9, th0=np.radians(30.0), L=1.0, m=0.2):
    N = Noyau()
    N.corps("masse", m, _rot9(np.diag([1e-8] * 3)), [L * sin(th0), 0.0, -L * cos(th0)])
    N.liaison("rotule", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
    E0 = N.energie()
    tr = N.simule(t_end, h, rho=rho)
    return N, tr, E0


def _pendule_ref(t_end, th0=np.radians(30.0), L=1.0):
    sol = solve_ivp(lambda t, y: [y[1], -G / L * sin(y[0])], (0, t_end), [th0, 0.0],
                    rtol=1e-12, atol=1e-14, dense_output=True)
    return sol.sol


def convergence():
    L, m, th0 = 1.0, 0.2, np.radians(30.0)
    T = 2 * pi * sqrt(L / G) * 1.07                      # ~ période à 30°
    t_end = 2 * T
    ref = _pendule_ref(t_end, th0, L)
    hs = [T / 100, T / 200, T / 400, T / 800]
    err_q, err_v, err_lam, viol = [], [], [], []
    for h in hs:
        N, tr, _ = _pendule(h, t_end, th0=th0)
        t = tr[-1][0]
        th, thd = ref(t)
        r_ref = np.array([L * sin(th), 0.0, -L * cos(th)])
        v_ref = np.array([L * thd * cos(th), 0.0, L * thd * sin(th)])
        r, v = np.array(tr[-1][1][0]), np.array(tr[-1][3][0])
        err_q.append(np.linalg.norm(r - r_ref))
        err_v.append(np.linalg.norm(v - v_ref))
        # tension du fil : multiplicateur de la rotule (3 composantes = force sur la masse)
        lam = np.array(tr[-1][5])
        tension_ref = m * (G * cos(th) + L * thd**2)
        err_lam.append(abs(np.linalg.norm(lam) - tension_ref))
        # contrainte en vitesse : r·v = 0 exactement pour |r| = L
        viol.append(max(abs(np.dot(e[1][0], e[3][0])) / (L * max(np.linalg.norm(e[3][0]), 1e-12))
                        for e in tr))
    def ordres(e):
        return [np.log2(e[i] / e[i + 1]) for i in range(len(e) - 1)]
    oq, ov, ol, ovi = ordres(err_q), ordres(err_v), ordres(err_lam), ordres(viol)
    print("  1. ordre de convergence (pendule 30°, 2 périodes, contre EDO minimale à 1e-12)")
    print(f"     h = T/{'  T/'.join(str(int(round(T / h))) for h in hs)}")
    print(f"     erreur position   {'  '.join(f'{e:.2e}' for e in err_q)}   ordres {', '.join(f'{o:.2f}' for o in oq)}")
    print(f"     erreur vitesse    {'  '.join(f'{e:.2e}' for e in err_v)}   ordres {', '.join(f'{o:.2f}' for o in ov)}")
    print(f"     erreur tension λ  {'  '.join(f'{e:.2e}' for e in err_lam)}   ordres {', '.join(f'{o:.2f}' for o in ol)}")
    print(f"  2. violation r·v/(L|v|) {'  '.join(f'{e:.2e}' for e in viol)}   ordres {', '.join(f'{o:.2f}' for o in ovi)}")
    assert min(oq[-2:]) > 1.8 and min(ov[-2:]) > 1.8, ("ordre 2 non tenu", oq, ov)
    assert min(ol[-2:]) > 1.5, ("multiplicateur : ordre", ol)
    # MESURÉ : la violation en vitesse est O(h²) — 5,7e-4 à T/100, 9,0e-6 à
    # T/800. C'est la limite CONNUE de l'index 3 direct (Φ tenu, Φ̇ pas) et
    # l'argument des formulations GGL / index 2 stabilisé. On exige l'ordre 2
    # et une borne à T/800 ; on ne prétend pas qu'elle est nulle.
    assert min(ovi) > 1.8 and viol[-1] < 2e-5, ("contrainte en vitesse", viol, ovi)
    return dict(h=hs, err_q=err_q, err_v=err_v, err_lam=err_lam, viol=viol)


def amortissement():
    """ρ∞ = 1 (trapèze) est INSTABLE en index 3 — mesuré : Newton diverge à
    t ≈ 16 s (résidu 1,9e5) après des oscillations croissantes de λ. C'est le
    résultat de Cardona–Géradin (1989) qui a motivé HHT-α puis l'α-généralisé.
    Le noyau le REFUSE ; en dessous, la dissipation doit croître quand ρ∞ baisse."""
    L, th0 = 1.0, np.radians(30.0)
    T = 2 * pi * sqrt(L / G) * 1.07
    print("  5. amortissement numérique (pendule, 10 périodes, h = T/400)")
    try:
        _pendule(T / 400, 10 * T, rho=1.0, th0=th0)
        raise AssertionError("ρ∞ = 1 doit être refusé (instable en index 3)")
    except (ValueError, RuntimeError) as e:
        # ValueError depuis le 3 sept. : un domaine d'ARGUMENT invalide se
        # refuse comme les gardes d'invariant. RuntimeError reste attrapé —
        # le solveur refuse aussi, si l'appel le contourne.
        print(f"     ρ∞ = 1     refusé : {e}")
    out = {}
    for rho in (0.95, 0.9, 0.5):
        N, tr, E0 = _pendule(T / 400, 10 * T, rho=rho, th0=th0)
        dE = (N.energie() - E0) / (0.2 * G * L * (1 - cos(th0)))
        out[rho] = dE
        print(f"     ρ∞ = {rho:<4}  ΔE/E = {dE:+.2e}")
    assert out[0.5] < out[0.9] < out[0.95] < 0, ("la dissipation doit croître quand ρ∞ baisse", out)
    return out


# ── 3 : toupie contre l'EDO de Lagrange ──────────────────────────────────────
def _toupie_ref(t_end, m, l, Ja, Jt, ws, th0):
    """Toupie symétrique, angles d'Euler (θ nutation, φ précession, ψ spin),
    Jt' = Jt + m l² au pivot. Lagrangien classique ; on intègre θ̈, φ̇ et ψ̇
    depuis les intégrales premières p_φ et p_ψ."""
    Jtp = Jt + m * l**2
    # à t = 0 : ω = ws·e_axe (pas de précession, pas de nutation) → ψ̇ = ws, φ̇ = 0, θ̇ = 0
    p_psi = Ja * ws                                   # = Ja (ψ̇ + φ̇ cos θ)
    p_phi = Ja * ws * cos(th0)                        # = Jt' sin²θ φ̇ + p_psi cos θ

    def f(t, y):
        th, thd = y
        s, c = sin(th), cos(th)
        phid = (p_phi - p_psi * c) / (Jtp * s * s)
        thdd = (phid**2 * s * c * Jtp - p_psi * phid * s + m * G * l * s) / Jtp
        return [thd, thdd]
    sol = solve_ivp(f, (0, t_end), [th0, 0.0], rtol=1e-11, atol=1e-13, dense_output=True)

    def axe(t):
        th = sol.sol(t)[0]
        # φ par intégration de φ̇ : on le reconstruit par quadrature fine
        return th
    # φ(t) par quadrature de φ̇(θ(t))
    ts = np.linspace(0, t_end, 200001)
    ths = sol.sol(ts)[0]
    phid = (p_phi - p_psi * np.cos(ths)) / (Jtp * np.sin(ths) ** 2)
    phis = np.concatenate([[0.0], np.cumsum(0.5 * (phid[1:] + phid[:-1]) * np.diff(ts))])
    return ts, ths, phis


def toupie():
    m, l, Ja, Jt, ws, th0 = 0.1, 0.05, 1e-4, 3e-4, 300.0, np.radians(30.0)
    t_end = 1.0
    ts, ths, phis = _toupie_ref(t_end, m, l, Ja, Jt, ws, th0)
    print("  3. toupie contre l'EDO de Lagrange (axe complet, 1 s, nutation comprise)")
    errs = []
    for h in (4e-4, 2e-4, 1e-4):
        R0 = _expm(np.array([th0, 0.0, 0.0]))
        axe0 = R0 @ np.array([0.0, 0.0, 1.0])
        N = Noyau()
        N.corps("toupie", m, _rot9(np.diag([Jt, Jt, Ja])), list(l * axe0), _rot9(R0), w=list(ws * axe0))
        N.liaison("pointe", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
        tr = N.simule(t_end, h, rho=0.9, tous=max(1, int(round(2e-3 / h))))
        pire = 0.0
        for e in tr:
            t = e[0]
            R = np.array(e[2][0]).reshape(3, 3)
            ax = R[:, 2]
            # l'axe de référence : θ(t) autour de l'axe x initial, φ(t) autour de z
            th = np.interp(t, ts, ths); ph = np.interp(t, ts, phis)
            ax_ref = _expm(np.array([0.0, 0.0, ph])) @ _expm(np.array([th, 0.0, 0.0])) @ np.array([0.0, 0.0, 1.0])
            pire = max(pire, float(np.linalg.norm(ax - ax_ref)))
        errs.append(pire)
        print(f"     h = {h:.0e}  écart max de l'axe {pire:.2e} rad")
    o = [np.log2(errs[i] / errs[i + 1]) for i in range(2)]
    print(f"     ordres {', '.join(f'{x:.2f}' for x in o)}")
    assert errs[-1] < 5e-4, ("toupie : l'axe s'écarte de l'EDO", errs)
    assert min(o) > 1.5, ("toupie : ordre", o)
    return errs


# ── 4 : déterminisme ─────────────────────────────────────────────────────────
def determinisme():
    """Même calcul, RAYON_NUM_THREADS=1 puis tous les fils : identique au bit.

    ⚠ SUR UN CAS OÙ LE PARALLÉLISME S'ACTIVE VRAIMENT. La première version ne
    testait qu'un pendule à deux corps — très en dessous de `PARALLELE_MIN` et
    de `CREUX_PAR = 2000` : elle vérifiait le déterminisme là où rien n'est
    parallèle, ce qui ne prouve rien. On ajoute une chaîne de 7 200 inconnues,
    au-delà du seuil où faer parallélise ses supernœuds, et un robot public.
    """
    import subprocess
    code = ("import numpy as np, json\nfrom vinkulum.verification import _pendule\n"
            "N, tr, _ = _pendule(0.005, 1.0)\nprint(json.dumps([tr[-1][1][0], tr[-1][3][0], tr[-1][5]]))")
    outs = []
    for env in ({"RAYON_NUM_THREADS": "1"}, {}):
        e = dict(os.environ, **env)
        outs.append(subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=e).stdout.strip())
    # le cas PARALLÉLISÉ : 1 200 corps, 7 200 inconnues > CREUX_PAR
    gros = ("import json, numpy as np, vinkulum as V\n"
            "J=[1e-3,0,0,0,1e-3,0,0,0,1e-3]\n"
            "n=V.Noyau([0,0,-9.81]); nb=1200\n"
            "ix=[n.corps(f'c{k}',1.0,J,[0.2*k,0,0]) for k in range(nb)]\n"
            "for k in range(nb-1): n.liaison(f'l{k}',ix[k],ix[k+1],pa=[0.2,0,0],bloque_r=[0,2])\n"
            "n.liaison('b',None,ix[0],pa=[0,0,0],bloque_r=[0,2])\n"
            "n.simule(0.005,1e-3,tous=10**9)\n"
            "p=np.asarray(n.etat()[1]).ravel()\n"
            "print(json.dumps([repr(float(x)) for x in p[-9:]]))")
    gros_outs = []
    for env in ({"RAYON_NUM_THREADS": "1"}, {}):
        e = dict(os.environ, **env)
        gros_outs.append(subprocess.run([sys.executable, "-c", gros], capture_output=True,
                                        text=True, env=e).stdout.strip())
    ok1 = outs[0] == outs[1]
    ok2 = gros_outs[0] == gros_outs[1] and gros_outs[0] != ""
    print(f"  4. déterminisme : pendule {'identique' if ok1 else 'DIFFÉRENT'} · "
          f"chaîne de 7 200 inconnues (LU creux PARALLÉLISÉ) "
          f"{'identique' if ok2 else 'DIFFÉRENT'} — au bit")
    assert ok1, (outs[0][:80], outs[1][:80])
    assert ok2, ("le cas parallelise n'est pas deterministe",
                 gros_outs[0][:80], gros_outs[1][:80])


def gyroscopique_tangent():
    """LE GYROSCOPIQUE N'EST COMPTÉ QU'UNE FOIS dans les tangentes exposées.

    Trouvé le 5 sept. par le pont adjoint : `raideur` retirait ω × J_sω « à la
    main » en plus de celui que `forces()` porte, et `amortissement` ajoutait
    de même [ω]×J_s − [J_sω]× — K_θθ et C_ωω valaient EXACTEMENT 2,000 fois
    l'analytique. Le résidu de Newton n'a jamais doublé (il ne passe pas par
    `forces()`) ; seules les tangentes exposées étaient fausses en rotation.
    Deux contrôles : C_ωω d'un corps libre contre [ω]×J_s − [J_sω]× à 1e-8,
    et la NUTATION d'un corps libre symétrique — `modes_complexes` doit rendre
    ω₃·J_a/J_t, la précession du vecteur vitesse de rotation autour du moment
    cinétique (cinématique d'Euler, aucun paramètre).
    """
    jt, ja, w3 = 1e-3, 2e-3, 30.0
    n = Noyau([0.0, 0.0, 0.0])
    n.corps("t", 1.0, [jt, 0, 0, 0, jt, 0, 0, 0, ja], [0.0] * 3, w=[0.0, 0.0, w3])
    _, c, m, _, _ = (np.asarray(x, float) for x in n.k_c_m_z())
    js, w = m[3:, 3:], np.array([0.0, 0.0, w3])
    sk = lambda v: np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    c_th = sk(w) @ js - sk(js @ w)
    ec = np.max(np.abs(c[3:, 3:] - c_th)) / max(np.max(np.abs(c_th)), 1e-30)
    f, z, sig = n.modes_complexes(combien=1)[0]
    f_th = w3 * ja / jt / (2 * np.pi)
    # ET LA TOUPIE RAPIDE CONTRAINTE — le cas que `modes_complexes` refusait :
    # ses deux modes sont les racines du trinôme Jt'Ω² − Jaω₃Ω + mgl = 0 (Jt'
    # au pivot), précession lente et nutation ; les deux doivent sortir.
    m, l, w3 = 1.0, 0.1, 300.0
    n2 = Noyau([0.0, 0.0, -9.81])
    c = n2.corps("t", m, [jt, 0, 0, 0, jt, 0, 0, 0, ja], [0.0, 0.0, l], w=[0.0, 0.0, w3])
    n2.liaison("piv", None, c, pa=[0.0, 0.0, 0.0], bloque_r=[])
    jtp = jt + m * l * l
    th = np.sort(np.abs(np.roots([jtp, -ja * w3, m * 9.81 * l]))) / (2 * np.pi)
    mc = n2.modes_complexes(combien=4)
    f2 = np.array([x[0] for x in mc])
    s2 = max(abs(x[2]) for x in mc)
    print(f"  29. gyroscopique compté UNE fois : C_ωω à {ec:.1e} de [ω]×J − [Jω]× · "
          f"nutation libre {f:.5f} Hz pour ω₃Ja/Jt = {f_th:.5f} Hz")
    print(f"      toupie rapide CONTRAINTE : {f2.round(6).tolist()} Hz pour les racines exactes "
          f"{th.round(6).tolist()} · |σ| max {s2:.1e}")
    assert ec < 1e-6, ("C_ωω n'est pas la tangente gyroscopique (double compte ?)", ec)
    assert abs(f / f_th - 1.0) < 1e-6 and abs(sig) < 1e-6, ("nutation fausse", f, f_th, sig)
    assert len(f2) == 2 and np.max(np.abs(f2 / th - 1.0)) < 1e-6 and s2 < 1e-6, (f2, th, s2)
    return dict(ec=float(ec), f=float(f), f_th=float(f_th), toupie=f2.tolist())


def audit_jacobien():
    """LE JACOBIEN DE NEWTON CONTRE UNE DIFFÉRENCE FINIE DU RÉSIDU, élément par élément.

    C'est l'instrument qui manquait au noyau, et il a rendu SIX défauts en une
    heure le 5 sept. au soir — tous « quasi-Newton déclarés » ou non déclarés,
    tous mesurés d'abord, puis corrigés :
      · la raideur de POSITION des contacts par pénalité (∂F/∂r) absente
        (7e-4 relatif à h = 1e-3, croît en k·h²) ;
      · ∂F_t/∂λ du frottement non lisse absent (0,22 relatif — le plus gros) ;
      · le terme ∂(G·u)/∂q des lignes GGL et non holonomes laissé de côté
        (5e-2 sur un cardan sous GGL ; srscm passait de 4,6 à 2,9 Newton/pas) ;
      · G de l'engrenage et de la vis « à la main », faux hors de l'axe (7e-4) ;
      · la colonne u̇ des lignes en VITESSE (nh, contact) passée par J_l comme
        une ligne en pose (4e-4) ;
      · ∂Φ_t/∂q d'une cible rhéonome sous GGL (5e-4).
    Le juge : `Noyau.audit_jacobien`, bloc par bloc (dyn | con | ggl) ×
    (u̇ | λ | ζ), écart max relatif au max du bloc, pas de DF calé sur ce que
    chaque inconnue DÉPLACE. Un bloc dont le max est sous 1e-3 du max global
    n'est pas jugé (bruit sur du négligeable). Ce qui reste laissé de côté,
    déclaré : la raideur λ·∂G/∂q du contact non lisse (7,7e-6 mesuré).
    """
    from vinkulum.rotor import polaire_lineaire
    from scipy.spatial.transform import Rotation as R
    J9 = [2e-3, 0, 0, 0, 3e-3, 0, 0, 0, 4e-3]
    rot = lambda v: list(R.from_rotvec(v).as_matrix().ravel())
    g = [0, 0, -9.81]

    def deux(N, sep=0.5):
        a = N.corps("a", 1.0, J9, [0.3, 0.1, 0.2], rot=rot([0.3, -0.2, 0.5]), v=[0.2, -0.1, 0.3], w=[1.0, -2.0, 0.5])
        b = N.corps("b", 0.7, J9, [0.3 + sep, 0.15, 0.1], rot=rot([-0.4, 0.1, 0.2]), v=[-0.1, 0.2, 0.1], w=[0.5, 1.5, -1.0])
        return a, b

    def m_pendules():
        N = Noyau(g); a, b = deux(N)
        N.liaison("p1", None, a, pa=[0, 0, 0], bloque_r=[0, 2]); N.liaison("p2", a, b, pa=[0.25, 0, 0], bloque_r=[1, 2]); return N

    def m_cardan():
        N = Noyau(g); a, b = deux(N)
        N.liaison("rot", None, a, pa=[0, 0, 0], bloque_r=[]); N.liaison("rot2", a, b, pa=[0.25, 0, 0], bloque_r=[])
        N.cardan("c", a, b, [1, 0, 0], [0, 1, 0]); return N

    def m_engrenage():
        N = Noyau(g); a, b = deux(N, 0.3)
        N.liaison("pa", None, a, pa=[0, 0, 0], bloque_r=[0, 1]); N.liaison("pb", None, b, pa=[0, 0, 0], bloque_r=[0, 1])
        N.engrenage("e", a, b, [0, 0, 1], [0, 0, 1], 3.0); return N

    def m_vis():
        N = Noyau(g); a, b = deux(N, 0.1)
        N.liaison("pa", None, a, pa=[0, 0, 0], bloque_r=[0, 1])
        N.liaison("pb", None, b, pa=[0, 0, 0], bloque_t=[0, 1], bloque_r=[0, 1, 2]); N.vis("v", a, b, [0, 0, 1], 0.01); return N

    def m_couples():
        N = Noyau(g); a, b = deux(N)
        N.liaison("pa", None, a, pa=[0, 0, 0], bloque_r=[2]); N.liaison("pb", a, b, pa=[0.25, 0, 0], bloque_r=[0])
        N.couple("k", None, a, [0, 0, 1], ("ressort", [3.0, 0.2, 0.1]))
        N.couple("pd", a, b, [1, 0, 0], ("pd", [2.0, 0.1, 5.0, 50.0]), cible=("table", [0, 0.1, 1, 0.3]))
        N.couple("gov", None, a, [0, 0, 1], ("gouverneur", [0.5, 2.0]), cible=("table", [0, 3.0, 1, 3.0]))
        N.couple("but", a, b, [1, 0, 0], ("butee", [50.0, 1.0, -0.05, 0.05])); return N

    def m_poutres():
        N = Noyau(g); a, b = deux(N, 0.3)
        c = N.corps("c", 0.3, J9, [0.9, 0.1, 0.15], rot=rot([0.1, 0.3, -0.2]), v=[0.1, 0.1, -0.2], w=[0.3, 0.2, 0.1])
        N.liaison("enc", None, a, pa=[0, 0, 0]); N.poutre("p1", a, b, 2e4, 1e4, 3.0, 5.0); N.poutre("p2", b, c, 2e4, 1e4, 3.0, 5.0); return N

    def m_contacts():
        N = Noyau(g); a, b = deux(N, 0.05)
        N.contact("sol", a, [0, 0, -0.1], 0.15, k=1e4, expo=1.5, c=0.3, mu=0.4, v_eps=1e-2)
        N.contact("caps", a, [0.1, 0, 0], 0.15, b=b, pb=[-0.1, 0, 0], rayon_b=0.1, p1=[0.3, 0, 0], p1b=[0.1, 0, 0], k=2e3, mu=0.3, v_eps=1e-2)
        return N

    def m_boite_cyl():
        N = Noyau(g); a, b = deux(N, 0.3)
        N.contact("box", a, [0.3, 0, 0], 0.12, b=b, pb=[0, 0, 0], demi=[0.2, 0.1, 0.1], k=3e3, mu=0.2, v_eps=1e-2)
        N.contact("cyl", a, [0.0, 0.3, 0], 0.15, b=b, pb=[0, 0.1, 0], cylindre=([0, 0, 1], 0.1, 0.2), k=3e3, mu=0.2, v_eps=1e-2)
        N.contact("cable", a, [0, 0, 0], 0.45, b=b, pb=[0, 0, 0], rayon_b=0.0, cable=True, k=5e3); return N

    def m_nonlisse():
        N = Noyau(g); a = N.corps("a", 1.0, J9, [0, 0, 0.1], v=[0.1, 0.0, -0.5], w=[0.3, 0.2, 0.1])
        N.contact("nl", a, [0, 0, 0], 0.1, nonlisse=True, mu=0.3, v_eps=1e-2, restitution=0.0); return N

    def m_pale(**kw):
        N = Noyau([0, 0, 0]); a, b = deux(N)
        N.liaison("piv", None, a, pa=[0, 0, 0], bloque_r=[2]); N.liaison("pb", a, b, pa=[0.25, 0, 0], bloque_r=[1])
        i = N.inflow([0, 0, 1], 0.5, 0.05); N.vent([2.0, 0.5, 0.3])
        N.pale("p", b, [0, 0, 0], [1, 0, 0], [0, 1, 0], 0.4, 0.05, polaire_lineaire(), inflow=i, **kw)
        N.pose_inflow(i, 1.5); return N

    def m_nh():
        N = Noyau(g); a, b = deux(N)
        N.liaison("nh", None, a, pa=[0, 0, 0], bloque_t=[0], bloque_r=[], nh=True); N.liaison("p", a, b, pa=[0.25, 0, 0], bloque_r=[1, 2]); return N

    def m_cible():
        N = Noyau(g); a, b = deux(N)
        N.liaison("mot", None, a, pa=[0, 0, 0], bloque_r=[0, 1, 2], cible_r=([0, 0, 1], ("table", [0, 0, 1, 2.0])))
        N.liaison("p", a, b, pa=[0.25, 0, 0], bloque_r=[1, 2]); return N

    cas = (("pendules", m_pendules, True), ("cardan", m_cardan, True), ("engrenage", m_engrenage, True),
           ("vis", m_vis, False), ("couples", m_couples, False), ("poutres", m_poutres, False),
           ("contacts", m_contacts, False), ("boîte/cylindre/câble", m_boite_cyl, False),
           ("non lisse", m_nonlisse, False), ("pale LB", lambda: m_pale(lb=True), False),
           ("non holonome", m_nh, False), ("cible rhéonome", m_cible, True))
    h = 1e-3
    print("  30. jacobien AD contre différence finie du résidu, par bloc (écart relatif max)")
    pire = 0.0
    for nom, f, avec_ggl in cas:
        for ggl in ((False, True) if avec_ggl else (False,)):
            N = f()
            N.simule(3 * h, h, tous=10 ** 9)
            for sigma in (0., .6, 1.):
                blocs, _ = N.audit_jacobien(h, ggl=ggl, sigma_lie=sigma)
                gmax = max(mx for _, _, mx in blocs)
                juges = [(b, e) for b, e, mx in blocs if mx > 1e-3 * gmax]
                e_max = max(e for _, e in juges)
                pire = max(pire, e_max)
                print(f"     {nom:22} {'GGL ' if ggl else '    '} σ={sigma:g}: " + " ".join(f"{b} {e:.0e}" for b, e in juges))
                assert e_max < 1e-5, (nom, ggl, sigma, juges)
    print(f"     pire écart sur le corpus : {pire:.1e} — un terme laissé de côté s'y lirait")
    return dict(pire=float(pire))


def validation_scores():
    """Une hypothèse qui passe ne doit pas masquer sa voisine qui échoue."""
    from contextlib import redirect_stdout
    from io import StringIO
    from vinkulum import validation
    avant = validation.LIGNES[:]
    try:
        validation.LIGNES.clear()
        fms = {False: {5.0: 0.59, 9.0: 0.8}, True: {5.0: 0.7, 9.0: float("nan")}}
        notes = {False: {5.0: "", 9.0: ""}, True: {5.0: "", 9.0: "ERREUR témoin"}}
        with redirect_stdout(StringIO()):
            validation._maryland_lignes("témoin", 0.59, fms, notes, 24000)
        lignes = validation.LIGNES[:]
        assert len(lignes) == 4, "un N_crit a été sélectionné sur la mesure"
        assert [l[6] for l in lignes] == [True, False, False, False]
        # Changer la mesure change les verdicts, jamais les prédictions publiées.
        validation.LIGNES.clear()
        with redirect_stdout(StringIO()):
            validation._maryland_lignes("témoin", 0.8, fms, notes, 24000)
        np.testing.assert_allclose([l[4] for l in validation.LIGNES], [l[4] for l in lignes], equal_nan=True)
        assert validation.LIGNES[-1][7].startswith("ERREUR")
    finally:
        validation.LIGNES[:] = avant
    print("     validation : toutes les hypothèses N_crit publiées, erreur technique conservée")


def _regressions():
    import unittest
    from .test_noyau import AuditNoyau
    from .test_analyses import AnalysesFiables
    from .test_campagnes import Campagnes
    suite = unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(c)
                               for c in (AuditNoyau, AnalysesFiables, Campagnes))
    resultat = unittest.TextTestRunner(verbosity=1).run(suite)
    assert resultat.wasSuccessful(), "régressions de l'audit du noyau"


def _taches():
    from ._campagnes import tache
    noms = ('convergence', 'toupie', 'determinisme', 'amortissement', 'robustesse',
            '_regressions', 'validation_scores', 'assemblage', 'liaisons_unilaterales',
            'vis_ecrou', 'superelement', 'bennett', 'kapitza', 'statique', 'six_barres',
            'srscm', 'lateral_buckling', 'arbre_tournant', 'invariance_unites',
            'quatre_barres_flexible', 'treillis_de_boucles', 'invariances',
            'modeles_mal_poses', 'echelle', 'non_holonome', 'coulomb_pente', 'accords',
            'lecteurs_casses', 'sparsite', 'refus_gardes', 'gyroscopique_tangent',
            'audit_jacobien')
    cas = [tache(n, 'vinkulum.verification', n, exclusif=n in ('echelle', 'refus_gardes'))
           for n in noms]
    for mod in ('trim', 'floquet', 'reduction', 'adjoint_temps', 'convergence',
                'urdf', 'mjcf', 'domaines', 'campagne'):
        kw = {'rapide': True} if mod in ('floquet', 'reduction', 'adjoint_temps') else {}
        cas.append(tache(mod+'.demo', 'vinkulum.'+mod, 'demo', **kw))
    return cas


def main():
    from ._campagnes import executer
    executer(_taches(), 'vérification')



def assemblage():
    """ASSEMBLER UN MÉCANISME depuis une pose APPROCHÉE — Φ(q) = 0, puis Φ̇ = 0.

    Tout code multicorps a cette fonction (ADAMS l'appelle *position
    analysis*, Simscape *assembly*) et le noyau ne l'avait pas : il fallait
    poser chaque corps à une configuration déjà cohérente. Sur un
    quatre-barres c'est une équation à résoudre à la main ; sur une tête de
    rotor à 71 contraintes, ce n'est pas raisonnable — et le piège des
    VITESSES incohérentes a été payé trois fois dans ce dépôt : le premier
    pas passe à réparer Φ̇, le solveur « ne finit jamais », et rien ne le dit.

    Newton sur la variété avec correction de NORME MINIMALE, δq = −Gᵀ(GGᵀ)⁻¹Φ.
    La norme minimale n'est pas un détail : elle garantit qu'une pose déjà
    cohérente ne bouge pas, et qu'une pose approchée est corrigée par le plus
    petit déplacement — pas par un saut vers l'AUTRE branche du mécanisme (un
    quatre-barres en a deux, et l'utilisateur veut celle qu'il a dessinée).
    """
    ji = [1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4]

    def quatre_barres(redondant=False):
        n = Noyau([0.0, 0.0, -9.81])
        n.corps("a", 1.0, ji, [0.5, 0.0, 0.0])
        n.corps("b", 1.0, ji, [1.1, 0.45, 0.0])
        n.corps("c", 1.0, ji, [0.0, 0.45, 0.0])
        n.liaison("p0", None, 0, pa=[-0.5, 0, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.liaison("p1", 0, 1, pa=[0.5, 0, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.liaison("p2", 1, 2, pa=[0.1, 0.45, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.liaison("p3", None, 2, pa=[0, -0.45, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
        if redondant:                      # une liaison de trop : GGᵀ singulière
            n.liaison("p4", None, 2, pa=[0, -0.45, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1])
        return n

    def deplace(n, d, dv=0.0):
        """Les corps bougent APRÈS la construction — `liaison` déduit la
        contrainte de la pose courante, donc perturber avant ne fait rien."""
        t, r, rot, v, w, vi = n.etat()
        r = [[x + d * (1.0 if k == 0 else -0.5) for k, x in enumerate(p)] for p in r]
        n.pose_etat(t, r, rot, [[dv, 0.0, 0.0] for _ in v], w, vi)

    print("  7. assemblage — Φ(q) = 0 depuis une pose approchée, puis Φ̇ = 0")
    for d in (0.01, 0.05, 0.20):
        n = quatre_barres()
        deplace(n, d, dv=0.3)
        p0 = max(abs(x) for x in n.phi())
        res, it = n.assemble()
        pd = max(abs(x) for x in n.phi_dot())
        print(f"     écart {d:.2f} m : |Φ| {p0:.2e} → {res:.1e} en {it} it · |Φ̇| 3.0e-01 → {pd:.1e}")
        assert res < 1e-11 and pd < 1e-13, (d, res, pd)
        assert it <= 8, ("l'assemblage traîne", d, it)

    # · NORME MINIMALE : une pose cohérente ne bouge pas — zéro itération.
    n = quatre_barres()
    res, it = n.assemble()
    assert it == 0 and res < 1e-14, ("une pose cohérente a bougé", res, it)

    # · REDONDANT : GGᵀ singulière, le repli SVD rend la norme minimale.
    n = quatre_barres(redondant=True)
    deplace(n, 0.05)
    res, it = n.assemble()
    assert res < 1e-12, ("l'assemblage échoue sur des contraintes redondantes", res)
    print(f"     redondant (GGᵀ singulière) : |Φ| → {res:.1e} par repli SVD")

    # · IMPOSSIBLE : deux corps encastrés à 10 m, une distance imposée de 1.
    #   Doit ÉCHOUER, pas boucler — c'est le pendant du `h = 0` de `robustesse`.
    m = Noyau([0.0, 0.0, 0.0])
    m.corps("a", 1.0, ji, [0.0, 0.0, 0.0])
    m.corps("b", 1.0, ji, [10.0, 0.0, 0.0])
    m.liaison("ga", None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    m.liaison("gb", None, 1, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    m.distance("d", 0, 1, [0, 0, 0], [0, 0, 0], l=1.0)
    try:
        m.assemble(iters=40)
        raise AssertionError("un mécanisme non assemblable a été « assemblé »")
    except ValueError:
        pass
    print("     mécanisme non assemblable : refusé avec son ‖Φ‖, pas de boucle")
    return dict(ok=True)


def liaisons_unilaterales():
    """BUTÉE et CÂBLE — les deux liaisons unilatérales que tout mécanisme a.

    Une butée d'articulation (vérin en fin de course, palonnier contre son
    arrêt, genou) et un câble (qui tire mais ne pousse jamais) manquaient au
    noyau. Les deux sont UNILATÉRAUX, donc ni une liaison ni une force
    ordinaire : ils n'agissent que d'un côté.

    · **Butée** : loi de couple `("butee", [k, c, θ_min, θ_max])`, nulle dans
      la plage. L'amortissement suit la PÉNÉTRATION, pas ω seul — sinon un
      amortisseur linéaire TIRE au décollement et arrache la pièce de sa
      butée. C'est la faute déjà payée sur le contact (Hunt–Crossley).
    · **Câble** : le contact au signe près. Un contact repousse quand deux
      corps se rapprochent, un câble tire quand ils s'éloignent — même
      élément, `cable=True`, et le « rayon » devient la longueur au repos.
      Il donne aussi le ressort unilatéral et la butée de TRANSLATION.

    Le juge de la butée est √(2E/k) : un ressort unilatéral laisse dépasser
    d'autant que l'énergie qu'il doit absorber, et c'est vérifiable sur
    plusieurs décades de raideur — un « mur » ne le serait pas.
    """
    jz = 2e-3
    print("  8. liaisons unilatérales — butée d'articulation et câble")
    ecarts = []
    for k in (20.0, 200.0, 2000.0, 20000.0):
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("d", 1.0, [jz, 0, 0, 0, jz, 0, 0, 0, jz], [0.0] * 3, w=[0.0, 0.0, 2.0])
        n.liaison("piv", None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.couple("but", None, 0, [0, 0, 1], ("butee", [k, 0.0, radians(-15.0), radians(15.0)]))
        mx = 0.0
        for j in range(2000):
            n.simule((j + 1) * 1e-4, 1e-4, tous=10 ** 9)
            mx = max(mx, abs(n.couples()[0][1]))
        th = sqrt(2.0 * (0.5 * jz * 4.0) / k)
        ecarts.append(abs((mx - radians(15.0)) / th - 1.0))
    print(f"     butée : le dépassement suit √(2E/k) sur 4 décades de raideur — "
          f"pire écart {100 * max(ecarts):.2f} %")
    assert max(ecarts) < 0.05, ("la butée n'est pas un ressort unilatéral", ecarts)
    # · contrôle NÉGATIF : sans butée le disque tourne librement.
    n = Noyau([0.0, 0.0, 0.0])
    n.corps("d", 1.0, [jz, 0, 0, 0, jz, 0, 0, 0, jz], [0.0] * 3, w=[0.0, 0.0, 2.0])
    n.liaison("piv", None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1])
    n.couple("nul", None, 0, [0, 0, 1], ("constant", [0.0]))
    n.simule(0.2, 1e-4, tous=10 ** 9)
    assert abs(n.couples()[0][1] - 0.4) < 1e-6, ("le contrôle négatif est faux", n.couples())

    # · CÂBLE : une masse suspendue tombe, puis EST RETENUE à la longueur.
    ji = [1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4]
    fin = {}
    for cab in (False, True):
        n = Noyau([0.0, 0.0, -9.81])
        n.corps("m", 1.0, ji, [0.0, 0.0, -0.1])
        n.corps("anc", 1.0, ji, [0.0] * 3)
        n.liaison("f", None, 1, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
        n.contact("c", 0, [0.0] * 3, 0.25, b=1, pb=[0.0] * 3, rayon_b=0.25,
                  k=1e4, c=2.0, expo=1.0, cable=cab)
        tr = n.simule(2.0, 1e-4, tous=200)
        fin[cab] = min(f[1][0][2] for f in tr)
    print(f"     câble : sans lui la masse tombe à {fin[False]:+.2f} m, avec lui "
          f"elle est retenue à {fin[True]:+.4f} m pour 0,50 de longueur")
    assert fin[False] < -15.0, ("le contrôle négatif ne tombe pas", fin)
    assert -0.53 < fin[True] < -0.50, ("le câble ne retient pas à sa longueur", fin)
    return dict(butee=float(max(ecarts)), cable=fin[True])


def vis_ecrou():
    """VIS–ÉCROU (glissière hélicoïdale) — translation liée à la rotation.

    Un vérin à vis, une table de machine-outil, un actionneur linéaire
    électrique en sont faits, et aucun ne se modélise par une glissière plus
    un pivot : il faut la CONTRAINTE qui lie les deux.

        Φ = n·(r_b − r_a) − pas·θ − c₀ ,  θ déroulé, pas en m/rad

    LE JUGE N'EST PAS LA CINÉMATIQUE mais l'INERTIE RAMENÉE. Vérifier que
    « N tours avancent de N × pas » ne teste que Φ, que le solveur tient de
    toute façon à 1e-17. Ce qui prouve que la liaison transmet la DYNAMIQUE
    est que l'écrou freine la vis exactement comme une inertie de plus :

        (I + m·(pas/2π)²)·ω = I·ω₀        (aucun effort extérieur)

    C'est une conservation, elle ne dépend d'aucun réglage, et elle croît avec
    le pas — donc elle discrimine sur deux décades.
    """
    iz, m = 1e-3, 1.0
    print("  9. vis–écrou — la contrainte, et l'inertie qu'elle ramène")
    ecarts, ecarts_v, phimax = [], [], 0.0
    for pas in (0.002, 0.010, 0.050, 0.100):
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("e", m, [iz, 0, 0, 0, iz, 0, 0, 0, iz], [0.0] * 3)
        n.liaison("gl", None, 0, bloque_t=[0, 1], bloque_r=[0, 1, 2])
        n.corps("v", 1.0, [iz, 0, 0, 0, iz, 0, 0, 0, iz], [0.0, 0.0, 0.2], w=[0.0, 0.0, 10.0])
        n.liaison("pv", None, 1, bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.vis("vis", 0, 1, [0, 0, 1], pas)
        n.simule(0.02, 1e-5, tous=10 ** 9)
        e = n.etat()
        w, vz = e[4][1][2], e[3][0][2]
        pr = pas / (2.0 * pi)
        w_th = iz * 10.0 / (iz + m * pr * pr)
        ecarts.append(abs(w / w_th - 1.0))
        ecarts_v.append(abs(vz + pr * w))
        phimax = max(phimax, max(abs(x) for x in n.phi()))
    print(f"     inertie ramenée (I + m(pas/2π)²)ω = Iω₀ sur 2 décades de pas — "
          f"pire écart {100 * max(ecarts):.4f} %")
    print(f"     |Φ| max {phimax:.1e} m · |v_z − pas·ω/2π| max {max(ecarts_v) * 1e9:.0f} nm/s")
    assert max(ecarts) < 1e-5, ("la vis ne ramène pas la bonne inertie", ecarts)
    assert phimax < 1e-14 and max(ecarts_v) < 1e-6, (phimax, ecarts_v)
    # · un pas NUL dégénérerait en glissière : refusé.
    try:
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("a", 1.0, [iz, 0, 0, 0, iz, 0, 0, 0, iz], [0.0] * 3)
        n.corps("b", 1.0, [iz, 0, 0, 0, iz, 0, 0, 0, iz], [0.0, 0.0, 0.1])
        n.vis("v", 0, 1, [0, 0, 1], 0.0)
        raise AssertionError("une vis de pas nul a été acceptée")
    except ValueError:
        pass
    return dict(inertie=float(max(ecarts)), phi=float(phimax))


def superelement():
    """SUPERÉLÉMENT — la raideur d'un maillage EF quelconque, dans le multicorps.

    C'est ainsi qu'un noyau multicorps devient GÉNÉRALISTE côté flexible : on
    ne maille pas une coque ou un solide dans le solveur (ADAMS et RecurDyn ne
    le font pas non plus), on prend le K d'un code EF, on le réduit —
    Craig–Bampton, que ce dépôt a déjà — et on l'attache à des nœuds qui sont
    des corps ordinaires. Liaisons, contact, modes, adjoints suivent.

    DEUX CONTRÔLES, et le second est celui qui peut réfuter :

    1. **le pont** : un K de poutre analytique doit rendre la même dynamique
       que les poutres géométriquement exactes du noyau, elles-mêmes validées
       contre Euler–Bernoulli. Deux discrétisations différentes convergent
       vers la même physique ;
    2. **le corotationnel** : une ROTATION RIGIDE ne doit produire AUCUNE
       force. Un K constant en repère global en produirait — il fabriquerait
       du rappel là où il n'y a pas de déformation, et rien dans le contrôle 1
       ne le montrerait (la console ne tourne pas).

    ⚠ Non couvert, déclaré : la raideur GÉOMÉTRIQUE (précontrainte,
    raidissement centrifuge). Pour ça les poutres GE restent le bon outil —
    c'est pourquoi elles ont été faites avant, et pourquoi Craig–Bampton seul
    manque la raideur centrifuge d'une pale.
    """
    e_mod, g_mod = 2.1e11, 8.1e10
    b, h = 0.02, 0.004
    a_s, iy, iz = b * h, b * h ** 3 / 12, h * b ** 3 / 12
    jt, rho, lg, ne = iy + iz, 7800.0, 0.5, 8
    le = lg / ne

    def k_poutre(l):
        k = np.zeros((12, 12))
        ea = e_mod * a_s / l
        k[0, 0] = k[6, 6] = ea
        k[0, 6] = k[6, 0] = -ea
        gj = g_mod * jt / l
        k[3, 3] = k[9, 9] = gj
        k[3, 9] = k[9, 3] = -gj
        for i1, i2, r1, r2, ii, sg in ((1, 7, 5, 11, iz, 1.0), (2, 8, 4, 10, iy, -1.0)):
            aa, bb = 12 * e_mod * ii / l ** 3, 6 * e_mod * ii / l ** 2
            cc, dd = 4 * e_mod * ii / l, 2 * e_mod * ii / l
            k[i1, i1] += aa; k[i2, i2] += aa; k[i1, i2] -= aa; k[i2, i1] -= aa
            k[r1, r1] += cc; k[r2, r2] += cc; k[r1, r2] += dd; k[r2, r1] += dd
            k[i1, r1] += sg * bb; k[r1, i1] += sg * bb
            k[i1, r2] += sg * bb; k[r2, i1] += sg * bb
            k[i2, r1] -= sg * bb; k[r1, i2] -= sg * bb
            k[i2, r2] -= sg * bb; k[r2, i2] -= sg * bb
        return k

    def console(par_superelement):
        n = Noyau([0.0, 0.0, 0.0])
        mn = rho * a_s * le
        idx = [n.corps(f"n{i}", mn if 0 < i < ne else mn / 2,
                       [mn * le * le / 12, 0, 0, 0, mn * le * le / 12, 0, 0, 0, mn * le * le / 12],
                       [i * le, 0.0, 0.0]) for i in range(ne + 1)]
        n.liaison("enc", None, idx[0], bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
        if par_superelement:
            k = k_poutre(le)
            for i in range(ne):
                n.superelement(f"s{i}", [idx[i], idx[i + 1]], list(k.ravel()))
        else:
            for i in range(ne):
                n.poutre(f"e{i}", idx[i], idx[i + 1], e_mod * a_s, g_mod * a_s,
                         g_mod * jt, e_mod * iy, ei3=e_mod * iz)
        return n

    f_ge = [x for x, _ in console(False).modes(combien=2)][:2]
    f_se = [x for x, _ in console(True).modes(combien=2)][:2]
    f_th = 1.875 ** 2 / (2 * pi * lg * lg) * sqrt(e_mod * iy / (rho * a_s))
    print("  10. superélément — la raideur d'un maillage EF, dans le multicorps")
    print(f"     console 8 él. : poutre GE {f_ge[0]:.4f} / {f_ge[1]:.2f} Hz · "
          f"superélément {f_se[0]:.4f} / {f_se[1]:.2f} Hz · théorie {f_th:.4f} Hz")
    assert abs(f_se[0] / f_ge[0] - 1) < 0.01, ("le K importé ne rend pas la dynamique de l'élément natif", f_ge, f_se)
    assert abs(f_se[0] / f_th - 1) < 0.02, ("ni celle de la théorie", f_se, f_th)

    # ── le COROTATIONNEL : une rotation rigide ne produit rien ────────────
    k = np.zeros((12, 12))
    for i in range(12):
        k[i, i] = 1e6
    for i, j in ((0, 6), (1, 7), (2, 8)):
        k[i, j] = k[j, i] = -1e6
    ji = [1e-5, 0, 0, 0, 1e-5, 0, 0, 0, 1e-5]
    pire = 0.0
    for ang in (15.0, 90.0, 180.0):
        c, s2 = cos(radians(ang)), sin(radians(ang))
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("a", 0.1, ji, [0.0] * 3)
        n.corps("b", 0.1, ji, [0.1, 0.0, 0.0])
        n.superelement("se", [0, 1], list(k.ravel()))
        t, r, _, v, w, vi = n.etat()
        rm = np.array([[c, -s2, 0.0], [s2, c, 0.0], [0.0, 0.0, 1.0]])
        n.pose_etat(t, [list(rm @ np.asarray(p)) for p in r], [list(rm.ravel())] * 2, v, w, vi)
        n.simule(1e-7, 1e-7, tous=10 ** 9)
        pire = max(pire, max(np.linalg.norm(n.etat()[3][i]) for i in (0, 1)))
    print(f"     corotationnel : rotation rigide jusqu'à 180° → vitesse induite {pire:.1e} m/s")
    assert pire < 1e-12, ("une rotation rigide produit une force : K n'est pas corotationnel", pire)
    # · et une VRAIE déformation en produit une, sinon le contrôle ne dit rien
    n = Noyau([0.0, 0.0, 0.0])
    n.corps("a", 0.1, ji, [0.0] * 3)
    n.corps("b", 0.1, ji, [0.1, 0.0, 0.0])
    n.superelement("se", [0, 1], list(k.ravel()))
    t, _, rot, v, w, vi = n.etat()
    n.pose_etat(t, [[0.0] * 3, [0.101, 0.0, 0.0]], rot, v, w, vi)
    n.simule(1e-5, 1e-6, tous=10 ** 9)
    assert np.linalg.norm(n.etat()[3][1]) > 1e-3, "une déformation réelle ne produit pas de force"

    # · K non symétrique refusée : une raideur qui ne dérive d'aucune énergie
    try:
        kb = k.copy()
        kb[0, 3] += 1.0
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("a", 0.1, ji, [0.0] * 3)
        n.corps("b", 0.1, ji, [0.1, 0.0, 0.0])
        n.superelement("se", [0, 1], list(kb.ravel()))
        raise AssertionError("un K non symétrique a été accepté")
    except ValueError:
        pass
    print("     K non symétrique refusée · rotation rigide neutre · déformation active")
    return dict(f_ge=f_ge, f_se=f_se, corot=float(pire))


def bennett():
    """LE MÉCANISME DE BENNETT (1903) — un modèle que l'auteur du code n'a PAS
    conçu, et dont le comportement contredit le calcul naïf.

    Tous les autres bancs de ce noyau ont été écrits par celui qui a écrit le
    solveur : c'est le biais qui reste quand la liste de fonctionnalités est
    remplie. On ne le lève pas en ajoutant une brique — on le lève en prenant
    un problème POSÉ AILLEURS, dont la réponse est publiée et qu'on ne peut
    pas faire passer par construction.

    Bennett est ce problème. Un 4R spatial : trois corps mobiles (18 ddl),
    quatre pivots (20 contraintes). **Grübler prédit une mobilité de −2**,
    c'est-à-dire une structure RIGIDE. Le mécanisme bouge quand même — il a
    une mobilité de 1 — pourvu que ses dimensions vérifient

        a / sin α  =  b / sin β        (liens opposés égaux)

    Trois contraintes y sont donc REDONDANTES, et un solveur qui compte ses
    contraintes au lieu d'en mesurer le RANG déclare le mécanisme bloqué.

    LE JUGE est direct et n'implique aucune dynamique : on projette une
    vitesse ARBITRAIRE sur le noyau des contraintes (`assemble`). Elle ne
    survit que si la mobilité est ≥ 1.

    ET LE CONTRE-CONTRÔLE EST LA MOITIÉ DU TEST : avec la condition de Bennett
    VIOLÉE, la même construction doit rendre un mécanisme rigide — vitesse
    projetée à zéro, et la boucle qui ne se ferme plus. Sans lui, « ça bouge »
    ne prouverait rien : un solveur trop permissif bougerait aussi.
    """
    from scipy.optimize import least_squares

    def dh(th, a, al):
        ct, st, ca, sa = cos(th), sin(th), cos(al), sin(al)
        return np.array([[ct, -st * ca, st * sa, a * ct],
                         [st, ct * ca, -ct * sa, a * st],
                         [0.0, sa, ca, 0.0], [0.0, 0.0, 0.0, 1.0]])

    def log3(r):
        c = np.clip((np.trace(r) - 1.0) / 2.0, -1.0, 1.0)
        th = np.arccos(c)
        v = np.array([r[2, 1] - r[1, 2], r[0, 2] - r[2, 0], r[1, 0] - r[0, 1]])
        return 0.5 * v if abs(th) < 1e-12 else th / (2.0 * sin(th)) * v

    def poses(fac, th1=radians(80.0)):
        al, be, a = radians(50.0), radians(80.0), 0.30
        b = a * sin(be) / sin(al) * fac
        par = [(a, al), (b, be), (a, al), (b, be)]

        def r6(x):
            t = np.eye(4)
            for th, (aa, alp) in zip([th1] + list(x), par):
                t = t @ dh(th, aa, alp)
            return np.concatenate([t[:3, 3], log3(t[:3, :3])])

        best = None
        for g in ([2.0, -1.0, 2.0], [1.0, 1.0, 1.0], [2.5, -2.0, 2.5]):
            r = least_squares(r6, g, xtol=1e-15, ftol=1e-15)
            if best is None or r.cost < best.cost:
                best = r
        ths, ts = [th1] + list(best.x), [np.eye(4)]
        for th, (aa, alp) in zip(ths, par):
            ts.append(ts[-1] @ dh(th, aa, alp))
        return ts, float(np.linalg.norm(best.fun))

    ji = [1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4]

    def monte(fac):
        ts, res = poses(fac)
        rng = np.random.default_rng(5)
        n = Noyau([0.0, 0.0, 0.0])
        for i in (1, 2, 3):
            n.corps(f"l{i}", 0.5, ji, list(ts[i][:3, 3]))
        t0, _, _, _, _, vi = n.etat()
        n.pose_etat(t0, [list(ts[i][:3, 3]) for i in (1, 2, 3)],
                    [list(ts[i][:3, :3].ravel()) for i in (1, 2, 3)],
                    [list(rng.normal(size=3)) for _ in range(3)],
                    [list(rng.normal(size=3)) for _ in range(3)], vi)
        pa_ir = [(None, 0), (0, 1), (1, 2), (2, None)]
        for k in range(4):
            a_, b_ = pa_ir[k]
            rk, pk = ts[k][:3, :3], ts[k][:3, 3]
            if a_ is None:
                pa, ra = list(pk), list(rk.ravel())
            else:
                ra_m, pa_v = ts[a_ + 1][:3, :3], ts[a_ + 1][:3, 3]
                pa, ra = list(ra_m.T @ (pk - pa_v)), list((ra_m.T @ rk).ravel())
            n.liaison(f"p{k + 1}", a_, b_, pa=pa, ra=ra, bloque_t=[0, 1, 2], bloque_r=[0, 1])
        return n, res

    print("  11. Bennett (1903) — un mécanisme POSÉ AILLEURS, que Grübler déclare rigide")
    out = {}
    for fac, nom in ((1.0, "condition vérifiée"), (1.3, "condition violée  ")):
        n, res = monte(fac)
        v0 = max(np.linalg.norm(n.etat()[3][i]) for i in range(3))
        n.assemble()
        v1 = max(np.linalg.norm(n.etat()[3][i]) for i in range(3))
        w1 = max(np.linalg.norm(n.etat()[4][i]) for i in range(3))
        pd = max(abs(x) for x in n.phi_dot())
        out[fac] = (res, v1, w1, pd)
        print(f"     {nom} : fermeture {res:.1e} · vitesse arbitraire {v0:.3f} → "
              f"{v1:.6f} m/s (ω {w1:.6f}) · |Φ̇| {pd:.1e}")
    # · Bennett BOUGE : 18 ddl, 20 contraintes, et pourtant mobilité ≥ 1.
    assert out[1.0][1] > 1e-3, ("le Bennett est déclaré rigide : la redondance "
                                "n'est pas détectée, le solveur COMPTE ses contraintes "
                                "au lieu d'en mesurer le rang", out)
    # · et la condition violée donne bien un mécanisme RIGIDE — sans ce
    #   contre-contrôle, « ça bouge » ne prouverait rien.
    assert out[1.3][1] < 1e-9, ("un mécanisme qui ne se ferme pas bouge quand même : "
                                "le solveur est trop permissif", out)
    assert out[1.0][0] < 1e-12 < out[1.3][0], ("la construction ne discrimine pas", out)
    assert max(out[f][3] for f in out) < 1e-12, ("Φ̇ n'est pas tenu", out)
    print("     ⇒ mobilité 1 retrouvée là où le compte des contraintes dit −2,")
    print("       et rigidité retrouvée quand la condition de Bennett tombe.")
    return {str(k): v for k, v in out.items()}


def kapitza():
    """LE PENDULE DE KAPITZA (1951) — second problème posé ailleurs, d'une
    nature tout autre que Bennett : non plus cinématique, mais STABILITÉ.

    Un pendule INVERSÉ tombe. Faites vibrer sa base assez vite, et il tient
    debout — c'est le résultat de Kapitza, et il est contre-intuitif au point
    qu'on ne peut pas le faire sortir d'un solveur par hasard. Le critère est
    analytique :

        a²ω² > 2 g L        (a amplitude de la base, ω sa pulsation)

    Rien de ce nombre n'est écrit dans le noyau : il doit ÉMERGER de la
    simulation. Et le contrôle discrimine des deux côtés — sous le seuil le
    pendule doit TOMBER, ce qui interdit à un solveur trop rigide (ou trop
    amorti) de « réussir » en bloquant tout.
    """
    lg, mp, g = 0.20, 0.5, 9.81
    seuil = 2.0 * g * lg

    def essai(a, f_hz, t_end=3.0, th0=radians(12.0)):
        w = 2.0 * pi * f_hz
        n = Noyau([0.0, 0.0, -g])
        base = n.corps("base", 1.0, [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3], [0.0] * 3)
        tab = [(t, a * cos(w * t)) for t in np.linspace(0.0, t_end, 4000)]
        n.liaison("vert", None, base, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2],
                  cible_t=([0, 0, 1], ("table", [x for pr in tab for x in pr])))
        pend = n.corps("p", mp, [1e-5, 0, 0, 0, 1e-5, 0, 0, 0, 1e-5],
                       [lg * sin(th0), 0.0, lg * cos(th0)])
        n.liaison("piv", base, pend, pa=[0.0] * 3, bloque_t=[0, 1, 2], bloque_r=[0, 2])
        tr = n.simule(t_end, 2e-5, tous=200)
        return max(abs(degrees(np.arctan2(f[1][1][0], f[1][1][2] - f[1][0][2]))) for f in tr)

    print(f"  12. Kapitza (1951) — un pendule INVERSÉ que la vibration stabilise "
          f"(critère a²ω² > 2gL = {seuil:.2f})")
    out = {}
    for a, f in ((0.005, 20), (0.010, 30), (0.015, 50), (0.020, 80)):
        w = 2.0 * pi * f
        r = (a * w) ** 2 / seuil
        mx = essai(a, f)
        out[(a, f)] = (r, mx)
        print(f"     a {a * 1e3:4.0f} mm · f {f:3d} Hz : a²ω²/(2gL) = {r:6.2f} → "
              f"|θ| max {mx:6.1f}°   {'TIENT DEBOUT' if mx < 40 else 'tombe'}")
    sous = [v for k, v in out.items() if v[0] < 1.0]
    sur = [v for k, v in out.items() if v[0] > 5.0]
    assert sous and all(v[1] > 90.0 for v in sous), ("sous le critère, le pendule devrait TOMBER — "
                                                     "un solveur qui le tient debout est trop rigide", out)
    assert sur and all(v[1] < 40.0 for v in sur), ("au-dessus du critère, le pendule devrait tenir "
                                                   "debout : la stabilisation paramétrique n'est pas reproduite", out)
    print("     ⇒ la bascule se fait AU CRITÈRE, et ce nombre n'est écrit nulle part dans le noyau.")
    return {f"{k[0]}_{k[1]}": v for k, v in out.items()}


def statique():
    """ANALYSE STATIQUE et POUTRE DE PRINCETON — le benchmark expérimental.

    Deux choses d'un coup, parce que la seconde a exigé la première.

    **L'analyse statique.** Tout code multicorps en a une (ADAMS
    *equilibrium*, MBDyn `initial assembly`, Simscape *steady state*) ;
    vinkulum atteignait l'équilibre en amortissant une trajectoire. On résout
    désormais le vrai système par Newton — f(q) + Gᵀλ = 0 et Φ(q) = 0,
    jacobien [[K, Gᵀ], [G, 0]] — avec continuation par paliers de charge
    (`pose_effort`), ce que fait tout code non linéaire.

    **La poutre de Princeton** (Dowell & Traybar 1975) est le juge : une
    poutre en porte-à-faux chargée en bout, à sept orientations de charge, du
    plan fort (EJY = 36,3) au plan faible (EJZ = 2,43). C'est LE banc de
    référence des poutres géométriquement exactes, et le seul des huit
    benchmarks livrés avec MBDyn dont la référence soit **expérimentale**.

    QUATRE HYPOTHÈSES ONT ÉTÉ PAYÉES avant que le cas passe, et les trois
    premières étaient fausses :
    · « c'est le bruit de relaxation » — le statique ne relaxe pas ;
    · « c'est le conditionnement » (six ordres entre EA et GJ) — l'équilibrage
      de Jacobi ne déplace pas le plafond d'une décimale ;
    · « un seul module de cisaillement pour une section qui en a deux » — GAY
      et GAZ donnent le même résultat ;
    · **la bonne** : le banc n'amortissait qu'autour de **z**. Un mouvement
      hors de ce plan n'était donc pas amorti du tout, et la relaxation
      s'arrêtait sur un état qui n'était PAS un équilibre — résidu 0,2 à
      0,7 N. Sur les trois axes, tout tombe à 3e-7 et les sept angles passent.

    *Un banc qui n'amortit qu'un axe ne relaxe qu'un plan* — et il rend des
    chiffres d'apparence plausible partout ailleurs. Le seul contrôle qui l'a
    dit est le RÉSIDU, pas l'écart à la référence.
    """
    lg, ea, ejy, ejz, gj, gay = 0.508, 2.84191e6, 3.62794e1, 2.42873e0, 3.10338e0, 6.40131e5
    p_charge, ne = 8.896, 10
    # MBDyn sur le même modèle (jeu de tests livré), (u₂, u₃) du bout par θ
    ref = {0: (0.0, 0.01072), 15: (0.04213, 0.01063), 30: (0.07946, 0.01015),
           45: (0.10904, 0.00889), 60: (0.12992, 0.00665), 75: (0.14222, 0.00357),
           90: (0.14626, 0.0)}
    le = lg / ne
    mn = 2700.0 * 3.175e-3 * 12.7e-3 * le

    def poutre_chargee(theta, amorti_3d=True):
        n = Noyau([0.0, 0.0, 0.0])
        idx = [n.corps(f"n{i}", max(mn if 0 < i < ne else mn / 2, 1e-9),
                       [1e-9, 0, 0, 0, mn * le * le / 12, 0, 0, 0, mn * le * le / 12],
                       [i * le, 0.0, 0.0]) for i in range(ne + 1)]
        n.liaison("enc", None, idx[0], bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
        for i in range(ne):
            n.poutre(f"e{i}", idx[i], idx[i + 1], ea, gay, gj, ejy, ei3=ejz)
        th = radians(theta)
        n.effort(idx[ne], [0.0, p_charge * sin(th), p_charge * cos(th)], [0.0, 0.0, 0.0])
        axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] if amorti_3d else [[0, 0, 1]]
        for i in idx:
            for k, ax in enumerate(axes):
                n.couple(f"a{i}_{k}", None, i, ax, ("ressort", [0.0, 5e-3, 0.0]))
        n.simule(8.0, 5e-4, tous=10 ** 9)
        # GARDE du Newton modifié inter-pas (7 sept.) : sur cette poutre RELAXÉE un
        # jacobien sert ~40 pas (384 pour 16 000, mesuré). S'il redevenait un par pas,
        # ce contrôle repasserait de 3 à 10 s sans qu'aucun assert de physique ne le
        # dise. Pas sur le cas à un seul axe : il n'atteint pas l'équilibre, son état
        # bouge à chaque pas et le jacobien avec (0,91 par pas, mesuré — c'est juste).
        if amorti_3d:
            assert n.stats()[2] < 0.1 * round(8.0 / 5e-4), ("jacobien recalculé à chaque pas", n.stats())
        e = n.etat()
        res = max(abs(x) for x in np.asarray(n.residu_statique()).reshape(-1, 6)[ne])
        return e[1][idx[ne]][1], e[1][idx[ne]][2], res

    print("  13. Princeton (Dowell & Traybar 1975) — poutre chargée, contre MBDyn")
    ecarts, residus = [], []
    for th in sorted(ref):
        u2, u3, res = poutre_chargee(th)
        r2, r3 = ref[th]
        d = float(np.hypot(u2 - r2, u3 - r3))
        ecarts.append(d)
        residus.append(res)
        print(f"     θ {th:3d}° : MBDyn ({r2:+.5f},{r3:+.5f}) · vinkulum ({u2:+.5f},{u3:+.5f})"
              f" · écart {100 * d / 0.146:5.2f} % · résidu {res:.1e} N")
    print(f"     les 7 angles convergés (résidu < {max(residus):.0e} N), "
          f"écart max {100 * max(ecarts) / 0.146:.2f} % du déplacement principal")
    assert max(residus) < 1e-5, ("un de ces états n'est pas un équilibre : ne rien conclure", residus)
    assert max(ecarts) < 0.146 * 0.01, ("vinkulum s'écarte de MBDyn sur Princeton", ecarts)

    # · LE CONTRÔLE QUI A TOUT DÉBLOQUÉ, gardé : n'amortir qu'un axe laisse
    #   des états NON équilibrés, avec des déplacements d'apparence plausible.
    _, _, res_1ax = poutre_chargee(45, amorti_3d=False)
    assert res_1ax > 1e-2, ("amortir un seul axe ne laisse plus de résidu : "
                            "ce contrôle ne prouve plus rien, le rejouer", res_1ax)
    print(f"     (à un seul axe amorti, le même cas s'arrête à {res_1ax:.2f} N de résidu :")
    print(f"      un banc qui n'amortit qu'un axe ne relaxe qu'un plan.)")

    # · l'analyse STATIQUE sur le même cas, sans relaxation du tout
    n = Noyau([0.0, 0.0, 0.0])
    idx = [n.corps(f"s{i}", max(mn if 0 < i < ne else mn / 2, 1e-9),
                   [1e-9, 0, 0, 0, mn * le * le / 12, 0, 0, 0, mn * le * le / 12],
                   [i * le, 0.0, 0.0]) for i in range(ne + 1)]
    n.liaison("enc", None, idx[0], bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    for i in range(ne):
        n.poutre(f"e{i}", idx[i], idx[i + 1], ea, gay, gj, ejy, ei3=ejz)
    n.effort(idx[ne], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    tot = 0
    for k in range(1, 11):
        n.pose_effort(0, [0.0, 0.0, p_charge * k / 10.0], [0.0, 0.0, 0.0])
        r, it = n.statique(tol=1e-8, iters=60)
        tot += it
    u3 = n.etat()[1][idx[ne]][2]
    ec = abs(u3 / ref[0][1] - 1.0)
    print(f"     et par analyse STATIQUE (θ=0, 10 paliers, {tot} it) : {100 * ec:.3f} % — "
          f"sans avoir à choisir un temps d'arrêt")
    assert ec < 0.01 and r < 1e-8, ("la statique s'écarte de MBDyn", u3, r)
    return dict(ecarts=[float(x) for x in ecarts], residus=[float(x) for x in residus],
                statique=float(ec))


def six_barres():
    """6barmech — le mécanisme à six barres du jeu de tests de MBDyn.

    Cinquième problème posé ailleurs, et le plus simple à énoncer : cinq
    corps, sept liaisons, contraintes choisies COMPOSANTE PAR COMPOSANTE par
    ses auteurs (29 pour 30 ddl, mobilité 1), lâché en rotation sous gravité.
    Rien n'est de moi : ni la géométrie, ni les inerties, ni le choix des
    contraintes, ni la référence.

    LA LEÇON DE MÉTHODE, et elle vaut pour toute comparaison de solveurs :
    **on ne compare pas un solveur RAFFINÉ à un solveur qui ne l'est pas.**
    À h = 8e-3 des deux côtés, l'écart valait 1,7e-3 ; en raffinant vinkulum
    seul il MONTAIT à 6,0e-3, ce qui ressemblait à une erreur de modèle. En
    raffinant MBDyn aussi, on voit qu'il se déplace exactement d'autant — et
    que les deux convergent vers le MÊME point à 2,8e-6. Comparer à pas
    différents mélange l'erreur de modèle et celle d'intégration, et fait
    accuser le modèle.
    """
    lg, mm, om = 1.0, 1.0, 1.0
    j1, j2 = mm * lg * lg / 12.0, 1e-9
    # MBDyn au pas fin (h = 5e-4) : position du corps 101 à t = 3 s
    ref = np.array([-0.049323, -0.835610])

    def run(h, t_end=3.0):
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
        n.simule(t_end, h, tous=10 ** 9)
        p = n.etat()[1][idx[101]]
        return np.array([p[0], p[2]]), phi0, phid0, max(abs(x) for x in n.phi())

    print("  14. 6barmech — 5 corps, 7 liaisons (29 contraintes), contre MBDyn")
    ds = []
    for h in (8e-3, 1e-3, 5e-4):
        v, phi0, phid0, phif = run(h)
        d = float(np.linalg.norm(v - ref))
        ds.append((h, d, phi0, phid0, phif))
        print(f"     h {h:.0e} : ({v[0]:+.6f},{v[1]:+.6f}) contre MBDyn ({ref[0]:+.6f},{ref[1]:+.6f}) "
              f"— écart {d:.2e} m")
    # · l'état initial du benchmark est ADMISSIBLE : Φ et Φ̇ nuls. Sans les
    #   vitesses ANGULAIRES (que j'avais d'abord oubliées) Φ̇ vaut 0,5 et la
    #   trajectoire n'a plus rien à voir.
    assert ds[0][2] < 1e-12 and ds[0][3] < 1e-12, ("état initial non admissible", ds[0])
    assert ds[-1][4] < 1e-12, ("les contraintes ne sont pas tenues", ds[-1])
    # · au MÊME pas fin, les deux solveurs sont au même point
    assert ds[-1][1] < 1e-5, ("vinkulum s'écarte de MBDyn sur 6barmech", ds)
    # · et le raffinement RAPPROCHE : c'est ce qui distingue une erreur
    #   d'intégration d'une erreur de modèle
    assert ds[0][1] > 10 * ds[-1][1], ("le raffinement ne rapproche pas", ds)
    print(f"     au même pas fin, les deux solveurs sont à {ds[-1][1]:.1e} m sur une course de 1 m")
    print("     (comparer à pas DIFFÉRENTS mélange erreur de modèle et d'intégration :")
    print("      vinkulum seul raffiné donnait 6,0e-3 — MBDyn se déplace exactement d'autant.)")
    return dict(ecarts=[(h, d) for h, d, *_ in ds])


def srscm():
    """srskm — le *spatial slider-crank* du jeu de tests de MBDyn.

    Sixième problème posé ailleurs, et le premier qui demande une liaison que
    le noyau n'avait pas : un **joint de Cardan**. Trois corps, quatre
    liaisons (17 contraintes pour 18 ddl, mobilité 1), vilebrequin lancé à
    6 rad/s puis lâché sous gravité — aucun moteur, aucune commande. Rien
    n'est de moi : géométrie, inerties, choix des contraintes et référence
    viennent du deck `srscm`.

    POURQUOI LE CARDAN EST UNE BRIQUE ET NON UN RÉGLAGE. On serait tenté de
    l'écrire `liaison(bloque_r = [k])` — deux rotations libres, une bloquée,
    c'est bien le compte. Mais cette liaison-là annule la composante k de la
    partie ANTISYMÉTRIQUE de uₐᵀu_b, soit ½(uₐ₃·u_b₂ − uₐ₂·u_b₃) : le second
    terme n'a rien à faire dans un cardan. Les deux ne coïncident qu'au
    premier ordre, et ce mécanisme balaye des angles pleins. Le joint réel
    est `Φ = (R_a·n_a)·(R_b·n_b) = 0` — les deux axes de la croix restent
    perpendiculaires — et c'est ce que `Noyau.cardan` pose.

    LE CONTRÔLE QUI NE DOIT RIEN AU SOLVEUR : les vitesses initiales. Le deck
    lâche le vilebrequin à ω = 6 rad/s et laisse bielle et coulisseau à
    l'arrêt — un état INADMISSIBLE, que MBDyn répare par son assemblage. On
    les résout ici à la main (4 équations linéaires : la rotule au maneton,
    la glissière, et la dérivée de la contrainte de cardan), et elles
    retombent sur celles de MBDyn au dernier chiffre. Le montage est donc
    juste avant qu'aucun pas de temps n'ait été fait.
    """
    L_C, L_R, OM = 0.08, 0.30, 6.0
    A = np.array([0.0, 0.1, 0.12])
    B = A + np.array([0.0, 0.0, L_C])
    C = np.array([np.sqrt(L_R ** 2 - 0.1 ** 2 - (0.12 + L_C) ** 2), 0.0, 0.0])
    # MBDyn au pas fin (h = 1e-4) : position des noeuds a t = 5 s
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

    # vitesses initiales admissibles : inconnues (vx du coulisseau, omega bielle)
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

    def run(h, cardan=True, ggl=False):
        n, (ro, bl) = monte(cardan)
        phi0 = max(abs(x) for x in n.phi()), max(abs(x) for x in n.phi_dot())
        n.simule(5.0, h, tous=10 ** 9, ggl=ggl)
        p = n.etat()[1]
        e = max(float(np.max(np.abs(np.array(p[i]) - REF[k])))
                for k, i in (("rod", ro), ("block", bl)))
        return e, phi0, max(abs(x) for x in n.phi()), max(abs(x) for x in n.phi_dot())

    print("  15. srskm — bielle-manivelle SPATIALE, joint de Cardan, contre MBDyn")
    dw = float(np.max(np.abs(w_ro - REF_W)))
    print(f"     vitesses initiales resolues a la main : omega bielle "
          f"({w_ro[0]:.2f},{w_ro[1]:.2f},{w_ro[2]:.2f}) — MBDyn a {dw:.1e} pres")
    ds = []
    for h in (4e-3, 5e-4):
        e, phi0, phif, _ = run(h)
        ds.append((h, e, phi0, phif))
        print(f"     h {h:.0e} : ecart max aux noeuds {e:.2e} m")
    sans, _, _, _ = run(4e-3, cardan=False)
    # GGL TRAVERSE LE POINT MORT (5 sept., soir). Jusqu'à ce soir `simule(ggl=True)`
    # cassait ici à t = 2,724 s, où cond(GGᵀ) fait un pic ×4 800 — et le README
    # disait « GGL ne traverse pas une configuration quasi singulière ». La cause
    # n'était pas la configuration : c'était la PRÉCISION des colonnes ζ du
    # jacobien, des différences finies à ~1e-8 relatif, que le conditionnement
    # multipliait en un incrément faux au pour-cent. Colonnes exactes
    # (`tangent::jacobien_ad`, produit ∂r/∂q · βh²G₀ᵀ) : le pic ne gêne plus.
    e_g, _, _, phid_g = run(4e-3, ggl=True)

    # · l'assemblage en vitesse est un fait de CINEMATIQUE : il ne doit rien
    #   au solveur, et il tombe sur MBDyn au dernier chiffre
    assert dw < 1e-9, ("les vitesses initiales ne sont pas celles de MBDyn", w_ro)
    # · etat initial admissible, contraintes tenues jusqu'au bout
    assert max(ds[0][2]) < 1e-12 and ds[-1][3] < 1e-12, ("contraintes", ds)
    # · au pas fin, les deux solveurs sont au meme point
    assert ds[-1][1] < 1e-5, ("vinkulum s'ecarte de MBDyn sur srskm", ds)
    # · ordre 2 observe entre les deux pas (facteur 8 => facteur ~64)
    r = ds[0][1] / ds[-1][1]
    assert r > 20, ("le raffinement ne rapproche pas assez", ds)
    # · CONTROLE NEGATIF : sans la croix, le mecanisme a deux ddl et part
    #   ailleurs. Sans lui, un vert ne prouverait pas que la brique sert.
    assert sans > 20 * ds[0][1], ("le cardan ne change rien — brique inerte", sans)
    print(f"     ordre observe {np.log(r) / np.log(8):.2f} · sans la croix, l'ecart "
          f"passe a {sans:.2e} m (controle negatif)")
    print(f"     GGL au meme pas : ecart {e_g:.2e} m, |Φ̇| {phid_g:.1e} — le point mort "
          f"de t = 2,72 s est traverse")
    # · le point mort est traversé, Φ̇ au plancher, et l'écart reste du même
    #   ordre que l'index 3 (un facteur sur la constante, pas sur l'ordre)
    assert phid_g < 1e-9 and e_g < 10 * ds[0][1], ("GGL ne traverse plus le point mort", e_g, phid_g)
    return dict(ecarts=[(h, e) for h, e, *_ in ds], sans_cardan=sans, ggl=(e_g, phid_g))

def lateral_buckling():
    """Flambement latéral d'une poutre mince — le banc publié de Multibody
    System Dynamics (2016) 37:29–48, implémenté par MBDyn dans son jeu de tests.

    Septième problème posé ailleurs, et le premier qui juge la
    **RAIDEUR GÉOMÉTRIQUE** : une poutre encastrée de 1 m, mince (EI_z/EI_y =
    1/100), tirée en bout par une bielle qu'une manivelle fait tourner d'un
    demi-tour en 0,4 s. La poutre plie dans son plan raide, la compression
    monte, et au-delà d'un seuil elle DÉVERSE latéralement. Rien de tout cela
    n'existe dans un modèle linéarisé : c'est la précontrainte qui rend le
    mode latéral instable.

    CE QU'IL FAUT JUGER, ET CE QU'IL NE FAUT PAS. Après le déversement la
    poutre oscille latéralement vite et fort : une position ponctuelle à
    t = 0,8 s ne dit rien (la phase dérive, MBDyn lui-même passe de −0,022 à
    +0,023 en 0,05 s). Les deux grandeurs robustes sont l'INSTANT du
    déversement et l'AMPLITUDE de l'oscillation qui suit. C'est aussi ce que
    compare le papier.

    LE CONTRÔLE NÉGATIF EST ICI LE PLUS FORT DE TOUS LES BANCS. Le benchmark
    porte un défaut initial `D = 0,1 mm` — sans lui le modèle est
    rigoureusement symétrique et ne PEUT pas déverser. Retiré, le noyau rend
    y ≡ 0 au bit près sur 1 600 pas : il ne fabrique aucune asymétrie. Un
    solveur qui en fabriquerait flamberait quand même, et son bel accord sur
    l'amplitude ne vaudrait rien.

    CE QUI DIFFÈRE DE MBDyn, DÉCLARÉ : lui discrétise en `beam3` (trois
    nœuds, deux points de Gauss), nous en poutres à deux nœuds, et la masse
    est concentrée aux nœuds des deux côtés mais pas selon la même règle. Le
    verdict se lit donc en CONVERGENCE DE MAILLAGE, comme sur Princeton.
    La commande de manivelle est un cosinus que MBDyn a en primitive et que
    le noyau approche par une table dense : de 2 000 à 32 000 points, les
    deux chiffres publiés ne bougent pas d'un millionième — l'approximation
    n'est pas dans le résultat.
    """
    T, D = 0.4, 0.1e-3
    SEG = dict(
        BEAM=dict(L=1.0, EA=73e6, GAY=5.025e6, GAZ=23.40e6, GJ=877.2, EJY=60830,
                  EJZ=608.3, M=2.68, J22=2233e-6, J33=22.33e-6),
        LINK=dict(L=0.25, EA=33.02e6, GAY=10.81e6, GAZ=10.81e6, GJ=914.5, EJY=1189,
                  EJZ=1189, M=1.212, J22=43.65e-6, J33=43.65e-6),
        CRANK=dict(L=0.05, EA=132.1e6, GAY=43.22e6, GAZ=43.22e6, GJ=14630, EJY=19020,
                   EJZ=19020, M=4.85, J22=698.3e-6, J33=698.3e-6))
    for d in SEG.values():
        d["J11"] = d["J22"] + d["J33"]
    # MBDyn sur le meme deck, h = T/800 : point milieu de la poutre
    REF_Y, REF_T = 0.044984, 0.1290

    def rep(v1, v2):
        e1 = np.asarray(v1, float) / np.linalg.norm(v1)
        e2 = np.asarray(v2, float) - np.dot(v2, e1) * e1
        return np.column_stack([e1, e2 / np.linalg.norm(e2), np.cross(e1, e2 / np.linalg.norm(e2))])

    R_B, R_L, R_K = np.eye(3), rep([0, 0, -1], [0, 1, 0]), rep([-1, 0, 0], [0, 1, 0])

    def marche(nb, nl, nk, d=D, npts=8000):
        R0, Tp = np.zeros(3), np.array([1.0, 0.0, 0.0])
        Cp = Tp + [0, d, 0]
        Bp = Cp + [0, 0, -0.25]
        Gp = Bp + [-0.05, 0, 0]
        n = Noyau([0.0, 0.0, 0.0])                 # ce benchmark est sans gravite
        idx = {}
        for lab, ne, o, rm in (("BEAM", nb, R0, R_B), ("LINK", nl, Cp, R_L),
                               ("CRANK", nk, Bp, R_K)):
            s, dl = SEG[lab], SEG[lab]["L"] / ne
            for i in range(ne + 1):
                la = dl if 0 < i < ne else dl / 2.0
                tige = s["M"] * la ** 3 / 12.0     # masse concentree : le troncon reste une tige
                idx[(lab, i)] = n.corps(
                    f"{lab}{i}", s["M"] * la,
                    [s["J11"] * la, 0, 0, 0, s["J22"] * la + tige, 0, 0, 0, s["J33"] * la + tige],
                    list(o + rm[:, 0] * (i * dl)), rot=list(rm.flatten()))
            for i in range(ne):
                n.poutre(f"{lab}p{i}", idx[(lab, i)], idx[(lab, i + 1)], s["EA"], s["GAY"],
                         s["GJ"], s["EJY"], ei3=s["EJZ"], ga3=s["GAZ"])
        n.liaison("clamp", None, idx[("BEAM", 0)], pa=list(R0))
        n.liaison("C", idx[("BEAM", nb)], idx[("LINK", 0)], pa=list(Cp - Tp),
                  bloque_t=[0, 1, 2], bloque_r=[])
        n.liaison("B", idx[("LINK", nl)], idx[("CRANK", 0)], pa=[0, 0, 0],
                  ra=list(R_L.T.flatten()), bloque_t=[0, 1, 2], bloque_r=[0, 2])
        # « cosine, 0, pi/T, pi/2, half » de MBDyn : (pi/2)(1 - cos(pi t/T)), figee a pi
        tb = []
        for k in range(npts + 1):
            u = 2 * T * k / npts
            tb += [u, np.pi / 2 * (1 - np.cos(np.pi * min(u, T) / T))]
        n.liaison("G", None, idx[("CRANK", nk)], pa=list(Gp), bloque_t=[0, 1, 2],
                  bloque_r=[0, 1, 2], cible_r=([0.0, -1.0, 0.0], ("table", tb)))
        im = idx[("BEAM", nb // 2)]
        h = T / 800
        ys = [(fr[0], abs(fr[1][im][1])) for fr in n.simule(2 * T, h, tous=1)]
        ym = max(y for _, y in ys)
        td = next((t for t, y in ys if y > 1e-3), float("nan"))
        return ym, td

    print("  16. flambement lateral (MSD 2016 37:29-48, deck MBDyn) — raideur geometrique")
    ds = []
    for m in ((10, 2, 2), (20, 4, 4)):
        ym, td = marche(*m)
        ds.append((m, ym, td))
        print(f"     {m[0]:2d} elements de poutre : |y| max {ym:.6f} m (MBDyn {REF_Y:.6f}) "
              f"· deversement a t = {td:.4f} s (MBDyn {REF_T:.4f})")
    sym, _ = marche(10, 2, 2, d=0.0)

    e = abs(ds[-1][1] - REF_Y) / REF_Y
    # · l'amplitude post-flambement, au maillage fin
    assert e < 5e-3, ("amplitude de deversement hors tolerance", ds)
    # · et le raffinement RAPPROCHE : erreur de discretisation, pas de modele
    assert abs(ds[0][1] - REF_Y) > 3 * abs(ds[-1][1] - REF_Y), ("le raffinement ne rapproche pas", ds)
    # · l'INSTANT du deversement, a un pas de temps pres
    assert abs(ds[-1][2] - REF_T) <= 1.5 * T / 800, ("instant de deversement", ds)
    # · CONTROLE NEGATIF : sans le defaut initial, le modele est symetrique et
    #   ne PEUT pas deverser. Un solveur qui fabrique de l'asymetrie flamberait
    #   quand meme, et son accord sur l'amplitude ne vaudrait rien.
    assert sym == 0.0, ("le noyau fabrique de l'asymetrie : y != 0 sans defaut", sym)
    print(f"     ecart au maillage fin {100 * e:.3f} % · sans le defaut initial, "
          f"|y| = {sym:.1e} sur 1600 pas (controle negatif)")
    return dict(maillages=[(m[0], y, t) for m, y, t in ds], sans_defaut=sym)

def arbre_tournant():
    """rotatingshaft — « Stability of a rotating shaft », le second banc publié
    de Multibody System Dynamics (2016) 37:29–48, dans l'implémentation du jeu
    de tests MBDyn.

    Huitième problème posé ailleurs. Un arbre tubulaire de 6 m bi-encastré en
    flexion, 32 éléments, portant en son milieu un disque de 70,6 kg EXCENTRÉ
    de 50 mm — le balourd. La rotation est imposée en deux rampes (0 → 48 puis
    48 → 72 rad/s), soit 22 tours, et il n'y a AUCUN amortissement : rien ne
    ramène la réponse vers un régime établi, tout ce qu'on lit est du
    transitoire.

    IL JUGE DEUX CHOSES QUE LES AUTRES BANCS NE JUGENT PAS.

    · **L'ANALYSE STATIQUE contre un tiers.** MBDyn atteint l'équilibre sous
      gravité en rampant la pesanteur sur dix secondes avec un amortisseur ;
      on le résout directement par `statique`. Les deux tombent sur la même
      flèche, et la théorie des poutres bi-encastrées aussi — trois chemins
      indépendants sur un même millimètre.
    · **Une pièce EXCENTRÉE.** Le disque ne peut ni être fondu dans le nœud
      (le balourd disparaîtrait) ni y être recentré (l'axe de la poutre ferait
      un coude de 48 mm) : c'est un corps à part, soudé au nœud. C'est ce que
      ferait n'importe quel utilisateur, et rien ne l'avait jamais exercé.

    CE QU'IL A TROUVÉ, ET C'EST LE PLUS UTILE. En le montant, j'ai bloqué la
    rotation du nœud d'encastrement ET imposé son angle par une seconde
    liaison. Deux contraintes qui se CONTREDISENT. Le noyau en a neutralisé
    une comme redondante — le mécanisme même qui fait passer un quatre-barres
    surcontraint ou un Bennett — et a rendu, sans un mot, un arbre qui ne
    tournait pas du tout. Φ valait 0,96.

    > La différence entre une redondance LICITE et une contradiction se lit
    > sur Φ, pas sur le rang. `Modele.tol_phi` la lit désormais après chaque
    > pas accepté, et le contrôle négatif ci-dessous rejoue la faute.
    """
    import math
    NB, L, rho, rI, rO = 16, 6.0, 7800.0, 0.045, 0.05
    ml = math.pi * rho * (rO ** 2 - rI ** 2)
    jxx = math.pi / 4 * rho * (rO ** 4 - rI ** 4)
    mD, rD, tD, dz = 70.573, 0.24, 0.05, 0.05
    jd11, jdxx = mD * rD ** 2 / 2, mD * (rD ** 2 / 4 + tD ** 2 / 12)
    EA, GA, GJ, EJ = 313.4e6, 60.5e6, 272.7e3, 354.5e3
    REF_Z, REF_A = -3.3100e-3, 0.26862         # MBDyn : fleche a t=0, excursion max

    def theta(t):                              # l'angle : les deux rampes de MBDyn, integrees
        a = (24 * t - 24 * math.sin(2 * math.pi * t) / (2 * math.pi)
             if t <= 0.5 else 12 + 48 * (t - 0.5))
        if t <= 1.0:
            b = 0.0
        elif t <= 1.25:
            b = 12 * (t - 1) - 12 * math.sin(4 * math.pi * (t - 1)) / (4 * math.pi)
        else:
            b = 3 + 24 * (t - 1.25)
        return a + b

    def monte(ne, spin, fautif=False):
        n = Noyau([0.0, 0.0, -9.81])
        dl, idx = L / ne, []
        for i in range(ne + 1):
            la = dl if 0 < i < ne else dl / 2.0
            dm = la * ml
            djxx = dm * la ** 2 / 12 + la * jxx
            idx.append(n.corps(f"s{i}", dm, [2 * jxx * la, 0, 0, 0, djxx, 0, 0, 0, djxx],
                               [i * dl, 0.0, 0.0]))
        for i in range(ne):
            n.poutre(f"p{i}", idx[i], idx[i + 1], EA, GA, GJ, EJ)
        c = ne // 2
        dq = n.corps("disque", mD, [jd11, 0, 0, 0, jdxx, 0, 0, 0, jdxx], [c * dl, 0.0, dz])
        n.liaison("soudure", idx[c], dq, pa=[0, 0, dz])
        tb = []
        for k in range(25001):
            u = 2.5 * k / 25000
            tb += [u, theta(u)]
        if fautif:                             # LA FAUTE, rejouee : deux liaisons pour un axe
            n.liaison("R", None, idx[0], pa=[0, 0, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
            n.liaison("spin", None, idx[0], pa=[0, 0, 0], bloque_t=[], bloque_r=[0],
                      cible_r=([1.0, 0.0, 0.0], ("table", tb)))
        else:
            n.liaison("R", None, idx[0], pa=[0, 0, 0], bloque_t=[0, 1, 2], bloque_r=[0, 1, 2],
                      **({"cible_r": ([1.0, 0.0, 0.0], ("table", tb))} if spin else {}))
        n.liaison("T", None, idx[-1], pa=[L, 0, 0], bloque_t=[1, 2], bloque_r=[1, 2])
        return n, idx

    print("  17. rotatingshaft (MSD 2016 37:29-48, deck MBDyn) — statique, balourd, 22 tours")
    n, idx = monte(32, spin=False)
    r, it = n.statique()
    zc = n.etat()[1][idx[16]][2]
    ana = -(ml * 9.81 * L ** 4 / (384 * EJ) + mD * 9.81 * L ** 3 / (192 * EJ))
    print(f"     STATIQUE : fleche milieu {zc * 1000:+.4f} mm · MBDyn {REF_Z * 1000:+.4f} · "
          f"poutre bi-encastree {ana * 1000:+.4f} ({it} it, ‖r‖ {r:.1e} N)")

    ds = []
    for ne in (16, 32):
        n, idx = monte(ne, spin=True)
        n.statique()
        h = n.simule(2.5, 1e-3, tous=1)
        c = idx[ne // 2]
        am = max(math.hypot(f[1][c][1], f[1][c][2]) for f in h)
        ds.append((ne, am, max(abs(x) for x in n.phi())))
        print(f"     {ne:2d} elements : excursion max du milieu {am:.5f} m "
              f"(MBDyn {REF_A:.5f}) · |Φ| {ds[-1][2]:.1e}")

    try:
        monte(16, spin=True, fautif=True)[0].simule(0.01, 1e-3, tous=10 ** 9)
        dit = ""
    except (RuntimeError, ValueError) as ex:
        dit = str(ex)

    # · la statique tombe sur MBDyn ET sur la theorie des poutres
    assert abs(zc - REF_Z) / abs(REF_Z) < 0.01, ("fleche statique", zc, REF_Z)
    assert abs(zc - ana) / abs(ana) < 0.01, ("la theorie ne confirme pas", zc, ana)
    # · la dynamique : 22 tours sans amortissement, aucun regime etabli
    assert abs(ds[-1][1] - REF_A) / REF_A < 0.03, ("excursion max hors tolerance", ds)
    assert ds[-1][2] < 1e-9, ("les contraintes ne sont pas tenues", ds)
    # · CONTROLE NEGATIF : la faute que ce banc a servi a trouver. Deux liaisons
    #   contradictoires sur le meme axe doivent etre DITES, pas neutralisees en
    #   silence -- sinon l'arbre ne tourne pas et rien ne le signale.
    assert "Φ" in dit and "redondante" in dit, ("la contradiction passe en silence", dit)
    print(f"     ecart dynamique {100 * abs(ds[-1][1] - REF_A) / REF_A:.2f} % · "
          f"la contradiction encastrement/rotation est refusee (controle negatif)")
    return dict(fleche=zc, excursions=[(ne, a) for ne, a, _ in ds])

def invariance_unites():
    """LE MÊME MÉCANISME, COTÉ DANS CINQ SYSTÈMES D'UNITÉS.

    Ce banc n'a pas de référence externe et n'en a pas besoin : rien de
    physique ne change d'un jeu à l'autre, seuls les NOMBRES. Un noyau
    agnostique doit rendre la même trajectoire une fois ramenée aux mêmes
    unités. Toute grandeur écrite en dur — une tolérance en newtons, une borne
    en mètres — le fait échouer.

    POURQUOI CE BANC EXISTE. Les huit modèles étrangers ont trouvé, en tout,
    ZÉRO erreur de physique et quatre défauts, dont deux de la même famille :
    une grandeur ABSOLUE dans un code censé ne pas connaître l'unité de son
    utilisateur. Les chercher un par un, au hasard des bancs, est la mauvaise
    méthode — celle-ci les trouve toutes d'un coup, et elle a immédiatement
    payé : le clamp d'accélération initiale valait 0,05 **mètre**, ce qui
    décalait la réponse de 3,7 % à l'échelle d'un MEMS, en silence. Il vaut
    désormais 5 % de la taille du modèle.

    ET IL A TROUVÉ UN DÉFAUT QUE J'AVAIS INTRODUIT UNE HEURE PLUS TÔT : le
    garde-fou de stagnation de `statique`, écrit « moins de 1 % de progrès
    cinq fois de suite », coupait un Newton AMORTI qui convergeait très bien —
    il descend par paliers avant d'entrer en régime quadratique. Sur une
    fenêtre de dix itérations, la même statique tient jusqu'à un facteur
    d'échelle 1e4 au lieu de 3e3.

    LA LIMITE EST PUBLIÉE, PAS EXPLIQUÉE. Au-delà d'un facteur ~1e5 sur les
    LONGUEURS, `statique` ne converge plus ; son résidu relatif se dégrade en
    L². Deux causes ont été supposées et RÉFUTÉES par la mesure —
    l'équilibrage des lignes de contrainte (essayé, sans effet, retiré) et le
    conditionnement de K (mesuré à 5e17 *aussi* à l'échelle SI, où tout
    converge). On publie donc la borne mesurée et on ne fabrique pas de
    troisième hypothèse. La dynamique, elle, est invariante sur les cinq
    systèmes.
    """
    # (nom, facteur de longueur, de masse, de temps)
    JEUX = [("SI            ", 1.0, 1.0, 1.0),
            ("mm-t-s        ", 1e3, 1e-3, 1.0),
            ("pouce-livre-ms", 39.3701, 2.20462, 1e3),
            ("micro (MEMS)  ", 1e6, 1e9, 1e6),
            ("kilo (civil)  ", 1e-3, 1e-6, 1e-3)]

    def cas(al, am, at, stat, iters=60):
        af = am * al / at ** 2                      # force
        ae = af * al ** 2                           # EI, GJ
        n = Noyau([0.0, 0.0, -9.81 * al / at ** 2])
        lg, ms, ei, ne = 1.0 * al, 2.0 * am, 50.0 * ae, 6
        dl, idx = lg / ne, []
        for i in range(ne + 1):
            la = dl if 0 < i < ne else dl / 2
            dm = ms / lg * la
            j = dm * la ** 2 / 12
            idx.append(n.corps(f"c{i}", dm, [j, 0, 0, 0, j, 0, 0, 0, j], [i * dl, 0, 0]))
        for i in range(ne):
            n.poutre(f"p{i}", idx[i], idx[i + 1], 1e4 * af, 5e3 * af, 0.5 * ei, ei)
        n.liaison("enc", None, idx[0], pa=[0, 0, 0])
        if stat:
            # paliers_max=1 : on mesure la statique NUE. La continuation
            # automatique masquerait justement ce que ce banc cherche — la
            # borne d'echelle au-dela de laquelle Newton ne demarre plus.
            r, it = n.statique(iters=iters, paliers_max=1)
            return n.etat()[1][idx[-1]][2] / al, r
        n.simule(0.4 * at, 1e-3 * at, tous=10 ** 9)
        return n.etat()[1][idx[-1]][2] / al, max(abs(x) for x in n.phi())

    print("  18. invariance d'unites — le meme mecanisme dans cinq systemes")
    out = {}
    for quoi in ("dynamique", "statique"):
        vals, ref = [], None
        for nom, al, am, at in JEUX:
            try:
                z, r = cas(al, am, at, quoi == "statique")
                ref = z if ref is None else ref
                vals.append((nom, z, abs(z - ref) / abs(ref), None))
            except (ValueError, RuntimeError) as ex:
                vals.append((nom, float("nan"), float("nan"), type(ex).__name__))
        ok = [v for v in vals if v[3] is None]
        pire = max(v[2] for v in ok)
        ref = ok[0][1]
        print(f"     {quoi:10s} : {len(ok)}/{len(JEUX)} systemes · fleche {ref:+.9f} m "
              f"· ecart max {pire:.1e}"
              + ("" if len(ok) == len(JEUX)
                 else "  (refus : " + ", ".join(v[0].strip() for v in vals if v[3]) + ")"))
        out[quoi] = (len(ok), pire)

    # borne d'echelle de la statique, mesuree
    bornes = []
    for al in (1e2, 1e4, 1e5):
        try:
            z, r = cas(al, 1.0, 1.0, True, iters=400)
            bornes.append((al, z, r, True))
        except (ValueError, RuntimeError):
            bornes.append((al, float("nan"), float("nan"), False))
    tenus = [b for b in bornes if b[3]]
    print(f"     statique : exacte au 9e chiffre jusqu'a un facteur de longueur "
          f"{max(b[0] for b in tenus):.0e}, refus a {min((b[0] for b in bornes if not b[3]), default=float('inf')):.0e} "
          f"(residu relatif en L^2 — limite MESUREE, cause non identifiee)")

    # · la dynamique est invariante sur les CINQ systemes
    assert out["dynamique"][0] == len(JEUX), ("un systeme d'unites refuse en dynamique", vals)
    assert out["dynamique"][1] < 1e-8, ("la dynamique depend des unites", out)
    # · la statique l'est sur ceux qu'elle tient
    assert out["statique"][1] < 1e-8, ("la statique depend des unites", out)
    # · et elle est EXACTE, pas seulement invariante, sur toute la plage tenue
    assert max(abs(b[1] - tenus[0][1]) for b in tenus) < 1e-9, ("resultat different", bornes)
    return dict(dynamique=out["dynamique"], statique=out["statique"],
                bornes=[(b[0], b[3]) for b in bornes])

def quatre_barres_flexible():
    """fourbar — le quatre-barres de Bauchau, troisième banc publié de
    Multibody System Dynamics (2016) 37:29–48, dans l'implémentation du jeu de
    tests MBDyn.

    Neuvième problème posé ailleurs, et **le seul qui ne peut pas bouger sans
    la flexibilité**. Trois barres, quatre pivots — et les axes sont
    VOLONTAIREMENT désalignés : celui du pivot C est incliné de 5°. Le
    mécanisme rigide correspondant est cinématiquement BLOQUÉ ; ce qui le fait
    tourner, c'est la déformation des barres. Un solveur qui traiterait les
    corps comme rigides ne rendrait pas un résultat approché : il ne rendrait
    rien du tout.

    La manivelle est entraînée à 0,6 rad/s pendant deux tours avant t = 0,
    puis on regarde jusqu'à t = 12 s — trois tours de plus.

    ÉTAT INITIAL : c'est une DONNÉE du deck, pas un résultat, et il est
    embarqué tel que MBDyn le construit (`donnees/fourbar.json`). Le
    reconstruire demandait de réimplémenter la sémantique des repères MBDyn en
    cascade — une source d'erreur qui n'apprend rien sur le noyau. Ce qui se
    vérifie, en revanche, c'est qu'il est ADMISSIBLE pour NOS contraintes :
    Φ et Φ̇ y sont nuls à la précision machine, ce qui n'irait pas de soi si
    l'un des quatre repères de liaison était mal posé.

    CE QU'IL A TROUVÉ. `phi_dot` rendait **Φ̇ = 0,6 sur cet état parfaitement
    admissible** : il calculait G·u et oubliait ∂Φ/∂t, qui n'est nul que si
    aucune liaison n'est PILOTÉE. Ici la manivelle l'est. Un état initial
    licite était donc déclaré faux — le pire que puisse faire un contrôle, et
    il m'a fait chercher une heure du mauvais côté.
    """
    import json
    import os
    d = json.load(open(os.path.join(os.path.dirname(__file__), "donnees", "fourbar.json")))
    NB, OM, TH = 5, 0.6, np.radians(5.0)
    P = {1: dict(L=0.12, EA=5.29920e7, GAY=1.68803e7, GAZ=1.68803e7, GJ=7.33488e2,
                 EJY=1.13050e3, EJZ=1.13050e3, m=1.99680, J1=8.51968e-5,
                 J2=4.25984e-5, J3=4.25984e-5),
         2: dict(L=0.24, EA=5.29920e7, GAY=1.68803e7, GAZ=1.68803e7, GJ=7.33488e2,
                 EJY=1.13050e3, EJZ=1.13050e3, m=1.99680, J1=8.51968e-5,
                 J2=4.25984e-5, J3=4.25984e-5),
         3: dict(L=0.12, EA=1.32480e7, GAY=4.22008e6, GAZ=4.22008e6, GJ=4.58430e1,
                 EJY=7.06560e1, EJZ=7.06560e1, m=4.99200e-1, J1=5.32480e-6,
                 J2=2.66240e-6, J3=2.66240e-6)}

    def rc(v1, v2):
        e1 = np.asarray(v1, float) / np.linalg.norm(v1)
        e2 = np.asarray(v2, float) - np.dot(v2, e1) * e1
        e2 /= np.linalg.norm(e2)
        return np.column_stack([e1, e2, np.cross(e1, e2)])

    def rodrigues(p):
        th = float(np.linalg.norm(p))
        if th < 1e-14:
            return np.eye(3)
        k = np.asarray(p) / th
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)

    R1, R2, R3 = rc([0, 0, 1], [-1, 0, 0]), rc([1, 0, 0], [0, 0, 1]), rc([0, 0, -1], [1, 0, 0])
    R23 = R2 @ rc([np.cos(TH), 0, -np.sin(TH)], [0, 1, 0])
    A = np.zeros(3)
    B = np.array([0.0, 0.0, P[1]["L"]])
    C = np.array([P[2]["L"], 0.0, P[1]["L"]])
    D = np.array([P[2]["L"], 0.0, 0.0])

    def monte(theta_c=TH):
        ne = 2 * NB
        r23 = R2 @ rc([np.cos(theta_c), 0, -np.sin(theta_c)], [0, 1, 0])
        n = Noyau([0.0, 0.0, 0.0])                    # ce benchmark est sans gravite
        bar = {}
        for j in (1, 2, 3):
            s, dl, bar[j] = P[j], P[j]["L"] / ne, []
            for i in range(ne + 1):
                la = dl if 0 < i < ne else dl / 2
                dm = s["m"] * la
                tige = dm * la ** 2 / 12
                ci = d["ci"][str(j * 1000 + i)]
                bar[j].append(n.corps(
                    f"b{j}n{i}", dm,
                    [s["J1"] * la, 0, 0, 0, s["J2"] * la + tige, 0, 0, 0, s["J3"] * la + tige],
                    ci[:3], rot=list(rodrigues(ci[3:6]).flatten()), v=ci[6:9], w=ci[9:12]))
            for i in range(ne):
                n.poutre(f"p{j}{i}", bar[j][i], bar[j][i + 1], s["EA"], s["GAY"], s["GJ"],
                         s["EJY"], ei3=s["EJZ"], ga3=s["GAZ"])
        rb1 = np.array(d["ci"][str(1000 + ne)][:3])
        rb2 = np.array(d["ci"][str(2000 + ne)][:3])
        n.liaison("A", None, bar[1][0], pa=list(A), ra=list(R1.flatten()),
                  bloque_t=[0, 1, 2], bloque_r=[0, 1, 2],
                  cible_r=([0.0, 0.0, 1.0], ("lineaire", [0.0, -OM])))
        n.liaison("B", bar[1][ne], bar[2][0], pa=list(R1.T @ (B - rb1)),
                  bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.liaison("C", bar[2][ne], bar[3][0], pa=list(R2.T @ (C - rb2)),
                  ra=list((R2.T @ r23).flatten()), bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.liaison("D", None, bar[3][ne], pa=list(D), ra=list(R3.flatten()),
                  bloque_t=[0, 1, 2], bloque_r=[0, 1])
        return n, bar

    print("  19. fourbar (MSD 2016 37:29-48, deck MBDyn) — il ne tourne QUE parce qu'il plie")
    n, bar = monte()
    phi0 = max(abs(x) for x in n.phi())
    phid0 = max(abs(x) for x in n.phi_dot())
    t0 = float(d["ci"]["1000"][0]) * 0 - 2 * (2 * np.pi / OM)
    n.simule(d["t_fin"], 4e-3, tous=10 ** 9)
    p = n.etat()[1]
    ecarts = []
    for j in (1, 2, 3):
        for i in (0, NB, 2 * NB):
            ref = np.array(d["fin"][str(j * 1000 + i)])
            ecarts.append(float(np.max(np.abs(np.array(p[bar[j][i]]) - ref))))
    e = max(ecarts)
    print(f"     etat initial : |Φ| {phi0:.1e} · |Φ̇| {phid0:.1e} (admissible pour NOS liaisons)")
    print(f"     apres {d['t_fin'] - t0:.1f} s et 3 tours : ecart max sur 9 noeuds {e:.2e} m")

    # controle NEGATIF : sans les 5 degres, le mecanisme rigide est bloque.
    # On ne demande pas qu'il plante -- on demande qu'il fasse AUTRE CHOSE.
    n0, bar0 = monte(theta_c=0.0)
    try:
        n0.simule(t0 + 2.0, 4e-3, tous=10 ** 9)
        pd = float(np.max(np.abs(np.array(n0.etat()[1][bar0[2][NB]])
                                 - np.array(d["fin"]["2005"]))))
    except (RuntimeError, ValueError):
        pd = float("inf")

    assert phi0 < 1e-12, ("l'etat initial du deck viole nos contraintes", phi0)
    # · l'etat initial est ADMISSIBLE : c'est ce que phi_dot ne savait pas dire
    assert phid0 < 1e-9, ("phi_dot ignore la derivee de la cible pilotee", phid0)
    # · et la dynamique suit MBDyn sur trois tours
    assert e < 5e-5, ("vinkulum s'ecarte de MBDyn sur fourbar", ecarts)
    # · sans le desalignement, ce n'est pas le meme mecanisme
    assert pd > 100 * e, ("les 5 degres ne changent rien -- le banc ne juge pas ce qu'il croit", pd)
    print(f"     sans les 5 deg d'inclinaison en C (mecanisme rigide bloque) : {pd:.2e} m "
          f"d'ecart au meme instant (controle negatif)")
    return dict(ecart=e, phi0=phi0, phid0=phid0)

def treillis_de_boucles():
    """multibarmech — le dernier cas du jeu de tests MBDyn, et le plus GROS :
    un treillis de 5 × 5 boucles fermées, 55 corps, 190 contraintes, lâché
    sous gravité avec une vitesse de rotation initiale.

    Dixième problème posé ailleurs. Il n'apporte aucune physique neuve : ce
    qu'il exerce, c'est l'ÉCHELLE et la REDONDANCE. Ses auteurs choisissent
    les contraintes composante par composante — un pivot y bloque deux
    translations et zéro rotation — pour que le treillis ne soit pas
    surdéterminé. Un noyau qui traiterait chaque liaison comme un bloc de six
    n'en résoudrait aucune.

    Il est passé du premier coup, état initial RECONSTRUIT compris : Φ = 0
    exactement et Φ̇ = 1,8e-15 sans rien lire de la sortie MBDyn.

    ⚠ ET IL DONNE L'HORIZON DE VALIDITÉ D'UNE COMPARAISON DE SOLVEURS. Le
    mécanisme a ~140 degrés de liberté et aucun amortissement : sa trajectoire
    diverge. Mesuré, et c'est le chiffre qui compte :

    ```text
        a t = 10 s, MBDyn contre LUI-MEME au pas raffine x5 : 4,34 m d'ecart
        au meme instant, vinkulum contre MBDyn au meme pas  : 0,0158 m
    ```

    **Le solveur de référence se déplace 280 fois plus en changeant son propre
    pas que ne l'en sépare le nôtre.** Publier un écart à t = 10 s serait
    publier du bruit — dans un sens comme dans l'autre. On compare donc sur la
    fenêtre où la mesure a un sens, et on publie la borne au lieu de la taire.
    """
    import json
    import os
    d = json.load(open(os.path.join(os.path.dirname(__file__), "donnees", "multibarmech.json")))
    NY, NZ, LG, MS = 5, 5, 1.0, 1.0
    DP = 2 * np.pi / 3

    def monte():
        n = Noyau([0.0, 0.0, -9.81])
        jv = [MS * LG ** 2 / 12, 0, 0, 0, MS * LG ** 2 / 12, 0, 0, 0, 1.0]
        jh = [MS * LG ** 2 / 12, 0, 0, 0, 1.0, 0, 0, 0, MS * LG ** 2 / 12]
        v, h = {}, {}
        for iz in range(NZ):
            for iy in range(NY + 1):
                v[(iz, iy)] = n.corps(f"v{iz}_{iy}", MS, jv, [0.0, iy * LG, -iz * LG - LG / 2],
                                      v=[0.0, DP * (iz * LG + LG / 2), 0.0], w=[DP, 0.0, 0.0])
            for iy in range(1, NY + 1):
                h[(iz, iy)] = n.corps(f"h{iz}_{iy}", MS, jh,
                                      [0.0, (iy - 1) * LG + LG / 2, -(iz + 1) * LG],
                                      v=[0.0, DP * (iz + 1) * LG, 0.0])
        k = [0]

        def li(a, b, pa, bt, br):
            k[0] += 1
            n.liaison(f"j{k[0]}", a, b, pa=pa, bloque_t=bt, bloque_r=br)

        for iz in range(NZ):
            if iz == 0:
                li(None, v[(0, 0)], [0.0, 0.0, 0.0], [0, 1, 2], [1, 2])
            else:
                li(v[(iz, 0)], v[(iz - 1, 0)], [0.0, 0.0, LG / 2], [0, 1, 2], [1, 2])
            for iy in range(1, NY + 1):
                if iz == 0:
                    li(None, v[(0, iy)], [0.0, iy * LG, 0.0], [1, 2], [])
                else:
                    li(v[(iz, iy)], v[(iz - 1, iy)], [0.0, 0.0, LG / 2], [1, 2], [])
            li(v[(iz, 0)], h[(iz, 1)], [0.0, 0.0, -LG / 2], [0, 1, 2], [1, 2])
            for iy in range(1, NY):
                li(v[(iz, iy)], h[(iz, iy + 1)], [0.0, 0.0, -LG / 2], [1, 2], [])
            for iy in range(1, NY + 1):
                li(v[(iz, iy)], h[(iz, iy)], [0.0, 0.0, -LG / 2], [1, 2], [])
        return n, v, h

    print("  20. multibarmech (jeu de tests MBDyn) — 5x5 boucles, 55 corps, 190 contraintes")
    n, v, h = monte()
    phi0 = max(abs(x) for x in n.phi())
    phid0 = max(abs(x) for x in n.phi_dot())
    H = n.simule(2.0, 1e-3, tous=1)
    ecarts = {}
    for t in ("0.5", "1.0", "2.0"):
        k = min(int(round(float(t) / 1e-3)) - 1, len(H) - 1)
        p = H[k][1]
        e = 0.0
        for (iz, iy), i in v.items():
            e = max(e, max(abs(a - b) for a, b in zip(p[i], d["ref"][t][str(100 * iz + iy)])))
        for (iz, iy), i in h.items():
            e = max(e, max(abs(a - b) for a, b in zip(p[i], d["ref"][t][str(10000 + 100 * iz + iy)])))
        ecarts[t] = e
    c = d["chaos"]
    print(f"     etat initial RECONSTRUIT : |Φ| {phi0:.1e} · |Φ̇| {phid0:.1e}")
    print("     ecart max sur les 55 noeuds : "
          + " · ".join(f"t={t} s {e:.2e} m" for t, e in ecarts.items()))
    print(f"     horizon : a t={c['t']:.0f} s MBDyn contre LUI-MEME (pas /5) vaut {c['mbdyn_h1e3_contre_h2e4_max']:.2f} m,"
          f" nous {c['vinkulum_contre_mbdyn_max']:.1e} — {c['mbdyn_h1e3_contre_h2e4_max']/c['vinkulum_contre_mbdyn_max']:.0f}x plus")

    # · l'etat initial se RECONSTRUIT, il n'a pas ete lu chez MBDyn
    assert phi0 < 1e-12 and phid0 < 1e-12, ("etat initial non admissible", phi0, phid0)
    # · et les 190 contraintes tiennent jusqu'au bout
    assert max(abs(x) for x in n.phi()) < 1e-12, "les contraintes derivent"
    # · accord sur la fenetre ou la mesure a un sens
    assert max(ecarts.values()) < 1e-4, ("vinkulum s'ecarte de MBDyn", ecarts)
    # · et la borne au-dela de laquelle comparer n'a plus de sens est PUBLIEE,
    #   pas tue : le solveur de reference bouge 280x plus en raffinant son
    #   propre pas que ne l'en separe le notre.
    assert c["mbdyn_h1e3_contre_h2e4_max"] > 100 * c["vinkulum_contre_mbdyn_max"], (
        "l'argument d'horizon ne tient plus — rejouer la mesure", c)
    return dict(ecarts=ecarts, phi0=phi0, phid0=phid0)

def invariances():
    """CE QU'UN UTILISATEUR SUPPOSE SANS LE VÉRIFIER.

    Dix modèles étrangers ont été pris, et le critère qui reste — « n'importe
    qui pose SON modèle dedans et obtient un résultat juste » — ne se fabrique
    pas de l'intérieur. Ce banc en est le plus proche substitut honnête : il
    n'exige pas un résultat, il exige les propriétés qu'un utilisateur tient
    pour acquises **sans jamais les tester**, et dont la violation ne se voit
    pas sur un banc de physique bien orienté.

    Quatre, chacune tuant une classe entière de fautes :

    · **repère** — le modèle entier tourné et translaté doit rendre la même
      chose. Attrape toute hypothèse cachée d'alignement aux axes : « z est
      vertical », « la gravité est en −z », un produit vectoriel écrit à la
      main dans un seul plan ;
    · **ordre de déclaration** — permuter les corps ne doit rien changer.
      Attrape les dépendances à l'ordre d'assemblage et une détection de
      redondance qui neutraliserait « la dernière » ligne plutôt que la bonne ;
    · **origine des temps** — démarrer l'horloge à t = −3 ou t = 1000 avec la
      même loi de commande recalée doit donner la même trajectoire. Attrape un
      `t` supposé nul quelque part ;
    · **invariance galiléenne** — un système LIBRE auquel on ajoute une
      vitesse d'ensemble doit se déformer exactement pareil. C'est le test
      d'objectivité de la formulation de poutre, et beaucoup de formulations
      le ratent ; il est plus fort que la rotation statique du repère.

    LE CONTRÔLE NÉGATIF EST OBLIGATOIRE ICI. Un banc d'invariance passe
    toujours du premier coup si le code est correct — donc il ne prouve rien
    tant qu'on n'a pas vérifié qu'il MORD. On rejoue la faute exacte qu'un
    code aligné aux axes commettrait : tourner le modèle et laisser la
    gravité en −z. L'écart passe de 1e-14 à 2e-2.

    CE QUE ÇA NE PROUVE PAS : qu'un modèle inconnu passera. Une invariance
    respectée est nécessaire, pas suffisante. Elle dit seulement qu'aucune de
    ces quatre familles de fautes n'est présente.
    """
    def rot(ax, th):
        k = np.asarray(ax, float)
        k = k / np.linalg.norm(k)
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)

    J = [2e-3, 0, 0, 0, 2e-3, 0, 0, 0, 2e-3]
    PT = [np.array([0.3 * i, 0.0, 0.0]) for i in range(4)]

    def bati(R=None, T=None, ordre=None, t0=0.0, g_tournee=True):
        R = np.eye(3) if R is None else R
        T = np.zeros(3) if T is None else T
        g = R @ np.array([0.0, 0.0, -9.81]) if g_tournee else np.array([0.0, 0.0, -9.81])
        n = Noyau(list(g))
        idx = {}
        for i in (range(4) if ordre is None else ordre):
            idx[i] = n.corps(f"c{i}", 0.5, J, list(R @ PT[i] + T), rot=list(R.flatten()))
        for i in range(3):
            n.poutre(f"p{i}", idx[i], idx[i + 1], 1e5, 5e4, 20.0, 40.0)
        n.liaison("pivot", None, idx[0], pa=list(R @ PT[0] + T), ra=list(R.flatten()),
                  bloque_t=[0, 1, 2], bloque_r=[0, 1, 2],
                  cible_r=([0.0, 0.0, 1.0], ("lineaire", [-t0 * 0.5, 0.5])))
        return n, idx, R, T

    def bout(n, idx, R, T):
        return R.T @ (np.array(n.etat()[1][idx[3]]) - T)

    n, idx, R, T = bati()
    n.simule(1.0, 1e-3, tous=10 ** 9)
    ref = bout(n, idx, R, T)

    print("  21. invariances — ce qu'un utilisateur suppose sans le verifier")
    res = {}
    e = 0.0
    for nom, R_, T_ in (("rot 90 deg", rot([0, 0, 1], np.pi / 2), np.zeros(3)),
                        ("rot quelconque", rot([1, 2, 3], 0.6458), np.zeros(3)),
                        ("translation", np.eye(3), np.array([10.0, -5.0, 3.0])),
                        ("les deux", rot([1, 2, 3], 0.6458), np.array([10.0, -5.0, 3.0]))):
        m, i, a, b = bati(R_, T_)
        m.simule(1.0, 1e-3, tous=10 ** 9)
        e = max(e, float(np.max(np.abs(bout(m, i, a, b) - ref))))
    res["repere"] = e

    e = 0.0
    for o in ([3, 2, 1, 0], [1, 3, 0, 2], [2, 0, 3, 1]):
        m, i, a, b = bati(ordre=o)
        m.simule(1.0, 1e-3, tous=10 ** 9)
        e = max(e, float(np.max(np.abs(bout(m, i, a, b) - ref))))
    res["ordre"] = e

    e = 0.0
    for t0 in (5.0, -3.0, 1000.0):
        m, i, a, b = bati(t0=t0)
        st = m.etat()
        m.pose_etat(t0, st[1], st[2], st[3], st[4])
        m.simule(t0 + 1.0, 1e-3, tous=10 ** 9)
        e = max(e, float(np.max(np.abs(bout(m, i, a, b) - ref))))
    res["temps"] = e

    def libre(vd):
        n = Noyau([0.0, 0.0, 0.0])
        idx = [n.corps(f"c{i}", 0.5, J, [0.3 * i, 0.0, 0.0], v=list(vd)) for i in range(4)]
        for i in range(3):
            n.poutre(f"p{i}", idx[i], idx[i + 1], 1e5, 5e4, 20.0, 40.0)
        n.effort(idx[3], [0.0, 30.0, 0.0], [0.0, 0.0, 0.0])
        n.simule(0.5, 1e-3, tous=10 ** 9)
        p = np.array(n.etat()[1])
        return p - p[0]                      # la DEFORMATION, vue du premier corps

    d0 = libre(np.zeros(3))
    res["galilee"] = max(float(np.max(np.abs(libre(np.array(v)) - d0)))
                         for v in ([5.0, 0, 0], [0, -7.0, 3.0], [100.0, 100.0, 100.0]))

    # CONTROLE NEGATIF : le modele tourne, la gravite laissee en -z
    m, i, a, b = bati(rot([1, 2, 3], 0.6458), np.zeros(3), g_tournee=False)
    m.simule(1.0, 1e-3, tous=10 ** 9)
    faute = float(np.max(np.abs(bout(m, i, a, b) - ref)))

    for k, val in res.items():
        print(f"     {k:9s} : ecart max {val:.2e} m")
    print(f"     controle negatif (gravite non tournee) : {faute:.2e} m — le banc mord")

    for k, val in res.items():
        assert val < 1e-12, (f"invariance {k} violee", val)
    # · sans ceci, un banc d'invariance qui passe ne prouve RIEN
    assert faute > 1e6 * max(res.values()), ("le banc ne mord pas", faute, res)
    return res

def modeles_mal_poses():
    """MODÈLES PLAUSIBLES MAL POSÉS — pas absurdes, JUSTE FAUX.

    `robustesse` couvre ce qu'un noyau fait devant un modèle ABSURDE (h = 0,
    indice hors bornes, liaison d'un corps à lui-même). Ce n'est pas ce qui
    fait échouer un utilisateur nouveau. Ce qui le fait échouer, c'est un
    modèle qui a l'air parfaitement normal et qui est faux : dix cas sondés le
    3 sept., cinq passaient en SILENCE.

    LA DISTINCTION QUI DÉCIDE : refuser, ou publier ?

    Deux sont des erreurs de saisie sans aucun usage légitime, donc REFUSÉES :
    un tenseur d'inertie **non symétrique** (le noyau lisait la moitié de ce
    que l'utilisateur avait écrit et jetait l'autre), et une rotation initiale
    **non orthogonale** (ré-orthonormalisée en douce à chaque pas, donc le
    modèle intégré n'était plus celui écrit).

    La troisième ne peut PAS être un refus, et c'est un banc étranger qui l'a
    montré : `multibarmech` donne à ses barres un tenseur `diag(mL²/12,
    mL²/12, 1)` — un `1.` de remplissage sur l'axe qui ne travaille pas, qui
    **viole l'inégalité triangulaire** des moments principaux. Aucun solide
    réel n'a ce tenseur ; un benchmark publié l'utilise quand même. Refuser
    aurait cassé le banc 20.

    ⇒ `Noyau.controle()` : la question « mon modèle est-il bien posé ? », que
    l'utilisateur nouveau ne sait pas encore poser. Il rend la liste des
    anomalies — inertie hors inégalité triangulaire, contraintes non
    satisfaites, vitesses incompatibles (que le premier pas corrigerait en
    silence), corps qu'aucun élément ne touche — et n'impose rien.

    Ce qui RESTE silencieux, nommé : un pas de temps mille fois trop grand
    rend un résultat lissé sans rien dire. Il n'existe pas de critère
    universel là-dessus ; le pas adaptatif est la réponse du noyau, et elle
    est facultative.
    """
    J = [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3]
    refus, passe = [], []

    def sonde(nom, f, doit_refuser):
        try:
            f()
            (passe if not doit_refuser else refus).append((nom, None))
            return None
        except (ValueError, RuntimeError) as ex:
            (refus if doit_refuser else passe).append((nom, type(ex).__name__))
            return str(ex)

    def non_sym():
        Noyau([0, 0, -9.81]).corps("c", 1.0, [1e-3, 5e-4, 0, 0, 1e-3, 0, 0, 0, 1e-3], [0, 0, 0])

    def non_orth():
        Noyau([0, 0, -9.81]).corps("c", 1.0, J, [0, 0, 0],
                                   rot=[1, 0, 0, 0, 1, 0.01, 0, 0, 1])

    def sym_ok():
        Noyau([0, 0, -9.81]).corps("c", 1.0, J, [0, 0, 0],
                                   rot=[0, -1, 0, 1, 0, 0, 0, 0, 1])

    m1 = sonde("inertie non symetrique", non_sym, True)
    m2 = sonde("rotation non orthogonale", non_orth, True)
    sonde("rotation licite (90 deg)", sym_ok, False)

    # ce que controle() voit, sur un modele qui a l'air normal
    n = Noyau([0.0, 0.0, -9.81])
    a = n.corps("a", 1.0, J, [0, 0, 0])
    b = n.corps("b", 1.0, [1.0, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3], [1, 0, 0], v=[0, 5, 0])
    n.corps("orphelin", 1.0, J, [5, 0, 0])
    n.liaison("l", None, a, pa=[0, 0, 0])
    n.liaison("m", a, b, pa=[1, 0, 0])
    vus = n.controle()
    sain = Noyau([0.0, 0.0, -9.81])
    s = sain.corps("s", 1.0, J, [0, 0, 0])
    sain.liaison("l", None, s, pa=[0, 0, 0])

    print("  22. modeles PLAUSIBLES mal poses — pas absurdes, juste faux")
    print(f"     refuses a la construction : {len(refus)} · "
          + " · ".join(n_ for n_, _ in refus))
    for s_ in vus:
        print(f"     controle() : {s_[:96]}")
    print(f"     modele sain : controle() rend {sain.controle()}")

    assert len(refus) == 2, ("un modele mal pose passe encore", refus, passe)
    assert "SYMÉTRIQUE" in m1 and "orthogonale" in m2, ("message peu utile", m1, m2)
    # · une rotation LICITE doit passer : un refus trop large est pire que rien
    assert not [x for x in passe if x[1]], ("un modele licite est refuse", passe)
    # · controle() voit les trois anomalies, et se tait sur un modele sain
    assert len(vus) == 3, ("controle() manque une anomalie", vus)
    assert any("triangulaire" in x for x in vus), "inegalite triangulaire non vue"
    assert any("vitesses" in x for x in vus), "vitesses incompatibles non vues"
    assert any("aucune liaison, force ou interaction" in x for x in vus), "corps orphelin non vu"
    assert sain.controle() == [], ("controle() crie au loup sur un modele sain",
                                   sain.controle())
    return dict(refus=[n_ for n_, _ in refus], controle=vus)

def _treillis_boucles(ny, nz, lg=1.0):
    """Le treillis de `multibarmech`, agrandi — des BOUCLES FERMEES, structure
    bien moins favorable qu'une chaine."""
    dp = 2 * np.pi / 3
    n = Noyau([0, 0, -9.81])
    jv = [lg * lg / 12, 0, 0, 0, lg * lg / 12, 0, 0, 0, 1.0]
    jh = [lg * lg / 12, 0, 0, 0, 1.0, 0, 0, 0, lg * lg / 12]
    v, h = {}, {}
    for iz in range(nz):
        for iy in range(ny + 1):
            v[(iz, iy)] = n.corps(f"v{iz}_{iy}", 1.0, jv, [0, iy * lg, -iz * lg - lg / 2],
                                  v=[0, dp * (iz * lg + lg / 2), 0], w=[dp, 0, 0])
        for iy in range(1, ny + 1):
            h[(iz, iy)] = n.corps(f"h{iz}_{iy}", 1.0, jh,
                                  [0, (iy - 1) * lg + lg / 2, -(iz + 1) * lg],
                                  v=[0, dp * (iz + 1) * lg, 0])
    k = [0]

    def li(a, b, pa, bt, br):
        k[0] += 1
        n.liaison(f"j{k[0]}", a, b, pa=pa, bloque_t=bt, bloque_r=br)

    for iz in range(nz):
        if iz == 0:
            li(None, v[(0, 0)], [0, 0, 0], [0, 1, 2], [1, 2])
        else:
            li(v[(iz, 0)], v[(iz - 1, 0)], [0, 0, lg / 2], [0, 1, 2], [1, 2])
        for iy in range(1, ny + 1):
            if iz == 0:
                li(None, v[(0, iy)], [0, iy * lg, 0], [1, 2], [])
            else:
                li(v[(iz, iy)], v[(iz - 1, iy)], [0, 0, lg / 2], [1, 2], [])
        li(v[(iz, 0)], h[(iz, 1)], [0, 0, -lg / 2], [0, 1, 2], [1, 2])
        for iy in range(1, ny):
            li(v[(iz, iy)], h[(iz, iy + 1)], [0, 0, -lg / 2], [1, 2], [])
        for iy in range(1, ny + 1):
            li(v[(iz, iy)], h[(iz, iy)], [0, 0, -lg / 2], [1, 2], [])
    return n


def echelle():
    """JUSQU'OÙ ÇA MONTE — mesuré, pas supposé.

    J'ai écrit le 3 sept. qu'un modèle industriel de 10⁴–10⁵ ddl était
    « inatteignable » en algèbre dense. **DEUX croyances dans une phrase, et
    les deux sont fausses.** Paul a relevé la première ; en mesurant, la
    seconde est tombée aussi.

    · **le noyau n'est PAS dense.** `DENSE_MAX = 160` : au-delà, il bascule
      sur un LU CREUX (faer, factorisation symbolique gardée, supernœuds
      parallélisés au-delà de `CREUX_PAR = 2000`). J'avais décrit
      l'architecture sans la lire, et proposé d'écrire un chantier qui
      existait déjà ;
    · l'exposant observé est **~1,07** — quasi-linéaire ;
    · **480 000 inconnues (80 000 corps) à 1,34 s/pas dans 816 Mo.** Il n'y a
      pas de mur là où j'en voyais un.

    CE QUE LA CHASSE A TROUVÉ, en suivant la mesure et non l'intuition. Tout
    le coût venait de **quatre assemblages DENSES de G (m × n)** que rien
    n'obligeait :

    | où | pourquoi c'était inutile | gain à 30 000 inconnues |
    |---|---|---|
    | `tol_phi` (posé le matin même) | ne veut que Φ | 1,3 s **par pas** |
    | `detecte_redondance` | assemblait AVANT de tester son seuil de saut, puis retournait 0 | 1,3 s + 1,7 Go |
    | `acc_init` ×2 | n'en tire que deux produits `G·u` | 2,6 s + 3,4 Go |
    | triplets de `acc_init` | balayait `for j in 0..n` pour 12 termes par ligne | — |

    `phi_seul`, `g_u` et une construction par élément les remplacent.
    **1 727 ms/pas et 4 177 Mo → 41 ms/pas et 96 Mo** : ×42 en temps, ×44 en
    mémoire, sur le même cas.

    ⚠ ET J'AI MESURÉ MA PROPRE RÉGRESSION COMME SI C'ÉTAIT LE SOLVEUR. Le
    premier de ces quatre était à moi, posé le matin même ; j'en avais tiré
    deux bancs et un plan de chantier « solveur creux » — pour un solveur qui
    était déjà creux. **Avant de conclure d'une mesure de performance,
    vérifier ce qu'on a soi-même changé dans la journée.**

    Le chrono manquant est posé au passage : `chronos()` rend maintenant la
    fin de pas et le TOTAL, et le non-instrumenté est tombé de 82 % à 1 %.

    Ce banc garde l'exposant : s'il repassait au-dessus de 2,5, quelque chose
    aurait changé de nature dans le solveur et il faudrait le savoir.
    """
    import time

    J = [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3]

    def chaine(nb):
        n = Noyau([0, 0, -9.81])
        idx = [n.corps(f"c{i}", 1.0, J, [0.2 * i, 0, 0]) for i in range(nb)]
        for i in range(nb - 1):
            n.liaison(f"l{i}", idx[i], idx[i + 1], pa=[0.2, 0, 0], bloque_r=[0, 2])
        n.liaison("bati", None, idx[0], pa=[0, 0, 0], bloque_r=[0, 2])
        return n

    print("  23. echelle — jusqu'ou ca monte, mesure et non suppose")
    pts = []
    for nb in (120, 240, 480):
        n = chaine(nb)
        n.simule(0.005, 1e-3, tous=10 ** 9)          # amorce jetee
        t0 = time.time()
        n.simule(0.045, 1e-3, tous=10 ** 9)
        dt = (time.time() - t0) / 40 * 1000
        inc = 6 * nb + 4 * (nb - 1) + 4
        pts.append((nb, inc, dt, inc ** 2 * 8 / 1e6))
        print(f"     {nb:4d} corps · {inc:5d} inconnues · {dt:6.1f} ms/pas "
              f"(une matrice DENSE aurait pese {inc ** 2 * 8 / 1e6:5.0f} Mo)")
    ex = [float(np.log(pts[i + 1][2] / pts[i][2]) / np.log(pts[i + 1][0] / pts[i][0]))
          for i in range(len(pts) - 1)]
    # LES POUTRES DOIVENT RESTER LINEAIRES. `Poutre::tangente` copiait TOUT le
    # vecteur des corps 24 fois par poutre et par jacobien — O(N^2) sur une
    # structure. Mesure du 3 sept. : 2 000 poutres a 4,6 s/pas et 8 000 a
    # 69 SECONDES par pas, dont 99,5 % dans le jacobien. Deux corps suffisent.
    def _poutres(nb):
        n = Noyau([0, 0, -9.81])
        ix = [n.corps(f"c{i}", 1.0, [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3],
                      [0.2 * i, 0, 0]) for i in range(nb)]
        for i in range(nb - 1):
            n.poutre(f"p{i}", ix[i], ix[i + 1], 1e6, 5e5, 200.0, 400.0)
        n.liaison("b", None, ix[0], pa=[0, 0, 0])
        return n

    pp = []
    for nb in (250, 500, 1000):
        n = _poutres(nb)
        n.simule(0.002, 1e-3, tous=10 ** 9)
        t0 = time.time()
        n.simule(0.007, 1e-3, tous=10 ** 9)
        pp.append((6 * nb + 6, (time.time() - t0) / 5 * 1000))
    ep = float(np.log(pp[-1][1] / pp[0][1]) / np.log(pp[-1][0] / pp[0][0]))
    print(f"     chaine de POUTRES : "
          + " · ".join(f"{a_} inc {b_:.1f} ms" for a_, b_ in pp)
          + f"  exposant {ep:.2f}")

    # LES AUTRES BRIQUES : contacts et aero. RESULTAT NEGATIF, garde comme tel
    # — le balayage de topologies avait sorti deux O(N^2) (tangente de poutre,
    # critere de parallelisation) ; ici il n'en sort aucun, et c'est ce qu'on
    # veut pouvoir affirmer plus tard sans le re-mesurer.
    def _spheres(nb):
        n = Noyau([0, 0, -9.81])
        k = int(np.ceil(nb ** (1 / 3)))
        for i in range(nb):
            a_, b_, c_ = i % k, (i // k) % k, i // (k * k)
            j = n.corps(f"s{i}", 1.0, [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3],
                        [0.15 * a_, 0.15 * b_, 0.15 * c_ + 1.0])
            n.sphere(j, [0, 0, 0], 0.1)
        n.appariement(1e5, 10.0)
        return n

    sp = []
    for nb in (200, 800, 1600):
        n = _spheres(nb)
        n.simule(0.002, 1e-3, tous=10 ** 9)
        t0 = time.time()
        n.simule(0.007, 1e-3, tous=10 ** 9)
        sp.append((nb, (time.time() - t0) / 5 * 1000, n.appariement_stats()[0]))
    es = float(np.log(sp[-1][1] / sp[0][1]) / np.log(sp[-1][0] / sp[0][0]))
    print(f"     spheres en CONTACT : "
          + " · ".join(f"{a_} corps {b_:.1f} ms" for a_, b_, _ in sp)
          + f"  ({sp[-1][2]} paires actives)  exposant {es:.2f}")

    # LE CAS QUI PANIQUAIT — treillis a boucles fermees, 820 corps, 7 480
    # inconnues. faer debordait ses tableaux et PANIQUAIT au travers de PyO3
    # (« range end index 95682 out of range for slice of length 95634 »)
    # parce que la symbolique GARDEE etait reutilisee sur un motif qui avait
    # change. Defaut ANTERIEUR, invisible sous 2 000 inconnues.
    nt = _treillis_boucles(20, 20)
    nt.simule(0.003, 1e-3, tous=10 ** 9)
    print(f"     treillis a boucles fermees, {len(nt.etat()[1])} corps · "
          f"{6 * len(nt.etat()[1]) + len(nt.phi())} inconnues : passe "
          f"(il PANIQUAIT avant le 3 sept.)")

    glob = float(np.log(pts[-1][2] / pts[0][2]) / np.log(pts[-1][0] / pts[0][0]))
    print(f"     exposant local : {' · '.join(f'{e:.2f}' for e in ex)}"
          f"   ·   global {glob:.2f}")
    print("     (une discontinuite reproductible s'observait ici entre 2400 et 4800")
    print("      inconnues, publiee comme un fait de machine — elle a DISPARU avec")
    print("      les assemblages denses. C'etait encore la meme regression.)")

    # · le cout n'est PAS cubique — c'etait ma croyance, et elle etait fausse.
    #   On juge sur l'exposant GLOBAL : les locaux portent la marche de cache,
    #   qui est un fait de machine et non une loi du solveur.
    assert glob < 2.5, ("l'exposant a change de nature", glob, ex, pts)
    # · et le cout MONTE quand meme avec la taille : un exposant nul dirait que
    #   la mesure est dominee par un forfait, pas par l'algebre
    # borne BASSE lâche à dessein : depuis que les assemblages denses sont
    # partis, ces tailles-là coûtent si peu que la mesure est dominée par les
    # forfaits d'allocation, et la gate les mesure sous charge parallèle
    # (observé : 0,73 en gate contre 0,90 à vide). Sous 0,5, ce ne serait plus
    # du bruit mais un vrai forfait dominant.
    assert glob > 0.5, ("mesure dominee par un forfait — allonger la fenetre", glob, ex)
    # · une structure de POUTRES reste lineaire : c'est la regression a ne pas
    #   laisser revenir (elle valait x151 sur 8 000 poutres)
    assert ep < 1.5, ("le jacobien de poutre est redevenu quadratique", pp, ep)
    # · l'appariement de contacts reste lineaire, paires actives comprises
    assert es < 1.5, ("l'appariement de contacts est devenu quadratique", sp, es)
    assert sp[-1][2] > 100, ("les spheres ne se touchent pas — le banc ne juge rien", sp)
    print("     (le noyau bascule en LU CREUX au-dela de DENSE_MAX = 160 inconnues.)")
    print("     mesure a part, apres suppression de quatre assemblages DENSES de G :")
    print("        60 000 inconnues  135 ms/pas   127 Mo")
    print("       480 000 inconnues 1339 ms/pas   816 Mo   exposant 1.07")
    print("     Le meme cas a 30 000 valait 1727 ms et 4177 Mo le matin meme : x42.")
    return dict(points=pts, exposants=ex)

def sparsite():
    """LA STRUCTURE DU PROBLÈME, contre le remplissage réel.

    ⚠ **CE BANC A ÉTÉ ÉCRIT SUR UNE FAUSSE PRÉMISSE, et c'est ce qu'il garde
    de plus utile.** Je le voulais pour chiffrer un chantier « solveur creux »
    — or le noyau EST déjà creux au-delà de `DENSE_MAX = 160`. Ce qu'il mesure
    reste vrai et sert à autre chose : la structure du problème.

    Mesuré sur des modèles réels : le nombre de termes non nuls de la raideur
    croît **linéairement** en n — environ 1 par degré de liberté sur une
    chaîne — et la demi-bande reste petite (8 sur une chaîne, 71 sur le
    treillis 5×5 de `multibarmech`, 4,7 à 15 % sur les deux robots URDF qui
    sont petits). La MATRICE est donc creuse en O(n).

    ET LE REMPLISSAGE EST DÉJÀ BON — mais c'est une TRACE, pas une mesure de
    ce banc, et il faut le dire : une factorisation creuse de la matrice de pas
    [[M+K, Gᵀ], [G, 0]] (SuperLU) donne **6 à 7 non-nuls par ligne avec COLAMD,
    9,3 avec l'ordering NATUREL**, exposant ~1,0 dans les deux cas — mesuré à
    part le 3 septembre, jamais rejoué ici. À 30 000 inconnues le LU tient donc
    dans quelques mégaoctets. Ce qui consomme 4,2 Go est ailleurs, et `echelle`
    dit où on en est de cette chasse.

    ⚠ **RECTIFIÉ LE 5 SEPTEMBRE — « le LU n'est PAS le point dur » était un
    argument de MÉMOIRE publié comme un verdict général.** Il l'est resté
    jusqu'à ce qu'un chrono soit lu : la RÉSOLUTION pèse **32 / 65 / 52 %** du
    pas à 18 / 180 / 540 inconnues (`refus_gardes` la publie). Les deux
    affirmations sont vraies et elles ne parlent pas de la même chose ; les
    confondre ferait écarter à tort une sous-structuration (note d'architecture
    §6, Peiret 2019) au motif d'un remplissage. Ce qui écarte Schur est
    ailleurs, et c'est mesuré : **zéro repli SVD**.
    """
    J = [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3]

    def chaine(nb):
        n = Noyau([0, 0, -9.81])
        idx = [n.corps(f"c{i}", 1.0, J, [0.2 * i, 0, 0]) for i in range(nb)]
        for i in range(nb - 1):
            n.liaison(f"l{i}", idx[i], idx[i + 1], pa=[0.2, 0, 0], bloque_r=[0, 2])
        n.liaison("b", None, idx[0], pa=[0, 0, 0], bloque_r=[0, 2])
        return n

    print("  24. sparsite — la structure, et l'ecart au remplissage reel")
    pts = []
    for nb in (20, 40, 80):
        n = chaine(nb)
        k, _, _ = n.k_m_z()
        N = 6 * nb
        K = np.asarray(k).reshape(N, N)
        nz = int(np.count_nonzero(K))
        ij = np.argwhere(K != 0)
        bande = int(np.max(np.abs(ij[:, 0] - ij[:, 1]))) if len(ij) else 0
        pts.append((nb, N, nz, nz / N, bande))
        print(f"     {nb:3d} corps · n={N:4d} · {nz:6d} non nuls sur {N * N:7d} "
              f"({100 * nz / (N * N):5.2f} %) · {nz / N:4.1f} par ddl · demi-bande {bande}")
    ratio = [p[3] for p in pts]

    # · le nombre de non nuls croit LINEAIREMENT : le ratio par ddl est borne
    assert max(ratio) < 2 * min(ratio), ("les non nuls ne croissent pas en O(n)", pts)
    # · et la demi-bande ne grandit pas avec le modele
    assert pts[-1][4] <= 2 * pts[0][4], ("la bande grandit avec n", pts)
    g = pts[-1][1] ** 2 / pts[-1][2]
    print(f"     la MATRICE est creuse en O(n) : facteur {g:.0f} deja a n={pts[-1][1]}")
    print("     TRACE du 3 sept., mesure faite A PART et que ce banc NE REJOUE PAS :")
    print("     SuperLU sur [[M+K,Gt],[G,0]] remplit a 6-7 non nuls par ligne avec")
    print("     COLAMD, 9.3 avec l'ordering NATUREL, exposant ~1.0 dans les deux.")
    print("     ⇒ le LU ne coute pas de MEMOIRE : quelques Mo a 30 000 inconnues.")
    print("       Il ne dit rien du TEMPS — voir `refus_gardes` (t_sol 32-65 % du pas).")
    return dict(points=pts, gain=g)

def non_holonome():
    """ROULEMENT SANS GLISSEMENT — la brique qui manquait pour qu'un véhicule
    entre dans ce noyau.

    Jusqu'au 3 sept. vinkulum ne connaissait que des contraintes **holonomes**,
    `Φ(q) = 0`. Or la condition de roulement ne s'écrit PAS comme une fonction
    de la position : une roue peut atteindre n'importe quelle pose en roulant.
    Elle s'écrit sur les VITESSES — `G·u = 0`, la vitesse du point de contact
    est nulle — et c'est une classe entière de mécanismes (véhicules, robots
    mobiles, patins, galets) qui restait dehors.

    Le jacobien est le MÊME que celui de la liaison holonome correspondante ;
    seul le NIVEAU auquel on l'impose change. D'où `liaison(..., nh=True)`
    plutôt qu'un élément séparé.

    TROIS CHOSES QU'IL A FALLU CORRIGER, chacune mesurée :

    · **la mise à l'échelle**. Les lignes de position sont divisées par βh²,
      ce qui rend leur jacobien `G·c`. Une ligne de vitesse divisée par la
      même chose donne un gradient **2 000 fois trop petit** : Newton diverge
      (résidu 1e38 mesuré). Elle se divise par γh ;
    · **le point de contact est GÉOMÉTRIQUE, pas matériel**. En gardant le
      `pb` figé dans la roue, on suit un point de la jante qui s'éloigne du
      sol : la roue décélérait et repartait en arrière (−0,28 m au lieu de
      +1,20). Le bras vers B se recalcule à chaque évaluation ;
    · **le contrôle `‖Φ‖` doit SAUTER ces lignes**. Une contrainte non
      holonome n'a pas de Φ à annuler — celui de sa liaison mesure la position
      du point de contact, qui dérive légitimement puisque la roue avance.
      Sans ça, le contrôle posé le matin même refusait le roulement.
    """
    rw, ms, om, tf = 0.3, 2.0, 4.0, 1.0
    jd, ja = ms * rw * rw / 4, ms * rw * rw / 2

    def monte(nh):
        n = Noyau([0.0, 0.0, -9.81])
        w = n.corps("roue", ms, [jd, 0, 0, 0, ja, 0, 0, 0, jd], [0.0, 0.0, rw],
                    v=[om * rw, 0.0, 0.0], w=[0.0, om, 0.0])
        n.liaison("essieu", None, w, pa=[0, 0, rw], bloque_t=[1, 2], bloque_r=[0, 2])
        if nh:
            n.liaison("contact", None, w, pa=[0.0, 0.0, 0.0],
                      bloque_t=[0], bloque_r=[], nh=True)
        return n, w

    n, w = monte(True)
    phi0 = max(abs(x) for x in n.phi())
    gu0 = max(abs(x) for x in n.phi_dot())
    n.simule(tf, 1e-3, tous=10 ** 9)
    x = n.etat()[1][w][0]
    v = n.etat()[3][w][0]
    e = abs(x - om * rw * tf)

    # CONTROLE NEGATIF : sans la contrainte, rien ne lie la rotation a
    # l'avance. La roue tourne sur place et ne parcourt que ce que son elan
    # initial lui donne — la contrainte doit changer la trajectoire.
    n2, w2 = monte(False)
    n2.pose_etat(0.0, n2.etat()[1], n2.etat()[2],
                 [[0.0, 0.0, 0.0]], [[0.0, om, 0.0]])
    n2.simule(tf, 1e-3, tous=10 ** 9)
    x2 = n2.etat()[1][w2][0]

    print("  25. non holonome — roulement sans glissement")
    print(f"     etat initial admissible : |Φ| {phi0:.1e} · |G·u| {gu0:.1e}")
    print(f"     apres {tf:.0f} s : x = {x:+.6f} m · theorie ωRt = {om * rw * tf:+.6f} · "
          f"ecart {e:.1e}")
    print(f"     vitesse finale {v:+.6f} m/s (theorie {om * rw:+.6f})")
    print(f"     sans la contrainte, meme rotation : x = {x2:+.6f} m (controle negatif)")

    assert phi0 < 1e-12 and gu0 < 1e-12, ("etat initial non admissible", phi0, gu0)
    # · la roue avance EXACTEMENT de ωRt : c'est la definition du roulement
    assert e < 1e-9, ("le roulement ne tient pas la relation ωRt", x, e)
    assert abs(v - om * rw) < 1e-9, ("la vitesse n'est pas maintenue", v)
    # · et sans la contrainte, ce n'est pas le meme mouvement
    assert abs(x2 - om * rw * tf) > 0.1, ("la contrainte ne change rien", x2)
    return dict(x=x, ecart=e, sans=x2)

def coulomb_pente():
    """FROTTEMENT DE COULOMB — les deux régimes, et le seuil qui les sépare.

    Une sphère lâchée sur une pente. Selon μ, deux physiques différentes, et
    **chacune a sa solution analytique** :

    ```text
        μ < μ_c   glissement :  a = g(sin α − μ cos α)
        μ > μ_c   roulement  :  a = g sin α / (1 + I/mR²)   — indépendant de μ
        μ_c = tan α · I/(I + mR²)
    ```

    Ce banc est né d'une erreur de lecture, et c'est ce qui le rend bon. J'ai
    d'abord cru à un frottement SATURÉ : μ = 0,20 et μ = 0,60 donnaient
    exactement la même trajectoire. C'était le régime de ROULEMENT — au-delà
    de μ_c l'accélération ne dépend plus de μ, c'est la physique qui le dit.
    Le noyau avait raison, la lecture non.

    ⇒ Le banc juge donc les DEUX régimes contre leur formule, et le seuil qui
    les sépare. C'est plus fort qu'un seul point : une erreur de signe ou
    d'échelle sur μ casserait le régime glissant, une erreur sur le bras de
    levier casserait le régime roulant.
    """
    al = np.radians(20.0)
    nrm = np.array([np.sin(al), 0.0, np.cos(al)])
    tan = np.array([np.cos(al), 0.0, -np.sin(al)])       # descente
    rs, ms, tf = 0.05, 1.0, 1.5
    inert = 1e-3                                          # = 2/5 m R² : sphère pleine
    mu_c = np.tan(al) * inert / (inert + ms * rs ** 2)
    a_roule = 9.81 * np.sin(al) / (1 + inert / (ms * rs ** 2))

    def glisse(mu):
        n = Noyau([0.0, 0.0, -9.81])
        c = n.corps("bille", ms, [inert, 0, 0, 0, inert, 0, 0, 0, inert], list(nrm * rs))
        n.contact("sol", c, [0, 0, 0], rs, normale=list(nrm), origine=[0, 0, 0],
                  k=1e6, expo=1.5, c=50.0, mu=mu, v_eps=1e-4)
        n.simule(tf, 2e-4, tous=10 ** 9)
        p = np.array(n.etat()[1][c]) - nrm * rs
        return float(np.dot(p, tan)), float(np.dot(np.array(n.etat()[3][c]), tan))

    print("  26. Coulomb sur une pente — glissement, roulement, et le seuil")
    print(f"     pente {np.degrees(al):.0f}° · tan α = {np.tan(al):.4f} · "
          f"μ_c = {mu_c:.4f} · a_roulement = {a_roule:.4f} m/s²")
    res = {}
    for mu in (0.05, 0.20, 0.60):
        x, v = glisse(mu)
        a_th = (9.81 * (np.sin(al) - mu * np.cos(al)) if mu < mu_c else a_roule)
        a_me = v / tf
        res[mu] = (x, v, a_me, a_th)
        reg = "glisse " if mu < mu_c else "roule  "
        print(f"     μ={mu:.2f} {reg}: a mesure {a_me:+.4f} · theorie {a_th:+.4f} · "
              f"ecart {100 * abs(a_me - a_th) / abs(a_th):.2f} %")

    # · regime GLISSANT : la formule de Coulomb, au pourcent
    assert abs(res[0.05][2] - res[0.05][3]) / res[0.05][3] < 0.02, ("glissement", res[0.05])
    # · regime ROULANT : independant de mu, et egal a g sin/(1+I/mR²)
    for mu in (0.20, 0.60):
        assert abs(res[mu][2] - a_roule) / a_roule < 0.02, ("roulement", mu, res[mu])
    # · et les deux regimes roulants coincident : c'est CE qui prouve que
    #   l'independance en mu n'est pas une saturation numerique
    assert abs(res[0.20][0] - res[0.60][0]) < 1e-3, ("les deux mu ne coincident pas", res)
    print(f"     les deux μ roulants coincident a {abs(res[0.20][0] - res[0.60][0]):.1e} m — "
          f"c'est la physique, pas une saturation")
    return dict(mu_c=mu_c, a_roule=a_roule, res={k: v[2] for k, v in res.items()})

def accords():
    """DEUX BRIQUES POUR LA MÊME PHYSIQUE DOIVENT S'ACCORDER.

    Chaque banc précédent juge UNE brique contre une référence. Celui-ci juge
    les briques ENTRE ELLES, là où deux chemins indépendants du noyau
    décrivent le même phénomène. C'est ce qui reste vérifiable sans utilisateur
    tiers, et c'est là que dorment les incohérences qu'un tiers rencontrerait :
    un code peut être juste sur chaque brique prise à part et se contredire dès
    qu'on en combine deux.

    · **le roulement** s'obtient de DEUX façons sans rapport l'une avec
      l'autre — par contact pénalisé plus frottement de Coulomb, ou par
      contrainte NON HOLONOME. La première est une force, la seconde un
      multiplicateur ; elles ne partagent pas une ligne de code ;
    · **une liaison rigide** est la limite d'une poutre dont on augmente la
      raideur. La liaison est algébrique (un multiplicateur), la poutre est
      élastique (une force interne). Elles doivent converger l'une vers
      l'autre, et la vitesse de convergence dit que la limite est la bonne.
    """
    al = np.radians(15.0)
    nrm = np.array([np.sin(al), 0.0, np.cos(al)])
    tng = np.array([np.cos(al), 0.0, -np.sin(al)])
    rs, ms, tf = 0.05, 1.0, 1.0
    inert = 0.4 * ms * rs * rs                       # sphère pleine
    jj = [inert, 0, 0, 0, inert, 0, 0, 0, inert]

    def par_frottement():
        n = Noyau([0.0, 0.0, -9.81])
        c = n.corps("b", ms, jj, list(nrm * rs))
        n.contact("sol", c, [0, 0, 0], rs, normale=list(nrm), origine=[0, 0, 0],
                  k=1e7, expo=1.5, c=50.0, mu=0.8, v_eps=1e-5)
        n.simule(tf, 2e-4, tous=10 ** 9)
        return float(np.dot(np.array(n.etat()[1][c]) - nrm * rs, tng))

    def par_contrainte():
        n = Noyau([0.0, 0.0, -9.81])
        c = n.corps("b", ms, jj, list(nrm * rs))
        n.liaison("plan", None, c, pa=list(nrm * rs),
                  ra=list(np.column_stack([nrm, tng, np.cross(nrm, tng)]).flatten()),
                  bloque_t=[0, 2], bloque_r=[0, 1])
        n.liaison("roule", None, c, pa=[0, 0, 0],
                  ra=list(np.column_stack([tng, nrm, np.cross(tng, nrm)]).flatten()),
                  bloque_t=[0], bloque_r=[], nh=True)
        n.simule(tf, 2e-4, tous=10 ** 9)
        return float(np.dot(np.array(n.etat()[1][c]) - nrm * rs, tng))

    xf, xc = par_frottement(), par_contrainte()
    th = 0.5 * 9.81 * np.sin(al) / (1 + inert / (ms * rs * rs)) * tf * tf

    jl = [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3]

    def rigide():
        n = Noyau([0.0, 0.0, -9.81])
        a = n.corps("a", 1.0, jl, [0, 0, 0])
        b = n.corps("b", 1.0, jl, [0.5, 0, 0])
        n.liaison("enc", None, a, pa=[0, 0, 0])
        n.liaison("rig", a, b, pa=[0.5, 0, 0])
        n.simule(0.3, 1e-4, tous=10 ** 9)
        return np.array(n.etat()[1][b])

    def souple(ech):
        n = Noyau([0.0, 0.0, -9.81])
        a = n.corps("a", 1.0, jl, [0, 0, 0])
        b = n.corps("b", 1.0, jl, [0.5, 0, 0])
        n.liaison("enc", None, a, pa=[0, 0, 0])
        n.poutre("p", a, b, 1e6 * ech, 5e5 * ech, 1e4 * ech, 2e4 * ech)
        n.simule(0.3, 1e-4, tous=10 ** 9)
        return np.array(n.etat()[1][b])

    pl = rigide()
    conv = [(e, float(np.max(np.abs(souple(e) - pl)))) for e in (1, 10, 100, 1000)]

    # 3. STATIQUE a charge croissante : la continuation doit se faire SEULE
    def console(fz, iters=100, tol=1e-10):
        n = Noyau([0.0, 0.0, 0.0])
        ix = [n.corps(f"c{i}", 0.5, jl, [i / 6.0, 0, 0]) for i in range(7)]
        for i in range(6):
            n.poutre(f"p{i}", ix[i], ix[i + 1], 1e6, 5e5, 50.0, 100.0)
        n.liaison("enc", None, ix[0], pa=[0, 0, 0])
        n.effort(ix[-1], [0.0, 0.0, -fz], [0.0, 0.0, 0.0])
        r, it = n.statique(tol=tol, iters=iters)
        return float(n.etat()[1][ix[-1]][2]), (it // 1000 if it >= 1000 else 1)

    stat = [(f,) + console(f) for f in (0.5, 5.0, 20.0)]
    # Le nouveau mérite peut résoudre 20 N en un palier : ce progrès ne
    # doit pas faire échouer le banc. Un budget Newton court exerce encore
    # la continuation et doit retrouver le même équilibre.
    z_cont, p_cont = console(20.0, iters=8, tol=1e-8)

    # 4. MODAL vs TEMPOREL : la frequence propre doit etre celle qu'on mesure.
    # Deux chemins sans rapport — un probleme aux valeurs propres sur le
    # systeme contraint LINEARISE, et une integration temporelle NON lineaire.
    # On excite exactement le premier mode : le mouvement devient
    # mono-frequentiel et le pic FFT est net. (Une perturbation quelconque
    # excite plusieurs modes et la mesure ne veut plus rien dire — essaye,
    # elle donnait 285 % d'ecart.)
    nm = Noyau([0.0, 0.0, 0.0])
    im = [nm.corps(f"m{i}", 0.5, jl, [i / 6.0, 0, 0]) for i in range(7)]
    for i in range(6):
        nm.poutre(f"q{i}", im[i], im[i + 1], 1e6, 5e5, 50.0, 100.0)
    nm.liaison("enc", None, im[0], pa=[0, 0, 0])
    md = nm.modes(4)
    f1, v1 = md[0][0], np.asarray(md[0][1])
    st = nm.etat()
    amp = 1e-3 / max(abs(v1).max(), 1e-30) * 2 * np.pi * f1
    vit = [list(amp * v1[6 * i:6 * i + 3]) for i in range(len(im))]
    omg = [list(amp * v1[6 * i + 3:6 * i + 6]) for i in range(len(im))]
    nm.pose_etat(0.0, st[1], st[2], vit, omg)
    hh = nm.simule(4.0 / f1, 2e-4, tous=1)
    tt = np.array([fr[0] for fr in hh])
    zz = np.array([fr[1][im[-1]][2] for fr in hh])
    zz = zz - zz.mean()
    sp = np.abs(np.fft.rfft(zz * np.hanning(len(zz))))
    frq = np.fft.rfftfreq(len(zz), tt[1] - tt[0])
    fpic = float(frq[1:][np.argmax(sp[1:])])

    # 5. SUPERELEMENT vs MAILLAGE COMPLET, en STATIQUE. La condensation de
    # Guyan est EXACTE pour la raideur ; c'est donc la seule ou l'accord doit
    # etre parfait, et c'est ce qui en fait un bon controle. En DYNAMIQUE elle
    # ne l'est pas (la masse reste concentree) — `coque.demo` mesure ce
    # residu-la, ~2 % sur f1, et le declare.
    ne, lg, fz = 6, 1.0, 2.0
    mm = Noyau([0.0, 0.0, 0.0])
    jxx = [mm.corps(f"c{i}", 0.5, jl, [i * lg / ne, 0, 0]) for i in range(ne + 1)]
    for i in range(ne):
        mm.poutre(f"p{i}", jxx[i], jxx[i + 1], 1e6, 5e5, 50.0, 100.0)
    kk, _, _ = mm.k_m_z()
    nn = 6 * (ne + 1)
    kmat = np.asarray(kk).reshape(nn, nn)
    gard = list(range(6)) + list(range(6 * ne, 6 * ne + 6))
    inte = [x for x in range(nn) if x not in gard]
    kr = (kmat[np.ix_(gard, gard)]
          - kmat[np.ix_(gard, inte)]
          @ np.linalg.solve(kmat[np.ix_(inte, inte)], kmat[np.ix_(inte, gard)]))

    nc = Noyau([0.0, 0.0, 0.0])
    ic = [nc.corps(f"c{i}", 0.5, jl, [i * lg / ne, 0, 0]) for i in range(ne + 1)]
    for i in range(ne):
        nc.poutre(f"p{i}", ic[i], ic[i + 1], 1e6, 5e5, 50.0, 100.0)
    nc.liaison("enc", None, ic[0], pa=[0, 0, 0])
    nc.effort(ic[-1], [0.0, 0.0, -fz], [0.0, 0.0, 0.0])
    nc.statique(tol=1e-10, iters=100)
    zc = float(nc.etat()[1][ic[-1]][2])

    ns = Noyau([0.0, 0.0, 0.0])
    sa = ns.corps("a", 0.5 * ne, jl, [0, 0, 0])
    sb = ns.corps("b", 0.5, jl, [lg, 0, 0])
    ns.superelement("se", [sa, sb], list(kr.flatten()), 0.0, 0.0)
    ns.liaison("enc", None, sa, pa=[0, 0, 0])
    ns.effort(sb, [0.0, 0.0, -fz], [0.0, 0.0, 0.0])
    ns.statique(tol=1e-10, iters=100)
    zs = float(ns.etat()[1][sb][2])

    print("  27. accords — deux briques pour la meme physique")
    print(f"     roulement : par FROTTEMENT {xf:+.6f} m · par CONTRAINTE {xc:+.6f} m · "
          f"theorie {th:+.6f}")
    print(f"                 accord entre briques {abs(xf - xc):.1e} m "
          f"(une force contre un multiplicateur)")
    print("     liaison RIGIDE = limite d'une poutre raidie : "
          + " · ".join(f"x{e} {d:.1e}" for e, d in conv))
    print(f"     modal vs temporel : f1 propre {f1:.4f} Hz · pic FFT de l'integration "
          f"{fpic:.4f} Hz · ecart {100 * abs(fpic - f1) / f1:.2f} %")
    print(f"     superelement (Guyan {nn}→{len(gard)} ddl) vs maillage complet, en "
          f"statique : {zs:+.8f} contre {zc:+.8f} m · {100 * abs(zc - zs) / abs(zc):.3f} %")
    print("     statique a charge croissante (continuation AUTOMATIQUE) : "
          + " · ".join(f"{f:.1f} N {z:+.4f} m en {p} palier(s)" for f, z, p in stat))

    # · les deux chemins du roulement donnent le meme resultat, et le bon
    assert abs(xf - xc) < 1e-4, ("frottement et contrainte se contredisent", xf, xc)
    assert abs(xc - th) / th < 1e-4, ("la contrainte ne rend pas la theorie", xc, th)
    # · la poutre raidie CONVERGE vers la liaison : d'un facteur au moins 100
    #   sur trois decades de raideur
    assert conv[-1][1] < conv[0][1] / 100, ("la poutre ne converge pas vers la liaison", conv)
    # · la statique subdivise SEULE quand Newton stagne. Sans continuation elle
    #   tenait jusqu'a 0,5 N et echouait a 2 N sur cette meme console —
    #   l'utilisateur devait deviner qu'il fallait monter la charge par paliers.
    assert p_cont > 1, ("la continuation ne s'active pas avec le budget court", p_cont)
    assert abs(z_cont - stat[-1][1]) < 1e-8, ("la continuation change l'équilibre", z_cont, stat)
    # · et la fleche reste proportionnelle a la charge tant qu'on est petit
    r1, r2 = stat[0][1] / 0.5, stat[1][1] / 5.0
    assert abs(r1 - r2) / abs(r1) < 0.02, ("la reponse n'est plus lineaire a 5 N", stat)
    # · l'analyse MODALE et l'INTEGRATION donnent la meme frequence. Le pas de
    #   la FFT vaut ici f1/4 : un ecart d'un bin serait 25 %, donc 2 % est un
    #   accord franc et pas une coincidence de resolution.
    assert abs(fpic - f1) / f1 < 0.02, ("modal et temporel se contredisent", f1, fpic)
    # · Guyan est EXACT pour la raideur : l'accord doit etre a la precision de
    #   la resolution, pas « proche ». Une erreur de partition ou de signe dans
    #   la condensation se verrait tout de suite ici.
    assert abs(zc - zs) / abs(zc) < 1e-3, ("le superelement ne rend pas le maillage", zc, zs)
    return dict(roulement=(xf, xc, th), convergence=conv, statique=stat,
                modal=(f1, fpic), superelement=(zc, zs))

def _chaine_tournante(nb=2, w=300.0):
    """Le témoin des gardes de refus — une chaîne de corps EN ROTATION.

    La rotation est le point, pas un décor : sans elle le terme gyroscopique
    est nul, la tangente d'amortissement redevient symétrique, et la ligne
    LDLᵀ ne prouverait plus rien.
    """
    J = [1e-3, 0, 0, 0, 2e-3, 0, 0, 0, 3e-3]
    N = Noyau([0.0, 0.0, -9.81])
    idx = [N.corps(f"c{i}", 1.0, J, [0.2 * i, 0, 0], w=[0, 0, w]) for i in range(nb)]
    N.liaison("b", None, idx[0], pa=[0, 0, 0], bloque_r=[0, 1])
    for i in range(nb - 1):
        N.liaison(f"l{i}", idx[i], idx[i + 1], pa=[0.2, 0, 0], bloque_r=[0, 1])
    return N


def refus_gardes():
    """CE QU'ON A REFUSÉ — la RAISON du refus, asserée et non recopiée.

    Écrit le 5 septembre 2026, après avoir confronté au noyau une note
    d'architecture (« Avancées mathématiques pour la prochaine génération de
    solveurs multicorps industriels »). Sur une vingtaine de recommandations,
    la moitié était déjà là et quatre étaient refusées par des mesures
    ANTÉRIEURES. Or le dépôt a un idiome pour ça — l'assert inversé, qui casse
    quand la raison d'un refus cesse d'être vraie (`bancs.adaptatif`,
    `bancs.py:722`, `verification.py:864`, `campagne.py:208`) — et **cinq
    refus n'en avaient aucun**. Un refus mesuré une fois puis jamais rejoué
    n'est plus une mesure : c'est une croyance datée.

    Ce contrôle ne réimplémente rien. Il rejoue la PRÉMISSE de chaque refus,
    sur un modèle témoin de deux corps tournants, et casse si elle tombe.

    ── 1. LDLᵀ (note §6 : « factorisation creuse bloc-structurée, LDLᵀ »)

    La recommandation suppose un point-selle SYMÉTRIQUE. Le nôtre ne l'est
    pas : la tangente d'amortissement porte le terme gyroscopique
    `[ω]×J_s − [J_sω]×` (`src/lib.rs`, `amortissement`), qui n'est pas égal à
    son transposé — mesuré ici à ~1,5 en relatif, c'est-à-dire de l'ordre de
    la matrice elle-même, pas un résidu. Le LU n'est donc pas un choix de
    confort, c'est le seul choix correct.

    Un SECOND motif existe, nommé et non mesuré ici parce que le jacobien de
    pas assemblé n'est pas exposé : ce jacobien porte `Gᵀ` en (1,2) et
    `G·J_l·c` en (2,1) (`src/tangent.rs`, en-tête) — asymétrique par `J_l`
    même si C l'était.

    ── 2. Précision mixte (note §6 : « simple précision + raffinement itératif »)

    Le noyau exige de Newton une tolérance que f32 ne sait pas atteindre. On
    ne recopie pas cette tolérance : on la LIT sur la signature de `simule`.
    Et ce refus n'est pas théorique — le défaut du 5 septembre (le vol S2 qui
    cassait à 50,84 s) était un problème d'**ulp** sur un angle déroulé, en
    double précision. En simple, ce budget n'existe pas.

    ── 3. Complément de Schur / sous-structuration (note §6, Peiret 2019)

    Il vise les assemblages localement RAIDES dont le conditionnement dégrade
    l'ensemble. Le signal existe déjà dans le noyau et personne ne le
    regardait : `stats()` compte les **replis SVD**, c'est-à-dire les pas où
    le LU n'a pas tenu (`Modele::resout` bascule en moindres carrés quand le
    résidu de la résolution dépasse son seuil). Tant qu'il vaut zéro, le
    conditionnement n'est pas le poste dur et Schur répare un mal absent.

    ⚠ Le témoin d'ici est un JOUET — deux corps, aucune raideur locale — et
    Schur vise précisément les assemblages localement raides. La mesure qui
    compte est faite AILLEURS, sur le cas que la recommandation vise : le vol
    S2 de FRELON (`cad/vinkulum_s2.py demo`), ligne d'engrenage à l'échelle
    1/(βh²) = 4e8, servos PD, aéro — **0 repli SVD sur 60 000 pas**, asserté
    là-bas. Ici on garde le cas simple, dans la gate ; là-bas le cas dur, hors
    gate.

    ── 4. Cayley (note §1.3) — PUBLIÉ, PAS ASSERTÉ SUR LE BON CHIFFRE, ET
    C'EST UNE RÉTRACTATION

    La note dit que l'application de Cayley réduit « le coût par itération de
    Newton par rapport à la formulation exponentielle ». J'ai d'abord refusé
    en citant « 1,07 jacobien par pas » — chiffre juste, mais qui décrit le
    **vol S2 de FRELON** et non le noyau : sur des chaînes tournantes il vaut
    1,44 · 2,10 · 4,88 à 18 · 180 · 540 inconnues. Et le jacobien pèse 26 à
    42 % du pas, ce qui n'est pas rien.

    (Le vol S2 lui-même rend 1,05 jac/pas, mesuré le 5 sept. : le chiffre
    était juste pour LUI, et faux comme propriété du noyau.)

    Le bon chiffre n'existait pas, alors il a été posé (`chronos()` rend
    désormais `t_exp` en sixième) : la part du jacobien passée à construire la
    tangente de exp. Elle vaut **0,4 à 1,0 %**, et elle DÉCROÎT avec la taille
    — parce que `jac_gauche` est appelée une fois par CORPS et par jacobien,
    pas par élément ni par itération de Newton. La prémisse de la note ne
    s'applique pas à cette implémentation.

    `t_jac/t_pas` et `jac/pas` sont publiés à côté SANS assert : ils varient
    avec la taille du modèle, donc les borner ici fabriquerait un cliquet.

    ── 5. SE(3) (refusé le 4 septembre, `docs/JOURNAL.md`)

    Le refus s'appuie sur un fait mesuré ailleurs : « l'invariance par
    changement de repère est déjà mesurée à 1e-14 par `verification` ». Rien
    ne reliait le refus à ce contrôle — si `invariances` disparaissait, le
    refus perdrait sa base en silence. Sentinelle d'oubli, patron des
    sentinelles de FRELON : ce qu'il faut interdire n'est pas une valeur
    fausse, c'est l'OUBLI.

    ── CE QUE CE CONTRÔLE NE FAIT PAS, dit tout de suite

    Il ne juge PAS les refus dont la mesure n'est plus rejouable parce que
    l'option a été retirée du code — `simule(projection=True)` (ordre 2,00 →
    ~0,8) est une TRACE de docstring, pas un contrôle, et la ressusciter pour
    la mesurer coûterait plus que ce qu'elle garde. Il ne juge pas non plus
    les refus déjà gardés ailleurs : le pas adaptatif (`bancs.adaptatif`), la
    métrique de masse pour GGL (`bancs.py`), ρ∞ haut (`bancs.py`), l'exposant
    d'échelle (`echelle`), les λ exacts des contacts multiples qui rendent
    SAP/convexe sans objet (`contact.multiples`).
    """
    import inspect

    print("  28. refus gardés — la raison d'un refus, rejouée")
    out = {}

    # ── 1. LDLᵀ : la tangente d'amortissement est-elle encore asymétrique ?
    N = _chaine_tournante()
    k, c, m, z, g = (np.asarray(x, float) for x in N.k_c_m_z())
    asym = float(np.abs(c - c.T).max())
    ech = float(np.abs(c).max())
    rel_c = asym / max(ech, 1e-30)
    sym_m = float(np.abs(m - m.T).max())
    out["asym_c"] = rel_c
    print(f"     LDLᵀ : asym(C)/‖C‖ = {rel_c:.3f} (gyroscopique) · asym(M) = {sym_m:.1e}")
    assert rel_c > 0.1, (
        "la tangente d'amortissement est redevenue symétrique — soit le terme "
        "gyroscopique a disparu (physique fausse), soit LDLᵀ devient licite : "
        "rejouer le refus de la note §6",
        rel_c,
    )

    # ── 2. précision mixte : la tolérance se LIT, elle ne se recopie pas
    tol = inspect.signature(Noyau.simule).parameters["tol"].default
    eps32 = 2.0 ** -24
    out["tol_defaut"], out["eps_f32"] = float(tol), eps32
    print(f"     précision mixte : tol par défaut {tol:.0e} contre eps(f32) {eps32:.2e} "
          f"— {eps32 / tol:.0e}× trop grossier")
    assert tol < eps32, (
        "la tolérance par défaut est passée au-dessus de l'epsilon f32 : la "
        "précision mixte redevient discutable, rejouer le refus de la note §6",
        tol,
    )

    # ── 3 et 4 : un seul vol du témoin
    h, ns = 1e-4, 200
    N.simule(ns * h, h, tous=10 ** 9)
    t_jac, t_sol, t_res, t_fin, t_pas, t_exp = N.chronos()
    _, svd, jac, _ = N.stats()
    out["svd"] = svd
    print(f"     Schur : replis SVD sur {ns} pas = {svd} "
          f"(le LU tient, le conditionnement n'est pas le poste dur)")
    assert svd == 0, (
        "le LU ne tient plus : le conditionnement devient le poste dur et le "
        "complément de Schur (note §6, Peiret 2019) redevient un candidat",
        svd,
    )

    part_exp = t_exp / max(t_jac, 1e-30)
    out.update(part_exp=part_exp, part_jac=t_jac / max(t_pas, 1e-30),
               part_sol=t_sol / max(t_pas, 1e-30), jac_par_pas=jac / ns)
    print(f"     Cayley : t_exp/t_jac = {100 * part_exp:.2f} % "
          f"(publié, non asserté : t_jac/t_pas {100 * out['part_jac']:.0f} % · "
          f"t_sol/t_pas {100 * out['part_sol']:.0f} % · {out['jac_par_pas']:.2f} jac/pas)")
    assert part_exp < 0.05, (
        "la tangente de exp est passée à plus de 5 % du jacobien — elle a dû "
        "entrer dans une boucle interne : rejouer le refus de Cayley (note §1.3)",
        part_exp,
    )

    # ── 5. SE(3) : le contrôle sur lequel le refus s'appuie existe-t-il encore ?
    appui = globals().get("invariances")
    print("     SE(3) : le refus s'appuie sur `invariances` (repère à 1e-14) — "
          f"{'présent' if callable(appui) else 'ABSENT'}")
    assert callable(appui), (
        "le refus de SE(3) (JOURNAL du 4 sept.) cite l'invariance de repère "
        "mesurée par `verification.invariances` — ce contrôle a disparu, le "
        "refus n'a plus de base"
    )

    print("     ⇒ cinq refus portent leur garde ; ils tombent si leur raison tombe.")
    return out


def robustesse():
    """CE QUE LE NOYAU FAIT DEVANT UN MODÈLE ABSURDE — un solveur mature ne
    plante pas, et surtout n'accepte pas en SILENCE.

    Ce banc est né d'un sondage (3 sept.) qui a trouvé trois défauts en dix
    minutes, chacun d'une nature différente :

    · `h = 0` ne faisait pas AVANCER t, donc la boucle tournait à l'infini —
      pas une erreur, un BLOCAGE : aucune trace, aucun diagnostic, un
      processus à tuer. Le pire des trois ;
    · un indice de corps hors bornes faisait PANIQUER Rust au travers de PyO3
      (`index out of bounds`), au lieu de lever une exception Python ;
    · `liaison(a, a)` — et `liaison(bâti, bâti)` — était ACCEPTÉE en silence,
      en produisant une contrainte identiquement nulle. Le corps tombait
      librement sous une liaison qu'on croyait avoir posée. Silencieux, donc
      pire qu'un plantage.

    On vérifie les DEUX sens : ce qui doit être refusé l'est, et ce qui est
    licite passe encore.
    """
    J = [1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4]
    nan = float("nan")
    PLAQUE = ([-1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1, 0], [0, 1, 2, 0, 2, 3])

    def mk():
        n = Noyau([0.0, 0.0, -9.81])
        n.corps("a", 1.0, J, [0.0, 0.0, 0.0])
        n.corps("b", 1.0, J, [0.1, 0.0, 0.0])
        return n

    def restaure_invalide(champ, valeur):
        n = mk()
        avant = n.etat()
        etat = list(n.etat())
        etat[1][0][0] += 0.2
        etat[champ] = valeur
        try:
            n.pose_etat(*etat)
        finally:
            assert n.etat() == avant, "une restauration refusée a modifié le modèle"

    refus = [
        ("sphère : rayon négatif", lambda: mk().sphere(0, [0, 0, 0], -1), ValueError),
        ("sphère : rayon infini", lambda: mk().sphere(0, [0, 0, 0], float("inf")), ValueError),
        ("sphère : centre NaN", lambda: mk().sphere(0, [nan, 0, 0], 1), ValueError),
        ("pas infini", lambda: mk().simule(0.01, float("inf")), ValueError),
        ("loi constante vide", lambda: mk().liaison("x", None, 0,
         cible_t=([1, 0, 0], ("constante", []))), ValueError),
        ("loi linéaire tronquée", lambda: mk().liaison("x", None, 0,
         cible_t=([1, 0, 0], ("lineaire", [0]))), ValueError),
        ("table impaire", lambda: mk().liaison("x", None, 0,
         cible_t=([1, 0, 0], ("table", [0, 0, 1]))), ValueError),
        ("loi non finie", lambda: mk().liaison("x", None, 0,
         cible_t=([1, 0, 0], ("constante", [nan]))), ValueError),
        ("axe de translation hors bornes", lambda: mk().liaison("x", None, 0, bloque_t=[3]), IndexError),
        ("axe de rotation hors bornes", lambda: mk().liaison("x", None, 0, bloque_r=[3]), IndexError),
        ("orientation réfléchie", lambda: mk().corps("miroir", 1, J, [0, 0, 0],
         rot=[-1, 0, 0, 0, 1, 0, 0, 0, 1]), ValueError),
        ("restauration : temps NaN", lambda: restaure_invalide(0, nan), ValueError),
        ("restauration : position NaN", lambda: restaure_invalide(1, [[0, 0, 0], [nan, 0, 0]]), ValueError),
        ("restauration : vitesse infinie", lambda: restaure_invalide(3, [[0, 0, 0], [float("inf"), 0, 0]]), ValueError),
        ("restauration : vitesse angulaire NaN", lambda: restaure_invalide(4, [[0, 0, 0], [nan, 0, 0]]), ValueError),
        ("restauration : orientation dégénérée au dernier corps",
         lambda: restaure_invalide(2, [np.eye(3).ravel().tolist(), [0] * 9]), ValueError),
        ("restauration : orientation réfléchie", lambda: restaure_invalide(2,
         [np.eye(3).ravel().tolist(), [-1, 0, 0, 0, 1, 0, 0, 0, 1]]), ValueError),
        ("restauration : inflow NaN", lambda: restaure_invalide(5, [nan]), ValueError),
        ("pas h = 0 (bouclait à l'infini)", lambda: mk().simule(0.01, 0.0), ValueError),
        ("pas h < 0", lambda: mk().simule(0.01, -1e-4), ValueError),
        ("ρ∞ = 1 (instable en index 3)", lambda: mk().simule(0.01, 1e-4, rho=1.0), ValueError),
        ("liaison d'un corps à lui-même", lambda: mk().liaison("x", 0, 0, bloque_t=[0, 1, 2], bloque_r=[]), ValueError),
        ("liaison bâti–bâti", lambda: mk().liaison("x", None, None, bloque_t=[0, 1, 2], bloque_r=[]), ValueError),
        ("corps inexistant (paniquait)", lambda: mk().liaison("x", None, 7, bloque_t=[0, 1, 2], bloque_r=[]), IndexError),
        ("poutre a == b", lambda: mk().poutre("e", 0, 0, 1e5, 1e5, 1e2, 1e2), ValueError),
        ("poutre sur corps inexistant", lambda: mk().poutre("e", 0, 9, 1e5, 1e5, 1e2, 1e2), IndexError),
        ("distance a == b", lambda: mk().distance("d", 0, 0, [0, 0, 0], [0, 0, 0]), ValueError),
        ("couple bâti–bâti", lambda: mk().couple("c", None, None, [0, 0, 1], ("constant", [1.0])), ValueError),
        ("engrenage a == b", lambda: mk().engrenage("g", 0, 0, [0, 0, 1], [0, 0, 1], 2.0), ValueError),
        ("effort sur corps inexistant", lambda: mk().effort(5, [0, 0, 1], [0, 0, 0]), IndexError),
        ("masse nulle", lambda: Noyau([0, 0, 0]).corps("a", 0.0, J, [0, 0, 0]), ValueError),
        ("inertie non définie positive", lambda: Noyau([0, 0, 0]).corps(
            "a", 1.0, [-1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4], [0, 0, 0]), ValueError),
        # SONDAGE DU 7 SEPT. — trois BLOCAGES (la SVD de nalgebra itère sans
        # plafond et ne converge jamais sur un NaN : `acc_init` figé, aucune
        # trace) et treize acceptations silencieuses, dont un maillage porté par
        # un autre corps que b — la bille tombait À TRAVERS le sol.
        ("masse NaN (BLOQUAIT dans la SVD)", lambda: Noyau([0, 0, 0]).corps("a", nan, J, [0, 0, 0]), ValueError),
        ("gravité NaN (bloquait)", lambda: Noyau([0, 0, nan]), ValueError),
        ("contact : normale nulle (bloquait — normalize() → NaN)",
         lambda: mk().contact("c", 0, [0, 0, 0], 0.05, normale=[0, 0, 0]), ValueError),
        ("contact : maillage porté par un autre corps que b (traversait)",
         lambda: (lambda n: (n.corps("c", 1.0, J, [0.2, 0, 0]), n.maillage("m", 1, *PLAQUE),
                            n.contact("c", 0, [0, 0, 0], 0.05, b=2, maille=0)))(mk()), ValueError),
        ("contact : a == b", lambda: mk().contact("c", 0, [0, 0, 0], 0.05, b=0, rayon_b=0.05), ValueError),
        ("contact : rayon négatif", lambda: mk().contact("c", 0, [0, 0, 0], -0.05), ValueError),
        ("contact : k = 0 (fantôme)", lambda: mk().contact("c", 0, [0, 0, 0], 0.05, k=0.0), ValueError),
        ("contact : barrière IPC à d_hat = 0 (ln 0)",
         lambda: mk().contact("c", 0, [0, 0, 0], 0.05, expo=-1.0, d_hat=0.0), ValueError),
        ("contact : frottement à v_eps = 0 (tanh(v/0))",
         lambda: mk().contact("c", 0, [0, 0, 0], 0.05, mu=0.5, v_eps=0.0), ValueError),
        ("contact : deux primitives sur le second corps",
         lambda: mk().contact("c", 0, [0, 0, 0], 0.05, b=1, demi=[.1, .1, .1], cylindre=([0, 0, 1], .1, .1)), ValueError),
        ("contact : restitution > 1", lambda: mk().contact("c", 0, [0, 0, 0], 0.05, nonlisse=True, restitution=1.5), ValueError),
        ("maillage : triangle dégénéré", lambda: mk().maillage("m", 1, PLAQUE[0], [0, 0, 1]), ValueError),
        ("maillage : sommet NaN", lambda: mk().maillage("m", 1, [nan] * 12, PLAQUE[1]), ValueError),
        ("simule : fin NaN (rendait une liste vide)", lambda: mk().simule(nan, 1e-4), ValueError),
        ("simule : fin avant t courant (liste vide)", lambda: mk().simule(-1.0, 1e-4), ValueError),
        ("simule : tous = 0 (liste vide)", lambda: mk().simule(0.01, 1e-4, tous=0), ValueError),
    ]
    print("  6. robustesse — un modèle absurde doit être refusé, JAMAIS accepté en silence")
    for nom, f, attendu in refus:
        try:
            f()
            raise AssertionError(f"« {nom} » a été ACCEPTÉ — c'est le défaut, pas le contrôle")
        except attendu:
            pass
        except AssertionError:
            raise
        except BaseException as e:                  # un panic PyO3 n'est pas une Exception
            raise AssertionError(f"« {nom} » : attendu {attendu.__name__}, reçu "
                                 f"{type(e).__name__} — {str(e)[:60]}") from None
    print(f"     {len(refus)} modèles absurdes, {len(refus)} refus typés (ValueError / IndexError)")

    # LES DEUX SENS : ce qui est licite doit encore passer, et donner le même
    # résultat qu'avant les gardes. Une garde qui refuse trop est un défaut.
    def poutre_de_longueur_nulle():
        n = Noyau([0, 0, 0])
        n.corps("a", 1.0, J, [0, 0, 0])
        n.corps("b", 1.0, J, [0, 0, 0])
        n.poutre("e", 0, 1, 1e5, 1e5, 1e2, 1e2)
    try:
        poutre_de_longueur_nulle()
        raise AssertionError("poutre de longueur nulle acceptée")
    except ValueError:
        pass
    n = mk()
    assert n.sphere(0, [0, 0, 0], 0.01) == 0
    n.liaison("bon", 0, 1, bloque_t=[0, 1, 2], bloque_r=[])
    n.simule(0.02, 1e-4, tous=10 ** 9)
    z = n.etat()[1][0][2]
    assert z < -1e-5, ("le cas LICITE ne tombe plus : une garde refuse trop", z)
    # et le pas de fin absorbé : t_end/h qui ne tombe pas juste ne doit plus
    # fabriquer un micro-pas ininversible (mesuré : Newton divergeait)
    m = mk()
    m.liaison("bon", 0, 1, bloque_t=[0, 1, 2], bloque_r=[])
    m.simule(0.02, 1.25e-5 * 3, tous=10 ** 9)
    # Un grand pas nominal est réduit à la fin demandée : son accélération
    # initiale doit être celle du pas effectivement intégré, même en reprise.
    for t0 in (0.0, 0.02):
        for h in (0.001, 0.01, 1.0):
            libre = mk()
            if t0:
                libre.simule(t0, 0.001)
            libre.simule(t0 + 0.01, h)
            etat = libre.etat()
            np.testing.assert_allclose(etat[1][0][2], -0.5 * 9.81 * (t0 + 0.01) ** 2, rtol=1e-12)
            np.testing.assert_allclose(etat[3][0][2], -9.81 * (t0 + 0.01), rtol=1e-12)
    # Les perturbations utilisées par Floquet restent ré-orthonormalisables.
    libre = mk()
    etat = list(libre.etat())
    etat[2][1][1] += 1e-6
    libre.pose_etat(*etat)
    rot = np.array(libre.etat()[2][1]).reshape(3, 3)
    np.testing.assert_allclose(rot.T @ rot, np.eye(3), atol=1e-14)
    assert np.linalg.det(rot) > 0
    print(f"     et le licite passe encore : chute libre {1e3 * z:.4f} mm, "
          f"pas de fin absorbé (t/h non entier)")
    return dict(refus=len(refus))


def lecteurs_casses():
    """UN UTILISATEUR TIERS APPORTE AUSSI DES FICHIERS CASSÉS.

    Le patron `robustesse` — un modèle mal posé doit être REFUSÉ AVEC MOTIF,
    jamais accepté en silence — n'avait jamais été appliqué aux LECTEURS. Il
    l'est ici en corrompant des modèles publics du lot, une mutation à la
    fois, et en exigeant un refus motivé.

    Sorti à sa première exécution (3 sept.), trois silences que le noyau
    laissait passer : une masse `nan` (le modèle s'intégrait SANS ERREUR en
    produisant du `nan` partout — le pire des silences), un axe de joint NUL
    (le degré libéré aurait été arbitraire), et un `quat="0 0 0 0"` pris pour
    une rotation. Plus une incohérence de forme : le lecteur MJCF laissait
    fuir un `ParseError` d'implémentation là où l'URDF disait déjà
    « XML invalide ».

    **Un cas reste ACCEPTÉ, et c'est une décision** : l'inertie non physique
    (J1 > J2 + J3). `multibarmech`, cas de référence MBDyn, la viole
    délibérément ; le refuser interdirait un modèle publié. C'est `controle()`
    qui le dit, en lint, sans bloquer la charge.

    Le contrôle négatif est dans le banc : le fichier NON muté doit charger.
    Sans lui, une garde trop large rendrait tout ce vert sans rien prouver —
    « un garde-fou qui crie au loup finit ignoré ».
    """
    import tempfile
    from pathlib import Path

    from vinkulum import mjcf, urdf

    print("  lecteurs devant des fichiers CASSÉS (patron robustesse)")
    d = Path(__file__).parent / "donnees"
    src_u = (d / "urdfs" / "panda.urdf").read_bytes()
    src_m = (d / "menagerie" / "panda.xml").read_bytes()

    def essai(nom, octets, charge, suffixe):
        with tempfile.NamedTemporaryFile("wb", suffix=suffixe, delete=False) as f:
            f.write(octets)
            p = f.name
        try:
            _, i = charge(p)
            return nom, None, i["n_corps"]
        except ValueError as ex:
            return nom, str(ex), None
        finally:
            os.unlink(p)

    # une mutation à la fois, sur un modèle qui charge sans elle
    cas = [
        ("URDF tronqué", src_u[: len(src_u) // 2], urdf.charge, ".urdf"),
        ("URDF masse NaN",
         re.sub(rb'<mass value="[^"]*"', b'<mass value="nan"', src_u, count=1),
         urdf.charge, ".urdf"),
        ("URDF axe nul",
         re.sub(rb'<axis xyz="[^"]*"', b'<axis xyz="0 0 0"', src_u, count=1),
         urdf.charge, ".urdf"),
        ("URDF parent inexistant",
         src_u.replace(b'<parent link="panda_link0"',
                       b'<parent link="pas_la"', 1), urdf.charge, ".urdf"),
        ("MJCF tronqué", src_m[: len(src_m) // 2], mjcf.charge, ".xml"),
        ("MJCF masse négative",
         re.sub(rb'mass="([0-9.]+)"', rb'mass="-\1"', src_m, count=1),
         mjcf.charge, ".xml"),
        ("MJCF inertie NaN",
         re.sub(rb'diaginertia="[^"]*"', b'diaginertia="nan nan nan"',
                src_m, count=1), mjcf.charge, ".xml"),
        ("MJCF quaternion nul",
         re.sub(rb'quat="[^"]*"', b'quat="0 0 0 0"', src_m, count=1),
         mjcf.charge, ".xml"),
    ]
    accepte = []
    for nom, octets, charge, sfx in cas:
        # la mutation doit avoir mordu, sinon le cas ne teste rien
        assert octets != (src_u if sfx == ".urdf" else src_m), \
            f"{nom} : la mutation n'a rien changé — le cas ne teste rien"
        n2, motif, nb = essai(nom, octets, charge, sfx)
        if motif is None:
            accepte.append((nom, nb))
        else:
            assert len(motif) > 20, (nom, "refus sans motif lisible", motif)
            print(f"    {nom:24s} refusé — {motif.split(' : ')[-1][:56]}")
    assert not accepte, ("acceptés EN SILENCE — le pire des défauts", accepte)

    # contrôle négatif : sans mutation, les deux fichiers chargent
    n4, iu = urdf.charge(str(d / "urdfs" / "panda.urdf"))
    _, im = mjcf.charge(str(d / "menagerie" / "panda.xml"))
    # et les DEUX lecteurs publient la même clé — l'incohérence
    # n_links/n_corps est sortie par ce banc même
    assert iu["n_corps"] == 13 and im["n_corps"] == 11, (iu, im)

    # l'acceptation DÉLIBÉRÉE : inertie non physique, dite par controle()
    cass = re.sub(rb'ixx="([0-9.eE+-]+)"', rb'ixx="9\1"', src_u, count=1)
    _, motif, nb = essai("inertie non physique", cass, urdf.charge, ".urdf")
    assert motif is None, ("refusée alors que multibarmech la viole exprès", motif)
    with tempfile.NamedTemporaryFile("wb", suffix=".urdf", delete=False) as f:
        f.write(cass)
        p = f.name
    try:
        n5, _ = urdf.charge(p)
        av = [x for x in n4.controle() if "inégalité" in x or "triangle" in x]
        ap = [x for x in n5.controle() if "inégalité" in x or "triangle" in x]
    finally:
        os.unlink(p)
    assert not av and ap, ("controle() ne voit pas l'inertie non physique", av, ap)
    print(f"    inertie non physique     ACCEPTÉE (multibarmech la viole exprès) "
          f"— dite par controle() : {ap[0][:44]}")
    print(f"    contrôle négatif : les deux fichiers NON mutés chargent "
          f"({iu['n_corps']} et {im['n_corps']} corps), même clé dans les deux")
    return dict(refus=len(cas))

if __name__ == "__main__":
    main()
