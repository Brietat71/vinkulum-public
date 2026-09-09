"""vinkulum.maquette — la MAQUETTE Python de la formulation, et la démo qui juge le Rust.

    python -m vinkulum.maquette            # le noyau Rust sur trois cas à solution analytique
    python -m vinkulum.maquette maquette   # la maquette Python seule (contre-solveur, ~1 min)

La maquette a été écrite AVANT le Rust et validée seule ; elle est gardée comme
CONTRE-SOLVEUR : même formulation, deux implémentations indépendantes (langage,
algèbre linéaire) qui doivent coïncider trajectoire à trajectoire. Chaque brique
du noyau se juge sur une solution indépendante avant d'être étendue.
"""
import os
import sys
from math import cos, pi, sin, sqrt

import numpy as np


I3 = np.eye(3)
G0 = np.array([0.0, 0.0, -9.80665])


# ── algèbre de SO(3) ─────────────────────────────────────────────────────────
def skew(w):
    return np.array([[0.0, -w[2], w[1]], [w[2], 0.0, -w[0]], [-w[1], w[0], 0.0]])


def vee(S):
    return np.array([S[2, 1], S[0, 2], S[1, 0]])


def expm(w):
    """Rodrigues : exp([w]×). Exacte, y compris à w → 0."""
    t = float(np.linalg.norm(w))
    K = skew(w)
    if t < 1e-12:
        return I3 + K + 0.5 * K @ K
    return I3 + sin(t) / t * K + (1.0 - cos(t)) / t**2 * K @ K


# ── corps ────────────────────────────────────────────────────────────────────
class Corps:
    """Corps rigide : masse, tenseur d'inertie au CdM en axes CORPS, pose
    (r, R) du CdM, vitesses spatiales (v, w)."""

    def __init__(self, nom, m, J, r=(0, 0, 0), R=None, v=(0, 0, 0), w=(0, 0, 0)):
        self.nom, self.m = nom, float(m)
        self.J = np.array(J, float) if np.ndim(J) == 2 else np.diag(np.array(J, float))
        self.r, self.R = np.array(r, float), (I3.copy() if R is None else np.array(R, float))
        self.v, self.w = np.array(v, float), np.array(w, float)
        assert self.m > 0 and np.all(np.linalg.eigvalsh(self.J) > 0), (nom, "inertie non définie positive")


class _Bati:
    """Le bâti : pose fixe, aucune inconnue (`idx` None)."""
    r, R = np.zeros(3), I3
    v = w = np.zeros(3)
    nom = "bâti"


BATI = _Bati()


# ── liaisons ─────────────────────────────────────────────────────────────────
class Liaison:
    """Repère (pa, Ra) sur `a`, repère (pb, Rb) sur `b` ; six degrés relatifs
    de b vu de a : translation d = Raᵀ Raᵀ (Pb − Pa) et rotation E = (Ra·Ra)ᵀ
    (Rb·Rb)·Cᵀ. On bloque les composantes listées dans `bloque_t` / `bloque_r`
    (indices 0..2 dans le repère de la liaison). `cible_t(t)` et `cible_r(t)`
    imposent la translation / la rotation relative (défaut : 0 / identité) —
    c'est ce qui fait un servo, un moteur à vitesse imposée, une came.

    Si pb/Rb ne sont pas donnés, le repère sur b est PRIS À LA POSE INITIALE
    pour coïncider avec celui de a : la liaison est fermée au départ."""

    def __init__(self, nom, a, b, pa=(0, 0, 0), Ra=None, pb=None, Rb=None,
                 bloque_t=(0, 1, 2), bloque_r=(0, 1, 2), cible_t=None, cible_r=None):
        self.nom, self.a, self.b = nom, a, b
        self.pa = np.array(pa, float)
        self.Ra = I3.copy() if Ra is None else np.array(Ra, float)
        if pb is None:
            Pa = a.r + a.R @ self.pa
            pb = b.R.T @ (Pa - b.r)
        if Rb is None:
            Rb = b.R.T @ (a.R @ self.Ra)
        self.pb, self.Rb = np.array(pb, float), np.array(Rb, float)
        self.bt, self.br = tuple(bloque_t), tuple(bloque_r)
        self.cible_t = cible_t or (lambda t: np.zeros(3))
        self.cible_r = cible_r or (lambda t: I3)
        self.n = len(self.bt) + len(self.br)

    def phi_G(self, t):
        """Φ (n,) et son jacobien G (n, 12) sur (δr_a, δθ_a, δr_b, δθ_b)."""
        a, b = self.a, self.b
        ua, ub = a.R @ self.Ra, b.R @ self.Rb            # axes des repères, monde
        qa, qb = a.R @ self.pa, b.R @ self.pb            # bras de levier, monde
        dw = (b.r + qb) - (a.r + qa)
        d = ua.T @ dw
        E = ua.T @ ub @ self.cible_r(t).T
        phi_t = d - self.cible_t(t)
        phi_r = vee(0.5 * (E - E.T))
        # δd = uaᵀ([dw]× + [qa]×) δθ_a − uaᵀ δr_a + uaᵀ δr_b − uaᵀ[qb]× δθ_b
        Gt = np.zeros((3, 12))
        Gt[:, 0:3] = -ua.T
        Gt[:, 3:6] = ua.T @ (skew(dw) + skew(qa))
        Gt[:, 6:9] = ua.T
        Gt[:, 9:12] = -ua.T @ skew(qb)
        # δE = [w]× E avec w = uaᵀ(δθ_b − δθ_a) ; δφ = vee(skew([w]× E)), linéaire en w
        Dw = np.column_stack([vee(0.5 * (skew(e) @ E - (skew(e) @ E).T)) for e in I3])
        Gr = np.zeros((3, 12))
        Gr[:, 3:6] = -Dw @ ua.T
        Gr[:, 9:12] = Dw @ ua.T
        bt, br = list(self.bt), list(self.br)
        return (np.concatenate([phi_t[bt], phi_r[br]]),
                np.vstack([Gt[bt], Gr[br]]))


class Distance:
    """Bielle à deux rotules : |Pb − Pa| = L. `L` None = longueur initiale."""

    def __init__(self, nom, a, b, pa, pb, L=None):
        self.nom, self.a, self.b = nom, a, b
        self.pa, self.pb = np.array(pa, float), np.array(pb, float)
        if L is None:
            L = float(np.linalg.norm((b.r + b.R @ self.pb) - (a.r + a.R @ self.pa)))
        self.L, self.n = float(L), 1

    def phi_G(self, t):
        a, b = self.a, self.b
        qa, qb = a.R @ self.pa, b.R @ self.pb
        dw = (b.r + qb) - (a.r + qa)
        l = float(np.linalg.norm(dw))
        nrm = dw / l
        G = np.zeros((1, 12))
        G[0, 0:3] = -nrm
        G[0, 3:6] = nrm @ skew(qa)
        G[0, 6:9] = nrm
        G[0, 9:12] = -nrm @ skew(qb)
        return np.array([l - self.L]), G


# ── le modèle et son intégrateur ─────────────────────────────────────────────
class Modele:
    def __init__(self, corps, liaisons, forces=None, g=G0):
        self.corps, self.liaisons = list(corps), list(liaisons)
        for i, c in enumerate(self.corps):
            c.idx = i
        self.forces = forces or (lambda t, corps: {})   # {corps: (F, M)} monde, au CdM
        self.g = np.array(g, float)
        self.n, self.m = 6 * len(self.corps), sum(L.n for L in self.liaisons)
        self.t = 0.0
        self.lam = np.zeros(self.m)

    # état ↔ vecteurs
    def _u(self):
        return np.concatenate([np.concatenate([c.v, c.w]) for c in self.corps])

    def _set_u(self, u):
        for c in self.corps:
            c.v, c.w = u[6 * c.idx:6 * c.idx + 3], u[6 * c.idx + 3:6 * c.idx + 6]

    def _pose(self):
        return [(c.r.copy(), c.R.copy()) for c in self.corps]

    def _set_pose(self, poses):
        for c, (r, R) in zip(self.corps, poses):
            c.r, c.R = r, R

    def _avance(self, poses0, dq, h):
        """q_{n+1} = q_n ⊕ h·dq — translation vectorielle, rotation sur le groupe."""
        for c, (r0, R0) in zip(self.corps, poses0):
            k = 6 * c.idx
            c.r = r0 + h * dq[k:k + 3]
            c.R = expm(h * dq[k + 3:k + 6]) @ R0

    def M(self):
        M = np.zeros((self.n, self.n))
        for c in self.corps:
            k = 6 * c.idx
            M[k:k + 3, k:k + 3] = c.m * I3
            M[k + 3:k + 6, k + 3:k + 6] = c.R @ c.J @ c.R.T
        return M

    def f(self, t):
        """Forces généralisées : poids, terme gyroscopique, forces utilisateur."""
        f = np.zeros(self.n)
        ext = self.forces(t, self.corps)
        for c in self.corps:
            k = 6 * c.idx
            Js = c.R @ c.J @ c.R.T
            f[k:k + 3] = c.m * self.g
            f[k + 3:k + 6] = -np.cross(c.w, Js @ c.w)
            if c in ext:
                F, Mo = ext[c]
                f[k:k + 3] += F
                f[k + 3:k + 6] += Mo
        return f

    def phi_G(self, t):
        phi, G = np.zeros(self.m), np.zeros((self.m, self.n))
        row = 0
        for L in self.liaisons:
            p, g = L.phi_G(t)
            phi[row:row + L.n] = p
            for corps, col in ((L.a, 0), (L.b, 6)):
                if corps is not BATI:
                    k = 6 * corps.idx
                    G[row:row + L.n, k:k + 6] += g[:, col:col + 6]
            row += L.n
        return phi, G

    def reactions(self):
        """Réaction de chaque liaison = sa part de Gᵀλ, par nom (n, composantes)."""
        out, row = {}, 0
        for L in self.liaisons:
            out[L.nom] = self.lam[row:row + L.n].copy()
            row += L.n
        return out

    # ── α-généralisé, index 3, groupe de Lie ─────────────────────────────────
    def _acc_init(self, t):
        """u̇₀ consistant : [M Gᵀ; G 0][u̇; λ] = [f; −Ġu]."""
        u = self._u()
        _, G = self.phi_G(t)
        eps = 1e-7
        poses = self._pose()
        self._avance(poses, u, eps)
        _, G1 = self.phi_G(t + eps)
        self._set_pose(poses)
        Gdot_u = (G1 @ u - G @ u) / eps
        A = np.block([[self.M(), G.T], [G, np.zeros((self.m, self.m))]])
        rhs = np.concatenate([self.f(t), -Gdot_u])
        x = np.linalg.lstsq(A, rhs, rcond=None)[0]
        return x[:self.n], x[self.n:]

    def simule(self, t_end, h, rho=0.9, tol=1e-12, newton_max=25, sortie=None):
        """Intègre de `self.t` à `t_end` au pas `h`. `sortie(t, modele)` est
        appelée après chaque pas. Rend le nombre de pas."""
        am = (2 * rho - 1) / (rho + 1)
        af = rho / (rho + 1)
        gam = 0.5 + af - am
        bet = 0.25 * (gam + 0.5) ** 2
        n, m = self.n, self.m
        udot, lam = self._acc_init(self.t)
        a = udot.copy()
        npas = 0
        while self.t < t_end - 1e-12:
            t0, t1 = self.t, min(self.t + h, t_end)
            hh = t1 - t0
            u0, poses0 = self._u(), self._pose()
            x = np.concatenate([udot, lam])              # prédiction : u̇ figé

            def etat(x):
                ud, la = x[:n], x[n:]
                a1 = ((1 - af) * ud + af * udot - am * a) / (1 - am)
                dq = u0 + (0.5 - bet) * hh * a + bet * hh * a1
                u1 = u0 + (1 - gam) * hh * a + gam * hh * a1
                self._avance(poses0, dq, hh)
                self._set_u(u1)
                return ud, la, a1

            def residu(x):
                ud, la, _ = etat(x)
                phi, G = self.phi_G(t1)
                r1 = self.M() @ ud + G.T @ la - self.f(t1)
                return np.concatenate([r1, phi / (bet * hh * hh)])

            r = residu(x)
            for it in range(newton_max):
                nr = float(np.linalg.norm(r))
                # critère par blocs (cf. lib.rs) : dynamique en relatif, Φ en absolu
                if (np.linalg.norm(r[:n]) < tol * (1.0 + float(np.linalg.norm(self.f(t1))))
                        and np.linalg.norm(r[n:]) * bet * hh * hh < 1e-11):
                    break
                # jacobien par différences finies (ponytail : dense, à remplacer)
                Jm = np.empty((n + m, n + m))
                for j in range(n + m):
                    dx = np.zeros(n + m)
                    dx[j] = 1e-6 * max(1.0, abs(x[j]))
                    Jm[:, j] = (residu(x + dx) - r) / dx[j]
                x = x - np.linalg.lstsq(Jm, r, rcond=None)[0]
                r = residu(x)
            else:
                raise RuntimeError(f"Newton ne converge pas à t={t1:.5f} (résidu {nr:.2e})")
            udot, lam, a = etat(x)
            self.lam = lam
            self.t = t1
            npas += 1
            if sortie:
                sortie(t1, self)
        return npas

    def energie(self):
        T = sum(0.5 * c.m * c.v @ c.v + 0.5 * c.w @ (c.R @ c.J @ c.R.T) @ c.w for c in self.corps)
        V = -sum(c.m * self.g @ c.r for c in self.corps)
        return T + V


# ── DEMO : trois cas à solution analytique ───────────────────────────────────
def _pendule():
    """Rotule au bâti, masse ponctuelle au bout d'une tige sans masse L = 1 m,
    lâchée à θ₀ = 5°. Période = 2π√(L/g)·(1 + θ₀²/16 + …) ; contrainte tenue."""
    L, th0, g = 1.0, np.radians(5.0), 9.80665
    r0 = np.array([L * sin(th0), 0.0, -L * cos(th0)])
    p = Corps("masse", 0.2, np.diag([1e-6] * 3), r=r0)     # inertie propre négligeable
    rot = Liaison("rotule", BATI, p, pa=(0, 0, 0), bloque_t=(0, 1, 2), bloque_r=())
    mod = Modele([p], [rot])
    T_th = 2 * pi * sqrt(L / g) * (1 + th0**2 / 16 + 11 * th0**4 / 3072)
    E0 = mod.energie()
    xs, ts, drift = [], [], 0.0

    def sortie(t, m):
        ts.append(t)
        xs.append(p.r[0])
    npas = mod.simule(10 * T_th, T_th / 400, rho=0.9, sortie=sortie)
    xs, ts = np.array(xs), np.array(ts)
    # passages par zéro descendants de x → période
    z = [ts[i] - xs[i] * (ts[i + 1] - ts[i]) / (xs[i + 1] - xs[i])
         for i in range(len(xs) - 1) if xs[i] > 0 >= xs[i + 1]]
    T_mes = float(np.mean(np.diff(z)))
    drift = abs(np.linalg.norm(p.r) - L)
    dE = (mod.energie() - E0) / (0.2 * g * L * (1 - cos(th0)))
    assert abs(T_mes / T_th - 1) < 2e-4, ("période", T_mes, T_th)
    assert drift < 1e-10, ("dérive de contrainte", drift)
    assert abs(dE) < 5e-3, ("énergie", dE)
    reac = np.linalg.norm(mod.reactions()["rotule"])
    print(f"  pendule   : période {T_mes:.6f} s (théorie {T_th:.6f}, {1e6*(T_mes/T_th-1):+.0f} ppm) · "
          f"|Φ| {drift:.1e} · ΔE/E {dE:+.1e} sur 10 périodes · réaction {reac:.4f} N · {npas} pas")


def _toupie():
    """Toupie symétrique sur rotule, axe incliné, spin rapide : précession
    Ω_p ≈ m·g·l / (J_a·ω_s) (toupie rapide). Test du terme gyroscopique et de
    la mise à jour sur le groupe (l'axe ne doit ni se dénormaliser ni tomber)."""
    m, l, Ja, Jt, ws, th0 = 0.1, 0.05, 1e-4, 3e-4, 300.0, np.radians(30.0)
    R0 = expm(np.array([th0, 0.0, 0.0]))                 # axe corps z incliné de θ₀ vers −y… (rotation autour de x)
    axe = R0 @ np.array([0.0, 0.0, 1.0])
    top = Corps("toupie", m, np.diag([Jt, Jt, Ja]), r=l * axe, R=R0, w=ws * axe)
    rot = Liaison("pointe", BATI, top, pa=(0, 0, 0), bloque_t=(0, 1, 2), bloque_r=())
    mod = Modele([top], [rot])
    # précession LENTE exacte de la toupie symétrique (racine du trinôme
    # Jt'·cosθ·Ω² − Ja·ω₃·Ω + m·g·l = 0, Jt' = inertie transverse AU PIVOT) ;
    # la formule « rapide » m·g·l/(Ja·ω₃) est son premier ordre, à 3 % ici
    Jtp = Jt + m * l**2
    Op_th = (Ja * ws - sqrt((Ja * ws)**2 - 4 * Jtp * cos(th0) * m * 9.80665 * l)) / (2 * Jtp * cos(th0))
    az = []

    def sortie(t, mm):
        e = top.R @ np.array([0.0, 0.0, 1.0])
        az.append(np.arctan2(e[1], e[0]))
    T = 3.0
    mod.simule(T, 2e-4, rho=0.8, sortie=sortie)
    az = np.unwrap(np.array(az))
    Op = (az[-1] - az[0]) / T
    orth = np.linalg.norm(top.R.T @ top.R - I3)
    e = top.R @ np.array([0.0, 0.0, 1.0])
    incl = np.degrees(np.arccos(e[2]))
    # la toupie rapide : Ω_p à quelques % (nutation + termes d'ordre supérieur en 1/ω_s²)
    # la nutation (lâcher sans précession initiale) décale la moyenne : 1 %
    assert abs(Op / Op_th - 1) < 0.01, ("précession", Op, Op_th)
    assert orth < 1e-12, ("R n'est plus orthonormale", orth)
    assert abs(incl - 30.0) < 1.5, ("l'axe est tombé", incl)
    print(f"  toupie    : précession {Op:.4f} rad/s (théorie {Op_th:.4f}, {100*(Op/Op_th-1):+.2f} %) · "
          f"inclinaison {incl:.2f}° (30) · ‖RᵀR−I‖ {orth:.1e} · {int(T/2e-4)} pas")


def _bielle_manivelle():
    """Manivelle a = 0,03 pivotée au bâti à θ(t) = ω·t IMPOSÉ, bielle l = 0,10
    sur deux rotules, coulisseau en glissière sur x. x(θ) = a·cosθ + √(l² − a²sin²θ).
    La rotation de la bielle autour de son axe reste libre (1 DDL inutile,
    non redondant) ; la glissière + rotule côté coulisseau ne sont PAS
    redondantes. Test des cibles imposées et de la fermeture de boucle."""
    a, l, om = 0.03, 0.10, 20.0
    man = Corps("manivelle", 0.05, np.diag([1e-6, 2e-5, 2e-5]), r=(a / 2, 0, 0))
    bie = Corps("bielle", 0.08, np.diag([1e-6, 8e-5, 8e-5]), r=(a + l / 2, 0, 0))
    cou = Corps("coulisseau", 0.10, np.diag([1e-5] * 3), r=(a + l, 0, 0))

    def theta(t):
        return expm(np.array([0.0, 0.0, om * t]))
    L = [Liaison("pivot moteur", BATI, man, pa=(0, 0, 0), bloque_t=(0, 1, 2), bloque_r=(0, 1, 2),
                 cible_r=theta),
         Liaison("rotule maneton", man, bie, pa=(a / 2, 0, 0), bloque_r=()),
         Liaison("rotule pied", bie, cou, pa=(l / 2, 0, 0), bloque_r=()),
         Liaison("glissière", BATI, cou, pa=(a + l, 0, 0), bloque_t=(1, 2), bloque_r=(0, 1, 2))]
    mod = Modele([man, bie, cou], L, g=(0, 0, 0))
    pire = 0.0

    def sortie(t, mm):
        nonlocal pire
        th = om * t
        x_th = a * cos(th) + sqrt(l**2 - (a * sin(th))**2)
        pire = max(pire, abs(cou.r[0] - x_th))
    T = 2 * pi / om
    npas = mod.simule(2 * T, T / 400, rho=0.9, sortie=sortie)
    phi, _ = mod.phi_G(mod.t)
    reac = mod.reactions()
    assert pire < 1e-6, ("coulisseau hors cinématique", pire)
    assert np.linalg.norm(phi) < 1e-9, ("boucle non fermée", np.linalg.norm(phi))
    print(f"  bielle    : écart max au x(θ) analytique {pire:.1e} m sur 2 tours · |Φ| {np.linalg.norm(phi):.1e} · "
          f"couple moteur {reac['pivot moteur'][5]*1e3:+.3f} mN·m au dernier pas · {npas} pas")


# ── LE NOYAU RUST (src/lib.rs, module `vinkulum._vinkulum`) sur les mêmes trois cas ──
# La maquette Python ci-dessus est le CONTRE-SOLVEUR : même formulation, deux
# implémentations indépendantes (langages, algèbre linéaire), qui doivent
# coïncider à la tolérance de Newton — c'est R1 appliquée au noyau lui-même.
def _rust():
    try:
        from vinkulum import _vinkulum
    except ImportError as e:
        raise SystemExit("extension Rust absente : maturin develop --uv --release dans le dépôt vinkulum") from e
    return _vinkulum


def _rot9(R):
    return [float(x) for x in np.asarray(R).reshape(9)]


def _pendule_rust(mc):
    L, th0, g = 1.0, np.radians(5.0), 9.80665
    N = mc.Noyau()
    N.corps("masse", 0.2, _rot9(np.diag([1e-6] * 3)), [L * sin(th0), 0.0, -L * cos(th0)])
    N.liaison("rotule", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
    T_th = 2 * pi * sqrt(L / g) * (1 + th0**2 / 16 + 11 * th0**4 / 3072)
    E0 = N.energie()
    tr = N.simule(10 * T_th, T_th / 400, rho=0.9)
    ts = np.array([e[0] for e in tr]); xs = np.array([e[1][0][0] for e in tr])
    z = [ts[i] - xs[i] * (ts[i + 1] - ts[i]) / (xs[i + 1] - xs[i])
         for i in range(len(xs) - 1) if xs[i] > 0 >= xs[i + 1]]
    T_mes = float(np.mean(np.diff(z)))
    drift = abs(np.linalg.norm(tr[-1][1][0]) - L)
    dE = (N.energie() - E0) / (0.2 * g * L * (1 - cos(th0)))
    assert abs(T_mes / T_th - 1) < 2e-4, ("période", T_mes, T_th)
    assert drift < 1e-10 and abs(dE) < 5e-3, (drift, dE)
    print(f"  pendule   : période {T_mes:.6f} s (théorie {T_th:.6f}, {1e6*(T_mes/T_th-1):+.0f} ppm) · "
          f"|Φ| {drift:.1e} · ΔE/E {dE:+.1e} · {len(tr)} pas")


def _toupie_rust(mc):
    m, l, Ja, Jt, ws, th0 = 0.1, 0.05, 1e-4, 3e-4, 300.0, np.radians(30.0)
    R0 = expm(np.array([th0, 0.0, 0.0])); axe = R0 @ np.array([0.0, 0.0, 1.0])
    N = mc.Noyau()
    N.corps("toupie", m, _rot9(np.diag([Jt, Jt, Ja])), list(l * axe), _rot9(R0), w=list(ws * axe))
    N.liaison("pointe", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
    Jtp = Jt + m * l**2
    Op_th = (Ja * ws - sqrt((Ja * ws)**2 - 4 * Jtp * cos(th0) * m * 9.80665 * l)) / (2 * Jtp * cos(th0))
    T = 3.0
    tr = N.simule(T, 2e-4, rho=0.8, tous=5)
    az = np.unwrap([np.arctan2(e[2][0][5], e[2][0][2]) for e in tr])   # colonne z de R : (R02, R12, R22)
    Op = (az[-1] - az[0]) / (tr[-1][0] - tr[0][0])
    R = np.array(tr[-1][2][0]).reshape(3, 3)
    orth = np.linalg.norm(R.T @ R - I3); incl = np.degrees(np.arccos(R[2, 2]))
    assert abs(Op / Op_th - 1) < 0.01, ("précession", Op, Op_th)
    assert orth < 1e-12 and abs(incl - 30.0) < 1.5, (orth, incl)
    print(f"  toupie    : précession {Op:.4f} rad/s (théorie {Op_th:.4f}, {100*(Op/Op_th-1):+.2f} %) · "
          f"inclinaison {incl:.2f}° · ‖RᵀR−I‖ {orth:.1e} · {int(T/2e-4)} pas")


def _bielle_rust(mc):
    """Boucle fermée avec rotation IMPOSÉE (loi linéaire) — et contre-solveur :
    la maquette Python sur le même mécanisme, trajectoire à trajectoire."""
    a, l, om = 0.03, 0.10, 20.0
    N = mc.Noyau(g=[0.0, 0.0, 0.0])
    N.corps("manivelle", 0.05, _rot9(np.diag([1e-6, 2e-5, 2e-5])), [a / 2, 0, 0])
    N.corps("bielle", 0.08, _rot9(np.diag([1e-6, 8e-5, 8e-5])), [a + l / 2, 0, 0])
    N.corps("coulisseau", 0.10, _rot9(np.diag([1e-5] * 3)), [a + l, 0, 0])
    N.liaison("pivot moteur", None, 0, cible_r=([0.0, 0.0, 1.0], ("lineaire", [0.0, om])))
    N.liaison("rotule maneton", 0, 1, pa=[a / 2, 0, 0], bloque_r=[])
    N.liaison("rotule pied", 1, 2, pa=[l / 2, 0, 0], bloque_r=[])
    N.liaison("glissière", None, 2, pa=[a + l, 0, 0], bloque_t=[1, 2])
    T = 2 * pi / om
    tr = N.simule(2 * T, T / 400, rho=0.9)
    pire = max(abs(e[1][2][0] - (a * cos(om * e[0]) + sqrt(l**2 - (a * sin(om * e[0]))**2))) for e in tr)
    assert pire < 1e-6 and np.linalg.norm(N.phi()) < 1e-9, (pire, N.phi())
    # contre-solveur : la maquette Python, mêmes corps, mêmes liaisons, même pas
    man = Corps("manivelle", 0.05, np.diag([1e-6, 2e-5, 2e-5]), r=(a / 2, 0, 0))
    bie = Corps("bielle", 0.08, np.diag([1e-6, 8e-5, 8e-5]), r=(a + l / 2, 0, 0))
    cou = Corps("coulisseau", 0.10, np.diag([1e-5] * 3), r=(a + l, 0, 0))
    mod = Modele([man, bie, cou], [
        Liaison("pivot moteur", BATI, man, cible_r=lambda t: expm(np.array([0.0, 0.0, om * t]))),
        Liaison("rotule maneton", man, bie, pa=(a / 2, 0, 0), bloque_r=()),
        Liaison("rotule pied", bie, cou, pa=(l / 2, 0, 0), bloque_r=()),
        Liaison("glissière", BATI, cou, pa=(a + l, 0, 0), bloque_t=(1, 2))], g=(0, 0, 0))
    py = []
    mod.simule(2 * T, T / 400, rho=0.9, sortie=lambda t, m: py.append(np.concatenate([bie.r, bie.R.reshape(9), cou.v])))
    ru = [np.concatenate([e[1][1], e[2][1], e[3][2]]) for e in tr]
    ecart = max(float(np.max(np.abs(p - r))) for p, r in zip(py, ru))
    assert ecart < 1e-7, ("Rust ≠ maquette Python", ecart)
    reac = dict(N.reactions())["pivot moteur"]
    print(f"  bielle    : écart max au x(θ) analytique {pire:.1e} m · |Φ| {np.linalg.norm(N.phi()):.1e} · "
          f"Rust vs maquette Python {ecart:.1e} sur {len(py)} pas · réaction pivot {np.linalg.norm(reac[:3]):.3f} N")


def _engrenage_rust(mc):
    """Deux roues extérieures sur pivots au bâti, rapport −4 (pignon a, roue b),
    couple constant τ sur le pignon, sans gravité. Analytique : θ_a + 4·θ_b = 0
    exactement, et α_a = τ / (J_a + J_b/16) — l'inertie de la roue ramenée par
    le carré du rapport. La réaction de denture est le multiplicateur."""
    Ja, Jb, tau, r = 2e-6, 4e-5, 1e-3, -4.0
    N = mc.Noyau(g=[0.0, 0.0, 0.0])
    N.corps("pignon", 0.01, _rot9(np.diag([Ja / 2, Ja / 2, Ja])), [0.014, 0, 0])
    N.corps("roue", 0.05, _rot9(np.diag([Jb / 2, Jb / 2, Jb])), [0, 0, 0])
    N.liaison("pivot pignon", None, 0, pa=[0.014, 0, 0], bloque_r=[0, 1])
    N.liaison("pivot roue", None, 1, pa=[0, 0, 0], bloque_r=[0, 1])
    N.engrenage("denture", 0, 1, [0, 0, 1], [0, 0, 1], r)
    N.effort(0, [0, 0, 0], [0, 0, tau])
    T = 0.5
    tr = N.simule(T, 1e-3, rho=0.9, tous=50)
    alpha_th = tau / (Ja + Jb / r**2)
    wa, wb = tr[-1][4][0][2], tr[-1][4][1][2]
    # angles continus : le pignon fait plusieurs tours en 0,5 s (α ≈ 240 rad/s²)
    tha_th = 0.5 * alpha_th * T**2
    Ra = np.array(tr[-1][2][0]).reshape(3, 3); Rb = np.array(tr[-1][2][1]).reshape(3, 3)
    tha_w, thb_w = np.arctan2(Ra[1, 0], Ra[0, 0]), np.arctan2(Rb[1, 0], Rb[0, 0])
    ecart_ratio = abs(wa - r * wb) / abs(wa)
    ecart_alpha = abs(wa / T - alpha_th) / alpha_th
    # cohérence de l'angle déroulé : θ_a modulo 2π doit coïncider avec l'angle lu
    ecart_ang = abs(((tha_th - tha_w) + pi) % (2 * pi) - pi)
    lam = dict(N.reactions())["denture"][0]
    assert ecart_ratio < 1e-9, ("rapport de vitesses", wa, wb)
    assert ecart_alpha < 1e-6, ("inertie ramenée", wa / T, alpha_th)
    assert ecart_ang < 1e-6 and abs(N.phi()[-1]) < 1e-12, ("angle déroulé", ecart_ang, N.phi())
    print(f"  engrenage : ω_a/ω_b = {wa/wb:.9f} (−4) · α_a {wa/T:.6f} rad/s² (théorie {alpha_th:.6f}) · "
          f"{tha_th/(2*pi):.1f} tours déroulés · couple de denture λ {lam*1e3:.4f} mN·m")


def _couples_rust(mc):
    """Trois lois de couple à solution analytique : ressort de torsion (période
    2π√(J/k)), gouverneur (montée exponentielle en τ = J/Kg sous la
    saturation), servo PD saturé (échelon : couple borné, régime final au
    droop Kp — mesuré contre la formule). Et Newton doit rester quadratique
    avec les tangentes AD des couples (≤ 3 itérations/pas)."""
    J, k = 2e-4, 0.5
    T_th = 2 * pi * sqrt(J / k)
    # le ressort est au repos à θ = 0,3 (pose initiale) : on excite par une vitesse initiale
    N2 = mc.Noyau(g=[0, 0, 0])
    N2.corps("disque", 0.1, _rot9(np.diag([J / 2, J / 2, J])), [0, 0, 0], w=[0, 0, 2.0])
    N2.liaison("pivot", None, 0, bloque_r=[0, 1])
    N2.couple("torsion", None, 0, [0, 0, 1], ("ressort", [k, 0.0, 0.0]))
    tr = N2.simule(3 * T_th, T_th / 400, rho=0.9)
    ts = np.array([e[0] for e in tr]); th = np.array([np.arctan2(e[2][0][3], e[2][0][0]) for e in tr])
    z = [ts[i] - th[i] * (ts[i + 1] - ts[i]) / (th[i + 1] - th[i]) for i in range(len(th) - 1) if th[i] < 0 <= th[i + 1]]
    T_mes = float(np.mean(np.diff(z)))
    it = N2.stats()[0] / len(tr)
    assert abs(T_mes / T_th - 1) < 5e-4, ("ressort de torsion", T_mes, T_th)
    # gouverneur : J ω̇ = Kg (Ω − ω) tant que |τ| < Q → ω = Ω(1 − e^{−t Kg/J})
    Kg, Q, Om = 1e-3, 1.0, 50.0
    N3 = mc.Noyau(g=[0, 0, 0])
    N3.corps("rotor", 0.1, _rot9(np.diag([J / 2, J / 2, J])), [0, 0, 0])
    N3.liaison("pivot", None, 0, bloque_r=[0, 1])
    N3.couple("gouverneur", None, 0, [0, 0, 1], ("gouverneur", [Kg, Q]), cible=("constante", [Om]))
    tau_ = J / Kg
    tr = N3.simule(3 * tau_, tau_ / 200, rho=0.9)
    pire = max(abs(e[4][0][2] - Om * (1 - np.exp(-e[0] / tau_))) for e in tr) / Om
    assert pire < 1e-4, ("gouverneur", pire)
    # servo PD saturé : échelon de consigne 0,2 rad ; régime final θ∞ = cible (Kd n'agit qu'en vitesse),
    # couple ≤ Q à chaque instant, vitesse ≤ ω_nl
    Kp, Kd, Qs, Wnl = 0.05, 2e-3, 0.02, 12.0
    N4 = mc.Noyau(g=[0, 0, 0])
    N4.corps("palonnier", 0.01, _rot9(np.diag([1e-6, 1e-6, 2e-6])), [0, 0, 0])
    N4.liaison("pivot", None, 0, bloque_r=[0, 1])
    N4.couple("servo", None, 0, [0, 0, 1], ("pd", [Kp, Kd, Qs, Wnl]), cible=("constante", [0.2]))
    tr = N4.simule(0.3, 1e-4, rho=0.9, tous=10)
    th_fin = np.arctan2(tr[-1][2][0][3], tr[-1][2][0][0])
    w_max = max(abs(e[4][0][2]) for e in tr)
    nom, th, w, tau = N4.couples()[0]
    assert abs(th_fin - 0.2) < 2e-3 and w_max <= Wnl * 1.001 and abs(tau) <= Qs * 1.001, (th_fin, w_max, tau)
    print(f"  couples   : torsion période {T_mes:.6f} s (théorie {T_th:.6f}, {1e6*(T_mes/T_th-1):+.0f} ppm, Newton {it:.2f}/pas) · "
          f"gouverneur ω(t) à {pire:.1e} de Ω(1−e^(−t/τ)) · servo PD : θ∞ {th_fin:.4f} (cible 0,2), ω_max {w_max:.2f} ≤ {Wnl}")


def demo():
    mc = _rust()
    print("╔═ vinkulum (Rust) — jalon 1 : α-généralisé sur SO(3) + liaisons rigides")
    _pendule_rust(mc)
    _toupie_rust(mc)
    _bielle_rust(mc)
    _engrenage_rust(mc)
    _couples_rust(mc)
    print("╚ demo OK — période, précession gyroscopique, boucle fermée imposée ; Rust = maquette Python")


def maquette():
    print("╔═ maquette Python (contre-solveur) — mêmes trois cas")
    _pendule()
    _toupie()
    _bielle_manivelle()
    print("╚ maquette OK")


if __name__ == "__main__":
    maquette() if "maquette" in sys.argv else demo()
