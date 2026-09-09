"""vinkulum.adjoint_temps — le gradient d'une TRAJECTOIRE, par l'adjoint discret.

    python -m vinkulum.adjoint_temps      # oscillateur : gradient exact contre DF

`vinkulum.sensibilite` fait l'adjoint sur un ÉQUILIBRE. Ici c'est le long du
temps : l'objectif est une fonctionnelle de la trajectoire, et son gradient
par rapport aux paramètres se calcule en remontant le temps UNE fois, quel
que soit le nombre de paramètres.

L'intégrateur définit implicitement, à chaque pas,

    R(x_{n+1} ; z_n, p) = 0      puis      z_{n+1} = Φ(x_{n+1}, z_n)

avec z = (q, u, u̇, a) l'état complet du schéma et x = (u̇, λ) l'inconnue de
Newton. La chaîne se différencie par le théorème des fonctions implicites, et
l'adjoint remonte :

    μ_N = ∂J/∂z_N ,   μ_{n} = A_nᵀ μ_{n+1} ,   dJ/dp = Σ B_nᵀ μ_{n+1}

où A_n = ∂z_{n+1}/∂z_n et B_n = ∂z_{n+1}/∂p suivent du jacobien du pas
et des formules de Newmark. PontNoyau calcule leurs produits transposés
sans former ces matrices ; il réassemble et factorise le jacobien au recul.

CE QUE CE MODULE COUVRE, ET CE QU'IL NE COUVRE PAS. Trois étages :
  · `Lineaire` — M, C, K constants (modèle NumPy séparé, référence analytique) ;
  · `NonLineaire` — K(t), C(t) réassemblés à chaque pas, même formule ; la
    trajectoire est gardée entière, ou en √N par le CHECKPOINTING de Griewank
    (`checkpoints=K`, rejeu identique au bit — le `revolve` binomial optimal
    n'est pas fait) ;
  · `PontNoyau` — l'adjoint sur la trajectoire du NOYAU lui-même (SO(3),
    liaisons, λ), depuis `simule` + `enregistre_schema` + `k_c_m_g_creux`. Sans
    contact, sans liaison non holonome, sans GGL, sans projection d'invariant. Jugé contre une DF sur le
    noyau : 3e-8 (double pendule), 3e-10 (pendule sphérique tournant), et
    compatible avec la DF, à la marge de dispersion près, sur un pendule
    flexible dont le juge lui-même n'est bon qu'à ~1 % (marches de Newton).
Ce qu'on a trouvé en le construisant, et qui était un défaut du NOYAU :
`raideur` et `amortissement` comptaient le gyroscopique DEUX FOIS (docstring
de `PontNoyau`). Corrigé dans le noyau le soir même ; le pont garde le
doublement comme CONTRÔLE NÉGATIF, et `verification.gyroscopique_tangent`
garde le noyau.
"""
import sys

import numpy as np


def coeffs(rho):
    """Les quatre paramètres de l'α-généralisé, comme le noyau les pose."""
    am = (2.0 * rho - 1.0) / (rho + 1.0)
    af = rho / (rho + 1.0)
    gam = 0.5 + af - am
    bet = 0.25 * (gam + 0.5) ** 2
    return am, af, gam, bet


class Lineaire:
    """M ü + C u̇ + K u = f(t, p), intégré par l'α-généralisé du noyau.

    Volontairement réécrit ici, en petit : l'adjoint doit se juger contre une
    trajectoire dont on contrôle CHAQUE terme, et sur laquelle une référence
    analytique existe. Le pont vers le noyau utilise plus loin les mêmes
    champs linéarisés, via `k_c_m_g_creux`.
    """

    def __init__(self, m, c, k, rho=0.9):
        self.m, self.c, self.k = (np.atleast_2d(np.asarray(x, float)) for x in (m, c, k))
        self.n = self.m.shape[0]
        self.am, self.af, self.gam, self.bet = coeffs(rho)

    def _pas(self, h, q, u, ud, a, f1):
        """Un pas : rend (q1, u1, ud1, a1) et les blocs de sensibilité."""
        am, af, gam, bet = self.am, self.af, self.gam, self.bet
        # a1 = ((1-af)ud1 + af·ud - am·a)/(1-am) ; q1 = q + h[u + (½-β)h a + βh a1]
        ca = (1.0 - af) / (1.0 - am)                    # ∂a1/∂ud1
        cq_u, cq_a, cq_ud = h, h * (0.5 - bet) * h - h * bet * h * am / (1 - am), \
            h * bet * h * af / (1 - am)
        cu_u, cu_a, cu_ud = 1.0, (1 - gam) * h - gam * h * am / (1 - am), \
            gam * h * af / (1 - am)
        # R(ud1) = M ud1 + C u1 + K q1 - f1 = 0, avec q1 et u1 affines en ud1
        dq_dud1 = h * bet * h * ca
        du_dud1 = gam * h * ca
        j = self.m + self.c * du_dud1 + self.k * dq_dud1
        rhs = f1 - self.c @ (u * cu_u + a * cu_a + ud * cu_ud) \
            - self.k @ (q + u * cq_u + a * cq_a + ud * cq_ud)
        ud1 = np.linalg.solve(j, rhs)
        a1 = ((1 - af) * ud1 + af * ud - am * a) / (1 - am)
        q1 = q + u * cq_u + a * cq_a + ud * cq_ud + dq_dud1 * ud1
        u1 = u * cu_u + a * cu_a + ud * cu_ud + du_dud1 * ud1
        # ∂ud1/∂(q, u, ud, a) — par le théorème des fonctions implicites
        ji = np.linalg.inv(j)
        d_ud1 = dict(
            q=-ji @ self.k,
            u=-ji @ (self.c * cu_u + self.k * cq_u),
            ud=-ji @ (self.c * cu_ud + self.k * cq_ud),
            a=-ji @ (self.c * cu_a + self.k * cq_a),
            f=ji)
        return (q1, u1, ud1, a1), d_ud1, (cq_u, cq_a, cq_ud, cu_u, cu_a, cu_ud,
                                          dq_dud1, du_dud1, ca)

    def _bloc_a(self, d_ud1, co):
        """A = ∂(q1,u1,ud1,a1)/∂(q,u,ud,a), en blocs (4n × 4n)."""
        cq_u, cq_a, cq_ud, cu_u, cu_a, cu_ud, dq, du, ca = co
        n, i = self.n, np.eye(self.n)
        am, af = self.am, self.af
        a = np.zeros((4 * n, 4 * n))
        for j, cle in enumerate(("q", "u", "ud", "a")):
            d = d_ud1[cle]
            base_q = {"q": i, "u": i * cq_u, "ud": i * cq_ud, "a": i * cq_a}[cle]
            base_u = {"q": np.zeros((n, n)), "u": i * cu_u, "ud": i * cu_ud,
                      "a": i * cu_a}[cle]
            a[0:n, j * n:(j + 1) * n] = base_q + dq * d
            a[n:2 * n, j * n:(j + 1) * n] = base_u + du * d
            a[2 * n:3 * n, j * n:(j + 1) * n] = d
            da = (1 - af) / (1 - am) * d
            if cle == "ud":
                da = da + af / (1 - am) * i
            if cle == "a":
                da = da - am / (1 - am) * i
            a[3 * n:4 * n, j * n:(j + 1) * n] = da
        return a

    def integre(self, h, nsteps, force, q0=None, u0=None):
        """Trajectoire + les blocs A_n et ∂z_{n+1}/∂f_{n+1}, gardés pour l'adjoint."""
        n = self.n
        q = np.zeros(n) if q0 is None else np.array(q0, float)
        u = np.zeros(n) if u0 is None else np.array(u0, float)
        ud = np.linalg.solve(self.m, force(0.0) - self.c @ u - self.k @ q)
        a = ud.copy()
        traj, blocs = [(0.0, q.copy(), u.copy())], []
        for i in range(nsteps):
            t1 = (i + 1) * h
            (q, u, ud, a), d, co = self._pas(h, q, u, ud, a, force(t1))
            aa = self._bloc_a(d, co)
            bb = np.zeros((4 * n, n))
            bb[0:n] = co[6] * d["f"]
            bb[n:2 * n] = co[7] * d["f"]
            bb[2 * n:3 * n] = d["f"]
            bb[3 * n:4 * n] = (1 - self.af) / (1 - self.am) * d["f"]
            blocs.append((aa, bb, t1))
            traj.append((t1, q.copy(), u.copy()))
        return traj, blocs

    def gradient(self, h, nsteps, force, d_force, dj_dz_final):
        """dJ/dp par l'adjoint : UNE remontée, quel que soit le nombre de p.

        `d_force(t)` rend ∂f/∂p à l'instant t, en (n × n_p).
        `dj_dz_final` : ∂J/∂z au dernier instant, en 4n.
        """
        _, blocs = self.integre(h, nsteps, force)
        mu = np.array(dj_dz_final, float)
        n_p = d_force(0.0).shape[1]
        grad = np.zeros(n_p)
        for aa, bb, t1 in reversed(blocs):
            grad += (bb.T @ mu) @ d_force(t1)
            mu = aa.T @ mu
        # ⚠ L'ÉTAT INITIAL DÉPEND DU PARAMÈTRE, et l'oublier ne se voit PAS
        # partout : l'accélération de départ est ü₀ = M⁻¹(f(0) − Cu₀ − Kq₀),
        # donc elle porte ∂f(0)/∂p. Sur un forçage en sin(ωt) le terme est nul
        # à t = 0 et le gradient sort juste ; sur un cos(ωt) il ne l'est pas, et
        # le gradient est faux de 1,3 % — un gradient PARTIELLEMENT juste, ce
        # qui est le pire cas pour être repéré.
        n = self.n
        d_ud0 = np.linalg.solve(self.m, d_force(0.0))
        grad += mu[2 * n:3 * n] @ d_ud0 + mu[3 * n:4 * n] @ d_ud0
        return grad


class NonLineaire(Lineaire):
    """M ü + g(q, u, t, p) = 0 — la RAIDEUR ET L'AMORTISSEMENT DÉPENDENT DE L'ÉTAT.

    C'est ce qui manquait pour optimiser sur une trajectoire complète : une
    machine réelle est non linéaire (grandes rotations, contact, aéro), et un
    adjoint qui suppose K constant ne peut l'accompagner que près d'un point
    de fonctionnement.

    LA FORMULE NE CHANGE PAS — c'est ce que le cas linéaire annonçait, et
    c'est vérifié ici : la récurrence μ_n = A_nᵀ μ_{n+1} tient telle quelle,
    seule A_n est évaluée avec les tangentes AU POINT COURANT. Ce qui change
    est le COÛT et la mémoire : il faut la trajectoire pour remonter, donc on
    la stocke (le *checkpointing* de Griewank la rachèterait en √N, non fait
    et nommé).

    `forces(q, u, t, p)` rend g ; `tangentes(q, u, t, p)` rend (K, C) =
    (∂g/∂q, ∂g/∂u) ; `d_forces(q, u, t, p)` rend ∂g/∂p. La masse reste
    constante — une masse variable est un autre sujet (systèmes ouverts).
    """

    def __init__(self, m, forces, tangentes, d_forces, rho=0.9):
        n = np.atleast_2d(np.asarray(m, float)).shape[0]
        super().__init__(m, np.zeros((n, n)), np.zeros((n, n)), rho=rho)
        self.g, self.dg, self.dgp = forces, tangentes, d_forces

    def _pas_nl(self, h, q, u, ud, a, t1, p, iters=30, tol=1e-13):
        """Un pas de Newton sur R(ud1) = M ud1 + g(q1, u1, t1) = 0.

        Le linéaire résolvait en UNE fois ; ici il faut itérer, et c'est la
        seule différence de fond. Les coefficients du schéma sont les mêmes.
        """
        am, af, gam, bet = self.am, self.af, self.gam, self.bet
        ca = (1.0 - af) / (1.0 - am)
        cq_u = h
        cq_a = h * (0.5 - bet) * h - h * bet * h * am / (1 - am)
        cq_ud = h * bet * h * af / (1 - am)
        cu_u, cu_a = 1.0, (1 - gam) * h - gam * h * am / (1 - am)
        cu_ud = gam * h * af / (1 - am)
        dq_dud1, du_dud1 = h * bet * h * ca, gam * h * ca
        base_q = q + u * cq_u + a * cq_a + ud * cq_ud
        base_u = u * cu_u + a * cu_a + ud * cu_ud
        ud1 = ud.copy()
        for _ in range(iters):
            q1, u1 = base_q + dq_dud1 * ud1, base_u + du_dud1 * ud1
            r = self.m @ ud1 + self.g(q1, u1, t1, p)
            if np.linalg.norm(r) <= tol * (1.0 + np.linalg.norm(self.m @ ud1)):
                break
            kk, cc = self.dg(q1, u1, t1, p)
            j = self.m + cc * du_dud1 + kk * dq_dud1
            ud1 = ud1 - np.linalg.solve(j, r)
        q1, u1 = base_q + dq_dud1 * ud1, base_u + du_dud1 * ud1
        kk, cc = self.dg(q1, u1, t1, p)
        j = self.m + cc * du_dud1 + kk * dq_dud1
        ji = np.linalg.inv(j)
        a1 = ((1 - af) * ud1 + af * ud - am * a) / (1 - am)
        d_ud1 = dict(q=-ji @ kk,
                     u=-ji @ (cc * cu_u + kk * cq_u),
                     ud=-ji @ (cc * cu_ud + kk * cq_ud),
                     a=-ji @ (cc * cu_a + kk * cq_a),
                     f=ji)
        return (q1, u1, ud1, a1), d_ud1, (cq_u, cq_a, cq_ud, cu_u, cu_a, cu_ud,
                                          dq_dud1, du_dud1, ca), (q1, u1)

    def _marche(self, h, i0, n_pas, p, etat):
        """`n_pas` pas depuis l'état complet `etat` = (q, u, u̇, a) au pas i0.

        Rend (blocs, état final). C'est la brique que le checkpointing rejoue
        segment par segment — la MÊME arithmétique que l'aller, donc les blocs
        recalculés sont identiques au bit."""
        n = self.n
        q, u, ud, a = (np.array(x, float) for x in etat)
        blocs = []
        for i in range(i0, i0 + n_pas):
            t1 = (i + 1) * h
            (q, u, ud, a), d, co, pt = self._pas_nl(h, q, u, ud, a, t1, p)
            aa = self._bloc_a(d, co)
            # ∂z₁/∂p passe par ∂g/∂p AU POINT COURANT — c'est là que le non
            # linéaire diffère : la dérivée du forçage se lit sur la
            # trajectoire, pas sur une loi figée.
            dgp = self.dgp(pt[0], pt[1], t1, p)
            bb = np.zeros((4 * n, dgp.shape[1]))
            df = -d["f"] @ dgp
            bb[0:n] = co[6] * df
            bb[n:2 * n] = co[7] * df
            bb[2 * n:3 * n] = df
            bb[3 * n:4 * n] = (1 - self.af) / (1 - self.am) * df
            blocs.append((aa, bb, t1))
        return blocs, (q, u, ud, a)

    def _etat0(self, p, q0=None, u0=None):
        n = self.n
        q = np.zeros(n) if q0 is None else np.array(q0, float)
        u = np.zeros(n) if u0 is None else np.array(u0, float)
        ud = np.linalg.solve(self.m, -self.g(q, u, 0.0, p))
        return q, u, ud, ud.copy()

    def integre(self, h, nsteps, p, q0=None, u0=None):
        etat = self._etat0(p, q0, u0)
        traj = [(0.0, etat[0].copy(), etat[1].copy())]
        blocs = []
        for i in range(nsteps):
            b, etat = self._marche(h, i, 1, p, etat)
            blocs += b
            traj.append(((i + 1) * h, etat[0].copy(), etat[1].copy()))
        return traj, blocs

    def gradient(self, h, nsteps, p, dj_dz_final, checkpoints=None):
        """dJ/dp par l'adjoint, sur une trajectoire NON linéaire.

        `checkpoints=None` : la trajectoire entière est gardée (N blocs A_n,
        B_n). `checkpoints=K` : CHECKPOINTING de Griewank à deux niveaux — en
        aller on ne garde que K états complets (q, u, u̇, a) aux frontières de
        segments ; en remontée chaque segment est REJOUÉ depuis son point de
        contrôle, ses blocs servent puis sont jetés. Mémoire : K états +
        ⌈N/K⌉ blocs au lieu de N blocs — en √N avec K = ⌈√N⌉. Le rejeu est la
        même arithmétique que l'aller, donc le gradient est identique au bit
        (mesuré : écart 0.0). Coût : un aller de plus.

        Ce n'est PAS le `revolve` binomial optimal de Griewank–Walther (1997),
        qui hiérarchise les points de contrôle en log N : √N est le premier
        étage, et il suffit tant que la mémoire n'est pas le mur. À passer
        en binomial le jour où elle le redevient.

        `self.pic_blocs` et `self.pic_etats` : le pic de blocs et d'états
        simultanément gardés au dernier appel — ce que le banc asserte.
        """
        n = self.n
        mu = np.array(dj_dz_final, float)
        if checkpoints is None:
            _, blocs = self.integre(h, nsteps, p)
            self.pic_blocs, self.pic_etats = len(blocs), 1
            grad = np.zeros(blocs[0][1].shape[1])
            for aa, bb, _ in reversed(blocs):
                grad += bb.T @ mu
                mu = aa.T @ mu
        else:
            k = max(1, int(checkpoints))
            lg = -(-nsteps // k)                       # longueur de segment ⌈N/K⌉
            bornes = list(range(0, nsteps, lg)) + [nsteps]
            etat = self._etat0(p)
            pts = []                                   # les K points de contrôle
            for i0, i1 in zip(bornes[:-1], bornes[1:]):
                pts.append((i0, tuple(x.copy() for x in etat)))
                _, etat = self._marche(h, i0, i1 - i0, p, etat)
            self.pic_etats = len(pts)
            self.pic_blocs = 0
            grad = None
            for (i0, e0), i1 in zip(reversed(pts), reversed(bornes[1:])):
                blocs, _ = self._marche(h, i0, i1 - i0, p, e0)   # rejeu du segment
                self.pic_blocs = max(self.pic_blocs, len(blocs))
                if grad is None:
                    grad = np.zeros(blocs[0][1].shape[1])
                for aa, bb, _ in reversed(blocs):
                    grad += bb.T @ mu
                    mu = aa.T @ mu
        # l'état initial dépend du paramètre, ici aussi : ü₀ = −M⁻¹g(q₀,u₀,0)
        q0, u0 = np.zeros(n), np.zeros(n)
        d_ud0 = -np.linalg.solve(self.m, self.dgp(q0, u0, 0.0, p))
        grad += mu[2 * n:3 * n] @ d_ud0 + mu[3 * n:4 * n] @ d_ud0
        return grad


def demo_nonlineaire():
    """L'adjoint NON LINÉAIRE contre une différence finie, sur un DUFFING.

    Le juge est une DF sur la trajectoire complète : elle ne partage avec
    l'adjoint ni la remontée, ni les tangentes, ni la formule. Et l'oscillateur
    de Duffing (raideur cubique) est non linéaire pour de bon — à l'amplitude
    retenue, le terme cubique pèse la moitié du terme linéaire, donc un
    adjoint qui figerait K se tromperait franchement.
    """
    m, k1, k3, c = 1.0, 40.0, 3.0e3, 0.30
    amp, om = 6.0, 5.0
    h, ns = 2e-4, 3000

    def g(q, u, t, p):
        return np.array([k1 * q[0] + k3 * q[0] ** 3 + c * u[0]
                         - (p[0] * amp * np.sin(om * t) + p[1] * amp * np.cos(om * t))])

    def dg(q, u, t, p):
        return np.array([[k1 + 3.0 * k3 * q[0] ** 2]]), np.array([[c]])

    def dgp(q, u, t, p):
        return np.array([[-amp * np.sin(om * t), -amp * np.cos(om * t)]])

    md = NonLineaire([[m]], g, dg, dgp)
    p0 = np.array([1.0, 0.6])

    def j_de(p):
        tr, _ = md.integre(h, ns, p)
        return tr[-1][1][0]                      # J = q(T)

    # ∂J/∂z_N : J ne dépend que de q ⇒ (1, 0, 0, 0)
    grad = md.gradient(h, ns, p0, np.array([1.0, 0.0, 0.0, 0.0]))
    df = np.zeros(2)
    for i in range(2):
        e = np.zeros(2)
        e[i] = 1e-6
        df[i] = (j_de(p0 + e) - j_de(p0 - e)) / 2e-6
    # amplitude du non-linéaire : combien le terme cubique pèse
    tr, _ = md.integre(h, ns, p0)
    qmax = max(abs(x[1][0]) for x in tr)
    part = k3 * qmax ** 2 / (k1 + k3 * qmax ** 2)
    print("╔═ vinkulum — adjoint NON LINÉAIRE en temps (Duffing, raideur cubique)")
    print(f"║ |q| max {qmax:.4f} — le terme cubique pèse {100 * part:.0f} % de la raideur")
    for i in range(2):
        print(f"║   p{i}   adjoint {grad[i]:+.9e}   DF {df[i]:+.9e}")
    ec = max(abs(grad[i] / df[i] - 1) for i in range(2))
    print(f"║ écart à la DF {100 * ec:.2e} %   ·   UNE remontée pour les 2 paramètres, "
          f"la DF en demande 4 trajectoires")
    assert ec < 1e-4, ("l'adjoint non linéaire s'écarte de la DF", grad, df)
    assert part > 0.3, ("le cas n'est pas assez non linéaire pour juger", part)
    # CHECKPOINTING DE GRIEWANK : même gradient, mémoire en √N. Le rejeu est
    # la même arithmétique, donc l'écart attendu est ZÉRO — pas « petit ».
    pic_plein = md.pic_blocs
    kc = int(np.ceil(np.sqrt(ns)))
    g_ck = md.gradient(h, ns, p0, np.array([1.0, 0.0, 0.0, 0.0]), checkpoints=kc)
    e_ck = float(np.max(np.abs(g_ck - grad) / np.maximum(np.abs(grad), 1e-30)))
    print(f"║ Griewank √N : {kc} points de contrôle, pic {md.pic_blocs} blocs + {md.pic_etats} états "
          f"contre {pic_plein} blocs ; écart au plein-stockage {e_ck:.1e}")
    assert e_ck < 1e-12, ("le rejeu depuis un point de contrôle change le gradient", e_ck)
    assert md.pic_blocs <= kc + 1 and md.pic_etats <= kc + 1, (md.pic_blocs, md.pic_etats, kc)
    # CONTRÔLE NÉGATIF du compteur : sans checkpointing il doit valoir N
    assert pic_plein == ns, pic_plein
    print("╚═ adjoint non linéaire OK")
    return dict(grad=grad.tolist(), df=df.tolist(), ecart=float(ec), part=float(part),
                griewank=dict(k=kc, pic_blocs=md.pic_blocs, pic_etats=md.pic_etats, ecart=e_ck))


# ── LE PONT VERS LE NOYAU ─────────────────────────────────────────────────────
def _sk(v):
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def _expm(th):
    k = _sk(th)
    f = np.linalg.norm(th)
    if f < 1e-8:
        return np.eye(3) + k + 0.5 * k @ k
    return np.eye(3) + np.sin(f) / f * k + (1 - np.cos(f)) / f ** 2 * k @ k


def _jl(th):
    """Jacobien GAUCHE de exp sur SO(3) : δ(exp(θ)) = [J_l δθ]ₓ exp(θ)."""
    k = _sk(th)
    f = np.linalg.norm(th)
    if f < 1e-6:
        return np.eye(3) + 0.5 * k + k @ k / 6.0
    return np.eye(3) + (1 - np.cos(f)) / f ** 2 * k + (f - np.sin(f)) / f ** 3 * k @ k


def _sk_blocs(v):
    """Matrices antisymétriques de vecteurs (nombre de corps × 3)."""
    out = np.zeros((len(v), 3, 3))
    out[:, 0, 1], out[:, 0, 2] = -v[:, 2], v[:, 1]
    out[:, 1, 0], out[:, 1, 2] = v[:, 2], -v[:, 0]
    out[:, 2, 0], out[:, 2, 1] = -v[:, 1], v[:, 0]
    return out


class PontNoyau:
    """L'ADJOINT EN TEMPS SUR LA VRAIE MACHINE — feuille de route n°2, premier pas.

    Jusqu'ici `Lineaire` et `NonLineaire` tournaient sur un modèle NumPy
    séparé (M ü + g, q vectoriel, sans Φ, sans λ). Ici la trajectoire est
    celle du NOYAU : `simule(tous=1)` rend (q, u, λ) par pas et
    `enregistre_schema` + `schema()` rendent les deux accélérations (u̇, a)
    que l'α-généralisé traîne d'un pas à l'autre. L'état du schéma est
    z = (q, u, u̇, a) avec q sur SO(3) — les rotations sont vécues en
    perturbations SPATIALES à gauche (δR = [δθ]ₓR), comme dans le noyau.

    Ce pont utilise exclusivement sigma_lie=0 (schéma historique). La
    correction sigma optionnelle de la version 0.8 n'est pas couverte par
    ses dérivées temporelles ; `aller` appelle `simule` avec le défaut zéro.

    Le pas du noyau à sigma=0, tel que `Pas::a1_dq`, `etat`, `avance` l'écrivent :

        a₁ = c u̇₁ + k_a u̇₀ − k_b a₀      c = (1−α_f)/(1−α_m), k_a = α_f/(1−α_m), k_b = α_m/(1−α_m)
        dq = u₀ + (½−β)h a₀ + βh a₁ ;   r₁ = r₀ + h dq_t ;   R₁ = exp(h dq_r) R₀
        u₁ = u₀ + (1−γ)h a₀ + γh a₁
        R(u̇₁, λ₁) = M(q₁)u̇₁ + ω₁×J_s(q₁)ω₁ − f(q₁, u₁) + G(q₁)ᵀλ₁ = 0 ,  Φ(q₁) = 0

    et sa règle de chaîne sur SO(3) : pour R₁ = exp(θ̂)R₀,
    δθ₁ = J_l(θ) δθ_vec + exp(θ̂) δθ₀. A_n = ∂z₁/∂z₀ et B_n = ∂z₁/∂p sortent
    du théorème des fonctions implicites sur (u̇₁, λ₁) — le système linéaire
    est EXACTEMENT le jacobien de Newton du noyau (en-tête de `tangent.rs`).

    LES TANGENTES VIENNENT DU NOYAU (`k_c_m_g_creux` au point posé, λ posé par
    `pose_lam`), AVEC UNE CORRECTION faite ici en Python : il MANQUE à K le
    terme d'inertie en orientation ∂(J_s u̇)/∂θ = −[J_s u̇]ₓ + J_s[u̇]ₓ —
    `raideur` dérive un champ de forces à u̇ = 0.
    ET UN DÉFAUT DU NOYAU TROUVÉ ICI, corrigé le soir même : `raideur` et
    `amortissement` comptaient le gyroscopique ω×J_sω DEUX FOIS (`forces()`
    le porte, et les deux l'ajoutaient encore « à la main » — mesuré sur un
    corps libre, K_θθ et C_ωω à exactement 2,000× l'analytique). Le pont l'a
    d'abord compensé chez lui ; depuis la correction il ne touche plus à K et
    C, et garde le DOUBLEMENT comme contrôle négatif (`gyro=True` remet la
    copie en trop : l'écart à la DF doit exploser).
    Le contrôle négatif du banc altère chacune des trois pièces (J_l, K_M,
    gyroscopique) et montre que l'écart à la DF GRANDIT.

    DOMAINE : pas de contact (le redémarrage de schéma après un pas à contact
    actif — `acc_init` dans `simule` — casse la forme lisse de A_n : refusé
    avec sa raison). Les candidats de contact automatique et les éléments
    aérodynamiques et les liaisons non holonomes sont également refusés.
    Pas de GGL, pas de projection d'invariant, pas de pas
    adaptatif. L'état initial dépend du paramètre par l'accélération
    consistante (`acc_init`, non exposée) : ∂(u̇₀, a₀)/∂p est pris par
    DIFFÉRENCE FINIE sur un pas de `simule` — dit, pas caché.

    Le recul résout Jᵀψ = Φ_xᵀμ, puis applique Φ_zᵀμ − R_zᵀψ et −R_pᵀψ.
    Il ne forme ni Z, ni A (4n × 4n), ni B. K et C sont assemblés localement.
    `creux=None` choisit SciPy à partir de 180 degrés physiques lorsqu'il est
    disponible ; `False` garde NumPy, `True` exige SciPy. Le LU est recalculé
    à chaque pas, avec équilibrage et contrôle de l'erreur arrière.

    `d_residu_p(N,p)` garde son contrat historique (n × n_p), pour des
    paramètres de forces avec masse et contraintes indépendantes de p.
    En alternative, `d_residu_transpose(N,p,ud,psi)` fournit R_pᵀψ pour le
    résidu complet (forces, contraintes) à q,u,ud,λ fixés. psi porte n+m
    composantes, ud les n accélérations physiques. Ce rappel peut éviter
    aussi la matrice n × n_p ; le signe moins est appliqué par le pont.

    `d_initial_transpose(p,mu)` peut fournir (∂z₀/∂p)ᵀmu directement.
    En son absence, le chemin historique suppose q₀ et u₀ indépendants de p
    et évalue les dérivées de l'accélération initiale par différences finies.
    Ces rappels doivent inclure toutes les dépendances de leur résidu ou
    de leur état initial. Une topologie et un rang de contraintes stables
    sont nécessaires ; les historiques internes non enregistrés et les
    changements de branche ne sont pas différenciés par ce pont.

    `n_kcmz` reste le compteur historique des linéarisations, bien qu'aucun
    appel à k_c_m_z ne soit nécessaire au gradient. `residu_adjoint_max`
    mesure l'erreur arrière en norme infinie du système équilibré ; il ne
    certifie ni le conditionnement ni l'erreur du gradient physique.
    """

    def __init__(self, construit, d_residu_p=None, rho=0.9, jl=True, k_m=True, gyro=False,
                 creux=None, d_residu_transpose=None, d_initial_transpose=None):
        self.construit, self.d_residu_p = construit, d_residu_p
        self.rho = rho
        self.am, self.af, self.gam, self.bet = coeffs(rho)
        self.jl, self.k_m, self.gyro = jl, k_m, gyro   # contrôles négatifs
        self.n_kcmz = 0
        if creux is not None and type(creux) is not bool:
            raise ValueError("creux doit être None, True ou False")
        self.creux = creux
        self.d_residu_transpose = d_residu_transpose
        self.d_initial_transpose = d_initial_transpose
        self.residu_adjoint_max = 0.
        self.n_resolutions_adjoint = 0

    def _utilise_creux(self, n):
        if self.creux is False or (self.creux is None and n < 180):
            return False
        try:
            import scipy.sparse  # optionnel ; le chemin NumPy reste disponible
        except ImportError:
            if self.creux:
                raise ImportError("PontNoyau(creux=True) nécessite SciPy") from None
            return False
        return True

    def aller(self, p, t_end, h):
        N = self.construit(p)
        if N.contacts():
            raise ValueError("PontNoyau : contact présent — le redémarrage de schéma après "
                             "impact (acc_init dans simule) casse la forme lisse de A_n ; "
                             "hors domaine, à traiter à part")
        N._verifie_domaine_adjoint()
        t0, r0, rot0, v0, w0, _ = N.etat()
        N.enregistre_schema(True)
        fr = N.simule(t_end, h, tous=1, rho=self.rho)
        sch = N.schema()
        assert len(sch) == len(fr) + 1, (len(sch), len(fr))
        fr = [(t0, r0, rot0, v0, w0, list(sch[0][2]))] + list(fr)
        return N, fr, sch

    def _tangentes(self, N, frame, ud, p, creux=False, parametres=True, historique=False):
        t, r, rot, v, w, lam = frame
        N.pose_etat(t, [list(x) for x in r], [list(x) for x in rot],
                    [list(x) for x in v], [list(x) for x in w])
        N.pose_lam(list(lam))
        self.n_kcmz += 1
        if historique:
            k, c, m, _, g = (np.array(x, float) for x in N.k_c_m_z(t))
            g = g.reshape(-1, k.shape[0]) if k.shape[0] else g.reshape(0, 0)
        elif creux:
            from scipy.sparse import csc_matrix, coo_matrix
            k, c, m, g = (csc_matrix((data, indices, pointers), shape=(nr, nc))
                         for nr, nc, pointers, indices, data in N.k_c_m_g_creux(t))
        else:
            def dense(export):
                nr, nc, pointers, indices, data = export
                a = np.zeros((nr, nc))
                if len(data):
                    cols = np.repeat(np.arange(nc), np.diff(pointers))
                    a[np.asarray(indices), cols] = data
                return a
            k, c, m, g = map(dense, N.k_c_m_g_creux(t))
        nb = len(r)
        indices = 6*np.arange(nb)[:, None] + np.arange(3, 6)
        rows, cols = np.repeat(indices, 3, axis=1).ravel(), np.tile(indices, (1, 3)).ravel()
        js = np.asarray(m[rows, cols]).reshape(nb, 3, 3)
        wi = np.asarray(w, float).reshape(nb, 3)
        udi = np.asarray(ud, float).reshape(nb, 6)[:, 3:]
        dk, dc = np.zeros_like(js), np.zeros_like(js)
        if self.gyro:  # contrôle négatif : double compte historique
            sw = _sk_blocs(wi)
            sjw = _sk_blocs((js @ wi[..., None])[..., 0])
            dk += sw @ (js @ sw - sjw)
            dc += sw @ js - sjw
        if self.k_m:
            dk += -_sk_blocs((js @ udi[..., None])[..., 0]) + js @ _sk_blocs(udi)
        if creux:
            k = k + coo_matrix((dk.ravel(), (rows, cols)), shape=k.shape).tocsc()
            if self.gyro:
                c = c + coo_matrix((dc.ravel(), (rows, cols)), shape=c.shape).tocsc()
        else:
            k[rows, cols] += dk.ravel()
            c[rows, cols] += dc.ravel()
        if not parametres:
            return k, c, m, g, None
        if self.d_residu_p is None:
            raise ValueError("fournir d_residu_p ou d_residu_transpose pour le gradient")
        rp = np.atleast_2d(np.asarray(self.d_residu_p(N, p), float))
        if rp.shape[0] != k.shape[0]:
            rp = rp.T
        if rp.shape != (k.shape[0], len(p)) or not np.isfinite(rp).all():
            raise ValueError("d_residu_p doit fournir un tableau fini de forme (6n, nombre de paramètres)")
        return k, c, m, g, rp

    def _resout_adjoint(self, systeme, rhs, creux):
        """Une résolution transposée, contrôlée par erreur arrière équilibrée.

        Le raffinement réutilise la factorisation creuse. Il ne remplace pas
        un contrôle du conditionnement ni une validation du gradient physique.
        """
        if len(rhs) == 0:
            return rhs.copy()
        data = systeme.data if creux else systeme
        if not np.isfinite(data).all() or not np.isfinite(rhs).all():
            raise ValueError("adjoint : système ou second membre non fini")
        try:
            if creux:
                from scipy.sparse.linalg import splu
                eq = systeme.tocsc(copy=True)
                rows = abs(eq).max(axis=1).toarray().ravel()
                rows[rows == 0] = 1.
                eq.data /= rows[eq.indices]
                cols = abs(eq).max(axis=0).toarray().ravel()
                cols[cols == 0] = 1.
                eq.data /= np.repeat(cols, np.diff(eq.indptr))
                lu = splu(eq)
                solve = lambda b: lu.solve(b, trans='T')
            else:
                rows = np.max(abs(systeme), axis=1)
                rows[rows == 0] = 1.
                eq = systeme / rows[:, None]
                cols = np.max(abs(eq), axis=0)
                cols[cols == 0] = 1.
                eq = eq / cols[None, :]
                solve = lambda b: np.linalg.solve(eq.T, b)
            # Aeq = Dr A Dc ; Aeqᵀ y = Dc b et x = Dr y.
            balanced_rhs = rhs / cols
            if not np.isfinite(balanced_rhs).all():
                raise FloatingPointError("adjoint : équilibrage du second membre non représentable")
            sol = solve(balanced_rhs)
            matrix_norm = float(np.max(np.asarray(abs(eq).sum(axis=0))))
            self.n_resolutions_adjoint += 1
            for iteration in range(4):
                residual = balanced_rhs - eq.T @ sol
                if not np.isfinite(sol).all() or not np.isfinite(residual).all():
                    raise FloatingPointError("adjoint : résolution non finie")
                nx, nb = float(np.max(abs(sol))), float(np.max(abs(balanced_rhs)))
                scale = max(nx, nb)
                error = (float(np.max(abs(residual)))/scale / (matrix_norm*(nx/scale)+nb/scale)) if scale else 0.
                if error <= 1e-11:
                    self.residu_adjoint_max = max(self.residu_adjoint_max, error)
                    out = sol / rows
                    if not np.isfinite(out).all():
                        raise FloatingPointError("adjoint : solution physique non représentable")
                    return out
                if iteration < 3:
                    sol += solve(residual)
                    self.n_resolutions_adjoint += 1
        except (RuntimeError, np.linalg.LinAlgError) as exc:
            raise RuntimeError("adjoint : Jacobien du pas non inversible ; vérifier le rang des contraintes et le domaine du modèle") from exc
        raise FloatingPointError(f"adjoint : erreur arrière {error:.3g} > 1e-11")

    def _recul(self, N, fr0, fr1, ud0, a0, ud1, a1, p, mu):
        """(Aᵀμ, Bᵀμ) par le résidu implicite, sans former A ou B."""
        am, af, gam, bet = self.am, self.af, self.gam, self.bet
        c, ka, kb = (1-af)/(1-am), af/(1-am), am/(1-am)
        h = fr1[0]-fr0[0]
        if not np.isfinite(h) or h <= 0:
            raise ValueError("adjoint : pas temporel non positif ou non fini")
        nb = len(fr0[1]); n = 6*nb
        mu = np.asarray(mu, float)
        if mu.shape != (4*n,) or not np.isfinite(mu).all():
            raise ValueError("adjoint : covecteur final fini de forme (24n,) requis")
        creux = self._utilise_creux(n)
        u0 = np.concatenate((np.asarray(fr0[3]).reshape(nb,3),np.asarray(fr0[4]).reshape(nb,3)),axis=1).reshape(n)
        dq = u0 + (.5-bet)*h*a0 + bet*h*a1
        es = np.tile(np.eye(6), (nb, 1, 1))
        ds = es.copy()
        for i in range(nb):
            th = h*dq[6*i+3:6*i+6]
            es[i, 3:, 3:] = _expm(th)
            ds[i, 3:, 3:] = _jl(th) if self.jl else np.eye(3)
        def transpose(blocs, v):
            return (blocs.transpose(0,2,1) @ v.reshape(nb,6,1)).reshape(n)
        k, cc, m, g, rp = self._tangentes(N, fr1, ud1, p, creux, self.d_residu_transpose is None)
        pq, pu = bet*h*h*c, gam*h*c
        if creux:
            from scipy.sparse import bsr_matrix, bmat
            d = bsr_matrix((ds,np.arange(nb),np.arange(nb+1)),shape=(n,n)).tocsc()
            d.eliminate_zeros()
            sy = bmat([[pq*(k@d)+pu*cc+m,g.T],[pq*(g@d),None]],format='csc')
        else:
            # Produits par les blocs cinématiques, sans multiplication dense n³.
            kd = np.einsum('ibj,bjk->ibk', k.reshape(n,nb,6), ds).reshape(n,n)
            mm = g.shape[0]
            gd = np.einsum('ibj,bjk->ibk', g.reshape(mm,nb,6), ds).reshape(mm,n)
            sy = np.block([[pq*kd+pu*cc+m,g.T],[pq*gd,np.zeros((mm,mm))]])
        mq, mu_u, mud, ma = mu.reshape(4,n)
        rhs = np.concatenate((pq*transpose(ds,mq)+pu*mu_u+mud+c*ma,np.zeros(g.shape[0])))
        psi = self._resout_adjoint(sy,rhs,creux)
        vq = mq - k.T@psi[:n] - g.T@psi[n:]
        vu = mu_u - cc.T@psi[:n]
        dtq = transpose(ds,vq)
        previous = np.concatenate((transpose(es,vq), h*dtq+vu,
            bet*h*h*ka*dtq+gam*h*ka*vu+ka*ma,
            h*h*(.5-bet-bet*kb)*dtq+h*(1-gam-gam*kb)*vu-kb*ma))
        if self.d_residu_transpose is None:
            grad = -rp.T@psi[:n]
        else:
            grad = -np.asarray(self.d_residu_transpose(N, p, ud1, psi), float)
            if grad.shape != (len(p),):
                raise ValueError("d_residu_transpose doit fournir un vecteur de longueur nombre de paramètres")
        if not np.isfinite(previous).all() or not np.isfinite(grad).all():
            raise FloatingPointError("adjoint : produit transposé non fini")
        return previous, grad

    def _pas(self, N, fr0, fr1, ud0, a0, ud1, a1, p):
        """A = ∂z₁/∂z₀ (4n × 4n) et B = ∂z₁/∂p (4n × n_p) du pas fr0 → fr1."""
        am, af, gam, bet = self.am, self.af, self.gam, self.bet
        c, ka, kb = (1 - af) / (1 - am), af / (1 - am), am / (1 - am)
        h = fr1[0] - fr0[0]
        nb = len(fr0[1])
        n = 6 * nb
        u0 = np.concatenate([np.concatenate([fr0[3][i], fr0[4][i]]) for i in range(nb)])
        dq = u0 + (0.5 - bet) * h * a0 + bet * h * a1
        e = np.eye(n)
        d = np.eye(n)
        for i in range(nb):
            th = h * dq[6 * i + 3:6 * i + 6]
            s = slice(6 * i + 3, 6 * i + 6)
            e[s, s] = _expm(th)
            d[s, s] = _jl(th) if self.jl else np.eye(3)
        k, cc, m, g, rp = self._tangentes(N, fr1, ud1, p, historique=True)
        # δq₁ = Q0·X₀ + Pq·δu̇₁ ; δu₁ = U0·X₀ + Pu·δu̇₁ ; X₀ = (δq₀, δu₀, δu̇₀, δa₀)
        q0m = np.zeros((n, 4 * n))
        q0m[:, 0:n] = e
        q0m[:, n:2 * n] = h * d
        q0m[:, 2 * n:3 * n] = bet * h * h * ka * d
        q0m[:, 3 * n:4 * n] = h * h * (0.5 - bet - bet * kb) * d
        u0m = np.zeros((n, 4 * n))
        u0m[:, n:2 * n] = np.eye(n)
        u0m[:, 2 * n:3 * n] = gam * h * ka * np.eye(n)
        u0m[:, 3 * n:4 * n] = ((1 - gam) * h - gam * h * kb) * np.eye(n)
        pq, pu = bet * h * h * c * d, gam * h * c * np.eye(n)
        mm = g.shape[0]
        sy = np.zeros((n + mm, n + mm))
        sy[:n, :n] = k @ pq + cc @ pu + m
        sy[:n, n:] = g.T
        sy[n:, :n] = g @ pq
        rhs = np.zeros((n + mm, 4 * n + rp.shape[1]))
        rhs[:n, :4 * n] = -(k @ q0m + cc @ u0m)
        rhs[n:, :4 * n] = -(g @ q0m)
        rhs[:n, 4 * n:] = -rp
        sol = np.linalg.solve(sy, rhs)
        sx, sp = sol[:n, :4 * n], sol[:n, 4 * n:]
        a_ = np.zeros((4 * n, 4 * n))
        a_[0:n] = q0m + pq @ sx
        a_[n:2 * n] = u0m + pu @ sx
        a_[2 * n:3 * n] = sx
        a_[3 * n:4 * n] = c * sx
        a_[3 * n:4 * n, 2 * n:3 * n] += ka * np.eye(n)
        a_[3 * n:4 * n, 3 * n:4 * n] -= kb * np.eye(n)
        b_ = np.zeros((4 * n, rp.shape[1]))
        b_[0:n] = pq @ sp
        b_[n:2 * n] = pu @ sp
        b_[2 * n:3 * n] = sp
        b_[3 * n:4 * n] = c * sp
        return a_, b_

    def gradient(self, p, t_end, h, dj_dz_final, eps0=1e-6):
        """dJ/dp par UNE remontée ; `dj_dz_final(N, frame)` → ∂J/∂z_N (4n)."""
        p = np.atleast_1d(np.asarray(p, float))
        if p.ndim != 1 or not np.isfinite(p).all():
            raise ValueError("les paramètres doivent former un vecteur fini")
        if self.d_initial_transpose is None and (not np.isfinite(eps0) or eps0 <= 0):
            raise ValueError("eps0 doit être strictement positif et fini")
        N, fr, sch = self.aller(p, t_end, h)
        mu = np.asarray(dj_dz_final(N, fr[-1]), float)
        n = 6 * len(fr[0][1])
        if mu.shape != (4*n,) or not np.isfinite(mu).all():
            raise ValueError("dj_dz_final doit fournir un covecteur fini de forme (24n,)")
        grad = np.zeros(len(p))
        for i in range(len(fr) - 1, 0, -1):
            ud0, a0 = np.asarray(sch[i - 1][0]), np.asarray(sch[i - 1][1])
            ud1, a1 = np.asarray(sch[i][0]), np.asarray(sch[i][1])
            mu, contribution = self._recul(N, fr[i - 1], fr[i], ud0, a0, ud1, a1, p, mu)
            grad += contribution
        if self.d_initial_transpose is not None:
            initial = np.asarray(self.d_initial_transpose(p, mu), float)
            if initial.shape != (len(p),) or not np.isfinite(initial).all():
                raise ValueError("d_initial_transpose doit fournir un gradient initial fini de longueur nombre de paramètres")
            return grad + initial
        # État initial : q₀ et u₀ indépendants de p dans ce chemin historique.
        # u̇₀ = a₀ = acc_init(p), par différence finie sur le premier pas.
        for j in range(len(p)):
            uds = []
            for sg in (+1.0, -1.0):
                pp = p.copy()
                pp[j] += sg * eps0 * max(1.0, abs(p[j]))
                M_ = self.construit(pp)
                M_.enregistre_schema(True)
                M_.simule(h, h, tous=10 ** 9)
                uds.append(np.asarray(M_.schema()[0][0]))
            d_ud0 = (uds[0] - uds[1]) / (2 * eps0 * max(1.0, abs(p[j])))
            grad[j] += mu[2 * n:3 * n] @ d_ud0 + mu[3 * n:4 * n] @ d_ud0
        return grad


def _double_pendule(g):
    """Deux barres pivotées (axe y), lâchées à l'horizontale : rotations
    finies ET liaisons, le cas que le modèle NumPy ne pouvait pas porter."""
    from vinkulum import Noyau
    N = Noyau([0.0, 0.0, -float(np.atleast_1d(g)[0])])
    ji = [1e-3, 0, 0, 0, 2e-2, 0, 0, 0, 2e-2]
    a = N.corps("b1", 1.0, ji, [0.5, 0.0, 0.0])
    b = N.corps("b2", 0.7, ji, [1.5, 0.0, 0.0])
    N.liaison("p1", None, a, pa=[0.0, 0.0, 0.0], bloque_r=[0, 2])
    N.liaison("p2", a, b, pa=[0.5, 0.0, 0.0], bloque_r=[0, 2])
    return N


def _pendule_spherique(g):
    """Un corps en ROTULE, inertie asymétrique, lâché avec un spin hors du plan
    du balancement : rotations finies dans les TROIS axes. C'est le cas où
    J_l, ∂(J_s u̇)/∂θ et le gyroscopique pèsent — le pendule plan y est
    aveugle (rotation d'axe fixe, principal : les trois termes y sont nuls)."""
    from vinkulum import Noyau
    N = Noyau([0.0, 0.0, -float(np.atleast_1d(g)[0])])
    a = N.corps("b", 1.0, [1e-2, 0, 0, 0, 3e-2, 0, 0, 0, 2e-2], [0.5, 0.0, 0.0],
                w=[4.0, 0.0, 2.0])
    N.liaison("rotule", None, a, pa=[0.0, 0.0, 0.0], bloque_r=[])
    return N


def _pendule_flexible(ei, n_elem=3):
    """Chaîne de poutres pivotée à une extrémité, lâchée à l'horizontale :
    le paramètre est la raideur de flexion, ∂R/∂EI est EXACT (`d_residu_poutre`)."""
    from vinkulum import Noyau
    L, b, hh, rho, e_mod = 0.6, 0.02, 0.004, 2700.0, 70e9
    a_ = b * hh
    g_mod = e_mod / 2.6
    jt = b * hh ** 3 / 12 + hh * b ** 3 / 12
    dl = L / n_elem
    mel = rho * a_ * dl
    N = Noyau([0.0, 0.0, -9.81])
    idx = []
    for i in range(n_elem + 1):
        mm = mel * (0.5 if i in (0, n_elem) else 1.0)
        j = np.diag([rho * jt * dl, rho * (b * hh ** 3 / 12) * dl, rho * (hh * b ** 3 / 12) * dl])
        idx.append(N.corps(f"n{i}", mm, list(j.ravel()), [i * dl, 0.0, 0.0]))
    N.liaison("pivot", None, idx[0], pa=[0.0, 0.0, 0.0], bloque_r=[0, 2])
    for i in range(n_elem):
        N.poutre(f"b{i}", idx[i], idx[i + 1], e_mod * a_, g_mod * a_, g_mod * jt, float(ei))
    return N


def demo_pont(rapide=False):
    """Le pont contre une différence finie CENTRÉE sur le noyau lui-même."""
    import time as _t
    print("╔═ vinkulum — ADJOINT EN TEMPS SUR LE NOYAU (pont, premier pas)")
    out = {}

    def j_final(fr, corps, comp):
        n = 6 * len(fr[1])
        d = np.zeros(4 * n)
        d[6 * corps + comp] = 1.0
        return d

    # 1. double pendule, paramètre g : ∂R/∂g = −m sur les lignes z
    t_end, h = (1.5, 4e-3) if rapide else (3.0, 2e-3)

    def rp_g(N, p):
        r = np.zeros(12)
        for i, mm in enumerate((1.0, 0.7)):
            r[6 * i + 2] = mm            # R_t = m(u̇ − g⃗), g⃗ = −g e_z ⇒ ∂R/∂g = +m e_z
        return r[:, None]

    def j_dp(g):
        N = _double_pendule(g)
        N.simule(t_end, h, tous=10 ** 9)
        return N.etat()[1][1][0]                      # x du second corps à T

    t0 = _t.time()
    pont = PontNoyau(_double_pendule, rp_g)
    g_adj = pont.gradient([9.81], t_end, h, lambda N, fr: j_final(fr, 1, 0))[0]
    t_adj = _t.time() - t0
    t0 = _t.time()
    e = 1e-5
    g_df = (j_dp(9.81 + e) - j_dp(9.81 - e)) / (2 * e)
    t_df = _t.time() - t0
    ec = abs(g_adj / g_df - 1)
    print(f"║ double pendule, p = g, {int(round(t_end / h))} pas : adjoint {g_adj:+.9e}  DF {g_df:+.9e}"
          f"  écart {ec:.2e}")
    print(f"║   coût : {pont.n_kcmz} linéarisations, {t_adj:.1f} s, contre 2 simule = {t_df:.2f} s "
          f"pour la DF à 1 paramètre")
    out["double_pendule"] = dict(adjoint=float(g_adj), df=float(g_df), ecart=float(ec),
                                 t_adjoint=t_adj, t_df=t_df, n_kcmz=pont.n_kcmz)
    # 1b. pendule SPHÉRIQUE tournant : rotations finies sur trois axes. Le
    # pendule plan est AVEUGLE aux trois termes (mesuré : les retirer ne
    # change pas le 9ᵉ chiffre) — un contrôle négatif s'y serait cru inerte.
    def rp_g1(N, p):
        r = np.zeros(6)
        r[2] = 1.0
        return r[:, None]

    def j_sp(g):
        N = _pendule_spherique(g)
        N.simule(t_end, h, tous=10 ** 9)
        return N.etat()[1][0][0]

    ps = PontNoyau(_pendule_spherique, rp_g1)
    gs = ps.gradient([9.81], t_end, h, lambda N, fr: j_final(fr, 0, 0))[0]
    dfs = (j_sp(9.81 + e) - j_sp(9.81 - e)) / (2 * e)
    ecs = abs(gs / dfs - 1)
    print(f"║ pendule sphérique tournant, p = g : adjoint {gs:+.9e}  DF {dfs:+.9e}  écart {ecs:.2e}")
    out["spherique"] = dict(adjoint=float(gs), df=float(dfs), ecart=float(ecs))
    # CONTRÔLES NÉGATIFS : chaque pièce de la tangente retirée doit ÉLOIGNER de la DF
    neg = {}
    for nom, kw in (("sans J_l", dict(jl=False)), ("sans K_M", dict(k_m=False)),
                    ("gyro doublé (défaut d'avant)", dict(gyro=True))):
        pn = PontNoyau(_pendule_spherique, rp_g1, **kw)
        gn = pn.gradient([9.81], t_end, h, lambda N, fr: j_final(fr, 0, 0))[0]
        neg[nom] = abs(gn / dfs - 1)
        print(f"║   contrôle négatif {nom:26} : écart {neg[nom]:.2e}")
    out["negatifs"] = {k: float(v) for k, v in neg.items()}

    # 2. pendule flexible, paramètre EI : ∂R/∂EI exact par linéarité
    t_end2, h2 = (0.6, 2e-3) if rapide else (1.0, 1e-3)
    ei0 = 70e9 * 0.02 * 0.004 ** 3 / 12.0
    ne = 3

    def rp_ei(N, p):
        return sum(np.array(N.d_residu_poutre(j, "ei"), float) for j in range(ne))[:, None]

    def j_fl(ei):
        N = _pendule_flexible(ei, ne)
        N.simule(t_end2, h2, tous=10 ** 9)
        return N.etat()[1][ne][2]                     # z du bout à T

    pont2 = PontNoyau(lambda p: _pendule_flexible(p[0], ne), rp_ei)
    g2 = pont2.gradient([ei0], t_end2, h2, lambda N, fr: j_final(fr, ne, 2))[0]
    # LE JUGE A SA PROPRE RÉSOLUTION, et elle se mesure : J(EI) porte des
    # marches de ~1e-8 relatif là où le compte d'itérations de Newton change
    # (critère sur l'incrément, 1e-11), et ΔJ vaut ~2e-7 par 1e-4 de EI. La DF
    # varie donc de ~1 % avec ε (mesuré : 3,24e-4 → 3,31e-4 de ε 3e-3 à 3e-5),
    # et ce n'est pas le pont qui bouge. On publie la DF à plusieurs ε et on
    # juge l'adjoint DANS leur étendue — pas contre l'une d'elles.
    dfs2 = []
    for rel in (3e-4, 1e-4, 3e-5):
        e2 = rel * ei0
        dfs2.append((j_fl(ei0 + e2) - j_fl(ei0 - e2)) / (2 * e2))
    lo, hi = min(dfs2), max(dfs2)
    print(f"║ pendule flexible, p = EI, {int(round(t_end2 / h2))} pas : adjoint {g2:+.9e}  "
          f"DF {lo:+.6e} … {hi:+.6e} (ε de 3e-4 à 3e-5, étendue {100 * (hi - lo) / abs(hi):.2f} %)")
    marge = (hi - lo)
    ec2 = max(0.0, lo - g2, g2 - hi) / hi
    out["flexible"] = dict(adjoint=float(g2), df=[float(x) for x in dfs2], hors_etendue=float(ec2))
    assert ec < 1e-6, ("le pont s'écarte de la DF sur le double pendule", g_adj, g_df)
    assert ecs < 1e-4, ("le pont s'écarte de la DF sur le pendule sphérique", gs, dfs)
    assert lo - marge <= g2 <= hi + marge, ("le pont sort de l'étendue de la DF (flexible)", g2, dfs2)
    assert all(v > 10 * ecs for v in neg.values()), ("un terme retiré ne change rien — inerte", neg, ecs)
    print("╚═ pont adjoint ↔ noyau OK")
    return out


def demo(rapide=False):
    print("╔═ vinkulum — adjoint EN TEMPS : une remontée pour tous les paramètres")
    # oscillateur forcé : m ü + c u̇ + k u = p₀·sin(ωt) + p₁·cos(ωt)
    m, c, k = 1.0, 0.4, 100.0
    w = 7.0
    sys_ = Lineaire([[m]], [[c]], [[k]])
    p = np.array([2.0, -1.0])
    h, nst = 1e-3, 2000 if not rapide else 800

    def force(t, pp=p):
        return np.array([pp[0] * np.sin(w * t) + pp[1] * np.cos(w * t)])

    def d_force(t):
        return np.array([[np.sin(w * t), np.cos(w * t)]])

    # objectif : le DÉPLACEMENT final. ∂J/∂z = (1, 0, 0, 0) sur les 4n
    dj = np.zeros(4)
    dj[0] = 1.0
    g_adj = sys_.gradient(h, nst, force, d_force, dj)

    # référence : différences finies sur la trajectoire ENTIÈRE
    g_df = np.zeros(2)
    for i in range(2):
        vals = []
        for sens in (+1.0, -1.0):
            pp = p.copy()
            pp[i] += sens * 1e-5
            traj, _ = sys_.integre(h, nst, lambda t, q=pp: force(t, q))
            vals.append(traj[-1][1][0])
        g_df[i] = (vals[0] - vals[1]) / 2e-5

    # référence ANALYTIQUE : le système est linéaire, donc la réponse est
    # linéaire en p — le gradient EST la réponse au forçage unitaire
    g_an = np.zeros(2)
    for i in range(2):
        e = np.zeros(2)
        e[i] = 1.0
        traj, _ = sys_.integre(h, nst, lambda t, q=e: force(t, q))
        g_an[i] = traj[-1][1][0]

    print(f"║ oscillateur forcé · {nst} pas · h {h:.0e} · 2 paramètres d'amplitude")
    for i in range(2):
        print(f"║   p{i}   adjoint {g_adj[i]:+.9e}   DF {g_df[i]:+.9e}   "
              f"linéarité {g_an[i]:+.9e}")
    e_df = np.max(np.abs(g_adj - g_df) / np.maximum(np.abs(g_df), 1e-30))
    e_an = np.max(np.abs(g_adj - g_an) / np.maximum(np.abs(g_an), 1e-30))
    print(f"║ écart à la DF {100 * e_df:.2e} %   ·   à la référence analytique "
          f"{100 * e_an:.2e} %")
    print(f"║ coût : UNE remontée de {nst} pas pour les 2 paramètres ; "
          f"la DF en demande {2 * 2} trajectoires complètes")
    assert e_df < 1e-6, (g_adj, g_df)
    assert e_an < 1e-9, (g_adj, g_an)
    demo_nonlineaire()
    demo_pont(rapide)
    print("╚═ adjoint en temps OK")
    return dict(grad=g_adj.tolist(), ecart_df=float(e_df), ecart_analytique=float(e_an))


if __name__ == "__main__":
    demo(rapide="rapide" in sys.argv)
