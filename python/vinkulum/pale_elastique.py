"""vinkulum.pale_elastique — la PALE ÉLASTIQUE couplée à l'aérodynamique.

    python -m vinkulum.pale_elastique            # bancs complets
    python -m vinkulum.pale_elastique rapide

Jusqu'ici la pale du noyau était RIGIDE sur charnière équivalente
(`vinkulum.rotor`) : c'est le modèle de mécanique du vol, pas une pale. Ici
la pale est une chaîne de poutres géométriquement exactes (Simo–Reissner,
`Noyau.poutre`) dont chaque nœud porte son tronçon d'aérodynamique
(`Noyau.pale`) : la section suit la ROTATION du nœud, donc la torsion
élastique change l'incidence, le battement change la vitesse vue — le
couplage aéroélastique est dans la cinématique, pas ajouté.

CE QUI EST JUGÉ, contre des références RECONNUES :
· le diagramme en éventail d'une poutre encastrée TOURNANTE (raidissement
  centrifuge), contre la table de Wright, Smith, Thresher & Wang (1982) —
  identique aux valeurs de Hodges & Rutkowski (1981) — jusqu'à Ω̄ = 12 ;
· l'amortissement aérodynamique du premier mode de battement en
  stationnaire contre le nombre de Lock GÉNÉRALISÉ au mode (Johnson,
  *Helicopter Theory*, ch. 9) — le « γ/16 » de la pale rigide, écrit sur la
  forme modale ;
· la cohérence poussée / flèche : la flèche du bout de pale du noyau contre
  une poutre intégrée en Python sous les MÊMES charges (aéro lues sur le
  noyau + centrifuge) — contrôle d'élément, pas de physique nouvelle ;
· en avancement, la réponse 1/rev de la pale élastique contre la pale rigide
  de même inertie — publiée, périodicité assertée.

MULTITHREAD, mesuré : un `ThreadPoolExecutor` sur des `Noyau` indépendants
rend 0,88× sous le GIL de 3.13 quand le noyau garde le GIL ; depuis le
5 sept. `simule`, `modes_complexes` et `k_c_m_z` le LIBÈRENT (`allow_threads`)
— le gain se mesure dans `bancs`, et le parallélisme interne (rayon) reste.

CE QUI N'EST PAS LÀ, déclaré : pas de couplage flexion-torsion par décalage
du centre de gravité (sections symétriques ici, le noyau le porterait par
`p0`) ; masse ajoutée non circulatoire absente (module `pale`) ; inflow
uniforme ou imposé, pas de sillage libre dans cette boucle ; la référence
« Wright 1982 » est prise de mémoire de la table classique — les valeurs
sont dans `WRIGHT_1982` avec leur provenance, à confronter au papier.
"""
import sys
import time
from math import cos, pi, radians, sin, sqrt

import numpy as np

from vinkulum import Noyau
from vinkulum.rotor import polaire_lineaire

# Premier mode de flexion d'une poutre uniforme encastrée tournante, sans
# déport de pied, ω̄ = ω√(ρA L⁴/EI) en fonction de Ω̄ = Ω√(ρA L⁴/EI).
# Wright, Smith, Thresher, Wang, J. Appl. Mech. 49 (1982) ; Hodges &
# Rutkowski, AIAA J. 19 (1981). Ω̄ = 0 : 1,8751² = 3,5160 (exact).
WRIGHT_1982 = {0.0: 3.5160, 1.0: 3.6816, 2.0: 4.1373, 3.0: 4.7973, 4.0: 5.5850,
               5.0: 6.4495, 6.0: 7.3604, 10.0: 11.2023, 12.0: 13.1702}


class PaleElastique:
    """Rotor à pales ÉLASTIQUES : nœuds en corps, poutres GE, aéro par nœud."""

    def __init__(self, *, n_pales=1, R=1.0, r0=0.0, corde=0.08, omega=50.0,
                 rho_l=0.6, ei_flap=40.0, ei_lag=800.0, gj=60.0, ea=2e6,
                 n_el=12, theta0=6.0, rho=1.225, polaire=None, a_lift=5.73, cl0=0.0):
        self.__dict__.update(n_pales=n_pales, R=R, r0=r0, corde=corde, omega=omega,
                             rho_l=rho_l, ei_flap=ei_flap, ei_lag=ei_lag, gj=gj, ea=ea,
                             n_el=n_el, theta0=theta0, rho=rho, a_lift=a_lift, cl0=cl0)
        self.polaire = polaire or polaire_lineaire(a_lift=a_lift, cl0=cl0)
        self.L = R - r0
        self.le = self.L / n_el
        self.aire = pi * R * R
        self.t_tour = 2 * pi / omega

    def _es(self, psi):
        es = np.array([-cos(psi), -sin(psi), 0.0])
        ec = np.array([sin(psi), -cos(psi), 0.0])
        return es, ec, np.cross(es, ec)

    def monte(self, aero=True, mu=0.0, rigide=False, precontrainte=True, omega=None, v_i=None):
        """Construit le rotor. `rigide` multiplie les raideurs par 1e4 (pale de
        référence). `precontrainte` place les nœuds à leur allongement
        centrifuge, pour partir d'un ÉQUILIBRE (sinon le mode axial sonne)."""
        om = self.omega if omega is None else omega
        f = 1e4 if rigide else 1.0
        N = Noyau([0.0, 0.0, 0.0])
        w0 = [0.0, 0.0, om]
        moyeu = N.corps("moyeu", 1.0, list(np.diag([0.01, 0.01, 0.02]).ravel()), [0.0] * 3, w=w0)
        N.liaison("arbre", None, moyeu, cible_r=([0.0, 0.0, 1.0], ("lineaire", [0.0, om])))
        inflow = N.inflow([0.0, 0.0, 1.0], self.aire, 0.15 * self.t_tour, self.rho) if aero else None
        if aero and v_i is not None:
            N.pose_inflow(inflow, float(v_i))
        th = radians(self.theta0)
        self.noeuds = []
        m_el = self.rho_l * self.le
        for i in range(self.n_pales):
            psi0 = 2 * pi * i / self.n_pales
            es, ec, en = self._es(psi0)
            # repère de nœud : es, puis (ec, en) tournés du pas θ₀ autour de es
            ecp = cos(th) * ec + sin(th) * en
            enp = -sin(th) * ec + cos(th) * en
            rot = np.column_stack([es, ecp, enp])
            noeuds = []
            for k in range(self.n_el + 1):
                r = self.r0 + k * self.le
                m = m_el if 0 < k < self.n_el else 0.5 * m_el
                jl = m * self.le * self.le / 12.0
                jj = np.diag([m * self.corde ** 2 / 12.0, jl, jl])
                v0 = np.cross(np.array(w0), r * es)
                noeuds.append(N.corps(f"p{i}n{k}", m, list(jj.ravel()), list(r * es),
                                      list(rot.ravel()), v=list(v0), w=w0))
            # encastrement au moyeu (hingeless), pas θ₀ dans la géométrie
            N.liaison(f"enc{i}", moyeu, noeuds[0], pa=[0.0] * 3, ra=list(rot.ravel()))
            for k in range(self.n_el):
                N.poutre(f"p{i}e{k}", noeuds[k], noeuds[k + 1], self.ea * f, 1e2 * self.ea * f,
                         self.gj * f, self.ei_flap * f, ei3=self.ei_lag * f, ga3=1e2 * self.ea * f)
            if aero:
                for k in range(self.n_el + 1):
                    a0 = -0.5 * self.le if k > 0 else 0.0
                    a1 = 0.5 * self.le if k < self.n_el else 0.0
                    N.pale(f"a{i}n{k}", noeuds[k], [a0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                           a1 - a0, self.corde, self.polaire, inflow=inflow, rho=self.rho, gauss=3)
            self.noeuds.append(noeuds)
        if precontrainte:
            # LA POUTRE PREND SA LONGUEUR AU REPOS À LA CRÉATION : les nœuds sont
            # créés non étirés, puis POSÉS à leur allongement centrifuge —
            # créés déjà étirés, la tension était nulle et l'éventail restait
            # plat (mesuré : ω̄₁ 3,50 à Ω̄ 12 pour 13,17). N(r) = ½ρΩ²(R² − r²)
            # ⇒ δ(r) = ρΩ²/(EA)·(R²r/2 − r³/6), rapporté au pied.
            t, rr, rot, v, w, vi = N.etat()
            rr = [list(x) for x in rr]
            d = lambda x: self.rho_l * om * om / (self.ea * f) * (self.R ** 2 * x / 2 - x ** 3 / 6)
            for i in range(self.n_pales):
                es = self._es(2 * pi * i / self.n_pales)[0]
                for k, n in enumerate(self.noeuds[i]):
                    r = self.r0 + k * self.le
                    rr[n] = list((r + d(r) - d(self.r0)) * es)
            N.pose_etat(t, rr, [list(x) for x in rot], [list(x) for x in v], [list(x) for x in w], list(vi))
        if mu:
            N.vent([-mu * om * self.R, 0.0, 0.0])
        self.inflow_idx = inflow
        return N

    # ── mesures ─────────────────────────────────────────────────────────────
    def fleche_bout(self, N, i=0):
        """Déplacement hors plan (en) du bout de pale, repère du moyeu."""
        t, r, rot, v, w, vi = N.etat()
        rh = np.array(rot[0]).reshape(3, 3)
        p = rh.T @ np.array(r[self.noeuds[i][-1]])
        return float(p[2])

    def poussee(self, N):
        return sum(p[1] for p in N.aero()[0])


# ── bancs ───────────────────────────────────────────────────────────────────
def eventail(rapide=False):
    """Diagramme en éventail : poutre encastrée tournante contre Wright (1982)."""
    n_el = 12 if rapide else 20
    pe = PaleElastique(n_el=n_el, rho_l=0.6, ei_flap=40.0, ei_lag=40.0, gj=60.0, ea=2e6, R=1.0)
    ref = sqrt(pe.rho_l * pe.L ** 4 / pe.ei_flap)
    print("╔═ vinkulum — PALE ÉLASTIQUE : éventail d'une poutre tournante contre Wright et al. (1982)")
    print(f"║ {n_el} éléments GE, masse concentrée, précontrainte centrifuge posée à la main")
    out, pires = {}, 0.0
    cles = sorted(WRIGHT_1982) if not rapide else [0.0, 2.0, 6.0, 12.0]
    for ob in cles:
        om = ob / ref
        N = pe.monte(aero=False, omega=om if om > 0 else 1e-9)
        mc = N.modes_complexes(combien=6)
        f = [x[0] for x in mc if x[0] > 1e-3]
        # le premier mode HORS PLAN : à Ω = 0 flap et traînée coïncident (EI égaux), on prend le premier
        w_nd = 2 * pi * f[0] * ref
        e = w_nd / WRIGHT_1982[ob] - 1.0
        pires = max(pires, abs(e))
        out[ob] = (w_nd, WRIGHT_1982[ob], e)
        print(f"║ Ω̄ {ob:4.1f} : ω̄₁ noyau {w_nd:8.4f}  table {WRIGHT_1982[ob]:8.4f}  écart {100 * e:+.3f} %")
    print(f"╚═ pire écart {100 * pires:.3f} %")
    return out, pires


def amortissement_lock(rapide=False):
    """Amortissement aéro du premier mode de battement contre le Lock généralisé."""
    pe = PaleElastique(n_el=10 if rapide else 16, omega=60.0, theta0=0.0, cl0=0.0)
    N = pe.monte(aero=False)
    mo = N.modes(combien=3)
    f0, phi = mo[0]
    # forme modale en déplacement hors plan (z monde) des nœuds de la pale
    nz = [phi[6 * n + 2] for n in pe.noeuds[0]]
    eta = np.abs(np.array(nz)) / max(abs(x) for x in nz)
    r = np.array([pe.r0 + k * pe.le for k in range(pe.n_el + 1)])
    m = np.array([pe.rho_l * pe.le * (1.0 if 0 < k < pe.n_el else 0.5) for k in range(pe.n_el + 1)])
    i_eta = float(np.sum(m * eta ** 2))
    c_aero = 0.5 * pe.rho * pe.a_lift * pe.corde * pe.omega * float(np.trapezoid(r * eta ** 2, r))
    zeta_th = c_aero / (2 * i_eta * 2 * pi * f0)
    # avec l'aéro : modes complexes autour de l'état tournant
    Na = pe.monte(aero=True, v_i=0.0)
    Na.pose_inflow(pe.inflow_idx, 0.0, impose=True)
    mc = Na.modes_complexes(combien=6)
    # le mode de battement : celui dont la fréquence est la plus proche de f0
    fa, za, sa = min(mc, key=lambda x: abs(x[0] - f0))
    print("╔═ vinkulum — PALE ÉLASTIQUE : amortissement aéro du battement contre le Lock généralisé")
    print(f"║ mode de battement {f0:.3f} Hz (sans aéro) → {fa:.3f} Hz avec aéro, ζ noyau {za:.4f}, "
          f"ζ théorie (γ modal) {zeta_th:.4f}, écart {100 * (za / zeta_th - 1):+.1f} %")
    print("╚═ (forme modale du noyau, ∫rη² dr trapèzes, portance quasi-stationnaire)")
    return dict(f0=f0, fa=fa, zeta=za, zeta_th=zeta_th)


def stationnaire(rapide=False):
    """Poussée et flèche : élastique contre rigide, et contre une poutre intégrée."""
    pe = PaleElastique(n_el=10 if rapide else 14, omega=60.0, theta0=8.0)
    tours = 6.0 if rapide else 10.0
    h = pe.t_tour / 400
    res = {}
    for nom, rig in (("rigide", True), ("élastique", False)):
        N = pe.monte(aero=True, rigide=rig)
        N.simule(tours * pe.t_tour, h, tous=10 ** 9)
        res[nom] = (pe.poussee(N), pe.fleche_bout(N), N.etats_inflow(pe.inflow_idx)[0], N)
    # poutre intégrée sous les charges DU NOYAU (aéro par nœud + centrifuge)
    N = res["élastique"][3]
    fz = np.zeros(pe.n_el + 1)
    for nom, t, q in N.aero()[0]:
        k = int(nom.split("n")[-1])
        fz[k] += t
    r = np.array([pe.r0 + k * pe.le for k in range(pe.n_el + 1)])
    # charge linéique : la poussée de chaque nœud sur SON tronçon (l_e, moitié aux bouts)
    lk = np.array([pe.le if 0 < k < pe.n_el else 0.5 * pe.le for k in range(pe.n_el + 1)])
    # flèche statique d'une poutre tendue : EI w'''' − (N w')' = f_z, N(r) = ½ρΩ²(R²−r²)
    n_c = 0.5 * pe.rho_l * pe.omega ** 2 * (pe.R ** 2 - r ** 2)
    nf = 200
    x = np.linspace(pe.r0, pe.R, nf + 1)
    dx = x[1] - x[0]
    fzl = np.interp(x, r, fz / lk)
    ncl = np.interp(x, r, n_c)
    # différences finies : A w = f, encastrée en 0, libre en L
    A = np.zeros((nf + 1, nf + 1))
    b = np.zeros(nf + 1)
    for i in range(2, nf - 1):
        A[i, i - 2:i + 3] += pe.ei_flap / dx ** 4 * np.array([1, -4, 6, -4, 1])
        A[i, i - 1:i + 2] -= np.array([ncl[i] - 0.5 * (ncl[i + 1] - ncl[i - 1]) / 2,
                                       -2 * ncl[i], ncl[i] + 0.5 * (ncl[i + 1] - ncl[i - 1]) / 2]) / dx ** 2
        b[i] = fzl[i]
    A[0, 0] = 1.0; A[1, 0], A[1, 1] = -1.0, 1.0                      # w = 0, w' = 0
    A[nf - 1, nf - 3:nf + 1] = np.array([-1, 3, -3, 1]) / dx ** 3   # w''' = 0 (effort tranchant)
    A[nf, nf - 2:nf + 1] = np.array([1, -2, 1]) / dx ** 2           # w'' = 0 (moment)
    w = np.linalg.solve(A, b)
    fl_int = float(w[-1])
    tr, te = res["rigide"][0], res["élastique"][0]
    print("╔═ vinkulum — PALE ÉLASTIQUE en stationnaire : poussée et flèche")
    print(f"║ poussée rigide {tr:.4f} N · élastique {te:.4f} N ({100 * (te / tr - 1):+.2f} %) · "
          f"v_i {res['élastique'][2]:.3f} m/s")
    print(f"║ flèche du bout : noyau {1e3 * res['élastique'][1]:.3f} mm · poutre tendue intégrée "
          f"{1e3 * fl_int:.3f} mm ({100 * (res['élastique'][1] / fl_int - 1):+.2f} %) · rigide {1e3 * res['rigide'][1]:.4f} mm")
    print("╚═")
    return dict(t_rigide=tr, t_elastique=te, fleche=res["élastique"][1], fleche_int=fl_int)


def avancement(rapide=False):
    """μ = 0,2 : réponse 1/rev du battement élastique, périodicité."""
    pe = PaleElastique(n_el=8 if rapide else 12, omega=60.0, theta0=8.0)
    tours = 8.0 if rapide else 12.0
    h = pe.t_tour / 400
    out = {}
    for nom, rig in (("rigide", True), ("élastique", False)):
        N = pe.monte(aero=True, rigide=rig, mu=0.2)
        traj = N.simule(tours * pe.t_tour, h, tous=4)
        t = np.array([f[0] for f in traj])
        z = np.array([f[1][pe.noeuds[0][-1]][2] for f in traj])
        m1 = z[(t > (tours - 2) * pe.t_tour) & (t <= (tours - 1) * pe.t_tour)]
        m2 = z[t > (tours - 1) * pe.t_tour]
        n = min(len(m1), len(m2))
        per = float(np.max(np.abs(m1[:n] - m2[:n]))) / max(np.ptp(m2), 1e-12)
        out[nom] = (float(np.mean(m2)), 0.5 * float(np.ptp(m2)), per, pe.poussee(N))
    print("╔═ vinkulum — PALE ÉLASTIQUE en avancement μ = 0,2")
    for nom, (mz, amp, per, T) in out.items():
        print(f"║ {nom:10} : flèche moyenne {1e3 * mz:+.3f} mm · amplitude 1/rev {1e3 * amp:.3f} mm · "
              f"écart entre les deux derniers tours {100 * per:.2f} % · poussée {T:.4f} N")
    print("╚═")
    return out


def demo(rapide=False):
    t0 = time.time()
    ev, pire = eventail(rapide)
    # · l'éventail : la table classique au demi pour-cent (n_el = 20) ; 1 % en rapide
    assert pire < (0.01 if rapide else 0.005), ev
    am = amortissement_lock(rapide)
    # · l'amortissement aéro : le Lock généralisé au mode, à 15 % (forme modale
    #   et intégration trapèzes ; la pale rigide donnerait γ/16 exactement)
    assert abs(am["zeta"] / am["zeta_th"] - 1.0) < 0.15, am
    st = stationnaire(rapide)
    # · la flèche du noyau et la poutre tendue intégrée sous les mêmes charges
    assert abs(st["fleche"] / st["fleche_int"] - 1.0) < 0.05, st
    # · la poussée ne bouge PAS : section symétrique sans décalage de CG, pas
    #   de couplage flexion-torsion, et un battement ÉTABLI ne change pas
    #   l'incidence — mesuré −0,11 % à 10 éléments, +0,23 % à 14 : c'est du
    #   bruit de discrétisation autour de zéro, pas un signe (le premier
    #   assert disait « porte moins » ; c'était une croyance)
    assert abs(st["t_elastique"] / st["t_rigide"] - 1.0) < 0.01, st
    av = avancement(rapide)
    # · périodique sur l'élastique (la rigide bat de 0,01 mm : sa périodicité
    #   relative est du bruit, publiée)
    assert av["élastique"][2] < 0.05, av
    print(f"╚═ pale élastique OK — {time.time() - t0:.0f} s")
    return dict(eventail={str(k): v for k, v in ev.items()}, amortissement=am, stationnaire=st,
                avancement={k: list(v) for k, v in av.items()})


if __name__ == "__main__":
    demo(rapide="rapide" in sys.argv)
