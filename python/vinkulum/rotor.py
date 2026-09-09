"""vinkulum.rotor — rotor d'essai en articulé équivalent, et son TRIM.

    python -m vinkulum.rotor        # banc : stationnaire contre la théorie de
                                    # l'élément de pale, puis trim 3×3 en avancement

Premier banc AÉROMÉCANIQUE du noyau : jusqu'ici les bancs étaient mécaniques
(pendule, toupie, bielle, engrenage, poutre) ou venaient de FRELON. Celui-ci
est le squelette de tout cas rotorcraft opposable — même structure de commande
que HART-II ou n'importe quel essai de rotor isolé : trois commandes
(θ₀, θ₁c, θ₁s), trois cibles (poussée, moment de roulis, moment de tangage).

CONVENTIONS, écrites parce qu'elles décident du signe de tout :
· monde : +x = sens de vol, +z = haut ; rotor tournant autour de +z ;
· azimut ψ mesuré depuis l'ARRIÈRE dans le sens de rotation (usage
  hélicoptère) : l'envergure vaut e_s(ψ) = (−cos ψ, −sin ψ, 0), donc ψ = 90°
  est la pale AVANÇANTE ;
· le rotor est ENCASTRÉ et c'est l'air qui bouge : « avancer à V » s'écrit
  « vent de face à V » (`Noyau.vent`), ce qui donne μ = V/(ΩR) ;
· pas commandé θ(ψ) = θ₀ + θ₁c cos ψ + θ₁s sin ψ, appliqué par la CIBLE de la
  liaison de pale — donc par la cinématique, comme une timonerie sans jeu.

CE QUE LE MODÈLE EST, ET CE QU'IL N'EST PAS :
· pale RIGIDE sur charnière de battement centrale + ressort — l'articulé
  équivalent, réglé sur ν_β. C'est le modèle de rotor hingeless de la
  mécanique du vol, pas une pale flexible ; la poutre géométriquement exacte
  du noyau la remplacera, et c'est un chantier nommé, pas une omission ;
· pas de traînée (lag) : sans effet sur le trim en poussée et moments au
  premier ordre, déclaré ;
· inflow UNIFORME (Glauert en avancement). HOST classe ce niveau tout en bas
  de sa colonne dynamique — au-dessus viennent Pitt & Peters à 3 états, puis
  Peters–He à 15 ou 153. C'est le premier écart connu à un code de métier ;
· l'ORDRE des rotations n'est pas choisi à la main : la liaison bloque les
  composantes de pas et de traînée du résidu et laisse le battement libre, ce
  qui ne peut être satisfait exactement que par « battement PUIS pas » —
  l'ordre du palier réel, où le pitch bearing est en dehors de la charnière.
"""
import sys
import time
from math import cos, degrees, pi, radians, sin

import numpy as np

from vinkulum import Noyau
from vinkulum.trim import trim


def polaire_lineaire(a_lift=5.73, cl0=0.1, cd0=0.0079, k_cd=0.4, alpha_dec=14.0):
    """Polaire c81 mono-Mach : linéaire jusqu'au décrochage, puis plateau.

    Volontairement ANALYTIQUE et non tabulée d'un profil réel : le contrôle du
    banc est la théorie de l'élément de pale, qui suppose cl = a(α − α₀). Une
    table mesurée y ajouterait un écart de modèle qu'on ne saurait pas séparer
    de l'erreur du solveur. Le profil réel viendra avec le cas de corrélation.
    """
    al = np.concatenate([np.array([-180.0, -90.0, -30.0]),
                         np.arange(-alpha_dec, alpha_dec + 0.5, 1.0),
                         np.array([30.0, 90.0, 180.0])])
    a0 = -degrees(cl0 / a_lift)
    lin = lambda d: a_lift * radians(d - a0)
    cl_dec = lin(alpha_dec)
    cl = np.array([-0.0 if abs(d) > 90 else
                   (np.sign(d) * cl_dec * max(0.0, 1.0 - (abs(d) - alpha_dec) / 60.0)
                    if abs(d) > alpha_dec else lin(d)) for d in al])
    cd = np.array([cd0 + k_cd * radians(min(abs(d), 90.0)) ** 2 for d in al])
    cm = np.full_like(al, -0.01)
    return (al.tolist(), cl.tolist(), cd.tolist(), cm.tolist())


class Rotor:
    """Rotor d'essai : Nb pales, battement articulé équivalent, inflow uniforme."""

    def __init__(self, *, n_pales=4, R=2.0, r0=0.44, corde=0.121, omega=109.0,
                 vrillage=-8.0, r_vrillage=0.75, rho=1.2055, nu_beta=1.12,
                 lock=6.0, zeta=0.02, segments=5, gauss=5, polaire=None,
                 a_lift=5.73, cl0=0.1):
        self.__dict__.update(
            n_pales=n_pales, R=R, r0=r0, corde=corde, omega=omega, rho=rho,
            vrillage=vrillage, r_vrillage=r_vrillage, nu_beta=nu_beta,
            zeta=zeta, segments=segments, gauss=gauss)
        self.polaire = polaire or polaire_lineaire(a_lift=a_lift, cl0=cl0)
        self.a_lift, self.cl0 = a_lift, cl0
        self.alpha0 = -degrees(cl0 / a_lift)      # incidence de portance nulle
        self.sigma = n_pales * corde / (pi * R)
        self.aire = pi * R ** 2
        self.v_bout = omega * R
        # inertie de battement par le nombre de Lock (la MASSE de pale n'est
        # pas une donnée du banc : c'est γ qui pilote la réponse en battement)
        self.i_b = rho * self.a_lift * corde * R ** 4 / lock
        self.lock = lock
        self.m_b = self.i_b * 3.0 * (R - r0) / (R ** 3 - r0 ** 3)
        self.k_beta = self.i_b * omega ** 2 * (nu_beta ** 2 - 1.0)
        self.c_beta = 2.0 * zeta * self.i_b * omega * nu_beta
        self.t_tour = 2 * pi / omega

    # ── construction ────────────────────────────────────────────────────────
    def _reperes(self, psi):
        es = np.array([-cos(psi), -sin(psi), 0.0])
        ec = np.array([sin(psi), -cos(psi), 0.0])
        return es, ec, np.cross(es, ec)

    def monte(self, cmd, mu=0.0, alpha_disque=0.0, montee=1.0, t_end=None, v_i=None, harm=None):
        """Construit le rotor ; `cmd` = (θ₀, θ₁c, θ₁s) en degrés.

        `montee` = nombre de tours sur lesquels la commande passe de 0 à sa
        valeur (démarrage doux : la contrainte de pas est cinématique, une
        marche y injecterait une accélération infinie).
        """
        th0, th1c, th1s = (radians(v) for v in cmd)
        t_end = t_end or 4.0 * self.t_tour
        N = Noyau([0.0, 0.0, 0.0])          # rotor isolé : pas de pesanteur
        L = self.R - self.r0

        # LE ROTOR TOURNE DÉJÀ : les vitesses initiales doivent satisfaire la
        # contrainte en VITESSE (Φ̇ = 0), sinon le premier pas doit rattraper un
        # saut de 109 rad/s et Newton n'y arrive pas. Piège payé une fois.
        w0 = [0.0, 0.0, self.omega]
        moyeu = N.corps("moyeu", 1.0, list(np.diag([0.01, 0.01, 0.02]).ravel()), [0.0] * 3,
                        w=w0)
        N.liaison("arbre", None, moyeu, cible_r=([0.0, 0.0, 1.0], ("lineaire", [0.0, self.omega])))

        inflow = N.inflow([0.0, 0.0, 1.0], self.aire, 0.15 * self.t_tour, self.rho)
        if harm is not None:
            # inflow IMPOSÉ avec ses harmoniques (sillage libre, `vinkulum.sillage`) :
            # λ = λ₀ + r̄(λ1c cos ψ + λ1s sin ψ), ψ depuis l'arrière ⇒ e₁ = −x
            N.pose_inflow_harmoniques(inflow, float(harm[0]), float(harm[1]), float(harm[2]),
                                      [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], self.v_bout)
        elif v_i is not None:
            N.pose_inflow(inflow, float(v_i))

        self.idx = []
        for i in range(self.n_pales):
            psi0 = 2 * pi * i / self.n_pales
            es, ec, en = self._reperes(psi0)
            rot = np.column_stack([es, ec, en])
            # barre mince le long de es ; l'inertie de battement est celle du
            # banc (i_b), rapportée au CdM par Steiner
            r_cg = 0.5 * (self.r0 + self.R)
            j_env = self.m_b * self.corde ** 2 / 12.0
            j_bat = max(self.i_b - self.m_b * r_cg ** 2, 1e-6)
            j = np.diag([j_env, j_bat, j_bat + j_env])
            v0 = np.cross(np.array(w0), r_cg * es)
            b = N.corps(f"pale{i}", self.m_b, list(j.ravel()), list(r_cg * es), list(rot.ravel()),
                        v=list(v0), w=w0)
            self.idx.append(b)
            # pas imposé autour de l'envergure, traînée bloquée, battement libre
            N.liaison(f"pale{i}", moyeu, b, pa=[0.0] * 3, ra=list(rot.ravel()),
                      bloque_t=[0, 1, 2], bloque_r=[0, 2],
                      cible_r=([1.0, 0.0, 0.0], self._table_pas(psi0, th0, th1c, th1s, t_end, montee)))
            N.couple(f"kbeta{i}", moyeu, b, list(ec),
                     ("ressort", [self.k_beta, self.c_beta, 0.0]))
            # vrillage : un élément de pale par segment, repère de section TOURNÉ
            for k in range(self.segments):
                s0 = self.r0 + k * L / self.segments
                s1 = s0 + L / self.segments
                d = radians(self.vrillage * (0.5 * (s0 + s1) / self.R - self.r_vrillage))
                ec_k = [0.0, cos(d), sin(d)]        # repère CORPS : es=x, ec=y, en=z
                N.pale(f"p{i}s{k}", b, [(s0 - r_cg), 0.0, 0.0], [1.0, 0.0, 0.0], ec_k,
                       s1 - s0, self.corde, self.polaire, inflow=inflow,
                       rho=self.rho, gauss=self.gauss)

        v = mu * self.v_bout
        N.vent([-v * cos(alpha_disque), 0.0, -v * sin(alpha_disque)])
        self.t_end, self.inflow_idx = t_end, inflow
        return N

    def _table_pas(self, psi0, th0, th1c, th1s, t_end, montee):
        n = max(4, int(round(t_end / self.t_tour * 72)))
        out = []
        for k in range(n + 1):
            t = t_end * k / n
            psi = psi0 + self.omega * t
            amp = min(1.0, t / (montee * self.t_tour)) if montee > 0 else 1.0
            out += [t, amp * (th0 + th1c * cos(psi) + th1s * sin(psi))]
        return ("table", out)

    # ── mesure ──────────────────────────────────────────────────────────────
    def torseur(self, cmd, mu=0.0, alpha_disque=0.0, h=2e-4, tours=4.0, moyenne=1.0,
                v_i=None, harm=None):
        """(T, Mx, My) moyens sur le DERNIER tour, plus l'état du rotor.

        La moyenne sur un tour entier est ce qui fait de la sortie une
        constante : les efforts instantanés d'un rotor à Nb pales portent des
        harmoniques Nb/rev que le trim ne doit pas chasser.
        """
        t_end = tours * 2 * pi / self.omega
        N = self.monte(cmd, mu, alpha_disque, t_end=t_end, v_i=v_i, harm=harm)
        N.moyenne_aero(t_end - moyenne * self.t_tour)
        N.simule(t_end, h, tous=10 ** 9)
        f, m, _ = N.torseur_moyen()
        beta = []
        for nom, th, w, tau in N.couples():
            if nom.startswith("kbeta"):
                beta.append(degrees(th))
        return dict(T=f[2], Mx=m[0], My=m[1], Fx=f[0], Fy=f[1],
                    v_i=N.aero()[1][self.inflow_idx][0], beta=beta, N=N)

    def trim(self, cibles, x0=(8.0, 0.0, 0.0), mu=0.0, alpha_disque=0.0, harm=None, **kw):
        """Trim 3×3 : (θ₀, θ₁c, θ₁s) en degrés pour (T, Mx, My) visés.
        `harm` = (v_i, v1c, v1s) en m/s impose l'inflow d'un sillage libre."""
        cib = np.array([cibles["T"], cibles.get("Mx", 0.0), cibles.get("My", 0.0)])
        etat = {}

        def residu(x):
            r = self.torseur(tuple(x), mu, alpha_disque, v_i=etat.get("v_i"), harm=harm, **kw)
            etat["v_i"] = r["v_i"]          # repart de l'inflow trouvé : un
            etat["dernier"] = r             # transitoire de moins par évaluation
            return np.array([r["T"], r["Mx"], r["My"]]) - cib

        ech = np.array([max(abs(cib[0]), 1.0), max(abs(cib[0]) * self.R * 0.05, 1.0),
                        max(abs(cib[0]) * self.R * 0.05, 1.0)])
        x, r, info = trim(residu, np.array(x0, float), pas=0.05, tol=1e-3, echelle=ech)
        return x, info, etat["dernier"]

    # ── référence indépendante ──────────────────────────────────────────────
    def ct_bet(self, theta0_deg, lam, mu=0.0, th1s_deg=0.0):
        """C_T par la théorie de l'élément de pale, inflow uniforme λ.

        C_T = (σ a / 2) ∫ ((θ(r̄) − α₀) r̄² − λ r̄) dr̄ sur [r₀/R, 1], vrillage
        compris. C'est la référence ANALYTIQUE du banc en stationnaire : elle ne
        partage avec le solveur ni la cinématique, ni l'intégration, ni le
        battement. L'incidence de portance nulle α₀ EN FAIT PARTIE — l'oublier
        décalait le verdict d'exactement 1,0° de collectif, c'est-à-dire de la
        totalité de l'écart qu'on croyait mesurer.
        """
        x0 = self.r0 / self.R
        th0 = radians(theta0_deg - self.alpha0)
        tw = radians(self.vrillage)
        # θ(r̄) = θ0 + tw·(r̄ − r_vrillage)
        # ⟨u_T²⟩ = r̄² + μ²/2 et ⟨u_P u_T⟩ = λ r̄ sur un tour (le coning ne
        # contribue pas en moyenne, ⟨cos ψ⟩ = ⟨cos ψ sin ψ⟩ = 0).
        # ⚠ LE CYCLIQUE, LUI, CONTRIBUE : ⟨sin ψ · u_T²⟩ = μ r̄ ≠ 0, donc le
        # LATÉRAL θ₁s entre dans la poussée MOYENNE, à hauteur de μ. L'oublier
        # donnait +8,5 % d'écart en avancement — le solveur était juste.
        m2 = mu ** 2 / 2
        i_th = (th0 * ((1 - x0 ** 3) / 3 + m2 * (1 - x0))
                + radians(th1s_deg) * mu * (1 - x0 ** 2) / 2
                + tw * ((1 - x0 ** 4) / 4 - self.r_vrillage * (1 - x0 ** 3) / 3
                        + m2 * ((1 - x0 ** 2) / 2 - self.r_vrillage * (1 - x0))))
        return 0.5 * self.sigma * self.a_lift * (i_th - lam * (1 - x0 ** 2) / 2)

    def theta_bet(self, T, mu=0.0, lam=None, th1s_deg=0.0):
        """θ₀ analytique pour une poussée T (inverse de `ct_bet`).

        En stationnaire λ vient de Froude et le contrôle est TOTALEMENT
        indépendant du solveur. En avancement, λ est celui que le solveur a
        trouvé : le contrôle porte alors sur la cinématique et l'intégration,
        pas sur le modèle d'inflow — et c'est dit, parce qu'un contrôle dont on
        ignore ce qu'il partage avec le calculé ne juge rien.
        """
        ct = T / (self.rho * self.aire * self.v_bout ** 2)
        if lam is None:
            lam = (ct / 2.0) ** 0.5
        x0 = self.r0 / self.R
        tw = radians(self.vrillage)
        m2 = mu ** 2 / 2
        i_tw = (tw * ((1 - x0 ** 4) / 4 - self.r_vrillage * (1 - x0 ** 3) / 3
                      + m2 * ((1 - x0 ** 2) / 2 - self.r_vrillage * (1 - x0)))
                + radians(th1s_deg) * mu * (1 - x0 ** 2) / 2)
        th0 = (ct / (0.5 * self.sigma * self.a_lift) + lam * (1 - x0 ** 2) / 2 - i_tw) \
            / ((1 - x0 ** 3) / 3 + m2 * (1 - x0))
        return degrees(th0) + self.alpha0, lam, ct


def demo(rapide=False):
    r = Rotor()
    print(f"rotor d'essai : {r.n_pales} pales · R {r.R} m · c {r.corde} m · "
          f"Ω {r.omega} rad/s · σ {r.sigma:.4f} · γ {r.lock} · ν_β {r.nu_beta}")
    t_cible = 3300.0

    # ── 1. STATIONNAIRE : le trim atteint sa cible, et θ₀ tombe sur la théorie
    t0 = time.time()
    x_h, info_h, det_h = r.trim(dict(T=t_cible), x0=(9.0, 0.0, 0.0), mu=0.0,
                                tours=3.0 if rapide else 5.0)
    dt_h = time.time() - t0
    th_bet, lam_bet, ct = r.theta_bet(t_cible)
    ecart = (x_h[0] - th_bet) / th_bet
    print(f"\n  STATIONNAIRE, cible T = {t_cible:.0f} N   (C_T {ct:.5f} · C_T/σ {ct / r.sigma:.4f})")
    print(f"    trim      θ₀ {x_h[0]:7.3f}°   θ₁c {x_h[1]:6.3f}°   θ₁s {x_h[2]:6.3f}°")
    print(f"    obtenu    T {det_h['T']:8.1f} N   Mx {det_h['Mx']:7.2f}   My {det_h['My']:7.2f} N·m")
    print(f"    v_i {det_h['v_i']:.3f} m/s (Froude {(t_cible / (2 * r.rho * r.aire)) ** 0.5:.3f})"
          f"   ·   battement {min(det_h['beta']):.3f}…{max(det_h['beta']):.3f}°")
    print(f"    THÉORIE de l'élément de pale : θ₀ {th_bet:.3f}° (λ {lam_bet:.5f}"
          f" · α₀ {r.alpha0:+.3f}°)   ·   écart {100 * ecart:+.2f} %")
    print(f"    coût : {info_h['evals']} évaluations, {info_h['matrices']} matrice(s), {dt_h:.1f} s")

    assert info_h["converge"], info_h
    assert abs(det_h["T"] - t_cible) / t_cible < 2e-3, det_h["T"]
    # le stationnaire est AXISYMÉTRIQUE : ni cyclique, ni battement, ni moment
    assert abs(x_h[1]) < 0.05 and abs(x_h[2]) < 0.05, x_h
    assert max(det_h["beta"]) - min(det_h["beta"]) < 0.02, det_h["beta"]
    # référence indépendante : le BET est une borne de MODÈLE, pas une vérité —
    # il linéarise cl, néglige l'inclinaison de la portance par φ et la
    # composante de traînée. Quelques pour-cent sont l'écart attendu ; au-delà,
    # c'est le solveur qu'il faut regarder, et non le seuil.
    assert abs(ecart) < 0.03, (x_h[0], th_bet, ecart)

    if rapide:
        print("\nrotor OK (rapide)")
        return

    # ── 2. AVANCEMENT : trim 3×3, la dissymétrie doit être compensée
    for mu in (0.15,):
        t0 = time.time()
        x, info, det = r.trim(dict(T=t_cible, Mx=0.0, My=0.0), x0=(x_h[0], 0.0, 0.0),
                              mu=mu, tours=5.0)
        dt = time.time() - t0
        print(f"\n  AVANCEMENT μ = {mu}   (V {mu * r.v_bout:.1f} m/s)")
        print(f"    trim      θ₀ {x[0]:7.3f}°   θ₁c {x[1]:6.3f}°   θ₁s {x[2]:6.3f}°")
        print(f"    obtenu    T {det['T']:8.1f} N   Mx {det['Mx']:7.2f}   My {det['My']:7.2f} N·m")
        print(f"    battement {min(det['beta']):.2f}…{max(det['beta']):.2f}°   ·   "
              f"v_i {det['v_i']:.3f} m/s   ·   H {det['Fx']:.1f} N")
        print(f"    coût : {info['evals']} évaluations, {info['matrices']} matrice(s), {dt:.1f} s")
        lam_mes = det["v_i"] / r.v_bout
        th_av, _, _ = r.theta_bet(t_cible, mu=mu, lam=lam_mes, th1s_deg=x[2])
        ec_av = (x[0] - th_av) / th_av
        print(f"    THÉORIE de l'élément de pale à μ (λ mesuré {lam_mes:.5f}, "
              f"θ₁s repris) : θ₀ {th_av:.3f}°   ·   écart {100 * ec_av:+.2f} %")
        assert info["converge"], info
        assert abs(det["T"] - t_cible) / t_cible < 3e-3, det["T"]
        assert abs(det["Mx"]) < 1.0 and abs(det["My"]) < 1.0, (det["Mx"], det["My"])
        assert abs(ec_av) < 0.05, (x[0], th_av, ec_av)
        # ce que le stationnaire ne pouvait pas montrer : il FAUT du cyclique
        assert max(abs(x[1]), abs(x[2])) > 0.5, x
        # et le disque reste PLAT — ν_β ≈ 1 lie le moment de moyeu au battement
        # 1/rev (M = K_β·β₁), donc trimer les moments à zéro annule β₁ : la
        # dispersion azimutale du battement doit rester sous le coning
        disp = max(det["beta"]) - min(det["beta"])
        assert disp < 0.1 * abs(np.mean(det["beta"])), (disp, det["beta"])
        print(f"    disque plat : dispersion du battement {disp:.3f}° pour un coning "
              f"de {abs(np.mean(det['beta'])):.2f}° — c'est M = K_β·β₁ avec M ≈ 0")
    print("\nrotor OK")


if __name__ == "__main__":
    demo(rapide="rapide" in sys.argv)
