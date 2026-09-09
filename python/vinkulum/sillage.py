"""vinkulum.sillage — le SILLAGE LIBRE d'un rotor en stationnaire.

    python -m vinkulum.sillage            # banc complet (~3 min)
    python -m vinkulum.sillage rapide     # banc court (~30 s)

Le noyau connaît deux inflows : l'uniforme (théorie du disque) et Pitt–Peters
(trois états). Tous deux POSENT la forme du sillage. Ici le sillage est
CALCULÉ : le tourbillon marginal de chaque pale est un filament libre, convecté
par la vitesse que TOUT le sillage induit sur lui — et la contraction, le pas
hélicoïdal, l'inflow au disque en SORTENT au lieu d'y entrer.

CE QUE C'EST :
· N_b filaments (un tourbillon marginal par pale), discrétisés en âge ζ, N_lib
  tours LIBRES relaxés + N_loin tours PRESCRITS (hélice rigide prolongeant le
  dernier tour libre, même rayon, même pas) — la troncature est déclarée, et le
  prolongement rend au disque le downwash du sillage lointain ;
· Biot–Savart par segment droit, noyau de Vatistas n = 2 (r_c = fraction de
  corde, `r_core` — c'est un BOUTON, pas une donnée : aucune mesure ici ne
  le fixe) ;
· équation d'âge dans le repère TOURNANT : x(ζ+Δζ) = R_z(−Δζ)·(x(ζ) +
  v̄·Δζ/Ω), résolue comme un POINT FIXE par Newton amorti sur les 3K
  inconnues (r, z, δφ) d'un filament, les N_b filaments étant identiques à
  une rotation près. Aucune marche en temps, aucune relaxation de Picard —
  celle-ci diverge (instabilité d'appariement, mesurée, cf. `resout`) ;
· la circulation liée est UNIFORME sur l'envergure et vient de la poussée :
  Γ = 2T/(ρ N_b Ω R²) (pale rectangulaire non vrillée, portance ∝ r).

MESURÉ SANS NAPPE (rotor FRELON Ø350, 2 pales, C_T 0,0069, Δζ 15°, 3 tours
libres — le modèle du 5 sept., gardé en contrôle négatif) :
    Newton 6 itérations, 1,1e-2 → 3,9e-10 ; λ +1,8 % du momentum ; r_v(2π)
    0,800 contre Landgrebe 0,807 (−0,9 %), asymptote 0,72 (Landgrebe 0,78) ;
    k₂ à −48 % de Landgrebe, k₁ de MAUVAIS SIGNE ; Δζ/2 : λ +2,3 %, r_v
    −0,2 % ; 2 tours libres au lieu de 3 : +0,3 % ; couplage au noyau : v_i à
    +0,7 % de la théorie du disque en 10 allers-retours.
AVEC NAPPE (le défaut depuis le 6 sept.) : k₁ −0,0198 (bon signe, Landgrebe
    −0,0238), k₂ −14 %, r_v(2π) 0,837 (+3,6 %), asymptote 0,797, λ +4,2 % ;
    Δζ/2 : λ +4,0 %, r_v −1,5 % ; 2 tours : −0,2 % ; couplage −0,6 %.

LA NAPPE INTERNE (6 sept.) : Γ(r) de la théorie de l'élément de pale au pas
qui rend la poussée (inflow de Froude), `n_nappe` éléments, filaments
internes ΔΓ_k et tourbillon de pied PRESCRITS en hélice de rayon fixe (pas
d'âge doublé, tours libres + 3), au pas de l'INFLOW LOCAL v(r_k) lu au disque
après un premier Newton au pas de Froude — deux passages. C'est le compromis
CAMRAD, nappe prescrite + marginal libre. Ce qu'elle rend, mesuré (15°,
3 tours) : k₁ +0,0026 → −0,0198 (Landgrebe −0,0238, facteur 1,2), k₂
−0,043 → −0,071 (−0,083), contraction 0,800 → 0,837 (Landgrebe 0,807),
asymptote 0,777 → 0,797 (0,78), λ 1,018 → 1,042 × momentum, couplage noyau
−0,6 % ; 20 → ~80 s. Essayé et refusé, mesuré : la nappe au pas de FROUDE
seul (k₁ −0,016, mais sur un rotor chargé — Maryland, C_T/σ 0,145 — ses
filaments restaient au ras du disque et y induisaient un upwash intérieur :
FM à +15/+19 % de la mesure au lieu de −3/−13, cf. `validation`) ; un
quasi-Newton au jacobien SANS la nappe (40 it, plancher 2e-7 : le marginal
jeune est trop près de la nappe) ; la nappe FINE (même pas d'âge, toute la
longueur : 142 s pour k₁ −0,011, moins bon que la grossière) ; 4 éléments
(contraction à +8 % de celle sans nappe). Libérer la nappe n'a pas été
essayé (×6 sur les inconnues). `nappe=False` rend le modèle d'avant, au
chiffre près.

CE QUE CE N'EST PAS, déclaré :
· la nappe est PRESCRITE (pas de contraction, pas de déformation) et Γ(r)
  vient d'une pale non vrillée à pas uniforme, pas de la polaire ;
· pas de vrillage, pas d'avancement (μ = 0 seulement), pas de dissipation du
  noyau avec l'âge ;
· le disque n'est pas un domaine de calcul : l'inflow « moyen » est une
  moyenne pondérée par l'aire sur un échantillon d'anneaux.
"""
import sys
import time

import numpy as np

RHO = 1.225


def _rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def biot_savart(pts, seg_a, seg_b, gam, r_core):
    """Vitesse induite en `pts` (P,3) par des segments droits A→B (S,3) de
    circulation `gam` (S,), noyau de Vatistas n = 2.

    Formule exacte du segment (r₁ + r₂)/(r₁r₂(r₁r₂ + r₁·r₂))·(r₁ × r₂), fois le
    facteur de noyau h²/√(r_c⁴ + h⁴), h = distance au segment. Un point sur le
    prolongement d'un segment (r₁ × r₂ = 0) rend zéro, sans division par zéro.
    """
    r1 = pts[:, None, :] - seg_a[None, :, :]
    r2 = pts[:, None, :] - seg_b[None, :, :]
    n1 = np.linalg.norm(r1, axis=2)
    n2 = np.linalg.norm(r2, axis=2)
    cr = np.cross(r1, r2)
    dl = np.linalg.norm(seg_b - seg_a, axis=1)[None, :]
    h2 = np.sum(cr * cr, axis=2) / np.maximum(dl * dl, 1e-30)
    den = n1 * n2 * (n1 * n2 + np.sum(r1 * r2, axis=2))
    coef = np.where(den > 1e-30, (n1 + n2) / np.maximum(den, 1e-30), 0.0)
    rc = np.asarray(r_core, float).reshape(1, -1) if np.ndim(r_core) else r_core
    noyau = h2 / np.sqrt(rc ** 4 + h2 * h2)
    w = gam[None, :] * coef * noyau / (4.0 * np.pi)
    return np.einsum("ps,psk->pk", w, cr)


class Sillage:
    """Sillage libre d'un rotor à `n_b` pales, rayon `R`, régime `omega`."""

    def __init__(self, R, n_b, omega, corde, r_pied=0.0, dzeta_deg=12.0,
                 tours_libres=4, tours_loin=6, r_core=None, delta=100.0,
                 nappe=True, n_nappe=6, a_lift=5.73):
        self.R, self.n_b, self.omega, self.corde = R, n_b, omega, corde
        self.r_pied = r_pied
        # NAPPE INTERNE (5–6 sept.) : Γ(r) de la théorie de l'élément de pale au
        # pas qui donne la poussée (inflow de Froude), discrétisée en `n_nappe`
        # éléments ; les filaments internes (ΔΓ aux frontières, tourbillon de
        # pied) sont PRESCRITS en hélice à rayon fixe et au pas de Froude —
        # le compromis des codes de métier (CAMRAD : nappe prescrite, marginal
        # libre). `nappe=False` = Γ uniforme, un seul filament : le modèle
        # d'avant, exactement.
        self.nappe, self.n_nappe, self.a_lift = nappe, n_nappe, a_lift
        self.r_st = np.array([r_pied, R])      # stations de la nappe
        self.g_el = np.array([0.0])           # Γ par élément (le dernier = marginal)
        self.dz = np.radians(dzeta_deg)
        self.k_lib = int(round(tours_libres * 2 * np.pi / self.dz))
        self.k_loin = int(round(tours_loin * 2 * np.pi / self.dz))
        # ponytail: r_c = 0,2 c est un bouton de robustesse numérique (0,05–0,4 c
        # dans la littérature) ; il n'est calibré sur aucune mesure ici
        self.r_core = 0.2 * corde if r_core is None else r_core
        # CROISSANCE DU NOYAU AVEC L'ÂGE (Squire ; Bhagwat–Leishman 2002) :
        # r_c² = r_c0² + 4 α δ ν ζ/Ω, α = 1,25643, ν = 1,5e-5, δ = coefficient
        # de viscosité turbulente (10–1000 dans la littérature). C'est ce qui
        # amortit l'instabilité d'appariement des tours lointains ; δ est un
        # BOUTON déclaré, non mesuré ici.
        self.delta = delta
        self.psi_b = 2 * np.pi * np.arange(n_b) / n_b
        self.gamma = 0.0
        self._v_nappe = 0.0
        self.x = None          # (n_b, k_lib+1, 3) nœuds libres, repère tournant
        self.residus = []

    # ── géométrie ─────────────────────────────────────────────────────────
    def circulation(self, poussee):
        """Γ depuis la poussée. Sans nappe : uniforme, T = ρ Γ Ω N_b R²/2.
        Avec nappe : Γ(r) = ½ a c Ω (θ r − λR)⁺, λ = √(C_T/2), θ résolu pour
        que ρ Ω N_b ∫ Γ r dr = T ; le marginal porte le Γ du dernier élément."""
        if not self.nappe or poussee <= 0.0:
            self.gamma = 2.0 * max(poussee, 0.0) / (RHO * self.n_b * self.omega * self.R ** 2)
            self.r_st, self.g_el = np.array([self.r_pied, self.R]), np.array([self.gamma])
            self.theta = 0.0
            return self.gamma
        vt = self.omega * self.R
        lam = np.sqrt(max(poussee, 0.0) / (RHO * np.pi * self.R ** 2 * vt * vt) / 2.0)
        r = np.linspace(self.r_pied, self.R, 400)
        g = lambda th: np.maximum(0.5 * self.a_lift * self.corde * self.omega * (th * r - lam * self.R), 0.0)
        t_de = lambda th: RHO * self.omega * self.n_b * np.trapezoid(g(th) * r, r)
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if t_de(mid) < poussee else (lo, mid)
        self.theta = 0.5 * (lo + hi)
        self.r_st = np.linspace(self.r_pied, self.R, self.n_nappe + 1)
        rm = 0.5 * (self.r_st[:-1] + self.r_st[1:])
        self.g_el = np.maximum(0.5 * self.a_lift * self.corde * self.omega * (self.theta * rm - lam * self.R), 0.0)
        self.gamma = float(self.g_el[-1])
        return self.gamma

    def _nappe(self, v_i):
        """Filaments internes PRESCRITS (A, B, Γ, r_c) : à chaque frontière r_k
        (pied compris, bout exclu — il est le marginal libre) une hélice de
        rayon r_k, pas v_i/Ω, force Γ_{k−1} − Γ_k (Kelvin ; au pied 0 − Γ₀)."""
        if not self.nappe:
            return None
        # segments deux fois plus longs que le marginal et longueur bornée à
        # (tours libres + 3) : la nappe pèse 6 filaments, c'est elle qui coûte
        dzn = 2.0 * self.dz
        n_seg = int(round((self.k_lib * self.dz + 3 * 2 * np.pi) / dzn))
        z = np.arange(n_seg + 1) * dzn
        rc = np.sqrt(self.r_core ** 2 + 4 * 1.25643 * self.delta * 1.5e-5 * (z[:-1] + 0.5 * dzn) / self.omega)
        A, B, G, C = [], [], [], []
        g_ext = np.concatenate([[0.0], self.g_el])
        v_k = np.broadcast_to(np.asarray(v_i, float), (len(self.r_st),))
        # LA NAPPE SE CONTRACTE AVEC LE MARGINAL : rayon r_k · r_v(ζ)/R, r_v lu sur
        # le marginal LIBRE (premier passage). À rayon fixe, le filament externe
        # de la nappe (0,84 R) se retrouvait DEHORS du marginal contracté (0,80 R
        # à 2π) et sa vorticité, de sens opposé, soufflait vers le HAUT au bord du
        # disque — mesuré sur Maryland : FM à +25 % au lieu de −3.
        if self.x is not None and self.x.shape[1] == self.k_lib + 1:
            _, rv, _ = self.marginal()
            ratio = np.interp(z, np.arange(self.k_lib + 1) * self.dz, rv)
        else:
            ratio = np.ones_like(z)
        for k in range(len(self.r_st) - 1):
            dg = g_ext[k] - g_ext[k + 1]
            if dg == 0.0:
                continue
            for b in range(self.n_b):
                az = self.psi_b[b] - z
                pts = np.column_stack([self.r_st[k] * ratio * np.cos(az), self.r_st[k] * ratio * np.sin(az),
                                       -v_k[k] * z / self.omega])
                A.append(pts[:-1]); B.append(pts[1:]); G.append(np.full(n_seg, dg)); C.append(rc)
        if not A:
            return None
        return np.vstack(A), np.vstack(B), np.concatenate(G), np.concatenate(C)

    def helice_rigide(self, v_i):
        """Hélice non déformée : rayon R, pas axial v_i/Ω par radian d'âge."""
        z = np.arange(self.k_lib + 1) * self.dz
        x = np.zeros((self.n_b, self.k_lib + 1, 3))
        for b in range(self.n_b):
            az = self.psi_b[b] - z
            x[b, :, 0] = self.R * np.cos(az)
            x[b, :, 1] = self.R * np.sin(az)
            x[b, :, 2] = -v_i * z / self.omega
        return x

    def _prolonge(self, x):
        """Sillage lointain PRESCRIT : chaque filament continue en hélice rigide
        au rayon et au pas de son dernier tour libre."""
        n_t = int(round(2 * np.pi / self.dz))
        out = []
        for b in range(self.n_b):
            last, prev = x[b, -1], x[b, -1 - n_t]
            r_last = np.hypot(last[0], last[1])
            dz_tour = last[2] - prev[2]
            az0 = np.arctan2(last[1], last[0])
            k = np.arange(1, self.k_loin + 1)
            az = az0 - k * self.dz
            far = np.column_stack([r_last * np.cos(az), r_last * np.sin(az),
                                   last[2] + dz_tour * k * self.dz / (2 * np.pi)])
            out.append(np.vstack([x[b], far]))
        return np.array(out)

    def segments(self, x=None):
        """Tous les segments (A, B, Γ) : sillage libre + lointain + tourbillons liés."""
        x = self.x if x is None else x
        full = self._prolonge(x)
        a = full[:, :-1, :].reshape(-1, 3)
        b = full[:, 1:, :].reshape(-1, 3)
        # tourbillon LIÉ : du pied au bout, par élément de la nappe, Γ_k (Kelvin :
        # le dernier se prolonge dans le marginal, qui porte Γ dans le sens de l'âge)
        lie_a, lie_b, lie_g = [], [], []
        for p in self.psi_b:
            e = np.array([np.cos(p), np.sin(p), 0.0])
            for k in range(len(self.r_st) - 1):
                lie_a.append(self.r_st[k] * e); lie_b.append(self.r_st[k + 1] * e); lie_g.append(self.g_el[k])
        lie_a, lie_b, lie_g = np.array(lie_a), np.array(lie_b), np.array(lie_g)
        z_seg = (np.arange(full.shape[1] - 1) + 0.5) * self.dz
        rc = np.sqrt(self.r_core ** 2 + 4 * 1.25643 * self.delta * 1.5e-5 * z_seg / self.omega)
        rc = np.concatenate([np.tile(rc, self.n_b), np.full(len(lie_g), self.r_core)])
        A, B = np.vstack([a, lie_a]), np.vstack([b, lie_b])
        G = np.concatenate([np.full(len(a), self.gamma), lie_g])
        nap = self._nappe(self._v_nappe)
        if nap is not None:
            A, B, G, rc = (np.vstack([A, nap[0]]), np.vstack([B, nap[1]]),
                           np.concatenate([G, nap[2]]), np.concatenate([rc, nap[3]]))
        return A, B, G, rc

    def vitesse(self, pts, x=None):
        """Vitesse induite (P,3) en repère monde par tout le système tourbillonnaire."""
        a, b, g, rc = self.segments(x)
        return biot_savart(np.atleast_2d(pts), a, b, g, rc)

    def _marche(self, x):
        """Une convection complète depuis le bout, sur le champ induit par `x` :
        x(ζ+Δζ) = R_z(−Δζ)·(x(ζ) + v̄·Δζ/Ω), v̄ = moyenne des vitesses aux deux
        nœuds (trapèze sur le champ figé)."""
        v = self.vitesse(x.reshape(-1, 3), x).reshape(x.shape)
        new = np.empty_like(x)
        new[:, 0] = x[:, 0]
        rot = _rz(-self.dz)
        for k in range(self.k_lib):
            vk = 0.5 * (v[:, k] + v[:, k + 1])
            new[:, k + 1] = (x[:, k] + vk * self.dz / self.omega) @ rot.T
        return new

    def _de_u(self, u):
        """Géométrie SYMÉTRIQUE depuis les inconnues (r, z, δφ) par âge k ≥ 1 :
        les N_b filaments sont identiques à une rotation près — exact en
        stationnaire, et c'est ce qui réduit les inconnues à 3K."""
        u = u.reshape(-1, 3)
        x = np.empty((self.n_b, self.k_lib + 1, 3))
        k = np.arange(1, self.k_lib + 1)
        for b in range(self.n_b):
            x[b, 0] = [self.R * np.cos(self.psi_b[b]), self.R * np.sin(self.psi_b[b]), 0.0]
            a = self.psi_b[b] - k * self.dz + u[:, 2]
            x[b, 1:, 0] = u[:, 0] * np.cos(a)
            x[b, 1:, 1] = u[:, 0] * np.sin(a)
            x[b, 1:, 2] = u[:, 1]
        return x

    def resout(self, v_i0=None, tol=1e-9, iters=40, libre=True):
        """Le point fixe x = marche(x), résolu par NEWTON AMORTI (jacobien par
        différences finies sur les 3K inconnues (r, z, δφ) par âge, moindres
        carrés, rebroussement sur le résidu). Convergence QUADRATIQUE mesurée :
        1,1e-2 → 1,8e-3 → 7e-6 → 1,6e-10 en 7 itérations.

        LA RELAXATION DE PICARD A ÉTÉ ESSAYÉE ET REFUSÉE, mesurée : sous-relaxée
        à 0,1–0,4, avec ou sans croissance du noyau (δ 100–500), avec 2 à 4 tours
        libres, symétrisée ou non, lissée ou non, le résidu MONTE (1,4e-2 →
        2e-2…1e-1 en 300 itérations) et les tours 2–4 remontent — c'est
        l'instabilité d'appariement des hélices, physique en stationnaire, que
        Picard reproduit au lieu de la traverser. Le stationnaire est une
        ÉQUATION, pas une marche en temps : Newton la résout que la solution
        soit stable ou non. (`scipy.optimize.root` hybr, essayé aussi : n'a pas
        quitté le point de départ — non instruit, le Newton maison suffit.)
        `libre=False` garde l'hélice rigide (contrôle négatif)."""
        if v_i0 is None:
            v_i0 = np.sqrt(self.gamma * self.n_b * self.omega / (4 * np.pi)) if self.gamma > 0 else 0.0
        self._v_nappe = v_i0                  # à chaque appel : Froude, puis le pas LOCAL (deux passages)
        if self.x is None or self.x.shape[1] != self.k_lib + 1:
            self.x = self.helice_rigide(v_i0)
        if not libre or self.gamma == 0.0:
            self.residus, self.succes = [0.0], True
            return self.x
        if self.nappe:
            # LE PAS DE LA NAPPE EST CELUI DE L'INFLOW LOCAL, PAS DE FROUDE (6 sept.) :
            # au pas de Froude, les filaments internes d'un rotor chargé (Maryland,
            # C_T/σ 0,145) restaient trop près du plan du disque et y induisaient un
            # upwash intérieur — FM à +15/+19 % de la mesure au lieu de −3/−13.
            # Deux passages : Newton au pas de Froude, lecture de v(r_k) au disque
            # (moyenne en azimut), Newton au pas local depuis la géométrie trouvée.
            self._resout_marginal(tol, iters)
            az = (np.arange(24) + 0.5) / 24 * 2 * np.pi
            rk = np.maximum(self.r_st, 1e-3 * self.R)
            pts = np.column_stack([np.outer(rk, np.cos(az)).ravel(), np.outer(rk, np.sin(az)).ravel(),
                                   np.zeros(rk.size * az.size)])
            vz = self.vitesse(pts)[:, 2].reshape(rk.size, az.size).mean(axis=1)
            self._v_nappe = np.maximum(-vz, 0.0)
        return self._resout_marginal(tol, iters)

    def _resout_marginal(self, tol, iters):
        """Newton sur le marginal, nappe (si présente) figée."""
        x0 = self.x[0, 1:]
        k = np.arange(1, self.k_lib + 1)
        u = np.column_stack([np.hypot(x0[:, 0], x0[:, 1]), x0[:, 2],
                             np.angle(np.exp(1j * (np.arctan2(x0[:, 1], x0[:, 0]) - self.psi_b[0]
                                                    + k * self.dz)))]).ravel()
        ech = np.tile([self.R, self.R, 1.0], self.k_lib)
        u = u / ech

        def res(v):
            x = self._de_u(v * ech)
            return ((self._marche(x) - x)[0, 1:] / self.R).ravel()

        n = len(u)
        r = res(u)
        self.residus = [float(np.abs(r).max())]
        self.succes = False
        for _ in range(iters):
            if self.residus[-1] < tol:
                self.succes = True
                break
            # (un quasi-Newton au jacobien SANS la nappe a été essayé : 40 itérations,
            # plancher 2e-7, 87 s contre 142 — le marginal jeune est trop près de
            # la nappe pour qu'on l'ignore dans la dérivée. Jacobien complet.)
            jac = np.empty((n, n))
            for j in range(n):
                e = np.zeros(n)
                e[j] = 1e-7
                jac[:, j] = (res(u + e) - r) / e[j]
            du = np.linalg.lstsq(jac, -r, rcond=None)[0]
            al = 1.0
            for _ in range(8):
                r_e = res(u + al * du)
                if np.abs(r_e).max() < self.residus[-1]:
                    break
                al *= 0.5
            u, r = u + al * du, r_e
            self.residus.append(float(np.abs(r).max()))
        self.x = self._de_u(u * ech)
        return self.x

    # ── sorties ───────────────────────────────────────────────────────────
    def marginal(self, b=0):
        """(ζ, r_v/R, z_v/R) du tourbillon marginal de la pale `b`."""
        z = np.arange(self.k_lib + 1) * self.dz
        x = self.x[b]
        return z, np.hypot(x[:, 0], x[:, 1]) / self.R, x[:, 2] / self.R

    def inflow_disque(self, n_r=12, n_az=24):
        """Vitesse induite moyenne à travers le disque (m/s, > 0 vers le bas),
        pondérée par l'aire sur `n_r` anneaux × `n_az` azimuts."""
        r = (np.arange(n_r) + 0.5) / n_r * self.R
        az = (np.arange(n_az) + 0.5) / n_az * 2 * np.pi
        rr, aa = np.meshgrid(r, az, indexing="ij")
        pts = np.column_stack([(rr * np.cos(aa)).ravel(), (rr * np.sin(aa)).ravel(),
                               np.zeros(rr.size)])
        w = np.repeat(r, n_az)
        vz = self.vitesse(pts)[:, 2]
        return float(-np.sum(vz * w) / np.sum(w))


def landgrebe(zeta, c_t, sigma, n_b, theta_tw=0.0):
    """Landgrebe (1972) : position du tourbillon marginal, empirique.

    r_v/R = A + (1 − A)e^{−λζ}, A = 0,78, λ = 0,145 + 27 C_T ;
    z_v/R = k₁ζ pour ζ ≤ 2π/N_b, puis k₁(2π/N_b) + k₂(ζ − 2π/N_b),
    k₁ = −0,25(C_T/σ + 0,001 θ_tw), k₂ = −(1,41 + 0,0141 θ_tw)√(C_T/2).
    """
    lam = 0.145 + 27.0 * c_t
    r = 0.78 + 0.22 * np.exp(-lam * zeta)
    k1 = -0.25 * (c_t / sigma + 0.001 * theta_tw)
    k2 = -(1.41 + 0.0141 * theta_tw) * np.sqrt(c_t / 2.0)
    zb = 2 * np.pi / n_b
    z = np.where(zeta <= zb, k1 * zeta, k1 * zb + k2 * (zeta - zb))
    return r, z, k1, k2


# ── couplage au noyau ─────────────────────────────────────────────────────────
def rotor_noyau(R, n_b, omega, corde, theta_deg, r_pied, tau=1e9):
    """Rotor rigide dans le noyau : un moyeu entraîné à Ω, `n_b` pales à pas
    fixe, un inflow uniforme dont le retard `tau` (énorme) GÈLE la relaxation
    interne — c'est le sillage qui posera v_i, le noyau ne doit pas le défaire."""
    from vinkulum import Noyau
    from vinkulum.rotor import polaire_lineaire
    N = Noyau([0.0, 0.0, 0.0])
    moyeu = N.corps("moyeu", 1.0, list(np.diag([0.01, 0.01, 0.02]).ravel()), [0.0] * 3,
                    w=[0.0, 0.0, omega])
    N.liaison("arbre", None, moyeu, cible_r=([0.0, 0.0, 1.0], ("lineaire", [0.0, omega])))
    i = N.inflow([0.0, 0.0, 1.0], np.pi * R * R, tau, RHO)
    th = np.radians(theta_deg)
    pol = polaire_lineaire()
    for b in range(n_b):
        p = 2 * np.pi * b / n_b
        es = np.array([np.cos(p), np.sin(p), 0.0])
        e_psi = np.array([-np.sin(p), np.cos(p), 0.0])
        ec = np.cos(th) * e_psi + np.sin(th) * np.array([0.0, 0.0, 1.0])
        N.pale(f"pale{b}", moyeu, list(es * r_pied), list(es), list(ec), R - r_pied, corde,
               pol, inflow=i, rho=RHO)
    return N, i


def couple(N, i, R, n_b, omega, corde, r_pied=0.0, h=1e-4, tol=0.01, iters=12, **kw):
    """Couple le sillage au noyau : poussée du noyau → Γ → sillage → v_i moyen
    → `pose_inflow` → un tour de rotor → … jusqu'à ce que v_i bouge de moins
    de `tol` en relatif. Rend (v_i sillage, v_i théorie du disque à la même
    poussée, poussée, itérations)."""
    s = Sillage(R, n_b, omega, corde, r_pied=r_pied, **kw)
    t_tour = 2 * np.pi / omega
    # départ à la théorie du disque sur la poussée du premier tour : trois
    # allers-retours de moins (chacun un Newton de sillage, 5 s avec la nappe)
    N.simule(N.t() + t_tour, h, tous=10 ** 9)
    p0 = sum(p[1] for p in N.aero()[0])
    N.pose_inflow(i, float(np.sqrt(max(p0, 0.0) / (2 * RHO * np.pi * R * R))))
    v_i, k = 0.0, 0
    for k in range(1, iters + 1):
        N.simule(N.t() + t_tour, h, tous=10 ** 9)
        poussee = sum(p[1] for p in N.aero()[0])
        s.circulation(poussee)
        s.resout()                       # Newton repart de la géométrie précédente
        v_new = s.inflow_disque()
        N.pose_inflow(i, v_new)
        if abs(v_new - v_i) < tol * max(v_new, 1e-9):
            v_i = v_new
            break
        v_i = v_new
    v_disque = np.sqrt(max(poussee, 0.0) / (2 * RHO * np.pi * R * R))
    return v_i, v_disque, poussee, k, s


# ── AVANCEMENT ────────────────────────────────────────────────────────────────
class SillageAvancement:
    """Sillage libre en VOL D'AVANCEMENT (μ > 0), en repère du moyeu.

    En vol stabilisé la géométrie est PÉRIODIQUE : l'élément lâché par la pale
    à l'azimut ψ_r se retrouve, à l'âge ζ, en r(ψ_r, ζ), et le filament de la
    pale k à l'azimut de pale ψ_b est {r(ψ_b − 2πk/N_b − ζ, ζ)}. L'inconnue est
    donc la GRILLE r(ψ_r, ζ) (N_ψ lâchers × N_ζ âges), et non plus les 3K
    inconnues d'un seul filament comme en stationnaire. Convection :

        r(ψ_r, ζ+Δζ) = r(ψ_r, ζ) + [V_air + ½(v(ψ_r,ζ) + v(ψ_r,ζ+Δζ))]·Δζ/Ω

    v étant la vitesse induite par TOUT le sillage tel qu'il est quand la
    pale est à ψ_b = ψ_r + ζ. Résolution par RELAXATION (marche en âge sur le
    champ figé, sous-relaxée) : c'est la relaxation qui DIVERGE en
    stationnaire (`Sillage.resout`), et qui converge ici — l'avancement
    emporte les tours l'un loin de l'autre et tue l'instabilité
    d'appariement. Mesuré, cf. `demo_avancement`. Circulation Γ(ψ_r)
    = Γ₀(1 + g_c cos ψ_r + g_s sin ψ_r) : uniforme par défaut (rotor trimé ⇒
    portance ~constante en azimut), une 1/rev déclarée en option.

    Conventions = celles de `vinkulum.rotor` : ψ compté depuis l'ARRIÈRE, bout
    de pale en (−R cos ψ, −R sin ψ, 0), rotation +z, air en (−μΩR, 0, 0).
    Harmoniques d'inflow : λ = λ₀ + r̄(λ1c cos ψ + λ1s sin ψ), e₁ = arrière.
    """

    def __init__(self, R, n_b, omega, corde, mu, r_pied=0.0, dpsi_deg=15.0,
                 tours_libres=3, tours_loin=4, r_core=None, delta=100.0,
                 gamma_harm=(0.0, 0.0), alpha_disque=0.0):
        self.R, self.n_b, self.omega, self.corde, self.mu = R, n_b, omega, corde, mu
        self.r_pied = r_pied
        self.n_psi = int(round(360.0 / dpsi_deg))
        if self.n_psi % n_b:
            raise ValueError("N_ψ doit être multiple de N_b (décalage entier entre pales)")
        self.dz = 2 * np.pi / self.n_psi
        self.n_z = int(round(tours_libres * self.n_psi))
        self.n_loin = int(round(tours_loin * self.n_psi))
        self.r_core = 0.2 * corde if r_core is None else r_core
        self.delta = delta
        self.gamma0 = 0.0
        self.g_harm = tuple(gamma_harm)
        v = mu * omega * R
        self.v_air = np.array([-v * np.cos(alpha_disque), 0.0, -v * np.sin(alpha_disque)])
        self.x = None
        self.residus = []
        self.succes = False
        self._vi_moy = 0.0

    # ── géométrie ─────────────────────────────────────────────────────────
    def circulation(self, poussee):
        self.gamma0 = 2.0 * poussee / (RHO * self.n_b * self.omega * self.R ** 2)
        return self.gamma0

    def _psi_r(self):
        return np.arange(self.n_psi) * self.dz

    def _gamma(self, psi):
        return self.gamma0 * (1.0 + self.g_harm[0] * np.cos(psi) + self.g_harm[1] * np.sin(psi))

    def _bout(self, psi):
        return np.column_stack([-self.R * np.cos(psi), -self.R * np.sin(psi), np.zeros_like(psi)])

    def rigide(self, v_i):
        """Sillage RIGIDE : chaque élément convecté en ligne droite par
        V_air − v_i ẑ (hélice inclinée, non déformée)."""
        psi = self._psi_r()
        z = np.arange(self.n_z + 1) * self.dz
        conv = self.v_air + np.array([0.0, 0.0, -v_i])
        return self._bout(psi)[:, None, :] + conv[None, None, :] * z[None, :, None] / self.omega

    def _idx(self, k, j, b=0):
        """indice de lâcher du nœud d'âge j du filament de la pale b, pale 0 à ψ_b = k·Δζ"""
        return (k - j - b * self.n_psi // self.n_b) % self.n_psi

    def _filaments(self, x, k):
        """Les N_b filaments (nœuds libres + lointain prescrit) et les tourbillons
        liés quand la pale 0 est à l'azimut k·Δζ. Rend (A, B, Γ, r_c)."""
        j = np.arange(self.n_z + 1)
        conv = self.v_air + np.array([0.0, 0.0, -self._vi_moy])
        seg_a, seg_b, gam, rc = [], [], [], []
        for b in range(self.n_b):
            i = self._idx(k, j, b)
            fil = x[i, j]
            # lointain : prolongement en ligne droite au taux de convection moyen
            m = np.arange(1, self.n_loin + 1)
            far = fil[-1][None, :] + conv[None, :] * m[:, None] * self.dz / self.omega
            full = np.vstack([fil, far])
            seg_a.append(full[:-1]); seg_b.append(full[1:])
            # Γ de chaque segment = Γ au lâcher de son nœud amont
            psi_r = self._idx(k, np.arange(full.shape[0] - 1), b) * self.dz
            gam.append(self._gamma(psi_r))
            zs = (np.arange(full.shape[0] - 1) + 0.5) * self.dz
            rc.append(np.sqrt(self.r_core ** 2 + 4 * 1.25643 * self.delta * 1.5e-5 * zs / self.omega))
            # tourbillon lié, pied → bout, de la pale b à son azimut
            pb = (k - b * self.n_psi // self.n_b) % self.n_psi * self.dz
            e = np.array([-np.cos(pb), -np.sin(pb), 0.0])
            seg_a.append((self.r_pied * e)[None, :]); seg_b.append((self.R * e)[None, :])
            gam.append(np.array([self._gamma(pb)])); rc.append(np.array([self.r_core]))
        return (np.vstack(seg_a), np.vstack(seg_b), np.concatenate(gam), np.concatenate(rc))

    def vitesse(self, pts, k, x=None):
        """Vitesse induite en `pts` (repère moyeu) quand la pale 0 est à k·Δζ."""
        x = self.x if x is None else x
        a, b, g, rc = self._filaments(x, k)
        return biot_savart(np.atleast_2d(pts), a, b, g, rc)

    def _vitesses(self, x):
        """v(ψ_r, ζ) pour toute la grille : pour chaque azimut de pale k, les
        nœuds du filament 0 (âges 0..n_z, lâchers k − j)."""
        v = np.empty_like(x)
        j = np.arange(self.n_z + 1)
        for k in range(self.n_psi):
            i = self._idx(k, j)
            v[i, j] = self.vitesse(x[i, j], k, x)
        return v

    def _marche(self, x, v):
        new = np.empty_like(x)
        new[:, 0] = x[:, 0]
        for jj in range(self.n_z):
            vk = 0.5 * (v[:, jj] + v[:, jj + 1])
            new[:, jj + 1] = x[:, jj] + (self.v_air[None, :] + vk) * self.dz / self.omega
        return new

    def glauert(self):
        """λ de Glauert à la circulation courante : λ = C_T/(2√(μ² + λ²))."""
        ct = self.gamma0 * self.n_b / (2 * np.pi * self.omega * self.R ** 2) if self.gamma0 > 0 else 0.0
        lam = np.sqrt(ct / 2.0) if ct > 0 else 0.0
        for _ in range(60):
            lam = ct / (2.0 * np.sqrt(self.mu ** 2 + lam ** 2)) if ct > 0 else 0.0
        return lam, ct

    def resout(self, v_i0=None, tol=1e-4, iters=500, relax=0.5, libre=True):
        """Relaxation : x ← x + ω(marche(x) − x), résidu = max|marche(x) − x|/R.
        `libre=False` garde le sillage rigide (contrôle négatif).

        DOMAINE MESURÉ : converge à μ ≥ 0,15 — à 2 tours libres jusqu'à 1e-6
        (μ 0,15 : 312 it ; μ 0,30 : 194 it), à 3 tours le troisième freine
        (μ 0,15 : 6e-5 après 500 it, soit 10 µm sur R = 175 mm ; μ 0,30 :
        2e-6 en 328 it). `tol` = 1e-4 R est le plancher UTILE : λ₀ et k_x n'y
        bougent plus (Δψ/2 les déplace de 0,4 %). En dessous de μ ≈ 0,12 elle
        retrouve l'instabilité d'appariement
        du stationnaire : μ 0,10 stagne à 1e-4, μ 0,08 oscille à 3e-2, μ 0,05
        plafonne à 8e-5 (ω 0,25, 800 it) — et ω 0,15 ne fait pas mieux. Le
        régime μ < 0,12 est la zone où ni cette relaxation ni le Newton du
        stationnaire (μ = 0) ne s'appliquent tels quels ; déclaré."""
        if v_i0 is None:
            v_i0 = self.glauert()[0] * self.omega * self.R
        self._vi_moy = v_i0
        if self.x is None or self.x.shape[:2] != (self.n_psi, self.n_z + 1):
            self.x = self.rigide(v_i0)
        self.residus, self.succes = [], False
        if not libre or self.gamma0 == 0.0:
            self.residus, self.succes = [0.0], True
            return self.x
        x = self.x
        for _ in range(iters):
            v = self._vitesses(x)
            new = self._marche(x, v)
            r = float(np.abs(new - x).max() / self.R)
            self.residus.append(r)
            x = x + relax * (new - x)
            # la convection moyenne du lointain suit le sillage
            self._vi_moy = max(-float(np.mean(v[:, :, 2])), 0.0)
            if r < tol:
                self.succes = True
                break
        self.x = x
        return x

    # ── sorties ───────────────────────────────────────────────────────────
    def inflow_disque(self, n_r=10, n_az=24):
        """λ(r̄, ψ) sur le disque (vers le bas > 0, adimensionné par ΩR), moyenné
        sur toutes les phases du sillage (les N_ψ azimuts de pale). Rend
        (r̄, ψ, λ[n_r, n_az])."""
        r = (np.arange(n_r) + 0.5) / n_r
        az = (np.arange(n_az) + 0.5) / n_az * 2 * np.pi
        rr, aa = np.meshgrid(r, az, indexing="ij")
        pts = np.column_stack([(-rr * np.cos(aa)).ravel() * self.R,
                               (-rr * np.sin(aa)).ravel() * self.R, np.zeros(rr.size)])
        vz = np.zeros(pts.shape[0])
        for k in range(self.n_psi):
            vz += self.vitesse(pts, k)[:, 2]
        lam = -vz / self.n_psi / (self.omega * self.R)
        return r, az, lam.reshape(n_r, n_az)

    def harmoniques(self, n_r=10, n_az=24, rmax=1.0):
        """(λ₀, λ1c, λ1s) par moindres carrés de λ = λ₀ + r̄(λ1c cos ψ + λ1s sin ψ)
        sur r̄ ≤ `rmax`, pondérés par l'aire ; k_x = λ1c/λ₀, k_y = λ1s/λ₀.

        DEUX k_x, ET ILS NE DISENT PAS LA MÊME CHOSE (mesuré à μ = 0,15) : sur
        tout le disque le fit rend 1,60, sur r̄ ≤ 0,5 il rend 1,10. Le champ
        n'est pas linéaire — upwash de −0,013 au bord avant, downwash de +0,058
        au bord arrière contre λ₀ = 0,025 — et Coleman comme Drees définissent
        k_x comme le GRADIENT AU CENTRE. C'est `rmax=0.5` qui se compare à eux ;
        `rmax=1.0` est la projection de tout le champ sur le modèle linéaire du
        noyau, celle qu'on lui impose."""
        r, az, lam = self.inflow_disque(n_r, n_az)
        rr, aa = np.meshgrid(r, az, indexing="ij")
        m = (rr <= rmax).ravel()
        w = np.sqrt(rr.ravel())[m]
        A = np.column_stack([np.ones(rr.size), (rr * np.cos(aa)).ravel(), (rr * np.sin(aa)).ravel()])[m]
        c = np.linalg.lstsq(A * w[:, None], lam.ravel()[m] * w, rcond=None)[0]
        return dict(lam0=float(c[0]), lam1c=float(c[1]), lam1s=float(c[2]),
                    kx=float(c[1] / c[0]) if c[0] else float("nan"),
                    ky=float(c[2] / c[0]) if c[0] else float("nan"))

    def geometrie(self):
        """Inclinaison χ du sillage (dérive longitudinale sur descente, premier
        tour) et contraction r_v/R à ζ = 2π (rayon transverse du tube autour de
        son centre, moyenné sur les lâchers)."""
        j2pi = self.n_psi
        d = self.x[:, j2pi] - self.x[:, 0]
        chi = float(np.arctan2(-np.mean(d[:, 0]), -np.mean(d[:, 2])))
        c = np.mean(self.x[:, j2pi], axis=0)
        rv = np.hypot(self.x[:, j2pi, 0] - c[0], self.x[:, j2pi, 1] - c[1])
        return dict(chi=chi, r_2pi=float(np.mean(rv) / self.R))

    def impose(self, N, i):
        """Impose au noyau les harmoniques du sillage (`pose_inflow_harmoniques`)."""
        h = self.harmoniques()
        vt = self.omega * self.R
        N.pose_inflow_harmoniques(i, h["lam0"] * vt, h["lam1c"] * vt, h["lam1s"] * vt,
                                  [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], vt)
        return h


def continuation(R, n_b, omega, corde, poussee, mus, r_pied=0.0, dpsi_deg=15.0, tours_libres=3,
                 tours_loin=4, tol=1e-4, iters=300, relaxs=(0.5, 0.25, 0.1)):
    """CONTINUATION EN μ (6 sept.) : on descend `mus` (décroissants) en repartant
    à chaque pas de la géométrie convergée du précédent, la sous-relaxation
    étant réduite (0,5 → 0,25 → 0,1) jusqu'à convergence. Rend une liste de
    dicts (mu, succes, residu, it, relax, lam0, glauert, kx, coleman, drees).

    MESURÉ (rotor FRELON, 2 tours, 18°) depuis μ 0,15 :
        μ 0,12  ω 0,25   73 it   1e-4   λ₀/Glauert 1,13   k_x 1,03 (Drees 0,99, Coleman 0,83)
        μ 0,10  ω 0,25  184 it   1e-4   λ₀/Glauert 1,10   k_x 0,72 (Drees 0,91, Coleman 0,70)
        μ 0,08  ω 0,1   300 it   3e-2   NON
        μ 0,06  ω 0,1   300 it   2e-2   NON
        μ 0,04  ω 0,1   222 it   1e-4   « convergé » sur un k_x NÉGATIF (−0,40) : géométrie
                                        non physique — un résidu qui descend ne suffit pas
        μ 0,02  ω 0,1   300 it   5e-4   NON
    Et à 3 tours, 15° (le mode complet) :
        μ 0,12  ω 0,1   293 it   1e-4   λ₀/Glauert 1,09   k_x 1,02 (Drees 1,00)
        μ 0,10  ω 0,1   300 it   4e-2   NON — k_x 0,89 (Drees 0,90) : la géométrie
                                        est là, le troisième tour ne se pose pas
        μ 0,08  ω 0,1   300 it   5e-3   NON — k_x 0,67 (Drees 0,77)
    LE PLANCHER DÉPEND DES TOURS LIBRES : 0,10 à 2 tours, 0,12 à 3 — c'est le
    troisième tour, celui qui stagnait déjà à μ 0,15, qui ne converge plus
    quand le sillage s'empile sous le rotor. Sous ce plancher tous les tours
    interagissent et 2–3 tours libres + un lointain rectiligne ne suffisent
    plus : c'est un Newton sur la grille périodique entière (3·N_ψ·N_ζ
    inconnues) qu'il faudrait, comme en stationnaire — chantier nommé, pas
    fait. La marche en pseudo-temps (PC2B) n'a pas été essayée. Un k_x qui
    tombe dans la bande Coleman–Drees sans résidu convergé (μ 0,10 à 3 tours)
    est PUBLIÉ, pas asserté.
    """
    out, prev = [], None
    for mu in mus:
        s = SillageAvancement(R, n_b, omega, corde, mu, r_pied=r_pied, dpsi_deg=dpsi_deg,
                              tours_libres=tours_libres, tours_loin=tours_loin)
        s.circulation(poussee)
        if prev is not None:
            s.x = prev.x.copy()
            s._vi_moy = prev._vi_moy
        relax = relaxs[0]
        for relax in relaxs:
            s.resout(tol=tol, iters=iters, relax=relax)
            if s.succes:
                break
        hm, h5 = s.harmoniques(), s.harmoniques(rmax=0.5)
        _, kx_c, kx_d = coleman(mu, hm["lam0"])
        out.append(dict(mu=mu, succes=s.succes, residu=s.residus[-1], it=len(s.residus), relax=relax,
                        lam0=hm["lam0"], glauert=s.glauert()[0], kx=h5["kx"], coleman=float(kx_c),
                        drees=float(kx_d), sillage=s))
        prev = s
    return out


def coleman(mu, lam):
    """k_x de Coleman–Feingold–Stempin (1945) et de Drees (1949), χ = atan(μ/λ)."""
    chi = np.arctan2(mu, lam)
    kx_c = np.tan(0.5 * chi)
    kx_d = 4.0 / 3.0 * (1.0 - np.cos(chi) - 1.8 * mu ** 2) / np.sin(chi) if mu > 0 else 0.0
    return chi, kx_c, kx_d


def couple_avancement(rotor, mu, cmd, tours=3, h=2e-4, iters=6, tol=0.02, **kw):
    """Couple le sillage d'avancement au rotor du noyau (`vinkulum.rotor.Rotor`) :
    poussée → Γ → sillage → harmoniques → `pose_inflow_harmoniques` → un
    torseur → … Rend (harmoniques, sillage, dernier torseur, itérations)."""
    s = SillageAvancement(rotor.R, rotor.n_pales, rotor.omega, rotor.corde, mu,
                          r_pied=rotor.r0, **kw)
    vt = rotor.omega * rotor.R
    harm, hn, k = None, None, 0
    for k in range(1, iters + 1):
        r = rotor.torseur(cmd, mu=mu, h=h, tours=tours, harm=harm)
        s.circulation(max(r["T"], 0.0))
        s.resout()
        hn = s.harmoniques()
        fin = harm is not None and abs(hn["lam0"] * vt - harm[0]) < tol * max(abs(hn["lam0"] * vt), 1e-9)
        harm = (hn["lam0"] * vt, hn["lam1c"] * vt, hn["lam1s"] * vt)
        if fin:
            break
    return hn, s, r, k, harm


def demo_avancement(rapide=False):
    """Le sillage en avancement jugé par Drees et Coleman, par sa continuité
    avec le stationnaire, et par le noyau.

    MESURÉ (rotor FRELON, 2 pales, C_T 0,0069, Δψ 15°, 3 tours libres) :
    · k_x AU CENTRE (r̄ ≤ 0,5) : 1,08 à μ 0,15 (Drees 1,08, Coleman 0,85) et
      1,12 à μ 0,30 (Drees 1,07, Coleman 0,96). Le sillage RIGIDE (hélice
      inclinée, même Γ) rend 0,92 et 1,06 : la déformation libre ajoute +17 %
      à μ 0,15 et +5 % à 0,30 — c'est peu, et c'est attendu : Coleman EST un
      tube incliné non déformé, et le gradient longitudinal est d'abord une
      affaire d'inclinaison ;
    · λ₀ contre Glauert : +4 % (0,15), −2 % (0,30) ;
    · Γ à 1/rev (−30 % en sin ψ) : k_y −0,28 (uniforme +0,02) ET k_x au centre
      −20 % — la 1/rev de Γ fabrique le gradient latéral, et déplace le
      longitudinal : la circulation en azimut n'est pas un détail du modèle ;
    · continuité : à μ 0,05 la relaxation ne converge PAS (résidu 1e-2 à
      3 tours, 8e-5 à 2, déclaré) mais λ₀/Glauert vaut 1,06 contre 1,02 au
      stationnaire (Newton, μ = 0) et 1,04 à μ 0,15 — le rapport est continu
      à ±6 % à travers les deux régimes, sur un état QUASI stationnaire ;
    · CONTINUATION en μ (6 sept., `continuation`) : depuis 0,15, la géométrie
      convergée du μ précédent fait converger 0,12 (2 et 3 tours) puis 0,10
      à 2 tours seulement (k_x 1,03 et 0,72 : entre Coleman et Drees) ; à
      3 tours 0,10 reste à 4e-2 sur un k_x pourtant à 1 % de Drees ; 0,08 et
      0,06 ne convergent pas, et 0,04 « converge » sur un k_x négatif. Le
      plancher est μ 0,12 (3 tours) ou 0,10 (2 tours), mesuré, publié ;
    · Δψ ÷ 2 : λ₀ +0,4 %, k_x −0,3 % ;
    · couplage noyau (rotor 2 pales, θ₀ 8°, μ 0,15) : λ₀ sillage à +5,5 % de
      Glauert en 3 allers-retours ; trim à T, Mx = My = 0 : le sillage
      déplace le cyclique LATÉRAL de +2,46° (θ1c 0,83° → 3,29°) et ne touche
      pas au longitudinal (Δθ1s +0,01°) — le gradient longitudinal d'inflow
      se paie en latéral, comme le dit toute la mécanique du vol.
    Chiffres à 3 tours libres, Δψ 15° ; en rapide (2 tours, 18°) : k_x 1,11
    et 1,07, λ₀ +9 % et −3 %.
    Tolérances : k_x au centre contre Drees à 15 % (Coleman et Drees diffèrent
    entre eux de 10–25 % sur 0,1 < μ < 0,3, et les sillages libres publiés
    tombent entre les deux — Leishman, *Principles*, ch. 10) ; λ₀ contre
    Glauert à 15 % ; Δψ ÷ 2 à 5 %.
    """
    from vinkulum.rotor import Rotor
    t0 = time.time()
    R, n_b, omega, corde = 0.175, 2, 334.9, 0.020          # FRELON
    r_pied = 0.15 * R
    poussee = 2.8
    vtip = omega * R
    c_t = poussee / (RHO * np.pi * R * R * vtip * vtip)
    dpsi, tl = (18.0, 2) if rapide else (15.0, 3)
    print("╔═ vinkulum — SILLAGE LIBRE en AVANCEMENT : grille périodique, relaxation")
    print(f"║ rotor Ø{2 * R * 1e3:.0f} · {n_b} pales · C_T {c_t:.5f} · Δψ {dpsi:.0f}° · {tl} tours libres + 4 prescrits")
    out = {}

    def monte(mu, **kw):
        s = SillageAvancement(R, n_b, omega, corde, mu, r_pied=r_pied, dpsi_deg=dpsi,
                              tours_libres=tl, tours_loin=4, **kw)
        s.circulation(poussee)
        return s

    for mu in (0.15, 0.30):
        s = monte(mu)
        s.resout()
        hm, hc5 = s.harmoniques(), s.harmoniques(rmax=0.5)
        g = s.geometrie()
        chi, kx_c, kx_d = coleman(mu, hm["lam0"])
        lam_gl = s.glauert()[0]
        sr = monte(mu)
        sr.resout(libre=False)
        hr5 = sr.harmoniques(rmax=0.5)
        out[mu] = dict(h=hm, kx=hc5["kx"], kx_rigide=hr5["kx"], g=g, coleman=float(kx_c), drees=float(kx_d),
                       glauert=float(lam_gl), it=len(s.residus), res=s.residus[-1], succes=s.succes)
        print(f"║ μ {mu:.2f} : relaxation {len(s.residus)} it, résidu {s.residus[0]:.1e} → {s.residus[-1]:.1e} · "
              f"λ₀ {hm['lam0']:.5f} (Glauert {lam_gl:.5f}, {100 * (hm['lam0'] / lam_gl - 1):+.1f} %) · "
              f"χ {np.degrees(g['chi']):.1f}° (atan(μ/λ) {np.degrees(chi):.1f}°) · r_v(2π) {g['r_2pi']:.3f}")
        print(f"║        k_x au centre {hc5['kx']:+.3f} contre Drees {kx_d:+.3f} ({100 * (hc5['kx'] / kx_d - 1):+.1f} %) "
              f"et Coleman {kx_c:+.3f} ({100 * (hc5['kx'] / kx_c - 1):+.1f} %) · rigide {hr5['kx']:+.3f} · "
              f"sur tout le disque {hm['kx']:+.3f} · k_y {hc5['ky']:+.3f}")
    # Γ à 1/rev : ce que ça change
    s1 = monte(0.15, gamma_harm=(0.0, -0.3))
    s1.resout()
    h1 = s1.harmoniques(rmax=0.5)
    print(f"║ μ 0,15 avec Γ(ψ) = Γ₀(1 − 0,3 sin ψ) : k_x {h1['kx']:+.3f} (uniforme {out[0.15]['kx']:+.3f}), "
          f"k_y {h1['ky']:+.3f} (uniforme {out[0.15]['h']['ky']:+.3f}) — la 1/rev de Γ fait le latéral et déplace le longitudinal")
    # continuité avec le stationnaire, par le rapport à Glauert
    mu_c = 0.05
    sc = monte(mu_c)
    sc.resout(relax=0.25, iters=800)              # 36 s : le prix de la continuité
    hc = sc.harmoniques()
    q_c = hc["lam0"] / sc.glauert()[0]
    s0 = Sillage(R, n_b, omega, corde, r_pied=r_pied, dzeta_deg=dpsi, tours_libres=tl, tours_loin=10 - tl)
    s0.circulation(poussee)
    s0.resout()
    q_0 = s0.inflow_disque() / vtip / np.sqrt(c_t / 2.0)
    q_15 = out[0.15]["h"]["lam0"] / out[0.15]["glauert"]
    print(f"║ continuité λ₀/Glauert : μ 0 (Newton) {q_0:.3f} · μ {mu_c} {q_c:.3f} (relaxation {'convergée' if sc.succes else 'NON convergée'}, "
          f"résidu {sc.residus[-1]:.1e}) · μ 0,15 {q_15:.3f}")
    # CONTINUATION EN μ sous 0,12 (cf. `continuation`) : jusqu'où ça descend
    mus = (0.15, 0.12) if rapide else (0.15, 0.12, 0.10, 0.08)
    cont = continuation(R, n_b, omega, corde, poussee, mus, r_pied=r_pied, dpsi_deg=dpsi,
                        tours_libres=tl, tours_loin=4)
    for c in cont:
        print(f"║ continuation μ {c['mu']:.2f} : {'convergé' if c['succes'] else 'NON convergé'} (ω {c['relax']}, "
              f"{c['it']} it, résidu {c['residu']:.1e}) · λ₀/Glauert {c['lam0'] / c['glauert']:.3f} · "
              f"k_x {c['kx']:+.3f} (Drees {c['drees']:+.3f}, Coleman {c['coleman']:+.3f})")
    # contrôles négatifs
    sz = monte(0.15)
    sz.gamma0 = 0.0
    sz.resout(v_i0=out[0.15]["h"]["lam0"] * vtip)
    hz = sz.harmoniques()
    sr = monte(0.15)
    sr.resout(libre=False)
    hr = sr.harmoniques()
    gr = sr.geometrie()
    print(f"║ contrôle négatif Γ = 0 : λ₀ {hz['lam0']:.1e} — rien")
    print(f"║ contrôle négatif RIGIDE (même Γ) : λ₀ {hr['lam0']:.5f} ({100 * (hr['lam0'] / out[0.15]['h']['lam0'] - 1):+.1f} % du libre), "
          f"k_y {hr['ky']:+.3f}, r_v(2π) {gr['r_2pi']:.3f} (pas de contraction)")
    e_dz = None
    if not rapide:
        s2 = SillageAvancement(R, n_b, omega, corde, 0.15, r_pied=r_pied, dpsi_deg=dpsi / 2,
                               tours_libres=tl, tours_loin=4)
        s2.circulation(poussee)
        s2.resout()
        h2, h25 = s2.harmoniques(), s2.harmoniques(rmax=0.5)
        e_dz = (h2["lam0"] / out[0.15]["h"]["lam0"] - 1.0, h25["kx"] / out[0.15]["kx"] - 1.0)
        print(f"║ Δψ/2 à μ 0,15 : λ₀ {100 * e_dz[0]:+.2f} %, k_x {100 * e_dz[1]:+.2f} %")
    # couplage noyau
    rot = Rotor(n_pales=2, R=R, r0=r_pied, corde=corde, omega=omega, vrillage=0.0,
                rho=RHO, lock=6.0, segments=4, gauss=3)
    mu_k = 0.15
    cmd0 = (8.0, 0.0, 0.0)
    hn, sk, rk, k_it, harm = couple_avancement(rot, mu_k, cmd0, tours=2 if rapide else 3, dpsi_deg=dpsi,
                                                tours_libres=tl, tours_loin=4)
    ct_k = rk["T"] / (RHO * np.pi * R * R * vtip * vtip)
    lam_u = ct_k / (2.0 * np.sqrt(mu_k ** 2 + hn["lam0"] ** 2))
    print(f"║ couplage noyau (μ {mu_k}, θ₀ {cmd0[0]}°) : T {rk['T']:.3f} N, λ₀ sillage {hn['lam0']:.5f} contre Glauert "
          f"{lam_u:.5f} ({100 * (hn['lam0'] / lam_u - 1):+.1f} %), imposé au noyau λ1c/λ₀ {hn['kx']:+.3f}, {k_it} allers-retours")
    trims = None
    if not rapide:
        cible = dict(T=rk["T"], Mx=0.0, My=0.0)
        xu, _, _ = rot.trim(cible, x0=cmd0, mu=mu_k, tours=3)
        xs, _, _ = rot.trim(cible, x0=tuple(xu), mu=mu_k, tours=3, harm=harm)
        trims = ([float(v) for v in xu], [float(v) for v in xs])
        print(f"║ trim (T {cible['T']:.3f} N, Mx = My = 0) : inflow uniforme θ₀ {xu[0]:.2f}° θ1c {xu[1]:.2f}° θ1s {xu[2]:.2f}° · "
              f"avec sillage θ₀ {xs[0]:.2f}° θ1c {xs[1]:.2f}° θ1s {xs[2]:.2f}° — Δθ1c {xs[1] - xu[1]:+.2f}°, Δθ1s {xs[2] - xu[2]:+.2f}°")
    print(f"╚═ sillage en avancement — {time.time() - t0:.1f} s")

    for mu in (0.15, 0.30):
        o = out[mu]
        assert o["succes"], ("la relaxation ne converge pas", mu, o["res"])
        # k_x AU CENTRE contre Drees : Coleman et Drees diffèrent entre eux de
        # 10–25 % sur 0,1 < μ < 0,3, les sillages libres publiés sont entre les deux
        assert abs(o["kx"] / o["drees"] - 1.0) < 0.15, ("k_x n'est pas celui de Drees", mu, o["kx"], o["drees"])
        assert o["kx"] > o["coleman"], ("k_x devrait être au-dessus de Coleman (tube non déformé)", mu, o["kx"], o["coleman"])
        assert abs(o["h"]["lam0"] / o["glauert"] - 1.0) < 0.15, ("λ₀ s'écarte de Glauert", mu, o["h"])
    # continuité : le rapport à Glauert reste dans ±15 % des deux côtés du régime
    # (le résidu à μ 0,05 est PUBLIÉ, pas asserté : la relaxation n'y converge pas)
    assert abs(q_c - 1.0) < 0.15 and abs(q_0 - 1.0) < 0.15, (q_0, q_c, sc.residus[-1])
    assert hz["lam0"] == 0.0, ("Γ = 0 induit quelque chose", hz)
    assert abs(hr["ky"]) < 0.05 and abs(gr["r_2pi"] - 1.0) < 1e-6, ("le rigide ne doit ni contracter ni pencher latéralement", hr, gr)
    assert abs(hr["lam0"] / out[0.15]["h"]["lam0"] - 1.0) > 0.02, ("le libre devrait différer du rigide", hr)
    assert abs(h1["ky"]) > 5 * abs(out[0.15]["h"]["ky"]) + 0.05, ("la 1/rev de Γ devrait faire un k_y", h1)
    if e_dz is not None:
        assert max(abs(e) for e in e_dz) < 0.05, ("dépend du pas d'azimut", e_dz)
    assert abs(hn["lam0"] / lam_u - 1.0) < 0.20, ("couplage noyau hors de Glauert", hn, lam_u)
    # · la continuation tient jusqu'à μ 0,12 (2 ET 3 tours) : convergée, k_x ENTRE
    #   les deux modèles publiés (Coleman et Drees diffèrent de 10–25 % ; ±10 % de
    #   la bande), λ₀ à Glauert. Au-dessous : publié, le plancher est mesuré
    for c in cont:
        lo, hi = min(c["coleman"], c["drees"]), max(c["coleman"], c["drees"])
        if c["mu"] >= 0.12:
            assert c["succes"], ("la continuation devrait converger jusqu'à μ 0,12", c["mu"], c["residu"])
        if c["succes"]:
            assert 0.9 * lo <= c["kx"] <= 1.1 * hi, ("k_x hors de la bande Coleman–Drees", c["mu"], c["kx"], lo, hi)
            assert abs(c["lam0"] / c["glauert"] - 1.0) < 0.15, ("λ₀ hors de Glauert", c["mu"], c["lam0"])
    return dict(out={str(k): v for k, v in out.items()}, gamma_1rev=h1, continuite=(q_0, q_c, q_15),
                rigide=hr, dpsi=e_dz, couplage=hn, trims=trims,
                continuation=[{k: v for k, v in c.items() if k != "sillage"} for c in cont])

# ── banc ──────────────────────────────────────────────────────────────────────
def demo(rapide=False):
    """Le sillage libre jugé par la théorie du disque et par Landgrebe.

    CE QUE LA MESURE A DIT, ET QUI N'ÉTAIT PAS ÉCRIT D'AVANCE : k₁, la pente
    axiale du tourbillon marginal AVANT le passage de la pale suivante, sort
    avec le MAUVAIS SIGNE (+0,003 contre −0,024 chez Landgrebe). Ce n'est pas
    un bug, c'est le modèle : sans nappe interne, le tourbillon jeune n'est
    porté que par celui de la pale précédente, presque à sa hauteur et juste
    à l'intérieur — il le SOULÈVE. La nappe interne, de vorticité opposée, est
    ce qui le fait descendre dans la réalité. Le chiffre est publié, il n'est
    PAS asserté ; ce qui est asserté est k₂ (après le passage), de bon signe
    et à ~50 % de Landgrebe. Et il dépend du noyau : δ = 1000 rend
    k₁ = −0,012 — un bouton, pas une prédiction.

    Tolérances, et d'où elles viennent :
    · inflow moyen contre √(C_T/2) : 15 % — un sillage libre en stationnaire
      retombe sur la théorie du disque à ~10 % (Leishman, *Principles of
      Helicopter Aerodynamics*, ch. 10 ; le lointain prescrit et la
      circulation uniforme en coûtent quelques points) ;
    · contraction à ζ = 2π contre Landgrebe : 10 % — c'est la dispersion de
      la corrélation elle-même sur les rotors qui l'ont produite ;
    · k₂ (convection axiale après le premier passage de pale) : signe, et
      ordre de grandeur à un facteur 2 — Landgrebe le tire de photos de fumée
      sur des rotors à 2–8 pales, et le nôtre n'a ni vrillage ni nappe
      interne. C'est le chiffre le plus faible du banc, et il est dit.
    """
    R, n_b, omega, corde = 0.175, 2, 334.9, 0.020          # FRELON : Ø350, 3 197 tr/min
    r_pied = 0.15 * R
    poussee = 2.8
    a, vtip = np.pi * R * R, omega * R
    c_t = poussee / (RHO * a * vtip * vtip)
    sigma = n_b * corde / (np.pi * R)
    lam_mom = np.sqrt(c_t / 2.0)
    dz, tl = (18.0, 2) if rapide else (15.0, 3)
    t0 = time.time()
    print("╔═ vinkulum — SILLAGE LIBRE en stationnaire : tourbillons marginaux, point fixe par Newton")
    print(f"║ rotor Ø{2 * R * 1e3:.0f} · {n_b} pales · C_T {c_t:.5f} · σ {sigma:.4f} · "
          f"Δζ {dz:.0f}° · {tl} tours libres + {10 - tl} prescrits · noyau r_c {0.2 * corde * 1e3:.1f} mm, δ 100")

    def monte(dzeta, tours, **kw):
        # 6 éléments de nappe, en rapide aussi : à 4 la contraction sort à +8 %
        # de celle sans nappe (mesuré), et le gain n'est que 10 s
        s = Sillage(R, n_b, omega, corde, r_pied=r_pied, dzeta_deg=dzeta, tours_libres=tours,
                    tours_loin=10 - tours, **kw)
        s.circulation(poussee)
        return s

    def mesure(s):
        lam = s.inflow_disque() / vtip
        z, rv, zv = s.marginal()
        i2pi = int(round(2 * np.pi / s.dz))
        kb = int(round(2 * np.pi / n_b / s.dz))
        k1 = np.polyfit(z[: kb + 1], zv[: kb + 1], 1)[0]
        k2 = np.polyfit(z[kb: i2pi + 1], zv[kb: i2pi + 1], 1)[0]
        return lam, z, rv, zv, i2pi, kb, k1, k2

    s = monte(dz, tl)
    s.resout()
    lam, z, rv, zv, i2pi, kb, k1_m, k2_m = mesure(s)
    r_l, z_l, k1_l, k2_l = landgrebe(z, c_t, sigma, n_b)
    e_lam = lam / lam_mom - 1.0
    e_r = rv[i2pi] / r_l[i2pi] - 1.0
    e_k2 = k2_m / k2_l - 1.0
    print(f"║ Newton : {len(s.residus) - 1} itérations, résidu {s.residus[0]:.1e} → {s.residus[-1]:.1e} "
          f"(déplacement / R) — " + " → ".join(f"{r:.1e}" for r in s.residus))
    print(f"║ nappe interne : {s.n_nappe} éléments, θ {np.degrees(s.theta):.2f}°, Γ marginal {s.gamma:.4f} "
          f"(uniforme équivalent {2 * poussee / (RHO * n_b * omega * R * R):.4f}) m²/s")
    print(f"║ inflow moyen λ {lam:.5f} contre √(C_T/2) = {lam_mom:.5f} — écart {100 * e_lam:+.1f} %")
    print(f"║ contraction r_v/R à ζ = 2π : {rv[i2pi]:.4f} contre Landgrebe {r_l[i2pi]:.4f} — "
          f"écart {100 * e_r:+.1f} % · asymptote {rv[-1]:.3f} (Landgrebe 0,78)")
    print(f"║ convection axiale k₂ (π < ζ < 2π) : {k2_m:+.4f} contre Landgrebe {k2_l:+.4f} — "
          f"écart {100 * e_k2:+.1f} %")
    print(f"║ k₁ (ζ < π) : {k1_m:+.4f} contre Landgrebe {k1_l:+.4f} — bon signe, facteur {k1_l / k1_m:.2f}")
    # ── contrôle négatif : SANS nappe, le modèle d'avant au chiffre près
    sn = monte(dz, tl, nappe=False)
    sn.resout()
    lam_n, _, rv_n, _, _, _, k1_n, k2_n = mesure(sn)
    print(f"║ contrôle négatif SANS nappe (Γ uniforme, le modèle d'avant le 6 sept.) : k₁ {k1_n:+.4f} "
          f"(MAUVAIS signe : le tourbillon jeune est soulevé par le précédent), r_v(2π) {rv_n[i2pi]:.4f}, "
          f"λ {lam_n / lam_mom:.3f}×momentum — la nappe déplace la contraction de {100 * (rv[i2pi] / rv_n[i2pi] - 1):+.1f} %")

    # ── contrôles négatifs : Γ = 0 ; hélice rigide ────────────────────────
    s0 = monte(dz, tl)
    s0.circulation(0.0)                  # Γ = 0 partout : marginal ET nappe
    s0.resout(v_i0=lam_mom * vtip)
    _, rv0, _ = s0.marginal()
    lam0 = s0.inflow_disque() / vtip
    sr = monte(dz, tl)
    sr.resout(libre=False)
    _, rvr, _ = sr.marginal()
    lam_r = sr.inflow_disque() / vtip
    print(f"║ contrôle négatif Γ = 0 : inflow {lam0:.1e}, contraction {abs(rv0 - 1).max():.1e} — rien")
    print(f"║ contrôle négatif hélice RIGIDE (même Γ) : r_v/R {rvr[i2pi]:.4f} (pas de contraction), "
          f"λ {lam_r:.5f} ({100 * (lam_r / lam_mom - 1):+.1f} % du momentum)")

    # ── sensibilité au pas d'âge et au nombre de tours libres ──────────────
    e_dz = None
    if not rapide:
        s2 = monte(dz / 2, 2)
        s2.resout()
        sb = monte(dz, 2)
        sb.resout()
        _, rv2, _ = s2.marginal()
        _, rvb, _ = sb.marginal()
        e_dz = (s2.inflow_disque() / sb.inflow_disque() - 1.0, rv2[2 * i2pi] / rvb[i2pi] - 1.0)
        e_tl = (sb.inflow_disque() / vtip / lam - 1.0, rvb[i2pi] / rv[i2pi] - 1.0)
        print(f"║ Δζ/2 : λ bouge de {100 * e_dz[0]:+.2f} %, r_v(2π) de {100 * e_dz[1]:+.2f} % · "
              f"2 tours libres au lieu de {tl} : {100 * e_tl[0]:+.2f} % et {100 * e_tl[1]:+.2f} %")
    else:
        print("║ (mode rapide : sans le contrôle de pas d'âge — Δζ/2 coûte 60 s)")

    # ── couplage au noyau ─────────────────────────────────────────────────
    N, i = rotor_noyau(R, n_b, omega, corde, 7.0, r_pied)
    v_s, v_d, t_n, k_c, _ = couple(N, i, R, n_b, omega, corde, r_pied=r_pied, dzeta_deg=dz,
                                   tours_libres=tl, tours_loin=10 - tl)
    e_c = v_s / v_d - 1.0
    print(f"║ couplage noyau : poussée {t_n:.3f} N, v_i sillage {v_s:.3f} m/s contre théorie du "
          f"disque {v_d:.3f} — écart {100 * e_c:+.1f} %, {k_c} allers-retours")
    print(f"╚═ sillage libre — {time.time() - t0:.1f} s")

    assert s.succes and s.residus[-1] < 1e-9, ("Newton ne converge pas", s.residus)
    assert abs(e_lam) < 0.15, ("l'inflow du sillage libre s'écarte de la théorie du disque", lam, lam_mom)
    assert abs(e_r) < 0.10, ("la contraction n'est pas celle de Landgrebe", rv[i2pi], r_l[i2pi])
    assert k2_m < 0 and 0.5 < k2_m / k2_l < 2.0, ("la convection axiale après passage", k2_m, k2_l)
    # · LA NAPPE REND À k₁ SON SIGNE (Landgrebe) et son ordre de grandeur ; sans
    #   elle il est positif — c'est le contrôle négatif, et la raison du chantier
    assert k1_m < 0 and 0.5 < k1_m / k1_l < 2.0, ("k₁ avec nappe", k1_m, k1_l)
    assert k1_n > 0, ("sans nappe k₁ devrait rester du mauvais signe (contrôle négatif)", k1_n)
    assert abs(rv[i2pi] / rv_n[i2pi] - 1.0) < 0.06, ("la nappe déplace trop la contraction", rv[i2pi], rv_n[i2pi])
    assert lam0 == 0.0 and abs(rv0 - 1).max() < 1e-12, ("Γ = 0 doit ne rien induire (inflow exactement nul, géométrie à l'arrondi près)", lam0)
    assert abs(rvr[i2pi] - 1.0) < 1e-12, "une hélice rigide ne se contracte pas"
    assert abs(rv[i2pi] - 1.0) > 0.1, ("le sillage LIBRE doit se contracter", rv[i2pi])
    if e_dz is not None:
        assert max(abs(e) for e in e_dz) < 0.05, ("le résultat dépend du pas d'âge", e_dz)
        assert max(abs(e) for e in e_tl) < 0.02, ("le résultat dépend du nombre de tours libres", e_tl)
    assert abs(e_c) < 0.15, ("le couplage noyau s'écarte de la théorie du disque", v_s, v_d)
    return dict(lam=lam, lam_momentum=lam_mom, r_2pi=float(rv[i2pi]), landgrebe=float(r_l[i2pi]),
                k1=float(k1_m), k1_landgrebe=float(k1_l), k2=float(k2_m), k2_landgrebe=float(k2_l),
                newton=s.residus, dzeta=e_dz, couplage=(v_s, v_d, t_n, k_c),
                avancement=demo_avancement(rapide))


if __name__ == "__main__":
    demo(rapide="rapide" in sys.argv)
