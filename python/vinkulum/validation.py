"""vinkulum.validation — la CAMPAGNE contre des essais et des modèles publiés.

    python -m vinkulum.validation            # campagne complète → docs/VALIDATION.md
    python -m vinkulum.validation rapide

Un banc juge le solveur contre l'analytique ; une VALIDATION le juge contre
ce que d'autres ont mesuré ou calculé et publié. La différence est celle
d'un code de recherche à un code de métier, et elle se mesure ici cas par
cas, avec la source, le chiffre publié, le chiffre du noyau, l'écart, et ce
que la comparaison ne couvre pas. Rien n'est recoté pour passer ; un écart
expliqué se publie avec sa cause.

CE QUE LA CAMPAGNE JUGE (référence → grandeur) :
  A. Wright, Smith, Thresher & Wang 1982 (= Hodges & Rutkowski 1981) —
     éventail d'une poutre encastrée tournante, ω̄₁(Ω̄) jusqu'à Ω̄ = 12 ;
  B. Johnson, *Helicopter Theory* (1980) §5-3 / Bramwell — coefficients de
     battement β₀, β1c, β1s d'un rotor articulé (ν = 1) en avancement,
     inflow uniforme, forme fermée ;
  C. Benedict, Winslow, Hayden & Chopra, *J. Aircraft* 52 (2015) — figure de
     mérite MESURÉE de micro-rotors (R 41 mm, 10 000 tr/min), quatre sections,
     polaires NeuralFoil aux Reynolds de l'essai. RÉSULTAT (5 sept.) : la
     plaque mince passe à −1 % ; les trois sections épaisses sortent de +12 à
     +35 % avec l'inflow UNIFORME du noyau, et retombent vers la mesure avec
     le ki effectif 1,34 que la gate V-2 de FRELON a calibré sur ces rotors —
     c'est le modèle d'inflow en stationnaire (pas de perte de bout, pas de
     distribution) qui est jugé, et il est jugé insuffisant : sillage libre
     (`vinkulum.sillage`) ou perte de bout à coupler, chantier nommé ;
  D. Lock généralisé au mode (Johnson ch. 9) — amortissement aéro du premier
     mode de battement d'une pale élastique ;
  E. Dowell & Traybar 1975 (poutre de Princeton), Multibody System Dynamics
     37 (2016) 29–48 (flambement latéral, arbre tournant, quatre-barres
     flexible), Bennett 1903, Kapitza 1951, Andrews (Schiehlen 1990),
     srscm / 6barmech (decks MBDyn) — déjà jugés dans `verification`,
     rejoués ici pour figurer au même tableau ;
  F. Drees 1949 et Coleman–Feingold–Stempin 1945 — gradient d'inflow
     longitudinal k_x du sillage libre en avancement, au centre du disque ;
  G. Elliott, Althoff & Sailey 1988 (NASA TM 100541/100543) — inflow MESURÉ
     au vélocimètre laser sur un rotor 4 pales à μ 0,15 et 0,30, 121 et
     146 points lus par OCR (`vinkulum.elliott1988`), contre Glauert, Drees
     et le sillage libre aux mêmes points ;
  H. Ramsay, Hoffman & Gregorek 1995 (OSU/NREL TP-442-7817) — S809 en
     tangage oscillant MESURÉ (table B3 transcrite, repères de la fig. C25),
     contre Leishman–Beddoes et Øye.

CE QUE LA CAMPAGNE NE COUVRE PAS, dit tout de suite : McCroskey 1982
(NASA TM 84245) est un scan de courbes sans table — le S809 d'OSU le
remplace, avec des repères lus à ±0,05 ; aucun rotor grandeur nature en
avancement mesuré avec charges de pale (HART II, UH-60 Airloads) — les
données ne sont pas publiques hors ligne ; les constantes de LB ne sont
calibrées sur aucune des sections jugées.

BAS REYNOLDS (consigne du 6 sept.) : FRELON vole à Re 5e4–8e4, et à cette
décade la polaire dépend du N_crit (bulle laminaire), le décrochage part du
bord d'attaque, et les constantes de Leishman–Beddoes — mesurées sur NACA
0012 à Re ≥ 1e6 — n'ont aucune raison de tenir. Dans cette campagne, SEUL
le cas C (Maryland, Re 3e4, mesuré) est à bas Reynolds ; Elliott (Re ~1e6)
et le S809 (Re 1e6) n'y sont pas, et le tableau le marque ligne par ligne.
Un essai de DÉCROCHAGE DYNAMIQUE à Re < 1e5 n'a pas été trouvé sous forme
de données (les études existantes — Kim & Chang, LES à 6,6e4, plaques en
tangage — publient des courbes, pas des tables) : LB reste hors domaine
mesuré à bas Re, et c'est écrit ici plutôt que caché derrière le S809.

LE CORPUS MARS DE LA NASA, EXAMINÉ ET REFUSÉ (6 sept., question de Paul) :
· Koning, Johnson & Allan 2019 (AIAA J., NTRS 20180006746) et Dull et al.
  2022 (NTRS 20210026235) publient les mesures d'Ingenuity au simulateur
  25 ft du JPL (CO2, ρ 0,0175) — FM contre C_T/σ, figures VECTORIELLES donc
  digitalisables. Refusé : rotor COAXIAL (le noyau n'a pas d'interaction de
  sillage entre rotors) et M_bout ≈ 0,7 (la table de pale du noyau est à UN
  Mach — pas de correction de compressibilité). Deux absences, pas une.
· Young et al. 2002 (AHS Forum, rotor ISOLÉ Ø2,44 m en chambre Ames) :
  Re_bout 55 000 au point de conception — la décade de FRELON — mais les
  runs mesurés sont à Re 30–37 000 ET M 0,40–0,50, et leur propre table 4
  montre la pente de portance effective tomber de 5,03 à 2,92 /rad avec le
  Mach : c'est la compressibilité qui domine ce qu'ils mesurent, pas le
  Reynolds. Figures rastérisées, hystérésis « puzzling » déclarée par les
  auteurs. Refusé.
· Winslow, Otsuka, Govindarajan & Chopra 2018 (J. Aircraft) : annoncé ici
  comme mesuré, c'est une étude RANS Spalart–Allmaras — même rang que
  NeuralFoil, pas un juge. Rétracté.
· Granlund, Ol & Bernal 2013 (JFM, plaques en tangage, Re 1e4) et Kim &
  Chang 2014 (NACA 0012 en tangage sinusoïdal, Re 2–5e4) : les seules
  mesures de décrochage dynamique à notre décade ; PDF inaccessibles
  (Cambridge, DTIC et ResearchGate rendent une page HTML au téléchargement,
  ScienceDirect payant). Chantier ouvert : obtenir les fichiers.
Ce que l'examen a appris : le corpus Mars est à bas Reynolds ET à haut
Mach — il ne sépare pas les deux, et le noyau ne modélise pas le second.
"""
import sys
import time
from math import cos, degrees, pi, radians, sin, sqrt

import numpy as np

from vinkulum import Noyau, depot_docs

LIGNES = []      # (cas, source, grandeur, référence, noyau, écart, verdict, note, Re)

# LE REYNOLDS DE CHAQUE LIGNE EST ÉCRIT (consigne de Paul, 6 sept. : « attention
# aux cas de bas Reynolds »). FRELON vole à Re 5e4–8e4 ; une validation à Re 1e6
# ne le couvre pas, et une ligne « OK » à haut Re ne dit rien de la bulle
# laminaire, du N_crit ni des constantes de décrochage dynamique à bas Re.
RE_FRELON = (5e4, 8e4)


def _ligne(cas, source, grandeur, ref, val, tol, note="", re=None):
    """`re` : Reynolds de chaque cas — None = sans Reynolds (structure, forme
    fermée), un nombre sinon. « bas Re » = dans la décade de FRELON."""
    e = (val - ref) / abs(ref) if ref else float("nan")
    ok = abs(e) <= tol
    LIGNES.append((cas, source, grandeur, ref, val, e, ok, note, re))
    bas = "" if re is None else (" · BAS Re" if re < 2e5 else " · haut Re")
    print(f"║ {cas:<44} {grandeur:<18} réf {ref:10.4f}  noyau {val:10.4f}  {100 * e:+7.2f} %  "
          f"{'OK' if ok else 'ÉCART'}{bas}")
    return ok


# ── A. éventail de Wright 1982 ───────────────────────────────────────────────
def cas_wright(rapide):
    from vinkulum import pale_elastique as P
    ev, pire = P.eventail(rapide)
    for ob, (w, ref, e) in ev.items():
        _ligne("A · poutre tournante (Wright 1982)", "Wright et al. 1982, table",
               f"ω̄₁ à Ω̄ = {ob:g}", ref, w, 0.01 if rapide else 0.005)


# ── B. battement en avancement, forme fermée ─────────────────────────────────
def cas_johnson(rapide):
    """β₀, β1c, β1s d'un rotor articulé à ν = 1 sans vrillage, inflow uniforme.

    Bramwell (1976) / Johnson (1980, §5-3), pale sans excentricité ni ressort :
        a₀ = γ/8·[θ₀(1 + μ²) − 4λ/3]        β = a₀ − a₁ cos ψ − b₁ sin ψ
        a₁ = (8/3·μ·θ₀ − 2μλ)/(1 − μ²/2)     b₁ = (4/3)·μ·a₀/(1 + μ²/2)
    λ est l'inflow que le NOYAU trouve (Glauert) : le contrôle porte sur la
    cinématique du battement et l'intégration, pas sur le modèle d'inflow.
    Négligé par la forme fermée, donc écart attendu de quelques % : flux
    inversé, découpe de pied (r₀/R = 0,05 ici), plateau de décrochage.
    """
    from vinkulum.rotor import Rotor
    r = Rotor(n_pales=4, R=2.0, r0=0.1, corde=0.121, omega=109.0, vrillage=0.0,
              nu_beta=1.0, zeta=0.0, lock=6.0, segments=8, cl0=0.0)
    th0, mu = 8.0, (0.15 if rapide else 0.20)
    t_end = (5.0 if rapide else 7.0) * r.t_tour
    N = r.monte((th0, 0.0, 0.0), mu=mu, t_end=t_end, montee=1.0)
    h = r.t_tour / 360
    traj = N.simule(t_end, h, tous=2)
    v_i = N.aero()[1][r.inflow_idx][0]
    lam = v_i / r.v_bout
    # β(ψ) de la pale 0 : sin β = composante z de son envergure (colonne 0 de R)
    t = np.array([f[0] for f in traj])
    b = r.idx[0]
    beta = np.array([np.arcsin(np.clip(f[2][b][6], -1, 1)) for f in traj])
    psi = r.omega * t  # pale 0 lâchée à ψ₀ = 0
    k = t > t_end - r.t_tour
    A = np.column_stack([np.ones(k.sum()), np.cos(psi[k]), np.sin(psi[k])])
    b0, b1c, b1s = np.linalg.lstsq(A, beta[k], rcond=None)[0]
    th = radians(th0 - r.alpha0)
    a0 = r.lock / 8 * (th * (1 + mu ** 2) - 4 * lam / 3)
    a1 = (8 / 3 * mu * th - 2 * mu * lam) / (1 - mu ** 2 / 2)
    bb1 = 4 / 3 * mu * a0 / (1 + mu ** 2 / 2)
    src = "Bramwell 1976 / Johnson 1980 §5-3"
    _ligne(f"B · battement articulé μ = {mu:.2f}", src, "β₀ coning (°)", degrees(a0), degrees(b0), 0.06)
    _ligne(f"B · battement articulé μ = {mu:.2f}", src, "β1c (°)", -degrees(a1), degrees(b1c), 0.10)
    _ligne(f"B · battement articulé μ = {mu:.2f}", src, "β1s (°)", -degrees(bb1), degrees(b1s), 0.12,
           note=f"λ noyau {lam:.4f}, γ {r.lock}")


# ── C. micro-rotors Maryland mesurés ─────────────────────────────────────────
def _naca4(m, p, t, n=120):
    x = 0.5 * (1 - np.cos(np.linspace(0, pi, n)))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x - x ** 2)) if p > 0 else 0 * x
    dyc = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x)) if p > 0 else 0 * x
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.vstack([np.column_stack([xu[::-1], yu[::-1]]), np.column_stack([xl[1:], yl[1:]])])


def _plaque_arc(camber, t_rel, n=120):
    x = 0.5 * (1 - np.cos(np.linspace(0, pi, n)))
    yc = 4 * camber * x * (1 - x)
    yt = 0.5 * t_rel * np.ones_like(x)
    yt[:3] *= np.array([0.2, 0.6, 0.9]); yt[-3:] *= np.array([0.9, 0.6, 0.2])
    return np.vstack([np.column_stack([x[::-1], (yc + yt)[::-1]]), np.column_stack([x[1:], (yc - yt)[1:]])])


MARYLAND = [   # (nom, coords, CT/σ à l'optimum, FM mesurée) — Benedict et al. 2015, fig. 5-7
    ("plaque cambrée 6,1 % (t/c 2,3)", lambda: _plaque_arc(0.061, 0.023), 0.145, 0.59),
    ("NACA 4504", lambda: _naca4(0.045, 0.5, 0.04), 0.115, 0.575),
    ("NACA 4512", lambda: _naca4(0.045, 0.5, 0.12), 0.135, 0.53),
    ("NACA 0012 (contre-exemple)", lambda: _naca4(0.0, 0.5, 0.12), 0.085, 0.45),
]


def _polaire_nf(coords, re, ncrit=9.0):
    """Polaire c81 depuis NeuralFoil au Reynolds de l'essai, prolongée en plateau."""
    import neuralfoil as nf
    al = np.arange(-8.0, 16.01, 0.5)
    out = nf.get_aero_from_coordinates(coordinates=coords, alpha=al, Re=re, model_size="xxlarge", n_crit=ncrit)
    cl, cd = np.asarray(out["CL"]), np.asarray(out["CD"])
    ok = np.asarray(out["analysis_confidence"]) >= 0.7
    al, cl, cd = al[ok], cl[ok], cd[ok]
    a_ext = np.concatenate([[-180, -90, -30], al, [30, 90, 180]])
    cl_ext = np.concatenate([[0, 0, cl[0] * 0.5], cl, [cl[-1] * 0.5, 0, 0]])
    cd_ext = np.concatenate([[0.02, 1.5, 0.6], cd, [0.6, 1.5, 0.02]])
    return (list(a_ext), list(cl_ext), list(cd_ext), [0.0] * len(a_ext))


def _maryland_lignes(nom, fm_mes, fms, notes, re):
    """Publie chaque hypothèse calculée, sans sélection selon la mesure."""
    for nappe, lib in ((False, "sans nappe"), (True, "AVEC nappe")):
        for ncrit in (5.0, 9.0):
            _ligne(f"C · micro-rotor Maryland, {nom}", "Benedict et al. 2015 (mesuré)",
                   f"FM (sillage {lib}, N_crit {ncrit:g})", fm_mes, fms[nappe][ncrit], 0.10, re=re,
                   note=notes[nappe][ncrit] + " — hypothèse publiée séparément, sans sélection sur la mesure")


def cas_maryland(rapide):
    """FM mesurée des micro-rotors Maryland, rotor rigide du noyau, inflow uniforme.

    ATTENDU AVANT MESURE : l'inflow uniforme du noyau (théorie du disque, sans
    perte de bout) rend une puissance induite PLUS BASSE que le réel — ki
    effectif de ces rotors mesuré 1,31–1,38 (gate V-2 de FRELON) contre 1,0
    ici — donc une FM trop HAUTE, de l'ordre de +10 à +25 %. Le tableau
    publie l'écart et le ki qu'il faudrait ; c'est le module `sillage` (mean
    inflow du sillage libre) et une perte de bout qui le fermeraient.

    MESURÉ (5 sept.) : avec l'inflow du sillage libre SANS nappe (tourbillons
    marginaux, Γ uniforme), FM à −2,9 / −3,9 / −12,7 / −5,5 % de la mesure —
    dans la bande ±10 % pour trois rotors sur quatre.

    ET LA NAPPE INTERNE DÉFAIT L'ACCORD (6 sept.) : le sillage AVEC nappe
    (`nappe=True`, le défaut de `Sillage` — k₁ au bon signe, contraction à
    +3 % de Landgrebe, c'est-à-dire la MEILLEURE géométrie contre la mesure
    indépendante de Landgrebe) rend FM +9,9 / +15,8 / +19,8 / −14,5 % :
    trois rotors sur quatre HORS bande, et l'écart change de signe. Sur ces
    rotors très chargés (C_T/σ 0,085–0,145, Re 3e4) la nappe prescrite au
    pas de l'inflow local reste près du disque et y induit un upwash
    intérieur : moins de puissance induite, FM trop haute. Les deux lignes
    sont publiées, aucune tolérance n'a bougé, et le désaccord est déclaré
    NON RÉSOLU : l'accord sans nappe peut être une compensation (pas de
    nappe = pas d'upwash intérieur = un ki qui tombe juste), et la nappe
    prescrite n'est pas une nappe libre. Ce qui trancherait : une nappe
    LIBRE (×6 sur les inconnues, non fait) ou une mesure d'inflow au disque
    sur un rotor de cette classe — il n'y en a pas de publiée à ce Reynolds.
    """
    from vinkulum.rotor import Rotor
    from vinkulum.sillage import Sillage
    R, nb, c, rho = 0.041, 2, 0.0113, 1.225
    om = 10000.0 * 2 * pi / 60
    re = om * 0.75 * R * c / 1.5e-5
    sig = nb * c / (pi * R)
    rbar = (np.arange(10) + 0.5) / 10
    n_az = 16
    az = (np.arange(n_az) + 0.5) / n_az * 2 * pi
    pts = np.column_stack([(np.outer(rbar * R, np.cos(az))).ravel(), (np.outer(rbar * R, np.sin(az))).ravel(),
                           np.zeros(rbar.size * n_az)])
    # LA POLAIRE EST L'INCONNUE, PAS LE SOLVEUR : à Re 30 000 NeuralFoil rend, pour
    # les sections ÉPAISSES, une traînée qui dépend fortement de N_crit (bulle
    # laminaire) — la gate V-2 de FRELON le publie en bande {5, 9} et juge la
    # bande contre ±10 % de la FM mesurée. Ici chaque N_crit est désormais
    # jugé séparément : retenir le meilleur après lecture de la mesure
    # confondait compatibilité de la bande et pouvoir prédictif du modèle.
    for nom, coords, cts, fm_mes in (MARYLAND if not rapide else MARYLAND[:2]):
        fms, notes = {False: {}, True: {}}, {False: {}, True: {}}
        for ncrit in (5.0, 9.0):
            try:
                pol = _polaire_nf(coords(), re, ncrit)
            except Exception as ex:  # noqa: BLE001 — branche N_crit morte (confiance)
                for nappe in (False, True):
                    fms[nappe][ncrit] = float("nan")
                    notes[nappe][ncrit] = f"ERREUR N_crit {ncrit:g} : {ex}"
                continue
            r = Rotor(n_pales=nb, R=R, r0=0.025 * R, corde=c, omega=om, vrillage=0.0, nu_beta=3.0,
                      zeta=0.05, lock=4.0, segments=8, polaire=pol, rho=rho)
            t_cible = cts * sig * rho * pi * R ** 2 * (om * R) ** 2
            h = r.t_tour / 300

            def torseur(th0, profil=None):
                N = r.monte((th0, 0.0, 0.0), t_end=4 * r.t_tour)
                if profil is not None:
                    # l'inflow du SILLAGE LIBRE remplace le disque uniforme : v(r̄) figé
                    N.pose_inflow(r.inflow_idx, 0.0, impose=True)
                    N.pose_inflow_profil(r.inflow_idx, list(rbar), list(profil), [0.0, 0.0, 0.0],
                                         [1.0, 0.0, 0.0], om * R)
                N.moyenne_aero(3 * r.t_tour)
                N.simule(4 * r.t_tour, h, tous=10 ** 9)
                f, m, _ = N.torseur_moyen()
                return f[2], -m[2]

            def fm_a_poussee(profil=None):
                lo, hi = 2.0, 20.0
                for _ in range(14):
                    mid = 0.5 * (lo + hi)
                    T, Q = torseur(mid, profil)
                    lo, hi = (mid, hi) if T < t_cible else (lo, mid)
                T, Q = torseur(0.5 * (lo + hi), profil)
                ct = T / (rho * pi * R ** 2 * (om * R) ** 2)
                cp = abs(Q) * om / (rho * pi * R ** 2 * (om * R) ** 3)
                return ct ** 1.5 / (sqrt(2) * cp), 0.5 * (lo + hi), T, ct, cp
            fm_u, th_u, T0, ct, cp = fm_a_poussee()
            # SILLAGE LIBRE STATIONNAIRE : Γ depuis la poussée, géométrie par Newton,
            # profil radial v(r̄) moyenné en azimut, imposé au noyau ; deux passes
            # (la poussée à θ₀ trimé ne bouge plus qu'au pour-cent). Les DEUX
            # sillages sont jugés — sans nappe (5 sept.) et avec nappe (6 sept.) —
            # parce qu'ils ne s'accordent PAS sur ces rotors, cf. le docstring.
            for nappe in (False, True):
                sil = Sillage(R, nb, om, c, r_pied=0.025 * R, dzeta_deg=15.0 if rapide else 12.0,
                              tours_libres=3, nappe=nappe)
                prof, T = None, T0
                for _ in range(2):
                    sil.circulation(T)
                    sil.resout()
                    vz = sil.vitesse(pts)[:, 2].reshape(rbar.size, n_az)
                    prof = -vz.mean(axis=1)
                    fm_w, th_w, T, ct_w, cp_w = fm_a_poussee(prof)
                fms[nappe][ncrit] = fm_w
                cp_i = ct_w ** 1.5 / sqrt(2)
                notes[nappe][ncrit] = (f"N_crit {ncrit:g} : FM uniforme {fm_u:.3f} → sillage {fm_w:.3f} "
                                      f"(ki sillage {(cp_w - (cp - cp_i)) / cp_i:.2f}), θ₀ {th_w:.1f}°")
        # Aucun choix après comparaison : les deux hypothèses, déclarées
        # avant calcul, gardent chacune leur valeur et leur verdict.
        _maryland_lignes(nom, fm_mes, fms, notes, re)


# ── D. Lock modal ────────────────────────────────────────────────────────────
def cas_lock(rapide):
    from vinkulum import pale_elastique as P
    am = P.amortissement_lock(rapide)
    _ligne("D · amortissement aéro du battement élastique", "Johnson 1980 ch. 9 (Lock modal)", "ζ premier mode",
           am["zeta_th"], am["zeta"], 0.15)


# ── E. les références déjà jugées, rejouées ──────────────────────────────────
def cas_references(rapide):
    from vinkulum import verification as V
    for nom, f, src in (("Kapitza 1951", V.kapitza, "Kapitza 1951 (critère a²ω² > 2gL)"),
                        ("Bennett 1903", V.bennett, "Bennett 1903 (mobilité 1)"),
                        ("srscm (deck MBDyn)", V.srscm, "MBDyn, jeu de tests"),
                        ("6barmech (deck MBDyn)", V.six_barres, "MBDyn, jeu de tests"),
                        ("flambement latéral", V.lateral_buckling, "MSD 37 (2016) 29–48"),
                        ("arbre tournant", V.arbre_tournant, "MSD 37 (2016) 29–48"),
                        ("quatre-barres flexible", V.quatre_barres_flexible, "MSD 37 (2016) 29–48")):
        if rapide and nom in ("arbre tournant", "quatre-barres flexible"):
            continue
        t0 = time.time()
        try:
            f()
            LIGNES.append((f"E · {nom}", src, "banc `verification`", float("nan"), float("nan"), 0.0, True,
                           f"rejoué, vert, {time.time() - t0:.0f} s", None))
        except Exception as ex:  # noqa: BLE001 — préserver les autres références
            _erreur(f"E · {nom}", src, "banc `verification`", ex)
    try:
        from vinkulum import andrews
        andrews.demo()
        LIGNES.append(("E · mécanisme d'Andrews", "Schiehlen 1990 (benchmark IUTAM)", "banc", float("nan"),
                       float("nan"), 0.0, True, "rejoué, vert", None))
    except Exception as ex:  # noqa: BLE001
        _erreur("E · mécanisme d'Andrews", "Schiehlen 1990 (benchmark IUTAM)", "banc", ex)


def _erreur(cas, source, grandeur, ex):
    """Une référence non exécutée ne doit être ni perdue ni recotée en écart."""
    LIGNES.append((cas, source, grandeur, float("nan"), float("nan"), 1.0, False,
                   f"ERREUR {type(ex).__name__} : {ex}", None))


# ── F. Coleman, si le sillage en avancement est là ───────────────────────────
def cas_coleman(rapide):
    """k_x du sillage libre en avancement contre Drees (1949) et Coleman (1945).

    Rotor FRELON (R 0,175, 2 pales, C_T 0,0069). Les deux sont des MODÈLES,
    pas des essais ; Drees et Coleman diffèrent entre eux de 10 à 25 %, et le
    k_x jugé est le gradient AU CENTRE (r̄ ≤ 0,5), comme leurs formules. Drees
    est asserté à 15 % ; Coleman publié.
    """
    from vinkulum.sillage import SillageAvancement, coleman
    R, n_b, omega, corde, r_pied, poussee = 0.175, 2, 334.9, 0.020, 0.15 * 0.175, 2.8
    dpsi, tl = (18.0, 2) if rapide else (15.0, 3)
    for mu in ((0.15,) if rapide else (0.15, 0.30)):
        s = SillageAvancement(R, n_b, omega, corde, mu, r_pied=r_pied, dpsi_deg=dpsi, tours_libres=tl, tours_loin=4)
        s.circulation(poussee)
        s.resout()
        hm, h5 = s.harmoniques(), s.harmoniques(rmax=0.5)
        chi, kx_c, kx_d = coleman(mu, hm["lam0"])
        _ligne(f"F · sillage libre μ = {mu:.2f}, gradient d'inflow", "Drees 1949", "k_x au centre", float(kx_d),
               h5["kx"], 0.15, re=omega * 0.75 * R * corde / 1.5e-5, note=f"Coleman 1945 : {kx_c:.3f} ({100 * (h5['kx'] / kx_c - 1):+.0f} %) · λ₀/Glauert {hm['lam0'] / s.glauert()[0]:.3f}")

def cas_elliott(rapide):
    """L'inflow MESURÉ d'Elliott et al. 1988 contre le sillage libre en avancement.

    Le seul essai de l'inflow en avancement de la campagne, et le classique du
    domaine (Bagai–Leishman l'ont utilisé pour leur sillage). Rotor 4 pales
    articulé, C_T 0,0064, μ 0,15 et 0,30 (`vinkulum.elliott1988`). On compare
    la vitesse induite NORMALE mesurée (une corde au-dessus du plan de
    l'arbre, moyennée en azimut) à trois modèles évalués aux mêmes points :
    Glauert uniforme, Drees linéaire (λ₀(1 + k_x r̄ cos ψ)), et le sillage
    libre (Γ uniforme par pale, moyenné sur ses phases, évalué une corde
    au-dessus du disque). Métrique : RMS de l'écart rapporté à la moyenne
    mesurée sur r̄ ≤ 1. Asserté : le sillage fait MIEUX que Drees ; les RMS
    absolus sont publiés, pas assertés — aucune tolérance de la littérature
    n'est assez précise pour en faire un seuil honnête.

    DEUX QUESTIONS TRANCHÉES PAR LA MESURE (6 sept.) : le SENS d'azimut — le
    miroir (ψ → −ψ) dégrade la RMS (0,285 → 0,299 à μ 0,15) et retourne k_y
    du mauvais côté : la convention directe est la bonne, et le k_y positif
    du sillage contre −0,06 mesuré reste un écart de MODÈLE (Γ uniforme) ;
    la HAUTEUR — au plan du disque la RMS vaut 0,46, une corde au-dessus
    0,285, deux cordes 0,26 : l'écart entre le λ₀ mesuré (0,002 à μ 0,30)
    et Glauert (0,011) est pour partie la sonde au-dessus du plan, que le
    sillage reproduit en tendance (0,0107 → 0,0088 → 0,0075) mais pas en
    amplitude. On évalue là où la sonde était : une corde.
    """
    from vinkulum import elliott1988 as E
    from vinkulum.sillage import SillageAvancement, coleman
    om = E.VTIP / E.R_M
    rho = 1.225
    for mu_cle in ((0.15,) if rapide else (0.15, 0.30)):
        cd = E.CONDITIONS[mu_cle]
        tab = E.table(mu_cle)
        tab = tab[tab[:, 1] <= 1.0]
        psi_m, r_m, v_m = np.radians(tab[:, 0]), tab[:, 1], -tab[:, 3]   # λ_i mesuré, > 0 vers le bas
        T = cd["ct"] * rho * pi * E.R_M ** 2 * E.VTIP ** 2
        s = SillageAvancement(E.R_M, E.N_PALES, om, E.CORDE_M, cd["mu"], r_pied=0.24 * E.R_M,
                              dpsi_deg=18.0 if rapide else 15.0, tours_libres=2 if rapide else 3, tours_loin=4,
                              alpha_disque=radians(-cd["alpha_arbre_deg"]))
        s.circulation(T)
        s.resout()
        # points de mesure dans le repère du sillage (ψ depuis l'aval, x = −r cos ψ), une corde au-dessus
        pts = np.column_stack([-r_m * np.cos(psi_m) * E.R_M, -r_m * np.sin(psi_m) * E.R_M,
                               np.full(r_m.size, E.CORDE_M)])
        vz = np.zeros(pts.shape[0])
        for k in range(s.n_psi):
            vz += s.vitesse(pts, k)[:, 2]
        lam_w = -vz / s.n_psi / E.VTIP
        lam0, ct = s.glauert()
        chi, kx_c, kx_d = coleman(cd["mu"], lam0)
        lam_u = np.full(r_m.size, lam0)
        lam_d = lam0 * (1.0 + kx_d * r_m * np.cos(psi_m))
        # échelle : la RMS de la MESURE elle-même — sa moyenne vaut 0,002 à μ 0,30
        # (upwash sur tout l'avant, une corde au-dessus du plan), rapporter à
        # la moyenne y ferait exploser la métrique sans rien juger
        sig = float(np.sqrt(np.mean(v_m ** 2)))
        rms = lambda x: float(np.sqrt(np.mean((x - v_m) ** 2)) / sig)
        r_u, r_d, r_w = rms(lam_u), rms(lam_d), rms(lam_w)
        # le gradient longitudinal mesuré, par le même ajustement que pour le sillage
        A = np.column_stack([np.ones(r_m.size), r_m * np.cos(psi_m), r_m * np.sin(psi_m)])
        c_m = np.linalg.lstsq(A, v_m, rcond=None)[0]
        c_w = np.linalg.lstsq(A, lam_w, rcond=None)[0]
        print(f"║   Elliott μ {cd['mu']:.3f} : λ₀ mesuré {c_m[0]:.4f} (Glauert {lam0:.4f}) · k_x mesuré {c_m[1] / c_m[0]:+.3f} "
              f"(Drees {kx_d:+.3f}, sillage {c_w[1] / c_w[0]:+.3f}) · k_y mesuré {c_m[2] / c_m[0]:+.3f} (sillage {c_w[2] / c_w[0]:+.3f})")
        re_e = E.VTIP * 0.75 * E.CORDE_M / 1.5e-5
        # PETERS–HE (7 sept.) : états finis M 4 / Q 6, charge de la théorie des
        # tranches d'un rotor RIGIDE trimé à C_T (θ₀ ajusté, θ1c/θ1s mesurés),
        # évalué AU DISQUE — la sonde d'Elliott est une corde au-dessus, et ce
        # décalage vaut 0,46 → 0,285 de RMS sur le sillage libre : PH est donc
        # jugé avec un handicap déclaré, et publié à côté des trois autres
        from vinkulum import peters_he as PH
        ph = PH.PetersHe(M=4, Q=6)
        th0, ct_ph = PH.trim_ct(ph, cd["ct"], E.R_M, E.N_PALES, om, E.CORDE_M, cd["mu"], radians(-cd["alpha_arbre_deg"]),
                                radians(cd["theta1c_deg"]), radians(cd["theta1s_deg"]), radians(E.VRILLAGE_DEG),
                                r_pied=0.24, rho=rho)
        lam_ph = ph.w(r_m, psi_m)
        r_ph = rms(lam_ph)
        c_ph = np.linalg.lstsq(A, lam_ph, rcond=None)[0]
        print(f"║   Peters–He μ {cd['mu']:.3f} : θ₀ trimé {degrees(th0):.2f}° (mesuré {cd['a0_deg']}), C_T {ct_ph:.5f}, "
              f"λ₀ {c_ph[0]:.4f}, k_x {c_ph[1] / c_ph[0]:+.3f}, k_y {c_ph[2] / c_ph[0]:+.3f}, RMS/σ {r_ph:.3f}")
        _ligne(f"G · inflow mesuré (Elliott 1988) μ = {cd['mu']:.2f}", "NASA TM 100541/100543 (mesuré, LV)",
               "RMS/σ_mes Peters–He (M 4, Q 6, au disque)", r_d, r_ph, 10.0, re=re_e,
               note=f"Glauert {r_u:.3f}, Drees {r_d:.3f}, sillage {r_w:.3f}, Peters–He {r_ph:.3f} ; "
                    f"θ₀ trimé {degrees(th0):.2f}° pour {cd['a0_deg']} mesuré ; k_x PH {c_ph[1] / c_ph[0]:+.2f} — rotor rigide, au disque")
        assert r_ph < r_u, ("Peters–He devrait battre Glauert uniforme", r_ph, r_u)
        _ligne(f"G · inflow mesuré (Elliott 1988) μ = {cd['mu']:.2f}", "NASA TM 100541/100543 (mesuré, LV)",
               "RMS/σ_mes sillage libre", r_d, r_w, 10.0, re=re_e,
               note=f"RMS/σ_mes : Glauert {r_u:.3f}, Drees {r_d:.3f}, sillage {r_w:.3f} ({tab.shape[0]} points, r̄ ≤ 1) ; "
                    f"λ₀ mesuré {c_m[0]:.4f}, Glauert {lam0:.4f} — asserté : sillage < Drees")
        assert r_w < r_d, (r_w, r_d)
        if mu_cle == 0.15:
            # à μ 0,30 le λ₀ mesuré une corde au-dessus vaut 0,002 (Glauert 0,011) :
            # le gradient rapporté à λ₀ n'a plus de sens, il est publié dans la note
            _ligne(f"G · inflow mesuré (Elliott 1988) μ = {cd['mu']:.2f}", "NASA TM 100541/100543 (mesuré, LV)",
                   "k_x mesuré vs sillage", c_m[1] / c_m[0], c_w[1] / c_w[0], 0.30, re=re_e,
                   note=f"Drees {kx_d:+.3f} ({100 * (kx_d / (c_m[1] / c_m[0]) - 1):+.0f} %)")


def cas_s809(rapide):
    """Décrochage dynamique MESURÉ : S809 en tangage (OSU/NREL 1995) contre Leishman–Beddoes et Øye.

    Polaire statique = table B3 transcrite ; tangage ±5,5° à k = 0,077 autour
    de 8°, 14°, 20° ; on compare le CL MAX de la boucle et le CL MIN de la
    descente lus sur la figure C25 (±0,05). Depuis le 6 sept. la section
    TANGUE autour de c/4 (cible de rotation en vitesse, `_tangage`), avec la
    masse ajoutée et l'incidence au 3/4 de corde — le harnais à vent tournant
    n'avait ni α̇ ni α̈. Les constantes de LB sont celles de la littérature (NACA 0012,
    Leishman–Beddoes 1989), non calibrées sur le S809 — un profil de 21 %
    d'épaisseur à décrochage de bord de fuite : l'écart publié est celui du
    modèle générique, et c'est ce qu'un utilisateur verra.
    """
    from math import cos, radians, sin
    from vinkulum import osu_s809 as S
    c, V = 0.1, 20.0
    om = 2.0 * S.K_RED * V / c
    T = 2 * pi / om
    pol = S.polaire_c81()

    def boucle(am, aa, wagner=False, **kw):
        # section qui TANGUE autour de c/4 (l'axe de l'essai OSU), cible en
        # vitesse ; montée douce 0 → α_m sur 0,05 s puis n_cyc cycles
        n_cyc = 2 if rapide else 3
        dt = T / 400
        t_pre = 0.05
        tt = np.arange(0.0, t_pre + n_cyc * T + dt, dt)
        al = np.where(tt < t_pre, am * 0.5 * (1 - np.cos(pi * tt / t_pre)), am + aa * np.sin(om * (tt - t_pre)))
        tab = list(zip(tt.tolist(), np.radians(al).tolist())) + [(tt[-1] + 1.0, float(np.radians(al[-1])))]
        N, q, _ = _tangage(c, V, pol, tab, 0.25, **kw)
        if wagner:
            N.instationnaire(0)          # le module attaché de LB = Wagner : les deux ensemble
        t, pts = 0.0, []
        while t < t_pre + n_cyc * T - 1e-12:
            t += dt
            N.simule(t, dt, tous=10 ** 9)
            if t > t_pre:
                a = am + aa * sin(om * (t - t_pre))
                pts.append((a, _cl_tangage(N, q), om * (t - t_pre)))
        pts = np.array(pts[-400:])
        desc = pts[np.cos(pts[:, 2]) < 0]        # descente : dα/dt < 0
        return float(pts[:, 1].max()), float(desc[:, 1].min())

    # LE HARNAIS D'ABORD : quasi-statique à 8° ± 5,5°, le CL max doit être celui de
    # la table B3 à 13,5° (1,03 ± 0,03 interpolé) — sinon on ne juge pas un modèle
    mx_qs, _ = boucle(8.0, S.AMPLITUDE_DEG)
    cl_tab = float(np.interp(13.5, S.STATIQUE_B3[:, 0], S.STATIQUE_B3[:, 1]))
    assert abs(mx_qs - cl_tab) < 0.03, ("harnais S809 : le quasi-statique ne rend pas la table", mx_qs, cl_tab)
    for am, ref in ((S.DYNAMIQUE_C25.items()) if not rapide else [(20.0, S.DYNAMIQUE_C25[20.0])]):
        mx_lb, mn_lb = boucle(am, S.AMPLITUDE_DEG, lb=True, wagner=True, masse_ajoutee=True, alpha_34=True)
        mx_oy, mn_oy = boucle(am, S.AMPLITUDE_DEG, oye=True, wagner=True, masse_ajoutee=True, alpha_34=True)
        mx_st, mn_st = boucle(am, S.AMPLITUDE_DEG)
        _ligne(f"H · S809 tangage {am:.0f}° ± 5,5°, k 0,077 (OSU 1995)", "NREL/TP-442-7817 fig. C25 (mesuré)",
               "CL max de la boucle", ref["cl_max"], mx_lb, 0.12, re=1.01e6,
               note=f"Øye+Wagner {mx_oy:.2f}, quasi-statique {mx_st:.2f} — LB+Wagner+masse ajoutée+3/4, pivot c/4, constantes NACA 0012, ±0,05 de lecture")
        _ligne(f"H · S809 tangage {am:.0f}° ± 5,5°, k 0,077 (OSU 1995)", "NREL/TP-442-7817 fig. C25 (mesuré)",
               "CL min en descente", ref["cl_min"], mn_lb, 0.20, re=1.01e6,
               note=f"Øye+Wagner {mn_oy:.2f}, quasi-statique {mn_st:.2f}")


# ── I. plaque en rampe 0 → 45°, Re 2e4 (Jantzen–Granlund–Ol 2014) ───────────
def _section(c, V, pol, rapide_dt, lb=False, wagner=False, oye=False):
    """Section fixe dont le VENT tourne (patron du S809) ; rend (N, vent, cl_de)."""
    from math import cos, radians, sin
    N = Noyau([0, 0, 0])
    b = N.corps("s", 1.0, list(np.eye(3).ravel()), [0, 0, 0])
    N.liaison("fixe", None, b)
    kw = {}
    if lb:
        kw["lb"] = True
    if oye:
        kw["oye"] = True
    N.pale("p", b, [0, 0, 0], [0, 1, 0], [1, 0, 0], 1.0, c, pol, **kw)
    if wagner:
        N.instationnaire(0)
    q = 0.5 * 1.225 * V * V * c

    def vent(a):
        N.vent([-V * cos(radians(a)), 0.0, -V * sin(radians(a))])

    def cl_de(fr, a):
        f = -np.array(fr[5][:3])                      # force aéro = −réaction
        return float(np.dot(f, [-sin(radians(a)), 0.0, cos(radians(a))])) / q   # ⟂ au vent
    return N, vent, cl_de


def _tangage(c, V, pol, table, pivot_xc, **pale_kw):
    """Section qui TANGUE (cible de rotation du corps, table (t, α_rad)), vent fixe.

    Le pivot est à `pivot_xc` de la corde depuis le bord d'attaque ; la ligne
    de référence de la pale (c/4) est posée à (0,25 − x_p)·c derrière lui.
    Rend (N, q) ; la portance se lit sur la sonde d'inflow (axe +z, τ = 1e9 :
    v_i reste nul), `N.aero()[0][0][1]` — c'est la force ⟂ au vent, exacte à
    tout angle. Le harnais du S809 (vent tournant) n'avait ni α̇ ni α̈ ; celui-ci
    les a, et c'est ce que la masse ajoutée et le 3/4 de corde demandent.
    """
    N = Noyau([0, 0, 0])
    N.vent([-V, 0.0, 0.0])
    b = N.corps("s", 1.0, [1e-3, 0, 0, 0, 1e-3, 0, 0, 0, 1e-3], [0, 0, 0])
    # NON HOLONOME (nh) : la cible est imposée au niveau VITESSE — ω est
    # EXACTEMENT la pente de la table sur chaque pas. Imposée en position
    # (index 3), la vitesse du corps est une quantité dérivée qui oscille d'un
    # pas à l'autre, et la masse ajoutée — une différence arrière de ẇ — en
    # fait un pic par pas : mesuré sur la rampe de Jantzen, C_L max 2 974 pour
    # 6. En vitesse, la différence arrière rend α̈ en moyenne, ce qu'elle doit.
    # AXE +y et lecture −sonde : avec es = y, ec = x la normale en = es × ec
    # pointe vers −z, donc une rotation nose-up autour de +y donne +α à la
    # section et sa portance sort le long de −z. Le contraire (axe −y, +sonde)
    # passait tous les contrôles sur polaire SYMÉTRIQUE — la section voyait −α
    # et la sonde lisait −L, deux signes qui s'annulent — et rendait 0,257 pour
    # 0,447 sur le S809 cambré (6 sept.). D'où `_tangage_signe` sur le S809.
    N.liaison("pivot", None, b, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2],
              cible_r=([0, 1, 0], ("table", [x for pr in table for x in pr])), nh=True)
    j = N.inflow([0.0, 0.0, 1.0], 1.0, 1e9, 1.225)
    N.pale("p", b, [-(0.25 - pivot_xc) * c, -0.5, 0.0], [0, 1, 0], [1, 0, 0], 1.0, c, pol,
           inflow=j, rho=1.225, gauss=5, **pale_kw)
    return N, 0.5 * 1.225 * V * V * c, j


def _cl_tangage(N, q):
    """Portance ⟂ au vent, depuis la sonde (axe +z) — signe : la normale de la
    section est −z, cf. `_tangage`."""
    return -N.aero()[0][0][1] / q


def _tangage_signe(c, V, pol=None):
    """Le signe de la rotation et de la lecture, jugés sur une polaire CAMBRÉE
    (S809) à deux angles — une polaire symétrique laisse passer deux signes
    faux qui s'annulent, mesuré le 6 sept."""
    from vinkulum import osu_s809 as S
    pol = S.polaire_c81()
    for a in (3.0, 13.5):
        tab = [(0.0, 0.0), (0.02, radians(a)), (1.0, radians(a))]
        N, q, _ = _tangage(c, V, pol, tab, 0.25)
        N.simule(0.2, 1e-4, tous=10 ** 9)
        cl = _cl_tangage(N, q)
        ref = float(np.interp(a, pol[0], pol[1]))
        assert abs(cl - ref) < 0.03, ("harnais de tangage : le quasi-statique ne rend pas la table", a, cl, ref)
    return 1.0


def cas_jantzen(rapide):
    """Plaque plane en RAMPE 0 → 45° à Re 20 000 (Jantzen, Taira, Granlund & Ol 2014, fig. 13, 2D, mesuré).

    Le noyau tourne la section (vent tournant, patron S809) suivant la rampe
    d'Eldredge, LB + Wagner sur une polaire statique de plaque CONSTRUITE
    (`jantzen2014.polaire_plaque_c81`, HYPOTHESE : le papier n'en publie pas).
    Jugé : la MOYENNE de C_L sur une fenêtre après la rampe — [2 ; 8] pour C1,
    [3 ; 6] (le pic pendant la rampe) et [6 ; 8] pour C6 — parce qu'à 45° la
    plaque lâche des tourbillons et qu'un instant ne veut rien dire. Le pic
    NON CIRCULATOIRE de C1 (C_L 6 à t 0,85, masse ajoutée) est exclu : le noyau
    n'a pas ce terme, déclaré. Le pivot est au bord d'attaque dans l'essai ; la
    vitesse induite par le taux de tangage n'est pas dans le harnais (le vent
    tourne, la section ne bouge pas) — à K π/8 c'est un biais déclaré, pas un
    détail. La sensibilité au niveau post-décrochage de la polaire (±15 % sur
    CN) est publiée dans la note.
    """
    from vinkulum import jantzen2014 as J
    c, V = 0.1, 20.0
    conv = c / V
    pol = J.polaire_plaque_c81()
    _tangage_signe(c, V, pol)

    def rampe(cas, pol, **kw):
        t0, dt = 0.02, conv / 200
        # un nœud de table PAR PAS : la cible est linéaire par morceaux, donc la
        # vitesse imposée est constante entre deux nœuds et SAUTE à chaque nœud ;
        # avec la masse ajoutée (différence arrière de ẇ) un nœud tous les cinq
        # pas fait un pic par nœud — mesuré : C_L max 0,72 pour 0,49 attendu à
        # k 0,76 sur polaire linéaire ; un nœud par pas rend 0,51. Le banc
        # Theodorsen ne le voit pas : il lit une amplitude de Fourier, qui filtre.
        tt = np.arange(0.0, t0 + 8.3 * conv, dt)
        al = np.where(tt < t0 + 0.25 * conv, 0.0,
                      np.radians(J.alpha_deg((tt - t0 - 0.25 * conv) / conv, cas)))
        tab = list(zip(tt.tolist(), al.tolist()))
        tab.append((tt[-1] + 1.0, float(al[-1])))
        N, q, _ = _tangage(c, V, pol, tab, 0.0, **kw)
        if kw.get("lb") or kw.get("retard_lineaire") or kw.get("masse_ajoutee"):
            N.instationnaire(0)                # Wagner : le module attaché de LB
        N.simule(t0, 1e-4, tous=10 ** 9)
        tc, cl = [], []
        t = t0
        while t < t0 + 8.25 * conv - 1e-12:
            t += dt
            N.simule(t, dt, tous=10 ** 9)
            tc.append((t - t0 - 0.25 * conv) / conv)
            cl.append(_cl_tangage(N, q))
        tc, cl = np.array(tc), np.array(cl)
        return (lambda a_, b_: float(cl[(tc >= a_) & (tc <= b_)].mean())), (tc, cl)

    def pol_scal(f):
        a, cl, cd, cm = pol
        a = np.array(a)
        return list(a), [v * (f if abs(x) > 12 else 1.0) for x, v in zip(a, cl)], cd, cm
    complet = dict(lb=True, masse_ajoutee=True, alpha_34=True)
    for cas, fenetres in (("C1", ((2.0, 8.0),)), ("C6", ((3.0, 6.0), (6.0, 8.0)))):
        if rapide and cas == "C6":
            continue
        m_lb, (tc, cl_t) = rampe(cas, pol, **complet)
        m_st, _ = rampe(cas, pol)
        m_sans, _ = rampe(cas, pol, lb=True)
        m_lo, _ = rampe(cas, pol_scal(0.85), **complet)
        m_hi, _ = rampe(cas, pol_scal(1.15), **complet)
        for a_, b_ in fenetres:
            ref = J.moyenne(cas, a_, b_)
            _ligne(f"I · plaque plane rampe 0→45°, {cas} (K {J.K_RED[cas]:.3f}), Re 2e4 (AFRL 2014)",
                   "Jantzen, Taira, Granlund & Ol 2014 fig. 13 (mesuré, 2D)", f"C_L moyen sur t ∈ [{a_:g} ; {b_:g}]",
                   ref, m_lb(a_, b_), 0.25, re=J.RE,
                   note=f"quasi-statique {m_st(a_, b_):.2f} ; LB+Wagner sans masse ajoutée ni 3/4 {m_sans(a_, b_):.2f} ; "
                        f"CN post-décrochage ±15 % : {m_lo(a_, b_):.2f}–{m_hi(a_, b_):.2f} — LB+Wagner+masse ajoutée+3/4, "
                        "pivot au BA, polaire statique CONSTRUITE")
        if cas == "C1":
            # le pic NON CIRCULATOIRE pendant l'accélération (t < 1,2) : c'est la
            # masse ajoutée, et elle est maintenant dans le modèle
            m = (tc > 0.0) & (tc < 1.2)
            pic = float(cl_t[m].max())
            mex = (J.T_CONV > 0.0) & (J.T_CONV < 1.2)
            _ligne("I · plaque plane rampe 0→45°, C1 (K 0.393), Re 2e4 (AFRL 2014)",
                   "Jantzen, Taira, Granlund & Ol 2014 fig. 13 (mesuré, 2D)", "pic non circulatoire C_L max, t < 1,2",
                   float(J.CL_C1[mex].max()), pic, 0.30, re=J.RE,
                   note="masse ajoutée πρb²ẇ, différence arrière, rampe d'Eldredge a = 11 — le pic dépend du lissage a (déclaré)")


# ── J. NACA 0012 en tangage sinusoïdal, Re 2,3e4 (Kim, Chang & Sohn 2017) ───
def cas_kim(rapide):
    """NACA 0012, α = 6°·sin(ωt) autour de c/4, Re 2,3e4, k 0,1 / 0,2 / 0,4 / 0,76 (mesuré, IJASS 2017).

    Polaire statique NeuralFoil au Reynolds de l'essai (N_crit 9 — la bulle
    laminaire y est, c'est elle qui fait le saut de 4° à 6°) ; LB + Wagner.
    Jugés : C_L max de la boucle (imprimé sur la figure, ±0,03) et le SENS de
    rotation de la boucle — anti-horaire jusqu'à k 0,4, horaire à 0,76 —
    qu'on mesure comme le signe de l'aire orientée ∮ C_L dα. Le sens est
    publié comme une ligne booléenne (réf 1 = même sens). Comme au S809 le
    vent tourne et la section ne bouge pas : à k 0,76 la vitesse de tangage
    manque au modèle, déclaré.

    MESURÉ AVANT DE JUGER — le harnais, polaire LINÉAIRE 2π, α 3°·sin ωt,
    Wagner du noyau contre Theodorsen (|C(k)|, phase) :
        k 0,05 : 0,910 / −8,8°  (0,918 / −8,2°)     k 0,4  : 0,648 / −16,2°  (0,646 / −14,8°)
        k 0,1  : 0,845 / −11,3° (0,850 / −11,7°)    k 0,76 : 0,559 / −12,9°  (0,571 / −12,2°)
        k 0,2  : 0,764 / −14,6° (0,752 / −14,5°)
    Le circulatoire est juste à 2 % à tous les k de l'essai. Ce qui casse sur
    Kim est donc ailleurs, et c'est mesuré : (1) la polaire NeuralFoil à Re
    2,3e4 porte le SAUT de la bulle laminaire (C_L 0,075 à 4°, 0,535 à 6°), et
    le noyau lit la table statique à l'incidence EFFECTIVE de Wagner — à k 0,76
    l'incidence effective plafonne à 0,57 × 6° = 3,4°, sous le saut, C_L ~0,05 ;
    l'essai dit l'inverse (« les événements de couche limite disparaissent à
    k ≥ 0,4 », C_L max 0,41 ≈ 3,7/rad × 6°) ; (2) la MASSE AJOUTÉE (non
    circulatoire, ∝ α̇) n'est pas modélisée : à k 0,76 elle vaut ~π·k·α_a ≈ 0,25,
    en quadrature — la moitié du C_L max mesuré. Les deux sont des limites du
    MODÈLE à bas Re, publiées.

    CLOS LA MÊME NUIT, en partie : masse ajoutée + 3/4 de corde + retard
    linéaire (a₀ 2π) — k 0,76 passe de −69 % à −3 %, k 0,1 reste à −34 %
    (la bulle porte plus en descente qu'en montée : aucun modèle de section
    ne l'a). Le SENS de la boucle : 3 sur 4 — le modèle complet renverse
    entre k 0,2 et 0,4 (la masse ajoutée mène en phase, comme Theodorsen),
    la mesure entre 0,4 et 0,76 ; le Wagner à table retardée rendait 4/4
    par un retard qui n'était pas le bon. Ce décalage visqueux du
    renversement EST la bulle, encore elle.
    """
    from vinkulum import kim2017 as K
    c, V = 0.1, 20.0
    pol = _polaire_nf(_naca4(0.0, 0.5, 0.12), K.RE, 9.0)
    _tangage_signe(c, V, pol)

    def boucle(k, **kw):
        om = 2.0 * k * V / c
        T = 2 * pi / om
        n_cyc = 3
        dt = T / 400
        tt = np.arange(0.0, n_cyc * T + dt, dt)      # un nœud de table PAR PAS (cf. cas_jantzen)
        al = np.radians(K.ALPHA_M_DEG + K.ALPHA_A_DEG * np.sin(om * tt))
        tab = list(zip(tt.tolist(), al.tolist())) + [(tt[-1] + 1.0, float(al[-1]))]
        instat = kw.pop("wagner", False)
        N, q, _ = _tangage(c, V, pol, tab, K.PIVOT_X_C, **kw)
        if instat:
            N.instationnaire(0)
        t, pts = 0.0, []
        while t < n_cyc * T - 1e-12:
            t += dt
            N.simule(t, dt, tous=10 ** 9)
            a = K.ALPHA_M_DEG + K.ALPHA_A_DEG * sin(om * t)
            pts.append((a, _cl_tangage(N, q)))
        p = np.array(pts[-400:])
        aire = 0.5 * float(np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1]))
        return float(p[:, 1].max()), (1 if aire > 0 else -1), aire

    complet = dict(wagner=True, masse_ajoutee=True, alpha_34=True, retard_lineaire=True, a0_rad=2 * pi)
    ks = (0.10, 0.76) if rapide else sorted(K.FIG8)
    for k in ks:
        clmax, sens, aire = boucle(k, **dict(complet))
        mx_lb, s_lb, _ = boucle(k, wagner=True, lb=True, masse_ajoutee=True, alpha_34=True)
        mx_w, s_w, _ = boucle(k, wagner=True)
        mx_st, s_st, _ = boucle(k)
        ref_max, _, ref_sens = K.FIG8[k]
        _ligne(f"J · NACA 0012 tangage 0° ± 6°, k {k:g}, Re 2,3e4 (IJASS 2017)", "Kim, Chang & Sohn 2017 fig. 8 (mesuré)",
               "C_L max de la boucle", ref_max, clmax, 0.12, re=K.RE,
               note=f"LB+Wagner+masse ajoutée+3/4 {mx_lb:.2f} ; Wagner seul (table à α retardée) {mx_w:.2f} ; quasi-statique {mx_st:.2f} — "
                    "jugé : Wagner+masse ajoutée+3/4+retard linéaire (a₀ 2π déclaré), polaire NeuralFoil N_crit 9, ±0,03 de lecture")
        _ligne(f"J · NACA 0012 tangage 0° ± 6°, k {k:g}, Re 2,3e4 (IJASS 2017)", "Kim, Chang & Sohn 2017 fig. 8 (flèche)",
               "sens de la boucle (1 = mesuré)", 1.0, float(sens * ref_sens), 0.5, re=K.RE,
               note=f"aire orientée ∮C_L dα {aire:+.3f} ; mesuré {'anti-horaire' if ref_sens > 0 else 'HORAIRE'} ; "
                    f"LB {'anti-horaire' if s_lb > 0 else 'horaire'} ; Wagner seul {'anti-horaire' if s_w > 0 else 'horaire'}")


def rapport(chemin):
    lignes = ["# Campagne de validation — vinkulum contre des essais et modèles publiés", "",
              "Généré par `python -m vinkulum.validation` ; chaque ligne porte sa source, le chiffre publié, ",
              "le chiffre du noyau et l'écart. Un ÉCART est publié, jamais recoté.", "",
              "| cas | source | grandeur | référence | noyau | écart | verdict | Re | note |", "|---|---|---|---|---|---|---|---|---|"]
    for cas, src, gr, ref, val, e, ok, note, re in LIGNES:
        fr = lambda x: "—" if x != x else f"{x:.4f}"
        sre = "—" if re is None else (f"{re:.1e}" + (" **bas**" if re < 2e5 else ""))
        lignes.append(f"| {cas} | {src} | {gr} | {fr(ref)} | {fr(val)} | {'—' if e != e else f'{100 * e:+.2f} %'} | "
                      f"{'OK' if ok else 'ÉCART'} | {sre} | {note} |")
    n_ok = sum(1 for l in LIGNES if l[6])
    n_bas = sum(1 for l in LIGNES if l[8] is not None and l[8] < 2e5)
    lignes += ["", f"**{n_ok} / {len(LIGNES)} lignes dans leur tolérance.** "
                   f"{n_bas} lignes à BAS Reynolds (décade de FRELON, 5e4–8e4) ; les autres sont à haut Re ou sans "
                   "Reynolds (structure, forme fermée) et ne couvrent pas la bulle laminaire, le N_crit ni les "
                   "constantes de décrochage dynamique à bas Re — consigne du 6 sept.", ""]
    open(chemin, "w").write("\n".join(lignes) + "\n")


def main(rapide=False):
    LIGNES.clear()
    t0 = time.time()
    print("╔═ vinkulum — CAMPAGNE DE VALIDATION contre des références publiées")
    for f in (cas_wright, cas_johnson, cas_lock, cas_maryland, cas_coleman, cas_elliott, cas_s809, cas_jantzen, cas_kim, cas_references):
        try:
            f(rapide)
        except Exception as ex:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            _erreur(f.__name__, "", "", ex)
    d = depot_docs()
    if d:
        import os
        rapport(os.path.join(d, "docs", "VALIDATION.md"))
    n_ok = sum(1 for l in LIGNES if l[6])
    print(f"╚═ {n_ok}/{len(LIGNES)} lignes dans leur tolérance — {time.time() - t0:.0f} s")
    return LIGNES


if __name__ == "__main__":
    lignes = main(rapide="rapide" in sys.argv)
    # Un écart physique reste un résultat ; une campagne incomplète à cause
    # d'une erreur d'exécution doit, elle, faire échouer la commande.
    sys.exit(int(any(l[7].startswith("ERREUR") for l in lignes)))
