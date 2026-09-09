"""vinkulum.floquet — stabilité d'un système PÉRIODIQUE, par la monodromie.

    python -m vinkulum.floquet        # contrôles analytiques, puis le rotor

POURQUOI CE MODULE EXISTE, et quelle dette il solde.

`Noyau.modes_complexes` linéarise autour d'un état et rend (f, ζ, σ). Il est
exact à 2e-16 sur un oscillateur amorti — et il REFUSE dès qu'un corps
contraint tourne. La mesure du 3 septembre a dit pourquoi, et ce n'était pas
ce que j'avais d'abord écrit : la raideur tangente en repère fixe autour d'un
état tournant est dominée par des termes en ω²(J_t − J_a) valant mille fois la
raideur cherchée, et le résultat n'émerge que d'une annulation quasi exacte.

**La sortie n'est pas de mieux projeter, c'est de ne pas linéariser en repère
fixe.** HOST passe aux composantes périodiques de l'état — « In HOST, the
stability is done using the periodic components of the state variables » — et
en tire un système à coefficients CONSTANTS, ce qui est une approximation
commode et standard. Floquet en est la version exacte : on ne linéarise aucun
champ, on propage la dynamique RÉELLE sur une période.

    δx(T) = Φ · δx(0)        λ = eig(Φ)        σ = ln λ / T
    f = |Im σ| / 2π          ζ = −Re σ / |σ|   stable ⟺ |λ| ≤ 1

CE QUE ÇA EXIGE, et le module le VÉRIFIE au lieu de le supposer :
· la trajectoire de référence doit être PÉRIODIQUE de période T (sinon Φ ne
  veut rien dire) — contrôlé en mesurant l'écart entre l'état à t₀ et à t₀+T ;
· les directions perturbées doivent être ADMISSIBLES (dans le noyau des
  contraintes), sinon le premier pas dépense son énergie à réparer Φ et la
  colonne mesure la projection, pas la dynamique — contrôlé sur ‖Φ‖ et ‖Φ̇‖.

CE QUE ÇA N'EST PAS : une analyse modale. Un système à coefficients
périodiques n'a pas de fréquences propres uniques — l'exposant est défini
modulo 2π/T (repliement), ce que le module publie au lieu de le taire.
"""
import sys

import numpy as np

from vinkulum import Noyau


# ── mesure d'une perturbation d'état ─────────────────────────────────────────
def _log(m):
    """vecteur rotation de R ∈ SO(3), branche principale."""
    c = max(-1.0, min(1.0, (np.trace(m) - 1.0) / 2.0))
    a = np.arccos(c)
    w = 0.5 * np.array([m[2, 1] - m[1, 2], m[0, 2] - m[2, 0], m[1, 0] - m[0, 1]])
    return w if a < 1e-9 else (a / np.sin(a)) * w


def _expm(th):
    a = np.linalg.norm(th)
    if a < 1e-14:
        return np.eye(3)
    n = th / a
    k = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(a) * k + (1 - np.cos(a)) * (k @ k)


def _rots(etat):
    return [np.array(r, float).reshape(3, 3) for r in etat[2]]


def ecart(etat, ref):
    """Perturbation de `etat` par rapport à `ref`, en 12·n composantes.

    La rotation se mesure par log(R·R_refᵀ) — la perturbation À GAUCHE, celle
    que le noyau applique lui-même (R ← exp([δθ]×)·R). La mesurer dans une
    carte globale (log R − log R_ref) y ajouterait la distorsion de la carte,
    qui n'est pas petite dès que la référence tourne.
    """
    nb = len(etat[1])
    x = np.zeros(12 * nb)
    ra, rb = _rots(etat), _rots(ref)
    for i in range(nb):
        x[6 * i:6 * i + 3] = np.array(etat[1][i]) - np.array(ref[1][i])
        x[6 * i + 3:6 * i + 6] = _log(ra[i] @ rb[i].T)
        x[6 * nb + 6 * i:6 * nb + 6 * i + 3] = np.array(etat[3][i]) - np.array(ref[3][i])
        x[6 * nb + 6 * i + 3:6 * nb + 6 * i + 6] = np.array(etat[4][i]) - np.array(ref[4][i])
    return x


def _perturbe(n, d, eps):
    """Applique eps·d à l'état courant de `n` (d en 12·nb composantes)."""
    t, r, rot, v, w, vi = n.etat()
    nb = len(r)
    r = [list(np.array(r[i]) + eps * d[6 * i:6 * i + 3]) for i in range(nb)]
    rot = [list((_expm(eps * d[6 * i + 3:6 * i + 6]) @ np.array(rot[i], float).reshape(3, 3)).ravel())
           for i in range(nb)]
    v = [list(np.array(v[i]) + eps * d[6 * nb + 6 * i:6 * nb + 6 * i + 3]) for i in range(nb)]
    w = [list(np.array(w[i]) + eps * d[6 * nb + 6 * i + 3:6 * nb + 6 * i + 6]) for i in range(nb)]
    n.pose_etat(t, r, rot, v, w, list(vi))


def monodromie(construit, periode, h, base, eps=1e-6, tol_phi=1e-7, tol_per=1e-6,
               journal=False):
    """Matrice de transition sur une période, réduite à `base`.

    `construit()` rend un `Noyau` FRAIS posé sur la trajectoire de référence
    (on ne réutilise pas un noyau : l'inflow et les états d'éléments ne sont
    pas tous dans `etat()`).
    `base` : liste de directions admissibles, chacune un vecteur de 12·nb.
    Rend (Φ réduite m×m, diagnostics).
    """
    n = construit()
    t0 = n.t()
    ref0 = n.etat()
    n.simule(t0 + periode, h, tous=10 ** 9)
    refT = n.etat()

    # LA RÉFÉRENCE EST-ELLE PÉRIODIQUE ? sans ça, Φ ne veut rien dire — mais
    # l'écart mesuré porte AUSSI l'erreur d'intégration sur un tour, qui n'est
    # pas un défaut de périodicité. Un seuil absolu confond les deux. Le seul
    # contrôle qui les sépare est de RAFFINER : l'erreur de schéma tombe en h²,
    # une vraie dérive ne bouge pas. (Elle s'annule de toute façon dans les
    # différences centrées, qui la partagent ; ce contrôle protège du cas où
    # la trajectoire n'est pas périodique DU TOUT.)
    d_per = float(np.max(np.abs(ecart(refT, ref0))))
    if d_per > tol_per:
        n2 = construit()
        n2.simule(t0 + periode, h / 2, tous=10 ** 9)
        d2 = float(np.max(np.abs(ecart(n2.etat(), ref0))))
        if d2 > 0.4 * d_per:
            raise ValueError(f"la trajectoire de référence n'est pas périodique à {periode} s : "
                             f"écart {d_per:.3e} à h, {d2:.3e} à h/2 — il ne décroît pas comme "
                             f"l'erreur de schéma, c'est une vraie dérive")
        d_per = d2

    b = np.column_stack([np.asarray(d, float) / np.linalg.norm(d) for d in base])
    m = b.shape[1]
    # référence de contrainte, à comparer aux valeurs perturbées. ⚠ `phi_dot()`
    # rend G·u, SANS le terme ∂c/∂t d'une cible mobile : sur un arbre à
    # rotation imposée il vaut Ω et non zéro. Ce n'est donc pas une violation,
    # c'est la convention — d'où la comparaison à la RÉFÉRENCE et non à zéro.
    n_ref = construit()
    phi_ref = np.array(n_ref.phi() or [0.0])
    phid_ref = np.array(n_ref.phi_dot() or [0.0])
    nb2 = len(ref0[1])
    phi_viol = 0.0
    cols = []
    for k in range(m):
        out = []
        for sens in (+1.0, -1.0):
            nk = construit()
            _perturbe(nk, b[:, k], sens * eps)
            # LA DIRECTION EST-ELLE ADMISSIBLE ? une perturbation hors du noyau
            # des contraintes fait dépenser le premier pas à réparer Φ, et la
            # colonne mesure la projection au lieu de la dynamique.
            # ADMISSIBILITÉ, et il faut savoir ce qu'on mesure. Perturber la
            # POSITION change G(q), donc Φ̇ = G·u varie au PREMIER ordre en ε
            # même pour une direction parfaitement admissible (mesuré : 1e-4
            # pour ε = 1e-5, soit ‖u‖·ε). Ce n'est pas une faute, c'est que le
            # noyau de G tourne avec q. On contrôle donc Φ pour toute
            # direction, et Φ̇ seulement pour les directions purement en
            # VITESSE, où G ne bouge pas et où G·δu doit être nul.
            phi_viol = max(phi_viol,
                           float(np.max(np.abs(np.array(nk.phi() or [0.0]) - phi_ref))))
            if np.max(np.abs(b[:6 * nb2, k])) < 1e-12:
                phi_viol = max(phi_viol,
                               float(np.max(np.abs(np.array(nk.phi_dot() or [0.0]) - phid_ref))))
            nk.simule(t0 + periode, h, tous=10 ** 9)
            out.append(ecart(nk.etat(), refT))
        cols.append(b.T @ ((out[0] - out[1]) / (2.0 * eps)))
        if journal:
            print(f"    direction {k} : |δx(T)| {np.linalg.norm(cols[-1]):.4e}")
    # une direction admissible laisse Φ à O(ε²) : le seuil suit ε², pas ε
    if phi_viol > tol_phi * (eps / 1e-6) ** 2:
        raise ValueError(f"base non admissible : Φ atteint {phi_viol:.3e} après une "
                         f"perturbation de {eps:.1e} — c'est du premier ordre, donc la "
                         f"direction sort du noyau des contraintes")
    return np.column_stack(cols), dict(periodicite=d_per, phi=phi_viol, m=m)


def exposants(phi, periode, f_propre=None):
    """Table triée par |λ| décroissant.

    ⚠ CE QUI EST AMBIGU ET CE QUI NE L'EST PAS. La partie RÉELLE de l'exposant,
    σ = ln|λ|/T, est un fait : c'est le taux de croissance ou de décroissance,
    et le verdict de stabilité (|λ| ≤ 1) en découle sans hypothèse. La partie
    IMAGINAIRE ne l'est pas : elle n'est connue que modulo 2π/T, donc `f` est la
    fréquence REPLIÉE dans la bande [0, 1/2T] et rien de plus.

    Conséquence directe, et c'est un piège qu'on a payé : **ζ ne se calcule pas
    sur la fréquence repliée**. Sur la pale battante à ν_β = 1,12, la fréquence
    propre est 19,43 Hz, la repliée 2,08 — et ζ sorti de la seconde vaut 0,427
    au lieu de 0,050. Passer `f_propre` (venue d'ailleurs : une théorie, une
    analyse modale, un balayage en T) rend le ζ vrai ; sans elle, le champ
    `zeta_bande` est publié pour ce qu'il est, et `zeta` reste None.
    """
    lam = np.linalg.eigvals(phi)
    lam = lam[np.argsort(-np.abs(lam))]
    out = []
    for l in lam:
        s = np.log(complex(l)) / periode if abs(l) > 0 else complex(-np.inf, 0)
        f = abs(s.imag) / (2 * np.pi)
        w = 2 * np.pi * f_propre if f_propre else None
        out.append(dict(mult=l, module=abs(l), sigma=s.real, f=f,
                        zeta=(-s.real / np.hypot(s.real, w) if w else None),
                        zeta_bande=(-s.real / abs(s) if abs(s) > 0 else 0.0)))
    return out


# ── cas de contrôle ──────────────────────────────────────────────────────────
def _pendule_amorti(f_hz=2.0, zeta=0.15):
    """Oscillateur en rotation autour de x, à l'équilibre : σ exacts connus.

    L'équilibre est trivialement périodique, donc Floquet s'y applique et doit
    retomber sur `modes_complexes`, par un chemin qui ne partage avec lui NI la
    raideur tangente, NI la forme d'état, NI le solveur aux valeurs propres.
    """
    j = 0.01
    w0 = 2 * np.pi * f_hz
    k, c = j * w0 ** 2, 2 * zeta * j * w0

    def construit():
        n = Noyau([0.0, 0.0, 0.0])
        b = n.corps("b", 1.0, list(np.diag([j, j, j]).ravel()), [0.0] * 3)
        n.liaison("pivot", None, b, bloque_t=[0, 1, 2], bloque_r=[1, 2])
        n.couple("ressort", None, b, [1.0, 0.0, 0.0], ("ressort", [k, c, 0.0]))
        return n
    d_q, d_u = np.zeros(12), np.zeros(12)
    d_q[3] = 1.0            # rotation autour de x
    d_u[6 + 3] = 1.0        # vitesse angulaire autour de x
    return construit, [d_q, d_u], dict(f=f_hz * np.sqrt(1 - zeta ** 2), zeta=zeta)


def _pale_battante(omega=109.0, nu_beta=1.12, zeta=0.05):
    """UNE pale sur charnière de battement centrale, arbre à Ω imposé, sans aéro.

    C'est le cas que `modes_complexes` ne sait pas traiter — un corps CONTRAINT
    qui TOURNE — et il a une solution analytique : le battement obéit à
        I(β̈ + Ω²β) + C β̇ + K_β β = 0
    (le terme Ω²β est le raidissement centrifuge, exact pour une charnière
    centrale), donc ω_β = Ω·ν_β avec ν_β² = 1 + K_β/(IΩ²), et ζ = C/(2Iω_β).
    """
    i_b, l = 2.0, 1.0
    k_beta = i_b * omega ** 2 * (nu_beta ** 2 - 1.0)
    w_beta = omega * nu_beta
    c_beta = 2.0 * zeta * i_b * w_beta

    def construit():
        n = Noyau([0.0, 0.0, 0.0])
        moyeu = n.corps("moyeu", 1.0, list(np.diag([0.01] * 3).ravel()), [0.0] * 3,
                        w=[0.0, 0.0, omega])
        # pale le long de +x, battement autour de +y, masse ponctuelle à r telle
        # que I = m r² soit l'inertie de battement voulue
        m_p = i_b / l ** 2
        j_p = np.diag([1e-4, 1e-4, 1e-4])
        p = n.corps("pale", m_p, list(j_p.ravel()), [l, 0.0, 0.0],
                    w=[0.0, 0.0, omega], v=[0.0, omega * l, 0.0])
        n.liaison("arbre", None, moyeu, cible_r=([0.0, 0.0, 1.0], ("lineaire", [0.0, omega])))
        n.liaison("flap", moyeu, p, pa=[0.0] * 3, bloque_t=[0, 1, 2], bloque_r=[0, 2])
        n.couple("kbeta", moyeu, p, [0.0, 1.0, 0.0], ("ressort", [k_beta, c_beta, 0.0]))
        return n

    # direction admissible : le battement (rotation de la pale autour de y du
    # MOYEU, qui à t=0 est y du monde) et sa vitesse. La pale bat autour de la
    # charnière : δθ = e_y, et le point matériel à distance l monte de δθ×r.
    nb = 2
    d_q, d_u = np.zeros(12 * nb), np.zeros(12 * nb)
    # ⚠ LE SIGNE EST UNE CONTRAINTE, PAS UNE CONVENTION : une rotation δθ
    # autour de +y envoie le CdM (à +l·x̂) vers δθ·(ŷ × l x̂) = −δθ·l·ẑ. Le
    # prendre positif fait sortir la direction du noyau des contraintes, et le
    # contrôle d'admissibilité le dit au lieu de laisser passer une colonne
    # qui mesurerait une projection.
    d_q[6 + 4] = 1.0                     # δθ_y de la pale (corps 1)
    d_q[6 + 2] = -l                      # sa translation induite : δz = −δθ_y · l
    d_u[6 * nb + 6 + 4] = 1.0            # δω_y
    d_u[6 * nb + 6 + 2] = -l             # δv_z
    return construit, [d_q, d_u], dict(f=w_beta / (2 * np.pi), zeta=zeta,
                                       periode=2 * np.pi / omega, nu=nu_beta)


def _rotor_trime(mu=0.0, tours_etabli=6.0, h=None):
    """Le rotor de `vinkulum.rotor`, amené au régime établi, puis figé.

    La référence doit être périodique : on intègre plusieurs tours (l'inflow
    est un état explicite à retard, il lui faut le temps de s'établir), on
    relève l'état COMPLET, et `construit()` le repose tel quel.
    """
    from vinkulum.rotor import Rotor
    r = Rotor()
    th0 = r.theta_bet(3300.0)[0]
    h = h or r.t_tour / 300
    n0 = r.monte((th0, 0.0, 0.0), mu, t_end=(tours_etabli + 2) * r.t_tour, montee=1.0)
    n0.simule(tours_etabli * r.t_tour, h, tous=10 ** 9)
    fige = n0.etat()

    def construit():
        n = r.monte((th0, 0.0, 0.0), mu, t_end=(tours_etabli + 2) * r.t_tour, montee=1.0)
        n.pose_etat(fige[0], [list(x) for x in fige[1]], [list(x) for x in fige[2]],
                    [list(x) for x in fige[3]], [list(x) for x in fige[4]], list(fige[5]))
        return n

    # BASE CONSTRUITE SUR LA POSE RÉELLE, pas sur la géométrie nominale.
    # En avancement la pale BAT : à l'instant figé son axe de battement est
    # tourné de β et son CdM n'est plus sur l'envergure nominale. Une base
    # nominale sort alors du noyau des contraintes, et le contrôle
    # d'admissibilité le refuse — mesuré (3,8e-5 pour une perturbation de 1e-5,
    # soit un écart d'ordre 1). L'axe de battement se lit dans la rotation du
    # corps (colonne 1 = e_c en repère corps), et le déplacement induit du CdM
    # est δθ × r puisque la charnière est à l'origine.
    nb = 1 + r.n_pales
    base = []
    for i in range(r.n_pales):
        rot_i = np.array(fige[2][1 + i], float).reshape(3, 3)
        axe = rot_i[:, 1]                          # e_c du repère de la pale
        pos = np.array(fige[1][1 + i], float)
        for bloc in (0, 1):                        # q puis u
            d = np.zeros(12 * nb)
            k = 6 * nb * bloc + 6 * (1 + i)
            d[k:k + 3] = np.cross(axe, pos)
            d[k + 3:k + 6] = axe
            base.append(d)
    return r, construit, base, h


def demo(rapide=False):
    print("╔═ vinkulum — stabilité par la monodromie de Floquet")

    # 0. SENSIBILITÉ des exposants — contre une dérivée ANALYTIQUE
    #    σ = −ζω₀ ± i ω₀√(1−ζ²) avec ω₀ = 2π f, donc ∂σ/∂ζ et ∂σ/∂f sont
    #    connues exactement : c'est le seul contrôle qui puisse réfuter la
    #    formule de perturbation, et il ne partage rien avec elle.
    f0, z0 = 2.0, 0.15
    cons0, base0, ref0 = _pendule_amorti(f0, z0)
    tp = 0.3 / ref0["f"]
    ph0, _ = monodromie(cons0, tp, tp / 2000, base0)
    dph = dmonodromie(lambda z: _pendule_amorti(f0, z)[0](), tp, tp / 2000, base0, "zeta", z0)
    lam, ds = sensibilite(ph0, dph, tp)
    w0 = 2 * np.pi * f0
    ex = np.array([-w0 - 1j * s * w0 * z0 / np.sqrt(1 - z0 ** 2) for s in (1, -1)])
    # on apparie par la partie imaginaire de σ, uniquement pour COMPARER
    sig = np.log(lam.astype(complex)) / tp
    ordre = np.argsort(-sig.imag)
    e = max(abs(ds[ordre[i]] - ex[i]) / abs(ex[i]) for i in (0, 1))
    print(f"║ ∂σ/∂ζ mesuré {ds[ordre[0]]:.4f} contre {ex[0]:.4f} exact — écart {100 * e:.3f} %")
    assert e < 2e-3, (ds[ordre], ex)
    print("║   la formule δλ = Wᵢ·δΦ·xᵢ donne la dérivée DU MODE, identifié par son")
    print("║   vecteur propre : aucun appariement, donc rien à casser à un croisement")

    # 1. équilibre amorti : régression croisée contre modes_complexes
    construit, base, ref = _pendule_amorti()
    # ⚠ PAS la période propre : à T = 1/f_d la phase vaut exactement 2π et
    # l'exposant se replie sur zéro — Floquet rendrait « f = 0, ζ = 1 » sur un
    # oscillateur parfaitement sain. On sonde sur un tiers de période.
    t_per = 0.3 / ref["f"]
    phi, diag = monodromie(construit, t_per, t_per / 2000, base)
    tab = exposants(phi, t_per, f_propre=ref["f"])
    f_m = max(t["f"] for t in tab)
    z_m = np.mean([t["zeta"] for t in tab])
    n = construit()
    mc = n.modes_complexes(2)
    print(f"║ 1. équilibre amorti — Floquet f {f_m:.5f} Hz · ζ {z_m:.5f}")
    print(f"║      théorie          f {ref['f']:.5f} Hz · ζ {ref['zeta']:.5f}"
          f"   ·   modes_complexes f {mc[0][0]:.5f} · ζ {mc[0][1]:.5f}")
    assert abs(f_m - ref["f"]) / ref["f"] < 2e-3, (f_m, ref["f"])
    assert abs(z_m - ref["zeta"]) < 3e-3, (z_m, ref["zeta"])

    # 2. LE CAS QUE modes_complexes REFUSE : un corps contraint qui tourne
    construit, base, ref = _pale_battante()
    per = ref["periode"]
    phi, diag = monodromie(construit, per, per / 2000, base)
    f_th, z_th = ref["f"], ref["zeta"]
    tab = exposants(phi, per, f_propre=f_th)
    print(f"║ 2. pale battante Ω 109 rad/s · ν_β {ref['nu']} — référence périodique "
          f"à {diag['periodicite']:.2e}")
    for t in tab:
        print(f"║      λ {t['mult'].real:+.6f}{t['mult'].imag:+.6f}i  |λ| {t['module']:.6f}"
              f"   σ {t['sigma']:+.4f} 1/s   f repliée {t['f']:.3f} Hz   ζ {t['zeta']:+.5f}")
    # REPLIEMENT : ω_β = ν_β·Ω et la période EST celle du rotor, donc
    # f·T = ν_β = 1,12 > 1/2. L'exposant se voit à (ν_β − 1)·Ω, et c'est une
    # propriété de Floquet, pas une erreur — la fréquence propre ne se retrouve
    # que par un renseignement extérieur.
    f_att = abs(1.0 / per - f_th) if f_th * per > 0.5 else f_th
    sig_att = -z_th * 2 * np.pi * f_th
    f_mes = max(t["f"] for t in tab)
    s_mes = np.mean([t["sigma"] for t in tab])
    z_mes = np.mean([t["zeta"] for t in tab])
    print(f"║      théorie : f propre {f_th:.3f} Hz (ν_β·Ω) → repliée {f_att:.3f} Hz"
          f"   ·   σ {sig_att:+.4f} 1/s   ·   ζ {z_th:.5f}")
    # σ est le fait NON AMBIGU : c'est lui qu'on asserte en premier
    assert abs(s_mes - sig_att) / abs(sig_att) < 0.02, (s_mes, sig_att)
    assert abs(f_mes - f_att) / f_att < 0.02, (f_mes, f_att)
    assert abs(z_mes - z_th) / z_th < 0.03, (z_mes, z_th)
    assert all(t["module"] < 1.0 for t in tab), tab
    print(f"║      écarts : σ {100 * abs(s_mes / sig_att - 1):.2f} %   "
          f"f {100 * abs(f_mes / f_att - 1):.2f} %   ζ {100 * abs(z_mes / z_th - 1):.2f} %")
    print(f"║      STABLE : |λ| max {max(t['module'] for t in tab):.6f} < 1"
          f"   ·   `modes_complexes` REFUSE ce cas, Floquet le tranche")

    if rapide:
        print("╚═ floquet OK (rapide — le rotor complet est dans les bancs)")
        return

    # 3. LE ROTOR COMPLET, trimé, avec son aérodynamique : le produit
    import time
    for mu, nom in ((0.0, "stationnaire"), (0.15, "μ = 0,15")):
        t_dep = time.time()
        r, construit, base, h = _rotor_trime(mu)
        try:
            phi, diag = monodromie(construit, r.t_tour, h, base, eps=1e-5, tol_per=2e-4)
        except ValueError as e:
            print(f"║ 3. rotor {nom} : REFUSÉ — {e}")
            continue
        nu = r.nu_beta
        f_prop = nu * r.omega / (2 * np.pi)
        tab = exposants(phi, r.t_tour, f_propre=f_prop)
        # amortissement AÉRODYNAMIQUE du battement : le nombre de Lock donne
        # ζ ≈ γ/16 en première approximation (moyeu articulé, inflow figé) —
        # c'est un ORDRE DE GRANDEUR classique, pas une identité, et il sert
        # ici de recoupement et non de référence.
        z_lock = r.lock / 16.0
        z_mes = np.mean([t["zeta"] for t in tab[:4]])
        print(f"║ 3. rotor {nom} — {r.n_pales} pales, aéro et inflow actifs, "
              f"{diag['m']} directions, {time.time() - t_dep:.1f} s")
        print(f"║      |λ| de {min(t['module'] for t in tab):.4f} à "
              f"{max(t['module'] for t in tab):.4f}   ·   périodicité {diag['periodicite']:.1e}")
        for t in tab[:4]:
            print(f"║      σ {t['sigma']:+8.2f} 1/s   f repliée {t['f']:6.2f} Hz   "
                  f"ζ {t['zeta']:+.4f}")
        print(f"║      ζ moyen des quatre modes de battement {z_mes:.3f}   "
              f"contre γ/16 = {z_lock:.3f} (ordre de grandeur)")
        assert all(t["module"] < 1.0 for t in tab), tab
        assert 0.3 * z_lock < z_mes < 3.0 * z_lock, (z_mes, z_lock)

    print("╚═ floquet OK")



def sensibilite(phi, dphi, periode):
    """∂σ/∂p de CHAQUE mode, sans re-diagonaliser et SANS APPARIEMENT.

    C'est le point : une étude paramétrique de stabilité classique recalcule
    la monodromie à chaque valeur du paramètre, la rediagonalise, puis doit
    APPARIER les valeurs propres d'un calcul à l'autre — par tri en module, en
    fréquence, ou à l'œil. L'appariement est ambigu dès que deux modes se
    croisent ou se frôlent (*veering*), et c'est là que la stabilité se joue :
    la résonance sol/air d'un rotor EST une coalescence de modes.

    La perturbation de valeur propre l'évite. Pour Φx = λx et W = X⁻¹ (dont
    les LIGNES sont les vecteurs propres gauches, normalisés à Wᵢ·xᵢ = 1) :

        δλᵢ = Wᵢ · δΦ · xᵢ        puis      δσᵢ = δλᵢ / (λᵢ · T)

    puisque σ = ln λ / T. Le mode est identifié par SON vecteur propre, pas
    par son rang dans une liste triée : la dérivée suit le mode à travers un
    croisement.

    `dphi` est ∂Φ/∂p (même forme que `phi`). Rend (λ, ∂σ/∂p), l'ordre étant
    celui de `numpy.linalg.eig`, c'est-à-dire celui des colonnes de X.
    """
    lam, x = np.linalg.eig(np.asarray(phi, float))
    w = np.linalg.inv(x)
    d = np.asarray(dphi, float)
    dlam = np.array([w[i] @ d @ x[:, i] for i in range(len(lam))])
    return lam, dlam / (lam * periode)


def dmonodromie(construit, periode, h, base, param, p0, rel=1e-4, **kw):
    """∂Φ/∂p par différence centrée sur DEUX monodromies.

    `construit` prend le paramètre ; `param` n'est là que pour la lisibilité
    des messages. Deux monodromies suffisent quel que soit le nombre de modes
    — c'est `sensibilite` qui les distribue ensuite sur tous les modes à la
    fois, alors qu'une DF sur les exposants en demanderait deux PAR mode et
    poserait le problème d'appariement à chaque fois.
    """
    d = abs(p0) * rel if p0 else rel
    a, _ = monodromie(lambda: construit(p0 + d), periode, h, base, **kw)
    b, _ = monodromie(lambda: construit(p0 - d), periode, h, base, **kw)
    return (np.asarray(a, float) - np.asarray(b, float)) / (2.0 * d)

if __name__ == "__main__":
    demo(rapide="rapide" in sys.argv)
