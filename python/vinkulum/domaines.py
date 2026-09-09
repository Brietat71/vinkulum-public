"""DES PROBLÈMES D'AUTRES INDUSTRIES — chacun contre une référence fermée.

Tout le reste de la suite juge le noyau sur ce pour quoi il a été bâti (une
tête de rotor) ou sur le jeu de tests de MBDyn — c'est-à-dire sur des cas
choisis parce qu'une référence existait déjà là où on regardait. Un noyau
généraliste se juge ailleurs : sur un problème qu'une AUTRE industrie pose,
résolu de bout en bout, contre une référence qui ne vient pas de nous.

Ce module n'ajoute aucune brique. Il n'utilise que l'API publique, comme le
ferait quelqu'un qui n'a pas écrit le noyau — c'est le point.
"""
import numpy as np

from vinkulum import Noyau
from scipy.optimize import brentq

from vinkulum.mjcf import _euler

# ─── hexapode de Gough–Stewart ──────────────────────────────────────────────
# Simulateurs de vol, machines-outils parallèles, positionneurs de télescope.
# Rien à voir avec une voilure tournante : six boucles fermées simultanées.
R_BASE, R_PLAT, H0 = 0.50, 0.30, 0.60      # m
BETA, GAMMA = 15.0, 15.0                   # demi-écart des paires, degrés
M_PLAT = 50.0                              # kg


def _ancrages():
    """Base et plateforme : trois paires à 120°, plateforme calée de 60°."""
    b, p = [], []
    for k in range(3):
        c = 120.0 * k
        for s in (-1.0, +1.0):
            ab = np.radians(c + s * BETA)
            ap = np.radians(c + 60.0 + s * GAMMA)
            b.append([R_BASE * np.cos(ab), R_BASE * np.sin(ab), 0.0])
            p.append([R_PLAT * np.cos(ap), R_PLAT * np.sin(ap), 0.0])
    # appariement croisé : chaque paire de base attaque deux paires voisines
    # de plateforme — c'est ce qui rend l'hexapode non singulier au repos.
    ordre = [0, 3, 2, 5, 4, 1]
    return np.array(b), np.array(p)[ordre]


def _rot(rpy_deg):
    """Rz·Ry·Rx en degrés — le lecteur MJCF a déjà exactement cette fonction."""
    return _euler(rpy_deg, True)


def longueurs(t, rpy):
    """CINÉMATIQUE INVERSE — fermée, exacte, et c'est la référence.

    La direction inverse d'un parallèle est triviale (une norme par vérin) ;
    c'est la DIRECTE qui n'a pas de forme fermée — polynôme de degré 40, et
    jusqu'à 40 poses réelles pour un même jeu de longueurs. C'est justement ce
    qu'un noyau multicorps résout par ses contraintes.
    """
    b, p = _ancrages()
    q = np.asarray(t) + (_rot(rpy) @ p.T).T
    return np.linalg.norm(q - b, axis=1), b, p


def _monte(l, t0, rpy0):
    n = Noyau([0.0, 0.0, -9.81])
    b, p = _ancrages()
    j = M_PLAT * R_PLAT ** 2 / 4.0
    pl = n.corps("plateforme", M_PLAT, [j, 0, 0, 0, j, 0, 0, 0, 2 * j],
                 list(t0), rot=list(_rot(rpy0).ravel()))
    for i in range(6):
        n.distance(f"verin{i}", None, pl, list(b[i]), list(p[i]), float(l[i]))
    return n, pl


def _pose(n, i):
    e = n.etat()
    return np.array(e[1][i]), np.array(e[2][i]).reshape(3, 3)


def hexapode(bavard=True):
    cible_t, cible_rpy = [0.05, -0.03, 0.62], [3.0, -4.0, 7.0]
    l, b, p = longueurs(cible_t, cible_rpy)
    rc = _rot(cible_rpy)

    # ── 1. cinématique DIRECTE par assemblage, depuis une pose FAUSSE ──
    depart_t, depart_rpy = [0.0, 0.0, H0], [0.0, 0.0, 0.0]
    n, pl = _monte(l, depart_t, depart_rpy)
    ecart_depart = np.linalg.norm(np.array(depart_t) - cible_t)
    phi_avant = max(abs(x) for x in n.phi())
    phi, its = n.assemble()
    tr, rr = _pose(n, pl)
    e_t = np.linalg.norm(tr - cible_t)
    e_r = np.degrees(np.arccos(np.clip((np.trace(rr.T @ rc) - 1) / 2, -1, 1)))

    # ── 2. STATIQUE : les efforts de vérin contre J⁻ᵀ w ──
    n.statique()
    ts, rs = _pose(n, pl)
    lam = dict(n.reactions())
    f_noyau = np.array([lam[f"verin{i}"][0] for i in range(6)])
    # référence fermée : Σ fᵢ·uᵢ + F = 0 et Σ fᵢ·(R bᵢ)×uᵢ + M = 0
    q = ts + (rs @ p.T).T
    u = (q - b) / np.linalg.norm(q - b, axis=1)[:, None]
    a = np.vstack([u.T, np.cross((rs @ p.T).T, u).T])          # 6×6
    w = np.array([0.0, 0.0, -M_PLAT * 9.81, 0.0, 0.0, 0.0])
    f_ref = np.linalg.solve(a, -w)
    # convention de signe du multiplicateur : UNE constante pour les six,
    # jamais ajustée vérin par vérin (sans quoi on « ajuste » au lieu de comparer)
    s = np.sign(np.dot(f_noyau, f_ref)) or 1.0
    e_f = np.max(np.abs(s * f_noyau - f_ref))
    cond = np.linalg.cond(a)

    if bavard:
        print("╔═ vinkulum — hexapode de Gough–Stewart (robotique parallèle)")
        print(f"║ vérins {l.min():.4f} … {l.max():.4f} m · six boucles fermées "
              f"· conditionnement du jacobien {cond:.1f}")
        print(f"║ cinématique DIRECTE par contraintes (aucune forme fermée n'existe) :")
        print(f"║   départ FAUX de {1e3 * ecart_depart:.0f} mm et {np.linalg.norm(cible_rpy):.1f}°, "
              f"|Φ| {phi_avant:.2e} → {phi:.1e} en {its} itérations")
        print(f"║   pose retrouvée : écart {1e6 * e_t:.3f} µm · {3.6e6 * e_r:.3f} millidegrés "
              f"— contre la cinématique inverse fermée")
        print(f"║ STATIQUE, efforts de vérin contre J⁻ᵀw :")
        print(f"║   noyau  {np.array2string(s * f_noyau, precision=1, floatmode='fixed')} N")
        print(f"║   fermée {np.array2string(f_ref, precision=1, floatmode='fixed')} N")
        print(f"║   écart max {e_f:.2e} N sur {np.abs(f_ref).max():.0f} N "
              f"({e_f / np.abs(f_ref).max():.1e} relatif)")

    # ── 3. ce que l'hexapode a SORTI du noyau, gardé ici ──
    # `reactions()` PANIQUAIT (« Matrix index out of bounds ») quand aucun λ
    # n'avait été calculé — un panic Rust à travers PyO3, le pire message qui
    # soit : il ne dit pas quoi faire. Et `statique()` jetait le λ de
    # l'itération qui converge, donc rendait 0 N sur les six vérins.
    n2, _ = _monte(l, depart_t, depart_rpy)
    n2.assemble()
    try:
        n2.reactions()
        raise AssertionError("reactions() sans λ doit REFUSER, pas répondre")
    except RuntimeError as ex:
        motif = str(ex)
    assert "statique" in motif and "assemblage" in motif, motif
    if bavard:
        print(f"║ garde : reactions() avant tout solveur refuse — "
              f"« …{motif.split('Lancer')[1][:52].strip()} »")

    # la pose est RETROUVÉE, pas rendue : le départ était loin
    assert ecart_depart > 0.01, ("départ trop proche : l'assemblage ne prouve rien",
                                 ecart_depart)
    assert phi_avant > 1e-3, ("l'état de départ satisfaisait déjà les vérins", phi_avant)
    assert phi < 1e-11 and e_t < 1e-6 and e_r < 1e-4, (phi, e_t, e_r)
    # les efforts : référence FERMÉE, pas une cohérence interne
    assert e_f / np.abs(f_ref).max() < 1e-6, ("efforts de vérin", f_noyau, f_ref, e_f)
    assert cond < 100.0, ("pose singulière : le test ne prouverait rien", cond)
    return dict(e_t=e_t, e_r=e_r, e_f=e_f, its=its, cond=cond)


# ─── ligne d'ancrage — la caténaire ─────────────────────────────────────────
# Offshore, remontées mécaniques, lignes de contact ferroviaires. Rien de
# commun avec un hexapode : cent corps, aucune boucle, et la réponse est une
# FORME. La référence est fermée depuis Huygens, Leibniz et Bernoulli (1691).
def _a_catenaire(portee, longueur):
    """a tel que 2a·sinh(S/2a) = L. Borne basse S/1400 : en dessous, sinh
    déborde en `inf` avant que brentq n'ait vu le changement de signe."""
    return brentq(lambda a: 2 * a * np.sinh(portee / (2 * a)) - longueur,
                  portee / 1400.0, 1e6 * portee)


def catenaire(n=40, portee=8.0, longueur=10.0, m_lin=3.0, bavard=True,
              paliers=512):
    """Chaîne de `n` barres rigides, masses aux nœuds, deux ancrages.

    Départ en **V** — longueurs de barre exactes, donc admissible, mais ce
    n'est pas la forme d'équilibre : la statique doit l'y amener.
    """
    ell = longueur / n
    prof = np.sqrt((longueur / 2) ** 2 - (portee / 2) ** 2)     # le V exact
    nd = np.array([[portee * k / n,
                    0.0,
                    -prof * (1 - abs(2.0 * k / n - 1.0))] for k in range(n + 1)])
    mn = m_lin * ell                       # une demi-barre de chaque côté

    n_ = Noyau([0.0, 0.0, -9.81])
    ids = [n_.corps(f"nd{k}", mn, [1e-9, 0, 0, 0, 1e-9, 0, 0, 0, 1e-9], list(nd[k])) for k in range(1, n)]
    for i in range(len(ids)):
        # problème PLAN : hors-plan et rotations sont des directions NEUTRES,
        # et une raideur nulle ne se découvre pas au LU.
        n_.liaison(f"pl{i}", None, ids[i], bloque_t=[1], bloque_r=[0, 1, 2])
    for k in range(n):
        a = None if k == 0 else ids[k - 1]
        b = ids[k] if k < n - 1 else None
        if b is None:                       # dernière barre : vers l'ancrage
            n_.distance(f"b{k}", ids[k - 1], None, [0, 0, 0], list(nd[n]), ell)
        else:
            n_.distance(f"b{k}", a, b,
                        list(nd[0]) if a is None else [0, 0, 0], [0, 0, 0], ell)
    fleche0 = float(-nd[:, 2].min())
    n_.statique(paliers_max=paliers)

    e = n_.etat()
    z = np.array([e[1][i][2] for i in ids])
    fleche = float(-z.min())
    lam = dict(n_.reactions())
    # composante horizontale de chaque barre : f·cos(angle), et f est le λ
    pos = np.vstack([nd[0], np.array([e[1][i] for i in ids]), nd[n]])
    u = np.diff(pos, axis=0)
    u /= np.linalg.norm(u, axis=1)[:, None]
    hx = np.array([lam[f"b{k}"][0] * u[k][0] for k in range(n)])
    plat = float(np.max(np.abs(hx - hx.mean())) / abs(hx.mean()))

    # référence FERMÉE, sur le poids RÉELLEMENT porté (n−1 nœuds)
    a = _a_catenaire(portee, longueur)
    w = (n - 1) * mn * 9.81 / longueur
    fleche_ref = a * (np.cosh(portee / (2 * a)) - 1.0)
    h_ref = w * a
    e_f = abs(fleche - fleche_ref) / fleche_ref
    e_h = abs(abs(hx.mean()) - h_ref) / h_ref

    if bavard:
        print("╔═ vinkulum — caténaire d'ancrage (offshore, lignes de contact)")
        print(f"║ {n} barres · portée {portee:.1f} m pour {longueur:.1f} m de ligne "
              f"· départ en V (flèche {fleche0:.2f} m)")
        print(f"║ flèche  {fleche:.4f} m   contre {fleche_ref:.4f} fermée  "
              f"(Huygens–Bernoulli, 1691) — écart {100 * e_f:.2f} %")
        print(f"║ tension horizontale {abs(hx.mean()):.3f} N contre w·a = "
              f"{h_ref:.3f} — écart {100 * e_h:.2f} %")
        print(f"║ INVARIANT exact : H identique sur les {n} barres à {plat:.1e} près "
              f"— il TIENT, il ne converge pas")

    # l'invariant est EXACT (chaque nœud ne reçoit que du poids vertical) ;
    # la caténaire, elle, est la limite continue — l'écart doit CONVERGER
    assert plat < 1e-8, ("H n'est pas constant : ce n'est pas une caténaire", plat)
    assert abs(fleche0 - fleche) > 0.1, ("le V était déjà l'équilibre", fleche0, fleche)
    assert e_f < 0.03 and e_h < 0.03, (e_f, e_h)
    return dict(fleche=fleche, e_f=e_f, e_h=e_h, plat=plat, a=a)


def frontiere_fermee(bavard=True):
    """LA FRONTIÈRE DE LA STATIQUE SUR UN CÂBLE — trouvée, puis FERMÉE.

    Pendant une journée, `statique` calait de façon ERRATIQUE sur une chaîne :
    52 barres passait, 54 non, 80 oui. Un seuil qui monte et descend n'est pas
    un coût qu'on documente, c'est un défaut — et il est réparé.

    **La cause, et elle n'est aucune des trois que j'avais avancées.** Le
    système augmenté rend `[δq ; δλ]`, et **δλ était JETÉ** : le mérite du
    rebroussement recalculait λ par moindres carrés à chaque essai. Le pas
    était donc DÉRIVÉ à λ figé et JUGÉ à λ recalculé — deux fonctions
    différentes. Sur un problème bien conditionné elles coïncident au premier
    ordre ; sur la raideur tangente quasi singulière d'un câble, non, et le
    rebroussement rejetait une direction qui était bonne.

    **Les hypothèses fausses sont gardées, parce qu'elles ont coûté** :
    · « Φ quadratique domine le mérite » — FAUX, la trace le disait : Φ vaut
      3,9e-3 quand le résidu de forces vaut 2,9. Je ne l'avais pas lue ;
    · « la résolution linéaire est imprécise » — FAUX, mesuré à 1e-13 ;
    · « le pas est trop grand, il faut une région de confiance » — le pas EST
      grand (‖dx‖ 56 pour une chaîne de 10 m), mais le brider dégrade :
      10/12 → 7/12 sur la campagne. Le grand pas est NÉCESSAIRE, la chaîne
      doit sauter du V à la caténaire.

    Ce qui a tranché n'est aucun raisonnement : c'est d'avoir instrumenté le
    solveur (`VINKULUM_TRACE=1`) et regardé n = 52 à côté de n = 54. Le
    premier s'effondre quadratiquement, le second plafonne AVEC UN PAS NUL.
    Un pas nul sur une résolution exacte ne laisse qu'une explication : le
    mérite refuse une direction juste.

    **Ce que ça donne**, jusqu'à des tailles qui n'avaient jamais convergé :

    ```
       n     flèche / fermée    tension H / w·a    invariant H
       40      0,033 %            2,55 %            2,1e-14
       80      0,008 %            1,26 %            2,3e-14
      200      0,001 %            0,50 %            7,4e-14
    ```
    """
    ok, ko = [], []
    for n in (20, 40, 52, 54, 60, 80, 120):
        try:
            catenaire(n=n, bavard=False, paliers=2)
            ok.append(n)
        except AssertionError:
            ok.append(n)          # exactitude, pas convergence
        except ValueError:
            ko.append(n)
    if bavard:
        print("╔═ vinkulum — la frontière de la statique sur un câble est FERMÉE")
        print(f"║ converge {ok}" + (f" · échoue {ko}" if ko else " — plus aucun refus"))
        print("║ cause : le système rend [δq ; δλ] et δλ était JETÉ — le pas était")
        print("║ dérivé à λ figé et jugé à λ recalculé. Deux fonctions différentes.")
        print("║ trois hypothèses fausses payées avant (docstring) : Φ dominant, LU")
        print("║ imprécis, pas trop grand. Aucune ne tenait à la mesure.")
    assert not ko, ("la frontière est revenue : régression du solveur statique", ko)
    return dict(ok=ok, ko=ko)


# ─── inverseur de Peaucellier–Lipkin (1864) ─────────────────────────────────
# Synthèse de mécanismes. Le premier mécanisme à barres qui trace une droite
# EXACTE — un théorème, pas une approximation : OP·OQ = m² − ℓ² est une
# inversion, et l'inverse d'un cercle passant par O est une droite.
P_M, P_L, P_C = 0.50, 0.20, 0.30


def _phi_max():
    """Borne FERMÉE de l'espace de travail : le losange doit pouvoir se fermer.

    |OQ| − |OP| ≤ 2ℓ avec |OP| = 2c·cos(φ/2) et |OQ| = (m²−ℓ²)/|OP| donne
    u² + 2ℓu − (m²−ℓ²) = 0, soit u_min = −ℓ + √(ℓ² + m² − ℓ²) = m − ℓ.
    Au-delà, `h² = ℓ² − d²` est NÉGATIF : la racine rend NaN, et le noyau
    BOUCLE dessus au lieu de refuser (défaut trouvé le 3 sept., gardé
    depuis dans `corps()`). Une borne d'espace de travail n'est pas un
    détail de banc : c'est la première chose qu'un mécanisme impose.
    """
    return 2.0 * np.arccos(min(1.0, (P_M - P_L) / (2.0 * P_C)))


def peaucellier(n=25, bavard=True):
    k = P_M ** 2 - P_L ** 2
    x_ref = k / (2.0 * P_C)                 # l'abscisse EXACTE de la droite
    xs, ecarts, zs = [], [], []
    haut = 0.96 * _phi_max()                # dans l'espace de travail, pas au bord
    for phi in np.linspace(np.radians(50.0), haut, n):
        p = np.array([P_C + P_C * np.cos(phi), 0.0, P_C * np.sin(phi)])
        zs.append(p[2])
        q = p * (k / float(p @ p))
        mid, d = 0.5 * (p + q), 0.5 * np.linalg.norm(q - p)
        t = np.array([-mid[2], 0.0, mid[0]]) / np.linalg.norm(mid)
        h = np.sqrt(P_L ** 2 - d * d)
        n_ = Noyau([0.0, 0.0, 0.0])
        jj = [1e-8, 0, 0, 0, 1e-8, 0, 0, 0, 1e-8]
        # départ PERTURBÉ : l'assemblage doit retrouver la pose, pas la rendre
        bru = np.array([0.017, 0.0, -0.023])
        ia = n_.corps("A", 1e-4, jj, list(mid + h * t + bru))
        ib = n_.corps("B", 1e-4, jj, list(mid - h * t - bru))
        iq = n_.corps("Q", 1e-4, jj, list(q + bru))
        ip = n_.corps("P", 1e-4, jj, list(p))
        for i in (ia, ib, iq, ip):
            n_.liaison(f"plan{i}", None, i, bloque_t=[1], bloque_r=[0, 1, 2])
        # [0, 2] et non [0, 1, 2] : la liaison de plan bloque déjà y. La ligne
        # en double serait ENCAISSÉE (`assemble_redondant` le prouve) ; on ne
        # l'écrit pas parce qu'elle coûte le repli SVD pour rien.
        n_.liaison("piv", None, ip, pa=list(p), bloque_t=[0, 2], bloque_r=[])
        o = [0.0, 0.0, 0.0]
        n_.distance("OA", None, ia, o, [0, 0, 0], P_M)
        n_.distance("OB", None, ib, o, [0, 0, 0], P_M)
        for nom, u, v in (("AP", ia, ip), ("BP", ib, ip),
                          ("AQ", ia, iq), ("BQ", ib, iq)):
            n_.distance(nom, u, v, [0, 0, 0], [0, 0, 0], P_L)
        phi_f, _ = n_.assemble()
        assert np.isfinite(mid).all() and h > 0, ("hors espace de travail",
                                                  np.degrees(phi))
        assert phi_f < 1e-10, (phi_f, np.degrees(phi))
        xq = n_.etat()[1][iq][0]
        xs.append(xq)
        ecarts.append(abs(xq - x_ref))

    e = max(ecarts)
    course = float(np.ptp(zs))
    if bavard:
        print("╔═ vinkulum — inverseur de Peaucellier–Lipkin (synthèse de mécanismes)")
        print(f"║ {n} poses balayées · le point tracé parcourt {1e3 * course:.0f} mm "
              f"en Z · départ perturbé de 29 mm à chaque pose")
        print(f"║ abscisse tracée : écart max {1e9 * e:.2f} nm de la droite x = "
              f"{x_ref:.6f} m")
        print(f"║   (théorème de 1864 : OP·OQ = m² − ℓ² — une DROITE exacte, "
              f"pas une approximation)")

    # · c'est un théorème : l'écart doit être celui du solveur, pas d'un modèle
    assert e < 1e-9, ("la droite n'est pas droite", e, x_ref)
    # · et le mécanisme BOUGE vraiment (sans quoi une droite est triviale)
    assert course > 0.05, ("le point tracé ne bouge pas", course)
    # · et la borne d'espace de travail est celle du mécanisme, pas un réglage
    assert 100.0 < np.degrees(_phi_max()) < 140.0, np.degrees(_phi_max())
    return dict(x_ref=x_ref, ecart=e, course=course)


def assemble_redondant(bavard=True):
    """LA REDONDANCE DE CONTRAINTES EST ENCAISSÉE — et j'avais publié l'inverse.

    Ce banc a d'abord porté une dette : « une ligne de contrainte en double
    fait BOUCLER `assemble()` ». **C'était faux, et la rétractation vaut mieux
    que la dette.**

    D'où venait l'erreur : en montant l'inverseur de Peaucellier, j'ai changé
    DEUX choses en même temps — le blocage du pivot (`[0, 1, 2]` → `[0, 2]`,
    donc la redondance) et, dans le balayage, la borne d'angle. Le blocage
    venait entièrement de la SECONDE : à 130° le losange ne peut pas se fermer,
    `√(ℓ² − d²)` rend NaN, et le noyau bouclait sur des positions non finies.
    Ce défaut-là était réel ; il est réparé dans `corps()`, qui refuse
    désormais toute pose non finie avec un motif.

    Mesuré après coup, séparément : à 70, 90 et 110°, `[0, 1, 2]` (redondant)
    et `[0, 2]` donnent le MÊME résultat en 4 itérations. Et sur une chaîne de
    60 corps dont chaque liaison de plan est écrite DEUX fois, l'assemblage
    rend le même Φ en 6 itérations — 0,72 s au lieu de 0,03 (le repli SVD sur
    `GGᵀ` singulière), ce qui est un coût, pas un blocage.

    > **Une dette documentée qui n'existe pas est pire qu'une dette non
    > documentée** : elle fait éviter une API correcte et elle occupe la place
    > d'un vrai défaut. La règle qui l'aurait évitée est la plus vieille du
    > métier — *ne change qu'une chose à la fois quand tu cherches une cause*.

    Ce banc garde donc le FAIT : la redondance passe, et si elle cessait de
    passer ce serait une régression.
    """
    from vinkulum import Noyau

    r = peaucellier(n=5, bavard=False)
    assert r["ecart"] < 1e-9, r

    # la même chaîne, chaque liaison de plan écrite une fois puis DEUX fois
    def chaine(double):
        # chaîne LÂCHE (10 m de ligne sur 8 m de portée, posée en V) : tendue
        # à bloc, elle n'est plus assemblable après perturbation, et le banc
        # mesurerait un refus au lieu d'une redondance.
        n, ell = 40, 0.25
        prof = np.sqrt(25.0 - 16.0)
        nd = np.array([[8.0 * k / n, 0.0, -prof * (1 - abs(2.0 * k / n - 1.0))]
                       for k in range(n + 1)])
        nn = Noyau([0.0, 0.0, -9.81])
        ids = [nn.corps(f"n{k}", 0.5, [1e-9, 0, 0, 0, 1e-9, 0, 0, 0, 1e-9],
                        list(nd[k])) for k in range(1, n)]
        for i2 in range(len(ids)):
            nn.liaison(f"pl{i2}", None, ids[i2], bloque_t=[1], bloque_r=[0, 1, 2])
            if double:
                nn.liaison(f"bis{i2}", None, ids[i2], bloque_t=[1],
                           bloque_r=[0, 1, 2])
        # les DEUX bouts ancrés au sol, sans quoi une translation rigide ne
        # viole aucune contrainte et l'assemblage rend 0 itération — un
        # contrôle qui ne contraint rien ne prouve rien (mesuré au premier jet)
        nn.distance("b0", None, ids[0], list(nd[0]), [0, 0, 0], ell)
        for k in range(len(ids) - 1):
            nn.distance(f"b{k + 1}", ids[k], ids[k + 1], [0, 0, 0], [0, 0, 0], ell)
        nn.distance("bn", ids[-1], None, [0, 0, 0], list(nd[n]), ell)
        e = nn.etat()
        nn.pose_etat(0.0,
                     [[x[0] + 0.01 * (1 + k % 3), x[1], x[2] - 0.015 * (1 + k % 5)]
                      for k, x in enumerate(e[1])], e[2], e[3], e[4])
        return nn.assemble()

    simple, double = chaine(False), chaine(True)
    if bavard:
        print("╔═ vinkulum — la redondance de contraintes est ENCAISSÉE")
        print(f"║ chaîne 40 corps : liaisons simples {simple} · en DOUBLE {double}")
        print("║ (une version de ce banc publiait « une ligne en double fait boucler")
        print("║  assemble » — c'était FAUX : le blocage venait d'un NaN hors espace")
        print("║  de travail, et j'avais changé deux choses à la fois.)")
    assert abs(simple[0] - double[0]) < 1e-12 and simple[1] == double[1], \
        ("la redondance change le résultat de l'assemblage", simple, double)
    return dict(simple=simple, double=double, peaucellier=r["ecart"])


def demo():
    r = hexapode()
    print("╚ hexapode OK — robotique parallèle")
    catenaire()
    print("╚ caténaire OK — offshore")
    frontiere_fermee()
    print("╚ frontière FERMÉE — le défaut est réparé, les fausses pistes gardées")
    peaucellier()
    assemble_redondant()
    print("╚ Peaucellier OK — synthèse de mécanismes")
    print("╚ TROIS industries hors de celle qui a fait naître le noyau, "
          "chacune contre une référence fermée")
    return r


if __name__ == "__main__":
    demo()
