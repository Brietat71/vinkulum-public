"""vinkulum.contact — le contact par pénalité régularisée, et ce qu'il vaut.

    python -m vinkulum.contact        # enfoncement de Hertz, puis restitution

Le manque le plus sérieux du noyau face à un code de contact (RecurDyn en a
fait son métier). Voici la brique : une sphère contre un demi-espace, loi de
Hertz avec amortissement de **Hunt–Crossley** et frottement de Coulomb lissé.

    F_n = k δ^e (1 + c δ̇)        F_t = −μ F_n tanh(‖v_t‖/v_ε) · v_t/‖v_t‖

POURQUOI LA PÉNALITÉ, ET PAS UN MULTIPLICATEUR. Les conditions d'unilatéralité
(δ ≥ 0, F ≥ 0, δ·F = 0) rendent le résidu NON DÉRIVABLE à l'instant du choc, ce
qui effondre le pas de tout intégrateur implicite. La pénalité régularisée
garde tout C¹ — et donc, ici, différentiable par l'AD du noyau, qui fournit la
tangente exacte du contact comme celle de n'importe quel autre élément.

POURQUOI HUNT–CROSSLEY, ET PAS UN AMORTISSEUR LINÉAIRE. Avec c·δ̇, la force
d'amortissement est NON NULLE à δ = 0 : elle saute à l'impact et TIRE au
décollement, ce qui est non physique et fabrique une discontinuité. En la
rendant proportionnelle à δ^e, elle part de zéro et y revient. C'est la même
fin que la « step function » d'ADAMS et de RecurDyn, obtenue sans paramètre
de seuil supplémentaire.
"""
import sys

import numpy as np

from vinkulum import Noyau

G = 9.81


def bille(k=1e6, expo=1.5, c=0.0, mu=0.0, m=0.5, R=0.02, z0=None, vz=0.0, nonlisse=False, e=0.0):
    N = Noyau([0.0, 0.0, -G])
    b = N.corps("bille", m, list((2.0 / 5.0 * m * R * R * np.eye(3)).ravel()),
                [0.0, 0.0, R if z0 is None else z0], v=[0.0, 0.0, vz])
    N.contact("sol", b, [0, 0, 0], R, [0, 0, 1.0], [0, 0, 0], k, expo=expo, c=c, mu=mu,
              nonlisse=nonlisse, restitution=e)
    return N, b, dict(m=m, R=R, k=k)


def restitution(c, v0=0.5, k=1e6, m=0.5, R=0.02, h=2e-6):
    """Rebond : vitesse après / vitesse avant, sans gravité (choc pur)."""
    N = Noyau([0.0, 0.0, 0.0])
    b = N.corps("bille", m, list((2.0 / 5.0 * m * R * R * np.eye(3)).ravel()),
                [0.0, 0.0, R], v=[0.0, 0.0, -v0])
    N.contact("sol", b, [0, 0, 0], R, [0, 0, 1.0], [0, 0, 0], k, expo=1.5, c=c)
    # durée : le contact hertzien dure ~2,94·(m/(k√R))^0,4 / v0^0,2 — on prend
    # large et on s'arrête quand la bille est repartie
    t, dt = 0.0, 2e-4
    while t < 0.05:
        N.simule(t + dt, h, tous=10 ** 9)
        t += dt
        if N.etat()[3][0][2] > 0 and N.contacts()[0][1] <= 0:
            break
    return N.etat()[3][0][2] / v0


def choc(m1, m2, v1, c=0.0, R=0.02, k=1e6, h=2e-6, t_max=0.2):
    """Choc central de deux billes libres, sans pesanteur.

    Le contrôle le plus exigeant du contact à DEUX corps, et il ne demande
    aucune donnée : à c = 0 le choc est élastique, donc la quantité de
    mouvement ET l'énergie se conservent, et les vitesses sortantes sont
    celles de la théorie — v₁' = (m₁−m₂)/(m₁+m₂)·v₁, v₂' = 2m₁/(m₁+m₂)·v₁.
    """
    N = Noyau([0.0, 0.0, 0.0])
    jj = lambda m: list((2.0 / 5.0 * m * R * R * np.eye(3)).ravel())
    a = N.corps("a", m1, jj(m1), [-0.05, 0, 0], v=[v1, 0, 0])
    b = N.corps("b", m2, jj(m2), [0.05, 0, 0])
    N.contact("cc", a, [0, 0, 0], R, k=k, c=c, b=b, pb=[0, 0, 0], rayon_b=R)
    p0, e0 = m1 * v1, 0.5 * m1 * v1 ** 2
    t = 0.0
    while t < t_max:
        N.simule(t + 2e-3, h, tous=10 ** 9)
        t += 2e-3
        if N.contacts()[0][1] <= 0 and N.etat()[3][1][0] > 1e-9 and t > 0.02:
            break
    va, vb = N.etat()[3][0][0], N.etat()[3][1][0]
    p = m1 * va + m2 * vb
    e = 0.5 * m1 * va ** 2 + 0.5 * m2 * vb ** 2
    return va, vb, p / p0 - 1.0, e / e0 - 1.0


def rebond(h_libre, h_contact=None, marge=0.0, t_end=1.2, k=1e10, c=0.6,
           R=0.02, m=0.5, z0=0.30):
    """Bille lâchée qui rebondit : le vol libre est long, le choc est bref.

    C'est le cas où le PAS PILOTÉ PAR LE CONTACT paie — la décomposition en
    temps sous sa forme la plus simple. L'estimateur n'a rien d'empirique :
    on connaît la distance à la collision, donc on sait exactement quand
    resserrer.
    """
    import time
    N = Noyau([0.0, 0.0, -G])
    b = N.corps("b", m, list((2.0 / 5.0 * m * R * R * np.eye(3)).ravel()), [0, 0, z0])
    N.contact("sol", b, [0, 0, 0], R, [0, 0, 1.0], [0, 0, 0], k, expo=1.5, c=c)
    pc = (h_libre, h_contact, marge) if h_contact else None
    t0 = time.time()
    N.simule(t_end, h_libre, tous=10 ** 9, pas_contact=pc)
    return N.etat()[1][0][2], time.time() - t0


def appariement(rapide=False):
    """L'APPARIEMENT automatique : exactitude, puis complexité.

    Déclarer chaque paire à la main ne passe pas l'échelle — c'est le fossé
    qui restait face à un code de contact. La grille de hachage spatial ne
    teste que les voisines des 27 cellules adjacentes : O(N) tant que la
    densité est bornée, contre O(N²) pour le balayage naïf. Même fin qu'une
    hiérarchie de volumes englobants, avec une structure plate reconstruite à
    chaque pas au lieu d'être mise à jour.
    """
    R, m, k = 0.02, 0.5, 1e6

    def deux(auto):
        N = Noyau([0.0, 0.0, 0.0])
        jj = lambda mm: list((2.0 / 5.0 * mm * R * R * np.eye(3)).ravel())
        a = N.corps("a", m, jj(m), [-0.05, 0, 0], v=[1.0, 0, 0])
        b = N.corps("b", m, jj(m), [0.05, 0, 0])
        if auto:
            N.sphere(a, [0, 0, 0], R)
            N.sphere(b, [0, 0, 0], R)
            N.appariement(k=k, marge=5e-3)
        else:
            N.contact("cc", a, [0, 0, 0], R, k=k, b=b, pb=[0, 0, 0], rayon_b=R)
        t = 0.0
        while t < 0.2:
            N.simule(t + 2e-3, 2e-6, tous=10 ** 9)
            t += 2e-3
            if N.etat()[3][1][0] > 1e-9 and t > 0.02:
                break
        return N.etat()[3][0][0], N.etat()[3][1][0]

    d, au = deux(False), deux(True)
    print(f"  appariement AUTO contre paire déclarée : v_a {au[0]:+.6f} contre {d[0]:+.6f}"
          f"   ·   v_b {au[1]:+.6f} contre {d[1]:+.6f}")
    assert abs(au[0] - d[0]) < 1e-9 and abs(au[1] - d[1]) < 1e-9, (au, d)

    rng = np.random.default_rng(3)
    tailles = (50, 200, 800) if rapide else (50, 200, 800, 3200)
    res = []
    for n in tailles:
        N = Noyau([0.0, 0.0, 0.0])
        L = (n / 300.0) ** (1 / 3) * 0.5           # densité constante
        for i in range(n):
            N.corps(f"b{i}", 0.01, list((1e-6 * np.eye(3)).ravel()),
                    list(rng.uniform(-L, L, 3)))
            N.sphere(i, [0, 0, 0], 0.01)
        N.appariement(k=1e5, marge=2e-3)
        N.simule(20 * 1e-4, 1e-4, tous=10 ** 9)
        pa, ta = N.appariement_stats()
        res.append((n, ta / 20 * 1e3, pa))
    print("  complexité (densité constante, 20 pas) :")
    for n, t, pa in res:
        print(f"      N = {n:5d}   {pa:3d} paires   {t:8.4f} ms/pas   "
              f"naïf O(N²) : {(n / res[0][0]) ** 2 * res[0][1]:9.4f} ms")
    x = np.log([a for a, _, _ in res])
    y = np.log([b for _, b, _ in res])
    expo = float(np.polyfit(x, y, 1)[0])
    print(f"      exposant mesuré {expo:.3f}   (le balayage naïf vaut 2,0)")
    assert expo < 1.3, expo
    return dict(exposant=expo, points=res)


def pente(th_deg, mu, t_end=0.4, h=2e-5, m=0.5, k=1e6, R=0.02, nonlisse=False):
    """Bille sur un plan incliné PORTÉ PAR UN CORPS : accélération le long de
    la pente.

    Deux régimes, deux formules exactes, et une transition connue :

        μ ≥ (2/7)·tan θ   roulement sans glissement   a = (5/7)·g·sin θ
        μ <  (2/7)·tan θ   glissement                 a = g(sin θ − μ cos θ)

    C'est le contrôle le plus complet du contact : il vérifie d'un coup la
    force normale, le frottement de Coulomb régularisé, le COUPLE de
    frottement (sans lui la bille glisserait au lieu de rouler, donc les bras
    de levier sont jugés aussi), et le plan attaché à un corps.
    """
    from math import cos, radians, sin
    th = radians(th_deg)
    N = Noyau([0.0, 0.0, -G])
    sol = N.corps("sol", 1.0, list((1e-3 * np.eye(3)).ravel()), [0, 0, 0])
    N.liaison("fixe", None, sol)
    n = [-sin(th), 0.0, cos(th)]
    b = N.corps("bille", m, list((2.0 / 5.0 * m * R * R * np.eye(3)).ravel()),
                [R * n[0], 0.0, R * n[2]])
    N.contact("c", b, [0, 0, 0], R, n, [0, 0, 0], k, expo=1.5, c=10.0, mu=mu,
              v_eps=1e-4, b=sol, rayon_b=0.0, nonlisse=nonlisse)
    N.simule(t_end, h, tous=10 ** 9)
    # ⚠ la plus grande PENTE est −(cos θ, 0, sin θ) : la projection de −ẑ sur
    # le plan. Se tromper de signe en x rendait six cas identiques et faux.
    t_hat = -np.array([cos(th), 0.0, sin(th)])
    return float(np.array(N.etat()[3][1]) @ t_hat) / t_end


def nonlisse():
    """CONTACT NON LISSE — multiplicateur sous complémentarité, loi d'impact en vitesse.

    Ce que la pénalité ne sait pas faire : une pénétration NULLE (à l'arrondi,
    pas à k^{−1}), une restitution CHOISIE plutôt que résultant d'un
    amortissement, et une force normale qui est un multiplicateur. Contrôles :
    la bille lâchée de H touche à √(2H/g), ne pénètre pas, ne rebondit pas à
    e = 0 et remonte à e²·H puis e⁴·H à e = 0,5 ; posée, λ = mg ; sur la
    pente, glissement et roulement suivent leurs formules exactes avec le
    Coulomb sur λ.
    """
    from math import cos, radians, sin, sqrt
    m, R, H, h = 0.5, 0.02, 0.1, 1e-4
    out = {}
    for e in (0.0, 0.5):
        N, b, _ = bille(c=10.0, mu=0.3, z0=R + H, nonlisse=True, e=e)
        traj = N.simule(0.6, h, tous=1)
        t = np.array([f[0] for f in traj])
        gap = np.array([f[1][0][2] for f in traj]) - R
        vz = np.array([f[3][0][2] for f in traj])
        # premier instant où la vitesse cesse de descendre : l'impact
        i_imp = int(np.argmax(np.diff(vz) > 0.5 * G * h))
        apex = []                                 # sommets après chaque impact
        for i in range(i_imp, len(vz) - 1):
            if vz[i] > 0 >= vz[i + 1]:
                apex.append(gap[i])
        out[e] = dict(t_imp=t[i_imp], pen=-gap.min(), apex=apex[:2], lam=N.reactions()[0][1][0])
        print(f"║ e = {e}: impact à {t[i_imp]:.4f} s (√(2H/g) = {sqrt(2 * H / G):.4f}), pénétration max "
              f"{-gap.min():.1e} m, sommets {', '.join(f'{a:.4e}' for a in apex[:2]) or 'aucun'}"
              f" (e²H = {e * e * H:.4e}, e⁴H = {e ** 4 * H:.4e}), λ final {out[e]['lam']:.4f} N (mg {m * G:.4f})")
    assert abs(out[0.0]["t_imp"] - sqrt(2 * H / G)) < 2 * h, out
    assert all(v["pen"] < 1e-5 for v in out.values()), out
    assert out[0.0]["apex"] == [] or max(out[0.0]["apex"]) < 1e-6, out
    assert abs(out[0.5]["apex"][0] / (0.25 * H) - 1) < 0.03, out
    assert abs(out[0.5]["apex"][1] / (0.0625 * H) - 1) < 0.06, out
    assert all(abs(v["lam"] / (m * G) - 1) < 0.05 for v in out.values()), out
    for th, mu in ((20.0, 0.05), (20.0, 0.5)):
        gl = mu < 2.0 / 7.0 * np.tan(radians(th))
        ref = G * (sin(radians(th)) - mu * cos(radians(th))) if gl else 5.0 / 7.0 * G * sin(radians(th))
        a = pente(th, mu, nonlisse=True)
        print(f"║ pente {th:.0f}° μ = {mu} ({'glisse' if gl else 'roule'}) : a = {a:.5f}  théorie {ref:.5f}"
              f"  écart {100 * (a / ref - 1):+.3f} %")
        assert abs(a / ref - 1) < 1e-3, (th, mu, a, ref)
    return out


def multiples(h=1e-4):
    """CONTACTS MULTIPLES QUI BASCULENT — ce qu'une bille seule ne teste pas.

    Deux scènes sur des lignes K DIFFÉRENTES qui s'activent et se relâchent
    dans le même vol : (1) trois billes empilées, lâchées avec un jeu, e = 0 —
    au repos les multiplicateurs valent 3mg, 2mg, mg (chaque contact porte ce
    qui est au-dessus) ; (2) trois billes alignées sans pesanteur, e = 1 — la
    première frappe la deuxième qui frappe la troisième : à masses égales la
    vitesse se TRANSMET entièrement (v, 0, 0 → 0, 0, v), quantité de mouvement
    et énergie conservées à e = 1. Chaque impact active une ligne, en relâche
    une autre, et redémarre le schéma.
    """
    m, R = 0.5, 0.02
    jj = list((2.0 / 5.0 * m * R * R * np.eye(3)).ravel())
    # (1) la pile
    N = Noyau([0.0, 0.0, -G])
    ids = [N.corps(f"b{i}", m, jj, [0.0, 0.0, R + i * (2 * R + 0.002)]) for i in range(3)]
    N.contact("sol", ids[0], [0, 0, 0], R, [0, 0, 1.0], [0, 0, 0], 1e6, nonlisse=True, restitution=0.0)
    for i in (0, 1):
        N.contact(f"b{i}b{i+1}", ids[i], [0, 0, 0], R, b=ids[i + 1], pb=[0, 0, 0], rayon_b=R,
                  k=1e6, nonlisse=True, restitution=0.0)
    N.simule(0.5, h, tous=10 ** 9)
    lam = [N.reactions()[k][1][0] for k in range(3)]
    z = [N.etat()[1][i][2] for i in range(3)]
    print(f"║ pile e = 0 : λ sol {lam[0]:.4f} (3mg {3*m*G:.4f}) · λ 0–1 {lam[1]:.4f} (2mg {2*m*G:.4f}) · "
          f"λ 1–2 {lam[2]:.4f} (mg {m*G:.4f}) · écarts {z[1]-z[0]-2*R:.1e}, {z[2]-z[1]-2*R:.1e} m")
    assert all(abs(l / (n * m * G) - 1) < 0.05 for l, n in zip(lam, (3, 2, 1))), lam
    assert all(abs(z[i + 1] - z[i] - 2 * R) < 1e-5 for i in (0, 1)), z
    # (2) le pendule de Newton à plat
    N = Noyau([0.0, 0.0, 0.0])
    v0 = 0.5
    ids = [N.corps(f"c{i}", m, jj, [i * (2 * R + 0.01), 0.0, 0.0], v=[v0 if i == 0 else 0.0, 0, 0]) for i in range(3)]
    for i in (0, 1):
        N.contact(f"c{i}c{i+1}", ids[i], [0, 0, 0], R, b=ids[i + 1], pb=[0, 0, 0], rayon_b=R,
                  k=1e6, nonlisse=True, restitution=1.0)
    N.simule(0.2, h, tous=10 ** 9)
    v = [N.etat()[3][i][0] for i in range(3)]
    p, e = m * sum(v), 0.5 * m * sum(x * x for x in v)
    print(f"║ trois billes e = 1 : v = {v[0]:+.4f}, {v[1]:+.4f}, {v[2]:+.4f} (théorie 0, 0, {v0}) · "
          f"p {p/(m*v0):.4f}·p₀ · E {e/(0.5*m*v0*v0):.4f}·E₀")
    assert abs(v[2] / v0 - 1) < 0.02 and abs(v[0]) < 0.02 * v0 and abs(v[1]) < 0.02 * v0, v
    assert abs(p / (m * v0) - 1) < 1e-3 and abs(e / (0.5 * m * v0 ** 2) - 1) < 0.03, (p, e)
    return dict(lam=lam, v=v)


def barriere_ipc(k=1e3, d_hat=1e-3, c=5.0, m=0.5, R=0.02, t_end=1.5, h=2e-5):
    """Bille posée, loi de BARRIÈRE IPC : elle ne touche jamais le sol.

    Li & al. 2020 (Incremental Potential Contact) remplacent la force qui
    croît AVEC l'enfoncement par un potentiel qui DIVERGE quand la distance
    tend vers zéro :

        B(d) = −(d − d̂)² ln(d/d̂)   sur 0 < d < d̂,   0 au-delà

    Trois propriétés que la pénalité de Hertz n'a pas : l'interpénétration
    devient IMPOSSIBLE au lieu d'être petite ; le support est COMPACT, donc le
    conditionnement n'est pas dégradé loin du contact ; et B est C², donc la
    force est C¹ et l'AD passe au travers.

    Ce qui N'EST PAS repris d'IPC, et il faut le dire : la formulation
    complète fait de chaque pas une MINIMISATION d'un potentiel incrémental,
    avec une recherche linéaire filtrée par détection de collision continue —
    c'est ELLE qui donne la garantie topologique. Ici la barrière est une loi
    de force dans notre schéma existant : elle rend l'interpénétration
    extrêmement improbable, elle ne la rend pas impossible par construction.
    """
    N = Noyau([0.0, 0.0, -G])
    b = N.corps("b", m, list((2.0 / 5.0 * m * R * R * np.eye(3)).ravel()), [0, 0, R + 5e-4])
    N.contact("sol", b, [0, 0, 0], R, [0, 0, 1.0], [0, 0, 0], k, expo=-1.0,
              d_hat=d_hat, c=c)
    N.simule(t_end, h, tous=10 ** 9)
    _, d, f = N.contacts()[0]
    return d, f


def tas(n, R=0.01, k=1e5, mu=0.3, t_end=0.6, h=2e-4, m=0.01, graine=11):
    """N billes lâchées en pluie sur un sol, contacts découverts et IPC.

    Le cas où tout joue ensemble : appariement automatique, contacts multiples
    simultanés, frottement, barrière. C'est aussi le terrain que les codes de
    contact mettent en avant (sols granulaires), et il pose une question à
    laquelle un enfoncement moyen ne répond pas : **une seule paire
    a-t-elle pénétré ?**
    """
    import time
    rng = np.random.default_rng(graine)
    N = Noyau([0.0, 0.0, -G])
    ids = []
    par = int(np.ceil(n ** (1 / 3)))
    c = 0
    for i in range(par):
        for j in range(par):
            for kk in range(par):
                if c >= n:
                    break
                p = [(i - par / 2) * 2.2 * R, (j - par / 2) * 2.2 * R, 0.03 + kk * 2.2 * R]
                b = N.corps(f"b{c}", m, list((2.0 / 5.0 * m * R * R * np.eye(3)).ravel()),
                            p, v=list(rng.normal(0, 0.01, 3)))
                ids.append(b)
                N.sphere(b, [0, 0, 0], R)
                c += 1
    kw = dict(k=k, expo=-1.0, d_hat=2e-4, c=3.0, mu=mu, v_eps=1e-3)
    for b in ids:
        N.contact(f"s{b}", b, [0, 0, 0], R, [0, 0, 1.0], [0, 0, 0], **kw)
    N.appariement(marge=3e-3, **kw)
    t0 = time.time()
    N.simule(t_end, h, tous=10 ** 9)
    dt = time.time() - t0
    d_max = max((d for _, d, _ in N.contacts()), default=-1.0)
    return dict(s=dt, npas=int(t_end / h), paires=N.appariement_stats()[0], d_max=d_max)


def tir(v0=50.0, h=1e-3, ccd=False, k=1e8, R=0.005, m=0.01, d_hat=1e-3):
    """Bille lancée à grande vitesse contre une autre : traverse-t-elle ?

    À pas grossier, une sphère rapide passe de l'autre côté SANS qu'aucun
    instant échantillonné ne la voie — les deux extrémités du pas sont saines,
    la trajectoire ne l'est pas. C'est le défaut que la détection de collision
    CONTINUE corrige, et c'est ce qui donne à IPC sa garantie topologique.

    Pour nos primitives le test est analytique : la distance entre centres est
    |Δc₀ + t·Δv|, dont le minimum sur [0,1] se calcule en fermé.
    """
    N = Noyau([0.0, 0.0, 0.0])
    jj = lambda mm: list((2.0 / 5.0 * mm * R * R * np.eye(3)).ravel())
    a = N.corps("a", m, jj(m), [-0.05, 0, 0], v=[v0, 0, 0])
    b = N.corps("b", m, jj(m), [0.05, 0, 0])
    N.contact("cc", a, [0, 0, 0], R, k=k, expo=-1.0, d_hat=d_hat, c=1.0,
              b=b, pb=[0, 0, 0], rayon_b=R)
    N.simule(0.3 / v0 * 2, h, tous=10 ** 9, ccd=ccd)
    e = N.etat()
    return e[1][0][0], e[1][1][0], e[3][0][0], e[3][1][0]


def tir_paroi(forme, ccd, v0=50.0, h=1e-3, k=1e8, R=0.005, m=0.01, d_hat=1e-3):
    """Bille rapide contre une PAROI mince fixée au bâti — boîte, cylindre ou
    maillage. Retourne x de la bille à la fin.

    Le pas grossier saute par-dessus la paroi : à 50 m/s et h = 1e-3 la bille
    va de −0,03 à +0,02 en un pas, et aux deux bouts elle est hors de portée de
    la barrière. Le CCD juge la CAPSULE balayée [x₀, x₁] de rayon R contre la
    primitive — pour la boîte et le cylindre par le minimum exact de la
    distance le long du segment (convexe), pour le maillage par la distance
    segment/triangle en forme fermée sous le BVH — et rejoue le pas plus
    court jusqu'à ce que la barrière voie la paroi.
    """
    N = Noyau([0.0, 0.0, 0.0])
    a = N.corps("a", m, list((2.0 / 5.0 * m * R * R * np.eye(3)).ravel()), [-0.03, 0, 0],
                v=[v0, 0, 0])
    b = N.corps("b", 1.0, [1e-2, 0, 0, 0, 1e-2, 0, 0, 0, 1e-2], [0.0] * 3)
    N.liaison("f", None, b, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    kw = dict(k=k, expo=-1.0, d_hat=d_hat, c=1.0, b=b, pb=[0, 0, 0])
    if forme == "boîte":
        N.contact("cc", a, [0, 0, 0], R, demi=[0.01, 0.1, 0.1], **kw)
    elif forme == "cylindre":
        N.contact("cc", a, [0, 0, 0], R, cylindre=([0, 1, 0], 0.01, 0.1), **kw)
    else:
        mi = N.maillage("paroi", b, [0, -.1, -.1, 0, .1, -.1, 0, .1, .1, 0, -.1, .1],
                        [0, 1, 2, 0, 2, 3])
        N.contact("cc", a, [0, 0, 0], R, maille=mi, **kw)
    N.simule(0.06 / v0, h, tous=10 ** 9, ccd=ccd)
    return N.etat()[1][0][0]


def capsules():
    """LA CAPSULE — un segment dilaté, et la sphère en est le cas dégénéré.

    Le noyau n'avait que des sphères et des plans ; une tige, une barre, un
    doigt, un axe, une biellette ne s'en approchent pas. La capsule est la
    primitive qui couvre le plus de mécanismes réels **sans maillage**, et
    elle ne coûte qu'une distance segment–segment en forme fermée.

    Trois régimes à vérifier, parce que l'algorithme a trois branches et que
    seule la mesure dit qu'on est dans la bonne :
      · les points les plus proches à l'INTÉRIEUR des deux segments ;
      · un point le plus proche à une EXTRÉMITÉ (segments décalés) ;
      · le cas DÉGÉNÉRÉ p₁ = p₀, qui doit rendre exactement la sphère
        d'avant — c'est la non-régression de tout le contact existant.

    ⚠ Déclaré : les `clamp` rendent la distance C⁰ mais pas C¹ aux TRANSITIONS
    de branche. L'AD y prend la dérivée de la branche active, ce qui est
    correct partout sauf sur ces surfaces de mesure nulle.
    """
    ji = [1e-4, 0, 0, 0, 1e-4, 0, 0, 0, 1e-4]

    def enfoncement(pa0, pa1, pb0, pb1, ra, rb, dz):
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("a", 1.0, ji, [0.0, 0.0, 0.0])
        n.corps("b", 1.0, ji, [0.0, 0.0, dz])
        n.contact("c", 0, list(pa0), ra, b=1, pb=list(pb0), rayon_b=rb,
                  p1=list(pa1), p1b=list(pb1), k=1.0)
        n.simule(1e-6, 1e-6, tous=10 ** 9)
        return n.contacts()[0][1]

    print("╔═ vinkulum — CAPSULES (segment dilaté), distance en forme fermée")
    pires = []
    # 1. segments CROISÉS : le contact est à l'intérieur des deux
    for dz in (0.15, 0.20, 0.30):
        d = enfoncement([-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], 0.1, 0.1, dz)
        pires.append(abs(d - (0.2 - dz)))
        print(f"║ croisées   dz {dz:.2f} : enfoncement {d:+.6f}  théorie {0.2 - dz:+.6f}")
    # 2. segments PARALLÈLES décalés : le contact passe par une EXTRÉMITÉ
    for dx in (0.0, 1.5, 2.5):
        d = enfoncement([-1, 0, 0], [1, 0, 0], [dx - 1, 0, 0], [dx + 1, 0, 0], 0.1, 0.1, 0.15)
        th = 0.2 - np.hypot(max(0.0, dx - 2.0), 0.15)
        pires.append(abs(d - th))
        print(f"║ parallèles dx {dx:.1f} : enfoncement {d:+.6f}  théorie {th:+.6f}")
    # 3. DÉGÉNÉRÉ p₁ = p₀ : la sphère d'avant, au dernier chiffre
    for dz in (0.15, 0.25):
        d = enfoncement([0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], 0.1, 0.1, dz)
        pires.append(abs(d - (0.2 - dz)))
    print(f"║ sphère (p₁ = p₀) : non-régression du contact existant")
    print(f"╚═ capsules OK — pire écart à la géométrie exacte {max(pires):.1e}")
    assert max(pires) < 1e-9, ("la distance capsule s'écarte de la géométrie", pires)
    return dict(pire=float(max(pires)))


def capsules_primitives():
    """LA CAPSULE contre BOÎTE, CYLINDRE et MAILLAGE — ce que le contact ne
    couvrait pas (il lisait p0 : une tige contre un carter n'était qu'une
    bille à son bout).

    Le point de l'axe le plus proche de la primitive est FIGÉ sur le pas,
    comme le point du maillage (`q_maille`) : abscisse s ∈ [0, 1] par le
    minimum exact de la distance convexe le long de l'axe (boîte, cylindre)
    ou par la distance segment/triangle en forme fermée sous le BVH
    (maillage, Ericson 5.1.10). Le résidu voit alors la sphère de l'axe qui
    s'y trouve — exact aux bouts, d'ordre 2 à l'intérieur, puisque la
    distance est stationnaire en s au minimum. Chaque cas est jugé contre la
    géométrie exacte, minimisée indépendamment (scipy, Brent borné).
    """
    from scipy.optimize import minimize_scalar
    ji = [1e-2, 0, 0, 0, 1e-2, 0, 0, 0, 1e-2]
    hb = np.array([0.3, 0.2, 0.1])
    R, L, r = 0.2, 0.3, 0.05
    PL = ([-1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1, 0], [0, 1, 2, 0, 2, 3])

    def enf(p0, p1, forme):
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("a", 1.0, ji, [0.0] * 3)
        n.corps("b", 1.0, ji, [0.0] * 3)
        kw = dict(b=1, pb=[0.0] * 3, p1=list(p1), k=1.0)
        if forme == "boîte":
            n.contact("c", 0, list(p0), r, demi=list(hb), **kw)
        elif forme == "cylindre":
            n.contact("c", 0, list(p0), r, cylindre=([0, 0, 1], R, L), **kw)
        else:
            n.maillage("m", 1, *PL)
            n.contact("c", 0, list(p0), r, maille=0, **kw)
        n.simule(1e-6, 1e-6, tous=10 ** 9)   # un pas : s est figé en fin de pas, puis lu
        return n.contacts()[0][1]

    def d_point(p, forme):
        p = np.asarray(p, float)
        if forme == "boîte":
            return np.linalg.norm(np.maximum(np.abs(p) - hb, 0.0))
        if forme == "cylindre":
            rho, z = np.hypot(p[0], p[1]), p[2]
            return np.hypot(max(rho - R, 0.0), max(abs(z) - L, 0.0))
        return np.linalg.norm(np.maximum(np.abs(p[:2]) - 1.0, 0.0).tolist() + [p[2]])

    def theorie(p0, p1, forme):
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
        f = lambda s: d_point(p0 + s * (p1 - p0), forme)
        res = minimize_scalar(f, bounds=(0.0, 1.0), method="bounded", options=dict(xatol=1e-13))
        return r - min(res.fun, f(0.0), f(1.0))

    cas = {
        "boîte": [("parallèle à la face", [-0.2, 0, 0.13], [0.2, 0, 0.13]),
                  ("bout le plus proche  ", [0, 0, 0.12], [0, 0, 0.5]),
                  ("au-dessus d'un coin  ", [0.4, 0.3, 0.2], [0.2, 0.1, 0.2]),
                  ("le long d'une arête  ", [0.5, 0.35, -0.3], [0.5, 0.35, 0.3])],
        "cylindre": [("parallèle à l'axe   ", [0.23, 0, -0.1], [0.23, 0, 0.1]),
                     ("perpendiculaire      ", [0.3, -0.5, 0], [0.3, 0.5, 0]),
                     ("bout sur le fond     ", [0, 0, 0.33], [0.1, 0, 0.6]),
                     ("par-dessus l'arête   ", [0.25, 0, 0.35], [0.15, 0, 0.35])],
        "maillage": [("parallèle au plan   ", [-0.5, 0, 0.03], [0.5, 0, 0.03]),
                     ("bout le plus proche  ", [0, 0, 0.02], [0.3, 0, 0.4]),
                     ("traverse le plan     ", [0, 0, -0.1], [0, 0, 0.1]),
                     ("au-delà du bord      ", [1.2, 0, 0.1], [1.6, 0, 0.1])],
    }
    print("╔═ vinkulum — CAPSULE contre boîte, cylindre et maillage (s figé sur le pas)")
    pires = []
    for forme, liste in cas.items():
        for nom, p0, p1 in liste:
            d, th = enf(p0, p1, forme), theorie(p0, p1, forme)
            pires.append(abs(d - th))
            print(f"║ {forme:8s} {nom} : enfoncement {d:+.6f}  théorie {th:+.6f}")
    # ACCORD : une tige qui tombe à plat sur une plaque, boîte contre maillage
    def chute(forme):
        n = Noyau([0.0, 0.0, -9.81])
        n.corps("t", 0.5, ji, [0.0, 0.0, 0.3])
        n.corps("sol", 1.0, ji, [0.0] * 3)
        n.liaison("f", None, 1, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
        kw = dict(b=1, p1=[0.2, 0, 0], k=1e6, expo=1.5, c=0.5)
        if forme == "boîte":
            n.contact("c", 0, [-0.2, 0, 0], 0.05, demi=[1, 1, 0.05], pb=[0, 0, -0.05], **kw)
        else:
            n.maillage("sol", 1, *PL)
            n.contact("c", 0, [-0.2, 0, 0], 0.05, maille=0, **kw)
        tr = n.simule(1.0, 2e-5, tous=100)
        z = [f[1][0][2] for f in tr]
        return min(z), max(z[len(z) // 2:])
    bo, ma = chute("boîte"), chute("maillage")
    print(f"║ tige à plat qui tombe — boîte z_min {bo[0]:+.5f} rebond {bo[1]:+.5f} · maillage {ma[0]:+.5f} / {ma[1]:+.5f}")
    assert abs(bo[0] - ma[0]) < 1e-6 and abs(bo[1] - ma[1]) < 1e-6, ("boîte et maillage ne s'accordent pas", bo, ma)
    print(f"╚═ capsules contre primitives OK — pire écart à la géométrie exacte {max(pires):.1e}")
    assert max(pires) < 1e-8, ("la capsule s'écarte de la géométrie", pires)
    return dict(pire=float(max(pires)))


def boites():
    """LA BOÎTE (OBB) — troisième primitive, et la dernière du plan v0.5.

    Avec la sphère et la capsule, elle couvre l'essentiel d'un mécanisme sans
    maillage : un bâti, une glissière, un carter, une came plane. Le point le
    plus proche est un CLAMP dans le repère de la boîte, donc exact — et les
    trois régimes (face, arête, coin) sont les trois façons dont le clamp
    mord, ce qui les rend vérifiables séparément.

    ⚠ NON COUVERT, déclaré : le centre de la primitive À L'INTÉRIEUR de la
    boîte. La distance y change de nature, et une pénalité ne doit jamais y
    arriver — si elle y arrive c'est la raideur qui est trop faible, et c'est
    ça qu'il faut corriger.
    """
    ji = [1e-2, 0, 0, 0, 1e-2, 0, 0, 0, 1e-2]
    h = np.array([0.3, 0.2, 0.1])
    r = 0.05

    def enf(p, ang=0.0):
        c, s2 = np.cos(ang), np.sin(ang)
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("s", 1.0, ji, list(p))
        n.corps("b", 1.0, ji, [0.0] * 3)
        n.contact("c", 0, [0.0] * 3, r, b=1, pb=[0.0] * 3, demi=list(h), k=1.0)
        if ang:
            t, rr, rot, v, w, vi = n.etat()
            n.pose_etat(t, rr, [rot[0], [c, -s2, 0, s2, c, 0, 0, 0, 1]], v, w, vi)
        n.simule(1e-6, 1e-6, tous=10 ** 9)
        return n.contacts()[0][1]

    def theorie(p, ang=0.0, hh=None):
        c, s2 = np.cos(ang), np.sin(ang)
        rm = np.array([[c, -s2, 0.0], [s2, c, 0.0], [0.0, 0.0, 1.0]])
        loc = rm.T @ np.asarray(p, float)
        return r - np.linalg.norm(np.maximum(np.abs(loc) - (h if hh is None else hh), 0.0))

    print("╔═ vinkulum — BOÎTE (OBB), point le plus proche par clamp")
    pires = []
    for nom, p in (("face ", [0.34, 0.0, 0.0]), ("face ", [0.0, 0.0, 0.13]),
                   ("arête", [0.33, 0.23, 0.0]), ("arête", [0.0, 0.22, 0.12]),
                   ("coin ", [0.32, 0.22, 0.12]), ("loin ", [0.5, 0.4, 0.3])):
        m, t = enf(np.array(p)), theorie(p)
        pires.append(abs(m - t))
        print(f"║ {nom} ({p[0]:+.2f},{p[1]:+.2f},{p[2]:+.2f}) : enfoncement {m:+.6f}  théorie {t:+.6f}")
    # · LA BOÎTE TOURNE AVEC SON CORPS : sans ça elle ne serait qu'une AABB
    #   déguisée, et rien dans les cas ci-dessus ne l'aurait montré.
    hp = np.array([0.3, 0.1, 0.1])
    for ang in (30.0, 90.0, 225.0):
        c, s2 = np.cos(np.radians(ang)), np.sin(np.radians(ang))
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("s", 1.0, ji, [0.0, 0.34, 0.0])
        n.corps("b", 1.0, ji, [0.0] * 3)
        n.contact("c", 0, [0.0] * 3, r, b=1, pb=[0.0] * 3, demi=list(hp), k=1.0)
        t0, rr, rot, v, w, vi = n.etat()
        n.pose_etat(t0, rr, [rot[0], [c, -s2, 0, s2, c, 0, 0, 0, 1]], v, w, vi)
        n.simule(1e-6, 1e-6, tous=10 ** 9)
        m = n.contacts()[0][1]
        th = theorie([0.0, 0.34, 0.0], np.radians(ang), hp)
        pires.append(abs(m - th))
        print(f"║ tournée {ang:5.0f}° : enfoncement {m:+.6f}  théorie {th:+.6f}")
    print(f"╚═ boîtes OK — pire écart à la géométrie exacte {max(pires):.1e}")
    assert max(pires) < 1e-12, ("la distance à la boîte s'écarte de la géométrie", pires)
    return dict(pire=float(max(pires)))


def cylindres():
    """LE CYLINDRE FINI — la primitive qui manquait, et ses QUATRE régimes.

    Un arbre, un galet, un rouleau, un tourillon : ni la capsule (bouts ronds)
    ni la boîte ne les décrivent. Le point le plus proche sépare l'axial
    (clampé à ±L) du radial (ramené à R), et les régimes se vérifient un par
    un contre la géométrie exacte : FLANC (|z| ≤ L, ρ ≥ R), FOND (|z| > L,
    ρ < R), ARÊTE circulaire (|z| > L, ρ ≥ R), et LOIN (aucun contact). Puis
    le cylindre TOURNE avec son corps — sans quoi ce ne serait qu'un cylindre
    d'axe fixe déguisé — et l'axe peut être quelconque dans le repère du corps.

    ⚠ NON COUVERT, déclaré : le centre de la primitive À L'INTÉRIEUR du
    cylindre (même limite que la boîte). Le CCD, lui, est là depuis le
    7 sept. (`tir_paroi`) : la sphère balayée sur le pas est une capsule, et
    la traversée se juge par la distance de son segment au cylindre —
    jusqu'au 5 sept. boîte et maillage tombaient DANS la branche « plan » du
    CCD avec la normale par défaut, puis ils en étaient exclus.
    """
    ji = [1e-2, 0, 0, 0, 1e-2, 0, 0, 0, 1e-2]
    R, L, r = 0.2, 0.3, 0.05

    def enf(p, axe=(0.0, 0.0, 1.0), rot=None):
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("s", 1.0, ji, list(p))
        n.corps("b", 1.0, ji, [0.0] * 3)
        n.contact("c", 0, [0.0] * 3, r, b=1, pb=[0.0] * 3, cylindre=(list(axe), R, L), k=1.0)
        if rot is not None:
            t, rr, ro, v, w, vi = n.etat()
            n.pose_etat(t, rr, [ro[0], list(np.asarray(rot).ravel())], v, w, vi)
        n.simule(1e-6, 1e-6, tous=10 ** 9)
        return n.contacts()[0][1]

    def theorie(p, axe=(0.0, 0.0, 1.0), rot=None):
        rm = np.eye(3) if rot is None else np.asarray(rot, float)
        loc = rm.T @ np.asarray(p, float)
        e = np.asarray(axe, float) / np.linalg.norm(axe)
        z = loc @ e
        rv = loc - z * e
        rho = np.linalg.norm(rv)
        dz = max(abs(z) - L, 0.0)
        dr = max(rho - R, 0.0)
        return r - np.hypot(dz, dr)

    print("╔═ vinkulum — CYLINDRE FINI, point le plus proche en quatre régimes")
    pires = []
    for nom, p in (("flanc", [0.24, 0.0, 0.1]), ("flanc", [0.0, -0.23, -0.25]),
                   ("fond ", [0.05, 0.05, 0.34]), ("fond ", [-0.1, 0.0, -0.33]),
                   ("arête", [0.23, 0.0, 0.33]), ("arête", [0.0, 0.22, -0.34]),
                   ("loin ", [0.5, 0.5, 0.6])):
        m, t = enf(p), theorie(p)
        pires.append(abs(m - t))
        print(f"║ {nom} ({p[0]:+.2f},{p[1]:+.2f},{p[2]:+.2f}) : enfoncement {m:+.6f}  théorie {t:+.6f}")
    # · LE CYLINDRE TOURNE AVEC SON CORPS, et son axe est quelconque
    def rz(ang):
        c, s2 = np.cos(np.radians(ang)), np.sin(np.radians(ang))
        return np.array([[c, -s2, 0.0], [s2, c, 0.0], [0.0, 0.0, 1.0]])
    for ang in (30.0, 90.0, 225.0):
        p = [0.24, 0.0, 0.32]                       # arête, axe x du corps
        m = enf(p, axe=(1.0, 0.0, 0.0), rot=rz(ang))
        th = theorie(p, axe=(1.0, 0.0, 0.0), rot=rz(ang))
        pires.append(abs(m - th))
        print(f"║ tourné {ang:5.0f}°, axe x : enfoncement {m:+.6f}  théorie {th:+.6f}")
    # · CONTRÔLE NÉGATIF : un point à ρ < R et |z| < L (intérieur) est HORS
    #   domaine — le noyau rend ~r (distance nulle), jamais un signe faux ; et
    #   une sphère exactement sur le flanc (ρ = R + r) a un enfoncement nul.
    m0 = enf([R + r, 0.0, 0.0])
    pires.append(abs(m0))
    print(f"║ sphère affleurant le flanc : enfoncement {m0:+.2e} (doit être 0)")
    print(f"╚═ cylindres OK — pire écart à la géométrie exacte {max(pires):.1e}")
    assert max(pires) < 1e-12, ("la distance au cylindre s'écarte de la géométrie", pires)
    # ⚠ la fonction refuse un cylindre sans corps porteur : c'est le domaine
    try:
        n = Noyau([0.0] * 3)
        n.corps("s", 1.0, ji, [0.0] * 3)
        n.contact("c", 0, [0.0] * 3, r, cylindre=([0, 0, 1], R, L), k=1.0)
        raise AssertionError("un cylindre sans porteur a été accepté")
    except ValueError:
        pass
    return dict(pire=float(max(pires)))


def _icosphere(n):
    """Sphère triangulée par subdivision de l'octaèdre — géométrie de test
    dont la distance exacte est connue analytiquement."""
    v = [np.array(x, float) for x in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                                      (0, -1, 0), (0, 0, 1), (0, 0, -1))]
    f = [(0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
         (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5)]
    for _ in range(n):
        nf, mid = [], {}

        def m(a, b):
            k = (min(a, b), max(a, b))
            if k not in mid:
                p = v[a] + v[b]
                p /= np.linalg.norm(p)
                v.append(p)
                mid[k] = len(v) - 1
            return mid[k]

        for a, b, c in f:
            ab, bc, ca = m(a, b), m(b, c), m(c, a)
            nf += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
        f = nf
    return np.array(v), np.array(f)


def maillages():
    """GÉOMÉTRIE QUELCONQUE — maillage triangulaire et BVH.

    Sphères, capsules et boîtes couvrent beaucoup ; elles ne couvrent pas une
    pièce usinée, un carter, une came. C'est ce maillage qui fait passer le
    contact de « trois primitives » à « n'importe quelle forme », et c'est le
    dernier écart fonctionnel de fond avec un solveur généraliste.

    Trois contrôles, aucun interne :
      1. une SPHÈRE maillée doit rendre sa distance ANALYTIQUE, à l'erreur de
         facettage près — laquelle doit décroître en O(h²) quand on raffine ;
      2. le BVH doit ÉLAGUER : l'exposant du coût en nombre de triangles doit
         être franchement sous 1 (le balayage naïf vaut 1) ;
      3. un contact DYNAMIQUE sur maillage doit donner le même mouvement que
         le même contact sur le plan analytique — sinon la géométrie est
         juste et la mécanique non.
    """
    import time
    ji = [1e-5, 0, 0, 0, 1e-5, 0, 0, 0, 1e-5]
    print("╔═ vinkulum — MAILLAGE triangulaire et BVH (géométrie quelconque)")

    # 1. distance analytique et convergence du facettage
    rng = np.random.default_rng(7)
    pts = []
    for _ in range(1500):
        q = rng.normal(size=3)
        pts.append(list(q * (1.5 + rng.random()) / np.linalg.norm(q)))
    errs = {}
    for sub in (1, 2, 3, 4):
        v, f = _icosphere(sub)
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("m", 1.0, ji, [0.0] * 3)
        mi = n.maillage("s", 0, list(v.ravel()), [int(x) for x in f.ravel()])
        errs[len(f)] = max(abs(n.distance_maillage(mi, p)[0] - (np.linalg.norm(p) - 1.0))
                           for p in pts)
    cles = sorted(errs)
    print("║ sphère maillée contre sa distance exacte : " +
          " · ".join(f"{k} tri {errs[k]:.4f}" for k in cles))
    # le facettage tombe en O(h²) : ×4 de triangles ⇒ ÷~4 d'erreur
    rap = [errs[cles[i]] / errs[cles[i + 1]] for i in range(len(cles) - 1)]
    assert min(rap) > 2.5, ("le facettage ne converge pas en O(h²)", errs, rap)

    # 2. le BVH élague : exposant du coût en N
    ts = []
    for sub in (3, 4, 5, 6):
        v, f = _icosphere(sub)
        n = Noyau([0.0, 0.0, 0.0])
        n.corps("m", 1.0, ji, [0.0] * 3)
        mi = n.maillage("s", 0, list(v.ravel()), [int(x) for x in f.ravel()])
        n.distance_maillage(mi, pts[0])
        t0 = time.perf_counter()
        for p in pts:
            n.distance_maillage(mi, p)
        ts.append((len(f), (time.perf_counter() - t0) / len(pts) * 1e6))
    ex = [np.log(ts[i + 1][1] / ts[i][1]) / np.log(ts[i + 1][0] / ts[i][0])
          for i in range(len(ts) - 1)]
    print(f"║ coût de requête {ts[0][1]:.2f} → {ts[-1][1]:.2f} µs pour {ts[0][0]} → {ts[-1][0]} "
          f"triangles : exposant {np.mean(ex):+.2f} (naïf +1,00)")
    assert np.mean(ex) < 0.7, ("le BVH n'élague pas", ts, ex)

    # 3. contact DYNAMIQUE : maillage contre plan analytique
    def chute(par_maillage):
        n = Noyau([0.0, 0.0, -9.81])
        n.corps("b", 0.5, ji, [0.0, 0.0, 0.3])
        n.corps("sol", 1.0, ji, [0.0] * 3)
        n.liaison("f", None, 1, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
        if par_maillage:
            n.maillage("sol", 1, [-1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1, 0], [0, 1, 2, 0, 2, 3])
            n.contact("c", 0, [0.0] * 3, 0.05, b=1, maille=0, k=1e6, expo=1.5, c=0.5)
        else:
            n.contact("c", 0, [0.0] * 3, 0.05, normale=[0, 0, 1], origine=[0.0] * 3,
                      k=1e6, expo=1.5, c=0.5)
        tr = n.simule(1.2, 2e-5, tous=100)
        z = [f[1][0][2] for f in tr]
        return min(z), max(z[len(z) // 2:])

    pl, ma = chute(False), chute(True)
    print(f"║ chute sur le sol — plan analytique z_min {pl[0]:+.5f} rebond {pl[1]:+.5f} · "
          f"maillage {ma[0]:+.5f} / {ma[1]:+.5f}")
    assert abs(ma[0] - pl[0]) < 1e-6 and abs(ma[1] - pl[1]) < 1e-6, ("le maillage ne reproduit pas le plan", pl, ma)
    print("╚═ maillages OK — la géométrie quelconque entre dans le contact")
    return dict(facettage=errs, exposant=float(np.mean(ex)))


def _penalite_ccd(rapide=False):
    print("╔═ vinkulum — contact par pénalité : Hertz statique, puis restitution")

    # 1. ENFONCEMENT STATIQUE : k δ^1,5 = mg, donc δ = (mg/k)^(2/3). EXACT.
    m, k, R = 0.5, 1e6, 0.02
    d_th = (m * G / k) ** (2.0 / 3.0)
    N, b, p = bille(k=k, c=10.0, m=m, R=R)
    N.simule(1.0, 2e-5, tous=10 ** 9)
    _, d, f = N.contacts()[0]
    print(f"║ enfoncement statique {1e6 * d:9.4f} µm   théorie {1e6 * d_th:9.4f} µm   "
          f"écart {100 * abs(d / d_th - 1):.4f} %")
    print(f"║ force de contact     {f:9.4f} N    poids   {m * G:9.4f} N")
    assert abs(d / d_th - 1) < 1e-3, (d, d_th)
    assert abs(f / (m * G) - 1) < 1e-4, (f, m * G)

    # 2. RESTITUTION : à c = 0 le contact est CONSERVATIF, donc e = 1
    #    exactement — c'est le contrôle qui ne demande aucune corrélation.
    e0 = restitution(0.0)
    print(f"║ restitution à c = 0 : e = {e0:.6f}   (conservatif : 1 exactement)")
    assert abs(e0 - 1.0) < 2e-3, e0

    # 3. et e doit DÉCROÎTRE avec c, de façon monotone — c'est ce que
    #    l'amortissement de Hunt–Crossley dissipe
    cs = (0.0, 0.5, 1.0, 2.0) if not rapide else (0.0, 1.0)
    es = [restitution(cc) for cc in cs]
    print("║ restitution contre amortissement : " + "  ".join(
        f"c={cc:.1f} → e={ee:.4f}" for cc, ee in zip(cs, es)))
    assert all(es[i + 1] < es[i] for i in range(len(es) - 1)), es
    # la pente initiale : e ≈ 1 − α·c·v₀ — on publie α mesuré plutôt que de
    # recopier une constante de la littérature qu'on ne peut pas vérifier ici
    # LA PENTE SE PREND À PETIT c : e n'est pas linéaire en c sur toute la
    # plage, donc le dernier point sous-estime la pente initiale de 30 %.
    a0 = (1.0 - es[1]) / (cs[1] * 0.5)
    af = (1.0 - es[-1]) / (cs[-1] * 0.5)
    print(f"║ pente à petit c : e ≈ 1 − {a0:.3f}·c·v₀   (8/15 = 0,533 dans la"
          f" littérature pour cette forme)   ·   à c = {cs[-1]:.0f} elle vaut {af:.3f}")
    assert 0.4 < a0 < 0.75, a0

    # 4. CHOC À DEUX CORPS MOBILES — le contact n'est plus contre un plan fixe
    print("║ choc central de deux billes libres, c = 0 (élastique) :")
    for m1, m2, lbl in ((0.5, 0.5, "masses égales"), (0.5, 1.5, "m₂ = 3 m₁")):
        va, vb, dp, de = choc(m1, m2, 1.0)
        u1, u2 = (m1 - m2) / (m1 + m2), 2 * m1 / (m1 + m2)
        print(f"║   {lbl:14} v_a {va:+.6f} v_b {vb:+.6f}   théorie {u1:+.6f} {u2:+.6f}"
              f"   ΔP/P {dp:+.1e}  ΔE/E {de:+.1e}")
        assert abs(va - u1) < 1e-4 and abs(vb - u2) < 1e-4, (va, vb, u1, u2)
        assert abs(dp) < 1e-12, dp
        assert abs(de) < 1e-6, de
    # 5. PAS PILOTÉ PAR LE CONTACT — la décomposition en temps, version simple
    z_ref, t_ref = rebond(1e-5)
    z_p, t_p = rebond(1e-4, 1e-5, 5e-3)
    print(f"║ bille qui rebondit (k = 1e10, acier) sur 1,2 s :")
    print(f"║   pas fixe 1e-5            z {z_ref:.6f} m   {t_ref:.3f} s")
    print(f"║   piloté 1e-4/1e-5         z {z_p:.6f} m   {t_p:.3f} s   ×{t_ref / t_p:.1f}"
          f"   écart {abs(z_p - z_ref) * 1e6:.1f} µm")
    assert abs(z_p - z_ref) < 1e-6, (z_p, z_ref)
    assert t_p < 0.6 * t_ref, (t_p, t_ref)
    print(f"║   l'estimateur est EXACT et gratuit : on connaît la distance à la"
          f" collision, on resserre avant.")
    # 6. PLAN PORTÉ PAR UN CORPS, et le frottement de Coulomb au complet
    from math import cos, radians, sin, tan
    print("║ bille sur plan incliné porté par un CORPS — roulement et glissement :")
    cas = ((0.30, 15.0), (0.30, 35.0), (0.10, 35.0), (0.02, 15.0)) if rapide else \
          ((0.30, 15.0), (0.30, 35.0), (0.10, 15.0), (0.10, 35.0), (0.02, 15.0), (0.02, 35.0))
    for mu, thd in cas:
        th = radians(thd)
        a = pente(thd, mu)
        roule = mu >= 2.0 / 7.0 * tan(th)
        att = 5.0 / 7.0 * G * sin(th) if roule else G * (sin(th) - mu * cos(th))
        print(f"║   μ {mu:.2f}  θ {thd:4.1f}°   a {a:7.4f}   attendu {att:7.4f}"
              f"   {'roule ' if roule else 'glisse'}   {100 * abs(a / att - 1):5.2f} %")
        assert abs(a / att - 1) < 5e-3, (mu, thd, a, att)

    # 7. BARRIÈRE IPC — l'interpénétration devient impossible
    print("║ barrière IPC contre pénalité de Hertz, bille posée (poids 4,9050 N) :")
    for lbl, kk in (("Hertz k = 1e6", None), ("Hertz k = 1e10", None)):
        kw = dict(k=1e6 if "1e6" in lbl else 1e10, c=10.0)
        N2, b2, _ = bille(**kw)
        N2.simule(1.5, 2e-5, tous=10 ** 9)
        _, dd, ff = N2.contacts()[0]
        print(f"║   {lbl:22} ENFONCE de {dd * 1e6:8.3f} µm   Fn {ff:.4f} N")
        assert dd > 0, dd
    for kk, dh in ((1e3, 1e-3), (1e5, 1e-4)):
        d, f = barriere_ipc(k=kk, d_hat=dh)
        print(f"║   IPC κ = {kk:.0e} d̂ = {dh:.0e}   SÉPARÉE de {-d * 1e6:8.3f} µm"
              f"   Fn {f:.4f} N")
        # LES DEUX FAITS : la bille ne touche JAMAIS, et la force est exacte
        assert d < 0, ("la barrière a laissé pénétrer", d)
        assert abs(f / (m * G) - 1) < 1e-3, (f, m * G)

    print("║ appariement automatique :")
    # 8. TOUT ENSEMBLE, À L'ÉCHELLE : appariement + contacts multiples + IPC
    print("║ N billes lâchées sur un sol, contacts découverts, barrière IPC :")
    tailles = (30, 100) if rapide else (30, 100, 300)
    pts = []
    for n in tailles:
        r = tas(n)
        pts.append((n, r["s"] * 1e3 / r["npas"]))
        print(f"║   {n:4d} billes   {r['paires']:4d} paires simultanées   "
              f"{r['s'] * 1e3 / r['npas']:6.3f} ms/pas   "
              f"enfoncement max {r['d_max'] * 1e6:+8.2f} µm")
        # LE FAIT QUI COMPTE : pas UNE seule paire n'a pénétré
        assert r["d_max"] < 0, ("une paire a pénétré", n, r["d_max"])
    ex = float(np.polyfit(np.log([a for a, _ in pts]), np.log([b for _, b in pts]), 1)[0])
    print(f"║   exposant du temps par pas {ex:.2f} — il suit le nombre de PAIRES, qui")
    print(f"║   croît plus vite que N pendant l'empilement ; l'appariement, lui, est"
          f" mesuré à 1,00")
    # 9. CCD — la garantie topologique : une sphère rapide ne traverse plus
    xa0, xb0, _, _ = tir(ccd=False)
    xa1, xb1, va1, vb1 = tir(ccd=True)
    print(f"║ bille à 50 m/s contre une autre, pas grossier 1e-3 :")
    print(f"║   sans CCD   x_a {xa0:+.4f}  x_b {xb0:+.4f}   → elle a TRAVERSÉ")
    print(f"║   avec CCD   x_a {xa1:+.4f}  x_b {xb1:+.4f}   v_a {va1:+6.2f}  v_b {vb1:+6.2f}"
          f"   → elle rebondit")
    assert xa0 > xb0, ("le cas ne discrimine pas : rien ne traverse", xa0, xb0)
    assert xa1 < xb1, ("le CCD a laissé traverser", xa1, xb1)
    assert vb1 > 0.9 * 50.0, ("le choc n'a pas transmis la vitesse", vb1)
    # 9b. CCD sur boîte, cylindre et maillage — la capsule balayée. La paroi
    # est mince (±0,01 ou plane) et un pas de 0,05 l'enjambe sans la voir.
    print(f"║ bille à 50 m/s contre une paroi mince fixe, pas 1e-3 (face à x = −0,01 / 0) :")
    for forme, face in (("boîte", -0.01), ("cylindre", -0.01), ("maillage", 0.0)):
        x0, x1 = tir_paroi(forme, False), tir_paroi(forme, True)
        print(f"║   {forme:9s} sans CCD x {x0:+.4f} → TRAVERSÉE · avec CCD x {x1:+.4f} → rebond")
        assert x0 > 0.01, ("le cas ne discrimine pas : rien ne traverse", forme, x0)
        assert x1 < face - 0.005, ("le CCD a laissé traverser", forme, x1)
    print("╚═ contact OK")
    return dict(enfoncement=(d, d_th), restitution=dict(zip(map(str, cs), es)))


def _taches(rapide=False):
    from ._campagnes import tache
    cas = [tache(n, 'vinkulum.contact', n, exclusif=n == 'maillages') for n in
           ('multiples', 'capsules', 'capsules_primitives', 'boites', 'cylindres', 'maillages', 'nonlisse')]
    cas.append(tache('penalite_ccd', 'vinkulum.contact', '_penalite_ccd', rapide=rapide, garder=True))
    cas.append(tache('appariement', 'vinkulum.contact', 'appariement', rapide=rapide, exclusif=True))
    return cas


def demo(rapide=False):
    from ._campagnes import executer
    return executer(_taches(rapide), 'contacts')['penalite_ccd']


if __name__ == "__main__":
    demo(rapide="rapide" in sys.argv)
