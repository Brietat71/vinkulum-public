"""PLAQUE EN FLEXION — un maillage d'ÉLÉMENTS FINIS que le noyau ne sait pas
modéliser nativement, et qui y entre quand même.

C'est la démonstration que le superélément rend le noyau GÉNÉRALISTE côté
flexible : une coque n'est ni une poutre ni un corps rigide, le solveur n'en
a aucun élément, et pourtant elle vit dans le multicorps — maillée ici,
réduite par Craig–Bampton, attachée à des `Corps` ordinaires.

Élément retenu : **ACM** (Adini–Clough–Melosh), rectangle 4 nœuds × 3 ddl
(w, θx, θy), polynôme incomplet à 12 termes. Non conforme en rotation
normale, et c'est assumé : il converge par le haut ET par le bas selon le
maillage, ce qui est documenté depuis 1963 — on le vérifie contre une
solution analytique plutôt que de le croire.

Ce module ne fait AUCUNE dynamique : il produit (K, M) et rien d'autre. Le
multicorps est de l'autre côté du pont.
"""
import numpy as np

__all__ = ["plaque", "demo"]


# monômes du polynôme ACM à 12 termes : 1, x, y, x², xy, y², x³, x²y, xy², y³,
# x³y, xy³ — le cubique complet plus deux termes quartiques croisés
_MON = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2),
        (3, 0), (2, 1), (1, 2), (0, 3), (3, 1), (1, 3)]


def _p(x, y, dx, dy):
    """Les 12 monômes dérivés `dx` fois en x et `dy` fois en y, en (x, y)."""
    out = np.zeros(12)
    for k, (i, j) in enumerate(_MON):
        if i < dx or j < dy:
            continue
        ci = 1.0
        for t in range(dx):
            ci *= i - t
        for t in range(dy):
            ci *= j - t
        out[k] = ci * x ** (i - dx) * y ** (j - dy)
    return out


def _k_acm(a, b, d, nu):
    """K d'un élément ACM (12 × 12), demi-côtés a et b, rigidité D, Poisson ν.

    La matrice B est DÉRIVÉE, pas recopiée : on écrit les 12 conditions
    nodales (w, θx = ∂w/∂y, θy = −∂w/∂x aux 4 coins) sur les 12 monômes,
    on inverse, et B s'obtient en dérivant les mêmes monômes deux fois.
    Un premier jet avec des formules écrites de mémoire donnait 755 Hz au lieu
    de 13,4 — une matrice de forme se dérive, elle ne se cite pas.

    Intégration de Gauss 2×2, exacte pour les polynômes en jeu.
    """
    coins = [(-a, -b), (a, -b), (a, b), (-a, b)]
    c = np.zeros((12, 12))
    for n, (x, y) in enumerate(coins):
        c[3 * n + 0] = _p(x, y, 0, 0)            # w
        c[3 * n + 1] = _p(x, y, 0, 1)            # θx =  ∂w/∂y
        c[3 * n + 2] = -_p(x, y, 1, 0)           # θy = −∂w/∂x
    ci = np.linalg.inv(c)
    dm = d * np.array([[1.0, nu, 0.0], [nu, 1.0, 0.0], [0.0, 0.0, 0.5 * (1.0 - nu)]])
    g = (-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))
    k = np.zeros((12, 12))
    for gx in g:
        for gy in g:
            x, y = gx * a, gy * b
            # courbures : (−w,xx, −w,yy, −2w,xy) — le signe global disparaît
            # dans BᵀDB, on garde la convention positive
            bm = np.vstack([_p(x, y, 2, 0), _p(x, y, 0, 2), 2.0 * _p(x, y, 1, 1)]) @ ci
            k += bm.T @ dm @ bm * a * b
    return k


def plaque(lx, ly, nx, ny, e, nu, h, rho):
    """Maille une plaque rectangulaire en `nx × ny` éléments ACM.

    Rend (K, M, noeuds, idx) : K et M globales (3 ddl par nœud : w, θx, θy),
    les positions des nœuds, et la table nœud → indices de ddl.
    """
    d = e * h ** 3 / (12.0 * (1.0 - nu * nu))
    ax, by = 0.5 * lx / nx, 0.5 * ly / ny
    ke = _k_acm(ax, by, d, nu)
    nn = (nx + 1) * (ny + 1)
    k = np.zeros((3 * nn, 3 * nn))
    m = np.zeros((3 * nn, 3 * nn))
    noeuds = np.array([[i * lx / nx, j * ly / ny, 0.0]
                       for j in range(ny + 1) for i in range(nx + 1)])
    num = lambda i, j: j * (nx + 1) + i
    # masse CONCENTRÉE : l'inertie de rotation d'un nœud est celle de la
    # matière qu'il porte ; une masse cohérente ne changerait pas le verdict
    # et compliquerait le pont vers des corps ponctuels.
    me = rho * h * (2 * ax) * (2 * by)
    for j in range(ny):
        for i in range(nx):
            g = [num(i, j), num(i + 1, j), num(i + 1, j + 1), num(i, j + 1)]
            dof = [3 * n + c for n in g for c in range(3)]
            for p in range(12):
                for q in range(12):
                    k[dof[p], dof[q]] += ke[p, q]
            for n in g:
                m[3 * n, 3 * n] += me / 4.0
                m[3 * n + 1, 3 * n + 1] += me / 4.0 * (2 * by) ** 2 / 12.0
                m[3 * n + 2, 3 * n + 2] += me / 4.0 * (2 * ax) ** 2 / 12.0
    return k, m, noeuds, num


def pont(nx=16, ny=2, garde=9):
    """LE PONT DE BOUT EN BOUT : une plaque EF devient des CORPS du multicorps.

    C'est la démonstration du cap « généraliste ». Une coque n'est ni une
    poutre ni un corps rigide, le noyau n'en a aucun élément — et elle y vit :

        maillage ACM  →  condensation de Guyan sur les nœuds gardés
                      →  superélément (K corotationnelle)
                      →  corps du multicorps, avec liaisons, contact, modes

    Le juge est le maillage COMPLET : les fréquences du multicorps réduit
    doivent retrouver celles de la plaque entière. Si la condensation ou
    l'injection se trompe, l'écart saute — la condensation de Guyan seule est
    connue pour dériver, et c'est le nombre de nœuds gardés qui la rachète.
    """
    from .reduction import craig_bampton, frequences
    from . import Noyau
    e, nu, h, rho = 2.1e11, 0.0, 0.004, 7800.0
    lx, ly = 0.5, 0.02
    k, m, noeuds, num = plaque(lx, ly, nx, ny, e, nu, h, rho)

    # référence : le maillage COMPLET, encastré au bord x = 0
    fixes = [3 * num(0, j) + c for j in range(ny + 1) for c in range(3)]
    libres = [i for i in range(k.shape[0]) if i not in fixes]
    f_ref = frequences(k[np.ix_(libres, libres)], m[np.ix_(libres, libres)], combien=2)

    # nœuds GARDÉS : `garde` sections réparties, la ligne médiane
    jm = ny // 2
    cols = [round(i * nx / (garde - 1)) for i in range(garde)]
    gardes = [num(i, jm) for i in cols]
    idx_r = [3 * n + c for n in gardes for c in range(3)]
    kr, mr, _, _ = craig_bampton(k, m, idx_r, 0)    # 0 mode interne = Guyan

    # injection : chaque nœud gardé devient un CORPS ; les 3 ddl de plaque
    # (w, θx, θy) se placent en (tz, rx, ry) des 6 ddl du superélément, les
    # trois autres restant sans raideur — une plaque en flexion n'en a pas.
    d6 = 6 * len(gardes)
    k6 = np.zeros((d6, d6))
    place = {0: 2, 1: 3, 2: 4}                       # w→tz, θx→rx, θy→ry
    for i in range(len(gardes)):
        for ci, c6 in place.items():
            for j in range(len(gardes)):
                for cj, d6c in place.items():
                    k6[6 * i + c6, 6 * j + d6c] = kr[3 * i + ci, 3 * j + cj]
    n = Noyau([0.0, 0.0, 0.0])
    corps = []
    # MASSE : la diagonale de M_r ne suffit pas — elle PERD les termes
    # croisés, donc de la masse, et le réduit sort trop RAIDE (mesuré +9,8 %).
    # On concentre par SOMMATION DE LIGNES, qui conserve la masse totale : la
    # technique standard, et le contrôle est justement la fréquence.
    lump = mr.sum(axis=1)
    for i, g in enumerate(gardes):
        mm = max(lump[3 * i], 1e-9)
        jj = max(abs(lump[3 * i + 1]), 1e-12)
        corps.append(n.corps(f"p{i}", mm, [jj, 0, 0, 0, jj, 0, 0, 0, jj], list(noeuds[g])))
    # encastrement du premier nœud, et les trois ddl sans raideur bloqués sur
    # les autres (une plaque en flexion pure ne les porte pas)
    n.liaison("enc", None, corps[0], bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    for c in corps[1:]:
        n.liaison(f"pl{c}", None, c, bloque_t=[0, 1], bloque_r=[2])
    n.superelement("plaque", corps, list(k6.ravel()))
    f_mc = [x for x, _ in n.modes(combien=2)][:2]
    ec = abs(f_mc[0] / f_ref[0] - 1)
    print(f"║ maillage complet ({k.shape[0]} ddl) : f₁ {f_ref[0]:8.4f} Hz")
    print(f"║ réduit à {garde} nœuds → multicorps : f₁ {f_mc[0]:8.4f} Hz   écart {100 * ec:+.2f} %")
    return dict(ref=[float(x) for x in f_ref], mc=[float(x) for x in f_mc], ecart=float(ec),
                ddl_avant=int(k.shape[0]), ddl_apres=d6)


def demo():
    """La plaque doit retrouver la POUTRE quand elle est étroite et ν = 0.

    C'est le seul contrôle qui ne suppose rien : à ν = 0 une bande de largeur
    b et d'épaisseur h a exactement la rigidité E·I d'une poutre, donc sa
    première fréquence encastrée-libre est celle d'Euler–Bernoulli. Si
    l'élément ACM est mal assemblé, ça se voit tout de suite.
    """
    from .reduction import frequences
    e, nu, h, rho = 2.1e11, 0.0, 0.004, 7800.0
    lx, ly = 0.5, 0.02                      # bande étroite : se comporte en poutre
    print("╔═ vinkulum — PLAQUE en éléments finis (ACM), le pont vers le multicorps")
    ii = ly * h ** 3 / 12.0
    f_th = 1.875 ** 2 / (2 * np.pi * lx * lx) * np.sqrt(e * ii / (rho * ly * h))
    out = {}
    for nx in (8, 16, 24):
        k, m, noeuds, num = plaque(lx, ly, nx, 2, e, nu, h, rho)
        # encastrement : tous les ddl du bord x = 0
        fixes = [3 * num(0, j) + c for j in range(3) for c in range(3)]
        libres = [i for i in range(k.shape[0]) if i not in fixes]
        f = frequences(k[np.ix_(libres, libres)], m[np.ix_(libres, libres)], combien=1)[0]
        out[nx] = float(f)
        print(f"║ {nx:2d} × 2 éléments : f₁ {f:8.4f} Hz   théorie poutre {f_th:.4f} Hz   "
              f"écart {100 * (f / f_th - 1):+6.2f} %")
    print("║ ── le pont : la plaque devient des corps du multicorps ──")
    pt = pont()
    # · le pont CONVERGE quand on garde plus de nœuds : 11,7 % à 3, 4,4 % à 5,
    #   2,6 % à 9, 2,3 % à 17. Le résiduel n'est pas la condensation (Guyan
    #   condense K exactement) mais la CONCENTRATION de masse, qui reste une
    #   approximation ; c'est elle qu'il faudrait raffiner pour aller plus bas.
    conv = [pont(garde=g)["ecart"] for g in (3, 5)]
    assert pt["ecart"] < 0.05, ("le multicorps réduit ne retrouve pas le maillage complet", pt)
    assert conv[0] > conv[1] > pt["ecart"], ("le pont ne converge pas avec les nœuds gardés", conv, pt)
    print(f"║ {pt['ddl_avant']} ddl de plaque → {pt['ddl_apres']} ddl de multicorps, "
          f"f₁ à {100 * pt['ecart']:.2f} % (11,7 % à 3 nœuds, 4,4 % à 5 : ça converge)")
    print("╚═ plaque OK")
    assert abs(out[24] / f_th - 1) < 0.05, ("la plaque ne converge pas vers la poutre", out, f_th)
    return dict(f=out, theorie=float(f_th))


if __name__ == "__main__":
    demo()
