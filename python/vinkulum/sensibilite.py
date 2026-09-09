"""vinkulum.sensibilite — le gradient par l'ADJOINT : une résolution pour tous
les paramètres.

    python -m vinkulum.sensibilite      # console : gradient adjoint contre DF

L'AXE DÉCISIF DU PROJET, et il se juge sur un fait simple : combien coûte le
gradient d'une sortie par rapport à N paramètres.

    différences finies globales   N + 1 résolutions complètes, 3 chiffres
    mode DIRECT (ABAQUS/Design)   N résolutions du système tangent
    mode ADJOINT (ici)            UNE résolution, quel que soit N

À l'équilibre R(q, p) = 0 et pour un objectif J(q, p) :

    dJ/dp = ∂J/∂p − ψᵀ · ∂R/∂p       avec     Kᵀ ψ = (∂J/∂q)ᵀ

ψ ne dépend PAS de p : on le résout une fois, puis chaque paramètre coûte un
produit scalaire. C'est ce qui rend l'optimisation de conception praticable à
cent paramètres, là où les différences finies en demandent cent une.

CE QUE FAIT LE CONCURRENT, écrit dans son propre manuel (ABAQUS Theory
Manual §2.17.1) : `ABAQUS/Design` couvre « nonperturbation, **static** stress
problems », « only solid elements with elastic or hyperelastic properties »,
et sa méthode est **semi-analytique** — les vecteurs élémentaires viennent de
différences finies, dont le manuel nomme lui-même le prix (« if the interval
is too small, round-off or cancellation errors occur […] if too large,
truncation errors »). Et son DSA est en mode DIRECT : une résolution par
paramètre.

CE QUE CE MODULE FAIT, ET CE QU'IL NE FAIT PAS. Il est en mode ADJOINT, donc
indépendant de N. Le jacobien K vient de la raideur TANGENTE du noyau, exacte
par duaux. `∂R/∂p`, lui, est encore obtenu par différences finies sur le
paramètre — semi-analytique, comme ABAQUS : chaque paramètre coûte une
ASSEMBLAGE, pas une résolution. Rendre ce terme exact demande que les
paramètres traversent les éléments en duaux, et c'est le chantier suivant ;
il est nommé, pas caché. La sensibilité le long d'une TRAJECTOIRE (adjoint en
temps) n'est pas non plus ici.
"""
import sys
import time

import numpy as np


def equilibre(construit, t_end=5.0, h=2e-3, tol_v=1e-6):
    """Amène le modèle à son équilibre statique par relaxation amortie.

    ⚠ **CE DOCSTRING A DIT « le noyau n'a pas de solveur statique » jusqu'au
    5 septembre. C'EST FAUX** — `Noyau.statique` existe (Newton sur R(q) = 0,
    système augmenté [[K, Gᵀ], [G, 0]], équilibrage de Jacobi, Newton amorti),
    et il a été écrit APRÈS cette ligne sans que personne ne revienne la
    corriger. Une phrase qui annonce une lacune comblée envoie le lecteur
    écrire ce qui existe déjà : c'est le genre de trace qui coûte un chantier.

    La relaxation reste ici, et c'est désormais un CHOIX et non une nécessité :
    elle marche sur ce corpus, elle est vérifiée par son propre contrôle de
    vitesse résiduelle, et la remplacer demanderait de rejouer les gradients.
    Ce qui n'est PAS mesuré, et qu'on n'affirmera donc pas : laquelle des deux
    est la plus robuste ici. La brique `statique` est là, c'est l'appel qui ne
    l'est pas.
    """
    N = construit()
    N.simule(t_end, h, tous=10 ** 9)
    v = max(np.abs(np.array(N.etat()[3], float)).max(),
            np.abs(np.array(N.etat()[4], float)).max())
    # L'ÉQUILIBRE DOIT ÊTRE ATTEINT, et ça se VÉRIFIE : à t_end trop court la
    # flèche vaut 2,89 mm au lieu de 4,60 et le gradient adjoint sort 10× trop
    # petit — sans que rien ne crie, parce que la DF globale, elle, refait la
    # relaxation à chaque point et voit la même erreur des deux côtés.
    if v > tol_v:
        raise ValueError(f"équilibre non atteint : vitesse résiduelle {v:.2e} > {tol_v:.0e} "
                         f"après {t_end} s — allonger t_end ou amortir plus")
    return N, v


def adjoint(construit, objectif, params, eps=None, t_end=5.0, h=2e-3, exact=None):
    """Gradient de `objectif` par rapport à `params`, par l'adjoint.

    `construit(vals)` -> Noyau, où `vals` est le vecteur des paramètres.
    `objectif(N)` -> (valeur, dJ/dq) avec dJ/dq sur les 6n degrés physiques.
    `params` : valeurs nominales.
    `exact` : liste de (indice de poutre, "ea"|"ga"|"gj"|"ei") — un par
    paramètre. Fournie, ∂R/∂p est pris EXACT par linéarité de l'énergie de
    poutre (aucune différence finie) ; absente, il est semi-analytique.
    """
    p0 = np.asarray(params, float)
    e = eps if eps is not None else np.maximum(np.abs(p0), 1.0) * 1e-6

    N, v_res = equilibre(lambda: construit(p0), t_end, h)
    j0, dj_dq = objectif(N)
    k, c, m, z, g = (np.array(x, float) for x in N.k_c_m_z())

    # L'ÉQUILIBRE VIT DANS LE NOYAU DES CONTRAINTES : le système tangent y est
    # Z'KZ, et l'adjoint aussi. Hors de là, K est singulière (les modes de
    # liaison), et une résolution directe rendrait n'importe quoi.
    kz = z.T @ k @ z
    psi_z = np.linalg.solve(kz.T, z.T @ dj_dq)

    # ∂R/∂p : semi-analytique. On reconstruit le modèle au paramètre perturbé
    # SANS le rééquilibrer, et on lit la variation du résidu à q figé — c'est
    # ce que ABAQUS appelle « holding the incremental displacement constant ».
    etat = N.etat()
    grad = np.zeros(len(p0))
    n_ass = 0
    for i in range(len(p0)):
        if exact is not None:
            # EXACT : coefficients condensés dérivés analytiquement, puis
            # gradient de l'énergie dérivée par duaux, sans pas à régler.
            j, quoi = exact[i]
            dr_dp = np.array(N.d_residu_poutre(j, quoi), float)
            n_ass += 1
        else:
            dr = []
            for sens in (+1.0, -1.0):
                pp = p0.copy()
                pp[i] += sens * e[i]
                Np = construit(pp)
                Np.pose_etat(etat[0], [list(x) for x in etat[1]], [list(x) for x in etat[2]],
                             [list(x) for x in etat[3]], [list(x) for x in etat[4]], list(etat[5]))
                dr.append(np.array(Np.residu_statique(), float))
            dr_dp = (dr[0] - dr[1]) / (2 * e[i])
            n_ass += 2
        grad[i] = -psi_z @ (z.T @ dr_dp)
    return j0, grad, dict(residu_vitesse=v_res, n_resolutions=1, n_assemblages=n_ass,
                          exact=exact is not None)


def gradient_df(construit, objectif_val, params, eps=None, t_end=5.0, h=2e-3):
    """Gradient par différences finies GLOBALES : 2N équilibres complets."""
    p0 = np.asarray(params, float)
    e = eps if eps is not None else np.maximum(np.abs(p0), 1.0) * 1e-4
    g = np.zeros(len(p0))
    for i in range(len(p0)):
        vals = []
        for sens in (+1.0, -1.0):
            pp = p0.copy()
            pp[i] += sens * e[i]
            N, _ = equilibre(lambda: construit(pp), t_end, h)
            vals.append(objectif_val(N))
        g[i] = (vals[0] - vals[1]) / (2 * e[i])
    return g


def gradient_frequence(N, mode, exact):
    """dλ/dp et df/dp EXACTS pour un mode propre, sans re-résolution.

    Pour φ normalisé en masse (φᵀMφ = 1) et λ = ω² simple :

        dλ/dp = φᵀ (∂K/∂p − λ ∂M/∂p) φ

    C'est la formule de perturbation d'une valeur propre simple. Elle est
    EXACTE, elle ne demande ni pas ni résolution supplémentaire, et le seul
    terme qui manquait — ∂K/∂p — est lui aussi exact par linéarité de
    l'énergie de poutre. La sensibilité d'une fréquence à N paramètres coûte
    donc N produits, contre N+1 extractions modales en différences finies.

    `mode` : rang du mode (0 = premier). `exact` : liste de (poutre, "ei"…).
    Rend (f_Hz, [df/dp]).
    """
    from scipy.linalg import eigh
    k, c, m, z, g = (np.array(x, float) for x in N.k_c_m_z())
    kz = 0.5 * (z.T @ k @ z + (z.T @ k @ z).T)
    mz = 0.5 * (z.T @ m @ z + (z.T @ m @ z).T)
    lam, vec = eigh(kz, mz)
    garde = lam > 1e-9
    lam, vec = lam[garde], vec[:, garde]
    l0 = lam[mode]
    phi = z @ vec[:, mode]                 # ramené en degrés physiques
    d = np.zeros(len(exact))
    for i, (j, quoi) in enumerate(exact):
        dk = np.array(N.d_raideur_poutre(j, quoi), float)
        dk = 0.5 * (dk + dk.T)
        # ∂M/∂p = 0 : une raideur ne porte pas de masse
        d[i] = float(phi @ dk @ phi)
    f = np.sqrt(l0) / (2 * np.pi)
    # f = √λ/2π ⇒ df/dλ = 1/(4π√λ) = 1/(8π²f). Le contrôle contre la DF a
    # sorti un facteur 2π exactement : la conversion se vérifie, elle ne se
    # pose pas.
    return f, d / (8 * np.pi ** 2 * f)


# ── cas de contrôle : console sous son propre poids ──────────────────────────
def _console(ei_par_elem, n_elem=6, L=0.6, b=0.02, hh=0.004, rho=2700.0, e_mod=70e9):
    from vinkulum import Noyau
    a = b * hh
    i3 = hh * b ** 3 / 12.0
    jt = b * hh ** 3 / 12.0 + i3
    g_mod = e_mod / (2 * 1.3)
    dl = L / n_elem
    mel = rho * a * dl
    N = Noyau([0.0, 0.0, -9.81])
    idx = []
    for i in range(n_elem + 1):
        mm = mel * (0.5 if i in (0, n_elem) else 1.0)
        j = np.diag([rho * jt * dl, rho * (b * hh ** 3 / 12) * dl, rho * i3 * dl])
        idx.append(N.corps(f"n{i}", max(mm, 1e-9), list(j.ravel()), [i * dl, 0.0, 0.0]))
    N.liaison("encastre", None, idx[0])
    for i in range(n_elem):
        # amortissement visqueux global : la relaxation vers l'équilibre
        N.poutre(f"b{i}", idx[i], idx[i + 1], e_mod * a, g_mod * a, g_mod * jt,
                 float(ei_par_elem[i]))
    for i in idx[1:]:
        N.couple(f"amo{i}", None, i, [0.0, 1.0, 0.0], ("ressort", [0.0, 4e-2, 0.0]))
    return N, idx


def demo(rapide=False):
    print("╔═ vinkulum — gradient par l'ADJOINT : une résolution pour tous les paramètres")
    n_elem = 4 if rapide else 6
    e_mod, b, hh = 70e9, 0.02, 0.004
    ei0 = np.full(n_elem, e_mod * b * hh ** 3 / 12.0)

    def build(p):
        return _console(p, n_elem=n_elem)[0]

    def obj(N):
        """flèche verticale du nœud du bout, et son gradient en q."""
        r = np.array(N.etat()[1], float)
        nb = r.shape[0]
        d = np.zeros(6 * nb)
        d[6 * (nb - 1) + 2] = 1.0
        return float(r[-1, 2]), d

    def obj_val(N):
        return float(np.array(N.etat()[1], float)[-1, 2])

    t0 = time.time()
    j0, g_adj, info = adjoint(build, obj, ei0,
                              exact=[(i, "ei") for i in range(n_elem)])
    t_adj = time.time() - t0
    _, g_semi, info_semi = adjoint(build, obj, ei0)     # ∂R/∂p par DF
    t0 = time.time()
    g_df = gradient_df(build, obj_val, ei0)
    t_df = time.time() - t0

    print(f"║ console {n_elem} éléments · flèche à l'équilibre {1e3 * j0:.4f} mm "
          f"(vitesse résiduelle {info['residu_vitesse']:.1e})")
    print(f"║ paramètres : les {n_elem} raideurs de flexion EI, nominal {ei0[0]:.4f} N·m²")
    print("║   i   adjoint EXACT (m/(N·m²))   DF globale        écart")
    for i in range(n_elem):
        ec = abs(g_adj[i] - g_df[i]) / max(abs(g_df[i]), 1e-30)
        print(f"║  {i:2d}    {g_adj[i]:+.6e}       {g_df[i]:+.6e}   {100 * ec:6.2f} %")
    print(f"║ coût : adjoint {info['n_resolutions']} résolution + "
          f"{info['n_assemblages']} assemblages, {t_adj:.2f} s   ·   "
          f"DF globale {2 * n_elem} équilibres complets, {t_df:.2f} s   "
          f"(×{t_df / max(t_adj, 1e-9):.1f})")

    ec = np.abs(g_adj - g_df) / np.maximum(np.abs(g_df), 1e-30)
    ec_semi = np.abs(g_semi - g_df) / np.maximum(np.abs(g_df), 1e-30)
    print(f"║ ∂R/∂p EXACT (règle de chaîne et duaux) : "
          f"{info['n_assemblages']} évaluations, écart max {100 * np.max(ec):.3f} %")
    print(f"║ ∂R/∂p semi-analytique (comme ABAQUS)           : "
          f"{info_semi['n_assemblages']} évaluations, écart max {100 * np.max(ec_semi):.3f} %")
    # LE POINT QUI DÉCIDE, et il se mesure : la dérivée semi-analytique a un
    # PAS à régler, et le manuel ABAQUS nomme lui-même les deux murs — « if the
    # interval is too small, round-off or cancellation errors occur […] if too
    # large, truncation errors ». L'exact n'a pas de pas du tout.
    print("║ balayage du pas de la différence finie sur ∂R/∂p :")
    pire = 0.0
    # LE MUR D'ANNULATION COMMENCE À 1e-13 sur 4 éléments (0,04 % à 1e-12,
    # 3 % à 1e-13, 14 % à 1e-14 — mesuré le 5 sept.) ; un balayage arrêté à
    # 1e-12 le voyait à peine et l'assert « le semi-analytique se dégrade »
    # tombait en mode rapide (rapport 4,4 pour 10). On va jusqu'au mur.
    for rel in (1e-2, 1e-4, 1e-6, 1e-9, 1e-12, 1e-13, 1e-14):
        _, g_r, _ = adjoint(build, obj, ei0, eps=np.abs(ei0) * rel)
        e_r = float(np.max(np.abs(g_r - g_df) / np.maximum(np.abs(g_df), 1e-30)))
        pire = max(pire, e_r)
        print(f"║      pas relatif {rel:.0e} : écart max {100 * e_r:10.4f} %")
    print(f"║      EXACT (aucun pas)       : écart max {100 * np.max(ec):10.4f} %  — invariant")
    print("║      Les très petits pas exposent l'annulation ; les grands pas peuvent")
    print("║      aussi tronquer la dépendance non linéaire des coefficients condensés.")
    # le semi-analytique DOIT se dégrader quelque part sur ce balayage, sinon
    # le cas ne discrimine rien et l'argument ne vaut pas
    assert pire > 10 * max(np.max(ec), 1e-12), (pire, np.max(ec))
    assert np.max(ec) < 0.05, (g_adj, g_df, ec)
    # l'exact coûte deux fois moins d'évaluations que le semi-analytique, et il
    # n'a ni pas à régler ni annulation — c'est le fait qu'on garde
    assert info["n_assemblages"] * 2 == info_semi["n_assemblages"], (info, info_semi)
    assert np.max(ec) <= np.max(ec_semi) * 1.5, (ec, ec_semi)
    # LE FAIT QUI COMPTE : le coût de l'adjoint ne suit PAS le nombre de
    # paramètres, celui des DF si. Une seule résolution, quel que soit N.
    assert info["n_resolutions"] == 1, info
    # et le signe est physique : raidir un élément REMONTE le bout
    assert np.all(g_adj > 0), g_adj
    # le premier élément pèse le plus : le moment y est maximal
    assert g_adj[0] == np.max(g_adj), g_adj
    print(f"║ le gradient est physique : raidir remonte le bout, et l'élément "
          f"d'ENCASTREMENT pèse ×{g_adj[0] / g_adj[-1]:.1f} le dernier")

    # ── sensibilité d'une FRÉQUENCE PROPRE, exacte et sans re-résolution ──
    exact = [(i, "ei") for i in range(n_elem)]
    N_eq, _ = equilibre(lambda: build(ei0), t_end=5.0, h=2e-3)
    f0, df = gradient_frequence(N_eq, 0, exact)
    # Pas de DF 1e-2 et non 1e-4 : c'est la RÉFÉRENCE qui est bruitée, pas
    # l'analytique. La fréquence sort d'une extraction modale, donc à un
    # plancher d'arrondi près ; divisée par 2δp elle rend un bruit en 1/δp,
    # qui écrase le signal des éléments de bout (leur influence sur f₁ est
    # mille fois moindre que celle de l'encastrement). Mesuré à 6 éléments,
    # écart de l'élément de bout : 57 % à δ 1e-4 · 7,6 % à 1e-3 · 0,53 % à
    # 1e-2 · 1,1 % à 1e-1. Un écart qui DÉCROÎT quand le pas grandit est la
    # signature du bruit ; la troncature ferait l'inverse, et c'est elle qui
    # remonte à 1e-1. On se pose au creux.
    d_rel = 1e-2
    df_num = np.zeros(n_elem)
    for i in range(n_elem):
        vals = []
        for sens in (+1.0, -1.0):
            pp = ei0.copy()
            pp[i] += sens * ei0[i] * d_rel
            Np, _ = equilibre(lambda: build(pp), t_end=5.0, h=2e-3)
            vals.append(gradient_frequence(Np, 0, exact)[0])
        df_num[i] = (vals[0] - vals[1]) / (2 * ei0[i] * d_rel)
    ecf = np.abs(df - df_num) / np.maximum(np.abs(df_num), 1e-30)
    print(f"║ FRÉQUENCE du 1er mode {f0:.5f} Hz — sa sensibilité aux mêmes EI :")
    for i in range(n_elem):
        print(f"║  {i:2d}    {df[i]:+.6e} Hz/(N·m²)   contre DF {df_num[i]:+.6e}"
              f"   {100 * ecf[i]:6.2f} %")
    print(f"║ formule de perturbation : {n_elem} produits, AUCUNE extraction modale "
          f"de plus ; la DF en demande {2 * n_elem}")
    assert np.max(ecf) < 0.01, (df, df_num, ecf)
    assert np.all(df > 0), df
    print("╚═ sensibilité OK")
    return dict(grad=g_adj.tolist(), ecart=float(np.max(ec)))


if __name__ == "__main__":
    demo(rapide="rapide" in sys.argv)
