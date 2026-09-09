"""vinkulum.trim — Newton sur les commandes, à matrice d'influence GARDÉE.

    python -m vinkulum.trim        # autotest (cas algébrique à solution connue)

HOST (Benoit & al., AHS 2000) énonce que le noyau d'un code d'aéromécanique a
TROIS fonctions et que tout le reste en découle : « The functions of HOST are
based on three main utilities: the trim calculation, the time domain
simulation and the calculation of linear equivalent system. » Vinkulum avait
la deuxième et l'analyse modale ; voici la première.

Deux choses seulement, et elles viennent de la même page :

1. **La matrice d'influence se GARDE.** « The convergence scheme iterates with
   the same matrix to reach the solution. As the matrix calculation is time
   consuming, it is only done again when the matrix becomes inaccurate. »
   C'est exactement l'invariant du Newton simplifié de l'intégrateur — un
   jacobien par pas, gardé tant que le résidu chute — appliqué un cran plus
   haut. Ici une matrice coûte n+1 simulations complètes, donc le gain est du
   même ordre que le nombre d'itérations.

2. **La LOI de trim appartient à l'appelant.** « A trim law defines the way in
   which the trim has to be found. It gives the parameters to be imposed […]
   and the parameters to set free ». Ce module ne sait donc RIEN d'un rotor :
   il reçoit une fonction résidu R^n → R^n. Le rotor, le vol en palier, la
   pesée d'un mécanisme ou l'identification d'un paramètre s'écrivent tous
   dedans sans que le trim change d'une ligne.

Ce que ce module N'A PAS et que HOST a : la représentation harmonique de
l'état gérée par le noyau (trim sur les composantes X0, Xic, Xis d'un état
périodique plutôt que sur une moyenne mesurée), et la réduction par isotropie
qui fait tomber le nombre d'inconnues. Déclaré, pas caché.
"""
import numpy as np


def trim(residu, x0, *, pas=1e-3, tol=1e-6, iters=30, seuil_refaire=0.5,
         echelle=None, journal=None, broyden=True):
    """Résout residu(x) = 0 par Newton à matrice d'influence gardée.

    `residu(x) -> array(n)`   ·   `x0` : point de départ (n,)
    `pas`     : perturbation pour la matrice d'influence (scalaire ou (n,))
    `tol`     : sur ‖r‖∞ rapporté à `echelle` (défaut : 1 par composante)
    `seuil_refaire` : la matrice est REFAITE quand une itération ne divise pas
                le résidu par au moins 1/seuil_refaire — c'est la seule
                définition opérationnelle de « the matrix becomes inaccurate ».
    `broyden` : entre deux réfections, la matrice est CORRIGÉE de rang 1 pour
                satisfaire la condition sécante A·δx = δr. Gratuit (le couple
                (δx, δr) est déjà en main), et il rend la convergence
                superlinéaire là où la matrice figée est linéaire.
    Rend (x, r, info) avec info = {evals, matrices, iters, norme, convergé}.

    Sur la mise à jour : le manuel de théorie ABAQUS décrit BFGS avec noyau
    gardé (§2.2.2, d'après Matthies & Strang) — mais BFGS exige un jacobien
    SYMÉTRIQUE et défini positif, ce qu'une matrice d'influence de trim n'est
    ni l'un ni l'autre. La mise à jour qui vaut ici est celle de Broyden, qui
    ne suppose rien de tel. Reprendre BFGS sans le voir aurait donné une
    matrice qui n'a plus de sens.
    """
    x = np.asarray(x0, dtype=float).copy()
    n = x.size
    dx_pert = np.full(n, float(pas)) if np.isscalar(pas) else np.asarray(pas, float)
    ech = np.ones(n) if echelle is None else np.asarray(echelle, float)
    evals = [0]

    def r_de(z):
        evals[0] += 1
        r = np.asarray(residu(z), dtype=float)
        if r.size != n:
            raise ValueError(f"trim : residu rend {r.size} composantes pour {n} commandes")
        return r

    def influence(z, rz):
        # matrice d'influence par perturbation : n simulations de plus.
        # (différences AVANT et non centrées — le coût d'une évaluation est une
        # simulation entière ; le pas se règle, la précision n'est pas le juge,
        # c'est la convergence de Newton qui l'est.)
        a = np.empty((n, n))
        for k in range(n):
            zk = z.copy()
            zk[k] += dx_pert[k]
            a[:, k] = (r_de(zk) - rz) / dx_pert[k]
        return a

    r = r_de(x)
    norme = float(np.max(np.abs(r / ech)))
    a = influence(x, r)
    mats = 1
    bloque = False
    for it in range(iters):
        if norme < tol:
            return x, r, dict(evals=evals[0], matrices=mats, iters=it, norme=norme, converge=True)
        try:
            dx = np.linalg.solve(a, -r)
        except np.linalg.LinAlgError:
            dx = -np.linalg.lstsq(a, r, rcond=None)[0]
        # rebroussement : un pas qui empire le résidu est coupé en deux
        lam, x_n, r_n, n_n = 1.0, None, None, None
        for _ in range(12):
            x_e = x + lam * dx
            r_e = r_de(x_e)
            n_e = float(np.max(np.abs(r_e / ech)))
            if n_e < norme:
                x_n, r_n, n_n = x_e, r_e, n_e
                break
            lam *= 0.5
        if x_n is None:
            # même amorti, ça n'avance plus. Refaire la matrice AU MÊME POINT
            # rendrait la même matrice : on ne la refait qu'une fois, puis on
            # sort en déclarant l'échec. (Sans cette garde, un résidu qui plafonne
            # — plateau saturé, commande sans autorité — consommait une matrice
            # par itération pour ne rien changer : 31 matrices mesurées.)
            if bloque:
                break
            bloque = True
            a = influence(x, r)
            mats += 1
            continue
        bloque = False
        progres = n_n / max(norme, 1e-300)
        if broyden:
            d_x, d_r = x_n - x, r_n - r
            den = float(d_x @ d_x)
            if den > 0.0:
                a = a + np.outer(d_r - a @ d_x, d_x) / den
        x, r, norme = x_n, r_n, n_n
        if journal:
            journal.append(dict(it=it, norme=norme, lam=lam, x=x.copy()))
        if norme >= tol and (progres > seuil_refaire or lam < 1.0):
            a = influence(x, r)               # la matrice a cessé d'être juste
            mats += 1
    return x, r, dict(evals=evals[0], matrices=mats, iters=iters, norme=norme,
                      converge=bool(norme < tol))


def demo():
    """Autotest : cas algébriques à solution connue, et le PRIX de la stratégie.

    Le contrôle qui compte n'est pas « ça converge » — un Newton qui refait sa
    matrice converge aussi, et plus vite en itérations. C'est le COÛT, compté
    en évaluations du résidu, qui est une simulation entière dans un vrai trim.
    """
    rng = np.random.default_rng(7)

    def cas(n, force):
        a = np.eye(n) * 3.0 + 0.3 * rng.standard_normal((n, n))
        cible = rng.standard_normal(n) * (1.0 + 4.0 * force)
        return (lambda x: a @ x + force * x ** 3 - cible)

    # 1. les deux stratégies rendent la MÊME solution
    res = cas(3, 0.15)
    x_g, _, i_g = trim(res, np.zeros(3), pas=1e-4, tol=1e-10, seuil_refaire=0.5)
    x_r, _, i_r = trim(res, np.zeros(3), pas=1e-4, tol=1e-10, seuil_refaire=0.0)
    assert i_g["converge"] and i_r["converge"], (i_g, i_r)
    assert np.max(np.abs(x_g - x_r)) < 1e-7, x_g - x_r
    assert np.max(np.abs(res(x_g))) < 1e-9

    # 2. LE PRIX DES TROIS STRATÉGIES, sur quatre régimes. Le coût est compté
    #    en évaluations du résidu — dans un vrai trim, chacune est une
    #    simulation entière, et la matrice en vaut n.
    strat = (dict(seuil_refaire=0.0, broyden=False),   # refaite à chaque tour
             dict(seuil_refaire=0.5, broyden=False),   # gardée, figée
             dict(seuil_refaire=0.5, broyden=True))    # gardée, corrigée rang 1
    couts = {}
    for n in (3, 12):
        for force in (0.15, 1.2):
            f = cas(n, force)
            out = [trim(f, np.zeros(n), pas=1e-4, tol=1e-10, iters=200, **k)[2]
                   for k in strat]
            assert all(o["converge"] for o in out), out
            couts[(n, force)] = [o["evals"] for o in out]

    # ce que la mesure établit, et rien de plus :
    # · Broyden gagne dans les QUATRE régimes — c'est lui le défaut ;
    assert all(c[2] < min(c[0], c[1]) for c in couts.values()), couts
    # · la matrice figée bat la matrice refaite en non-linéarité MODÉRÉE, et
    #   l'écart s'ouvre avec n (une matrice coûte n évaluations) ;
    assert couts[(3, 0.15)][1] < couts[(3, 0.15)][0]
    assert (couts[(12, 0.15)][0] / couts[(12, 0.15)][1]
            > couts[(3, 0.15)][0] / couts[(3, 0.15)][1]), couts
    # · en non-linéarité FORTE à petit n, elle la perd : la matrice de départ
    #   décrit une autre machine, et les itérations linéaires coûtent plus que
    #   les quadratiques. C'est la raison d'être du seuil de réfection.
    assert couts[(3, 1.2)][1] >= couts[(3, 1.2)][0] * 0.9, couts

    # 4. le rebroussement sauve un départ où le pas de Newton complet sort du
    #    domaine (tanh : la tangente vise 18 fois trop loin)
    plateau = lambda x: np.array([np.tanh(4.0 * x[0]) - 0.9])
    xs, _, i_s = trim(plateau, np.array([1.0]), pas=1e-5, tol=1e-9)
    assert i_s["converge"], i_s
    assert abs(np.tanh(4.0 * xs[0]) - 0.9) < 1e-8

    # 5. et quand la commande n'a PLUS d'autorité (plateau saturé à 1e-6 de
    #    pente), l'échec est DÉCLARÉ, pas dépensé : deux matrices, pas trente
    _, _, i_p = trim(plateau, np.array([2.0]), pas=1e-5, tol=1e-9)
    assert not i_p["converge"], i_p
    assert i_p["matrices"] <= 2, i_p

    # 6. une loi de trim mal posée est REFUSÉE, pas résolue de travers
    try:
        trim(lambda z: np.zeros(4), np.zeros(3))
        raise AssertionError("un résidu de taille != n aurait du lever")
    except ValueError:
        pass

    print(f"  meme solution des deux cotes : ecart {np.max(np.abs(x_g - x_r)):.1e}")
    print("  cout en evaluations du residu :")
    print("      n   non-lin. |  refaite   figee   Broyden")
    for (n, force), c in couts.items():
        print(f"     {n:2d}   {force:5.2f}   |   {c[0]:4d}    {c[1]:4d}     {c[2]:4d}"
              f"    (Broyden x{min(c[0], c[1]) / c[2]:.2f} sur le meilleur des deux autres)")
    print(f"  depart hors domaine (tanh) : converge en {i_s['iters']} iterations par rebroussement")
    print(f"  commande sans autorite (pente 1.8e-6) : echec DECLARE en "
          f"{i_p['evals']} evaluations, {i_p['matrices']} matrices")
    print("trim OK")


if __name__ == "__main__":
    demo()
