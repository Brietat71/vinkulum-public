"""vinkulum.reduction — Craig–Bampton : un corps flexible réduit à ses interfaces.

    python -m vinkulum.reduction        # console encastrée, réduction et son prix

Une sous-structure flexible maillée finement coûte cher à intégrer alors que
son mouvement utile tient dans peu de coordonnées. Craig & Bampton (1968) la
remplacent par ses degrés d'INTERFACE plus quelques coordonnées modales :

    u_E = Ψ · u_R + Φ · q,      Ψ = −K_EE⁻¹ K_ER      (modes de CONTRAINTE)
                                Φ  = modes normaux à interfaces BLOQUÉES

    K_r = [[K_RR + K_RE Ψ, 0], [0, Λ]]        M_r = Tᵀ M T,   T = [[I, 0], [Ψ, Φ]]

ABAQUS le décrit en §2.14.1 (« Substructuring and superelement analysis ») :
la partie Ψ seule est la réduction de Guyan, celle que NASTRAN-95 fait déjà en
1972 ; ce qui la rend dynamiquement juste est l'ajout de Φ, « an additional,
and generally more effective, technique ».

CE QUE CETTE RÉDUCTION N'EST PAS, et il faut le dire avant de s'en servir :
elle est LINÉAIRE autour d'un état. Le manuel ABAQUS l'écrit — « the response
within a substructure […] is considered to be a linear perturbation about the
state of the substructure at the time it is made into a superelement ». Un
superélément **ne tourne pas** : réduire puis mettre en rotation perd le
raidissement centrifuge, qui dépend du régime. C'est pourquoi une PALE garde
ses poutres géométriquement exactes et pourquoi Craig–Bampton vise ici la
CELLULE raide, pas le rotor.

Ce module fait la RÉDUCTION et la juge ; il ne l'intègre pas au solveur comme
un élément — ce serait des degrés de liberté généralisés dans la formulation,
et ça ne se fait qu'une fois qu'on sait ce que la réduction vaut.
"""
import sys
import time

import numpy as np


def craig_bampton(k, m, idx_r, n_modes):
    """Réduit (K, M) en gardant les DDL `idx_r` et `n_modes` modes internes.

    Rend (K_r, M_r, T, freqs_bloquees) où T est la transformation
    (n × (len(idx_r) + n_modes)) qui ramène l'espace complet.
    """
    k = 0.5 * (np.asarray(k, float) + np.asarray(k, float).T)
    m = 0.5 * (np.asarray(m, float) + np.asarray(m, float).T)
    n = k.shape[0]
    r = list(idx_r)
    e = [i for i in range(n) if i not in set(r)]
    if not e:
        raise ValueError("rien à réduire : tous les degrés sont retenus")
    k_rr, k_re = k[np.ix_(r, r)], k[np.ix_(r, e)]
    k_ee = k[np.ix_(e, e)]
    # modes de CONTRAINTE : la déformée statique quand on impose 1 sur un DDL
    # d'interface, les autres bloqués. C'est Guyan, et c'est exact en statique.
    psi = -np.linalg.solve(k_ee, k_re.T)
    # modes NORMAUX à interfaces bloquées : ce que Guyan ignore, et ce qui rend
    # la réduction juste en dynamique
    from scipy.linalg import eigh
    lam, phi = eigh(k_ee, m[np.ix_(e, e)])
    lam, phi = lam[:n_modes], phi[:, :n_modes]
    # normalisation en masse (φᵀ M φ = I), ce que `eigh` généralisé donne déjà
    t = np.zeros((n, len(r) + n_modes))
    for j, i in enumerate(r):
        t[i, j] = 1.0
    t[np.ix_(e, range(len(r)))] = psi
    t[np.ix_(e, range(len(r), len(r) + n_modes))] = phi
    k_red = t.T @ k @ t
    m_red = t.T @ m @ t
    # le bloc de couplage de K est nul par construction (Ψ annule K_ER) : on le
    # SYMÉTRISE et on l'ASSERTE plutôt que de le supposer
    if n_modes:
        couplage = np.abs(k_red[:len(r), len(r):]).max()
        ech = max(abs(k_red).max(), 1e-30)
        assert couplage < 1e-8 * ech, f"couplage K_Rq non nul : {couplage:.2e} pour ‖K‖ {ech:.2e}"
    f_bloq = np.sqrt(np.maximum(lam, 0.0)) / (2 * np.pi)
    return 0.5 * (k_red + k_red.T), 0.5 * (m_red + m_red.T), t, f_bloq


def frequences(k, m, combien=6):
    """Fréquences propres (Hz) d'un couple (K, M) symétrique, M définie positive."""
    from scipy.linalg import eigh
    lam = eigh(np.asarray(k, float), np.asarray(m, float), eigvals_only=True)
    lam = lam[lam > 1e-9]
    return np.sqrt(lam[:combien]) / (2 * np.pi)


def _console(n_elem, L=1.0, b=0.02, hh=0.004, e_mod=70e9, rho=2700.0):
    """Console encastrée en poutres géométriquement exactes, comme `bancs`."""
    from vinkulum import Noyau
    a = b * hh
    i2 = b * hh ** 3 / 12.0
    i3 = hh * b ** 3 / 12.0
    jt = i2 + i3
    g = e_mod / (2 * 1.3)
    dl = L / n_elem
    mel = rho * a * dl
    N = Noyau([0.0, 0.0, 0.0])
    idx = []
    for i in range(n_elem + 1):
        mm = mel * (0.5 if i in (0, n_elem) else 1.0)
        j = np.diag([rho * jt * dl, rho * i2 * dl, rho * i3 * dl])
        idx.append(N.corps(f"n{i}", max(mm, 1e-9), list(j.ravel()), [i * dl, 0.0, 0.0]))
    N.liaison("encastre", None, idx[0])
    for i in range(n_elem):
        N.poutre(f"b{i}", idx[i], idx[i + 1], e_mod * a, g * a, g * jt, e_mod * i2)
    return N, idx


def demo(rapide=False):
    print("╔═ vinkulum — Craig–Bampton : ce que la réduction garde et ce qu'elle coûte")
    n_elem = 10 if rapide else 20
    N, idx = _console(n_elem)
    k, m, z = N.k_m_z()
    k, m, z = np.array(k), np.array(m), np.array(z)
    # on travaille dans la base ADMISSIBLE : les contraintes (l'encastrement)
    # y sont déjà éliminées, sinon K_EE est singulière
    ka, ma = z.T @ k @ z, z.T @ m @ z
    ka = 0.5 * (ka + ka.T)
    ma = 0.5 * (ma + ma.T)
    f_ref = frequences(ka, ma, 4)
    print(f"║ console {n_elem} éléments · {ka.shape[0]} degrés admissibles")
    print(f"║ complet   f = " + "  ".join(f"{x:8.3f}" for x in f_ref) + "  Hz")

    # INTERFACE : les degrés du nœud du BOUT, exprimés dans la base admissible.
    # Z est orthonormale, donc la projection d'un degré physique j est sa ligne.
    # On prend les p directions de Z qui portent le plus le nœud du bout.
    bout = list(range(6 * idx[-1], 6 * idx[-1] + 6))
    poids = np.linalg.norm(z[bout, :], axis=0)
    ordre = np.argsort(-poids)
    idx_r = sorted(int(i) for i in ordre[:6])
    print(f"║ interface : les 6 directions qui portent le nœud du bout "
          f"(poids {poids[ordre[:6]].min():.3f} à {poids[ordre[0]]:.3f})")

    t0 = time.time()
    lignes = []
    for nm in (0, 1, 2, 3, 5, 8):
        kr, mr, tr, f_bl = craig_bampton(ka, ma, idx_r, nm)
        f = frequences(kr, mr, 4)
        err = [abs(f[i] / f_ref[i] - 1) for i in range(min(len(f), len(f_ref)))]
        lignes.append((nm, kr.shape[0], f, max(err)))
        print(f"║ {('Guyan seul' if nm == 0 else f'+ {nm} modes'):12} "
              f"{kr.shape[0]:3d} ddl   f = " + "  ".join(f"{x:8.3f}" for x in f)
              + f"   écart max {100 * max(err):6.2f} %")
    print(f"║ {time.time() - t0:.2f} s")

    # CE QUE LA MESURE DOIT ÉTABLIR, et qui peut échouer :
    # · Guyan seul est MAUVAIS en dynamique — c'est sa limite connue ;
    guyan = lignes[0][3]
    assert guyan > 0.05, f"Guyan seul à {100 * guyan:.2f} % : le cas ne discrimine rien"
    # · l'ajout de modes internes améliore, et de façon MONOTONE ;
    errs = [l[3] for l in lignes]
    assert all(errs[i + 1] <= errs[i] * 1.001 for i in range(len(errs) - 1)), errs
    # · quelques modes suffisent : le mécanisme n'a d'intérêt que s'il CONVERGE
    assert lignes[-1][3] < 1e-3, lignes[-1]
    n_util = next(l[0] for l in lignes if l[3] < 1e-2)
    n_ddl = next(l[1] for l in lignes if l[3] < 1e-2)
    print(f"║ VERDICT : {n_ddl} degrés au lieu de {ka.shape[0]} (×{ka.shape[0] / n_ddl:.1f}) "
          f"pour 1 % sur les quatre premières fréquences — {n_util} modes internes")
    print(f"║ Guyan SEUL : {100 * guyan:.1f} % d'écart. C'est ce que NASTRAN-95 fait "
          f"en 1972, et c'est ce que les modes internes réparent.")
    print("╚═ réduction OK")
    return dict(f_ref=[float(x) for x in f_ref], lignes=[(l[0], l[1], float(l[3])) for l in lignes])


if __name__ == "__main__":
    demo(rapide="rapide" in sys.argv)
