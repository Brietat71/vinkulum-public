"""Pilote de la routine HCB officielle Exudyn 1.11, sur matrices fournies.

Cette adaptation est ALGEBRIQUE : chaque corps à six coordonnées devient
deux groupes de trois indices Position. Il ne s'agit pas d'un maillage solide
et ce modèle ne doit pas être exporté vers ObjectFFRFreducedOrder/RBE2/RBE3.
Le mode allBoundaryNodes conserve chacune des six coordonnées terminales.

Documentation uniquement, sans lecture de FEM.py :
https://exudyn.readthedocs.io/en/v1.11.0/docs/RST/pythonUtilities/FEM.html
https://jgerstmayr.github.io/EXUDYN/docs/RST/ModelOrderReductionAndComponentModeSynthesis.html
"""
import argparse
import json
import operator
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh, solve
from scipy.sparse import csr_matrix, diags


def prepare(D, M, nbModes, *, projection="standard", positions=None):
    """Retourne base physique, projections et diagnostics HCB officiels.

    Les six dernières coordonnées sont les ports physiques. Construction de
    K=D.T@D, équilibrage symétrique S=diag(K)^(-1/2), appel officiel HCB,
    puis Tphys=S@basis. Projection standard : Tphys.T K Tphys. Variante
    énergie : (D Tphys).T (D Tphys), avec exactement la même routine HCB.
    La masse est toujours projetée par Tphys.T M Tphys.
    Une diagonalisation équilibrée de la petite paire est également incluse
    dans la préparation. Les diagnostics de
    référence et l'évaluation des réponses sont séparés par l'appelant.
    """
    import exudyn
    from exudyn.FEM import FEMinterface, HCBstaticModeSelection
    debut = time.perf_counter()
    if isinstance(nbModes, (bool, np.bool_)) or operator.index(nbModes) < 0:
        raise ValueError("nombre de modes entier non négatif requis")
    nbModes = operator.index(nbModes)
    if projection not in ("standard", "energie"):
        raise ValueError("projection standard ou energie requise")
    d, m = csr_matrix(D, dtype=float), csr_matrix(M, dtype=float)
    n = d.shape[1]
    if n < 12 or n % 6 or m.shape != (n, n) or nbModes >= n-6:
        raise ValueError("6 coordonnées par corps, six ports terminaux, nbModes < dimension intérieure requis")
    if not np.all(np.isfinite(d.data)) or not np.all(np.isfinite(m.data)):
        raise ValueError("D et M finis requis")
    k = (d.T @ d).tocsr()
    diag = k.diagonal()
    if np.any(diag <= 0):
        raise ValueError("diagonale de K strictement positive requise")
    echelles = 1/np.sqrt(diag)
    s = diags(echelles)
    khat, mhat = (s @ k @ s).tocsr(), (s @ m @ s).tocsr()
    preparation_entree = time.perf_counter()-debut
    fem = FEMinterface()
    xy = np.zeros((n//3, 3)) if positions is None else np.asarray(positions, dtype=float)
    if xy.shape != (n//3, 3):
        raise ValueError("positions auxiliaires de forme (nombre de triplets,3) requises")
    fem.nodes = {"Position": xy.copy()}
    fem.massMatrix, fem.stiffnessMatrix = mhat, khat
    # Ni RBE, ni suppression des six coordonnées terminales comme mouvement
    # rigide : l'encastrement est déjà supprimé des matrices physiques.
    avant_hcb = time.perf_counter()
    fem.ComputeHurtyCraigBamptonModes(
        boundaryNodesList=[[n//3-2, n//3-1]], nEigenModes=nbModes,
        useSparseSolver=True, computationMode=HCBstaticModeSelection.allBoundaryNodes,
        excludeRigidBodyMotion=False, verboseMode=False)
    apres_hcb = time.perf_counter()
    base = np.asarray(fem.modeBasis["matrix"])
    if base.shape != (n, 6+nbModes):
        raise ArithmeticError("dimension HCB incompatible avec les six ports conservés")
    trace_statique = float(np.linalg.norm(base[-6:, :6]-np.eye(6)))
    trace_modes = float(np.linalg.norm(base[-6:, 6:]))
    if trace_statique > 1e-12 or trace_modes > 1e-12:
        raise ArithmeticError("la trace HCB ne conserve pas les six ports indépendants")
    tphys = echelles[:, None]*base
    if projection == "energie":
        dt = d @ tphys
        kr_brut = dt.T @ dt
    else:
        kr_brut = tphys.T @ (k @ tphys)
    mr_brut = tphys.T @ (m @ tphys)
    kr, mr = .5*(kr_brut+kr_brut.T), .5*(mr_brut+mr_brut.T)
    if kr.shape != (6+nbModes, 6+nbModes) or mr.shape != kr.shape:
        raise ArithmeticError("projections réduites de dimensions invalides")
    termine_projection = time.perf_counter()
    eq_reduit = 1/np.sqrt(np.diag(kr))
    valeurs_reduites, vecteurs = eigh(eq_reduit[:, None]*kr*eq_reduit[None, :],
                                     eq_reduit[:, None]*mr*eq_reduit[None, :])
    q = eq_reduit[:, None]*vecteurs
    projection_force = q.T @ tphys[-6:].T
    termine_spectre = time.perf_counter()
    # Diagnostics hors coût de préparation, explicitement chronométrés.
    if projection == "energie":
        energie = kr_brut
        standard = tphys.T @ (k @ tphys)
        standard = .5*(standard+standard.T)
    else:
        dt = d @ tphys
        energie = dt.T @ dt
        standard = kr
    statique = khat[:-6] @ base[:, :6]
    norme_statique = float(np.linalg.norm(statique) /
        max(np.linalg.norm(khat.data)*np.linalg.norm(base[:, :6]), np.finfo(float).tiny))
    eigenvalues = np.asarray(fem.eigenValues)
    residu_modal = 0.
    if nbModes:
        v = base[:-6, 6:]
        kv = khat[:-6, :-6] @ v
        mv = (mhat[:-6, :-6] @ v)*eigenvalues[None, :]
        residu_modal = float(np.linalg.norm(kv-mv)/(np.linalg.norm(kv)+np.linalg.norm(mv)))
    diagnostics = dict(exudyn_version=exudyn.__version__, dimensions=[n, 6+nbModes], projection=projection,
        modes_internes=nbModes, ports=6, type_base=fem.modeBasis.get("type"),
        erreur_trace_statique=trace_statique, erreur_trace_modale=trace_modes,
        residu_statique_equilibre_relatif=norme_statique,
        residu_modal_equilibre_relatif=residu_modal,
        valeurs_propres=eigenvalues.tolist(),
        ecart_projection_K_vers_D_relatif=float(np.linalg.norm(standard-energie)/np.linalg.norm(energie)),
        ecart_projection_statique_relatif=float(np.linalg.norm(standard[:6, :6]-energie[:6, :6])
                                                /np.linalg.norm(energie[:6, :6])),
        asymetrie_projection_K=float(np.linalg.norm(kr_brut-kr_brut.T)),
        asymetrie_projection_M=float(np.linalg.norm(mr_brut-mr_brut.T)),
        convention="triplets algébriques ; allBoundaryNodes ; aucune interface RBE ni FFRF physique",
        preparation_entree_s=preparation_entree, routine_HCB_officielle_s=apres_hcb-avant_hcb,
        projections_s=termine_projection-apres_hcb,
        spectre_reduit_s=termine_spectre-termine_projection,
        preparation_totale_s=termine_spectre-debut,
        diagnostics_s=time.perf_counter()-termine_spectre,
        cout_total_avec_diagnostics_s=time.perf_counter()-debut)
    return dict(Tphys=tphys, Kred=kr, Mred=mr, K=k, M=m, Q=q,
                valeurs_reduites=valeurs_reduites, projection_force=projection_force,
                diagnostics=diagnostics)


def reponse(preparation, omega, force_physique, *, resolution="spectrale"):
    """Projection/reconstruction de l'enveloppe, pas solveur Exudyn.

    Accepte un vecteur physique de six coefficients ou plusieurs RHS (6,r).
    resolution='directe' sert seulement au contrôle du calcul spectral.
    """
    t = preparation["Tphys"]
    force_physique = np.asarray(force_physique, dtype=float)
    if force_physique.ndim not in (1, 2) or force_physique.shape[0] != 6:
        raise ValueError("force physique de dimension (6,) ou (6,r) requise")
    if resolution == "spectrale":
        denom = preparation["valeurs_reduites"]-float(omega)**2
        charge = preparation["projection_force"] @ force_physique
        coef = charge/denom if charge.ndim == 1 else charge/denom[:, None]
        return t @ (preparation["Q"] @ coef)
    if resolution != "directe":
        raise ValueError("résolution spectrale ou directe requise")
    f = np.zeros((t.shape[0],)+force_physique.shape[1:])
    f[-6:] = force_physique
    a = preparation["Kred"]-omega**2*preparation["Mred"]
    rhs = t.T @ f
    # Équilibrage du petit système pour que les différentes échelles des
    # colonnes de base HCB ne dégradent pas artificiellement sa résolution.
    eq = 1/np.sqrt(np.abs(np.diag(a)))
    if rhs.ndim == 1:
        q = eq*solve(eq[:, None]*a*eq[None, :], eq*rhs, assume_a="sym")
    else:
        q = eq[:, None]*solve(eq[:, None]*a*eq[None, :], eq[:, None]*rhs, assume_a="sym")
    return t @ q


def transfert_ports(preparation, omega):
    """Six RHS physiques, uniquement les six déplacements de port."""
    b = preparation["projection_force"]
    return b.T @ (b/(preparation["valeurs_reduites"]-float(omega)**2)[:, None])


def _lire_sparse(fichier, nom):
    return csr_matrix((fichier[nom+"_data"], fichier[nom+"_indices"],
                       fichier[nom+"_indptr"]), shape=tuple(fichier[nom+"_shape"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture", type=Path, help="NPZ sans pickle : d/m CSR, force, omega, references")
    parser.add_argument("--modes", type=int, nargs="+", default=[0, 6, 12])
    parser.add_argument("--projection", choices=["standard", "energie"], default="standard")
    parser.add_argument("--controle-geometrie", action="store_true")
    args = parser.parse_args()
    with np.load(args.fixture, allow_pickle=False) as f:
        d, m = _lire_sparse(f, "d"), _lire_sparse(f, "m")
        force, omegas, references = f["force"], f["omega"], f["references"]
    sorties = []
    for modes in args.modes:
        np.random.seed(3081)
        p = prepare(d, m, modes, projection=args.projection)
        reponses = [reponse(p, w, force) for w in omegas]
        erreurs = []
        for w, x, ref in zip(omegas, reponses, references):
            e = x-ref
            erreurs.append(dict(omega=float(w), ports=x[-6:].tolist(),
                erreur_relative_masse=float(np.sqrt((e @ (m @ e))/(ref @ (m @ ref)))),
                erreur_relative_deformation=float(np.linalg.norm(d @ e)/np.linalg.norm(d @ ref))))
        sortie = dict(fixture=str(args.fixture), diagnostic=p["diagnostics"], reponses=erreurs)
        if args.controle_geometrie:
            positions = np.random.default_rng(23).normal(size=(d.shape[1]//3, 3))
            np.random.seed(3081)
            autre = prepare(d, m, modes, projection=args.projection, positions=positions)
            sortie["ecart_reponses_positions_auxiliaires"] = max(
                float(np.linalg.norm(reponse(autre, w, force)-x)/np.linalg.norm(x))
                for w, x in zip(omegas, reponses))
        sorties.append(sortie)
    print(json.dumps(sorties, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
