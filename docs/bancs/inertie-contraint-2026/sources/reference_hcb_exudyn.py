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


def _preparer_reponse_reduite(kr, mr, tphys):
    """Prépare la réponse sans dépendre de la génération de la base HCB.

    La paire équilibrée cible est conservée avant diagonalisation. Les
    valeurs/vecteurs arrondis fournissent un préconditionneur approché et
    ne redéfinissent pas cette paire. Tous ces calculs sont comptabilisés
    dans la préparation de ``prepare``.
    """
    eq = 1/np.sqrt(np.diag(kr))
    ke = eq[:, None]*kr*eq[None, :]
    me = eq[:, None]*mr*eq[None, :]
    valeurs, v = eigh(ke, me)
    q = eq[:, None]*v
    te = tphys*eq[None, :]
    return dict(Q=q, valeurs_reduites=valeurs,
                projection_force=q.T@tphys[-6:].T,
                echelles_reduites=eq, Kred_equilibree=ke, Mred_equilibree=me,
                vecteurs_reduits_equilibres=v, Tphys_equilibree=te,
                projection_force_equilibree=te[-6:].T)


def prepare(D, M, nbModes, *, projection="standard", positions=None,
            diagnostics_complets=True):
    """Retourne base physique, projections et diagnostics HCB officiels.

    Les six dernières coordonnées sont les ports physiques. Construction de
    K=D.T@D, équilibrage symétrique S=diag(K)^(-1/2), appel officiel HCB,
    puis Tphys=S@basis. Projection standard : Tphys.T K Tphys. Variante
    énergie : (D Tphys).T (D Tphys), avec exactement la même routine HCB.
    La masse est toujours projetée par Tphys.T M Tphys.
    La diagonalisation équilibrée de la petite paire et les données des
    deux corrections résiduelles sont incluses dans la préparation.
    ``diagnostics_complets=False`` conserve les contrôles de dimensions et
    trace, sans calculer la projection alternative ni les résidus d'audit.
    Le coût des diagnostics optionnels est identifié séparément.
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
    donnees_reponse = _preparer_reponse_reduite(kr, mr, tphys)
    termine_spectre = time.perf_counter()
    # L'audit complet calcule notamment l'autre projection O(n*r²).
    # Il reste disponible pour les sondes, sans être imposé à la mesure.
    norme_statique = residu_modal = ecart_projection = ecart_statique = None
    asymetrie_k = asymetrie_m = None
    eigenvalues = np.asarray(fem.eigenValues)
    if diagnostics_complets:
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
        residu_modal = 0.
        if nbModes:
            v = base[:-6, 6:]
            kv = khat[:-6, :-6] @ v
            mv = (mhat[:-6, :-6] @ v)*eigenvalues[None, :]
            residu_modal = float(np.linalg.norm(kv-mv)/(np.linalg.norm(kv)+np.linalg.norm(mv)))
        ecart_projection = float(np.linalg.norm(standard-energie)/np.linalg.norm(energie))
        ecart_statique = float(np.linalg.norm(standard[:6, :6]-energie[:6, :6])
                               /np.linalg.norm(energie[:6, :6]))
        asymetrie_k = float(np.linalg.norm(kr_brut-kr_brut.T))
        asymetrie_m = float(np.linalg.norm(mr_brut-mr_brut.T))
    diagnostics = dict(complet=bool(diagnostics_complets),
        exudyn_version=exudyn.__version__, dimensions=[n, 6+nbModes], projection=projection,
        modes_internes=nbModes, ports=6, type_base=fem.modeBasis.get("type"),
        erreur_trace_statique=trace_statique, erreur_trace_modale=trace_modes,
        residu_statique_equilibre_relatif=norme_statique,
        residu_modal_equilibre_relatif=residu_modal,
        valeurs_propres=eigenvalues.tolist(),
        ecart_projection_K_vers_D_relatif=ecart_projection,
        ecart_projection_statique_relatif=ecart_statique,
        asymetrie_projection_K=asymetrie_k,
        asymetrie_projection_M=asymetrie_m,
        convention="triplets algébriques ; allBoundaryNodes ; aucune interface RBE ni FFRF physique",
        resolutions_disponibles=["spectrale", "spectrale_corrigee", "directe"],
        reponse_spectrale_corrigee=dict(
            nombre_corrections=2,
            residu="b - (Kred_equilibree - omega² Mred_equilibree) q",
            preconditionneur="V diag(valeurs_reduites - omega²)^(-1) Vᵀ",
            matrices_et_projections_preparees=True,
            refus_pole="dénominateur modal exactement nul ; aucune marge certifiée",
            garantie="aucune contraction ou précision machine garantie a priori"),
        preparation_entree_s=preparation_entree, routine_HCB_officielle_s=apres_hcb-avant_hcb,
        projections_s=termine_projection-apres_hcb,
        spectre_reduit_s=termine_spectre-termine_projection,
        preparation_totale_s=termine_spectre-debut,
        diagnostics_s=time.perf_counter()-termine_spectre,
        cout_total_avec_diagnostics_s=time.perf_counter()-debut)
    return dict(Tphys=tphys, Kred=kr, Mred=mr, K=k, M=m,
                diagnostics=diagnostics, **donnees_reponse)


def reponse(preparation, omega, force_physique, *, resolution="spectrale"):
    """Projection/reconstruction de l'enveloppe, pas solveur Exudyn.

    Accepte un vecteur physique de six coefficients ou plusieurs RHS (6,r).
    ``spectrale_corrigee`` effectue deux corrections dans la petite paire
    équilibrée avec l'inverse modal comme préconditionneur. Cela ne corrige
    pas l'erreur de réduction ni celle de formation/projection de K.
    ``spectrale`` et ``directe`` restent disponibles pour l'audit.
    """
    t = preparation["Tphys"]
    if (np.iscomplexobj(force_physique) or np.iscomplexobj(omega)
            or np.ndim(omega) != 0):
        raise ValueError("pulsation scalaire réelle et forces réelles requises")
    force_physique = np.asarray(force_physique, dtype=float)
    omega = float(omega)
    if (force_physique.ndim not in (1, 2) or force_physique.shape[0] != 6
            or (force_physique.ndim == 2 and force_physique.shape[1] == 0)):
        raise ValueError("force physique de dimension (6,) ou (6,r) requise")
    if not np.isfinite(omega) or omega < 0 or not np.all(np.isfinite(force_physique)):
        raise ValueError("pulsation positive ou nulle et forces finies requises")
    if resolution in ("spectrale", "spectrale_corrigee"):
        denom = preparation["valeurs_reduites"]-omega**2
        if np.any(denom == 0):
            raise ArithmeticError("préconditionneur spectral singulier : dénominateur modal nul")
        charge = preparation["projection_force"] @ force_physique
        coef = charge/denom if charge.ndim == 1 else charge/denom[:, None]
        if resolution == "spectrale":
            return t @ (preparation["Q"] @ coef)
        v = preparation["vecteurs_reduits_equilibres"]
        qb = v@coef
        rhs = preparation["projection_force_equilibree"]@force_physique
        a = preparation["Kred_equilibree"]-omega**2*preparation["Mred_equilibree"]
        for _ in range(2):
            residu = rhs-a@qb
            modal = v.T@residu
            correction = modal/denom if modal.ndim == 1 else modal/denom[:, None]
            qb += v@correction
        return preparation["Tphys_equilibree"]@qb
    if resolution != "directe":
        raise ValueError("résolution spectrale, spectrale_corrigee ou directe requise")
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
