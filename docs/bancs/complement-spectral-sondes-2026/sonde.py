"""Sonde flottante du complément spectral, sans prétention de certification.

Le modèle énergie reste le D binary64 fourni. Les directions approchées
proviennent du facteur QR et sont auditées dans D. Aucun solveur concurrent
n'intervient. Les temps isolés sont indicatifs, sans protocole comparatif.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import time

THREADS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS")
for key in THREADS:
    if os.environ.get(key) != "1":
        raise RuntimeError(key + " doit valoir 1 avant le lancement")

import numpy as np
import scipy
from scipy.linalg import cholesky, null_space, solve, solve_triangular, svdvals
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import LinearOperator, eigsh
import vinkulum
from vinkulum._ports.condensation_energie import CondensationEnergie
from vinkulum._ports.inverse_selectionnee import InverseSelectionnee
from vinkulum._ports.certificat_spectral import certifier_spectral


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        while data := f.read(1024 * 1024):
            h.update(data)
    return h.hexdigest()


def lire(path):
    with np.load(path, allow_pickle=False) as z:
        matrices = [csr_matrix((z[n+"_data"], z[n+"_indices"], z[n+"_indptr"]),
                               shape=tuple(z[n+"_shape"])) for n in ("d", "m")]
        return (*matrices, z["metrique"].copy())


def trace_complement(qr, m, trace, b, eta):
    debut = time.perf_counter()
    u = qr.solve(b)
    h = (b.T @ u + u.T @ b) / 2
    g = u.T @ (m @ u)
    g = (g + g.T) / 2
    sigma = svdvals(h)
    # L'équilibrage ne change que la résolution du petit système. B et les
    # contraintes ciblées restent exactement les valeurs stockées à l'entrée.
    e = 1 / np.sqrt(np.diag(h))
    he = e[:, None] * h * e[None, :]
    ge = e[:, None] * g * e[None, :]
    chi = float(np.trace(solve(he, ge, assume_a="pos")))
    tau = float(trace - chi)
    borne = (1 - eta) / tau if tau > 0 else None
    bn = b / np.linalg.norm(b, axis=0)
    bs = svdvals(bn)
    residu_q = qr.produit(u) - b
    residu_d = qr.di.T @ (qr.di @ u) - b
    dual_q = qr.dual(residu_q)
    dual_d = qr.dual(residu_d)
    norme_b = np.linalg.norm(qr.dual(b), axis=0)
    return dict(
        retenues=b.shape[1], trace_totale=float(trace), trace_retiree=chi,
        trace_complement=tau, facteur_annulation=(abs(trace)+abs(chi))/abs(tau),
        lambda_complement_estimee=borne,
        frequence_limite_estimee_hz=np.sqrt(borne)/(2*np.pi) if borne else None,
        bande_40_hz_flottante=bool(borne and borne > (2*np.pi*40)**2),
        condition_H=float(sigma[0]/sigma[-1]),
        condition_H_equilibre=float(np.linalg.cond(he)),
        sigma_min_B_colonnes_normalisees=float(bs[-1]),
        residu_resolution_QR_dual_relatif=float(np.max(np.linalg.norm(dual_q, axis=0)/norme_b)),
        residu_resolution_D_dual_relatif=float(np.max(np.linalg.norm(dual_d, axis=0)/norme_b)),
        calcul_s=time.perf_counter()-debut,
        certification_machine=False,
    )


def cas(n, entrees, sortie):
    path = entrees / f"n{n}-f40.npz"
    d, m, metrique = lire(path)
    total = m.shape[0]
    interieur, ports = np.arange(total-6), np.arange(total-6, total)
    t = time.perf_counter()
    qr = CondensationEnergie(d, interieur, ports, metrique)
    temps_qr = time.perf_counter()-t
    mi = m[interieur][:, interieur].tocsr()
    masse = mi.diagonal()
    if (mi - diags(masse)).nnz or np.any(masse <= 0):
        raise ValueError("Cette sonde exige la masse diagonale native positive")
    t = time.perf_counter()
    inv = InverseSelectionnee(qr.r, mi, qr.echelles)
    temps_trace = time.perf_counter()-t
    cert = certifier_spectral(qr.di, inv)
    eta = cert["eta_superieur"]
    racine = np.sqrt(masse)
    appels = 0
    def produit_inverse(x):
        nonlocal appels
        appels += 1
        vectoriel = x.ndim == 1
        xx = x[:, None] if vectoriel else x
        y = racine[:, None]*qr.solve(racine[:, None]*xx)
        return y[:, 0] if vectoriel else y
    op = LinearOperator((len(interieur), len(interieur)), matvec=produit_inverse,
                        matmat=produit_inverse, dtype=float)
    t = time.perf_counter()
    rng = np.random.default_rng(107+n)
    valeurs_inverse, z = eigsh(op, k=12, which="LM", tol=1e-11,
                              v0=rng.standard_normal(len(interieur)), ncv=40)
    ordre = np.argsort(valeurs_inverse)[::-1]
    valeurs_inverse, z = valeurs_inverse[ordre], z[:, ordre]
    phi = z/racine[:, None]
    temps_modes = time.perf_counter()-t
    dp = qr.di @ phi
    kproj = (dp.T@dp)
    mproj = phi.T@(mi@phi)
    valeurs_d = np.diag(kproj)/np.diag(mproj)
    residu = qr.di.T@dp-(mi@phi)*valeurs_d[None, :]
    residu_dual = qr.dual(residu)
    residus = np.linalg.norm(residu_dual, axis=0)/np.sqrt(np.diag(kproj))
    references = dict(
        n_poutres=n, dimension_interieure=len(interieur), entree_sha256=sha(path),
        reference_sha256=sha(path.with_suffix(".ref.npy")),
        certificat_spectral_initial=cert,
        trace_QR_float=inv.trace_masse,
        frequence_limite_initiale_certifiee_hz=np.sqrt(cert["lambda_min"])/(2*np.pi),
        frequences_modes_D_rayleigh_hz=(np.sqrt(valeurs_d)/(2*np.pi)).tolist(),
        frequences_modes_QR_inverse_hz=(1/np.sqrt(valeurs_inverse)/(2*np.pi)).tolist(),
        orthogonalite_massique=float(np.linalg.norm(mproj-np.eye(12), 2)),
        residus_modes_D_dual_energie_relatifs=residus.tolist(),
        appels_inverse_eigsh=appels,
        couts_indicatifs_s=dict(qr=temps_qr, trace_selectionnee=temps_trace,
            certificat_spectral_initial=cert["preparation_s"], modes_12=temps_modes),
        certification_complement_machine=False,
    )
    b_all = mi @ phi  # Cible de contraintes : ces valeurs binary64 exactes.
    np.savez_compressed(sortie/f"n{n}-directions.npz", phi=phi, b=b_all,
        r_data=qr.r.data, r_indices=qr.r.indices, r_indptr=qr.r.indptr,
        r_shape=qr.r.shape, echelles=qr.echelles)
    references["directions_sha256"] = sha(sortie/f"n{n}-directions.npz")
    references["selection"] = [trace_complement(qr, mi, inv.trace_masse, b_all[:, :k], eta)
                              for k in (1, 2, 4, 6, 12)]
    if n == 32:
        # Contre-calcul explicite de petite taille seulement : une base
        # massiquement orthonormale du complément, puis une trace positive.
        # Toutes ces opérations restent flottantes et non certifiées.
        references["controle_complement_explicite"] = []
        for k in (1, 2, 4, 6, 12):
            zz = null_space((b_all[:, :k]/racine[:, None]).T)
            cc = qr.r@(zz/(racine*qr.echelles)[:, None])
            kc = cc.T@cc
            l = cholesky(kc, lower=True)
            il = solve_triangular(l, np.eye(l.shape[0]), lower=True)
            tau_positive = float(np.sum(il*il))
            references["controle_complement_explicite"].append(dict(
                retenues=k, trace_positive_QR=tau_positive,
                lambda_min_QR_flottante=float(np.linalg.eigvalsh(kc)[0])))
    references["sensibilite"] = []
    # Mélange déterministe vers l'orthogonal massique de toutes les 12
    # directions. Aucune direction perturbée n'est présumée mode propre.
    bruit = rng.standard_normal(z.shape)
    bruit -= z @ (z.T @ bruit)
    bruit, _ = np.linalg.qr(bruit, mode="reduced")
    for k in (2, 4):
        for amplitude in (1e-8, 1e-4, 1e-2, .1):
            zp, _ = np.linalg.qr(z[:, :k]+amplitude*bruit[:, :k], mode="reduced")
            pp = zp/racine[:, None]
            r = trace_complement(qr, mi, inv.trace_masse, mi@pp, eta)
            r["perturbation_sous_espace"] = amplitude
            references["sensibilite"].append(r)
    references["rotation_vers_modes_suivants"] = []
    for k in (1, 2, 4):
        for angle in (.01, .1, .5, 1., 1.4):
            zp = np.cos(angle)*z[:, :k]+np.sin(angle)*z[:, k:2*k]
            r = trace_complement(qr, mi, inv.trace_masse, mi@(zp/racine[:, None]), eta)
            r["angle_principal_radians"] = angle
            references["rotation_vers_modes_suivants"].append(r)
    references["reechelonnement"] = []
    for k in (2, 4, 12):
        facteurs = np.geomspace(1e-6, 1e6, k)
        r = trace_complement(qr, mi, inv.trace_masse, b_all[:, :k]*facteurs[None, :], eta)
        r["facteurs_colonnes"] = facteurs.tolist()
        references["reechelonnement"].append(r)
    return references


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--entrees", type=Path, default=Path("/tmp/vinkulum-confrontation-ports-0.10.0-corrigee"))
    p.add_argument("--sortie", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--n", type=int, nargs="+", default=[32, 128, 512])
    args = p.parse_args()
    rapport = dict(
        statut="sonde flottante : aucun certificat du complément ni réponse 40 Hz produite",
        date_utc=datetime.now(timezone.utc).isoformat(), python=sys.version,
        vinkulum=vinkulum.__version__, numpy=np.__version__, scipy=scipy.__version__,
        fils={key: os.environ[key] for key in THREADS},
        cpu=sorted(os.sched_getaffinity(0)), script_sha256=sha(__file__), cas=[])
    for n in args.n:
        print("début", n, flush=True)
        result = cas(n, args.entrees, args.sortie)
        rapport["cas"].append(result)
        (args.sortie/"rapport.json").write_text(json.dumps(rapport, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
        print("terminé", n, "Hz", [x["frequence_limite_estimee_hz"] for x in result["selection"]], flush=True)


if __name__ == "__main__":
    main()
