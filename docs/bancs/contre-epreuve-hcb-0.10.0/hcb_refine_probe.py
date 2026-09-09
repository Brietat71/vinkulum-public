#!/usr/bin/env python3
"""Contre-épreuve de l'enveloppe harmonique HCB, sans modifier la campagne.

Exemple, à lancer seulement lorsque la campagne concurrente est arrêtée :
  python /tmp/hcb_refine_probe.py n32-f20.npz n32-f20.ref.npy \
      --modes 185 --output /tmp/hcb-refinement-n32-f20.json --cpu 8

La base vient de l'API officielle via reference_hcb_exudyn. Les corrections
portent exclusivement sur la petite paire projetée, pas sur son exactitude
vis-à-vis du modèle D.T D. Les temps sont une observation, sans répétitions.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import warnings


def empreinte(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("entree", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--modes", type=int, default=185)
    parser.add_argument("--projection", nargs="+", choices=("standard", "energie"),
                        default=["standard", "energie"])
    parser.add_argument("--ci", type=Path,
                        default=Path("/tmp/vinkulum-energie-native-travail/ci"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cpu", type=int)
    args = parser.parse_args()
    for nom in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "RAYON_NUM_THREADS"):
        os.environ[nom] = "1"
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    sys.path.insert(0, str(args.ci.resolve()))

    import numpy as np
    import scipy
    import exudyn
    from scipy.linalg import solve
    from confronte_ports_exudyn import lire, juger, normes
    from reference_hcb_exudyn import prepare, reponse

    d, m, metrique, forces, omega = lire(args.entree)
    reference = np.load(args.reference, allow_pickle=False)
    if len(omega) != 257 or forces.shape != (6, 6):
        raise ValueError("La contre-épreuve exige les 257 fréquences et six charges de la campagne")
    if reference.shape != (len(omega), d.shape[1], 6):
        raise ValueError("Dimensions de la référence incompatibles")
    rapport = dict(
        schema=1, date_utc=datetime.now(timezone.utc).isoformat(),
        portee="Contre-épreuve de notre enveloppe harmonique, base HCB officielle inchangée",
        temps="Une observation par méthode, imports/lecture/juge exclus ; aucun classement de vitesse",
        python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
        exudyn=exudyn.__version__, cpu=sorted(os.sched_getaffinity(0)),
        modes=args.modes, dimension_physique=d.shape[1],
        frequences_hz=(omega/(2*np.pi)).tolist(), charges=6,
        sources_sha256={
            "sonde": empreinte(__file__),
            "pont_hcb": empreinte(args.ci/"reference_hcb_exudyn.py"),
            "juge": empreinte(args.ci/"confronte_ports_exudyn.py")},
        entree_sha256=empreinte(args.entree), reference_sha256=empreinte(args.reference),
        projections=[])

    for projection in args.projection:
        np.random.seed(42123)
        debut = time.perf_counter()
        preparation = prepare(d, m, args.modes, projection=projection)
        temps_preparation = time.perf_counter()-debut
        kr, mr = preparation["Kred"], preparation["Mred"]
        eq = 1/np.sqrt(np.diag(kr))
        ke = eq[:, None]*kr*eq[None, :]
        me = eq[:, None]*mr*eq[None, :]
        t_equilibre = preparation["Tphys"]*eq[None, :]
        # prepare stocke Q = diag(eq) V, où V diagonalise la paire équilibrée.
        v = preparation["Q"]/eq[:, None]
        valeurs = preparation["valeurs_reduites"]
        rhs = t_equilibre[-6:].T@forces
        rhs_normes = np.linalg.norm(rhs, axis=0)
        if np.any(rhs_normes <= 0):
            raise ArithmeticError("Charge réduite nulle")
        modal_rhs = v.T@rhs

        def spectrale(w):
            return reponse(preparation, w, forces, resolution="spectrale"), None

        def directe(w):
            # Même petite paire ; équilibrage diagonal dépendant de la
            # fréquence dans le pont public déjà utilisé par les sondages.
            return reponse(preparation, w, forces, resolution="directe"), None

        def directe_equilibre_fixe(w):
            qb = solve(ke-w*w*me, rhs, assume_a="sym")
            return t_equilibre@qb, None

        def spectrale_corrigee(w):
            denominateurs = valeurs-w*w
            a = ke-w*w*me
            qb = v@(modal_rhs/denominateurs[:, None])
            residus = []
            for iteration in range(3):
                residu = rhs-a@qb
                residus.append(float(np.max(np.linalg.norm(residu, axis=0)/rhs_normes)))
                if iteration < 2:
                    qb += v@((v.T@residu)/denominateurs[:, None])
            return t_equilibre@qb, residus

        cas = dict(projection=projection, preparation_s=temps_preparation,
                   diagnostic_hcb=preparation["diagnostics"],
                   spectre_reduit_min=float(valeurs.min()),
                   spectre_reduit_max=float(valeurs.max()),
                   methodes=[])
        rapport["projections"].append(cas)
        methodes = (("spectrale", spectrale),
                    ("directe", directe),
                    ("spectrale_deux_corrections", spectrale_corrigee),
                    ("directe_equilibre_fixe", directe_equilibre_fixe))
        for nom, calculer in methodes:
            champs = np.empty_like(reference)
            nm = np.empty((len(omega), 6))
            nd = np.empty_like(nm)
            historiques = []
            resultat = dict(methode=nom)
            try:
                with warnings.catch_warnings(record=True) as avertissements:
                    warnings.simplefilter("always")
                    debut = time.perf_counter()
                    for j, w in enumerate(omega):
                        champs[j], historique = calculer(w)
                        nm[j], nd[j] = normes(d, m, champs[j])
                        if historique is not None:
                            historiques.append(historique)
                    temps = time.perf_counter()-debut
                    resultat.update(
                        reponses_champs_et_normes_s=temps,
                        avertissements=sorted(set(str(a.message) for a in avertissements)))
                juge = juger(d, m, metrique, champs, reference)
                resultat.update(statut="termine", accepte=juge["accepte"],
                    maxima=juge["maxima"], maxima_operateurs=juge["maxima_operateurs"],
                    operateurs=juge["operateurs"],
                    frequences_pires_hz={n: float(omega[np.argmax(x)]/(2*np.pi))
                                          for n, x in juge["operateurs"].items()},
                    champ_sha256=hashlib.sha256(champs.tobytes()).hexdigest())
                if historiques:
                    resultat["residus_reduits_par_frequence_initial_apres_1_apres_2"] = historiques
                    resultat["residus_reduits_max_initial_apres_1_apres_2"] = \
                        np.max(historiques, axis=0).tolist()
            except Exception as exc:
                resultat.update(statut="refus_ou_echec", erreur=type(exc).__name__+": "+str(exc))
            cas["methodes"].append(resultat)
            print(json.dumps(dict(projection=projection, methode=nom,
                    statut=resultat["statut"], accepte=resultat.get("accepte"),
                    maxima_operateurs=resultat.get("maxima_operateurs"),
                    reponses_s=resultat.get("reponses_champs_et_normes_s")),
                    ensure_ascii=False, allow_nan=False), flush=True)

    contenu = json.dumps(rapport, ensure_ascii=False, indent=2, allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(contenu)
        print("Rapport : "+str(args.output), flush=True)
    else:
        print(contenu)


if __name__ == "__main__":
    main()
