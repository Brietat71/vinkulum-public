"""Témoin reproductible de l'API native des ports, sans comparaison concurrente.

Exemple : python ci/experience_reduction_native.py --n 32 128 512
Vérification de l'archive : python ci/experience_reduction_native.py --verifier

L'oracle dynamique reconstruit la flexion de Timoshenko depuis les paramètres
physiques et résout une chaîne de blocs 2x2 en Decimal, sans lire D ni K natif.
Les deux précisions de l'oracle contrôlent sa convergence, pas ses arrondis.
"""
import argparse
import copy
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time


THREADS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS")
for _variable in THREADS:
    os.environ.setdefault(_variable, "1")

import numpy as np
import scipy
import vinkulum
from vinkulum import Noyau
from vinkulum.reduction_ports import reduire_poutres


ROOT = Path(__file__).resolve().parents[1]
SORTIE = ROOT / "docs/bancs/reduction-native-0.10.0.json"
SOURCES = ("reduction_ports.py", "_ports/condensation_energie.py",
           "_ports/ports_krylov.py", "_ports/ports_releves.py",
           "_ports/champ_interieur.py", "_ports/inverse_selectionnee.py",
           "_ports/certificat_spectral.py")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def serialisable(valeur):
    if isinstance(valeur, dict):
        return {str(k): serialisable(v) for k, v in valeur.items()}
    if isinstance(valeur, (list, tuple)):
        return [serialisable(v) for v in valeur]
    if isinstance(valeur, np.ndarray):
        return serialisable(valeur.tolist())
    if isinstance(valeur, np.generic):
        return serialisable(valeur.item())
    if isinstance(valeur, Decimal):
        return str(valeur)
    return valeur


def ecrire(path, document):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(serialisable(document), ensure_ascii=False,
                               indent=2, allow_nan=False) + "\n")


def parametres(n):
    longueur, largeur, hauteur = 1., .04, .006
    young, poisson, rho, cisaillement = 70e9, .3, 2700., 5/6
    section = largeur * hauteur
    iy, iz = largeur * hauteur**3 / 12, hauteur * largeur**3 / 12
    # Approximation usuelle de Saint-Venant pour un rectangle mince,
    # fixée comme donnée de ce modèle ; l'oracle sollicite la flexion seule.
    jt = largeur * hauteur**3 / 3 * (1 - .63*hauteur/largeur + .052*(hauteur/largeur)**5)
    g = young / (2 * (1 + poisson))
    pas = longueur / n
    masse = rho * section * pas
    inerties = [masse*(largeur**2 + hauteur**2)/12,
                masse*(pas**2 + hauteur**2)/12,
                masse*(pas**2 + largeur**2)/12]
    return dict(n=n, longueur_m=longueur, largeur_y_m=largeur, hauteur_z_m=hauteur,
                young_pa=young, poisson=poisson, rho_kg_m3=rho,
                facteur_cisaillement=cisaillement, pas_m=pas,
                ea_n=young*section, gay_n=cisaillement*g*section,
                gaz_n=cisaillement*g*section, gj_nm2=g*jt,
                eiy_nm2=young*iy, eiz_nm2=young*iz,
                masse_noeud_interieur_kg=masse,
                inerties_noeud_interieur_kg_m2=inerties,
                masse_totale_avant_encastrement_kg=rho*section*longueur,
                force_port_physique=[0., 0., .1, 0., 0., 0.],
                masse_discrete=("Poids nodaux 1/2 aux deux extrémités, 1 ailleurs. "
                                 "Chaque inertie est ce même poids de celle d'un "
                                 "parallélépipède pas×largeur×hauteur. La masse de "
                                 "la racine est retirée avec ses six coordonnées. "
                                 "Masse concentrée diagonale, sans matrice de masse consistante."))


def modele(p):
    n = Noyau([0., 0., 0.])
    for j in range(p["n"] + 1):
        poids = .5 if j in (0, p["n"]) else 1.
        inertie = np.diag(np.asarray(p["inerties_noeud_interieur_kg_m2"]) * poids)
        n.corps(str(j), poids*p["masse_noeud_interieur_kg"], inertie.ravel().tolist(),
                [j*p["pas_m"], 0., 0.])
    n.liaison("encastrement", None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    for j in range(p["n"]):
        n.poutre(str(j), j, j+1, p["ea_n"], p["gay_n"], p["gj_nm2"], p["eiy_nm2"],
                 ei3=p["eiz_nm2"], ga3=p["gaz_n"], formulation="integree")
    libres = np.arange(6, 6*(p["n"]+1))
    ports = np.arange(6*p["n"], 6*(p["n"]+1))
    l = p["longueur_m"]
    # Les deux blocs de flexion comprennent la compliance de translation
    # et rotation, sans mélanger leurs unités par une norme arbitraire.
    metrique = np.zeros((6, 6))
    metrique[0, 0], metrique[3, 3] = p["ea_n"]/l, p["gj_nm2"]/l
    for ids, ei, ga, signes in (([1, 5], p["eiz_nm2"], p["gay_n"], [1., 1.]),
                                ([2, 4], p["eiy_nm2"], p["gaz_n"], [1., -1.])):
        compliance = np.array([[l**3/(3*ei)+l/ga, l*l/(2*ei)], [l*l/(2*ei), l/ei]])
        metrique[np.ix_(ids, ids)] = np.linalg.inv(compliance) * np.outer(signes, signes)
    return n, libres, ports, metrique


def _decimal(x):
    return Decimal.from_float(float(x))


def reference_statique(p, precision=90):
    with localcontext() as contexte:
        contexte.prec = precision
        l, h, ei, ga, force = map(_decimal, (p["longueur_m"], p["pas_m"], p["eiy_nm2"],
                                            p["gaz_n"], p["force_port_physique"][2]))
        champ = np.zeros((p["n"], 6))
        for j in range(1, p["n"]+1):
            x = j*h
            champ[j-1, 2] = float(force*(x*x*(3*l-x)/(6*ei)+x/ga))
            champ[j-1, 4] = -float(force*x*(2*l-x)/(2*ei))
        return champ.ravel()


def _solve2(a, b):
    det = a[0][0]*a[1][1] - a[0][1]*a[1][0]
    if det == 0:
        raise ArithmeticError("pivot nul de l'oracle dynamique Decimal")
    return [[(a[1][1]*b[0][j]-a[0][1]*b[1][j])/det for j in range(len(b[0]))],
            [(-a[1][0]*b[0][j]+a[0][0]*b[1][j])/det for j in range(len(b[0]))]]


def reference_dynamique(p, omega, precision):
    """Flexion (uz, -ry), matrice nodale exacte de Timoshenko et masses données."""
    with localcontext() as contexte:
        contexte.prec = precision
        zero = Decimal(0)
        h, ga, ei, masse, inertie = map(_decimal, (p["pas_m"], p["gaz_n"], p["eiy_nm2"],
            p["masse_noeud_interieur_kg"], p["inerties_noeud_interieur_kg_m2"][1]))
        z = _decimal(omega)**2
        ceff = 1/(1/ga+h*h/(12*ei))
        a, b = ceff/h, ceff/2
        c, dd = ei/h+ceff*h/4, -ei/h+ceff*h/4
        liaison = [[-a, b], [-b, dd]]
        transferts = []
        for j in range(p["n"]):
            dernier = j == p["n"]-1
            poids = Decimal("0.5") if dernier else Decimal(1)
            pivot = [[(a if dernier else 2*a)-z*poids*masse, -b if dernier else zero],
                     [-b if dernier else zero, (c if dernier else 2*c)-z*poids*inertie]]
            if j:
                g = transferts[-1]
                pivot = [[pivot[r][s]-sum((liaison[k][r]*g[k][s] for k in range(2)), zero)
                          for s in range(2)] for r in range(2)]
            if dernier:
                x = _solve2(pivot, [[_decimal(p["force_port_physique"][2])], [zero]])
            else:
                transferts.append(_solve2(pivot, liaison))
        solution = [[x[0][0], x[1][0]]]
        for g in reversed(transferts):
            suivant = solution[-1]
            solution.append([-sum((g[r][k]*suivant[k] for k in range(2)), zero) for r in range(2)])
        champ = np.zeros((p["n"], 6))
        for j, v in enumerate(reversed(solution)):
            champ[j, 2], champ[j, 4] = float(v[0]), -float(v[1])
        return champ.ravel()


def bilan_reponse(r, rep, reference, omega):
    champ = np.asarray(rep["champ_physique"])
    defaut = champ-reference
    m, d = r.m, r.qr.d
    erreur_m = float(np.sqrt(max(float(defaut @ (m @ defaut)), 0.)))
    erreur_d = float(np.linalg.norm(d @ defaut))
    norme_m = float(np.sqrt(reference @ (m @ reference)))
    norme_d = float(np.linalg.norm(d @ reference))
    ports_reference = reference[r.qr.s]
    erreur_port = float(np.linalg.norm(np.linalg.solve(rep["normalisation_ports"],
                                  rep["deplacement_ports_physique"]-ports_reference)))
    return dict(omega_rad_s=float(omega), frequence_hz=float(omega/(2*math.pi)),
                deplacement_ports_physique=rep["deplacement_ports_physique"].copy(),
                reference_ports_physique=ports_reference.copy(),
                coordonnees_ports_normalisees=rep["coordonnees_ports_normalisees"].copy(),
                borne_ports_normalises=rep["borne_ports_normalises"],
                erreur_ports_normalises_vers_reference=erreur_port,
                erreur_relative_champ_euclidienne=float(np.linalg.norm(defaut)/np.linalg.norm(reference)),
                masse=dict(rep["masse"], erreur_vers_reference=erreur_m,
                           erreur_relative_vers_reference=erreur_m/norme_m),
                deformation=dict(rep["deformation"], erreur_vers_reference=erreur_d,
                                 erreur_relative_vers_reference=erreur_d/norme_d),
                certification_machine=rep["certification_machine"], marge=rep["marge"],
                borne_schur=rep["borne_schur"],
                champ_float64_sha256=hashlib.sha256(champ.tobytes()).hexdigest())


def cas(n, args):
    p = parametres(n)
    sorties = dict(parametres=p, statut="en_cours", couts_observes_s={})
    debut = time.perf_counter()
    phase = "construction_noyau"
    try:
        noyau, libres, ports, metric = modele(p)
        sorties["couts_observes_s"][phase] = time.perf_counter()-debut
        sorties["metrique_ports"] = metric.copy()
        omegas = 2*math.pi*np.asarray(args.frequences_hz)
        phase, demarre = "reduction_native_et_certificat", time.perf_counter()
        r = reduire_poutres(noyau, libres, ports, metric, float(max(omegas)),
                            tolerance=args.tolerance_schur, max_blocs=1,
                            max_directions=args.max_directions,
                            precision_spectrale=args.precision_spectrale,
                            budget_spectral=args.budget_spectral, budget_qr=args.budget_qr)
        sorties["couts_observes_s"][phase] = time.perf_counter()-demarre
        sorties.update(dimension_physique=len(libres), dimension_interieure=len(r.qr.i), ports=r.p,
                       directions_initiales=r.taille_interieure, certificat_spectral=copy.deepcopy(r.certificat_spectral),
                       audit_facteur=copy.deepcopy(r.audit_modele()),
                       domaine_natif=r.origine["donnees_natives"]["infos_domaine"],
                       residu_contraintes_initial=r.origine["residu_contraintes_initial"],
                       tolerance_contraintes=r.origine["tolerance_contraintes"],
                       lambda_min=r.lambda_min, nnz_D=r.qr.d.nnz, nnz_M=r.m.nnz)
        phase, demarre = "reponses_initiales", time.perf_counter()
        initiales = copy.deepcopy([r.reponse(w, p["force_port_physique"]) for w in omegas])
        sorties["couts_observes_s"][phase] = time.perf_counter()-demarre
        facteur, certificat = r.qr, r.certificat_spectral
        phase, demarre = "controle_sans_budget_enrichissement", time.perf_counter()
        dimension_initiale = r.taille_interieure
        sans_budget = r.adapter(omegas, p["force_port_physique"],
                                tolerance_relative=args.tolerance_relative, max_etapes=0)
        sorties["controle_sans_enrichissement"] = dict(
            max_etapes=0, statut=sans_budget["statut"],
            accepte=bool(sans_budget["historique"][-1]["accepte"]),
            directions_avant=dimension_initiale, directions_apres=r.taille_interieure,
            historique=copy.deepcopy(sans_budget["historique"]))
        sorties["couts_observes_s"][phase] = time.perf_counter()-demarre
        phase, demarre = "adaptation_et_reponses_finales", time.perf_counter()
        adapte = r.adapter(omegas, p["force_port_physique"],
                            tolerance_relative=args.tolerance_relative, max_etapes=args.max_etapes)
        finales = copy.deepcopy(adapte["reponses"])
        sorties["couts_observes_s"][phase] = time.perf_counter()-demarre
        sorties.update(statut=adapte["statut"], historique=copy.deepcopy(adapte["historique"]),
                       directions_finales=r.taille_interieure, facteur_conserve=r.qr is facteur,
                       certificat_conserve=(r.certificat_spectral is certificat
                                            and certificat == sorties["certificat_spectral"]),
                       certification_champ_machine=adapte["certification_machine"])
        phase, demarre = "oracles_independants", time.perf_counter()
        references = [reference_dynamique(p, float(w), 90) for w in omegas]
        ecarts_precision = [float(np.linalg.norm(reference_dynamique(p, float(w), 70)-ref)
                                  /np.linalg.norm(ref)) for w, ref in zip(omegas, references)]
        statique = reference_statique(p)
        dynamique_zero = reference_dynamique(p, 0., 90)
        sorties["references"] = dict(precision_decimal=[70, 90],
            ecarts_relatifs_champ_70_90=ecarts_precision,
            statique_analytique_ports=statique[-6:].copy(),
            ecart_relatif_statique_analytique_dynamique_zero=float(np.linalg.norm(statique-dynamique_zero)
                                                                  /np.linalg.norm(statique)))
        sorties["initiales"] = [bilan_reponse(r, rep, ref, w)
                                 for rep, ref, w in zip(initiales, references, omegas)]
        sorties["finales"] = [bilan_reponse(r, rep, ref, w)
                               for rep, ref, w in zip(finales, references, omegas)]
        sorties["couts_observes_s"][phase] = time.perf_counter()-demarre
        sorties["precision_reference_compatible_tolerance"] = all(
            rep[norme]["erreur_relative_vers_reference"] <= args.tolerance_relative
            for rep in sorties["finales"] for norme in ("masse", "deformation"))
        sorties["acceptation_api"] = bool(adapte["historique"][-1]["accepte"])
    except Exception as erreur:
        sorties.update(statut="refus_ou_echec", phase_refus=phase,
                       erreur=dict(type=type(erreur).__name__, message=str(erreur),
                                   bilan=getattr(erreur, "bilan", None)), acceptation_api=False)
    sorties["couts_observes_s"]["total_temoignage_avec_oracles"] = time.perf_counter()-debut
    return sorties


def environnement():
    paquet = Path(vinkulum.__file__).parent
    extension = next(paquet.glob("*.so"))
    def git(*options):
        result = subprocess.run(["git", *options], cwd=ROOT, capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None
    return dict(date_utc=datetime.now(timezone.utc).isoformat(), python=sys.version,
                executable=sys.executable, vinkulum=vinkulum.__version__, numpy=np.__version__,
                scipy=scipy.__version__, plateforme=platform.platform(),
                cpu_affinite=sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
                fils={k: os.environ.get(k) for k in THREADS},
                script_sha256=sha(__file__), extension_sha256=sha(extension),
                sources_python_sha256={nom: sha(paquet/nom) for nom in SOURCES},
                git_head=git("rev-parse", "HEAD"), git_statut=git("status", "--short"))


def verifier(path):
    archive = json.loads(path.read_text())
    actuel = environnement()
    assert archive["schema"] == 1
    for cle in ("script_sha256", "extension_sha256", "sources_python_sha256"):
        assert archive["environnement"][cle] == actuel[cle], "empreinte différente : " + cle
    protocole = archive["protocole"]
    assert [c["parametres"]["n"] for c in archive["cas"]] == protocole["tailles_demandees"]
    frequences = protocole["frequences_hz"]
    tolerance = protocole["tolerance_relative"]
    controlees = 0
    for resultat in archive["cas"]:
        if resultat["statut"] == "refus_ou_echec":
            assert not resultat["acceptation_api"] and resultat["erreur"]["message"]
            continue
        p = resultat["parametres"]
        np.testing.assert_array_equal(reference_statique(p)[-6:],
                                      resultat["references"]["statique_analytique_ports"])
        assert resultat["facteur_conserve"] and resultat["certificat_conserve"]
        assert resultat["certificat_spectral"]["certification_machine"]
        assert not resultat["certification_champ_machine"]
        assert len(resultat["references"]["ecarts_relatifs_champ_70_90"]) == len(frequences)
        assert len(resultat["initiales"]) == len(resultat["finales"]) == len(frequences)
        for hz, initial, rep in zip(frequences, resultat["initiales"], resultat["finales"]):
            assert initial["omega_rad_s"] == rep["omega_rad_s"] == 2*math.pi*hz
            ref = reference_dynamique(p, rep["omega_rad_s"], 90)
            np.testing.assert_array_equal(ref[-6:], rep["reference_ports_physique"])
            np.testing.assert_array_equal(ref[-6:], initial["reference_ports_physique"])
            assert resultat["lambda_min"] > rep["omega_rad_s"]**2
            assert not rep["certification_machine"]
        assert resultat["acceptation_api"] == resultat["historique"][-1]["accepte"]
        precision_reference = all(math.isfinite(rep[norme]["erreur_relative_vers_reference"])
            and 0 <= rep[norme]["erreur_relative_vers_reference"] <= tolerance
            for rep in resultat["finales"] for norme in ("masse", "deformation"))
        assert precision_reference == resultat["precision_reference_compatible_tolerance"]
        if resultat["acceptation_api"]:
            assert resultat["statut"] == "tolerance_aux_frequences_demandees"
            assert precision_reference, "acceptation API contredite par la précision de l'oracle physique"
            for rep in resultat["finales"]:
                for norme in ("masse", "deformation"):
                    borne = rep[norme]["borne_relative"]
                    assert borne is not None and math.isfinite(borne) and 0 <= borne <= tolerance
        sans_budget = resultat["controle_sans_enrichissement"]
        assert sans_budget["directions_avant"] == sans_budget["directions_apres"]
        assert sans_budget["max_etapes"] == 0
        assert sans_budget["accepte"] == sans_budget["historique"][-1]["accepte"]
        if not sans_budget["accepte"]:
            assert sans_budget["statut"] == "budget_etapes"
        controlees += 1
    print(json.dumps(dict(archive=str(path), cas_verifies=controlees,
                          empreintes_et_oracles="conformes"), ensure_ascii=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, nargs="+", default=[32, 128, 512])
    parser.add_argument("--output", type=Path, default=SORTIE)
    parser.add_argument("--frequences-hz", type=float, nargs="+", default=[0., 1., 2.])
    parser.add_argument("--tolerance-relative", type=float, default=1e-6)
    parser.add_argument("--tolerance-schur", type=float, default=1e-10)
    parser.add_argument("--max-directions", type=int, default=64)
    parser.add_argument("--max-etapes", type=int, default=6)
    parser.add_argument("--budget-qr", type=int, default=10_000_000)
    parser.add_argument("--budget-spectral", type=int, default=10_000_000)
    parser.add_argument("--precision-spectrale", type=int, default=80)
    parser.add_argument("--charge-concurrente", action="store_true")
    parser.add_argument("--verifier", action="store_true")
    args = parser.parse_args()
    if args.verifier:
        verifier(args.output)
        return
    if any(n < 2 for n in args.n) or len(set(args.n)) != len(args.n):
        parser.error("n distincts et supérieurs ou égaux à 2 requis")
    if not all(math.isfinite(w) and w >= 0 for w in args.frequences_hz):
        parser.error("fréquences réelles finies et non négatives requises")
    document = dict(schema=1, objet="Témoin natif de réduction matérielle par ports, version 0.10.0",
        environnement=environnement(),
        protocole=dict(commandes=sys.argv, tailles_demandees=args.n, frequences_hz=args.frequences_hz,
            repetitions=1, echauffement=False, charge_concurrente=args.charge_concurrente,
            note_couts="Durées murales observées une seule fois ; imports exclus, oracles séparés. "
                       "Pas de mesure isolée ni de comparaison de vitesse avec un autre solveur.",
            tolerance_relative=args.tolerance_relative, tolerance_schur=args.tolerance_schur,
            bloc_initial=1, max_directions=args.max_directions, max_etapes=args.max_etapes,
            budget_qr=args.budget_qr, budget_spectral=args.budget_spectral,
            precision_spectrale=args.precision_spectrale),
        limites=["Petites perturbations matérielles au repos : K=DᵀD, sans amortissement ni charges intérieures.",
                 "Masse nodale concentrée explicitement définie ; pas de revendication de convergence vers un continuum dynamique.",
                 "Certificat spectral pour D et M binary64 d'entrée ; contrôles de champ évalués en doubles non certifiés.",
                 "L'oracle reconstruit les coefficients physiques avant les racines du facteur natif : les modèles peuvent différer par arrondis.",
                 "Decimal 70/90 contrôle la convergence observée de l'oracle, sans certificat machine de cet oracle.",
                 "Acceptation adaptative aux seules fréquences demandées, sans garantie d'acceptation entre ces points.",
                 "Un seul passage par taille, sans classement de performances ni campagne exhaustive."], cas=[])
    for n in args.n:
        print(f"console {n} : construction, réduction native, adaptation, oracles", flush=True)
        resultat = cas(n, args)
        document["cas"].append(resultat)
        ecrire(args.output, document)
        print(json.dumps(serialisable({k: resultat[k] for k in
            ("statut", "acceptation_api", "couts_observes_s")}), ensure_ascii=False), flush=True)
    verifier(args.output)


if __name__ == "__main__":
    main()
