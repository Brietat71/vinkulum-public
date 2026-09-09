"""36 processus frais : inertie compilée, Decimal témoin, LU corrigée.

Le contrôleur par facteurs reste identique. Toute la permutation physique
est reconstruite et chronométrée ; aucun certificat ni direction offert.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time
import warnings
import experience_controle_facteurs as precedent
from experience_controle_facteurs import (CIBLE, MAILLAGES, THREADS, FILS_SUPPLEMENTAIRES,
    GAMMA_INERTIE, OMEGA_MAX, BUDGETS_INERTIE, RANGS, DONNEES_FIGEES,
    maintenant, sha, jsonable, ecrire, fichier, donnees, sauver_champs,
    controle_oracle, qualifier, qualifier_controle, verifier_diagnostic_facteurs,
    verifier_reference)

CI=Path(__file__).resolve().parent
FACTEURS=CANDIDATS=INERTIE=CONTROLES=("binaire_controle","decimal_controle")
VARIANTES=CANDIDATS+("lu_corrigee",)
SCHEMA=1
TYPE_CAMPAGNE="inertie_binaire_36_essais"
EXTENSION_FACTEURS="gram"
SOURCES=precedent.SOURCES+("experience_inertie_binaire.py","inertie_binaire_compile.py",
    "inertie_binaire_native.cpp","matrices_certificat.py","ordre_separateurs.py")


def sources():
    return {nom:sha(CI/nom) for nom in SOURCES}


def environnement(variante):
    return precedent.environnement("facteurs_controle" if variante in CANDIDATS else variante)


def parametres_certificat(variante):
    if variante=="binaire_controle":
        return dict(type="inertie_kkt_binary64_compilee",arithmetique="binary64_nextafter_cpp",
            gamma_binary64=GAMMA_INERTIE,gamma_hex=GAMMA_INERTIE.hex(),frequence_seuil_hz=80.,
            permutation="separateurs_bfs",taille_feuille=8,**BUDGETS_INERTIE)
    return precedent.parametres_certificat("facteurs_controle" if variante=="decimal_controle" else variante)


def parametres_reduction(variante):
    return precedent.parametres_reduction("facteurs_controle" if variante in CANDIDATS else variante)


def parametres_controle(variante):
    return precedent.parametres_controle("facteurs_controle" if variante in CANDIDATS else variante)


def worker_candidat(entree, reference, variante, n, dossier):
    import numpy as np
    from scipy.sparse.linalg import LinearOperator, eigsh
    from condensation_energie import CondensationEnergie
    from krylov_contraint import KrylovContraint
    from controle_complement import ControleComplement
    from controle_facteurs import ControleFacteurs
    from inertie_complement_dirigee import certifier_inertie_complement
    from confronte_ports_exudyn import juger, normes

    from inertie_binaire_compile import BibliothequeInertie, certifier_inertie_compilee
    from ordre_separateurs import ordre_separateurs
    bib = BibliothequeInertie(os.environ['VINKULUM_INERTIE_BINARY64_LIB']) if variante == 'binaire_controle' else None
    d, m, metrique, forces, omega = donnees(entree, reference, n)
    parametres = parametres_certificat(variante)
    result = dict(variante=variante, n=n, max_hz=40., modes_demandes=1,
                  type_certificat=parametres["type"], parametres_certificat=parametres,
                  parametres_reduction=parametres_reduction(variante),
                  parametres_controle=parametres_controle(variante),
                  omega_max_binary64=float(omega[-1]),
                  date_utc=maintenant(), environnement=environnement(variante),
                  sources_sha256=sources(), entree_sha256=sha(entree),
                  reference_sha256=sha(reference), nombre_frequences=257, charges=6,
                  cible=CIBLE, certification_machine_reponses=False,
                  dimension_physique=d.shape[1], nnz_d=d.nnz, nnz_m=m.nnz,
                  memoire_comparable=False,
                  memoire_note="VmHWM candidat avant oracle ; témoins historiques avec ru_maxrss")
    champs = np.full((257, d.shape[1], 6), np.nan)
    nm, nd = np.full((257, 6), np.nan), np.full((257, 6), np.nan)
    retours, valides, echecs = [], [], []
    b = phi = reduction = controle = None
    phases = {}
    phase = "condensation"
    prepare = None
    debut = precedent = time.perf_counter()

    def marquer(nom):
        nonlocal precedent
        fin = time.perf_counter()
        phases[nom] = fin-precedent
        precedent = fin

    with warnings.catch_warnings(record=True) as avertissements:
        warnings.simplefilter("always")
        try:
            qr = CondensationEnergie(d, np.arange(d.shape[1]-6),
                                     np.arange(d.shape[1]-6, d.shape[1]), metrique)
            marquer(phase)
            phase = "selection_direction"
            mi = m[qr.i][:, qr.i].tocsr()
            racine = np.sqrt(mi.diagonal())
            appels = 0

            def inverse_masse(v):
                nonlocal appels
                appels += 1
                return (racine[:, None]*qr.solve(racine[:, None]*np.asarray(v).reshape(-1, 1))).ravel()

            op = LinearOperator(mi.shape, matvec=inverse_masse, dtype=np.float64)
            valeur, vecteur = eigsh(op, k=1, which="LA", tol=1e-11, ncv=20,
                v0=np.random.default_rng(107+n).normal(size=len(qr.i)))
            phi = vecteur/racine[:, None]
            b = mi@phi  # Ce produit arrondi définit B, puis le certificat cible ce B exact.
            if valeur[0] <= 0 or not np.all(np.isfinite(b)):
                raise ArithmeticError("direction modale inverse invalide")
            result["selection"] = dict(k=1, ncv=20, tolerance=1e-11, graine=107+n,
                which="LA", appels_inverse=appels, valeurs_inverse=valeur,
                valeur_propre_qr=1/valeur[0], definition="M_ii^1/2 K_Q^-1 M_ii^1/2")
            marquer(phase)
            phase = "certification_complement"
            # Même seuil et mêmes données exactes ; permutation entièrement comptée.
            gamma = float((2*np.pi*80)**2)
            if gamma != GAMMA_INERTIE:
                raise ArithmeticError("seuil binary64 différent du protocole")
            if variante == "binaire_controle":
                debut_ordre = time.perf_counter()
                permutation = ordre_separateurs(d[:,qr.i],mi,taille_feuille=8)
                result['ordre_separateurs_s'] = time.perf_counter()-debut_ordre
                cert = certifier_inertie_compilee(d[:,qr.i],mi,b,gamma,
                    bibliotheque=bib,permutation=permutation,**BUDGETS_INERTIE)
            else:
                cert = certifier_inertie_complement(d[:,qr.i],mi,b,
                    gamma=gamma,precision=32,**BUDGETS_INERTIE)
            result["certificat_complement"] = cert
            marge_exacte = Fraction(cert["lambda_min"])-Fraction(float(omega[-1]))**2
            result["marge_bande_exacte"] = dict(numerateur=str(marge_exacte.numerator),
                                               denominateur=str(marge_exacte.denominator))
            if cert["certification_machine"] is not True or marge_exacte <= 0:
                raise ArithmeticError("bande non démontrée : lambda doit dépasser omega_max² exactement")
            marquer(phase)
            phase = "construction_krylov"
            reduction = KrylovContraint(qr, m, b, phi, cert["lambda_min"], omega[-1],
                                         blocs=4, max_directions=128)
            marquer(phase)
            phase = "construction_controle"
            if variante in FACTEURS:
                controle = ControleFacteurs(reduction, profondeur=8,
                                           extension=EXTENSION_FACTEURS)
            elif variante in CONTROLES:
                controle = ControleComplement(reduction, profondeur=8)
            marquer(phase)
            prepare = precedent
            phase = "reponses"
            solveur = controle if controle is not None else reduction
            for j, w in enumerate(omega):
                try:
                    rep = solveur.reponses(w, forces)
                    champs[j] = rep["champ"]
                    if controle is None:
                        nm[j], nd[j] = normes(d, m, rep["champ"])
                    else:
                        # Normes déjà calculées sur M et D originaux par le contrôleur.
                        nm[j] = rep["bornes"]["masse"]["normes"]
                        nd[j] = rep["bornes"]["deformation"]["normes"]
                    if not all(np.all(np.isfinite(a)) for a in (champs[j], nm[j], nd[j])):
                        raise ArithmeticError("champ ou norme non fini")
                    valides.append(j)
                    retours.append({k: v for k, v in rep.items() if k not in ("champ", "coordonnees")})
                except Exception as exc:
                    echecs.append(dict(frequence=j, omega=float(w), type=type(exc).__name__, message=str(exc)))
                    retours.append(dict(refus=echecs[-1]))
            marquer(phase)
            result["statut"] = "termine" if not echecs else "refus_ou_echec"
        except Exception as exc:
            marquer(phase)
            result.update(statut="refus_ou_echec", phase=phase,
                          erreur=dict(type=type(exc).__name__, message=str(exc),
                                      diagnostic=getattr(exc, "diagnostic", getattr(exc, "bilan", None))))
    termine = precedent
    result.update(phases_s=phases, temps_calcul_s=termine-debut,
                  avertissements=sorted(set(str(w.message) for w in avertissements)),
                  frequences_calculees=len(valides), refus_calcul=echecs)
    if prepare is not None:
        result.update(preparation_s=prepare-debut, reponses_s=termine-prepare,
                      total_s=termine-debut)
    else:
        result["temps_jusquau_refus_s"] = termine-debut
    result["memoire_avant_oracle_kib"] = {
        ligne.split(":")[0]: int(ligne.split()[1])
        for ligne in Path("/proc/self/status").read_text().splitlines()
        if ligne.startswith(("VmHWM:", "VmRSS:"))}
    # Toutes les sauvegardes et diagnostics supplémentaires commencent ici.
    hors = time.perf_counter()
    result["champs"] = sauver_champs(dossier, champs, valides)
    if b is not None and phi is not None:
        bases = dossier/"bases.npz"
        autres = {} if reduction is None else dict(K_reduit=reduction.kr, M_reduit=reduction.mr,
            base_physique=reduction.base, indices_interieur=reduction.i, indices_ports=reduction.s)
        with bases.open("xb") as f:
            np.savez_compressed(f, B=b, Phi=phi, **autres)
        result["bases"] = dict(fichier=bases.name, sha256=sha(bases),
            b_sha256=hashlib.sha256(np.asarray(b, dtype="<f8").tobytes()).hexdigest(),
            phi_sha256=hashlib.sha256(np.asarray(phi, dtype="<f8").tobytes()).hexdigest(),
            forme=list(b.shape), sauvegarde_hors_chronometre=True)
    if reduction is not None:
        result["reduction"] = {nom: getattr(reduction, nom) for nom in
            ("taille_reduite", "taille_complement_reduit", "defaut_contrainte",
             "defaut_canonique", "resolutions_statiques", "statut", "historique")}
        result["directions"] = reduction.taille_complement_reduit
    if controle is not None:
        result["controle_uniforme"] = {nom: getattr(controle, nom) for nom in
            ("rho", "resolvante_marge", "ymax", "delta", "delta_colonnes", "termes",
             "defaut_masse", "defaut_deformation", "cm", "cd",
             "defaut_fonctionnel", "borne_schur", "preparation_s",
             "norme_w_m", "norme_w_d", "norme_h_m", "norme_h_d", "norme_theta")}
        if hasattr(controle, "facteurs"):
            result["controle_facteurs"] = controle.facteurs
    result["retours_par_frequence"] = retours
    # NaN des champs incomplets restent dans NPY, jamais remplacés par zéro dans le juge.
    result["normes_candidates"] = dict(indices=valides, masse=nm[valides], deformation=nd[valides])
    if len(valides) == 257:
        try:
            ref = np.load(reference, mmap_mode="r", allow_pickle=False)
            result["juge"] = juger(d, m, metrique, champs, ref)
            if controle is not None:
                result["controle_oracle"] = controle_oracle(d, m, champs, ref, retours)
        except Exception as exc:
            result.update(statut="refus_ou_echec", phase="juge",
                          erreur=dict(type=type(exc).__name__, message=str(exc)))
    result["sauvegarde_et_juge_hors_mesure_s"] = time.perf_counter()-hors
    if bib is not None:
        result['bibliotheque_avant'] = bib.identite
        result['bibliotheque_apres_sha256'] = sha(bib.identite['chemin'])
        if result['bibliotheque_apres_sha256'] != bib.identite['sha256']:
            raise RuntimeError('bibliothèque modifiée pendant le worker')
    result["sources_apres_sha256"] = sources()
    if result["sources_apres_sha256"] != result["sources_sha256"]:
        raise RuntimeError("sources modifiées pendant le worker")
    ecrire(dossier/"resultat.json", result)


def worker_temoin(entree, reference, variante, n, dossier):
    """Le worker historique reste intact, avec son propre chronomètre.

    Le crochet du juge retient seulement le tableau déjà calculé, après
    l'arrêt de ce chronomètre et la mesure mémoire. Aucun solveur remplacé.
    Le JSON historique est conservé séparément, sans réécriture.
    """
    import confronte_ports_exudyn as historique
    donnees(entree, reference, n)
    identite, empreintes = environnement(variante), sources()
    capture = []
    original = historique.juger

    def conserver(d, m, metrique, candidat, ref):
        capture.append(candidat)
        return original(d, m, metrique, candidat, ref)

    brut = dossier/"temoin-brut.json"
    historique.juger = conserver
    try:
        historique.worker(entree, reference, variante,
                          RANGS[n] if variante.startswith("hcb") else 0, brut)
    finally:
        historique.juger = original
    result = json.loads(brut.read_text())
    result.update(n=n, max_hz=40., type_certificat="aucun",
                  parametres_certificat=parametres_certificat(variante),
                  environnement=identite, sources_sha256=empreintes,
                  sources_apres_sha256=sources(), memoire_comparable=False,
                  resultat_historique=dict(fichier=brut.name, sha256=sha(brut)))
    if capture:
        result["champs"] = sauver_champs(dossier, capture[0], range(len(capture[0])))
    if result["sources_apres_sha256"] != empreintes:
        raise RuntimeError("sources modifiées pendant le témoin")
    ecrire(dossier/"resultat.json", result)


def analyser(rapport):
    if rapport.get("schema") != SCHEMA or rapport.get("type_campagne") != TYPE_CAMPAGNE:
        raise ValueError("schéma de campagne inattendu")
    if len(rapport["essais"]) != 36:
        raise ValueError("campagne incomplète : 36 essais requis")
    protocole = rapport["protocole"]
    if (protocole["variantes"] != list(VARIANTES)
            or protocole["certificats"] != {v: parametres_certificat(v) for v in CANDIDATS}
            or protocole["controles"] != {v: parametres_controle(v) for v in CONTROLES}
            or protocole["max_directions"] != 128
            or protocole["blocs"] != 4 or protocole["profondeur_controle"] != 8
            or protocole["directions_retenues"] != 1):
        raise ValueError("paramètres de campagne différents du protocole")
    attendus = {(n, v, p) for n in MAILLAGES for v in VARIANTES for p in range(4)}
    observes, versions, identites, lignes = set(), set(), {}, []
    for entree in rapport["essais"]:
        e = entree["resultat"]
        cle = (e["n"], e["variante"], entree["passage"])
        if cle in observes or cle not in attendus:
            raise ValueError("essai dupliqué ou imprévu")
        observes.add(cle)
        if (entree["role"], entree["repetition"]) != (
                "chauffe" if entree["passage"] == 0 else "mesure", entree["passage"]-1):
            raise ValueError("rôle incorrect")
        env = e["environnement"]
        if env["cpu"] != [8] or env["fils"] != THREADS:
            raise ValueError("CPU ou nombre de fils différents")
        if e["sources_sha256"] != rapport["sources_sha256"] or e["sources_apres_sha256"] != rapport["sources_sha256"]:
            raise ValueError("source différente pendant la campagne")
        cas = rapport["cas"][str(e["n"])]
        if e["entree_sha256"] != cas["entree_sha256"] or e["reference_sha256"] != cas["reference"]["fichier_sha256"]:
            raise ValueError("données ou oracle différents")
        if (e["entree_sha256"], e["reference_sha256"]) != DONNEES_FIGEES[e["n"]]:
            raise ValueError("identités historiques des données perdues")
        modes = RANGS[e["n"]] if e["variante"].startswith("hcb") else 1 if e["variante"] in CANDIDATS else 0
        if e["modes_demandes"] != modes:
            raise ValueError("rang demandé différent du protocole")
        parametres = parametres_certificat(e["variante"])
        if (e["type_certificat"] != parametres["type"]
                or e["parametres_certificat"] != parametres):
            raise ValueError("type ou paramètres du certificat différents du protocole")
        if e["variante"] in CANDIDATS:
            if (e["parametres_reduction"] != parametres_reduction(e["variante"])
                    or e["parametres_controle"] != parametres_controle(e["variante"])
                    or e["omega_max_binary64"] != OMEGA_MAX):
                raise ValueError("réduction ou bande différentes du protocole")
            if e["statut"] == "termine":
                cert = e["certificat_complement"]
                if (cert["certification_machine"] is not True
                        or not isinstance(cert["lambda_min"], (int, float))
                        or isinstance(cert["lambda_min"], bool)
                        or not math.isfinite(cert["lambda_min"])):
                    raise ValueError("certificat différent du protocole")
                dimension = 6*e["n"]-6
                if (cert["dimension_interieure"], cert["rang_contraintes"],
                        cert["dimension_complement"]) != (dimension, 1, dimension-1):
                    raise ValueError("dimensions du certificat incohérentes")
                if e["variante"] in INERTIE:
                    if (cert["lambda_min"] != GAMMA_INERTIE
                            or cert["preuve_kkt"]["signature"] != [dimension, 1, 0]
                            or cert["preuve_masse"]["signature"] != [dimension, 0, 0]
                            or sorted(cert["permutation_physique"]) != list(range(dimension))):
                        raise ValueError("seuil ou inertie différents du protocole")
                if e['variante']=='decimal_controle':
                    if cert['precision_decimal']!=32 or cert['permutation_physique']!=list(range(dimension)):
                        raise ValueError('témoin Decimal différent du protocole')
                else:
                    if (cert['arithmetique']!='binary64_nextafter_cpp'
                            or cert['encodage_intervalles']!='hex_binary64'
                            or e['bibliotheque_avant']['sha256']!=rapport['compilation']['bibliotheque_sha256']
                            or e['bibliotheque_apres_sha256']!=rapport['compilation']['bibliotheque_sha256']
                            or cert['bibliotheque']!=e['bibliotheque_avant']):
                        raise ValueError('bibliothèque ou arithmétique différente du protocole')
                marge = e["marge_bande_exacte"]
                fraction = Fraction(int(marge["numerateur"]), int(marge["denominateur"]))
                attendue = Fraction(cert["lambda_min"])-Fraction(e["omega_max_binary64"])**2
                if fraction != attendue or fraction <= 0:
                    raise ValueError("bande non démontrée exactement")
                selection = e["selection"]
                if (selection["k"], selection["graine"], selection["ncv"],
                        selection["tolerance"], selection["which"]) != (1, 107+e["n"], 20, 1e-11, "LA"):
                    raise ValueError("sélection modale différente du protocole")
                if e["variante"] in FACTEURS:
                    verifier_diagnostic_facteurs(e)
                elif "controle_facteurs" in e:
                    raise ValueError("diagnostic de facteurs dans le contrôle historique")
        famille = "candidat" if e["variante"] in CANDIDATS else "hcb" if e["variante"].startswith("hcb") else "lu"
        identite = json.dumps(env, sort_keys=True)
        if famille in identites and identites[famille] != identite:
            raise ValueError("environnement modifié dans une famille")
        identites[famille] = identite
        versions.add((env["python"], env["numpy"], env["scipy"]))
        qualifier(e)
    if observes != attendus or len(versions) != 1:
        raise ValueError("couverture ou versions numériques non comparables")
    for n in MAILLAGES:
        for v in VARIANTES:
            essais = [e for e in rapport["essais"] if e["resultat"]["n"] == n and e["resultat"]["variante"] == v]
            mesures = [e["resultat"] for e in essais if e["role"] == "mesure"]
            acceptes = [qualifier(e["resultat"]) for e in essais]
            ligne = dict(n=n, max_hz=40., variante=v, modes_demandes=(RANGS[n] if v.startswith("hcb") else 1 if v in CANDIDATS else 0),
                eligible_champs=all(acceptes), acceptations_champs=sum(qualifier(e) for e in mesures),
                chauffe_acceptee=next(qualifier(e["resultat"]) for e in essais if e["role"] == "chauffe"))
            if all(e["statut"] == "termine" for e in mesures):
                for nom in ("preparation_s", "reponses_s", "total_s"):
                    valeurs = [e[nom] for e in mesures]
                    ligne[nom] = statistics.median(valeurs)
                    ligne[nom+"_plage"] = [min(valeurs), max(valeurs)]
                ligne["erreur_mediane"] = statistics.median(max(*e["juge"]["maxima"].values(), *e["juge"]["maxima_operateurs"].values()) for e in mesures)
                ligne["erreur_pire"] = max(max(*e["juge"]["maxima"].values(), *e["juge"]["maxima_operateurs"].values()) for e in mesures)
            ligne["type_certificat"] = parametres_certificat(v)["type"]
            ligne["parametres_certificat"] = parametres_certificat(v)
            if v in CANDIDATS:
                ligne["parametres_reduction"] = parametres_reduction(v)
                ligne["parametres_controle"] = parametres_controle(v)
                if all(e["statut"] == "termine" for e in mesures):
                    for phase in ("condensation", "selection_direction", "certification_complement",
                                  "construction_krylov", "construction_controle", "reponses"):
                        valeurs = [e["phases_s"][phase] for e in mesures]
                        ligne.setdefault("phases_s", {})[phase] = statistics.median(valeurs)
                        ligne.setdefault("phases_s_plages", {})[phase] = [min(valeurs), max(valeurs)]
                    ligne["lambda_min_certifie"] = min(e["certificat_complement"]["lambda_min"] for e in mesures)
                    ligne["tailles_reduites"] = sorted({e["reduction"]["taille_reduite"] for e in mesures})
            if v in CONTROLES:
                controles = [qualifier_controle(e["resultat"]) for e in essais]
                ligne["controle_accepte_tous_essais"] = all(c[0] for c in controles)
                ligne["majorants_couvrent_oracles"] = all(c[1] for c in controles)
                ligne["eligible_controle"] = (ligne["eligible_champs"]
                    and ligne["controle_accepte_tous_essais"] and ligne["majorants_couvrent_oracles"])
                ligne["reponses_certifiees_machine"] = False
            lignes.append(ligne)
    return dict(cible=CIBLE, nombre_essais=36, type_campagne=TYPE_CAMPAGNE, configurations=lignes,
                portee="consoles natives, contrôle par facteurs et LU corrigée ; aucun nouveau classement Exudyn/MBDyn/Simpack")


def plan():
    return [(n, v, p) for n in MAILLAGES for p in range(4)
            for v in (VARIANTES if p % 2 == 0 else VARIANTES[::-1])]


def ecrire_manifest(dossier):
    ecrire(dossier/"manifest.json", dict(schema=1, date_utc=maintenant(),
        fichiers_sha256={str(p.relative_to(dossier)): sha(p) for p in sorted(dossier.rglob("*"))
                         if p.is_file() and p.name != "manifest.json"}))


def verifier(dossier):
    """SHA, protocole et qualifications enregistrées ; aucun rejeu de solveur.

    La preuve d'inertie n'est pas recalculée ici. Le vérificateur d'archive
    compacte séparé pourra la rejouer depuis D_i, M_ii et le B conservé.
    Les transformations des facteurs et les normes de réponse ne sont
    pas rejouées : leurs diagnostics permettent de qualifier les essais.
    """
    import numpy as np
    dossier = Path(dossier)
    manifest = json.loads((dossier/"manifest.json").read_text())["fichiers_sha256"]
    presents = {str(p.relative_to(dossier)) for p in dossier.rglob("*") if p.is_file() and p.name != "manifest.json"}
    if set(manifest) != presents:
        raise ValueError("inventaire incohérent")
    for nom, empreinte in manifest.items():
        if sha(dossier/nom) != empreinte:
            raise ValueError("fichier altéré : "+nom)
    rapport = json.loads((dossier/"rapport.json").read_text())
    bilan = analyser(rapport)
    if rapport["statut"] != "termine" or bilan != rapport["analyse"] or bilan != json.loads((dossier/"bilan.json").read_text()):
        raise ValueError("bilan incohérent")
    if [(e["resultat"]["n"], e["resultat"]["variante"], e["passage"]) for e in rapport["essais"]] != plan():
        raise ValueError("ordre de passage incorrect")
    for nom, empreinte in rapport["sources_sha256"].items():
        if sha(dossier/"sources"/nom) != empreinte:
            raise ValueError("source incohérente")
    if sha(dossier/'build/inertie_binaire.so') != rapport['compilation']['bibliotheque_sha256']:
        raise ValueError('bibliothèque archivée différente')
    for n, cas in rapport["cas"].items():
        verifier_reference(cas["reference"])
        if sha(dossier/"entrees"/f"n{n}-f40.npz") != cas["entree_sha256"]:
            raise ValueError("donnée incohérente")
        oracle = dossier/"oracles"/f"n{n}-f40.ref.npy"
        if oracle.exists() and sha(oracle) != cas["reference"]["fichier_sha256"]:
            raise ValueError("oracle archivé altéré")
    for e in rapport["essais"]:
        rep = dossier/"essais"/e["nom"]
        if sha(rep/"resultat.json") != e["resultat_sha256"] or json.loads((rep/"resultat.json").read_text()) != e["resultat"]:
            raise ValueError("essai incohérent")
        if sha(dossier/"essais"/(e["nom"]+".log")) != e["journal_sha256"]:
            raise ValueError("journal incohérent")
        terminal = json.loads((dossier/"essais"/(e["nom"]+".terminal.json")).read_text())
        if terminal != e["terminal"] or terminal["retour"] != 0 or terminal["attente_terminee"] is not True:
            raise ValueError("handle de processus non terminal")
        for cle in ("bases", "champs", "resultat_historique"):
            if cle in e["resultat"]:
                artefact = e["resultat"][cle]
                if sha(rep/artefact["fichier"]) != artefact["sha256"]:
                    raise ValueError("artefact de calcul incohérent")
        if "bases" in e["resultat"] and "certificat_complement" in e["resultat"]:
            if e["resultat"]["bases"]["b_sha256"] != e["resultat"]["certificat_complement"]["contraintes_sha256"]:
                raise ValueError("B archivé différent du B certifié")
        if "bases" in e["resultat"]:
            meta = e["resultat"]["bases"]
            with np.load(rep/meta["fichier"], allow_pickle=False) as z:
                for nom, cle in (("B", "b_sha256"), ("Phi", "phi_sha256")):
                    tableau = z[nom]
                    if tableau.shape != (6*e["resultat"]["n"]-6, 1) or not np.all(np.isfinite(tableau)):
                        raise ValueError("B/Phi archivé invalide")
                    if hashlib.sha256(np.asarray(tableau, dtype="<f8").tobytes()).hexdigest() != meta[cle]:
                        raise ValueError("coefficients de B/Phi incohérents")
        if "champs" in e["resultat"]:
            meta = e["resultat"]["champs"]
            champ = np.load(rep/meta["fichier"], mmap_mode="r", allow_pickle=False)
            valides = meta["indices_frequences_valides"]
            if champ.shape != (257, 6*e["resultat"]["n"], 6) or len(set(valides)) != len(valides):
                raise ValueError("champ archivé de dimensions incorrectes")
            if any(j not in range(257) for j in valides) or not np.all(np.isfinite(champ[valides])):
                raise ValueError("champ déclaré valide non fini")
            if e["resultat"]["statut"] == "termine" and valides != list(range(257)):
                raise ValueError("champ terminé incomplet")
    return bilan



def campagne(dossier,entrees,python):
    from inertie_binaire_compile import compiler
    dossier,entrees=Path(dossier).absolute(),Path(entrees).absolute()
    if dossier.exists():raise FileExistsError('dossier présent : aucune reprise')
    rapport=dict(schema=SCHEMA,type_campagne=TYPE_CAMPAGNE,date_utc=maintenant(),statut='en_cours',cpu=8,
        sources_sha256=sources(),cas={},essais=[],python=fichier(python),
        protocole=dict(maillages=list(MAILLAGES),max_hz=40.,frequences=257,charges=6,variantes=list(VARIANTES),
            certificats={v:parametres_certificat(v) for v in CANDIDATS},
            controles={v:parametres_controle(v) for v in CONTROLES},
            max_directions=128,blocs=4,profondeur_controle=8,directions_retenues=1,
            repetitions=3,echauffements=1,processus_frais=True,ordre_inverse_par_passage=True,
            comparaison_bande='Fraction(lambda)>Fraction(omega_max)**2',
            cout='QR, sélection, permutation complète si binaire, certificat, Krylov, contrôle par facteurs, champs et normes physiques',
            exclus='compilation unique et chargement du module, imports, lectures, sauvegardes, provenance et juge',
            reponses_certifiees_machine=False,reprise_automatique=False))
    for n in MAILLAGES:
        nom=f'n{n}-f40'
        cas=json.loads((entrees/(nom+'.npz.json')).read_text())
        cas['reference']=json.loads((entrees/(nom+'.ref.npy.json')).read_text());verifier_reference(cas['reference'])
        if (sha(entrees/(nom+'.npz')),sha(entrees/(nom+'.ref.npy')))!=DONNEES_FIGEES[n]:
            raise ValueError('entrée ou oracle différent du protocole')
        cas['oracle_original']=str(entrees/(nom+'.ref.npy'));rapport['cas'][str(n)]=cas
    dossier.mkdir(parents=True,exist_ok=False)
    for name in ('sources','entrees','essais'):(dossier/name).mkdir()
    for name,expected in rapport['sources_sha256'].items():
        shutil.copy2(CI/name,dossier/'sources'/name)
        if sha(dossier/'sources'/name)!=expected:raise ValueError('source modifiée pendant la capture')
    for n in MAILLAGES:
        for suffix in ('.npz','.npz.json','.ref.npy.json'):
            name=f'n{n}-f40{suffix}';shutil.copy2(entrees/name,dossier/'entrees'/name)
    lib=compiler(dossier/'build')
    rapport['compilation']=json.loads((dossier/'build/identite.json').read_text())
    if rapport['compilation']['source_sha256']!=rapport['sources_sha256']['inertie_binaire_native.cpp']:
        raise ValueError('source compilée différente de la capture')
    ecrire(dossier/'rapport.json',rapport)
    env=dict(os.environ,**THREADS,**FILS_SUPPLEMENTAIRES,PYTHONDONTWRITEBYTECODE='1',VINKULUM_INERTIE_BINARY64_LIB=str(lib))
    env.pop('PYTHONPATH',None);env.pop('PYTHONHOME',None)
    for n,variante,passage in plan():
        name=f'n{n}-f40-{variante}-passage{passage}';sortie=dossier/'essais'/name
        command=[str(python),str(dossier/'sources'/Path(__file__).name),'--worker',str(dossier/'entrees'/f'n{n}-f40.npz'),
            rapport['cas'][str(n)]['oracle_original'],variante,str(n),str(sortie)]
        print(name,flush=True);start=time.perf_counter();log=dossier/'essais'/(name+'.log')
        with log.open('x') as stream:
            process=subprocess.Popen(command,env=env,stdout=stream,stderr=subprocess.STDOUT);code=process.wait()
        terminal=dict(pid=process.pid,retour=code,attente_terminee=True,date_utc=maintenant(),processus_s=time.perf_counter()-start,commande=command)
        ecrire(dossier/'essais'/(name+'.terminal.json'),terminal)
        if code!=0 or not (sortie/'resultat.json').is_file():
            rapport.update(statut='interrompu',essai_interrompu=name,terminal=terminal)
            ecrire(dossier/'rapport.json',rapport,remplacer=True)
            raise RuntimeError('processus terminal en échec ; aucune reprise')
        resultat=json.loads((sortie/'resultat.json').read_text())
        rapport['essais'].append(dict(nom=name,passage=passage,role='chauffe' if passage==0 else 'mesure',repetition=passage-1,
            resultat=resultat,resultat_sha256=sha(sortie/'resultat.json'),journal_sha256=sha(log),terminal=terminal))
        ecrire(dossier/'rapport.json',rapport,remplacer=True)
    rapport.update(statut='termine',analyse=analyser(rapport),fin_utc=maintenant())
    ecrire(dossier/'rapport.json',rapport,remplacer=True);ecrire(dossier/'bilan.json',rapport['analyse'])
    ecrire_manifest(dossier);verifier(dossier)
    print('36 essais terminés et vérifiés : '+str(dossier),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--campagne',nargs=3,metavar=('DOSSIER_NEUF','ENTREES','PYTHON_GELE'))
    group.add_argument('--worker',nargs=5,metavar=('ENTREE','REFERENCE','VARIANTE','N','DOSSIER_NEUF'))
    group.add_argument('--verifier',type=Path)
    args=p.parse_args()
    if args.campagne or args.worker:
        if 8 not in os.sched_getaffinity(0):p.error('CPU 8 requis')
        os.sched_setaffinity(0,{8});os.environ.update(THREADS);os.environ.update(FILS_SUPPLEMENTAIRES)
    if args.campagne:campagne(*args.campagne)
    elif args.worker:
        entree,reference,variante,n,dossier=args.worker;n=int(n)
        if variante not in VARIANTES or n not in MAILLAGES:p.error('configuration imprévue')
        dossier=Path(dossier).absolute();dossier.mkdir(parents=True,exist_ok=False)
        (worker_candidat if variante in CANDIDATS else worker_temoin)(entree,reference,variante,n,dossier)
    else:verifier(args.verifier)


if __name__=='__main__':main()
