"""Archive compacte des 84 essais, avec rejeu des certificats d'inertie.

Rapport original gzip, résultats reconstruits octet pour octet, NPY de
champs/oracles exclus avec SHA. Aucun champ ni solveur physique n'est
rejoué. Les congruences dirigées d'inertie sont recalculées depuis D/M/B
pour chaque identité distincte ; jamais à partir de sources exécutées
dans l'archive. Les décisions de réponse restent des contrôles flottants.
"""
import argparse
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil

import numpy as np
from scipy.sparse import csr_matrix

from archive_krylov_contraint import (
    _candidat as _candidat_trace, _chemin, _ecrire, _environnement,
    _entete_npy, _inventaire, _json, _manifest, _original, _tableau,
    ecrire_manifest, verifier_controle)
from experience_krylov_contraint import DONNEES_FIGEES, MAILLAGES, RANGS, verifier_reference
from test_archive_complement_spectral import _decimal, _exiger, _grille, _hash, _nombre, _sha, _juge

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT/"docs/bancs/inertie-contraint-2026"
INERTIES = ("inertie_champs", "inertie_controle")
CANDIDATS = INERTIES+("krylov_champs", "krylov_controle")
VARIANTES = CANDIDATS+("lu_corrigee", "hcb_standard", "hcb_energie")
PILOTE_SHA256 = "26c08cdd54ff3ee338147a7e4620ba155f11f096b92c467727c8115cc0d89ee1"


def _experience():
    # Module local connu ; aucune importation depuis destination/sources.
    import experience_inertie_contrainte
    return experience_inertie_contrainte


def _verifier_sources_rejeu(rapport):
    import inspect
    import inertie_complement_dirigee as inertie
    import vinkulum._ports.certificat_spectral as arithmetique
    import vinkulum._ports.inverse_selectionnee as stockage
    # Le pilote courant peut recevoir une correction de journalisation ;
    # son instantané mesuré reste épinglé séparément. Le code effectivement
    # rejoué pour la preuve doit en revanche rester identique.
    for module, nom in ((inertie, "inertie_complement_dirigee.py"),):
        _exiger(_sha(inspect.getfile(module)) == rapport["sources_sha256"][nom],
                 "module local du rejeu différent du module mesuré")
    e = next(e["resultat"] for e in rapport["essais"] if e["resultat"]["variante"] in INERTIES)
    hashes = e["environnement"]["sources_python_sha256"]
    for module, nom in ((arithmetique, "_ports/certificat_spectral.py"),
                        (stockage, "_ports/inverse_selectionnee.py")):
        _exiger(_sha(inspect.getfile(module)) == hashes[nom],
                 "dépendance locale du rejeu différente de la dépendance mesurée")


def plan():
    return [(n, v, p) for n in MAILLAGES for p in range(4)
            for v in (VARIANTES if p % 2 == 0 else VARIANTES[::-1])]


def _intervalle(x):
    _exiger(isinstance(x, list) and len(x) == 2, "intervalle de pivot incomplet")
    a, b = (_decimal(v, "pivot dirigé") for v in x)
    _exiger(a <= b, "intervalle de pivot inversé")
    return a, b


def verifier_pivots(preuve, dimension, precision):
    """Recalcule les signatures et déterminants, indépendamment du LDL."""
    bas = Context(prec=precision, rounding=ROUND_FLOOR)
    haut = Context(prec=precision, rounding=ROUND_CEILING)
    indices, somme = [], [0, 0, 0]
    _exiger(isinstance(preuve["signature"], list) and len(preuve["signature"]) == 3
             and all(type(v) is int and v >= 0 for v in preuve["signature"]),
             "signature entière requise")
    for p in preuve["pivots"]:
        ids = p["indices"]
        _exiger(isinstance(ids, list) and len(ids) in (1, 2)
                 and all(type(i) is int and 0 <= i < dimension for i in ids), "indices de pivot invalides")
        indices.extend(ids)
        _exiger(isinstance(p["signature"], list) and len(p["signature"]) == 2
                 and all(type(v) is int and v >= 0 for v in p["signature"]),
                 "signature entière de pivot requise")
        if len(ids) == 1:
            a, b = _intervalle(p["diagonal"])
            attendu = [1, 0] if a > 0 else [0, 1] if b < 0 else None
        else:
            a, b, c = (_intervalle(p[k]) for k in ("a", "b", "c"))
            prod = (min(bas.multiply(x, y) for x in a for y in c),
                    max(haut.multiply(x, y) for x in a for y in c))
            carre = (Decimal(0) if b[0] <= 0 <= b[1]
                     else min(bas.multiply(x, x) for x in b),
                     max(haut.multiply(x, x) for x in b))
            det = (bas.subtract(prod[0], carre[1]), haut.subtract(prod[1], carre[0]))
            _exiger(_intervalle(p["determinant"]) == det, "déterminant dirigé du pivot incohérent")
            tr = bas.add(a[0], c[0]), haut.add(a[1], c[1])
            attendu = ([1, 1] if det[1] < 0 else [2, 0] if det[0] > 0 and tr[0] > 0
                       else [0, 2] if det[0] > 0 and tr[1] < 0 else None)
        _exiger(attendu is not None and p["signature"] == attendu,
                 "signature de pivot non démontrée")
        somme[0] += attendu[0]
        somme[1] += attendu[1]
    _exiger(len(indices) == dimension and sorted(indices) == list(range(dimension)),
             "élimination incomplète ou indice réutilisé")
    _exiger(preuve["signature"] == somme, "somme des signatures incohérente")
    for cle in ("largeur_max", "pic_coefficients", "coefficients_crees"):
        _exiger(type(preuve[cle]) is int and preuve[cle] >= 0, "compteur de preuve invalide")
    return somme


def _hors_temps(cert):
    return {k: v for k, v in cert.items() if k not in ("preparation_s", "phases_s")}


def verifier_inertie(cert, d, m, b, omega, parametres, cache):
    """Audite les congruences réelles ; cache par D/M/B et paramètres exacts."""
    from inertie_complement_dirigee import certifier_inertie_complement, _empreinte
    n, s = b.shape
    _exiger(parametres == _experience().parametres_certificat("inertie_champs"),
             "paramètres du certificat d'inertie différents du protocole")
    _exiger(all(type(cert[k]) is int for k in ("precision_decimal", "dimension_interieure",
             "rang_contraintes", "dimension_complement")), "dimensions/précision entières requises")
    _exiger(cert["certification_machine"] is True and cert["precision_decimal"] == 32
             and cert["dimension_interieure"] == n and cert["rang_contraintes"] == s
             and cert["dimension_complement"] == n-s, "dimensions/précision du certificat d'inertie incohérentes")
    _exiger(cert["portee"] == "coercivité stricte de D_i.T D_i-gamma M_ii sur ker(B.T) ; aucune réponse certifiée",
             "portée du certificat d'inertie modifiée")
    hb = hashlib.sha256(np.asarray(b, dtype="<f8").tobytes()).hexdigest()
    _exiger(cert["contraintes_sha256"] == hb and cert["d_sha256"] == _empreinte(d)
             and cert["masse_sha256"] == _empreinte(m), "empreinte B/D/M différente des entrées")
    gamma = _nombre(cert["lambda_min"], "gamma", positif=True)
    _exiger(gamma == parametres["gamma_binary64"] and gamma == float((2*np.pi*80.)**2),
             "seuil d'inertie différent de 80 Hz")
    _exiger(Fraction(gamma) > max(Fraction(float(w))**2 for w in omega), "bande non couverte en Fraction")
    _exiger(cert["permutation_physique"] == list(range(n)), "permutation physique différente du protocole")
    _exiger(verifier_pivots(cert["preuve_kkt"], n+s, 32) == [n, s, 0],
             "signature KKT différente de la coercivité requise")
    masse = cert["preuve_masse"]
    if masse["methode"] == "diagonale_binary64_positive":
        _nombre(masse["minimum_binary64"], "minimum massique", positif=True)
        _exiger(isinstance(masse["signature"], list) and all(type(v) is int for v in masse["signature"]),
                 "signature entière de masse requise")
        _exiger(m.nnz == n and np.array_equal(m.indices, np.arange(n))
                 and np.array_equal(m.indptr, np.arange(n+1)) and np.all(m.data > 0)
                 and masse["minimum_binary64"] == float(np.min(m.data))
                 and masse["signature"] == [n, 0, 0], "preuve de masse diagonale invalide")
    else:
        _exiger(masse["methode"] == "congruences_dirigees"
                 and verifier_pivots(masse, n, 32) == [n, 0, 0], "preuve de masse positive invalide")
    budgets = {k: parametres[k] for k in ("budget_operations", "budget_coefficients", "budget_rectangulaire")}
    _exiger(budgets == dict(budget_operations=100_000_000, budget_coefficients=2_000_000,
                           budget_rectangulaire=2_000_000), "budgets d'inertie différents du protocole")
    _exiger(type(cert["operations_decimal"]) is int and 0 < cert["operations_decimal"] <= budgets["budget_operations"]
             and cert["preuve_kkt"]["pic_coefficients"] <= budgets["budget_coefficients"]
             and cert["coefficients_contraintes"] == n*s <= budgets["budget_rectangulaire"],
             "compteurs du certificat hors budget")
    _nombre(cert["preparation_s"], "temps certificat", positif=True)
    _exiger(set(cert["phases_s"]) == {"masse", "assemblage", "elimination"}, "phases du certificat incomplètes")
    for v in cert["phases_s"].values():
        _nombre(v, "phase certificat")
    _exiger(abs(sum(cert["phases_s"].values())-cert["preparation_s"]) < 1e-9,
             "total des phases du certificat incohérent")
    cle = (cert["d_sha256"], cert["masse_sha256"], hb, gamma, 32, tuple(budgets.values()))
    if cle not in cache:
        calcule = certifier_inertie_complement(d, m, b, gamma, precision=32, **budgets)
        cache[cle] = _hors_temps(calcule)
    _exiger(_hors_temps(cert) == cache[cle], "certificat différent du rejeu dirigé depuis D/M/B")


def _uniforme(e, omega):
    u = e["controle_uniforme"]
    for cle, valeur in u.items():
        if cle not in ("termes", "delta_colonnes"):
            _nombre(valeur, "contrôle uniforme "+cle)
    termes = u["termes"]
    _exiger([v["degre"] for v in termes] == list(range(1, 8)), "degrés d'enveloppe incomplets")
    for terme in termes:
        for cle in ("polynome", "queue", "majorant"):
            _nombre(terme[cle], "enveloppe "+cle)
        _grille(terme["majorants_colonnes"], (7,), "enveloppe par colonne")
        _exiger(terme["majorant"] == terme["polynome"]+terme["queue"], "somme polynôme/queue incohérente")
    _exiger(u["delta"] == min(t["majorant"] for t in termes)
             and u["delta_colonnes"] == np.min([t["majorants_colonnes"] for t in termes], axis=0).tolist(),
             "minimum des enveloppes incohérent")
    lam = e["certificat_complement"]["lambda_min"]
    alpha = lam-float(omega[-1])**2
    attendus = dict(rho=float(omega[-1])**2/lam,
        resolvante_marge=1-u["rho"]*u["norme_theta"],
        cm=u["delta"]/alpha+u["defaut_masse"],
        cd=math.sqrt(lam)*u["delta"]/alpha+u["defaut_deformation"],
        borne_schur=u["delta"]**2/alpha+u["defaut_fonctionnel"])
    _exiger(u["resolvante_marge"] > 0 and all(math.isclose(u[k], v, rel_tol=64*2.**-52, abs_tol=0)
             for k, v in attendus.items()), "identités du contrôle uniforme incohérentes")
    verifier_controle(e)


def _candidat(destination, e, nom, omega, donnees, cache):
    if e["variante"] not in INERTIES:
        _exiger(e["type_certificat"] == "trace_inverse_selectionnee_dirigee", "type trace incohérent")
        _candidat_trace(destination, e, nom, omega)
        return
    _exiger(e["type_certificat"] == "inertie_kkt_dirigee", "type inertie incohérent")
    commun = dict(e)
    commun.pop("certificat_complement", None)
    if e["variante"] == "inertie_champs":
        commun["variante"] = "krylov_champs"
    # Réutilisation des contrôles de fréquences, normes et chronos ; la
    # preuve spectrale d'inertie est auditée séparément ci-dessous.
    _candidat_trace(destination, commun, nom, omega)
    if "certificat_complement" not in e:
        return
    ni = 6*e["n"]-6
    meta = e["bases"]
    with np.load(destination/"essais"/nom/"bases.npz", allow_pickle=False) as z:
        b, phi = z["B"], z["Phi"]
        _tableau(b, (ni, 1), "B")
        _tableau(phi, (ni, 1), "Phi")
        for cle, a in (("b_sha256", b), ("phi_sha256", phi)):
            _exiger(hashlib.sha256(np.asarray(a, dtype="<f8").tobytes()).hexdigest() == _hash(meta[cle]),
                     "empreinte B/Phi différente")
        _exiger(meta["forme"] == [ni, 1] and meta["sauvegarde_hors_chronometre"] is True,
                 "métadonnées des bases incohérentes")
        if "reduction" in e:
            nr = e["reduction"]["taille_reduite"]
            _exiger(nr == 7+e["directions"] and e["directions"] == e["reduction"]["taille_complement_reduit"],
                     "dimensions réduites incohérentes")
            _exiger(set(z.files) == {"B", "Phi", "K_reduit", "M_reduit", "base_physique",
                     "indices_interieur", "indices_ports"}, "inventaire des bases incomplet")
            for cle in ("K_reduit", "M_reduit"):
                _tableau(z[cle], (nr, nr), cle)
            _tableau(z["base_physique"], (ni+6, nr), "base physique")
            _exiger(np.array_equal(z["indices_interieur"], np.arange(ni))
                     and np.array_equal(z["indices_ports"], np.arange(ni, ni+6)), "partition physique différente")
        cert = e["certificat_complement"]
        marge = Fraction(cert["lambda_min"])-Fraction(float(omega[-1]))**2
        declaree = e["marge_bande_exacte"]
        _exiger(marge > 0 and marge == Fraction(int(declaree["numerateur"]), int(declaree["denominateur"])),
                 "marge de bande Fraction incohérente")
        d, m = donnees
        verifier_inertie(cert, d, m, b, omega, e["parametres_certificat"], cache)
    if e["variante"] == "inertie_controle" and e["statut"] == "termine":
        _uniforme(e, omega)


def verifier_complet(destination=ARCHIVE, *, audit=None):
    """Vérifie l'archive autonome ; ne lit aucun chemin des fichiers exclus."""
    destination = Path(destination)
    rejeu, donnees = {}, {}
    fichiers = _manifest(destination)
    provenance = _json((destination/"provenance-fichiers.json").read_bytes())
    brut = _json((destination/"manifeste-campagne-brute.json").read_bytes())
    exclus = _json((destination/"exclus.json").read_bytes())
    r = _json(_original(destination, provenance["rapport.json"]))
    _exiger([(e["resultat"]["n"], e["resultat"]["variante"], e["passage"])
             for e in r["essais"]] == plan(), "couverture/ordre des 84 essais incorrect")
    essais = {e["nom"]: e for e in r["essais"]}
    originaux = {}
    for nom, meta in provenance.items():
        _chemin(nom)
        originaux[nom] = _original(destination, meta, essais)
        _exiger(meta["sha256"] == brut["fichiers_sha256"][nom],
                 "empreinte différente du manifeste de campagne")
    champs_exclus = {n for n, v in exclus.items() if v["type"] == "champ"}
    _exiger(set(provenance).isdisjoint(champs_exclus)
             and set(provenance) | champs_exclus == set(brut["fichiers_sha256"]),
             "inventaire brut incomplet ou exclusion non déclarée")
    attendus = {v["archive"] for v in provenance.values()} | {
        "manifeste-campagne-brute.json", "provenance-fichiers.json", "exclus.json"}
    _exiger(set(fichiers)-attendus <= {"README.md", "journal-parent.log.gz"}
             and attendus <= set(fichiers), "fichiers compacts inattendus/manquants")
    _exiger(r["schema"] == 2 and r["type_campagne"] == "inertie_contrainte_84_essais" and r["statut"] == "termine" and r["cpu"] == 8,
             "campagne non terminée")
    _exiger(set(r["sources_sha256"]) == set(_experience().SOURCES)
             and r["sources_sha256"]["experience_inertie_contrainte.py"] == PILOTE_SHA256,
             "pilote différent du protocole figé")
    for nom, h in r["sources_sha256"].items():
        _exiger(hashlib.sha256(originaux["sources/"+nom]).hexdigest() == _hash(h),
                 "source capturée différente")
    _verifier_sources_rejeu(r)
    p = r["protocole"]
    strict = dict(maillages=list(MAILLAGES), variantes=list(VARIANTES), max_hz=40.,
        frequences=257, charges=6, blocs=4, directions_retenues=1, profondeur_controle=8,
        max_directions=128, repetitions=3, echauffements=1, processus_frais=True,
        ordre_inverse_par_passage=True, reponses_certifiees_machine=False, reprise_automatique=False)
    _exiger(all(type(p[k]) is type(v) and p[k] == v for k, v in strict.items()),
             "protocole compact modifié")
    _exiger([(e["resultat"]["n"], e["resultat"]["variante"], e["passage"])
             for e in r["essais"]] == plan(), "couverture/ordre des 84 essais incorrect")
    _exiger(set(r["cas"]) == {str(n) for n in MAILLAGES}, "maillages manquants")
    historique = _json(originaux["historique-hcb.json"])
    hm = _json(originaux["historique-hcb-manifest.json"])
    for fichier, cle in (("historique-hcb.json", "historique_hcb"),
                         ("historique-hcb-manifest.json", "manifeste_historique_hcb")):
        _exiger(hashlib.sha256(originaux[fichier]).hexdigest() == r[cle]["sha256"]
                 and len(originaux[fichier]) == r[cle]["octets"], "provenance HCB incohérente")
    _exiger(hm["fichiers_sha256"]["bilan.json"] == r["historique_hcb"]["sha256"],
             "bilan HCB différent du manifeste historique")
    rangs = {(x["n"], x["variante"]): x["modes"] for x in historique["configurations"] if x["max_hz"] == 40}
    for n in MAILLAGES:
        _exiger(all(rangs.get((n, v)) == RANGS[n] for v in ("hcb_standard", "hcb_energie")),
                 "rangs historiques HCB incohérents")
    omegas = {}
    exclusions_attendues = set()
    for ns, c in r["cas"].items():
        n = int(ns)
        prefixe = f"entrees/n{n}-f40"
        sidecar = _json(originaux[prefixe+".npz.json"])
        ref = _json(originaux[prefixe+".ref.npy.json"])
        _exiger({k: v for k, v in c.items() if k not in ("reference", "oracle_original")} == sidecar
                 and c["reference"] == ref, "métadonnées d'entrée différentes")
        verifier_reference(ref)
        _exiger((c["entree_sha256"], ref["fichier_sha256"]) == DONNEES_FIGEES[n],
                 "entrée/oracle historique différent")
        _exiger(hashlib.sha256(originaux[prefixe+".npz"]).hexdigest() == c["entree_sha256"],
                 "empreinte entrée NPZ différente")
        with np.load(destination/(prefixe+".npz"), allow_pickle=False) as z:
            omega = z["omega"]
            _tableau(omega, (257,), "omega")
            _exiger(np.array_equal(omega, 2*np.pi*np.linspace(0., 40., 257)), "bande NPZ modifiée")
            omegas[n] = omega.copy()
            matrices = []
            for cle in ("d", "m"):
                a = csr_matrix((z[cle+"_data"], z[cle+"_indices"], z[cle+"_indptr"]),
                               shape=tuple(z[cle+"_shape"]))
                _exiger(a.shape == (6*n, 6*n), "dimensions des modèles incorrectes")
                a.eliminate_zeros()
                matrices.append(a)
            ni = 6*n-6
            donnees[n] = (matrices[0][:, :ni].tocsr(), matrices[1][:ni, :ni].tocsr())
        nom = f"oracles/n{n}-f40.ref.npy"
        exclusions_attendues.add(nom)
        x = exclus[nom]
        _exiger(x["type"] == "oracle" and x["sha256"] == ref["fichier_sha256"]
                 and x["original"] == c["oracle_original"] and x["forme"] == [257, 6*n, 6],
                 "oracle exclu incohérent")
    for entree in r["essais"]:
        e, nom = entree["resultat"], entree["nom"]
        _exiger(nom == f"n{e['n']}-f40-{e['variante']}-passage{entree['passage']}", "nom d'essai incohérent")
        prefixe = "essais/"+nom
        _exiger(_json(originaux[prefixe+"/resultat.json"]) == e
                 and provenance[prefixe+"/resultat.json"]["sha256"] == entree["resultat_sha256"],
                 "résultat original différent du rapport")
        _exiger(provenance[prefixe+".log"]["sha256"] == entree["journal_sha256"], "journal différent")
        t = _json(originaux[prefixe+".terminal.json"])
        _exiger(t == entree["terminal"] and t["retour"] == 0 and t["attente_terminee"] is True,
                 "processus non terminal")
        _nombre(t["processus_s"], "temps processus", positif=True)
        _environnement(e)
        _exiger(e["cible"] == 1e-6 and e["max_hz"] == 40., "cible/bande modifiée")
        if e["statut"] == "termine":
            _juge(e["juge"])
        else:
            _exiger("erreur" in e or "refus" in e or e.get("refus_calcul"), "échec sans diagnostic")
        for cle in ("bases", "resultat_historique"):
            if cle in e:
                a = e[cle]
                _exiger(provenance[prefixe+"/"+_chemin(a["fichier"])]["sha256"] == a["sha256"],
                         "artefact différent du résultat")
        if "champs" in e:
            a = e["champs"]
            chemin = prefixe+"/"+_chemin(a["fichier"])
            exclusions_attendues.add(chemin)
            x = exclus[chemin]
            _exiger(x["type"] == "champ" and x["sha256"] == a["sha256"]
                     == brut["fichiers_sha256"][chemin]
                     and x["forme"] == a["forme"] == [257, 6*e["n"], 6]
                     and x["indices_frequences_valides"] == a["indices_frequences_valides"],
                     "champ exclu différent du résultat")
            indices = a["indices_frequences_valides"]
            _exiger(indices == sorted(set(indices)) and all(type(j) is int and 0 <= j < 257 for j in indices),
                     "indices de champ invalides")
            if e["statut"] == "termine":
                _exiger(indices == list(range(257)), "champ terminé incomplet")
        if e["variante"] in CANDIDATS:
            _candidat(destination, e, nom, omegas[e["n"]], donnees[e["n"]], rejeu)
    _exiger(set(exclus) == exclusions_attendues, "inventaire des exclusions incohérent")
    for x in exclus.values():
        _hash(x["sha256"])
        _exiger(x["dtype"] == "float64" and x["octets"] >= math.prod(x["forme"])*8,
                 "taille/type de fichier exclu incohérent")
    bilan = _experience().analyser(r)
    _exiger(bilan == r["analyse"] == _json(originaux["bilan.json"]), "bilan requalifié différent")
    if audit is not None:
        audit.update(certificats_inertie_distincts_rejoues=len(rejeu),
                     methodologie="rejeu dirigé depuis D/M/B ; aucun champ/oracle recalculé")
    return bilan


def archiver(source, destination=ARCHIVE, journal_parent=None, *, audit=None):
    """Crée une archive neuve après vérification SHA de la campagne terminale."""
    source, destination = Path(source), Path(destination)
    _exiger(not destination.exists(), "destination déjà présente")
    brut = _json((source/"manifest.json").read_bytes())
    r = _json((source/"rapport.json").read_bytes())
    _exiger(r["statut"] == "termine" and len(r["essais"]) == 84, "campagne incomplète")
    fichiers = _inventaire(source)
    _exiger(set(fichiers) == set(brut["fichiers_sha256"]), "inventaire brut incohérent")
    for n, p in fichiers.items():
        _chemin(n)
        _exiger(_sha(p) == brut["fichiers_sha256"][n], "campagne brute altérée : "+n)
    destination.mkdir(parents=True)
    shutil.copyfile(source/"manifest.json", destination/"manifeste-campagne-brute.json")
    provenance, exclus = {}, {}
    essais = {e["nom"]: e for e in r["essais"]}
    champs = {"essais/"+e["nom"]+"/"+e["resultat"]["champs"]["fichier"]: e["resultat"]["champs"]
              for e in r["essais"] if "champs" in e["resultat"]}
    for n, p in fichiers.items():
        if p.suffix == ".npy":
            _exiger(n in champs, "NPY inattendu dans la campagne")
            shape, fortran, dtype = _entete_npy(p)
            exclus[n] = dict(type="champ", original=str(p.absolute()),
                sha256=brut["fichiers_sha256"][n], octets=p.stat().st_size,
                forme=list(shape), dtype=str(dtype),
                indices_frequences_valides=champs[n]["indices_frequences_valides"])
            continue
        compresser = n == "rapport.json" or n.endswith(("/resultat.json", "/temoin-brut.json", ".log"))
        if n.endswith("/resultat.json"):
            nom_essai = p.parent.name
            data = (json.dumps(essais[nom_essai]["resultat"], ensure_ascii=False,
                               indent=2, allow_nan=False)+"\n").encode()
            _exiger(hashlib.sha256(data).hexdigest() == brut["fichiers_sha256"][n],
                     "reconstruction du résultat non identique")
            provenance[n] = dict(archive="rapport.json.gz", encodage="extrait_rapport",
                essai=nom_essai, sha256=brut["fichiers_sha256"][n], octets=len(data))
            continue
        cible = n+".gz" if compresser else n
        out = destination/cible
        out.parent.mkdir(parents=True, exist_ok=True)
        data = p.read_bytes()
        out.write_bytes(gzip.compress(data, compresslevel=9, mtime=0) if compresser else data)
        provenance[n] = dict(archive=cible, encodage="gzip" if compresser else "identique",
                             sha256=brut["fichiers_sha256"][n], octets=len(data))
    for ns, c in r["cas"].items():
        path = Path(c["oracle_original"])
        h = _sha(path)
        _exiger(h == c["reference"]["fichier_sha256"], "oracle original altéré")
        shape, fortran, dtype = _entete_npy(path)
        exclus[f"oracles/n{ns}-f40.ref.npy"] = dict(type="oracle", original=str(path),
            sha256=h, octets=path.stat().st_size, forme=list(shape), dtype=str(dtype))
    _ecrire(destination/"provenance-fichiers.json", provenance)
    _ecrire(destination/"exclus.json", exclus)
    if journal_parent is not None:
        (destination/"journal-parent.log.gz").write_bytes(
            gzip.compress(Path(journal_parent).read_bytes(), mtime=0))
    ecrire_manifest(destination)
    return verifier_complet(destination, audit=audit)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", nargs="?", type=Path, default=ARCHIVE)
    parser.add_argument("--creer-depuis", type=Path)
    parser.add_argument("--verifier", type=Path)
    parser.add_argument("--journal-parent", type=Path)
    args = parser.parse_args()
    if args.verifier is not None:
        args.destination = args.verifier
    audit = {}
    resultat = (archiver(args.creer_depuis, args.destination, args.journal_parent, audit=audit)
                if args.creer_depuis else verifier_complet(args.destination, audit=audit))
    print(json.dumps(dict(archive=str(args.destination), nombre_essais=resultat["nombre_essais"],
                         statut="archive compacte vérifiée", audit=audit), ensure_ascii=False))
