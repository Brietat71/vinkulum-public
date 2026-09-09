"""Archive compacte des 60 essais de contrôle par facteurs.

Les preuves d'inertie sont rejouées depuis les modèles et B exacts, avec
cache par identité ; les décisions physiques sont recalculées depuis les
détails sauvegardés. Aucun champ, oracle ou contrôleur n'est rejoué et
aucun contenu de sources/ dans l'archive n'est exécuté.
"""
import argparse
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
    ecrire_manifest,
)
from archive_inertie_contrainte import (
    _uniforme, _verifier_sources_rejeu, verifier_inertie,
)
from experience_krylov_contraint import DONNEES_FIGEES, MAILLAGES, RANGS, verifier_reference
from test_archive_complement_spectral import _exiger, _grille, _hash, _nombre, _sha, _juge

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT/"docs/bancs/controle-facteurs-2026"
CANDIDATS = ("facteurs_controle", "inertie_controle")
VARIANTES = CANDIDATS+("lu_corrigee", "hcb_standard", "hcb_energie")
PILOTE_SHA256 = "508efe7f731cfd3c368229b62706b16e32456ae22ae881571e0ad6a3a1f30a7b"
CONTROLE_SHA256 = "fed375d6b120df072d056a7029e84916dd666ab8bbfc926b3ccae849dd769bbe"


def _experience():
    # Module local connu ; aucune importation depuis destination/sources.
    import experience_controle_facteurs
    return experience_controle_facteurs


def plan():
    return [(n, v, p) for n in MAILLAGES for p in range(4)
            for v in (VARIANTES if p % 2 == 0 else VARIANTES[::-1])]


def comparaisons_admises(bilan):
    """Ratios du seul candidat nouveau, entre configurations toutes admises."""
    lignes = {(e["n"], e["variante"]): e for e in bilan["configurations"]}
    comparaisons = []
    for n in MAILLAGES:
        nouveau = lignes[(n, "facteurs_controle")]
        if nouveau.get("eligible_controle") is not True:
            continue
        for variante in VARIANTES[1:]:
            reference = lignes[(n, variante)]
            cle = "eligible_controle" if variante in CANDIDATS else "eligible_champs"
            if reference.get(cle) is not True:
                continue
            ratios = {k: _nombre(reference[k], "temps référence", positif=True)
                      / _nombre(nouveau[k], "temps candidat", positif=True)
                      for k in ("preparation_s", "reponses_s", "total_s")}
            comparaisons.append(dict(n=n, candidat="facteurs_controle",
                reference=variante, ratios_reference_sur_candidat=ratios))
    return comparaisons


def verifier_facteurs(e):
    """Audite les diagnostics déclarés ; les images physiques ne sont pas rejouées."""
    _experience().verifier_diagnostic_facteurs(e)
    f = e["controle_facteurs"]
    _exiger(set(f) == {"extension", "gamma_gram", "gamma_sommes",
             "certification_machine", "reparation"}, "inventaire des diagnostics de facteurs modifié")
    u = 2.**-53
    for cle, n in (("gamma_gram", 6*e["n"]), ("gamma_sommes", 7)):
        _exiger(f[cle] == n*u/(1-n*u), "budget d'arrondi de Gram différent des dimensions")
    for image in f["reparation"].values():
        _exiger(set(image) == {"rang_facteur", "colonnes_image", "defaut_colonnes",
                 "kappa", "certification_machine"}, "inventaire de réparation modifié")
        _exiger(type(image["rang_facteur"]) is int and type(image["colonnes_image"]) is int,
                 "dimensions entières des facteurs requises")
        _grille(image["defaut_colonnes"], (e["directions"],), "défauts de réparation")


def _candidat(destination, e, nom, omega, donnees, cache):
    _exiger(e["type_certificat"] == "inertie_kkt_dirigee", "type inertie incohérent")
    commun = dict(e)
    commun.pop("certificat_complement", None)
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
    if e["statut"] == "termine":
        _uniforme(e, omega)
        if e["variante"] == "facteurs_controle":
            verifier_facteurs(e)
        else:
            _exiger("controle_facteurs" not in e, "facteurs présents dans le contrôle historique")


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
             for e in r["essais"]] == plan(), "couverture/ordre des 60 essais incorrect")
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
    _exiger(r["schema"] == 3 and r["type_campagne"] == "controle_facteurs_60_essais" and r["statut"] == "termine" and r["cpu"] == 8,
             "campagne non terminée")
    _exiger(set(r["sources_sha256"]) == set(_experience().SOURCES)
             and r["sources_sha256"]["experience_controle_facteurs.py"] == PILOTE_SHA256
             and r["sources_sha256"]["controle_facteurs.py"] == CONTROLE_SHA256,
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
             for e in r["essais"]] == plan(), "couverture/ordre des 60 essais incorrect")
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
        champs = {(e["resultat"]["n"], e["resultat"]["variante"], e["passage"]):
                  e["resultat"].get("champs", {}).get("sha256") for e in r["essais"]}
        identiques = all(champs[(n, CANDIDATS[0], p)] is not None
                         and champs[(n, CANDIDATS[0], p)] == champs[(n, CANDIDATS[1], p)]
                         for n in MAILLAGES for p in range(4))
        audit.update(certificats_inertie_distincts_rejoues=len(rejeu),
                     champs_candidats_sha256_identiques=identiques,
                     comparaisons_admises=comparaisons_admises(bilan),
                     methodologie="rejeu dirigé depuis D/M/B ; aucun champ/oracle recalculé")
    return bilan


def archiver(source, destination=ARCHIVE, journal_parent=None, *, audit=None):
    """Crée une archive neuve après vérification SHA de la campagne terminale."""
    source, destination = Path(source), Path(destination)
    _exiger(not destination.exists(), "destination déjà présente")
    brut = _json((source/"manifest.json").read_bytes())
    r = _json((source/"rapport.json").read_bytes())
    _exiger(r["statut"] == "termine" and len(r["essais"]) == 60, "campagne incomplète")
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
