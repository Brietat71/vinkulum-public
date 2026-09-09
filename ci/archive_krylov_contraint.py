"""Archive compacte et validation de la campagne Krylov contrainte.

Le rapport original est comprimé sans modification ; les résultats qui y
figurent sont reconstruits avec le sérialiseur JSON standard, puis leur SHA
est comparé à celui des octets originaux. Le manifeste brut lie ces
empreintes aux NPY exclus. La validation ne charge aucun champ, ne
rejoue aucun solveur et n'exécute jamais les sources capturées. Elle contrôle
les détails déclarés, les décisions, les identités et l'arithmétique finale
des certificats ; elle ne recrée ni les intervalles U/H/J ni les oracles.
"""
import argparse
from copy import deepcopy
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import shutil

import numpy as np

from experience_krylov_contraint import (CANDIDATS, DONNEES_FIGEES,
    FILS_SUPPLEMENTAIRES, MAILLAGES, RANGS, SOURCES, THREADS, VARIANTES,
    analyser, plan, qualifier_controle, verifier_reference)
from test_archive_complement_spectral import (
    _certificat, _exiger, _grille, _hash, _juge, _nombre, _sha)

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "docs/bancs/krylov-contraint-2026"
PILOTE_SHA256 = "eca0ad55fcc208fa16e5b9a701be34782f79f318040544b1c415a296ba4f744e"


def _json(data):
    def constant(x):
        raise ValueError("constante JSON non finie : " + x)
    def paires(items):
        d = {}
        for k, v in items:
            _exiger(k not in d, "clé JSON dupliquée : " + k)
            d[k] = v
        return d
    return json.loads(data, parse_constant=constant, object_pairs_hook=paires)


def _ecrire(path, objet):
    Path(path).write_text(json.dumps(objet, ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n")


def _chemin(nom):
    _exiger(isinstance(nom, str) and nom and "\\" not in nom,
             "chemin d'archive invalide")
    p = PurePosixPath(nom)
    _exiger(not p.is_absolute() and ".." not in p.parts and str(p) == nom,
             "chemin d'archive non relatif/canonique")
    return nom


def _inventaire(destination):
    fichiers = {}
    for p in sorted(Path(destination).rglob("*")):
        _exiger(not p.is_symlink(), "lien symbolique interdit dans l'archive")
        if p.is_file() and p != Path(destination)/"manifest.json":
            fichiers[str(p.relative_to(destination))] = p
    return fichiers


def ecrire_manifest(destination):
    """Recalcule uniquement le manifeste compact (utile après ajout du README)."""
    fichiers = _inventaire(destination)
    _ecrire(Path(destination)/"manifest.json", dict(schema=1,
        fichiers_sha256={n: _sha(p) for n, p in fichiers.items()},
        fichiers_octets={n: p.stat().st_size for n, p in fichiers.items()}))


def _manifest(destination):
    m = _json((destination/"manifest.json").read_bytes())
    fichiers = _inventaire(destination)
    _exiger(m["schema"] == 1 and set(m["fichiers_sha256"]) == set(fichiers)
             == set(m["fichiers_octets"]), "inventaire compact incohérent")
    for n, p in fichiers.items():
        _chemin(n)
        _exiger(_sha(p) == _hash(m["fichiers_sha256"][n])
                 and p.stat().st_size == m["fichiers_octets"][n],
                 "fichier compact altéré : " + n)
    return fichiers


def _original(destination, meta, essais=None):
    path = destination/_chemin(meta["archive"])
    _exiger(meta["encodage"] in ("identique", "gzip", "extrait_rapport"), "encodage inconnu")
    if meta["encodage"] == "extrait_rapport":
        _exiger(meta["archive"] == "rapport.json.gz" and essais is not None,
                 "source de l'extrait incorrecte")
        data = (json.dumps(essais[meta["essai"]]["resultat"], ensure_ascii=False,
                           indent=2, allow_nan=False)+"\n").encode()
    else:
        data = gzip.decompress(path.read_bytes()) if meta["encodage"] == "gzip" else path.read_bytes()
    _exiger(len(data) == meta["octets"] and hashlib.sha256(data).hexdigest()
             == _hash(meta["sha256"]), "octets originaux non conservés")
    return data


def _tableau(a, forme, nom):
    _exiger(a.shape == forme and a.dtype == np.dtype("float64")
             and bool(np.all(np.isfinite(a))), "NPZ " + nom + " invalide")


def _entete_npy(path):
    with Path(path).open("rb") as f:
        version = np.lib.format.read_magic(f)
        _exiger(version in ((1, 0), (2, 0)), "version NPY inattendue")
        lecteur = (np.lib.format.read_array_header_1_0 if version == (1, 0)
                   else np.lib.format.read_array_header_2_0)
        return lecteur(f)


def verifier_controle(e):
    """Requalifie les bornes et recoupe les deux juges, sans champ brut."""
    qualifier_controle(e)
    c = e["controle_oracle"]
    refus = []
    for j, retour in enumerate(e["retours_par_frequence"]):
        sigma = _nombre(retour["sigma_min_bloc_conserve"], "sigma")
        schur = _nombre(retour["borne_schur"], "borne Schur")
        marge = retour["marge"]
        _exiger(type(marge) in (int, float) and math.isfinite(marge)
                 and marge == sigma-schur, "marge globale incohérente")
        _exiger(retour["certification_machine"] is False
                 and retour["taille_conservee"] == 7,
                 "portée du contrôle modifiée")
        for cle in ("residu_coefficients", "erreur_coefficients"):
            _nombre(retour[cle], cle)
        for nom in ("masse", "deformation"):
            b = retour["bornes"][nom]
            normes = _grille(b["normes"], (6,), "normes candidates")
            _exiger(e["normes_candidates"][nom][j] == b["normes"],
                     "normes candidates différentes du retour")
            _exiger((b["absolues"] is not None) == (marge > 0),
                     "majorant absolu incompatible avec la marge")
            if b["absolues"] is not None:
                bornes = _grille(b["absolues"], (6,), "majorants absolus")
                relatives = [float(v/(n-v)) if n > v else None
                             for v, n in zip(bornes, normes, strict=True)]
                _exiger(b["relatives"] == relatives, "majorant relatif incohérent")
            else:
                _exiger(b["relatives"] == [None]*6, "majorants relatifs sans marge")
        if not c["acceptations"][j]:
            refus.append(j)
    _exiger([r["frequence"] for r in c["refus"]] == refus
             and all(isinstance(r["raisons"], list) and r["raisons"] for r in c["refus"]),
             "refus du contrôle incomplets")
    for nom in ("masse", "deformation"):
        err = _grille(c["erreurs_absolues"][nom], (257, 6), "erreurs absolues")
        ref = _grille(c["normes_reference"][nom], (257, 6), "normes référence", positif=True)
        rel = np.asarray(e["juge"]["erreurs"][nom])
        # Les deux juges emploient des ordres d'opérations distincts pour M.
        _exiger(np.allclose(err/ref, rel, rtol=256*2.**-52, atol=0),
                 "erreurs absolues et juge relatif incompatibles")


def _environnement(e):
    env = e["environnement"]
    _exiger(env["fils"] == THREADS and env["fils_supplementaires"] == FILS_SUPPLEMENTAIRES
             and env["cpu"] == [8] and e["memoire_comparable"] is False,
             "protocole environnement/mémoire incohérent")
    _hash(env["executable"]["sha256"])
    _nombre(env["executable"]["octets"], "taille Python", positif=True)
    for d in env["distributions"].values():
        _hash(d["record_sha256"])
        if "roue_sha256" in d:
            _hash(d["roue_sha256"])
        for h in d.get("extensions_sha256", {}).values():
            _hash(h)
    for h in env["bibliotheques_numeriques_sha256"].values():
        _hash(h)
    if env["threadpoolctl"] is not None:
        _exiger(all(x["num_threads"] == 1 for x in env["threadpoolctl"]),
                 "bibliothèque numérique multithread")
    if e["variante"] != "lu_corrigee":
        for cle in ("record_sha256", "roue_sha256"):
            _hash(env["distribution"][cle])
        _exiger(env["distribution"]["version"] == env["moteur"], "version moteur incohérente")
    for cle in ("sources_python_sha256", "extensions_sha256"):
        for h in env.get(cle, {}).values():
            _hash(h)
    if "extension_sha256" in env:
        _hash(env["extension_sha256"])


def _candidat(destination, e, nom, omega):
    _exiger(e["certification_machine_reponses"] is False, "réponses déclarées certifiées")
    for v in e["phases_s"].values():
        _nombre(v, "phase chronométrée")
    _exiger(abs(sum(e["phases_s"].values())-e["temps_calcul_s"]) < 1e-9,
             "somme des phases incohérente")
    if e["statut"] == "termine":
        _exiger(len(e["retours_par_frequence"]) == 257
                 and e["frequences_calculees"] == 257 and e["refus_calcul"] == []
                 and e["normes_candidates"]["indices"] == list(range(257)),
                 "fréquences candidates incomplètes")
        _exiger(e["reponses_s"] == e["phases_s"]["reponses"]
                 and e["total_s"] == e["temps_calcul_s"], "temps candidat incohérents")
        for nom_norme in ("masse", "deformation"):
            _grille(e["normes_candidates"][nom_norme], (257, 6), "normes candidates")
        if e["variante"] == "krylov_champs":
            for retour in e["retours_par_frequence"]:
                _exiger(retour["certification_machine"] is False
                         and retour["taille_reduite"] == e["reduction"]["taille_reduite"]
                         and retour["taille_complement_reduit"] == e["directions"],
                         "dimensions/portée du retour Krylov incohérentes")
                residu = np.asarray(retour["residu_reduit"], dtype=object)
                _exiger(residu.shape == (retour["taille_reduite"], 6)
                         and all(type(v) in (int, float) and math.isfinite(v) for v in residu.flat),
                         "résidu réduit incomplet/non fini")
                _nombre(retour["defaut_contrainte"], "défaut de contrainte")
    if "certificat_complement" in e:
        meta = e["bases"]
        with np.load(destination/"essais"/nom/"bases.npz", allow_pickle=False) as z:
            b, phi = z["B"], z["Phi"]
            ni = 6*e["n"]-6
            _tableau(b, (ni, 1), "B")
            _tableau(phi, (ni, 1), "Phi")
            for cle, a in (("b_sha256", b), ("phi_sha256", phi)):
                _exiger(hashlib.sha256(np.asarray(a, dtype="<f8").tobytes()).hexdigest()
                         == _hash(meta[cle]), "empreinte B/Phi différente")
            _exiger(meta["forme"] == [ni, 1] and meta["sauvegarde_hors_chronometre"] is True,
                     "métadonnées des bases incohérentes")
            if "reduction" in e:
                nr = e["reduction"]["taille_reduite"]
                _exiger(nr == 7+e["directions"] and e["directions"]
                         == e["reduction"]["taille_complement_reduit"], "dimensions réduites incohérentes")
                _exiger(set(z.files) == {"B", "Phi", "K_reduit", "M_reduit", "base_physique",
                         "indices_interieur", "indices_ports"}, "inventaire des bases incomplet")
                for cle in ("K_reduit", "M_reduit"):
                    _tableau(z[cle], (nr, nr), cle)
                _tableau(z["base_physique"], (6*e["n"], nr), "base physique")
                _exiger(np.array_equal(z["indices_interieur"], np.arange(ni))
                         and np.array_equal(z["indices_ports"], np.arange(ni, ni+6)),
                         "partition des coordonnées physiques incorrecte")
            cert = deepcopy(e["certificat_complement"])
            _exiger(cert["precision_decimal"] == 16, "précision du certificat différente")
            marge = Fraction(cert["lambda_min"])-Fraction(float(omega[-1]))**2
            m = e["marge_bande_exacte"]
            _exiger(marge > 0 and marge == Fraction(int(m["numerateur"]), int(m["denominateur"])),
                     "marge de bande Fraction incohérente")
            # Le certificat frais laisse la décision de bande au pilote.
            cert.update(bande_demandee_binary64_couverte=True, omega_max_binary64=float(omega[-1]))
            _certificat(cert, b, omega, 1)
    if e["variante"] == "krylov_controle" and e["statut"] == "termine":
        u = e["controle_uniforme"]
        for cle, valeur in u.items():
            if cle not in ("termes", "delta_colonnes"):
                _nombre(valeur, "contrôle uniforme "+cle)
        termes = u["termes"]
        _exiger([v["degre"] for v in termes] == list(range(1, 8)),
                 "degrés de l'enveloppe incomplets")
        for terme in termes:
            for cle in ("polynome", "queue", "majorant"):
                _nombre(terme[cle], "enveloppe "+cle)
            _grille(terme["majorants_colonnes"], (7,), "enveloppe par colonne")
            _exiger(terme["majorant"] == terme["polynome"]+terme["queue"],
                     "somme polynôme/queue incohérente")
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


def verifier_complet(destination=ARCHIVE):
    """Vérifie l'archive autonome ; ne lit aucun chemin des fichiers exclus."""
    destination = Path(destination)
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
    _exiger(r["schema"] == 1 and r["statut"] == "termine" and r["cpu"] == 8,
             "campagne non terminée")
    _exiger(set(r["sources_sha256"]) == set(SOURCES)
             and r["sources_sha256"]["experience_krylov_contraint.py"] == PILOTE_SHA256,
             "pilote différent du protocole figé")
    for nom, h in r["sources_sha256"].items():
        _exiger(hashlib.sha256(originaux["sources/"+nom]).hexdigest() == _hash(h),
                 "source capturée différente")
    p = r["protocole"]
    strict = dict(maillages=list(MAILLAGES), variantes=list(VARIANTES), max_hz=40.,
        frequences=257, charges=6, blocs=4, directions_retenues=1, profondeur_controle=8,
        precision_decimal=16, repetitions=3, echauffements=1, processus_frais=True,
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
            _candidat(destination, e, nom, omegas[e["n"]])
    _exiger(set(exclus) == exclusions_attendues, "inventaire des exclusions incohérent")
    for x in exclus.values():
        _hash(x["sha256"])
        _exiger(x["dtype"] == "float64" and x["octets"] >= math.prod(x["forme"])*8,
                 "taille/type de fichier exclu incohérent")
    bilan = analyser(r)
    _exiger(bilan == r["analyse"] == _json(originaux["bilan.json"]), "bilan requalifié différent")
    return bilan


def archiver(source, destination=ARCHIVE, journal_parent=None):
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
    return verifier_complet(destination)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", nargs="?", type=Path, default=ARCHIVE)
    parser.add_argument("--creer-depuis", type=Path)
    parser.add_argument("--verifier", type=Path)
    parser.add_argument("--journal-parent", type=Path)
    args = parser.parse_args()
    if args.verifier is not None:
        args.destination = args.verifier
    resultat = (archiver(args.creer_depuis, args.destination, args.journal_parent)
                if args.creer_depuis else verifier_complet(args.destination))
    print(json.dumps(dict(archive=str(args.destination), nombre_essais=resultat["nombre_essais"],
                         statut="archive compacte vérifiée"), ensure_ascii=False))
