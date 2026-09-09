"""Validation sémantique additionnelle de l'archive du complément spectral.

Aucun solveur n'est chronométré, aucun champ oracle n'est recalculé. Les
hashes, données, décisions et arithmétiques finales des certificats sont
recoupés. Les intervalles H/U/J et les champs bruts ne sont pas archivés :
ce contrôle ne rejoue donc ni la preuve entière ni le juge physique.
"""
import argparse
from copy import deepcopy
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT/"docs/bancs/complement-spectral-2026"
NORMES = {"masse", "deformation", "port"}
SOURCES = {"trace_complement_dirigee.py", "retention_complement.py",
           "experience_complement_spectral.py", "confronte_ports_exudyn.py"}
DIAGNOSTICS = {"sigma_min_bloc_conserve", "defaut_contrainte",
               "residu_kkt", "ecart_fonctionnel"}
FILS = {k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                        "MKL_NUM_THREADS", "RAYON_NUM_THREADS")}


def _exiger(condition, message):
    if not condition:
        raise ValueError(message)


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for bloc in iter(lambda: f.read(1024*1024), b""):
            h.update(bloc)
    return h.hexdigest()


def _hash(h):
    _exiger(isinstance(h, str) and re.fullmatch("[0-9a-f]{64}", h), "hash SHA256 invalide")
    return h


def _nombre(x, nom, positif=False):
    _exiger(type(x) in (int, float) and math.isfinite(x)
             and (x > 0 if positif else x >= 0), nom+" fini non négatif requis")
    return x


def _entier(x, nom, minimum=0):
    _exiger(type(x) is int and x >= minimum, nom+" entier requis")
    return x


def _decimal(x, nom):
    _exiger(isinstance(x, str), nom+" doit être un décimal textuel")
    try:
        d = Decimal(x)
    except Exception as exc:
        raise ValueError(nom+" décimal invalide") from exc
    _exiger(d.is_finite(), nom+" décimal fini requis")
    return d


def _intervalle(r, bas, haut, nom):
    lo, hi = _decimal(r[bas], nom), _decimal(r[haut], nom)
    _exiger(lo <= hi, nom+" intervalle inversé")
    return lo, hi


def _conversion(x, decimal, haut, nom):
    _nombre(x, nom, positif=True)
    exact = Fraction(decimal)
    courant = Fraction(x)
    voisin = float(np.nextafter(float(x), -np.inf if haut else np.inf))
    _exiger(math.isfinite(voisin), nom+" voisin binary64 non fini")
    if haut:
        bon = Fraction(voisin) < exact <= courant
    else:
        bon = courant <= exact < Fraction(voisin)
    _exiger(bon, nom+" conversion binary64 non dirigée")


def _grille(x, forme, nom, positif=False):
    _exiger(isinstance(x, list), nom+" liste requise")
    a = np.asarray(x, dtype=object)
    _exiger(a.shape == forme, nom+" couverture incorrecte")
    for v in a.flat:
        _nombre(v, nom, positif=positif)
    return np.asarray(x, dtype=float)


@lru_cache(maxsize=1)
def _ancres():
    """Lecture unique des preuves historiques ; seuls les en-têtes restent en mémoire."""
    ancien = ROOT/"docs/bancs/confrontation-ports-exudyn-0.10.0"
    manifest = json.loads((ancien/"manifest.json").read_text())
    for nom in ("essais.json.gz", "sources/oracle_champs_ports.py",
                "sources/reference_ports_precision.py", "sources/confronte_ports_exudyn.py"):
        _exiger(_sha(ancien/nom) == manifest["fichiers_sha256"][nom],
                "provenance historique altérée")
    with gzip.open(ancien/"essais.json.gz", "rt") as f:
        historique = json.load(f)
    cas = {c["n"]: {"entree_sha256": c["entree_sha256"], "reference": c["reference"]}
           for c in historique["cas"] if c["max_hz"] == 40}
    _exiger(set(cas) == {32, 128, 512}, "références historiques incomplètes")
    version = json.loads((ROOT/"docs/bancs/version-0.11.0.json").read_text())
    return cas, version


def _csr(z, nom, forme):
    shape = z[nom+"_shape"]
    _exiger(shape.dtype.kind in "iu" and shape.shape == (2,)
             and tuple(shape) == forme, nom+" dimensions incorrectes")
    vals, ids, ptr = z[nom+"_data"], z[nom+"_indices"], z[nom+"_indptr"]
    _exiger(vals.dtype.kind == "f" and vals.dtype.itemsize == 8
             and vals.ndim == 1 and np.all(np.isfinite(vals)), nom+" données binary64 finies requises")
    _exiger(ids.dtype.kind in "iu" and ptr.dtype.kind in "iu"
             and ids.shape == vals.shape and ptr.shape == (forme[0]+1,)
             and ptr[0] == 0 and ptr[-1] == len(vals)
             and np.all(ptr >= 0) and np.all(ptr <= len(vals))
             and np.all(ptr[1:] >= ptr[:-1]) and np.all(ids >= 0)
             and np.all(ids < forme[1]), nom+" stockage CSR invalide")
    for i in range(forme[0]):
        ligne = ids[ptr[i]:ptr[i+1]]
        _exiger(np.all(ligne[1:] > ligne[:-1]), nom+" CSR non canonique")
    return csr_matrix((vals, ids, ptr), shape=forme)


def _tableau(a, forme, nom):
    _exiger(a.shape == forme and a.dtype.kind == "f" and a.dtype.itemsize == 8
             and np.all(np.isfinite(a)), nom+" tableau binary64 fini de dimensions exactes requis")


def _donnees(destination, c):
    n, ni = c["n"], 6*c["n"]-6
    with np.load(destination/f"entrees/n{n}-f40.npz", allow_pickle=False) as z:
        attendus = {f"{nom}_{cle}" for nom in ("d", "m")
                    for cle in ("data", "indices", "indptr", "shape")}
        _exiger(set(z.files) == attendus | {"metrique", "forces", "omega"},
                 "inventaire des données physiques incorrect")
        d, m = _csr(z, "d", (6*n, 6*n)), _csr(z, "m", (6*n, 6*n))
        _tableau(z["metrique"], (6, 6), "métrique")
        # La métrique native provient d'une inversion flottante de
        # compliance. Le contrat QR accepte sa très petite asymétrie ;
        # imposer ici une égalité exacte rejetterait les données mesurées.
        _exiger(np.allclose(z["metrique"], z["metrique"].T, rtol=0.,
                           atol=1e-14*np.max(abs(z["metrique"]))), "métrique non symétrique")
        try:
            np.linalg.cholesky(z["metrique"])
        except np.linalg.LinAlgError as exc:
            raise ValueError("métrique non positive") from exc
        _tableau(z["forces"], (6, 6), "forces")
        _exiger(np.array_equal(z["forces"], np.diag([.1]*3+[.01]*3)), "charges natives modifiées")
        _tableau(z["omega"], (257,), "omega")
        omega = z["omega"].copy()
        _exiger(np.array_equal(omega, 2*np.pi*np.linspace(0., 40., 257)),
                 "fréquences NPZ différentes de la grille annoncée")
        _exiger(m.nnz == 6*n and np.all(m.diagonal() > 0)
                 and (m-m.T).nnz == 0, "masse native diagonale positive requise")
    with np.load(destination/f"directions/n{n}-directions.npz", allow_pickle=False) as z:
        _exiger(set(z.files) == {"phi", "b", "echelles", "r_data", "r_indices",
                                "r_indptr", "r_shape"}, "inventaire des directions incorrect")
        r = _csr(z, "r", (ni, ni))
        _exiger(np.all(r.diagonal() != 0), "facteur R singulier")
        for i in range(ni):
            _exiger(np.all(r.indices[r.indptr[i]:r.indptr[i+1]] >= i),
                     "facteur R non triangulaire supérieur")
        _tableau(z["echelles"], (ni,), "échelles")
        _exiger(np.all(z["echelles"] > 0), "échelles non positives")
        b, phi = z["b"].copy(), z["phi"].copy()
        _exiger(b.ndim == 2 and b.shape[0] == ni and b.shape[1] >= 2,
                 "B de dimensions incorrectes")
        _tableau(b, b.shape, "B")
        _tableau(phi, b.shape, "Phi")
        # La réponse n'emploie que le rang1. Cette vérification exacte du
        # couplage de coordonnées ne suppose pas B=M Phi.
        couplage = sum((Fraction(float(x))*Fraction(float(y))
                        for x, y in zip(b[:, 0], phi[:, 0])), Fraction(0))
        _exiger(couplage != 0, "B.T Phi singulier pour les coordonnées retenues")
    return b, omega


def _certificat(v, b, omega, rang):
    ni = b.shape[0]
    _exiger(v["rang_contraintes"] == rang and type(v["rang_contraintes"]) is int
             and v["dimension_interieure"] == ni and type(v["dimension_interieure"]) is int
             and v["dimension_complement"] == ni-rang and type(v["dimension_complement"]) is int,
             "dimensions du certificat incohérentes")
    _exiger(v["certification_machine"] is True
             and v["bande_demandee_binary64_couverte"] is True, "certificat ou couverture absent")
    _exiger(v["contraintes_cibles"] ==
             "noyau de B.T ; valeurs binary64 de B interprétées exactement",
             "convention des contraintes différente")
    _exiger(v["portee"] == "coercivité sur le complément seulement ; ni Schur global ni réponse certifiés",
             "portée du certificat modifiée")
    h = hashlib.sha256(np.asarray(b[:, :rang], dtype="<f8").tobytes()).hexdigest()
    _exiger(_hash(v["contraintes_sha256"]) == h, "hash exact de B différent des données NPZ")
    _exiger(v["coefficients_rectangulaires"] == ni*rang+rang*rang
             and type(v["coefficients_rectangulaires"]) is int, "compteur rectangulaire incohérent")
    precision = _entier(v["precision_decimal"], "précision", 16)
    _exiger(precision <= 4096, "précision déclarée excessive pour ce format")
    base = v["certificat_total"]
    _exiger(base["precision_decimal"] == precision
             and base["certification_machine"] is True
             and base["minoration_D_original_etablie"] is True, "certificat énergétique initial absent")
    _exiger(base["modele"] == "D_I.T D_I et M_II : valeurs binary64 d'entrée interprétées exactement",
             "modèle spectral différent")
    total_lo, total_hi = _intervalle(base, "trace_inferieure_decimal",
                                    "trace_superieure_decimal", "trace totale")
    _exiger(total_lo > 0, "trace totale non positive")
    eta = _decimal(base["eta_superieur_decimal"], "eta")
    _exiger(0 <= eta < 1, "perturbation énergétique hors domaine")
    chi_lo, chi_hi = _intervalle(v, "correction_inferieure_decimal",
                                "correction_superieure_decimal", "correction")
    _exiger(chi_hi > 0, "correction supérieure non positive")
    lo, hi = _intervalle(v, "trace_complement_inferieure_decimal",
                        "trace_complement_superieure_decimal", "trace complémentaire")
    _exiger(hi > 0, "trace complémentaire supérieure non positive")
    bas, haut = Context(prec=precision, rounding=ROUND_FLOOR), Context(prec=precision, rounding=ROUND_CEILING)
    _exiger(lo == bas.subtract(total_lo, chi_hi) and hi == haut.subtract(total_hi, chi_lo),
             "soustraction dirigée des traces incohérente")
    pivots = v["pivots_h_inferieurs_decimal"]
    _exiger(isinstance(pivots, list) and len(pivots) == rang
             and all(_decimal(p, "pivot H") > 0 for p in pivots), "pivots H non strictement positifs")
    for objet, trace, nom in ((base, total_hi, "lambda totale"), (v, hi, "lambda complémentaire")):
        lam = _decimal(objet["lambda_inferieur_decimal"], nom)
        _exiger(lam > 0 and lam == bas.divide(bas.subtract(Decimal(1), eta), trace),
                 nom+" quotient dirigé incohérent")
        _conversion(objet["lambda_min"], lam, False, nom)
        _nombre(objet["preparation_s"], "préparation certificat", positif=True)
        _entier(objet["operations_decimal"], "opérations Decimal", 1)
    _conversion(base["trace_superieure"], total_hi, True, "trace supérieure")
    if eta == 0:
        _exiger(type(base["eta_superieur"]) in (int, float) and base["eta_superieur"] == 0,
                 "conversion eta incohérente")
    else:
        _conversion(base["eta_superieur"], eta, True, "eta supérieur")
    _exiger(v["operations_decimal"] >= base["operations_decimal"],
             "compteur Decimal total inférieur au certificat initial")
    for cle in ("coefficients_gram_ecart", "produits_gram", "coefficients_inverse"):
        _entier(base[cle], cle)
    _entier(v["operations_inverse_float"], "opérations inverse", 1)
    _nombre(v["omega_max_binary64"], "omega du certificat", positif=True)
    _exiger(v["omega_max_binary64"] == float(omega[-1]), "bande du certificat différente des omega NPZ")
    _exiger(Fraction(v["lambda_min"]) > max(Fraction(float(w))**2 for w in omega),
             "bande NPZ non couverte en Fraction")


def _juge(juge):
    _exiger(set(juge) == {"maxima", "maxima_operateurs", "erreurs", "operateurs", "accepte"},
             "structure du juge incomplète")
    maxima = []
    for cle in ("maxima", "maxima_operateurs", "erreurs", "operateurs"):
        _exiger(set(juge[cle]) == NORMES, "normes du juge incomplètes")
    for nom in NORMES:
        e = _grille(juge["erreurs"][nom], (257, 6), "erreurs "+nom)
        o = _grille(juge["operateurs"][nom], (257,), "opérateurs "+nom)
        for valeurs, cle in ((e, "maxima"), (o, "maxima_operateurs")):
            maxi = _nombre(juge[cle][nom], "maximum "+nom)
            _exiger(maxi == float(np.max(valeurs)), "maximum du juge incohérent")
            maxima.append(maxi)
    _exiger(type(juge["accepte"]) is bool
             and juge["accepte"] == (max(maxima) <= 1e-6), "qualification du juge incohérente")


def verifier_complet(destination=ARCHIVE):
    """Vérifie l'archive et ses liens historiques, sans rejouer une expérience."""
    destination = Path(destination)
    try:
        manifest = json.loads((destination/"manifest.json").read_text())
        _exiger(manifest["version_support"] == "0.11.0"
                 and manifest["prototype"] == "complément spectral, hors API de production",
                 "version ou portée du manifeste incorrecte")
        hashes = manifest["fichiers_sha256"]
        necessaires = {"rapport.json.gz"} | {"sources/"+n for n in SOURCES}
        necessaires |= {f"entrees/n{n}-f40.npz" for n in (32, 128, 512)}
        necessaires |= {f"directions/n{n}-directions.npz" for n in (32, 128, 512)}
        autorises = necessaires | {"README.md"}
        fichiers = [p for p in destination.rglob("*") if p.is_file() and p != destination/"manifest.json"]
        _exiger(all(not p.is_symlink() for p in destination.rglob("*")), "lien symbolique dans l'archive")
        presents = {str(p.relative_to(destination)) for p in fichiers}
        _exiger(necessaires <= presents <= autorises and set(hashes) == presents,
                 "inventaire d'archive incomplet ou inattendu")
        for nom, h in hashes.items():
            _exiger(_sha(destination/nom) == _hash(h), "archive altérée")
        rapport = json.loads(gzip.decompress((destination/"rapport.json.gz").read_bytes()))
        _exiger(len(rapport["cas"]) == 3 and [c["n"] for c in rapport["cas"]] == [32, 128, 512],
                 "couverture des maillages incomplète")
        _exiger(set(rapport["sources_sha256"]) == SOURCES, "sources du rapport incomplètes")
        for nom, h in rapport["sources_sha256"].items():
            _exiger(hashes["sources/"+nom] == h, "provenance source différente")
        historique, version = _ancres()
        livraison = version["controle_roue"]
        wheel = livraison["concordance_mesure_livraison"]["roue_mesuree"]["sha256"]
        _exiger(rapport["vinkulum"] == rapport["distribution"]["version"] == "0.11.0"
                 and rapport["distribution"]["roue_sha256"] == wheel
                 and rapport["extension_sha256"] == version["roue"]["extension_sha256"],
                 "installation différente de la version mesurée")
        _hash(rapport["distribution"]["record_sha256"])
        sources_attendues = {k.removeprefix("vinkulum/"): v
                            for k, v in livraison["sources_python_sha256"].items()}
        _exiger(rapport["sources_python_sha256"] == sources_attendues,
                 "sources Python différentes de la roue qualifiée")
        _exiger(rapport["cpu"] == [8] and rapport["fils"] == FILS,
                 "affinité ou budget de fils incorrect")
        for nom in ("python", "numpy", "scipy"):
            _exiger(rapport[nom] == version["environnement"][nom], "version numérique différente")
        _exiger(rapport["protocole"] ==
                 "une exécution, durées indicatives ; aucune comparaison de vitesse",
                 "protocole transformé en classement de performances")
        for c in rapport["cas"]:
            _entier(c["n"], "maillage", 1)
            n, ni = c["n"], 6*c["n"]-6
            _exiger(c["frequences"] == 257 and type(c["frequences"]) is int
                     and c["charges"] == 6 and type(c["charges"]) is int
                     and c["bande_hz"] == [0, 40], "couverture des charges ou fréquences incorrecte")
            _exiger(c["taille_conservee"] == 7 and type(c["taille_conservee"]) is int
                     and c["taille_kkt"] == ni+1 and type(c["taille_kkt"]) is int,
                     "dimensions du témoin incohérentes")
            _exiger(c["certification_champ_machine"] is False, "champ improprement déclaré certifié")
            for dossier, fichier, cle in (
                    ("entrees", f"n{n}-f40.npz", "entree_sha256"),
                    ("directions", f"n{n}-directions.npz", "directions_sha256")):
                _exiger(hashes[f"{dossier}/{fichier}"] == _hash(c[cle]), "provenance donnée différente")
            _exiger(c["entree_sha256"] == historique[n]["entree_sha256"],
                     "entrée différente de l'oracle historique")
            ref = c["reference"]
            _exiger(ref["precision"] == [70, 90]
                     and _hash(ref["fichier_sha256"]) == _hash(c["reference_sha256"]),
                     "référence 70/90 ou hash incohérent")
            maximum = []
            for cle in ("maxima", "maxima_operateurs"):
                _exiger(set(ref[cle]) == NORMES, "qualification de référence incomplète")
                maximum += [_nombre(ref[cle][nom], "erreur de référence") for nom in NORMES]
            _exiger(_nombre(ref["ecart_relatif_max"], "écart référence") == max(maximum)
                     and max(maximum) <= 1e-8, "référence non qualifiée")
            _nombre(ref["temps_hors_mesures_s"], "temps référence", positif=True)
            _exiger(ref == historique[n]["reference"], "référence différente de l'archive historique")
            _exiger(set(c["diagnostics"]) == DIAGNOSTICS, "diagnostics incomplets")
            for nom in DIAGNOSTICS:
                _grille(c["diagnostics"][nom], (257,), nom,
                        positif=(nom == "sigma_min_bloc_conserve"))
            for nom in ("preparation_temoin_indicative_s", "reponses_temoin_indicative_s"):
                _nombre(c[nom], nom, positif=True)
            _juge(c["juge"])
            b, omega = _donnees(destination, c)
            _exiger(len(c["certificats"]) == 2, "couverture des certificats incomplète")
            for rang, cert in enumerate(c["certificats"], 1):
                _certificat(cert, b, omega, rang)
        return rapport
    except (KeyError, TypeError, IndexError, OverflowError) as exc:
        raise ValueError("structure de l'archive invalide") from exc


class ArchiveComplementSpectral(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads(gzip.decompress((ARCHIVE/"rapport.json.gz").read_bytes()))

    def alterer(self, mutation, motif, donnees=None):
        """Mutations cohérentes avec un manifeste recalculé : tester le sens."""
        with TemporaryDirectory() as tmp:
            dest = Path(tmp)/"archive"
            shutil.copytree(ARCHIVE, dest)
            r = deepcopy(self.original)
            mutation(r)
            if donnees:
                donnees(dest, r)
            (dest/"rapport.json.gz").write_bytes(gzip.compress(
                json.dumps(r, ensure_ascii=False).encode(), mtime=0))
            m = json.loads((dest/"manifest.json").read_text())
            m["fichiers_sha256"] = {str(p.relative_to(dest)): _sha(p)
                for p in dest.rglob("*") if p.is_file() and p.name != "manifest.json"}
            (dest/"manifest.json").write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError, motif):
                verifier_complet(dest)

    def test_archive_reelle_complete(self):
        r = verifier_complet(ARCHIVE)
        self.assertEqual(sum(len(c["certificats"]) for c in r["cas"]), 6)
        self.assertTrue(all(c["juge"]["accepte"] for c in r["cas"]))

    def test_erreurs_nan_inf_negatives_et_booleennes(self):
        for v in (float("nan"), float("inf"), -1., True):
            with self.subTest(v=v):
                self.alterer(lambda r: r["cas"][0]["juge"]["erreurs"]["masse"][0].__setitem__(0, v),
                             "fini non négatif")

    def test_maxima_et_qualification(self):
        self.alterer(lambda r: r["cas"][0]["juge"]["maxima"].__setitem__("masse", .1), "maximum")
        def fausse_acceptation(r):
            j = r["cas"][0]["juge"]
            j["operateurs"]["port"][0] = j["maxima_operateurs"]["port"] = .1
        self.alterer(fausse_acceptation, "qualification")
        self.alterer(lambda r: r["cas"][0]["juge"].__setitem__("accepte", 1), "qualification")

    def test_diagnostics_couverture_et_signe(self):
        self.alterer(lambda r: r["cas"][0]["diagnostics"]["residu_kkt"].pop(), "couverture")
        for v in (float("nan"), -1.):
            self.alterer(lambda r: r["cas"][0]["diagnostics"]["defaut_contrainte"].__setitem__(0, v),
                         "fini non négatif")
        self.alterer(lambda r: r["cas"][0]["diagnostics"]["sigma_min_bloc_conserve"].__setitem__(0, 0),
                     "fini non négatif")

    def test_reference_precisions_erreurs_et_hash(self):
        self.alterer(lambda r: r["cas"][0]["reference"].__setitem__("precision", [50, 70]), "70/90")
        self.alterer(lambda r: r["cas"][0]["reference"]["maxima"].__setitem__("masse", 1e-4),
                     "non qualifiée")
        def hash_coherent_mais_autre(r):
            r["cas"][0]["reference"]["fichier_sha256"] = r["cas"][0]["reference_sha256"] = "0"*64
        self.alterer(hash_coherent_mais_autre, "historique")

    def test_dimensions_rang_et_statuts(self):
        self.alterer(lambda r: r["cas"][0].__setitem__("taille_conservee", 6), "dimensions")
        self.alterer(lambda r: r["cas"][0]["certificats"][0].__setitem__("dimension_complement", 186),
                     "dimensions")
        self.alterer(lambda r: r["cas"][0]["certificats"][0].__setitem__("certification_machine", 1),
                     "certificat")
        self.alterer(lambda r: r["cas"][0].__setitem__("certification_champ_machine", True), "champ")

    def test_pivots_eta_et_intervalles(self):
        self.alterer(lambda r: r["cas"][0]["certificats"][0]["pivots_h_inferieurs_decimal"].__setitem__(0, "0"),
                     "pivots")
        self.alterer(lambda r: r["cas"][0]["certificats"][0]["certificat_total"].__setitem__("eta_superieur_decimal", "1"),
                     "perturbation")
        self.alterer(lambda r: r["cas"][0]["certificats"][0].__setitem__("trace_complement_inferieure_decimal", "1"),
                     "intervalle")
        self.alterer(lambda r: r["cas"][0]["certificats"][0].__setitem__("trace_complement_superieure_decimal", "0.1"),
                     "soustraction")

    def test_quotient_et_conversion_diriges(self):
        self.alterer(lambda r: r["cas"][0]["certificats"][0].__setitem__("lambda_inferieur_decimal", "200000"),
                     "quotient")
        def voisin_trop_haut(r):
            v = r["cas"][0]["certificats"][0]
            v["lambda_min"] = float(np.nextafter(v["lambda_min"], np.inf))
        self.alterer(voisin_trop_haut, "conversion binary64")
        self.alterer(lambda r: r["cas"][0]["certificats"][0].__setitem__("lambda_min", True),
                     "fini non négatif")

    def test_bande_est_celle_du_npz(self):
        self.alterer(lambda r: r["cas"][0]["certificats"][0].__setitem__("omega_max_binary64", 1.),
                     "omega NPZ")
        self.alterer(lambda r: r["cas"][0].__setitem__("bande_hz", [0, 20]), "couverture")

    def test_npz_b_hash_dimensions_et_valeurs(self):
        def modifier(dest, rapport, mode):
            path = dest/"directions/n32-directions.npz"
            with np.load(path, allow_pickle=False) as f:
                z = {k: f[k].copy() for k in f.files}
            if mode == "hash":
                z["b"][0, 0] += .125
            elif mode == "nan":
                z["phi"][0, 0] = np.nan
            elif mode == "forme":
                z["phi"] = z["phi"][:-1]
            elif mode == "coordonnees":
                z["phi"][:, 0] = 0
            np.savez_compressed(path, **z)
            rapport["cas"][0]["directions_sha256"] = _sha(path)
        for mode, motif in (("hash", "hash exact de B"), ("nan", "Phi tableau"),
                            ("forme", "Phi tableau"), ("coordonnees", "B.T Phi")):
            self.alterer(lambda r: None, motif,
                         lambda dest, r: modifier(dest, r, mode))

    def test_installation_sources_et_temps(self):
        self.alterer(lambda r: r["distribution"].__setitem__("roue_sha256", "0"*64), "installation")
        self.alterer(lambda r: r["sources_python_sha256"].__setitem__("_ports/certificat_spectral.py", "0"*64),
                     "sources Python")
        self.alterer(lambda r: r["cas"][0].__setitem__("reponses_temoin_indicative_s", -1.), "fini non négatif")

    def test_source_modifiee_meme_avec_manifest_recalcule(self):
        def modifier(dest, r):
            p = dest/"sources/retention_complement.py"
            p.write_text(p.read_text()+"\n# modification après mesure\n")
        self.alterer(lambda r: None, "provenance source", modifier)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verifier", type=Path)
    options, suite = parser.parse_known_args()
    if options.verifier:
        verifier_complet(options.verifier)
        print("Archive sémantique conforme : 6 certificats, 3 bandes et provenance historique")
    else:
        unittest.main(argv=[__file__]+suite)
