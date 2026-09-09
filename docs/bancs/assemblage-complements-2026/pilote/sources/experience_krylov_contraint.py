"""Campagne 40 Hz : Krylov contraint frais, LU corrigée et HCB historique.

Les 60 processus comprennent une chauffe et trois mesures par variante et
maillage. Aucun B, Phi ou certificat préparé n'entre dans un worker. Les
champs, y compris refusés, sont conservés après le chronomètre. Ce pilote
ne lance aucune reprise et refuse tout dossier de sortie déjà existant.
"""
import argparse
from contextlib import redirect_stdout
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import time
import warnings

CI = Path(__file__).resolve().parent
CIBLE = 1e-6
MAILLAGES = (32, 128, 512)
CANDIDATS = ("krylov_champs", "krylov_controle")
VARIANTES = CANDIDATS + ("lu_corrigee", "hcb_standard", "hcb_energie")
RANGS = {32: 185, 128: 761, 512: 768}
# Identités communes aux archives corrigées 0.10.0 et réponses 0.11.0.
# Le manifeste historique est également conservé par la campagne.
DONNEES_FIGEES = {
    32: ("df335d97b6979ea280016d019e4a7fffdc7675cd61c2a5a1fb7aa68d1d4d6b8a",
         "5c19a294a18346328aa2046c2d8d98107010579b01ee0ca81f3de2c8ee9ca50b"),
    128: ("a1549264160cd4210d2e3e60ee31ebc0bd32e214b0925ce03ee1edb1741c740a",
          "2053db09eb87ab41cf368f08303c11a79ed4c8cda8ef2840d5b3abbf7103bb93"),
    512: ("8934fa28cfca143c836b0ec0d78f3ffe16049b1769974c92c9c4f2f27a616905",
          "a0449dfdc2717f19ceec9d359bf848666a91fbdd37fbf31ad9e76771f0f7c05f"),
}
THREADS = {k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                            "MKL_NUM_THREADS", "RAYON_NUM_THREADS")}
FILS_SUPPLEMENTAIRES = {k: "1" for k in ("BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                                       "VECLIB_MAXIMUM_THREADS")}
SOURCES = ("experience_krylov_contraint.py", "krylov_contraint.py",
           "controle_complement.py", "inverse_contrainte_energie.py",
           "racine_masse_blocs.py", "trace_complement_dirigee.py",
           "condensation_energie.py", "confronte_ports_exudyn.py",
           "reference_hcb_exudyn.py", "oracle_champs_ports.py",
           "reference_ports_precision.py", "experience_reduction_native.py")


def maintenant():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for bloc in iter(lambda: f.read(1024*1024), b""):
            h.update(bloc)
    return h.hexdigest()


def jsonable(obj):
    """Conversion hors chronomètre ; non-finis conservés comme diagnostics."""
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if hasattr(obj, "tolist"):
        return jsonable(obj.tolist())
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, float) and not math.isfinite(obj):
        return dict(valeur_non_finie=str(obj))
    return obj


def ecrire(path, contenu, *, remplacer=False):
    path = Path(path)
    texte = json.dumps(jsonable(contenu), ensure_ascii=False, indent=2,
                       allow_nan=False)+"\n"
    if remplacer:
        # Checkpoint appartenant à ce processus ; aucune reprise sur un
        # dossier préexistant. Un lecteur ne voit jamais de JSON tronqué.
        temporaire = path.with_name(path.name+".tmp")
        with temporaire.open("x") as f:
            f.write(texte)
        os.replace(temporaire, path)
    else:
        with path.open("x") as f:
            f.write(texte)


def fichier(path):
    path = Path(path).absolute()
    return dict(path=str(path), sha256=sha(path), octets=path.stat().st_size)


def sources():
    return {nom: sha(CI/nom) for nom in SOURCES}


def environnement(variante):
    """Identités seulement : pas d'inspection de l'implémentation Exudyn."""
    import numpy as np
    import scipy
    import scipy.linalg
    from importlib.metadata import distribution
    from confronte_ports_exudyn import distribution_figee

    def config(module):
        flux = io.StringIO()
        with redirect_stdout(flux):
            module.show_config()
        return flux.getvalue()

    result = dict(python=sys.version, executable=fichier(sys.executable),
                  numpy=np.__version__, scipy=scipy.__version__,
                  plateforme=platform.platform(), machine=platform.machine(),
                  cpu=sorted(os.sched_getaffinity(0)),
                  fils={k: os.environ.get(k) for k in THREADS},
                  fils_supplementaires={k: os.environ.get(k) for k in FILS_SUPPLEMENTAIRES},
                  numpy_config=config(np), scipy_config=config(scipy),
                  distributions={})
    result["cpu_modele"] = next((ligne.split(":", 1)[1].strip()
        for ligne in Path("/proc/cpuinfo").read_text().splitlines()
        if ligne.startswith("model name")), "non disponible")
    for nom in ("numpy", "scipy"):
        dist = distribution(nom)
        result["distributions"][nom] = distribution_figee(nom)
        result["distributions"][nom]["extensions_sha256"] = {
            str(p): sha(dist.locate_file(p)) for p in dist.files or ()
            if str(p).endswith(".so") and Path(dist.locate_file(p)).is_file()}
    if variante in CANDIDATS:
        import vinkulum
        paquet = Path(vinkulum.__file__).parent
        result.update(moteur=vinkulum.__version__,
                      distribution=distribution_figee("vinkulum"),
                      extensions_sha256={p.name: sha(p) for p in paquet.glob("*.so")},
                      sources_python_sha256={str(p.relative_to(paquet)): sha(p)
                                              for p in sorted(paquet.rglob("*.py"))})
        if "roue_sha256" not in result["distribution"]:
            raise ValueError("le candidat requiert une roue locale dont le SHA est accessible")
    elif variante.startswith("hcb"):
        import exudyn
        from exudyn.FEM import FEMinterface, HCBstaticModeSelection  # imports du worker historique
        module = sys.modules[exudyn.SystemContainer.__module__]
        result.update(moteur=exudyn.__version__,
                      distribution=distribution_figee("exudyn"),
                      module_natif=module.__name__, extension_sha256=sha(module.__file__))
        if "roue_sha256" not in result["distribution"]:
            raise ValueError("Exudyn requiert une roue locale dont le SHA est accessible")
    bibliotheques = set()
    for ligne in Path("/proc/self/maps").read_text().splitlines():
        chemin = ligne.split()[-1]
        if chemin.startswith("/") and any(s in chemin.lower() for s in
                ("openblas", "libblas", "liblapack", "libmkl", "libblis", "libgomp", "libiomp")):
            if Path(chemin).is_file():
                bibliotheques.add(chemin)
    result["bibliotheques_numeriques_sha256"] = {p: sha(p) for p in sorted(bibliotheques)}
    try:
        from threadpoolctl import threadpool_info
    except ImportError:
        result["threadpoolctl"] = None
    else:
        result["threadpoolctl"] = threadpool_info()
        if any(p.get("num_threads") != 1 for p in result["threadpoolctl"]):
            raise ValueError("bibliothèque numérique effectivement multithread")
    return result


def donnees(entree, reference, n):
    import numpy as np
    from confronte_ports_exudyn import lire
    if (sha(entree), sha(reference)) != DONNEES_FIGEES[n]:
        raise ValueError("entrée ou oracle différent des témoins historiques figés")
    d, m, metrique, forces, omega = lire(entree)
    masses = m.diagonal()
    coord = m.tocoo()
    if (d.shape[1] != 6*n or m.shape != (6*n, 6*n)
            or metrique.shape != (6, 6) or forces.shape != (6, 6)
            or omega.shape != (257,) or not np.all(np.isfinite(omega))
            or np.any(np.diff(omega) <= 0) or omega[0] != 0
            or not np.array_equal(omega, 2*np.pi*np.linspace(0., 40., 257))
            or np.any((coord.row != coord.col) & (coord.data != 0))
            or np.any(masses <= 0)):
        raise ValueError("console native à masse diagonale, 6 charges et grille 0–40 Hz/257 requise")
    # Le fichier oracle n'est pas chargé avant le chronométrage.
    if not Path(reference).is_file():
        raise ValueError("oracle absent")
    return d, m, metrique, forces, omega


def sauver_champs(dossier, champs, valides):
    import numpy as np
    path = dossier/"champs.npy"
    with path.open("xb") as f:
        np.save(f, champs, allow_pickle=False)
    return dict(fichier=path.name, sha256=sha(path), forme=list(champs.shape),
                indices_frequences_valides=list(valides),
                convention="fréquences × coordonnées physiques × six charges ; NaN si calcul absent")


def controle_oracle(d, m, champs, reference, retours):
    """Diagnostic flottant des majorants ; aucune preuve machine ajoutée."""
    import numpy as np
    result = {nom: [] for nom in ("masse", "deformation")}
    erreurs_absolues = {nom: [] for nom in result}
    normes_reference = {nom: [] for nom in result}
    epsilon = float(np.finfo(float).eps)
    acceptations, refus = [], []
    for j, (x, y, rep) in enumerate(zip(champs, reference, retours, strict=True)):
        delta = x-y
        bon = math.isfinite(rep["marge"]) and rep["marge"] > 0
        raisons = [] if bon else ["marge_bloc_conserve_non_positive"]
        for nom, erreur, norme in (
                ("masse", np.sqrt(np.maximum(np.sum(delta*(m@delta), axis=0), 0.)),
                 np.sqrt(np.maximum(np.sum(y*(m@y), axis=0), 0.))),
                ("deformation", np.linalg.norm(d@delta, axis=0), np.linalg.norm(d@y, axis=0))):
            erreurs_absolues[nom].append(erreur.tolist())
            normes_reference[nom].append(norme.tolist())
            b = rep["bornes"][nom]
            if b["absolues"] is None:
                result[nom].append([False]*6)
            else:
                result[nom].append((erreur <= b["absolues"]+64*epsilon*norme).tolist())
            if not all(result[nom][-1]):
                bon = False
                raisons.append(dict(norme=nom, motif="majorant_non_confirme_par_oracle",
                                    charges=[c for c, ok in enumerate(result[nom][-1]) if not ok]))
            for c, borne in enumerate(b["relatives"]):
                if borne is None or not math.isfinite(borne) or not 0 <= borne <= CIBLE:
                    bon = False
                    raisons.append(dict(norme=nom, charge=c, borne_relative=borne,
                                        motif="majorant_indisponible_ou_hors_tolerance"))
        acceptations.append(bool(bon))
        if raisons:
            refus.append(dict(frequence=j, raisons=raisons))
    return dict(accepte_toutes_frequences=all(acceptations), acceptations=acceptations,
                refus=refus, majorants_confrontes=result,
                erreurs_absolues=erreurs_absolues, normes_reference=normes_reference,
                facteur_marge_arrondi=64, epsilon_binary64=epsilon,
                marge_comparaison="64 epsilon × norme de référence, flottant non certifié")


def worker_candidat(entree, reference, variante, n, dossier):
    import numpy as np
    from scipy.sparse.linalg import LinearOperator, eigsh
    from condensation_energie import CondensationEnergie
    from krylov_contraint import KrylovContraint
    from controle_complement import ControleComplement
    from trace_complement_dirigee import certifier_complement
    from confronte_ports_exudyn import juger, normes

    d, m, metrique, forces, omega = donnees(entree, reference, n)
    result = dict(variante=variante, n=n, max_hz=40., modes_demandes=1,
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
            cert = certifier_complement(d[:, qr.i], qr.r, qr.echelles, mi, b, precision=16)
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
            if variante == "krylov_controle":
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
                                      diagnostic=getattr(exc, "diagnostic", None)))
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
    result.update(n=n, max_hz=40., environnement=identite, sources_sha256=empreintes,
                  sources_apres_sha256=sources(), memoire_comparable=False,
                  resultat_historique=dict(fichier=brut.name, sha256=sha(brut)))
    if capture:
        result["champs"] = sauver_champs(dossier, capture[0], range(len(capture[0])))
    if result["sources_apres_sha256"] != empreintes:
        raise RuntimeError("sources modifiées pendant le témoin")
    ecrire(dossier/"resultat.json", result)


def qualifier(essai):
    """Recalcule la qualification de champ, jamais à partir du seul booléen."""
    if essai["statut"] != "termine":
        return False
    if (essai["nombre_frequences"], essai["charges"]) != (257, 6):
        raise ValueError("couverture incomplète")
    juge = essai["juge"]
    accepte = True
    for nom in ("masse", "deformation", "port"):
        valeurs, op = juge["erreurs"][nom], juge["operateurs"][nom]
        if len(valeurs) != 257 or any(len(ligne) != 6 for ligne in valeurs) or len(op) != 257:
            raise ValueError("juge incomplet")
        nombres = [v for ligne in valeurs for v in ligne]
        if not all(math.isfinite(v) and v >= 0 for v in nombres+op):
            raise ValueError("erreur du juge invalide")
        if max(nombres) != juge["maxima"][nom] or max(op) != juge["maxima_operateurs"][nom]:
            raise ValueError("maximum du juge incohérent")
        accepte &= max(nombres+op) <= CIBLE
    if bool(accepte) != juge["accepte"]:
        raise ValueError("acceptation incohérente")
    for cle in ("preparation_s", "reponses_s", "total_s"):
        if not math.isfinite(essai[cle]) or essai[cle] <= 0:
            raise ValueError("temps invalide")
    if abs(essai["total_s"]-essai["preparation_s"]-essai["reponses_s"]) > 1e-9:
        raise ValueError("total des temps incohérent")
    return bool(accepte)


def qualifier_controle(essai):
    if essai["statut"] != "termine" or "controle_oracle" not in essai:
        return False, False
    diagnostic = essai["controle_oracle"]
    retours = essai["retours_par_frequence"]
    if len(retours) != 257 or len(diagnostic["acceptations"]) != 257:
        raise ValueError("contrôle incomplet")

    def fini(x):
        return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)

    acceptations = []
    for retour in retours:
        bon = fini(retour["marge"]) and retour["marge"] > 0
        if retour["certification_machine"] is not False:
            raise ValueError("certification machine des réponses non autorisée")
        for nom in ("masse", "deformation"):
            valeurs = retour["bornes"][nom]["relatives"]
            if len(valeurs) != 6:
                raise ValueError("majorants par charge incomplets")
            bon &= all(fini(v) and 0 <= v <= CIBLE for v in valeurs)
        acceptations.append(bool(bon))
    couverture = True
    if diagnostic["facteur_marge_arrondi"] != 64 or diagnostic["epsilon_binary64"] != 2.**-52:
        raise ValueError("marge flottante du contrôle différente du protocole")
    for nom in ("masse", "deformation"):
        valeurs = diagnostic["majorants_confrontes"][nom]
        if len(valeurs) != 257 or any(len(ligne) != 6 for ligne in valeurs):
            raise ValueError("confrontation des majorants incomplète")
        if not all(isinstance(v, bool) for ligne in valeurs for v in ligne):
            raise ValueError("confrontation des majorants non booléenne")
        erreurs, normes = diagnostic["erreurs_absolues"][nom], diagnostic["normes_reference"][nom]
        if len(erreurs) != 257 or len(normes) != 257:
            raise ValueError("données de confrontation des majorants incomplètes")
        for j, (err, ref) in enumerate(zip(erreurs, normes, strict=True)):
            if len(err) != 6 or len(ref) != 6 or not all(fini(v) and v >= 0 for v in err+ref):
                raise ValueError("erreurs absolues ou normes de référence invalides")
            absolues = retours[j]["bornes"][nom]["absolues"]
            if absolues is None:
                attendus = [False]*6
            else:
                if len(absolues) != 6 or not all(fini(v) and v >= 0 for v in absolues):
                    raise ValueError("majorant absolu invalide")
                attendus = [e <= b+64*2.**-52*r for e, b, r in zip(err, absolues, ref, strict=True)]
            if attendus != valeurs[j]:
                raise ValueError("majorant confronté incohérent avec les erreurs absolues")
        couverture &= all(v for ligne in valeurs for v in ligne)
        acceptations = [a and all(ligne) for a, ligne in zip(acceptations, valeurs, strict=True)]
    if acceptations != diagnostic["acceptations"] or all(acceptations) != diagnostic["accepte_toutes_frequences"]:
        raise ValueError("qualification du contrôle incohérente")
    return all(acceptations), bool(couverture)


def analyser(rapport):
    if len(rapport["essais"]) != 60:
        raise ValueError("campagne incomplète : 60 essais requis")
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
        if e["variante"] in CANDIDATS and e["statut"] == "termine":
            cert = e["certificat_complement"]
            if cert["precision_decimal"] != 16 or cert["certification_machine"] is not True:
                raise ValueError("certificat différent du protocole")
            marge = e["marge_bande_exacte"]
            if Fraction(int(marge["numerateur"]), int(marge["denominateur"])) <= 0:
                raise ValueError("bande non démontrée")
            if e["selection"]["k"] != 1 or e["selection"]["graine"] != 107+e["n"]:
                raise ValueError("sélection modale différente du protocole")
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
            if v == "krylov_controle":
                controles = [qualifier_controle(e["resultat"]) for e in essais]
                ligne["controle_accepte_tous_essais"] = all(c[0] for c in controles)
                ligne["majorants_couvrent_oracles"] = all(c[1] for c in controles)
                ligne["eligible_controle"] = (ligne["eligible_champs"]
                    and ligne["controle_accepte_tous_essais"] and ligne["majorants_couvrent_oracles"])
                ligne["reponses_certifiees_machine"] = False
            lignes.append(ligne)
    return dict(cible=CIBLE, nombre_essais=60, configurations=lignes,
                portee="consoles natives, réponses matricielles ; HCB officiel avec pont de réponse local, aucune simulation FFRF")


def plan():
    return [(n, v, p) for n in MAILLAGES for p in range(4)
            for v in (VARIANTES if p % 2 == 0 else VARIANTES[::-1])]


def verifier_reference(meta):
    valeurs = [meta[cle][nom] for cle in ("maxima", "maxima_operateurs")
               for nom in ("masse", "deformation", "port")]
    if (meta["precision"] != [70, 90] or max(valeurs) != meta["ecart_relatif_max"]
            or not all(math.isfinite(v) and 0 <= v <= CIBLE/100 for v in valeurs)):
        raise ValueError("référence Decimal70/90 non qualifiée")


def campagne(dossier, entrees, py_candidat, py_lu, py_hcb, historique):
    dossier, entrees, historique = Path(dossier).absolute(), Path(entrees).absolute(), Path(historique).absolute()
    if dossier.exists():
        raise FileExistsError("dossier de sortie déjà présent ; aucune reprise ou réécriture automatique")
    # Validation de provenance avant création des fichiers de campagne.
    archive = json.loads(historique.read_text())
    manifeste_historique = historique.parent/"manifest.json"
    manifeste = json.loads(manifeste_historique.read_text())["fichiers_sha256"]
    if manifeste.get(historique.name) != sha(historique):
        raise ValueError("bilan historique différent de son manifeste")
    rangs = {(e["n"], e["variante"]): e["modes"] for e in archive["configurations"] if e["max_hz"] == 40.}
    for n in MAILLAGES:
        for v in ("hcb_standard", "hcb_energie"):
            if rangs.get((n, v)) != RANGS[n]:
                raise ValueError("rangs HCB historiques inattendus")
    rapport = dict(schema=1, date_utc=maintenant(), statut="en_cours", cpu=8,
        sources_sha256=sources(), historique_hcb=fichier(historique),
        manifeste_historique_hcb=fichier(manifeste_historique), cas={}, essais=[],
        protocole=dict(maillages=list(MAILLAGES), max_hz=40., frequences=257, charges=6,
            variantes=list(VARIANTES), blocs=4, directions_retenues=1, profondeur_controle=8,
            precision_decimal=16, comparaison_bande="Fraction(lambda)>Fraction(omega_max)**2",
            eigsh=dict(k=1, tol=1e-11, ncv=20, graine="107+n_poutres", which="LA"),
            cout="QR, sélection eigsh, B=fl(M Phi), certificat dirigé frais, Krylov, contrôle éventuel, champs et normes M/D",
            exclus="imports, lecture des entrées, empreintes, sauvegardes, oracle et juge",
            repetitions=3, echauffements=1, processus_frais=True, ordre_inverse_par_passage=True,
            hcb="rangs historiques offerts, pont spectral corrigé deux fois dans le worker inchangé",
            reponses_certifiees_machine=False, reprise_automatique=False))
    for n in MAILLAGES:
        nom = f"n{n}-f40"
        cas = json.loads((entrees/(nom+".npz.json")).read_text())
        cas["reference"] = json.loads((entrees/(nom+".ref.npy.json")).read_text())
        verifier_reference(cas["reference"])
        if (cas["n"], cas["max_hz"], cas["nombre_frequences"], cas["charges"]) != (n, 40., 257, 6):
            raise ValueError("métadonnées d'entrée inattendues")
        if sha(entrees/(nom+".npz")) != cas["entree_sha256"] or sha(entrees/(nom+".ref.npy")) != cas["reference"]["fichier_sha256"]:
            raise ValueError("entrée ou oracle altéré")
        if (cas["entree_sha256"], cas["reference"]["fichier_sha256"]) != DONNEES_FIGEES[n]:
            raise ValueError("provenance différente des archives historiques")
        cas["oracle_original"] = str(entrees/(nom+".ref.npy"))
        rapport["cas"][str(n)] = cas
    pythons = {"candidat": str(Path(py_candidat).absolute()), "lu": str(Path(py_lu).absolute()), "hcb": str(Path(py_hcb).absolute())}
    for py in pythons.values():
        if not os.access(py, os.X_OK):
            raise ValueError("Python exécutable absent : "+py)
    rapport["pythons"] = pythons
    dossier.mkdir(parents=True, exist_ok=False)
    (dossier/"sources").mkdir()
    (dossier/"entrees").mkdir()
    (dossier/"essais").mkdir()
    for nom, empreinte in rapport["sources_sha256"].items():
        shutil.copy2(CI/nom, dossier/"sources"/nom)
        if sha(dossier/"sources"/nom) != empreinte:
            raise ValueError("source modifiée pendant sa capture")
    shutil.copy2(historique, dossier/"historique-hcb.json")
    shutil.copy2(manifeste_historique, dossier/"historique-hcb-manifest.json")
    for n in MAILLAGES:
        for suffixe in (".npz", ".npz.json", ".ref.npy.json"):
            nom = f"n{n}-f40{suffixe}"
            shutil.copy2(entrees/nom, dossier/"entrees"/nom)
    ecrire(dossier/"rapport.json", rapport)
    env = dict(os.environ, **THREADS, **FILS_SUPPLEMENTAIRES, PYTHONDONTWRITEBYTECODE="1")
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    for n, v, passage in plan():
        nom = f"n{n}-f40-{v}-passage{passage}"
        sortie = dossier/"essais"/nom
        famille = "candidat" if v in CANDIDATS else "hcb" if v.startswith("hcb") else "lu"
        commande = [pythons[famille], str(dossier/"sources"/Path(__file__).name), "--worker",
                    str(dossier/"entrees"/f"n{n}-f40.npz"), rapport["cas"][str(n)]["oracle_original"],
                    v, str(n), str(sortie)]
        print(nom, flush=True)
        debut = time.perf_counter()
        # Le parent possède le journal ; le worker crée exclusivement son dossier.
        journal = dossier/"essais"/(nom+".log")
        with journal.open("x") as flux:
            processus = subprocess.Popen(commande, env=env, stdout=flux, stderr=subprocess.STDOUT)
            retour = processus.wait()  # Handle terminal vérifié avant l'essai suivant.
        terminal = dict(pid=processus.pid, retour=retour, attente_terminee=True,
                        date_utc=maintenant(), processus_s=time.perf_counter()-debut, commande=commande)
        ecrire(dossier/"essais"/(nom+".terminal.json"), terminal)
        if retour != 0 or not (sortie/"resultat.json").is_file():
            rapport.update(statut="interrompu", essai_interrompu=nom, terminal=terminal)
            ecrire(dossier/"rapport.json", rapport, remplacer=True)
            raise RuntimeError(f"essai {nom} arrêté ({retour}) ; aucune relance automatique")
        resultat = json.loads((sortie/"resultat.json").read_text())
        rapport["essais"].append(dict(nom=nom, passage=passage,
            role="chauffe" if passage == 0 else "mesure", repetition=passage-1,
            resultat=resultat, resultat_sha256=sha(sortie/"resultat.json"),
            journal_sha256=sha(journal), terminal=terminal))
        ecrire(dossier/"rapport.json", rapport, remplacer=True)
    rapport.update(statut="termine", analyse=analyser(rapport), fin_utc=maintenant())
    ecrire(dossier/"rapport.json", rapport, remplacer=True)
    ecrire(dossier/"bilan.json", rapport["analyse"])
    ecrire_manifest(dossier)
    verifier(dossier)
    print("Campagne terminée : "+str(dossier), flush=True)


def ecrire_manifest(dossier):
    ecrire(dossier/"manifest.json", dict(schema=1, date_utc=maintenant(),
        fichiers_sha256={str(p.relative_to(dossier)): sha(p) for p in sorted(dossier.rglob("*"))
                         if p.is_file() and p.name != "manifest.json"}))


def verifier(dossier):
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
    for nom, cle in (("historique-hcb.json", "historique_hcb"),
                     ("historique-hcb-manifest.json", "manifeste_historique_hcb")):
        if sha(dossier/nom) != rapport[cle]["sha256"]:
            raise ValueError("historique HCB incohérent")
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


def archiver(dossier, destination):
    dossier, destination = Path(dossier), Path(destination)
    if destination.exists():
        raise FileExistsError("destination d'archive déjà présente")
    verifier(dossier)
    rapport = json.loads((dossier/"rapport.json").read_text())
    # Données autonomes, y compris oracles et champs hors tolérance ; les
    # installations restent identifiées par leurs empreintes et versions.
    for cas in rapport["cas"].values():
        if sha(cas["oracle_original"]) != cas["reference"]["fichier_sha256"]:
            raise ValueError("oracle d'origine modifié avant archivage")
    shutil.copytree(dossier, destination, ignore=shutil.ignore_patterns("manifest.json"))
    (destination/"oracles").mkdir()
    for n, cas in rapport["cas"].items():
        shutil.copy2(cas["oracle_original"], destination/"oracles"/f"n{n}-f40.ref.npy")
    ecrire_manifest(destination)
    verifier(destination)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    choix = p.add_mutually_exclusive_group(required=True)
    choix.add_argument("--worker", nargs=5, metavar=("ENTREE", "REFERENCE", "VARIANTE", "N", "DOSSIER_NEUF"))
    choix.add_argument("--campagne", nargs=5, metavar=("DOSSIER_NEUF", "ENTREES", "PY_CANDIDAT", "PY_LU", "PY_HCB"))
    choix.add_argument("--archiver", nargs=2, metavar=("DOSSIER", "DESTINATION_NEUVE"))
    choix.add_argument("--verifier", type=Path)
    p.add_argument("--historique", type=Path,
                   default=CI.parent/"docs/bancs/confrontation-ports-exudyn-0.10.0/bilan.json")
    args = p.parse_args()
    if args.worker or args.campagne:
        if 8 not in os.sched_getaffinity(0):
            p.error("CPU 8 indisponible ; protocole inchangé requis")
        os.sched_setaffinity(0, {8})
        os.environ.update(THREADS)
        os.environ.update(FILS_SUPPLEMENTAIRES)
    if args.worker:
        entree, reference, variante, n, dossier = args.worker
        n = int(n)
        if variante not in VARIANTES or n not in MAILLAGES:
            p.error("variante ou maillage non prévu")
        dossier = Path(dossier).absolute()
        dossier.mkdir(parents=True, exist_ok=False)
        if variante in CANDIDATS:
            worker_candidat(entree, reference, variante, n, dossier)
        else:
            worker_temoin(entree, reference, variante, n, dossier)
    elif args.campagne:
        campagne(*args.campagne, args.historique)
    elif args.archiver:
        archiver(*args.archiver)
    else:
        verifier(args.verifier)
        print("Campagne/Archive conforme : "+str(args.verifier))


if __name__ == "__main__":
    main()
