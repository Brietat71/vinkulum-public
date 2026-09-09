"""Confrontation matricielle : ports Vinkulum, HCB officiel Exudyn et LU.

Les modèles sont les mêmes D/M natifs figés. Le pont de réponse HCB est
notre enveloppe des bases officielles, pas une simulation FFRF Exudyn.
Chaque processus frais calcule les six champs de réponse à 257 fréquences.
Les références Decimal et le juge physique sont hors chronométrage.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
THREADS = {k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                            "MKL_NUM_THREADS", "RAYON_NUM_THREADS")}
VARIANTES = ("vinkulum_api", "vinkulum_reduit", "lu", "lu_corrigee",
            "hcb_standard", "hcb_energie")
CIBLE = 1e-6


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def distribution_figee(nom):
    """Empreinte du manifeste installé ; aucune consultation du code concurrent."""
    from importlib.metadata import distribution
    dist = distribution(nom)
    record = dist.read_text("RECORD")
    if record is None:
        raise ValueError("distribution sans manifeste RECORD : "+nom)
    resultat = dict(version=dist.version,
                    record_sha256=hashlib.sha256(record.encode()).hexdigest())
    provenance = dist.read_text("direct_url.json")
    if provenance:
        from urllib.parse import unquote, urlparse
        url = urlparse(json.loads(provenance)["url"])
        roue = Path(unquote(url.path))
        if url.scheme == "file" and not url.netloc and roue.suffix == ".whl" and roue.is_file():
            resultat["roue_sha256"] = sha(roue)
    return resultat


def ecrire(path, valeur):
    Path(path).write_text(json.dumps(valeur, ensure_ascii=False, indent=2,
                                   allow_nan=False)+"\n")


def sparse_pack(nom, matrice):
    a = matrice.tocsr()
    return {nom+"_data": a.data, nom+"_indices": a.indices,
            nom+"_indptr": a.indptr, nom+"_shape": a.shape}


def lire(path):
    import numpy as np
    from scipy.sparse import csr_matrix
    with np.load(path, allow_pickle=False) as z:
        matrices = [csr_matrix((z[n+"_data"], z[n+"_indices"], z[n+"_indptr"]),
                               shape=tuple(z[n+"_shape"])) for n in ("d", "m")]
        return (*matrices, *(z[n].copy() for n in ("metrique", "forces", "omega")))


def entrees(n, max_hz, sortie, nombre=257):
    import numpy as np
    import vinkulum
    from scipy.sparse import csc_matrix
    from experience_reduction_native import parametres, modele
    p = parametres(n)
    noyau, libres, ports, metrique = modele(p)
    export = noyau.facteurs_materiels_poutres()
    def convertir(nom):
        nl, nc, ptr, ids, data = export[nom]
        return csc_matrix((data, ids, ptr), shape=(nl, nc))
    d = convertir("d")[:, libres].tocsr()
    m = convertir("masse")[libres][:, libres].tocsr()
    d.sort_indices(); m.sort_indices()
    forces = np.diag([.1, .1, .1, .01, .01, .01])
    np.savez_compressed(sortie, **sparse_pack("d", d), **sparse_pack("m", m),
                        metrique=metrique, forces=forces,
                        omega=2*np.pi*np.linspace(0., max_hz, nombre))
    return dict(n=n, max_hz=max_hz, nombre_frequences=nombre, charges=6,
                modele=p, vinkulum=vinkulum.__version__, entree_sha256=sha(sortie),
                extension_sha256=sha(next(Path(vinkulum.__file__).parent.glob("*.so"))))


def references(entree, sortie):
    import numpy as np
    from oracle_champs_ports import OracleChamps
    d, m, metrique, forces, omega = lire(entree)
    debut = time.perf_counter()
    a = OracleChamps(d, m, p=6, dps=70)
    b = OracleChamps(d, m, p=6, dps=90)
    champs, champs70 = [], []
    for w in omega:
        x, y = a.reponse(w, forces), b.reponse(w, forces)
        champs70.append(x)
        champs.append(y)
    juge = juger(d, m, metrique, np.asarray(champs70), np.asarray(champs))
    ecart = max(*juge["maxima"].values(), *juge["maxima_operateurs"].values())
    if ecart > CIBLE/100:
        raise ArithmeticError("les deux précisions de référence ne concordent pas au seuil requis")
    np.save(sortie, np.asarray(champs), allow_pickle=False)
    return dict(precision=[70, 90], ecart_relatif_max=ecart,
                controle="colonnes et opérateurs masse/déformation/port, après conversion binary64",
                maxima=juge["maxima"], maxima_operateurs=juge["maxima_operateurs"],
                fichier_sha256=sha(sortie), temps_hors_mesures_s=time.perf_counter()-debut)


def normes(d, m, x):
    import numpy as np
    return (np.sqrt(np.maximum(np.sum(x*(m@x), axis=0), 0.)),
            np.linalg.norm(d@x, axis=0))


def juger(d, m, metrique, candidat, reference):
    """Erreurs physiques par charge et sur toutes leurs combinaisons réelles.

    Les six colonnes de référence sont normalisées avant QR pour éviter
    qu'une différence d'amplitude entre axial et flexion soit prise pour
    une perte de rang. La même transformation s'applique à l'erreur :
    ||B_erreur R^-1|| reste la norme opérateur relative exacte. Le rang
    numérique insuffisant est refusé ; aucune direction n'est supprimée.
    """
    import numpy as np
    from scipy.linalg import cholesky, qr, solve_triangular, svdvals
    from scipy.sparse import issparse
    if any(np.iscomplexobj(a) for a in (d, m, metrique, candidat, reference)):
        raise ValueError("le juge physique requiert des données réelles")
    candidat, reference = np.asarray(candidat), np.asarray(reference)
    metrique = np.asarray(metrique, dtype=float)
    if (candidat.ndim != 3 or candidat.shape != reference.shape
            or candidat.shape[0] == 0 or candidat.shape[1] < 6 or candidat.shape[2] != 6):
        raise ValueError("champs de mêmes dimensions fréquences×coordonnées×6 requis")
    total = candidat.shape[1]
    if (len(d.shape) != 2 or d.shape[1] != total or m.shape != (total, total)
            or metrique.shape != (6, 6)):
        raise ValueError("dimensions incompatibles de D, M ou de la métrique de port")
    for nom, a in (("D", d), ("M", m), ("métrique", metrique),
                   ("candidat", candidat), ("référence", reference)):
        valeurs = a.data if issparse(a) else np.asarray(a)
        if not np.all(np.isfinite(valeurs)):
            raise ValueError(nom+" contient des valeurs non finies")
    masses = np.asarray(m.diagonal()).ravel()
    if issparse(m):
        triplets = m.tocoo()
        hors_diagonale = np.any((triplets.row != triplets.col) & (triplets.data != 0))
    else:
        hors_diagonale = np.any(np.asarray(m) != np.diag(masses))
    if hors_diagonale or np.any(masses <= 0):
        raise ValueError("le juge actuel requiert une masse diagonale strictement positive")
    eps = np.finfo(float).eps
    if np.max(np.abs(metrique-metrique.T)) > 64*eps*np.max(np.abs(metrique)):
        raise ValueError("métrique de port non symétrique")
    # Même convention que la normalisation des ports : Cholesky du
    # triangle inférieur ; une dissymétrie significative est refusée.
    l = cholesky(metrique, lower=True)
    racine_masse = np.sqrt(masses)[:, None]
    erreurs = {n: [] for n in ("masse", "deformation", "port")}
    operateurs = {n: [] for n in erreurs}
    for frequence, (x, y) in enumerate(zip(candidat, reference, strict=True)):
        delta = x-y
        couples = (("masse", racine_masse*delta, racine_masse*y),
                   ("deformation", d@delta, d@y),
                   ("port", l.T@delta[-6:], l.T@y[-6:]))
        for nom, be, br in couples:
            if not np.all(np.isfinite(be)) or not np.all(np.isfinite(br)):
                raise ArithmeticError(f"champ pondéré non fini ({nom}, fréquence {frequence})")
            numerateur = np.linalg.norm(be, axis=0)
            denominateur = np.linalg.norm(br, axis=0)
            if (not np.all(np.isfinite(numerateur)) or not np.all(np.isfinite(denominateur))
                    or np.any(denominateur <= 0)):
                raise ArithmeticError(f"référence de norme nulle ou non représentable ({nom}, fréquence {frequence})")
            colonnes = numerateur/denominateur
            reference_normalisee = br/denominateur[None, :]
            erreur_normalisee = be/denominateur[None, :]
            if br.shape[0] < 6:
                raise ArithmeticError(f"référence de rang inférieur à six ({nom}, fréquence {frequence})")
            _, r = qr(reference_normalisee, mode="economic", check_finite=False)
            valeurs = svdvals(r, check_finite=False)
            seuil_rang = eps*max(reference_normalisee.shape)*valeurs[0]
            if valeurs[-1] <= seuil_rang:
                raise ArithmeticError(f"référence numériquement de rang inférieur à six ({nom}, fréquence {frequence})")
            # Résolution à droite ; aucune inverse explicite ni Gram formé.
            relatif = solve_triangular(r.T, erreur_normalisee.T, lower=True,
                                       check_finite=False).T
            operateur = float(svdvals(relatif, check_finite=False)[0])
            if not np.all(np.isfinite(colonnes)) or not np.isfinite(operateur):
                raise ArithmeticError(f"erreur relative non finie ({nom}, fréquence {frequence})")
            erreurs[nom].append(colonnes.tolist())
            operateurs[nom].append(operateur)
    maxima = {n: float(np.max(v)) for n, v in erreurs.items()}
    maxima_operateurs = {n: float(np.max(v)) for n, v in operateurs.items()}
    return dict(erreurs=erreurs, maxima=maxima, operateurs=operateurs,
                maxima_operateurs=maxima_operateurs,
                accepte=all(v <= CIBLE for v in (*maxima.values(), *maxima_operateurs.values())))


def worker(entree, reference, variante, modes, sortie):
    import resource
    import warnings
    import numpy as np
    import scipy
    from scipy.linalg import solve
    from scipy.sparse import diags
    from scipy.sparse.linalg import splu
    d, m, metrique, forces, omega = lire(entree)
    n, ports = d.shape[1], np.arange(d.shape[1]-6, d.shape[1])
    result = dict(variante=variante, modes_demandes=modes, entree_sha256=sha(entree),
        reference_sha256=sha(reference), python=sys.version, numpy=np.__version__,
        scipy=scipy.__version__, date_utc=datetime.now(timezone.utc).isoformat(),
        cpu=sorted(os.sched_getaffinity(0)), fils={k: os.environ.get(k) for k in THREADS},
        nombre_frequences=len(omega), charges=6, cible=CIBLE, statut="en_cours")
    # Imports propres à la méthode exclus, mais leurs bibliothèques restent
    # dans le RSS. L'assemblage D.T D et toutes les préparations sont comptés.
    if variante.startswith("vinkulum"):
        import vinkulum
        from vinkulum.reduction_ports import ReductionMaterielle
        from vinkulum._ports import certificat_spectral, inverse_selectionnee  # imports différés exclus
        result.update(moteur=vinkulum.__version__, extension_sha256=sha(
            next(Path(vinkulum.__file__).parent.glob("*.so"))))
        paquet = Path(vinkulum.__file__).parent
        result.update(distribution=distribution_figee("vinkulum"),
                      sources_python_sha256={str(p.relative_to(paquet)): sha(p)
                                              for p in sorted(paquet.rglob("*.py"))})
    elif variante.startswith("hcb"):
        import exudyn
        from exudyn.FEM import FEMinterface, HCBstaticModeSelection  # import différé exclu
        from reference_hcb_exudyn import prepare, reponse
        module = sys.modules[exudyn.SystemContainer.__module__]
        result.update(moteur=exudyn.__version__, extension_sha256=sha(module.__file__),
                      module_natif=module.__name__,
                      resolution_hcb="spectrale_corrigee", corrections_hcb=2,
                      distribution=distribution_figee("exudyn"))
        np.random.seed(42123)
    champs = np.empty((len(omega), n, 6))
    nm, nd = np.empty((len(omega), 6)), np.empty((len(omega), 6))
    bornes = {nom: [] for nom in ("masse", "deformation")}
    debut = time.perf_counter()
    phase = "preparation"
    try:
        with warnings.catch_warnings(record=True) as avertissements:
            warnings.simplefilter("always")
            if variante.startswith("vinkulum"):
                r = ReductionMaterielle(d, m, np.arange(n-6), ports, metrique, omega[-1],
                     tolerance=1e-10, max_blocs=12, max_directions=128)
                f_normal = r.qr.w.T@forces
                result.update(directions=r.taille_interieure, certificat_spectral=r.certificat_spectral,
                              statut_schur=r.statut, borne_schur=r.borne_uniforme)
            elif variante.startswith("hcb"):
                r = prepare(d, m, modes, projection=variante.split("_")[1],
                            diagnostics_complets=False)
                result.update(directions=modes, diagnostic_hcb=r["diagnostics"])
            else:
                k = (d.T@d).tocsc()
                eq = 1/np.sqrt(k.diagonal())
                s = diags(eq)
                ke, me = (s@k@s).tocsc(), (s@m@s).tocsc()
                f = np.zeros((n, 6)); f[-6:] = forces
                fe = eq[:, None]*f
                result["directions"] = n-6
            prepare = time.perf_counter()
            phase = "reponses"
            for j, w in enumerate(omega):
                if variante == "vinkulum_api":
                    bm, bd = [], []
                    for col in range(6):
                        rep = r.reponse(w, forces[:, col])
                        champs[j, :, col] = rep["champ_physique"]
                        nm[j, col] = rep["masse"]["norme_candidate"]
                        nd[j, col] = rep["deformation"]["norme_candidate"]
                        bm.append(rep["masse"]["borne_relative"])
                        bd.append(rep["deformation"]["borne_relative"])
                    bornes["masse"].append(bm); bornes["deformation"].append(bd)
                else:
                    if variante == "vinkulum_reduit":
                        y = solve(r.schur(w), f_normal, assume_a="sym")
                        x = r.reconstruire(w, y)
                    elif variante.startswith("hcb"):
                        x = reponse(r, w, forces, resolution="spectrale_corrigee")
                    else:
                        facteur = splu(ke-w*w*me, permc_spec="MMD_AT_PLUS_A")
                        x = eq[:, None]*facteur.solve(fe)
                        if variante == "lu_corrigee":
                            for _ in range(2):
                                residu = f-d.T@(d@x)+w*w*(m@x)
                                x += eq[:, None]*facteur.solve(eq[:, None]*residu)
                    champs[j] = x
                    nm[j], nd[j] = normes(d, m, x)
            termine = time.perf_counter()
            result.update(statut="termine", preparation_s=prepare-debut,
                          reponses_s=termine-prepare, total_s=termine-debut,
                          rss_avant_oracle_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          avertissements=sorted(set(str(w.message) for w in avertissements)))
        if not np.all(np.isfinite(champs)):
            raise ArithmeticError("champ non fini")
        result["champ_sha256"] = hashlib.sha256(champs.tobytes()).hexdigest()
        phase = "juge"
        ref = np.load(reference, mmap_mode="r", allow_pickle=False)
        result["juge"] = juger(d, m, metrique, champs, ref)
        result["bornes_api"] = bornes if variante == "vinkulum_api" else None
        result["normes_candidates"] = dict(masse=nm.tolist(), deformation=nd.tolist())
    except Exception as erreur:
        result.update(statut="refus_ou_echec", phase=phase,
                      erreur=dict(type=type(erreur).__name__, message=str(erreur)),
                      temps_jusquau_refus_s=time.perf_counter()-debut)
    ecrire(sortie, result)


def analyse(rapport):
    import math
    attendus = {(n, float(hz)) for n in (32, 128, 512) for hz in (2, 20, 40)}
    observes = [(c["n"], c["max_hz"]) for c in rapport["cas"]]
    if len(observes) != len(attendus) or set(observes) != attendus:
        raise ValueError("familles de cas incomplètes ou dupliquées")
    configurations = []
    environnements = set()
    binaires = {}
    for cas in rapport["cas"]:
        if (not math.isfinite(cas["reference"]["ecart_relatif_max"])
                or not 0 <= cas["reference"]["ecart_relatif_max"] <= CIBLE/100):
            raise ValueError("référence non qualifiée")
        if rapport["schema"] >= 2:
            controles = [cas["reference"][cle][nom] for cle in ("maxima", "maxima_operateurs")
                         for nom in ("masse", "deformation", "port")]
            if (not all(math.isfinite(x) and 0 <= x <= CIBLE/100 for x in controles)
                    or max(controles) != cas["reference"]["ecart_relatif_max"]):
                raise ValueError("contrôle opérateur de référence incohérent")
        noms = [c["variante"] for c in cas["configurations"]]
        if len(noms) != len(VARIANTES) or set(noms) != set(VARIANTES):
            raise ValueError("configurations incomplètes ou dupliquées")
        for config in cas["configurations"]:
            essais = config["essais"]
            roles = [(e["role"], e["repetition"]) for e in essais]
            if len(set(roles)) != len(roles):
                raise ValueError("répétition dupliquée")
            if any(role not in ("chauffe", "mesure") for role, _ in roles):
                raise ValueError("rôle de répétition inconnu")
            finaux = [e for e in essais if e["role"] == "mesure"]
            complet = (sorted(e["repetition"] for e in finaux) == [0, 1, 2]
                        and sorted(rep for role, rep in roles if role == "chauffe") == [-1])
            eligible = complet
            for e in essais:
                if e["variante"] != config["variante"] or e["modes_demandes"] != config["modes"]:
                    raise ValueError("configuration de processus incohérente")
                if e["cpu"] != [rapport["cpu"]] or e["fils"] != THREADS:
                    raise ValueError("budget de calcul différent")
                environnements.add((e["python"], e["numpy"], e["scipy"]))
                if e["variante"].startswith(("vinkulum", "hcb")):
                    moteur = "vinkulum" if e["variante"].startswith("vinkulum") else "hcb"
                    identite = (e["moteur"], e["extension_sha256"])
                    if rapport["schema"] >= 2:
                        identite += (json.dumps(e["distribution"], sort_keys=True),)
                        if moteur == "vinkulum":
                            identite += (json.dumps(e["sources_python_sha256"], sort_keys=True),)
                    if moteur in binaires and binaires[moteur] != identite:
                        raise ValueError("binaire différent dans une campagne")
                    binaires[moteur] = identite
                if e["variante"].startswith("hcb") and rapport["protocole"].get("resolution_hcb"):
                    if (e.get("resolution_hcb") != rapport["protocole"]["resolution_hcb"]
                            or e.get("corrections_hcb") != rapport["protocole"]["corrections_hcb"]):
                        raise ValueError("résolution HCB différente du protocole")
                if e["entree_sha256"] != cas["entree_sha256"] or e["reference_sha256"] != cas["reference"]["fichier_sha256"]:
                    raise ValueError("entrées ou référence différentes")
                if e["statut"] != "termine":
                    eligible = False
                    continue
                if e["nombre_frequences"] != cas["nombre_frequences"] or e["charges"] != 6:
                    raise ValueError("couverture de réponses incomplète")
                qualifie = True
                for nom in ("masse", "deformation", "port"):
                    valeurs = e["juge"]["erreurs"][nom]
                    if len(valeurs) != cas["nombre_frequences"] or any(len(v) != 6 for v in valeurs):
                        raise ValueError("couverture du juge incomplète")
                    maxi = max(x for ligne in valeurs for x in ligne)
                    if not all(math.isfinite(x) and x >= 0 for ligne in valeurs for x in ligne):
                        raise ValueError("erreur non finie ou négative")
                    if maxi != e["juge"]["maxima"][nom]:
                        raise ValueError("maximum du juge incohérent")
                    operateurs = e["juge"]["operateurs"][nom]
                    if (len(operateurs) != cas["nombre_frequences"]
                            or not all(math.isfinite(v) and v >= 0 for v in operateurs)):
                        raise ValueError("couverture opérateur incomplète ou invalide")
                    if max(operateurs) != e["juge"]["maxima_operateurs"][nom]:
                        raise ValueError("maximum opérateur incohérent")
                    qualifie &= maxi <= CIBLE and max(operateurs) <= CIBLE
                if e["juge"]["accepte"] != bool(qualifie):
                    raise ValueError("qualification du juge incohérente")
                eligible &= qualifie
                if not all(math.isfinite(e[t]) and e[t] > 0 for t in ("preparation_s", "reponses_s", "total_s")):
                    raise ValueError("temps invalides")
                if abs(e["preparation_s"]+e["reponses_s"]-e["total_s"]) > 1e-9*max(e["total_s"], 1.):
                    raise ValueError("total de temps incohérent")
            summary = dict(n=cas["n"], max_hz=cas["max_hz"], variante=config["variante"],
                           modes=config["modes"], eligible=bool(eligible), complet=complet)
            if complet and all(e["statut"] == "termine" for e in finaux):
                summary.update({t: statistics.median(e[t] for e in finaux)
                                for t in ("preparation_s", "reponses_s", "total_s", "rss_avant_oracle_kib", "directions")})
                if rapport["schema"] >= 2:
                    summary["dispersion_s"] = {t: [min(e[t] for e in finaux), max(e[t] for e in finaux)]
                                               for t in ("preparation_s", "reponses_s", "total_s")}
                summary["erreurs_max"] = {nom: max(e["juge"]["maxima"][nom] for e in finaux)
                                           for nom in ("masse", "deformation", "port")}
                summary["erreurs_operateur_max"] = {nom: max(e["juge"]["maxima_operateurs"][nom] for e in finaux)
                                                     for nom in ("masse", "deformation", "port")}
            configurations.append(summary)
    if len(environnements) != 1:
        raise ValueError("versions de Python/NumPy/SciPy différentes")
    return dict(configurations=configurations, cible=CIBLE,
                portee="routines matricielles sur consoles ; aucune simulation FFRF ou comparaison généraliste")


def campagne(dossier, python_vinkulum, python_exudyn, cpu):
    """Sondages HCB puis quatre passages alternés, le premier hors médianes."""
    dossier = Path(dossier)
    env = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE="1")
    rapport = dict(schema=2, date_utc=datetime.now(timezone.utc).isoformat(),
        cpu=cpu, plateforme=platform.platform(), cible=CIBLE, cas=[],
        protocole=dict(repetitions=3, echauffements=1, frequences=257, charges=6,
            imports_et_entrees_exclus=True, assemblage_gram_comptabilise=True,
            reponses_comprennent_champs_et_normes=True,
            diagnostic_hcb_comptabilise="trace et dimensions ; audits complets exclus des runs",
            resolution_hcb="spectrale_corrigee", corrections_hcb=2,
            selection_modes_hcb_hors_cout=True,
            description="Ordre inversé un passage sur deux, processus frais, un CPU et un fil. "
                        "HCB : premier compte de la grille satisfaisant le juge commun au sondage ; "
                        "réponse spectrale avec deux corrections résiduelles dans la petite paire équilibrée, "
                        "tous leurs coûts de préparation et de réponse sont comptés ; "
                        "les projections alternatives et résidus modaux de diagnostic sont désactivés ; "
                        "cette recherche avec oracle est offerte à HCB. Vinkulum : construction "
                        "automatique au seuil Schur1e-10, avec certificat spectral compté. "
                        "Les bornes de champ internes sont rapportées séparément de la qualification."),
        sources_sha256={nom: sha(ROOT/"ci"/nom) for nom in
            ("confronte_ports_exudyn.py", "reference_hcb_exudyn.py", "oracle_champs_ports.py",
             "reference_ports_precision.py", "experience_reduction_native.py")})
    def sauvegarder():
        ecrire(dossier/"rapport.json", rapport)
    def lancer(cas, config, role, repetition):
        nom = f"n{cas['n']}-f{cas['max_hz']:g}"
        variante, modes = config["variante"], config["modes"]
        base = f"{nom}-{variante}-r{modes}-{role}{repetition}"
        sortie = dossier/(base+".json")
        py = python_exudyn if variante.startswith("hcb") else python_vinkulum
        commande = [str(py), str(Path(__file__).resolve()), "--worker", str(dossier/(nom+".npz")),
                    str(dossier/(nom+".ref.npy")), variante, str(modes), str(sortie), str(cpu)]
        print(base, flush=True)
        debut = time.perf_counter()
        with (dossier/(base+".log")).open("w") as journal:
            retour = subprocess.run(commande, env=env, stdout=journal, stderr=subprocess.STDOUT)
        if retour.returncode != 0 or not sortie.exists():
            raise RuntimeError(f"processus {base} en échec ({retour.returncode}), voir son journal")
        resultat = json.loads(sortie.read_text())
        resultat.update(role=role, repetition=repetition, processus_s=time.perf_counter()-debut,
                         journal_sha256=sha(dossier/(base+".log")))
        return resultat
    for n in (32, 128, 512):
        for hz in (2., 20., 40.):
            nom = f"n{n}-f{hz:g}"
            cas = json.loads((dossier/(nom+".npz.json")).read_text())
            cas.update(reference=json.loads((dossier/(nom+".ref.npy.json")).read_text()),
                       sondages=[], configurations=[])
            if cas["reference"]["ecart_relatif_max"] > CIBLE/100:
                raise ValueError("référence non qualifiée")
            if sha(dossier/(nom+".npz")) != cas["entree_sha256"]:
                raise ValueError("entrée altérée après génération")
            rapport["cas"].append(cas)
            for variante in VARIANTES:
                modes = 0
                if variante.startswith("hcb"):
                    grille = sorted(set([r for r in (0, 12, 24, 48, 96, 192, 384, 768)
                                         if r < 6*(n-1)] + [min(6*(n-1)-1, 768)]))
                    for rang in grille:
                        modes = rang
                        essai = lancer(cas, dict(variante=variante, modes=rang), "sondage", 0)
                        cas["sondages"].append(essai)
                        sauvegarder()
                        if essai["statut"] == "termine" and essai["juge"]["accepte"]:
                            break
                cas["configurations"].append(dict(variante=variante, modes=modes, essais=[]))
            for passage in range(4):
                ordre = cas["configurations"] if passage % 2 == 0 else cas["configurations"][::-1]
                for config in ordre:
                    config["essais"].append(lancer(cas, config, "chauffe" if passage == 0 else "mesure",
                                                   -1 if passage == 0 else passage-1))
                    sauvegarder()
    rapport["analyse"] = analyse(rapport)
    sauvegarder()
    print("Campagne terminée : "+str(dossier/"rapport.json"), flush=True)


def archiver(dossier, destination):
    """Conserve les mesures, entrées natives et sources ; oracles recalculables."""
    import shutil
    dossier, destination = Path(dossier), Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    rapport = json.loads((dossier/"rapport.json").read_text())
    bilan = analyse(rapport)
    if len(rapport["cas"]) != 9:
        raise ValueError("campagne incomplète")
    contenu = json.dumps(rapport, ensure_ascii=False, allow_nan=False).encode()
    (destination/"essais.json.gz").write_bytes(gzip.compress(contenu, mtime=0))
    ecrire(destination/"bilan.json", bilan)
    for nom, empreinte in rapport["sources_sha256"].items():
        if sha(ROOT/"ci"/nom) != empreinte:
            raise ValueError("source modifiée depuis les mesures : "+nom)
        cible = destination/"sources"/nom
        cible.parent.mkdir(exist_ok=True)
        shutil.copy2(ROOT/"ci"/nom, cible)
    for source in sorted(dossier.glob("*.npz")):
        cible = destination/"entrees"/source.name
        cible.parent.mkdir(exist_ok=True)
        shutil.copy2(source, cible)
    journaux = {p.name: p.read_text() for p in dossier.glob("*.log")}
    (destination/"journaux.json.gz").write_bytes(gzip.compress(
        json.dumps(journaux, ensure_ascii=False).encode(), mtime=0))
    ecrire(destination/"manifest.json", dict(schema=1,
        date_utc=datetime.now(timezone.utc).isoformat(), version_vinkulum="0.10.0",
        reference="DᵀD exact des coefficients binary64, Decimal70/90 ; champs complets recalculables",
        fichiers_sha256={str(p.relative_to(destination)): sha(p)
                         for p in sorted(destination.rglob("*")) if p.is_file() and p.name != "manifest.json"}))


def verifier(destination):
    destination = Path(destination)
    manifest = json.loads((destination/"manifest.json").read_text())
    declares = set(manifest["fichiers_sha256"])
    presents = {str(p.relative_to(destination)) for p in destination.rglob("*")
                if p.is_file() and p.name != "manifest.json"}
    if declares != presents:
        raise ValueError("inventaire de l'archive incomplet ou incohérent")
    indispensables = {"essais.json.gz", "bilan.json", "journaux.json.gz"}
    indispensables.update("sources/"+nom for nom in
        ("confronte_ports_exudyn.py", "reference_hcb_exudyn.py", "oracle_champs_ports.py",
         "reference_ports_precision.py", "experience_reduction_native.py"))
    indispensables.update(f"entrees/n{n}-f{hz}.npz" for n in (32, 128, 512) for hz in (2, 20, 40))
    if not indispensables <= declares:
        raise ValueError("fichiers indispensables de l'archive absents")
    for nom, empreinte in manifest["fichiers_sha256"].items():
        if sha(destination/nom) != empreinte:
            raise ValueError("archive altérée : "+nom)
    rapport = json.loads(gzip.decompress((destination/"essais.json.gz").read_bytes()))
    if analyse(rapport) != json.loads((destination/"bilan.json").read_text()):
        raise ValueError("bilan incohérent")
    sources = {"confronte_ports_exudyn.py", "reference_hcb_exudyn.py", "oracle_champs_ports.py",
               "reference_ports_precision.py", "experience_reduction_native.py"}
    if set(rapport["sources_sha256"]) != sources:
        raise ValueError("sources de reproduction incomplètes")
    for nom, empreinte in rapport["sources_sha256"].items():
        if sha(destination/"sources"/nom) != empreinte:
            raise ValueError("provenance de source incohérente : "+nom)
    journaux = json.loads(gzip.decompress((destination/"journaux.json.gz").read_bytes()))
    for cas in rapport["cas"]:
        nom = f"n{cas['n']}-f{cas['max_hz']:g}"
        if sha(destination/"entrees"/(nom+".npz")) != cas["entree_sha256"]:
            raise ValueError("provenance d'entrée incohérente : "+nom)
        essais = cas["sondages"]+[e for c in cas["configurations"] for e in c["essais"]]
        for e in essais:
            journal = (f"{nom}-{e['variante']}-r{e['modes_demandes']}-"
                       f"{e['role']}{e['repetition']}.log")
            if journal not in journaux or hashlib.sha256(journaux[journal].encode()).hexdigest() != e["journal_sha256"]:
                raise ValueError("provenance de journal incohérente : "+journal)
            if e["variante"].startswith("vinkulum") and (e["moteur"], e["extension_sha256"]) != (cas["vinkulum"], cas["extension_sha256"]):
                raise ValueError("binaire de génération différent des essais")
    print("Archive de confrontation HCB conforme : "+str(destination))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--worker", nargs=6, metavar=("ENTREE", "REFERENCE", "VARIANTE", "MODES", "SORTIE", "CPU"))
    p.add_argument("--entree", nargs=3, metavar=("N", "MAX_HZ", "SORTIE"))
    p.add_argument("--reference", nargs=2, metavar=("ENTREE", "SORTIE"))
    p.add_argument("--campagne", nargs=3, metavar=("DOSSIER", "PY_VINKULUM", "PY_EXUDYN"))
    p.add_argument("--cpu", type=int, default=8)
    p.add_argument("--archiver", nargs=2, metavar=("DOSSIER", "DESTINATION"))
    p.add_argument("--verifier", type=Path)
    args = p.parse_args()
    if args.worker:
        entree, reference, variante, modes, sortie, cpu = args.worker
        os.sched_setaffinity(0, {int(cpu)})
        if variante not in VARIANTES:
            p.error("variante inconnue")
        worker(entree, reference, variante, int(modes), sortie)
    elif args.entree:
        n, hz, sortie = args.entree
        ecrire(str(sortie)+".json", entrees(int(n), float(hz), sortie))
    elif args.reference:
        entree, sortie = args.reference
        ecrire(str(sortie)+".json", references(entree, sortie))
    elif args.campagne:
        campagne(*args.campagne, args.cpu)
    elif args.archiver:
        archiver(*args.archiver)
    elif args.verifier:
        verifier(args.verifier)
    else:
        p.error("choisir une opération")


if __name__ == "__main__":
    main()
