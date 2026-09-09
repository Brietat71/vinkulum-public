"""Comparaison locale Krylov adaptatif / HCB creux à tolérance de Schur commune.

Préparation, enveloppe uniforme, requêtes et une reconstruction sont comptées.
Imports, assemblage de la fixture et audit indépendant sont exclus des durées.
HCB reçoit un nombre de modes présélectionné : le coût de sa recherche n'est
pas compté. Aucun temps n'est attribué à Exudyn, MBDyn ou Simpack.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
THREADS = {k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS")}
SOURCES = ("mesure_ports_krylov.py", "ports_krylov.py", "reference_ports_hcb.py",
           "modeles_ports.py", "test_ports_krylov.py", "test_reference_ports_hcb.py")
# Les comptes HCB sont figés après sondage, avant la campagne archivée.
CASES = [
    dict(nom="chaine32", famille="chaine", n=32, fraction=.8, modes=28, tolerance=1e-6),
    dict(nom="chaine256", famille="chaine", n=256, fraction=.8, modes=88, tolerance=1e-6),
    dict(nom="chaine2048", famille="chaine", n=2048, fraction=.8, modes=96, tolerance=1e-6),
    dict(nom="chaine2048_pole", famille="chaine", n=2048, fraction=.98, modes=256, tolerance=1e-6),
    dict(nom="console16", famille="console", n=16, fraction=.7, modes=32, tolerance=1e-6),
    dict(nom="console128", famille="console", n=128, fraction=.7, modes=56, tolerance=1e-6),
    dict(nom="console128_arrondi", famille="console", n=128, fraction=.7, modes=762, tolerance=1e-10),
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ecrit(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+"\n")


def worker(nom, variante, output):
    import resource
    import numpy as np
    import scipy
    from scipy.sparse.linalg import spsolve
    import vinkulum
    from modeles_ports import chaine, console, borne_lambda_trace_console, metrique_ports, reference_chaine
    from ports_krylov import reduire as krylov
    from reference_ports_hcb import reduire as hcb
    cas = next(c for c in CASES if c["nom"] == nom)
    modele = chaine(cas["n"]) if cas["famille"] == "chaine" else console(cas["n"])
    metrique = metrique_ports(modele)
    if cas["famille"] == "chaine":
        lam = modele.metadata["mu_min_interieur_s_moins_2"]
        origine_borne = "spectre analytique chaîne Dirichlet"
    else:
        lam = borne_lambda_trace_console(modele.metadata)
        origine_borne = "inverse trace flexibilité massique analytique"
    omega_max = cas["fraction"]*np.sqrt(lam)
    lam *= .999  # marge conventionnelle ; ne certifie pas les arrondis.
    frequences = np.linspace(0., omega_max, 257)
    debut = time.perf_counter()
    if variante == "krylov":
        r = krylov(modele.k, modele.m, modele.interieur, modele.interface, metrique,
                   lambda_min=lam, omega_max=omega_max, tolerance=cas["tolerance"],
                   max_blocs=12, max_directions=128)
        borne = r.interieur.borne_uniforme
        diagnostic = dict(statut=r.interieur.statut, historique=r.interieur.historique,
                          resolutions=r.interieur.resolutions, seconds_membres=r.interieur.seconds_membres,
                          nnz_facteurs=r.interieur.nnz_facteurs,
                          residu_resolution_max=r.interieur.residu_resolution_max,
                          defaut_gram=r.interieur.defaut_gram)
        w = r.normalisation
    else:
        r = hcb(modele.k, modele.m, modele.interieur, modele.interface, cas["modes"], metrique)
        borne = r.borne_uniforme(lam, omega_max)
        diagnostic = dict(r.compteurs)
        diagnostic["bilan_borne"] = r.bilan_borne
        w = r.w
    prepare = time.perf_counter()
    schurs = np.array([r.schur(om) for om in frequences])
    requetes = time.perf_counter()
    om = float(frequences[-1])
    port = np.eye(len(modele.interface))
    champ = r.reconstruire(om, port) if variante == "krylov" else r.reconstruit(om, port)
    termine = time.perf_counter()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Audit indépendant après la zone chronométrée. Erreur en norme spectrale,
    # avec métrique statique physique, pas une erreur relative au Schur résonant.
    i, s = modele.interieur, modele.interface
    erreurs = []
    ecarts_loewner = []
    for om, sr in zip(frequences, schurs):
        if cas["famille"] == "chaine":
            ref = np.array([[reference_chaine(cas["n"], om)["schur"]*cas["n"]]])
        else:
            a = modele.k-om*om*modele.m
            b = a[i][:, s].toarray()
            ref = w.T@(a[s][:, s].toarray()-b.T@spsolve(a[i][:, i], b))@w
        diff = (sr-ref+sr.T-ref.T)*.5
        erreurs.append(float(np.linalg.norm(diff, 2)))
        ecarts_loewner.append(float(np.linalg.eigvalsh(diff)[0]))
    a = modele.k-omega_max**2*modele.m
    residual = (a@champ)[i]
    force = (a@champ)[s]
    fermeture = float(np.linalg.norm(residual, np.inf)/max(np.linalg.norm(force, np.inf), 1e-300))
    erreur = max(erreurs)
    accepte_borne = borne <= cas["tolerance"]
    accepte_audit = erreur <= cas["tolerance"]
    extension = next(Path(vinkulum.__file__).parent.glob("_vinkulum*.so"), None)
    if extension is None:
        extension = next(Path(vinkulum.__file__).parent.glob("*.so"))
    result = dict(cas=cas, variante=variante, date=datetime.now(timezone.utc).isoformat(),
                  fils={k: os.environ.get(k) for k in THREADS}, cpu=sorted(os.sched_getaffinity(0)),
                  python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                  extension_sha256=sha(extension), metadata=modele.metadata,
                  lambda_min=lam, origine_borne=origine_borne, omega_max=omega_max,
                  directions=r.taille_interieure, ports=len(s), nombre_requetes=len(frequences),
                  preparation_s=prepare-debut, requetes_s=requetes-prepare,
                  reconstruction_s=termine-requetes, total_s=termine-debut,
                  rss_avant_oracle_kib=rss, diagnostic=diagnostic, borne_uniforme=borne,
                  erreur_schur_audit=erreur, minimum_ecart_loewner=min(ecarts_loewner),
                  residu_champ_relatif=fermeture, accepte_borne=accepte_borne,
                  accepte_audit=accepte_audit, accepte=accepte_borne and accepte_audit,
                  certification_machine=False)
    ecrit(output, result)


def campaign(path, cpu):
    if path.exists() or cpu not in os.sched_getaffinity(0):
        raise SystemExit("sortie existante ou CPU indisponible")
    (path/"sources").mkdir(parents=True)
    hashes = {}
    for name in SOURCES:
        (path/"sources"/name).write_bytes((ROOT/"ci"/name).read_bytes())
        hashes[name] = sha(path/"sources"/name)
    essais = []
    for cas in CASES:
        for rep in range(4):
            variantes = ("krylov", "hcb") if rep % 2 == 0 else ("hcb", "krylov")
            for variante in variantes:
                name = f'{cas["nom"]}-{variante}-{rep}.json'
                subprocess.run(["taskset", "-c", str(cpu), sys.executable, __file__, "--worker",
                                cas["nom"], variante, str(path/name)],
                               env=dict(os.environ, **THREADS), check=True, timeout=180)
                data = json.loads((path/name).read_text())
                essais.append(dict(fichier=name, sha256=sha(path/name), cas=cas["nom"],
                                   variante=variante, repetition=rep, echauffement=rep == 0))
                print(cas["nom"], variante, rep, data["directions"], data["accepte"],
                      f'{data["total_s"]:.5g}s', flush=True)
    ecrit(path/"manifest.json", dict(schema=1, cas=CASES, repetitions=4, essais=essais,
                                   sources=hashes, date=datetime.now(timezone.utc).isoformat(),
                                   cpu_info=Path("/proc/cpuinfo").read_text().split("model name\t: ")[1].split("\n")[0]))
    verify(path)


def verify(path):
    manifest = json.loads((path/"manifest.json").read_text())
    for name, expected in manifest["sources"].items():
        assert sha(path/"sources"/name) == expected
    expected = {(c["nom"], v, r) for c in manifest["cas"]
                for v in ("krylov", "hcb") for r in range(manifest["repetitions"])}
    found, raws = set(), []
    for entry in manifest["essais"]:
        key = (entry["cas"], entry["variante"], entry["repetition"])
        assert key in expected and key not in found
        found.add(key)
        assert sha(path/entry["fichier"]) == entry["sha256"]
        data = json.loads((path/entry["fichier"]).read_text())
        assert data["cas"] in manifest["cas"] and data["cas"]["nom"] == entry["cas"]
        assert data["variante"] == entry["variante"]
        assert data["fils"] == THREADS and len(data["cpu"]) == 1
        assert entry["echauffement"] == (entry["repetition"] == 0)
        assert data["certification_machine"] is False
        assert 0 <= data["omega_max"]**2 < data["lambda_min"]
        assert 0 <= data["directions"] <= data["metadata"]["ddl"]-data["ports"]
        for name in ("total_s", "preparation_s", "requetes_s", "reconstruction_s",
                     "borne_uniforme", "erreur_schur_audit", "residu_champ_relatif"):
            assert math.isfinite(data[name]) and data[name] >= 0
        assert math.isclose(sum(data[k] for k in ("preparation_s", "requetes_s", "reconstruction_s")),
                            data["total_s"], rel_tol=1e-12)
        assert data["accepte_borne"] == (data["borne_uniforme"] <= data["cas"]["tolerance"])
        assert data["accepte_audit"] == (data["erreur_schur_audit"] <= data["cas"]["tolerance"])
        assert data["accepte"] == (data["accepte_borne"] and data["accepte_audit"])
        raws.append((entry, data))
    assert found == expected
    for c in manifest["cas"]:
        ligne = []
        for v in ("krylov", "hcb"):
            trials = [d for e, d in raws if e["cas"] == c["nom"] and e["variante"] == v and not e["echauffement"]]
            ligne.append(dict(variante=v, directions=[d["directions"] for d in trials],
                              accepte=[d["accepte"] for d in trials],
                              total_ms=1000*statistics.median(d["total_s"] for d in trials),
                              borne=max(d["borne_uniforme"] for d in trials),
                              erreur=max(d["erreur_schur_audit"] for d in trials)))
        print(c["nom"], ligne)
    print(f"Archive vérifiée : {len(found)} essais, échauffements et refus conservés.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--worker", nargs=3)
    group.add_argument("--sortie", type=Path)
    group.add_argument("--verifier", type=Path)
    parser.add_argument("--cpu", type=int, default=8)
    args = parser.parse_args()
    if args.worker:
        worker(*args.worker)
    elif args.verifier:
        verify(args.verifier)
    else:
        campaign(args.sortie, args.cpu)
