"""Audit de précision des ports relevés et de leur assemblage.

Les oracles Decimal sont préparés avant les mesures et archivés séparément.
Chaque mesure part dans un processus frais : construction, contrôle uniforme,
257 requêtes et un champ ; imports, modèles et oracles sont exclus du temps.
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
CASES = [dict(nom="chaine2048", famille="chaine", n=2048)] + [
    dict(nom=f"console{n}", famille="console", n=n) for n in (32, 128, 512)]
VARIANTES = ("direct", "qr", "lu_energie")
SOURCES = ("mesure_ports_releves.py", "ports_releves.py", "ports_krylov.py",
           "condensation_energie.py", "condensation_energie_lu.py", "energie_ports.py",
           "modeles_ports.py", "reference_ports_precision.py", "chaine_sous_structures.py",
           "test_energie_ports.py", "test_ports_releves.py", "test_reference_ports_precision.py")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ecrit(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+"\n")


def modele(cas):
    import numpy as np
    from modeles_ports import chaine, console, metrique_ports, borne_lambda_trace_console
    from energie_ports import energie_chaine, energie_console
    a = chaine(cas["n"]) if cas["famille"] == "chaine" else console(cas["n"])
    e = energie_chaine(cas["n"]) if cas["famille"] == "chaine" else energie_console(a.metadata)
    lam = (a.metadata["mu_min_interieur_s_moins_2"] if cas["famille"] == "chaine"
           else borne_lambda_trace_console(a.metadata))
    return a, e, metrique_ports(a), .999*lam, .7*np.sqrt(lam)


def empreinte_matrice(a):
    a = a.tocsr()
    h = hashlib.sha256()
    for x in (a.indptr, a.indices, a.data):
        h.update(x.tobytes())
    h.update(str(a.shape).encode())
    return h.hexdigest()


def reference(cas, path):
    import numpy as np
    from scipy.linalg import cholesky, solve_triangular
    from reference_ports_precision import schur_console, schur_facteur
    from modeles_ports import reference_chaine
    a, e, metric, lam, om = modele(cas)
    w = solve_triangular(cholesky(metric, lower=True).T, np.eye(len(e.interface)), lower=False)
    result = []
    for omega in np.linspace(0., om, 17):
        sd = w.T@schur_facteur(e.d, e.m, e.interieur, e.interface, omega)@w
        if cas["famille"] == "console":
            sp = w.T@schur_console(a.metadata, omega)@w
        else:
            sp = np.array([[cas["n"]*reference_chaine(cas["n"], omega)["schur"]]])
        result.append(dict(omega=float(omega), schur_D=sd.tolist(), schur_physique=sp.tolist()))
    s90 = w.T@schur_facteur(e.d, e.m, e.interieur, e.interface, om, dps=90)@w
    ecart_precision = float(np.linalg.norm(s90-np.array(result[-1]["schur_D"]), 2))
    ecrit(path, dict(cas=cas, lambda_min=lam, omega_max=om, d_sha256=empreinte_matrice(e.d),
                     m_sha256=empreinte_matrice(e.m), metrique=metric.tolist(), metadata=a.metadata,
                     dps=70, dps_verification=90, ecart_70_90=ecart_precision, frequences=result))


def worker(nom, variante, oracle, output):
    import resource
    import numpy as np
    import scipy
    import vinkulum
    from ports_releves import PortsReleves
    from ports_krylov import reduire
    cas = next(c for c in CASES if c["nom"] == nom)
    a, e, metric, lam, om = modele(cas)
    frequences = np.linspace(0., om, 257)
    debut = time.perf_counter()
    if variante == "direct":
        r = reduire(a.k, a.m, a.interieur, a.interface, metric, lambda_min=lam,
                    omega_max=om, tolerance=1e-10)
        borne = r.interieur.borne_uniforme
        diagnostic = dict(statut=r.interieur.statut, historique=r.interieur.historique)
    else:
        r = PortsReleves(e.d, e.m, e.interieur, e.interface, metric, lam, om,
                        tolerance=1e-10, methode=variante)
        borne = r.borne_uniforme
        diagnostic = dict(statut=r.statut, historique=r.interieur.historique,
                          audit_bande=r.audit_bande, borne_reduction=r.borne_qr)
    prepare = time.perf_counter()
    _ = [r.schur(omega) for omega in frequences]
    requetes = time.perf_counter()
    champ = r.reconstruire(om, np.eye(len(e.interface)))
    termine = time.perf_counter()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ref = json.loads(Path(oracle).read_text())
    assert ref["cas"] == cas and ref["lambda_min"] == lam and ref["omega_max"] == om
    assert ref["d_sha256"] == empreinte_matrice(e.d) and ref["m_sha256"] == empreinte_matrice(e.m)
    erreurs_d, erreurs_p = [], []
    for point in ref["frequences"]:
        sr = r.schur(point["omega"])
        erreurs_d.append(float(np.linalg.norm(sr-np.array(point["schur_D"]), 2)))
        erreurs_p.append(float(np.linalg.norm(sr-np.array(point["schur_physique"]), 2)))
    dx, mx = e.d@champ, e.m@champ
    energie = dx.T@dx-om**2*(champ.T@mx)
    ecart_energie = float(np.linalg.norm(energie-r.schur(om), 2))
    if variante != "direct":
        diagnostic["audit_modele"] = r.audit_modele()
        au = r.audit_energie(om)
        diagnostic["audit_ponctuel_total"] = au["majorant_total"]
    extension = next(Path(vinkulum.__file__).parent.glob("*.so"))
    resultat = dict(cas=cas, variante=variante, oracle_sha256=sha(oracle),
                     date=datetime.now(timezone.utc).isoformat(), python=sys.version,
                     numpy=np.__version__, scipy=scipy.__version__, extension_sha256=sha(extension),
                     fils={k: os.environ.get(k) for k in THREADS}, cpu=sorted(os.sched_getaffinity(0)),
                     tolerance=1e-10, lambda_min=lam, omega_max=om, directions=r.taille_interieure,
                     ports=len(e.interface), nombre_requetes=257, nombre_frequences_audit=17,
                     preparation_s=prepare-debut, requetes_s=requetes-prepare,
                     reconstruction_s=termine-requetes, total_s=termine-debut,
                     rss_avant_oracle_kib=rss, borne_uniforme=borne,
                     erreur_D=max(erreurs_d), erreur_physique=max(erreurs_p),
                     ecart_energie=ecart_energie, diagnostic=diagnostic,
                     accepte=bool(borne <= 1e-10 and max(erreurs_d) <= 1e-10 and max(erreurs_p) <= 1e-10),
                     certification_machine=False)
    ecrit(output, resultat)


def assemblage(n, segments, output):
    import numpy as np
    from chaine_sous_structures import ChaineAssemblee
    from modeles_ports import chaine, reference_chaine
    lam = chaine(n).metadata["mu_min_interieur_s_moins_2"]
    om = .8*np.sqrt(lam)
    debut = time.perf_counter()
    r = ChaineAssemblee(n, segments, om, tolerance=1e-10)
    prepare = time.perf_counter()
    frequences = np.linspace(0., om, 129)
    reponses = [r.reponse(omega) for omega in frequences]
    requetes = time.perf_counter()
    champ = r.reconstruire(om, reponses[-1]["deplacement"])
    termine = time.perf_counter()
    erreurs, rapports, relatifs, non_bornes = [], [], [], 0
    indices = np.arange(r.longueur-1, n, r.longueur)
    for omega, sol in zip(frequences, reponses):
        ref = reference_chaine(n, omega, force=1/r.echelle)["deplacement"]
        exact = ref[indices]/r.echelle
        err = float(np.linalg.norm(sol["deplacement"]-exact))
        erreurs.append(err)
        relatifs.append(err/max(float(np.linalg.norm(exact)), 1e-300))
        if sol["borne_erreur_norme"] is None:
            non_bornes += 1
        else:
            rapports.append(err/max(sol["borne_erreur_norme"], 1e-300))
    ref = reference_chaine(n, om, force=1/r.echelle)["deplacement"]
    ecrit(output, dict(n=n, segments=segments, variante="assemblage_qr",
                        fils={k: os.environ.get(k) for k in THREADS}, cpu=sorted(os.sched_getaffinity(0)),
                        preparation_s=prepare-debut, requetes_s=requetes-prepare,
                        reconstruction_s=termine-requetes, total_s=termine-debut,
                        directions=r.taille_interieure, ports=segments,
                        borne_uniforme=r.borne_uniforme, erreur_port_max=max(erreurs),
                        erreur_relative_port_max=max(relatifs), rapport_erreur_borne_max=max(rapports, default=0.),
                        nombre_frequences=129, reponses_non_bornees=non_bornes,
                        erreur_champ_final=float(np.linalg.norm(champ-ref)/np.linalg.norm(ref)),
                        accepte=bool(r.borne_uniforme <= 1e-10 and non_bornes == 0
                                      and max(rapports, default=0.) <= 1.),
                        certification_machine=False))


def campaign(path, cpu):
    if path.exists() or cpu not in os.sched_getaffinity(0):
        raise SystemExit("sortie existante ou CPU indisponible")
    (path/"sources").mkdir(parents=True)
    (path/"references").mkdir()
    hashes = {}
    for name in SOURCES:
        (path/"sources"/name).write_bytes((ROOT/"ci"/name).read_bytes())
        hashes[name] = sha(path/"sources"/name)
    references = {}
    for c in CASES:
        name = c["nom"]+".json"
        reference(c, path/"references"/name)
        references[name] = sha(path/"references"/name)
    essais = []
    env = dict(os.environ, **THREADS)
    for c in CASES:
        for rep in range(4):
            variantes = VARIANTES[rep % 3:]+VARIANTES[:rep % 3]
            for variante in variantes:
                name = f'{c["nom"]}-{variante}-{rep}.json'
                subprocess.run(["taskset", "-c", str(cpu), sys.executable, __file__, "--worker",
                                c["nom"], variante, str(path/"references"/(c["nom"]+".json")), str(path/name)],
                               env=env, check=True, timeout=180)
                essais.append(dict(fichier=name, sha256=sha(path/name), cas=c["nom"],
                                   variante=variante, repetition=rep, echauffement=rep == 0))
                d = json.loads((path/name).read_text())
                print(c["nom"], variante, rep, d["directions"], d["accepte"], d["erreur_physique"], flush=True)
    for n, seg in ((64, 4), (2048, 8), (2048, 32)):
        for rep in range(4):
            name = f"assemblage-{n}-{seg}-{rep}.json"
            subprocess.run(["taskset", "-c", str(cpu), sys.executable, __file__, "--assemblage",
                            str(n), str(seg), str(path/name)], env=env, check=True, timeout=180)
            essais.append(dict(fichier=name, sha256=sha(path/name), cas=f"assemblage-{n}-{seg}",
                               variante="assemblage_qr", repetition=rep, echauffement=rep == 0))
            print(name, flush=True)
    ecrit(path/"manifest.json", dict(schema=1, date=datetime.now(timezone.utc).isoformat(),
                                    cas=CASES, variantes=VARIANTES, repetitions=4,
                                    assemblages=[[64,4],[2048,8],[2048,32]], sources=hashes,
                                    references=references, essais=essais,
                                    cpu_info=Path("/proc/cpuinfo").read_text().split("model name\t: ")[1].split("\n")[0]))
    verify(path)


def verify(path):
    m = json.loads((path/"manifest.json").read_text())
    for name, h in m["sources"].items():
        assert sha(path/"sources"/name) == h
    for name, h in m["references"].items():
        assert sha(path/"references"/name) == h
        d = json.loads((path/"references"/name).read_text())
        assert d["dps"] == 70 and d["dps_verification"] == 90 and d["ecart_70_90"] <= 1e-14
    expected = {(c["nom"], v, r) for c in m["cas"] for v in m["variantes"] for r in range(m["repetitions"])}
    expected |= {(f"assemblage-{n}-{s}", "assemblage_qr", r) for n, s in m["assemblages"] for r in range(m["repetitions"])}
    seen, donnees = set(), []
    for e in m["essais"]:
        key = (e["cas"], e["variante"], e["repetition"])
        assert key in expected and key not in seen
        seen.add(key)
        assert sha(path/e["fichier"]) == e["sha256"]
        d = json.loads((path/e["fichier"]).read_text())
        assert d["variante"] == e["variante"] and d["fils"] == THREADS and len(d["cpu"]) == 1
        assert e["echauffement"] == (e["repetition"] == 0) and d["certification_machine"] is False
        for name in ("total_s", "preparation_s", "requetes_s", "reconstruction_s", "borne_uniforme"):
            assert math.isfinite(d[name]) and d[name] >= 0
        assert math.isclose(sum(d[k] for k in ("preparation_s", "requetes_s", "reconstruction_s")), d["total_s"], rel_tol=1e-12)
        if d["variante"] == "assemblage_qr":
            assert d["accepte"] == (d["borne_uniforme"] <= 1e-10 and d["reponses_non_bornees"] == 0 and d["rapport_erreur_borne_max"] <= 1)
        else:
            assert d["cas"]["nom"] == e["cas"] and d["tolerance"] == 1e-10
            assert d["oracle_sha256"] == m["references"][e["cas"]+".json"]
            assert d["accepte"] == (d["borne_uniforme"] <= d["tolerance"] and d["erreur_D"] <= d["tolerance"] and d["erreur_physique"] <= d["tolerance"])
        donnees.append((e, d))
    assert seen == expected
    for cas, variante in sorted({(e["cas"], e["variante"]) for e in m["essais"]}):
        ds = [d for e, d in donnees if e["cas"] == cas and e["variante"] == variante and not e["echauffement"]]
        print(cas, variante, "r", [d["directions"] for d in ds], "accepte", [d["accepte"] for d in ds],
              "ms", 1000*statistics.median(d["total_s"] for d in ds),
              "erreur", max(d.get("erreur_physique", d.get("erreur_relative_port_max", 0.)) for d in ds))
    print(f"Archive vérifiée : {len(seen)} essais, références haute précision et refus conservés.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--worker", nargs=4)
    g.add_argument("--assemblage", nargs=3)
    g.add_argument("--sortie", type=Path)
    g.add_argument("--verifier", type=Path)
    p.add_argument("--cpu", type=int, default=8)
    a = p.parse_args()
    if a.worker:
        worker(*a.worker)
    elif a.assemblage:
        assemblage(int(a.assemblage[0]), int(a.assemblage[1]), a.assemblage[2])
    elif a.verifier:
        verify(a.verifier)
    else:
        campaign(a.sortie, a.cpu)
