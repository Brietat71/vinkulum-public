"""Expériences de précision physique ; aucune comparaison chronométrée.

Les fréquences d'acceptation sont déclarées avant calcul. Une solution
analytique contrôle a posteriori chaque champ ; elle ne décide pas l'arrêt.
Les échecs près de la résonance restent archivés avec les réussites.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil

import numpy as np
import scipy

from chaine_sous_structures import ChaineAssemblee
from champ_interieur import ControleChamp, adapter_champ
from modeles_ports import reference_chaine

SOURCES = ("experience_champ_interieur.py", "champ_interieur.py", "ports_releves.py",
           "ports_krylov.py", "chaine_sous_structures.py", "condensation_energie.py",
           "condensation_energie_lu.py", "energie_ports.py", "reference_ports_precision.py",
           "modeles_ports.py", "test_ports_releves.py", "test_champ_interieur.py")
CAS = [(64, 4, 3), (2048, 8, 3), (2048, 32, 3), (2048, 8, 17), (2048, 32, 17)]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ecrire(p, valeur):
    Path(p).write_text(json.dumps(valeur, indent=2, ensure_ascii=False, allow_nan=False)+"\n")


def precision(a, frequences, reponses):
    valeurs = []
    for omega, r in zip(frequences, reponses):
        u = a.reconstruire(omega, r["deplacement"])
        oracle = reference_chaine(a.n, omega, force=1/a.echelle)["deplacement"]
        erreur = u-oracle
        du, de = np.diff(np.r_[0., oracle]), np.diff(np.r_[0., erreur])
        resultat = dict(omega=float(omega), norme_oracle=float(np.linalg.norm(oracle)),
                        norme_deformation_oracle=float(np.linalg.norm(du)),
                        erreur_masse_relative=float(np.linalg.norm(erreur)/np.linalg.norm(oracle)),
                        erreur_deformation_relative=float(np.linalg.norm(de)/np.linalg.norm(du)),
                        borne_masse_relative=r["masse"]["borne_relative"],
                        borne_deformation_relative=r["deformation"]["borne_relative"],
                        marge=r["marge"])
        valeurs.append(resultat)
    return valeurs


def experience(n, segments, nombre_frequences):
    om = .8*2*np.sin(np.pi/(2*n))
    frequences = (np.array([0., .37*om, om]) if nombre_frequences == 3
                  else np.linspace(0., om, nombre_frequences))
    a = ChaineAssemblee(n, segments, om)
    facteurs = [id(r.qr) for r in a.structures]
    f = np.eye(segments)[-1]
    avant = ControleChamp(a.reduction)
    initial = precision(a, frequences, [avant.reponse(w, f) for w in frequences])
    resultat = adapter_champ(a.reduction, frequences, f, tolerance_relative=1e-10, max_etapes=3)
    final = precision(a, frequences, resultat["reponses"])
    return dict(n=n, segments=segments, nombre_frequences=nombre_frequences,
                frequences=frequences.tolist(), tolerance_relative=1e-10,
                omega_max=float(om), force_physique=1/a.echelle,
                omega_premier_mode_global=float(2*np.sin(np.pi/(4*n+2))),
                historique=resultat["historique"], initial=initial, final=final,
                statut=resultat["statut"], facteurs_reutilises=facteurs == [id(r.qr) for r in a.structures],
                certification_machine=False)


def campagne(destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    (destination/"sources").mkdir()
    sources = {}
    for nom in SOURCES:
        shutil.copyfile(Path(__file__).with_name(nom), destination/"sources"/nom)
        sources[nom] = sha(destination/"sources"/nom)
    essais = []
    for n, s, nf in CAS:
        nom = f"chaine-{n}-{s}-{nf}.json"
        resultat = experience(n, s, nf)
        ecrire(destination/nom, resultat)
        essais.append(dict(fichier=nom, sha256=sha(destination/nom), cas=[n, s, nf]))
        print(nom, resultat["statut"], flush=True)
    ecrire(destination/"manifest.json", dict(date=datetime.now(timezone.utc).isoformat(),
            python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
            essais=essais, sources=sources, chronometrage=False,
            oracle="solution analytique de la chaîne discrète, force physique 1/echelle"))
    verifier(destination)


def verifier(destination):
    d = Path(destination)
    m = json.loads((d/"manifest.json").read_text())
    assert m["chronometrage"] is False
    for nom, h in m["sources"].items():
        assert sha(d/"sources"/nom) == h
    assert [e["cas"] for e in m["essais"]] == [list(c) for c in CAS]
    for e in m["essais"]:
        assert sha(d/e["fichier"]) == e["sha256"]
        r = json.loads((d/e["fichier"]).read_text())
        assert r["certification_machine"] is False and r["facteurs_reutilises"]
        assert [r[k] for k in ("n", "segments", "nombre_frequences")] == e["cas"]
        for stade in ("initial", "final"):
            assert len(r[stade]) == r["nombre_frequences"]
            for point, omega in zip(r[stade], r["frequences"]):
                assert point["omega"] == omega
                for norme in ("masse", "deformation"):
                    erreur = point[f"erreur_{norme}_relative"]
                    borne = point[f"borne_{norme}_relative"]
                    assert np.isfinite(erreur) and erreur >= 0
                    if borne is not None:
                        assert np.isfinite(borne) and erreur <= borne+1e-13
        accepte = all(p[f"borne_{norme}_relative"] is not None
                      and p[f"borne_{norme}_relative"] <= r["tolerance_relative"]
                      for p in r["final"] for norme in ("masse", "deformation"))
        assert accepte == (r["statut"] == "tolerance_aux_frequences_demandees")
        assert accepte == r["historique"][-1]["accepte"]
        assert accepte == (r["nombre_frequences"] == 3)
    print("Archive champ intérieur : 5 expériences cohérentes, dont 2 refus de tolérance.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    choix = parser.add_mutually_exclusive_group(required=True)
    choix.add_argument("--destination", type=Path)
    choix.add_argument("--verifier", type=Path)
    args = parser.parse_args()
    campagne(args.destination) if args.destination else verifier(args.verifier)
