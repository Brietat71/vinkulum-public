"""Archive et vérifie les mesures du prototype mixte, sans exécutable tiers.

Le vérificateur recalcule les erreurs à partir des sorties brutes. Il contrôle
les empreintes des programmes mesurés, les répétitions et les équilibres.
Une archive cohérente peut contenir des échecs : ils restent dans le bilan.
"""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import statistics
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def exige(condition, message):
    if not condition:
        raise ValueError(message)


def proche(a, b, message, tolerance=1e-16):
    exige(np.allclose(a, b, rtol=1e-13, atol=tolerance), message)


def profil(data):
    line = np.asarray(data["positions"])
    exige(line.shape == (961, 3) and np.isfinite(line).all(), "profil invalide")
    proche(line[-1], data["position"], "bout différent du profil", 1e-14)
    for key, limit in (("residu", 1e-8), ("residu_max_paliers", 1e-8), ("residu_interne", 1e-9)):
        exige(0 <= data[key] < limit, f"équilibre insuffisant : {key}")
    exige(data["paliers"] == 50, "nombre de paliers modifié")
    return line


def verifie(archive):
    exige(archive["schema"] == 1, "schéma inconnu")
    report = json.loads(archive["campagne_json"])
    reference = json.loads(archive["reference_json"])
    diagnostics = json.loads(archive["diagnostics_json"])
    sources = archive["sources"]
    for metadata in (report["metadata"], diagnostics):
        for name, expected in metadata["sources_sha256"].items():
            exige(sha(sources[name].encode()) == expected, f"source mesurée différente : {name}")
    exige(
        sha(archive["reference_json"].encode()) == report["metadata"]["reference_sha256"],
        "référence différente de celle de la campagne",
    )
    exige(
        sha(sources["ci/reference_cosserat.py"].encode()) == reference["programme_sha256"],
        "générateur de référence différent",
    )
    exige(diagnostics["executable_sha256"] == report["metadata"]["executable_sha256"], "binaires différents")
    refs = reference["references"]
    exige(len(refs) == 2, "deux références continues requises")
    force = np.array([0.0, 8.896 / np.sqrt(2), 8.896 / np.sqrt(2)])
    for ref in refs:
        line = np.asarray(ref["positions"])
        rotations = np.asarray(ref["rotations"])
        moments = np.asarray(ref["moments"])
        exige(line.shape == (961, 3) and rotations.shape == (961, 3, 3), "dimensions de référence")
        exige(np.isfinite(line).all() and np.isfinite(rotations).all(), "référence non finie")
        proche(ref["coordonnee"], np.linspace(0.0, 0.508, 961), "grille continue différente")
        orthogonal = np.max(np.abs(np.einsum("nji,njk->nik", rotations, rotations) - np.eye(3)))
        balance = np.max(np.abs(moments + np.cross(line, force) - ref["moment_racine"]))
        proche(orthogonal, ref["defaut_orthogonalite"], "défaut de rotation incohérent")
        proche(balance, ref["defaut_bilan_moment_nm"], "bilan continu incohérent")
        exige(orthogonal < 1e-8 and balance < 1e-9, "invariant continu insuffisant")
        proche(line[0], np.zeros(3), "racine déplacée")
        proche(rotations[0], np.eye(3), "racine tournée")
        exige(np.max(np.abs(moments[-1])) < 1e-10, "moment terminal non nul")
        exige(len(ref["paliers"]) == 16, "continuation continue incomplète")
        exige(max(row["residu_moment_nm"] for row in ref["paliers"]) < 1e-10, "tir non convergé")
    change = float(np.max(np.abs(np.asarray(refs[0]["positions"]) - refs[1]["positions"])))
    proche(change, reference["ecart_raffinement_m"], "raffinement continu incohérent")
    exige(change < 1e-9, "référence insuffisamment raffinée")
    fine = np.asarray(refs[-1]["positions"])
    expected = {(k, n, r) for k, ns in ((1, (2, 4, 8, 16)), (2, (1, 2, 3, 4, 8))) for n in ns for r in (-1, 0, 1, 2)}
    rows = report["essais"]
    exige(len(rows) == len(expected), "essai omis ou ajouté")
    exige({(s["ordre"], s["elements"], s["repetition"]) for s in rows} == expected, "plan d'essais différent")
    for row in rows:
        key = f'ordre-{row["ordre"]}-elements-{row["elements"]}-essai-{row["repetition"]}'
        raw = archive["sorties"][key]
        for channel in ("stdout", "stderr"):
            exige(sha(raw[channel].encode()) == row[f"{channel}_sha256"], f"sortie différente : {key}/{channel}")
        if not row["ok"]:
            continue
        exige(row["code"] == 0 and not row["delai_depasse"], "succès contredit par le processus")
        data = json.loads(raw["stdout"])
        exige(data == row["resultat"], "résultat différent du JSON brut")
        exige(data["cas"] == "princeton" and data["ordre"] == row["ordre"] and data["elements"] == row["elements"], "cas incohérent")
        line = profil(data)
        proche(np.max(np.abs(line[-1] - fine[-1])), row["erreur_bout_m"], "erreur au bout incohérente")
        proche(np.max(np.abs(line - fine)), row["erreur_ligne_m"], "erreur sur la ligne incohérente")
        reaction = np.asarray(data["reaction_racine"])
        proche(reaction[:3] + force, np.zeros(3), "bilan global des forces", 1e-8)
        proche(reaction[3:] + np.cross(line[-1], force), np.zeros(3), "bilan global des moments", 2e-8)
        proche(np.max(np.abs(reaction[3:] + refs[-1]["moment_racine"])), row["erreur_reaction_moment_nm"], "erreur de réaction incohérente")
        exige(row["temps_total_s"] > data["temps_s"] > 0, "temps incohérents")
        exige(int(raw["rss"].splitlines()[-1]) == row["rss_kib"], "RSS différent")
    configs = report["configurations"]
    exige(len(configs) == 9 and len({(c["ordre"], c["elements"]) for c in configs}) == 9, "configurations incomplètes")
    for config in configs:
        samples = [s for s in rows if s["repetition"] >= 0 and (s["ordre"], s["elements"]) == (config["ordre"], config["elements"])]
        good = [s for s in samples if s["ok"]]
        exige(config["essais"] == len(samples) == 3 and config["reussites"] == len(good), "répétitions incohérentes")
        if not good:
            continue
        proche(statistics.median(s["temps_total_s"] for s in good), config["temps_median_s"], "temps médian incohérent")
        proche(statistics.median(s["resultat"]["temps_s"] for s in good), config["temps_calcul_median_s"], "temps de calcul incohérent")
        for key in ("erreur_bout_m", "erreur_ligne_m"):
            proche(max(s[key] for s in good), config[key], f"maximum incohérent : {key}")
        exige(max(s["rss_kib"] for s in good) == config["rss_max_kib"], "RSS maximal incohérent")
        delta = max(float(np.max(np.abs(np.asarray(s["resultat"]["positions"]) - good[0]["resultat"]["positions"]))) for s in good)
        proche(delta, config["ecart_repetitions_m"], "répétabilité incohérente")
    s = np.linspace(0.0, 1.0, 961)
    exact = np.column_stack((np.sin(s), 1.0 - np.cos(s), np.zeros_like(s)))
    moment = diagnostics["moment"]
    exige(len(moment) == 8 and {(s["ordre"], s["elements"]) for s in moment} == {(k, n) for k in (1, 2) for n in (1, 2, 4, 8)}, "sondages d'arc incomplets")
    for row in moment:
        if row["code"] != 0:
            continue
        data = json.loads(row["stdout"])
        line = profil(data)
        proche(np.max(np.abs(line[-1] - exact[-1])), row["erreur_bout_m"], "erreur d'arc au bout incohérente")
        proche(np.max(np.abs(line - exact)), row["erreur_ligne_m"], "erreur d'arc incohérente")
        proche(data["angle_z"], 1.0, "angle final de l'arc", 2e-10)
        proche(data["reaction_racine"], [0, 0, 0, 0, 0, -100], "réaction du moment pur", 2e-8)
    linear = diagnostics["linearisation"]
    exige(linear["code"] == 0, "sonde linéaire non exécutée")
    probes = [json.loads(line) for line in linear["stdout"].splitlines()]
    exige(len(probes) == 20, "sondes linéaires incomplètes")
    exige("8 passed; 0 failed" in archive["tests"], "journal des contrôles physiques absent")
    return {
        "essais_princeton": len(rows),
        "succes_princeton": sum(s["ok"] for s in rows),
        "echauffements": sum(s["repetition"] < 0 for s in rows),
        "sondages_moment": len(moment),
        "sondages_lineaires": len(probes),
        "refus_condensation": sum(not p["condensation"] for p in probes),
        "raideurs_negatives": sum(p.get("valeur_propre_min", 0) < 0 for p in probes),
        "ecart_reference_m": change,
        "configurations": configs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--verifier", type=Path, help="préfixe des deux fichiers d'archive")
    group.add_argument("--creer", type=Path, help="préfixe de sortie, sans écrasement")
    parser.add_argument("--campagne", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--diagnostics", type=Path)
    parser.add_argument("--tests", type=Path)
    args = parser.parse_args()
    prefix = args.verifier or args.creer
    payload_path = Path(str(prefix) + "-donnees.json.gz")
    manifest_path = Path(str(prefix) + ".json")
    if args.verifier:
        payload = payload_path.read_bytes()
        manifest = json.loads(manifest_path.read_bytes())
        exige(sha(payload) == manifest["archive_sha256"], "empreinte d'archive incorrecte")
        archive = json.loads(gzip.decompress(payload))
        bilan = verifie(archive)
        exige(bilan == manifest["bilan"], "bilan du manifeste différent des mesures")
    else:
        exige(not payload_path.exists() and not manifest_path.exists(), "archive déjà présente")
        exige(all((args.campagne, args.reference, args.diagnostics, args.tests)), "entrées manquantes")
        report = json.loads((args.campagne / "resultats.json").read_bytes())
        diagnostics = json.loads(args.diagnostics.read_bytes())
        paths = set(report["metadata"]["sources_sha256"]) | set(diagnostics["sources_sha256"])
        paths |= {"Cargo.toml", "Cargo.lock", "ci/archive_poutre_mixte.py"}
        sorties = {}
        for row in report["essais"]:
            key = f'ordre-{row["ordre"]}-elements-{row["elements"]}-essai-{row["repetition"]}'
            directory = args.campagne / key
            sorties[key] = {
                "stdout": (directory / "stdout.json").read_text(),
                "stderr": (directory / "stderr.txt").read_text(),
                "rss": (directory / "rss.txt").read_text() if (directory / "rss.txt").exists() else "",
            }
        archive = {
            "schema": 1,
            "campagne_json": (args.campagne / "resultats.json").read_text(),
            "reference_json": args.reference.read_text(),
            "diagnostics_json": args.diagnostics.read_text(),
            "sources": {name: (ROOT / name).read_text() for name in sorted(paths)},
            "sorties": sorties,
            "tests": args.tests.read_text(),
        }
        bilan = verifie(archive)
        payload = gzip.compress((json.dumps(archive, ensure_ascii=False) + "\n").encode(), mtime=0)
        manifest = {
            "description": "Prototype statique autonome d'ordres 1 et 2 ; aucune nouvelle API distribuée.",
            "base_git": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "archive": payload_path.name,
            "archive_sha256": sha(payload),
            "executable_sha256": report["metadata"]["executable_sha256"],
            "sondes_executable_sha256": diagnostics["sondes_sha256"],
            "sources_sha256": {name: sha(text.encode()) for name, text in archive["sources"].items()},
            "bilan": bilan,
        }
        payload_path.write_bytes(payload)
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in bilan.items() if k != "configurations"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
