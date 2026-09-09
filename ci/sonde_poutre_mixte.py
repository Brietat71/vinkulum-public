"""Sonde l'arc exact et la flexibilité linéaire du prototype de poutre mixte.

Ces exécutions sont des diagnostics physiques, sans classement en temps.
Les sorties et les refus sont conservés, y compris la limite sans flexion.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def empreinte(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(command):
    process = subprocess.Popen(
        [str(arg) for arg in command],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    timeout = False
    try:
        stdout, stderr = process.communicate(timeout=60)
    except subprocess.TimeoutExpired:
        timeout = True
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
    return {
        "commande": [str(arg) for arg in command],
        "code": process.returncode,
        "delai_depasse": timeout,
        "stdout": stdout.decode(),
        "stderr": stderr.decode(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--sondes", type=Path, required=True)
    parser.add_argument("--sortie", type=Path, required=True)
    args = parser.parse_args()
    executable, sondes = args.executable.resolve(), args.sondes.resolve()
    sources = [
        Path(__file__).resolve().relative_to(ROOT),
        Path("examples/poutre_mixte_sondes.rs"),
        Path("examples/poutre_mixte/element.rs"),
        Path("src/ad.rs"),
    ]
    report = {
        "definition": __doc__,
        "executable_sha256": empreinte(executable),
        "sondes_sha256": empreinte(sondes),
        "sources_sha256": {str(path): empreinte(ROOT / path) for path in sources},
        "moment": [],
        "linearisation": execute([sondes]),
    }
    s = np.linspace(0.0, 1.0, 961)
    exact = np.column_stack((np.sin(s), 1.0 - np.cos(s), np.zeros_like(s)))
    for ordre in (1, 2):
        for elements in (1, 2, 4, 8):
            row = execute([executable, ordre, elements, "moment", 50, "profil"])
            row.update(ordre=ordre, elements=elements)
            if row["code"] == 0:
                data = json.loads(row["stdout"])
                line = np.asarray(data["positions"])
                if line.shape != exact.shape or not np.isfinite(line).all():
                    raise ValueError("profil du moment pur invalide")
                row.update(
                    erreur_bout_m=float(np.max(np.abs(line[-1] - exact[-1]))),
                    erreur_ligne_m=float(np.max(np.abs(line - exact))),
                )
            report["moment"].append(row)
            print(ordre, elements, row.get("erreur_ligne_m", row["stderr"]), flush=True)
    args.sortie.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
