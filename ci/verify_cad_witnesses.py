"""Adversarial replay checks on a retained CAD-witness mesh capture."""

import argparse
import json
import shutil
import tempfile
from dataclasses import asdict
from pathlib import Path

from vinkulum_studio.meshing import CadWitnessMeshedSolid, load_mesh


def verify(mesh, legacy):
    solid, _ = load_mesh(mesh)
    if not isinstance(solid, CadWitnessMeshedSolid) or len(solid.cad_faces) < 2:
        raise ValueError("Use a capture with at least two CAD witness faces.")
    files = [f"cad-face-{identifier}.brep" for identifier, _ in solid.cad_faces]
    checks = ["unaltered witness capture replays"]
    for mode in ("missing", "modified", "swapped", "extra", "missing_payload"):
        with tempfile.TemporaryDirectory(prefix="vinkulum-witness-test-") as temporary:
            copy = Path(temporary) / "mesh"
            shutil.copytree(mesh, copy)
            first, second = copy / files[0], copy / files[1]
            if mode == "missing":
                first.unlink()
            elif mode == "modified":
                first.write_bytes(first.read_bytes() + b"\n")
            elif mode == "swapped":
                a, b = first.read_bytes(), second.read_bytes()
                first.write_bytes(b)
                second.write_bytes(a)
            elif mode == "extra":
                (copy / "cad-face-999999.brep").write_bytes(first.read_bytes())
            else:
                payload = json.loads((copy / "result.json").read_text())
                del payload["result"]["cad_faces"]
                (copy / "result.json").write_text(json.dumps(payload))
            try:
                load_mesh(copy)
            except ValueError:
                checks.append(mode + " witness refused")
            else:
                raise AssertionError(mode + " witness was admitted")
    if legacy is not None:
        old, _ = load_mesh(legacy)
        payload = json.loads((legacy / "result.json").read_text())["result"]
        if json.loads(json.dumps(asdict(old))) != payload:
            raise AssertionError("Legacy mesh payload changed on replay")
        checks.append("legacy payload replays without added fields")
    return {"status": "passed", "checks": checks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--legacy-mesh", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = verify(args.mesh, args.legacy_mesh)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
