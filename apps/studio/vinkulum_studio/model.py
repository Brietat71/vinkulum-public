"""G0 data contracts. No Qt, mutable solver state, or executable expressions."""

import csv
import io
import json
import math
import os
from dataclasses import asdict, dataclass, fields
from pathlib import Path
import tempfile

MAX_STEPS = 20_000
MAX_RESULT_BYTES = 16 * 1024 * 1024
G = 9.80665
INERTIA = 1e-8
RHO = 0.9
UNITS = {"time": "s", "x": "m", "z": "m", "angle": "rad"}
PARAMETER_UNITS = {"length": "m", "mass": "kg", "angle_deg": "deg",
                   "duration": "s", "step": "s"}


def finite_number(value):
    if type(value) not in (int, float):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


@dataclass(frozen=True)
class Parameters:
    length: float = 1.0
    mass: float = 0.2
    angle_deg: float = 30.0
    duration: float = 2.0
    step: float = 0.005

    def __post_init__(self):
        bounds = {"length": (0.05, 20), "mass": (0.01, 100),
                  "angle_deg": (-170, 170), "duration": (0.001, 600),
                  "step": (0.000001, 1)}
        for name, (low, high) in bounds.items():
            value = getattr(self, name)
            if not finite_number(value):
                raise ValueError(f"{name} : nombre fini requis.")
            if not low <= value <= high:
                raise ValueError(f"{name} : valeur attendue entre {low} et {high}.")
        if self.step > self.duration:
            raise ValueError("Le pas doit être inférieur ou égal à la durée.")
        if self.duration / self.step > MAX_STEPS:
            raise ValueError(f"Limite G0 : {MAX_STEPS:,} pas par calcul. Augmentez le pas.")

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict) or set(data) != {f.name for f in fields(cls)}:
            raise ValueError("Liste de paramètres absente, inconnue ou incomplète.")
        return cls(**data)


def _unique_object(pairs):
    data = {}
    for key, value in pairs:
        if key in data:
            raise ValueError(f"Clé JSON dupliquée : {key}")
        data[key] = value
    return data


def read_json(path, max_bytes=16_384):
    with open(path, "rb") as stream:
        raw = stream.read(max_bytes + 1)
    if len(raw) > max_bytes:
        raise ValueError("Fichier trop volumineux.")
    def invalid(value):
        raise ValueError(f"Nombre JSON interdit : {value}")
    try:
        return json.loads(raw, object_pairs_hook=_unique_object, parse_constant=invalid)
    except (UnicodeError, json.JSONDecodeError, RecursionError) as exc:
        raise ValueError("Fichier JSON illisible.") from exc


def atomic_text(path, content):
    """Replace only after a complete write; leave the previous file on failure."""
    target = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=target.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def write_json(path, data):
    atomic_text(path, json.dumps(data, ensure_ascii=False, allow_nan=False, indent=2) + "\n")


def save_parameters(path, parameters):
    write_json(path, {"format": "vinkulum-studio-parameters", "schema_version": 1,
                      "parameters": asdict(parameters)})


def load_parameters(path):
    data = read_json(path)
    if (not isinstance(data, dict)
            or set(data) != {"format", "schema_version", "parameters"}
            or data["format"] != "vinkulum-studio-parameters"
            or type(data["schema_version"]) is not int or data["schema_version"] != 1):
        raise ValueError("Format ou version de paramètres non pris en charge (version 1 requise).")
    return Parameters.from_dict(data["parameters"])


@dataclass(frozen=True)
class Result:
    run_id: str
    parameters: Parameters
    samples: tuple[tuple[float, float, float, float], ...]
    manifest_json: str

    @classmethod
    def from_dict(cls, data, run_id, parameters):
        if (not isinstance(data, dict) or type(data.get("schema_version")) is not int
                or data["schema_version"] != 1):
            raise ValueError("Version de résultat incorrecte.")
        if data.get("run_id") != run_id or Parameters.from_dict(data.get("parameters")) != parameters:
            raise ValueError("Le résultat ne correspond pas aux entrées capturées.")
        if data.get("status") != "completed":
            raise ValueError("Le worker n'a pas terminé le calcul.")
        rows = data.get("samples")
        if not isinstance(rows, list) or not 2 <= len(rows) <= MAX_STEPS + 2:
            raise ValueError("Nombre d'échantillons incorrect.")
        previous = -1.0
        for row in rows:
            if (not isinstance(row, list) or len(row) != 4
                    or any(not finite_number(v) for v in row)):
                raise ValueError("Échantillon non fini ou mal formé.")
            t, x, z, angle = row
            if t <= previous or t < 0 or t > parameters.duration + 1e-9:
                raise ValueError("Temps des échantillons incohérents.")
            if abs(math.hypot(x, z) - parameters.length) > 1e-5 * max(1, parameters.length):
                raise ValueError("Géométrie de pendule incohérente.")
            if abs(math.atan2(math.sin(angle - math.atan2(x, -z)),
                              math.cos(angle - math.atan2(x, -z)))) > 1e-8:
                raise ValueError("Angle et position incohérents.")
            previous = t
        if rows[0][0] != 0 or not math.isclose(rows[-1][0], parameters.duration, abs_tol=1e-9):
            raise ValueError("Trajectoire incomplète.")
        initial = math.radians(parameters.angle_deg)
        if abs(rows[0][3] - initial) > 1e-12:
            raise ValueError("État initial incohérent.")
        manifest = data.get("manifest")
        if (not isinstance(manifest, dict) or manifest.get("units") != UNITS
                or manifest.get("parameter_units") != PARAMETER_UNITS
                or manifest.get("scientific_status") != "NotAssessed"):
            raise ValueError("Manifeste de provenance absent ou incorrect.")
        return cls(run_id, parameters, tuple(tuple(r) for r in rows),
                   json.dumps(manifest, ensure_ascii=False, allow_nan=False))

    def export_csv(self, path):
        buffer = io.StringIO(newline="")
        metadata = {"format": "vinkulum-studio-samples", "schema_version": 1,
                    "run_id": self.run_id, "parameters": asdict(self.parameters),
                    "manifest": json.loads(self.manifest_json)}
        buffer.write("# " + json.dumps(metadata, ensure_ascii=False, allow_nan=False) + "\n")
        writer = csv.writer(buffer)
        writer.writerow(["time [s]", "x [m]", "z [m]", "angle [rad]"])
        writer.writerows(self.samples)
        atomic_text(path, buffer.getvalue())
