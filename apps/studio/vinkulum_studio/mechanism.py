"""Captured-project adapter and versioned scientific result, without Qt."""

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import zipfile

import numpy as np

from . import __version__
from .document import MAX_BODY_SAMPLES, Project
from .model import atomic_text

MAX_ARCHIVE_BYTES = 64 * 1024 * 1024
ARRAY_NAMES = {
    "time",
    "position",
    "rotation",
    "velocity",
    "angular_velocity",
    "joint_coordinate",
}


def build_solver(project):
    from vinkulum import Noyau

    issues = project.diagnostics()
    if issues:
        raise ValueError("\n".join(d.message for d in issues))
    solver = Noyau(g=list(project.gravity))
    indices = {}
    for body in project.bodies:
        indices[body.id] = solver.corps(
            body.name,
            body.mass,
            list(body.inertia()),
            list(body.position),
            rot=list(body.orientation),
        )
    for joint in project.joints:
        bt = [0, 1] if joint.kind == "glissiere" else [0, 1, 2]
        br = (
            []
            if joint.kind == "rotule"
            else [0, 1]
            if joint.kind == "pivot"
            else [0, 1, 2]
        )
        options = {}
        if joint.motion is not None:
            target = ([0.0, 0.0, 1.0], joint.motion.native())
            if joint.kind == "pivot":
                br = [0, 1, 2]
                options["cible_r"] = target
            else:
                bt = [0, 1, 2]
                options["cible_t"] = target
        solver.liaison_reperes(
            joint.name,
            indices.get(joint.a),
            indices.get(joint.b),
            list(joint.pa),
            list(joint.ra),
            list(joint.pb),
            list(joint.rb),
            bloque_t=bt,
            bloque_r=br,
            **options,
        )
    for load in project.loads:
        solver.effort_temporel(
            load.name,
            indices[load.body],
            [law.native() for law in load.force],
            [law.native() for law in load.moment],
            point=list(load.point),
        )
    return solver


def joint_coordinates(project, positions, rotations):
    rows = positions.shape[0]
    indices = {b.id: i for i, b in enumerate(project.bodies)}
    joints = [j for j in project.joints if j.kind in {"pivot", "glissiere"}]
    coordinates = np.empty((rows, len(joints)))

    def frame(ref, point, orientation):
        if ref is None:
            R = np.broadcast_to(np.eye(3), (rows, 3, 3))
            r = np.zeros((rows, 3))
        else:
            index = indices[ref]
            r, R = positions[:, index], rotations[:, index].reshape(-1, 3, 3)
        return r + np.einsum("tij,j->ti", R, point), R @ np.array(orientation).reshape(
            3, 3
        )

    for column, joint in enumerate(joints):
        pa, A = frame(joint.a, joint.pa, joint.ra)
        pb, B = frame(joint.b, joint.pb, joint.rb)
        if joint.kind == "glissiere":
            coordinates[:, column] = np.einsum("ti,ti->t", A[:, :, 2], pb - pa)
        else:
            relative = A.transpose(0, 2, 1) @ B
            coordinates[:, column] = np.arctan2(relative[:, 1, 0], relative[:, 0, 0])
    return coordinates


def simulate_project(project, run_id, directory):
    import vinkulum
    import vinkulum._vinkulum as native

    solver = build_solver(project)
    before = solver.etat()
    if project.joints:
        solver.assemble(t=0.0, tol=1e-12, iters=60, vitesses=True)
    initial = solver.etat()
    if not np.allclose(before[1], initial[1], rtol=0.0, atol=1e-10) or not np.allclose(
        before[2], initial[2], rtol=0.0, atol=1e-10
    ):
        raise ValueError(
            "L'initialisation déplacerait les corps : corrigez les ancrages du document."
        )
    trajectory = [initial, *solver.simule(project.duration, project.step, rho=0.9)]
    arrays = {
        "time": np.array([r[0] for r in trajectory], dtype=np.float64),
        "position": np.array([r[1] for r in trajectory], dtype=np.float64),
        "rotation": np.array([r[2] for r in trajectory], dtype=np.float64),
        "velocity": np.array([r[3] for r in trajectory], dtype=np.float64),
        "angular_velocity": np.array([r[4] for r in trajectory], dtype=np.float64),
    }
    arrays["joint_coordinate"] = joint_coordinates(
        project, arrays["position"], arrays["rotation"]
    )
    # No compression: the on-disk budget also bounds uncompressed arrays.
    path = Path(directory) / "trajectory.npz"
    np.savez(path, **arrays)
    if path.stat().st_size > MAX_ARCHIVE_BYTES:
        raise ValueError("Budget du résultat dépassé.")
    metadata = {
        "schema_version": 2,
        "status": "completed",
        "kind": "mechanism",
        "run_id": run_id,
        "project": asdict(project),
        "archive_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "body_ids": [b.id for b in project.bodies],
        "manifest": {
            "app_version": __version__,
            "kernel_version": vinkulum.__version__,
            "kernel_sha256": hashlib.sha256(
                Path(native.__file__).read_bytes()
            ).hexdigest(),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "method": "generalized-alpha / SO(3)",
            "rho_infinity": 0.9,
            "units": {
                "time": "s",
                "position": "m",
                "rotation": "1",
                "velocity": "m/s",
                "angular_velocity": "rad/s",
            },
            "frame": "world; row-major body-to-world rotation matrices",
            "initial_velocities": arrays["velocity"][0].tolist(),
            "initial_angular_velocities": arrays["angular_velocity"][0].tolist(),
            "joint_angles": "principal angles [-pi,pi]; no turn-count inference",
            "scientific_status": "NotAssessed",
            "sampling": "native steps plus initialized state",
            "limitations": "No trajectory error bound or checkpoint; rigid mechanisms only.",
        },
    }
    # Validate before publishing metadata; the parent independently repeats checks.
    MechanicalResult.read(metadata, run_id, project, directory)
    return metadata


def immutable(array):
    array = np.ascontiguousarray(array, dtype=np.float64)
    return np.frombuffer(array.tobytes(), dtype=np.float64).reshape(array.shape)


@dataclass(frozen=True)
class MechanicalResult:
    run_id: str
    project: Project
    time: np.ndarray
    position: np.ndarray
    rotation: np.ndarray
    velocity: np.ndarray
    angular_velocity: np.ndarray
    joint_coordinate: np.ndarray
    manifest_json: str

    @classmethod
    def read(cls, data, run_id, project, directory):
        if (
            not isinstance(data, dict)
            or type(data.get("schema_version")) is not int
            or data["schema_version"] != 2
            or data.get("kind") != "mechanism"
            or data.get("status") != "completed"
            or data.get("run_id") != run_id
            or Project.from_dict(data.get("project")) != project
            or data.get("body_ids") != [b.id for b in project.bodies]
        ):
            raise ValueError("Résultat incompatible avec le projet capturé.")
        manifest = data.get("manifest")
        if (
            not isinstance(manifest, dict)
            or manifest.get("scientific_status") != "NotAssessed"
        ):
            raise ValueError("Manifeste scientifique incorrect.")
        path = Path(directory) / "trajectory.npz"
        if path.stat().st_size > MAX_ARCHIVE_BYTES:
            raise ValueError("Archive de résultat trop volumineuse.")
        if hashlib.sha256(path.read_bytes()).hexdigest() != data.get("archive_sha256"):
            raise ValueError("Empreinte de résultat incorrecte.")
        with zipfile.ZipFile(path) as archive:
            infos = archive.infolist()
            if (
                len(infos) != len(ARRAY_NAMES)
                or {i.filename for i in infos} != {n + ".npy" for n in ARRAY_NAMES}
                or sum(i.file_size for i in infos) > MAX_ARCHIVE_BYTES
            ):
                raise ValueError(
                    "Contenu ou taille décompressée du résultat incorrects."
                )
            for info in infos:
                with archive.open(info) as stream:
                    version = np.lib.format.read_magic(stream)
                    if version == (1, 0):
                        shape, _, dtype = np.lib.format.read_array_header_1_0(stream)
                    elif version == (2, 0):
                        shape, _, dtype = np.lib.format.read_array_header_2_0(stream)
                    else:
                        raise ValueError("Format de tableau non pris en charge.")
                    size = math.prod(shape) * dtype.itemsize
                    if (
                        dtype != np.dtype("float64")
                        or size != info.file_size - stream.tell()
                        or size > MAX_ARCHIVE_BYTES
                    ):
                        raise ValueError(
                            "Taille ou type déclaré du tableau incohérent."
                        )
        with np.load(path, allow_pickle=False) as source:
            arrays = {name: source[name] for name in ARRAY_NAMES}
        times = arrays["time"]
        bodies = len(project.bodies)
        if times.ndim != 1 or len(times) < 2 or len(times) * bodies > MAX_BODY_SAMPLES:
            raise ValueError("Nombre d'échantillons incorrect.")
        for name, array in arrays.items():
            if array.dtype != np.float64 or not np.isfinite(array).all():
                raise ValueError("Données scientifiques non finies ou type incorrect.")
            if name in {"position", "rotation", "velocity", "angular_velocity"}:
                if array.shape != (len(times), bodies, 9 if name == "rotation" else 3):
                    raise ValueError("Dimensions des données incorrectes.")
        if (
            times[0] != 0
            or not np.all(np.diff(times) > 0)
            or not math.isclose(times[-1], project.duration, abs_tol=1e-9)
        ):
            raise ValueError("Trajectoire incomplète ou temps incohérents.")
        rotations = arrays["rotation"].reshape(-1, 3, 3)
        if (
            np.max(np.abs(rotations.transpose(0, 2, 1) @ rotations - np.eye(3))) > 1e-8
            or np.max(np.abs(np.linalg.det(rotations) - 1)) > 1e-8
        ):
            raise ValueError("Rotations hors SO(3).")
        if not np.allclose(
            arrays["position"][0],
            [b.position for b in project.bodies],
            rtol=0.0,
            atol=1e-10,
        ) or not np.allclose(
            arrays["rotation"][0],
            [b.orientation for b in project.bodies],
            rtol=0.0,
            atol=1e-10,
        ):
            raise ValueError("État initial différent du projet capturé.")
        expected = joint_coordinates(project, arrays["position"], arrays["rotation"])
        if arrays["joint_coordinate"].shape != expected.shape or not np.allclose(
            arrays["joint_coordinate"], expected, atol=1e-10, rtol=0.0
        ):
            raise ValueError("Coordonnées de liaison incohérentes.")
        return cls(
            run_id,
            project,
            **{k: immutable(v) for k, v in arrays.items()},
            manifest_json=json.dumps(manifest, ensure_ascii=False, allow_nan=False),
        )

    def poses(self, index):
        return {
            body.id: (self.position[index, i], self.rotation[index, i])
            for i, body in enumerate(self.project.bodies)
        }

    def series(self):
        """Label, unit, immutable sample array; all labels identify their object."""
        output = []
        for index, body in enumerate(self.project.bodies):
            for axis, name in enumerate("XYZ"):
                output.append(
                    (
                        f"{body.name} [{body.id[:8]}] · position {name}",
                        "m",
                        self.position[:, index, axis],
                    )
                )
                output.append(
                    (
                        f"{body.name} [{body.id[:8]}] · vitesse {name}",
                        "m/s",
                        self.velocity[:, index, axis],
                    )
                )
        index = 0
        for joint in self.project.joints:
            if joint.kind in {"pivot", "glissiere"}:
                suffix = "angle principal" if joint.kind == "pivot" else "déplacement"
                output.append(
                    (
                        f"{joint.name} [{joint.id[:8]}] · {suffix}",
                        "rad" if joint.kind == "pivot" else "m",
                        self.joint_coordinate[:, index],
                    )
                )
                index += 1
        return output

    def export_csv(self, path):
        buffer = io.StringIO(newline="")
        buffer.write(
            "# "
            + json.dumps(
                {
                    "format": "vinkulum-studio-mechanism-samples",
                    "schema_version": 1,
                    "run_id": self.run_id,
                    "project": asdict(self.project),
                    "manifest": json.loads(self.manifest_json),
                },
                ensure_ascii=False,
                allow_nan=False,
            )
            + "\n"
        )
        writer = csv.writer(buffer)
        header = ["time [s]"]
        for body in self.project.bodies:
            for prefix, unit, components in (
                ("position", "m", "xyz"),
                ("rotation", "1", range(9)),
                ("velocity", "m/s", "xyz"),
                ("angular_velocity", "rad/s", "xyz"),
            ):
                header.extend(f"{body.id}.{prefix}.{c} [{unit}]" for c in components)
        joints = [j for j in self.project.joints if j.kind in {"pivot", "glissiere"}]
        header.extend(
            f"{j.id}.coordinate [{'rad principal' if j.kind == 'pivot' else 'm'}]"
            for j in joints
        )
        writer.writerow(header)
        for i, t in enumerate(self.time):
            row = [t]
            for body in range(len(self.project.bodies)):
                for array in (
                    self.position,
                    self.rotation,
                    self.velocity,
                    self.angular_velocity,
                ):
                    row.extend(array[i, body])
            row.extend(self.joint_coordinate[i])
            writer.writerow(row)
        atomic_text(path, buffer.getvalue())
