"""Mechanical authoring data, independent of Qt and mutable solver instances.

All lengths are metres, rotations are row-major body-to-world matrices,
and joint coordinates are expressed in each body's own frame. None is ground.
"""

from dataclasses import asdict, dataclass, replace
import math
import uuid

import numpy as np

from .cad_data import CadGeometry

from .model import Parameters, finite_number, load_parameters, read_json, write_json

IDENTITY = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
ZERO = (0.0, 0.0, 0.0)
MAX_BODIES = 32
MAX_JOINTS = 64
MAX_LOADS = 128
MAX_BODY_SAMPLES = 100_000
MAX_PROJECT_BYTES = 32 * 1024 * 1024


def new_id():
    return str(uuid.uuid4())


def vector(value, size, name):
    if isinstance(value, (tuple, list)):
        value = tuple(v.item() if isinstance(v, np.generic) else v for v in value)
    if (
        not isinstance(value, (tuple, list))
        or len(value) != size
        or not all(finite_number(v) for v in value)
    ):
        raise ValueError(f"{name} : {size} composantes finies requises.")
    return tuple(float(v) for v in value)


def rotation(value):
    values = vector(value, 9, "Rotation")
    matrix = np.array(values).reshape(3, 3)
    if (
        np.linalg.norm(matrix.T @ matrix - np.eye(3)) > 1e-9
        or abs(np.linalg.det(matrix) - 1) > 1e-9
    ):
        raise ValueError("Rotation hors SO(3).")
    return values


def named(identifier, name):
    if not isinstance(identifier, str) or str(uuid.UUID(identifier)) != identifier:
        raise ValueError("Identité UUID invalide.")
    if not isinstance(name, str) or not name.strip() or len(name) > 128:
        raise ValueError("Nom requis, limité à 128 caractères.")


def strict(data, cls):
    from dataclasses import fields

    if not isinstance(data, dict) or set(data) != {f.name for f in fields(cls)}:
        raise ValueError(f"Champs inconnus ou manquants pour {cls.__name__}.")
    return data


@dataclass(frozen=True)
class Law:
    kind: str = "constante"
    values: tuple = (0.0,)

    def __post_init__(self):
        values = self.values
        if not isinstance(values, tuple) or not all(finite_number(v) for v in values):
            raise ValueError("Loi : coefficients finis requis.")
        if self.kind == "constante" and len(values) == 1:
            return
        if self.kind == "lineaire" and len(values) == 2:
            return
        if self.kind == "table" and 2 <= len(values) <= 20_000 and len(values) % 2 == 0:
            if all(a < b for a, b in zip(values[::2], values[2::2])):
                return
        raise ValueError(
            "Loi invalide : constante, linéaire ou table à temps strictement croissants requise."
        )

    def value(self, t):
        if self.kind == "constante":
            return self.values[0]
        if self.kind == "lineaire":
            return self.values[0] + self.values[1] * t
        return float(np.interp(t, self.values[::2], self.values[1::2]))

    def native(self):
        return self.kind, list(self.values)

    @classmethod
    def from_dict(cls, data):
        strict(data, cls)
        return cls(data["kind"], tuple(data["values"]))


@dataclass(frozen=True)
class Body:
    id: str
    name: str
    shape: str = "box"
    dimensions: tuple = (0.2, 0.2, 0.2)
    mass: float = 1.0
    position: tuple = ZERO
    orientation: tuple = IDENTITY
    inertia_mode: str = "homogeneous"
    explicit_inertia: tuple = IDENTITY
    cad: CadGeometry | None = None

    def __post_init__(self):
        named(self.id, self.name)
        if self.shape not in {"box", "cylinder", "sphere", "cad"}:
            raise ValueError("Primitive inconnue.")
        object.__setattr__(
            self,
            "dimensions",
            vector(
                self.dimensions,
                {"box": 3, "cylinder": 2, "sphere": 1, "cad": 3}[self.shape],
                "Dimensions",
            ),
        )
        if any(v <= 0 or v > 1000 for v in self.dimensions):
            raise ValueError("Dimensions strictement positives, au plus 1000 m.")
        if not finite_number(self.mass) or not 0 < self.mass <= 1e9:
            raise ValueError("Masse strictement positive, au plus 10⁹ kg.")
        object.__setattr__(self, "position", vector(self.position, 3, "Position"))
        object.__setattr__(self, "orientation", rotation(self.orientation))
        object.__setattr__(
            self, "explicit_inertia", vector(self.explicit_inertia, 9, "Inertie")
        )
        if self.inertia_mode not in {"homogeneous", "explicit"}:
            raise ValueError("Mode d'inertie inconnu.")
        if (self.shape == "cad" and not isinstance(self.cad, CadGeometry)) or (
            self.shape != "cad" and self.cad is not None
        ):
            raise ValueError("Géométrie CAD incohérente avec le corps.")
        if self.inertia_mode == "explicit":
            j = np.array(self.explicit_inertia).reshape(3, 3)
            if (
                not np.allclose(j, j.T, atol=1e-14, rtol=1e-12)
                or np.linalg.eigvalsh(j).min() <= 0
            ):
                raise ValueError("Inertie symétrique définie positive requise.")
            eigen = np.linalg.eigvalsh(j)
            if eigen[-1] > sum(eigen[:2]) + 1e-12 * eigen[-1]:
                raise ValueError(
                    "Inertie non physique : inégalité triangulaire violée."
                )

    def inertia(self):
        if self.inertia_mode == "explicit":
            return self.explicit_inertia
        if self.cad is not None:
            return tuple(self.mass * v for v in self.cad.unit_inertia_m2)
        if self.shape == "box":
            x, y, z = self.dimensions
            diagonal = (
                self.mass * (y * y + z * z) / 12,
                self.mass * (x * x + z * z) / 12,
                self.mass * (x * x + y * y) / 12,
            )
        elif self.shape == "sphere":
            diagonal = (2 * self.mass * self.dimensions[0] ** 2 / 5,) * 3
        else:
            radius, length = self.dimensions  # cylinder axis is local Z
            transverse = self.mass * (3 * radius**2 + length**2) / 12
            diagonal = (transverse, transverse, self.mass * radius**2 / 2)
        return tuple(float(v) for v in np.diag(diagonal).flat)

    @classmethod
    def from_dict(cls, data):
        # Schema 1 bodies predate the optional BREP payload.
        if isinstance(data, dict) and "cad" not in data:
            data = {**data, "cad": None}
        data = dict(strict(data, cls))
        for key in ("dimensions", "position", "orientation", "explicit_inertia"):
            data[key] = tuple(data[key])
        data["cad"] = (
            None if data["cad"] is None else CadGeometry.from_dict(data["cad"])
        )
        return cls(**data)


@dataclass(frozen=True)
class Joint:
    id: str
    name: str
    kind: str
    a: str | None
    b: str | None
    pa: tuple = ZERO
    ra: tuple = IDENTITY
    pb: tuple = ZERO
    rb: tuple = IDENTITY
    motion: Law | None = None

    def __post_init__(self):
        named(self.id, self.name)
        if self.kind not in {"pivot", "rotule", "glissiere", "encastrement"}:
            raise ValueError("Type de liaison inconnu.")
        if self.a == self.b:
            raise ValueError("La liaison doit relier deux objets distincts.")
        for reference in (self.a, self.b):
            if reference is not None:
                named(reference, "référence")
        object.__setattr__(self, "pa", vector(self.pa, 3, "Ancrage A"))
        object.__setattr__(self, "pb", vector(self.pb, 3, "Ancrage B"))
        object.__setattr__(self, "ra", rotation(self.ra))
        object.__setattr__(self, "rb", rotation(self.rb))
        if self.motion is not None and (
            not isinstance(self.motion, Law) or self.kind not in {"pivot", "glissiere"}
        ):
            raise ValueError("Mouvement imposé réservé au pivot et à la glissière.")

    @classmethod
    def from_dict(cls, data):
        data = dict(strict(data, cls))
        for key in ("pa", "ra", "pb", "rb"):
            data[key] = tuple(data[key])
        data["motion"] = (
            None if data["motion"] is None else Law.from_dict(data["motion"])
        )
        return cls(**data)


@dataclass(frozen=True)
class Load:
    id: str
    name: str
    body: str
    point: tuple = ZERO
    force: tuple = (Law(), Law(), Law())
    moment: tuple = (Law(), Law(), Law())

    def __post_init__(self):
        named(self.id, self.name)
        named(self.body, "référence")
        object.__setattr__(self, "point", vector(self.point, 3, "Point d'application"))
        for laws in (self.force, self.moment):
            if (
                not isinstance(laws, tuple)
                or len(laws) != 3
                or not all(isinstance(v, Law) for v in laws)
            ):
                raise ValueError("Trois lois scalaires requises par vecteur de charge.")

    @classmethod
    def from_dict(cls, data):
        data = dict(strict(data, cls))
        data["point"] = tuple(data["point"])
        for key in ("force", "moment"):
            data[key] = tuple(Law.from_dict(v) for v in data[key])
        return cls(**data)


@dataclass(frozen=True)
class Diagnostic:
    object_id: str
    message: str


@dataclass(frozen=True)
class Project:
    id: str
    name: str = "Mécanisme"
    revision: int = 0
    bodies: tuple = ()
    joints: tuple = ()
    loads: tuple = ()
    deleted: tuple = ()
    gravity: tuple = (0.0, 0.0, -9.80665)
    duration: float = 2.0
    step: float = 0.005

    def __post_init__(self):
        named(self.id, self.name)
        if type(self.revision) is not int or self.revision < 0:
            raise ValueError("Révision entière positive ou nulle requise.")
        object.__setattr__(self, "gravity", vector(self.gravity, 3, "Gravité"))
        Parameters(duration=self.duration, step=self.step)
        for objects, kind, limit in (
            (self.bodies, Body, MAX_BODIES),
            (self.joints, Joint, MAX_JOINTS),
            (self.loads, Load, MAX_LOADS),
        ):
            if (
                not isinstance(objects, tuple)
                or len(objects) > limit
                or not all(isinstance(v, kind) for v in objects)
            ):
                raise ValueError(
                    f"Collection {kind.__name__} invalide ou limite {limit} dépassée."
                )
        ids = [o.id for o in (*self.bodies, *self.joints, *self.loads)]
        cad_parts = [b.cad for b in self.bodies if b.cad is not None]
        if (
            sum(len(c.brep_mm) for c in cad_parts) > 4_000_000
            or sum(len(c.triangles) + len(c.vertices_m) for c in cad_parts) > 100_000
        ):
            raise ValueError(
                "Budget CAD du document dépassé (4 Mo BREP, 100 000 éléments de maillage)."
            )
        if len(ids) != len(set(ids)):
            raise ValueError("Identités dupliquées.")
        if not isinstance(self.deleted, tuple) or len(set(self.deleted)) != len(
            self.deleted
        ):
            raise ValueError("Liste des objets supprimés invalide.")
        for identifier in self.deleted:
            named(identifier, "objet supprimé")
            if identifier in ids:
                raise ValueError("Un objet ne peut être présent et supprimé.")
        if (
            len(self.bodies) * (math.ceil(self.duration / self.step) + 1)
            > MAX_BODY_SAMPLES
        ):
            raise ValueError(
                f"Budget dépassé : {MAX_BODY_SAMPLES} échantillons-corps au maximum."
            )

    def pose(self, identifier):
        if identifier is None:
            return np.zeros(3), np.eye(3)
        body = next(b for b in self.bodies if b.id == identifier)
        return np.array(body.position), np.array(body.orientation).reshape(3, 3)

    def diagnostics(self):
        issues = []
        ids = {body.id for body in self.bodies}
        if not ids:
            issues.append(Diagnostic(self.id, "Ajoutez au moins un corps."))
        for obj in (*self.joints, *self.loads):
            references = (obj.a, obj.b) if isinstance(obj, Joint) else (obj.body,)
            missing = [r for r in references if r is not None and r not in ids]
            for r in missing:
                state = "supprimé" if r in self.deleted else "introuvable"
                issues.append(Diagnostic(obj.id, f"Corps {state} : {r}."))
            if missing or isinstance(obj, Load):
                continue
            a, A = self.pose(obj.a)
            b, B = self.pose(obj.b)
            QA, QB = (
                A @ np.array(obj.ra).reshape(3, 3),
                B @ np.array(obj.rb).reshape(3, 3),
            )
            delta = QA.T @ (b + B @ obj.pb - a - A @ obj.pa)
            limit = 1e-8 * max(1.0, np.linalg.norm(obj.pa), np.linalg.norm(obj.pb))
            transl = delta[:2] if obj.kind == "glissiere" else delta
            if np.linalg.norm(transl) > limit:
                issues.append(
                    Diagnostic(
                        obj.id,
                        "Ancrages incompatibles : corrigez le placement ou les repères.",
                    )
                )
            rotational_error = 0.0
            if obj.kind == "pivot":
                rotational_error = np.linalg.norm(QA[:, 2] - QB[:, 2])
            elif obj.kind in {"glissiere", "encastrement"}:
                rotational_error = np.linalg.norm(QA - QB)
            if rotational_error > 1e-8:
                issues.append(
                    Diagnostic(obj.id, "Axes ou orientations de liaison incompatibles.")
                )
            if obj.motion is not None:
                relative = QA.T @ QB
                coordinate = (
                    delta[2]
                    if obj.kind == "glissiere"
                    else math.atan2(relative[1, 0], relative[0, 0])
                )
                if abs(coordinate - obj.motion.value(0)) > 1e-8:
                    issues.append(
                        Diagnostic(
                            obj.id, "La consigne initiale ne correspond pas à la pose."
                        )
                    )
        return tuple(issues)

    def replace_object(self, obj):
        key = (
            "bodies"
            if isinstance(obj, Body)
            else "joints"
            if isinstance(obj, Joint)
            else "loads"
        )
        objects = getattr(self, key)
        if not any(o.id == obj.id for o in objects):
            objects = (*objects, obj)
        else:
            objects = tuple(obj if old.id == obj.id else old for old in objects)
        return replace(self, **{key: objects})

    def remove(self, identifier):
        if not any(
            o.id == identifier for o in (*self.bodies, *self.joints, *self.loads)
        ):
            raise ValueError("Objet absent.")
        return replace(
            self,
            bodies=tuple(b for b in self.bodies if b.id != identifier),
            joints=tuple(j for j in self.joints if j.id != identifier),
            loads=tuple(l for l in self.loads if l.id != identifier),
            deleted=(*self.deleted, identifier),
        )

    @classmethod
    def from_dict(cls, data):
        data = dict(strict(data, cls))
        for key, kind in (("bodies", Body), ("joints", Joint), ("loads", Load)):
            data[key] = tuple(kind.from_dict(v) for v in data[key])
        data["gravity"] = tuple(data["gravity"])
        data["deleted"] = tuple(data["deleted"])
        return cls(**data)


def joint_at(project, kind, a, b, point=ZERO, axis=(0.0, 0.0, 1.0)):
    point = np.array(vector(point, 3, "Point"))
    z = np.array(vector(axis, 3, "Axe"))
    norm = np.linalg.norm(z)
    if not math.isfinite(norm) or norm < 1e-12:
        raise ValueError("Axe non nul requis.")
    z /= norm
    seed = np.eye(3)[int(np.argmin(abs(z)))]
    x = np.cross(seed, z)
    x /= np.linalg.norm(x)
    frame = np.column_stack((x, np.cross(z, x), z))
    ar, A = project.pose(a)
    br, B = project.pose(b)
    return Joint(
        new_id(),
        kind.capitalize(),
        kind,
        a,
        b,
        tuple(A.T @ (point - ar)),
        tuple((A.T @ frame).flat),
        tuple(B.T @ (point - br)),
        tuple((B.T @ frame).flat),
    )


class History:
    def __init__(self, project):
        self.current = project
        self.past = []
        self.future = []

    def commit(self, project):
        if project == self.current:
            return
        self.past.append(self.current)
        self.past = self.past[-100:]
        self.future.clear()
        self.current = replace(project, revision=self.current.revision + 1)

    def undo(self):
        if self.past:
            self.future.append(self.current)
            self.current = replace(self.past.pop(), revision=self.current.revision + 1)

    def redo(self):
        if self.future:
            self.past.append(self.current)
            self.current = replace(
                self.future.pop(), revision=self.current.revision + 1
            )


def save_project(path, project):
    data = asdict(project)
    for body in data["bodies"]:
        if body["cad"] is None:
            del body["cad"]
    write_json(
        path,
        {
            "format": "vinkulum-studio-project",
            "schema_version": 2
            if any(b.cad is not None for b in project.bodies)
            else 1,
            "project": data,
        },
    )


def load_project(path):
    data = read_json(path, MAX_PROJECT_BYTES)
    if (
        not isinstance(data, dict)
        or set(data) != {"format", "schema_version", "project"}
        or data["format"] != "vinkulum-studio-project"
        or type(data["schema_version"]) is not int
        or data["schema_version"] not in (1, 2)
    ):
        raise ValueError(
            "Format de projet inconnu. Utilisez l'import G0 pour un ancien pendule."
        )
    return Project.from_dict(data["project"])


def pendulum(parameters=Parameters()):
    angle = math.radians(parameters.angle_deg)
    body = Body(
        new_id(),
        "Masse du pendule",
        "sphere",
        (0.06,),
        parameters.mass,
        (
            parameters.length * math.sin(angle),
            0.0,
            -parameters.length * math.cos(angle),
        ),
        inertia_mode="explicit",
        explicit_inertia=(1e-8, 0.0, 0.0, 0.0, 1e-8, 0.0, 0.0, 0.0, 1e-8),
    )
    project = Project(
        new_id(),
        "Pendule G0",
        bodies=(body,),
        duration=parameters.duration,
        step=parameters.step,
    )
    return project.replace_object(joint_at(project, "rotule", None, body.id))


def import_g0(path):
    return pendulum(load_parameters(path))


def replace_cad_body(project, body):
    """Rebase attachment coordinates when a CAD edit moves the mass centre."""
    old = next((b for b in project.bodies if b.id == body.id), None)
    if old is None:
        return project.replace_object(body)
    old_r = np.array(old.orientation).reshape(3, 3)
    new_r = np.array(body.orientation).reshape(3, 3)
    delta = np.array(old.position) - body.position

    def point(p):
        return tuple(new_r.T @ (old_r @ p + delta))

    def frame(r):
        return tuple((new_r.T @ old_r @ np.array(r).reshape(3, 3)).flat)

    joints = []
    for joint in project.joints:
        values = {}
        for side in ("a", "b"):
            if getattr(joint, side) == body.id:
                values["p" + side] = point(getattr(joint, "p" + side))
                values["r" + side] = frame(getattr(joint, "r" + side))
        joints.append(replace(joint, **values))
    loads = tuple(
        replace(load, point=point(load.point)) if load.body == body.id else load
        for load in project.loads
    )
    return replace(project.replace_object(body), joints=tuple(joints), loads=loads)
