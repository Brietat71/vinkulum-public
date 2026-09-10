"""Immutable solid-feature graph; no geometry engine is imported by its reader."""

import hashlib
import json
import math
import uuid
from dataclasses import asdict, dataclass, fields, replace

import numpy as np

from .sketch import Sketch

IDENTITY = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
PARAMETERS = {
    "box": ("Length", "Width", "Height"),
    "cylinder": ("Radius", "Height"),
    "sphere": ("Radius",),
    "extrude_rectangle": ("Length", "Width", "Height"),
    "extrude_circle": ("Radius", "Height"),
    "extrude_sketch": ("Height",),
    "fillet": ("Radius",),
    "cut": (),
    "fuse": (),
    "common": (),
    "snapshot": (),
}


def _identifier(value):
    if not isinstance(value, str) or str(uuid.UUID(value)) != value:
        raise ValueError("A canonical CAD feature UUID is required.")


def _vector(value, count, name, limit):
    if (
        not isinstance(value, tuple)
        or len(value) != count
        or not all(
            type(v) in (int, float) and math.isfinite(v) and abs(v) <= limit
            for v in value
        )
    ):
        raise ValueError(f"Invalid CAD {name}.")


def _fields(data, cls):
    if not isinstance(data, dict) or set(data) != {f.name for f in fields(cls)}:
        raise ValueError(f"Unknown or missing {cls.__name__} fields.")
    return dict(data)


@dataclass(frozen=True)
class CadFeature:
    id: str
    name: str
    kind: str
    inputs: tuple = ()
    dimensions_mm: tuple = ()
    position_mm: tuple = (0.0, 0.0, 0.0)
    orientation: tuple = IDENTITY
    brep_mm: str | None = None
    source_body_id: str | None = None
    source_body_sha256: str | None = None

    def __post_init__(self):
        _identifier(self.id)
        if not isinstance(self.name, str) or not 1 <= len(self.name.strip()) <= 128:
            raise ValueError("A CAD feature name of 1–128 characters is required.")
        if self.kind not in PARAMETERS:
            raise ValueError("Unsupported CAD feature kind.")
        if self.kind == "extrude_sketch" and not isinstance(self, SketchExtrusion):
            raise ValueError("A sketch extrusion requires its persistent profile.")
        count = (
            2
            if self.kind in ("cut", "fuse", "common")
            else 1
            if self.kind == "fillet"
            else 0
        )
        if not isinstance(self.inputs, tuple) or len(self.inputs) != count:
            raise ValueError("Invalid number of CAD feature inputs.")
        for value in self.inputs:
            _identifier(value)
        if len(set(self.inputs)) != count:
            raise ValueError("CAD feature inputs must be distinct.")
        _vector(self.dimensions_mm, len(PARAMETERS[self.kind]), "dimensions", 1e6)
        if any(v < 0.001 for v in self.dimensions_mm):
            raise ValueError("CAD dimensions must be at least 0.001 mm.")
        _vector(self.position_mm, 3, "placement", 1e6)
        _vector(self.orientation, 9, "rotation", 1.00000001)
        R = np.array(self.orientation, dtype=float).reshape(3, 3)
        if not np.allclose(R.T @ R, np.eye(3), atol=1e-10, rtol=0) or not math.isclose(
            np.linalg.det(R), 1, abs_tol=1e-10
        ):
            raise ValueError("CAD feature orientation must be a proper rotation.")
        if count and (
            self.position_mm != (0.0, 0.0, 0.0) or self.orientation != IDENTITY
        ):
            raise ValueError(
                "Operation features inherit the design frame from their inputs."
            )
        if self.kind == "snapshot":
            if (
                not isinstance(self.brep_mm, str)
                or not 1 <= len(self.brep_mm) <= 2_000_000
                or not self.brep_mm.isascii()
                or "CASCADE Topology" not in self.brep_mm[:200]
            ):
                raise ValueError("A captured CAD input requires a bounded ASCII BREP.")
        elif self.brep_mm is not None:
            raise ValueError("Only a captured solid input can contain BREP.")
        if self.source_body_id is not None:
            _identifier(self.source_body_id)
            if (
                self.kind != "snapshot"
                or not isinstance(self.source_body_sha256, str)
                or len(self.source_body_sha256) != 64
                or any(v not in "0123456789abcdef" for v in self.source_body_sha256)
            ):
                raise ValueError(
                    "A captured tool requires its source body fingerprint."
                )
        elif self.source_body_sha256 is not None:
            raise ValueError("A tool fingerprint requires a source body identity.")

    @classmethod
    def from_dict(cls, data):
        if (
            cls is CadFeature
            and isinstance(data, dict)
            and data.get("kind") == "extrude_sketch"
        ):
            return SketchExtrusion.from_dict(data)
        data = _fields(data, cls)
        for key in ("inputs", "dimensions_mm", "position_mm", "orientation"):
            data[key] = tuple(data[key])
        if cls is SketchExtrusion:
            data["profile"] = Sketch.from_dict(data["profile"])
        return cls(**data)


@dataclass(frozen=True)
class SketchExtrusion(CadFeature):
    profile: Sketch | None = None

    def __post_init__(self):
        super().__post_init__()
        if self.kind != "extrude_sketch" or not isinstance(self.profile, Sketch):
            raise ValueError(
                "A sketch extrusion requires its persistent planar profile."
            )


@dataclass(frozen=True)
class CadRecipe:
    features: tuple
    root: str
    # Design-frame origin, expressed in the centred mechanical body's axes [m].
    origin_in_body_m: tuple = (0.0, 0.0, 0.0)

    def __post_init__(self):
        _identifier(self.root)
        _vector(self.origin_in_body_m, 3, "design origin", 1000)
        if (
            not isinstance(self.features, tuple)
            or not 1 <= len(self.features) <= 100
            or not all(isinstance(f, CadFeature) for f in self.features)
        ):
            raise ValueError("A CAD recipe requires 1–100 immutable features.")
        if len({f.id for f in self.features}) != len(self.features):
            raise ValueError("Duplicate CAD feature identities.")
        if self.embedded_brep_bytes > 2_000_000:
            raise ValueError("Captured CAD inputs exceed the 2 MB recipe budget.")
        self.ordered()

    @property
    def embedded_brep_bytes(self):
        return sum(len(f.brep_mm or "") for f in self.features)

    def ordered(self):
        """Dependency order is determined by identities, never array positions."""
        by_id = {f.id: f for f in self.features}
        active, visited, result = set(), set(), []

        def visit(identifier):
            if identifier in active:
                raise ValueError("Cycle in CAD feature dependencies.")
            if identifier not in by_id:
                raise ValueError("A referenced CAD feature is missing.")
            if identifier in visited:
                return
            active.add(identifier)
            feature = by_id[identifier]
            for source in feature.inputs:
                visit(source)
            active.remove(identifier)
            visited.add(identifier)
            result.append(feature)

        visit(self.root)
        if len(visited) != len(by_id):
            raise ValueError("The CAD recipe contains disconnected features.")
        return tuple(result)

    def replace_feature(self, feature):
        if not any(f.id == feature.id for f in self.features):
            raise ValueError("The CAD feature to edit is missing.")
        return replace(
            self,
            features=tuple(feature if f.id == feature.id else f for f in self.features),
        )

    def fingerprint(self):
        data = asdict(self)
        data["features"] = sorted(data["features"], key=lambda f: f["id"])
        profiles = {
            f.id: f.profile for f in self.features if isinstance(f, SketchExtrusion)
        }
        for feature in data["features"]:
            if feature["id"] in profiles:
                feature["profile"] = profiles[feature["id"]].canonical()
        return hashlib.sha256(
            json.dumps(
                data, sort_keys=True, separators=(",", ":"), allow_nan=False
            ).encode()
        ).hexdigest()

    @classmethod
    def from_dict(cls, data):
        data = _fields(data, cls)
        data["features"] = tuple(CadFeature.from_dict(f) for f in data["features"])
        data["origin_in_body_m"] = tuple(data["origin_in_body_m"])
        return cls(**data)


def recipe_for_body(body):
    """Older solids are captured explicitly; their prose journal is not replayed."""
    if body.cad and body.cad.recipe:
        return body.cad.recipe
    feature = CadFeature(
        str(uuid.uuid4()),
        body.name,
        "snapshot" if body.cad else body.shape,
        dimensions_mm=() if body.cad else tuple(v * 1000 for v in body.dimensions),
        brep_mm=body.cad.brep_mm if body.cad else None,
    )
    return CadRecipe((feature,), feature.id)
