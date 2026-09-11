"""Captured OCCT-8 solid to a bounded, face-labelled Gmsh tetrahedral mesh.

Gmsh is a separate executable. No CAD or mesher library is imported here.
Face identities belong to this captured mesh; they are not persistent CAD names.
"""

import hashlib
import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass, fields
from pathlib import Path

import numpy as np

from .document import Body
from .mesh_data import MAX_ELEMENTS, MAX_NODES, FiniteMesh
from .model import finite_number, read_json, write_json

MAX_MESH_BYTES = 32 * 1024 * 1024


def fingerprint(value):
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def _hash(value):
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError("Invalid mesh fingerprint.")


def _fields(data, cls):
    if not isinstance(data, dict) or set(data) != {f.name for f in fields(cls)}:
        raise ValueError(f"Unknown or missing {cls.__name__} fields.")
    return dict(data)


def read_bounded(path):
    with Path(path).open("rb") as stream:
        raw = stream.read(MAX_MESH_BYTES + 1)
    if len(raw) > MAX_MESH_BYTES:
        raise ValueError("Meshing file exceeded its 32 MiB budget.")
    return raw


@dataclass(frozen=True)
class MeshRequest:
    body: Body
    size_mm: float
    order: int = 2
    threads: int = 2

    def __post_init__(self):
        if not isinstance(self.body, Body) or self.body.cad is None:
            raise ValueError("Select one CAD solid before generating a volume mesh.")
        if int(self.body.cad.occt_version.split(".")[0]) < 8:
            raise ValueError("The captured CAD solid requires OCCT 8 or newer.")
        if not finite_number(self.size_mm) or not 0.001 <= self.size_mm <= 1e6:
            raise ValueError("Mesh size must be between 0.001 and 1,000,000 mm.")
        if type(self.order) is not int or self.order not in (1, 2):
            raise ValueError("Tetrahedron order must be 1 or 2.")
        if type(self.threads) is not int or self.threads < 1:
            raise ValueError("The mesher requires a positive integer thread count.")

    @property
    def element_type(self):
        return "C3D10" if self.order == 2 else "C3D4"

    @classmethod
    def from_dict(cls, data):
        if isinstance(data, dict) and "cad_face_witnesses" in data:
            values = _fields(data, CadWitnessMeshRequest)
            values["body"] = Body.from_dict(values["body"])
            return CadWitnessMeshRequest(**values)
        if isinstance(data, dict) and "algorithm" in data:
            values = _fields(data, HxtMeshRequest)
            values["body"] = Body.from_dict(values["body"])
            return HxtMeshRequest(**values)
        data = _fields(data, cls)
        data["body"] = Body.from_dict(data["body"])
        return cls(**data)

    def geo_script(self):
        # BREP is unitless: Studio's captured BREP uses millimetres. Mesh in
        # that local frame, then convert coordinates to world SI exactly once.
        return f"""SetFactory("OpenCASCADE");
General.AbortOnError = 2;
General.NumThreads = {self.threads};
Merge "solid.brep";
v() = Volume{{:}};
If (#v() != 1)
  Error("Vinkulum requires exactly one imported CAD solid");
  Exit;
EndIf
s() = Surface{{:}};
For i In {{0:#s()-1}}
  Physical Surface(Sprintf("Face %g", s(i)), s(i)) = {{s(i)}};
EndFor
Physical Volume("Solid", 1) = {{v()}};
Mesh.MeshSizeMin = {self.size_mm / 5:.17g};
Mesh.MeshSizeMax = {self.size_mm:.17g};
Mesh.MeshSizeFromCurvature = 12;
Mesh.Algorithm3D = 1;
Mesh.RandomSeed = 1;
Mesh.ElementOrder = {self.order};
Mesh.HighOrderOptimize = {2 if self.order == 2 else 0};
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
Mesh.SaveAll = 0;
"""


@dataclass(frozen=True)
class HxtMeshRequest(MeshRequest):
    """Parallel 3D Delaunay request; legacy requests retain their exact payload."""

    algorithm: str = "hxt"

    def __post_init__(self):
        super().__post_init__()
        if self.algorithm != "hxt":
            raise ValueError("Unsupported parallel meshing algorithm.")

    def geo_script(self):
        return (
            super()
            .geo_script()
            .replace(
                "Mesh.Algorithm3D = 1;",
                f"Mesh.Algorithm3D = 10;\nMesh.MaxNumThreads1D = {self.threads};\n"
                f"Mesh.MaxNumThreads2D = {self.threads};\nMesh.MaxNumThreads3D = {self.threads};",
            )
        )


@dataclass(frozen=True)
class CadWitnessMeshRequest(HxtMeshRequest):
    """Export CAD witnesses in the same Gmsh process that generates the mesh."""

    cad_face_witnesses: str = "occt-brep-witness-1"

    def __post_init__(self):
        super().__post_init__()
        if self.cad_face_witnesses != "occt-brep-witness-1":
            raise ValueError("Unsupported CAD face witness format.")

    def geo_script(self):
        return (
            super().geo_script()
            + """Geometry.OCCExportOnlyVisible = 1;
Mesh.Format = 10;
Hide { Volume{:}; Surface{:}; Curve{:}; Point{:}; }
For i In {0:#s()-1}
  Show { Surface{s(i)}; }
  Save Sprintf("cad-face-%g.brep", s(i));
  Hide { Surface{s(i)}; }
EndFor
Show { Volume{:}; Surface{:}; Curve{:}; Point{:}; }
Mesh.Format = 1;
"""
        )


@dataclass(frozen=True)
class MeshSurface:
    id: int
    name: str
    triangles: tuple

    def __post_init__(self):
        if (
            type(self.id) is not int
            or self.id < 1
            or not isinstance(self.name, str)
            or not 1 <= len(self.name) <= 128
        ):
            raise ValueError("Invalid captured mesh face identity.")
        if not isinstance(self.triangles, (tuple, list)) or not self.triangles:
            raise ValueError("A mesh face requires boundary triangles.")
        object.__setattr__(
            self, "triangles", tuple(tuple(row) for row in self.triangles)
        )

    @classmethod
    def from_dict(cls, data):
        return cls(**_fields(data, cls))


@dataclass(frozen=True)
class MeshedSolid:
    request: MeshRequest
    mesh: FiniteMesh
    surfaces: tuple
    gmsh_version: str
    occ_version: str
    executable_sha256: str
    input_sha256: str
    raw_mesh_sha256: str
    log_sha256: str

    def __post_init__(self):
        if not isinstance(self.request, MeshRequest) or not isinstance(
            self.mesh, FiniteMesh
        ):
            raise TypeError(
                "Captured geometry and a validated volume mesh are required."
            )
        if isinstance(self.request, CadWitnessMeshRequest) and not hasattr(
            self, "cad_faces"
        ):
            raise ValueError("CAD witness payload is missing.")
        if self.mesh.element_type != self.request.element_type:
            raise ValueError("Mesh element family differs from the captured request.")
        if (
            not isinstance(self.surfaces, tuple)
            or not self.surfaces
            or len(self.surfaces) > 10000
            or not all(isinstance(s, MeshSurface) for s in self.surfaces)
        ):
            raise ValueError("Invalid mesh surface groups.")
        if len({s.id for s in self.surfaces}) != len(self.surfaces):
            raise ValueError("Duplicate mesh face identity.")
        for version in (self.gmsh_version, self.occ_version):
            if (
                not isinstance(version, str)
                or re.fullmatch(r"\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?", version)
                is None
            ):
                raise ValueError("Invalid meshing engine version.")
        if int(self.occ_version.split(".")[0]) < 8:
            raise ValueError("Gmsh must use OCCT 8 or newer.")
        for value in (
            self.executable_sha256,
            self.input_sha256,
            self.raw_mesh_sha256,
            self.log_sha256,
        ):
            _hash(value)
        if self.input_sha256 != fingerprint(asdict(self.request)):
            raise ValueError(
                "Mesh input fingerprint differs from the captured solid/parameters."
            )
        expected = {tuple(sorted(row)): row for row in self.mesh.boundary_faces}
        found = []
        for surface in self.surfaces:
            for row in surface.triangles:
                if any(type(i) is not int for i in row) or tuple(row) != expected.get(
                    tuple(sorted(row))
                ):
                    raise ValueError(
                        "A mesh face differs from the outward volume boundary."
                    )
                found.append(tuple(sorted(row)))
        if len(found) != len(set(found)) or set(found) != set(expected):
            raise ValueError(
                "Surface groups must cover the volume boundary exactly once."
            )

    @classmethod
    def from_dict(cls, data):
        kind = (
            CadWitnessMeshedSolid
            if isinstance(data, dict) and "cad_faces" in data
            else cls
        )
        data = _fields(data, kind)
        data["request"] = MeshRequest.from_dict(data["request"])
        data["mesh"] = FiniteMesh(**_fields(data["mesh"], FiniteMesh))
        data["surfaces"] = tuple(MeshSurface.from_dict(s) for s in data["surfaces"])
        return kind(**data)

    @property
    def sha256(self):
        return fingerprint(asdict(self))


@dataclass(frozen=True)
class CadWitnessMeshedSolid(MeshedSolid):
    cad_faces: tuple

    def __post_init__(self):
        super().__post_init__()
        if not isinstance(self.request, CadWitnessMeshRequest):
            raise TypeError("CAD witnesses require a captured witness request.")
        rows = tuple(tuple(row) for row in self.cad_faces)
        if any(len(row) != 2 or type(row[0]) is not int for row in rows):
            raise ValueError("Invalid CAD face witness record.")
        if [row[0] for row in rows] != sorted(surface.id for surface in self.surfaces):
            raise ValueError("CAD witnesses must cover each mesh surface exactly once.")
        for _, digest in rows:
            _hash(digest)
        object.__setattr__(self, "cad_faces", rows)


def checked_cad_witnesses(root, surfaces):
    expected = {f"cad-face-{surface.id}.brep" for surface in surfaces}
    if {p.name for p in root.glob("cad-face-*.brep")} != expected:
        raise ValueError("CAD face witness files differ from the mesh surface set.")
    total = 0
    rows = []
    for surface in sorted(surfaces, key=lambda item: item.id):
        raw = read_bounded(root / f"cad-face-{surface.id}.brep")
        total += len(raw)
        if not raw.startswith(
            (b"\nCASCADE Topology V", b"DBRep_DrawableShape\n\nCASCADE Topology V")
        ):
            raise ValueError("CAD witness is not an ASCII OCCT BREP export.")
        if total > MAX_MESH_BYTES:
            raise ValueError("CAD face witnesses exceed the capture budget.")
        rows.append((surface.id, sha256(raw)))
    return tuple(rows)


def parse_msh(raw, request, cancellation=None):
    """Read only the captured MSH 2.2 ASCII triangle/tetrahedron contract."""
    if len(raw) > MAX_MESH_BYTES:
        raise ValueError("Raw mesh exceeded its budget.")
    lines = raw.decode("ascii").splitlines()
    sections = {}
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        key = lines[i].strip()
        if (
            key not in ("$MeshFormat", "$PhysicalNames", "$Nodes", "$Elements")
            or key in sections
        ):
            raise ValueError("Unsupported or duplicate Gmsh section.")
        end = "$End" + key[1:]
        i += 1
        start = i
        while i < len(lines) and lines[i].strip() != end:
            i += 1
        if i == len(lines):
            raise ValueError("Incomplete Gmsh section.")
        sections[key] = lines[start:i]
        i += 1
    if set(sections) != {
        "$MeshFormat",
        "$PhysicalNames",
        "$Nodes",
        "$Elements",
    } or sections["$MeshFormat"] != ["2.2 0 8"]:
        raise ValueError("MSH 2.2 ASCII with named boundary groups is required.")

    def counted(key, maximum):
        section = sections[key]
        if not section:
            raise ValueError("Missing Gmsh section count.")
        count = int(section[0])
        if not 1 <= count <= maximum or len(section) != count + 1:
            raise ValueError("Invalid or oversized Gmsh section count.")
        return section[1:]

    names = {}
    for line in counted("$PhysicalNames", 10001):
        match = re.fullmatch(r'(2|3)\s+(\d+)\s+"([^"\r\n]+)"', line)
        if not match:
            raise ValueError("Invalid Gmsh physical name.")
        key = tuple(map(int, match.group(1, 2)))
        if key in names:
            raise ValueError("Duplicate Gmsh physical name.")
        names[key] = match[3]
    coordinates = {}
    for line in counted("$Nodes", MAX_NODES):
        row = line.split()
        if len(row) != 4:
            raise ValueError("Malformed Gmsh node.")
        identifier = int(row[0])
        point = tuple(map(float, row[1:]))
        if (
            identifier < 1
            or identifier in coordinates
            or not all(finite_number(v) for v in point)
        ):
            raise ValueError("Invalid or duplicate Gmsh node.")
        coordinates[identifier] = point
    node_ids = sorted(coordinates)
    remap = {old: i + 1 for i, old in enumerate(node_ids)}
    rotation = np.array(request.body.orientation).reshape(3, 3)
    world = (
        np.array([coordinates[n] for n in node_ids]) * 0.001 @ rotation.T
        + request.body.position
    )
    volumes = []
    boundary = []
    element_ids = set()
    volume_type = 11 if request.order == 2 else 4
    surface_type = 9 if request.order == 2 else 2
    for line in counted("$Elements", MAX_ELEMENTS * 5):
        row = list(map(int, line.split()))
        if len(row) < 6:
            raise ValueError("Malformed Gmsh element.")
        identifier, kind, ntags = row[:3]
        if identifier < 1 or identifier in element_ids or not 2 <= ntags <= 8:
            raise ValueError("Duplicate element or invalid Gmsh tags.")
        element_ids.add(identifier)
        group, entity = row[3:5]
        try:
            nodes = tuple(remap[n] for n in row[3 + ntags :])
        except KeyError as error:
            raise ValueError("Gmsh element references a missing node.") from error
        if kind == volume_type and group == 1:
            if len(nodes) != (10 if request.order == 2 else 4):
                raise ValueError("Wrong tetrahedron node count.")
            # Gmsh edges 9 and 10 are (3,2) and (3,1); CalculiX reverses them.
            if request.order == 2:
                nodes = nodes[:8] + (nodes[9], nodes[8])
            volumes.append((identifier, nodes))
        elif kind == surface_type and (2, group) in names and group == entity:
            if len(nodes) != (6 if request.order == 2 else 3):
                raise ValueError("Wrong boundary triangle node count.")
            boundary.append((group, nodes))
        else:
            raise ValueError("Unsupported or unlabelled mesh element.")
    mesh = FiniteMesh(
        tuple(tuple(map(float, p)) for p in world),
        tuple(nodes for _, nodes in sorted(volumes)),
        request.element_type,
        cancellation=cancellation,
    )
    outer = {tuple(sorted(row)): row for row in mesh.boundary_faces}
    groups = {}
    for group, nodes in boundary:
        canonical = outer.get(tuple(sorted(nodes)))
        if canonical is None:
            raise ValueError("Gmsh triangle is not a conforming exterior face.")
        # Also validate midside-to-edge assignments before normalising orientation.
        if request.order == 2:
            edges = {
                frozenset((canonical[a], canonical[b])): canonical[3 + j]
                for j, (a, b) in enumerate(((0, 1), (1, 2), (2, 0)))
            }
            if any(
                edges.get(frozenset((nodes[a], nodes[b]))) != nodes[3 + j]
                for j, (a, b) in enumerate(((0, 1), (1, 2), (2, 0)))
            ):
                raise ValueError("Invalid quadratic boundary node ordering.")
        groups.setdefault(group, []).append(canonical)
    surfaces = tuple(
        MeshSurface(i, names[(2, i)], tuple(rows)) for i, rows in sorted(groups.items())
    )
    return mesh, surfaces


def identify_mesher(executable):
    found = shutil.which(str(executable))
    if not found:
        raise FileNotFoundError("Install a Gmsh executable built with OCCT 8 or newer.")
    path = Path(found).resolve()
    return path, sha256(path.read_bytes())


def parse_info(text, request=None):
    version = re.search(r"^Version\s*:\s*(\S+)", text, re.MULTILINE)
    occ = re.search(r"^OCC version\s*:\s*(\S+)", text, re.MULTILINE)
    if not version or not occ or int(occ[1].split(".")[0]) < 8:
        raise ValueError("Gmsh must report an OCCT 8 or newer geometry kernel (-info).")
    if isinstance(request, HxtMeshRequest):
        options = re.search(r"^Build options\s*:\s*(.+)$", text, re.MULTILINE)
        if not options or not {"hxt", "openmp"} <= set(options[1].lower().split()):
            raise ValueError(
                "Parallel HXT meshing requires a Gmsh build with Hxt and OpenMP."
            )
    return version[1], occ[1]


@dataclass(frozen=True)
class PreparedMesh:
    request: MeshRequest
    root: Path
    executable: Path
    executable_sha256: str
    gmsh_version: str
    occ_version: str


def prepare_mesh(request, directory, executable, executable_sha256, info):
    gmsh, occ = parse_info(info, request)
    root = Path(directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    write_json(root / "input.json", asdict(request))
    (root / "solid.brep").write_text(request.body.cad.brep_mm, encoding="ascii")
    (root / "mesh.geo").write_text(request.geo_script(), encoding="ascii")
    (root / "engine-info.txt").write_text(info)
    return PreparedMesh(request, root, executable, executable_sha256, gmsh, occ)


def finish_mesh(run, returncode, *, publish=True, cancellation=None):
    if returncode != 0:
        raise RuntimeError(
            f"Gmsh exited with code {returncode}; inspect {run.root / 'mesher.log'}."
        )
    log = read_bounded(run.root / "mesher.log")
    if re.search(rb"^Error\s*:", log, re.MULTILINE):
        raise RuntimeError("Gmsh reported an error; mesh not accepted.")
    if read_bounded(run.root / "solid.brep") != run.request.body.cad.brep_mm.encode(
        "ascii"
    ) or read_bounded(run.root / "mesh.geo") != run.request.geo_script().encode(
        "ascii"
    ):
        raise ValueError("CAD or meshing script changed during execution.")
    if read_json(run.root / "input.json", MAX_MESH_BYTES) != json.loads(
        json.dumps(asdict(run.request))
    ):
        raise ValueError("Captured meshing parameters changed during execution.")
    if sha256(run.executable.read_bytes()) != run.executable_sha256:
        raise ValueError("The meshing executable changed during execution.")
    raw = read_bounded(run.root / "mesh.msh")
    mesh, surfaces = parse_msh(raw, run.request, cancellation)
    kind = (
        CadWitnessMeshedSolid
        if isinstance(run.request, CadWitnessMeshRequest)
        else MeshedSolid
    )
    extra = (
        {"cad_faces": checked_cad_witnesses(run.root, surfaces)}
        if kind is CadWitnessMeshedSolid
        else {}
    )
    result = kind(
        run.request,
        mesh,
        surfaces,
        run.gmsh_version,
        run.occ_version,
        run.executable_sha256,
        fingerprint(asdict(run.request)),
        sha256(raw),
        sha256(log),
        **extra,
    )
    if cancellation is not None and cancellation():
        raise InterruptedError("Mesh validation cancelled.")
    if publish:
        save_mesh_result(run, result)
    return result


def save_mesh_result(run, result, *, filename="result.json"):
    write_json(
        run.root / filename,
        {"format": "vinkulum-solid-mesh", "schema": 1, "result": asdict(result)},
    )


def load_mesh(directory):
    root = Path(directory)
    if root.is_file() and root.name == "result.json":
        root = root.parent
    data = read_json(root / "result.json", MAX_MESH_BYTES)
    if (
        not isinstance(data, dict)
        or set(data) != {"format", "schema", "result"}
        or data["format"] != "vinkulum-solid-mesh"
        or type(data["schema"]) is not int
        or data["schema"] != 1
    ):
        raise ValueError("Unsupported solid mesh archive.")
    result = MeshedSolid.from_dict(data["result"])
    raw = read_bounded(root / "mesh.msh")
    log = read_bounded(root / "mesher.log")
    if sha256(raw) != result.raw_mesh_sha256 or sha256(log) != result.log_sha256:
        raise ValueError("Captured mesher output fingerprint does not match.")
    if re.search(rb"\bError\s*:", log, re.IGNORECASE):
        raise ValueError("Archived mesher log contains an error.")
    if (
        read_bounded(root / "solid.brep")
        != result.request.body.cad.brep_mm.encode("ascii")
        or read_bounded(root / "mesh.geo")
        != result.request.geo_script().encode("ascii")
        or read_json(root / "input.json", MAX_MESH_BYTES)
        != json.loads(json.dumps(asdict(result.request)))
    ):
        raise ValueError("Captured CAD/meshing inputs are inconsistent.")
    if parse_info(read_bounded(root / "engine-info.txt").decode(), result.request) != (
        result.gmsh_version,
        result.occ_version,
    ):
        raise ValueError("Captured mesher versions are inconsistent.")
    mesh, surfaces = parse_msh(raw, result.request)
    if mesh != result.mesh or surfaces != result.surfaces:
        raise ValueError("Mesh payload does not reproduce the raw mesher output.")
    if (
        isinstance(result, CadWitnessMeshedSolid)
        and checked_cad_witnesses(root, surfaces) != result.cad_faces
    ):
        raise ValueError("Captured CAD face witness fingerprint does not match.")
    return result, root


def run_mesh(request, directory, *, executable="gmsh", timeout=120):
    if not finite_number(timeout) or not 0 < timeout <= 600:
        raise ValueError("Meshing timeout must be between 0 and 600 seconds.")
    engine, identity = identify_mesher(executable)
    from .engine_threads import process_environment

    environment = process_environment(request.threads, engine="gmsh")
    probe = subprocess.run(
        [str(engine), "-info"],
        capture_output=True,
        text=True,
        timeout=10,
        check=True,
        env=environment,
    )
    run = prepare_mesh(
        request, directory, engine, identity, probe.stdout + probe.stderr
    )
    with (run.root / "mesher.log").open("wb") as log:
        process = subprocess.run(
            [
                str(engine),
                "-nt",
                str(request.threads),
                "mesh.geo",
                "-3",
                "-format",
                "msh2",
                "-o",
                "mesh.msh",
            ],
            cwd=run.root,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
            env=environment,
        )
    return finish_mesh(run, process.returncode)
