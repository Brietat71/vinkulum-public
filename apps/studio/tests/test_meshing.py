"""Independent tetrahedral boundary identities, physical loads and archive checks."""

import hashlib
import json
import math
import sys
import tempfile
import unittest
import uuid
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from vinkulum_studio.calculix import load_study, study_document
from vinkulum_studio.document import load_project
from vinkulum_studio.mesh_binding import MeshBinding, MeshCondition, mesh_measurements
from vinkulum_studio.meshing import (
    MeshedSolid,
    MeshRequest,
    fingerprint,
    finish_mesh,
    load_mesh,
    parse_info,
    parse_msh,
    prepare_mesh,
)

ROOT = Path(__file__).resolve().parents[3]


def request(order=2):
    body = load_project(ROOT / "examples/studio/platine-percee.vinkulum.json").bodies[0]
    return MeshRequest(replace(body, position=(0.0, 0.0, 0.0)), 10.0, order)


def tetra_msh(order=2, curve=0.0):
    # Artificial 1 m simplex, independent of the source CAD geometry: these
    # tests qualify parsing and boundary integration, not CAD fidelity.
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    mids = {(1, 2): 5, (2, 3): 6, (1, 3): 7, (1, 4): 8, (2, 4): 9, (3, 4): 10}
    if order == 2:
        points += [
            (0.5, 0.0, curve),
            (0.5, 0.5, 0.0),
            (0.0, 0.5, 0.0),
            (0.0, 0.0, 0.5),
            (0.5, 0.0, 0.5),
            (0.0, 0.5, 0.5),
        ]
    node = lambda i: i * 10 + 3
    rows = []
    for i, corners in enumerate(((1, 3, 2), (1, 2, 4), (2, 3, 4), (3, 1, 4)), 1):
        ids = corners
        if order == 2:
            ids += tuple(
                mids[tuple(sorted((corners[a], corners[b])))]
                for a, b in ((0, 1), (1, 2), (2, 0))
            )
        rows.append(
            f"{i} {9 if order == 2 else 2} 2 {i} {i} "
            + " ".join(str(node(n)) for n in ids)
        )
    ids = (1, 2, 3, 4, 5, 6, 7, 8, 10, 9) if order == 2 else (1, 2, 3, 4)
    rows += [
        f"5 {11 if order == 2 else 4} 2 1 1 " + " ".join(str(node(n)) for n in ids)
    ]
    names = [f'2 {i} "Face {i}"' for i in range(1, 5)] + ['3 1 "Solid"']
    nodes = [
        f"{node(i)} " + " ".join(format(v * 1000, ".17g") for v in p)
        for i, p in enumerate(points, 1)
    ]
    text = [
        "$MeshFormat",
        "2.2 0 8",
        "$EndMeshFormat",
        "$PhysicalNames",
        "5",
        *names,
        "$EndPhysicalNames",
        "$Nodes",
        str(len(nodes)),
        *nodes,
        "$EndNodes",
        "$Elements",
        "5",
        *rows,
        "$EndElements",
    ]
    return ("\n".join(text) + "\n").encode("ascii")


def solid(order=2, curve=0.0):
    req, raw = request(order), tetra_msh(order, curve)
    mesh, surfaces = parse_msh(raw, req)
    return MeshedSolid(
        req,
        mesh,
        surfaces,
        "5.0.0-test",
        "8.0.1",
        "a" * 64,
        fingerprint(asdict(req)),
        hashlib.sha256(raw).hexdigest(),
        "b" * 64,
    )


def condition(kind, surfaces, *, axes=(), values=()):
    return MeshCondition(str(uuid.uuid4()), kind, kind, tuple(surfaces), axes, values)


class MeshContracts(unittest.TestCase):
    def test_units_ordering_world_transform_and_boundary_coverage(self):
        s = solid()
        self.assertEqual(s.mesh.elements, (tuple(range(1, 11)),))
        self.assertEqual(s.mesh.nodes[1], (1.0, 0.0, 0.0))
        self.assertEqual(len(s.surfaces), 4)
        req = s.request
        rotation = (0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
        changed = replace(
            req, body=replace(req.body, orientation=rotation, position=(2.0, 3.0, 4.0))
        )
        transformed, _ = parse_msh(tetra_msh(), changed)
        np.testing.assert_array_equal(
            transformed.nodes,
            np.array(s.mesh.nodes) @ np.array(rotation).reshape(3, 3).T + (2, 3, 4),
        )

    def test_missing_counts_ids_faces_orders_and_legacy_occt_are_rejected(self):
        req, raw = request(), tetra_msh()
        invalid = (
            raw.replace(b"$Nodes\n10", b"$Nodes"),
            raw.replace(b"$Nodes\n10", b"$Nodes\n6001"),
            raw.replace(b"103 ", b"13 "),
            raw.replace(b"5 11 2 1 1", b"5 4 2 1 1"),
            raw.replace(
                b"5 11 2 1 1 13 23 33 43 53 63 73 83 103 93",
                b"5 11 2 1 1 13 23 33 43 53 63 73 83 93 103",
            ),
            raw.replace(
                b"$EndElements", b"$EndElements\n$Comment\nignored\n$EndComment"
            ),
        )
        for data in invalid:
            self.assertNotEqual(data, raw)
            with self.assertRaises(ValueError):
                parse_msh(data, req)
        empty_count = raw.replace(
            b"$Nodes\n10\n" + raw.split(b"$Nodes\n10\n")[1].split(b"$EndNodes")[0],
            b"$Nodes\n",
        )
        with self.assertRaisesRegex(ValueError, "count"):
            parse_msh(empty_count, req)
        with self.assertRaisesRegex(ValueError, "OCCT 8"):
            parse_info("Version : 4.12.1\nOCC version : 7.6.3\n")
        with self.assertRaises(InterruptedError):
            parse_msh(raw, req, cancellation=lambda: True)

    def test_curved_pressure_closed_surface_force_moment_and_volume(self):
        for order, curve in ((1, 0.0), (2, 0.0), (2, 0.05)):
            s = solid(order, curve)
            binding = MeshBinding(
                s, (condition("pressure", (1, 2, 3, 4), values=(120.0,)),)
            )
            forces = np.zeros((len(s.mesh.nodes), 3))
            for i, *values in binding.arrays[1]:
                forces[i - 1] = values
            np.testing.assert_allclose(forces.sum(axis=0), 0.0, atol=2e-13)
            np.testing.assert_allclose(
                np.cross(s.mesh.nodes, forces).sum(axis=0), 0.0, atol=2e-13
            )
            # x=r, y=s, z=t+4*c*r*(1-r-s-t): detJ=1-4*c*r.
            reference = (1 - curve) / 6
            self.assertTrue(
                math.isclose(
                    mesh_measurements(s)["signed_mesh_volume_m3"],
                    reference,
                    rel_tol=3e-14,
                )
            )
        planar = solid()
        binding = MeshBinding(planar, (condition("pressure", (1,), values=(120.0,)),))
        resultant = np.array([r[1:] for r in binding.arrays[1]]).sum(axis=0)
        # Face z=0 has outward normal -z. Positive pressure is compressive +z.
        np.testing.assert_allclose(resultant, (0, 0, 60), atol=2e-13)

    def test_total_force_supports_scaling_and_physical_intent_are_preserved(self):
        s = solid(curve=0.05)
        conditions = (
            condition("support", (1,), axes=(1, 2, 3)),
            condition("total_force", (2, 3), values=(10.0, -20.0, 30.0)),
        )
        binding = MeshBinding(s, conditions)
        np.testing.assert_allclose(
            np.array([r[1:] for r in binding.arrays[1]]).sum(axis=0),
            (10.0, -20.0, 30.0),
            atol=2e-13,
        )
        study = binding.study(210e9, 0.3)
        scaled = study.scaled_loads(-2.0)
        np.testing.assert_allclose(
            np.array([r[1:] for r in scaled.forces]).sum(axis=0),
            (-20, 40, -60),
            atol=5e-13,
        )
        self.assertEqual(scaled.fixed_dofs, study.fixed_dofs)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "study.json"
            data = study_document(study)
            path.write_text(json.dumps(data))
            self.assertEqual(load_study(path), study)
            changed = {**data, "schema": 2}
            path.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "schema 3"):
                load_study(path)
            data["study"]["forces"] = ((1, 11.0, -20.0, 30.0),)
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "physical intent"):
                load_study(path)

    def test_archive_rechecks_raw_mesh_and_rejects_error_logs(self):
        req = request()
        engine = Path(sys.executable).resolve()
        identity = hashlib.sha256(engine.read_bytes()).hexdigest()
        with tempfile.TemporaryDirectory() as directory:
            run = prepare_mesh(
                req,
                Path(directory) / "run",
                engine,
                identity,
                "Version : 5.0.0-test\nOCC version : 8.0.1\n",
            )
            (run.root / "mesh.msh").write_bytes(tetra_msh())
            (run.root / "mesher.log").write_text("Synthetic parser fixture\n")
            captured = finish_mesh(run, 0)
            self.assertEqual(load_mesh(run.root)[0], captured)
            # Even matching fingerprints must not admit an explicit error log.
            error = b"Error : failed meshing\n"
            (run.root / "mesher.log").write_bytes(error)
            data = json.loads((run.root / "result.json").read_text())
            data["result"]["log_sha256"] = hashlib.sha256(error).hexdigest()
            (run.root / "result.json").write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "log contains an error"):
                load_mesh(run.root)


if __name__ == "__main__":
    unittest.main()
