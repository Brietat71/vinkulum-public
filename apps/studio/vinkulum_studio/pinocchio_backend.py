"""Optional Pinocchio worker: operators at one captured articulated state.

This module has no Qt, VTK, OCCT or native Vinkulum dependency. Import Pinocchio
only in its dedicated environment. No trajectory integrator is implied.
"""

import argparse
import hashlib
import json
import platform
import sys
from importlib import metadata
from pathlib import Path

import numpy as np

from .articulated import operator_conventions, state_vector, tree_links
from .document import load_project, save_project
from .model import finite_number, read_json, write_json

PINOCCHIO_VERSION = "4.1.0"
ADAPTER_VERSION = "0.1.0"
STATE_KEYS = {"q", "velocity", "acceleration", "effort", "time_s"}


def _engine_fingerprints():
    files = {}
    for name in ("pin", "libpinocchio"):
        distribution = metadata.distribution(name)
        for file in distribution.files or ():
            if ".so" in file.name:
                path = distribution.locate_file(file)
                files[f"{name}/{file}"] = hashlib.sha256(path.read_bytes()).hexdigest()
    if not files:
        raise ValueError("Could not identify the installed Linux Pinocchio binaries.")
    return files


def _symmetric_upper(matrix):
    upper = np.triu(np.array(matrix, copy=True))
    return upper + np.triu(upper, 1).T


class ArticulatedModel:
    def __init__(self, project):
        self.project = project
        self.links = tree_links(project)
        import pinocchio as pin

        if pin.__version__ != PINOCCHIO_VERSION:
            raise ValueError(
                f"The adapter requires qualified Pinocchio {PINOCCHIO_VERSION}; found {pin.__version__}."
            )
        self.pin = pin
        self.model = model = pin.Model()
        model.gravity.linear = np.array(project.gravity)
        model.gravity.angular = np.zeros(3)
        # Keep kinematic Jacobian lever arms near the assembly, even when its
        # world coordinates are kilometres from zero. Output poses are restored.
        self.world_origin = np.array(self.links[0].parent_anchor[0])
        self.joint_ids = {}
        self.frame_ids = {}
        body_placements = {None: pin.SE3.Identity()}
        joint_ids = {None: 0}
        frames = {None: 0}
        bodies = {body.id: body for body in project.bodies}
        for link in self.links:
            parent_point, parent_rotation = link.parent_anchor
            child_point, child_rotation = link.child_anchor
            mount = pin.SE3(
                np.array(parent_rotation).reshape(3, 3), np.array(parent_point)
            )
            if link.parent is None:
                mount.translation -= self.world_origin
            body_placement = pin.SE3(
                np.array(child_rotation).reshape(3, 3), np.array(child_point)
            ).inverse()
            # Re-rooting reverses the axis, not the user's A→B coordinate.
            axis = np.array((0.0, 0.0, float(link.sign)))
            joint_model = (
                pin.JointModelRevoluteUnaligned(axis)
                if link.joint.kind == "pivot"
                else pin.JointModelPrismaticUnaligned(axis)
            )
            joint = model.addJoint(
                joint_ids[link.parent],
                joint_model,
                body_placements[link.parent] * mount,
                link.joint.id,
            )
            body = bodies[link.child]
            inertia = np.array(body.inertia()).reshape(3, 3)
            # Document origins are body centres of mass, including centred CAD.
            model.appendBodyToJoint(
                joint, pin.Inertia(body.mass, np.zeros(3), inertia), body_placement
            )
            joint_frame = model.addJointFrame(joint, frames[link.parent])
            frame = model.addBodyFrame(body.id, joint, body_placement, joint_frame)
            self.joint_ids[link.joint.id] = joint
            self.frame_ids[body.id] = frame
            joint_ids[body.id], frames[body.id] = joint, frame
            body_placements[body.id] = body_placement
        if model.nq != len(self.links) or model.nv != len(self.links):
            raise ValueError("Unexpected coordinate dimensions from Pinocchio.")
        self.initial_q = state_vector(None, self.links, "q", initial=True)
        initial = self.kinematics(self.initial_q)
        for body in project.bodies:
            pose, _ = initial[body.id]
            # Account for subtraction/addition at a translated world origin.
            tolerance = 1e-8 + 128 * np.finfo(float).eps * max(
                1.0, np.max(np.abs(body.position))
            )
            if not np.allclose(
                pose.translation, body.position, rtol=0, atol=tolerance
            ) or not np.allclose(
                pose.rotation,
                np.array(body.orientation).reshape(3, 3),
                rtol=0,
                atol=1e-8,
            ):
                raise ValueError(
                    f"{body.name}: Pinocchio conversion did not preserve the captured pose."
                )

    def kinematics(self, q):
        pin = self.pin
        data = self.model.createData()
        pin.computeJointJacobians(self.model, data, q)
        pin.updateFramePlacements(self.model, data)
        result = {}
        for identifier, frame in self.frame_ids.items():
            pose = data.oMf[frame].copy()
            pose.translation += self.world_origin
            # EigenPy squeezes a 6×1 matrix into a vector for a one-DOF model.
            jacobian = np.array(
                pin.getFrameJacobian(self.model, data, frame, pin.LOCAL_WORLD_ALIGNED),
                copy=True,
            ).reshape(6, self.model.nv)
            result[identifier] = pose, jacobian
        return result

    def analyse(
        self, *, q=None, velocity=None, acceleration=None, effort=None, time_s=0.0
    ):
        if not finite_number(time_s) or not 0 <= time_s <= self.project.duration:
            raise ValueError(
                "Load evaluation time must lie within the captured project duration."
            )
        q = state_vector(q, self.links, "q", initial=True)
        v = state_vector(velocity, self.links, "velocity")
        a = state_vector(acceleration, self.links, "acceleration")
        tau = state_vector(effort, self.links, "effort")
        pin, model = self.pin, self.model
        data = model.createData()
        kinematics = self.kinematics(q)
        external = np.zeros(model.nv)
        load_values = []
        for load in self.project.loads:
            pose, jacobian = kinematics[load.body]
            arm = pose.rotation @ load.point
            force = np.array([law.value(time_s) for law in load.force])
            moment = np.array([law.value(time_s) for law in load.moment])
            point_jacobian = jacobian[:3] + np.cross(jacobian[3:].T, arm).T
            external += point_jacobian.T @ force + jacobian[3:].T @ moment
            load_values.append(
                {
                    "id": load.id,
                    "body": load.body,
                    "point_world_m": (pose.translation + arm).tolist(),
                    "force_world_N": force.tolist(),
                    "moment_world_Nm": moment.tolist(),
                }
            )
        mass = _symmetric_upper(pin.crba(model, data, q))
        np.linalg.cholesky(mass)
        bias = pin.rnea(model, data, q, v, np.zeros(model.nv)).copy()
        inverse = pin.rnea(model, data, q, v, a).copy()
        forward = pin.aba(model, data, q, v, tau + external).copy()
        recovered = pin.rnea(model, data, q, v, forward).copy()
        derivatives = [
            np.array(part, copy=True).reshape(model.nv, model.nv)
            for part in pin.computeRNEADerivatives(model, data, q, v, a)
        ]
        derivatives[2] = _symmetric_upper(derivatives[2])
        # These intrinsic derivatives deliberately exclude configuration-dependent
        # applied world loads; their contribution is reported separately.
        residual = inverse - (mass @ a + bias)
        forward_residual = recovered - (tau + external)
        for error, scale in (
            (residual, np.abs(inverse) + np.abs(mass) @ np.abs(a) + np.abs(bias)),
            (forward_residual, np.abs(recovered) + np.abs(tau) + np.abs(external)),
        ):
            if not np.all(np.abs(error) <= 1e-12 + 1e-10 * scale):
                raise ValueError("Pinocchio operator consistency check failed.")
        body_values = []
        potential = 0.0
        for body in self.project.bodies:
            pose, jacobian = kinematics[body.id]
            potential -= body.mass * np.array(self.project.gravity) @ pose.translation
            body_values.append(
                {
                    "id": body.id,
                    "name": body.name,
                    "position_m": pose.translation.tolist(),
                    "body_to_world": pose.rotation.tolist(),
                    "jacobian": jacobian.tolist(),
                    "velocity_world": (jacobian @ v).tolist(),
                }
            )
        result = {
            "format": "vinkulum-pinocchio-operators",
            "schema": 1,
            "engine": "Pinocchio",
            "engine_version": pin.__version__,
            "adapter_version": ADAPTER_VERSION,
            "scientific_status": "NotAssessed",
            "nq": model.nq,
            "nv": model.nv,
            "coordinate_order": [
                {
                    "joint_id": link.joint.id,
                    "name": link.joint.name,
                    "kind": link.joint.kind,
                    "parent_body": link.parent,
                    "child_body": link.child,
                    "coordinate_unit": link.coordinate_unit,
                    "effort_unit": link.effort_unit,
                    "positive_direction": "B relative to A along joint A local +Z",
                }
                for link in self.links
            ],
            "state": {
                "q": q.tolist(),
                "velocity": v.tolist(),
                "acceleration": a.tolist(),
                "effort": tau.tolist(),
                "time_s": time_s,
            },
            "mass_matrix": mass.tolist(),
            "intrinsic_bias": bias.tolist(),
            "external_effort": external.tolist(),
            "inverse_effort": (inverse - external).tolist(),
            "forward_acceleration": forward.tolist(),
            "intrinsic_inverse_derivatives": {
                key: value.tolist()
                for key, value in zip(("q", "velocity", "acceleration"), derivatives)
            },
            "kinetic_energy_J": float(0.5 * v @ mass @ v),
            "potential_energy_J": float(potential),
            "bodies": body_values,
            "loads": load_values,
            "checks": {
                "inverse_residual": residual.tolist(),
                "forward_effort_residual": forward_residual.tolist(),
                "absolute_tolerance_effort_units": 1e-12,
                "relative_tolerance": 1e-10,
            },
            "conventions": operator_conventions(),
        }
        # Reject nonfinite channels before crossing the worker/JSON boundary.
        json.dumps(result, allow_nan=False)
        return result


def run_analysis(project, directory, state=None):
    state = {} if state is None else state
    if not isinstance(state, dict) or not set(state) <= STATE_KEYS:
        raise ValueError("Unknown or invalid articulated state fields.")
    root = Path(directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    save_project(root / "project.vinkulum.json", project)
    write_json(root / "requested-state.json", state)
    engine_files = _engine_fingerprints()
    result = ArticulatedModel(project).analyse(**state)
    if _engine_fingerprints() != engine_files:
        raise ValueError("Pinocchio binary files changed during the calculation.")
    result["engine_files_sha256"] = engine_files
    result["adapter_files_sha256"] = {
        name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in (
            "articulated.py",
            "pinocchio_backend.py",
            "document.py",
            "model.py",
            "cad_data.py",
        )
    }
    result["project_sha256"] = hashlib.sha256(
        (root / "project.vinkulum.json").read_bytes()
    ).hexdigest()
    result["requested_state_sha256"] = hashlib.sha256(
        (root / "requested-state.json").read_bytes()
    ).hexdigest()
    result["environment"] = dict(
        sorted(
            (distribution.metadata["Name"], distribution.version)
            for distribution in metadata.distributions()
        )
    )
    result["worker_python"] = sys.version
    result["platform"] = platform.platform()
    write_json(root / "result.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project", type=Path)
    parser.add_argument("--state", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = run_analysis(
            load_project(args.project),
            args.output,
            None if args.state is None else read_json(args.state),
        )
    except (
        ValueError,
        OSError,
        TypeError,
        ImportError,
        np.linalg.LinAlgError,
    ) as error:
        print(str(error), file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "engine_version": result["engine_version"],
                "nq": result["nq"],
                "nv": result["nv"],
                "scientific_status": result["scientific_status"],
                "output": str(args.output),
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
