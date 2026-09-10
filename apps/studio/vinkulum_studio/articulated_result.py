"""Read a captured articulated result without loading a numerical engine."""

import hashlib
import math
import re
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from .applied_loads import ancestor_coordinates, world_load_derivative
from .articulated import operator_conventions, state_vector, tree_links
from .document import Law, Project, load_project, rotation
from .model import finite_number, read_json
from .pinocchio_backend import ADAPTER_VERSION, PINOCCHIO_VERSION

MAX_OPERATOR_BYTES = 4 * 1024 * 1024


def _array(value, shape, name):
    raw = np.asarray(value, dtype=object)
    if raw.shape != shape or not all(finite_number(v) for v in raw.flat):
        raise ValueError(
            f"Invalid Pinocchio channel {name}: expected finite values with shape {shape}."
        )
    return np.array(raw, dtype=float)


def _close(actual, expected, name, *, atol=1e-12, rtol=1e-10):
    if (
        not np.isfinite(actual).all()
        or not np.isfinite(expected).all()
        or not np.allclose(actual, expected, rtol=rtol, atol=atol)
    ):
        raise ValueError(f"Inconsistent Pinocchio channel: {name}.")


@dataclass(frozen=True)
class OperatorResult:
    project: Project
    state: dict
    report: dict
    directory: Path
    display_project: Project


def load_operators(directory, *, expected_project=None, expected_state=None):
    """Validate stored data and operator identities, without rerunning Pinocchio.

    This checks internal consistency; it does not authenticate the producer or
    prove the correctness of a finite-precision dynamics implementation.
    """
    root = Path(directory).resolve()
    if root.name == "result.json" and root.is_file():
        root = root.parent
    report = read_json(root / "result.json", MAX_OPERATOR_BYTES)
    if not isinstance(report, dict):
        raise TypeError("A Pinocchio result object is required.")
    schema = report.get("schema")
    versions = {1: "0.1.0", 2: ADAPTER_VERSION}
    if type(schema) is not int or schema not in versions:
        raise ValueError("Unsupported Pinocchio result contract: schema.")
    for key, expected in {
        "format": "vinkulum-pinocchio-operators",
        "engine": "Pinocchio",
        "engine_version": PINOCCHIO_VERSION,
        "adapter_version": versions[schema],
        "scientific_status": "NotAssessed",
    }.items():
        if type(report.get(key)) is not type(expected) or report[key] != expected:
            raise ValueError(f"Unsupported Pinocchio result contract: {key}.")
    project = load_project(root / "project.vinkulum.json")
    requested = read_json(root / "requested-state.json")
    if expected_project is not None and project != expected_project:
        raise ValueError("The worker returned a different captured project.")
    if expected_state is not None and requested != expected_state:
        raise ValueError("The worker returned a different requested state.")
    if not isinstance(requested, dict) or not set(requested) <= {
        "q",
        "velocity",
        "acceleration",
        "effort",
        "time_s",
    }:
        raise ValueError("Invalid requested articulated state.")
    for name, key in (
        ("project.vinkulum.json", "project_sha256"),
        ("requested-state.json", "requested_state_sha256"),
    ):
        if report.get(key) != hashlib.sha256((root / name).read_bytes()).hexdigest():
            raise ValueError(f"Pinocchio captured file fingerprint mismatch: {name}.")
    for key in ("engine_files_sha256", "adapter_files_sha256"):
        fingerprints = report.get(key)
        if (
            not isinstance(fingerprints, dict)
            or not fingerprints
            or not all(
                isinstance(name, str)
                and isinstance(value, str)
                and re.fullmatch(r"[0-9a-f]{64}", value)
                for name, value in fingerprints.items()
            )
        ):
            raise ValueError(f"Missing or invalid provenance: {key}.")
    environment = report.get("environment")
    if (
        not isinstance(environment, dict)
        or environment.get("pin") != PINOCCHIO_VERSION
        or not all(
            isinstance(name, str) and name and isinstance(value, str) and value
            for name, value in environment.items()
        )
        or any(
            not isinstance(report.get(key), str) or not report[key]
            for key in ("worker_python", "platform")
        )
    ):
        raise ValueError("Incomplete Pinocchio worker provenance.")
    links = tree_links(project)
    n = len(links)
    if any(
        type(report.get(key)) is not int or report[key] != n for key in ("nq", "nv")
    ):
        raise ValueError("Pinocchio coordinate dimensions do not match the project.")
    order = report.get("coordinate_order")
    if not isinstance(order, list) or len(order) != n:
        raise ValueError("Missing Pinocchio coordinate identities.")
    for entry, link in zip(order, links):
        expected = {
            "joint_id": link.joint.id,
            "name": link.joint.name,
            "kind": link.joint.kind,
            "parent_body": link.parent,
            "child_body": link.child,
            "coordinate_unit": link.coordinate_unit,
            "effort_unit": link.effort_unit,
            "positive_direction": "B relative to A along joint A local +Z",
        }
        if entry != expected:
            raise ValueError("Pinocchio coordinate identity, sign or unit mismatch.")
    state = report.get("state")
    if not isinstance(state, dict) or set(state) != {
        "q",
        "velocity",
        "acceleration",
        "effort",
        "time_s",
    }:
        raise ValueError("Incomplete Pinocchio state.")
    vectors = {}
    for key in ("q", "velocity", "acceleration", "effort"):
        vectors[key] = _array(state[key], (n,), key)
        expected = state_vector(requested.get(key), links, key, initial=key == "q")
        _close(vectors[key], expected, key, atol=0, rtol=0)
    t = state["time_s"]
    if (
        not finite_number(t)
        or not 0 <= t <= project.duration
        or t != requested.get("time_s", 0.0)
    ):
        raise ValueError("Inconsistent Pinocchio load time.")
    mass = _array(report.get("mass_matrix"), (n, n), "mass matrix")
    _close(mass, mass.T, "mass symmetry", rtol=1e-12, atol=1e-14)
    try:
        np.linalg.cholesky(mass)
    except np.linalg.LinAlgError as error:
        raise ValueError(
            "The Pinocchio mass matrix is not positive definite."
        ) from error
    channels = {
        key: _array(report.get(key), (n,), key)
        for key in (
            "intrinsic_bias",
            "external_effort",
            "inverse_effort",
            "forward_acceleration",
        )
    }
    for acceleration, effort, name in (
        (
            vectors["acceleration"],
            channels["inverse_effort"] + channels["external_effort"],
            "inverse dynamics",
        ),
        (
            channels["forward_acceleration"],
            vectors["effort"] + channels["external_effort"],
            "forward dynamics",
        ),
    ):
        calculated = mass @ acceleration + channels["intrinsic_bias"]
        scale = (
            np.abs(mass) @ np.abs(acceleration)
            + np.abs(channels["intrinsic_bias"])
            + np.abs(effort)
        )
        if (
            not np.isfinite(calculated).all()
            or not np.isfinite(scale).all()
            or not np.all(np.abs(calculated - effort) <= 1e-12 + 1e-10 * scale)
        ):
            raise ValueError(f"Inconsistent Pinocchio {name}.")
    checks = report.get("checks")
    if (
        not isinstance(checks, dict)
        or checks.get("absolute_tolerance_effort_units") != 1e-12
        or checks.get("relative_tolerance") != 1e-10
    ):
        raise ValueError("Invalid captured Pinocchio consistency checks.")
    for key in ("inverse_residual", "forward_effort_residual"):
        _array(checks.get(key), (n,), key)
    derivatives = report.get("intrinsic_inverse_derivatives")
    if not isinstance(derivatives, dict) or set(derivatives) != {
        "q",
        "velocity",
        "acceleration",
    }:
        raise ValueError("Missing intrinsic derivative channels.")
    for key, value in derivatives.items():
        matrix = _array(value, (n, n), f"intrinsic derivative/{key}")
        if key == "acceleration":
            _close(matrix, mass, "acceleration derivative")
    load_derivatives = {}
    for channel in ("external_effort_derivatives", "loaded_inverse_derivatives"):
        if schema == 1:
            if channel in report:
                raise ValueError("Schema 1 does not define applied-load derivatives.")
            continue
        values = report.get(channel)
        if not isinstance(values, dict) or set(values) != set(derivatives):
            raise ValueError(f"Missing {channel} channels.")
        load_derivatives[channel] = {
            key: _array(value, (n, n), f"{channel}/{key}")
            for key, value in values.items()
        }
    if schema == 2:
        external_derivatives = load_derivatives["external_effort_derivatives"]
        for key in derivatives:
            _close(
                load_derivatives["loaded_inverse_derivatives"][key],
                np.array(derivatives[key]) - external_derivatives[key],
                f"loaded derivative/{key}",
            )
            if key != "q":
                _close(
                    external_derivatives[key],
                    np.zeros((n, n)),
                    f"external derivative/{key}",
                    atol=0,
                    rtol=0,
                )
    captured_bodies = report.get("bodies")
    if not isinstance(captured_bodies, list) or len(captured_bodies) != len(
        project.bodies
    ):
        raise ValueError("Missing captured body values.")
    displayed, jacobians, poses = [], {}, {}
    for body, value in zip(project.bodies, captured_bodies):
        if (
            not isinstance(value, dict)
            or value.get("id") != body.id
            or value.get("name") != body.name
        ):
            raise ValueError(
                "Pinocchio body identities do not match the captured model."
            )
        position = _array(value.get("position_m"), (3,), "body position")
        R = _array(value.get("body_to_world"), (3, 3), "body rotation")
        rotation(tuple(R.flat))
        J = _array(value.get("jacobian"), (6, n), "body Jacobian")
        _close(
            _array(value.get("velocity_world"), (6,), "body velocity"),
            J @ vectors["velocity"],
            "body velocity",
        )
        jacobians[body.id], poses[body.id] = J, (position, R)
        displayed.append(
            replace(body, position=tuple(position), orientation=tuple(R.flat))
        )
    display_project = replace(
        project,
        bodies=tuple(displayed),
        loads=tuple(
            replace(
                load,
                force=tuple(Law(values=(law.value(t),)) for law in load.force),
                moment=tuple(Law(values=(law.value(t),)) for law in load.moment),
            )
            for load in project.loads
        ),
    )
    # Recheck joint compatibility and represented coordinates from returned poses.
    for i, link in enumerate(tree_links(display_project)):
        error = link.initial_coordinate - vectors["q"][i]
        if link.joint.kind == "pivot":
            error = math.remainder(error, 2 * math.pi)
        if abs(error) > 1e-8:
            raise ValueError(
                "Returned body poses do not represent the requested joint coordinates."
            )
    load_values = report.get("loads")
    if not isinstance(load_values, list) or len(load_values) != len(project.loads):
        raise ValueError("Missing captured load values.")
    external = np.zeros(n)
    external_dq = np.zeros((n, n))
    paths = ancestor_coordinates(links)
    for load, value in zip(project.loads, load_values):
        if (
            not isinstance(value, dict)
            or value.get("id") != load.id
            or value.get("body") != load.body
        ):
            raise ValueError("Captured load identity mismatch.")
        position, R = poses[load.body]
        J = jacobians[load.body]
        point_J = J[:3] + np.cross(J[3:].T, R @ load.point).T
        force = np.array([law.value(t) for law in load.force])
        moment = np.array([law.value(t) for law in load.moment])
        for key, expected in (
            ("point_world_m", position + R @ load.point),
            ("force_world_N", force),
            ("moment_world_Nm", moment),
        ):
            _close(_array(value.get(key), (3,), key), expected, key)
        external += point_J.T @ force + J[3:].T @ moment
        if schema == 2:
            external_dq += world_load_derivative(
                point_J, J[3:], force, moment, paths[load.body]
            )
    _close(channels["external_effort"], external, "applied world loads")
    if schema == 2:
        _close(
            external_derivatives["q"], external_dq, "applied world load derivative/q"
        )
    expected_energy = {
        "kinetic_energy_J": float(
            0.5 * vectors["velocity"] @ mass @ vectors["velocity"]
        ),
        "potential_energy_J": float(
            -sum(
                body.mass * np.array(project.gravity) @ poses[body.id][0]
                for body in project.bodies
            )
        ),
    }
    for key, value in expected_energy.items():
        if not finite_number(report.get(key)):
            raise ValueError(f"Missing or invalid {key}.")
        _close(report[key], value, key)
    if report.get("conventions") != operator_conventions(schema):
        raise ValueError("Unsupported operator units, frames or derivative scope.")
    return OperatorResult(project, state, report, root, display_project)
