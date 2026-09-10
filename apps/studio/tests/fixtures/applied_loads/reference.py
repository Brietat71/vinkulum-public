"""Original rotation/virtual-work references for fixed-world applied loads."""

import math
from uuid import NAMESPACE_URL, uuid5
from dataclasses import replace

import numpy as np
from vinkulum_studio.document import Body, Law, Load, Project, joint_at, new_id


def fixture_id(name):
    return str(uuid5(NAMESPACE_URL, "vinkulum-applied-loads/v1/" + name))


def spatial_pair(
    *,
    mixed=False,
    centre=(0.15, -0.12, 0.2),
    point=(0.11, 0.21, 0.34),
    force=(2, -3, 5),
    moment=(7, 11, 13),
):
    first = Body(fixture_id("first"), "First hinge", dimensions=(0.1, 0.1, 0.1))
    second = Body(
        fixture_id("second"), "Loaded body", position=centre, dimensions=(0.1, 0.2, 0.3)
    )
    project = Project(
        fixture_id("mixed" if mixed else "spatial"), bodies=(first, second), duration=2
    )
    joints = (
        joint_at(
            project, "pivot", None, first.id, axis=(0, 0, 1) if mixed else (1, 0, 0)
        ),
        joint_at(
            project,
            "glissiere" if mixed else "pivot",
            first.id,
            second.id,
            axis=(1, 0, 0) if mixed else (0, 1, 0),
        ),
    )
    joints = tuple(
        replace(joint, id=fixture_id(f"joint-{i}")) for i, joint in enumerate(joints)
    )
    load = Load(
        fixture_id("load"),
        "World wrench",
        second.id,
        point=point,
        force=tuple(Law(values=(value,)) for value in force),
        moment=tuple(Law(values=(value,)) for value in moment),
    )
    return replace(project, joints=joints, loads=(load,))


def rotation_reference(q, arm, force, moment, *, mixed=False):
    """Differentiate explicit Rx Ry or Rz/translation matrices; never use a solver J."""
    x, y = q
    c, s, C, S = math.cos(x), math.sin(x), math.cos(y), math.sin(y)
    if mixed:
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        dR = np.array([[-s, -c, 0], [c, -s, 0], [0, 0, 0]])
        ddR = np.array([[-c, s, 0], [-s, -c, 0], [0, 0, 0]])
        u = np.array(arm) + np.array([y, 0, 0])
        position = R @ u
        J = np.column_stack((dR @ u, R[:, 0]))
        H = np.zeros((3, 2, 2))
        H[:, 0, 0], H[:, 0, 1], H[:, 1, 0] = ddR @ u, dR[:, 0], dR[:, 0]
        W = np.array([[0, 0], [0, 0], [1, 0]])
        dW = np.zeros((3, 2, 2))
    else:
        X = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
        dX = np.array([[0, 0, 0], [0, -s, -c], [0, c, -s]])
        ddX = np.array([[0, 0, 0], [0, -c, s], [0, -s, -c]])
        Y = np.array([[C, 0, S], [0, 1, 0], [-S, 0, C]])
        dY = np.array([[-S, 0, C], [0, 0, 0], [-C, 0, -S]])
        ddY = np.array([[-C, 0, -S], [0, 0, 0], [S, 0, -C]])
        position = X @ Y @ arm
        J = np.column_stack((dX @ Y @ arm, X @ dY @ arm))
        H = np.empty((3, 2, 2))
        H[:, 0, 0], H[:, 1, 1] = ddX @ Y @ arm, X @ ddY @ arm
        H[:, 0, 1] = H[:, 1, 0] = dX @ dY @ arm
        W = np.array([[1, 0], [0, c], [0, s]])
        dW = np.zeros((3, 2, 2))
        dW[:, 1, 0] = [0, -s, c]
    return (
        position,
        J,
        W,
        H,
        dW,
        J.T @ force + W.T @ moment,
        np.einsum("kij,k->ij", H, force) + np.einsum("kij,k->ij", dW, moment),
    )


def move_world(project, R, t):
    def rotate_laws(laws):
        # The fixtures use constants or affine time laws. Rotate both coefficients.
        samples = np.array([[law.value(time) for law in laws] for time in (0, 1)])
        moved = samples @ R.T
        return tuple(Law("lineaire", (float(a), float(b - a))) for a, b in zip(*moved))

    joints = []
    for joint in project.joints:
        data = {}
        for side in ("a", "b"):
            if getattr(joint, side) is None:
                data["p" + side] = tuple(R @ getattr(joint, "p" + side) + t)
                data["r" + side] = tuple(
                    (R @ np.array(getattr(joint, "r" + side)).reshape(3, 3)).flat
                )
        joints.append(replace(joint, **data))
    return replace(
        project,
        bodies=tuple(
            replace(
                body,
                position=tuple(R @ body.position + t),
                orientation=tuple((R @ np.array(body.orientation).reshape(3, 3)).flat),
            )
            for body in project.bodies
        ),
        joints=tuple(joints),
        gravity=tuple(R @ project.gravity),
        loads=tuple(
            replace(
                load, force=rotate_laws(load.force), moment=rotate_laws(load.moment)
            )
            for load in project.loads
        ),
    )


def reframe_bodies(project, S):
    """Rotate body coordinates about each COM; preserve physical inertia and point."""
    return replace(
        project,
        bodies=tuple(
            replace(
                body,
                orientation=tuple((np.array(body.orientation).reshape(3, 3) @ S).flat),
                inertia_mode="explicit",
                explicit_inertia=tuple(
                    (S.T @ np.array(body.inertia()).reshape(3, 3) @ S).flat
                ),
            )
            for body in project.bodies
        ),
        joints=tuple(
            replace(
                joint,
                **{
                    key
                    + side: (
                        tuple(S.T @ getattr(joint, key + side))
                        if key == "p"
                        else tuple(
                            (
                                S.T @ np.array(getattr(joint, key + side)).reshape(3, 3)
                            ).flat
                        )
                    )
                    for side in ("a", "b")
                    if getattr(joint, side) is not None
                    for key in ("p", "r")
                },
            )
            for joint in project.joints
        ),
        loads=tuple(
            replace(load, point=tuple(S.T @ load.point)) for load in project.loads
        ),
    )
