"""Pinocchio-independent admission and coordinate contract for rigid trees."""

import math
from collections import deque
from dataclasses import dataclass

import numpy as np

from .document import Joint, Project, vector


@dataclass(frozen=True)
class TreeLink:
    joint: Joint
    parent: str | None
    child: str
    sign: int
    initial_coordinate: float

    @property
    def coordinate_unit(self):
        return "rad" if self.joint.kind == "pivot" else "m"

    @property
    def effort_unit(self):
        return "N m" if self.joint.kind == "pivot" else "N"

    @property
    def parent_anchor(self):
        return (
            (self.joint.pa, self.joint.ra)
            if self.sign == 1
            else (self.joint.pb, self.joint.rb)
        )

    @property
    def child_anchor(self):
        return (
            (self.joint.pb, self.joint.rb)
            if self.sign == 1
            else (self.joint.pa, self.joint.ra)
        )


def tree_links(project):
    """Orient each declared edge away from ground; preserve its A→B coordinate.

    Ground is None. Multiple branches attached to ground form one fixed-base
    tree. A loop, detached component or unsupported law is rejected explicitly.
    """
    if not isinstance(project, Project):
        raise TypeError("A validated Studio project is required.")
    issues = project.diagnostics()
    if issues:
        raise ValueError(
            "Invalid mechanism: "
            + "; ".join(f"{issue.object_id}: {issue.message}" for issue in issues)
        )
    adjacency = {body.id: [] for body in project.bodies}
    adjacency[None] = []
    for joint in project.joints:
        if joint.kind not in {"pivot", "glissiere"}:
            raise ValueError(
                f"{joint.name}: the first Pinocchio adapter supports revolute and prismatic joints."
            )
        if joint.motion is not None:
            raise ValueError(
                f"{joint.name}: prescribed motion requires a separate constrained-dynamics contract."
            )
        adjacency[joint.a].append(joint)
        adjacency[joint.b].append(joint)
    queue = deque([None])
    seen = {None}
    edges = set()
    links = []
    while queue:
        parent = queue.popleft()
        for joint in sorted(adjacency[parent], key=lambda item: item.id):
            if joint.id in edges:
                continue
            edges.add(joint.id)
            sign = 1 if joint.a == parent else -1
            child = joint.b if sign == 1 else joint.a
            if child in seen:
                raise ValueError(
                    f"{joint.name}: a closed loop is outside the first tree adapter."
                )
            seen.add(child)
            queue.append(child)
            a, A = project.pose(joint.a)
            b, B = project.pose(joint.b)
            QA = A @ np.array(joint.ra).reshape(3, 3)
            QB = B @ np.array(joint.rb).reshape(3, 3)
            relative = QA.T @ QB
            coordinate = (
                math.atan2(relative[1, 0], relative[0, 0])
                if joint.kind == "pivot"
                else float(QA[:, 2] @ (b + B @ joint.pb - a - A @ joint.pa))
            )
            if not math.isfinite(coordinate):
                raise ValueError(f"{joint.name}: nonfinite initial coordinate.")
            links.append(TreeLink(joint, parent, child, sign, coordinate))
    if len(seen) != len(project.bodies) + 1:
        raise ValueError(
            "Every body must be connected to ground; floating or detached components are unsupported."
        )
    return tuple(links)


def state_vector(value, links, name, *, initial=False):
    if value is None:
        return np.array([link.initial_coordinate if initial else 0.0 for link in links])
    return np.array(vector(value, len(links), name))
