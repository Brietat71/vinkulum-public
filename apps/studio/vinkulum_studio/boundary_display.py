"""Sampled boundary directions; display data never replaces physical conditions."""

from dataclasses import dataclass

import numpy as np

from .mesh_binding import triangle_shape
from .model import finite_number


@dataclass(frozen=True)
class SurfaceFrame:
    point: tuple
    normal: tuple


@dataclass(frozen=True)
class BoundaryGlyph:
    condition: str
    surface: int
    kind: str
    point: tuple
    normal: tuple
    direction: tuple


def unit_direction(values):
    """Normalize finite vectors without squaring their physical magnitude."""
    values = np.asarray(values, dtype=float)
    scale = np.max(np.abs(values))
    if not np.isfinite(values).all() or scale == 0:
        return None
    scaled = values / scale
    return tuple(map(float, scaled / np.linalg.norm(scaled)))


def spread_indices(points, limit):
    """Bounded, deterministic spatial sampling, with no load-weight meaning."""
    if len(points) <= limit:
        return tuple(range(len(points)))
    local = np.asarray(points, dtype=float) - points[0]
    scale = np.max(np.abs(local))
    if scale == 0:
        return tuple(range(limit))
    local /= scale
    index = int(np.argmin(np.sum((local - np.mean(local, axis=0)) ** 2, axis=1)))
    chosen = []
    distance = np.full(len(points), np.inf)
    for _ in range(limit):
        chosen.append(index)
        distance = np.minimum(distance, np.sum((local - local[index]) ** 2, axis=1))
        distance[chosen] = -1
        index = int(np.argmax(distance))
    return tuple(chosen)


def boundary_frames(solid, *, per_face=9):
    """Evaluate P1/P2 faces at triangle centres, using their actual curved map.

    Run with the mesh's other geometric measurements in the background worker.
    A direction that floating-point display arithmetic cannot resolve is omitted;
    this does not change admission of the already captured physical study.
    """
    if type(per_face) is not int or not 1 <= per_face <= 32:
        raise ValueError("Use 1–32 direction samples per face.")
    points = np.asarray(solid.mesh.nodes)
    shapes = {n: triangle_shape(1 / 3, 1 / 3, n == 6) for n in (3, 6)}
    result = {}
    for surface in solid.surfaces:
        candidates = []
        for triangle in surface.triangles:
            p = points[np.asarray(triangle) - 1]
            local = p - p[0]
            scale = np.max(np.abs(local))
            if scale == 0:
                continue
            basis, gradient = shapes[len(triangle)]
            tangent = (local / scale).T @ gradient
            normal = unit_direction(np.cross(tangent[:, 0], tangent[:, 1]))
            if normal is not None:
                candidates.append(
                    SurfaceFrame(tuple(map(float, p[0] + basis @ local)), normal)
                )
        indices = spread_indices([f.point for f in candidates], per_face)
        result[surface.id] = tuple(candidates[i] for i in indices)
    return result


def effective_values(condition, load_factor):
    if condition.kind == "support":
        return ()
    if not finite_number(load_factor):
        return None
    values = tuple(v * load_factor for v in condition.values)
    return values if all(finite_number(v) for v in values) else None


def boundary_glyphs(frames, conditions, load_factor, *, per_condition=24):
    """Directions of physical intent, not nodal-force arrows or force magnitudes.

    Positive pressure is inward. Total-force components and constrained axes
    are expressed in the world frame. Supports do not scale with the load factor.
    """
    glyphs, invalid = [], []
    for condition in conditions:
        values = effective_values(condition, load_factor)
        if values is None:
            invalid.append(condition.id)
            continue
        candidates = [
            (identifier, frame)
            for identifier in condition.surfaces
            for frame in frames.get(identifier, ())
        ]
        # One support anchor per face; up to 24 sampled anchors for a load.
        if condition.kind == "support":
            candidates = [(i, f[0]) for i in condition.surfaces if (f := frames.get(i))]
        indices = spread_indices([f.point for _, f in candidates], per_condition)
        for index in indices:
            identifier, frame = candidates[index]
            if condition.kind == "support":
                directions = [
                    tuple(float(j + 1 == axis) for j in range(3))
                    for axis in condition.axes
                ]
            elif condition.kind == "pressure":
                directions = (
                    [tuple(-np.sign(values[0]) * n for n in frame.normal)]
                    if values[0]
                    else []
                )
            else:
                direction = unit_direction(values)
                directions = [direction] if direction is not None else []
            for direction in directions:
                glyphs.append(
                    BoundaryGlyph(
                        condition.id,
                        identifier,
                        condition.kind,
                        frame.point,
                        frame.normal,
                        direction,
                    )
                )
    return tuple(glyphs), tuple(invalid)
