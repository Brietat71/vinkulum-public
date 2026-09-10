"""Pure CAD response admission and immutable document preparation; no Qt imports."""

import math
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from .cad_history import CadRecipe
from .document import MAX_PROJECT_BYTES, Body, Project, replace_cad_body
from .execution import ExecutionPlan
from .model import read_json
from .sketch import Sketch


class CadAdmissionError(ValueError):
    def __init__(self, kind, message):
        super().__init__(message)
        self.kind = kind


@dataclass(frozen=True)
class CadAdmission:
    payload: dict
    body: Body | None
    source_project: Project | None
    project: Project | None


def admit_cad_response(
    path, request, input_sha256, plan, code, source_project=None, service_id=None
):
    path = Path(path)
    body = None
    result = read_json(path, MAX_PROJECT_BYTES) if path.is_file() else {}
    if not isinstance(result, dict):
        raise TypeError("Invalid CAD worker response.")
    if result.get("status") == "failed":
        raise CadAdmissionError(
            "operation_failed",
            "CAD operation failed: " + str(result.get("error", "No detail returned.")),
        )
    if code != 0:
        raise CadAdmissionError(
            "worker_failed",
            f"CAD worker exited with code {code}. The previous model is preserved.",
        )
    if (
        result.get("status") != "completed"
        or result.get("request_sha256") != input_sha256
        or (service_id is not None and result.get("service_request_id") != service_id)
    ):
        raise ValueError("CAD response does not match the captured operation.")
    if "body" in result:
        body = Body.from_dict(result["body"])
        if body.cad is None:
            raise ValueError("The CAD worker returned no solid geometry.")
        density = request.get("density", 7800)
        if body.inertia_mode != "homogeneous" or not math.isclose(
            body.mass, body.cad.volume_m3 * density, rel_tol=1e-12, abs_tol=0
        ):
            raise ValueError("CAD mass properties do not match the captured density.")
        operation = request.get("operation")
        if operation == "extrude_sketch":
            features = body.cad.recipe.features if body.cad.recipe else ()
            if (
                len(features) != 1
                or features[0].kind != "extrude_sketch"
                or features[0].profile != Sketch.from_dict(request["profile"])
                or features[0].dimensions_mm != tuple(request["dimensions_mm"])
            ):
                raise ValueError(
                    "The CAD worker changed the captured sketch or extrusion height."
                )
        if (
            operation in ("regenerate", "fillet", "cut", "fuse", "common")
            and body.id != request["a"]["id"]
        ):
            raise ValueError("The CAD worker changed the target body identity.")
        if operation == "regenerate":
            expected = CadRecipe.from_dict(request["recipe"])
            if (
                body.cad.recipe is None
                or replace(body.cad.recipe, origin_in_body_m=expected.origin_in_body_m)
                != expected
            ):
                raise ValueError("The CAD worker returned a different feature graph.")
            source = Body.from_dict(request["a"])
            if body.name != source.name or body.orientation != source.orientation:
                raise ValueError("CAD regeneration changed the captured body frame.")
            expected_position = np.array(source.position) + np.array(
                source.orientation
            ).reshape(3, 3) @ (
                np.array(expected.origin_in_body_m) - body.cad.recipe.origin_in_body_m
            )
            tolerance = 1e-10 + 16 * np.spacing(
                np.maximum(np.abs(expected_position), np.abs(source.position))
            )
            if not np.all(
                np.abs(np.array(body.position) - expected_position) <= tolerance
            ):
                raise ValueError(
                    "CAD regeneration changed the design frame's world origin."
                )
    elif request.get("operation") == "export_step":
        if not isinstance(result.get("step"), str) or not result["step"].startswith(
            "ISO-10303-21;"
        ):
            raise ValueError("Invalid CAD STEP output.")
    elif request.get("operation") != "self_check":
        raise ValueError("The CAD worker returned no part.")
    execution = result.get("execution")
    if (
        not isinstance(execution, dict)
        or ExecutionPlan.from_dict(execution.get("allocation")) != plan
    ):
        raise ValueError("CAD response changed the captured CPU allocation.")
    if execution.get("backend") != "occt-native" or any(
        type(execution.get(key)) is not int or execution[key] != plan.threads
        for key in ("threads", "default_threads")
    ):
        raise ValueError(
            "CAD response does not identify the allocated native OCCT pool."
        )
    prepared = (
        replace_cad_body(source_project, body)
        if source_project is not None and body is not None
        else None
    )
    return CadAdmission(result, body, source_project, prepared)
