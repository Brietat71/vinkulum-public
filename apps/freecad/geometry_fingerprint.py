"""Versioned geometry representation identity for native static analyses.

SPDX-License-Identifier: Apache-2.0
Signed-zero geometry tokens are normalized and only the Checked metadata bit is masked in OCCT V1 topology records emitted by
Part. All geometry, tolerances, topology order and other flags remain bytewise.
This is representation identity, not a proof of mathematical shape equality.
"""

import json
import re

import FreeCAD as App
import Part

if __package__:
    from . import bridge
else:
    import bridge

LEGACY = "brep-world-raw-v1"
CURRENT = "brep-v1-zero-checked0-local-world-v1"


def normalize(text):
    if not text.startswith("\nCASCADE Topology V1, (c) Matra-Datavision\n"):
        raise ValueError("Only the observed V1 writer is supported")
    tables = list(re.finditer(r"^TShapes ([0-9]+)\n", text, re.MULTILINE))
    if len(tables) != 1:
        raise ValueError("Expected one topology table")
    table = tables[0]
    body = text[table.end() :]
    tags = list(re.finditer(r"^(?:Ve|Ed|Wi|Fa|Sh|So|CS|Co)\n", body, re.MULTILINE))
    if len(tags) != int(table[1]) or not tags or tags[0].start() != 0:
        raise ValueError("Topology record count mismatch")
    suffix = re.compile(
        r"\n([01]{7})\n((?:[+\-ie][0-9]+ [0-9]+\s+)*\*)\n(?:\n[+\-ie][0-9]+ [0-9]+\s*)?\Z"
    )
    out = []
    for index, tag in enumerate(tags):
        end = tags[index + 1].start() if index + 1 < len(tags) else len(body)
        record = body[tag.start() : end]
        match = suffix.search(record)
        if match is None:
            raise ValueError(f"Unrecognized topology suffix {index}")
        flags = match[1]
        # Numeric -0 is the same real coordinate as 0. Never rewrite signs in
        # topology links: their +/- prefix encodes orientation, not a number.
        geometry = re.sub(r"(?<!\S)-0(?!\S)", "0", record[: match.start(1)])
        out.append(geometry + flags[:2] + "0" + flags[3:] + record[match.end(1) :])
    prefix = re.sub(r"(?<!\S)-0(?!\S)", "0", text[: table.end()])
    return prefix + "".join(out)


def signature(source, kind=CURRENT):
    if kind == LEGACY:
        return bridge.signature(source)
    if kind != CURRENT:
        raise ValueError("Unsupported geometry fingerprint version.")
    shape = Part.Shape(source.Shape)
    shape.Placement = App.Placement()
    shape = shape.cleaned()
    geometry = normalize(shape.exportBrepToString())
    matrix = source.getGlobalPlacement().toMatrix()
    pose = [
        0.0 if getattr(matrix, f"A{i}{j}") == 0 else getattr(matrix, f"A{i}{j}")
        for i in range(1, 5)
        for j in range(1, 5)
    ]
    return bridge.sha(
        json.dumps(
            {"kind": kind, "brep": geometry, "world": pose},
            sort_keys=True,
            allow_nan=False,
        ).encode("ascii")
    )
