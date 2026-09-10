"""Original planar STEP construction, v1.0; Python standard library only.

SPDX-License-Identifier: Apache-2.0
No CAD kernel, mesh exporter or Vinkulum code generates this fixture.
"""

from fractions import Fraction
from pathlib import Path


def render():
    # Counterclockwise L outline in mm, extruded from z=0 to z=8.
    outline = ((0, 0), (60, 0), (60, 12), (12, 12), (12, 48), (0, 48))
    rotation = ((-3, -4, 12), (12, 3, 4), (-4, 12, 3))  # Divide by 13.
    translation = (17, -23, 31)
    points = [
        tuple(
            Fraction(sum(a * b for a, b in zip(row, (x, y, z))), 13) + t
            for row, t in zip(rotation, translation)
        )
        for z in (0, 8)
        for x, y in outline
    ]
    # Every face loop has outward orientation, including the concave end faces.
    faces = [tuple(reversed(range(6))), tuple(range(6, 12))]
    faces += [(i, (i + 1) % 6, (i + 1) % 6 + 6, i + 6) for i in range(6)]
    entities = []

    def entity(value):
        number = len(entities) + 1
        entities.append(f"#{number}={value};")
        return f"#{number}"

    def numbers(values):
        return "(" + ",".join(f"{float(v):.17E}" for v in values) + ")"

    app = entity("APPLICATION_CONTEXT('automotive design')")
    entity(
        f"APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,{app})"
    )
    context = entity(f"PRODUCT_CONTEXT('',{app},'mechanical')")
    product = entity(
        f"PRODUCT('independent-l-bracket','Independent L bracket','Original Apache-2.0 fixture',({context}))"
    )
    entity(f"PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,({product}))")
    formation = entity(
        f"PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('1','',{product},.NOT_KNOWN.)"
    )
    definition_context = entity(
        f"PRODUCT_DEFINITION_CONTEXT('part definition',{app},'design')"
    )
    definition = entity(
        f"PRODUCT_DEFINITION('design','',{formation},{definition_context})"
    )
    definition_shape = entity(f"PRODUCT_DEFINITION_SHAPE('','',{definition})")
    length = entity("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
    angle = entity("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
    solid_angle = entity("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())")
    uncertainty = entity(
        f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-8),{length},'distance_accuracy_value','mm')"
    )
    geometry_context = entity(
        f"(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(({uncertainty}))"
        f"GLOBAL_UNIT_ASSIGNED_CONTEXT(({length},{angle},{solid_angle}))REPRESENTATION_CONTEXT('','3D'))"
    )
    point_ids = [entity(f"CARTESIAN_POINT('',{numbers(p)})") for p in points]
    face_ids = []
    for face in faces:
        a, b, c = (points[i] for i in face[:3])
        u, v = tuple(y - x for x, y in zip(a, b)), tuple(y - x for x, y in zip(a, c))
        normal = tuple(
            u[(i + 1) % 3] * v[(i + 2) % 3] - u[(i + 2) % 3] * v[(i + 1) % 3]
            for i in range(3)
        )
        axis = entity(f"DIRECTION('',{numbers(normal)})")
        reference = entity(f"DIRECTION('',{numbers(u)})")
        placement = entity(
            f"AXIS2_PLACEMENT_3D('',{point_ids[face[0]]},{axis},{reference})"
        )
        plane = entity(f"PLANE('',{placement})")
        loop = entity(f"POLY_LOOP('',({','.join(point_ids[i] for i in face)}))")
        bound = entity(f"FACE_OUTER_BOUND('',{loop},.T.)")
        face_ids.append(entity(f"FACE_SURFACE('',({bound}),{plane},.T.)"))
    shell = entity(f"CLOSED_SHELL('',({','.join(face_ids)}))")
    solid = entity(f"FACETED_BREP('Independent L bracket',{shell})")
    shape = entity(
        f"FACETED_BREP_SHAPE_REPRESENTATION('',({solid}),{geometry_context})"
    )
    entity(f"SHAPE_DEFINITION_REPRESENTATION({definition_shape},{shape})")
    return (
        "ISO-10303-21;\nHEADER;\n"
        "FILE_DESCRIPTION(('Original planar L bracket; Apache-2.0'),'2;1');\n"
        "FILE_NAME('l-bracket.step','2026-09-10T00:00:00',('Vinkulum contributors'),"
        "('Vinkulum'),'independent_step/generate.py 1.0','Python standard library','Apache-2.0');\n"
        "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n"
        + "\n".join(entities)
        + "\nENDSEC;\nEND-ISO-10303-21;\n"
    )


if __name__ == "__main__":
    Path(__file__).with_name("l-bracket.step").write_text(render(), encoding="ascii")
