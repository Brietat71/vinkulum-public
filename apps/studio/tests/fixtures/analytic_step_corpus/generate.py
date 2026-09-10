"""Original analytic AP214 fixtures, v1.0; standard library only.

SPDX-License-Identifier: Apache-2.0
Geometry is written directly from lines, circles and elementary surfaces.
No CAD kernel, exported CAD geometry or numerical-reference module is used.
"""

from fractions import Fraction
from pathlib import Path

ROTATION = ((-10, 2, 11), (10, -5, 10), (5, 14, 2))  # Divide by 15.


def numbers(values):
    return "(" + ",".join(f"{float(value):.17E}" for value in values) + ")"


def vector(a, b):
    return tuple(y - x for x, y in zip(a, b))


def cross(a, b):
    return tuple(
        a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3]
        for i in range(3)
    )


class Writer:
    def __init__(self, name, *, metres, translation):
        self.name, self.translation = name, translation
        self.scale = Fraction(1, 1000) if metres else Fraction(1)
        self.entities, self.points, self.vertices, self.lines = [], {}, {}, {}
        app = self.add("APPLICATION_CONTEXT('automotive design')")
        self.add(
            f"APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,{app})"
        )
        context = self.add(f"PRODUCT_CONTEXT('',{app},'mechanical')")
        product = self.add(
            f"PRODUCT('{name}','{name}','Original Apache-2.0 analytic fixture',({context}))"
        )
        self.add(f"PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,({product}))")
        formation = self.add(
            f"PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('1','',{product},.NOT_KNOWN.)"
        )
        definition_context = self.add(
            f"PRODUCT_DEFINITION_CONTEXT('part definition',{app},'design')"
        )
        definition = self.add(
            f"PRODUCT_DEFINITION('design','',{formation},{definition_context})"
        )
        self.definition_shape = self.add(
            f"PRODUCT_DEFINITION_SHAPE('','',{definition})"
        )
        unit = "$" if metres else ".MILLI."
        length = self.add(f"(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT({unit},.METRE.))")
        angle = self.add("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
        solid_angle = self.add(
            "(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())"
        )
        uncertainty = self.add(
            f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({float(self.scale) * 1e-7:.17E}),{length},'distance_accuracy_value','declared length unit')"
        )
        self.context = self.add(
            f"(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(({uncertainty}))GLOBAL_UNIT_ASSIGNED_CONTEXT(({length},{angle},{solid_angle}))REPRESENTATION_CONTEXT('','3D'))"
        )

    def add(self, value):
        identifier = f"#{len(self.entities) + 1}"
        self.entities.append(f"{identifier}={value};")
        return identifier

    def rotate(self, values):
        return tuple(
            Fraction(sum(a * b for a, b in zip(row, values)), 15) for row in ROTATION
        )

    def point(self, local):
        local = tuple(local)
        if local not in self.points:
            transformed = tuple(
                self.scale * (value + offset)
                for value, offset in zip(self.rotate(local), self.translation)
            )
            self.points[local] = self.add(f"CARTESIAN_POINT('',{numbers(transformed)})")
        return self.points[local]

    def vertex(self, local):
        local = tuple(local)
        if local not in self.vertices:
            self.vertices[local] = self.add(f"VERTEX_POINT('',{self.point(local)})")
        return self.vertices[local]

    def direction(self, local):
        return self.add(f"DIRECTION('',{numbers(self.rotate(local))})")

    def placement(self, origin, normal=(0, 0, 1), x=(1, 0, 0)):
        return self.add(
            f"AXIS2_PLACEMENT_3D('',{self.point(origin)},{self.direction(normal)},{self.direction(x)})"
        )

    def line(self, start, end):
        if (end, start) in self.lines:
            return self.lines[end, start], False
        if (start, end) not in self.lines:
            direction = self.direction(vector(start, end))
            magnitude = self.add(f"VECTOR('',{direction},1.)")
            geometry = self.add(f"LINE('',{self.point(start)},{magnitude})")
            self.lines[start, end] = self.add(
                f"EDGE_CURVE('',{self.vertex(start)},{self.vertex(end)},{geometry},.T.)"
            )
        return self.lines[start, end], True

    def bound(self, edges, *, outer=True, forward=True):
        oriented = [
            self.add(f"ORIENTED_EDGE('',*,*,{edge},{'.T.' if sense else '.F.'})")
            for edge, sense in edges
        ]
        loop = self.add(f"EDGE_LOOP('',({','.join(oriented)}))")
        return self.add(
            f"{'FACE_OUTER_BOUND' if outer else 'FACE_BOUND'}('',{loop},{'.T.' if forward else '.F.'})"
        )

    def face(self, bounds, surface, *, forward=True):
        return self.add(
            f"ADVANCED_FACE('',({','.join(bounds)}),{surface},{'.T.' if forward else '.F.'})"
        )

    def box_shell(self, low, high):
        x0, y0, z0 = low
        x1, y1, z1 = high
        points = (
            (x0, y0, z0),
            (x1, y0, z0),
            (x1, y1, z0),
            (x0, y1, z0),
            (x0, y0, z1),
            (x1, y0, z1),
            (x1, y1, z1),
            (x0, y1, z1),
        )
        loops = [(3, 2, 1, 0), (4, 5, 6, 7)] + [
            (i, (i + 1) % 4, (i + 1) % 4 + 4, i + 4) for i in range(4)
        ]
        faces = []
        for loop in loops:
            vertices = [points[index] for index in loop]
            a, b, c = vertices[:3]
            u, v = vector(a, b), vector(a, c)
            surface = self.add(f"PLANE('',{self.placement(a,cross(u,v),u)})")
            bound = self.bound(
                [
                    self.line(start, end)
                    for start, end in zip(vertices, vertices[1:] + vertices[:1])
                ]
            )
            faces.append(self.face([bound], surface))
        return self.add(f"CLOSED_SHELL('',({','.join(faces)}))")

    def rim(self, x, y, z, radius):
        circle = self.add(
            f"CIRCLE('',{self.placement((x,y,z))},{float(radius*self.scale):.17E})"
        )
        points = tuple(
            (x + dx, y + dy, z)
            for dx, dy in ((radius, 0), (0, radius), (-radius, 0), (0, -radius))
        )
        edges = [
            self.add(
                f"EDGE_CURVE('',{self.vertex(start)},{self.vertex(end)},{circle},.T.)"
            )
            for start, end in zip(points, points[1:] + points[:1])
        ]
        return points, edges

    def cylindrical_wall(self, x, y, radius, height, *, inward=False):
        low, bottom = self.rim(x, y, 0, radius)
        high, top = self.rim(x, y, height, radius)
        upright = [self.line(a, b) for a, b in zip(low, high)]
        surface = self.add(
            f"CYLINDRICAL_SURFACE('',{self.placement((x,y,0))},{float(radius*self.scale):.17E})"
        )
        faces = []
        for i in range(4):
            j = (i + 1) % 4
            edges = [
                (bottom[i], True),
                upright[j],
                (top[i], False),
                (upright[i][0], not upright[i][1]),
            ]
            faces.append(self.face([self.bound(edges)], surface, forward=not inward))
        return faces, bottom, top

    def finish(self, solids):
        shape = self.add(
            f"ADVANCED_BREP_SHAPE_REPRESENTATION('',({','.join(solids)}),{self.context})"
        )
        self.add(f"SHAPE_DEFINITION_REPRESENTATION({self.definition_shape},{shape})")
        return (
            "ISO-10303-21;\nHEADER;\n"
            "FILE_DESCRIPTION(('Original analytic corpus; Apache-2.0'),'2;1');\n"
            f"FILE_NAME('{self.name}.step','2026-09-10T00:00:00',('Vinkulum contributors'),('Vinkulum'),'analytic_step_corpus/generate.py 1.0','Python standard library','Apache-2.0');\n"
            "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n"
            + "\n".join(self.entities)
            + "\nENDSEC;\nEND-ISO-10303-21;\n"
        )


def eccentric_bore():
    writer = Writer("eccentric-bore-metres", metres=True, translation=(300, -140, 220))
    outer, outer_bottom, outer_top = writer.cylindrical_wall(0, 0, 25, 18)
    inner, inner_bottom, inner_top = writer.cylindrical_wall(8, -5, 6, 18, inward=True)
    faces = outer + inner
    for z, normal, outer_edges, inner_edges, forward in (
        (18, (0, 0, 1), outer_top, inner_top, True),
        (0, (0, 0, -1), outer_bottom, inner_bottom, False),
    ):
        surface = writer.add(f"PLANE('',{writer.placement((0,0,z),normal)})")
        bounds = [
            writer.bound([(edge, True) for edge in outer_edges], forward=forward),
            writer.bound(
                [(edge, True) for edge in inner_edges], outer=False, forward=not forward
            ),
        ]
        faces.append(writer.face(bounds, surface))
    shell = writer.add(f"CLOSED_SHELL('',({','.join(faces)}))")
    solid = writer.add(f"MANIFOLD_SOLID_BREP('Eccentric through bore',{shell})")
    return writer.finish([solid])


def offset_void(*, extra_solid=False):
    writer = Writer("offset-closed-void-mm", metres=False, translation=(-80, 25, 41))
    outer = writer.box_shell((0, 0, 0), (60, 44, 30))
    inner = writer.box_shell((25, 8, 10), (47, 26, 22))
    oriented = writer.add(f"ORIENTED_CLOSED_SHELL('',*,{inner},.F.)")
    solid = writer.add(f"BREP_WITH_VOIDS('Offset closed cavity',{outer},({oriented}))")
    solids = [solid]
    if extra_solid:
        other = writer.box_shell((80, 0, 0), (90, 10, 10))
        solids.append(writer.add(f"MANIFOLD_SOLID_BREP('Separate part',{other})"))
    return writer.finish(solids)


if __name__ == "__main__":
    root = Path(__file__).parent
    (root / "eccentric-bore-metres.step").write_text(eccentric_bore(), encoding="ascii")
    (root / "offset-closed-void-mm.step").write_text(offset_void(), encoding="ascii")
