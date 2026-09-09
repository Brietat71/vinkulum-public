"""Small editable mechanical examples; all physics goes through the public API."""

import math
from dataclasses import replace

from .document import Body, Law, Load, Project, joint_at, new_id, pendulum


def ry(angle):
    c, s = math.cos(angle), math.sin(angle)
    return (c, 0.0, s, 0.0, 1.0, 0.0, -s, 0.0, c)


def rz(angle):
    c, s = math.cos(angle), math.sin(angle)
    return (c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0)


def double_pendulum():
    a, b = math.radians(30), math.radians(-15)
    end = (math.sin(a), 0.0, -math.cos(a))
    first = Body(
        new_id(),
        "Upper arm",
        dimensions=(0.06, 0.06, 1.0),
        mass=1.0,
        position=tuple(v / 2 for v in end),
        orientation=ry(-a),
    )
    second = Body(
        new_id(),
        "Lower arm",
        dimensions=(0.06, 0.06, 1.0),
        mass=1.0,
        position=(end[0] + math.sin(b) / 2, 0.0, end[2] - math.cos(b) / 2),
        orientation=ry(-b),
    )
    p = Project(
        new_id(),
        "Rigid double pendulum",
        bodies=(first, second),
        duration=2.0,
        step=0.005,
    )
    p = p.replace_object(joint_at(p, "pivot", None, first.id, axis=(0.0, 1.0, 0.0)))
    p = p.replace_object(
        joint_at(p, "pivot", first.id, second.id, point=end, axis=(0.0, 1.0, 0.0))
    )
    return p


def slider_crank():
    angle = 0.5
    crank_length = 0.3
    rod_length = 0.8
    pin = (crank_length * math.cos(angle), crank_length * math.sin(angle), 0.0)
    x = pin[0] + math.sqrt(rod_length**2 - pin[1] ** 2)
    rod_angle = math.atan2(-pin[1], x - pin[0])
    crank = Body(
        new_id(),
        "Crank",
        dimensions=(crank_length, 0.04, 0.04),
        mass=0.5,
        position=(pin[0] / 2, pin[1] / 2, 0.0),
        orientation=rz(angle),
    )
    rod = Body(
        new_id(),
        "Connecting rod",
        dimensions=(rod_length, 0.04, 0.04),
        mass=0.4,
        position=((pin[0] + x) / 2, pin[1] / 2, 0.0),
        orientation=rz(rod_angle),
    )
    slider = Body(
        new_id(),
        "Slider",
        dimensions=(0.14, 0.16, 0.12),
        mass=0.3,
        position=(x, 0.0, 0.0),
    )
    p = Project(
        new_id(),
        "Driven slider-crank",
        bodies=(crank, rod, slider),
        gravity=(0.0, 0.0, 0.0),
        duration=1.0,
        step=0.005,
    )
    motor = joint_at(p, "pivot", None, crank.id)
    p = p.replace_object(
        replace(motor, name="Driven revolute joint", motion=Law("lineaire", (0.0, 2.0)))
    )
    p = p.replace_object(joint_at(p, "pivot", crank.id, rod.id, point=pin))
    p = p.replace_object(joint_at(p, "pivot", rod.id, slider.id, point=(x, 0.0, 0.0)))
    p = p.replace_object(
        joint_at(
            p, "glissiere", None, slider.id, point=(x, 0.0, 0.0), axis=(1.0, 0.0, 0.0)
        )
    )
    return p


def temporal_force():
    b = Body(new_id(), "Driven mass", dimensions=(0.2, 0.2, 0.2), mass=2.0)
    p = Project(
        new_id(),
        "Increasing force",
        bodies=(b,),
        gravity=(0.0, 0.0, 0.0),
        duration=1.0,
        step=0.005,
    )
    return p.replace_object(
        Load(
            new_id(),
            "Linear X force",
            b.id,
            force=(Law("lineaire", (0.0, 6.0)), Law(), Law()),
        )
    )


EXAMPLES = {
    "G0 pendulum": pendulum,
    "Double pendulum": double_pendulum,
    "Slider-crank": slider_crank,
    "Time-dependent force": temporal_force,
}
