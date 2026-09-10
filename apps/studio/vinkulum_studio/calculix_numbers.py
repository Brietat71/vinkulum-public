"""Bounded decimal transport for CalculiX's ``read(..., '(f20.0)')``.

The reader in nodes.f, cloads.f and elastics.f consumes only 20 characters.
Fortran permits an exponent without E, so ``123-5`` denotes 0.00123 here.
Prefer an exact binary64 round trip. Otherwise retain at least 15 significant
digits and enforce a 64-ULP bound, including subnormals. The caller must capture
the requested inputs and revalidate the transported geometry and material.
"""

import math
import re
from decimal import Context, Decimal, localcontext
from fractions import Fraction
from functools import lru_cache


def read_number_field(field):
    """Decode our restricted F20.0 grammar; never emulate silent truncation."""
    if not isinstance(field, str) or not 1 <= len(field) <= 20:
        raise ValueError("CalculiX numeric fields must fit within 20 characters.")
    match = re.fullmatch(r"(-?(?:\d+(?:\.\d*)?|\.\d+))(?:(?:[eE])?([+-]?\d+))?", field)
    if not match:
        raise ValueError("Invalid CalculiX numeric field.")
    return float(match[1] + ("e" + match[2] if match[2] else ""))


def _compact(decimal):
    sign, digits, exponent = decimal.as_tuple()
    digits = "".join(map(str, digits))
    while len(digits) > 1 and digits.endswith("0"):
        digits, exponent = digits[:-1], exponent + 1
    prefix = "-" if sign else ""
    candidates = []
    fixed = format(decimal, "f")
    fixed = re.sub(r"^(-?)0\.", r"\1.", fixed)
    if len(fixed) <= 20:
        candidates.append(fixed)
    for split in range(len(digits) + 1):
        power = exponent + len(digits) - split
        mantissa = digits[:split] + (
            "." + digits[split:] if split < len(digits) else ""
        )
        text = prefix + mantissa + (f"{power:+d}" if power else "")
        if len(text) <= 20:
            candidates.append(text)
    return min(candidates, key=lambda text: (len(text), text)) if candidates else None


def number_field(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("CalculiX inputs must be finite.")
    if value == 0:
        # Python's numeric cache keys equate +0.0 and -0.0.
        return "-0" if math.copysign(1.0, value) < 0 else "0"
    return _nonzero_field(value)


@lru_cache(maxsize=32768)
def _nonzero_field(value):
    original = format(value, ".17g")
    if len(original) <= 20:
        return original
    exact = _compact(Decimal(repr(value)))
    if exact is not None and read_number_field(exact) == value:
        return exact
    candidates = []
    for precision in (17, 16, 15):
        rounded = Decimal(format(value, f".{precision}g"))
        # Include adjacent decimals: nearest decimal may overflow at DBL_MAX.
        with localcontext(Context(prec=32)):
            quantum = Decimal(1).scaleb(rounded.as_tuple().exponent)
            neighbours = (rounded, rounded - quantum, rounded + quantum)
        for decimal in neighbours:
            field = _compact(decimal)
            if field is None:
                continue
            actual = read_number_field(field)
            if math.isfinite(actual):
                error = abs(Fraction(actual) - Fraction(value))
                candidates.append((error, len(field), field))
    if not candidates:
        raise ValueError("No bounded CalculiX decimal representation is available.")
    error, _, field = min(candidates)
    if error > 64 * Fraction(math.ulp(value)):
        raise ValueError("CalculiX input conversion exceeds its 64-ULP budget.")
    return field


def transported_number(value):
    return read_number_field(number_field(value))
