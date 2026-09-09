"""Identité de construction embarquée, indépendante du dépôt."""
from functools import lru_cache
import hashlib
from vinkulum import _vinkulum

@lru_cache(maxsize=1)
def _construction():
    compiler, options, inertia, mass = _vinkulum._construction_certificat()
    return dict(liaison="extension_rust", compilateur=compiler, options=options,
                inertie_source_sha256=hashlib.sha256(inertia.encode()).hexdigest(),
                comparaison_source_sha256=hashlib.sha256(mass.encode()).hexdigest())

def identite():
    return dict(_construction())
