"""Compatibilité des outils de recherche ; implémentation distribuée dans vinkulum."""
import sys
from vinkulum._ports import condensation_energie_lu as _implementation
sys.modules[__name__] = _implementation
