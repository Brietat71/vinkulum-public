"""Compatibilité des outils de recherche ; implémentation distribuée dans vinkulum."""
import sys
from vinkulum._ports import ports_releves as _implementation
sys.modules[__name__] = _implementation
