"""Compatibilité des outils de recherche ; implémentation distribuée dans vinkulum."""
import sys
from vinkulum._ports import certificat_spectral as _implementation
sys.modules[__name__] = _implementation
