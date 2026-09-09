"""Compatibilité des contre-épreuves du prototype avec le paquet qualifié."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'python/vinkulum'))
from _oracle_assemblage import *
