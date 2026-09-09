"""Compatibilité des contre-épreuves du prototype avec le paquet qualifié."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'python/vinkulum'))
from _verification_assemblage import *

if __name__ == '__main__':
    import json
    bornes = verifier_document(json.loads(Path(sys.argv[1]).read_text()))
    print(json.dumps({'inclusion_stricte': True, 'bornes': list(map(str, bornes))}))
