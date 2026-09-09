"""Capture les trajectoires du défaut historique pour deux roues isolées."""
import hashlib
import json
from pathlib import Path
import sys

import vinkulum
from vinkulum import _vinkulum
from mesure_sigma import CASES, model

out = {'version': vinkulum.__version__,
       'binaire_sha256': hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
       'modeles': {}}
for name, (end, h) in CASES.items():
    n = model(name)
    # Appel sans le nouvel argument, disponible aussi en 0.7.2.
    tr = n.simule(end, h, rho=.9, tol=1e-12, newton_max=25)
    raw = json.dumps(tr, separators=(',', ':'), allow_nan=False).encode()
    out['modeles'][name] = {'trajectoire_sha256': hashlib.sha256(raw).hexdigest(),
                            'sorties': len(tr), 'stats': list(n.stats())}
Path(sys.argv[1]).write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
