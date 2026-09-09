"""Contrôle d'usage d'une roue installée hors du dépôt."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import vinkulum
from vinkulum import _vinkulum


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--readme', type=Path, required=True)
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    root = args.readme.resolve().parent
    module = Path(_vinkulum.__file__).resolve()
    assert not module.is_relative_to(root)
    assert not Path.cwd().resolve().is_relative_to(root)
    code = re.search(r'```python\n(.*?)\n```', args.readme.read_text(), re.S).group(1)
    ns = {}
    exec(compile(code, str(args.readme), 'exec'), ns)
    assert ns['bilan']['statut'] == 'tolerance' and ns['bilan']['strict']
    assert abs(ns['angle']-.02) < 1e-10
    ns['n'].pose_etat(*ns['sauvegarde'])
    assert ns['n'].etat_precis() == ns['sauvegarde']
    out = dict(version=vinkulum.__version__, module=str(module),
               extension_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),
               import_hors_depot=True, exemple_readme_execute=True,
               readme_sha256=hashlib.sha256(args.readme.read_bytes()).hexdigest(),
               angle_readme_rad=ns['angle'], statique_stricte=ns['bilan'],
               etat_precis_restaure=True)
    args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n')
    print('Roue indépendante : version, exemple README et restauration vérifiés.')
