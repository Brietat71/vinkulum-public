"""Un corps immobile contraint selon deux directions presque parallèles."""
import hashlib
import json
import math
from pathlib import Path
from vinkulum import Noyau, _vinkulum


def modele(delta):
    n = Noyau([0., 0., 0.])
    n.corps('support', 1., [1., 0., 0., 0., 1., 0., 0., 0., 1.], [0., 0., 0.])
    for name, angle in [('x', 0.), ('proche', math.atan(delta))]:
        c, s = math.cos(angle), math.sin(angle)
        n.liaison(name, None, 0, bloque_t=[0], bloque_r=[],
                  ra=[c, -s, 0., s, c, 0., 0., 0., 1.])
    n.liaison('reste', None, 0, bloque_t=[2], bloque_r=[0, 1, 2])
    n.effort(0, [2., 1., 3.], [0., 0., 0.])
    return n


def main():
    output = dict(extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
                  programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), cas=[])
    for delta in [1., 1e-3, 1e-6, 1e-9, 1e-12, 0.]:
        n = modele(delta)
        initial = n.etat_precis()
        erreur = None
        try:
            n.statique(strict=True, tol=1e-8, iters=10, paliers_max=1)
        except Exception as e:
            erreur = str(e)
        try:
            reactions = n.reactions()
        except RuntimeError:
            reactions = None
        output['cas'].append(dict(delta=delta, erreur=erreur, rapport=n.statique_info(),
                                 reactions=reactions, etat_initial=initial, etat_final=n.etat_precis()))
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
