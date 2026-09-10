"""Vinkulum — noyau multicorps (Rust) : intégrateur α-généralisé sur SO(3),
liaisons rigides par contraintes, multithread par défaut.

    from vinkulum import Noyau
    N = Noyau()                      # gravité −z par défaut
    i = N.corps("masse", m, J9, r)   # tenseur 3×3 en ligne, axes corps, au CdM
    N.liaison("rotule", None, i, bloque_t=[0, 1, 2], bloque_r=[])
    traj = N.simule(t_end, h)        # (t, r, R, v, w, λ) par pas
"""
import sys as _sys

if _sys.version_info < (3, 14):
    raise ImportError("Vinkulum nécessite Python 3.14 ou plus ; recréez le venv.")

from vinkulum._vinkulum import ExecutionPool, LB_TF, LB_TP, LB_TV, LB_TVL, Noyau

try:                                    # la version vient du paquet installé,
    from importlib.metadata import version as _v   # jamais recopiée ici
    __version__ = _v("vinkulum")
except Exception:                       # noqa: BLE001  (source non installée)
    __version__ = "0.0.0+source"

def depot_docs():
    """Le `docs/` du DÉPÔT, ou None si on tourne depuis un paquet INSTALLÉ.

    `dirname(__file__)/../../docs` suppose l'arborescence des sources
    (`python/vinkulum/`). Depuis un `site-packages`, ce chemin sort de la
    bibliothèque et désigne un répertoire qui n'existe pas — `doc` accusait
    alors la référence d'API d'être « périmée » (elle est simplement absente
    du paquet) et `bancs` écrivait ses résultats DANS la bibliothèque de
    l'utilisateur. Mesuré le 3 sept. en installant la roue dans un venv neuf.

    Un seul marqueur suffit et il ne peut pas se retrouver dans un paquet
    installé : le `Cargo.toml` du noyau.
    """
    import os
    r = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    return r if os.path.exists(os.path.join(r, "Cargo.toml")) else None


__all__ = ["Noyau", "ExecutionPool", "demo", "depot_docs", "__version__", "LB_TP", "LB_TF", "LB_TV", "LB_TVL"]


def demo():
    from vinkulum.maquette import demo as _d
    _d()
