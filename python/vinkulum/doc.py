"""RÉFÉRENCE D'API GÉNÉRÉE — `docs/API.md` sort du CODE, jamais de la main.

Une documentation recopiée diverge : c'est la règle que ce projet a payée dix
fois côté FRELON, et il n'y a aucune raison qu'une API y échappe. Ce module
introspecte le noyau compilé et les modules Python, et rend le Markdown.

Deux modes, comme `cad/etat.py` chez FRELON :

    python -m vinkulum.doc            CONTRÔLE — rouge si docs/API.md a dérivé
    python -m vinkulum.doc --ecrire   régénère

Le mode par défaut est le CONTRÔLE : la CI ne répare pas le fichier qu'elle
juge, elle rougit et donne la commande.
"""
import importlib
import inspect
import os
import sys

__all__ = ["rendu", "controle", "ecrire"]

_MODULES = ["certification", "trim", "rotor", "floquet", "reduction", "reduction_ports", "reduction_contrainte", "sensibilite", "convergence",
            "adjoint_temps", "rotation", "contact", "coque", "andrews",
            "verification", "bancs"]


def _resume(doc):
    """Première phrase utile d'un docstring, sur une ligne."""
    if not doc:
        return ""
    for ligne in doc.strip().splitlines():
        l = ligne.strip()
        if l:
            return l.rstrip(".")
    return ""


def rendu():
    """Le Markdown complet de la référence."""
    import vinkulum
    from vinkulum import _vinkulum

    out = ["# Référence d'API — vinkulum",
           "",
           "**Ce fichier est GÉNÉRÉ** par `python -m vinkulum.doc --ecrire`.",
           "Toute édition à la main est écrasée et fait rougir la CI — une",
           "documentation recopiée finit par décrire une autre bibliothèque.",
           "",
           f"Version : `{getattr(vinkulum, '__version__', '?')}`",
           "",
           "## Le noyau — `vinkulum.Noyau`",
           "",
           _resume(_vinkulum.Noyau.__doc__) or "Le modèle multicorps et son intégrateur.",
           ""]

    membres = [(n, m) for n, m in inspect.getmembers(_vinkulum.Noyau)
               if not n.startswith("_") and callable(m)]
    out += ["| méthode | ce qu'elle fait |", "|---|---|"]
    for nom, m in sorted(membres):
        out.append(f"| `{nom}` | {_resume(m.__doc__)} |")
    out += ["", "## Les modules d'analyse", ""]
    for nm in _MODULES:
        try:
            mod = importlib.import_module(f"vinkulum.{nm}")
        except Exception as e:                      # noqa: BLE001
            out += [f"### `vinkulum.{nm}`", "", f"*non importable : {e}*", ""]
            continue
        out += [f"### `vinkulum.{nm}`", "", _resume(mod.__doc__), ""]
        fns = [(n, f) for n, f in inspect.getmembers(mod, inspect.isfunction)
               if not n.startswith("_") and f.__module__ == mod.__name__]
        if fns:
            out += ["| fonction | ce qu'elle fait |", "|---|---|"]
            for n, f in sorted(fns):
                try:
                    sig = str(inspect.signature(f))
                except (ValueError, TypeError):
                    sig = "(…)"
                out.append(f"| `{n}{sig}` | {_resume(f.__doc__)} |")
            out.append("")
    return "\n".join(out).rstrip() + "\n"


def _chemin():
    from vinkulum import depot_docs
    r = depot_docs()
    return os.path.join(r, "docs", "API.md") if r else None


def ecrire():
    p = _chemin()
    if p is None:
        raise SystemExit("vinkulum.doc --ecrire vise le DÉPÔT ; ici le paquet "
                         "est installé (pas de sources). Lancer depuis un clone.")
    p = os.path.abspath(p)
    with open(p, "w", encoding="utf-8") as f:
        f.write(rendu())
    print(f"→ {p}")
    return p


def controle():
    """Rouge si `docs/API.md` ne correspond plus au code."""
    p = _chemin()
    if p is None:
        # paquet INSTALLÉ : la référence n'y est pas, et ce n'est pas une
        # dérive. Accuser ici ferait crier au loup chez l'utilisateur.
        print(f"╚═ référence d'API non embarquée dans le paquet installé — "
              f"ce contrôle vise le dépôt ({rendu().count(chr(10) + '| `')} "
              f"entrées générées depuis le code)")
        return 0
    p = os.path.abspath(p)
    neuf = rendu()
    vieux = open(p, encoding="utf-8").read() if os.path.exists(p) else ""
    if vieux != neuf:
        n_meth = neuf.count("\n| `")
        print(f"╔═ RÉFÉRENCE D'API PÉRIMÉE — {n_meth} entrées attendues")
        print(f"║ {p}")
        print("║ ne correspond plus au code. Régénérer :")
        print("║     python -m vinkulum.doc --ecrire")
        print("╚═ (la CI ne répare pas le fichier qu'elle juge)")
        return 1
    print(f"╚═ référence d'API à jour ({neuf.count(chr(10) + '| `')} entrées)")
    return 0


if __name__ == "__main__":
    if "--ecrire" in sys.argv:
        ecrire()
    else:
        sys.exit(controle())
