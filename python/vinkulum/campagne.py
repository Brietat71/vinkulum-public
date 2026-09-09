"""CAMPAGNE DE NON-RÉGRESSION DU SOLVEUR STATIQUE.

Le 3 sept., trois réparations de `statique` ont été essayées sur la frontière
erratique des structures à précontrainte seule. Chacune faisait passer des cas
et en cassait d'autres, et il a fallu conclure « l'arbitrage demande une
campagne de non-régression ». Cette campagne-là n'était pas une décision à
prendre : c'était du code à écrire. Le voici.

Ce que ça donne, et qui manquait : **un changement du solveur se juge sur une
table, pas sur le cas qu'on avait sous la main.** Chaque ligne porte son
verdict et, quand une référence existe, son écart.

L'axe du corpus est la RAIDEUR, parce que c'est là que vit le défaut :

    précontrainte SEULE   chaîne de barres — aucune raideur matérielle
    raideur matérielle    la MÊME chaîne en poutres
    raideur + boucles     hexapode (six boucles, zéro degré libre)
    grands déplacements   console en flexion

La paire chaîne-barres / chaîne-poutres est le contrôle **falsifiant** : si la
seconde échouait comme la première, le diagnostic « c'est la précontrainte
seule » serait faux et devrait être rejoué.

    python -m vinkulum.campagne            la table
    python -m vinkulum.campagne --json     de quoi comparer deux versions
"""
import json
import sys
import time

import numpy as np

from vinkulum import Noyau

# ─── les cas ────────────────────────────────────────────────────────────────
PORTEE, LONGUEUR, M_LIN = 8.0, 10.0, 3.0


def _noeuds(n):
    prof = np.sqrt((LONGUEUR / 2) ** 2 - (PORTEE / 2) ** 2)
    return np.array([[PORTEE * k / n, 0.0, -prof * (1 - abs(2.0 * k / n - 1.0))]
                     for k in range(n + 1)])


def _chaine(n, poutres):
    """La même chaîne, tenue par des BARRES ou par des POUTRES.

    Les deux portent le même poids sur la même portée depuis le même V. Seule
    la nature de la raideur change — et c'est la variable du corpus.
    """
    ell, nd = LONGUEUR / n, _noeuds(n)
    mn = M_LIN * ell
    nn = Noyau([0.0, 0.0, -9.81])
    ids = [nn.corps(f"nd{k}", mn, [1e-9, 0, 0, 0, 1e-9, 0, 0, 0, 1e-9], list(nd[k]))
           for k in range(1, n)]
    for i in range(len(ids)):
        nn.liaison(f"pl{i}", None, ids[i], bloque_t=[1], bloque_r=[0, 1, 2])
    if poutres:
        # ⚠ LE JUMEAU DOIT RÉSOUDRE LE MÊME PROBLÈME, et c'est plus subtil
        # qu'il n'y paraît. Premier jet : poutres créées SUR LE V. Une poutre
        # prend sa configuration de référence à la pose de création, donc elle
        # y est NON DÉFORMÉE — le V devient un équilibre à EI près, et la
        # chaîne « convergeait » à sa flèche de départ sans avoir rien résolu.
        # Un contrôle qui ne bouge pas ne contrôle rien (l'assert de vacuité
        # l'a sorti). Ici les poutres naissent DROITES, sur toute la longueur
        # de la ligne : l'ancrage à `PORTEE` raccourcit la portée de 2 m et
        # force le fléchissement, comme un vrai câble.
        ell = LONGUEUR / n
        nn = Noyau([0.0, 0.0, -9.81])
        # DÉFAUT INITIAL, sans quoi rien ne fléchit : une ligne parfaitement
        # droite sous compression est un équilibre (instable), et la statique
        # s'y arrête — mesuré, flèche 3,7 mm au lieu de 2,65 m. C'est le cas
        # d'école du flambement, et il vaut ici comme partout : un solveur
        # statique ne bifurque pas tout seul.
        ids = [nn.corps(f"nd{k}", mn, [1e-9, 0, 0, 0, 1e-9, 0, 0, 0, 1e-9],
                        [k * ell, 0.0,
                         -0.05 * LONGUEUR * np.sin(np.pi * k / n)])
               for k in range(1, n)]
        for i2 in range(len(ids)):
            nn.liaison(f"pl{i2}", None, ids[i2], bloque_t=[1], bloque_r=[0, 1, 2])
        a0 = nn.corps("a0", 1e-9, [1e-12, 0, 0, 0, 1e-12, 0, 0, 0, 1e-12],
                      [0.0, 0.0, 0.0])
        a1 = nn.corps("a1", 1e-9, [1e-12, 0, 0, 0, 1e-12, 0, 0, 0, 1e-12],
                      [n * ell, 0.0, 0.0])
        ch = [a0] + ids + [a1]
        for k in range(len(ch) - 1):
            # (nom, a, b, EA, GA, GJ, EI) — un CÂBLE : quasi inextensible,
            # flexion négligeable. La raideur vient de la MATIÈRE.
            nn.poutre(f"p{k}", ch[k], ch[k + 1], 1e5, 1e5, 1e-6, 1e-6)
        nn.liaison("ea", None, a0)
        # l'ancrage aval RAMÈNE le bout de 10 m à 8 m : c'est lui qui fait
        # fléchir, et il est violé de 2 m au départ (l'assemblage s'en charge)
        nn.liaison("eb", None, a1, pa=[PORTEE, 0.0, 0.0])
        nn.assemble()
        return nn, ids
    else:
        for k in range(n):
            if k == 0:
                nn.distance("b0", None, ids[0], list(nd[0]), [0, 0, 0], ell)
            elif k == n - 1:
                nn.distance(f"b{k}", ids[k - 1], None, [0, 0, 0], list(nd[n]), ell)
            else:
                nn.distance(f"b{k}", ids[k - 1], ids[k], [0, 0, 0], [0, 0, 0], ell)
    return nn, ids


def _fleche(nn, ids):
    e = nn.etat()
    return float(-min(e[1][i][2] for i in ids))


def _cas_chaine(n, poutres):
    def f():
        nn, ids = _chaine(n, poutres)
        nn.statique()
        return dict(valeur=_fleche(nn, ids))
    return f


def _console(ne=8):
    """Console en flexion — raideur matérielle, grands déplacements."""
    def f():
        nn = Noyau([0.0, 0.0, 0.0])
        ix = [nn.corps(f"c{i}", 0.05, [1e-6, 0, 0, 0, 1e-6, 0, 0, 0, 1e-6],
                       [i * 1.0 / ne, 0.0, 0.0]) for i in range(ne + 1)]
        for i in range(ne):
            nn.poutre(f"p{i}", ix[i], ix[i + 1], 1e6, 5e5, 50.0, 100.0)
        nn.liaison("enc", None, ix[0], pa=[0, 0, 0])
        nn.effort(ix[-1], [0.0, 0.0, -30.0], [0.0, 0.0, 0.0])
        nn.statique(tol=1e-10, iters=100)
        return dict(valeur=float(nn.etat()[1][ix[-1]][2]))
    return f


def _hexapode():
    def f():
        from vinkulum import domaines
        r = domaines.hexapode(bavard=False)
        return dict(valeur=r["e_f"])
    return f


def cas():
    """Le corpus. Nom, fonction, et la référence quand il en existe une."""
    c = [("console 8 poutres", _console(), None)]
    # 120 barres est RETIRÉ du corpus par défaut : il coûte 110 s pour une
    # information que 54 et 60 donnent déjà (la frontière est erratique).
    # `domaines.frontiere_fermee` le garde pour qui veut la carte complète.
    for n in (20, 40, 52, 54, 60, 80):
        c.append((f"chaîne {n} BARRES (précontrainte seule)",
                  _cas_chaine(n, False), None))
    for n in (20, 40, 60, 120):
        c.append((f"chaîne {n} poutres (raideur matérielle)",
                  _cas_chaine(n, True), None))
    c.append(("hexapode (6 boucles, 0 ddl)", _hexapode(), None))
    return c


def campagne(bavard=True):
    res = {}
    for nom, f, _ in cas():
        t0 = time.time()
        try:
            r = f()
            res[nom] = dict(ok=True, valeur=r["valeur"], s=time.time() - t0)
        except Exception as ex:                       # noqa: BLE001
            res[nom] = dict(ok=False, motif=str(ex)[:70], s=time.time() - t0)
    if bavard:
        print("╔═ vinkulum — CAMPAGNE de non-régression du solveur statique")
        for nom, r in res.items():
            v = f"{r['valeur']:+.6f}" if r["ok"] else "—"
            print(f"║ {'OK   ' if r['ok'] else 'REFUS'} {nom:44s} {v:>12s} "
                  f"{r['s']:5.1f} s")
        n_ok = sum(1 for r in res.values() if r["ok"])
        print(f"║ {n_ok}/{len(res)} convergent")

    # ── CE QUE LE CORPUS NE PROUVE PAS ENCORE, et pourquoi ──
    #
    # Le contrôle qu'il faudrait est un JUMEAU à raideur matérielle du câble :
    # même poids, même portée, même forme d'équilibre, la raideur venant de la
    # matière au lieu de la précontrainte. Il départagerait « le défaut vient
    # de la précontrainte seule » d'une simple affaire de taille.
    #
    # TROIS FORMES ESSAYÉES LE 3 SEPT., AUCUNE NE RÉSOUT LE MÊME PROBLÈME :
    #
    #   poutres créées SUR LE V   flèche 3,0037 — elles y sont NON DÉFORMÉES,
    #                             donc le V est un équilibre : ça « converge »
    #                             sans avoir rien résolu
    #   poutres créées DROITES    flèche 0,0037 — l'ancrage COMPRIME au lieu
    #                             de faire fléchir ; une ligne droite sous
    #                             compression est un équilibre (instable)
    #   + défaut initial en sinus flèche 0,5037 — elle reste EXACTEMENT sur le
    #                             défaut imposé (0,05 × 10 m), à tout n
    #
    # ⇒ fabriquer le jumeau est un problème de CONTINUATION (charge ou
    # raccourcissement par paliers, avec suivi de branche), pas un choix de
    # raideurs. C'est un chantier à part, et il est nommé ici plutôt
    # qu'esquivé. Les lignes « poutres » du corpus restent : elles mesurent ce
    # qu'elles mesurent — une chaîne raide qui converge — et pas davantage.
    fl = [r["valeur"] for k, r in res.items()
          if "poutres" in k and "chaîne" in k and r["ok"]]
    if bavard:
        print("║ ce que le corpus NE prouve pas : le jumeau à raideur matérielle")
        print("║ d'un câble reste à construire — trois formes essayées, aucune ne")
        print(f"║ résout le même problème (flèches {min(fl):.4f} contre 2,6546 aux barres)")
    # garde du CONSTAT : tant que les poutres ne descendent pas près de la
    # caténaire, la ligne ci-dessus est vraie et doit rester écrite.
    assert fl and min(fl) > 0.4, ("le jumeau à poutres atteint enfin la "
                                  "caténaire : le contrôle falsifiant devient "
                                  "possible, rejouer ce bloc", fl)
    # ── LA FRONTIÈRE EST FERMÉE, ET C'EST CE BANC QUI L'A DIT ──
    # Il a rougi le 3 sept. sur son propre garde-fou de constat (« la
    # frontière des barres est toujours là ») quand elle a disparu. Le
    # garde-fou change donc de sens : toutes les chaînes de barres doivent
    # converger, et si l'une repassait au rouge ce serait une RÉGRESSION.
    bar = [(k, r["ok"]) for k, r in res.items() if "BARRES" in k]
    assert all(o for _, o in bar), ("RÉGRESSION : la frontière erratique de la "
                                    "statique sur les câbles est revenue", bar)
    return res


def demo():
    r = campagne()
    print("╚ campagne OK — le corpus qui permet de JUGER un changement du solveur")
    return r


if __name__ == "__main__":
    r = campagne(bavard="--json" not in sys.argv)
    if "--json" in sys.argv:
        print(json.dumps(r, indent=1, sort_keys=True))
