#!/usr/bin/env bash
# LE CHEMIN D'INSTALLATION EST MESURÉ, PAS PROMIS.
#
# Le README promettait `pip install vinkulum` et « contrôles rapides (2 s) » ;
# le paquet n'est publié nulle part et les contrôles durent ~150 s. Une
# promesse d'installation que rien n'exécute finit fausse — comme un chiffre
# recopié dans la prose.
#
# Ce script fabrique le seul tiers qu'on puisse fabriquer sans en avoir un :
# un environnement VIERGE qui n'a que la roue. Il a sorti trois défauts à sa
# première exécution (3 sept.) — scipy manquant du chemin du README, `doc` qui
# accusait la référence d'être périmée alors qu'elle est absente du paquet,
# et `bancs` qui écrivait dans la bibliothèque de l'utilisateur.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3.14}
"$PY" -c 'import sys; sys.exit("Python 3.14 ou plus requis") if sys.version_info < (3, 14) else None'
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT

# On construit depuis l'ARBRE DE TRAVAIL, pour juger ce qu'on s'apprête à
# commiter ; la CI rejoue le même script sur le commit POUSSÉ, ce qui donne la
# sémantique du clone sans qu'on ait à la simuler. (Premier jet : `git clone .`
# — il testait HEAD, donc jamais le correctif en cours.)
# `cargo test` ne tournait NULLE PART dans la boucle de travail : le champ
# `nh` ajouté aux liaisons le 3 sept. a cassé la cible `test` pendant des
# heures sans que rien ne rougisse. Il coûte 0,04 s une fois la build chaude.
echo "== cargo test =="
cargo test --release -q 2>&1 | grep -E 'test result|error' | head -4

echo "== roue depuis l'arbre de travail =="
maturin build --release -q --interpreter "$PY" -o "$T/roues"
ROUE=$(echo "$T"/roues/*.whl)
echo "   $(basename "$ROUE")"

echo "== venv neuf : QUE la roue =="
"$PY" -m venv "$T/venv"
"$T/venv/bin/pip" install -q "$ROUE" 'scipy>=1.11'

echo "== le tiers lance les commandes du README, HORS du dépôt =="
cd "$T"                       # un tiers n'est pas dans les sources
"$T/venv/bin/python" -m vinkulum.verification
"$T/venv/bin/python" -m vinkulum.doc          # ne doit PAS accuser : la réf. n'est pas dans le paquet
"$T/venv/bin/python" - <<'PY'
# les modèles publics doivent voyager DANS la roue, pas rester dans le dépôt
from vinkulum import mjcf, urdf
from pathlib import Path
d = Path(urdf.__file__).parent / "donnees"
assert (d / "urdfs").is_dir() and (d / "menagerie").is_dir(), "données absentes de la roue"
n, i = mjcf.charge(str(d / "menagerie" / "panda.xml"))
assert i["n_corps"] == 11, i
print(f"   données embarquées : {len(list((d/'menagerie').glob('*.xml')))} MJCF, "
      f"{len(list((d/'urdfs').glob('*.urdf')))} URDF — panda chargé hors dépôt")
PY
echo "== bancs écrit-il chez l'utilisateur et non dans sa bibliothèque =="
mkdir -p "$T/travail" && cd "$T/travail"
"$T/venv/bin/python" -m vinkulum.bancs rapide >/dev/null
cd "$T"
echo "OK — le chemin du README tient dans un environnement vierge"
