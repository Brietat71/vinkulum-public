#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
: "${PY:?Définir PY vers le Python contenant le noyau Vinkulum à confronter}"
command -v lake >/dev/null || { echo 'Lean 4.19.0 requis : voir preuves/README.md' >&2; exit 1; }
case "$(lean --version)" in
  'Lean (version 4.19.0,'*) ;;
  *) echo 'Version Lean différente de 4.19.0' >&2; exit 1 ;;
esac
lake build
lake build garde
AUDIT_VINKULUM=$(mktemp)
trap 'rm -f "$AUDIT_VINKULUM"' EXIT
lake env lean Audit.lean > "$AUDIT_VINKULUM"
"$PY" - "$AUDIT_VINKULUM" <<'PY'
import re
import sys
from pathlib import Path
s = Path(sys.argv[1]).read_text()
attendus = {'garde_exact', 'decode_exact', 'decode_bornes', 'bit_implicite',
            'boucle_exacte', 'decalage_exact', 'produit_u128', 'exposants_rust',
            'tolerance_admise', 'fraction_masque', 'champ_masque', 'signe_decalage',
            'refus_non_fini'}
lignes = re.findall(r"'Vinkulum\.(\w+)' depends on axioms: \[([^\]]*)\]", s)
if {n for n, _ in lignes} != attendus or len(lignes) != len(attendus):
    raise SystemExit('audit incomplet ou format inattendu : ' + s)
for nom, axiomes in lignes:
    if set(axiomes.split(', ')) - {'propext', 'Classical.choice', 'Quot.sound'}:
        raise SystemExit(f'{nom}: axiome non autorisé : {axiomes}')
print(s, end='')
PY
"$PY" contre_epreuves.py
