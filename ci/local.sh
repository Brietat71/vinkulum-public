#!/usr/bin/env bash
# CI LOCALE — également appelée par .github/workflows/ci.yml. Branchée en `pre-push` par
# `git config core.hooksPath ci/hooks` (fait par ci/installe_hook.sh).
#
#   ci/local.sh            fmt · clippy -D · cargo test · roue → venv · verification · doc
#   ci/local.sh --bancs    … puis bancs rapide et contact
# Les campagnes indépendantes utilisent jusqu'à 8 processus, dans le budget
# CPU disponible. VINKULUM_BANCS_JOBS=1 impose le séquentiel ;
# VINKULUM_BANCS_CPUS borne le budget total, VINKULUM_BANCS_LOGS les journaux.
#
# PY = l'interpréteur qui reçoit la roue (défaut : venv actif, sinon .venv du dépôt).
# Un échec rend non zéro : le push est refusé. `git push --no-verify` passe outre,
# et c'est écrit ici pour que ce soit un choix, pas un oubli.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH"
PY=${PY:-${VIRTUAL_ENV:-$PWD/.venv}/bin/python}
if [[ ! -x "$PY" ]]; then
  printf 'Python introuvable : %s. Créez .venv ou définissez PY.\n' "$PY" >&2
  exit 1
fi
VINKULUM_VENV=$("$PY" -c '
import sys
if sys.version_info < (3, 14):
    raise SystemExit("Python 3.14 ou plus requis ; recréez le venv ou définissez PY")
if sys.prefix == sys.base_prefix:
    raise SystemExit("PY doit désigner un venv")
print(sys.prefix)
')
MATURIN=${MATURIN:-$VINKULUM_VENV/bin/maturin}
if [[ ! -x "$MATURIN" ]]; then
  MATURIN=$(command -v "$MATURIN" || command -v maturin) || {
    printf 'maturin introuvable : installez-le ou définissez MATURIN.\n' >&2
    exit 1
  }
fi
INSTALL=()
if command -v uv >/dev/null 2>&1; then
  INSTALL+=(--uv)
elif ! "$PY" -c 'import pip' >/dev/null 2>&1; then
  printf 'Installation : uv ou pip dans le venv est requis.\n' >&2
  exit 1
fi
etape() { printf '\n══ %s\n' "$*"; }

etape "format et lint"
cargo fmt --check
cargo clippy --release --all-targets -- -D warnings

etape "tests du noyau"
cargo test --release -q

etape "contrôles du prototype de poutre mixte"
cargo test --release --example poutre_mixte -q

etape "roue → $PY"
VIRTUAL_ENV="$VINKULUM_VENV" "$MATURIN" develop "${INSTALL[@]}" --release --extras verification

etape "contrôles rapides (campagne parallèle)"
"$PY" -m vinkulum.verification

etape "charges temporelles et liaisons à deux repères : références indépendantes"
"$PY" -m unittest vinkulum.test_charges_temporelles

etape "régressions Python : noyau, analyses, campagnes et cinématique sigma"
"$PY" -m unittest vinkulum.test_modes_creux vinkulum.test_certification_creuse vinkulum.test_certification vinkulum.test_certification_assemblage vinkulum.test_certification_quotient vinkulum.test_noyau vinkulum.test_assemblage vinkulum.test_analyses vinkulum.test_campagnes vinkulum.test_schema_lie vinkulum.test_operateurs vinkulum.test_reduction_ports vinkulum.test_reponses_groupees vinkulum.test_reduction_contrainte vinkulum.test_fiabilite_lisse vinkulum.test_validation

etape "épreuves physiques : covariance des projections et convergence du pendule"
VINKULUM_AUDIT_ASSEMBLAGE=$(mktemp)
"$PY" ci/audit_assemblage_physique.py --sortie "$VINKULUM_AUDIT_ASSEMBLAGE"
rm "$VINKULUM_AUDIT_ASSEMBLAGE"

etape "diagnostic énergie–moment : trace non mutante et attribution exacte"
"$PY" ci/test_trace_moments_em.py
"$PY" ci/test_attribution_moments_em.py
"$PY" ci/test_attribution_rotation_em.py
"$PY" ci/test_archive_attribution_reperes.py

etape "prototype de certificat géométrique : intervalles, export, oracle et documents"
"$PY" ci/test_prototype_assemblage_certifie.py
"$PY" ci/test_export_geometrie_certifiee.py
"$PY" ci/test_oracle_polynomial_assemblage.py
"$PY" ci/test_archive_assemblage_certifie.py
"$PY" ci/test_archive_assemblage_public.py

etape "garantie temporelle du pendule : intervalles, référence et archives"
"$PY" ci/test_intervalle_temporel.py
"$PY" ci/test_taylor_temporel.py
"$PY" ci/test_reference_temporelle.py
"$PY" ci/test_archive_trajectoire_certifiee.py
"$PY" ci/test_trace_pendule_courante.py

etape "certificats creux : facteurs, perturbations, archives et matrices natives"
"$PY" ci/test_lineaire_creux.py
"$PY" ci/test_dependances_structurelles.py
"$PY" ci/test_archive_lineaire_creux.py
"$PY" ci/test_lineaire_creux_natif.py
"$PY" ci/test_archive_certification_creuse.py
"$PY" ci/test_archive_modes_creux.py
"$PY" ci/test_confronte_modes.py
"$PY" ci/test_archive_confrontation_modale.py

etape "certification : garde exact des vitesses et obligations ouvertes"
VINKULUM_DOSSIER_CERTIFICATION=$(mktemp)
"$PY" ci/certifie_noyau.py --sortie "$VINKULUM_DOSSIER_CERTIFICATION"
rm "$VINKULUM_DOSSIER_CERTIFICATION"

etape "preuves Lean et correspondance du garde exact"
PY="$PY" preuves/verifier.sh

etape "la référence d'API est-elle à jour"
"$PY" -m vinkulum.doc

etape "archive du prototype de poutre mixte"
"$PY" ci/archive_poutre_mixte.py --verifier docs/bancs/poutre-mixte-2026

etape "quotient orthogonal local : contre-épreuves et archive"
"$PY" ci/test_contraintes_orthogonales.py
"$PY" ci/mesure_quotient_orthogonal.py --verifier docs/bancs/quotient-orthogonal-2026

etape "réduction par ports : identités, enveloppes et confrontation HCB"
"$PY" ci/test_ports_krylov.py
"$PY" ci/test_reference_ports_hcb.py
"$PY" ci/mesure_ports_krylov.py --verifier docs/bancs/ports-krylov-2026

etape "ports relevés : énergie, précision indépendante et assemblage"
"$PY" ci/test_energie_ports.py
"$PY" ci/test_reference_ports_precision.py
"$PY" ci/test_ports_releves.py
"$PY" ci/mesure_ports_releves.py --verifier docs/bancs/ports-releves-2026

etape "minoration spectrale : inverse sélectionnée et arithmétique dirigée"
"$PY" ci/test_inverse_selectionnee.py
"$PY" ci/test_trace_complement.py

etape "complément contraint : certificat dirigé, composition et inertie"
"$PY" ci/test_trace_complement_dirigee.py
"$PY" ci/test_retention_complement.py
"$PY" ci/test_inertie_complement.py
"$PY" ci/test_bornes_produits_creux.py
"$PY" ci/test_archive_complement_spectral.py

etape "Krylov contraint : projecteurs, enveloppes et réponses physiques"
"$PY" ci/test_inverse_contrainte_energie.py
"$PY" ci/test_racine_masse_blocs.py
"$PY" ci/test_krylov_contraint_preuves.py
"$PY" ci/test_krylov_contraint.py
"$PY" ci/test_archive_krylov_contraint.py

etape "inertie dirigée : congruences exactes et certificat creux"
"$PY" ci/test_inertie_dirigee_preuves.py
"$PY" ci/test_inertie_complement_dirigee.py
"$PY" ci/test_archive_inertie_contrainte.py

etape "contrôle par facteurs : majorations et champs physiques"
"$PY" ci/test_controle_facteurs_preuves.py
"$PY" ci/test_controle_facteurs.py
"$PY" ci/test_archive_controle_facteurs.py

etape "inertie binary64 : intervalles compilés, séparateurs et archives"
"$PY" ci/test_inertie_binaire_compile.py
"$PY" ci/test_archive_inertie_binaire.py

etape "masses couplées : comparaison, quotient, oracle et archives"
"$PY" ci/test_comparaison_masse.py
"$PY" ci/test_controle_masse_comparee.py
"$PY" ci/test_oracle_champs_couples.py
"$PY" ci/test_juge_masse_couplee.py
"$PY" ci/test_modeles_masse_consistante.py
"$PY" ci/test_archive_masse_comparee.py

etape "assemblages contraints : travail, modes privés, marges globales et archives"
"$PY" ci/test_assemblage_complements_preuves.py
"$PY" ci/test_assemblage_complements.py
"$PY" ci/test_oracle_assemblage_complements.py
"$PY" ci/test_archive_assemblage_complements.py
"$PY" ci/test_certificat_paquet.py
"$PY" ci/test_livraison_contrainte.py

etape "champ intérieur : observables, enrichissement et refus près des résonances"
"$PY" ci/test_champ_interieur.py
"$PY" ci/experience_champ_interieur.py --verifier docs/bancs/champ-interieur-2026

etape "confrontation des ports : oracle, juge, pont HCB et archive"
"$PY" ci/test_oracle_champs_ports.py
"$PY" ci/test_juge_ports_exudyn.py
# NumPy/SciPy suffisent aux régressions algébriques du pont ; seuls les
# tests d'intégration officiels sont ignorés si Exudyn n'est pas installé.
"$PY" ci/test_reference_hcb_exudyn.py
# Vérification des mesures conservées ; aucune campagne n'est relancée.
"$PY" ci/test_archive_ports_exudyn.py
"$PY" ci/test_reponses_groupees_archive.py

etape "archive de l'expérience sigma"
"$PY" ci/bilan_sigma.py docs/bancs/integration-sigma-2026.json

etape "archive de globalisation statique"
"$PY" ci/bilan_globalisation.py docs/bancs/globalisation-statique-0.8.1.json

etape "archive de la tangente analytique des poutres"
"$PY" ci/bilan_globalisation.py docs/bancs/tangente-poutre-0.8.2.json

etape "confrontations externes et contrôles du classement"
"$PY" ci/archive_confrontation.py --verifier docs/bancs/confrontation-statique-mbdyn-0.8.2
"$PY" ci/confronte_exudyn.py --verifier docs/bancs/confrontation-exudyn-0.8.2
"$PY" ci/test_confrontation_exudyn.py

etape "adjoints par produits : archives et contrôles de comparaison"
"$PY" ci/mesure_adjoint_operateurs.py --verifier docs/bancs/adjoint-operateurs-0.9.0-essais.json.gz
"$PY" ci/mesure_adjoint_operateurs.py --verifier docs/bancs/adjoint-operateurs-0.9.0-livraison-essais.json.gz
"$PY" ci/test_adjoint_operateurs.py

etape "dossier de fiabilité : inventaire, décisions et contre-épreuves"
"$PY" ci/test_qualifie_validation.py
"$PY" ci/verifie_fiabilite.py docs/bancs/version-0.14.1.json
"$PY" ci/test_verifie_fiabilite.py

if [[ "${1:-}" == "--bancs" ]]; then
  etape "bancs rapide"
  "$PY" -m vinkulum.bancs rapide
  etape "bancs de contact"
  "$PY" -m vinkulum.contact
fi
printf '\n╚ CI locale OK\n'
