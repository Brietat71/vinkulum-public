# Audit de réalisation du plan de fiabilité

Cet audit suit les six livrables de [PLAN_FIABILITE.md](PLAN_FIABILITE.md).
Il ne redéfinit pas leur réalisation comme une certification du noyau entier.
État du lot final : neuf attributions admises et archivées, roue finale
qualifiée et CI étendue réussie. La livraison est identifiée par `v0.18.0`.

| Exigence | Éléments contrôlables | Portée et état |
|---|---|---|
| 1. Référence 0.14.0 identifiée et reconstruite | `docs/bancs/version-0.14.1.json`, archive de preuves associée, `ci/verifie_fiabilite.py` | Identité des paquets, deux campagnes de 62 métriques, 44 admises et 18 écarts conservés ; contrôles d'archive rejoués lors du lot final |
| 1. Suivi des écarts historiques | Registre dans l'archive `fiabilite-0.14.1-preuves.json.gz` ; `VERSION_0.14.1.md` | Les écarts physiques restent nommés ; aucune clôture générale fondée sur la seule CI |
| 2. Newton énergie–moment sur petites inerties | Correction 0.14.1, corpus avant/après, `python/vinkulum/test_fiabilite_lisse.py` | 32 violations historiques éliminées sur le corpus ; échelles d'inertie de 1e-20 à 1e20 et convergence indépendante |
| 2. Refus et restitution | `python/vinkulum/test_noyau.py`, `test_validation.py`, `ci/test_trace_moments_em.py` | Contrôles des chemins testés, y compris échec après modification de la copie de diagnostic ; aucune preuve de toutes les sorties Rust |
| 2. Diagnostic des repères à 20 s | `ATTRIBUTION_REPERES_EM.md`, traces natives et vérificateur rationnel | Neuf attributions admises, 420 000 pas comparés ; archive complète et 36 pics contrôlés ; qualification 0.18.0 conservée |
| 3. Décodage binary64 et garde entier | `preuves/README.md`, sources Lean, `preuves/verifier.sh` | 13 théorèmes et audit des axiomes ; 13 314 confrontations Lean/Rust/Fraction ; pas de preuve du compilateur ou de tout le programme Rust |
| 4. Assemblage géométrique certifié | `CERTIFICATION_ASSEMBLAGE.md`, `ci/test_archive_assemblage_public.py`, API `certifier_assemblage` | Existence et unicité locales par Krawczyk, jauges explicites, enclosures géométriques, contraintes originales contrôlées ; petits mécanismes admis |
| 5. Première garantie temporelle | `GARANTIE_TEMPORELLE_PENDULE.md`, traces et certificats archivés, `ci/test_archive_trajectoire_certifiee.py` | Pendule sur 1 s, trois grilles, Picard et Taylor avec reste ; bornes aux instants enregistrés, sans interpolation certifiée ni généralisation aux contacts |
| 6. Opérateurs creux et dépendances perturbées | `CERTIFICATION_CREUSE.md`, API 0.16.0, tests d'archives creuses et de dépendances structurelles | Certificats sur les opérateurs et perturbations structurées admis, KKT qualifiés jusqu'à 2 816 inconnues ; rang arbitrairement perturbé hors domaine |
| 6. Qualification modale | `MODES_CREUX.md`, `ci/test_archive_modes_creux.py`, API 0.17.0 | 18 cas, trois familles, jusqu'à 1 024 éléments ; résidus physiques et références indépendantes ; le rang modal numérique n'est pas déclaré certifié |
| 6. Confrontation MBDyn et Exudyn | `CONFRONTATION_MODALE.md`, protocole et archive de 140 essais, `ci/test_archive_confrontation_modale.py` | Précision commune avant classement, temps/mémoire limités, configurations en arbre, échecs et résultats défavorables conservés ; aucun classement général des solveurs |
| Livraison finale | Roue isolée, notes de version, CI locale étendue, commit, tag et remote | Roue finale qualifiée, modules identiques au candidat, 209 tests isolés, 20 tests de diagnostic, 87 tests Rust, 46 bancs et 9 contacts ; tag `v0.18.0`, journaux dans `docs/bancs/version-0.18.0` |

## Ce que les contrôles établissent

Une vérification d'archive contrôle l'intégrité, les décisions et les
contre-épreuves prévues par son vérificateur. Elle ne recrée pas une mesure
de performance historique ni une nouvelle campagne expérimentale. Les
rejeux mathématiques exécutent les calculs exacts explicitement décrits.
La CI courante ajoute les régressions du paquet effectivement reconstruit.
Ces niveaux d'évidence ne sont pas interchangeables.

Le lot énergie–moment doit préserver les trajectoires ordinaires mesurées :
le dossier de capture compare tous les bits des dates, rotations et vitesses
à ceux de la roue 0.17.0. Les moments portés sont lus dans l'intégrateur,
et non reconstruits à partir des sorties. Chaque attribution doit refuser
un horizon incomplet ou un budget excessif. Le diagnostic peut être fermé
comme origine quantifiée des écarts du corpus ; les écarts numériques eux-mêmes
ne deviennent pas nuls et aucune précision universelle n'en découle.

## Obligations générales maintenues

La preuve du programme compilé, les trajectoires des mécanismes généraux,
les événements et contacts non lisses, les erreurs de modèle et la validation
expérimentale des domaines d'usage restent des travaux distincts. Le
[registre de certification](CERTIFICATION_NOYAU.md) les conserve. Terminer
ce plan de premiers livrables ne permet pas de les déclarer réalisés.
