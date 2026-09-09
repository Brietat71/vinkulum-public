# Priorités de fiabilité du noyau

Choix retenu : fiabilité démontrable, en commençant par les mécanismes
rigides à contraintes holonomes sans contact. Les certificats algébriques
existants restent locaux. Aucun résultat ci-dessous ne certifie le noyau entier.

| Étape | Livrable et critère de sortie | État |
|---|---|---|
| 1. Référence reproductible | Roue 0.14.0 identifiée, reconstruction comparée, CI et campagne physique rejouées ; suivi individuel des lignes historiques | Dossier de la 0.14.1 |
| 2. Défauts du noyau lisse | Contre-épreuves avant/après, conservation des seuils physiques, état restitué sur les refus testés | Correction du Newton énergie–moment et de la gestion des erreurs de validation ; diagnostic des repères quantifié par neuf attributions exactes encadrées sur 420 000 pas comparés ; voir `ATTRIBUTION_REPERES_EM.md` |
| 3. Consolider les preuves | Formaliser le décodage binary64 et le garde entier dans Lean 4, documenter la correspondance Rust et la base de confiance | 13 théorèmes audités sous Lean 4 ; 13 314 contre-épreuves Lean/Rust/Fraction sans divergence ; contrôle obligatoire en CI et base de confiance dans `preuves/README.md` |
| 4. Assemblage géométrique | Intervalles de contraintes/Jacobien et inclusion de Krawczyk, petits mécanismes carrés avec jauge explicite, contrôle de toutes les contraintes originales | API 0.15.0 en lecture seule, schéma public compatible, bornes de distance recalculées et cinq documents produits par la roue isolée ; domaine et preuves dans `CERTIFICATION_ASSEMBLAGE.md` |
| 5. Première garantie temporelle | Pendule à pas fixe, horizon fini, référence indépendante validée par intervalles ; bornes aux instants contrôlés | Premier banc réalisé sur la roue 0.15.0 : trois traces, horizon 1 s, intervalles rationnels, Picard et Taylor avec reste ; certificats autonomes et contre-épreuves dans `GARANTIE_TEMPORELLE_PENDULE.md` ; pas d'API générale |
| 6. Portée et comparaison | Opérateurs creux, dépendances perturbées, qualification du chantier modal puis confrontation MBDyn/Exudyn à précision commune | Réalisé dans les domaines déclarés : API 0.16.0, KKT qualifiés jusqu'à 2 816 inconnues et quotient structurel autonome ; voir `CERTIFICATION_CREUSE.md`. Analyse modale 0.17.0 qualifiée sur 18 cas jusqu’à 1 024 éléments (`MODES_CREUX.md`) ; confrontation modale à précision commune réalisée : 140 essais, configurations Exudyn en arbre incluses, résultats défavorables et budgets conservés (`CONFRONTATION_MODALE.md`) |

## Premier lot : 0.14.1

Le [dossier de livraison](VERSION_0.14.1.md) décrit la correction, les
contre-épreuves et les résultats. Le registre distingue défaut logiciel,
limite de modèle, conditionnement, limitation du contrat et résultat
inexpliqué. Une anomalie ne devient pas « close » parce qu'un test voisin passe.

Le travail modal non commité de `/root/vinkulum` est préservé. Les travaux
présents partent du commit livré `e90f9ac`, dans un clone distinct.

## Règles des prochains lots

- Ne pas modifier les seuils après lecture d'un résultat pour le faire passer.
- Associer chaque affirmation à ses hypothèses, une preuve ou mesure,
  un artefact vérifiable et un cas qui détecte une violation.
- Distinguer refus du solveur, refus de preuve et désaccord avec une mesure.
- Qualifier la roue isolée, exécuter la CI locale obligatoire avant push et
  conserver les contre-exemples avec leurs résultats défavorables.
- Publier une fonction de certification seulement après qualification du
  prototype ; préserver les schémas des certificats déjà livrés.

Une preuve de l'algorithme, une propriété de coefficients stockés, une
borne de trajectoire et une validation expérimentale sont quatre objets
différents. L'étape 3 est documentée dans `preuves/` ; l'étape 4 est intégrée
en 0.15.0. L'étape 5 qualifie un premier banc temporel sur cette roue, avec
des scripts et documents séparés. L'étape 6 dispose désormais d'une confrontation modale commune ; son domaine
est explicite et ne permet pas un classement général. Le diagnostic de
repères de l’étape 2 est désormais quantifié sur les douze trajectoires
archivées. Les écarts sont conservés, leur propagation est attribuée ;
aucune borne d’erreur à la solution continue n’en est déduite.
Les six livrables sont réalisés dans leurs domaines déclarés et réunis
par la livraison 0.18.0. La CI étendue et la roue isolée sont qualifiées ;
l’[audit des exigences](AUDIT_PLAN_FIABILITE.md) et le
[dossier de livraison](bancs/version-0.18.0/qualification.json) donnent les
évidences correspondantes. Les obligations générales de certification du
noyau demeurent explicitement ouvertes.
