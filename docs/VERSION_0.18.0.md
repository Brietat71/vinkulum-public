# Vinkulum 0.18.0 — attribution vérifiée des écarts de repères

Livraison identifiée par le tag annoté `v0.18.0`. Le [dossier de qualification](bancs/version-0.18.0/qualification.json) conserve les empreintes et les journaux.

Cette version ajoute des outils d'attribution a posteriori pour la toupie
libre de l'intégrateur énergie–moment. Ils expliquent et encadrent les
contributions à la différence entre deux trajectoires discrètes stockées.
La nouvelle fonctionnalité de diagnostic motive une version mineure ;
les équations et opérations numériques du pas restent celles de la 0.17.0.

## Changement concret

Une entrée privée, `Noyau._trace_em`, expose sur une copie du modèle les
moments réellement portés par le pas, les rotations, les vitesses spatiales
et les inerties stockées. Les outils dans `ci/` relisent ces valeurs comme
des rationnels exacts. L'état source n'est pas modifié, même si le calcul
échoue après avoir avancé un premier corps de la copie.

Les écarts sont séparés en contributions initiales, contributions des
données transformées, défauts algébriques propagés et post-traitement des
observations. Les identités quadratiques et de résolvante sont vérifiées
exactement. Les enclosures sont dirigées et doivent recomposer les
observations à chaque pas. Les transformations stockées ne sont jamais
supposées exactement orthogonales.

Voir la [dérivation et le protocole](ATTRIBUTION_REPERES_EM.md). Les outils
refusent les horizons incomplets, les grilles incohérentes et les budgets
excessifs. Une limite de nombre de pas reste explicitement un calcul partiel.

## Résultats et traçabilité

La capture compare les dates, rotations et vitesses aux sorties ordinaires
de la roue candidate et de la roue publiée 0.17.0 : **560 012 points,
aucune composante binary64 différente** sur les douze trajectoires.
Les douze couples de maxima historiques de la 0.14.1 sont aussi retrouvés
exactement, y compris le raffinement où l'écart augmente.

Les neuf attributions couvrent trois changements de repères et trois pas,
sur 20 s : **420 000 pas comparés, tous admis**. L’archive complète est
vérifiée, ainsi que ses 36 pics. À chacun des neuf pics de vitesse, les
bornes établissent une contribution dominante des défauts de pas propagés,
sous la décomposition sécante déclarée.
Les écarts numériques observés sont conservés ; aucun seuil physique
historique n'est modifié pour les reclasser.

Le dossier de capture identifie la roue candidate et son extension native.
L’[archive](bancs/attribution-reperes-0.18.0/manifest.json) conserve les douze traces complètes, les neuf résultats,
les sources exactes des vérificateurs et les journaux de calcul. Le contrôle
d'intégrité ne se substitue pas au rejeu rationnel complet.

## Validation de livraison

Les 209 tests du paquet isolé et les 20 tests de trace, d’attribution et
d’archive passent. Ces derniers incluent les altérations de résultats, les
refus de réflexions et le rattachement historique. La CI locale étendue passe :
87 tests Rust, huit tests du prototype de poutre, 41 vérifications rapides,
46 bancs et neuf bancs de contact. Les 13 théorèmes Lean et 13 314
contre-épreuves exactes sont contrôlés. Trois intégrations du pont Exudyn
sont ignorées dans cet environnement sans Exudyn ; les archives de la
confrontation indépendante restent vérifiées.

La roue finale est installée dans un environnement neuf ; ses 209 tests
Python et les cinq tests de trace passent aussi hors du dépôt. Tous ses
modules Python et son extension sont identiques au candidat mesuré. Seuls
le README dans METADATA, la date et l’UUID du SBOM, puis RECORD diffèrent.
La comparaison détaillée est conservée dans le dossier de qualification.

Une première CI a été invalidée par la modification de son script pendant
que Bash le lisait : erreur de shell après des suites déjà terminées.
Ce journal est conservé ; il ne vaut pas une CI de livraison réussie.

La CI étendue réussie a été suivie d’une CI complète après correction du
lancement du rejeu : `-B` empêche la création de bytecode dans l’archive.
Une contre-épreuve exécute réellement un pas sur une copie, vérifie son
incomplétude déclarée et toutes ses empreintes avant/après. Les deux
journaux de CI restent distincts dans le dossier de qualification.

## Portée maintenue

Cette attribution concerne les différences entre les trajectoires discrètes
contrôlées. Elle n'est pas une borne d'erreur à la solution continue, une
preuve de toute l'implémentation ni une certification des contacts ou de
tous les mécanismes. Elle n'établit pas non plus un classement général
contre MBDyn ou Exudyn. La [confrontation modale](CONFRONTATION_MODALE.md)
conserve son domaine, ses résultats défavorables et ses budgets.

La formalisation Lean déjà livrée est désormais correctement mentionnée
dans le dossier général de certification. Les obligations géométriques
et de rang y distinguent les domaines 0.15.0/0.16.0 réalisés de leurs
extensions encore ouvertes. L'[audit du plan](AUDIT_PLAN_FIABILITE.md)
réunit les preuves de chaque livrable.
