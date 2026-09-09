# Première publication publique — 0.18.1

La publication est un instantané des sources livrées en 0.18.0, complété par
Apache-2.0, les notices tierces, les métadonnées 0.18.1 et les instructions
publiques. Le code numérique reste inchangé. Le dépôt de travail historique
reste privé : des messages de commit et une ancienne pull request contiennent
des liens de sessions de travail qui n'ont pas vocation à être publiés.

## Périmètre contrôlé

L'audit préalable a porté sur les branches distantes, les tags et la référence
de pull request récupérables, ainsi que sur les fichiers présents dans leurs
historiques. Gitleaks 8.30.1 a examiné les changements Git et une extraction
des 3 013 blobs uniques : textes, 458 contenus gzip, 50 archives ZIP/NPZ,
quatre PDF et chaînes lisibles de 13 fichiers binaires. Les tableaux NPY
numériques ont été traités comme des données, avec inspection des en-têtes.

Les deux scans ont chacun produit 25 alertes `generic-api-key`, toutes
identifiées comme des empreintes SHA-256 de fichiers dans les manifestes.
Aucun secret d'accès n'a été identifié par ces contrôles. Ce constat n'est
pas une preuve formelle de l'absence de toute information confidentielle.

Le dépôt public ne reprend ni les messages de commits privés, ni les
discussions, ni les configurations locales, ni les fichiers non suivis du
répertoire de travail. Il conserve les dossiers scientifiques suivis, dont
les références techniques à FRELON, premier cas d'application. Aucun dossier
CAO ni fichier de configuration privé de cette application n'est ajouté.

Les 36 fichiers XML/URDF tiers ont une provenance publique vérifiée et leurs
notices sont distribuées. Le graphe Cargo verrouillé fait l'objet d'un
registre des licences. Voir [les attributions](../THIRD_PARTY_NOTICES.md).

## Reproduction et limites

Suivre les instructions du README depuis un clone neuf. Le projet requiert
Python 3.14+, Rust/Cargo et un compilateur C++17. Les campagnes contre des
solveurs externes nécessitent leur installation ; elles ne font pas partie
de l'exemple minimal.

Les dossiers historiques conservent leurs chemins locaux et leurs identités
de commits à titre de provenance. Certaines archives externes volumineuses
ne sont pas incluses ; leurs chemins ne doivent pas être interprétés comme
des ressources téléchargeables. Les certificats et sources figées présents
dans le dépôt restent consultables avec leurs vérificateurs documentés.

La publication ne change aucune garantie scientifique. Le noyau demeure
alpha, les 18 écarts physiques historiques restent documentés, et les
comparaisons MBDyn/Exudyn conservent leur périmètre d'origine. Le
[registre de certification](CERTIFICATION_NOYAU.md) fait foi sur les garanties.

## Validation de la livraison

La construction de la roue a été effectuée depuis un répertoire de compilation
neuf. Son installation dans un environnement indépendant du dépôt a exécuté
l'exemple du README et les 41 vérifications mécaniques avec succès. Les
36 modèles et les notices des 150 dépendances ont été contrôlés ; les notices
de la roue sont identiques octet par octet à celles du dépôt.

La CI locale complète est passée : 87 tests Rust, huit contrôles de poutre
mixte, 209 tests Python, les contre-épreuves et vérificateurs d'archives,
13 théorèmes Lean audités et 13 314 confrontations sans divergence. Trois
tests d'intégration Exudyn sont ignorés faute d'installation de ce solveur ;
les archives des comparaisons sont vérifiées séparément par la même CI.
Les bancs étendus de la 0.18.0 ne sont pas réattribués à cette livraison.

Le premier essai de CI a été arrêté par l'absence de Lean dans le `PATH`.
Le passage complet a été relancé avec Lean 4.19.0 et le cache Mathlib figé.
Il ne s'agit pas d'une reconstruction de Mathlib depuis ses sources.

Les journaux et l'empreinte de la roue sont conservés dans
[`bancs/publication-0.18.1/`](bancs/publication-0.18.1/).
