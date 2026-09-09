# Vinkulum 0.15.0 — certificat géométrique local

La 0.15.0 ajoute `certifier_assemblage` à `vinkulum.certification`.
L'appel établit, dans le domaine admis, l'existence et l'unicité d'une pose
admissible dans une boîte donnée. Il borne sa distance à la pose native,
y compris l'écart de conversion de la rotation en quaternion.
L'état du noyau reste inchangé, y compris lorsque la certification refuse.

Le [contrat complet](CERTIFICATION_ASSEMBLAGE.md) décrit les coordonnées,
les jauges et les refus. Le domaine initial est limité à quatre corps,
liaisons holonomes à repères matériels stockés finis sans loi imposée et
distances positives. Toutes les contraintes originales restent présentes.
Cette capacité locale ne constitue pas une certification du noyau entier.

## Compatibilité

Nouvelle API compatible, donc augmentation mineure selon les règles du
projet. Les schémas linéaire et quotient conservent leurs vérificateurs.
`verifier_certificat` reconnaît désormais `vinkulum.assemblage.1` et refait
le calcul polynomial exact ainsi que les bornes de distance annoncées.
La lecture autonome fonctionne avec la bibliothèque standard seule.

Le calcul dynamique et les intégrateurs ne sont pas modifiés dans ce lot.
Le code Rust ajouté exporte une copie de la géométrie dans le domaine admis.
Les [preuves Lean du garde](../preuves/README.md), livrées après 0.14.1,
restent distinctes de ce certificat géométrique. Leurs artefacts compilés
Mathlib font explicitement partie de la base de confiance documentée.

## Qualification isolée

Roue : `vinkulum-0.15.0-cp314-cp314-linux_x86_64.whl`.
SHA-256 : `057f89c7256b5bfcb15abee35c831778c64edd3c167e6500ce01280d6fb3a89c`.
Installation dans un environnement séparé ; tests lancés avec `python -I`
depuis `/tmp`, sans import des modules Python du dépôt.

- 188 tests Python de régression réussis sur cette roue.
- Cinq certificats produits par l'API installée : trois échelles de
  longueur, un repère tourné et un double pendule. Les états natifs restent inchangés.
- Chaque certificat est relu par le vérificateur installé dans un processus
  Python `-S`, sans chargement de Vinkulum ni de ses dépendances numériques.
- Les empreintes des cinq modules concernés correspondent aux sources.
- 46/46 bancs rapides et 9/9 cas de contact réussis sur cette roue.
- 24 tests du prototype vérifient intervalles, dérivées, correspondance
  mécanique, cas singuliers, branches parasites et documents altérés.

Les [documents et leur qualification](bancs/assemblage-certifie-0.15.0/qualification.json)
conservent la roue identifiée, les empreintes des modules, les modèles,
les jauges, les poses initiales et les bornes. La relecture des cinq
documents fait partie de la CI obligatoire.

La campagne physique historique des 62 mesures reste celle de la 0.14.1 ;
elle n'est pas présentée comme rejouée dans ce lot. Le diagnostic de repères
à long terme, les garanties temporelles, les dépendances perturbées et les
comparaisons actualisées restent ouverts dans le [plan](PLAN_FIABILITE.md).
