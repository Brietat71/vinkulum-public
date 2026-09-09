# Version 0.8.0

La version 0.8.0 ajoute une correction géométrique expérimentale à
l'intégrateur generalized-alpha de Lie. Elle vise les rotations spatiales
à axe variable. Le choix du schéma historique reste inchangé par défaut.

## API et compatibilité

- `Noyau.simule(..., sigma_lie=0.0)` et
  `Noyau.audit_jacobien(..., sigma_lie=0.0)` acceptent un argument final,
  fini, compris entre 0 et 1.
- Rust expose `Modele.sigma_lie`, initialisé à zéro par le constructeur.
- `sigma_lie != 0` avec un estimateur adaptatif est refusé avant intégration.
  Le multi-rythme Python fixe σ=0 ; le multi-rythme Rust refuse σ non nul.
  Le pont adjoint temporel conserve le schéma historique.
- Les sorties, unités et conventions de repère ne changent pas.

Il s'agit d'une extension de l'API : une version mineure, suivant les
[règles du projet](VERSIONNEMENT.md). Aucun tag antérieur n'est déplacé.

## Résultat scientifique

La [note d'intégration](INTEGRATION_SIGMA.md) dérive la relation cinématique
en rotations spatiales et sa tangente implicite, y compris le traitement
distinct des colonnes GGL. Elle précise la différence avec le correcteur
approximé publié en 2025. L'implémentation ne reprend aucun code de solveur
concurrent et ne revendique pas sa simplification « sans tangente ».

La campagne compare trois réglages, quatre pas et sept cas, avec trois
répétitions et contrôle des références. Le gain d'orientation de la toupie
coexiste avec une dégradation des réactions au même pas. Le témoin plan,
la branche flexible et la boucle spatiale empêchent de conclure à un gain
général. Le choix σ reste donc explicite et expérimental.

La [confrontation externe avec MBDyn](CONFRONTATION_MBDYN_0.7.2.md)
reste attachée à la roue 0.7.2. Cette version ne fournit pas une nouvelle
mesure comparative avec MBDyn ou Simpack.

## Vérification et livraison

Les nouvelles régressions contrôlent les équations de rotation indépendantes,
les dérivées sur SO(3), les petites rotations, les changements de repère,
GGL, les valeurs invalides, le rejeu d'un pas hors carte et la restauration
après échec. L'audit des douze familles d'éléments est étendu à trois σ.

Les nombres de contrôles, les empreintes des sources et de la roue, ainsi
que les journaux de livraison figurent dans le
[relevé 0.8.0](bancs/version-0.8.0.json). Le programme de contrôle de
l'archive scientifique est intégré à la CI locale et à sa définition
GitHub conservée dans le dépôt.
