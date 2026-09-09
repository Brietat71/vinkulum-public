# Vinkulum 0.12.1 — assemblage et initialisation dynamique

Cette corrective traite des défauts reproduits dans le noyau 0.12.0.
Les signatures publiques sont conservées. Elle porte sur l'assemblage,
la projection des vitesses, leur diagnostic temporel et l'accélération
initiale de l'intégrateur α-généralisé.

## Défauts corrigés

| Situation | Comportement 0.12.0 reproduit | Comportement corrigé |
|---|---|---|
| Deux commandes imposent 1 et 2 m/s au même degré de liberté | Succès avec 1,5 m/s, résidus ±0,5 m/s | Refus explicite et restitution de l'état avant l'appel |
| Liaison à bras de 10⁵ m, translation et rotation bloquées | Vitesse parasite de −2,14 × 10⁻⁶ m/s après projection d'une vitesse de 1 m/s | Contraintes satisfaites à la tolérance demandée |
| Même liaison avec un bras de 10⁸ m | Translation presque inchangée et rotation parasite de −10⁻⁸ rad/s | Les deux vitesses bloquées sont corrigées |
| Une correction de position suffit, `iters=1` | Faux échec ; la dernière pose n'était pas contrôlée | Succès avec une correction |
| Commande par table, pentes 1 puis 3 m/s | 2 m/s au nœud intérieur ; demi-pente aux extrémités | Convention analytique de la loi, commune au diagnostic et à la projection |
| Translation non holonome pilotée, assemblage à `t=1` | Échec en cherchant une position qui n'est pas imposée | Correction des seules positions holonomes, puis de toutes les vitesses requises |
| Diagnostic d'une vitesse non holonome déjà satisfaite | Fausse alerte sur une position non imposée | Le contrôle de position porte sur les seules contraintes holonomes |
| Débordement d'une dérivée de commande | Possibilité de succès avec une vitesse non finie | Refus avant écriture des vitesses et restitution de la pose |

La dérivée analytique supprime aussi les erreurs de différence temporelle
pour les rotations rapides et les grandes dates. Le diagnostic d'une table
constante hors de son domaine reste défini jusqu'à la plus grande date
finie ; il ne requiert plus une date perturbée au-delà de cette limite.

## Résolution et contrat d'échec

### Accélération initiale : défaut physique révélé par la campagne

Une campagne indépendante, ajoutée après les premières régressions et une
CI étendue verte, a trouvé une erreur supplémentaire. Un plafonnement de
l'accélération initiale en fonction du pas altérait l'état consistant calculé
par le noyau. Sur un pendule de 1 000 m, la dérive relative d'énergie atteignait
0,23 % en 256 pas, contre une référence d'EDO indépendante. Ce comportement
est reproduit sur la roue 0.12.0 et sur la première candidate corrective.

La contre-épreuve minimale est la chute libre depuis le repos, pendant une
seconde sous 9,81 m/s². Avec un pas de 1 s, la 0.12.0 donnait −2,482819 m
au lieu de −4,905 m ; avec un pas de 0,1 s, −4,676626 m. Le schéma à
accélération constante doit reproduire cette solution quadratique.

Ce plafonnement est supprimé. Les accélérations physique et algorithmique
initiales restent consistantes ; les difficultés de Newton sont traitées
par les reprises et subdivisions existantes. Les régressions couvrent les
forces constantes, la gravité et un couple constant autour d'un axe principal.
Cette correction concerne bien un aspect critique de l'intégration temporelle.

### Assemblage et projection

La correction résout directement le système rectangulaire `G δ = r`,
avec équilibrage des lignes et la factorisation orthogonale QR/SVD déjà
présente dans le noyau. La matrice normale `GGᵀ`, qui mettait le
conditionnement au carré, n'est plus formée dans ces deux chemins.
Pour des équations compatibles, l'équilibrage préserve l'ensemble des
solutions et la correction de norme minimale dans les coordonnées utilisées.

La projection contrôle **chaque équation originale**, y compris les lignes
redondantes écartées de la factorisation dynamique. Une solution de moindres
carrés ne suffit plus à annoncer un assemblage réussi. Le contrôle utilise
un résidu accumulé avec compensation et, pour chaque ligne, la limite

```text
|G_i u + Φ_t,i| ≤ tol + 64 ε (|Φ_t,i| + Σ_j |G_ij u_j|).
```

Au plus trois corrections sont tentées. Un échec numérique, une précision
insuffisante ou des vitesses incompatibles produisent une exception.
La sauvegarde de l'appel Python couvre désormais **les deux phases** :
positions, rotations, petites translations compensées, vitesses, horloge,
schéma et multiplicateurs sont restitués après un refus. `vitesses=False`
conserve son usage de résolution des seules positions.

Les corrections ne déplacent plus les corps gelés d'une partition ; les
équations et contrôles de bassin des éléments gelés sont réservés à leur
propre sous-pas. Le masque de redondance n'est pas un droit à violer une
contrainte physique.

## Validation et limites

### Portée scientifique des défauts

Les faux succès remettent en cause le **critère d'admissibilité d'un état**
dans les chemins concernés. Ils peuvent affecter la préparation d'une
simulation et sa reprise après changement du masque de redondance
(`Modele::masque_rejoue`). L'assemblage est également appelé par les imports
URDF/MJCF et plusieurs campagnes, dont le mécanisme d'Andrews. Le problème
n'est donc pas limité à un exemple utilisateur isolé.

Le défaut de grand bras est explicable sans comparaison entre solveurs.
Sur les deux coordonnées `(vy, wz)`, les équations donnent

```text
G = [[1, L], [0, 1]],      G Gᵀ = [[1 + L², L], [L, 1]].
```

Le déterminant exact de la seconde matrice vaut 1. Pour `L = 10⁸`,
`1 + L²` est arrondi à `L²` en binary64 : la formation de cette matrice
efface l'information qui sépare les deux contraintes. Le conditionnement
de G croît comme L² et celui de GGᵀ comme L⁴. Équilibrer les lignes puis
travailler sur G conserve ici une séparation exploitable ; cela ne rend
pas toutes les géométries arbitrairement bien conditionnées.

L'erreur d'accélération initiale remet aussi en cause les résultats temporels
pour lesquels le plafonnement était actif. Les observations ne démontrent
pas une erreur générale des équations du mouvement ou des tangentes. Elles ne suffisent pas non
plus à certifier les trajectoires et réactions antérieures. Les résultats
qui utilisent ces chemins doivent être rejugés sur leurs références : une
CI verte couvre les cas exécutés, pas tous les modèles admissibles. Les
archives historiques conservent leur version d'origine ; aucun résultat
de classement n'est réattribué à la corrective.

### Contrôles exécutables

Les régressions cinématiques sont dans
[`test_assemblage.py`](../python/vinkulum/test_assemblage.py) : solutions
fermées de projection, grands bras de levier, commandes temporelles,
échecs répétés et restauration après une simulation ayant calculé son schéma
et ses réactions. Deux tests Rust contrôlent les corps et éléments gelés.
Le contrôle des grandes dates de `test_analyses.py` est adapté à la dérivée
analytique ; les autres attentes mécaniques existantes sont conservées.

La campagne [`audit_assemblage_physique.py`](../ci/audit_assemblage_physique.py)
fixe ses critères avant les mesures : 2 304 projections dans 24 repères
orthogonaux exacts, six longueurs, deux origines, permutations et duplications
de contraintes ; deux explorations au-delà du domaine principal ; six
trajectoires de pendule sur trois échelles, comparées à DOP853 à deux
tolérances. Les seuils initiaux ne sont pas relevés après l'échec du pendule.

Les 2 304 projections ont une erreur maximale de 4,45 × 10⁻¹⁶ sur les
coordonnées contrôlées. Les deux explorations à des bras supérieurs à 10¹² m
sont refusées, avec restitution de l'état. Sur le pendule de 1 000 m à
256 pas, l'erreur adimensionnelle finale passe de 9,83 × 10⁻⁴ à
3,93 × 10⁻⁵, et la dérive relative d'énergie de 2,30 × 10⁻³ à
5,79 × 10⁻⁸. Aux trois longueurs, doubler le nombre de pas divise l'erreur
par environ quatre. Les références resserrées diffèrent de moins de
3,6 × 10⁻¹³. Ces observations portent sur les problèmes et grilles déclarés.

Le [bilan de livraison](bancs/version-0.12.1.json) conserve les résultats,
les empreintes et les journaux de validation, avec une comparaison des
nouvelles régressions à la roue 0.12.0 figée.

Cette livraison ne démontre pas la fiabilité de tout modèle multicorps.
L'assemblage reste dense, le seuil de rang des équations équilibrées vaut
10⁻¹² et une géométrie trop mal conditionnée peut être refusée. La tolérance
continue de porter sur des composantes exprimées dans les unités du modèle ;
elle ne constitue pas un critère universel sans dimension. Aux nœuds d'une
table, la pente gauche est une convention de la loi, pas une dérivée
mathématique bilatérale. Aucun gain de vitesse face à un autre solveur
n'est revendiqué pour cette corrective.
