# Réactions de norme minimale par QR de la transposée

La projection statique résout désormais ses grands blocs larges de rang
plein sans SVD dense. Les réactions restent calculées en norme minimale,
avec les mêmes tolérances physiques et le même contrôle des pivots.
Cette étape prolonge la [projection par composantes](CONTRAINTES_CREUSES.md).

Sur la boucle soudée spatiale à 128 corps, la médiane passe de **2,663 s
à 2,313 ms**, soit **1 151×**, et le pic mémoire de **53,5 à 41,0 Mio**.
Il s'agit d'un calcul de réactions sur un modèle déjà immobile, sans
itération de correction de pose.

Le témoin est le commit `6a09df2079885da8540f4e8e7aaa8aa3969c7c5b`, avec
son extension conservée. Les comparaisons portent sur deux versions de
Vinkulum ; aucune nouvelle exécution MBDyn ou Simpack n'est incluse.

## Calcul

Pour un bloc `A = Gᵀ` ayant davantage de réactions inconnues que d'équations
indépendantes, on factorise sa transposée : `Aᵀ P = Q₁ R`. La solution est
`λ = Q₁ R⁻ᵀ Pᵀ f`. Les colonnes orthonormales de Q₁ donnent la solution
de norme minimale. Le principe QR/LQ pour les systèmes sous-déterminés
de rang plein est décrit dans le [guide LAPACK](https://www.netlib.org/lapack/lug/node27.html).

Le noyau effectue une seule résolution triangulaire, puis applique Q₁ par
les réflecteurs de Householder stockés par faer. Il parcourt les panneaux
et leurs sous-panneaux en ordre inverse et injecte chaque segment de la
solution aux pivots correspondants. Il ne forme ni Q, ni les équations
normales, ni un produit de la matrice originale avec une double résolution
triangulaire. Ce dernier calcul réintroduirait des soustractions de grands
termes ; un test conserve une composante exacte égale à 2 en présence
d'autres réactions de l'ordre de `1e12`.

Le QR symbolique reste réutilisé seulement à dimensions et motif identiques.
Le calcul utilise explicitement un fil, sans changer le réglage global.
Les blocs dont la plus petite dimension dépasse 32 utilisent le chemin
creux ; les petits blocs gardent le traitement précédent.

Un pivot inférieur au seuil `ε × max(dimensions) × ‖A‖F` entraîne le repli
antérieur. Une composante réellement déficiente en rang peut donc encore
utiliser une SVD dense. Le contrôle des diagonales n'est pas une garantie
universelle de révélation du rang ou de précision des réactions.

Avant de tenter le QR transposé, le noyau compte les colonnes exactement
distinctes. S'il y en a moins que de lignes, le rang plein est impossible
et la tentative est évitée. Cette borne ne déclare jamais un rang plein ;
elle n'emploie aucune tolérance sur les coefficients et ne fusionne pas
les réactions. Les autres dépendances restent contrôlées par la factorisation.

## Mesures physiques

Trois répétitions par version et par taille alternent sur le même CPU Linux,
dans des processus isolés. Le temps exclut les imports, la construction du
modèle et les contrôles indépendants. Le pic `VmHWM`, lu juste après la
statique, inclut Python et ses imports. Le mode est strict à `tol=1e-8`.

La boucle droite reprend les chaînes soudées de l'étape précédente, avec
une fermeture au sol. La boucle spatiale ajoute des orientations différentes,
des positions non colinéaires, des masses variables et une pesanteur à
trois composantes. Une référence indépendante de réactions utilise le
bilan des poids en aval et seulement six auto-contraintes de fermeture.
Les forces et moments de chaque corps sont aussi contrôlés séparément.

| Boucle | Corps | Avant | Après | Gain | Pic avant → après |
|---|---:|---:|---:|---:|---:|
| Droite | 16 | 0,420 ms | 0,427 ms | 0,98× | 39,2 → 39,2 Mio |
| Droite | 32 | 1,807 ms | 0,554 ms | 3,26× | 39,2 → 39,3 Mio |
| Droite | 64 | 13,700 ms | 0,538 ms | 25,48× | 39,7 → 39,4 Mio |
| Droite | 128 | 199,872 ms | 1,026 ms | 194,81× | 40,5 → 39,6 Mio |
| Droite | 256 | 2,868 s | 1,920 ms | 1 493,57× | 45,6 → 40,1 Mio |
| Droite | 512 | 22,539 s | 3,877 ms | 5 813,62× | 64,4 → 40,5 Mio |
| Spatiale | 16 | 2,346 ms | 0,307 ms | 7,65× | 40,3 → 40,3 Mio |
| Spatiale | 32 | 16,764 ms | 0,584 ms | 28,71× | 40,5 → 40,4 Mio |
| Spatiale | 64 | 145,015 ms | 1,179 ms | 123,01× | 43,0 → 40,5 Mio |
| Spatiale | 128 | 2,663 s | 2,313 ms | 1 151,35× | 53,5 → 41,0 Mio |

Les boucles convergent en une évaluation, sans déplacer leurs corps. Ces
gains mesurent la suppression d'une SVD dense dans la projection des
réactions ; ils ne classent pas les simulations multicorps en général.
À 128 corps spatiaux, les trois durées sont dans `[2,652 ; 2,673] s` avant
et `[2,310 ; 2,315] ms` après.

Sur cette boucle spatiale, les défauts des bilans indépendants restent
sous `6,92e-12 N` et `6,54e-12 N·m` après changement. L'écart maximal
de réaction à la référence à six inconnues est `6,78e-11` en unités SI,
contre `2,11e-10` avant ; le défaut relatif d'orthogonalité aux
auto-contraintes passe de `6,63e-15` à `2,16e-15`.

Les grandes tailles suivantes ne sont mesurées qu'avec la nouvelle version :

| Corps, boucle droite | Temps médian | Pic mémoire | Défaut d'équilibre SI | Orthogonalité relative |
|---:|---:|---:|---:|---:|
| 1 024 | 8,336 ms | 42,5 Mio | 1,12e-9 | 2,89e-12 |
| 2 048 | 17,272 ms | 46,8 Mio | 8,44e-9 | 3,09e-13 |

L'orthogonalité relative est `‖Zᵀλ‖₂/(‖Z‖F ‖λ‖₂)`. Le contrôle à
`1e-12`, tenu jusqu'à 512 corps dans la comparaison, n'est pas tenu
à 1 024 corps. Le défaut d'équilibre à 2 048 corps approche aussi la
tolérance de `1e-8`. Ces résultats ne prouvent donc pas une précision
uniforme à toutes les tailles. Ils n'ont pas de gain avant/après associé.

La chaîne articulée avec liaisons dupliquées possède de vrais degrés libres.
Les sommes des réactions sont comparées aux poids en aval, leur partage
égal entre doublons est vérifié, et chaque ressort équilibre le moment
des poids. Elle expose le coût restant du projecteur déficient en rang
et de la résolution de Newton singulière.

| Chaîne articulée avec doublons | Avant | Après |
|---:|---:|---:|
| 8 corps | 14,370 ms | 14,374 ms |
| 16 corps | 158,245 ms | 160,011 ms |
| 32 corps | 2,075 s | 2,076 s |
| 64 corps | 16,447 s | 16,560 s |

Toutes demandent quatre évaluations, avec positions compensées et
orientations identiques entre versions. Les moments et réactions sont
vérifiés indépendamment, à moins de `2,11e-11 N·m` et `5,52e-12 N`.
Cette famille **ne bénéficie pas du grand gain** des boucles soudées.
Elle reste un problème de coût prioritaire.

La version intermédiaire, avant la borne sur les colonnes distinctes,
coûtait `48,240 ms` à huit corps, contre `14,590 ms` pour son témoin.
Un diagnostic séparé avec `VINKULUM_TRACE=1` situait le surcoût dans la première tentative de QR du bloc
dupliqué. Le contrôle exact supprime cette tentative et rétablit le coût
du tableau final ; la campagne intermédiaire et son empreinte de binaire
restent archivées.

Les témoins sans ce nouveau chemin sont conservés : à 128 corps, la chaîne
soudée passe de `0,791` à `0,781 ms`, les pivots indépendants de `1,870`
à `1,880 ms` et la chaîne articulée sans doublons de `9,528` à `9,882 ms`
(+3,7 % sur cette série). Leurs états exportés sont identiques. Ces petites
variations sur trois répétitions ne suffisent pas à établir un changement
général de coût hors des blocs visés.

## Validation et reproduction

`ci/local.sh --bancs` passe : format, Clippy, **38 tests Rust**, **61 tests
Python**, **41 groupes de vérification**, **46 bancs rapides**, **9 bancs
de contact**, référence d'API à **185 entrées**.

Les tests ajoutés couvrent le QR transposé sur matrices creuses et denses,
les permutations, le changement des valeurs et du motif, une direction
libre conservant son résidu, les directions presque dépendantes et la
répartition analytique d'une boucle mécanique redondante.

Les **108 équilibres** des mesures finales tiennent le mode strict et les
contrôles physiques indépendants. Les positions compensées et orientations
comparées entre versions sont identiques. Les **400 paliers Princeton**
passent pour les deux formulations, ainsi que les trois témoins
d'allongement exact `2^-60`. Le câble condensé extrême reste en échec
depuis ses deux initialisations ; l'élément historique conserve ses succès.

```bash
.venv314/bin/python ci/mesure_redondances.py \
  --ancien-pythonpath /chemin/vers/le/temoin \
  --sortie docs/bancs/redondances-mesures.json
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_positions_compensees.py \
  --sortie docs/bancs/redondances-princeton-strict.json
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_cable.py \
  --sortie docs/bancs/redondances-cable.json
# Après journalisation de ci/local.sh --bancs :
.venv314/bin/python ci/valide_redondances.py --journal-ci /tmp/redondances-ci.log
```

[Échantillons complets](bancs/redondances-mesures.json),
[comparaison calculée](bancs/redondances-comparaison.json),
[400 paliers Princeton](bancs/redondances-princeton-strict.json),
[câble](bancs/redondances-cable.json),
[version intermédiaire](bancs/redondances-intermediaire.json),
[trace du témoin](bancs/redondances-double8-avant.trace),
[trace intermédiaire](bancs/redondances-double8-intermediaire.trace),
[trace après correction](bancs/redondances-double8-apres.trace),
[journal CI](bancs/redondances-ci.log),
[tests Python](bancs/redondances-tests-python.log),
[empreintes et validation](bancs/redondances-validation.json).

Le travail restant comprend les rangs déficients, le Newton singulier,
les autres projecteurs et les matrices publiques. Les limites de la
formulation intégrée du câble et les comparaisons industrielles restent
ouvertes. L'[objectif généraliste](OBJECTIF_MBDYN.md) demeure actif.
