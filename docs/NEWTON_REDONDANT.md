# Newton statique sur les contraintes équivalentes

La chaîne mobile de 64 corps avec pivots dupliqués passe de **16,702 s à
7,159 ms**, soit **2 333×**, avec quatre évaluations de Newton dans les
deux versions. Le pic mémoire passe de **72,5 à 40,5 Mio**. Les réactions
individuelles sont conservées et vérifiées par des bilans indépendants.

La réduction corrige aussi une perte de convergence : la poutre Princeton
à charge pleine converge maintenant avec un encastrement dupliqué. Le
témoin échoue pour les deux formulations et les deux maillages mesurés.

Le témoin est `d93593cc751d3aeaf523c0c11507d6a4caf4143f`, avec son
extension conservée. Ces résultats comparent deux versions de Vinkulum.
Aucune nouvelle exécution MBDyn ou Simpack n'entre dans cette campagne.

## Fondement algébrique et portée de la preuve

Pour un groupe de `k` contraintes, on exige exactement
`G_i = s_i g` et `Φ_i = s_i h`, avec `s_i ∈ {−1,+1}`. On définit une
colonne de U par `U_i = s_i/√k`, nulle en dehors du groupe. Les groupes
sont disjoints ; les lignes non regroupées ont leur propre colonne.
Ainsi, en arithmétique réelle :

```text
UᵀU = I             G = U G_r             Φ = U Φ_r
(G_r)_groupe = √k g                        (Φ_r)_groupe = √k h
```

Toute réaction s'écrit `λ = U μ + z`, avec `Uᵀz = 0`. Alors
`Gᵀλ = G_rᵀμ` et `‖λ‖₂² = ‖μ‖₂² + ‖z‖₂²`. La réaction de norme
minimale a donc `z = 0`. Résoudre les moindres carrés réduits et reconstruire
`λ = U μ` conserve à la fois la projection des forces et la norme minimale.
Cette propriété reste vraie si le gradient réduit possède d'autres
dépendances, à condition que son solveur fournisse lui aussi la solution
de norme minimale.

Le même changement de coordonnées s'applique à Newton. Avec
`T = diag(I,U)`, son système complet et son système réduit vérifient :

```text
A = [[K, Gᵀ], [G, 0]] = T A_r Tᵀ
A_r = [[K, G_rᵀ], [G_r, 0]]
b = [r_f; −Φ] = T [r_f; −Φ_r]
```

Comme `TᵀT = I`, les conditions de Moore–Penrose donnent
`A⁺ = T A_r⁺ Tᵀ`. La résolution réduite, suivie du relèvement des
incréments de réaction, produit donc la même solution de norme minimale
que la SVD du système complet. Cette identité ne suppose pas K symétrique.
L'équilibrage de Jacobi utilisé ici ne change que les coordonnées
physiques ; il commute avec cette injection des multiplicateurs.

K est calculé **après reconstruction de toutes les réactions**, car il
contient la dérivée de `Gᵀλ`. Une réaction affectée uniquement au premier
exemplaire modifierait sa norme et pourrait modifier la précontrainte
si des fonctions n'ont la même tangente qu'à la configuration courante.

Cette preuve repose sur l'orthogonalité et la décomposition en somme
directe, des outils classiques d'algèbre linéaire. La conservation de la
norme par transformation orthogonale est rappelée dans le
[guide LAPACK](https://www.netlib.org/lapack/lug/node39.html).
Les tests du noyau comparent notamment un système augmenté non symétrique
et équilibré à sa SVD complète. Ils contrôlent aussi une réaction
manufacturée appartenant à l'image de G, des multiplicités inégales,
des signes opposés et une dépendance supplémentaire non supprimée.
Ces contrôles quantifient l'effet des arrondis ; la preuve en arithmétique
réelle ne garantit pas à elle seule une précision uniforme en flottants.

## Traitement numérique

La détection construit temporairement les lignes du gradient au format
CSR. Les clés empruntent ces lignes et comparent les coefficients exacts,
avec signe canonique, ainsi que l'écart de contrainte signé. Une collision
de hachage n'entraîne pas de fusion sans égalité des coefficients.
Le coût attendu de cette détection et sa mémoire temporaire sont linéaires
en `nnz(G)+m`. La mise en ordre du gradient réduit utilise ensuite le tri
habituel des triplets.

Il n'y a aucun seuil de fusion. Deux directions séparées de `1e-12`,
ou deux écarts différents, restent distincts. Les lignes nulles restent
également distinctes. Un dépassement de capacité lors de la mise à
l'échelle abandonne la compression et conserve le chemin antérieur.
Les groupes sont recalculés à chaque évaluation de Newton. Les caches
de factorisation continuent de contrôler dimensions et motif.

Toutes les contraintes originales restent dans le modèle, les réactions
publiques et le critère final. Le mérite physique de la recherche de pas
les utilise aussi. Pour le mérite naturel, un essai dont les écarts
ne respectent plus le regroupement courant est refusé : deux fonctions
ayant initialement la même valeur et tangente peuvent ensuite diverger.
Ce refus conservateur ne constitue pas une preuve générale de convergence.

Le solveur réduit peut réutiliser la factorisation LU pour le mérite
naturel. Le système complet singulier empruntait auparavant le repli SVD,
qui désactivait ce mérite. Cela explique le changement de verdict du
témoin Princeton avec encastrement double.

## Protocole

Trois répétitions alternent les versions dans des processus séparés,
sur le même CPU Linux, avec un fil demandé à Rayon, OpenBLAS et OpenMP.
Les campagnes générales et les compilations sont terminées avant ces
mesures. Le temps exclut les imports, la construction des modèles et
les contrôles indépendants. `VmHWM` est relevé immédiatement après la
statique ; il inclut Python et les imports déjà effectués.

Les essais utilisent `strict=True`, `tol=1e-8`, `iters=100` et
`paliers_max=1`. L'échelle de force reste celle du déséquilibre libre
initial, avec plancher 1. Les échecs sont archivés avec leur état restauré.
Les grandes tailles ne disposant que d'une mesure actuelle sont indiquées
comme telles. Aucun gain n'est calculé entre un échec et un succès.

## Mécanismes mobiles

La chaîne mesure 1 m pour 1 kg ; chaque pivot est associé à un ressort
de torsion de raideur `k = 10 × N`, en N·m/rad, où N est le nombre de corps. La
famille « double » possède deux exemplaires de chaque pivot. La famille
« mixte » alterne un, deux et trois exemplaires, avec des repères opposés.

| Famille | Corps | Avant | Après | Gain | Pic avant → après |
|---|---:|---:|---:|---:|---:|
| Double | 8 | 14,492 ms | 0,965 ms | 15,02× | 39,8 → 39,8 Mio |
| Double | 16 | 158,780 ms | 1,833 ms | 86,63× | 41,5 → 40,1 Mio |
| Double | 32 | 2,096 s | 3,568 ms | 587,44× | 48,0 → 40,2 Mio |
| Double | 64 | 16,702 s | 7,159 ms | 2 333,03× | 72,5 → 40,5 Mio |
| Mixte | 16 | 106,652 ms | 1,834 ms | 58,14× | 42,1 → 40,2 Mio |
| Mixte | 32 | 627,233 ms | 3,570 ms | 175,71× | 47,6 → 40,3 Mio |
| Mixte | 64 | 4,502 s | 7,075 ms | 636,32× | 72,0 → 40,8 Mio |

Toutes demandent quatre évaluations et une tentative. À 64 corps doubles,
les trois durées sont dans `[16,652 ; 16,745] s` avant et environ
`[7,134 ; 7,170] ms` après. Les données exactes sont dans la comparaison
calculée. Les plus grands écarts de pose entre versions, sur ces chaînes,
sont inférieurs à `1,30e-13 m` et `1,22e-13` sur les matrices de rotation.

Les références additionnent les poids en aval dans le repère de chaque
pivot et calculent leur moment à partir des positions compensées. À
64 corps doubles, le défaut des réactions totales est `8,00e-15 N`, celui
des moments `2,77e-13 N·m`, et la fermeture tient `8,28e-14 m`. Les deux
copies reçoivent exactement la même réaction dans les sorties mesurées.
Le test de multiplicité mixte vérifie chaque réaction individuelle,
après changement de repère, et chaque ressort de torsion.

| Grandes tailles, version actuelle seule | Temps médian | Pic mémoire |
|---|---:|---:|
| Double, 128 corps | 14,857 ms | 42,0 Mio |
| Double, 256 corps | 30,265 ms | 44,4 Mio |
| Double, 512 corps | 62,605 ms | 49,2 Mio |
| Double, 1 024 corps | 128,207 ms | 59,5 Mio |
| Mixte, 128 corps | 14,697 ms | 42,0 Mio |
| Mixte, 256 corps | 30,636 ms | 44,5 Mio |
| Mixte, 512 corps | 61,649 ms | 49,2 Mio |

À 1 024 corps doubles, les défauts indépendants restent sous `2,71e-13 N`
pour les réactions, `1,20e-12 N·m` pour les ressorts et `3,87e-15 m` pour
les fermetures. Cette croissance mesurée du coût ne garantit pas la même
pente pour d'autres topologies.

## Convergence de Princeton

La charge de 8,896 N, inclinée entre les axes principaux de la poutre,
est appliquée entièrement dès le premier appel. Le seul changement de
modèle est un second exemplaire de l'encastrement. Le témoin échoue après
111 évaluations cumulées sur deux tentatives, puis restaure l'état initial.
Les temps jusqu'à cet échec restent archivés, sans ratio de vitesse.

| Formulation | Intervalles | Témoin | Version actuelle | Évaluations actuelles |
|---|---:|---|---:|---:|
| Milieu | 10 | Échec | Tolérance en 17,815 ms | 47 |
| Milieu | 20 | Échec | Tolérance en 40,983 ms | 48 |
| Intégrée | 10 | Échec | Tolérance en 19,058 ms | 45 |
| Intégrée | 20 | Échec | Tolérance en 42,833 ms | 45 |

Les réactions d'encastrement sont comparées à la force extérieure et à
son moment au bout, avec partage égal entre les deux copies. Le test
automatique compare aussi les poses avec celles de la poutre à un seul
encastrement. Rejoué avec l'ancienne extension, ce test échoue par absence
de convergence ; son journal est conservé.

## Témoins et limites mesurées

Chaque parallélogramme fermé possède trois corps mobiles, quatre pivots,
vingt équations et une mobilité. Les trois dépendances de son gradient
ne sont pas des lignes équivalentes. Sa référence de pose provient de
l'équation scalaire d'énergie
`100 δ + 2 × 9,81 cos(π/3 + δ) = 0` ; les forces et moments sont vérifiés
séparément sur chacun des corps. Les états exportés restent identiques
entre les deux versions.

| Parallélogrammes indépendants | Avant | Après |
|---:|---:|---:|
| 1 | 0,577 ms | 0,605 ms |
| 4 | 8,325 ms | 8,460 ms |
| 8 | 120,588 ms | 119,917 ms |

Cette famille garde le coût du système augmenté singulier. Elle fournit
un témoin concret pour généraliser le traitement du rang. Le plus grand
écart à la pose analytique est inférieur à `1e-11 m` ; les bilans de force
et de moment restent sous `2,08e-12 N` et `8,22e-11 N·m`.

La détection conserve un surcoût sur les familles sans lignes équivalentes :

| Témoin | Avant | Après | Variation |
|---|---:|---:|---:|
| Chaîne articulée, 8 corps | 0,703 ms | 0,781 ms | +11,0 % |
| Chaîne articulée, 128 corps | 9,838 ms | 10,454 ms | +6,3 % |
| Chaîne soudée, 128 corps | 0,787 ms | 0,893 ms | +13,4 % |
| Rotors, 128 corps | 1,840 ms | 1,977 ms | +7,4 % |
| Boucle soudée spatiale, 128 corps | 2,308 ms | 2,622 ms | +13,6 % |
| Boucle soudée droite, 512 corps | 3,799 ms | 4,401 ms | +15,9 % |

Leurs poses sont identiques entre versions. Les mesures portent sur trois
répétitions ; leurs étendues sont publiées. Le coût de détection hors des
cas accélérés demeure à réduire.

## Vérification et reproduction

La campagne générale `ci/local.sh --bancs` passe : format, Clippy,
**42 tests Rust**, **41 groupes de vérification**, **46 bancs rapides**,
**9 bancs de contact** et référence d'API à **185 entrées**. Le test
Princeton à encastrement double a été ajouté ensuite ; la suite Python
complète a été rejouée et ses **64 tests** passent.

Les **141 essais comparatifs** comprennent **129 équilibres stricts** et
les **12 échecs Princeton du témoin**. Les **400 paliers Princeton** sans
doublons passent encore pour les deux formulations, ainsi que les trois
témoins d'allongement `2^-60`. Le câble condensé extrême reste en échec
depuis ses deux initialisations ; la formulation historique réussit.

```bash
.venv314/bin/python ci/mesure_newton_redondant.py \
  --ancien-pythonpath /chemin/vers/le/temoin \
  --sortie docs/bancs/newton-redondant-mesures.json
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_positions_compensees.py \
  --sortie docs/bancs/newton-redondant-princeton-strict.json
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_cable.py \
  --sortie docs/bancs/newton-redondant-cable.json
.venv314/bin/python ci/valide_newton_redondant.py \
  --journal-ci /chemin/vers/ci.log \
  --journal-tests /chemin/vers/tests-python.log \
  --journal-negatif /chemin/vers/princeton-negatif.log
```

[Échantillons complets](bancs/newton-redondant-mesures.json),
[comparaison calculée](bancs/newton-redondant-comparaison.json),
[400 paliers stricts](bancs/newton-redondant-princeton-strict.json),
[câble](bancs/newton-redondant-cable.json),
[journal CI](bancs/newton-redondant-ci.log),
[64 tests Python](bancs/newton-redondant-tests-python.log),
[échec sur le témoin](bancs/newton-redondant-princeton-negatif.log),
[empreintes et validation](bancs/newton-redondant-validation.json).

La réduction générale des dépendances, les changements de rang en cours
de mouvement, les autres projecteurs et les comparaisons industrielles
restent à traiter. L'[objectif généraliste](OBJECTIF_MBDYN.md), avec recours
massif aux mathématiques fondamentales, reste actif.
