# Dépendances générales : factorisation orthogonale et contrôle de la SVD

Le noyau sait désormais calculer les moindres carrés de norme minimale sur
des dépendances combinant plusieurs contraintes. Cette étape concerne la
projection statique des forces et le repli du système augmenté de Newton.
Elle complète la [réduction des lignes équivalentes](NEWTON_REDONDANT.md).
Les méthodes employées sont classiques ; leur intégration, les décisions
de rang et leurs effets mécaniques sont vérifiés ci-dessous.

## Calcul dans les deux espaces orthogonaux

Pour le système transmis au solveur, `A P = Q R`. Après sélection du rang
`r`, on factorise la transposée des `r` premières lignes : `R_rᵀ = Z T`.
La solution est alors

```text
x = P Z [ T⁻ᵀ (Qᵀ b)_r ; 0 ].
```

En arithmétique exacte et au rang exact, la première transformation
sépare la partie du second membre atteignable par `A`. La seconde sépare
l'espace des lignes du noyau droit ; annuler la composante dans ce noyau
donne la solution de norme minimale. Un QR avec pivotage seul ne suffit
pas à cette dernière propriété : voir le
[guide LAPACK sur la factorisation orthogonale complète](https://www.netlib.org/lapack/lug/node43.html).

Les deux transformations sont appliquées par réflecteurs de Householder,
sans former `Q` ou `Z`. Aucun produit `AᵀA` n'est construit. La matrice est
normalisée par un scalaire commun à `A`, `b` et au seuil. Les opérations
faer utilisent explicitement un seul fil, sans modifier le parallélisme
global de l'application. La mémoire et le calcul restent denses dans ce
repli ; le QR avec pivotage de la bibliothèque n'est pas présenté comme
un algorithme intégralement par blocs.

La norme minimale concerne les coordonnées du système fourni : les
réactions pour la projection, les variables équilibrées par Jacobi pour
Newton, comme dans son ancien repli SVD. Toutes les réactions originales
sont reconstruites avant la raideur de précontrainte et les résidus.

## Rang numérique et limites des garanties

Le seuil de chaque appel est conservé. Le rang proposé par la diagonale
du QR n'est accepté que si deux contrôles le séparent nettement du seuil
`τ` :

- la norme de Frobenius majorée du triangle restant est au plus `τ/4` ;
- une borne inférieure de `σ_min(R_11)` dépasse `4τ`.

Pour cette seconde borne, les résolutions triangulaires avec diagonale
`|R_ii|` et coefficients hors diagonale `−|R_ij|` majorent les sommes
absolues des lignes et colonnes de `R_11⁻¹`. Leur produit fournit
`σ_min(R_11) ≥ 1/√(B_1 B_inf)`. Les opérations positives de la majoration
sont arrondies vers le haut, puis la borne finale vers le bas.

Cette borne concerne le **triangle calculé**. Elle ne constitue pas une
preuve par intervalles du rang de la matrice d'entrée. Elle peut être très
pessimiste ; une séparation insuffisante conserve le repli SVD. Le
[guide LAPACK sur le pivotage](https://www.netlib.org/lapack/lug/node42.html)
distingue lui aussi la factorisation de la décision du rang numérique.

La stationnarité `Aᵀ(Ax−b)` est ensuite contrôlée sur la matrice d'entrée
normalisée, sans former de matrice normale. Le seuil de ce contrôle tient
compte de la précision machine, des dimensions et des normes des données.
Les sous-flux lors de la normalisation et les valeurs non finies entraînent
un repli ou un refus explicite. Les résidus physiques complets décident
toujours de la convergence statique ; les tolérances mécaniques restent
inchangées.

## Erreur détectée dans l'ancien repli SVD

Une matrice fabriquée `7×5`, de rang numérique 3, a révélé une erreur
dans l'appel nalgebra 0.33.3 utilisé précédemment. On choisit `x_ref=Aᵀy`
puis `b=A x_ref` : la référence de norme minimale est connue sans SVD.
La troisième valeur singulière était erronée, et pas seulement une valeur
proche du seuil de coupure.

`svd_sure` utilise maintenant la SVD mince de faer, avec refus préalable
des entrées non finies, vérification des facteurs et propagation des
erreurs de convergence. Son type de retour reste compatible avec les
appels existants, dont la dynamique et la remise des rotations sur SO(3).
La reconstruction, l'orthogonalité, la solution fabriquée, les matrices
vides et nulles, et les échelles `10⁻²⁰⁰` à `10²⁰⁰` sont testées.

Le programme `examples/diagnostic_svd.rs` exporte les matrices et les
facteurs des deux méthodes. `ci/valide_svd.py` les confronte à la référence
fabriquée et au pilote LAPACK `gesvd` de SciPy.

| Contrôle sur ce cas | Ancien appel nalgebra | Nouvel appel faer |
|---|---:|---:|
| Erreur de solution, norme 2 | 2,384 × 10⁻³ | 1,920 × 10⁻¹⁵ |
| Erreur relative de reconstruction | 2,390 × 10⁻² | 6,775 × 10⁻¹⁶ |
| Troisième valeur singulière | 1,8264947021 | 1,8227188893 |

LAPACK retrouve `1,8227188893` et une erreur de solution de
`3,152 × 10⁻¹⁶`. Les données et facteurs complets sont dans
[l'archive numérique](bancs/factorisations-svd-lapack.json). Ces chiffres
décrivent ce cas fabriqué, sans constituer une précision universelle.

## Préparation des petits calculs

Le premier prototype préparait aussi des réflecteurs par blocs pour les
petites matrices. Sur un seul parallélogramme, le premier appel coûtait
environ 35 ms.
Le traçage des ouvertures de fichiers a révélé 2 880 accès aux informations
de cache des processeurs lors de cette première préparation.

Jusqu'à 64 lignes et colonnes, le calcul applique maintenant les
réflecteurs individuellement. Il conserve le QR avec pivotage, les bornes
de rang et la seconde transformation orthogonale. Les accès aux fichiers
de cache disparaissent sur ce petit cas. Le
[rapport de traçage](bancs/factorisations-demarrage.json) et ses journaux
compressés conservent cette vérification ; les temps instrumentés ne
servent pas à comparer les vitesses. Les mesures principales incluent le
premier appel dans un processus neuf.

## Référence mécanique indépendante

La cascade comporte `N` parallélogrammes reliés entre eux, soit `3N`
corps, `20N` contraintes, `N` mobilités et `3N` dépendances. Le coupleur
d'une cellule porte la suivante ; le graphe est connexe. Chaque bras a
une longueur `L=1/N`, chaque corps une masse `m=1/(3N)` et chaque cellule
un ressort de torsion `k=100/N`. L'angle initial est `α=π/3`.

L'énergie réduite vaut, avec les cellules numérotées de zéro à `N−1`,

```text
V = Σ_i [ m g L (3(N−i)−1) sin(α+δ_i) + k δ_i²/2 ].
```

Les angles d'équilibre sont donc les racines indépendantes de
`kδ_i + m g L(3(N−i)−1)cos(α+δ_i)=0`. L'énergie est strictement convexe
sur l'intervalle testé `[-0,4 ; 0]`. Le banc contrôle les poses obtenues
par ces racines, la fermeture des articulations, les efforts hors du plan
et l'équilibre des forces et moments de chaque corps aux points d'ancrage.
Ces calculs n'utilisent ni les gradients ni les tangentes du noyau.

## Résultats et surcoûts conservés

La comparaison finale utilise le commit `6d38aa4` comme témoin, avec la
même tolérance stricte, dans des processus neufs alternés. Les temps
ci-dessous sont les médianes de trois appels ; leurs étendues, réactions
et poses sont dans le [bilan complet](bancs/factorisations-mesures-bilan.json).

| Modèle | Taille | Avant, ms | Après, ms | Évaluations avant → après |
|---|---:|---:|---:|---:|
| Cascade | 1 cellule | 0,614 | 0,519 | 4 → 4 |
| Cascade | 4 cellules | 36,430 | 40,181 | 20 → 5 |
| Parallélogrammes indépendants | 1 cellule | 0,606 | 0,505 | 4 → 4 |
| Parallélogrammes indépendants | 8 cellules | 119,667 | 57,404 | 4 → 4 |
| Chaîne à liaisons doubles | 64 corps | 7,006 | 7,105 | 4 → 4 |
| Chaîne à multiplicité variable | 64 corps | 7,118 | 7,191 | 4 → 4 |
| Chaîne articulée | 128 corps | 10,146 | 11,006 | 4 → 4 |
| Boucle spatiale | 128 corps | 2,545 | 2,521 | 1 → 1 |
| Chaîne soudée | 128 corps | 0,908 | 0,895 | 1 → 1 |
| Rotors | 128 corps | 2,022 | 2,133 | 2 → 2 |
| Boucle soudée | 512 corps | 4,422 | 4,403 | 1 → 1 |
| Princeton milieu, encastrement double | 10 intervalles | 17,807 | 19,074 | 47 → 47 |
| Princeton intégrée, encastrement double | 10 intervalles | 19,094 | 20,290 | 45 → 45 |

Le gain sur huit parallélogrammes vaut **2,085×**. La correction des petits
calculs ramène le premier parallélogramme sous le temps du témoin. Des
surcoûts restent mesurés : **10,3 %** sur la cascade à quatre cellules,
**8,5 %** sur la chaîne articulée, et **6,3–7,1 %** sur Princeton. La
préparation des produits par blocs reste utilisée sur les grandes
matrices ; l'accélération n'est pas universelle.

Les échecs du témoin sont conservés dans la
[campagne initiale](bancs/factorisations-blocs-bilan.json). Les cas concernés
ont ensuite été vérifiés avec le binaire final, sans recalculer un ratio
de vitesse à partir d'un échec :

| Modèle | Ancien témoin | Binaire final : temps / évaluations / RSS |
|---|---|---|
| Cascade, 8 cellules / 24 corps | Échec, 111 évaluations | 70,56 ms / 5 / 41,5 Mio |
| Cascade, 16 cellules / 48 corps | Échec | 249,90 ms / 6 / 50,5 Mio |
| Cascade, 32 cellules / 96 corps | Non exécuté dans cette comparaison | 2,038 s / 9 / 86,0 Mio |
| 16 parallélogrammes indépendants | Échec | 156,43 ms / 4 / 52,0 Mio |
| 32 parallélogrammes indépendants | Délai de processus de 90 s dépassé | 801,81 ms / 4 / 87,0 Mio |

Sur la cascade à 32 cellules, l'écart de position à la référence d'énergie
est inférieur à `6,4 × 10⁻¹³ m`, le bilan des forces à `2,3 × 10⁻¹³ N`
et celui des moments à `2,2 × 10⁻¹² N·m`. Les
[mesures aux grandes tailles](bancs/factorisations-grands-bilan.json)
conservent tous les autres contrôles. La croissance du coût reste nette
avec la taille : cette étape ne fournit pas encore un repli creux pour
les dépendances générales.

## Contribution de chaque changement

L'[ablation](bancs/factorisations-ablation-bilan.json) compare deux prototypes
qui possèdent tous deux la nouvelle SVD. La SVD seule échoue sur les
cascades à 4, 8, 16 et 32 cellules dans le budget fixé. L'association avec
la factorisation orthogonale converge sur ces quatre tailles. Sur les
parallélogrammes indépendants à 8, 16 et 32 cellules, les gains médians
supplémentaires valent respectivement **1,48×, 2,05× et 2,57×**. Pour
32 cellules, le RSS passe de **171,5 à 87,0 Mio** sur ces essais convergés.

La correction de la SVD ne suffit donc pas à expliquer la convergence
obtenue. Cette comparaison ne constitue pas une garantie de convergence
pour tous les mécanismes. L'étendue temporelle plus large de la troisième
répétition est conservée, avec les rapports calculés pour chaque paire.

Les quatre campagnes totalisent **255 essais** : **231 équilibres stricts**,
**21 échecs avec état restauré** et **3 délais de processus dépassés**.
Les deux campagnes portant sur le binaire final totalisent 93 équilibres
stricts, dont 54 calculés avec ce binaire et 39 avec le témoin.
Les 48 tests Rust, 65 tests Python, 41 cas généraux, 46 bancs mécaniques,
9 cas de contact et 400 paliers Princeton passent. Le câble condensé
extrême conserve son échec connu. Le
[manifeste](bancs/factorisations-validation.json) lie les sources, binaires,
journaux et résultats.

## Reproduction

```bash
cargo run --release --example diagnostic_svd > /tmp/svd-rust.json
python ci/valide_svd.py /tmp/svd-rust.json /tmp/svd-lapack.json
python ci/diagnostic_cascade.py 8
ci/local.sh --bancs
python ci/valide_factorisations.py
```

Les comparaisons de binaires utilisent `ci/mesure_svd.py`, avec des
paquets Python contenant une copie distincte de chaque extension. Les
processus sont alternés, fixés sur le même CPU, à un fil, avec trois
répétitions, la même tolérance stricte `1e-8` et le même budget. Les états,
réactions, erreurs, empreintes et étendues temporelles sont conservés.
Un échec de convergence n'est jamais converti en gain de vitesse.

Les patches `factorisations-ablation-svd.patch`,
`factorisations-ablation-cod.patch` et `factorisations-blocs.patch`
reconstituent les sources des prototypes depuis cette version. Ils sont
accompagnés de leurs empreintes et s'appliquent avec
`git apply --unidiff-zero`. Le prototype COD précède une annotation
de type dans un test ; sa reproduction concerne la compilation de
production. Les tests complets portent sur le binaire final. Chaque
variante doit être compilée et copiée dans son propre paquet avant de
lancer les mesures ; une extension déjà chargée ne doit pas être écrasée.

La couverture universelle, le passage à grande taille sur ces dépendances
générales et la supériorité sur MBDyn ou Simpack restent à établir. Cette
étape n'inclut aucune nouvelle exécution de ces deux logiciels.
