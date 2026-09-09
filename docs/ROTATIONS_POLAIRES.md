# Rotations polaires et SVD normalisée

Cette étape traite les remises sur SO(3) en assemblage, en recherche de pas
statique et dans l'intégrateur énergie–moment. Elle corrige aussi une erreur
de SVD sur des matrices denses très petites ou très grandes. Les comparaisons
portent sur des versions de Vinkulum ; aucune nouvelle exécution MBDyn ou
Simpack ne permet encore de conclure à une supériorité externe.

## Projection près de SO(3)

Pour une matrice `A` à déterminant positif, on pose `E = AᵀA − I`.
Lorsque `||E||F ≤ 1/8`, la projection utilise l'itération de Newton–Schulz :

```text
X₀ = A
Eₖ = XₖᵀXₖ − I
Xₖ₊₁ = Xₖ − ½ XₖEₖ.
```

L'identité `Eₖ₊₁ = −¾Eₖ² + ¼Eₖ³` donne
`||Eₖ₊₁||F ≤ (¾ + ¼||Eₖ||F)||Eₖ||F²`. À partir de `1/8`, quatre
itérations suffisent en arithmétique exacte pour passer sous epsilon.
Les valeurs singulières initiales sont alors comprises entre `√(7/8)` et
`√(9/8)` ; les facteurs de correction restent positifs et conservent
l'orientation. La limite est le facteur polaire propre.

Le noyau contrôle le résidu **recalculé en f64**, avec arrêt à `8ε` et
au plus cinq corrections. Il applique au moins une correction, même si
l'entrée est déjà proche du groupe. Les comparaisons utilisent les normes
au carré pour éviter les racines carrées. Le chemin proche n'alloue aucune
matrice sur le tas : ses produits sont des matrices 3 × 3 de taille fixe.

Cette utilisation locale s'appuie sur les itérations polaires étudiées par
[Higham et Schreiber](https://eprints.maths.manchester.ac.uk/340/).
Leur intérêt est ici de remplacer une factorisation générale répétée par
des produits très petits. La stabilité de Newton–Schulz est conditionnelle,
comme le montrent
[Nakatsukasa et Higham](https://eprints.maths.manchester.ac.uk/1839/).
La borne locale et le résidu calculé ne constituent pas un certificat par
intervalles, ni une garantie globale d'erreur rétrograde.

Pour des changements de repères orthogonaux propres `L` et `B`, l'itération
commute, en arithmétique exacte, avec `A ↦ LAB`. Elle ne choisit aucune colonne
privilégiée. Les tests fabriquent `A = Q P D Pᵀ` avec `D` positive : le
facteur connu est `Q`. Ils vérifient cette référence, les changements de
repères, les zéros d'une rotation plane et 40 000 projections composées
dans les deux sens.

## Repli et changement d'échelle

Hors du voisinage, ou après un résidu insuffisant, le noyau calcule une SVD.
Pour `A = UΣVᵀ`, il forme `U diag(1, 1, det(UVᵀ)) Vᵀ` : la correction
d'une réflexion porte sur la plus petite valeur singulière, selon l'ordre
décroissant de faer. Le résultat est une rotation propre la plus proche au
sens de Frobenius. Elle peut être non unique pour une matrice déficiente
en rang ou un spectre dégénéré ; le test nul vérifie son appartenance à
SO(3), sans lui imposer une orientation particulière.

Ce contrôle a révélé un défaut que les anciens tests diagonaux masquaient.
Sur `A = Q diag(2s, s, −s/2) Pᵀ`, avec deux rotations denses indépendantes,
la SVD faer brute retrouve, pour `s = 10⁻²⁰⁰`, les valeurs normalisées
`[1,526823528 ; 0,993094333 ; 0,659508888]` au lieu de `[2 ; 1 ; 0,5]`.
L'erreur du facteur propre atteint **1,513** en norme de Frobenius.

| Échelle uniforme `s` | faer brut | Adaptateur normalisé |
|---|---|---|
| `10⁻³⁰⁰`, `10⁻²⁰⁰` | Spectre et rotation erronés | Spectre et reconstruction validés |
| `1` | Validé | Validé |
| `10²⁰⁰`, `10³⁰⁰` | `NoConvergence` | Spectre et reconstruction validés |

`svd_sure` divise maintenant l'entrée par une puissance de deux tirée du
plus grand coefficient absolu, puis remet les valeurs singulières à
l'échelle. Un aller-retour coefficient par coefficient vérifie que cette
normalisation n'a perdu aucune donnée représentée. Il couvre aussi les
coefficients sous-normaux, dont le plus petit f64 positif.

Une très grande **étendue interne** reste une limite : par exemple,
`diag(MAX_f64, min_subnormal_f64)` est explicitement refusée. Une perte par
arrondi vers un sous-normal est également refusée, même sans passage à
zéro. Les valeurs singulières non représentables entraînent un refus. Le
correctif n'offre donc pas une SVD générale de toutes les matrices à
coefficients finis et ne remplace pas un traitement des matrices très
inégalement mises à l'échelle.

Le [diagnostic Rust](../examples/diagnostic_echelle_svd.rs) exporte les
matrices et tous les facteurs, y compris les échecs bruts. Le
[contrôle indépendant](../ci/valide_echelle_svd.py) utilise les rotations de
SciPy pour la référence fabriquée et LAPACK `gesvd` sur l'entrée originale.
Les erreurs sont évaluées après remise à l'échelle pour ne pas faire
déborder les normes. Les données sont conservées dans
[le bilan numérique](bancs/rotations-echelles-lapack.json).

## Dynamique et contrôles physiques

Le domaine déjà admis par `simule_em` reste celui des corps libres sous la
seule gravité. La projection polaire remplace son orthogonalisation de
Gram–Schmidt. La mise à jour principale en alpha généralisé conserve son
produit d'exponentielles ; cette étape n'en change pas la projection.

Le schéma au point milieu et la transformation de Cayley restent ceux
du noyau. Avec `Π̄ = (Π₀ + Π₁)/2`, l'équation
`Π₁ − Π₀ = h Π̄ × J⁻¹Π̄` conserve l'énergie quadratique en arithmétique
exacte. La mise à jour cohérente de l'orientation conserve alors le moment
spatial. Le Newton, les produits et la projection flottants introduisent
une erreur d'arrondi ; les commentaires du noyau précisent désormais
cette distinction.

Le [banc indépendant](../ci/diagnostic_rotation.py) part d'une rotation
dense et d'une toupie anisotrope proche de son axe intermédiaire :
`J = diag(0,001 ; 0,002 ; 0,003) kg·m²`,
`ωmat = (0,05 ; 20 ; 0,03) rad/s`. Il compare quatre repères matériels et
spatiaux, avec `h = 0,004 ; 0,002 ; 0,001 s` jusqu'à une seconde, puis
`h = 0,001 s` jusqu'à vingt secondes. La référence temporelle à une seconde
intègre séparément les équations d'Euler et un quaternion par DOP853.

Le validateur recalcule énergie, moment et orthogonalité sur au moins
200 instantanés par trajectoire. Ce contrôle échantillonné ne couvre pas
chaque pas exportable. Une régression supplémentaire emploie trois corps
d'inerties distinctes, des rotations matérielles denses, une chute sous
gravité et deux appels successifs. Elle confronte les translations à la
solution quadratique et vérifie l'ordre deux face aux rotations analytiques
autour des axes principaux.

## Mesures et limites

Le témoin est le commit `f445fc3`. Les processus sont neufs, alternés,
sur le même CPU avec un seul fil demandé. Chaque cas est répété trois
fois, sans compilation ni campagne concurrente pendant les chronométrages.
Les médianes ci-dessous incluent la résolution et ses contrôles internes ;
la construction des modèles et les références indépendantes sont hors du
chronométrage.

| Statique stricte | Taille | Avant, ms | Après, ms |
|---|---:|---:|---:|
| Princeton milieu, encastrement double | 10 intervalles | 19,055 | 16,486 |
| Princeton intégrée, encastrement double | 10 intervalles | 20,427 | 17,625 |
| Chaîne articulée | 128 corps | 11,062 | 10,396 |
| Rotors | 128 corps | 2,140 | 1,964 |
| Chaîne à liaisons doubles | 64 corps | 7,162 | 6,837 |
| Boucle spatiale | 128 corps | 2,500 | 2,585 |
| Chaîne soudée | 128 corps | 0,910 | 0,935 |
| Boucle soudée | 512 corps | 4,432 | 4,494 |
| Cascade connexe | 96 corps | 2 041,855 | 2 040,020 |

La baisse sur Princeton vaut **13,5–13,7 %**, avec le même nombre
d'évaluations et des écarts de position inférieurs à `4,8 × 10⁻¹⁷ m`.
Les réactions diffèrent au plus de `1,05 × 10⁻¹⁰` dans leurs unités SI.
Les grandes cascades restent dominées par la factorisation dense : le
changement de temps y est inférieur à la variation entre répétitions.
Les surcoûts observés atteignent **3,4 %** sur la boucle spatiale.
Le [bilan statique complet](bancs/rotations-statique-bilan.json) conserve
les 18 cas, les étendues, les pics RSS et les contrôles mécaniques.
Les deux campagnes statiques totalisent **216 équilibres stricts**,
dont 108 dans la comparaison du binaire final.

Sur la toupie à vingt secondes, le repère initial passe de **5,012 à
5,109 ms**, soit **1,95 % de surcoût**. Sur les quatre repères, ce surcoût
est de **1,9–2,3 %**. La première projection avec racines carrées demandait
5,233 ms sur le même repère. Les temps et les autres pas sont dans le
[bilan dynamique](bancs/rotations-dynamique-carres-bilan.json).

La dérive relative maximale du moment sur les instantanés passe de
`1,022 × 10⁻¹³` à `2,685 × 10⁻¹⁴`. Celle de l'énergie reste sous
`2,360 × 10⁻¹⁴`. Le défaut d'orthogonalité reste sous `8,1 × 10⁻¹⁶` dans
la nouvelle version. Les ordres observés sont compris entre **1,9984 et
1,9998**, par rapport à la référence indépendante à une seconde. La
précision temporelle reste celle du schéma au point milieu.

**L'équivalence numérique entre repères à long terme reste imparfaite.**
À vingt secondes, les écarts atteignent `1,722 × 10⁻⁷` en norme de
Frobenius sur la rotation et `7,243 × 10⁻⁷ rad/s` sur la vitesse angulaire,
dans les deux versions. La nouvelle projection ne résout pas cet écart
de trajectoire. L'origine et son amplification restent à isoler ; la
conservation des invariants ne suffit pas à prouver la précision de phase.

Un prototype supplémentaire stockait l'inverse de l'inertie matérielle
avec le moment, pour éviter de l'inverser à chaque pas. Malgré la réduction
du nombre d'inversions, il ralentit ce banc de **23–24 %** : 6,225 ms contre
5,034 ms dans le repère initial de cette campagne. Il a été retiré.
La cause du surcoût au niveau du code machine n'a pas été isolée. Le
[bilan du prototype](bancs/rotations-dynamique-inertie-bilan.json) conserve
aussi ses contrôles physiques réussis et les mesures des autres versions.

Les trois campagnes dynamiques totalisent **432 simulations**, avec
les trajectoires sauvegardées, les références et les résultats de chaque
répétition. Les patchs de reproduction des deux prototypes s'appliquent
à cette version avec `git apply --unidiff-zero`, dans une copie de travail
distincte ; le validateur vérifie les empreintes des sources reconstruites.

La validation finale comprend **52 tests Rust, 66 tests Python, 41 cas
généraux, 46 bancs mécaniques, 9 cas de contact et 400 paliers Princeton
stricts**. Le câble condensé extrême conserve son échec connu. Les
journaux et empreintes sont liés par le
[manifeste](bancs/rotations-validation.json). Les coûts mesurés pendant la
CI et la campagne des 400 paliers ne servent pas à comparer les vitesses.

## Reproduction

```bash
cargo run --release --example diagnostic_echelle_svd > /tmp/echelles.json
python ci/valide_echelle_svd.py /tmp/echelles.json /tmp/echelles-lapack.json
python ci/mesure_rotations.py --version avant=/tmp/avant --version apres=/tmp/apres \
  --sortie /tmp/rotations.json
python ci/bilan_rotations.py /tmp/rotations.json /tmp/rotations-bilan.json
ci/local.sh --bancs
python ci/valide_rotations.py
```

Les répertoires fournis comme `PYTHONPATH` doivent contenir chacun un paquet
`vinkulum` avec son propre fichier d'extension copié, et les modules Python
du dépôt. Il faut terminer les compilations avant les mesures isolées.
