# Princeton strict : précision des positions

Cette page archive le diagnostic avant modification du stockage. Le
[correctif dans le noyau](POSITIONS_COMPENSEES.md) et ses résultats complets
sont décrits séparément.

Le refus du premier petit palier à quarante intervalles révèle une limite
des positions absolues stockées en double précision. Recalculer seulement
les forces avec davantage de chiffres ne lève pas cette limite. Une
expérience séparée, conservant les déplacements autour de la géométrie
initiale, réduit fortement ce défaut.

**Ce diagnostic ne corrige pas encore le solveur complet.** Il reconstruit
l'équilibre en translation à rotations fixées. Les équations de moment
ne sont pas résolues de nouveau, et le résultat n'est pas chargé dans le
noyau. Aucune tolérance ni règle d'acceptation n'est modifiée.

## Expérience reproductible

On reprend Princeton, formulation `milieu`, au premier des cinquante
paliers externes : fraction de charge `0.0009866357858642205`, soit
environ 0,0987 % de la charge finale. Les réglages sont `tol=1e-8`,
100 itérations et un seul palier interne autorisé.

Le mode habituel fournit l'état étudié. Le mode strict est rejoué
séparément : il réussit à dix et vingt intervalles, et échoue à quarante
et soixante, avec restauration vérifiée. L'état analysé n'est donc pas
la dernière tentative du calcul strict avant sa restauration.

On conserve les rotations nodales du résultat habituel, puis :

1. On recalcule les forces en double précision avec une implémentation
   NumPy séparée. L'écart maximal avec les translations du résidu du
   noyau reste inférieur à 6,5e-14 N sur les quatre maillages.
2. On recalcule aux **mêmes positions** en `longdouble` : le résidu reste
   du même ordre. Sur cette machine, ce type offre 64 bits significatifs,
   contre 53 pour `float64` ; ce n'est pas de la quadruple précision.
3. À rotations fixées, on reconstruit les positions en précision étendue
   pour que chaque poutre porte le même effort au bout. On mesure le
   résidu, puis celui obtenu après le seul arrondi des positions en f64.
4. On conserve plutôt les déplacements par rapport aux positions
   initiales, chacun en f64, et on évalue la déformation sans reformer
   d'abord les grandes positions absolues.

Norme infinie du résidu de **translation des nœuds libres**, en newtons :

| Intervalles | Noyau, état habituel | Mêmes positions, calcul étendu | Reconstruction étendue | Positions reconstruites arrondies en f64, calcul étendu | Déplacements conservés en f64 |
|---|---:|---:|---:|---:|---:|
| 10 | 3,16e-9 | 2,61e-9 | 8,32e-13 | 2,74e-9 | 4,11e-10 |
| 20 | 5,18e-9 | 5,83e-9 | 3,40e-12 | 7,61e-9 | 3,65e-10 |
| 40 | 2,10e-8 | 2,11e-8 | 4,88e-12 | 1,99e-8 | 4,73e-10 |
| 60 | 3,91e-8 | 3,88e-8 | 9,62e-12 | 2,97e-8 | 5,74e-10 |

À quarante intervalles, la reconstruction déplace les positions de moins
de 5,1e-17 m. La raideur axiale par intervalle vaut environ 2,24e8 N/m :
de si petites différences peuvent donc peser sur un critère de force à
1e-8 N. C'est un ordre de grandeur, pas une borne d'erreur certifiée.

Extraire les déplacements **après** avoir arrondi les positions ne
récupère pas les chiffres perdus : le résidu reste à 2,02e-8 N pour
quarante intervalles. Il faut conserver ces termes pendant le calcul.

## Reconstruction et périmètre

Pour la poutre `milieu`, à rotation moyenne matérielle `R` fixée :

```
F = R C (Rᵀ d/L - e1)
d = L R⁻ᵀ (e1 + C⁻¹ R⁻¹ F)
```

Les cordes `d` s'accumulent depuis l'encastrement. Le calcul utilise les
inverses de `R`, sans supposer que les matrices stockées sont parfaitement
orthogonales ni les réorthonormaliser en précision étendue. Les longueurs
et constantes matérielles restent les mêmes valeurs f64 que dans le
modèle initial.

Pour cette géométrie droite suivant x, écrire `r = r_initial + u` donne
la même déformation mathématique sous la forme :

```
gamma = (Rᵀ e1 - e1) + Rᵀ (u_B - u_A)/L
```

Le second terme peut rester petit sans disparaître dans une position de
l'ordre du demi-mètre. L'expérience établit une piste pour le stockage
des translations ; elle ne démontre ni la convergence des cinquante
paliers, ni celle de la formulation intégrée, ni un avantage sur MBDyn.

Une intégration dans le noyau doit préserver cette représentation dans
les mises à jour, les dérivées, les sauvegardes et les restaurations. Le
contrat public d'état doit permettre de reproduire le même état mécanique ;
un complément caché perdu par `etat()` / `pose_etat()` ne suffirait pas.
Les autres interactions et la dynamique devront être contrôlées aussi.

## Contrôles et artefacts

Les quatre nouveaux tests du diagnostic passent : chaîne anisotrope à
solution affine exacte, même solution tournée et translatée, moyenne de
rotations planes, et barre dont l'allongement `2^-60` disparaît dans la
position f64 `1 + 2^-60` mais reste conservable comme déplacement f64.

Le noyau et l'extension sont inchangés depuis l'étape
[statut statique](STATUT_STATIQUE.md) ; leurs empreintes sont contrôlées.
Les campagnes générales de cette étape précédente ne sont pas rejouées
pour cet ajout de diagnostic.

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/test_precision_statique.py
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_precision_statique.py \
  --sortie /tmp/precision-statique.json
```

Le programme refuse une plateforme où `longdouble` n'est pas plus précis
que f64. Il archive les positions étendues sous forme de chaînes décimales
afin de ne pas les arrondir lors de la sérialisation JSON.

[États, résidus et empreintes](bancs/precision-statique-diagnostic.json).
[Contrôles et journaux](bancs/precision-statique-bilan.json).
L'[objectif face à MBDyn](OBJECTIF_MBDYN.md) reste actif.
