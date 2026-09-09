# Globalisation de Newton : corrections et fermeture des liaisons

La version corrective 0.8.1 réduit le travail de Newton statique. La
campagne ci-dessous compare les roues isolées 0.8.0 et 0.8.1, à paramètres
identiques, avec trois répétitions par configuration.

## Problème et comportement recherché

Sur les cinquante paliers de la console anisotrope Princeton, la 0.8.0
répète une tentative en mérite de forces avant de restaurer le palier et
de réussir en mérite de corrections : 879 évaluations principales et
100 tentatives. Une faible erreur de configuration peut produire une
forte erreur de force dans une direction très raide. La norme de force
fait alors rejeter un pas qui réduit efficacement la correction de pose.

Le candidat utilise cette correction pendant la première tentative.
Il autorise un pas refusé par les forces lorsqu'il réduit suffisamment
la correction prédite **et** respecte un contrôle séparé des contraintes
géométriques. La fermeture des liaisons est indispensable : accepter une
petite correction préconditionnée seule a fait diverger deux chaînettes.

Le verdict final conserve les critères physiques de `statique` : résidu
libre et contraintes sous leur tolérance, ou stagnation explicitement
signalée dans le mode non strict. Une petite correction ne constitue
jamais à elle seule un verdict de convergence.

## Dérivation de la règle implémentée

À une itération donnée, la raideur tangente avec précontrainte et les
contraintes forment le système augmenté `A δx = b`. Les premières
composantes de `δx` sont les corrections de pose, les suivantes les
corrections de réaction. Le noyau construit déjà un équilibrage diagonal
positif `D` et résout `A_s = D⁻¹ A D⁻¹`, avec `b_s = D⁻¹ b`.

Pour un essai `α`, les rotations suivent l'exponentielle spatiale et les
translations utilisent les incréments compensés. Le résidu d'essai
emploie la réaction `λ + α δλ`, cohérente avec la direction calculée.
Reprojeter les réactions à chaque essai changerait la fonction évaluée.

Notons `P_q` l'extraction des coordonnées de pose et `D_q` leur métrique.
La norme de référence est `C₀ = ‖D_q δq‖₂`. La correction simplifiée de
l'essai est `C(α) = ‖P_q A_s⁻¹ D⁻¹ b(α)‖₂` : elle réutilise la
factorisation de l'itération, sans former une nouvelle tangente. Le modèle
linéaire prédit une correction résiduelle `(1−α) δq`.

La nouvelle acceptation exige simultanément :

1. `C₀` strictement positif et fini, `C(α)` fini ;
2. `C(α) ≤ (1−α/4) C₀` ;
3. `‖Φ(q_α)‖∞ ≤ max(tol, (1−α/4) ‖Φ(q)‖∞)`.

La deuxième condition exige un quart de la diminution prédite par le
modèle linéaire. La troisième utilise les contraintes originales, avant
réduction des dépendances, dans les unités et à la tolérance de l'API.
Elle n'interdit pas les pas déjà acceptés par le mérite de forces : elle
encadre l'acceptation supplémentaire par les corrections.

La correction supplémentaire n'est calculée que si les forces refusent
l'essai et si le filtre géométrique l'autorise. Elle nécessite une
résolution directe ou raffinée contrôlée ; les pas obtenus par le repli
orthogonal conservent la stratégie existante. Le pas effectivement retenu
est celui qui satisfait les deux contrôles, même si un essai précédent
avait une norme de force plus petite.

La fenêtre d'abandon fondée sur dix valeurs de force ne condamne pas un
progrès accepté par les corrections dans cette fenêtre. Le budget
d'itérations reste borné. Après un échec, la restauration du palier et la
tentative historique en mérite naturel restent disponibles, puis la
subdivision en charge. Le mode sans amortissement conserve son chemin.

## Sources, transfert et limites théoriques

Les méthodes Newton orientées vers l'erreur et l'invariance affine sont
un socle exposé par Deuflhard dans *Newton Methods for Nonlinear Problems*
(2004). La notice de l'éditeur a été consultée ; nous ne revendiquons pas
une lecture intégrale du livre ni la reproduction d'un de ses algorithmes.
[Springer, ouvrage de 2004](https://link.springer.com/book/9783540210993).

L'étude d'Aparicio-Estrems, Gargallo-Peiró et Roca (2024) associe recherche
de pas, préconditionnement et contrôle du travail pour une optimisation
de maillage. Elle motive le contrôle des évaluations et de la prédiction
locale. Son optimisation et son Newton–CG ne constituent pas notre système
multicorps augmenté, généralement indéfini.
[Publication et texte intégral](https://arxiv.org/html/2403.13654v1).

**La règle ci-dessus est une adaptation propre, avec le coefficient 1/4
choisi pour cette expérience.** Ce coefficient n'est pas déduit de la
constante homonyme de l'article de 2024. Le filtre n'est pas une
implémentation complète d'une méthode de filtre d'optimisation publiée.
Nous ne revendiquons ni convergence globale sur tout mécanisme, ni
invariance affine de tout le solveur. Les changements de métrique,
contraintes redondantes, bifurcations et arrondis restent déterminants.

## Expériences discriminantes

Les variantes précédentes ont été conservées dans des roues distinctes :

- Une acceptation en corrections sans contrôle géométrique améliore
  Princeton, mais fait échouer les chaînettes à 80 et 120 segments.
- Vérifier la contraction avec la tangente de l'itération suivante ne
  rétablit pas les deux chaînettes. Avec le filtre géométrique, ce contrôle
  ajoute aussi des reprises inutiles sur les poutres et boucles tournées.
- Le candidat courant conserve le filtre et retire cette confirmation.

Les nouveaux tests vérifient la petite charge anisotrope sans reprise
répétée, avec contrôle indépendant de la résultante et du moment à
l'encastrement. Pour les chaînettes, une référence **discrète** est calculée
par un problème scalaire indépendant : la tension horizontale est telle
que la somme des projections horizontales des segments vaut la portée.
Les positions, longueurs, tensions et bilans nodaux sont contrôlés.
Les contrôles négatifs détectent la reprise dans la roue 0.8.0 et la
divergence dans la variante sans filtre.

## Campagne et résultats

`ci/diagnostic_globalisation.py` construit les modèles et leurs contrôles.
`ci/mesure_globalisation.py` compare deux binaires identifiés par empreinte,
en alternant leur ordre et en conservant les refus. Le corpus comprend
51 configurations : poutres de deux formulations, rampes de cinquante
paliers, directions de charge, changements d'unités, rotations du modèle,
moment pur, chaînettes, boucles en cascade, contacts Hertz et absence
d'équilibre. Les bilans sont indépendants lorsque les références du cas
le permettent ; les états finaux sont aussi comparés entre binaires.

L'[archive des 306 résolutions chronométrées](bancs/globalisation-statique-0.8.1.json)
et son [bilan recalculable](bancs/globalisation-statique-0.8.1-bilan.json)
conservent trois répétitions pour chaque version et chacune des
51 configurations. Chaque processus effectue un échauffement identique
avant sa résolution mesurée ; le résultat détaillé de cet échauffement
n'est pas conservé. Les nombres d'évaluations sont identiques sur les
trois répétitions de chaque couple version–configuration :

| Cas | 0.8.0 | Candidat |
|---|---:|---:|
| Princeton, cinquante paliers, deux formulations | 879 évaluations / 100 tentatives | 334 / 50 |
| Princeton intégré, charge directe, 8 éléments | 45 évaluations | 34 |
| Six chaînettes, 40 à 120 segments | 8 évaluations chacune | 8 chacune |
| Cascade tournée, 32 boucles | 10 évaluations | 8 |

Temps médians en millisecondes, AMD EPYC 7543, processeur logique 0,
un fil demandé par bibliothèque ; Python 3.14.7, NumPy 2.5.3 et SciPy 1.18.1.
Les deux moteurs sont exécutés successivement et alternent leur ordre.
Aucune compilation ni campagne de tests n'accompagne ces mesures.

| Configuration | 0.8.0 (ms) | 0.8.1 (ms) | Temps 0.8.1 / 0.8.0 |
|---|---:|---:|---:|
| Princeton milieu, 8 éléments, rampe | 250.181 | 76.397 | 0.305 |
| Princeton milieu, 20 éléments, rampe | 723.631 | 227.913 | 0.315 |
| Princeton milieu, 60 éléments, rampe | 1971.134 | 611.319 | 0.310 |
| Princeton intégré, 8 éléments, rampe | 283.283 | 85.977 | 0.304 |
| Princeton intégré, 20 éléments, rampe | 809.409 | 253.772 | 0.314 |
| Princeton intégré, 60 éléments, rampe | 2243.155 | 682.134 | 0.304 |
| Chaînette, 80 segments | 7.394 | 7.472 | 1.011 |
| Chaînette, 120 segments | 11.086 | 11.173 | 1.008 |
| Cascade, 32 boucles | 35.416 | 30.355 | 0.857 |
| Cascade tournée, 32 boucles | 556.122 | 485.982 | 0.874 |

Les rampes sont **3.17 à 3.29 fois plus rapides** dans ce protocole interne.
Les six chaînettes varient entre −0.5 % et +1.8 % de temps ; aucune avance
de vitesse n'y est revendiquée. Les deux contacts Hertz varient de moins
de 1 %. Les petits écarts restent sensibles au bruit de chronométrage.

Les 49 configurations satisfaisant leurs contrôles dans la roue 0.8.0
les satisfont aussi dans la roue 0.8.1, sur les trois répétitions. Ce total
comprend le contrôle négatif « absence d'équilibre », dont le résultat
attendu est un refus avec restauration : il ne signifie donc pas
49 équilibres calculés. Sur les cas convergés, l'écart maximal de position
entre versions est `5.801e-10 m` et celui des composantes des matrices de
rotation `1.997e-9`. Ces écarts entre binaires ne remplacent pas les
références physiques indépendantes du corpus.

Deux poutres tournées, avec unité de longueur multipliée par `10⁻³`,
échouent avec **les deux versions** au réglage du diagnostic. La tolérance
scalaire y est resserrée à `10⁻¹¹`, à la fois pour les forces et pour les
contraintes ; ce n'est pas une comparaison à tolérance physique identique
dans toutes les composantes. La normalisation et le plancher numérique
méritent un chantier distinct. La 0.8.1 consacre davantage d'itérations à
ces échecs, sa fenêtre de forces étant moins restrictive : 197 ou
200 évaluations contre 111. Leur temps de refus passe respectivement de
48.810 à 78.599 ms et de 55.874 à 89.911 ms, soit **environ 61 % de plus**.
Les restaurations passent. Ces cas restent inclus et ne sont pas comptés
comme des équilibres acquis.

Le chronomètre comprend l'application des charges, `statique`, la collecte
des rapports et des sauvegardes cinématiques. Imports, construction du
modèle, échauffement et contrôles indépendants sont hors chronomètre.
Un seul effort est mis à jour pendant les rampes. Ce coût interne diffère
du protocole d'usage complet de la confrontation MBDyn. **Il ne permet
pas de recalculer par division le rapport de vitesse face à MBDyn ou à
Simpack.**

La [confrontation externe de la 0.8.1](CONFRONTATION_STATIQUE_MBDYN_0.8.1.md)
a depuis été exécutée séparément : 150 calculs et 3 800 paliers Vinkulum
stricts passent. Au seuil de 10 µm au bout, MBDyn conserve un avantage
de temps de 44.6 fois face à la formulation milieu et 11.4 fois face à
l'intégrée. Les coûts totaux, marges des références et mémoires y sont
présentés avec leurs limites.

## Vérification du lot

Les contrôles couvrent 70 tests Rust, 8 du prototype flexible, 78 tests
Python, 41 groupes de vérification, 46 bancs mécaniques et 9 bancs de
contact. La roue finale installée hors du dépôt passe les 78 tests Python,
les 41 vérifications et l'exemple du README avec restauration précise.
Les étapes et empreintes sont précisées dans le
[relevé de livraison](bancs/version-0.8.1.json).

Le vérificateur de campagne garde explicitement les deux limites connues :

```bash
python ci/bilan_globalisation.py docs/bancs/globalisation-statique-0.8.1.json
```

Il relit les résultats, contrôle les répétitions, les statuts, les
restaurations et les écarts d'état, puis recalcule le bilan. Sa réussite
atteste la cohérence de cette archive et l'absence de nouvelle perte de
convergence sur ce corpus ; elle ne démontre pas une robustesse générale.
