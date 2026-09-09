# Newton statique singulier : équilibrage et raffinement creux

La 0.7.1 traite davantage de systèmes statiques redondants sans densifier
la matrice de Newton. Elle associe un équilibrage par blocs de coordonnées
spatiales à une factorisation auxiliaire creuse, dont les corrections sont
jugées sur le système original. Le modèle mécanique et le contrat public
de `statique()` sont conservés.

Les [mesures alternées](bancs/raffinement-mesures.json) comparent la roue
0.7.0 archivée au nouveau noyau, sur 26 configurations et trois répétitions
par version. Le [bilan](bancs/raffinement-bilan.json) publie aussi les coûts
défavorables, les plages de temps, la mémoire et les contrôles physiques.
Cette comparaison porte sur deux versions de Vinkulum ; elle ne constitue
pas une nouvelle confrontation exécutée avec MBDyn ou Simpack.

| Cas | 0.7.0 | 0.7.1 | Rapport des temps 0.7.0 / 0.7.1 |
|---|---:|---:|---:|
| Cascade, 8 cellules | 70,245 ms | 5,017 ms | 14,00× |
| Cascade, 16 cellules | 250,198 ms | 12,544 ms | 19,95× |
| Cascade, 32 cellules / 96 corps | 2,0384 s | 37,652 ms | 54,14× |
| Même cascade tournée, 96 corps | 2,1339 s | 639,634 ms | 3,34× |
| Parallélogrammes indépendants, 32 cellules | 802,278 ms | 10,074 ms | 79,64× |
| Princeton intégré, 40 intervalles, encastrement double | 81,621 ms | 76,843 ms | 1,062× |
| Chaîne articulée, 128 corps | 10,073 ms | 10,755 ms | 0,937× |
| Princeton milieu, 10 intervalles, encastrement double | 16,418 ms | 17,471 ms | 0,940× |
| Rotors indépendants, 32 corps | 0,559 ms | 0,595 ms | 0,940× |

Ce sont les médianes de trois essais, avec une tolérance statique stricte
de `1e-8`. La cascade de 96 corps conserve neuf évaluations de Newton et
son pic mémoire passe de **86,0 à 40,9 Mio**. Le cas tourné utilise dix
évaluations et son pic augmente de **107,3 à 108,8 Mio**. Le pire surcoût
des 26 configurations atteint **6,8 %**, sur la chaîne articulée à 128 corps.
Les faibles écarts sur les cas sous la milliseconde ne sont pas présentés
comme une amélioration ou une régression universelle ; les plages de
mesure restent consultables dans le bilan.

## Défaut observé

La cascade connexe à 32 cellules contient 96 corps, 576 coordonnées,
640 contraintes et 96 redondances. Son système augmenté a une dimension
de 1 216. L'ancien équilibrage prenait la racine de la diagonale physique,
avec une valeur de secours seulement sous `1e-30`. Une diagonale presque
nulle produisait donc une échelle de `2,24e-9`, bien que les termes de liaison
soient ordinaires.

Le [diagnostic indépendant avec LAPACK](bancs/raffinement-diagnostic.json)
rejoue la première matrice effectivement assemblée :

| Matrice | Plus grand coefficient | Norme spectrale | Rang calculé au seuil absolu `1e-12` |
|---|---:|---:|---:|
| Avant équilibrage, reconstruite | 4,415 | 12,900 | 1 120 |
| Ancien Jacobi | 700 291 923 | 1 337 483 275 | 1 156 |
| Nouveaux blocs | 1,104 | 3,139 | 1 120 |

Le rang calculé après Jacobi illustre la sensibilité d'un seuil absolu à
l'amplification numérique. Il ne prouve pas l'apparition de directions
physiques supplémentaires. Le spectre complet est conservé ; aucun de ces
seuils n'est présenté comme un certificat de rang des données.

La matrice brute n'est pas symétrique : `max|A−Aᵀ| = 0,05157`. Les termes
de Newton en rotations hors équilibre imposent de garder ce cas dans le
solveur général. [MINRES-QLP](https://web.stanford.edu/group/SOL/software/minresqlp/)
exige une matrice symétrique ou hermitienne et ne peut donc pas être choisi
ici sans transformation supplémentaire.

## Algorithme retenu

Pour `A δx = b`, le noyau construit `Aₛ = D⁻¹ A D⁻¹`, `bₛ = D⁻¹ b`,
puis restitue `δx = D⁻¹ y`. Chaque translation et chaque rotation partagent
une échelle pour leurs trois composantes ; les multiplicateurs sont
équilibrés individuellement. Les normes de Frobenius des blocs de lignes
et de colonnes sont invariantes sous rotation orthogonale de ces vecteurs.
Cela évite de fonder l'échelle sur une seule composante diagonale.

La démarche s'inspire de l'équilibrage simultané de
[Ruiz](https://cerfacs.fr/wp-content/uploads/2017/06/14_DanielRuiz.pdf).
La variante implémentée arrondit les facteurs à des puissances de deux et
borne le travail à huit passes. Elle réutilise les échelles du Newton
précédent ; si celles-ci perdraient un coefficient du nouveau système,
elle repart de l'identité. Elle ne promet ni conditionnement optimal ni
invariance bit à bit près des frontières de quantification.

Les changements d'exposant sont contrôlés par reconstruction exacte des
coefficients. Le second membre et la correction rendue sont également
contrôlés. Une perte par sous-flux ou une valeur non représentable provoque
un refus explicite, avec la restauration prévue par l'appel statique.

Le solveur direct reste le premier essai. Quand il échoue sur un grand
système avec contraintes, le noyau factorise

```text
B = Aₛ − μ diag(0_coordonnées, I_multiplicateurs)
μ = √ε × max|Aₛᵢⱼ|
y₀ = B⁻¹ bₛ
rⱼ = bₛ − Aₛ yⱼ
yⱼ₊₁ = yⱼ + B⁻¹ rⱼ
```

Le décalage appartient uniquement à la factorisation auxiliaire. Les
résidus portent toujours sur `Aₛ`, sans ajouter de raideur aux coordonnées
physiques. Le motif de la factorisation symbolique est réutilisé seulement
s'il est identique. Le refus du LU creux ne déclenche plus un LU dense
supplémentaire avant cet essai.

Les produits et sommes de `bₛ−Aₛ y` conservent leurs petits termes avec
`TwoSum` et le reste du produit par FMA. Ces transformations s'appuient sur
les travaux de [Ogita, Rump et Oishi](https://www.tuhh.de/ti3/paper/rump/OgRuOi05.pdf).
La tentative est acceptée si `‖r‖₂ ≤ 1e-12 ‖bₛ‖₂`, avec au plus huit
contrôles ; elle est abandonnée si le résidu ne baisse pas d'au moins 25 %,
devient non fini ou dépasse ce budget. Le repli orthogonal dense reste
disponible. La même factorisation sert à la recherche de pas naturelle.

Les travaux récents d'[Epperly, Greenbaum et Nakatsukasa (2025)](https://arxiv.org/abs/2502.17767)
montrent le rôle du raffinement et du redémarrage dans la stabilité de
méthodes préconditionnées pour systèmes non symétriques. Ils motivent le
contrôle du résidu recalculé, mais leurs résultats ne sont pas transposés
en garantie générale pour les systèmes singuliers de Vinkulum. Ici, le
décalage est un choix heuristique contrôlé à l'exécution.

## Essais écartés et coûts de préparation

[LSMR](https://stanford.edu/group/SOL/software/lsmr/LSMR-SISC-2011.pdf)
sans préconditionneur adapté atteint sa limite de 9 728 itérations sur la
matrice brute comme sur celle équilibrée par Jacobi. Le diagnostic archive
les résidus recalculés. La seule suppression du mauvais Jacobi ne suffit
donc pas à fournir un solveur itératif rapide.

Un décalage auxiliaire de `1e-3` diverge sur la première matrice équilibrée
par blocs : le résidu passe de 0,260 à 6,147 dès la correction suivante.
Les essais `1e-5`, `1e-7` et `1e-9` atteignent le seuil. Cette expérience
justifie le rejet sur croissance et le contrôle sur la matrice originale.
Le diagnostic SciPy sert à l'algèbre ; les temps annoncés viennent des
calculs mécaniques du noyau Rust.

Le [premier prototype](bancs/raffinement-prototype-bilan.json) accélérait
déjà les cascades, mais ajoutait jusqu'à 25 % au temps de certains témoins.
Le code final construit directement les puissances de deux, réutilise les
inverses dans l'accumulation des normes et conserve les échelles entre
itérations. Les [mesures du prototype](bancs/raffinement-prototype-mesures.json)
et son [patch contre le commit 5d0559c](bancs/raffinement-prototype.patch)
sont conservés ; les deux candidats ne sont pas confondus.

## Validation et reproduction

Les cascades sont contrôlées contre l'équilibre dérivé de leur énergie
réduite, avec bilans indépendants des forces, moments et fermetures. La
même mécanique est aussi tournée par le vecteur de rotation
`[0,31, −0,47, 0,22]`. Les chaînes et poutres vérifient leurs torseurs de
réaction, notamment le partage entre liaisons dupliquées. Les configurations
déjà immobiles sont identifiées par leur nombre d'évaluations dans le bilan.

La CI étendue couvre 57 tests Rust, 67 tests Python, 41 vérifications,
46 bancs mécaniques et 9 bancs de contact. Les 400 paliers Princeton
restent stricts, à 10, 20, 40 et 60 intervalles pour les deux formulations.
Trois barres conservent exactement l'allongement de `2⁻⁶⁰ m` dans les deux
composantes exportées de la position. Le
[relevé de validation](bancs/raffinement-validation.json) lie ces résultats
aux sources et aux empreintes des extensions.

```bash
# AVANT contient le paquet extrait de la roue 0.7.0 ; CPU est disponible.
python ci/mesure_raffinement.py --version avant="$AVANT" \
  --version apres="$PWD/python" --cpu "$CPU" --repetitions 3 \
  --sortie docs/bancs/raffinement-mesures.json
python ci/bilan_raffinement.py docs/bancs/raffinement-mesures.json \
  docs/bancs/raffinement-bilan.json
OPENBLAS_NUM_THREADS=1 python ci/diagnostic_raffinement.py \
  --jacobi docs/bancs/raffinement-matrices-jacobi.tar.gz \
  --blocs docs/bancs/raffinement-matrices-blocs.tar.gz \
  --sortie docs/bancs/raffinement-diagnostic.json
python ci/valide_raffinement.py
```

La mesure utilise des processus neufs alternés, le même CPU et les mêmes
versions Python, NumPy et SciPy, avec un seul fil demandé. Le temps couvre
l'appel statique ; la mémoire est le pic du processus. Les compilations et
autres campagnes sont terminées avant les mesures de performance. Les
imports propres au repère tourné interdisent de comparer directement son
pic mémoire à celui du cas plan ; la comparaison entre versions reste
faite au sein du même cas.

Les archives contiennent les triplets, seconds membres et échelles des
itérations réelles. Le [patch d'export](bancs/raffinement-export.patch)
décrit l'instrumentation temporaire du noyau 0.7.0 ; elle ne fait pas partie
de l'extension distribuée. L'archive par blocs provient du prototype
d'équilibrage, avant le raffinement creux.

## Limites et suite

Le résidu linéaire accepté ne constitue pas une borne générale de l'erreur
de déplacement. La correction intermédiaire des multiplicateurs n'est pas
promise de norme minimale dans toutes les métriques. Les réactions finales
continuent d'être calculées par le projecteur QR/COD existant ; ses
composantes réellement déficientes restent denses, notamment dans le cas
tourné. Les modes physiques libres peuvent laisser la factorisation
auxiliaire singulière ; les contraintes contradictoires conservent le repli.

Le calcul public `k_c_m_z()` échoue encore sur la cascade initiale à
32 cellules (`SymbolicSingular`, indice 58). Cette limite des matrices
d'analyse est distincte du calcul statique corrigé ici. Les prochains
travaux doivent traiter les projecteurs creux de rang déficient et ce
chemin d'analyse, puis rejouer les confrontations externes à précision
commune. L'objectif généraliste reste ouvert.
