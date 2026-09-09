# Krylov contraint : réduction rapide au-delà du premier mode intérieur

Cette note conserve la campagne initiale de 60 essais. Le
[certificat par inertie dirigée](INERTIE_CONTRAINTE_PROTOTYPE.md) apporte
ensuite une campagne de 84 essais : à 512 poutres, les champs passent de
0,937 à 0,451 s et le service avec contrôle de 1,305 à 0,816 s dans cette
nouvelle confrontation. Les nombres ci-dessous restent ceux du premier lot.

Travail expérimental du 8 septembre 2026, après la livraison 0.11.0.
Le prototype conserve une direction intérieure et réduit son complément,
au lieu de résoudre un grand KKT à chaque fréquence. Dans la campagne de
**60 essais**, les **24 essais Krylov** passent le seuil physique 10⁻⁶
sur 0–40 Hz ; les 12 essais avec contrôle passent aussi leurs majorants.
Sur 512 poutres, les champs seuls sont 7,3 fois plus rapides pendant la
phase de réponse que la LU corrigée, mais le coût total reste 29 % plus
élevé. Les sources sont dans `ci/` ; ce travail ne modifie pas encore
l'API publique.

## Le mécanisme mathématique

La réduction matérielle publique refuse une bande dès que son minorant
spectral intérieur ne garantit plus la coercivité. Le
[certificat du complément](COMPLEMENT_SPECTRAL_DIRIGE.md) contourne cette
limite en gardant une direction intérieure parmi les inconnues. Son
covecteur effectivement stocké, `B = fl(M_ii Phi)`, définit le complément
`ker(Bᵀ)`. Le certificat vise exactement ce B, même si la direction est
approchée. Les vrais pôles du système complet restent présents.

La nouvelle construction applique l'inverse statique **contrainte** à
tous les couplages du relèvement conservé. Elle réutilise un facteur QR
et des réflecteurs de Householder ; elle ne forme ni l'inverse complète,
ni une base dense du noyau des contraintes. Quatre blocs de Krylov
construisent ensuite une base fixe. Les produits énergétiques emploient
le facteur physique D d'entrée, et les produits massiques M d'entrée.

Deux services sont mesurés séparément :

- **Champs** : un petit système projeté rend les six champs physiques et
  leurs normes à chaque fréquence.
- **Champs avec contrôle** : une enveloppe résiduelle préparée sur la
  bande, puis une marge sur le bloc conservé entier, donnent des majorants
  par charge en masse et en déformation. Ce service reconstruit aussi les
  champs et évalue leurs normes physiques.

Les [preuves de l'enveloppe](KRYLOV_CONTRAINT_PREUVES.md) et de sa
[variante par direction](KRYLOV_CONTRAINT_ANISOTROPIE.md) détaillent la
composition. En particulier, si δ majore la norme du résidu et η(q) son
action sur la réponse conservée, le défaut de Schur vérifie

\[
\|(\bar S-S)q\|\le
\frac{\delta\,\eta(q)}{\lambda_c-\Omega^2}.
\]

Remplacer le numérateur par η(q)² serait faux. Cette distinction permet
un contrôle plus précis pour les directions peu affectées par le résidu,
sans ignorer son amplification par le système conservé.

## Ce qui est prouvé et ce qui est évalué

La minoration spectrale du complément utilise des arrondis dirigés. La
séparation `lambda_c > omega_max²` est aussi comparée exactement sur les
nombres binary64 stockés. Le pilote choisit **16 chiffres décimaux avant
la campagne**, après des sondes à 16, 24, 32 et 80 chiffres ; le défaut
du certificat entraîne un refus, sans reprise à une autre précision.
La valeur par défaut de la fonction de certification reste 80 chiffres.

Les identités du contrôleur sont démontrées en arithmétique exacte. Les
normes, résidus, rangs, réparations et petits systèmes qui les évaluent
restent calculés en binary64 : **les réponses ne sont pas certifiées
machine**. Une enveloppe résiduelle uniforme n'est pas non plus une
certification de réponse entre les fréquences testées, car la marge du
système conservé est évaluée à chaque fréquence.

La construction répare les contraintes après normalisation des nouvelles
directions. Un test reproduit l'amplification du défaut que provoquerait
l'omission de cette dernière réparation. Le contrôle conserve ensuite
explicitement son défaut résiduel et l'écart du fonctionnel énergétique.
La projection duale retire les réactions par réflecteurs avant les normes,
ce qui évite la soustraction de deux grands Grams presque égaux.

Les 37 tests algébriques et physiques comprennent 12 contre-épreuves
rationnelles, des masses couplées, des ports sans masse, des directions
obliques, des charges qui se compensent, l'enrichissement et des pôles
globaux. Ces tests ne constituent pas une preuve de tous les arrondis.

## Protocole comparatif

La campagne contient 60 processus frais : trois consoles natives de
32, 128 et 512 poutres, cinq variantes, une chauffe puis trois mesures.
L'ordre des variantes s'inverse à chaque passage. Chaque essai traite
257 fréquences de 0 à 40 Hz et six charges terminales. CPU 8, bibliothèques
numériques à un fil, Python 3.14.7, NumPy 2.5.3 et SciPy 1.18.1 sont figés.

Chaque candidat recommence le QR, la sélection de la direction intérieure,
le produit définissant B, le certificat dirigé et la base de Krylov.
Le coût total comprend cette préparation, le contrôle éventuel, les
réponses et leurs normes physiques. Les imports, lectures, empreintes,
sauvegardes et le juge indépendant sont hors chronomètre pour tous.
Aucune direction ou certification préparée n'est offerte au candidat.

Les témoins sont la LU corrigée et les bases HCB de l'API officielle
d'Exudyn 1.11.0, avec le pont de réponse historique. Les rangs HCB
185, 761 et 768 proviennent de la confrontation antérieure ; la recherche
de ces rangs leur est offerte. Leur préparation est recomptée. Les deux
projections historiques, standard et énergétique, restent distinctes.
Ce protocole ne mesure pas une simulation FFRF ni le moteur Exudyn entier.

Le juge impose une erreur relative de 10⁻⁶ sur les champs en masse,
les déformations et les ports, par charge et sur toutes les combinaisons
des six charges. L'oracle Decimal à 70 et 90 chiffres donne les mêmes
champs après conversion en binary64 ; ce constat ne vaut pas un oracle
à intervalles certifié. Les majorants sont aussi confrontés aux erreurs
absolues observées, avec la marge flottante documentée `64 eps × norme
de référence`. Leurs booléens d'acceptation sont recalculables.

Les mesures mémoire des pilotes n'ont pas une convention commune ;
aucun classement mémoire n'est tiré de cette campagne. Tous les échecs
de calcul, de précision ou de contrôle font partie des résultats.

## Résultats à précision commune

Secondes, médianes de trois mesures après une chauffe. Chaque total est
la médiane des totaux par essai ; les médianes des phases peuvent donc
ne pas s'additionner exactement. La qualification inclut la chauffe.

| Poutres | Méthode | Préparation | Réponses et normes | Total | Champs à 10⁻⁶ |
|---:|---|---:|---:|---:|:---:|
| 32 | Krylov, champs | 0,0692 | 0,0294 | 0,0985 | 4/4 |
| 32 | Krylov, contrôle | 0,0777 | 0,1226 | 0,2006 | 4/4 |
| 32 | LU corrigée | 0,0011 | 0,1011 | 0,1022 | 4/4 |
| 32 | HCB standard | 0,0307 | 0,0781 | 0,1088 | 4/4 |
| 32 | HCB énergie | 0,0303 | 0,0780 | 0,1084 | 4/4 |
| 128 | Krylov, champs | 0,2193 | 0,0437 | 0,2634 | 4/4 |
| 128 | Krylov, contrôle | 0,2408 | 0,1753 | 0,4161 | 4/4 |
| 128 | LU corrigée | 0,0012 | 0,2219 | 0,2231 | 4/4 |
| 128 | HCB standard | 0,6627 | 0,9725 | 1,6401 | 2/4 |
| 128 | HCB énergie | 0,6541 | 0,9857 | 1,6363 | 4/4 |
| 512 | Krylov, champs | 0,8382 | 0,0988 | 0,9373 | 4/4 |
| 512 | Krylov, contrôle | 0,9210 | 0,3817 | 1,3026 | 4/4 |
| 512 | LU corrigée | 0,0018 | 0,7236 | 0,7254 | 4/4 |
| 512 | HCB standard | 4,9743 | 1,5347 | 6,5273 | 0/4 |
| 512 | HCB énergie | 4,9644 | 1,5292 | 6,4797 | 0/4 |

À 128 poutres, le service avec contrôle prend **0,416 s contre 1,636 s**
pour HCB énergie, soit 3,93 fois plus vite au total ; les deux passent.
Il reste 1,87 fois plus coûteux que la LU corrigée. À 32 poutres, les
champs seuls gagnent environ 3,6 % au total face à la LU dans ces trois
mesures ; ce petit écart n'établit pas un avantage robuste sur d'autres
machines. Le service avec contrôle est plus lent sur ce petit cas.

À 512 poutres, HCB dépasse le seuil avec une erreur maximale d'environ
3,23 × 10⁻⁵, pour les deux projections au rang historique 768. **Aucun
facteur d'accélération à précision commune n'est calculé contre ces
échecs.** Cela ne démontre pas qu'un autre rang ou réglage Exudyn échouerait.
À 128 poutres, HCB standard passe deux mesures sur trois et échoue à la
chauffe ; ce réglage reste exclu du classement admissible.

Le complément de Krylov compte 24, 32 et 40 directions : avec les six
ports et une direction intérieure retenue, les systèmes des champs ont
**31, 39 et 47 inconnues**, contre 192, 768 et 3 072 physiques. La variante
contrôlée résout en outre un bloc conservé de taille 7. Sur les trois
mesures, les pires erreurs physiques des champs seuls sont respectivement
1,30 × 10⁻¹⁰, 5,75 × 10⁻¹⁰ et 5,86 × 10⁻⁹ ; celles de la variante
contrôlée restent sous 2,20 × 10⁻⁹ sur l'ensemble des cas.

Tous les contrôles passent les 257 fréquences et les six charges dans
leurs quatre essais. Les maxima de leurs majorants relatifs, chauffe
comprise et arrondis vers le haut, sont :

| Poutres | Masse | Déformation | Limite spectrale du complément |
|---:|---:|---:|---:|
| 32 | 1,30 × 10⁻⁹ | 2,82 × 10⁻⁹ | ≈ 65,46 Hz |
| 128 | 1,23 × 10⁻⁹ | 1,35 × 10⁻⁹ | ≈ 65,70 Hz |
| 512 | 1,77 × 10⁻⁸ | 1,07 × 10⁻⁸ | ≈ 65,71 Hz |

Ces majorants concernent les charges demandées. Le juge indépendant
éprouve aussi l'opérateur sur toutes leurs combinaisons ; il ne transforme
pas les majorants par charge en certificat d'opérateur uniforme.

## Priorité après cette mesure

Sur 512 poutres, la préparation contrôlée se répartit approximativement
en 0,085 s de QR, 0,022 s de sélection de direction, **0,685 s de
certification**, 0,046 s de Krylov et 0,083 s d'enveloppe de contrôle.
Le certificat représente environ 74 % de la préparation et 53 % du total.
Les champs seuls restent plus lents au total que la LU sur 128 et 512
poutres. Optimiser seulement le petit système de réponse ne suffira pas.

La priorité issue de cette mesure était un encadrement creux du seuil spectral moins coûteux,
notamment par le [critère d'inertie déjà dérivé](INERTIE_COMPLEMENT_PREUVES.md),
avec preuve des arrondis et mesure complète. Le
[prototype d'inertie dirigée](INERTIE_CONTRAINTE_PROTOTYPE.md) réalise
désormais cette étape. La
[compression commune des résidus](KRYLOV_CONTRAINT_COMPRESSION.md) offre
ensuite une réduction du coût des normes, à évaluer sans changer la
physique. Une réutilisation de préparation pourra être mesurée sur un
usage qui la justifie ; elle ne doit pas être offerte dans un comparatif
qui impose une préparation fraîche.

## Reproduire et auditer

L'[archive compacte](bancs/krylov-contraint-2026/README.md) conserve les
sources figées, les modèles, les bases physiques, B/Phi, les petits K/M,
les résultats détaillés, journaux et empreintes des grands champs exclus.
Les [sondes préalables](bancs/krylov-contraint-sondes-2026/README.md)
conservent les échecs qui ont conduit aux quatre blocs et à la réparation
finale. Leurs anciens états de code ne sont pas tous figés : elles ne
servent pas de classement formel.

Pour auditer l'archive depuis le dépôt avec NumPy/SciPy disponibles :

```bash
python ci/archive_krylov_contraint.py --verifier docs/bancs/krylov-contraint-2026
```

Pour reproduire la campagne, fournir un dossier neuf, les modèles et les
oracles aux SHA historiques, deux environnements figés Vinkulum 0.11.0
et Exudyn 1.11.0, puis utiliser le pilote archivé :

```bash
"$PY_VINKULUM" docs/bancs/krylov-contraint-2026/sources/experience_krylov_contraint.py \
  --campagne "$SORTIE_NEUVE" "$ENTREES_AVEC_ORACLES" \
  "$PY_VINKULUM" "$PY_VINKULUM" "$PY_EXUDYN" \
  --historique docs/bancs/confrontation-ports-exudyn-0.10.0/bilan.json
```

Le pilote refuse une sortie déjà présente et ne relance pas les essais
interrompus. Les installations restent identifiées par leurs versions,
roues et fichiers exécutables ; les résultats d'une reproduction doivent
être qualifiés avant de comparer ses durées.

## Limites de transfert

Les consoles sont linéaires, non amorties, à masse diagonale. Les tests
du contrôleur couvrent aussi des masses intérieures SPD par composantes
connexes de taille au plus six et une masse complète PSD. Une masse
consistante plus largement couplée exige une autre racine creuse.

L'assemblage de plusieurs sous-structures, les grandes rotations,
la précontrainte, les contacts et les trajectoires restent à éprouver.
L'usage d'une matrice réduite beaucoup plus petite ne suffit pas à établir
un avantage de temps total. Cette expérience ne classe pas Vinkulum face
aux moteurs multicorps complets MBDyn, Exudyn ou Simpack.
