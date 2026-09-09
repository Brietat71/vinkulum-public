# Inertie dirigée : certificat moins coûteux et réduction à 40 Hz

Cette note conserve la campagne de 84 essais. Le
[contrôle par facteurs](CONTROLE_FACTEURS_PROTOTYPE.md) apporte ensuite
une campagne indépendante de 60 essais : à 512 poutres, le service avec
majorants descend de 0,819 à 0,671 s, contre 0,696 s pour la LU du même lot.
Les résultats ci-dessous restent ceux de la campagne d'inertie initiale.

Point expérimental du 8 septembre 2026. Le certificat creux du complément
est désormais calculé par congruences à arrondis dirigés. Sur les trois
consoles mesurées, les champs passent la précision relative 10⁻⁶ et leur
calcul complet devient **1,46 à 1,58 fois plus rapide que la LU corrigée**.
À 512 poutres, le total passe de **0,937 à 0,451 s**. Avec le contrôle
d'erreur, il passe de **1,305 à 0,816 s**, mais la LU reste plus rapide
(0,711 s). Toutes ces durées comprennent une préparation fraîche.

La campagne contient **84 essais**, dont 24 avec le nouveau certificat,
24 avec l'ancien et 36 témoins. Les 48 essais des prototypes passent le
juge physique ; les 24 avec contrôle passent aussi leurs majorants.
Au total, 75/84 essais passent : les neuf refus HCB restent archivés.
Cette expérience n'étend pas encore l'API publique de la **0.11.0**.

## Changement mathématique

Soient les facteurs matériels intérieurs D, la masse M et les covecteurs
stockés B. Chaque nombre binary64 est interprété exactement. On assemble
par intervalles le système symétrique

```text
Aγ = D.T D − γ M,
Cγ = [[Aγ, B], [B.T, 0]],       B de taille n × s.
```

Une élimination complète dont les pivots prouvent la signature `(n,s,0)`
établit simultanément le rang s de B et la coercivité stricte

```text
x.T D.T D x > γ x.T M x    pour x ≠ 0 et B.T x = 0.
```

L'identité d'inertie par congruence conserve les directions négatives
associées aux contraintes et isole celles du complément. Elle reste
valable lorsque Aγ est singulière. Aucun inverse complet, soustraction
de traces ou base dense du noyau n'est nécessaire.
La [preuve algébrique](INERTIE_COMPLEMENT_PREUVES.md) et le
[carnet d'arrondis](INERTIE_DIRIGEE_PREUVES.md) précisent l'invariant de
Schur, les pivots doubles et les refus.

Ce principe n'est pas un nouveau théorème. La
[recherche ciblée dans cinq sources primaires](INERTIE_DIRIGEE_SOURCES.md)
relie les vérifications creuses de 1995 aux travaux de 2025–2026 sur
l'inertie et les résidus de factorisation. Notre transfert concret est
le certificat de **D.T D exact sur le complément contraint**, composé
avec la réduction énergétique existante et confronté aux champs physiques.
Le repli par décalage unilatéral dérivé dans cette recherche reste une
piste : il n'est pas utilisé dans les durées publiées ici.

## Contrat de l'implémentation

[ci/inertie_complement_dirigee.py](../ci/inertie_complement_dirigee.py)
certifie d'abord M strictement positive. Une masse exactement diagonale
positive est reconnue sans opération arithmétique ; une masse couplée
passe une élimination dirigée séparée. La seule positivité des diagonales
ne suffirait pas. Cette certification de masse n'impose pas de blocs de
taille six ; cette restriction appartient encore à la racine de masse
employée par le contrôleur de réponse.

Le Gram est assemblé directement depuis D, avec produits et sommes
dirigés. Certifier un Gram déjà arrondi en double changerait le problème.
Les pivots scalaires doivent avoir un signe strictement séparé ; sinon
un pivot 2×2 voisin est essayé, avec déterminant et trace encadrés.
La mise à jour encadre le Schur exact, même lorsque des occurrences
partagent les mêmes variables. Seul un intervalle exactement `[0,0]`
peut être retiré. Aucun pivot régularisé ni coefficient simplement petit
n'est substitué à une donnée d'entrée.

Une signature complète différente réfute la coercivité au seuil demandé.
Un pivot indécidable ou un budget dépassé produit un refus distinct,
sans conclusion négative sur le problème exact. La campagne fixe avant
son lancement **32 chiffres Decimal**, l'ordre physique naturel, puis
les multiplicateurs, et les budgets suivants : 100 millions d'opérations
dirigées, deux millions de coefficients symétriques et deux millions
d'entrées rectangulaires de B. Ces budgets sont des compteurs logiques,
pas une limite de mémoire en octets : les dictionnaires conservent les
deux triangles et la validation des entrées peut aussi allouer.

Le seuil fixé est `γ = 252661.87266788757`, soit le binary64
`0x1.ed7aefb394d2bp+17`, construit à partir de 80 Hz. La comparaison
`Fraction(γ) > Fraction(Ω)**2`, avec Ω le maximum stocké de la bande
0–40 Hz, est vérifiée exactement. « 80 Hz » décrit ce choix ; la preuve
porte sur le γ stocké. L'ancien minorant correspondait à environ 65,7 Hz.
Le certificat est donc moins coûteux et son minorant est plus élevé.

Les [sondes préalables](bancs/inertie-contrainte-sondes-2026/README.md)
conservent notamment le refus à 16 chiffres sur 512 poutres. Elles ont
aussi éprouvé trois blocs de Krylov : les majorants échouaient sur deux
maillages. Les **quatre blocs**, une direction intérieure retenue et une
profondeur de contrôle de huit sont donc fixés avant les 84 essais.
Aucune reprise ni changement de précision pendant la campagne.

## Protocole et résultats

Trois consoles matérielles linéaires de 32, 128 et 512 poutres, soit
192, 768 et 3 072 coordonnées libres ; 257 fréquences de 0 à 40 Hz ;
six charges terminales. Chaque variante utilise un processus neuf pour
une chauffe et trois mesures, avec ordre inversé à chaque passage.
CPU 8, bibliothèques à un fil, Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1,
roues Vinkulum 0.11.0 et Exudyn 1.11.0 identifiées par leurs empreintes.

Les deux nouveaux candidats diffèrent des deux candidats Krylov par le
certificat. Les douze sources de l'ancienne campagne restent identiques.
Chaque processus candidat recalcule QR, direction intérieure, B, certificat,
base réduite, contrôle éventuel, réponses et normes physiques. Les imports,
lectures, captures de provenance, sauvegardes et juge indépendant restent
hors chronomètre. Le chronomètre extérieur du certificat compte son appel
entier, y compris les empreintes de son résultat.

Les bases HCB proviennent de l'API officielle Exudyn avec le pont de
réponse historique. Les rangs 185, 761 et 768, recherchés antérieurement,
sont offerts aux témoins ; leur construction est recomptée. Les projections
standard et énergétique sont conservées. Le banc n'exécute pas une
simulation FFRF complète. Aucun code d'implémentation concurrent n'est lu.

Le juge contrôle les champs en masse, déformations et ports, par charge
et sur toutes les combinaisons des six charges, au seuil relatif 10⁻⁶.
Les champs de référence Decimal à 70 et 90 chiffres coïncident après
conversion en double ; ce constat ne certifie pas l'oracle par intervalles.

Durées en secondes, médianes des trois mesures. Les médianes de phases
peuvent ne pas s'additionner exactement à celle du total. L'admission
compte aussi la chauffe ; une seule occurrence refusée exclut un ratio
de vitesse à précision commune pour cette configuration.

| Poutres | Variante | Préparation | Réponses | Total | Admissions |
|---:|---|---:|---:|---:|---:|
| 32 | Inertie, champs | 0,0415 | 0,0294 | **0,0708** | 4/4 |
| 32 | Inertie, contrôle | 0,0495 | 0,1215 | 0,1710 | 4/4 |
| 32 | Trace, champs | 0,0696 | 0,0297 | 0,0991 | 4/4 |
| 32 | Trace, contrôle | 0,0774 | 0,1224 | 0,1998 | 4/4 |
| 32 | LU corrigée | 0,0010 | 0,1024 | 0,1035 | 4/4 |
| 32 | HCB standard | 0,0308 | 0,0781 | 0,1083 | 4/4 |
| 32 | HCB énergie | 0,0304 | 0,0775 | 0,1084 | 4/4 |
| 128 | Inertie, champs | 0,1022 | 0,0432 | **0,1454** | 4/4 |
| 128 | Inertie, contrôle | 0,1233 | 0,1741 | 0,2980 | 4/4 |
| 128 | Trace, champs | 0,2201 | 0,0436 | 0,2637 | 4/4 |
| 128 | Trace, contrôle | 0,2407 | 0,1742 | 0,4162 | 4/4 |
| 128 | LU corrigée | 0,0012 | 0,2206 | 0,2218 | 4/4 |
| 128 | HCB standard | 0,6599 | 0,9742 | 1,6397 | 3/4 |
| 128 | HCB énergie | 0,6554 | 0,9715 | 1,6269 | 4/4 |
| 512 | Inertie, champs | 0,3534 | 0,0980 | **0,4514** | 4/4 |
| 512 | Inertie, contrôle | 0,4348 | 0,3819 | 0,8159 | 4/4 |
| 512 | Trace, champs | 0,8382 | 0,0982 | 0,9367 | 4/4 |
| 512 | Trace, contrôle | 0,9244 | 0,3813 | 1,3049 | 4/4 |
| 512 | LU corrigée | 0,0018 | 0,7096 | 0,7114 | 4/4 |
| 512 | HCB standard | 4,9993 | 1,5242 | 6,5236 | 0/4 |
| 512 | HCB énergie | 4,9477 | 1,5436 | 6,4880 | 0/4 |

Les champs avec inertie gagnent respectivement **1,46×, 1,53× et 1,58×**
sur la LU au total. À 128 poutres, HCB énergie est admissible : les gains
sont **11,19× pour les champs** et **5,46× avec contrôle**. À 32 poutres,
les champs devancent HCB, mais le service avec contrôle reste plus lent.
À 512 poutres, aucun ratio HCB n'est recevable : ses deux configurations
dépassent 10⁻⁶, avec une erreur maximale voisine de 3,23 × 10⁻⁵.
À 128 poutres, une des trois mesures HCB standard atteint 1,034 × 10⁻⁶ ;
les deux autres et la chauffe passent. Aucune occurrence n'est écartée.

Les tailles réduites restent 31, 39 et 47. La pire erreur du nouveau
candidat champs vaut respectivement 9,94 × 10⁻¹¹, 6,89 × 10⁻¹⁰ et
5,63 × 10⁻⁹ (arrondis vers le haut). Le candidat avec contrôle reste
en dessous de 2,17 × 10⁻⁹ sur ces trois cas. Les majorants couvrent les
erreurs observées avec la marge flottante documentée `64 eps × norme
de référence` ; ils ne sont pas des certificats machine des réponses.

## Coût gagné, coût restant

Sur le candidat champs, la certification passe de 0,0401/0,1669/0,6845 s
à **0,0121/0,0490/0,1992 s**, soit environ **3,3 à 3,4 fois plus rapide**.
Le grand KKT n'a que quatre voisins futurs au maximum dans cet ordre
d'élimination. Il crée 12 255 coefficients symétriques cumulés et utilise
204 206 opérations Decimal. Ce faible remplissage est une propriété de
ces consoles et de cet ordre, pas une complexité générale des KKT.
Les 24 certificats d'inertie de la campagne utilisent uniquement des
pivots scalaires ; les pivots 2×2 sont éprouvés par les tests séparés.

Le contrôle reste la limite du gain total : à 512 poutres, environ
0,198 s de certificat, 0,083 s de préparation du contrôle et 0,382 s de
réponses contrôlées. Ces dernières deviennent la plus grosse phase.
La [compression des résidus](KRYLOV_CONTRAINT_COMPRESSION.md) doit être
éprouvée avec ses défauts d'approximation ; économiser seulement sa
préparation ne suffit pas à établir une victoire sur la LU. Le coût des
normes physiques et du contrôle par fréquence doit aussi être traité.

Le certificat spectral couvre le complément pour toute la bande sous
son seuil. Les réponses, la marge du système conservé et les réparations
sont évaluées en doubles aux fréquences demandées. Les majorants par
charge ne certifient pas uniformément toutes les combinaisons ; celles-ci
sont éprouvées séparément par le juge. Aucune garantie de réponse entre
les points échantillonnés ni sur le continuum mécanique n'est ajoutée.

Le transfert public, les masses consistantes plus largement couplées
dans le contrôleur, l'amortissement, la précontrainte, les contacts et les
trajectoires non linéaires restent ouverts. Aucun classement mémoire
commun ni nouvelle mesure MBDyn/Simpack n'est produit par ce banc.
L'objectif généraliste reste actif.

## Preuves, archive et reproduction

Les 18 tests rationnels indépendants éprouvent les congruences et les
frontières ; 12 autres tests confrontent l'implémentation à des calculs
Fraction et à des cas hostiles. Ils incluent masse couplée 7×7, permutations,
contraintes obliques et déficientes, données extrêmes, budgets, pivots 2×2
et un faux positif qu'aurait produit un Gram arrondi.

L'[archive des 84 essais](bancs/inertie-contraint-2026/README.md) conserve
modèles, bases physiques, B, résultats, journaux et sources figées. Les
24 certificats d'inertie correspondent à trois identités distinctes D/M/B ;
le validateur rejoue chacune et compare toute la preuve, hormis ses temps.
Les grands champs exclus restent identifiés par empreinte ; leur jugement
est requalifié depuis les diagnostics conservés, sans reconstruire ces
champs. Consulter le README de l'archive pour la portée exacte de chaque
vérification.

Après la campagne, le pilote courant reçoit un correctif limité à la
sérialisation des refus : il conserve aussi les détails des exceptions
qui exposent `.bilan` au lieu de `.diagnostic`. Le pilote archivé reste
celui des 84 essais ; aucun calcul numérique mesuré n'est modifié.

```bash
python ci/archive_inertie_contrainte.py --verifier docs/bancs/inertie-contraint-2026
```

Pour une nouvelle campagne, les environnements et données doivent respecter
les identités publiées ; fournir un dossier de sortie neuf :

```bash
"$PY_VINKULUM" docs/bancs/inertie-contraint-2026/sources/experience_inertie_contrainte.py \
  --campagne "$SORTIE_NEUVE" "$ENTREES_AVEC_ORACLES" \
  "$PY_VINKULUM" "$PY_VINKULUM" "$PY_EXUDYN" \
  --historique docs/bancs/confrontation-ports-exudyn-0.10.0/bilan.json
```

Le pilote refuse un dossier existant et ne reprend pas une campagne
interrompue. Les temps d'une reproduction ne deviennent comparables
qu'après sa propre qualification physique.
