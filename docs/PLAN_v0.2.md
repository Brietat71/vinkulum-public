# Vinkulum v0.2 — le noyau, pas le squelette

Décision de Paul, 4 septembre 2026 : « préparer une prochaine version beaucoup
plus crédible et sérieuse, qui dépasse déjà la concurrence en performances
fondationnelles ». Ce plan dit ce que « fondationnel » veut dire, comment on le
MESURE, et dans quel ordre on le construit. Rien n'y est un slogan : chaque
ligne a un chiffre de sortie et un juge.

## 1. Ce qu'on juge — les performances fondationnelles

Un noyau se compare sur des courbes **travail–précision** (erreur contre temps
CPU), pas sur des temps seuls ni sur des listes de fonctionnalités. Les cinq
axes, avec la référence à battre et l'état v0.1 :

| axe | mesure | v0.1 (3 sept.) | cible v0.2 | référence |
|---|---|---|---|---|
| coût par pas à précision donnée | ms/pas pour 1e-6 rad sur la tête S2 (143 inconnues) | 12 | **≤ 0,3** | MBDyn même deck : à mesurer, ordre 0,1–0,5 |
| passage à l'échelle | temps/pas en fonction de N, de 10² à 10⁴ inconnues | dense O(N³) | **O(N^1,2)** creux, mesuré | Exudyn / Chrono |
| exactitude du jacobien | itérations de Newton par pas ; convergence quadratique | 2–4, linéaire (DF) | **≤ 2**, quadratique (AD) | analytique chez tous |
| parallélisme | accélération forte 1 → 64 cœurs sur 10⁴ inconnues | colonnes DF seulement | **≥ 20×** sur l'assemblage et le résidu — *mesuré v0.2.0 : ×6 à 16 fils sur l'assemblage ; le solveur direct ne parallélise pas sur une CHAÎNE (arbre d'élimination = chemin), cible à rejuger sur un arbre large* | MBDyn : quasi nul ; Exudyn : partiel |
| robustesse | pas rejetés / échecs Newton sur la suite de bancs, 10⁵ pas | non mesuré | **0 échec**, pas adaptatif avec estimateur d'erreur | — |

Et deux axes de **crédibilité** qui ne sont pas de la vitesse :

- **exactitude** : ordre 2 en q, v, λ tenu sur TOUS les bancs ; violation en
  vitesse publiée (O(h²)) et une option index 2 (GGL) mesurée à côté ;
- **contre-solveurs** : chaque banc a une référence indépendante (analytique,
  EDO minimale, MBDyn sur le même deck, valeurs publiées).

## 2. La suite de bancs — le juge de tout le reste

Elle existe AVANT le code qu'elle juge. Trois familles :

1. **Analytiques** (v0.1, gardés) : pendule, toupie, bielle-manivelle imposée,
   engrenage.
2. **Bancs publiés du métier** (IFToMM multibody benchmark, valeurs de
   référence publiées) : mécanisme d'Andrews (boucles fermées, raide en
   temps), **mécanisme de Bricard** (contraintes REDONDANTES — celles qu'on
   n'a jamais jugées), régulateur à boules (gyroscopie + ressorts), pendule
   double chaotique (sensibilité, ordre), rotor de Jeffcott (raideur).
3. **La tête S2 de FRELON** en dynamique LIBRE : gouverneur, servos PD en
   couple, puis aéro c81 — contre le vol MBDyn image par image. C'est le seul
   banc à 10² inconnues et boucles fermées multiples ; on y mesure le coût par
   pas et la précision en même temps.

Sortie : `python -m vinkulum.bancs` produit les courbes travail–précision
(JSON + figure) pour vinkulum ET pour MBDyn quand le deck existe. « Dépasser »
se lit sur ces courbes, nulle part ailleurs.

## 3. Ce qu'on construit, dans l'ordre — chaque brique jugée avant la suivante

### 3.1 Jacobien exact par différentiation automatique (le poste n°1)

Les différences finies coûtent N résidus par jacobien et bornent Newton à une
convergence linéaire. Le résidu est une composition de petites fonctions par
élément (Φ, G, M, f) : **mode avant par nombres duaux** (`ad-trait`, tangentes
SIMD, sur Rust stable — Enzyme quand il sera stable), par ÉLÉMENT : chaque
liaison rend son bloc jacobien local (n_e × 12 ou 18) avec les dérivées de G
et de f incluses, exactes à l'arrondi. Coût : O(nnz) au lieu de O(N·coût du
résidu). Juge : convergence quadratique de Newton (résidu ÷ 100 par
itération), et coïncidence avec les DF à 1e-7.

### 3.2 Assemblage creux et solveur creux (le poste n°2)

Système de Newton bloc-creux : M bloc-diagonal, Gᵀ par éléments. Assemblage
**par fils** (rayon sur les éléments, tampons COO par fil → CSR), solveur
direct creux **faer** (LU/LDLᵀ symétrique indéfini avec pivotage), mise à
l'échelle des lignes de contrainte (déjà là). Contraintes redondantes :
détection de rang à l'initialisation (SVD de G₀) et **régularisation** des
lignes dépendantes — jugée sur Bricard. Juge : O(N^1,2) mesuré sur une chaîne
de 10⁴ inconnues, résultat identique au dense à 1e-12.

### 3.3 Données en SoA, résidu et assemblage parallèles

Corps et états en tableaux contigus (positions, R, v, ω, masses, tenseurs),
éléments par type en tableaux ; boucles `par_iter` sur les éléments. Juge :
accélération forte mesurée à 1/4/16/64 fils, et déterminisme au bit conservé
(réduction ordonnée).

### 3.4 Intégrateur : pas adaptatif, initialisation consistante, option index 2

Estimateur d'erreur local (différence entre le prédicteur et la solution ;
Arnold–Brüls) → contrôle du pas, rejet, redémarrage ; initialisation
consistante des accélérations et de λ (déjà) ; **option GGL** (index 2,
contrainte en vitesse tenue) mesurée sur les mêmes bancs pour publier ce
qu'elle coûte et ce qu'elle rend. Juge : 0 échec sur 10⁵ pas de la suite,
ordre 2 conservé sous pas variable.

### 3.5 Éléments de force et de commande

Ressort-amortisseur (translation, rotation), **servo PD saturé en couple et
en vitesse** (le patron du n0 FRELON), gouverneur (couple asservi à une
consigne de vitesse), couple/force suiveur. Juge : la tête S2 en dynamique
libre contre MBDyn — régime, droop, T/W.

### 3.6 Aérodynamique de pale : élément c81 + inflow

Élément de pale (stations, c81, inflow uniforme puis Pitt-Peters) sur un corps
rigide — le niveau 1 de MBDyn, chez nous. Juge : le vol MBDyn (T/W 0,89 à
inflow uniforme), puis le BEMT calibré V-2 de FRELON.

### 3.7 Poutre géométriquement exacte (Simo–Reissner sur SO(3))

Le premier corps flexible, sur le même groupe que l'intégrateur — la pale.
Juge : cas de Bathe (cantilever 45°), fréquences propres d'une poutre en
rotation (raidissement centrifuge — Southwell), puis la pale FRELON contre
`battement.py` (ν de battement 1,28–1,40).

### 3.8 Interfaces : zéro-copie, artefacts signés, format de sortie

`rust-numpy` (vues sans copie), replay signé par les cotes consommées
(patron FRELON), sortie colonnaire (Arrow) pour les longues trajectoires.

## 4. Ce qu'on NE fait pas en v0.2, et pourquoi

Contact, corps flexibles 3D (solides), co-simulation, GPU, géométrie B-Rep :
aucun banc de FRELON ne les exige encore, et chacun est une brique à juger à
part. Le SDF comme représentation maîtresse : refusé (README). Enzyme : quand
stable.

## 5. Ordre de marche et jalons de sortie

| jalon | contenu | sortie (mesurée) |
|---|---|---|
| **v0.2.0** | suite de bancs + AD + creux + SoA | courbes travail–précision, S2 ≤ 0,3 ms/pas, Bricard tenu |
| **v0.2.1** | pas adaptatif + GGL en option + éléments de force | S2 en dynamique libre = MBDyn, 0 échec sur 10⁵ pas |
| **v0.2.2** | aéro c81 + inflow | vol S2 = MBDyn image par image ; T/W, régime, droop |
| **v0.3.0** | poutre géométriquement exacte | pale FRELON flexible contre `battement.py` |

La règle ne change pas : **une brique n'entre que jugée**, sa référence
indépendante nommée, ses chiffres écrits ici et dans le README. Ce que la
concurrence a et que ce plan n'a pas se dit ; ce que ce plan rend et qu'elle
n'a pas se mesure.
