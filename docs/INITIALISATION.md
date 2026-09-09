# Initialisation cohérente : courbure, redondances et décollement

Travail de la corrective **0.7.2**, comparé à la roue figée **0.7.1**.
Le périmètre est l'accélération initiale utilisée au départ de `simule`,
au redémarrage après un impact et dans le calcul des réactions de
précontrainte de `k_c_m_z()` lorsque les multiplicateurs manquent.

## Défauts reproduits

Une bielle de longueur `L` impose à une masse animée d'une vitesse
tangentielle `v` une accélération radiale `−v²/L`. L'ancien calcul de
courbure déplaçait les poses pendant `ε = 10⁻⁷ s` puis différenciait `G u`.
Dans ce cas, il donnait `−v²/√(L² + ε²v²)`. À **L = 1 µm** et **v = 10 m/s**,
le résultat était `−7,071067811865476e7 m/s²`, contre `−1e8 m/s²` attendu :
**29,289 % d'erreur**. La nouvelle évaluation donne le résultat analytique
exactement dans l'arrondi f64 de ce cas. Cela ne promet pas une exactitude
universelle à toutes les échelles.

L'analyse d'une cascade connexe de parallélogrammes pouvait échouer dès
**24 corps**, avec `LU creux : SymbolicSingular`. Ses dépendances sont
compatibles : huit paramètres angulaires décrivent les mouvements du
modèle. Le défaut venait de l'initialisation des multiplicateurs utilisée
par la raideur de précontrainte. Le noyau exigeait un LU inversible au-delà
du seuil dense.

Après un rebond sans frottement, l'ancienne initialisation imposait encore
une accélération normale nulle à une bille dont la vitesse l'éloignait du
plan. Le redémarrage doit retrouver `−g`. Le diagnostic à un pas mesure
**0 → −9,81 m/s²**. Une première combinaison de la courbure exacte avec
l'ancien ensemble de contacts a aussi révélé une accélération attractive
au contact rasant de deux sphères : ce prototype n'est pas livré.

Les [références antérieures](bancs/initialisation-0.7.1-references.json)
identifient la roue témoin par l'empreinte de son extension. Le diagnostic
commun est [exécutable](../ci/diagnostic_initialisation.py) avec chaque roue.

## Dérivée directionnelle sur les poses

Pour une contrainte holonome, le terme convectif est

```text
γ = d²/ds² Φ(q ⊕ s u, t+s), évalué en s = 0,
G u̇ = −γ.
```

Les positions suivent `r(s) = r + s v`, les orientations
`R(s) = exp(s[ω]×) R`. Deux nombres duaux emboîtés, chacun avec un seul
axe, calculent cette dérivée mixte. Les petites translations compensées
restent dans les poses. Aucune pose du modèle n'est déplacée pour l'essai
et aucune Hessienne globale n'est assemblée.

Pour les contraintes non holonomes linéaires en vitesse, la quantité
différenciée est directement `G(q ⊕ s u,t+s) u + Φ_t` : elle correspond à
la loi de vitesse de l'élément. Les commandes sont différenciées avec la
géométrie, sur le segment courant de leur loi. Deux rotations commandées
`R_b(t) = R_x(2t) R_z(3t)` donnent ainsi `α_b(0) = −6 e_y rad/s²` sans les
composantes parasites de l'ancien déplacement d'essai.

Les lois tabulées sont affines par morceaux. À une cassure de pente, il
n'existe pas de dérivée seconde classique ; le calcul conserve la branche
choisie par `Loi::derivee`. Il ne résout pas une impulsion de commande.
Pour le contact, les changements de primitive ou de point le plus proche
conservent leurs limites de différentiabilité existantes.

## Projection dans la métrique de masse

La formulation idéale de Gauss choisit l'accélération compatible la plus
proche de l'accélération libre dans la métrique de masse. Le cadre, ses
hypothèses et une extension aux contraintes inconsistantes sont présentés
par [Udwadia, 2023](https://ruk.usc.edu/bio/udwadia/papers/The%20General%20Gauss%20Principle%20of%20Least%20Constraint%20JAM%202023.pdf).
Vinkulum conserve ici des contraintes compatibles : il ne remplace pas
silencieusement une incompatibilité par un ajustement des contraintes.

Pour `M = L Lᵀ`, on pose `W = L⁻¹ Gᵀ`, `z₀ = L⁻¹ f` et `c = −γ`.
Avec la décomposition mince `W = U Σ Vᵀ` et son rang numérique `r` :

```text
z_particulier = U_r Σ_r⁻¹ V_rᵀ c,
z_libre      = (I − U_r U_rᵀ) z₀,
u̇           = L⁻ᵀ (z_libre + z_particulier),
λ            = V_r Σ_r⁻¹ (U_rᵀ z₀ − Σ_r⁻¹ V_rᵀ c).
```

Le solveur ne forme pas `G M⁻¹ Gᵀ`. Il calcule séparément la composante
libre et la solution particulière, pour éviter de retrouver une petite
accélération par soustraction de grandes réactions. Les facteurs de masse
sont des blocs 3×3. Le graphe sépare les composantes indépendantes avant
les SVD ; les corps gelés suivent la réduction multirythme existante.
Les lignes explicitement inactives gardent un multiplicateur nul.

Le seuil de rang est `ε_machine max(n_bloc,m_bloc) ‖W_bloc‖F`. C'est un
seuil numérique, dépendant des unités et des échelles ; il n'est pas une
certification par intervalles. Le repli fournit la réaction de norme
euclidienne minimale dans le rang retenu. Le chemin LU conservé ne vérifie
pas indépendamment cette norme minimale pour une dépendance que ses
contrôles numériques laisseraient passer : l'unification avec les
projecteurs orthogonaux reste un travail à faire.

Le LU reste prioritaire lorsqu'aucune singularité structurelle immédiate
ne l'exclut. Un contrôle compensé sépare le résidu dynamique des équations
d'accélération. Jusqu'à trois corrections réutilisent les mêmes facteurs
avant le repli. Cette étape conserve le pendule vertical URDF : le premier
prototype le refusait parce que son accélération presque nulle contenait
un résidu d'arrondi. Le [journal de cet échec](bancs/initialisation-ci-premier-refus.log)
est conservé ; le critère n'a pas été relâché pour le faire passer.

## Sélection des contacts au redémarrage

Un contact dont la vitesse normale indique un décollement est libéré.
Pour les contacts pouvant rester fermés, le problème local impose
`λ ≥ 0`, `δ̈ ≤ 0` et `λ δ̈ = 0`, avec `δ = −gap`. Un premier calcul avec
tous les contacts admissibles est conservé si ses réactions sont positives.

Sinon, une sélection duale augmente la réaction d'une contrainte violée
jusqu'à sa fermeture ou jusqu'à l'annulation d'une autre réaction. Cette
seconde possibilité permet de relâcher un appui. Le point de départ libre
et la sélection duale sont des principes classiques de programmation
quadratique, décrits par
[Goldfarb et Idnani, 1983](https://link.springer.com/article/10.1007/BF02591962).
L'implémentation présente recalcule ses directions par projection de masse ;
elle n'implémente pas leurs mises à jour incrémentales de facteurs QR.

Les directions faibles sont contrôlées sans former de matrice normale.
Un budget de pivots borne le travail ; un échec de compatibilité, de signe
ou de résidu produit une erreur. Les tests comparent les accélérations
à l'énumération indépendante de toutes les faces sur 32 problèmes, avec
inertie anisotrope, ainsi qu'à des contraintes presque parallèles à `10⁻¹²`.
Un contact excentré sur un corps tournant vérifie aussi la courbure et la
réaction par une formule analytique.

Ce calcul concerne l'accélération lisse du redémarrage. La loi d'impact et
son schéma temporel restent ceux du noyau. Le frottement non lisse reste
évalué dans les forces avec le dernier multiplicateur disponible ; sa
résolution simultanée avec la nouvelle réaction n'est pas réalisée ici.

## Validation et mesures

La CI étendue passe : **66 tests Rust**, **72 tests Python**, **41/41 groupes
de vérification**, **46/46 bancs mécaniques**, **9/9 contacts**, formatage,
Clippy et référence des **185 entrées d'API**. Les
[journaux](bancs/initialisation-ci.log) et
[tests Python](bancs/initialisation-tests-python.log) sont conservés.

Les [mesures brutes](bancs/initialisation-mesures.json) comprennent **138
appels**, soit 23 configurations × deux versions × trois répétitions,
dans l'ordre AB / BA / AB. Chaque appel utilise un processus neuf, fixé
au CPU indiqué dans l'archive, avec un fil demandé. Le chronomètre couvre
l'appel public et, pour l'initialisation, la lecture de l'accélération du
schéma. Les imports et la construction du modèle précèdent le chronomètre ;
le contre-calcul indépendant le suit. Le pic RSS du processus est relevé
avant ce contre-calcul. Les coûts du premier appel restent inclus.

Les **69 essais candidats** passent ; **30 des 69 essais témoins** échouent
par singularité. Les références de cascade utilisent un paramétrage indépendant par les
angles : matrice de masse réduite, Hessienne du potentiel, mobilité,
orthogonalité et absence de mutation de l'état public.

| Appel, cascade | Corps | Temps médian 0.7.2 | Pic RSS maximal |
|---|---:|---:|---:|
| Analyse, repère initial | 24 | 43,17 ms | 73,4 Mio |
| Analyse, repère initial | 96 | 312,0 ms | 133,8 Mio |
| Analyse, repère initial | 144 | 828,2 ms | 214,3 Mio |
| Analyse, repère tourné | 144 | 848,6 ms | 214,6 Mio |
| Initialisation, repère initial | 144 | 328,1 ms | 138,5 Mio |
| Initialisation, repère tourné | 144 | 326,4 ms | 139,0 Mio |

Les versions témoins échouent sur ces appels : aucun rapport de vitesse
n'est calculé avec un échec. Sur les treize configurations valides dans
les deux versions, la variation des médianes va de **−4,1 % à +1,5 %**.
Le cas de 2 048 corps libres passe de **5,373 à 5,235 ms**, avec accélération
analytique exacte dans l'arrondi du test. Trois répétitions ne suffisent
pas à établir une absence générale de surcoût.

Sur les composantes d'accélération exprimées en SI, le plus grand écart
absolu aux coordonnées réduites est **7,42e-12** (translations en m/s²,
rotations en rad/s²). Le défaut maximal `G Z` est **1,27e-15**, celui de
`Zᵀ Z − I` **3,11e-15**. Les écarts de la raideur réduite à la Hessienne
du potentiel restent sous **2,85e-14 N·m/rad**. Ces seuils correspondent
aux modèles et états testés, pas à une borne générale d'erreur.

Les [références 0.7.2](bancs/initialisation-0.7.2-references.json), la roue
vérifiée hors du dépôt et ses empreintes figurent dans le
[relevé de livraison](bancs/version-0.7.2.json).

Les matrices d'analyse publiques et les SVD de composantes connectées
restent denses. Leur passage à l'échelle, les autres projecteurs, les
sensibilités aux changements de rang et les contacts avec frottement
couplé restent ouverts. Ces corrections ne constituent pas une nouvelle
confrontation exécutée avec MBDyn ou Simpack, ni une preuve de supériorité
générale sur ces solveurs.
