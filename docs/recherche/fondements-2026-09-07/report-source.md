# Vinkulum : les fondements d’une rupture de performance

## Recherche et décisions d’architecture · 7 septembre 2026

**Destinataire : conception scientifique du noyau généraliste.** Approfondissement de l’étude 2001–2026, au-delà du catalogue de méthodes. Référence logicielle : Vinkulum 0.9.0, commit 7850f6d. Les propositions de ce rapport sont des programmes de recherche ; elles ne sont pas des fonctionnalités livrées.

### La décision

Il faut déplacer le centre de gravité de Vinkulum vers **la dimension réellement nécessaire du calcul, la structure des couplages et des bornes d’erreur exploitables**. Une tangente plus rapide améliore le calcul existant. Une formulation qui évite de construire ou de résoudre une grande partie du problème peut changer son ordre de coût.

Trois programmes méritent de passer avant une nouvelle optimisation isolée : l’élimination locale des contraintes avec contrôle du rang ; les espaces de transfert optimaux entre composants ; une formulation mixte et un préconditionnement dont la stabilité résiste aux paramètres physiques. Les variétés adaptatives et la fermeture thermodynamique constituent un second horizon, plus risqué.

Le résultat de théorie des algorithmes de **Fürer–Hoppen–Trevisan, 2025**, donne une raison précise d’étudier la largeur des graphes. Les travaux **Smetana–Patera, 2016**, puis **Buhr–Smetana, 2018**, donnent un moyen de choisir et de contrôler l’information transmise entre composants. Ces résultats structurent les obligations de preuve détaillées dans les pages suivantes. Ils ne démontrent pas une accélération du multicorps général.

### Ce que « souverain » doit signifier pour le projet

Pouvoir redériver les équations, connaître les hypothèses qui rendent l’algorithme fiable, maîtriser son implémentation et reproduire ses mesures. Le prestige d’une théorie ou sa date de publication ne suffit pas. Les mathématiques fondamentales peuvent fournir le mécanisme de rupture ; l’architecture du logiciel, les données et l’exécution sur machine restent nécessaires pour le réaliser.

Le document joint sur la « stack théorique idéale » et la précédente étude de 58 références ont été pris en compte. Cette nouvelle recherche examine **21 sources primaires ou documentations officielles**, avec lecture ciblée des théorèmes structurants. Elle distingue quatre niveaux : résultat publié ; conséquence algébrique redérivée ici ; hypothèse de transfert ; performance à mesurer.

**Conclusion scientifique actuelle :** une avance généraliste reste à construire. Une cible crédible serait un solveur qui adapte localement la représentation au niveau de précision demandé, conserve les lois mécaniques et sait détecter quand sa réduction devient invalide. Aucune exclusivité mondiale de cette combinaison n’est établie.

---

## 1. Définir mathématiquement le gain recherché

### Le coût à précision physique imposée

Pour un modèle, des charges, une durée et des observables fixés, la quantité à minimiser est le coût total nécessaire pour respecter une erreur ε. Un comptage utile est :

**C = Cpréparation + Σpas ΣNewton (Cassemblage + Crésolution + Ccontrôle) + Csorties.**

Cette décomposition est notre cadre d’évaluation. Diminuer les inconnues peut réduire plusieurs termes simultanément. À l’inverse, un modèle réduit peut perdre tout son avantage si sa construction, sa mise à jour ou l’évaluation de ses forces visite encore tout le modèle complet.

L’erreur doit porter sur les grandeurs utilisées : positions, orientations, efforts de liaison, fréquences, phase, énergie, impulsions ou sensibilités. Ces exigences donnent des problèmes d’approximation différents. Une bonne fréquence propre ne certifie ni les réactions ni le comportement après impact.

### Les trois paramètres cachés

**Dimension effective r(ε).** Combien de variables sont nécessaires pour représenter les réponses admissibles dans une norme physique ? Cette dimension dépend des charges et du domaine d’exploitation, pas seulement du nombre de pièces.

**Largeur des couplages w.** Combien de variables restent simultanément connectées lors d’une élimination ? Deux mécanismes ayant le même nombre de corps et de liaisons peuvent avoir des coûts très différents.

**Constante de stabilité β.** Quelle amplification relie une erreur dans les équations à une erreur dans la solution ? Près d’une singularité, d’un flambement ou d’un changement de contact, un petit résidu peut masquer une grande erreur.

Ces trois paramètres doivent être instrumentés séparément. La recherche de rupture consiste à agir sur r et w tout en empêchant β de s’effondrer artificiellement. Une singularité physique réelle doit cependant rester visible.

### Une limite incontournable

Si l’on exige d’écrire N valeurs indépendantes à chaque pas, le seul coût de sortie est au moins proportionnel à N dans le modèle usuel de calcul. Une promesse de coût inférieur à N doit donc préciser quelles sorties restent réduites, quels prétraitements sont amortis et quelles hypothèses limitent les entrées. C’est une observation de comptage, pas une limitation propre à un concurrent.

Le qualificatif « généraliste » doit porter sur les modèles acceptés et les contrôles disponibles. Il n’impose pas que tous les modèles utilisent la même représentation interne.

---

## 2. Graphes : un résultat fondamental de 2025

### L’énoncé qui peut changer l’architecture

Fürer, Hoppen et Trevisan traitent une matrice m × n sur un corps quelconque. Avec une décomposition arborescente fournie du graphe biparti lignes–colonnes, de largeur k et de taille O(k(m+n)), leur résultat permet notamment rang, déterminant et résolution en **O(k²(m+n)) opérations de corps**. Le calcul explicite de toutes les colonnes d’une forme échelonnée peut coûter davantage. La décomposition du graphe est une entrée, et aucune borne de stabilité en virgule flottante n’est fournie par cet énoncé. [Fast Gaussian Elimination, ESA 2025](https://drops.dagstuhl.de/storage/00lipics/lipics-vol351-esa2025/LIPIcs.ESA.2025.116/LIPIcs.ESA.2025.116.pdf).

### Ce que le multicorps récent exploite déjà

LCABA reformule des éliminations par programmation dynamique non séquentielle. L’analyse de sa première résolution donne O(n+m+m_c²n), où m_c mesure le voisinage couplé maximal de l’algorithme. Un voisinage borné donne un coût linéaire ; les itérations proximales nécessaires à la convergence restent à compter. **m_c n’est pas la largeur arborescente k du théorème précédent.** [Sathya–Carpentier, version auteur associée à TRO 2026](https://www.ajaysathya.com/assets/pdf/LCABA.pdf).

Le manuscrit 2026 sur les opérateurs de Delassus calcule des produits sans construire J M⁻¹ Jᵀ. Il distingue le produit, annoncé linéaire dans son modèle, de l’inverse amorti, dont le pire cas inclut O(n+m²d+m³) opérations. **Divergence non résolue :** l’introduction annonce O(n+m²) mémoire, la section V-C indique O(n+dm²). Aucune garantie mémoire générale n’est retenue ici. Il faut aussi compter la taille des contraintes si leur arité croît. Un inverse amorti exact n’est pas l’inverse du problème non amorti. [Matrix-Free Delassus Operations](https://www.ajaysathya.com/assets/pdf/26_matrix-free_delassus.pdf).

### Programme Vinkulum

Construire un graphe des inconnues et équations avant condensation ; mesurer les tailles des fronts réellement produits. Éliminer les blocs internes à une sous-structure lorsque cela ne crée pas un front plus coûteux. Garder une représentation implicite des directions admissibles : une base globale dense peut annuler le bénéfice de cette structure.

**Verrou de recherche :** trouver une élimination numériquement stable qui exploite la petite largeur sans supposer des pivots favorables. Le rang structurel, le rang exact et le rang numérique doivent être distingués. Une dépendance presque exacte peut imposer un pivotement qui augmente le front.

**Expérience réfutante :** comparer chaînes, boucles locales répétées et boucles globalement couplées à taille égale, puis approcher une singularité. Si les fronts ou les corrections de stabilité croissent fortement, une revendication de coût linéaire doit être retirée pour cette famille.

---

## 3. Topologie et analyse fonctionnelle : éviter les faux degrés de liberté

### Le véritable apport du calcul extérieur discret

Dans le cadre FEEC, le théorème 3.8 d’Arnold–Falk–Winther relie la stabilité du problème mixte de Hodge à un sous-complexe discret et à des projections de cochaînes uniformément bornées. La relation de commutation s’écrit **d π_h = π_h d**. Les estimations dépendent notamment des constantes de Poincaré et des projections. Cette structure permet de préserver au niveau discret les relations qui rendent le problème continu bien posé. Ce n’est pas un théorème de stabilité déjà applicable à une poutre SE(3) non linéaire. [FEEC, 2010, théorèmes 3.8–3.9](https://arxiv.org/pdf/0906.4325).

Mardal et Winther organisent le préconditionnement autour des espaces fonctionnels et de leurs applications de Riesz. Pour les systèmes selles dépendant de paramètres, les normes doivent refléter ces paramètres ; une mise à l’échelle diagonale générique ne remplace pas cette construction. Leur cadre indique comment relier stabilité de l’opérateur et stabilité de la résolution discrète. [Préconditionnement opérateur, 2011](https://web-backend.simula.no/sites/default/files/publications/Simula.SC.570.pdf).

### La cible propre aux poutres et aux liaisons

Vinkulum doit dériver une formulation mixte dans laquelle efforts, moments, déformations et cinématique ont des espaces compatibles. La condition à établir est une borne du type :

**inf_λ sup_v b(v,λ)/(∥v∥V ∥λ∥Q) ≥ β₀ > 0.**

Elle doit être accompagnée de la coercivité sur le noyau, d’une continuité uniforme et d’un traitement explicite des mouvements rigides. L’objectif est que les constantes restent contrôlées lorsque le maillage change et lorsque le rapport cisaillement/flexion augmente, dans un domaine physique précisé. Il faut utiliser un quotient pour les multiplicateurs non uniques.

Ce serait une preuve de robustesse de formulation. Ajouter des inconnues indépendantes ou condenser une matrice ne la fournit pas. En particulier, les difficultés du prototype de poutre mixte à contraste extrême ne seront pas résolues par le seul choix d’un factoriseur plus rapide.

### Pourquoi cela peut accélérer

Une formulation qui ne verrouille pas peut atteindre la précision demandée avec moins d’éléments ; une norme adaptée peut éviter la croissance du nombre d’itérations. Ces deux effets sont distincts et doivent être mesurés séparément.

**Expérience réfutante :** vérifier simultanément convergence de déplacement, efforts et spectre sur une suite de maillages et d’élancements. Si la borne discrète s’effondre ou si apparaissent des modes parasites, le transfert reste invalide. Un essai réussi au bout d’une console ne suffit pas.

---

## 4. Gamblets : résoudre dans les bonnes échelles

### Une construction issue de la complexité de l’information

Owhadi construit des fonctions adaptées à l’opérateur par récupération optimale à partir de mesures hiérarchiques. Pour le cadre elliptique étudié, elles donnent des espaces orthogonaux en énergie entre niveaux, des blocs bien conditionnés et une décroissance permettant la localisation. Les théorèmes 5.5–5.6 contrôlent localisation et résolutions approchées ; la complexité devient quasi linéaire sous les hypothèses de l’article. Ce sont des bases déterministes ; l’interprétation probabiliste ne signifie pas qu’il faille entraîner un réseau. [Gamblets, SIAM Review 2017](https://arxiv.org/pdf/1503.03467).

Owhadi–Zhang adaptent ensuite ces espaces à des systèmes implicites contenant masse et raideur, avec un opérateur de type **4M/Δt² + K**. L’article distingue préparation en O(N log^(3d)N) et propagation en O(N log^(d+1)N) pour son objectif de précision à l’échelle du maillage. Le théorème 3.1 borne l’écart à une intégration au point milieu résolue exactement. Le cadre utilise notamment une ellipticité uniforme et un opérateur linéaire déterminé. [Gamblets temporelles, 2017](https://arxiv.org/pdf/1606.07686).

### Transfert proposé

Tester ces idées sur les blocs flexibles coercifs d’une tangente effective, avant le KKT général. Pour un pas mécanique implicite, une matrice candidate peut combiner M, Δt C et Δt² K. Les coefficients précis dépendent du schéma ; sa symétrie et sa positivité doivent être vérifiées, pas supposées.

Une hiérarchie construite une fois peut servir à plusieurs résolutions si l’opérateur varie modérément. Il faut borner cette variation dans la métrique de l’opérateur de référence. Une rotation rigide peut être traitée par transport exact lorsque la formulation le permet ; une déformation qui change fortement la tangente exige un contrôle différent.

### La difficulté de composition

Projeter d’abord toutes les contraintes peut transformer un opérateur local en opérateur dense. Construire séparément la hiérarchie et la projection peut donc détruire l’hypothèse de localisation. Le problème scientifique intéressant est de conserver simultanément contraintes locales, modes rigides et décomposition en échelles.

**Limites :** les constantes de contraste ne disparaissent pas parce que les coefficients sont peu réguliers. Une tangente proche du flambement ou non symétrique sort du cadre coercif utilisé ici. Des mises à jour fréquentes peuvent rendre la préparation plus chère que les gains.

**Expérience réfutante :** publier coût de construction, coût par résolution et coût de reconstruction quand Δt ou la configuration change, contre une factorisation creuse et un multigrille convenablement réglés. Le gain doit survivre au coût total.

---

## 5. Interfaces optimales : réduire ce qui traverse les composants

### Un critère optimal, pas un nombre de modes arbitraire

Smetana–Patera construisent un opérateur compact T : une donnée sur le bord externe d’un ensemble de composants produit une réponse sur une interface interne, après résolution locale. Les vecteurs singuliers dominants de T définissent un espace optimal au sens de Kolmogorov : **d_r = σ_(r+1)(T)** dans les espaces et normes spécifiés. La proposition 3.8 donne une borne d’erreur de condensation avec une constante géométrique. Le résultat paramétrique du spectral greedy concerne le domaine d’apprentissage considéré et une borne sur le nombre d’interfaces, sans garantie universelle hors de cet ensemble. [Espaces optimaux, 2016](https://doi.org/10.1137/15M1009603).

Les noyaux rigides et les solutions associées aux charges internes doivent être ajoutés à l’espace de transfert. L’exemple elliptique séparable permet de voir la décroissance des composantes qui atteignent l’interface. Cette décroissance ne se transpose pas sans preuve aux ondes propagatives, aux résonances ou aux contraintes de contact appliquées à l’intérieur du composant. [Smetana, développement et estimation, 2018](https://arxiv.org/pdf/1808.02946).

### Rendre la construction abordable

Buhr–Smetana remplacent le calcul complet de l’espace singulier par des résolutions locales avec données aléatoires et enrichissement adaptatif. La proposition 3.2 donne une approximation presque optimale en espérance ; la proposition 3.7 donne un majorant probabiliste de norme. Les métriques interviennent dans les constantes. La probabilité d’échec doit être cumulée sur les contrôles ; un seuil par test ne vaut pas le même seuil pour une longue campagne. [Réduction locale randomisée, 2018](https://arxiv.org/pdf/1706.09179).

### Application propre à Vinkulum

Sur une structure ramifiée, conserver les mouvements rigides exacts de chaque composant et choisir les traces élastiques par le transfert. Ajouter les directions nécessaires aux efforts de liaison, puis enrichir localement quand le résidu indique une perte de précision.

La distinction avec une base de modes propres est importante : une faible fréquence ne mesure pas à elle seule la capacité d’un mode à transmettre une charge donnée entre interfaces.

**Limite souvent oubliée :** l’interface d’une poutre à deux nœuds contient déjà peu de variables. Y appliquer une réduction de port supplémentaire peut ne rien apporter. La cible est une sous-structure de nombreuses poutres, une interface surfacique ou un composant flexible complexe.

**Expérience réfutante :** changer l’emplacement et la direction des charges, déplacer une liaison et exciter une fréquence hors construction. Publier aussi le cas où l’enrichissement revient presque au modèle complet.

---

## 6. Largeurs stables : savoir quand la compression est impossible

### La question de dimension intrinsèque

Cohen, DeVore, Petrova et Wojtaszczyk définissent des largeurs de variétés stables : encodeur et décodeur doivent satisfaire des bornes de Lipschitz. Leur théorème 3.3 relie ces largeurs aux nombres d’entropie du compact à approximer. Cela limite les performances d’une représentation non linéaire stable : on ne peut pas invoquer sa non-linéarité pour ignorer la richesse de la famille des solutions. Le théorème ne fournit ni une base Vinkulum ni sa dimension nécessaire. [Optimal Stable Nonlinear Approximation, 2022](https://arxiv.org/pdf/2009.09907).

### Conséquence pour l’apprentissage et les variétés mécaniques

Pour un encodeur E et un décodeur D, mesurer seulement ∥u−D(E(u))∥ sur les états d’apprentissage ne suffit pas. Il faut aussi étudier la sensibilité des coordonnées et de leur reconstruction. Une petite erreur d’état peut coexister avec un mauvais effort ou un mauvais gradient.

La famille à approximer doit donc inclure les variations de charge, de géométrie et d’état de contact réellement visées. Si elle contient de nombreux événements localisés indépendants, sa dimension effective peut augmenter fortement. Le choix raisonnable peut être plusieurs représentations locales, avec transitions contrôlées, plutôt qu’une seule très petite variété.

### Des bases qui suivent la dynamique

Pagliantini construit une réduction dynamique de systèmes hamiltoniens paramétrés sur une variété de matrices de rang fixé. La base évolue avec la solution et conserve une structure orthosymplectique ; la géométrie du fibré fixe les ambiguïtés de représentation. C’est une alternative à une base globale figée. La dimension complète intervient encore dans l’évolution de la base et l’évaluation du champ ; le cadre concerne des systèmes hamiltoniens et des collections de paramètres. [Dynamical reduced basis methods, 2021](https://arpi.unipi.it/bitstream/11568/1169930/2/s00211-021-01211-w-15.pdf).

### Décision

Pour Vinkulum, une variété adaptative est une hypothèse de recherche, avec quatre contrôles : défaut d’invariance dynamique ; conditionnement de la représentation ; coût de mise à jour ; erreur des observables hors construction. La conservation d’une structure symplectique ne certifie pas à elle seule la précision de phase ou des contraintes.

**Expérience réfutante :** forcer un transfert d’énergie entre modes ou déplacer une zone de flexion importante. Si le rang nécessaire croît presque comme la dimension complète, ou si la mise à jour coûte davantage que la résolution directe, la réduction ne convient pas à cette famille.

---

## 7. Physique fondamentale : que deviennent les variables supprimées ?

### La raideur ne suffit pas à éliminer un mode

Le théorème 1 de Haller–Ponsioen justifie une réduction lente–rapide sous trois conditions : prolongement régulier dans la limite singulière ; existence d’une variété critique ; stabilité asymptotique de la dynamique rapide. Pour un paramètre assez petit, il existe une variété lente attirante et une expansion de la dynamique réduite. Un oscillateur rapide non amorti ne satisfait pas, par sa seule fréquence élevée, la condition d’attraction. L’existence exacte de la variété ne rend pas exacte une expansion tronquée. [Slow–Fast Decomposition, 2017](https://georgehaller.com/reprints/slow_fast.pdf).

### Une élimination exacte peut créer de la mémoire

La dérivation suivante est élémentaire et indépendante. Pour deux blocs linéaires :

**ẋ = A₁₁x + A₁₂y ; ẏ = A₂₁x + A₂₂y.**

La variation des constantes donne :

**y(t) = exp(tA₂₂)y(0) + ∫₀ᵗ exp((t−s)A₂₂)A₂₁x(s) ds.**

En remplaçant y dans la première équation, on obtient une force dépendant de l’histoire de x et de l’état initial supprimé. Si exp(tA₂₂) oscille sans décroître, rien ne justifie de remplacer cette mémoire par une force instantanée. La dimension a été réduite, mais le coût a pu migrer vers l’histoire.

### Une avancée de 2025 délimite la fermeture thermodynamique

Mielke, Peletier et Zimmer partent d’un système hamiltonien fini couplé linéairement à un bain infini linéaire. L’hypothèse 1.2 impose que la compression de l’évolution du bain soit un semi-groupe dissipatif fini ; c’est ce qui autorise une fermeture markovienne. Avec une variable d’énergie supplémentaire, le théorème 4.10 construit une structure GENERIC : énergie conservée et entropie non décroissante dans son cadre. Ce résultat ne dit pas qu’une structure flexible finie quelconque possède cette propriété. [Deriving a GENERIC system, 2025](https://link.springer.com/article/10.1007/s00205-025-02119-7).

### Conséquence d’architecture

Tout composant réduit doit déclarer les états cachés qu’il conserve : variables de mémoire, modes rapides, dissipation et initialisation. Une simple interface effort–position ne suffit pas toujours. L’adjoint doit voir ces états ; les ignorer peut produire le gradient exact d’un mauvais modèle réduit.

**Expérience réfutante :** initialiser deux modèles avec les mêmes variables conservées mais des états internes différents. Si les sorties futures diffèrent sensiblement, une fermeture déterministe sans mémoire dans ces seules variables est inadéquate.

---

## 8. Contact : la régularité est une ressource algorithmique

### Ce que démontre réellement Newton semismooth*

Gfrerer, Mandlmayr, Outrata et Valdman travaillent sur des équations généralisées. Leur théorème 4.4 garantit une convergence locale superlinéaire lorsque l’application est SCD semismooth* et SCD régulière autour de la solution, avec les étapes d’approximation prévues par l’algorithme. L’application mécanique examinée est un contact statique avec frottement de Coulomb. Ce résultat ne garantit ni convergence globale depuis n’importe quel état ni différentiabilité classique d’une trajectoire d’impacts. [SCD semismooth* Newton, 2023](https://arxiv.org/pdf/2205.15129).

### Trois questions qui doivent rester distinctes

**Résoudre la loi.** La complémentarité normale et la loi tangentielle doivent être satisfaites pour le modèle choisi. Une modification convexe ou une régularisation peuvent changer le problème physique.

**Suivre la branche.** Au voisinage d’une solution régulière, une méthode de Newton adaptée peut être très rapide. À une bifurcation ou une transition dégénérée, la même garantie peut disparaître. Le solveur doit exposer cette perte, pas simplement durcir un seuil numérique.

**Différencier l’événement.** Un gradient de la résolution locale à contact fixé ne comprend pas automatiquement la variation de l’instant d’impact. Il faut une dérivation séparée de la trajectoire et de ses transitions. La généralité de cette composition reste une obligation de preuve pour Vinkulum.

### Réduire aussi la mécanique de l’interface

Monjaraz Tec, Gross et Krack associent synthèse modale et frontière sans masse, avec des variantes de MacNeal et Craig–Bampton, pour le contact élastodynamique. Les exemples de l’article montrent une diminution des oscillations parasites et des coûts. Cette voie modifie la représentation discrète de la masse à la frontière ; elle ne se réduit pas à accélérer le solveur de contact existant. [Massless boundary CMS, 2022](https://arxiv.org/pdf/2111.07693).

### Programme et frontière de validité

Conserver une résolution physique des zones susceptibles de contact et réduire l’intérieur. Étudier la stabilité du système mixte ainsi obtenu, notamment lorsque la masse devient singulière. Une métrique fondée sur M⁻¹ ne peut alors être reprise sans adaptation.

**Expérience réfutante :** impact oblique avec frottement, impacts presque simultanés, séparation et retour en contact. Contrôler impulsions, bilan d’énergie, chronologie et réactions. Une réduction qui passe les vibrations libres mais change ces observables doit être enrichie ou refusée.

---

## 9. Composer les briques : une borne que l’on peut dériver

### Un résultat conditionnel simple

Voici une conséquence d’algèbre linéaire redérivée pour cadrer Vinkulum ; ce n’est pas un nouveau théorème revendiqué. Soit un problème linéaire symétrique contraint, à contraintes homogènes de rang fixé. Dans des coordonnées admissibles analytiques Z, poser **H = Zᵀ A Z**, supposé défini positif, et résoudre Hu = f. Z sert ici à la preuve ; il n’est pas nécessaire de le construire en production.

Pour une approximation admissible û, poser r = f−Hû et définir ∥v∥H² = vᵀHv. Alors :

**∥u−û∥H² = rᵀ H⁻¹r.**

Si un opérateur B vérifie **cH⁻¹ ≤ B ≤ CH⁻¹** avec c et C positifs connus, dans l’ordre des formes quadratiques, on obtient :

**(rᵀBr)/C ≤ ∥u−û∥H² ≤ (rᵀBr)/c.**

La preuve est immédiate en remplaçant r par H(u−û), puis en appliquant les inégalités de formes. Pour une observable linéaire ℓ, Cauchy–Schwarz donne aussi **|ℓ(u−û)| ≤ ∥ℓ∥H⁻¹ ∥u−û∥H**.

### Pourquoi cette borne compte

Elle relie le préconditionneur, le résidu et l’approximation réduite dans la même norme. Un enrichissement local peut être décidé par l’erreur de solution visée, et non par un seuil de résidu arbitraire. Mais si c est seulement estimé, le résultat calculé est un indicateur, pas un certificat garanti.

Si l’approximation n’est pas admissible, il faut contrôler la correction de contrainte par une borne inf-sup. Si l’opérateur est indéfini ou non symétrique, la preuve ci-dessus ne s’applique plus telle quelle. Si le résidu est lui-même hyperréduit, son erreur doit être majorée ; on ne peut pas certifier un modèle uniquement avec les équations qu’il a choisi de conserver.

### Les quatre obligations de composition

1. **Approximation :** relier erreur de transfert, erreur intérieure et erreur de sortie, avec constantes d’assemblage explicites.

2. **Stabilité :** maintenir les bornes nécessaires malgré changement d’unités, élancement, configuration et contraintes redondantes.

3. **Temps :** propager ces erreurs sur l’horizon simulé ; la borne d’un solveur statique ne devient pas une borne de trajectoire par addition informelle.

4. **Sensibilités :** contrôler dérivées de bases, changements de représentation et événements. Une erreur d’état petite ne borne pas, seule, l’erreur de gradient.

**La contribution scientifique à rechercher est cette chaîne vérifiable.** Partager SO(3) ou une notation énergétique ne prouve pas que les briques la satisfont ensemble.

---

## 10. Ce que les concurrents ont déjà

| Référence | Capacité établie ou documentée | Conséquence pour Vinkulum |
|---|---|---|
| Exudyn 1.11.0 | FFRF, Hurty–Craig–Bampton, modes creux ; arbres en coordonnées minimales avec connexions de fermeture | Une comparaison doit inclure ses représentations efficaces, pas seulement un modèle redondant |
| MBDyn | Poutres géométriquement exactes et synthèse modale de composants | Réduction et géométrie ne sont pas des nouveautés suffisantes |
| Simpack | Réduction linéaire et non linéaire, SIMBEAM, couplage Abaqus ; amélioration de sélection automatique des modes annoncée en 2025x | Ni réduction non linéaire ni automatisation ne peuvent être revendiquées comme terrain vierge |

Ces capacités proviennent des [guides FFRF/CMS Exudyn](https://exudyn.readthedocs.io/en/v1.11.0/docs/RST/ModelOrderReductionAndComponentModeSynthesis.html), de la [documentation KinematicTree](https://exudyn.readthedocs.io/en/v1.11.0/docs/RST/items/ObjectKinematicTree.html), de la [FAQ MBDyn](https://www.mbdyn.org/Documentation/FAQ.html), des [modules flexibles Simpack](https://www.3ds.com/products/simulia/simpack/flexible-body) et des [nouveautés Simpack 2025x](https://3dswym.3dexperience.3ds.com/wiki/simulia-community/simpack-release-2025x-news_6CSq0Y9uTB-kTqrJ9h-qzQ). Elles décrivent un périmètre ; elles ne classent pas les performances et ne rendent pas publics tous les algorithmes internes.

### Où placer une différenciation défendable

Le pari proposé est l’association d’une représentation locale adaptable, de contraintes traitées sans densification inutile et d’un contrôle d’erreur utile aux efforts et sensibilités. L’objectif serait de réduire l’expertise manuelle nécessaire pour choisir une base tout en gagnant du temps total sur des modèles connectés difficiles.

Cette combinaison est une hypothèse de positionnement. La documentation consultée ne suffit ni à établir son absence chez les concurrents ni à démontrer la supériorité d’une future réalisation.

### État propre à ne pas confondre avec la cible

Vinkulum 0.9.0 a livré des opérateurs adjoints creux dans son domaine lisse holonome. Les améliorations mesurées contre ses propres versions précédentes ne sont pas un classement d’adjoints contre Exudyn, MBDyn ou Simpack. Le prototype de modes sélectifs reste expérimental ; aucune performance de ses sondages n’est promue dans cette étude.

Les confrontations archivées du dépôt restent les preuves disponibles pour les cas déjà exécutés. **Cette recherche n’ajoute aucun nouveau benchmark externe.** Simpack n’a pas été exécuté dans ce travail. Le niveau atteint par un concurrent doit être apprécié sur un modèle physique et un protocole comparables, avec ses options pertinentes.

---

## 11. Programme de recherche et critères d’abandon

### Priorité 1 — Structure stable des contraintes

Livrable : dérivation des éliminations par sous-structure, coût en fonction des fronts et politique de rang numérique. Cas : chaîne, réseau ramifié, boucles locales, boucles globales, singularité approchée. Témoins : KKT creux actuel, représentation minimale lorsque disponible, Exudyn KinematicTree pour les cas compatibles. Mesurer remplissage, mémoire, temps, contraintes et réactions.

**Critère scientifique :** le coût doit suivre la structure annoncée sans sacrifier la précision. Une explosion des fronts ou des corrections invalide la revendication sur la famille concernée.

### Priorité 2 — Réduction locale avec contrôle d’erreur

Livrable : transfert local, noyaux exacts, estimateur et enrichissement. Comparer base modale, espace optimal et construction randomisée sur des composants réellement volumineux. Déplacer les charges et les interfaces ; inclure des paramètres hors construction et des cas proches de résonances.

**Critère scientifique :** publier erreur effective et majorant, puis le coût complet. Pour Q utilisations, une préparation supplémentaire C₀ n’est amortie que si Q(Ccomplet−Créduit) > C₀. Si le dénominateur est nul ou négatif, il n’existe pas d’amortissement favorable. Le contrôle d’erreur doit être inclus dans Créduit.

### Priorité 3 — Robustesse de formulation

Livrable : choix des espaces mixtes, preuve ou borne vérifiable de stabilité et préconditionnement cohérent. La première étape est linéaire ; l’extension à grandes rotations doit conserver les propriétés déjà établies. Une convergence correcte avec un maillage excessif ne répond pas au problème de formulation.

**Critère scientifique :** distinguer mauvais conditionnement physique, défaut discret et perte due à l’arithmétique. La perte de stabilité près d’un flambement réel ne doit pas être masquée pour faire réussir le test.

### Priorité 4 — Dynamique réduite, mémoire et contact

À engager après les trois précédentes : évolution de base, fermeture avec mémoire, traitement des frontières de contact et sensibilités des transitions. Chaque extension introduit ses propres hypothèses. Le saut direct à un solveur universel réduit, non lisse et différentiable n’est pas justifié par les sources.

### Exemple concret de décision attendue

Pour un mécanisme déployable formé de nombreux segments flexibles, le solveur conserverait les mouvements rigides, réduirait les déformations internes et ne transmettrait entre segments que les directions nécessaires à la précision des efforts. Lorsqu’une butée entre en contact, il enrichirait la région touchée et contrôlerait la trajectoire résultante. **Ce comportement est la cible**, pas une capacité disponible ni une accélération chiffrée promise.

La prochaine intégration doit être choisie par une dérivation et une expérience discriminante. C’est ainsi que cette recherche peut transformer l’architecture plutôt qu’allonger sa bibliographie.

---

## Sources et accès

Les sources suivantes ont été consultées le 7 septembre 2026. Les textes mathématiques ont fait l’objet de lectures ciblées des définitions, résultats et hypothèses indiqués dans le rapport ; ce document ne prétend pas à un audit intégral de toutes leurs preuves. Les versions auteur peuvent différer des versions finales. Les sources logicielles désignent des documentations, jamais du code d’implémentation concurrent.

**Fürer, Hoppen, Trevisan (2025).** Fast Gaussian Elimination for Low Treewidth Matrices. ESA 2025, 116:1–15. Texte primaire, théorème 1 et conséquences. [Article Dagstuhl, DOI 10.4230/LIPIcs.ESA.2025.116](https://drops.dagstuhl.de/storage/00lipics/lipics-vol351-esa2025/LIPIcs.ESA.2025.116/LIPIcs.ESA.2025.116.pdf).

**Sathya, Carpentier (2026 ; préprint 2025).** Constrained Articulated Body Dynamics Algorithms for Closed-Loop Mechanisms. IEEE Transactions on Robotics, DOI 10.1109/TRO.2026.3651683. PDF auteur 20 pages, analyses IV-D et V-D ; préprint antérieur recoupé. [Manuscrit auteur](https://www.ajaysathya.com/assets/pdf/LCABA.pdf).

**Sathya, Montaut, de Mont-Marin, Carpentier (2026).** Matrix-Free Delassus Operations: Scalable and Memory-Efficient Algorithms. Manuscrit auteur, 8 pages, encore présenté avant fin d’évaluation dans le fichier consulté. Le statut d’une version finale distincte n’est pas utilisé pour les conclusions. [Texte auteur](https://www.ajaysathya.com/assets/pdf/26_matrix-free_delassus.pdf).

**Arnold, Falk, Winther (2010).** Finite element exterior calculus: from Hodge theory to numerical stability. Bulletin AMS 47, 281–354. Prépublication complète, théorèmes 3.8–3.9. [Texte primaire](https://arxiv.org/pdf/0906.4325).

**Mardal, Winther (2011).** Preconditioning discretizations of systems of partial differential equations. Numerical Linear Algebra with Applications 18, 1–40, DOI 10.1002/nla.716. Manuscrit institutionnel ; date de publication recoupée, distincte du pied de page générique du manuscrit. [Texte Simula](https://web-backend.simula.no/sites/default/files/publications/Simula.SC.570.pdf).

**Owhadi (2017).** Multigrid with Rough Coefficients and Multiresolution Operator Decomposition from Hierarchical Information Games. SIAM Review 59, 99–149, DOI 10.1137/15M1013894. Texte primaire, localisation et complexité. [Prépublication complète](https://arxiv.org/pdf/1503.03467).

**Owhadi, Zhang (2017).** Gamblets for opening the complexity-bottleneck of implicit schemes for hyperbolic and parabolic ODEs/PDEs with rough coefficients. Journal of Computational Physics 347, 99–128. Texte primaire, opérateur 2.1 et théorème 3.1. [Prépublication complète](https://arxiv.org/pdf/1606.07686).

**Smetana, Patera (2016).** Optimal Local Approximation Spaces for Component-Based Static Condensation Procedures. SIAM Journal on Scientific Computing 38, A3318–A3356. Métadonnées éditeur et texte auteur déposé sur ResearchGate ; propositions 3.5 et 3.8. [Publication, DOI 10.1137/15M1009603](https://doi.org/10.1137/15M1009603).

**Smetana (2018).** Static condensation optimal port/interface reduction and error estimation for structural health monitoring. Prépublication arXiv 1808.02946 ; section 3 et proposition 3.4. [Texte primaire](https://arxiv.org/pdf/1808.02946).

**Buhr, Smetana (2018).** Randomized Local Model Order Reduction. SIAM Journal on Scientific Computing 40, A2120–A2151, DOI 10.1137/17M1138480. Algorithme 1, propositions 3.2 et 3.7–3.8. [Prépublication complète](https://arxiv.org/pdf/1706.09179).

**Cohen, DeVore, Petrova, Wojtaszczyk (2022).** Optimal Stable Nonlinear Approximation. Foundations of Computational Mathematics 22, DOI 10.1007/s10208-021-09494-z. Définitions des largeurs stables et théorème 3.3. [Prépublication complète](https://arxiv.org/pdf/2009.09907).

**Pagliantini (2021).** Dynamical reduced basis methods for Hamiltonian systems. Numerische Mathematik 148, 409–448, DOI 10.1007/s00211-021-01211-w. Géométrie des espaces de rang fixé et projections. [Texte institutionnel](https://arpi.unipi.it/bitstream/11568/1169930/2/s00211-021-01211-w-15.pdf).

**Haller, Ponsioen (2017).** Exact model reduction by a slow–fast decomposition of nonlinear mechanical systems. Nonlinear Dynamics 90, 617–647, DOI 10.1007/s11071-017-3685-9. Hypothèses A1–A3 et théorème 1, prépublication et version publiée recoupées. [Article auteur](https://georgehaller.com/reprints/slow_fast.pdf).

**Mielke, Peletier, Zimmer (22 septembre 2025).** Deriving a GENERIC system from a Hamiltonian system. Archive for Rational Mechanics and Analysis 249, 62. Intégral ouvert ; hypothèse 1.2 et théorème 4.10. [Article éditeur](https://link.springer.com/article/10.1007/s00205-025-02119-7).

**Gfrerer, Mandlmayr, Outrata, Valdman (2023 ; en ligne 2022).** On the SCD semismooth* Newton method for generalized equations with application to a class of static contact problems with Coulomb friction. Computational Optimization and Applications 86, 1159–1191, DOI 10.1007/s10589-022-00429-0. Théorème 4.4. [Prépublication complète](https://arxiv.org/pdf/2205.15129).

**Monjaraz Tec, Gross, Krack (2022).** A massless boundary component mode synthesis method for elastodynamic contact problems. Computers & Structures 260, 106698, DOI 10.1016/j.compstruc.2021.106698. Formulation et exemples, gains propres à l’article. [Prépublication complète](https://arxiv.org/pdf/2111.07693).

**Exudyn, version 1.11.0.** Model order reduction and component mode synthesis. Guide théorique FFRF, HCB et calcul modal. [Documentation officielle](https://exudyn.readthedocs.io/en/v1.11.0/docs/RST/ModelOrderReductionAndComponentModeSynthesis.html).

**Exudyn, version 1.11.0.** ObjectKinematicTree. Description officielle accessible par index de recherche ; l’ouverture directe de la page a échoué. [Documentation officielle](https://exudyn.readthedocs.io/en/v1.11.0/docs/RST/items/ObjectKinematicTree.html).

**MBDyn.** FAQ, éléments flexibles et synthèse modale. [Documentation officielle](https://www.mbdyn.org/Documentation/FAQ.html).

**Dassault Systèmes / SIMULIA.** Flexible Body Simulation Modules. Capacités déclarées de Simpack, sans accès aux algorithmes internes. [Documentation produit](https://www.3ds.com/products/simulia/simpack/flexible-body).

**Dassault Systèmes / SIMULIA (2025x).** Simpack Release 2025x News. Annonce d’améliorations de sélection automatique des modes. [Annonce officielle](https://3dswym.3dexperience.3ds.com/wiki/simulia-community/simpack-release-2025x-news_6CSq0Y9uTB-kTqrJ9h-qzQ).

### Limites et portée de la livraison

L’étude porte sur des transferts mathématiques possibles vers Vinkulum. Elle ne revendique ni exhaustivité bibliographique mondiale ni validation de toutes les compositions proposées. Les inconnues décisives sont désormais des preuves de stabilité, des extensions de domaine et des expériences à exécuter. Aucun résultat commercial propriétaire n’est inféré d’une absence de publication.
