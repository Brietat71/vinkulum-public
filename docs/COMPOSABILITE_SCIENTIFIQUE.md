# Composabilité scientifique : lecture critique et décisions pour Vinkulum

Relecture du 7 septembre 2026 du document transmis par l'utilisateur,
*Stack théorique idéale en dynamique multicorps : état de l'art mathématique
brique par brique et composabilité*. Les corrections ci-dessous s'appuient
sur des articles primaires ; elles ne constituent pas un nouvel inventaire
exhaustif des solveurs mondiaux. Aucun code source concurrent n'a été consulté.

**Décision : conserver la vocation généraliste, avec plusieurs méthodes
dont les domaines et les interfaces sont vérifiés.** Le document identifie
utilement les objectifs distincts de précision structurelle, d'optimisation
et de contact massif. En déduire un choix d'architecture irréversible serait
excessif. Le partage de SO(3) ou SE(3) ne suffit pas, à lui seul, à prouver
la conservation, la stabilité ou la justesse des gradients d'un assemblage.

## Corrections qui changent les choix de développement

| Affirmation du document | Lecture retenue et conséquence pour Vinkulum |
|---|---|
| Craig–Bampton et grandes rotations seraient incompatibles | Dans un repère flottant, la rotation rigide peut être grande tandis que la déformation locale reste petite. Distinguer amplitude du mouvement et amplitude de déformation ; conserver cette voie de réduction. [1] |
| L'ANCF manquerait intrinsèquement d'objectivité | Cette exclusion générale est incorrecte : une formulation ANCF cohérente peut donner une déformation nulle sous mouvement rigide arbitraire. Le verrouillage dépend de l'élément et de son intégration. Comparer des éléments définis, pas des étiquettes. [2] |
| Les dérivées RNEA/ABA seraient toutes O(n) | Distinguer une direction dérivée, une matrice complète et la profondeur du graphe. L'analyse de Singh et al. donne O(Nd) pour les dérivées de dynamique inverse, avec d la profondeur de l'arbre. Une chaîne peut donc coûter O(N²). Un opérateur appliqué à un vecteur reste une cible distincte. [3] |
| Contact dur et gradients corrects s'excluraient | Des sensibilités locales existent pour des transitions régulières : la saltation inclut la variation de l'instant d'impact. Il faut des gardes et resets différentiables, une traversée transverse et une séquence d'événements localement stable. Cela ne résout pas les impacts rasants, simultanés ou les bifurcations de contact. [4] |
| SSM et ECSW seraient composables par simple géométrie commune | Une SSM exige notamment des conditions spectrales et de non-résonance ; sa propriété d'invariance n'est pas une preuve automatique de conservation de l'énergie du modèle hyper-réduit. Vérifier séparément variété, projection, pondérations et intégration temporelle. [5] |
| Ajouter GGL donnerait un intégrateur symplectique | Le principe variationnel GGL de 2023 permet de construire un schéma spécifique. La simple stabilisation GGL d'un autre intégrateur n'hérite pas de cette preuve. Notre option GGL doit garder ses propres contrôles. [6] |
| La géométrie de Lie supprimerait toutes les singularités | Elle évite une carte globale d'angles, mais le logarithme principal garde une coupure à π. Les incréments, transports et changements de branche doivent être définis. Vinkulum teste explicitement cette limite pour les dérivées des poutres. |

Deux autres distinctions sont nécessaires. Un algorithme conservatif et un
algorithme dissipatif ne peuvent préserver et diminuer la même énergie au
même pas sur le même problème ; ils peuvent pourtant coexister dans un
moteur. De même, plusieurs représentations de rotations peuvent communiquer
si les conversions, espaces tangents et efforts conjugués sont cohérents.
Les incompatibilités concernent les propriétés revendiquées d'un calcul,
pas nécessairement l'existence de plusieurs méthodes dans un logiciel.

Les affirmations « aucun code au monde » et « SOTA incontesté » ne sont
pas retenues comme conclusions vérifiées. L'absence d'un exemple dans une
revue ne prouve pas l'absence mondiale d'une capacité. Les coûts GPU, les
lois de contact et les précisions industrielles exigent des confrontations
exécutées, à problème physique et erreur comparables.

## Ce que cela apporte au programme existant

| Brique | État vérifié dans Vinkulum | Prochaine preuve utile |
|---|---|---|
| Intégration géométrique | Option σ de 2025 livrée en 0.8.0 ; gain ciblé et réactions parfois dégradées | Choisir sur plusieurs observables à erreur commune ; conserver le défaut tant qu'aucun gain général n'est démontré |
| Flexibilité non linéaire | Deux énergies de poutre intégrées au noyau ; prototype mixte 2026 autonome | Efforts et tangente exacts moins coûteux ; puis condensation stable, inertie et dynamique du prototype |
| Dérivées et résolution | Duaux, adjoints et assemblages locaux ; globalisation statique livrée en 0.8.1 | Mesurer produits directs/transposés, résidus adjoints et passage à l'échelle sur chaînes, arbres et boucles |
| Réduction | La réduction non linéaire générale n'est pas démontrée | Commencer par un domaine lisse défini, publier coût hors ligne, erreur hors charges de construction et repli vers le modèle complet |
| Contact | Plusieurs traitements présents, sans garantie globale de gradient aux transitions | Loi normale/tangentielle couplée, bilan de dissipation, différentiation implicite locale et saltation sous hypothèses explicites |
| Couplages | Aucun certificat universel de composabilité | Contrôler travail aux interfaces, inertie transportée et erreur de couplage avant de multiplier les sous-solveurs |

Le chantier de tangente de poutre applique ce programme à une dépense
mesurée du noyau. Le travail virtuel donne les moments analytiques ; un
passage de dérivation automatique d'ordre un fournit leur tangente. Le
contre-calcul reste la dérivation emboîtée de l'énergie avec rotations
successives. Ce transfert repose sur la mécanique variationnelle et les
transports de Lie ; il ne transforme pas notre interpolation en l'élément
SE(3) précis de Sonneville–Cardona–Brüls [7]. Une simplification analytique
utile n'est pas, à elle seule, une nouveauté scientifique.

## Conditions à partager entre les briques

Ces conditions sont une proposition d'architecture et de validation,
pas un certificat déjà obtenu pour tout le noyau.

1. **Cinématique et travail.** Définir pose, perturbation gauche ou droite,
   vitesses, unités, repères et efforts conjugués. Une conversion doit
   préserver la puissance ; une rotation rigide superposée doit transporter
   efforts et tangentes correctement.
2. **Équations et dérivées.** Fournir résidu, action de sa dérivée et action
   transposée cohérents avec les mêmes états internes. Préciser si le
   gradient porte sur le problème convergé, le pas discret ou le modèle
   continu ; contrôler les erreurs de résolution associées.
3. **Énergie et événements.** Séparer énergie physique, travail externe,
   dissipation constitutive et dissipation numérique. Documenter les
   conditions où un événement possède une dérivée, plusieurs dérivées
   directionnelles ou aucune sensibilité classique utilisable.
4. **Erreur et coût.** Comparer à tolérance physique commune, compter
   préparation, réduction, résolution et sorties, puis conserver aussi les
   cas défavorables et les refus. La compatibilité doit être testée sur les
   assemblages, au-delà des tests de chaque élément isolé.

L'avantage recherché est une composition plus économique et plus fiable
de méthodes appropriées à chaque problème. Il reste à mesurer face aux
concurrents : la [confrontation statique 0.8.1](CONFRONTATION_STATIQUE_MBDYN_0.8.1.md)
laisse encore MBDyn 11,4 fois plus rapide que l'option intégrée sur Princeton
au seuil choisi. Aucune comparaison exécutée à Simpack n'est disponible.

## Sources primaires consultées

1. Zwölfer et Gerstmayr, *The nodal-based floating frame of reference
   formulation with modal reduction*, Acta Mechanica, 2021.
   [Article, introduction et réduction](https://link.springer.com/article/10.1007/s00707-020-02886-2).
2. *Absolute nodal coordinate formulation – Multilevel finite element
   framework for the nonlinear multi-scale multibody dynamic analysis of
   composite structures*, Computers & Structures, 2023.
   [Introduction consultable](https://www.sciencedirect.com/science/article/pii/S0045794923002171).
3. Singh, Russell et Wensing, *Efficient Analytical Derivatives of Rigid-Body
   Dynamics using Spatial Vector Algebra*, IEEE RA-L, 2022.
   [Texte intégral, introduction et complexité](https://arxiv.org/pdf/2105.05102).
   Le texte distingue explicitement la complexité de l'algorithme présenté
   en 2018 et celle de son implémentation ; notre lecture reste celle des
   équations et analyses publiées.
4. Kong et al., *Saltation Matrices: The Essential Tool for Linearizing
   Hybrid Dynamical Systems*, prépublication 2023, publication 2024.
   [Texte intégral, §III et V](https://arxiv.org/html/2306.06862v2).
5. Haller et Ponsioen, *Nonlinear normal modes and spectral submanifolds:
   existence, uniqueness and use in model reduction*, 2016.
   [Article et conditions d'existence](https://www.georgehaller.com/reprints/NNM.pdf).
6. *The GGL variational principle for constrained mechanical systems*, 2023.
   [Texte intégral, §3–4](https://link.springer.com/article/10.1007/s11044-023-09889-6).
7. Sonneville, Cardona et Brüls, *Geometrically exact beam finite element
   formulated on the special Euclidean group SE(3)*, CMAME, 2014.
   [Postprint institutionnel](https://orbi.uliege.be/bitstream/2268/159471/1/pa_SonnevilleCardonaBruls2014_GeometricallyExactBeamFEOnSE3.pdf).
