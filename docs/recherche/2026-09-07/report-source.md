# Vinkulum : les avancées scientifiques capables de faire la différence

## Recherche et décisions d’architecture · 2001–2026

**7 septembre 2026 · Étude en français · Vinkulum généraliste**

Destinataire : conception et développement de Vinkulum. Périmètre : travaux publiquement accessibles de 2001 au 7 septembre 2026, confrontés aux capacités documentées des principaux solveurs et à l’état mesuré du projet. État de référence : paquet **0.7.2**, documentation et prototype au commit **723112b**.

### Réponse directe

**Les avancées étudiées justifient une ambition forte, mais la supériorité de Vinkulum reste à construire et à mesurer.** Le projet possède déjà de la géométrie des rotations, de la différentiation automatique et des adjoints. Les limites actuelles de la flexibilité, des contacts, du passage à l’échelle et des couplages empêchent d’en déduire une avance générale. [État du noyau Vinkulum](https://github.com/Brietat71/vinkulum/blob/723112bd51fc71b81d2757bdadc854529b055831/README.md)

**L’avantage le plus crédible à rechercher est une simulation précise de grands mécanismes connectés, rigides et flexibles, accompagnée de sensibilités fiables, à faible coût total.** C’est notre hypothèse stratégique. Elle concerne aussi bien une transmission, une machine, une éolienne ou un mécanisme aéronautique qu’un robot. Sa valeur dépend de la combinaison des méthodes et de leur validation sur plusieurs familles de problèmes.

| Priorité proposée | Apport recherché | Situation réelle |
|---|---|---|
| Résolution exploitant le graphe et dérivation implicite | Réduire les inconnues, les factorisations et le coût des gradients | Méthodes publiées ; transfert complet à démontrer |
| Flexibilité d’ordre élevé, réduction et hyperréduction structurées | Conserver la précision avec moins de travail global | Un prototype local existe ; coût et conditionnement encore défavorables |
| Contact couplé, trajectoires de collision contrôlées, dérivées d’événements | Fiabiliser les cas difficiles et l’optimisation | Plusieurs verrous restent ouverts dans le noyau |
| Intégration géométrique récente et contrôle des erreurs | Améliorer précision temporelle et couplages | Expériences proches du noyau, sans gain acquis |

Le critère de réussite est le **coût pour obtenir une réponse physique vérifiée**, puis pour identifier ou optimiser le système. Ces priorités sont une synthèse des preuves détaillées ci-après. Elles ne constituent pas un classement de performances ni une promesse de domination générale.

---

## 1. Se mesurer aux meilleures références actuelles

Il n’existe pas de classement mondial unique applicable à tous les usages multicorps. Les familles suivantes couvrent les références industrielles, la mécanique flexible, la robotique de contact et le calcul massif. Les capacités ci-dessous proviennent des éditeurs ou des équipes : **elles ne constituent pas une comparaison chronométrée indépendante**.

| Références | Capacités documentées utiles à la comparaison | Conséquence pour Vinkulum |
|---|---|---|
| Simpack | Analyses temporelles, fréquentielles et temps réel ; SIMBEAM, flexibilité et réduction non linéaires, FlexContact. [Modules](https://www.3ds.com/products/simulia/simpack/core) ; [flexibilité](https://www.3ds.com/products/simulia/simpack/flexible-body) | Une poutre non linéaire ou une réduction seule ne crée pas d’avance. |
| Adams, RecurDyn, Simcenter | Chaînes de simulation industrielles ; RecurDyn décrit notamment MFBD, contact et couplages MPS/DEM ; Simcenter intègre le mouvement rigide et flexible à la CAO/CAE. [Adams](https://nexus.hexagon.com/home/product/adams/) ; [RecurDyn](https://support.functionbay.com/en/page/single/232/recurdyn-all-products) ; [Simcenter](https://www.siemens.com/en-us/products/simcenter/simulation-test/motion-simulation/) | Évaluer également assemblage, échanges de modèles et couverture métier. |
| MBDyn | Rigide, flexible, poutres, coques, réduction modale et couplages aérodynamiques, hydrauliques, électriques. [MBDyn](https://www.mbdyn.org/) | Témoin généraliste accessible, à dépasser sur des cas plus variés. |
| Exudyn, Chrono | Exudyn associe Python et C++ pour le multicorps flexible ; Chrono poursuit une couverture multiphysique et parallèle, avec une version 10.0 annoncée en mars 2026. [Exudyn](https://exudyn.readthedocs.io/en/latest/docs/RST/Exudyn.html) ; [Chrono](https://projectchrono.org/news/) | Comparer à des architectures modernes proches de nos ambitions. |
| Drake | Contact hydroélastique distribué ; le modèle décrit n’intègre pas un état complet de déformation ni les ondes élastiques rapides. [Guide Drake](https://drake.mit.edu/doxygen_cxx/group__hydroelastic__user__guide.html) | Comparer des lois physiques équivalentes avant de comparer les temps. |
| MuJoCo, MJX, MJWarp, Newton | Calcul CPU/GPU et simulation en lots ; fonctionnalités et différentiabilité varient selon l’implémentation ou le solveur choisi. [MJX](https://mujoco.readthedocs.io/en/stable/mjx.html) ; [MJWarp](https://mujoco.readthedocs.io/en/latest/mjwarp/) ; [Newton](https://newton-physics.github.io/newton/stable/solvers/index.html) | Séparer latence d’un mécanisme, débit des lots et calcul des gradients. |

### Un indice stratégique, avec une limite claire

La documentation MJWarp distingue explicitement débit et latence : un pas d’un seul environnement peut être plus lent sur GPU. Elle indique aussi que les grands mécanismes connectés restent un axe d’optimisation. Cela soutient notre choix d’un banc consacré au mécanisme connecté ; cela ne prouve pas que Vinkulum puisse y gagner. Il ne faut pas transformer la mention d’environ 60 degrés de liberté en limite universelle de capacité. [Documentation MJWarp, consultée le 7 septembre 2026](https://mujoco.readthedocs.io/en/latest/mjwarp/)

Les méthodes internes des produits propriétaires sont partiellement publiques. Une capacité non documentée ne doit pas être déclarée absente. **Aucun essai Simpack, Adams, RecurDyn, Simcenter, Exudyn, Chrono, Drake ou MuJoCo n’a été exécuté pour ce rapport.**

---

## 2. Le point de départ mesuré et le problème à gagner

### Ce que les mesures établissent

La confrontation de la roue Vinkulum 0.7.2 à MBDyn conserve 326 exécutions et compare les meilleurs réglages admissibles parmi les grilles essayées. Sur le mécanisme plan, les temps médians sont **0,269 contre 0,317 s** au seuil de 100 µm ; sur le spatial, **0,749 contre 0,735 s**. MBDyn reste plus rapide sur Andrews. Sur Princeton, Vinkulum intégré prend **0,387 s contre 0,0186 s**, soit un rapport de **20,8**, au seuil de 10 µm au bout. Les démarrages, runtimes et sorties diffèrent : ce sont des coûts d’usage du protocole, avec erreurs estimées, sans classement isolé des noyaux. [Confrontation exécutée](https://github.com/Brietat71/vinkulum/blob/723112bd51fc71b81d2757bdadc854529b055831/docs/CONFRONTATION_MBDYN_0.7.2.md)

L’effort prioritaire doit donc traiter un coût réel de formulation et de résolution. La présence de méthodes modernes dans le dépôt ne suffit pas à compenser ce retard.

### Exemple concret : une machine flexible en boucle fermée

Considérons une machine dont plusieurs branches flexibles portent un outil, avec une transmission, des appuis frottants et des actionneurs. Le concepteur souhaite modifier sections, matériaux, amortissements et paramètres de commande pour réduire la vibration de l’outil et les pics d’effort.

Le calcul utile comporte quatre résultats : la trajectoire et la phase de l’outil ; les efforts et vibrations dans les branches ; les transitions d’adhérence et de glissement ; la variation de l’objectif lorsque les paramètres changent. Une bonne position terminale peut masquer une mauvaise courbure, une réaction erronée ou un gradient incorrect.

**Expérience proposée, non exécutée :** faire croître indépendamment le nombre de branches, les boucles, les modes flexibles, les contacts et les paramètres. Fixer des seuils sur les observables avant de choisir les réglages. Mesurer d’abord une simulation, puis une optimisation complète, en comptant chaque rejet et chaque recalcul.

### Une architecture généraliste

Les mêmes opérations doivent servir une chaîne cinématique, un rotor flexible et un système couplé à un réseau électrique ou hydraulique : éliminer localement ce qui peut l’être, conserver les interfaces nécessaires, résoudre les contraintes, contrôler les échanges d’énergie, calculer les sensibilités valides.

Notre proposition n’est pas de spécialiser le produit sur cet exemple. Elle consiste à utiliser un cas exigeant comme point d’entrée, puis à imposer des contre-épreuves planes, spatiales, flexibles, non lisses et multiphysiques. Le résultat recherché est un **ensemble de compromis mesurés entre erreur, temps, mémoire et robustesse**, accompagné de limites explicites.

---

## 3. Le levier central : exploiter la structure mécanique

### De 2001 à 2026 : ce qui mérite réellement d’être repris

Les avancées utiles forment une continuité. Les travaux de réduction mécanique préservant la structure et les méthodes Newton–Krylov étaient déjà établis au début des années 2000. Les développements des années 2010 ont renforcé géométrie, algèbre creuse et réduction non linéaire. Les travaux récents apportent des formulations et garanties plus ciblées. La date d’un article ne suffit donc pas à fixer sa priorité. [Lall et al., 2003](https://doi.org/10.1016/S0167-2789(03)00227-6) ; [Knoll et Keyes, 2004](https://doi.org/10.1016/j.jcp.2003.08.010)

### Réduire le travail global avant de changer de processeur

Un système contraint couple inertie, accélérations et réactions. Les inconnues internes d’un élément flexible, les branches d’un arbre cinématique et les interfaces de boucles n’ont pas le même rôle. Notre architecture proposée est :

**Éléments et lois locales → éliminations locales → interfaces et contraintes → résolution globale → contrôles physiques → sensibilités.**

La résolution proximale et creuse de Carpentier, Budhiraja et Mansard exploite le graphe cinématique, y compris avec contraintes redondantes ; elle est déjà mise en œuvre dans Pinocchio. Le lien entre dynamique contrainte et régulation quadratique optimale fournit aussi des éliminations structurées. Les complexités linéaires de certaines variantes ne couvrent toutefois pas les boucles internes générales. [Carpentier et al., 2021](https://www.roboticsproceedings.org/rss17/p017.pdf) ; [Sathya et al., 2024](https://doi.org/10.1109/TRO.2023.3335652)

**Transfert proposé :** réutiliser la structure symbolique tant que la topologie reste valide, condenser les champs internes, traiter les boucles aux interfaces et conserver un repli général. Instrumenter le nombre de factorisations, le remplissage, les produits matrice-vecteur et les réassemblages. Une optimisation pour les arbres doit annoncer son domaine et rester correcte en présence de boucles.

### Le rang et la métrique font partie du problème

LSMR permet de résoudre des moindres carrés par produits avec une matrice et sa transposée. MINRES-QLP traite des systèmes symétriques indéfinis ou singuliers. Le second exige la symétrie effective de l’opérateur ; aucun de ces noms ne justifie une symétrisation arbitraire des équations. Un préconditionnement peut changer la norme minimisée. [Fong et Saunders, 2011](https://stanford.edu/group/SOL/software/lsmr/LSMR-SISC-2011.pdf) ; [Choi et al., 2011](https://doi.org/10.1137/100787921)

Il faut distinguer l’accélération mécaniquement déterminée et les réactions individuelles éventuellement indéterminées. Notre choix doit rendre explicites unités, pondérations et convention de réaction. Une petite norme d’un résidu préconditionné ne suffit pas : les équations physiques originales et les contraintes doivent être contrôlées séparément.

**Décision : priorité structurante.** Elle sert simultanément les grandes assemblées, la robustesse et les adjoints. Le gain devra survivre à la duplication et à la permutation des contraintes, aux changements d’unités et aux rapports de raideurs élevés.

---

## 4. Géométrie et intégration : une amélioration proche du noyau

### Ce qui est déjà un socle partagé

Intégrer les rotations sur leur groupe évite de traiter les orientations comme des vecteurs ordinaires. Les schémas generalized-α sur groupes de Lie pour multicorps flexibles contraints étaient publiés en 2012 ; les poutres géométriquement exactes formulées sur SE(3) en 2014. Ces bases apportent cohérence géométrique et objectivité, mais leur adoption seule ne différencie plus Vinkulum. [Brüls et al., 2012](https://orbi.uliege.be/handle/2268/119869) ; [Sonneville et al., 2014](https://orbi.uliege.be/bitstream/2268/159471/1/pa_SonnevilleCardonaBruls2014_GeometricallyExactBeamFEOnSE3.pdf)

### L’apport précis de la modification σ de 2025

Holzinger, Arnold et Gerstmayr modifient l’étape géométrique de generalized-α. β et γ désignent les paramètres usuels du schéma. Un choix **σ = γ / (3β)** annule une composante du terme d’erreur dominant liée au crochet de Lie ; **σ = 1** simplifie le correcteur. Il ne s’agit pas d’annuler toute l’erreur ni d’augmenter automatiquement l’ordre. L’ordre deux suppose une initialisation adaptée. La publication indique une implémentation dans Exudyn : c’est une amélioration à rattraper et à évaluer. Certaines comparaisons du papier mêlent Python et C++, ce qui limite l’interprétation de leurs temps. [Holzinger et al., 2025](https://doi.org/10.1016/j.mechmachtheory.2025.106236)

**Expérience proposée :** comparer σ = 0, σ = 1 et σ = γ / (3β) dans le même Rust, avec les mêmes assemblages, tolérances, sorties et politique de Newton. Inclure un rotor incliné, une toupie, une boucle spatiale, une chaîne flexible et un témoin plan. Suivre trajectoire, phase, réactions, refus et coût à erreur commune. Le témoin plan permet de vérifier si le bénéfice se limite aux mouvements où les termes géométriques concernés sont significatifs.

### Conservation et précision ont des critères distincts

Les intégrateurs variationnels et les schémas énergie–moment répondent à des propriétés différentes. Une énergie bien conservée ne démontre ni une phase correcte ni la convergence du solveur non linéaire. Les comparaisons de Betsch et ses coauteurs illustrent la nécessité de contrôler plusieurs indicateurs. [Betsch et al., 2010](https://doi.org/10.1115/1.4001388)

Notre proposition est de choisir explicitement l’objectif temporel : vibrations à long terme, amortissement des hautes fréquences, événement rapide ou réponse d’intérêt. Toute correction d’énergie ou de moment doit respecter les symétries et forces réellement présentes. L’extension d’un cas conservatif aux actionneurs, au frottement ou à un couplage externe exige de nouvelles équations de bilan.

**Décision : première expérience temporelle.** Elle est assez proche du noyau pour donner rapidement une réponse utile, avec une contribution scientifique précise et des témoins défavorables inclus.

---

## 5. Flexibilité : précision, condensation et physique des sections

### Ce qu’apportent les formulations récentes

La formulation mixte de Humer, Steinbrecher et Pechstein introduit rotations internes discontinues, moments indépendants et courbure discrète. L’hybridation permet de conserver des poses nodales globales et d’éliminer des champs internes. Cette prépublication de mai 2026 fournit une piste pour gagner en ordre d’approximation sans multiplier autant les inconnues globales. [Humer et al., 2026](https://arxiv.org/html/2605.04573v1)

D’autres travaux enrichissent le modèle physique de section : directeurs extensibles, déformations de gauchissement et formulations mixtes. Choi et ses coauteurs étudient notamment une discrétisation objective dont certaines propriétés exactes d’énergie–moment dépendent du matériau choisi. Les B-splines généralisées sur SE(3) constituent une autre piste d’ordre élevé ; nous n’avons consulté que le résumé et des extraits de ce dernier article. [Choi et al., 2024](https://arxiv.org/html/2407.14637v2) ; [Ren et al., 2025](https://www.sciencedirect.com/science/article/pii/S0045782525002518)

### Ce que notre prototype démontre, et où il échoue

Le prototype autonome Rust d’ordres un et deux a été exécuté avec huit tests et 36 calculs Princeton. À l’ordre deux, trois éléments donnent **2,906 µm au bout**, mais **102,52 µm sur la grille de la ligne moyenne**. Huit éléments donnent **5,453 µm sur cette grille**, en **3,823 s**. La référence continue indépendante est évaluée sur 961 positions ; ce maximum échantillonné n’est pas une borne entre les points. [Prototype et archive Vinkulum](https://github.com/Brietat71/vinkulum/blob/723112bd51fc71b81d2757bdadc854529b055831/docs/POUTRE_MIXTE_PROTOTYPE.md)

Les sondages détectent aussi une raideur condensée négative artificielle à un rapport de raideurs de **10¹⁶**. Ce défaut concerne notre condensation en double précision, sans invalider la formulation en arithmétique exacte. La dynamique et les adjoints ne sont pas implémentés dans cet exemple. Aucun gain de vitesse externe n’est établi. [Limites mesurées du prototype](https://github.com/Brietat71/vinkulum/blob/723112bd51fc71b81d2757bdadc854529b055831/docs/POUTRE_MIXTE_PROTOTYPE.md)

### Travail scientifique à mener

Nous proposons de comparer la condensation actuelle à une résolution mixte complète correctement mise à l’échelle, puis à des éliminations numériquement plus stables. Les contrôles porteront sur la ligne entière, les courbures, les efforts, les modes rigides et la tangente. La précision terminale peut bénéficier d’une convergence plus rapide que celle du champ intérieur.

Il faudra ensuite dériver inerties et dynamique internes, charges distribuées, flambement, jonctions et limite du câble tendu. La condensation statique seule ne définit pas une réduction dynamique valide. Le choix entre ordre élevé, interpolation géométrique et enrichissement de section devra se faire sur **l’erreur physique obtenue par unité de coût**.

**Décision : maintenir un chantier de formulation, avec un verrou de stabilité avant intégration.** Le résultat actuel fournit un contre-exemple utile à l’idée qu’une publication récente suffit à produire un moteur supérieur.

---

## 6. Réduire les modèles flexibles sans perdre leurs propriétés

### Réduction non linéaire et hyperréduction

Réduire un modèle non linéaire à quelques coordonnées ne supprime pas nécessairement le coût de l’évaluation de ses forces. Deux opérations doivent être distinguées : réduire le nombre d’inconnues, puis réduire le travail d’évaluation du modèle complet.

Carlberg, Tuminaro et Boggs construisent la réduction à partir de la métrique cinétique, du potentiel, de la dissipation et des efforts pour conserver la structure lagrangienne. L’intérêt porte notamment sur des calculs répétés, qui amortissent la préparation du modèle. [Carlberg et al., 2015](https://arxiv.org/pdf/1401.8044)

ECSW échantillonne et pondère des contributions élémentaires afin de diminuer le coût non linéaire tout en conservant certaines propriétés de structure, de stabilité et de précision. Le résumé publié contient des garanties à paramètres fixés ; il ne fournit pas un certificat universel pour des chargements hors du domaine de construction. L’article intégral n’a pas été accessible pendant cette étude. [Farhat et al., 2015](https://doi.org/10.1002/nme.4820)

**Proposition Vinkulum :** conserver les coordonnées nécessaires aux articulations et interfaces, réduire les déformations internes, puis contrôler le modèle réduit avec un résidu indépendant. En cas d’insuffisance, enrichir la base ou revenir au modèle complet. Cette chaîne avec contrôle et repli est une proposition d’architecture ; les publications citées ne démontrent pas son fonctionnement général pour nos contacts et couplages.

### Variétés spectrales : un apport de mathématiques fondamentales

Les sous-variétés spectrales, ou SSM, donnent un cadre d’existence et d’unicité pour certaines réductions non linéaires locales autour d’un mouvement de référence. Les hypothèses de régularité et de non-résonance comptent ; une série formelle n’est pas toujours convergente. [Haller et Ponsioen, 2016](https://www.georgehaller.com/reprints/NNM.pdf)

Les méthodes de calcul ultérieures atteignent de grands modèles d’éléments finis sans calculer tout le spectre. Elles ouvrent l’analyse des résonances et des réponses forcées à coût réduit, au prix d’une préparation parfois importante. Cela ne constitue pas une réduction générale des impacts ou des changements de contact. [Jain et Haller, 2022](https://arxiv.org/pdf/2103.10264)

**Usage proposé :** identifier les résonances, interactions modales et zones dangereuses d’une transmission ou d’une structure tournante ; confronter les branches de réponse à une continuation du modèle complet. Il faut mesurer le coût hors ligne, le domaine d’amplitude couvert et le nombre d’analyses nécessaires pour amortir le modèle.

**Décision : potentiel élevé pour un avantage durable.** La réduction flexible devra servir la simulation, l’analyse et les sensibilités, avec un domaine de validité visible et des contrôles sur des charges non utilisées pour construire le modèle.

---

## 7. Contact : résoudre une loi physique clairement définie

### Les formulations ne sont pas interchangeables

Non-pénétration, cône de frottement circulaire et dissipation maximale imposent des conditions distinctes. La comparaison de Le Lidec et ses coauteurs montre pourquoi des relaxations apparemment proches peuvent produire des comportements différents. Une réaction hyperstatique sélectionnée numériquement n’est pas automatiquement la réaction physiquement pertinente. [Le Lidec et al., 2024](https://simple-robotics.github.io/publications/contact-models/static/paper/lelidec2024contacts.pdf)

Les formulations convexes d’Anitescu apportent des résultats de convergence sous hypothèses, sans être un Coulomb exact à pas fini. La formulation d’Acary et ses coauteurs associe un sous-problème conique à un point fixe extérieur : la convexité du sous-problème ne règle pas toute la convergence. [Anitescu, 2006](https://doi.org/10.1007/s10107-005-0590-7) ; [Acary et al., 2011](https://doi.org/10.1002/zamm.201000073)

### Trois familles modernes à confronter

**Proximal et complémentarité.** La méthode de Carpentier, Le Lidec et Montaut combine ADMM, régularisation proximale et correction de De Saxcé. Ses garanties convexes s’appliquent lorsque cette correction est retirée ; elles ne deviennent pas un théorème général pour le contact Coulomb exact. [Carpentier et al., 2024](https://www.roboticsproceedings.org/rss20/p108.pdf)

**Analyse variationnelle non lisse.** Le Newton semismooth* SCD utilise la structure locale d’une équation multivoque. Les résultats étudiés donnent une convergence locale superlinéaire sous régularité appropriée ; leurs expériences portent sur du contact élastique statique et la globalisation reste heuristique. C’est une piste pour les régimes difficiles, à dériver pour notre dynamique. [Gfrerer et al., 2023](https://link.springer.com/article/10.1007/s10589-022-00429-0)

**Lois régularisées cohérentes.** Les champs de contact irrotationnels analysent plusieurs approximations. La variante « Lagged » évite un décollement artificiel présent dans certaines formulations convexes, tout en restant découplée et régularisée. Elle est déjà utilisée dans Drake. Il faut comparer sa loi effective, y compris raideur et dissipation, avant de la choisir. [Castro et al., version 2025](https://arxiv.org/html/2312.03908v3)

### Une piste très récente qui exige encore des preuves

La décomposition de Song et ses coauteurs, publiée en préprint en juillet 2026, combine un sous-problème fortement convexe avec une correction extérieure non monotone. Les auteurs ne disposent pas du théorème monotone usuel pour cette correction. Sur leur arche, 47 % des sous-pas atteignent 10⁻⁶ avec le budget considéré ; la médiane est proche, à 1,1 × 10⁻⁶. Une incohérence entre texte et tableau d’annexe reste non résolue. Nous retenons l’idée de décomposition pour expérimentation, sans reprendre une promesse de robustesse générale. [Song et al., 2026](https://arxiv.org/html/2607.19599v1)

**Décision : résoudre ensemble normale et tangente, puis juger avec un résidu physique indépendant.** Distinguer la compliance qui modélise un matériau et la régularisation introduite pour calculer. Les paradoxes de Painlevé rappellent que le modèle parfaitement rigide peut lui-même perdre existence ou unicité ; changer de solveur ne supprime pas ce problème physique. [Hogan et Uldall Kristiansen, 2017](https://arxiv.org/pdf/1610.00143)

---

## 8. Collisions et événements : des garanties utiles à l’utilisateur

### Contrôler la trajectoire réellement parcourue

La détection continue par inclusion propose des tests conservatifs sur les primitives et trajectoires étudiées, avec une évaluation explicite des faux négatifs et des faux positifs. Son intérêt est de transformer « la collision semble détectée » en un contrat géométrique testable. [Wang et al., 2021](https://continuous-collision-detection.github.io/tight_inclusion/CCD-benchmark-paper-350ppi.pdf)

Rigid IPC traite les mouvements rigides courbes en bornant l’écart entre trajectoire et corde, puis en adaptant le test ou la subdivision. Cette idée doit être redérivée pour l’interpolation temporelle exacte de Vinkulum. La preuve ne s’étend pas automatiquement à toutes nos capsules ; les auteurs eux-mêmes mentionnent les capsules parmi les extensions. [Ferguson et al., 2021](https://ipc-sim.github.io/rigid-ipc/assets/rigid_ipc_paper_350ppi.pdf)

**Contre-épreuve proposée :** une capsule tourne autour d’un support. Ses poses de début et de fin sont libres, mais son arc traverse un obstacle. Il faut tester les positions intermédiaires pertinentes avec une borne conservative. Ajouter rasance, primitives presque parallèles, très petits jeux et grandes coordonnées. Publier les collisions manquées, les fausses alertes et le coût de subdivision.

### Une trajectoire sans intersection peut encore être mécaniquement fausse

Les travaux sur la convergence d’IPC étudient aussi le raffinement de surface, de barrière et de temps. Un autre article montre des forces tangentielles parasites capables de freiner un objet dans un cas sans frottement. La garantie géométrique doit donc être accompagnée de tests sur les forces et le travail. [Li et al., 2023](https://arxiv.org/html/2307.15908v1) ; [Du et al., 2023](https://arxiv.org/pdf/2308.01696)

### Dériver aussi l’instant de l’événement

Lorsqu’un paramètre change, l’impact peut survenir plus tôt ou plus tard. Différencier seulement la loi de remise à jour à instant fixe omet cet effet. La matrice de saltation ajoute cette correction de temps à la linéarisation d’un système hybride. Les hypothèses classiques exigent des lois suffisamment régulières, une traversée transverse et une séquence locale d’événements inchangée. Rasance, simultanéité et accumulation d’impacts sortent de cette garantie. [Kong et al., 2024](https://arxiv.org/html/2306.06862v3)

**Proposition Vinkulum :** fournir un gradient accompagné de son statut : régime régulier vérifié, dérivée directionnelle lorsque disponible, ou sensibilité non fiable dans le régime rencontré. Les différences finies de contrôle devront recalculer indépendamment les instants d’impact. Une dérivée formellement obtenue par AD ne devra pas être présentée comme une sensibilité physique validée.

**Décision : traiter géométrie, mécanique et sensibilité comme trois contrats à contrôler ensemble.** Cette capacité intéresse la conception, l’identification et le contrôle de tous les mécanismes avec butées, chocs ou frottement.

---

## 9. Différentiation : obtenir le gradient au coût du problème convergé

### Dériver la solution, avec les hypothèses nécessaires

Pour un état convergé z et des paramètres θ, écrivons **F(z, θ) = 0**. La sensibilité dans une direction résout **Fz δz = −Fθ δθ**. Pour un objectif scalaire L, un adjoint résout **Fzᵀ λ = Lzᵀ**, puis donne **dL/dθ = Lθ − λᵀ Fθ**. Fz, Fθ, Lz et Lθ désignent les dérivées partielles correspondantes. La formulation modulaire de Blondel et ses coauteurs organise ces opérations autour des conditions du problème. Elle suppose la régularité et l’inversibilité locales appropriées ; près d’une singularité, une petite erreur d’état peut produire une grande erreur de gradient. [Blondel et al., 2022](https://papers.neurips.cc/paper_files/paper/2022/file/228b9279ecf9bbafe582406850c57115-Paper-Conference.pdf)

Le théorème des fonctions implicites n’est pas nouveau. L’apport pratique récent réside dans l’organisation des opérateurs, de la différentiation et des résolutions. **Notre proposition** est de réutiliser les blocs et préconditionneurs du calcul primal lorsque cela est mathématiquement valide, avec un contrôle séparé du résidu adjoint.

### Ne calculer que les dérivées nécessaires

Un produit Jacobienne-vecteur, ou JVP, décrit l’effet d’une perturbation sans construire toute la Jacobienne. Son opération transposée, VJP, sert à remonter un objectif. Les Hessiennes complètes restent utiles pour certains blocs locaux ; les produire partout peut imposer un coût inutile.

Enzyme montre comment différencier une représentation LLVM déjà optimisée. Le support Rust documenté demeure expérimental, sur certaines plateformes et chaînes nightly. Une expérience isolée peut comparer ce calcul aux jets actuels ; la production ne doit pas dépendre d’une disponibilité stable qui n’est pas établie. [Moses et Churavy, 2020](https://papers.neurips.cc/paper_files/paper/2020/hash/9332c513ef44b682e9347822c2e457ac-Abstract.html) ; [Guide officiel Rust](https://rustc-dev-guide.rust-lang.org/autodiff/installation.html)

### Expérience qui mesure une valeur produit

Comparer un objectif dépendant de 1, 10 puis 1 000 paramètres, sur un modèle régulier, une boucle redondante et un cas avec impact isolé. Vérifier les développements de Taylor et l’identité entre produit direct et produit transposé ; mesurer mémoire, reconstruction éventuelle de trajectoire et temps complet.

Sur les contraintes redondantes, définir la convention et les grandeurs identifiables avant de différencier. Pour les impacts, appliquer les conditions de la page précédente. Mesurer enfin le nombre de simulations nécessaires pour atteindre un objectif de conception, avec des initialisations identiques.

**Décision : priorité forte, après stabilisation des opérateurs primaux.** Le noyau possède déjà des adjoints ; l’avance à rechercher concerne leur portée, leur coût et la validité de leurs résultats dans les régimes où l’utilisateur en a besoin.

---

## 10. Informatique de pointe : stabilité, matériel et apprentissage

### Accélérer une résolution en gardant un contrôle numérique

Le raffinement itératif à plusieurs précisions sépare précision de factorisation, du calcul du résidu et des itérations de Krylov. Il peut permettre une arithmétique moins coûteuse dans une partie du calcul, sous ses hypothèses de convergence. Il ne garantit pas la résolution d’un système KKT singulier quelconque. [Amestoy et al., 2024](https://doi.org/10.1137/23M1549079)

Les travaux d’Epperly, Greenbaum et Nakatsukasa sur les équations normales préconditionnées établissent des résultats de stabilité avec rang numérique plein et préconditionnement adapté. Les solveurs aléatoires SPIR/FOSSILS apportent aussi des garanties pour certains moindres carrés surdéterminés. Ces résultats récents ne lèvent pas automatiquement les difficultés de projection redondante de Vinkulum. [Epperly et al., 2026, systèmes linéaires](https://arxiv.org/html/2502.17767v2) ; [Epperly et al., 2026, moindres carrés](https://arxiv.org/html/2406.03468v3)

**Transfert proposé :** commencer par les blocs réguliers bien identifiés, recalculer le résidu original en précision suffisante et conserver un repli robuste. Ne modifier le solveur des cas singuliers qu’après une analyse propre de leur rang et de leur métrique.

### Globaliser Newton et compter son travail réel

Les stratégies récentes Newton–CG pour maillages courbes associent préconditionnement, adaptation de précision et recherche de pas. Leur domaine d’origine est une optimisation de maillage, pas un système multicorps contraint général. Nous proposons d’en reprendre les principes de contrôle du travail, en dérivant la globalisation pour nos équations. [Aparicio-Estrems et al., 2024](https://arxiv.org/html/2403.13654v1)

### Deux bancs matériels, puis un arbitrage

Le banc CPU évaluera la latence d’un mécanisme connecté et sa mémoire. Le banc GPU fera varier séparément la taille de chaque mécanisme et le nombre de simulations simultanées. Préparation, compilation, transferts, précision, rejets et coût des gradients seront comptés. Une accélération de mille environnements indépendants ne mesure pas la latence d’une grande assemblée.

### Apprentissage : une aide contrôlée à la résolution

Les préconditionneurs fondés sur DeepONet apprennent un opérateur approché ou un espace grossier. Les propriétés de symétrie et de positivité dépendent de la construction ; une application neuronale non linéaire ne se traite pas comme un préconditionneur linéaire fixe. [Kopaničáková et Karniadakis, 2025](https://arxiv.org/html/2401.02016v2)

Nous proposons d’apprendre initialisations, bases réduites ou préconditionneurs, tout en laissant les équations physiques décider de l’acceptation. Les échecs connus des PINN montrent qu’une pénalisation de la physique dans une fonction de coût ne garantit pas son respect. [Krishnapriyan et al., 2021](https://www.stat.berkeley.edu/~mmahoney/pubs/failure-modes-neurips21.pdf)

**Décision : accélérations conditionnelles.** Compter apprentissage, hors-domaine et repli ; comparer au meilleur préconditionneur classique disponible dans le même protocole.

---

## 11. Physique et multiphysique : préserver les échanges utiles

### Des ports énergétiques pour composer le système

La formulation port-hamiltonienne de Latussek, Kinon et Betsch organise les interconnexions multicorps autour des échanges d’énergie et d’une représentation par directeurs. Les propriétés exactes du point milieu présentées dépendent notamment du caractère au plus quadratique du Hamiltonien et des contraintes. La commutation entre interconnexion et discrétisation exige le même schéma aux sous-systèmes. Ce préprint de mars 2026 ne certifie donc pas toutes les lois constitutives ni tous les couplages multi-rythmes. [Latussek et al., 2026](https://arxiv.org/html/2603.12841v1)

Le bénéfice recherché pour Vinkulum est de disposer d’une description cohérente des puissances aux interfaces mécanique, électrique, hydraulique et commande. Les bilans doivent distinguer travail imposé, stockage, dissipation physique et erreur numérique.

### Contrôler la cosimulation avec les informations disponibles

ECCO utilise le défaut d’énergie échangée pour estimer l’erreur et adapter le macro-pas de cosimulation sans accès aux Jacobiennes internes. Cette méthode fournit un indicateur opérationnel, sans établir une stabilité universelle en présence de couplages raides ou de boucles algébriques. [Sadjina et al., 2017](https://arxiv.org/pdf/1602.06434)

Les estimateurs a posteriori fondés sur les adjoints GARK répartissent l’erreur spatiale et temporelle selon une quantité d’intérêt. Le cadre consulté et ses exemples ne démontrent pas directement le même résultat pour les DAE multicorps avec impacts. [Narayanamurthi et al., 2020](https://arxiv.org/html/2001.08824v1)

**Proposition :** commencer par des interfaces de puissance mesurables et confronter une résolution monolithique à des partitions. Suivre simultanément énergie, phase, contraintes et observable métier. Adapter ensuite pas, sous-pas ou résolution locale selon l’erreur qui affecte réellement cette observable. Un bilan énergétique presque fermé peut accompagner une phase incorrecte.

### Incertitudes et multi-fidélité

La gestion optimale de modèles multi-fidélité répartit les évaluations en fonction de leurs coûts et corrélations. Elle peut réduire le coût d’une statistique d’un modèle de haute fidélité en combinant plusieurs niveaux. L’absence de biais concerne cette sortie de haute fidélité, pas automatiquement la réalité physique. Pour ce travail, seuls le résumé et la bibliographie de l’auteur ont été accessibles. [Peherstorfer et al., 2016](https://cims.nyu.edu/~pehersto/)

Cette piste devient intéressante pour la dispersion de fabrication, la fiabilité ou l’identification répétée. Elle réduit le nombre de calculs coûteux nécessaires à une étude ; elle n’accélère pas à elle seule un pas du noyau.

**Décision : construire des bilans d’interface et des modèles réduits contrôlés avant une stratégie multi-rythme ou multi-fidélité générale.** Le progrès physique recherché tient à la fidélité des modèles et de leurs interactions, avec une erreur numérique maîtrisée.

---

## 12. Programme d’expériences qui peut confirmer ou invalider les choix

Les expériences ci-dessous sont **proposées**. Leur ordre suit les dépendances techniques ; aucun gain chiffré de ces transferts n’est annoncé. Elles complètent le prototype mixte déjà archivé.

| Expérience | Contrôle décisif | Motif d’abandon ou de révision |
|---|---|---|
| Variantes temporelles géométriques | Même code, erreurs de phase, position et réaction ; témoins plans et spatiaux | Bénéfice limité à un seul rotor ou annulé par Newton |
| Poutre mixte, condensation et ordre | Modèle mixte complet de référence ; ligne, efforts, modes rigides, anisotropie et raideurs extrêmes | Instabilité artificielle, coût supérieur sans gain utile |
| Graphe, redondance et résolution | Arbres, boucles et maillages ; duplication, permutation, unités ; petite référence SVD | Résultat physique modifié par préconditionnement ou sélection implicite des réactions |
| Réduction et hyperréduction | Charges inédites, inversion de charge, préflambement ; résidu indépendant et repli | Erreur non détectée ou préparation non amortie |
| Contact et trajectoires | Adhérence/glissement, repère tangent tournant, rapports de masses, capsule en rotation | Loi effective différente, travail parasite ou collision manquée |
| Sensibilités | Tests de Taylor, identité JVP/VJP, impacts recalculés ; singularités et rasance | Gradient déclaré valide sans convergence du contrôle |
| Couplage et multi-rythme | Référence monolithique ; interfaces raides, retard et réponse directe | Bon bilan d’énergie mais mauvaise phase ou instabilité |
| CPU, GPU et apprentissage | Taille connectée et lots indépendants ; préparation, mémoire, précision et repli inclus | Gain disparaissant dans le coût complet ou hors du domaine appris |

---

## 13. Protocole de comparaison et critères d’intégration

Un banc doit permettre de choisir une méthode et de reproduire ce choix. Les critères suivants constituent notre proposition de protocole commun ; ils ne sont pas des résultats scientifiques déjà obtenus.

1. **Fixer la question physique.** Géométrie, matériau, frottement, chargement, unités, horizon et observables doivent être équivalents. Un contact régularisé et un contact rigide exact constituent des modèles différents.
2. **Fixer la référence et sa marge.** Solution analytique si disponible, modèle continu indépendant ou raffinements concordants. Distinguer erreur estimée sur une grille et borne certifiée.
3. **Explorer les réglages admissibles.** Publier les compromis temps–erreur et mémoire–erreur, les tentatives échouées et les configurations retenues. Séparer initialisation, préparation, calcul, sorties et coût total.
4. **Exiger plusieurs familles.** Un gain doit survivre à des témoins qui n’ont pas guidé son développement. Les seuils d’intégration sont décidés avant de regarder les résultats finaux.
5. **Comparer l’usage complet.** Pour une étude répétée ou une optimisation, compter préparation, gradients, nombre d’appels et coût de validation. Conserver les versions, matériels, modèles et données brutes.

Les comparateurs devront être choisis par capacité et par accès effectif : MBDyn et Exudyn pour la mécanique flexible ; Chrono pour certains couplages ; Drake ou MuJoCo pour les contacts et sensibilités pris en charge ; produits industriels pour leurs cas pertinents dès qu’une exécution reproductible est possible.

### Accepter une amélioration sans masquer son domaine

Une méthode sera candidate à l’intégration lorsqu’elle apporte un gain utile à erreur commune, avec les contrôles physiques requis et sans régression inexpliquée sur les témoins généralistes. Les gains et pertes seront publiés ensemble. Une stratégie spécialisée pourra rester disponible sur son domaine démontré, avec une sélection explicite et un repli contrôlé.

Pour les cas raides ou singuliers, un refus diagnostiqué vaut une information distincte d’une solution acceptée. Les statistiques doivent conserver les deux. Pour les modèles réduits ou appris, le coût amorti sera donné en fonction du nombre de réutilisations, afin que l’utilisateur puisse juger sa propre étude.

Une revendication de supériorité générale exigerait également des preuves sur la couverture des modèles, la qualité des importations, le diagnostic, le temps réel et les chaînes de travail industrielles. Le présent rapport ne mesure pas ces dimensions.

---

## 14. Décision d’investissement et limites de l’étude

### Ordre de travail recommandé

**Première étape : un noyau mesurable et une formulation stable.** Instrumenter les coûts de Newton et des dérivées, expérimenter la modification temporelle récente, traiter le coût et la stabilité de la poutre mixte. Maintenir les témoins de la livraison 0.7.2.

**Deuxième étape : structure du système et sensibilités.** Unifier les blocs physiques et leurs produits directs/transposés, exploiter le graphe, rendre explicites rang et métrique, puis dériver les solutions convergées avec des contrôles dédiés.

**Troisième étape : flexibilité réduite et contact exigeant.** Construire une réduction contrôlée sur charges nouvelles ; résoudre normale et tangente ensemble ; traiter les trajectoires de rotation et les événements. Cette étape doit produire une amélioration du coût de conception sur plusieurs familles généralistes.

**Étapes conditionnelles : multiphysique avancée, GPU, précision mixte et apprentissage.** Les engager lorsque les mesures identifient un coût qu’elles peuvent réellement réduire. Les variétés spectrales et la multi-fidélité ajoutent des capacités d’analyse au-dessus du noyau.

### Ce qui constituerait le véritable avantage

Notre proposition est une architecture où le modèle, la résolution et les sensibilités partagent les mêmes structures mécaniques, avec une erreur physique contrôlée. La différenciation durable viendrait de cette cohérence, de la diversité des cas validés et du coût complet obtenu. Il reste à démontrer qu’elle procure un avantage face aux implémentations concurrentes, qui continuent elles aussi de progresser.

Les publications justifient les mécanismes locaux et certaines garanties. **Elles ne prouvent pas la performance de leur combinaison dans Vinkulum.** Le moteur généraliste devra conserver plusieurs stratégies lorsque leurs domaines de validité diffèrent. Une sélection automatique n’est souhaitable qu’avec des critères mesurables et un repli vérifié.

### Limites et niveau de preuve

L’étude s’appuie sur **58 notices de provenance** : articles, prépublications, documentation officielle et mesures propres. Les trois volets — géométrie et physique, contact et contraintes, calcul et différentiation — ont été explorés séparément puis recoupés. Les affirmations les plus déterminantes ont été relues dans les sources primaires ; les liens renvoient aux sources ou à leurs pages de publication, et les notes précisent la version consultée.

L’accès intégral a manqué pour certains articles, notamment ECSW, les B-splines de 2025 et la multi-fidélité. Leurs conclusions sont limitées au niveau consulté. Les prépublications de 2026 sont identifiées comme telles. Un article annoncé dans un numéro de novembre 2026 a été exclu du corpus consolidé, faute de date d’accès public vérifiée avant le 7 septembre.

Les comparaisons industrielles reposent sur les capacités documentées. Le code des solveurs concurrents n’a pas été consulté. Les résultats locaux existants ne couvrent pas encore contact, multiphysique, gradients et temps réel face à ces concurrents.

La recherche a été arrêtée lorsque chaque décision importante disposait d’un support primaire ou d’une limite explicite, et que les contradictions restantes étaient circonscrites. **La prochaine preuve utile doit venir des expériences proposées.** Une collecte supplémentaire de publications ne permettrait pas d’annoncer les gains encore non mesurés.

---

## Sources et accès

Les références suivent leur première apparition. Les liens proches des affirmations et ceux ci-dessous sont cliquables. Les pages de documentation évoluent ; elles sont datées par leur consultation du 7 septembre 2026. Les détails d’accès indiquent le niveau effectivement examiné, sans prétendre à une lecture intégrale uniforme du corpus.

**Projet Vinkulum · 2026-09-07.** [Vinkulum — état et limites du noyau](https://github.com/Brietat71/vinkulum/blob/723112bd51fc71b81d2757bdadc854529b055831/README.md). README et sources au commit 723112b. Accès : Lecture README et points d'entrée du code propre ; aucun code concurrent.

**Dassault Systèmes · s.d. ; consulté le 7 septembre 2026.** [Simpack Base Modules and Wizard](https://www.3ds.com/products/simulia/simpack/core). Documentation produit officielle. Accès : Page officielle.

**Dassault Systèmes · s.d. ; consulté le 7 septembre 2026.** [Simpack Flexible Body Modules](https://www.3ds.com/products/simulia/simpack/flexible-body). Documentation produit officielle. Accès : Page officielle.

**Hexagon · 2025.1 ; page consultée le 7 septembre 2026.** [Adams — Products and 2025.1 update](https://nexus.hexagon.com/home/product/adams/). Documentation produit officielle. Accès : Page officielle.

**FunctionBay · s.d. ; consulté le 7 septembre 2026.** [RecurDyn — All Product Components](https://support.functionbay.com/en/page/single/232/recurdyn-all-products). Documentation produit officielle. Accès : Page officielle.

**Siemens · s.d. ; consulté le 7 septembre 2026.** [Motion simulation](https://www.siemens.com/en-us/products/simcenter/simulation-test/motion-simulation/). Documentation produit officielle. Accès : Page officielle.

**Équipe MBDyn, Politecnico di Milano · s.d. ; consulté le 7 septembre 2026.** [MBDyn — Introduction and capabilities](https://www.mbdyn.org/). Documentation officielle. Accès : Page officielle.

**Gerstmayr et équipe Exudyn · publication 2024 ; documentation consultée le 7 septembre 2026.** [Exudyn — flexible multibody dynamics with Python and C++](https://exudyn.readthedocs.io/en/latest/docs/RST/Exudyn.html). Multibody System Dynamics et documentation officielle. Accès : Documentation officielle ; versions de documentation et d'essais à figer séparément.

**Project Chrono · 2026-03-27.** [Chrono 10.0.0 released](https://projectchrono.org/news/). Actualités officielles. Accès : Annonce et documentation des modules consultées.

**Équipe Drake · s.d. ; consulté le 7 septembre 2026.** [Hydroelastic Contact User Guide](https://drake.mit.edu/doxygen_cxx/group__hydroelastic__user__guide.html). Documentation officielle. Accès : Guide officiel.

**Équipe MuJoCo · s.d. ; consulté le 7 septembre 2026.** [MuJoCo XLA — Implementations and Feature Parity](https://mujoco.readthedocs.io/en/stable/mjx.html). Documentation officielle. Accès : Tableau des fonctionnalités consulté.

**Équipes MuJoCo, Google DeepMind et NVIDIA · s.d. ; consulté le 7 septembre 2026.** [MuJoCo Warp — When To Use MJWarp?](https://mujoco.readthedocs.io/en/latest/mjwarp/). Documentation officielle. Accès : Sections débit/latence/mécanismes/différentiabilité relues par le coordinateur.

**Newton Developers · 2026 ; consulté le 7 septembre 2026.** [Newton Physics — Solvers](https://newton-physics.github.io/newton/stable/solvers/index.html). Documentation officielle. Accès : Guide et tableau des fonctionnalités consultés.

**Projet Vinkulum · 2026-09-07.** [Vinkulum 0.7.2 face à MBDyn — confrontation exécutée](https://github.com/Brietat71/vinkulum/blob/723112bd51fc71b81d2757bdadc854529b055831/docs/CONFRONTATION_MBDYN_0.7.2.md). Dépôt Vinkulum, mesures locales archivées. Accès : Sources locales et archives vérifiées ; commit poussé.

**Lall, Krysl, Marsden · 2003-10-01.** [Structure-preserving model reduction for mechanical systems](https://doi.org/10.1016/S0167-2789(03)00227-6). Physica D 184, 304–318. Accès : Résumé et introduction éditeur.

**Knoll, Keyes · 2004-01-20.** [Jacobian-free Newton–Krylov methods: a survey of approaches and applications](https://doi.org/10.1016/j.jcp.2003.08.010). Journal of Computational Physics. Accès : Résumé éditeur.

**Carpentier, Budhiraja, Mansard · 2021-07.** [Proximal and Sparse Resolution of Constrained Dynamic Equations](https://www.roboticsproceedings.org/rss17/p017.pdf). Robotics: Science and Systems. Accès : Texte intégral ; lecture coordonnateur de la formulation et de l'introduction.

**Sathya, Bruyninckx, Decré, Pipeleers · 2024.** [Efficient Constrained Dynamics Algorithms Based on an Equivalent LQR Formulation Using Gauss’ Principle of Least Constraint](https://doi.org/10.1109/TRO.2023.3335652). IEEE Transactions on Robotics 40. Accès : Texte intégral institutionnel ; complexité et limites §IX.E.

**Fong, Saunders · 2011-10-27.** [LSMR: An Iterative Algorithm for Sparse Least-Squares Problems](https://stanford.edu/group/SOL/software/lsmr/LSMR-SISC-2011.pdf). SIAM Journal on Scientific Computing. Accès : Texte intégral et documentation auteurs.

**Choi, Paige, Saunders · 2011-08-04.** [MINRES-QLP: A Krylov Subspace Method for Indefinite or Singular Symmetric Systems](https://doi.org/10.1137/100787921). SIAM Journal on Scientific Computing. Accès : Article et [analyse complémentaire du préconditionnement](https://www.mcs.anl.gov/papers/P3027-0812.pdf).

**Brüls, Cardona, Arnold · 2012.** [Lie group generalized-alpha time integration of constrained flexible multibody systems](https://orbi.uliege.be/handle/2268/119869). Mechanism and Machine Theory 48, 121–137. Accès : Résumé institutionnel et extraits consultés.

**Sonneville, Cardona, Brüls · 2014.** [Geometrically exact beam finite element formulated on the special Euclidean group SE(3)](https://orbi.uliege.be/bitstream/2268/159471/1/pa_SonnevilleCardonaBruls2014_GeometricallyExactBeamFEOnSE3.pdf). CMAME 268, 451–474. Accès : Postprint intégral institutionnel.

**Holzinger, Arnold, Gerstmayr · 2025-10.** [σ-modified Lie group generalized-α methods for constrained multibody systems](https://doi.org/10.1016/j.mechmachtheory.2025.106236). Mechanism and Machine Theory 217, 106236. Accès : Texte intégral déposé par le premier auteur sur ResearchGate, relu par le coordinateur ; éditeur inaccessible.

**Betsch, Hesch, Sänger, Uhlar · 2010.** [Variational Integrators and Energy-Momentum Schemes for Flexible Multibody Dynamics](https://doi.org/10.1115/1.4001388). Journal of Computational and Nonlinear Dynamics 5, 031001. Accès : Texte auteur consulté.

**Humer, Steinbrecher, Pechstein · 2026-05-06.** [Mixed Finite Elements for Geometrically Exact Beams using Discontinuous Rotations and Discrete Curvature](https://arxiv.org/html/2605.04573v1). arXiv 2605.04573v1, prépublication. Accès : Théorie et résultats lus ; dérivation propre et prototype exécuté dans Vinkulum.

**Choi, Klinkel, Klarmann, Sauer · 2024-11-19 (v2).** [An objective isogeometric mixed finite element formulation for nonlinear elastodynamic beams with incompatible warping strains](https://arxiv.org/html/2407.14637v2). arXiv 2407.14637v2. Accès : Texte intégral.

**Ren, Yuan, Liu · 2025.** [Geometrically exact beam finite element with generalized B-spline interpolation on the special Euclidean group SE(3)](https://www.sciencedirect.com/science/article/pii/S0045782525002518). CMAME 441, 117979. Accès : Aperçu éditeur ; dérivation complète non auditée.

**Projet Vinkulum · 2026-09-07.** [Poutre mixte : transfert d'une formulation de 2026](https://github.com/Brietat71/vinkulum/blob/723112bd51fc71b81d2757bdadc854529b055831/docs/POUTRE_MIXTE_PROTOTYPE.md). Prototype autonome, commit 723112b. Accès : Prototype exécuté, 36 runs Princeton, 28 sondages et 8 tests ; archive vérifiée.

**Carlberg, Tuminaro, Boggs · 2015 ; prépublication 2014.** [Preserving Lagrangian Structure in Nonlinear Model Reduction with Application to Structural Dynamics](https://arxiv.org/pdf/1401.8044). SIAM Journal on Scientific Computing. Accès : Texte intégral.

**Farhat, Chapman, Avery · 2015-03-03.** [Structure-preserving, stability, and accuracy properties of the energy-conserving sampling and weighting method for the hyper reduction of nonlinear finite element dynamic models](https://doi.org/10.1002/nme.4820). IJNME 102, 1077–1110. Accès : Résumé éditeur détaillé ; intégral inaccessible.

**Haller, Ponsioen · 2016-08-01.** [Nonlinear normal modes and spectral submanifolds: existence, uniqueness and use in model reduction](https://www.georgehaller.com/reprints/NNM.pdf). Nonlinear Dynamics 86, 1493–1534. Accès : Texte intégral auteur.

**Jain, Haller · 2022 ; en ligne 2021-10-12.** [How to compute invariant manifolds and their reduced dynamics in high-dimensional finite element models](https://arxiv.org/pdf/2103.10264). Nonlinear Dynamics. Accès : Texte intégral.

**Le Lidec, Jallet, Montaut, Laptev, Schmid, Carpentier · 2024.** [Contact Models in Robotics: a Comparative Analysis](https://simple-robotics.github.io/publications/contact-models/static/paper/lelidec2024contacts.pdf). IEEE Transactions on Robotics. Accès : Texte intégral.

**Anitescu · 2006 ; en ligne 2005.** [Optimization-based simulation of nonsmooth rigid multibody dynamics](https://doi.org/10.1007/s10107-005-0590-7). Mathematical Programming. Accès : Texte intégral ; théorème 3 et hypothèses A1–A5.

**Acary, Cadoux, Lemaréchal, Malick · 2011 ; en ligne 2010.** [A formulation of the linear discrete Coulomb friction problem via convex optimization](https://doi.org/10.1002/zamm.201000073). ZAMM. Accès : Texte intégral §§3–5.

**Carpentier, Le Lidec, Montaut · 2024-07.** [From Compliant to Rigid Contact Simulation: a Unified and Efficient Approach](https://www.roboticsproceedings.org/rss20/p108.pdf). Robotics: Science and Systems. Accès : Texte intégral ; §VI.

**Gfrerer, Mandlmayr, Outrata, Valdman · 2023 ; en ligne 2022-11-07.** [On the SCD semismooth* Newton method for generalized equations with application to a class of static contact problems with Coulomb friction](https://link.springer.com/article/10.1007/s10589-022-00429-0). Computational Optimization and Applications. Accès : Texte intégral éditeur.

**Castro, Han, Masterjohn · 2025-07-15 (v3) ; prépublication initiale 2023.** [Irrotational Contact Fields](https://arxiv.org/html/2312.03908v3). arXiv 2312.03908v3. Accès : Texte intégral ; tableau I recoupé par le coordinateur.

**Song, Fan, Ascher, Pai · 2026-07-21.** [A Splitting Architecture for Exact Reduced Coulomb Friction](https://arxiv.org/html/2607.19599v1). arXiv 2607.19599v1, prépublication. Accès : Texte intégral ; non-monotonie et tableaux de convergence vérifiés par le coordinateur.

**Hogan, Uldall Kristiansen · 2017.** [On the regularization of impact without collision: the Painlevé paradox and compliance](https://arxiv.org/pdf/1610.00143). Proceedings of the Royal Society A. Accès : Texte intégral.

**Wang, Ferguson, Schneider, Jiang, Attene, Panozzo · 2021-10.** [A Large Scale Benchmark and an Inclusion-Based Algorithm for Continuous Collision Detection](https://continuous-collision-detection.github.io/tight_inclusion/CCD-benchmark-paper-350ppi.pdf). ACM Transactions on Graphics, DOI 10.1145/3460775. Accès : Texte intégral.

**Ferguson et al. · 2021-08.** [Intersection-free Rigid Body Dynamics](https://ipc-sim.github.io/rigid-ipc/assets/rigid_ipc_paper_350ppi.pdf). ACM Transactions on Graphics, DOI 10.1145/3450626.3459802. Accès : Texte intégral ; §4.3 et limites §7.

**Li, Ferguson, Schneider, Langlois, Zorin, Panozzo, Jiang, Kaufman · 2023-07-29.** [Convergent Incremental Potential Contact](https://arxiv.org/html/2307.15908v1). arXiv 2307.15908v1. Accès : Texte intégral.

**Du, Li, Coros, Thomaszewski · 2023-08.** [No Free Slide: Spurious Contact Forces in Incremental Potential Contact](https://arxiv.org/pdf/2308.01696). arXiv 2308.01696. Accès : Texte intégral ; contre-exemple §6.

**Kong, Payne, Zhu, Johnson · 2024.** [Saltation Matrices: The Essential Tool for Linearizing Hybrid Dynamical Systems](https://arxiv.org/html/2306.06862v3). Proceedings of the IEEE, DOI 10.1109/JPROC.2024.3440211. Accès : Texte intégral ; hypothèses et §III.

**Blondel et al. · 2022.** [Efficient and Modular Implicit Differentiation](https://papers.neurips.cc/paper_files/paper/2022/file/228b9279ecf9bbafe582406850c57115-Paper-Conference.pdf). NeurIPS. Accès : Texte intégral ; équations et erreur de dérivation consultées.

**Moses, Churavy · 2020.** [Instead of Rewriting Foreign Code for Machine Learning, Automatically Synthesize Fast Gradients](https://papers.neurips.cc/paper_files/paper/2020/hash/9332c513ef44b682e9347822c2e457ac-Abstract.html). NeurIPS. Accès : Publication et article.

**Rust Compiler Development Guide · s.d. ; consulté le 7 septembre 2026.** [Automatic differentiation: Installation](https://rustc-dev-guide.rust-lang.org/autodiff/installation.html). Documentation officielle Rust. Accès : Page officielle consultée et recoupée par le coordinateur.

**Amestoy, Buttari, Higham, L’Excellent, Mary, Vieublé · 2024-02-08.** [Five-Precision GMRES-Based Iterative Refinement](https://doi.org/10.1137/23M1549079). SIAM Journal on Matrix Analysis and Applications. Accès : Résumé publié et prépublication consultés ; hypothèses détaillées à relire avant portage.

**Epperly, Greenbaum, Nakatsukasa · 2026-05-26 ; prépublication 2025.** [Stable algorithms for general linear systems by preconditioning the normal equations](https://arxiv.org/html/2502.17767v2). Numerische Mathematik, DOI 10.1007/s00211-026-01547-1. Accès : Prépublication intégrale ; date de publication recoupée institutionnellement ; final éditeur sous abonnement.

**Epperly, Meier, Nakatsukasa · 2026-02 ; en ligne 2025-09-30.** [Fast randomized least-squares solvers can be just as accurate and stable as classical direct solvers](https://arxiv.org/html/2406.03468v3). Communications on Pure and Applied Mathematics 79(2). Accès : Texte intégral.

**Aparicio-Estrems, Gargallo-Peiró, Roca · 2024.** [A Globalized and Preconditioned Newton-CG Solver for Metric-Aware Curved High-Order Mesh Optimization](https://arxiv.org/html/2403.13654v1). Computer-Aided Design 168. Accès : Texte intégral.

**Kopaničáková, Karniadakis · 2025-02-20.** [DeepONet Based Preconditioning Strategies for Solving Parametric Linear Systems of Equations](https://arxiv.org/html/2401.02016v2). SIAM Journal on Scientific Computing. Accès : Texte intégral ; §4 et conclusion.

**Krishnapriyan, Gholami, Zhe, Kirby, Mahoney · 2021.** [Characterizing possible failure modes in physics-informed neural networks](https://www.stat.berkeley.edu/~mmahoney/pubs/failure-modes-neurips21.pdf). NeurIPS. Accès : Texte intégral.

**Latussek, Kinon, Betsch · 2026-03-13.** [Port-Hamiltonian multibody dynamics: Lagrangian formulation, consistent interconnection, structure-preserving simulation and index-reduction](https://arxiv.org/html/2603.12841v1). arXiv 2603.12841v1, prépublication. Accès : Texte intégral ; hypothèses recoupées par le coordinateur.

**Sadjina, Kyllingstad, Pedersen, Skjong · 2017 ; prépublication 2016.** [Energy Conservation and Power Bonds in Co-Simulations: Non-Iterative Adaptive Step Size Control and Error Estimation](https://arxiv.org/pdf/1602.06434). Engineering with Computers. Accès : Texte intégral.

**Narayanamurthi, Römer, Sandu · 2020-01-23.** [Goal-oriented a posteriori estimation of numerical errors in the solution of multiphysics systems](https://arxiv.org/html/2001.08824v1). arXiv 2001.08824v1. Accès : Texte intégral : cadre, adjoints et description des expériences.

**Peherstorfer, Willcox, Gunzburger · 2016.** [Optimal model management for multifidelity Monte Carlo estimation](https://cims.nyu.edu/~pehersto/). SIAM Journal on Scientific Computing 38(5), A3163–A3194. Accès : Résumé et bibliographie auteur ; PDF MIT et DOI inaccessibles pendant cette recherche.
