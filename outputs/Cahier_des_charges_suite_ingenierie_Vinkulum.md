# Cahier des charges — Suite d’ingénierie fondée sur Vinkulum

**Version :** 1.1 — interface graphique prioritaire  
**Date :** 9 septembre 2026  
**Nom de travail :** Suite d’ingénierie Vinkulum  
**Socle étudié :** dépôt privé `Brietat71/vinkulum`, branche `main`, version déclarée 0.18.0  
**Référence Git relevée :** `20803a2712efc4e50b9388bf5d7846cfa1660e6b`

## 1. Objet et portée du document

Ce cahier des charges définit une application d’ingénierie locale et collaborative permettant de concevoir, paramétrer, simuler et comparer des systèmes mécaniques, puis d’étendre ces usages à la CAO et à des couplages multiphysiques.

Vinkulum constitue le moteur mécanique initial. Le projet doit réutiliser ses capacités avant d’en développer de nouvelles. L’application ajoute un document métier persistant, des références stables, une validation asynchrone, une gestion des calculs et une interface utilisateur.

Le document définit une cible globale et une première livraison limitée. Il ne constitue ni une certification du moteur existant, ni une déclaration de disponibilité des fonctions futures. Les tests de Vinkulum n’ont pas été rejoués lors de sa rédaction ; les performances et qualifications mentionnées par le dépôt restent celles de ses dossiers de référence.

**Vocabulaire normatif :**

- **DOIT / NE DOIT PAS :** exigence obligatoire pour le lot concerné.
- **DEVRAIT :** préférence justifiée, avec dérogation documentée possible.
- **CIBLE :** objectif à mesurer et à qualifier avant engagement de livraison.
- **RECHERCHE :** expérimentation sans promesse de généralisation.

Les identifiants d’exigences doivent être repris dans les tickets, les tests et les rapports de recette. Une exigence d’un lot ultérieur ne bloque pas la livraison d’un lot antérieur, sauf dépendance explicitement indiquée.

### Priorité de livraison : interface graphique immédiatement utilisable

**Décision utilisateur : disposer rapidement d’une interface graphique, même temporaire.** Le premier livrable est donc une application de démonstration autour de l’API Vinkulum existante. Il ne doit pas attendre le modèle documentaire complet, la collaboration, la CAO, ni un nouvel évaluateur Rust.

La rapidité vient du périmètre réduit : un pendule, des champs numériques, un calcul séparé de l’interface et des résultats visibles. Les contrats de capture des entrées, de gestion des erreurs et de séparation des calculs restent obligatoires sous une forme minimale.

## 2. Finalité et utilisateurs

### 2.1 Finalité

Permettre à un ingénieur de passer d’un modèle paramétrique à un résultat physique traçable, sans perdre l’association entre hypothèses, géométrie, configuration numérique et résultats.

Le fonctionnement local est prioritaire : ouverture, édition, validation et calcul doivent rester possibles sans connexion. La collaboration synchronise les intentions de conception ; elle n’impose pas la réplication permanente de tous les états du solveur.

### 2.2 Utilisateurs et responsabilités

| Profil | Besoin principal | Capacités attendues |
|---|---|---|
| Concepteur | Construire et modifier un système | Paramètres, géométrie, composants, références et assemblages |
| Ingénieur calcul | Évaluer un comportement physique | Hypothèses, solveurs, tolérances, convergence et résultats |
| Collaborateur | Travailler sur un projet partagé | Édition concurrente, commentaires de diagnostic, résolution des ambiguïtés |
| Relecteur | Examiner un résultat sans altérer le modèle | Lecture du snapshot, provenance, diagnostics et export |
| Développeur scientifique | Ajouter ou qualifier une capacité | API, cas indépendants, tests, mesures et documentation du domaine |

Les profils sont des usages. Ils ne nécessitent pas cinq systèmes de permissions distincts. Pour le partage, trois droits suffisent initialement : propriétaire, éditeur et lecteur.

### 2.3 Parcours de référence

1. Créer un projet depuis un exemple ou ouvrir un projet existant.
2. Modifier des paramètres avec unités et expressions.
3. Observer les diagnostics et corriger les références ou contraintes invalides.
4. Choisir une analyse et ses paramètres numériques.
5. Lancer le calcul sur un instantané identifié et validé pour cette analyse.
6. Continuer à éditer sans modifier les entrées du calcul en cours.
7. Consulter les résultats, leur statut, leurs limites et leur provenance.
8. Comparer deux exécutions et exporter un dossier reproductible.

## 3. État de départ et réutilisation de Vinkulum

L’inspection du manifeste, du README et de portions du code confirme un noyau Rust avec API Python via PyO3. Le dépôt comprend des corps rigides, contraintes, poutres, calculs dynamiques et statiques, ainsi que des modules et dossiers d’analyse, de comparaison et de certification limitée.

`Modele` contient des données mécaniques et un état d’exécution mutable : temps, multiplicateurs, compteurs et caches notamment. `Noyau` expose ce modèle à Python. Ces objets ne doivent pas devenir directement le document collaboratif.

| Élément | Position retenue |
|---|---|
| Noyau mécanique Rust | Réutilisation prioritaire, corrections motivées par cas reproductibles |
| API Python et construction Maturin | Réutilisation ; compatibilité des usages existants à préserver |
| Cas analytiques et bancs du dépôt | Point de départ de la qualification, sans étendre leurs garanties |
| Certifications existantes | Présenter leur domaine précis ; ne jamais les résumer en « noyau certifié » |
| Document paramétrique et collaboration | Nouvelle couche ; absence à confirmer par inventaire ciblé avant ajout |
| CAO interactive et interface de bureau | Nouvelle intégration, non considérée comme disponible |
| Co-simulation généralisée et IGA | Lots futurs ou recherche selon le périmètre |

**BASE-01.** L’intégration DOIT conserver la possibilité d’utiliser Vinkulum indépendamment de l’application graphique et des services collaboratifs.

**BASE-02.** Toute modification du moteur DOIT respecter ses règles de contribution, de versionnement et de validation applicables au checkout utilisé.

**BASE-03.** Les défauts et limitations connus du dépôt DOIVENT rester visibles. Une nouvelle interface ne transforme pas un calcul hors domaine en calcul qualifié.

## 4. Périmètre fonctionnel et exclusions

### 4.1 Première livraison : interface temporaire G0

- Fenêtre PySide6 avec exemple de pendule préchargé depuis un cas existant de Vinkulum.
- Champs numériques avec unités affichées : paramètres effectivement exposés par le cas retenu, par exemple longueur, masse, angle initial, durée et pas de temps.
- Boutons Lancer et Arrêter, statut visible et message d’erreur compréhensible.
- Calcul dans un processus séparé ; les paramètres sont copiés au clic sur Lancer et rattachés au résultat.
- Animation simple du pendule après calcul, avec lecture/pause et déplacement dans le temps.
- Courbes de résultat, au minimum angle en fonction du temps ; autres courbes uniquement si disponibles sans nouvelle physique.
- Enregistrement et rechargement des paramètres dans un petit fichier JSON versionné ; export des échantillons avec unités.

Une représentation 2D du pendule est suffisante pour G0. Une vue 3D peut être ajoutée si une intégration existante la rend immédiate ; elle ne doit pas retarder la première fenêtre utilisable. Aucun streaming scientifique temps réel n’est requis : l’animation peut lire la trajectoire retournée par Vinkulum.

**G0 ne requiert pas :** expressions, graphe nodal générique, UUID pour chaque entité, CRDT, validation géométrique, historique complet, conteneur `.engproj` ou nouveau moteur de calcul. Un objet de paramètres indépendant par exécution suffit à figer les entrées.

### 4.2 Premier produit local structuré

- Projet local contenant un pendule mécanique fondé sur un exemple Vinkulum existant.
- Paramètres scalaires avec unités et expressions limitées.
- Références stables, détection des dépendances cycliques et diagnostics.
- Instantanés, validation asynchrone et exécutions indépendantes.
- Sauvegarde, réouverture, historique des calculs et courbes simples.
- Interface minimale permettant de modifier, lancer, observer et comprendre un échec.

La CLI sert aux tests et à l’automatisation ; elle ne remplace pas le premier livrable graphique demandé. Le produit local structuré étend ensuite G0 avec les exigences documentaires complètes.

### 4.3 Cible produit étendue

- Édition paramétrique de mécanismes rigides puis flexibles.
- Vue 3D, arbre de construction, propriétés et résultats synchronisés.
- Collaboration locale d’abord, puis partage entre machines avec travail hors ligne.
- CAO paramétrique et références topologiques explicites.
- Ateliers mécanique, structures, modal et analyses complémentaires qualifiées.
- Import de modèles et intégration de solveurs externes selon besoins démontrés.
- Calculs par lots, comparaisons et export scientifique.

### 4.4 Hors première livraison

CAO générale, contact complexe, CFD industrielle, IGA, réécriture des solveurs, P2P sans infrastructure, modèle universel de plugins, calcul distribué, boutique d’extensions, rendu 16K à 360 Hz et garanties bit à bit entre architectures.

Ces exclusions limitent le premier lot ; elles ne réduisent pas l’ambition globale. Aucun composant ne doit être ajouté uniquement pour anticiper un usage non décrit.

## 5. Architecture logique

La chaîne de responsabilité est : **document éditable → capture cohérente → validation → configuration d’exécution → instance de solveur → artefacts de résultats**.

| Composant | Responsabilité | Interdiction principale |
|---|---|---|
| Document métier | Intentions, paramètres, références et dépendances | Contenir les caches mutables du solveur |
| Réplication | Synchronisation des modifications autorisées | Déclarer le modèle physiquement valide |
| Capture | Produire des entrées cohérentes et immuables | Lire un assemblage de versions différentes |
| Validation | Établir les diagnostics pour une analyse donnée | Réparer silencieusement le document partagé |
| Adaptateur Vinkulum | Traduire des entrées validées en modèle de calcul | Réutiliser l’état mutable d’un autre calcul |
| Exécuteur | Lancer, surveiller, annuler et recueillir les sorties | Confondre succès du processus et acceptation physique |
| Worker géométrique | Reconstruire et résoudre les références CAO | Modifier le document ou publier des fichiers finaux directement |
| Stockage | Persister, vérifier et restaurer | Annoncer une sauvegarde avant sa réussite effective |
| Interface | Éditer et présenter l’état | Constituer le seul contrôle d’autorisation du lancement |

**ARCH-01.** La logique de document et de validation DOIT être testable sans interface graphique.

**ARCH-02.** Les calculs lourds et appels natifs bloquants NE DOIVENT PAS s’exécuter sur le thread d’interface.

**ARCH-03.** Une instance de calcul DOIT posséder son état mutable. Le partage de caches ne sera admis qu’avec un contrat d’immuabilité ou de synchronisation testé.

**ARCH-04.** Le processus de géométrie DOIT être isolé de l’application lorsque la CAO est introduite. L’exécution Vinkulum DEVRAIT également utiliser un worker pour permettre annulation, confinement d’un crash et maîtrise des ressources.

La séparation logique ne prescrit pas des microservices. G0 doit réutiliser l’API Python publique de Vinkulum et garder une structure minimale : interface, préparation des paramètres, worker de calcul et affichage. Le produit structuré peut ensuite regrouper document, validation légère, stockage et interface dans une application, avec un processus de calcul séparé.

## 6. Modèle de document

### 6.1 Objets persistants

| Objet | Champs ou relations minimaux |
|---|---|
| `Project` | Identité, version de schéma, titre, documents et paramètres de stockage |
| `Document` | Identité stable, objets métier, expressions, connexions et analyses déclarées |
| `Node` | Identité stable, type métier versionné, nom affiché, propriétés et ports |
| `Port` | Identité stable dans son contexte, direction, type, cardinalité et caractère requis |
| `Connection` | Identité, source et destination explicites |
| `Parameter` | Identité, nature physique, unité d’affichage, valeur littérale ou expression |
| `Frame` | Identité, origine, orientation, parent éventuel et conventions explicites |
| `TopologyReference` | Producteur, provenance ou sélection, cardinalité attendue et politique de résolution |
| `Asset` | Identité de contenu, taille, type et données ou référence locale résoluble |
| `AnalysisSpec` | Type d’analyse, sous-modèle ciblé, charges, hypothèses et options |

**DOC-01.** Les noms affichés et indices de tableaux NE DOIVENT PAS servir d’identité persistante. Renommer un objet ne doit pas rompre ses références.

**DOC-02.** Le modèle DOIT distinguer objets supprimés, références manquantes et types inconnus. Un objet inconnu à l’ouverture ne doit pas être supprimé silencieusement.

**DOC-03.** Chaque type de nœud DOIT déclarer ses entrées, sorties, contraintes et version de schéma. La migration d’un port référencé doit être explicite.

**DOC-04.** Les données dérivées — maillage généré, matrice, tessellation, trajectoire — DOIVENT référencer leurs entrées ; elles ne doivent pas être traitées comme paramètres collaboratifs par défaut.

### 6.2 Types et connexions

Le typage doit distinguer forme de donnée, nature physique et contexte. Un vecteur de force et un vecteur de position ne sont pas interchangeables. Deux vecteurs exprimés dans des repères différents nécessitent une transformation explicite ou une règle de conversion connue.

`Expression` désigne une représentation d’entrée, pas une nature physique. `TopoDS_Shape` reste un détail du worker OCCT ; le document utilise une référence géométrique indépendante de cette structure native.

**TYPE-01.** Une entrée simple NE DOIT PAS recevoir plusieurs sources sans règle métier définie. Les entrées multiples doivent déclarer leur ordre ou leur caractère non ordonné.

**TYPE-02.** Les valeurs non finies, références invalides et tailles incompatibles DOIVENT être refusées avant construction du modèle numérique.

### 6.3 Graphes

Le graphe de dépendances évaluables est construit à partir des connexions, références des expressions et dépendances de production géométrique. Il doit être acyclique dans le périmètre initial. Les expressions implicites simultanées sont hors périmètre initial.

Le graphe des contraintes mécaniques est distinct. Il peut comporter des boucles légitimes ; leur admissibilité relève des contrôles mécaniques et numériques.

**GRAPH-01.** Toute référence contribuant à une évaluation DOIT participer à l’analyse des dépendances. Les dépendances cachées dans du code arbitraire ne sont pas admises dans l’évaluateur initial.

**GRAPH-02.** Un cycle DOIT produire un diagnostic identifiant au moins un chemin cyclique et les objets concernés.

**GRAPH-03.** Le recalcul initial peut porter sur tout le sous-modèle requis. L’incrémentalité ne sera introduite qu’après mesure, avec preuve qu’elle produit les mêmes résultats que le recalcul complet sur les cas de recette.

## 7. Expressions, unités et conventions physiques

**UNIT-01.** L’évaluateur DOIT prendre en charge, au premier lot, constantes avec unités, références de paramètres, opérations arithmétiques et parenthèses. Toute fonction supplémentaire doit déclarer domaine et règles dimensionnelles.

**UNIT-02.** Les références de l’AST DOIVENT utiliser des identifiants stables. Le texte original peut être conservé pour l’édition ; son statut de source ou de représentation dérivée doit être unique et explicite.

**UNIT-03.** L’évaluation DOIT distinguer dimensions SI et nature de grandeur lorsque nécessaire : angle, énergie/couple, température absolue/écart notamment. Une égalité de dimensions n’autorise pas toutes les substitutions métier.

**UNIT-04.** Les conversions internes vers SI DOIVENT préserver l’unité choisie pour l’affichage. Les conversions affines ne doivent pas être traitées comme de simples facteurs multiplicatifs.

**UNIT-05.** Les règles de domaine DOIVENT être attachées au modèle physique utilisé : masse positive, inertie admissible, raideur et paramètres constitutifs selon formulation. Une règle d’un modèle ne doit pas devenir une interdiction universelle non justifiée.

**UNIT-06.** Les conventions de repère, rotation, ordre des composantes, signe des efforts et point d’application DOIVENT être définies aux frontières des adaptateurs.

**UNIT-07.** Le parseur DOIT borner taille, profondeur et coût d’évaluation. Aucun `eval` de code utilisateur général ne doit servir d’évaluateur d’expressions.

Recette minimale : conversion mm/m, incompatibilité longueur/température, renommage d’un paramètre, cycle uniquement dans les expressions, division par zéro, valeur non finie et expression trop profonde.

## 8. Révisions, instantanés et provenance

### 8.1 Identités distinctes

- **Fraîcheur locale :** `(DocumentID, SessionID, local_revision)`.
- **Identité persistante :** `SnapshotID`, unique pour un instantané enregistré.
- **Intégrité et déduplication :** empreintes des contenus et manifeste versionné.

Une révision locale n’est ni une horloge mondiale, ni une preuve d’égalité entre répliques. Un hash ne prouve ni validité physique ni provenance de confiance ; il identifie un contenu selon une procédure définie.

**SNAP-01.** Toute transaction métier effectivement appliquée, locale ou distante, DOIT invalider le statut courant et avancer la révision avant le debounce. Caméra, sélection et présence ne modifient pas cette révision.

**SNAP-02.** La capture DOIT associer atomiquement un contenu cohérent et sa révision : transaction de lecture, propriétaire unique du document ou copie protégée. Une édition pendant la capture ne doit pas produire un mélange.

**SNAP-03.** Un instantané DOIT figer ses dépendances nécessaires. Un chemin externe mutable seul est insuffisant. Les gros blobs peuvent être partagés par adresse de contenu, sans duplication à chaque capture.

**SNAP-04.** Les snapshots nécessaires à un calcul ou un rapport conservé NE DOIVENT PAS être supprimés par le nettoyage automatique.

**SNAP-05.** Si des empreintes sémantiques sont utilisées pour le cache, leur canonicalisation DOIT être spécifiée : ordre, nombres, unités, versions et dépendances. À défaut, hasher l’archive exacte et ne pas prétendre reconnaître tous les documents équivalents.

### 8.2 Séparation des états

| Objet | États minimaux |
|---|---|
| Validation | `Pending`, `Running`, `Valid`, `Invalid`, `Failed`, `Cancelled` |
| Exécution | `Queued`, `Running`, `Succeeded`, `Failed`, `Cancelled`, `Interrupted` |
| Issue numérique | `Converged`, `MaxIterReached`, `Diverged`, `NotAssessed` selon analyse |
| Résultats | `Partial`, `Complete`, `Unavailable` avec manifeste |

Les checkpoints sont des artefacts d’exécution ; ils ne remplacent pas le statut principal du calcul. Un calcul terminé peut conserver des avertissements ou sortir du domaine qualifié. L’interface doit distinguer ces situations.

## 9. Validateur métier asynchrone

### 9.1 Déroulement

1. Invalidation immédiate après transaction métier.
2. Regroupement des modifications rapprochées ; debounce initial proposé de 300 ms, configurable.
3. Capture cohérente de la version à contrôler.
4. Validation de structure et de types.
5. Validation des expressions, unités, références et règles métier.
6. Reconstruction géométrique si requise pour l’analyse.
7. Publication d’un rapport immuable lié au snapshot et à la configuration de validation.

**VAL-01.** Le validateur DOIT lire uniquement ses entrées figées. Les contrôles dépendant d’une bibliothèque ou d’un environnement doivent enregistrer leurs versions.

**VAL-02.** Un rapport ne peut mettre à jour le statut courant que si sa clé de fraîcheur correspond au document courant et si la demande de validation est encore pertinente.

**VAL-03.** La file automatique DOIT être bornée. Les demandes intermédiaires peuvent être remplacées par la plus récente. Un snapshot explicitement soumis à calcul doit rester identifiable et recevoir une issue.

**VAL-04.** L’annulation est une optimisation. La protection contre les rapports périmés DOIT rester correcte même si le travail ancien termine après son annulation demandée.

**VAL-05.** Un rapport DOIT contenir : identités, version du validateur, analyse ciblée, état de chaque étape requise, diagnostics typés et statut final. Absence d’erreur ne vaut pas achèvement.

**VAL-06.** Chaque diagnostic DOIT fournir code stable, objets concernés, gravité, caractère bloquant et explication. Une proposition de réparation doit devenir une transaction explicite.

**VAL-07.** Un calcul DOIT être refusé si une étape requise est incomplète, échouée ou non prise en charge. Un crash du validateur ne doit pas être présenté comme une erreur de conception.

**VAL-08.** L’autorisation de lancement DOIT être vérifiée dans la couche d’exécution, y compris pour CLI et scripts ; le bouton grisé n’est pas une garantie suffisante.

### 9.2 Tests de concurrence obligatoires

- Rapport V1 favorable reçu après rapport V2 défavorable.
- Rapport V1 reçu après édition V2 mais avant expiration du debounce V2.
- Réouverture du même document pendant une ancienne validation.
- Capture simultanée à une modification.
- Annulation non honorée par un worker, puis retour tardif.
- Rapport vide produit après interruption d’une étape obligatoire.

Ces tests doivent contrôler l’ordre des événements, sans dépendre de longues attentes réelles.

## 10. Exécution des calculs et intégration Vinkulum

**RUN-01.** L’adaptateur DOIT construire une nouvelle instance mécanique depuis le snapshot accepté. Les identifiants documentaires doivent être reliés explicitement aux indices internes du solveur.

**RUN-02.** Chaque exécution DOIT enregistrer snapshot, configuration numérique, version et empreinte des composants pertinents, plateforme et options influençant le calcul. Le binaire du solveur relève du manifeste d’exécution, pas obligatoirement du document métier.

**RUN-03.** Les modifications du document après lancement NE DOIVENT PAS modifier les entrées du calcul. L’interface doit indiquer si les résultats appartiennent à une version antérieure.

**RUN-04.** La file des calculs DOIT avoir une limite de concurrence configurable. Les budgets mémoire et de threads doivent tenir compte des bibliothèques numériques internes.

**RUN-05.** L’utilisateur DOIT pouvoir demander l’annulation. Le système doit proposer un arrêt forcé après délai si le worker ne répond pas, marquer les sorties partielles et conserver le diagnostic.

**RUN-06.** Les sorties DOIVENT être écrites dans un espace propre au calcul. Leur publication finale nécessite un manifeste cohérent ; une interruption ne doit pas faire passer un fichier incomplet pour un résultat complet.

**RUN-07.** Si la reprise par checkpoint est proposée, l’état sauvegardé DOIT couvrir les historiques nécessaires au solveur. La reprise doit être testée contre une exécution continue dans les tolérances annoncées.

**RUN-08.** Les exceptions, refus de paramètres, non-convergences et crashes DOIVENT rester distincts. Les résultats partiels peuvent être consultés mais doivent être identifiés comme tels.

Le premier cas est un pendule dans le domaine lisse admis par Vinkulum. Sa recette s’appuie sur un cas existant et une référence indépendante applicable au régime choisi. Il n’est pas demandé de réécrire ce solveur.

## 11. Interface et ateliers

### 11.1 Prototype graphique prioritaire G0

**GUI-01.** L’application DOIT ouvrir directement un exemple fonctionnel. L’utilisateur ne doit pas devoir créer un document nodal avant son premier calcul.

**GUI-02.** Les contrôles de G0 portent sur les champs exposés : valeurs finies, unités annoncées et domaine du cas retenu. Les validations existantes de Vinkulum restent actives.

**GUI-03.** Lancer DOIT capturer les entrées, créer une instance dédiée et démarrer un worker. Pendant ce calcul, toute modification des champs prépare le prochain lancement ; les paramètres de la trajectoire affichée restent consultables.

**GUI-04.** G0 peut limiter la concurrence à un calcul. Le bouton Lancer est alors désactivé pendant l’exécution ; l’édition et le bouton Arrêter restent disponibles.

**GUI-05.** Un worker qui échoue ou est arrêté NE DOIT PAS remplacer la dernière trajectoire réussie par un résultat incomplet présenté comme réussi. La trajectoire conservée doit rester clairement identifiée comme issue du calcul précédent.

**GUI-06.** L’animation DOIT offrir lecture/pause et sélection du temps. La géométrie affichée doit correspondre aux paramètres du calcul affiché, même si les champs ont changé depuis.

**GUI-07.** Le chargement JSON DOIT valider les champs et la version ; il ne doit exécuter aucun code. Les sorties exportées doivent inclure les paramètres du calcul et les unités.

**GUI-08.** La recette DOIT être réalisable en ouvrant la fenêtre : lancer l’exemple, voir l’animation et la courbe, modifier un paramètre, relancer, arrêter un calcul, enregistrer puis recharger les paramètres.

Choix proposé pour aller vite : PySide6 et API Python existante. Le dessin du pendule et une courbe simple peuvent utiliser Qt directement ; aucune bibliothèque 3D n’est obligatoire pour G0. Le code numérique reste dans Vinkulum. La future interface peut remplacer cette présentation sans réécrire le solveur ni perdre le format des paramètres documenté.

### 11.2 Organisation cible

L’interface cible comporte un explorateur de projet, un arbre de construction, un panneau de propriétés avec unités, une vue 3D, un panneau de diagnostics et un espace de résultats. Un éditeur nodal peut être ajouté lorsque les usages le justifient ; il n’est pas indispensable à la première livraison.

PySide6 est la préférence initiale pour l’application de bureau. VTK est un candidat pour la visualisation scientifique. La sélection du backend graphique doit être qualifiée sur les plateformes ciblées avant d’être figée.

### 11.3 Exigences

**UI-01.** L’interface DOIT distinguer brouillon, validation en cours, modèle invalide, validation indisponible et version validée.

**UI-02.** Les diagnostics DOIVENT permettre de retrouver l’objet concerné et expliquer l’action attendue sans afficher uniquement une exception technique.

**UI-03.** La navigation et l’édition courante DOIVENT rester disponibles pendant un calcul. Les champs en édition non encore soumis ne doivent pas être présentés comme déjà intégrés au snapshot courant.

**UI-04.** Une animation interpolée DOIT être distinguable des échantillons calculés. L’interpolation visuelle ne doit pas créer de nouvelles valeurs scientifiques exportées comme calculées.

**UI-05.** Courbes et tableaux DOIVENT afficher unités, repères pertinents, temps et identité du calcul. Les comparaisons doivent signaler les différences de modèles et configurations.

**UI-06.** Les fonctions essentielles DOIVENT être accessibles au clavier ; les états ne doivent pas reposer sur la couleur seule. Les panneaux doivent rester utilisables avec mise à l’échelle de l’affichage.

### 11.4 Ateliers cibles

| Atelier | Contenu | Condition d’introduction |
|---|---|---|
| Mécanismes | Corps, liaisons, efforts et commandes | Premier parcours local qualifié |
| Structures | Poutres, matériaux, charges et réponses | Domaines des formulations documentés |
| Analyse modale | Modes, fréquences et normalisation | Réutilisation des API existantes qualifiées |
| CAO | Esquisses, opérations, références et tessellation | Worker géométrique et contrat topologique qualifiés |
| Aérodynamique | Modèles réduits et échanges mécaniques | Domaine de chaque modèle et couplage établis |
| Thermique / fluides | Solveurs adaptés et transferts | Besoin utilisateur et adaptateur vérifiable |

L’ajout d’un atelier ne nécessite pas un système de plugins universel au préalable.

## 12. CAO et références topologiques

**CAD-01.** OCCT est le candidat initial pour la géométrie exacte. Build123d peut être évalué comme couche de construction, sans imposer son usage si l’accès existant suffit.

**CAD-02.** Les entités d’esquisse DOIVENT posséder des identifiants stables. Un réordonnancement ne doit pas modifier leur identité.

**CAD-03.** Une référence topologique DOIT décrire producteur, entités sources, relation de provenance et cardinalité attendue. Les noms affichés ne doivent pas être les clés de résolution.

**CAD-04.** Le résultat de résolution DOIT distinguer cible résolue, absente et ambiguë. Les évolutions de type division/fusion doivent être traitées selon une politique explicite, jamais par choix arbitraire de la première face.

**CAD-05.** L’intention de charge DOIT être distincte du sélecteur : pression, force totale, distribution et conservation attendue. Une réassociation géométrique ne doit pas décider implicitement de la physique.

**CAD-06.** L’intégration DOIT alimenter et tester l’historique topologique nécessaire à TNaming. L’existence de cette bibliothèque n’est pas une garantie autonome de persistance des références.

**CAD-07.** Le worker DOIT produire des artefacts distincts pour B-Rep et tessellation, identifiés par leurs entrées et paramètres. La tessellation d’affichage ne remplace pas la géométrie exacte.

**CAD-08.** Crash, délai dépassé et erreur d’opération DOIVENT produire des diagnostics distincts. Le parent doit rester utilisable et les écritures temporaires ne doivent pas altérer le projet.

Recette : extrusion simple, renommage et réordonnancement d’esquisse, suppression d’entité, division de face, fusion, changement amont, worker interrompu et worker bloqué.

## 13. Collaboration et fonctionnement hors ligne

**COL-01.** Le modèle partagé DOIT couvrir les intentions : paramètres, expressions, objets et connexions. Présence, caméra et curseurs relèvent d’un canal éphémère distinct.

**COL-02.** La représentation répliquée DOIT préciser la granularité des écritures. Pour une expression courte, le remplacement atomique du champ est acceptable au premier lot, avec règle de concurrence visible ; une fusion de fragments d’AST ne doit pas être supposée valide.

**COL-03.** Après reconnexion et réception de toutes les mises à jour nécessaires, les répliques DOIVENT converger vers le même document partagé. La validité métier doit être réévaluée localement.

**COL-04.** Les cycles, suppressions concurrentes et réécritures incompatibles DOIVENT produire des diagnostics. Le système doit proposer une édition corrective plutôt qu’effacer silencieusement une intention.

**COL-05.** L’annulation collaborative DOIT être spécifiée et testée : une annulation locale ne doit pas restaurer aveuglément l’ensemble du document en écrasant les modifications d’autrui.

**COL-06.** Les résultats volumineux sont partagés explicitement ou téléchargés à la demande. Leurs manifestes doivent identifier le calcul et les contenus disponibles.

**COL-07.** Les permissions de lecture et d’écriture DOIVENT être contrôlées au service de synchronisation. Une identité transportée par le client ne suffit pas à autoriser une modification.

**COL-08.** La reconnexion DOIT prévoir resynchronisation et reprise des assets interrompus. Toute opération persistante doit disposer d’une livraison fiable ou d’une récupération vérifiable.

Yrs/Yjs est un candidat de réplication à qualifier. Un relais central avec transport fiable est suffisant au premier lot distribué. QUIC, HTTP/3 et P2P ne sont pas des prérequis. La révocation d’accès doit empêcher les nouveaux accès autorisés par le service ; elle ne peut pas effacer les copies déjà détenues hors ligne.

## 14. Persistance, sauvegarde et formats

**STORE-01.** Le projet DOIT être utilisable depuis le disque local sans base serveur obligatoire.

**STORE-02.** Le format DOIT comporter une version de schéma, le document, un manifeste d’assets, les snapshots conservés et les références vers les calculs. L’extension de travail peut être `.engproj` ; le format physique sera choisi après le premier parcours.

**STORE-03.** Pour le prototype, un répertoire de projet avec manifeste lisible et fichiers de données est acceptable. Une archive transportable peut être ajoutée sans imposer sa réécriture complète à chaque sauvegarde.

**STORE-04.** La sauvegarde DOIT préserver la dernière version récupérable en cas d’interruption. La publication atomique d’un manifeste, les contrôles d’intégrité et le comportement sur disque plein doivent être testés sur chaque plateforme qualifiée.

**STORE-05.** L’ouverture DOIT vérifier schéma, tailles, références et intégrité des assets. Un asset manquant doit être signalé ; il ne doit pas être remplacé silencieusement par un fichier homonyme.

**STORE-06.** Toute migration DOIT conserver une copie récupérable de la version d’origine et fournir son résultat. Une version future inconnue doit être refusée proprement ou ouverte selon un mode limité explicite.

**STORE-07.** Le nettoyage DOIT respecter les références actives des snapshots, calculs et lecteurs. Il doit informer l’utilisateur avant suppression de résultats conservés explicitement.

**STORE-08.** Un export reproductible DOIT inclure les entrées nécessaires, la configuration d’exécution, les identités des composants et les limites de disponibilité des binaires ou dépendances.

Arrow, HDF5 ou un autre format scientifique seront choisis selon accès partiel, volumes, compression et compatibilité. PostgreSQL peut devenir un index de service ; Memgraph ne sera ajouté qu’en présence de requêtes mesurées qui le justifient. Les index serveur doivent rester distingués de la source documentaire.

## 15. Cohérence mémoire et flux de résultats

**MEM-01.** Aucun tampon publié NE DOIT être modifié ou réutilisé tant qu’un lecteur le détient. La durée de détention inclut les vues NumPy et objets de rendu.

**MEM-02.** Un échange atomique de pointeur ou un `Arc` ne suffit pas à démontrer cette propriété. Le protocole de propriété, acquisition et libération doit être testé.

**MEM-03.** Le flux visuel DOIT être borné et peut omettre des publications. Seules les trames en attente non détenues peuvent être remplacées ; une trame détenue ne doit jamais être écrasée.

**MEM-04.** Le flux scientifique NE DOIT PAS perdre silencieusement de données. Il peut ralentir le calcul ou échouer explicitement au plafond configuré, en conservant les sorties partielles identifiées.

**MEM-05.** Les limites mémoire DOIVENT être configurables et observables. Une croissance non bornée pour compenser un lecteur bloqué est interdite.

**MEM-06.** Le zéro-copie reste une optimisation locale. Il n’est pas exigé entre tous les processus, bibliothèques et appareils. Toute optimisation doit préserver durée de vie, cohérence et diagnostic des erreurs.

SoA, SIMD, alignement et parallélisme doivent être décidés sur les noyaux mesurés. Aucun changement global de disposition mémoire de Vinkulum n’est requis pour livrer la couche document.

## 16. Validité numérique et qualification scientifique

**NUM-01.** Chaque analyse DOIT déclarer son domaine : modèles physiques, conditions admissibles, approximations, tolérances, régimes exclus et références.

**NUM-02.** Les critères d’acceptation DOIVENT distinguer erreur de discrétisation, résidu numérique, erreur de modèle et comparaison expérimentale. Un petit résidu ne prouve pas une petite erreur physique.

**NUM-03.** Les normes mixtes DOIVENT utiliser des échelles explicites. Les seuils sont choisis avant mesure et justifiés pour le cas ; aucun `10⁻⁹` universel n’est prescrit.

**NUM-04.** La qualification DOIT inclure cas analytiques ou références indépendantes adaptées, convergence temporelle et spatiale lorsque pertinente, bilans et contre-exemples.

**NUM-05.** La reproductibilité initiale est définie par configuration enregistrée et tolérances. Une identité bit à bit ne doit être revendiquée que dans un périmètre testé et documenté.

**NUM-06.** Les résultats défavorables et refus DOIVENT être conservés dans les campagnes. Un refus de certification ne signifie pas automatiquement invalidité du résultat numérique.

**NUM-07.** Les analyses hors domaine qualifié DOIVENT être bloquées ou marquées expérimentales selon une politique explicite ; aucun avertissement ne doit disparaître à l’export.

**NUM-08.** Une comparaison de performance DOIT utiliser erreurs contrôlées, mêmes hypothèses, ressources décrites et coûts d’entrée/sortie explicités.

## 17. Co-simulation et extensions multiphysiques

### 17.1 Contrat d’adaptateur

Chaque participant doit déclarer variables échangées, unités, repères, signes, maillage d’interface, horloge, pas acceptés, état sauvegardable et modes d’échec. Les capacités de restauration et d’évaluation répétée ne doivent pas être supposées.

**COUP-01.** Commencer par deux sous-systèmes 1D masse-ressort-amortisseur avec référence indépendante. Qualifier échange, erreur énergétique et convergence avant couplage géométrique complexe.

**COUP-02.** Les interfaces DOIVENT expliciter conservation, stockage et dissipation. Pour une interface idéale sans stockage ni dissipation, le défaut de somme des travaux doit être mesuré ; une interface dissipative doit inclure son terme physique.

**COUP-03.** Les transferts DOIVENT préciser propriétés recherchées : conservation des efforts, cohérence des déplacements et compatibilité du travail discret selon formulation.

**COUP-04.** Un couplage implicite DOIT définir résidus, seuils absolus/relatifs, maximum d’itérations et stratégie d’échec. IQN-ILS est une option d’accélération, pas une garantie de stabilité.

**COUP-05.** Le rejet d’un pas DOIT restaurer tous les états nécessaires des participants avant nouvelle tentative. Nombre maximal de reprises et pas minimal doivent être bornés.

**COUP-06.** La séparation diagnostique entre intégration, transfert et couplage DOIT indiquer ce qui est mesuré, estimé ou non identifiable séparément ; aucun registre ne doit être présenté comme une borne sans justification.

### 17.2 Priorités d’intégration

Réutiliser les fonctions mécaniques Vinkulum qualifiées. Évaluer ensuite, selon cas métier, un adaptateur vers MBDyn/Chrono, un maillage Gmsh, des analyses CalculiX, puis aérodynamique réduite ou OpenFOAM. Pinocchio peut être évalué pour besoins de cinématique et dynamique inverse insuffisamment couverts.

Cette liste décrit des candidats, pas des dépendances imposées. preCICE doit être examiné avant de développer une infrastructure de couplage équivalente. Aucun solveur externe n’est programmé pour être remplacé uniquement en raison de son langage.

## 18. Sécurité, confidentialité et tolérance aux pannes

**SEC-01.** Les fichiers importés et mises à jour réseau DOIVENT être validés aux frontières : tailles, types, références, chemins, profondeur et coûts raisonnablement bornés.

**SEC-02.** L’ouverture d’un projet NE DOIT PAS exécuter automatiquement ses scripts. Le premier lot privilégie expressions et données déclaratives. L’exécution de scripts futurs nécessite une action explicite et un environnement aux capacités définies.

**SEC-03.** Un processus séparé confine certains crashes mais n’est pas, à lui seul, une sandbox de sécurité. Les workers de contenus non fiables doivent disposer d’accès fichiers et réseau limités selon leurs besoins.

**SEC-04.** Les workers DOIVENT écrire dans un espace temporaire dédié. Le parent valide les résultats et contrôle leur publication dans le projet.

**SEC-05.** Les communications partagées DOIVENT authentifier les participants et protéger le transport. Les secrets ne doivent pas être inscrits dans les snapshots ni les logs.

**SEC-06.** Les logs exportés DOIVENT être inspectables avant partage et ne doivent pas embarquer automatiquement des fichiers ou chemins sensibles inutiles.

**SEC-07.** Aucun changement de licence ou publication du dépôt n’est implicite dans ce projet. Les conditions de redistribution des composants retenus doivent être examinées avant livraison d’un installateur.

## 19. Performances et plateformes

Les valeurs ci-dessous sont des **cibles initiales de recette**, proposées pour rendre les mesures concrètes. Elles ne décrivent pas les performances déjà obtenues. La machine, le système, la version et le jeu de données de référence doivent être figés au lot 0.

| Mesure | Charge de référence proposée | Cible initiale |
|---|---|---|
| Retour visuel après édition | Document mécanique de 100 nœuds | p95 ≤ 100 ms, hors reconstruction et calcul |
| Validation structure/types | 1 000 nœuds, 5 000 dépendances, expressions bornées | p95 ≤ 200 ms, hors debounce et géométrie |
| Navigation 3D simple | Scène du pendule, viewport 1080p | 60 images/s visées sur machine qualifiée |
| Prise en compte annulation | Worker répondant au protocole | Accusé visible ≤ 1 s ; arrêt selon frontière sûre documentée |
| Mémoire sous lecteur bloqué | Publications répétées et lecteur volontairement retenu | Respect du plafond configuré, aucune croissance non bornée |
| Surcoût orchestration | Cas dont calcul domine le lancement | Mesurer capture, démarrage et sérialisation séparément ; budget fixé après baseline |

La performance d’un solveur n’est pas une propriété du seul nombre de corps. Les benchmarks doivent préciser contraintes, contacts, formulation, pas, durée et tolérances.

**PERF-01.** Les mesures DOIVENT présenter distributions ou répétitions pertinentes, pas seulement le meilleur temps.

**PERF-02.** L’application DOIT limiter la surallocation de threads entre processus, Rayon et bibliothèques numériques.

**PLAT-01.** La qualification initiale vise Linux x86_64 pour continuité avec les dossiers du dépôt, et macOS ARM64 pour l’usage de bureau ciblé. La disponibilité effective des dépendances doit être vérifiée. Windows constitue un lot de portage ultérieur, sans promesse préalable de compatibilité.

**PLAT-02.** Le paquet doit figer les versions compatibles de Rust, Python, bibliothèques natives et interface. Le dépôt inspecté exige Python ≥ 3.14 ; cette exigence ne doit pas être modifiée sans qualification.

## 20. Recherche indépendante

### 20.1 Résolution monolithique

Objectif : étudier les avantages d’une résolution unifiée sur un couplage défini. Le prototype doit expliciter inconnues, résidus, discrétisation, Jacobien, contacts et stratégie linéaire. Les choix direct/Krylov, Schur ou autres dépendent de la structure observée.

Critère de transfert produit : gain mesuré de robustesse, précision ou coût sur cas représentatifs, domaine explicite, tests de non-régression et coût de maintenance acceptable. « Tout Rust » n’est pas un critère suffisant.

### 20.2 Analyse isogéométrique

Commencer par Poisson ou élasticité plane sur géométries contrôlées. Comparer erreur à coût total donné, puis coût à erreur donnée, en documentant approximation géométrique et régularité de la solution. Tester raccordements, puis trimming dans des campagnes distinctes.

Le nombre de géométries du catalogue doit être justifié par sa couverture, pas fixé arbitrairement. Une victoire sur un cas ne justifie pas le remplacement d’OCCT. Une recherche négative reste un livrable utile.

### 20.3 Modèles réduits et apprentissage

Évaluer les modèles existants avant remplacement de backend. Documenter domaine d’entraînement ou de qualification, comportement hors domaine et dérivabilité requise. Un modèle rapide ne doit pas être présenté comme une référence physique indépendante.

Les recherches peuvent avancer indépendamment du produit. Elles ne deviennent des dépendances de livraison qu’après décision explicite et qualification.

## 21. Découpage de livraison et critères de sortie

| Lot | Livrable | Dépendances | Critère de sortie |
|---|---|---|---|
| L0 — Vérification ciblée | Checkout identifié, installation et cas pendule existant rejoué | Aucune | Cas retenu exécutable ; limites pertinentes connues ; aucun audit global préalable imposé |
| G0 — Interface immédiate | Fenêtre, paramètres, worker, animation, courbe, arrêt et JSON | L0 | Parcours GUI-08 réalisé ; erreur et annulation sans blocage de la fenêtre |
| L1 — Document et exécution structurés | Identités, expressions, snapshots, validation et Runs intégrés à G0 | G0 | Courses critiques refusées ; entrées et résultats traçables |
| L2 — Produit local enrichi | Historique, sauvegarde robuste, diagnostics et visualisation étendue | L1 | Parcours complet utilisable hors ligne et récupération après interruption |
| L3 — Collaboration | Deux clients, synchronisation, permissions, reconnexion, conflits | L2 | Convergence après partition ; aucun modèle invalide lancé par erreur |
| L4 — CAO et références | Worker OCCT, opérations limitées, sélection et charges | L2 | Références stables sur cas admis ; ambiguïtés et crashes explicites |
| L5 — Ateliers mécaniques étendus | Structures et analyses choisies dans les capacités Vinkulum | L2 ; L4 si géométrie requise | Recettes physiques et usages documentés par atelier |
| L6 — Couplage | Adaptateurs et cas d’interface limités | L1, capacité de checkpoint des participants | Bilans, convergence et rollback qualifiés |
| R — Recherche | Monolithe, IGA, modèles réduits et optimisations | Selon expérience | Rapport reproductible ; intégration conditionnelle |

L3, L4 et certains ateliers peuvent être menés indépendamment après stabilisation des contrats communs. Aucun délai ni budget n’est inventé : l’estimation sera produite après L0, selon équipe, plateformes et profondeur des ateliers retenus.

## 22. Recette et traçabilité des exigences

| Test | Scénario | Exigences principales | Résultat attendu |
|---|---|---|---|
| T00 | Ouvrir G0, calculer, animer, modifier, arrêter et recharger | GUI-01 à GUI-08 | Première fenêtre utilisable avec Vinkulum réel ; entrées du calcul inchangées |
| T01 | Renommer puis réordonner un objet référencé | DOC-01, CAD-02 | Références inchangées |
| T02 | Cycle créé uniquement dans deux expressions | GRAPH-01/02, UNIT-02 | Diagnostic de cycle ; lancement refusé |
| T03 | Rapport ancien reçu pendant debounce | SNAP-01, VAL-02 | Statut courant reste invalidé |
| T04 | Mutation pendant capture | SNAP-02 | Snapshot entièrement avant ou après transaction |
| T05 | Validation interrompue avec liste d’erreurs vide | VAL-05/07/08 | Lancement refusé |
| T06 | Modifier pendant simulation | RUN-01/03 | Entrées du Run inchangées ; résultats rattachés à l’ancienne version |
| T07 | Pendule et raffinement temporel | NUM-01/04, RUN-01 | Accord et convergence selon référence choisie |
| T08 | Arrêter worker puis rouvrir le projet | RUN-06/08, STORE-04 | Projet récupérable ; résultats partiels identifiés |
| T09 | Lecteur conserve un tampon longtemps | MEM-01/03/05 | Pas de mutation observable ; mémoire bornée |
| T10 | Coupure disque ou espace insuffisant | STORE-04 | Sauvegarde précédente récupérable ; échec explicite |
| T11 | Deux clients créent un cycle après fusion | COL-03/04, VAL-07 | Même document convergé ; calcul bloqué |
| T12 | Suppression concurrente d’un objet référencé | COL-04/05, DOC-02 | Référence manquante visible ; aucune cible substituée |
| T13 | Face divisée portant une force totale | CAD-04/05 | Politique explicite ou résolution requise |
| T14 | Crash et blocage du worker géométrique | CAD-08, SEC-04 | Hôte vivant ; statut Failed ; document conservé |
| T15 | Rejet puis reprise d’un pas couplé | COUP-04/05 | Restauration complète ; limites de reprise respectées |
| T16 | Fichier contenant script ou chemin sortant du projet | SEC-01/02/04 | Aucune exécution automatique ni écriture non autorisée |
| T17 | Installation dans environnement neuf | BASE-01, PLAT-02 | Exemple, import et tests requis exécutables hors checkout |
| T18 | Reprise après checkpoint | RUN-07 | Écart à exécution continue dans contrat documenté |

Chaque test doit identifier version, données, résultat attendu et artefact observé. Les tests numériques fixent leurs seuils avant exécution. Les tests de concurrence contrôlent l’ordonnancement. Les tests de performance utilisent un environnement isolé des campagnes concurrentes.

## 23. Qualité, exploitation et définition de terminé

**QUAL-01.** Une fonctionnalité livrée DOIT disposer d’un cas d’usage démontrable, d’un domaine déclaré, de diagnostics et d’au moins une vérification pertinente de ses risques.

**QUAL-02.** Les contrôles existants requis par Vinkulum DOIVENT être exécutés pour les changements qui le concernent. Les campagnes lourdes doivent être sélectionnées selon impact ; leurs absences et motifs sont consignés.

**QUAL-03.** Les journaux DOIVENT corréler document, snapshot, validation, worker et Run. Les erreurs doivent rester consultables après redémarrage lorsque le projet a pu être sauvegardé.

**QUAL-04.** Les formats publics, schémas de résultats et certificats DOIVENT être versionnés. Une rupture exige migration ou refus explicite, jamais interprétation silencieuse.

**QUAL-05.** Les agents de codage peuvent aider à implémenter et analyser, mais leurs affirmations ne remplacent pas une mesure, un test ou une preuve vérifiable. Les changements numériques requièrent revue de formulation et contre-épreuves adaptées.

Un lot est terminé lorsque ses exigences obligatoires sont vérifiées, ses défauts bloquants fermés, ses limites publiées et son installation reproduite sur les plateformes annoncées. Un rapport favorable sur un sous-ensemble ne clôt pas les autres obligations.

## 24. Risques principaux et réponses prévues

| Risque | Conséquence | Réponse |
|---|---|---|
| Confondre document et état de calcul | Résultats modifiés par édition | Adaptateur et instance dédiée par Run |
| Références géométriques ambiguës | Charge appliquée à mauvaise cible | Identités stables, provenance, résolution explicite |
| Validation périmée | Calcul lancé sur modèle non contrôlé | Révision immédiate, capture cohérente, contrôle côté exécuteur |
| Crash natif | Perte de session | Workers, écritures temporaires, diagnostic et récupération |
| Explosion mémoire | Blocage ou arrêt du système | Flux bornés, budgets et politique distincte visuel/scientifique |
| Réplication valide mais modèle incohérent | Résultats non exploitables | Validation métier après fusion |
| Extension excessive du périmètre | Retard sans parcours utilisable | Lots verticaux et critères de sortie |
| Garanties physiques exagérées | Mauvaise décision d’ingénierie | Domaines, références et limites présents dans l’interface et exports |
| Dépendances natives incompatibles | Installation impossible | Qualification plateforme dès L0 et installateurs testés |
| Couplage sans restauration complète | Trajectoire faussée après reprise | Contrat de checkpoint par participant |

## 25. Décisions retenues et décisions restantes

### Retenues par ce cahier des charges

- Vinkulum est le moteur mécanique initial ; pas de réécriture de principe.
- Le document métier est distinct du modèle numérique mutable.
- Le premier livrable est une interface graphique temporaire réutilisant un pendule existant ; il précède le modèle documentaire complet.
- Les calculs consomment des snapshots validés pour leur analyse.
- Le fonctionnement local précède les services distribués.
- Les ambiguïtés métier et erreurs techniques ont des statuts distincts.
- CAO native isolée, résultats persistants traçables, flux mémoire bornés.
- Recherche et objectifs produit possèdent des critères de sortie distincts.

### À trancher au moment où elles deviennent nécessaires

| Décision | Moment | Évidence attendue |
|---|---|---|
| Emplacement de la couche document dans le dépôt ou paquet associé | L0/L1 | Inventaire du code, maintien API et contraintes de build |
| Implémentation Rust du domaine et frontière Python | L1 | G0 réutilise Python ; introduction de Rust selon bénéfice et contrats du domaine, sans réécriture préalable |
| Format physique `.engproj` | L2 | Tests de sauvegarde, volumes et récupération |
| Backend 3D exact | L2 | Prototype sur plateformes qualifiées ; G0 peut utiliser un dessin 2D Qt |
| Granularité CRDT et transport | L3 | Scénarios concurrents, taille d’historique et reconnexion |
| API OCCT directe ou couche de construction | L4 | Opérations nécessaires et historique topologique disponible |
| Solveur externe à intégrer en premier | L6 | Cas utilisateur non couvert par Vinkulum |
| Budgets et calendrier | Après L0 | Capacité équipe, résultats de baseline et lots retenus |

Ces décisions ne nécessitent pas de suspendre les travaux indépendants déjà définis. La solution la plus simple satisfaisant les contrats doit être privilégiée.

## 26. Sources et limites de l’état des lieux

### Sources du projet

- [Dépôt Vinkulum](https://github.com/Brietat71/vinkulum)
- [Référence Git relevée](https://github.com/Brietat71/vinkulum/commit/20803a2712efc4e50b9388bf5d7846cfa1660e6b)
- [README : capacités et limites déclarées](https://github.com/Brietat71/vinkulum/blob/main/README.md)
- [Manifeste Rust](https://github.com/Brietat71/vinkulum/blob/main/Cargo.toml)
- [Manifeste Python](https://github.com/Brietat71/vinkulum/blob/main/pyproject.toml)
- [Modèle et API du noyau](https://github.com/Brietat71/vinkulum/blob/main/src/lib.rs)
- [Cas et contre-solveur Python](https://github.com/Brietat71/vinkulum/blob/main/python/vinkulum/maquette.py)
- [Plan de fiabilité](https://github.com/Brietat71/vinkulum/blob/main/docs/PLAN_FIABILITE.md)

Les fichiers ont été consultés sur `main` pendant l’inspection. La référence Git ci-dessus doit être figée et les contrôles rejoués au lot 0. Les liens `main` peuvent évoluer. L’état des lieux n’est pas un audit exhaustif du dépôt.

### Références techniques mobilisées dans la définition des contrats

- [Yjs : convergence sous réception des mises à jour](https://docs.yjs.dev/)
- [Rust : limites de catch_unwind](https://doc.rust-lang.org/std/panic/fn.catch_unwind.html)
- [Rust : propriété partagée avec Arc](https://doc.rust-lang.org/std/sync/struct.Arc.html)
- [Rayon : ordre des réductions et flottants](https://docs.rs/rayon/latest/rayon/iter/trait.ParallelIterator.html)
- [OCCT : mécanismes TNaming](https://dev.opencascade.org/doc/refman/html/package_tnaming.html)
- [BIPM : brochure du SI](https://www.bipm.org/fr/publications/si-brochure)
- [preCICE : checkpoint et couplage implicite](https://precice.org/couple-your-code-implicit-coupling)
- [preCICE : transferts entre maillages](https://precice.org/configuration-mapping)
- [preCICE : accélération du couplage](https://precice.org/configuration-acceleration)

Les exigences du présent document sont des choix de conception proposés pour ce projet. Les références techniques éclairent leurs contraintes ; elles ne certifient pas l’implémentation future.
