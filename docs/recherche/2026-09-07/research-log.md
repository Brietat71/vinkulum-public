# Journal de recherche — 7 septembre 2026

## Périmètre et décision

Étudier les avancées accessibles de 2001 au 7 septembre 2026 susceptibles
d'améliorer Vinkulum généraliste face aux références mondiales. Le
destinataire est la conception du projet ; le livrable retenu est un PDF
français avec priorités d'architecture, limites et expériences décisives.
Référence : paquet Vinkulum 0.7.2 ; documentation et prototype au commit
`723112b`. Sources primaires, prépublications datées, documentation officielle
et mesures propres sont distinguées. Aucun code de solveur concurrent lu.

L'outil `update_plan` a été appelé puis réessayé ; il était indisponible.
Le plan cadrage, découverte, vérifications, synthèse et contrôle a été
maintenu dans la conversation, avec une seule étape en cours.

## Recherches et lectures effectuées

Trois voies substantielles ont été confiées à trois chercheurs. Les familles
ci-dessous regroupent les requêtes et lectures ; elles ne constituent pas
une transcription exhaustive. Les sources et identifiants natifs sont
conservés dans `claim-source-ledger.json`.

1. **Géométrie et physique :** generalized-alpha sur Lie et modification
   sigma ; poutres SE(3), mixtes, gauchissement, B-splines ; réduction
   lagrangienne, ECSW, SSM ; énergie–moment, port-hamiltonien et ECCO.
2. **Contact et contraintes :** Coulomb et relaxations convexes ; proximal,
   ADMM, SCD semismooth*, champs irrotationnels, splitting 2026, Painlevé ;
   CCD par inclusion, trajectoires rigides, IPC, forces parasites ; saltation.
3. **Calcul et différentiation :** Newton–Krylov, LSMR, MINRES-QLP, graphe
   proximal, LQR ; globalisation, précision mixte, stabilité et randomisation ;
   différentiation implicite, Enzyme/Rust, préconditionnement appris, PINN.

La coordination a étudié l'état propre du projet et les sources officielles
MBDyn, Simpack, Adams, RecurDyn, Simcenter, Chrono, Exudyn, Drake,
MuJoCo/MJX/MJWarp et Newton. Elle a complété la recherche par les
estimateurs d'erreur GARK et la multi-fidélité.

Une première fusion a précédé la vague de vérifications ciblées. La
seconde vague a recherché hypothèses, contre-exemples, dates et différences
de lois, sans relancer une découverte générale.

## Matrice de lacunes après vérification

| Décision | Preuve et confiance | Limite ou contradiction | Suite utile |
|---|---|---|---|
| SO(3), AD ou GPU comme avantage | Plusieurs adoptions documentées ; élevée | Algorithmes propriétaires partiellement publics | Mesurer la combinaison ; ne pas inférer une absence |
| Modification sigma | Intégral auteur 2025 relu par la coordination | Gain partiel, initialisation, déjà dans Exudyn, temps Python/C++ mêlés | Même Rust, plusieurs familles, erreur commune |
| Poutre mixte | Théorie, prototype et archive vérifiés | Coût, instabilité à 10^16, statique seulement | Condensation stable puis dynamique |
| Rang et résolution par graphe | LSMR, MINRES-QLP, proximal et LQR | Symétrie, métrique, topologie ; réactions non uniques | Boucles, unités, duplication et résidu original |
| Réduction et SSM | Intégraux structurés/SSM, résumé ECSW | Préparation, garanties locales et domaine des charges | Hors construction, enrichissement, repli, amortissement |
| Coulomb exact et convexité | Sources primaires 2006–2026 | Sous-problème convexe et boucle extérieure distincts | Lois équivalentes, contrôle normale/tangente |
| Splitting juillet 2026 | Décomposition et convergence relues | Arche : 47 % sous 10^-6 au budget étudié ; incohérence texte/tableau non résolue | Prototype éventuel, aucune robustesse générale déduite |
| CCD et rotations | Inclusion et borne trajectoire–corde | Primitives et trajectoires limitées ; capsules à dériver | Rotation intermédiaire, rasance, petits jeux |
| Dérivées d'événements | Saltation et hypothèses explicites | Rasance/simultanéité hors garantie classique | Recalcul des impacts, dérivée directionnelle ou statut invalide |
| GPU et grands mécanismes | Guide MJWarp actuel relu | Débit différent de latence ; travail connecté encore prioritaire | Bancs mécanisme unique/lots ; aucun plafond universel à 60 ddl |
| Stabilité et AD | RSS 2021, Blondel 2022, Epperly 2026, Rust relus | Rang plein pour certains résultats ; AD Rust expérimentale | Produits et adjoints contrôlés, repli |
| Couplages | PH 2026, ECCO et GARK | Quadraticité, même schéma ; aucune garantie générale DAE avec impacts | Monolithique, énergie et phase, puis adaptation |
| Multi-fidélité et apprentissage | Résumé primaire MFMC ; DeepONet et contre-exemples PINN | MFMC intégral inaccessible ; hors-domaine et préparation | Résidu indépendant et coût amorti |
| Supériorité mondiale | Documentation officielle et banc MBDyn propre | Aucun autre moteur exécuté dans cette étude | Protocoles reproductibles et accès exécutable |

## Recoupements et arrêt

La coordination a relu les sources sigma 2025, Song 2026, PH 2026, ICF
version 2025, proximal RSS 2021, Blondel 2022, Epperly et ses coauteurs,
ainsi que les guides MJWarp et Rust. Les versions, niveaux d'accès et
dates sont conservés dans le registre de provenance.

ECSW, les B-splines de 2025 et MFMC n'ont pas été intégralement accessibles.
Les conclusions correspondantes restent au niveau consulté. Le travail
annoncé dans un numéro de novembre 2026 (DOI 10.1016/j.cma.2026.119227)
est exclu du corpus consolidé : disponibilité publique avant le 7 septembre
non vérifiée.

Arrêt : chaque décision dispose d'une preuve primaire ou d'une limite
explicite ; les contradictions susceptibles de changer l'architecture sont
résolues ou circonscrites. Les inconnues restantes portent principalement
sur le transfert à Vinkulum et exigent des expériences. Les trois voies
sont terminées ; la synthèse précède la création du rapport canonique.

## Production et contrôle

`report-source.md` est la source canonique. Les 58 références sont citées
près des affirmations et décrites en fin de rapport avec auteur, date,
publication, lien et niveau d'accès. Les identifiants natifs restent dans
le registre interne.

Le PDF vectoriel est produit sans réseau par `ci/rend_recherche.py`, avec
ReportLab 5.0.1 et DejaVu. Aucun graphique de données ne requiert un
tableur compagnon. Le premier rendu avait trois débordements vers des
pages presque vides. Le contrôle structurel et les vues des pages denses
ont motivé une synthèse plus courte, des colonnes adaptées et une page
propre au protocole. Le rendu corrigé compte 18 pages ; les métadonnées
portent la date éditoriale du 7 septembre 2026.

Les empreintes et résultats du contrôle final figurent dans
`verification.json`. La version du paquet et le tag 0.7.2 restent
inchangés : ce livrable documente la recherche et ses décisions.
