# Objectif : faire considérablement mieux que MBDyn, Simpack et les autres références

Objectif actif. La [confrontation statique 0.8.2](CONFRONTATION_EXUDYN_0.8.2.md)
compare désormais Vinkulum, MBDyn et Exudyn 1.11.0 sur Princeton ;
la [confrontation 0.7.2](CONFRONTATION_MBDYN_0.7.2.md)
reste la référence des dynamiques. Celle de la 0.6.2 reste le point de départ
archivé, avec les paramètres, les références et les échecs de chaque moteur.
Aucun résultat sur une seule famille ne vaut achèvement de cet objectif.

La [confrontation harmonique 0.10.0](CONFRONTATION_PORTS_EXUDYN_0.10.0.md)
ajoute une mesure des réductions sur les mêmes facteurs matériels,
avec HCB officiel Exudyn et deux témoins LU. Le juge contrôle les
champs, déformations et ports, y compris les combinaisons de charges.
Ce front expose le coût de la préparation et les refus de bande de
Vinkulum ; il ne remplace pas les comparaisons statiques et dynamiques.
Le [nouveau carnet de preuve](TRACE_COMPLEMENT_SPECTRAL_PREUVES.md)
étudie une élimination du complément des directions résonantes,
avec trace contrainte et obligations d'encadrement numérique.

La [livraison 0.11.0](REPONSES_GROUPEES_0.11.0.md) supprime une partie
du coût des bornes répétées : six charges groupées partagent leurs calculs
à fréquence fixée. Sur 128 poutres à 20 Hz, son API avec bornes est 5,7 fois
plus rapide que la 0.10.0 et 4,5 fois plus rapide que le témoin HCB retenu,
préparation comprise. La LU corrigée reste plus rapide sur les six cas
acceptés par Vinkulum. Le coût du certificat et la restriction de bande
restent à traiter dans l'API publique.

Le [complément spectral expérimental](COMPLEMENT_SPECTRAL_DIRIGE.md) possède
maintenant un certificat à arrondis dirigés. Une direction intérieure
conservée porte la limite certifiée du complément à environ 65,7 Hz sur
512 poutres. Le témoin KKT passe les 257 fréquences de 0–40 Hz sur les trois
maillages, toutes combinaisons des six charges comprises. Il résout encore
un grand système à chaque fréquence : aucune avance de vitesse n'en découle.
Le [prototype Krylov contraint](KRYLOV_CONTRAINT_PROTOTYPE.md) réalise
cette réduction avec enveloppe de résidu, marge sur tout le bloc conservé
et 31, 39 et 47 inconnues. La nouvelle
[campagne par inertie dirigée](INERTIE_CONTRAINTE_PROTOTYPE.md) contient
84 essais : les 24 avec inertie et les 24 avec trace passent sur 0–40 Hz ;
les 24 avec contrôle passent aussi leurs majorants. Sur 512 poutres,
le certificat passe d'environ 0,685 à 0,199 s. Les champs deviennent
plus rapides que la LU au total : **0,451 contre 0,711 s**. Avec contrôle,
le total descend de 1,305 à **0,816 s**, encore derrière la LU.
À 128 poutres, ce dernier service devance HCB énergie (0,298 contre
1,627 s), mais reste derrière la LU (0,222 s). Les neuf refus HCB sont
conservés ; leurs configurations sont exclues des ratios de vitesse.
Le [contrôle par facteurs](CONTROLE_FACTEURS_PROTOTYPE.md) traite ensuite
les normes d'extension et les images de réparation. Ses 60 essais conservent
24/24 contrôles admis et dix refus HCB. Sur 512 poutres, le total avec
contrôle descend de **0,819 à 0,671 s**, légèrement devant la LU du même
lot (**0,696 s**, environ 3,6 % de temps en moins). Les champs et leurs
normes physiques sont identiques entre les deux contrôleurs. À 128 poutres,
le total de **0,253 s** devance HCB énergie (**1,654 s**), mais reste
derrière la LU (**0,222 s**). Les petits cas et le coût de préparation
restent à améliorer. L'assemblage et le transfert dans l'API ont ensuite
été réalisés, comme indiqué ci-dessous.
La [campagne suivante de 36 essais](INERTIE_BINAIRE_SEPARATEURS.md) remplace
le certificat Decimal par des congruences binary64 compilées, après dissection
du graphe physique. Tous les essais passent, et les 24 contrôles conservent
champs, normes et majorants identiques. À 512 poutres, le total devient
**0,516 contre 0,673 s** pour le témoin Decimal et **0,715 s pour la LU**.
Le gain sur la LU est de 27,8 % ; à 128, il est de 4,1 %. À 32, le prototype
reste 39,6 % plus lent. La permutation est comprise dans ces durées. Aucun
nouveau classement Exudyn/MBDyn/Simpack ni résultat mémoire commun n'est établi.
Les réponses restent évaluées en doubles, sans certificat machine de leurs arrondis.
Le [transport des normes par comparaison massique](MASSE_COMPAREE_PREUVES.md)
retire ensuite la limite de six coordonnées par composante connexe de masse
du contrôleur expérimental. Les couplages physiques restent dans les équations.
Sur la nouvelle masse consistante, les 24 champs et 12 contrôles passent ;
le total à 512 éléments est **0,559 contre 0,729 s pour la LU corrigée**.
À 128, l'écart est inférieur à 1 % ; à 32, le prototype reste plus lent.
Les certificats de comparaison sont rejoués sur les données exactes stockées.
Cette campagne expérimentale ne renouvelle aucun classement face aux
moteurs multicorps externes ; son transfert dans le paquet arrive en 0.12.0.
L'[assemblage des compléments contraints](ASSEMBLAGE_COMPLEMENTS_PREUVES.md)
relie désormais plusieurs sous-structures tout en gardant leurs modes
intérieurs privés. La marge porte sur le bloc global complet. Le nouveau
lot de trois branches orientées, avec bras de levier et masses consistantes,
admet 24/24 champs et 12/12 contrôles. À 512 éléments par branche,
le total est **1,665 s contre 2,288 s pour la LU corrigée** ; à 128,
**0,657 contre 0,722 s**. Le petit cas reste plus lent. La normalisation
par somme des métriques de port limite la perte des bornes lors de la
composition. Les raccordements mobiles et les grands blocs d'interface
restent ouverts. La [0.12.0](VERSION_0.12.0.md) rend la réduction contrainte,
les certificats natifs et l'assemblage accessibles dans une roue installable,
avec contrôle supplémentaire de la masse complète. Cette livraison
ne transforme pas les majorants numériques des réponses en certificats.
L'[audit des résultats creux de 2026](VERIFICATION_CREUSE_2026.md)
et la [preuve d'inertie dirigée](INERTIE_DIRIGEE_PREUVES.md) étayent cette
voie de certification, avec hypothèses et contre-exemples explicites.

L'objectif couvre tous les tableaux grâce aux avancées mathématiques et
algorithmiques : précision, robustesse, temps, mémoire, couverture et usage.
L'exigence explicite est d'utiliser massivement les avancées en
mathématiques fondamentales, en informatique de pointe et en physique.
Chaque méthode doit préciser ses hypothèses,
les invariants qu'elle conserve et les contrôles qui confrontent sa théorie
au calcul flottant et aux modèles physiques.
Le [point de recherche mathématique](RECHERCHE_MATHEMATIQUE.md) distingue
les méthodes déjà livrées des pistes récentes étudiées.
L'[approfondissement des fondements mathématiques](FONDEMENTS_MATHEMATIQUES_VINKULUM_2026.pdf)
réoriente le programme vers la dimension effective, la largeur des couplages
et la stabilité : élimination locale avec contrôle du rang, réduction des
interfaces avec contrôle d'erreur, formulation mixte et préconditionnement
robustes aux paramètres. Ses 21 sources et obligations de preuve ne
constituent ni une nouvelle capacité livrée ni un classement de performances.
L'[étude scientifique 2001–2026](ETUDE_SCIENTIFIQUE_VINKULUM_2001_2026.pdf)
hiérarchise désormais quatre axes : structure des contraintes et opérateurs
creux, flexibles et réduction préservant la structure, contacts couplés et
événements, intégration géométrique et échanges d'énergie.
Simpack fait désormais explicitement partie des références à confronter.
Ses analyses statiques, temporelles et fréquentielles, ses fonctions de temps
réel et son environnement de modélisation sont décrits par
[Dassault Systèmes](https://www.3ds.com/products/simulia/simpack/core).
Ses formulations de corps flexibles, ses contacts déformables et ses
interfaces avec les éléments finis doivent également entrer dans la
[comparaison de couverture](https://www.3ds.com/products/simulia/simpack/flexible-body).
Ces descriptions de l'éditeur ne constituent pas des mesures comparatives.
Aucune exécution Simpack n'est encore archivée dans ce dépôt. Les protocoles
doivent comparer des problèmes physiques équivalents, avec références,
tolérances, sorties et coûts de préparation explicités pour chaque moteur.

## Opérateurs adjoints : version 0.9.0

La [livraison 0.9.0](OPERATEURS_ADJOINTS.md) supprime les matrices de
transition denses du recul temporel et fournit les tangentes CSC et les
produits de dérivées des poutres. Le gradient partagé à 120 poutres gagne
**32,47×** face à la roue 0.8.2 avec l'interface existante. Le gradient
des 120 rigidités gagne **43,61×** avec les produits locaux et la dérivée
initiale exacte du modèle de contrôle ; le pic RSS passe de 335 à 77 Mio.

La campagne finale conserve 156 essais, en plus d'un premier lot distinct,
et les tests confrontent les gradients aux différences finies. Le produit
isolé à 30 poutres ralentit à cause de l'import de SciPy ; les limitations
aux contraintes holonomes lisses, aux historiques et au stockage de
trajectoire sont explicites. Les calculs modaux restent denses. Ce résultat
interne avance le front optimisation sans établir un classement externe
face aux adjoints d'Exudyn, MBDyn ou Simpack.

## Première application : cinématique σ de 2025

La [version 0.8.0](INTEGRATION_SIGMA.md) ajoute une option σ au noyau,
avec élimination locale et tangente implicite. Sept familles mesurées
montrent un gain ciblé sur l'orientation de la toupie, mais une dégradation
des réactions au même pas et aucun gain général sur les flexibles ou
les boucles. Le schéma historique reste le défaut. Cette expérience
concrétise une piste de la recherche ; elle ne clôt aucun des fronts
ci-dessous et ne remplace pas les comparaisons externes.

## Fronts à couvrir

| Front | Preuve attendue | État au début du travail |
|---|---|---|
| Précision | Convergence temporelle et spatiale, références indépendantes, réactions et événements | Bons accords sélectionnés ; poutres à raffiner davantage que les `beam3` |
| Robustesse | Modèles difficiles terminés avec équilibre et contraintes vérifiés ; restaurations après échec | Princeton statique en échec ; initialisations difficiles des deux moteurs |
| Temps de calcul | Comparaisons répétées au même seuil d'erreur, coûts du noyau et de l'usage distingués | Avantages partagés selon le problème et le réglage |
| Temps réel | Latence par pas, maximum et quantiles, échéances manquées, plateformes et erreurs documentées | Aucune supériorité démontrée sur les fonctions temps réel de Simpack |
| Mémoire et échelle | Pics RSS et pente en taille, familles de topologies différentes | Coût des trajectoires Python et des assemblages globaux à réduire |
| Dérivées et optimisation | Gradients indépendamment vérifiés et coût selon le nombre de paramètres, limites aux discontinuités | Sensibilités et adjoints présents ; comparaison externe à effectuer |
| Couverture généraliste | Modèles vérifiés en robotique, transmissions, structures, contacts et multiphysique | Couverture multiphysique moins étendue que MBDyn |
| Usage et confiance | Modèles reproductibles, diagnostics, documentation et utilisateurs indépendants | API Python et CI disponibles ; validation externe encore limitée |

## Globalisation de Newton : version 0.8.1

Le [travail de globalisation](GLOBALISATION_STATIQUE.md) réduit les cinquante
paliers Princeton de 879 à 334 évaluations principales et de 100 à
50 tentatives. Une contraction des corrections prédites autorise certains
pas refusés par les forces, avec un filtre séparé de fermeture des liaisons.
Les six chaînettes du diagnostic retrouvent leur convergence en huit
évaluations, y compris les tailles que la variante sans filtre faisait
diverger. Les critères physiques et les restaurations sont conservés.

Le corpus de comparaison couvre aussi des directions de charge, des
changements d'unités, des modèles tournés, des boucles redondantes et des
contacts. Les 306 résolutions chronométrées montrent un gain de 3.17 à
3.29 fois sur les rampes, entre roues Vinkulum isolées. Deux réglages
resserrés à petite échelle restent en échec dans
les deux versions ; le candidat peut consacrer davantage d'itérations à
leur refus, avec environ 61 % de temps supplémentaire. Le résultat ne
clôt donc ni le travail de robustesse, ni celui des tolérances. La mise
à jour de la confrontation externe est maintenant exécutée : 150 calculs,
3 800 paliers Vinkulum stricts, aucun échec. Au seuil de 10 µm au bout,
MBDyn reste 44.6 fois plus rapide que la formulation milieu et 11.4 fois
plus rapide que l'intégrée, avec moins de mémoire dans ce protocole.
Les écarts se réduisent mais restent importants. Les temps complets ne
doivent pas être confondus avec les temps internes de globalisation.

## Efforts analytiques et tangente des poutres : version 0.8.2

Le [gradient par travail virtuel](TANGENTE_POUTRE_ANALYTIQUE.md) supprime
les différences finies du bloc de rotation. Les six rampes Princeton
prennent 1,98 à 2,56 fois moins de temps que la roue 0.8.1, avec les
mêmes équilibres et les mêmes deux refus connus à petite échelle. Les
sensibilités aux rigidités partagent les nouvelles expressions.

Le bénéfice ne s'étend pas automatiquement à tous les appels : l'analyse
complète `k_c_m_z` ne gagne que 1,04–1,09 fois à 120 éléments. Le coût des
matrices et projections globales doit être traité. La
[lecture critique de la stack scientifique](COMPOSABILITE_SCIENTIFIQUE.md)
renforce les priorités d'opérateurs creux, de réduction avec domaine
vérifié et de sensibilités aux contacts sous hypothèses explicites.
La confrontation externe est maintenant actualisée : 232 calculs avec
MBDyn et Exudyn 1.11.0, Newton complet/modifié et binaires standard/fast
pour ce dernier. Au seuil commun de 10 µm au bout, l'intégrée Vinkulum
prend 131,10 ms contre 271,62 ms pour le meilleur Exudyn testé et
15,88 ms pour MBDyn. Le défaut milieu reste plus lent qu'Exudyn.

Exudyn devient une référence explicite pour les flexibles réduits,
les modes creux et l'intégration géométrique : sa documentation propose
déjà la cinématique σ de 2025 et une chaîne FFRF/Hurty–Craig–Bampton.
Le [rapport comparatif](CONFRONTATION_EXUDYN_0.8.2.md) distingue ces
capacités documentées du seul résultat statique exécuté. Ni les contacts,
ni la dynamique Exudyn, ni Simpack ne sont classés par cette campagne.

## Confrontation de la roue 0.7.2

La campagne exécute 326 calculs, dont 89 échauffements, et conserve dix
sondages préalables. Le juge vérifie les deux références par moteur, les
trois répétitions des candidats et chaque palier statique strict.

Vinkulum prend 15,3 % de temps en moins sur le mécanisme plan ; les temps
spatiaux sont proches ; MBDyn est 1,86× plus rapide sur Andrews. Sur
Princeton, les 3 800 paliers Vinkulum passent, mais MBDyn reste 109,8× plus
rapide face à la poutre milieu et 20,8× face à l'option intégrée, au seuil
de 10 µm. Les temps complets incluent les runtimes et les politiques de
sortie propres à chaque outil. Les écarts de mémoire restent défavorables
à Vinkulum dans ce protocole.

Une priorité ressort des diagnostics : les deux configurations Princeton
retenues consomment chacune 879 évaluations principales pour 50 paliers.
Chaque palier échoue avec le mérite en forces puis réussit avec le mérite
en corrections. Il faut réduire ce travail répété tout en conservant les
garanties de restauration, les tolérances physiques et les modèles qui
bénéficient du mérite en forces. La formulation flexible et les coûts de
sortie restent à améliorer. Ces constats ne ferment aucun des autres
fronts, notamment le temps réel et la confrontation exécutée à Simpack.

## Prototype de formulation mixte de 2026

Le [prototype mixte de 2026](POUTRE_MIXTE_PROTOTYPE.md) ouvre maintenant un
travail de formulation indépendant de ces correctifs. Les ordres 1 et 2
ont été dérivés, exécutés et comparés à une référence continue de Cosserat.
Le gain de convergence spatiale est mesuré ; le temps et la condensation
aux rapports de raideurs extrêmes restent défavorables. Les sources et
les échecs sont archivés. Ce prototype n'ajoute aucune capacité à la roue
0.7.2 ; la dynamique et les couplages généralistes restent à construire.

## Premier défaut traité : statique anisotrope

`statique()` essaie toujours le mérite en forces. En cas d'échec du palier,
il restaure l'état initial de ce palier puis essaie un mérite fondé sur la
correction de pose prédite par la même factorisation. Les conditions finales
sur l'équilibre et les contraintes restent celles du solveur.

L'abandon prématuré fondé sur la seule évolution des forces n'est pas
appliqué à cette seconde stratégie. Le budget d'itérations reste borné ;
la continuation en charge et la restauration en cas d'échec sont conservées.

La substitution globale du premier mérite avait dégradé une chaînette.
La stratégie de reprise a donc été vérifiée contre cette régression, en
plus de Princeton. Un contrôle négatif rejoue le nouveau test avec le binaire
antérieur : il échoue à la charge pleine, comme la confrontation le prédit.

Les mesures finales de cette étape sont dans le
[rapport du correctif](CORRECTIF_STATIQUE_PRINCETON.md).

## Assemblage local des poutres réalisé

[Mesures et contrôles](RAIDEUR_LOCALE.md) : gain de 6,76× sur Princeton à
60 intervalles, sans changement de précision. Les matrices d'analyse sont
également accélérées. Les tangentes locales des poutres remplacent leurs
réévaluations globales ; à cette étape, les autres contributions restaient
globales.

## Résolution statique optimisée

La [factorisation creuse en statique](STATIQUE_CREUSE.md) réduit ensuite le
temps à 60 intervalles de 0,406 s à 0,231 s, à précision inchangée. Le gain
cumulé depuis le premier correctif est de 11,87×. La résolution vérifie son
résidu et conserve un repli pour les systèmes singuliers. Les 25 tests
Rust, 44 tests Python, 41 groupes de vérification, 46 bancs rapides et
9 bancs de contact passent. Dans le profil à 60 intervalles, la raideur
représente désormais environ 71 % du temps, la résolution linéaire 8 %.

## Nouvelle formulation de poutre : gain de précision, limite identifiée

L'option [poutre intégrée](POUTRE_INTEGREE.md) reproduit la raideur linéaire
de Timoshenko et l'arc sous moment pur. Sur Princeton à dix intervalles,
l'erreur au bout passe de 254,87 à 5,33 µm. Sous le seuil de 7,1 µm, les
configurations mesurées passent de 60 intervalles / 0,2322 s à
10 intervalles / 0,0266 s, soit 8,73× sur ce critère, entre formulations
de Vinkulum. Cela ne classe pas encore Vinkulum et MBDyn.

La condensation modifie les sensibilités : leur règle de chaîne et
l'adjoint dynamique ont été contrôlés. La campagne générale passe avec
47 tests Python et les nouveaux cas spécifiques à l'option intégrée.

**Le remplacement général n'est pas acquis.** Sur la chaîne de vingt
éléments du corpus, avec EI presque nul et rotations nodales bloquées,
l'option dépasse le délai de vingt secondes tandis que l'élément historique
termine. Le défaut `formulation="milieu"` est donc conservé et
`formulation="integree"` reste explicite et expérimentale. Le prochain
travail sur cette option doit traiter la flexion interne non linéaire et
ce cas limite, sans faire disparaître son résultat de la campagne.

## Tangente mixte exacte et diagnostic du câble

Les [blocs de translation et de couplage](TANGENTE_MIXTE.md) sont désormais
calculés sans différences finies, pour les deux formulations. Le bloc entre
rotations reste approché. Le temps statique Princeton diminue de 18 à 25 %
à précision spatiale inchangée. Les 26 tests Rust, 47 tests Python et les
campagnes complètes passent.

Le diagnostic montre que le modèle du câble bloque les rotations : son
équilibre en translation est linéaire et la formulation intégrée y prédit
une flèche de 320 km, avec un conditionnement de 8,72e10. Ce n'est pas une
caténaire physique. Newton échoue encore depuis l'état original dans le
budget borné de huit itérations.

**Défaut identifié à cette étape : l'échelle du critère statique.** Depuis la
solution linéaire calculée séparément, le solveur annonce un succès mais
normalise par 1,24e8, contre 14,715 au départ original. Le résidu de 4,12 µN
ne tient pas la tolérance rapportée à la première échelle. Le diagnostic
enregistre désormais ce contrôle à échelle fixe ; ce résultat ne peut donc
pas servir à déclarer le câble résolu. Il faut corriger cette dépendance et
rejouer les cas de précontrainte et de charges faibles, puis poursuivre la
validation physique de la formulation intégrée.

## Réactions exclues de la normalisation statique

Le [correctif d'échelle libre](ECHELLE_STATIQUE_LIBRE.md) normalise désormais
par le résidu initial hors contraintes. Sur le test analytique du ressort,
une charge de 1e12 N reprise par le guide ne masque plus la force libre de
1 N : la flèche passe du faux résultat nul à la valeur exacte de 0,1 m.
Le test échoue sur l'ancien binaire et vérifie aussi la restauration après
un refus pour budget insuffisant.

Le câble extrême ne bénéficie plus du succès trompeur depuis sa référence
linéaire ; il reste non résolu dans le budget du diagnostic. Les 26 tests
Rust, 48 tests Python, 41 groupes de vérification, 46 bancs rapides et
9 bancs de contact passent. Princeton à charge pleine conserve exactement
ses positions finales ; le chemin en paliers peut demander davantage de
travail avec le critère resserré, ce qui reste publié.

**Limite identifiée du contrat de convergence :** le critère reste relatif
au déséquilibre libre initial, avec plancher 1, et le contrôle de stagnation
peut encore accepter sous 1e-6 fois cette échelle sans tenir `tol`.
Il faut distinguer ce statut d'une convergence à la tolérance demandée et
clarifier les contrôles absolus, avant de revendiquer une robustesse
générale supérieure. La validation physique de l'option intégrée reste
également ouverte.

## Verdict statique explicite et mode strict

Le [rapport de convergence](STATUT_STATIQUE.md) distingue désormais
`tolerance`, `stagnation` et `echec`. Le secours historique émet un
`RuntimeWarning` ; `strict=True` le refuse. Les reprises, les résidus et
la restauration sont consultables par `statique_info()`. Un avertissement
transformé en exception restaure également le modèle, y compris après un
déplacement effectif.

Le témoin sans équilibre, auparavant accepté sans avertissement, est
maintenant signalé comme stagnation et refusé en mode strict. Sur Princeton
à quarante intervalles, le parcours en cinquante paliers utilise ce secours
49 fois à `tol=1e-8` ; le mode strict refuse dès le premier palier, avec un
résidu relatif de 1,90e-8. La charge pleine passe dans les deux modes.
Les diagnostics conservent ces résultats distincts.

Les 26 tests Rust, 50 tests Python, 41 groupes de vérification, 46 bancs
rapides et 9 bancs de contact passent ; l'API comporte 184 entrées. La
campagne générale utilise encore le mode habituel et ne prouve donc pas
la convergence stricte de tous ses modèles.

**Travail suivant :** améliorer la précision et la convergence des cas
stricts refusés, et préciser les contrôles absolus et les échelles de
grandeurs. Le nouveau rapport améliore le contrat et sa vérifiabilité ;
il ne résout ni ces cas, ni la validité physique du câble condensé, ni les
autres fronts de comparaison avec MBDyn.

## Précision des positions : cause isolée, intégration encore ouverte

Le [diagnostic numérique](PRECISION_STATIQUE.md) montre que recalculer les
forces en précision étendue aux positions déjà arrondies conserve le
résidu du premier petit palier Princeton, environ 2,11e-8 N à quarante
intervalles. À rotations fixées, une reconstruction indépendante donne
4,88e-12 N ; arrondir ses positions en f64 ramène le résidu à 1,99e-8 N.

Conserver les déplacements autour de la géométrie initiale et évaluer
séparément leur contribution donne 4,73e-10 N en f64. Les extraire après
l'arrondi des positions ne récupère pas cette précision. Quatre tests
analytiques contrôlent le diagnostic, dont un allongement exact de 2^-60
perdu dans la position f64 1 + 2^-60.

Cette expérience résout seulement les translations à rotations fixées.
Le noyau est inchangé, et le mode strict échoue encore à quarante et
soixante intervalles au premier palier. Il faut porter une représentation
précise des translations dans les mises à jour, les interactions, les
dérivées et les sauvegardes, avec un contrat public d'état reproductible,
puis vérifier de nouveau le système complet et la dynamique. Aucun gain
de convergence complète ni de classement face à MBDyn n'est déclaré ici.

## Translations compensées intégrées au noyau

Le [stockage compensé](POSITIONS_COMPENSEES.md) conserve les petits
incréments dans les mises à jour statiques et dynamiques, les sous-pas
et les sauvegardes. Les poutres, superéléments et contraintes emploient
la partie basse ; la géométrie de contact est recentrée. `etat_precis()`
expose cette information et permet une restauration cinématique exacte
des données exportées, y compris des orientations.

Les cinquante paliers Princeton passent maintenant en mode strict à
10, 20, 40 et 60 intervalles pour les deux formulations : **400 paliers
à la tolérance**, aucun secours de stagnation. Un contrôle séparé des
forces aux positions conservées confirme le premier palier. La barre
analytique d'allongement 2^-60 passe aussi avec les deux poutres et un
superélément ; l'ancien binaire échoue dans les trois cas.

Les 28 tests Rust, 55 tests Python, 41 groupes de vérification,
46 bancs rapides et 9 bancs de contact passent. La vis–écrou conserve ses
contrôles d'inertie et de vitesse ; l'adjoint intégré est rejoué. Le câble
condensé extrême reste en échec, et sa validité physique reste ouverte.
Les coûts sont mesurés contre le binaire antérieur, sans en déduire un
nouveau classement chronométrique face à MBDyn.

Les vues historiques et les trajectoires restent arrondies. L'état précis
exporté n'est pas une sauvegarde de tous les historiques du modèle ; les
tolérances restent celles du contrat statique existant. Cette étape ne
prouve pas une précision universelle de toutes les interactions, ni un
avantage sur les autres fronts de l'objectif.

## Assemblage local étendu et système statique creux

Les [tangentes locales](ASSEMBLAGE_LOCAL.md) couvrent maintenant les couples,
pales, contacts, superéléments, effets gyroscopiques et réactions de liaison.
Le système statique est assemblé en triplets ; le chemin creux évite la
raideur et la matrice augmentée denses. Sur Princeton à 480 intervalles,
le contrôle alterné passe de 11,31 à 0,943 s, soit 11,99×, et de 191,4 à
66,5 Mio de pic mémoire médian. La position du bout et le travail de Newton
sont identiques, à la même tolérance stricte. Ce gain compare deux versions
de Vinkulum.

Le contrôle des dérivées corrige aussi le roulement chargé, les translations
dans les inflows spatiaux et l'amortissement faible sous forte précharge.
Les petits témoins dynamiques publient leurs coûts : environ 10 % de gain
sur le roulement, mais 9 à 12 % de surcoût pour les champs d'inflow spatiaux.
Les 31 tests Rust, 59 tests Python et les campagnes complètes passent ;
les 400 paliers stricts restent à la tolérance, le câble condensé en échec.

## Projection statique par composantes et QR creux

La [nouvelle projection des contraintes](CONTRAINTES_CREUSES.md) résout
directement les moindres carrés sur `Gᵀ`, sans `G` dense ni `GGᵀ` dans
Newton statique et sa recherche de pas. Les composantes indépendantes
sont séparées ; les facteurs symboliques creux sont réutilisés à motif
identique. Les rangs déficients conservent une SVD de norme minimale sur
la composante concernée.

Sur une chaîne articulée non linéaire à 256 corps, le contrôle alterné
mesure **1,076 → 0,01991 s**, soit **54,05×**, et **112,4 → 42,7 Mio**,
avec le même travail de Newton. Les forces, moments et fermetures sont
vérifiés indépendamment. Le petit cas à huit corps coûte 8,9 % de plus.
La boucle soudée redondante à 256 corps progresse de 3,53×, mais garde un
coût dense important. Les grands gains des modèles déjà immobiles ou
à pivots indépendants sont publiés avec cette portée limitée.

Le témoin de guides presque parallèles atteint maintenant la tolérance
jusqu'à un écart de directions de `1e-12`, contre un échec dès `1e-6`
auparavant. Les **35 tests Rust**, **60 tests Python**, campagnes complètes,
**132 équilibres alternés** et **400 paliers Princeton stricts** passent.
Ces progrès comparent deux versions de Vinkulum ; ils ne prouvent pas
la supériorité globale sur MBDyn ou Simpack.

## Réactions redondantes par QR de la transposée

Les [grands blocs larges de rang plein](REDONDANCES_CREUSES.md) utilisent
maintenant une seule résolution triangulaire et l'application implicite
de Q pour calculer les réactions de norme minimale. Le noyau évite la SVD
dense sur ce chemin, contrôle toujours les pivots, et écarte la tentative
quand le nombre de colonnes exactement distinctes prouve un défaut de rang.

Les mesures comprennent une boucle soudée droite, une boucle avec
orientations et chargement en trois dimensions, et une chaîne articulée
dont les liaisons sont dupliquées. Cette dernière conserve le coût dense
du système de Newton singulier et distingue clairement ce qu'il reste
à résoudre. Les références de réactions et les bilans des corps sont
calculés indépendamment du gradient du noyau.

Sur la boucle spatiale immobile à 128 corps, le contrôle alterné donne
**2,663 s → 2,313 ms**, soit **1 151×**, à la même tolérance stricte,
avec **53,5 → 41,0 Mio**. Les 108 équilibres comparés, 400 paliers Princeton,
38 tests Rust, 61 tests Python et les campagnes complètes passent. Le
surcoût de la petite chaîne à doublons a été corrigé par la borne de rang.
La chaîne mobile à doublons de 64 corps coûte encore environ 16,6 s ;
la précision de norme minimale à certaines grandes tailles reste publiée.

## Newton des contraintes équivalentes

La [réduction orthonormale](NEWTON_REDONDANT.md) conserve la projection des
forces, la norme minimale des réactions et le système augmenté de Newton
par une injection isométrique. La preuve s'applique aussi à une raideur
non symétrique. Les réactions complètes restent disponibles pour la
précontrainte et pour tous les contrôles physiques. Les contraintes
contradictoires et les directions presque dépendantes restent distinctes.

Sur la chaîne mobile à 64 corps et liaisons dupliquées, les mesures alternées
donnent **16,702 s → 7,159 ms**, soit **2 333×**, avec quatre évaluations
et **72,5 → 40,5 Mio**. Princeton à encastrement double converge désormais
en mode strict pour les deux formulations ; le témoin échouait dans le
même budget. Les **141 essais** conservent ces **12 échecs du témoin** et
vérifient **129 équilibres stricts**. Les **42 tests Rust**, **64 tests
Python**, campagnes générales et **400 paliers Princeton** passent.

Les parallélogrammes fermés exposent les dépendances linéaires qui ne sont
pas des lignes équivalentes et conservent leur coût dense. Le surcoût de
détection sur les autres témoins, entre environ 6 et 16 %, reste publié.

## Factorisations orthogonales et dépendances générales

Le [nouveau repli orthogonal](FACTORISATIONS_ORTHOGONALES.md) traite les
dépendances combinant plusieurs contraintes. Deux transformations de
Householder calculent la solution de norme minimale sans former les
matrices orthogonales ni les équations normales. Une séparation du rang
et un contrôle de stationnarité conditionnent son emploi ; les résidus
mécaniques complets gardent la décision de convergence.

Ce travail a aussi révélé une erreur de l'ancienne SVD sur une matrice
fabriquée de rang déficient. Le remplacement est contrôlé par reconstruction,
orthogonalité, référence algébrique et LAPACK indépendant. Une cascade
connexe de parallélogrammes fournit désormais une régression physique
dont les angles d'équilibre dérivent de l'énergie réduite.

La cascade à 8 cellules converge en 5 évaluations et celle à 16 cellules
en 6 ; le témoin échoue dans le budget fixé. À 32 cellules, les 96 corps
atteignent l'équilibre strict en 9 évaluations. Huit parallélogrammes
indépendants passent de 119,667 à 57,404 ms, soit **2,085×**. Le premier
petit calcul passe aussi sous le temps du témoin après suppression de la
préparation par blocs. Des surcoûts de **1 à 10,3 %** restent publiés sur
plusieurs témoins. Les 48 tests Rust, 65 tests Python, campagnes complètes
et 400 paliers Princeton passent.

Le calcul reste dense sur ce repli. La borne de rang concerne le triangle
calculé, sans certificat par intervalles du rang des données. Le rapport
publie les mesures alternées, les échecs conservés et la comparaison entre
la nouvelle SVD seule et son association à la factorisation orthogonale.

## Géométrie des rotations et échelles numériques

La [projection polaire](ROTATIONS_POLAIRES.md) remplace les SVD répétées
près de SO(3) par une itération locale sur matrices 3 × 3, avec résidu
recalculé et repli. Elle remplace aussi Gram–Schmidt dans le domaine des
corps libres de l'intégrateur énergie–moment. Les références fabriquées
vérifient le facteur polaire connu, les changements de repères et les
compositions réversibles. Une normalisation exacte par puissance de deux
corrige en parallèle une erreur de SVD faer sur des matrices denses entre
`10⁻³⁰⁰` et `10³⁰⁰` ; les étendues internes incompatibles avec la
normalisation f64 restent explicitement refusées.

Les temps Princeton baissent de **13,5–13,7 %** sur les deux témoins à dix
intervalles. Le corps libre conserve son ordre deux et présente une dérive
relative du moment sous `2,7 × 10⁻¹⁴` aux instantanés contrôlés. Un surcoût
dynamique de **1,9–2,3 %** reste mesuré. Les écarts entre repères à vingt
secondes persistent dans les deux versions : ils demandent un diagnostic
de phase et de sensibilité, au-delà du seul contrôle des invariants.

Une tentative de préparation des inverses d'inertie a ralenti le banc
de 23–24 % ; ses données sont conservées et son code a été retiré.
Les 52 tests Rust, 66 tests Python, campagnes complètes et 400 paliers
Princeton passent. Cette étape ne démontre toujours aucune supériorité
globale ni nouvelle comparaison directe avec MBDyn ou Simpack.

## Équilibrage spatial et raffinement creux de Newton

Le [raffinement statique](RAFFINEMENT_STATIQUE.md) corrige un équilibrage
diagonal qui amplifiait les coefficients de certaines cascades jusqu'à
700 millions. Les normes par blocs de trois coordonnées respectent la
rotation des vecteurs spatiaux, avec des échelles en puissances de deux.
Une factorisation auxiliaire creuse permet ensuite de résoudre davantage
de systèmes singuliers par résidus compensés sur la matrice originale.
La tentative est rejetée si elle ne converge pas ; les réactions finales
conservent leur projecteur de norme minimale.

La comparaison alternée de la **0.7.1 avec la roue 0.7.0** donne
**2,038 s → 37,65 ms**, soit **54,1×**, sur la cascade connexe de 96 corps,
avec **86,0 → 40,9 Mio**. Le même mécanisme tourné gagne **3,34×** ; les
32 parallélogrammes indépendants gagnent **79,6×**. Les **156 équilibres**
passent en mode strict avec références mécaniques indépendantes.

Les optimisations de préparation ramènent le pire surcoût observé de
25 % dans le prototype à **6,8 %** dans la version retenue. Les échecs des
essais LSMR sans préconditionneur adapté et la divergence d'un décalage
auxiliaire trop grand sont archivés. Les **57 tests Rust**, **67 tests
Python**, campagnes générales et **400 paliers Princeton** passent.

Le projecteur de rang déficient et le chemin public `k_c_m_z()` restent
à traiter ; ce dernier échoue sur la cascade initiale à 96 corps. Cette
étape ne comprend aucune nouvelle exécution MBDyn ou Simpack.

## Prochains développements mathématiques

Les résultats précédents orientent le travail suivant vers trois axes,
à développer et à mesurer dans le cadre de l'objectif complet :

1. **Décompositions orthogonales et rang.** Étendre le traitement des
   dépendances générales aux grandes matrices creuses, avec reconstruction
   de la solution de norme minimale et contrôle de l'erreur rétrograde.
   La factorisation orthogonale complète fournit le cadre algébrique ; le
   [guide LAPACK](https://www.netlib.org/lapack/lug/node43.html) précise le
   rôle de la transformation supplémentaire à droite. Les parallélogrammes,
   les directions presque dépendantes et les changements de rang devront
   rester dans la validation.
2. **Géométrie de l'espace admissible et analyse spectrale.** Construire une
   correction particulière des contraintes, puis résoudre dans leur noyau
   tangent avec `ZᵀKZ`, ou appliquer implicitement ce projecteur avec un
   préconditionneur approprié. Le cadre et ses hypothèses sont présentés par
   [Benzi, Golub et Liesen, §6](https://page.math.tu-berlin.de/~liesen/Publicat/BenGolLie05.pdf).
   Il faudra mesurer le remplissage des facteurs et conserver le traitement
   des raideurs non symétriques, des modes libres et des bifurcations.
3. **Géométrie différentielle et dynamique à plusieurs échelles.** Exploiter
   plus largement la structure des variétés de contraintes, les groupes de
   Lie et l'analyse des équations modifiées pour relier erreurs locales,
   dérive des invariants et coût à long terme. Les fondements de cette
   démarche sont exposés par
   [Hairer, Lubich et Wanner](https://www.unige.ch/~hairer/preprints/gniverlet.html).
   Les garanties devront préciser le régime considéré, notamment pour les
   oscillations rapides, la dissipation et les contacts. Les bancs d'énergie,
   de moment, de sous-cyclage et de contact existants seront des points de
   départ, à compléter par des confrontations externes.

## Limites suivantes

La [corrective 0.7.2](INITIALISATION.md) traite l'initialisation et l'analyse
des mécanismes redondants. La dérivée directionnelle des contraintes
supprime l'erreur de 29,3 % du cas centripète à 1 µm. Le repli dans la
métrique de masse permet l'analyse et l'initialisation de la cascade de
144 corps ; le décollement après impact retrouve son accélération libre.
Les 69 essais candidats passent, contre 30 échecs parmi les 69 témoins.
La CI étendue passe avec 66 tests Rust, 72 Python, 41 vérifications,
46 bancs mécaniques et 9 contacts. L'analyse de 144 corps consomme encore
214–215 Mio et demande 0,83–0,85 s : le coût des matrices denses reste à
traiter. Cette étape ne rejoue pas la confrontation externe.

Les grandes composantes réellement déficientes en rang, les matrices d'analyse publiques,
les autres projecteurs et les replis singuliers restent denses. Il faut
traiter ces coûts et réduire celui des dérivées aérodynamiques spatiales. La poutre intégrée
réduit le besoin de maillage sur Princeton, mais sa robustesse générale,
sa validité sur le câble et les comparaisons externes à précision commune
restent à établir. Aucune mesure Simpack n'est encore disponible.

Les autres fronts du tableau restent ouverts ; cet objectif n'est pas achevé.
