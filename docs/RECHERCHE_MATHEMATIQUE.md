# Recherche mathématique et transfert dans Vinkulum

Point du 8 septembre 2026, après la confrontation de la roue 0.8.2 à
MBDyn et Exudyn 1.11.0, puis les opérateurs adjoints 0.9.0 et la réduction
matérielle 0.10.0 et les réponses groupées 0.11.0.

La [confrontation harmonique des ports](CONFRONTATION_PORTS_EXUDYN_0.10.0.md)
emploie maintenant les bases HCB de l'API officielle Exudyn, des
références Decimal et un contrôle de toutes les combinaisons des six
charges terminales. Les coûts de préparation, de réponse groupée et des
bornes publiques sont distingués. Cette expérience concerne des consoles
matérielles linéaires ; elle ne classe pas les solveurs multicorps complets.
Le [carnet sur la trace complémentaire](TRACE_COMPLEMENT_SPECTRAL_PREUVES.md)
dérive une voie pour conserver les directions résonantes et borner
l'opérateur restant. Le [prototype dirigé](COMPLEMENT_SPECTRAL_DIRIGE.md)
certifie maintenant ce complément, puis un témoin KKT contrôle les champs
sur 0–40 Hz : trois maillages acceptés au seuil relatif 10⁻⁶, toutes
combinaisons des six charges comprises. Les [preuves de composition](RETENTION_INTERIEURE_PREUVES.md)
conservent les couplages et les véritables résonances globales. Cette
extension reste expérimentale. Le [Krylov contraint](KRYLOV_CONTRAINT_PROTOTYPE.md)
remplace désormais ce témoin par une petite base fixe, avec
[enveloppe résiduelle uniforme](KRYLOV_CONTRAINT_PREUVES.md) et
[action du défaut par direction](KRYLOV_CONTRAINT_ANISOTROPIE.md).
Les 37 tests comprennent douze contre-épreuves rationnelles. Une réparation
des contraintes après normalisation supprime une amplification effectivement
observée ; les résidus et les majorants gardent leurs limites d'arrondi.
Dans la campagne de 60 essais, les 24 essais Krylov passent sur les trois
maillages ; les 12 avec contrôle passent aussi leurs majorants.
La [campagne suivante par inertie dirigée](INERTIE_CONTRAINTE_PROTOTYPE.md)
confronte 84 essais, dont 48 des prototypes tous acceptés physiquement
et 24 avec contrôle tous acceptés par leurs majorants. Les neuf refus
HCB sont conservés. Le certificat assemble le produit exact D.T D et
encadre chaque congruence, sans soustraction de traces. Trente tests
supplémentaires éprouvent la preuve et l'implémentation, dont 18 tests
rationnels indépendants. L'archive rejoue les trois identités distinctes
des 24 certificats d'inertie, pivots et compteurs compris.
Le [contrôle par facteurs](CONTROLE_FACTEURS_PROTOTYPE.md) ajoute ensuite
une campagne de 60 essais : les 24 essais contrôlés passent, avec champs
et normes physiques identiques entre ancien et nouveau contrôleur.
La plus grande norme d'extension est majorée depuis le Gram déjà calculé ;
les images de réparation sont factorisées avec leur défaut explicite.
Les [dix-sept contre-épreuves de preuve](CONTROLE_FACTEURS_PREUVES.md) et onze
tests d'implémentation distinguent norme d'opérateur et petite combinaison,
rang théorique et rang stocké, sous-flux et vraie nullité. Deux réfutations
du premier prototype sont conservées avec leurs corrections. La
[lecture de six sources primaires](CONTROLE_FACTEURS_SOURCES.md) situe les
estimateurs stables de 2014 et les vérifications de facteurs de 2024–2026.

La [0.11.0](REPONSES_GROUPEES_0.11.0.md) transfère la séparation entre
opérateur de réponse et second membre dans l'API publique avec bornes.
Le gain est algorithmique ; aucune nouveauté de théorème n'est revendiquée.
La contre-épreuve des petites combinaisons impose de calculer les normes
sur les champs physiques, malgré la tentation d'un Gram plus petit.
La préparation représente encore environ 78 % du coût groupé sur 512
poutres à 20 Hz. Sur le prototype contraint à 40 Hz et 512 poutres,
l'inertie réduit le certificat d'environ 0,685 à 0,199 s et permet aux
champs de devancer la LU au total : **0,451 contre 0,711 s**. Le contrôle
par facteurs abaisse ensuite les réponses contrôlées de 0,382 à **0,230 s**.
Le total devient **0,671 s contre 0,819 s pour l'ancien contrôle et
0,696 s pour la LU**, soit un gain encore modeste sur ce dernier témoin.
La préparation reprend alors environ 65 % du coût.
Les [séparateurs et congruences binary64 compilées](INERTIE_BINAIRE_SEPARATEURS.md)
abaissent ensuite le total à **0,516 s**, contre **0,673 s** pour le témoin
Decimal et **0,715 s** pour la LU du nouveau lot. Les 36 essais passent.
La permutation change les dépendances de l'élimination et permet ici une
preuve qui refusait en ordre naturel à cette précision ; le gain général
sur les largeurs d'intervalles n'est pas démontré. Les onze contre-épreuves
du noyau et sept tests d'archive couvrent les arrondis, les congruences,
les modes FTZ/DAZ et les certificats falsifiés. Les sources de 2020 et 2026
sur les permutations motivent l'axe, sans constituer une méthode nouvelle
revendiquée ni fournir les accélérations mesurées dans Vinkulum. La
[compression commune des résidus](KRYLOV_CONTRAINT_COMPRESSION.md) et la
compression des métriques par fréquence restent à éprouver avec leurs défauts
et les petites combinaisons de charges ; l'intégration publique reste ouverte.
Le calcul de trace reste sensible aux annulations et au conditionnement
des contraintes ; les refus et les bornes valides mais inutiles sont testés.

L'[audit ciblé de vérification creuse de 2026](VERIFICATION_CREUSE_2026.md)
examine les manuscrits récents de Rump et persiste une contre-épreuve exacte
d'un lemme général tel qu'imprimé. La branche symétrique positive n'est pas
réfutée par cet exemple. Une [dérivation par inertie](INERTIE_COMPLEMENT_PREUVES.md)
permet de viser directement la coercivité du complément, sans soustraction
de traces ni base dense. Dix tests rationnels éprouvent l'identité, y compris
à la frontière singulière. Le [certificat machine creux](INERTIE_DIRIGEE_PREUVES.md)
est maintenant éprouvé, avec pivots doubles et rang déduit de la régularité.
La [lecture de cinq sources primaires supplémentaires](INERTIE_DIRIGEE_SOURCES.md)
situe l'antériorité (1995), les apports de 2025–2026 et un repli par décalage
unilatéral encore non implémenté. Le principe d'inertie n'est pas revendiqué
comme un théorème nouveau.

L'[approfondissement des fondements mathématiques](FONDEMENTS_MATHEMATIQUES_VINKULUM_2026.pdf)
ajoute une lecture ciblée de 21 sources : complexité paramétrée par la
largeur du graphe (résultat de 2025), FEEC et stabilité inf-sup,
gamblets, espaces de transfert optimaux, largeurs de variétés stables,
réduction lente–rapide, fermeture GENERIC de 2025 et opérateurs de
Delassus de 2026. Les théorèmes, hypothèses, limites et expériences
réfutantes y sont distingués. Il ne constitue pas une nouvelle campagne
de performances ni une capacité livrée.

**Priorité scientifique révisée :** réduire le coût du problème à précision
physique imposée, en agissant sur la dimension nécessaire, les couplages
et la stabilité. Les trois programmes prioritaires sont l'élimination
locale des contraintes avec contrôle du rang, la réduction de
sous-structures par leurs interfaces avec contrôle d'erreur, puis une
formulation mixte et un préconditionnement robustes aux paramètres.
Les modes sélectifs en cours de prototypage restent un outil de comparaison.
La réduction adaptative avec mémoire et contact est une extension à démontrer.

Le [prototype de quotient orthogonal local](CONTRAINTES_ORTHOGONALES_PROTOTYPE.md)
éprouve maintenant le premier axe : projection et efforts sans base dense,
avec contre-calcul cinématique des cascades redondantes. La lecture
complémentaire de Foster–Davis (2013), Scott–Tůma (2022) et spaQR (2020)
précise l'antériorité et le rôle de la compression des interfaces.
Le remplissage des facteurs, la sensibilité au repère et un contre-exemple
de mobilité au seuil de rang empêchent une généralisation au noyau à ce stade.
Un [carnet de preuves](QUOTIENT_LOCAL_PREUVES.md) dérive la composition
par quotients, la masse induite, le report exact des pivots faibles, une
borne d'erreur physique et une limite de dimension dynamique. Il montre
aussi pourquoi une dépendance du jacobien à une pose singulière peut
cacher une compatibilité non linéaire d'ordre deux.

La comparaison doit inclure les arbres cinématiques et FFRF/HCB d'Exudyn,
la CMS de MBDyn et les réductions non linéaires déjà proposées par Simpack.
La combinaison envisagée ne bénéficie d'aucune exclusivité mondiale établie.

Le [prototype de réduction par ports](PORTS_KRYLOV_PROTOTYPE.md) éprouve
maintenant le deuxième axe avec une référence HCB creuse et 56 essais
archivés. Son [carnet de preuves](PORTS_KRYLOV_PREUVES.md) établit une
enveloppe résiduelle uniforme, conserve les défauts de rang et explique
la décroissance algébrique de l'erreur modale sur la chaîne. Les gains
locaux atteignent 4,90× et 17,76× sur les grandes chaînes testées,
construction et contrôle inclus ; le petit cas régresse et la tolérance
1e−10 échoue en flottants. Les [sources complémentaires](PORTS_KRYLOV_SOURCES.md)
précisent le rattachement aux encadrements matriciels de 2025. La réduction
n'est pas encore intégrée aux trajectoires du noyau. La confrontation
harmonique 0.10.0 teste désormais sa déclinaison matérielle native face
aux bases HCB officielles d'Exudyn.

Le [relèvement en énergie](PORTS_RELEVES_PROTOTYPE.md) traite ensuite
la perte de précision du Schur : factorisation QR du facteur élémentaire,
ou LU corrigée par des résidus énergétiques, avec audit uniforme dans
le modèle d'entrée. Les oracles à 70 chiffres et 60 essais confirment
le seuil 1e−10 sur les consoles testées jusqu'à 512 éléments. Cette
précision coûte plus cher que la voie directe lorsque celle-ci suffit.
L'assemblage de sous-structures et son contrôle de réponse sont éprouvés ;
le passage aux tangentes non linéaires reste une exigence à traiter.
Le contrôle du champ a depuis été ajouté, comme décrit ci-dessous.
Le [carnet de preuves](PORTS_RELEVEMENT_PREUVES.md)
distingue erreur de réduction, perturbation du modèle et amplification
globale ; les [six sources supplémentaires](PORTS_RELEVEMENT_SOURCES.md)
précisent les antériorités et les limites de précision flottante.

Le [contrôle du champ intérieur](CHAMP_INTERIEUR_PROTOTYPE.md) traite
maintenant les déplacements en norme massique, les déformations en norme
énergétique et les observables linéaires. Le [nouveau carnet de preuves](CHAMP_INTERIEUR_PREUVES.md)
établit le passage du résidu au champ, son amplification après assemblage
et une contraction locale sous coercivité. L'enrichissement conserve les
facteurs et la base existants ; trois expériences satisfont le seuil relatif
1e−10 aux fréquences demandées, deux le refusent près d'une résonance.
Cette acceptation ne couvre pas toute la bande globale. L'[extraction native
du facteur matériel des poutres](FACTEURS_ENERGIE_NOYAU.md) fournit aussi
les déformations pondérées et leur jacobien, sans former K. La 0.10.0
raccorde cette extraction aux [ports Python publics](REDUCTION_PORTS.md).
Elle ne remplace pas la tangente précontrainte et ne fournit pas encore
de trajectoires réduites.

Le [carnet de minoration spectrale](INVERSE_SELECTIONNEE_PREUVES.md)
démontre maintenant une borne obtenue par trace d'inverse sélectionnée,
puis contrôle l'écart énergétique de la factorisation. Les deux calculs
sont encadrés par arithmétique à arrondis dirigés : la minoration vise
les nombres binary64 de D et M fournis, sans supposer le QR exact.
Le procédé évite une inverse complète et un calcul modal, sous des
budgets explicites de remplissage et d'opérations. Le minorant peut être
conservateur ; le coût n'est pas linéaire pour tout graphe et les autres
briques de contrôle du champ restent en flottants non certifiés.

L'[étude scientifique 2001–2026](ETUDE_SCIENTIFIQUE_VINKULUM_2001_2026.pdf)
réunit désormais 58 références de mathématiques, informatique, physique
et documentation des solveurs. Elle précise les hypothèses des méthodes,
leurs limites, leur adoption chez les concurrents et les expériences
requises pour décider des transferts dans le moteur généraliste.
La [lecture critique du document sur la composabilité](COMPOSABILITE_SCIENTIFIQUE.md)
précise maintenant les exclusions à corriger, les hypothèses des gradients
de contact et les conditions à partager entre les briques du noyau.

**Le recours aux avancées récentes n'est pas encore à la hauteur de
l'objectif.** Les corrections récentes de Vinkulum emploient principalement
des outils éprouvés : différentiation automatique, géométrie de SO(3),
projections orthogonales et raffinement. Leur utilité est vérifiée, mais
elle ne justifie pas de qualifier tout le noyau d'état de l'art.

La [confrontation mesurée](CONFRONTATION_EXUDYN_0.8.2.md) fournit les
témoins à dépasser. Sur Princeton au seuil de 10 µm, l'option intégrée
est 2,07× plus rapide que le meilleur réglage Exudyn testé, mais reste
8,26× plus lente que MBDyn en temps total. Le défaut milieu reste plus
lent qu'Exudyn. La globalisation et les tangentes analytiques améliorent
donc le noyau sans établir une avance généraliste.

Exudyn dispose déjà d'une cinématique σ issue de 2025, de modes creux
et d'une chaîne FFRF/Hurty–Craig–Bampton. Le rapport donne les sources
primaires et les domaines comparés. Ces méthodes sont des références à
composer et à dépasser par des mesures, pas des exclusivités de Vinkulum.

## Priorités issues de l'étude scientifique

L'hypothèse stratégique est de réduire le coût des grands mécanismes
connectés rigides et flexibles, avec des sensibilités physiquement
vérifiées. Elle reste à démontrer sur plusieurs familles de problèmes.

- **Intégration géométrique.** La modification sigma de 2025 a été
  [comparée et livrée en option dans la 0.8.0](INTEGRATION_SIGMA.md).
  Son gain reste ciblé, avec des réactions parfois dégradées. Continuer
  les comparaisons à erreur commune sur plusieurs observables.
  [Holzinger, Arnold et Gerstmayr](https://doi.org/10.1016/j.mechmachtheory.2025.106236).
- **Structure et sensibilités.** Exploiter le graphe, les éliminations
  locales et les produits directs/transposés, puis différencier les
  équations convergées. Rang, métrique et résidu adjoint doivent être
  contrôlés. [Proximal creux, 2021](https://www.roboticsproceedings.org/rss17/p017.pdf)
  et [différentiation implicite, 2022](https://papers.neurips.cc/paper_files/paper/2022/file/228b9279ecf9bbafe582406850c57115-Paper-Conference.pdf).
- **Flexibilité réduite.** Préserver la structure mécanique, compter le
  coût hors ligne et vérifier les charges hors construction, avec un
  enrichissement ou un repli à dériver. [Carlberg, Tuminaro et Boggs](https://arxiv.org/pdf/1401.8044).
- **Contact et événements.** Résoudre normale et tangente ensemble,
  contrôler la trajectoire courbe réellement suivie et dériver l'instant
  d'événement quand les hypothèses le permettent. [Comparaison des lois](https://simple-robotics.github.io/publications/contact-models/static/paper/lelidec2024contacts.pdf)
  et [matrices de saltation](https://arxiv.org/html/2306.06862v3).

Les ports énergétiques, le multi-rythme, la précision mixte, le GPU et les
préconditionneurs appris complètent ces axes selon les coûts mesurés.
Le rapport contient le programme commun de comparaison et les motifs de
refus d'une intégration. Les transferts livrés sont identifiés explicitement ;
les autres propositions ne sont pas des fonctionnalités disponibles.

## Différentiation implicite par produits — livraison 0.9.0

Le [pont temporel](OPERATEURS_ADJOINTS.md) compose maintenant les tangentes
locales du noyau, un système contraint creux et les produits transposés de
sensibilité. Il n'exporte plus une base dense inutile et ne construit plus
les matrices de transition `A` et `B`. C'est une application mécanique du
principe modulaire de différentiation implicite de Blondel et al. (2022),
fondé sur un théorème classique, avec les transports de covecteurs sur SO(3).

Sur la chaîne de 120 poutres et vingt pas, le gradient d'une rigidité
partagée gagne **32,47×** contre la roue 0.8.2 avec l'interface existante.
Le gradient de 120 rigidités gagne **43,61×** avec produits locaux et
dérivée initiale exacte pour ce modèle. Les différences finies indépendantes
et le contre-calcul dense contrôlent les gradients ; le résidu adjoint est
équilibré et vérifié. Le coût du premier import de SciPy dégrade cependant
le produit isolé à 30 poutres, et le stockage de trajectoire reste complet.

Cette brique rend l'optimisation des grands modèles moins coûteuse dans
son domaine holonome lisse. Elle ne traite pas les gradients d'événements
ni les historiques cachés et ne fournit pas une borne de conditionnement.
Les modes creux sélectifs et les sensibilités statiques restent un travail
distinct. Aucune comparaison d'adjoints Exudyn, MBDyn ou Simpack n'est
déduite de ces gains internes.

## Efforts et tangente des poutres — livraison 0.8.2

La [dérivation par travail virtuel](TANGENTE_POUTRE_ANALYTIQUE.md) donne
les moments nodaux analytiques des deux énergies existantes. Leur dérivée
par duaux d'ordre un remplace les différences finies du bloc de rotation.
Le contre-calcul par énergie et dérivation emboîtée reste indépendant.
Les tests couvrent notamment l'objectivité, les changements d'unités,
les rigidités signées des sensibilités et la coupure du logarithme à π.
La comparaison entre roues figées donne un gain de 1,98 à 2,56 fois
sur les six rampes Princeton ; les 49 contrôles précédemment réussis
sont conservés. Les deux échecs à petite échelle restent présents.

## Poutre mixte avec moments indépendants — priorité de formulation

Le préprint de **Humer, Steinbrecher et Pechstein, mai 2026**, introduit des
moments et des rotations internes indépendants, avec courbure discrète
aux interfaces. Sa forme hybride conserve positions et rotations nodales
comme inconnues globales et permet une condensation des champs internes.
Les sections théoriques 2–3 et les résultats numériques ont été consultés.
La construction part d'une transformation de Legendre de l'énergie de
flexion ; au premier ordre, la rotation interne est constante par élément.
[Article et équations](https://arxiv.org/html/2605.04573v1).

**Transfert réalisé en prototype autonome, pas encore intégré au noyau :**
les ordres 1 et 2 sont implémentés en Rust, avec élimination analytique des
moments, résolution interne et tangente condensée différenciée. Huit tests
passent, dont l'objectivité, la flexibilité anisotrope, la tangente et une
jonction ramifiée. Une référence continue indépendante et 36 exécutions
Princeton mesurent la précision sur toute la ligne.
[Dérivation, résultats et sources exécutées](POUTRE_MIXTE_PROTOTYPE.md).

L'ordre 2 donne 2,9 µm au bout avec trois éléments, mais il faut huit
éléments pour atteindre 5,5 µm sur les 961 points de la ligne. Son coût
reste défavorable et une raideur négative artificielle apparaît dans la
condensation quadratique au rapport de raideurs 10¹⁶. Les sondages conservent
ce défaut et le refus à flexion nulle. La stabilité numérique, les branches
de rotation, la dynamique et les couplages doivent être traités avant
l'intégration. Une méthode récente n'est pas automatiquement une méthode
plus fiable dans le noyau.

## Approximation géométrique d'ordre supérieur — alternative à mesurer

Une publication de **2025** étudie les B-splines généralisées sur SE(3)
pour des poutres géométriquement exactes. Son résumé décrit une
interpolation plus régulière et des gains de convergence sur ses exemples.
Le résumé éditeur a été consulté ; la dérivation complète reste à examiner
avant toute implémentation.
[Publication](https://www.sciencedirect.com/science/article/pii/S0045782525002518).

**Transfert envisagé :** comparer le raffinement en ordre au raffinement
en nombre d'éléments, à erreur commune, pour des chaînes et des structures
ramifiées. Il faut compter les inconnues, les quadratures, les dérivées et
les sorties. Un gain de convergence dans un article ne donne pas directement
un gain de temps dans le noyau.

## Résolution non linéaire — associer préconditionnement et globalisation

Un travail de **2024** sur Newton–CG pour l'optimisation de maillages
courbes associe stratégie de globalisation, préconditionnement et précision
variable des résolutions. Son résumé a été consulté : c'est une piste de
méthode, pas une preuve sur les équilibres multicorps.
[Prépublication](https://arxiv.org/abs/2403.13654).

**Transfert à dériver :** agir sur les variables ou les équations qui
déséquilibrent le problème non linéaire, et adapter le travail linéaire au
progrès effectivement obtenu. Un CG destiné à une optimisation ne se
branche pas directement sur les matrices KKT de Vinkulum : leur
indéfinition, les contraintes redondantes, les forces non conservatives
et les contacts imposent d'autres hypothèses. Le système physique original
doit rester celui du contrôle final.

## Preuves requises pour intégrer une méthode

Chaque piste doit aboutir à une chaîne vérifiable : **référence précise,
hypothèses, dérivation propre, prototype exécuté, comparaison au témoin
figé, intégration et contrôles généralistes**. Une référence récente ou
un résultat de prototype ne vaut pas une capacité livrée.

Pour la poutre mixte, les premiers contrôles exécutés couvrent la limite
linéaire anisotrope, l'arc sous moment pur, les mouvements rigides
superposés, les dérivées, une jonction et Princeton. La limite du câble
tendu reste ouverte. La dynamique exigera ses propres bilans d'énergie
et de quantité de mouvement.
Pour Newton, il faudra conserver les restaurations après refus et rejouer
les mécanismes redondants, les échelles extrêmes et les cas où le mérite
en forces fonctionne bien. Les mesures devront publier les coûts
défavorables autant que les gains.

Ces axes complètent les autres fronts de l'[objectif généraliste](OBJECTIF_MBDYN.md).
Les méthodes encore au stade de piste ou de prototype ne constituent pas
des capacités livrées. Aucun de ces résultats ne démontre une supériorité
générale sur MBDyn ou Simpack.
