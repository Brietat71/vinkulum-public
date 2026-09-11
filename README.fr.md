# Vinkulum — documentation technique du noyau

[Accueil du projet](README.md) · [Démarrer dans FreeCAD](docs/FREECAD_FIRST_RUN.md) ·
[Contribuer](docs/CONTRIBUTOR_PROJECTS.md)

Cette référence conserve les résultats historiques avec leurs versions et limites.
**L'interface maison Vinkulum Studio est abandonnée. FreeCAD sous Linux est
l'interface de référence pour les nouveaux développements et livraisons.**
Le code et les recettes de l'ancienne GUI restent archivés pour reproduire les
résultats historiques. Ses adaptateurs de calcul utilisés par FreeCAD restent
maintenus ; aucune nouvelle fonctionnalité ni livraison de la GUI Studio n'est prévue.

Vinkulum est un moteur de simulation mécanique écrit en **Rust**, piloté par
une **API Python**. Il calcule les mouvements, les équilibres, les réactions
de liaison et les vibrations de systèmes de corps rigides et flexibles.

Sa vocation est généraliste : robotique, transmissions, suspensions, machines
industrielles, structures et mécanismes aéronautiques. FRELON, un projet
d'hélicoptère nano-UAV, est son premier cas d'application ; les modèles
aérodynamiques prolongent ce socle mécanique.

**Version courante du noyau : 0.20.0** · Python **3.14 ou plus** · Phase **alpha** ·
[Apache-2.0](LICENSE), avec [licences tierces distinctes](THIRD_PARTY_NOTICES.md).

Première publication publique des sources : voir le [dossier d'ouverture](docs/PUBLICATION_PUBLIQUE.md).
Le noyau entier n'est pas certifié ; les garanties et limites sont précisées
ci-dessous. Aucun paquet n'est publié sur PyPI.

[Versions et fonctionnalités actuelles](README.md) ·
[Certification : garanties et obligations](docs/CERTIFICATION_NOYAU.md) ·
[Référence d'API](docs/API.md) ·
[Objectif et travaux du noyau](docs/OBJECTIF_MBDYN.md) ·
[Règles de versionnement](docs/VERSIONNEMENT.md)

## Démarrer

Pour l'interface graphique, suivre le [premier lancement dans FreeCAD](docs/FREECAD_FIRST_RUN.md)
et installer son [moteur de calcul séparé](docs/FREECAD_ENGINE.md). La conception
reste dans le document FreeCAD ; les calculs s'exécutent dans des processus séparés.
Le guide distingue l'extension publiée des fonctionnalités de développement et
précise leurs limites de qualification.
Le [cahier des charges v1.1](outputs/Cahier_des_charges_suite_ingenierie_Vinkulum.md)
décrit aussi les lots futurs de la suite généraliste ; ils ne sont pas tous livrés.

Les commandes suivantes s'exécutent dans un clone du dépôt, sous un shell
Bash, avec **Python 3.14+, Rust/Cargo, un compilateur C++17 et un éditeur de liens
système, ainsi que uv** disponibles. L'installation compile l'extension en
mode optimisé et ajoute SciPy pour les références numériques indépendantes.

```bash
git clone https://github.com/Brietat71/vinkulum-public.git
cd vinkulum-public
uv venv --python 3.14 .venv
source .venv/bin/activate
uv pip install 'maturin>=1.15,<2'
maturin develop --uv --release --extras verification
python -c 'import vinkulum; print(vinkulum.__version__)'
python -m vinkulum.verification
```

Cette installation de développement utilise les modules Python du dépôt.
Après une modification du Rust, relancer la commande `maturin develop`.
Le tag Git annoté `v0.18.1` identifie la première publication ; les commits
et tags de version restent distincts des publications de paquets.

Pour construire une roue installable indépendamment du dépôt :

```bash
maturin build --release --interpreter "$VIRTUAL_ENV/bin/python"
```

La roue est écrite dans `target/wheels/` et s'installe avec `uv pip install`
suivi de son chemin. Ajouter SciPy pour lancer les vérifications. La
[livraison 0.14.1](docs/bancs/version-0.14.1.json) a été validée avec une roue
CPython 3.14 / Linux x86_64 (`linux_x86_64`) dans un environnement neuf,
hors du dépôt. La compatibilité d'une roue dépend de sa plateforme et de
son interpréteur.

## Certification du noyau

La **0.18.0** ajoute une attribution vérifiée des écarts de repères de la
toupie libre : neuf comparaisons sur 20 s, 420 000 pas comparés, moments
portés natifs et enclosures rationnelles. Les trajectoires ordinaires
mesurées restent identiques bit à bit à la 0.17.0. Cette garantie porte
sur les différences entre trajectoires discrètes ; elle ne borne pas
l’erreur à la solution continue. Voir le [diagnostic et ses preuves](docs/ATTRIBUTION_REPERES_EM.md)
et l’[audit du plan de fiabilité](docs/AUDIT_PLAN_FIABILITE.md).

La **0.14.0** étend `certifier_initialisation` avec `redondances=True`.
Elle vérifie exactement la dépendance et la compatibilité des contraintes,
établit leur rang, et certifie l'accélération unique ainsi que la réaction
généralisée unique. Les multiplicateurs sont bornés pour le représentant
choisi par le noyau ; leur répartition n'est pas déclarée unique lorsque
les contraintes sont redondantes. Voir le [contrat du quotient](docs/CERTIFICATION_QUOTIENT.md).

La **0.13.0** ajoute `certifier_systeme` et `certifier_initialisation` dans
`vinkulum.certification`. Pour les systèmes admis, un calcul exact établit
l'unicité de la solution et borne l'erreur de chaque composante. Le second
appel contrôle les accélérations et multiplicateurs réellement calculés par
le noyau, sur une copie du modèle. Les preuves JSON sont vérifiables avec
la bibliothèque standard Python seule.

Cette certification dense est facultative, limitée à 128 inconnues et à
un périmètre sans contacts non lisses ni partitions gelées pour l'initialisation.
Les contraintes écartées exigent l'option explicite `redondances=True` et
une preuve de dépendance exacte ; une proximité numérique ne suffit pas. Un refus ne démontre pas la
singularité du système. Voir le [contrat, la démonstration et l'exemple](docs/CERTIFICATION_LINEAIRE.md).

La **0.12.2** ajoute un contrôle entier exact avant l'acceptation des
vitesses assemblées. Il décide, sans arrondi, si chaque contrainte respecte
le budget annoncé sur les coefficients binary64 utilisés par le noyau.
Un cas de frontière accepté à tort par la 0.12.1 est désormais corrigé.

**Le noyau entier n'est pas certifié.** Ces garanties portent sur le résidu
d'assemblage et les systèmes linéaires admis ; la géométrie, les redondances perturbées,
les trajectoires et les contacts ont
encore des obligations de preuve. Le [dossier de certification](docs/CERTIFICATION_NOYAU.md)
donne la démonstration, les hypothèses, les contre-épreuves et les limites.
La CI confronte le garde compilé à un oracle rationnel indépendant.
Le [dossier Lean 4](preuves/README.md), ajouté après la 0.14.1, formalise
le décodage binary64 et l’équivalence du garde entier au critère rationnel.
La confrontation Lean/Rust/Fraction couvre 13 314 cas sans divergence ;
elle reste distincte d’une preuve du programme Rust compilé.

La **0.15.0** ajoute `certifier_assemblage`, en lecture seule : un certificat
local d'existence et d'unicité géométriques pour de petits mécanismes admis,
avec jauges explicites et bornes de distance à la pose native. Le
[domaine, les refus et la base de confiance](docs/CERTIFICATION_ASSEMBLAGE.md)
sont explicites ; cette garantie ne porte pas sur une trajectoire.

Un [premier dossier temporel](docs/GARANTIE_TEMPORELLE_PENDULE.md), ajouté
après la livraison 0.15.0, borne l'erreur de trois traces d'un pendule sur
une seconde. La référence utilise des intervalles rationnels, une inclusion
de Picard et un reste de Taylor explicite. Les certificats se relisent sans
le module natif ; leurs garanties concernent les nœuds enregistrés de ce
banc, avec les incertitudes initiales annoncées.

La **0.17.0** ajoute `Noyau.modes_creux`, une [analyse modale locale](docs/MODES_CREUX.md)
avec formes normalisées en masse, réactions et contrôles des équations
originales. Dix-huit cas jusqu’à 1 024 éléments sont qualifiés face à des
références indépendantes ; le rang et le spectre restent non certifiés.
La [confrontation modale](docs/CONFRONTATION_MODALE.md) conserve 140 essais
face à MBDyn et Exudyn : à 1 024 barres, Vinkulum est 1,81× plus rapide
qu’Exudyn arbre creux ; à 256 barres, cette variante Exudyn est 1,22× plus
rapide. MBDyn domine les petits cas. Ces mesures ne sont pas un classement
général des solveurs multicorps.

La **0.16.0** ajoute les [API de certification creuse](docs/CERTIFICATION_CREUSE.md),
avec enveloppes de perturbation, facteurs vérifiés sans inverse dense et
quotients à dépendances structurées. Huit certificats autonomes sont
qualifiés, dont cinq KKT natifs jusqu'à 2 816 inconnues.

## Fiabilisation du noyau existant

La **0.14.1** corrige un arrêt prématuré de Newton dans `simule_em` :
un seuil absolu pouvait figer la vitesse angulaire matérielle des corps
à petite inertie, hors axes principaux. Le résidu est maintenant relatif
au moment cinétique. Des contre-épreuves comparent huit facteurs d'inertie,
de `1e-20` à `1e20`, dans deux repères matériels, aux équations d'Euler
indépendantes. Cette correction ne constitue pas une certification temporelle.

La campagne physique conserve les erreurs d'exécution de chaque référence,
y compris Andrews, et poursuit les références suivantes. Un nouveau pilote
qualifie une roue identifiée sans écraser le rapport historique, en conservant
les tolérances et le suivi de chaque ligne. Voir les [résultats et limites](docs/VERSION_0.14.1.md)
et les [étapes prioritaires](docs/PLAN_FIABILITE.md).

La **0.12.1** corrige l'assemblage et l'initialisation dynamique : suppression
d'un plafonnement de l'accélération qui faussait même la chute libre,
refus des commandes incompatibles avec restitution de l'état, projection
robuste aux grands bras de levier et dérivées analytiques des lois.
Les contraintes non holonomes agissent au
niveau des vitesses ; les partitions gelées restent préservées.
Voir les [cas reproduits et le contrat](docs/VERSION_0.12.1.md).

## Réduction contrainte installable

La **0.12.0** ajoute `ReductionContrainte` et `AssemblageContraint` dans
`vinkulum.reduction_contrainte` : directions basses conservées, complément
spectral vérifié par intervalles, masses couplées et ports partagés.
Les certificats natifs sont embarqués dans la roue. Les bornes des réponses
restent numériques ; cette capacité concerne le modèle matériel linéaire
`K=DᵀD`, sans amortissement. Voir le
[contrat et l'exemple de résonance intérieure](docs/REDUCTION_CONTRAINTE.md).
La qualification de livraison admet six cas sur 257 fréquences et six
charges, face à des références indépendantes. Les temps historiques des
prototypes ci-dessous ne sont pas réattribués à cette nouvelle API.

## Exemple : une articulation rappelée par un ressort

Un bras tourne autour d'un pivot vertical. Un ressort de torsion de
**100 N·m/rad** s'oppose à un couple appliqué de **2 N·m**. À l'équilibre,
l'angle attendu est `couple / raideur = 0,02 rad`, soit environ **1,146°**.
Le modèle utilise les unités SI et désactive la gravité pour isoler cet effet.

```python
from math import atan2, degrees, isclose
from vinkulum import Noyau

n = Noyau(g=[0.0, 0.0, 0.0])
# Tenseur d'inertie au centre de masse, en kg·m², matrice 3 × 3 en ligne.
j = [0.02, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.01]
bras = n.corps("bras", 1.0, j, [0.0, 0.0, 0.0])
n.liaison("pivot_z", None, bras, bloque_t=[0, 1, 2], bloque_r=[0, 1])
n.couple("rappel", None, bras, [0.0, 0.0, 1.0],
         ("ressort", [100.0, 0.0, 0.0]))  # raideur, amortissement, angle au repos
n.effort(bras, [0.0, 0.0, 0.0], [0.0, 0.0, 2.0])

n.statique(tol=1e-10, strict=True)
bilan = n.statique_info()
_, rotation = n.pose(bras)
angle = atan2(rotation[3], rotation[0])
assert bilan["statut"] == "tolerance"
assert isclose(angle, 2.0 / 100.0, abs_tol=1e-10)
print(f"Angle : {angle:.6f} rad ({degrees(angle):.3f}°)")
print(f"Convergence : {bilan['statut']}")

# Conserver et restaurer l'état cinématique avec ses petites translations.
sauvegarde = n.etat_precis()
n.pose_etat(*sauvegarde)

# Retirer le couple appliqué, puis simuler les oscillations autour du repos.
n.pose_effort(0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
trajectoire = n.simule(t_end=0.1, h=0.0001)
```

Sortie attendue :

```text
Angle : 0.020000 rad (1.146°)
Convergence : tolerance
```

`None` désigne le bâti. Une liaison se ferme sur la pose courante : ses axes
bloqués sont exprimés dans son repère. Chaque entrée de `trajectoire`
contient `(t, r, R, v, w, λ)` ; l'état initial n'est pas inclus.
`etat_precis()` conserve la cinématique, mais ne sauvegarde pas tous les
historiques de contact et d'aérodynamique.

## Apport de la 0.11.0 : réponses groupées avec leurs bornes

`ReductionMaterielle.reponses(omega, forces)` traite plusieurs charges à
une même pulsation. Chaque colonne reçoit son champ physique, ses normes
et ses bornes. Les appels successifs à `reponse` réutilisent également les
calculs de la dernière fréquence ; un enrichissement invalide cet état.

Sur 128 poutres et 257 fréquences jusqu'à 20 Hz, six charges avec leurs
bornes prennent **0,358 s**, préparation comprise, contre **2,029 s** pour
l'API 0.10.0 et **1,619 s** pour le témoin HCB Exudyn le plus rapide retenu.
La LU corrigée reste devant à **0,220 s**. Les bandes à 40 Hz restent refusées.
Ces mesures concernent la réduction matérielle linéaire, à erreur physique
commune sur les fréquences testées ; elles ne classent pas les moteurs
multicorps complets. Voir le [protocole et les neuf cas](docs/REPONSES_GROUPEES_0.11.0.md)
et le [guide de l'API](docs/REDUCTION_PORTS.md).

## Apport de la 0.10.0 : réduction matérielle avec borne spectrale vérifiée

`Noyau.facteurs_materiels_poutres()` exporte directement les facteurs
d'énergie, masses et contraintes en stockage creux. La nouvelle API
`vinkulum.reduction_ports` construit un modèle réduit, reconstruit le champ
physique et enrichit sa base aux fréquences demandées. Elle conserve le
facteur statique pendant cet enrichissement.

La borne spectrale intérieure peut être calculée automatiquement : une
inverse sélectionnée fournit une trace, puis l'arithmétique à arrondis
dirigés vérifie cette trace et l'écart de factorisation. Cela évite une
résolution de valeurs propres ou l'inverse complète. Le [carnet de preuves](docs/INVERSE_SELECTIONNEE_PREUVES.md)
établit les inégalités et les conditions de refus ; le
[guide avec exemple exécutable](docs/REDUCTION_PORTS.md) décrit les unités,
forces, tolérances et diagnostics. SciPy est requis (`vinkulum[reduction]`).

Le certificat concerne la minoration spectrale du modèle matériel
`K = DᵀD` défini par les nombres d'entrée. Les autres bornes sur le champ
restent évaluées en doubles. La réduction ne couvre pas encore les
trajectoires non linéaires, les contacts ou la tangente précontrainte.
Cette livraison n'établit pas une avance généraliste sur MBDyn,
Exudyn ou Simpack.

La [confrontation harmonique avec les bases HCB officielles d'Exudyn](docs/CONFRONTATION_PORTS_EXUDYN_0.10.0.md)
mesure désormais neuf cas, à précision physique commune sur les champs,
les déformations et les ports. Elle distingue préparation, réponses
groupées et API avec bornes, et conserve les refus de bande et les
contre-calculs. Le [prototype du complément spectral](docs/COMPLEMENT_SPECTRAL_DIRIGE.md)
certifie maintenant l'espace éliminé après conservation d'une direction
intérieure : sur 512 poutres, sa limite passe d'environ 28,3 à **65,7 Hz**.
Le témoin mécanique passe les 257 fréquences de 0–40 Hz sur trois maillages,
avec une erreur maximale d'opérateur en déformation de **2,22 × 10⁻⁹**.
Le [prototype Krylov contraint](docs/KRYLOV_CONTRAINT_PROTOTYPE.md) remplace
ensuite ce KKT par des systèmes de 31, 39 et 47 inconnues. Le nouveau
[certificat par inertie dirigée](docs/INERTIE_CONTRAINTE_PROTOTYPE.md)
réduit son coût de certification d'environ **3,4×**. Dans une campagne
de 84 essais, les 24 essais avec inertie et les 24 avec trace passent
le seuil physique 10⁻⁶ ; les 24 avec contrôle passent aussi leurs majorants.
Sur 512 poutres, les champs avec inertie prennent **0,451 s contre
0,711 s pour la LU corrigée**, préparation comprise.
Le [contrôle par facteurs](docs/CONTROLE_FACTEURS_PROTOTYPE.md) réduit
ensuite le coût des majorations en conservant les champs physiques.
Dans sa campagne de 60 essais, les 24 essais des deux contrôleurs passent
la précision et leurs majorants. À 512 poutres, le nouveau service avec
contrôle prend **0,671 s contre 0,819 s auparavant et 0,696 s pour la LU**.
L'avance sur la LU est modeste et ne s'étend pas aux deux petits maillages.
À 128 poutres, il prend **0,253 s contre 1,654 s pour HCB énergie Exudyn**,
les deux au seuil demandé. Les dix refus HCB de ce nouveau lot sont
conservés sans ratio pour leurs configurations.
Le [certificat binaire avec séparateurs](docs/INERTIE_BINAIRE_SEPARATEURS.md)
réduit ensuite la préparation : **0,516 s contre 0,673 s pour le contrôle
Decimal et 0,715 s pour la LU**, au total sur 512 poutres. Les 36 essais
du nouveau lot passent, dont 24 contrôles ; champs, normes et majorants
restent identiques entre contrôleurs. Le prototype devance aussi légèrement
la LU à 128 poutres, mais reste derrière à 32. L'ordre et la preuve sont
entièrement comptés ; aucun nouveau ratio Exudyn n'est tiré de ce lot.
Le [contrôle par comparaison de masses](docs/MASSE_COMPAREE_PREUVES.md)
étend maintenant le prototype aux couplages cinétiques entre éléments,
avec la matrice de masse complète. Sur une nouvelle discrétisation à masse
consistante, les 24 calculs et les 12 contrôles passent : à 512 éléments,
le total médian est de **0,559 s contre 0,729 s pour la LU corrigée**.
À 128, les temps sont presque égaux ; à 32, le prototype reste plus lent.
Les nouvelles références sont indépendantes des anciennes masses concentrées.
L'[assemblage des compléments](docs/ASSEMBLAGE_COMPLEMENTS_PREUVES.md)
compose ensuite plusieurs pièces avec modes privés et bornes globales.
Sur trois branches de 512 éléments raccordées à une jonction rigide
linéarisée, le total médian est de **1,665 s contre 2,288 s pour la LU corrigée**.
Les 24 champs et les 12 contrôles passent. Une résonance d'une pièce isolée
peut être traversée si l'assemblage reste régulier ; une résonance globale
est refusée. Cette capacité reste expérimentale.
La [0.12.0](docs/VERSION_0.12.0.md) intègre ces capacités dans le paquet
installable ; les temps ci-dessus restent ceux des prototypes mesurés. Les [preuves de composition](docs/RETENTION_INTERIEURE_PREUVES.md)
et l'[audit de vérification creuse de 2026](docs/VERIFICATION_CREUSE_2026.md)
précisent les couplages, les singularités et les résultats retenus.

À 128 poutres sur 0–20 Hz, les champs groupés prennent **0,274 s** contre
**1,623 s** pour le meilleur réglage HCB retenu, préparation incluse.
Mais la LU corrigée prend **0,220 s**, et l'API Vinkulum avec ses bornes
**2,034 s**. Réduire le coût de préparation et de contrôle reste nécessaire.

## Apport de la 0.9.0 : gradients sans matrices de transition

Le pont adjoint temporel calcule directement les produits de dérivées :
un système transposé remplace la construction des grandes matrices `A`
et `B`. L'export `k_c_m_g_creux` assemble les tangentes sans base dense `Z`.
Des rappels optionnels peuvent aussi éviter la matrice des dérivées par
paramètre et les différences finies de l'état initial.

À **120 poutres**, le gradient d'un EI partagé prend **0,334 s contre
10,846 s** avec la roue 0.8.2, soit **32,5× plus rapide** avec l'interface
existante. Pour **120 EI distincts**, les produits locaux et une dérivée
initiale exacte donnent **0,264 s contre 11,532 s**, soit **43,6×** ; le pic
mémoire passe de **335 à 77 Mio**. Ces mesures concernent vingt pas d'un
pendule flexible, avec **156 essais finaux** et gradients contrôlés.
Un produit isolé à 30 poutres est toutefois **3,63× plus lent**, à cause
du premier chargement de SciPy inclus dans la mesure.

La [dérivation et les mesures](docs/OPERATEURS_ADJOINTS.md) précisent les
contrôles de gradient, l'équilibrage du système et les limites. Le pont
reste limité à la dynamique lisse à pas fixe, avec contraintes holonomes,
sans contact, aérodynamique, GGL ni correction σ. La trajectoire complète
reste en mémoire. L’analyse locale `modes_creux` est disponible depuis la 0.17.0.

## Apport de la 0.8.2 : efforts de poutre analytiques

Les moments nodaux des deux formulations de poutre sont calculés par
travail virtuel. Leur tangente exacte est obtenue en un passage de dérivation
automatique, ce qui supprime les différences finies du bloc de rotation et
leurs douze réévaluations d'efforts. Les sensibilités aux rigidités utilisent
les mêmes expressions. L'API et les critères de convergence sont conservés.

La [dérivation et les contrôles](docs/TANGENTE_POUTRE_ANALYTIQUE.md) couvrent
les repères, les changements d'unités, les très petits allongements et la
limite du logarithme à π. Le [point sur la composabilité scientifique](docs/COMPOSABILITE_SCIENTIFIQUE.md)
précise comment ces briques s'inscrivent dans la vocation généraliste.

## Apport de la 0.8.1 : moins de reprises en statique

Newton peut accepter un pas qui réduit la correction de pose prédite,
avec un contrôle séparé de fermeture des liaisons. Sur les cinquante
paliers Princeton mesurés, le travail passe de **879 à 334 évaluations
principales**, et de **100 à 50 tentatives**, pour les deux formulations
de poutre. Les chaînettes de 40 à 120 segments conservent leur convergence.
Les tolérances physiques et la restauration après refus restent contrôlées.

Les [mesures, équations et limites](docs/GLOBALISATION_STATIQUE.md) incluent
les témoins généralistes et les échecs : deux configurations à petite
échelle restent non résolues au réglage resserré du diagnostic. Ces
résultats internes n'établissent pas une supériorité sur MBDyn ou Simpack.

## Apport de la 0.8.0 : correction géométrique des rotations

La nouvelle option `N.simule(..., sigma_lie=0.665)` applique, à `rho=0.9`,
une cinématique issue des travaux de Holzinger, Arnold et Gerstmayr (2025).
Elle réduit nettement l'erreur d'orientation de la toupie rapide mesurée.
Le bénéfice dépend du problème et de la grandeur visée : la réaction au
pivot est moins précise au même pas, le mouvement plan ne gagne rien et
la boucle spatiale ne s'améliore pas. **Le défaut reste `sigma_lie=0`.**

L'option est expérimentale, disponible à pas fixe avec une tangente
implicite incluant GGL. L'adaptatif et le pont adjoint temporel ne la
prennent pas en charge. Les
[équations, mesures et limites](docs/INTEGRATION_SIGMA.md) rendent compte
des sept cas et des 252 essais chronométrés ; ils prolongent
l'[étude scientifique 2001–2026](docs/ETUDE_SCIENTIFIQUE_VINKULUM_2001_2026.pdf).

## Apports de la série 0.7

- **0.7.2 : accélérations initiales et analyse corrigées.** La courbure des
  contraintes est différenciée directement : l'erreur de **29,3 %** sur
  l'accélération centripète d'une bielle de 1 µm disparaît dans le cas
  analytique testé. La projection dans la métrique de masse permet
  l'analyse de mécanismes redondants auparavant refusés. Après un rebond,
  les contacts en décollement sont libérés. Voir les
  [références, mesures et limites](docs/INITIALISATION.md).
- **0.7.1 : accélération des équilibres redondants.** Un équilibrage des
  coordonnées spatiales et un raffinement creux évitent davantage de
  factorisations denses dans Newton. Les critères physiques restent
  contrôlés sur le système original. Les
  [mesures et limites](docs/RAFFINEMENT_STATIQUE.md) comprennent des
  mécanismes tournés dans l'espace et les coûts sur les autres familles.

- **Un verdict statique exploitable.** `statique(strict=True)` exige la
  tolérance demandée. `statique_info()` distingue `tolerance`, `stagnation`
  et `echec`, avec résidus, tentatives et restauration. Le mode habituel
  reste disponible ; une acceptation sur stagnation émet un `RuntimeWarning`.
- **Des petits déplacements préservés.** Les translations compensées
  évitent de perdre les incréments sous l'arrondi des positions. Les
  cinquante paliers Princeton passent en mode strict à 10, 20, 40 et
  60 intervalles, pour les deux formulations de poutre : **400 paliers**.
- **Une poutre optionnelle plus précise sur les cas mesurés.**
  `poutre(..., formulation="integree")` ajoute une formulation expérimentale
  à flexibilités condensées. Le défaut reste `formulation="milieu"`.
- **Une statique qui exploite davantage la structure des modèles.**
  Assemblage local des tangentes, résolution creuse, projection des
  contraintes par QR, traitement des liaisons équivalentes et repli
  orthogonal pour les dépendances générales.
- **Des rotations et factorisations corrigées.** Projection polaire près
  de SO(3), contrôles de rang et corrections de SVD, notamment aux échelles
  uniformes extrêmes. Les domaines de validité et les refus sont documentés.

Les [notes de la 0.7.0](docs/VERSION_0.7.0.md) détaillent la compatibilité
Python et les adaptations nécessaires aux utilisateurs directs des
structures Rust. Les [0.7.1](docs/VERSION_0.7.1.md) et
[0.7.2](docs/VERSION_0.7.2.md) conservent ces interfaces.

## Capacités disponibles

| Domaine | Capacités |
|---|---|
| Statique et assemblage | Fermeture des contraintes depuis une pose approchée, équilibre non linéaire, continuation en charge, réactions et rapport de convergence. |
| Dynamique | α-généralisé sur SO(3), option géométrique σ à pas fixe, stabilisation GGL, pas adaptatif, sous-cyclage à deux partitions ; schéma énergie–moment pour les corps libres. |
| Liaisons et actionnement | Liaison générique à degrés de liberté sélectionnables, distance, engrenage, vis–écrou, Cardan, roulement non holonome, butées, câbles, ressorts et commandes en couple. |
| Structures flexibles | Poutres géométriquement exactes anisotropes, superéléments corotationnels, plaque ACM, réduction de Guyan et Craig–Bampton. |
| Contact | Pénalité Hertz/Hunt–Crossley, barrière IPC, frottement de Coulomb régularisé, contact non lisse avec impact ; sphères et capsules face à plusieurs primitives et maillages triangulaires avec BVH. |
| Analyse | Modes réels et complexes, spectre local, bilan de stabilité, Floquet, études de convergence et audit des jacobiens. |
| Sensibilités | Dérivées de poutres et produits transposés locaux, sensibilités statiques et modales, adjoint temporel par résolution transposée dans son domaine pris en charge. |
| Entrées | API Python, lecteurs URDF et MJCF ; corpus de 21 fichiers MJCF et 9 URDF embarqué dans la roue 0.7.2, avec restrictions des lecteurs. |
| Aérodynamique | Pales à polaires c81, Wagner/Jones, Øye, Leishman–Beddoes, inflow uniforme ou Pitt–Peters, sillage libre et pale élastique couplée ; module Peters–He encore séparé du noyau. |

Les contrôles généralistes incluent notamment un hexapode de
Gough–Stewart, une caténaire et l'inverseur de Peaucellier–Lipkin
(`python -m vinkulum.domaines`), ainsi que des mécanismes plans, spatiaux,
redondants et flexibles. Les signatures, conventions de repère et restrictions
sont accessibles par `help(Noyau)` et dans la [référence d'API](docs/API.md).

## Résultats mesurés et comparaison externe

La **0.7.2 face à la roue 0.7.1** : les **69 essais candidats** d'analyse et
d'initialisation passent ; **30 essais témoins** échouent par singularité.
L'analyse de la cascade de **144 corps** termine en **0,83–0,85 s**, selon
le repère, avec **214–215 Mio** de pic mémoire. Ces valeurs décrivent un
succès nouveau ; le temps d'un échec antérieur ne constitue pas une base
de gain de vitesse. Sur les cas déjà valides, les médianes varient de
**−4,1 % à +1,5 %**. Les [138 essais alternés](docs/INITIALISATION.md)
conservent les erreurs, contrôles physiques et coûts complets mesurés.

La **0.7.1 face à la roue 0.7.0**, sur trois essais alternés par version :

| Cas statique | Temps médian 0.7.0 → 0.7.1 | Gain |
|---|---:|---:|
| Cascade connexe, 96 corps | 2,038 s → 37,65 ms | 54,1× |
| Même cascade tournée dans l'espace | 2,134 s → 639,6 ms | 3,34× |
| 32 parallélogrammes indépendants | 802,3 ms → 10,07 ms | 79,6× |

Les **156 équilibres** de cette comparaison atteignent la même tolérance
stricte et passent les contrôles mécaniques indépendants. Sur la cascade
de 96 corps, le pic mémoire passe de **86,0 à 40,9 Mio**. Des surcoûts
subsistent sur d'autres familles, jusqu'à **6,8 %** dans cette campagne.
Le [rapport du raffinement](docs/RAFFINEMENT_STATIQUE.md) publie le
protocole, les 26 configurations et les limites.

Les mesures ci-dessous proviennent **d'étapes de développement intégrées à
la 0.7.0**. Chaque rapport identifie ses binaires témoins, ses paramètres et
ses contrôles. Les bases de comparaison diffèrent entre lignes ; les gains
ne se multiplient pas entre eux.

| Cas | Résultat mesuré | Rapport |
|---|---|---|
| Princeton, 10 intervalles | Erreur au bout de 254,87 à 5,33 µm entre les deux formulations de Vinkulum, contre la référence MBDyn raffinée. | [Poutre intégrée](docs/POUTRE_INTEGREE.md) |
| Princeton, 480 intervalles | Temps statique de 11,31 à 0,943 s ; pic mémoire du processus de 191,4 à 66,5 Mio, même position finale et tolérance stricte. | [Assemblage local](docs/ASSEMBLAGE_LOCAL.md) |
| Chaîne articulée, 256 corps | Temps statique de 1,076 à 0,01991 s, à la même tolérance stricte et avec bilans mécaniques indépendants. | [Contraintes creuses](docs/CONTRAINTES_CREUSES.md) |
| Chaîne, 64 corps et liaisons dupliquées | Temps statique de 16,702 s à 7,159 ms, avec réactions individuelles de norme minimale conservées. | [Newton redondant](docs/NEWTON_REDONDANT.md) |
| Princeton, 10 intervalles et encastrement double | Temps réduit de 13,5–13,7 % par la projection polaire, selon la formulation. | [Rotations polaires](docs/ROTATIONS_POLAIRES.md) |

Les coûts défavorables sont également publiés : l'étape des rotations polaires ajoute
environ **2 %** au temps du corps libre en énergie–moment et **3,4 %** sur
le témoin de boucle spatiale à 128 corps. Ces mesures historiques restent
attachées à leurs versions ; les coûts de la 0.7.1 sont publiés dans le
[bilan du raffinement statique](docs/RAFFINEMENT_STATIQUE.md).

La **confrontation dynamique avec MBDyn utilise la roue 0.7.2** :
[protocole, trajectoires, mémoire et échecs](docs/CONFRONTATION_MBDYN_0.7.2.md).
Les configurations retenues satisfont le même seuil d'erreur estimée,
avec deux références raffinées par moteur et trois répétitions.

| Problème | Seuil | Temps total médian Vinkulum / MBDyn |
|---|---|---|
| Mécanisme plan | 100 µm | 0,269 / 0,317 s |
| Mécanisme spatial | 100 µm | 0,749 / 0,735 s |
| Andrews, initialisation corrigée | 1 mrad | 0,170 / 0,091 s |

La [confrontation statique 0.8.2](docs/CONFRONTATION_EXUDYN_0.8.2.md)
ajoute **Exudyn 1.11.0** et remesure les trois moteurs sur Princeton,
avec le même seuil de **10 µm au bout** :

| Réglage retenu | Intervalles | Temps total médian | Erreur estimée avec marge |
|---|---:|---:|---:|
| MBDyn, `beam3` | 6 | 15,88 ms | 5,740 µm |
| Vinkulum intégré, expérimental | 8 | 131,10 ms | 9,107 µm |
| Exudyn, Newton complet | 60 | 271,62 ms | 8,277 µm |
| Vinkulum milieu, défaut | 60 | 361,47 ms | 7,925 µm |

Les **232 calculs**, dont **62 échauffements**, terminent. L'intégrée
est **2,07× plus rapide qu'Exudyn** dans les réglages testés, mais reste
**8,26× plus lente que MBDyn**. Les variantes Exudyn avec Newton modifié
et binaire fast sont également mesurées ; elles ne changent pas ce
classement. Chaque variante possède deux références raffinées, et
l'effet des tolérances Exudyn est contrôlé. La formulation milieu de
Vinkulum reste plus lente qu'Exudyn sur ce cas.

Ces coûts comprennent les runtimes et les sorties ; ils ne classent pas
les seuls noyaux numériques. Cette console ne prouve aucune supériorité
généraliste, notamment sur les corps flexibles réduits, les contacts et
la dynamique. Les [mesures 0.8.1](docs/CONFRONTATION_STATIQUE_MBDYN_0.8.1.md)
et les campagnes antérieures restent archivées.
**Aucune exécution Simpack n'est archivée.**

L'objectif est de dépasser ces références en précision, robustesse, temps,
mémoire, couverture et usage. Son
[suivi](docs/OBJECTIF_MBDYN.md) distingue les progrès établis des preuves
encore attendues, notamment pour le temps réel et la validation externe.

Un [prototype de poutre mixte issu d'une formulation de mai 2026](docs/POUTRE_MIXTE_PROTOTYPE.md)
est désormais exécuté : huit éléments quadratiques donnent 5,5 µm d'erreur
sur 961 positions de Princeton. Les tests, les mesures et les limites de
conditionnement sont archivés. Il reste autonome : son coût et sa stabilité
aux rapports de raideurs extrêmes doivent être améliorés avant intégration.

L'[étude scientifique 2001–2026](docs/ETUDE_SCIENTIFIQUE_VINKULUM_2001_2026.pdf)
confronte **58 références** à cet état mesuré et aux capacités documentées
des principaux solveurs. Elle propose les priorités de résolution par
graphe, réduction flexible, sensibilités et contact, avec les expériences
nécessaires pour démontrer leurs gains. Les résultats de prototype et les
transferts à réaliser sont distingués dans le [suivi de recherche](docs/RECHERCHE_MATHEMATIQUE.md).

## Limites à connaître

- **Convergence statique.** Le mode strict contrôle le déséquilibre libre
  relatif à son échelle initiale, avec plancher 1, et les contraintes. Il ne
  fournit pas une borne absolue de l'erreur de déplacement ni une invariance
  à toutes les unités. Voir le [contrat statique](docs/STATUT_STATIQUE.md).
- **Poutres et matrices difficiles.** La poutre intégrée reste expérimentale ;
  le câble condensé extrême conserve son échec connu et sa limite physique.
  Les dépendances générales peuvent imposer un calcul dense coûteux. Une
  matrice dont la normalisation SVD perdrait des coefficients est refusée.
- **Initialisation et analyse.** Le repli par composante reste dense. L’analyse
  `modes_creux` peut revenir à une QR dense des contraintes et ne certifie
  pas le rang ; un export creux des tangentes est disponible. Les cassures de commandes ou de géométrie
  n'ont pas de dérivée seconde classique. Au redémarrage après impact,
  le frottement dépend encore de la dernière réaction connue ; il n'est
  pas résolu simultanément avec la nouvelle réaction normale.
- **Invariants et couplages.** `simule_em` est limité aux corps libres sous
  la seule gravité. La projection de moment exige une symétrie SO(3)
  globale ; celle d'énergie reste limitée aux corps libres. Le multi-rythme
  couple deux partitions par les forces, sans liaison traversante, aéro,
  contact non lisse, GGL ni adaptatif. Le correcteur de couplage reste à faire.
  Des écarts entre repères à long terme restent à diagnostiquer dans le
  [banc de rotations](docs/ROTATIONS_POLAIRES.md).
- **Géométrie de contact.** Le CCD suppose que le porteur ne tourne pas sur
  le pas. Le balayage d'une capsule utilise trois sphères et peut manquer
  un obstacle entre elles ; le contact capsule–plan traite encore le point
  `p0` comme une sphère. Une primitive déjà à l'intérieur d'une boîte n'est
  pas couverte. Les points les plus proches sont figés sur le pas et leurs
  changements de branche limitent la différentiabilité.
- **Flexibilité et adjoints.** Le noyau n'intègre pas directement de maillage
  de coques ou de solides 3D ; les superéléments importent leur raideur,
  sans raideur géométrique de précontrainte. Le pont adjoint temporel exclut
  contact, aérodynamique, liaisons non holonomes et GGL, ainsi que les
  historiques non enregistrés et changements de branche. Il garde toute
  la trajectoire. La sauvegarde en √N du modèle NumPy séparé ne réalise
  pas encore l'ordonnancement binomial de Revolve.
- **Validation aérodynamique.** Pitt–Peters est intégré explicitement hors
  du Newton mécanique ; Peters–He n'y est pas couplé. La masse ajoutée et
  l'incidence aux trois quarts de corde sont optionnelles ; CD et CM restent
  statiques dans Leishman–Beddoes. La relaxation du sillage libre échoue
  encore à faible avancement. Les désaccords Maryland et bas Reynolds
  restent publiés dans la [campagne physique](docs/VALIDATION.md).

Les priorités portent sur le passage à l'échelle des systèmes redondants,
les formulations flexibles et contacts difficiles, les couplages et dérivées,
les garanties temporelles et les comparaisons à erreur contrôlée sur des
modèles apportés par des utilisateurs indépendants.

## Vérifier et contribuer

La livraison **0.12.0** ajoute 14 contre-épreuves exécutables dans la roue,
trois tests de conservation des preuves et calculs historiques et deux
tests de l'archive de qualification. Les 129 tests Python du paquet passent
hors du dépôt ; les six cas de qualification respectent le seuil physique
10⁻⁶ et leurs majorants. Voir le [bilan](docs/bancs/version-0.12.0.json).

La livraison historique **0.11.0** a passé : **83 tests Rust**, **8 tests du prototype
flexible**, **115 tests Python publics**,
**41/41 groupes de vérification**, **46/46 bancs mécaniques**, **9/9 bancs de
contact**, formatage, Clippy et contrôle des **189 entrées d'API**. Les
contre-épreuves spécialisées comprennent **12 tests de minoration spectrale**,
avec arithmétique rationnelle exacte et contrôle des refus, **8 régressions
des réponses groupées** et **11 contrôles de leur archive**. Ces nombres
décrivent la livraison ; les campagnes spécialisées conservent leurs propres
critères et leurs échecs. Les journaux et empreintes sont dans le
[relevé de version](docs/bancs/version-0.11.0.json).

Depuis le dépôt, avec le venv de développement activé :

Le [workflow GitHub Actions](.github/workflows/ci.yml) est configuré pour les
pushes sur `main`, les pull requests et les déclenchements manuels. Il appelle
`ci/local.sh --bancs` avec Lean 4.19.0, puis un second job construit et teste
une roue CPython 3.14 hors du dépôt. Les journaux et la roue vérifiée sont
conservés comme artefacts pendant 14 jours. Consulter les
[exécutions GitHub](https://github.com/Brietat71/vinkulum-public/actions/workflows/ci.yml)
pour leur état réel ; la réussite locale ne vaut pas réussite sur GitHub.

```bash
ci/installe_hook.sh          # branche la CI locale avant chaque push
ci/local.sh                 # formatage, Clippy, tests, installation, vérification, API
ci/local.sh --bancs          # ajoute les bancs mécaniques et de contact
python -m vinkulum.domaines  # robotique parallèle, caténaire, mécanisme inverseur
python -m vinkulum.campagne  # corpus de statique et limites des stratégies
python -m vinkulum.bancs     # campagne complète et courbes travail–précision
```

`ci/local.sh` utilise le venv actif, ou `.venv` ; `PY` permet de choisir un
autre interpréteur de venv. Les campagnes parallélisables emploient jusqu'à
huit processus dans leur budget CPU. `VINKULUM_BANCS_CPUS=8` borne ce budget ;
`VINKULUM_BANCS_JOBS=1` impose une exécution séquentielle. Les mesures de
performance doivent être isolées, avec versions, nombre de fils et précision
contrôlés, conformément au protocole du rapport concerné.

`PY="$VIRTUAL_ENV/bin/python" ci/installe_neuf.sh` construit une roue et
rejoue les vérifications hors du dépôt dans un venv neuf ; ce contrôle
utilise aussi `pip` et nécessite l'accès aux dépendances.

La [référence d'API](docs/API.md) est générée depuis le paquet installé :
`python -m vinkulum.doc` la vérifie et `python -m vinkulum.doc --ecrire` la
régénère. Les [règles de versionnement](docs/VERSIONNEMENT.md) précisent la
synchronisation Rust/Python, les migrations et l'immuabilité des tags.
Les rapports techniques conservent les mesures détaillées ; le
[journal](docs/JOURNAL.md) et les [revues](docs/REVUES_EXTERNES.md) retracent
les choix d'architecture.

## Développement indépendant et licence

Le protocole clean-room du projet s'appuie sur la littérature scientifique,
les manuels, les fichiers de modèle et le comportement observable des
solveurs de référence. **Les sources des solveurs de référence ne sont pas
lues pour écrire Vinkulum.** Les modèles et résultats utilisés comme témoins
sont identifiés dans leurs rapports.

Le code et les documents originaux de cette publication sont distribués sous
[Apache-2.0](LICENSE). Les modèles tiers conservent leurs licences : leur
provenance, leurs empreintes et leurs notices figurent dans
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Pour contribuer, consulter [CONTRIBUTING.md](CONTRIBUTING.md). Les rapports
historiques conservent leurs versions, chemins de travail et hypothèses ;
ils ne constituent pas une promesse de certification générale ou de
supériorité sur tous les solveurs multicorps.
