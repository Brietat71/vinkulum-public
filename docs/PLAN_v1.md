# Vinkulum v1.0 — dépasser Abaqus, Nastran, Simpack, HOST et CAMRAD II

Ambition posée par Paul le 3 septembre 2026. Ce document dit ce qu'elle
signifie **exactement**, ce qui est atteignable, ce qui ne l'est pas, et
comment on le PROUVE. Rien ici n'est une prétention : chaque ligne est une
mesure à produire sur un banc que d'autres peuvent rejouer.

## 1. Ce qu'on ne dépassera pas, et il faut le dire d'abord

| code | ce qu'il a et qu'on n'aura pas | ordre de grandeur |
|---|---|---|
| **Nastran** | 60 ans de bibliothèque d'éléments, de solutions et de corrélation | ~10⁶ lignes, des milliers d'éléments |
| **Abaqus** | matériaux, contact, rupture, mise en forme, thermique, explicite | ~10⁷ lignes, 40 ans |
| **Simpack** | bibliothèques métier (rail, automobile), interfaces CAO, co-simulation | 30 ans, écosystème |
| **CAMRAD II / HOST** | corrélation essais en vol sur des dizaines d'appareils réels | 35 ans, données propriétaires |
| **RecurDyn** | **contact** (BVH sur Parasolid, pénalité régularisée), MFBD sur solides et coques, bibliothèques métier | 25 ans, cœur de métier |

**La largeur fonctionnelle n'est pas rattrapable, et la viser serait la
garantie de tout rater.** Ce qu'on peut faire, c'est être **strictement
meilleur sur des axes où ces codes sont structurellement faibles** — faibles
par leur âge et leur architecture, pas par manque de talent.

## 2. Les cinq axes de dépassement, et pourquoi ils sont tenables

### A. Différentiabilité de bout en bout — l'axe décisif

Aucun de ces cinq ne rend le **gradient exact** de n'importe quelle sortie par
rapport à n'importe quelle entrée. Nastran a des sensibilités sur un
sous-ensemble (SOL 200), Abaqus des dérivées de conception limitées, Simpack
et CAMRAD II essentiellement rien. Tous obligent aux différences finies :
n+1 simulations pour n paramètres, avec le bruit d'arrondi qu'on a mesuré
ici même (une composante nulle sort à 5e-3 de bruit à ε = 1e-7).

Vinkulum a déjà l'AD dans le jacobien. L'étendre au MODE INVERSE sur la
trajectoire entière donne ∂(n'importe quoi)/∂(tout) **en un seul passage** :
optimisation de conception, identification de paramètres, assimilation de
données d'essai, contrôle optimal. C'est un saut de nature, pas de degré.

- **Métrique** : gradient de l'endurance FRELON par rapport aux 20 cotes de la
  tête, exact à 1e-10, en < 5× le coût d'une simulation. Les DF demandent 21×
  et rendent 3 chiffres.
- **Banc** : identification des raideurs d'une poutre à partir de sa réponse.

### B. Couplage monolithique mécanisme + flexible + aéro + commande

Chacun de ces codes couvre une part et **co-simule** le reste : CAMRAD II fait
le rotor mais pas la structure détaillée de la cellule ; Abaqus fait la
structure mais pas le rotor en vol ; Simpack appelle un aéro externe. La
co-simulation coûte de la stabilité, de la précision d'interface et du temps.

Vinkulum résout **un seul système DAE** : liaisons, engrenage, poutres,
pales, inflow, servos et gouverneur y sont des éléments du même résidu, avec
le même jacobien. C'est déjà vrai aujourd'hui — la tête S2 en vol le prouve.

- **Métrique** : un cas où le couplage change le résultat de plus de 5 % par
  rapport à la co-simulation (le droop servo sous moment aéro en est un
  candidat déjà observé au n1).

### C. Vérifiabilité — une supériorité méthodologique, pas technique

Ces cinq codes sont des boîtes noires. Leur validation est interne et, pour
l'essentiel, non rejouable par l'utilisateur.

Vinkulum publie sa **suite de bancs** (`docs/bancs/<version>.json`), chaque
brique contre une référence indépendante, l'historique des versions étant la
preuve. Deux diagnostics (`fuite_transport`, `diag_modal`) existent pour
qu'une explication puisse être **réfutée**.

- **Métrique** : tout chiffre publié est rejouable par un tiers en une
  commande, et chaque banc nomme la référence qui pourrait le contredire.

### D. Déterminisme au bit

Le parallélisme de ces codes n'est pas reproductible au bit (réductions en
ordre variable). Vinkulum l'est, et c'est mesuré (1 fil / tous les fils :
identiques). Pour la certification et le débogage, c'est décisif.

### E. Coût par pas à précision donnée

Pas « plus rapide » : **la courbe travail–précision**. La tête S2 (143
inconnues, boucles fermées, engrenage, aéro) tourne à 0,23 ms/pas ; six
secondes de vol en 16 s de calcul.

- **Métrique** : sur un cas modélisable par les deux, produire la courbe
  erreur/temps de vinkulum et celle du concurrent. Sans cette courbe, aucune
  affirmation de supériorité n'a de sens.

## 3. Les bancs de la barre — publics, rejouables, non choisis par nous

Prétendre dépasser CAMRAD II exige de se comparer sur SES cas, pas sur les
nôtres. Cibles, par ordre de dureté :

| banc | ce qu'il juge | référence |
|---|---|---|
| **Bathe–Bolourchi** cantilever 45° | poutre GE en grands déplacements 3D | solution publiée |
| **Princeton beam** (Dowell–Traybar) | flexion-torsion couplée, mesures | essais publiés |
| **IFToMM** (Andrews, Bricard, régulateur) | multicorps raide et redondant | valeurs de référence |
| **Southwell** poutre en rotation | raidissement centrifuge, ν_β | analytique |
| **HART-II** | rotor en soufflerie, charges et BVI | données publiques |
| **UH-60A Airloads** | corrélation essais en vol — la barre de CAMRAD II/HOST | données publiques NASA |
| **NREL 5 MW / IEA 15 MW** | pale flexible en rotation, aéroélasticité | OpenFAST, public |

**Tant que HART-II et UH-60A ne sont pas passés, « faire mieux que CAMRAD II »
n'est pas une affirmation défendable** — c'est une intention. Le dire
autrement serait exactement la faute de méthode corrigée le 3 septembre.

## 4. Ce qui manque pour y aller, dans l'ordre

1. **AD en mode inverse sur la trajectoire** (axe A) — le différenciateur.
2. **Pas adaptatif et contrôle d'erreur** : aucun code sérieux n'impose un pas
   fixe à l'utilisateur.
3. **Poutre en rotation validée** (Southwell) puis **pale flexible** — sans
   quoi rien de rotorcraft n'est crédible.
4. **La limite gyroscopique** : formulation en repère tournant ou Floquet.
   C'est le verrou identifié et mesuré (‖Kr‖ = 50,85 contre mgl = 0,049 :
   soustraction catastrophique).
5. **Inflow dynamique** (Pitt–Peters, puis sillage prescrit) — CAMRAD II vit
   là ; sans modèle de sillage, pas de HART-II.
6. **Craig–Bampton** pour importer une cellule maillée ailleurs : c'est ainsi
   qu'on cohabite avec Nastran plutôt que de le refaire.
7. Contact et frottement — pour Simpack et Abaqus, pas pour le rotor.

## 5. La règle qui vaut plus que le plan

Aucune de ces sept briques n'entre sans sa référence indépendante ET sans le
diagnostic capable de la réfuter. Le 3 septembre, une explication publiée sans
mesure s'est révélée fausse en une exécution ; c'est le genre d'erreur qui,
répétée, transformerait ce plan en prospectus. **On mesure, ou on se tait.**

## 6. Ce que le rapport SABRE D1.6 apporte (lu le 3 septembre)

Livrable D1.6 du projet européen SABRE (TUM, mai 2018) : couplage faible
**CAMRAD II ↔ TAU** (CFD DLR) sur un rotor **HART-II**. Il documente de
l'intérieur la façon dont l'état de l'art rotorcraft est réellement construit.

### 6.1 La citation qui valide l'axe B, écrite par la partie adverse

> « A third approach has also been suggested in which a single code is
> developed to solve all of the differential equations governing aerodynamics
> and structural dynamics. This is known as the **monolithic approach**.
> However, for rotorcraft applications, **such a monolithic approach has not
> been developed and successfully tested**. […] it offers relatively few
> advantages compared to the amount of work required to implement it, namely a
> huge increase in development time. » (§1.1, p. 3)

Le monolithique n'est pas écarté parce qu'il serait inférieur : il est écarté
parce qu'il **coûte trop cher à développer** quand on part de codes existants.
Vinkulum le fait par construction et l'a déjà démontré sur la tête S2
(liaisons, engrenage, poutres, pales, inflow, servos et gouverneur dans un
seul résidu). C'est l'axe où la position est structurellement favorable — et
la seule chose qui manque est l'échelle, pas la méthode.

### 6.2 Ce que fait vraiment CAMRAD II, et ce que vinkulum a déjà

| CAMRAD II (§2.1) | vinkulum aujourd'hui |
|---|---|
| aéro **lifting line 2D + tables** (Cl, Cd, Cm vs α et Mach) | même chose : élément `Pale`, tables c81 |
| trim par **jacobien des sorties vs commandes + Newton–Raphson** | absent — **c'est la brique manquante n° 1** |
| capteurs structuraux/aéro en sortie | `aero()`, `couples()`, `poutres()`, `modes()` |
| pales flexibles par sections/nœuds | poutres géométriquement exactes (validées) |

Le niveau aérodynamique de CAMRAD II **seul** est celui que vinkulum a déjà.
Ce qui manque pour se comparer à lui n'est pas la physique de section : c'est
le **trim**. Et son jacobien de trim est obtenu par différences finies — là où
notre AD le rendrait exact (axe A).

### 6.3 Le coût réel du couplage, et ce qu'il ouvre

45 heures sur **280 cœurs** pour UNE révolution de rotor en CFD ; 4 itérations
de couplage ; la limite de 48 h par job du supercalculateur ne suffit même pas
à une solution couplée, d'où un enchaînement manuel. Le couplage faible
n'échange qu'**une fois par tour** — c'est ce qui le rend praticable, et c'est
aussi sa limite de fidélité.

### 6.4 Le banc HART-II, entièrement paramétré par ce rapport

| paramètre | valeur |
|---|---|
| μ | 0,15 |
| ρ · T∞ | 1,2055 kg/m³ · 17,3 °C |
| profil · corde | NACA 23012 · 0,121 m |
| vrillage | −8 °/R, nul à r = 0,75 |
| pales · R · σ | 4 · 2,0 m · 0,077 |
| ω | 109 rad/s |
| θ_shaft | −4,5° |
| cibles de trim | T = 3300 N, Mx = My = −20 N·m |

**Résultats de trim de CAMRAD II à reproduire** (Table 4.3, itération 1, non
couplée — le cas comparable à vinkulum sans CFD) :

```
   Θcoll 6,45°   Θlat −2,03°   Θlong 2,69°      T_CAMRAD 3298,4 N
```

et la trajectoire de convergence du couplage (Θcoll 6,45 → 5,96 → 5,95 →
5,93 ; T_TAU 3805,8 → 3297,2 → 3313,2), utile comme cible d'un futur couplage.

**C'est le premier banc rotorcraft directement opposable**, et il ne demande
pas de CFD : reproduire le trim de CAMRAD II en lifting line sur le même
rotor, aux mêmes cibles. Réserve du rapport à ne pas oublier : **ses pales y
sont supposées RIGIDES** et un seul capteur structural est utilisé — c'est
précisément là que nos poutres géométriquement exactes peuvent faire mieux, et
c'est mesurable.

### 6.5 Ce que ça change à l'ordre de marche

Le **trim** monte en tête, devant l'AD inverse : sans lui, aucune comparaison
rotorcraft n'est possible, et il est le premier consommateur naturel de l'AD
(jacobien de trim exact au lieu de différences finies). Ordre révisé :

1. **trim** (Newton sur les commandes, jacobien par AD) → banc HART-II §6.4 ;
2. AD inverse sur la trajectoire ;
3. poutre en rotation (Southwell) puis pale flexible → l'axe où le rapport
   déclare lui-même sa faiblesse ;
4. le reste inchangé.

## 7. Ce que le papier HOST apporte (lu le 3 septembre)

*« HOST, a General Helicopter Simulation Tool for Germany and France »*,
Benoit, Dequin, Kampa, von Grünhagen, Basset, Gimonet — AHS 56ᵗʰ Annual Forum,
2000. C'est le code commun Eurocopter / ONERA / DLR, et le papier décrit son
ARCHITECTURE, ce qui en fait la source la plus directement utile des cinq :
SABRE disait ce que CAMRAD II calcule, HOST dit comment un noyau
d'aéromécanique est organisé et pourquoi.

### 7.1 Le noyau a TROIS fonctions, et nous n'en avions qu'une et demie

> « The functions of HOST are based on three main utilities: the trim
> calculation, the time domain simulation and the calculation of linear
> equivalent system. The others are pre- and post-processing routines and are
> a combination of that three ones. »

| fonction du noyau | HOST | vinkulum au 3 septembre |
|---|---|---|
| domaine temporel | oui | **oui** — α-généralisé sur groupe de Lie, index 3 |
| trim | oui, harmonique, loi définie par l'utilisateur | **fait ce jour** (`vinkulum.trim`) |
| système linéaire équivalent | oui, périodique → coefficients constants | modes propres et complexes seulement |

La troisième est celle qui manque encore, et le papier en donne la recette
exacte : l'état est remplacé par ses composantes harmoniques (X₀, X_ic, X_is),
les dérivées temporelles font apparaître les termes ±2iΩ et −(iΩ)², et on
obtient un système à **coefficients constants** sur lequel les modes
collectif, progressif et régressif se lisent directement. Deux traitements
précèdent : correction de couplage inertiel (les perturbations sont appliquées
à accélération moyenne nulle) et réduction quasi-statique.

> **C'est la réponse à la limite gyroscopique qu'on a mesurée le 3 septembre
> au matin.** `modes_complexes` refuse dès qu'un corps contraint tourne, parce
> que la raideur tangente en repère fixe est dominée par des termes en
> ω²(J_t − J_a) valant mille fois la raideur cherchée. HOST ne linéarise pas
> en repère fixe autour d'un état tournant : il passe aux composantes
> périodiques de l'état. « In HOST, the stability is done using the periodic
> components of the state variables. » Le refus reste juste ; la voie de sortie
> est nommée par un code qui la pratique depuis trente ans.

### 7.2 Le trim, et ce qu'on lui a repris

Deux choses, écrites dans `python/vinkulum/trim.py` :

- **la matrice d'influence se garde** — « The convergence scheme iterates with
  the same matrix to reach the solution. As the matrix calculation is time
  consuming, it is only done again when the matrix becomes inaccurate. » C'est
  l'invariant du Newton simplifié de notre intégrateur, un cran plus haut ;
- **la loi de trim appartient à l'appelant** — « A trim law defines the way in
  which the trim has to be found. It gives the parameters to be imposed […]
  and the parameters to set free. » Notre trim ne sait donc rien d'un rotor.

Ce qu'on a **mesuré et qui nuance le premier point** : garder la matrice n'est
pas gratuit, c'est un échange d'itérations quadratiques contre des itérations
linéaires. Sur nos cas algébriques, elle gagne d'un facteur 1,2 à n = 3 et
2,0 à n = 12 en non-linéarité modérée — et elle **perd** (×1,27) quand la
non-linéarité est forte, parce que la matrice initiale décrit une autre
machine. Le seuil de réfection existe pour ça, et l'autotest garde les deux
faits. HOST a raison pour son cas, où le trim porte sur « several hundred
influencing parameters » réduits par isotropie : à grand n la matrice domine
tout.

Ce que HOST a et que nous n'avons pas : la représentation harmonique de l'état
gérée par le noyau (trim sur X₀/X_ic/X_is au lieu d'une moyenne mesurée sur un
tour), et la réduction par isotropie. Déclaré dans le module.

### 7.3 Trois fonctions qui découlent du trim, et qu'on aura presque gratuitement

- **simulation inverse** : objectifs sur les SORTIES, commandes ajustées par
  matrice d'influence sur un pas ΔT_inv > Δt, avec mémorisation et
  restauration de l'état. HOST l'a développée et finalisée **en trois
  semaines** grâce à la structure — c'est la mesure de ce que l'architecture
  achète. Notre `Modele` est déjà clonable, donc snapshot/restore existe ;
- **identification paramétrique** : Newton-Raphson du second ordre sur un
  critère d'erreur de sortie, sensibilités par perturbation, conditionnement
  du second gradient surveillé par le rapport des valeurs propres extrêmes,
  paramètres insensibles détectés. C'est notre axe A (l'AD rend ces
  sensibilités exactes au lieu de les approcher) ;
- **système linéaire équivalent** : §7.1.

### 7.4 L'échelle d'inflow, et où nous sommes dessus

HOST publie sa hiérarchie complète (fig. 9), et elle sert de règle graduée :

```
   Meijer-Drees · Coleman · Blake & White      analytiques, quasi-statiques (TRIM)
   Pitt & Peters                    3 états    dynamique — le seuil des réponses hors-axe
   Peters–He / FiSuW          15 à 153 états   harmoniques × polynômes de Legendre (POD)
   anneaux tourbillonnaires · METAR            vortex lattice (trim seulement)
```

Vinkulum a l'**inflow uniforme dynamique** — le premier barreau de la colonne
dynamique. Le saut suivant est Pitt & Peters : trois états (λ₀, λ_s, λ_c), une
matrice de masse apparente [M] et une matrice de gains [L], trois EDO de plus.
HOST le désigne comme ce qui « améliore nettement la simulation des réponses
sur l'axe, en particulier pour le tangage en avancement », et le DLR montre
que sans lui les couplages roulis-tangage sont mal prédits. C'est peu de code
pour un effet nommé — il monte dans l'ordre de marche, juste après le trim.

> Réserve utile que le papier écrit lui-même sur ses modèles à états finis :
> « there is no representation of the vortices. Therefore the Blade Vortex
> Interactions can not be captured, thus this model is not suited for the
> aero-acoustic studies. » Le sillage n'est donc pas un raffinement de
> l'inflow : c'est une autre physique, et c'est pour ça que HART-II est un cas
> de BVI.

### 7.5 Ce que ça change à l'ordre de marche

1. ~~**trim**~~ — **fait le 3 septembre**, avec son banc rotor (`vinkulum.rotor`) ;
2. ~~**système linéaire équivalent périodique**~~ — **fait le 3 septembre**, et
   par la voie EXACTE plutôt que par l'approximation à coefficients constants :
   `vinkulum.floquet` propage la dynamique réelle sur une période. La limite
   gyroscopique déclarée le matin même est levée (§7.6) ;
3. **Pitt & Peters (3 états)** — le premier écart d'inflow nommé par HOST ;
4. ~~AD inverse~~ — **mode ADJOINT fait le 3 septembre sur l'équilibre**
   (`vinkulum.sensibilite`) : une résolution pour tous les paramètres, exact à
   0,00 % contre les différences finies globales, ×11,7 moins cher à six
   paramètres et le rapport monte avec N. **∂R/∂p est exact** depuis le même
   jour (linéarité de l'énergie de poutre), et **l'adjoint EN TEMPS aussi**
   (`vinkulum.adjoint_temps`, exact à 2e-11 % contre une référence
   analytique) ;
5. ~~poutre en rotation~~ — **faite le 3 septembre** (`vinkulum.rotation`) :
   courbe de Southwell, et le cas limite EXACT qui ne demande aucune donnée
   extérieure — raideur élastique → 0 fait de la poutre une chaîne, donc une
   pale articulée, dont ν vaut 1 (mesuré 1,0023). Reste la pale flexible
   branchée sur FRELON ;
6. ~~**pas adaptatif par résidu de demi-pas**~~ — **fait et MESURÉ le
   3 septembre : ×2,8 plus cher que le pas fixe à précision égale sur un cas
   oscillatoire.** Implémenté, opt-in, désactivé par défaut ; le résultat
   négatif est gardé par un assert (`bancs.adaptatif`). Voir
   `CAPITALISATION.md` §2 — l'estimateur est local, l'erreur d'un mécanisme
   peu dissipatif est de phase.

## 7bis. RecurDyn : où on peut le battre, où on ne le peut pas (3 septembre)

Cible ajoutée par Paul, avec la consigne « sans forcer ». Voici donc ce que
la mesure dit, et rien de plus.

### Le O(N) récursif : même ordre, par un autre chemin — MESURÉ

RecurDyn tient sa réputation d'une formulation RÉCURSIVE : le mouvement d'un
corps est écrit relativement à son parent, ce qui remplace une grande matrice
par une suite de petites. Sur une topologie ARBORESCENTE, le coût est O(N).

Vinkulum est en coordonnées ABSOLUES — l'approche que la récursive prétend
dépasser — mais avec un jacobien creux et un LU qui exploite sa structure
bande. Mesuré sur une chaîne de maillons, hors initialisation :

| corps | inconnues | contraintes | ms/pas |
|---|---|---|---|
| 30 | 270 | 90 | 0,375 |
| 300 | 2 700 | 900 | 4,117 |
| 1 200 | 10 800 | 3 600 | 12,93 |
| 2 400 | 21 600 | 7 200 | 24,13 |

**Exposant 0,93 de 270 à 21 600 inconnues, monotone.** L'ordre est donc le
même ; ce qui sépare les deux approches est la CONSTANTE, pas la complexité.
Et la récursive a une limite que la note omet : elle est O(N) sur un ARBRE ;
une boucle fermée doit être coupée et refermée par des contraintes, ce qui
ramène un système à résoudre. Un rotor à plateau cyclique — 71 contraintes en
boucles multiples — est exactement ce cas.

> Une correction au passage, parce qu'elle circule : le coût d'un solveur
> global n'est pas « exponentiel ». Il est O(N³) au pire en dense, et O(N) à
> O(N^1,2) dès qu'on exploite le creux — ce que les codes modernes en
> coordonnées absolues font tous.

### Ce que la mesure a trouvé en chemin

Le temps par pas n'était pas monotone : 69,6 ms à 5 400 inconnues, 17,7 à
10 800. **Un balayage non monotone sur un paramètre monotone n'est pas un
résultat.** La cause n'était pas le solveur : c'était la détection A PRIORI
des contraintes redondantes, un QR pivoté DENSE en O(m²n), payé une fois et
amorti — 5,9 s à 1 800 contraintes, soit sept fois le coût d'un pas. Son
seuil était posé à 2 000, ce qui le faisait sauter juste APRÈS le pic.
Ramené à sa mesure (900), le total tombe de 69,6 à 9,4 ms/pas, **×7,4**, et
le banc de redondance passe toujours (quatre-barres, 1,7e-13 m).

### MFBD : nos poutres font déjà ce que la note décrit — pour des poutres

« Sous-domaines nodaux, déformation non linéaire locale, forces réinjectées
dans l'équation du mouvement » : c'est exactement une poutre géométriquement
exacte, et c'est ce que vinkulum a. Le plan l'écrit depuis le premier jour et
la mesure du 3 septembre le confirme (raidissement centrifuge, cas limite
ν = 1 exact). **Craig–Bampton est chez nous la voie de la cellule RAIDE, pas
des pales** — pour la raison même que la note donne : réduire linéarise, et
une pièce qui tourne perd son raidissement.

Ce qui manque de notre côté et qu'il faut dire : des éléments **coque et
solide 3D**. Nous avons des poutres, rien d'autre.

### Le contact : la brique existe depuis le 3 septembre, le fossé reste

**Fait le jour même** (`vinkulum.contact`) : sphère contre demi-espace, loi de
Hertz F_n = k δ^e (1 + c δ̇) avec l'amortissement de Hunt–Crossley et un
Coulomb lissé par tanh. Tout est C¹, donc l'AD du noyau en donne la tangente
exacte comme pour n'importe quel élément.

| contrôle | mesuré |
|---|---|
| enfoncement statique, k δ^1,5 = mg | **288,6862 µm contre 288,6862** — 0,0000 % |
| force de contact contre le poids | 4,9050 N contre 4,9050 |
| restitution à c = 0 (conservatif) | **e = 1,000000** |
| pente à petit c | e ≈ 1 − 0,573·c·v₀, contre 8/15 = 0,533 dans la littérature |
| choc central de deux corps mobiles | vitesses de la théorie élastique, ΔP/P 4e-16 |
| bille sur plan incliné porté par un CORPS | **0,00 % sur six cas**, roulement ET glissement |

Le dernier contrôle est le plus complet du contact : il juge d'un coup la
force normale, le frottement de Coulomb régularisé, le COUPLE de frottement
(sans lui la bille glisserait au lieu de rouler, donc les bras de levier sont
jugés aussi) et le plan attaché à un corps mobile. Ses deux formules sont
exactes et leur transition est connue —
μ ≥ (2/7)·tan θ donne a = (5/7)·g·sin θ, sinon a = g(sin θ − μ cos θ).

Le choix de Hunt–Crossley plutôt qu'un amortisseur linéaire est le même que
celui de la « step function » d'ADAMS et de RecurDyn, obtenu sans paramètre de
seuil : en rendant l'amortissement proportionnel à δ^e, la force part de zéro
ET y revient, là où c·δ̇ saute à l'impact et TIRE au décollement.

### La BARRIÈRE IPC : l'interpénétration cesse d'être petite, elle disparaît

Li & al. 2020 remplacent la force qui croît AVEC l'enfoncement par un
potentiel qui DIVERGE quand la distance tend vers zéro :

    B(d) = −(d − d̂)² ln(d/d̂)   sur 0 < d < d̂,   0 au-delà

**Implémenté le 3 septembre** (`expo < 0` bascule la loi), et mesuré sur une
bille posée sous son poids :

| loi | résultat | force |
|---|---|---|
| Hertz k = 1e6 | **enfoncée** de 288,686 µm | 4,9050 N |
| Hertz k = 1e10 | **enfoncée** de 0,622 µm | 4,9050 N |
| IPC κ = 1e3, d̂ = 1e-3 | **séparée** de 227,701 µm | 4,9050 N |
| IPC κ = 1e5, d̂ = 1e-4 | **séparée** de 64,985 µm | 4,9050 N |

La force à l'équilibre est exacte des deux côtés ; ce qui change est le SIGNE
de la distance. Trois propriétés que la pénalité n'a pas : l'interpénétration
devient impossible au lieu d'être petite ; le support est COMPACT, donc le
conditionnement n'est pas dégradé loin du contact ; et B est C², donc la force
est C¹ et l'AD passe au travers — la barrière ne coûte rien à l'axe A.

### Le CCD complète la barrière : la traversée devient impossible

Une barrière empêche l'APPROCHE ; elle n'empêche pas qu'un pas grossier fasse
passer une sphère rapide de l'autre côté, les deux extrémités du pas restant
saines pendant que la trajectoire, elle, ne l'est pas.

**Fait le 3 septembre.** Pour nos primitives le test est analytique : la
distance entre centres vaut |Δc₀ + t·Δv|, dont le minimum sur [0,1] se calcule
en fermé. Un pas qui traverse est REJOUÉ plus court.

```
   bille à 50 m/s contre une autre, pas grossier 1e-3
   sans CCD   x_a +0,5500  x_b +0,0500              → elle a TRAVERSÉ
   avec CCD   x_a +0,0568  x_b +0,5432              → elle rebondit
              v_a +1,59    v_b +48,41                 (choc quasi élastique)
```

Et quand la barrière est trop faible pour arrêter la bille, le CCD **REFUSE**
au lieu de laisser passer — un arrêt vaut mieux qu'une traversée silencieuse.

**Ce qui n'est toujours PAS repris d'IPC** : la formulation complète fait de
chaque pas une MINIMISATION d'un potentiel incrémental, et le filtre CCD y
agit dans la recherche linéaire de cette minimisation. Ici il agit au niveau
du PAS, après coup, avec rejeu. La garantie est la même sur le résultat, le
chemin pour l'obtenir ne l'est pas — et un pas rejoué coûte, là où une
recherche linéaire filtrée ne rejoue rien.

**L'APPARIEMENT est fait aussi** (3 septembre) : grille de hachage spatial,
les paires se découvrent et disparaissent seules, et la structure symbolique du
LU creux est invalidée quand la liste change.

| N sphères | paires | appariement | naïf O(N²) |
|---|---|---|---|
| 50 | 0 | 0,038 ms/pas | 0,038 |
| 800 | 4 | 0,613 | 9,65 |
| 3 200 | 22 | 2,466 | 154,4 |

**Exposant mesuré 1,005** — linéaire, contre 2,0 pour le balayage naïf, soit
×62 à 3 200 sphères. Et l'appariement automatique rend le même résultat qu'une
paire déclarée à la main, au dernier chiffre.

### Tout ensemble, à l'échelle

Billes lâchées sur un sol, contacts découverts par la grille, barrière IPC,
frottement actif :

| billes | paires simultanées | ms/pas | enfoncement max |
|---|---|---|---|
| 30 | 18 | 0,065 | **−267,38 µm** |
| 100 | 108 | 0,449 | **−246,92 µm** |
| 300 | 394 | 3,002 | **−191,43 µm** |

**Pas une seule paire n'a pénétré**, sur 394 contacts simultanés — c'est la
propriété d'IPC, et un enfoncement MOYEN ne l'aurait pas dit : c'est le
maximum sur toutes les paires qui est asserté. L'exposant 1,67 du temps par
pas suit le nombre de PAIRES, qui croît plus vite que N pendant l'empilement ;
l'appariement lui-même est mesuré à 1,00.

**Ce qui reste du fossé** : nos primitives sont des SPHÈRES. RecurDyn apparie
des géométries quelconques — surfaces Parasolid exactes, hiérarchies de volumes
englobants — là où nous demandons qu'un solide soit approché par des sphères.
Et FRELON fait toujours ses collisions **hors du solveur**, par booléens OCCT
sur des poses résolues : le pont entre les deux reste à faire. RecurDyn en a fait son cœur de métier —
hiérarchies de volumes englobants sur géométrie Parasolid exacte, force de
pénalité régularisée avec amortissement et frottement de Coulomb, choisie
précisément parce qu'un multiplicateur strict rend le résidu discontinu et
fait diverger le solveur au choc.

C'est le manque le plus sérieux, et il est structurel : FRELON fait
aujourd'hui ses collisions **hors du solveur**, par booléens OCCT sur des
poses résolues. C'est correct pour un contrôle de non-interpénétration, ce
n'est pas de la dynamique de contact.

### Le compte honnête

| | vinkulum | RecurDyn |
|---|---|---|
| coût par pas, grands systèmes | exposant 0,93 mesuré | O(N) revendiqué |
| poutres géométriquement exactes | oui, validées | oui (MFBD) |
| coques et solides flexibles | **non** | oui |
| contact | **sphère/plan, Hertz + Hunt–Crossley + Coulomb lissé**, tangente exacte par AD | oui, sur géométrie Parasolid exacte, BVH, tout appariement |
| gradient exact (AD, adjoint) | oui, mesuré à 2e-11 % | non |
| déterminisme au bit | oui, mesuré | non documenté |
| bancs publics rejouables | oui | non |
| couplage aéro monolithique | oui | co-simulation |

**On ne le bat pas sur son terrain, et le dire est la condition pour être
cru sur le reste.** Ce qui se tient : l'ordre de complexité, les poutres, et
tout l'axe de la différentiabilité — où il n'a rien.

## 7ter. Les trois axes de recherche 2026, passés à la mesure (3 septembre)

Notes de recherche transmises par Paul : intégrateurs variationnels
symplectiques, contact non-lisse par complémentarité (Project Chrono NSC),
décomposition de domaine asynchrone. Chacun est réel ; chacun demande une
nuance que la mesure fournit.

### Noether discret : la critique porte, et nous la chiffrons contre nous

Un corps libre sans couple extérieur a un moment cinétique EXACTEMENT
constant. Un intégrateur variationnel le conserve à la précision machine ; le
nôtre le conserve à l'ordre du schéma. Mesuré sur 5,6 tours d'une toupie
libre :

| ρ∞ | h | dérive de \|L\| | de sa direction |
|---|---|---|---|
| 0,9 | 1e-3 | −1,42e-08 | 9,11e-04° |
| 0,9 | 1e-4 | −1,66e-10 | 9,16e-06° |
| 0,5 | 1e-3 | **+7,63e-07** | 1,33e-03° |

Deux faits, et le second nous est défavorable : la dérive suit l'ordre 2
(h/10 la divise par 86), **et elle dépend de ρ∞ — ρ∞ = 0,5 dérive 54 fois
plus que 0,9.** La dissipation numérique coûte donc bien en conservation
d'invariant, exactement comme la littérature variationnelle le reproche.

À notre échelle l'effet reste minuscule : 1,4e-08 sur 5,6 tours, soit ~0,3 %
sur un MILLION de tours. C'est là — mécanique céleste, dynamique moléculaire,
transmissions sur des horizons énormes — que le variationnel gagne, et pas
sur un vol de quelques minutes.

> **Une note de méthode, parce qu'elle a failli coûter une conclusion
> inverse** : un premier essai concluait « ρ∞ n'y est pour rien », sur une
> mesure où `rho` n'était pas passé au solveur. Les quatre lignes identiques
> auraient dû alerter avant la conclusion. L'assert du banc vérifie désormais
> que ρ∞ CHANGE quelque chose — un contrôle qui ne peut pas distinguer ses
> paramètres ne contrôle rien.

### Symplectique : la conservation « au bit près » n'existe pas

Un intégrateur symplectique **ne conserve pas l'énergie exactement**. Il
conserve une énergie MODIFIÉE, et l'énergie vraie oscille dans une bande
bornée — c'est le résultat d'analyse rétrograde qui fait toute sa valeur, et
c'est autre chose que « préservée au bit près ». Le moment cinétique, lui, est
exactement conservé si le schéma est équivariant (Noether discret).

Et pour NOUS, la question ne se pose pas dans ces termes, parce que deux
mesures la ferment :

- **la dissipation n'est pas un défaut à supprimer, c'est une nécessité** : en
  index 3 avec contraintes, le trapèze non dissipatif est INSTABLE
  (Cardona–Géradin 1989), et le noyau refuse ρ∞ = 1 sur cette base — mesuré,
  divergence à t = 16 s ;
- **la dérive qu'un symplectique corrigerait est déjà négligeable** : ΔE/E =
  9,4e-06 sur cent périodes du pendule contraint (banc `projection_vitesse`).

Ce qui resterait à gagner est donc borné par 1e-5, contre une reformulation
complète. Le rapport n'y est pas.

### Contact non-lisse : le gain est réel, et son prix est notre axe A

Le reproche à la pénalité — « elle force h à 1e-7 » — est **exagéré d'un
facteur cent, et vrai dans son principe**. Mesuré sur le choc élastique de
deux billes, en faisant varier la raideur de contact :

| k (N/m^1,5) | durée du choc | h admissible | pas dans le choc |
|---|---|---|---|
| 1e6 | 19 397 µs | 1e-4 | 194 |
| 1e8 | 3 074 µs | 1e-5 | 307 |
| 1e10 (acier) | 487 µs | **1e-5** | 49 |
| 1e12 | 77 µs | 1e-6 | 77 |

Ce qui commande n'est pas k mais la **durée du choc**, qui décroît en
k^−0,4 : il faut cinquante à trois cents pas dedans, quelle que soit la
raideur. À de l'acier sur acier, h = 1e-5 suffit — pas 1e-7.

**Mais le reproche tient sur le fond** : un mécanisme à 50 Hz se contente de
h = 1e-3 sans contact et exige 1e-5 avec, soit **un facteur cent** de pas
imposé par le contact et non par la dynamique. C'est précisément ce que la
complémentarité (Moreau, LCP/CCP, Chrono::NSC) supprime.

> **Son prix touche l'axe A, et j'ai d'abord écrit cet arbitrage trop
> tranché.** J'avais posé « un LCP ne donne pas de gradient ». C'est faux :
> le sous-différentiel de Clarke et la différentiation implicite des
> conditions KKT actives donnent un gradient, valide **presque partout**, et
> c'est un domaine actif (physique différentiable). Ce qui reste vrai, et
> qu'il faut dire à sa juste mesure :
>
> · le gradient d'un LCP est **discontinu aux changements d'ensemble actif** —
>   il existe presque partout, mais « presque partout » n'est pas partout, et
>   une descente qui traverse une transition ne voit pas venir le saut ;
> · le sous-différentiel de Clarke est un **ensemble**, pas un vecteur : un
>   algorithme doit en choisir un élément, et les garanties de convergence du
>   cas lisse tombent ;
> · la littérature de la simulation différentiable rapporte que les gradients
>   à travers un choc **rigide** sont souvent peu informatifs, au point que
>   plusieurs travaux réintroduisent délibérément du lissage.
>
> **NOTRE contact, lui, est mesuré différentiable à travers un choc complet.**
> Balayage de la raideur autour de 1e10 sur un rebond entier :
>
> ```
>    k de 9,80e9 à 1,02e10 (±2 %)   ·   z final de 0,073813743 à 0,073816641 m
>    dérivée 6,98e-15 … 7,43e-15 m/(N/m^1,5)   monotone, variation 6,4 %
> ```
>
> La sortie est C¹ et bien conditionnée **de l'autre côté du rebond** — la
> variation de 6,4 % est la courbure, pas du bruit. C'est ce qui autorise
> l'adjoint à traverser le contact.
>
> L'arbitrage honnête est donc : la complémentarité achète de grands pas et
> paie en régularité du gradient ; la pénalité fait l'inverse. Aucune des deux
> n'interdit la différentiation — l'une la rend franche, l'autre la rend
> conditionnelle.

Ce qui se ferait sans rien perdre : garder la pénalité et lui donner un **pas
adaptatif local** — c'est-à-dire l'axe 3, la décomposition en temps.

### Décomposition en temps : la version simple est faite, et elle paie

C'est l'idée qui rend la pénalité économique : micro-pas là où il y a choc,
macro-pas ailleurs. **Fait le 3 septembre** sous sa forme la plus directe —
le pas global suit l'ÉCHELLE ACTIVE, resserré dès qu'une paire approche du
contact, relâché sinon.

```
   bille qui rebondit, k = 1e10 (acier), 1,2 s
   pas fixe 1e-5        z 0,021372 m   0,151 s
   piloté 1e-4 / 1e-5   z 0,021372 m   0,053 s   ×2,8   écart 0,0 µm
```

**Et l'estimateur n'a rien d'empirique** — c'est ce qui le sépare du résidu de
demi-pas, refusé le même jour : on CONNAÎT la distance à la collision, donc on
sait exactement quand resserrer. Le pas adaptatif global coûtait ×2,8 ; celui-ci
en rapporte 2,8, sur le même noyau, parce que son critère est exact au lieu
d'être heuristique.

Le gain décroît quand le pas libre grandit (×2,8 à 1e-4, ×1,7 à 1e-3) : la
marge d'approche doit grandir avec lui, donc la zone à pas fin s'élargit.
C'est un compromis explicite, pas un réglage caché.

**Ce qui reste de l'axe** : le sous-cyclage vrai, où chaque sous-domaine a son
propre pas et où les cœurs n'attendent pas. Ici tout le système suit le pas du
contact ; sur un mécanisme où le contact est local, c'est encore trop.

### Une correction sur l'intégration intrinsèque

La note affirme qu'un intégrateur de Lie sur SE(3) élimine « les DAE d'indice
3 au profit d'un schéma d'indice 0 ou 1 ». **Cela vaut pour la contrainte de
NORMALISATION** — un quaternion dérive, une rotation mise à jour sur le
groupe ne dérive pas, et c'est exact et c'est ce que vinkulum fait depuis le
premier jour.

Mais **cela ne supprime aucune contrainte de LIAISON.** Les 71 contraintes de
la tête S2 — pivots, distances, engrenage — restent des équations algébriques
d'indice 3, quelle que soit la façon dont on paramètre les rotations. Le
raccourci confond deux choses : la géométrie de la variété de configuration
d'un corps, et les contraintes qui relient les corps entre eux.

## 7quater. Cinq axes de plus, et ce que la mesure en retient (3 septembre)

### S-MOR : la critique porte, et nous la CHIFFRONS

Une réduction de Galerkin fige la raideur de l'ÉTAT où elle est faite. Nous
l'avions déclaré ; voici ce que ça coûte, mesuré sur la poutre en rotation :

| Ω (rad/s) | f complet | f du réduit-à-l'arrêt | erreur |
|---|---|---|---|
| 20 | 4,723 | 3,266 | **−30,9 %** |
| 100 | 16,984 | 3,266 | **−80,8 %** |
| 200 | 34,234 | 3,266 | **−90,5 %** |

Sur le rotor FRELON (Ω = 335 rad/s) l'erreur serait pire encore. **C'est
exactement pourquoi une pale garde ses poutres géométriquement exactes** et
pourquoi Craig–Bampton vise la cellule raide. La réduction symplectique sur
variété NON LINÉAIRE, qui lèverait le verrou, n'est pas faite — et elle ne
lève pas le problème du raidissement centrifuge, qui n'est pas un problème de
symplecticité mais de dépendance de la raideur au régime.

### DDG / Discrete Elastic Rods : une alternative, pas un dépassement

Repère de Bishop, courbure discrète aux sommets, torsion réduite à un scalaire
par arête : c'est élégant et c'est la bonne formulation pour une tige MINCE.
Le gain annoncé « ×100 à ×1000 » se compare à des éléments finis 3D dégénérés,
pas à une poutre géométriquement exacte — dont c'est justement l'argument
depuis Simo–Vu Quoc. Notre poutre coûte 1,57 ms/pas à 20 éléments et rend le
premier mode à +0,01 %. À rouvrir le jour où un câble entre au modèle.

### AVMI : notre pas piloté en est la version dégénérée

Le multi-rythme variationnel donne à chaque sous-domaine son propre pas, sans
perte d'énergie à l'interface. Notre pas piloté par le contact fait la même
chose avec UN pas global, ce qui est plus grossier — et rapporte quand même
×2,8. Le vrai AVMI est le chantier de performance suivant.

### Koopman : réserve, et elle est de fond

Pour un système hamiltonien, le générateur est bien anti-auto-adjoint et le
semi-groupe unitaire : le spectre est sur l'axe imaginaire. Mais la
**diagonalisation exacte vit en dimension INFINIE** ; toute mise en œuvre
(EDMD) tronque sur une base finie d'observables, et devient une
approximation dont l'erreur dépend du choix de cette base. « Résoudre par
déphasages spectraux sans assembler de matrice tangente » décrit le cas
linéaire ou quasi-intégrable, pas un mécanisme à contacts et butées. À
surveiller, pas à adopter.

### Transport optimal / JKO : hors périmètre, et pour une raison de fond

Le flot de gradient dans l'espace des mesures rend une **densité**, pas des
trajectoires individuelles. C'est le bon outil pour du régolithe à un million
de grains dont on veut le comportement d'ensemble ; c'est le mauvais pour un
rotor dont on veut le pas de chaque pale.

## 8. La stabilité en rotation est tranchée, par Floquet (3 septembre)

`modes_complexes` refusait tout état en rotation, et la mesure avait dit
pourquoi : la raideur tangente en repère fixe y est dominée par des termes en
ω²(J_t − J_a) valant mille fois la raideur cherchée. **La sortie n'était pas
de mieux projeter, c'était de ne pas linéariser en repère fixe.**

HOST fait la version approchée de cette idée — composantes périodiques de
l'état, puis coefficients constants. Floquet en est la version exacte : on
propage la dynamique réelle sur une période et on lit les valeurs propres de
la matrice de transition. Aucune linéarisation de champ, donc aucune
annulation catastrophique.

| cas | ce qui est mesuré |
|---|---|
| équilibre amorti | f 1,97737 Hz · ζ 0,15000 — **identique** à la théorie et à `modes_complexes` |
| pale battante, Ω 109 rad/s, ν_β 1,12 | σ à **0,01 %** · ζ à **0,13 %** · f repliée à 1,19 % |
| rotor 4 pales, aéro et inflow, μ = 0 et 0,15 | **stable**, \|λ\| max 0,16 ; ζ de battement 0,29 et 0,31 contre γ/16 = 0,375 |

Deux propriétés de Floquet que le module publie au lieu de les taire, parce
qu'elles se paient sinon : **l'exposant n'est connu que modulo 2π/T**, donc à
T égal à la période propre la phase vaut exactement 2π et un oscillateur sain
sort à « f = 0, ζ = 1 » ; et **ζ ne se calcule pas sur la fréquence repliée**
— sur la pale à ν_β = 1,12 il sortait à 0,427 au lieu de 0,050. Ce qui reste
non ambigu est la partie RÉELLE de l'exposant, et le verdict de stabilité n'en
dépend pas.

Ce que ça coûte : 2m intégrations d'une période pour m directions. Le rotor
complet, huit directions, aéro active : **1,1 s**.

## 9. Invariants — ce qui est fait, et le chemin SOTA pour la suite (3 septembre)

Mesuré ici, dans les deux sens, et gardé par des bancs :

| | domaine | verdict |
|---|---|---|
| projection sur le **moment** | symétrie SO(3) globale (corps libres, liaisons entre corps) | **retenue** — dérive 3,6e-10 → 4,1e-20 /s, ordre 2,00 préservé (libre ET contraint), dissipation HF identique au 8ᵉ chiffre, \|Φ\| inchangé |
| projection sur l'**énergie** | systèmes rigides **sans liaison** | retenue sur ce domaine seul — E à 1e-17, ordre 2,02 ; sous liaison `G δu = 0` rend la correction transverse et l'erreur cesse de converger (ordre ~0) |
| les **deux ensemble** | — | **impossible par projection** : pour un corps rigide `∂E/∂ω = ωᵀ·∂L/∂ω` exactement, rang 3 sur 4 |
| bouton **ρ∞** | — | ne répare rien : monté il dégrade les systèmes contraints, Newton lâche à 0,999 (index 3, Cardona–Géradin) |

### Le chemin qui reste, et il est nommé dans la littérature

La limite est structurelle : imposer E **et** L demande de corriger aussi la
configuration, donc de changer la FORMULATION. Trois familles, par ordre de
pertinence pour ce projet :

1. **Schémas à décroissance d'énergie contrôlée** — Bauchau & Bottasso, *On the
   design of energy preserving and decaying schemes for flexible, nonlinear
   multi-body systems* (CMAME 169, 1999) et la série qui suit. C'est
   **exactement** le dilemme mesuré ici — dissiper les hautes fréquences des
   contraintes SANS perdre les invariants — et c'est traité pour les
   **voilures tournantes**, notre cas. À lire en clean-room avant toute
   implémentation.
2. **Schémas énergie-moment par gradient discret** — Simo & Tarnow (1992),
   Gonzalez (1996), Betsch & Steinmann : remplacer ∇V(q) par un gradient
   discret vérifiant `G·(q₁−q₀) = V₁−V₀`. Conservation exacte des deux, mais
   **aucune** dissipation HF : à combiner avec (1), pas à substituer.
3. **Intégrateurs variationnels de Lie** — Lee, Leok & McClamroch : Noether
   discret exact par construction, symplectiques. Même réserve : sans
   dissipation HF, l'index 3 direct est instable (mesuré ici à ρ∞ 0,999).

⇒ **La cible est (1)**, et (2) en fournit la brique conservative. La projection
livrée aujourd'hui n'est pas un substitut : c'est la réparation qui coûte
zéro changement de formulation, et son domaine est écrit.


## 10. État au 3 septembre (soir) — ce que la passe de maturité a livré

| brique | verdict mesuré |
|---|---|
| projection sur le **moment** | dérive 3,6e-10 → 4,1e-20 /s ; ordre 2,00 et dissipation HF intacts |
| projection sur l'**énergie** | domaine restreint aux systèmes sans liaison, mesuré des deux côtés |
| **sensibilité de Floquet** ∂σ/∂p | 0,004 % contre l'analytique, **sans appariement de modes** |
| **robustesse** | 14 modèles absurdes, 14 refus typés ; 3 défauts trouvés (dont un blocage infini) |
| **Pitt–Peters** 3 états | Froude 0,0000 %, τ théorique −0,70 %, Glauert −0,02 % |

### Ce qui reste, par ordre de valeur

1. **Schéma énergie-moment dans la formulation** (Bauchau–Bottasso 1999) — la
   seule voie pour conserver E *et* L, et elle est écrite pour les voilures
   tournantes. Le cadre est au §9.
2. **Aérodynamique instationnaire de section** (Theodorsen linéaire, puis
   Leishman–Beddoes) : le décrochage dynamique n'est pas couvert, et c'est ce
   qui borne un rotor en manœuvre.
3. **Contact sur géométrie quelconque** — aujourd'hui sphères et plans. C'est
   le plus gros écart fonctionnel avec un solveur généraliste du commerce.
4. **Sillage libre** (vortex particles / lifting line) : Pitt–Peters donne le
   gradient 1/rev, pas l'interaction pale-tourbillon.
5. Adjoint en temps **non linéaire** (le cas linéaire est fait), coques et
   solides 3D, sous-cyclage multi-rythme.


## 11. Plan v0.5 — ACHEVÉ (3 septembre)

Les six briques nommées au départ sont livrées, chacune avec son banc et un
juge qui ne partage rien avec l'implémentation :

| brique | juge indépendant | écart |
|---|---|---|
| **capsule** | géométrie segment–segment exacte, 3 régimes | 1e-10 |
| **boîte (OBB)** | clamp exact, face/arête/coin + rotation du corps | 4e-15 |
| **câble** | chute libre contre retenue à la longueur | −0,513 pour 0,50 |
| **butée** | √(2E/k) sur 4 décades de raideur | 0,08 % |
| **vis–écrou** | inertie ramenée (I + m(pas/2π)²)ω = Iω₀ | 0,0001 % |
| **assemblage** | Φ et Φ̇ depuis une pose approchée, norme minimale | 1e-12 en 3 it |

Ce que le plan nommait et qui reste **délibérément dehors**, avec son motif :

- **came** — une came n'est pas une primitive, c'est un PROFIL. Elle demande
  une courbe paramétrée et sa normale, donc une représentation géométrique que
  le noyau n'a pas encore ; l'ajouter en dur pour un profil circulaire serait
  une capsule déguisée. À faire avec le pont B-Rep, pas avant.
- **cylindre** — la capsule le couvre à ses extrémités près (un cylindre est
  une capsule à bouts plats). Le gain est un cas de bord, le coût une branche
  de plus dans une distance déjà à trois régimes : mesuré comme non rentable.
- pénétration profonde dans une boîte, hiérarchie de volumes englobants,
  géométrie quelconque.


## 12. Cap « noyau généraliste » — les six manques comblés (3 septembre)

Le 3 septembre, à la question « le noyau est-il mature ? », la réponse fut
« sur son domaine, oui ; sans qualificatif, non », avec six manques nommés.
Ils sont comblés, chacun avec un juge qui ne partage rien avec lui :

| manque nommé | ce qui le comble | juge |
|---|---|---|
| ni coques ni solides 3D | **superélément** corotationnel | rotation rigide 180° → 7e-18 |
| — | **coque ACM** + pont vers le multicorps | poutre analytique −0,10 % ; pont 2,6 % |
| pas de géométrie quelconque | **maillage + BVH** | facettage O(h²) ; exposant 0,44 |
| adjoint en temps linéaire | **adjoint non linéaire** | Duffing 66 % cubique, 9,6e-7 % |
| pas de doc d'API | **`docs/API.md` générée**, 123 entrées | contrôle de fraîcheur dans la CI |
| ni CI ni packaging | **CI** 7 étapes + roues 3.11→3.13 | fmt, clippy -D warnings, tests, bancs |

### Ce qui reste, et qu'aucune brique ne lève

1. **Un utilisateur tiers**, sur un modèle qu'il apporte. Trois problèmes
   posés ailleurs ont été ajoutés depuis — Andrews (IFToMM, jugé par MBDyn sur
   la même définition, 0,0050 %), Bennett, Kapitza — et **Andrews a trouvé une
   faute qu'aucun banc interne n'aurait vue** (un ressort modélisé en câble
   alors qu'il est comprimé : 2,4× trop vite, avec des contraintes tenues à
   1e-17). Le biais est donc mesuré et entamé, pas supprimé : les modèles
   viennent d'ailleurs, la mise en œuvre reste de l'auteur du solveur.
2. **Schéma énergie-moment dans la formulation** (Bauchau–Bottasso) — §9.
3. **Cylindre**, pénétration profonde dans une boîte, came (§11).
4. **Checkpointing** de l'adjoint non linéaire (mémoire en √N).
5. Sillage libre, décrochage dynamique, sous-cyclage multi-rythme.


## 13. Modèles étrangers — ce qui est pris, ce qui est ouvert (3 septembre)

`~/src/mbdyn/tests/benchmarks/` (livré avec MBDyn) contient huit cas posés par
des tiers. Ils sont la source de modèles étrangers la plus proche, et lire
leurs DONNÉES n'est pas copier du code — la distinction clean-room du projet
s'applique telle quelle.

**Pris** : `andrewssqueezer` (Schiehlen 1990, IFToMM) — 7 corps, 41
contraintes, ressort raide, jugé par MBDyn sur la même définition à
**0,0050 %** après 1941° parcourus, et il a trouvé une faute qu'aucun banc
interne n'aurait vue.

**Commencé, NON publié — et voici exactement où il en est** : `princeton`,
l'expérience de poutre de Dowell & Traybar (1975), le juge de référence des
poutres géométriquement exactes et le seul des huit dont la référence soit
**expérimentale**.

Deux obstacles franchis :
· Octave est absent : `.set`, `.nod`, `.elm` régénérés à la main depuis
  `princeton_gen.m` ;
· **la convention du repère était mal lue** — dans le générateur, la ligne
  `2, 0., cos(THETA), sin(THETA)` est COMMENTÉE et la ligne active est
  `2, 0., 1., 0.` : le repère de la poutre NE tourne pas, seule la CHARGE
  tourne. Avec le bon `.ref`, MBDyn donne le comportement attendu (u₂ de 0 à
  0,146 m quand la charge passe du plan fort au plan faible, rapport 13,6 ≈
  EJY/EJZ).

Ce qui reste, et pourquoi ce n'est pas publié :

| θ | MBDyn (u₂, u₃) | vinkulum | écart |
|---|---|---|---|
| 15° | (+0,04213, +0,01063) | (+0,04220, +0,00738) | 3,3 mm |
| 45° | (+0,10904, +0,00889) | (+0,10885, +0,00510) | 3,8 mm |
| 90° | (+0,14626, 0) | (+0,14597, 0) | 0,3 mm |

**u₂ — le déplacement dominant — est retrouvé à 0,2 %.** Tout l'écart est sur
u₃, le petit.

Et le raffinement de maillage tranche, mais pas comme je l'avais d'abord
écrit. Sur deux points je concluais « raffiner AGGRAVE » (2,62 % à 10
éléments, 3,82 % à 20) ; le troisième dit **2,66 % à 40**. La suite est donc
**NON MONOTONE**, et un balayage non monotone sur un paramètre monotone n'est
pas un résultat — c'est du bruit. (Le dépôt a écrit cette leçon une première
fois sur la poutre de queue de FRELON ; c'est la seconde.)

Ce que ça élimine et ce que ça désigne : ni la discrétisation (elle
convergerait), ni la poutre (u₂ est juste à 0,2 %). Il reste la RELAXATION —
le cas est résolu ici en amortissant jusqu'à l'immobilité, donc chaque
maillage s'arrête à un endroit différent de sa décroissance. **La piste est un
solveur STATIQUE** (`residu_statique` existe déjà), pas un élément.

**Ouverts** : `6barmech` (surcontraint), `fourbar`, `lateralbuckling`,
`multibarmech`, `rotatingshaft`, `srskm`.


## 14. Princeton hors plan — deux explications réfutées, et ce que la mesure dit

Le solveur statique (§13 de `verification`) résout le cas plan à **0,27 %** de
MBDyn, contre 2,6 % par relaxation. Hors plan, il plafonne à ‖r‖ ≈ 0,7 N.

**Réfuté** — « c'est le bruit de relaxation » : le solveur statique ne relaxe
pas, et le cas plan y gagne cent fois.
**Réfuté** — « c'est le conditionnement » : l'équilibrage de Jacobi a été
implémenté (il reste, il ne nuit pas) et ne déplace pas le plafond d'une
décimale.

**Réfuté** — « un seul module de cisaillement pour une section qui en a deux
(GAY ≠ GAZ) » : les deux valeurs donnent le MÊME résultat à θ = 90°.

**LA QUATRIÈME HYPOTHÈSE ÉTAIT LA BONNE, et le défaut était dans mon BANC** :
il n'amortissait qu'autour de **z**. Un mouvement hors de ce plan n'était donc
pas amorti du tout, et la relaxation s'arrêtait sur un état qui n'était pas un
équilibre — résidu 0,2 à 0,7 N — en rendant des déplacements d'apparence
plausible. Amorti sur les trois axes, **les sept angles passent** :

| θ | 0° | 15° | 30° | 45° | 60° | 75° | 90° |
|---|---|---|---|---|---|---|---|
| écart à MBDyn | 0,02 % | 0,09 % | 0,14 % | 0,17 % | 0,19 % | 0,20 % | 0,20 % |
| résidu (N) | 2e-7 | 4e-8 | 5e-7 | 3e-8 | 5e-7 | 1e-7 | 9e-8 |

**Princeton est donc validé de bout en bout**, et le noyau n'avait aucun
défaut. *Un banc qui n'amortit qu'un axe ne relaxe qu'un plan* — et le seul
contrôle qui l'a dit est le RÉSIDU, jamais l'écart à la référence. Le
contre-contrôle est gardé : à un seul axe, le même cas doit encore s'arrêter à
0,65 N.
