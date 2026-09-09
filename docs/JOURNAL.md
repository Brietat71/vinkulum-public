# Journal de vinkulum

Le récit brique par brique, du jalon 1 à v0.3.2. Sorti du README le 4 septembre :
un README décrit ce qu'un noyau EST, un journal ce par quoi il est passé — et
les deux ne vieillissent pas au même rythme. Les chiffres de ce fichier sont des
**traces datées** ; ce qui fait foi est la sortie de `verification`, `campagne`
et `bancs`, qui les recalculent.

## Vérification approfondie du jalon 1 (3 septembre, soir) — `python -m vinkulum.verification`, 9 s

La démo juge des grandeurs intégrales ; ici on juge ce que la théorie PROMET,
contre des références indépendantes (EDO minimales intégrées à 1e-12), et
deux défauts sont sortis de la mesure, corrigés le soir même.

```
   1. ordre de convergence (pendule 30°, 2 périodes)   h = T/100 … T/800
      position  2,13e-3 → 3,35e-5   ordres 1,99 · 2,00 · 2,00
      vitesse   9,21e-3 → 1,44e-4   ordres 2,00 · 2,00 · 2,00
      tension λ 5,48e-3 → 8,58e-5   ordres 2,00 · 2,00 · 2,00     (Arnold–Brüls 2007 : ordre 2 en q ET λ — tenu)
   2. contrainte en VITESSE r·v/(L|v|)   5,7e-4 → 9,5e-6   ordre 2  — NON nulle : Φ tenu, Φ̇ pas
   3. toupie contre l'EDO de Lagrange, axe complet, nutation comprise, 1 s
      h 4e-4 / 2e-4 / 1e-4 : 4,8e-3 / 1,2e-3 / 3,0e-4 rad   ordres 1,99 · 2,00
   4. déterminisme : 1 fil et tous les fils IDENTIQUES au bit
   5. ρ∞ = 1 REFUSÉ ; 0,95 / 0,9 / 0,5 : ΔE/E −1,4e-8 / −9,5e-8 / −2,3e-5 sur 10 périodes
```

**Ce que la mesure a corrigé.**
- *Tolérance de Newton.* À `tol` 1e-9 l'ordre du multiplicateur tombait à
  1,77 à h = T/800 : l'erreur de Newton passait devant celle du schéma. À 1e-13
  il vaut 2,00. Défaut porté à **1e-12**. Et le critère devient **par blocs** :
  les lignes de contrainte sont mises à l'échelle 1/(β·h²) — 4e8 à h = 1e-4 —
  et un Φ au plancher d'arrondi (1e-17 m) pesait 4e-9 dans la norme globale :
  Newton « ne convergeait pas » sur un résidu nul. Dynamique jugée en relatif
  aux forces, Φ en absolu dans ses unités (1e-11).
- *ρ∞ = 1 est instable.* Pendule à 30°, h = T/400 : Newton diverge à t ≈ 16 s
  (résidu 1,9e5) après des oscillations croissantes de λ. C'est le résultat de
  **Cardona–Géradin (1989)** : « la règle du trapèze de Newmark est
  inconditionnellement instable en présence de contraintes » — ce qui a motivé
  HHT-α puis l'α-généralisé. Le noyau **refuse** ρ∞ = 1 ; la note Aether qui
  vantait « symplectique, plus d'amortissement numérique » proposait donc
  exactement le schéma instable en index 3.

**Ce que la mesure confirme, et qu'il faut savoir.** L'index 3 direct laisse
la contrainte en vitesse violée à O(h²) (9,5e-6 à T/800) — c'est l'argument
réel des formulations GGL / index 2 stabilisé (Adams SI2). Sur une tête rotor
à 3 000 tr/min ça se mesure à h ; si un cas le demande, GGL entre en option,
sur mesure. La rotation relative de π reste une solution parasite du résidu
vee(skew(E)), hors bassin de Newton — documentée, pas gardée par un test.

## Jalon 2 (3 septembre, nuit) — l'engrenage, et la tête rotor S2 de FRELON comme premier banc

**Engrenage** (`Noyau.engrenage`) : θ_a − rapport·θ_b = 0 autour de deux axes
du porteur, angles DÉROULÉS (la branche d'atan2 la plus proche du dernier pas
accepté ; plus d'un quart de tour par pas → erreur nette). C'est la liaison
que MBDyn n'a pas ; sa réaction est le couple de denture, en sortie. Jugé :
deux roues sur pivots, rapport −4, couple constant sur le pignon —
ω_a/ω_b = −4,000000000, α_a = τ/(J_a + J_b/16) à 1e-9 (l'inertie de la roue
ramenée par le carré du rapport), 4,4 tours déroulés, λ = la part du couple
qu'absorbe la roue.

**La tête S2 de FRELON** (`cad/vinkulum_s2.py` dans FRELON — le modèle est à
FRELON, le noyau est ici) : 12 corps, 71 contraintes, dont l'étage 4:1 en
engrenage avec le rotor moteur SUR le pignon. À commande tenue, mât en
rotation imposée, servos en rotation imposée, les pas de pale retombent sur
la chaîne analytique `timonerie_s2` (7 inconnues) à **< 1e-5°** sur 96 images,
|Φ| 1,5e-12 ; collectif pur idem. 8 000 pas de 2e-4 s en 98 s (143 inconnues,
jacobien par différences finies en parallèle — c'est le poste que l'AD
ramènera).

Ce que ce banc ne fait pas encore, déclaré : gouverneur et servos en COUPLE
(PD saturé), aéro c81 — le jalon suivant, jugé contre le vol MBDyn.

## v0.2.0 (4 septembre) — jacobien exact par AD, Newton simplifié, redondance : ×20 à ×34

Le poste n°1 du plan est fait. `src/ad.rs` : nombres duaux à N tangentes,
emboîtables, algèbre 3×3 générique, exp(SO(3)) en série près de zéro, testés
contre gradient et Hessienne analytiques (1e-14). `src/tangent.rs` : le
jacobien de Newton assemblé ANALYTIQUEMENT — M, le terme gyroscopique et sa
tangente en (θ, ω) par duaux sur chaque corps, G analytique, et **K = ∂(Gᵀλ)/∂δ
par duaux du premier ordre sur le gradient analytique** (12 ou 18 tangentes),
composés avec la tangente gauche de l'exponentielle J_l(θ) parce que
R = exp(θ)R₀. Contrôles : G par duaux = G analytique (1e-12) ; K = différences
finies centrées du gradient (1e-7) sur liaison, bielle, engrenage ; trajectoires
AD = DF à 1e-14.

**Deux choses apprises en le faisant.**
- La Hessienne de λ·Φ en coordonnées exponentielles (duaux emboîtés) n'est PAS
  la dérivée du champ de gradient sous perturbation gauche — elle porte en plus
  la courbure de exp. Newton exige la seconde ; les deux ont été confondues une
  heure, et la version emboîtée était en plus 2,7× plus lente (144 tangentes).
- **Newton simplifié** : le jacobien se calcule à la première itération du pas
  et se garde tant que le résidu chute de 10× par itération — mesuré : 2 à 2,3
  itérations par pas, UN jacobien par pas, convergence quadratique conservée
  (résidu 5e-4 → 1,5e-10 → 6e-15).

**Redondance structurelle** : QR à pivotage de colonnes de Gᵀ au départ, les
lignes dépendantes sont neutralisées (λ = 0). Le quatre-barres plan en 3D
(3 redondantes) rend le même mouvement que sa modélisation non redondante à
1e-14 m, sans plus jamais retomber sur la SVD.

| banc (inconnues) | v0.1.1 DF, ms/pas | v0.2.0 AD, ms/pas | gain |
|---|---|---|---|
| pendule (9) | 0,24 | 0,009 | ×27 |
| toupie (9) | 0,24 | 0,007 | ×34 |
| double pendule (18) | 0,36 | 0,018 | ×20 |
| bielle-manivelle (35) | 1,03 | 0,044 | ×23 |
| quatre-barres redondant (38) | 0,9 | 0,04 | ×22 |
| **tête S2 de FRELON (143)** | **12** | **0,25** | **×48** |

Erreurs et ordres inchangés (`docs/bancs/0.2.0.json` contre `0.1.1.json`).
**La cible S2 ≤ 0,3 ms/pas est atteinte** — attribué poste par poste
(`chronos()`) : jacobien 0,075, factorisation + résolution 0,136, résidus
0,033 ms/pas. Deux leçons de mesure : la résolution prenait 81 % du pas quand
le LU de nalgebra était refactorisé à chaque itération — la factorisation
(faer, pivotage partiel) se fait UNE fois par jacobien, gardée avec lui ; et
faer parallélise ses factorisations sur rayon par défaut — sur une matrice de
143 la synchronisation coûte 30× le calcul (1,5 ms au lieu de 0,05) :
séquentiel, le parallélisme du noyau est ailleurs. `simule(...,
jacobien="df")` garde les différences finies comme contrôle.

## v0.2.0, suite (4 septembre) — creux, résidu par blocs, éléments en parallèle : O(N^1,1)

Le banc de passage à l'échelle (chaîne de N maillons en rotules, 9N
inconnues) a d'abord rendu un exposant **2,3 à 2,7** : le LU dense (O(N³)),
mais aussi un résidu qui construisait M et G DENSES à chaque évaluation (16 ms
par pas à 1 800 inconnues, plus que le jacobien), et — trouvé en attribuant
le temps — une **SVD dense de tout le système à l'initialisation** : O(N³) une
fois, 670 ms par pas amortis à 3 600 inconnues, plus que tout le reste.

Corrigé : résidu en O(n + nnz) par blocs de corps et d'éléments ; jacobien
assemblé en **triplets**, factorisation **dense (faer, pivotage partiel) sous
160 inconnues, creuse (faer, structure symbolique gardée d'un jacobien à
l'autre) au-delà** ; initialisation consistante par la même factorisation ;
blocs d'éléments et de corps calculés **en parallèle** (rayon) au-delà de 48
éléments — en dessous, le dispatch coûtait ×4 sur la bielle, mesuré.

| N maillons (inconnues) | ms/pas | jacobien | factorisation + résolution | résidus |
|---|---|---|---|---|
| 10 (90) | 0,12 | 0,044 | 0,057 | 0,011 |
| 30 (270) | 0,35 | 0,13 | 0,18 | 0,034 |
| 100 (900) | 1,6 | 0,64 | 0,69 | 0,11 |
| 200 (1 800) | 3,5 | 0,97 | 1,4 | 0,22 |
| 400 (3 600) | 14,9* | 2,6 | 4,8 | 0,69 |
| 800 (7 200) | 13,0 | 2,9 | 7,7 | 1,1 |

\\* la détection de redondance (QR dense à pivotage, O(n·m²)) s'amortit encore
sur les pas à 1 200 contraintes ; elle est sautée au-delà de 2 000 (rang plein
supposé — les grands systèmes de ce noyau sont des chaînes et des maillages)
et c'est une QR creuse à pivotage qui la remplacera. **De 900 à 7 200
inconnues : ×8 d'inconnues pour ×8,3 de temps — exposant 1,05**, la cible
O(N^1,2) est tenue ; la tête S2 reste à 0,23 ms/pas, les petits bancs à leur
temps d'avant.

**Accélération forte, mesurée à 7 200 inconnues (chaîne de 800), et ce qu'elle
dit.** Assemblage du jacobien : 9,9 ms à 1 fil → 1,7 ms à 16 fils (×6), 2,7 à
64 (le découpage coûte). Résolution creuse : **6,9 → 7,3 ms, aucune
accélération, même autorisée à paralléliser** — et c'est structurel, pas un
défaut d'implémentation : l'arbre d'élimination d'une CHAÎNE est un chemin, il
n'offre aucun parallélisme à un solveur direct, quel qu'il soit. La cible
« ≥ 20× sur 64 cœurs » du plan est donc mal posée sur ce banc ; elle se
mesurera sur une structure à arbre large (plusieurs mécanismes, un maillage).
Total sur la chaîne : ×1,6 (Amdahl : le solveur pèse 45 % à 1 fil, 55 % à
16). Réglage : `RAYON_NUM_THREADS=16` suffit à cette taille.

## v0.2.1 (4 septembre) — couples à loi, et la tête S2 en dynamique libre = MBDyn

`Noyau.couple` : couple entre deux corps autour d'un axe (angle relatif
déroulé, vitesse relative) — constant, **ressort-amortisseur**, **servo PD
saturé en couple ET en vitesse** (u = tanh((Kp·(cible − θ) − Kd·ω)/Q), τ =
max(u,0)·Q·clamp(1 − ω/ω_nl) + min(u,0)·Q·clamp(1 + ω/ω_nl) — le patron du n0
FRELON), **gouverneur** (τ = sat(Kg·(Ω(t) − ω), ±Q)). Tangentes ∂τ/∂(θ_a, θ_b,
ω_a, ω_b) par duaux (`tanh`, `max`, `min` par branche dans `Scalar`),
assemblées dans le jacobien avec J_l et les coefficients du schéma. Jugé :
torsion +21 ppm de 2π√(J/k), gouverneur à 8e-7 de Ω(1 − e^{−t·Kg/J}), servo
PD : cible atteinte, couple et vitesse bornés.

**Le banc qui compte** : la tête S2 de FRELON en dynamique LIBRE — gouverneur
sur le mât, trois servos PD, gravité, engrenage — contre le `.mov` MBDyn du
même scénario (`mbdyn_s2.regression`), image par image sur 1,6 s :

```
   ψ du mât : 0,002°   ·   pas de pale : 0,0002°   ·   0,22 ms/pas, Newton 1,01/pas
```

Deux intégrateurs indépendants (multistep MBDyn à tolérance 1e-5 ;
α-généralisé de Lie) sur le même modèle : l'écart est celui des solveurs, et
il est au niveau de la tolérance de MBDyn. Piège payé : le `.mov` sur disque
était celui du dernier run (collectif pur), pas du scénario comparé — 6° d'écart
sur les pas, 0,005° sur ψ ; on fait rejouer MBDyn dans la comparaison.

## v0.2.2 (4 septembre) — l'aéro : le vol de la tête S2 sur vinkulum = le vol MBDyn

`Noyau.pale` : théorie des tranches quasi-stationnaire — stations de Gauss le
long de l'envergure, repère de section (avance, normale, envergure) porté par
le corps, incidence α = −atan2(u_P, u_T), Cl/Cd/Cm interpolés sur la polaire
c81 mono-Mach de FRELON, portance ⟂ au vent relatif, traînée le long, moment
de tangage. `Noyau.inflow` : v_i = √(T/2ρA) suivi avec un retard du premier
ordre (état explicite mis à jour au pas accepté — pas de couplage dans le
jacobien). Tangentes ∂(F, M)/∂(θ, v, ω) par duaux ; AD = DF au dernier chiffre.

**Le vol de la tête S2** (6 s : spin-up sous gouverneur, collectif de
stationnaire, cyclique tournant, servos PD, gravité, aéro), **16 s de calcul,
0,27 ms/pas**, contre le vol MBDyn du même deck :

```
                  vinkulum      MBDyn
   régime         3 161,4       3 162,0 tr/min      (−0,02 %)
   couple gouv.   51,74         51,27 mN·m          (+0,9 %)
   cyclique |A−B| 9,02          9,02 °
   T aéro − poids de tête 2,47 N contre réaction du bâti 2,35 N (grandeurs voisines, pas identiques)
```

**Et le contre-solveur a trouvé un défaut du deck MBDyn de FRELON.** Le premier
vol vinkulum, pales au rayon PHYSIQUE (48–175 mm), rendait un couple 19 % sous
MBDyn ; posées comme le deck les posait (`yc − GRIP_Y0` = 9 mm trop loin, une
formule héritée du deck balancier où le nœud était au pied), l'écart tombait à
3 %. Le deck est corrigé, MBDyn rejoué : son T/W passe de 0,89 à **0,76**, son
couple de 63,6 à 51,3 mN·m — les chiffres publiés le 2 sept. portaient 9 mm
d'envergure en trop. Deux codes qui s'accordent sur un modèle FAUX ne valent
rien ; deux codes qui divergent disent où regarder.

Pièges payés : pale à l'arrêt, |V| = 0 rend la direction du vent indéterminée
et √ non dérivable — Newton stagnait, 10 min perdues avant le plancher à 1 mm/s.
Déclaré, pas modélisé : le vent (air properties), Pitt–Peters, le sillage.

## v0.3.0 (4 septembre) — la poutre géométriquement exacte

`Noyau.poutre` : élément à deux nœuds, déformation et courbure CONSTANTES
(un point d'intégration — pas de verrouillage en cisaillement), énergie
U = ½L(γᵀC_Nγ + κᵀC_Mκ) avec γ = R̄ᵀ(r_B − r_A)/L − e₁ et κ = log(R_AᵀR_B)/L,
`log` sur SO(3) en série près de zéro. Les forces nodales sont −∇U **par
duaux** (12 tangentes) : rien n'est dérivé à la main. Les nœuds sont des corps
ordinaires — la poutre se monte dans un mécanisme par les mêmes liaisons.

| banc | mesuré | référence |
|---|---|---|
| console 10 él., 1ᵉʳ mode | 4,1051 Hz | 4,1126 (1,875⁴EI/ρAL⁴) — **−0,18 %** |
| console 20 él., 1ᵉʳ mode | 4,1131 Hz | **+0,01 %** |
| flèche de bout | 11,426 mm | 11,429 (FL³/3EI) — **−0,03 %** |
| enroulement à 180° (M = πEI/L) | κ 3,124 /m, bout à 0,642 m | π/L = 3,1416, 2L/π = 0,6366 — **−0,6 % / +0,8 %** |

L'enroulement complet est le vrai juge : 180° de rotation, aucune
linéarisation ne le passe. Convergence en maillage d'ordre 3–4 sur le mode.

**Quatre défauts trouvés par ce banc, tous dans l'INTÉGRATEUR, pas dans la
poutre** — et c'est ce qu'un élément raide apporte :
- le critère de Newton ignorait le **plancher d'arrondi des forces internes**
  (EA·ε ≈ 1e-9 N pour EA = 7e6) : « ne converge pas » sur un résidu déjà nul ;
- l'**accélération initiale consistante** vaut M/I sur un nœud d'inertie
  minuscule (7e7 rad/s²), et le prédicteur en gardait 5 rad — bornée à 0,3 rad
  sur un pas ;
- le **prédicteur** extrapolait 12 rad de rotation — borné pareil ;
- Newton part maintenant en **pas amorti** (rebroussement ×½, six fois) et sort
  sur **stagnation** du résidu au plancher.
Aucun de ces quatre ne se voyait sur un mécanisme rigide ; tous les bancs
antérieurs sont inchangés.

**Dette soldée le lendemain, et l'accusé était le juge.** La tangente par duaux
emboîtés rendait ~0 là où les différences finies « mesuraient » 5,8e-4 : je
l'ai déclarée fausse. Elle ne l'était pas. Le champ de forces internes porte
un plancher d'arrondi EA·ε ≈ 1e-9 N ; à ε = 1e-7 une composante NULLE de K
sort à 1e-9/2e-7 = 5e-3 de **bruit pur**, et les deux « écarts » observés
avaient les mêmes chiffres significatifs à un facteur 10 près — la signature
d'un quotient de bruit. Le patron emboîté est prouvé isolément
(`ad::tests_emboite`, contre DF d'un cas scalaire sur SO(3)), la DF est passée
à ε = 1e-5, et les deux tangentes coïncident.

> **Une différence finie ne peut pas vérifier une composante nulle** d'un
> champ qui porte un plancher d'arrondi : elle mesure le plancher divisé par
> le pas. Le seuil de comparaison doit porter ce bruit, sinon c'est le juge
> qui accuse.

**Ce qui tourne en production reste les différences finies** — et c'est la
mesure qui l'a décidé, pas la théorie : l'emboîtée est **2,9× plus lente**
(0,69 contre 0,24 ms/pas, 144 tangentes par multiplication contre 24
évaluations de gradient), exactement comme pour les liaisons. L'exacte reste
comme contre-calcul dans le test. Contrôles en place : **modes rigides
annulés**, et **exacte = DF** sur les douze colonnes. Non symétrique : c'est
établi hors équilibre (Simo–Vu Quoc), pas un défaut.

## v0.3.1 (4 septembre) — analyse modale : le noyau sait dire ses fréquences propres

`Noyau.modes(combien)` : linéarisation autour de l'état courant et modes du
système CONTRAINT. K = −∂f/∂q par différences finies du champ de forces TOTAL
(poutres, couples, pales, gravité) **précontrainte comprise** (Gᵀλ suit la
configuration — c'est ce qui raidit un mât chargé) ; M spatiale ; base du
noyau de G par les vecteurs propres de GᵀG ; problème réduit ZᵀKZ v = ω²ZᵀMZ v
résolu par Cholesky de la masse. Rend (fréquence, forme), triées.

**Le contrôle le plus dense du noyau, et le moins cher.** Trois modes d'une
console d'un coup, en 12 ms, contre 1,875⁴ / 4,694⁴ / 7,855⁴ :

```
   4,11   25,75   72,17 Hz      théorie 4,11  25,77  72,17      écart max 0,07 %
```

Deux secondes de simulation temporelle donnaient UN mode à 0,18 % ; 12 ms en
donnent trois à 0,07 %. Piège payé : la SVD **compacte** de nalgebra ne rend
que min(m, n) lignes de Vᵀ — elle ne contient pas le noyau dès que m < n, le
cas général. On passe par les vecteurs propres de GᵀG.

Déclaré : pas de terme gyroscopique ni de raidissement centrifuge en ω (il
faudrait un problème quadratique aux valeurs propres) — valable au repos ou à
rotation lente. K est symétrisée : la tangente d'une poutre géométriquement
exacte ne l'est pas hors équilibre, et un mode propre n'a de sens que pour sa
partie symétrique.

## v0.3.2 (4 septembre) — modes COMPLEXES : l'amortissement et la stabilité

Une fréquence ne dit pas si un mode diverge. `Noyau.modes_complexes` rend
(f, ζ, σ) : valeurs propres de la forme d'état [[0, I], [−M⁻¹K, −M⁻¹C]] sur le
noyau de G, avec **C = −∂f/∂u** par différences finies (Kd des servos, pente
des saturateurs, traînée aéro) **plus le terme gyroscopique** ∂(ω×Jω)/∂ω =
[ω]×J − [Jω]×, qui vit dans le résidu et pas dans f. **σ > 0 = instable.**

```
   oscillateur amorti   f 7,86771 Hz · ζ 0,150000      écart 2e-16 sur les deux
   boucle servo PD      f 109,438 Hz · ζ 0,700000      écart 0,000 %
   amortissement < 0    σ = +2,50                      INSTABILITÉ détectée
```

### L'état tournant : trois corrections justes, une limite structurelle, un refus

La toupie rapide sortait INSTABLE (σ = ∓145) alors qu'elle est stable. Trois
défauts trouvés en creusant, chacun réel et corrigé :

1. **la raideur ignorait ω × J_s(q) ω**, qui vit dans le RÉSIDU et pas dans f.
   Sa dérivée en orientation vaut ∝ ω² : c'est le **raidissement centrifuge**,
   celui qui donne son ν_β ≥ 1 à une pale de rotor ;
2. **λ n'était pas calculé** hors simulation, donc la précontrainte était nulle
   — et le couple de rappel d'un poids ne vient PAS de ∂f/∂q (une force
   constante), il vient de ∂(Gᵀλ)/∂q. Sans λ, la toupie n'avait plus de gravité ;
3. **le bloc (0,0) de la forme d'état n'est pas nul** autour d'un état tournant :
   R = exp(s)R₀(t) donne ṡ = δω + ω₀ × s (et non ṡ = δω). Le signe compte : à
   l'envers, les fréquences doublent.

Avec les trois, la toupie devient **stable — σ ~ 1e-24, purement imaginaire**,
ce qui est le verdict juste. Ses fréquences restent fausses de 47 et 92 %.

**La cause a été MESURÉE, et ma première explication était fausse.** J'avais
écrit — dans le code, ici, et dans un commit — que « le transport sort de
l'espace admissible, donc la projection perd de l'information ». Mesuré
(`fuite_transport`) : la fuite vaut **5,0 %**, bien trop peu pour expliquer
92 %. Ce que dit la mesure (`diag_modal`, toupie à 300 rad/s) :

```
   ‖Kr‖ = 50,85 N·m/rad        la raideur PHYSIQUE cherchée : mgl = 0,049
   ‖T‖  = 423                  fuite du transport : 5,0 %
   contrôle à ω = 0 : ‖Kr‖ = 0,0692 = mgl·√2 — la raideur est EXACTE
```

La raideur tangente est dominée par les termes en ω²·(Jt − Ja) = 18, soit
**mille fois** la raideur cherchée. Le résultat physique n'apparaît que par
**annulation quasi exacte** de ces grands termes avec ceux du changement de
variables : une soustraction catastrophique. 5 % d'erreur sur le transport,
c'est 2,5 N·m/rad — **cinquante fois mgl**. Ce n'est pas la fuite qui est
grande, c'est le terme sur lequel elle porte.

> Le correctif n'est donc pas « mieux projeter » — ce que ma fausse
> explication aurait fait chercher. C'est une formulation où ces grands termes
> n'ont pas à se soustraire : repère tournant, coordonnées relatives, Floquet.

**Et la leçon de méthode, qui vaut plus que le correctif** : une explication
de défaut ne se publie pas sans mesure. Les deux diagnostics
(`fuite_transport`, `diag_modal`) font quinze lignes chacun, et le premier a
réfuté mon explication en une exécution. Un contrôle qui ne peut que confirmer
ne vaut rien.

> **`modes_complexes` REFUSE donc de répondre dès qu'un corps CONTRAINT
> tourne**, en disant pourquoi. Un noyau qui rendrait ces chiffres serait pire
> qu'un noyau qui refuse.

Les modes NON tournants sont exacts — c'est ce que les trois contrôles
prouvent, et rien de plus.

## L'ordre des corps flexibles, corrigé (3 septembre, soir)

**Vérifié le 3 septembre au soir, et corrigé.** J'avais écrit « corps
flexibles par réduction modale (Craig–Bampton) AVANT les poutres
géométriquement exactes ». C'est l'ordre INVERSE pour une voilure tournante :
les codes de référence du métier (CAMRAD II, RCAS, Dymore) modélisent la pale
en **poutre non linéaire / géométriquement exacte**, parce que la réduction
modale linéaire manque la raideur géométrique (centrifuge) et les grandes
rotations — il lui faut des « dérivées modales » pour s'en approcher, et
encore. Ordre retenu :

1. **poutres géométriquement exactes** (Simo–Reissner sur SO(3), le même
   groupe que l'intégrateur) — pales, mât ;
2. **Craig–Bampton** pour ce qui est raide et linéaire — cellule, berceau —
   avec dérivées modales si l'interface bouge beaucoup.

## Revue de complexité du cœur (4 septembre) — ce qui a été coupé, et ce qui a résisté

Une revue « qu'est-ce qui peut disparaître » sur les fondements (SO(3),
α-généralisé, Newton, contraintes, projections). Les fondements tiennent ;
ce qui est parti est du doublon et de l'option morte, jamais de la physique.
`verification` (166 s), `campagne` 12/12, `cargo test`, clippy et la doc
générée sont verts sur le résultat.

**Coupé :**
- `expm` et `log_so3` de `lib.rs` étaient des doublons de `ad::expm` et de
  `tangent::logm` génériques — une seule Rodrigues, un seul log, appelés à
  `T = f64` ;
- la bissection à la main d'`interp` (pales) → `partition_point` ;
- **`simule(projection=True)`** — la projection de vitesse PAR PAS, refusée
  par défaut depuis sa mesure (Φ̇ annulé, ordre 2,00 → ~0,8, énergie dissipée
  ×100 sur 100 périodes), n'était passée par aucun appelant. Retirée du
  binding et de la boucle ; `projette_vitesse` SURVIT parce qu'`assemble`
  l'appelle à l'initialisation, où elle est licite. Le banc
  `projection_vitesse` devient `violation_vitesse` : il garde ce qui reste
  vrai sans l'option — Φ̇ à O(h²), borné, sans coût d'énergie (9e-6 sur 100
  périodes) — et porte la mesure de l'option en TRACE dans son docstring ;
- **`simule(jacobien="df")`** — le jacobien par différences finies, jamais
  passé par un appelant (« contrôle » de l'AD que personne ne lançait) ;
- `Pas::etat` et `Pas::dtheta` recopiaient a₁ et dq : `a1_dq` une fois ;
- le hachage FNV maison du motif creux → `std::hash::DefaultHasher` ;
- la condition de rebroussement à quatre termes → `k_bt >= 6 || (fini &&
  ‖r‖ ≤ ‖r₀‖)`, même logique, lisible ;
- la permutation du QR pivoté reconstruite via une matrice 1×m et une ligne
  « no-op de typage » → `from_fn` et une boucle.

**A résisté, et c'est mesuré :** la sortie « stagnation » de Newton (résidu
qui ne descend plus sous 1e-6·(1 + ‖f‖)). Retirée au motif que le critère
sur l'incrément (`dx_rel < 1e-11`) la couvrirait, elle a fait casser
`verification` au premier passage : « Newton ne converge pas à t=2.36332
(résidu 1.14e-7) ». Les deux critères ne se recouvrent pas — remise, avec la
mesure écrite à côté. **Retiré aussi de la liste, avant d'être touché :**
`energie=True`, que la revue avait pris pour une option morte ; `bancs`
l'exerce sur un corps libre (ordre 1,9+, erreur à 5 % de celle sans
projection) et vérifie ses deux refus. Une revue de complexité se juge
comme le solveur : sur la table, pas sur l'intuition.

Net : **−110 lignes** sur `lib.rs`, `tangent.rs` et `bancs.py`, API réduite
de deux arguments que personne ne passait.

## Trois fragilités du cœur, fermées ou mesurées (4 septembre, suite)

La revue de complexité avait aussi rendu un jugement de fond : formulation
saine, mais **pas de stabilisation de contrainte**, **masque de redondance
figé à t = 0**, et **sept seuils de Newton posés à la main et mesurés sur le
seul corpus de l'auteur**. Les trois sont traités ; ce qui ne l'est pas est
dit en fin de bloc.

**GGL (`simule(ggl=True)`).** Gear–Gupta–Leimkuhler sur l'α-généralisé de Lie
(Arnold–Brüls 2007) : m multiplicateurs ζ de plus, la pose reçoit βh·G₀ᵀζ
(G au début du pas ; le facteur βh rend ∂Φ/∂ζ d'ordre 1, comme ∂Φ/∂u̇), et la
ligne G·u₁ + Φ_t = 0 entre au système, mise à l'échelle 1/(γh) comme le non
holonome. Φ_t est pris à q₀ par différence finie — exact pour une cible à
taux constant, O(h) sinon, déclaré. Les lignes inactives ou déjà en vitesse
portent ζ = 0.

```
   pendule 30°, 2 périodes        index 3                 GGL
   erreur T/100 → T/400           2,14e-3 → 1,44e-4       1,92e-3 → 1,30e-4
   ordres                         1,97 · 1,92             1,97 · 1,91
   |Φ̇| à T/400                    2,6e-5                  3,3e-16
   itérations de Newton           1 590                   2 128   (×1,3)
```

L'ordre est gardé, Φ̇ tombe au plancher, l'erreur en position est même un
peu plus petite. La projection de vitesse après coup, retirée le matin,
faisait tomber l'ordre à 0,8 : la différence est que ζ est dans la
FORMULATION du pas, Newton le résout avec u̇ et λ, et l'accélération
algorithmique reste cohérente.

*Le jacobien, et ce qui y manque.* Les colonnes ζ sont EXACTES par
différences finies sur le résidu — la dépendance en ζ passe par q, donc par
tout ce que le résidu contient (Gᵀλ, forces, gyroscopique), et une copie du
modèle par colonne les capture toutes ; O(m·(n + nnz)) par jacobien, c'est le
plafond déclaré (`ponytail:` dans le code) et la raison pour laquelle GGL
reste **opt-in** : l'index 3 direct tient l'échelle, GGL sert là où la
cohérence des vitesses compte (réactions, états initiaux de Floquet). Les
lignes G·u, elles, n'ont que G·c en colonnes u̇ : le terme (∂G/∂q·u₁)·βh² est
laissé de côté, O(h·|u|/L) relatif — Newton converge linéairement à ce taux,
d'où les 1,3× d'itérations. C'est mesuré, borné par un assert (< 3×).

**Le masque de redondance n'est plus figé.** Banc : le parallélogramme parti
À PLAT — A, B, D, C alignés, point de changement où les branches
parallélogramme et antiparallélogramme se croisent et où G perd un rang de
plus qu'ailleurs. Détecté là, le QR neutralisait quatre lignes au lieu de
trois ; la quatrième redevient structurelle dès que le mécanisme sort de
l'alignement, et le noyau rendait « Newton ne converge pas à t=0.011 ».
Désormais : quand ‖Φ‖ dépasse `tol_phi` OU quand Newton cale, on re-détecte
à la pose courante ; si le masque change, l'état d'avant le pas est ramené
sur Φ = Φ̇ = 0 sous le nouveau masque (`assemble` + projection de vitesse —
la seule place où elle est licite) et le pas se rejoue. Référence fermée,
sans gravité : la branche parallélogramme tourne uniformément (énergie
cinétique indépendante de l'angle), centre du coupleur (1 + cos φ, sin φ).

```
   φ₀ = 5° (hors bifurcation)   index 3 8,5e-8 m   GGL 3,7e-8 m   3 redondantes
   φ₀ = 0° (à plat)             index 3 1,8e-6 m   GGL 1,9e-6 m   3 redondantes
```

Le départ à plat paie l'assemblage (×20 sur l'écart, à 1e-6), et c'est
asserté dans ce sens : si un jour il ne le payait plus, le mécanisme de
rejeu aurait changé de nature.

*Deux essais avant celui-là.* Rejouer le pas depuis l'état dérivé sans le
ramener sur Φ = 0 : Newton devait résorber Φ/(βh²) d'un coup et divergeait.
Ne re-détecter que sur la tolérance de Φ : en GGL, la dérive de la ligne
inactive est plus lente et Newton calait AVANT que Φ ne la trahisse
(r_dyn 5e-9 sans descendre, 24 itérations) — d'où le second déclencheur.

**Les seuils sont mesurés, pas seulement posés.** `VINKULUM_SEUILS=k`
multiplie d'un coup la borne du prédicteur (0,3 rad, 5 % de longueur), le
rafraîchissement du jacobien (0,1), le plancher de stagnation (1e-6) et le
critère d'incrément (1e-11). `bancs.seuils()` rejoue la maquette, l'ordre du
pendule et la bifurcation à ×0,5, ×1, ×2 : **les trois passent**. Ce n'est
pas une preuve que les seuils sont bons ; c'est la preuve que ces trois bancs
ne les mesurent pas — et le jour où l'un cassera à ×2, on saura qu'il
mesurait un seuil.

**Ce qui reste, dit tel quel :** contact non lisse (Moreau–Jean, LCP) — une
semaine, pas un correctif ; décrochage dynamique et sillage — feuille de
route physique, pas fragilité du solveur ; SE(3) au lieu de SO(3)×R³ —
l'invariance est mesurée à 1e-14, pas garantie par construction ; estimateur
d'erreur sur λ — le demi-pas ne regarde que la dynamique ; validation contre
une EXPÉRIENCE — il n'y a pas de données ici, et on n'en fabrique pas.

## Contact non lisse, décrochage dynamique, Φ à mi-pas — et pourquoi pas SE(3) (4 septembre, soir)

Demande : contact non lisse, décrochage dynamique + sillage, SE(3), estimateur
d'erreur sur λ. Trois sur quatre sont là ; le quatrième est refusé pour une
raison mesurable, écrite en fin de bloc.

**Contact non lisse** (`contact(..., nonlisse=True, restitution=e)`). La
normale cesse d'être une pénalité : c'est un multiplicateur λ ≥ 0 sous
complémentarité, porté par un élément `Elem::K` — une ligne de plus dans λ,
écrite avec Fischer–Burmeister φ(a, b) = a + b − √(a² + b²), semi-lisse, que
le Newton existant résout (la ligne reçoit un facteur ∂φ/∂a sur G·c et ∂φ/∂b
sur sa diagonale). Le frottement reste le Coulomb lissé, appliqué sur λ.
Pour les analyses (assemblage, redondance, modes, statique) le contact est
ABSENT : Φ et G nuls, la ligne jamais masquée.

*Trois versions avant celle qui tient, chacune mesurée.*

1. **Ligne de POSITION** (gap ≥ 0 ⊥ λ). Pente à 0,001 % — et la bille
   lâchée de 10 cm **remontait à 25 cm**. Une contrainte de position sous
   Newmark rebondit d'elle-même : q₁ bloqué ⇒ a₁ impulsif ⇒
   u₁ ≈ u₀(1 − γ/β) ≈ −u₀. Restitution parasite ≈ 1, c'est le résultat
   classique qui a fait naître Moreau–Jean.
2. **Ligne de VITESSE** (ġap₁ + e·ġap₀ + gap₀/h ⊥ λ, contact actif si le gap
   prédit ferme). Encore un rebond à e = 0, à 0,8 fois la vitesse
   d'arrivée : l'α-généralisé GARDE l'accélération impulsive du pas d'impact
   et la réinjecte au suivant par (1 − γ)h·a. Et le rappel gap₀/h « en un
   pas » dictait la vitesse de sortie à la place de la loi d'impact.
3. **Retenu** : loi de Newton pure ġap₁ = −e·ġap₀, rappel de pénétration
   doux (κ = 0,2, l'ERP d'ODE, ne s'applique qu'à gap < 0, converge en κβ/γ
   par pas sans dépasser), et **redémarrage du schéma** après tout pas à
   contact actif — l'accélération lisse, consistante avec les contacts
   actifs, remplace la mémoire impulsive (Chen–Acary–Virlez–Brüls 2013
   séparent de même les parts lisse et impulsive). Deux défauts de plus sont
   sortis en mesurant : l'activation par le seul gap prédit BATTAIT à
   gap = 0⁺ (inactif → chute de gh²/2 → impulsion ×3 ; λ à −3,5 % au repos,
   +7,8 % sur l'accélération de la pente) — un contact déjà fermé reste
   actif, λ ≥ 0 le laisse partir ; et `forces()` ignorait le Coulomb sur λ,
   donc l'accélération lisse du redémarrage n'avait pas de frottement
   (+7,9 % sur la pente, exactement (1 − γ)·μg cos θ).

```
   bille lâchée de H = 0,1 m, h = 1e-4
     impact à 0,1427 s (√(2H/g) = 0,1428) · pénétration max < 1e-6 m
     e = 0   : aucun rebond, λ final = mg à 1e-5
     e = 0,5 : sommets 2,509e-2 et 6,30e-3 pour e²H = 2,5e-2 et e⁴H = 6,25e-3
   pente 20°, Coulomb sur λ : glissement +0,001 %, roulement +0,001 %
```

**Ce que ce contact n'est pas** : le pas de contact est d'ordre 1 (le
redémarrage) ; le frottement est lissé (tanh), pas un cône du second ordre ;
la raideur de contact λ·∂G/∂q est hors jacobien (quasi-Newton, comme la
pénalité qui ignore déjà sa raideur en position) ; ni IPC ni CCD ne
s'appliquent à la ligne non lisse.

**Décrochage dynamique d'Øye** (`pale(..., oye=True)`). Un état f par
station : CL = f·CL_att + (1 − f)·CL_fs, CL_att = a₀(α − α₀) par moindres
carrés sur ±5° autour du zéro de portance, CL_fs la table « entièrement
décollée » déduite de la polaire statique par f_st = (2√(CL/CL_att) − 1)²,
et ḟ = (f_st(α) − f)/τ avec τ = 4c/|V|, exponentielle exacte sur le pas
(patron des états de Wagner). Le c81 statique est la limite t → ∞ par
construction. Échelon 5° → 20° sur une section fixe, polaire synthétique :
portance ×2,185 à 0⁺, 1,0001× la statique à 0,2 s, **τ mesuré 20,00 ms =
4c/V**. Pas de Leishman–Beddoes (retard de pression, tourbillon de bord
d'attaque) — déclaré. **Sillage** : rien de plus que Pitt–Peters, qui est
déjà le sillage dynamique à trois états ; un sillage libre est un projet,
pas un correctif — non fait, dit.

**Estimateur de demi-pas : Φ aussi.** Le résidu à t + h/2 ne regardait que
les lignes dynamiques ; il regarde désormais |Φ(q_{h/2})|, relatif à la
taille du modèle — l'erreur locale que λ paie, invisible au résidu
d'équilibre — et rejoue le pas si elle dépasse la tolérance. Sur le pendule
adaptatif : 2,7e-8, 15 rejeux ; nulle sur un modèle sans liaison, ce que le
banc vérifie aussi. `adapt_stats()` la publie en cinquième.

**SE(3) — refusé, et voici pourquoi.** Les inconnues du noyau sont les
vitesses SPATIALES du centre de masse. La mise à jour SE(3) (r₁ = r₀ +
h·V(hω)·v, V le jacobien de exp) intègre exactement un mouvement HÉLICOÏDAL
à torseur constant dans le corps ; pour un corps libre, dont le centre de
masse va DROIT, elle courberait la trajectoire — c'est faux, pas plus
précis. Le gain de Sonneville–Brüls vient d'une formulation entière en
torseurs matériels (unknowns, inertie, forces, toutes les tangentes), soit
un autre noyau. Ce qu'on en attendrait — l'invariance par changement de
repère — est déjà mesuré à 1e-14 par `verification`. On ne réécrit pas un
solveur pour une propriété qu'on a.

> **Complété le 5 septembre : ce refus ne répondait qu'à MOITIÉ.** Il réfute
> la mise à jour hélicoïdale, ce qui est juste. Il ne répond pas à l'autre
> argument, celui que la littérature met en avant (note d'architecture de
> sept. 2026, §1.2) : **en vitesses de repère CORPS, la matrice de masse ne
> dépend plus de la configuration**, ce qui simplifie le jacobien et améliore
> le conditionnement. C'est vrai, et c'est le prix que ce noyau paie —
> `J_s = R·J·Rᵀ` se recompose partout, et le gyroscopique `ω × J_s ω` vit dans
> le résidu avec sa tangente exacte `[ω]×J_s − [J_sω]×`.
>
> La réponse tient en deux lignes, et la seconde est neuve. **Un :** avoir M
> constante exige que les inconnues, l'inertie, les efforts ET toutes les
> tangentes soient en torseurs matériels — c'est-à-dire un autre noyau, pas
> une option. **Deux :** le conditionnement n'est pas le poste dur, et ce
> n'est plus une opinion — `refus_gardes` assert **zéro repli SVD** (le LU
> tient à chaque pas ; quand il ne tient pas, `Modele::resout` bascule en
> moindres carrés et `stats()` le compte). On ne réécrit pas un solveur pour
> améliorer un conditionnement dont rien ne se plaint.

## Les quick wins, du moins cher au plus cher — et le vol long a cassé à 50,8 s (5 septembre)

Six candidats classés par coût ; cinq faits, un bloqué. Chacun a rendu
quelque chose qu'on ne savait pas.

**Seuils de Newton sur la MACHINE** (`vinkulum_s2.py seuils`, FRELON).
`bancs.seuils` ne rejouait que deux jouets ; le vol S2 sous
`VINKULUM_SEUILS` ×0,5 / ×1 / ×2 rend régime 3160,83 tr/min et couple
53,312 mN·m **identiques au 5ᵉ chiffre**, Newton 2,472 / 2,471 / 2,468 par
pas. Seuls les jacobiens bougent (1,18 / 1,07 / 1,03 par pas) — les seuils
mordent, la physique non.

**Contacts multiples qui basculent** (`contact.multiples`). Pile de trois
billes non lisses, e = 0 : λ = 3mg / 2mg / mg à 1e-5, jeux < 1e-6 m. Trois
billes alignées, e = 1 : v → (0, 0, 0,5000), p et E à 1e-4. Chaque impact
active une ligne K, en relâche une autre, redémarre le schéma — ce que la
bille seule ne testait pas.

**La garde du bassin** (`Liaison::bassin`, `Modele::verifie_bassin`). La
branche θ = π du résidu vee(skew E) était « documentée, pas gardée ». Mesuré
avant de garder : une cible de rotation qui saute de **0,3 rad en UN pas**
(h = 0,01) y tombe en silence — l'encastrement accepte, |Φ| 1e-16, le corps
est retourné de π − 0,3. Même saut réparti sur trois pas ou plus : juste.
Le seuil coïncide avec la borne du prédicteur (0,3 rad). Garde : sur un pivot
(deux lignes bloquées) E_kk sur l'axe LIBRE, insensible à la rotation
propre ; sur un encastrement tr E ; < 0 = refus, à l'acceptation du pas et
à l'assemblage. Une seule ligne bloquée (cardan) : indécidable, non gardé.

**Φ_t analytique, à la pose COURANTE** (`Liaison::phi_dt`). Cherché comme
raffinement, trouvé comme défaut : sur une rotation imposée à 8 rad/s, GGL
rendait |Φ̇| **2,6e-4 contre 4,3e-6 en index 3** — 60× pire — parce que Φ_t
était pris par différence finie à la pose du début du pas alors que G·u est
pris à la pose de Newton. Analytique (−dir·θ̇ ; ∂E/∂t = −E·[axe]×·θ̇), à la
pose courante : **1,3e-9**. Et `Loi::derivee` rend la pente des tables.

**Vol long — LE résultat.** Aucun banc ne dépassait 3 s hors le corps libre.
`vinkulum_s2.py long 60` : **Newton ne converge pas à t = 50,84 s**, résidu
5,25e-3 ; à h = 2e-5 dès 50,4 s : 50,83 s, résidu 1,31e-1 ; sur un vol de
100 s (phase de cyclique différente) : 50,89 s. Rien ne dérive avant — |Φ|
1e-16, Newton 3,5/pas, réactions et couples sains à 50,82. Le seul chiffre
qui ne cesse de croître est l'angle DÉROULÉ du pignon : −4·16 388 =
**−65 552 rad**. À 65 536 rad l'ulp du double passe de 7,3e-12 à 1,5e-11 ;
la ligne d'engrenage θ_a − 4θ_b est mise à l'échelle 1/(βh²) = 4e8 à
h = 1e-4, soit 6e-3 de résidu incompressible — c'est le 5,25e-3 mesuré, et
×25 à h = 2e-5 (1e10 × 1,5e-11 = 0,15 ≈ 0,131). La prédiction du temps de
casse : 65 536/(4 × 331,1) + 1,3 s de spin-up = 50,8 s. **Correction** :
Φ = θ_a − r·θ_b est invariant sous (θ_a − 2πk·r, θ_b − 2πk) — on ramène θ_b
dans un tour à chaque pas accepté, θ_a suit. La vis–écrou garde son angle
déroulé (sa translation porte la même origine) — même mur à ~10 000 tours,
déclaré.
Rejoué : **60 s de vol tenus**, 600 000 pas, 0,300 ms/pas, Newton 3,02/pas,
régime par fenêtre de 5 s entre 3161,39 et 3161,43 tr/min — dérive
**0,001 %** — Q 52,72 mN·m. Suite complète verte derrière (verification,
bancs, contact, campagne).

**CI** : `ci/github-actions.yml` existe ; le jeton de ce dépôt n'a pas le
scope `workflow`, elle ne peut pas être poussée d'ici. Rien à coder.

**Trouvé en route, fermé le lendemain** — voir la section suivante.

Deux corrections de dépôt au passage : clippy 1.98 refuse `impl` après le
module de tests (déplacé en fin de fichier) ; `nproc` = 1 sur cette
machine, deux simulations en parallèle se partagent le cœur.


## L'ÉCHEC GGL EST INSTRUIT — c'était l'incrément d'une différence finie (5 septembre)

GGL passait partout sauf sur `srscm`, `six_barres` et la **tête S2 de FRELON**.
Instruit sur ordre de Paul : mesurer d'abord, juger ensuite la pertinence de la
forme moderne complète que décrit une revue externe (projection mise à
l'échelle par la masse, élimination de μ dans le pas de Newton, ∂(Gv)/∂q par
différentiation automatique).

### Ce qui NE sépare pas les deux groupes

Relevé modèle par modèle avant de toucher au code : ni `h` (3,93e-4 passe,
2e-4 échoue ; 8e-3 échoue, 2,15e-2 passe), ni le type d'élément (l'engrenage
seul passe, une cible de rotation à 20 rad/s passe, le cardan seul passe), ni
la planéité, ni le nombre de lignes, ni la durée simulée. Cinq pistes écartées
au tableau, sans un essai de correction.

Le seul séparateur : les trois rouges ont **au moins deux boucles fermées
indépendantes** et un **G de rang plein**, quand les modèles à boucle qui
passent portent des lignes masquées.

### L'ÉCHELLE DE BOUCLES — reproduire hors de FRELON (`bancs.ggl_boucles`)

Chaîne de parallélogrammes, 1/2/3 boucles × pivot (plan, redondant) ou rotule
(rang plein) × deux **échelles séparées d'un facteur 100**. Le temps est
adimensionné (h et t_end suivent τ = √(L/g)), donc deux échelles ne diffèrent
QUE par les unités. **Six des douze cas cassaient** — dont, à L/100, une seule
boucle. L'échec était reproduit dans un modèle de trente lignes, sans FRELON.

### LA CAUSE

Les colonnes ζ du jacobien sont des différences finies, et leur incrément
était `dxj = 1e-6·max(|ζ|,1)` — calé sur **ζ**, alors que ζ n'agit sur le
résidu que par la pose, qu'il déplace de **βh²·(ligne j de G₀)**. À
h = 1e-3, βh² ≈ 3e-7 : la pose bougeait de 3e-13 en relatif, soit **trois
chiffres significatifs** avant division par l'incrément. La colonne était du
bruit, et plus le système couple (boucles) ou se raidit en conditionnement
(petite échelle), moins Newton le tolère.

Calé sur la sensibilité — `dxj = 1e-8·L/(βh²·‖G₀ⱼ‖)`, avec un repli pour les
lignes masquées, dont G₀ ne porte aucun terme :

```
   échelle de boucles     6 échecs / 12   →   1 / 12
   six_barres             ÉCHEC t=1,23 s  →   OK
   tête S2 (regression)   ÉCHEC 2e pas    →   OK, pas de pale à 3,5e-6°, |Φ| 1,8e-15
   tête S2 (dynamique)    ÉCHEC 2e pas    →   OK
```

*Quatre chiffres significatifs sur une colonne de jacobien ne se voient nulle
part : le résidu ne ment pas, il stagne.*

### CE QUI RESTE, ET IL EST MESURÉ AUSSI

`srscm` casse toujours, au même instant qu'avant correction. Ce n'est ni le
cardan (retiré : même casse, à la milliseconde et au chiffre près), ni la
vitesse (gravité ÷10 : passe ; ÷100 : casse à 2,38 s — non monotone, donc pas
un effet de vitesse), ni le masque de redondance (0 ligne masquée, 0 rejeu de
bout en bout). C'est la **configuration** :

```
   t          2,600     2,700     2,720     2,723     2,724     2,730     2,800
   cond(GGᵀ)  5,9e2     1,2e3     2,7e3     5,6e4    2,8e6     1,3e3     4,0e2
```

Un **pic de conditionnement d'un facteur 4 800 à l'instant exact de la
casse**, rang inchangé à 17 : le mécanisme spatial passe par un point mort.
L'index 3 le traverse ; GGL, qui ajoute m inconnues couplées par ce même G,
non. Même motif pour la dernière cellule de l'échelle (L/100, 3 boucles,
redondante, cond 2,7e5).

> **Portée exacte de GGL, à partir d'aujourd'hui** : stabilisation d'index 2
> qui tient tant que le conditionnement de G reste modéré ; elle ne traverse
> pas une configuration quasi singulière. Le contre-exemple est nommé et gardé
> par un banc, il n'est pas caché derrière un ε de régularisation — qui aurait
> été un paramètre de réglage de plus, et la revue elle-même le dit.

### LA MÉTRIQUE DE MASSE : PRÉDICTION RÉFUTÉE, MESURÉE

La revue attribue le mal à une projection dimensionnellement incohérente
(q̇ = v − Gᵀμ mêle mètres et radians) et propose q̇ = v − M⁻¹Gᵀμ. **La moitié
du diagnostic est confirmée** : `cond(GGᵀ)` suit **1/L²** — ×1e4 quand la
machine rétrécit de 100, asserté sur les six paires de l'échelle. C'est bien
un conditionnement qui dépend du choix d'unités.

**Le remède proposé, lui, ne marche pas ici** : `cond(G·M⁻¹·Gᵀ)` est **PIRE
que `cond(G·Gᵀ)` sur les douze cellules**, d'un facteur 5 à 20 — parce que
pour une barre élancée M⁻¹ porte lui-même un rapport 12/L² entre rotation et
translation, et l'ajoute au lieu de le retirer. Le banc **asserte le résultat
négatif** : si M⁻¹ améliorait un jour, la conclusion devrait être rejouée.

Restent à juger, non faits : l'élimination de μ par un solve GM⁻¹Gᵀ (qui
supprimerait les m clones du modèle par jacobien) et ∂(Gv)/∂q par AD (le terme
que `tangent.rs` omet). Ni l'un ni l'autre n'est nécessaire pour ce qui vient
d'être réparé. Déjà en place et confirmées par la revue : α-généralisé sur
groupe de Lie, redondance par rang révélé, contact hors GGL au niveau vitesse,
initialisation cohérente.

## UNE NOTE D'ARCHITECTURE CONFRONTÉE AU NOYAU — et cinq refus qui n'avaient pas de garde (5 septembre)

Une note de synthèse (« Avancées mathématiques pour la prochaine génération de
solveurs multicorps industriels », sept. 2026 : groupes de Lie, mécanique
discrète, contact convexe + ADMM, simulation différentiable, réduction
structure-preserving, algèbre linéaire) a été lue contre le code.

Sur une vingtaine de recommandations : **neuf sont déjà là**, quatre étaient
refusées par des mesures antérieures, une est **fausse ici**, une est
**contre-indiquée**, et une seule mérite d'être lue sans être écrite.

### Le chantier n'est pas venu des recommandations, il est venu de l'état des refus

Ce dépôt a un idiome pour garder un résultat négatif : l'**assert inversé**,
qui casse quand la raison d'un refus cesse d'être vraie — `bancs.adaptatif`
(« si un jour cet assert casse, le verdict devra être rejoué »), la métrique
de masse réfutée, ρ∞ haut, l'exposant d'échelle, `campagne`. **Cinq refus n'en
avaient aucun**, et un refus mesuré une fois puis jamais rejoué n'est plus une
mesure : c'est une croyance datée. `verification.refus_gardes` les porte.

```
   LDLᵀ              asym(C)/‖C‖ = 1,500   le gyroscopique n'est pas symétrique
   précision mixte   tol 1e-12 < eps(f32) 5,96e-8   —   6e4× trop grossier
   Schur             replis SVD = 0        le LU tient à chaque pas
   Cayley            t_exp/t_jac = 0,4–1,0 %  (chrono POSÉ pour ça)
   SE(3)             `invariances` présent (le refus cite ce contrôle)
```

Les cinq contrôles négatifs sont faits : chacun **mord** quand on lui retire
sa prémisse (C symétrisée, tolérance desserrée, repli SVD forcé, `t_exp` porté
à `t_jac`, `invariances` supprimée). Sans eux, on aurait cinq asserts dont on
ne saurait pas s'ils peuvent casser.

### DEUX AFFIRMATIONS DÉJÀ ÉCRITES SONT TOMBÉES, ET C'EST LE VRAI ACQUIS

**Le refus de Cayley n'était pas mesuré.** Je l'avais refusé en citant
« 1,07 jacobien par pas » — chiffre juste, mais qui décrit le **vol S2 de
FRELON**, pas le noyau : sur des chaînes tournantes il vaut **1,44 · 2,10 ·
4,88** à 18 · 180 · 540 inconnues, et le jacobien pèse **26 à 42 %** du pas.
Le bon chiffre n'existait pas. Il a donc été posé — `chronos()` rend
maintenant `t_exp` en sixième, la part du jacobien passée à construire les
`J_l` — et il vaut **0,4 à 1,0 %**, décroissant avec la taille. La raison est
structurelle et se lit dans le code : `jac_gauche` est appelée une fois par
CORPS et par jacobien, **pas par élément ni par itération de Newton**, alors
que l'argument de la note (§1.3, « le coût par itération de Newton est
réduit ») suppose l'inverse. La prémisse ne s'applique pas à cette
implémentation.

**« Le LU n'est PAS le point dur » est un argument de MÉMOIRE**, publié depuis
le 3 septembre comme un verdict général. Le remplissage est bon (6–7 non-nuls
par ligne avec COLAMD) et il en découle quelques mégaoctets à 30 000 inconnues
— rien de plus. Le **temps** dit autre chose : la résolution pèse **32 / 65 /
52 %** du pas aux trois tailles. Les deux affirmations sont vraies et ne
parlent pas de la même chose ; les confondre ferait écarter une
sous-structuration au motif d'un remplissage. Ce qui écarte Schur est
ailleurs, et c'est maintenant asserté : zéro repli SVD.

*Une mesure qui décrit une machine et sert à en juger une autre : c'est la
famille de défauts que ce dépôt a payée le plus souvent, et elle vient de se
reproduire sur mon propre verdict, une heure après l'avoir écrit.*

### Trois autres traces périmées, trouvées en croisant la note et le code

- la puce « **GGL** non fait, par mesure » du README contredisait la table du
  même fichier depuis le 4 septembre. Ce qui reste vrai — le défaut que GGL
  répare ne fait pas de mal — est gardé ; le reste est daté ;
- la ligne « sensibilité » du README fondait deux choses qui ne se ressemblent
  pas : l'adjoint sur **ÉQUILIBRE** est branché au noyau (`k_c_m_z`,
  `d_residu_poutre`), l'adjoint **EN TEMPS** tourne sur un modèle NumPy
  **séparé** (`M ü + g`, pas de Φ, pas de λ, pas de contact). Le code le
  déclarait, le README non ;
- `sensibilite.equilibre` annonçait « le noyau n'a pas de solveur statique ».
  `Noyau.statique` existe. La phrase envoyait le lecteur écrire ce qui est
  déjà là.

### Ce qu'on n'implémente pas, et pourquoi

| recommandation | motif |
|---|---|
| **SE(3)** (§1) | inconnues = vitesses spatiales ; M constante exige une formulation entière en torseurs matériels = un autre noyau. Et le conditionnement n'est pas le poste dur : 0 repli SVD, asserté |
| **Cayley** (§1.3) | `jac_gauche` est hors de la boucle interne — mesuré, 0,4–1,0 % du jacobien |
| **RK semi-explicite sur groupe de Lie** (§1.3) | vise les Cosserat NON RAIDES à contraintes internes. **Motif d'abord mal posé** : j'avais écrit « ici tout est raide », ce qui juge le NOYAU avec le critère de FRELON — or le point n°6 de la feuille de route (le jumeau de câble) est exactement un problème non raide de ce genre. Reste non fait, pour la bonne raison : la note ne chiffre aucun gain, et le câble n'a pas encore de modèle qui converge sur lequel le mesurer |
| **ADMM** (§3.3) | ses trois motifs sont mesurés absents : LU creux à remplissage borné, Fischer–Burmeister semi-lisse sans lissage dans le même Newton, warm-start déjà là. ADMM paie à 10⁴–10⁵ contacts simultanés, pas à dix |
| **SAP / formulation convexe** (§3.2) | achèterait l'UNICITÉ, que le noyau n'a pas — mais le défaut ne se manifeste pas : `contact.multiples` rend λ = 3mg/2mg/mg exacts et trois billes à e = 1 transmettent v à 0,5000. Même argument que celui qui a différé GGL |
| **LDLᵀ** (§6) | **faux ici** : la tangente d'amortissement est asymétrique à 1,50 en relatif, et le jacobien de pas porte `Gᵀ` en (1,2) contre `G·J_l·c` en (2,1) |
| **précision mixte** (§6) | contre-indiquée : le défaut du 5 septembre (vol S2 cassé à 50,84 s) était un problème d'**ulp** en DOUBLE précision |
| **GPU, ML** (§5–6) | pas la taille, pas les données |

Et la seule qui mérite d'être **lue sans être écrite** : le **principe
variationnel GGL** de Kinon, Betsch & Schneider (*Multibody Syst. Dyn.* 57,
2023, pour le principe ; *Nonlinear Dynamics* 111, 2023, pour les schémas).
Son point porte exactement sur la jointure de deux mécanismes qui vivent ici
SÉPARÉMENT — la stabilisation GGL d'index 2 et la projection sur la variété du
moment : appliquer un intégrateur à structure préservée à GGL *tel quel* ne
donne pas une méthode à structure préservée. Mais l'adopter, c'est **changer
d'intégrateur** — l'α-généralisé n'est ni symplectique ni conservatif, et le
code le dit — donc un fork, pas une brique. La note rappelle aussi
Ge–Marsden, que ce noyau a établi par le calcul et un cran plus bas : pour un
corps rigide `∂E/∂ω = ωᵀ·∂L/∂ω` exactement, la pile est de rang 3 sur 4.

### Ce que la note recommande et qui est DÉJÀ l'architecture

Sa §4.2 conseille une formulation **à un seul niveau** — dynamique et contact
dans un système KKT différentié implicitement, plutôt que dériver à travers
les itérations d'un solveur de contact. C'est la construction du noyau : le
contact non lisse EST une ligne de λ dans le même Newton. Avec une réserve
qu'il faut écrire, parce qu'elle borne le pont adjoint à venir : après tout
pas à contact actif, **le schéma REDÉMARRE** sur l'accélération lisse — donc
`∂z₁/∂z₀` n'y a pas la forme lisse que l'adjoint remonte. Le pont
adjoint-en-temps ↔ noyau (`k_c_m_z`, annoncé dans son propre docstring) devra
traiter ce cas ou déclarer qu'il s'arrête au premier impact.

### Addendum du même jour — deux mesures que le bloc ci-dessus devait à son lecteur

Paul : « pourquoi rien d'implémenté ? ». Réponse honnête : deux de mes refus
étaient posés sans la mesure qui les tranche, et une troisième raison
était mal posée (RK semi-explicite, corrigé dans la table). Les deux
mesures :

**GGL et la projection de Noether COMPOSENT.** Kinon–Betsch–Schneider
établissent qu'un intégrateur à structure préservée appliqué à GGL tel
quel ne préserve plus la structure ; ce noyau porte les deux mécanismes et
personne n'avait vérifié. Modèle de `bancs.moment_projete` (deux corps
libres, rotule, sans gravité, 40 s) :

```
   ggl    moment   dérive L /s   |Φ| max    |Φ̇| max   newton/pas
   False  False    +1,08e-11     9,2e-16    1,07e-06     2,21
   True   True     −5,29e-18     1,8e-15    3,15e-11     2,36
```

L au plancher machine ET Φ̇ au niveau GGL, ensemble. La raison est
simple une fois dite : la projection vit dans le noyau de G, et G ne
change pas sous GGL. Le problème que Kinon résout est celui d'un schéma
qui REVENDIQUE la symplecticité — celui-ci projette. Gardé par la
5ᵉ réfutation de `moment_projete` ; la référence descend de « à lire »
à « sans objet ici, mesuré ».

**Schur, sur le vrai cas raide.** Le garde de `refus_gardes` tourne sur
deux corps sans raideur locale — là où la recommandation ne mord pas.
Mesuré sur le vol S2 de FRELON (ligne d'engrenage à 4e8, servos PD,
aéro, 60 000 pas) : **0 repli SVD**, `jac/pas` 1,05, 0,277 ms/pas.
Asserté dans `vinkulum_s2 demo` côté FRELON. Le refus tient — et il est
maintenant tenu là où il devait l'être.


## HUIT BRIQUES EN UNE SOIRÉE — et deux défauts du noyau trouvés en les jugeant (5 septembre, soir)

Ordre de Paul, dans ses mots : « Implémenté Leishman–Beddoes, sillage libre,
cylindre au contact, checkpointing Griewank, multi-rythme, énergie-moment
Simo–Tarnow, pont adjoint-en-temps, GGL ne traverse pas un point mort. »
C'est la feuille de route du matin, entière, y compris le pont arbitré
« pas maintenant » quelques heures plus tôt. Tout est entré, chaque brique
avec son banc et son contrôle négatif ; rien n'a été recoté pour passer.

| brique | où | ce que le banc mesure |
|---|---|---|
| **GGL au point mort** | `tangent::jacobien_ad`, `verification.srscm`, `bancs.ggl_boucles` | `srscm` sous GGL passe t = 2,724 s (Φ̇ 5,5e-12, écart 1,9e-3 m au pas grossier contre 3,6e-4 en index 3) ; échelle de boucles **12/12** (11 le matin) ; assert resserré, pas desserré |
| **cylindre fini** | `Contact.cylindre`, `contact.cylindres` | quatre régimes (flanc, fond, arête, loin) et trois rotations à **1,2e-14** de la géométrie ; refus sans corps porteur |
| **Leishman–Beddoes** | `Pale.lb`, `bancs.leishman_beddoes` | échelon 5°→20° : surcroît intégré 30,8 ms contre 23,9 pour Øye, C_N^v max 0,53 puis éteint ; tangage 10° ± 10° à k = 0,1 : **CL max 1,684 pour 1,316 statique**, aire ∮CL dα 0,098 ; statique 3e-15, régime attaché 1e-17 |
| **énergie-moment** | `Modele::simule_em`, `bancs.energie_moment` | E et L à **1,0e-14** sur 20 000 pas d'un corps qui culbute ; α-gen 6,7e-6 / 6,3e-4 ; `moment=True` 1,3e-5 / 3,5e-17 ; `energie=True` 4,9e-15 / 6,3e-4 ; ordre **2,00 / 2,00** ; liaison refusée |
| **multi-rythme** | `Modele::simule_multi`, `bancs.multirythme` | chaîne de 8 pendules à h = 1e-3, ensemble vibrant à 300 Hz à h/100 : amplitude à **0,959×** le fin (grossier 0,806×), bout de chaîne à 1,9 mm sur 2 m, **0,81×** le temps du fin, 1 Newton par pas ; liaison qui enjambe refusée |
| **pont adjoint ↔ noyau** | `adjoint_temps.PontNoyau` | double pendule **1,6e-9** de la DF, pendule sphérique tournant **1,9e-10**, flexible (EI) dans l'étendue de la DF (1,4 %) ; sans J_l 5e-3, sans ∂(J_su̇)/∂θ 4e-2, gyroscopique doublé 1,3 |
| **Griewank √N** | `adjoint_temps.NonLineaire.gradient(checkpoints=)` | 55 points de contrôle, pic 55 blocs contre 3 000, gradient identique **au bit** |
| **sillage libre** | `vinkulum.sillage` | λ à **+1,8 %** du momentum, r_v/R(2π) 0,800 contre Landgrebe 0,807, k₂ à −48 %, **k₁ de mauvais signe** (+0,0026 contre −0,0238), Γ = 0 → rien, hélice rigide → pas de contraction, couplage noyau +0,7 % |

### GGL ne traversait pas un point mort — et ce n'était pas le point mort

Le matin, `srscm` cassait sous GGL à l'instant exact d'un pic de
conditionnement ×4 800, et j'ai écrit « GGL ne traverse pas une
configuration quasi singulière ». C'était une lecture. La cause : les
colonnes ζ du jacobien étaient des différences finies, exactes à ~1e-8
relatif, et le conditionnement multipliait ce bruit en un incrément de
Newton faux au pour-cent. Les colonnes sont désormais le produit creux
∂r/∂q · βh²G₀ᵀ, tiré du même jacobien AD (une liste `kqt` des blocs « par
unité de pose », remplie seulement quand GGL est armé). Plus une copie du
modèle par colonne, plus d'incrément à régler. Le README et le banc portent
la rétractation ; l'assert de l'échelle est resserré à 12/12.

### Le gyroscopique était compté DEUX fois dans K et C

Trouvé par le pont adjoint, qui ne fermait pas : `Modele::raideur`
retirait ω × J_sω « à la main » alors que `forces()` le porte déjà, et
`amortissement` ajoutait de même [ω]×J_s − [J_sω]×. Mesuré sur un corps
libre : K_θθ et C_ωω à **exactement 2,000×** l'analytique. Le résidu de
Newton n'a jamais doublé (il ne passe pas par `forces()`) — donc aucune
trajectoire n'était fausse ; seules les tangentes EXPOSÉES l'étaient :
`k_c_m_z`, `modes_complexes`, `diag_modal`, Floquet, l'adjoint. Corrigé ;
`verification.gyroscopique_tangent` (n° 29) garde la nutation libre
ω₃J_a/J_t à 1e-6 — 9,54930 Hz retrouvés. La « limite gyroscopique » que
`modes_complexes` déclare sur la toupie rapide (σ = ∓145) en était très
probablement la cause : à rejouer, écrit à la feuille de route.

### Ce que les bancs ont refusé avant de passer

- **Énergie-moment** : relire Π de (R, ω) à chaque pas faisait dériver E
  et L en **N²** (4e-13 à 2 000 pas, 7,5e-11 à 20 000, 5,6e-9 à 200 000) —
  l'aller-retour J·Rᵀ·R·J⁻¹ sur une R orthonormale à N·ε près. Π est porté
  d'un pas à l'autre et R re-orthonormalisée : 1e-14 à 200 000 pas. Puis la
  lecture de L elle-même dérivait avec un corps qui s'éloigne à 2 m/s (le
  terme m·r×v à 400 m) : mesure de la LECTURE, pas du schéma — le banc
  tourne en rotation pure, la translation est jugée sur E sous gravité.
- **Multi-rythme**, deux fois. Le premier banc mettait la poutre RAIDE
  entre les partitions : un partenaire gelé sur 35 kN/m pendant un pas
  lent, c'est 9 N de force parasite, le pendule ne bougeait plus — la
  raideur doit vivre DANS la partition rapide, c'est la condition du
  multi-rythme et pas une limite du code. Puis le bout de chaîne sortait à
  **4 cm** : l'α-généralisé REDÉMARRÉ à chaque sous-pas (a₀ = u̇₀) n'est
  plus que d'ordre 1. L'état (u̇, a, λ) survit désormais par partition
  (`schema_init`/`schema_fin` sur `simule`) : 1,9 mm.
- **Leishman–Beddoes** : le premier banc lisait C_N^v = 0 et un surcroît
  égal à Øye — la section voyait α = −a (signe de la normale), le
  tourbillon était négatif et le `max` le cachait. Signe retourné, métrique
  changée : le pic à 0⁺ est le saut attaché, identique pour les deux
  modèles ; ce que LB allonge est le surcroît INTÉGRÉ.
- **Sillage** (rapport de l'agent) : la relaxation de Picard DIVERGE en
  stationnaire (instabilité d'appariement des hélices, physique), remplacée
  par un Newton amorti sur le point fixe ; `scipy.optimize.root` n'a pas
  quitté le point de départ.
- **Cylindre** : en passant, boîte et maillage tombaient dans la branche
  « plan » du CCD avec la normale par défaut — exclus explicitement.

### Ce qui reste hors domaine, et c'est dit

Énergie-moment avec liaisons (gradient discret de Φ), multi-rythme sans
restriction du système (pas plus rapide) ni corrector ni aéro, sillage en
avancement et nappe interne, LB sur CD/CM et constantes calibrées, pont
adjoint sans contact et sans Rust, revolve binomial. Tout est à la feuille
de route du README, avec le chiffre qui le motive.

Vérification : 29 contrôles, 174 s ; `cargo test` 6/6, clippy à zéro ;
`multicorps demo` et `vinkulum_s2` verts côté FRELON (0 repli SVD sur 6 666
images).


## L'AUDIT DU JACOBIEN — sept termes manquants ou faux, et le refus de `modes_complexes` qui tombe (5 septembre, nuit)

Ordre de Paul : « Analyse le noyau pour corriger ses défauts. » On ne
corrige pas ce qu'on ne mesure pas : le premier geste a été un INSTRUMENT,
`Noyau.audit_jacobien` — le jacobien AD de Newton contre une différence
finie centrée du MÊME résidu, bloc par bloc (dyn | con | ggl) × (u̇ | λ | ζ),
sur un corpus de dix-sept modèles qui couvre chaque type d'élément. Le pas
de la DF se cale sur ce que chaque inconnue DÉPLACE (u̇ et ζ n'agissent sur
Φ que par la pose, à βh²c près) : à 1e-6 sur u̇ la pose bougeait de 1e-13 et
les colonnes de contrainte étaient du bruit à 1e-3 — la leçon des colonnes ζ,
rejouée sur l'instrument lui-même. Plancher atteint : 1e-10.

Ce qu'il a rendu, dans l'ordre où c'est tombé, chaque chiffre AVANT → APRÈS :

| défaut | mesure | correction |
|---|---|---|
| raideur de POSITION des contacts par pénalité absente du jacobien (« quasi-Newton, comme la pénalité qui ignore déjà sa raideur ») | 7e-4 à h = 1e-3, en k·h² au-delà | `Contact::tangentes` à 24 duaux (translation comprise) → 1e-7 |
| ∂F_t/∂λ du frottement non lisse absent | **0,22** relatif | colonne λ = force par unité de λ → 1e-11 |
| ∂(G·u)/∂q des lignes GGL et non holonomes « laissé de côté, O(h) » | 5e-2 sur un cardan sous GGL | `Elem::d_gu`, duaux emboîtés composés exp(sû)·exp(δ̂) ; `srscm` sous GGL : 4,6 → 2,9 Newton/pas, Φ̇ 5,5e-12 → 1,7e-16 |
| G de l'engrenage et de la vis « à la main » (n_w), faux hors de l'axe | 7e-4 | forme fermée dθ/dε = [½c(tr E·n − Eᵀn) + s·a]/(s²+c²) → 1e-10 ; les duaux emboîtés 18×18 essayés d'abord coûtaient ×15 sur le pas |
| colonne u̇ des lignes en VITESSE (nh, contact) passée par J_l comme une ligne en pose | 4e-4 | sans J_l → 5e-10 |
| ∂Φ_t/∂q d'une cible rhéonome sous GGL | 5e-4 | `Elem::d_phit` → 7e-12 |
| **Hessienne des éléments par duaux emboîtés en UNE exponentielle** (vis, cardan, puis engrenage) : il manquait ½[Gᵀλ]× (BCH) | K[4,3] 0 contre −0,40 en DF (test Rust) ; cardan 8e-7 | `eval_k` compose exp(δ₁)·exp(δ₂)·R, δ₁ à gauche → 4e-14 |

Reste laissé de côté, mesuré et déclaré : la raideur λ·∂G/∂q du contact non
lisse (7,7e-6). Le contrôle 30 de `verification` garde le corpus à 1e-5.

### `modes_complexes` refusait l'état tournant — la cause n'était pas celle écrite

Le refus du 4 sept. accusait « le transport d'une perturbation qui ne reste
pas admissible » et mesurait des fréquences fausses de 47 à 92 % sur la
toupie rapide. Le double compte gyroscopique corrigé ce soir était le
suspect ; forcé, le calcul rendait encore 0,004 et 51,6 Hz pour 0,27 et
8,41. La cause : le système d'état posait δq̇ = δu. Pour R = exp(δθ)·R₀(t)
qui tourne, δω = δθ̇ − ω₀×δθ, et la translation d'un corps contraint suit
δr = δθ×r₀. Écrit en Python d'abord (`Zᵀ(M δu̇ + C δu + K δq)` avec
δu = Zξ̇ + Sξ, S = Ż − T_L Z) : **0,268524 et 8,412654 Hz, les deux racines
exactes du trinôme, σ 1e-10** ; le corps libre à 95,49297 Hz = ω₃Ja/Jt ;
l'état au repos inchangé. Porté dans le noyau, la base admissible dérivée le
long du mouvement par différence finie alignée. Le refus tombe ; le contrôle
29 garde les deux racines à 1e-6.

### Le rouge de `bancs rapide` avait une cause, et c'était le balayage

`sensibilite.demo(rapide)` exigeait que le semi-analytique se dégrade d'un
facteur 10 sur un balayage du pas qui s'arrêtait à 1e-12 — là où, sur
4 éléments, le mur d'annulation COMMENCE (0,04 % ; 3 % à 1e-13, 14 % à
1e-14, mesuré). Le balayage va jusqu'au mur ; rapport 1 370 pour 10. Aucun
seuil n'a bougé.

### Ce que ça change côté FRELON

Le vol S2 est identique (régime, couple, cyclique au centième, 0,279 ms/pas)
et les 153 Hz de la boucle de commande sont **rejoués, inchangés** — ce
modèle-là est linéarisé au repos, le double compte n'y entrait pas.


## PALES ÉLASTIQUES, SILLAGE EN AVANCEMENT, ET LA PREMIÈRE CAMPAGNE CONTRE DES ESSAIS (5–6 septembre, nuit)

Ordre de Paul : « Ajoute les pales élastiques couplées à l'aéro et le sillage
en avancement. Une fois ces deux points implémentés avec succès, lance une
campagne de validation contre des essais et modèles publiés et reconnus. »
Les deux briques sont entrées, chacune avec ses juges, et la campagne
(`python -m vinkulum.validation` → `docs/VALIDATION.md`) rend **26 lignes
sur 27 dans leur tolérance**, la 27ᵉ publiée avec sa cause.

### La pale élastique (`vinkulum.pale_elastique`)

Chaîne de poutres géométriquement exactes dont chaque nœud porte son tronçon
de `pale` : la section suit la rotation du nœud, donc torsion et battement
élastiques entrent dans l'incidence par la cinématique. Jugée :

```
   éventail de Wright et al. 1982 (Ω̄ 0 → 12, 20 éléments)     −0,12 à −0,18 %
   amortissement aéro du 1er mode contre le Lock modal          −1,7 %
   flèche du bout contre une poutre TENDUE intégrée (mêmes charges)  +0,8 %
   μ = 0,2 : périodique à 1,9 %, amplitude 1/rev 52 mm sur 1 m (pale souple)
```

Trois pièges payés : `poutre` prend sa longueur au repos à la création (les
nœuds créés déjà étirés donnaient une poutre SANS tension et un éventail
plat — ω̄₁ 3,50 à Ω̄ 12) ; la charge linéique de la poutre intégrée devait
diviser par la demi-longueur aux bouts (+12 % → +0,3 %) ; et « la pale
élastique porte moins » était une croyance — mesuré −0,11 % puis +0,23 % selon
la discrétisation, c'est-à-dire zéro : section symétrique, pas de couplage.

### Le sillage libre en avancement (agent, `sillage.SillageAvancement`)

Grille périodique lâcher × âge, relaxation sous-relaxée (elle CONVERGE en
avancement, là où elle divergeait en stationnaire — mesuré des deux côtés),
harmoniques d'inflow imposées au noyau (`pose_inflow_harmoniques`, neuf).

```
   k_x au centre, μ 0,15 : 1,082  (Drees 1,084 → −0,2 % ; Coleman 0,854 → +27 %)
   k_x au centre, μ 0,30 : 1,118  (Drees 1,068 → +4,7 % ; Coleman 0,963 → +16 %)
   trim T = 3,33 N, Mx = My = 0 : θ1c 0,83° (uniforme) → 3,29° (sillage), Δ +2,46°
```

Le trim latéral qui bouge de 2,5° est le fait qu'un code de métier veut
voir. Deux k_x existent et ne disent pas la même chose : sur tout le disque
1,60 (champ non linéaire, upwash au bord avant), au centre 1,08 — Coleman et
Drees sont des gradients AU CENTRE. μ < 0,12 est hors des deux solveurs,
déclaré.

### La campagne (`docs/VALIDATION.md`)

| référence | ce qui est jugé | résultat |
|---|---|---|
| Wright/Hodges 1982 | ω̄₁(Ω̄) poutre tournante, 9 points | −0,12 à −0,18 % |
| Bramwell/Johnson (forme fermée) | β₀, β1c, β1s articulé μ 0,2 | +0,03 / −0,08 / −2,2 % |
| Johnson ch. 9, Lock modal | ζ du battement élastique | −1,7 % |
| **Benedict et al. 2015 (MESURÉ)** | FM de 4 micro-rotors, polaires NeuralFoil, **inflow du sillage libre** | −2,9 / −3,9 / **−12,7** / −5,5 % |
| Drees 1949 (Coleman publié) | k_x du sillage en avancement | −0,2 / +4,7 % |
| MBDyn, MSD 2016, Bennett, Kapitza, Andrews | rejoués | verts |

**Ce que Maryland a dit, et c'est le résultat de la nuit** : avec l'inflow
UNIFORME du noyau, les quatre FM sortaient à +12, +12, +15 et +35 % de la
mesure (la plaque mince à −1 %) ; le ki effectif 1,34 que la gate V-2 de
FRELON avait calibré sur ces rotors disait déjà que c'est l'inflow qui
manque. Un **profil radial imposé** (`pose_inflow_profil`, neuf) porté par le
sillage libre stationnaire ramène trois rotors sur quatre sous 6 % d'un
essai publié. Le quatrième (NACA 4512, −12,7 %) reste sur la polaire :
NeuralFoil à Re 30 000 sur 12 % d'épaisseur, sensible à N_crit, bande {5, 9}
publiée.

### Python 3.14, free-threaded, et PyO3 0.29

Consignes de Paul dans la nuit : multithread sauf exception justifiée,
3.14 sauf obstacle, PyO3 dernière version stable. Fait : PyO3 0.23 → 0.29
(`allow_threads` → `detach`), module `gil_used = false`, l'interpréteur
détaché dans `simule` et `modes_complexes`. Trois roues (3.13, 3.14, 3.14t),
`verification` verte sur les trois, quatre fils Python sur des `Noyau`
indépendants : **×3,5 / ×3,0 / ×3,4**. Piège : sous free-threaded, un objet
Python partagé entre fils (une liste mutée dans un `monte`) se court-circuite
— le harnais, pas le noyau.

## LE NON-COUVERT INSTRUIT — la nappe rend k₁ à Landgrebe et DÉFAIT Maryland (6 septembre)

Ordre : « implémente le non couvert », c'est-à-dire les quatre trous que la
campagne de la nuit déclarait — nappe interne du sillage, μ < 0,12, une
mesure d'inflow en avancement, un essai de décrochage dynamique. Les quatre
sont instruits ; deux ferment, deux rendent un résultat qui n'est pas celui
qu'on attendait, et c'est celui-là qui compte.

### La nappe interne (`Sillage(nappe=True)`, le défaut)

Compromis CAMRAD : Γ(r) de l'élément de pale au pas qui rend la poussée,
filaments internes ΔΓ_k et tourbillon de pied **prescrits** en hélice au
pas de l'inflow LOCAL lu au disque (deux passages), marginal libre par
Newton. Mesuré sur le rotor FRELON, Δζ 15°, 3 tours :

```
   k₁ (ζ < π)        +0,0026  →  −0,0195     Landgrebe −0,0238   bon signe, facteur 1,2
   k₂ (π < ζ < 2π)   −0,043   →  −0,076      Landgrebe −0,083    −8 %
   r_v(2π)            0,800   →   0,833      Landgrebe 0,807     +3 %
   asymptote          0,777   →   0,781      Landgrebe 0,78
   λ / momentum       1,018   →   1,083
```

Trois contrôles négatifs : sans nappe, le modèle du 5 sept. au chiffre près
(k₁ redevient positif — c'est asserté) ; Γ = 0, rien ; hélice rigide, pas de
contraction. Essayés et refusés, mesurés : la nappe au pas de Froude seul
(k₁ −0,016 mais upwash intérieur sur les rotors chargés), un quasi-Newton au
jacobien sans la nappe (plancher 2e-7 : le marginal jeune est trop près),
la nappe fine (142 s pour un k₁ moins bon que la grossière), 4 éléments.
Prix : 20 → 80 s.

### Et Maryland, rejugé avec la nappe, SORT de la bande

```
                        sans nappe (5 sept.)   AVEC nappe (6 sept.)
   plaque cambrée           −2,9 %                 +9,9 %
   NACA 4504                −3,9 %                +15,8 %
   NACA 4512               −12,7 %                +19,8 %
   NACA 0012                −5,5 %                −14,5 %
```

**La géométrie la plus proche de Landgrebe rend la FM la plus loin de
Benedict.** Sur ces rotors très chargés (C_T/σ 0,085–0,145) la nappe
prescrite reste au ras du disque et y induit un upwash intérieur : moins de
puissance induite, FM trop haute — et l'écart change de signe. Je ne choisis
pas la ligne qui passe : `validation` juge désormais les DEUX sillages,
aucune tolérance ne bouge, trois lignes sont rouges, et le désaccord est
déclaré non résolu dans le docstring de `cas_maryland`. Deux lectures
restent ouvertes : l'accord sans nappe était une compensation (pas de
nappe = pas d'upwash = un ki qui tombe juste), ou la nappe prescrite n'est
pas une nappe libre. Ce qui trancherait — une nappe libre (×6 sur les
inconnues) ou une mesure d'inflow au disque à Re 3e4 — n'existe pas encore.

### μ < 0,12 par continuation (`sillage.continuation`)

Depuis μ 0,15, la géométrie convergée du pas précédent, relaxation réduite
0,5 → 0,25 → 0,1. Le plancher est **mesuré** : 0,12 à trois tours libres,
0,10 à deux ; 0,08 et 0,06 ne convergent pas ; 0,04 « converge » sur un k_x
NÉGATIF — un résidu qui descend ne suffit pas, et c'est asserté que la
continuation tienne jusqu'à 0,12 avec k_x dans la bande Coleman–Drees. Sous
le plancher tous les tours interagissent : c'est un Newton sur la grille
périodique entière qu'il faudrait, comme en stationnaire. Chantier nommé.

### Elliott 1988 et le S809 (Re ~1e6 — PAS le Reynolds de FRELON)

Elliott, Althoff & Sailey (NASA TM 100541/100543), inflow mesuré au vélocimètre
laser sur un rotor à 4 pales, tables 3 parsées depuis les scans (121 et
146 points) : RMS/σ_mesure 0,285 contre 0,333 pour Drees à μ 0,15 (k_x −11 %),
0,744 contre 0,911 à μ 0,30. Le sillage bat le modèle linéaire là où il est
jugé. Ramsay, Hoffman & Gregorek 1995 (OSU/NREL), S809 en tangage 5,5° à
k 0,077, Leishman–Beddoes + Wagner : CL max à −9 / +2 / −31 % (8 / 14 / 20°),
CL min à −5 / +10 / +0,3 %. Le 20° est hors bande : les constantes de LB sont
celles du NACA 0012, jamais calibrées sur la section, et le tourbillon de
bord d'attaque du S809 en décrochage profond ne se règle pas avec.

**Bas Reynolds, dit ligne par ligne.** Paul l'a rappelé : les validations à
haut Re ne couvrent pas FRELON (Re 5e4–8e4). Le tableau porte désormais le
Reynolds de chaque ligne ; seul Maryland (3e4) est à bas Re, et **aucun
essai publié de décrochage dynamique ni d'inflow mesuré n'existe à ce
Reynolds** — cherché (McCroskey 1982 est un scan de courbes sans table). Le
non-couvert a changé de nature : ce n'est plus un modèle absent, c'est une
donnée absente, et ça s'écrit tel quel.

L'agent qui portait la nappe est mort sur la limite de session avant de
lire son propre résultat Maryland ; c'est sa sortie stockée qui l'a donné.

## DEUX MESURES DE DÉCROCHAGE DYNAMIQUE À BAS REYNOLDS, ET LE HARNAIS QUI MENTAIT À GRAND ANGLE (6 septembre, soir)

Paul : « fais un effort ». Le corpus Mars de la NASA refusé le matin (bas
Re mais haut Mach, coaxial), restait à trouver ce qu'on avait dit
introuvable — une mesure de section instationnaire à Re 1e4–5e4. Deux
existent, libres d'accès une fois qu'on cesse de demander aux éditeurs :

- **Jantzen, Taira, Granlund & Ol 2014** (Phys. Fluids, DTIC ADA603996 via
  la Wayback Machine) : plaque plane 2D en rampe 0 → 45° au bord d'attaque,
  tunnel hydraulique AFRL, **Re 20 000**, C_L(t) en tirets sur une figure
  VECTORIELLE — les deux courbes (C1 : K π/8, C6 : K π/48) sont extraites du
  SVG et calibrées sur les graduations, 642 et 676 points, `jantzen2014.py` ;
- **Kim, Chang & Sohn 2017** (IJASS, koreascience) : NACA 0012 en tangage
  sinusoïdal 0° ± 6° autour de c/4, **Re 23 000**, k 0,1 / 0,2 / 0,4 / 0,76,
  C_L par pression instationnaire. Figure rastérisée mais repères IMPRIMÉS
  (C_L max, pente, flèche de sens), `kim2017.py`.

### Le harnais projetait la portance sur z, pas ⟂ au vent

En regardant la rampe à 45° : quasi-statique à 0,67 là où la polaire construite
donne 0,95. Le harnais de section (celui du S809) lisait `l_vec · z`, c'est-à-dire
L·cos α — −3 % à 13,5°, **−10 % à 25,5°, −29 % à 45°**. Le sondage du S809
(« |L|/q = table à 1e-3 ») vérifiait le MODULE, pas la projection. Corrigé
(e_L = (−sin α, 0, cos α)) ; les lignes S809 bougent (8° : −9 → −7 % ;
14° : +2 → +6 % ; 20° : C_L max 1,12 → 1,17, toujours −28 %). Faute de harnais,
pas de noyau : le banc LB de `bancs.py` projette juste.

### Wagner est juste ; ce qui casse à bas Re est ailleurs

Avant de juger, le circulatoire du noyau contre Theodorsen sur polaire
linéaire, aux k de l'essai : |C| et phase à **≤ 2 %** de 0,05 à 0,76 (table
dans le docstring de `cas_kim`). Puis Kim :

```
   k 0,10 : mesuré 0,65 · quasi-statique 0,54 · Wagner 0,37 · LB+Wagner 0,14
   k 0,76 : mesuré 0,41 · quasi-statique 0,54 · Wagner 0,05 · LB+Wagner 0,12
   sens de la boucle : anti-horaire à 0,1, HORAIRE à 0,76 — RENDU aux deux
```

Deux causes, mesurées. **(1)** NeuralFoil à Re 2,3e4 porte le saut de la
bulle laminaire (C_L 0,075 à 4°, 0,535 à 6°) et le noyau lit la table à
l'incidence EFFECTIVE de Wagner : à k 0,76 elle plafonne à 3,4°, sous le
saut. L'essai dit l'inverse — « les événements de couche limite disparaissent
à k ≥ 0,4 », la portance redevient quasi linéaire. Passer une non-linéarité
statique de bas Re à travers un retard circulatoire l'AMPLIFIE au lieu de la
lisser. **(2)** La masse ajoutée (∝ α̇) n'est pas modélisée : à k 0,76 elle
vaut ~π·k·α_a ≈ 0,25 — la moitié du C_L max mesuré. Rien de ces deux-là ne se
règle par une constante de LB. Le sens de la boucle, lui — le renversement
entre 0,4 et 0,76 que Theodorsen place entre 0,1 et 0,2 — est rendu.

Jantzen C1 (K π/8) : C_L moyen sur [2 ; 8] **2,02 pour 1,73 mesuré (+16 %)**,
quasi-statique 0,95, CN post-décrochage ±15 % → 1,90–2,11. C6 (K π/48) :
**1,83 pour 2,35 sur [3 ; 6] (−22 %)** — le pic pendant la rampe lente, LB
sous-porte — puis **2,05 pour 1,43 sur [6 ; 8] (+44 %, ÉCART)** : la plaque
mesurée a fini de lâcher son tourbillon et retombe, le modèle tient un
plateau. Le pic non circulatoire (C_L 6 à t 0,85) est exclu et déclaré ; la
polaire statique de plaque est CONSTRUITE (le papier n'en publie pas) et
c'est l'entrée faible. Kim aux quatre k : C_L max à −69 … −78 %, sens rendu
4/4. **41 / 51 lignes** dans `docs/VALIDATION.md`, 21 à bas Reynolds.

**Ce que ça change à la priorité** : la donnée bas-Re existe. Ce qui manque
est dans le noyau — la masse ajoutée et la façon de composer retard et
polaire non linéaire — et c'est un chantier de modèle, pas de données.

## CLORE L'ÉQUATION — masse ajoutée, 3/4 de corde, retard linéaire, et une cible en vitesse (nuit du 6 au 7 septembre)

Paul : « il faut clore l'équation ». Les deux manques nommés le soir —
masse ajoutée absente, retard × polaire non linéaire — sont dans le noyau,
avec ce que leur mise au banc a exigé de découvrir sur le harnais.

### Trois briques dans `Pale`

- **`masse_ajoutee`** : L_nc = πρb²·ẇ, w la vitesse normale de l'air vue au
  MILIEU de corde (= ḧ + Vα̇ − b·a·α̈ chez Theodorsen), ẇ par différence
  arrière sur le pas — ordre 1, implicite en vitesse (l'AD porte πρb²/h dans
  la tangente). Pas de moment non circulatoire, déclaré.
- **`alpha_34`** : l'incidence circulatoire lue au 3/4 de corde — le
  b(½ − a)α̇ de Theodorsen. Sans lui, une section qui tangue ne voit pas son
  taux de tangage ; c'est vrai aussi d'une pale de rotor en cyclique, et
  c'est un chantier FRELON nommé (défaut à False : le vol S2 ne bouge pas).
- **`retard_lineaire`** : Wagner ne retarde que a₀(α − α₀) ; Δ(α) = c81 −
  linéaire est lu à l'incidence instantanée. a₀ par la table ou `a0_rad`
  imposé (2π sur Kim : la table NeuralFoil à Re 2,3e4 n'a pas de partie
  linéaire, c'est une bulle).

Banc `theodorsen_complet` — pivot c/4, C_L/(2πα) = C(k)(1 + ik) + ik/2 − k²/4 :

```
   k 0,10 : |H| 0,843 pour 0,848 (−0,5 %)   phase −2,3° pour −2,6°
   k 0,40 : 0,709 pour 0,711 (−0,3 %)       +22,4° pour +23,6°
   k 0,76 : 0,828 pour 0,850 (−2,6 %)       +53,0° pour +53,6°
   masse ajoutée seule 0,474 (0,488) · 3/4 seul 0,701 (0,717) · ni l'un ni l'autre 0,559 (|C| 0,571)
```

Chaque terme rend SA part, et leur somme la forme fermée : c'est ce qui dit
que ce sont deux termes et non un réglage.

### Le harnais a menti deux fois, et les deux fois la polaire symétrique le couvrait

**Une cible en position fait osciller la vitesse.** Sur la rampe de Jantzen
la masse ajoutée rendait un pic de C_L **2 974** pour 6 mesuré : la cible de
rotation était imposée en POSITION (index 3), la vitesse du corps est alors
une quantité dérivée qui oscille d'un pas à l'autre, et une différence
arrière de ẇ en fait un pic par pas. Un nœud tous les cinq pas n'y change
rien (0,72 pour 0,49 sur polaire linéaire). Imposée en VITESSE — liaison
non holonome `nh=True` avec cible — ω est exactement la pente de la table
sur le pas et la différence arrière rend α̈ en moyenne. Sauf que la ligne
non holonome ignorait sa cible (résidu G·u = 0, sans Φ_t) : la section ne
bougeait pas. Corrigé dans le noyau, docstring de `liaison` à jour.

**Un axe faux et une lecture fausse s'annulent sur une polaire symétrique.**
Avec es = y, ec = x la normale de section est −z : une rotation nose-up
autour de −y donne −α à la section, et la sonde d'inflow (+z) lit −L. Sur
NACA 0012, plaque, polaire linéaire : rien ne se voit, tous les contrôles
passent. Sur le S809 cambré : 0,257 pour 0,447 à 3°. Axe +y, lecture −sonde,
et `_tangage_signe` juge désormais sur le S809 à 3° ET 13,5° — une polaire
symétrique ne peut plus servir de contrôle de signe.

### Ce que ça ferme, et ce que ça ne ferme pas

```
   Kim k 0,76   0,41 mesuré :  0,05 (Wagner, table retardée)  →  0,399  (−3 %)
   Kim k 0,10   0,65 mesuré :  0,36                            →  0,43   (−34 %)
   Jantzen C1   pic non circulatoire 6,05 :                        4,90   (−19 %)
   Jantzen C1   moyenne [2 ; 8] 1,73 :                             1,92   (+11 %)
   sens de boucle Kim : 3 sur 4 — le modèle complet RENVERSE entre k 0,2 et
   0,4 (comme Theodorsen, la masse ajoutée mène en phase), la mesure entre
   0,4 et 0,76. Le Wagner « table à α retardée » rendait 4/4 : c'était le
   retard artificiel de la bulle lue en retard qui tombait juste.
```

Campagne complète : **42 / 52**, 22 à bas Re (Kim −34 / −27 / −28 / −3 % de
k 0,1 à 0,76 ; Jantzen C1 +11 %, pic −19 %, C6 −6 % sur [3 ; 6] et +31 % sur
[6 ; 8] — le modèle tient un plateau que la plaque mesurée quitte ; S809
−6 / +7 / −27 %). À k 0,76 l'équation est close : la moitié du C_L était la
masse ajoutée, et elle y est. À k 0,1 le manque est la BULLE laminaire — la mesure porte PLUS
en descente qu'en montée (réattachement retardé), un mécanisme que ni
Wagner, ni LB, ni le retard linéaire ne contiennent — et c'est déclaré, pas
maquillé par une constante. LB reste à 0,14 sur Kim : son f_st, inversé de
Kirchhoff sur une polaire à bulle, dit « décollé » à 4° ; ce n'est pas le
modèle qu'il faut à cette décade.

## PETERS–HE, LE NIVEAU QUI MANQUAIT — et le ½ que le rapport ne disait pas (7 septembre, petit matin)

Paul a transmis un rapport de master (S. Li, UC Davis 2020) : l'inflow
dynamique à états finis de Peters–He, validé sur le rotor 2MRTS d'Elliott —
le même essai que la ligne G, avec ses réglages de trim mesurés (table 5 :
θ1c, θ1s, α_arbre), que `elliott1988.py` avait lus à l'identique. Le
rapport donne φ, H, L̃, Γ, V, τ (éq. 28–41) ; il manquait deux choses que
Ferreira et al. 2021 (TEMA, libre) fournissent : le « 2 » au dénominateur du
Γ impair, et le **½ devant τ**. Sans ce ½ la convergence à N états rend
exactement 2× le momentum — mesuré ici (λ̄ → 1,006·C_T/V à 8 états) AVANT
de le lire, et c'est ce qui a fait chercher le texte.

`vinkulum.peters_he` : états (m, n) par la méthode des tables, K = 2H/π, L̃
en forme fermée, `stationnaire` (point fixe sur les Ṽ de He, sous-relaxé —
en stationnaire pur λ ← c/λ oscille), `avance` (Euler implicite sur
[K]α̇ + Ṽ L̃⁻¹α = ½τ), une théorie des tranches de rotor rigide trimée à C_T.
Trois formes fermées SORTENT du modèle et sont assertées :

```
   1 état, m = 0            λ̄ = 9/8 · C_T/(2V)        (la troncature, pas une erreur)
   8 états, m = 0           λ̄ = 1,006 · C_T/(2V)      momentum
   3 états, charge uniforme k_x = (2π/3)·X = 2,094 X   Pitt–Peters 15π/32 = 1,473 X, Coleman X
```

Sur Elliott à μ 0,15 (au disque, rotor rigide, θ₀ trimé 6,72° pour 9,37
mesuré — pas de battement) : RMS/σ **0,386** — Glauert 0,677, Drees 0,333,
sillage libre 0,285 une corde au-dessus ; à μ 0,30 : 1,030 (Drees 0,911,
sillage 0,744) — au disque le λ₀ de PH vaut 0,0092 là où la sonde mesure
0,002, et c'est la hauteur de sonde qui pèse, déclarée. k_x ajusté 1,05 pour
1,67 mesuré : à 3 états et charge uniforme le modèle rend 1,65, ce sont les
modes radiaux supérieurs et les harmoniques de charge qui le rabattent —
vérifié en les coupant un à un (`M1Q1` uniforme 1,65 · `M4Q6` uniforme 1,10 ·
`M4Q6` complet 1,05). Publié, pas asserté au-delà de « bat Glauert ».

**Non couvert, déclaré** : le couplage au noyau. `Inflow` n'accepte qu'un
profil radial et une harmonique 1/rev ; imposer un champ w(r̄, ψ) complet à
une pale est la brique suivante — c'est elle qui mettrait Peters–He sous la
tête S2, à la place du Pitt–Peters à trois états.

### Et la carte d'inflow, pour que ça se branche (même matin)

`pose_inflow_carte(i, rbar, psi, w, centre, e1, vtip)` : un champ w(r̄, ψ)
complet imposé à une pale, bilinéaire, périodique en ψ — là où le noyau
n'acceptait qu'un profil radial et une harmonique 1/rev. Banc `inflow_carte`
sur le rotor de `vinkulum.rotor` : une carte constante rend la poussée de
`pose_inflow` uniforme à 1e-6 ; une carte v₀ + r̄(v1c cos ψ + v1s sin ψ) rend
les moments de `pose_inflow_harmoniques` à 0,06 % (l'interpolation linéaire
d'un cosinus sur 5°, (Δψ)²/8) ; tournée de 180° elle les retourne. Ce qui
manque encore pour fermer la boucle Peters–He ↔ noyau : la projection des
charges de pale sur les φ_n^m côté Rust (le noyau n'expose que poussée et
couple par disque). Nommé.

## LE JACOBIEN SURVIT AU PAS — Newton modifié inter-pas, et `expm` à valeur nulle (7 septembre, matin)

Ordre : « optimise les performances du noyau ». Instrumenté AVANT de toucher
(règle du 3 sept.) : `chronos()` et `perf` sur les deux charges qui comptent.

```
   Princeton (11 corps, 10 poutres, 33 amortisseurs, 66 inconnues, h 5e-4)
      jacobien 81 % du pas — m_mul<Dual<12>> 42 %, expm<Dual<12>> 14 %
   tête S2 FRELON (12 corps, 71 contraintes, 143 inconnues, h 1e-4)
      jacobien 34 % · LU dense + descentes 50 % · résidus 11 %
```

Le jacobien se calculait UNE fois par pas (Newton simplifié du 4 sept.) et
mourait avec lui. Il porte h et la dimension du système ; tant qu'ils ne
bougent pas, rien ne l'oblige à mourir.

**Deux leviers, tous deux mesurés, le second bit-identique.**

- **Le jacobien factorisé est GARDÉ d'un pas à l'autre**, avec le juge qui
  existait déjà à l'intérieur du pas : il sert tant que le résidu chute d'au
  moins 10× par descente, sinon il se rafraîchit. Il tombe si h ou `dim`
  changent (adaptatif, contact activé, masque rejoué) ; en dessous de
  `JAC_HERITE_MIN = 48` inconnues il n'est jamais gardé — le jacobien y est
  bon marché et le quadratique y achète l'invariance au plancher (banc
  `invariances`, 30 inconnues : 5e-13 m ; sous jacobien hérité 1e-12).
- **`expm` sur duaux du premier ordre à valeur nulle** — le cas de TOUTES les
  tangentes du noyau (`pose_pert`) — rend `I + [w]×` sans racine ni k² : c'est
  ce que la série calculait déjà (a = 1, b·k² = 0), 27 produits de duaux en
  moins. Refusé pour les duaux emboîtés, où le terme mixte ½s[u, δ] vit.

**Trois règles qu'il a fallu ajouter, chacune sortie d'un banc qui a rougi.**

1. *Le quatre-barres flexible à h 4e-3 lâchait à t = 0,016.* Un premier pas
   hérité accepté à 2× de chute sortait du bassin. **Un pas sous jacobien
   hérité n'est pris que s'il fait les 10× d'un Newton sain** ; sinon on RESTE
   au point courant et on rafraîchit — ni rebroussement (halver un pas faux ne
   le rend pas juste), ni pas forcé. Partout où l'hérité n'est pas aussi bon
   que le frais, la trajectoire de Newton est celle du frais.
2. *Invariance de repère 3e-10 m au lieu de 5e-13 ; DF de sensibilité modale
   sur le 6ᵉ EI à 3,5 % au lieu de 0,6.* Sous jacobien hérité la convergence
   est linéaire : le dernier itéré s'arrête AU seuil, là où le quadratique le
   dépassait de plusieurs décades. **Sous le seuil on continue de descendre
   avec le même jacobien** (une descente coûte un résidu, pas un jacobien),
   toute descente est prise, et un pas qui ne descend plus EST le plancher.
   Résultat : sensibilité à 0,26 % — mieux que l'ancien 0,64.
3. *Princeton remonté à 5 s.* La règle 1 (10× ou refus) refusait la descente
   vers le plancher de la règle 2 et jetait le jacobien. Les deux régimes sont
   séparés par `sous_tol`.

Essayé et **refusé, mesuré** : « la stagnation (rapport > 0,9) ne compte que
sous jacobien frais » — sur la tête S2 le plancher est atteint sous jacobien
hérité, la règle y forçait un jacobien par pas (0,38 ms/pas, ×1,7 PIRE) ; et
`lto = "fat"` + `codegen-units = 1` : −3,7 % sur S2, −2 % sur Princeton, build
11 → 65 s. Pas pour ça.

```
                              avant        après       jacobiens/pas   Newton/pas
   Princeton 8 s              9,22 s       1,98 s      1,00 → 0,024    2,57 → 2,33
   tête S2 (vol 6 s)          0,281 ms/pas 0,219       1,05 → 0,50     2,57 → 4,9
   verification               171 s        104 s
   gate FRELON multicorps     173 s        ~130 s
   invariances (dim 30)       4,8e-13 m    1,2e-14     (chemin inchangé + expm exact)
   S2 contre MBDyn            régime 3161,4 / 3162,0 · Q 51,74 / 51,27 — inchangés
```

Garde exécutable : `statique` asserte que Princeton consomme moins de 0,1
jacobien par pas. Non couvert, nommé : sur les petits modèles (Kapitza, 18
inconnues, 35 µs/pas) le temps est dans `phi_g` et la mécanique du pas, pas
dans le jacobien ; et sur S2, 10 % du mur est l'interpréteur Python de la
télémétrie par pas, hors noyau.

Rebasé sur `ebe63c9` (PR #1, paniques PyO3 et contacts manuels) : `bancs
rapide` y rougit sur `contact.tas` (« Newton ne converge pas à t = 0,0634 »)
**avec ou sans ce commit** — mesuré sur l'amont seul, résidu 5,65e-9. Il
n'est pas de ce chantier ; il est nommé ici pour qu'on ne le lui impute pas.

## LA RÉSOLUTION, PIÈCE PAR PIÈCE — et vinkulum passe en 3.14 (7 septembre, matin, suite)

Le chrono grossier disait « résolution 50 % du pas S2 » ; il fallait le
couper avant de toucher. Trois compteurs fins (`VINKULUM_CHRONO=1`) :

```
   tête S2, 60 000 pas : facto 3,16 s (105 µs pièce) · descente 2,13 s (7,2 µs)
                         · contrôle Jx−b 1,10 s (3,7 µs, À CHAQUE descente)
```

- **Le contrôle Jx − b juge la factorisation, pas le second membre** : une
  fois par factorisation suffit (la tête fait dix descentes par jacobien).
  1,10 → 0,12 s.
- **Factorisation assemblée directement dans la matrice faer** (plus de
  détour par nalgebra puis `from_fn`) : 3,26 → 2,88 s. Descente en place sur
  une copie de b : sans effet mesurable (2,12 → 2,10) — la copie n'était pas
  le coût, gardée parce que plus courte.
- **Le LU creux à 143 inconnues est PLUS LENT** que le dense : 138 µs contre
  105 (descente 3,2 µs contre 7,2, mais c'est la factorisation qui compte ; le
  hachage du motif n'y est pour rien, 0,29 s sur 4,4). `DENSE_MAX = 160` reste,
  mesuré. nnz 3 036 sur 20 449 : 15 % de remplissage, trop pour un creux.
- **Princeton, côté résidu** (53 000 évaluations, 19 µs) : `pose_pert` sur
  duaux à valeur nulle sans expm ni produit de duaux (R₀ + [δθ]×·R₀, 18
  produits par constante), et `m_mul_cst` pour les repères matériels
  constants (27 produits par constante au lieu de 27 produits de duaux —
  `Scalar::scale`). 1,99 → 1,74 s.

```
                              tour 1        tour 2
   tête S2 (vol 6 s)          0,219 ms/pas  0,192      MBDyn : 3161,4 / 3162,0 · 51,74 / 51,27 inchangés
   Princeton 8 s              1,98 s        1,74
   gate FRELON multicorps     129 s         123
   verification (3.14)        104 s         92
```

**Python 3.14** (consigne de Paul, 7 sept.) : la vérification et les bancs
de vinkulum tournent sur `.venv314` et `.venv314t` (roues cp314 / cp314t,
`maturin build --release -i .venv314/bin/python -i .venv314t/bin/python`,
`uv pip install --force-reinstall`). Princeton 1,74 s en 3.14, 1,73 en 3.14t
— le noyau est en Rust, l'interpréteur n'y pèse rien. La tête S2 reste en
3.13 : `vinkulum_s2` importe `build123d`, donc OCP. **La roue cp314 de
cadquery-ocp existe maintenant** (7.9.3.1.1, vérifié sur PyPI) ; ce qui
retient FRELON est son binding OCCT 8 compilé localement pour 3.13 — un
chantier à part, pas un obstacle.

Non couvert, nommé : `bancs rapide` reste rouge sur `contact.tas` (amont
`ebe63c9`, PR #1), donc les bancs après lui ne sont pas rejoués ici.

## PARTIR DE 2a_n − a_{n−1} — le prédicteur d'ordre 1, et le recul en pas fixe qu'il a fait sortir (7 septembre, matin, troisième tour)

Profil de la tête S2 devenu plat (aucun symbole au-dessus de 6 %) : les
micro-leviers sont épuisés, restait le NOMBRE de descentes, 4,94 par pas.
Il dépend du point de départ de Newton, qui était l'accélération du pas
précédent (ordre 0). **Extrapolation linéaire** (a et λ : 2a_n − a_{n−1}),
réinitialisée sur impact et si h change, réservée aux modèles ≥ 48
inconnues (même seuil que le jacobien hérité : Kapitza, 18 inconnues sous
une consigne en table à cassures, y perdait 6 %).

```
   tête S2   Newton 4,94 → 4,01 par pas · 0,192 → 0,173 ms/pas (MBDyn inchangé)
   Princeton 384 → 362 jacobiens, 1,74 s (stable)
```

**Trois bancs ont rougi, et chacun a donné une règle.**

1. *Quatre-barres flexible (h 4e-3)* : le départ extrapolé sort du bassin.
   **Repli** : si Newton lâche sur un départ extrapolé, le pas se rejoue à
   l'ordre 0 — l'échec n'est un échec qu'à l'ordre 0.
2. *KUKA iiwa sous gravité avec butées (h 2e-3)* : l'extrapolation à travers
   une butée change de branche, et le rejeu ne retrouve pas la trajectoire
   (chaotique). **On n'extrapole qu'après un pas LISSE** : ni rebroussement de
   Newton, ni pas hérité qui fait MONTER le résidu. Un pas hérité qui descend
   de moins de 10× n'est pas un accroc — le compter coûtait la moitié du gain
   (0,185 au lieu de 0,172).
3. *Le même KUKA, sondé* : à h 1,8 · 2,1 · 2,2e-3 l'ordre 0 SANS extrapolation
   lâchait aussi ; il ne passait qu'à 2,0e-3 exactement. **Le banc tenait par
   chance.** Ce n'est pas le prédicteur, c'est un Newton à pas fixe sans
   aucun recul sur un événement. **Recul en pas fixe** : un Newton qui lâche
   se rejoue au demi-pas jusqu'au point de grille manqué, puis h reprend
   (plancher h/64, sinon l'erreur d'origine). Sondé : 2,0 / 1,8 / 2,2e-3
   passent désormais, ordre 0 ou 1. Ce qui échouait passe, le reste est
   inchangé au bit.

Non couvert, nommé : le KUKA à 2,1e-3 lâche encore, même à h/64 — résidu de
départ 1e7 à chaque sous-pas, puis divergence : c'est la butée en pénalité
(couple `effort` sur `garde_deg` 1°) qui est impulsive, pas le pas. Et
l'écart avec/sans butées (0,438 m aujourd'hui, 0,494 ce matin, 0,677 entre
les deux) dit que cette chute de 3 s est CHAOTIQUE : le banc n'en juge que le
signe, c'est bien.

`VINKULUM_PRED0=1` rend l'ordre 0 partout (A/B, robustesse).
Régressions : cargo test/clippy/fmt, `verification` 92 s (3.14), gate FRELON
123 s, démo S2 régime 3161,4 / 3162,0 · Q 51,74 / 51,27.

## LE JACOBIEN TOURNE AVEC LES CORPS — la factorisation gardée se réutilise sans refaire J (7 septembre, quatrième tour)

Paul : « fais-le » — le chantier qui rend le jacobien moins périssable en
rotation. Trois voies, deux refusées à la mesure, la troisième sortie d'une
invariance relue dans le code.

**Refusé 1 — Broyden.** Quasi-Newton par sécantes de rang 1 sur le LU gardé
(Woodbury). jac/pas 0,47 → 0,30 mais Newton 4 → 5,3 et 0,242 ms/pas : la
dérive en rotation n'est pas une direction, c'est un sous-espace ; la
sécante d'un essai refusé (qui vise précisément la direction où J est faux)
n'en rattrape qu'une à la fois. Code retiré.

**Refusé 2 — paralléliser l'assemblage.** À 22 éléments (S2), rayon coûte
plus qu'il ne rend : jac 2,3 → 4,9 s à 8 fils, 10 s à 64. `PARALLELE_MIN = 48`
confirmé ; `VINKULUM_PAR_MIN` le règle pour mesurer.

**Retenu — J(q) ≈ T_r·J₀·T_cᵀ.** Les lignes de Φ sont MATÉRIELLES (d = uaᵀ·dw,
phi_r dans le repère de a) : elles ne tournent pas. Les colonnes (δr, δθ, ẇ
spatiaux) et les lignes dynamiques d'un corps tournent avec lui, Qᵢ = Rᵢ·Rᵢ⁰ᵀ
depuis la pose d'assemblage. Exact pour tout élément entre corps co-rotatifs
ou vers le bâti, approché sur une liaison mixte (tournant ↔ fixe : le bloc des
colonnes du fixe porte uaᵀ, qui a tourné). Et J⁻¹ = T_c·J₀⁻¹·T_rᵀ : **le LU
gardé sert tel quel**, on tourne b avant et x après (12 rotations 3×3, rien).
La règle des 10× reste seule juge de la péremption.

```
   S2   Newton 4,01 → 3,08 par pas         1ʳᵉ descente héritée : ×2 000 000 au lieu de ×1 600
   Princeton  362 → 178 jacobiens          (les nœuds de la poutre fléchie tournent aussi)
```

Mais jac/pas restait à 0,49 : la trace a montré le résidu dynamique refusé à
**3,3× du seuil** (2,47e-10 contre 7,5e-11) parce qu'il ne faisait plus 10×
sur les lignes de la liaison mixte — un jacobien encore bon, condamné à un
pas du but. **À moins de 100× du seuil, toute descente héritée est prise.**

```
                              tour 3        tour 4
   tête S2 (vol 6 s)          0,173 ms/pas  0,150      jac/pas 0,47 → 0,34 · MBDyn inchangé
   Princeton 8 s              1,74 s        1,65       jacobiens 362 → 156
   verification (3.14)        92 s          86
   gate FRELON multicorps     123 s         116
   KUKA à butées              2,0 / 1,9 / 1,8e-3 passent
```

Ce qui borne maintenant : le jacobien vit ~3 pas ; au-delà, c'est la liaison
mixte du plateau (couronne tournante ↔ bague fixe) dont le bloc dérive de 1,9°
par pas. Le pas suivant, s'il en faut un, est de rendre ces blocs exacts — ce
qui n'est PAS une rotation de lignes/colonnes (démontré : les colonnes de a
et celles de b exigeraient deux transformations de ligne incompatibles) mais
une ré-évaluation des 2–3 éléments mixtes, avec le LU à refaire. Ou le repère
tournant pour le rotor entier, où la liaison mixte devient un roulement
axisymétrique dont le bloc ne bouge pas.

## CE QUI VALAIT LE PLUS — contact.tas, la CI locale, et le jacobien mesuré bloc par bloc (7 septembre, soir)

Paul : « fais ce qui vaudrait le plus avant d'attaquer le maillage et la
robustesse ». Trois choses, dans l'ordre annoncé.

**`contact.tas` est vert — et ce qu'il jugeait avant était faux.** La PR #1
(6 sept.) a rendu les contacts manuels que l'appariement automatique effaçait
« au premier pas, sans un mot » : le tas de billes tombait SANS SOL, et le banc
passait. Avec le sol, c'est un vrai tas — Newton à pas fixe y lâchait à
t = 0,0634. Le recul au demi-pas (quatrième tour du matin) le tient : **25/25
bancs, `bancs OK` en 370 s.**

**La CI est locale.** Pas de crédit GitHub (décision Paul) : `ci/local.sh`
rejoue l'enchaînement de `ci/github-actions.yml` — fmt, clippy -D, test,
roue installée, `verification`, fraîcheur de l'API ; `--bancs` ajoute les
bancs et le contact. Le hook versionné `ci/hooks/pre-push` la lance à chaque
push (`git config core.hooksPath ci/hooks`). 1 min 42 sur 64 cœurs, et elle a
mordu à son premier passage : `docs/API.md` était périmé depuis la PR #1.

**Le jacobien, mesuré bloc par bloc.** Le repère tournant promis le matin
n'aurait rien donné de plus que le J tourné : les blocs entre corps
co-rotatifs sont déjà exacts, et un roulement entre un corps qui tourne et un
qui ne tourne pas garde un bloc qui dépend de l'angle relatif dans TOUT
repère. Restait à savoir ce qui tue vraiment J à 3 pas. Deux mesures :

- le taux de descente du J tourné ne bouge PAS avec son âge (1e-7 à l'âge 1,
  2e-8 à l'âge 11) — le rafraîchissement venait de la fin du pas, sur des
  lignes de Φ à l'arrondi (échelle 1/βh² = 4e8) que `r_tot` prenait pour de
  la péremption. **À moins de 100× du seuil, toute descente héritée est
  prise, un pas qui ne descend plus est le plancher, et `r_tot` ne juge plus.**
  S2 0,150 → **0,141 ms/pas**, jac/pas 0,34 → 0,29, Princeton 156 → 99
  jacobiens ;
- `VINKULUM_DIAG_J=1` compare T·J₀·Tᵀ au J frais au premier rafraîchissement
  après 3 s, bloc par bloc et case par case. Ce qui dérive : le pivot du
  PIGNON au bâti (34 % — Gᵀ vaut I à frais et Q tourné : ses lignes de corps
  tournent, pas ses colonnes de λ), le 6701ZZ (7 % et 3 %, mixte), le mât au
  bâti (7 %). Tout le reste sous 1e-4.

Essayé et **refusé, mesuré** : tourner aussi les lignes de Φ des pivots au
bâti comme leur corps b — dérive 6,2 → 4,1 %, jac/pas 0,29 → 0,26, mais
Newton 4,0 → 4,4 et **0,141 → 0,153 ms/pas**. Un J plus proche n'est pas un
Newton plus court ; la case du pignon qui reste (les lignes de rotation, 0,40
à frais contre 0,66 tourné) n'est pas une rotation.

**Le KUKA à butées est sorti du pile ou face.** À h = 2e-3 la chute avec
butées en pénalité passait ou lâchait selon le chemin de Newton (ordre 0 ou 1,
règles près du seuil), et à 1,8 / 2,1 / 2,2e-3 elle lâchait avec le solveur
d'origine ; de 1,0 à 1,9e-3 elle passe dans tous les modes essayés. Le banc
tourne à **1e-3** et le dit. Ce n'est pas un seuil qu'on déplace : c'est un
tirage qu'on cesse de jouer.

```
   tête S2   0,281 (matin) → 0,141 ms/pas (−50 %) · MBDyn 3161,4 / 3162,0 · 51,74 / 51,27 inchangés
   verification 82 s · gate FRELON multicorps ~100 s · bancs rapide 25/25 · CI locale 1 min 42
```

Ce qui reste nommé : la correction exacte des 2–3 blocs qui dérivent (pignon,
6701ZZ) par Woodbury sur le LU gardé — la seule voie qui ferait vivre J
dix pas ; le mécanisme est dessiné (lignes de l'élément comme sélecteur,
W₀ = J₀⁻¹P calculé une fois par factorisation), pas écrit.

## LE CCD SUR BOÎTE, CYLINDRE ET MAILLAGE — la capsule balayée (7 septembre, sur le ThinkPad)

Paul : « avance », puis « boîte et cylindre : capsule ». Le CCD ne jugeait
que sphère/sphère et sphère/plan ; boîte, cylindre et maillage en étaient
EXCLUS depuis le 5 (avant, ils tombaient dans la branche « plan » avec la
normale par défaut). Une bille de 5 mm à 50 m/s, h 1e-3, enjambe une paroi
de ±0,01 sans qu'aucun bout du pas ne soit à portée de la barrière (d̂ 1e-3) :
elle finit à +0,03, de l'autre côté, et le banc ne le voyait pas.

**Une idée, pas trois distances.** La sphère balayée sur le pas est une
CAPSULE [pa0, pa1] de rayon r, et elle traverse si la distance de son
segment à la primitive est ≤ r. La distance d'un point à un convexe
(triangle, boîte, cylindre fini) est convexe ; composée avec le mouvement
affine du centre elle est convexe en t, donc son minimum sur [0, 1] est
GLOBAL, et une section dorée (80 itérations, 0,618⁸⁰ ≈ 2e-17) le rend à
l'arrondi. Aucune distance segment/boîte ni segment/cylindre à écrire :
`min_convexe_segment` sur la distance point/primitive qu'on avait déjà
(`dist_boite`, `dist_cylindre`, `point_triangle`). Le maillage se parcourt
sous le BVH dilaté de r (slabs sur l'AABB — le cube contient la boule,
conservatif), triangle par triangle. Le segment est pris dans le repère du
PORTEUR, chaque bout à sa pose du moment : exact si le porteur ne tourne
pas sur le pas, approché sinon — la limite du plan porté, déclarée.

```
   bille 5 mm à 50 m/s, h 1e-3, paroi fixe (face à x = −0,01 ; maillage à 0)
                  sans CCD       avec CCD
   boîte          +0,0300        −0,0547      traversée → rebond
   cylindre       +0,0300        −0,0547
   maillage       +0,0300        −0,0361
```

Test Rust `ccd_capsule_balayee` : triangle traversé (0), parallèle à 0,3
(0,3 à 1e-12), coin de boîte (√2), flanc (0,5) et fond (0,25) de cylindre,
maillage touché / non touché / hors du triangle à x = 0. Banc
`contact.tir_paroi`, section 9b de `contact`.

Non couvert, nommé : la capsule (p1 ≠ p0) contre boîte, cylindre ou
maillage — le contact lui-même ne la traite pas (il lit p0) ; le porteur qui
tourne vite sur le pas ; et une sphère DÉJÀ dedans à t₀ — le CCD la voit
traverser, divise jusqu'à h/256 et lève, ce qui est juste pour une barrière
et faux pour une pénalité qui pénètre par construction (`ccd=True` est fait
pour IPC, comme avant). Coût : 84 évaluations de `point_triangle` par
triangle candidat — un maillage fin sous une sphère au repos le sentira, et
la distance segment/triangle en fermé est la relève, nommée dans le code.

Première construction sur le ThinkPad (16 cœurs, 14 Go) : `cargo build`
1 min 14, roue maturin en 3.14, `cargo test` 12/12. Le banc `maillages` a
rougi UNE fois sur l'exposant du BVH (0,83 : 43 µs à 32 768 triangles)
pendant que clippy compilait à côté ; seul, 0,37 deux fois. C'est une
mesure de temps sur un portable, pas le noyau — et c'est la deuxième fois
qu'un chrono juge une machine au lieu du code (§ Princeton).

## LE SECOND SONDAGE — trois blocages dans la SVD, et treize silences (7 septembre, suite)

Paul : « continue ». Après le CCD, la robustesse — au sens du 3 septembre :
un modèle absurde est REFUSÉ, jamais accepté en silence, jamais bloqué.
Quarante entrées absurdes jetées sur `contact`, `maillage`, `corps`, `Noyau`
et `simule`, chaque cas d'exécution dans son propre processus sous
`timeout` — une alarme Python ne peut pas interrompre une boucle Rust.

**Trois BLOCAGES, une seule cause.** Masse NaN, gravité NaN, et une normale
de plan NULLE : `simule` ne rendait jamais la main, aucune trace, un
processus à tuer — la classe de `h = 0`. La pile d'appels (gdb sur le
processus figé) : `nalgebra::SVD::new` ← `Modele::lstsq` ← `acc_init`.
`Matrix::svd` de nalgebra itère SANS plafond, et sur un NaN il ne converge
jamais. La normale nulle y arrive par `normalize()` (0/0 = NaN) ; la masse
et la gravité, par le système de l'accélération initiale. Racine, pas
symptôme : `svd_sure` refuse une matrice non finie avant d'itérer et
plafonne `try_svd` — sur les SIX sites SVD du noyau (`lstsq`, assemblage,
projection de vitesse, masque de redondance, statique ×2), et `acc_init`
refuse une accélération non finie. Puis les gardes à la déclaration, pour
que le refus dise QUOI : masse `!(m > 0)` (le `m <= 0` laissait passer NaN),
gravité finie, normale non nulle et finie.

**Treize acceptations silencieuses**, dont une qui simulait FAUX : un
contact `maille=` avec `b` un autre corps que le porteur du maillage lisait
le point le plus proche dans le mauvais repère — la bille tombait à travers
le sol, z = −0,93 à 0,5 s. Désormais `b` absent = le corps du maillage, `b`
autre = refus. Les autres : `a == b`, rayon négatif, `k = 0` (contact
fantôme), `d_hat = 0` sous barrière (ln 0), `v_eps = 0` avec frottement
(tanh(v/0)), deux primitives sur le second corps (la priorité tranchait en
silence), restitution > 1, triangle dégénéré (sommets répétés ou aire nulle
— `point_triangle` y divise par zéro), sommet NaN, et `simule` à fin NaN ou
passée, ou `tous = 0` : une trajectoire VIDE rendue sans un mot.

Ce que le sondage a jugé sain, et qui reste : la barrière IPC clampe d à
1e-4·d̂ (une pénétration initiale ne fait pas NaN, elle repousse) ; un
sphère déjà dedans sous CCD lève « traversée même au pas minimal », ce qui
est juste ; une normale non unitaire est normalisée ; les indices hors
bornes lèvent `IndexError` partout.

`verification.robustesse` : 14 → **30 modèles absurdes, 30 refus typés**,
dans la gate. Les seize nouveaux sont écrits avec ce qu'ils faisaient avant
(bloquait / traversait / liste vide), pour qu'on sache ce que chaque garde
tient.

Méthode, pour la prochaine fois : `timeout 20` par cas et par processus, et
`gdb -batch -ex run -ex bt` sous `timeout --foreground -s INT` pour lire où
ça boucle — le débogage du Newton (`VINKULUM_DEBUG`) ne disait rien parce
que le blocage était AVANT le premier pas.

## LA CAPSULE CONTRE LES TROIS PRIMITIVES, ET LE SEGMENT/TRIANGLE EN FORME FERMÉE (7 septembre, suite)

Paul : « fais la capsule contre les trois primitives et la distance
segment/triangle fermée ». Une tige contre un carter n'était qu'une bille à
son bout : le contact lisait p0 pour boîte, cylindre et maillage.

**Le point de l'axe est figé sur le pas, comme le point du maillage.**
L'abscisse s ∈ [0, 1] du point de l'axe [p0, p1] le plus proche de la
primitive est calculée en fin de pas avec `q_maille` — minimum exact de la
distance convexe le long de l'axe pour la boîte et le cylindre
(`min_convexe_segment` rend désormais l'abscisse), distance segment/triangle
sous le BVH pour le maillage — et le résidu AD voit la sphère de l'axe qui
s'y trouve (`s_axe`, nul partout ailleurs : p0 au bit près). Exact aux bouts ;
à l'intérieur, la distance est stationnaire en s au minimum, donc l'erreur
d'un s figé est d'ordre 2 en ce que le pas fait bouger. Pas de golden
section en duaux, pas de nouvelle branche AD : les trois branches existantes
servent telles quelles.

**Segment/triangle en forme fermée** (Ericson 5.1.10) : le minimum est une
intersection (0), ou entre le segment et une des trois arêtes (`seg_seg`,
5.1.9, dégénérés et parallèles compris), ou entre un bout et le triangle
(`point_triangle`). `Maillage::plus_proche_segment` parcourt le BVH avec le
slab dilaté de la meilleure distance courante et rend (d, s, q) ; le CCD s'y
réduit (`traverse` = d ≤ r), et les 84 évaluations par triangle candidat de
la section dorée sont parties avec leur commentaire `ponytail:`.

```
   contact.capsules_primitives — 12 cas contre la géométrie exacte (scipy, Brent borné)
   boîte     face ∥ +0,020 · bout +0,030 · au-dessus d'un coin −0,050 · le long d'une arête −0,200
   cylindre  ∥ axe +0,020 · ⊥ −0,050 · bout sur le fond +0,020 · par-dessus l'arête 0,000
   maillage  ∥ plan +0,020 · bout +0,030 · traverse +0,050 · au-delà du bord −0,1736
   pire écart 1,0e-12 · tige à plat qui tombe : boîte et maillage z_min +0,04502, rebond +0,07278, identiques
```

Test Rust : `seg_tri` traversée (0, s ½), parallèle (0,3), au-delà d'une
arête (0,5 — et le minimum est PLAT : tout s de [⅓, ⅔] est juste, le test
l'avait d'abord exigé à ½), segment dégénéré ; `plus_proche_segment` = `seg_tri`
sur un triangle seul.

Non couvert, nommé : la capsule contre le PLAN (porté ou fixe) reste une
bille à p0 — le point le plus proche y est un bout, et il saute d'un bout à
l'autre quand la tige bascule ; le CCD d'une capsule contre ces primitives
balaie trois sphères (les deux bouts, le point figé), une tige fine peut
enjamber un petit obstacle entre eux ; et le premier pas d'une simulation
lit s = 0 et q = origine, comme le maillage depuis le 5 sept.

## Robustesse de l'API et validation sans sélection du meilleur accord — 6 septembre 2026

L'analyse de `a2bf879` a reproduit une chute libre faussée d'environ 49 %
lorsque le pas demandé dépassait la durée restante, des paniques sur des lois
et indices mal formés, une réflexion acceptée comme rotation, un temps NaN
restauré sans erreur et un rayon de sphère négatif accepté. Les corrections,
leurs cas de reproduction et les changements de contrat sont détaillés dans
la [note d'intégration](ROBUSTESSE_2026-09-06.md).

La limitation de l'accélération initiale utilise désormais le pas réduit à
la durée disponible. Les arguments mal formés sont refusés à la déclaration.
`pose_etat` valide l'état entier avant de le modifier, tout en conservant la
ré-orthonormalisation des perturbations valides. La suite de robustesse passe
de 30 à **48 refus typés** et vérifie aussi les cas valides, la chute libre
analytique et l'absence de restauration partielle.

La CI locale utilise le venv actif ou celui du dépôt, avec `PY` et `MATURIN`
surchargeables, et affiche les diagnostics complets. Son hook est activé dans
ce clone ; cette configuration doit être refaite après un nouveau clone.

Le score Maryland choisissait après calcul le `N_crit` donnant le meilleur
accord avec la mesure. Chaque hypothèse est maintenant publiée séparément.
La campagne complète a été rejouée sous Python 3.13.12 en **1 959 s**, avant
le relèvement du minimum à 3.14, sans erreur d'exécution :
**44/62 lignes dans leur tolérance**, au lieu de 44/54 dans l'ancienne
présentation. Les huit nouvelles lignes rendent visibles les hypothèses
écartées auparavant ; les tolérances et les paramètres physiques n'ont pas
été ajustés. `verification.validation_scores` garde cette propriété.

Contrôles : 13 tests Rust, formatage, Clippy et CI complète réussis ;
vérification Python en 76 s et API générée de 178 entrées à jour. La dernière
garde sur les sphères a été vérifiée dans une roue reconstruite, avec les
48 cas de robustesse, puis installée. Les écarts physiques restent ouverts
dans [le rapport généré](VALIDATION.md).

**Consigne de compatibilité : Python 3.14 ou plus.** Le minimum est maintenant
imposé par les métadonnées du paquet, l'import et les scripts de CI ;
`.python-version` sélectionne 3.14 et le workflow préparé cible 3.14/3.14t.
Le venv local a été recréé en **3.14.7**. La CI complète y passe en **71 s**
pour la vérification Python, avec 13 tests Rust et 48 refus typés. Les cas
capsule/primitives et CCD/parois ont aussi été rejoués, ainsi qu'un calcul de
polaire NeuralFoil. La roue `cp314` exige `Requires-Python: >=3.14`.
Les anciens chiffres sous 3.13 restent des traces historiques ; la campagne
physique complète de 1 959 s n'a pas été rejouée sous 3.14. La
[note d'intégration](ROBUSTESSE_2026-09-06.md) distingue ces périmètres.
