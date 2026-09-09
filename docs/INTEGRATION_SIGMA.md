# Cinématique σ : première application de l'étude scientifique

La version 0.8.0 ajoute `sigma_lie` à `Noyau.simule` et
`Noyau.audit_jacobien`. **Zéro reste le défaut.** Cette option expérimentale
vise l'erreur géométrique des rotations spatiales dont l'axe varie pendant
un pas. Elle ne remplace ni la formulation des corps flexibles ni le
solveur de contact.

Cette réalisation applique une piste de l'[étude scientifique 2001–2026](ETUDE_SCIENTIFIQUE_VINKULUM_2001_2026.pdf).
Sa source principale est Holzinger, Arnold et Gerstmayr,
*σ-modified Lie Group Generalized-α Methods for Constrained Multibody Systems*,
Mechanism and Machine Theory 217 (2025), 106236,
[DOI : 10.1016/j.mechmachtheory.2025.106236](https://doi.org/10.1016/j.mechmachtheory.2025.106236).
Le code résulte d'une dérivation des équations cinématiques (12)–(16),
avec les conventions de Vinkulum. Aucun code de solveur concurrent n'a
été utilisé.

## Résultats et décision

**252 essais chronométrés sur 252 terminent**, sans réduction du pas. Le
tableau compare le plus grand pas de chaque cas ; les erreurs portent sur
tous les pas du rejeu de validation. Les temps sont les médianes de trois
répétitions sur CPU 0, avec un fil demandé par bibliothèque.

| Cas | Pas (s) | Erreur σ=0 | Erreur σ*=0.665 | Coût σ*/σ=0 au même pas |
|---|---:|---:|---:|---:|
| Solide sphérique, vitesse affine | 0.04 | 1.352e-4 rad | 7.744e-9 rad | 1.42× |
| Solide libre anisotrope | 0.04 | 9.778e-3 rad | 8.133e-3 rad | 1.49× |
| Pendule plan | 0.02 | 4.593e-4 rad | 4.593e-4 rad | 1.09× |
| Toupie rapide, 300 rad/s | 0.001 | 4.028e-3 rad | 6.150e-4 rad | 1.27× |
| Deux axes commandés | 0.02 | 2.000e-3 rad/s | 8.040e-4 rad/s | 1.15× |
| Branche flexible | 0.002 | 2.76029e-4 rad | 2.76030e-4 rad | 1.07× |
| Boucle spatiale avec Cardan | 0.004 | 1.56657e-2 rad | 1.59007e-2 rad | 1.03× |

Le solide sphérique isole précisément le terme géométrique visé : ce très
grand gain ne se transpose pas aux autres problèmes. Sur la toupie,
l'erreur d'orientation baisse de **6.55×**, mais celle de la norme de
réaction passe de **0.01883 N à 0.03944 N**. Le réglage σ=1 la porte à
0.06154 N. Le gain dépend donc aussi de la grandeur observée.

Le [bilan calculable](bancs/integration-sigma-bilan-2026.json) retient les
meilleurs temps mesurés satisfaisant différents seuils. Ces seuils servent
à lire les résultats après exploration ; ce ne sont ni des critères
préenregistrés de supériorité, ni une optimisation continue du pas.

| Objectif sur la grille mesurée | σ=0 | σ*=0.665 | Conséquence |
|---|---:|---:|---|
| Toupie : orientation ≤ 1 mrad | 4.585 ms à h=0.0005 | 3.012 ms à h=0.001 | σ* prend 34 % de temps en moins |
| Même toupie : orientation ≤ 1 mrad **et** réaction ≤ 5 mN | 4.585 ms à h=0.0005 | 10.655 ms à h=0.00025 | σ* prend 2.32× plus de temps |
| Axes commandés : vitesse angulaire ≤ 1e-4 rad/s | 8.361 ms à h=0.0025 | 4.954 ms à h=0.005 | σ* prend 41 % de temps en moins |
| Branche flexible : position ≤ 1 µm | 74.258 ms à h=0.00025 | 77.327 ms au même pas | Aucun gain |
| Boucle spatiale : position ≤ 100 µm | 143.922 ms à h=0.001 | 146.928 ms au même pas | Aucun gain |

Les coûts inférieurs à la milliseconde sont sensibles au bruit du système ;
les minimums et maximums sont conservés, sans suppression d'échantillons.
L'écart entre les deux références internes reste sous 2 % de la plus petite
erreur mesurée, dans les quatre grandeurs cinématiques. Le contrôle par
resserrement de tolérance des références EDO est également archivé.

**Décision : livrer comme option expérimentale, garder zéro par défaut.**
L'intégration géométrique apporte un gain réel et ciblé. Elle ne résout
pas le retard de coût des poutres ni la couverture multiphysique. Le banc
ne contient aucune nouvelle exécution MBDyn ou Simpack : les résultats
externes restent ceux de la [confrontation 0.7.2](CONFRONTATION_MBDYN_0.7.2.md).

Les appels sans le nouvel argument reproduisent **bit à bit**, y compris
les statistiques Newton, les sept trajectoires comparées à la roue 0.7.2
gelée ([contrôle des deux binaires](bancs/integration-sigma-zero-2026.json)).
Cela vérifie ce corpus, pas une identité garantie sur tout modèle.

## Utilisation et domaine

```python
rho = 0.9
gamma = 0.5 + (1-rho)/(1+rho)
beta = 0.25*(gamma+0.5)**2
sigma = gamma/(3*beta)  # 0.665 pour rho=0.9
trajectoire = N.simule(1.0, 0.001, rho=rho, sigma_lie=sigma)
blocs, coefficients = N.audit_jacobien(0.001, rho=rho, sigma_lie=sigma)
```

`sigma_lie` doit être fini et compris entre 0 et 1 ; `rho` conserve son
domaine `[0,1[`. Chaque appel Python fixe explicitement son réglage, même
après un appel utilisant une valeur différente. L'audit restaure l'état
physique et le réglage σ qu'il a trouvés. Côté Rust, le champ est
`Modele.sigma_lie`, initialisé à zéro par `Modele::new`.

La résolution locale utilise la carte `||θ|| < π`. Un essai hors carte
rend un résidu non admissible, sans modifier le modèle. Le mécanisme
existant de réduction du pas peut le rejouer ; l'échec au pas minimal
rend une erreur et restaure le dernier état accepté. Ce domaine concerne
l'incrément d'un pas, pas l'orientation cumulée du solide.

Le couplage GGL possède sa propre dérivée de configuration et est testé.
`sigma_lie != 0` avec `adaptatif` est refusé : l'estimateur de demi-pas
historique n'est pas encore validé pour cette cinématique. Le multi-rythme
Python conserve σ=0 ; le multi-rythme Rust refuse un champ σ non nul.
`PontNoyau`, le pont adjoint temporel, utilise exclusivement σ=0.
`simule_em` reste un intégrateur distinct. Les tests des tangentes incluent
les contacts, mais le gain d'ordre des trajectoires lisses ne constitue
pas une garantie aux impacts ni un gain de précision du CCD.

## Équations effectivement résolues

Les vitesses angulaires de Vinkulum sont spatiales, avec composition à
gauche :

\[
R_1=\exp(\widehat\theta)R_0,\qquad
\omega_1=J_l(\theta)\dot\theta_1.
\]

Soit `b` l'incrément Newmark classique, incluant la correction GGL lorsqu'elle
est active. La vitesse et l'accélération algorithmique gardent leurs
équations generalized-alpha. Seule la configuration change :

\[
\theta=b+sD(\theta)\omega_1,\qquad
s=\sigma h\beta/\gamma,\qquad D(\theta)=J_l(\theta)^{-1}-I.
\]

Le code évalue directement

\[
D(\theta)\omega=-\tfrac12\theta\times\omega+
\frac{1-(t/2)\cot(t/2)}{t^2}\,
\theta\times(\theta\times\omega),\quad t=\|\theta\|.
\]

Une série paire évite les soustractions imprécises près de zéro. Le Newton
local porte sur `δ=θ-b`, avec vérification du résidu et recherche de
descente. Il est borné à douze itérations, chacune avec au plus dix essais.
Il élimine ces trois inconnues localement, sans agrandir le système global.

Avec

\[
A=I-s\,\partial_\theta[D(\theta)\omega_1],\qquad
B=I+\sigma D(\theta),
\]

la dérivée spatiale de configuration vis-à-vis de l'accélération physique
est `β h² c J_l A⁻¹ B`, où `c=(1-α_f)/(1-α_m)`. Celle d'une variation
indépendante de `b`, notamment GGL, utilise `J_l A⁻¹`. Le Jacobien global
conserve cette distinction, ainsi que les dérivées physiques déjà présentes.

**Ce n'est pas une reproduction du correcteur approximé de l'article.**
Vinkulum reconstruit l'état à partir de chaque essai du Newton global,
indépendamment du chemin parcouru par les itérations. Il résout donc la
relation implicite ci-dessus et garde la dérivée de `D`. La simplification
« sans tangente » du correcteur publié pour σ=1 ne s'applique pas directement :
ici `J_l A⁻¹ B` n'est généralement pas l'identité. Aucun gain de temps de
l'article n'est repris comme résultat de Vinkulum.

Pour une vitesse spatiale affine, le développement de Magnus contient
`−h³(ω₀×α₀)/12`. Le choix `σ*=γ/(3β)` restitue ce terme. Il annule une
composante de l'erreur locale ; il ne rend pas tous les problèmes d'ordre
supérieur. L'initialisation algorithmique et les autres sources d'erreur
restent déterminantes.

## Protocole reproductible

Le [programme de mesure](../ci/mesure_sigma.py) compare σ=0, σ*=0.665 et
σ=1 dans **le même binaire Rust**, à ρ=0.9, tolérance Newton `1e-12` et
plafond de 25 itérations. Sept cas, quatre pas et trois répétitions donnent
252 calculs mesurés, après une chauffe de chaque variante. L'ordre des
variantes est inversé une répétition sur deux. Les sorties ont la même
cadence physique pour les quatre pas d'un cas. La construction du modèle,
les références et le calcul des erreurs sont hors du temps d'intégration.

Chaque calcul chronométré est rejoué hors du chronométrage conservé, avec
une sortie à chaque pas : 252 rejouements de validation, dont l'état final
doit être identique à celui du calcul chronométré. Les erreurs retenues
sont des maxima sur **tous ces pas**, y compris le transitoire initial :
position, distance géodésique entre orientations, vitesse linéaire et angulaire.
Le pendule et la toupie ajoutent la norme de la réaction au pivot. Ces
maxima échantillonnés ne certifient pas les extrema entre pas. L'archive
garde aussi les erreurs sur la grille de sortie commune aux chronométrages.

Les quatre modèles à un solide utilisent les équations d'Euler dans le
repère matériel, intégrées avec DOP853 à deux tolérances. Le solide sphérique
sous couple constant possède en outre une vitesse spatiale affine exacte.
Les deux axes commandés utilisent des rotations analytiques. La branche
flexible et la boucle spatiale utilisent une référence **interne** affinée,
contrôlée par division du pas et comparaison de deux valeurs de σ ; elle
ne prouve pas à elle seule la justesse du modèle physique.

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python ci/mesure_sigma.py --sortie /tmp/integration-sigma.json
```

La campagne fixe l'affinité au premier CPU autorisé, journalise le binaire,
le programme, Python, les fils demandés et les statistiques Newton. Les
[données brutes](bancs/integration-sigma-2026.json) conservent les échecs
éventuels et les dispersions ; elles sont liées à leur empreinte de binaire.

## Contrôles de correction

Quatre tests Rust couvrent l'inverse géométrique, les petites rotations,
le signe du terme de Magnus, la covariance et les tangentes implicites
contre des perturbations sur SO(3). Quatre tests Python couvrent les
équations indépendantes à trois valeurs de ρ, les options invalides et
leur atomicité, le défaut zéro, le mouvement plan, une trajectoire GGL,
le rejeu d'un pas hors carte et la restauration après échec.

L'audit global couvre douze familles d'éléments, dont les poutres, les
contacts, l'aérodynamique, les contraintes non holonomes et les commandes
rhéonomes. Il est exécuté à σ=0, 0.6 et 1, avec les variantes GGL existantes.
Les limites de cet audit restent celles de son protocole : les blocs trop
petits face au maximum global sont filtrés, et la raideur géométrique du
contact non lisse n'est pas entièrement incluse dans le Jacobien.
