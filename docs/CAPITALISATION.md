# Ce que le métier a déjà payé — et que nous n'avons pas à repayer

Abaqus n'a pas inventé la mécanique : il a rétro-conçu ce que Nastran avait
appris, et Nastran ce que la littérature avait établi. Chaque génération
capitalise sur les erreurs publiées de la précédente. **Leur expérience est
publique ; seul leur code ne l'est pas.** C'est l'asymétrie qui rend
l'objectif atteignable : on ne rattrape pas 40 ans de développement, mais on
part au niveau où 40 ans de publications nous déposent.

Ce fichier liste les pièges DÉJÀ PAYÉS par d'autres, avec le remède et l'état
dans vinkulum. Il n'est utile que s'il reste opérationnel : chaque ligne dit
quoi faire, pas seulement quoi savoir.

## Intégration en temps

| piège payé | par qui | remède | état |
|---|---|---|---|
| Le trapèze de Newmark est **inconditionnellement instable** avec contraintes | Cardona–Géradin 1989 | dissipation haute fréquence obligatoire (HHT-α, α-généralisé) | **acquis** — ρ∞ = 1 refusé, retrouvé par la mesure |
| L'index 3 direct **dérive** si Φ n'est pas tenu à chaque pas | Baumgarte, puis toute la littérature DAE | résoudre Φ = 0 par Newton, pas stabiliser | **acquis** — Φ à 1e-15 |
| L'index 3 laisse Φ̇ violé à O(h²) | Gear–Gupta–Leimkuhler 1985 | GGL / index 2 stabilisé | **mesuré sans conséquence** : Φ̇ ne DÉRIVE pas (même ordre à 2 et 100 périodes) et l'énergie tient à 9e-6 |
| La projection de vitesse annule Φ̇ mais DISSIPE | classique (et le tableau de Paul le dit) | ne pas la prendre pour gratuite | **payé par nous** : ordre 1,92 → 0,76 et ΔE/E ×4 500 |
| Les angles d'Euler dégénèrent aux grandes rotations | tout le domaine depuis les années 80 | rester sur SO(3), mettre à jour sur le groupe | **acquis** |
| Ordre 2 en q ET en λ n'est pas automatique | Arnold–Brüls 2007 | schéma Lie-groupe index 3, preuve publiée | **acquis** — 2,00 mesuré sur les trois |
| Le critère de Newton mis à l'échelle sur les forces EXTÉRIEURES échoue là où la dynamique est inertielle | — | y ajouter ‖M·u̇‖ + ‖ω × J_s ω‖ | **payé par nous** le 3 sept. (rotor : `forces()` = 0) |
| Linéariser en repère FIXE autour d'un état tournant | Coleman–Feingold, puis tout le rotorcraft | coordonnées multipales (HOST) ou Floquet exact | **acquis** — `vinkulum.floquet` |
| L'exposant de Floquet est ambigu modulo 2π/T, et ζ n'est pas invariant par repliement | classique | publier σ (non ambigu) ; ζ demande la fréquence propre | **payé par nous** — ζ 0,427 au lieu de 0,050 |

## Éléments et formulation

| piège payé | par qui | remède | état |
|---|---|---|---|
| **Verrouillage en cisaillement** des poutres | des centaines d'articles depuis 1970 | intégration réduite (un point) ou formulation mixte | **acquis** — un point d'intégration |
| Réduction modale linéaire : rate la raideur géométrique | littérature rotorcraft (CAMRAD, Dymore) | poutres **géométriquement exactes** pour les pales, Craig–Bampton pour le raide | **acquis** (poutre GE), CB à faire |
| La tangente d'une poutre GE **n'est pas symétrique** hors équilibre | Simo–Vu Quoc 1986 | ne pas symétriser dans Newton, ne pas s'en alarmer | **acquis** — écrit dans le test |
| Contraintes **redondantes** → système singulier | mécanismes plans en 3D, Bricard | détection de rang, ou moindres carrés | **acquis** — QR pivoté au départ |

## Aéromécanique du rotor

| piège payé | par qui | remède | état |
|---|---|---|---|
| Lifting line 2D + tables suffit pour le **trim**, pas pour les charges | CAMRAD II, HOST, RCAS | garder la LL pour le trim, coupler CFD ensuite | LL acquise, trim à faire |
| Le trim se fait par **jacobien commandes→sorties + Newton** | CAMRAD II (§2.1 SABRE D1.6) | même méthode, mais jacobien par **AD** au lieu de DF | à faire — brique n° 1 |
| Le couplage CFD/CSD **faible** (une fois par tour) suffit et se trime par construction | Potsdam 2004, Altmikus 2002, SABRE 2018 | delta-method ΔF, historique cumulé | à faire, si CFD un jour |
| Le sillage domine les charges en vol d'avancement | HART-II, UH-60A Airloads | inflow dynamique puis sillage prescrit | inflow uniforme seulement |
| Linéariser en **repère fixe autour d'un état tournant** ne marche pas | Floquet, littérature rotor | repère tournant ou théorie de Floquet | **payé par nous** — refus explicite, mesuré |

## Méthode et outillage

| piège payé | par qui | remède | état |
|---|---|---|---|
| Une **différence finie** ne vérifie pas une composante nulle d'un champ à plancher d'arrondi | classique en AD | seuil portant le bruit plancher/2ε | **payé par nous** le 3 sept. |
| Les sensibilités par DF coûtent n+1 simulations et rendent 3 chiffres | tous les codes du métier | AD, mode avant puis inverse | avant acquis, inverse à faire |
| Une matrice d'influence de trim FIGÉE converge linéairement ; la refaire coûte n évaluations | ABAQUS §2.2.2 (BFGS à noyau), HOST | mise à jour de rang 1 entre deux réfections | **acquis** — Broyden, ×1,3 à ×1,6 mesuré |
| Reprendre **BFGS** sur un jacobien non symétrique | classique | Broyden pour le non symétrique, BFGS pour le SPD | **acquis** — écrit dans `trim.py` |
| Un solveur qui **converge** ne prouve pas qu'il a simulé le bon modèle | leçon n1 MBDyn (FRELON, 31 août) | chaque grandeur publiée est relue en sortie | **acquis** |
| Deux codes qui s'accordent sur un modèle **faux** ne valent rien | FRELON, 4 sept. (9 mm d'envergure) | contre-solveur indépendant, et divergence = signal | **acquis** |

## Ce que le manuel de théorie ABAQUS ajoute (lu le 3 septembre)

Paul a fourni le *ABAQUS Theory Manual* entier. C'est de la THÉORIE — des
équations publiées, pas du code — donc lisible sans toucher au clean-room qui
nous interdit les sources de MBDyn. Quatre choses en sortent, dont une qui
valide notre axe décisif **par écrit du concurrent**.

### 1. L'axe A est confirmé par le manuel d'Abaqus lui-même

`ABAQUS/Design` (§2.17.1) est leur module de sensibilités. Le manuel écrit ce
qu'il couvre, et c'est plus étroit qu'on ne l'imaginait :

> « ABAQUS/Design supports design sensitivity analysis (DSA) for
> **nonperturbation, static stress problems** […]. In the current capability
> **only solid elements with elastic or hyperelastic properties** may be made
> design dependent. »

Ni dynamique, ni poutres, ni contact frottant, ni rotor. Et la méthode est
**semi-analytique** : les vecteurs élémentaires sont obtenus par différences
finies, seule la résolution est analytique. Le manuel en nomme lui-même le
prix, dans les mots exacts qu'on a mesurés ici le matin même :

> « If the interval is too small, round-off or **cancellation errors** occur
> due to loss of precision during the differencing operations. On the other
> hand, if the interval is too large, **truncation errors** may occur. ABAQUS
> will automatically choose a perturbation size that provides the best
> compromise. »

**Un compromis que l'AD n'a pas à faire** — il n'y a ni annulation ni
troncature dans un nombre dual. Notre axe A n'est donc pas une prétention
contre un concurrent supposé faible : c'est un fait que son propre manuel
écrit. Et leur DSA est en mode DIRECT (une résolution par paramètre), donc son
coût croît avec le nombre de paramètres — ce que le mode inverse supprime.

### 2. Le pas adaptatif : recette suivie, mesurée, et REFUSÉE PAR DÉFAUT

> **Fait le 3 septembre, puis mesuré. Le verdict n'est pas celui qu'on
> attendait, et c'est pour ça qu'on le garde.** À précision ÉGALE, le contrôle
> de pas par résidu de demi-pas coûte **×2,8** le pas fixe sur un problème
> oscillatoire à deux échelles (banc `adaptatif`, gardé par assert). Il est donc
> implémenté, opt-in, et **désactivé par défaut** — comme chez ABAQUS, qui
> offre les deux.
>
> La raison n'est pas un défaut d'implémentation, c'est le DOMAINE de
> l'estimateur : le résidu de demi-pas mesure le déséquilibre LOCAL, et sur un
> mécanisme peu dissipatif l'erreur est dominée par la PHASE accumulée, qu'il
> ne voit pas. Le manuel le dit lui-même : l'algorithme est « purely
> empirical » et économique « in initially excited problems with **high
> dissipation**, such as impulsively loaded problems, with extensive
> plasticity ». Une voilure tournante n'est ni l'un ni l'autre.
>
> Ce que la brique apporte quand même : le pas se REJOUE quand Newton ne
> converge pas, au lieu d'échouer. Et le jour où un cas dissipatif entre au
> banc — contact, impact — elle est là, avec sa mesure à côté.

La recette, pour mémoire :

C'est le **résidu de demi-pas** (Hibbitt & Karlsson 1979, §2.4.1), et l'idée
tient en une phrase :

> « Satisfaction […] at the end of each time step ensures equilibrium at these
> points in time but **does not say anything about the quality of equilibrium
> at intermediate time points**. »

L'accélération étant supposée linéaire sur le pas (c'est la base de Newmark),
l'état à mi-pas s'intègre exactement — avec ā(τ) = a₀ + (τ/h)(a₁ − a₀) :

```
   ü(h/2) = (a₀ + a₁)/2
   u̇(h/2) = u̇₀ + (h/8)(3a₀ + a₁)
   Δu(h/2) = (h/2)u̇₀ + h²(5a₀ + a₁)/48
```

On évalue le résidu d'équilibre à cet état, on prend son plus grand terme, et
on le compare à une force typique P du problème. Le manuel donne même
l'échelle, ce qui évite de la calibrer à l'aveugle :

| R(t+h/2) | qualité |
|---|---|
| ≈ 0,1 P | haute précision |
| ≈ P | précision moyenne |
| ≈ 10 P | grossier |

Coût : une évaluation de résidu par pas, soit une fraction d'une itération de
Newton. **Non fait** — c'est le chantier suivant côté intégrateur, et il est
désormais spécifié ligne à ligne au lieu d'être un vœu.

Deux mises en garde du même paragraphe, à ne pas perdre : changer le pas
**injecte du bruit haute fréquence**, d'où la dissipation légère (α = −0,05
chez eux ; ρ∞ = 0,9 chez nous) qui devient nécessaire et non plus optionnelle ;
et l'algorithme de décision lui-même est « **purely empirical** ».

### 3. Ce que leurs superéléments ne savent pas faire, et qui nous concerne

Craig–Bampton est là (§2.14.1) : modes de contrainte statiques plus modes
normaux à interfaces bloquées. Mais le manuel pose une limite décisive pour
une voilure tournante :

> « the response within a substructure, once it has been reduced to a
> superelement, is considered to be a **linear perturbation about the state**
> of the substructure at the time it is made into a superelement. »

Un superélément Abaqus **ne tourne pas**. Pour une pale, la raideur
centrifuge et le raidissement géométrique dépendent du régime : réduire puis
faire tourner perd exactement ce qui compte. C'est pourquoi les codes
rotorcraft gardent des poutres géométriquement exactes là où un code de
structure réduirait — et c'est ce que vinkulum fait déjà.

### 4. Une règle de méthode qu'on appliquait sans l'avoir écrite

> « the intermediate, nonconverged solutions obtained during the iteration
> process are usually **not on the actual solution path**; thus, the
> integration of history-dependent variables must be performed completely
> over the increment at each iteration and **not obtained as the sum** of
> integrations associated with each Newton iteration. »

C'est la raison pour laquelle notre `Inflow` est un état explicite mis à jour
**après** le pas accepté, et jamais dans le résidu. On l'avait fait par
prudence ; c'est une règle, et elle a un nom.

## NASTRAN-95 : ce qu'il donne, et une CORRECTION DE LICENCE (3 septembre)

Paul a transmis `github.com/nasa/NASTRAN-95`. Deux constats, dont un qui
touche la doctrine.

**Ce n'est PAS le domaine public.** Le dépôt est sous **NASA Open Source
Agreement 1.3**, une licence de réciprocité : les œuvres dérivées se
redistribuent sous NOSA, avec obligation de documenter les modifications.
Nos notes internes disaient « NASA, domaine public » — c'était faux, et sur
une question de licence ça ne se laisse pas passer. **Le régime clean-room
s'applique à NASTRAN-95 exactement comme à MBDyn (GPL-2)** : on lit les
formulations, on ne copie pas une ligne. Les manuels, eux, sont de la
documentation NASA — les lire est du même ordre que lire la théorie ABAQUS.

**Et il s'arrête à Guyan.** Le Programmer's Manual (SP-223, 1972) décrit la
réduction par transformation de contraintes — `MCE2` partitionne [Kgg] et
forme [Knn] = [K̄nn] + [Gm]ᵀ[Kmn] + [K̄mn]ᵀ[Gm] + [Gm]ᵀ[K̄mm][Gm], la même
opération étant faite sur [Mgg], [Bgg] — c'est-à-dire la réduction STATIQUE,
appliquée aussi à la masse. **Craig–Bampton n'y est pas** : les modes
normaux à interfaces bloquées, qui sont l'autre moitié de la méthode, sont
arrivés dans les NASTRAN commerciaux plus tard. Pour cette brique, la source
utile reste le manuel ABAQUS §2.14.1, qui l'a en entier.

Ce qu'il apporte quand même : trois méthodes d'extraction de valeurs propres
(puissance inverse, déterminant, Givens) et une doctrine de partitionnement
des degrés de liberté par « sets » qui a survécu cinquante ans dans tous les
codes de structure.

## Ce que ce tableau dit du calendrier

Sur vingt et un pièges recensés, **seize sont déjà derrière nous** — non parce
qu'on est malin, mais parce qu'ils étaient publiés et qu'on les a lus. Les
quatre qui restent (sillage, adjoint en temps NON LINÉAIRE, contact, corps
flexibles 3D) sont des chantiers identifiés, pas des inconnues. Le trim, la
stabilité en rotation, Craig–Bampton, l'AD inverse et le raidissement
centrifuge en sont sortis le 3 septembre ; le pas adaptatif et GGL en sont
sortis autrement — **mesurés, et refusés sur mesure**.

**C'est exactement la thèse de cette capitalisation** : le temps de
développement d'un noyau moderne n'est pas le temps qu'ont mis les anciens.
Ce qui reste incompressible n'est pas l'algorithmique — c'est la
**validation** : les données d'essais, les cas industriels, la confrontation
au réel. C'est là, et seulement là, que 40 ans comptent encore.
