# Projection statique des contraintes sans équations normales

La statique assemble maintenant son gradient de contraintes en triplets
et calcule les réactions directement à partir de `Gᵀ`. Elle évite la
construction de `G` dense et de `GGᵀ`, y compris pendant la recherche du pas.
Les contributions locales et la factorisation du système de Newton
complètent l'[assemblage local précédent](ASSEMBLAGE_LOCAL.md).

Sur la chaîne articulée non linéaire à 256 corps, le temps médian passe de
**1,076 s à 0,01991 s**, soit **54,05×**, et le pic mémoire de **112,4 à
42,7 Mio**, avec le même travail de Newton et la même tolérance stricte.

Le témoin est Vinkulum au commit
`4b79a761c47d3de2d74018808258b38d9c5958db`, avec son extension conservée.
Les mesures comparent deux versions de Vinkulum. Aucun nouveau classement
face à MBDyn ou Simpack n'en découle.

## Méthode et contrat numérique

Le projecteur minimise `‖Gᵀλ − f‖₂`. Il décompose le graphe biparti entre
lignes de contraintes et degrés physiques en composantes indépendantes.
Les lignes nulles ont un multiplicateur nul. Les coefficients partagés
sont sommés dans l'ordre de dispersion ; seuls les zéros exacts sont omis.
Le champ reste celui de `phi_g()`, sans lui appliquer le masque d'activité
de la dynamique.

Les blocs scalaires se résolvent par division. Les petits blocs utilisent
un QR dense avec pivotage des colonnes. Au-delà de 32 multiplicateurs,
les blocs ayant au moins autant de degrés physiques utilisent un QR
supernodal creux avec ordre COLAMD et réflecteurs de Householder, fourni
par faer 0.24.4. L'analyse symbolique est réutilisée lorsque les dimensions
et le motif local sont identiques. La factorisation numérique est refaite
avec les valeurs courantes ; aucun réglage global des fils n'est modifié.

Le QR applique les transformations orthogonales sans construire Q ni les
équations normales. Ce choix repose sur des méthodes établies de calcul
numérique, décrites dans le [cours de David Bindel, Cornell](https://www.cs.cornell.edu/courses/cs6210/2025fa/lec/2025-09-29.html).
La structure creuse de l'entrée ne garantit pas celle des facteurs : le
remplissage dépend du graphe.

La diagonale de R est contrôlée **avant** la résolution triangulaire.
Un bloc douteux repasse par le QR à pivotage, puis si nécessaire par une
SVD dense de `Gᵀ` local, qui donne la solution de norme minimale dans le
rang retenu. Les blocs avec davantage de multiplicateurs que de degrés
physiques utilisent directement cette SVD. La norme de multiplicateurs
reste la norme euclidienne des coordonnées actuelles.

Le seuil est `ε × max(n_lignes, n_colonnes) × ‖Gᵀ_local‖F`.
Il remplace le seuil absolu `1e-12` de l'ancienne SVD de `GGᵀ` ; le choix
du rang numérique peut donc changer. Le contrôle des diagonales ne
constitue pas une estimation générale du conditionnement. Les critères
physiques finaux, le mode strict, les reprises et les restaurations
restent ceux du solveur. Un essai non fini pendant la recherche linéaire
peut encore être raccourci.

## Témoins physiques

Un corps soumis à `(2, 1, 3) N` est maintenu par deux guides de directions
`x` et `(cos(a), sin(a), 0)`, plus les quatre blocages complémentaires.
Ses deux réactions valent exactement `2 − cot(a)` et `1/sin(a)`.
Avec `a = atan(δ)`, l'ancienne projection échoue dès `δ = 1e-6` dans le
diagnostic. Le nouveau solveur atteint la tolérance jusqu'à `δ = 1e-12`,
sans déplacer ce corps déjà en équilibre. À `δ = 0`, la force selon y
n'a plus de réaction possible : les deux versions refusent l'équilibre
et restaurent l'état.

Les quatre familles chronométrées sont entièrement décrites par les
programmes archivés : chaîne soudée sous poids, pivots indépendants à
ressort, chaîne soudée refermée au sol et chaîne articulée sous gravité.
Cette dernière mesure 1 m, pèse 1 kg et porte des ressorts de raideur
`10N N·m/rad` pour N éléments ; elle part horizontale et se déforme.
Un bilan indépendant compare chaque couple de ressort au moment de tous
les poids en aval, ainsi que les réactions et les fermetures de pivots.
La boucle redondante vérifie l'équilibre des corps et l'orthogonalité
des réactions aux six auto-contraintes analytiques.

## Mesures

Les versions alternent sur le même CPU Linux, dans un processus neuf par
échantillon, avec un fil demandé et trois répétitions. Les temps excluent
les imports, la construction du modèle et les contrôles Python. `VmHWM`
est lu immédiatement après la statique ; il inclut Python et ses imports.
Les poses complètes, réactions, résidus et étendues sont conservés.

Chaîne articulée, `strict=True`, `tol=1e-8`, quatre évaluations et une
tentative dans tous les échantillons :

| Corps | Avant | Après | Gain | Pic avant | Pic après |
|---:|---:|---:|---:|---:|---:|
| 8 | 0,658 ms | 0,716 ms | 0,92× | 39,3 Mio | 39,5 Mio |
| 16 | 1,358 ms | 1,294 ms | 1,05× | 39,3 Mio | 39,7 Mio |
| 32 | 5,632 ms | 2,441 ms | 2,31× | 40,0 Mio | 39,6 Mio |
| 64 | 26,065 ms | 4,685 ms | 5,56× | 43,4 Mio | 40,1 Mio |
| 128 | 152,531 ms | 9,843 ms | 15,50× | 57,4 Mio | 41,0 Mio |
| 256 | 1 076,203 ms | 19,910 ms | 54,05× | 112,4 Mio | 42,7 Mio |

Le surcoût médian du petit cas à huit corps est de **8,9 %**. Le coût fixe
de décomposition et de QR n'est donc pas avantageux partout. À 256 corps,
les trois durées sont dans `[1,063 ; 1,076] s` avant et
`[0,01981 ; 0,01996] s` après.

Les positions compensées de tous les corps diffèrent de moins de
`2,23e-16 m` entre versions ; les matrices d'orientation de moins de
`4,45e-16` par composante. Le contrôle indépendant des moments reste sous
`1,22e-11 N·m`, celui des réactions sous `1,26e-11 N` pour les deux versions.
Sur la version actuelle, les réactions sont à moins de `1,34e-13 N` de la
référence analytique.

Les autres topologies, à leur plus grande taille mesurée :

| Famille | Corps | Avant | Après | Gain | Pic avant → après |
|---|---:|---:|---:|---:|---:|
| Chaîne soudée | 512 | 2,974 s | 2,957 ms | 1 005,86× | 257,0 → 40,5 Mio |
| Pivots indépendants | 512 | 4,004 s | 7,580 ms | 528,22× | 270,8 → 43,0 Mio |
| Boucle soudée redondante | 256 | 10,138 s | 2,872 s | 3,53× | 129,1 → 45,1 Mio |

Les chaînes soudées ont zéro degré libre et sont déjà en équilibre :
leur gain mesure principalement la projection des réactions. Les pivots
indépendants convergent en deux évaluations. Ces gains ne décrivent pas
une accélération générale de mille fois des simulations multicorps.
Les poses exportées de ces trois familles sont identiques entre versions.

La boucle à 256 corps conserve le critère commun d'équilibre de `1e-8`
en unités SI. Son résidu passe de `2,92e-9` à `6,09e-11`. Le défaut relatif
d'orthogonalité aux auto-contraintes, `‖Zᵀλ‖₂/(‖Z‖F ‖λ‖₂)`, passe de
`1,09e-10` à `3,61e-14`. Le contrôle de norme minimale à `1e-12` est exigé
de la nouvelle projection ; l'erreur du témoin est conservée. Le coût
augmente encore fortement entre 128 et 256 corps sur ce repli dense.

Princeton, dont les seules contraintes sont à l'encastrement, apporte un
contrôle sur une famille peu concernée par le changement : à 480 intervalles,
la médiane passe de `0,9421` à `0,9328 s`, avec des étendues respectives
`[0,9391 ; 0,9432]` et `[0,9294 ; 0,9392] s`. Les coordonnées du bout sont
identiques. Ce petit écart n'établit pas un gain supplémentaire général.

## Validation et limites

`ci/local.sh --bancs` passe : format, Clippy, **35 tests Rust**, **60 tests
Python**, **41 groupes de vérification**, **46 bancs rapides**, **9 bancs
de contact** et référence d'API à **185 entrées**. Les **132 équilibres**
alternés tiennent la tolérance stricte et leurs bilans physiques séparés.
Les **400 paliers Princeton** passent pour les deux formulations, ainsi
que les trois témoins d'allongement exact `2^-60`.

Les tests Rust couvrent la quasi-dépendance, les multiplicateurs de norme
minimale, les colonnes redondantes d'un grand bloc, les permutations,
le changement des coefficients et du motif, les alias de corps et le
masque dynamique. Le test mécanique des guides échoue avec l'ancienne
extension et passe avec la nouvelle. Le petit exemple
[`qr_contraintes`](../examples/qr_contraintes.rs) est un diagnostic de
l'API QR de faer sur une matrice singulière ; il illustre pourquoi la
diagonale doit être contrôlée avant résolution.

Une grande composante redondante peut encore imposer une SVD dense.
Les matrices publiques d'analyse, les autres projecteurs et les replis
singuliers de Newton restent à traiter. La normalisation physique et le
choix des unités ne sont pas changés par ce travail. Le diagnostic du câble
condensé extrême reste en échec ; sa validité physique n'est pas rétablie.
La couverture généraliste, le temps réel et la comparaison industrielle
restent des fronts ouverts de l'[objectif complet](OBJECTIF_MBDYN.md).

## Reproduction

Le paquet témoin utilise les mêmes fichiers Python, avec la seule
extension compilée au commit de référence. Les empreintes des extensions,
des programmes et des sources sont liées aux résultats par le validateur.
Le JSON des échantillons est compacté pour l'archivage sans changer ses
valeurs ; les espaces terminaux du journal négatif sont retirés.

```bash
.venv314/bin/python ci/mesure_contraintes_alternee.py \
  --ancien-pythonpath /chemin/vers/le/temoin \
  --sortie docs/bancs/contraintes-alternees.json
.venv314/bin/python ci/bilan_contraintes.py \
  docs/bancs/contraintes-alternees.json docs/bancs/contraintes-comparaison.json

# Exécuter avec chaque extension ; conserver aussi les échecs du témoin.
.venv314/bin/python ci/diagnostic_contraintes_proches.py

RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_positions_compensees.py \
  --sortie docs/bancs/contraintes-princeton-strict.json
.venv314/bin/python ci/mesure_statique_alternee.py \
  --ancien-pythonpath /chemin/vers/le/temoin \
  --sortie docs/bancs/contraintes-princeton-alterne.json
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/diagnostic_cable.py \
  --sortie docs/bancs/contraintes-cable.json

# Après collecte des archives et journalisation de ci/local.sh --bancs :
.venv314/bin/python ci/valide_contraintes.py --journal-ci /tmp/contraintes-ci.log
```

[Échantillons alternés](bancs/contraintes-alternees.json),
[comparaison calculée](bancs/contraintes-comparaison.json),
[directions proches avant](bancs/contraintes-proches-avant.json),
[directions proches après](bancs/contraintes-proches-apres.json),
[Princeton strict](bancs/contraintes-princeton-strict.json),
[temps Princeton](bancs/contraintes-princeton-alterne.json),
[câble](bancs/contraintes-cable.json),
[CI](bancs/contraintes-ci.log),
[tests Python](bancs/contraintes-tests-python.log),
[contrôle négatif](bancs/contraintes-negatif-avant.log),
[empreintes et validation](bancs/contraintes-validation.json).
