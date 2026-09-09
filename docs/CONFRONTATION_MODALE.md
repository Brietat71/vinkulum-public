# Confrontation modale Vinkulum, Exudyn et MBDyn

Vinkulum 0.17.0 dispose d'un avantage mesuré sur le grand cas de cette
campagne, avec des résultats favorables aux concurrents sur d'autres tailles.
Ces mesures ne démontrent pas une supériorité générale en dynamique multicorps.

Le [protocole](PROTOCOLE_CONFRONTATION_MODALE.md) fixe le modèle, les
configurations, la précision et les budgets. La chaîne articulée est décrite
par les interfaces publiques, avec une référence indépendante issue de ses
énergies. Les **140 exécutions** couvrent cinq tailles, sept configurations,
un échauffement et trois répétitions. Toutes les observations, y compris
les refus, sont dans [l'archive](bancs/confrontation-modale-0.17.0/manifest.json).

## Temps à précision commune

Médiane du temps de processus complet, en secondes. Un chiffre n'est donné
que lorsque les quatre exécutions satisfont le seuil relatif de fréquence
`1e-8` et le contrôle de partie réelle applicable à MBDyn.

| Barres | Vinkulum creux | Exudyn arbre creux | MBDyn LAPACK | MBDyn ARPACK élargi |
|---:|---:|---:|---:|---:|
| 4 | 0,260 | 0,240 | 0,013 | 0,012 |
| 16 | 0,228 | 0,240 | 0,040 | 0,046 |
| 64 | 0,236 | 0,241 | 1,170 | 1,289 |
| 256 | 0,321 | 0,263 | délai > 60 s | délai > 60 s |
| 1 024 | 0,333 | 0,602 | échec d'exécution | échec d'exécution |

À **1 024 barres**, Vinkulum prend 0,3318 à 0,3329 s contre 0,5969 à
0,6027 s pour Exudyn arbre creux : rapport des médianes **1,81×**. Le maximum
RSS des répétitions est **79,3 Mio contre 156,9 Mio**, soit un rapport de
**1,98×**. Les erreurs relatives maximales de fréquence, échauffement compris,
sont respectivement `7,99e-15` et `7,84e-14` ; il s'agit d'écarts mesurés à la
référence, pas de bornes certifiées.

À **256 barres**, Exudyn arbre creux est **1,22× plus rapide** : ses temps
sont compris entre 0,2624 et 0,2648 s, contre 0,3191 à 0,3292 s pour Vinkulum.
À **4 et 16 barres**, MBDyn est nettement plus rapide sur le temps complet.
À 64 barres, Vinkulum et Exudyn arbre creux ont des temps proches. Les imports
Python et les différences d'interface comptent dans ces mesures ; ce tableau
ne compare pas uniquement des factorisations déjà préparées.

## Résultats non classés

Les trois variantes Exudyn et les variantes MBDyn LAPACK et ARPACK élargi
passent les quatre essais aux tailles 4, 16 et 64. La recherche ARPACK courte
ne restitue pas le spectre demandé, sur toute la grille.

À 256 barres, Exudyn cartésien dense dépasse le seuil avec une erreur relative
maximale d'environ `1,10e-6`, et Exudyn arbre dense avec `1,70e-7`. L'arbre
creux reste qualifié. MBDyn LAPACK et ARPACK élargi dépassent chacun le budget
de 60 secondes lors des quatre exécutions.

À 1 024 barres, Exudyn arbre dense atteint environ `6,29e-5` d'erreur relative.
Exudyn cartésien dense et les variantes MBDyn LAPACK et ARPACK élargi signalent
des échecs d'allocation sous le plafond commun de **4 Gio d'espace virtuel**.
Leurs temps ne sont pas classés. Cette observation ne dit pas qu'ils seraient
incapables de traiter le cas avec davantage de ressources.

Le [bilan complet](bancs/confrontation-modale-0.17.0/bilan.json) conserve les
sept configurations, leurs quatre décisions et, lorsqu'elles sont qualifiées,
la médiane, l'étendue temporelle et le pic RSS. Aucun échec n'est transformé
en ratio de vitesse favorable à Vinkulum.

## Portée scientifique

Ce cas teste un spectre local conservatif de chaîne ouverte, sans amortissement,
contact, gyroscopie, boucle fermée ni dépendance de contraintes. MBDyn calcule
un problème général du premier ordre ; Vinkulum recherche ici le spectre de
la partie symétrique de la raideur, et Exudyn dispose d'une formulation en
arbre adaptée à cette topologie. Les différentes portées des interfaces
empêchent d'extrapoler le classement à la dynamique générale.

La campagne montre une capacité utile de Vinkulum sur ce grand problème,
et indique où les interfaces concurrentes sont plus rapides. Elle ne prouve
ni un avantage universel, ni une certification spectrale, ni une garantie
sur les performances d'autres mécanismes. Les preuves algébriques livrées
par Vinkulum restent des garanties distinctes à hypothèses explicites.

## Reproductibilité

Machine : AMD EPYC 7543, CPU 0 autorisé, un fil demandé, processus séquentiels.
Python 3.14.7 ; NumPy 2.5.3 ; SciPy 1.18.1 ; Vinkulum 0.17.0 ; Exudyn 1.11.0.
MBDyn develop, sources `eb3bb5796e99c70e9ee0af072c38aada13a199f6`, compilation
`-O3 -DNDEBUG`, OpenBLAS, ARPACK et UMFPACK activés.
Le binaire a pour SHA-256
`65e45ef27c1399f5249891cc7e20aceeb65f81b48c81d9dde19fbae5b4a544d5`.
Les empreintes des extensions Python et les paramètres de compilation sont
conservés dans `campagne.json.gz`, avec la charge du système pour chaque essai.

```sh
OPENBLAS_NUM_THREADS=1 python ci/confronte_modes.py \
  --sortie /tmp/nouvelle-campagne \
  --mbdyn /chemin/mbdyn \
  --python-vinkulum /chemin/venv-vinkulum/bin/python \
  --python-exudyn /chemin/venv-exudyn/bin/python

OPENBLAS_NUM_THREADS=1 python ci/archive_confrontation_modale.py \
  --verifier docs/bancs/confrontation-modale-0.17.0
python ci/test_confronte_modes.py
python ci/test_archive_confrontation_modale.py
```

La relecture recalcule les références, le spectre demandé et le classement
sans les exécutables concurrents. Les instantanés des producteurs, les modèles
MBDyn et le protocole sont archivés. Les tests détectent notamment un mode
manquant, répété ou imprécis, une partie réelle excessive, une grille incomplète
et un faux succès. Les temps et RSS restent des observations instrumentées.

La relecture distingue la métrique recalculée avec sa référence archivée de
la décision vérifiée avec une nouvelle référence indépendante. Le parallélisme
BLAS peut changer les derniers bits de cette dernière : un écart de `2,79e-14`
a été observé sur la métrique du cas Vinkulum à 256 barres, sans changer la
décision. Le vérificateur exige que les deux références conduisent à la même
décision au seuil inchangé de `1e-8`.
