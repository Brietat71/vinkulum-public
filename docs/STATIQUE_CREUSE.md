# Factorisation creuse en statique

Étape après l'assemblage local des poutres, sur la branche de travail après
0.6.2. Sur Princeton à 60 intervalles, le calcul direct passe de **0,406 s
à 0,231 s**, soit **1,75× plus rapide**, avec le même écart de **7,08 µm**
à la référence MBDyn raffinée. Depuis le premier correctif de robustesse,
le gain cumulé est de **11,87×**. Ces rapports comparent des versions de
Vinkulum, pas les deux moteurs.

## Résolution

Au-delà de 160 inconnues, le Newton statique utilise la factorisation
creuse existante lorsque moins d'un quart des coefficients sont non nuls,
et la factorisation dense faer sinon. Les petits systèmes conservent leur
LU nalgebra. Seuls les zéros exacts sont omis : aucun petit coefficient
physique n'est filtré.

L'analyse symbolique est réutilisée dans un appel Newton tant que la
dimension et le motif des coefficients ne changent pas. Le cache reste
local à cet appel, indépendant de celui du solveur dynamique.

Chaque résolution optimisée contrôle la finitude et le résidu du système
équilibré : `norm(A*x-b) <= 1e-9*(1+norm(b))`. Un échec de factorisation
creuse revient au LU nalgebra ; une solution initiale inutilisable conduit
au repli SVD existant. Pendant la recherche amortie, une correction
inutilisable entraîne le rejet de l'essai. Les critères finaux d'équilibre
physique et de contraintes restent inchangés.

L'équilibrage parcourt les matrices par colonnes et réutilise leur
allocation. L'assemblage reste dense avant conversion : cette étape ne
résout donc pas encore le coût mémoire quadratique des grands modèles.

## Mesures

Même modèle et protocole que l'étape précédente : un fil demandé, un
échauffement puis trois répétitions, construction et démarrage Python hors
chronomètre. Charge pleine de 8,896 N à 45°, sans continuation imposée.

| Intervalles | Assemblage local seul | Avec résolution optimisée | Gain | Écart au bout |
|---|---|---|---|---|
| 10 | 0,0264 s | 0,0263 s | 1,00× | 254,87 µm |
| 20 | 0,0637 s | 0,0633 s | 1,01× | 63,72 µm |
| 40 | 0,1847 s | 0,1400 s | 1,32× | 15,93 µm |
| 60 | 0,4062 s | 0,2315 s | 1,75× | 7,08 µm |

Le chemin de 50 paliers est aussi vérifié à 10, 20 et 40 intervalles.
À 40, il passe de 3,97 s à 3,03 s ; il s'agit d'essais uniques de
validation, pas de médianes de performance.

- [Mesures avant](bancs/statique-princeton-local.json).
- [Mesures après et empreintes du binaire et du source](bancs/statique-princeton-creux.json).

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/mesure_statique.py --sortie /tmp/statique-creuse.json
```

## Contrôles et limites

Le test Princeton couvre désormais 40 intervalles pour exercer le chemin
creux avec contraintes et grands déplacements. Un nouveau test à 180
inconnues vérifie le repli singulier : trente corps soumis à un rappel en
rotation atteignent leur équilibre sans déplacement dans leurs cinq
directions neutres respectives.

Formatage, Clippy et vérification du diff passent ; **25 tests Rust**, **44
tests Python**, **41/41 groupes de vérification**, **46/46 bancs rapides**
et **9/9 bancs de contact** passent. La référence d'API reste à jour à
183 entrées. L'extension a été construite en release et chargée directement
dans le venv Python 3.14 ; l'installation d'une roue n'a pas été rejouée.

Le [bilan archivé](bancs/statique-creuse-bilan.json) contient les journaux,
les empreintes et les gains recalculés. Les positions finales changent de
moins de 6e-17 m entre les deux étapes sur les quatre maillages mesurés.

Un profil séparé à 60 intervalles, avec `VINKULUM_TRACE=1`, attribue environ
71 % du temps instrumenté à la raideur et 8 % à la résolution linéaire,
contre respectivement 39 % et 44 % avant cette étape. La résolution inclut
ici la conversion en triplets et le contrôle du résidu. Ces proportions
servent au diagnostic ; les médianes ci-dessus sont mesurées sans trace.

La précision spatiale de l'élément de poutre reste le prochain obstacle
sur Princeton. Accélérer sa factorisation ne réduit pas l'erreur de
discrétisation. Les autres fronts de l'[objectif généraliste](OBJECTIF_MBDYN.md)
restent ouverts.
