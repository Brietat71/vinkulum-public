# Assemblage local des tangentes de poutre

Étape archivée. La [résolution statique creuse](STATIQUE_CREUSE.md) traite
ensuite une partie du coût de factorisation décrit dans les limites.

Étape suivant le correctif statique de Princeton, sur la branche de travail
après 0.6.2. Le calcul direct à 60 intervalles passe de **2,75 s à 0,41 s**,
soit **6,76× plus rapide**, avec une erreur au bout de **7,08 µm** inchangée.
Ces gains comparent deux versions de Vinkulum ; ils ne démontrent pas une
supériorité sur MBDyn.

## Changement

La raideur retirait auparavant deux champs globaux de forces par coordonnée
perturbée. Cela recalculait toutes les poutres pour chaque colonne.
Désormais, leur contribution est assemblée à partir des tangentes locales
12×12, avec leurs véritables indices globaux. Le reste du champ conserve sa
différentiation globale. Les réactions sont initialisées avec le **modèle
complet**, avant cette séparation, pour préserver la précontrainte.

Les forces constantes dans le repère monde sont retirées des différences
finies où leur dérivée est nulle. Dans `amortissement()`, les forces
conservatives de poutres sont également retirées : dans `forces_a`, elles
ne dépendent pas des vitesses. Les contributions des couples, du mouvement
gyroscopique, des contacts, de l'aérodynamique et des superéléments restent
évaluées.

Le champ des forces physiques, le schéma temporel et les critères
physiques de convergence ne changent pas.

## Statique à précision inchangée

Même programme, modèle, références et protocole que le correctif précédent :
trois répétitions après échauffement, un fil demandé, temps de chargement et
de calcul hors construction et démarrage Python.

| Intervalles | Correctif précédent | Assemblage local | Gain | Écart au bout / MBDyn raffiné |
|---|---|---|---|---|
| 10 | 0.085 s | 0.026 s | 3.24× | 254.87 µm |
| 20 | 0.320 s | 0.064 s | 5.02× | 63.72 µm |
| 40 | 1.240 s | 0.185 s | 6.72× | 15.93 µm |
| 60 | 2.747 s | 0.406 s | 6.76× | 7.08 µm |

Les positions entre les deux versions varient de moins de 2e-16 m sur ces
cas. Le chemin des 50 paliers reste vérifié aux trois maillages 10/20/40 ;
à 40 intervalles, son temps passe de 25,69 s à 3,97 s (essais uniques de
validation, pas des médianes de performance).

## Matrices d'analyse

`k_c_m_z()` est mesuré à un état déformé et tournant imposé, pour ne pas
mélanger ce résultat avec le solveur statique. Trois répétitions après
échauffement ; copie vers Python incluse.

| Intervalles | Avant | Après | Gain |
|---|---|---|---|
| 10 | 3.36 ms | 0.86 ms | 3.89× |
| 30 | 27.31 ms | 4.31 ms | 6.33× |
| 60 | 110.34 ms | 19.29 ms | 5.72× |
| 120 | 458.92 ms | 93.90 ms | 4.89× |

Les quatre projections déterministes archivées de K changent de moins de
1e-16 en relatif sur ces modèles ; celles de C sont identiques. Cela ne
remplace pas le contrôle matriciel indépendant : le nouveau test compare
**toutes les entrées** de K et C à des différences finies du champ complet,
sur deux poutres déformées partageant un nœud, des indices non contigus,
un couple de rappel amorti et des corps tournants.

## Validation et reproduction

Formatage et Clippy passent ; **25 tests Rust**, **43 tests Python**,
**41/41 groupes de vérification**, **46/46 bancs rapides** et **9/9 bancs de
contact** passent. L'API générée reste à jour, 183 entrées.

- [Bilan calculé, empreintes et journaux](bancs/raideur-locale-bilan.json).
- [Statique avant](bancs/statique-princeton-correctif.json) et [après](bancs/statique-princeton-local.json).
- [Tangentes avant](bancs/tangentes-avant-local.json) et [après](bancs/tangentes-local.json).

```bash
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/mesure_statique.py --sortie /tmp/statique-locale.json
RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv314/bin/python ci/mesure_tangentes.py --sortie /tmp/tangentes-locales.json
```

## Ce qui reste

Les matrices retournées sont denses ; le calcul des autres contributions
et la factorisation statique gardent encore des coûts globaux. L'erreur
spatiale des poutres à deux nœuds demande également davantage d'intervalles
que les `beam3` de MBDyn sur Princeton. Les prochains gains doivent traiter
ces coûts et cette précision, puis être jugés au même seuil d'erreur.

[L'objectif généraliste complet reste actif](OBJECTIF_MBDYN.md).
