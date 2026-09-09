# Sondes préalables et audit du complément spectral

Cette archive conserve les expériences exploratoires du 8 septembre 2026
ayant précédé la qualification des réponses à 40 Hz. Les fichiers originaux
ont été copiés sans modification depuis
`/tmp/vinkulum-sonde-complement-2026`. Aucun cas, résultat refusé ou diagnostic
défavorable n'a été retiré. Elle distingue deux expériences de portées
différentes :

- `sonde.py`, `rapport.json`, `NOTE.md` et `execution.log` : sélection flottante
  de 1/2/4/6/12 directions sur les 32/128/512 poutres natives. Les estimations
  du complément restent non certifiées, et les temps sont indicatifs.
- `audit_dirige.py`, `audit_dirige.json`, `audit_dirige.log` et
  `AUDIT_DIRIGE.md` : 70 cas adversariaux du certificat dirigé expérimental,
  avec 68 résultats acceptés vérifiés en fractions exactes et deux refus sûrs.
  `audit_entrees.py` et `audit_entrees.json` complètent cet audit par les
  entrées sparse, mutations entre appels et budgets. `environnement.json`
  identifie le manifeste RECORD et l'extension native installés.

Les notes sont des instantanés de ces étapes. Leur mention « aucune réponse
à 40 Hz produite » décrit la sonde initiale ; la qualification physique
ultérieure relève de l'expérience distincte
[`ci/experience_complement_spectral.py`](../../../ci/experience_complement_spectral.py).
Cette archive ne constitue ni une campagne comparative de performances ni
une certification machine des champs physiques.

## Provenance

La base du dépôt était le commit
`9062db807da2a406ec8681193a06ae046418cf92` (version 0.11.0). Le module
`ci/trace_complement_dirigee.py` était un ajout expérimental non encore
commité ; son SHA256 audité figure dans `audit_dirige.json` et le manifeste.
Le calcul emploie la roue figée de Vinkulum 0.11.0, Python 3.14.7,
NumPy 2.5.3 et SciPy 1.18.1. Les scripts de référence exacts proviennent du
dépôt Vinkulum ; aucun code de solveur concurrent ni document tiers n'est
contenu ici.

Les entrées natives et les oracles Decimal existants viennent de
`/tmp/vinkulum-confrontation-ports-0.10.0-corrigee`. Le rapport en conserve
les empreintes pour les trois maillages. Les fichiers
`n32-directions.npz`, `n128-directions.npz` et `n512-directions.npz` sont
référencés par leur SHA256 dans `rapport.json` et `manifest.json`. Ils sont
conservés avec l'archive principale de qualification du complément et ne
sont pas dupliqués dans cette archive de sondes.

## Reproduction et intégrité

Les scripts archivés conservent leurs chemins absolus d'origine : ils
documentent exactement les exécutions rapportées. Pour une autre machine,
adapter ces chemins dans une copie de travail des scripts, avec la même
version de Vinkulum et les dépendances indiquées. `sonde.py` accepte
`--entrees`, `--sortie` et `--n`; les audits importent le certificat
expérimental et les auxiliaires rationnels depuis le répertoire `ci`.
Ne pas écrire les nouvelles mesures dans cette archive.

Les exécutions d'origine utilisaient :

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 \
  /tmp/vinkulum-release-0.11.0-final/venv/bin/python \
  /tmp/vinkulum-sonde-complement-2026/sonde.py
```

Les deux audits s'exécutent avec le même environnement en remplaçant le
dernier argument par `audit_dirige.py` ou `audit_entrees.py`. Aucun autre
solveur n'intervient et aucun temps de ces audits ne sert à un classement.

`manifest.json` donne l'inventaire exact, les tailles et SHA256 de tous les
autres fichiers de ce répertoire, ainsi que les empreintes des sources
auditées et des directions externes. Il ne s'inclut pas dans son propre
inventaire. Les noms, les tailles et les empreintes ont été revérifiés après
la copie ; les 70 cas et leurs deux refus sont conservés intégralement.
