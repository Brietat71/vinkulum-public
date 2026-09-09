# Qualification du complément spectral — 8 septembre 2026

Trois cas natifs, deux certificats spectraux par cas (une puis deux
directions), et trois contrôles de champs sur 0–40 Hz avec une direction
retenue. Chaque contrôle couvre 257 fréquences, six charges terminales
et la norme d'opérateur sur leurs combinaisons. Les trois cas satisfont
le seuil relatif 10⁻⁶. Le [compte rendu](../../COMPLEMENT_SPECTRAL_DIRIGE.md)
précise les résultats, les preuves et leurs limites.

`rapport.json.gz` conserve les intervalles décimaux, pivots, budgets,
empreintes, diagnostics, erreurs par charge et erreurs d'opérateur à
chaque fréquence. Les temps proviennent d'une seule exécution et sont
indicatifs ; ils excluent la sélection préalable des directions. Aucun
classement de vitesse n'en est déduit. Le certificat couvre la coercivité
du complément ; les champs ne possèdent pas de certificat machine.

## Provenance

L'expérience emploie la roue figée candidate 0.11.0 dont les fichiers
exécutables ont été confrontés à la livraison définitive ; le
[relevé de version](../version-0.11.0.json) documente cette correspondance.
Le rapport identifie RECORD, extension native et sources Python réellement
chargés, Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, un fil par bibliothèque
et le processeur logique 8. La base Git est
`9062db807da2a406ec8681193a06ae046418cf92` ; les quatre sources
expérimentales effectivement exécutées sont copiées dans `sources/`.

`entrees/` contient les trois jeux natifs binary64 de la confrontation
0.10.0 corrigée. `directions/` conserve les douze directions sélectionnées
et le facteur QR ; seules les premières colonnes servent à cette
qualification. La [sonde et son audit](../complement-spectral-sondes-2026/README.md)
conservent le choix de ces directions et sa provenance.

Les oracles à 70 et 90 chiffres sont réutilisés de la
[confrontation historique](../confrontation-ports-exudyn-0.10.0/manifest.json).
Leurs empreintes et le contrôle de convergence figurent dans le rapport.
Les gros tableaux `.ref.npy` ne sont pas dupliqués ; ils se recalculent
depuis les entrées et les sources de l'oracle déjà archivées dans cette
confrontation (`oracle_champs_ports.py`, `reference_ports_precision.py`,
`confronte_ports_exudyn.py`). Une concordance à deux précisions ne constitue
pas une preuve d'arrondi de l'oracle.

## Reproduction

Utiliser une installation isolée de la version 0.11.0 et les dépendances
indiquées. Depuis la racine du dépôt, préparer un dossier de travail avec
les trois NPZ de `entrees/`, puis reconstruire chaque référence :

```sh
python ci/confronte_ports_exudyn.py --reference \
  /tmp/complement-entrees/n32-f40.npz /tmp/complement-entrees/n32-f40.ref.npy
```

Répéter pour 128 et 512, puis comparer les empreintes des trois tableaux
aux `reference_sha256` du rapport. Lancer une nouvelle qualification dans
un répertoire distinct, en choisissant un CPU disponible :

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 \
python ci/experience_complement_spectral.py --cpu 8 --mesurer \
  /tmp/complement-entrees docs/bancs/complement-spectral-2026/directions \
  /tmp/complement-nouvelle-qualification
```

Les scripts de `sources/` sont les instantanés de l'exécution conservée.
Pour une reproduction historique après évolution du code, utiliser leurs
copies à la place des scripts correspondants dans un checkout de travail
isolé. Conserver l'archive existante et écrire les nouveaux résultats à part.

## Contrôles de l'archive

```sh
python ci/experience_complement_spectral.py --verifier docs/bancs/complement-spectral-2026
python ci/test_archive_complement_spectral.py
```

Le premier contrôle relit inventaire, empreintes et résultats. Le second
ajoute la cohérence sémantique des intervalles, pivots, conversions,
dimensions, contraintes NPZ, références et résultats, avec des mutations
adversariales du rapport. Il ne relance pas le solveur ni les oracles.
Le manifeste couvre tous les fichiers de cette archive sauf lui-même ;
l'ajout de ce README n'a modifié aucun résultat ni source mesurée.
