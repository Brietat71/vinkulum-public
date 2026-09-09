# Dimension à retenir pour ouvrir le complément au-delà de 40 Hz

La sonde flottante indique qu'**une seule direction intérieure retenue suffit
sur les trois maillages** pour porter le minorant estimé du complément vers
65,5–65,7 Hz. Deux directions portent cette limite vers 100,4–101,0 Hz.
Le refus spectral actuel est donc compatible avec un complément de codimension
très faible. Ces résultats ne sont pas un certificat du complément et aucune
réponse harmonique nouvelle à 40 Hz n'est produite par cette expérience.

| Poutres | Limite initiale certifiée | 1 direction, estimée | 2 | 4 | 6 | 12 |
|---:|---:|---:|---:|---:|---:|---:|
| 32 | 28,2969 Hz | 65,4639 Hz | 100,3827 Hz | 156,5188 Hz | 212,9584 Hz | 384,5003 Hz |
| 128 | 28,3267 Hz | 65,7001 Hz | 100,9976 Hz | 158,3146 Hz | 216,1096 Hz | 397,0312 Hz |
| 512 | 28,3286 Hz | 65,7149 Hz | 101,0366 Hz | 158,4297 Hz | 216,3135 Hz | 397,8737 Hz |

Les colonnes estimées évaluent `(1-eta)/tau_c` en binary64, puis la racine
divisée par `2*pi`. `eta` provient du certificat spectral dirigé existant sur
le D natif et le facteur QR exacts. **L'emploi d'un eta certifié ne certifie
pas le quotient : tau_c et la soustraction sont ici flottants.**

## Données et construction

L'environnement est la roue figée Vinkulum 0.11.0, Python 3.14.7, NumPy
2.5.3 et SciPy 1.18.1. Les fichiers natifs `n32-f40.npz`, `n128-f40.npz` et
`n512-f40.npz` proviennent de la confrontation 0.10.0 corrigée. Le rapport
contient leur SHA256 et celui des références Decimal 90 chiffres existantes.
Les références sont seulement identifiées ; aucune nouvelle réponse n'est
comparée à elles ici.

`CondensationEnergie` produit R et E. Le facteur énergétique analysé pour la
trace est `K_Q=E^-1 R.T R E^-1`. `InverseSelectionnee` calcule sa trace
massique. `certifier_spectral` fournit le contrôle dirigé de l'écart entre
ce facteur et le `D_I.T D_I` interprété exactement.

`eigsh` agit sur l'opérateur symétrique `M_II^1/2 K_Q^-1 M_II^1/2`, appliqué
par deux résolutions triangulaires, sans former d'inverse dense. Le vecteur
initial est pseudoaléatoire avec graine fixée. Les 12 directions sont
ordonnées par valeur propre inverse décroissante, puis transportées en
coordonnées physiques. Leurs résidus sont ensuite évalués avec **D original**.
La plus grande erreur résiduelle duale relative parmi les 12 directions vaut
environ 1,04e-13, 7,61e-13 et 3,51e-12 selon le maillage ; ces résidus restent
des diagnostics flottants.

La contrainte ciblée est `B.T x=0` pour **les valeurs binary64 stockées de B**,
calculées par `B=M_II @ Phi`. Il ne faut pas identifier exactement cette
contrainte à celle du produit mathématique non arrondi `M_II Phi`.
Les archives `n*-directions.npz` contiennent `phi`, `b`, R au format CSR et E.
La première colonne fournit la sélection à une direction ; les premières
k colonnes fournissent les autres sélections.

Pour chaque sélection, la sonde résout `U=K_Q^-1 B`, forme `H=B.T U` et
`G=U.T M_II U`, puis calcule `tau_c=trace(M_II K_Q^-1)-trace(H^-1 G)`.
Le petit système H est équilibré diagonalement avant sa résolution. Le B
stocké n'est pas modifié par cet équilibrage algébrique interne.

## Conditionnement et vérifications indépendantes

À une direction, le facteur d'annulation `(tau+chi)/tau_c` reste entre
9,70 et 9,77 ; à douze directions, il monte entre 368 et 394. Ce cas natif
ne présente donc pas l'annulation catastrophique de la contre-épreuve
`diag(1e-30,1,2,3)`, mais la future preuve doit toujours encadrer la
soustraction et refuser une borne supérieure non exploitable.

Pour les modes non perturbés, H brut a un conditionnement qui monte jusqu'à
environ 1006 à douze directions ; l'équilibrage ramène son conditionnement à
environ 1. La multiplication des colonnes de B par des facteurs de 1e-6 à
1e6 préserve les limites estimées à environ 1e-10 Hz près après équilibrage,
alors que le conditionnement brut peut atteindre environ 1e24. Le rapport
conserve les deux conditionnements ; le résultat ne justifie pas d'ignorer
le rang ou les erreurs d'arrondi du petit système dans une preuve dirigée.

Pour n=32 seulement, une base explicite massiquement orthonormale du
complément est construite indépendamment. La trace inverse du Gram réduit
est évaluée comme somme positive de carrés d'un facteur triangulaire inverse.
Elle concorde avec la soustraction des traces à moins de 1,7e-12 relatif
pour les cinq sélections. Ce contre-calcul est lui aussi flottant.

La sensibilité examine deux perturbations déterministes : mélange avec des
directions aléatoires orthogonales aux douze modes, puis rotation vers les
modes immédiatement suivants. Cette dernière est plus informative : avec
une direction retenue, une rotation de 1 radian vers le deuxième mode donne
encore une limite estimée de 47,3–47,4 Hz ; à 1,4 radian elle tombe à
32,4 Hz et ne suffit plus pour 40 Hz. Il faut donc certifier le complément
effectivement utilisé, et non se fier au nombre de directions.

## Coûts indicatifs et suite

Une seule exécution, un fil BLAS/OpenMP/Rayon, sans affinité dédiée ni
répétitions : les temps ne sont **pas** des mesures comparatives de performance.
Les 12 modes prennent environ 26/30/43 ms ; la trace complémentaire à une
direction prend environ 2,0/2,2/3,2 ms. Le certificat spectral initial prend
environ 33/137/560 ms. Le coût d'une certification dirigée supplémentaire
de U, H et de la correction de trace n'est pas encore mesuré par cette sonde.

La suite concrète est de certifier d'abord le B à une direction sur les trois
maillages, puis de garder cette direction dans le système couplé et de
condenser seulement son complément. Tous les couplages de masse et de
raideur doivent rester présents. Les résonances globales du modèle physique
subsistent ; une coercivité certifiée du complément ne suffit ni à garantir
le Schur restant ni à valider les champs à 40 Hz.

Reproduction depuis ce répertoire :

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 \
  /tmp/vinkulum-release-0.11.0-final/venv/bin/python sonde.py
```

`rapport.json` conserve les résultats complets. `sonde.py` est autonome
pour le chargement des entrées et ne lit aucun code de solveur concurrent.
