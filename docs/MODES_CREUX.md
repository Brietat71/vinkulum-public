# Analyse modale creuse locale

`Noyau.modes_creux(combien=6, t=None, decalage=None, tol=1e-8,
maxiter=1000)` recherche les valeurs propres signées proches du décalage du
problème local

\[
K_s x + G^T\eta = \lambda Mx,\qquad Gx=0,\qquad
K_s=(K+K^T)/2.
\]

Les formes sont normalisées dans la masse originale : `X.T @ M @ X ≈ I`.
Les résultats contiennent les formes, réactions, valeurs propres en s⁻²,
fréquences positives en Hz et taux non amortis associés aux valeurs négatives.
L'appel conserve l'état du modèle. SciPy est requis.

```python
from vinkulum import Noyau
n = Noyau([0, 0, 0])
n.corps('masse', 1., [1.,0.,0.,0.,1.,0.,0.,0.,1.], [0.,0.,0.])
r = n.modes_creux(3)
print([m['valeur_propre'] for m in r['modes']])  # modes rigides nuls
```

## Algorithme et coût

La masse native doit être définie positive et diagonale par blocs de corps
6 × 6. Une masse couplée entre corps est refusée. La normalisation utilise
les facteurs de Cholesky de ces blocs. Une projection des contraintes et
une factorisation KKT creuse fournissent l'opérateur inverse utilisé par
ARPACK. Le décalage par défaut est légèrement négatif, à l'échelle de la
raideur normalisée ; un décalage explicite est exprimé en s⁻².

Le petit problème de Rayleigh–Ritz est recalculé dans les matrices physiques
originales, avec accumulation `longdouble` lorsqu'elle est prise en charge.
La résolution du petit problème reste en binary64. `precision_ritz_bits`
annonce la précision effectivement utilisée (64 bits significatifs sur la
plateforme qualifiée ; elle peut rester à 53 sur une autre plateforme).
Ce raffinement limite la cancellation des petites valeurs propres sans
changer les seuils de validation.

Il n'y a pas de construction d'une base dense du sous-espace admissible.
Cependant, les petites familles de contraintes utilisent une QR dense de
leurs normales. Les grandes familles peuvent aussi y revenir lorsque la
projection creuse ne permet pas la décision de rang : `repli_qr_dense` le
signale. Demander tout le spectre sans contraintes utilise une résolution
dense, signalée par `spectre_complet_dense`. « Creux » ne garantit donc pas
un coût mémoire linéaire pour tout modèle.

## Portée des contrôles

`tol` porte sur les résidus calculés et l'orthogonalité, **pas sur une borne
certifiée d'erreur de fréquence**. Les résidus des équations physiques sont
recalculés dans les matrices originales, en plus des contrôles normalisés.
Le champ `residu_lambda_s2` est une estimation numérique issue du résidu,
non une inclusion spectrale validée par intervalles.

Le rang est numérique : `rang_certifie=False`. Une dépendance supprimée est
signalée par `reduction_contraintes_numerique` et son défaut calculé par
`defaut_reduction_contraintes`. Par exemple, les lignes `(1,0,0)` et
`(1,10⁻¹⁸,0)` sont exactement indépendantes mais peuvent être déclarées
numériquement dépendantes. Le spectre exact peut changer dans ce voisinage.
Un test conserve ce contre-exemple. Les certificats de
[quotient structuré](CERTIFICATION_CREUSE.md) ont un autre contrat, avec
hypothèses explicites ; ils ne sont pas automatiquement produits ici.

`certification_spectrale=False` s'applique aussi au signe annoncé. Une
recherche près d'un décalage peut manquer une valeur négative plus éloignée.
Ce calcul n'établit pas la stabilité complète : il utilise la partie
symétrique de la raideur locale, sans amortissement ni gyroscopie. Les
contacts, l'aérodynamique et les contraintes non holonomes sont hors domaine.
Un état qui n'est pas un équilibre ne devient pas un équilibre par cet appel.

## Qualification reproductible de la roue 0.17.0

Les références sont assemblées indépendamment dans `ci/modeles_modes.py` :
flexibilité analytique du modèle discrétisé pour la console, énergie
cinétique et potentiel en coordonnées angulaires pour les chaînes. La
console emploie la même raideur de flexion dans les deux plans et des inerties
de masse distinctes. Elle ne mesure pas l'erreur de discrétisation vis-à-vis
d'une poutre continue.

La grille comprend les trois familles `console`, `articulee`, `gravite`,
chacune aux tailles 4, 16, 64, 256, 512 et 1 024. Les 18 cas passent le seuil
commun `1e-8` sur l'erreur relative de fréquence, l'orthogonalité massique et
les résidus relatifs recalculés. Le maximum d'erreur de fréquence est environ
`4,9e-11`. Le maximum de résidu physique relatif recalculé est `6,26e-9`,
celui de fermeture relative des contraintes `2,13e-12` et le défaut
d’orthogonalité massique `1,78e-15`. Ces indicateurs ont des normalisations
différentes ; ils ne sont pas des bornes d’erreur de fréquence.
Les matrices, formes et réactions sont conservées avec leurs empreintes dans [le manifeste](bancs/modes-creux-0.17.0/qualification.json).

Le brouillon protégé antérieur à cette qualification est figé dans
`ci/modes_creux_avant_qualification.py` et exécuté sur le même module natif.
Sur la console à 1 024 éléments, il donne environ `2,46e-8` d'erreur relative
de fréquence, au-delà du seuil commun. Ce résultat défavorable est conservé
pour attribuer le gain au raffinement modal. Il ne s'agit pas d'une régression
d'une API livrée en 0.16.0 : ce brouillon n'était pas publié.

```sh
OPENBLAS_NUM_THREADS=1 python ci/qualifie_modes_creux.py \
  --verifier docs/bancs/modes-creux-0.17.0
python ci/test_archive_modes_creux.py
```

La relecture recalcule références, résidus et orthogonalité sans rappeler le
solveur modal. Elle reste une vérification flottante, distincte d'une preuve
formelle. Les contre-épreuves modifient fréquences, formes, réactions et
statut de certification. Quatorze tests du paquet couvrent également les
modes rigides, instables, les répétitions, redondances, rotations, changements
d'unités, décalages et refus. Aucun gain de temps face à un concurrent n'est
inféré de cette campagne de précision.
