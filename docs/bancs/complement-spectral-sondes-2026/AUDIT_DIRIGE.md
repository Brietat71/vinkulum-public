# Audit indépendant de trace_complement_dirigee.py

**Aucun faux certificat trouvé.** La relecture algébrique et 70 cas adversariaux
concordent : 68 certificats acceptés vérifiés en fractions exactes, deux refus
sûrs lorsque la précision 16 ne démontre pas le rang. Les données, résultats
et SHA256 du module audité figurent dans `audit_dirige.json`. Le module de
production expérimental et sa suite de tests n'ont pas été modifiés.

## Vérifications mathématiques

- `_resoudre_qr` encadre successivement `R.T y=E B`, `R z=y`, puis `U=E z`.
  Il calcule donc bien `K_Q^-1 B` pour `K_Q=E^-1 R.T R E^-1`.
- Les triangles de H et J peuvent être copiés par symétrie : leurs cibles
  exactes `B.T K_Q^-1 B` et `U.T M U` sont symétriques, même si deux évaluations
  flottantes indépendantes auraient fourni des nombres différents.
- La LDL intervalle encadre les opérations de la LDL exacte. Chaque pivot
  inférieur strictement positif prouve un pivot exact positif ; elle établit
  H SPD et donc le rang colonne plein de B. `_resoudre_ldl` résout bien les
  trois étapes L, D, L.T ; le terme diagonal est divisé avant la remontée.
- La soustraction dirigée `tau-chi` fournit un majorant de la trace positive
  exacte du complément. Sa borne inférieure peut être négative : la positivité
  exacte découle déjà de K_Q/M SPD, du rang de B et de `s<n`. Diviser par le
  majorant positif reste alors sûr, éventuellement très conservateur.
- La relation dirigée existante `K_D >= (1-eta) K_Q`, avec eta<1, se restreint
  au même noyau de B.T. Le quotient par le majorant de tau_c, calculé vers
  le bas et converti vers le bas, donne un minorant valable dans D original.

Chaque certificat accepté a été comparé à `Fraction(float(valeur))`, sans
arrondir auparavant R/E : les fractions exactes représentent séparément R
et E. Les intervalles de tau, chi et tau_c contiennent leurs valeurs exactes.
Le quotient est inférieur à `(1-eta)/tau_c` exact. Une base rationnelle du
noyau de B.T est obtenue par Gauss ; les mineurs principaux dominants de
`Z.T (D.T D-lambda M) Z` vérifient la positivité requise.

## Cas adversariaux et limites de précision

Le lot comprend une masse couplée SPD, des colonnes de B rééchelonnées par
`2^±50`, `2^±200` et `2^±400`, des facteurs R anisotropes jusqu'à
`diag(2^-100,1,2^100,2^200)`, des contraintes quasi redondantes, une annulation
de grande trace et 40 cas pseudoaléatoires reproductibles. Les précisions
16, 40 et 80 sont examinées.

Le contre-exemple le plus pertinent pour la qualité de la borne est

```
B = [[1, 1], [0, 2^-52], [0, 0], [0, 0]]
```

avec R, E et la masse couplée du cas oblique. À 16 chiffres, le rang est
refusé. À 40 chiffres, le rang est démontré mais la borne supérieure de la
trace est environ **1,34e7 fois sa valeur exacte** ; à 80 elle redevient
précise. Cela ne compromet pas la sûreté du certificat ; cela exclut de
confondre « certificat positif obtenu » et « bande cible démontrée ».

À 16 chiffres, l'annulation `R=diag(2^-50,1,2,3), B=e1, M=I` fournit une
borne de trace environ 4,41e15 fois trop grande, toujours sûre. Plusieurs
cas ont une borne inférieure négative de tau_c et ont été vérifiés
exactement sans fausse conclusion.

## Entrées, mutations et budgets

`audit_entrees.py` vérifie les refus de CSR contenant des doublons malgré
des drapeaux `has_canonical_format=True` mensongers, pour D, R et M ; le
triangle inférieur non nul de R ; une masse asymétrique ; les paramètres
de budget booléens. Les seuils exacts de budget Decimal et rectangulaire
sont acceptés, leurs prédécesseurs sont refusés. Une mutation ultérieure
de R est réévaluée dans eta et dans le minorant ; une mutation de B change
l'empreinte, sans réutiliser de cache antérieur.

Deux restrictions restent explicites : les budgets sont des compteurs
logiques, pas des garanties mémoire en octets ; aucune garantie de cohérence
de plusieurs entrées modifiées concurremment pendant l'appel n'est établie.
Les copies internes évitent un cache périmé entre appels, pas une transaction
atomique sur un modèle mutable.

## Couverture conseillée

Les sept tests existants sont cohérents avec leur portée. Deux compléments
seraient particulièrement utiles : le B quasi redondant ci-dessus à 16/40/80
chiffres, et le rééchelonnement exact des colonnes par `2^±400` avec masse
couplée. Ces cas distinguent rang, conditionnement et largeur effective du
certificat. Le fichier `audit_dirige.py` contient leurs données complètes.

Une limite de la composition actuelle est volontaire : eta porte sur tout
l'intérieur. Une erreur QR située entièrement dans les directions retenues
peut donc dégrader le minorant du complément. Par exemple, D=I, R diagonal
`(10,1,1,1)`, M=I, B=e1 donne eta=0,99 et lambda≈1/300, alors que le complément
exact est inchangé et la borne de trace intrinsèque vaut 1/3. C'est une perte
de précision de la borne, pas une violation de sa sûreté. Une perturbation
énergétique directement restreinte serait une amélioration distincte.
