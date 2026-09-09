# Inertie du KKT par élimination dirigée

Cette note justifie un certificat spectral expérimental sur `ker(B.T)`.
Les contre-épreuves indépendantes sont dans
[ci/test_inertie_dirigee_preuves.py](../ci/test_inertie_dirigee_preuves.py).
Elles utilisent des fractions exactes et un petit témoin d'intervalles ;
elles ne mesurent ni ne valident une implémentation creuse de production.

## 1. Objet exact et conclusion recherchée

Soient `D ∈ R^(m×n)`, `M=M.T ∈ R^(n×n)`, `B ∈ R^(n×s)`, avec `0<s<n`.
Chaque entrée binary64 finie est interprétée comme son nombre dyadique
**exact**, y compris B. On ne remplace pas B par un nouveau produit `M Phi`.
Le paramètre γ est lui aussi un nombre exact spécifié par le contrat.
Posons

```text
Aγ = D.T D − γ M,
Cγ = [[Aγ, B], [B.T, 0]].
```

On suppose M strictement positive, avec une preuve séparée portant sur ces
mêmes entrées. La positivité de la diagonale de M ne suffit pas. Si M est
stockée par blocs découplés, chacun de ces blocs doit être certifié positif.
Une partition prétendue doit aussi exclure exactement tout couplage entre
blocs.

L'inertie `(p,q,z)` compte respectivement les directions positives,
négatives et nulles. Une preuve

```text
Inertia(Cγ) = (n,s,0)
```

établit que B a rang s et que, pour tout `x≠0` tel que `B.T x=0`,

```text
x.T D.T D x > γ x.T M x.
```

Elle donne donc `λ_min((D.T D,M)|ker(B.T)) > γ`. Retourner `lambda_min=γ`
est un minorant valide, sans convertir ce résultat en certification des
réponses, de leur bloc de ports, ou du modèle physique continu.

En effet, si `B y=0`, alors `Cγ [0;y]=0`. La régularité du KKT prouve donc
le rang colonne de B, sans Gram `B.T B`. Pour une base Z de `ker(B.T)`,
l'identité exacte, valable même si la restriction est singulière, est

```text
Inertia(Cγ) = Inertia(Z.T Aγ Z) + (s,s,0).
```

La congruence qui la démontre est donnée dans
[INERTIE_COMPLEMENT_PREUVES.md](INERTIE_COMPLEMENT_PREUVES.md).
Elle n'inverse ni Aγ ni `Z.T Aγ Z`. La construction de Z n'est pas requise
par le certificat numérique.

## 2. Assemblage : certifier D.T D, pas un Gram arrondi

Pour chaque paire `(i,j)`, l'entrée exacte est

```text
a_ij = Σ_k d_ki d_kj − γ m_ij.
```

Chaque produit et chaque somme doivent être encadrés par arrondis dirigés.
`Decimal.from_float(x)` fournit le nombre binary64 exact ; `Decimal(str(x))`
définit généralement un autre nombre. Une précision de travail finie
intervient dans les opérations, pas dans la définition du point d'entrée.
Les conversions rationnelles éventuelles doivent elles-mêmes être
dirigées. Les non-finis, dépassements et divisions non admissibles causent
un refus ; ils ne produisent pas de certificat partiel utilisable.

Il suffit d'assembler un triangle, puis de réutiliser le même intervalle
pour son symétrique. M doit être exactement symétrique dans la convention
d'entrée. Une moyenne arrondie de M et M.T n'est pas une preuve de cette
symétrie. Pour des données creuses avec doublons, il faut soit refuser les
doublons, soit les sommer dans l'arithmétique exacte/dirigée du contrat.
Une canonicalisation préalable en binary64 peut changer le problème.

Si une ligne de D a `t_k` colonnes présentes, elle engendre au plus
`t_k(t_k+1)/2` produits dans le triangle. Le coût d'assemblage dépend de
`Σ_k t_k²`, pas seulement du nombre de valeurs stockées. Les zéros
structurels sont exacts ; aucune entrée simplement petite ne peut être
supprimée sans une preuve supplémentaire d'erreur.

Contre-exemple : avec `ε=2⁻³⁰`, les colonnes `(1,ε).T` et `(1,0).T` ont
un Gram exact strictement positif de déterminant `ε²`. Le Gram binary64
peut devenir `[[1,1],[1,1]]`, de rang un. Certifier ce dernier certifierait
un objet différent, même si la factorisation ultérieure était parfaite.

## 3. Invariant d'élimination

Après une permutation symétrique, partitionnons une matrice exacte restante

```text
W = [[P, F.T], [F, Q]],   S = Q − F P⁻¹ F.T.
```

Si P est inversible, une congruence inversible donne

```text
Inertia(W) = Inertia(P) + Inertia(S).
```

L'invariant numérique est : **chaque coefficient de la matrice de Schur
exacte, déterminée par les entrées initiales et les pivots déjà choisis,
appartient à son intervalle stocké**. Une arithmétique d'intervalles avec
arrondis dirigés préserve cet invariant pour l'expression de S.

Les facteurs d'intervalles ne constituent pas une factorisation exacte
avec des milieux flottants. La preuve porte sur la suite de matrices
exactes cachées dans les intervalles. Il est permis de choisir le prochain
pivot d'après les données approchées ou les intervalles, puis de vérifier
ses conditions. Les propriétés d'un facteur approché ne remplacent pas
ces vérifications.

Les occurrences répétées de P, de F ou d'une variable sont dépendantes.
Les évaluer comme des occurrences indépendantes dans une extension
naturelle d'intervalles **élargit** généralement l'encadrement ; cela ne
retire pas le point exact. En revanche, resserrer arbitrairement cet
encadrement en supposant une dépendance non démontrée n'est pas valide.
On peut intersecter deux encadrements seulement si chacun contient déjà
le même scalaire exact. Une division par le milieu du pivot n'encadre pas
son inverse exact.

## 4. Règles de signature des pivots

Pour un pivot scalaire `[d]`, sa signature vaut `(1,0,0)` si `inf[d]>0`,
et `(0,1,0)` si `sup[d]<0`. Sinon ce pivot n'est pas admissible.

Pour un pivot symétrique `P=[[a,b],[b,c]]`, calculer un encadrement

```text
[Δ] = [a][c] − square([b]),    [t] = [a]+[c].
```

`square([b])` doit contenir tous les carrés : si l'intervalle traverse
zéro, sa borne inférieure est zéro. Le produit générique `[b]*[b]` reste
un encadrement sûr, mais souvent plus large. Calculer seulement les
carrés des deux extrémités et prendre leur minimum serait faux lorsque
l'intervalle contient zéro.

Les règles suffisantes sont :

| Condition dirigée | Inertie prouvée de P |
|---|---|
| `sup[Δ]<0` | `(1,1,0)` |
| `inf[Δ]>0` et `inf[t]>0` | `(2,0,0)` |
| `inf[Δ]>0` et `sup[t]<0` | `(0,2,0)` |

Avec `inf[Δ]>0`, le signe strictement prouvé de a ou c peut aussi
remplacer celui de la trace. Ces règles proviennent du produit et de la
somme des deux valeurs propres réelles. Un déterminant positif seul
ne distingue pas `I₂` de `−I₂`. Une trace positive seule ne distingue pas
une matrice positive d'une matrice de signature `(1,1,0)`.

Le calcul du Schur utilise l'inverse encadré

```text
P⁻¹ = (1/Δ) [[c,−b],[−b,a]].
```

Toutes les divisions exigent `0∉[Δ]`. Si une règle de signature n'est
pas séparée de zéro, il faut choisir un autre pivot, augmenter la
précision dans une nouvelle tentative explicitement comptée, ou refuser.
Un signe de milieu, un seuil de positivité ou une petite valeur propre
flottante ne lève pas ce refus. Une singularité exacte doit rester visible.

La règle ne suppose pas la diagonale de P inversible : `[[0,b],[b,0]]`
avec `b≠0` est hyperbolique et admissible. Par exemple, `D=diag(1,2)`,
`M=I`, `γ=1`, `B=e₁` donnent Aγ singulière mais Cγ régulière de signature
`(2,1,0)`. Le pivot physique/multiplicateur permet une preuve là où une
élimination scalaire obstinée sur le premier zéro échoue.

## 5. Complétude, refus et bande de fréquences

Si tous les `n+s` indices sont éliminés par pivots inversibles dont la
signature est prouvée, leur somme est l'inertie exacte de Cγ et sa nullité
est zéro. Il est alors légitime d'en déduire le rang de B. On ne peut pas
en déduire ce rang d'un simple préfixe de pivots, même tous bien séparés.
Si B est redondant, un vecteur nul multiplicateur subsiste : une
élimination complète ne peut pas réussir correctement.

Avec une signature complète différente de `(n,s,0)`, la coercivité stricte
au seuil γ est réfutée. Avec un pivot non séparé ou un budget dépassé,
le résultat est **inconnu**, pas une preuve de non-coercivité. Ces cas
doivent conserver des motifs distincts dans le diagnostic.

Pour couvrir toutes les fréquences ω stockées, on vérifie exactement

```text
Fraction(γ) > max(Fraction(ω)**2).
```

L'égalité peut aussi suffire mathématiquement puisque la coercivité au
seuil γ est stricte ; exiger `>` conserve une marge rationnelle positive
et le contrat actuel. Il ne faut pas conclure à partir de la seule
comparaison de γ à `fl(ω²)`. Pour `ω=1+2⁻⁵²`, le carré binary64 usuel
est inférieur au carré exact de `2⁻¹⁰⁴`. Un arrondi supérieur du carré
suivi de la comparaison Fraction évite cette ambiguïté.

## 6. Graphe, budgets et coût

Un pivot scalaire crée potentiellement toutes les arêtes entre ses
voisins restants. Pour un pivot de taille deux, la réunion des voisins
des deux indices devient potentiellement une clique. Le budget doit
couvrir ces entrées et les opérations avant leur création. Un Schur
non nul ne peut pas être omis parce que l'arête n'existait pas dans Cγ.

Le petit nombre s de contraintes ne suffit pas à borner le remplissage.
Avec `A=I` et `B=1`, éliminer tôt un couple physique/multiplicateur dont
le pivot vaut `[[1,1],[1,0]]` laisse `I+11.T` sur les indices physiques
restants : le graphe devient complet. En éliminant les physiques avant
le multiplicateur, la même étoile possède un ordre beaucoup plus petit.
Une permutation et une stratégie de secours à pivots doubles sont donc
des composantes du coût à mesurer, pas des détails de la preuve.

Ne retirer que les intervalles exactement `[0,0]` est sûr. Garder une
entrée dont l'intervalle contient zéro est parfois coûteux mais nécessaire
sans preuve d'annulation exacte. Il faut compter l'assemblage, les
intervalles stockés, les opérations dirigées, les recherches/permutations,
les blocs de masse et tous les essais de précision. Un résultat obtenu
après abandon d'un budget ne constitue pas un résultat sous ce budget.

Les preuves ci-dessus n'établissent à elles seules aucun gain par rapport
aux 0,685 s observées sur le certificat de trace à 512 éléments. La
[campagne de 84 essais](INERTIE_CONTRAINTE_PROTOTYPE.md) mesure séparément
ce gain sur trois consoles. Le remplissage et la largeur des intervalles
peuvent rendre cette alternative plus chère ou indécidable sur d'autres
problèmes au budget choisi.
