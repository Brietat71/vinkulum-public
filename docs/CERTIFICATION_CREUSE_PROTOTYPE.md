# Certificats linéaires creux et dépendances structurées

Ce dossier conserve la qualification historique du prototype. Son
[intégration publique en 0.16.0](CERTIFICATION_CREUSE.md) ajoute les API et
le format autonome du quotient structurel, avec une nouvelle qualification.

Ce prototype de l'étape 6 vérifie des systèmes carrés creux et des familles
de systèmes perturbés, sans construire leur inverse. Cinq certificats portent
sur des KKT construits avec les matrices de masse et de contraintes exportées
par la roue **0.15.0**. Le plus grand comporte **2 816 inconnues**.
Le prototype ne modifie ni le noyau ni les API et schémas publics déjà livrés.

## Contrat linéaire

Les coefficients stockés de `A`, `b`, du candidat `x`, des facteurs `L,U`
et des enveloppes sont interprétés comme rationnels exacts. On veut borner
l'écart entre `x` et la solution de chaque système :

```text
(A + ΔA) x_exact = b + Δb,
|ΔA| ≤ D,   |Δb| ≤ d.
```

Les inégalités sont composante par composante ; les coefficients absents de
`D` ont une incertitude nulle. Des poids strictement positifs `w` définissent
`||v||_w = max_i |v_i|/w_i`. Ils sont fournis dans les unités des coordonnées,
sans prétendre définir automatiquement une norme physique.

Le vérificateur impose des facteurs triangulaires à diagonale non nulle.
Ils peuvent être inexacts, non symétriques et comporter des pivots négatifs.
La positivité définie de la matrice n'est pas requise pour ce certificat
d'inversibilité. Une factorisation proposée par SuperLU n'est jamais tenue
pour une preuve.

Les permutations sont archivées explicitement :
`A_tilde[i,j]=A[pr[i],pc[j]]`, `b_tilde[i]=b[pr[i]]`,
`x_tilde[j]=x[pc[j]]`. La même transformation est appliquée aux enveloppes
et aux poids. Les bornes sont remises dans l'ordre original des coordonnées.
Le document contient la matrice originale ; le générateur vérifie aussi
sa reconstruction depuis la proposition permutée.

## Démonstration

On note formellement `R=U⁻¹L⁻¹`, sans jamais le matérialiser. Le produit
`E=LU−A` est accumulé en rationnels, une ligne creuse à la fois.

Pour un vecteur positif `v`, deux substitutions de comparaison calculent
un majorant `B(v) ≥ |U⁻¹| |L⁻¹| v`. La première est :

```text
y_i = arrondi_sup((v_i + sum_(j<i) |L_ij| y_j) / |L_ii|),
```

puis la substitution arrière analogue utilise `U`. Ces inégalités suivent
directement de la résolution triangulaire et de l'inégalité triangulaire,
par induction. Elles peuvent être pessimistes lorsque des annulations sont
essentielles dans les inverses.

Posons `z = B((|E|+D)w)` et `q=max_i z_i/w_i`. Alors, pour toute perturbation
admise :

```text
|I − R(A+ΔA)|w ≤ z,
||I − R(A+ΔA)||_w ≤ q.
```

Le vérificateur exige **`q < 1`**. La série de Neumann établit l'inversibilité
de `R(A+ΔA)`, donc celle de `A+ΔA`, puisque `L,U` sont inversibles.

Le résidu du candidat vérifie
`|(A+ΔA)x−(b+Δb)| ≤ r = |Ax−b|+D|x|+d`. On calcule `s=B(r)` et :

```text
beta ≥ max_i(s_i/w_i)/(1−q),
|x_i−x_exact_i| ≤ min(w_i beta, s_i+z_i beta).
```

Ces deux bornes suivent respectivement de la norme de Neumann et de
`e = R résidu + (I−R(A+ΔA))e`. Chaque arrondi des majorants est dirigé
vers le haut. Le vérificateur recalcule les bornes, puis refuse celles
annoncées en dessous du résultat démontré.

### Arrondis relatifs rationnels

La fonction `haut` n'utilise aucun flottant. Pour `q=n/d>0`, soit
`k=bit_length(n)−bit_length(d)` et `h=2^(k−128)`. Le calcul entier renvoie
`ceil(q/h)h`. On a `2^(k−1)<q<2^(k+1)`, donc :

```text
q ≤ haut(q) ≤ q (1+2⁻¹²⁷).
```

Zéro reste exactement zéro. L'exposant variable évite un plancher absolu
dépendant des unités. Les produits et résidus de coefficients sont exacts ;
les substitutions et bornes positives utilisent cet arrondi supérieur.
Les tests couvrent notamment des systèmes remis à l'échelle par `10⁻²⁰⁰`
et `10²⁰⁰`, sans changer la décision de preuve.

## Dépendances perturbées : ce qui peut être garanti

Une faible perturbation arbitraire ne préserve pas le rang. Considérons :

```text
G0 = [1 0; 1 0],   G_epsilon = [1 0; 1 epsilon].
```

Pour tout `epsilon != 0`, le noyau de `G_epsilon` est nul, alors que celui
de `G0` contient `(0,1)`. Avec masse identité et force `(0,1)`,
l'accélération contrainte passe de `(0,1)` à `(0,0)` ; les multiplicateurs
valent `(-1/epsilon, 1/epsilon)`. Les contre-épreuves rationnelles vérifient
ces faits jusqu'à `epsilon=10⁻¹⁰⁰`. Une certification uniforme conservant
le degré de liberté de `G0` serait fausse pour cette famille.

Le contrat implémenté est donc **structurel**, fourni par le modèle :

```text
G' = T' C',   second_membre = T' d',
T'[base,:] = I exactement.
```

`C'`, `d'`, les forces, la masse et les autres lignes de `T'` disposent
d'enveloppes explicites. Les lignes de base de `T'` ne sont pas perturbables.
Ce contrat n'est jamais déduit d'une proximité de coefficients et n'affirme
pas qu'une géométrie native quelconque le respecte.

Le certificat précédent couvre uniformément le KKT réduit
`[M' C'ᵀ; C' 0]`. Son inversibilité implique le plein rang de `C'` :
une dépendance de ses lignes fournirait un vecteur nul du KKT.
`T'` a plein rang colonne grâce à son sous-bloc identité. Ainsi toutes
les contraintes originales `G'a=T'd'` sont équivalentes à `C'a=d'`.
Aucune ligne originale n'est écartée par une décision de seuil.

L'accélération et la réaction généralisée sont uniques. Pour les
multiplicateurs complets, un représentant met les multiplicateurs réduits
sur les lignes de base et zéro ailleurs ; aucune norme minimale n'est
revendiquée. La réaction proposée est `Cᵀ nu` et sa borne inclut
`|C|ᵀ erreur_nu + D_Cᵀ(|nu|+erreur_nu)`. La positivité physique de la masse
reste une obligation distincte de l'inversibilité algébrique démontrée ici.

## Qualification native et archives

Le banc crée des chaînes verticales à pivots y, avec masse et inertie
unitaires, longueurs unitaires et gravité stockée `9.81`. Les matrices
`M,G` proviennent de `k_c_m_g_creux`, dont la non-mutation est contrôlée.
Le système qualifié est `[M Gᵀ;G 0]` et son second membre contient les
forces de gravité ainsi qu'une accélération de contrainte nulle à cette
pose immobile. Il s'agit d'un système construit à partir de l'export ;
le certificat ne prétend pas tracer une résolution interne du noyau.

| Corps | Inconnues | Majorant de contraction (lecture supérieure) | Maximum numérique des bornes (lecture supérieure) |
|---:|---:|---:|---:|
| 1 | 11 | 6,939 × 10⁻¹⁷ | 0 |
| 4 | 44 | 1,353 × 10⁻¹⁵ | 6,751 × 10⁻¹⁴ |
| 16 | 176 | 7,517 × 10⁻¹⁴ | 1,338 × 10⁻¹² |
| 64 | 704 | 4,998 × 10⁻¹² | 5,591 × 10⁻¹⁰ |
| 256 | 2 816 | 2,841 × 10⁻¹⁰ | 9,682 × 10⁻⁹ |

Le maximum numérique mélange accélérations et multiplicateurs dans leurs
unités fournies, avec `w=1` ; ce n'est pas une norme mécanique homogène.
Les bornes composante par composante dans l'ordre original sont archivées.
Une élimination rationnelle dense indépendante contrôle le cas à 44 inconnues,
notamment ses accélérations exactement nulles. Le test courant à 704 inconnues
reconstruit un certificat depuis le noyau installé, au-delà du budget dense
public de 128 inconnues.

Roue : `vinkulum-0.15.0-cp314-cp314-linux_x86_64.whl`, SHA-256
`057f89c7256b5bfcb15abee35c831778c64edd3c167e6500ce01280d6fb3a89c`.
Les cinq fichiers comprimés et leur manifeste figurent dans
`bancs/lineaire-creux-prototype/`. Le manifeste distingue construction/export
du modèle, proposition, construction du document et relecture dans un
nouveau processus. Ces temps sont des observations de qualification uniques,
pas un classement de performances entre solveurs. Pour le plus grand cas,
la relecture autonome observée vaut environ 0,35 s et le document comprimé
occupe 133 389 octets ; les coûts varient avec le remplissage et les facteurs.

## Vérification et limites

```sh
python -S ci/archive_lineaire_creux.py docs/bancs/lineaire-creux-prototype/chaine-256.json.gz
python ci/test_lineaire_creux.py
python ci/test_dependances_structurelles.py
python ci/test_archive_lineaire_creux.py
python ci/test_lineaire_creux_natif.py
```

Les 17 tests incluent les familles perturbées aux sommets contre une
élimination rationnelle, les pivots négatifs, les permutations, les refus
de rang, les données falsifiées et un cas diagonal à 4 096 inconnues.
Le format autonome `vinkulum.lineaire.creux.prototype.1` couvre le système
linéaire et ses enveloppes ; le contrat structurel dispose ici d'une fonction
de prototype et de contre-épreuves, pas encore d'un format autonome dédié.

Limites explicites : 4 096 inconnues, 200 000 termes par matrice et dans
le défaut accumulé, 2 millions de produits scalaires de coefficients pour
`LU`, coefficients d'entrée sur 2 048 bits, exposant de majoration limité
à ±16 384. Le produit peut densifier : son budget est contrôlé, pas supposé
linéaire. Le lecteur limite le document décomprimé à 16 Mo. Les facteurs
peuvent rendre les comparaisons trop pessimistes ; un refus de preuve ne
démontre pas une singularité.

La base de confiance comprend Python, ses entiers/Fraction, le code de
comparaison et l'argument mathématique ci-dessus. Le générateur et le
vérificateur partagent ce code ; les petits oracles d'élimination constituent
des contre-épreuves distinctes. Aucune preuve Lean du certificat creux n'est
revendiquée. Les empreintes identifient les artefacts sans authentifier leur
provenance mécanique. L'extension des API, la qualification du chantier modal
et les nouvelles confrontations MBDyn/Exudyn restent à réaliser.
