# Inertie vérifiée : changer les dépendances avant de réduire la précision

Campagne du 8 septembre 2026, après le [contrôle par facteurs](CONTROLE_FACTEURS_PROTOTYPE.md).
Le nouveau prototype calcule les mêmes champs, normes et majorants, avec
une préparation moins coûteuse. Sur 512 poutres, le service complet passe
de **0,673 à 0,516 s**, contre **0,715 s pour la LU corrigée du même lot**.
Les 36 essais physiques et les 24 contrôles passent. La version publique
reste **0.11.0** : le transfert de ce prototype dans l'API reste à réaliser.

## Résultats à service et tolérance communs

Médianes de trois mesures, précédées d'une chauffe par configuration.
Chaque passage emploie un processus neuf, CPU 8, bibliothèques à un fil.
Chaque essai traite 257 fréquences entre 0 et 40 Hz et six charges terminales,
avec comparaison physique par charge et sur toutes leurs combinaisons.
La tolérance relative est 10⁻⁶, sur masse, déformations et ports.

| Poutres | Variante | Préparation (s) | Réponses (s) | Total (s) | Admis |
|---:|---|---:|---:|---:|---:|
| 32 | Inertie binaire, contrôle | 0,0414 | 0,1024 | **0,1439** | 4/4 |
| 32 | Inertie Decimal, contrôle | 0,0503 | 0,1031 | 0,1534 | 4/4 |
| 32 | LU corrigée | 0,0010 | 0,1021 | **0,1031** | 4/4 |
| 128 | Inertie binaire, contrôle | 0,0853 | 0,1285 | **0,2134** | 4/4 |
| 128 | Inertie Decimal, contrôle | 0,1245 | 0,1299 | 0,2545 | 4/4 |
| 128 | LU corrigée | 0,0012 | 0,2212 | 0,2225 | 4/4 |
| 512 | Inertie binaire, contrôle | 0,2870 | 0,2287 | **0,5160** | 4/4 |
| 512 | Inertie Decimal, contrôle | 0,4442 | 0,2292 | 0,6734 | 4/4 |
| 512 | LU corrigée | 0,0018 | 0,7134 | 0,7151 | 4/4 |

Les médianes de phases ne s'additionnent pas nécessairement à celle du total.
Le gain total sur le contrôle Decimal est **6,2 %, 16,1 % et 23,4 %**.
Sur la LU, le nouveau candidat économise **4,1 % à 128** et **27,8 % à 512**,
mais coûte encore **39,6 % de plus à 32 poutres**. À 512, les trois mesures
vont de 0,5150 à 0,5171 s pour le nouveau candidat et de 0,6731 à 0,6756 s
pour le contrôle Decimal. Aucun nouveau ratio Exudyn n'est calculé : ses
configurations HCB n'ont pas été relancées dans cette campagne. Les résultats
[Exudyn précédents](CONTROLE_FACTEURS_PROTOTYPE.md) gardent leur propre lot.

Les huit champs de contrôle par maillage ont une seule empreinte SHA256.
Les normes et tous les retours de majorants sont également identiques entre
les deux variantes et leurs quatre passages. Les pires erreurs physiques
des contrôleurs restent inférieures à 5,58 × 10⁻¹¹, 4,17 × 10⁻¹⁰ et
2,17 × 10⁻⁹. Ces observations ne prouvent pas l'identité sur toute plateforme.

## Le problème mathématique conservé

Avec D, M et B stockés en binary64, interprétés exactement, on vise

```text
Aγ = DᵀD − γM,
Cγ = [[Aγ, B], [Bᵀ, 0]],       B ∈ ℝⁿˣˢ.
```

Le seuil reste `γ = 0x1.ed7aefb394d2bp+17`, associé au choix de 80 Hz.
M est d'abord prouvée définie positive. Une suite complète de congruences
qui donne l'inertie `(n,s,0)` de Cγ prouve le rang s de B et

```text
xᵀDᵀDx > γ xᵀMx,    pour x ≠ 0 et Bᵀx = 0.
```

La [preuve d'inertie](INERTIE_COMPLEMENT_PREUVES.md) reste applicable,
y compris lorsque Aγ est singulière. Une permutation physique P transforme
le KKT par la congruence diag(P,I) ; elle ne change donc pas le problème
certifié. B, Phi et les facteurs de réponse ne sont pas permutés dans le
solveur de réponse : cette permutation appartient seulement au certificat.

### Lemme d'encadrement binary64

Si t est le résultat correctement arrondi d'une opération exacte finie z,
alors son prédécesseur et son successeur binary64 encadrent z. Ce choix
est volontairement conservateur, même si l'opération est exacte. Lorsque
z sous-flue vers zéro, les voisins ±2⁻¹⁰⁷⁴ conservent l'encadrement sous
l'hypothèse de sous-flux graduel. Un voisin infini entraîne un refus.

Pour deux intervalles finis a et b, leurs opérations sont encadrées par
les extrémités monotones pour + et −, les produits extrêmes pour ×, et les
quotients extrêmes pour ÷ lorsque b ne contient pas zéro. Les identités
avec zéro ou un exact évitent des opérations sans élargissement artificiel.
Un carré utilise sa dépendance particulière et une borne inférieure ≥ 0.
Le noyau calcule ces expressions, puis `nextafter` vers l'extérieur.

Le contrat requiert IEEE binary64, arrondi au plus proche et sous-flux
graduel. Le noyau vérifie le mode d'arrondi et les bits de deux résultats
subnormaux ; une comparaison flottante seule pourrait être trompée par DAZ.
Les modes FTZ, DAZ et les trois modes dirigés sont refusés par les tests.
Les modes de contrôle ne sont jamais changés par le noyau. La compilation
interdit fast-math et les contractions FMA. Ce contrat vise le compilateur
et la plateforme éprouvés ; une option de compilation seule ne prouve pas
la conformité de tous les environnements possibles.

### Invariant des Schur

Si le bloc courant H et ses liaisons V sont encadrés, les intervalles
calculés pour `S − V H⁻¹ Vᵀ` contiennent le Schur exact. Les dépendances
entre occurrences peuvent élargir le résultat, mais ne rendent pas son
inclusion fausse. Un pivot scalaire doit exclure strictement zéro. Pour
un pivot symétrique 2×2, un déterminant négatif prouve un signe de chaque ;
un déterminant positif et une trace signée prouvent deux signes identiques.
La formule d'inverse et les mises à jour sont elles aussi encadrées.

Par induction, la somme des signatures des pivots est celle du KKT exact.
Seuls les coefficients exactement encadrés par `[0,0]` peuvent être retirés.
Un pivot ambigu, un dépassement, un budget épuisé ou un défaut de plateforme
refuse le certificat. Une signature complète différente réfute la coercivité
au seuil choisi. Il n'y a ni régularisation, ni remplacement du Gram exact
par le Gram arrondi, ni repli automatique en précision supérieure.

## Pourquoi changer l'ordre

Les [sondes conservées](bancs/inertie-binaire-sondes-2026/README.md) montrent
que l'ordre naturel échoue sur 512 poutres à 16 chiffres Decimal et en
intervalles binary64. Une dissection par séparateurs passe dans les deux
cas. Cette observation ne constitue pas un théorème de stabilité universelle.

L'algorithme construit le motif booléen de DᵀD, uni à celui de M. Il peut
conserver des arêtes dont la valeur physique exacte se compense ; il n'en
supprime aucune nécessaire. Dans chaque composante, deux parcours BFS
choisissent un départ puis une couche proche de la moitié des sommets.
Cette couche est éliminée après les deux côtés ; les sous-problèmes sont
traités de même jusqu'à huit sommets. Les multiplicateurs restent à la fin.

Une arête d'un graphe non orienté relie deux niveaux BFS dont les distances
diffèrent d'au plus un. Retirer une couche sépare donc les niveaux inférieurs
et supérieurs. Tant qu'un pivot demeure dans l'un des côtés, ses Schur
n'introduisent pas de liaison avec l'autre côté. Le remplissage reste permis
sur les séparateurs et les multiplicateurs. Un pivot double traversant cette
partition pourrait changer ce comportement ; aucune borne générale de
remplissage n'est déduite de l'heuristique. Le certificat contrôle toujours
les pivots réellement obtenus, quelle que soit cette structure.

Sur les consoles du lot, la largeur maximale passe de quatre à six voisins.
L'ordre est donc plus coûteux en opérations, mais il permet ici les intervalles
binary64. La suppression de longues dépendances est une explication plausible
de cette réussite, pas une estimation démontrée des largeurs d'intervalles.
Une analyse quantitative exigerait les marges des pivots et la croissance
des mises à jour ; la profondeur de l'arbre seule ne suffit pas.

## Coûts réellement comptés

| Poutres | Certificat Decimal (s) | Ordre + certificat binaire (s) | Dont ordre (s) |
|---:|---:|---:|---:|
| 32 | 0,01223 | 0,00305 | 0,00173 |
| 128 | 0,04919 | 0,01000 | 0,00675 |
| 512 | 0,19774 | 0,04388 | 0,03215 |

À 512, le noyau C++ lui-même prend environ 0,00578 s pour masse, assemblage
et élimination ; ce chiffre exclut la permutation, le pont Python, les copies,
les empreintes et la sérialisation de la preuve. Le tableau compte ces coûts.
La préparation du contrôleur par facteurs prend encore 0,091 s et les
réponses 0,229 s ; l'élimination n'est plus le coût dominant.

La validation CSR vérifie les indices réels et les frontières de lignes
avec des tableaux, au lieu d'un appel NumPy par ligne. Une paire d'indices
adjacents doit être strictement croissante sauf à une frontière de ligne.
Les lignes vides répètent des frontières sans supprimer une paire intérieure.
Les drapeaux de format canonique en cache ne sont jamais considérés comme
une preuve. Les doublons, même accompagnés de faux drapeaux, sont refusés.

Le C++ est compilé une fois, avant les processus, puis chargé avant chaque
chronomètre. Sa commande, son source et son binaire sont capturés. Il ne
s'agit pas de recompilation à chaque modèle. Chaque candidat reconstruit QR,
sélection modale, B, permutation éventuelle, certificat, quatre blocs Krylov,
contrôle de profondeur huit et les 257 réponses. La génération de la direction
utilise la même graine, le même eigsh et les mêmes opérations qu'auparavant.
Le worker LU est inchangé. Aucune occurrence n'a été écartée ou relancée.

## Vérification et limites

Les 11 tests du noyau incluent plus de mille confrontations arithmétiques
rationnelles, 72 KKT entiers/rationnels avec signes admis et réfutés,
les pivots 2×2, les masses couplées, plusieurs contraintes obliques,
les permutations, les frontières singulières, le faux positif d'un Gram
arrondi, les budgets, les extrêmes et les modes flottants hostiles.
Les congruences rapportées sont rejouées en Fraction, indépendamment de
l'élimination numérique. Les tests comparent aussi Python et C++.

La gamme binary64 a des limites visibles : des contraintes multipliées
par 2⁻³⁰⁰ ou 2³⁰⁰ peuvent faire sous-fluer ou déborder un carré intermédiaire,
alors que la coercivité exacte demeure inchangée. Ces cas sont refusés,
jamais classés coercifs à tort. Un équilibrage exact par puissances de deux
avec contrôle de représentabilité reste à étudier ; il n'est pas caché dans
le candidat. La réussite sur trois consoles ne démontre pas celle sur des
structures 3D générales, des contraintes denses nombreuses ou des seuils
presque singuliers.

Les sept tests d'archive requalifient les diagnostics, vérifient les copies
et rejouent les six identités distinctes de certificats : trois Decimal et
trois binaires. Ils réfutent notamment un pivot falsifié conservant la bonne
signature globale et une permutation complète différente de celle prévue.
Les grands champs NPY exclus restent identifiés ; leur jugement est conservé,
mais il n'est pas recalculable depuis cette seule archive compacte.

La masse intérieure du contrôleur de réponse reste limitée aux composantes
connexes de taille au plus six. Le certificat n'a pas cette restriction.
Le continuum physique, l'amortissement, la précontrainte, les contacts, les
trajectoires non linéaires, la mémoire comparable et la supériorité générale
sur MBDyn/Simpack restent à traiter. Les réponses et leurs majorants restent
**non certifiés machine** ; seule la coercivité discrète est certifiée.

## Sources primaires et antériorité

- **Zarebavami et al., 2026, SIGGRAPH**, [notice et version du 4 juin](https://arxiv.org/abs/2602.00898v2),
  [texte, sections 2–5](https://arxiv.org/html/2602.00898v2).
  L'article exploite des groupes de mailles, leur graphe quotient et la
  hiérarchie des séparateurs pour réduire le coût des permutations. Il motive
  l'étude conjointe ordre/coût et la réutilisation de la hiérarchie. Notre BFS
  ne reproduit pas son algorithme et ne reprend pas ses facteurs d'accélération.
  Son contexte de Cholesky sur maillages triangulaires ne certifie pas notre
  KKT indéfini ni la propagation des intervalles. Aucun code de l'article lu.
- **Ost, Schulz, Strash, 2020**, [Engineering Data Reduction for Nested Dissection](https://arxiv.org/abs/2004.11315).
  La notice décrit des réductions de graphe avant dissection pour améliorer
  remplissage et temps. Cette piste concerne les 0,032 s d'ordre qui restent
  sur le grand cas ; aucune de ces réductions n'est implémentée ici. L'effet
  sur un graphe routier ne donne pas un gain garanti pour une structure mécanique.
- **Rump, 2026**, [Verified error bounds for sparse systems, Part II](https://www.tuhh.de/ti3/paper/rump/sparselss_II_final.pdf).
  Le théorème 1.1 utilise inertie, décalages et résidus pour vérifier des systèmes
  creux ; il rappelle une antériorité de 1995. Notre récurrence par intervalles
  vérifie chaque congruence, sans appliquer ce théorème de résidu. Le
  [dossier antérieur](INERTIE_DIRIGEE_SOURCES.md) distingue ces voies et leurs
  hypothèses. Ni l'inertie de Sylvester ni la dissection emboîtée ne sont de
  nouveaux théorèmes revendiqués par Vinkulum.
- **GCC**, [options d'optimisation](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html),
  notamment `-ffp-contract`, `-frounding-math` et `-ffast-math`.
  La commande mesurée emploie GCC 13.3, sans fast-math, avec contraction
  désactivée. Les gardes et contre-épreuves restent nécessaires. Le
  [projet de norme C++, bibliothèque mathématique](https://eel.is/c++draft/c.math)
  documente les fonctions flottantes disponibles ; l'encadrement ci-dessus
  ajoute ses hypothèses explicites de plateforme.

## Reproduction

Depuis le dépôt, avec le Python gelé utilisé pour les mesures :

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 \
  /tmp/vinkulum-release-0.11.0-final/venv/bin/python \
  ci/experience_inertie_binaire.py --campagne /tmp/nouveau-lot-inertie \
  /tmp/vinkulum-confrontation-ports-0.10.0-corrigee \
  /tmp/vinkulum-release-0.11.0-final/venv/bin/python
```

Le dossier de sortie doit être neuf. La [campagne compacte](bancs/inertie-binaire-2026/README.md)
contient les sources, modèles, bases, preuves et provenance de compilation.
`PYTHONPATH=ci python ci/archive_inertie_binaire.py` vérifie les archives et
rejoue uniquement les certificats avec le code courant connu, sans exécuter
les sources ou le binaire capturés. Un compilateur C++17 est nécessaire.
