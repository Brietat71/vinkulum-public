# Première garantie temporelle : pendule rigide

Trois traces de la roue **0.15.0** sont confrontées à une intégration
indépendante par intervalles rationnels, avec reste de Taylor explicite.
Les documents de `bancs/trajectoire-certifiee-pendule/` permettent de
recalculer des **majorants d'erreur à tous les nœuds enregistrés**, y compris
l'instant initial. Cette qualification réalise le premier cas de l'étape 5.
Elle ne constitue ni une certification du noyau entier ni une API générale.

## Objet mathématique et observables

Le modèle de référence est un solide de masse 1 kg, de matrice d'inertie
identité en kg·m², dont le centre est à 1 m d'un pivot fixe d'axe y.
La gravité est `-g ez`. En posant `r = (-sin(theta), 0, -cos(theta))` et
`R = Ry(theta)`, l'inertie autour du pivot vaut `Jyy + m L² = 2` :

```text
theta' = omega
omega' = -k sin(theta),   k = g/2.
```

Les conditions initiales sont `theta ∈ [1/2 − 10⁻¹², 1/2 + 10⁻¹²]`,
`omega ∈ [−10⁻¹², 10⁻¹²]`. Le coefficient `k` couvre à la fois `981/200`
et la moitié du binary64 stocké pour `9.81`, soit
`5522539043063071 / 1125899906842624`. Ces intervalles sont arrondis
vers l'extérieur sur la grille dyadique à 128 bits. Les bornes couvrent
toutes les solutions de cette famille, pour un coefficient constant
quelconque dans l'intervalle annoncé.

Les observables exactes de référence sont :

```text
r = (-sin(theta), 0, -cos(theta))
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)]
v = (-cos(theta) omega, 0, sin(theta) omega)
w = (0, omega, 0).
```

La trace native conserve exactement les valeurs binary64 retournées par
`Noyau.simule` : position **partie haute**, matrice ligne par ligne, vitesse
linéaire et vitesse angulaire. Les parties basses compensées ne sont pas
exposées dans cette trace ; la garantie concerne les valeurs effectivement
retournées. La matrice et la position initiales calculées par `sin`/`cos`
flottants ne sont pas supposées appartenir exactement à la variété des
contraintes. Leur écart à la famille idéale est inclus dans les comparaisons,
dès `t = 0`.

Le générateur fait **une seule simulation continue** par pas, avec `rho=0.9`,
`tous=1` et les autres options par défaut : il ne réinitialise pas l'histoire
algorithmique entre deux observations. Il refuse tout nombre de sorties ou
toute date qui diffère de la grille dyadique attendue. Horizon : **1 s** ;
pas : `1/64`, `1/128` et `1/256` s. Le certificat prouve l'écart de la trace
archivée à la référence, sans avoir besoin de supposer la justesse du schéma
numérique qui a produit cette trace.

## Preuve de l'encadrement

### 1. Arithmétique dirigée

Pour chaque rationnel `q = a/b`, avec `b > 0`, les opérations définissent :

```text
bas(q)  = floor(a 2^128 / b) / 2^128
haut(q) = -bas(-q).
```

Les comparaisons et divisions entières établissent `bas(q) ≤ q ≤ haut(q)`.
Chaque opération d'intervalle arrondit ses extrémités dans ces directions.
Les quatre produits d'extrémités couvrent la multiplication. La division
utilisée est seulement la division par un scalaire rationnel non nul.
Les décisions d'acceptation ne font intervenir aucun flottant.

Pour `|theta| ≤ 2`, les polynômes de Taylor de degré 40 de sinus et cosinus
sont évalués par Horner en intervalles. Le reste de Lagrange est borné par
`max|theta|^41 / 41!`, car les dérivées correspondantes sont bornées par 1.
Aucune approximation de pi ni réduction d'argument n'est utilisée. Le
calcul refuse les domaines qui sortent de cet intervalle angulaire.

### 2. Existence et inclusion sur chaque pas

Pour une boîte initiale `X`, un pas positif `h` et une boîte proposée `B`,
le vérificateur exige `X ⊂ B` et :

```text
X + [0,h] f(B) ⊂ intérieur(B),
h max(1, max|k|) < 1.
```

Pour chaque état initial et chaque coefficient constant admissibles,
l'opérateur de Picard agit sur les chemins continus à valeurs dans `B`.
La première inclusion assure qu'il préserve cet ensemble fermé. La seconde
est une contraction en norme uniforme infinie : une constante de Lipschitz
de `(omega, -k sin(theta))` est `max(1, |k|)`. Le théorème du point fixe
donne un chemin solution unique sur le pas, entièrement contenu dans `B`.
La proposition de boîte peut échouer ; seul le contrôle exact autorise
la poursuite.

### 3. Taylor avec reste sur le domaine démontré

Notons `a_n`, `b_n`, `s_n`, `c_n` les coefficients de Taylor normalisés de
`theta`, `omega`, `sin(theta)`, `cos(theta)`. Les identités différentielles
donnent :

```text
(n+1) a_(n+1) = b_n
(n+1) b_(n+1) = -k s_n
(n+1) s_(n+1) = sum_(j=0..n) (j+1) a_(j+1) c_(n-j)
(n+1) c_(n+1) = -sum_(j=0..n) (j+1) a_(j+1) s_(n-j).
```

L'évaluation par intervalles sur `X` encadre les coefficients jusqu'à
l'ordre 7. L'évaluation sur **tout le domaine de Picard `B`** encadre le
coefficient d'ordre 8 au point intermédiaire du reste de Taylor. Pour
chacune des deux composantes, la boîte suivante est donc :

```text
X_suivant = sum_(j=0..7) A_j(X) h^j + A_8(B) h^8.
```

Le coefficient `k` est évalué comme intervalle à chaque opération ; la
perte de corrélation élargit les boîtes et n'invalide pas l'inclusion.
L'induction sur les pas couvre toute la famille initiale jusqu'à 1 s.

### 4. Écart aux sorties natives

Les formules des observables sont évaluées sur chaque boîte démontrée.
Pour une sortie native exacte `z` et son intervalle de référence `[a,b]`,
`max(|a−z|, |b−z|)` majore l'erreur. Le maximum est pris sur toutes les
composantes et tous les nœuds, initial inclus. Le vérificateur recalcule
ce maximum et refuse une borne annoncée plus petite.

## Résultats archivés

Les lectures ci-dessous sont arrondies **vers le haut** ; les rationnels
exacts figurent dans les certificats. L'erreur de rotation est le maximum
des erreurs des neuf coefficients de matrice, sans conversion en angle.

| Pas (s) | Nœuds | Position (m) | Rotation | Vitesse (m/s) | Vitesse angulaire (rad/s) |
|---|---:|---:|---:|---:|---:|
| 1/64 | 65 | 1,084 × 10⁻⁴ | 1,084 × 10⁻⁴ | 1,580 × 10⁻⁴ | 1,438 × 10⁻⁴ |
| 1/128 | 129 | 2,711 × 10⁻⁵ | 2,711 × 10⁻⁵ | 3,950 × 10⁻⁵ | 3,595 × 10⁻⁵ |
| 1/256 | 257 | 6,777 × 10⁻⁶ | 6,777 × 10⁻⁶ | 9,874 × 10⁻⁶ | 8,987 × 10⁻⁶ |

La largeur maximale des boîtes `(theta, omega)` reste inférieure à
`3,013 × 10⁻¹¹` dans les trois cas. La réduction proche d'un facteur quatre
des majorants est une observation sur ce banc ; ce n'est pas une preuve
générale de l'ordre de convergence du solveur.

Roue qualifiée : `vinkulum-0.15.0-cp314-cp314-linux_x86_64.whl`, SHA-256
`057f89c7256b5bfcb15abee35c831778c64edd3c167e6500ce01280d6fb3a89c`.
Le manifeste conserve l'environnement, les empreintes des scripts et des
trois documents, ainsi que les résultats de leur relecture autonome.

## Reproduction et contre-épreuves

Depuis la racine du dépôt, bibliothèque standard seule pour la relecture :

```sh
python -S ci/archive_trajectoire_certifiee.py docs/bancs/trajectoire-certifiee-pendule/pendule-256.json
python ci/test_intervalle_temporel.py
python ci/test_taylor_temporel.py
python ci/test_reference_temporelle.py
python ci/test_archive_trajectoire_certifiee.py
```

Le test `python ci/test_trace_pendule_courante.py` exécute en plus le noyau
installé au pas `1/64`, reconstruit son certificat et contrôle les quatre
limites publiées dans le tableau, gelées pour les régressions ultérieures.
Ces cinq suites sont obligatoires dans la CI locale.

Le générateur `ci/qualifie_trajectoire_certifiee.py --sortie DOSSIER_NEUF
--roue CHEMIN_ROUE` utilise le paquet installé dans son interpréteur.
Les contre-épreuves couvrent les arrondis signés, des sommes trigonométriques
rationnelles plus longues, le mouvement libre exact et la séparatrice
`omega(t) = ±2/cosh(t)`, évaluée par une exponentielle rationnelle encadrée.
Un test d'ordre bas montre que retirer le reste ne couvre plus cette
solution analytique. Les tests de documents refusent notamment les bornes
sous-estimées, une trace modifiée avec des bornes devenues insuffisantes,
les nœuds supprimés, les témoins Picard invalides, les boîtes falsifiées,
les paramètres hors contrat et les clés JSON dupliquées.

## Domaine et base de confiance

Le format `vinkulum.trajectoire.pendule.prototype.1` impose ce modèle,
ces incertitudes, ces trois grilles, l'ordre 8 et l'arithmétique à 128 bits.
Il ne passe pas par le vérificateur public des certificats algébriques.
La relecture autonome n'importe ni Vinkulum, ni NumPy, ni SciPy.
Le lecteur limite les documents à 4 Mo, les rationnels à 2 048 bits par
numérateur/dénominateur et leur représentation canonique à 1 300 caractères.

La base de confiance comprend Python et ses entiers/Fraction, les formules
du modèle, l'implémentation des intervalles et de Taylor, et l'argument
mathématique ci-dessus. Le générateur et le vérificateur partagent ces
modules ; les contrôles analytiques servent de contre-épreuves distinctes.
**Cette preuve temporelle n'est pas formalisée dans Lean.** Les théorèmes
Lean du garde de vitesse ne lui transfèrent aucune garantie.

Les empreintes assurent l'identité des artefacts contrôlés, sans prouver
que le compilateur ou le générateur ont produit la trace à partir du modèle
annoncé. Un autre document satisfaisant les mêmes contrôles pourrait aussi
être valide : la vérification n'est pas une authentification de provenance.

Aucune borne n'est annoncée sur une interpolation native entre nœuds,
sur un horizon plus long, sur d'autres mécanismes ou sur le contact,
l'adaptatif, GGL, les changements de rang ou les autres intégrateurs.
Le désaccord de repères à temps long du solide libre reste ouvert.
