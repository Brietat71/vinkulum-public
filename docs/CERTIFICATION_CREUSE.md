# Certification creuse — API 0.16.0

`vinkulum.certification` expose deux fonctions facultatives, vérifiables par
`verifier_certificat`. Elles copient leurs entrées et ne modifient pas les
solveurs du noyau. La [démonstration du prototype](CERTIFICATION_CREUSE_PROTOTYPE.md)
s'applique au même algorithme par facteurs triangulaires, défaut exact,
substitutions positives et contraction stricte. Aucun inverse dense n'est construit.

## Systèmes linéaires perturbés

```python
from scipy.sparse import csc_matrix
from vinkulum.certification import certifier_systeme_creux, verifier_certificat

preuve = certifier_systeme_creux(
    csc_matrix([[4., 1.], [1., 3.]]), [1., 2.],
    delta_a=csc_matrix([[1e-6, 0.], [0., 1e-6]]),
    delta_b=[1e-7, 1e-7],
)
resultat = verifier_certificat(preuve)
assert resultat['unicite_uniforme']
print(resultat['bornes_composantes'])
```

La garantie couvre chaque système `(A+ΔA)y=b+Δb` tel que
`|ΔA|≤delta_a`, `|Δb|≤delta_b`, composante par composante. Les entrées
absentes de `delta_a` ont une incertitude nulle. Les poids positifs optionnels
définissent `max_i |erreur_i|/poids_i`, avec 1 par défaut. Les bornes restent
exprimées dans les unités des coordonnées fournies.

Les matrices sont creuses SciPy et réelles. Leurs contributions sont converties
en binary64, puis les doublons sont sommés exactement : `1e16, 1, -1e16`
au même emplacement donnent 1. La proposition flottante peut approcher cette
somme rationnelle ; le défaut `LU−A` est toujours vérifié contre la somme
exacte archivée.

SuperLU propose solution, facteurs et permutations. Un candidat `x` peut
remplacer la solution proposée. Avec `facteurs=(L,U)`, un candidat `x` est
obligatoire ; `permutations=(pr,pc)` indique `A_permute[i,j]=A[pr[i],pc[j]]`,
sinon l'ordre reste inchangé. Les poids et enveloppes suivent les mêmes
permutations ; les bornes sont remises dans l'ordre original. Les facteurs
fournis utilisent des permutations sous forme de listes d'entiers Python.
Ils
peuvent avoir des pivots négatifs ; leur triangularité et leurs diagonales
non nulles sont vérifiées. Les calculs flottants ne constituent pas une preuve.

## Dépendances structurées

```python
from scipy.sparse import eye, csc_matrix
from vinkulum.certification import certifier_quotient_structurel, verifier_certificat

preuve = certifier_quotient_structurel(
    eye(2, format='csc'), csc_matrix([[1., 0.]]),
    csc_matrix([[1.], [2.]]), [0], [0., 1.], [0.],
    delta_c=csc_matrix([[.001, .001]]),
    delta_t=csc_matrix([[0.], [.1]]),
)
resultat = verifier_certificat(preuve)
assert resultat['acceleration_unique']
assert resultat['reaction_generalisee_unique']
assert not resultat['multiplicateurs_uniques']
```

Les arguments sont `M,C,T,base,force,d`. La famille est explicitement
`G'=T'C'`, second membre complet `T'd'`, avec `T'[base,:]=I` exactement.
Toutes les contraintes originales sont conservées par cette relation.
Le vérificateur contrôle les lignes de base et leurs incertitudes nulles,
les deux couplages `C',C'ᵀ` du KKT et leurs enveloppes, puis la preuve
linéaire imbriquée. Celle-ci couvre aussi des variations indépendantes plus
larges que les variations corrélées de C : cette majoration est conservative.

L'inversibilité uniforme du KKT implique le plein rang de C'. Le sous-bloc
identité de T' maintient l'équivalence avec les contraintes complètes. Les
autres coefficients de T' peuvent varier, changer de signe ou s'annuler
dans leur enveloppe. `delta_m`, `delta_force` et `delta_d` encadrent aussi
la masse, la force et le second membre réduit.

La réaction généralisée `C'ᵀ nu'` est bornée, ainsi que les accélérations
et multiplicateurs réduits. Un représentant des multiplicateurs complets
porte les valeurs réduites sur base et zéro ailleurs ; aucune norme minimale
n'est revendiquée. La positivité physique de M reste une obligation distincte.

Cette structure n'est jamais déduite de contraintes numériquement proches.
Une perturbation arbitrairement petite peut changer le rang : le contre-exemple
rationnel du prototype reste publié. Les paramètres mécaniques doivent
justifier la structure annoncée ; la fonction ne certifie pas cette provenance.

## Formats, refus et budgets

Les formats `vinkulum.lineaire.creux.1` et `vinkulum.quotient.structurel.1`
contiennent des rationnels canoniques. Le quotient contient le certificat
du KKT réduit et les données structurelles et de réaction. Les anciens
formats publics et le prototype creux historique restent lisibles.
Après écriture JSON, la relecture utilise la bibliothèque standard seule :

```sh
python -S /chemin/du/paquet/vinkulum/_verification_lineaire.py preuve.json
```

SciPy, fourni par l'extra `verification`, est nécessaire à la génération.
`CertificationImpossible` signale un domaine ou budget refusé, une proposition
impossible ou une contraction non établie. Un refus ne démontre pas une singularité.

Budgets : 4 096 inconnues, 200 000 contributions par matrice, 200 000 termes
du défaut et 2 millions de produits pour LU. Les coefficients de la preuve
linéaire sont limités à 2 048 bits, les exposants de majoration à ±16 384.
Les majorations rationnelles sont dirigées vers le haut avec erreur relative
au plus `2⁻¹²⁷`. Dans le quotient, `1≤r≤n`, `n+r≤4096` et le nombre total
de contraintes est au plus 4 096. Le remplissage peut provoquer un refus ;
aucun coût linéaire universel n'est revendiqué.

## Qualification et base de confiance

Huit documents de la roue isolée figurent dans
`bancs/certification-creuse-0.16.0/` : cinq chaînes natives jusqu'à 2 816
inconnues, une somme de contributions dupliquées et deux familles de rangs
1 et 2, dont une base non contiguë et permutée. Chaque document est relu
dans un processus `python -S`. La conservation des entrées est contrôlée.
Les tests couvrent solutions analytiques aux extrémités des enveloppes,
permutations, refus et falsifications ; les archives historiques sont relues.

La base de confiance comprend Python, ses entiers/Fraction et les modules
de preuve. Générateur et vérificateur partagent certains calculs ; les
oracles analytiques et rationnels sont des contre-épreuves distinctes.
Les empreintes identifient les artefacts, sans authentifier leur origine
mécanique. Aucune preuve Lean des nouveaux certificats n'est revendiquée.
Ils ne certifient ni géométrie, ni trajectoire générale, ni contact, ni stabilité spectrale.

`certifier_initialisation` conserve son contrat dense existant. Les API
creuses reçoivent des matrices explicites ; elles ne remplacent pas
automatiquement la certification des systèmes internes exportés par cette fonction.
