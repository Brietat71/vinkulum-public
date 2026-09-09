# Certificat géométrique d'assemblage

`vinkulum.certification.certifier_assemblage` établit l'existence et l'unicité
d'une pose admissible dans une boîte locale, pour le domaine ci-dessous.
Il conserve la pose native dans le document et borne sa distance à la pose
admissible. L'appel lit un instantané ; il ne modifie pas l'état du noyau.
Cette fonction requiert SciPy, disponible via l'extra `verification` du paquet.

```python
from fractions import Fraction
from vinkulum import Noyau
from vinkulum.certification import certifier_assemblage, verifier_certificat

n = Noyau([0., 0., -9.81])
n.corps('pendule', 1., [1.,0.,0., 0.,1.,0., 0.,0.,1.], [0.,0.,-1.])
n.liaison('pivot', None, 0, bloque_r=[0, 2])
certificat = certifier_assemblage(
    n, jauges=[(5, 0)], rayons=[Fraction(1, 10**6)]*7)
preuve = verifier_certificat(certificat)
assert preuve['existence_unique_locale']
```

## Contrat

Un à quatre corps rigides, liaisons holonomes sans loi imposée à repères
matériels stockés finis et distances de longueur strictement positive. Les
contacts, maillages, poutres, superéléments, corps gelés et les autres
contraintes sont refusés. Les contraintes originales sont toutes conservées,
y compris celles que la résolution dynamique aurait désactivées.

Chaque corps apporte `(rx, ry, rz, qw, qx, qy, qz)`. Les jauges sont des
couples `(indice global, valeur fixée)` ; elles complètent les contraintes
mécaniques et une équation de norme de quaternion par corps pour obtenir un
système carré. Dans l'exemple, l'indice 5 fixe `qy` et donc la liberté du
pivot au voisinage de la pose considérée. Une jauge dépendante ou une
contrainte redondante peut empêcher de prouver l'inclusion : aucun retrait
silencieux de ligne n'est effectué.

Les rayons sont strictement positifs et explicites, dans les unités des
coordonnées : trois translations puis quatre composantes de quaternion.
Une même valeur pour ces sept composantes n'est pas une norme physique
invariante au changement d'unités. Les données entières, `Fraction` et
binary64 finies sont admises ; les booléens, complexes et chaînes sont
refusés par le producteur. Les rationnels encodés sont limités à 8192 bits
pour leur numérateur et leur dénominateur. Un calcul ou un encodage hors
budget est refusé.

Le centre de translation conserve exactement la somme des parties haute
et basse de la position. SciPy propose un quaternion à partir de la matrice
native. Cette conversion n'est pas supposée exacte : les distances de
rotation sont bornées relativement à la matrice native conservée, en
incluant cet écart. NumPy propose un inverse du Jacobien ; son résultat est
soumis aux vérifications rationnelles.

## Garantie et vérification

Le [prototype qualifié](ASSEMBLAGE_CERTIFIE_PROTOTYPE.md) décrit les
équations, les hypothèses et les contre-épreuves antérieures à l'API.
Les fonctions de rotation sont polynomiales en quaternion, avec équation de
norme unitaire. Une inclusion stricte de Krawczyk sur toute la boîte et un
contrôle de branche des liaisons établissent existence et unicité locales.
Cela ne prouve ni une unicité globale, ni la convergence du Newton natif.

Le générateur utilise des jets d'intervalles rationnels. Le vérificateur
reconstruit des polynômes, dérive leurs monômes et refait une majoration par
rayons. Les deux partagent la définition mécanique et les primitives
d'intervalle ; la seconde dérivation ne constitue pas une seconde
spécification indépendante de la mécanique.

Le document `vinkulum.assemblage.1` contient le modèle, les jauges, le
centre, les rayons, l'inverse proposé, la pose native initiale et les bornes
annoncées. La vérification recalcule l'inclusion et les distances. Elle
refuse une borne annoncée trop petite, même si le reste du certificat est
valide. Les bornes de distance sont :

- valeurs absolues des trois écarts de translation par corps ;
- carré de la norme de Frobenius de l'écart des matrices de rotation.

La seconde quantité n'est pas un angle ni une erreur de quaternion.
Les documents peuvent être sérialisés en JSON. La relecture autonome ne
requiert ni NumPy, ni SciPy, ni le binaire de Vinkulum :

```sh
python -S chemin/vers/vinkulum/_verification_lineaire.py certificat.json
```

Les modules `_verification_assemblage.py`, `_oracle_assemblage.py` et
`_geometrie_certifiee.py` doivent rester à côté de ce vérificateur. Les
schémas linéaire et quotient précédents sont toujours reconnus.

## Limites de la preuve

Les points d'attache et matrices de repère stockés sont considérés comme
des rationnels exacts. Les matrices de repère sont utilisées telles quelles ;
l'appartenance exacte à SO(3) concerne les rotations inconnues des corps.
Le certificat ne propage pas l'erreur de leur construction à partir des
paramètres utilisateur et n'authentifie pas la provenance du document.
Il certifie les données conservées, pas un objet natif qui aurait changé
après la photographie. Il ne garantit ni forces, ni accélérations, ni
trajectoires, ni contacts.

La base de confiance inclut CPython/Fraction, le modèle de contraintes,
l'arithmétique d'intervalle, le vérificateur et l'export Rust avec ses outils
de compilation. Le théorème géométrique utilisé n'est pas formalisé dans
Lean. Une inclusion refusée signifie que ces calculs n'ont pas établi la
preuve demandée ; elle ne démontre pas l'absence de solution physique.
