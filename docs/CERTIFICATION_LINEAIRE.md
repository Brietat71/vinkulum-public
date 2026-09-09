# Certificats linéaires : unicité et erreur en avant

Disponible depuis **0.13.0** dans `vinkulum.certification`. La garantie
porte sur le système de nombres fourni, pas sur une trajectoire mécanique.
La certification est facultative et dense : 64 inconnues par défaut,
128 au plus. Un refus signifie que cette méthode n'établit pas la garantie ;
il ne démontre pas à lui seul la singularité du système.

Depuis 0.14.0, l'option `redondances=True` de l'initialisation utilise un
[certificat de quotient distinct](CERTIFICATION_QUOTIENT.md). Le contrat
linéaire ci-dessous et le comportement sans cette option sont conservés.

## 1. Contrat mathématique

Soient une matrice carrée réelle A et des vecteurs b, x dont les
coefficients binary64 sont considérés comme des rationnels exacts.
Le vecteur x est une solution approchée quelconque ; R est une matrice
proposée comme inverse approchée. R n'est pas supposée correcte.

Le générateur calcule exactement :

    r = b − Ax
    E = I − RA
    z = Rr
    s_i = Σ_j |E_ij|
    η = max_i s_i
    ζ = max_i |z_i|

Si `η < 1`, alors **A est inversible** et la solution exacte x* de Ax*=b
est unique. On obtient :

    ||x* − x||_∞ ≤ ζ / (1 − η)

La sortie donne un binary64 B supérieur ou égal à cette borne rationnelle.
Pour chaque coordonnée, elle donne aussi un binary64 h_i tel que :

    |x*_i − x_i| ≤ |z_i| + s_i B ≤ h_i

Ces bornes sont absolues dans les unités des coordonnées fournies. Une
accélération linéaire, une accélération angulaire et un multiplicateur
n'ont pas la même unité physique : le maximum de leurs valeurs numériques
n'est pas une norme physique invariante aux changements d'unités. La
garantie composante par composante reste explicite.

## 2. Démonstration et vérificateur

Comme la norme infinie est sous-multiplicative et `||E||∞=η<1`, la série
`I + E + E² + …` converge et est l'inverse de `I−E=RA`. Ainsi RA est
inversible. Pour A carrée, cela impose que A et R soient inversibles.
Il n'est donc pas nécessaire de faire confiance à l'inversion numérique.

Pour un système KKT dont les lignes inférieures sont de la forme `[G, 0]`,
cette inversibilité implique aussi le plein rang des lignes du G stocké :
une dépendance de ces lignes serait une dépendance des lignes de A. Ce
résultat ne classe pas le rang des systèmes pour lesquels la certification
est refusée, et ne borne pas l'erreur de calcul de G.

Avec `e=x*−x`, l'identité `(I−E)e=z` donne `e=z+Ee`. D'où :

    ||e||∞ ≤ ζ + η ||e||∞
    (1−η) ||e||∞ ≤ ζ

Le facteur `1−η` est strictement positif, ce qui justifie la division.
La majoration de chaque composante suit de
`|e_i| ≤ |z_i| + Σ_j |E_ij| |e_j|`.

Le vérificateur autonome emploie un critère d'inclusion : il recalcule
E et z en rationnels, exige `η<1`, puis vérifie pour chaque ligne :

    |z_i| + s_i B ≤ B
    |z_i| + s_i B ≤ h_i

La première inégalité garantit que l'application `v ↦ z+Ev` conserve la
boule fermée de rayon B. Les itérés depuis zéro y restent, et convergent
par `η<1` vers l'unique solution. Celle-ci est donc dans cette boule.
La seconde inégalité majore alors ses composantes. Ce vérificateur ne
recopie pas la division utilisée pour construire B et n'accepte aucune
déclaration enregistrée de « succès » sans refaire les calculs.

Les produits doubles formant E et r sont évalués par le produit scalaire
entier Rust de la 0.12.2, via son entrée privée ; sa décision de tolérance
n'est pas utilisée. Les produits triples formant z et les bornes utilisent
`fractions.Fraction`. La conversion des majorants en binary64 est
contrôlée par une comparaison rationnelle après arrondi vers le haut.
Un majorant non représentable en binary64 fini est refusé.

Le vérificateur `python/vinkulum/_verification_lineaire.py` n'utilise que
la bibliothèque standard Python, avec un chemin de calcul rationnel
indépendant du produit scalaire Rust. Les nombres sont encodés en
hexadécimal binary64 canonique : leur valeur n'est pas modifiée par une
sérialisation décimale.

**Base de confiance :** opérations entières Rust et `num-bigint`,
arithmétique rationnelle et conversion binary64 de CPython, compilateurs,
matériel, fidélité de l'export natif. La démonstration est mathématique ;
elle n'est pas une preuve formelle du code dans un assistant de preuve.
La relecture indépendante d'un certificat réduit la confiance nécessaire
au générateur. Elle n'authentifie pas la provenance physique de ses données.

## 3. Accélérations réellement calculées par le noyau

`certifier_initialisation(noyau)` prend une copie du modèle et appelle
le même `acc_init` que l'intégrateur. Un témoin facultatif capture les
contributions de A, b et le vecteur x **effectivement résolu**. Il n'est
pas reconstruit depuis une trajectoire et la solution n'est pas remplacée
par celle d'un solveur externe. NumPy propose seulement R dans ce chemin.

La matrice cible est la somme exacte des contributions capturées. Lorsque
des contributions partagent un indice, elles sont additionnées en
rationnels. Si cette somme n'est pas représentable exactement en binary64,
la certification est refusée : aucune erreur d'assemblage supplémentaire
n'est dissimulée par cette conversion. Les contributions figurent dans le
document ; le vérificateur contrôle aussi leur égalité à la matrice cible.

Les variables sont, pour chaque corps, les trois accélérations linéaires
en m/s² puis les trois accélérations angulaires en rad/s², dans le repère
spatial ; viennent ensuite les multiplicateurs de liaison, dans l'ordre
des lignes du noyau. L'indice de séparation est fourni dans `origine`.

Périmètre actuellement refusé :

- matrice vide, au-delà du budget de dimension, ou sans contraction établie ;
- contraintes désactivées ou redondantes sans l'option `redondances=True` ;
- contacts non lisses et sélection de leurs ensembles actifs ;
- partitions gelées ou réduites, et schéma de redémarrage imposé ;
- données non finies, agrégation non représentable, précision demandée non établie.

Le modèle initial reste inchangé, y compris si la certification échoue.
Le certificat ne prouve pas l'exactitude mécanique de M, G, des forces ou
de leur biais d'accélération. Il ne prouve pas non plus que la pose est
admissible : assembler le mécanisme reste une opération distincte.
La certification de la géométrie, des redondances perturbées, des événements et de
l'erreur globale de temps reste ouverte.

## 4. Exemple exécutable et preuve transportable

```python
import json
import numpy as np
from vinkulum import Noyau
from vinkulum.certification import certifier_initialisation, verifier_certificat

n = Noyau([0., -9.81, 0.])
n.corps("pendule", 1., (.2*np.eye(3)).ravel().tolist(),
        [0., -1., 0.], v=[.3, 0., 0.], w=[0., 0., .3])
n.liaison("pivot", None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1])
n.assemble()
preuve = certifier_initialisation(n, erreur_max=1e-11)
controle = verifier_certificat(preuve)
print(float(controle["borne_erreur_inf"]))
with open("preuve.json", "w") as f:
    json.dump(preuve, f, indent=2)
```

Pour vérifier la preuve sans charger Vinkulum, NumPy ni les paquets du site :

```bash
python -I -S python/vinkulum/_verification_lineaire.py preuve.json
```

Une [preuve de ce pendule est conservée](bancs/certificat-pendule-0.13.0.json).
Les données du document définissent le système contrôlé ; pour rattacher
la preuve à un calcul particulier, conserver aussi l'identité de la roue,
les entrées du modèle et la provenance de cet export.

## 5. Contre-épreuves et limites de performance

`test_certification.py` compare les inclusions à une élimination de Gauss
rationnelle indépendante sur des systèmes denses, permutés, indéfinis et
mal conditionnés. Il altère des certificats (solution, matrice, inverse,
second membre, contributions et bornes) et exige le refus. Une borne
abaissée d'un seul ulp sous 1/3 est détectée. Le chemin natif est confronté
à l'accélération enregistrée au départ du premier pas réel.

Deux cas distinguent cette garantie d'un contrôle de résidu :

- Pour `A=diag(1,2⁻⁶⁰)`, `b=(1,2⁻⁶⁰)`, `x=(1,0)`, le résidu est
  inférieur à 10⁻¹², mais l'erreur exacte vaut 1. Le certificat donne 1
  et refuse une exigence de 10⁻¹².
- Pour `A=(2⁻¹⁰⁷⁴)`, `b=(2⁻¹⁰⁷⁴)`, `x=(1)`, `R=(1)`, on a
  `η=1−2⁻¹⁰⁷⁴<1` et une erreur nulle. Arrondir η en binary64 détruirait
  cette décision ; la comparaison rationnelle conserve la preuve.

La [campagne](../ci/audit_certification_lineaire.py) contrôle six
initialisations de pendules, quatre systèmes denses et ces deux cas.
Ses critères de précision sont déclarés avant exécution. Les certificats
et journaux de livraison sont archivés avec leurs empreintes.

La construction d'une inverse dense et les vérifications ont un coût
cubique. La copie du modèle ajoute un coût qui dépend de ses données.
Cette capacité s'exécute à la demande ; aucun certificat linéaire dense
n'est ajouté automatiquement aux pas de simulation. Elle ne constitue
pas encore une méthode de certification adaptée aux très grands systèmes.

## 6. Position scientifique

La certification a posteriori par préconditionnement appartient au calcul
validé établi. Les travaux de
[Minamihata, Ogita, Rump et Oishi (2020)](https://tore.tuhh.de/entities/publication/850096c1-ad54-4348-b236-6c37f555b742)
étudient des bornes plus fines pour les systèmes denses, notamment près
du mauvais conditionnement critique. La présente version implémente le
critère suffisant élémentaire démontré ci-dessus ; elle ne revendique
ni leur algorithme plus avancé ni une rupture de performance.

La suite devra étendre les quotients exacts aux redondances perturbées,
aux opérateurs creux, aux incertitudes des coefficients et aux trajectoires,
avec des certificats dont les hypothèses peuvent elles-mêmes être contrôlées.
