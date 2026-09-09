# Prototype de certificat géométrique

Ce document conserve la qualification du prototype antérieure à l'API.
L'intégration publique est décrite dans [le contrat d'assemblage](CERTIFICATION_ASSEMBLAGE.md).
Elle étend aussi l'export aux repères matériels stockés finis ; le domaine
identité décrit ci-dessous est celui du prototype initial. Les certificats linéaires et de
quotient existants gardent leurs schémas et leur vérificateur.

## Objet mathématique

Pour chaque corps, les inconnues sont trois translations et un quaternion
`(w,x,y,z)`. Une équation `w²+x²+y²+z²=1` est ajoutée. Les matrices de
rotation sont des polynômes ; elles appartiennent à SO(3) à toute racine
admise. Les repères matériels de liaison du premier domaine sont l'identité.
Les points d'attache stockés sont interprétés comme des rationnels binary64
exacts. Les parties basse et haute des translations natives sont additionnées
rationnellement, sans perte lors de l'export.

Les résidus de liaison reprennent les translations dans le repère A et la
partie antisymétrique de `R_Aᵀ R_B`. Pour une distance de longueur strictement
positive, le résidu est la distance au carré moins le carré de la longueur :
ses zéros sont ceux de la contrainte de norme. Les normes de quaternion et
les jauges explicites complètent le système carré. Aucune ligne mécanique
originale n'est supprimée. Une dépendance peut donc provoquer un refus.

Le générateur évalue fonctions et Jacobien en intervalles rationnels par
règles de différentiation automatique. L'oracle collecte des polynômes,
les dérive symboliquement et refait l'inclusion par des rayons et magnitudes.
Les deux calculs partagent la définition des fonctions mécaniques et les
primitives d'intervalle ; ils ne constituent pas deux spécifications
indépendantes de la mécanique.

Pour un centre `c`, des rayons positifs `r`, une matrice proposée `R` et
`X=[-r,r]`, les opérations exactes vérifient

```
-R f(c) + (I - R J(c+X)) X ⊂ intérieur(X).
```

Le programme vérifie également une contraction en norme infinie pondérée
par les rayons. Ces conditions établissent l'existence et l'unicité locale,
pas une unicité globale ni un bassin de convergence du Newton natif.
La matrice proposée par NumPy n'est jamais tenue pour une preuve.
Le [théorème 13.3 de Rump, Acta Numerica 2010](https://www.tuhh.de/ti3/rump/intlab/ActaNumerica2010.pdf)
donne le résultat d'inclusion utilisé. Ce théorème géométrique n'est pas
encore formalisé dans le projet Lean.

Pour les liaisons bloquant deux ou trois rotations, le critère de branche
utilisé par le noyau doit être strictement positif sur la boîte entière.
Cela exclut la racine parasite à 180° dans le domaine qualifié.

## Domaine et limites

L'export privé est une copie en lecture seule : un à quatre corps, liaisons
holonomes sans loi imposée et distances positives. Les contacts, maillages,
corps gelés, poutres, superéléments, repères matériels non identité et autres
types de contraintes sont refusés explicitement. Une interface généraliste
ne doit pas cacher ces refus derrière une sélection de lignes.

Les documents expérimentaux certifient leurs données. Ils n'authentifient
pas leur provenance native et ne propagent pas l'erreur de construction des
points d'attache à partir des paramètres initiaux de l'utilisateur. La
preuve ne concerne ni des forces, ni une accélération, ni une trajectoire.
La base de confiance inclut CPython/Fraction, les opérations d'intervalle,
la définition des contraintes, le vérificateur, l'export Rust et sa chaîne
de compilation. NumPy et SciPy proposent un inverse et un centre quaternion ;
ils ne décident pas l'inclusion.

Une fonction privée borne aussi la distance aux translations et matrices de
rotation natives. Elle inclut l'écart introduit par la conversion de matrice
en quaternion ; sa borne de rotation est un carré de norme de Frobenius,
pas une distance angulaire. Elle est désormais intégrée au schéma public, qui conserve aussi la pose native.

## Vérifications conservées

- Produits d'intervalles contre tous les sommets et milieux de 784 paires.
- Solutions analytiques, poses et normes perturbées, refus de singularité,
  branche retournée, jauge absente, inverse altéré et boîte insuffisante.
- Comparaison des dérivées symboliques et des jets sur vingt états rationnels
  d'un mécanisme à deux corps.
- Comparaison avec `phi` et `phi_dot` natifs sur 32 poses hors équilibre,
  avec seuil fixé à `1e-12`. C'est une contre-épreuve, pas une preuve du Rust.
- Pendules à trois échelles de longueur, double pendule, trilatération,
  conservation des lignes redondantes et contrôle de l'état natif inchangé.
- Documents altérés refusés ; relecture dans un processus Python `-S`.

Quatre documents sont conservés dans
[assemblage-certifie-prototype](bancs/assemblage-certifie-prototype/).
Exemple de vérification autonome avec la bibliothèque standard :

```sh
python -S ci/archive_assemblage_certifie.py docs/bancs/assemblage-certifie-prototype/double-pendule.json
```

Les quatre suites `ci/test_*assemblage_certifie*.py` et
`ci/test_export_geometrie_certifiee.py` sont raccordées à la CI locale.
La qualification de livraison sur roue isolée est consignée séparément.
Les limitations historiques de ce prototype restent explicites ci-dessus.
