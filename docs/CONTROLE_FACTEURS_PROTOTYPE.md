# Contrôle par facteurs : accélérer les majorations sans modifier les champs

Cette note conserve la campagne de 60 essais. Le
[certificat binaire avec séparateurs](INERTIE_BINAIRE_SEPARATEURS.md)
apporte ensuite 36 essais : à 512 poutres, le même service avec majorants
prend 0,516 s contre 0,673 s avec le certificat Decimal et 0,715 s pour
la LU de ce nouveau lot. Les résultats ci-dessous ne sont pas réécrits.

Expérience du 8 septembre 2026, après le
[certificat par inertie dirigée](INERTIE_CONTRAINTE_PROTOTYPE.md).
Le contrôle par fréquence devenait la plus grosse phase du service avec
majorants. Le nouveau prototype conserve sa résolution et ses normes
physiques ; il réduit le coût de deux majorations d'opérateur.
Le transfert dans l'API publique reste à réaliser : la version est 0.11.0.

## De la mesure au choix mathématique

Le [profil exploratoire](bancs/controle-facteurs-sondes-2026/profil/README.md)
attribue l'essentiel du surcoût aux normes et aux images physiques, plutôt
qu'aux petits solveurs. Sur 512 poutres, les deux évaluations répétées de
la norme spectrale de DX représentent environ 0,097 s sous instrumentation ;
l'application des deux images de réparation et leurs normes, environ 0,080 s.
Le solve de Y prend environ 0,009 s. Ces attributions instrumentées ne sont
pas des gains comparatifs et ne s'additionnent pas au chronométrage ordinaire.

On note `X=T−E_i W Y`, où Y dépend de la fréquence. Le Schur conservé
emploie déjà `G=(DX).T DX`. Deux observations guident le changement.

**La norme d'opérateur est une grande valeur propre.** Pour une matrice
symétrique G, `lambda_max(G) <= ||G||_inf`. Une majoration du défaut de
formation du Gram permet d'utiliser cette inégalité pour majorer `||DX||₂`,
sans SVD physique supplémentaire. On peut aussi utiliser sa trace, car
le Gram exact est positif. Le prototype prend le minimum des deux
majorants évalués, avec des budgets d'arrondi pour les produits et sommes.
Il ne remplace aucune petite norme de charge par `sqrt(q.T G q)` : cette
dernière formule peut perdre une combinaison presque nulle.

**Le correctif provient de peu de contraintes.** Avec
`H=C J`, `J=(B.T C)^−1 B.T W`, les images massique et énergétique de H
possèdent des facteurs associés au nombre de contraintes. On prépare
une QR de ces seuls facteurs gauches, puis une petite matrice L.
Pour l'image stockée A, on évalue explicitement le défaut `A−Q L`.
Il couvre aussi les arrondis qui empêchent A d'être exactement de faible
rang. Aucune petite singularité n'est supprimée au moyen d'un seuil.

Si `||Q||₂ <= kappa` et les normes des colonnes du défaut sont majorées
par le vecteur e, alors, pour chaque v,

```text
||A v||₂ <= kappa ||L v||₂ + e.T |v|.
```

Les défauts de calcul du petit produit et de l'ancien produit physique
s'ajoutent. La multiplication précède la norme ; le défaut reste présent
lorsque les facteurs prédisent une compensation exacte. Le calcul en
ligne dépend alors du nombre de contraintes et de directions réduites,
au lieu du nombre de coordonnées physiques.

Le [carnet de preuves](CONTROLE_FACTEURS_PREUVES.md) distingue ces
identités exactes de leur évaluation. Le prototype
[ci/controle_facteurs.py](../ci/controle_facteurs.py) conserve aussi une
option QR pour la norme d'extension. Le
[pilote à neuf essais](bancs/controle-facteurs-sondes-2026/README.md)
compare les deux options ; Gram est choisi avant la campagne formelle.
Les sources et résultats de cette alternative sont conservés.

La revue du premier prototype a également produit deux réfutations :
une norme d'audit pouvait perdre des petits carrés lors de sa sommation,
et un Gram nul ou subnormal pouvait cacher une image non nulle. Leurs
[données et ancien code](bancs/controle-facteurs-sondes-2026/refutations-v1/README.md)
restent archivés. Avant la campagne, les normes d'audit ont reçu une mise
à l'échelle, un budget relatif et absolu des opérations, puis des arrondis
scalaires vers l'extérieur. Une diagonale Gram nulle ou subnormale entraîne
un repli sur une norme de Frobenius physique robuste. Le second pilote
recompte le coût de ces corrections ; il conserve les neuf champs identiques
et les contrôles acceptés. Cela ne certifie pas tous les autres calculs.

## Ce qui reste identique, ce qui peut changer

La construction du complément, le certificat d'inertie, Y, X, le Schur,
sa factorisation, les coordonnées et la reconstruction des champs gardent
les opérations historiques. Les normes par charge sont toujours calculées
sur les champs physiques retournés avec D et M originaux. La norme du
résidu de Y est simplement réutilisée au lieu d'être calculée deux fois.

Les normes d'extension et de réparation deviennent des majorants
supplémentaires. Elles peuvent donc augmenter les bornes de réponse ou
provoquer un refus de contrôle que l'ancienne évaluation aurait évité.
Une proximité de résonance globale ne peut pas être effacée par ce changement.
Le candidat conserve toutes les décisions de marge et toutes les sorties
de refus ; une borne indisponible reste indisponible.

Les champs identiques dans les sondes et contre-épreuves ne constituent
pas à eux seuls une preuve universelle d'identité entre bibliothèques et
plateformes. La campagne compare aussi indépendamment les champs aux
oracles, avec la tolérance physique de 10⁻⁶.

## Protocole comparatif

Le lot comporte 60 processus frais : trois consoles natives de 32, 128 et
512 poutres, cinq variantes, une chauffe et trois mesures par configuration.
L'ordre s'inverse à chaque passage. Chaque essai traite 257 fréquences de
0 à 40 Hz et six charges terminales. CPU 8 et bibliothèques numériques à
un fil ; environnements et sources sont figés et identifiés par empreintes.

Les variantes sont le nouveau contrôle par facteurs, le précédent contrôle
avec inertie, la LU corrigée et les deux projections HCB historiques
d'Exudyn 1.11.0. Les deux candidats recalculent QR, direction intérieure,
B, certificat d'inertie à 32 chiffres et seuil gamma associé à 80 Hz,
quatre blocs de Krylov et contrôle de profondeur huit. Le nouveau candidat
emploie explicitement `extension="gram"`. Aucun de ces coûts n'est offert.

Le chronomètre comprend préparation, réponses et normes physiques.
Imports, lectures, sauvegardes et juge indépendant sont hors chronomètre.
Les empreintes calculées à l'intérieur de l'appel du certificat restent
comptées par son chronomètre extérieur. Les rangs HCB 185, 761 et 768 sont
ceux de la confrontation antérieure ; leur recherche est offerte aux
témoins, mais la construction de leurs bases est recomptée.

Le juge évalue champs en masse, déformations et ports, par charge et sur
toutes les combinaisons des six charges. Les majorants par charge sont
confrontés aux erreurs absolues, avec la marge flottante documentée
`64 eps × norme de référence`. Les références Decimal 70/90 coïncident
après conversion en double ; elles ne sont pas des oracles par intervalles.
Toutes les occurrences, y compris les refus de précision HCB, sont retenues.

## Résultats de la campagne

Les **24 essais des deux contrôleurs** passent le juge et leurs majorants.
Au total, **50/60 essais** et **12/15 configurations** sont admis. Les dix
refus sont conservés : deux mesures HCB standard à 128 poutres et les
huit occurrences HCB à 512 poutres. Une configuration avec un seul refus
est exclue des ratios de vitesse à précision commune.

Durées en secondes, médianes des trois mesures ; l'admission comprend la
chauffe. Les médianes des phases ne s'additionnent pas nécessairement à
celle du total.

| Poutres | Variante | Préparation | Réponses | Total | Admissions |
|---:|---|---:|---:|---:|---:|
| 32 | Facteurs, contrôle | 0,0509 | 0,1039 | 0,1548 | 4/4 |
| 32 | Ancien contrôle | 0,0495 | 0,1220 | 0,1712 | 4/4 |
| 32 | LU corrigée | 0,0011 | 0,1026 | 0,1037 | 4/4 |
| 32 | HCB standard | 0,0305 | 0,0775 | 0,1079 | 4/4 |
| 32 | HCB énergie | 0,0308 | 0,0779 | 0,1087 | 4/4 |
| 128 | Facteurs, contrôle | 0,1245 | 0,1288 | 0,2533 | 4/4 |
| 128 | Ancien contrôle | 0,1233 | 0,1748 | 0,2980 | 4/4 |
| 128 | LU corrigée | 0,0012 | 0,2210 | 0,2222 | 4/4 |
| 128 | HCB standard | 0,6643 | 0,9698 | 1,6350 | 2/4 |
| 128 | HCB énergie | 0,6603 | 0,9956 | 1,6542 | 4/4 |
| 512 | Facteurs, contrôle | 0,4394 | 0,2297 | 0,6712 | 4/4 |
| 512 | Ancien contrôle | 0,4371 | 0,3820 | 0,8191 | 4/4 |
| 512 | LU corrigée | 0,0019 | 0,6939 | 0,6959 | 4/4 |
| 512 | HCB standard | 5,0153 | 1,5513 | 6,5536 | 0/4 |
| 512 | HCB énergie | 4,9559 | 1,5307 | 6,5088 | 0/4 |

Par rapport à l'ancien contrôle, les réponses accélèrent de **1,17×,
1,36× et 1,66×** ; les totaux, de **1,11×, 1,18× et 1,22×**. La nouvelle
préparation coûte un peu plus cher : la factorisation et ses audits sont
comptés. Sur 512 poutres, les réponses contrôlées prennent **0,230 s**,
contre **0,382 s** auparavant.

Le total à 512 poutres devient légèrement meilleur que la LU :
**0,671 contre 0,696 s**, soit environ **3,6 % de temps en moins**.
Les plages des trois mesures sont respectivement 0,6690–0,6727 s et
0,6951–0,7002 s. Ce petit gain sur ce lot ne démontre pas une avance
universelle. À 32 et 128 poutres, le nouveau service reste plus lent que
la LU : 0,155 contre 0,104 s, puis 0,253 contre 0,222 s.

À 128 poutres, HCB énergie satisfait la précision : le nouveau contrôle
est **6,53× plus rapide au total** (0,253 contre 1,654 s). À 32 poutres,
HCB reste plus rapide. À 512 poutres, ses deux configurations dépassent
la tolérance : aucun ratio de vitesse HCB n'est publié. HCB standard
à 128 poutres atteint 1,662 × 10⁻⁶ au pire, arrondi vers le haut, et ne
passe qu'une des trois mesures ; sa chauffe passe.

Les pires erreurs des contrôleurs restent respectivement inférieures à
5,58 × 10⁻¹¹, 4,17 × 10⁻¹⁰ et 2,17 × 10⁻⁹. Les champs ont les mêmes
empreintes entre les deux contrôleurs pour chaque maillage et passage.
Les facteurs changent les majorants, pas la réponse observée. Les normes
physiques restent identiques. Les bornes conservent leur statut numérique,
même lorsqu'elles couvrent tous les oracles du lot.

La préparation redevient le coût dominant sur le grand cas : environ
0,439 s sur 0,671 s, dont 0,197 s de certificat d'inertie, 0,086 s de
condensation et 0,090 s de préparation du contrôle. Supprimer encore du
travail répété dans les réponses et réduire la préparation sont les deux
axes à comparer, avec les mêmes services et budgets. La compression complète
des métriques demanderait une nouvelle preuve et une nouvelle campagne.

## Portée et prochaines étapes

Le certificat d'inertie continue de démontrer la coercivité du complément
discret exact défini par D, M et B stockés. Les nouveaux budgets d'arrondi
du contrôleur sont évalués en doubles, avec les hypothèses usuelles
d'absence de dépassement et de sous-flux non couvert. Les normes, sommes,
racines, QR et majorants de ces budgets ne sont pas eux-mêmes encadrés par
intervalles : **les réponses ne sont pas certifiées machine**.

Les contre-épreuves couvrent les deux variantes, plusieurs contraintes,
des masses couplées et une masse complète semi-définie positive avec ports
sans masse. La masse intérieure du contrôleur reste limitée aux composantes
connexes de taille au plus six. Les champs et leurs petites normes ne sont
pas remplacés par des métriques réduites ; la reconstruction physique garde
donc un coût proportionnel à la taille du modèle.

La [lecture de six sources primaires](CONTROLE_FACTEURS_SOURCES.md) précise
l'antériorité des estimateurs réduits stables de 2014, les vérifications de
QR de 2024 et les travaux de 2025–2026 sur les résidus. Elle expose aussi
une compression plus ambitieuse des métriques et ses défauts de Schur.
Cette extension ne fait pas partie du candidat mesuré. L'interpolation
empirique et les estimateurs aléatoires ne remplacent pas le contrat
déterministe actuel.

Les contre-épreuves Fraction calculent indépendamment les champs exacts,
les erreurs et les majorations d'opérateur. Un test perturbe l'image de
réparation hors des facteurs proposés : une compensation qui annule leur
action doit garder un majorant non nul. Un autre rend le minorant spectral
valide mais inutilisable pour le contrôle : les bornes restent refusées.
Dix-sept tests de preuve et onze tests d'implémentation couvrent ces contrats ;
l'un vérifie 48 normes contre leurs carrés exacts, jusqu'aux valeurs
subnormales et à l'échelle 2¹⁰⁰⁰.

Le banc reste une confrontation de réponses matérielles linéaires avec
bases HCB officielles et pont local de réponse. Il ne mesure ni une
simulation FFRF complète, ni une dynamique non linéaire avec contacts,
ni Simpack. Aucun classement mémoire commun n'en est déduit.
L'objectif de supériorité généraliste reste actif.

## Reproduction

L'[archive formelle](bancs/controle-facteurs-2026/README.md) conserve les
résultats, sources, bases, diagnostics des facteurs et empreintes des
grands champs exclus. Son validateur requalifie les diagnostics physiques
et les majorants, puis rejoue les certificats d'inertie distincts depuis
D/M/B. Il ne reconstitue pas les grands champs ni les transformations du
contrôleur à partir des seuls diagnostics archivés.

```bash
python ci/archive_controle_facteurs.py --verifier docs/bancs/controle-facteurs-2026
```

Une reproduction exige les environnements et données aux identités
publiées, puis un dossier neuf ; aucun essai interrompu n'est relancé
pour remplacer son résultat :

```bash
"$PY_VINKULUM" docs/bancs/controle-facteurs-2026/sources/experience_controle_facteurs.py \
  --campagne "$SORTIE_NEUVE" "$ENTREES_AVEC_ORACLES" \
  "$PY_VINKULUM" "$PY_VINKULUM" "$PY_EXUDYN" \
  --historique docs/bancs/confrontation-ports-exudyn-0.10.0/bilan.json
```
